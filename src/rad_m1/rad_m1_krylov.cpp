//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_krylov.cpp
//! \brief the implicit-M1 Krylov loop on SEVERAL ranks (tests_m1/runs_3w_krylov).
//! implicit_krylov_pipe defaults OFF; implicit_halo_mpi defaults ON where its
//! preconditions hold (ImplicitInit; bitwise the ordinary exchange either way).
//!
//!  * <rad_m1>/implicit_halo_mpi (needs implicit_halo_direct): the ghost zones of the
//!    implicit exchanges on a mesh whose neighbours are NOT all on this rank.  The
//!    on-rank neighbours are filled by the copy kernel of implicit_halo_direct (hd_src
//!    holds -1 for an off-rank neighbour), the off-rank ones by ONE message per
//!    neighbour rank: one pack kernel over a precomputed cell list, one MPI message pair
//!    per neighbour rank, one unpack kernel.  The same numbers the pack / exchange /
//!    unpack chain of the general boundary machinery delivers (pure copies), so the
//!    solve is bitwise the implicit_halo_direct = true multi-rank path.
//!  * <rad_m1>/implicit_krylov_pipe (needs implicit_krylov_fuse = 3): the pipelined
//!    right-preconditioned BiCGStab of Cools & Vanroose (Parallel Computing 65, 2017,
//!    Alg. 4, "p-BiCGStab").  Two reductions per iteration as before, but NEITHER blocks
//!    the device: each is reduced into pinned host memory and the host waits only for
//!    that kernel (an event), while the preconditioner and operator application that
//!    follow are already queued; on several ranks the partial sums then go into ONE
//!    MPI_Iallreduce, posted before the halo exchange and completed after the operator
//!    is launched.  Same iterates as BiCGStab in exact arithmetic; the recurrences for
//!    the auxiliary vectors (w = A M^-1 r, s = A M^-1 p, z = A M^-1 s, ...) accumulate
//!    their own round-off, so convergence is always confirmed on the TRUE residual.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#if defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace radm1 {

namespace {
//! five sums and one max (|r| >= 0, so 0 is the identity of the max slot)
struct M1PipeVal {
  Real s[5];
  Real mx;
};

template <class Space>
struct M1PipeRed {
 public:
  using reducer = M1PipeRed<Space>;
  using value_type = M1PipeVal;
  using result_view_type = Kokkos::View<value_type, Space, Kokkos::MemoryUnmanaged>;

 private:
  result_view_type value;
  bool refs;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit M1PipeRed(value_type &v) : value(&v), refs(true) {}
  KOKKOS_INLINE_FUNCTION
  explicit M1PipeRed(const result_view_type &v) : value(v), refs(false) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    for (int q = 0; q < 5; ++q) {d.s[q] += s.s[q];}
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    for (int q = 0; q < 5; ++q) {v.s[q] = 0.0;}
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return refs;}
};

//! the flattened (m,k,j,i) index, as par_for flattens it
KOKKOS_INLINE_FUNCTION
void M1PipeIdx(const int idx, const int nkji, const int nji, const int ni,
               int &m, int &k, int &j, int &i) {
  m = idx/nkji;
  int r = idx - m*nkji;
  k = r/nji;
  r -= k*nji;
  j = r/ni;
  i = r - j*ni;
}

//! wait for the work queued on the device up to record(), not for what came after it
class M1Evt {
 public:
#if defined(KOKKOS_ENABLE_HIP)
  void record() {
    if (ev_ == nullptr) {(void)hipEventCreateWithFlags(&ev_, hipEventDisableTiming);}
    (void)hipEventRecord(ev_, DevExeSpace().hip_stream());
  }
  void wait() {(void)hipEventSynchronize(ev_);}

 private:
  hipEvent_t ev_ = nullptr;
#elif defined(KOKKOS_ENABLE_CUDA)
  // m1-perf-0928: the CUDA version (was a full device fence, which also waited for the
  // interior operator queued behind the pack and so serialised implicit_halo_overlap)
  void record() {
    if (ev_ == nullptr) {(void)cudaEventCreateWithFlags(&ev_, cudaEventDisableTiming);}
    (void)cudaEventRecord(ev_, DevExeSpace().cuda_stream());
  }
  void wait() {(void)cudaEventSynchronize(ev_);}

 private:
  cudaEvent_t ev_ = nullptr;
#else
  void record() {}
  void wait() {DevExeSpace().fence();}
#endif
};

M1Evt &M1EvtA() {static M1Evt e; return e;}
M1Evt &M1EvtB() {static M1Evt e; return e;}
//! implicit_halo_mpi: the unpack kernel of the last exchange (it reads hm_rbuf, which the
//! next exchange's receives overwrite) and whether one is still to be waited for
M1Evt &M1EvtU() {static M1Evt e; return e;}
bool &M1HmUnpackPending() {static bool b = false; return b;}
#if MPI_PARALLEL_ENABLED
//! implicit_halo_mpi: the requests of the exchange in flight (Post -> Finish)
std::vector<MPI_Request> &M1HmReqR() {static std::vector<MPI_Request> v; return v;}
std::vector<MPI_Request> &M1HmReqS() {static std::vector<MPI_Request> v; return v;}
#endif

#if MPI_PARALLEL_ENABLED
void M1PipeOpFn(void *in, void *inout, int *len, MPI_Datatype *) {
  Real *a = static_cast<Real *>(in);
  Real *b = static_cast<Real *>(inout);
  for (int q = 0; q < *len; ++q) {
    for (int c = 0; c < 5; ++c) {b[6*q + c] += a[6*q + c];}
    b[6*q + 5] = (a[6*q + 5] > b[6*q + 5]) ? a[6*q + 5] : b[6*q + 5];
  }
}

//! the datatype, the operation and a private communicator of the pipelined reductions
struct M1PipeMPI {
  MPI_Datatype typ = MPI_DATATYPE_NULL;
  MPI_Op op = MPI_OP_NULL;
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Request req = MPI_REQUEST_NULL;
  Real loc[6], glb[6];
};
M1PipeMPI &M1Pipe() {
  static M1PipeMPI p;
  if (p.op == MPI_OP_NULL) {
    MPI_Type_contiguous(6, MPI_ATHENA_REAL, &p.typ);
    MPI_Type_commit(&p.typ);
    MPI_Op_create(&M1PipeOpFn, 1, &p.op);
    MPI_Comm_dup(MPI_COMM_WORLD, &p.comm);
  }
  return p;
}
#endif

//! start the global reduction of v (nothing to do on one rank)
void M1PipePost(const M1PipeVal &v) {
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks == 1) {return;}
  M1PipeMPI &p = M1Pipe();
  for (int c = 0; c < 5; ++c) {p.loc[c] = v.s[c];}
  p.loc[5] = v.mx;
  MPI_Iallreduce(p.loc, p.glb, 1, p.typ, p.op, p.comm, &p.req);
#else
  (void)v;
#endif
}

//! ...and finish it into v
void M1PipeFinish(M1PipeVal &v) {
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks == 1) {return;}
  M1PipeMPI &p = M1Pipe();
  MPI_Wait(&p.req, MPI_STATUS_IGNORE);
  for (int c = 0; c < 5; ++c) {v.s[c] = p.glb[c];}
  v.mx = p.glb[5];
#else
  (void)v;
#endif
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloMPIInit
//! \brief implicit_halo_mpi: the cell lists of the off-rank part of the implicit
//! exchanges.  Receiver ghost region (G, o) of MeshBlock G towards direction o is filled
//! from the active cells of the same-level neighbour g, cell (k,j,i) <- g's cell
//! (k - o3 nx3, j - o2 nx2, i - o1 nx1) (the map of ImplicitHaloDirect).  On each rank
//! pair both sides enumerate the regions sorted by (receiver gid, o) and each region in
//! (k,j,i) order, so the message needs no header.  hm_state = -1 (the ordinary
//! exchange keeps running) unless EVERY rank has same-level neighbours only, outside the
//! cubed-sphere and polar transforms.

void RadiationM1::ImplicitHaloMPIInit() {
  hm_state = -1;
#if MPI_PARALLEL_ENABLED
  auto *pm = pmy_pack->pmesh;
  int ok = (pm->multilevel || pm->use_cubed_sphere || pm->use_polar_boundary) ? 0 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int n1 = nx1 + 2*ng;
  const int n2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  const bool md = (nx2 > 1), td = (nx3 > 1);
  const int e2 = md ? 1 : 0;
  const int e3 = td ? 1 : 0;
  auto &nb = pmy_pack->pmb->nghbr;
  auto &lev = pmy_pack->pmb->mb_lev;
  const int me = global_variable::my_rank;
  // (neighbour rank, receiver gid, receiver direction, local m, o1, o2, o3)
  using Reg = std::tuple<int, int, int, int, int, int, int>;
  std::vector<Reg> sr, rr;
  int hm_mpif[6] = {0, 0, 0, 0, 0, 0};   // faces with an off-rank ghost (the pack)
  for (int m = 0; m < nmb; ++m) {
    for (int o3 = -e3; o3 <= e3; ++o3) {
      for (int o2 = -e2; o2 <= e2; ++o2) {
        for (int o1 = -1; o1 <= 1; ++o1) {
          if (o1 == 0 && o2 == 0 && o3 == 0) continue;
          int n = NeighborIndex(o1, o2, o3, 0, 0);
          if (n < 0 || n >= pmy_pack->pmb->nnghbr) {ok = 0; continue;}
          const NeighborBlock &q = nb.h_view(m,n);
          if (q.gid < 0) continue;
          if (q.lev != lev.h_view(m)) {ok = 0; continue;}
          if (q.rank == me) continue;
          const int oi = (o1+1) + 3*(o2+1) + 9*(o3+1);
          const int oin = (1-o1) + 3*(1-o2) + 9*(1-o3);
          rr.emplace_back(q.rank, pmy_pack->gids + m, oi, m, o1, o2, o3);
          if (o1 != 0) {hm_mpif[(o1 < 0) ? 0 : 1] = 1;}
          if (o2 != 0) {hm_mpif[(o2 < 0) ? 2 : 3] = 1;}
          if (o3 != 0) {hm_mpif[(o3 < 0) ? 4 : 5] = 1;}
          sr.emplace_back(q.rank, q.gid, oin, m, o1, o2, o3);
        }
      }
    }
  }
  {int g = ok;
  MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  ok = g;}
  if (ok != 1) {
    if (me == 0) {
      std::cout << "<rad_m1> implicit_halo_mpi: mesh not supported (levels / seams), "
                << "the ordinary exchange is used" << std::endl;
    }
    return;
  }
  std::sort(sr.begin(), sr.end());
  std::sort(rr.begin(), rr.end());
  // the neighbour ranks (the same set both ways on a same-level mesh)
  hm_rank.clear();
  for (auto &t : sr) {hm_rank.push_back(std::get<0>(t));}
  for (auto &t : rr) {hm_rank.push_back(std::get<0>(t));}
  std::sort(hm_rank.begin(), hm_rank.end());
  hm_rank.erase(std::unique(hm_rank.begin(), hm_rank.end()), hm_rank.end());
  const int nseg = static_cast<int>(hm_rank.size());
  hm_soff.assign(nseg, 0); hm_slen.assign(nseg, 0);
  hm_roff.assign(nseg, 0); hm_rlen.assign(nseg, 0);
  auto rng = [&](int o, int s, int e, int nx, bool on, int &lo, int &hi) {
    if (!on) {lo = s; hi = e; return;}
    (void)nx;
    if (o < 0) {
      lo = s - ng; hi = s - 1;
    } else if (o > 0) {
      lo = e + 1; hi = e + ng;
    } else {
      lo = s; hi = e;
    }
  };
  std::vector<int> sm, skji, sseg, rm, rkji, rseg;
  // one pass per list: sender cells of the receiver's ghost region, and ghost cells
  for (int pass = 0; pass < 2; ++pass) {
    auto &lst = (pass == 0) ? sr : rr;
    auto &vm = (pass == 0) ? sm : rm;
    auto &vk = (pass == 0) ? skji : rkji;
    auto &vs = (pass == 0) ? sseg : rseg;
    auto &off = (pass == 0) ? hm_soff : hm_roff;
    auto &len = (pass == 0) ? hm_slen : hm_rlen;
    int seg = -1;
    for (auto &t : lst) {
      int r = std::get<0>(t);
      int sg = static_cast<int>(std::lower_bound(hm_rank.begin(), hm_rank.end(), r)
                                - hm_rank.begin());
      if (sg != seg) {
        seg = sg;
        off[seg] = static_cast<int>(vm.size());
      }
      const int m = std::get<3>(t);
      const int o1 = std::get<4>(t), o2 = std::get<5>(t), o3 = std::get<6>(t);
      // the receiver's direction towards the sender
      const int r1 = (pass == 0) ? -o1 : o1;
      const int r2 = (pass == 0) ? -o2 : o2;
      const int r3 = (pass == 0) ? -o3 : o3;
      int il, iu, jl, ju, kl, ku;
      rng(r1, is, ie, nx1, true, il, iu);
      rng(r2, js, je, nx2, md, jl, ju);
      rng(r3, ks, ke, nx3, td, kl, ku);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            int kk = k, jj = j, ii = i;
            if (pass == 0) {   // the sender's own cell
              kk = k - r3*nx3;
              jj = j - r2*nx2;
              ii = i - r1*nx1;
            }
            vm.push_back(m);
            vk.push_back((kk*n2 + jj)*n1 + ii);
            vs.push_back(seg);
          }
        }
      }
      len[seg] = static_cast<int>(vm.size()) - off[seg];
    }
  }
  hm_nsend = static_cast<int>(sm.size());
  hm_nrecv = static_cast<int>(rm.size());
  hm_nq = std::max(std::max(M1_NHALO_T, M1_NHALO_Q), M1_NVIMP_X);
  auto mk = [](const char *nm, const std::vector<int> &v) {
    DualArray1D<int> a(nm, std::max(static_cast<int>(v.size()), 1));
    for (size_t q = 0; q < v.size(); ++q) {a.h_view(q) = v[q];}
    a.modify_host();
    a.sync_device();
    return a;
  };
  hm_sm = mk("m1_hm_sm", sm);
  hm_skji = mk("m1_hm_skji", skji);
  hm_sseg = mk("m1_hm_sseg", sseg);
  hm_rm = mk("m1_hm_rm", rm);
  hm_rkji = mk("m1_hm_rkji", rkji);
  hm_rseg = mk("m1_hm_rseg", rseg);
  std::vector<int> segs(4*std::max(nseg, 1), 0);
  for (int s = 0; s < nseg; ++s) {
    segs[4*s] = hm_soff[s];
    segs[4*s + 1] = hm_slen[s];
    segs[4*s + 2] = hm_roff[s];
    segs[4*s + 3] = hm_rlen[s];
  }
  hm_segs = mk("m1_hm_segs", segs);
  hm_sbuf = DvceArray1D<Real>("m1_hm_sbuf", std::max(hm_nsend, 1)*hm_nq);
  hm_rbuf = DvceArray1D<Real>("m1_hm_rbuf", std::max(hm_nrecv, 1)*hm_nq);
  MPI_Comm *c = new MPI_Comm;
  MPI_Comm_dup(MPI_COMM_WORLD, c);
  hm_comm = static_cast<void *>(c);
  hm_state = 1;
  for (int f = 0; f < 6; ++f) {hm_face[f] = hm_mpif[f];}
  if (impl_halo_ipc) {ImplicitHaloIPCInit();}
  if (me == 0) {
    std::cout << "<rad_m1> implicit_halo_mpi: ON, rank 0 has " << nseg
              << " neighbour ranks, " << hm_nsend << " cells sent per component"
              << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloMPI
//! \brief implicit_halo_mpi: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list).  Receives posted first, then the pack kernel,
//! the on-rank copy kernel queued behind it, the sends as soon as the pack is done, the
//! unpack kernel once every message is in.  Buffer layout per neighbour rank: component
//! major, i.e. entry e of component n at off*nq + n*len + (e - off).

void RadiationM1::ImplicitHaloMPI(int nq, int c0) {
  ImplicitHaloMPIPost(nq, c0);
  ImplicitHaloMPIFinish(nq, c0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloMPIPost
//! \brief implicit_halo_mpi, first half: receives posted, pack kernel, event, on-rank
//! copy kernel.  The host does not wait here, so the caller may queue more work
//! (implicit_halo_overlap: the interior operator) before ImplicitHaloMPIFinish.

void RadiationM1::ImplicitHaloMPIPost(int nq, int c0) {
#if MPI_PARALLEL_ENABLED
  if (hi_state == 1) {ImplicitHaloIPCPost(nq, c0); return;}
  MPI_Comm comm = *static_cast<MPI_Comm *>(hm_comm);
  const int nseg = static_cast<int>(hm_rank.size());
  auto &rq = M1HmReqR();
  auto &sq = M1HmReqS();
  rq.assign(nseg, MPI_REQUEST_NULL);
  sq.assign(nseg, MPI_REQUEST_NULL);
  auto sbuf = hm_sbuf;
  auto rbuf = hm_rbuf;
  constexpr int tag = 4242;
  // hrdet-1009: the receive buffer is still being read by the PREVIOUS exchange's unpack
  // kernel until that kernel has run; a receive posted now lets a fast neighbour's
  // message (GPU-aware MPI, IPC) overwrite it first, so that unpack hands the ghosts the
  // NEXT component.  Back-to-back exchanges with no host sync in between (the three
  // single-component halo exchanges of the half-range lateral faces) hit this:
  // run-to-run differences from cycle 1 on the 16-block BSG layout.
  if (M1HmUnpackPending()) {
    M1EvtU().wait();
    M1HmUnpackPending() = false;
  }
  for (int s = 0; s < nseg; ++s) {
    if (hm_rlen[s] > 0) {
      MPI_Irecv(rbuf.data() + static_cast<size_t>(hm_roff[s])*nq, hm_rlen[s]*nq,
                MPI_ATHENA_REAL, hm_rank[s], tag, comm, &rq[s]);
    }
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nc0 = c0;
  auto iw_ = iw;
  auto segs = hm_segs.d_view;
  if (hm_nsend > 0) {
    auto em = hm_sm.d_view;
    auto ek = hm_skji.d_view;
    auto es = hm_sseg.d_view;
    const int ns = hm_nsend;
    par_for("m1_impl_hm_pack", DevExeSpace(), 0, nq-1, 0, ns-1,
    KOKKOS_LAMBDA(const int n, const int e) {
      const int sg = es(e);
      const int off = segs(4*sg), len = segs(4*sg + 1);
      const int kji = ek(e);
      const int k = kji/(n1*n2);
      const int j = (kji - k*n1*n2)/n1;
      const int i = kji - (k*n2 + j)*n1;
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      sbuf(off*nq + n*len + (e - off)) = iw_(em(e),nc,k,j,i);
    });
  }
  M1EvtB().record();
  ImplicitHaloDirect(nq, c0);   // the on-rank neighbours, behind the pack
#else
  (void)nq;
  (void)c0;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloMPIFinish
//! \brief implicit_halo_mpi, second half: wait for the pack, sends, wait for every
//! message, unpack kernel, wait for the sends.

void RadiationM1::ImplicitHaloMPIFinish(int nq, int c0) {
#if MPI_PARALLEL_ENABLED
  if (hi_state == 1) {ImplicitHaloIPCFinish(nq, c0); return;}
  MPI_Comm comm = *static_cast<MPI_Comm *>(hm_comm);
  const int nseg = static_cast<int>(hm_rank.size());
  auto &rq = M1HmReqR();
  auto &sq = M1HmReqS();
  auto sbuf = hm_sbuf;
  auto rbuf = hm_rbuf;
  constexpr int tag = 4242;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nc0 = c0;
  auto iw_ = iw;
  auto segs = hm_segs.d_view;
  M1EvtB().wait();   // the pack kernel is done: the send buffer may go out
  for (int s = 0; s < nseg; ++s) {
    if (hm_slen[s] > 0) {
      MPI_Isend(sbuf.data() + static_cast<size_t>(hm_soff[s])*nq, hm_slen[s]*nq,
                MPI_ATHENA_REAL, hm_rank[s], tag, comm, &sq[s]);
    }
  }
  MPI_Waitall(nseg, rq.data(), MPI_STATUSES_IGNORE);
  if (hm_nrecv > 0) {
    auto em = hm_rm.d_view;
    auto ek = hm_rkji.d_view;
    auto es = hm_rseg.d_view;
    const int nr = hm_nrecv;
    par_for("m1_impl_hm_unpack", DevExeSpace(), 0, nq-1, 0, nr-1,
    KOKKOS_LAMBDA(const int n, const int e) {
      const int sg = es(e);
      const int off = segs(4*sg + 2), len = segs(4*sg + 3);
      const int kji = ek(e);
      const int k = kji/(n1*n2);
      const int j = (kji - k*n1*n2)/n1;
      const int i = kji - (k*n2 + j)*n1;
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      iw_(em(e),nc,k,j,i) = rbuf(off*nq + n*len + (e - off));
    });
    M1EvtU().record();
    M1HmUnpackPending() = true;
  }
  MPI_Waitall(nseg, sq.data(), MPI_STATUSES_IGNORE);
#else
  (void)nq;
  (void)c0;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloOp
//! \brief the halo of x followed by y = A x (+ the reductions `red` of
//! ImplicitOffDiagOpC into out[0..3]).  Off (the default), exactly ImplicitKrylovHalo
//! then ImplicitOpX.  implicit_halo_overlap (with implicit_halo_mpi active and the
//! stencil operator): the exchange is split around the operator on the interior cells,
//!   receives + pack + on-rank copy | interior operator | [host: wait for the pack,
//!   sends, wait for the messages] unpack | shell operator,
//! so the GPU computes the interior while the messages fly.  Every cell gets the same
//! arithmetic as ImplicitStencilOp, so y is bitwise; with red > 0 the sums are summed
//! in two parts (interior, shell), i.e. round-off.

void RadiationM1::ImplicitHaloOp(int xc, int yc, int red, Real *out) {
  bool ovl = impl_halo_ovl && impl_stencil && !halo_direct_on && impl_halo_mpi;
  if (ovl && hm_state == 0) {ImplicitHaloMPIInit();}
  ovl = ovl && (hm_state == 1);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int w = vimp_now ? 2 : 1;       // the reach of the operator (M1VimpRow: 2)
  if (ovl) {
    // an interior box that is not empty in every non-degenerate direction
    ovl = (indcs.nx1 > 2*w) && ((indcs.nx2 == 1) || (indcs.nx2 > 2*w)) &&
          ((indcs.nx3 == 1) || (indcs.nx3 > 2*w));
  }
  if (!ovl) {
    ImplicitKrylovHalo(xc);
    ImplicitOpX(xc, yc, red, out);
    return;
  }
  if (ho_h.extent_int(0) != 8) {
    ho_h = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_ho_h", 8);
  }
  Real *h0 = ho_h.data(), *h1 = ho_h.data() + 4;
  if (impl_opsplit && red != 0) {
    // implicit_op_split_red (m1-fast3): both parts as plain kernels, then ONE read-only
    // reduction over every active cell (ImplicitStencilOp, red_only): the sums are those
    // of the one-rank-style flat reduction, not interior + shell (round-off).
    // implicit_op_team_red alone keeps the fused interior + shell parts below (bitwise
    // the default at 2+ ranks): the split path measured 2-3 % slower there (runs_5l)
    ImplicitHaloMPIPost(1, xc);
    ImplicitStencilOpPart(xc, yc, 0, 1, w, h0);
    ImplicitHaloMPIFinish(1, xc);
    ImplicitStencilOpPart(xc, yc, 0, 2, w, h1);
    ImplicitStencilOp(xc, yc, red, out, true);
    return;
  }
  ImplicitHaloMPIPost(1, xc);
  ImplicitStencilOpPart(xc, yc, red, 1, w, h0);
  ImplicitHaloMPIFinish(1, xc);
  ImplicitStencilOpPart(xc, yc, red, 2, w, h1);
  if (red != 0) {
    DevExeSpace().fence();
    out[0] = h0[0] + h1[0];
    out[1] = h0[1] + h1[1];
    out[2] = h0[2] + h1[2];
    out[3] = (h1[3] > h0[3]) ? h1[3] : h0[3];
  }
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabPipe
//! \brief implicit_krylov_pipe: preconditioned p-BiCGStab (Cools & Vanroose 2017,
//! Alg. 4), right preconditioning x = M^-1 u as in ImplicitBiCGStabTwo.  Vectors
//! (iw slots / kpw slots):
//!   x KX, r KR, rtilde KRH (shadow), b KB,
//!   w = A rhat KP, what = M^-1 w KY, t = A what KTT,
//!   z = A shat KS, zhat = M^-1 z KZ, v = A zhat KV,
//!   phat, s = A phat, shat, rhat = M^-1 r in kpw.
//! Iteration i (alpha_i, beta_{i-1}, omega_{i-1} known):
//!   K1: phat = rhat + beta (phat - omega shat); s = w + beta (s - omega z);
//!       shat = what + beta (shat - omega zhat); z = t + beta (z - omega v);
//!       q = r - alpha s, y = w - alpha z; reduce (q,y), (y,y)          [R1]
//!   P1: zhat = M^-1 z; halo; v = A zhat                    (R1 in flight)
//!   omega = (q,y)/(y,y)
//!   K2: qhat = rhat - alpha shat; x += alpha phat + omega qhat; r = q - omega y;
//!       rhat = qhat - omega (what - alpha zhat); w = y - omega (t - alpha v);
//!       reduce (rt,r), (rt,w), (rt,s), (rt,z), max|r|                  [R2]
//!   P2: what = M^-1 w; halo; t = A what                    (R2 in flight)
//!   beta = (alpha/omega) (rt,r)/rho, alpha = (rt,r)/((rt,w) + beta (rt,s)
//!          - beta omega (rt,z)), rho = (rt,r).
//! The convergence test on max|r_{i+1}| is one P2 late (P2 is wasted on the last
//! iteration) and is confirmed on the true residual; a failed confirmation, and every
//! breakdown, restarts from the true residual of x (restart and fallback rules of the
//! fused loop).

int RadiationM1::ImplicitBiCGStabPipe(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nmb1 = nmb - 1;
  if (kpw.extent_int(0) != nmb) {
    kpw = DvceArray5D<Real>("m1_kpw", nmb, 4, iw.extent(2), iw.extent(3), iw.extent(4));
  }
  if (kp_h.extent_int(0) != 12) {
    kp_h = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_kp_h", 12);
  }
  auto iw_ = iw;
  auto kp_ = kpw;
  constexpr int PH = 0, SV = 1, SH = 2, RHT = 3;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, nmb*nkji);
  using HRed = M1PipeRed<Kokkos::HostSpace>;
  using PRed = M1PipeRed<Kokkos::SharedHostPinnedSpace>;
  M1PipeVal *h1 = reinterpret_cast<M1PipeVal *>(kp_h.data());
  M1PipeVal *h2 = reinterpret_cast<M1PipeVal *>(kp_h.data() + 6);
  typename PRed::result_view_type v1(h1), v2(h2);
  auto blocking_sum = [&](M1PipeVal &v) {
    M1PipePost(v);
    M1PipeFinish(v);
    bcg_nred += 1.0;
  };
  // r = b - A x (true residual) into KR (and into KRH with `shadow`), (r,r), max|r|
  auto true_res = [&](bool shadow) -> M1PipeVal {
    ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
    M1PipeVal tr;
    const bool sh = shadow;
    Kokkos::parallel_reduce("m1_impl_pipe_res", pol,
    KOKKOS_LAMBDA(const int idx, M1PipeVal &v) {
      int m, k, j, i;
      M1PipeIdx(idx, nkji, nji, ni, m, k, j, i);
      k += ks; j += js; i += is;
      Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = r;
      if (sh) {iw_(m,M1_IW_KRH,k,j,i) = r;}
      v.s[0] += r*r;
      Real a = fabs(r);
      v.mx = (a > v.mx) ? a : v.mx;
    }, HRed(tr));
    blocking_sum(tr);
    return tr;
  };

  // x0 = the Picard iterate
  par_for("m1_impl_pipe_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  M1PipeVal r0 = true_res(true);
  Real rnorm = r0.mx;
  Real rho = r0.s[0];
  bcg_r0rel = rnorm/bscale;
  const bool ew = (impl_ew_max > 0.0) && !ew_tight;   // ew_tight: runs_4a_accel
  Real tabs = tol*bscale;
  if (ew) {
    Real eta = impl_ew_max;
    if (ew_fprev > 0.0) {
      eta = impl_ew_gam*SQR(rnorm/ew_fprev);
      const Real sg = impl_ew_gam*SQR(ew_etaprev);
      if (sg > 0.1) {eta = fmax(eta, sg);}
      eta = fmin(eta, impl_ew_max);
    }
    ew_fprev = rnorm;
    ew_etaprev = eta;
    tabs = fmax(tabs, eta*rnorm);
  }
  auto lin_done = [=](const Real r) {
    return ew ? (r < tabs) : (r/bscale < tol);
  };
  // the start-up from r (KR) and rtilde (KRH): rhat = M^-1 r, w = A rhat, what = M^-1 w,
  // t = A what; returns (rtilde, w)
  auto start = [&]() -> Real {
    ImplicitPrecondX(M1_IW_KR, M1_IW_KY, 0, 0.0, 0.0);
    Real o4[4];
    ImplicitHaloOp(M1_IW_KY, M1_IW_KP, 1, o4);
    M1PipeVal a;
    for (int c = 0; c < 5; ++c) {a.s[c] = 0.0;}
    a.s[0] = o4[0];
    a.mx = 0.0;
    blocking_sum(a);
    par_for("m1_impl_pipe_st", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      kp_(m,RHT,k,j,i) = iw_(m,M1_IW_KY,k,j,i);
      kp_(m,PH,k,j,i) = 0.0;
      kp_(m,SV,k,j,i) = 0.0;
      kp_(m,SH,k,j,i) = 0.0;
      iw_(m,M1_IW_KS,k,j,i) = 0.0;
      iw_(m,M1_IW_KZ,k,j,i) = 0.0;
      iw_(m,M1_IW_KV,k,j,i) = 0.0;
    });
    ImplicitPrecondX(M1_IW_KP, M1_IW_KY, 0, 0.0, 0.0);
    ImplicitHaloOp(M1_IW_KY, M1_IW_KTT, 0, nullptr);
    return a.s[0];
  };

  int nit = 0;
  int nrestart = 0;
  Real alpha = 1.0, beta = 0.0, omega = 1.0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  bool need_start = true;
  while (!done && nit < impl_lin_maxit) {
    bool breakdown = false;
    if (need_start) {
      need_start = false;
      const Real rw = start();
      if (!(fabs(rho) > M1_BCG_EPS) || !(fabs(rw) > M1_BCG_EPS)) {
        breakdown = true;
      } else {
        alpha = rho/rw;
        beta = 0.0;
        omega = 1.0;
      }
    }
    if (!breakdown) {
      ++nit;
      const Real al = alpha, be = beta, om = omega;
      // K1 + R1
      Kokkos::parallel_reduce("m1_impl_pipe_k1", pol,
      KOKKOS_LAMBDA(const int idx, M1PipeVal &v) {
        int m, k, j, i;
        M1PipeIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        const Real w = iw_(m,M1_IW_KP,k,j,i);
        const Real z = iw_(m,M1_IW_KTT,k,j,i)
                       + be*(iw_(m,M1_IW_KS,k,j,i) - om*iw_(m,M1_IW_KV,k,j,i));
        const Real s = w + be*(kp_(m,SV,k,j,i) - om*iw_(m,M1_IW_KS,k,j,i));
        const Real shold = kp_(m,SH,k,j,i);
        kp_(m,PH,k,j,i) = kp_(m,RHT,k,j,i) + be*(kp_(m,PH,k,j,i) - om*shold);
        kp_(m,SH,k,j,i) = iw_(m,M1_IW_KY,k,j,i)
                          + be*(shold - om*iw_(m,M1_IW_KZ,k,j,i));
        kp_(m,SV,k,j,i) = s;
        iw_(m,M1_IW_KS,k,j,i) = z;
        const Real q = iw_(m,M1_IW_KR,k,j,i) - al*s;
        const Real y = w - al*z;
        v.s[0] += q*y;
        v.s[1] += y*y;
      }, PRed(v1));
      M1EvtA().record();
      // P1: zhat = M^-1 z, halo, v = A zhat, with R1 in flight
      ImplicitPrecondX(M1_IW_KS, M1_IW_KZ, 0, 0.0, 0.0);
      M1EvtA().wait();
      M1PipeVal a = *h1;
      M1PipePost(a);
      ImplicitHaloOp(M1_IW_KZ, M1_IW_KV, 0, nullptr);
      M1PipeFinish(a);
      bcg_nred += 1.0;
      omega = (a.s[1] > 0.0) ? (a.s[0]/a.s[1]) : 0.0;
      if (!(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}   // x, r untouched so far
    }
    if (!breakdown) {
      const Real al = alpha, om = omega;
      // K2 + R2
      Kokkos::parallel_reduce("m1_impl_pipe_k2", pol,
      KOKKOS_LAMBDA(const int idx, M1PipeVal &v) {
        int m, k, j, i;
        M1PipeIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        const Real s = kp_(m,SV,k,j,i);
        const Real z = iw_(m,M1_IW_KS,k,j,i);
        const Real w = iw_(m,M1_IW_KP,k,j,i);
        const Real q = iw_(m,M1_IW_KR,k,j,i) - al*s;
        const Real y = w - al*z;
        const Real qh = kp_(m,RHT,k,j,i) - al*kp_(m,SH,k,j,i);
        iw_(m,M1_IW_KX,k,j,i) += al*kp_(m,PH,k,j,i) + om*qh;
        const Real r = q - om*y;
        iw_(m,M1_IW_KR,k,j,i) = r;
        kp_(m,RHT,k,j,i) = qh - om*(iw_(m,M1_IW_KY,k,j,i) - al*iw_(m,M1_IW_KZ,k,j,i));
        const Real wn = y - om*(iw_(m,M1_IW_KTT,k,j,i) - al*iw_(m,M1_IW_KV,k,j,i));
        iw_(m,M1_IW_KP,k,j,i) = wn;
        const Real rt = iw_(m,M1_IW_KRH,k,j,i);
        v.s[0] += rt*r;
        v.s[1] += rt*wn;
        v.s[2] += rt*s;
        v.s[3] += rt*z;
        const Real ar = fabs(r);
        v.mx = (ar > v.mx) ? ar : v.mx;
      }, PRed(v2));
      M1EvtA().record();
      // P2: what = M^-1 w, halo, t = A what, with R2 in flight
      ImplicitPrecondX(M1_IW_KP, M1_IW_KY, 0, 0.0, 0.0);
      M1EvtA().wait();
      M1PipeVal b = *h2;
      M1PipePost(b);
      ImplicitHaloOp(M1_IW_KY, M1_IW_KTT, 0, nullptr);
      M1PipeFinish(b);
      bcg_nred += 1.0;
      if (lin_done(b.mx)) {
        const M1PipeVal tr = true_res(false);
        if (lin_done(tr.mx)) {
          done = true;
          break;
        }
        breakdown = true;   // the recursive residual drifted: restart from the true one
      } else {
        const Real rn = b.s[0];
        const Real bt = (alpha/omega)*(rn/rho);
        const Real den = b.s[1] + bt*b.s[2] - bt*omega*b.s[3];
        if (!(fabs(rn) > M1_BCG_EPS) || !(fabs(den) > M1_BCG_EPS)) {
          breakdown = true;
        } else {
          beta = bt;
          alpha = rn/den;
          rho = rn;
        }
      }
    }
    if (breakdown && !done) {
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        fell_back = true;
        break;
      }
      const M1PipeVal tr = true_res(true);
      rho = tr.s[0];
      if (lin_done(tr.mx)) {
        done = true;
        break;
      }
      need_start = true;
    }
  }
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
// implicit_halo_ipc (m1-perf-0928): node-local halo exchange over CUDA IPC.
//
// Buffers (one per rank, raw cudaMalloc so the IPC handle's base IS the pointer):
//   hi_rbuf  2 x hi_half Reals, parity p = seq & 1; segment s of component n of a
//            call with nq components at p*hi_half + roff[s]*nq + n*rlen[s] + e (the
//            layout of the MPI path's hm_rbuf, so the unpack is the same arithmetic)
//   hi_sig   64 long long: [s] = the last exchange whose data from segment s arrived,
//            [32 + s] = the last exchange of MY data that segment s has unpacked
// Protocol of exchange seq (every rank runs the same sequence of exchanges):
//   Post:   wait until every neighbour acked seq-2 (it is done with this parity's
//           buffer), pack straight into the neighbour's buffer, system fence, set the
//           neighbour's arrival flag [pslot] = seq
//   Finish: spin until every arrival flag >= seq, unpack, ack into the sender's
//           [32 + pslot] = seq
// Every step is a kernel on the one device stream, so the host never waits.  A spin
// gives up after ~20 s of GPU clock and raises hi_err; the host checks it every 4096
// exchanges and at the end of the run's solves, and aborts with a message.
namespace {
constexpr int M1_HI_NSIG = 32;
constexpr long long M1_HI_SPIN = 40000000000LL;   // GPU cycles (~20 s at ~2 GHz)

KOKKOS_INLINE_FUNCTION
void M1HiWaitGE(const long long *f, const long long target, int *err) {
#if defined(__CUDA_ARCH__)
  const volatile long long *vf = f;
  const long long t0 = clock64();
  while (*vf < target) {
    __nanosleep(256);
    if (clock64() - t0 > M1_HI_SPIN) {*err = 1; return;}
  }
  __threadfence_system();
#else
  (void)f; (void)target; (void)err;
#endif
}

KOKKOS_INLINE_FUNCTION
void M1HiSet(long long *f, const long long v) {
#if defined(__CUDA_ARCH__)
  __threadfence_system();
  volatile long long *vf = f;
  *vf = v;
  __threadfence_system();
#else
  (void)f; (void)v;
#endif
}
} // namespace

void RadiationM1::ImplicitHaloIPCInit() {
  hi_state = -1;
#if MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
  const int me = global_variable::my_rank;
  const int nseg = static_cast<int>(hm_rank.size());
  MPI_Comm comm = *static_cast<MPI_Comm *>(hm_comm);
  // every neighbour rank on this node?
  int ok = (nseg <= M1_HI_NSIG) ? 1 : 0;
  MPI_Comm node;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, me, MPI_INFO_NULL, &node);
  int nn;
  MPI_Comm_size(node, &nn);
  std::vector<int> nodeg(nn);
  MPI_Allgather(&me, 1, MPI_INT, nodeg.data(), 1, MPI_INT, node);
  MPI_Comm_free(&node);
  for (int s = 0; s < nseg; ++s) {
    if (std::find(nodeg.begin(), nodeg.end(), hm_rank[s]) == nodeg.end()) {ok = 0;}
  }
  hi_half = static_cast<size_t>(std::max(hm_nrecv, 1))*hm_nq;
  if (ok) {
    if (cudaMalloc(reinterpret_cast<void **>(&hi_rbuf), 2*hi_half*sizeof(Real))
        != cudaSuccess) {ok = 0;}
    if (ok && cudaMalloc(reinterpret_cast<void **>(&hi_sig),
                         2*M1_HI_NSIG*sizeof(long long)) != cudaSuccess) {ok = 0;}
    if (ok && cudaMalloc(reinterpret_cast<void **>(&hi_err), sizeof(int))
        != cudaSuccess) {ok = 0;}
    if (ok) {
      (void)cudaMemset(hi_sig, 0, 2*M1_HI_NSIG*sizeof(long long));
      (void)cudaMemset(hi_err, 0, sizeof(int));
      (void)cudaDeviceSynchronize();
    }
  }
  {int g = ok;
  MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  ok = g;}
  if (ok) {
    // swap with every neighbour: IPC handles, my segment index for it, my receive
    // offset/length for its data, my parity-half size
    cudaIpcMemHandle_t hb, hs;
    (void)cudaIpcGetMemHandle(&hb, hi_rbuf);
    (void)cudaIpcGetMemHandle(&hs, hi_sig);
    struct Info {cudaIpcMemHandle_t b, s; long long half; int slot, roff, rlen;};
    std::vector<Info> mine(nseg), theirs(nseg);
    std::vector<MPI_Request> rq(2*nseg);
    for (int s = 0; s < nseg; ++s) {
      mine[s] = Info{hb, hs, static_cast<long long>(hi_half), s, hm_roff[s], hm_rlen[s]};
      MPI_Irecv(&theirs[s], sizeof(Info), MPI_BYTE, hm_rank[s], 4243, comm, &rq[s]);
    }
    for (int s = 0; s < nseg; ++s) {
      MPI_Isend(&mine[s], sizeof(Info), MPI_BYTE, hm_rank[s], 4243, comm, &rq[nseg + s]);
    }
    MPI_Waitall(2*nseg, rq.data(), MPI_STATUSES_IGNORE);
    hi_prb = DualArray1D<unsigned long long>("m1_hi_prb", std::max(nseg, 1));
    hi_psig = DualArray1D<unsigned long long>("m1_hi_psig", std::max(nseg, 1));
    hi_phalf = DualArray1D<long long>("m1_hi_phalf", std::max(nseg, 1));
    hi_proff = DualArray1D<int>("m1_hi_proff", std::max(nseg, 1));
    hi_pslot = DualArray1D<int>("m1_hi_pslot", std::max(nseg, 1));
    for (int s = 0; s < nseg; ++s) {
      void *pb = nullptr, *ps = nullptr;
      if (theirs[s].rlen != hm_slen[s]) {ok = 0;}
      if (cudaIpcOpenMemHandle(&pb, theirs[s].b, cudaIpcMemLazyEnablePeerAccess)
          != cudaSuccess) {ok = 0;}
      if (cudaIpcOpenMemHandle(&ps, theirs[s].s, cudaIpcMemLazyEnablePeerAccess)
          != cudaSuccess) {ok = 0;}
      hi_prb.h_view(s) = reinterpret_cast<unsigned long long>(pb);
      hi_psig.h_view(s) = reinterpret_cast<unsigned long long>(ps);
      hi_phalf.h_view(s) = theirs[s].half;
      hi_proff.h_view(s) = theirs[s].roff;
      hi_pslot.h_view(s) = theirs[s].slot;
    }
    (void)cudaGetLastError();
    for (auto *a : {&hi_prb, &hi_psig}) {a->modify_host(); a->sync_device();}
    hi_phalf.modify_host(); hi_phalf.sync_device();
    hi_proff.modify_host(); hi_proff.sync_device();
    hi_pslot.modify_host(); hi_pslot.sync_device();
    int g = ok;
    MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    ok = g;
    MPI_Barrier(MPI_COMM_WORLD);   // every peer's flags are zero before any exchange
  }
  hi_state = ok ? 1 : -1;
  hi_seq = 0;
  hi_nchk = 0;
  if (me == 0) {
    std::cout << "<rad_m1> implicit_halo_ipc: "
              << (ok ? "ON (node-local CUDA IPC, no MPI in the implicit halo)"
                     : "not possible (neighbour off-node, or IPC failed): MPI path")
              << std::endl;
  }
#endif
}

void RadiationM1::ImplicitHaloIPCPost(int nq, int c0) {
#if MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
  const int nseg = static_cast<int>(hm_rank.size());
  const long long seq = ++hi_seq;
  const long long par = seq & 1;
  long long *sig = hi_sig;
  int *err = hi_err;
  auto psig = hi_psig.d_view;
  auto pslot = hi_pslot.d_view;
  // the neighbours are done with this parity's buffers (exchange seq-2)
  Kokkos::parallel_for("m1_hi_wack", Kokkos::RangePolicy<DevExeSpace>(0, 1),
  KOKKOS_LAMBDA(const int) {
    for (int s = 0; s < nseg; ++s) {M1HiWaitGE(sig + M1_HI_NSIG + s, seq - 2, err);}
  });
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nc0 = c0;
  auto iw_ = iw;
  auto segs = hm_segs.d_view;
  if (hm_nsend > 0) {
    auto em = hm_sm.d_view;
    auto ek = hm_skji.d_view;
    auto es = hm_sseg.d_view;
    auto prb = hi_prb.d_view;
    auto phalf = hi_phalf.d_view;
    auto proff = hi_proff.d_view;
    const int ns = hm_nsend;
    par_for("m1_impl_hm_pack", DevExeSpace(), 0, nq-1, 0, ns-1,
    KOKKOS_LAMBDA(const int n, const int e) {
      const int sg = es(e);
      const int off = segs(4*sg), len = segs(4*sg + 1);
      const int kji = ek(e);
      const int k = kji/(n1*n2);
      const int j = (kji - k*n1*n2)/n1;
      const int i = kji - (k*n2 + j)*n1;
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      Real *dst = reinterpret_cast<Real *>(prb(sg)) + par*phalf(sg)
                  + static_cast<long long>(proff(sg))*nq + n*len + (e - off);
      *dst = iw_(em(e),nc,k,j,i);
    });
  }
  // arrival flags of this exchange in every neighbour
  Kokkos::parallel_for("m1_hi_flag", Kokkos::RangePolicy<DevExeSpace>(0, 1),
  KOKKOS_LAMBDA(const int) {
    for (int s = 0; s < nseg; ++s) {
      M1HiSet(reinterpret_cast<long long *>(psig(s)) + pslot(s), seq);
    }
  });
  ImplicitHaloDirect(nq, c0);   // the on-rank neighbours
#else
  (void)nq;
  (void)c0;
#endif
}

void RadiationM1::ImplicitHaloIPCFinish(int nq, int c0) {
#if MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
  const int nseg = static_cast<int>(hm_rank.size());
  const long long seq = hi_seq;
  const long long par = seq & 1;
  long long *sig = hi_sig;
  int *err = hi_err;
  auto psig = hi_psig.d_view;
  auto pslot = hi_pslot.d_view;
  Kokkos::parallel_for("m1_hi_warr", Kokkos::RangePolicy<DevExeSpace>(0, 1),
  KOKKOS_LAMBDA(const int) {
    for (int s = 0; s < nseg; ++s) {M1HiWaitGE(sig + s, seq, err);}
  });
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nc0 = c0;
  auto iw_ = iw;
  auto segs = hm_segs.d_view;
  if (hm_nrecv > 0) {
    auto em = hm_rm.d_view;
    auto ek = hm_rkji.d_view;
    auto es = hm_rseg.d_view;
    const Real *rb = hi_rbuf + par*static_cast<long long>(hi_half);
    const int nr = hm_nrecv;
    par_for("m1_impl_hm_unpack", DevExeSpace(), 0, nq-1, 0, nr-1,
    KOKKOS_LAMBDA(const int n, const int e) {
      const int sg = es(e);
      const int off = segs(4*sg + 2), len = segs(4*sg + 3);
      const int kji = ek(e);
      const int k = kji/(n1*n2);
      const int j = (kji - k*n1*n2)/n1;
      const int i = kji - (k*n2 + j)*n1;
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      iw_(em(e),nc,k,j,i) = rb[static_cast<long long>(off)*nq + n*len + (e - off)];
    });
  }
  // tell every sender its buffer half is free again
  Kokkos::parallel_for("m1_hi_ack", Kokkos::RangePolicy<DevExeSpace>(0, 1),
  KOKKOS_LAMBDA(const int) {
    for (int s = 0; s < nseg; ++s) {
      M1HiSet(reinterpret_cast<long long *>(psig(s)) + M1_HI_NSIG + pslot(s), seq);
    }
  });
  if (++hi_nchk >= 4096) {
    hi_nchk = 0;
    int h = 0;
    (void)cudaMemcpy(&h, hi_err, sizeof(int), cudaMemcpyDeviceToHost);
    if (h != 0) {
      std::cout << "### FATAL ERROR in rad_m1 implicit_halo_ipc: a halo spin timed out "
                << "(rank " << global_variable::my_rank << ", exchange " << seq << ")"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
#else
  (void)nq;
  (void)c0;
#endif
}

} // namespace radm1
