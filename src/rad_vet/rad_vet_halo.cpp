//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_halo.cpp
//! \brief the ghost exchange of the work array on several ranks: implicit_halo_mpi (the
//! on-rank neighbours by the implicit_halo_direct copy kernel, the off-rank ones by one
//! message per neighbour rank; bitwise the ordinary exchange) and implicit_halo_ipc
//! (node-local, over CUDA IPC).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#if defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace radvet {

namespace {
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

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitHaloMPIInit
//! \brief implicit_halo_mpi: the cell lists of the off-rank part of the implicit
//! exchanges.  Receiver ghost region (G, o) of MeshBlock G towards direction o is filled
//! from the active cells of the same-level neighbour g, cell (k,j,i) <- g's cell
//! (k - o3 nx3, j - o2 nx2, i - o1 nx1) (the map of ImplicitHaloDirect).  On each rank
//! pair both sides enumerate the regions sorted by (receiver gid, o) and each region in
//! (k,j,i) order, so the message needs no header.  hm_state = -1 (the ordinary
//! exchange keeps running) unless EVERY rank has same-level neighbours only, outside the
//! cubed-sphere and polar transforms.

void RadiationVet::ImplicitHaloMPIInit() {
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
      std::cout << "<rad_vet> implicit_halo_mpi: mesh not supported (levels / seams), "
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
  // G5: the ghost depth of the lists = <rad_vet>/implicit_halo_depth where it
  // applies (the spherical-polar wedge: the same layers the implicit kernels read,
  // poisoned-ghost gate), else all ng layers as before
  const int dg = (impl_hdepth > 0 && impl_hdepth < ng) ? impl_hdepth : ng;
  auto rng = [&](int o, int s, int e, int nx, bool on, int &lo, int &hi) {
    if (!on) {lo = s; hi = e; return;}
    (void)nx;
    if (o < 0) {
      lo = s - dg; hi = s - 1;
    } else if (o > 0) {
      lo = e + 1; hi = e + dg;
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
  hm_sm = mk("vet_hm_sm", sm);
  hm_skji = mk("vet_hm_skji", skji);
  hm_sseg = mk("vet_hm_sseg", sseg);
  hm_rm = mk("vet_hm_rm", rm);
  hm_rkji = mk("vet_hm_rkji", rkji);
  hm_rseg = mk("vet_hm_rseg", rseg);
  std::vector<int> segs(4*std::max(nseg, 1), 0);
  for (int s = 0; s < nseg; ++s) {
    segs[4*s] = hm_soff[s];
    segs[4*s + 1] = hm_slen[s];
    segs[4*s + 2] = hm_roff[s];
    segs[4*s + 3] = hm_rlen[s];
  }
  hm_segs = mk("vet_hm_segs", segs);
  hm_sbuf = DvceArray1D<Real>("vet_hm_sbuf", std::max(hm_nsend, 1)*hm_nq);
  hm_rbuf = DvceArray1D<Real>("vet_hm_rbuf", std::max(hm_nrecv, 1)*hm_nq);
#if defined(KOKKOS_ENABLE_HIP)
  hm_shost = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("vet_hm_shost",
                                                               hm_sbuf.extent(0));
  hm_rhost = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("vet_hm_rhost",
                                                               hm_rbuf.extent(0));
#endif
  MPI_Comm *c = new MPI_Comm;
  MPI_Comm_dup(MPI_COMM_WORLD, c);
  hm_comm = static_cast<void *>(c);
  hm_state = 1;
  for (int f = 0; f < 6; ++f) {hm_face[f] = hm_mpif[f];}
  if (impl_halo_ipc) {ImplicitHaloIPCInit();}
  if (me == 0) {
    std::cout << "<rad_vet> implicit_halo_mpi: ON, rank 0 has " << nseg
              << " neighbour ranks, " << hm_nsend << " cells sent per component"
              << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitHaloMPI
//! \brief implicit_halo_mpi: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list).  Receives posted first, then the pack kernel,
//! the on-rank copy kernel queued behind it, the sends as soon as the pack is done, the
//! unpack kernel once every message is in.  Buffer layout per neighbour rank: component
//! major, i.e. entry e of component n at off*nq + n*len + (e - off).

void RadiationVet::ImplicitHaloMPI(int nq, int c0) {
  ImplicitHaloMPIPost(nq, c0);
  ImplicitHaloMPIFinish(nq, c0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitHaloMPIPost
//! \brief implicit_halo_mpi, first half: receives posted, pack kernel, event, on-rank
//! copy kernel.  The host does not wait here, so the caller may queue more work
//! (implicit_halo_overlap: the interior operator) before ImplicitHaloMPIFinish.

void RadiationVet::ImplicitHaloMPIPost(int nq, int c0) {
#if MPI_PARALLEL_ENABLED
  if (hi_state == 1) {ImplicitHaloIPCPost(nq, c0); return;}
  MPI_Comm comm = *static_cast<MPI_Comm *>(hm_comm);
  const int nseg = static_cast<int>(hm_rank.size());
  auto &rq = M1HmReqR();
  auto &sq = M1HmReqS();
  rq.assign(nseg, MPI_REQUEST_NULL);
  sq.assign(nseg, MPI_REQUEST_NULL);
#if defined(KOKKOS_ENABLE_HIP)
  // the pack kernel writes the pinned host send buffer itself
  auto sbuf = hm_shost;
#else
  auto sbuf = hm_sbuf;
  auto rbuf = hm_rbuf;
#endif
  constexpr int tag = 4242;
  // the receive buffer is still being read by the PREVIOUS exchange's unpack
  // kernel until that kernel has run; a receive posted now lets a fast neighbour's
  // message (GPU-aware MPI, IPC) overwrite it first, so that unpack hands the ghosts the
  // NEXT component.  Back-to-back exchanges with no host sync in between (the three
  // single-component halo exchanges of the half-range lateral faces) hit this:
  // run-to-run differences from cycle 1 on the 16-block BSG layout.
  if (M1HmUnpackPending()) {
    M1EvtU().wait();
    M1HmUnpackPending() = false;
  }
  // on HIP the messages go through host copies (ImplicitHaloMPIFinish)
#if defined(KOKKOS_ENABLE_HIP)
  Real *rdst = hm_rhost.data();
#else
  Real *rdst = rbuf.data();
#endif
  for (int s = 0; s < nseg; ++s) {
    if (hm_rlen[s] > 0) {
      MPI_Irecv(rdst + static_cast<size_t>(hm_roff[s])*nq, hm_rlen[s]*nq,
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
    par_for("vet_impl_hm_pack", DevExeSpace(), 0, nq-1, 0, ns-1,
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
//! \fn void RadiationVet::ImplicitHaloMPIFinish
//! \brief implicit_halo_mpi, second half: wait for the pack, sends, wait for every
//! message, unpack kernel, wait for the sends.

void RadiationVet::ImplicitHaloMPIFinish(int nq, int c0) {
#if MPI_PARALLEL_ENABLED
  if (hi_state == 1) {ImplicitHaloIPCFinish(nq, c0); return;}
  MPI_Comm comm = *static_cast<MPI_Comm *>(hm_comm);
  const int nseg = static_cast<int>(hm_rank.size());
  auto &rq = M1HmReqR();
  auto &sq = M1HmReqS();
#if defined(KOKKOS_ENABLE_HIP)
  // the unpack kernel reads the pinned host receive buffer itself
  auto rbuf = hm_rhost;
#else
  auto sbuf = hm_sbuf;
  auto rbuf = hm_rbuf;
#endif
  constexpr int tag = 4242;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nc0 = c0;
  auto iw_ = iw;
  auto segs = hm_segs.d_view;
  M1EvtB().wait();   // the pack kernel is done: the send buffer may go out
  // on HIP (MI300A) the GPU-aware MPI path delivered these device buffers
  // corrupted (942 of 2140 messages of a 2-rank accretor run had a receive checksum !=
  // the sender's, with the send buffer fenced and complete; wrong and run-to-run random
  // ghosts), so MPI is handed host copies there.  CUDA builds keep the device buffers.
#if defined(KOKKOS_ENABLE_HIP)
  // the pack kernel wrote hm_shost (pinned host) and is complete (event
  // above), so the messages go out of it as they did out of the staging copy (was a
  // D2H deep_copy here and an H2D one after the receives: d6fast2-1112 had cut them to
  // the [0, n*nq) prefix); the same values cross, the ghosts are bitwise
  Real *ssrc = hm_shost.data();
#else
  Real *ssrc = sbuf.data();
#endif
  for (int s = 0; s < nseg; ++s) {
    if (hm_slen[s] > 0) {
      MPI_Isend(ssrc + static_cast<size_t>(hm_soff[s])*nq, hm_slen[s]*nq,
                MPI_ATHENA_REAL, hm_rank[s], tag, comm, &sq[s]);
    }
  }
  MPI_Waitall(nseg, rq.data(), MPI_STATUSES_IGNORE);
  if (hm_nrecv > 0) {
    auto em = hm_rm.d_view;
    auto ek = hm_rkji.d_view;
    auto es = hm_rseg.d_view;
    const int nr = hm_nrecv;
    par_for("vet_impl_hm_unpack", DevExeSpace(), 0, nq-1, 0, nr-1,
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

void RadiationVet::ImplicitHaloIPCInit() {
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
    hi_prb = DualArray1D<unsigned long long>("vet_hi_prb", std::max(nseg, 1));
    hi_psig = DualArray1D<unsigned long long>("vet_hi_psig", std::max(nseg, 1));
    hi_phalf = DualArray1D<long long>("vet_hi_phalf", std::max(nseg, 1));
    hi_proff = DualArray1D<int>("vet_hi_proff", std::max(nseg, 1));
    hi_pslot = DualArray1D<int>("vet_hi_pslot", std::max(nseg, 1));
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
    std::cout << "<rad_vet> implicit_halo_ipc: "
              << (ok ? "ON (node-local CUDA IPC, no MPI in the implicit halo)"
                     : "not possible (neighbour off-node, or IPC failed): MPI path")
              << std::endl;
  }
#endif
}

void RadiationVet::ImplicitHaloIPCPost(int nq, int c0) {
#if MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
  const int nseg = static_cast<int>(hm_rank.size());
  const long long seq = ++hi_seq;
  const long long par = seq & 1;
  long long *sig = hi_sig;
  int *err = hi_err;
  auto psig = hi_psig.d_view;
  auto pslot = hi_pslot.d_view;
  // the neighbours are done with this parity's buffers (exchange seq-2)
  Kokkos::parallel_for("vet_hi_wack", Kokkos::RangePolicy<DevExeSpace>(0, 1),
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
    par_for("vet_impl_hm_pack", DevExeSpace(), 0, nq-1, 0, ns-1,
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
  Kokkos::parallel_for("vet_hi_flag", Kokkos::RangePolicy<DevExeSpace>(0, 1),
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

void RadiationVet::ImplicitHaloIPCFinish(int nq, int c0) {
#if MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
  const int nseg = static_cast<int>(hm_rank.size());
  const long long seq = hi_seq;
  const long long par = seq & 1;
  long long *sig = hi_sig;
  int *err = hi_err;
  auto psig = hi_psig.d_view;
  auto pslot = hi_pslot.d_view;
  Kokkos::parallel_for("vet_hi_warr", Kokkos::RangePolicy<DevExeSpace>(0, 1),
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
    par_for("vet_impl_hm_unpack", DevExeSpace(), 0, nq-1, 0, nr-1,
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
  Kokkos::parallel_for("vet_hi_ack", Kokkos::RangePolicy<DevExeSpace>(0, 1),
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
      std::cout << "### FATAL ERROR in rad_vet implicit_halo_ipc: a halo spin timed out "
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

} // namespace radvet
