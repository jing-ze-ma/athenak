#ifndef RAD_VET_RAYS_SC_CART_EXCHANGE_HPP_
#define RAD_VET_RAYS_SC_CART_EXCHANGE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_exchange.hpp
//! \brief Cartesian short characteristics: exchange of bands and planes between
//! MeshBlocks.

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

static inline void VetFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in rad_vet/rays/sc_cart" << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------
//! \struct ScCartMBState
//! \brief state of the several-MeshBlock sweep (ScCartSweepMB): banded plane buffers,
//! neighbour tables and message buffers.  A banded plane is (kk, jj), kk in [0, n3w),
//! jj in [0, n2w), active cells at kk = k - ks + b3, jj = j - js + b2.

struct ScCartMBState {
  int w2 = 0, w3 = 0;          // reach of a ray per layer in x2, x3 (w3 = 0 in 2-D)
  int b2 = 0, b3 = 0;          // band widths = w2, w3 (the ghost band of a plane)
  int nring = 2;               // intensity planes kept
  DvceArray2D<int> rdep;       // (nray, 4) cells of ghost band a ray reads per layer:
                               // below / above in x2, below / above in x3
  DvceArray2D<int> offs, offr;   // (nof, nray+1) per-ray offsets of the plane messages
  std::vector<int> cnts, cntr;   // (nof) plane message sizes (sent / received)
  int n2w = 1, n3w = 1;        // banded plane extents
  int nbx1 = 1, nx1g = 0;      // MeshBlocks along x1, global layers along x1
  int nof = 0;                 // horizontal neighbour slots: 2 (2-D) or 8 (3-D)
  int maxreg = 0;              // largest band region (cells) of one slot
  Kokkos::Array<int, 8> oj, ok;  // slot offsets in (lx2, lx3); opposite slot = nof-1-o
  DualArray1D<int> lx1;        // (nmb) x1 logical position of the block
  DualArray1D<int> hloc;       // (8 nmb) local index of the slot's neighbour, -1 remote
  DualArray1D<int> xloc;       // (2 nmb) x1 neighbour lo/hi: local index, -1 remote,
                               // -2 none (physical x1 face)
  std::vector<int> hrank, xrank;
  bool hmpi = false, xmpi = false;   // some horizontal / x1 neighbour is remote
  DvceArray5D<Real> csw;       // (nmb, nx1+4, 2, n3w, n2w) extinction and source, banded,
                               // layers 0, 1 and nx1+2, nx1+3 = the x1 neighbours' two
                               // face layers (0 beyond a physical x1 face)
  DvceArray5D<Real> ipl;       // (nmb, nring, nray, n3w, n2w) I of the last layers
  DvceArray5D<Real> spl;       // (nmb, 4, NFS, n3w, n2w) closure_consistency: face sums
                               // of the x1 neighbours' cells, for the lower (0) and the
                               // upper (1) one, received from the upper (2), lower (3)
  MeshBoundaryValuesCC *pbd = nullptr;  // vet_tensor = full: ghost exchange of D
  DvceArray5D<Real> dsc, dsc_c;         // its scratch array (nmb, 6, k, j, i) and dummy
#if MPI_PARALLEL_ENABLED
  MPI_Comm comm;
  std::vector<MPI_Request> req;  // the band exchange in flight (VetBandPost/VetBandWait)
#endif
  int pa0 = 0, pna = 0, pb0 = 0, pnb = 0;  // its plane range
  // sweep blocks and rays: the rank's own blocks and all rays (nsb = nmb, nrl = nray)
  int nsb = 0, gsb0 = 0;       // sweep blocks, gid of the first one
  int nrl = 0, nhl = 0;        // rays, of which up (mu1 > 0; they come first)
  std::vector<int> hlid, xlid;  // receiver's sweep index of a remote slot neighbour
  // one message per neighbour rank and exchange, pieces packed in the receiver's
  // (block, slot) order.  Plan P: intensity plane, plan B: csw band.
  struct Agg {
    DvceArray1D<int> soff, roff;     // (8 nsb) piece offset in the flat buffers, -1
    std::vector<int> prk, sdsp, scnt, rdsp, rcnt;   // per peer rank
    int stot = 0, rtot = 0;
  } aP, aB;
  DvceArray1D<Real> sfl, rfl;  // flat send / receive buffers
  int npb = 1;                 // band passes: a band wider than a block takes cells of
                               // the neighbours' band, filled by the pass before
};

// the band region of slot (dj, dk) on the RECEIVER: first plane index and extent
KOKKOS_INLINE_FUNCTION
static void VetRegion(const int d, const int nx, const int w, int &r0, int &nr) {
  r0 = (d < 0) ? 0 : ((d == 0) ? w : (w + nx));
  nr = (d == 0) ? nx : w;
}

// the pass of the band exchange that fills band cell (kk, jj): one per block width
// between the cell and the active cells (0 for an active cell)
static KOKKOS_INLINE_FUNCTION
int VetBandPass(const int kk, const int jj, const int nx2, const int nx3, const int b2,
                const int b3) {
  auto pass = [](const int c, const int nx, const int b) {
    const int d = (c < b) ? b - c : ((c >= b + nx) ? c - b - nx + 1 : 0);
    return (d + nx - 1)/nx;
  };
  const int p2 = pass(jj, nx2, b2), p3 = pass(kk, nx3, b3);
  return (p2 > p3) ? p2 : p3;
}

//----------------------------------------------------------------------------------------
//! \fn VetFor
//! \brief flat device loop over [0, n) with the functor as a kernel argument
//! (HintLightWeight): no host wait between the sweep's per-layer launches.  Functors
//! must capture Views (.d_view of DualViews), never DualViews or host structs.

template <class F>
static inline void VetFor(const char *name, const int n, const F &f) {
  static_assert(sizeof(F) <= 3072, "vet_sc kernel functor too large for a kernel-"
                "argument launch: capture Views, not DualViews or host structs");
  if (n <= 0) return;
  Kokkos::parallel_for(name, Kokkos::Experimental::require(
                       Kokkos::RangePolicy<>(DevExeSpace(), 0, n),
                       Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
}

//----------------------------------------------------------------------------------------
//! \fn VetAggPost
//! \brief post the one message per neighbour rank of plan A (st.sfl / st.rfl, filled by
//! the caller's completed pack kernel)

static inline void VetAggPost(ScCartMBState &st, const ScCartMBState::Agg &A) {
#if MPI_PARALLEL_ENABLED
  const int np = static_cast<int>(A.prk.size());
  for (int p = 0; p < np; ++p) {
    if (A.rcnt[p] == 0) continue;
    st.req.push_back(MPI_REQUEST_NULL);
    MPI_Irecv(st.rfl.data() + A.rdsp[p], A.rcnt[p], MPI_ATHENA_REAL, A.prk[p], 5,
              st.comm, &st.req.back());
  }
  for (int p = 0; p < np; ++p) {
    if (A.scnt[p] == 0) continue;
    st.req.push_back(MPI_REQUEST_NULL);
    MPI_Isend(st.sfl.data() + A.sdsp[p], A.scnt[p], MPI_ATHENA_REAL, A.prk[p], 5,
              st.comm, &st.req.back());
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn VetBandPost / VetBandWait
//! \brief fill the horizontal ghost band (b2 x b3) of the banded planes A(m, a0..+na,
//! b0..+nb, kk, jj) from the 2 or 8 horizontal neighbours: periodic wrap, device copy on
//! the rank, MPI between ranks.  A band wider than a block takes st.npb passes, pass p
//! filling the cells p block widths out (no race).  Post launches the local copy and
//! pack and posts the messages; Wait completes and unpacks.

static inline
void VetBandPost(ScCartMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                 const int b0, const int nb, const int nmb, const int nx2,
                 const int nx3, const int pass) {
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto a_ = a;
  st.pa0 = a0;
  st.pna = na;
  st.pb0 = b0;
  st.pnb = nb;
  const int nab = na*nb;
  VetFor("vet_vet_band_loc", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int ab = q%nab;
    q /= nab;
    const int o = q%nof, m = q/nof;
    const int n = hl_(8*m + o);
    if (n < 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(dj, nx2, b2, j0, nj);
    VetRegion(dk, nx3, b3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int jj = j0 + jr, kk = k0 + kr;
    if (VetBandPass(kk, jj, nx2, nx3, b2, b3) != pass) return;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    a_(m,aa,bb,kk,jj) = a_(n,aa,bb,kk - dk*nx3,jj - dj*nx2);
  });
#if MPI_PARALLEL_ENABLED
  st.req.clear();
  if (!st.hmpi) return;
  auto sf_ = st.sfl;
  auto so_ = st.aB.soff;
  VetFor("vet_vet_band_packa", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int ab = q%nab;
    q /= nab;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(-dj, nx2, b2, j0, nj);
    VetRegion(-dk, nx3, b3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    sf_(so_(8*m + o) + ab*nj*nk + idx) = a_(m,aa,bb,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
  });
  Kokkos::fence();
  VetAggPost(st, st.aB);
#endif
}

static inline
void VetBandWait(ScCartMBState &st, DvceArray5D<Real> &a, const int nmb, const int nx2,
                 const int nx3) {
#if MPI_PARALLEL_ENABLED
  if (!st.hmpi || st.req.empty()) return;
  MPI_Waitall(static_cast<int>(st.req.size()), st.req.data(), MPI_STATUSES_IGNORE);
  st.req.clear();
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg;
  const int a0 = st.pa0, b0 = st.pb0, nb = st.pnb, nab = st.pna*st.pnb;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto a_ = a;
  auto rf_ = st.rfl;
  auto ro_ = st.aB.roff;
  VetFor("vet_vet_band_unpka", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int ab = q%nab;
    q /= nab;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(dj, nx2, b2, j0, nj);
    VetRegion(dk, nx3, b3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    a_(m,aa,bb,k0 + kr,j0 + jr) = rf_(ro_(8*m + o) + ab*nj*nk + idx);
  });
#endif
}

static inline
void VetBandExchange(ScCartMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                     const int b0, const int nb, const int nmb, const int nx2,
                     const int nx3) {
  for (int p = 1; p <= st.npb; ++p) {
    VetBandPost(st, a, a0, na, b0, nb, nmb, nx2, nx3, p);
    VetBandWait(st, a, nmb, nx2, nx3);
  }
}

//----------------------------------------------------------------------------------------
//! \fn VetRayRegion
//! \brief band region of slot direction d for one ray: only the ghost cells it reads
//! (dlo below / dhi above the active cells)

static KOKKOS_INLINE_FUNCTION
void VetRayRegion(const int d, const int nx, const int b, const int dlo, const int dhi,
                  int &r0, int &nr) {
  if (d < 0) {
    r0 = b - dlo;
    nr = dlo;
  } else if (d == 0) {
    r0 = b;
    nr = nx;
  } else {
    r0 = b + nx;
    nr = dhi;
  }
}

//----------------------------------------------------------------------------------------
//! \fn VetPlanePost / VetPlaneWait
//! \brief band exchange of intensity plane p cut down to what each ray reads (one side,
//! floor(reach) + 1 cells); same values in every cell read

static inline
void VetPlanePost(ScCartMBState &st, DvceArray5D<Real> &a, const int p, const int nray,
                  const int nmb, const int nx2, const int nx3, const int pass = 1) {
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto rd_ = st.rdep;
  auto a_ = a;
  VetFor("vet_vet_plane_loc", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    const int n = hl_(8*m + o);
    if (n < 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(dj, nx2, b2, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(dk, nx3, b3, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int jj = j0 + jr, kk = k0 + kr;
    if (VetBandPass(kk, jj, nx2, nx3, b2, b3) != pass) return;
    a_(m,p,r,kk,jj) = a_(n,p,r,kk - dk*nx3,jj - dj*nx2);
  });
#if MPI_PARALLEL_ENABLED
  st.req.clear();
  st.pa0 = p;
  if (!st.hmpi) return;
  auto os_ = st.offs;
  auto sf_ = st.sfl;
  auto so_ = st.aP.soff;
  VetFor("vet_vet_plane_packa", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(-dj, nx2, b2, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(-dk, nx3, b3, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    sf_(so_(8*m + o) + os_(o,r) + idx) = a_(m,p,r,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
  });
  Kokkos::fence();
  VetAggPost(st, st.aP);
#endif
}

static inline
void VetPlaneWait(ScCartMBState &st, DvceArray5D<Real> &a, const int nray, const int nmb,
                  const int nx2, const int nx3) {
#if MPI_PARALLEL_ENABLED
  if (!st.hmpi || st.req.empty()) return;
  MPI_Waitall(static_cast<int>(st.req.size()), st.req.data(), MPI_STATUSES_IGNORE);
  st.req.clear();
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg;
  const int p = st.pa0;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto rd_ = st.rdep;
  auto or_ = st.offr;
  auto a_ = a;
  auto rf_ = st.rfl;
  auto ro_ = st.aP.roff;
  VetFor("vet_vet_plane_unpka", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(dj, nx2, b2, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(dk, nx3, b3, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    a_(m,p,r,k0 + kr,j0 + jr) = rf_(ro_(8*m + o) + or_(o,r) + idx);
  });
#endif
}

//----------------------------------------------------------------------------------------
//! \fn VetX1Move
//! \brief copy whole banded planes A(nbr, as, b0..+nb, :, :) -> A(m, ad, ...) across x1
//! block faces.  s = 0: m receives from its lower x1 neighbour, s = 1 from the upper;
//! only blocks with lx1 == rcv receive (-2: all).  Planes are contiguous (LayoutRight).

static inline
void VetX1Move(ScCartMBState &st, DvceArray5D<Real> &a, const int s, const int rcv,
               const int ad, const int as, const int b0, const int nb, const int nmb) {
  const int n2w = st.n2w, n3w = st.n3w;
  auto xl_ = st.xloc;
  auto lx_ = st.lx1;
  auto a_ = a;
  par_for("vet_vet_x1mv", DevExeSpace(), 0, nmb-1, 0, nb-1, 0, n3w-1, 0, n2w-1,
  KOKKOS_LAMBDA(const int m, const int b, const int kk, const int jj) {
    const int n = xl_.d_view(2*m + s);
    if (n < 0) return;
    if (rcv != -2 && lx_.d_view(m) != rcv) return;
    a_(m,ad,b0 + b,kk,jj) = a_(n,as,b0 + b,kk,jj);
  });
#if MPI_PARALLEL_ENABLED
  if (!st.xmpi) return;
  Kokkos::fence();
  const size_t cnt = static_cast<size_t>(nb)*n3w*n2w;
  auto plane = [&](const int m, const int aa) {
    return a.data() + (((static_cast<size_t>(m)*a.extent(1) + aa)*a.extent(2) + b0)
                       *a.extent(3))*a.extent(4);
  };
  std::vector<MPI_Request> req;
  for (int m = 0; m < nmb; ++m) {
    const int lx = st.lx1.h_view(m);
    // receive from my side-s neighbour
    if (st.xloc.h_view(2*m + s) == -1 && (rcv == -2 || lx == rcv)) {
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(plane(m, ad), static_cast<int>(cnt), MPI_ATHENA_REAL, st.xrank[2*m + s],
                16*m + 8 + s, st.comm, &req.back());
    }
    // send to my opposite-side neighbour, which receives on its side s
    const int so = 1 - s;
    const int lxn = (lx + ((so == 1) ? 1 : -1) + st.nbx1) % st.nbx1;
    if (st.xloc.h_view(2*m + so) == -1 && (rcv == -2 || lxn == rcv)) {
      const int rk = st.xrank[2*m + so];
      const int lidn = st.xlid[2*m + so];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(plane(m, as), static_cast<int>(cnt), MPI_ATHENA_REAL, rk,
                16*lidn + 8 + s, st.comm, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
#endif
}

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_CART_EXCHANGE_HPP_
