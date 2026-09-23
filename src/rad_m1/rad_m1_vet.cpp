//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vet.cpp
//! \brief <rad_m1>/closure = vet_sc: a VARIABLE EDDINGTON TENSOR for the multi-D implicit
//! M1 solve from a grey formal solution of the transfer equation by SHORT
//! CHARACTERISTICS on the M1 grid.
//!
//! WHAT IS SOLVED.  Once per hydro step, at the start-of-step state,
//!     mu . grad I = chi (S - I),   chi = rho (kappa_F + kappa_s)  (= opac M1_OP_T)
//! with the grey source (E units, eps = 4 pi I / c, so that E = <eps>)
//!     S = eps_th a T^4 + (1 - eps_th) E,   eps_th = min(kappa_P / chi, 1),
//! i.e. (kappa_a B + kappa_s J)/chi with J taken from the current M1 E (no lambda
//! iteration: E is what the moment solve delivers and it is lagged with the tensor).
//! The stellar tables carry no separate scattering opacity (kappa_s = 0, kappa_F =
//! Rosseland incl. electron scattering), so the absorption fraction is estimated by the
//! Planck/extinction ratio and CAPPED: S is a convex blend of a T^4 and E, never
//! negative.
//! (The form E + (kappa_P a T^4 - kappa_E E)/chi, the module's own net emission, is NOT
//! usable: kappa_P >> kappa_R in the iron bump amplifies T_gas - T_rad by kappa_P/chi and
//! drove S to 0 at tau ~ 2-9, K/J = 0.48 there, in the first 200 s test.)
//! The Eddington tensor D = K/J of the formal solution is projected onto the uniaxial
//! form the solver takes (vet_axis = flux, the default):
//!     D ~ (1-chi)/2 delta + (3 chi - 1)/2 n n,   n = H/|H|,   chi = n.D.n
//! with H the flux of the FORMAL SOLUTION (not the M1 flux of the cell): the trace and
//! the component along the flux are kept exactly, the two perpendicular ones averaged
//! (exact for a plane-parallel or axisymmetric field).  chi is clamped to [1/3, 1].
//! vet_axis = eigen takes chi = the largest eigenvalue of D and n its vector instead;
//! in the nearly isotropic interior (|D - 1/3| ~ 1e-4) that axis FLIPS between x1, x2
//! and x3 from cell to cell, D_22 jumps by 1.5 (chi - 1/3) between neighbours, and
//! the huge deep E turns the jump into a horizontal force (first 200 s test: bottom-cell
//! rms v2 1.8e4 cm/s, KE_2 3.5e4 x the T3_tau arm).  The flux axis is smooth.
//!
//! ANGLES.  Double-Gauss in mu_1 = Omega.x1 (vet_nmu nodes on each hemisphere) times
//! vet_nphi uniform azimuths about x1; no ray is horizontal.
//!
//! SWEEP.  Layer by layer in x1: the ray arriving at a cell centre of layer i is traced
//! back to the PLANE of cell centres of the upwind layer i -+ 1, where I, chi and S are
//! (bi)linearly interpolated (periodic wrap in x2/x3); along the segment the source is
//! linear in tau (the positive, first-order short-characteristics weights).  Upward and
//! downward rays are swept in the same launches (layer is+l up, ie-l down), each launch
//! a par_for over the (k, j) plane with the rays looped inside, so the moments need no
//! atomics.  3-D is the same code with the x3 shift switched on.
//!
//! BOUNDARIES.  Bottom (x1min face), upward rays: the diffusion intensity
//! eps = E + 3 (F . Omega) / c of the bottom cell's M1 state, carried half a cell with
//! S of that cell.  Top (x1max face), downward rays: vacuum.
//!
//! SEVERAL MESHBLOCKS / MPI RANKS (VetSweepMB).  The same sweep over the GLOBAL layers,
//! EXACT: the planes are banded by w = floor(max dx1 |mu_perp/mu_1|/dx_perp) + 1 cells
//! in x2 (x3), the band of each new intensity plane is exchanged with the horizontal
//! neighbours after every launch (device copy on a rank, MPI between ranks), and at an
//! x1 block face the upwind plane is handed to the next block of the stack.  A split in
//! x2/x3 runs every block of a layer in the same launch (no serialisation); a split in
//! x1 sweeps the stack in turn, as the gathered x1 line solve does.  The alternative,
//! block-Jacobi with the neighbours' intensities lagged by one sweep, is kept as the
//! DIAGNOSTIC vet_mb_lag = K (VetMBSweeps).  The six D_ab of vet_tensor = full get their
//! ghosts from the ordinary cell-centred exchange.
//!
//! MULTI-GPU OPTIONS (tests_m1/runs_3q_scscale/README.md; keys read only when given).
//! vet_mb_agg (default true, exact): one MPI message per neighbour rank and exchange.
//! vet_mb_mom_fuse (exact): the moments summed inside the ray launch (team scratch).
//! vet_mb_agroup = G (round-off): HYBRID decomposition, the rays split over the G
//! consecutive ranks of a group (one node), the space over the groups; the moments are
//! summed over the group in rank order.
//!
//! LIMITS (checked at start-up).  Periodic x2 (and x3), uniform mesh, every MeshBlock at
//! least w cells wide in x2 (x3) and 2 cells in x1; vet_milne needs one block along x1.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {

void VetFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in rad_m1_vet.cpp" << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

// Gauss-Legendre nodes and weights on [0,1] (weights sum to 1)
void GaussLegendre01(const int n, std::vector<double> &x, std::vector<double> &w) {
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(M_PI*(i + 0.75)/(n + 0.5));
    double pp = 1.0;
    for (int it = 0; it < 200; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int jj = 1; jj <= n; ++jj) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*jj - 1.0)*z*p2 - (jj - 1.0)*p3)/jj;
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;
      if (std::fabs(z - z1) < 1.0e-15) break;
    }
    x[i] = 0.5*(1.0 + z);
    w[i] = 1.0/((1.0 - z*z)*pp*pp);
  }
}

// periodic index in [0, n)
KOKKOS_INLINE_FUNCTION
int VetWrap(const int a, const int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

// largest eigenvalue and its unit eigenvector of a symmetric 3x3 matrix (cyclic Jacobi)
KOKKOS_INLINE_FUNCTION
void VetEigMax(Real a[3][3], Real &lmax, Real &e1, Real &e2, Real &e3) {
  Real v[3][3];
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 6; ++sw) {
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  int im = 0;
  if (a[1][1] > a[im][im]) im = 1;
  if (a[2][2] > a[im][im]) im = 2;
  lmax = a[im][im];
  e1 = v[0][im];
  e2 = v[1][im];
  e3 = v[2][im];
}

// all eigenvalues lam[p] and unit eigenvectors v[.][p] of a symmetric 3x3 matrix (the
// cyclic Jacobi of VetEigMax, run to round-off); a is destroyed
KOKKOS_INLINE_FUNCTION
void VetEigSym(Real a[3][3], Real lam[3], Real v[3][3]) {
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 12; ++sw) {
    Real off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off <= 1.0e-18*(fabs(a[0][0]) + fabs(a[1][1]) + fabs(a[2][2]))) break;
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  for (int p = 0; p < 3; ++p) {lam[p] = a[p][p];}
}

} // namespace

//----------------------------------------------------------------------------------------
//! \struct VetMBState
//! \brief the state of the SEVERAL-MESHBLOCK sweep (VetSweepMB): banded plane buffers,
//! neighbour tables, message buffers, and the block-Jacobi emulation of the diagnostic
//! <rad_m1>/vet_mb_lag.
//!
//! Plane index convention: a banded plane is (kk, jj), kk in [0, n3w), jj in [0, n2w),
//! with the active cells at kk = k - ks + w3, jj = j - js + w2; w2 (w3) is the band a
//! ray reaches per layer, floor(max dx1 |mu_2|/(|mu_1| dx2)) + 1 cells.

struct VetMBState {
  int w2 = 0, w3 = 0;          // reach of a ray per layer in x2, x3 (w3 = 0 in 2-D)
  int hk = 1;                  // vet_mb_halo: layers per band exchange
  int b2 = 0, b3 = 0;          // band widths = hk w2, hk w3 (the ghost band of a plane)
  bool ovl = false;            // vet_mb_overlap: interior of layer l+1 under the MPI
  bool raypar = true;          // vet_mb_kernel = ray (thread per ray) | cell (per column)
  int nbat = 1;                // vet_mb_mom_batch: launches per moment launch
  int nring = 2;               // intensity planes kept (max(nbat, 2))
  DvceArray2D<int> rdep;       // (nray, 4) cells of ghost band a ray reads per layer:
                               // below / above in x2, below / above in x3
  DvceArray2D<int> offs, offr;   // (nof, nray+1) per-ray offsets of the plane messages
  std::vector<int> cnts, cntr;   // (nof) plane message sizes (sent / received)
  int n2w = 1, n3w = 1;        // banded plane extents
  int nbx1 = 1, nx1g = 0;      // MeshBlocks along x1, global layers along x1
  int nof = 0;                 // horizontal neighbour slots: 2 (2-D) or 8 (3-D)
  int maxreg = 0;              // largest band region (cells) of one slot
  int maxcnt = 0;              // largest message (Reals) of one slot
  Kokkos::Array<int, 8> oj, ok;  // slot offsets in (lx2, lx3); opposite slot = nof-1-o
  DualArray1D<int> lx1;        // (nmb) x1 logical position of the block
  DualArray1D<int> hloc;       // (8 nmb) local index of the slot's neighbour, -1 remote
  DualArray1D<int> hself;      // (8 nmb) 1 if the slot's neighbour is the block itself
  DualArray1D<int> xloc;       // (2 nmb) x1 neighbour lo/hi: local index, -1 remote,
                               // -2 none (physical x1 face)
  std::vector<int> hgid, hrank, xgid, xrank;
  std::vector<int> gid0;       // (nranks) first global id of each rank
  bool hmpi = false, xmpi = false;   // some horizontal / x1 neighbour is remote
  DvceArray5D<Real> csw;       // (nmb, nx1+2, 2, n3w, n2w) extinction and source, banded,
                               // layers 0 and nx1+1 = the x1 neighbours' face layers
  DvceArray5D<Real> ipl;       // (nmb, nring, nray, n3w, n2w) I of the last layers
  DvceArray3D<Real> sbuf, rbuf;  // (nmb, nof, maxcnt) message buffers
  MeshBoundaryValuesCC *pbd = nullptr;  // vet_tensor = full: ghost exchange of D
  DvceArray5D<Real> dsc, dsc_c;         // its scratch array (nmb, 6, k, j, i) and dummy
  // DIAGNOSTIC vet_mb_lag = K > 0: block-Jacobi emulation (see VetMBSweeps)
  int nlag = 0;
  DvceArray5D<Real> hist;      // (nmb, nx1g, nray, n3w, n2w) the last sweep's bands
  DvceArray5D<Real> histx;     // (nmb, 2, nray, n3w, n2w) the last sweep's x1 planes
  DvceArray5D<Real> kref;      // (nmb, 7, k, j, i) J, K_ab of the exact sweep
  std::vector<double> e0max, emax, esum, erms;  // per sweep index: first call max; max,
  double ncalls = 0.0;                           // sum of max, sum of rms over calls >= 1
#if MPI_PARALLEL_ENABLED
  MPI_Comm comm;
  std::vector<MPI_Request> req;  // the band exchange in flight (VetBandPost/VetBandWait)
#endif
  int pa0 = 0, pna = 0, pb0 = 0, pnb = 0;  // its plane range
  // PROTOTYPE vet_mb_angles: the angle-decomposed sweep (VetSweepAng)
  bool ang = false;
  int nmbt = 0, nbx2 = 1, n2g = 1, n3g = 1, nrm = 0;
  std::vector<int> gcnt, gdsp, mcnt, mdsp;  // gather / moment message counts, offsets
  DvceArray2D<Real> rrcv;             // (nranks, my moment count) received partials
  DvceArray1D<Real> gsnd, gbuf;       // (nmb chunk), (nmbt chunk) gather buffers
  DvceArray4D<Real> gcs;              // (nx1g, 2, n3g, n2g) global extinction, source
  DvceArray3D<Real> gbot;             // (4, n3g, n2g) E, F1..3 of the bottom cells
  DualArray1D<int> glx;               // (3 nmbt) lx1, lx2, lx3 of every block
  DualArray1D<int> gofd;              // (nbx1 nbx2 nbx3) gid of each logical location
  DvceArray2D<Real> gdx;              // (nmbt, 3) dx1, dx2, dx3 of every block
  DvceArray1D<int> rlist;             // (nrm) my rays
  DvceArray4D<Real> gip;              // (2, nrm, n3g, n2g) I of my rays, last two layers
  DvceArray5D<Real> gmom;             // (nmbt, 10, nx3, nx2, nx1) partial moments
  // SWEEP BLOCKS and LOCAL RAYS.  The banded sweep runs over nsb blocks and nrl rays:
  // by default the rank's own blocks and all rays (nsb = nmb, nrl = nray).  With the
  // hybrid decomposition vet_mb_agroup = G > 1 the G consecutive ranks of a group (one
  // node) all hold the extinction and source of the GROUP's blocks (gathered once per
  // call) and each sweeps the rays r = p, p + G, ... (p = rank in the group) over them:
  // the per-layer band exchange stays inside the rank for group neighbours and goes to
  // the rank with the same p in the neighbour group otherwise; the partial moments go
  // to the owners of the blocks and are summed in group-rank order (round-off).
  int ng = 1, gp = 0;          // vet_mb_agroup: ranks per group, my rank in the group
  int nsb = 0, gsb0 = 0;       // sweep blocks, gid of the first one
  int sown = 0;                // sweep index of my first own block
  int nrl = 0, nhl = 0;        // local rays, of which up (mu1 > 0; they come first)
  DvceArray2D<Real> angl;      // (nrl, 4) mu1, mu2, mu3, weight of the local rays
  DvceArray4D<Real> bwe;       // (nsb, 4, nx3, nx2) E, F1..3 of the bottom (is) cells
  DvceArray2D<Real> sdx;       // (nsb, 3) dx1, dx2, dx3 of the sweep blocks
  DvceArray5D<Real> mom;       // (nsb, 10, nx1, nx3, nx2) the sweep's moments, j fastest
                               // (coalesced additions), into vet_cell by VetMomOut
  std::vector<int> hlid, xlid;  // receiver's sweep index of a remote slot neighbour
  std::vector<int> gpr;        // (G) world ranks of my group
  std::vector<int> gso, gsn;   // (G) first sweep index and block count of each
  // vet_mb_agg (default on): ONE message per neighbour rank and exchange; the pieces of
  // all (block, slot) pairs to one rank are packed into one buffer in the receiver's
  // (block, slot) order.  Plan P: the intensity plane, plan B: the csw band.
  bool agg = true;
  struct Agg {
    DvceArray1D<int> soff, roff;     // (8 nsb) piece offset in the flat buffers, -1
    std::vector<int> prk, sdsp, scnt, rdsp, rcnt;   // per peer rank
    int stot = 0, rtot = 0;
  } aP, aB;
  DvceArray1D<Real> sfl, rfl;  // flat send / receive buffers
  // vet_mb_mom_fuse: the moments are summed in the ray launch (team scratch)
  bool fuse = false;
};

namespace {

// the band region of slot (dj, dk) on the RECEIVER: first plane index and extent
KOKKOS_INLINE_FUNCTION
void VetRegion(const int d, const int nx, const int w, int &r0, int &nr) {
  r0 = (d < 0) ? 0 : ((d == 0) ? w : (w + nx));
  nr = (d == 0) ? nx : w;
}

//----------------------------------------------------------------------------------------
//! \fn VetFor
//! \brief a flat device loop over [0, n) launched with the functor as a KERNEL ARGUMENT
//! (HintLightWeight), i.e. without the host wait of the constant-memory path that a
//! 0.5-32 kB functor takes by default (see BvalsTeamFor in bvals.hpp).  The sweep issues
//! a few launches per x1 layer, so every one of them must be free of host syncs:
//! functors capture Views (.d_view of DualViews), never DualViews or host structs.

template <class F>
inline void VetFor(const char *name, const int n, const F &f) {
  static_assert(sizeof(F) <= 3072, "vet_sc kernel functor too large for a kernel-"
                "argument launch: capture Views, not DualViews or host structs");
  if (n <= 0) return;
  Kokkos::parallel_for(name, Kokkos::Experimental::require(
                       Kokkos::RangePolicy<>(DevExeSpace(), 0, n),
                       Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
}

//----------------------------------------------------------------------------------------
//! \fn VetAggPost
//! \brief vet_mb_agg: post the ONE message per neighbour rank of plan A (flat buffers
//! st.sfl / st.rfl, filled by the caller's pack kernel, which must have completed).

void VetAggPost(VetMBState &st, const VetMBState::Agg &A) {
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
//! \brief fill the horizontal ghost band (width b2 x b3) of the banded planes A(m,
//! a0..a0+na-1, b0..b0+nb-1, kk, jj) from the (2 or 8) horizontal neighbours' active
//! cells: periodic wrap, a device copy between blocks of this rank, MPI (device buffers,
//! like the boundary machinery) between ranks.  Ghost and source regions never overlap,
//! so the local copy has no race.  VetBandPost launches the local copy and the pack and
//! posts the messages; VetBandWait completes them and unpacks.  Work that does not read
//! the remote part of the band may be launched in between.

void VetBandPost(VetMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                 const int b0, const int nb, const int nmb, const int nx2,
                 const int nx3) {
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
  VetFor("m1_vet_band_loc", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
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
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    a_(m,aa,bb,kk,jj) = a_(n,aa,bb,kk - dk*nx3,jj - dj*nx2);
  });
#if MPI_PARALLEL_ENABLED
  st.req.clear();
  if (!st.hmpi) return;
  if (st.agg) {
    auto sf_ = st.sfl;
    auto so_ = st.aB.soff;
    VetFor("m1_vet_band_packa", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
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
    return;
  }
  auto sb_ = st.sbuf;
  // pack: slot o sends my cells that fill the neighbour's ghost region of slot nof-1-o
  // (offset -o), i.e. its ghost index + o * nx
  VetFor("m1_vet_band_pack", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
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
    sb_(m,o,ab*nj*nk + idx) = a_(m,aa,bb,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
  });
  Kokkos::fence();
  auto rb_ = st.rbuf;
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      if (st.hrank[8*m + o] == global_variable::my_rank) continue;
      int nj = (st.oj[o] == 0) ? nx2 : b2, nk = (st.ok[o] == 0) ? nx3 : b3;
      st.req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(rb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, nab*nj*nk,
                MPI_ATHENA_REAL, st.hrank[8*m + o], 16*m + o, st.comm, &st.req.back());
    }
  }
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      int rk = st.hrank[8*m + o];
      if (rk == global_variable::my_rank) continue;
      int nj = (st.oj[o] == 0) ? nx2 : b2, nk = (st.ok[o] == 0) ? nx3 : b3;
      int lidn = st.hlid[8*m + o];
      st.req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(sb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, nab*nj*nk,
                MPI_ATHENA_REAL, rk, 16*lidn + (nof - 1 - o), st.comm, &st.req.back());
    }
  }
#endif
}

void VetBandWait(VetMBState &st, DvceArray5D<Real> &a, const int nmb, const int nx2,
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
  auto rb_ = st.rbuf;
  auto a_ = a;
  if (st.agg) {
    auto rf_ = st.rfl;
    auto ro_ = st.aB.roff;
    VetFor("m1_vet_band_unpka", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
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
    return;
  }
  VetFor("m1_vet_band_unpk", nmb*nof*nab*mreg, KOKKOS_LAMBDA(const int t) {
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
    a_(m,aa,bb,k0 + kr,j0 + jr) = rb_(m,o,ab*nj*nk + idx);
  });
#endif
}

void VetBandExchange(VetMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                     const int b0, const int nb, const int nmb, const int nx2,
                     const int nx3) {
  VetBandPost(st, a, a0, na, b0, nb, nmb, nx2, nx3);
  VetBandWait(st, a, nmb, nx2, nx3);
}

//----------------------------------------------------------------------------------------
//! \fn VetRayRegion
//! \brief the band region of slot direction d for ONE ray: only the ghost cells that ray
//! reads (dlo below / dhi above the active cells per layer, times hk layers)

KOKKOS_INLINE_FUNCTION
void VetRayRegion(const int d, const int nx, const int b, const int hk, const int dlo,
                  const int dhi, int &r0, int &nr) {
  if (d < 0) {
    r0 = b - hk*dlo;
    nr = hk*dlo;
  } else if (d == 0) {
    r0 = b;
    nr = nx;
  } else {
    r0 = b + nx;
    nr = hk*dhi;
  }
}

//----------------------------------------------------------------------------------------
//! \fn VetPlanePost / VetPlaneWait
//! \brief the band exchange of the intensity plane p (VetBandPost/VetBandWait) cut down
//! to what each ray reads: a ray going up in x2 reads the band below only, and only
//! floor(reach) + 1 cells of it (1 cell for most rays; the band w is set by the most
//! oblique one).  Same values in every cell that is read.

void VetPlanePost(VetMBState &st, DvceArray5D<Real> &a, const int p, const int nray,
                  const int nmb, const int nx2, const int nx3) {
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg, hk = st.hk;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto rd_ = st.rdep;
  auto a_ = a;
  VetFor("m1_vet_plane_loc", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    const int n = hl_(8*m + o);
    if (n < 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(dj, nx2, b2, hk, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(dk, nx3, b3, hk, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int jj = j0 + jr, kk = k0 + kr;
    a_(m,p,r,kk,jj) = a_(n,p,r,kk - dk*nx3,jj - dj*nx2);
  });
#if MPI_PARALLEL_ENABLED
  st.req.clear();
  st.pa0 = p;
  if (!st.hmpi) return;
  auto os_ = st.offs;
  if (st.agg) {
    auto sf_ = st.sfl;
    auto so_ = st.aP.soff;
    VetFor("m1_vet_plane_packa", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
      const int idx = t%mreg;
      int q = t/mreg;
      const int r = q%nray;
      q /= nray;
      const int o = q%nof, m = q/nof;
      if (hl_(8*m + o) >= 0) return;
      const int dj = oj[o], dk = ok[o];
      int j0, nj, k0, nk;
      VetRayRegion(-dj, nx2, b2, hk, rd_(r,0), rd_(r,1), j0, nj);
      VetRayRegion(-dk, nx3, b3, hk, rd_(r,2), rd_(r,3), k0, nk);
      if (idx >= nj*nk) return;
      const int kr = idx/nj, jr = idx - kr*nj;
      sf_(so_(8*m + o) + os_(o,r) + idx) = a_(m,p,r,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
    });
    Kokkos::fence();
    VetAggPost(st, st.aP);
    return;
  }
  auto sb_ = st.sbuf;
  VetFor("m1_vet_plane_pack", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(-dj, nx2, b2, hk, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(-dk, nx3, b3, hk, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    sb_(m,o,os_(o,r) + idx) = a_(m,p,r,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
  });
  Kokkos::fence();
  auto rb_ = st.rbuf;
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      if (st.hrank[8*m + o] == global_variable::my_rank) continue;
      if (st.cntr[o] == 0) continue;
      st.req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(rb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, st.cntr[o],
                MPI_ATHENA_REAL, st.hrank[8*m + o], 16*m + o, st.comm, &st.req.back());
    }
  }
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      int rk = st.hrank[8*m + o];
      if (rk == global_variable::my_rank) continue;
      if (st.cnts[o] == 0) continue;
      int lidn = st.hlid[8*m + o];
      st.req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(sb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, st.cnts[o],
                MPI_ATHENA_REAL, rk, 16*lidn + (nof - 1 - o), st.comm, &st.req.back());
    }
  }
#endif
}

void VetPlaneWait(VetMBState &st, DvceArray5D<Real> &a, const int nray, const int nmb,
                  const int nx2, const int nx3) {
#if MPI_PARALLEL_ENABLED
  if (!st.hmpi || st.req.empty()) return;
  MPI_Waitall(static_cast<int>(st.req.size()), st.req.data(), MPI_STATUSES_IGNORE);
  st.req.clear();
  const int nof = st.nof, b2 = st.b2, b3 = st.b3, mreg = st.maxreg, hk = st.hk;
  const int p = st.pa0;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc.d_view;
  auto rd_ = st.rdep;
  auto or_ = st.offr;
  auto rb_ = st.rbuf;
  auto a_ = a;
  if (st.agg) {
    auto rf_ = st.rfl;
    auto ro_ = st.aP.roff;
    VetFor("m1_vet_plane_unpka", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
      const int idx = t%mreg;
      int q = t/mreg;
      const int r = q%nray;
      q /= nray;
      const int o = q%nof, m = q/nof;
      if (hl_(8*m + o) >= 0) return;
      const int dj = oj[o], dk = ok[o];
      int j0, nj, k0, nk;
      VetRayRegion(dj, nx2, b2, hk, rd_(r,0), rd_(r,1), j0, nj);
      VetRayRegion(dk, nx3, b3, hk, rd_(r,2), rd_(r,3), k0, nk);
      if (idx >= nj*nk) return;
      const int kr = idx/nj, jr = idx - kr*nj;
      a_(m,p,r,k0 + kr,j0 + jr) = rf_(ro_(8*m + o) + or_(o,r) + idx);
    });
    return;
  }
  VetFor("m1_vet_plane_unpk", nmb*nof*nray*mreg, KOKKOS_LAMBDA(const int t) {
    const int idx = t%mreg;
    int q = t/mreg;
    const int r = q%nray;
    q /= nray;
    const int o = q%nof, m = q/nof;
    if (hl_(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRayRegion(dj, nx2, b2, hk, rd_(r,0), rd_(r,1), j0, nj);
    VetRayRegion(dk, nx3, b3, hk, rd_(r,2), rd_(r,3), k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    a_(m,p,r,k0 + kr,j0 + jr) = rb_(m,o,or_(o,r) + idx);
  });
#endif
}

//----------------------------------------------------------------------------------------
//! \fn VetX1Move
//! \brief copy whole banded planes A(nbr, as, b0..b0+nb-1, :, :) -> A(m, ad, ...) across
//! the x1 faces between MeshBlocks.  Side s = 0: m receives from its LOWER x1 neighbour,
//! s = 1 from its UPPER one; only the blocks with lx1 == rcv receive (rcv = -2: every
//! block that has such a neighbour).  The planes are contiguous in A (LayoutRight), so
//! MPI moves them in place.

void VetX1Move(VetMBState &st, DvceArray5D<Real> &a, const int s, const int rcv,
               const int ad, const int as, const int b0, const int nb, const int nmb) {
  const int n2w = st.n2w, n3w = st.n3w;
  auto xl_ = st.xloc;
  auto lx_ = st.lx1;
  auto a_ = a;
  par_for("m1_vet_x1mv", DevExeSpace(), 0, nmb-1, 0, nb-1, 0, n3w-1, 0, n2w-1,
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
    const int lxn = (so == 1) ? (lx + 1) : (lx - 1);
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

//----------------------------------------------------------------------------------------
//! \struct VetRayK
//! \brief ONE RAY AT ONE CELL of launch L of the banded sweep: the intensity of local ray
//! r at the banded plane cell (kk, jj) of sweep block m, written to ip(m, pw, r, kk, jj).
//! The arithmetic is that of the per-cell kernel of VetSweepMB term for term; c0, s0
//! come from the banded copy csw (the same numbers as vet_cell), the bottom E, F from
//! bwe (the same numbers as iw at is), dx from sdx (the same numbers as mb_size).
//! Returns false (nothing written) where the ray does not run: outside its valid
//! overlap (ek of its own per-layer reaches, rdep), or its hemisphere has no layer of
//! block m in this launch.

struct VetRayK {
  int l, nx1, nx1g, b2, b3, ek, nx2, nx3, pw, pr;
  bool thrd, milne;
  Real fmil, cl, efl;
  DvceArray2D<int> rd;
  DvceArray4D<Real> bw;
  DvceArray5D<Real> cs, ip;
  DvceArray1D<int> lx;
  DvceArray2D<Real> ang, dx;

  KOKKOS_INLINE_FUNCTION
  bool operator()(const int m, const int r, const int kk, const int jj, Real &ivo) const {
    if (ek > 0) {
      // the overlap a ray still has valid upwind data for: ek of its own reaches
      if (jj < b2 - ek*rd(r,0) || jj >= b2 + nx2 + ek*rd(r,1)) return false;
      if (thrd && (kk < b3 - ek*rd(r,2) || kk >= b3 + nx3 + ek*rd(r,3))) return false;
    }
    const int lb = lx(m);
    const int lu = l - lb*nx1;
    const int ld = (nx1g - 1 - l) - lb*nx1;
    const bool ua = (lu >= 0 && lu < nx1), da = (ld >= 0 && ld < nx1);
    const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2);
    const bool up = (m1 > 0.0);
    if (up ? !ua : !da) return false;
    const Real dx1 = dx(m,0);
    const Real dx2 = dx(m,1);
    const Real dx3 = dx(m,2);
    const int li0 = (up ? lu : ld) + 1;   // csw layer of the cell
    const Real am1 = fabs(m1);
    const Real c0 = cs(m,li0,0,kk,jj);
    const Real s0 = cs(m,li0,1,kk,jj);
    Real iv;
    if (l == 0) {
      // the physical x1 face (active cells only: no overlap at l = 0)
      const Real su = fmax(1.5*s0 - 0.5*cs(m,up ? (li0 + 1) : (li0 - 1),1,kk,jj), 0.0);
      Real ib = 0.0;
      if (up) {
        if (milne) {
          ib = su + 3.0*fmil/cl*m1;
        } else {
          const int kb = kk - b3, jb = jj - b2;
          Real e = fmax(bw(m,0,kb,jb), efl);
          Real fo = bw(m,1,kb,jb)*m1 + bw(m,2,kb,jb)*m2 + bw(m,3,kb,jb)*m3;
          ib = fmax(e + 3.0*fo/cl, 0.0);
        }
      }
      const Real dtau = 0.5*c0*dx1/am1;
      const Real ex = exp(-dtau);
      Real w0, wu;
      if (dtau < 1.0e-3) {
        w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
        wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
      } else {
        const Real g = (1.0 - ex)/dtau;
        w0 = 1.0 - g;
        wu = g - ex;
      }
      iv = ib*ex + wu*su + w0*s0;
    } else {
      const int li = up ? (li0 - 1) : (li0 + 1);
      const Real sh2 = -dx1*m2/(am1*dx2);
      const Real fl2 = floor(sh2);
      const Real a2 = sh2 - fl2;
      const int o2 = static_cast<int>(fl2);
      const int ja = jj + o2;
      const int jb = ja + 1;
      int ka = kk, kb = kk;
      Real a3 = 0.0;
      if (thrd) {
        const Real sh3 = -dx1*m3/(am1*dx3);
        const Real fl3 = floor(sh3);
        a3 = sh3 - fl3;
        const int o3 = static_cast<int>(fl3);
        ka = kk + o3;
        kb = ka + 1;
      }
      const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
      const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
      const Real iup_v = waa*ip(m,pr,r,ka,ja) + wab*ip(m,pr,r,ka,jb)
                         + wba*ip(m,pr,r,kb,ja) + wbb*ip(m,pr,r,kb,jb);
      const Real cu = waa*cs(m,li,0,ka,ja) + wab*cs(m,li,0,ka,jb)
                      + wba*cs(m,li,0,kb,ja) + wbb*cs(m,li,0,kb,jb);
      const Real su = waa*cs(m,li,1,ka,ja) + wab*cs(m,li,1,ka,jb)
                      + wba*cs(m,li,1,kb,ja) + wbb*cs(m,li,1,kb,jb);
      const Real dtau = 0.5*(cu + c0)*dx1/am1;
      const Real ex = exp(-dtau);
      Real w0, wu;
      if (dtau < 1.0e-3) {
        w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
        wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
      } else {
        const Real g = (1.0 - ex)/dtau;
        w0 = 1.0 - g;
        wu = g - ex;
      }
      iv = iup_v*ex + wu*su + w0*s0;
    }
    ivo = fmax(iv, 0.0);
    ip(m,pw,r,kk,jj) = ivo;
    return true;
  }
};

//----------------------------------------------------------------------------------------
//! \fn VetRayLaunch
//! \brief launch L of the banded sweep with ONE THREAD PER (block, ray, cell): VetRayK at
//! the plane cells jj in [jlo, jlo+nj), kk in [klo, klo+nk) (banded indices; outside
//! the active area = the redundant overlap of vet_mb_halo), skipping the box [sj0, sj1)
//! x [sk0, sk1).  The moments are summed by VetMomLaunch.

void VetRayLaunch(const char *name, const int nmb, const int nray, const int jlo,
                  const int nj, const int klo, const int nk, const int sj0,
                  const int sj1, const int sk0, const int sk1, const VetRayK &rk) {
  VetFor(name, nmb*nray*nk*nj, KOKKOS_LAMBDA(const int t) {
    const int jx = t%nj;
    int q = t/nj;
    const int kx = q%nk;
    q /= nk;
    const int r = q%nray, m = q/nray;
    const int jj = jlo + jx, kk = klo + kx;
    if (jj >= sj0 && jj < sj1 && kk >= sk0 && kk < sk1) return;
    Real iv;
    rk(m, r, kk, jj, iv);
  });
}

//----------------------------------------------------------------------------------------
//! \fn VetMomAdd
//! \brief the ten moments (J, K11 K22 K33 K12 K13 K23, H1 H2 H3) of the rays r0..r1-1
//! whose hemisphere is `upw` (up: mu1 > 0) at one cell of intensity plane p, ADDED to
//! ac[]: the rays in the order and with the arithmetic of the per-cell kernel; the
//! intensities are loaded eight at a time so that the loads of a thread overlap.  rl_
//! maps a ray slot to its angle (the identity unless rl = true).

template <class IP>
KOKKOS_INLINE_FUNCTION
void VetMomAdd(Real ac[10], const bool upw, const int r0, const int r1, const IP &ipr,
               const DvceArray2D<Real> &ang_, const DvceArray1D<int> &rl_,
               const bool rl) {
  for (int q0 = r0; q0 < r1; q0 += 8) {
    Real iv[8];
    for (int u = 0; u < 8; ++u) {
      iv[u] = (q0 + u < r1) ? ipr(q0 + u) : 0.0;
    }
    for (int u = 0; u < 8 && q0 + u < r1; ++u) {
      const int r = rl ? rl_(q0 + u) : (q0 + u);
      const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
      if ((m1 > 0.0) != upw) continue;
      const Real a = wr*iv[u];
      ac[0] += a;
      ac[1] += a*m1*m1;
      ac[2] += a*m2*m2;
      ac[3] += a*m3*m3;
      ac[4] += a*m1*m2;
      ac[5] += a*m1*m3;
      ac[6] += a*m2*m3;
      ac[7] += a*m1;
      ac[8] += a*m2;
      ac[9] += a*m3;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn VetRayMomLaunch
//! \brief vet_mb_mom_fuse: VetRayLaunch and the moments in ONE launch.  A team takes one
//! row kk of block m and a chunk of <= 512/nray cells jj; one thread per (ray, cell) of
//! the chunk runs VetRayK and keeps the intensity in team scratch; after a barrier one
//! thread per (active cell, hemisphere) sums the rays of the hemisphere with VetMomAdd
//! (the arithmetic and order of VetMomLaunch) and adds the ten sums to the target cell
//! (up rays: layer lu, down rays: layer ld of the launch).  A cell gets its two sums
//! from two launches onto 0, or both from the up thread when lu == ld: the same numbers
//! as VetMomLaunch.  Target vt_(m, n0 + n, i0 + i, k0 + k, j0 + j) (block-local k, j, i;
//! st.mom, j fastest).

template <class VT>
void VetRayMomLaunch(const char *name, const int nmb, const int nray, const int nh,
                     const int jlo, const int nj, const int klo, const int nk,
                     const int sj0, const int sj1, const int sk0, const int sk1,
                     const VetRayK &rk, const VT &vt_, const int i0, const int j0,
                     const int k0, const int n0) {
  if (nj <= 0 || nk <= 0 || nmb <= 0) return;
  // one thread per (ray, cell) of the chunk: chunks of <= 512/nray cells
  const int chm = std::max(1, 512/nray);
  const int nch = (nj + chm - 1)/chm, ch = (nj + nch - 1)/nch;
  const int lg = nmb*nk*nch;
  const size_t scr = ScrArray1D<Real>::shmem_size(static_cast<size_t>(nray)*ch);
  const int b2 = rk.b2, b3 = rk.b3, nx2 = rk.nx2, nx3 = rk.nx3;
  const int l = rk.l, nx1 = rk.nx1, nx1g = rk.nx1g;
  auto lx_ = rk.lx;
  auto ang_ = rk.ang;
  DvceArray1D<int> nol_;   // VetMomAdd's ray list: not used (rl = false)
  // launch bounds 512: the default 1024-thread bound caps the registers of the ray
  // arithmetic (spills)
  using Pol = Kokkos::TeamPolicy<Kokkos::LaunchBounds<512, 1>>;
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
  const int tsz = std::min(512, nray*ch);
  Pol pol(DevExeSpace(), lg, tsz);
#else
  Pol pol(DevExeSpace(), lg, Kokkos::AUTO);
#endif
  Kokkos::parallel_for(name, Kokkos::Experimental::require(
                       pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
                       Kokkos::Experimental::WorkItemProperty::HintLightWeight),
  KOKKOS_LAMBDA(const Pol::member_type &tm) {
    const int lgi = tm.league_rank();
    const int c = lgi%nch;
    int q = lgi/nch;
    const int kx = q%nk, m = q/nk;
    const int jb = jlo + c*ch;
    const int nc = (jlo + nj - jb < ch) ? (jlo + nj - jb) : ch;
    const int kk = klo + kx;
    ScrArray1D<Real> si(tm.team_scratch(0), nray*ch);
    const bool skr = (kk >= sk0 && kk < sk1);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nray*nc), [&](const int u) {
      const int jx = u%nc, r = u/nc;
      const int jj = jb + jx;
      if (skr && jj >= sj0 && jj < sj1) return;
      Real iv = 0.0;
      rk(m, r, kk, jj, iv);
      si(r*ch + jx) = iv;
    });
    tm.team_barrier();
    if (kk < b3 || kk >= b3 + nx3) return;
    const int lb = lx_(m);
    const int lu = l - lb*nx1;
    const int ld = (nx1g - 1 - l) - lb*nx1;
    const bool ua = (lu >= 0 && lu < nx1), da = (ld >= 0 && ld < nx1);
    const bool same = ua && da && (lu == ld);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, 2*nc), [&](const int u) {
      const int jx = u%nc, h = u/nc;
      const int jj = jb + jx;
      if (jj < b2 || jj >= b2 + nx2) return;
      if (skr && jj >= sj0 && jj < sj1) return;
      if ((h == 0) ? !ua : !da) return;
      if (h == 1 && same) return;
      const int k = k0 + kk - b3, j = j0 + jj - b2;
      auto ipr = [&](const int r) { return si(r*ch + jx); };
      Real ac[10];
      for (int n = 0; n < 10; ++n) {
        ac[n] = 0.0;
      }
      VetMomAdd(ac, h == 0, (h == 0) ? 0 : nh, (h == 0) ? nh : nray, ipr, ang_, nol_,
                false);
      const int i = i0 + ((h == 0) ? lu : ld);
      for (int n = 0; n < 10; ++n) {
        vt_(m,n0 + n,i,k,j) += ac[n];
      }
      if (h == 0 && same) {
        for (int n = 0; n < 10; ++n) {
          ac[n] = 0.0;
        }
        VetMomAdd(ac, false, nh, nray, ipr, ang_, nol_, false);
        for (int n = 0; n < 10; ++n) {
          vt_(m,n0 + n,i,k,j) += ac[n];
        }
      }
    });
  });
}

//----------------------------------------------------------------------------------------
//! \fn VetMomLaunch
//! \brief the moments of launches L0..L1 (a batch of vet_mb_mom_batch) from the planes
//! VetRayLaunch wrote (ring slot L % nring): one thread per (block, launch, hemisphere,
//! active cell), the ten moments of that hemisphere's rays.  A cell receives exactly two
//! contributions per sweep, the up rays of launch G and the down rays of launch
//! nx1g-1-G, added to 0; two additions onto 0 are exact in either order, and when both
//! fall in one batch the up thread adds both (no race).  So the moments are the numbers
//! of the per-cell kernel.  The rays of the up hemisphere are 0..nh-1.  Target
//! vc_(m, n0 + n, i0 + i, k0 + k, j0 + j): st.mom (block, n, i, k, j), j fastest, so that
//! the additions are coalesced (VetMomOut hands it to vet_cell).

void VetMomLaunch(const int nmb, const int nray, const int l0, const int l1,
                  const int nring, const int nx1, const int nx1g, const int is,
                  const int js, const int ks, const int nx2, const int nx3, const int b2,
                  const int b3, const DvceArray5D<Real> &vc_,
                  const DvceArray5D<Real> &ip_, const DvceArray1D<int> &lx_,
                  const DvceArray2D<Real> &ang_, const DvceArray1D<int> &rl_,
                  const int nh, const int n0) {
  const int nl2 = 2*(l1 - l0 + 1);
  VetFor("m1_vet_mb_mom", nmb*nl2*nx3*nx2, KOKKOS_LAMBDA(const int t) {
    const int jx = t%nx2;
    int q = t/nx2;
    const int kx = q%nx3;
    q /= nx3;
    const int s = q%nl2, m = q/nl2;
    const int lc = l0 + s/2;
    const bool upt = ((s & 1) == 0);
    const int g = upt ? lc : (nx1g - 1 - lc);   // the target global layer
    const int li = g - lx_(m)*nx1;
    if (li < 0 || li >= nx1) return;
    const int lo = upt ? (nx1g - 1 - g) : g;    // the launch of its other hemisphere
    const bool both = (lo >= l0 && lo <= l1);
    if (!upt && both) return;
    const int kc = b3 + kx, jc = b2 + jx;
    Real ac[10];
    for (int n = 0; n < 10; ++n) {
      ac[n] = 0.0;
    }
    const int pc = lc%nring;
    auto ipr = [&](const int r) { return ip_(m,pc,r,kc,jc); };
    VetMomAdd(ac, upt, upt ? 0 : nh, upt ? nh : nray, ipr, ang_, rl_, false);
    Real ad[10];
    if (both) {
      for (int n = 0; n < 10; ++n) {
        ad[n] = 0.0;
      }
      const int po = lo%nring;
      auto ipo = [&](const int r) { return ip_(m,po,r,kc,jc); };
      VetMomAdd(ad, false, nh, nray, ipo, ang_, rl_, false);
    }
    const int k = ks + kx, j = js + jx, i = is + li;
    for (int n = 0; n < 10; ++n) {
      vc_(m,n0+n,i,k,j) += ac[n];
      if (both) {
        vc_(m,n0+n,i,k,j) += ad[n];
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn VetMomOut
//! \brief the moments of the banded sweep, accumulated in mom (sweep block, n, i, k, j)
//! with j fastest (so that the per-launch additions are coalesced), into vet_cell
//! (m, M1_VET_J + n, k, j, i): a tiled transpose through team scratch, one team per
//! (block, moment, k, chunk of 16 j).  ng > 1 (vet_mb_agroup): the value is the sum of
//! the group's partials in group-rank order (mine from mom(so + m), the others from rv_,
//! laid out like my blocks' slice of mom); otherwise mom(so + m) itself, which is then
//! zeroed for the next call (zero = true).  A copy: the numbers are those the moment
//! kernels would have added into vet_cell directly.

void VetMomOut(const DvceArray5D<Real> &mo_, const DvceArray2D<Real> &rv_, const int so,
               const int ng, const int gp, const bool zero, const DvceArray5D<Real> &vc_,
               const int nmb, const int nx1, const int nx2, const int nx3, const int is,
               const int js, const int ks) {
  const int cj = 16;
  const int nch = (nx2 + cj - 1)/cj;
  const int lg = nmb*10*nx3*nch;
  const size_t scr = ScrArray1D<Real>::shmem_size(static_cast<size_t>(nx1)*cj);
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
  Kokkos::TeamPolicy<> pol(DevExeSpace(), lg, 256);
#else
  Kokkos::TeamPolicy<> pol(DevExeSpace(), lg, Kokkos::AUTO);
#endif
  Kokkos::parallel_for("m1_vet_mb_momout", Kokkos::Experimental::require(
                       pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
                       Kokkos::Experimental::WorkItemProperty::HintLightWeight),
  KOKKOS_LAMBDA(const TeamMember_t &tm) {
    const int lgi = tm.league_rank();
    const int c = lgi%nch;
    int q = lgi/nch;
    const int k = q%nx3;
    q /= nx3;
    const int n = q%10, m = q/10;
    const int j0 = c*cj;
    const int ncj = (nx2 - j0 < cj) ? (nx2 - j0) : cj;
    ScrArray1D<Real> tl(tm.team_scratch(0), nx1*cj);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nx1*ncj), [&](const int u) {
      const int jj = u%ncj, i = u/ncj;
      const int j = j0 + jj;
      Real v;
      if (ng > 1) {
        const int t = (((m*10 + n)*nx1 + i)*nx3 + k)*nx2 + j;
        Real sum = 0.0;
        for (int p = 0; p < ng; ++p) {
          const Real w = (p == gp) ? mo_(so + m,n,i,k,j) : rv_(p,t);
          sum = (p == 0) ? w : (sum + w);
        }
        v = sum;
      } else {
        v = mo_(so + m,n,i,k,j);
        if (zero) mo_(so + m,n,i,k,j) = 0.0;
      }
      tl(i*cj + jj) = v;
    });
    tm.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nx1*ncj), [&](const int u) {
      const int i = u%nx1, jj = u/nx1;
      vc_(m,M1_VET_J + n,ks + k,js + j0 + jj,is + i) = tl(i*cj + jj);
    });
  });
}

//----------------------------------------------------------------------------------------
//! \fn VetSweepAng
//! \brief vet_mb_angles = true (PROTOTYPE): ANGLE DECOMPOSITION of the sweep.  Every rank
//! holds the extinction and source of the WHOLE mesh (one all-to-all gather per call)
//! and sweeps its share of the rays (r = p, p + P, ...) over the whole mesh as one
//! periodic block, one thread per (ray, cell), with no exchange inside the sweep; the
//! partial moments of all ranks are sent to the owners of the cells and summed there in
//! rank order.  The arithmetic of a ray is that of the banded sweep (with the dx
//! of the block that holds the cell); only the order of the sum over rays changes (by
//! rank), i.e. D equals the exact sweep to round-off, and on ONE rank the moments are
//! the very same numbers.

void VetSweepAng(VetMBState &st, const int nmb, const int is, const int js,
                 const int ks, const int nx1, const int nx2, const int nx3,
                 const bool thrd, const bool milne, const Real fmil, const Real cl,
                 const Real efl, const DvceArray5D<Real> &iw_,
                 const DvceArray5D<Real> &vc_, const DvceArray2D<Real> &ang_) {
  const int nx1g = st.nx1g, n2g = st.n2g, n3g = st.n3g, nbx2 = st.nbx2;
  const int nbx1 = st.nbx1;
  const int ncb = nx3*nx2*nx1, npb = nx3*nx2;
  const int chunk = 2*ncb + 4*npb;
  auto gs_ = st.gsnd;
  auto gb_ = st.gbuf;
  // (1) pack this rank's extinction, source and bottom-plane E, F; gather everywhere
  VetFor("m1_vet_ang_pack", nmb*chunk, KOKKOS_LAMBDA(const int t) {
    const int m = t/chunk;
    int q = t - m*chunk;
    if (q < 2*ncb) {
      const int c = q/ncb;
      q -= c*ncb;
      const int k = q/(nx2*nx1);
      q -= k*nx2*nx1;
      const int j = q/nx1, i = q - j*nx1;
      gs_(t) = vc_(m,(c == 0) ? M1_VET_CHX : M1_VET_SRC,ks + k,js + j,is + i);
    } else {
      q -= 2*ncb;
      const int c = q/npb;
      q -= c*npb;
      const int k = q/nx2, j = q - k*nx2;
      const int n = (c == 0) ? M1_IW_EN : ((c == 1) ? M1_IW_F1
                                       : ((c == 2) ? M1_IW_F2 : M1_IW_F3));
      gs_(t) = iw_(m,n,ks + k,js + j,is);
    }
  });
  Kokkos::fence();
  // the gather as point-to-point messages (MPI collectives on device buffers stage
  // through the host: 0.1 s per call measured with MPI_Allgatherv/Reduce_scatter)
  const int np = global_variable::nranks, me = global_variable::my_rank;
  {
    auto dst = Kokkos::subview(gb_, std::make_pair(st.gdsp[me], st.gdsp[me] + nmb*chunk));
    Kokkos::deep_copy(dst, gs_);
  }
#if MPI_PARALLEL_ENABLED
  if (np > 1) {
    std::vector<MPI_Request> rq;
    for (int p = 0; p < np; ++p) {
      if (p == me) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(gb_.data() + st.gdsp[p], st.gcnt[p], MPI_ATHENA_REAL, p, 1, st.comm,
                &rq.back());
    }
    Kokkos::fence();
    for (int p = 0; p < np; ++p) {
      if (p == me) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Isend(gs_.data(), nmb*chunk, MPI_ATHENA_REAL, p, 1, st.comm, &rq.back());
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
  }
#endif
  // (2) the global planes cs(i, c, k, j) and the bottom plane (c, k, j)
  auto gc_ = st.gcs;
  auto gbt_ = st.gbot;
  auto gl_ = st.glx.d_view;
  const int nmbt = st.nmbt;
  VetFor("m1_vet_ang_scat", nmbt*chunk, KOKKOS_LAMBDA(const int t) {
    const int g = t/chunk;
    int q = t - g*chunk;
    const int l1 = gl_(3*g), l2 = gl_(3*g + 1), l3 = gl_(3*g + 2);
    if (q < 2*ncb) {
      const int c = q/ncb;
      q -= c*ncb;
      const int k = q/(nx2*nx1);
      q -= k*nx2*nx1;
      const int j = q/nx1, i = q - j*nx1;
      gc_(l1*nx1 + i,c,l3*nx3 + k,l2*nx2 + j) = gb_(t);
    } else {
      if (l1 != 0) return;
      q -= 2*ncb;
      const int c = q/npb;
      q -= c*npb;
      const int k = q/nx2, j = q - k*nx2;
      gbt_(c,l3*nx3 + k,l2*nx2 + j) = gb_(t);
    }
  });
  // (3) the sweep over the global layers, my rays only: launch l does layer l (up rays)
  // and nx1g-1-l (down rays), as the banded sweep
  auto gp_ = st.gip;
  auto gm_ = st.gmom;
  auto rl_ = st.rlist;
  auto gof_ = st.gofd.d_view;
  auto gdx_ = st.gdx;
  const int nrm = st.nrm, nring = st.nring, nbat = st.nbat;
  for (int l = 0; l < nx1g; ++l) {
    const int pw = l%nring, pr = (l + nring - 1)%nring;
    VetFor("m1_vet_ang_ray", nrm*n3g*n2g, KOKKOS_LAMBDA(const int t) {
      const int j = t%n2g;
      int q = t/n2g;
      const int k = q%n3g, rr = q/n3g;
      const int r = rl_(rr);
      const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2);
      const bool up = (m1 > 0.0);
      const int i = up ? l : (nx1g - 1 - l);
      const int g = gof_(((k/nx3)*nbx2 + j/nx2)*nbx1 + i/nx1);
      const Real dx1 = gdx_(g,0), dx2 = gdx_(g,1), dx3 = gdx_(g,2);
      const Real am1 = fabs(m1);
      const Real c0 = gc_(i,0,k,j);
      const Real s0 = gc_(i,1,k,j);
      Real iv;
      if (l == 0) {
        const Real su = fmax(1.5*s0 - 0.5*gc_(up ? (i + 1) : (i - 1),1,k,j), 0.0);
        Real ib = 0.0;
        if (up) {
          if (milne) {
            ib = su + 3.0*fmil/cl*m1;
          } else {
            Real e = fmax(gbt_(0,k,j), efl);
            Real fo = gbt_(1,k,j)*m1 + gbt_(2,k,j)*m2 + gbt_(3,k,j)*m3;
            ib = fmax(e + 3.0*fo/cl, 0.0);
          }
        }
        const Real dtau = 0.5*c0*dx1/am1;
        const Real ex = exp(-dtau);
        Real w0, wu;
        if (dtau < 1.0e-3) {
          w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
          wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
        } else {
          const Real gg = (1.0 - ex)/dtau;
          w0 = 1.0 - gg;
          wu = gg - ex;
        }
        iv = ib*ex + wu*su + w0*s0;
      } else {
        const int iup = up ? (i - 1) : (i + 1);
        const Real sh2 = -dx1*m2/(am1*dx2);
        const Real fl2 = floor(sh2);
        const Real a2 = sh2 - fl2;
        const int o2 = static_cast<int>(fl2);
        const int ja = VetWrap(j + o2, n2g);
        const int jb = VetWrap(j + o2 + 1, n2g);
        int ka = k, kb = k;
        Real a3 = 0.0;
        if (thrd) {
          const Real sh3 = -dx1*m3/(am1*dx3);
          const Real fl3 = floor(sh3);
          a3 = sh3 - fl3;
          const int o3 = static_cast<int>(fl3);
          ka = VetWrap(k + o3, n3g);
          kb = VetWrap(k + o3 + 1, n3g);
        }
        const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
        const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
        const Real iup_v = waa*gp_(pr,rr,ka,ja) + wab*gp_(pr,rr,ka,jb)
                           + wba*gp_(pr,rr,kb,ja) + wbb*gp_(pr,rr,kb,jb);
        const Real cu = waa*gc_(iup,0,ka,ja) + wab*gc_(iup,0,ka,jb)
                        + wba*gc_(iup,0,kb,ja) + wbb*gc_(iup,0,kb,jb);
        const Real su = waa*gc_(iup,1,ka,ja) + wab*gc_(iup,1,ka,jb)
                        + wba*gc_(iup,1,kb,ja) + wbb*gc_(iup,1,kb,jb);
        const Real dtau = 0.5*(cu + c0)*dx1/am1;
        const Real ex = exp(-dtau);
        Real w0, wu;
        if (dtau < 1.0e-3) {
          w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
          wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
        } else {
          const Real gg = (1.0 - ex)/dtau;
          w0 = 1.0 - gg;
          wu = gg - ex;
        }
        iv = iup_v*ex + wu*su + w0*s0;
      }
      gp_(pw,rr,k,j) = fmax(iv, 0.0);
    });
    if ((l + 1)%nbat != 0 && l != nx1g - 1) continue;
    // the partial moments of my rays of launches l0..l, into the block-ordered sum
    // buffer: two contributions per cell and sweep, onto 0 (see VetMomLaunch)
    const int l0 = l - l%nbat, nl2 = 2*(l - l0 + 1), l1 = l;
    VetFor("m1_vet_ang_mom", nl2*n3g*n2g, KOKKOS_LAMBDA(const int t) {
      const int j = t%n2g;
      int q = t/n2g;
      const int k = q%n3g, s = q/n3g;
      const int lc = l0 + s/2;
      const bool upt = ((s & 1) == 0);
      const int g = upt ? lc : (nx1g - 1 - lc);
      const int lo = upt ? (nx1g - 1 - g) : g;
      const bool both = (lo >= l0 && lo <= l1);
      if (!upt && both) return;
      Real ac[10], ad[10];
      for (int n = 0; n < 10; ++n) {
        ac[n] = 0.0;
        ad[n] = 0.0;
      }
      const int pc = lc%nring, po = lo%nring;
      auto ipr = [&](const int r) { return gp_(pc,r,k,j); };
      VetMomAdd(ac, upt, 0, nrm, ipr, ang_, rl_, true);
      if (both) {
        auto ipo = [&](const int r) { return gp_(po,r,k,j); };
        VetMomAdd(ad, false, 0, nrm, ipo, ang_, rl_, true);
      }
      const int kb = k/nx3, jb = j/nx2;
      const int gg = gof_((kb*nbx2 + jb)*nbx1 + g/nx1);
      const int kl = k - kb*nx3, jl = j - jb*nx2, il = g - (g/nx1)*nx1;
      for (int n = 0; n < 10; ++n) {
        gm_(gg,n,kl,jl,il) += ac[n];
        if (both) {
          gm_(gg,n,kl,jl,il) += ad[n];
        }
      }
    });
  }
  // (4) the partial moments of every rank's blocks go to their owner (point-to-point),
  // which sums them in rank order into vet_cell
  Kokkos::fence();
  auto rv_ = st.rrcv;
  const int mc = st.mcnt[me];
#if MPI_PARALLEL_ENABLED
  if (np > 1) {
    std::vector<MPI_Request> rq;
    for (int p = 0; p < np; ++p) {
      if (p == me) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(rv_.data() + static_cast<size_t>(p)*mc, mc, MPI_ATHENA_REAL, p, 2,
                st.comm, &rq.back());
    }
    for (int p = 0; p < np; ++p) {
      if (p == me) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Isend(gm_.data() + st.mdsp[p], st.mcnt[p], MPI_ATHENA_REAL, p, 2, st.comm,
                &rq.back());
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
  }
#endif
  const int g0 = st.mdsp[me]/(10*ncb);
  VetFor("m1_vet_ang_unpk", nmb*10*ncb, KOKKOS_LAMBDA(const int t) {
    const int i = t%nx1;
    int q = t/nx1;
    const int j = q%nx2;
    q /= nx2;
    const int k = q%nx3;
    q /= nx3;
    const int n = q%10, m = q/10;
    Real sum = 0.0;
    for (int p = 0; p < np; ++p) {
      const Real v = (p == me) ? gm_(g0 + m,n,k,j,i) : rv_(p,t);
      sum = (p == 0) ? v : (sum + v);
    }
    vc_(m,M1_VET_J + n,ks + k,js + j,is + i) = sum;
  });
  Kokkos::deep_copy(gm_, 0.0);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetInit

void RadiationM1::VetInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  if (!trans_on) {
    VetFatal("<rad_m1>/closure = vet_sc needs transport = implicit on a MULTI-D mesh");
  }
  // several MeshBlocks or ranks: the exact banded sweep of VetSweepMB.  vet_mb_force
  // (DIAGNOSTIC, read only when given) runs it on one MeshBlock too, where it must
  // reproduce the single-block sweep bit for bit.
  bool mbpath = (pm->nmb_total != 1 || global_variable::nranks != 1);
  if (pin->DoesParameterExist("rad_m1", "vet_mb_force")) {
    mbpath = mbpath || pin->GetBoolean("rad_m1", "vet_mb_force");
  }
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
      (trans_x3 && pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic)) {
    VetFatal("<rad_m1>/closure = vet_sc needs periodic x2 (and x3) boundaries");
  }
  vet_nmu = pin->GetOrAddInteger("rad_m1", "vet_nmu", 4);
  vet_nphi = pin->GetOrAddInteger("rad_m1", "vet_nphi", 8);
  vet_milne = pin->GetOrAddBoolean("rad_m1", "vet_milne", false);
  vet_dump = pin->GetOrAddString("rad_m1", "vet_dump", "none");
  if (vet_dump.compare("none") == 0) {vet_dump = "";}
  vet_dump_every = pin->GetOrAddInteger("rad_m1", "vet_dump_every", 0);
  {std::string ax = pin->GetOrAddString("rad_m1", "vet_axis", "flux");
  if (ax.compare("flux") == 0) {
    vet_axis_flux = true;
  } else if (ax.compare("eigen") == 0) {
    vet_axis_flux = false;
  } else {
    VetFatal("<rad_m1>/vet_axis = '" + ax + "' not implemented (flux | eigen)");
  }
  }
  // vet_tensor = uniaxial (default) | full.  Read only when given, so that an input
  // without it writes the restart file the uniaxial path always wrote.
  vet_full = false;
  vet_eig_min = 0.0;
  if (pin->DoesParameterExist("rad_m1", "vet_tensor")) {
    std::string vt = pin->GetString("rad_m1", "vet_tensor");
    if (vt.compare("full") == 0) {
      vet_full = true;
    } else if (vt.compare("uniaxial") != 0) {
      VetFatal("<rad_m1>/vet_tensor = '" + vt + "' not implemented (uniaxial | full)");
    }
  }
  if (pin->DoesParameterExist("rad_m1", "vet_eig_min")) {
    vet_eig_min = pin->GetReal("rad_m1", "vet_eig_min");
  }
  if (vet_eig_min < 0.0 || vet_eig_min > 0.33) {
    VetFatal("<rad_m1>/vet_eig_min must lie in [0, 0.33]");
  }
  if (vet_nmu < 1 || vet_nphi < 1) {VetFatal("<rad_m1>/vet_nmu, vet_nphi must be >= 1");}
  // vet_x1_periodic = true (default false; tests_m1/runs_3r_radwave): a PERIODIC x1
  // sweep for a periodic x1 mesh.  The sweep is repeated vet_x1_npass times, each pass
  // entering layer is (ie) from the upward (downward) intensity the previous pass left
  // in layer ie (is); the first pass uses the ordinary boundary intensities and only
  // the last pass accumulates the moments.  The inflow error decays as exp(-tau) along
  // each ray per box crossing.  Single MeshBlock only.  Read only when named.
  vet_x1per = false;
  vet_x1npass = 1;
  if (pin->DoesParameterExist("rad_m1", "vet_x1_periodic")) {
    vet_x1per = pin->GetBoolean("rad_m1", "vet_x1_periodic");
  }
  if (vet_x1per) {
    if (mbpath) {VetFatal("<rad_m1>/vet_x1_periodic needs a single MeshBlock");}
    if (vet_milne) {VetFatal("<rad_m1>/vet_x1_periodic excludes vet_milne");}
    if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
      VetFatal("<rad_m1>/vet_x1_periodic needs periodic x1 boundaries");
    }
    vet_x1npass = pin->GetOrAddInteger("rad_m1", "vet_x1_npass", 64);
    if (vet_x1npass < 1) {VetFatal("<rad_m1>/vet_x1_npass must be >= 1");}
  }
  vet_nray = 2*vet_nmu*vet_nphi;

  // the angle set: (mu1, mu2, mu3, weight), upward hemisphere first
  std::vector<double> xg, wg;
  GaussLegendre01(vet_nmu, xg, wg);
  Kokkos::realloc(vet_ang, vet_nray, 4);
  auto ah = Kokkos::create_mirror_view(vet_ang);
  int r = 0;
  for (int h = 0; h < 2; ++h) {
    for (int a = 0; a < vet_nmu; ++a) {
      double m1 = (h == 0) ? xg[a] : -xg[a];
      double st = std::sqrt(1.0 - m1*m1);
      for (int l = 0; l < vet_nphi; ++l) {
        double ph = 2.0*M_PI*(l + 0.5)/vet_nphi;
        ah(r, 0) = m1;
        ah(r, 1) = st*std::cos(ph);
        ah(r, 2) = st*std::sin(ph);
        ah(r, 3) = 0.5*wg[a]/vet_nphi;
        ++r;
      }
    }
  }
  Kokkos::deep_copy(vet_ang, ah);

  auto &indcs = pm->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pm->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(vet_cell, nmb, M1_VET_NC, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(vet_cell, 0.0);
  if (mbpath) {
    VetMBInit(pin);   // the banded plane buffers live in vet_mbs
  } else {
    Kokkos::realloc(vet_ipl, nmb, 2, vet_nray, ncells3, ncells2);
    Kokkos::deep_copy(vet_ipl, 0.0);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> closure = vet_sc: short characteristics, " << vet_nray
              << " rays (double-Gauss nmu=" << vet_nmu << " x nphi=" << vet_nphi
              << "), lagged per hydro step" << (vet_milne ? ", MILNE source (diag)" : "")
              << (vet_full ? ", FULL tensor D = K/J (eigenvalue floor " : "")
              << (vet_full ? std::to_string(vet_eig_min) + ")" : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetMBInit
//! \brief several MeshBlocks / ranks: band widths, neighbour tables, buffers.  Uniform
//! mesh (the implicit solve already fatals on SMR/AMR), periodic x2 (x3), every block at
//! least as wide as the band.

void RadiationM1::VetMBInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  vet_mbs = new VetMBState;
  VetMBState &st = *vet_mbs;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const bool thrd = trans_x3;
  if ((mindcs.nx1 % nx1) != 0 || (mindcs.nx2 % nx2) != 0 ||
      (thrd && (mindcs.nx3 % nx3) != 0)) {
    VetFatal("<rad_m1>/closure = vet_sc on several MeshBlocks needs a uniform mesh");
  }
  st.nbx1 = mindcs.nx1/nx1;
  const int nbx2 = mindcs.nx2/nx2;
  const int nbx3 = thrd ? (mindcs.nx3/nx3) : 1;
  st.nx1g = mindcs.nx1;
  if (vet_milne && st.nbx1 > 1) {
    VetFatal("<rad_m1>/vet_milne needs ONE MeshBlock along x1 (tau from the top)");
  }
  if (nx1 < 2) {
    VetFatal("<rad_m1>/closure = vet_sc needs <meshblock>/nx1 >= 2");
  }
  const int me = global_variable::my_rank;
  const int nrk = global_variable::nranks;
  // vet_mb_agg (exact, default on): one message per neighbour rank and exchange;
  // vet_mb_mom_fuse (exact, default off): the moments summed inside the ray launch;
  // vet_mb_agroup = G (round-off, default 1 = off): the hybrid decomposition, rays split
  // over the G consecutive ranks of a group (one node), space over the groups.  All
  // three are read only when given.
  if (pin->DoesParameterExist("rad_m1", "vet_mb_agg")) {
    st.agg = pin->GetBoolean("rad_m1", "vet_mb_agg");
  }
  if (pin->DoesParameterExist("rad_m1", "vet_mb_mom_fuse")) {
    st.fuse = pin->GetBoolean("rad_m1", "vet_mb_mom_fuse");
  }
  if (pin->DoesParameterExist("rad_m1", "vet_mb_agroup")) {
    st.ng = pin->GetInteger("rad_m1", "vet_mb_agroup");
  }
  if (st.ng < 1 || (nrk % st.ng) != 0 || st.ng > vet_nray) {
    VetFatal("<rad_m1>/vet_mb_agroup must be >= 1, divide the number of ranks and not "
             "exceed the number of rays");
  }
  st.gp = me % st.ng;
  {
    const int r0 = me - st.gp;
    st.gsb0 = pm->gids_eachrank[r0];
    st.nsb = 0;
    st.gpr.assign(st.ng, 0);
    st.gso.assign(st.ng, 0);
    st.gsn.assign(st.ng, 0);
    for (int q = 0; q < st.ng; ++q) {
      st.gpr[q] = r0 + q;
      st.gso[q] = pm->gids_eachrank[r0 + q] - st.gsb0;
      st.gsn[q] = pm->nmb_eachrank[r0 + q];
      st.nsb += st.gsn[q];
    }
    st.sown = st.gso[st.gp];
  }
  // the local rays r = p, p + G, ... (all of them without groups); up rays come first
  {
    std::vector<int> rl;
    for (int r = st.gp; r < vet_nray; r += st.ng) {
      rl.push_back(r);
    }
    st.nrl = static_cast<int>(rl.size());
    st.nhl = 0;
    for (int r : rl) {
      if (r < vet_nray/2) ++st.nhl;
    }
    Kokkos::realloc(st.rlist, st.nrl);
    auto rlh = Kokkos::create_mirror_view(st.rlist);
    for (int r = 0; r < st.nrl; ++r) {
      rlh(r) = rl[r];
    }
    Kokkos::deep_copy(st.rlist, rlh);
    auto agh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
    Kokkos::realloc(st.angl, st.nrl, 4);
    auto alh = Kokkos::create_mirror_view(st.angl);
    for (int r = 0; r < st.nrl; ++r) {
      for (int c = 0; c < 4; ++c) {
        alh(r,c) = agh(rl[r],c);
      }
    }
    Kokkos::deep_copy(st.angl, alh);
  }

  // the band: the largest horizontal reach of a ray per layer, in cells
  auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
  const Real dx1 = pmy_pack->pmb->mb_size.h_view(0).dx1;
  const Real dx2 = pmy_pack->pmb->mb_size.h_view(0).dx2;
  const Real dx3 = pmy_pack->pmb->mb_size.h_view(0).dx3;
  Real s2 = 0.0, s3 = 0.0;
  for (int r = 0; r < vet_nray; ++r) {
    s2 = std::max(s2, dx1*std::fabs(ah(r,1))/(std::fabs(ah(r,0))*dx2));
    s3 = std::max(s3, dx1*std::fabs(ah(r,2))/(std::fabs(ah(r,0))*dx3));
  }
  st.w2 = static_cast<int>(std::floor(s2)) + 1;
  st.w3 = thrd ? (static_cast<int>(std::floor(s3)) + 1) : 0;
  if (st.w2 > nx2 || (thrd && st.w3 > nx3)) {
    VetFatal("<rad_m1>/closure = vet_sc: a ray crosses " + std::to_string(st.w2) + " x "
             + std::to_string(st.w3) + " cells per layer, more than a MeshBlock ("
             + std::to_string(nx2) + " x " + std::to_string(nx3)
             + "): use wider MeshBlocks");
  }
  // scaling options (read only when given, so that an input without them writes the
  // restart file it always wrote): vet_mb_halo = K exchanges a band K w wide every K
  // layers and recomputes the overlap redundantly (same arithmetic, same numbers);
  // vet_mb_overlap sweeps the interior of the next layer while the MPI part of the band
  // is in flight; vet_mb_kernel = ray (one thread per ray and cell) | cell (one thread
  // per cell looping over the rays, the original kernel).  None changes a single bit.
  if (pin->DoesParameterExist("rad_m1", "vet_mb_halo")) {
    st.hk = pin->GetInteger("rad_m1", "vet_mb_halo");
  }
  if (pin->DoesParameterExist("rad_m1", "vet_mb_overlap")) {
    st.ovl = pin->GetBoolean("rad_m1", "vet_mb_overlap");
  }
  if (pin->DoesParameterExist("rad_m1", "vet_mb_kernel")) {
    std::string kn = pin->GetString("rad_m1", "vet_mb_kernel");
    if (kn.compare("ray") == 0) {
      st.raypar = true;
    } else if (kn.compare("cell") == 0) {
      st.raypar = false;
    } else {
      VetFatal("<rad_m1>/vet_mb_kernel = '" + kn + "' not implemented (ray | cell)");
    }
  }
  // DEFAULT vet_mb_halo = 3 (tests_m1/runs_3p_fastdefault; exact, 8.2 -> 5.6 ms per
  // sweep on 2 GPUs, tests_m1/runs_3m_vetmb/README_SCALING.md) when the input does not
  // name it, with the ray kernel and no vet_mb_lag, narrowed until the band fits a block.
  if (!pin->DoesParameterExist("rad_m1", "vet_mb_halo") && st.raypar &&
      !(pin->DoesParameterExist("rad_m1", "vet_mb_lag") &&
        pin->GetInteger("rad_m1", "vet_mb_lag") > 0)) {
    st.hk = 3;
    while (st.hk > 1 && (st.hk*st.w2 > nx2 || (thrd && st.hk*st.w3 > nx3))) {--st.hk;}
  }
  if (st.hk < 1) {
    VetFatal("<rad_m1>/vet_mb_halo must be >= 1");
  }
  st.b2 = st.hk*st.w2;
  st.b3 = st.hk*st.w3;
  if (st.b2 > nx2 || (thrd && st.b3 > nx3)) {
    VetFatal("<rad_m1>/vet_mb_halo = " + std::to_string(st.hk) + ": the band ("
             + std::to_string(st.b2) + " x " + std::to_string(st.b3)
             + " cells) is wider than a MeshBlock");
  }
  if (!st.raypar && st.hk > 1) {
    VetFatal("<rad_m1>/vet_mb_halo > 1 needs vet_mb_kernel = ray");
  }
  // vet_mb_mom_batch = B: the moments of B launches are summed by one launch (the
  // intensity planes of the last max(B, 2) launches are kept)
  st.nbat = st.raypar ? 4 : 1;
  if (pin->DoesParameterExist("rad_m1", "vet_mb_mom_batch")) {
    st.nbat = pin->GetInteger("rad_m1", "vet_mb_mom_batch");
  }
  if (pin->DoesParameterExist("rad_m1", "vet_mb_lag") &&
      pin->GetInteger("rad_m1", "vet_mb_lag") > 0) {
    st.nbat = 1;
  }
  if (st.nbat < 1 || (!st.raypar && st.nbat > 1)) {
    VetFatal("<rad_m1>/vet_mb_mom_batch must be >= 1 (1 with vet_mb_kernel = cell)");
  }
  st.nring = st.fuse ? 2 : std::max(st.nbat, 2);
  st.n2w = nx2 + 2*st.b2;
  st.n3w = thrd ? (nx3 + 2*st.b3) : 1;
  const int nx3e = thrd ? nx3 : 1;
  if (thrd) {
    st.nof = 8;
    for (int t = 0, o = 0; t < 9; ++t) {
      if (t == 4) continue;
      st.oj[o] = t%3 - 1;
      st.ok[o] = t/3 - 1;
      ++o;
    }
  } else {
    st.nof = 2;
    st.oj[0] = -1;
    st.oj[1] = 1;
    st.ok[0] = st.ok[1] = 0;
  }
  for (int o = st.nof; o < 8; ++o) {
    st.oj[o] = st.ok[o] = 0;
  }
  st.maxreg = thrd ? std::max(std::max(st.b2*nx3, nx2*st.b3), st.b2*st.b3) : st.b2;
  st.maxcnt = std::max(st.nrl, 2*(nx1 + 2))*st.maxreg;
  // the per-ray band of the intensity planes: cells a ray reads beyond the active area
  // per layer, below and above in x2 (x3): floor of the foot offset (one more cell when
  // the offset is within round-off of an integer, where blocks might floor differently).
  // Local rays only (index rl, angle rlist(rl)).
  {
    auto rlh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), st.rlist);
    Kokkos::realloc(st.rdep, st.nrl, 4);
    auto rdh = Kokkos::create_mirror_view(st.rdep);
    for (int rl = 0; rl < st.nrl; ++rl) {
      const int r = rlh(rl);
      for (int d = 0; d < 2; ++d) {
        rdh(rl,2*d) = 0;
        rdh(rl,2*d + 1) = 0;
        if (d == 1 && !thrd) continue;
        const Real dxp = (d == 0) ? dx2 : dx3;
        const Real sh = -dx1*ah(r,1 + d)/(std::fabs(ah(r,0))*dxp);
        const int o = static_cast<int>(std::floor(sh));
        const int eps = (std::fabs(sh - std::round(sh)) < 1.0e-9) ? 1 : 0;
        rdh(rl,2*d) = std::max(0, -o) + eps;
        rdh(rl,2*d + 1) = std::max(0, o + 1) + eps;
        if (rdh(rl,2*d) > ((d == 0) ? st.w2 : st.w3) ||
            rdh(rl,2*d + 1) > ((d == 0) ? st.w2 : st.w3)) {
          VetFatal("<rad_m1>/closure = vet_sc: per-ray band wider than the band");
        }
      }
    }
    Kokkos::deep_copy(st.rdep, rdh);
    Kokkos::realloc(st.offs, st.nof, st.nrl + 1);
    Kokkos::realloc(st.offr, st.nof, st.nrl + 1);
    auto osh = Kokkos::create_mirror_view(st.offs);
    auto orh = Kokkos::create_mirror_view(st.offr);
    st.cnts.assign(st.nof, 0);
    st.cntr.assign(st.nof, 0);
    for (int o = 0; o < st.nof; ++o) {
      for (int sd = 0; sd < 2; ++sd) {   // sd = 0: sent (direction -o), 1: received
        const int dj = (sd == 0) ? -st.oj[o] : st.oj[o];
        const int dk = (sd == 0) ? -st.ok[o] : st.ok[o];
        int off = 0;
        for (int r = 0; r < st.nrl; ++r) {
          ((sd == 0) ? osh : orh)(o,r) = off;
          int j0, nj, k0, nk;
          VetRayRegion(dj, nx2, st.b2, st.hk, rdh(r,0), rdh(r,1), j0, nj);
          VetRayRegion(dk, nx3e, st.b3, st.hk, rdh(r,2), rdh(r,3), k0, nk);
          off += nj*nk;
        }
        ((sd == 0) ? osh : orh)(o,st.nrl) = off;
        ((sd == 0) ? st.cnts : st.cntr)[o] = off;
        st.maxcnt = std::max(st.maxcnt, off);
      }
    }
    Kokkos::deep_copy(st.offs, osh);
    Kokkos::deep_copy(st.offr, orh);
  }

  // neighbour tables from the global LogicalLocation list, over the SWEEP blocks s =
  // 0..nsb-1 (gid gsb0 + s): a neighbour in my group is local (its sweep index), any
  // other one belongs to the rank with my group position p in its owner's group
  const int nbt = st.nbx1*nbx2*nbx3;
  std::vector<int> gof(nbt, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + l.lx1] = g;
  }
  for (int t = 0; t < nbt; ++t) {
    if (gof[t] < 0) {
      VetFatal("<rad_m1>/closure = vet_sc: incomplete MeshBlock grid");
    }
  }
  st.gid0.assign(nrk, 0);
  for (int rk = 0; rk < nrk; ++rk) {
    st.gid0[rk] = pm->gids_eachrank[rk];
  }
  const int nsb = st.nsb, g0 = st.gsb0, ngr = st.ng;
  const int mygrp = me/ngr;
  // the counterpart of a block's owner: rank, and the block's sweep index there
  auto cpart = [&](const int g, int &rk, int &lid) {
    const int gr = pm->rank_eachmb[g]/ngr;
    rk = gr*ngr + st.gp;
    lid = g - pm->gids_eachrank[gr*ngr];
  };
  Kokkos::realloc(st.lx1, nsb);
  Kokkos::realloc(st.hloc, 8*nsb);
  Kokkos::realloc(st.hself, 8*nsb);
  Kokkos::realloc(st.xloc, 2*nsb);
  st.hgid.assign(8*nsb, -1);
  st.hrank.assign(8*nsb, me);
  st.hlid.assign(8*nsb, -1);
  st.xgid.assign(2*nsb, -1);
  st.xrank.assign(2*nsb, -1);
  st.xlid.assign(2*nsb, -1);
  for (int m = 0; m < nsb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    st.lx1.h_view(m) = l.lx1;
    for (int o = 0; o < 8; ++o) {
      st.hloc.h_view(8*m + o) = -1;
      st.hself.h_view(8*m + o) = 0;
      if (o >= st.nof) continue;
      const int l2 = (l.lx2 + st.oj[o] + nbx2) % nbx2;
      const int l3 = (l.lx3 + st.ok[o] + nbx3) % nbx3;
      const int g = gof[(l3*nbx2 + l2)*st.nbx1 + l.lx1];
      st.hgid[8*m + o] = g;
      st.hself.h_view(8*m + o) = (g == g0 + m) ? 1 : 0;
      if (pm->rank_eachmb[g]/ngr == mygrp) {
        st.hloc.h_view(8*m + o) = g - g0;
      } else {
        cpart(g, st.hrank[8*m + o], st.hlid[8*m + o]);
        st.hmpi = true;
      }
    }
    for (int s = 0; s < 2; ++s) {
      const int lx = l.lx1 + ((s == 0) ? -1 : 1);
      st.xloc.h_view(2*m + s) = -2;
      if (lx < 0 || lx >= st.nbx1) continue;
      const int g = gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + lx];
      st.xgid[2*m + s] = g;
      if (pm->rank_eachmb[g]/ngr == mygrp) {
        st.xrank[2*m + s] = me;
        st.xloc.h_view(2*m + s) = g - g0;
      } else {
        cpart(g, st.xrank[2*m + s], st.xlid[2*m + s]);
        st.xloc.h_view(2*m + s) = -1;
        st.xmpi = true;
      }
    }
  }
  st.lx1.modify_host();
  st.lx1.sync_device();
  st.hloc.modify_host();
  st.hloc.sync_device();
  st.hself.modify_host();
  st.hself.sync_device();
  st.xloc.modify_host();
  st.xloc.sync_device();

  // vet_mb_agg: the message plans.  A piece is one (sweep block, slot) pair with a
  // remote neighbour; the pieces to one rank are ordered by the RECEIVER's (block, slot)
  // = (hlid, nof-1-o) on the sender and (m, o) on the receiver, the same order.
  auto plan = [&](VetMBState::Agg &A, const std::vector<int> &csz,
                  const std::vector<int> &crz) {
    std::vector<std::array<int, 4>> ps, pr;   // (rank, key1, key2, index 8m+o)
    for (int m = 0; m < nsb; ++m) {
      for (int o = 0; o < st.nof; ++o) {
        if (st.hloc.h_view(8*m + o) >= 0) continue;
        ps.push_back({st.hrank[8*m + o], st.hlid[8*m + o], st.nof - 1 - o, 8*m + o});
        pr.push_back({st.hrank[8*m + o], m, o, 8*m + o});
      }
    }
    std::sort(ps.begin(), ps.end());
    std::sort(pr.begin(), pr.end());
    Kokkos::realloc(A.soff, 8*nsb);
    Kokkos::realloc(A.roff, 8*nsb);
    auto soh = Kokkos::create_mirror_view(A.soff);
    auto roh = Kokkos::create_mirror_view(A.roff);
    for (int t = 0; t < 8*nsb; ++t) {
      soh(t) = -1;
      roh(t) = -1;
    }
    A.prk.clear();
    A.sdsp.clear();
    A.scnt.clear();
    A.rdsp.clear();
    A.rcnt.clear();
    for (const auto &x : ps) {
      if (A.prk.empty() || A.prk.back() != x[0]) {
        A.prk.push_back(x[0]);
        A.sdsp.push_back(A.stot);
        A.scnt.push_back(0);
      }
      soh(x[3]) = A.stot;
      const int c = csz[x[3]%8];
      A.stot += c;
      A.scnt.back() += c;
    }
    // the receive side: the same peers (the neighbour relation is symmetric per slot)
    A.rdsp.assign(A.prk.size(), 0);
    A.rcnt.assign(A.prk.size(), 0);
    size_t ip = 0;
    int last = -1;
    for (const auto &x : pr) {
      while (ip < A.prk.size() && A.prk[ip] != x[0]) ++ip;
      if (ip == A.prk.size()) {
        VetFatal("<rad_m1>/vet_mb_agg: asymmetric neighbour ranks");
      }
      if (static_cast<int>(ip) != last) {
        A.rdsp[ip] = A.rtot;
        last = static_cast<int>(ip);
      }
      roh(x[3]) = A.rtot;
      const int c = crz[x[3]%8];
      A.rtot += c;
      A.rcnt[ip] += c;
    }
    Kokkos::deep_copy(A.soff, soh);
    Kokkos::deep_copy(A.roff, roh);
  };
  if (st.agg) {
    std::vector<int> csP(8, 0), crP(8, 0), csB(8, 0), crB(8, 0);
    const int nab = 2*(nx1 + 2);
    for (int o = 0; o < st.nof; ++o) {
      csP[o] = st.cnts[o];
      crP[o] = st.cntr[o];
      const int nj = (st.oj[o] == 0) ? nx2 : st.b2, nk = (st.ok[o] == 0) ? nx3e : st.b3;
      csB[o] = nab*nj*nk;
      crB[o] = nab*nj*nk;
    }
    plan(st.aP, csP, crP);
    plan(st.aB, csB, crB);
    const int ns = std::max(std::max(st.aP.stot, st.aB.stot), 1);
    const int nr = std::max(std::max(st.aP.rtot, st.aB.rtot), 1);
    Kokkos::realloc(st.sfl, ns);
    Kokkos::realloc(st.rfl, nr);
  }

  // the bottom E, F and dx of the sweep blocks
  Kokkos::realloc(st.bwe, nsb, 4, nx3e, nx2);
  Kokkos::deep_copy(st.bwe, 0.0);
  {
    std::vector<Real> dxl(3*nmb), dxg(3*pm->nmb_total);
    for (int m = 0; m < nmb; ++m) {
      dxl[3*m] = pmy_pack->pmb->mb_size.h_view(m).dx1;
      dxl[3*m + 1] = pmy_pack->pmb->mb_size.h_view(m).dx2;
      dxl[3*m + 2] = pmy_pack->pmb->mb_size.h_view(m).dx3;
    }
#if MPI_PARALLEL_ENABLED
    if (ngr > 1) {
      std::vector<int> c3(nrk), d3(nrk);
      for (int rk = 0; rk < nrk; ++rk) {
        c3[rk] = 3*pm->nmb_eachrank[rk];
        d3[rk] = 3*pm->gids_eachrank[rk];
      }
      MPI_Allgatherv(dxl.data(), 3*nmb, MPI_ATHENA_REAL, dxg.data(), c3.data(),
                     d3.data(), MPI_ATHENA_REAL, MPI_COMM_WORLD);
    } else {
      for (int t = 0; t < 3*nmb; ++t) {
        dxg[3*g0 + t] = dxl[t];
      }
    }
#else
    dxg = dxl;
#endif
    Kokkos::realloc(st.sdx, nsb, 3);
    auto sdh = Kokkos::create_mirror_view(st.sdx);
    for (int m = 0; m < nsb; ++m) {
      for (int c = 0; c < 3; ++c) {
        sdh(m,c) = dxg[3*(g0 + m) + c];
      }
    }
    Kokkos::deep_copy(st.sdx, sdh);
  }
  if (ngr > 1) {
    // the group gather of extinction, source and bottom E, F (chunk per block) and
    // the partial moments of the sweep blocks
    const int ncb = nx3e*nx2*nx1, chunk = 2*ncb + 4*nx3e*nx2;
    Kokkos::realloc(st.gsnd, nmb*chunk);
    Kokkos::realloc(st.gbuf, nsb*chunk);
    Kokkos::realloc(st.rrcv, ngr, nmb*10*ncb);
  }
  Kokkos::realloc(st.mom, nsb, 10, nx1, nx3e, nx2);
  Kokkos::deep_copy(st.mom, 0.0);

  Kokkos::realloc(st.csw, nsb, nx1 + 2, 2, st.n3w, st.n2w);
  Kokkos::deep_copy(st.csw, 0.0);
  Kokkos::realloc(st.ipl, nsb, st.nring, st.nrl, st.n3w, st.n2w);
  Kokkos::deep_copy(st.ipl, 0.0);
  if (!st.agg) {
    Kokkos::realloc(st.sbuf, nsb, st.nof, st.maxcnt);
    Kokkos::realloc(st.rbuf, nsb, st.nof, st.maxcnt);
  }
  Kokkos::realloc(vet_ipl, 1, 1, 1, 1, 1);
  if (vet_full) {
    const int n1 = nx1 + 2*indcs.ng;
    const int n2 = nx2 + 2*indcs.ng;
    const int n3 = thrd ? (nx3 + 2*indcs.ng) : 1;
    Kokkos::realloc(st.dsc, nmb, 6, n3, n2, n1);
    Kokkos::deep_copy(st.dsc, 0.0);
    Kokkos::realloc(st.dsc_c, nmb, 6, 1, 1, 1);
    st.pbd = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    st.pbd->InitializeBuffers(6);
  }
  // DIAGNOSTIC: vet_mb_lag = K (read only when given) emulates the BLOCK-JACOBI sweep
  if (pin->DoesParameterExist("rad_m1", "vet_mb_lag")) {
    st.nlag = pin->GetInteger("rad_m1", "vet_mb_lag");
  }
  // PROTOTYPE vet_mb_angles: the angle-decomposed sweep (VetSweepAng)
  if (pin->DoesParameterExist("rad_m1", "vet_mb_angles")) {
    st.ang = pin->GetBoolean("rad_m1", "vet_mb_angles");
  }
  if (st.ang) {
    if (vet_milne || st.nlag > 0) {
      VetFatal("<rad_m1>/vet_mb_angles excludes vet_milne and vet_mb_lag");
    }
    const int np = global_variable::nranks;
    st.nmbt = pm->nmb_total;
    st.nbx2 = nbx2;
    st.n2g = mindcs.nx2;
    st.n3g = thrd ? mindcs.nx3 : 1;
    const int ncb = nx3e*nx2*nx1, chunk = 2*ncb + 4*nx3e*nx2;
    st.gcnt.assign(np, 0);
    st.gdsp.assign(np, 0);
    st.mcnt.assign(np, 0);
    st.mdsp.assign(np, 0);
    for (int rk = 0; rk < np; ++rk) {
      st.gcnt[rk] = pm->nmb_eachrank[rk]*chunk;
      st.gdsp[rk] = pm->gids_eachrank[rk]*chunk;
      st.mcnt[rk] = pm->nmb_eachrank[rk]*10*ncb;
      st.mdsp[rk] = pm->gids_eachrank[rk]*10*ncb;
    }
    Kokkos::realloc(st.gsnd, nmb*chunk);
    Kokkos::realloc(st.gbuf, st.nmbt*chunk);
    Kokkos::realloc(st.gcs, st.nx1g, 2, st.n3g, st.n2g);
    Kokkos::realloc(st.gbot, 4, st.n3g, st.n2g);
    Kokkos::realloc(st.glx, 3*st.nmbt);
    Kokkos::realloc(st.gofd, nbt);
    for (int g = 0; g < st.nmbt; ++g) {
      LogicalLocation &l = pm->lloc_eachmb[g];
      st.glx.h_view(3*g) = l.lx1;
      st.glx.h_view(3*g + 1) = l.lx2;
      st.glx.h_view(3*g + 2) = l.lx3;
    }
    for (int t = 0; t < nbt; ++t) {
      st.gofd.h_view(t) = gof[t];
    }
    st.glx.modify_host();
    st.glx.sync_device();
    st.gofd.modify_host();
    st.gofd.sync_device();
    // dx of every block (the banded sweep reads the block's own)
    std::vector<Real> dxl(3*nmb), dxg(3*st.nmbt);
    for (int m = 0; m < nmb; ++m) {
      dxl[3*m] = pmy_pack->pmb->mb_size.h_view(m).dx1;
      dxl[3*m + 1] = pmy_pack->pmb->mb_size.h_view(m).dx2;
      dxl[3*m + 2] = pmy_pack->pmb->mb_size.h_view(m).dx3;
    }
#if MPI_PARALLEL_ENABLED
    std::vector<int> c3(np), d3(np);
    for (int rk = 0; rk < np; ++rk) {
      c3[rk] = 3*pm->nmb_eachrank[rk];
      d3[rk] = 3*pm->gids_eachrank[rk];
    }
    MPI_Allgatherv(dxl.data(), 3*nmb, MPI_ATHENA_REAL, dxg.data(), c3.data(), d3.data(),
                   MPI_ATHENA_REAL, MPI_COMM_WORLD);
#else
    dxg = dxl;
#endif
    Kokkos::realloc(st.gdx, st.nmbt, 3);
    auto gdh = Kokkos::create_mirror_view(st.gdx);
    for (int g = 0; g < st.nmbt; ++g) {
      for (int c = 0; c < 3; ++c) {
        gdh(g,c) = dxg[3*g + c];
      }
    }
    Kokkos::deep_copy(st.gdx, gdh);
    // my rays: r = p, p + P, ... (both hemispheres on every rank)
    std::vector<int> rl;
    for (int r = me; r < vet_nray; r += np) {
      rl.push_back(r);
    }
    st.nrm = static_cast<int>(rl.size());
    if (st.nrm < 1) {
      VetFatal("<rad_m1>/vet_mb_angles: more ranks than rays");
    }
    Kokkos::realloc(st.rlist, st.nrm);
    auto rlh = Kokkos::create_mirror_view(st.rlist);
    for (int r = 0; r < st.nrm; ++r) {
      rlh(r) = rl[r];
    }
    Kokkos::deep_copy(st.rlist, rlh);
    Kokkos::realloc(st.gip, st.nring, st.nrm, st.n3g, st.n2g);
    Kokkos::realloc(st.gmom, st.nmbt, 10, nx3e, nx2, nx1);
    Kokkos::realloc(st.rrcv, np, nmb*10*ncb);
    Kokkos::deep_copy(st.gmom, 0.0);
    if (me == 0) {
      const double mb = 8.0e-6*(static_cast<double>(st.gsnd.size()) + st.gbuf.size()
                                + st.gcs.size() + st.gbot.size() + st.gip.size()
                                + st.gmom.size() + st.rrcv.size());
      std::cout << "<rad_m1> vet_mb_angles (PROTOTYPE): " << st.nrm << " of " << vet_nray
                << " rays on rank 0, whole mesh per rank: " << mb
                << " MB of angle-sweep buffers on rank 0" << std::endl;
    }
  }
  if (st.nlag > 0 && (st.hk > 1 || st.ovl)) {
    VetFatal("<rad_m1>/vet_mb_lag needs vet_mb_halo = 1 and vet_mb_overlap = false");
  }
  if (st.ng > 1 && (st.ang || st.nlag > 0 || !st.raypar || vet_milne)) {
    VetFatal("<rad_m1>/vet_mb_agroup > 1 excludes vet_mb_angles, vet_mb_lag, "
             "vet_mb_kernel = cell and vet_milne");
  }
  if (st.fuse && (!st.raypar || st.nlag > 0 || st.nrl > 512)) {
    VetFatal("<rad_m1>/vet_mb_mom_fuse needs vet_mb_kernel = ray, no vet_mb_lag and at "
             "most 512 rays per rank");
  }
  if (st.nlag > 0) {
    const int n1 = nx1 + 2*indcs.ng;
    const int n2 = nx2 + 2*indcs.ng;
    const int n3 = thrd ? (nx3 + 2*indcs.ng) : 1;
    Kokkos::realloc(st.hist, nmb, st.nx1g, vet_nray, st.n3w, st.n2w);
    Kokkos::deep_copy(st.hist, 0.0);
    Kokkos::realloc(st.histx, nmb, 2, vet_nray, st.n3w, st.n2w);
    Kokkos::deep_copy(st.histx, 0.0);
    Kokkos::realloc(st.kref, nmb, 7, n3, n2, n1);
    st.e0max.assign(st.nlag, 0.0);
    st.emax.assign(st.nlag, 0.0);
    st.esum.assign(st.nlag, 0.0);
    st.erms.assign(st.nlag, 0.0);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Comm_dup(MPI_COMM_WORLD, &st.comm);
#endif
  if (global_variable::my_rank == 0) {
    const double mb = 8.0e-6*(static_cast<double>(st.csw.size()) + st.ipl.size()
                              + 2.0*st.sbuf.size() + st.dsc.size() + st.sfl.size()
                              + st.rfl.size() + st.gsnd.size() + st.gbuf.size()
                              + st.gmom.size() + st.rrcv.size());
    std::cout << "<rad_m1> vet_sc on " << pm->nmb_total << " MeshBlocks (" << st.nbx1
              << " x " << nbx2 << " x " << nbx3 << "), " << global_variable::nranks
              << " rank(s): exact banded sweep, reach " << st.w2 << " x " << st.w3
              << " cells, band " << st.b2 << " x " << st.b3 << " (exchange every "
              << st.hk << " layer(s)" << (st.ovl ? ", overlapped" : "") << "), "
              << (st.raypar ? "ray" : "cell") << "-parallel kernel, moments every "
              << (st.fuse ? std::string("FUSED in the ray launch")
                          : std::to_string(st.nbat) + " layer(s)") << ", "
              << (st.agg ? "one message per rank" : "one message per slot") << ", "
              << ((st.ng > 1) ? ("HYBRID: rays over groups of " + std::to_string(st.ng)
                                 + " ranks (" + std::to_string(st.nrl)
                                 + " rays and " + std::to_string(st.nsb)
                                 + " sweep blocks on rank 0), ")
                              : std::string(""))
              << st.nof << " horizontal slots, " << mb
              << " MB of sweep buffers on this rank"
              << ((st.nlag > 0) ? (", DIAGNOSTIC block-Jacobi lag with "
                                   + std::to_string(st.nlag) + " sweep(s) per call")
                                : std::string(""))
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetFree

void RadiationM1::VetFree() {
  if (vet_mbs == nullptr) return;
  if (vet_mbs->pbd != nullptr) {
    delete vet_mbs->pbd;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Comm_free(&(vet_mbs->comm));
#endif
  delete vet_mbs;
  vet_mbs = nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetSweepMB
//! \brief the short-characteristics sweep of VetShortChar over SEVERAL MeshBlocks, exact.
//!
//! Launch L (0 .. nx1g-1, the GLOBAL layers) does global layer L for the upward rays and
//! nx1g-1-L for the downward ones in whichever block holds them; the other blocks return
//! at once.  Before the sweep the extinction and source planes are copied into banded
//! planes and their ghost band (x2/x3 neighbours) and x1 face layers are exchanged; after
//! each launch the new intensity plane's ghost band is exchanged (VetBandExchange), and
//! when the sweep crosses an x1 block face (L a multiple of nx1) the upwind plane is
//! handed to the next block (VetX1Move).  Every value a block reads is the very number
//! the neighbour computed, and the arithmetic is that of the single-block kernel term
//! for term: with one MeshBlock (vet_mb_force) the moments are bit-identical to it.
//! A split along x2/x3 costs no serialisation (all blocks of a layer run together);
//! along x1 the stack is swept in turn, which the x1 line solve (implicit_partition =
//! gather) already imposes.
//!
//! GPU cost (tests_m1/runs_3m_vetmb/README_SCALING.md).  A launch holds only one x1
//! layer, so the default kernel runs one thread per (block, ray, cell) (VetRayLaunch)
//! and sums the moments in a second launch (VetMomLaunch, same order, same numbers);
//! vet_mb_kernel = cell keeps the per-cell ray loop.  vet_mb_halo = K exchanges a band
//! K reaches wide after every K-th launch only, and the K-1 launches in between also
//! sweep the part of the band still valid (K-1, K-2, .. reaches) redundantly: the same
//! arithmetic on the same numbers.  vet_mb_overlap sweeps the interior cells of the next
//! launch while the MPI part of the band is in flight.  Every launch passes its functor
//! as a kernel argument (VetFor): no host wait between launches.
//!
//! `lagged` (DIAGNOSTIC vet_mb_lag): after every exchange the ghost band (and the handed
//! x1 plane) is SWAPPED with the one stored by the previous sweep -- each block then
//! sweeps with the neighbours' intensities of the previous sweep, the block-Jacobi
//! alternative, while the self-wrap of a block that spans the whole period stays exact.

void RadiationM1::VetSweepMB(bool lagged) {
  VetMBState &st = *vet_mbs;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nx3 = trans_x3 ? indcs.nx3 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const int nray = vet_nray, nh = vet_nray/2;
  // the sweep blocks and local rays (= nmb, nray, nh without vet_mb_agroup)
  const int nsb = st.nsb, nrl = st.nrl, nhl = st.nhl;
  const bool hyb = (st.ng > 1);
  const int nx1g = st.nx1g, nbx1 = st.nbx1;
  // w2, w3 here = the BAND (plane offset of the active cells); rw2, rw3 = the reach
  const int w2 = st.b2, w3 = st.b3, n2w = st.n2w, n3w = st.n3w;
  const int rw2 = st.w2, rw3 = st.w3, hk = st.hk;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto cs_ = st.csw;
  auto ip_ = st.ipl;
  auto lx_ = st.lx1;
  auto ang_ = vet_ang;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real efl = e_floor;
  const bool milne = vet_milne;
  const Real fmil = iflux_x1min;
  auto bw_ = st.bwe;
  const int ncb = nx3*nx2*nx1, npb = nx3*nx2;
#if MPI_PARALLEL_ENABLED
  const int ngr = st.ng, gp = st.gp;
#endif

  // (a) the moments start from zero (VetShortChar zeroes them; a repeated sweep of
  // vet_mb_lag zeroes them here); the banded extinction/source planes, layer i - is + 1,
  // and the bottom E, F
  if (!hyb) {
    par_for("m1_vet_mb_csw", DevExeSpace(), 0, nmb-1, is, ie, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
      if (lagged) {
        for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {
          vc_(m,n,k,j,i) = 0.0;
        }
      }
      cs_(m,i-is+1,0,k-ks+w3,j-js+w2) = vc_(m,M1_VET_CHX,k,j,i);
      cs_(m,i-is+1,1,k-ks+w3,j-js+w2) = vc_(m,M1_VET_SRC,k,j,i);
      if (i == is) {
        bw_(m,0,k-ks,j-js) = iw_(m,M1_IW_EN,k,j,i);
        bw_(m,1,k-ks,j-js) = iw_(m,M1_IW_F1,k,j,i);
        bw_(m,2,k-ks,j-js) = iw_(m,M1_IW_F2,k,j,i);
        bw_(m,3,k-ks,j-js) = iw_(m,M1_IW_F3,k,j,i);
      }
    });
  } else {
    // vet_mb_agroup: every rank of the group gets the extinction, source and bottom
    // E, F of all the group's blocks (one message per rank pair)
    const int chunk = 2*ncb + 4*npb;
    auto gs_ = st.gsnd;
    auto gb_ = st.gbuf;
    VetFor("m1_vet_hy_pack", nmb*chunk, KOKKOS_LAMBDA(const int t) {
      const int m = t/chunk;
      int q = t - m*chunk;
      if (q < 2*ncb) {
        const int c = q/ncb;
        q -= c*ncb;
        const int k = q/(nx2*nx1);
        q -= k*nx2*nx1;
        const int j = q/nx1, i = q - j*nx1;
        gs_(t) = vc_(m,(c == 0) ? M1_VET_CHX : M1_VET_SRC,ks + k,js + j,is + i);
      } else {
        q -= 2*ncb;
        const int c = q/npb;
        q -= c*npb;
        const int k = q/nx2, j = q - k*nx2;
        const int n = (c == 0) ? M1_IW_EN : ((c == 1) ? M1_IW_F1
                                         : ((c == 2) ? M1_IW_F2 : M1_IW_F3));
        gs_(t) = iw_(m,n,ks + k,js + j,is);
      }
    });
    {
      auto dst = Kokkos::subview(gb_, std::make_pair(st.sown*chunk,
                                                     (st.sown + nmb)*chunk));
      Kokkos::deep_copy(dst, gs_);
    }
    Kokkos::fence();
#if MPI_PARALLEL_ENABLED
    {
      std::vector<MPI_Request> rq;
      for (int q = 0; q < ngr; ++q) {
        if (q == gp) continue;
        rq.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(gb_.data() + static_cast<size_t>(st.gso[q])*chunk, st.gsn[q]*chunk,
                  MPI_ATHENA_REAL, st.gpr[q], 6, st.comm, &rq.back());
      }
      for (int q = 0; q < ngr; ++q) {
        if (q == gp) continue;
        rq.push_back(MPI_REQUEST_NULL);
        MPI_Isend(gs_.data(), nmb*chunk, MPI_ATHENA_REAL, st.gpr[q], 6, st.comm,
                  &rq.back());
      }
      MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    }
#endif
    VetFor("m1_vet_hy_scat", nsb*chunk, KOKKOS_LAMBDA(const int t) {
      const int s = t/chunk;
      int q = t - s*chunk;
      if (q < 2*ncb) {
        const int c = q/ncb;
        q -= c*ncb;
        const int k = q/(nx2*nx1);
        q -= k*nx2*nx1;
        const int j = q/nx1, i = q - j*nx1;
        cs_(s,i + 1,c,k + w3,j + w2) = gb_(t);
      } else {
        q -= 2*ncb;
        const int c = q/npb;
        q -= c*npb;
        const int k = q/nx2, j = q - k*nx2;
        bw_(s,c,k,j) = gb_(t);
      }
    });
    Kokkos::deep_copy(st.mom, 0.0);
  }
  VetX1Move(st, st.csw, 0, -2, 0, nx1, 0, 2, nsb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 1, 1, 0, 2, nsb);
  VetBandExchange(st, st.csw, 0, nx1 + 2, 0, 2, nsb, nx2, nx3);

  auto hs_ = st.hist;
  auto hx_ = st.histx;
  auto hsf_ = st.hself;
  auto xl_ = st.xloc;
  const int nof = st.nof, mreg = st.maxreg;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto lxd_ = st.lx1.d_view;
  const bool ovl = st.ovl && st.hmpi && st.raypar;
  bool pend = false;   // the band exchange of the last plane is in flight
  const int nring = st.nring, nbat = st.nbat;
  auto rd_ = st.rdep;
  // the moment target: mom (sweep block, n, i, k, j), handed to vet_cell by VetMomOut
  auto mt_ = st.mom;
  const int mi0 = 0, mj0 = 0, mk0 = 0, mn0 = 0;
  VetRayK rk;
  rk.nx1 = nx1;
  rk.nx1g = nx1g;
  rk.b2 = w2;
  rk.b3 = w3;
  rk.nx2 = nx2;
  rk.nx3 = nx3;
  rk.thrd = thrd;
  rk.milne = milne;
  rk.fmil = fmil;
  rk.cl = cl;
  rk.efl = efl;
  rk.rd = rd_;
  rk.bw = bw_;
  rk.cs = cs_;
  rk.ip = ip_;
  rk.lx = lxd_;
  rk.ang = st.angl;
  rk.dx = st.sdx;
  for (int l = 0; l < nx1g; ++l) {
    // ring slot of this launch's plane and of the upwind one (l & 1, pw ^ 1 for 2)
    const int pw = l%nring, pr = (l + nring - 1)%nring;
    if (l > 0 && (l % nx1) == 0) {
      // the x1 move copies whole planes: the band must have arrived
      if (pend) {
        VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
        pend = false;
      }
      // the sweep enters the next block of the stack: up rays from below, down from above
      const int q = l/nx1;
      VetX1Move(st, st.ipl, 0, q, pr, pr, 0, nhl, nsb);
      VetX1Move(st, st.ipl, 1, nbx1 - 1 - q, pr, pr, nhl, nrl - nhl, nsb);
      if (lagged) {
        par_for("m1_vet_mb_lagx", DevExeSpace(), 0, nmb-1, 0, nh-1, 0, n3w-1, 0, n2w-1,
        KOKKOS_LAMBDA(const int m, const int r, const int kk, const int jj) {
          for (int s = 0; s < 2; ++s) {
            if (xl_.d_view(2*m + s) == -2) continue;
            if (lx_.d_view(m) != ((s == 0) ? q : (nbx1 - 1 - q))) continue;
            const int rr = r + s*nh;
            const Real t = ip_(m,pr,rr,kk,jj);
            ip_(m,pr,rr,kk,jj) = hx_(m,s,rr,kk,jj);
            hx_(m,s,rr,kk,jj) = t;
          }
        });
      }
    }
    // vet_mb_halo = K: launch l = l0 + tt after the exchange of plane l0 (tt = 1..K; the
    // plane of launch 0 is exchanged too) sweeps (K - tt) reaches beyond the active cells
    const int tt = (l == 0) ? hk : ((l - 1) % hk) + 1;
    if (st.raypar) {
      const int ek = hk - tt;
      const int e2 = ek*rw2, e3 = thrd ? ek*rw3 : 0;
      const int jlo = w2 - e2, nj = nx2 + 2*e2, klo = w3 - e3, nk = nx3 + 2*e3;
      // the interior: cells whose upwind stencil (reach rw) holds no ghost cell
      const int sj0 = w2 + rw2, sj1 = w2 + nx2 - rw2;
      const int sk0 = thrd ? (w3 + rw3) : 0, sk1 = thrd ? (w3 + nx3 - rw3) : 1;
      rk.l = l;
      rk.ek = ek;
      rk.pw = pw;
      rk.pr = pr;
      if (pend && sj1 > sj0 && sk1 > sk0) {
        // vet_mb_overlap: the interior under the MPI part of the band, then the rest
        if (st.fuse) {
          VetRayMomLaunch("m1_vet_mb_raym_in", nsb, nrl, nhl, sj0, sj1 - sj0, sk0,
                          sk1 - sk0, 0, 0, 0, 0, rk, mt_, mi0, mj0, mk0, mn0);
        } else {
          VetRayLaunch("m1_vet_mb_ray_in", nsb, nrl, sj0, sj1 - sj0, sk0, sk1 - sk0, 0,
                       0, 0, 0, rk);
        }
        VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
        pend = false;
        if (st.fuse) {
          VetRayMomLaunch("m1_vet_mb_raym_fr", nsb, nrl, nhl, jlo, nj, klo, nk, sj0,
                          sj1, sk0, sk1, rk, mt_, mi0, mj0, mk0, mn0);
        } else {
          VetRayLaunch("m1_vet_mb_ray_fr", nsb, nrl, jlo, nj, klo, nk, sj0, sj1, sk0,
                       sk1, rk);
        }
      } else {
        if (pend) {
          VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
          pend = false;
        }
        if (st.fuse) {
          VetRayMomLaunch("m1_vet_mb_raym", nsb, nrl, nhl, jlo, nj, klo, nk, 0, 0, 0, 0,
                          rk, mt_, mi0, mj0, mk0, mn0);
        } else {
          VetRayLaunch("m1_vet_mb_ray", nsb, nrl, jlo, nj, klo, nk, 0, 0, 0, 0, rk);
        }
      }
      // the moments of the batch that ends with this launch
      if (!st.fuse && ((l + 1)%nbat == 0 || l == nx1g - 1)) {
        VetMomLaunch(nsb, nrl, l - l%nbat, l, nring, nx1, nx1g, mi0, mj0, mk0, nx2, nx3,
                     w2, w3, mt_, ip_, lxd_, st.angl, st.rlist, nhl, mn0);
      }
    } else {
      par_for("m1_vet_mb_sweep", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        const int lb = lx_.d_view(m);
        const int lu = l - lb*nx1;
        const int ld = (nx1g - 1 - l) - lb*nx1;
        const bool ua = (lu >= 0 && lu < nx1), da = (ld >= 0 && ld < nx1);
        if (!ua && !da) return;
        const Real dx1 = mbsize.d_view(m).dx1;
        const Real dx2 = mbsize.d_view(m).dx2;
        const Real dx3 = mbsize.d_view(m).dx3;
        const int iu = is + lu, id = is + ld;
        const int jc = j - js + w2, kc = k - ks + w3;
        Real acu[10], acd[10];
        for (int n = 0; n < 10; ++n) {
          acu[n] = 0.0;
          acd[n] = 0.0;
        }
        for (int r = 0; r < nray; ++r) {
          const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
          const bool up = (m1 > 0.0);
          if (up ? !ua : !da) continue;
          const int i = up ? iu : id;
          const Real am1 = fabs(m1);
          const Real c0 = vc_(m,M1_VET_CHX,k,j,i);
          const Real s0 = vc_(m,M1_VET_SRC,k,j,i);
          Real iv;
          if (l == 0) {
            // the physical x1 face: VetShortChar, verbatim
            const int iin = up ? (i + 1) : (i - 1);
            const Real su = fmax(1.5*s0 - 0.5*vc_(m,M1_VET_SRC,k,j,iin), 0.0);
            Real ib = 0.0;
            if (up) {
              if (milne) {
                ib = su + 3.0*fmil/cl*m1;
              } else {
                Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
                Real fo = iw_(m,M1_IW_F1,k,j,i)*m1 + iw_(m,M1_IW_F2,k,j,i)*m2
                          + iw_(m,M1_IW_F3,k,j,i)*m3;
                ib = fmax(e + 3.0*fo/cl, 0.0);
              }
            }
            const Real dtau = 0.5*c0*dx1/am1;
            const Real ex = exp(-dtau);
            Real w0, wu;
            if (dtau < 1.0e-3) {
              w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
              wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
            } else {
              const Real g = (1.0 - ex)/dtau;
              w0 = 1.0 - g;
              wu = g - ex;
            }
            iv = ib*ex + wu*su + w0*s0;
          } else {
            // the upwind plane: banded layer i -+ 1 (0 and nx1+1 = the neighbour's face)
            const int li = up ? (i - is) : (i - is + 2);
            const Real sh2 = -dx1*m2/(am1*dx2);
            const Real fl2 = floor(sh2);
            const Real a2 = sh2 - fl2;
            const int o2 = static_cast<int>(fl2);
            const int ja = jc + o2;
            const int jb = ja + 1;
            int ka = kc, kb = kc;
            Real a3 = 0.0;
            if (thrd) {
              const Real sh3 = -dx1*m3/(am1*dx3);
              const Real fl3 = floor(sh3);
              a3 = sh3 - fl3;
              const int o3 = static_cast<int>(fl3);
              ka = kc + o3;
              kb = ka + 1;
            }
            const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
            const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
            const Real iup_v = waa*ip_(m,pr,r,ka,ja) + wab*ip_(m,pr,r,ka,jb)
                               + wba*ip_(m,pr,r,kb,ja) + wbb*ip_(m,pr,r,kb,jb);
            const Real cu = waa*cs_(m,li,0,ka,ja) + wab*cs_(m,li,0,ka,jb)
                            + wba*cs_(m,li,0,kb,ja) + wbb*cs_(m,li,0,kb,jb);
            const Real su = waa*cs_(m,li,1,ka,ja) + wab*cs_(m,li,1,ka,jb)
                            + wba*cs_(m,li,1,kb,ja) + wbb*cs_(m,li,1,kb,jb);
            const Real dtau = 0.5*(cu + c0)*dx1/am1;
            const Real ex = exp(-dtau);
            Real w0, wu;
            if (dtau < 1.0e-3) {
              w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
              wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
            } else {
              const Real g = (1.0 - ex)/dtau;
              w0 = 1.0 - g;
              wu = g - ex;
            }
            iv = iup_v*ex + wu*su + w0*s0;
          }
          iv = fmax(iv, 0.0);
          ip_(m,pw,r,kc,jc) = iv;
          Real *ac = up ? acu : acd;
          const Real a = wr*iv;
          ac[0] += a;
          ac[1] += a*m1*m1;
          ac[2] += a*m2*m2;
          ac[3] += a*m3*m3;
          ac[4] += a*m1*m2;
          ac[5] += a*m1*m3;
          ac[6] += a*m2*m3;
          ac[7] += a*m1;
          ac[8] += a*m2;
          ac[9] += a*m3;
        }
        for (int n = 0; n < 10; ++n) {
          if (ua) {
            vc_(m,M1_VET_J+n,k,j,iu) += acu[n];
          }
          if (da) {
            vc_(m,M1_VET_J+n,k,j,id) += acd[n];
          }
        }
      });
    }
    if (tt == hk && (l < nx1g - 1 || lagged)) {
      VetPlanePost(st, st.ipl, pw, nrl, nsb, nx2, nx3);
      if (ovl && l < nx1g - 1 && ((l + 1) % nx1) != 0) {
        pend = true;
      } else {
        VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
      }
    }
    if (lagged) {
      // block-Jacobi: the band a block reads is the one of the PREVIOUS sweep
      par_for("m1_vet_mb_lagb", DevExeSpace(), 0, nmb-1, 0, nof-1, 0, nray-1, 0, mreg-1,
      KOKKOS_LAMBDA(const int m, const int o, const int r, const int idx) {
        if (hsf_.d_view(8*m + o) != 0) return;
        int j0, nj, k0, nk;
        VetRegion(oj[o], nx2, w2, j0, nj);
        VetRegion(ok[o], nx3, w3, k0, nk);
        if (idx >= nj*nk) return;
        const int kk = k0 + idx/nj, jj = j0 + idx%nj;
        const Real t = ip_(m,pw,r,kk,jj);
        ip_(m,pw,r,kk,jj) = hs_(m,l,r,kk,jj);
        hs_(m,l,r,kk,jj) = t;
      });
    }
  }
  if (!st.raypar) return;   // the per-cell kernel added into vet_cell itself
  if (!hyb) {
    VetMomOut(st.mom, st.rrcv, 0, 1, 0, true, vc_, nmb, nx1, nx2, nx3, is, js, ks);
    return;
  }
  // vet_mb_agroup: the partial moments of the group's blocks go to their owners (one
  // message per rank pair), which sum them in group-rank order into vet_cell
  Kokkos::fence();
  auto gm_ = st.mom;
  auto rv_ = st.rrcv;
  const int mc = nmb*10*ncb;
#if MPI_PARALLEL_ENABLED
  {
    std::vector<MPI_Request> rq;
    for (int q = 0; q < ngr; ++q) {
      if (q == gp) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(rv_.data() + static_cast<size_t>(q)*mc, mc, MPI_ATHENA_REAL, st.gpr[q],
                7, st.comm, &rq.back());
    }
    for (int q = 0; q < ngr; ++q) {
      if (q == gp) continue;
      rq.push_back(MPI_REQUEST_NULL);
      MPI_Isend(gm_.data() + static_cast<size_t>(st.gso[q])*10*ncb, st.gsn[q]*10*ncb,
                MPI_ATHENA_REAL, st.gpr[q], 7, st.comm, &rq.back());
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
  }
#endif
  VetMomOut(gm_, rv_, st.sown, st.ng, st.gp, false, vc_, nmb, nx1, nx2, nx3, is, js, ks);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetMBSweeps
//! \brief the sweep(s) of one VetShortChar call on several MeshBlocks.  Default: ONE
//! exact sweep.  DIAGNOSTIC vet_mb_lag = K: the exact sweep is kept as the reference
//! (J, K_ab), then K block-Jacobi sweeps run, each with the neighbours' band of the
//! sweep before (the last one of the previous call for the first); after each, max and
//! rms over the mesh of max_ab |K_ab/J - (K_ab/J)_exact| are recorded.  The solve reads
//! the K-th lagged sweep, i.e. the run IS the block-Jacobi scheme with K sweeps per step.

void RadiationM1::VetMBSweeps() {
  VetMBState &st = *vet_mbs;
  if (st.ang) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    VetSweepAng(st, pmy_pack->nmb_thispack, indcs.is, indcs.js, indcs.ks, indcs.nx1,
                indcs.nx2, trans_x3 ? indcs.nx3 : 1, trans_x3, vet_milne, iflux_x1min,
                c_light, e_floor, iw, vet_cell, vet_ang);
    return;
  }
  VetSweepMB(false);
  if (st.nlag <= 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  auto vc_ = vet_cell;
  auto kr_ = st.kref;
  par_for("m1_vet_mb_kref", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int n = 0; n < 7; ++n) {
      kr_(m,n,k,j,i) = vc_(m,M1_VET_J+n,k,j,i);
    }
  });
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nc = nmb*nk*nj*ni;
  for (int s = 0; s < st.nlag; ++s) {
    VetSweepMB(true);
    Real emx = 0.0, e2 = 0.0;
    Kokkos::parallel_reduce("m1_vet_mb_err", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
    KOKKOS_LAMBDA(const int idx, Real &smx, Real &s2) {
      const int m = idx/(nk*nj*ni);
      int r = idx - m*(nk*nj*ni);
      const int k = ks + r/(nj*ni);
      r -= (k - ks)*(nj*ni);
      const int j = js + r/ni;
      const int i = is + r - (j - js)*ni;
      const Real ja = vc_(m,M1_VET_J,k,j,i), jb = kr_(m,0,k,j,i);
      Real d = 0.0;
      if (ja > 0.0 && jb > 0.0) {
        for (int n = 1; n < 7; ++n) {
          d = fmax(d, fabs(vc_(m,M1_VET_J+n,k,j,i)/ja - kr_(m,n,k,j,i)/jb));
        }
      } else {
        d = 1.0;
      }
      smx = fmax(smx, d);
      s2 += d*d;
    }, Kokkos::Max<Real>(emx), Kokkos::Sum<Real>(e2));
    Real cnt = static_cast<Real>(nc);
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &emx, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &e2, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const double rms = std::sqrt(e2/cnt);
    if (vet_ncall == 0.0) {
      st.e0max[s] = emx;
    } else {
      st.emax[s] = std::max(st.emax[s], static_cast<double>(emx));
      st.esum[s] += emx;
      st.erms[s] += rms;
    }
  }
  if (vet_ncall > 0.0) {
    st.ncalls += 1.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetShortChar

void RadiationM1::VetShortChar() {
  Kokkos::fence();
  Kokkos::Timer timer;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const int nray = vet_nray;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto ip_ = vet_ipl;
  auto ang_ = vet_ang;
  auto opac_ = opac;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real ar = arad;
  const Real efl = e_floor;
  const bool thermal = (pmy_pack->phydro != nullptr) && coupling && !opac_zero;
  const bool milne = vet_milne;
  const Real fmil = iflux_x1min;

  // (1) extinction and source: cell by cell, or column by column for the Milne
  // diagnostic, which needs tau (the same arithmetic; a column loop per thread leaves
  // the GPU nearly idle)
  if (!milne) {
    auto opd_ = opac;
    auto iwd_ = iw;
    auto vcd_ = vet_cell;
    const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
    VetFor("m1_vet_srcc", (nmb1 + 1)*nk*nj*ni, KOKKOS_LAMBDA(const int t) {
      const int i = is + t%ni;
      int q = t/ni;
      const int j = js + q%nj;
      q /= nj;
      const int k = ks + q%nk, m = q/nk;
      Real chx = fmax(opd_(m,M1_OP_T,k,j,i), 1.0e-300);
      Real e = fmax(iwd_(m,M1_IW_EN,k,j,i), efl);
      Real s = e;
      if (thermal) {
        Real tg = iwd_(m,M1_IW_TP,k,j,i);
        Real t2 = tg*tg;
        Real eth = fmin(opd_(m,M1_OP_P,k,j,i)/chx, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
      vcd_(m,M1_VET_CHX,k,j,i) = chx;
      vcd_(m,M1_VET_SRC,k,j,i) = s;
      for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {
        vcd_(m,n,k,j,i) = 0.0;
      }
    });
  } else {
    par_for("m1_vet_src", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real tau = 0.0;
      for (int i = ie; i >= is; --i) {
        Real chx = fmax(opac_(m,M1_OP_T,k,j,i), 1.0e-300);
        Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
        Real s = e;
        if (milne) {
          // exact grey Milne problem: S = J = 3 H (tau + q(tau)), Hopf q fitted
          Real tc = tau + 0.5*chx*dx1;
          Real qh = 0.710446 - 0.133054*exp(-3.4488*tc);
          s = 3.0*fmil/cl*(tc + qh);
        } else if (thermal) {
          Real t = iw_(m,M1_IW_TP,k,j,i);
          Real t2 = t*t;
          // thermal fraction: the Planck (emission) mean over the extinction, capped at 1
          Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chx, 1.0);
          s = eth*ar*t2*t2 + (1.0 - eth)*e;
        }
        tau += chx*dx1;
        vc_(m,M1_VET_CHX,k,j,i) = chx;
        vc_(m,M1_VET_SRC,k,j,i) = s;
        for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {vc_(m,n,k,j,i) = 0.0;}
      }
    });
  }

  // (2) the sweep: launch l does layer is+l for the upward rays and ie-l for the
  // downward ones, reading the intensity of launch l-1 from the other plane buffer.
  // Several MeshBlocks: the banded sweep over the global layers (VetSweepMB).
  if (vet_mbs != nullptr) {
    VetMBSweeps();
  }
  const int nlaunch = (vet_mbs != nullptr) ? 0 : nx1;   // the banded sweep ran instead
  const int npass = vet_x1per ? vet_x1npass : 1;
  for (int p = 0; p < npass; ++p) {
  const bool wrap = (p > 0);            // layer is (ie) reads the periodic image
  const bool accum = (p == npass - 1);  // only the last pass adds to the moments
  for (int l = 0; l < nlaunch; ++l) {
    const int gl = p*nlaunch + l;
    const int pw = gl & 1, pr = pw ^ 1;
    par_for("m1_vet_sweep", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real dx1 = mbsize.d_view(m).dx1;
      const Real dx2 = mbsize.d_view(m).dx2;
      const Real dx3 = mbsize.d_view(m).dx3;
      const int iu = is + l, id = ie - l;
      Real acu[10], acd[10];
      for (int n = 0; n < 10; ++n) {acu[n] = 0.0; acd[n] = 0.0;}
      for (int r = 0; r < nray; ++r) {
        const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
        const bool up = (m1 > 0.0);
        const int i = up ? iu : id;
        const Real am1 = fabs(m1);
        const Real c0 = vc_(m,M1_VET_CHX,k,j,i);
        const Real s0 = vc_(m,M1_VET_SRC,k,j,i);
        Real iv;
        if (l == 0 && !wrap) {
          // half a cell from the boundary face; S at the face linearly extrapolated
          // from the two boundary cells (exact for the deep linear S of a diffusion
          // regime), and linear in tau along the segment
          const int iin = up ? (i + 1) : (i - 1);
          const Real su = fmax(1.5*s0 - 0.5*vc_(m,M1_VET_SRC,k,j,iin), 0.0);
          Real ib = 0.0;
          if (up) {
            if (milne) {
              // the exact deep Milne intensity 3 H (tau + q + mu) = S + 3 H mu
              ib = su + 3.0*fmil/cl*m1;
            } else {
              // diffusion: eps = E + 3 (F . Omega)/c from the bottom cell's M1 state
              Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
              Real fo = iw_(m,M1_IW_F1,k,j,i)*m1 + iw_(m,M1_IW_F2,k,j,i)*m2
                        + iw_(m,M1_IW_F3,k,j,i)*m3;
              ib = fmax(e + 3.0*fo/cl, 0.0);
            }
          }
          const Real dtau = 0.5*c0*dx1/am1;
          const Real ex = exp(-dtau);
          Real w0, wu;
          if (dtau < 1.0e-3) {
            w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
            wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
          } else {
            const Real g = (1.0 - ex)/dtau;
            w0 = 1.0 - g;
            wu = g - ex;
          }
          iv = ib*ex + wu*su + w0*s0;
        } else {
          const int iup = (l == 0) ? (up ? ie : is) : (up ? (i - 1) : (i + 1));
          // the foot of the characteristic on the upwind plane, in cell units
          const Real sh2 = -dx1*m2/(am1*dx2);
          const Real fl2 = floor(sh2);
          const Real a2 = sh2 - fl2;
          const int o2 = static_cast<int>(fl2);
          const int ja = js + VetWrap(j - js + o2, nx2);
          const int jb = js + VetWrap(j - js + o2 + 1, nx2);
          int ka = k, kb = k;
          Real a3 = 0.0;
          if (thrd) {
            const Real sh3 = -dx1*m3/(am1*dx3);
            const Real fl3 = floor(sh3);
            a3 = sh3 - fl3;
            const int o3 = static_cast<int>(fl3);
            ka = ks + VetWrap(k - ks + o3, nx3);
            kb = ks + VetWrap(k - ks + o3 + 1, nx3);
          }
          const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
          const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
          const Real iup_v = waa*ip_(m,pr,r,ka,ja) + wab*ip_(m,pr,r,ka,jb)
                             + wba*ip_(m,pr,r,kb,ja) + wbb*ip_(m,pr,r,kb,jb);
          const Real cu = waa*vc_(m,M1_VET_CHX,ka,ja,iup)
                          + wab*vc_(m,M1_VET_CHX,ka,jb,iup)
                          + wba*vc_(m,M1_VET_CHX,kb,ja,iup)
                          + wbb*vc_(m,M1_VET_CHX,kb,jb,iup);
          const Real su = waa*vc_(m,M1_VET_SRC,ka,ja,iup)
                          + wab*vc_(m,M1_VET_SRC,ka,jb,iup)
                          + wba*vc_(m,M1_VET_SRC,kb,ja,iup)
                          + wbb*vc_(m,M1_VET_SRC,kb,jb,iup);
          const Real dtau = 0.5*(cu + c0)*dx1/am1;
          const Real ex = exp(-dtau);
          Real w0, wu;
          if (dtau < 1.0e-3) {
            w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
            wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
          } else {
            const Real g = (1.0 - ex)/dtau;
            w0 = 1.0 - g;
            wu = g - ex;
          }
          iv = iup_v*ex + wu*su + w0*s0;
        }
        iv = fmax(iv, 0.0);
        ip_(m,pw,r,k,j) = iv;
        Real *ac = up ? acu : acd;
        const Real a = wr*iv;
        ac[0] += a;
        ac[1] += a*m1*m1;
        ac[2] += a*m2*m2;
        ac[3] += a*m3*m3;
        ac[4] += a*m1*m2;
        ac[5] += a*m1*m3;
        ac[6] += a*m2*m3;
        ac[7] += a*m1;
        ac[8] += a*m2;
        ac[9] += a*m3;
      }
      if (accum) {
        for (int n = 0; n < 10; ++n) {
          vc_(m,M1_VET_J+n,k,j,iu) += acu[n];
          vc_(m,M1_VET_J+n,k,j,id) += acd[n];
        }
      }
    });
  }
  }

  // (3) the uniaxial projection (chi, n) of D = K/J
  const bool axflux = vet_axis_flux;
  par_for("m1_vet_proj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real jj = vc_(m,M1_VET_J,k,j,i);
    Real chi = 1.0/3.0, n1 = 1.0, n2 = 0.0, n3 = 0.0;
    if (jj > 0.0) {
      Real ij = 1.0/jj;
      Real a[3][3];
      a[0][0] = vc_(m,M1_VET_K11,k,j,i)*ij;
      a[1][1] = vc_(m,M1_VET_K11+1,k,j,i)*ij;
      a[2][2] = vc_(m,M1_VET_K11+2,k,j,i)*ij;
      a[0][1] = a[1][0] = vc_(m,M1_VET_K11+3,k,j,i)*ij;
      a[0][2] = a[2][0] = vc_(m,M1_VET_K11+4,k,j,i)*ij;
      a[1][2] = a[2][1] = vc_(m,M1_VET_K11+5,k,j,i)*ij;
      // axis = the formal solution's OWN flux direction, chi = n.D.n (vet_axis = flux)
      // or the principal axis of D (vet_axis = eigen, see the file header)
      Real h1 = vc_(m,M1_VET_H1,k,j,i), h2 = vc_(m,M1_VET_H1+1,k,j,i);
      Real h3 = vc_(m,M1_VET_H1+2,k,j,i);
      Real hm = sqrt(h1*h1 + h2*h2 + h3*h3);
      if (axflux && hm > 1.0e-8*jj) {
        n1 = h1/hm;
        n2 = h2/hm;
        n3 = h3/hm;
        chi = n1*(a[0][0]*n1 + a[0][1]*n2 + a[0][2]*n3)
              + n2*(a[1][0]*n1 + a[1][1]*n2 + a[1][2]*n3)
              + n3*(a[2][0]*n1 + a[2][1]*n2 + a[2][2]*n3);
      } else {
        VetEigMax(a, chi, n1, n2, n3);
      }
      Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
      if (nn > 0.0) {
        n1 /= nn;
        n2 /= nn;
        n3 /= nn;
      } else {
        n1 = 1.0;
        n2 = 0.0;
        n3 = 0.0;
      }
      chi = fmin(fmax(chi, 1.0/3.0), 1.0);
    }
    vc_(m,M1_VET_CHI,k,j,i) = chi;
    vc_(m,M1_VET_N1,k,j,i) = n1;
    vc_(m,M1_VET_N1+1,k,j,i) = n2;
    vc_(m,M1_VET_N1+2,k,j,i) = n3;
  });

  // (4) vet_tensor = full: the GUARDED D = K/J the solve reads in place of (chi, n)
  if (vet_full) {VetFullTensor();}
  Kokkos::fence();
  vet_time += timer.seconds();

  // column dumps (diagnostic, outside the timer)
  if (!vet_dump.empty()) {
    int nc = static_cast<int>(vet_ncall);
    if (nc == 0 || (vet_dump_every > 0 && (nc % vet_dump_every) == 0)) {
      char buf[32];
      std::snprintf(buf, sizeof(buf), ".%06d.txt", nc);
      VetDump(vet_dump + std::string(buf));
    }
  }
  vet_ncall += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetFullTensor
//! \brief vet_tensor = full: D_ab = K_ab/J of the last formal solution, REALIZABILITY
//! GUARDED, into vet_cell(M1_VET_D11 .. +5), active cells and ghosts.
//!
//! By construction K = sum_r w_r I_r Omega_r Omega_r^T with w_r > 0 and I_r >= 0 (the
//! sweep clips I at 0) and sum_r w_r Omega_r Omega_r^T has trace sum_r w_r = J/J, so the
//! raw D is symmetric, has trace 1 to round-off and is positive semidefinite with
//! eigenvalues in [0, 1].  The guard makes that EXACT and enforces a floor:
//!   (1) J <= 0 or a non-positive trace: D = delta/3 (counted);
//!   (2) D /= tr D (a round-off correction, not counted);
//!   (3) eigen-decomposition (cyclic Jacobi); if an eigenvalue lies below vet_eig_min
//!       (default 0) the spectrum is mapped affinely, lam -> emin + (lam' - emin)
//!       (1 - 3 emin)/(sum lam' - 3 emin) with lam' = max(lam, emin), which keeps the
//!       eigenvectors, the trace 1 and every eigenvalue in [emin, 1] (counted).
//! Consequences for the solve: D_dd >= emin >= 0 on every diagonal, so each x1 line row
//! gains nu df (D_11,i + D_11,i) >= 0 on its diagonal and loses nu df D_11,i+-1 <= 0 off
//! it -- the column-diagonally-dominant M-matrix of the uniaxial form, and the diagonal
//! stays >= 1.  The transverse diagonal TDIA >= 0 likewise.
//! Ghosts: periodic wrap in x2/x3 (the SC sweep already requires ONE MeshBlock with
//! periodic x2/x3), copy of the edge cell in x1 (never read at a physical x1 face).
//! Per-cell |D_guarded - D_raw| goes to M1_VET_GD for the run statistics.

void RadiationM1::VetFullTensor() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*(indcs.ng);
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  const bool thrd = trans_x3;
  const Real emin = vet_eig_min;
  auto vc_ = vet_cell;

  par_for("m1_vet_full", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jj = vc_(m,M1_VET_J,k,j,i);
    Real d[6];
    Real raw[6];
    for (int n = 0; n < 6; ++n) {raw[n] = 0.0;}
    Real tr = 0.0;
    if (jj > 0.0) {
      const Real ij = 1.0/jj;
      for (int n = 0; n < 6; ++n) {raw[n] = vc_(m,M1_VET_K11+n,k,j,i)*ij;}
      tr = raw[0] + raw[1] + raw[2];
    }
    bool hit = false;
    if (!(jj > 0.0) || !(tr > 0.0)) {
      d[0] = d[1] = d[2] = 1.0/3.0;
      d[3] = d[4] = d[5] = 0.0;
      hit = true;
    } else {
      const Real itr = 1.0/tr;
      for (int n = 0; n < 6; ++n) {d[n] = raw[n]*itr;}
      Real a[3][3], v[3][3], lam[3];
      a[0][0] = d[0];
      a[1][1] = d[1];
      a[2][2] = d[2];
      a[0][1] = a[1][0] = d[3];
      a[0][2] = a[2][0] = d[4];
      a[1][2] = a[2][1] = d[5];
      VetEigSym(a, lam, v);
      if (lam[0] < emin || lam[1] < emin || lam[2] < emin) {
        Real sl = 0.0;
        for (int p = 0; p < 3; ++p) {
          lam[p] = fmax(lam[p], emin);
          sl += lam[p];
        }
        const Real sc = (1.0 - 3.0*emin)/fmax(sl - 3.0*emin, 1.0e-300);
        for (int p = 0; p < 3; ++p) {lam[p] = emin + (lam[p] - emin)*sc;}
        Real b[3][3];
        for (int p = 0; p < 3; ++p) {
          for (int q = 0; q < 3; ++q) {
            b[p][q] = lam[0]*v[p][0]*v[q][0] + lam[1]*v[p][1]*v[q][1]
                      + lam[2]*v[p][2]*v[q][2];
          }
        }
        d[0] = b[0][0];
        d[1] = b[1][1];
        d[2] = b[2][2];
        d[3] = 0.5*(b[0][1] + b[1][0]);
        d[4] = 0.5*(b[0][2] + b[2][0]);
        d[5] = 0.5*(b[1][2] + b[2][1]);
        hit = true;
      }
    }
    Real dev = 0.0;
    for (int n = 0; n < 6; ++n) {
      vc_(m,M1_VET_D11+n,k,j,i) = d[n];
      dev = fmax(dev, fabs(d[n] - raw[n]));
    }
    // <= -1 marks a cell the guard left alone (trace normalisation only): -1 - dev
    vc_(m,M1_VET_GD,k,j,i) = hit ? dev : (-1.0 - dev);
  });

  // ghosts: periodic in x2 (x3), edge copy in x1.  Several MeshBlocks: the six D_ab go
  // through the ordinary cell-centred exchange (periodic wrap, edges, corners, MPI), and
  // only a PHYSICAL x1 face takes the edge copy.
  if (vet_mbs != nullptr) {
    VetMBState &st = *vet_mbs;
    auto ds_ = st.dsc;
    auto lx_ = st.lx1;   // over the SWEEP blocks: own block m is sweep block sown + m
    const int nbx1 = st.nbx1, so = st.sown;
    par_for("m1_vet_full_dsc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = 0; n < 6; ++n) {
        ds_(m,n,k,j,i) = vc_(m,M1_VET_D11+n,k,j,i);
      }
    });
    MeshBoundaryValuesCC *pb = st.pbd;
    while (pb->InitRecv(6) == TaskStatus::incomplete) {}
    while (pb->PackAndSendCC(st.dsc, st.dsc_c) == TaskStatus::incomplete) {}
    while (pb->RecvAndUnpackCC(st.dsc, st.dsc_c) == TaskStatus::incomplete) {}
    while (pb->ClearSend() == TaskStatus::incomplete) {}
    while (pb->ClearRecv() == TaskStatus::incomplete) {}
    par_for("m1_vet_full_ghmb", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (act) return;
      const int lb = lx_.d_view(so + m);
      int ii = i;
      if (i < is && lb == 0) ii = is;
      if (i > ie && lb == nbx1 - 1) ii = ie;
      for (int n = 0; n < 6; ++n) {
        vc_(m,M1_VET_D11+n,k,j,i) = ds_(m,n,k,j,ii);
      }
    });
  } else {
    par_for("m1_vet_full_gh", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (act) return;
      const int ii = (i < is) ? is : ((i > ie) ? ie : i);
      const int jj = js + VetWrap(j - js, nx2);
      const int kk = thrd ? (ks + VetWrap(k - ks, nx3)) : k;
      for (int n = 0; n < 6; ++n) {
        vc_(m,M1_VET_D11+n,k,j,i) = vc_(m,M1_VET_D11+n,kk,jj,ii);
      }
    });
  }

  // run statistics of the guard (host reductions over the active cells)
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nc = (nmb1 + 1)*nk*nj*ni;
  Real nhit = 0.0, gmax = 0.0, dmin = 1.0;
  Kokkos::parallel_reduce("m1_vet_full_st", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
  KOKKOS_LAMBDA(const int idx, Real &sh, Real &sg, Real &sd) {
    const int m = idx/(nk*nj*ni);
    int r = idx - m*(nk*nj*ni);
    const int k = ks + r/(nj*ni);
    r -= (k - ks)*(nj*ni);
    const int j = js + r/ni;
    const int i = is + r - (j - js)*ni;
    const Real g = vc_(m,M1_VET_GD,k,j,i);
    if (g >= 0.0) {sh += 1.0;}
    sg = fmax(sg, (g >= 0.0) ? g : (-1.0 - g));
    sd = fmin(sd, vc_(m,M1_VET_D11,k,j,i));
  }, Kokkos::Sum<Real>(nhit), Kokkos::Max<Real>(gmax), Kokkos::Min<Real>(dmin));
  vet_nguard += nhit;
  vet_ncell += static_cast<Real>(nc);
  vet_guard_max = std::max(vet_guard_max, gmax);
  vet_d11_min = std::min(vet_d11_min, dmin);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetDump

void RadiationM1::VetDump(const std::string &fname) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  auto vh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_cell);
  auto ih = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
  Real dx1 = pmy_pack->pmb->mb_size.h_view(0).dx1;
  Real x1m = pmy_pack->pmb->mb_size.h_view(0).x1min;
  // the whole (j, i) plane at k = ks: raw K/J, the uniaxial projection and, under
  // vet_tensor = full, the guarded D the solve reads.  Several MeshBlocks: one file per
  // local block, `.plane.g<gid>`, with GLOBAL (i, j) and 17 digits.
  const int nmbd = (vet_mbs != nullptr) ? pmy_pack->nmb_thispack : 1;
  for (int m = 0; m < nmbd; ++m) {
    const int je = indcs.je;
    std::string gname = fname + ".plane";
    int i0 = 0, j0 = 0;
    if (vet_mbs != nullptr) {
      const int gid = pmy_pack->gids + m;
      char gb[16];
      std::snprintf(gb, sizeof(gb), ".g%05d", gid);
      gname += gb;
      LogicalLocation &ll = pmy_pack->pmesh->lloc_eachmb[gid];
      i0 = ll.lx1*indcs.nx1;
      j0 = ll.lx2*indcs.nx2;
    }
    std::ofstream g(gname);
    g << "# closure = vet_sc plane dump, m = " << m << ", k = ks; call "
      << static_cast<int>(vet_ncall) << (vet_full ? " (full)" : " (uniaxial)") << "\n"
      << "# 1 i  2 j  3 J  4-9 K/J 11 22 33 12 13 23  10-12 H/J  13 chi  14-16 n  "
      << "17-22 D_guarded 11 22 33 12 13 23 (0 unless full)  23 E_m1  24 chi_ext\n";
    g << std::setprecision((vet_mbs != nullptr) ? 17 : 10);
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real jj = vh(m,M1_VET_J,ks,j,i);
        Real ij = (jj > 0.0) ? 1.0/jj : 0.0;
        g << i - is + i0 << " " << j - js + j0 << " " << jj;
        for (int n = 0; n < 6; ++n) {g << " " << vh(m,M1_VET_K11+n,ks,j,i)*ij;}
        for (int n = 0; n < 3; ++n) {g << " " << vh(m,M1_VET_H1+n,ks,j,i)*ij;}
        for (int n = 0; n < 4; ++n) {g << " " << vh(m,M1_VET_CHI+n,ks,j,i);}
        for (int n = 0; n < 6; ++n) {
          g << " " << (vet_full ? vh(m,M1_VET_D11+n,ks,j,i) : 0.0);
        }
        g << " " << ih(m,M1_IW_EN,ks,j,i) << " " << vh(m,M1_VET_CHX,ks,j,i) << "\n";
      }
    }
  }
  if (global_variable::my_rank != 0) return;
  std::ofstream f(fname);
  f << "# closure = vet_sc column dump, m = 0, k = ks, j = js; call "
    << static_cast<int>(vet_ncall) << "\n"
    << "# 1 i  2 x1  3 tau_top(centre)  4 chi_ext  5 S  6 J_sc  7 K11/J  8 K22/J  "
    << "9 K33/J  10 K12/J  11 H1/J (=F1/cE)  12 chi_proj  13 n1  14 n2  15 E_m1  "
    << "16 F1/(c E)_m1  17 K/J Milne-exact (tau+q_inf)/(3(tau+q(tau)))\n";
  f << std::setprecision(10);
  Real tau = 0.0;
  std::vector<Real> tc(ie + 1);
  for (int i = ie; i >= is; --i) {
    Real chx = vh(0,M1_VET_CHX,ks,js,i);
    tc[i] = tau + 0.5*chx*dx1;
    tau += chx*dx1;
  }
  for (int i = is; i <= ie; ++i) {
    Real jj = vh(0,M1_VET_J,ks,js,i);
    Real ij = (jj > 0.0) ? 1.0/jj : 0.0;
    Real t = tc[i];
    Real qh = 0.710446 - 0.133054*std::exp(-3.4488*t);
    Real exact = (t + 0.710446)/(3.0*(t + qh));
    Real e = ih(0,M1_IW_EN,ks,js,i);
    f << i - is << " " << x1m + (i - is + 0.5)*dx1 << " " << t << " "
      << vh(0,M1_VET_CHX,ks,js,i) << " " << vh(0,M1_VET_SRC,ks,js,i) << " " << jj
      << " " << vh(0,M1_VET_K11,ks,js,i)*ij << " " << vh(0,M1_VET_K11+1,ks,js,i)*ij
      << " " << vh(0,M1_VET_K11+2,ks,js,i)*ij << " " << vh(0,M1_VET_K11+3,ks,js,i)*ij
      << " " << vh(0,M1_VET_H1,ks,js,i)*ij << " " << vh(0,M1_VET_CHI,ks,js,i)
      << " " << vh(0,M1_VET_N1,ks,js,i) << " " << vh(0,M1_VET_N1+1,ks,js,i)
      << " " << e << " " << ih(0,M1_IW_F1,ks,js,i)/(c_light*e) << " " << exact << "\n";
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetReport

void RadiationM1::VetReport() {
  if (global_variable::my_rank != 0 || vet_ncall <= 0.0) return;
  std::cout << "<rad_m1> vet_sc: formal solutions=" << vet_ncall
            << " SC seconds=" << vet_time << " (" << vet_time/vet_ncall << " per call)"
            << " implicit-solve seconds incl. SC=" << vet_itime
            << " SC/(implicit solve without SC)="
            << vet_time/std::fmax(vet_itime - vet_time, 1.0e-300) << std::endl;
  if (vet_full) {
    std::cout << "<rad_m1> vet_tensor=full: guard triggered in " << vet_nguard << " of "
              << vet_ncell << " cell-calls (eigenvalue floor " << vet_eig_min
              << "); max |D_guarded - K/J|=" << vet_guard_max
              << " (incl. trace normalisation); min D_11=" << vet_d11_min << std::endl;
  }
  if (vet_mbs != nullptr && vet_mbs->nlag > 0) {
    const VetMBState &st = *vet_mbs;
    const double nc = std::max(st.ncalls, 1.0);
    for (int s = 0; s < st.nlag; ++s) {
      std::cout << "<rad_m1> vet_mb_lag sweep " << s + 1 << "/" << st.nlag
                << ": max_ab |K/J - exact| call 0 max=" << st.e0max[s]
                << "; calls >= 1 (" << st.ncalls << "): max=" << st.emax[s]
                << " mean of max=" << st.esum[s]/nc << " mean rms=" << st.erms[s]/nc
                << std::endl;
    }
  }
}

} // namespace radm1
