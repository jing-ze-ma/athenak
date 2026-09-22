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
//! LIMITS (checked at start-up).  Periodic x2 (and x3), uniform mesh, every MeshBlock at
//! least w cells wide in x2 (x3) and 2 cells in x1; vet_milne needs one block along x1.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
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
  int w2 = 0, w3 = 0;          // band widths in x2, x3 (w3 = 0 in 2-D)
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
  DvceArray5D<Real> ipl;       // (nmb, 2, nray, n3w, n2w) I of the last two layers
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
#endif
};

namespace {

// the band region of slot (dj, dk) on the RECEIVER: first plane index and extent
KOKKOS_INLINE_FUNCTION
void VetRegion(const int d, const int nx, const int w, int &r0, int &nr) {
  r0 = (d < 0) ? 0 : ((d == 0) ? w : (w + nx));
  nr = (d == 0) ? nx : w;
}

//----------------------------------------------------------------------------------------
//! \fn VetBandExchange
//! \brief fill the horizontal ghost band of the banded planes A(m, a0..a0+na-1,
//! b0..b0+nb-1, kk, jj) from the (2 or 8) horizontal neighbours' active cells: periodic
//! wrap, a device copy between blocks of this rank, MPI (device buffers, like the
//! boundary machinery) between ranks.  Ghost and source regions never overlap, so the
//! local copy has no race.

void VetBandExchange(VetMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                     const int b0, const int nb, const int nmb, const int nx2,
                     const int nx3) {
  const int nof = st.nof, w2 = st.w2, w3 = st.w3, mreg = st.maxreg;
  const auto oj = st.oj;
  const auto ok = st.ok;
  auto hl_ = st.hloc;
  auto a_ = a;
  par_for("m1_vet_band_loc", DevExeSpace(), 0, nmb-1, 0, nof-1, 0, na*nb-1, 0, mreg-1,
  KOKKOS_LAMBDA(const int m, const int o, const int ab, const int idx) {
    const int n = hl_.d_view(8*m + o);
    if (n < 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(dj, nx2, w2, j0, nj);
    VetRegion(dk, nx3, w3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int jj = j0 + jr, kk = k0 + kr;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    a_(m,aa,bb,kk,jj) = a_(n,aa,bb,kk - dk*nx3,jj - dj*nx2);
  });
#if MPI_PARALLEL_ENABLED
  if (!st.hmpi) return;
  auto sb_ = st.sbuf;
  auto rb_ = st.rbuf;
  // pack: slot o sends my cells that fill the neighbour's ghost region of slot nof-1-o
  // (offset -o), i.e. its ghost index + o * nx
  par_for("m1_vet_band_pack", DevExeSpace(), 0, nmb-1, 0, nof-1, 0, na*nb-1, 0, mreg-1,
  KOKKOS_LAMBDA(const int m, const int o, const int ab, const int idx) {
    if (hl_.d_view(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(-dj, nx2, w2, j0, nj);
    VetRegion(-dk, nx3, w3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    sb_(m,o,ab*nj*nk + idx) = a_(m,aa,bb,k0 + kr + dk*nx3,j0 + jr + dj*nx2);
  });
  Kokkos::fence();
  std::vector<MPI_Request> req;
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      if (st.hrank[8*m + o] == global_variable::my_rank) continue;
      int nj = (st.oj[o] == 0) ? nx2 : w2, nk = (st.ok[o] == 0) ? nx3 : w3;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(rb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, na*nb*nj*nk,
                MPI_ATHENA_REAL, st.hrank[8*m + o], 16*m + o, st.comm, &req.back());
    }
  }
  for (int m = 0; m < nmb; ++m) {
    for (int o = 0; o < nof; ++o) {
      int rk = st.hrank[8*m + o];
      if (rk == global_variable::my_rank) continue;
      int nj = (st.oj[o] == 0) ? nx2 : w2, nk = (st.ok[o] == 0) ? nx3 : w3;
      int lidn = st.hgid[8*m + o] - st.gid0[rk];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(sb_.data() + (static_cast<size_t>(m)*nof + o)*st.maxcnt, na*nb*nj*nk,
                MPI_ATHENA_REAL, rk, 16*lidn + (nof - 1 - o), st.comm, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  par_for("m1_vet_band_unpk", DevExeSpace(), 0, nmb-1, 0, nof-1, 0, na*nb-1, 0, mreg-1,
  KOKKOS_LAMBDA(const int m, const int o, const int ab, const int idx) {
    if (hl_.d_view(8*m + o) >= 0) return;
    const int dj = oj[o], dk = ok[o];
    int j0, nj, k0, nk;
    VetRegion(dj, nx2, w2, j0, nj);
    VetRegion(dk, nx3, w3, k0, nk);
    if (idx >= nj*nk) return;
    const int kr = idx/nj, jr = idx - kr*nj;
    const int aa = a0 + ab/nb, bb = b0 + ab%nb;
    a_(m,aa,bb,k0 + kr,j0 + jr) = rb_(m,o,ab*nj*nk + idx);
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
      const int lidn = st.xgid[2*m + so] - st.gid0[rk];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(plane(m, as), static_cast<int>(cnt), MPI_ATHENA_REAL, rk,
                16*lidn + 8 + s, st.comm, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
#endif
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
  st.n2w = nx2 + 2*st.w2;
  st.n3w = thrd ? (nx3 + 2*st.w3) : 1;
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
  st.maxreg = thrd ? std::max(std::max(st.w2*nx3, nx2*st.w3), st.w2*st.w3) : st.w2;
  st.maxcnt = std::max(vet_nray, 2*(nx1 + 2))*st.maxreg;

  // neighbour tables from the global LogicalLocation list
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
  st.gid0.assign(global_variable::nranks, 0);
  for (int rk = 0; rk < global_variable::nranks; ++rk) {
    st.gid0[rk] = pm->gids_eachrank[rk];
  }
  const int me = global_variable::my_rank;
  const int g0 = pmy_pack->gids;
  Kokkos::realloc(st.lx1, nmb);
  Kokkos::realloc(st.hloc, 8*nmb);
  Kokkos::realloc(st.hself, 8*nmb);
  Kokkos::realloc(st.xloc, 2*nmb);
  st.hgid.assign(8*nmb, -1);
  st.hrank.assign(8*nmb, me);
  st.xgid.assign(2*nmb, -1);
  st.xrank.assign(2*nmb, -1);
  for (int m = 0; m < nmb; ++m) {
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
      st.hrank[8*m + o] = pm->rank_eachmb[g];
      st.hself.h_view(8*m + o) = (g == g0 + m) ? 1 : 0;
      if (pm->rank_eachmb[g] == me) {
        st.hloc.h_view(8*m + o) = g - g0;
      } else {
        st.hmpi = true;
      }
    }
    for (int s = 0; s < 2; ++s) {
      const int lx = l.lx1 + ((s == 0) ? -1 : 1);
      st.xloc.h_view(2*m + s) = -2;
      if (lx < 0 || lx >= st.nbx1) continue;
      const int g = gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + lx];
      st.xgid[2*m + s] = g;
      st.xrank[2*m + s] = pm->rank_eachmb[g];
      if (pm->rank_eachmb[g] == me) {
        st.xloc.h_view(2*m + s) = g - g0;
      } else {
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

  Kokkos::realloc(st.csw, nmb, nx1 + 2, 2, st.n3w, st.n2w);
  Kokkos::deep_copy(st.csw, 0.0);
  Kokkos::realloc(st.ipl, nmb, 2, vet_nray, st.n3w, st.n2w);
  Kokkos::deep_copy(st.ipl, 0.0);
  Kokkos::realloc(st.sbuf, nmb, st.nof, st.maxcnt);
  Kokkos::realloc(st.rbuf, nmb, st.nof, st.maxcnt);
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
                              + 2.0*st.sbuf.size() + st.dsc.size());
    std::cout << "<rad_m1> vet_sc on " << pm->nmb_total << " MeshBlocks (" << st.nbx1
              << " x " << nbx2 << " x " << nbx3 << "), " << global_variable::nranks
              << " rank(s): exact banded sweep, band " << st.w2 << " x " << st.w3
              << " cells, " << st.nof << " horizontal slots, " << mb
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
  const int nx1g = st.nx1g, nbx1 = st.nbx1;
  const int w2 = st.w2, w3 = st.w3, n2w = st.n2w, n3w = st.n3w;
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

  // (a) the moments start from zero (a repeated sweep of vet_mb_lag); the banded
  // extinction/source planes, layer i - is + 1
  par_for("m1_vet_mb_csw", DevExeSpace(), 0, nmb-1, is, ie, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
    for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {
      vc_(m,n,k,j,i) = 0.0;
    }
    cs_(m,i-is+1,0,k-ks+w3,j-js+w2) = vc_(m,M1_VET_CHX,k,j,i);
    cs_(m,i-is+1,1,k-ks+w3,j-js+w2) = vc_(m,M1_VET_SRC,k,j,i);
  });
  VetX1Move(st, st.csw, 0, -2, 0, nx1, 0, 2, nmb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 1, 1, 0, 2, nmb);
  VetBandExchange(st, st.csw, 0, nx1 + 2, 0, 2, nmb, nx2, nx3);

  auto hs_ = st.hist;
  auto hx_ = st.histx;
  auto hsf_ = st.hself;
  auto xl_ = st.xloc;
  const int nof = st.nof, mreg = st.maxreg;
  const auto oj = st.oj;
  const auto ok = st.ok;
  for (int l = 0; l < nx1g; ++l) {
    const int pw = l & 1, pr = pw ^ 1;
    if (l > 0 && (l % nx1) == 0) {
      // the sweep enters the next block of the stack: up rays from below, down from above
      const int q = l/nx1;
      VetX1Move(st, st.ipl, 0, q, pr, pr, 0, nh, nmb);
      VetX1Move(st, st.ipl, 1, nbx1 - 1 - q, pr, pr, nh, nh, nmb);
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
    VetBandExchange(st, st.ipl, pw, 1, 0, nray, nmb, nx2, nx3);
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
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetMBSweeps
//! \brief the sweep(s) of one VetShortChar call on several MeshBlocks.  Default: ONE
//! exact sweep.  DIAGNOSTIC vet_mb_lag = K: the exact sweep is kept as the reference
//! (J, K_ab), then K block-Jacobi sweeps run, each with the neighbours' band of the
//! sweep before (the last one of the previous call for the first); after each, max and
//! rms over the mesh of max_ab |K_ab/J - (K_ab/J)_exact| are recorded.  The solve reads the K-th
//! lagged sweep, i.e. the run IS the block-Jacobi scheme with K sweeps per step.

void RadiationM1::VetMBSweeps() {
  VetMBState &st = *vet_mbs;
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

  // (1) extinction and source, column by column (the Milne diagnostic needs tau)
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

  // (2) the sweep: launch l does layer is+l for the upward rays and ie-l for the
  // downward ones, reading the intensity of launch l-1 from the other plane buffer.
  // Several MeshBlocks: the banded sweep over the global layers (VetSweepMB).
  if (vet_mbs != nullptr) {
    VetMBSweeps();
  }
  const int nlaunch = (vet_mbs != nullptr) ? 0 : nx1;   // the banded sweep ran instead
  for (int l = 0; l < nlaunch; ++l) {
    const int pw = l & 1, pr = pw ^ 1;
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
        if (l == 0) {
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
          const int iup = up ? (i - 1) : (i + 1);
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
      for (int n = 0; n < 10; ++n) {
        vc_(m,M1_VET_J+n,k,j,iu) += acu[n];
        vc_(m,M1_VET_J+n,k,j,id) += acd[n];
      }
    });
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
    auto lx_ = st.lx1;
    const int nbx1 = st.nbx1;
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
      const int lb = lx_.d_view(m);
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
