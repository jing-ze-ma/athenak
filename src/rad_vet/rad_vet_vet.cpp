//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vet.cpp
//! \brief <rad_vet>/closure = vet_sc: a VARIABLE EDDINGTON TENSOR for the multi-D
//! implicit M1 solve from a grey formal solution of the transfer equation by SHORT
//! CHARACTERISTICS on the M1 grid.
//!
//! WHAT IS SOLVED.  Once per hydro step, at the start-of-step state,
//!     mu . grad I = chi (S - I),   chi = rho (kappa_F + kappa_s)  (= opac M1_OP_T)
//! with the grey source (E units, eps = 4 pi I / c, so that E = <eps>)
//!     S = eps_th a T^4 + (1 - eps_th) E,   eps_th = min(kappa_P / chi, 1),
//! i.e. (kappa_a B + kappa_s J)/chi with J taken from the current M1 E (no lambda
//! iteration: E is what the moment solve delivers and it is lagged with the tensor).
//! (vet_scatter = true, opt-in, read only when named: the table split kappa_R = kappa_a
//! + kappa_e and S = VetScatterMix, rad_vet.hpp: Jiang 2021 eq. 6 or true absorption.)
//! Without it the stellar tables carry no separate scattering opacity (kappa_s = 0,
//! kappa_F = Rosseland incl. electron scattering), so the absorption fraction is
//! estimated by the Planck/extinction ratio and CAPPED: S is a convex blend of a T^4
//! and E, never negative.
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
//! S of that cell (vet_bc_diff = face, DEFAULT: E moved
//! to the face by the source's linear extrapolation; cell: the old bottom-cell E).
//! Top (x1max face), downward rays: vacuum.  With vet_bc_bath (default true, m1-h2div)
//! an x1 end whose implicit_bc is marshak gives its entering rays the isotropic incident
//! bath eps = implicit_ebath of that end instead (0: vacuum).
//! (#2, VetBcInit in rad_vet.cpp): the inner end follows vet_bc_inner
//! (shared with vet_gd): emit = the Marshak bath or the black sphere 4 F/c of a flux
//! face, diff = the diffusion intensity above, auto (default) = emit at a Marshak face;
//! a REFLECTING x1max takes the mirrored leaving intensity (vet_bc_reflect = auto, the
//! default; mirror: both ends; open = the old diffusion / vacuum).
//!
//! SEVERAL MESHBLOCKS / MPI RANKS (VetSweepMB).  The same sweep over the GLOBAL layers,
//! EXACT: the planes are banded by w = floor(max dx1 |mu_perp/mu_1|/dx_perp) + 1 cells
//! in x2 (x3), the band of each new intensity plane is exchanged with the horizontal
//! neighbours after every launch (device copy on a rank, MPI between ranks), and at an
//! x1 block face the upwind plane is handed to the next block of the stack.  A split in
//! x2/x3 runs every block of a layer in the same launch (no serialisation); a split in
//! x1 sweeps the stack in turn, as the gathered x1 line solve does; one MPI message per
//! neighbour rank and exchange.  The six D_ab of vet_tensor = full get their ghosts from
//! the ordinary cell-centred exchange.
//!
//! LIMITS (checked at start-up).  Periodic x2 (and x3), uniform mesh, every MeshBlock at
//! least w cells wide in x2 (x3) and 2 cells in x1.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_vetseg.hpp"
#include "rad_vet/rad_vet_vetsrc.hpp"
#include "utils/recipe.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {

// the per-cell sums of one vet_sc launch (VetShortCharSweep): the moments of the up
// (u = 0) and down (u = 1) rays of the cells iu, id and the face sums into them (f2) and
// into their upwind cells (fn), added ray by ray in the ray order (the per-cell and the
// ray-parallel kernel give the same bits)
struct VetScAcc {
  Real ac[2][10] = {}, hp[2][3] = {}, kp[2][M1_VET_NKP] = {}, nl[2] = {};
  Real f2[2][M1_VET_NFS] = {}, fn[2][M1_VET_NFS] = {};
  // ray r of intensity iv in a cell of source s0; g the scalars of its face sums
  // (entry, x2 and x3 interpolation flux, face crossing, exit), each times (1, mu)
  template <class A>
  KOKKOS_INLINE_FUNCTION void Add(const A &ang, const int r, const Real iv, const Real s0,
                                  const Real g[5], const bool hrs, const bool kps,
                                  const bool nls) {
    const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2), wr = ang(r,3);
    const bool up = (m1 > 0.0);
    const int u = up ? 0 : 1;
    auto add4 = [&](Real *fx, const Real q) {
      fx[0] += q;
      fx[1] += q*m1;
      fx[2] += q*m2;
      fx[3] += q*m3;
    };
    add4(f2[u] + 4, g[0]);
    add4(fn[u], g[1]);
    add4(fn[u] + 12, g[2]);
    add4(up ? f2[0] + 8 : fn[1] + 8, g[3]);
    add4(f2[u] + 4, g[4]);
    Real *c = ac[u];
    const Real a = wr*iv;
    c[0] += a;
    c[1] += a*m1*m1;
    c[2] += a*m2*m2;
    c[3] += a*m3*m3;
    c[4] += a*m1*m2;
    c[5] += a*m1*m3;
    c[6] += a*m2*m3;
    c[7] += a*m1;
    c[8] += a*m2;
    c[9] += a*m3;
    if (hrs) {
      hp[u][0] += a*fmax(m1, 0.0);
      hp[u][1] += a*fmax(m2, 0.0);
      hp[u][2] += a*fmax(m3, 0.0);
    }
    if (kps) {
      const Real mu[3] = {m1, m2, m3};
      for (int d = 0; d < 3; ++d) {
        if (!(mu[d] > 0.0)) continue;
        for (int b = 0; b < 3; ++b) kp[u][3*d+b] += a*mu[b]*mu[d];
      }
    }
    if (nls) {nl[u] += wr*fabs(iv - s0);}
  }
  template <class V>
  KOKKOS_INLINE_FUNCTION void Write(const V &vc, const int m, const int k, const int j,
                                    const int iu, const int id, const int is,
                                    const int ie, const bool mrlo, const bool mrhi,
                                    const bool hrs, const bool kps, const bool nls,
                                    const bool f2s) const {
    for (int n = 0; n < 10; ++n) {
      vc(m,M1_VET_J+n,k,j,iu) += ac[0][n];
      vc(m,M1_VET_J+n,k,j,id) += ac[1][n];
    }
    if (hrs) {
      for (int n = 0; n < 3; ++n) {
        vc(m,M1_VET_HP1+n,k,j,iu) += hp[0][n];
        vc(m,M1_VET_HP1+n,k,j,id) += hp[1][n];
      }
    }
    if (kps) {
      for (int n = 0; n < M1_VET_NKP; ++n) {
        vc(m,M1_VET_KP+n,k,j,iu) += kp[0][n];
        vc(m,M1_VET_KP+n,k,j,id) += kp[1][n];
      }
    }
    if (nls) {
      vc(m,M1_VET_NL,k,j,iu) += nl[0];
      vc(m,M1_VET_NL,k,j,id) += nl[1];
    }
    if (f2s) {
      const int ium = (iu > is) ? iu - 1 : (mrlo ? is : ie);
      const int idp = (id < ie) ? id + 1 : (mrhi ? ie : is);
      for (int n = 0; n < M1_VET_NFS; ++n) {
        vc(m,M1_VET_F2+n,k,j,iu) += f2[0][n];
        vc(m,M1_VET_F2+n,k,j,id) += f2[1][n];
        vc(m,M1_VET_F2+n,k,j,ium) += fn[0][n];
        vc(m,M1_VET_F2+n,k,j,idp) += fn[1][n];
      }
    }
  }
};

void VetFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in rad_vet_vet.cpp" << std::endl << msg << std::endl;
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
//! neighbour tables and message buffers.
//!
//! Plane index convention: a banded plane is (kk, jj), kk in [0, n3w), jj in [0, n2w),
//! with the active cells at kk = k - ks + w3, jj = j - js + w2; w2 (w3) is the band a
//! ray reaches per layer, floor(max dx1 |mu_2|/(|mu_1| dx2)) + 1 cells.

struct VetMBState {
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
  DvceArray5D<Real> spl;       // (nmb, 4, NFS, n3w, n2w) clean_consistency: face sums
                               // of the x1 neighbours' cells, for the lower (0) and the
                               // upper (1) one, received from the upper (2), lower (3)
  MeshBoundaryValuesCC *pbd = nullptr;  // vet_tensor = full: ghost exchange of D
  DvceArray5D<Real> dsc, dsc_c;         // its scratch array (nmb, 6, k, j, i) and dummy
#if MPI_PARALLEL_ENABLED
  MPI_Comm comm;
  std::vector<MPI_Request> req;  // the band exchange in flight (VetBandPost/VetBandWait)
#endif
  int pa0 = 0, pna = 0, pb0 = 0, pnb = 0;  // its plane range
  // SWEEP BLOCKS and RAYS: the rank's own blocks and all rays (nsb = nmb, nrl = nray)
  int nsb = 0, gsb0 = 0;       // sweep blocks, gid of the first one
  int nrl = 0, nhl = 0;        // rays, of which up (mu1 > 0; they come first)
  std::vector<int> hlid, xlid;  // receiver's sweep index of a remote slot neighbour
  // ONE message per neighbour rank and exchange; the pieces of all (block, slot) pairs
  // to one rank are packed into one buffer in the receiver's (block, slot) order.  Plan
  // P: the intensity plane, plan B: the csw band.
  struct Agg {
    DvceArray1D<int> soff, roff;     // (8 nsb) piece offset in the flat buffers, -1
    std::vector<int> prk, sdsp, scnt, rdsp, rcnt;   // per peer rank
    int stot = 0, rtot = 0;
  } aP, aB;
  DvceArray1D<Real> sfl, rfl;  // flat send / receive buffers
  int npb = 1;                 // band passes: a band wider than a block takes cells of
                               // the neighbours' band, filled by the pass before
};

namespace {

// the band region of slot (dj, dk) on the RECEIVER: first plane index and extent
KOKKOS_INLINE_FUNCTION
void VetRegion(const int d, const int nx, const int w, int &r0, int &nr) {
  r0 = (d < 0) ? 0 : ((d == 0) ? w : (w + nx));
  nr = (d == 0) ? nx : w;
}

// the pass of the band exchange that fills band cell (kk, jj): one per block width
// between the cell and the active cells (0 for an active cell)
KOKKOS_INLINE_FUNCTION
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
//! like the boundary machinery) between ranks.  A band wider than a block takes st.npb
//! passes: pass p fills the cells p block widths out from the neighbours' cells one
//! width nearer (active cells, or band cells of the pass before), so the local copy,
//! which writes the cells of its pass only, has no race.  VetBandPost launches the local
//! copy and the pack and
//! posts the messages; VetBandWait completes them and unpacks.  Work that does not read
//! the remote part of the band may be launched in between.

void VetBandPost(VetMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
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

void VetBandExchange(VetMBState &st, DvceArray5D<Real> &a, const int a0, const int na,
                     const int b0, const int nb, const int nmb, const int nx2,
                     const int nx3) {
  for (int p = 1; p <= st.npb; ++p) {
    VetBandPost(st, a, a0, na, b0, nb, nmb, nx2, nx3, p);
    VetBandWait(st, a, nmb, nx2, nx3);
  }
}

//----------------------------------------------------------------------------------------
//! \fn VetRayRegion
//! \brief the band region of slot direction d for ONE ray: only the ghost cells that ray
//! reads (dlo below / dhi above the active cells)

KOKKOS_INLINE_FUNCTION
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
//! \brief the band exchange of the intensity plane p (VetBandPost/VetBandWait) cut down
//! to what each ray reads: a ray going up in x2 reads the band below only, and only
//! floor(reach) + 1 cells of it (1 cell for most rays; the band w is set by the most
//! oblique one).  Same values in every cell that is read.

void VetPlanePost(VetMBState &st, DvceArray5D<Real> &a, const int p, const int nray,
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

void VetPlaneWait(VetMBState &st, DvceArray5D<Real> &a, const int nray, const int nmb,
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

} // namespace

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetBathLo
//! \brief the isotropic intensity entering the vet_sc sweep at the
//! inner x1 face when vet_bc_inner resolves to emit (VetBcEmit): the Marshak bath
//! implicit_ebath_x1min, the black sphere 4 F/c of a flux face, nothing at an outflow
//! face (scheme = clean); -1 = the diffusion
//! intensity.  The old rule (vet_bc_bath at a
//! Marshak face) is vet_bc_inner = auto, bitwise.

Real RadiationVet::VetBathLo() const {
  if (!VetBcEmit()) {return -1.0;}
  if (ibc_x1min == M1_IBC_OUTFLOW) return 0.0;
  return (ibc_x1min == M1_IBC_FLUX) ? 4.0*iflux_x1min/c_light : iebath_x1min;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetInit

void RadiationVet::VetInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  if (!trans_on) {
    VetFatal("<rad_vet>/closure = vet_sc needs transport = implicit on a MULTI-D mesh");
  }
  // several MeshBlocks or ranks: the exact banded sweep of VetSweepMB
  const bool mbpath = (pm->nmb_total != 1 || global_variable::nranks != 1);
  // a non-periodic x2 (x3) mesh boundary is a vacuum face of the formal
  // solution (rays entering through it carry I = 0), on the single-MeshBlock sweep; the
  // banded multi-block sweep still needs the periodic wrap
  vet_x2open = (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
                pm->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic);
  vet_x3open = trans_x3 &&
               (pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic ||
                pm->mesh_bcs[BoundaryFace::outer_x3] != BoundaryFlag::periodic);
  // the moment solve's transverse faces must be the same vacuum (implicit_bc_x2min|max,
  // x3min|max = marshak): a reflecting (zero-flux) face has no formal-solution
  // counterpart here
  auto vac = [](const int b) {return b == M1_IBC_MARSHAK || b == M1_IBC_OUTFLOW;};
  if ((vet_x2open && !(vac(ibc_x2min) && vac(ibc_x2max))) ||
      (vet_x3open && !(vac(ibc_x3min) && vac(ibc_x3max)))) {
    VetFatal("<rad_vet>/closure = vet_sc on a non-periodic x2 (x3) mesh needs "
             "implicit_bc_x2min = implicit_bc_x2max (x3min, x3max) = marshak: the "
             "vacuum face of both the formal solution and the moment solve");
  }
  if ((vet_x2open || vet_x3open) && mbpath) {
    VetFatal("<rad_vet>/closure = vet_sc with non-periodic x2 (x3) boundaries (vacuum "
             "faces of the formal solution) needs a single MeshBlock on one rank");
  }
  vet_nmu = pin->GetOrAddInteger("rad_vet", "vet_nmu", 4);
  vet_nphi = pin->GetOrAddInteger("rad_vet", "vet_nphi", 8);
  // m1-h2div (tests_m1/runs_5j_h2div): at an x1 end with implicit_bc = marshak the rays
  // entering the box carry the incident bath, eps = implicit_ebath (0 = the vacuum of a
  // free surface), instead of the diffusion intensity E + 3 F.Omega/c of the end cell's
  // M1 state.  That intensity closes a lagged loop through the tensor: in the thin
  // diffuse-wall shadow the wall cell went to f = 1, E = 1.5-1.6 E_bath (exact 0.5) at
  // every dt and time scheme, and under hesdirk2 its stages failed every step.  Default
  // true; a restart whose file lacks the key keeps the old boundary (false). vet_interp =
  // auto | lin | log | corner, the upwind stencil of the SC segment (VetUpCorner,
  // rad_vet_vetseg.hpp).  auto = corner at a material edge (stencil chi contrast > 4),
  // log elsewhere (as vet_gd, which reads the same key).  DEFAULT auto on a fresh run; a
  // restart whose file lacks the key keeps lin (the old vet_sc form, bitwise).
  {
    const std::string vi = pin->GetOrAddString("rad_vet", "vet_interp",
                               global_variable::restart_run ? "lin" : "auto");
    if (vi.compare("lin") == 0) {
      vet_interp = 0;
    } else if (vi.compare("log") == 0) {
      vet_interp = 1;
    } else if (vi.compare("corner") == 0) {
      vet_interp = 2;
    } else if (vi.compare("auto") == 0) {
      vet_interp = 3;
    } else {
      VetFatal("<rad_vet>/vet_interp = '" + vi + "' is not a choice (auto | lin | log | "
               "corner)");
    }
  }
  // vet_source = relaxed | gas, the VET source of a lagged build
  // (and of predict / rebuild unless time2_vet_source overrides it), as vet_col / vet_gd
  // (vet_col_source there).  DEFAULT gas on vet_sc (fresh and restart: its old source).
  // NOT unified (#3, OPEN): relaxed / blend force S -> B where c dt
  // rho kappa_a >> 1 in a scattering-dominated atmosphere with held gas -- the Cartesian
  // non-LTE tests nltek eps 1e-1..1e-5 degrade (E_top/ex 0.9994 -> 1.079 / 0.858, L1
  // x2) and nlte2 loses 5 % in the top tau bands (vet_gd shows the same with relaxed,
  // which is why its nlte2 deck names vet_col_source = gas); gas on vet_gd costs bsgR
  // Picard 2.8 -> 5.0.  The key and the source function are shared; the defaults stay
  // per closure until the relaxation is weighted by the escape stiffness.
  {
    // (USER 10-20): DEFAULT relaxed_esc under the new default set
    // (utils/recipe.hpp: a fresh run, or <job>/recipe = prod1019; it equals gas on the
    // non-LTE tests and keeps relaxed where the exchange is stiff); a restart whose input
    // lacks the key, or recipe = none, keeps the old vet_sc default gas
    // vet_col_source (the vet_col spelling) has always been inert on vet_sc: a deck that
    // names only it was written for the old vet_sc default and keeps gas (noted), so an
    // old spelling never changes the vet_sc method
    if (!pin->DoesParameterExist("rad_vet", "vet_source")) {
      const bool alias = pin->DoesParameterExist("rad_vet", "vet_col_source");
      const bool nd = recipe::NewDefaults(pin) && !alias;
      pin->GetOrAddString("rad_vet", "vet_source", nd ? "relaxed_esc" : "gas");
      recipe::Log(std::string("<rad_vet>/vet_source = ") + (nd ? "relaxed_esc (default)"
          : (alias ? "gas (the input names only vet_col_source, inert on vet_sc: the old "
                     "vet_sc default; name vet_source to choose)"
          : (recipe::GetMode(pin) == recipe::Mode::none
             ? "gas (recipe none: the old vet_sc default)"
             : "gas (restart input has no key: the old vet_sc default)"))));
    }
    const std::string vs = pin->GetString("rad_vet", "vet_source");
    if (vs.compare("relaxed") == 0) {
      vcol_srelax = true;
    } else if (vs.compare("relaxed_esc") == 0) {
      // (opt-in): relaxed with the cell's photon escape (rad_vet_vetsrc.hpp)
      vcol_srelax = true;
      vcol_sesc = true;
    } else if (vs.compare("gas") == 0) {
      vcol_srelax = false;
    } else {
      VetFatal("<rad_vet>/vet_source = '" + vs +
               "' is not a choice (relaxed | relaxed_esc | gas)");
    }
  }
  // vet_bc_inner / vet_bc_reflect (VetBcInit, rad_vet.cpp; read
  // before vet_bc_bath is recorded, whose old inner meaning a restart keeps).  The inner
  // end now follows vet_bc_inner (VetBathLo); vet_bc_bath keeps the top bath.
  VetBcInit(pin, false);
  vet_bc_bath = pin->GetOrAddBoolean("rad_vet", "vet_bc_bath",
                                     !global_variable::restart_run);
  // the bottom (x1min) diffusion intensity of the formal solution,
  // eps = E + 3 F.Omega/c, takes E of the bottom CELL (cell, the default) or E moved to
  // the bottom FACE, E + (S_f - S_0) with S_f = 1.5 S_0 - 0.5 S_1 the face source the
  // half-cell carry already uses (face).  The cell value is low by dE/dtau dtau_cell/2,
  // a first-order boundary error that the sweep carries ~1 mean free path into the
  // atmosphere (deep atm order 1.5).  DEFAULT face; a restart whose
  // file lacks the key keeps cell (the old boundary); cell reproduces older runs bitwise.
  {
    std::string bd = pin->GetOrAddString("rad_vet", "vet_bc_diff",
                                         global_variable::restart_run ? "cell" : "face");
    if (bd.compare("face") == 0) {
      vet_bc_face = true;
    } else if (bd.compare("cell") == 0) {
      vet_bc_face = false;
    } else {
      VetFatal("<rad_vet>/vet_bc_diff must be cell or face");
    }
  }
  // vet_seg = mean | split | emit | jump (rad_vet_vetseg.hpp).  DEFAULT emit
  // (fresh runs: jump); a restart whose file
  // lacks the key keeps mean (bitwise older runs).
  {
    const char *sgd = global_variable::restart_run ? "mean" : "jump";
    std::string sg = pin->GetOrAddString("rad_vet", "vet_seg", sgd);
    if (sg.compare("split") == 0) {
      vet_seg_split = true;
      vet_seg_emit = false;
    } else if (sg.compare("emit") == 0) {
      vet_seg_split = true;
      vet_seg_emit = true;
    } else if (sg.compare("jump") == 0) {
      vet_seg_split = true;
      vet_seg_emit = true;
      vet_seg_jump = true;
      vet_seg_plateau = pin->GetOrAddBoolean("rad_vet", "vet_seg_plateau", false);
    } else if (sg.compare("mean") == 0) {
      vet_seg_split = false;
    } else {
      VetFatal("<rad_vet>/vet_seg must be mean, split, emit or jump");
    }
  }
  // vet_seg_cross (default false): where the ray crosses a row or layer boundary of the
  // upwind column before the x1 face and that cell is a material edge to the stencil,
  // the last part of the upwind half-segment takes the chi and S of the cell it crosses
  // (VetSegCross; with vet_interp != lin)
  vet_xcr = pin->GetOrAddBoolean("rad_vet", "vet_seg_cross", false);
  {std::string ax = pin->GetOrAddString("rad_vet", "vet_axis", "flux");
  if (ax.compare("flux") == 0) {
    vet_axis_flux = true;
  } else if (ax.compare("eigen") == 0) {
    vet_axis_flux = false;
  } else {
    VetFatal("<rad_vet>/vet_axis = '" + ax + "' not implemented (flux | eigen)");
  }
  }
  // vet_tensor = uniaxial | full.  DEFAULT full since defaults-0927 (2026-09-27;
  // docs/dev/switch_inventory_2026-09-24.md: every vet_sc gate since runs_3j ran full,
  // NC 0, +0.03 ms per SC call).  A restart whose file lacks the key (written when the
  // key was read only when named) keeps uniaxial; the resolved value is recorded.
  vet_full = false;
  {
    std::string vt = pin->GetOrAddString("rad_vet", "vet_tensor",
                         global_variable::restart_run ? "uniaxial" : "full");
    if (vt.compare("full") == 0) {
      vet_full = true;
    } else if (vt.compare("uniaxial") != 0) {
      VetFatal("<rad_vet>/vet_tensor = '" + vt + "' not implemented (uniaxial | full)");
    }
  }
  if (vet_nmu < 1 || vet_nphi < 1) {VetFatal("<rad_vet>/vet_nmu, vet_nphi must be >= 1");}
  // vet_x1_periodic = true (default false; tests_m1/runs_3r_radwave): a PERIODIC x1
  // sweep for a periodic x1 mesh.  The sweep is repeated vet_x1_npass times, each pass
  // entering layer is (ie) from the upward (downward) intensity the previous pass left
  // in layer ie (is); the first pass uses the ordinary boundary intensities and only
  // the last pass accumulates the moments.  The inflow error decays as exp(-tau) along
  // each ray per box crossing.  Single MeshBlock only.  Read only when named.
  vet_x1per = false;
  vet_x1npass = 1;
  // NOT a slip (rad_vet.hpp vet_mirror_npass): vet_x1_npass is the one
  // pass-count key of both multi-pass x1 sweeps, the periodic one (read again below) and
  // the cavity mirrored at both x1 ends (vet_bc_reflect = mirror, default 4, >= 2)
  if (pin->DoesParameterExist("rad_vet", "vet_x1_npass")) {
    vet_mirror_npass = std::max(2, pin->GetInteger("rad_vet", "vet_x1_npass"));
  }
  //   on a periodic x1 mesh the open sweep (upward rays enter layer is with
  // the boundary intensity, downward rays enter ie with 0) is NOT periodic: an x1-uniform
  // state gets an x1-dependent tensor, i.e. a ramp along x1 with a jump at the wrap
  // (hot_blob slab facing x2: 22 % in E, 5 % in rho in 60 cycles).  So a fresh run on a
  // periodic x1 mesh defaults to the periodic sweep (on several MeshBlocks with a fixed
  // number of passes).  A restart whose file lacks the key keeps the open sweep (bitwise
  // continuation).
  const bool x1meshper =
      (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic);
  if (pin->DoesParameterExist("rad_vet", "vet_x1_periodic")) {
    vet_x1per = pin->GetBoolean("rad_vet", "vet_x1_periodic");
  } else if (x1meshper && !global_variable::restart_run) {
    vet_x1per = true;
    pin->SetBoolean("rad_vet", "vet_x1_periodic", vet_x1per);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_sc: periodic x1 mesh, vet_x1_periodic defaulted to true"
                << std::endl;
    }
  }
  if (vet_x1per) {
    if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
      VetFatal("<rad_vet>/vet_x1_periodic needs periodic x1 boundaries");
    }
    const bool np_named = pin->DoesParameterExist("rad_vet", "vet_x1_npass");
    vet_x1npass = pin->GetOrAddInteger("rad_vet", "vet_x1_npass", 64);
    if (vet_x1npass < 1) {VetFatal("<rad_vet>/vet_x1_npass must be >= 1");}
    // vet_x1_ptol: ADAPTIVE passes.  After pass p >= 1 the change of the
    // exit plane (= the next pass's entry), d_p = max|I_p - I_{p-1}| / max I_p over all
    // rays and (k,j), with the contraction q = d_p/d_{p-1}, bounds the error of the entry
    // pass p used by d_p/(1 - q); the sweep stops when that is < vet_x1_ptol, at most
    // vet_x1_npass passes.  Default 1e-6 (the tensor moves by <~ 1e-6 relative, far
    // below the SC closure's own angular/spatial error, 1e-3..1e-2).  0 = exactly
    // vet_x1_npass passes, the default when vet_x1_npass is named without vet_x1_ptol
    // (an existing deck keeps its fixed count) and on a restart whose file lacks the key.
    // From the second build on an adaptive sweep starts from the last exit plane; a fine
    // wave on a large background needs vet_x1_ptol ~ 1e-14.
    const bool fixed = np_named || global_variable::restart_run || mbpath;
    vet_x1ptol = pin->GetOrAddReal("rad_vet", "vet_x1_ptol", fixed ? 0.0 : 1.0e-6);
    if (vet_x1ptol < 0.0) {VetFatal("<rad_vet>/vet_x1_ptol must be >= 0");}
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_sc periodic x1 sweep: "
                << ((vet_x1ptol > 0.0) ? "adaptive, vet_x1_ptol=" : "fixed passes")
                << ((vet_x1ptol > 0.0) ? std::to_string(vet_x1ptol) : std::string())
                << ", vet_x1_npass=" << vet_x1npass
                << ((vet_x1ptol > 0.0) ? " (cap)" : "") << std::endl;
    }
    // vet_x1_periodic_inflow: the first pass's entry intensity.
    //   source   : BOTH directions enter with the face source S (the intensity of an
    //              infinite uniform medium): a uniform state gives I = S on every ray and
    //              H = 0 exactly, whatever npass (default on fresh runs)
    //   boundary : the ordinary boundary intensities, upward E + 3 F.Omega/c, downward 0
    //              (the old form; a restart whose file lacks the key).  One-sided: at
    //              npass = 2 a uniform state far from S (E >> S) gets a spurious net
    //              upward H ~ J (predfix_1011 eq2: H1/J = 0.77, which the thin switch and
    //              the half-range flux turned into F1 = 0.77 c E in one step), and an
    //              existing F feeds back through the upward inflow
    std::string vin = pin->GetOrAddString("rad_vet", "vet_x1_periodic_inflow",
                                          global_variable::restart_run ? "boundary"
                                                                       : "source");
    if (vin.compare("source") == 0) {
      vet_x1src = true;
    } else if (vin.compare("boundary") == 0) {
      vet_x1src = false;
    } else {
      VetFatal("<rad_vet>/vet_x1_periodic_inflow = '" + vin
               + "' is not a choice (source | boundary)");
    }
  }
  vet_nray = 2*vet_nmu*vet_nphi;
  // vet_ring = N (read only when given; 2-D tests) replaces the angle set by
  // the IDORT / Jiang 2021 2-D ring: N in-plane directions theta = (l + 1/2) 2 pi/N from
  // +x1 toward +x2 (N a multiple of 4: no ray along x2), equal weights.  vet_ring_groups
  // = 2 adds the tilted ring sqrt(1/3) (cos, sin) theta, +-sqrt(2/3) in x3 (weights 1/2 :
  // 1/4 : 1/4), which gives f = 1/3 for an isotropic field (Ma et al. 2025 Sect. 3.1).
  int vring = 0, vrgrp = 1;
  if (pin->DoesParameterExist("rad_vet", "vet_ring")) {
    vring = pin->GetInteger("rad_vet", "vet_ring");
    vrgrp = pin->GetOrAddInteger("rad_vet", "vet_ring_groups", 1);
    if (vring < 4 || vring % 4 != 0)
      VetFatal("<rad_vet>/vet_ring must be a multiple of 4");
    if (vrgrp != 1 && vrgrp != 2) VetFatal("<rad_vet>/vet_ring_groups must be 1 or 2");
    vet_nray = vring*((vrgrp == 2) ? 3 : 1);
  }

  // vet_sc_set = gd (read only when named): the vet_gd direction set (the
  // product vet_gd_gl_nmu x vet_gd_gl_nphi [_stagger], or vet_gd_ls / vet_gd_geo) with
  // its pole on x1.  On a Cartesian block this x1-layer sweep is vet_gd's shell sweep
  // (the layers are its shells) and vet_gd's z rotation turns the set about x1, so the
  // three geometries share the set, the rotation and the closure moments
  if (pin->DoesParameterExist("rad_vet", "vet_sc_set") &&
      pin->GetString("rad_vet", "vet_sc_set").compare("gd") == 0) {
    if (vring > 0) VetFatal("<rad_vet>/vet_sc_set = gd replaces vet_ring");
    vsc_gd = true;
    vgd_gl_nmu = pin->GetOrAddInteger("rad_vet", "vet_gd_gl_nmu", 0);
    vgd_gl_nphi = pin->GetOrAddInteger("rad_vet", "vet_gd_gl_nphi", vgd_gl_nphi);
    vgd_gl_stag = pin->GetOrAddBoolean("rad_vet", "vet_gd_gl_stagger", vgd_gl_stag);
    vgd_ls = pin->GetOrAddInteger("rad_vet", "vet_gd_ls", 0);
    // the vet_gd default set, the geodesic vet_gd_geo = 4 (162 directions; HEALPix is not
    // offered: its equator ring would lie in the layers)
    const int geo = pin->GetOrAddInteger("rad_vet", "vet_gd_geo",
                                         (vgd_gl_nmu > 0 || vgd_ls != 0) ? 0 : 4);
    if (geo > 0) vgd_ls = -geo;
    // the period of the rotation about x1: C_nphi of the product set, else a full turn
    vet_nphi = (vgd_gl_nmu > 0) ? vgd_gl_nphi : 1;
  }
  vet_ring_n = vring;
  vet_ring_g = vrgrp;
  VetBuildAngles(0.0);
  if (vsc_gd && vsc_gdmir < 0 && ((vet_bc_mirror && ibc_x1min == M1_IBC_REFLECT) ||
                                  (vet_bc_mirror_hi && ibc_x1max == M1_IBC_REFLECT))) {
    VetFatal("<rad_vet>/vet_sc_set = gd: a mirrored x1 end needs a set with x1 mirror "
             "pairs (e.g. vet_gd_gl_stagger = false)");
  }
  if (vsc_gd && global_variable::my_rank == 0) {
    auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
    Real mn = 1.0;
    for (int r = 0; r < vet_nray; ++r) {mn = std::min(mn, std::fabs(ah(r,0)));}
    std::cout << "<rad_vet> vet_sc_set = gd: " << vet_nray << " rays, min |mu1| " << mn
              << ", x1 mirror pairs " << ((vsc_gdmir > 0) ? "yes" : "no") << std::endl;
  }

  auto &indcs = pm->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pm->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(vet_cell, nmb, M1_VET_NC + M1_VET_NKP +
                  (cln_consist ? M1_VET_NFS : 0), ncells3, ncells2, ncells1);
  Kokkos::deep_copy(vet_cell, 0.0);
  if (mbpath) {
    VetMBInit(pin);   // the banded plane buffers live in vet_mbs
  } else {
    Kokkos::realloc(vet_ipl, nmb, 2, vet_nray, ncells3, ncells2);
    Kokkos::deep_copy(vet_ipl, 0.0);
    // vet_sc_kernel = ray (the default on a device) | cell (on the host): the
    // launch of a layer as one thread per (block, ray, cell) and the per-cell sums, or
    // one thread per cell over the rays; the same bits (read only when named)
    vsc_rayk = !std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value;
    if (pin->DoesParameterExist("rad_vet", "vet_sc_kernel")) {
      const std::string kn = pin->GetString("rad_vet", "vet_sc_kernel");
      if (kn.compare("ray") != 0 && kn.compare("cell") != 0) {
        VetFatal("<rad_vet>/vet_sc_kernel = '" + kn + "' is not a choice (ray | cell)");
      }
      vsc_rayk = (kn.compare("ray") == 0);
    }
    if (vsc_rayk) {Kokkos::realloc(vet_rayg, nmb, 5*vet_nray, ncells3, ncells2);}
    if (vet_x1per && vet_x1ptol > 0.0) {
      Kokkos::realloc(vet_x1ent, nmb, vet_nray, ncells3, ncells2);
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> closure = vet_sc: short characteristics, " << vet_nray
              << (vsc_gd ? " rays (vet_sc_set = gd: the vet_gd set, pole on x1"
                         : " rays (double-Gauss nmu=" + std::to_string(vet_nmu) +
                           " x nphi=" + std::to_string(vet_nphi))
              << "), lagged per hydro step"
              << (vet_full ? ", FULL tensor D = K/J (eigenvalue floor " : "")
              << (vet_full ? std::to_string(vet_eig_min) + ")" : "")
              << ", vet_bc_diff = " << (vet_bc_face ? "face" : "cell")
              << ", vet_seg = "
              << (vet_seg_jump ? "jump" : (vet_seg_emit ? "emit"
                  : (vet_seg_split ? "split" : "mean")))
              << ((vet_seg_jump && !vet_seg_plateau) ? " without the plateau weights"
                                                     : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetBuildAngles
//! \brief the vet_sc angle set (mu1, mu2, mu3, weight), upward hemisphere first, turned
//! by beta (rotation; beta = 0 is the fixed set, bit for bit): double-Gauss
//! azimuths about x1 + beta, ring angles + beta (|beta| < pi/N keeps N/2 rays up).

void RadiationVet::VetBuildAngles(const Real beta) {
  if (vsc_gd) {
    VetGdAngles(beta);
    return;
  }
  const int vring = vet_ring_n, vrgrp = vet_ring_g;
  // the angle set: (mu1, mu2, mu3, weight), upward hemisphere first
  std::vector<double> xg, wg;
  GaussLegendre01(vet_nmu, xg, wg);
  if (vet_ang.extent_int(0) != vet_nray) {Kokkos::realloc(vet_ang, vet_nray, 4);}
  auto ah = Kokkos::create_mirror_view(vet_ang);
  int r = 0;
  for (int h = 0; h < ((vring > 0) ? 0 : 2); ++h) {
    for (int a = 0; a < vet_nmu; ++a) {
      double m1 = (h == 0) ? xg[a] : -xg[a];
      double st = std::sqrt(1.0 - m1*m1);
      for (int l = 0; l < vet_nphi; ++l) {
        double ph = 2.0*M_PI*(l + 0.5)/vet_nphi;
        if (beta != 0.0) {ph += beta;}
        ah(r, 0) = m1;
        ah(r, 1) = st*std::cos(ph);
        ah(r, 2) = st*std::sin(ph);
        ah(r, 3) = 0.5*wg[a]/vet_nphi;
        ++r;
      }
    }
  }
  if (vring > 0) {   // the 2-D ring (up rays first, ray r + nray/2 its
                     // x1 mirror, as the double-Gauss set)
    std::vector<std::array<double, 4>> up, dn;
    for (int g = 0; g < ((vrgrp == 2) ? 3 : 1); ++g) {
      const double sz = (g == 0) ? 0.0 : std::sqrt(2.0/3.0)*((g == 1) ? 1.0 : -1.0);
      const double sp = (g == 0) ? 1.0 : std::sqrt(1.0/3.0);
      const double wg0 = (vrgrp == 1) ? 1.0 : ((g == 0) ? 0.5 : 0.25);
      for (int l = 0; l < vring; ++l) {
        double th = 2.0*M_PI*(l + 0.5)/vring;
        if (beta != 0.0) {th += beta;}
        if (std::cos(th) > 0.0) {
          up.push_back({sp*std::cos(th), sp*std::sin(th), sz, wg0/vring});
        } else {
          dn.push_back({sp*std::cos(th), sp*std::sin(th), sz, wg0/vring});
        }
      }
    }
    const int nhu = static_cast<int>(up.size());
    if (beta != 0.0 && static_cast<int>(dn.size()) != nhu) {
      VetFatal("<rad_vet>/vet_ring: the rotated ring lost its up/down balance");
    }
    for (int q = 0; q < nhu; ++q) {
      for (int c = 0; c < 4; ++c) {
        ah(q, c) = up[q][c];
        // the fixed ring: the x1 mirror of the up ray (bitwise as before); a rotated
        // ring is not x1-symmetric, its own down rays
        ah(q + nhu, c) = (beta != 0.0) ? dn[q][c] : ((c == 0) ? -up[q][c] : up[q][c]);
      }
    }
  }
  Kokkos::deep_copy(vet_ang, ah);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdAngles
//! \brief gdcart-1117 vet_sc_set = gd: the vet_gd set (VetGdInit's builders) with its
//! pole on x1, (mu1, mu2, mu3) = (z, x, y), turned about x1 by beta; the up rays
//! (mu1 > 0) first, then ray r + nray/2 the x1 mirror of ray r where the set has them
//! (vsc_gdmir = 1), else the down rays in the builder's order (vsc_gdmir = -1)

void RadiationVet::VetGdAngles(const Real beta) {
  std::vector<double> d;
  VgdBuildSet(d);
  const int n = static_cast<int>(d.size()/4);
  // the pole: the set's z axis, or where a ray lies in the plane z = 0 (a geodesic set
  // in its built orientation) the axis p of a 4096-point Fibonacci grid that keeps every
  // ray furthest from the layers (max over p of min |n . p|), so that no ray runs along
  // the x1 layers (vet_sc's double-Gauss set: |mu1| >= its first Gauss node)
  double e[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}, best = 0.0;
  bool tilt = false;
  for (int q = 0; q < n; ++q) {tilt = tilt || (std::fabs(d[4*q+2]) < 1.0e-9);}
  for (int g = 0; g < 4096 && tilt; ++g) {
    const double z = 1.0 - (g + 0.5)/4096.0, rz = std::sqrt(1.0 - z*z);
    const double ph = g*M_PI*(3.0 - std::sqrt(5.0));
    const double pg[3] = {rz*std::cos(ph), rz*std::sin(ph), z};
    double mn = 1.0;
    for (int q = 0; q < n; ++q) {
      mn = std::min(mn, std::fabs(d[4*q]*pg[0] + d[4*q+1]*pg[1] + d[4*q+2]*pg[2]));
    }
    if (mn <= best + 1.0e-12) continue;
    best = mn;
    const double a = (std::fabs(pg[0]) < 0.9) ? 1.0 : 0.0, b = 1.0 - a;
    const double h = a*pg[0] + b*pg[1], u[3] = {a - h*pg[0], b - h*pg[1], -h*pg[2]};
    const double un = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    for (int c = 0; c < 3; ++c) {
      e[0][c] = u[c]/un;
      e[2][c] = pg[c];
    }
    for (int c = 0; c < 3; ++c) {
      e[1][c] = e[2][(c+1)%3]*e[0][(c+2)%3] - e[2][(c+2)%3]*e[0][(c+1)%3];
    }
  }
  for (int q = 0; q < n && tilt; ++q) {
    double t[3];
    for (int b = 0; b < 3; ++b) {
      t[b] = e[b][0]*d[4*q] + e[b][1]*d[4*q+1] + e[b][2]*d[4*q+2];
    }
    for (int b = 0; b < 3; ++b) {d[4*q+b] = t[b];}
  }
  const double cb = std::cos(beta), sb = std::sin(beta);
  std::vector<std::array<double, 4>> up, dn;
  for (int q = 0; q < n; ++q) {
    const double *u = &d[4*q];
    const std::array<double, 4> a = {u[2], cb*u[0] - sb*u[1], sb*u[0] + cb*u[1], u[3]};
    if (std::fabs(a[0]) < 1.0e-9) {
      VetFatal("<rad_vet>/vet_sc_set = gd: a ray along the x1 layers");
    }
    ((a[0] > 0.0) ? up : dn).push_back(a);
  }
  const int nh = static_cast<int>(up.size());
  if (2*nh != n) VetFatal("<rad_vet>/vet_sc_set = gd: the set is not half up, half down");
  vsc_gdmir = 1;
  std::vector<std::array<double, 4>> dm;
  for (int q = 0; q < nh && vsc_gdmir > 0; ++q) {
    for (const auto &b : dn) {
      if (std::fabs(b[0] + up[q][0]) + std::fabs(b[1] - up[q][1]) +
          std::fabs(b[2] - up[q][2]) < 1.0e-12) {dm.push_back(b);}
    }
    if (static_cast<int>(dm.size()) != q + 1) vsc_gdmir = -1;
  }
  if (vsc_gdmir > 0) dn = dm;
  vet_nray = n;
  if (vet_ang.extent_int(0) != n) {Kokkos::realloc(vet_ang, n, 4);}
  auto ah = Kokkos::create_mirror_view(vet_ang);
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < 4; ++c) {ah(r, c) = (r < nh) ? up[r][c] : dn[r - nh][c];}
  }
  Kokkos::deep_copy(vet_ang, ah);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotUpdate
//! \brief the vet_sc set of the current rotation block (VetGdBlock, the
//! vet_gd clock: the block of the step's midpoint, the same for every stage and build of
//! a step).  q = VetGdAngleBlk/(pi/2) in [0, 1); double-Gauss: azimuth (2 pi/nphi) q;
//! ring: 0.8 (q - 1/2) (2 pi/N) (rays stay at least 0.1 of a ring step from +-x2).
//! Rebuilt only when the angle changes; the banded sweep's local copy follows.

void RadiationVet::VetRotUpdate() {
  if (!vsc_rot) {return;}
  const Real beta = VetRotBeta(VetGdBlock(pmy_pack->pmesh->ncycle));
  if (beta == vsc_beta) {return;}
  VetRotSet(beta);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetRotBeta
//! \brief the vet_sc set rotation of block blk (VetRotUpdate's formula)

Real RadiationVet::VetRotBeta(const int64_t blk) const {
  // a mirrored x1 end takes ray r + nray/2 as the x1 mirror of ray r: a turned ring has
  // no mirror pairs (nor any ray along a mirror image), so it stays unturned there
  if (vet_ring_n > 0 && !vet_x1per && ((vet_bc_mirror && ibc_x1min == M1_IBC_REFLECT) ||
      (vet_bc_mirror_hi && ibc_x1max == M1_IBC_REFLECT))) {
    return 0.0;
  }
  // the stratified phase with vet_rotate_avg (RotPhase); else as before
  const Real q = vra_on ? static_cast<Real>(RotPhase(blk))
                        : VetGdAngleBlk(blk)/(0.5*M_PI);
  return (vet_ring_n > 0) ? 0.8*(q - 0.5)*2.0*M_PI/vet_ring_n : q*2.0*M_PI/vet_nphi;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotSet
//! \brief the vet_sc set turned by beta

void RadiationVet::VetRotSet(const Real beta) {
  vsc_beta = beta;
  VetBuildAngles(beta);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetScRotAvgSetup
//! \brief the vet_sc fields averaged by vet_rotate_avg (rad_vet_rotavg.cpp): J (ratio)
//! and, normalised by J (kind 3), K_ab, H_a, the half-range sums (implicit_flux_ beam =
//! halfrange) and the C1 sum (thin switch fs*); the source, extinction and eps are
//! set-independent.  The uniaxial projection and the full tensor follow from the averaged
//! K/J in VetShortCharSweep.

void RadiationVet::VetScRotAvgSetup() {
  RotAvgTab &t = vra_tab;
  t.nf = 0;
  vra_a0 = vet_cell;
  vra_a1 = vet_cell;
  vra_c = DvceArray3D<Real>();
  auto add = [&](const int comp, const int kind, const Real lo, const Real hi) {
    t.arr[t.nf] = 0; t.comp[t.nf] = comp; t.kind[t.nf] = kind;
    t.lo[t.nf] = lo; t.hi[t.nf] = hi; ++t.nf;
  };
  t.jf = 0;
  t.jw = vra_eq ? -1 : 0;   // the sets weighted by their J (union of the sets)
  add(M1_VET_J, 2, 0.0, 0.0);
  for (int c = 0; c < 3; ++c) {add(M1_VET_K11 + c, 3, 0.0, 1.0);}
  for (int c = 3; c < 6; ++c) {add(M1_VET_K11 + c, 3, -1.0, 1.0);}
  for (int c = 0; c < 3; ++c) {add(M1_VET_H1 + c, 3, -1.0, 1.0);}
  if (impl_beam_hr) {
    for (int c = 0; c < 3; ++c) {add(M1_VET_HP1 + c, 3, 0.0, 1.0);}
  }
  if (thsw_c1) {add(M1_VET_NL, 3, 0.0, 1.0e300);}
  t.farr = -1;
  t.fcomp = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetMBInit
//! \brief several MeshBlocks / ranks: band widths, neighbour tables, buffers.  Uniform
//! mesh (the implicit solve already fatals on SMR/AMR), periodic x2 (x3), every block at
//! least as wide as the band.

void RadiationVet::VetMBInit(ParameterInput *pin) {
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
    VetFatal("<rad_vet>/closure = vet_sc on several MeshBlocks needs a uniform mesh");
  }
  st.nbx1 = mindcs.nx1/nx1;
  const int nbx2 = mindcs.nx2/nx2;
  const int nbx3 = thrd ? (mindcs.nx3/nx3) : 1;
  st.nx1g = mindcs.nx1;
  if (nx1 < 2) {
    VetFatal("<rad_vet>/closure = vet_sc needs <meshblock>/nx1 >= 2");
  }
  const int me = global_variable::my_rank;
  // the sweep blocks: the rank's own; the local rays: all, the up rays (mu1 > 0) first
  st.gsb0 = pm->gids_eachrank[me];
  st.nsb = nmb;
  st.nrl = vet_nray;
  st.nhl = vet_nray/2;

  // the band: the largest horizontal reach of a ray per layer, in cells
  auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
  const Real dx1 = MeshBlock0Dx(0), dx2 = MeshBlock0Dx(1), dx3 = MeshBlock0Dx(2);
  Real s2 = 0.0, s3 = 0.0;
  for (int r = 0; r < vet_nray; ++r) {
    s2 = std::max(s2, dx1*std::fabs(ah(r,1))/(std::fabs(ah(r,0))*dx2));
    s3 = std::max(s3, dx1*std::fabs(ah(r,2))/(std::fabs(ah(r,0))*dx3));
  }
  if (vsc_rot) {
    // the band covers every rotated set: the transverse reach
    // |mu_perp|/|mu1| bounds both axes (double-Gauss, rotation about x1); the ring's
    // steepest ray, 0.1 of a ring step from x2, at either extreme of the rotation
    for (int r = 0; r < vet_nray; ++r) {
      const Real tp = std::sqrt(ah(r,1)*ah(r,1) + ah(r,2)*ah(r,2))/std::fabs(ah(r,0));
      if (vet_ring_n == 0) {
        s2 = std::max(s2, dx1*tp/dx2);
        s3 = std::max(s3, dx1*tp/dx3);
      }
    }
    if (vet_ring_n > 0) {
      const Real tm = std::tan(0.5*M_PI - 0.1*2.0*M_PI/vet_ring_n);
      s2 = std::max(s2, dx1*tm/dx2);
      s3 = std::max(s3, dx1*tm/dx3);
    }
  }
  // scheme = clean: one cell more for the slopes of the second-order foot (VetPlmFoot);
  // a reach within round-off of an integer one more (the per-ray band below)
  const int wp = 1;
  auto wband = [wp](const Real s) {
    return static_cast<int>(std::floor(s)) + 1 + wp +
           ((std::fabs(s - std::round(s)) < 1.0e-9) ? 1 : 0);
  };
  st.w2 = wband(s2);
  st.w3 = thrd ? wband(s3) : 0;
  st.b2 = st.w2;
  st.b3 = st.w3;
  st.npb = std::max((st.b2 + nx2 - 1)/nx2, thrd ? (st.b3 + nx3 - 1)/nx3 : 1);
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
  // the per-ray band of the intensity planes: cells a ray reads beyond the active area
  // per layer, below and above in x2 (x3): floor of the foot offset (one more cell when
  // the offset is within round-off of an integer, where blocks might floor differently).
  // with the vet_sc set rotation (vsc_rot, VetRotSet) a ray slot's foot
  // offset moves with the angle while this band is built once, so a turned ray read
  // ghost cells that were never exchanged (stale intensity: the banded sweep depended on
  // the MeshBlock layout, J up to 6 % of its maximum on oneo).  The band of each slot is
  // now the union over the WHOLE rotation range (sampled finely, padded by the largest
  // offset change between samples, capped at the band w, which already bounds every
  // rotated set).  Without rotation: as before (bitwise); one MeshBlock: not used.
  {
    Kokkos::realloc(st.rdep, st.nrl, 4);
    auto rdh = Kokkos::create_mirror_view(st.rdep);
    for (int rl = 0; rl < st.nrl; ++rl) {
      const int r = rl;
      for (int d = 0; d < 2; ++d) {
        rdh(rl,2*d) = 0;
        rdh(rl,2*d + 1) = 0;
        if (d == 1 && !thrd) continue;
        const Real dxp = (d == 0) ? dx2 : dx3;
        const Real sh = -dx1*ah(r,1 + d)/(std::fabs(ah(r,0))*dxp);
        const int o = static_cast<int>(std::floor(sh));
        const int eps = (std::fabs(sh - std::round(sh)) < 1.0e-9) ? 1 : 0;
        rdh(rl,2*d) = std::max(0, -o) + eps + wp;
        rdh(rl,2*d + 1) = std::max(0, o + 1) + eps + wp;
        if (rdh(rl,2*d) > ((d == 0) ? st.w2 : st.w3) ||
            rdh(rl,2*d + 1) > ((d == 0) ? st.w2 : st.w3)) {
          VetFatal("<rad_vet>/closure = vet_sc: per-ray band wider than the band");
        }
      }
    }
    if (vsc_rot) {
      // every set VetRotBeta can give: double-Gauss beta = q 2 pi/nphi, ring beta =
      // 0.8 (q - 1/2) 2 pi/N, q in [0, 1]
      const int nsam = 512;
      // a real copy (create_mirror, not _view: on a host build the mirror view IS
      // vet_ang)
      auto keep = Kokkos::create_mirror(Kokkos::HostSpace(), vet_ang);
      Kokkos::deep_copy(keep, vet_ang);
      std::vector<Real> smn(2*st.nrl, 1.0e300), smx(2*st.nrl, -1.0e300);
      std::vector<Real> sprev(2*st.nrl, 0.0), dmax(2*st.nrl, 0.0);
      for (int qs = 0; qs <= nsam; ++qs) {
        const Real q = static_cast<Real>(qs)/nsam;
        const Real beta = (vet_ring_n > 0) ? 0.8*(q - 0.5)*2.0*M_PI/vet_ring_n
                                           : q*2.0*M_PI/vet_nphi;
        VetBuildAngles(beta);
        auto bh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
        for (int rl = 0; rl < st.nrl; ++rl) {
          const int r = rl;
          for (int d = 0; d < 2; ++d) {
            if (d == 1 && !thrd) continue;
            const Real dxp = (d == 0) ? dx2 : dx3;
            const Real sh = -dx1*bh(r,1 + d)/(std::fabs(bh(r,0))*dxp);
            const int x = 2*rl + d;
            smn[x] = std::min(smn[x], sh);
            smx[x] = std::max(smx[x], sh);
            if (qs > 0) {dmax[x] = std::max(dmax[x], std::fabs(sh - sprev[x]));}
            sprev[x] = sh;
          }
        }
      }
      Kokkos::deep_copy(vet_ang, keep);   // the current set, bit for bit
      for (int rl = 0; rl < st.nrl; ++rl) {
        for (int d = 0; d < 2; ++d) {
          if (d == 1 && !thrd) continue;
          const int x = 2*rl + d;
          const int w = (d == 0) ? st.w2 : st.w3;
          const int olo = static_cast<int>(std::floor(smn[x] - dmax[x] - 1.0e-9));
          const int ohi = static_cast<int>(std::floor(smx[x] + dmax[x] + 1.0e-9));
          rdh(rl,2*d) = std::min(w, std::max(rdh(rl,2*d), std::max(0, -olo) + wp));
          rdh(rl,2*d + 1) = std::min(w, std::max(rdh(rl,2*d + 1),
                                                 std::max(0, ohi + 1) + wp));
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
          VetRayRegion(dj, nx2, st.b2, rdh(r,0), rdh(r,1), j0, nj);
          VetRayRegion(dk, nx3e, st.b3, rdh(r,2), rdh(r,3), k0, nk);
          off += nj*nk;
        }
        ((sd == 0) ? osh : orh)(o,st.nrl) = off;
        ((sd == 0) ? st.cnts : st.cntr)[o] = off;
      }
    }
    Kokkos::deep_copy(st.offs, osh);
    Kokkos::deep_copy(st.offr, orh);
  }

  // neighbour tables from the global LogicalLocation list, over the SWEEP blocks s =
  // 0..nsb-1 (gid gsb0 + s): a neighbour on my rank is local (its sweep index)
  const int nbt = st.nbx1*nbx2*nbx3;
  std::vector<int> gof(nbt, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + l.lx1] = g;
  }
  for (int t = 0; t < nbt; ++t) {
    if (gof[t] < 0) {
      VetFatal("<rad_vet>/closure = vet_sc: incomplete MeshBlock grid");
    }
  }
  const int nsb = st.nsb, g0 = st.gsb0;
  // the owner of a block: rank, and the block's sweep index there
  auto cpart = [&](const int g, int &rk, int &lid) {
    rk = pm->rank_eachmb[g];
    lid = g - pm->gids_eachrank[rk];
  };
  Kokkos::realloc(st.lx1, nsb);
  Kokkos::realloc(st.hloc, 8*nsb);
  Kokkos::realloc(st.xloc, 2*nsb);
  st.hrank.assign(8*nsb, me);
  st.hlid.assign(8*nsb, -1);
  st.xrank.assign(2*nsb, -1);
  st.xlid.assign(2*nsb, -1);
  for (int m = 0; m < nsb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    st.lx1.h_view(m) = l.lx1;
    for (int o = 0; o < 8; ++o) {
      st.hloc.h_view(8*m + o) = -1;
      if (o >= st.nof) continue;
      const int l2 = (l.lx2 + st.oj[o] + nbx2) % nbx2;
      const int l3 = (l.lx3 + st.ok[o] + nbx3) % nbx3;
      const int g = gof[(l3*nbx2 + l2)*st.nbx1 + l.lx1];
      if (pm->rank_eachmb[g] == me) {
        st.hloc.h_view(8*m + o) = g - g0;
      } else {
        cpart(g, st.hrank[8*m + o], st.hlid[8*m + o]);
        st.hmpi = true;
      }
    }
    for (int s = 0; s < 2; ++s) {
      int lx = l.lx1 + ((s == 0) ? -1 : 1);
      st.xloc.h_view(2*m + s) = -2;
      if (vet_x1per) lx = (lx + st.nbx1) % st.nbx1;   // the periodic image
      if (lx < 0 || lx >= st.nbx1) continue;
      const int g = gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + lx];
      if (pm->rank_eachmb[g] == me) {
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
        VetFatal("<rad_vet>/vet_mb_agg: asymmetric neighbour ranks");
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
  {
    std::vector<int> csP(8, 0), crP(8, 0), csB(8, 0), crB(8, 0);
    const int nab = 2*(nx1 + 4);
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

  Kokkos::realloc(st.csw, nsb, nx1 + 4, 2, st.n3w, st.n2w);
  Kokkos::deep_copy(st.csw, 0.0);
  Kokkos::realloc(st.ipl, nsb, st.nring, st.nrl, st.n3w, st.n2w);
  Kokkos::deep_copy(st.ipl, 0.0);
  if (cln_consist) Kokkos::realloc(st.spl, nsb, 4, M1_VET_NFS, st.n3w, st.n2w);
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
    // the D-tensor halo follows <mesh>/cs_seam_blockend (as pbval_u;
    // VetInit runs from ImplicitInit, before pbval_u exists, so the key is read here)
    st.pbd->cs_blockend_ext_ = CsSeamBlockendExtend(pmy_pack, pin);
    // the six components of the symmetric tensor D: no slot is a tangential VECTOR
    // member (the default rule would treat slots 2,3 as one).  A tensor transform at
    // poles/seams is a later stage; rad_vet is a fatal on sp/cs until then.
    st.pbd->SetVectorPairs(6, {});
  }
  if (vet_x1per && vet_x1ptol > 0.0) {
    VetFatal("<rad_vet>/vet_x1_periodic on several MeshBlocks needs a fixed "
             "vet_x1_npass");
  }
#if MPI_PARALLEL_ENABLED
  MPI_Comm_dup(MPI_COMM_WORLD, &st.comm);
#endif
  if (global_variable::my_rank == 0) {
    const double mb = 8.0e-6*(static_cast<double>(st.csw.size()) + st.ipl.size()
                              + st.dsc.size() + st.sfl.size() + st.rfl.size());
    std::cout << "<rad_vet> vet_sc on " << pm->nmb_total << " MeshBlocks (" << st.nbx1
              << " x " << nbx2 << " x " << nbx3 << "), " << global_variable::nranks
              << " rank(s): exact banded sweep, reach " << st.w2 << " x " << st.w3
              << " cells, cell-parallel kernel, one message per rank, " << st.nof
              << " horizontal slots, " << mb << " MB of sweep buffers on this rank"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetFree

void RadiationVet::VetFree() {
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
//! \fn void RadiationVet::VetSweepMB
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
//! for term: on one MeshBlock the moments are bit-identical to it.
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

void RadiationVet::VetSweepMB() {
  // vet_source_noesrc for vet_sc (as vet_col, fixbundle-1009 F1): the
  // source and the boundary intensity read the PHYSICAL E^n = u0 under be, not the
  // solve's old vector E^n + dt (chat/c) esrc (M1_IW_EN)
  const bool noes_ = VsrcNoes();
  auto esv_ = noes_ ? u0 : iw;
  const int esc_ = noes_ ? M1_E : M1_IW_EN;
  VetMBState &st = *vet_mbs;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nx3 = trans_x3 ? indcs.nx3 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const bool hrs_ = impl_beam_hr;   // the half-range sums (cell kernel)
  const bool nls_ = thsw_c1;   // sum w |I - S| (M1_VET_NL)
  const int nray = vet_nray;
  const int nsb = st.nsb, nrl = st.nrl, nhl = st.nhl;
  const int nx1g = st.nx1g, nbx1 = st.nbx1;
  // w2, w3 here = the BAND (plane offset of the active cells)
  const int w2 = st.b2, w3 = st.b3;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto cs_ = st.csw;
  auto ip_ = st.ipl;
  auto lx_ = st.lx1;
  auto ang_ = vet_ang;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real efl = e_floor;
  const Real bblo = VetBathLo();   // vet_bc_inner
  const Real bbhi = (vet_bc_bath && ibc_x1max == M1_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool bfc = vet_bc_face;
  const bool spl = vet_seg_split;
  const bool xcr = vet_xcr;
  const int spe = VetSegMode();
  const int vim = vet_interp;   // (#6)
  const bool plm_ = true;

  // (a) the banded extinction/source planes, layer i - is + 2 (VetShortChar zeroes the
  // moments)
  par_for("vet_vet_mb_csw", DevExeSpace(), 0, nmb-1, is, ie, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
    cs_(m,i-is+2,0,k-ks+w3,j-js+w2) = vc_(m,M1_VET_CHX,k,j,i);
    cs_(m,i-is+2,1,k-ks+w3,j-js+w2) = vc_(m,M1_VET_SRC,k,j,i);
  });
  // the neighbours' two face layers, the nearer first (with nx1 = 1 the farther one is
  // the neighbour's own face layer from ITS neighbour)
  VetX1Move(st, st.csw, 0, -2, 1, nx1 + 1, 0, 2, nsb);
  VetX1Move(st, st.csw, 0, -2, 0, nx1, 0, 2, nsb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 2, 2, 0, 2, nsb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 3, 3, 0, 2, nsb);
  VetBandExchange(st, st.csw, 0, nx1 + 4, 0, 2, nsb, nx2, nx3);

  const int nring = st.nring;
  // vet_x1_periodic: npass sweeps, the first layer of a pass entering with the exit plane
  // of the pass before (the periodic image); the moments are those of the last pass.
  // clean_consistency: the face sums of a neighbour's cell go to spl (VetX1Move after)
  const int npass = vet_x1per ? vet_x1npass : 1, ngl = npass*nx1g;
  const bool x1src = vet_x1per && vet_x1src;
  const bool kps_ = true, f2s_ = cln_consist, x1per = vet_x1per;
  const int bmn_ = ibm_nr;
  const int bmr0 = ibm_ray[0], bmr1 = ibm_ray[1], bmr2 = ibm_ray[2], bmr3 = ibm_ray[3];
  auto bmi_ = ibm_i;
  auto sp_ = st.spl;
  if (f2s_) Kokkos::deep_copy(st.spl, 0.0);
  for (int gl = 0; gl < ngl; ++gl) {
    const int l = gl%nx1g;
    const bool wrap = (gl >= nx1g), accum = (gl >= ngl - nx1g);
    // ring slot of this launch's plane and of the upwind one (l & 1, pw ^ 1 for 2)
    const int pw = gl%nring, pr = (gl + nring - 1)%nring;
    if (gl > 0 && (l % nx1) == 0) {
      // the sweep enters the next block of the stack: up rays from below, down from above
      const int q = l/nx1;
      VetX1Move(st, st.ipl, 0, q, pr, pr, 0, nhl, nsb);
      VetX1Move(st, st.ipl, 1, nbx1 - 1 - q, pr, pr, nhl, nrl - nhl, nsb);
    }
    par_for("vet_vet_mb_sweep", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
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
      Real hpu[3] = {0.0, 0.0, 0.0}, hpd[3] = {0.0, 0.0, 0.0};
      Real nlu = 0.0, nld = 0.0;
      Real kpu[M1_VET_NKP], kpd[M1_VET_NKP];
      for (int n = 0; n < M1_VET_NKP; ++n) {kpu[n] = 0.0; kpd[n] = 0.0;}
      // face fluxes into the cells iu, id (F1, FX) and their upwind cells (F2, FX, F3)
      Real f2u[M1_VET_NFS] = {0.0}, f2d[M1_VET_NFS] = {0.0};
      Real fnu[M1_VET_NFS] = {0.0}, fnd[M1_VET_NFS] = {0.0};
      for (int r = 0; r < nray; ++r) {
        const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
        const bool up = (m1 > 0.0);
        if (up ? !ua : !da) continue;
        const int i = up ? iu : id;
        const Real am1 = fabs(m1);
        const Real c0 = vc_(m,M1_VET_CHX,k,j,i);
        const Real s0 = vc_(m,M1_VET_SRC,k,j,i);
        Real iv;
        auto add4 = [&](Real *fx, const Real g) {
          fx[0] += g;
          fx[1] += g*m1;
          fx[2] += g*m2;
          fx[3] += g*m3;
        };
        Real ifc = 0.0;   // clean_consistency: the intensity at the face crossing
        if (l == 0 && !wrap) {
          // the physical x1 face: VetShortChar, verbatim
          const int iin = up ? (i + 1) : (i - 1);
          const Real s1 = vc_(m,M1_VET_SRC,k,j,iin);
          const Real su = fmax(1.5*s0 - 0.5*s1, 0.0);
          Real ib = up ? 0.0 : bbhi;
          if (x1src) {
            ib = su;
          } else if (up) {
            if (bblo >= 0.0) {
              ib = bblo;
            } else {
              Real e = esv_(m,esc_,k,j,i);
              if (bfc) e += 0.5*(s0 - s1);   // vet_bc_diff = face
              e = fmax(e, efl);
              Real fo = iw_(m,M1_IW_F1,k,j,i)*m1 + iw_(m,M1_IW_F2,k,j,i)*m2
                        + iw_(m,M1_IW_F3,k,j,i)*m3;
              ib = fmax(e + 3.0*fo/cl, 0.0);
            }
            const int q = (r == bmr0) ? 0 : ((r == bmr1) ? 1
                          : ((r == bmr2) ? 2 : ((r == bmr3) ? 3 : -1)));
            if (q >= 0 && q < bmn_) {ib += bmi_(m,q,k-ks,j-js);}
          }
          if (f2s_) add4((up ? f2u : f2d) + 4, wr*m1*ib);
          const Real dtau = 0.5*c0*dx1/am1;
          Real ex, w0, wu;
          VetSegW(dtau, ex, w0, wu);
          iv = ib*ex + wu*su + w0*s0;
        } else {
          // the upwind plane: banded layer i -+ 1 (1 and nx1+2 = the neighbour's face)
          const int li = up ? (i - is + 1) : (i - is + 3);
          const Real sh2 = -dx1*m2/(am1*dx2);
          const Real fl2 = floor(sh2);
          const Real a2 = sh2 - fl2;
          const int o2 = static_cast<int>(fl2);
          const int ja = jc + o2;
          const int jb = ja + 1;
          int ka = kc, kb = kc, o3 = 0;
          Real a3 = 0.0;
          const Real sh3 = thrd ? -dx1*m3/(am1*dx3) : 0.0;
          if (thrd) {
            const Real fl3 = floor(sh3);
            a3 = sh3 - fl3;
            o3 = static_cast<int>(fl3);
            ka = kc + o3;
            kb = ka + 1;
          }
          Real dv = 0.0, fl[2] = {0.0, 0.0};
          if (plm_) {
            dv = VetPlmFoot([&](const int dk, const int dj) {
              return ip_(m,pr,r,ka + dk,ja + dj);
            }, a2, a3, thrd, fl);
          }
          if (f2s_) {
            // the upwind plane's intensity carried across the low x2 (x3) face by the
            // shift of the foot interpolation (VetShortCharSweep)
            auto g2 = [&](const int d) {
              if (!thrd) return ip_(m,pr,r,kc,jc + d);
              return 0.5*(ip_(m,pr,r,kc,jc + d) + (1.0 - a3)*ip_(m,pr,r,ka,jc + d)
                          + a3*ip_(m,pr,r,kb,jc + d));
            };
            auto g3 = [&](const int d) {
              return 0.5*(ip_(m,pr,r,kc + d,jc) + (1.0 - a2)*ip_(m,pr,r,kc + d,ja)
                          + a2*ip_(m,pr,r,kc + d,jb));
            };
            auto phi = [&](const Real sh, const auto &g) {
              const Real sa = fabs(sh), fa = sa - floor(sa);
              const int na = static_cast<int>(floor(sa)), sg = (sh < 0.0) ? 1 : -1;
              Real p = fa*g((sg > 0) ? -1 - na : na);
              for (int q = 0; q < na; ++q) p += g((sg > 0) ? -1 - q : q);
              return sg*p;
            };
            add4(up ? fnu : fnd, wr*am1*dx2/dx1*(phi(sh2, g2) + fl[0]));
            if (thrd) {
              add4((up ? fnu : fnd) + 12, wr*am1*dx3/dx1*(phi(sh3, g3) + fl[1]));
            }
          }
          const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
          const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
          const Real iup_v = waa*ip_(m,pr,r,ka,ja) + wab*ip_(m,pr,r,ka,jb)
                             + wba*ip_(m,pr,r,kb,ja) + wbb*ip_(m,pr,r,kb,jb) + dv;
          const Real cu = waa*cs_(m,li,0,ka,ja) + wab*cs_(m,li,0,ka,jb)
                          + wba*cs_(m,li,0,kb,ja) + wbb*cs_(m,li,0,kb,jb);
          const Real su = waa*cs_(m,li,1,ka,ja) + wab*cs_(m,li,1,ka,jb)
                          + wba*cs_(m,li,1,kb,ja) + wbb*cs_(m,li,1,kb,jb);
          // vet_seg = jump: the corners' S one x1 layer further upwind (VetSegFace; 0
          // beyond a physical x1 face: q = 1, the boundary cell keeps jump)
          const bool juu = spl && spe == 2;
          const int lu2 = up ? (li - 1) : (li + 1);
          auto s2 = [&](const int kq, const int jq) {
            return juu ? cs_(m,lu2,1,kq,jq) : 0.0;
          };
          const Real squ[4] = {s2(ka,ja), s2(ka,jb), s2(kb,ja), s2(kb,jb)};
          const Real *suq = juu ? squ : nullptr;
          // and the own cell's S one x1 layer downstream (VetSegSplit)
          const Real sdd = juu ? cs_(m,up ? (li + 2) : (li - 2),1,kc,jc) : -1.0;
          if (vim != 0) {   // vet_interp: log | corner
            const Real wq[4] = {waa, wab, wba, wbb};
            const Real iq[4] = {ip_(m,pr,r,ka,ja) + dv, ip_(m,pr,r,ka,jb) + dv,
                                ip_(m,pr,r,kb,ja) + dv, ip_(m,pr,r,kb,jb) + dv};
            const Real cq[4] = {cs_(m,li,0,ka,ja), cs_(m,li,0,ka,jb), cs_(m,li,0,kb,ja),
                                cs_(m,li,0,kb,jb)};
            const Real sq[4] = {cs_(m,li,1,ka,ja), cs_(m,li,1,ka,jb), cs_(m,li,1,kb,ja),
                                cs_(m,li,1,kb,jb)};
            iv = VetUpCorner(vim, wq, iq, cq, sq, s0, c0, dx1/am1, spl, spe,
                             f2s_ ? &ifc : nullptr, suq, sdd);
            int r2 = 0, r3 = 0;
            const Real f = xcr ? fmin(VetCrossFrac(sh2, r2),
                                      thrd ? VetCrossFrac(sh3, r3) : 1.0) : 1.0;
            if (f < 1.0 && VetCrossEdge(wq, cq, cs_(m,li,0,kc + r3,jc + r2))) {
              // vet_seg_cross: the stencil over 1 - f, then the cell crossed
              Real ib;
              (void) VetUpCorner(vim, wq, iq, cq, sq, s0, c0, (1.0 - f)*dx1/am1, spl,
                                 spe, &ib, suq, sdd);
              iv = VetSegCross(ib, cs_(m,li,0,kc + r3,jc + r2),
                               cs_(m,li,1,kc + r3,jc + r2), s0, c0, dx1/am1, f, spe,
                               f2s_ ? &ifc : nullptr,
                               juu ? s2(kc + r3,jc + r2) : -1.0, sdd);
            }
          } else if (spl) {
            // vet_seg = split, the two halves at the x1 face
            iv = VetSegSplit(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1, spe,
                             f2s_ ? &ifc : nullptr,
                             juu ? (waa*squ[0] + wab*squ[1] + wba*squ[2] + wbb*squ[3])
                                 : -1.0, sdd);
          } else {
            const Real dtau = 0.5*(cu + c0)*dx1/am1;
            Real ex, w0, wu;
            VetSegW(dtau, ex, w0, wu);
            iv = iup_v*ex + wu*su + w0*s0;
            // the face crossing at the face (0.5 chi_u hl from u), not at half the tau
            if (f2s_) ifc = VetSegFaceI(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1,
                                        0);
          }
          if (f2s_) add4(up ? f2u + 8 : fnd + 8, wr*m1*fmax(ifc, 0.0));
        }
        iv = fmax(iv, 0.0);
        ip_(m,pw,r,kc,jc) = iv;
        if (f2s_ && !x1per && l == nx1g - 1) {
          // the ray leaving through the x1 boundary face: the half cell at chi_0, S_0
          const Real ex = exp(-0.5*c0*dx1/am1);
          add4((up ? f2u : f2d) + 4, wr*m1*(iv*ex + (1.0 - ex)*s0));
        }
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
        if (hrs_) {   // the half-range sums, as VetShortChar
          Real *hp = up ? hpu : hpd;
          hp[0] += a*fmax(m1, 0.0);
          hp[1] += a*fmax(m2, 0.0);
          hp[2] += a*fmax(m3, 0.0);
        }
        if (kps_) {
          Real *kp = up ? kpu : kpd;
          const Real mu[3] = {m1, m2, m3};
          for (int d = 0; d < 3; ++d) {
            if (!(mu[d] > 0.0)) continue;
            for (int b = 0; b < 3; ++b) kp[3*d+b] += a*mu[b]*mu[d];
          }
        }
        if (nls_) {
          if (up) {nlu += wr*fabs(iv - s0);} else {nld += wr*fabs(iv - s0);}
        }
      }
      if (!accum) return;
      for (int n = 0; n < 10; ++n) {
        if (ua) {
          vc_(m,M1_VET_J+n,k,j,iu) += acu[n];
        }
        if (da) {
          vc_(m,M1_VET_J+n,k,j,id) += acd[n];
        }
      }
      if (hrs_) {
        for (int n = 0; n < 3; ++n) {
          if (ua) {vc_(m,M1_VET_HP1+n,k,j,iu) += hpu[n];}
          if (da) {vc_(m,M1_VET_HP1+n,k,j,id) += hpd[n];}
        }
      }
      if (nls_) {
        if (ua) {vc_(m,M1_VET_NL,k,j,iu) += nlu;}
        if (da) {vc_(m,M1_VET_NL,k,j,id) += nld;}
      }
      if (kps_) {
        for (int n = 0; n < M1_VET_NKP; ++n) {
          if (ua) {vc_(m,M1_VET_KP+n,k,j,iu) += kpu[n];}
          if (da) {vc_(m,M1_VET_KP+n,k,j,id) += kpd[n];}
        }
      }
      if (f2s_) {
        for (int n = 0; n < M1_VET_NFS; ++n) {
          if (ua) {
            vc_(m,M1_VET_F2+n,k,j,iu) += f2u[n];
            if (lu > 0) {
              vc_(m,M1_VET_F2+n,k,j,iu-1) += fnu[n];
            } else {
              sp_(m,0,n,kc,jc) += fnu[n];
            }
          }
          if (da) {
            vc_(m,M1_VET_F2+n,k,j,id) += f2d[n];
            if (ld < nx1 - 1) {
              vc_(m,M1_VET_F2+n,k,j,id+1) += fnd[n];
            } else {
              sp_(m,1,n,kc,jc) += fnd[n];
            }
          }
        }
      }
    });
    if (gl < ngl - 1) {
      for (int p = 1; p < st.npb; ++p) {
        VetPlanePost(st, st.ipl, pw, nrl, nsb, nx2, nx3, p);
        VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
      }
      VetPlanePost(st, st.ipl, pw, nrl, nsb, nx2, nx3, st.npb);
      VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
    }
  }
  if (f2s_) {
    // the face sums of the x1 neighbours' boundary cells
    VetX1Move(st, st.spl, 1, -2, 2, 0, 0, M1_VET_NFS, nsb);
    VetX1Move(st, st.spl, 0, -2, 3, 1, 0, M1_VET_NFS, nsb);
    par_for("vet_vet_mb_spill", DevExeSpace(), 0, nmb-1, 0, M1_VET_NFS-1, ks, ke, js,
            je,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j) {
      vc_(m,M1_VET_F2+n,k,j,ie) += sp_(m,2,n,k-ks+w3,j-js+w2);
      vc_(m,M1_VET_F2+n,k,j,is) += sp_(m,3,n,k-ks+w3,j-js+w2);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetShortChar

void RadiationVet::VetShortChar() {
  // (vet_rotate_avg; as VetGdBuild): the first build of a new block moves
  // the ring on by the Delta of the old block's last build (else re-init); the first
  // build of a run sweeps the K-1 previous sets; the last build of a block also sweeps
  // the NEXT set at the same state (Delta), the set then stays on the next one and
  // vet_cell of this build is put back
  const bool ra = vra_on && vsc_rot;
  const int ilo = pmy_pack->pmesh->mb_indcs.is;
  const int64_t b = VetGdBlock(pmy_pack->pmesh->ncycle);
  if (ra && vra_state == 1 && b != vra_blk) {
    if (vra_pend == 3) {
      VetScRotAvgSetup();
      RotAvgShift(vra_k - 1);
      vra_pend = 0;
      vra_nrot += 1.0;
    } else {
      vra_state = 0;
    }
  }
  if (ra && vra_state == 0) {
    Kokkos::Timer tr;
    vra_head = 0;
    for (int jb = 1; jb < vra_k; ++jb) {
      VetRotSet(VetRotBeta(b - jb));
      vra_raw = true;
      VetShortCharSweep();
      vra_raw = false;
      VetScRotAvgSetup();
      RotAvgSnap(jb - 1, ilo);
      vra_nsweep += 1.0;
    }
    vra_pend = 1;
    vra_state = 1;
    vra_ninit += 1.0;
    vra_time += tr.seconds();
  }
  VetRotUpdate();   // the set of this step's rotation block
  vra_last = ra && vra_state == 1 && RotLastStep();
  VetShortCharSweep();
  if (vra_last) {
    Kokkos::Timer tr;
    if (vra_bk_vc.extent(0) != vet_cell.extent(0) ||
        vra_bk_vc.extent(1) != vet_cell.extent(1)) {
      Kokkos::realloc(vra_bk_vc, vet_cell.extent(0), vet_cell.extent(1),
                      vet_cell.extent(2), vet_cell.extent(3), vet_cell.extent(4));
    }
    Kokkos::deep_copy(vra_bk_vc, vet_cell);
    VetRotSet(VetRotBeta(b + 1));
    vra_raw = true;
    VetShortCharSweep();
    vra_raw = false;
    VetScRotAvgSetup();
    RotAvgDelta(vra_k - 1, ilo);
    Kokkos::deep_copy(vet_cell, vra_bk_vc);
    vra_pend = 3;
    vra_nsweep += 1.0;
    vra_time += tr.seconds();
  }
  vra_last = false;
  vra_blk = b;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetShortCharSweep
//! \brief the vet_sc build with the current set (split from VetShortChar)

void RadiationVet::VetShortCharSweep() {
  // vet_source_noesrc for vet_sc (as vet_col, fixbundle-1009 F1): the
  // source and the boundary intensity read the PHYSICAL E^n = u0 under be, not the
  // solve's old vector E^n + dt (chat/c) esrc (M1_IW_EN)
  const bool noes_ = VsrcNoes();
  auto esv_ = noes_ ? u0 : iw;
  const int esc_ = noes_ ? M1_E : M1_IW_EN;
  // the half-range sums sum w I max(mu_a, 0) (M1_VET_HP1..3);
  // every sweep variant carries them (the banded MB sweep in its moment arrays)
  const bool hrs_ = impl_beam_hr;
  // scheme = clean: the half-range stress columns K+^d_bd (M1_VET_KP, boundary fluxes)
  const bool kps_ = true;
  // (implicit_thin_switch = fs | fskn): the non-locality sum w |I - S| of the
  // rays at each cell (M1_VET_NL, C1 of M1ThinCell)
  const bool nls_ = thsw_c1;
  // clean_consistency: the rays' face fluxes (M1_VET_F2, F1, FX): the x1 flux of the
  // intensity where the segment crosses the face, the transverse part of the linear foot
  // interpolation as an x2 flux of the upwind plane's layer, so that the cell balance of
  // the face fluxes is the exchange of the two half-segments inside the cell
  const bool f2s_ = cln_consist;
  const bool x1per = vet_x1per;
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
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real ar = arad;
  const Real efl = e_floor;
  const bool thermal = fl_on && coupling && !opac_zero;
  const int vsf = vscat ? vscat_form : 0;   // vet_scatter
  const Real bblo = VetBathLo();   // vet_bc_inner
  const Real bbhi = (vet_bc_bath && ibc_x1max == M1_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool bfc = vet_bc_face;
  const bool x1src = vet_x1per && vet_x1src;
  const bool spl = vet_seg_split;
  const bool xcr = vet_xcr;
  const int spe = VetSegMode();
  const bool x2o = vet_x2open, x3o = vet_x3open;   // vacuum x2 (x3) faces
  const int vim = vet_interp;   // (#6)
  const bool plm_ = true;
  // the boundary beam rays (ImplBeamInit)
  const int bmn_ = ibm_nr;
  const int bmr0 = ibm_ray[0], bmr1 = ibm_ray[1], bmr2 = ibm_ray[2], bmr3 = ibm_ray[3];
  auto bmi_ = ibm_i;

  // (1) extinction and source, cell by cell
  auto opd_ = opac;
  auto iwd_ = iw;
  auto vcd_ = vet_cell;
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  // the VET source rule of vet_col / vet_gd (VetLatSweep): vet_source = relaxed (srx:
  // the state relaxed by a local backward-Euler exchange over the step,
  // VetRelaxedSource) or gas (the state's own source); under time2_vet = predict |
  // rebuild time2_vet_source = state | relaxed | blend (Time2VetColAt sets vcol_srelax
  // / vcol_sbl for the build).  vet_sc before: always the state's source.
  const bool srx = thermal && vcol_srelax;
  const bool sbl = vcol_sbl;
  const Real sx0 = vcol_sx0;
  // (vet_source = relaxed_esc): the x1 cell width for the escape
  const bool sesc = srx && vcol_sesc;
  const Real clq = gas_feedback ? c_light : 0.0;   // held gas: T* = T_b
  auto &mbs_ = pmy_pack->pmb->mb_size;
  const int vsr = (vsf != 0 && vscat_jrel) ? (vsf + 2) : vsf;   // relaxed J
  FluidRef flv = FluidRef::Get(pmy_pack);
  auto eos = flv.eos;
  const bool cvx = cln_srccvx && !srx;   // clean_source_convex
  auto uh = flv.u0;
  const Real chs = chat;
  const Real dts = mr_on ? mr_dt : pmy_pack->pmesh->dt;
  VetFor("vet_vet_srcc", (nmb1 + 1)*nk*nj*ni, KOKKOS_LAMBDA(const int t) {
    const int i = is + t%ni;
    int q = t/ni;
    const int j = js + q%nj;
    q /= nj;
    const int k = ks + q%nk, m = q/nk;
    Real chx = fmax(opd_(m,M1_OP_T,k,j,i), 1.0e-300);
    Real e = fmax(esv_(m,esc_,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iwd_(m,M1_IW_TP,k,j,i);
      Real t2 = tg*tg;
      if (srx) {
        const Real eth = fmin(opd_(m,M1_OP_P,k,j,i)/chx, 1.0);
        const Real rks = (vsf != 0) ? opd_(m,M1_OP_S,k,j,i) : 0.0;
        // relaxed_esc: the escape rate of the cell (0 = relaxed)
        const Real qe = VetEscapeRate(sesc, chx, mbs_.d_view(m).dx1);
        s = sesc ? (sbl ? VetRelaxedSourceBQ(eos, uh(m,IDN,k,j,i), tg, e,
                                    opd_(m,M1_OP_P,k,j,i),
                                    opd_(m,M1_OP_E,k,j,i), eth, ar, clq, chs, dts, vsr,
                                    chx, rks, sx0, qe)
                : VetRelaxedSourceQ(eos, uh(m,IDN,k,j,i), tg, e, opd_(m,M1_OP_P,k,j,i),
                                   opd_(m,M1_OP_E,k,j,i), eth, ar, clq, chs, dts, vsr,
                                   chx, rks, qe))
             : (sbl ? VetRelaxedSourceB(eos, uh(m,IDN,k,j,i), tg, e,
                                    opd_(m,M1_OP_P,k,j,i),
                                    opd_(m,M1_OP_E,k,j,i), eth, ar, cl, chs, dts, vsr,
                                    chx, rks, sx0)
                : VetRelaxedSource(eos, uh(m,IDN,k,j,i), tg, e, opd_(m,M1_OP_P,k,j,i),
                                   opd_(m,M1_OP_E,k,j,i), eth, ar, cl, chs, dts, vsr,
                                   chx, rks));
      } else if (vsf != 0) {   // vet_scatter (VetScatterMix, rad_vet.hpp)
        const Real rkp = cvx ? fmin(opd_(m,M1_OP_P,k,j,i), chx) : opd_(m,M1_OP_P,k,j,i);
        s = VetScatterMix(vsf, rkp, chx, opd_(m,M1_OP_S,k,j,i), ar*t2*t2, e);
      } else {
        Real eth = fmin(opd_(m,M1_OP_P,k,j,i)/chx, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
    }
    vcd_(m,M1_VET_CHX,k,j,i) = chx;
    vcd_(m,M1_VET_SRC,k,j,i) = s;
    for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {
      vcd_(m,n,k,j,i) = 0.0;
    }
    if (hrs_) {
      for (int n = 0; n < 3; ++n) {vcd_(m,M1_VET_HP1+n,k,j,i) = 0.0;}
    }
    if (kps_) {
      for (int n = 0; n < M1_VET_NKP; ++n) {vcd_(m,M1_VET_KP+n,k,j,i) = 0.0;}
    }
    if (f2s_) {
      for (int n = 0; n < M1_VET_NFS; ++n) {vcd_(m,M1_VET_F2+n,k,j,i) = 0.0;}
    }
    if (nls_) {
      vcd_(m,M1_VET_NL,k,j,i) = 0.0;
      // vacuum (no extinction): no lagged scattering source, eps = 1
      const Real eps = thermal ? fmin(opd_(m,M1_OP_P,k,j,i)/chx, 1.0) : 0.0;
      vcd_(m,M1_VET_EPS,k,j,i) = !(opd_(m,M1_OP_T,k,j,i) > 0.0) ? 1.0 : eps;
    }
  });

  // (2) the sweep: launch l does layer is+l for the upward rays and ie-l for the
  // downward ones, reading the intensity of launch l-1 from the other plane buffer.
  // Several MeshBlocks: the banded sweep over the global layers (VetSweepMB).
  if (vet_mbs != nullptr) {
    VetSweepMB();
  }
  const int nlaunch = (vet_mbs != nullptr) ? 0 : nx1;   // the banded sweep ran instead
  // vet_bc_reflect = mirror: at a REFLECTING x1 end the entering
  // rays carry the leaving intensity of their x1 mirror ray (ray r +- nray/2: both angle
  // sets list the upward half first, the mirror of up ray r at r + nray/2).  A second
  // pass reads it from the plane buffer the first pass left at the boundary layer, as the
  // periodic image is read (the mirror image of the boundary cell is the cell itself):
  // the interior SC segment of length dx1/|mu1| through the boundary cell and its image.
  // One reflecting end: exact in 2 passes (the leaving rays do not depend on it); both:
  // a closed cavity, vet_x1_npass passes (default 4; each pass adds one round trip).
  const bool mlo = vet_bc_mirror && (ibc_x1min == M1_IBC_REFLECT) && !vet_x1per;
  const bool mhi = vet_bc_mirror_hi && (ibc_x1max == M1_IBC_REFLECT) && !vet_x1per;
  if ((mlo || mhi) && vet_mbs != nullptr) {
    VetFatal("<rad_vet>/vet_bc_reflect = mirror (the default with a reflecting x1 end) "
             "needs a single MeshBlock on vet_sc (the banded sweep has no mirror); name "
             "vet_bc_reflect = open");
  }
  const int nhray = nray/2;
  const int npass = vet_x1per ? vet_x1npass : ((mlo && mhi) ? vet_mirror_npass
                                               : ((mlo || mhi) ? 2 : 1));
  // vet_x1_ptol > 0: adaptive periodic passes.  Every pass accumulates into
  // moments zeroed at its start, so the moments kept are those of the last pass run (the
  // same values as accumulating the last pass only); after each pass the exit plane is
  // compared with the previous one (one device reduction, one host read per pass)
  const bool adapt = vet_x1per && (vet_x1ptol > 0.0) && (nlaunch > 0);
  Real dprev = -1.0;
  int npdone = 0;
  bool conv = false;
  // an adaptive sweep after the first build starts from the last build's exit plane (its
  // periodic image then), not from the boundary inflow
  const bool warm = adapt && (vet_x1_nbuild > 0.0);
  if (warm) {
    auto ip_ = vet_ipl;
    auto ent_ = vet_x1ent;
    par_for("vet_vet_x1warm", DevExeSpace(), 0, nmb1, 0, nray - 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
      ip_(m,1,r,k,j) = ent_(m,r,k,j);
    });
  }
  for (int p = 0; p < npass; ++p) {
  const bool wrap = vet_x1per && (p > 0 || warm);   // layer is (ie) reads the image
  const bool mrlo = mlo && (p > 0), mrhi = mhi && (p > 0);   // ... or the mirror image
  // only the last pass adds to the moments (adaptive: each pass, on zeroed moments)
  const bool accum = adapt || (p == npass - 1);
  if (adapt && p > 0) {
    par_for("vet_vet_x1zero", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {vc_(m,n,k,j,i) = 0.0;}
      if (hrs_) {for (int n = 0; n < 3; ++n) {vc_(m,M1_VET_HP1+n,k,j,i) = 0.0;}}
      if (kps_) {for (int n = 0; n < M1_VET_NKP; ++n) {vc_(m,M1_VET_KP+n,k,j,i) = 0.0;}}
      if (f2s_) {for (int n = 0; n < M1_VET_NFS; ++n) {vc_(m,M1_VET_F2+n,k,j,i) = 0.0;}}
      if (nls_) {vc_(m,M1_VET_NL,k,j,i) = 0.0;}
    });
  }
  for (int l = 0; l < nlaunch; ++l) {
    const int gl = p*nlaunch + l;
    const int pw = gl & 1, pr = pw ^ 1;
    const int iu = is + l, id = ie - l;
    // ray r of launch l in the cell (k, j): its intensity (into the plane buffer) and
    // the scalars g of its face sums (VetScAcc::Add), shared by the per-cell kernel and
    // the ray-parallel one (vet_sc_kernel = ray: one thread per (block, ray, cell), then
    // the sums per cell in the same ray order, the same bits)
    auto rayf = KOKKOS_LAMBDA(const int m, const int k, const int j, const int r,
                              Real g[5]) {
      const Real dx1 = mbsize.d_view(m).dx1;
      const Real dx2 = mbsize.d_view(m).dx2;
      const Real dx3 = mbsize.d_view(m).dx3;
      const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
      const bool up = (m1 > 0.0);
      const int i = up ? iu : id;
      const Real am1 = fabs(m1);
      const Real c0 = vc_(m,M1_VET_CHX,k,j,i);
      const Real s0 = vc_(m,M1_VET_SRC,k,j,i);
      Real iv;
      const bool mimg = (l == 0) && (up ? mrlo : mrhi);   // vet_bc_reflect = mirror
      for (int q = 0; q < 5; ++q) {g[q] = 0.0;}
      Real ifc = 0.0;   // clean_consistency: the intensity at the face crossing
      if (l == 0 && !wrap && !mimg) {
        // half a cell from the boundary face; S at the face linearly extrapolated
        // from the two boundary cells (exact for the deep linear S of a diffusion
        // regime), and linear in tau along the segment
        const int iin = up ? (i + 1) : (i - 1);
        const Real s1 = vc_(m,M1_VET_SRC,k,j,iin);
        const Real su = fmax(1.5*s0 - 0.5*s1, 0.0);
        Real ib = up ? 0.0 : bbhi;
        if (x1src) {
          ib = su;   // vet_x1_periodic_inflow = source: symmetric entry
        } else if (up) {
          if (bblo >= 0.0) {
            ib = bblo;
          } else {
            // diffusion: eps = E + 3 (F . Omega)/c from the bottom cell's M1 state;
            // vet_bc_diff = face: E moved to the FACE by the source's
            // linear extrapolation, E_f = E + (S_f - S_0) (= 1.5 E_0 - 0.5 E_1 when
            // S = E), the face value the half-cell carry with linear S expects
            Real e = esv_(m,esc_,k,j,i);
            if (bfc) e += 0.5*(s0 - s1);
            e = fmax(e, efl);
            Real fo = iw_(m,M1_IW_F1,k,j,i)*m1 + iw_(m,M1_IW_F2,k,j,i)*m2
                      + iw_(m,M1_IW_F3,k,j,i)*m3;
            ib = fmax(e + 3.0*fo/cl, 0.0);
          }
          if (bmn_ > 0) {   // boundary beams
            const int q = (r == bmr0) ? 0 : ((r == bmr1) ? 1 : ((r == bmr2) ? 2
                          : ((r == bmr3) ? 3 : -1)));
            if (q >= 0 && q < bmn_) {ib += bmi_(m,q,k-ks,j-js);}
          }
        }
        if (f2s_) g[0] = wr*m1*ib;
        const Real dtau = 0.5*c0*dx1/am1;
        Real ex, w0, wu;
        VetSegW(dtau, ex, w0, wu);
        iv = ib*ex + wu*su + w0*s0;
      } else {
        const int iup = mimg ? i : ((l == 0) ? (up ? ie : is) : (up ? (i - 1) : (i + 1)));
        const int rb = mimg ? ((r < nhray) ? (r + nhray) : (r - nhray)) : r;
        // the foot of the characteristic on the upwind plane, in cell units
        const Real sh2 = -dx1*m2/(am1*dx2);
        const Real fl2 = floor(sh2);
        const Real a2 = sh2 - fl2;
        const int o2 = static_cast<int>(fl2);
        int ja = js + VetWrap(j - js + o2, nx2);
        int jb = js + VetWrap(j - js + o2 + 1, nx2);
        // vacuum x2 (x3) faces: a foot point outside the box carries
        // I = 0 (the ray entered through the face) and the edge cell's chi and S
        Real ga = 1.0, gb = 1.0;
        if (x2o) {
          const int qa = j - js + o2, qb = qa + 1;
          ga = (qa >= 0 && qa < nx2) ? 1.0 : 0.0;
          gb = (qb >= 0 && qb < nx2) ? 1.0 : 0.0;
          ja = js + ((qa < 0) ? 0 : ((qa >= nx2) ? (nx2 - 1) : qa));
          jb = js + ((qb < 0) ? 0 : ((qb >= nx2) ? (nx2 - 1) : qb));
        }
        int ka = k, kb = k, o3 = 0;
        Real a3 = 0.0, ka_g = 1.0, kb_g = 1.0;
        const Real sh3 = thrd ? -dx1*m3/(am1*dx3) : 0.0;
        if (thrd) {
          const Real fl3 = floor(sh3);
          a3 = sh3 - fl3;
          o3 = static_cast<int>(fl3);
          ka = ks + VetWrap(k - ks + o3, nx3);
          kb = ks + VetWrap(k - ks + o3 + 1, nx3);
          if (x3o) {
            const int qa = k - ks + o3, qb = qa + 1;
            ka_g = (qa >= 0 && qa < nx3) ? 1.0 : 0.0;
            kb_g = (qb >= 0 && qb < nx3) ? 1.0 : 0.0;
            ka = ks + ((qa < 0) ? 0 : ((qa >= nx3) ? (nx3 - 1) : qa));
            kb = ks + ((qb < 0) ? 0 : ((qb >= nx3) ? (nx3 - 1) : qb));
          }
        }
        // the upwind plane at offsets (dk, dj) from the cell, 0 beyond a vacuum face
        auto ipq = [&](const int dk, const int dj) {
          const int qj = j - js + dj, qk = k - ks + dk;
          if ((x2o && (qj < 0 || qj >= nx2)) || (x3o && (qk < 0 || qk >= nx3))) {
            return 0.0;
          }
          return ip_(m,pr,rb,thrd ? ks + VetWrap(qk, nx3) : k,js + VetWrap(qj, nx2));
        };
        Real dv = 0.0, fl[2] = {0.0, 0.0};
        if (plm_) {
          dv = VetPlmFoot([&](const int dk, const int dj) {
            return ipq(o3 + dk, o2 + dj);
          }, a2, a3, thrd, fl);
        }
        if (f2s_) {
          // Phi: the upwind plane's intensity carried across the low x2 (x3) face by
          // the shift -sh2 (-sh3) cells of the foot interpolation.  The x2 and x3 parts
          // S2, S3 of the interpolation commute: the x2 face carries (I + S3 I)/2, the
          // x3 face (I + S2 I)/2, and the cell balance of the two is S3 S2 I - I.
          auto g2 = [&](const int d) {
            return thrd ? 0.5*(ipq(0,d) + (1.0 - a3)*ipq(o3,d) + a3*ipq(o3+1,d))
                        : ipq(0,d);
          };
          auto g3 = [&](const int d) {
            return 0.5*(ipq(d,0) + (1.0 - a2)*ipq(d,o2) + a2*ipq(d,o2+1));
          };
          auto phi = [&](const Real sh, const auto &g) {
            const Real s = fabs(sh), f = s - floor(s);
            const int n = static_cast<int>(floor(s)), sg = (sh < 0.0) ? 1 : -1;
            Real p = f*g(sg > 0 ? -1 - n : n);
            for (int q = 0; q < n; ++q) p += g(sg > 0 ? -1 - q : q);
            return sg*p;
          };
          // (vet_bc_reflect = mirror: the upwind plane is the wall column's mirror
          // image, so its shift crosses the wall column's own x2 faces: the reflected
          // ray changes rows there, and the wall face carries no net flux)
          g[1] = wr*am1*dx2/dx1*(phi(sh2, g2) + fl[0]);
          if (thrd) {
            g[2] = wr*am1*dx3/dx1*(phi(sh3, g3) + fl[1]);
          }
        }
        const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
        const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
        // the vacuum-face gates on the mirror-aware ray rb
        // (#2); a gate is exactly 1 off the open faces
        const Real gaa = ga*ka_g, gab = gb*ka_g, gba = ga*kb_g, gbb = gb*kb_g;
        Real iup_v;
        if (x2o || x3o) {
          iup_v = waa*gaa*ip_(m,pr,rb,ka,ja) + wab*gab*ip_(m,pr,rb,ka,jb)
                  + wba*gba*ip_(m,pr,rb,kb,ja) + wbb*gbb*ip_(m,pr,rb,kb,jb);
        } else {
          iup_v = waa*ip_(m,pr,rb,ka,ja) + wab*ip_(m,pr,rb,ka,jb)
                  + wba*ip_(m,pr,rb,kb,ja) + wbb*ip_(m,pr,rb,kb,jb);
        }
        iup_v += dv;
        const Real cu = waa*vc_(m,M1_VET_CHX,ka,ja,iup)
                        + wab*vc_(m,M1_VET_CHX,ka,jb,iup)
                        + wba*vc_(m,M1_VET_CHX,kb,ja,iup)
                        + wbb*vc_(m,M1_VET_CHX,kb,jb,iup);
        const Real su = waa*vc_(m,M1_VET_SRC,ka,ja,iup)
                        + wab*vc_(m,M1_VET_SRC,ka,jb,iup)
                        + wba*vc_(m,M1_VET_SRC,kb,ja,iup)
                        + wbb*vc_(m,M1_VET_SRC,kb,jb,iup);
        // vet_seg = jump: the corners' S one x1 cell further upwind (VetSegFace), its
        // periodic or mirror image beyond an x1 end; none (q = 1: the boundary cell keeps
        // jump) beyond an open one (vet_cell holds no ghost CHX, SRC)
        int iuu = up ? (i - 2) : (i + 2);
        if (x1per) {
          iuu = is + VetWrap(iuu - is, nx1);
        } else if (up ? (iuu < is && mrlo) : (iuu > ie && mrhi)) {
          iuu = up ? (2*is - 1 - iuu) : (2*ie + 1 - iuu);
        }
        const bool juu = spl && spe == 2 && iuu >= is && iuu <= ie;
        auto s2 = [&](const int kq, const int jq) {
          return juu ? vc_(m,M1_VET_SRC,kq,jq,iuu) : 0.0;
        };
        const Real squ[4] = {s2(ka,ja), s2(ka,jb), s2(kb,ja), s2(kb,jb)};
        const Real *suq = juu ? squ : nullptr;
        // and the own cell's S one x1 cell downstream (VetSegSplit; active cells only)
        int idn = up ? (i + 1) : (i - 1);
        if (x1per) idn = is + VetWrap(idn - is, nx1);
        const Real sdd = (spl && spe == 2 && idn >= is && idn <= ie)
                         ? vc_(m,M1_VET_SRC,k,j,idn) : -1.0;
        if (vim != 0) {   // vet_interp: log | corner
          const Real wq[4] = {waa, wab, wba, wbb};
          // a corner outside a vacuum face carries I = 0 (gate exactly 1
          // elsewhere, bitwise), its segment the edge cell's chi and S, as lin above
          const Real iq[4] = {gaa*ip_(m,pr,rb,ka,ja) + dv, gab*ip_(m,pr,rb,ka,jb) + dv,
                              gba*ip_(m,pr,rb,kb,ja) + dv, gbb*ip_(m,pr,rb,kb,jb) + dv};
          const Real cq[4] = {vc_(m,M1_VET_CHX,ka,ja,iup), vc_(m,M1_VET_CHX,ka,jb,iup),
                              vc_(m,M1_VET_CHX,kb,ja,iup), vc_(m,M1_VET_CHX,kb,jb,iup)};
          const Real sq[4] = {vc_(m,M1_VET_SRC,ka,ja,iup), vc_(m,M1_VET_SRC,ka,jb,iup),
                              vc_(m,M1_VET_SRC,kb,ja,iup), vc_(m,M1_VET_SRC,kb,jb,iup)};
          iv = VetUpCorner(vim, wq, iq, cq, sq, s0, c0, dx1/am1, spl, spe,
                           f2s_ ? &ifc : nullptr, suq, sdd);
          int r2 = 0, r3 = 0;
          Real f = xcr ? fmin(VetCrossFrac(sh2, r2),
                              thrd ? VetCrossFrac(sh3, r3) : 1.0) : 1.0;
          const int q2 = j - js + r2, q3 = k - ks + r3;
          if ((x2o && (q2 < 0 || q2 >= nx2)) || (x3o && (q3 < 0 || q3 >= nx3))) f = 1.0;
          const int jq = js + VetWrap(q2, nx2), kq = thrd ? ks + VetWrap(q3, nx3) : k;
          if (f < 1.0 && VetCrossEdge(wq, cq, vc_(m,M1_VET_CHX,kq,jq,iup))) {
            // vet_seg_cross: the stencil over 1 - f, then the cell crossed
            Real ib;
            (void) VetUpCorner(vim, wq, iq, cq, sq, s0, c0, (1.0 - f)*dx1/am1, spl, spe,
                               &ib, suq, sdd);
            iv = VetSegCross(ib, vc_(m,M1_VET_CHX,kq,jq,iup),
                             vc_(m,M1_VET_SRC,kq,jq,iup), s0, c0, dx1/am1, f, spe,
                             f2s_ ? &ifc : nullptr, juu ? s2(kq,jq) : -1.0, sdd);
          }
        } else if (spl) {
          // vet_seg = split, the two halves at the x1 face
          iv = VetSegSplit(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1, spe,
                           f2s_ ? &ifc : nullptr,
                           juu ? (waa*squ[0] + wab*squ[1] + wba*squ[2] + wbb*squ[3])
                               : -1.0, sdd);
        } else {
          const Real dtau = 0.5*(cu + c0)*dx1/am1;
          Real ex, w0, wu;
          VetSegW(dtau, ex, w0, wu);
          iv = iup_v*ex + wu*su + w0*s0;
          // the face crossing at the face (0.5 chi_u hl from u), not at half the tau
          if (f2s_) ifc = VetSegFaceI(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1,
                                      0);
        }
        // the face crossing of the single segment: a face flux
        if (f2s_ && !mimg) g[3] = wr*m1*fmax(ifc, 0.0);
      }
      iv = fmax(iv, 0.0);
      ip_(m,pw,r,k,j) = iv;
      if (f2s_ && !x1per && i == (up ? ie : is)) {
        // the ray leaving through the x1 boundary face: the half cell at chi_0, S_0
        const Real ex = exp(-0.5*c0*dx1/am1);
        g[4] = wr*m1*(iv*ex + (1.0 - ex)*s0);
      }
      return iv;
    };
    if (vsc_rayk) {
      auto rg_ = vet_rayg;
      par_for("vet_vet_sweep_ray", DevExeSpace(), 0, nmb1, 0, nray - 1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
        Real g[5];
        rayf(m, k, j, r, g);
        for (int q = 0; q < 5; ++q) {rg_(m,q*nray + r,k,j) = g[q];}
      });
      par_for("vet_vet_sweep_sum", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        VetScAcc a;
        for (int r = 0; r < nray; ++r) {
          const Real g[5] = {rg_(m,r,k,j), rg_(m,nray + r,k,j), rg_(m,2*nray + r,k,j),
                             rg_(m,3*nray + r,k,j), rg_(m,4*nray + r,k,j)};
          const int i = (ang_(r,0) > 0.0) ? iu : id;
          a.Add(ang_, r, ip_(m,pw,r,k,j), vc_(m,M1_VET_SRC,k,j,i), g, hrs_, kps_, nls_);
        }
        if (accum) {
          a.Write(vc_, m, k, j, iu, id, is, ie, mrlo, mrhi, hrs_, kps_, nls_, f2s_);
        }
      });
    } else {
      par_for("vet_vet_sweep", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        VetScAcc a;
        for (int r = 0; r < nray; ++r) {
          Real g[5];
          const Real iv = rayf(m, k, j, r, g);
          const int i = (ang_(r,0) > 0.0) ? iu : id;
          a.Add(ang_, r, iv, vc_(m,M1_VET_SRC,k,j,i), g, hrs_, kps_, nls_);
        }
        if (accum) {
          a.Write(vc_, m, k, j, iu, id, is, ie, mrlo, mrhi, hrs_, kps_, nls_, f2s_);
        }
      });
    }
  }
  ++npdone;
  if (adapt) {
    // the exit plane of this pass (launch nlaunch-1 wrote plane (gl & 1)): up rays at ie,
    // down rays at is, i.e. the entry of the next pass; d = max|I_p - I_{p-1}|, a = max I
    const int pl = (p*nlaunch + nlaunch - 1) & 1;
    auto ent_ = vet_x1ent;
    const int nj = je - js + 1, nk = ke - ks + 1;
    const int ntot = (nmb1 + 1)*nray*nk*nj;
    const bool cmp = (p > 0);
    Real dmax = 0.0, amax = 0.0;
    Kokkos::parallel_reduce("vet_vet_x1conv",
      Kokkos::RangePolicy<DevExeSpace>(0, ntot),
      KOKKOS_LAMBDA(const int t, Real &dl, Real &al) {
        int q = t;
        const int j = js + q % nj; q /= nj;
        const int k = ks + q % nk; q /= nk;
        const int r = q % nray;
        const int m = q / nray;
        const Real v = ip_(m,pl,r,k,j);
        if (cmp) {dl = fmax(dl, fabs(v - ent_(m,r,k,j)));}
        al = fmax(al, v);
        ent_(m,r,k,j) = v;
      }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(amax));
    if (p > 0) {
      const Real d = (amax > 0.0) ? (dmax/amax) : 0.0;
      // error of the entry this pass used: d + d q + d q^2 ... = d/(1 - q)
      const Real q = (dprev > 0.0) ? fmin(d/dprev, 0.999) : 0.999;
      conv = (d == 0.0) || (d/(1.0 - q) < vet_x1ptol);
      dprev = d;
      if (conv) {break;}
    }
  }
  }
  if (adapt) {
    vet_x1_nbuild += 1.0;
    vet_x1_npsum += static_cast<Real>(npdone);
    vet_x1_npmax = std::max(vet_x1_npmax, npdone);
    if (!conv) {
      if (vet_x1_ncap == 0.0 && global_variable::my_rank == 0) {
        std::cout << "<rad_vet> WARNING vet_sc periodic x1 sweep: vet_x1_npass="
                  << vet_x1npass << " reached before vet_x1_ptol=" << vet_x1ptol
                  << " (first time; counted in the final report)" << std::endl;
      }
      vet_x1_ncap += 1.0;
    }
  }

  // (vet_rotate_avg): the pending offsets from this (current-set) sweep, then
  // the mean over the last K sets of J and the J-normalised moments (not on a snapshot)
  if (vra_on && vsc_rot && !vra_raw && vra_state == 1) {VetScRotAvgSetup();}
  // a history loaded from the restart file must fit the table
  if (vra_on && vsc_rot && !vra_raw && vra_state == 1 && RotAvgRstCheck()) {
    Kokkos::Timer tr;
    if (vra_last) {RotAvgSnap(vra_k - 1, is);}   // the current set's raw moments (Delta)
    RotAvgUpdate(is);
    RotAvgApply(is);
    Kokkos::fence();
    vra_time += tr.seconds();
  }

  // (3) the uniaxial projection (chi, n) of D = K/J
  const bool axflux = vet_axis_flux;
  par_for("vet_vet_proj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
  vet_ncall += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetFullTensor
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
//!       (= 0) the spectrum is mapped affinely, lam -> emin + (lam' - emin)
//!       (1 - 3 emin)/(sum lam' - 3 emin) with lam' = max(lam, emin), which keeps the
//!       eigenvectors, the trace 1 and every eigenvalue in [emin, 1] (counted).
//! Consequences for the solve: D_dd >= emin >= 0 on every diagonal, so each x1 line row
//! gains nu df (D_11,i + D_11,i) >= 0 on its diagonal and loses nu df D_11,i+-1 <= 0 off
//! it -- the column-diagonally-dominant M-matrix of the uniaxial form, and the diagonal
//! stays >= 1.  The transverse diagonal TDIA >= 0 likewise.
//! Ghosts: periodic wrap in x2/x3 (the SC sweep already requires ONE MeshBlock with
//! periodic x2/x3), copy of the edge cell in x1 (never read at a physical x1 face).
//! Per-cell |D_guarded - D_raw| goes to M1_VET_GD for the run statistics.

void RadiationVet::VetFullTensor() {
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
  const bool x2o = vet_x2open, x3o = vet_x3open;
  auto vc_ = vet_cell;

  par_for("vet_vet_full", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
    par_for("vet_vet_full_dsc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
    par_for("vet_vet_full_ghmb", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
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
    par_for("vet_vet_full_gh", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (act) return;
      const int ii = (i < is) ? is : ((i > ie) ? ie : i);
      // edge copy at a vacuum (non-periodic) x2 (x3) face
      const int jj = x2o ? ((j < js) ? js : ((j > je) ? je : j))
                         : (js + VetWrap(j - js, nx2));
      const int kk = thrd ? (x3o ? ((k < ks) ? ks : ((k > ke) ? ke : k))
                                 : (ks + VetWrap(k - ks, nx3))) : k;
      for (int n = 0; n < 6; ++n) {
        vc_(m,M1_VET_D11+n,k,j,i) = vc_(m,M1_VET_D11+n,kk,jj,ii);
      }
    });
  }

  // run statistics of the guard (host reductions over the active cells)
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nc = (nmb1 + 1)*nk*nj*ni;
  Real nhit = 0.0, gmax = 0.0, dmin = 1.0;
  Kokkos::parallel_reduce("vet_vet_full_st", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
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
//! \fn void RadiationVet::VetReport

void RadiationVet::VetReport() {
  if (global_variable::my_rank != 0 || vet_ncall <= 0.0) return;
  std::cout << "<rad_vet> vet_sc: formal solutions=" << vet_ncall
            << " SC seconds=" << vet_time << " (" << vet_time/vet_ncall << " per call)"
            << " implicit-solve seconds incl. SC=" << vet_itime
            << " SC/(implicit solve without SC)="
            << vet_time/std::fmax(vet_itime - vet_time, 1.0e-300) << std::endl;
  if (vet_x1_nbuild > 0.0) {   // adaptive periodic x1 sweep
    std::cout << "<rad_vet> vet_sc periodic x1 sweep (vet_x1_ptol=" << vet_x1ptol
              << ", cap vet_x1_npass=" << vet_x1npass << "): builds=" << vet_x1_nbuild
              << " passes mean=" << vet_x1_npsum/vet_x1_nbuild << " max=" << vet_x1_npmax
              << " builds that hit the cap=" << vet_x1_ncap << std::endl;
  }
  if (vra_on && vsc_rot) {
    std::cout << "<rad_vet> vet_rotate_avg = " << vra_k << ": snapshot sweeps "
              << vra_nsweep << " (rotations " << vra_nrot << ", inits " << vra_ninit
              << "), " << vra_time << " s incl. offsets (rank 0, fenced)" << std::endl;
  }
  if (impl_hr_sf) {
    std::cout << "<rad_vet> implicit_hr_face=sf: cells with Jh > J (sf face ratios) "
              << hsf_nsf << " of " << hsf_ncell << " cell-solves (rank 0)" << std::endl;
  }
  if (vet_full) {
    std::cout << "<rad_vet> vet_tensor=full: guard triggered in " << vet_nguard << " of "
              << vet_ncell << " cell-calls (eigenvalue floor " << vet_eig_min
              << "); max |D_guarded - K/J|=" << vet_guard_max
              << " (incl. trace normalisation); min D_11=" << vet_d11_min << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplBeamInit
//! \brief boundary-injected collimated beams (the Jiang et al. 2012 /
//! Jiang 2021 beam and shadow tests).  <rad_vet>/implicit_beam_n = 1 or 2 (read only
//! when given); beam b (b = 0, 1): implicit_beam<b>_e (energy density E_b of the beam),
//! implicit_beam<b>_ang (in-plane angle in degrees from +x1 toward +x2, |ang| < 90),
//! implicit_beam<b>_y0 / _y1 (the x2 patch of the x1min face).  The beam is carried by
//! the up-going SC ray nearest (largest Omega . n_b) the in-plane direction n_b (ties:
//! the smallest |mu3|); in 2-D its mu3 mirror gets half of it, so F3 = 0.  The entering
//! intensity frac E_b / (n w_r) (frac = the covered fraction of the column's x2 extent,
//! n = 1 or 2 rays) is added to whatever the x1min face gives that ray (the Marshak bath,
//! 0 for vacuum), and the moment solve takes the same rays' face flux
//! F_in = c frac E_b mu1 as an imposed flux on those columns (Marshak elsewhere).
//! Needs closure = vet_sc, implicit_bc_x1min = marshak and ONE MeshBlock (VetShortChar;
//! several with scheme = clean).

void RadiationVet::ImplBeamInit(ParameterInput *pin) {
  ibm_n = 0;
  if (!pin->DoesParameterExist("rad_vet", "implicit_beam_n")) return;
  ibm_n = pin->GetInteger("rad_vet", "implicit_beam_n");
  if (ibm_n <= 0) {ibm_n = 0; return;}
  if (ibm_n > 2) VetFatal("<rad_vet>/implicit_beam_n must be 0, 1 or 2");
  if (!vet_sc) VetFatal("<rad_vet>/implicit_beam_n needs closure = vet_sc");
  if (ibc_x1min != M1_IBC_MARSHAK)
    VetFatal("<rad_vet>/implicit_beam_n needs implicit_bc_x1min = marshak");
  if (vet_x1per) VetFatal("<rad_vet>/implicit_beam_n is not for a periodic x1 sweep");
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx2 = indcs.nx2, nx3 = trans_x3 ? indcs.nx3 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
  const int nray = vet_nray;
  ibm_nr = 0;
  Kokkos::realloc(ibm_i, nmb, 4, nx3, nx2);
  Kokkos::realloc(ibm_f, nmb, nx3, nx2);
  Kokkos::realloc(ibm_s, nmb, 3, nx3, nx2);
  auto hi = Kokkos::create_mirror_view(ibm_i);
  auto hf = Kokkos::create_mirror_view(ibm_f);
  auto hs = Kokkos::create_mirror_view(ibm_s);
  Kokkos::deep_copy(hi, 0.0);
  Kokkos::deep_copy(hf, 0.0);
  Kokkos::deep_copy(hs, 0.0);
  auto &msz = pmy_pack->pmb->mb_size;
  for (int b = 0; b < ibm_n; ++b) {
    const std::string p = "implicit_beam" + std::to_string(b) + "_";
    const Real eb = pin->GetReal("rad_vet", p + "e");
    const Real ang = pin->GetReal("rad_vet", p + "ang")*M_PI/180.0;
    const Real y0 = pin->GetReal("rad_vet", p + "y0");
    const Real y1 = pin->GetReal("rad_vet", p + "y1");
    int rb = -1;
    Real best = 1.0e30;
    for (int r = 0; r < nray; ++r) {
      if (!(ah(r,0) > 0.0)) continue;
      // the ray nearest the in-plane beam direction (cos ang, sin ang, 0) in 3-D
      const Real d = -(ah(r,0)*std::cos(ang) + ah(r,1)*std::sin(ang))
                     + 1.0e-9*std::fabs(ah(r,2));
      if (d < best - 1.0e-12) {best = d; rb = r;}
    }
    int rs[2] = {rb, -1};
    int nr = 1;
    if (!trans_x3) {
      for (int r = 0; r < nray; ++r) {
        if (r != rb && std::fabs(ah(r,0) - ah(rb,0)) < 1.0e-12 &&
            std::fabs(ah(r,1) - ah(rb,1)) < 1.0e-12 &&
            std::fabs(ah(r,2) + ah(rb,2)) < 1.0e-12) {
          rs[1] = r;
          nr = 2;
          break;
        }
      }
    }
    for (int q = 0; q < nr; ++q) {
      const int qq = ibm_nr++;
      ibm_ray[qq] = rs[q];
      const Real wr = ah(rs[q],3);
      for (int m = 0; m < nmb; ++m) {
        // the beams enter through the x1min face of the mesh
        auto bf = pmy_pack->pmb->mb_bcs.h_view(m,BoundaryFace::inner_x1);
        if (bf == BoundaryFlag::block) continue;
        const Real x2a = msz.h_view(m).x2min, dx2 = msz.h_view(m).dx2;
        for (int k = 0; k < nx3; ++k) {
          for (int j = 0; j < nx2; ++j) {
            const Real ya = x2a + j*dx2, yb = ya + dx2;
            const Real fr = std::fmax(std::fmin(yb, y1) - std::fmax(ya, y0), 0.0)/dx2;
            if (!(fr > 0.0)) continue;
            hi(m,qq,k,j) += fr*eb/(nr*wr);
            hf(m,k,j) += c_light*wr*ah(rs[q],0)*fr*eb/(nr*wr);
            for (int a = 0; a < 3; ++a) {
              hs(m,a,k,j) += c_light*wr*ah(rs[q],a)*ah(rs[q],0)*fr*eb/(nr*wr);
            }
          }
        }
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "  implicit_beam" << b << ": E_b=" << eb << " requested ang="
                << ang*180.0/M_PI << " deg -> ray " << rb << " (mu1 mu2 mu3) = ("
                << ah(rb,0) << " " << ah(rb,1) << " " << ah(rb,2) << "), in-plane ang="
                << std::atan2(ah(rb,1), ah(rb,0))*180.0/M_PI << " deg, " << nr
                << " ray(s), x2 patch [" << y0 << ", " << y1 << "]" << std::endl;
    }
  }
  Kokkos::deep_copy(ibm_i, hi);
  Kokkos::deep_copy(ibm_f, hf);
  Kokkos::deep_copy(ibm_s, hs);
}

} // namespace radvet
