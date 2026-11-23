#ifndef RAD_M1_RAD_M1_CLEAN_HPP_
#define RAD_M1_RAD_M1_CLEAN_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean.hpp
//! \brief Data layout and device functions of the clean implicit radiation solve
//! (<rad_m1>/scheme = clean, see docs/dev/rad_clean_scheme_design.md).
//!
//! Each implicit stage solves for the radiation energy E in every cell, with the face
//! fluxes written in terms of E:
//!   F_f = (1-w) [theta (G_f - h c^2 grad(D E)) + (1-theta) A_f]
//!         + w c (h+_f E_L + h-_f E_R)
//! where theta = 1/(1 + h c chi_f) is the implicit factor of the flux source, A = a E the
//! advective (lab-frame) flux with a_d = v_d + D_db v_b, h+-_f the rays' half-range
//! fluxes at the face divided by J of the donor cell, and w a smooth function of the
//! rays' non-locality C1 = <|I - S|>/J.  The cell equation is
//!   E (1 + h c kE) + h div F = E_old + h c kP a T^4.

#include <cmath>
#include <functional>
#include <limits>
#include <string>

#include "athena.hpp"
#include "bvals/bvals.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "reconstruct/plm.hpp"

// the face flux is not inlined on the device: inlined into the kernels that call it
// many times, nvcc (cicc) compiles this layout for more than an hour.  Host builds
// inline it as before
#if defined(__CUDACC__)
#define CLEAN_DEVICE_NOINLINE __noinline__
#else
#define CLEAN_DEVICE_NOINLINE
#endif
// HIP: the face flux, stress and source always inlined.  hipcc called them (its inliner
// declined them): each kernel using them carried a 468 B/lane scratch frame (old
// production's M1_INL).  Where they are a path the pass does not take (the stored face
// terms of cln_rhs and cln_velocity) they stay calls (CleanFaceFluxCall, ...): inlined
// there they raised the register use of the path taken
#if defined(__HIPCC__)
#define CLEAN_HIP_INLINE __attribute__((always_inline))
#define CLEAN_HIP_NOINLINE __attribute__((noinline))
#else
#define CLEAN_HIP_INLINE
#define CLEAN_HIP_NOINLINE
#endif
// HIP: a kernel reads CleanArgs from its device copy (CleanDevArgs; a capture indexed at
// run time is copied to every lane's stack); CUDA keeps the capture (parameter space,
// faster on the A100)
#if defined(KOKKOS_ENABLE_HIP)
#define CLEAN_HIPX(...) __VA_ARGS__
#else
#define CLEAN_HIPX(...)
#endif

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

//! \fn void CleanMaxRanks
//! \brief the maxima of v[0..n-1] over all ranks in one reduction, a non-finite one as
//! inf (MPI_MAX of a NaN may differ between the ranks, whose Picard exits then disagree)
inline void CleanMaxRanks(Real *v, const int n) {
  for (int q = 0; q < n; ++q) {
    if (!std::isfinite(v[q])) v[q] = std::numeric_limits<Real>::infinity();
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
}

// CleanFor: a 4-D clean pass; on CUDA at kCleanForMb blocks of 256 threads per SM (a
// register cap of 65536/(256 kCleanForMb)), with the flattening and index arithmetic of
// athena's par_for (bitwise).  Measured on R136 (20 cycles, hst identical): GH200 MB 3
// -19 ms/cycle (MB 2 -7), A100 MB 2 -20 ms/cycle (MB 3 +14, MB 1 +29).  HIP and the host
// keep athena's par_for; on HIP a pass may ask for kHipMb blocks (kHipMb 4: a cap of 128
// VGPRs, 4 waves per SIMD, for the face passes that fit it without spilling)
#if defined(ATHENA_CLEAN_FOR_MB)   // build-time override (register cap only, bitwise)
constexpr int kCleanForMb = ATHENA_CLEAN_FOR_MB;
#elif defined(KOKKOS_ARCH_HOPPER90)
constexpr int kCleanForMb = 3;
#else
constexpr int kCleanForMb = 2;
#endif
// (kCudaMb: a pass's own CUDA blocks per SM, measured per kernel on A100 R136 (spdc):
// cln_fwc, cln_stencil, cln_beta at 3 (-2.3, -1.0, -0.4 ms/cycle))
template <int kHipMb = 0, int kCudaMb = kCleanForMb, typename F>
inline void CleanFor(const char *name, const int nmb1, const int ks, const int ke,
                     const int js, const int je, const int is, const int ie,
                     const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int mb = kHipMb;
#elif defined(KOKKOS_ENABLE_CUDA)
  constexpr int mb = kCudaMb;
#else
  constexpr int mb = 0;
#endif
  if constexpr (mb > 0) {
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    AthFor(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, mb>>(
           DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/ni;
      f(m, ks + k, js + j, is + idx - m*nkji - k*nji - j*ni);
    });
  } else {
    par_for(name, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie, f);
  }
}

// on HIP, CleanFor with k fastest, then j, then i (the order of the rows in a VgdRag slab
// of vgd_i: d fastest, then k, then j), so that adjacent lanes read adjacent intensity
// rows, as prod's m1_vgd_mom; MI300A R136 -62 ms/cycle, A100 +6 (CUDA and the host keep
// CleanFor).  For passes whose cells are independent (bitwise)
template <typename F>
inline void CleanForK(const char *name, const int nmb1, const int ks, const int ke,
                      const int js, const int je, const int is, const int ie,
                      const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  par_for(name, DevExeSpace(), 0, nmb1, is, ie, js, je, ks, ke,
  KOKKOS_LAMBDA(const int m, const int i, const int j, const int k) {f(m, k, j, i);});
#else
  CleanFor(name, nmb1, ks, ke, js, je, is, ie, f);
#endif
}

// a pass of columns (one thread per (m, k, j), marching in i): on HIP blocks
// of 64 with k fastest (blocks of 256 put the HG level's 12k columns on 48 of the 228
// MI300A CUs); CUDA and the host keep par_for.  The same arithmetic per column
template <typename F>
inline void CleanForCol(const char *name, const int nmb1, const int ks, const int ke,
                        const int js, const int je, const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  const int nk = ke - ks + 1, njk = (je - js + 1)*nk;
  Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<64,
                       1>>(DevExeSpace(), 0, (nmb1 + 1)*njk),
  KOKKOS_LAMBDA(const int idx) {
    const int m = idx/njk, j = (idx - m*njk)/nk;
    f(m, ks + idx - m*njk - j*nk, js + j);
  });
#else
  par_for(name, DevExeSpace(), 0, nmb1, ks, ke, js, je, f);
#endif
}

// closure variables (arrays cln_c and cln_s); cln_s holds the first CLN_FH (cubed
// sphere: CLN_F1)
constexpr int CLN_D11 = 0;   // Eddington tensor D11 D22 D33 D12 D13 D23
constexpr int CLN_HP = 6;    // h+_a = H+_a/J, a = 1,2,3
constexpr int CLN_HM = 9;    // h-_a = H-_a/J (<= 0)
constexpr int CLN_C1 = 12;   // non-locality C1 = <|I - S|>/J
constexpr int CLN_J = 13;    // J of the rays (energy-density units)
constexpr int CLN_F2 = 14;   // the rays' flux through the low x2 face (M1_VET_F2, not
                             // normalised; clean_consistency)
constexpr int CLN_FX = 18;   // their flux through the low x1 face (M1_VET_FX)
constexpr int CLN_F3 = 22;   // their flux through the low x3 face (M1_VET_F3)
constexpr int CLN_G = 26;    // curvature of their stress in the cell layout (y_r,
                             // y_theta)
constexpr int CLN_JR = 28;   // their J also below clean_ray_min E (CLN_J is 0 there): the
                             // staggered flux factors (clean_stag_gflux)
constexpr int CLN_FH = 29;   // cubed sphere: through the high x2, x3 face of the last
                             // cell of a MeshBlock (+4; a seam ghost's low face differs)
constexpr int CLN_F1 = 37;   // cell layout: their flux through an x1 boundary face
                             // (M1_VET_F1)
constexpr int CLN_KP = 41;   // cell layout: k+^d_b = K+^d_bd/J (outgoing stress of
                             // mu_d > 0), 3 d + b
constexpr int CLN_N = 50;
// clean_ring_fh (cubed sphere, clean_rays = sc): in place of k+ (not held there) each
// cell's half ranges about its own faces' normals, h+ out of its high faces (3), h- of
// its low faces (3) (CleanCellHalf), so that the ring carries them as cell moments
constexpr int CLN_HO = CLN_KP;
constexpr int CLN_NHO = CLN_HO + 6;

// the slots of cln_s
inline int CleanNs(const bool cub, const bool ho = false) {
  return cub ? (ho ? CLN_NHO : CLN_F1) : CLN_FH;
}

// the moments of the cell (D, h+-, C1, J, k+; the others are face sums)
KOKKOS_INLINE_FUNCTION
bool CleanCellMoment(const int n) {return n <= CLN_J || n >= CLN_KP;}

// the slot of channel n of the restart files, which keep push 16's order (D h+ h- k+ C1
// J F2 F1 FX F3 FH G JR)
inline int CleanRstSlot(const int n) {
  constexpr int o[9] = {0, 12, 21, 23, 27, 31, 39, 47, CLN_N};
  constexpr int s[8] = {CLN_D11, CLN_KP, CLN_C1, CLN_F2, CLN_F1, CLN_FX, CLN_FH, CLN_G};
  int g = 0;
  while (n >= o[g+1]) ++g;
  return s[g] + n - o[g];
}

// face variables (arrays cln_fc1, cln_fc2, cln_fc3)
constexpr int CLF_W = 0;     // blend weight w
constexpr int CLF_TH = 1;    // theta = 1/(1 + h c chi_f)
constexpr int CLF_G = 2;     // flux of the old vector G_f
constexpr int CLF_HP = 3;    // h+_f
constexpr int CLF_HM = 4;    // h-_f
constexpr int CLF_N = 5;     // (the shared ones)
// staggered face slots after CLF_N (clean_stag_fsolve: CLF_N .. CLF_N + 3 hold its F,
// the dissipation term, F^H and, with clean_stag_fkin, the energy flux's memory d)
constexpr int CLF_QL = CLF_N + 4;   // the operator's al, ar (CleanFaceWeights) and cross
                                    // factor -(1-w) theta h c^2 of the pass (3)
constexpr int CLF_PT = CLF_QL;      // outside the linear solve (cln_fwc sets QL after
                                    // the RHS): the faces' momentum terms of the velocity
                                    // update (CleanStagDrag) and the RHS's a E
                                    // (cln_rhsaf)
constexpr int CLF_SR = CLF_QL + 1;  //   and the RHS's face sources (CleanFaceSource)
constexpr int CLF_DS = CLF_N + 7;   // clean_stag_fkin: d of the stage-loop start
constexpr int CLF_DB = CLF_N + 8;   // clean_stag_vimp: drag share b = beta/(1 + beta) and
constexpr int CLF_DQ = CLF_N + 9;   //   the face momentum Q = G - c^2 dm_f of the iterate
constexpr int CLF_KD = CLF_N + 10;  // clean_stag_fkin: the slopes K1, K2 of d (2)
constexpr Real CLN_KPOS_TH = 0.1;   // clean_stag_kpos: E_old at least this share of E_s
constexpr int CLN_KPOS_RN = 12;     //   rounds at most, the inflow-free share from
constexpr int CLN_KPOS_RF = 6;      //   round CLN_KPOS_RF
constexpr int CLN_KPOS_REDO = 4;    // clean_stag_kpos = redo: rounds, the last all faces
constexpr Real CLN_MEMUP = 1.0/6.0;  // clean_stag_gflux: a memory outflow above this
                                     // share of the donor's E_old goes upwind

// work vectors (array cln_w): the first CLW_NS are shared, then the staggered layout's
// own (CLW_NS..CLW_N - 1)
constexpr int CLW_EO = 0;    // E of the old vector
constexpr int CLW_EK = 1;    // Newton iterate of E
constexpr int CLW_T = 2;     // gas temperature of the iterate
constexpr int CLW_EI = 3;    // gas internal energy of the old vector
// exchange linearisation r = s/(d cv + s), s = 4 h c kP a T^3
constexpr int CLW_RR = 4;
constexpr int CLW_A1 = 5;    // advective coefficients a_1, a_2, a_3
constexpr int CLW_DM = 8;    // momentum change of the iterate (3)
constexpr int CLW_B = 11;    // right-hand side
constexpr int CLW_X = 12;    // solution
constexpr int CLW_BT = 13;   // beta = v/c of the iterate (3)
constexpr int CLW_RE = 16;   // h S_E = hc kP G a T^4 - hc kE E + RE E + RY.y - W (Ma et
constexpr int CLW_GM = 17;   //   al. 2025, eq. 56, G = CLW_GM): the beta terms (RY: 3);
                             //   staggered: RE = hc kP a T^4, RY 1 - r
constexpr int CLW_EP = 18;   // -(RE E + RY.y) of the solution (eint of the gas)
constexpr int CLW_JV = 19;   // (1 - r) dA/dbeta/(c^2 rho), dS_p/dbeta/(c^2 rho) of the
constexpr int CLW_GX = 22;   //   iterate (3 + 1): implicit gas velocity;
                             //   staggered clean_stag_sbimp: GX = r (CleanRightHandSide);
                             //   staggered gflux: JV + 1 = hc kE dop (CleanNewtonUpdate)
constexpr int CLW_RY = 23;
constexpr int CLW_NS = 24;
constexpr int CLW_AM = 24;   // staggered layout: x1-line coefficient of E_{i-1}
constexpr int CLW_A0 = 25;   // diagonal
constexpr int CLW_AP = 26;   // x1-line coefficient of E_{i+1}
constexpr int CLW_CP = 27;   // tridiagonal solve scratch
constexpr int CLW_R = 28;    // BiCGStab vectors
constexpr int CLW_R0 = 29;
constexpr int CLW_P = 30;
constexpr int CLW_V = 31;
constexpr int CLW_S = 32;
constexpr int CLW_Q = 33;
constexpr int CLW_Y = 34;
constexpr int CLW_Z = 35;
constexpr int CLW_TM = 36;   // clean_stag_pc = mg: coefficients of E_{j-1}, E_{j+1},
                             //   E_{k-1}, E_{k+1} (4)
constexpr int CLW_L = 40;    // staggered, cubed sphere: the lateral stress sums Y of the
                             // face families q of a cell, CLW_L + 3 q
                             // (CleanStagCsLateral)
constexpr int CLW_SE = 49;   // staggered clean_conv = budget: the gross terms S of the E
constexpr int CLW_SG = 50;   //   and the gas rows (their scales), the weighted B and a
constexpr int CLW_BW = 51;   //   preconditioner input (CleanBiCGStab)
constexpr int CLW_PW = 52;
constexpr int CLW_N = 53;
constexpr int CLW_KF = 53;   // staggered clean_stag_kpos: E of the old vector as formed,
constexpr int CLW_KM = 54;   //   scratch (the new r of a round, the limited E_old),
constexpr int CLW_PD = 55;   //   the limiter r (setup); clean_stag_xpat: the rows'
                             //   Patankar diagonal
constexpr int CLW_NK = 56;   //   (allocated only with clean_stag_kpos / xpat)
constexpr int CLW_XA = 56;   // staggered clean_xt_exact: B less its rest (the row's
constexpr int CLW_XW = 57;   //   exchange as linearised), the row's work W
constexpr int CLW_NX = 58;   //   (allocated only with clean_xt_exact)
// staggered clean_stag_kpos / kreal: the momentum increment's
constexpr int CLW_KS = 58;
                             //   applied share s (CleanOldVector, kreal; 1 without),
constexpr int CLW_KT = 59;   //   a donor's faces moved: 1 (redo: 2 the new ones; with
                             //   ghosts; CleanStagKPos, CleanStagKRedo)
constexpr int CLW_NKP = 60;  //   (allocated with clean_stag_kpos)
constexpr int CLW_KR = 60;   // staggered clean_stag_kreal: the deferred kick (1 - s) D of
constexpr int CLW_NR = 63;   //   the stage's gas momentum in the frame of y (3)

//----------------------------------------------------------------------------------------
//! \struct CleanArgs
//! \brief everything the device functions of the clean solve read

struct CleanArgs {
  DvceArray5D<Real> cs;        // closure of the stage
  CleanW w;                    // work vectors
  DvceArray5D<Real> fc[3];     // face variables in x1, x2, x3
  DvceArray3D<Real> beam;      // boundary beam flux at x1min (m, k-ks, j-js)
  DvceArray4D<Real> bstress;   // beam stress column c K_b1 at x1min (m, b, k-ks, j-js)
  Real h, c, idx[3];           // stage step, light speed, inverse cell widths
  int bc[6];                   // boundary types x1min, x1max, ..., x3max
  DvceArray2D<BoundaryFlag> mbc;  // boundary flags of the six faces of each MeshBlock
  Real ebath[6];               // boundary bath energy densities
  Real iflux[2];               // imposed fluxes of x1 flux faces (M1_IBC_FLUX), x1min/max
  DvceArray3D<Real> ifx2[2];   // those of x2 flux faces (m, k, i), x2min/max
  Real marshak_q;
  bool marshak;                // Marshak boundary form (else the rays' half range)
  bool use_beam;
  bool feedback;               // gas coupling (advective flux, exchange linearisation)
  bool vzero;                  // the gas momentum is frozen (cln_vzero, steady start)
  bool fd;                     // staggered, clean_stag_fkin: the energy flux's
                               // G + F^H - F
  bool gf;                     // staggered, clean_stag_gflux: flux factors, low order
  // cubed sphere, clean_stag_csvec_radial (CleanStagCsVector)
  bool cvr;
  bool fre;                    // staggered, clean_stag_freal: a realisable comoving sb
  int fdg;                     // staggered bracket: face slot of its lagged flux DG (-1)
  bool fs;                     // staggered, clean_stag_fsolve: the face F of the F solve
  bool tf;                     //   and each face's momentum exchange T of its owner block
                               // (CLF_N) with its dissipation term (CLF_N + 1)
  bool hi;                     // staggered, clean_stag_split: the face momentum state
                               // F^H (w = 0), not the energy flux
  int xq;                      // staggered, clean_stag_xdonor: face slots qL, qR of the
                               // lagged fluxes' donor part, their explicit part, C(E_k)
                               // (-1: the central cross flux)
  bool pd;                     // staggered, kpos / xpat: the rows' diagonal CLW_PD
  Real kt;                     // staggered, clean_stag_kpos = delta: tau/h, the slopes'
                               // linear face flux tau (al E_L + ar E_R) implicit
  int kq;                      // staggered, clean_stag_kpos: face slots qL, qR of the
                               // moved increment flux (-1: none; the E rows only)
  bool xlag, xsrc;             //   lag: explicit up to the donor's share; xsrc: the
                               //   face sources too
  bool mcl;                    // staggered: a face's momentum term at most the force of
                               // a realisable flux, h chi_f E (the iterate's,
                               // CleanMomentumChange)
  DvceArray1D<Real> mclc;      // (mcl: the face terms the cap moved, counted)
  bool fa;                     // staggered: the faces' momentum terms of the pass stored
                               // (CLF_PT, CleanStagDrag)
  bool drag;                   // staggered, clean_stag_vimp: implicit gas velocity
  bool upw;                    // staggered, clean_stag_upwind: donor-cell advective flux
  bool ucen;                   //   clean_stag_upwind_central: its lagged rest central
  bool fwc;                    // staggered: al, ar and the cross factor from CLF_QL
  bool sbi;                    // staggered, clean_stag_sbimp: E-row beta terms implicit
  bool sar;                    // staggered: the streaming share's net flux in area form
  bool bar;                    // clean_bc_area: an open x1 end's h+- E in area form too
                               // (CleanStagArea)
  DvceArray5D<Real> op;        // opacities (M1_OP_*)
  int nf;                      // number of flux components (dimensions)
  int is, ie, js, je, ks, ke;
  bool multi_d, three_d;
  // spherical-polar geometry: x1 = r, x2 = theta, x3 = phi, vectors in the local
  // orthonormal basis
  bool sph;
  DvceArray4D<Real> vol;       // cell volumes
  DvceArray4D<Real> area[3];   // face areas normal to x1, x2, x3
  DvceArray4D<Real> dxf[3];    // distances between the cell centres across those faces
  DvceArray2D<Real> x1v, x1f;  // cell-centre (volume centroid) and face radii
  DvceArray2D<Real> x2v, x2f;  // cell-centre and face polar angles
  DvceArray3D<Real> tg;         // tan, tan, sin, sin of x2v, x2f, x2v, x2f by (m, q, j)
  Real cseam;                  // cot of a theta-periodic seam face: its two copies' mean
  // cubed sphere: y in the cell's orthonormal frame (r, theta, phi) about the hydro's
  // radius vector; vol, area, dxf and x1v, x1f as above, g the geometry (CLG_*) of
  // CleanCsGeometry
  bool cub;
  bool tdn = false;   // clean_thin_donor (CleanStagFRay)
  bool lwb = false;   // clean_cs_lat_wb (CleanStagCsLat)
  bool ys = false;    // clean_stag_lat_sums: the lateral sums Y off the cubed sphere
  bool lzbf = false;  // clean_lz_bnd_free on spherical polar (CleanFaceStress)
  DvceArray5D<Real> yl;   //   (cln_yl: CLW_L + 3 q on the cubed sphere)
  DvceArray4D<Real> ae0;   // (the cubed sphere: pcoord's radial face area, CleanDivW)
  CleanGeo g;
  // the conserved gas energy also holds d phi (<hydro>/etotgrav) and |B|^2/2 (<mhd>)
  bool etg, mhd;
  DvceArray4D<Real> phi, emag;
};

// the rows of the staggered operator (clean_stag_stencil): the coefficients of E at the
// cell and its 18 neighbours that share a face or an edge.  Slot of offset (dk, dj, di):
// 0 the cell; 1 + 2 a + (d > 0) the faces along axis a (i, j, k); 7 + 4 p + 2 (d1 > 0) +
// (d2 > 0) the edges in the plane p of axes (j, i), (k, i), (k, j), offsets d1, d2
constexpr int CLN_ST = 19;
// cubed sphere, clean_stag_stencil: the coefficients of the lateral stress sums Y_q,b of
// a cell (CleanStagCsLateral) on E of its neighbour across face s of family q (slot 3 (2
// q + s) + b) and on its own E (18 + 3 q + b)
constexpr int CLN_CSY = 27;

KOKKOS_INLINE_FUNCTION
constexpr int CleanStSlot(const int dk, const int dj, const int di) {
  return ((dk != 0) + (dj != 0) + (di != 0) == 0) ? 0 :
         ((dk != 0) + (dj != 0) + (di != 0) == 1) ?
         1 + 2*((dk != 0) ? 2 : (dj != 0)) + (dk + dj + di > 0) :
         7 + 4*((dk == 0) ? 0 : ((dj == 0) ? 1 : 2)) + 2*(((dk == 0) ? dj : dk) > 0) +
         (((di == 0) ? dj : di) > 0);
}

KOKKOS_INLINE_FUNCTION
void CleanStOff(const int s, int &dk, int &dj, int &di) {   // (the offset of slot s)
  dk = dj = di = 0;
  if (s > 0 && s < 7) {
    const int d = ((s - 1) % 2) ? 1 : -1;
    ((s < 3) ? di : ((s < 5) ? dj : dk)) = d;
  } else if (s >= 7) {
    const int p = (s - 7)/4, d1 = (((s - 7)/2) % 2) ? 1 : -1, d2 = ((s - 7) % 2) ? 1 : -1;
    ((p == 0) ? dj : dk) = d1;
    ((p == 2) ? dj : di) = d2;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanGasEint
//! \brief gas internal energy of cell (k,j,i) of the conserved state uh with kinetic
//! energy ekin: the total energy less ekin, the gravitational energy d phi
//! (<hydro>/etotgrav) and the magnetic energy (<mhd>), as RadiationM1::Opacity forms it

KOKKOS_INLINE_FUNCTION
Real CleanGasEint(const CleanArgs &x, const DvceArray5D<Real> &uh, const Real ekin,
                  const int m, const int k, const int j, const int i) {
  Real eint = uh(m,IEN,k,j,i) - ekin;
  if (x.etg) eint -= uh(m,IDN,k,j,i)*x.phi(m,k,j,i);
  if (x.mhd) eint -= x.emag(m,k,j,i);
  return eint;
}

//----------------------------------------------------------------------------------------
//! \fn int CleanBc
//! \brief boundary type of side s (x1min, x1max, ..., x3max) of MeshBlock m: the mesh's
//! on a physical boundary, M1_IBC_PERIODIC on a face shared with a MeshBlock (or a
//! cubed-sphere panel)

KOKKOS_INLINE_FUNCTION
int CleanBc(const CleanArgs &x, const int m, const int s) {
  const BoundaryFlag f = x.mbc(m,s);
  return (f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
          f == BoundaryFlag::panel) ? M1_IBC_PERIODIC : x.bc[s];
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanOpen
//! \brief true if the physical boundary of side s of MeshBlock m lets the rays' flux
//! through (bath, beam): neither a wall nor an imposed flux

KOKKOS_INLINE_FUNCTION
bool CleanOpen(const CleanArgs &x, const int m, const int s) {
  const int bc = CleanBc(x, m, s);
  return bc != M1_IBC_REFLECT && bc != M1_IBC_FLUX;
}

//----------------------------------------------------------------------------------------
//! \struct CleanShell
//! \brief the cells of an n3 x n2 x n1 array outside its box [k0,k1] x [j0,j1] x [i0,i1]
//! (CleanInner by g), counted plane by plane: the k planes below and above the box whole,
//! the planes through it without the box's rows

struct CleanShell {
  int n1, n2, n3, i0, i1, j0, j1, k0, k1;
  int nio, ncell;
  CleanShell(const int n3_, const int n2_, const int n1_, const int g, const int is,
             const int ie, const int js, const int je, const int ks, const int ke,
             const bool md, const bool td) :
      n1(n1_), n2(n2_), n3(n3_), i0(is + g), i1(ie - g), j0(md ? js + g : 0),
      j1(md ? je - g : n2_ - 1), k0(td ? ks + g : 0), k1(td ? ke - g : n3_ - 1) {
    nio = n1 - (i1 - i0 + 1);
    ncell = n3*n2*n1 - (k1 - k0 + 1)*(j1 - j0 + 1)*(i1 - i0 + 1);
  }
  bool Empty() const {return i1 < i0 || j1 < j0 || k1 < k0;}
};

//----------------------------------------------------------------------------------------
//! \fn void CleanShellCell
//! \brief cell (k,j,i) number q of the shell s

KOKKOS_INLINE_FUNCTION
void CleanShellCell(const CleanShell &s, int q, int &k, int &j, int &i) {
  const int np = s.n2*s.n1, nj = s.j1 - s.j0 + 1;
  const int pl = np - nj*(s.i1 - s.i0 + 1), nb = s.k0*np, nm = (s.k1 - s.k0 + 1)*pl;
  if (q < nb || q >= nb + nm) {
    if (q >= nb) q -= nb + nm - (s.k1 + 1)*np;
    k = q/np;
    j = (q - k*np)/s.n1;
    i = q - k*np - j*s.n1;
    return;
  }
  q -= nb;
  k = s.k0 + q/pl;
  int r = q%pl;
  const int rb = s.j0*s.n1, rm = nj*s.nio;
  if (r < rb || r >= rb + rm) {
    if (r >= rb) r -= rb + rm - (s.j1 + 1)*s.n1;
    j = r/s.n1;
    i = r%s.n1;
    return;
  }
  r -= rb;
  j = s.j0 + r/s.nio;
  const int t = r%s.nio;
  i = (t < s.i0) ? t : s.i1 + 1 + (t - s.i0);
}

//----------------------------------------------------------------------------------------
//! \struct CleanLevel
//! \brief one grid level of the staggered multigrid (CleanStagLevel): an array holding
//! its blocks and vectors, with their slots

struct CleanLevel {
  CleanW u;                    // blocks and vectors (m, slot, k, j, i)
  int o[3];                    // index of the first active cell along x1, x2, x3
  int n[3];                    // active cells along x1, x2, x3
  int a0;                      // diagonal block
  int sm[3], sp[3];            // blocks of the lower and upper neighbour along each dir
};

// slots of a coarse level of the staggered mg (arrays cln_mg): the x1-line coefficients,
// diagonal and Thomas scratch in the order of cln_w (CLW_AM..CLW_CP), then as CLW_TM
constexpr int CLS_A0 = 1;    // diagonal (CLS_A0 -+ 1: E_{i-1}, E_{i+1}; + 2: scratch)
constexpr int CLS_TM = 4;    // E_{j-1}, E_{j+1}, E_{k-1}, E_{k+1}
constexpr int CLS_B = 8;     // right-hand side
constexpr int CLS_X = 9;     // solution
constexpr int CLS_N = 10;

//----------------------------------------------------------------------------------------
//! \fn Real CleanDuu
//! \brief u.D.u of cln_s's Eddington tensor (Cartesian on the cubed sphere) along the
//! unit vector u = g(m, gs..gs+2, k,j,i) of cln_gn (CLG_R, CLG_E, CLG_E + 3: the hydro's
//! directions), for <hydro>/rad_signal_speed

KOKKOS_INLINE_FUNCTION
Real CleanDuu(const DvceArray5D<Real> &s, const CleanGeo &g, const int gs,
              const int m, const int k, const int j, const int i) {
  Real d = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      d += g(m,gs+a,k,j,i)*g(m,gs+b,k,j,i)*s(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSmoothStep
//! \brief 3x^2 - 2x^3 on [0,1], 0 below and 1 above

KOKKOS_INLINE_FUNCTION
Real CleanSmoothStep(Real x) {
  x = fmin(fmax(x, 0.0), 1.0);
  return x*x*(3.0 - 2.0*x);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanExchangeTemperature
//! \brief the gas temperature of a cell after an implicit local exchange with radiation:
//! the root of e(d,T) + A T^4 = rhs, with A = h c kP a/(1 + h c kE) folded in by the
//! caller.  The left side increases monotonically with T, so the root is unique and lies
//! below T(d, rhs).  Safeguarded Newton iteration, general EOS.  rtol > 0
//! (clean_troot_rtol): the root also stands once the residual is at that level of the
//! row's terms, |g| <= rtol (e + A T^4 + rhs), the noise of a tabulated e(T) lying above
//! a step of 1e-15 T; the bound T(d, rhs) (an EOS inversion) then only when a step
//! leaves the bracket

KOKKOS_INLINE_FUNCTION
Real CleanExchangeTemperature(const EOS_Data &eos, const Real d, const Real a,
                              const Real rhs, const Real tguess, const Real rtol = 0.0) {
  if (!(rhs > 0.0)) return 0.0;
  Real lo = 0.0, hi = (rtol > 0.0 && tguess > 0.0) ? -1.0 : eos.Temperature(d, rhs);
  Real t = (tguess > 0.0 && (hi < 0.0 || tguess < hi)) ? tguess : 0.5*hi;
  for (int n = 0; n < 100; ++n) {
    Real e, p, chi_rho, chi_t, cv;
    eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
    const Real t3 = t*t*t;
    const Real g = e + a*t3*t - rhs;
    if (rtol > 0.0 && fabs(g) <= rtol*(fabs(e) + a*t3*t + rhs)) {
      return t - g/(d*cv + 4.0*a*t3);
    }
    if (g > 0.0) {
      hi = t;
    } else {
      lo = t;
    }
    Real tn = t - g/(d*cv + 4.0*a*t3);
    if (!(tn > lo && (hi < 0.0 || tn < hi))) {
      if (hi < 0.0) hi = eos.Temperature(d, rhs);
      tn = 0.5*(lo + hi);
    }
    if (fabs(tn - t) <= 1.0e-15*t) return tn;
    t = tn;
  }
  return t;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagSph
//! \brief staggered layout on spherical polar (x1 = r, x2 = theta; vectors in the local
//! orthonormal basis): q = 0 the inverse centre distance across face (k,j,i) normal to
//! dir, q = 1 its inverse radius (r_f on r faces, r of the cell on theta faces), q = 2, 3
//! the curvature coefficient of the face momentum equation of cell L or R: the stress
//! divergence is d_n(D_nn E) + transverse + (h_L E_L + h_R E_R)/(2 r) with
//!   r faces:     h = 2 D_rr - D_thth - D_phph + cot(theta) D_rth  (theta of the cell)
//!   theta faces: h = cot(theta_f) (D_thth - D_phph) + 3 D_rth
//!   phi faces:   none: the stress divergence is in angular-momentum (torque) form,
//!                [sum over the dual faces of A R (D E)_.ph + A_ph ((D E)_phph)]/(V R),
//!                R = r sin(theta), so that L_z changes only by boundary torques
//! (a uniform isotropic state, D = I/3, is exactly steady); q = 4 the part of q = 2, 3
//! of the normal stress itself per D_nn E, s = 1/r_f (r faces), cot(theta_f)/(2 r)
//! (theta faces), so that a normal stress P of the cells enters as (1/dx +- s) P_R,L,
//! the discrete (1/A) d_n(A P)

KOKKOS_INLINE_FUNCTION
Real CleanStagSph(const CleanArgs &x, const int q, const int dir, const int m,
                  const int k, const int j, const int i) {
  if (q == 0) {
    return (dir == 2) ? x.area[2](m,k,j,i)/x.vol(m,k,j,i) : 1.0/x.dxf[dir](m,k,j,i);
  }
  if (dir == 2 && q >= 2) return 0.0;   // (phi faces: torque form, CleanFaceStress)
  // (a theta-periodic seam face takes one cot in both copies, as 32fab890f, so that its
  // flux is single valued)
  const BoundaryFlag pf = BoundaryFlag::periodic;
  const bool sm = (dir == 1) &&
                  ((j == x.js && x.mbc(m,BoundaryFace::inner_x2) == pf) ||
                   (j == x.je + 1 && x.mbc(m,BoundaryFace::outer_x2) == pf));
  if (q == 4) {
    return (dir == 0) ? 1.0/x.x1f(m,i) :
           (sm ? 0.5*x.cseam/x.x1v(m,i) : 0.5/(x.tg(m,1,j)*x.x1v(m,i)));
  }
  if (q == 1) return 1.0/((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
  const int kc = (q == 2) ? k - (dir == 2) : k, jc = (q == 2) ? j - (dir == 1) : j;
  const int ic = (q == 2) ? i - (dir == 0) : i;
  auto d = [&](const int n) {return x.cs(m,CLN_D11+n,kc,jc,ic);};
  if (dir == 0) {
    return 2.0*d(0) - d(1) - d(2) + (x.multi_d ? d(3)/x.tg(m,0,jc) : 0.0);
  }
  if (sm) return (d(1) - d(2))*x.cseam + 3.0*d(3);
  return (d(1) - d(2))/x.tg(m,1,j) + 3.0*d(3);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDivW
//! \brief weight of the flux through the low (side 0) or high (side 1) face along dir in
//! the divergence of cell (k,j,i): 1/dx_d, on spherical polar A_f/V; on the cubed sphere
//! a radial face's curved area r^2 Omega (pcoord's, the hydro's: 4 pi r^2 over a shell),
//! not that of the closed vector areas (r^2 |S|, CleanCsGeometry, which the momentum and
//! the rays' face sums take; 7e-4 less over the sphere at N = 16)

KOKKOS_INLINE_FUNCTION
Real CleanDivW(const CleanArgs &x, const int dir, const int m, const int k, const int j,
               const int i, const int side) {
  if (!x.sph && !x.cub) return x.idx[dir];
  const auto &a = (x.cub && dir == 0) ? x.ae0 : x.area[dir];
  return a(m,k+side*(dir == 2),j+side*(dir == 1),i+side*(dir == 0))/x.vol(m,k,j,i);
}

//! \brief CleanDivW's divergence of the low and high face values fl, fh along dir
//! (Cartesian (fh - fl)/dx), and the volume of cell (k,j,i) (Cartesian dx dy dz)

KOKKOS_INLINE_FUNCTION
Real CleanDivF(const CleanArgs &x, const int dir, const int m, const int k, const int j,
               const int i, const Real fl, const Real fh) {
  if (!x.sph && !x.cub) return x.idx[dir]*(fh - fl);
  return CleanDivW(x, dir, m, k, j, i, 1)*fh - CleanDivW(x, dir, m, k, j, i, 0)*fl;
}

KOKKOS_INLINE_FUNCTION
Real CleanVol(const CleanArgs &x, const int m, const int k, const int j, const int i) {
  return (x.sph || x.cub) ? x.vol(m,k,j,i) : 1.0/(x.idx[0]*x.idx[1]*x.idx[2]);
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanCsSeamFace
//! \brief cubed sphere: true for a panel seam face (k,j,i) normal to dir, except next to
//! a cube vertex (production's foot-only vertex)

KOKKOS_INLINE_FUNCTION
bool CleanCsSeamFace(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i) {
  if (!x.cub || dir == 0) return false;
  const BoundaryFlag pf = BoundaryFlag::panel;
  const int b = 3 - dir, p = (dir == 1) ? j : k, o = (dir == 1) ? k : j;
  const int pl = (dir == 1) ? x.js : x.ks, ph = ((dir == 1) ? x.je : x.ke) + 1;
  const int ol = (dir == 1) ? x.ks : x.js, oh = (dir == 1) ? x.ke : x.je;
  return ((p == pl && x.mbc(m,2*dir) == pf) || (p == ph && x.mbc(m,2*dir+1) == pf)) &&
         !((o == ol && x.mbc(m,2*b) == pf) || (o == oh && x.mbc(m,2*b+1) == pf));
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsNormA
//! \brief cubed sphere: the area A_n that weights the normal stress difference of the
//! face-normal momentum of face (k,j,i) normal to dir over V_f, A_n (P_R - P_L)/V_f.  It
//! is the face area, except on a panel seam face (CleanCsSeamFace): there the two cells
//! are mirror images, their O(h) volume variation adds in V_f instead of cancelling, and
//! A_n = V_f/d with d the chord of the mirror pair's cell centres (where the cell values
//! sit), d = 2 |X| |r_hat.n_f| of the cell on this side (production's two-point seam row
//! over r angm, M1CsRowGeo / M1CsSeamGeom)

KOKKOS_INLINE_FUNCTION
Real CleanCsNormA(const CleanArgs &x, const int dir, const int m, const int k,
                  const int j, const int i) {
  if (!CleanCsSeamFace(x, dir, m, k, j, i)) return x.area[dir](m,k,j,i);
  const int kl = k - (dir == 2), jl = j - (dir == 1);
  const bool lo = ((dir == 1) ? j : k) == ((dir == 1) ? x.js : x.ks);
  const int ko = lo ? k : kl, jo = lo ? j : jl;
  Real xr = 0.0, rn = 0.0;
  for (int c = 0; c < 3; ++c) {
    xr += SQR(x.g.At<CLG_X>(m,c,ko,jo,i));
    rn += x.g.At<CLG_R>(m,c,ko,jo,i)*x.g.At<CLG_N>(m,3*dir+c,k,j,i);
  }
  return 0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,i))/(2.0*std::sqrt(xr)*fabs(rn));
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagDnn
//! \brief n.D.n of cell (k,j,i), n the unit normal of face (kf,jf,iff) normal to dir,
//! the normal stress of the face-normal momentum: D_dir,dir on Cartesian and spherical
//! polar meshes (n = e_dir), on the cubed sphere with D Cartesian

KOKKOS_INLINE_FUNCTION
Real CleanStagDnn(const CleanArgs &x, const int dir, const int m, const int kf,
                  const int jf, const int iff, const int k, const int j, const int i) {
  if (!x.cub) return x.cs(m,CLN_D11+dir,k,j,i);
  Real d = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      d += x.g.At<CLG_N>(m,3*dir+a,kf,jf,iff)*x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff)*
           x.cs(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \struct CleanFaceMet
//! \brief the metric of the momentum row of face (k,j,i) normal to dir on every mesh: the
//! normal stress difference of the cells is a (P_R - P_L)/v + s (P_R + P_L), plus (curved
//! meshes) hr (h_L P_L + h_R P_R) with the curvature coefficients h of CleanStagSph:
//!   Cartesian:       a = 1/dx, v = 1, s = hr = 0;
//!   spherical polar: a = CleanStagSph q = 0, v = 1, s its q = 4, hr half its q = 1;
//!   cubed sphere:    a = A_n (CleanCsNormA), v = V_f of the owners' half cells, s = hr = 0

struct CleanFaceMet {Real a, v, s, hr;};

KOKKOS_INLINE_FUNCTION
CleanFaceMet CleanFaceMetric(const CleanArgs &x, const int dir, const int m, const int k,
                             const int j, const int i) {
  const int kl = k - (dir == 2), jl = j - (dir == 1);
  if (x.cub) {
    return {CleanCsNormA(x, dir, m, k, j, i),
            0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,i - (dir == 0))), 0.0, 0.0};
  }
  if (!x.sph) return {x.idx[dir], 1.0, 0.0, 0.0};
  return {CleanStagSph(x, 0, dir, m, k, j, i), 1.0, CleanStagSph(x, 4, dir, m, k, j, i),
          0.5*CleanStagSph(x, 1, dir, m, k, j, i)};
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagMemUp
//! \brief clean_stag_gflux: the memory (1 - w) theta (G + d) of the energy flux through
//! face (k,j,i) normal to dir taken implicitly upwind, as q E_donor with q = memory/E_old
//! of the donor, |q| <= c (the donor-cell form of production's): on a
//! low-order face, and wherever its explicit outflow over the stage, h A/V |memory|,
//! would take more than CLN_MEMUP of the donor's E_old.  Returns q (> 0: from L, < 0:
//! from R), 0 where the memory stays explicit.  At E = E_old it is the explicit memory.

KOKKOS_INLINE_FUNCTION
Real CleanStagMemUp(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i) {
  if (!x.gf || x.hi) return 0.0;
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  // (clean_stag_vimp: G takes theta (1 - b), CleanFaceSource)
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0;
  const Real g = (1.0 - x.fc[dir](m,CLF_W,k,j,i))*x.fc[dir](m,CLF_TH,k,j,i)*
                 ((1.0 - qd)*x.fc[dir](m,CLF_G,k,j,i) +
                  (x.fd ? x.fc[dir](m,CLF_N+3,k,j,i) : 0.0));
  if (g == 0.0) return 0.0;
  const bool up = (g > 0.0);   // (the donor: L for an outflow to R)
  const Real eo = up ? x.w(m,CLW_EO,kl,jl,il) : x.w(m,CLW_EO,k,j,i);
  // (cubed sphere: A_f over the owner's V, so that a seam face decides alike on both
  // sides; a seam ghost's own V is this panel's extension, as CleanStagXSplit)
  const Real out = x.h*fabs(g)*(!x.cub ? (up ? CleanDivW(x, dir, m, kl, jl, il, 1) :
                                               CleanDivW(x, dir, m, k, j, i, 0)) :
                   ((dir == 0) ? x.ae0 : x.area[dir])(m,k,j,i)/
                   (up ? x.g(m,CLG_F+3,kl,jl,il) : x.g(m,CLG_F+3,k,j,i)));
  if (out <= CLN_MEMUP*eo) return 0.0;
  const Real q = (eo > 0.0) ? fmin(fabs(g)/eo, x.c) : x.c;
  return up ? q : -q;
}

//----------------------------------------------------------------------------------------
//! \fn int CleanCsSlot
//! \brief cubed sphere, staggered: the face (kf,jf,iff) normal to dir of cell (k,j,i), or
//! with far set the other face of its family, in the axes of the cell's owner: returns
//! 2 q + hi (q its family, hi = 1 for the owner's upper face), s = +-1 the owner's
//! orientation of the face normal along ours (a seam ghost's axes may be swapped and
//! flipped; the face is in the owner's x2 family if its x2 normal is along n)

KOKKOS_INLINE_FUNCTION
int CleanCsSlotQ(const DvceArray5D<Real> &g, const int dir, const int m, const int kf,
                 const int jf, const int iff, const int k, const int j, const int i,
                 Real &s) {
  Real t2 = 0.0, t3 = 0.0;
  for (int c = 0; c < 3; ++c) {
    t2 += g(m,CLG_F+c,k,j,i)*g(m,CLG_N+3*dir+c,kf,jf,iff);
    t3 += g(m,CLG_F+4+c,k,j,i)*g(m,CLG_N+3*dir+c,kf,jf,iff);
  }
  const int q = (fabs(t2) > fabs(t3)) ? 1 : 2;
  s = (((q == 1) ? t2 : t3) > 0.0) ? 1.0 : -1.0;
  return q;
}

KOKKOS_INLINE_FUNCTION
int CleanCsSlot(const CleanArgs &x, const int dir, const int m, const int kf,
                const int jf, const int iff, const int k, const int j, const int i,
                const bool far, Real &s) {
  s = 1.0;
  const bool up = (dir == 0) ? (iff > i) : ((dir == 1) ? (jf > j) : (kf > k));
  if (dir == 0) return (up != far) ? 1 : 0;
  // (CleanCsSlotQ of the face's cell on this side, CLG_Q)
  const int b = static_cast<int>(x.g(m,CLG_Q,kf,jf,iff)) >> (4*(dir - 1) + 2*up);
  const int q = 1 + (b & 1);
  s = (b & 2) ? 1.0 : -1.0;
  return 2*q + (((up == (s > 0.0)) != far) ? 1 : 0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagCsLat
//! \brief cubed sphere: n.(Y_p + Y_p') of cell (k,j,i) (CleanStagCsLateral), n the
//! normal of its face (kf,jf,iff) normal to dir and p, p' the cell's two face families
//! other than the face's, in the owner's axes (CleanCsSlot).  clean_cs_lat_wb, a lateral
//! face: n R in place of n, R the rotation of the cell's radial direction a (its centroid
//! direction, halo-exchanged) onto b, a projected into the face plane: a radial sum adds
//! n.b = 0, so a spherically symmetric P gives the face no force (the equiangular cells
//! are not mirror images, n.a_L + n.a_R != 0); O(dtheta^2) for a smooth P

KOKKOS_INLINE_FUNCTION
Real CleanStagCsLat(const CleanArgs &x, const int dir, const int m, const int kf,
                    const int jf, const int iff, const int k, const int j, const int i) {
  Real s, d = 0.0, n[3];   // (off the cubed sphere, clean_stag_lat_sums: q = dir, e_q)
  const int q = x.cub ? CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2 : dir;
  const int p = (q == 0) ? 1 : 0, pp = (q == 2) ? 1 : 2;
  for (int c = 0; c < 3; ++c) n[c] = x.cub ? x.g.At<CLG_N>(m,3*dir+c,kf,jf,iff) : (c == q);
  if (x.lwb && dir > 0) {   // (n R = n - w x n + w x (w x n)/(1 + a.b), w = a x b)
    Real a[3], b[3], aa = 0.0, an = 0.0, bb = 0.0, ab = 0.0;
    for (int c = 0; c < 3; ++c) {
      a[c] = x.g.At<CLG_X>(m,c,k,j,i);
      aa += a[c]*a[c];
      an += a[c]*n[c];
    }
    for (int c = 0; c < 3; ++c) {b[c] = a[c] - an*n[c]; bb += b[c]*b[c];}
    aa = 1.0/std::sqrt(aa);
    bb = 1.0/std::sqrt(bb);
    for (int c = 0; c < 3; ++c) {a[c] *= aa; b[c] *= bb; ab += a[c]*b[c];}
    const Real w[3] = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                       a[0]*b[1] - a[1]*b[0]};
    const Real u[3] = {w[1]*n[2] - w[2]*n[1], w[2]*n[0] - w[0]*n[2],
                       w[0]*n[1] - w[1]*n[0]};
    const Real v[3] = {w[1]*u[2] - w[2]*u[1], w[2]*u[0] - w[0]*u[2],
                       w[0]*u[1] - w[1]*u[0]};
    for (int c = 0; c < 3; ++c) n[c] += v[c]/(1.0 + ab) - u[c];
  }
  for (int c = 0; c < 3; ++c) {
    d += n[c]*(x.cub ? x.w(m,CLW_L+3*p+c,k,j,i) + x.w(m,CLW_L+3*pp+c,k,j,i)
                     : x.yl(m,3*p+c,k,j,i) + x.yl(m,3*pp+c,k,j,i));
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagN
//! \brief the component along face (kf,jf,iff) normal to dir of the vector in work slots
//! v, v+1, v+2 of cell (k,j,i): slot v + dir; on the cubed sphere (Cartesian) n.v

KOKKOS_INLINE_FUNCTION
Real CleanStagN(const CleanArgs &x, const int v, const int dir, const int m,
                const int kf, const int jf, const int iff, const int k, const int j,
                const int i) {
  if (!x.cub) return x.w(m,v+dir,k,j,i);
  Real d = 0.0;
  for (int b = 0; b < 3; ++b) d += x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff)*x.w(m,v+b,k,j,i);
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagArea
//! \brief the net part s = h+ + h- of the streaming share w c (h+ E_L + h- E_R) through
//! face (k,j,i) normal to dir in area form, from the donor upwind of s: adds w c s
//! (A_d/A_f - 1) to the coefficient pl or pr of E_L or E_R, A_d the area at the donor's
//! centre (spherical polar: x1v^2 of r faces, sin(x2v) of theta faces; cubed sphere: the
//! mean area of its two faces of the family, the owner's, CLG_A), so that a free stream
//! keeps A_d E_d (an isotropic field, s = 0, is unchanged; the signs, and with them the
//! M-matrix, are kept).  Only where h+- are the cells' own, analytic or of point
//! intensities (clean_rays = sc): the flux factors of clean_stag_gflux are half-range
//! fluxes through the face, already in area form.  CleanAreaRatio: A_d/A_f of donor
//! (kd,jd,id) and face (k,j,i).

KOKKOS_INLINE_FUNCTION
Real CleanAreaRatio(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i, const int kd, const int jd, const int id) {
  if (x.cub) {
    Real sd;
    return x.g.At<CLG_A>(m,CleanCsSlot(x, dir, m, k, j, i, kd, jd, id, false, sd)/2,
                         kd,jd,id)/x.area[dir](m,k,j,i);
  }
  return (dir == 0) ? x.x1v(m,id)*x.x1v(m,id)/(x.x1f(m,i)*x.x1f(m,i))
                    : x.tg(m,2,jd)/x.tg(m,3,j);
}

KOKKOS_INLINE_FUNCTION
void CleanStagArea(const CleanArgs &x, const int dir, const int m, const int k,
                   const int j, const int i, const Real w, Real &pl, Real &pr) {
  if (!x.sar || x.hi || !(x.sph || x.cub) || (!x.cub && dir == 2)) return;
  const Real s = w*x.c*(x.fc[dir](m,CLF_HP,k,j,i) + x.fc[dir](m,CLF_HM,k,j,i));
  const int d = (s > 0.0) ? 1 : 0, id = i - d*(dir == 0), jd = j - d*(dir == 1);
  const int kd = k - d*(dir == 2);
  ((d == 1) ? pl : pr) += s*(CleanAreaRatio(x, dir, m, k, j, i, kd, jd, id) - 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceChi
//! \brief chi of face (k,j,i) normal to dir: the mean of its two cells; with
//! clean_face_opac = harmonic an x1 face takes their distance-weighted harmonic mean (old
//! production implicit_face_opac = harmonic): at an unresolved opacity jump the
//! transparent side's mean free path limits the flux.  A MeshBlock face takes the
//! neighbour's ghost like an interior face; a physical boundary face keeps the mean.

KOKKOS_INLINE_FUNCTION
Real CleanFaceChi(const DvceArray5D<Real> &op, const DvceArray2D<Real> &x1f,
                  const DvceArray2D<BoundaryFlag> &bc, const bool hrm, const int dir,
                  const int m, const int k, const int j, const int i, const int is,
                  const int ie) {
  const Real cl = op(m,M1_OP_T,k-(dir == 2),j-(dir == 1),i-(dir == 0));
  const Real cr = op(m,M1_OP_T,k,j,i);
  auto phys = [&](const int s) {
    return bc(m,s) != BoundaryFlag::block && bc(m,s) != BoundaryFlag::periodic;
  };
  if (!hrm || dir != 0 || (i == is && phys(0)) || (i == ie + 1 && phys(1)) ||
      !(cl > 0.0 && cr > 0.0)) {
    return 0.5*(cl + cr);
  }
  if (x1f.extent_int(1) <= i + 1) return 2.0*cl*cr/(cl + cr);   // (Cartesian: uniform)
  const Real hl = x1f(m,i) - x1f(m,i-1), hr = x1f(m,i+1) - x1f(m,i);
  return (hl + hr)/(hl/cl + hr/cr);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceWeights
//! \brief coefficients of E_L and E_R in the flux through face (k,j,i) normal to dir,
//! i.e. everything but the off-diagonal Eddington terms and the old-vector flux G.
//! Physical boundaries (non-periodic) take the boundary flux of their form (none at a
//! wall or an imposed flux, CleanFaceSource).  up (if given) is the donor-cell share phi
//! of clean_stag_upwind.

KOKKOS_INLINE_FUNCTION
bool CleanBoundaryFace(const CleanArgs &x, const int dir, const int m, const int k,
                       const int j, const int i);

KOKKOS_INLINE_FUNCTION
void CleanFaceWeights(const CleanArgs &x, const int dir, const int m, const int k,
                      const int j, const int i, Real &al, Real &ar, Real *up = nullptr) {
  if (x.fwc) {
    al = x.fc[dir](m,CLF_QL,k,j,i);
    ar = x.fc[dir](m,CLF_QL+1,k,j,i);
    return;
  }
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  al = 0.0;
  ar = 0.0;
  // (clean_bc_area: an open r end's outgoing c h+- E in area form as CleanStagArea, the
  // end cell the donor: its net stream s = h+ + h- out of the domain takes A_d/A_f, the
  // part the bath's isotropic inflow balances (and a net inflow, the bath donor) none)
  const bool ba = x.bar && !x.marshak && dir == 0 && (x.sph || x.cub);
  if (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC) {
    if (!CleanOpen(x, m, 2*dir)) return;   // (a wall; an imposed flux is all source)
    const bool beam = x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0;
    if (beam && x.marshak) return;
    ar = x.marshak ? -x.c*x.marshak_q : x.c*x.cs(m,CLN_HM+dir,k,j,i);
    if (ba) ar += x.c*fmin(x.cs(m,CLN_HP,k,j,i) + x.cs(m,CLN_HM,k,j,i), 0.0)*
                  (CleanAreaRatio(x, 0, m, k, j, i, k, j, i) - 1.0);
    // (clean_stag_upwind: an open end carries the advective flux a E of the energy, the
    // end cell's when it flows out, the bath's when in, CleanFaceSource)
    if (x.upw && x.feedback && !x.hi && CleanOpen(x, m, 2*dir) && !beam) {
      ar += fmin(CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i), 0.0);
    }
    return;
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC) {
    if (!CleanOpen(x, m, 2*dir+1)) return;
    al = x.marshak ? x.c*x.marshak_q : x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
    if (ba) al += x.c*fmax(x.cs(m,CLN_HP,kl,jl,il) + x.cs(m,CLN_HM,kl,jl,il), 0.0)*
                  (CleanAreaRatio(x, 0, m, k, j, i, kl, jl, il) - 1.0);
    if (x.upw && x.feedback && !x.hi && CleanOpen(x, m, 2*dir+1)) {
      al += fmax(CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il), 0.0);
    }
    return;
  }
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  // (clean_stag_vimp: the stress takes theta_d = theta + b (1 - theta), the advective
  // flux 1 - theta_d and the beta term theta (1 - b), CleanStagDrag)
  const Real t0 = x.fc[dir](m,CLF_TH,k,j,i);
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0, th = t0 + qd*(1.0 - t0);
  // (the face-normal momentum's stress a n.D.n E/v, curved: and its curvature terms)
  const CleanFaceMet fm = CleanFaceMetric(x, dir, m, k, j, i);
  const Real dl = CleanStagDnn(x, dir, m, k, j, i, kl, jl, il);
  const Real dr = CleanStagDnn(x, dir, m, k, j, i, k, j, i);
  const Real gs = (1.0 - w)*th*x.h*x.c*x.c;
  if (x.sph) {
    al = gs*(fm.a*dl - fm.hr*CleanStagSph(x, 2, dir, m, k, j, i));
    ar = -gs*(fm.a*dr + fm.hr*CleanStagSph(x, 3, dir, m, k, j, i));
  } else {
    al = gs*fm.a/fm.v*dl;
    ar = -(gs*fm.a/fm.v)*dr;
  }
  al += w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
  ar += w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  CleanStagArea(x, dir, m, k, j, i, w, al, ar);
  const Real q = CleanStagMemUp(x, dir, m, k, j, i);
  if (q > 0.0) al += q;
  if (q < 0.0) ar += q;
  if (x.feedback) {
    const Real s = 0.5*(1.0 - w)*(1.0 - th);
    const Real ul = CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il);
    const Real ur = CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i);
    const Real sl = s*ul, sr = s*ur;
    // the beta term th c beta_f (hc kP a T^4 - hc kE E)_f of eq. 57 (Ma et al. 2025)
    const Real b = -0.25*(1.0 - w)*t0*(1.0 - qd)*x.c*
                   (CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                    CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*x.h*x.c;
    // of (1 - r) kE E: the emission follows E through the gas T (r = s/(d cv + s)), as
    // the old production's linearised exchange g0 = -(SRCR - SRCB E), SRCB = (1 - r) hc
    // kE; the whole hc kE beta made the face an advection at ~ tau_cell v, which broke
    // the M-matrix of an opaque moving cell at any dt
    const Real rl = 1.0 - x.w(m,CLW_GX,kl,jl,il);
    const Real rr = 1.0 - x.w(m,CLW_GX,k,j,i);
    const Real bl = b*x.op(m,M1_OP_E,kl,jl,il)*rl, br = b*x.op(m,M1_OP_E,k,j,i)*rr;
    al += sl;
    ar += sr;
    al += bl;
    ar += br;
    // (clean_stag_upwind: where the central flux breaks the M-matrix, al < 0 or ar > 0
    // (cell Peclet number > 1), the share phi of the advective flux that goes donor cell
    // is the least that restores it, continuous in the coefficients, so that the operator
    // changes smoothly between passes; CleanStagAdvectCorr adds phi times the rest,
    // lagged)
    if (x.upw && !(al >= 0.0 && ar <= 0.0)) {
      auto need = [](const Real c, const Real d) {   // least p >= 0 with c + p d >= 0
        return (c >= 0.0) ? 0.0 : ((d > 0.0) ? -c/d : 1.0);
      };
      const Real o = (ul + ur > 0.0) ? 1.0 : -1.0;
      const Real phi = fmin(fmax(need(al, o*sl), need(-ar, o*sr)), 1.0);
      if (up) *up = phi;
      al += o*phi*sl;
      ar -= o*phi*sr;
    }
  }
  // (clean_stag_kpos = delta: the slopes' linear flux tau (al E_L + ar E_R) implicit on
  // an interior face, the E rows only (the momentum's args carry kq = -1))
  if (x.kq >= 0 && x.kt != 0.0 && !CleanBoundaryFace(x, dir, m, k, j, i)) {
    al *= 1.0 + x.kt;
    ar *= 1.0 + x.kt;
  }
  // (clean_stag_xdonor: the cross flux q E_donor, CleanStagXDonor; keeps al >= 0 >= ar)
  if (x.xq >= 0 && !x.hi) {   // (the energy flux; F^H keeps the central form)
    al += x.fc[dir](m,x.xq,k,j,i);
    ar += x.fc[dir](m,x.xq+1,k,j,i);
  }
  // (clean_stag_kpos: the increment's moved flux q E_donor, CleanStagKPos; the E rows
  // and the write-back's energy only: the momentum's args carry kq = -1)
  if (x.kq >= 0) {
    al += x.fc[dir](m,x.kq,k,j,i);
    ar += x.fc[dir](m,x.kq+1,k,j,i);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanKposFlux
//! \brief clean_stag_kpos: the moved increment flux q E_donor of face (k,j,i) normal to
//! dir in work vector e (0 without kpos); the stored Phi and the slopes leave it out

KOKKOS_INLINE_FUNCTION
Real CleanKposFlux(const CleanArgs &x, const int e, const int dir, const int m,
                   const int k, const int j, const int i) {
  if (x.kq < 0) return 0.0;
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real f = x.fc[dir](m,x.kq,k,j,i)*x.w(m,e,kl,jl,il) +
           x.fc[dir](m,x.kq+1,k,j,i)*x.w(m,e,k,j,i);
  // (delta: tau al/ar part)
  if (x.kt != 0.0 && !CleanBoundaryFace(x, dir, m, k, j, i)) {
    CleanArgs y = x;
    y.kq = -1;
    y.xq = -1;
    Real al, ar;
    CleanFaceWeights(y, dir, m, k, j, i, al, ar);
    f += x.kt*(al*x.w(m,e,kl,jl,il) + ar*x.w(m,e,k,j,i));
  }
  return f;
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanBoundaryFace
//! \brief true for a face on a physical (non-periodic) boundary

KOKKOS_INLINE_FUNCTION
bool CleanBoundaryFace(const CleanArgs &x, const int dir, const int m, const int k,
                       const int j, const int i) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  return (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC) ||
         (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC);
}

//----------------------------------------------------------------------------------------
//! \fn int CleanTransSide
//! \brief 1 if the transverse stress terms of cell (k,j,i) reach across its face on side
//! s (0 low, 1 high) along b: inside, across a MeshBlock face or a periodic face and a
//! reflecting wall (ghosts that carry the field: centred through them), 0 at an open
//! boundary, an imposed flux, a pole and a cubed-sphere seam, whose ghosts do not (the
//! pole's are a mirror at the same phi, not the cell across the axis; a seam ghost is
//! the neighbour panel's cell, off the line, as production's in-panel M1CsTanDer; one
//! sided, the cell itself in place of the ghost: production's clamped form, M1OffDiv)

KOKKOS_INLINE_FUNCTION
int CleanTransSide(const CleanArgs &x, const int b, const int m, const int k, const int j,
                   const int i, const int s) {
  const int pos = (b == 0) ? i : ((b == 1) ? j : k);
  const int end = s ? ((b == 0) ? x.ie : ((b == 1) ? x.je : x.ke))
                    : ((b == 0) ? x.is : ((b == 1) ? x.js : x.ks));
  const int t = CleanBc(x, m, 2*b + s);
  const BoundaryFlag f = x.mbc(m,2*b+s);
  return (pos != end || (f != BoundaryFlag::polar && f != BoundaryFlag::panel &&
                         (t == M1_IBC_PERIODIC || t == M1_IBC_REFLECT))) ? 1 : 0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCsSeamCoef
//! \brief cubed sphere: the ghost's P - P0 of the cell (k,j,i) next to a seam along q as
//! c0 P0 + c1 P1 + c2 P2 over its cells 0, si, 2 si in panel, P = D_n E: the in-panel
//! second-order extrapolant (d = (4 P1 - 3 P0 - P2)/2, the ghost P1 - 2 d), with d
//! clamped to [0, 2 s1] (s1 = P1 - P0) as production's M1CsTanDer with cs_cross_mono
//! (M1CsTanMono, k = 0).  The clamp's case is taken from the old vector E^n (CLW_EO), so
//! that the operator is one linear map for the whole step (its stencil rows, cln_csyb,
//! the same), continuous in E^n: at either bound the two cases coincide

KOKKOS_INLINE_FUNCTION
void CleanCsSeamCoef(const CleanArgs &x, const int n, const int q, const int m,
                     const int k, const int j, const int i, const int si, Real c[3]) {
  const int dk = (q == 2)*si, dj = (q == 1)*si, di = (q == 0)*si;
  auto p = [&](const int o) {
    return x.cs(m,n,k+o*dk,j+o*dj,i+o*di)*x.w(m,CLW_EO,k+o*dk,j+o*dj,i+o*di);
  };
  const Real s1 = p(1) - p(0), d = 0.5*(4.0*p(1) - 3.0*p(0) - p(2));
  const Real lo = fmin(0.0, 2.0*s1), hi = fmax(0.0, 2.0*s1);
  const bool cl = (d < lo || d > hi), z = (s1 > 0.0) ? (d < lo) : (d > hi);
  c[0] = !cl ? 2.0 : (z ? -1.0 : 3.0);   // (d: P1 - 2 d; d = 0: P1; d = 2 s1: P1 - 4 s1)
  c[1] = !cl ? -3.0 : (z ? 1.0 : -3.0);
  c[2] = !cl ? 1.0 : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceStress
//! \brief divergence of the radiation stress D E at face (k,j,i) normal to dir,
//! d_b(D_{dir,b} E): the normal derivative from the two adjacent cells (if normal is
//! set) and, away from physical boundaries, the transverse derivatives averaged over the
//! four cells around the face.  E is work vector e.

CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceStress(const CleanArgs &x, const int e, const int dir, const int m,
                     const int k, const int j, const int i, const bool normal) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real div = 0.0;
  const CleanFaceMet fm = CleanFaceMetric(x, dir, m, k, j, i);
  if (normal) {   // (the normal stress of the two (half) cells, CleanFaceMet)
    div = fm.a*(CleanStagDnn(x, dir, m, k, j, i, k, j, i)*x.w(m,e,k,j,i) -
                CleanStagDnn(x, dir, m, k, j, i, kl, jl, il)*x.w(m,e,kl,jl,il));
    if (x.sph) {
      div += fm.hr*(CleanStagSph(x, 2, dir, m, k, j, i)*x.w(m,e,kl,jl,il) +
                    CleanStagSph(x, 3, dir, m, k, j, i)*x.w(m,e,k,j,i));
    }
  }
  if (x.cub) {   // (cubed sphere: the lateral sums CleanStagCsLat, all over V_f; a seam
    // face takes none: its mirror chord is normal to the face, so the two-point term is
    // the whole n.div(D E), as production's seam row)
    if (!CleanBoundaryFace(x, dir, m, k, j, i) && !CleanCsSeamFace(x, dir, m, k, j, i)) {
      div += CleanStagCsLat(x, dir, m, k, j, i, kl, jl, il) +
             CleanStagCsLat(x, dir, m, k, j, i, k, j, i);
    }
    return div/fm.v;
  }
  if (!x.multi_d || CleanBoundaryFace(x, dir, m, k, j, i)) return div;
  // (clean_stag_lat_sums: the lateral sums Y of the two cells over V_f, as the cubed
  // sphere's; spherical polar phi faces keep the torque form below)
  if (x.ys && !(x.sph && dir == 2)) {
    const Real vf = x.sph ? 0.5*(x.vol(m,k,j,i) + x.vol(m,kl,jl,il)) : 1.0;
    return div + (CleanStagCsLat(x, dir, m, k, j, i, kl, jl, il) +
                  CleanStagCsLat(x, dir, m, k, j, i, k, j, i))/vf;
  }
  for (int b = 0; b < 3; ++b) {
    if (b == dir || (b == 1 && !x.multi_d) || (b == 2 && !x.three_d)) continue;
    const int ab = (dir + b == 1) ? 3 : ((dir + b == 2) ? 4 : 5);   // D12, D13, D23
    const int dk = (b == 2), dj = (b == 1), di = (b == 0);
    auto de = [&](const int kk, const int jj, const int ii) {
      return x.cs(m,ab,kk,jj,ii)*x.w(m,e,kk,jj,ii);
    };
    // (spherical polar: d_r on theta faces, (1/r) d_theta on r faces, central; one sided
    // where a side's ghosts do not carry the field, CleanTransSide)
    const int up = CleanTransSide(x, b, m, k, j, i, 1);
    const int dn = CleanTransSide(x, b, m, k, j, i, 0);
    // (torque form: dual r faces i, i+1 and theta faces j, j+1)
    if (x.sph && dir == 2) {
      const Real vr = x.vol(m,k,j,i)*x.x1v(m,i)*x.tg(m,2,j);
      auto pe = [&](const int jj, const int ii) {   // (D E)_{b,ph} at a dual face
        const int jh = (jj > j && !up) ? j : jj, ih = (ii > i && !up) ? i : ii;
        const int jo = jj - (b == 1), io = ii - (b == 0);
        const int jw = (jo < j && !dn) ? j : jo, iw = (io < i && !dn) ? i : io;
        return 0.25*(de(k,jh,ih) + de(k-1,jh,ih) + de(k,jw,iw) + de(k-1,jw,iw));
      };
      if (b == 0) {   // (clean_lz_bnd_free: no torque through a radial boundary)
        auto tq = [&](const int ii) {
          return (x.lzbf && CleanBoundaryFace(x, 0, m, k, j, ii)) ? 0.0 :
                 x.area[0](m,k,j,ii)*x.x1f(m,ii)*pe(j,ii);
        };
        div += (tq(i+1) - tq(i))*x.tg(m,2,j)/vr;
      } else {
        div += (x.area[1](m,k,j+1,i)*x.tg(m,3,j+1)*pe(j+1,i) -
                x.area[1](m,k,j,i)*x.tg(m,3,j)*pe(j,i))*x.x1v(m,i)/vr;
      }
      continue;
    }
    if (up + dn == 0) continue;
    const Real rb = !x.sph ? 1.0 : ((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
    const Real w = (up + dn == 2) ? 1.0 : 2.0;   // (the span: two cells or one)
    const Real ib = !x.sph ? w*x.idx[b] :
                    ((b == 0) ? 2.0/(x.x1v(m,i+up) - x.x1v(m,i-dn)) :
                    ((b == 1) ? 2.0/(rb*(x.x2v(m,j+up) - x.x2v(m,j-dn))) :
                    w*x.idx[2]/(rb*((dir == 1) ? x.tg(m,3,j) : x.tg(m,2,j)))));
    div += 0.25*ib*(de(k+up*dk,j+up*dj,i+up*di) + de(kl+up*dk,jl+up*dj,il+up*di)
                    - de(k-dn*dk,j-dn*dj,i-dn*di) - de(kl-dn*dk,jl-dn*dj,il-dn*di));
  }
  return div;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceCross
//! \brief flux of the off-diagonal Eddington terms, -(1-w) theta h c^2 d_b(D_{dir,b} E)
//! for b != dir.  E is work vector e.

KOKKOS_INLINE_FUNCTION
Real CleanFaceCross(const CleanArgs &x, const int e, const int dir, const int m,
                    const int k, const int j, const int i) {
  if (!x.multi_d || (x.xq >= 0 && !x.hi) || CleanBoundaryFace(x, dir, m, k, j, i)) {
    return 0.0;
  }
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  const Real t0 = x.fc[dir](m,CLF_TH,k,j,i);   // (clean_stag_vimp: theta_d, as
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0;   //  CleanFaceWeights)
  const Real th = t0 + qd*(1.0 - t0);
  return (x.fwc ? x.fc[dir](m,CLF_QL+2,k,j,i) : -(1.0 - w)*th*x.h*x.c*x.c)*
         CleanFaceStress(x, e, dir, m, k, j, i, false);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagAdvectCorr
//! \brief clean_stag_upwind: on a face whose central advective flux breaks the
//! M-matrix (cell Peclet number > 1, CleanFaceWeights) the share phi of the flux of face
//! (k,j,i) normal to dir is
//!   A_f = a_up E_up + [a_f E_f - a_up E_up]^k,  a_f = (a_L + a_R)/2,
//! the donor cell (upwind on a_f) in the matrix, which keeps the M-matrix at any
//! beta tau_cell, and the rest at the iterate E^k (CLW_EK, a deferred correction): E_f is
//! the van Leer value from the upwind side (src/reconstruct/plm.hpp).  At the fixed point
//! the face carries a_f E_f, second order in smooth flows and monotone at jumps; the
//! rest (all of it where the M-matrix holds) stays central.
//! clean_stag_upwind_central: the face carries a_f E_up + a_down m/2, m = minmod(2
//! (E_up - E_upup), E_down - E_up): the central mean of the cells' a E, (a_L E_L + a_R
//! E_R)/2 (the old production's cell flux a E), wherever (E_up - E_upup)/(E_down - E_up)
//! >= 1/2 (Sweby's TVD bound), and a_f E_up at an extremum of E.  The van Leer a_f E_f
//! is not central at r != 1, an upwind-biased O(dx^2) term whose up- and downflow parts
//! do not cancel in a convective zone: +10 % of the net advective flux at cell Peclet >
//! 1, the RSG/HG deep L +1.5-3.8 % of old's at 0.2-0.75 R_ph (both the energy flux and
//! F^H)

KOKKOS_INLINE_FUNCTION
Real CleanStagAdvectCorr(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  auto e = [&](const int o) {return x.w(m,CLW_EK,k+o*dk,j+o*dj,i+o*di);};
  const Real ul = CleanStagN(x, CLW_A1, dir, m, k, j, i, k-dk, j-dj, i-di);
  const Real ur = CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i);
  if (x.ucen) {   // (sl: the limited E_down - E_up; up, the upwind cell)
    const bool fw = (ul + ur > 0.0);
    const Real eu = fw ? e(-1) : e(0), df = (fw ? e(0) : e(-1)) - eu;
    const Real du = 2.0*(eu - (fw ? e(-2) : e(1)));
    const Real sl = (df*du > 0.0) ? ((fabs(df) < fabs(du)) ? df : du) : 0.0;
    return 0.5*((fw ? ur - ul : ul - ur)*eu + (fw ? ur : ul)*sl);
  }
  Real ef, dum;
  if (ul + ur > 0.0) {
    PLM(e(-2), e(-1), e(0), ef, dum);
    return 0.5*(ul + ur)*ef - ul*e(-1);
  }
  PLM(e(-1), e(0), e(1), dum, ef);
  return 0.5*(ul + ur)*ef - ur*e(0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanThinFaceDG
//! \brief staggered: the lagged additive flux DG of a face at thin weight w, the old
//! production's M1ThinFace (implicit_thin_corr = hc): w times CleanStagBracket's
//! q c (F_y - h+ J_L - h- J_R) (0 without the bracket)

KOKKOS_INLINE_FUNCTION
Real CleanThinFaceDG(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i, const Real w) {
  return (x.fdg >= 0) ? w*x.fc[dir](m,x.fdg,k,j,i) : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceSource
//! \brief the part of the face flux that does not depend on E: (1-w) theta G inside, the
//! incoming flux (beam, bath) or the imposed flux on physical boundaries

CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i, const Real *phw = nullptr) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  if (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir) == M1_IBC_REFLECT) return 0.0;
    if (CleanBc(x, m, 2*dir) == M1_IBC_FLUX) return dir ? x.ifx2[0](m,k,i) : x.iflux[0];
    if (x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0) {
      return x.beam(m,k-x.ks,j-x.js);
    }
    const bool up = x.upw && x.feedback && !x.hi;
    const Real ab = up ? fmax(CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i), 0.0) : 0.0;
    return (x.marshak ? x.c*x.marshak_q*x.ebath[2*dir] : 0.25*x.c*x.ebath[2*dir]) +
           ab*x.ebath[2*dir];
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir+1) == M1_IBC_REFLECT) return 0.0;
    if (CleanBc(x, m, 2*dir+1) == M1_IBC_FLUX) return dir ? x.ifx2[1](m,k,i) : x.iflux[1];
    const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
    const bool up = x.upw && x.feedback && !x.hi;
    const Real ab = up ? fmin(CleanStagN(x, CLW_A1, dir, m, k, j, i, kl, jl, il), 0.0)
                       : 0.0;
    return (x.marshak ? -x.c*x.marshak_q*x.ebath[2*dir+1] : -0.25*x.c*x.ebath[2*dir+1]) +
           ab*x.ebath[2*dir+1];
  }
  // (clean_stag_xdonor: the explicit share of the lagged fluxes of the pass,
  // CleanStagXSplit)
  if (x.xq >= 0 && !x.hi) return x.fc[dir](m,x.xq+2,k,j,i);
  const Real w = x.hi ? 0.0 : x.fc[dir](m,CLF_W,k,j,i);
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const Real th = (1.0 - w)*x.fc[dir](m,CLF_TH,k,j,i);
  // (clean_stag_vimp: G and the beta term take theta (1 - b), and b Q joins them)
  const Real qd = x.drag ? x.fc[dir](m,CLF_DB,k,j,i) : 0.0, tg = th*(1.0 - qd);
  Real s = tg*x.fc[dir](m,CLF_G,k,j,i);
  if (x.fd && !x.hi) s += th*x.fc[dir](m,CLF_N+3,k,j,i);   // (fkin: of F^H, the last F^H)
  if (CleanStagMemUp(x, dir, m, k, j, i) != 0.0) s = 0.0;   // (gflux: upwind, implicit)
  if (x.feedback) {   // (the emission part of the beta term, CleanFaceWeights)
    s += 0.25*tg*x.c*(CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                      CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*
         (x.w(m,CLW_RE,kl,jl,il) + x.w(m,CLW_RE,k,j,i));
    // (less the emission's linear response at E_k, CleanFaceWeights)
    {
      auto rk = [&](const int kk, const int jj, const int ii) {
        return x.w(m,CLW_GX,kk,jj,ii)*x.op(m,M1_OP_E,kk,jj,ii)*x.w(m,CLW_EK,kk,jj,ii);
      };
      s -= 0.25*tg*x.c*(CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                        CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*x.h*x.c*
           (rk(kl, jl, il) + rk(k, j, i));
    }
  }
  if (x.drag) s += (1.0 - w)*qd*x.fc[dir](m,CLF_DQ,k,j,i);
  if (!x.hi) s += CleanThinFaceDG(x, dir, m, k, j, i, w);   // (CleanStagBracket's DG)
  // (the advective share 1 - theta_d of an upwind face, CleanFaceWeights)
  if (x.upw && x.feedback) {
    Real al, ar, phi = 0.0;   // (phw: the caller's phi of these args, CleanFaceFlux)
    if (phw) phi = *phw; else CleanFaceWeights(x, dir, m, k, j, i, al, ar, &phi);
    if (phi > 0.0) {
      s += phi*((1.0 - w) - th)*(1.0 - qd)*CleanStagAdvectCorr(x, dir, m, k, j, i);
    }
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagXSplit
//! \brief clean_stag_xdonor: the lagged fluxes of face (k,j,i) normal to dir for the
//! pass, the cross flux C(E_k) (CleanStagXDonor, slot xq + 3) and the face source S (xr:
//! as it is), into a donor-cell part and an explicit part:
//!   p1:  C donor-implicit; with xsrc S as lag below, else S explicit;
//!   lag: L = C (+ S with xsrc) explicit up to its share s of the donor, where its
//!        outflow over the stage h (A/V) |L| exceeds CLN_MEMUP E_old,donor (>= 0): the
//!        explicit part s L, s = CLN_MEMUP E_old/outflow (continuous), the rest donor.
//! A donor part l enters as q E_donor, q = |l|/E_k,donor at most c (slots xq: qL >= 0
//! of E_L, xq + 1: qR <= 0 of E_R), the explicit part in slot xq + 2
//! (CleanFaceSource).  So the right-hand side keeps E_old - sum of the explicit
//! outflows >= 0 and the operator its M-matrix; at the fixed point the face carries
//! C + S wherever q < c.

KOKKOS_INLINE_FUNCTION
void CleanStagXSplit(const CleanArgs &x, const CleanArgs &xr, const int dir, const int m,
                     const int k, const int j, const int i) {
  const auto &f = x.fc[dir];
  const Real sr = CleanFaceSource(xr, dir, m, k, j, i);
  Real ql = 0.0, qr = 0.0, ex = 0.0;
  if (CleanBoundaryFace(x, dir, m, k, j, i)) {
    ex = sr;
  } else {
    const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
    // (the explicit share s of a lagged flux l: its outflow over the stage at most
    // CLN_MEMUP of the donor's E_old)
    auto share = [&](const Real l) {
      if (l == 0.0) return 1.0;
      const bool up = (l > 0.0);
      const Real eo = fmax(up ? x.w(m,CLW_EO,kl,jl,il) : x.w(m,CLW_EO,k,j,i), 0.0);
      // (cubed sphere: A_f over the owner's V, so a seam face takes one share on both
      // sides; the ghost's own V is this panel's extension, not the neighbour's cell)
      const Real out = x.h*fabs(l)*(!x.cub ? (up ? CleanDivW(x, dir, m, kl, jl, il, 1) :
                                                   CleanDivW(x, dir, m, k, j, i, 0)) :
                     ((dir == 0) ? x.ae0 : x.area[dir])(m,k,j,i)/
                     (up ? x.g(m,CLG_F+3,kl,jl,il) : x.g(m,CLG_F+3,k,j,i)));
      return (out > CLN_MEMUP*eo) ? CLN_MEMUP*eo/out : 1.0;
    };
    // (a donor part li: q = |li|/E_donor at most c; the cross flux of E_k over E_k (at
    // the fixed point exactly C), a face source, E-independent, over E_old as
    // CleanStagMemUp (the donor-cell form in time, one linear map for the step))
    auto put = [&](const Real li, const bool src) {
      if (li == 0.0) return;
      const bool up = (li > 0.0);
      const int e = src ? CLW_EO : CLW_EK;
      const Real ed = up ? x.w(m,e,kl,jl,il) : x.w(m,e,k,j,i);
      const Real q = (fabs(li) < x.c*ed) ? fabs(li)/ed : x.c;
      if (up) ql += q;
      if (!up) qr -= q;
    };
    const Real cf = f(m,x.xq+3,k,j,i);
    const Real sx = x.xsrc ? sr : 0.0;   // (the face source in the split)
    if (!x.xsrc) ex += sr;
    if (x.xlag) {
      const Real s = share(cf + sx);
      ex += s*(cf + sx);
      put((1.0 - s)*cf, false);
      put((1.0 - s)*sx, true);
    } else {
      put(cf, false);
      const Real s = share(sx);
      ex += s*sx;
      put((1.0 - s)*sx, true);
    }
  }
  f(m,x.xq,k,j,i) = ql;
  f(m,x.xq+1,k,j,i) = qr;
  f(m,x.xq+2,k,j,i) = ex;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFlux
//! \brief flux through face (k,j,i) normal to dir for E in work vector e; with_source
//! adds the E-independent part

CLEAN_DEVICE_NOINLINE CLEAN_HIP_INLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceFlux(const CleanArgs &x, const int e, const int dir, const int m,
                   const int k, const int j, const int i, const bool with_source) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real al, ar, phi = 0.0;   // (phi: the source's upwind share, weights once)
  CleanFaceWeights(x, dir, m, k, j, i, al, ar, &phi);
  Real f = CleanFaceCross(x, e, dir, m, k, j, i);
  if (al != 0.0) f += al*x.w(m,e,kl,jl,il);
  if (ar != 0.0) f += ar*x.w(m,e,k,j,i);
  if (with_source) f += CleanFaceSource(x, dir, m, k, j, i, &phi);
  return f;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFluxCall, CleanFaceStressCall
//! \brief the face flux and stress as calls on HIP (CLEAN_HIP_INLINE), for a path the
//! pass does not take

CLEAN_HIP_NOINLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceFluxCall(const CleanArgs &x, const int e, const int dir, const int m,
                       const int k, const int j, const int i, const bool with_source) {
  return CleanFaceFlux(x, e, dir, m, k, j, i, with_source);
}

CLEAN_HIP_NOINLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceStressCall(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i, const bool hi) {
  return CleanFaceStress(x, e, dir, m, k, j, i, hi);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFluxDivergence
//! \brief h div F of cell (k,j,i) for E in work vector e

KOKKOS_INLINE_FUNCTION
Real CleanFluxDivergence(const CleanArgs &x, const int e, const int m, const int k,
                         const int j, const int i, const bool with_source) {
  Real d = 0.0;
  for (int n = 0; n < 3; ++n) {   // (constant trip: x.fc[n], x.area[n] resolved)
    if (n >= x.nf) break;
    const int dk = (n == 2), dj = (n == 1), di = (n == 0);
    d += CleanDivF(x, n, m, k, j, i, CleanFaceFlux(x, e, n, m, k, j, i, with_source),
                   CleanFaceFlux(x, e, n, m, k+dk, j+dj, i+di, with_source));
  }
  return x.h*d;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStAdd
//! \brief adds f D_ab E of the cell at offset (OK, OJ, OI) from (k,j,i) to its slot of
//! row c

template <int OK, int OJ, int OI>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStAdd(const CleanArgs &x, const int ab, const int m, const int k, const int j,
                const int i, const Real f, Real (&c)[CLN_ST]) {
  c[CleanStSlot(OK, OJ, OI)] += f*x.cs(m,ab,k+OK,j+OJ,i+OI);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStCross
//! \brief the row terms of the transverse derivative along B of the off-diagonal stress
//! through face SD (0 low, 1 high) normal to D of cell (k,j,i), times cr: the terms of
//! CleanFaceStress (normal = false) with every slot a constant (as production's
//! ImplicitStencilBuild)

template <int D, int SD, int B>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStCross(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  const Real cr, Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  constexpr int bk = (B == 2), bj = (B == 1), bi = (B == 0);
  constexpr int ab = (D + B == 1) ? 3 : ((D + B == 2) ? 4 : 5);   // D12, D13, D23
  // (the cells L, R of the face)
  constexpr int lk = (SD - 1)*ek, lj = (SD - 1)*ej, li = (SD - 1)*ei;
  constexpr int rk = SD*ek, rj = SD*ej, ri = SD*ei;
  const int kf = k + rk, jf = j + rj, iff = i + ri;
  // (one sided where a side's ghosts do not carry the field, CleanTransSide)
  const int up = CleanTransSide(x, B, m, kf, jf, iff, 1);
  const int dn = CleanTransSide(x, B, m, kf, jf, iff, 0);
  if constexpr (D == 2) {
    if (x.sph) {   // (torque form, CleanFaceStress)
      const Real vr = x.vol(m,kf,jf,iff)*x.x1v(m,iff)*x.tg(m,2,jf);
      if constexpr (B == 0) {
        const Real s = x.tg(m,2,jf)/vr*cr*0.25;
        auto bf = [&](const int ii) {   // (clean_lz_bnd_free, CleanFaceStress)
          return x.lzbf && CleanBoundaryFace(x, 0, m, kf, jf, ii);
        };
        const Real fp = bf(iff+1) ? 0.0 : x.area[0](m,kf,jf,iff+1)*x.x1f(m,iff+1)*s;
        const Real fm = bf(iff) ? 0.0 : -x.area[0](m,kf,jf,iff)*x.x1f(m,iff)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        CleanStAdd<rk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<lk, 0, 1>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<rk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        CleanStAdd<lk, 0, -1>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
      } else {
        const Real s = x.x1v(m,iff)/vr*cr*0.25;
        const Real fp = x.area[1](m,kf,jf+1,iff)*x.tg(m,3,jf+1)*s;
        const Real fm = -x.area[1](m,kf,jf,iff)*x.tg(m,3,jf)*s;
        const Real f0 = (up ? fp : 2.0*fp) + (dn ? fm : 2.0*fm);
        CleanStAdd<rk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<lk, 1, 0>(x, ab, m, k, j, i, up ? fp : 0.0, c);
        CleanStAdd<rk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<lk, 0, 0>(x, ab, m, k, j, i, f0, c);
        CleanStAdd<rk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
        CleanStAdd<lk, -1, 0>(x, ab, m, k, j, i, dn ? fm : 0.0, c);
      }
      return;
    }
  }
  if (up + dn == 0) return;
  const Real rb = !x.sph ? 1.0 : ((D == 0) ? x.x1f(m,iff) : x.x1v(m,iff));
  const Real w = (up + dn == 2) ? 1.0 : 2.0;
  const Real ib = !x.sph ? w*x.idx[B] :
                  ((B == 0) ? 2.0/(x.x1v(m,iff+up) - x.x1v(m,iff-dn)) :
                  ((B == 1) ? 2.0/(rb*(x.x2v(m,jf+up) - x.x2v(m,jf-dn))) :
                  w*x.idx[2]/(rb*sin((D == 1) ? x.x2f(m,jf) : x.x2v(m,jf)))));
  const Real q = 0.25*ib*cr;
  if (up) {
    CleanStAdd<rk + bk, rj + bj, ri + bi>(x, ab, m, k, j, i, q, c);
    CleanStAdd<lk + bk, lj + bj, li + bi>(x, ab, m, k, j, i, q, c);
  } else {
    CleanStAdd<rk, rj, ri>(x, ab, m, k, j, i, q, c);
    CleanStAdd<lk, lj, li>(x, ab, m, k, j, i, q, c);
  }
  if (dn) {
    CleanStAdd<rk - bk, rj - bj, ri - bi>(x, ab, m, k, j, i, -q, c);
    CleanStAdd<lk - bk, lj - bj, li - bi>(x, ab, m, k, j, i, -q, c);
  } else {
    CleanStAdd<rk, rj, ri>(x, ab, m, k, j, i, -q, c);
    CleanStAdd<lk, lj, li>(x, ab, m, k, j, i, -q, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStFace
//! \brief the row terms of face SD (0 low, 1 high) normal to D of cell (k,j,i) in its h
//! div F: the face weights al, ar and the cross factor of the pass (x.fwc, CLF_QL) with
//! the off-diagonal stress of CleanFaceCross

template <int D, int SD>
KOKKOS_FORCEINLINE_FUNCTION
void CleanStFace(const CleanArgs &x, const int m, const int k, const int j, const int i,
                 Real (&c)[CLN_ST]) {
  constexpr int ek = (D == 2), ej = (D == 1), ei = (D == 0);
  const int kf = k + SD*ek, jf = j + SD*ej, iff = i + SD*ei;
  const Real wd = (SD ? x.h : -x.h)*CleanDivW(x, D, m, k, j, i, SD);
  c[CleanStSlot((SD - 1)*ek, (SD - 1)*ej, (SD - 1)*ei)] += wd*x.fc[D](m,CLF_QL,kf,jf,iff);
  c[CleanStSlot(SD*ek, SD*ej, SD*ei)] += wd*x.fc[D](m,CLF_QL+1,kf,jf,iff);
  // (the cubed sphere's lateral stress comes through its sums Y, CleanApply)
  if (!x.multi_d || x.cub || CleanBoundaryFace(x, D, m, kf, jf, iff)) return;
  const Real cr = wd*x.fc[D](m,CLF_QL+2,kf,jf,iff);
  if constexpr (D != 0) CleanStCross<D, SD, 0>(x, m, k, j, i, cr, c);
  if constexpr (D != 1) CleanStCross<D, SD, 1>(x, m, k, j, i, cr, c);
  if constexpr (D != 2) {
    if (x.three_d) CleanStCross<D, SD, 2>(x, m, k, j, i, cr, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStencilRow
//! \brief c = the row of cell (k,j,i) of the operator of CleanApply (x.fwc): the diagonal
//! and the coefficients of E in CleanFluxDivergence (clean_stag_stencil)

KOKKOS_INLINE_FUNCTION
void CleanStencilRow(const CleanArgs &x, const int m, const int k, const int j,
                     const int i, Real (&c)[CLN_ST]) {
  for (int s = 0; s < CLN_ST; ++s) c[s] = 0.0;
  const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
  c[0] = 1.0 + x.h*x.c*x.op(m,M1_OP_E,k,j,i)*qr + (x.pd ? x.w(m,CLW_PD,k,j,i) : 0.0);
  CleanStFace<0, 0>(x, m, k, j, i, c);
  CleanStFace<0, 1>(x, m, k, j, i, c);
  if (x.multi_d) {
    CleanStFace<1, 0>(x, m, k, j, i, c);
    CleanStFace<1, 1>(x, m, k, j, i, c);
  }
  if (x.three_d) {
    CleanStFace<2, 0>(x, m, k, j, i, c);
    CleanStFace<2, 1>(x, m, k, j, i, c);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSourceDivergence
//! \brief h div of the E-independent part of the face fluxes of cell (k,j,i): the face
//! sources of the pass (CLF_SR, cln_rhsaf)

KOKKOS_INLINE_FUNCTION
Real CleanSourceDivergence(const CleanArgs &x, const int m, const int k, const int j,
                           const int i) {
  Real d = 0.0;
  for (int n = 0; n < 3; ++n) {
    if (n >= x.nf) break;
    const int dk = (n == 2), dj = (n == 1), di = (n == 0);
    d += CleanDivF(x, n, m, k, j, i, x.fc[n](m,CLF_SR,k,j,i),
                   x.fc[n](m,CLF_SR,k+dk,j+dj,i+di));
  }
  return x.h*d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsArm
//! \brief cubed sphere: the L_z arm of the lateral face (kf,jf,iff) normal to q (1, 2),
//! lambda = (M x n)_z/A of its exact first moment M = int x dA (CLG_M): its momentum
//! p n has the L_z lambda p (clean_cs_lz_tel, CleanCsTorque)

KOKKOS_INLINE_FUNCTION
Real CleanCsArm(const CleanArgs &x, const int q, const int m, const int kf, const int jf,
                const int iff) {
  return (x.g.At<CLG_M>(m,3*(q-1),kf,jf,iff)*x.g.At<CLG_N>(m,3*q+1,kf,jf,iff) -
          x.g.At<CLG_M>(m,3*(q-1)+1,kf,jf,iff)*x.g.At<CLG_N>(m,3*q,kf,jf,iff))/
         x.area[q](m,kf,jf,iff);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagCsVector
//! \brief cubed sphere, staggered: the Cartesian vector u of cell (k,j,i) whose
//! projections p on the unit normals n of its six faces (low and high face of family q at
//! 2 q, 2 q + 1) fit in the least squares of weights w, given c = w p:
//! (sum w n n^T) u = sum c n.  Exact for a uniform vector (p = n.u) on any cell: the face
//! families are not orthogonal and the two faces of a family not parallel, so no face
//! projection is a component of u in a basis (nor are the mean normals' projections).
//! clean_stag_csvec_radial: the lateral faces' normals without their component along the
//! cell's r_hat (a radial field has none on them, a uniform one does), so that a radial
//! u is exact as well: else its fit is low by 1/(1 + 2 sin^2(dalpha/2)), 0.5 % at N = 16.

KOKKOS_INLINE_FUNCTION
void CleanStagCsVector(const CleanArgs &x, const int m, const int k, const int j,
                       const int i, const Real c[6], const Real w[6], Real u[3]) {
  Real a[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Real b[3] = {0.0, 0.0, 0.0};
  for (int q = 0; q < 3; ++q) {
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1), iff = i + s*(q == 0);
      Real n[3], nr = 0.0;
      for (int r = 0; r < 3; ++r) {
        n[r] = x.g.At<CLG_N>(m,3*q+r,kf,jf,iff);
        nr += n[r]*x.g(m,CLG_R+r,k,j,i);
      }
      for (int r = 0; r < 3; ++r)
        n[r] -= (x.cvr && q > 0) ? nr*x.g(m,CLG_R+r,k,j,i) : 0.0;
      for (int r = 0; r < 3; ++r) {
        b[r] += c[2*q+s]*n[r];
        for (int t = 0; t < 3; ++t) a[r][t] += w[2*q+s]*n[r]*n[t];
      }
    }
  }
  // (Cramer: the columns of the inverse are the cross products of the rows; written out,
  // no runtime-indexed local)
  Real ca[3][3];
  ca[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ca[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  ca[0][2] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  ca[1][0] = a[2][1]*a[0][2] - a[2][2]*a[0][1];
  ca[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
  ca[1][2] = a[2][0]*a[0][1] - a[2][1]*a[0][0];
  ca[2][0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ca[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  ca[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];
  const Real det = a[0][0]*ca[0][0] + a[0][1]*ca[0][1] + a[0][2]*ca[0][2];
  for (int t = 0; t < 3; ++t) u[t] = (ca[0][t]*b[0] + ca[1][t]*b[1] + ca[2][t]*b[2])/det;
}

//! \brief staggered: the cell vector u of cell (k,j,i) of its six face values c (as
//! CleanStagCsVector) in the components of CleanStagGasMomentum: the mean of the two
//! faces of each family; on the cubed sphere CleanStagCsVector of weights 1 (the same
//! where the normals are the cell's axes)

KOKKOS_INLINE_FUNCTION
void CleanStagCellVector(const CleanArgs &x, const int m, const int k, const int j,
                         const int i, const Real c[6], Real u[3]) {
  if (x.cub) {
    const Real w[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    CleanStagCsVector(x, m, k, j, i, c, w, u);
    return;
  }
  for (int q = 0; q < 3; ++q) u[q] = 0.5*(c[2*q] + c[2*q+1]);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagBetaE
//! \brief the E-dependent part of the E-row beta terms (Ma et al. 2025, eq. 56) of cell
//! (k,j,i), sb(E) = hc kE [(beta^2/2 - beta.a/c) E + 2 beta.y(E)], y the mean of the two
//! faces' momentum F^H/c without its E-independent part (xh: the face momentum state),
//! E in work vector e (CleanRightHandSide, clean_stag_sbimp)

KOKKOS_INLINE_FUNCTION
Real CleanStagBetaE(const CleanArgs &xh, const int e, const int m, const int k,
                    const int j, const int i) {
  Real b2 = 0.0, ba = 0.0, by = 0.0, fv[3], c[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int d = 0; d < 6; ++d) {   // (the cell's flux vector, CleanStagCellVector)
    const int q = d/2, s = d%2;
    if (q < xh.nf) {
      c[d] = CleanFaceFlux(xh, e, q, m, k+s*(q == 2), j+s*(q == 1), i+s*(q == 0), false);
    }
  }
  CleanStagCellVector(xh, m, k, j, i, c, fv);
  for (int n = 0; n < 3; ++n) {
    const Real bt = xh.w(m,CLW_BT+n,k,j,i);
    if (n < xh.nf) by += bt*fv[n]/xh.c;
    b2 += SQR(bt);
    ba += bt*xh.w(m,CLW_A1+n,k,j,i)/xh.c;
  }
  return xh.h*xh.c*xh.op(m,M1_OP_E,k,j,i)*((0.5*b2 - ba)*xh.w(m,e,k,j,i) + 2.0*by);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanMomentumChange
//! \brief the momentum the radiation gives to cell (k,j,i) along dir over the stage, in
//! conservative form: at each face, minus the left side of the radiation flux equation,
//! -[(F - G)/c^2 + h div(D E)], half of each face to each neighbour (clean_stag_split: in
//! proportion to their chi, as the two half cells absorb).  E is work vector e.  With cf,
//! wf set: the term of each face (low, high) and its weight wt; with tq the face's whole
//! term, before its weight (clean_cs_lz_tel, CleanCsTorque).

template <bool kFs = false>
KOKKOS_INLINE_FUNCTION
Real CleanMomentumChange(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i, Real *mag = nullptr,
                         Real *cf = nullptr, Real *wf = nullptr, Real *tq = nullptr) {
  // (kFs: x.fs and x.fa known set, the face fluxes and terms read, none computed)
  const bool fs = kFs || x.fs, fa = kFs || x.fa;
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  Real dm = 0.0;
  for (int side = 0; side < 2; ++side) {
    const int kf = k + side*dk, jf = j + side*dj, iff = i + side*di;
    const Real f = fs ? x.fc[dir](m,CLF_N,kf,jf,iff) :
                   CleanFaceFluxCall(x, e, dir, m, kf, jf, iff, true);
    const Real g = x.fc[dir](m,CLF_G,kf,jf,iff);
    // (clean_stag_split: shared as the two half cells absorb, chi_c/(chi_L + chi_R))
    const int o = 2*side - 1;
    const Real xc = x.op(m,M1_OP_T,k,j,i), xs = xc + x.op(m,M1_OP_T,k+o*dk,j+o*dj,i+o*di);
    Real wt = (x.hi && xs > 0.0) ? xc/xs : 0.5;
    if (x.hi && (x.sph || x.cub)) {   // (each half cell its absorption of V_f F; cubed
      const int ko = k + o*dk, jo = j + o*dj, io = i + o*di;   //  sphere: the owners' V)
      const Real vc = x.cub ? x.g(m,CLG_F+3,k,j,i) : x.vol(m,k,j,i);
      const Real vo = x.cub ? x.g(m,CLG_F+3,ko,jo,io) : x.vol(m,ko,jo,io);
      const Real xv = xc*vc + x.op(m,M1_OP_T,k+o*dk,j+o*dj,i+o*di)*vo;
      wt = (xv > 0.0) ? 0.5*xc*(vc + vo)/xv : 0.5;
    }
    if (x.hi && CleanBoundaryFace(x, dir, m, kf, jf, iff)) {
      // (a boundary flux is no F equation: the inner half cell's absorption of it;
      // gflux: of an open face the share 1 - exp(-chi dx/2) the half cell absorbs, not
      // chi dx/2; an imposed flux runs through the half cell, which takes chi dx/2 of it,
      // half the face as the production row, rad_m1_implicit.cpp:15272-15276)
      const Real dc = 1.0/CleanDivW(x, dir, m, k, j, i, side), q = 0.5*xc*dc;
      const bool gq = x.gf && CleanOpen(x, m, 2*dir + side);
      const Real db = gq ? -expm1(-q)*x.h*f/(x.c*dc) : 0.5*x.h*xc*f/x.c;
      dm += db;
      if (cf) {
        cf[side] = db;
        wf[side] = 0.5;
      }
      if (tq) tq[side] = 2.0*db;
      continue;
    }
    // (clean_stag_ftie: T of the face as its owner, the block whose low face it is,
    // computed it, CleanStagFaceTie)
    const Real tu = x.tf ? x.w(m,CLW_L+3*dir,kf,jf,iff) :
                    (fa ? x.fc[dir](m,CLF_PT,kf,jf,iff) : ((f - g)/(x.c*x.c) +
                    x.h*CleanFaceStressCall(x, e, dir, m, kf, jf, iff, true) -
                    (fs ? x.h*x.fc[dir](m,CLF_N+1,kf,jf,iff) : 0.0)));
    if (tq) tq[side] = tu;
    Real t = wt*tu;
    // (mcl, the Newton passes: the conservative term equals the face's absorption
    // h chi_f F/c only where the iterate's F satisfies its equation; capped there at the
    // force of a realisable flux, |F| <= c E, as the production deposit, so that a lagged
    // F does not kick floor gas to ~c; the write-back keeps the conservative term)
    if (x.mcl) {
      const Real cap = wt*(1.0/x.fc[dir](m,CLF_TH,kf,jf,iff) - 1.0)/x.c*
                       fmax(fmax(x.w(m,e,kf,jf,iff), x.w(m,e,kf-dk,jf-dj,iff-di)), 0.0);
      const Real tc = fmin(fmax(t, -cap), cap);
      if (tc != t) Kokkos::atomic_add(&x.mclc(0), 1.0);
      t = tc;
    }
    dm -= t;
    if (cf) {
      cf[side] = -t;
      wf[side] = wt;
    }
    // (mag: the size of its terms, for the round-off level of dm)
    auto de = [&](const int kk, const int jj, const int ii) {
      return fabs(x.cs(m,CLN_D11+dir,kk,jj,ii)*x.w(m,e,kk,jj,ii));
    };
    const Real id = (x.sph || x.cub) ? 1.0/x.dxf[dir](m,kf,jf,iff) : x.idx[dir];
    if (mag) *mag += 0.5*((fabs(f) + fabs(g))/(x.c*x.c) + x.h*id*
                          (de(kf,jf,iff) + de(kf-dk,jf-dj,iff-di)));
  }
  return dm;
}

// dissipation weight candidates (<rad_m1>/clean_alpha)
constexpr int CLA_TAU = 0;          // A(tau)^2, A = (1 - exp(-tau))/tau
constexpr int CLA_TAU_SENSOR = 1;   // max(A(tau)^2 psi, s_ray)
constexpr int CLA_TAU_RAY = 2;      // max(A(tau)^2, s_ray)
constexpr int CLA_ONE = 3;          // 1
constexpr int CLA_ZERO = 4;         // 0

//----------------------------------------------------------------------------------------
//! \fn Real CleanDissipationWeight
//! \brief the weight al of the numerical dissipation at a face, from the cell optical
//! depth tau = chi_f dx_d, the smoothness sensor psi and the rays' streaming weight sray,
//! and rc = (1 - al)/tau, with which the steady RC coefficient (1 - al) c/chi is c dx rc
//! (analytic to tau = 0 for clean_alpha = tau).  A^2 makes the dissipation O(1/tau^2)
//! where the physical diffusion flux is O(1/tau), also where the jumps are first order.
//! tau enters only here (numerical dissipation), never the closure.

KOKKOS_INLINE_FUNCTION
Real CleanDissipationWeight(const int mode, const Real tau, const Real psi,
                            const Real sray, Real &rc, const int p = 2) {
  // A = 1 - b tau
  const Real b = (tau > 1.0e-4) ? (1.0 + expm1(-tau)/tau)/tau
                                : 0.5 - tau*(1.0 - 0.25*tau)/6.0;
  const Real a = (tau > 1.0e-4) ? -expm1(-tau)/tau : 1.0 - b*tau;
  // A^p, p = <rad_m1>/clean_dissip_power (temporary test key)
  const Real ap = (p == 3) ? a*a*a : a*a;
  Real al = 0.0;
  switch (mode) {
    case CLA_TAU: al = ap; break;
    case CLA_TAU_SENSOR: al = fmax(ap*psi, sray); break;
    case CLA_TAU_RAY: al = fmax(ap, sray); break;
    case CLA_ONE: al = 1.0; break;
    default: break;
  }
  rc = (mode == CLA_TAU) ? b*(1.0 + a + ((p == 3) ? a*a : 0.0))
                         : (1.0 - al)/fmax(tau, 2.0);
  return al;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDiag
//! \brief component D_ab of the Eddington tensor of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanDiag(const CleanArgs &x, const int a, const int b, const int m, const int k,
               const int j, const int i) {
  return x.cs(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanNf
//! \brief component b of the unit normal of face (kf,jf,iff) normal to dir: e_dir, on the
//! cubed sphere the face's own normal

KOKKOS_INLINE_FUNCTION
Real CleanNf(const CleanArgs &x, const int b, const int dir, const int m, const int kf,
             const int jf, const int iff) {
  return x.cub ? x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff) : ((b == dir) ? 1.0 : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsFrame
//! \brief cubed sphere: the orthonormal frame of y in cell (k,j,i), rows q[0] = r_hat
//! (the hydro's), q[1] = theta_hat = phi_hat x r_hat, q[2] = phi_hat = z x r_hat/|.|, in
//! Cartesian components (y_cart = Q^T y, Q^T = Q^-1).  Returns R = x1v |z x r_hat|, the
//! arm of L_z (CleanInit refuses an odd panel resolution, which has R = 0 cells).

KOKKOS_INLINE_FUNCTION
Real CleanCsFrame(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  Real q[3][3]) {
  for (int c = 0; c < 3; ++c) q[0][c] = x.g.At<CLG_R>(m,c,k,j,i);
  const Real rho = std::sqrt(q[0][0]*q[0][0] + q[0][1]*q[0][1]);
  q[2][0] = -q[0][1]/rho;
  q[2][1] = q[0][0]/rho;
  q[2][2] = 0.0;
  q[1][0] = q[2][1]*q[0][2];
  q[1][1] = -q[2][0]*q[0][2];
  q[1][2] = q[2][0]*q[0][1] - q[2][1]*q[0][0];
  return x.x1v(m,i)*rho;
}

//! \brief the arm R of L_z = R p_phi of a cell vector p in the frame of y: spherical
//! polar x1v sin(theta), the cubed sphere CleanCsFrame's; 0 on Cartesian
KOKKOS_INLINE_FUNCTION
Real CleanLzArm(const CleanArgs &x, const int m, const int k, const int j, const int i) {
  Real q[3][3];
  return x.cub ? CleanCsFrame(x, m, k, j, i, q) : (x.sph ? x.x1v(m,i)*x.tg(m,2,j) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCsFamily
//! \brief cubed sphere: the normals n_d of the face families of cell (k,j,i) for the
//! half-range closure (CLN_HP, CLN_HM, CLN_KP): r_hat of the frame and the mean normals
//! of the two x2 and the two x3 faces

KOKKOS_INLINE_FUNCTION
void CleanCsFamily(const CleanArgs &x, const int m, const int k, const int j, const int i,
                   Real nd[3][3]) {
  for (int d = 0; d < 3; ++d) {
    const int dk = (d == 2), dj = (d == 1);
    Real s = 0.0;
    for (int c = 0; c < 3; ++c) {
      nd[d][c] = (d == 0) ? x.g.At<CLG_R>(m,c,k,j,i) :
                 x.g.At<CLG_N>(m,3*d+c,k,j,i) + x.g.At<CLG_N>(m,3*d+c,k+dk,j+dj,i);
      s += nd[d][c]*nd[d][c];
    }
    for (int c = 0; c < 3; ++c) nd[d][c] /= std::sqrt(s);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanRayFlux
//! \brief the rays' flux factor h = H/J of cell (k,j,i) in the components of the face
//! fluxes from its half-range sums s_b = h+_b + h-_b: s itself; on the cubed sphere
//! (Cartesian) the solution of n_b . h = s_b for the family normals n_b (CleanCsFamily)

KOKKOS_INLINE_FUNCTION
void CleanRayFlux(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  const Real s[3], Real h[3]) {
  for (int b = 0; b < 3; ++b) h[b] = s[b];
  if (!x.cub) return;
  Real nd[3][3], c[3][3];   // c_b = n_{b+1} x n_{b+2}
  CleanCsFamily(x, m, k, j, i, nd);
  for (int b = 0; b < 3; ++b) {
    const Real *p = nd[(b+1)%3], *q = nd[(b+2)%3];
    for (int e = 0; e < 3; ++e) c[b][e] = p[(e+1)%3]*q[(e+2)%3] - p[(e+2)%3]*q[(e+1)%3];
  }
  const Real det = nd[0][0]*c[0][0] + nd[0][1]*c[0][1] + nd[0][2]*c[0][2];
  for (int e = 0; e < 3; ++e) h[e] = (s[0]*c[0][e] + s[1]*c[1][e] + s[2]*c[2][e])/det;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanGasVector, CleanGasMomentum
//! \brief the gas momentum of cell (k,j,i) in the frame of y.  On the cubed sphere the
//! hydro stores m_r = p.r_hat and the COVARIANT components m_a = p.e_a on its unit panel
//! tangents (non-orthogonal, g_ab = e_a.e_b, hydro_am.hpp): p_t = g^ab m_b e_a.

KOKKOS_INLINE_FUNCTION
void CleanGasVector(const CleanArgs &x, const int m, const int k, const int j,
                    const int i, const Real u[3], Real p[3]) {
  for (int n = 0; n < 3; ++n) p[n] = u[n];
  if (!x.cub) return;
  Real q[3][3], pt[3], c = 0.0;
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) c += x.g.At<CLG_E>(m,b,k,j,i)*x.g.At<CLG_E>(m,3+b,k,j,i);
  const Real u1 = (p[1] - c*p[2])/(1.0 - c*c), u2 = (p[2] - c*p[1])/(1.0 - c*c);
  for (int b = 0; b < 3; ++b) {
    pt[b] = u1*x.g.At<CLG_E>(m,b,k,j,i) + u2*x.g.At<CLG_E>(m,3+b,k,j,i);
  }
  for (int a = 1; a < 3; ++a) p[a] = q[a][0]*pt[0] + q[a][1]*pt[1] + q[a][2]*pt[2];
}

//! (CleanGasVector: the same for hydro momentum components u, a linear map)
KOKKOS_INLINE_FUNCTION
void CleanGasMomentum(const CleanArgs &x, const DvceArray5D<Real> &uh, const int m,
                      const int k, const int j, const int i, Real p[3]) {
  const Real u[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
  CleanGasVector(x, m, k, j, i, u, p);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanAdvect
//! \brief a = v + D v of cell (k,j,i) for the gas velocity v in the frame of y (D in
//! Cartesian components on the cubed sphere)

KOKKOS_INLINE_FUNCTION
void CleanAdvect(const CleanArgs &x, const int m, const int k, const int j, const int i,
                 const Real v[3], Real a[3]) {
  if (!x.cub) {
    for (int b = 0; b < 3; ++b) {
      a[b] = v[b];
      for (int c = 0; c < 3; ++c) a[b] += CleanDiag(x, b, c, m, k, j, i)*v[c];
    }
    return;
  }
  Real q[3][3], vc[3], ac[3];
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) vc[b] = v[0]*q[0][b] + v[1]*q[1][b] + v[2]*q[2][b];
  for (int b = 0; b < 3; ++b) {
    ac[b] = vc[b];
    for (int c = 0; c < 3; ++c) ac[b] += CleanDiag(x, b, c, m, k, j, i)*vc[c];
  }
  for (int b = 0; b < 3; ++b) a[b] = q[b][0]*ac[0] + q[b][1]*ac[1] + q[b][2]*ac[2];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanGasChange
//! \brief the change dh of the hydro's momentum components of cell (k,j,i) for a change
//! dm of the gas momentum in the frame of y (the inverse map of CleanGasMomentum): on the
//! cubed sphere dh_r = dm_r and the covariant dh_a = dm_t.e_a

KOKKOS_INLINE_FUNCTION
void CleanGasChange(const CleanArgs &x, const int m, const int k, const int j,
                    const int i, const Real dm[3], Real dh[3]) {
  for (int n = 0; n < 3; ++n) dh[n] = dm[n];
  if (!x.cub) return;
  Real q[3][3], dt[3];
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) dt[b] = dm[1]*q[1][b] + dm[2]*q[2][b];
  for (int a = 1; a < 3; ++a) {
    dh[a] = dt[0]*x.g.At<CLG_E>(m,3*(a-1),k,j,i) +
            dt[1]*x.g.At<CLG_E>(m,3*(a-1)+1,k,j,i) +
            dt[2]*x.g.At<CLG_E>(m,3*(a-1)+2,k,j,i);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCellFaceDistance
//! \brief distance between the centres of the two cells across face (k,j,i) normal to
//! dir: dx_d on a Cartesian mesh; on spherical polar dr, r dtheta and r sin(theta) dphi
//! between the volume centroids, on the cubed sphere the centroids' distance

KOKKOS_INLINE_FUNCTION
Real CleanCellFaceDistance(const CleanArgs &x, const int dir, const int m, const int k,
                           const int j, const int i) {
  return (x.sph || x.cub) ? x.dxf[dir](m,k,j,i) : 1.0/x.idx[dir];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanHcTauW
//! \brief clean_face_hctau / clean_face_jse (old production implicit_thin_hcface = tau):
//! the weights of a cell value interpolated linearly in optical depth from the two cell
//! centres to face (k,j,i) normal to dir, wl = tau_R/(tau_L + tau_R), tau the half-cell
//! depths, x (r/r_f)^2 on the radial faces of sp / cs (exact for a free stream)

KOKKOS_INLINE_FUNCTION
void CleanHcTauW(const CleanArgs &x, const int dir, const int m, const int k, const int j,
                 const int i, Real &wl, Real &wr) {
  const int il = i - (dir == 0), jl = j - (dir == 1), kl = k - (dir == 2);
  const bool rad = (dir == 0 && (x.sph || x.cub));
  const Real dl = rad ? x.x1f(m,i) - x.x1v(m,il) : 1.0;
  const Real dr = rad ? x.x1v(m,i) - x.x1f(m,i) : 1.0;
  const Real tl = fmax(x.op(m,M1_OP_T,kl,jl,il), 0.0)*dl;
  const Real tr = fmax(x.op(m,M1_OP_T,k,j,i), 0.0)*dr;
  wl = (tl + tr > 0.0) ? tr/(tl + tr) : 0.5;
  wr = 1.0 - wl;
  if (rad) {
    wl *= SQR(x.x1v(m,il)/x.x1f(m,i));
    wr *= SQR(x.x1v(m,i)/x.x1f(m,i));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanOverlapFor
//! \brief old production's implicit_halo_overlap: f(m,k,j,i) over the active cells
//! CleanInner by ng while a ghost exchange is in flight, then fin() (its second half),
//! then f over the other active cells.  ov = false (or MeshBlocks too thin): fin(), then
//! f over all of them.  Each cell's value is the same either way.  fd >= 0: f over the
//! faces normal to fd, the inner ones those whose two neighbour faces along fd are the
//! block's (the F solve's, CleanStagFSolve).  post (the producer's side): f over the
//! other cells first, which the exchange sends, then post() (its first half), f over the
//! inner ones while it is in flight, then fin()

// (fin a std::function: nvcc takes no host lambda type in an extended lambda's parent)
template <typename F>
void CleanOverlapFor(const std::string &nm, MeshBlockPack *pp, const bool ov, const F &f,
                     const std::function<void()> &fin, const int fd = -1,
                     const std::function<void()> &post = nullptr) {
  auto &id = pp->pmesh->mb_indcs;
  const int g = id.ng, nmb1 = pp->nmb_thispack - 1, is = id.is, js = id.js, ks = id.ks;
  const bool md = pp->pmesh->multi_d, td = pp->pmesh->three_d;
  const int nx[3] = {id.nx1, id.nx2, id.nx3};
  int n[3], lo[3], hi[3];
  for (int d = 0; d < 3; ++d) {
    const bool ac = (d == 0) || (d == 1 && md) || (d == 2 && td);
    n[d] = nx[d] + (d == fd);
    lo[d] = (fd < 0) ? (ac ? g : 0) : (d == fd);
    hi[d] = (fd < 0) ? (ac ? nx[d] - 1 - g : nx[d] - 1) : nx[d] - 1 - (d == fd);
  }
  const CleanShell sh(n[2], n[1], n[0], 0, lo[0], hi[0], lo[1], hi[1], lo[2], hi[2],
                      true, true);
  if (!ov || sh.Empty()) {
    if (!post) fin();
    par_for(nm, DevExeSpace(), 0, nmb1, ks, id.ke + (fd == 2), js, id.je + (fd == 1), is,
            id.ie + (fd == 0), f);
    if (post) post(), fin();
    return;
  }
  auto inner = [&]() {
    par_for(nm, DevExeSpace(), 0, nmb1, ks + sh.k0, ks + sh.k1, js + sh.j0, js + sh.j1,
            is + sh.i0, is + sh.i1, f);
  };
  if (!post) inner(), fin();
  par_for(nm + "_sh", DevExeSpace(), 0, nmb1, 0, sh.ncell-1,
  KOKKOS_LAMBDA(const int m, const int q) {
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    f(m, ks + k, js + j, is + i);
  });
  if (post) post(), inner(), fin();
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanPredictE
//! \brief the predicted E = e0 + de kept positive: linear down to 0 and e0 exp(de/e0)
//! below (cut: a jump to e0/e at e0 + de = 0), or e0 exp(de/e0) for every de < 0
//! (smooth, clean_predict_pos: continuous, e0 + de to first order)

KOKKOS_INLINE_FUNCTION
Real CleanPredictE(const Real e0, const Real de, const bool smooth) {
  if (smooth) return (de >= 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
  return (e0 + de > 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
}

}  // namespace radm1
#endif  // RAD_M1_RAD_M1_CLEAN_HPP_
