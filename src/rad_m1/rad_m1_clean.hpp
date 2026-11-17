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

namespace radm1 {

// closure variables (arrays cln_c and cln_s)
constexpr int CLN_D11 = 0;   // Eddington tensor D11 D22 D33 D12 D13 D23
constexpr int CLN_HP = 6;    // h+_a = H+_a/J, a = 1,2,3
constexpr int CLN_HM = 9;    // h-_a = H-_a/J (<= 0)
constexpr int CLN_KP = 12;   // k+^d_b = K+^d_bd/J (outgoing stress of mu_d > 0), 3 d + b
constexpr int CLN_C1 = 21;   // non-locality C1 = <|I - S|>/J
constexpr int CLN_J = 22;    // J of the rays (energy-density units)
constexpr int CLN_F2 = 23;   // the rays' flux through the low x2 face (M1_VET_F2, not
                             // normalised; clean_consistency)
constexpr int CLN_F1 = 27;   // their flux through an x1 boundary face (M1_VET_F1)
constexpr int CLN_FX = 31;   // their flux through the low x1 face (M1_VET_FX)
constexpr int CLN_F3 = 35;   // their flux through the low x3 face (M1_VET_F3)
constexpr int CLN_FH = 39;   // cubed sphere: through the high x2, x3 face of the last
                             // cell of a MeshBlock (+4; a seam ghost's low face differs)
constexpr int CLN_G = 47;    // curvature of their stress in the cell layout (y_r, y_theta)
constexpr int CLN_JR = 49;   // their J also below clean_ray_min E (CLN_J is 0 there): the
                             // staggered flux factors (clean_stag_gflux)
constexpr int CLN_N = 50;

// face variables (arrays cln_fc1, cln_fc2, cln_fc3)
constexpr int CLF_W = 0;     // blend weight w
constexpr int CLF_TH = 1;    // theta = 1/(1 + h c chi_f)
constexpr int CLF_G = 2;     // flux of the old vector G_f
constexpr int CLF_HP = 3;    // h+_f
constexpr int CLF_HM = 4;    // h-_f
constexpr int CLF_RC = 5;    // cell layout: RC coefficient (1 - al) c/chi_f (h-free)
constexpr int CLF_N = 6;
// staggered face slots after CLF_N (clean_stag_fsolve: CLF_N .. CLF_N + 3 hold its F,
// the dissipation term, F^H and, with clean_stag_fkin, the energy flux's memory d)
constexpr int CLF_QL = CLF_N + 4;   // the operator's al, ar (CleanFaceWeights) and cross
                                    // factor -(1-w) theta h c^2 of the pass (3)
constexpr int CLF_DS = CLF_N + 7;   // clean_stag_fkin: d of the stage-loop start
constexpr int CLF_DB = CLF_N + 8;   // clean_stag_vimp: drag share b = beta/(1 + beta) and
constexpr int CLF_DQ = CLF_N + 9;   //   the face momentum Q = G - c^2 dm_f of the iterate
constexpr Real CLN_MEMUP = 1.0/6.0;  // clean_stag_gflux: a memory outflow above this
                                     // share of the donor's E_old goes upwind

// work vectors (array cln_w)
constexpr int CLW_EO = 0;    // E of the old vector
constexpr int CLW_EK = 1;    // Newton iterate of E
constexpr int CLW_T = 2;     // gas temperature of the iterate
constexpr int CLW_EI = 3;    // gas internal energy of the old vector
constexpr int CLW_RR = 4;    // exchange linearisation r = s/(d cv + s), s = 4 h c kP a T^3
constexpr int CLW_A1 = 5;    // advective coefficients a_1, a_2, a_3
constexpr int CLW_DM = 8;    // momentum change of the iterate (3)
constexpr int CLW_AM = 11;   // x1-line coefficient of E_{i-1}
constexpr int CLW_A0 = 12;   // diagonal
constexpr int CLW_AP = 13;   // x1-line coefficient of E_{i+1}
constexpr int CLW_CP = 14;   // tridiagonal solve scratch
constexpr int CLW_B = 15;    // right-hand side
constexpr int CLW_X = 16;    // solution
constexpr int CLW_R = 17;    // BiCGStab vectors
constexpr int CLW_R0 = 18;
constexpr int CLW_P = 19;
constexpr int CLW_V = 20;
constexpr int CLW_S = 21;
constexpr int CLW_Q = 22;
constexpr int CLW_Y = 23;
constexpr int CLW_Z = 24;
constexpr int CLW_WK = 25;   // work W of the right-hand side
constexpr int CLW_DK = 26;   // d(rho kappa_P)/dT, d(rho kappa_E)/dT of the iterate
constexpr int CLW_BT = 28;   // beta = v/c of the iterate (3, cell layout: frame of y)
constexpr int CLW_RE = 31;   // h S_E = hc kP G a T^4 - hc kE E + RE E + RY.y - W (Ma et
constexpr int CLW_RY = 32;   //   al. 2025, eq. 56, G = CLW_GM): the beta terms (RY: 3);
                             //   staggered: RE = hc kP a T^4, RY 1 - r
constexpr int CLW_GM = 35;   // G = 1 - beta^2/2
constexpr int CLW_EP = 36;   // -(RE E + RY.y) of the solution (eint of the gas)
constexpr int CLW_JV = 37;   // (1 - r) dA/dbeta/(c^2 rho), dS_p/dbeta/(c^2 rho) of the
constexpr int CLW_GX = 40;   //   iterate (3 + 1): implicit gas velocity, CleanCellDrag;
                             //   staggered clean_stag_sbimp: GX = r (CleanRightHandSide);
                             //   staggered gflux: JV + 1 = hc kE dop (CleanNewtonUpdate)
constexpr int CLW_TM = 41;   // clean_stag_pc = mg: coefficients of E_{j-1}, E_{j+1},
                             //   E_{k-1}, E_{k+1} (4)
constexpr int CLW_L = 45;    // staggered, cubed sphere: the lateral stress sums Y of the
                             // face families q of a cell, CLW_L + 3 q
                             // (CleanStagCsLateral)
constexpr int CLW_N = 54;

//----------------------------------------------------------------------------------------
//! \struct CleanArgs
//! \brief everything the device functions of the clean solve read

struct CleanArgs {
  DvceArray5D<Real> cs;        // closure of the stage
  DvceArray5D<Real> w;         // work vectors
  DvceArray5D<Real> fc[3];     // face variables in x1, x2, x3
  DvceArray3D<Real> beam;      // boundary beam flux at x1min (m, k-ks, j-js)
  DvceArray4D<Real> bstress;   // beam stress column c K_b1 at x1min (m, b, k-ks, j-js)
  DvceArray5D<Real> u;         // cell vectors of the cell layout (CLC_*)
  DvceArray5D<Real> ff;        // face fluxes, 4 per direction (ffon: CleanCellDivergence)
  bool ffon = false;
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
  bool recon;                  // cell layout: reconstructed face states in the jumps
  bool all_low;                // cell layout: every interior face low-order
  bool consist;                // cell layout: face fluxes carry the consistency bracket
  bool milne;                  // cell layout: sub-cell Eddington face model
  bool fvrays;                 // cell layout: FV rays (clean_rays = fv1 | fv)
  bool fd;                     // staggered, clean_stag_fkin: the energy flux's
                               // G + F^H - F
  bool gf;                     // staggered, clean_stag_gflux: flux factors, low order
  bool fs;                     // staggered, clean_stag_fsolve: the face F of the F solve
                               // (CLF_N) with its dissipation term (CLF_N + 1)
  bool hi;                     // staggered, clean_stag_split: the face momentum state
                               // F^H (w = 0), not the energy flux
  bool drag;                   // staggered, clean_stag_vimp: implicit gas velocity
  bool upw;                    // staggered, clean_stag_upwind: donor-cell advective flux
  bool fwc;                    // staggered: al, ar and the cross factor from CLF_QL
  bool sbi;                    // staggered, clean_stag_sbimp: E-row beta terms implicit
  bool sar;                    // staggered: the streaming share's net flux in area form
                               // (CleanStagArea)
  DvceArray5D<Real> op;        // opacities (M1_OP_*)
  int nf;                      // cell layout: number of flux components (dimensions)
  int is, ie, js, je, ks, ke;
  bool multi_d, three_d;
  // spherical-polar geometry (cell layout): x1 = r, x2 = theta, x3 = phi, vectors in the
  // local orthonormal basis
  bool sph;
  DvceArray4D<Real> vol;       // cell volumes
  DvceArray4D<Real> area[3];   // face areas normal to x1, x2, x3
  DvceArray4D<Real> dxf[3];    // distances between the cell centres across those faces
  DvceArray2D<Real> x1v, x1f;  // cell-centre (volume centroid) and face radii
  DvceArray2D<Real> x2v, x2f;  // cell-centre and face polar angles
  Real cseam;                  // cot of a theta-periodic seam face: its two copies' mean
  // cubed sphere (cell layout): y in the cell's orthonormal frame (r, theta, phi) about
  // the hydro's radius vector; vol, area, dxf and x1v, x1f as above, g the geometry
  // (CLG_*) of CleanCsGeometry
  bool cub;
  DvceArray5D<Real> g;
  // the conserved gas energy also holds d phi (<hydro>/etotgrav) and |B|^2/2 (<mhd>)
  bool etg, mhd;
  DvceArray4D<Real> phi, emag;
};

// the source moments of the FV rays (array cln_fvq; lld: all, else chi and S)
constexpr int CLQ_X = 0;     // chi
constexpr int CLQ_J = 1;     // J of the source
constexpr int CLQ_E = 2;     // the gas's emission q = S - (1 - eps) J
constexpr int CLQ_S = 3;     // slopes of S along x1, x2, x3 (+3)
constexpr int CLQ_W = 6;     // 1 - eps
constexpr int CLQ_C = 7;     // curvatures Q of S along x1, x2, x3 (+3)
constexpr int CLQ_S0 = 10;   // S (every FV mode)
constexpr int CLQ_N = 11;

// cubed-sphere geometry (array cln_gn), Cartesian components
constexpr int CLG_N = 0;     // unit normal of the lower face along dir: CLG_N + 3 dir
constexpr int CLG_X = 9;     // cell centroid
constexpr int CLG_R = 12;    // r_hat of the cell centre (the hydro's, hydro_am.hpp)
constexpr int CLG_E = 15;    // the hydro's unit panel tangents e_xi, e_eta (+3)
constexpr int CLG_M = 21;    // first moment int x dA of the lower x2, x3 face (+3)
constexpr int CLG_F = 27;    // staggered, exchanged (a seam ghost's are its owner's): the
                             // normal of the lower x2 face, the volume (+3), the normal
                             // of the lower x3 face (+4) and the MeshBlock gid (+7)
constexpr int CLG_A = 35;    // staggered, exchanged: the mean area of the cell's two
                             // faces of family x1, x2, x3 (+2), CleanStagArea
constexpr int CLG_NV = 38;

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
//! \fn bool CleanInner
//! \brief cell (k,j,i) lies inside the active cells, at least g from every block face
//! along the mesh dimensions (md: multi-D, td: 3-D)

KOKKOS_INLINE_FUNCTION
bool CleanInner(const int k, const int j, const int i, const int g, const int is,
                const int ie, const int js, const int je, const int ks, const int ke,
                const bool md, const bool td) {
  return i >= is + g && i <= ie - g && (!md || (j >= js + g && j <= je - g)) &&
         (!td || (k >= ks + g && k <= ke - g));
}

//----------------------------------------------------------------------------------------
//! \struct CleanLevel
//! \brief one grid level of the cell-layout preconditioner: an array holding the 4x4
//! blocks of the first-order operator (diagonal block and the blocks of the two
//! neighbours along each direction) and 4-component vectors, with their slots

struct CleanLevel {
  DvceArray5D<Real> u;         // blocks and vectors (m, slot, k, j, i)
  int o[3];                    // index of the first active cell along x1, x2, x3
  int n[3];                    // active cells along x1, x2, x3
  int a0;                      // diagonal block
  int sm[3], sp[3];            // blocks of the lower and upper neighbour along each dir
};

// slots of a coarse level (arrays cln_mg, no ghost cells)
constexpr int CLM_A0 = 0;    // diagonal block (row major)
constexpr int CLM_M = 16;    // block of the lower neighbour along dir: CLM_M + 32 dir
constexpr int CLM_P = 32;    // block of the upper neighbour along dir: CLM_P + 32 dir
constexpr int CLM_B = 112;   // right-hand side
constexpr int CLM_X = 116;   // solution
constexpr int CLM_T1 = 120;  // residual
constexpr int CLM_T2 = 124;  // line-solve correction
constexpr int CLM_N = 128;

// slots of a coarse level of the staggered mg (arrays cln_mg): the x1-line coefficients,
// diagonal and Thomas scratch in the order of cln_w (CLW_AM..CLW_CP), then as CLW_TM
constexpr int CLS_A0 = 1;    // diagonal (CLS_A0 -+ 1: E_{i-1}, E_{i+1}; + 2: scratch)
constexpr int CLS_TM = 4;    // E_{j-1}, E_{j+1}, E_{k-1}, E_{k+1}
constexpr int CLS_B = 8;     // right-hand side
constexpr int CLS_X = 9;     // solution
constexpr int CLS_N = 10;

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
//! below T(d, rhs).  Safeguarded Newton iteration, general EOS.

KOKKOS_INLINE_FUNCTION
Real CleanExchangeTemperature(const EOS_Data &eos, const Real d, const Real a,
                              const Real rhs, const Real tguess) {
  if (!(rhs > 0.0)) return 0.0;
  Real lo = 0.0, hi = eos.Temperature(d, rhs);
  Real t = (tguess > 0.0 && tguess < hi) ? tguess : 0.5*hi;
  for (int n = 0; n < 100; ++n) {
    Real e, p, chi_rho, chi_t, cv;
    eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
    const Real t3 = t*t*t;
    const Real g = e + a*t3*t - rhs;
    if (g > 0.0) {
      hi = t;
    } else {
      lo = t;
    }
    Real tn = t - g/(d*cv + 4.0*a*t3);
    if (!(tn > lo && tn < hi)) tn = 0.5*(lo + hi);
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
           (sm ? 0.5*x.cseam/x.x1v(m,i) : 0.5/(tan(x.x2f(m,j))*x.x1v(m,i)));
  }
  if (q == 1) return 1.0/((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
  const int kc = (q == 2) ? k - (dir == 2) : k, jc = (q == 2) ? j - (dir == 1) : j;
  const int ic = (q == 2) ? i - (dir == 0) : i;
  auto d = [&](const int n) {return x.cs(m,CLN_D11+n,kc,jc,ic);};
  if (dir == 0) {
    return 2.0*d(0) - d(1) - d(2) + (x.multi_d ? d(3)/tan(x.x2v(m,jc)) : 0.0);
  }
  if (sm) return (d(1) - d(2))*x.cseam + 3.0*d(3);
  return (d(1) - d(2))/tan(x.x2f(m,j)) + 3.0*d(3);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDivW
//! \brief weight of the flux through the low (side 0) or high (side 1) face along dir in
//! the divergence of cell (k,j,i): 1/dx_d, on spherical polar A_f/V

KOKKOS_INLINE_FUNCTION
Real CleanDivW(const CleanArgs &x, const int dir, const int m, const int k, const int j,
               const int i, const int side) {
  if (!x.sph && !x.cub) return x.idx[dir];
  return x.area[dir](m,k+side*(dir == 2),j+side*(dir == 1),i+side*(dir == 0))/
         x.vol(m,k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagDnn
//! \brief cubed sphere: n.D.n of cell (k,j,i), n the unit normal of face (kf,jf,iff)
//! normal to dir (D Cartesian), the normal stress of the face-normal momentum

KOKKOS_INLINE_FUNCTION
Real CleanStagDnn(const CleanArgs &x, const int dir, const int m, const int kf,
                  const int jf, const int iff, const int k, const int j, const int i) {
  Real d = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      d += x.g(m,CLG_N+3*dir+a,kf,jf,iff)*x.g(m,CLG_N+3*dir+b,kf,jf,iff)*
           x.cs(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagMemUp
//! \brief clean_stag_gflux: the memory (1 - w) theta (G + d) of the energy flux through
//! face (k,j,i) normal to dir taken implicitly upwind, as q E_donor with q = memory/E_old
//! of the donor, |q| <= c (the donor-cell form of production's postr-1031): on a
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
  const Real out = x.h*fabs(g)*(up ? CleanDivW(x, dir, m, kl, jl, il, 1) :
                                     CleanDivW(x, dir, m, k, j, i, 0));
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
int CleanCsSlot(const CleanArgs &x, const int dir, const int m, const int kf,
                const int jf, const int iff, const int k, const int j, const int i,
                const bool far, Real &s) {
  s = 1.0;
  const bool up = (dir == 0) ? (iff > i) : ((dir == 1) ? (jf > j) : (kf > k));
  if (dir == 0) return (up != far) ? 1 : 0;
  Real t2 = 0.0, t3 = 0.0;
  for (int c = 0; c < 3; ++c) {
    t2 += x.g(m,CLG_F+c,k,j,i)*x.g(m,CLG_N+3*dir+c,kf,jf,iff);
    t3 += x.g(m,CLG_F+4+c,k,j,i)*x.g(m,CLG_N+3*dir+c,kf,jf,iff);
  }
  const int q = (fabs(t2) > fabs(t3)) ? 1 : 2;
  s = (((q == 1) ? t2 : t3) > 0.0) ? 1.0 : -1.0;
  return 2*q + (((up == (s > 0.0)) != far) ? 1 : 0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagCsLat
//! \brief cubed sphere: n.(Y_p + Y_p') of cell (k,j,i) (CleanStagCsLateral), n the
//! normal of its face (kf,jf,iff) normal to dir and p, p' the cell's two face families
//! other than the face's, in the owner's axes (CleanCsSlot)

KOKKOS_INLINE_FUNCTION
Real CleanStagCsLat(const CleanArgs &x, const int dir, const int m, const int kf,
                    const int jf, const int iff, const int k, const int j, const int i) {
  Real s, d = 0.0;
  const int q = CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2;
  const int p = (q == 0) ? 1 : 0, pp = (q == 2) ? 1 : 2;
  for (int c = 0; c < 3; ++c) {
    d += x.g(m,CLG_N+3*dir+c,kf,jf,iff)*(x.w(m,CLW_L+3*p+c,k,j,i) +
                                        x.w(m,CLW_L+3*pp+c,k,j,i));
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
  for (int b = 0; b < 3; ++b) d += x.g(m,CLG_N+3*dir+b,kf,jf,iff)*x.w(m,v+b,k,j,i);
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
//! intensities (clean_rays = sc): the FV rays' (vet_gd, fv1 | fv | lld) are half-range
//! fluxes through the face, which their own finite-volume transport already carries in
//! area form, and so are the flux factors of clean_stag_gflux.

KOKKOS_INLINE_FUNCTION
void CleanStagArea(const CleanArgs &x, const int dir, const int m, const int k,
                   const int j, const int i, const Real w, Real &pl, Real &pr) {
  if (!x.sar || x.hi || !(x.sph || x.cub) || (!x.cub && dir == 2)) return;
  const Real s = w*x.c*(x.fc[dir](m,CLF_HP,k,j,i) + x.fc[dir](m,CLF_HM,k,j,i));
  const int d = (s > 0.0) ? 1 : 0, id = i - d*(dir == 0), jd = j - d*(dir == 1);
  const int kd = k - d*(dir == 2);
  Real sd, a;
  if (x.cub) {
    a = x.g(m,CLG_A + CleanCsSlot(x, dir, m, k, j, i, kd, jd, id, false, sd)/2,kd,jd,id)/
        x.area[dir](m,k,j,i);
  } else {
    a = (dir == 0) ? x.x1v(m,id)*x.x1v(m,id)/(x.x1f(m,i)*x.x1f(m,i))
                   : sin(x.x2v(m,jd))/sin(x.x2f(m,j));
  }
  ((d == 1) ? pl : pr) += s*(a - 1.0);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceWeights
//! \brief coefficients of E_L and E_R in the flux through face (k,j,i) normal to dir, i.e.
//! everything but the off-diagonal Eddington terms and the old-vector flux G.  Physical
//! boundaries (non-periodic) take the boundary flux of their form (none at a wall or an
//! imposed flux, CleanFaceSource).  up (if given) is the donor-cell share phi of
//! clean_stag_upwind.

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
  if (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC) {
    if (!CleanOpen(x, m, 2*dir)) return;   // (a wall; an imposed flux is all source)
    const bool beam = x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0;
    if (beam && x.marshak) return;
    ar = x.marshak ? -x.c*x.marshak_q : x.c*x.cs(m,CLN_HM+dir,k,j,i);
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
  const Real g = (1.0 - w)*th*x.h*x.c*x.c*x.idx[dir];
  al = g*x.cs(m,CLN_D11+dir,kl,jl,il) + w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
  ar = -g*x.cs(m,CLN_D11+dir,k,j,i) + w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  if (x.sph && !x.cub) {   // (centre distance and the curvature terms, CleanStagSph)
    const Real gs = (1.0 - w)*th*x.h*x.c*x.c, id = CleanStagSph(x, 0, dir, m, k, j, i);
    const Real hr = 0.5*CleanStagSph(x, 1, dir, m, k, j, i);
    al = gs*(id*x.cs(m,CLN_D11+dir,kl,jl,il) - hr*CleanStagSph(x, 2, dir, m, k, j, i)) +
         w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
    ar = -gs*(id*x.cs(m,CLN_D11+dir,k,j,i) + hr*CleanStagSph(x, 3, dir, m, k, j, i)) +
         w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  }
  if (x.cub) {   // (cubed sphere: the face-normal momentum, A n.D.n E over V_f)
    const Real gs = (1.0 - w)*th*x.h*x.c*x.c*x.area[dir](m,k,j,i)/
                    (0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,il)));
    al = gs*CleanStagDnn(x, dir, m, k, j, i, kl, jl, il) +
         w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
    ar = -gs*CleanStagDnn(x, dir, m, k, j, i, k, j, i) + w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  }
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
    // (clean_stag_sbimp: of (1 - r) kE E, the emission's response through the gas T)
    const Real rl = x.sbi ? 1.0 - x.w(m,CLW_GX,kl,jl,il) : 1.0;
    const Real rr = x.sbi ? 1.0 - x.w(m,CLW_GX,k,j,i) : 1.0;
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
//! \fn Real CleanFaceStress
//! \brief divergence of the radiation stress D E at face (k,j,i) normal to dir,
//! d_b(D_{dir,b} E): the normal derivative from the two adjacent cells (if normal is
//! set) and, away from physical boundaries, the transverse derivatives averaged over the
//! four cells around the face.  E is work vector e.

KOKKOS_INLINE_FUNCTION
Real CleanFaceStress(const CleanArgs &x, const int e, const int dir, const int m,
                     const int k, const int j, const int i, const bool normal) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real div = 0.0;
  if (x.cub) {   // (cubed sphere: n.div(D E) of the two half cells, the normal stress
    if (normal) {  //  and the lateral sums CleanStagCsLat, over V_f)
      div = x.area[dir](m,k,j,i)*(CleanStagDnn(x, dir, m, k, j, i, k, j, i)*x.w(m,e,k,j,i)
            - CleanStagDnn(x, dir, m, k, j, i, kl, jl, il)*x.w(m,e,kl,jl,il));
    }
    if (!CleanBoundaryFace(x, dir, m, k, j, i)) {
      div += CleanStagCsLat(x, dir, m, k, j, i, kl, jl, il) +
             CleanStagCsLat(x, dir, m, k, j, i, k, j, i);
    }
    return div/(0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,il)));
  }
  if (normal && x.sph) {
    div = CleanStagSph(x, 0, dir, m, k, j, i)*(x.cs(m,CLN_D11+dir,k,j,i)*x.w(m,e,k,j,i) -
          x.cs(m,CLN_D11+dir,kl,jl,il)*x.w(m,e,kl,jl,il)) +
          0.5*CleanStagSph(x, 1, dir, m, k, j, i)*
          (CleanStagSph(x, 2, dir, m, k, j, i)*x.w(m,e,kl,jl,il) +
           CleanStagSph(x, 3, dir, m, k, j, i)*x.w(m,e,k,j,i));
  } else if (normal) {
    div = x.idx[dir]*(x.cs(m,CLN_D11+dir,k,j,i)*x.w(m,e,k,j,i) -
                      x.cs(m,CLN_D11+dir,kl,jl,il)*x.w(m,e,kl,jl,il));
  }
  if (!x.multi_d || CleanBoundaryFace(x, dir, m, k, j, i)) return div;
  for (int b = 0; b < 3; ++b) {
    if (b == dir || (b == 1 && !x.multi_d) || (b == 2 && !x.three_d)) continue;
    const int ab = (dir + b == 1) ? 3 : ((dir + b == 2) ? 4 : 5);   // D12, D13, D23
    const int dk = (b == 2), dj = (b == 1), di = (b == 0);
    auto de = [&](const int kk, const int jj, const int ii) {
      return x.cs(m,ab,kk,jj,ii)*x.w(m,e,kk,jj,ii);
    };
    // (spherical polar: d_r on theta faces, (1/r) d_theta on r faces, central)
    // (torque form: dual r faces i, i+1 and theta faces j, j+1)
    if (x.sph && dir == 2) {
      const Real vr = x.vol(m,k,j,i)*x.x1v(m,i)*sin(x.x2v(m,j));
      auto pe = [&](const int jj, const int ii) {   // (D E)_{b,ph} at a dual face
        return 0.25*(de(k,jj,ii) + de(k-1,jj,ii) + de(k,jj-(b == 1),ii-(b == 0)) +
                     de(k-1,jj-(b == 1),ii-(b == 0)));
      };
      if (b == 0) {
        div += (x.area[0](m,k,j,i+1)*x.x1f(m,i+1)*pe(j,i+1) -
                x.area[0](m,k,j,i)*x.x1f(m,i)*pe(j,i))*sin(x.x2v(m,j))/vr;
      } else {
        div += (x.area[1](m,k,j+1,i)*sin(x.x2f(m,j+1))*pe(j+1,i) -
                x.area[1](m,k,j,i)*sin(x.x2f(m,j))*pe(j,i))*x.x1v(m,i)/vr;
      }
      continue;
    }
    const Real rb = !x.sph ? 1.0 : ((dir == 0) ? x.x1f(m,i) : x.x1v(m,i));
    const Real ib = !x.sph ? x.idx[b] : ((b == 0) ? 2.0/(x.x1v(m,i+1) - x.x1v(m,i-1)) :
                    ((b == 1) ? 2.0/(rb*(x.x2v(m,j+1) - x.x2v(m,j-1))) :
                    x.idx[2]/(rb*sin((dir == 1) ? x.x2f(m,j) : x.x2v(m,j)))));
    div += 0.25*ib*(de(k+dk,j+dj,i+di) + de(kl+dk,jl+dj,il+di)
                    - de(k-dk,j-dj,i-di) - de(kl-dk,jl-dj,il-di));
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
  if (!x.multi_d || CleanBoundaryFace(x, dir, m, k, j, i)) return 0.0;
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

KOKKOS_INLINE_FUNCTION
Real CleanStagAdvectCorr(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  auto e = [&](const int o) {return x.w(m,CLW_EK,k+o*dk,j+o*dj,i+o*di);};
  const Real ul = CleanStagN(x, CLW_A1, dir, m, k, j, i, k-dk, j-dj, i-di);
  const Real ur = CleanStagN(x, CLW_A1, dir, m, k, j, i, k, j, i);
  Real ef, dum;
  if (ul + ur > 0.0) {
    PLM(e(-2), e(-1), e(0), ef, dum);
    return 0.5*(ul + ur)*ef - ul*e(-1);
  }
  PLM(e(-1), e(0), e(1), dum, ef);
  return 0.5*(ul + ur)*ef - ur*e(0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceSource
//! \brief the part of the face flux that does not depend on E: (1-w) theta G inside, the
//! incoming flux (beam, bath) or the imposed flux on physical boundaries

KOKKOS_INLINE_FUNCTION
Real CleanFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i) {
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
    // (clean_stag_sbimp: less the emission's linear response at E_k, CleanFaceWeights)
    if (x.sbi) {
      auto rk = [&](const int kk, const int jj, const int ii) {
        return x.w(m,CLW_GX,kk,jj,ii)*x.op(m,M1_OP_E,kk,jj,ii)*x.w(m,CLW_EK,kk,jj,ii);
      };
      s -= 0.25*tg*x.c*(CleanStagN(x, CLW_BT, dir, m, k, j, i, kl, jl, il) +
                        CleanStagN(x, CLW_BT, dir, m, k, j, i, k, j, i))*x.h*x.c*
           (rk(kl, jl, il) + rk(k, j, i));
    }
  }
  if (x.drag) s += (1.0 - w)*qd*x.fc[dir](m,CLF_DQ,k,j,i);
  // (the advective share 1 - theta_d of an upwind face, CleanFaceWeights)
  if (x.upw && x.feedback) {
    Real al, ar, phi = 0.0;
    CleanFaceWeights(x, dir, m, k, j, i, al, ar, &phi);
    if (phi > 0.0) {
      s += phi*((1.0 - w) - th)*(1.0 - qd)*CleanStagAdvectCorr(x, dir, m, k, j, i);
    }
  }
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFlux
//! \brief flux through face (k,j,i) normal to dir for E in work vector e; with_source adds
//! the E-independent part

CLEAN_DEVICE_NOINLINE KOKKOS_INLINE_FUNCTION
Real CleanFaceFlux(const CleanArgs &x, const int e, const int dir, const int m,
                   const int k, const int j, const int i, const bool with_source) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  Real al, ar;
  CleanFaceWeights(x, dir, m, k, j, i, al, ar);
  Real f = CleanFaceCross(x, e, dir, m, k, j, i);
  if (al != 0.0) f += al*x.w(m,e,kl,jl,il);
  if (ar != 0.0) f += ar*x.w(m,e,k,j,i);
  if (with_source) f += CleanFaceSource(x, dir, m, k, j, i);
  return f;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFluxDivergence
//! \brief h div F of cell (k,j,i) for E in work vector e

KOKKOS_INLINE_FUNCTION
Real CleanFluxDivergence(const CleanArgs &x, const int e, const int m, const int k,
                         const int j, const int i, const bool with_source) {
  if (x.sph || x.cub) {   // (finite volume, sum A_f F_f/V)
    Real d = 0.0;
    for (int n = 0; n < 3; ++n) {   // (constant trip: x.fc[n], x.area[n] resolved)
      if (n >= x.nf) break;
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      d += CleanDivW(x, n, m, k, j, i, 1)*CleanFaceFlux(x, e, n, m, k+dk, j+dj, i+di,
                                                         with_source) -
           CleanDivW(x, n, m, k, j, i, 0)*CleanFaceFlux(x, e, n, m, k, j, i, with_source);
    }
    return x.h*d;
  }
  Real div = x.idx[0]*(CleanFaceFlux(x, e, 0, m, k, j, i+1, with_source) -
                       CleanFaceFlux(x, e, 0, m, k, j, i, with_source));
  if (x.multi_d) {
    div += x.idx[1]*(CleanFaceFlux(x, e, 1, m, k, j+1, i, with_source) -
                     CleanFaceFlux(x, e, 1, m, k, j, i, with_source));
  }
  if (x.three_d) {
    div += x.idx[2]*(CleanFaceFlux(x, e, 2, m, k+1, j, i, with_source) -
                     CleanFaceFlux(x, e, 2, m, k, j, i, with_source));
  }
  return x.h*div;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanSourceDivergence
//! \brief h div of the E-independent part of the face fluxes of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanSourceDivergence(const CleanArgs &x, const int m, const int k, const int j,
                           const int i) {
  if (x.sph || x.cub) {
    Real d = 0.0;
    for (int n = 0; n < 3; ++n) {
      if (n >= x.nf) break;
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      d += CleanDivW(x, n, m, k, j, i, 1)*CleanFaceSource(x, n, m, k+dk, j+dj, i+di) -
           CleanDivW(x, n, m, k, j, i, 0)*CleanFaceSource(x, n, m, k, j, i);
    }
    return x.h*d;
  }
  Real div = x.idx[0]*(CleanFaceSource(x, 0, m, k, j, i+1) -
                       CleanFaceSource(x, 0, m, k, j, i));
  if (x.multi_d) {
    div += x.idx[1]*(CleanFaceSource(x, 1, m, k, j+1, i) -
                     CleanFaceSource(x, 1, m, k, j, i));
  }
  if (x.three_d) {
    div += x.idx[2]*(CleanFaceSource(x, 2, m, k+1, j, i) -
                     CleanFaceSource(x, 2, m, k, j, i));
  }
  return x.h*div;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagCsVector
//! \brief cubed sphere, staggered: the Cartesian vector u of cell (k,j,i) whose
//! projections p on the unit normals n of its six faces (low and high face of family q at
//! 2 q, 2 q + 1) fit in the least squares of weights w, given c = w p:
//! (sum w n n^T) u = sum c n.  Exact for a uniform vector (p = n.u) on any cell: the face
//! families are not orthogonal and the two faces of a family not parallel, so no face
//! projection is a component of u in a basis (nor are the mean normals' projections).

KOKKOS_INLINE_FUNCTION
void CleanStagCsVector(const CleanArgs &x, const int m, const int k, const int j,
                       const int i, const Real c[6], const Real w[6], Real u[3]) {
  Real a[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Real b[3] = {0.0, 0.0, 0.0};
  for (int q = 0; q < 3; ++q) {
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1), iff = i + s*(q == 0);
      for (int r = 0; r < 3; ++r) {
        const Real nr = x.g(m,CLG_N+3*q+r,kf,jf,iff);
        b[r] += c[2*q+s]*nr;
        for (int t = 0; t < 3; ++t) a[r][t] += w[2*q+s]*nr*x.g(m,CLG_N+3*q+t,kf,jf,iff);
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

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagBetaE
//! \brief the E-dependent part of the E-row beta terms (Ma et al. 2025, eq. 56) of cell
//! (k,j,i), sb(E) = hc kE [(beta^2/2 - beta.a/c) E + 2 beta.y(E)], y the mean of the two
//! faces' momentum F^H/c without its E-independent part (xh: the face momentum state),
//! E in work vector e (CleanRightHandSide, clean_stag_sbimp)

KOKKOS_INLINE_FUNCTION
Real CleanStagBetaE(const CleanArgs &xh, const int e, const int m, const int k,
                    const int j, const int i) {
  Real b2 = 0.0, ba = 0.0, by = 0.0, fv[3];
  if (xh.cub) {   // (the cell's flux vector, CleanStagCsVector, as CleanRightHandSide)
    Real c[6], w[6];
    for (int d = 0; d < 6; ++d) {
      const int q = d/2, s = d%2;
      c[d] = CleanFaceFlux(xh, e, q, m, k+s*(q == 2), j+s*(q == 1), i+s*(q == 0), false);
      w[d] = 1.0;
    }
    CleanStagCsVector(xh, m, k, j, i, c, w, fv);
  }
  for (int n = 0; n < 3; ++n) {
    const Real bt = xh.w(m,CLW_BT+n,k,j,i);
    const int dk = (n == 2), dj = (n == 1), di = (n == 0);
    if (n < xh.nf && xh.cub) {
      by += bt*fv[n]/xh.c;
    } else if (n < xh.nf) {
      by += bt*0.5*(CleanFaceFlux(xh, e, n, m, k, j, i, false) +
                    CleanFaceFlux(xh, e, n, m, k+dk, j+dj, i+di, false))/xh.c;
    }
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
//! wf set: the term of each face (low, high) and its weight wt.

KOKKOS_INLINE_FUNCTION
Real CleanMomentumChange(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i, Real *mag = nullptr,
                         Real *cf = nullptr, Real *wf = nullptr) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  Real dm = 0.0;
  for (int side = 0; side < 2; ++side) {
    const int kf = k + side*dk, jf = j + side*dj, iff = i + side*di;
    const Real f = x.fs ? x.fc[dir](m,CLF_N,kf,jf,iff) :
                   CleanFaceFlux(x, e, dir, m, kf, jf, iff, true);
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
      // gflux: the share 1 - exp(-chi dx/2) the half cell absorbs, not chi dx/2)
      const Real dc = 1.0/CleanDivW(x, dir, m, k, j, i, side), q = 0.5*xc*dc;
      const Real db = x.gf ? -expm1(-q)*x.h*f/(x.c*dc) : 0.5*x.h*xc*f/x.c;
      dm += db;
      if (cf) {
        cf[side] = db;
        wf[side] = 0.5;
      }
      continue;
    }
    const Real t = wt*((f - g)/(x.c*x.c) + x.h*CleanFaceStress(x, e, dir, m, kf, jf, iff,
                   true) - (x.fs ? x.h*x.fc[dir](m,CLF_N+1,kf,jf,iff) : 0.0));
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

}  // namespace radm1
#endif  // RAD_M1_RAD_M1_CLEAN_HPP_
