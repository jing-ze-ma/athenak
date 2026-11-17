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
constexpr int CLN_N = 49;

// face variables (arrays cln_fc1, cln_fc2, cln_fc3)
constexpr int CLF_W = 0;     // blend weight w
constexpr int CLF_TH = 1;    // theta = 1/(1 + h c chi_f)
constexpr int CLF_G = 2;     // flux of the old vector G_f
constexpr int CLF_HP = 3;    // h+_f
constexpr int CLF_HM = 4;    // h-_f
constexpr int CLF_RC = 5;    // cell layout: RC coefficient (1 - al) c/chi_f (h-free)
constexpr int CLF_N = 6;

// work vectors (array cln_w); the cell layout allocates the first CLW_NC of them
constexpr int CLW_EO = 0;    // E of the old vector
constexpr int CLW_EK = 1;    // Newton iterate of E
constexpr int CLW_T = 2;     // gas temperature of the iterate
constexpr int CLW_EI = 3;    // gas internal energy of the old vector
constexpr int CLW_RR = 4;    // exchange linearisation r = s/(d cv + s), s = 4 h c kP a T^3
constexpr int CLW_A1 = 5;    // advective coefficients a_1, a_2, a_3
constexpr int CLW_DM = 8;    // momentum change of the iterate (3)
constexpr int CLW_B = 11;    // right-hand side
constexpr int CLW_X = 12;    // solution
constexpr int CLW_WK = 13;   // work W of the right-hand side
constexpr int CLW_DK = 14;   // d(rho kappa_P)/dT, d(rho kappa_E)/dT of the iterate
constexpr int CLW_BT = 16;   // cell layout: beta = v/c of the iterate (3, frame of y)
constexpr int CLW_RE = 19;   // h S_E = hc kP G a T^4 - hc kE E + RE E + RY.y - W (Ma et
constexpr int CLW_RY = 20;   //   al. 2025, eq. 56, G = CLW_GM): the beta terms (RY: 3)
constexpr int CLW_GM = 23;   // G = 1 - beta^2/2
constexpr int CLW_EP = 24;   // -(RE E + RY.y) of the solution (eint of the gas)
constexpr int CLW_JV = 25;   // (1 - r) dA/dbeta/(c^2 rho), dS_p/dbeta/(c^2 rho) of the
constexpr int CLW_GX = 28;   //   iterate (3 + 1): implicit gas velocity, CleanCellDrag
constexpr int CLW_NC = 29;
constexpr int CLW_AM = 29;   // staggered layout: x1-line coefficient of E_{i-1}
constexpr int CLW_A0 = 30;   // diagonal
constexpr int CLW_AP = 31;   // x1-line coefficient of E_{i+1}
constexpr int CLW_CP = 32;   // tridiagonal solve scratch
constexpr int CLW_R = 33;    // BiCGStab vectors
constexpr int CLW_R0 = 34;
constexpr int CLW_P = 35;
constexpr int CLW_V = 36;
constexpr int CLW_S = 37;
constexpr int CLW_Q = 38;
constexpr int CLW_Y = 39;
constexpr int CLW_Z = 40;
constexpr int CLW_N = 41;

//----------------------------------------------------------------------------------------
//! \struct CleanFvI
//! \brief the FV rays' intensities: cln_fvi, or cln_fvis in single precision
//! (clean_fv_float).  Reads give Real, Set stores in the precision of the array, Round
//! is a value as stored.

struct CleanFvI {
  DvceArray5D<Real> d;
  DvceArray5D<float> f;
  bool f32 = false;
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int n, const int k, const int j, const int i) const {
    return f32 ? static_cast<Real>(f(m,n,k,j,i)) : d(m,n,k,j,i);
  }
  KOKKOS_INLINE_FUNCTION
  Real Round(const Real v) const {
    return f32 ? static_cast<Real>(static_cast<float>(v)) : v;
  }
  KOKKOS_INLINE_FUNCTION
  void Set(const int m, const int n, const int k, const int j, const int i,
           const Real v) const {
    if (f32) {
      f(m,n,k,j,i) = static_cast<float>(v);
    } else {
      d(m,n,k,j,i) = v;
    }
  }
  int extent_int(const int r) const {return f32 ? f.extent_int(r) : d.extent_int(r);}
  // the intensities as they are held: with clean_fv_float, the doubles a restart read
  // until the next sweep converts them
  static CleanFvI Of(const DvceArray5D<Real> &d_, const DvceArray5D<float> &f_,
                     const bool f32_) {
    return CleanFvI{d_, f_, f32_ && d_.extent_int(1) == 0};
  }
};

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
  bool bracket;                // clean_consistency_bracket: with its bracket (else w only)
  bool posreal;                // clean_pos_real: kinetic y flux on low-order faces
  bool tilt;                   // clean_bracket_tilt: the E bracket in the factors g+-
  bool curvk;                  // clean_curv_iterate: the rays' curvature terms per E
  bool milne;                  // cell layout: sub-cell Eddington face model
  bool fvrays;                 // cell layout: FV rays (clean_rays = fv1 | fv)
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
constexpr int CLG_NV = 27;

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
//! \fn void CleanLineRange
//! \brief the lines [lo, hi) in steps st that team t of nt takes out of nline
//! (CleanCellLineRows): on the host the contiguous chunk [t nline/nt, (t+1) nline/nt), so
//! that a thread's lines lie next to each other in memory; on a GPU t, t + nt, ...

KOKKOS_INLINE_FUNCTION
void CleanLineRange(const int nline, const int t, const int nt, int &lo, int &hi,
                    int &st) {
  if constexpr (std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>) {
    lo = static_cast<int>(static_cast<int64_t>(nline)*t/nt);
    hi = static_cast<int>(static_cast<int64_t>(nline)*(t + 1)/nt);
    st = 1;
  } else {
    lo = t;
    hi = nline;
    st = nt;
  }
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
//! \fn void CleanFaceWeights
//! \brief coefficients of E_L and E_R in the flux through face (k,j,i) normal to dir, i.e.
//! everything but the off-diagonal Eddington terms and the old-vector flux G.  Physical
//! boundaries (non-periodic) take the boundary flux of their form.

KOKKOS_INLINE_FUNCTION
void CleanFaceWeights(const CleanArgs &x, const int dir, const int m, const int k,
                      const int j, const int i, Real &al, Real &ar) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  al = 0.0;
  ar = 0.0;
  if (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir) == M1_IBC_REFLECT) return;
    const bool beam = x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0;
    if (beam && x.marshak) return;
    ar = x.marshak ? -x.c*x.marshak_q : x.c*x.cs(m,CLN_HM+dir,k,j,i);
    return;
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir+1) == M1_IBC_REFLECT) return;
    al = x.marshak ? x.c*x.marshak_q : x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
    return;
  }
  const Real w = x.fc[dir](m,CLF_W,k,j,i);
  const Real th = x.fc[dir](m,CLF_TH,k,j,i);
  const Real g = (1.0 - w)*th*x.h*x.c*x.c*x.idx[dir];
  al = g*x.cs(m,CLN_D11+dir,kl,jl,il) + w*x.c*x.fc[dir](m,CLF_HP,k,j,i);
  ar = -g*x.cs(m,CLN_D11+dir,k,j,i) + w*x.c*x.fc[dir](m,CLF_HM,k,j,i);
  if (x.feedback) {
    const Real s = 0.5*(1.0 - w)*(1.0 - th);
    al += s*x.w(m,CLW_A1+dir,kl,jl,il);
    ar += s*x.w(m,CLW_A1+dir,k,j,i);
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
  if (normal) {
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
    div += 0.25*x.idx[b]*(de(k+dk,j+dj,i+di) + de(kl+dk,jl+dj,il+di)
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
  const Real w = x.fc[dir](m,CLF_W,k,j,i);
  const Real th = x.fc[dir](m,CLF_TH,k,j,i);
  return -(1.0 - w)*th*x.h*x.c*x.c*CleanFaceStress(x, e, dir, m, k, j, i, false);
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
    return x.marshak ? x.c*x.marshak_q*x.ebath[2*dir] : 0.25*x.c*x.ebath[2*dir];
  }
  if (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC) {
    if (CleanBc(x, m, 2*dir+1) == M1_IBC_REFLECT) return 0.0;
    if (CleanBc(x, m, 2*dir+1) == M1_IBC_FLUX) return dir ? x.ifx2[1](m,k,i) : x.iflux[1];
    return x.marshak ? -x.c*x.marshak_q*x.ebath[2*dir+1] : -0.25*x.c*x.ebath[2*dir+1];
  }
  const Real w = x.fc[dir](m,CLF_W,k,j,i);
  return (1.0 - w)*x.fc[dir](m,CLF_TH,k,j,i)*x.fc[dir](m,CLF_G,k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceFlux
//! \brief flux through face (k,j,i) normal to dir for E in work vector e; with_source adds
//! the E-independent part

KOKKOS_INLINE_FUNCTION
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
//! \fn Real CleanMomentumChange
//! \brief the momentum the radiation gives to cell (k,j,i) along dir over the stage, in
//! conservative form: at each face, minus the left side of the radiation flux equation,
//! -[(F - G)/c^2 + h div(D E)], half of each face to each neighbour.  E is work vector e.

KOKKOS_INLINE_FUNCTION
Real CleanMomentumChange(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  Real dm = 0.0;
  for (int side = 0; side < 2; ++side) {
    const int kf = k + side*dk, jf = j + side*dj, iff = i + side*di;
    const Real f = CleanFaceFlux(x, e, dir, m, kf, jf, iff, true);
    const Real g = x.fc[dir](m,CLF_G,kf,jf,iff);
    dm -= 0.5*((f - g)/(x.c*x.c) + x.h*CleanFaceStress(x, e, dir, m, kf, jf, iff, true));
  }
  return dm;
}

}  // namespace radm1
#endif  // RAD_M1_RAD_M1_CLEAN_HPP_
