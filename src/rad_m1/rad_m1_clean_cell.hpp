#ifndef RAD_M1_RAD_M1_CLEAN_CELL_HPP_
#define RAD_M1_RAD_M1_CLEAN_CELL_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_cell.hpp
//! \brief Device functions of the cell-centred layout of the clean implicit solve
//! (<rad_m1>/clean_layout = cell).
//!
//! The unknowns of every cell are E and y = F/c (all components in energy-density
//! units).  A stage of step h solves
//!   E (1 + h c kE) + h div fE = E_old + h c kP a T^4 - W
//!   y_b (1 + h c chi) - h chi a_b E + h div fY_b = y_old,b
//! with a_b = v_b + D_bc v_c (A = a E the advective flux) and W the work of the radiation
//! force on the gas.  At a face normal to d between cells L and R:
//!   fE   = c (y_d,L + y_d,R)/2 - al (s/2) (E_R - E_L) + RC
//!   fY_b = c (D_bd,L E_L + D_bd,R E_R)/2 - al (s/2) (y_b,R - y_b,L)
//!   RC   = (1 - al) th h c^2 [(g_L + g_R)/2 - ((D_dd E)_R - (D_dd E)_L)/dx_d]
//! with g_i = ((D_dd E)_{i+1} - (D_dd E)_{i-1})/(x_{i+1} - x_{i-1}), dx_d the distance
//! between the cell centres across the face, th = 1/(1 + h c chi_f),
//! s = c sqrt(max(D_dd,L, D_dd,R)) and al the dissipation weight.  RC is the momentum
//! interpolation of the flux to the face: where h c chi >> 1 it turns fE into the compact
//! diffusion flux.  With force_reference = clean (1 - al) th h c^2 is its steady value
//! (1 - al) c/chi_f (CleanDissipationWeight), so that a steady state of the scheme does
//! not depend on h.  The jumps
//! (E_R - E_L, y_R - y_L) are those of the face states: first order, or
//! (clean_recon = plm) the van Leer reconstruction of src/reconstruct/plm.hpp,
//! u_L^+ = u_L + dl dr/(dl + dr), whose weights dr/(dl + dr), dl/(dl + dr) are taken
//! from the Newton iterate (lagged) so that every pass stays a linear solve.
//! With clean_consistency every interior face flux also carries the fixed bracket
//! C = F_ray - f(Ur), the rays' kinetic flux less the face flux of their cell moments
//! Ur = (J, J h), so that Ur is a steady state of the discrete moment equations.  Its
//! streaming share w moves E and y with the rays' half-range flux factors,
//! c (h+_L U_L + h-_R U_R): the momentum deviation R = y - h E, which takes what the
//! closure's change leaves (y_old - h E), streams out as the rays do and is absorbed at
//! chi, so y relaxes to h E and the gas takes the exchange h chi (y - a E/c); on open
//! boundaries and walls R leaves with the inner cell's h+- (CleanCellDeviationFlux).
//!
//! On a spherical-polar mesh y holds the physical components (y_r, y_theta, y_phi), div
//! is the finite-volume divergence sum A_f f/V, with the angular-momentum conserving
//! weights of CleanCellRows, and the stress rows carry the curvature terms of
//! CleanCellSphSource.  E changes only by face fluxes, so the total energy is conserved
//! to round-off.  On the cubed sphere y holds (y_r, y_theta, y_phi) in the orthonormal
//! frame about the hydro's radius vector (CleanCsFrame).  The face fluxes are formed in
//! Cartesian components with the exact face vectors of CleanCsGeometry: d above is the
//! face normal n (y_d -> y.n, D_bd -> (D n)_b, D_dd -> n.D.n; D is Cartesian).  The rows
//! (CleanCellRows) are sum A_f f/V for E, its projections on r_hat and theta_hat, and the
//! torque about z for y_phi: E and L_z change only by face fluxes, a uniform isotropic
//! state is exact, and no curvature terms are needed.  The rays' face sums are Cartesian
//! and their kinetic flux takes the half ranges about the face normal (CleanCellKinetic).

#include <cmath>

#include "athena.hpp"
#include "rad_m1/rad_m1_clean.hpp"

namespace radm1 {

// face variables of the cell layout (arrays cln_fc1..3; CLF_TH: RC weight (1 - al) th)
constexpr int CLF_AL = 0;    // dissipation weight al
constexpr int CLF_AF = 2;    // dissipation speed s = c sqrt(max(D_dd,L, D_dd,R))
constexpr int CLF_LO = 3;    // 1: low-order face (CleanCellM1Jacobian)
constexpr int CLF_UG = 4;    // 1: physical boundary face of a marked cell, upwind stress
                             //    against the ghost (clean_pos_ghost)
constexpr int CLF_CB = CLF_N;  // consistency bracket of (fE, fY_1..3)
constexpr int CLF_WK = CLF_CB + 4;  // weight w of the rays' kinetic flux of U - Ur
constexpr int CLF_HK = CLF_WK + 1;  // cubed sphere: its h+ of L and h- of R about the
                                    // face normal (CleanCellConsistency)
constexpr int CLF_ML = CLF_HK + 2;  // clean_face_model = milne: mu and the sub-cell
                                    // terms of the derivatives (CleanCellMilneSetup)
constexpr int CLF_KK = CLF_ML + 9;  // cubed sphere, milne: k+ of L and k- of R about
                                    // the face normal (CleanCellConsistency)

// cell vectors (array cln_x), 4 components (E, y1, y2, y3) each
constexpr int CLC_O = 0;     // old vector
constexpr int CLC_K = 4;     // Newton iterate
constexpr int CLC_X = 8;     // solution
constexpr int CLC_B = 12;    // right-hand side
constexpr int CLC_R = 16;    // BiCGStab vectors
constexpr int CLC_R0 = 20;
constexpr int CLC_P = 24;
constexpr int CLC_V = 28;
constexpr int CLC_S = 32;
constexpr int CLC_Q = 36;
constexpr int CLC_Y = 40;
constexpr int CLC_Z = 44;
constexpr int CLC_AM = 48;   // preconditioner: 4x4 block of cell i-1 (row major)
constexpr int CLC_A0 = 64;   // diagonal block
constexpr int CLC_AP = 80;   // block of cell i+1
constexpr int CLC_CP = 96;   // block-tridiagonal solve scratch
constexpr int CLC_PSI = 112; // smoothness sensor of the old vector along x1, x2, x3
constexpr int CLC_PL = 115;  // lagged slope weights pl, pr of (dir, component):
constexpr int CLC_PR = 127;  //   slot + 4 dir + component
constexpr int CLC_JM = 139;  // preconditioner blocks of cells j-1, j+1, k-1, k+1
constexpr int CLC_JP = 155;
constexpr int CLC_KM = 171;
constexpr int CLC_KP = 187;
constexpr int CLC_T1 = 203;  // preconditioner scratch vectors
constexpr int CLC_T2 = 207;
constexpr int CLC_U = 211;   // the rays' cell moments (J, J h)
constexpr int CLC_MB = 215;  // B = kP a T^4/kE of the old vector (clean_face_model)
constexpr int CLC_N = 216;

// dissipation weight candidates (<rad_m1>/clean_alpha)
constexpr int CLA_TAU = 0;          // A(tau)^2, A = (1 - exp(-tau))/tau
constexpr int CLA_TAU_SENSOR = 1;   // max(A(tau)^2 psi, s_ray)
constexpr int CLA_TAU_RAY = 2;      // max(A(tau)^2, s_ray)
constexpr int CLA_ONE = 3;          // 1
constexpr int CLA_ZERO = 4;         // 0

//----------------------------------------------------------------------------------------
//! \fn int CleanLineSlotM
//! \brief slot of the preconditioner block of the lower neighbour along dir

KOKKOS_INLINE_FUNCTION
int CleanLineSlotM(const int dir) {
  return (dir == 0) ? CLC_AM : ((dir == 1) ? CLC_JM : CLC_KM);
}

//----------------------------------------------------------------------------------------
//! \fn int CleanLineSlotP
//! \brief slot of the preconditioner block of the upper neighbour along dir

KOKKOS_INLINE_FUNCTION
int CleanLineSlotP(const int dir) {
  return (dir == 0) ? CLC_AP : ((dir == 1) ? CLC_JP : CLC_KP);
}

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
                            const Real sray, Real &rc) {
  // A = 1 - b tau
  const Real b = (tau > 1.0e-4) ? (1.0 + expm1(-tau)/tau)/tau
                                : 0.5 - tau*(1.0 - 0.25*tau)/6.0;
  const Real a = (tau > 1.0e-4) ? -expm1(-tau)/tau : 1.0 - b*tau;
  Real al = 0.0;
  switch (mode) {
    case CLA_TAU: al = a*a; break;
    case CLA_TAU_SENSOR: al = fmax(a*a*psi, sray); break;
    case CLA_TAU_RAY: al = fmax(a*a, sray); break;
    case CLA_ONE: al = 1.0; break;
    default: break;
  }
  rc = (mode == CLA_TAU) ? b*(1.0 + a) : (1.0 - al)/fmax(tau, 2.0);
  return al;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDiag
//! \brief component D_ab of the Eddington tensor of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanDiag(const CleanArgs &x, const int a, const int b, const int m, const int k,
               const int j, const int i) {
  const int dab[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
  return x.cs(m,CLN_D11+dab[a][b],k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanNf
//! \brief component b of the unit normal of face (kf,jf,iff) normal to dir: e_dir, on the
//! cubed sphere the face's own normal

KOKKOS_INLINE_FUNCTION
Real CleanNf(const CleanArgs &x, const int b, const int dir, const int m, const int kf,
             const int jf, const int iff) {
  return x.cub ? x.g(m,CLG_N+3*dir+b,kf,jf,iff) : ((b == dir) ? 1.0 : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDn
//! \brief component b of D n of cell (k,j,i), n the unit normal of face (kf,jf,iff)
//! normal to dir (D_b,dir off the cubed sphere)

KOKKOS_INLINE_FUNCTION
Real CleanDn(const CleanArgs &x, const int b, const int dir, const int m, const int kf,
             const int jf, const int iff, const int k, const int j, const int i) {
  if (!x.cub) return CleanDiag(x, b, dir, m, k, j, i);
  Real d = 0.0;
  for (int c = 0; c < 3; ++c) {
    d += CleanDiag(x, b, c, m, k, j, i)*CleanNf(x, c, dir, m, kf, jf, iff);
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDnn
//! \brief n.D.n of cell (k,j,i) with the normal of face (kf,jf,iff) normal to dir

KOKKOS_INLINE_FUNCTION
Real CleanDnn(const CleanArgs &x, const int dir, const int m, const int kf, const int jf,
              const int iff, const int k, const int j, const int i) {
  if (!x.cub) return CleanDiag(x, dir, dir, m, k, j, i);
  Real d = 0.0;
  for (int b = 0; b < 3; ++b) {
    d += CleanNf(x, b, dir, m, kf, jf, iff)*CleanDn(x, b, dir, m, kf, jf, iff, k, j, i);
  }
  return d;
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
  for (int c = 0; c < 3; ++c) q[0][c] = x.g(m,CLG_R+c,k,j,i);
  const Real rho = std::sqrt(q[0][0]*q[0][0] + q[0][1]*q[0][1]);
  q[2][0] = -q[0][1]/rho;
  q[2][1] = q[0][0]/rho;
  q[2][2] = 0.0;
  q[1][0] = q[2][1]*q[0][2];
  q[1][1] = -q[2][0]*q[0][2];
  q[1][2] = q[2][0]*q[0][1] - q[2][1]*q[0][0];
  return x.x1v(m,i)*rho;
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
      nd[d][c] = (d == 0) ? x.g(m,CLG_R+c,k,j,i) :
                 x.g(m,CLG_N+3*d+c,k,j,i) + x.g(m,CLG_N+3*d+c,k+dk,j+dj,i);
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
//! \fn void CleanCellState
//! \brief E and y of cell (k,j,i) at slot v of x.u, y in Cartesian components on the
//! cubed sphere (the face fluxes are formed in them)

KOKKOS_INLINE_FUNCTION
void CleanCellState(const CleanArgs &x, const int v, const int m, const int k,
                    const int j, const int i, Real s[4]) {
  for (int c = 0; c < 4; ++c) s[c] = x.u(m,v+c,k,j,i);
  if (!x.cub) return;
  Real q[3][3];
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) {
    s[1+b] = x.u(m,v+1,k,j,i)*q[0][b] + x.u(m,v+2,k,j,i)*q[1][b] +
             x.u(m,v+3,k,j,i)*q[2][b];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellRows
//! \brief the contribution g of the flux f (fE and the stress column) through face
//! (kf,jf,iff) normal to dir to h div of the rows (E, y_1..3) of cell (k,j,i).
//! Cartesian: h f/dx_d.  Spherical polar: h A_f/V f, times r_f/r on the r faces of the
//! y_theta and y_phi rows and sin(theta_f)/sin(theta) on the theta faces of the y_phi
//! row: the angular-momentum conserving forms r^-3 d_r(r^3 P_ra) and
//! sin^-2 d_theta(sin^2 P_tp), which contain the curvature terms P_rt/r, P_rp/r and
//! cot(theta) P_tp/r.  Cubed sphere (f in Cartesian components): h/V times A fE,
//! A r_hat.f, A theta_hat.f and the torque (M x f)_z/R about z with M = int x dA of the
//! face (r_f A n on a radial face, no torque of an isotropic stress): the y_phi row is
//! the conservation form of L_z = V R y_phi, the y_r and y_theta rows the projected
//! divergence (a uniform isotropic stress cancels since sum_f A_f = 0).

KOKKOS_INLINE_FUNCTION
void CleanCellRows(const CleanArgs &x, const int dir, const int m, const int kf,
                   const int jf, const int iff, const int k, const int j, const int i,
                   const Real f[4], Real g[4]) {
  if (!x.sph && !x.cub) {
    for (int r = 0; r < 4; ++r) g[r] = x.h*x.idx[dir]*f[r];
    return;
  }
  const Real hv = x.h/x.vol(m,k,j,i), a = x.area[dir](m,kf,jf,iff);
  if (x.sph) {
    Real w[4] = {hv*a, hv*a, hv*a, hv*a};
    if (dir == 0) {
      const Real rf = x.x1f(m,iff)/x.x1v(m,i);
      w[2] *= rf;
      w[3] *= rf;
    } else if (dir == 1) {
      w[3] *= sin(x.x2f(m,jf))/sin(x.x2v(m,j));
    }
    for (int r = 0; r < 4; ++r) g[r] = w[r]*f[r];
    return;
  }
  Real q[3][3], mf[3];
  const Real rc = CleanCsFrame(x, m, k, j, i, q);
  for (int c = 0; c < 3; ++c) {
    mf[c] = (dir == 0) ? x.x1f(m,iff)*a*x.g(m,CLG_N+c,kf,jf,iff)
                       : x.g(m,CLG_M+3*(dir-1)+c,kf,jf,iff);
  }
  g[0] = hv*a*f[0];
  for (int b = 0; b < 2; ++b) {
    g[1+b] = hv*a*(q[b][0]*f[1] + q[b][1]*f[2] + q[b][2]*f[3]);
  }
  g[3] = hv*(mf[0]*f[2] - mf[1]*f[1])/rc;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellBlock
//! \brief the block out of the rows of cell (k,j,i) (CleanCellRows) in the states (E, y
//! in its frame on the cubed sphere) of cell (kn,jn,in), from the Jacobian jj of the face
//! flux through (kf,jf,iff) normal to dir

KOKKOS_INLINE_FUNCTION
void CleanCellBlock(const CleanArgs &x, const int dir, const int m, const int kf,
                    const int jf, const int iff, const int k, const int j, const int i,
                    const int kn, const int jn, const int in, const Real jj[4][4],
                    Real out[4][4]) {
  Real q[3][3];
  if (x.cub) CleanCsFrame(x, m, kn, jn, in, q);
  for (int c = 0; c < 4; ++c) {
    Real f[4], g[4];
    for (int r = 0; r < 4; ++r) {
      f[r] = (c == 0 || !x.cub) ? jj[r][c]
             : jj[r][1]*q[c-1][0] + jj[r][2]*q[c-1][1] + jj[r][3]*q[c-1][2];
    }
    CleanCellRows(x, dir, m, kf, jf, iff, k, j, i, f, g);
    for (int r = 0; r < 4; ++r) out[r][c] = g[r];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanGasMomentum
//! \brief the gas momentum of cell (k,j,i) in the frame of y.  On the cubed sphere the
//! hydro stores m_r = p.r_hat and the COVARIANT components m_a = p.e_a on its unit panel
//! tangents (non-orthogonal, g_ab = e_a.e_b, hydro_am.hpp): p_t = g^ab m_b e_a.

KOKKOS_INLINE_FUNCTION
void CleanGasMomentum(const CleanArgs &x, const DvceArray5D<Real> &uh, const int m,
                      const int k, const int j, const int i, Real p[3]) {
  for (int n = 0; n < 3; ++n) p[n] = uh(m,IM1+n,k,j,i);
  if (!x.cub) return;
  Real q[3][3], pt[3], c = 0.0;
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) c += x.g(m,CLG_E+b,k,j,i)*x.g(m,CLG_E+3+b,k,j,i);
  const Real u1 = (p[1] - c*p[2])/(1.0 - c*c), u2 = (p[2] - c*p[1])/(1.0 - c*c);
  for (int b = 0; b < 3; ++b) {
    pt[b] = u1*x.g(m,CLG_E+b,k,j,i) + u2*x.g(m,CLG_E+3+b,k,j,i);
  }
  for (int a = 1; a < 3; ++a) p[a] = q[a][0]*pt[0] + q[a][1]*pt[1] + q[a][2]*pt[2];
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
    dh[a] = dt[0]*x.g(m,CLG_E+3*(a-1),k,j,i) + dt[1]*x.g(m,CLG_E+3*(a-1)+1,k,j,i) +
            dt[2]*x.g(m,CLG_E+3*(a-1)+2,k,j,i);
  }
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
//! \fn Real CleanRefForce
//! \brief component b of the reference force density that force_reference withholds
//! from the gas: the stored force of the first coupled solve (clean), d arad_ref along x1
//! (wb_arad), else 0.  The clean one is in the restart file with the hesdirk2 slopes
//! (Time2RstPack).

KOKKOS_INLINE_FUNCTION
Real CleanRefForce(const int fref, const DvceArray5D<Real> &fr,
                   const DvceArray4D<Real> &aref, const Real d, const int b, const int m,
                   const int k, const int j, const int i) {
  if (fref == M1_FREF_CLEAN) return fr(m,b,k,j,i);
  return (fref == M1_FREF_WB_ARAD && b == 0) ? d*aref(m,k,j,i) : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellDrag
//! \brief with gas coupling: g_bc = h chi E_k (delta_bc + D_bc)/(c rho), the change of
//! the advective source -h chi a_b E of row y_b per unit change of y_c + h div fY_c, the
//! gas momentum the cell gives up (dv = -(dy + h div fY)/(c rho)): the gas velocity taken
//! implicit in the solve (radiation force and drag), so that the velocity iteration
//! contracts where the gas responds within h (rows and columns in the frame of y);
//! zero while the gas momentum is frozen (steady start)

KOKKOS_INLINE_FUNCTION
void CleanCellDrag(const CleanArgs &x, const Real chi, const Real rho, const int m,
                   const int k, const int j, const int i, Real g[3][3]) {
  const Real s = x.vzero ? 0.0 : x.h*chi*x.u(m,CLC_K,k,j,i)/(x.c*rho);
  Real q[3][3];
  if (x.cub) CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) {
    for (int c = 0; c < 3; ++c) {
      Real d = 0.0;
      if (x.cub) {   // D in the frame of y: Q D Q^T
        for (int e = 0; e < 3; ++e) {
          for (int f = 0; f < 3; ++f) d += q[b][e]*CleanDiag(x, e, f, m, k, j, i)*q[c][f];
        }
      } else {
        d = CleanDiag(x, b, c, m, k, j, i);
      }
      g[b][c] = s*((b == c ? 1.0 : 0.0) + d);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanJumpCoefficients
//! \brief coefficients cj of u_{L-1}, u_L, u_R, u_{R+1} in the jump u_R^- - u_L^+ of
//! component c at face (k,j,i) normal to dir: (0, -1, 1, 0) at first order, with the
//! lagged half-slopes s_i = (pl_i (u_i - u_{i-1}) + pr_i (u_{i+1} - u_i))/2 otherwise

KOKKOS_INLINE_FUNCTION
void CleanJumpCoefficients(const CleanArgs &x, const int dir, const int c, const int m,
                           const int k, const int j, const int i, Real cj[4]) {
  cj[0] = 0.0;
  cj[1] = -1.0;
  cj[2] = 1.0;
  cj[3] = 0.0;
  if (!x.recon) return;
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  const int n = 4*dir + c;
  const Real pll = x.u(m,CLC_PL+n,kl,jl,il), prl = x.u(m,CLC_PR+n,kl,jl,il);
  const Real plr = x.u(m,CLC_PL+n,k,j,i), prr = x.u(m,CLC_PR+n,k,j,i);
  cj[0] = 0.5*pll;
  cj[1] += 0.5*(plr - pll + prl);
  cj[2] += 0.5*(prr - plr - prl);
  cj[3] = -0.5*prr;
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanFaceDissipation
//! \brief dissipation weight al and al s/2 of the interior face (k,j,i) normal to dir;
//! a low-order face takes al = 1, s = c (its fY).  Returns true for a low-order face.

KOKKOS_INLINE_FUNCTION
bool CleanFaceDissipation(const CleanArgs &x, const int dir, const int m, const int k,
                          const int j, const int i, Real &al, Real &hs) {
  const bool lo = x.fc[dir](m,CLF_LO,k,j,i) > 0.0;
  al = lo ? 1.0 : x.fc[dir](m,CLF_AL,k,j,i);
  hs = 0.5*al*(lo ? x.c : x.fc[dir](m,CLF_AF,k,j,i));
  return lo;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellKinetic
//! \brief the rays' half-range factors of the kinetic flux through the interior face
//! (k,j,i) normal to dir: h+ of L and h- of R of the face family; on the cubed sphere
//! about the face's own normal (CLF_HK)

KOKKOS_INLINE_FUNCTION
void CleanCellKinetic(const CleanArgs &x, const int dir, const int m, const int k,
                      const int j, const int i, Real &hp, Real &hm) {
  const int kl = k - (dir == 2), jl = j - (dir == 1), il = i - (dir == 0);
  hp = x.cub ? x.fc[dir](m,CLF_HK,k,j,i) : x.cs(m,CLN_HP+dir,kl,jl,il);
  hm = x.cub ? x.fc[dir](m,CLF_HK+1,k,j,i) : x.cs(m,CLN_HM+dir,k,j,i);
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
//! \fn void CleanCellRhieChow
//! \brief coefficients q of (D_dd E) in the cells L-1, L, R, R+1 of the RC term of the
//! interior face (k,j,i) normal to dir (none at a low-order face).  The cell gradients
//! g_L, g_R are central over the centre distances, the face gradient compact.

KOKKOS_INLINE_FUNCTION
void CleanCellRhieChow(const CleanArgs &x, const int dir, const bool lo, const int m,
                       const int k, const int j, const int i, Real q[4]) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const Real a = lo ? 0.0 : x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c;
  const Real d0 = CleanCellFaceDistance(x, dir, m, k, j, i);
  const Real dm = CleanCellFaceDistance(x, dir, m, k-dk, j-dj, i-di);
  const Real dp = CleanCellFaceDistance(x, dir, m, k+dk, j+dj, i+di);
  const Real gl = 0.5*a/(dm + d0), gr = 0.5*a/(d0 + dp);
  q[0] = -gl;
  q[1] = a/d0 - gr;
  q[2] = gl - a/d0;
  q[3] = gr;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellSphSource
//! \brief spherical polar: h times the remaining curvature terms of the stress divergence
//! of cell (k,j,i) per unit E of the cell, g_b for the rows y_b: -(P_tt + P_pp)/r in y_r
//! and -cot(theta) P_pp/r in y_theta, in the face-area form -c (D_tt + D_pp)
//! (A_r+ - A_r-)/(2V) and -c D_pp (A_t+ - A_t-)/V that cancels the area terms of a
//! uniform isotropic stress exactly.  With clean_consistency the streaming share w of
//! the faces that carry the rays' flux F_ray (low-order or not) moves only y beyond the
//! rays' moments: its curvature terms are the rays' own, CLN_G (kin = true, independent
//! of E, shared by the r faces' w as the area terms), so that the rays' moments (J, J h)
//! are its steady state; the y_phi row keeps its angular-momentum flux form (no source),
//! so that L_z stays exact.

KOKKOS_INLINE_FUNCTION
void CleanCellSphSource(const CleanArgs &x, const int m, const int k, const int j,
                        const int i, Real g[3], const bool kin = false) {
  const Real hv = x.h*x.c/x.vol(m,k,j,i);
  // w of the faces that carry the rays' flux (interior and open x1 faces)
  auto wf = [&](const int dir, const int kk, const int jj, const int ii) {
    const bool on = !CleanBoundaryFace(x, dir, m, kk, jj, ii) ||
                    (dir == 0 && CleanBc(x, m, (ii == x.is) ? 0 : 1) != M1_IBC_REFLECT);
    return (x.consist && on) ? x.fc[dir](m,CLF_WK,kk,jj,ii) : 0.0;
  };
  const Real am = x.area[0](m,k,j,i), ap = x.area[0](m,k,j,i+1);
  if (kin) {
    const Real w = (wf(0, k, j, i+1)*ap - wf(0, k, j, i)*am)/(ap - am);
    for (int b = 0; b < 3; ++b) g[b] = (b < 2) ? w*hv*x.cs(m,CLN_G+b,k,j,i) : 0.0;
    return;
  }
  const Real dtt = CleanDiag(x, 1, 1, m, k, j, i), dpp = CleanDiag(x, 2, 2, m, k, j, i);
  auto d = [&](const int dir, const int kk, const int jj, const int ii, const Real dc) {
    return (1.0 - wf(dir, kk, jj, ii))*dc*x.area[dir](m,kk,jj,ii);
  };
  const Real dr = 0.5*(dtt + dpp);
  g[0] = -hv*(d(0, k, j, i+1, dr) - d(0, k, j, i, dr));
  g[1] = (x.nf > 1) ? -hv*(d(1, k, j+1, i, dpp) - d(1, k, j, i, dpp)) : 0.0;
  g[2] = 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellGhost
//! \brief incoming state outside the open physical boundary face (k,j,i) normal to dir:
//! the beam at x1min (E_g = F_b/(c sqrt(D_11)), y_g,b = D_1b F_b/(c D_11) with D of the
//! inner cell, exact for a collimated beam), else the isotropic bath (D = 1/3, y = 0).
//! Returns E_g, y_g and the stress column D_g n E_g (n the face normal).

KOKKOS_INLINE_FUNCTION
void CleanCellGhost(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i, const bool low, const Real n[3], Real &eg,
                    Real yg[3], Real pg[3]) {
  for (int b = 0; b < 3; ++b) {
    yg[b] = 0.0;
    pg[b] = 0.0;
  }
  if (low && dir == 0 && x.use_beam && x.beam(m,k-x.ks,j-x.js) > 0.0) {
    const Real yb = x.beam(m,k-x.ks,j-x.js)/x.c;
    const Real d11 = fmax(CleanDiag(x, 0, 0, m, k, j, i), 1.0e-12);
    eg = yb/std::sqrt(d11);
    for (int b = 0; b < x.nf; ++b) {
      yg[b] = CleanDiag(x, 0, b, m, k, j, i)*yb/d11;
      pg[b] = CleanDiag(x, b, 0, m, k, j, i)*eg;
    }
    return;
  }
  eg = x.ebath[2*dir + (low ? 0 : 1)];
  for (int b = 0; b < 3; ++b) pg[b] = n[b]*eg/3.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCellOutStress
//! \brief stress coefficient of the radiation leaving through the open boundary face
//! (kf,jf,iff) of inner cell (k,j,i) normal to dir, column b: the rays' half-range
//! k+^dir_b = K+_b,dir/J at a high face, (D n)_b - k+^dir_b at a low face

KOKKOS_INLINE_FUNCTION
Real CleanCellOutStress(const CleanArgs &x, const int dir, const int b, const int m,
                        const int kf, const int jf, const int iff, const int k,
                        const int j, const int i, const bool low) {
  const Real kp = x.cs(m,CLN_KP+3*dir+b,k,j,i);
  return low ? CleanDn(x, b, dir, m, kf, jf, iff, k, j, i) - kp : kp;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellM1Source
//! \brief the E-independent part of the moment scheme's face fluxes: the consistency
//! bracket and, on a physical boundary, the incoming E flux of CleanFaceSource and the
//! stress of the incoming radiation (beam: its stress column; isotropic bath: c E_bath/6
//! on the normal row); in a first-order solve or on a low-order face the ghost part of
//! the upwind stress

KOKKOS_INLINE_FUNCTION
void CleanCellM1Source(const CleanArgs &x, const int dir, const int m, const int k,
                       const int j, const int i, Real s[4]) {
  for (int r = 0; r < 4; ++r) s[r] = x.consist ? x.fc[dir](m,CLF_CB+r,k,j,i) : 0.0;
  if (!CleanBoundaryFace(x, dir, m, k, j, i)) return;
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const bool low = (pos == lo);
  if (CleanBc(x, m, 2*dir + (low ? 0 : 1)) == M1_IBC_REFLECT) return;
  Real n[3];
  for (int b = 0; b < 3; ++b) n[b] = CleanNf(x, b, dir, m, k, j, i);
  s[0] += CleanFaceSource(x, dir, m, k, j, i);
  if (x.all_low || x.fc[dir](m,CLF_LO,k,j,i) > 0.0 || x.fc[dir](m,CLF_UG,k,j,i) > 0.0) {
    // ghost part of the upwind stress of a first-order solve or a marked boundary face
    const int kc = low ? k : k - (dir == 2), jc = low ? j : j - (dir == 1);
    const int ic = low ? i : i - (dir == 0);
    Real eg, yg[3], pg[3];
    CleanCellGhost(x, dir, m, kc, jc, ic, low, n, eg, yg, pg);
    const Real sg = low ? 0.5*x.c : -0.5*x.c;
    for (int b = 0; b < x.nf; ++b) s[1+b] += 0.5*x.c*pg[b] + sg*yg[b];
    return;
  }
  if (low && dir == 0 && x.use_beam && x.beam(m,k-x.ks,j-x.js) > 0.0) {
    for (int b = 0; b < x.nf; ++b) s[1+b] += x.bstress(m,b,k-x.ks,j-x.js);
    return;
  }
  for (int b = 0; b < x.nf; ++b) s[1+b] += n[b]*x.c*x.ebath[2*dir + (low ? 0 : 1)]/6.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellDeviationFlux
//! \brief adds the flux wh P R of the momentum deviation R = y - h E of cell (k,j,i)
//! through the boundary face (kf,jf,iff) normal to dir to the flux rows of jj (wh = w c
//! h+-, the rays' half-range flux factor): P = I on an open face, n n^T at a wall (its
//! mirror image returns R_t and reverses R_n)

KOKKOS_INLINE_FUNCTION
void CleanCellDeviationFlux(const CleanArgs &x, const int dir, const int m, const int kf,
                            const int jf, const int iff, const int k, const int j,
                            const int i, const bool wall, const Real wh, Real jj[4][4]) {
  Real s[3], h[3];
  for (int b = 0; b < 3; ++b) s[b] = x.cs(m,CLN_HP+b,k,j,i) + x.cs(m,CLN_HM+b,k,j,i);
  CleanRayFlux(x, m, k, j, i, s, h);
  for (int b = 0; b < x.nf; ++b) {
    for (int c = 0; c < x.nf; ++c) {
      const Real p = wall ? CleanNf(x, b, dir, m, kf, jf, iff)*
                            CleanNf(x, c, dir, m, kf, jf, iff) : ((b == c) ? 1.0 : 0.0);
      jj[1+b][1+c] += wh*p;
      jj[1+b][0] -= wh*p*h[c];
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellM1Jacobian
//! \brief derivatives jl, jr of the moment fluxes (fE, fY_1..3) at face (k,j,i) normal to
//! dir with respect to the states (E, y_1..3) of the cells L and R.  The terms in the
//! states of L-1 and R+1 (RC, reconstruction) are left to CleanCellFaceFlux.  On
//! physical boundaries fE takes the boundary form of CleanFaceWeights and fY the stress
//! of the radiation leaving through the face, c k_out E of the inner cell
//! (CleanCellOutStress; the incoming part is in CleanCellFaceSource).  In a first-order
//! solve and on a low-order boundary face fY is the upwind (a = c) flux between the inner
//! cell and the ghost of CleanCellGhost; a reflecting wall carries only the normal stress
//! of the inner cell.  A low-order interior face takes its E flux from the rays' half
//! ranges and the compact diffusion flux, blended with the face's weight a,
//!   fE = a c (h+_L E_L + h-_R E_R) - (1 - a) th h c^2 ((D_dd E)_R - (D_dd E)_L)/dx_d,
//! free of y and monotone in E (as the boundary form), so that the E row of a cell with
//! every face low-order is an M-matrix (the work W aside); its fY has the cells' own
//! stress with first-order dissipation at speed c.  With clean_consistency the leaving
//! stress and the wall carry the momentum deviation of the inner cell out
//! (CleanCellDeviationFlux, weight w).

KOKKOS_INLINE_FUNCTION
void CleanCellM1Jacobian(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i, Real jl[4][4], Real jr[4][4]) {
  const int kl = k - (dir == 2), jl_ = j - (dir == 1), il = i - (dir == 0);
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      jl[r][c] = 0.0;
      jr[r][c] = 0.0;
    }
  }
  const bool blo = (pos == lo && CleanBc(x, m, 2*dir) != M1_IBC_PERIODIC);
  const bool bhi = (pos == hi && CleanBc(x, m, 2*dir+1) != M1_IBC_PERIODIC);
  if (blo || bhi) {
    const int bc = blo ? CleanBc(x, m, 2*dir) : CleanBc(x, m, 2*dir+1);
    Real (*jj)[4] = blo ? jr : jl;
    const int kc = blo ? k : kl, jc = blo ? j : jl_, ic = blo ? i : il;
    const Real wo = x.consist ? x.c*x.fc[dir](m,CLF_WK,k,j,i)*
                                x.cs(m,(blo ? CLN_HM : CLN_HP)+dir,kc,jc,ic) : 0.0;
    if (bc == M1_IBC_REFLECT) {
      const Real dnn = CleanDnn(x, dir, m, k, j, i, kc, jc, ic);
      for (int b = 0; b < x.nf; ++b) jj[1+b][0] = x.c*CleanNf(x, b, dir, m, k, j, i)*dnn;
      CleanCellDeviationFlux(x, dir, m, k, j, i, kc, jc, ic, true, 2.0*wo, jj);
      return;
    }
    const Real sg = blo ? -0.5*x.c : 0.5*x.c;
    Real al, ar;
    CleanFaceWeights(x, dir, m, k, j, i, al, ar);
    jj[0][0] = blo ? ar : al;
    const bool up = x.all_low || x.fc[dir](m,CLF_LO,k,j,i) > 0.0 ||
                    x.fc[dir](m,CLF_UG,k,j,i) > 0.0;
    for (int b = 0; b < x.nf; ++b) {
      if (up) {
        // first-order solve or marked boundary face: upwind (a = c) stress between the
        // inner cell and the ghost
        jj[1+b][0] = 0.5*x.c*CleanDn(x, b, dir, m, k, j, i, kc, jc, ic);
        jj[1+b][1+b] = sg;
      } else {
        jj[1+b][0] = x.c*CleanCellOutStress(x, dir, b, m, k, j, i, kc, jc, ic, blo);
      }
    }
    if (!up) CleanCellDeviationFlux(x, dir, m, k, j, i, kc, jc, ic, false, wo, jj);
    if (up || x.consist) return;
    // without the consistency term the moment system's outgoing wave (speed c sqrt(D_nn),
    // eigenvector (1, -+sqrt(D_nn) n)) carries the deviation R_n = (y - h E).n from the
    // rays' flux out, (c R_n/2, -+c sqrt(D_nn) R_n n/2): the boundary prescribes only the
    // incoming wave and a consistent state (R = 0) keeps the kinetic flux
    Real sh[3], h[3], n[3], hn = 0.0;
    for (int b = 0; b < 3; ++b) {
      sh[b] = x.cs(m,CLN_HP+b,kc,jc,ic) + x.cs(m,CLN_HM+b,kc,jc,ic);
      n[b] = CleanNf(x, b, dir, m, k, j, i);
    }
    CleanRayFlux(x, m, kc, jc, ic, sh, h);
    for (int b = 0; b < 3; ++b) hn += h[b]*n[b];
    const Real sd = (blo ? -0.5 : 0.5)*x.c*std::sqrt(fmax(CleanDnn(x, dir, m, k, j, i, kc,
                                                                    jc, ic), 0.0));
    jj[0][0] -= 0.5*x.c*hn;
    for (int b = 0; b < x.nf; ++b) {
      jj[0][1+b] += 0.5*x.c*n[b];
      jj[1+b][0] -= sd*n[b]*hn;
      for (int c = 0; c < x.nf; ++c) jj[1+b][1+c] += sd*n[b]*n[c];
    }
    return;
  }
  Real al, hs;
  const bool low = CleanFaceDissipation(x, dir, m, k, j, i, al, hs);
  Real q[4];
  CleanCellRhieChow(x, dir, low, m, k, j, i, q);
  Real cj[4] = {0.0, -1.0, 1.0, 0.0};
  if (!low) CleanJumpCoefficients(x, dir, 0, m, k, j, i, cj);
  jl[0][0] = -hs*cj[1] + q[1]*CleanDnn(x, dir, m, k, j, i, kl, jl_, il);
  jr[0][0] = -hs*cj[2] + q[2]*CleanDnn(x, dir, m, k, j, i, k, j, i);
  for (int b = 0; b < x.nf; ++b) {
    const Real nb = CleanNf(x, b, dir, m, k, j, i);
    jl[0][1+b] = 0.5*x.c*nb;
    jr[0][1+b] = 0.5*x.c*nb;
    if (!low) CleanJumpCoefficients(x, dir, 1+b, m, k, j, i, cj);
    jl[1+b][0] = 0.5*x.c*CleanDn(x, b, dir, m, k, j, i, kl, jl_, il);
    jr[1+b][0] = 0.5*x.c*CleanDn(x, b, dir, m, k, j, i, k, j, i);
    jl[1+b][1+b] = -hs*cj[1];
    jr[1+b][1+b] = -hs*cj[2];
  }
  if (low) {
    const Real a = x.fc[dir](m,CLF_AL,k,j,i);
    const Real g = x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c/
                   CleanCellFaceDistance(x, dir, m, k, j, i);
    Real hp, hm;
    CleanCellKinetic(x, dir, m, k, j, i, hp, hm);
    jl[0][0] = a*x.c*hp + g*CleanDnn(x, dir, m, k, j, i, kl, jl_, il);
    jr[0][0] = a*x.c*hm - g*CleanDnn(x, dir, m, k, j, i, k, j, i);
    for (int b = 0; b < 3; ++b) jl[0][1+b] = jr[0][1+b] = 0.0;
  }
  if (!x.consist) return;
  // clean_consistency: w of the face takes the rays' kinetic upwind flux, the cells'
  // half-range factors on E and y (c h+_L U_L + c h-_R U_R)
  const Real w = x.fc[dir](m,CLF_WK,k,j,i);
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      jl[r][c] *= 1.0 - w;
      jr[r][c] *= 1.0 - w;
    }
  }
  Real hp, hm;
  CleanCellKinetic(x, dir, m, k, j, i, hp, hm);
  for (int r = 0; r <= x.nf; ++r) {
    jl[r][r] += w*x.c*hp;
    jr[r][r] += w*x.c*hm;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCellMilneWeight
//! \brief clean_face_model = milne: weight mu = (1 - exp(-2 tau_h))^4 of the sub-cell
//! face fluxes at face (k,j,i) normal to dir, tau_h = chi d/2 the larger optical depth
//! from a cell centre to the face (2 E_3(tau_h) ~ exp(-2 tau_h) is the share of the
//! outgoing half-range flux that crosses the half cell).  The sub-cell flux, which does
//! not carry y, differs from the moment scheme's by O(1) in thin cells, so
//! mu = O(tau_h^4) leaves smooth solutions as they are; below 1e-12 it is zero.
//! Low-order faces, a first-order solve and reflecting boundaries keep the moment
//! scheme.

KOKKOS_INLINE_FUNCTION
Real CleanCellMilneWeight(const CleanArgs &x, const int dir, const int m, const int k,
                          const int j, const int i) {
  if (!x.milne || x.all_low || x.fc[dir](m,CLF_LO,k,j,i) > 0.0 ||
      x.fc[dir](m,CLF_UG,k,j,i) > 0.0) return 0.0;
  const bool bdry = CleanBoundaryFace(x, dir, m, k, j, i);
  const bool blo = bdry && (dir == 0 ? i == x.is : (dir == 1 ? j == x.js : k == x.ks));
  if (bdry && CleanBc(x, m, 2*dir + (blo ? 0 : 1)) == M1_IBC_REFLECT) return 0.0;
  const Real d = CleanCellFaceDistance(x, dir, m, k, j, i);
  const Real tl = blo ? 0.0 : x.op(m,M1_OP_T,k-(dir == 2),j-(dir == 1),i-(dir == 0));
  const Real tr = (bdry && !blo) ? 0.0 : x.op(m,M1_OP_T,k,j,i);
  const Real mu = SQR(SQR(expm1(-fmax(tl, tr)*d)));
  return (mu > 1.0e-12) ? mu : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellMilne
//! \brief clean_face_model = milne: the sub-cell Eddington fluxes of face (k,j,i) normal
//! to dir, blended with weight mu into the moment scheme's derivatives jl, jr (s null) or
//! into its E-independent part s (jl, jr then hold the moment scheme's derivatives).
//! Between its centre and the face each cell solves E'' = k^2 (E - B), y = -E'/3 in the
//! optical depth tau (k = sqrt(3 eps), eps = kE/chi, B = kP a T^4/kE linear in tau with
//! the limited slope of the cell centres), so that the outward flux at the face is
//!   y = (e - E_f)/r,  e = c E_c + B_f - c B_c - (B_f - B_c) t,  r = 3 tau_h t,
//! z = k tau_h, c = 1/cosh z, t = tanh(z)/z, linear in the face value E_f.  The kinetic
//! interface y = A_L E_fL - A_R E_fR + (th_L + th_R) y/2 couples the sides: a side emits
//! the share th = 1 - exp(-2 tau_h) of the Eddington half range E_f/4 + y/2, the rest by
//! the rays' (A = (1 - th) h + th/4).  Eliminating E_f,
//!   y = (A_L e_L - A_R e_R)/(1 - (th_L + th_R)/2 + A_L r_L + A_R r_R):
//! the rays' kinetic flux of thin cells, the diffusion flux (E_L - E_R)/(3 chi d) of
//! thick ones, the Milne/Marshak flux of a thick cell against a thin one or a boundary
//! (whose incoming flux is the moment scheme's; with clean_bc = marshak A = q and th = 0
//! in the denominator).  The stress is the kinetic one of the same half ranges,
//! K+-_dd = E_f/6 +- 3y/8 (Eddington) and the rays' k+- E_f.  The E and stress rows take
//! these with weight mu.  The E bracket of clean_consistency stays (only its divergence
//! is defined).  With the conservative FV rays (clean_rays = fv1 | fv) the stress rows
//! also replace the rays' kinetic flux of y and its bracket, which carry the rays' too
//! bright sub-cell photosphere into y; the short-characteristics face fluxes are defined
//! only up to a divergence-free part, so with them the y bracket stays.  On the cubed
//! sphere the normal is the face's own and an interior face takes the rays' half ranges
//! about it (CLF_HK, CLF_KK; a seam ghost's family slots are in the neighbour's axes).

KOKKOS_INLINE_FUNCTION
void CleanCellMilne(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i, Real jl[4][4], Real jr[4][4], Real *s) {
  const Real mu = CleanCellMilneWeight(x, dir, m, k, j, i);
  if (mu == 0.0) return;
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const bool bdry = CleanBoundaryFace(x, dir, m, k, j, i);
  const bool blo = bdry && (dir == 0 ? i == x.is : (dir == 1 ? j == x.js : k == x.ks));
  const Real d = 0.5*CleanCellFaceDistance(x, dir, m, k, j, i);
  // sides n = 0 (L), 1 (R): A, r, c, e - c E_c, th and the stress half ranges; a
  // boundary's ghost side has none
  Real a[2] = {0.0, 0.0}, r[2] = {0.0, 0.0}, cc[2] = {0.0, 0.0}, e0[2] = {0.0, 0.0};
  Real th[2] = {0.0, 0.0}, ks[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, den = 1.0;
  Real hk[2] = {0.0, 0.0}, nv[3];
  if (!bdry) CleanCellKinetic(x, dir, m, k, j, i, hk[0], hk[1]);
  for (int b = 0; b < 3; ++b) nv[b] = CleanNf(x, b, dir, m, k, j, i);
  for (int n = 0; n < 2; ++n) {
    if (bdry && (n == 0) == blo) continue;
    const int kc = k - dk*(1 - n), jc = j - dj*(1 - n), ic = i - di*(1 - n);
    const int sg = 2*n - 1;   // from the face into the cell
    const Real chi = x.op(m,M1_OP_T,kc,jc,ic);
    const Real eps = (chi > 0.0) ? fmin(x.op(m,M1_OP_E,kc,jc,ic)/chi, 1.0) : 0.0;
    const Real tau = chi*d, z = std::sqrt(3.0*eps)*tau;
    const Real t = (z > 1.0e-4) ? tanh(z)/z : 1.0 - z*z/3.0;
    th[n] = -expm1(-2.0*tau);
    cc[n] = 1.0/cosh(z);
    r[n] = 3.0*tau*t;
    const Real bc = x.u(m,CLC_MB,kc,jc,ic);
    const Real bl = x.u(m,CLC_MB,kc-sg*dk,jc-sg*dj,ic-sg*di) - bc;
    const Real br = bc - x.u(m,CLC_MB,kc+sg*dk,jc+sg*dj,ic+sg*di);
    const Real bf = bc + ((bl*br > 0.0) ? copysign(0.5*fmin(fabs(bl), fabs(br)), bl)
                                        : 0.0);
    e0[n] = bf - cc[n]*bc - (bf - bc)*t;
    Real h = n ? -hk[1] : hk[0];
    if (bdry) {   // the kinetic outgoing flux (CleanFaceWeights), not the full column
      Real fl, fr;
      CleanFaceWeights(x, dir, m, k, j, i, fl, fr);
      h = fabs(blo ? fr : fl)/x.c;
    }
    const Real tm = (bdry && x.marshak) ? 0.0 : th[n];
    a[n] = (1.0 - tm)*h + 0.25*tm;
    den += a[n]*r[n] - 0.5*tm;
    for (int b = 0; b < x.nf; ++b) {
      const Real kr = (x.cub && !bdry) ? x.fc[dir](m,CLF_KK+3*n+b,k,j,i) :
                      CleanCellOutStress(x, dir, b, m, k, j, i, kc, jc, ic, n == 1);
      ks[n][b] = (1.0 - th[n])*kr + th[n]/6.0*nv[b];
    }
  }
  // y = yl E_L + yr E_R + y0, E_fL = c_L E_L + e0_L - r_L y, E_fR = c_R E_R + e0_R
  // + r_R y
  const Real yl = a[0]*cc[0]/den, yr = -a[1]*cc[1]/den, q = 0.375*(th[0] - th[1]);
  if (s != nullptr) {
    // the incoming flux of a boundary face (inside, s[0] is the bracket)
    const Real cb = x.consist ? x.fc[dir](m,CLF_CB,k,j,i) : 0.0;
    const Real y0 = (a[0]*e0[0] - a[1]*e0[1] + (s[0] - cb)/x.c)/den;
    s[0] = (1.0 - mu)*s[0] + mu*(cb + x.c*y0);
    for (int b = 0; b < x.nf; ++b) {
      const Real cy = (x.consist && x.fvrays) ? x.fc[dir](m,CLF_CB+1+b,k,j,i) : 0.0;
      s[1+b] += mu*(x.c*(ks[0][b]*(e0[0] - r[0]*y0) + ks[1][b]*(e0[1] + r[1]*y0) +
                         q*nv[b]*y0) - cy);
    }
    return;
  }
  // the E row and the stress columns E (the dissipation of y stays)
  for (int c = 0; c < 4; ++c) {
    jl[0][c] *= 1.0 - mu;
    jr[0][c] *= 1.0 - mu;
  }
  jl[0][0] += mu*x.c*yl;
  jr[0][0] += mu*x.c*yr;
  for (int b = 0; b < x.nf; ++b) {
    const Real g = q*nv[b];
    jl[1+b][0] = (1.0 - mu)*jl[1+b][0] +
                 mu*x.c*(ks[0][b]*(cc[0] - r[0]*yl) + (ks[1][b]*r[1] + g)*yl);
    jr[1+b][0] = (1.0 - mu)*jr[1+b][0] +
                 mu*x.c*(ks[1][b]*(cc[1] + r[1]*yr) + (g - ks[0][b]*r[0])*yr);
  }
  if (!bdry || !x.consist) return;
  // the boundary's deviation flux keeps its E part (CleanCellDeviationFlux)
  const int kc = blo ? k : k - dk, jc = blo ? j : j - dj, ic = blo ? i : i - di;
  const Real wo = mu*x.c*x.fc[dir](m,CLF_WK,k,j,i)*
                  x.cs(m,(blo ? CLN_HM : CLN_HP)+dir,kc,jc,ic);
  Real sh[3], hr[3];
  for (int b = 0; b < 3; ++b) {
    sh[b] = x.cs(m,CLN_HP+b,kc,jc,ic) + x.cs(m,CLN_HM+b,kc,jc,ic);
  }
  CleanRayFlux(x, m, kc, jc, ic, sh, hr);
  for (int b = 0; b < x.nf; ++b) (blo ? jr : jl)[1+b][0] -= wo*hr[b];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceSource
//! \brief the E-independent part of the face fluxes (CleanCellM1Source, CleanCellMilne)

KOKKOS_INLINE_FUNCTION
void CleanCellFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i, Real s[4]) {
  CleanCellM1Source(x, dir, m, k, j, i, s);
  if (CleanCellMilneWeight(x, dir, m, k, j, i) == 0.0) return;
  Real jl[4][4], jr[4][4];
  CleanCellM1Jacobian(x, dir, m, k, j, i, jl, jr);
  CleanCellMilne(x, dir, m, k, j, i, jl, jr, s);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceJacobian
//! \brief derivatives jl, jr of the face fluxes: CleanCellM1Jacobian and the blend of
//! CleanCellMilne, whose weight and terms CleanCellMilneSetup stores per face (CLF_ML)

KOKKOS_INLINE_FUNCTION
void CleanCellFaceJacobian(const CleanArgs &x, const int dir, const int m, const int k,
                           const int j, const int i, Real jl[4][4], Real jr[4][4]) {
  CleanCellM1Jacobian(x, dir, m, k, j, i, jl, jr);
  const Real mu = x.milne ? x.fc[dir](m,CLF_ML,k,j,i) : 0.0;
  if (mu == 0.0) return;
  for (int c = 0; c < 4; ++c) {
    jl[0][c] *= 1.0 - mu;
    jr[0][c] *= 1.0 - mu;
  }
  // with clean_consistency and FV rays the stress rows drop the rays' kinetic y flux too
  const bool ky = x.consist && x.fvrays;
  for (int r = 1; r <= x.nf; ++r) {
    for (int c = 0; c < (ky ? 4 : 1); ++c) {
      jl[r][c] *= 1.0 - mu;
      jr[r][c] *= 1.0 - mu;
    }
  }
  for (int r = 0; r <= x.nf; ++r) {
    jl[r][0] += x.fc[dir](m,CLF_ML+1+r,k,j,i);
    jr[r][0] += x.fc[dir](m,CLF_ML+5+r,k,j,i);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceFlux
//! \brief face fluxes f = (fE, fY_1..3) at face (k,j,i) normal to dir of the cell vector
//! at slot v of x.u; with_source adds the incoming boundary flux and, with
//! clean_consistency, the bracket

KOKKOS_INLINE_FUNCTION
void CleanCellFaceFlux(const CleanArgs &x, const int v, const int dir, const int m,
                       const int k, const int j, const int i, const bool with_source,
                       Real f[4]) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const int kl = k - dk, jl_ = j - dj, il = i - di;
  Real jl[4][4], jr[4][4], ul[4], ur[4];
  CleanCellFaceJacobian(x, dir, m, k, j, i, jl, jr);
  CleanCellState(x, v, m, kl, jl_, il, ul);
  CleanCellState(x, v, m, k, j, i, ur);
  for (int r = 0; r < 4; ++r) {
    f[r] = 0.0;
    for (int c = 0; c < 4; ++c) f[r] += jl[r][c]*ul[c] + jr[r][c]*ur[c];
  }
  const bool bdry = CleanBoundaryFace(x, dir, m, k, j, i);
  if (with_source && (bdry || x.consist || x.milne)) {
    Real s[4];
    CleanCellFaceSource(x, dir, m, k, j, i, s);
    for (int r = 0; r < 4; ++r) f[r] += s[r];
  }
  if (bdry) return;
  // RC terms in E_{L-1} and E_{R+1}
  const int kp = k + dk, jp = j + dj, ip = i + di;
  const int km = kl - dk, jm = jl_ - dj, im = il - di;
  Real al, hs;
  const bool lo = CleanFaceDissipation(x, dir, m, k, j, i, al, hs);
  const Real s = x.consist ? 1.0 - x.fc[dir](m,CLF_WK,k,j,i) : 1.0;
  const Real mu = x.milne ? x.fc[dir](m,CLF_ML,k,j,i) : 0.0;
  Real q[4];
  CleanCellRhieChow(x, dir, lo, m, k, j, i, q);
  hs *= s;
  f[0] += (1.0 - mu)*s*(q[3]*CleanDnn(x, dir, m, k, j, i, kp, jp, ip)*x.u(m,v,kp,jp,ip) +
                        q[0]*CleanDnn(x, dir, m, k, j, i, km, jm, im)*x.u(m,v,km,jm,im));
  // reconstruction terms of the dissipation in the states of L-1 and R+1
  if (!x.recon || lo) return;
  CleanCellState(x, v, m, km, jm, im, ul);
  CleanCellState(x, v, m, kp, jp, ip, ur);
  for (int c = 0; c <= x.nf; ++c) {
    Real cj[4];
    CleanJumpCoefficients(x, dir, c, m, k, j, i, cj);
    f[c] -= (c == 0 ? 1.0 - mu : 1.0)*hs*(cj[0]*ul[c] + cj[3]*ur[c]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellDivergence
//! \brief h div of the face fluxes of cell (k,j,i) for the cell vector at slot v (adiv
//! adds the sum of their magnitudes)

KOKKOS_INLINE_FUNCTION
void CleanCellDivergence(const CleanArgs &x, const int v, const int m, const int k,
                         const int j, const int i, const bool with_source, Real div[4],
                         Real *adiv = nullptr) {
  for (int r = 0; r < 4; ++r) div[r] = 0.0;
  for (int dir = 0; dir < 3; ++dir) {
    if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
    Real fp[4], fm[4];
    CleanCellFaceFlux(x, v, dir, m, k+dk, j+dj, i+di, with_source, fp);
    CleanCellFaceFlux(x, v, dir, m, k, j, i, with_source, fm);
    if (x.sph || x.cub) {
      Real gp[4], gm[4];
      CleanCellRows(x, dir, m, k+dk, j+dj, i+di, k, j, i, fp, gp);
      CleanCellRows(x, dir, m, k, j, i, k, j, i, fm, gm);
      for (int r = 0; r < 4; ++r) div[r] += gp[r] - gm[r];
      for (int r = 0; adiv && r < 4; ++r) adiv[r] += fabs(gp[r]) + fabs(gm[r]);
    } else {
      for (int r = 0; r < 4; ++r) div[r] += x.h*x.idx[dir]*(fp[r] - fm[r]);
      for (int r = 0; adiv && r < 4; ++r) {
        adiv[r] += x.h*x.idx[dir]*(fabs(fp[r]) + fabs(fm[r]));
      }
    }
  }
  // the curvature terms of the stress divergence (linear in E of the cell)
  if (x.sph) {
    Real g[3];
    CleanCellSphSource(x, m, k, j, i, g);
    for (int b = 0; b < x.nf; ++b) div[1+b] += g[b]*x.u(m,v,k,j,i);
    if (with_source && x.consist) {
      CleanCellSphSource(x, m, k, j, i, g, true);
      for (int b = 0; b < x.nf; ++b) div[1+b] += g[b];
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellSourceDivergence
//! \brief h div of the E-independent face fluxes of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
void CleanCellSourceDivergence(const CleanArgs &x, const int m, const int k, const int j,
                               const int i, Real div[4]) {
  for (int r = 0; r < 4; ++r) div[r] = 0.0;
  for (int dir = 0; dir < 3; ++dir) {
    if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
    Real s[4], g[4];
    if (x.consist || x.milne || CleanBoundaryFace(x, dir, m, k+dk, j+dj, i+di)) {
      CleanCellFaceSource(x, dir, m, k+dk, j+dj, i+di, s);
      CleanCellRows(x, dir, m, k+dk, j+dj, i+di, k, j, i, s, g);
      for (int r = 0; r < 4; ++r) div[r] += g[r];
    }
    if (x.consist || x.milne || CleanBoundaryFace(x, dir, m, k, j, i)) {
      CleanCellFaceSource(x, dir, m, k, j, i, s);
      CleanCellRows(x, dir, m, k, j, i, k, j, i, s, g);
      for (int r = 0; r < 4; ++r) div[r] -= g[r];
    }
  }
  if (x.sph && x.consist) {
    Real g[3];
    CleanCellSphSource(x, m, k, j, i, g, true);
    for (int b = 0; b < x.nf; ++b) div[1+b] += g[b];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanBlockSolve
//! \brief solves a b' = b in place for nb columns of b (4x4, partial pivoting)

KOKKOS_INLINE_FUNCTION
void CleanBlockSolve(Real a[4][4], Real b[4][5], const int nb) {
  for (int p = 0; p < 4; ++p) {
    int piv = p;
    for (int r = p + 1; r < 4; ++r) {
      if (fabs(a[r][p]) > fabs(a[piv][p])) piv = r;
    }
    if (piv != p) {
      for (int c = 0; c < 4; ++c) {
        const Real t = a[p][c];
        a[p][c] = a[piv][c];
        a[piv][c] = t;
      }
      for (int n = 0; n < nb; ++n) {
        const Real t = b[p][n];
        b[p][n] = b[piv][n];
        b[piv][n] = t;
      }
    }
    const Real inv = 1.0/a[p][p];
    for (int r = p + 1; r < 4; ++r) {
      const Real f = a[r][p]*inv;
      for (int c = p; c < 4; ++c) a[r][c] -= f*a[p][c];
      for (int n = 0; n < nb; ++n) b[r][n] -= f*b[p][n];
    }
  }
  for (int p = 3; p >= 0; --p) {
    for (int n = 0; n < nb; ++n) {
      Real s = b[p][n];
      for (int c = p + 1; c < 4; ++c) s -= a[p][c]*b[c][n];
      b[p][n] = s/a[p][p];
    }
  }
}

}  // namespace radm1
#endif  // RAD_M1_RAD_M1_CLEAN_CELL_HPP_
