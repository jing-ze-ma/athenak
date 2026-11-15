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
//!   fE   = c (y_d,L + y_d,R)/2 - al (s/2) (E_R - E_L) + (1 - al) RC
//!   fY_b = c (D_bd,L E_L + D_bd,R E_R)/2 - al (s/2) (y_b,R - y_b,L)
//!   RC   = th h c^2 [(g_L + g_R)/2 - ((D_dd E)_R - (D_dd E)_L)/dx_d]
//! with g_i = ((D_dd E)_{i+1} - (D_dd E)_{i-1})/(2 dx_d), th = 1/(1 + h c chi_f),
//! s = c sqrt(max(D_dd,L, D_dd,R)) and al the dissipation weight.  RC is the momentum
//! interpolation of the flux to the face: where h c chi >> 1 it turns fE into the compact
//! diffusion flux.  The jumps (E_R - E_L, y_R - y_L) are those of the face states: first
//! order, or (clean_recon = plm) the van Leer reconstruction of src/reconstruct/plm.hpp,
//! u_L^+ = u_L + dl dr/(dl + dr), whose weights dr/(dl + dr), dl/(dl + dr) are taken
//! from the Newton iterate (lagged) so that every pass stays a linear solve.

#include <cmath>

#include "athena.hpp"
#include "rad_m1/rad_m1_clean.hpp"

namespace radm1 {

// face variables of the cell layout (arrays cln_fc1..3; CLF_TH as in the staggered one)
constexpr int CLF_AL = 0;    // dissipation weight al
constexpr int CLF_AF = 2;    // dissipation speed s = c sqrt(max(D_dd,L, D_dd,R))
constexpr int CLF_LO = 3;    // 1: low-order face (first-order jumps, al = 1, s = c,
                             //    normal stress factor 1)

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
constexpr int CLC_N = 211;

// dissipation weight candidates (<rad_m1>/clean_alpha)
constexpr int CLA_TAU = 0;          // A(tau) = (1 - exp(-tau))/tau
constexpr int CLA_TAU_SENSOR = 1;   // max(A(tau) psi, s_ray)
constexpr int CLA_TAU_RAY = 2;      // max(A(tau), s_ray)
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
//! depth tau = chi_f dx_d, the smoothness sensor psi and the rays' streaming weight sray.
//! tau enters only here (numerical dissipation), never the closure.

KOKKOS_INLINE_FUNCTION
Real CleanDissipationWeight(const int mode, const Real tau, const Real psi,
                            const Real sray) {
  const Real a = (tau > 1.0e-8) ? -expm1(-tau)/tau : 1.0 - 0.5*tau;
  switch (mode) {
    case CLA_TAU: return a;
    case CLA_TAU_SENSOR: return fmax(a*psi, sray);
    case CLA_TAU_RAY: return fmax(a, sray);
    case CLA_ONE: return 1.0;
    default: return 0.0;
  }
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
//! a low-order face takes al = 1, s = c.  Returns true for a low-order face.

KOKKOS_INLINE_FUNCTION
bool CleanFaceDissipation(const CleanArgs &x, const int dir, const int m, const int k,
                          const int j, const int i, Real &al, Real &hs) {
  const bool lo = x.fc[dir](m,CLF_LO,k,j,i) > 0.0;
  al = lo ? 1.0 : x.fc[dir](m,CLF_AL,k,j,i);
  hs = 0.5*al*(lo ? x.c : x.fc[dir](m,CLF_AF,k,j,i));
  return lo;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellGhost
//! \brief incoming state outside the open physical boundary face (k,j,i) normal to dir:
//! the beam at x1min (E_g = F_b/(c sqrt(D_11)), y_g,b = D_1b F_b/(c D_11) with D of the
//! inner cell, exact for a collimated beam), else the isotropic bath (D = 1/3, y = 0).
//! Returns E_g, y_g and the stress column D_g,bd E_g.

KOKKOS_INLINE_FUNCTION
void CleanCellGhost(const CleanArgs &x, const int dir, const int m, const int k,
                    const int j, const int i, const bool low, Real &eg, Real yg[3],
                    Real pg[3]) {
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
  pg[dir] = eg/3.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceSource
//! \brief the E-independent part of the face fluxes on a physical boundary: the
//! incoming E flux of CleanFaceSource and the ghost part of the upwind stress flux

KOKKOS_INLINE_FUNCTION
void CleanCellFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                         const int j, const int i, Real s[4]) {
  for (int r = 0; r < 4; ++r) s[r] = 0.0;
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const bool low = (pos == lo);
  if (x.bc[2*dir + (low ? 0 : 1)] == M1_IBC_REFLECT) return;
  s[0] = CleanFaceSource(x, dir, m, k, j, i);
  // the ghost sits outside the face; its closure is read from the inner cell
  const int kc = low ? k : k - (dir == 2), jc = low ? j : j - (dir == 1);
  const int ic = low ? i : i - (dir == 0);
  Real eg, yg[3], pg[3];
  CleanCellGhost(x, dir, m, kc, jc, ic, low, eg, yg, pg);
  for (int b = 0; b < x.nf; ++b) {
    s[1+b] = 0.5*x.c*pg[b] + (low ? 0.5 : -0.5)*x.c*yg[b];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceJacobian
//! \brief derivatives jl, jr of the face fluxes (fE, fY_1..3) at face (k,j,i) normal to
//! dir with respect to the states (E, y_1..3) of the cells L and R.  The terms in the
//! states of L-1 and R+1 (RC, reconstruction) are left to CleanCellFaceFlux.  On
//! physical boundaries fE takes the boundary form of CleanFaceWeights and fY the
//! upwind (a = c) flux between the inner cell and the incoming ghost state of
//! CleanCellGhost; a reflecting wall carries only the normal stress of the inner cell.

KOKKOS_INLINE_FUNCTION
void CleanCellFaceJacobian(const CleanArgs &x, const int dir, const int m, const int k,
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
  const bool blo = (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC);
  const bool bhi = (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC);
  if (blo || bhi) {
    const int bc = blo ? x.bc[2*dir] : x.bc[2*dir+1];
    Real (*jj)[4] = blo ? jr : jl;
    const int kc = blo ? k : kl, jc = blo ? j : jl_, ic = blo ? i : il;
    if (bc == M1_IBC_REFLECT) {
      jj[1+dir][0] = x.c*CleanDiag(x, dir, dir, m, kc, jc, ic);
      return;
    }
    Real al, ar;
    CleanFaceWeights(x, dir, m, k, j, i, al, ar);
    jj[0][0] = blo ? ar : al;
    for (int b = 0; b < x.nf; ++b) {
      jj[1+b][0] = 0.5*x.c*CleanDiag(x, b, dir, m, kc, jc, ic);
      jj[1+b][1+b] = blo ? -0.5*x.c : 0.5*x.c;
    }
    return;
  }
  Real al, hs;
  const bool low = CleanFaceDissipation(x, dir, m, k, j, i, al, hs);
  const Real th = x.fc[dir](m,CLF_TH,k,j,i);
  const Real q = 0.75*(1.0 - al)*th*x.h*x.c*x.c*x.idx[dir];
  Real cj[4] = {0.0, -1.0, 1.0, 0.0};
  if (!low) CleanJumpCoefficients(x, dir, 0, m, k, j, i, cj);
  jl[0][0] = -hs*cj[1] + q*CleanDiag(x, dir, dir, m, kl, jl_, il);
  jr[0][0] = -hs*cj[2] - q*CleanDiag(x, dir, dir, m, k, j, i);
  jl[0][1+dir] = 0.5*x.c;
  jr[0][1+dir] = 0.5*x.c;
  for (int b = 0; b < x.nf; ++b) {
    if (!low) CleanJumpCoefficients(x, dir, 1+b, m, k, j, i, cj);
    // a low-order face carries the normal stress with pressure factor 1 (upwind
    // characteristics E +- F/c with speeds +-c)
    const bool chr = low && b == dir;
    jl[1+b][0] = 0.5*x.c*(chr ? 1.0 : CleanDiag(x, b, dir, m, kl, jl_, il));
    jr[1+b][0] = 0.5*x.c*(chr ? 1.0 : CleanDiag(x, b, dir, m, k, j, i));
    jl[1+b][1+b] = -hs*cj[1];
    jr[1+b][1+b] = -hs*cj[2];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellFaceFlux
//! \brief face fluxes f = (fE, fY_1..3) at face (k,j,i) normal to dir of the cell vector
//! at slot v of x.u; with_source adds the incoming boundary flux

KOKKOS_INLINE_FUNCTION
void CleanCellFaceFlux(const CleanArgs &x, const int v, const int dir, const int m,
                       const int k, const int j, const int i, const bool with_source,
                       Real f[4]) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const int kl = k - dk, jl_ = j - dj, il = i - di;
  Real jl[4][4], jr[4][4];
  CleanCellFaceJacobian(x, dir, m, k, j, i, jl, jr);
  for (int r = 0; r < 4; ++r) {
    f[r] = 0.0;
    for (int c = 0; c < 4; ++c) {
      f[r] += jl[r][c]*x.u(m,v+c,kl,jl_,il) + jr[r][c]*x.u(m,v+c,k,j,i);
    }
  }
  if (CleanBoundaryFace(x, dir, k, j, i)) {
    if (with_source) {
      Real s[4];
      CleanCellFaceSource(x, dir, m, k, j, i, s);
      for (int r = 0; r < 4; ++r) f[r] += s[r];
    }
    return;
  }
  // RC terms in E_{L-1} and E_{R+1}
  const int kp = k + dk, jp = j + dj, ip = i + di;
  const int km = kl - dk, jm = jl_ - dj, im = il - di;
  Real al, hs;
  const bool lo = CleanFaceDissipation(x, dir, m, k, j, i, al, hs);
  const Real th = x.fc[dir](m,CLF_TH,k,j,i);
  const Real q = 0.25*(1.0 - al)*th*x.h*x.c*x.c*x.idx[dir];
  f[0] += q*(CleanDiag(x, dir, dir, m, kp, jp, ip)*x.u(m,v,kp,jp,ip) -
             CleanDiag(x, dir, dir, m, km, jm, im)*x.u(m,v,km,jm,im));
  // reconstruction terms of the dissipation in the states of L-1 and R+1
  if (!x.recon || lo) return;
  for (int c = 0; c <= x.nf; ++c) {
    Real cj[4];
    CleanJumpCoefficients(x, dir, c, m, k, j, i, cj);
    f[c] -= hs*(cj[0]*x.u(m,v+c,km,jm,im) + cj[3]*x.u(m,v+c,kp,jp,ip));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellDivergence
//! \brief h div of the face fluxes of cell (k,j,i) for the cell vector at slot v

KOKKOS_INLINE_FUNCTION
void CleanCellDivergence(const CleanArgs &x, const int v, const int m, const int k,
                         const int j, const int i, const bool with_source, Real div[4]) {
  for (int r = 0; r < 4; ++r) div[r] = 0.0;
  for (int dir = 0; dir < 3; ++dir) {
    if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
    Real fp[4], fm[4];
    CleanCellFaceFlux(x, v, dir, m, k+dk, j+dj, i+di, with_source, fp);
    CleanCellFaceFlux(x, v, dir, m, k, j, i, with_source, fm);
    for (int r = 0; r < 4; ++r) div[r] += x.h*x.idx[dir]*(fp[r] - fm[r]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCellSourceDivergence
//! \brief h div of the E-independent boundary fluxes of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
void CleanCellSourceDivergence(const CleanArgs &x, const int m, const int k, const int j,
                               const int i, Real div[4]) {
  for (int r = 0; r < 4; ++r) div[r] = 0.0;
  for (int dir = 0; dir < 3; ++dir) {
    if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
    Real s[4];
    if (CleanBoundaryFace(x, dir, k+dk, j+dj, i+di)) {
      CleanCellFaceSource(x, dir, m, k+dk, j+dj, i+di, s);
      for (int r = 0; r < 4; ++r) div[r] += x.h*x.idx[dir]*s[r];
    }
    if (CleanBoundaryFace(x, dir, k, j, i)) {
      CleanCellFaceSource(x, dir, m, k, j, i, s);
      for (int r = 0; r < 4; ++r) div[r] -= x.h*x.idx[dir]*s[r];
    }
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
