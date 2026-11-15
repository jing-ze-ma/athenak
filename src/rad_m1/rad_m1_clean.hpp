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
constexpr int CLN_N = 27;

// face variables (arrays cln_fc1, cln_fc2, cln_fc3)
constexpr int CLF_W = 0;     // blend weight w
constexpr int CLF_TH = 1;    // theta = 1/(1 + h c chi_f)
constexpr int CLF_G = 2;     // flux of the old vector G_f
constexpr int CLF_HP = 3;    // h+_f
constexpr int CLF_HM = 4;    // h-_f
constexpr int CLF_N = 5;

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
constexpr int CLW_N = 25;

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
  Real h, c, idx[3];           // stage step, light speed, inverse cell widths
  int bc[6];                   // boundary types x1min, x1max, ..., x3max
  Real ebath[6];               // boundary bath energy densities
  Real marshak_q;
  bool marshak;                // Marshak boundary form (else the rays' half range)
  bool use_beam;
  bool feedback;               // gas coupling (advective flux, exchange linearisation)
  bool recon;                  // cell layout: reconstructed face states in the jumps
  bool all_low;                // cell layout: every interior face low-order
  bool consist;                // cell layout: face fluxes carry the consistency bracket
  int nf;                      // cell layout: number of flux components (dimensions)
  int is, ie, js, je, ks, ke;
  bool multi_d, three_d;
};

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
  int cp;                      // block-tridiagonal solve scratch
};

// slots of a coarse level (arrays cln_mg, no ghost cells)
constexpr int CLM_A0 = 0;    // diagonal block (row major)
constexpr int CLM_M = 16;    // block of the lower neighbour along dir: CLM_M + 32 dir
constexpr int CLM_P = 32;    // block of the upper neighbour along dir: CLM_P + 32 dir
constexpr int CLM_CP = 112;  // block-tridiagonal solve scratch
constexpr int CLM_B = 128;   // right-hand side
constexpr int CLM_X = 132;   // solution
constexpr int CLM_T1 = 136;  // residual
constexpr int CLM_T2 = 140;  // line-solve correction
constexpr int CLM_N = 144;

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
  if (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) {
    if (x.bc[2*dir] == M1_IBC_REFLECT) return;
    const bool beam = x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0;
    if (beam && x.marshak) return;
    ar = x.marshak ? -x.c*x.marshak_q : x.c*x.cs(m,CLN_HM+dir,k,j,i);
    return;
  }
  if (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC) {
    if (x.bc[2*dir+1] == M1_IBC_REFLECT) return;
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
bool CleanBoundaryFace(const CleanArgs &x, const int dir, const int k, const int j,
                       const int i) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  return (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) ||
         (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC);
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
  if (!x.multi_d || CleanBoundaryFace(x, dir, k, j, i)) return div;
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
  if (!x.multi_d || CleanBoundaryFace(x, dir, k, j, i)) return 0.0;
  const Real w = x.fc[dir](m,CLF_W,k,j,i);
  const Real th = x.fc[dir](m,CLF_TH,k,j,i);
  return -(1.0 - w)*th*x.h*x.c*x.c*CleanFaceStress(x, e, dir, m, k, j, i, false);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFaceSource
//! \brief the part of the face flux that does not depend on E: (1-w) theta G inside, the
//! incoming flux (beam, bath) on physical boundaries

KOKKOS_INLINE_FUNCTION
Real CleanFaceSource(const CleanArgs &x, const int dir, const int m, const int k,
                     const int j, const int i) {
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie + 1 : ((dir == 1) ? x.je + 1 : x.ke + 1);
  if (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) {
    if (x.bc[2*dir] == M1_IBC_REFLECT) return 0.0;
    if (x.use_beam && dir == 0 && x.beam(m,k-x.ks,j-x.js) > 0.0) {
      return x.beam(m,k-x.ks,j-x.js);
    }
    return x.marshak ? x.c*x.marshak_q*x.ebath[2*dir] : 0.25*x.c*x.ebath[2*dir];
  }
  if (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC) {
    if (x.bc[2*dir+1] == M1_IBC_REFLECT) return 0.0;
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
