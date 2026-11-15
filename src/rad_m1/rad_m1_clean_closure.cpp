//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_closure.cpp
//! \brief Closure of the clean implicit radiation solve.  Once per step the formal
//! solution (short characteristics) is run on the state predicted at the end of the step,
//! and its moments give the Eddington tensor D, the half-range factors h+- and the
//! non-locality C1 of every cell.

#include <algorithm>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"
#include "utils/imex_sdirk2.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn void CleanStoreMoments
//! \brief the closure of a cell from its ray moments q (CLN_* layout, not normalised):
//! D = K/J, h+- = H+-/J, k+ = K+/J, C1 = <|I - S|>/J; without rays (J = 0) isotropic
//! and marked for the fill.  The face fluxes are copied.

KOKKOS_INLINE_FUNCTION
void CleanStoreMoments(const DvceArray5D<Real> &cs, const DvceArray5D<Real> &valid,
                       const Real *q, const int m, const int k, const int j,
                       const int i) {
  if (q[CLN_J] > 0.0) {
    const Real ij = 1.0/q[CLN_J];
    for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = q[CLN_D11+a]*ij;
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_HP+a,k,j,i) = fmax(q[CLN_HP+a]*ij, 0.0);
      cs(m,CLN_HM+a,k,j,i) = fmin(q[CLN_HM+a]*ij, 0.0);
    }
    for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = q[CLN_KP+a]*ij;
    cs(m,CLN_C1,k,j,i) = fmax(q[CLN_C1]*ij, 0.0);
    cs(m,CLN_J,k,j,i) = q[CLN_J];
    valid(m,0,k,j,i) = 1.0;
  } else {
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
      cs(m,CLN_D11+3+a,k,j,i) = 0.0;
      cs(m,CLN_HP+a,k,j,i) = 0.25;
      cs(m,CLN_HM+a,k,j,i) = -0.25;
    }
    for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    cs(m,CLN_C1,k,j,i) = 1.0;
    cs(m,CLN_J,k,j,i) = 0.0;
    valid(m,0,k,j,i) = 0.0;
  }
  for (int n = CLN_F2; n < CLN_N; ++n) cs(m,n,k,j,i) = q[n];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceVector
//! \brief spherical polar: the exact vector area a of the low face of cell (k,j,i) normal
//! to dir (r faces along r_hat, theta faces along theta_hat, phi faces along phi_hat)

KOKKOS_INLINE_FUNCTION
void CleanFaceVector(const CleanArgs &x, const DvceArray2D<Real> &x3f, const int dir,
                     const int m, const int k, const int j, const int i, Real a[3]) {
  const Real r1 = x.x1f(m,i), r2 = x.x1f(m,i+1), t1 = x.x2f(m,j), t2 = x.x2f(m,j+1);
  const Real p1 = x3f(m,k), p2 = x3f(m,k+1);
  const Real hr = 0.5*(r2*r2 - r1*r1);
  if (dir == 0) {
    const Real s2 = 0.5*(t2 - t1) - 0.25*(sin(2.0*t2) - sin(2.0*t1));
    a[0] = r1*r1*s2*(sin(p2) - sin(p1));
    a[1] = r1*r1*s2*(cos(p1) - cos(p2));
    a[2] = r1*r1*0.5*(SQR(sin(t2)) - SQR(sin(t1)))*(p2 - p1);
  } else if (dir == 1) {
    a[0] = sin(t1)*hr*cos(t1)*(sin(p2) - sin(p1));
    a[1] = sin(t1)*hr*cos(t1)*(cos(p1) - cos(p2));
    a[2] = -sin(t1)*hr*sin(t1)*(p2 - p1);
  } else {
    a[0] = -hr*(t2 - t1)*sin(p1);
    a[1] = hr*(t2 - t1)*cos(p1);
    a[2] = 0.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFvFace
//! \brief clean_rays = fv: the intensity of direction d on the low face of cell (k,j,i)
//! normal to dir for a beam along +e_dir (pos) or -e_dir: that of the upwind cell, with
//! o2 its minmod slope towards the face (second order in index space, positive)

KOKKOS_INLINE_FUNCTION
Real CleanFvFace(const DvceArray5D<Real> &fi, const bool o2, const int m, const int d,
                 const int dir, const int k, const int j, const int i, const bool pos) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const int su = pos ? -1 : 0, sd = pos ? 0 : -1, sf = pos ? -2 : 1;
  const Real u = fi(m,d,k+su*dk,j+su*dj,i+su*di);
  if (!o2) return u;
  const Real a = u - fi(m,d,k+sf*dk,j+sf*dj,i+sf*di);
  const Real b = fi(m,d,k+sd*dk,j+sd*dj,i+sd*di) - u;
  return u + ((a*b > 0.0) ? 0.5*((fabs(a) < fabs(b)) ? a : b) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFvFaceAxi
//! \brief CleanFvFace on a phi face (index k + o) of a 2-D axisymmetric mesh: the cell s
//! columns away is the own cell with the directions rotated by s phi extents (mp: the
//! seam maps cln_fvm, side 2 below, 3 above)

KOKKOS_INLINE_FUNCTION
Real CleanFvFaceAxi(const DvceArray5D<Real> &fi, const DvceArray4D<int> &mp,
                    const bool o2, const int m, const int d, const int k, const int j,
                    const int i, const int o, const bool pos) {
  auto at = [&](const int s) {
    int e = d;
    for (int q = 0; q < ((s < 0) ? -s : s); ++q) e = mp(m, (s < 0) ? 2 : 3, k, e);
    return fi(m,e,k,j,i);
  };
  const Real u = at(o - (pos ? 1 : 0));
  if (!o2) return u;
  const Real a = u - at(o + (pos ? -2 : 1));
  const Real b = at(o - (pos ? 0 : 1)) - u;
  return u + ((a*b > 0.0) ? 0.5*((fabs(a) < fabs(b)) ? a : b) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanLocalBasis
//! \brief the local orthonormal basis (r_hat, theta_hat, phi_hat) at (th, ph)

KOKKOS_INLINE_FUNCTION
void CleanLocalBasis(const Real th, const Real ph, Real e[3][3]) {
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  e[0][0] = st*cp; e[0][1] = st*sp; e[0][2] = ct;
  e[1][0] = ct*cp; e[1][1] = ct*sp; e[1][2] = -st;
  e[2][0] = -sp;   e[2][1] = cp;    e[2][2] = 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosure
//! \brief builds the closure for t^{n+1}.  Both stage solves and a backward-Euler redo end
//! at t^{n+1}, so one closure serves the whole step; it is rebuilt for a new step or a
//! new dt.  clean_closure = eddington takes the isotropic closure instead of the rays.

void RadiationM1::CleanClosure() {
  if (cln_sphere) {
    if (!cln_havec) CleanClosureSphere();
    cln_havec = true;
    return;
  }
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_havec && time == cln_tc && dt == cln_dtb) return;

  if (cln_edd) {
    CleanClosureEddington();
  } else {
    CleanPredictState(true);
    CleanRays();
    if (cln_opsaved) Kokkos::deep_copy(opac, cln_opb);
    CleanClosureFill();
  }
  CleanGhostFill(cln_s, 0, CLN_N);

  cln_tc = time;
  cln_dtb = dt;
  cln_havec = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureEddington
//! \brief clean_closure = eddington: the isotropic closure in every cell (D = I/3,
//! h+- = +-1/4, k+ = I/6, C1 = 0), with no formal solution

void RadiationM1::CleanClosureEddington() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto cs = cln_s;
  par_for("cln_edd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
      cs(m,CLN_D11+3+a,k,j,i) = 0.0;
      cs(m,CLN_HP+a,k,j,i) = 0.25;
      cs(m,CLN_HM+a,k,j,i) = -0.25;
    }
    for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    cs(m,CLN_C1,k,j,i) = 0.0;
    cs(m,CLN_J,k,j,i) = 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureIterate
//! \brief closure of the current Newton/Picard iterate (clean_closure = iterate): the
//! formal solution of the iterate's E, gas temperature and F, with the stage
//! opacities, so that the Eddington tensor belongs to the state being solved for

void RadiationM1::CleanClosureIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool cell = cln_cell;
  const Real cl = c_light;
  auto iw_ = iw;
  auto u0_ = u0;
  auto w = cln_w;
  auto u = cln_x;
  par_for("cln_clositer", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_EN,k,j,i) = fmax(cell ? u(m,CLC_K,k,j,i) : w(m,CLW_EK,k,j,i), 0.0);
    iw_(m,M1_IW_TP,k,j,i) = w(m,CLW_T,k,j,i);
    iw_(m,M1_IW_F1,k,j,i) = cell ? cl*u(m,CLC_K+1,k,j,i) : u0_(m,M1_F1,k,j,i);
    iw_(m,M1_IW_F2,k,j,i) = cell ? cl*u(m,CLC_K+2,k,j,i) : u0_(m,M1_F1+1,k,j,i);
    iw_(m,M1_IW_F3,k,j,i) = cell ? cl*u(m,CLC_K+3,k,j,i) : u0_(m,M1_F1+2,k,j,i);
  });
  CleanRays();
  CleanClosureFill();
  CleanGhostFill(cln_s, 0, CLN_N);
  cln_havec = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRayMoments
//! \brief formal solution of the current state (u0 and the gas, no prediction) and its
//! moments J, K = D J, H+, H-, <|I - S|> and the ray face fluxes per cell into dst (CLN_*
//! layout, not normalised; J = 0 marks a cell whose rays do not carry its radiation).
//! With solve = false the moments of the closure in place (the rays of the last formal
//! solution) are taken.

void RadiationM1::CleanRayMoments(DvceArray5D<Real> dst, const bool solve) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  if (solve) {
    CleanPredictState(false);
    CleanRays();
    if (cln_opsaved) Kokkos::deep_copy(opac, cln_opb);
  }

  auto cs = cln_s;
  par_for("cln_raymom", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = cs(m,CLN_J,k,j,i);
    for (int n = 0; n < CLN_J; ++n) dst(m,n,k,j,i) = jr*cs(m,n,k,j,i);
    for (int n = CLN_J; n < CLN_N; ++n) dst(m,n,k,j,i) = cs(m,n,k,j,i);
  });
  cln_havec = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureAfterSolve
//! \brief clean_closure = extrap(2): the ray moments of a converged solution.  The last
//! solve of a step gives the moments at t^{n+1} (used from the next step on); with
//! extrap2 the stage-1 solve gives those at the stage-1 time.  With extrap1 only a
//! backward-Euler step runs the formal solution here; after a stage step the closure of
//! the stage-1 solve (rays of its first iterate, a state at t^{n+1}) gives the moments.

void RadiationM1::CleanClosureAfterSolve() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_extrap == 3 && t2_solve == M1_T2S_STAGE2) return;
  if (cln_extrap == 3 && t2_solve == M1_T2S_STAGE1) {
    CleanRayMoments(cln_mc, false);
    cln_tm = time + dt;
    cln_have_c = true;
    cln_cand_be = false;
    return;
  }
  if (t2_solve == M1_T2S_STAGE1) {
    if (cln_extrap != 2) return;
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    CleanRayMoments(cln_m1);
    cln_t1 = time + (tb.a21 + tb.g)*dt;
    cln_have_1 = true;
    return;
  }
  CleanRayMoments(cln_mc);
  cln_tm = time + dt;
  cln_have_c = true;
  cln_cand_be = (t2_solve != M1_T2S_STAGE2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureExtrap
//! \brief clean_closure = extrap(2): the closure at the stage time from ray moments of
//! converged solutions, extrapolated linearly in time.
//!
//! The formal solution is linear in the source and the boundary intensity, so the
//! extrapolated moments are those of the extrapolated intensity, exact when S and the
//! inflow vary linearly in time; the closure D = K/J, h+- = H+-/J, C1 is then formed from
//! them.  The line passes through the moments at t^{n-1} and t^n (extrap), or, for stage
//! 2 with extrap2, through those at t^n and of the stage-1 solution.  A cell keeps the
//! moments at t^n where the extrapolated ones are not those of a non-negative intensity
//! (J <= 0, |H| > J or K - H H^T/J not positive semidefinite; counted).  The first step,
//! a restart and backward-Euler solves use the moments at t^n alone; without any (first
//! step, restart) the closure of the state at t^n is built (CleanClosure).  With extrap1
//! the moments are those of the stage-1 closures and stage 2 keeps the stage-1 closure.

void RadiationM1::CleanClosureExtrap() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  cln_nguard = 0;
  if (cln_extrap == 3 && t2_solve == M1_T2S_STAGE2) return;

  // a new step: the end of the last step becomes t^n
  if (cln_have_c && cln_tm == time) {
    std::swap(cln_mp, cln_mn);
    std::swap(cln_mn, cln_mc);
    cln_tp = cln_tn;
    cln_tn = time;
    cln_have_p = cln_have_n && !cln_cand_be;
    cln_have_n = true;
    cln_have_c = false;
    cln_have_1 = false;
  }
  if (!cln_have_n || cln_tn != time) {
    CleanClosure();
    return;
  }

  // the line M(t) = Mb + (t - tb)/(tb - ta) (Mb - Ma) at the stage time ts
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  Real ts = time + dt;
  if (t2_solve == M1_T2S_STAGE1) {
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    ts = time + (tb.a21 + tb.g)*dt;
  }
  const bool use1 = (cln_extrap == 2) && (t2_solve == M1_T2S_STAGE2) && cln_have_1 &&
                    (cln_t1 > cln_tn);
  const bool line = stage && (use1 || cln_have_p);
  auto ma = use1 ? cln_mn : cln_mp;
  auto mb = use1 ? cln_m1 : cln_mn;
  const Real tb = use1 ? cln_t1 : cln_tn;
  const Real ta = use1 ? cln_tn : cln_tp;
  const Real wx = line ? (ts - tb)/(tb - ta) : 0.0;
  auto mn = cln_mn;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  auto cs = cln_s;
  auto valid = cln_v;
  int nguard = 0;
  Kokkos::parallel_reduce("cln_extrap",
                          Kokkos::RangePolicy<DevExeSpace>(0, (nmb1 + 1)*nk*nj*ni),
  KOKKOS_LAMBDA(const int t, int &ng) {
    const int i = is + t%ni;
    const int j = js + (t/ni)%nj;
    const int k = ks + (t/(ni*nj))%nk;
    const int m = t/(ni*nj*nk);
    Real q[CLN_N];
    for (int n = 0; n < CLN_N; ++n) {
      q[n] = mb(m,n,k,j,i) + wx*(mb(m,n,k,j,i) - ma(m,n,k,j,i));
    }
    // the extrapolated moments must be those of a non-negative intensity
    bool ok = line && ma(m,CLN_J,k,j,i) > 0.0 && q[CLN_J] > 0.0;
    if (ok) {
      const Real ij = 1.0/q[CLN_J];
      Real h[3], a[3][3];
      for (int b = 0; b < 3; ++b) h[b] = (q[CLN_HP+b] + q[CLN_HM+b])*ij;
      a[0][0] = q[CLN_D11]*ij;
      a[1][1] = q[CLN_D11+1]*ij;
      a[2][2] = q[CLN_D11+2]*ij;
      a[0][1] = a[1][0] = q[CLN_D11+3]*ij;
      a[0][2] = a[2][0] = q[CLN_D11+4]*ij;
      a[1][2] = a[2][1] = q[CLN_D11+5]*ij;
      for (int b = 0; b < 3; ++b) {
        for (int c = 0; c < 3; ++c) a[b][c] -= h[b]*h[c];
      }
      // all principal minors of D - h h^T non-negative (to round-off)
      const Real eps = 1.0e-12;
      const Real m01 = a[0][0]*a[1][1] - a[0][1]*a[0][1];
      const Real m02 = a[0][0]*a[2][2] - a[0][2]*a[0][2];
      const Real m12 = a[1][1]*a[2][2] - a[1][2]*a[1][2];
      const Real det = a[0][0]*m12 - a[0][1]*(a[0][1]*a[2][2] - a[1][2]*a[0][2])
                       + a[0][2]*(a[0][1]*a[1][2] - a[1][1]*a[0][2]);
      ok = (h[0]*h[0] + h[1]*h[1] + h[2]*h[2] <= 1.0 + eps) &&
           (a[0][0] >= -eps) && (a[1][1] >= -eps) && (a[2][2] >= -eps) &&
           (m01 >= -eps) && (m02 >= -eps) && (m12 >= -eps) && (det >= -eps);
    }
    if (!ok) {
      if (line && mn(m,CLN_J,k,j,i) > 0.0) ++ng;
      for (int n = 0; n < CLN_N; ++n) q[n] = mn(m,n,k,j,i);
    }
    CleanStoreMoments(cs, valid, q, m, k, j, i);
  }, Kokkos::Sum<int>(nguard));
  cln_nguard = static_cast<int>(CleanReduce(nguard, false));
  cln_nguard_sum += cln_nguard;
  if (line) cln_next_sum += 1.0;

  CleanClosureFill();
  CleanGhostFill(cln_s, 0, CLN_N);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPredictState
//! \brief the state at t^{n+1} predicted from the physics of this step, as input of the
//! formal solution (iw: E, T, F) and, for temperature-dependent opacities, opac.
//!
//! The hydro's explicit stage already holds the gas at t^{n+1} without the radiation
//! coupling, and K1 = f(U^n) (the FSAL slope of the implicit part) gives the coupling's
//! rates at U^n:
//!  1. momentum: hydro + dt K1_M.  Where radiation supports the gas this cancels the
//!     explicit gravity kick; the kinetic energy is taken from this momentum.
//!  2. radiation transport: dt (K1_E + K1_EN) = -dt div F, less the work of the radiation
//!     force on the gas.  A loss larger than the content is applied at constant relative
//!     rate.
//!  3. exchange: the implicit local relaxation of (eint, E) over dt, from the hydro eint
//!     and the transported E (stiff and local, so solved as such; no overshoot).
//! Without gas feedback the gas is held fixed and E advances at its own rate.
//! Density is not changed by radiation within a step.  The prediction is within O(dt^2)
//! of U^{n+1}.  Without a valid K1 (first step) or with use_k1 = false the coupling
//! rates are zero and the formal solution sees the current state itself.

void RadiationM1::CleanPredictState(const bool use_k1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real dt = pmy_pack->pmesh->dt;

  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fluid = fl.on;
  const bool rates = use_k1 && cln_pred && t2_ok && (t2k1.extent_int(0) > nmb1);

  // opacities at the predicted state; the solve itself uses those of the stage
  cln_opsaved = fluid && !opac_zero && (opacity_type != M1_OPAC_CONST);
  if (cln_opsaved) {
    if (cln_opb.extent_int(1) != opac.extent_int(1)) {
      Kokkos::realloc(cln_opb, opac.extent_int(0), opac.extent_int(1),
                      opac.extent_int(2), opac.extent_int(3), opac.extent_int(4));
    }
    Kokkos::deep_copy(cln_opb, opac);
  }
  const bool new_opac = cln_opsaved;

  auto iw_ = iw;
  auto u0_ = u0;
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto k1 = t2k1;
  auto opac_ = opac;
  int otype = opacity_type;
  Real kp = kappa_p, ke_ = kappa_e, kf = kappa_f, ksc = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  M1OpacTab otab_ = otab;
  bool scat = vscat;
  Real arad_ = arad, cl = c_light;
  const bool feedback = cln_feedback;

  par_for("cln_predict", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // K1 is stored on active cells only; ghost cells take the explicit state
    const bool active = (k >= ks && k <= ke && j >= js && j <= je && i >= is && i <= ie);
    const bool use_rates = rates && active;
    const Real e0 = u0_(m,M1_E,k,j,i);
    Real e = e0, d = 1.0, t = 0.0;
    if (fluid) {
      d = uh(m,IDN,k,j,i);
      Real mom[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
      const Real ekin = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
      const Real eint_h = fmax(uh(m,IEN,k,j,i) - ekin, 0.0);
      t = eos.Temperature(d, eint_h);
      if (use_rates && !feedback) {
        // gas held fixed: E advances at its own rate, which includes the exchange
        const Real de = dt*k1(m,M1_T2_E,k,j,i);
        e = (e0 + de > 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
      } else if (use_rates) {
        // radiation transport over dt, less the work of the radiation force on the gas,
        // then the implicit local exchange of (eint, E) over dt
        for (int n = 0; n < 3; ++n) mom[n] += dt*k1(m,M1_T2_M1+n,k,j,i);
        const Real ekin_new = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
        const Real de = dt*(k1(m,M1_T2_E,k,j,i) + k1(m,M1_T2_EN,k,j,i))
                        - (ekin_new - ekin);
        e = (e0 + de > 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
        // E (1 + h c kE) = e + h c kP a T^4,  eint(T) + E = eint_h + e
        const Real hce = dt*cl*opac_(m,M1_OP_E,k,j,i);
        const Real hcp = dt*cl*opac_(m,M1_OP_P,k,j,i)*arad_;
        t = CleanExchangeTemperature(eos, d, hcp/(1.0 + hce), eint_h + e*hce/(1.0 + hce),
                                     t);
        e = (e + hcp*SQR(SQR(t)))/(1.0 + hce);
      }
    } else if (use_rates) {
      const Real de = dt*k1(m,M1_T2_E,k,j,i);
      e = (e0 + de > 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
    }

    if (active) {
      iw_(m,M1_IW_EN,k,j,i) = e;
      iw_(m,M1_IW_TP,k,j,i) = t;
      iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
      iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F1+1,k,j,i);
      iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F1+2,k,j,i);
    }
    if (new_opac) {
      Real op, oe, of, os;
      if (otype == M1_OPAC_TABLE) {
        M1TableOpacities(otab_, d, t, op, oe, of, os);
        M1ScatterEos(otab_, eos, d, t, of, os);
      } else {
        M1Opacities(otype, d, t, kp, ke_, kf, ksc, rref, tref, aa, bb, op, oe, of, os);
      }
      opac_(m,M1_OP_P,k,j,i) = d*op;
      opac_(m,M1_OP_E,k,j,i) = d*oe;
      opac_(m,M1_OP_T,k,j,i) = d*(of + os);
      if (scat) opac_(m,M1_OP_S,k,j,i) = d*os;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureSphere
//! \brief spherical polar: the exact ray moments of a uniformly bright sphere of radius
//! R = clean_sphere_radius about the origin seen through vacuum.  Outside it the rays
//! fill the cone mu > mu_c = sqrt(1 - R^2/r^2) about r_hat: D_rr = (1 + mu_c + mu_c^2)/3,
//! D_tt = D_pp = (1 - D_rr)/2, h+_r = (1 + mu_c)/2, h-_r = 0,
//! h+-_t = h+-_p = +-S/(pi (1 - mu_c)) with S = int_{mu_c}^1 sqrt(1 - mu^2) dmu, C1 = 1,
//! and the half-range stress columns k+^r = (D_rr, 0, 0), k+^t = (k, D_tt/2, 0),
//! k+^p = (k, 0, D_pp/2) with k = (R/r) (1 + mu_c)/(3 pi); R -> 0 is the point source.
//! Inside (r <= R) the radiation is isotropic: D = I/3, h+- = +-1/4, k+ = I/6, C1 = 0.

void RadiationM1::CleanClosureSphere() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real rs = cln_sph_r;
  auto x1v = pmy_pack->pcoord->x1v;
  auto cs = cln_s;

  par_for("cln_sphere", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    Real drr = 1.0/3.0, hpr = 0.25, hmr = -0.25, ht = 0.25, c1 = 0.0, kr = 0.0,
         krr = 1.0/6.0;
    if (r > rs) {
      const Real q = rs/r;
      const Real mc = std::sqrt(fmax(1.0 - q*q, 0.0));
      drr = (1.0 + mc + mc*mc)/3.0;
      hpr = 0.5*(1.0 + mc);
      hmr = 0.0;
      // S/(pi (1 - mu_c)) with 1 - mu_c = q^2/(1 + mu_c); S ~ 2 q^3/3 for small q
      const Real s = 0.5*(std::asin(q) - mc*q);
      ht = (q > 1.0e-4) ? s*(1.0 + mc)/(M_PI*q*q) : 2.0*q/(3.0*M_PI);
      c1 = 1.0;
      kr = q*(1.0 + mc)/(3.0*M_PI);
      krr = drr;
    }
    cs(m,CLN_D11,k,j,i) = drr;
    cs(m,CLN_D11+1,k,j,i) = 0.5*(1.0 - drr);
    cs(m,CLN_D11+2,k,j,i) = 0.5*(1.0 - drr);
    for (int a = 3; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = 0.0;
    cs(m,CLN_HP,k,j,i) = hpr;
    cs(m,CLN_HM,k,j,i) = hmr;
    for (int a = 1; a < 3; ++a) {
      cs(m,CLN_HP+a,k,j,i) = ht;
      cs(m,CLN_HM+a,k,j,i) = -ht;
    }
    for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = 0.0;
    cs(m,CLN_KP,k,j,i) = krr;
    cs(m,CLN_KP+3,k,j,i) = kr;
    cs(m,CLN_KP+4,k,j,i) = 0.25*(1.0 - drr);
    cs(m,CLN_KP+6,k,j,i) = kr;
    cs(m,CLN_KP+8,k,j,i) = 0.25*(1.0 - drr);
    cs(m,CLN_C1,k,j,i) = c1;
    cs(m,CLN_J,k,j,i) = 1.0;
  });
  CleanGhostFill(cln_s, 0, CLN_N);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRays
//! \brief the formal solution of the state in iw (E, T, F; chi from opac) and its
//! closure: vet_sc on a Cartesian mesh, vet_gd from the first shell on spherical polar.
//! The source function is the state's own (no relaxed source or relaxed scattering J,
//! which would make the closure depend on dt).

void RadiationM1::CleanRays() {
  vcol_srelax = false;
  vcol_sbl = false;
  vcol_sesc = false;
  vscat_jrel = false;
  if (!vgd_on) {
    VetShortChar();
    CleanClosureFromRays();
    return;
  }
  if (!tau_ready) VetColInit();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto iw_ = iw;
  auto opac_ = opac;
  par_for("cln_gd_chi", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
  });
  if (cln_rfv > 0) {
    VetGdSource();
    CleanFvSweep();
  } else {
    VetGdBuild();
  }
  CleanClosureFromGd();
  if (cln_consist) CleanGdFaceSums();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvSeams
//! \brief the ghosts of cln_fvi beyond a periodic theta or phi seam of the mesh (not a
//! whole circle in phi) take their directions in the frame of their image cell
//! (cln_fvm); the ghost exchange has copied the image cell's intensities

void RadiationM1::CleanFvSeams() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int n = vgd_n;
  const Real pp = pmy_pack->pmesh->mesh_size.x3max - pmy_pack->pmesh->mesh_size.x3min;
  const bool circ = std::abs(pp - 2.0*M_PI) < 1.0e-10;
  CleanArgs x = CleanMakeArgs();
  auto fi = cln_fvi;
  auto mp = cln_fvm;
  const int n1 = fi.extent_int(4), n2 = fi.extent_int(3), n3 = fi.extent_int(2);
  if (cln_fvt.extent_int(1) != n) Kokkos::realloc(cln_fvt, n3, n2, n1, n);
  auto t = cln_fvt;
  for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
  par_for("cln_fv_seam", DevExeSpace(), 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int k, const int j, const int i) {
    int sd = -1;
    if (j < js && x.mbc(m,2) == BoundaryFlag::periodic) sd = 0;
    if (j > je && x.mbc(m,3) == BoundaryFlag::periodic) sd = 1;
    if (!circ && k < ks && x.mbc(m,4) == BoundaryFlag::periodic) sd = 2;
    if (!circ && k > ke && x.mbc(m,5) == BoundaryFlag::periodic) sd = 3;
    if (sd < 0 || (sd < 2 && (k < ks || k > ke))) return;
    for (int d = 0; d < n; ++d) t(k,j,i,d) = fi(m,mp(m,sd,k,d),k,j,i);
    for (int d = 0; d < n; ++d) fi(m,d,k,j,i) = t(k,j,i,d);
  });
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::CleanFvWallInit
//! \brief clean_rays = fv | fv1: whether the low (side 0) or high (side 1) theta face of
//! the mesh is a reflecting wall; at the first call the specular map of its walls
//! (cln_fvwm, per MeshBlock and column k): each incoming set direction takes the leaving
//! direction nearest to its mirror image n - 2 (n . t) t about the wall's normal t at
//! the face; cln_fvwu holds |n . t| (t the face's unit normal, the flux factor of every
//! beam through it).

bool RadiationM1::CleanFvWallInit(const int side) {
  Mesh *pm = pmy_pack->pmesh;
  if (pm->use_polar_boundary || pm->mesh_bcs[side ? BoundaryFace::outer_x2
                                                  : BoundaryFace::inner_x2] !=
      BoundaryFlag::reflect) {
    return false;
  }
  if (cln_fvwm.extent_int(0) > 0) return true;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack, n = vgd_n, c3 = cln_s.extent_int(2);
  Kokkos::realloc(cln_fvwm, nmb, 2, c3, n);
  Kokkos::realloc(cln_fvwu, nmb, 2, c3, n);
  Kokkos::realloc(cln_fvw, nmb, 2, c3, cln_s.extent_int(4), n);
  Kokkos::deep_copy(cln_fvw, 0.0);
  auto m_h = Kokkos::create_mirror_view(cln_fvwm);
  auto u_h = Kokkos::create_mirror_view(cln_fvwu);
  auto d_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_dir);
  auto t_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx2f);
  auto p_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx3f);
  for (int m = 0; m < nmb; ++m) {
    for (int sd = 0; sd < 2; ++sd) {
      const Real th = t_h(m, sd ? indcs.je + 1 : indcs.js), sgn = sd ? 1.0 : -1.0;
      for (int k = indcs.ks; k <= indcs.ke; ++k) {
        // the unit normal of the wall face (CleanFaceVector: the mean theta_hat over the
        // face's phi range)
        const Real p1 = p_h(m,k), p2 = p_h(m,k+1);
        Real t[3] = {cos(th)*(sin(p2) - sin(p1)), cos(th)*(cos(p1) - cos(p2)),
                     -sin(th)*(p2 - p1)};
        const Real tn = std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
        for (int c = 0; c < 3; ++c) t[c] /= tn;
        for (int d = 0; d < n; ++d) {
          const Real nt = d_h(d,0)*t[0] + d_h(d,1)*t[1] + d_h(d,2)*t[2];
          u_h(m,sd,k,d) = std::abs(nt);
          m_h(m,sd,k,d) = -1;
          if (!(sgn*nt < 0.0)) continue;
          const Real v[3] = {d_h(d,0) - 2.0*nt*t[0], d_h(d,1) - 2.0*nt*t[1],
                             d_h(d,2) - 2.0*nt*t[2]};
          Real bd = -2.0;
          for (int e = 0; e < n; ++e) {
            const Real et = d_h(e,0)*t[0] + d_h(e,1)*t[1] + d_h(e,2)*t[2];
            const Real dt = v[0]*d_h(e,0) + v[1]*d_h(e,1) + v[2]*d_h(e,2);
            if (sgn*et > 0.0 && dt > bd) {bd = dt; m_h(m,sd,k,d) = e;}
          }
        }
      }
    }
  }
  Kokkos::deep_copy(cln_fvwm, m_h);
  Kokkos::deep_copy(cln_fvwu, u_h);
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvWalls
//! \brief the incoming intensities of the reflecting theta walls (cln_fvw): each incoming
//! direction e takes the intensity of the wall cell along its mirror direction, scaled
//! on each wall face so that it returns exactly the energy it receives,
//! sum_e w_e |n_e . t| I(e) = sum_out w_d |n_d . t| I_d (isotropic radiation stays
//! isotropic: the set integrates n exactly)

void RadiationM1::CleanFvWalls() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, n = vgd_n;
  CleanArgs x = CleanMakeArgs();
  auto fi = cln_fvi;
  auto fw = cln_fvw;
  auto mp = cln_fvwm;
  auto u = cln_fvwu;
  auto dir_ = vgd_dir;
  par_for("cln_fv_wall", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, 1, indcs.ks,
          indcs.ke, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int i) {
    const int jc = sd ? je : js;
    if (!CleanBoundaryFace(x, 1, m, k, sd ? je + 1 : js, i)) return;
    Real so = 0.0, si = 0.0;
    for (int e = 0; e < n; ++e) {
      const int d = mp(m,sd,k,e);
      const Real q = dir_(e,3)*u(m,sd,k,e)*fi(m,(d >= 0) ? d : e,k,jc,i);
      fw(m,sd,k,i,e) = (d >= 0) ? fi(m,d,k,jc,i) : 0.0;
      if (d >= 0) {si += q;} else {so += q;}
    }
    for (int e = 0; e < n; ++e) fw(m,sd,k,i,e) *= (si > 0.0) ? so/si : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvSweep
//! \brief clean_rays = fv | fv1: the conservative formal solution on the vet_gd
//! direction set (cln_fvi).  Per direction every face carries n . A_f times its upwind
//! intensity (CleanFvFace: the upwind cell's, with fv its minmod slope), A_f the exact
//! vector area, so that a cell passes on what enters it plus its exchange:
//! sum_out n . A_f I_f + V chi I = sum_in |n . A_f| I_f + V chi S.  A beam keeps its
//! energy; the cell intensities are those of the transported flux.  The
//! cells are taken in the order of n . x of their centres (Gauss-Seidel, repeated with
//! the ghost exchange until the largest change is below 1e-13 of the largest intensity).
//! With fv at most 8 passes per build (the slopes converge over the builds of a run).
//! The inner face emits the sweep's emission (vet_bc_inner emit) or the S of the first
//! cell, the top is vacuum, a reflecting theta wall is specular (CleanFvWalls), another
//! physical theta or phi face is closed.

void RadiationM1::CleanFvSweep() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ncell = nx1*nx2*nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  const bool axi = !pmy_pack->pmesh->three_d;
  CleanArgs x = CleanMakeArgs();
  auto x3f = pmy_pack->pcoord->xx3f;
  if (cln_fvi.extent_int(1) != n) {
    // the intensities, the outward face vectors of every cell and the cell order of
    // every direction (geometry and directions are fixed)
    Kokkos::realloc(cln_fvi, nmb, n, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
    Kokkos::realloc(cln_fva, nmb, 18, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
    auto fa = cln_fva;
    par_for("cln_fv_geo", DevExeSpace(), 0, nmb - 1, ks, indcs.ke, js, indcs.je, is,
            indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int f = 0; f < 6; ++f) {
        const int dir = f/2, up = f%2;
        Real a[3];
        CleanFaceVector(x, x3f, dir, m, k + (up && dir == 2), j + (up && dir == 1),
                        i + (up && dir == 0), a);
        for (int c = 0; c < 3; ++c) fa(m,3*f+c,k,j,i) = up ? a[c] : -a[c];
      }
    });
    Kokkos::realloc(cln_fvo, nmb, n, ncell);
    auto o_h = Kokkos::create_mirror_view(cln_fvo);
    auto d_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_dir);
    auto r_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x1v);
    auto t_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x2v);
    auto p_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x3v);
    std::vector<Real> key(ncell);
    std::vector<int> idx(ncell);
    for (int m = 0; m < nmb; ++m) {
      for (int d = 0; d < n; ++d) {
        for (int c = 0; c < ncell; ++c) {
          const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
          const Real st = std::sin(t_h(m,j)), r = r_h(m,i);
          key[c] = r*(d_h(d,0)*st*std::cos(p_h(m,k)) + d_h(d,1)*st*std::sin(p_h(m,k)) +
                      d_h(d,2)*std::cos(t_h(m,j)));
          idx[c] = c;
        }
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {return key[a] < key[b];});
        for (int c = 0; c < ncell; ++c) o_h(m,d,c) = idx[c];
      }
    }
    Kokkos::deep_copy(cln_fvo, o_h);
    Kokkos::deep_copy(cln_fvi, 0.0);
    // a periodic theta or phi seam of a wedge: the ghosts' directions in the frame of
    // their image cell, n' = sum_b (n . e_b(ghost)) e_b(image), the nearest set direction
    // (an exact permutation for a phi period of the set)
    const Real pt = pmy_pack->pmesh->mesh_size.x2max - pmy_pack->pmesh->mesh_size.x2min;
    const Real pp = pmy_pack->pmesh->mesh_size.x3max - pmy_pack->pmesh->mesh_size.x3min;
    const int c3 = cln_s.extent_int(2);
    Kokkos::realloc(cln_fvm, nmb, 4, c3, n);
    auto m_h = Kokkos::create_mirror_view(cln_fvm);
    for (int m = 0; m < nmb; ++m) {
      for (int sd = 0; sd < 4; ++sd) {
        for (int k = 0; k < c3; ++k) {
          // sd: theta low, theta high, phi low, phi high; the ghost one cell outside
          const bool th = (sd < 2), lo = (sd%2 == 0);
          const Real tg = th ? (lo ? t_h(m,js) - 1.0e-3 : t_h(m,indcs.je) + 1.0e-3)
                             : 0.5*M_PI;
          const Real pg = th ? p_h(m,k) : 0.0;
          const Real ti = th ? tg + (lo ? pt : -pt) : tg;
          const Real pi_ = th ? pg : pg + (lo ? pp : -pp);
          Real eg[3][3], ei[3][3];
          CleanLocalBasis(tg, pg, eg);
          CleanLocalBasis(ti, pi_, ei);
          for (int d = 0; d < n; ++d) {
            Real v[3] = {0.0, 0.0, 0.0};
            for (int b = 0; b < 3; ++b) {
              const Real nb = d_h(d,0)*eg[b][0] + d_h(d,1)*eg[b][1] + d_h(d,2)*eg[b][2];
              for (int c = 0; c < 3; ++c) v[c] += nb*ei[b][c];
            }
            int best = 0;
            Real bd = -2.0;
            for (int e = 0; e < n; ++e) {
              const Real dt = v[0]*d_h(e,0) + v[1]*d_h(e,1) + v[2]*d_h(e,2);
              if (dt > bd) {bd = dt; best = e;}
            }
            m_h(m,sd,k,d) = best;
          }
        }
      }
    }
    Kokkos::deep_copy(cln_fvm, m_h);
  }
  auto fi = cln_fvi;
  auto fa = cln_fva;
  auto o_ = cln_fvo;
  auto fm = cln_fvm;
  const bool o2 = (cln_rfv == 2);
  const bool wall[2] = {CleanFvWallInit(0), CleanFvWallInit(1)};
  auto fw = cln_fvw;
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  const bool emi = VetBcEmit();
  const Real iemit = (ibc_x1min == M1_IBC_FLUX) ? 4.0*iflux_x1min/c_light : iebath_x1min;
  // (fv: the slopes make the sweep nonlinear; a few passes per build, each starting from
  // the intensities of the last build)
  for (int it = 0; it < (o2 ? 8 : 100); ++it) {
    Real dmax = 0.0, imax = 0.0;
    Kokkos::parallel_reduce("cln_fv", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*n),
    KOKKOS_LAMBDA(const int md, Real &dm, Real &im) {
      const int m = md/n, d = md%n;
      const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
      for (int p = 0; p < ncell; ++p) {
        const int c = o_(m,d,p);
        const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
        // the outflow faces' slope parts are taken at the current iterate
        Real in = 0.0, out = 0.0, cor = 0.0;
        const Real ic = fi(m,d,k,j,i);
        for (int f = 0; f < 6; ++f) {
          const int dir = f/2, up = f%2;
          const int kf = k + (up && dir == 2), jf = j + (up && dir == 1);
          const int i_f = i + (up && dir == 0);
          const Real a = nv[0]*fa(m,3*f,k,j,i) + nv[1]*fa(m,3*f+1,k,j,i) +
                         nv[2]*fa(m,3*f+2,k,j,i);
          // (2-D axisymmetric: the phi faces open, their neighbours the own cell)
          const bool ax = axi && dir == 2;
          const bool bdry = !ax && CleanBoundaryFace(x, dir, m, kf, jf, i_f);
          // a reflecting theta wall: specular (its incoming beams from cln_fvw)
          const int sw = (dir == 1 && bdry && wall[up]) ? up : -1;
          if (bdry && dir != 0 && sw < 0) continue;
          if (sw >= 0 && a < 0.0) {
            in -= a*fw(m,sw,k,i,d);
            continue;
          }
          if (a >= 0.0) {
            out += a;
            cor += a*((ax ? CleanFvFaceAxi(fi, fm, o2, m, d, k, j, i, up, up)
                          : CleanFvFace(fi, o2, m, d, dir, kf, jf, i_f, up)) - ic);
          } else if (!bdry) {
            in -= a*(ax ? CleanFvFaceAxi(fi, fm, o2, m, d, k, j, i, up, !up)
                        : CleanFvFace(fi, o2, m, d, dir, kf, jf, i_f, !up));
          } else if (!up) {
            in -= a*(emi ? iemit : exp(cw_(m,1,k+ok,j+og,i)));
          }
        }
        const Real vc = x.vol(m,k,j,i)*exp(cw_(m,0,k+ok,j+og,i));
        const Real sc = vc*exp(cw_(m,1,k+ok,j+og,i));
        const Real iv = (out + vc > 0.0) ? fmax((in - cor + sc)/(out + vc), 0.0) : 0.0;
        dm = fmax(dm, fabs(iv - fi(m,d,k,j,i)));
        im = fmax(im, iv);
        fi(m,d,k,j,i) = iv;
      }
    }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(imax));
    CleanGhostFill(cln_fvi, 0, n);
    CleanFvSeams();
    if (wall[0] || wall[1]) CleanFvWalls();
    if (CleanReduce(dmax, true) <= 1.0e-13*CleanReduce(imax, true)) break;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFromGd
//! \brief the closure from the vet_gd intensities in the local basis (r, theta, phi) of
//! each cell, mu_b = n . e_b: J = <I>, K_bc = <I mu_b mu_c>, H+-_b = <I mu_b> over mu_b
//! >< 0, K+^d_b = <I mu_b mu_d> over mu_d > 0 and C1 = <|I - S|>/J (CLN_* layout, as
//! CleanClosureFromRays; no ray face fluxes).

void RadiationM1::CleanClosureFromGd() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
  const bool fv = (cln_rfv > 0);
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;
  par_for("cln_gd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real st = sin(x2v(m,j)), ct = cos(x2v(m,j));
    const Real sp = sin(x3v(m,k)), cp = cos(x3v(m,k));
    const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
    Real q[CLN_N] = {0.0};
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      const Real mu[3] = {nx*st*cp + ny*st*sp + nz*ct, nx*ct*cp + ny*ct*sp - nz*st,
                          -nx*sp + ny*cp};
      const Real iv = fv ? fi_(m,d,k,j,i) : vi_(m,d,k+ok,j+og,i);
      const Real a = dir_(d,3)*iv;
      q[CLN_J] += a;
      for (int b = 0; b < 3; ++b) {
        q[CLN_D11+b] += a*mu[b]*mu[b];
        q[((mu[b] > 0.0) ? CLN_HP : CLN_HM)+b] += a*mu[b];
        if (!(mu[b] > 0.0)) continue;
        for (int c = 0; c < 3; ++c) q[CLN_KP+3*b+c] += a*mu[c]*mu[b];
      }
      q[CLN_D11+3] += a*mu[0]*mu[1];
      q[CLN_D11+4] += a*mu[0]*mu[2];
      q[CLN_D11+5] += a*mu[1]*mu[2];
      q[CLN_C1] += dir_(d,3)*fabs(iv - s0);
    }
    if (q[CLN_J] < ray_min*iw_(m,M1_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    CleanStoreMoments(cs, valid, q, m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGdFaceSums
//! \brief the vet_gd rays' face fluxes for clean_consistency (CLN_FX, CLN_F2, CLN_F3,
//! CLN_F1: sum_d w (n_d . A_f) I_f / A_f and the stress columns in the face's local
//! basis), A_f the exact vector area of the face (closed for every cell).  A theta or
//! phi face carries the intensity of its upwind cell (none on a physical boundary).  The
//! r faces of a column take the cell balances from the inner face: the outer r face of a
//! cell gets what its inner one carries plus the exchange V chi (S - J, -H) less the
//! lateral outflow (the momentum in the column's basis), so that the energy balance of
//! every cell is the vet_gd exchange to round-off.  CLN_G is what the momentum balance of
//! the rays (in the cell's basis) differs by from the cell layout's divergence of these
//! face components: the curvature of their stress.  Several MeshBlocks along x1: a level
//! starts from the top face of the level below (vgd_x1b -> vgd_x1r), one level per pass;
//! at a reflecting top the remainder is spread over the outermost level (zero wall flux).

void RadiationM1::CleanGdFaceSums() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
  auto fm_ = cln_fvm;
  const bool fv = (cln_rfv > 0), o2 = (cln_rfv == 2);
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  auto x3f = pmy_pack->pcoord->xx3f;
  auto x3v = pmy_pack->pcoord->x3v;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  const bool wint = vgd_wint;
  auto cs = cln_s;
  CleanArgs x = CleanMakeArgs();
  // nx3 = 1 (axisymmetric): no phi face is a boundary; the upwind cell of the low one is
  // the own cell with the directions rotated back by the phi extent (the wall map), the
  // high one carries the low one's sums in its own basis
  const bool axi = !pmy_pack->pmesh->three_d;
  auto bnd = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i) {
    return !(axi && dir == 2) && CleanBoundaryFace(x, dir, m, k, j, i);
  };
  // the inner face: no flux at a reflecting one, else the outward beams of the sweep's
  // emission (vet_bc_inner emit) or of the first cell, the inward ones of the first cell
  const bool open = (ibc_x1min != M1_IBC_REFLECT);
  const bool emi = VetBcEmit();
  const Real iemit = (ibc_x1min == M1_IBC_FLUX) ? 4.0*iflux_x1min/c_light : iebath_x1min;
  // the short characteristics' intensities in the mesh layout with the ghost exchange
  // (the image cell's beyond a periodic seam, re-indexed by the wall map as the sweep)
  if (!fv) {
    if (cln_fvi.extent_int(1) != n) {
      Kokkos::realloc(cln_fvi, nmb1 + 1, n, cs.extent_int(2), cs.extent_int(3),
                      cs.extent_int(4));
    }
    fi_ = cln_fvi;
    par_for("cln_gd_copy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < n; ++d) fi_(m,d,k,j,i) = vi_(m,d,k+ok,j+og,i);
    });
    CleanGhostFill(cln_fvi, 0, n);
  }
  auto inten = [=] KOKKOS_FUNCTION (const int m, const int d, const int k, const int j,
                                    const int i) {
    return fi_(m,d,k,j,i);
  };
  // vector area of the low face of cell (k,j,i) normal to dir and the local basis there
  auto face = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                   const int i, Real a[3], Real e[3][3]) {
    CleanFaceVector(x, x3f, dir, m, k, j, i, a);
    CleanLocalBasis((dir == 1) ? x.x2f(m,j) : x.x2v(m,j),
                    (dir == 2) ? x3f(m,k) : x3v(m,k), e);
  };
  // beam of direction d through the low theta or phi face of cell (k,j,i)
  auto lat = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d) {
    if (bnd(dir, m, k, j, i)) return 0.0;
    Real a[3], e[3][3];
    face(dir, m, k, j, i, a, e);
    const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
    const int kl = k - (dir == 2), jl = j - (dir == 1);
    if (fv && axi && dir == 2) {
      return f*CleanFvFaceAxi(fi_, fm_, o2, m, d, k, j, i, 0, f > 0.0);
    }
    if (fv) return f*CleanFvFace(fi_, o2, m, d, dir, k, j, i, f > 0.0);
    if (!(f > 0.0)) return f*inten(m, d, k, j, i);
    if (axi && dir == 2) return f*inten(m, mp_(m,kl+ok,j+og,d), k, j, i);
    const bool seam = (dir == 1) ? (j == js && x.mbc(m,2) == BoundaryFlag::periodic)
                                 : (k == ks && x.mbc(m,4) == BoundaryFlag::periodic);
    if (!seam || wl_(m,kl+ok,jl+og) == 0) return f*inten(m, d, kl, jl, i);
    if (!wint) return f*inten(m, mp_(m,kl+ok,jl+og,d), kl, jl, i);
    Real v = 0.0;
    for (int t = 0; t < 3; ++t) {
      v += w3_(m,kl+ok,jl+og,d,t)*inten(m, m3_(m,kl+ok,jl+og,d,t), kl, jl, i);
    }
    return f*v;
  };
  // sums of the beams of every direction through the face, per unit area
  auto add = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d, const Real f, Real q[4]) {
    Real a[3], e[3][3];
    face(dir, m, k, j, i, a, e);
    const Real w = dir_(d,3)*f/x.area[dir](m,k,j,i);
    q[0] += w;
    for (int b = 0; b < 3; ++b) {
      q[1+b] += w*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
    }
  };
  // a mirror at the r face i, per unit area: no energy, the normal momentum of the
  // beams of cell ic that hit it (f of sign sg) doubled, no tangential part
  auto mirror = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                     const int ic, const Real sg, Real q[4]) {
    Real a[3], e[3][3], w[4] = {0.0};
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      if (sg*f > 0.0) add(0, m, k, j, i, d, 2.0*f*inten(m, d, k, j, ic), w);
    }
    q[0] = 0.0;
    q[1] = w[1];
    q[2] = 0.0;
    q[3] = 0.0;
  };
  for (int dir = 1; dir < 3; ++dir) {
    const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
    par_for("cln_gd_lat", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4] = {0.0};
      for (int d = 0; d < n; ++d) add(dir, m, k, j, i, d, lat(dir, m, k, j, i, d), q);
      for (int c = 0; c < 4; ++c) cs(m,nf+c,k,j,i) = q[c];
    });
  }
  CleanGhostFill(cln_s, CLN_F2, 4);
  CleanGhostFill(cln_s, CLN_F3, 4);
  // the conservative rays (clean_rays = fv | fv1): the r face i carries their own upwind
  // fluxes (area times the per-area sums), not the balance march, which loses the
  // diffusion flux in the difference V chi (S - J) of thick cells
  auto rflux = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                    Real q[4]) {
    Real a[3], e[3][3];
    for (int c = 0; c < 4; ++c) q[c] = 0.0;
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      Real iv = CleanFvFace(fi_, o2, m, d, 0, k, j, i, f > 0.0);
      if (CleanBoundaryFace(x, 0, m, k, j, i)) {   // the sweep's x1 boundary
        iv = (i == is) ? ((f > 0.0) ? (emi ? iemit : exp(cw_(m,1,k+ok,j+og,is))) : iv)
                       : ((f > 0.0) ? iv : 0.0);
      }
      add(0, m, k, j, i, d, f*iv, q);
    }
    for (int c = 0; c < 4; ++c) q[c] *= x.area[0](m,k,j,i);
  };
  auto xf_ = vgd_x1f;
  auto qb_ = vgd_x1b;
  auto qr_ = vgd_x1r;
  const int nj = je - js + 1, nkj = (ke - ks + 1)*nj;
  for (int lev = 0; lev < vgd_nb1; ++lev) {
    if (lev > 0) VetGdX1Pass(true);
    par_for("cln_gd_col", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      Real a[3], e[3][3], ef[3][3], q[4] = {0.0};
      const bool in = (xf_(m) & 1);
      const int kj = (k - ks)*nj + j - js;
      face(0, m, k, j, is, a, e);
      for (int d = 0; open && in && d < n; ++d) {
        const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
        const Real iv = (emi && f > 0.0) ? iemit : inten(m, d, k, j, is);
        add(0, m, k, j, is, d, f*iv, q);
      }
      if (in && !open) mirror(m, k, j, is, is, -1.0, q);
      for (int c = 0; c < 4; ++c) {
        q[c] = in ? q[c]*x.area[0](m,k,j,is) : qr_(m,c*nkj + kj);
      }
      if (fv) rflux(m, k, j, is, q);
      const Real vb = in ? 0.0 : qr_(m,4*nkj + kj);   // the column's volume below
      for (int i = is; i <= ie; ++i) {
        for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
        Real g[3] = {0.0, 0.0, 0.0}, qi[3] = {q[1], q[2], q[3]};
        // the exchange V chi (S - J, -H_b) and the lateral outflow of the cell
        const Real chi = x.vol(m,k,j,i)*exp(cw_(m,0,k+ok,j+og,i));
        const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
        for (int d = 0; d < n; ++d) {
          const Real di = dir_(d,3)*inten(m, d, k, j, i);
          q[0] += chi*(dir_(d,3)*s0 - di);
          for (int b = 0; b < 3; ++b) {
            q[1+b] -= chi*di*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
          }
        }
        for (int dir = 1; dir < 3; ++dir) {
          const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
          for (int up = 0; up < 2; ++up) {
            const int kk = k + (up && dir == 2), jj = j + (up && dir == 1);
            if (bnd(dir, m, kk, jj, i)) continue;
            face(dir, m, kk, jj, i, a, ef);
            const Real sg = (up ? -1.0 : 1.0)*x.area[dir](m,kk,jj,i);
            const int kc = (axi && dir == 2) ? k : kk;
            q[0] += sg*cs(m,nf,kc,jj,i);
            for (int b = 0; b < 3; ++b) {
              const Real wt = (dir == 2 && axi) ? 0.0 : ((dir == 1 && b == 2) ?
                              sin(x.x2f(m,jj))/sin(x.x2v(m,j)) : 1.0);
              g[b] += sg*wt*cs(m,nf+1+b,kc,jj,i);
              for (int c = 0; c < 3; ++c) {
                const Real r = sg*cs(m,nf+1+c,kc,jj,i)*(ef[c][0]*e[b][0] +
                               ef[c][1]*e[b][1] + ef[c][2]*e[b][2]);
                q[1+b] += r;
                g[b] -= r;
              }
            }
          }
        }
        if (fv) rflux(m, k, j, i + 1, q);
        // the angular-momentum weights r_f/r of the y_theta, y_phi rows on the r faces
        for (int b = 1; b < 3; ++b) {
          g[b] += q[1+b]*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) -
                  qi[b]*(1.0 - x.x1f(m,i)/x.x1v(m,i));
        }
        for (int b = 0; b < 3; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
      }
      if (xf_(m) & 4) {
        // a reflecting top: the rays' remainder at the wall (their non-conservation, all
        // but the mirror's momentum) spread over the cells of the column by volume, so
        // that no energy leaves through it (the levels below: in the pass down)
        Real qw[4];
        mirror(m, k, j, ie + 1, ie, 1.0, qw);
        for (int c = 0; c < 4; ++c) q[c] -= qw[c]*x.area[0](m,k,j,ie+1);
        Real vt = vb, vc = vb;
        for (int i = is; i <= ie; ++i) vt += x.vol(m,k,j,i);
        for (int i = is + (in ? 1 : 0); i <= ie; ++i) {
          if (i > is) vc += x.vol(m,k,j,i-1);
          for (int c = 0; c < 4; ++c) {
            cs(m,CLN_FX+c,k,j,i) -= q[c]*vc/(vt*x.area[0](m,k,j,i));
          }
        }
        for (int c = 0; c < 4; ++c) qb_(m,c*nkj + kj) = q[c];
        qb_(m,4*nkj + kj) = vt;
        qb_(m,5*nkj + kj) = vb;
        for (int c = 0; c < 4; ++c) q[c] = qw[c]*x.area[0](m,k,j,ie+1);
      }
      for (int c = 0; c < 4; ++c) {
        cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
        cs(m,CLN_F1+c,k,j,ie) = q[c]/x.area[0](m,k,j,ie+1);
        if (!(xf_(m) & 4)) qb_(m,c*nkj + kj) = q[c];
      }
      if (!(xf_(m) & 4)) {
        Real v = vb;
        for (int i = is; i <= ie; ++i) v += x.vol(m,k,j,i);
        qb_(m,4*nkj + kj) = v;
      }
    });
  }
  // a reflecting top over several levels: the remainder, the column volume and the
  // volume below the level above handed down; each level takes its share
  for (int lev = vgd_nb1 - 2; lev >= 0 && ibc_x1max == M1_IBC_REFLECT; --lev) {
    VetGdX1Pass(false);
    par_for("cln_gd_top", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if ((xf_(m) >> 3) != lev) return;
      const int kj = (k - ks)*nj + j - js;
      const Real vt = qr_(m,4*nkj + kj);
      Real vc = qr_(m,5*nkj + kj);
      for (int i = is; i <= ie; ++i) vc -= x.vol(m,k,j,i);
      for (int c = 0; c < 6; ++c) qb_(m,c*nkj + kj) = qr_(m,c*nkj + kj);
      qb_(m,5*nkj + kj) = vc;
      for (int i = is; i <= ie + 1; ++i) {
        if (i > is) vc += x.vol(m,k,j,i-1);
        for (int c = 0; c < 4; ++c) {
          const Real dq = qr_(m,c*nkj + kj)*vc/vt;
          if (i <= ie) cs(m,CLN_FX+c,k,j,i) -= dq/x.area[0](m,k,j,i);
          if (i == is) cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
          if (i == ie + 1) cs(m,CLN_F1+c,k,j,ie) -= dq/x.area[0](m,k,j,ie+1);
        }
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFromRays
//! \brief closure of every cell whose rays carry its radiation (J >= clean_ray_min E):
//! D = K/J, h+- = H+-/J, k+ = K+/J, C1 = <|I - S|>/J.  Other cells are marked for the
//! fill.  The rays' face fluxes (CLN_F2 ... CLN_F3) are copied for every cell.

void RadiationM1::CleanClosureFromRays() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto vc = vet_cell;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;

  const bool f2 = cln_consist;
  par_for("cln_rays", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = vc(m,M1_VET_J,k,j,i);
    for (int n = 0; n < M1_VET_NFS; ++n) {
      cs(m,CLN_F2+n,k,j,i) = f2 ? vc(m,M1_VET_F2+n,k,j,i) : 0.0;
    }
    if (jr > 0.0 && jr >= ray_min*iw_(m,M1_IW_EN,k,j,i)) {
      for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = vc(m,M1_VET_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,M1_VET_HP1+a,k,j,i)/jr;
        const Real h = vc(m,M1_VET_H1+a,k,j,i)/jr;
        cs(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cs(m,CLN_HM+a,k,j,i) = fmin(h - hp, 0.0);
      }
      for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = vc(m,M1_VET_KP+a,k,j,i)/jr;
      cs(m,CLN_C1,k,j,i) = fmax(vc(m,M1_VET_NL,k,j,i)/jr, 0.0);
      cs(m,CLN_J,k,j,i) = jr;
      valid(m,0,k,j,i) = 1.0;
    } else {
      // isotropic until the fill gives it the closure of its neighbours
      for (int a = 0; a < 3; ++a) {
        cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
        cs(m,CLN_D11+3+a,k,j,i) = 0.0;
        cs(m,CLN_HP+a,k,j,i) = 0.25;
        cs(m,CLN_HM+a,k,j,i) = -0.25;
      }
      for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
      cs(m,CLN_C1,k,j,i) = 1.0;
      cs(m,CLN_J,k,j,i) = 0.0;
      valid(m,0,k,j,i) = 0.0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFill
//! \brief cells without ray radiation take the mean closure (D, h+-, k+) of their valid
//! neighbours, one layer per pass (clean_fill_passes passes)

void RadiationM1::CleanClosureFill() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int dj = pmy_pack->pmesh->multi_d ? 1 : 0;
  const int dk = pmy_pack->pmesh->three_d ? 1 : 0;
  auto cs = cln_s;
  auto valid = cln_v;

  for (int pass = 0; pass < cln_fillpass; ++pass) {
    CleanGhostFill(cln_s, 0, CLN_N);
    CleanGhostFill(cln_v, 0, 2);
    // valid(1) marks the cells filled in this pass; only cells valid before it are read
    par_for("cln_fill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,1,k,j,i) = valid(m,0,k,j,i);
      if (valid(m,0,k,j,i) > 0.0) return;
      Real sum[CLN_C1] = {0.0};
      int nv = 0;
      for (int kk = k-dk; kk <= k+dk; ++kk) {
        for (int jj = j-dj; jj <= j+dj; ++jj) {
          for (int ii = i-1; ii <= i+1; ++ii) {
            if (valid(m,0,kk,jj,ii) <= 0.0) continue;
            ++nv;
            for (int n = 0; n < CLN_C1; ++n) sum[n] += cs(m,n,kk,jj,ii);
          }
        }
      }
      if (nv == 0) return;
      valid(m,1,k,j,i) = 1.0;
      for (int n = 0; n < CLN_C1; ++n) cs(m,n,k,j,i) = sum[n]/nv;
    });
    par_for("cln_fill_mark", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,0,k,j,i) = valid(m,1,k,j,i);
    });
  }
}

}  // namespace radm1
