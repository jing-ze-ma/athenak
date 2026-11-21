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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "utils/imex_sdirk2.hpp"
#include "utils/scratch_view.hpp"
// M1_TU_PART: this file is compiled as 3 translation units (rad_m1_clean_closure.cpp =
// part 0, rad_m1_clean_closure_p<k>.cpp = part k, see split_tu.py); every member function
// below is in exactly one part.
#ifndef M1_TU_PART_N
#define M1_TU_PART_N 0
#endif
#define M1_TU_PART(k) (M1_TU_PART_N == (k))

namespace radm1 {
extern bool g_vra_adv;   // vet_rotate_advance (rad_m1_vetgd.cpp)
//----------------------------------------------------------------------------------------
//! \fn void CleanStoreMoments
//! \brief the closure of a cell from its ray moments q (CLN_* layout, not normalised):
//! D = K/J, h+- = H+-/J, k+ = K+/J, C1 = <|I - S|>/J; without rays (J = 0) isotropic
//! and marked for the fill.  The face fluxes are copied.

KOKKOS_INLINE_FUNCTION
void CleanStoreMoments(const DvceArray5D<Real> &cs, const DvceArray5D<Real> &valid,
                       const Real *q, const int m, const int k, const int j,
                       const int i) {
  const int ns = cs.extent_int(1);   // (the staggered layout: no k+, CleanInit)
  if (q[CLN_J] > 0.0) {
    const Real ij = 1.0/q[CLN_J];
    for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = q[CLN_D11+a]*ij;
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_HP+a,k,j,i) = fmax(q[CLN_HP+a]*ij, 0.0);
      cs(m,CLN_HM+a,k,j,i) = fmin(q[CLN_HM+a]*ij, 0.0);
    }
    for (int a = 0; a < 9 && ns > CLN_KP; ++a) cs(m,CLN_KP+a,k,j,i) = q[CLN_KP+a]*ij;
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
    for (int a = 0; a < 9 && ns > CLN_KP; ++a) {
      cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    }
    cs(m,CLN_C1,k,j,i) = 1.0;
    cs(m,CLN_J,k,j,i) = 0.0;
    valid(m,0,k,j,i) = 0.0;
  }
  for (int n = CLN_F2; n < CLN_KP && n < ns; ++n) cs(m,n,k,j,i) = q[n];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceVector
//! \brief spherical polar: the exact vector area a of the low face of cell (k,j,i) normal
//! to dir (r faces along r_hat, theta faces along theta_hat, phi faces along phi_hat);
//! on the cubed sphere that of CleanCsGeometry; Cartesian the area along e_dir

KOKKOS_INLINE_FUNCTION
void CleanFaceVector(const CleanArgs &x, const DvceArray2D<Real> &x3f, const int dir,
                     const int m, const int k, const int j, const int i, Real a[3]) {
  if (x.cub) {   // the exact closed face vectors of CleanCsGeometry
    for (int c = 0; c < 3; ++c) a[c] = x.area[dir](m,k,j,i)*x.g.At<CLG_N>(m,3*dir+c,k,j,i);
    return;
  }
  if (!x.sph) {
    for (int c = 0; c < 3; ++c) a[c] = (c == dir) ? x.area[dir](m,k,j,i) : 0.0;
    return;
  }
  // (a theta or phi face reads only its own angle: the high face of a 1-D (2-D) cell has
  // no j+1 (k+1) edge)
  const Real r1 = x.x1f(m,i), r2 = x.x1f(m,i+1), t1 = x.x2f(m,j);
  const Real t2 = (dir == 1) ? t1 : x.x2f(m,j+1);
  const Real p1 = x3f(m,k), p2 = (dir == 2) ? p1 : x3f(m,k+1);
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
//! \fn Real CleanThickI
//! \brief clean_ray_cut: the intensity of direction n in a cell below the cut, the
//! diffusion limit I = J + 3 n . H of the cell's E en and F f (iw, the state the rays are
//! built for) in its frame e (CleanThickFrame); its moments are D = I/3 and h+- = +-1/4
//! + h/2 of the set's quadrature

KOKKOS_INLINE_FUNCTION
void CleanThickFrame(const CleanArgs &x, const DvceArray2D<Real> &x3v, const int m,
                     const int k, const int j, const int i, Real e[3][3]) {
  if (x.cub) {
    CleanCsFrame(x, m, k, j, i, e);
  } else {
    CleanLocalBasis(x.x2v(m,j), x3v(m,k), e);
  }
}

KOKKOS_INLINE_FUNCTION
Real CleanThickI(const Real e[3][3], const Real f[3], const Real en, const Real ic,
                 const Real n[3]) {
  Real nf = 0.0;
  for (int b = 0; b < 3; ++b) nf += (n[0]*e[b][0] + n[1]*e[b][1] + n[2]*e[b][2])*f[b];
  return fmax(en + 3.0*ic*nf, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosure
//! \brief builds the closure for t^{n+1}.  Both stage solves and a backward-Euler redo
//! end at t^{n+1}, so one closure serves the whole step; it is rebuilt for a new step or
//! a new dt.  clean_closure = eddington takes the isotropic closure instead of the rays.

#if M1_TU_PART(0)
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
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));

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
  CleanFor("cln_edd", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
      cs(m,CLN_D11+3+a,k,j,i) = 0.0;
      cs(m,CLN_HP+a,k,j,i) = 0.25;
      cs(m,CLN_HM+a,k,j,i) = -0.25;
    }
    for (int a = 0; a < 9 && cs.extent_int(1) > CLN_KP; ++a) {
      cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    }
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
  auto iw_ = iw;
  auto u0_ = u0;
  auto w = cln_w;
  par_for("cln_clositer", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_EN,k,j,i) = fmax(w(m,CLW_EK,k,j,i), 0.0);
    iw_(m,M1_IW_TP,k,j,i) = w(m,CLW_T,k,j,i);
    iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
    iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F1+1,k,j,i);
    iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F1+2,k,j,i);
  });
  CleanRays();
  CleanClosureFill();
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
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
  CleanFor("cln_raymom", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = cs(m,CLN_J,k,j,i);
    for (int n = 0; n < CLN_N; ++n) {   // (0: a slot the staggered layout drops)
      dst(m,n,k,j,i) = (n >= cs.extent_int(1)) ? 0.0 :
                       ((CleanCellMoment(n) && n != CLN_J) ? jr : 1.0)*cs(m,n,k,j,i);
    }
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
  // (clean_extrap_share: writing cln_mc overwrites the moments at t^{n-1})
  if (cln_extrap == 3 && t2_solve == M1_T2S_STAGE1) {
    if (cln_mshare) cln_have_p = false;
    CleanRayMoments(cln_mc, false);
    cln_tm = time + dt;
    cln_have_c = true;
    cln_cand_be = false;
    CleanRingJump();
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
  if (cln_mshare) cln_have_p = false;
  CleanRayMoments(cln_mc);
  cln_tm = time + dt;
  cln_have_c = true;
  cln_cand_be = (t2_solve != M1_T2S_STAGE2);
  CleanRingJump();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRingSolveSet
//! \brief clean_rot_ring = jump: the swept set whose intensities give the face half
//! ranges of the solve is the set of the step's rotation block;
//! the step after a block boundary reads the next set of the boundary build

void RadiationM1::CleanRingSolveSet() {
  if (cln_cssc) return;   // (the sc rays on the cubed sphere: CleanScUnion2's two sets)
  const int64_t bs = VetGdBlock(pmy_pack->pmesh->ncycle);
  const int ns = static_cast<int>(cln_fvdm.size());
  for (int q = 0; q < ns && ns > 1; ++q) {
    if (cln_fvblk[q] == bs) Kokkos::deep_copy(cln_fvd, cln_fvdm[q]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRingJump
//! \brief clean_rot_ring = jump: the moments of a step that closes a rotation block
//! become those of the next block's window at its end (the closure jumps there), and the
//! next step does not extrapolate through the jump

void RadiationM1::CleanRingJump() {
  if (!cln_rjump) return;
  Kokkos::deep_copy(cln_mc, cln_rnew);
  cln_cand_be = true;
  cln_rjump = false;
}
#endif  // M1_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn bool CleanRealizable
//! \brief whether ray moments q (CLN_* layout, not normalised, J > 0) are those of a
//! non-negative intensity: |H| <= J and K - H H^T/J positive semidefinite (to round-off)

KOKKOS_INLINE_FUNCTION
bool CleanRealizable(const CleanArgs &x, const Real *q, const int m, const int k,
                     const int j, const int i) {
  const Real ij = 1.0/q[CLN_J];
  Real s[3], h[3], a[3][3];
  for (int b = 0; b < 3; ++b) s[b] = (q[CLN_HP+b] + q[CLN_HM+b])*ij;
  CleanRayFlux(x, m, k, j, i, s, h);
  a[0][0] = q[CLN_D11]*ij;
  a[1][1] = q[CLN_D11+1]*ij;
  a[2][2] = q[CLN_D11+2]*ij;
  a[0][1] = a[1][0] = q[CLN_D11+3]*ij;
  a[0][2] = a[2][0] = q[CLN_D11+4]*ij;
  a[1][2] = a[2][1] = q[CLN_D11+5]*ij;
  for (int b = 0; b < 3; ++b) {
    for (int c = 0; c < 3; ++c) a[b][c] -= h[b]*h[c];
  }
  // all principal minors of D - h h^T non-negative
  const Real eps = 1.0e-12;
  const Real m01 = a[0][0]*a[1][1] - a[0][1]*a[0][1];
  const Real m02 = a[0][0]*a[2][2] - a[0][2]*a[0][2];
  const Real m12 = a[1][1]*a[2][2] - a[1][2]*a[1][2];
  const Real det = a[0][0]*m12 - a[0][1]*(a[0][1]*a[2][2] - a[1][2]*a[0][2])
                   + a[0][2]*(a[0][1]*a[1][2] - a[1][1]*a[0][2]);
  return (h[0]*h[0] + h[1]*h[1] + h[2]*h[2] <= 1.0 + eps) &&
         (a[0][0] >= -eps) && (a[1][1] >= -eps) && (a[2][2] >= -eps) &&
         (m01 >= -eps) && (m02 >= -eps) && (m12 >= -eps) && (det >= -eps);
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

#if M1_TU_PART(0)
void RadiationM1::CleanClosureExtrap() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  cln_nguard = 0;
  if (cln_extrap == 3 && t2_solve == M1_T2S_STAGE2) return;

  // a new step: the end of the last step becomes t^n
  if (cln_have_c && cln_tm == time) {
    std::swap(cln_mp, cln_mn);
    if (cln_mshare) {
      cln_mc = cln_mp;
    } else {
      std::swap(cln_mn, cln_mc);
    }
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
  CleanArgs x = CleanMakeArgs();
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
    const bool ok = line && ma(m,CLN_J,k,j,i) > 0.0 && q[CLN_J] > 0.0 &&
                    CleanRealizable(x, q, m, k, j, i);
    if (!ok) {
      if (line && mn(m,CLN_J,k,j,i) > 0.0) ++ng;
      // (the face sums of a MeshBlock face on the cubed sphere, held on both sides, stay
      // extrapolated alike)
      for (int n = 0; n < CLN_N; ++n) {
        if (!x.cub || CleanCellMoment(n)) q[n] = mn(m,n,k,j,i);
      }
    }
    CleanStoreMoments(cs, valid, q, m, k, j, i);
  }, Kokkos::Sum<int>(nguard));
  cln_nguard = static_cast<int>(CleanReduce(nguard, false));
  cln_nguard_sum += cln_nguard;
  if (line) cln_next_sum += 1.0;

  CleanClosureFill();
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
}
#endif  // M1_TU_PART(0)

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
//! of U^{n+1}.  Without a valid K1 (first step, backward Euler) the rates of 1. and 2.
//! are zero, but the exchange 3. still relaxes the gas, as old production's relaxed ray
//! source does in every step (a hot gas the step cools is no source of the closure).
//! With use_k1 = false the formal solution sees the current state itself.

#if M1_TU_PART(0)
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
  const bool exch = use_k1 && cln_pred;
  const bool rates = exch && t2_ok && (t2k1.extent_int(0) > nmb1);

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
  const bool ps = cln_psmooth;
  CleanArgs x = CleanMakeArgs();

  CleanFor("cln_predict", nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // K1 is stored on active cells only; ghost cells take the explicit state
    const bool active = (k >= ks && k <= ke && j >= js && j <= je && i >= is && i <= ie);
    const bool use_rates = rates && active;
    const bool relax = exch && feedback && active;
    const Real e0 = u0_(m,M1_E,k,j,i);
    Real e = e0, d = 1.0, t = 0.0;
    if (fluid) {
      d = uh(m,IDN,k,j,i);
      Real mom[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
      const Real ekin = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
      const Real eint_h = fmax(CleanGasEint(x, uh, ekin, m, k, j, i), 0.0);
      t = eos.Temperature(d, eint_h);
      if (use_rates && !feedback) {
        // gas held fixed: E advances at its own rate, which includes the exchange
        const Real de = dt*k1(m,M1_T2_E,k,j,i);
        e = CleanPredictE(e0, de, ps);
      } else if (relax) {
        // radiation transport over dt (a valid K1), less the work of the radiation force
        // on the gas, then the implicit local exchange of (eint, E) over dt
        if (use_rates) {
          for (int n = 0; n < 3; ++n) mom[n] += dt*k1(m,M1_T2_M1+n,k,j,i);
          const Real ekin_new = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
          const Real de = dt*(k1(m,M1_T2_E,k,j,i) + k1(m,M1_T2_EN,k,j,i))
                          - (ekin_new - ekin);
          e = CleanPredictE(e0, de, ps);
        }
        // E (1 + h c kE) = e + h c kP a T^4,  eint(T) + E = eint_h + e
        const Real hce = dt*cl*opac_(m,M1_OP_E,k,j,i);
        const Real hcp = dt*cl*opac_(m,M1_OP_P,k,j,i)*arad_;
        t = CleanExchangeTemperature(eos, d, hcp/(1.0 + hce), eint_h + e*hce/(1.0 + hce),
                                     t);
        e = (e + hcp*SQR(SQR(t)))/(1.0 + hce);
      }
    } else if (use_rates) {
      const Real de = dt*k1(m,M1_T2_E,k,j,i);
      e = CleanPredictE(e0, de, ps);
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
//! On the cubed sphere D and the stress columns take Cartesian components about
//! r_hat (that of the frame of y) and the mean normals n_d of the lateral families:
//! k+^r = k_rr r_hat, k+^d = k r_hat + D_tt n_d/2.

void RadiationM1::CleanClosureSphere() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real rs = cln_sph_r;
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  auto x1v = pmy_pack->pcoord->x1v;
  auto cs = cln_s;
  CleanArgs x = CleanMakeArgs();

  CleanFor("cln_sphere", nmb1, ks, ke, js, je, is, ie,
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
    cs(m,CLN_C1,k,j,i) = c1;
    cs(m,CLN_J,k,j,i) = 1.0;
    const bool kp = cs.extent_int(1) > CLN_KP;
    if (kp) {
      for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = 0.0;
      cs(m,CLN_KP,k,j,i) = krr;
      cs(m,CLN_KP+3,k,j,i) = kr;
      cs(m,CLN_KP+4,k,j,i) = 0.25*(1.0 - drr);
      cs(m,CLN_KP+6,k,j,i) = kr;
      cs(m,CLN_KP+8,k,j,i) = 0.25*(1.0 - drr);
    }
    if (!cub) return;
    const Real dt = 0.5*(1.0 - drr);
    Real nd[3][3];
    CleanCsFamily(x, m, k, j, i, nd);
    const int ab[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    for (int a = 0; a < 6; ++a) {
      cs(m,CLN_D11+a,k,j,i) = (drr - dt)*nd[0][ab[a][0]]*nd[0][ab[a][1]] +
                              ((a < 3) ? dt : 0.0);
    }
    for (int b = 0; b < 3 && kp; ++b) {
      cs(m,CLN_KP+b,k,j,i) = krr*nd[0][b];
      for (int d = 1; d < 3; ++d) {
        cs(m,CLN_KP+3*d+b,k,j,i) = kr*nd[0][b] + 0.5*dt*nd[d][b];
      }
    }
  });
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRays
//! \brief the formal solution of the state in iw (E, T, F; chi from opac) and its
//! closure: vet_sc on a Cartesian mesh, vet_gd from the first shell on spherical polar.
//! The source function is the state's own (no relaxed source or relaxed scattering J,
//! which would make the closure depend on dt).

void RadiationM1::CleanRays() {
  MemcutPoison(2);   // diagnostic (env MEMCUT_POISON)
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
  // (the closure's time: a step at least one rotation block long, the steady start's
  // large steps, takes the sets of its start, as VetGdBlock)
  const Real dt = pmy_pack->pmesh->dt, tc = pmy_pack->pmesh->time +
                  ((vgd_rot_dt > 0.0 && dt >= vgd_rot_dt) ? 0.0 : dt);
  par_for("cln_gd_chi", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
  });
  if (cln_ring) {   // (the sc rays: each set's sweep in CleanRingMoments)
    VetGdSource();
    CleanFvSets(tc, pmy_pack->pmesh->ncycle);
  } else {
    VetGdBuild();
    if (cln_scmem) VetGdTopCopy();
  }
  if (cln_ring) {
    CleanRingMoments();
  } else {
    CleanClosureFromGd();
    // (clean_rays = sc on the cubed sphere: also the intensities cln_fvi of the face half
    // ranges, CleanFaceHalf)
    if (cln_consist || cln_cssc) CleanGdFaceSums();
    if (cln_cssc) CleanFaceHalf(vgd_dir, 0, true, 0);
    if (cln_cssc) cln_fvi = DvceArray5D<Real>();   // (build scratch)
  }
  MemcutPoison(3);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanScUnion2
//! \brief clean_rot_ring, the sc rays on the cubed sphere: the union the face half ranges
//! read (cln_fvd, cln_fvn), the sets b and b+1 of the step in their cln_fvi slots
//! (blk mod 2) with the cross-fade weights 1 - s and s: two sets of intensities resident

void RadiationM1::CleanScUnion2() {
  const int n = vgd_n;
  const int64_t b = cln_rbld;
  if (cln_fvd.extent_int(0) != 2*n) cln_fvd = DvceArray2D<Real>("cln_fvd", 2*n, 4);
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  for (int64_t blk = b; blk <= b + 1; ++blk) {
    const int q = static_cast<int>(blk & 1);
    const std::vector<double> d = VgdRotate(VetGdAngleBlk(blk));
    const Real w = (blk > b) ? cln_rs : 1.0 - cln_rs;
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  cln_fvn = 2*n;
}
#endif  // M1_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn bool CleanRingMix
//! \brief clean_rot_ring: the raw moments q of the union of the current set (slot 0 of
//! rl, CleanRingL), the next one (slot qn, weight s) and the older sets r (weights 1,
//! the last 1 - s) from their normalised moments (rl, the offsets rr): J-weighted with
//! the ratios r_r = J_{b-r}/J_b, J their weighted mean.  Invalid offsets (r_r < 0) and a
//! current set without rays (false) are left out, as the production ring (RotAvgApply)

struct CleanRingL {
  DvceArray5D<Real> rl, cs;
  Kokkos::Array<int, CLN_N> fm;
  int nf, ol, ks, js, is;
  // (slot o's normalised moment n at the interior cell (k,j,i): the swept set ol from
  // cln_s, where its sweep left it (the last), as the raw copy would hold it)
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int o, const int n, const int k, const int j,
                  const int i) const {
    if (o == ol) {
      const Real jr = cs(m,CLN_J,k+ks,j+js,i+is), ij = (jr > 0.0) ? 1.0/jr : 0.0;
      return (CleanCellMoment(n) ? 1.0 : ij)*cs(m,n,k+ks,j+js,i+is);
    }
    return rl(m,((ol < 0) ? o*nf : 0)+fm[n],k,j,i);
  }
};

KOKKOS_INLINE_FUNCTION
bool CleanRingMix(const CleanRingL &rl, const DvceArray5D<Real> &rr,
                  const Kokkos::Array<int, CLN_N> &fr, const int nfr, const int qn,
                  const int nr, const Real s, const int m, const int k, const int j,
                  const int i, Real *q) {
  const Real jb = rl(m,0,CLN_J,k,j,i);
  if (!(jb > 0.0)) return false;
  Real a[RadiationM1::M1_RA_NF];
  const Real wn = s*rl(m,qn,CLN_J,k,j,i)/jb;
  Real na = 1.0 + s, nw = 1.0 + wn;
  for (int r = 0; r < nr; ++r) {
    const Real rj = rr(m,r*nfr+fr[int{CLN_J}],k,j,i), w = (r == nr - 1) ? 1.0 - s : 1.0;
    a[r] = (rj >= 0.0) ? w*rj : 0.0;
    na += (rj >= 0.0) ? w : 0.0;
    nw += a[r];
  }
  q[CLN_J] = jb*nw/na;
  for (int n = 0; n < CLN_N; ++n) {
    if (n == CLN_J) continue;
    if (rl.fm[n] < 0) {q[n] = 0.0; continue;}
    const Real gb = rl(m,0,n,k,j,i);
    if (fr[n] < 0) {   // (not on the ring: the swept sets', as the fallback below)
      q[n] = (1.0 - s)*(gb*jb) + s*(rl(m,qn,n,k,j,i)*rl(m,qn,CLN_J,k,j,i));
      continue;
    }
    Real v = gb + wn*rl(m,qn,n,k,j,i);
    for (int r = 0; r < nr; ++r) v += a[r]*(gb + rr(m,r*nfr+fr[n],k,j,i));
    q[n] = v/nw*q[CLN_J];
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanRingEma
//! \brief clean_rot_ring with vet_rotate_avg_ema (production's form, RotAvgEma*,
//! RotXfBlend): the raw moments q of (1 - s) A + s T in the normalised fields and J. A =
//! the current set with the EMA offsets (e = rr's first nfr, its J entry the ratio r =
//! J_avg/J_b), T = the next set with the EMA the coming boundary makes of this state's
//! Delta (RotAvgEmaShift): aw = max(1/N, a_ch, 1/(n + 1)), r' = (1 - aw) r/R + aw, e' =
//! (1 - aw)(r/R)/r' (e - D), R = J_{b+1}/J_b; a_ch is this build's (pd's entry nfr), n
//! rr's entry nfr + 2 (0: no history, A and T the sets themselves)

KOKKOS_INLINE_FUNCTION
bool CleanRingEma(const CleanRingL &rl, const DvceArray5D<Real> &rr,
                  const DvceArray5D<Real> &pd, const Kokkos::Array<int, CLN_N> &fr,
                  const int nfr, const int qn, const Real s, const Real a0, const int m,
                  const int k, const int j, const int i, Real *q) {
  const Real jb = rl(m,0,CLN_J,k,j,i), jn = rl(m,qn,CLN_J,k,j,i);
  if (!(jb > 0.0)) return false;
  const Real nset = rr(m,nfr+2,k,j,i), r = rr(m,fr[int{CLN_J}],k,j,i);
  const bool hist = nset > 0.5 && r > 0.0, sh = hist && jn > 0.0;
  Real jt = jn, fac = 0.0;
  if (sh) {
    const Real aw = (a0 < 0.0) ? fmax(0.0, fmin(pd(m,nfr,k,j,i), 1.0)) :
                     fmax(fmax(a0, fmin(pd(m,nfr,k,j,i), 1.0)), 1.0/(nset + 1.0));
    const Real rq = r*jb/jn, r1 = (1.0 - aw)*rq + aw;
    fac = (1.0 - aw)*rq/r1;
    jt = jn*r1;
  }
  q[CLN_J] = (1.0 - s)*jb*(hist ? r : 1.0) + s*jt;
  for (int n = 0; n < CLN_N; ++n) {
    if (n == CLN_J) continue;
    if (rl.fm[n] < 0) {q[n] = 0.0; continue;}
    const Real gb = rl(m,0,n,k,j,i), gn = rl(m,qn,n,k,j,i);
    if (fr[n] < 0) {   // (not on the ring: the swept sets', as CleanRingMix)
      q[n] = (1.0 - s)*(gb*jb) + s*(gn*jn);
      continue;
    }
    const Real e = rr(m,fr[n],k,j,i);
    const Real ga = gb + (hist ? e : 0.0), gt = sh ? gn + fac*(e - (gn - gb)) : gn;
    q[n] = ((1.0 - s)*ga + s*gt)*q[CLN_J];
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanRingFields
//! \brief clean_rot_ring: the fields the offsets and the Delta carry, fm(CLN index) =
//! their slot or -1, and their number.  The staggered layout leaves out those none of its
//! solves reads: those cln_s does not hold (CleanNs: CLN_KP, CLN_F1, off the cubed sphere
//! CLN_FH) and the stress parts of the face sums CLN_FX, CLN_F2, CLN_F3 and (cubed
//! sphere) CLN_FH: only their energy flux is read; on the cubed sphere also CLN_G (the
//! sp F solve reads it).  The mixture gives 0 there (CleanRingMix).  ring: those of the
//! offsets and the Delta (with clean_ring_fields = closure the cell moments and CLN_JR;
//! the others are mixed from the swept sets alone).  file: the fields of the restart
//! files' layout (CleanRingRstCopy): up to version 6 it keeps the stress parts and CLN_G
//! on the cubed sphere (0 in a file written since), from version 7 the ring's own

#if M1_TU_PART(0)
int RadiationM1::CleanRingFields(int *fm, const bool ring, const bool file) const {
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  const int ns = CleanNs(cub);   // (the slots of cln_s)
  int nf = 0;
  for (int n = 0; n < CLN_N; ++n) {
    const bool fs = (n >= CLN_F2 && n < CLN_G && (n - CLN_F2)%4 != 0) ||   // (stress
                    (n >= CLN_FH && n < CLN_F1 && (n - CLN_FH)%4 != 0);    // parts)
    const bool un = fs || (cub && (n == CLN_G || n == CLN_G + 1));
    const bool cm = !(ring && cln_rcell) || CleanCellMoment(n) || n == CLN_JR;
    const bool f6 = file && t2_rstv < 7;
    fm[n] = ((n < ns && !(f6 ? fs && !cub : un)) && cm) ? nf++ : -1;
  }
  return nf;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRingShift
//! \brief clean_rot_ring: the first build of block b shifts the ring by the pending
//! Delta (that of the boundary of block b when it closed the block): d_1 = -Delta,
//! d_j -= Delta, the J ratios r_j' = r_{j-1}/R, r_1' = 1/R

void RadiationM1::CleanRingShift(const int64_t b) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nr = vra_k - 1;
  Kokkos::Array<int, CLN_N> fm;
  const int nf = CleanRingFields(fm.data(), true);
  const Real tk = vgd_rot_t0 + static_cast<Real>(b)*vgd_rot_dt;
  if (std::abs(cln_rpt - tk) > 1.0e-9*vgd_rot_dt) cln_roff_sum += 1.0;
  vra_gen += 1;
  auto rr = cln_rr;
  auto pd = cln_rpd;
  if (vra_ema) {
    // (RotAvgEmaShift: the EMA moves on with the Delta and the change weight a_ch of the
    // boundary; the aux := the next set there.  No history or no Delta: the current set
    // alone)
    const Real a0 = vra_ema_wb ? -1.0 : 1.0/vra_ema_n;   // (< 0 = a_ch alone)
    par_for("cln_ring_eshift", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0,
            id.nx3 - 1, 0, id.nx2 - 1, 0, id.nx1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int fj = fm[int{CLN_J}];
      const Real R = pd(m,fj,k,j,i), nset = rr(m,nf+2,k,j,i), r = rr(m,fj,k,j,i);
      rr(m,nf,k,j,i) = pd(m,nf+1,k,j,i);
      rr(m,nf+1,k,j,i) = pd(m,nf+2,k,j,i);
      if (!(R > 0.0) || !(nset > 0.5) || !(r > 0.0)) {
        for (int f = 0; f < nf; ++f) rr(m,f,k,j,i) = (f == fj) ? 1.0 : 0.0;
        rr(m,nf+2,k,j,i) = (R > 0.0) ? 1.0 : 0.0;
        return;
      }
      const Real aw = (a0 < 0.0) ? fmax(0.0, fmin(pd(m,nf,k,j,i), 1.0)) :
                      fmax(fmax(a0, fmin(pd(m,nf,k,j,i), 1.0)), 1.0/(nset + 1.0));
      const Real r1 = (1.0 - aw)*r/R + aw, fac = (1.0 - aw)*(r/R)/r1;
      for (int f = 0; f < nf; ++f) {
        rr(m,f,k,j,i) = (f == fj) ? r1 : fac*(rr(m,f,k,j,i) - pd(m,f,k,j,i));
      }
      rr(m,nf+2,k,j,i) = (aw >= 1.0) ? 1.0 : fmin(nset + 1.0, 1.0e9);
    });
  } else {
  par_for("cln_ring_shift", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, id.nx3 - 1,
          0, id.nx2 - 1, 0, id.nx1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int fj = fm[int{CLN_J}];
    const Real rp = pd(m,fj,k,j,i);
    for (int f = 0; f < nf; ++f) {
      for (int r = nr - 1; r > 0; --r) {
        const Real o = rr(m,(r-1)*nf+f,k,j,i);
        rr(m,r*nf+f,k,j,i) = (f != fj) ? o - pd(m,f,k,j,i) :
                             ((rp > 0.0 && o >= 0.0) ? o/rp : -1.0);
      }
      rr(m,f,k,j,i) = (f != fj) ? -pd(m,f,k,j,i) : ((rp > 0.0) ? 1.0/rp : -1.0);
    }
  });
  }
  cln_rok = true;
  cln_rb = b;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRingPreShift
//! \brief dt_retry checkpoint (DtRetryQuiesce): with clean_rot_ring = live, the
//! cross-fade and vet_gd_rotate_align the build of a step lies in one block for every dt
//! of its retries (a step ends on each boundary): a step that starts on the boundary of
//! block cln_rb shifts the ring here, before the checkpoint, as its first build would.
//! The offsets then do not change within the cycle and the pending Delta is read by no
//! build of it: the checkpoint skips both (cln_rfix).  Otherwise the build shifts and
//! both are copied

int64_t RadiationM1::CleanRingBlk(const Real dt) const {
  const Real tau = (pmy_pack->pmesh->time + dt - vgd_rot_t0)/vgd_rot_dt;
  const int64_t b = static_cast<int64_t>(std::floor(tau));
  return (tau - static_cast<Real>(b) < 1.0e-9) ? b - 1 : b;   // (CleanFvSets)
}

// (live with the cross-fade: every build sweeps the sets b and b+1)
bool RadiationM1::CleanRingTwo() const {
  return clean && cln_ring == 1 && vgd_xf && !(vra_on && !g_vra_adv);
}

// (the sc rays on the staggered cubed sphere: their intensities are a build's scratch,
// folded into the face half ranges cln_fh; those are rebuilt before every solve where
// every solve builds the closure first and every build sweeps all the sets they read,
// else they are state: the dt_retry checkpoint and the restart files keep them)
bool RadiationM1::CleanFhBuilt() const {
  return cln_cssc && cln_extrap == 0 && !cln_closit && !cln_edd &&
         !cln_sphere && (!cln_ring || CleanRingTwo());
}

void RadiationM1::CleanRingPreShift() {
  const Real dt = pmy_pack->pmesh->dt;
  cln_rfix = CleanRingTwo() && vgd_rot_dt > 0.0 && vgd_rot_align && dt < vgd_rot_dt &&
             CleanRingBlk(dt) == CleanRingBlk(1.0e-6*dt);   // (any retry's dt)
  if (cln_rfix && cln_rok && cln_rpb == cln_rb && CleanRingBlk(dt) == cln_rb + 1) {
    CleanRingShift(cln_rb + 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::CleanRingRst*
//! \brief clean_rot_ring in the restart files (version 4): the offsets cln_rr, then the
//! Delta cln_rpd, of the fields the ring carries on the interior cells, in their own
//! layout flat per MeshBlock (the fields of the file layout, CleanRingFields: 0 in those
//! cln_rr does not hold; with vet_rotate_avg_ema each followed by its 3 entries,
//! CleanRingNr / Nd), dense in CleanRingRstN() ghosted slabs (the rest of the last one
//! 0).  CleanRingRstCopy packs (pk; 0 without have) or unpacks the slabs c ..
//! c + CleanRingRstN() - 1 of a (c < 0: a's channel 0 is slab -c) of MeshBlocks m0 ..

int RadiationM1::CleanRingRstN() const {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N];
  const int nff = CleanRingFields(fm, true, true);
  const int64_t nv = static_cast<int64_t>(CleanRingNr(nff) + CleanRingNd(nff))*id.nx1*
                     id.nx2*id.nx3;
  const int64_t ns = static_cast<int64_t>(id.nx1 + 2*id.ng)*
                     ((id.nx2 > 1) ? id.nx2 + 2*id.ng : 1)*
                     ((id.nx3 > 1) ? id.nx3 + 2*id.ng : 1);
  return static_cast<int>((nv + ns - 1)/ns);
}

void RadiationM1::CleanRingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb,
                                   const int m0, const bool pk, const bool have) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N], ff[CLN_N];
  const int nf = CleanRingFields(fm, true), nff = CleanRingFields(ff, true, true);
  Kokkos::Array<int, CLN_N> mp;   // (the file's field slot -> cln_rr's, -1: not held)
  for (int n = 0; n < CLN_N; ++n) {
    if (ff[n] >= 0) mp[ff[n]] = fm[n];
  }
  if (!pk) {   // (the writers of cln_rr bump vra_gen, DtRetry)
    vra_gen += 1;
    if (cln_rr.extent_int(0) != nmb) {
      Kokkos::realloc(cln_rr, nmb, CleanRingNr(nf), id.nx3, id.nx2, id.nx1);
      Kokkos::realloc(cln_rpd, nmb, CleanRingNd(nf), id.nx3, id.nx2, id.nx1);
    }
  }
  // (the field groups of cln_rr: the K - 1 offsets, or the EMA's 1; the entries after
  // the groups are copied as they are)
  const int gr = vra_ema ? 1 : vra_k - 1;
  const int64_t ni = static_cast<int64_t>(id.nx1)*id.nx2*id.nx3;
  const int64_t nr = CleanRingNr(nf)*ni, ndm = CleanRingNd(nf)*ni;
  const int64_t nro = CleanRingNr(nff)*ni, nd = nro + CleanRingNd(nff)*ni;
  const int64_t ns = static_cast<int64_t>(a.extent(2))*a.extent(3)*a.extent(4);
  const int64_t q0 = std::max(c, 0) - c, q1 = std::min<int64_t>(c + CleanRingRstN(),
                                                                 a.extent(1)) - c;
  if (q1 <= q0) return;
  const int64_t nq = (q1 - q0)*ns, sa = a.stride(0);
  Real *pa = a.data();
  Real *pr = (have && cln_rr.extent_int(0) >= m0 + nmb) ? cln_rr.data() : nullptr;
  Real *pd = (pr != nullptr) ? cln_rpd.data() : nullptr;
  Kokkos::parallel_for("m1_ring_rst",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>(DevExeSpace(), 0, nmb*nq),
  KOKKOS_LAMBDA(const int64_t p) {
    const int64_t m = p/nq, q = q0*ns + p - m*nq;
    Real *x = pa + m*sa + c*ns + q;
    Real *y = nullptr;
    if (pr != nullptr && q < nd) {
      const bool o = (q < nro);
      const int64_t u = o ? q : q - nro, ic = u%ni;
      const int g = static_cast<int>(u/ni), ng = o ? gr : 1;
      const int f = (g >= ng*nff) ? g - ng*(nff - nf) :
                    ((mp[g%nff] < 0) ? -1 : (g/nff)*nf + mp[g%nff]);
      if (f >= 0) y = (o ? pr + (m + m0)*nr : pd + (m + m0)*ndm) + f*ni + ic;
    }
    if (pk) {
      *x = (y != nullptr) ? *y : 0.0;
    } else if (y != nullptr) {
      *y = *x;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn RadiationM1::CleanFvRst*
//! \brief the sc rays' intensities on the cubed sphere in the restart files (version 5):
//! a slab of their count nv, then cln_fvi's interior cells (v, k, j, i), dense in the
//! ghosted slabs as CleanRingRstCopy.  The ghosts are those of the cell-centred exchange
//! of the interior (CleanGdFaceSums: CleanGhostFill), filled again after a restart has
//! read the interior (cln_fvgf, CleanRstFaces).  Unpacking (c <= 0: a's channel 0 is
//! record slab -c) slab 0 allocates cln_fvi

int RadiationM1::CleanFvRstN() const {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int64_t nv = static_cast<int64_t>(cln_fvi.extent_int(1))*id.nx1*id.nx2*id.nx3;
  const int64_t ns = static_cast<int64_t>(cln_fvi.extent_int(2))*cln_fvi.extent_int(3)*
                     cln_fvi.extent_int(4);
  return 1 + static_cast<int>((nv + ns - 1)/ns);
}

void RadiationM1::CleanFvRstCopy(DvceArray5D<Real> &a, const int c, const int nmb,
                                 const int m0, const bool pk) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  auto AL = Kokkos::ALL;
  if (!pk && c == 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                 Kokkos::subview(a, 0, 0, 0, 0, 0));
    Kokkos::realloc(cln_fvi, nmb, static_cast<int>(std::llround(h())), a.extent_int(2),
                    a.extent_int(3), a.extent_int(4));
    cln_fvgf = true;
    return;
  }
  if (pk) {
    Kokkos::deep_copy(Kokkos::subview(a, std::make_pair(0, nmb), c, AL, AL, AL),
                      static_cast<Real>(cln_fvi.extent_int(1)));
  }
  const int64_t ni = id.nx1, nj = id.nx2, nk = id.nx3, nc = ni*nj*nk;
  const int64_t nd = cln_fvi.extent_int(1)*nc;
  const int64_t ns = static_cast<int64_t>(a.extent(2))*a.extent(3)*a.extent(4);
  const int c1 = c + 1;   // (the data behind the count)
  const int64_t q0 = std::max(c1, 0) - c1, q1 = std::min<int64_t>(c + CleanFvRstN(),
                                                                   a.extent(1)) - c1;
  if (q1 <= q0) return;
  const int64_t nq = (q1 - q0)*ns, sa = a.stride(0);
  const int is = id.is, js = id.js, ks = id.ks;
  Real *pa = a.data();
  auto fi = cln_fvi;
  Kokkos::parallel_for("m1_fv_rst",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>(DevExeSpace(), 0, nmb*nq),
  KOKKOS_LAMBDA(const int64_t p) {
    const int64_t m = p/nq, q = q0*ns + p - m*nq;
    Real *x = pa + m*sa + c1*ns + q;
    const int v = static_cast<int>(q/nc), r = static_cast<int>(q - v*nc);
    const int i = is + r%ni, j = js + (r/ni)%nj, k = ks + r/(ni*nj);
    if (pk) {
      *x = (q < nd) ? fi(m + m0,v,k,j,i) : 0.0;
    } else if (q < nd) {
      fi(m,v,k,j,i) = *x;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRingMoments
//! \brief clean_rot_ring: the closure of the rotating union of vet_rotate_avg = K sets
//! with the cross-fade from a sweep of the current set b and the next set b+1 alone.
//! The union's moments (CLN_*, not normalised; the face sums included) are
//!   M = (1 - s)/K M_{b-K+1} + 1/K sum_{j=1}^{K-2} M_{b-j} + 1/K M_b + s/K M_{b+1}
//! with the older sets taken from the normalised moments g = M/J of the current one and
//! the offsets d_j = g_{b-j} - g_b, r_j = J_{b-j}/J_b (CleanRingMix) measured on the same
//! state, as production's ring (rad_m1_rotavg.cpp): the ray pattern is carried, the
//! current set's J and moments stay live, so the mixture follows the state's scale.
//! The build that closes block b (on its boundary t_b+1,
//! s = 1; vet_gd_rotate_align puts a step end there on every dt) gives Delta = M_{b+1} -
//! M_b, and the first build of block b+1 shifts the ring: d_1 = -Delta, d_j -= Delta.
//! The offsets are therefore piecewise constant in t with steps at the block boundaries:
//! the closure is a function of t, the state and the states at the boundaries, continuous
//! in t, so the integrator keeps its order in dt; against the union it lags by the
//! change of the ray-pattern differences over at most K rotation blocks.  A fresh start
//! (or a missed boundary) sweeps the whole union once to fill the ring.  A cell whose
//! mixture is not that of a non-negative intensity takes the two swept sets' (1 - s) M_b
//! + s M_{b+1} (the cell moments only on the cubed sphere, whose face sums are held on
//! both sides of a MeshBlock face).  Linear in the moments, so every face sum stays
//! single-valued and conservative.  The sc rays (clean_rays = sc, spherical polar) sweep
//! each set alone here, in its own tables (VgdXfTables keeps two).

void RadiationM1::CleanRingMoments() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nk = CleanRingK(), nr = nk - 1, ns = cln_rns;
  const int64_t b = cln_rbld;
  const Real s = cln_rs;
  const int c3 = ke - ks + 1, c2 = je - js + 1, c1 = ie - is + 1;
  // (the swept sets' moments: the current set in cln_rl slot 0, the next in slot 1; the
  // older ones of a fill straight into their offsets' slots of cln_rr.  All on the
  // interior cells; the offsets and the Delta only of the fields fm carries)
  // (fr: the fields of the offsets and the Delta, nfr of them)
  Kokkos::Array<int, CLN_N> fm, fr;
  const int nf = CleanRingFields(fm.data()), nfr = CleanRingFields(fr.data(), true);
  // (the sweep order: the sc rays on the cubed sphere the older sets in block order,
  // then the live b, b+1 in their slot order blk mod 2, as CleanFaceHalf folds the union
  // (the sets' sweeps are independent); elsewhere the slot order (on sp the block order
  // breaks the bitwise restart of q_def, reflecting top).  The set swept last keeps its
  // moments in cln_s: the live ring's rl holds the other one, ol the slot of the last
  // (0: b, 1: b + 1), or both where an older set of a fill is swept last)
  std::vector<int> ord(ns);
  for (int q = 0; q < ns; ++q) ord[q] = q;
  auto okey = [&](int u) {
    const int64_t bu = cln_fvblk[u];
    return (bu >= b) ? std::make_pair(int64_t{1}, bu & 1)
                     : std::make_pair(int64_t{0}, bu);
  };
  if (cln_cssc) {
    std::sort(ord.begin(), ord.end(), [&](int u, int v) {return okey(u) < okey(v);});
  }
  const int64_t pl = cln_fvblk[ord[ns - 1]] - (b - nk + 1);
  const int ol = (cln_ring != 1) ? -1 : ((pl == nr) ? 0 : ((pl == nk) ? 1 : -1));
  // (over the stencil coefficients cln_csy where they fit: every pass's
  // CleanPreconditionerSetup rewrites them before they are read, dead across a build)
  const int nrl = ((ol < 0) ? 2 : 1)*nf;
  if (cln_csy.span() >= static_cast<std::size_t>(nmb)*nrl*c3*c2*c1) {
    cln_rl = DvceArray5D<Real>(cln_csy.data(), nmb, nrl, c3, c2, c1);
  } else if (cln_rl.extent_int(0) != nmb || cln_rl.extent_int(1) != nrl) {
    Kokkos::realloc(cln_rl, nmb, nrl, c3, c2, c1);
  }
  if (cln_rr.extent_int(0) != nmb) {   // (a restart has read them)
    vra_gen += 1;   // (cln_rr's writers: the dt_retry checkpoint skips a clean copy)
    Kokkos::realloc(cln_rr, nmb, CleanRingNr(nfr), c3, c2, c1);
    Kokkos::realloc(cln_rpd, nmb, CleanRingNd(nfr), c3, c2, c1);
    cln_rok = false;
  }
  // vet_rotate_avg_ema (CleanRingEma): rr = the EMA offsets, the aux J and D11 (the
  // current set at the last boundary) and the set count n; pd also a_ch and the next
  // set's J and D11
  const bool ema = vra_ema, chj = vra_ema_chj, estep = vra_ema_chstep;
  const Real a0 = ema ? (vra_ema_wb ? -1.0 : 1.0/vra_ema_n) : 1.0, eps = vra_ema_eps;
  std::vector<int> slot(nk + 2, -1);   // slot of block b - nk + 1 + p
  for (int q = 0; q < ns; ++q) {
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    if (p >= 0 && p <= nk) slot[p] = q;
  }
  bool full = true;
  for (int p = 0; p < nr; ++p) full = full && slot[p] >= 0;
  const bool shift = cln_rok && b == cln_rb + 1 && cln_rpb == cln_rb;
  const bool fill = (!cln_rok || b != cln_rb) && !shift && full;
  // the moments of each swept set (the ray_min test on the mixture below); the sc rays on
  // the cubed sphere: the union's directions first (the face half ranges, CleanFaceHalf)
  if (cln_cssc) CleanScUnion2();
  auto dfull = cln_fvd;
  bool fh0 = true;
  const Real ray_min = cln_raymin;
  cln_raymin = 0.0;
  // (clean_ring_jmin: the J ratios of the offsets and the Delta only where both sets'
  // J >= clean_ray_min E: at a fresh start's ring fill a set can see a source hole, J ~
  // 1e-200 in all sets, and J_s/J_b ~ 1e50 there came back times the next build's J)
  const Real jmin = cln_ringjmin ? ray_min : 0.0;
  auto iwr = iw;
  auto cs = cln_s;
  const CleanRingL rl{cln_rl, cln_s, fm, nf, ol, ks, js, is};
  auto rr = cln_rr;
  // clean_rot_fuse: the live ring's two sets of a build (b, b+1) in one fused sweep
  // (VetGdSweep with vgd_xff_now: each set's own tables and intensities, vgd_i and
  // vgd_i_xb), their lists made first; the loop below then takes each set's moments, its
  // tables and intensities brought back by VgdXfTables (not with a reflecting top)
  const bool fuse = cln_rfuse && !cln_scmem && ns == 2 && !fill && VetGdXffOk();
  if (cln_rfuse && !fuse && ns == 2 && !fill && !cln_rfuse_said &&
      global_variable::my_rank == 0) {
    std::cout << "<rad_m1> clean_rot_fuse: not taken (needs the sc rays, no "
              << "reflecting top and a fused-sweep configuration, VetGdXffOk): two "
              << "sweeps per build" << std::endl;
    cln_rfuse_said = true;
  }
  if (fuse) {
    const int qa = (cln_fvblk[0] < cln_fvblk[1]) ? 0 : 1, qb = 1 - qa;
    vgd_scut = vlat_icut;
    VgdXfTables(cln_fvkey[qb]);
    vgd_xff_prep = true;
    VetGdSweep();
    VgdXfTables(cln_fvkey[qa]);
    VetGdSweep();
    vgd_xff_prep = false;
    vgd_xff_now = true;
    VetGdSweep();
    vgd_xff_now = false;
    vgd_xff_n += 1.0;   // (the fused builds, reported with vet_gd_rotate_xfade_fuse)
  }
  for (int qo = 0; qo < ns; ++qo) {
    const int q = ord[qo];
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    const bool older = (p >= 0 && p < nr);
    const bool keep = (older && fill) || p == nr || p == nk;
    if (older && fill) vra_gen += 1;
    cln_fvd = cln_fvdm[q];
    // the sc rays: set q swept alone (its tables, its own weights)
    cln_rq = q;
    cln_rqm = static_cast<int>(cln_fvblk[q] & 1);
    VgdXfTables(cln_fvkey[q]);
    vgd_scut = vlat_icut;
    if (!fuse) VetGdSweep();
    if (cln_scmem) VetGdTopCopy();
    CleanClosureFromGd();
    if (cln_consist || cln_cssc) CleanGdFaceSums();
    if (cln_cssc && (p == nr || p == nk)) {
      CleanFaceHalf(dfull, cln_rqm*vgd_n, fh0, 0);
      fh0 = false;
    }
    const int sl = (p == nk) ? 1 : 0;   // (the slot of a live set)
    if (!keep || (!older && sl == ol)) continue;
    const int o = static_cast<int>(older ? (nr - 1 - p)*nfr : ((ol < 0) ? sl*nf : 0));
    auto ra = older ? rr : cln_rl;
    const auto fa = older ? fr : fm;
    par_for("cln_ring_raw", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jr = cs(m,CLN_J,k,j,i), ij = (jr > 0.0) ? 1.0/jr : 0.0;
      for (int n = 0; n < CLN_N; ++n) {
        if (fa[n] < 0) continue;
        ra(m,o+fa[n],k-ks,j-js,i-is) = (CleanCellMoment(n) ? 1.0 : ij)*cs(m,n,k,j,i);
      }
    });
  }
  cln_fvd = dfull;
  cln_raymin = ray_min;
  // (jump: the next set only at the build that closes the block; its weight is s = 0)
  const int qb = 0, qn = (slot[nk] >= 0) ? 1 : 0;
  auto pd = cln_rpd;
  if (!cln_rok || b != cln_rb) {
    const bool f0 = !cln_insteady;   // (a fill outside the steady start: below)
    if (shift) {
      CleanRingShift(b);
    } else if (full && ema) {
      // RotAvgEmaInit: the J-weighted equal mix of the current and the previous set (its
      // moments in rr's offsets, above); the aux := the current set.  Outside the steady
      // start the fill0 form (as the window's fill below): the current set alone, count 1
      par_for("cln_ring_einit", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i), jp = rr(m,fr[int{CLN_J}],k,j,i);
        const Real jm = jmin*iwr(m,M1_IW_EN,k0,j0,i0);
        const bool v = !f0 && jb > 0.0 && jb >= jm && jp > 0.0 && jp >= jm;
        const Real ws = v ? jp/jb : 1.0;
        for (int n = 0; n < CLN_N; ++n) {
          const int f = fr[n];
          if (f < 0) continue;
          rr(m,f,k,j,i) = (n == CLN_J) ? (v ? 0.5*(1.0 + ws) : 1.0) :
                          (v ? ws*(rr(m,f,k,j,i) - rl(m,qb,n,k,j,i))/(1.0 + ws) : 0.0);
        }
        rr(m,nfr,k,j,i) = jb;
        rr(m,nfr+1,k,j,i) = rl(m,qb,CLN_D11,k,j,i);
        rr(m,nfr+2,k,j,i) = v ? 2.0 : ((f0 && jb > 0.0 && jb >= jm) ? 1.0 : 0.0);
      });
    } else if (full) {
      // a fill: within CleanSteadyStart (a fresh start's ring, built at the initial state
      // and time) the offsets of the older sets from the current one (their moments,
      // above); elsewhere the older sets start as the current one (zero offsets, J ratio
      // 1) and enter one per block through the shifts: their sweep at a step's build
      // would take a state one step on, at a time that depends on dt, into the window for
      // K blocks, an O(dt) error (deck D with rotation 1st order, 2nd with this form)
      par_for("cln_ring_fill", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i);
        const Real jm = jmin*iwr(m,M1_IW_EN,k0,j0,i0);
        for (int r = 0; r < nr; ++r) {
          for (int n = 0; n < CLN_N; ++n) {
            const int f = fr[n];
            if (f < 0) continue;
            const Real o = rr(m,r*nfr+f,k,j,i);
            rr(m,r*nfr+f,k,j,i) = f0 ? ((n != CLN_J) ? 0.0 : 1.0) : (n != CLN_J) ?
                                  o - rl(m,qb,n,k,j,i) :
                                  ((jb > 0.0 && jb >= jm && o >= jm) ? o/jb : -1.0);
          }
        }
      });
    } else {
      // a boundary missed: this build takes the swept sets, the next one the union
      cln_roff_sum += 1.0;
      cln_rok = false;
    }
    if (shift || full) {
      cln_rok = true;
      cln_rb = b;
    }
  }
  // the pending Delta of this build
  if (slot[nk] >= 0) {
    par_for("cln_ring_delta", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
      const int k = k0 - ks, j = j0 - js, i = i0 - is;
      const Real jb = rl(m,qb,CLN_J,k,j,i), jn = rl(m,qn,CLN_J,k,j,i);
      const Real jm = jmin*iwr(m,M1_IW_EN,k0,j0,i0);
      for (int n = 0; n < CLN_N; ++n) {
        if (fr[n] < 0) continue;
        pd(m,fr[n],k,j,i) = (n != CLN_J) ? rl(m,qn,n,k,j,i) - rl(m,qb,n,k,j,i) :
                                ((jb > 0.0 && jb >= jm && jn >= jm) ? jn/jb : -1.0);
      }
      if (!ema) return;
      // (RotAvgEmaChange: a_ch of the boundary, the current set now against the aux,
      // c = max |ln J_b/J_aux|, |D11_b - D11_aux|; the next set's J and D11 for the aux)
      Real ach = -1.0;
      if (eps > 0.0) {
        const Real jo = rr(m,nfr,k,j,i);
        Real c = !chj ? 0.0 : ((jb > 0.0 && jo > 0.0) ? fabs(log(jb/jo)) : 1.0e30);
        c = fmax(c, fabs(rl(m,qb,CLN_D11,k,j,i) - rr(m,nfr+1,k,j,i)));
        ach = estep ? fmin(1.0, fmax(0.0, c/eps - 1.0)) : fmin(1.0, c/eps);
      }
      pd(m,nfr,k,j,i) = ach;
      pd(m,nfr+1,k,j,i) = jn;
      pd(m,nfr+2,k,j,i) = rl(m,qn,CLN_D11,k,j,i);
    });
    cln_rpb = b;
    cln_rpt = vgd_rot_t0 + (static_cast<Real>(b) + (cln_rlast ? 1.0 : s))*vgd_rot_dt;
  }
  // the mixture; jump, at the build that closes the block: also the moments of the next
  // block's window there (cln_rnew, the closure history of the next step)
  const bool ring = cln_rok;
  const bool nxt = ring && cln_rlast && slot[nk] >= 0;
  if (nxt && cln_rnew.extent_int(0) != nmb) {   // (the layout of cln_mc it becomes)
    Kokkos::realloc(cln_rnew, nmb, CLN_N, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
  }
  auto rw = cln_rnew;
  auto iw_ = iw;
  auto valid = cln_v;
  CleanArgs x = CleanMakeArgs();
  int nguard = 0;
  const int ni = ie - is + 1, nj = je - js + 1, nkk = ke - ks + 1;
  Kokkos::parallel_reduce("cln_ring_mix",
                          Kokkos::RangePolicy<DevExeSpace>(0, nmb*nkk*nj*ni),
  KOKKOS_LAMBDA(const int t, int &ng) {
    const int i = is + t%ni;
    const int j = js + (t/ni)%nj;
    const int k = ks + (t/(ni*nj))%nkk;
    const int m = t/(ni*nj*nkk);
    // (the swept sets' raw moments M = g J)
    const int kr = k - ks, jr = j - js, ir = i - is;
    // (slot o; the fields the ring does not carry: 0, none of the solves reads them)
    auto raw = [&](const int o, const int n) {
      const Real jo = rl(m,o,CLN_J,kr,jr,ir);
      return (n == CLN_J) ? jo : ((fm[n] < 0) ? 0.0 : rl(m,o,n,kr,jr,ir)*jo);
    };
    Real q[CLN_N];
    const bool mix = ring && (ema ? CleanRingEma(rl, rr, pd, fr, nfr, qn, s, a0, m, kr,
                                                 jr, ir, q)
                                  : CleanRingMix(rl, rr, fr, nfr, qn, nr, s, m, kr, jr,
                                                 ir, q));
    if (!mix) {
      for (int n = 0; n < CLN_N; ++n) q[n] = (1.0 - s)*raw(0, n) + s*raw(qn, n);
    }
    if (mix && !(q[CLN_J] > 0.0 && CleanRealizable(x, q, m, k, j, i))) {
      if ((1.0 - s)*raw(0, CLN_J) + s*raw(qn, CLN_J) > 0.0) ++ng;
      for (int n = 0; n < CLN_N; ++n) {
        if (!x.cub || CleanCellMoment(n)) q[n] = (1.0 - s)*raw(0, n) + s*raw(qn, n);
      }
    }
    if (q[CLN_J] < ray_min*iw_(m,M1_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    CleanStoreMoments(cs, valid, q, m, k, j, i);
    if (!nxt) return;
    // (the next block's window: the union at s = 1)
    if (!(CleanRingMix(rl, rr, fr, nfr, qn, nr, 1.0, m, kr, jr, ir, q) &&
          q[CLN_J] > 0.0 && CleanRealizable(x, q, m, k, j, i))) {
      if (rl(m,qn,CLN_J,kr,jr,ir) > 0.0) ++ng;
      for (int n = 0; n < CLN_N; ++n) {
        if (!x.cub || CleanCellMoment(n)) q[n] = raw(qn, n);
      }
    }
    if (q[CLN_J] < ray_min*iw_(m,M1_IW_EN,k,j,i)) {
      for (int n = 0; n < CLN_N; ++n) {
        if (CleanCellMoment(n)) q[n] = 0.0;
      }
    }
    for (int n = 0; n < CLN_N; ++n) rw(m,n,k,j,i) = q[n];
  }, Kokkos::Sum<int>(nguard));
  cln_rjump = nxt;
  cln_rguard_sum += CleanReduce(nguard, false);
  cln_rl = DvceArray5D<Real>();   // (build scratch: freed, every build writes it first)
  if (cln_cssc) cln_fvi = DvceArray5D<Real>();   // (and so cln_fvi)
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvSets
//! \brief clean_rot_ring: the swept directions (cln_fvd, cln_fvn) for the closure of
//! time tc (of the step of cycle cyc).  The rotating vet_gd set (vet_gd_rotate_dt |
//! _every) is swept as the union of the sets the closure averages, set by set in slots of
//! n directions with the set's weight folded into the quadrature weights, so that every
//! moment and face sum is the weighted mean of the sets' (linear in I: single-valued,
//! conservative and realisable).  With vet_rotate_avg = K these are the last K sets
//! b-K+1 .. b of block b, 1/K each; vet_gd_rotate_xfade adds the next set b+1 and weights
//! the oldest (1 - s)/K and the next s/K, s the fraction of block b passed at tc, so that
//! the closure is continuous in time.  A function of t and the state alone (a restart
//! rebuilds the same union).  Block blk takes slot blk mod (K + xfade).  Without
//! vet_rotate_advance the K sets repeat: the union is fixed and has nothing to fade.

void RadiationM1::CleanFvSets(const Real tc, const int cyc) {
  const int n = vgd_n;
  const Mesh *pm = pmy_pack->pmesh;
  const int nk = vra_on ? CleanRingK() : 1;
  const bool xf = vgd_xf && !(vra_on && !g_vra_adv);
  int64_t b = VetGdBlock(cyc);
  Real s = 0.0;
  if (vgd_rot_dt > 0.0) {
    // (without the cross-fade the block of the step's midpoint, as VetGdBlock)
    const Real dt = (pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    const Real tau = (tc - (xf ? 0.0 : 0.5*dt) - vgd_rot_t0)/vgd_rot_dt;
    b = static_cast<int64_t>(std::floor(tau));
    if (xf) s = tau - static_cast<Real>(b);
    // clean_rot_ring = live: a build on a block boundary closes the block (s = 1)
    if (cln_ring == 1 && s < 1.0e-9) {
      b -= 1;
      s = 1.0;
    }
  }
  // clean_rot_ring with its history: the current set and the next one, (1 - s) and s;
  // jump: the next one (weight 0) only at the build that closes the block
  const bool two = cln_ring && cln_rok;
  cln_rlast = (cln_ring == 2) &&
              tc >= vgd_rot_t0 + (static_cast<Real>(b) + 1.0 - 1.0e-9)*vgd_rot_dt;
  const int nx = (xf || cln_rlast) ? 1 : 0;
  const int ns = (two ? 1 : nk) + nx;
  // (the sets' tables by their own count: CleanScUnion2 resets cln_fvn to two sets)
  if (cln_fvn != ns*n || static_cast<int>(cln_fvkey.size()) != ns) {
    cln_fvn = ns*n;
    cln_fvd = DvceArray2D<Real>("cln_fvd", ns*n, 4);
    cln_fvkey.assign(ns, std::numeric_limits<Real>::quiet_NaN());
    cln_fvblk.assign(ns, 0);
    if (cln_ring) {
      cln_fvdm.resize(ns);
      for (int q = 0; q < ns; ++q) cln_fvdm[q] = DvceArray2D<Real>("cln_fvdm", ns*n, 4);
    }
  }
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  for (int64_t blk = (two ? b : b - nk + 1); blk <= b + nx; ++blk) {
    const int q = static_cast<int>(((blk % ns) + ns) % ns);
    const Real key = VetGdAngleBlk(blk);
    Real w = ((xf && blk == b - nk + 1) ? 1.0 - s : ((blk > b) ? s : 1.0))/nk;
    if (two) w = (blk > b) ? s : 1.0 - s;
    const std::vector<double> d = VgdRotate(key);
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
    cln_fvkey[q] = key;
    cln_fvblk[q] = blk;
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  if (cln_ring) {
    // the directions of each slot alone, with the set's own weights
    for (int q = 0; q < ns; ++q) {
      auto m_h = Kokkos::create_mirror_view(cln_fvdm[q]);
      const std::vector<double> d = VgdRotate(cln_fvkey[q]);
      for (int e = 0; e < ns*n; ++e) {
        for (int c = 0; c < 3; ++c) m_h(e,c) = d_h(e,c);
        m_h(e,3) = (e/n == q) ? d[4*(e%n)+3] : 0.0;
      }
      Kokkos::deep_copy(cln_fvdm[q], m_h);
    }
    cln_rns = ns;
    cln_rs = s;
    cln_rbld = b;
  }
}
#endif  // M1_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn int CleanFvMirror
//! \brief the leaving direction (of d0 .. d0+nd-1) nearest to the mirror image
//! n - 2 (n . t) t of an incoming direction d at a wall of unit normal t, sg t outward
//! of the cell (-1: d is not incoming); u = |n . t|.  Of directions as near as each other
//! to round-off (a symmetric image: the cubed sphere's faces) the first, so that the map
//! does not depend on the last bits of the geometry (GPU against CPU)

KOKKOS_INLINE_FUNCTION
int CleanFvMirror(const DvceArray2D<Real> &dir_, const int d0, const int nd, const int d,
                  const Real t[3], const Real sg, Real &u) {
  const Real nt = dir_(d,0)*t[0] + dir_(d,1)*t[1] + dir_(d,2)*t[2];
  u = fabs(nt);
  int r = -1;
  if (!(sg*nt < 0.0)) return r;
  const Real v[3] = {dir_(d,0) - 2.0*nt*t[0], dir_(d,1) - 2.0*nt*t[1],
                     dir_(d,2) - 2.0*nt*t[2]};
  Real bd = -2.0;
  for (int e = d0; e < d0 + nd; ++e) {
    const Real et = dir_(e,0)*t[0] + dir_(e,1)*t[1] + dir_(e,2)*t[2];
    const Real dt = v[0]*dir_(e,0) + v[1]*dir_(e,1) + v[2]*dir_(e,2);
    if (sg*et > 0.0 && dt > bd + 1.0e-12) {bd = dt; r = e;}
  }
  return r;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRayCut
//! \brief clean_ray_cut, from the rays' source of this build (vlat_cs, VetLatCsFinish):
//! the first shell swept, vlat_icut (one global min over the columns).  Below it the
//! closure is the diffusion limit (CleanThickI), whose error is O(R^2) with R = |grad
//! S|/(chi S): a column's cut lies clean_ray_cut_margin optical depths below its first
//! cell from the bottom up with R >= clean_ray_cut_r (central differences along all three
//! axes, one-sided at a physical face) or thinner than vet_col_lat_cut_thin, and at
//! tau_top >= clean_ray_cut_taumin.  The cut moves down at once and up only after
//! clean_ray_cut_wait builds that would raise it; the cut and that count go into the
//! restart file (clean_ray_cut_now, _count).  The staggered layout's sc rays only
//! (clean_rays = sc: the vet_gd sweep starts at the cut with the diffusion inflow);
//! several MeshBlocks along x1: no cut.

#if M1_TU_PART(1)
void RadiationM1::CleanRayCut() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  vlat_icut = 0;
  if (!cln_rcut_on || !cln_init || pmy_pack->pmesh->mesh_indcs.nx1 != indcs.nx1) return;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, n1 = indcs.nx1;
  const int nk = ke - ks + 1, nj = je - js + 1, nf = pmy_pack->pmesh->three_d ? 3 :
                                                  (pmy_pack->pmesh->multi_d ? 2 : 1);
  const Real tcut = cln_rcut_tau, tthin = vlat_cut_thin, er = cln_rcut_r;
  const Real mg = cln_rcut_mg;
  CleanArgs x = CleanMakeArgs();
  auto cs_ = vlat_cs;
  int lmin = 2*n1;
  Kokkos::parallel_reduce("cln_raycut", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmy_pack->nmb_thispack*nk*nj),
  KOKKOS_LAMBDA(const int idx, int &lm) {
    const int m = idx/(nk*nj), k = ks + (idx/nj) % nk, j = js + idx % nj;
    auto tc = [&](const int i) {return exp(cs_(m,0,k,j,i))*(x.x1f(m,i+1) - x.x1f(m,i));};
    int lt = 0;   // (tau_top >= clean_ray_cut_taumin, as vet_col_lat_taucut)
    Real tau = 0.0;
    for (int i = ie; i >= is; --i) {
      if (tau + 0.5*tc(i) >= tcut) {lt = i - is; break;}
      tau += tc(i);
    }
    int i = is;
    for (; i <= ie && !(tthin > 0.0 && tc(i) < tthin); ++i) {
      Real r2 = 0.0;
      for (int dir = 0; dir < nf; ++dir) {
        const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
        Real g = 0.0;
        for (int up = 0; up < 2; ++up) {
          const int s = up ? 1 : -1;
          if (CleanBoundaryFace(x, dir, m, k + up*dk, j + up*dj, i + up*di)) continue;
          g = fmax(g, fabs(cs_(m,1,k+s*dk,j+s*dj,i+s*di) - cs_(m,1,k,j,i))/
                      x.dxf[dir](m,k+up*dk,j+up*dj,i+up*di));
        }
        r2 += g*g;
      }
      if (!(sqrt(r2) < er*exp(cs_(m,0,k,j,i)))) break;
    }
    int lc = i - is;
    for (Real t = 0.0; lc > 0 && t < mg;) t += tc(is + --lc);
    // (2 lc + 1 where tau_top binds the column)
    lc = (lc <= lt) ? 2*lc : 2*lt + 1;
    lm = (lc < lm) ? lc : lm;
  }, Kokkos::Min<int>(lmin));
  lmin = -static_cast<int>(CleanReduce(-lmin, true));
  const bool bt = (lmin % 2 == 1);
  lmin = std::max(0, std::min(lmin/2, n1 - 2));
  const int prev = cln_rcut;
  cln_rcut_n = (lmin > cln_rcut && cln_rcut >= 0) ? cln_rcut_n + 1 : 0;
  if (lmin < cln_rcut || cln_rcut < 0 || cln_rcut_n >= cln_rcut_wait) {
    cln_rcut = lmin;
    cln_rcut_n = 0;
  }
  if (global_variable::my_rank == 0 && prev != cln_rcut) {
    std::cout << "CLEAN ray cut cycle " << pmy_pack->pmesh->ncycle << " shell "
              << prev << " -> " << cln_rcut << " (want " << lmin << ", "
              << (bt ? "tau_top" : "R") << ")" << std::endl;
  }
  pin_rcut_->SetInteger("rad_m1", "clean_ray_cut_now", cln_rcut);
  pin_rcut_->SetInteger("rad_m1", "clean_ray_cut_count", cln_rcut_n);
  vlat_icut = cln_rcut;
}
#endif  // M1_TU_PART(1)

// the directions of one cell as a strided row, the accessor's address once
// per cell instead of per direction (the same values; out of the slab: its dummy)
struct CleanRow {
  Real *p;
  int64_t s;
  KOKKOS_INLINE_FUNCTION Real operator[](const int d) const {return p[d*s];}
};
KOKKOS_INLINE_FUNCTION
CleanRow VgdRow(const RadiationM1::VgdRag &vi, const int m, const int kb, const int jb,
                const int i) {
  Real *p = &vi(m,0,kb,jb,i);
  return CleanRow{p, (p == vi.d.data() + vi.dummy) ? 0 : 1};
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFromGd
//! \brief the closure from the vet_gd intensities in the local basis (r, theta, phi) of
//! each cell, mu_b = n . e_b: J = <I>, K_bc = <I mu_b mu_c>, H+-_b = <I mu_b> over mu_b
//! >< 0, K+^d_b = <I mu_b mu_d> over mu_d > 0 and C1 = <|I - S|>/J (CLN_* layout, as
//! CleanClosureFromRays; no ray face fluxes).  On the cubed sphere mu_b = n_b (Cartesian
//! D and k+ columns) and the half ranges are about the family normals (CleanCsFamily).

#if M1_TU_PART(1)
void RadiationM1::CleanClosureFromGd() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;
  CleanArgs x = CleanMakeArgs();
  // (in device memory, CleanGdFaceSums)
  cln_gdx = CleanDevArgs(x, 3);
  const CleanArgs *px = cln_gdx;
  // clean_ray_cut: below the cut the sc rays' intensities are the diffusion limit, in
  // vgd_i for the moments here and the face sums (CleanGdFaceSums)
  if (vgd_on && vlat_icut > 0) {
    const Real ic = 1.0/c_light;
    // (one team per cell, its lanes along the directions: the stores of a warp
    // contiguous; the same values)
    const int nc = vlat_icut, nj = je - js + 1, nk = ke - ks + 1;
    const int tsz = std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace> ? 1 : 32;
    Kokkos::parallel_for("cln_gd_thk",
                         Kokkos::TeamPolicy<>(DevExeSpace(), (nmb1 + 1)*nk*nj*nc, tsz),
    KOKKOS_LAMBDA(TeamMember_t tm) {
      const int t = tm.league_rank(), i = is + t % nc, j = js + (t/nc) % nj;
      const int k = ks + (t/(nc*nj)) % nk, m = t/(nc*nj*nk);
      const CleanArgs &x = *px;
      Real e[3][3];   // (the cell's frame, F and E once, not per direction)
      CleanThickFrame(x, x3v, m, k, j, i, e);
      const Real f[3] = {iw_(m,M1_IW_F1,k,j,i), iw_(m,M1_IW_F2,k,j,i),
                         iw_(m,M1_IW_F3,k,j,i)}, en = iw_(m,M1_IW_EN,k,j,i);
      const CleanRow r = VgdRow(vi_, m, k+ok, j+og, i);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, n), [&](const int d) {
        const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
        r.p[d*r.s] = CleanThickI(e, f, en, ic, nv);
      });
    });
  }
  CleanForK("cln_gd", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px;
    const Real st = sin(x2v(m,j)), ct = cos(x2v(m,j));
    const Real sp = sin(x3v(m,k)), cp = cos(x3v(m,k));
    const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
    // the cubed sphere: D in Cartesian components, half ranges about the families n_d
    Real nd[3][3];
    if (x.cub) CleanCsFamily(x, m, k, j, i, nd);
    Real q[CLN_N] = {0.0}, hp[3] = {0.0}, hm[3] = {0.0};   // (q in registers)
    const CleanRow r = VgdRow(vi_, m, k+ok, j+og, i);
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      Real mu[3] = {nx*st*cp + ny*st*sp + nz*ct, nx*ct*cp + ny*ct*sp - nz*st,
                    -nx*sp + ny*cp}, nu[3];
      for (int b = 0; b < 3; ++b) {
        if (!x.sph) mu[b] = dir_(d,b);
        nu[b] = x.cub ? nx*nd[b][0] + ny*nd[b][1] + nz*nd[b][2] : mu[b];
      }
      const Real nv[3] = {nx, ny, nz};
      const Real iv = r[d];
      const Real a = dir_(d,3)*iv;
      q[CLN_J] += a;
      for (int b = 0; b < 3; ++b) {
        q[CLN_D11+b] += a*mu[b]*mu[b];
        if (!(nu[b] > 0.0)) {
          hm[b] += a*nu[b];
          continue;
        }
        hp[b] += a*nu[b];
        for (int c = 0; c < 3; ++c) q[CLN_KP+3*b+c] += a*mu[c]*nu[b];
      }
      q[CLN_D11+3] += a*mu[0]*mu[1];
      q[CLN_D11+4] += a*mu[0]*mu[2];
      q[CLN_D11+5] += a*mu[1]*mu[2];
      q[CLN_C1] += dir_(d,3)*fabs(iv - s0);
    }
    for (int b = 0; b < 3; ++b) {
      q[CLN_HP+b] = hp[b];
      q[CLN_HM+b] = hm[b];
    }
    q[CLN_JR] = q[CLN_J];
    if (q[CLN_J] < ray_min*iw_(m,M1_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    CleanStoreMoments(cs, valid, q, m, k, j, i);
  });
}
#endif  // M1_TU_PART(1)

//----------------------------------------------------------------------------------------
//! \struct CleanGdI
//! \brief the sc rays' intensities in the mesh layout without a copy (CleanGdFaceSums,
//! CleanFaceHalf): vgd_i, the first ghost layer of the MeshBlock's faces from gl (side
//! 0: j = js - 1, 1: k = ks - 1, (r, t) = (k - ks, i), (j - js, i); the cubed sphere
//! (six) also 2: j = je + 1, 3: k = ke + 1 alike, 4: i = is - 1, 5: i = ie + 1 at
//! (k - ks, j - js))

struct CleanGdI {
  RadiationM1::VgdRag vi;
  DvceArray5D<Real> gl;
  int ok, og, ks, js;
  bool six = false;
  int ke = 0, je = 0, is = 0, ie = 0;
  // (the ghost side of (k, j, i) or -1, and its two indices)
  KOKKOS_INLINE_FUNCTION
  int Side(const int k, const int j, const int i, int &r, int &t) const {
    if (j < js || (six && j > je)) {r = k - ks; t = i; return (j < js) ? 0 : 2;}
    if (k < ks || (six && k > ke)) {r = j - js; t = i; return (k < ks) ? 1 : 3;}
    if (six && (i < is || i > ie)) {r = k - ks; t = j - js; return (i < is) ? 4 : 5;}
    return -1;
  }
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int d, const int k, const int j, const int i) const {
    int r, t;
    const int s = Side(k, j, i, r, t);
    return (s >= 0) ? gl(m,d,s,r,t) : vi(m,d,k+ok,j+og,i);
  }
  KOKKOS_INLINE_FUNCTION
  CleanRow At(const int m, const int k, const int j, const int i) const {
    int r, t;
    const int s = Side(k, j, i, r, t);
    if (s >= 0) return CleanRow{&gl(m,0,s,r,t), static_cast<int64_t>(gl.stride(1))};
    return VgdRow(vi, m, k+ok, j+og, i);
  }
};

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFaceHalf
//! \brief the sc rays on the cubed sphere, staggered layout: the half ranges of every
//! face about its own normal (CleanFaceCoefficients), as the sums of the union's beams of
//! the face's two cells, J_L, J_R, H+_L, H-_R (cln_fh, 4 per face family), folded set by
//! set right after its sweep (the set's intensities in vgd_i with the ghost layers of
//! CleanGdFaceSums, or a version-5 restart's in cln_fvi from row fo; its union rows o of
//! dv): no copy of the intensities.  zero: the first set of a build

#if M1_TU_PART(2)
void RadiationM1::CleanFaceHalf(const DvceArray2D<Real> &dv, const int o,
                                const bool zero, const int fo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, n = vgd_n, nmb1 = pmy_pack->nmb_thispack - 1;
  if (cln_fh.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(cln_fh, nmb1 + 1, 12, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
  }
  auto fh = cln_fh;
  auto fi = cln_fvi;
  auto gn = cln_gn;
  // (the swept set in vgd_i and its ghost layers, CleanGdFaceSums; a version-5 restart's
  // intensities in cln_fvi)
  const bool fv = cln_fvi.extent_int(0) > 0;
  const CleanGdI gi{vgd_i, cln_gdgl, vgd_w - ks, vgd_w - indcs.ng, ks, js, true, ke, je,
                    is, ie};
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    CleanForK("cln_facehalf", nmb1, ks, ke+dk, js, je+dj, is, ie+di,   // (HIP: k fastest)
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      // (J_L, J_R, H+_L, H-_R as scalars: nvcc -O3 miscompiles a runtime-indexed local)
      Real jl0 = zero ? 0.0 : fh(m,4*dir,k,j,i), jr0 = zero ? 0.0 : fh(m,4*dir+1,k,j,i);
      Real hp0 = zero ? 0.0 : fh(m,4*dir+2,k,j,i), hm0 = zero ? 0.0 : fh(m,4*dir+3,k,j,i);
      // (the two cells' intensity rows once, not an index per direction)
      const CleanRow gl_ = fv ? CleanRow{nullptr, 0} : gi.At(m, kl, jl, il);
      const CleanRow gr_ = fv ? CleanRow{nullptr, 0} : gi.At(m, k, j, i);
#pragma unroll 4   // (rays in flight: hipcc kept one at a time since 089efb501)
      for (int e = 0; e < n; ++e) {
        const int d = o + e;
        Real mu = 0.0;
        for (int b = 0; b < 3; ++b) mu += dv(d,b)*gn.At<CLG_N>(m,3*dir+b,k,j,i);
        const Real il_ = fv ? fi(m,fo+e,kl,jl,il) : gl_[e];
        const Real ir_ = fv ? fi(m,fo+e,k,j,i) : gr_[e];
        const Real wl = dv(d,3)*il_, wr = dv(d,3)*ir_;
        jl0 += wl;
        jr0 += wr;
        if (mu > 0.0) hp0 += wl*mu;
        if (!(mu > 0.0)) hm0 += wr*mu;
      }
      fh(m,4*dir,k,j,i) = jl0;
      fh(m,4*dir+1,k,j,i) = jr0;
      fh(m,4*dir+2,k,j,i) = hp0;
      fh(m,4*dir+3,k,j,i) = hm0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGdFaceSums
//! \brief the vet_gd rays' face fluxes for clean_consistency (CLN_FX, CLN_F2, CLN_F3,
//! CLN_F1: sum_d w (n_d . A_f) I_f / A_f and the stress columns in the face's local
//! basis), A_f the exact vector area of the face (closed for every cell).  A theta or
//! phi face carries the intensity of its upwind cell (none on a physical boundary), and
//! so does an r face of the staggered layout.  Else the
//! r faces of a column take the cell balances from the inner face: the outer r face of a
//! cell gets what its inner one carries plus the exchange V chi (S - J, -H) less the
//! lateral outflow (the momentum in the column's basis), so that the energy balance of
//! every cell is the vet_gd exchange to round-off.  CLN_G is what the momentum balance of
//! the rays (in the cell's basis) differs by from the cell layout's divergence of these
//! face components: the curvature of their stress.  Several MeshBlocks along x1: a level
//! starts from the top face of the level below (vgd_x1b -> vgd_x1r), one level per pass;
//! at a reflecting top the remainder is spread over the outermost level (zero wall flux).
//! On the cubed sphere the momentum is Cartesian and the r faces carry the rays' upwind
//! fluxes (no curvature terms).

void RadiationM1::CleanGdFaceSums() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
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
  const bool f1 = cln_s.extent_int(1) > CLN_F1;   // (the cell layout's CLN_F1)
  CleanArgs x = CleanMakeArgs();
  // (the kernels read the args in device memory: a light-weight launch without a host
  // wait)
  cln_gdx = CleanDevArgs(x, 3);
  const CleanArgs *px = cln_gdx;
  // nx3 = 1 (axisymmetric): no phi face is a boundary; the upwind cell of the low one is
  // the own cell with the directions rotated back by the phi extent (the wall map), the
  // high one carries the low one's sums in its own basis
  const bool axi = !pmy_pack->pmesh->three_d && x.sph;
  auto bnd = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i) {
    const CleanArgs &x = *px;
    return !(axi && dir == 2) && CleanBoundaryFace(x, dir, m, k, j, i);
  };
  // the inner face: no flux at a reflecting one, else the outward beams of the sweep's
  // emission (vet_bc_inner emit) or of the first cell, the inward ones of the first cell
  const bool open = (ibc_x1min != M1_IBC_REFLECT);
  const bool emi = VetBcEmit();
  const Real iemit = VetBathLo();
  const Real ihi = (ibc_x1max == M1_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool rwsc[2] = {!open, ibc_x1max == M1_IBC_REFLECT};
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
  // the short characteristics' intensities in the mesh layout with the ghost exchange
  // (the image cell's beyond a periodic seam, re-indexed by the wall map as the sweep)
  // (clean_rot_ring: two slots, of the sets blk mod 2, on the cubed sphere for the face
  // half ranges (CleanScUnion2), on sp with a reflecting top for its mirror)
  const int ns2 = (cln_ring && cln_scmem) ? 2 : 1;
  const int doff = (ns2 == 1) ? 0 : cln_rqm*n;
  // Staggered spherical polar (one radial MeshBlock, no reflecting top): vgd_i itself,
  // and only the low theta and phi ghost layers exchanged, in chunks of clean_halo_vars
  // directions (one exchange each)
  // The staggered cubed sphere (one radial MeshBlock) likewise, with the first ghost
  // layer of every MeshBlock face (CleanGdI six: the high lateral faces' sums CLN_FH and
  // the face half ranges, CleanFaceHalf, read them): no copy of the intensities
  const bool six = x.cub && vgd_nb1 == 1 && !cln_scmem;
  const bool gl = (x.sph && !axi && pmy_pack->pmesh->three_d &&
                   vgd_nb1 == 1 && !cln_scmem) || six;
  CleanGdI gi{vi_, DvceArray5D<Real>(), ok, og, ks, js, six, ke, je, is, ie};
  if (gl) {
    const int ng = indcs.ng;
    // (the cubed sphere: chunks of a quarter of clean_halo_vars, the chunk array made for
    // this call only, the largest scratch of a build there)
    const int nc = std::min(n, six ? std::max(cln_hmax/4, 1) : cln_hmax);
    const int nt = std::max(ke - ks + 1, je - js + 1);
    // (persistent scratch: no cudaMalloc / cudaFree per build; gl's unwritten rows stay
    // the zeros of its allocation, ch is written before it is read)
    cln_gdgl = ScratchView<DvceArray5D<Real>>("cln_gd_gl", nmb1 + 1, n, six ? 6 : 2, nt,
                                              std::max(cs.extent_int(4), je - js + 1));
    gi.gl = cln_gdgl;
    // (six: over the stencil rows cln_st where they fit, dead across a build as cln_csy)
    const bool ost = cln_st.span() >= static_cast<std::size_t>(nmb1 + 1)*nc*
                     cs.extent_int(2)*cs.extent_int(3)*cs.extent_int(4);
    auto ch = (six && ost) ? DvceArray5D<Real>(cln_st.data(), nmb1 + 1, nc,
                             cs.extent_int(2), cs.extent_int(3), cs.extent_int(4))
            : six ? DvceArray5D<Real>("cln_gd_ch", nmb1 + 1, nc, cs.extent_int(2),
                                      cs.extent_int(3), cs.extent_int(4))
                  : ScratchView<DvceArray5D<Real>>("cln_gd_ch", nmb1 + 1, nc,
                                                   cs.extent_int(2), cs.extent_int(3),
                                                   cs.extent_int(4));
    auto gg = gi.gl;
    for (int d0 = 0; d0 < n; d0 += nc) {
      const int nd = std::min(nc, n - d0);
      // (only the ng layers along the MeshBlock's theta and phi faces (six: and r faces):
      // the exchange reads no other cell)
      CleanFor("cln_gd_chc", nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (j >= js + ng && j <= je - ng && k >= ks + ng && k <= ke - ng &&
            !(six && (i < is + ng || i > ie - ng))) {
          return;
        }
        const CleanRow r = VgdRow(vi_, m, k+ok, j+og, i);
        for (int c = 0; c < nd; ++c) ch(m,c,k,j,i) = r[d0+c];
      });
      CleanGhostFill(ch, 0, nd, true, 0, 1);   // (glc, glr read the first layer)
      par_for("cln_gd_glc", DevExeSpace(), 0, nmb1, 0, nd - 1, 0, nt - 1, is, ie,
      KOKKOS_LAMBDA(const int m, const int c, const int t, const int i) {
        if (t <= ke - ks) gg(m,d0+c,0,t,i) = ch(m,c,ks+t,js-1,i);
        if (t <= je - js) gg(m,d0+c,1,t,i) = ch(m,c,ks-1,js+t,i);
        if (six && t <= ke - ks) gg(m,d0+c,2,t,i) = ch(m,c,ks+t,je+1,i);
        if (six && t <= je - js) gg(m,d0+c,3,t,i) = ch(m,c,ke+1,js+t,i);
      });
      if (!six) continue;
      par_for("cln_gd_glr", DevExeSpace(), 0, nmb1, 0, nd - 1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
        gg(m,d0+c,4,k-ks,j-js) = ch(m,c,k,j,is-1);
        gg(m,d0+c,5,k-ks,j-js) = ch(m,c,k,j,ie+1);
      });
    }
  } else {
    if (cln_fvi.extent_int(1) != ns2*n) {
      Kokkos::realloc(cln_fvi, nmb1 + 1, ns2*n, cs.extent_int(2),
                      cs.extent_int(3), cs.extent_int(4));
    }
    fi_ = cln_fvi;
    CleanFor("cln_gd_copy", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < n; ++d) fi_(m,doff+d,k,j,i) = vi_(m,d,k+ok,j+og,i);
    });
    CleanGhostFill(cln_fvi, doff, n);
  }
  auto inten = [=] KOKKOS_FUNCTION (const int m, const int d, const int k, const int j,
                                    const int i) {
    return gl ? gi(m,d,k,j,i) : fi_(m,doff+d,k,j,i);
  };
  // vector area of the low face of cell (k,j,i) normal to dir and the local basis there
  // (Cartesian on the cubed sphere)
  auto face = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                   const int i, Real a[3], Real e[3][3]) {
    const CleanArgs &x = *px;
    CleanFaceVector(x, x3f, dir, m, k, j, i, a);
    if (x.sph) {
      CleanLocalBasis((dir == 1) ? x.x2f(m,j) : x.x2v(m,j),
                      (dir == 2) ? x3f(m,k) : x3v(m,k), e);
      return;
    }
    for (int b = 0; b < 3; ++b) {
      for (int c = 0; c < 3; ++c) e[b][c] = (b == c) ? 1.0 : 0.0;
    }
  };
  // beam of direction d through the low theta or phi face of cell (k,j,i), of vector
  // area a (not a boundary)
  // (gl: ro, ru the rows of the cell and of its neighbour below along dir)
  auto lat = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d, const Real a[3],
                                  const CleanRow &ro, const CleanRow &ru) {
    const CleanArgs &x = *px;
    const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
    const int kl = k - (dir == 2), jl = j - (dir == 1);
    if (!(f > 0.0)) return f*(gl ? ro[d] : inten(m, d, k, j, i));
    if (axi && dir == 2) return f*inten(m, mp_(m,kl+ok,j+og,d), k, j, i);
    const bool seam = (dir == 1) ? (j == js && x.mbc(m,2) == BoundaryFlag::periodic)
                                 : (k == ks && x.mbc(m,4) == BoundaryFlag::periodic);
    if (!seam || wl_(m,kl+ok,jl+og) == 0) return f*(gl ? ru[d] : inten(m, d, kl, jl, i));
    if (!wint) {
      const int e = mp_(m,kl+ok,jl+og,d);
      return f*(gl ? ru[e] : inten(m, e, kl, jl, i));
    }
    Real v = 0.0;
    for (int t = 0; t < 3; ++t) {
      v += w3_(m,kl+ok,jl+og,d,t)*inten(m, m3_(m,kl+ok,jl+og,d,t), kl, jl, i);
    }
    return f*v;
  };
  // sums of the beams of every direction through the face (local basis e, area ar), per
  // unit area (the face's geometry once per face, not per direction: the same numbers)
  auto add = [=] KOKKOS_FUNCTION (const Real e[3][3], const Real ar, const int d,
                                  const Real f, Real q[4]) {
    const Real w = dir_(d,3)*f/ar;
    q[0] += w;
    for (int b = 0; b < 3; ++b) {
      q[1+b] += w*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
    }
  };
  // the same sums with the momentum in Cartesian components, projected on
  // the face's basis and divided by its area ar once per face (CleanGdFaceSums' lateral
  // and r faces: 8 flops per direction instead of 26 and a division)
  auto addc = [=] KOKKOS_FUNCTION (const int d, const Real f, Real s[4]) {
    const Real w = dir_(d,3)*f;
    s[0] += w;
    for (int c = 0; c < 3; ++c) s[1+c] += w*dir_(d,c);
  };
  auto proj = [=] KOKKOS_FUNCTION (const Real e[3][3], const Real ar, const Real s[4],
                                   Real q[4]) {
    q[0] = s[0]/ar;
    for (int b = 0; b < 3; ++b) q[1+b] = (e[b][0]*s[1] + e[b][1]*s[2] + e[b][2]*s[3])/ar;
  };
  // a mirror at the r face i, per unit area: no energy, the normal momentum of the
  // beams of cell ic that hit it (f of sign sg) doubled, no tangential part
  auto mirror = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                     const int ic, const Real sg, Real q[4]) {
    const CleanArgs &x = *px;
    Real a[3], e[3][3], w[4] = {0.0}, nf[3], wn = 0.0;
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      const Real tr = inten(m, d, k, j, ic);
      if (sg*f > 0.0) add(e, x.area[0](m,k,j,i), d, 2.0*f*tr, w);
    }
    for (int b = 0; b < 3; ++b) nf[b] = CleanNf(x, b, 0, m, k, j, i);
    for (int b = 0; b < 3; ++b) wn += w[1+b]*nf[b];
    q[0] = 0.0;
    for (int b = 0; b < 3; ++b) q[1+b] = wn*nf[b];
  };
  // sp, the raw rows (gl), a whole build: both lateral faces and the r face of
  // a cell in one pass over the directions (cln_gd_fs below; clean_gd_face_fuse = false:
  // the three kernels), each sum its own terms in the same order (bitwise)
  const bool fs = cln_gdfs && gl && x.sph && !x.cub && !axi;
  for (int dir = 1; dir < 3 && !fs; ++dir) {
    const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
    if (!x.sph && dir >= x.nf) continue;
    CleanForK("cln_gd_lat", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      // the cubed sphere: the high face of the MeshBlock too (CLN_FH of its last cell)
      const bool hi = x.cub && ((dir == 1) ? j == je : k == ke);
      for (int up = 0; up <= (hi ? 1 : 0); ++up) {
        const int kf = k + up*(dir == 2), jf = j + up*(dir == 1);
        Real q[4], sc[4] = {0.0}, a[3], e[3][3];
        face(dir, m, kf, jf, i, a, e);
        const bool bd = bnd(dir, m, kf, jf, i);
        const Real ar = x.area[dir](m,kf,jf,i);
        const CleanRow ro = gl ? gi.At(m, kf, jf, i) : CleanRow{nullptr, 0};
        const CleanRow ru = (gl && !bd) ? gi.At(m, kf - (dir == 2), jf - (dir == 1), i)
                                        : CleanRow{nullptr, 0};
        for (int d = 0; d < n; ++d) {
          addc(d, bd ? 0.0 : lat(dir, m, kf, jf, i, d, a, ro, ru), sc);
        }
        proj(e, ar, sc, q);
        for (int c = 0; c < 4; ++c) cs(m,(up ? CLN_FH + 4*(dir - 1) : nf)+c,k,j,i) = q[c];
      }
    });
  }
  if (!fs) CleanGhostFill(cln_s, CLN_F2, 4);
  if (!fs) CleanGhostFill(cln_s, CLN_F3, 4);
  // the staggered layout's sc rays on spherical polar (rr): the r face i carries their
  // own upwind fluxes (area times the per-area sums), not the balance march, which loses
  // the diffusion flux in the difference V chi (S - J) of thick cells and, with the
  // bracket, feeds back on E
  const bool rr = x.sph;
  auto rflux = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                    Real q[4]) {
    const CleanArgs &x = *px;
    Real a[3], e[3][3], sc[4] = {0.0};
    face(0, m, k, j, i, a, e);
    // (the sc rays at a reflecting x1 face: the mirror images of the cell's leaving
    // beams, CleanFvMirror, scaled so that the face returns the energy it receives)
    const int ic = (i == is) ? is : ie;
    const Real sg = (i == is) ? -1.0 : 1.0, an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const Real t[3] = {a[0]/an, a[1]/an, a[2]/an};
    const bool msc = CleanBoundaryFace(x, 0, m, k, j, i) && rwsc[i != is];
    Real sm = 1.0, ho = 0.0, hm = 0.0, u;
    for (int d = 0; d < n && msc; ++d) {
      const int r = CleanFvMirror(dir_, 0, n, d, t, sg, u);
      if (r < 0) ho += dir_(d,3)*u*inten(m, d, k, j, ic);
      if (r >= 0) hm += dir_(d,3)*u*inten(m, r, k, j, ic);
    }
    if (ho > 0.0 && hm > 0.0) sm = ho/hm;
    const CleanRow rl = gl ? gi.At(m, k, j, i - 1) : CleanRow{nullptr, 0};
    const CleanRow rh = gl ? gi.At(m, k, j, i) : CleanRow{nullptr, 0};
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      Real iv = gl ? ((f > 0.0) ? rl[d] : rh[d])
                   : inten(m, d, k, j, (f > 0.0) ? i - 1 : i);
      if (CleanBoundaryFace(x, 0, m, k, j, i)) {   // the sweep's x1 boundary
        const Real s0 = exp(cw_(m,1,k+ok,j+og,is));
        iv = (i == is) ? ((f > 0.0) ? (emi ? iemit : s0) : iv)
                       : ((f > 0.0) ? iv : ihi);
        for (int q = 0; q < bmn && i == is && f > 0.0; ++q) {
          if (d == bmr[q]) iv += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
        }
        const int r = msc ? CleanFvMirror(dir_, 0, n, d, t, sg, u) : -1;
        if (r >= 0) iv = sm*inten(m, r, k, j, ic);
      }
      addc(d, f*iv, sc);
    }
    proj(e, 1.0, sc, q);   // (the face's sums, not per area)
  };
  auto xf_ = vgd_x1f;
  auto qb_ = vgd_x1b;
  auto qr_ = vgd_x1r;
  const int nj = je - js + 1, nkj = (ke - ks + 1)*nj;
  // rr: the rays' r face fluxes, all faces at once (the columns only read them)
  DvceArray5D<Real> rf = ScratchView<DvceArray5D<Real>>("cln_gd_rf",
                       rr ? nmb1 + 1 : 1, 4, rr ? ke - ks + 1 : 1, rr ? nj : 1,
                       rr ? ie - is + 2 : 1);
  Kokkos::deep_copy(DevExeSpace(), rf, 0.0);   // (the zeros of a fresh one)
  if (fs) {
    // (cells is..ie: the low theta and phi faces as cln_gd_lat; faces is..ie+1: the r
    // face as cln_gd_rf, a boundary one by rflux after the loop)
    CleanForK("cln_gd_fs", nmb1, ks, ke, js, je, is, ie + 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      const bool cl = (i <= ie), rb = CleanBoundaryFace(x, 0, m, k, j, i);
      Real a1[3], a2[3], a0[3], e[3][3], q[4];
      Real s1[4] = {0.0}, s2[4] = {0.0}, s0[4] = {0.0};
      const bool b1 = !cl || bnd(1, m, k, j, i), b2 = !cl || bnd(2, m, k, j, i);
      const CleanRow ro = gi.At(m, k, j, i);
      const CleanRow u1 = b1 ? CleanRow{nullptr, 0} : gi.At(m, k, j - 1, i);
      const CleanRow u2 = b2 ? CleanRow{nullptr, 0} : gi.At(m, k - 1, j, i);
      const CleanRow rl = gi.At(m, k, j, i - 1);
      if (cl) {
        CleanFaceVector(x, x3f, 1, m, k, j, i, a1);
        CleanFaceVector(x, x3f, 2, m, k, j, i, a2);
      }
      if (!rb) CleanFaceVector(x, x3f, 0, m, k, j, i, a0);
      for (int d = 0; d < n; ++d) {
        if (cl) {
          addc(d, b1 ? 0.0 : lat(1, m, k, j, i, d, a1, ro, u1), s1);
          addc(d, b2 ? 0.0 : lat(2, m, k, j, i, d, a2, ro, u2), s2);
        }
        if (!rb) {
          const Real f = dir_(d,0)*a0[0] + dir_(d,1)*a0[1] + dir_(d,2)*a0[2];
          addc(d, f*((f > 0.0) ? rl[d] : ro[d]), s0);
        }
      }
      if (cl) {
        face(1, m, k, j, i, a1, e);
        proj(e, x.area[1](m,k,j,i), s1, q);
        for (int c = 0; c < 4; ++c) cs(m,CLN_F2+c,k,j,i) = q[c];
        face(2, m, k, j, i, a2, e);
        proj(e, x.area[2](m,k,j,i), s2, q);
        for (int c = 0; c < 4; ++c) cs(m,CLN_F3+c,k,j,i) = q[c];
      }
      if (rb) {
        rflux(m, k, j, i, q);
      } else {
        face(0, m, k, j, i, a0, e);
        proj(e, 1.0, s0, q);
      }
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
    CleanGhostFill(cln_s, CLN_F2, 4);
    CleanGhostFill(cln_s, CLN_F3, 4);
  } else if (rr) {
    CleanForK("cln_gd_rf", nmb1, ks, ke, js, je, is, ie + 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4];
      rflux(m, k, j, i, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
  }
  // the lateral faces' part of cell (k,j,i): their sums into q, the curvature into g
  auto latg = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                   const Real e[3][3], Real q[4], Real g[2]) {
    const CleanArgs &x = *px;
    Real a[3], ef[3][3];
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
          if (b < 2) g[b] += ((dir == 2 && axi) ? 0.0 : sg*cs(m,nf+1+b,kc,jj,i));
          for (int c = 0; c < 3; ++c) {
            const Real r = sg*cs(m,nf+1+c,kc,jj,i)*(ef[c][0]*e[b][0] +
                           ef[c][1]*e[b][1] + ef[c][2]*e[b][2]);
            q[1+b] += r;
            if (b < 2) g[b] -= r;
          }
        }
      }
    }
  };
  // rr on sp without a reflecting top: no march (every r face is the rays'), one thread
  // per cell (the same numbers); the columns then only close their ends
  const bool cel = rr && !x.cub && ibc_x1max != M1_IBC_REFLECT;
  if (cel) {
    CleanFor("cln_gd_cel", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      Real a[3], e[3][3], q[4], g[2] = {0.0, 0.0};
      face(0, m, k, j, is, a, e);
      for (int c = 0; c < 4; ++c) q[c] = rf(m,c,k-ks,j-js,i-is);
      for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
      const Real qi = q[2];
      latg(m, k, j, i, e, q, g);
      const Real q2 = rf(m,2,k-ks,j-js,i+1-is);
      g[1] += q2*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) - qi*(1.0 - x.x1f(m,i)/x.x1v(m,i));
      for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
    });
  }
  for (int lev = 0; lev < vgd_nb1; ++lev) {
    if (lev > 0) VetGdX1Pass(true);
    CleanForCol("cln_gd_col", nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const CleanArgs &x = *px;
      Real a[3], e[3][3], q[4] = {0.0};
      const bool in = (xf_(m) & 1);
      const int kj = (k - ks)*nj + j - js;
      face(0, m, k, j, is, a, e);
      for (int d = 0; open && in && !rr && d < n; ++d) {   // (rr: the rays' rf below)
        const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
        const Real iv = (emi && f > 0.0) ? iemit : inten(m, d, k, j, is);
        add(e, x.area[0](m,k,j,is), d, f*iv, q);
      }
      if (in && !open && !rr) mirror(m, k, j, is, is, -1.0, q);
      for (int c = 0; c < 4; ++c) {
        q[c] = in ? q[c]*x.area[0](m,k,j,is) : qr_(m,c*nkj + kj);
      }
      for (int c = 0; c < 4 && rr; ++c) q[c] = rf(m,c,k-ks,j-js,0);
      const Real vb = in ? 0.0 : qr_(m,4*nkj + kj);   // the column's volume below
      for (int i = is; i <= ie && !cel; ++i) {
        for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
        if (x.cub) {   // (no curvature terms: Cartesian momentum)
          rflux(m, k, j, i + 1, q);
          continue;
        }
        Real g[2] = {0.0, 0.0}, qi = q[2];
        // the exchange V chi (S - J, -H_b) and the lateral outflow of the cell
        const Real chi = x.vol(m,k,j,i)*exp(cw_(m,0,k+ok,j+og,i));
        const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
        for (int d = 0; d < n && !rr; ++d) {   // (rr: the face fluxes are the rays')
          const Real di = dir_(d,3)*inten(m, d, k, j, i);
          q[0] += chi*(dir_(d,3)*s0 - di);
          for (int b = 0; b < 3; ++b) {
            q[1+b] -= chi*di*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
          }
        }
        latg(m, k, j, i, e, q, g);
        for (int c = 0; c < 4 && rr; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
        // the angular-momentum weight r_f/r of the y_theta row on the r faces
        g[1] += q[2]*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) - qi*(1.0 - x.x1f(m,i)/x.x1v(m,i));
        for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
      }
      for (int c = 0; c < 4 && cel; ++c) q[c] = rf(m,c,k-ks,j-js,ie+1-is);
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
        if (f1) cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
        if (f1) cs(m,CLN_F1+c,k,j,ie) = q[c]/x.area[0](m,k,j,ie+1);
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
      const CleanArgs &x = *px;
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
          if (i == is && f1) cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
          if (i == ie + 1 && f1) cs(m,CLN_F1+c,k,j,ie) -= dq/x.area[0](m,k,j,ie+1);
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
  CleanFor("cln_rays", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = vc(m,M1_VET_J,k,j,i);
    for (int n = 0; n < M1_VET_NFS; ++n) {   // (F2, F1, FX, F3)
      const int c = (n < 4) ? CLN_F2 + n : ((n < 8) ? CLN_F1 + n - 4 : CLN_FX + n - 8);
      if (c < cs.extent_int(1)) cs(m,c,k,j,i) = f2 ? vc(m,M1_VET_F2+n,k,j,i) : 0.0;
    }
    cs(m,CLN_JR,k,j,i) = jr;
    if (jr > 0.0 && jr >= ray_min*iw_(m,M1_IW_EN,k,j,i)) {
      for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = vc(m,M1_VET_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,M1_VET_HP1+a,k,j,i)/jr;
        const Real h = vc(m,M1_VET_H1+a,k,j,i)/jr;
        cs(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cs(m,CLN_HM+a,k,j,i) = fmin(h - hp, 0.0);
      }
      for (int a = 0; a < 9 && cs.extent_int(1) > CLN_KP; ++a) {
        cs(m,CLN_KP+a,k,j,i) = vc(m,M1_VET_KP+a,k,j,i)/jr;
      }
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
      for (int a = 0; a < 9 && cs.extent_int(1) > CLN_KP; ++a) {
        cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
      }
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
  // (D, h+- and, where held, k+: the slots t < CLN_C1, then CLN_KP + t - CLN_C1)
  const int nfl = (cln_s.extent_int(1) > CLN_KP) ? CLN_C1 + 9 : CLN_C1;

  for (int pass = 0; pass < cln_fillpass; ++pass) {
    // (the last chunk of cln_s and the marks in one exchange: the same ghosts)
    // (only what cln_fill reads, D and h+-, valid(0), depth 1: every caller exchanges all
    // of cln_s after the fill; the last pass also valid(1), as the ghosts were left)
    const int ns = cln_s.extent_int(1), nf = (ns/cln_hmax)*cln_hmax;
    if (nfl == CLN_C1 && cln_hmax > CLN_C1) {
      CleanGhostDefer(cln_s, 0, CLN_C1);
      cln_hpend.back().dep = 1;
      CleanGhostFill(cln_v, 0, (pass + 1 < cln_fillpass) ? 1 : 2, true, 0, 1);
    } else {
      if (nf > 0) CleanGhostFill(cln_s, 0, nf);
      if (ns > nf) CleanGhostDefer(cln_s, nf, ns - nf);
      CleanGhostFill(cln_v, 0, 2);
    }
    // valid(1) marks the cells filled in this pass; only cells valid before it are read
    par_for("cln_fill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,1,k,j,i) = valid(m,0,k,j,i);
      if (valid(m,0,k,j,i) > 0.0) return;
      Real sum[CLN_C1 + 9] = {0.0};
      int nv = 0;
      for (int kk = k-dk; kk <= k+dk; ++kk) {
        for (int jj = j-dj; jj <= j+dj; ++jj) {
          for (int ii = i-1; ii <= i+1; ++ii) {
            if (valid(m,0,kk,jj,ii) <= 0.0) continue;
            ++nv;
            for (int t = 0; t < nfl; ++t) {
              sum[t] += cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,kk,jj,ii);
            }
          }
        }
      }
      if (nv == 0) return;
      valid(m,1,k,j,i) = 1.0;
      for (int t = 0; t < nfl; ++t) {
        cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,k,j,i) = sum[t]/nv;
      }
    });
    par_for("cln_fill_mark", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,0,k,j,i) = valid(m,1,k,j,i);
    });
  }
}
#endif  // M1_TU_PART(2)

}  // namespace radm1
