//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_clean_closure.cpp
//! \brief Closure of the clean implicit radiation solve: the formal solution on the
//! predicted end-of-step state gives each cell's Eddington tensor D, half-range factors
//! h+- and non-locality C1.

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
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "utils/imex_sdirk2.hpp"
#include "utils/scratch_view.hpp"
// VET_TU_PART: compiled as 3 TUs (this file part 0, rad_vet_clean_closure_p<k>.cpp part
// k) to cut nvcc/ptxas time; each member function is in exactly one part
#ifndef VET_TU_PART_N
#define VET_TU_PART_N 0
#endif
#define VET_TU_PART(k) (VET_TU_PART_N == (k))

namespace radvet {
extern bool g_vra_adv;   // vet_rotate_advance (rad_vet_vetgd.cpp)
//----------------------------------------------------------------------------------------
//! \fn void CleanStoreMoments
//! \brief a cell's closure from its ray moments q (CLN_*, not normalised): D = K/J,
//! h+- = H+-/J, k+ = K+/J, C1 = <|I - S|>/J; J = 0: isotropic, marked for the fill

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
    for (int a = 0; a < 9 && CLN_KP + a < ns; ++a) cs(m,CLN_KP+a,k,j,i) = q[CLN_KP+a]*ij;
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
    for (int a = 0; a < 9 && CLN_KP + a < ns; ++a) {
      cs(m,CLN_KP+a,k,j,i) = (ns == CLN_NHO) ? ((a < 3) ? 0.25 : -0.25) :
                             ((a%4 == 0) ? 1.0/6.0 : 0.0);
    }
    cs(m,CLN_C1,k,j,i) = 1.0;
    cs(m,CLN_J,k,j,i) = 0.0;
    valid(m,0,k,j,i) = 0.0;
  }
  // (a fixed trip count: q's indices known at compile time, q in registers)
  for (int n = CLN_F2; n < CLN_KP; ++n) {
    if (n < ns) cs(m,n,k,j,i) = q[n];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanFaceVector
//! \brief exact vector area a of the low face of cell (k,j,i) normal to dir (spherical
//! polar exact, cubed sphere CleanCsGeometry's, Cartesian along e_dir)

KOKKOS_INLINE_FUNCTION
void CleanFaceVector(const CleanArgs &x, const DvceArray2D<Real> &x3f, const int dir,
                     const int m, const int k, const int j, const int i, Real a[3]) {
  if (x.cub) {   // the exact closed face vectors of CleanCsGeometry
    for (int c = 0; c < 3; ++c)
      a[c] = x.area[dir](m,k,j,i)*x.g.At<CLG_N>(m,3*dir+c,k,j,i);
    return;
  }
  if (!x.sph) {
    for (int c = 0; c < 3; ++c) a[c] = (c == dir) ? x.area[dir](m,k,j,i) : 0.0;
    return;
  }
  // (a theta or phi face reads only its own angle: no j+1 (k+1) edge in 1-D (2-D))
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
//! \brief clean_ray_cut: intensity along n in a cell below the cut, the diffusion limit
//! I = J + 3 n . H of its E en and F f in the frame e (CleanThickFrame)

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
//! \fn void RadiationVet::CleanClosure
//! \brief builds the closure for t^{n+1}, once per step and dt (all stage solves end
//! there); clean_closure = eddington: isotropic

#if VET_TU_PART(0)
void RadiationVet::CleanClosure() {
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
  cln_cip = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureEddington
//! \brief clean_closure = eddington: D = I/3, h+- = +-1/4, k+ = I/6, C1 = 0

void RadiationVet::CleanClosureEddington() {
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
    for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
      cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    }
    cs(m,CLN_C1,k,j,i) = 0.0;
    cs(m,CLN_J,k,j,i) = 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureIterate
//! \brief clean_closure = iterate: closure from the formal solution of the current
//! iterate's E, T and F

void RadiationVet::CleanClosureIterate() {
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
    iw_(m,VET_IW_EN,k,j,i) = fmax(w(m,CLW_EK,k,j,i), 0.0);
    iw_(m,VET_IW_TP,k,j,i) = w(m,CLW_T,k,j,i);
    iw_(m,VET_IW_F1,k,j,i) = u0_(m,VET_F1,k,j,i);
    iw_(m,VET_IW_F2,k,j,i) = u0_(m,VET_F1+1,k,j,i);
    iw_(m,VET_IW_F3,k,j,i) = u0_(m,VET_F1+2,k,j,i);
  });
  CleanRays();
  CleanClosureFill();
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
  cln_havec = false;
  cln_cip = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRayMoments
//! \brief ray moments J, K, H+-, <|I - S|> and face fluxes (CLN_*, not normalised) into
//! dst: of a new formal solution of the current state, or (solve = false) of the closure
//! in place

void RadiationVet::CleanRayMoments(DvceArray5D<Real> dst, const bool solve) {
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
    for (int n = 0; n < CLN_N; ++n) {   // (0: a slot the staggered layout lacks)
      dst(m,n,k,j,i) = (n >= cs.extent_int(1)) ? 0.0 :
                       ((CleanCellMoment(n) && n != CLN_J) ? jr : 1.0)*cs(m,n,k,j,i);
    }
  });
  cln_havec = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureAfterSolve
//! \brief clean_closure = extrap(2): store the ray moments of a converged solution: the
//! step's last solve at t^{n+1}, extrap2 also stage 1 at its time.  extrap1: the stage-1
//! closure's moments at t^n + (a21 + g) dt; only backward Euler solves here

void RadiationVet::CleanClosureAfterSolve() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE2) return;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE1) {
    CleanRayMoments(cln_mc, false);
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    cln_tm = time + ((t2_tab == 0) ? dt : (tb.a21 + tb.g)*dt);
    cln_have_c = true;
    cln_cand_be = false;
    CleanRingJump();
    return;
  }
  if (t2_solve == VET_T2S_STAGE1) {
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
  cln_cand_be = (t2_solve != VET_T2S_STAGE2);
  CleanRingJump();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRingSolveSet
//! \brief clean_rot_ring = jump: the face half ranges use the swept set of the step's
//! rotation block

void RadiationVet::CleanRingSolveSet() {
  if (cln_cssc) return;   // (cubed-sphere sc rays: CleanScUnion2's two sets)
  const int64_t bs = VetGdBlock(pmy_pack->pmesh->ncycle);
  const int ns = static_cast<int>(cln_fvdm.size());
  for (int q = 0; q < ns && ns > 1; ++q) {
    if (cln_fvblk[q] == bs) Kokkos::deep_copy(cln_fvd, cln_fvdm[q]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRingJump
//! \brief clean_rot_ring = jump: at a rotation-block end the moments jump to the next
//! block's; the next step does not extrapolate through it

void RadiationVet::CleanRingJump() {
  if (!cln_rjump) return;
  Kokkos::deep_copy(cln_mc, cln_rnew);
  cln_cand_be = true;
  cln_rjump = false;
}
#endif  // VET_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn bool CleanRealizable
//! \brief true if moments q (J > 0) are realisable: |H| <= J, K - H H^T/J PSD (to 1e-12)

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
//! \fn void RadiationVet::CleanClosureExtrap
//! \brief clean_closure = extrap(2): the closure at the stage time from ray moments of
//! converged solutions, extrapolated linearly in time (exact for S and inflow linear in
//! t).  The line runs through t^{n-1} and t^n, or (stage 2, extrap2) t^n and stage 1.
//! Cells whose extrapolate is not realisable keep the latest moments (counted); first
//! step, restart and backward Euler use t^n alone (or build it).  extrap1 with stage 1
//! before t^{n+1} (sdirk2, trbdf2): stage 1 uses the last step's moments, stage 2 the
//! line through the last and this step's stage-1 sets, weighted 1/(1 + z), z = c dt/dx
//! (no over-relaxation where the field is quasi-static)

#if VET_TU_PART(0)
void RadiationVet::CleanClosureExtrap() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  cln_nguard = 0;
  const bool sx = (cln_extrap == 3 && t2_tab != 0);   // (stage-1 moments before t^n+1)
  // (stage 2 of sx: the line through two stage-1 sets, the last not from backward
  // Euler, a restart or a ring jump)
  const bool usec = sx && t2_solve == VET_T2S_STAGE2 && cln_have_c && !cln_cand_be &&
                    cln_have_n && cln_tn < time && cln_tm > cln_tn;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE2 && !usec) return;

  // a new step: the end of the last step (sx: its stage-1 time) becomes t^n
  if (cln_have_c && (cln_tm == time || (sx && cln_tm < time))) {
    std::swap(cln_mp, cln_mn);
    std::swap(cln_mn, cln_mc);
    cln_tp = cln_tn;
    cln_tn = cln_tm;
    cln_have_p = cln_have_n && !cln_cand_be;
    cln_have_n = true;
    cln_have_c = false;
    cln_have_1 = false;
  }
  if (!cln_have_n || !(cln_tn == time || (sx && cln_tn < time))) {
    CleanClosure();
    return;
  }

  // the line M(t) = Mb + (t - tb)/(tb - ta) (Mb - Ma) at the stage time ts
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  Real ts = time + dt;
  if (t2_solve == VET_T2S_STAGE1) {
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    ts = time + (tb.a21 + tb.g)*dt;
  }
  const bool use1 = (cln_extrap == 2) && (t2_solve == VET_T2S_STAGE2) && cln_have_1 &&
                    (cln_t1 > cln_tn);
  const bool line = stage && (use1 || usec || (cln_have_p && !sx));
  auto ma = (use1 || usec) ? cln_mn : cln_mp;
  auto mb = use1 ? cln_m1 : (usec ? cln_mc : cln_mn);
  const Real tb = use1 ? cln_t1 : (usec ? cln_tm : cln_tn);
  const Real ta = (use1 || usec) ? cln_tn : cln_tp;
  const Real wx = line ? (ts - tb)/(tb - ta) : 0.0;
  auto mn = usec ? cln_mc : cln_mn;   // (a guarded cell: the latest moments)

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
    // (stage 2 of sx: the line times 1/(1 + z), z = c dt/dx of the thinnest direction)
    Real z = 0.0;
    for (int d = 0; d < x.nf && usec; ++d) {
      z = fmax(z, x.c*dt*CleanDivF(x, d, m, k, j, i, 0.0, 1.0));
    }
    const Real wc = wx/(1.0 + z);
    Real q[CLN_N];
    for (int n = 0; n < CLN_N; ++n) {
      q[n] = mb(m,n,k,j,i) + wc*(mb(m,n,k,j,i) - ma(m,n,k,j,i));
    }
    // the extrapolated moments must be those of a non-negative intensity
    const bool ok = line && ma(m,CLN_J,k,j,i) > 0.0 && q[CLN_J] > 0.0 &&
                    CleanRealizable(x, q, m, k, j, i);
    if (!ok) {
      if (line && mn(m,CLN_J,k,j,i) > 0.0) ++ng;
      // (cubed sphere: face sums, held on both sides, stay extrapolated alike)
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
#endif  // VET_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanPredictState
//! \brief predict the t^{n+1} state for the formal solution (iw: E, T, F; opac if
//! T-dependent), O(dt^2), from the hydro's explicit stage and K1 = f(U^n) (FSAL):
//!  1. momentum: hydro + dt K1_M (cancels the explicit gravity kick where supported);
//!  2. transport: dt (K1_E + K1_EN), losses beyond the content at constant relative rate;
//!  3. exchange: implicit local relaxation of (eint, E) over dt.
//! No valid K1 (first step, backward Euler): only 3.  No feedback: gas fixed.
//! use_k1 = false: the current state itself

#if VET_TU_PART(0)
void RadiationVet::CleanPredictState(const bool use_k1) {
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
  cln_opsaved = fluid && !opac_zero && (opacity_type != VET_OPAC_CONST);
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
  CleanArgs x = CleanMakeArgs();

  CleanFor("cln_predict", nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // K1 is stored on active cells only; ghost cells take the explicit state
    const bool active = (k >= ks && k <= ke && j >= js && j <= je && i >= is && i <= ie);
    const bool use_rates = rates && active;
    const bool relax = exch && feedback && active;
    const Real e0 = u0_(m,VET_E,k,j,i);
    Real e = e0, d = 1.0, t = 0.0;
    if (fluid) {
      d = uh(m,IDN,k,j,i);
      Real mom[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
      const Real ekin = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
      const Real eint_h = fmax(CleanGasEint(x, uh, ekin, m, k, j, i), 0.0);
      t = eos.Temperature(d, eint_h);
      if (use_rates && !feedback) {
        // gas held fixed: E advances at its own rate, which includes the exchange
        const Real de = dt*k1(m,VET_T2_E,k,j,i);
        e = CleanPredictE(e0, de);
      } else if (relax) {
        // transport over dt (valid K1) less the radiation force's work, then the
        // implicit local exchange of (eint, E)
        if (use_rates) {
          for (int n = 0; n < 3; ++n) mom[n] += dt*k1(m,VET_T2_M1+n,k,j,i);
          const Real ekin_new = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
          const Real de = dt*(k1(m,VET_T2_E,k,j,i) + k1(m,VET_T2_EN,k,j,i))
                          - (ekin_new - ekin);
          e = CleanPredictE(e0, de);
        }
        // E (1 + h c kE) = e + h c kP a T^4,  eint(T) + E = eint_h + e
        const Real hce = dt*cl*opac_(m,VET_OP_E,k,j,i);
        const Real hcp = dt*cl*opac_(m,VET_OP_P,k,j,i)*arad_;
        t = CleanExchangeTemperature(eos, d, hcp/(1.0 + hce), eint_h + e*hce/(1.0 + hce),
                                     t);
        e = (e + hcp*SQR(SQR(t)))/(1.0 + hce);
      }
    } else if (use_rates) {
      const Real de = dt*k1(m,VET_T2_E,k,j,i);
      e = CleanPredictE(e0, de);
    }

    if (active) {
      iw_(m,VET_IW_EN,k,j,i) = e;
      iw_(m,VET_IW_TP,k,j,i) = t;
      iw_(m,VET_IW_F1,k,j,i) = u0_(m,VET_F1,k,j,i);
      iw_(m,VET_IW_F2,k,j,i) = u0_(m,VET_F1+1,k,j,i);
      iw_(m,VET_IW_F3,k,j,i) = u0_(m,VET_F1+2,k,j,i);
    }
    if (new_opac) {
      Real op, oe, of, os;
      if (otype == VET_OPAC_TABLE) {
        M1TableOpacities(otab_, d, t, op, oe, of, os);
        M1ScatterEos(otab_, eos, d, t, of, os);
      } else {
        M1Opacities(otype, d, t, kp, ke_, kf, ksc, rref, tref, aa, bb, op, oe, of, os);
      }
      opac_(m,VET_OP_P,k,j,i) = d*op;
      opac_(m,VET_OP_E,k,j,i) = d*oe;
      opac_(m,VET_OP_T,k,j,i) = d*(of + os);
      if (scat) opac_(m,VET_OP_S,k,j,i) = d*os;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureSphere
//! \brief exact closure of a uniformly bright sphere of radius R = clean_sphere_radius
//! in vacuum.  Outside, the cone mu > mu_c = sqrt(1 - R^2/r^2): D_rr = (1 + mu_c +
//! mu_c^2)/3, D_tt = D_pp = (1 - D_rr)/2, h+_r = (1 + mu_c)/2, h-_r = 0, h+-_t = h+-_p =
//! +-S/(pi (1 - mu_c)), S = int_{mu_c}^1 sqrt(1 - mu^2) dmu, C1 = 1; k+^r = (D_rr, 0, 0),
//! k+^t = (k, D_tt/2, 0), k+^p = (k, 0, D_pp/2), k = (R/r)(1 + mu_c)/(3 pi).  Inside:
//! isotropic.  Cubed sphere: Cartesian about r_hat and the lateral family normals n_d

void RadiationVet::CleanClosureSphere() {
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
    const bool kp = cs.extent_int(1) > CLN_KP + 8;   // (k+, not CLN_HO)
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
//! \fn void RadiationVet::CleanRays
//! \brief formal solution of the state in iw (chi from opac) and its closure: vet_sc,
//! or vet_gd on curved meshes; the state's own source (no dt-dependent relaxation)

void RadiationVet::CleanRays() {
  // (the closure's time; a step of at least one rotation block takes its start's sets)
  const Real dt = pmy_pack->pmesh->dt, tc = pmy_pack->pmesh->time +
                  ((vgd_rot_dt > 0.0 && dt >= vgd_rot_dt) ? 0.0 : dt);
  if (!vgd_on && cln_ring) {   // (the rotating vet_sc set: CleanRingMoments sweeps it)
    CleanFvSets(tc, pmy_pack->pmesh->ncycle);
    CleanRingMoments();
    return;
  }
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
    iw_(m,VET_IW_KT,k,j,i) = opac_(m,VET_OP_T,k,j,i);
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
    // (cubed-sphere sc rays: also the face half ranges from cln_fvi, CleanFaceHalf)
    if (cln_consist || cln_cssc) CleanGdFaceSums();
    if (cln_cssc && cln_rfh) CleanCellHalf(vgd_dir);
    if (cln_cssc && !cln_rfh) CleanFaceHalf(vgd_dir, 0, true, 0);
    if (cln_cssc) cln_fvi = DvceArray5D<Real>();   // (build scratch)
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanScUnion2
//! \brief clean_rot_ring, cubed-sphere sc rays: the direction union (cln_fvd) of sets
//! b, b+1 (slots blk mod 2) with cross-fade weights 1 - s, s

void RadiationVet::CleanScUnion2() {
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
#endif  // VET_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn bool CleanRingMix
//! \brief clean_rot_ring: raw moments q of the ring mix: current set (slot 0), next (qn,
//! weight s), older sets r (weights 1, the last 1 - s), J-weighted with r_r =
//! J_{b-r}/J_b.  Invalid offsets (r_r < 0) skipped; false if the current set has no rays

struct CleanRingL {
  DvceArray5D<Real> rl, cs;
  Kokkos::Array<int, CLN_N> fm;
  int nf, ol, ks, js, is;
  // (slot o's normalised moment n at interior cell (k,j,i); slot ol is read from cln_s)
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
  Real a[RadiationVet::VET_RA_NF];
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
//! rr's entry nfr + 2 (0: no history, A and T the sets themselves).  The face sums (not
//! cell moments) as (1 - s) A + s T in raw moments, linear in the faces' own sums: no
//! cell's J ratio weights its low face, so that a face and its mirror face mix alike

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
    q[n] = CleanCellMoment(n) ? ((1.0 - s)*ga + s*gt)*q[CLN_J]
                              : (1.0 - s)*ga*jb*(hist ? r : 1.0) + s*gt*jt;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CleanRingFields
//! \brief clean_rot_ring: the fields the ring carries, fm[n] = slot or -1; returns their
//! number.  Left out: fields cln_s lacks, the stress parts of the face sums and (cubed
//! sphere) CLN_G (mixed as 0).  file: the restart layout; files before version 7 keep the
//! stress parts (and CLN_G on the cubed sphere)

#if VET_TU_PART(0)
int RadiationVet::CleanRingFields(int *fm, const bool file) const {
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  const int ns = CleanNs(cub, cln_rfh);   // (the slots of cln_s)
  int nf = 0;
  for (int n = 0; n < CLN_N; ++n) {
    const bool fs = (n >= CLN_F2 && n < CLN_G && (n - CLN_F2)%4 != 0) ||   // (stress
                    (n >= CLN_FH && n < CLN_F1 && (n - CLN_FH)%4 != 0);    // parts)
    const bool un = fs || (cub && (n == CLN_G || n == CLN_G + 1));
    const bool f6 = file && t2_rstv < 7;
    fm[n] = (n < ns && !(f6 ? fs && !cub : un)) ? nf++ : -1;
  }
  return nf;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRingShift
//! \brief clean_rot_ring: block b's first build shifts the ring by the pending Delta:
//! d_1 = -Delta, d_j -= Delta, r_j' = r_{j-1}/R, r_1' = 1/R

void RadiationVet::CleanRingShift(const int64_t b) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nr = vra_k - 1;
  Kokkos::Array<int, CLN_N> fm;
  const int nf = CleanRingFields(fm.data());
  const Real tk = vgd_rot_t0 + static_cast<Real>(b)*vgd_rot_dt;
  if (std::abs(cln_rpt - tk) > 1.0e-9*vgd_rot_dt) cln_roff_sum += 1.0;
  vra_gen += 1;
  auto rr = cln_rr;
  auto pd = cln_rpd;
  if (vra_ema) {
    // (EMA update with the boundary's Delta and weight a_ch; no history or no Delta:
    // the current set alone)
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
//! \fn void RadiationVet::CleanRingPreShift
//! \brief before the dt_retry checkpoint: with clean_rot_ring = live, cross-fade and
//! vet_gd_rotate_align every retry's build lies in one block, so a step starting on a
//! block boundary shifts the ring here and the checkpoint skips offsets and Delta
//! (cln_rfix).  CleanRingBlk: the block of a build at time + dt

int64_t RadiationVet::CleanRingBlk(const Real dt) const {
  const Real tau = (pmy_pack->pmesh->time + dt - vgd_rot_t0)/vgd_rot_dt;
  const int64_t b = static_cast<int64_t>(std::floor(tau));
  return (tau - static_cast<Real>(b) < 1.0e-9) ? b - 1 : b;   // (CleanFvSets)
}

// (live with the cross-fade: every build sweeps the sets b and b+1)
bool RadiationVet::CleanRingTwo() const {
  return cln_ring == 1 && vgd_xf && !(vra_on && !g_vra_adv);
}

// (cubed-sphere sc rays: true if the face half ranges cln_fh are rebuilt before every
// solve; else they are state, kept by the dt_retry checkpoint and restart files)
bool RadiationVet::CleanFhBuilt() const {
  return cln_rfh || (cln_cssc && cln_extrap == 0 && !cln_closit && !cln_edd &&
                     !cln_sphere && (!cln_ring || CleanRingTwo()));
}

void RadiationVet::CleanRingPreShift() {
  const Real dt = pmy_pack->pmesh->dt;
  cln_rfix = CleanRingTwo() && vgd_rot_dt > 0.0 && vgd_rot_align && dt < vgd_rot_dt &&
             CleanRingBlk(dt) == CleanRingBlk(1.0e-6*dt);   // (any retry's dt)
  if (cln_rfix && cln_rok && cln_rpb == cln_rb && CleanRingBlk(dt) == cln_rb + 1) {
    CleanRingShift(cln_rb + 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::CleanRingRst*
//! \brief clean_rot_ring in restart files (version >= 4): offsets cln_rr then Delta
//! cln_rpd on interior cells in the file layout (CleanRingFields; 0 where cln_rr holds
//! none; ema: 3 extra entries), flat per MeshBlock in CleanRingRstN() ghosted slabs.
//! CleanRingRstCopy packs (pk; zeros without have) or unpacks slabs c .. of array a
//! (c < 0: channel 0 is slab -c) for MeshBlocks m0 .. m0 + nmb - 1

int RadiationVet::CleanRingRstN() const {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N];
  const int nff = CleanRingFields(fm, true);
  const int64_t nv = static_cast<int64_t>(CleanRingNr(nff) + CleanRingNd(nff))*id.nx1*
                     id.nx2*id.nx3;
  const int64_t ns = static_cast<int64_t>(id.nx1 + 2*id.ng)*
                     ((id.nx2 > 1) ? id.nx2 + 2*id.ng : 1)*
                     ((id.nx3 > 1) ? id.nx3 + 2*id.ng : 1);
  return static_cast<int>((nv + ns - 1)/ns);
}

void RadiationVet::CleanRingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb,
                                    const int m0, const bool pk, const bool have) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  int fm[CLN_N], ff[CLN_N];
  const int nf = CleanRingFields(fm), nff = CleanRingFields(ff, true);
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
  // (cln_rr's field groups: K - 1 offsets or the EMA's 1; later entries copied as is)
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
  Kokkos::parallel_for("vet_ring_rst",
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
//! \fn void RadiationVet::CleanRingMoments
//! \brief clean_rot_ring: closure of the rotating union of K = vet_rotate_avg sets with
//! cross-fade, sweeping only the current set b and the next b+1:
//!   M = (1 - s)/K M_{b-K+1} + 1/K sum_{j=1}^{K-2} M_{b-j} + 1/K M_b + s/K M_{b+1}
//! older sets from the current g = M/J plus offsets d_j = g_{b-j} - g_b,
//! r_j = J_{b-j}/J_b (CleanRingMix).  The build closing block b (s = 1) gives Delta =
//! M_{b+1} - M_b; block b+1's first build shifts the ring (CleanRingShift).  Offsets are
//! piecewise constant, so the closure stays continuous in t.  A fresh start fills the
//! ring with one full sweep.  Unrealisable mixtures fall back to (1 - s) M_b + s M_{b+1}
//! (counted); linear in the moments, so face sums stay single-valued and conservative

void RadiationVet::CleanRingMoments() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nk = CleanRingK(), nr = nk - 1, ns = cln_rns;
  const int64_t b = cln_rbld;
  const Real s = cln_rs;
  const int c3 = ke - ks + 1, c2 = je - js + 1, c1 = ie - is + 1;
  // (swept moments on interior cells: current set in cln_rl slot 0, next in slot 1; a
  // fill's older sets straight into cln_rr.  fr: fields of offsets and Delta)
  Kokkos::Array<int, CLN_N> fm, fr;
  const int nf = CleanRingFields(fm.data()), nfr = CleanRingFields(fr.data());
  // (sweep order: cubed-sphere sc rays the older sets by block, then b, b+1 by slot
  // blk mod 2, as CleanFaceHalf folds them; elsewhere slot order (keeps the restart
  // bitwise).  The last swept set stays in cln_s, ol its slot (-1: rl holds both))
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
  // (aliases cln_csy where it fits: rewritten before every read, dead across a build)
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
  // vet_rotate_avg_ema: rr = EMA offsets, aux J and D11 (last boundary), set count n;
  // pd also a_ch and the next set's J and D11
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
  // moments of each swept set (ray_min tested on the mixture); cubed-sphere sc rays: the
  // union's directions first
  if (cln_cssc) CleanScUnion2();
  auto dfull = cln_fvd;
  bool fh0 = true;
  const Real ray_min = cln_raymin;
  cln_raymin = 0.0;
  // (clean_ring_jmin: J ratios only where both sets' J >= clean_ray_min E, so a source
  // hole's J ~ 1e-200 gives no huge ratio)
  const Real jmin = cln_ringjmin ? ray_min : 0.0;
  auto iwr = iw;
  auto cs = cln_s;
  const CleanRingL rl{cln_rl, cln_s, fm, nf, ol, ks, js, is};
  auto rr = cln_rr;
  // clean_rot_fuse: sweep the live sets b, b+1 in one fused VetGdSweep (vgd_xff_now;
  // own tables and intensities each); the loop below then takes each set's moments
  const bool fuse = vgd_on && cln_rfuse && !cln_scmem && ns == 2 && !fill && VetGdXffOk();
  if (vgd_on && cln_rfuse && !fuse && ns == 2 && !fill && !cln_rfuse_said &&
      global_variable::my_rank == 0) {
    std::cout << "<rad_vet> clean_rot_fuse: not taken (needs the sc rays, no "
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
    vgd_xff_n += 1.0;   // (fused build count)
  }
  for (int qo = 0; qo < ns; ++qo) {
    const int q = ord[qo];
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    const bool older = (p >= 0 && p < nr);
    const bool keep = (older && fill) || p == nr || p == nk;
    if (older && fill) vra_gen += 1;
    if (!vgd_on) {   // (vet_sc: the set of block cln_fvblk[q] alone)
      VetRotSet(VetRotBeta(cln_fvblk[q]));
      VetShortCharSweep();
      CleanClosureFromRays();
    } else {
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
      if (cln_rfh) CleanCellHalf(vgd_dir);   // (every swept set: the ring's moments)
      if (cln_cssc && !cln_rfh && (p == nr || p == nk)) {
        CleanFaceHalf(dfull, cln_rqm*vgd_n, fh0, 0);
        fh0 = false;
      }
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
    // (f0: fill outside the steady start; clean_rot_ring_fill_rst: a restart's first
    // fill from an older-format file fills as the steady start does)
    const bool f0 = !cln_insteady && !(cln_rfrst && cln_rstold && !cln_rfdone && full);
    cln_rfdone = cln_rfdone || (cln_rstold && full);
    if (shift) {
      CleanRingShift(b);
    } else if (full && ema) {
      // EMA init: J-weighted equal mix of the current and previous set, aux := current.
      // f0: the current set alone, count 1; a restart with converged EMA history
      // (cln_rhist): count N, no warm-up
      const Real n0 = cln_rhist ? vra_ema_n : 1.0;
      cln_rhist = false;
      par_for("cln_ring_einit", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i), jp = rr(m,fr[int{CLN_J}],k,j,i);
        const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
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
        rr(m,nfr+2,k,j,i) = v ? fmax(2.0, n0) : ((f0 && jb > 0.0 && jb >= jm) ? n0 : 0.0);
      });
    } else if (full) {
      // a fill: in CleanSteadyStart the older sets' offsets from their sweeps; elsewhere
      // zero offsets (J ratio 1), the older sets entering one per block via the shifts
      // (sweeping them at a step's build would add an O(dt) error)
      par_for("cln_ring_fill", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb,CLN_J,k,j,i);
        const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
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
      const Real jm = jmin*iwr(m,VET_IW_EN,k0,j0,i0);
      for (int n = 0; n < CLN_N; ++n) {
        if (fr[n] < 0) continue;
        pd(m,fr[n],k,j,i) = (n != CLN_J) ? rl(m,qn,n,k,j,i) - rl(m,qb,n,k,j,i) :
                                ((jb > 0.0 && jb >= jm && jn >= jm) ? jn/jb : -1.0);
      }
      if (!ema) return;
      // (EMA change weight a_ch from c = max |ln J_b/J_aux|, |D11_b - D11_aux|; the next
      // set's J and D11 for the aux)
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
  // the mixture; jump: at the block-closing build also the next block's window (cln_rnew)
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
    // (slot o; fields the ring does not carry: 0, never read)
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
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) q[CLN_J] = 0.0;
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
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) {
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
//! \fn void RadiationVet::CleanFvSets
//! \brief the swept directions (cln_fvd, cln_fvn) for the closure at time tc: the union
//! of the rotating sets, the set weights folded into the quadrature weights (moments are
//! weighted means: single-valued, conservative, realisable).  vet_rotate_avg = K: sets
//! b-K+1 .. b, 1/K each; vet_gd_rotate_xfade adds b+1, weights (1 - s)/K oldest and s/K
//! next, s the fraction of block b passed.  Slot blk mod (K + xfade); a function of t

void RadiationVet::CleanFvSets(const Real tc, const int cyc) {
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
#endif  // VET_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn int CleanFvMirror
//! \brief leaving direction of d0 .. d0+nd-1 nearest the mirror image n - 2 (n . t) t of
//! incoming d at a wall (sg t outward; -1: d not incoming); u = |n . t|.  Ties within
//! 1e-12 take the first (CPU/GPU independent)

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
//! \fn void RadiationVet::CleanRayCut
//! \brief clean_ray_cut: first swept shell vlat_icut (global min over columns); below
//! it the closure is the diffusion limit (CleanThickI, error O(R^2),
//! R = |grad S|/(chi S)).
//! A column's cut: clean_ray_cut_margin optical depths below its lowest cell with R >=
//! clean_ray_cut_r or thinner than vet_col_lat_cut_thin, and at tau_top >=
//! clean_ray_cut_taumin.  Moves down at once, up after clean_ray_cut_wait builds (both
//! in the restart file).  sc rays only; no cut with several MeshBlocks along x1

#if VET_TU_PART(1)
void RadiationVet::CleanRayCut() {
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
  if (global_variable::my_rank == 0 && cln_verbose > 0 && prev != cln_rcut) {
    std::cout << "CLEAN ray cut cycle " << pmy_pack->pmesh->ncycle << " shell "
              << prev << " -> " << cln_rcut << " (want " << lmin << ", "
              << (bt ? "tau_top" : "R") << ")" << std::endl;
  }
  pin_rcut_->SetInteger("rad_vet", "clean_ray_cut_now", cln_rcut);
  pin_rcut_->SetInteger("rad_vet", "clean_ray_cut_count", cln_rcut_n);
  vlat_icut = cln_rcut;
}
#endif  // VET_TU_PART(1)

// one cell's directions as a strided row: the address once per cell (dummy: stride 0)
struct CleanRow {
  Real *p;
  int64_t s;
  KOKKOS_INLINE_FUNCTION Real operator[](const int d) const {return p[d*s];}
};
KOKKOS_INLINE_FUNCTION
CleanRow VgdRow(const RadiationVet::VgdRag &vi, const int m, const int kb, const int jb,
                const int i) {
  Real *p = &vi(m,0,kb,jb,i);
  return CleanRow{p, (p == vi.d.data() + vi.dummy) ? 0 : 1};
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureFromGd
//! \brief closure from the vet_gd intensities, mu_b = n . e_b in the cell's basis:
//! J = <I>, K_bc = <I mu_b mu_c>, H+-_b = <I mu_b> over mu_b >< 0, K+^d_b = <I mu_b mu_d>
//! over mu_d > 0, C1 = <|I - S|>/J.  Cubed sphere: Cartesian D, half ranges about the
//! family normals

#if VET_TU_PART(1)
void RadiationVet::CleanClosureFromGd() {
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
  // clean_ray_cut: below the cut vgd_i holds the diffusion limit
  if (vgd_on && vlat_icut > 0) {
    const Real ic = 1.0/c_light;
    // (one team per cell, lanes along the directions: contiguous stores)
    const int nc = vlat_icut, nj = je - js + 1, nk = ke - ks + 1;
    const int tsz =
        std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace> ? 1 : 32;
    Kokkos::parallel_for("cln_gd_thk",
                         Kokkos::TeamPolicy<>(DevExeSpace(), (nmb1 + 1)*nk*nj*nc, tsz),
    KOKKOS_LAMBDA(TeamMember_t tm) {
      const int t = tm.league_rank(), i = is + t % nc, j = js + (t/nc) % nj;
      const int k = ks + (t/(nc*nj)) % nk, m = t/(nc*nj*nk);
      const CleanArgs &x = *px;
      Real e[3][3];   // (the cell's frame, F and E once, not per direction)
      CleanThickFrame(x, x3v, m, k, j, i, e);
      const Real f[3] = {iw_(m,VET_IW_F1,k,j,i), iw_(m,VET_IW_F2,k,j,i),
                         iw_(m,VET_IW_F3,k,j,i)}, en = iw_(m,VET_IW_EN,k,j,i);
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
    // cubed sphere: Cartesian D, half ranges about the families n_d
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
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    CleanStoreMoments(cs, valid, q, m, k, j, i);
  });
}
#endif  // VET_TU_PART(1)

//----------------------------------------------------------------------------------------
//! \struct CleanGdI
//! \brief sc ray intensities in the mesh layout without a copy: vgd_i, the first ghost
//! layer from gl (side 0: j = js - 1, 1: k = ks - 1; six: also 2, 3 high j, k and 4, 5
//! low, high i)

struct CleanGdI {
  RadiationVet::VgdRag vi;
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
//! \fn void RadiationVet::CleanFaceHalf
//! \brief cubed-sphere sc rays: each face's half-range sums J_L, J_R, H+_L, H-_R about
//! its normal (cln_fh), accumulated set by set after each sweep from vgd_i (or, for a
//! version-5 restart, cln_fvi from row fo); union rows o of dv; zero: first set

#if VET_TU_PART(2)
void RadiationVet::CleanFaceHalf(const DvceArray2D<Real> &dv, const int o,
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
  // (the swept set in vgd_i plus ghost layers; a version-5 restart's in cln_fvi)
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
      // (the two cells' intensity rows once)
      const CleanRow gl_ = fv ? CleanRow{nullptr, 0} : gi.At(m, kl, jl, il);
      const CleanRow gr_ = fv ? CleanRow{nullptr, 0} : gi.At(m, k, j, i);
#pragma unroll 4   // (several rays in flight on hipcc)
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
//! \fn void RadiationVet::CleanCellHalf
//! \brief vet_rotate_avg_ema_wb on the cubed sphere: CleanFaceHalf's half ranges per
//! cell, each about the normal of the face it leaves through: CLN_HO + q = h+ out of the
//! high face of family q, CLN_HO + 3 + q = h- out of the low one (over the cell's J of
//! the set dv).  The ring carries them as cell moments, so a set change keeps a static
//! state's; a face takes its two cells' (a ghost's in its own axes, CleanCsSlot)

void RadiationVet::CleanCellHalf(const DvceArray2D<Real> &dv) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, n = vgd_n;
  auto cs = cln_s;
  auto gn = cln_gn;
  const CleanGdI gi{vgd_i, cln_gdgl, vgd_w - ks, vgd_w - indcs.ng, ks, js, true,
                    indcs.ke, indcs.je, is, indcs.ie};
  CleanForK("cln_cellhalf", pmy_pack->nmb_thispack - 1, ks, indcs.ke, js, indcs.je, is,
            indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real jc = 0.0, h[6] = {0.0};
    const CleanRow r = gi.At(m, k, j, i);
    for (int e = 0; e < n; ++e) {
      const Real w = dv(e,3)*r[e];
      jc += w;
      for (int q = 0; q < 3; ++q) {
        Real mh = 0.0, ml = 0.0;
        for (int b = 0; b < 3; ++b) {
          mh += dv(e,b)*gn.At<CLG_N>(m,3*q+b,k+(q == 2),j+(q == 1),i+(q == 0));
          ml += dv(e,b)*gn.At<CLG_N>(m,3*q+b,k,j,i);
        }
        if (mh > 0.0) h[q] += w*mh;
        if (!(ml > 0.0)) h[3+q] += w*ml;
      }
    }
    const bool ok = cs(m,CLN_J,k,j,i) > 0.0 && jc > 0.0;
    for (int a = 0; a < 6; ++a) {
      cs(m,CLN_HO+a,k,j,i) = ok ? h[a]/jc : ((a < 3) ? 0.25 : -0.25);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanGdFaceSums
//! \brief vet_gd face fluxes (CLN_FX, F2, F3, F1): sum_d w (n_d . A_f) I_f / A_f and the
//! stress columns in the face basis, A_f the exact closed vector area.  Lateral faces
//! (and staggered r faces) carry the upwind intensity.  Else r faces follow the cell
//! balances from the inner face (outer = inner + V chi (S - J, -H) - lateral outflow),
//! so each cell's energy balance is the vet_gd exchange to round-off; CLN_G is the stress
//! curvature left over.  Several MeshBlocks along x1: one level per pass from the one
//! below; a reflecting top spreads the remainder over the outermost level.  Cubed
//! sphere: Cartesian momentum, upwind r faces

void RadiationVet::CleanGdFaceSums() {
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
  // (args in device memory: light launches)
  cln_gdx = CleanDevArgs(x, 3);
  const CleanArgs *px = cln_gdx;
  // nx3 = 1 (axisymmetric): no phi boundary; the low phi face's upwind cell is the cell
  // itself through the wall map, the high face carries the low one's sums in its basis
  const bool axi = !pmy_pack->pmesh->three_d && x.sph;
  auto bnd = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i) {
    const CleanArgs &x = *px;
    return !(axi && dir == 2) && CleanBoundaryFace(x, dir, m, k, j, i);
  };
  // inner face: no flux if reflecting, else outward beams from the emission (vet_bc_inner
  // emit) or the first cell, inward beams from the first cell
  const bool open = (ibc_x1min != VET_IBC_REFLECT);
  const bool emi = VetBcEmit();
  const Real iemit = VetBathLo();
  const Real ihi = (ibc_x1max == VET_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool rwsc[2] = {!open, ibc_x1max == VET_IBC_REFLECT};
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
  // intensities in the mesh layout with ghosts (re-indexed by the wall map at a periodic
  // seam); clean_rot_ring: two slots (sets blk mod 2) with cln_scmem
  const int ns2 = (cln_ring && cln_scmem) ? 2 : 1;
  const int doff = (ns2 == 1) ? 0 : cln_rqm*n;
  // gl: one radial MeshBlock (sp: no reflecting top): vgd_i itself plus exchanged first
  // ghost layers, low theta, phi only (cubed sphere, six: all faces), in direction chunks
  const bool six = x.cub && vgd_nb1 == 1 && !cln_scmem;
  const bool gl = (x.sph && !axi && pmy_pack->pmesh->three_d &&
                   vgd_nb1 == 1 && !cln_scmem) || six;
  CleanGdI gi{vi_, DvceArray5D<Real>(), ok, og, ks, js, six, ke, je, is, ie};
  if (gl) {
    const int ng = indcs.ng;
    // (cubed sphere: chunks of clean_halo_vars/4, the build's largest scratch)
    const int nc = std::min(n, six ? std::max(cln_hmax/4, 1) : cln_hmax);
    const int nt = std::max(ke - ks + 1, je - js + 1);
    // (persistent scratch; gl's unwritten rows stay zero, ch is written before read)
    cln_gdgl = ScratchView<DvceArray5D<Real>>("cln_gd_gl", nmb1 + 1, n, six ? 6 : 2, nt,
                                              std::max(cs.extent_int(4), je - js + 1));
    gi.gl = cln_gdgl;
    // (six: aliases cln_st where it fits, dead across a build)
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
      // (only the ng layers the exchange reads)
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
  // vector area of the low face normal to dir and the local basis there
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
  // beam of direction d through the low lateral face of vector area a (not a boundary);
  // gl: ro, ru the rows of the cell and its lower neighbour
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
  // vet_gd_wall_sym / vet_gd_theta_seam_sym: a seam face's sums q averaged with those
  // seen from the other side (own beams exact, incoming via the mirrored wall map)
  const bool wsx = vgd_wsx, wst = vgd_wst;
  const Real wsc = cos(vgd_wsw), wss = sin(vgd_wsw);
  auto mi_ = vgd_mir;
  auto mz_ = vgd_mirz;
  auto seam2 = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                    const int i, const Real a[3], const Real e[3][3],
                                    const Real ar, Real q[4]) {
    const CleanArgs &x = *px;
    const int p = (dir == 2) ? 1 : 0;
    if (!(p ? wsx : wst) || (p ? (k != ks) : (j != js)) ||
        x.mbc(m, p ? 4 : 2) != BoundaryFlag::periodic ||
        wl_(m,k-p+ok,j-1+p+og) == 0) {return;}
    // the other side's face vector b and basis eb: phi, rotated about z by W; theta,
    // Z = diag(1, 1, -1): b = -Z a, eb = (Z e_r, -Z e_theta, Z e_phi)
    Real b[3], eb[3][3], s[4] = {0.0}, qb[4];
    for (int c = 0; c < 3; ++c) {
      const Real zc = (c == 2) ? -1.0 : 1.0;
      b[c] = p ? ((c == 0) ? wsc*a[0] - wss*a[1] : (c == 1) ? wss*a[0] + wsc*a[1] : a[2])
               : -zc*a[c];
      for (int t = 0; t < 3; ++t) {
        eb[t][c] = p ? ((c == 0) ? wsc*e[t][0] - wss*e[t][1]
                                 : (c == 1) ? wss*e[t][0] + wsc*e[t][1] : e[t][2])
                     : ((t == 1) ? -zc : zc)*e[t][c];
      }
    }
    const int kl = k - p, jl = j - 1 + p;
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*b[0] + dir_(d,1)*b[1] + dir_(d,2)*b[2];
      const int dm = p ? mi_(d) : mz_(d);
      Real v = (f > 0.0) ? inten(m, d, kl, jl, i) : 0.0;
      for (int t = 0; t < 3 && !(f > 0.0); ++t) {
        const int r = m3_(m,kl+ok,jl+og,dm,t);
        v += w3_(m,kl+ok,jl+og,dm,t)*inten(m, p ? mi_(r) : mz_(r), k, j, i);
      }
      const Real w = dir_(d,3)*f*v;
      s[0] += w;
      for (int t = 0; t < 3; ++t) s[1+t] += w*dir_(d,t);
    }
    for (int t = 0; t < 4; ++t) {
      qb[t] = (t == 0) ? s[0]/ar
                       : (eb[t-1][0]*s[1] + eb[t-1][1]*s[2] + eb[t-1][2]*s[3])/ar;
      q[t] = 0.5*(q[t] + qb[t]);
    }
  };
  // beam sums through the face per unit area (basis e, area ar)
  auto add = [=] KOKKOS_FUNCTION (const Real e[3][3], const Real ar, const int d,
                                  const Real f, Real q[4]) {
    const Real w = dir_(d,3)*f/ar;
    q[0] += w;
    for (int b = 0; b < 3; ++b) {
      q[1+b] += w*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
    }
  };
  // the same with Cartesian momentum, projected and divided by ar once per face (proj)
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
  // a mirror r face, per unit area: no energy, twice the normal momentum of cell ic's
  // beams that hit it (sign sg), no tangential part
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
  // sp with gl: both lateral faces and the r face in one pass over the directions
  // (cln_gd_fs), each sum in the same order as the separate kernels
  const bool fs = gl && x.sph && !x.cub && !axi;
  for (int dir = 1; dir < 3 && !fs; ++dir) {
    const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
    if (!x.sph && dir >= x.nf) continue;
    CleanForK("cln_gd_lat", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      // cubed sphere: also the MeshBlock's high face (CLN_FH of its last cell)
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
        if (!bd && !up) {seam2(dir, m, kf, jf, i, a, e, ar, q);}
        for (int c = 0; c < 4; ++c) cs(m,(up ? CLN_FH + 4*(dir - 1) : nf)+c,k,j,i) = q[c];
      }
    });
  }
  if (!fs) CleanGhostFill(cln_s, CLN_F2, 4);
  if (!fs) CleanGhostFill(cln_s, CLN_F3, 4);
  // rr (spherical polar): r faces carry the rays' own upwind fluxes, not the balance
  // march, which loses the diffusion flux of thick cells in V chi (S - J)
  const bool rr = x.sph;
  auto rflux = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                    Real q[4]) {
    const CleanArgs &x = *px;
    Real a[3], e[3][3], sc[4] = {0.0};
    face(0, m, k, j, i, a, e);
    // (reflecting x1 face: mirror images of the leaving beams, CleanFvMirror, scaled so
    // the face returns the energy it receives)
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
  Kokkos::deep_copy(DevExeSpace(), rf, 0.0);
  if (fs) {
    // (cells is..ie: low lateral faces as cln_gd_lat; faces is..ie+1: the r face as
    // cln_gd_rf, a boundary one by rflux)
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
        if (!b1) {seam2(1, m, k, j, i, a1, e, x.area[1](m,k,j,i), q);}
        for (int c = 0; c < 4; ++c) cs(m,CLN_F2+c,k,j,i) = q[c];
        face(2, m, k, j, i, a2, e);
        proj(e, x.area[2](m,k,j,i), s2, q);
        if (!b2) {seam2(2, m, k, j, i, a2, e, x.area[2](m,k,j,i), q);}
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
  // rr without a reflecting top: no march, one thread per cell; the columns then only
  // close their ends
  const bool cel = rr && !x.cub && ibc_x1max != VET_IBC_REFLECT;
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
        // reflecting top: the remainder at the wall (all but the mirror's momentum)
        // spread over the column by volume, so no energy leaves (lower levels: below)
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
  // reflecting top over several levels: remainder and volumes handed down, each level
  // takes its share
  for (int lev = vgd_nb1 - 2; lev >= 0 && ibc_x1max == VET_IBC_REFLECT; --lev) {
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
//! \fn void RadiationVet::CleanClosureFromRays
//! \brief closure from the vet_sc moments where J >= clean_ray_min E (else isotropic,
//! marked for the fill); the rays' face fluxes copied everywhere

void RadiationVet::CleanClosureFromRays() {
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
  DvceArray1D<Real> nb("cln_raynan", 1);   // (a ray J not finite: counted, no silent D)
  CleanFor("cln_rays", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = vc(m,VETRAY_J,k,j,i);
    if (!Kokkos::isfinite(jr)) Kokkos::atomic_add(&nb(0), 1.0);
    for (int n = 0; n < VETRAY_NFS; ++n) {   // (F2, F1, FX, F3)
      const int c = (n < 4) ? CLN_F2 + n : ((n < 8) ? CLN_F1 + n - 4 : CLN_FX + n - 8);
      if (c < cs.extent_int(1)) cs(m,c,k,j,i) = f2 ? vc(m,VETRAY_F2+n,k,j,i) : 0.0;
    }
    cs(m,CLN_JR,k,j,i) = jr;
    if (jr > 0.0 && jr >= ray_min*iw_(m,VET_IW_EN,k,j,i)) {
      for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = vc(m,VETRAY_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,VETRAY_HP1+a,k,j,i)/jr;
        const Real h = vc(m,VETRAY_H1+a,k,j,i)/jr;
        cs(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cs(m,CLN_HM+a,k,j,i) = fmin(h - hp, 0.0);
      }
      for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
        cs(m,CLN_KP+a,k,j,i) = vc(m,VETRAY_KP+a,k,j,i)/jr;
      }
      cs(m,CLN_C1,k,j,i) = fmax(vc(m,VETRAY_NL,k,j,i)/jr, 0.0);
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
      for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
        cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
      }
      cs(m,CLN_C1,k,j,i) = 1.0;
      cs(m,CLN_J,k,j,i) = 0.0;
      valid(m,0,k,j,i) = 0.0;
    }
  });
  const Real nn = CleanReduce(Kokkos::create_mirror_view_and_copy(HostMemSpace(), nb)(0),
                              false);
  if (nn > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: %.0f cells with a ray J not finite "
                "(isotropic)\n", pmy_pack->pmesh->ncycle, nn);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanClosureFill
//! \brief cells without ray radiation take the mean closure (D, h+-, k+) of their valid
//! neighbours, one layer per pass (clean_fill_passes passes).  clean_stag_thinf: the
//! first layer's, next to the rays, weighted by their E, so that a neighbour the rays
//! hardly light does not count as much as a beam; the deeper layers hold no radiation

void RadiationVet::CleanClosureFill() {
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
  const int nfl = CLN_C1 + std::max(std::min(cln_s.extent_int(1) - CLN_KP, 9), 0);
  // (clean_fill_cs_remap: a neighbour's h+- across a cubed-sphere seam remapped to our
  // axes (CleanCsSlot); families the seam turns, and diagonal neighbours, left out)
  const bool rmp = cln_fillcs && pmy_pack->pmesh->use_cubed_sphere;
  CleanArgs x = CleanMakeArgs();
  auto iw_ = iw;
  const bool thf = cln_thinf;
  const Real efl = fmax(e_floor, 1.0e-300);
  if (thf) {   // the weight, E
    par_for("cln_fill_w", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,2,k,j,i) = fmax(iw_(m,VET_IW_EN,k,j,i), efl);
    });
    CleanGhostFill(cln_v, 2, 1);
  }

  for (int pass = 0; pass < cln_fillpass; ++pass) {
    // (only what cln_fill reads, D, h+- and valid(0), depth 1, in one exchange; callers
    // exchange all of cln_s afterwards; the last pass also valid(1))
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
      Real sum[CLN_C1 + 9] = {0.0}, nh[6] = {0.0}, nv = 0.0;
      for (int kk = k-dk; kk <= k+dk; ++kk) {
        for (int jj = j-dj; jj <= j+dj; ++jj) {
          for (int ii = i-1; ii <= i+1; ++ii) {
            if (valid(m,0,kk,jj,ii) <= 0.0) continue;
            const Real wt = (thf && pass == 0) ? valid(m,2,kk,jj,ii) : 1.0;
            nv += wt;
            // (rmp: the families q[a] of its h+- slots in our axes, sign sg[a], -1: none)
            int q[3] = {0, 1, 2};
            Real sg[3] = {1.0, 1.0, 1.0};
            for (int d = 1; d < 3 && rmp; ++d) {
              const int o = (d == 1) ? jj - j : kk - k;
              if (o == 0) continue;
              Real sd;
              const int n = CleanCsSlot(x, d, m, k + (d == 2 && o > 0),
                                        j + (d == 1 && o > 0), i, kk, jj, ii, false,
                                        sd)/2;
              if (n == d && sd > 0.0) continue;
              const bool dg = (jj != j) && (kk != k);
              q[d] = dg ? -1 : n, sg[d] = sd, q[3 - d] = -1;
            }
            for (int t = 0; t < nfl; ++t) {
              const int a = (t >= CLN_HP && t < CLN_HP + 6) ? (t - CLN_HP)%3 : -1;
              if (a < 0) {
                sum[t] += wt*cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,kk,jj,ii);
                continue;
              }
              if (q[a] < 0) continue;
              const bool hp = t < CLN_HM;   // (a reversed axis: h+ = -h- of the other)
              const int u = ((hp == (sg[a] > 0.0)) ? CLN_HP : CLN_HM) + q[a];
              sum[t] += wt*sg[a]*cs(m,u,kk,jj,ii);
              nh[t - CLN_HP] += wt;
            }
          }
        }
      }
      if (nv == 0.0) return;
      valid(m,1,k,j,i) = 1.0;
      for (int t = 0; t < nfl; ++t) {
        const bool h = rmp && t >= CLN_HP && t < CLN_HP + 6;
        if (h && nh[t - CLN_HP] == 0.0) continue;
        cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,k,j,i) =
            sum[t]/(h ? nh[t - CLN_HP] : nv);
      }
    });
    par_for("cln_fill_mark", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,0,k,j,i) = valid(m,1,k,j,i);
    });
  }
}
#endif  // VET_TU_PART(2)

}  // namespace radvet
