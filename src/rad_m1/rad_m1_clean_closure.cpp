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
#include <limits>
#include <utility>
#include <vector>

#include <Kokkos_Sort.hpp>

#include "athena.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"
#include "utils/imex_sdirk2.hpp"

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
//! to dir (r faces along r_hat, theta faces along theta_hat, phi faces along phi_hat);
//! on the cubed sphere that of CleanCsGeometry; Cartesian the area along e_dir

KOKKOS_INLINE_FUNCTION
void CleanFaceVector(const CleanArgs &x, const DvceArray2D<Real> &x3f, const int dir,
                     const int m, const int k, const int j, const int i, Real a[3]) {
  if (x.cub) {   // the exact closed face vectors of CleanCsGeometry
    for (int c = 0; c < 3; ++c) a[c] = x.area[dir](m,k,j,i)*x.g(m,CLG_N+3*dir+c,k,j,i);
    return;
  }
  if (!x.sph) {
    for (int c = 0; c < 3; ++c) a[c] = (c == dir) ? x.area[dir](m,k,j,i) : 0.0;
    return;
  }
  const Real r1 = x.x1f(m,i), r2 = x.x1f(m,i+1), t1 = x.x2f(m,j), t2 = x.x2f(m,j+1);
  // (a phi face reads only its own phi: the high face of a 2-D cell has no k+1 edge)
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
//! \fn Real CleanFvFace
//! \brief clean_rays = fv: the intensity of direction d on the low face of cell (k,j,i)
//! normal to dir for a beam along +e_dir (pos) or -e_dir: that of the upwind cell, with
//! o2 its minmod slope towards the face (second order in index space, positive), the
//! downwind cell's intensity read from fd; lld (nl directions per block) its trace,
//! with the curvature Q of its source from fd (cln_fvq)

KOKKOS_INLINE_FUNCTION
Real CleanFvFace(const DvceArray5D<Real> &fi, const DvceArray5D<Real> &fd, const bool o2,
                 const int m, const int d, const int dir, const int k, const int j,
                 const int i, const bool pos, const int nl = 0) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  const int su = pos ? -1 : 0, sd = pos ? 0 : -1, sf = pos ? -2 : 1;
  const Real u = fi(m,d,k+su*dk,j+su*dj,i+su*di);
  if (nl > 0) {
    return u + (pos ? 1.0 : -1.0)*fi(m,(1+dir)*nl+d,k+su*dk,j+su*dj,i+su*di) +
           fd(m,CLQ_C+dir,k+su*dk,j+su*dj,i+su*di);
  }
  if (!o2) return u;
  const Real a = u - fi(m,d,k+sf*dk,j+sf*dj,i+sf*di);
  const Real b = fd(m,d,k+sd*dk,j+sd*dj,i+sd*di) - u;
  return u + ((a*b > 0.0) ? 0.5*((fabs(a) < fabs(b)) ? a : b) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFvNbr
//! \brief clean_rays = fv: the intensity of direction d in the cell s cells along dir
//! from (k,j,i) (2-D axisymmetric phi: the own cell rotated, as CleanFvFaceAxi); lld
//! block b at offset ob = b n

KOKKOS_INLINE_FUNCTION
Real CleanFvNbr(const DvceArray5D<Real> &fi, const DvceArray4D<int> &mp, const bool ax,
                const int m, const int d, const int dir, const int k, const int j,
                const int i, const int s, const int ob = 0) {
  if (!ax) return fi(m,ob+d,k+s*(dir == 2),j+s*(dir == 1),i+s*(dir == 0));
  int e = d;
  for (int q = 0; q < ((s < 0) ? -s : s); ++q) e = mp(m,(s < 0) ? 2 : 3,k,e);
  return fi(m,ob+e,k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFvFaceAxi
//! \brief CleanFvFace on a phi face (index k + o) of a 2-D axisymmetric mesh: the cell s
//! columns away is the own cell with the directions rotated by s phi extents (mp: the
//! seam maps cln_fvm, side 2 below, 3 above); the downwind intensity from fd

KOKKOS_INLINE_FUNCTION
Real CleanFvFaceAxi(const DvceArray5D<Real> &fi, const DvceArray5D<Real> &fd,
                    const DvceArray4D<int> &mp, const bool o2, const int m, const int d,
                    const int k, const int j, const int i, const int o, const bool pos,
                    const int nl = 0) {
  auto at = [&](const DvceArray5D<Real> &f, const int s, const int ob = 0) {
    int e = d;
    for (int q = 0; q < ((s < 0) ? -s : s); ++q) e = mp(m, (s < 0) ? 2 : 3, k, e);
    return f(m,ob+e,k,j,i);
  };
  const Real u = at(fi, o - (pos ? 1 : 0));
  if (nl > 0) {
    return u + (pos ? 1.0 : -1.0)*at(fi, o - (pos ? 1 : 0), 3*nl) + fd(m,CLQ_C+2,k,j,i);
  }
  if (!o2) return u;
  const Real a = u - at(fi, o + (pos ? -2 : 1));
  const Real b = at(fd, o - (pos ? 0 : 1)) - u;
  return u + ((a*b > 0.0) ? 0.5*((fabs(a) < fabs(b)) ? a : b) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanFvKey
//! \brief the sweep key n_d . x of cell c from its factors xc (CleanFvTables; cub: the
//! Cartesian position)

template <typename X, typename D>
KOKKOS_INLINE_FUNCTION
Real CleanFvKey(const X &xc, const D &dir, const bool cub, const int m, const int d,
                const int c) {
  if (cub) return dir(d,0)*xc(m,0,c) + dir(d,1)*xc(m,1,c) + dir(d,2)*xc(m,2,c);
  return xc(m,0,c)*(dir(d,0)*xc(m,1,c)*xc(m,2,c) + dir(d,1)*xc(m,1,c)*xc(m,3,c) +
                    dir(d,2)*xc(m,4,c));
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
//! \fn void RadiationM1::CleanFvSourceAlloc
//! \brief the FV rays' source moments in the mesh layout (cln_fvq, CLQ_*)

void RadiationM1::CleanFvSourceAlloc() {
  if (cln_fvq.extent_int(0) == pmy_pack->nmb_thispack) return;
  Kokkos::realloc(cln_fvq, pmy_pack->nmb_thispack, CLQ_N, cln_s.extent_int(2),
                  cln_s.extent_int(3), cln_s.extent_int(4));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvCartSource
//! \brief Cartesian FV rays: the direction set of vet_sc and the source of its sweep,
//! chi and S = eps B + (1 - eps) E (eps = min(kappa_P / chi, 1); the mix of vet_scatter)

void RadiationM1::CleanFvCartSource() {
  if (vgd_n != vet_nray) {
    vgd_dir = vet_ang;
    vgd_n = vet_nray;
  }
  CleanFvSourceAlloc();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool thermal = fl_on && coupling && !opac_zero;
  const int vsf = vscat ? vscat_form : 0;
  const Real ar = arad, efl = e_floor;
  auto iw_ = iw;
  auto op = opac;
  auto fq = cln_fvq;
  par_for("cln_fv_cart", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real chi = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
    const Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl), t2 = SQR(iw_(m,M1_IW_TP,k,j,i));
    const Real eth = fmin(op(m,M1_OP_P,k,j,i)/chi, 1.0);
    Real s = e;
    if (thermal) {
      s = (vsf != 0) ? VetScatterMix(vsf, op(m,M1_OP_P,k,j,i), chi, op(m,M1_OP_S,k,j,i),
                                     ar*t2*t2, e)
                     : eth*ar*t2*t2 + (1.0 - eth)*e;
    }
    fq(m,CLQ_X,k,j,i) = chi;
    fq(m,CLQ_S0,k,j,i) = fmax(s, 1.0e-300);
  });
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
    bool ok = line && ma(m,CLN_J,k,j,i) > 0.0 && q[CLN_J] > 0.0;
    if (ok) {
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
      // (the face sums of a MeshBlock face on the cubed sphere, held on both sides, stay
      // extrapolated alike)
      for (int n = 0; n < (x.cub ? CLN_F2 : CLN_N); ++n) q[n] = mn(m,n,k,j,i);
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
    if (!cub) return;
    const Real dt = 0.5*(1.0 - drr);
    Real nd[3][3];
    CleanCsFamily(x, m, k, j, i, nd);
    const int ab[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    for (int a = 0; a < 6; ++a) {
      cs(m,CLN_D11+a,k,j,i) = (drr - dt)*nd[0][ab[a][0]]*nd[0][ab[a][1]] +
                              ((a < 3) ? dt : 0.0);
    }
    for (int b = 0; b < 3; ++b) {
      cs(m,CLN_KP+b,k,j,i) = krr*nd[0][b];
      for (int d = 1; d < 3; ++d) {
        cs(m,CLN_KP+3*d+b,k,j,i) = kr*nd[0][b] + 0.5*dt*nd[d][b];
      }
    }
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
  if (!vgd_on && cln_rfv == 0) {
    VetShortChar();
    CleanClosureFromRays();
    return;
  }
  if (vgd_on && !tau_ready) VetColInit();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto iw_ = iw;
  auto opac_ = opac;
  par_for("cln_gd_chi", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
  });
  if (cln_rfv > 0) {
    if (vgd_on) {
      VetGdSource();
    } else {
      CleanFvCartSource();
    }
    CleanFvSets(pmy_pack->pmesh->time + pmy_pack->pmesh->dt, pmy_pack->pmesh->ncycle);
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
  const int n = cln_fvn, nb = cln_fvb;
  const Real pp = pmy_pack->pmesh->mesh_size.x3max - pmy_pack->pmesh->mesh_size.x3min;
  const bool circ = std::abs(pp - 2.0*M_PI) < 1.0e-10;
  CleanArgs x = CleanMakeArgs();
  auto fi = cln_fvi;
  auto mp = cln_fvm;
  const int n1 = fi.extent_int(4), n2 = fi.extent_int(3), n3 = fi.extent_int(2);
  if (cln_fvt.extent_int(3) != n) Kokkos::realloc(cln_fvt, n3, n2, n1, n);
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
    for (int b = 0; b < nb; ++b) {
      for (int d = 0; d < n; ++d) t(k,j,i,d) = fi(m,b*n+mp(m,sd,k,d),k,j,i);
      for (int d = 0; d < n; ++d) fi(m,b*n+d,k,j,i) = t(k,j,i,d);
    }
  });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvPanelSlopes
//! \brief cubed sphere, clean_rays = lld: the ghosts beyond a panel seam hold the slopes
//! (q: the source curvatures) along their own panel's x2, x3; they take them along ours,
//! swapped and reversed as the seam maps the indices (GetPanelBoundary)

void RadiationM1::CleanFvPanelSlopes(const bool q) {
  if (!pmy_pack->pmesh->use_cubed_sphere || cln_fvb == 1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng, nl = q ? 1 : cln_fvn, na = std::max(indcs.nx2, indcs.nx3);
  const int o2 = q ? CLQ_C + 1 : 2*cln_fvn, o3 = q ? CLQ_C + 2 : 3*cln_fvn;
  auto a = q ? cln_fvq : cln_fvi;
  auto nb = pmy_pack->pmb->nghbr.d_view;
  auto pnl = pmy_pack->pmb->mb_panel.d_view;
  auto bc = pmy_pack->pmb->mb_bcs.d_view;
  par_for("cln_fv_panel", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, 3, 0, ng - 1,
          0, na - 1, is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int f, const int g, const int t, const int i) {
    // f: the -x2, +x2, -x3, +x3 face of the MeshBlock; the ghost layer g, along it t
    if (bc(m,2+f) != BoundaryFlag::panel) return;
    const int s = (f%2) ? 1 : -1;
    const int nn = (f < 2) ? NeighborIndex(0,s,0,0,0) : NeighborIndex(0,0,s,0,0);
    const int pn = nb(m,nn).panel;
    if (pn == pnl(m)) return;
    const int j = (f < 2) ? ((s > 0) ? je + 1 + g : js - 1 - g) : js + t;
    const int k = (f < 2) ? ks + t : ((s > 0) ? ke + 1 + g : ks - 1 - g);
    if ((f < 2) ? (k > ke) : (j > je)) return;
    const PanelBoundaries pb = GetPanelBoundary(pn, pnl(m));
    for (int d = 0; d < nl; ++d) {
      Real s2 = a(m,o2+d,k,j,i), s3 = a(m,o3+d,k,j,i);
      if (pb.swap_ax) {const Real w = s2; s2 = s3; s3 = w;}
      a(m,o2+d,k,j,i) = (pb.rev_a && !q) ? -s2 : s2;
      a(m,o3+d,k,j,i) = (pb.rev_b && !q) ? -s3 : s3;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::CleanFvWall
//! \brief clean_rays = fv | fv1: whether the low (side 0) or high (side 1) theta face of
//! the mesh is a reflecting wall (specular: CleanFvWalls, maps in CleanFvTables)

bool RadiationM1::CleanFvWall(const int side) {
  Mesh *pm = pmy_pack->pmesh;
  return pm->use_spherical_polar && !pm->use_polar_boundary &&
         pm->mesh_bcs[side ? BoundaryFace::outer_x2 : BoundaryFace::inner_x2] ==
         BoundaryFlag::reflect;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvSets
//! \brief clean_rays = fv | fv1: the swept directions (cln_fvd, cln_fvn) and their tables
//! for the closure of time tc (of the step of cycle cyc).  A fixed vet_gd set is swept as
//! it is.  A rotating one (vet_gd_rotate_dt | _every) is
//! swept as the union of the sets the closure averages, set by set in slots of n
//! directions with the set's weight folded into the quadrature weights, so that every
//! moment and face sum is the weighted mean of the sets' (linear in I: single-valued,
//! conservative and realisable).  With vet_rotate_avg = K these are the last K sets
//! b-K+1 .. b of block b, 1/K each; vet_gd_rotate_xfade adds the next set b+1 and weights
//! the oldest (1 - s)/K and the next s/K, s the fraction of block b passed at tc, so that
//! the closure is continuous in time.  A function of t and the state alone (a
//! restart rebuilds the same union).  Block blk takes slot blk mod (K + xfade); a slot
//! whose set changes rebuilds its tables and starts its sweep from the old set's
//! intensities (time-dependent rays would remap them to the new set here).  Without
//! vet_rotate_advance the K sets repeat: the union is fixed and has nothing to fade.

void RadiationM1::CleanFvSets(const Real tc, const int cyc) {
  const int n = vgd_n;
  if (!VgdRotOn()) {
    if (cln_fvn == n) return;
    cln_fvd = vgd_dir;
    cln_fvn = n;
    CleanFvTables(0, n);
    return;
  }
  const Mesh *pm = pmy_pack->pmesh;
  const int nk = vra_on ? vra_k : 1;
  const bool xf = vgd_xf && !(vra_on && !g_vra_adv);
  const int ns = nk + (xf ? 1 : 0);
  int64_t b = VetGdBlock(cyc);
  Real s = 0.0;
  if (vgd_rot_dt > 0.0) {
    // (without the cross-fade the block of the step's midpoint, as VetGdBlock)
    const Real dt = (pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    const Real tau = (tc - (xf ? 0.0 : 0.5*dt) - vgd_rot_t0)/vgd_rot_dt;
    b = static_cast<int64_t>(std::floor(tau));
    if (xf) s = tau - static_cast<Real>(b);
  }
  if (cln_fvn != ns*n) {
    cln_fvn = ns*n;
    cln_fvd = DvceArray2D<Real>("cln_fvd", ns*n, 4);
    cln_fvkey.assign(ns, std::numeric_limits<Real>::quiet_NaN());
  }
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  std::vector<bool> chg(ns, false);
  for (int64_t blk = b - nk + 1; blk <= b + (xf ? 1 : 0); ++blk) {
    const int q = static_cast<int>(((blk % ns) + ns) % ns);
    const Real key = VetGdAngleBlk(blk);
    const Real w = ((xf && blk == b - nk + 1) ? 1.0 - s : ((blk > b) ? s : 1.0))/nk;
    const std::vector<double> d = VgdRotate(key);
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
    chg[q] = !(key == cln_fvkey[q]);
    cln_fvkey[q] = key;
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  for (int q = 0; q < ns; ++q) {
    if (chg[q]) CleanFvTables(q*n, n);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvTables
//! \brief the tables of the swept directions d0 .. d0+nd-1 (one set; every search stays
//! in it): the cell order of each direction by n . x of the cell centres (centroids on
//! the cubed sphere; cln_fvo), the direction maps of the periodic seams of a wedge
//! (cln_fvm: the ghost's direction in the frame of its image cell, n' = sum_b (n .
//! e_b(ghost)) e_b(image), the nearest set direction; an exact permutation for a phi
//! period of the set) and of the reflecting theta walls (cln_fvwm: an incoming direction
//! takes the leaving one nearest to its mirror image n - 2 (n . t) t about the wall's
//! unit normal t, the mean theta_hat over the face's phi range; cln_fvwu = |n . t|)

void RadiationM1::CleanFvTables(const int d0, const int nd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, ncell = nx1*nx2*indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack, nn = cln_fvn, c3 = cln_s.extent_int(2);
  if (cln_fvo.extent_int(1) != nn) {
    Kokkos::realloc(cln_fvo, nmb, nn, ncell);
    Kokkos::realloc(cln_fvm, nmb, 4, c3, nn);
  }
  cln_fvlv.clear();   // (the sweep's levels follow the order: CleanFvLevels)
  CleanArgs x = CleanMakeArgs();
  auto dir_ = cln_fvd;
  auto x1v = pmy_pack->pcoord->x1v;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto g = cln_gn;
  // the factors of n . x of every cell (r, sin th, cos ph, sin ph, cos th; the centroid)
  DvceArray3D<Real> xc("cln_fvxc", nmb, 5, ncell);
  par_for("cln_fv_xc", DevExeSpace(), 0, nmb - 1, 0, ncell - 1,
  KOKKOS_LAMBDA(const int m, const int c) {
    const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
    const Real f[5] = {x1v(m,i), sin(x2v(m,j)), cos(x3v(m,k)), sin(x3v(m,k)),
                       cos(x2v(m,j))};
    // (the cubed sphere: the centroid; Cartesian: the centre)
    const Real xp[3] = {x1v(m,i), x2v(m,j), x3v(m,k)};
    for (int a = 0; a < 5; ++a) {
      xc(m,a,c) = (a >= 3 || x.sph) ? f[a] : (x.cub ? g(m,CLG_X+a,k,j,i) : xp[a]);
    }
  });
  // the cell order of each direction: n . x ascending, equal keys by the cell index
  auto o_ = cln_fvo;
  const bool cub = !x.sph;
  DvceArray1D<Real> key("cln_fvkey", ncell);
  for (int m = 0; m < nmb; ++m) {
    for (int d = d0; d < d0 + nd; ++d) {
      auto ord = Kokkos::subview(o_, m, d, Kokkos::ALL);
      par_for("cln_fv_key", DevExeSpace(), 0, ncell - 1, KOKKOS_LAMBDA(const int c) {
        key(c) = CleanFvKey(xc, dir_, cub, m, d, c);
        ord(c) = c;
      });
      Kokkos::sort(DevExeSpace(), ord, KOKKOS_LAMBDA(const int a, const int b) {
        return key(a) < key(b) || (key(a) == key(b) && a < b);
      });
    }
  }
  // the seams (sd: theta low, theta high, phi low, phi high; the ghost one cell outside)
  const Real pt = pmy_pack->pmesh->mesh_size.x2max - pmy_pack->pmesh->mesh_size.x2min;
  const Real pp = pmy_pack->pmesh->mesh_size.x3max - pmy_pack->pmesh->mesh_size.x3min;
  auto fm = cln_fvm;
  par_for("cln_fv_seammap", DevExeSpace(), 0, nmb - 1, 0, 3, 0, c3 - 1, d0, d0 + nd - 1,
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int d) {
    if (!x.sph && !x.cub) {   // (Cartesian: the image cell's frame is the own)
      fm(m,sd,k,d) = d;
      return;
    }
    const bool th = (sd < 2), lo = (sd%2 == 0);
    const Real tg = th ? (lo ? x2v(m,js) - 1.0e-3 : x2v(m,je) + 1.0e-3) : 0.5*M_PI;
    const Real pg = th ? x3v(m,k) : 0.0;
    Real eg[3][3], ei[3][3], v[3] = {0.0, 0.0, 0.0};
    CleanLocalBasis(tg, pg, eg);
    CleanLocalBasis(th ? tg + (lo ? pt : -pt) : tg, th ? pg : pg + (lo ? pp : -pp), ei);
    for (int b = 0; b < 3; ++b) {
      const Real nb = dir_(d,0)*eg[b][0] + dir_(d,1)*eg[b][1] + dir_(d,2)*eg[b][2];
      for (int c = 0; c < 3; ++c) v[c] += nb*ei[b][c];
    }
    Real bd = -2.0;
    for (int e = d0; e < d0 + nd; ++e) {
      const Real dt = v[0]*dir_(e,0) + v[1]*dir_(e,1) + v[2]*dir_(e,2);
      if (dt > bd) {bd = dt; fm(m,sd,k,d) = e;}
    }
  });
  if (!CleanFvWall(0) && !CleanFvWall(1)) return;
  if (cln_fvwm.extent_int(3) != nn) {
    Kokkos::realloc(cln_fvwm, nmb, 2, c3, nn);
    Kokkos::realloc(cln_fvwu, nmb, 2, c3, nn);
    Kokkos::realloc(cln_fvw, nmb, 2, c3, cln_s.extent_int(4), cln_fvb*nn);
  }
  auto x2f = pmy_pack->pcoord->xx2f;
  auto x3f = pmy_pack->pcoord->xx3f;
  auto wm = cln_fvwm;
  auto wu = cln_fvwu;
  par_for("cln_fv_wallmap", DevExeSpace(), 0, nmb - 1, 0, 1, ks, ke, d0, d0 + nd - 1,
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int d) {
    const Real th = x2f(m, sd ? je + 1 : js), sgn = sd ? 1.0 : -1.0;
    const Real p1 = x3f(m,k), p2 = x3f(m,k+1);
    Real t[3] = {cos(th)*(sin(p2) - sin(p1)), cos(th)*(cos(p1) - cos(p2)),
                 -sin(th)*(p2 - p1)};
    const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    for (int c = 0; c < 3; ++c) t[c] /= tn;
    const Real nt = dir_(d,0)*t[0] + dir_(d,1)*t[1] + dir_(d,2)*t[2];
    wu(m,sd,k,d) = fabs(nt);
    wm(m,sd,k,d) = -1;
    if (!(sgn*nt < 0.0)) return;
    const Real v[3] = {dir_(d,0) - 2.0*nt*t[0], dir_(d,1) - 2.0*nt*t[1],
                       dir_(d,2) - 2.0*nt*t[2]};
    Real bd = -2.0;
    for (int e = d0; e < d0 + nd; ++e) {
      const Real et = dir_(e,0)*t[0] + dir_(e,1)*t[1] + dir_(e,2)*t[2];
      const Real dt = v[0]*dir_(e,0) + v[1]*dir_(e,1) + v[2]*dir_(e,2);
      if (sgn*et > 0.0 && dt > bd) {bd = dt; wm(m,sd,k,d) = e;}
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvWalls
//! \brief the incoming intensities of the reflecting theta walls (cln_fvw): each incoming
//! direction e takes the intensity of the wall cell along its mirror direction, scaled
//! on each wall face so that it returns exactly the energy it receives,
//! sum_e w_e |n_e . t| I(e) = sum_out w_d |n_d . t| I_d (isotropic radiation stays
//! isotropic: the set integrates n exactly); lld: the wall face traces, with the
//! slopes along x1 and x3 (blocks 1, 3)

void RadiationM1::CleanFvWalls() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, n = cln_fvn, nl = (cln_fvb > 1) ? n : 0;
  CleanArgs x = CleanMakeArgs();
  auto fi = cln_fvi;
  auto fq = cln_fvq;
  auto fw = cln_fvw;
  auto mp = cln_fvwm;
  auto u = cln_fvwu;
  auto dir_ = cln_fvd;
  par_for("cln_fv_wall", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, 1, indcs.ks,
          indcs.ke, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int i) {
    const int jc = sd ? je : js;
    if (!CleanBoundaryFace(x, 1, m, k, sd ? je + 1 : js, i)) return;
    Real so = 0.0, si = 0.0;
    for (int e = 0; e < n; ++e) {
      const int d = mp(m,sd,k,e), g = (d >= 0) ? d : e;
      const Real tr = fi(m,g,k,jc,i) + ((nl > 0) ? (sd ? 1.0 : -1.0)*fi(m,2*n+g,k,jc,i) +
                                                   fq(m,CLQ_C+1,k,jc,i) : 0.0);
      const Real q = dir_(e,3)*u(m,sd,k,e)*tr;
      fw(m,sd,k,i,e) = (d >= 0) ? tr : 0.0;
      for (int b = 1; b < 4 && nl > 0; b += 2) {
        fw(m,sd,k,i,b*n+e) = (d >= 0) ? fi(m,b*n+d,k,jc,i) : 0.0;
      }
      if (d >= 0) {si += q;} else {so += q;}
    }
    for (int e = 0; e < n; ++e) {
      for (int b = 0; b < 4; ++b) {
        if (b == 0 || (nl > 0 && b != 2)) fw(m,sd,k,i,b*n+e) *= (si > 0.0) ? so/si : 0.0;
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvLevels
//! \brief the wavefronts of the sweep: in each direction's order (cln_fvo) a cell takes
//! the level after every face neighbour placed before it across a face the direction
//! crosses (fv: every cell up to two along an axis), so that a level's cells of all
//! directions read nothing written in the level and a pass gives the sequential
//! sweep's intensities; cln_fvl lists (m n + d, cell) by level and cell, cln_fvlv the
//! level starts.  One team per direction takes the order in chunks; the levels of a
//! chunk are relaxed from below until they no longer change (the longest path).

void RadiationM1::CleanFvLevels() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3, nc = nx1*nx2*nx3;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nmb = pmy_pack->nmb_thispack, n = cln_fvn, r2 = (cln_rfv == 2) ? 2 : 1;
  const int nmd = nmb*n, nb = 1024;
  const bool axi = !pmy_pack->pmesh->three_d && pmy_pack->pmesh->use_spherical_polar;
  auto o_ = cln_fvo;
  auto fa = cln_fva;
  auto dir_ = cln_fvd;
  DvceArray2D<int> pos("cln_fvpos", nmd, nc), lev("cln_fvlev", nmd, nc);
  par_for("cln_fvl_pos", DevExeSpace(), 0, nmd - 1, 0, nc - 1,
  KOKKOS_LAMBDA(const int md, const int p) {
    pos(md,o_(md/n,md%n,p)) = p;
    lev(md,p) = 0;
  });
  int lmax = 0;
  Kokkos::parallel_reduce("cln_fvl_lev", Kokkos::TeamPolicy<>(DevExeSpace(), nmd,
                                                              Kokkos::AUTO),
  KOKKOS_LAMBDA(const TeamMember_t &tm, int &lm) {
    const int md = tm.league_rank(), m = md/n, d = md%n;
    const int nx[3] = {nx1, nx2, nx3};
    for (int p0 = 0; p0 < nc; p0 += nb) {
      const int p1 = (p0 + nb < nc) ? p0 + nb : nc;
      int chg = 1;
      while (chg > 0) {
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, p0, p1),
        [&](const int p, int &ch) {
          const int c = o_(m,d,p), x[3] = {c%nx1, (c/nx1)%nx2, c/(nx1*nx2)};
          int l = 0;
          for (int f = 0; f < 6; ++f) {
            const int e = f/2, s = (f%2) ? 1 : -1;
            const Real a = dir_(d,0)*fa(m,3*f,ks+x[2],js+x[1],is+x[0]) +
                           dir_(d,1)*fa(m,3*f+1,ks+x[2],js+x[1],is+x[0]) +
                           dir_(d,2)*fa(m,3*f+2,ks+x[2],js+x[1],is+x[0]);
            if ((axi && e == 2) || (a == 0.0 && r2 == 1)) continue;
            for (int r = 1; r <= r2; ++r) {
              int y[3] = {x[0], x[1], x[2]};
              y[e] += r*s;
              if (y[e] < 0 || y[e] >= nx[e]) break;
              const int g = y[0] + nx1*(y[1] + nx2*y[2]);
              if (pos(md,g) < p) l = (l > lev(md,g) + 1) ? l : lev(md,g) + 1;
            }
          }
          if (l != lev(md,c)) {
            lev(md,c) = l;
            ++ch;
          }
        }, chg);
      }
    }
    int l1 = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nc), [&](const int c, int &l) {
      l = (l > lev(md,c) + 1) ? l : lev(md,c) + 1;
    }, Kokkos::Max<int>(l1));
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {lm = (lm > l1) ? lm : l1;});
  }, Kokkos::Max<int>(lmax));
  // the list sorted by (level, m n + d, cell): one key per entry
  const int64_t nl = static_cast<int64_t>(nmd)*nc;
  DvceArray1D<int64_t> key("cln_fvl_key", nl);
  par_for("cln_fvl_key", DevExeSpace(), 0, nmd - 1, 0, nc - 1,
  KOKKOS_LAMBDA(const int md, const int c) {
    const int64_t q = md*static_cast<int64_t>(nc) + c;
    key(q) = lev(md,c)*nl + q;
  });
  pos = DvceArray2D<int>();
  lev = DvceArray2D<int>();
  Kokkos::sort(DevExeSpace(), key);
  Kokkos::realloc(cln_fvl, nl, 2);
  DvceArray1D<int> st("cln_fvl_start", lmax);
  auto fl = cln_fvl;
  Kokkos::parallel_for("cln_fvl_list", Kokkos::RangePolicy<>(DevExeSpace(), 0, nl),
  KOKKOS_LAMBDA(const int q) {
    const int l = static_cast<int>(key(q)/nl), mc = static_cast<int>(key(q)%nl);
    fl(q,0) = mc/nc;
    fl(q,1) = mc%nc;
    if (q == 0 || key(q-1)/nl != l) st(l) = q;
  });
  auto st_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), st);
  cln_fvlv.assign(lmax + 1, static_cast<int>(nl));
  for (int l = 0; l < lmax; ++l) cln_fvlv[l] = st_h(l);
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
//! cell and the boundary beams (implicit_beam), the top its Marshak bath, a reflecting
//! theta wall is specular (CleanFvWalls), another physical theta or phi (Cartesian x2,
//! x3) face is closed.  Cartesian meshes take vet_sc's direction set.

void RadiationM1::CleanFvSweep() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const int n = cln_fvn, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  const bool axi = !pmy_pack->pmesh->three_d && pmy_pack->pmesh->use_spherical_polar;
  CleanArgs x = CleanMakeArgs();
  auto x3f = pmy_pack->pcoord->xx3f;
  if (cln_fvi.extent_int(1) != cln_fvb*n) {
    Kokkos::realloc(cln_fvi, nmb, cln_fvb*n, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
  }
  if (cln_fva.extent_int(0) != nmb) {
    // the outward face vectors of every cell (CleanFvSets: the directions and tables)
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
  }
  if (cln_fvlv.empty()) CleanFvLevels();
  auto fi = cln_fvi;
  auto fa = cln_fva;
  auto fl = cln_fvl;
  auto fm = cln_fvm;
  const bool o2 = (cln_rfv == 2);
  const bool wall[2] = {CleanFvWall(0), CleanFvWall(1)};
  auto fw = cln_fvw;
  auto dir_ = cln_fvd;
  auto cw_ = vgd_cs;
  const bool emi = VetBcEmit();
  const Real iemit = VetBathLo();
  const Real ihi = (ibc_x1max == M1_IBC_MARSHAK) ? iebath_x1max : 0.0;
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
  // fv: the downwind cell of a slope is read from the last build (the own cell's outflow
  // slopes are implicit), so that a pass is a complete upwind sweep as in fv1; the slopes
  // settle with the run, exactly at a steady state
  if (o2 && cln_fvp.extent_int(1) != n) {
    Kokkos::realloc(cln_fvp, fi.extent_int(0), n, fi.extent_int(2), fi.extent_int(3),
                    fi.extent_int(4));
  }
  if (o2) Kokkos::deep_copy(cln_fvp, fi);
  auto fp = o2 ? cln_fvp : fi;
  // 2-D axisymmetric: a phi face's neighbour is the own cell in other directions, which
  // other threads write in the same pass: they are read from the last pass (cln_fvx)
  if (axi && cln_fvx.extent_int(1) != cln_fvb*n) {
    Kokkos::realloc(cln_fvx, fi.extent_int(0), cln_fvb*n, fi.extent_int(2),
                    fi.extent_int(3), fi.extent_int(4));
  }
  auto fx = axi ? cln_fvx : fi;
  // lld: the source moments.  S = q + (1 - eps) J, J and the emission q of the gas each
  // with face values interpolated linearly in optical depth between the cell centres
  // (single valued, exact for profiles linear in tau), the slope (S_f+ - S_f-)/2 and the
  // curvature Q = (S_f+ + S_f-)/2 - S that the outflow traces carry in the diffusion
  // limit, so that they match the face values: weighted by w = V chi / (V chi + sum
  // |A_f| / 4) (the outflow area of the mean direction) of the cell and of its thinner
  // neighbour along the axis (none at a photosphere, where the field decouples).  Where
  // the central slope of q exceeds its smaller one-sided difference d (an edge of the
  // gas) its slope and curvature take the factor d/|c| (none at an extremum); one sided
  // at a boundary
  const bool ld = (cln_rfv == 3);
  const int nl = ld ? n : 0;
  if (vgd_on) {   // chi and S in the mesh layout (Cartesian: CleanFvCartSource)
    CleanFvSourceAlloc();
    auto fq = cln_fvq;
    par_for("cln_fv_cs", DevExeSpace(), 0, nmb - 1, ks, indcs.ke, js, indcs.je, is,
            indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fq(m,CLQ_X,k,j,i) = exp(cw_(m,0,k+ok,j+og,i));
      fq(m,CLQ_S0,k,j,i) = exp(cw_(m,1,k+ok,j+og,i));
    });
  }
  auto fq = cln_fvq;
  if (ld) {
    auto iw_ = iw;
    auto op = opac;
    par_for("cln_fv_sc", DevExeSpace(), 0, nmb - 1, ks, indcs.ke, js, indcs.je, is,
            indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real chi = fq(m,CLQ_X,k,j,i), sj = fq(m,CLQ_S0,k,j,i);
      const Real ej = fmax(1.0 - op(m,M1_OP_P,k,j,i)/chi, 0.0);
      const Real jr = fmin(fmax(iw_(m,M1_IW_EN,k,j,i), 0.0), (ej > 0.0) ? sj/ej : 0.0);
      fq(m,CLQ_J,k,j,i) = jr;
      fq(m,CLQ_E,k,j,i) = sj - ej*jr;
      fq(m,CLQ_W,k,j,i) = ej;
    });
    CleanGhostFill(cln_fvq, CLQ_X, 3);
    par_for("cln_fv_src", DevExeSpace(), 0, nmb - 1, ks, indcs.ke, js, indcs.je, is,
            indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real ar = 0.0;
      for (int f = 0; f < 6; ++f) {
        ar += sqrt(SQR(fa(m,3*f,k,j,i)) + SQR(fa(m,3*f+1,k,j,i)) +
                   SQR(fa(m,3*f+2,k,j,i)));
      }
      const Real vc = x.vol(m,k,j,i)*fq(m,CLQ_X,k,j,i), wq = vc/(vc + 0.25*ar);
      for (int e = 0; e < 3; ++e) {
        const int di = (e == 0), dj = (e == 1), dk = (e == 2);
        Real sl = 0.0, cv = 0.0, wn = 1.0;
        for (int c = CLQ_J; c <= CLQ_E; ++c) {
          const Real v0 = fq(m,c,k,j,i), wc = (c == CLQ_J) ? fq(m,CLQ_W,k,j,i) : 1.0;
          Real vf[2] = {v0, v0}, dq[2] = {0.0, 0.0};
          int nf = 0;
          for (int up = 0; up < 2; ++up) {
            const int s = up ? 1 : -1, kn = k + s*dk, jn = j + s*dj, in = i + s*di;
            if ((e == 1 && !x.multi_d) || (e == 2 && !x.three_d) ||
                CleanBoundaryFace(x, e, m, k + up*dk, j + up*dj,
                                                      i + up*di)) continue;
            const Real d0 = (e == 0) ? fabs(x.x1f(m,i+up) - x.x1v(m,i)) : 1.0;
            const Real dn = (e == 0) ? fabs(x.x1v(m,in) - x.x1f(m,i+up)) : 1.0;
            const Real t0 = fq(m,CLQ_X,k,j,i)*d0, tn = fq(m,CLQ_X,kn,jn,in)*dn;
            vf[up] = (tn*v0 + t0*fq(m,c,kn,jn,in))/(t0 + tn);
            const Real vn = x.vol(m,k,j,i)*fq(m,CLQ_X,kn,jn,in);
            wn = fmin(wn, vn/(vn + 0.25*ar));
            dq[up] = s*(fq(m,c,kn,jn,in) - v0);
            ++nf;
          }
          Real q = (vf[1] - vf[0])*((nf == 2) ? 0.5 : 1.0);
          const Real dm = fmin(fabs(dq[0]), fabs(dq[1]));
          const Real ph = (nf < 2) ? 0.0 : (fabs(q) > dm) ? dm/fabs(q) : 1.0;
          cv += wc*((c == CLQ_J) ? 1.0 : ph)*(0.5*(vf[0] + vf[1]) - v0);
          if (c == CLQ_E && nf == 2) q *= (dq[0]*dq[1] > 0.0 || ph == 1.0) ? ph : 0.0;
          sl += wc*q;
        }
        fq(m,CLQ_S+e,k,j,i) = sl;
        fq(m,CLQ_C+e,k,j,i) = wq*wn*cv;
      }
    });
    CleanGhostFill(cln_fvq, CLQ_C, 3);
    CleanFvPanelSlopes(true);
  }
  // the walls' incoming beams of the intensities the sweep starts from (a restart has
  // no others)
  if (wall[0] || wall[1]) CleanFvWalls();
  // (each level reduces -(its largest change) and its largest intensity on the device:
  // no wait)
  const int nlev = static_cast<int>(cln_fvlv.size()) - 1;
  Kokkos::View<Kokkos::MinMaxScalar<Real>*, DevMemSpace> red("cln_fv_red", nlev);
  // fv1: a pass after the first takes only the cells at the faces of the MeshBlock and
  // those whose upwind neighbour changed after it took them (due: the pass that must
  // take a cell, done: the pass that has); the others would repeat their intensity
  const bool act = !o2 && !ld && !axi;
  const int ie = indcs.ie, je = indcs.je, ke = indcs.ke, nc = nx1*nx2*indcs.nx3;
  if (act && cln_fvdue.extent_int(1) != nc) {
    Kokkos::realloc(cln_fvdue, nmb*n, nc);
    Kokkos::realloc(cln_fvdone, nmb*n, nc);
  }
  auto due = cln_fvdue;
  auto done = cln_fvdone;
  for (int it = 0; it < 100; ++it) {
    Real dmax = 0.0, imax = 0.0;
    if (axi) Kokkos::deep_copy(cln_fvx, fi);
    if (act) Kokkos::deep_copy(done, -1);
    for (int l = 0; l < nlev; ++l) {
    Kokkos::parallel_reduce("cln_fv", Kokkos::RangePolicy<>(DevExeSpace(), cln_fvlv[l],
                                                            cln_fvlv[l+1]),
    KOKKOS_LAMBDA(const int q, Kokkos::MinMaxScalar<Real> &mm) {
      const int md = fl(q,0), m = md/n, d = md%n, c = fl(q,1);
      const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
      {
        const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
        if (act && it > 0 && due(md,c) != it && i != is && i != ie && j != js &&
            j != je && k != ks && k != ke) {
          done(md,c) = it;
          mm.max_val = fmax(mm.max_val, fi(m,d,k,j,i));
          return;
        }
        // fv: the outflow faces carry the cell's own minmod slope, I + mm(I - U, D - I)/2
        // (U, D the neighbours up- and downstream), solved with the cell: the residual
        // increases with I (slope >= (out + vc)/2), a few Newton steps on its pieces
        // lld: per axis the outflow weights cd, the inflow weights hd and inflow rd of
        // the slope equations, outflow sides po (+1, high face; 2, low face)
        Real in = 0.0, out = 0.0, ou[6], ud[6][2], cd[3] = {0.0, 0.0, 0.0};
        Real hd[3] = {0.0, 0.0, 0.0}, rd[3] = {0.0, 0.0, 0.0}, od[3] = {0.0, 0.0, 0.0};
        int no = 0, po[3] = {0, 0, 0}, om = 0;
        for (int f = 0; f < 6; ++f) {
          const int dir = f/2, up = f%2;
          const int kf = k + (up && dir == 2), jf = j + (up && dir == 1);
          const int i_f = i + (up && dir == 0);
          const Real a = nv[0]*fa(m,3*f,k,j,i) + nv[1]*fa(m,3*f+1,k,j,i) +
                         nv[2]*fa(m,3*f+2,k,j,i);
          // (2-D axisymmetric: the phi faces open, their neighbours the own cell)
          const bool ax = axi && dir == 2;
          // (Cartesian: no faces along an axis of one cell)
          if (!x.sph && !x.cub && (dir == 1 ? !x.multi_d : dir == 2 && !x.three_d)) {
            continue;
          }
          const bool bdry = !ax && CleanBoundaryFace(x, dir, m, kf, jf, i_f);
          // a reflecting theta wall: specular (its incoming beams from cln_fvw)
          const int sw = (dir == 1 && bdry && wall[up]) ? up : -1;
          if (bdry && dir != 0 && sw < 0) continue;
          const Real sg = up ? 1.0 : -1.0;
          if (sw >= 0 && a < 0.0) {
            in -= a*fw(m,sw,k,i,d);
            for (int e = 0; e < 3 && ld; ++e) {
              rd[e] -= a*((e == 1) ? sg*fw(m,sw,k,i,d) : fw(m,sw,k,i,(1+e)*n+d));
            }
            hd[1] -= sg*a;
            continue;
          }
          if (a >= 0.0) {
            if (!bdry) om |= 1 << f;
            out += a;
            cd[dir] += sg*a;
            od[dir] += a;
            po[dir] |= up ? 1 : 2;
            if (o2) {
              const int sg = up ? 1 : -1;
              ou[no] = a;
              ud[no][0] = CleanFvNbr(ax ? fx : fi, fm, ax, m, d, dir, k, j, i, -sg);
              ud[no][1] = CleanFvNbr(fp, fm, ax, m, d, dir, k, j, i, sg);
              ++no;
            }
          } else {
            // the upwind trace (lld: the neighbour's transverse slopes too)
            auto &fd = ld ? fq : fp;
            Real t = bdry ? (up ? ihi : (emi ? iemit : fq(m,CLQ_S0,k,j,i)))
                   : ax ? CleanFvFaceAxi(fx, fd, fm, o2, m, d, k, j, i, up, !up, nl)
                   : CleanFvFace(fi, fd, o2, m, d, dir, kf, jf, i_f, !up, nl);
            for (int q = 0; q < bmn && bdry && !up; ++q) {
              if (d == bmr[q]) t += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
            }
            in -= a*t;
            hd[dir] -= sg*a;
            for (int e = 0; e < 3 && ld; ++e) {
              rd[e] -= a*((e == dir) ? sg*t : bdry ? 0.0
                          : CleanFvNbr(ax ? fx : fi, fm, ax, m, d, dir, k, j, i, up ? 1 : -1,
                                       (1+e)*n));
            }
          }
        }
        const Real vc = x.vol(m,k,j,i)*fq(m,CLQ_X,k,j,i), sc = vc*fq(m,CLQ_S0,k,j,i);
        Real iv = (out + vc > 0.0) ? (in + sc)/(out + vc) : 0.0;
        // lld: the outflow traces I + s I_e + Q_e (Q_e the source's curvature along e),
        // the lumped slope equations hd I + T I_e = rd + V chi S_e - cd Q_e in the
        // balance T I + sum cd I_e = in + V chi S - sum od Q_e (T = out + V chi); where
        // the average or an outflow trace would be negative the slopes are scaled
        // towards the step value
        Real sl[3] = {0.0, 0.0, 0.0};
        if (ld && out + vc > 0.0) {
          const Real t = out + vc;
          Real num = in + sc, den = t, th = 1.0, qe[3];
          for (int e = 0; e < 3; ++e) {
            qe[e] = fq(m,CLQ_C+e,k,j,i);
            num -= od[e]*qe[e];
          }
          iv = num/t;
          for (int e = 0; e < 3; ++e) {
            rd[e] += vc*fq(m,CLQ_S+e,k,j,i) - cd[e]*qe[e];
            num -= cd[e]*rd[e]/t;
            den -= cd[e]*hd[e]/t;
          }
          const Real i0 = num/den;
          for (int e = 0; e < 3; ++e) sl[e] = (rd[e] - hd[e]*i0)/t;
          auto lim = [&](const Real b, const Real q) {
            if (q < 0.0) th = fmin(th, (b > 0.0) ? b/(b - q) : 0.0);
          };
          lim(iv, i0);
          for (int e = 0; e < 3; ++e) {
            if (po[e] & 1) lim(iv + qe[e], i0 + sl[e] + qe[e]);
            if (po[e] & 2) lim(iv + qe[e], i0 - sl[e] + qe[e]);
          }
          iv += th*(i0 - iv);
          for (int e = 0; e < 3; ++e) sl[e] *= th;
        }
        for (int q = 0; q < 8 && no > 0 && out + vc > 0.0; ++q) {
          Real g = (out + vc)*iv - in - sc, dg = out + vc;
          for (int o = 0; o < no; ++o) {
            const Real p1 = iv - ud[o][0], p2 = ud[o][1] - iv;
            if (p1*p2 > 0.0) {
              const bool sm = fabs(p1) < fabs(p2);
              g += 0.5*ou[o]*(sm ? p1 : p2);
              dg += 0.5*ou[o]*(sm ? 1.0 : -1.0);
            }
          }
          const Real dv = g/dg;
          iv -= dv;
          if (fabs(dv) <= 1.0e-15*fabs(iv)) break;
        }
        iv = fmax(iv, 0.0);
        mm.min_val = fmin(mm.min_val, -fabs(iv - fi(m,d,k,j,i)));
        mm.max_val = fmax(mm.max_val, iv);
        // (a change: the downstream neighbours are due in this pass, or the next one if
        // this pass has taken them)
        const bool chg = act && !(iv == fi(m,d,k,j,i));
        fi(m,d,k,j,i) = iv;
        if (act) done(md,c) = it;
        for (int f = 0; f < 6 && chg; ++f) {
          const int x0 = (f < 2) ? i - is : ((f < 4) ? j - js : k - ks);
          const int x1 = (f < 2) ? ie - is : ((f < 4) ? je - js : ke - ks);
          if (!(om & (1 << f)) || x0 == ((f%2) ? x1 : 0)) continue;   // (in the block)
          const int s = (f%2) ? 1 : -1, st = (f < 2) ? 1 : ((f < 4) ? nx1 : nx1*nx2);
          const int cn = c + s*st;
          due(md,cn) = (done(md,cn) == it) ? it + 1 : it;
        }
        for (int e = 0; e < 3 && ld; ++e) fi(m,(1+e)*n+d,k,j,i) = sl[e];
      }
    }, Kokkos::MinMax<Real, DevMemSpace>(Kokkos::subview(red, l)));
    }
    auto r_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), red);
    for (int l = 0; l < nlev; ++l) {
      dmax = fmax(dmax, -r_h(l).min_val);
      imax = fmax(imax, r_h(l).max_val);
    }
    CleanGhostFill(cln_fvi, 0, cln_fvb*n);
    CleanFvPanelSlopes(false);
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
//! CleanClosureFromRays; no ray face fluxes).  On the cubed sphere mu_b = n_b (Cartesian
//! D and k+ columns) and the half ranges are about the family normals (CleanCsFamily).

void RadiationM1::CleanClosureFromGd() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool fv = (cln_rfv > 0);
  const int n = fv ? cln_fvn : vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
  auto dir_ = fv ? cln_fvd : vgd_dir;
  auto cw_ = vgd_cs;
  auto fq = cln_fvq;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;
  CleanArgs x = CleanMakeArgs();
  par_for("cln_gd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real st = sin(x2v(m,j)), ct = cos(x2v(m,j));
    const Real sp = sin(x3v(m,k)), cp = cos(x3v(m,k));
    const Real s0 = fv ? fq(m,CLQ_S0,k,j,i) : exp(cw_(m,1,k+ok,j+og,i));
    // the cubed sphere: D in Cartesian components, half ranges about the families n_d
    Real nd[3][3];
    if (x.cub) CleanCsFamily(x, m, k, j, i, nd);
    Real q[CLN_N] = {0.0};
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      Real mu[3] = {nx*st*cp + ny*st*sp + nz*ct, nx*ct*cp + ny*ct*sp - nz*st,
                    -nx*sp + ny*cp}, nu[3];
      for (int b = 0; b < 3; ++b) {
        if (!x.sph) mu[b] = dir_(d,b);
        nu[b] = x.cub ? nx*nd[b][0] + ny*nd[b][1] + nz*nd[b][2] : mu[b];
      }
      const Real iv = fv ? fi_(m,d,k,j,i) : vi_(m,d,k+ok,j+og,i);
      const Real a = dir_(d,3)*iv;
      q[CLN_J] += a;
      for (int b = 0; b < 3; ++b) {
        q[CLN_D11+b] += a*mu[b]*mu[b];
        q[((nu[b] > 0.0) ? CLN_HP : CLN_HM)+b] += a*nu[b];
        if (!(nu[b] > 0.0)) continue;
        for (int c = 0; c < 3; ++c) q[CLN_KP+3*b+c] += a*mu[c]*nu[b];
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
//! On the cubed sphere the momentum is Cartesian and the r faces carry the rays' upwind
//! fluxes (no curvature terms).

void RadiationM1::CleanGdFaceSums() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool fv = (cln_rfv > 0), o2 = (cln_rfv == 2);
  const int n = fv ? cln_fvn : vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
  auto fm_ = cln_fvm;
  const int nl = (cln_rfv == 3) ? n : 0;
  auto fq_ = (nl > 0) ? cln_fvq : cln_fvi;
  auto sq_ = cln_fvq;
  auto dir_ = fv ? cln_fvd : vgd_dir;
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
  const bool axi = !pmy_pack->pmesh->three_d && x.sph;
  auto bnd = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i) {
    return !(axi && dir == 2) && CleanBoundaryFace(x, dir, m, k, j, i);
  };
  // the inner face: no flux at a reflecting one, else the outward beams of the sweep's
  // emission (vet_bc_inner emit) or of the first cell, the inward ones of the first cell
  const bool open = (ibc_x1min != M1_IBC_REFLECT);
  const bool emi = VetBcEmit();
  const Real iemit = VetBathLo();
  const Real ihi = (ibc_x1max == M1_IBC_MARSHAK) ? iebath_x1max : 0.0;
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
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
  // (Cartesian on the cubed sphere)
  auto face = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                   const int i, Real a[3], Real e[3][3]) {
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
  // beam of direction d through the low theta or phi face of cell (k,j,i)
  auto lat = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d) {
    if (bnd(dir, m, k, j, i)) return 0.0;
    Real a[3], e[3][3];
    face(dir, m, k, j, i, a, e);
    const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
    const int kl = k - (dir == 2), jl = j - (dir == 1);
    if (fv && axi && dir == 2) {
      return f*CleanFvFaceAxi(fi_, fq_, fm_, o2, m, d, k, j, i, 0, f > 0.0, nl);
    }
    if (fv) return f*CleanFvFace(fi_, fq_, o2, m, d, dir, k, j, i, f > 0.0, nl);
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
    Real a[3], e[3][3], w[4] = {0.0}, nf[3], wn = 0.0;
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      const Real tr = inten(m, d, k, j, ic) +
                      ((nl > 0) ? sg*fi_(m,nl+d,k,j,ic) + fq_(m,CLQ_C,k,j,ic) : 0.0);
      if (sg*f > 0.0) add(0, m, k, j, i, d, 2.0*f*tr, w);
    }
    for (int b = 0; b < 3; ++b) nf[b] = CleanNf(x, b, 0, m, k, j, i);
    for (int b = 0; b < 3; ++b) wn += w[1+b]*nf[b];
    q[0] = 0.0;
    for (int b = 0; b < 3; ++b) q[1+b] = wn*nf[b];
  };
  for (int dir = 1; dir < 3; ++dir) {
    const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
    if (!x.sph && dir >= x.nf) continue;
    par_for("cln_gd_lat", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // the cubed sphere: the high face of the MeshBlock too (CLN_FH of its last cell)
      const bool hi = x.cub && ((dir == 1) ? j == je : k == ke);
      for (int up = 0; up <= (hi ? 1 : 0); ++up) {
        const int kf = k + up*(dir == 2), jf = j + up*(dir == 1);
        Real q[4] = {0.0};
        for (int d = 0; d < n; ++d) {
          add(dir, m, kf, jf, i, d, lat(dir, m, kf, jf, i, d), q);
        }
        for (int c = 0; c < 4; ++c) cs(m,(up ? CLN_FH + 4*(dir - 1) : nf)+c,k,j,i) = q[c];
      }
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
      Real iv = CleanFvFace(fi_, fq_, o2, m, d, 0, k, j, i, f > 0.0, nl);
      if (CleanBoundaryFace(x, 0, m, k, j, i)) {   // the sweep's x1 boundary
        iv = (i == is) ? ((f > 0.0) ? (emi ? iemit : sq_(m,CLQ_S0,k,j,is)) : iv)
                       : ((f > 0.0) ? iv : ihi);
        for (int q = 0; q < bmn && i == is && f > 0.0; ++q) {
          if (d == bmr[q]) iv += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
        }
      }
      add(0, m, k, j, i, d, f*iv, q);
    }
    for (int c = 0; c < 4; ++c) q[c] *= x.area[0](m,k,j,i);
  };
  if (!x.sph && !x.cub) {   // Cartesian: every x1 face carries the rays' upwind flux
    par_for("cln_gd_x1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4];
      for (int up = 0; up <= ((i == ie) ? 1 : 0); ++up) {
        rflux(m, k, j, i + up, q);
        for (int c = 0; c < 4; ++c) {
          const Real f = q[c]/x.area[0](m,k,j,i+up);
          if (!up) cs(m,CLN_FX+c,k,j,i) = f;
          if ((!up && i == is) || up) cs(m,CLN_F1+c,k,j,i) = f;
        }
      }
      for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = 0.0;
    });
    return;
  }
  auto xf_ = vgd_x1f;
  auto qb_ = vgd_x1b;
  auto qr_ = vgd_x1r;
  const int nj = je - js + 1, nkj = (ke - ks + 1)*nj;
  // fv: the rays' r face fluxes, all faces at once (the columns only read them)
  DvceArray5D<Real> rf("cln_gd_rf", fv ? nmb1 + 1 : 1, 4, fv ? ke - ks + 1 : 1,
                       fv ? nj : 1, fv ? ie - is + 2 : 1);
  if (fv) {
    par_for("cln_gd_rf", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie + 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4];
      rflux(m, k, j, i, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
  }
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
      for (int c = 0; c < 4 && fv; ++c) q[c] = rf(m,c,k-ks,j-js,0);
      const Real vb = in ? 0.0 : qr_(m,4*nkj + kj);   // the column's volume below
      for (int i = is; i <= ie; ++i) {
        for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
        if (x.cub) {   // (no curvature terms: Cartesian momentum)
          if (!fv) rflux(m, k, j, i + 1, q);
          for (int c = 0; c < 4 && fv; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
          continue;
        }
        Real g[2] = {0.0, 0.0}, qi = q[2];
        // the exchange V chi (S - J, -H_b) and the lateral outflow of the cell
        const Real chi = x.vol(m,k,j,i)*exp(cw_(m,0,k+ok,j+og,i));
        const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
        for (int d = 0; d < n && !fv; ++d) {   // (fv: the face fluxes are the rays')
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
        for (int c = 0; c < 4 && fv; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
        // the angular-momentum weight r_f/r of the y_theta row on the r faces
        g[1] += q[2]*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) - qi*(1.0 - x.x1f(m,i)/x.x1v(m,i));
        for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
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
