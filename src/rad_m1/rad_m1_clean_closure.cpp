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

#include <Kokkos_Sort.hpp>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
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
    for (int c = 0; c < 3; ++c) a[c] = x.area[dir](m,k,j,i)*x.g(m,CLG_N+3*dir+c,k,j,i);
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
//! \fn Real CleanFvFace
//! \brief clean_rays = fv: the intensity of direction d on the low face of cell (k,j,i)
//! normal to dir for a beam along +e_dir (pos) or -e_dir: that of the upwind cell, with
//! o2 its minmod slope towards the face (second order in index space, positive), the
//! downwind cell's intensity read from fd; lld (nl directions per block) its trace,
//! with the curvature Q of its source from fd (cln_fvq)

template <typename FI, typename FD>
KOKKOS_INLINE_FUNCTION
Real CleanFvFace(const FI &fi, const FD &fd, const bool o2,
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

template <typename FI>
KOKKOS_INLINE_FUNCTION
Real CleanFvNbr(const FI &fi, const DvceArray4D<int> &mp, const bool ax,
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

template <typename FI, typename FD>
KOKKOS_INLINE_FUNCTION
Real CleanFvFaceAxi(const FI &fi, const FD &fd,
                    const DvceArray4D<int> &mp, const bool o2, const int m, const int d,
                    const int k, const int j, const int i, const int o, const bool pos,
                    const int nl = 0) {
  auto at = [&](const auto &f, const int s, const int ob = 0) {
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
//! \fn Real CleanThickI
//! \brief clean_ray_cut: the intensity of direction n in a cell below the cut, the
//! diffusion limit I = J + 3 n . H of the cell's E and F (iw, the state the rays are
//! built for); its moments are D = I/3 and h+- = +-1/4 + h/2 of the set's quadrature

KOKKOS_INLINE_FUNCTION
Real CleanThickI(const CleanArgs &x, const DvceArray5D<Real> &iw,
                 const DvceArray2D<Real> &x3v, const Real ic, const Real n[3],
                 const int m, const int k, const int j, const int i) {
  Real e[3][3];
  if (x.cub) {
    CleanCsFrame(x, m, k, j, i, e);
  } else {
    CleanLocalBasis(x.x2v(m,j), x3v(m,k), e);
  }
  const Real f[3] = {iw(m,M1_IW_F1,k,j,i), iw(m,M1_IW_F2,k,j,i), iw(m,M1_IW_F3,k,j,i)};
  Real nf = 0.0;
  for (int b = 0; b < 3; ++b) nf += (n[0]*e[b][0] + n[1]*e[b][1] + n[2]*e[b][2])*f[b];
  return fmax(iw(m,M1_IW_EN,k,j,i) + 3.0*ic*nf, 0.0);
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
  par_for("cln_edd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
  CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
  cln_havec = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvSourceAlloc
//! \brief the FV rays' source moments in the mesh layout (cln_fvq, CLQ_*)

void RadiationM1::CleanFvSourceAlloc() {
  if (cln_fvq.extent_int(0) == pmy_pack->nmb_thispack) return;
  Kokkos::realloc(cln_fvq, pmy_pack->nmb_thispack, (cln_rfv == 3) ? CLQ_N : CLQ_N0,
                  cln_s.extent_int(2), cln_s.extent_int(3), cln_s.extent_int(4));
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
//! ranges of the solve (CleanCellConsistency) is the set of the step's rotation block;
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

  par_for("cln_predict", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
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
  MemcutPoison(2);   // memcut-1113 diagnostic (env MEMCUT_POISON)
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
  } else if (cln_ring) {   // (the sc rays: each set's sweep in CleanRingMoments)
    VetGdSource();
    CleanFvSets(pmy_pack->pmesh->time + pmy_pack->pmesh->dt, pmy_pack->pmesh->ncycle);
  } else {
    VetGdBuild();
    if (cln_scmem) VetGdTopCopy();
  }
  if (cln_ring) {
    CleanRingMoments();
  } else {
    CleanClosureFromGd();
    // (clean_rays = sc on the cubed sphere: also the intensities cln_fvi of the face half
    // ranges, CleanCellConsistency)
    if (cln_consist || cln_cssc) CleanGdFaceSums();
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

//----------------------------------------------------------------------------------------
//! \fn bool CleanRingMix
//! \brief clean_rot_ring: the raw moments q of the union of the current set (slot 0 of
//! rl), the next one (slot qn, weight s) and the older sets r (weights 1, the last 1 - s)
//! from their normalised moments (rl, the offsets rr): J-weighted with the ratios r_r =
//! J_{b-r}/J_b, J their weighted mean.  Invalid offsets (r_r < 0) and a current set
//! without rays (false) are left out, as the production ring (RotAvgApply)

KOKKOS_INLINE_FUNCTION
bool CleanRingMix(const DvceArray5D<Real> &rl, const DvceArray5D<Real> &rr,
                  const Kokkos::Array<int, CLN_N> &fm, const int nf, const int qn,
                  const int nr, const Real s, const int m, const int k, const int j,
                  const int i, Real *q) {
  const int fj = fm[int{CLN_J}];
  const Real jb = rl(m,fj,k,j,i);
  if (!(jb > 0.0)) return false;
  Real a[RadiationM1::M1_RA_NF];
  const Real wn = s*rl(m,qn*nf+fj,k,j,i)/jb;
  Real na = 1.0 + s, nw = 1.0 + wn;
  for (int r = 0; r < nr; ++r) {
    const Real rj = rr(m,r*nf+fj,k,j,i), w = (r == nr - 1) ? 1.0 - s : 1.0;
    a[r] = (rj >= 0.0) ? w*rj : 0.0;
    na += (rj >= 0.0) ? w : 0.0;
    nw += a[r];
  }
  q[CLN_J] = jb*nw/na;
  for (int n = 0; n < CLN_N; ++n) {
    if (n == CLN_J) continue;
    const int f = fm[n];
    if (f < 0) {q[n] = 0.0; continue;}
    const Real gb = rl(m,f,k,j,i);
    Real v = gb + wn*rl(m,qn*nf+f,k,j,i);
    for (int r = 0; r < nr; ++r) v += a[r]*(gb + rr(m,r*nf+f,k,j,i));
    q[n] = v/nw*q[CLN_J];
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanRingFields
//! \brief clean_rot_ring: the fields the offsets and the Delta carry, fm(CLN index) =
//! their slot or -1, and their number.  The staggered layout leaves out those none of its
//! solves reads: those cln_s does not hold (CleanNs: CLN_KP, CLN_F1, off the cubed sphere
//! CLN_FH) and off the cubed sphere the stress parts of the face sums CLN_FX, CLN_F2,
//! CLN_F3 (only their energy flux is read; the sp F solve reads CLN_G).  The mixture
//! gives 0 there (CleanRingMix)

int RadiationM1::CleanRingFields(int *fm) const {
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  const int ns = CleanNs(cln_cell, cub);   // (the slots of cln_s)
  int nf = 0;
  for (int n = 0; n < CLN_N; ++n) {
    const bool fs = n >= CLN_F2 && n < CLN_G && (n - CLN_F2)%4 != 0;   // (stress parts)
    fm[n] = (cln_cell || (n < ns && !(fs && !cub))) ? nf++ : -1;
  }
  return nf;
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
  const int nmb = pmy_pack->nmb_thispack, nk = vra_k, nr = nk - 1, ns = cln_rns;
  const int64_t b = cln_rbld;
  const Real s = cln_rs;
  const int c3 = ke - ks + 1, c2 = je - js + 1, c1 = ie - is + 1;
  // (the swept sets' moments: the current set in cln_rl slot 0, the next in slot 1; the
  // older ones of a fill straight into their offsets' slots of cln_rr.  All on the
  // interior cells; the offsets and the Delta only of the fields fm carries)
  Kokkos::Array<int, CLN_N> fm;
  const int nf = CleanRingFields(fm.data());
  if (cln_rl.extent_int(0) != nmb || cln_rl.extent_int(1) != 2*nf) {
    Kokkos::realloc(cln_rl, nmb, 2*nf, c3, c2, c1);
  }
  if (cln_rr.extent_int(0) != nmb) {   // (a restart has read them)
    vra_gen += 1;   // (cln_rr's writers: the dt_retry checkpoint skips a clean copy)
    Kokkos::realloc(cln_rr, nmb, nr*nf, c3, c2, c1);
    Kokkos::realloc(cln_rpd, nmb, nf, c3, c2, c1);
    cln_rok = false;
  }
  std::vector<int> slot(nk + 2, -1);   // slot of block b - nk + 1 + p
  for (int q = 0; q < ns; ++q) {
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    if (p >= 0 && p <= nk) slot[p] = q;
  }
  bool full = true;
  for (int p = 0; p < nr; ++p) full = full && slot[p] >= 0;
  const bool shift = cln_rok && b == cln_rb + 1 && cln_rpb == cln_rb;
  const bool fill = (!cln_rok || b != cln_rb) && !shift && full;
  // the moments of each swept set (the ray_min test on the mixture below)
  auto dfull = cln_fvd;
  const Real ray_min = cln_raymin;
  cln_raymin = 0.0;
  auto cs = cln_s;
  auto rl = cln_rl;
  auto rr = cln_rr;
  // (the sc rays on the cubed sphere: the sets in block order, so that the last two
  // swept, b and b+1, hold the two cln_fvi slots blk mod 2; on sp the slot order: the
  // block order breaks the bitwise restart of q_def, reflecting top)
  std::vector<int> ord(ns);
  for (int q = 0; q < ns; ++q) ord[q] = q;
  std::sort(ord.begin(), ord.end(),
            [&](int u, int v) {return cln_fvblk[u] < cln_fvblk[v];});
  for (int qo = 0; qo < ns; ++qo) {
    const int q = (cln_rfv == 0 && cln_cssc) ? ord[qo] : qo;
    const int64_t p = cln_fvblk[q] - (b - nk + 1);
    const bool older = (p >= 0 && p < nr);
    const bool keep = (older && fill) || p == nr || p == nk;
    if (older && fill) vra_gen += 1;
    cln_fvd = cln_fvdm[q];
    if (cln_rfv == 0) {   // the sc rays: set q swept alone (its tables, its own weights)
      cln_rq = q;
      cln_rqm = static_cast<int>(cln_fvblk[q] & 1);
      VgdXfTables(cln_fvkey[q]);
      vgd_scut = vlat_icut;
      VetGdSweep();
      if (cln_scmem) VetGdTopCopy();
    }
    CleanClosureFromGd();
    if (cln_consist || cln_cssc) CleanGdFaceSums();
    if (!keep) continue;
    const int o = static_cast<int>(older ? (nr - 1 - p)*nf : ((p == nk) ? nf : 0));
    auto ra = older ? rr : rl;
    par_for("cln_ring_raw", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jr = cs(m,CLN_J,k,j,i), ij = (jr > 0.0) ? 1.0/jr : 0.0;
      for (int n = 0; n < CLN_N; ++n) {
        if (fm[n] < 0) continue;
        ra(m,o+fm[n],k-ks,j-js,i-is) = (CleanCellMoment(n) ? 1.0 : ij)*cs(m,n,k,j,i);
      }
    });
  }
  cln_fvd = dfull;
  if (cln_cssc) CleanScUnion2();
  cln_raymin = ray_min;
  // (jump: the next set only at the build that closes the block; its weight is s = 0)
  const int qb = 0, qn = (slot[nk] >= 0) ? 1 : 0;
  auto pd = cln_rpd;
  if (!cln_rok || b != cln_rb) {
    if (shift) {
      // the pending Delta is that of the boundary of block b when it closed the block
      const Real tk = vgd_rot_t0 + static_cast<Real>(b)*vgd_rot_dt;
      if (std::abs(cln_rpt - tk) > 1.0e-9*vgd_rot_dt) cln_roff_sum += 1.0;
      vra_gen += 1;
      par_for("cln_ring_shift", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is, fj = fm[int{CLN_J}];
        const Real rp = pd(m,fj,k,j,i);   // (J ratios: r_j' = r_{j-1}/R, r_1' = 1/R)
        for (int f = 0; f < nf; ++f) {
          for (int r = nr - 1; r > 0; --r) {
            const Real o = rr(m,(r-1)*nf+f,k,j,i);
            rr(m,r*nf+f,k,j,i) = (f != fj) ? o - pd(m,f,k,j,i) :
                                 ((rp > 0.0 && o >= 0.0) ? o/rp : -1.0);
          }
          rr(m,f,k,j,i) = (f != fj) ? -pd(m,f,k,j,i) : ((rp > 0.0) ? 1.0/rp : -1.0);
        }
      });
    } else if (full) {
      // the whole union swept: the offsets of the older sets (their moments, above) from
      // the current one
      par_for("cln_ring_fill", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k0, const int j0, const int i0) {
        const int k = k0 - ks, j = j0 - js, i = i0 - is;
        const Real jb = rl(m,qb*nf+fm[int{CLN_J}],k,j,i);
        for (int r = 0; r < nr; ++r) {
          for (int n = 0; n < CLN_N; ++n) {
            const int f = fm[n];
            if (f < 0) continue;
            const Real o = rr(m,r*nf+f,k,j,i);
            rr(m,r*nf+f,k,j,i) = (n != CLN_J) ? o - rl(m,qb*nf+f,k,j,i) :
                                 ((jb > 0.0) ? o/jb : -1.0);
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
      const int fj = fm[int{CLN_J}];
      const Real jb = rl(m,qb*nf+fj,k,j,i), jn = rl(m,qn*nf+fj,k,j,i);
      for (int n = 0; n < CLN_N; ++n) {
        const int f = fm[n];
        if (f < 0) continue;
        pd(m,f,k,j,i) = (n != CLN_J) ? rl(m,qn*nf+f,k,j,i) - rl(m,qb*nf+f,k,j,i) :
                            ((jb > 0.0) ? jn/jb : -1.0);
      }
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
      const Real jo = rl(m,o*nf+fm[int{CLN_J}],kr,jr,ir);
      return (n == CLN_J) ? jo : ((fm[n] < 0) ? 0.0 : rl(m,o*nf+fm[n],kr,jr,ir)*jo);
    };
    Real q[CLN_N];
    const bool mix = ring && CleanRingMix(rl, rr, fm, nf, qn, nr, s, m, kr, jr, ir, q);
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
    if (!(CleanRingMix(rl, rr, fm, nf, qn, nr, 1.0, m, kr, jr, ir, q) &&
          q[CLN_J] > 0.0 && CleanRealizable(x, q, m, k, j, i))) {
      if (rl(m,qn*nf+fm[int{CLN_J}],kr,jr,ir) > 0.0) ++ng;
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
  CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto mp = cln_fvm;
  const int n1 = fi.extent_int(4), n2 = fi.extent_int(3), ng = indcs.ng;
  // the ghost cells of the theta seams (rows j < js, j > je of ks <= k <= ke) and of the
  // phi seams (layers k < ks, k > ke), i fastest, in groups that permute their
  // intensities through a scratch array of at most 2^24 values
  const int64_t nt = 2*ng*(ke - ks + 1)*n1, np = circ ? 0 : 2*ng*n2*n1;
  const int64_t ne = pmy_pack->nmb_thispack*(nt + np);
  const int64_t nq = std::max(int64_t(1), std::min(ne, int64_t(1 << 24)/n));
  if (cln_fvt.extent_int(1) != n) Kokkos::realloc(cln_fvt, nq, n);
  auto t = cln_fvt;
  for (int64_t e0 = 0; e0 < ne; e0 += nq) {
    Kokkos::parallel_for("cln_fv_seam", Kokkos::RangePolicy<>(DevExeSpace(), e0,
                                                             std::min(ne, e0 + nq)),
    KOKKOS_LAMBDA(const int64_t e) {
      const int m = static_cast<int>(e/(nt + np)), i = static_cast<int>(e%n1);
      const int64_t r = e%(nt + np), q = (r < nt) ? r/n1 : (r - nt)/n1;
      const bool th = (r < nt);
      const int g = static_cast<int>(th ? q/(ke - ks + 1) : q/n2);
      const int gs = (g < ng) ? g + (th ? js : ks) - ng : g - ng + (th ? je : ke) + 1;
      const int k = th ? ks + static_cast<int>(q%(ke - ks + 1)) : gs;
      const int j = th ? gs : static_cast<int>(q%n2);
      const int sd = (th ? 0 : 2) + ((g < ng) ? 0 : 1);
      if (x.mbc(m,2+sd) != BoundaryFlag::periodic) return;
      for (int b = 0; b < nb; ++b) {
        for (int d = 0; d < n; ++d) t(e-e0,d) = fi(m,b*n+mp(m,sd,k,d),k,j,i);
        for (int d = 0; d < n; ++d) fi.Set(m, b*n+d, k, j, i, t(e-e0,d));
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
  // clean_fv_lag: the last step swept the same set unless it belongs to another rotation
  // block (a function of tc and dt alone, so that a restart decides as the run would)
  cln_fvfull = (cln_ring != 0);
  cln_fvtc = tc;
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
  int64_t b = VetGdBlock(cyc);
  if (vgd_rot_dt <= 0.0) cln_fvfull = cln_fvfull || (cyc < 1) || VetGdBlock(cyc - 1) != b;
  Real s = 0.0;
  if (vgd_rot_dt > 0.0) {
    // (without the cross-fade the block of the step's midpoint, as VetGdBlock)
    const Real dt = (pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    const Real tau = (tc - (xf ? 0.0 : 0.5*dt) - vgd_rot_t0)/vgd_rot_dt;
    b = static_cast<int64_t>(std::floor(tau));
    if (xf) s = tau - static_cast<Real>(b);
    cln_fvfull = cln_fvfull || std::floor(tau - pm->dt/vgd_rot_dt) != std::floor(tau);
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
  const int on = cln_fvn/n;
  const std::vector<int64_t> oblk = cln_fvblk;
  const CleanFvI ofv = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  const bool remap = cln_ring && cln_rfv > 0 && on > 0 && cln_fvn != ns*n &&
                     ofv.extent_int(1) == cln_fvb*on*n;
  // (the sets' tables by their own count: CleanScUnion2 resets cln_fvn to two sets)
  if (cln_fvn != ns*n || static_cast<int>(cln_fvkey.size()) != ns) {
    cln_fvn = ns*n;
    cln_fvd = DvceArray2D<Real>("cln_fvd", ns*n, 4);
    cln_fvkey.assign(ns, std::numeric_limits<Real>::quiet_NaN());
    cln_fvblk.assign(ns, 0);
    if (cln_ring) {
      // (the sweep's pass flags are sized for the swept directions)
      if (cln_fvdue.extent_int(0) < pmy_pack->nmb_thispack*ns*n) {
        cln_fvdue = DvceArray2D<int8_t>();
      }
      cln_fvdm.resize(ns);
      for (int q = 0; q < ns; ++q) cln_fvdm[q] = DvceArray2D<Real>("cln_fvdm", ns*n, 4);
    }
  }
  auto d_h = Kokkos::create_mirror_view(cln_fvd);
  std::vector<bool> chg(ns, false);
  for (int64_t blk = (two ? b : b - nk + 1); blk <= b + nx; ++blk) {
    const int q = static_cast<int>(((blk % ns) + ns) % ns);
    // (Cartesian: vet_sc's set of the block, turned about x1 by beta, VetRotBeta)
    const Real key = vgd_on ? VetGdAngleBlk(blk) : VetRotBeta(blk);
    Real w = ((xf && blk == b - nk + 1) ? 1.0 - s : ((blk > b) ? s : 1.0))/nk;
    if (two) w = (blk > b) ? s : 1.0 - s;
    std::vector<double> d;
    if (vgd_on) {
      d = VgdRotate(key);
    } else {
      VetBuildAngles(key);
      auto a_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vet_ang);
      d.assign(a_h.data(), a_h.data() + 4*n);
    }
    for (int e = 0; e < n; ++e) {
      for (int c = 0; c < 3; ++c) d_h(q*n+e,c) = d[4*e+c];
      d_h(q*n+e,3) = w*d[4*e+3];
    }
    chg[q] = !(key == cln_fvkey[q]);
    cln_fvkey[q] = key;
    cln_fvblk[q] = blk;
  }
  Kokkos::deep_copy(cln_fvd, d_h);
  if (!vgd_on) VetBuildAngles(vsc_beta);   // (vet_sc's own set back)
  if (remap) {
    // clean_rot_ring: a set kept keeps its intensities, a new one starts from the newest
    // old set's (the sweep's starting point)
    auto ofi = cln_fvi;
    auto ofs = cln_fvis;
    const int e0 = ofv.extent_int(0), e2 = ofv.extent_int(2), e3 = ofv.extent_int(3);
    const int e4 = ofv.extent_int(4);
    DvceArray5D<Real> fi;
    DvceArray5D<float> fs;
    if (ofv.f32) {
      fs = DvceArray5D<float>("cln_fvis", e0, cln_fvb*ns*n, e2, e3, e4);
    } else {
      fi = DvceArray5D<Real>("cln_fvi", e0, cln_fvb*ns*n, e2, e3, e4);
    }
    auto AL = Kokkos::ALL;
    for (int q = 0; q < ns; ++q) {
      int oq = 0;
      for (int r = 0; r < on; ++r) {
        if (oblk[r] == cln_fvblk[q] || (oblk[oq] != cln_fvblk[q] && oblk[r] > oblk[oq])) {
          oq = r;
        }
      }
      for (int c = 0; c < cln_fvb; ++c) {
        const int a = (c*ns + q)*n, o = (c*on + oq)*n;
        auto ra = std::make_pair(a, a + n), ro = std::make_pair(o, o + n);
        if (ofv.f32) {
          Kokkos::deep_copy(Kokkos::subview(fs, AL, ra, AL, AL, AL),
                            Kokkos::subview(ofs, AL, ro, AL, AL, AL));
        } else {
          Kokkos::deep_copy(Kokkos::subview(fi, AL, ra, AL, AL, AL),
                            Kokkos::subview(ofi, AL, ro, AL, AL, AL));
        }
      }
    }
    cln_fvi = fi;
    cln_fvis = fs;
  }
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
  for (int q = 0; q < ns && cln_rfv > 0; ++q) {   // (the sc rays: vet_gd's own tables)
    if (chg[q]) CleanFvTables(q*n, n);
  }
}

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
//! \fn void RadiationM1::CleanFvTables
//! \brief the tables of the swept directions d0 .. d0+nd-1 (one set; every search stays
//! in it): the direction maps of the periodic seams of a wedge (cln_fvm: the ghost's
//! direction in the frame of its image cell, n' = sum_b (n . e_b(ghost)) e_b(image), the
//! nearest set direction; an exact permutation for a phi period of the set) and of the
//! reflecting theta walls (cln_fvwm: an incoming direction takes the leaving one nearest
//! to its mirror image n - 2 (n . t) t about the wall's unit normal t, the mean
//! theta_hat over the face's phi range; cln_fvwu = |n . t|) and of the reflecting x1
//! faces (cln_fvrm, cln_fvru: t along the face's vector area); their sweep levels are due

void RadiationM1::CleanFvTables(const int d0, const int nd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb = pmy_pack->nmb_thispack, nn = cln_fvn;
  const int c3 = cln_s.extent_int(2), c2 = cln_s.extent_int(3);
  if (cln_fvm.extent_int(3) != nn) {
    Kokkos::realloc(cln_fvm, nmb, 4, c3, nn);
    cln_fvdo.assign(nn, true);
  }
  for (int d = d0; d < d0 + nd; ++d) cln_fvdo[d] = true;
  CleanArgs x = CleanMakeArgs();
  auto dir_ = cln_fvd;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
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
  const bool rw = (ibc_x1min == M1_IBC_REFLECT) || (ibc_x1max == M1_IBC_REFLECT);
  if (!CleanFvWall(0) && !CleanFvWall(1) && !rw) return;
  if (cln_fvwm.extent_int(3) != nn) {
    Kokkos::realloc(cln_fvwm, nmb, 2, c3, nn);
    Kokkos::realloc(cln_fvwu, nmb, 2, c3, nn);
    Kokkos::realloc(cln_fvrm, nmb, 2, c3, c2, nn);
    Kokkos::realloc(cln_fvru, nmb, 2, c3, c2, nn);
    Kokkos::realloc(cln_fvw, nmb, 4, c3, std::max(cln_s.extent_int(4), c2), cln_fvb*nn);
  }
  auto x2f = pmy_pack->pcoord->xx2f;
  auto x3f = pmy_pack->pcoord->xx3f;
  auto wm = cln_fvwm;
  auto wu = cln_fvwu;
  auto rm = cln_fvrm;
  auto ru = cln_fvru;
  if (rw) {
    par_for("cln_fv_x1map", DevExeSpace(), 0, nmb - 1, 0, 1, ks, ke, js, je, d0,
            d0 + nd - 1,
    KOKKOS_LAMBDA(const int m, const int sd, const int k, const int j, const int d) {
      Real t[3];
      CleanFaceVector(x, x3f, 0, m, k, j, sd ? ie + 1 : is, t);
      const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
      for (int c = 0; c < 3; ++c) t[c] /= tn;
      rm(m,sd,k,j,d) = CleanFvMirror(dir_, d0, nd, d, t, sd ? 1.0 : -1.0, ru(m,sd,k,j,d));
    });
  }
  if (!CleanFvWall(0) && !CleanFvWall(1)) return;
  par_for("cln_fv_wallmap", DevExeSpace(), 0, nmb - 1, 0, 1, ks, ke, d0, d0 + nd - 1,
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int d) {
    const Real th = x2f(m, sd ? je + 1 : js), sgn = sd ? 1.0 : -1.0;
    const Real p1 = x3f(m,k), p2 = x3f(m,k+1);
    Real t[3] = {cos(th)*(sin(p2) - sin(p1)), cos(th)*(cos(p1) - cos(p2)),
                 -sin(th)*(p2 - p1)};
    const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    for (int c = 0; c < 3; ++c) t[c] /= tn;
    wm(m,sd,k,d) = CleanFvMirror(dir_, d0, nd, d, t, sgn, wu(m,sd,k,d));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvWalls
//! \brief the incoming intensities of the reflecting theta walls and x1 faces (cln_fvw):
//! each incoming direction e takes the intensity of the wall cell along its mirror
//! direction, scaled on each wall face so that it returns exactly the energy it
//! receives, sum_e w_e |n_e . t| I(e) = sum_out w_d |n_d . t| I_d (isotropic radiation
//! stays isotropic: the set integrates n exactly); lld: the wall face traces, with the
//! slopes along the face (blocks 1, 3 on a theta wall, 2, 3 on an x1 face).  A
//! reflecting x1max (top) only mirrors the intensities a sweep starts from (the last
//! build, kept by the restart files: a one-build lag, as vet_gd's VetGdTopMirror), so
//! that the passes of a closed thin cavity converge

void RadiationM1::CleanFvWalls(const bool top) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int n = cln_fvn, nl = (cln_fvb > 1) ? n : 0;
  // (theta low, high, x1min, x1max)
  const bool on[4] = {CleanFvWall(0), CleanFvWall(1), ibc_x1min == M1_IBC_REFLECT,
                      top && ibc_x1max == M1_IBC_REFLECT};
  CleanArgs x = CleanMakeArgs();
  CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto fq = cln_fvq;
  auto fw = cln_fvw;
  auto mp = cln_fvwm;
  auto u = cln_fvwu;
  auto rm = cln_fvrm;
  auto ru = cln_fvru;
  auto dir_ = cln_fvd;
  par_for("cln_fv_wall", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, 0, 3, indcs.ks,
          indcs.ke, 0, std::max(ie, je),
  KOKKOS_LAMBDA(const int m, const int sd, const int k, const int c) {
    // (c: the wall cell's i on a theta wall, its j on an x1 face)
    const int s = sd%2, r = sd/2, b0 = 2 - r;   // (b0: the slope normal to the wall)
    if (!on[sd] || c < (r ? js : is) || c > (r ? je : ie)) return;
    const int jc = r ? c : (s ? je : js), ic = r ? (s ? ie : is) : c;
    if (!CleanBoundaryFace(x, 1 - r, m, k, jc + (!r && s), ic + (r && s))) return;
    Real so = 0.0, si = 0.0;
    for (int e = 0; e < n; ++e) {
      const int d = r ? rm(m,s,k,jc,e) : mp(m,s,k,e), g = (d >= 0) ? d : e;
      const Real tr = fi(m,g,k,jc,ic) +
                      ((nl > 0) ? (s ? 1.0 : -1.0)*fi(m,b0*n+g,k,jc,ic) +
                                  fq(m,CLQ_C+b0-1,k,jc,ic) : 0.0);
      const Real q = dir_(e,3)*(r ? ru(m,s,k,jc,e) : u(m,s,k,e))*tr;
      fw(m,sd,k,c,e) = (d >= 0) ? tr : 0.0;
      for (int b = 1; b < 4 && nl > 0; ++b) {
        if (b != b0) fw(m,sd,k,c,b*n+e) = (d >= 0) ? fi(m,b*n+d,k,jc,ic) : 0.0;
      }
      if (d >= 0) {si += q;} else {so += q;}
    }
    for (int e = 0; e < n; ++e) {
      for (int b = 0; b < 4; ++b) {
        if (b == 0 || (nl > 0 && b != b0)) fw(m,sd,k,c,b*n+e) *= (si > 0.0) ? so/si : 0.0;
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvLevels
//! \brief the wavefronts of the sweep for the directions that are due: in each
//! direction's cell order (n . x ascending; centroids on the cubed sphere) a cell takes
//! the level after every face neighbour placed before it across a face the direction
//! crosses (fv: every cell up to two along an axis), so that a level's cells of all
//! directions read nothing written in the level and a pass gives the sequential
//! sweep's intensities.  cln_fvc lists the cells of each m n + d by (level, cell),
//! cln_fvls their level starts; level l of the sweep takes the m n + d in turn
//! (cln_fvg, cln_fvo), cln_fvlv[l] cells.  One team per direction takes the order in
//! chunks; the levels of a chunk are relaxed from below until they no longer change
//! (the longest path).  The directions go in groups of at most 2^27 cells (scratch).

void RadiationM1::CleanFvLevels() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3, nc = nx1*nx2*nx3;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nmb = pmy_pack->nmb_thispack, n = cln_fvn, r2 = (cln_rfv == 2) ? 2 : 1;
  const int nmd = nmb*n, nb = 1024;
  const bool axi = !pmy_pack->pmesh->three_d && pmy_pack->pmesh->use_spherical_polar;
  CleanArgs ca = CleanMakeArgs();
  auto fa = cln_fva;
  auto dir_ = cln_fvd;
  auto x1v = pmy_pack->pcoord->x1v;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto gn = cln_gn;
  // the factors of n . x of every cell (r, sin th, cos ph, sin ph, cos th; the centroid)
  DvceArray3D<Real> xc("cln_fvxc", nmb, 5, nc);
  par_for("cln_fv_xc", DevExeSpace(), 0, nmb - 1, 0, nc - 1,
  KOKKOS_LAMBDA(const int m, const int c) {
    const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
    const Real f[5] = {x1v(m,i), sin(x2v(m,j)), cos(x3v(m,k)), sin(x3v(m,k)),
                       cos(x2v(m,j))};
    // (the cubed sphere: the centroid; Cartesian: the centre)
    const Real xp[3] = {x1v(m,i), x2v(m,j), x3v(m,k)};
    for (int a = 0; a < 5; ++a) {
      xc(m,a,c) = (a >= 3 || ca.sph) ? f[a] : (ca.cub ? gn(m,CLG_X+a,k,j,i) : xp[a]);
    }
  });
  if (cln_fvc.extent(0) != static_cast<size_t>(nmd)*nc) {
    Kokkos::realloc(cln_fvc, static_cast<size_t>(nmd)*nc);
    cln_fvls.assign(nmd, std::vector<int>());
  }
  // (clean_fv_joint: the MeshBlocks of a direction in one group)
  cln_fvjoint = cln_fvjoint && cln_rfv == 1 && !ca.sph && !pmy_pack->pmesh->multilevel &&
                !VgdRotOn();
  const bool jt = cln_fvjoint;
  std::vector<int> due;
  for (int q = 0; q < nmd; ++q) {
    const int md = jt ? (q%nmb)*n + q/nmb : q;
    if (cln_fvdo[md%n]) due.push_back(md);
  }
  int nw = std::max(1, std::min(static_cast<int>(due.size()), (1 << 27)/nc));
  if (jt) nw = std::max(nmb, nw - nw%nmb);
  cln_fvl0 = std::numeric_limits<int>::max();
  cln_fvjn[0] = cln_fvjn[1] = 0;
  if (jt) {
    if (cln_fvcore.extent_int(0) != nmd || cln_fvcore.extent_int(1) != nc) {
      Kokkos::realloc(cln_fvcore, nmd, nc);
    }
    Kokkos::deep_copy(cln_fvcore, 0);
    cln_fvseg.clear();
  }
  DvceArray1D<int> mds("cln_fvl_md", nw);
  DvceArray2D<int> ord("cln_fvord", nw, nc), pos("cln_fvpos", nw, nc);
  DvceArray2D<int> lev("cln_fvlev", nw, nc);
  DvceArray1D<Real> key("cln_fvkey", nc);
  DvceArray1D<int64_t> lk("cln_fvl_key", static_cast<size_t>(nw)*nc);
  auto mds_h = Kokkos::create_mirror_view(mds);
  auto fc = cln_fvc;
  const bool cub = !ca.sph;
  for (size_t w0 = 0; w0 < due.size(); w0 += nw) {
    const int nq = std::min(nw, static_cast<int>(due.size() - w0));
    for (int w = 0; w < nq; ++w) mds_h(w) = due[w0 + w];
    Kokkos::deep_copy(mds, mds_h);
    int lmax = 0;
    if (jt) {
      lmax = CleanFvJoint(mds, nq, lev, xc);
    } else {
      // the cell order of each direction: n . x ascending, equal keys by the cell index
      for (int w = 0; w < nq; ++w) {
        const int m = mds_h(w)/n, d = mds_h(w)%n;
        auto o = Kokkos::subview(ord, w, Kokkos::ALL);
        par_for("cln_fv_key", DevExeSpace(), 0, nc - 1, KOKKOS_LAMBDA(const int c) {
          key(c) = CleanFvKey(xc, dir_, cub, m, d, c);
          o(c) = c;
        });
        Kokkos::sort(DevExeSpace(), o, KOKKOS_LAMBDA(const int a, const int b) {
          return key(a) < key(b) || (key(a) == key(b) && a < b);
        });
      }
      par_for("cln_fvl_pos", DevExeSpace(), 0, nw - 1, 0, nc - 1,
      KOKKOS_LAMBDA(const int w, const int p) {
        if (w < nq) pos(w,ord(w,p)) = p;
        lev(w,p) = 0;
      });
      Kokkos::parallel_reduce("cln_fvl_lev", Kokkos::TeamPolicy<>(DevExeSpace(), nq,
                                                                  Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamMember_t &tm, int &lm) {
        const int w = tm.league_rank(), m = mds(w)/n, d = mds(w)%n;
        const int nx[3] = {nx1, nx2, nx3};
        for (int p0 = 0; p0 < nc; p0 += nb) {
          const int p1 = (p0 + nb < nc) ? p0 + nb : nc;
          int chg = 1;
          while (chg > 0) {
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, p0, p1),
            [&](const int p, int &ch) {
              const int c = ord(w,p), x[3] = {c%nx1, (c/nx1)%nx2, c/(nx1*nx2)};
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
                  if (pos(w,g) < p) l = (l > lev(w,g) + 1) ? l : lev(w,g) + 1;
                }
              }
              if (l != lev(w,c)) {
                lev(w,c) = l;
                ++ch;
              }
            }, chg);
          }
        }
        int l1 = 0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nc),
        [&](const int c, int &l) {
          l = (l > lev(w,c) + 1) ? l : lev(w,c) + 1;
        }, Kokkos::Max<int>(l1));
        Kokkos::single(Kokkos::PerTeam(tm), [&]() {lm = (lm > l1) ? lm : l1;});
      }, Kokkos::Max<int>(lmax));
    }
    // the cells of each direction sorted by (level, cell), and the count of each level
    DvceArray2D<int> cnt("cln_fvl_n", nq, lmax);
    par_for("cln_fvl_key", DevExeSpace(), 0, nw - 1, 0, nc - 1,
    KOKKOS_LAMBDA(const int w, const int c) {
      const int64_t wl = w*static_cast<int64_t>(lmax) + lev(w,c);
      lk(w*static_cast<int64_t>(nc) + c) = wl*nc + c;
      if (w < nq) Kokkos::atomic_add(&cnt(w,lev(w,c)), 1);
    });
    Kokkos::sort(DevExeSpace(), lk);
    par_for("cln_fvl_list", DevExeSpace(), 0, nq - 1, 0, nc - 1,
    KOKKOS_LAMBDA(const int w, const int p) {
      fc(mds(w)*static_cast<int64_t>(nc) + p) =
          static_cast<int>(lk(w*static_cast<int64_t>(nc) + p)%nc);
    });
    auto cnt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cnt);
    for (int w = 0; w < nq; ++w) {
      auto &ls = cln_fvls[mds_h(w)];
      ls.assign(lmax + 1, 0);
      for (int l = 0; l < lmax; ++l) ls[l+1] = ls[l] + cnt_h(w,l);
    }
  }
  if (jt) {
    const Real nd = CleanReduce(cln_fvjn[0], false), ne = CleanReduce(cln_fvjn[1], false);
    cln_fvexact = (ne == 0.0);
    if (global_variable::my_rank == 0 && !cln_fvjoint_said) {
      std::printf("<rad_m1> clean_fv_joint: %.0f dropped inflows (rank 0: %d cycle "
                  "segments), %.0f cells reading another rank%s\n", nd,
                  static_cast<int>(cln_fvseg.size()), ne,
                  cln_fvexact ? "" : ": the MeshBlocks' own levels instead");
      cln_fvjoint_said = true;
    }
    // (a rank's sweep waiting for another's would repeat its longer levels: off)
    if (!cln_fvexact) {
      cln_fvjoint = false;
      cln_fvdo.assign(n, true);
      CleanFvLevels();
      return;
    }
  }
  // the levels of the sweep: each takes the m n + d in turn
  int nlev = 0;
  for (int md = 0; md < nmd; ++md) {
    nlev = std::max(nlev, static_cast<int>(cln_fvls[md].size()) - 1);
  }
  Kokkos::realloc(cln_fvg, nlev, nmd + 1);
  Kokkos::realloc(cln_fvo, nlev, nmd);
  auto g_h = Kokkos::create_mirror_view(cln_fvg);
  auto o_h = Kokkos::create_mirror_view(cln_fvo);
  cln_fvlv.assign(nlev, 0);
  for (int l = 0; l < nlev; ++l) {
    int s = 0;
    for (int md = 0; md < nmd; ++md) {
      const auto &ls = cln_fvls[md];
      const int nl = static_cast<int>(ls.size()) - 1;
      g_h(l,md) = s;
      o_h(l,md) = md*static_cast<int64_t>(nc) + ((l < nl) ? ls[l] : 0) - s;
      s += (l < nl) ? ls[l+1] - ls[l] : 0;
    }
    g_h(l,nmd) = s;
    cln_fvlv[l] = s;
  }
  // the m n + d of every 64th entry of a level (where the sweep's search starts)
  cln_fvkl.assign(nlev + 1, 0);
  for (int l = 0; l < nlev; ++l) cln_fvkl[l+1] = cln_fvkl[l] + (cln_fvlv[l] + 63)/64;
  Kokkos::realloc(cln_fvk, cln_fvkl[nlev]);
  auto k_h = Kokkos::create_mirror_view(cln_fvk);
  for (int l = 0; l < nlev; ++l) {
    for (int md = 0; md < nmd; ++md) {
      for (int b = (g_h(l,md) + 63)/64; 64*b < g_h(l,md+1); ++b) {
        k_h(cln_fvkl[l] + b) = md;
      }
    }
  }
  Kokkos::deep_copy(cln_fvg, g_h);
  Kokkos::deep_copy(cln_fvo, o_h);
  Kokkos::deep_copy(cln_fvk, k_h);
  cln_fvdo.assign(n, false);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanFvJoint
//! \brief clean_fv_joint: the sweep levels of the directions mds(w), w < nq (each
//! direction's MeshBlocks consecutive, w = (direction) nmb + m), over all MeshBlocks of
//! the rank: a cell takes the level after its upwind face neighbours, across a MeshBlock
//! face the cell the ghost exchange copies (cln_fvsrc).  The levels are the rounds of a
//! topological order; where it stalls the remaining cells on and between cycles of the
//! upwind graph (those with a remaining downstream cell) drop their inflow from cells
//! later in n . x, and such a cell's upwind neighbour waits for it (the sweep reads the
//! last pass there).  Returns the number of levels; cln_fvl0 is the first level with a
//! dropped inflow or one from another rank, cln_fvexact whether there is none.

int RadiationM1::CleanFvJoint(const DvceArray1D<int> &mds, const int nq,
                              DvceArray2D<int> lev, const DvceArray3D<Real> &xc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nc = nx1*nx2*indcs.nx3;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, n = cln_fvn, nw = lev.extent_int(0);
  const int n1 = cln_s.extent_int(4), n2 = cln_s.extent_int(3), n3 = cln_s.extent_int(2);
  CleanArgs ca = CleanMakeArgs();
  if (cln_fvsrc.extent_int(0) != nmb) {
    // the source of every ghost: the exchange of the cells' global indices
    DvceArray5D<Real> id("cln_fvid", nmb, 1, n3, n2, n1);
    const int gs = pmy_pack->gids;
    par_for("cln_fvj_id", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool in = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      const int c = (i - is) + nx1*((j - js) + nx2*(k - ks));
      id(m,0,k,j,i) = in ? static_cast<Real>(gs + m)*nc + c : -1.0;
    });
    CleanGhostFill(id, 0, 1, false);
    Kokkos::realloc(cln_fvsrc, nmb, n3, n2, n1);
    auto src = cln_fvsrc;
    par_for("cln_fvj_src", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool in = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      const Real v = id(m,0,k,j,i);
      const int64_t g = static_cast<int64_t>(v), b = g/nc - gs;
      src(m,k,j,i) = (in || v < 0.0 || b < 0 || b >= nmb) ? -1 :
                     static_cast<int>(b*nc + g%nc);
    });
  }
  auto src = cln_fvsrc;
  auto fa = cln_fva;
  auto dir_ = cln_fvd;
  const bool cub = !ca.sph;
  // the face neighbour of cell c of w across face f: w' nc + c', -1 none (a physical
  // boundary, an axis of one cell), -2 a cell of another rank; up: an upwind one
  auto nbr = [=] KOKKOS_FUNCTION (const int w, const int c, const int f) {
    const int m = mds(w)/n, e = f/2, s = (f%2) ? 1 : -1;
    int x[3] = {c%nx1, (c/nx1)%nx2, c/(nx1*nx2)};
    const int nx[3] = {nx1, nx2, nc/(nx1*nx2)};
    if (!ca.cub && (e == 1 ? !ca.multi_d : e == 2 && !ca.three_d)) return -1;
    if (CleanBoundaryFace(ca, e, m, ks + x[2] + (f == 5), js + x[1] + (f == 3),
                          is + x[0] + (f == 1))) return -1;
    x[e] += s;
    if (x[e] >= 0 && x[e] < nx[e]) return w*nc + x[0] + nx1*(x[1] + nx2*x[2]);
    const int v = src(m,ks+x[2],js+x[1],is+x[0]);
    return (v < 0) ? -2 : (w - m + v/nc)*nc + v%nc;
  };
  auto up = [=] KOKKOS_FUNCTION (const int w, const int c, const int f) {
    const int m = mds(w)/n, d = mds(w)%n;
    const int k = ks + c/(nx1*nx2), j = js + (c/nx1)%nx2, i = is + c%nx1;
    const Real a = dir_(d,0)*fa(m,3*f,k,j,i) + dir_(d,1)*fa(m,3*f+1,k,j,i) +
                   dir_(d,2)*fa(m,3*f+2,k,j,i);
    return (a < 0.0) ? nbr(w, c, f) : -1;
  };
  auto key = [=] KOKKOS_FUNCTION (const int v) {
    return CleanFvKey(xc, dir_, cub, mds(v/nc)/n, mds(v/nc)%n, v%nc);
  };
  DvceArray2D<int8_t> drop("cln_fvj_drop", nw, nc), al("cln_fvj_al", nw, nc);
  DvceArray2D<int8_t> hd("cln_fvj_hd", nw, nc);
  par_for("cln_fvj_init", DevExeSpace(), 0, nw - 1, 0, nc - 1,
  KOKKOS_LAMBDA(const int w, const int c) {
    lev(w,c) = (w < nq) ? -1 : 0;
  });
  // the rounds: a cell is ready when its inflows and the cells it dropped an inflow
  // from wait for (dropped inflows of its neighbours from it) are of earlier rounds
  const auto pol = Kokkos::RangePolicy<>(DevExeSpace(), 0, nq*nc);
  auto core = cln_fvcore;
  int left = nq*nc, r = 0, ndrop = 0;
  auto round = [&](const int sg) {
    int na = 0;
    Kokkos::parallel_reduce("cln_fvj_round", pol,
    KOKKOS_LAMBDA(const int q, int &cnt) {
      const int w = q/nc, c = q%nc;
      if (lev(w,c) >= 0 || (sg > 0 && core(mds(w),c) != sg)) return;
      auto done = [&](const int v) {return lev(v/nc,v%nc) >= 0 && lev(v/nc,v%nc) < r;};
      for (int f = 0; f < 6; ++f) {
        const int v = up(w, c, f);
        if (v >= 0 && !(drop(w,c) & (1 << f)) && !done(v)) return;
        const int u = nbr(w, c, f);
        if (u < 0 || drop(u/nc,u%nc) == 0) continue;
        for (int g = 0; g < 6; ++g) {
          if ((drop(u/nc,u%nc) & (1 << g)) && up(u/nc, u%nc, g) == w*nc + c && !done(u)) {
            return;
          }
        }
      }
      lev(w,c) = r;
      ++cnt;
    }, na);
    return na;
  };
  while (left > 0) {
    const int na = round(0);
    if (na > 0) {
      left -= na;
      ++r;
      continue;
    }
    // a stall: the remaining cells al, less those with no remaining downstream cell
    par_for("cln_fvj_al", DevExeSpace(), 0, nw - 1, 0, nc - 1,
    KOKKOS_LAMBDA(const int w, const int c) {
      al(w,c) = (lev(w,c) < 0) ? 1 : 0;
    });
    for (int nt = 1; nt > 0;) {
      Kokkos::deep_copy(hd, 0);
      par_for("cln_fvj_hd", DevExeSpace(), 0, nq - 1, 0, nc - 1,
      KOKKOS_LAMBDA(const int w, const int c) {
        for (int f = 0; f < 6 && al(w,c); ++f) {
          const int v = up(w, c, f);
          if (v >= 0 && !(drop(w,c) & (1 << f)) && al(v/nc,v%nc)) hd(v/nc,v%nc) = 1;
        }
      });
      nt = 0;
      Kokkos::parallel_reduce("cln_fvj_trim", pol,
      KOKKOS_LAMBDA(const int q, int &cnt) {
        const int w = q/nc, c = q%nc;
        if (al(w,c) && !hd(w,c)) {
          al(w,c) = 0;
          ++cnt;
        }
      }, nt);
    }
    // (none left: every remaining cell)
    int nd = 0, all = 0;
    for (; all < 2 && nd == 0; ++all) {
      Kokkos::parallel_reduce("cln_fvj_drop", pol,
      KOKKOS_LAMBDA(const int q, int &cnt) {
        const int w = q/nc, c = q%nc;
        if (lev(w,c) >= 0 || !(all || al(w,c))) return;
        const int me = w*nc + c;
        for (int f = 0; f < 6; ++f) {
          const int v = up(w, c, f);
          if (v < 0 || (drop(w,c) & (1 << f)) || !(all || al(v/nc,v%nc)) ||
              lev(v/nc,v%nc) >= 0) continue;
          const Real kv = key(v), kc = key(me);
          if (kv > kc || (kv == kc && v > me)) {
            drop(w,c) |= (1 << f);
            ++cnt;
          }
        }
      }, nd);
    }
    if (nd == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/clean_fv_joint: the sweep levels do not close"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ndrop += nd;
    // the core is a segment of levels of its own, repeated by the sweep until it
    // converges (the cells downstream of it wait for it)
    const int sg = static_cast<int>(cln_fvseg.size()) + 1, fb = (all == 2);
    int ns = 0;
    Kokkos::parallel_reduce("cln_fvj_core", pol,
    KOKKOS_LAMBDA(const int q, int &cnt) {
      const int w = q/nc, c = q%nc;
      if (lev(w,c) < 0 && (fb || al(w,c))) {
        core(mds(w),c) = sg;
        ++cnt;
      }
    }, ns);
    const int r0 = r;
    for (int na = 0; ns > 0; ns -= na, left -= na, ++r) {
      na = round(sg);
      if (na == 0) break;
    }
    cln_fvseg.push_back({r0, r - 1, sg});
  }
  // the first level whose cells read another rank, and the cells that do
  auto ext = [=] KOKKOS_FUNCTION (const int w, const int c) {
    bool e = false;
    for (int f = 0; f < 6; ++f) e = e || (up(w, c, f) == -2);
    return e;
  };
  int l0 = 0, ne = 0;
  Kokkos::parallel_reduce("cln_fvj_l0", pol,
  KOKKOS_LAMBDA(const int q, int &l) {
    if (ext(q/nc, q%nc)) l = (lev(q/nc,q%nc) < l) ? lev(q/nc,q%nc) : l;
  }, Kokkos::Min<int>(l0));
  Kokkos::parallel_reduce("cln_fvj_ne", pol,
  KOKKOS_LAMBDA(const int q, int &cnt) {
    if (ext(q/nc, q%nc)) ++cnt;
  }, ne);
  cln_fvl0 = std::min(cln_fvl0, l0);
  cln_fvjn[0] += ndrop;
  cln_fvjn[1] += ne;
  return r;
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
//! the ghost exchange until the largest change is below clean_fv_tol (1e-10) of the
//! largest intensity).
//! With fv at most 8 passes per build (the slopes converge over the builds of a run).
//! The inner face emits the sweep's emission (vet_bc_inner emit) or the S of the first
//! cell and the boundary beams (implicit_beam), the top its Marshak bath, a reflecting
//! theta wall or x1 face is specular (CleanFvWalls), another physical theta or phi
//! (Cartesian x2, x3) face is closed.  Cartesian meshes take vet_sc's direction set.

void RadiationM1::CleanFvSweep() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const int n = cln_fvn, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  const bool axi = !pmy_pack->pmesh->three_d && pmy_pack->pmesh->use_spherical_polar;
  CleanArgs x = CleanMakeArgs();
  auto x3f = pmy_pack->pcoord->xx3f;
  // clean_fv_lag: one pass from the intensities of the last build (a restart's), the
  // ghosts refreshed from them first: the inflow across MeshBlock faces and periodic
  // seams lags by one build, the order inside a MeshBlock stays exact
  const int nfv = cln_fvb*n;
  const bool have = (cln_fvi.extent_int(1) == nfv) ||
                    (cln_fvf32 && cln_fvis.extent_int(1) == nfv);
  // (extrap: the history a restart left in the ghosts of the intensities it read, taken
  // before clean_fv_float converts them)
  if (cln_fvrst && cln_fvlag == 2 && cln_fvi.extent_int(1) == nfv) CleanFvRstRead();
  if (cln_fvf32) {
    if (axi) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/clean_fv_float does not take a 2-D "
                << "axisymmetric mesh" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // (a restart has read the intensities in double precision)
    if (cln_fvis.extent_int(1) != cln_fvb*n) {
      Kokkos::realloc(cln_fvis, nmb, cln_fvb*n, cln_s.extent_int(2), cln_s.extent_int(3),
                      cln_s.extent_int(4));
      if (cln_fvi.extent_int(1) == cln_fvb*n) {
        auto fs = cln_fvis;
        auto fd = cln_fvi;
        const int nv = cln_fvb*n, n2 = fs.extent_int(2) - 1;
        const int n3 = fs.extent_int(3) - 1, n4 = fs.extent_int(4) - 1;
        par_for("cln_fv_tof", DevExeSpace(), 0, nmb - 1, 0, nv - 1, 0, n2, 0, n3, 0, n4,
        KOKKOS_LAMBDA(const int m, const int v, const int k, const int j, const int i) {
          fs(m,v,k,j,i) = static_cast<float>(fd(m,v,k,j,i));
        });
      }
      cln_fvi = DvceArray5D<Real>();
    }
  } else if (cln_fvi.extent_int(1) != cln_fvb*n) {
    Kokkos::realloc(cln_fvi, nmb, cln_fvb*n, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
  }
  if (!have && ibc_x1max == M1_IBC_REFLECT) {
    // (a reflecting top mirrors the last build (CleanFvWalls): before the first one the
    // top cell's E)
    CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
    auto iw_ = iw;
    const int ie = indcs.ie;
    par_for("cln_fv_top0", DevExeSpace(), 0, nmb - 1, ks, indcs.ke, js, indcs.je, 0,
            n - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
      if (CleanBoundaryFace(x, 0, m, k, j, ie + 1)) {
        fi.Set(m, d, k, j, ie, iw_(m,M1_IW_EN,k,j,ie));
      }
    });
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
  if (std::find(cln_fvdo.begin(), cln_fvdo.end(), true) != cln_fvdo.end()) {
    CleanFvLevels();
  }
  // (the plain fv1 sweep only; extrap: a set new to the history restarts it)
  const bool lagon = (cln_fvlag > 0) && have && !cln_fvfull && !cln_fvjoint &&
                     cln_rfv == 1 && !axi;
  if (!lagon) cln_fvnv = 0;
  bool lag = lagon && (cln_fvlag == 1);
  CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto fa = cln_fva;
  auto fc = cln_fvc;
  auto fg = cln_fvg;
  auto fo = cln_fvo;
  auto fk = cln_fvk;
  const int nmd = nmb*n;
  auto fm = cln_fvm;
  const bool o2 = (cln_rfv == 2);
  const bool wall[2] = {CleanFvWall(0), CleanFvWall(1)};
  const bool rwl[2] = {ibc_x1min == M1_IBC_REFLECT, ibc_x1max == M1_IBC_REFLECT};
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
  if (o2) Kokkos::deep_copy(cln_fvp, cln_fvi);
  const CleanFvI fp = o2 ? CleanFvI{cln_fvp} : fi;
  // 2-D axisymmetric: a phi face's neighbour is the own cell in other directions, which
  // other threads write in the same pass: they are read from the last pass (cln_fvx)
  if (axi && cln_fvx.extent_int(1) != cln_fvb*n) {
    Kokkos::realloc(cln_fvx, fi.extent_int(0), cln_fvb*n, fi.extent_int(2),
                    fi.extent_int(3), fi.extent_int(4));
  }
  const CleanFvI fx = axi ? CleanFvI{cln_fvx} : fi;
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
    CleanGhostFill(cln_fvq, CLQ_X, 4);   // (chi, S, J, q)
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
  if (lagon) {
    if (cln_fvf32) {
      CleanGhostFill(cln_fvis, 0, cln_fvb*n);
    } else {
      CleanGhostFill(cln_fvi, 0, cln_fvb*n);
    }
    CleanFvPanelSlopes(false);
    CleanFvSeams();
    if (cln_fvlag == 2) lag = CleanFvExtrap();
  }
  cln_fvrst = false;
  // the walls' incoming beams of the intensities the sweep starts from (a restart has
  // no others)
  if (wall[0] || wall[1] || rwl[0] || rwl[1]) CleanFvWalls(true);
  // (each level reduces -(its largest change) and its largest intensity on the device:
  // no wait)
  const int nlev = static_cast<int>(cln_fvlv.size());
  Kokkos::View<Kokkos::MinMaxScalar<Real>*, DevMemSpace> red("cln_fv_red", nlev);
  // fv1: a pass after the first takes only the cells at the faces of the MeshBlock and
  // those whose upwind neighbour changed after it took them (due: the pass that must
  // take a cell, done: the pass that has); the others would repeat their intensity
  const bool act = !o2 && !ld && !axi;
  const int ie = indcs.ie, je = indcs.je, ke = indcs.ke, nc = nx1*nx2*indcs.nx3;
  // (clean_pc = x1line with clean_pc_factor, the factors rebuilt every pass: the flags
  // live in the factor store, which no solve reads before its next rebuild)
  bool shf = act && cln_cell && cln_pc == CLN_PC_X1LINE && cln_pcfac &&
             cln_pcreuse == 1;
  const size_t nf = static_cast<size_t>(nmb)*n*nc;
  if (shf) {
    CleanCellFactorStore();
    const size_t nb = cln_pcf32 ? cln_lfs.span()*sizeof(float)
                                : cln_lf.span()*sizeof(Real);
    shf = (nb >= 2*nf);   // (a store too small for the two flag arrays)
  }
  if (shf) {
    int8_t *pf = cln_pcf32 ? reinterpret_cast<int8_t *>(cln_lfs.data())
                           : reinterpret_cast<int8_t *>(cln_lf.data());
    cln_fvdue = DvceArray2D<int8_t>(pf, nmb*n, nc);
    cln_fvdone = DvceArray2D<int8_t>(pf + nf, nmb*n, nc);
    Kokkos::deep_copy(cln_fvdue, -1);
    cln_lf_ok = false;
  } else if (act && (cln_fvdue.extent_int(0) != nmb*n || cln_fvdue.extent_int(1) != nc ||
                     !cln_fvdue.use_count())) {
    Kokkos::realloc(cln_fvdue, nmb*n, nc);
    Kokkos::realloc(cln_fvdone, nmb*n, nc);
  }
  auto due = cln_fvdue;
  auto done = cln_fvdone;
  cln_fvnp[0] += 1;
  // clean_fv_joint: a later pass starts at the first level that reads the last pass
  const bool jt = cln_fvjoint;
  auto fsrc = cln_fvsrc;
  auto fcore = cln_fvcore;
  int lseg = 1;
  for (const auto &sg : cln_fvseg) lseg = std::max(lseg, sg[1] - sg[0] + 1);
  Kokkos::View<Kokkos::MinMaxScalar<Real>*, DevMemSpace> rsg("cln_fv_rsg", lseg);
  using RedView = decltype(Kokkos::subview(rsg, 0));
  for (int it = 0; it < 100; ++it) {
    Real dmax = 0.0, imax = 0.0;
    if (axi) Kokkos::deep_copy(cln_fvx, cln_fvi);
    if (act) Kokkos::deep_copy(done, -1);
    const int l1 = (jt && it > 0) ? std::min(cln_fvl0, nlev) : 0;
    // one level (sg > 0: only the cells of cycle segment sg) into the reduction rv; pl
    // std::true_type: the plain fv1 sweep (CUDA, which takes no device lambda inside a
    // generic one: a flag, the branches of the other forms kept)
#if defined(KOKKOS_ENABLE_CUDA)
    auto level = [&](const int l, const int sg, const RedView rv, const bool P) {
#else
    auto level = [&](const int l, const int sg, const RedView rv, auto pl) {
#endif
    const int64_t k0 = cln_fvkl[l];
    const int nkb = static_cast<int>(cln_fvkl[l+1] - k0);
    Kokkos::parallel_reduce("cln_fv", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                            cln_fvlv[l]),
    KOKKOS_LAMBDA(const int q, Kokkos::MinMaxScalar<Real> &mm) {
#if !defined(KOKKOS_ENABLE_CUDA)
      // (pl: the plain fv1 sweep, the branches of the other forms compiled out)
      constexpr bool P = decltype(pl)::value;
#endif
      const bool o2c = !P && o2, ldc = !P && ld, axc = !P && axi, acc = P || act;
      const int nlc = ldc ? nl : 0;
      // (the m n + d of entry q: the last whose start in the level is <= q, between
      // those of its 64-entry block and of the next)
      int md = fk(k0 + q/64);
      for (int h = (q/64 + 1 < nkb) ? fk(k0 + q/64 + 1) + 1 : nmd; h - md > 1;) {
        const int mid = (md + h)/2;
        if (fg(l,mid) <= q) {md = mid;} else {h = mid;}
      }
      const int m = md/n, d = md%n, c = fc(fo(l,md) + q);
      if (sg > 0 && fcore(md,c) != sg) return;
      const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
      {
        const int i = is + c%nx1, j = js + (c/nx1)%nx2, k = ks + c/(nx1*nx2);
        if (acc && it > 0 && due(md,c) != it && i != is && i != ie && j != js &&
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
          const bool ax = axc && dir == 2;
          // (Cartesian: no faces along an axis of one cell)
          if (!x.sph && !x.cub && (dir == 1 ? !x.multi_d : dir == 2 && !x.three_d)) {
            continue;
          }
          const bool bdry = !ax && CleanBoundaryFace(x, dir, m, kf, jf, i_f);
          // a reflecting theta wall or x1 face: specular (its incoming beams from
          // cln_fvw, by i or j)
          const int sw = !bdry ? -1 : (dir == 1 && wall[up]) ? up
                       : (dir == 0 && rwl[up]) ? 2 + up : -1;
          if (bdry && dir != 0 && sw < 0) continue;
          const Real sg = up ? 1.0 : -1.0;
          if (sw >= 0 && a < 0.0) {
            const int cw = dir ? i : j;
            in -= a*fw(m,sw,k,cw,d);
            for (int e = 0; e < 3 && ldc; ++e) {
              rd[e] -= a*((e == dir) ? sg*fw(m,sw,k,cw,d) : fw(m,sw,k,cw,(1+e)*n+d));
            }
            hd[dir] -= sg*a;
            continue;
          }
          if (a >= 0.0) {
            if (!bdry) om |= 1 << f;
            out += a;
            cd[dir] += sg*a;
            od[dir] += a;
            po[dir] |= up ? 1 : 2;
            if (o2c) {
              const int sg = up ? 1 : -1;
              ou[no] = a;
              ud[no][0] = CleanFvNbr(ax ? fx : fi, fm, ax, m, d, dir, k, j, i, -sg);
              ud[no][1] = CleanFvNbr(fp, fm, ax, m, d, dir, k, j, i, sg);
              ++no;
            }
          } else {
            // the upwind trace (lld: the neighbour's transverse slopes too)
            const CleanFvI fd = ldc ? CleanFvI{fq} : fp;
            // (clean_fv_joint: across a MeshBlock face the cell its ghost copies)
            const int s = up ? 1 : -1;
            const int sv = (jt && !bdry) ? fsrc(m,k+s*(dir == 2),j+s*(dir == 1),
                                                i+s*(dir == 0)) : -1;
            const int cv = sv%nc;
            Real t = bdry ? (up ? ihi : (emi ? iemit : fq(m,CLQ_S0,k,j,i)))
                   : ax ? CleanFvFaceAxi(fx, fd, fm, o2c, m, d, k, j, i, up, !up, nlc)
                   : (sv >= 0) ? fi(sv/nc,d,ks+cv/(nx1*nx2),js+(cv/nx1)%nx2,is+cv%nx1)
                   : CleanFvFace(fi, fd, o2c, m, d, dir, kf, jf, i_f, !up, nlc);
            for (int q = 0; q < bmn && bdry && !up; ++q) {
              if (d == bmr[q]) t += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
            }
            in -= a*t;
            hd[dir] -= sg*a;
            for (int e = 0; e < 3 && ldc; ++e) {
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
        if (ldc && out + vc > 0.0) {
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
        iv = fi.Round(fmax(iv, 0.0));
        mm.min_val = fmin(mm.min_val, -fabs(iv - fi(m,d,k,j,i)));
        mm.max_val = fmax(mm.max_val, iv);
        // (a change: the downstream neighbours are due in this pass, or the next one if
        // this pass has taken them)
        const bool chg = acc && !(iv == fi(m,d,k,j,i));
        fi.Set(m, d, k, j, i, iv);
        if (acc) done(md,c) = it;
        for (int f = 0; f < 6 && chg; ++f) {
          const int x0 = (f < 2) ? i - is : ((f < 4) ? j - js : k - ks);
          const int x1 = (f < 2) ? ie - is : ((f < 4) ? je - js : ke - ks);
          if (!(om & (1 << f)) || x0 == ((f%2) ? x1 : 0)) continue;   // (in the block)
          const int s = (f%2) ? 1 : -1, st = (f < 2) ? 1 : ((f < 4) ? nx1 : nx1*nx2);
          const int cn = c + s*st;
          due(md,cn) = (done(md,cn) == it) ? it + 1 : it;
        }
        for (int e = 0; e < 3 && ldc; ++e) fi.Set(m, (1+e)*n+d, k, j, i, sl[e]);
      }
    }, Kokkos::MinMax<Real, DevMemSpace>(rv));
    };
    auto lev = [&](const int l, const int sg, const RedView rv) {
      if (act) {
        level(l, sg, rv, std::true_type());
      } else {
        level(l, sg, rv, std::false_type());
      }
    };
    for (int l = l1; l < nlev; ++l) {
      lev(l, 0, Kokkos::subview(red, l));
      // clean_fv_joint: a cycle segment ending here is repeated until it converges
      for (const auto &sg : cln_fvseg) {
        for (int rep = 0; jt && sg[1] == l && rep < 100; ++rep) {
          for (int q = sg[0]; q <= l; ++q) {
            lev(q, sg[2], Kokkos::subview(rsg, q - sg[0]));
          }
          auto s_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rsg);
          Real dm = 0.0, im = 0.0;
          for (int q = 0; q <= l - sg[0]; ++q) {
            dm = fmax(dm, -s_h(q).min_val);
            im = fmax(im, s_h(q).max_val);
          }
          cln_fvnp[2] += l - sg[0] + 1;
          if (dm <= cln_fvtol*im) break;
        }
      }
    }
    auto r_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), red);
    for (int l = 0; l < nlev; ++l) {
      if (l >= l1) dmax = fmax(dmax, -r_h(l).min_val);
      imax = fmax(imax, r_h(l).max_val);
    }
    if (cln_fvf32) {
      CleanGhostFill(cln_fvis, 0, cln_fvb*n);
    } else {
      CleanGhostFill(cln_fvi, 0, cln_fvb*n);
    }
    CleanFvPanelSlopes(false);
    CleanFvSeams();
    if (wall[0] || wall[1] || rwl[0]) CleanFvWalls(false);
    cln_fvnp[1] += 1;
    cln_fvnp[2] += nlev - l1;
    if (jt && cln_fvexact && !rwl[0]) break;
    // (extrap: one pass if the inflow it extrapolated was close to the new one)
    if (lag && (cln_fvlag == 1 ||
                CleanFvExtrapMiss() <= cln_fvlagtol*CleanReduce(imax, true))) {
      break;
    }
    if (lag) {
      lag = false;
      continue;
    }
    if (CleanReduce(dmax, true) <= cln_fvtol*CleanReduce(imax, true)) break;
  }
  // extrap: this build's inflow is the newest of the history (a rebuild at the same time
  // replaces it)
  if (cln_fvlag == 2 && (cln_fvnv == 0 || cln_fvtc > cln_fvt1)) {
    cln_fvt1 = cln_fvtc;
    cln_fvnv = std::min(cln_fvnv + 1, 2);
  }
}

//----------------------------------------------------------------------------------------
//! \fn CleanShell CleanFvLayer
//! \brief clean_fv_lag = extrap: the first ghost layer of a MeshBlock (the cells an fv1
//! sweep reads across its faces) as the shell of the active cells in the box one larger;
//! cell q is (k,j,i) + (ks, js, is) - (td, md, 1) of the mesh indices

CleanShell CleanFvLayer(const RegionIndcs &ix, const bool md, const bool td) {
  return CleanShell(ix.nx3 + (td ? 2 : 0), ix.nx2 + (md ? 2 : 0), ix.nx1 + 2, 0, 1,
                    ix.nx1, md ? 1 : 0, md ? ix.nx2 : 0, td ? 1 : 0, td ? ix.nx3 : 0, md,
                    td);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::CleanFvExtrap
//! \brief clean_fv_lag = extrap, before the pass (the ghosts hold the inflow G1 of the
//! last build at t1): with the build before (G2 at t2, cln_fvgp) the ghosts take
//! G1 + (G1 - G2) (tc - t1)/(t1 - t2), at least 0, and G1 becomes the history.  Returns
//! whether it extrapolated (one pass); the first two builds of a set, a restart without
//! the record and a step taken back (dt_retry) sweep in full.

bool RadiationM1::CleanFvExtrap() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const CleanShell sh = CleanFvLayer(indcs, md, td);
  const int nmb = pmy_pack->nmb_thispack, n = cln_fvn;
  const int oi = indcs.is - 1, oj = indcs.js - (md ? 1 : 0), ok = indcs.ks - (td ? 1 : 0);
  if (cln_fvgp.extent_int(0) != nmb || cln_fvgp.extent_int(1) != n ||
      cln_fvgp.extent_int(2) != sh.ncell) {
    Kokkos::realloc(cln_fvgp, nmb, n, sh.ncell);
    cln_fvnv = std::min(cln_fvnv, 1);
  }
  const Real tc = cln_fvtc;
  if (tc < cln_fvt1) cln_fvnv = 0;
  const bool ex = (cln_fvnv >= 2) && (cln_fvt1 > cln_fvt2);
  const bool rec = (cln_fvnv >= 1) && (tc > cln_fvt1);
  if (!ex && !rec) return false;
  const Real w = ex ? (tc - cln_fvt1)/(cln_fvt1 - cln_fvt2) : 0.0;
  if (ex && (cln_fvgx.extent_int(1) != n || cln_fvgx.extent_int(2) != sh.ncell)) {
    Kokkos::realloc(cln_fvgx, nmb, n, sh.ncell);
  }
  const CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto gp = cln_fvgp;
  auto gx = cln_fvgx;
  par_for("cln_fv_extrap", DevExeSpace(), 0, nmb - 1, 0, n - 1, 0, sh.ncell - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int q) {
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    const Real g1 = fi(m,d,k+ok,j+oj,i+oi);
    if (ex) {
      const Real ge = fi.Round(fmax(g1 + w*(g1 - gp(m,d,q)), 0.0));
      fi.Set(m, d, k+ok, j+oj, i+oi, ge);
      gx(m,d,q) = ge;
    }
    if (rec) gp(m,d,q) = g1;
  });
  if (rec) cln_fvt2 = cln_fvt1;
  return ex;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanFvExtrapMiss
//! \brief clean_fv_lag = extrap, after the pass and its ghost exchange: the largest
//! difference (all ranks) between the new inflow and the extrapolated one the pass took,
//! over the ghosts beyond MeshBlock faces and periodic seams

Real RadiationM1::CleanFvExtrapMiss() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const CleanShell sh = CleanFvLayer(indcs, md, td);
  const int nmb = pmy_pack->nmb_thispack, n = cln_fvn;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int oi = is - 1, oj = js - (md ? 1 : 0), ok = ks - (td ? 1 : 0);
  const CleanFvI fi = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto gx = cln_fvgx;
  auto bc = pmy_pack->pmb->mb_bcs.d_view;
  auto open = [] KOKKOS_FUNCTION (const BoundaryFlag f) {
    return f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
           f == BoundaryFlag::panel;
  };
  Real miss = 0.0;
  const int nq = sh.ncell;
  Kokkos::parallel_reduce("cln_fv_miss", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*n*nq),
  KOKKOS_LAMBDA(const int idx, Real &mx) {
    const int m = idx/(n*nq), d = (idx/nq)%n, q = idx%nq;
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    k += ok;
    j += oj;
    i += oi;
    // (a ghost beyond a physical face is not inflow)
    if ((i < is && !open(bc(m,BoundaryFace::inner_x1))) ||
        (i > ie && !open(bc(m,BoundaryFace::outer_x1))) ||
        (md && j < js && !open(bc(m,BoundaryFace::inner_x2))) ||
        (md && j > je && !open(bc(m,BoundaryFace::outer_x2))) ||
        (td && k < ks && !open(bc(m,BoundaryFace::inner_x3))) ||
        (td && k > ke && !open(bc(m,BoundaryFace::outer_x3)))) {
      return;
    }
    mx = fmax(mx, fabs(fi(m,d,k,j,i) - gx(m,d,q)));
  }, Kokkos::Max<Real>(miss));
  return CleanReduce(miss, true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvRstWrite
//! \brief clean_fv_lag = extrap: the history in the restart's copy a of cln_fvi (channels
//! c .., MeshBlocks m0 .. m0+nmb-1): each first-layer ghost of cln_fvgp in the second
//! layer one cell further out (an unused copy of the ghosts a restart refills), and
//! (magic, t1, t2, nv) in the ghost (ks, js-1, is-2) of the first four channels.  Needs a
//! multi-D mesh with two ghost layers; else a restart sweeps its first builds in full.

constexpr Real kCleanFvRstMagic = -7.5e7;

void RadiationM1::CleanFvRstWrite(DvceArray5D<Real> &a, const int c, const int nmb,
                                  const int m0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const int n = cln_fvn;
  if (cln_fvlag != 2 || cln_fvnv < 1 || !md || indcs.ng < 2 || n < 4 ||
      cln_fvgp.extent_int(1) != n) {
    return;
  }
  const CleanShell sh = CleanFvLayer(indcs, md, td);
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int oi = is - 1, oj = js - 1, ok = ks - (td ? 1 : 0);
  auto gp = cln_fvgp;
  auto a_ = a;
  par_for("cln_fv_rstw", DevExeSpace(), 0, nmb - 1, 0, n - 1, 0, sh.ncell - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int q) {
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    k += ok;
    j += oj;
    i += oi;
    k += td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    j += (j < js) ? -1 : ((j > je) ? 1 : 0);
    i += (i < is) ? -1 : ((i > ie) ? 1 : 0);
    a_(m,c+d,k,j,i) = gp(m0+m,d,q);
  });
  const Real v[4] = {kCleanFvRstMagic, cln_fvt1, cln_fvt2, static_cast<Real>(cln_fvnv)};
  for (int e = 0; e < 4; ++e) {
    const Real ve = v[e];
    par_for("cln_fv_rsts", DevExeSpace(), 0, nmb - 1, KOKKOS_LAMBDA(const int m) {
      a_(m,c+e,ks,js-1,is-2) = ve;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFvRstRead
//! \brief clean_fv_lag = extrap: the history CleanFvRstWrite left in the ghosts of the
//! intensities a restart read (none: nv = 0, the first builds sweep in full)

void RadiationM1::CleanFvRstRead() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const int n = cln_fvn, nmb = pmy_pack->nmb_thispack;
  cln_fvnv = 0;
  if (!md || indcs.ng < 2 || n < 4) return;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  // (the record through a contiguous copy: a strided device subview has no copy to the
  // host on CUDA)
  auto fi = cln_fvi;
  DvceArray1D<Real> rec("cln_fv_rec", 4);
  par_for("cln_fv_rec", DevExeSpace(), 0, 3, KOKKOS_LAMBDA(const int e) {
    rec(e) = fi(0,e,ks,js-1,is-2);
  });
  auto v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rec);
  if (!(v(0) == kCleanFvRstMagic)) return;
  cln_fvt1 = v(1);
  cln_fvt2 = v(2);
  cln_fvnv = static_cast<int>(v(3));
  const CleanShell sh = CleanFvLayer(indcs, md, td);
  Kokkos::realloc(cln_fvgp, nmb, n, sh.ncell);
  const int oi = is - 1, oj = js - 1, ok = ks - (td ? 1 : 0);
  auto gp = cln_fvgp;
  par_for("cln_fv_rstr", DevExeSpace(), 0, nmb - 1, 0, n - 1, 0, sh.ncell - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int q) {
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    k += ok;
    j += oj;
    i += oi;
    k += td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    j += (j < js) ? -1 : ((j > je) ? 1 : 0);
    i += (i < is) ? -1 : ((i > ie) ? 1 : 0);
    gp(m,d,q) = fi(m,d,k,j,i);
  });
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

void RadiationM1::CleanRayCut() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  vlat_icut = 0;
  if (!cln_rcut_on || !cln_init || cln_rfv > 0 || cln_cell ||
      pmy_pack->pmesh->mesh_indcs.nx1 != indcs.nx1) return;
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

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFromGd
//! \brief the closure from the vet_gd intensities in the local basis (r, theta, phi) of
//! each cell, mu_b = n . e_b: J = <I>, K_bc = <I mu_b mu_c>, H+-_b = <I mu_b> over mu_b
//! >< 0, K+^d_b = <I mu_b mu_d> over mu_d > 0 and C1 = <|I - S|>/J (CLN_* layout, as
//! CleanClosureFromRays; no ray face fluxes).  On the cubed sphere mu_b = n_b (Cartesian
//! D and k+ columns) and the half ranges are about the family normals (CleanCsFamily).

void RadiationM1::CleanClosureFromGd(const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (vgd_pk && i0 < 0) return;   // (planb-1117: done per shell in the sweep)
  int is = (i0 < 0) ? indcs.is : i0, ie = (i0 < 0) ? indcs.ie : i1;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool fv = (cln_rfv > 0);
  const int n = fv ? cln_fvn : vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  const CleanFvI fi_ = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
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
  // (in device memory, CleanGdFaceSums; vgd_pk: one copy per sweep, at its first shell)
  if (i0 < 0 || i0 == indcs.is) cln_gdx = CleanDevArgs(x, 3);
  const CleanArgs *px = cln_gdx;
  // (clean_ray_cut: below the cut the sc rays' intensity is the diffusion limit)
  const int icut = indcs.is + (vgd_on ? vlat_icut : 0);
  const Real ic = 1.0/c_light;
  par_for("cln_gd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px;
    const Real st = sin(x2v(m,j)), ct = cos(x2v(m,j));
    const Real sp = sin(x3v(m,k)), cp = cos(x3v(m,k));
    const Real s0 = fv ? fq(m,CLQ_S0,k,j,i) : exp(cw_(m,1,k+ok,j+og,i));
    // the cubed sphere: D in Cartesian components, half ranges about the families n_d
    Real nd[3][3];
    if (x.cub) CleanCsFamily(x, m, k, j, i, nd);
    Real q[CLN_N] = {0.0}, hp[3] = {0.0}, hm[3] = {0.0};   // (q in registers)
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      Real mu[3] = {nx*st*cp + ny*st*sp + nz*ct, nx*ct*cp + ny*ct*sp - nz*st,
                    -nx*sp + ny*cp}, nu[3];
      for (int b = 0; b < 3; ++b) {
        if (!x.sph) mu[b] = dir_(d,b);
        nu[b] = x.cub ? nx*nd[b][0] + ny*nd[b][1] + nz*nd[b][2] : mu[b];
      }
      const Real nv[3] = {nx, ny, nz};
      const Real iv = fv ? fi_(m,d,k,j,i) : (i < icut)
                    ? CleanThickI(x, iw_, x3v, ic, nv, m, k, j, i) : vi_(m,d,k+ok,j+og,i);
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

//----------------------------------------------------------------------------------------
//! \struct CleanGdI
//! \brief the sc rays' intensities in the mesh layout without a copy (CleanGdFaceSums):
//! vgd_i, the low theta and phi ghost layers from gl (side 0: j = js - 1, side 1:
//! k = ks - 1); planb-1117 (pk): those from the sweep's raw band or the on-rank owner
//! (own), as the sweep reads them (rd0t), before the wall map

struct CleanGdI {
  RadiationM1::VgdRag vi;
  DvceArray5D<Real> gl;
  int ok, og, ks, js;
  bool pk;
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> own;
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int d, const int k, const int j, const int i) const {
    if (pk && (j < js || k < ks)) {
      const int kb = k + ok, jb = j + og, ow = own(m,kb,jb,0);
      return (ow < 0) ? vi(m,d,kb,jb,i) : vi(ow,d,own(m,kb,jb,1),own(m,kb,jb,2),i);
    }
    if (j < js) return gl(m,d,0,k-ks,i);
    if (k < ks) return gl(m,d,1,j-js,i);
    return vi(m,d,k+ok,j+og,i);
  }
};

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGdFaceSums
//! \brief the vet_gd rays' face fluxes for clean_consistency (CLN_FX, CLN_F2, CLN_F3,
//! CLN_F1: sum_d w (n_d . A_f) I_f / A_f and the stress columns in the face's local
//! basis), A_f the exact vector area of the face (closed for every cell).  A theta or
//! phi face carries the intensity of its upwind cell (none on a physical boundary), and
//! so does an r face of the FV rays and of the staggered layout.  Else the
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

void RadiationM1::CleanGdFaceSums(const int s0, const int s1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  // planb-1117 (vgd_pk): s0 >= 0 the lateral sums of shells s0..s1 and the r faces they
  // close (VetGdPkHook); s0 < 0 the rest
  const bool pk = vgd_pk, one = (s0 >= 0);
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool fv = (cln_rfv > 0), o2 = (cln_rfv == 2);
  const int n = fv ? cln_fvn : vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  CleanFvI fi_ = CleanFvI::Of(cln_fvi, cln_fvis, cln_fvf32);
  auto fm_ = cln_fvm;
  const int nl = (cln_rfv == 3) ? n : 0;
  const CleanFvI fq_ = (nl > 0) ? CleanFvI{cln_fvq} : fi_;
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
  const bool f1 = cln_s.extent_int(1) > CLN_F1;   // (the cell layout's CLN_F1)
  CleanArgs x = CleanMakeArgs();
  // (the kernels read the args in device memory: a light-weight launch without a host
  // wait, per shell too; vgd_pk: one copy per sweep, CleanClosureFromGd)
  if (!one) cln_gdx = CleanDevArgs(x, 3);
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
  const bool rwl[2] = {fv && !open, fv && ibc_x1max == M1_IBC_REFLECT};
  const bool rwsc[2] = {!fv && !open, !fv && ibc_x1max == M1_IBC_REFLECT};
  auto fw_ = cln_fvw;
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
  // the short characteristics' intensities in the mesh layout with the ghost exchange
  // (the image cell's beyond a periodic seam, re-indexed by the wall map as the sweep)
  // (clean_rot_ring: two slots, of the sets blk mod 2, on the cubed sphere for the face
  // half ranges (CleanScUnion2), on sp with a reflecting top for its mirror)
  const int ns2 = (cln_ring && (cln_cssc || cln_scmem)) ? 2 : 1;
  const int doff = (ns2 == 1) ? 0 : cln_rqm*n;
  // Staggered spherical polar (one radial MeshBlock, no reflecting top): vgd_i itself,
  // and only the low theta and phi ghost layers exchanged, in chunks of clean_halo_vars
  // directions (one exchange each)
  const bool gl = !fv && x.sph && !cln_cell && !axi && pmy_pack->pmesh->three_d &&
                  vgd_nb1 == 1 && !cln_scmem;
  if (gl && vgd_on && vlat_icut > 0) {   // (clean_ray_cut: below it the diffusion limit)
    auto iw_ = iw;
    const Real ic = 1.0/c_light;
    par_for("cln_gd_thk", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, is + vlat_icut - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      for (int d = 0; d < n; ++d) {
        const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
        vi_(m,d,k+ok,j+og,i) = CleanThickI(x, iw_, x3v, ic, nv, m, k, j, i);
      }
    });
  }
  CleanGdI gi{vi_, DvceArray5D<Real>(), ok, og, ks, js, pk, vgd_own};
  if (gl && pk) {
  } else if (gl) {
    const int ng = indcs.ng;
    const int nc = std::min(n, cln_hmax), nt = std::max(ke - ks + 1, je - js + 1);
    gi.gl = DvceArray5D<Real>("cln_gd_gl", nmb1 + 1, n, 2, nt, cs.extent_int(4));
    DvceArray5D<Real> ch(Kokkos::view_alloc(Kokkos::WithoutInitializing, "cln_gd_ch"),
                         nmb1 + 1, nc, cs.extent_int(2), cs.extent_int(3),
                         cs.extent_int(4));
    auto gg = gi.gl;
    for (int d0 = 0; d0 < n; d0 += nc) {
      const int nd = std::min(nc, n - d0);
      // (only the ng layers along the MeshBlock's theta and phi faces: the exchange reads
      // no other cell)
      par_for("cln_gd_chc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (j >= js + ng && j <= je - ng && k >= ks + ng && k <= ke - ng) return;
        for (int c = 0; c < nd; ++c) ch(m,c,k,j,i) = vi_(m,d0+c,k+ok,j+og,i);
      });
      CleanGhostFill(ch, 0, nd);
      par_for("cln_gd_glc", DevExeSpace(), 0, nmb1, 0, nd - 1, 0, nt - 1, is, ie,
      KOKKOS_LAMBDA(const int m, const int c, const int t, const int i) {
        if (t <= ke - ks) gg(m,d0+c,0,t,i) = ch(m,c,ks+t,js-1,i);
        if (t <= je - js) gg(m,d0+c,1,t,i) = ch(m,c,ks-1,js+t,i);
      });
    }
  } else if (!fv) {
    if (cln_fvi.extent_int(1) != ns2*n) {
      Kokkos::realloc(cln_fvi, nmb1 + 1, ns2*n, cs.extent_int(2),
                      cs.extent_int(3), cs.extent_int(4));
    }
    fi_ = CleanFvI{cln_fvi};
    par_for("cln_gd_copy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < n; ++d) fi_.Set(m, doff + d, k, j, i, vi_(m,d,k+ok,j+og,i));
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
  auto lat = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d, const Real a[3]) {
    const CleanArgs &x = *px;
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
  // a mirror at the r face i, per unit area: no energy, the normal momentum of the
  // beams of cell ic that hit it (f of sign sg) doubled, no tangential part
  auto mirror = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                     const int ic, const Real sg, Real q[4]) {
    const CleanArgs &x = *px;
    Real a[3], e[3][3], w[4] = {0.0}, nf[3], wn = 0.0;
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      const Real tr = inten(m, d, k, j, ic) +
                      ((nl > 0) ? sg*fi_(m,nl+d,k,j,ic) + fq_(m,CLQ_C,k,j,ic) : 0.0);
      if (sg*f > 0.0) add(e, x.area[0](m,k,j,i), d, 2.0*f*tr, w);
    }
    for (int b = 0; b < 3; ++b) nf[b] = CleanNf(x, b, 0, m, k, j, i);
    for (int b = 0; b < 3; ++b) wn += w[1+b]*nf[b];
    q[0] = 0.0;
    for (int b = 0; b < 3; ++b) q[1+b] = wn*nf[b];
  };
  for (int dir = 1; dir < 3; ++dir) {
    const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
    if ((!x.sph && dir >= x.nf) || (pk && !one)) continue;
    par_for("cln_gd_lat", DevExeSpace(), 0, nmb1, ks, ke, js, je, one ? s0 : is,
            one ? s1 : ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      // the cubed sphere: the high face of the MeshBlock too (CLN_FH of its last cell)
      const bool hi = x.cub && ((dir == 1) ? j == je : k == ke);
      for (int up = 0; up <= (hi ? 1 : 0); ++up) {
        const int kf = k + up*(dir == 2), jf = j + up*(dir == 1);
        Real q[4] = {0.0}, a[3], e[3][3];
        face(dir, m, kf, jf, i, a, e);
        const bool bd = bnd(dir, m, kf, jf, i);
        const Real ar = x.area[dir](m,kf,jf,i);
        for (int d = 0; d < n; ++d) {
          add(e, ar, d, bd ? 0.0 : lat(dir, m, kf, jf, i, d, a), q);
        }
        for (int c = 0; c < 4; ++c) cs(m,(up ? CLN_FH + 4*(dir - 1) : nf)+c,k,j,i) = q[c];
      }
    });
  }
  if (!one) CleanGhostFill(cln_s, CLN_F2, 4);
  if (!one) CleanGhostFill(cln_s, CLN_F3, 4);
  // the conservative rays (clean_rays = fv | fv1), and the staggered layout's sc rays on
  // spherical polar (rr): the r face i carries their own upwind fluxes (area times the
  // per-area sums), not the balance march, which loses the diffusion flux in the
  // difference V chi (S - J) of thick cells and, with the bracket, feeds back on E
  const bool rr = fv || (x.sph && !cln_cell);
  auto rflux = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                    Real q[4]) {
    const CleanArgs &x = *px;
    Real a[3], e[3][3];
    for (int c = 0; c < 4; ++c) q[c] = 0.0;
    face(0, m, k, j, i, a, e);
    // (the sc rays at a reflecting x1 face: the mirror images of the cell's leaving
    // beams, CleanFvMirror, scaled so that the face returns the energy it receives)
    const int ic = (i == is) ? is : ie;
    const Real sg = (i == is) ? -1.0 : 1.0, an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const Real t[3] = {a[0]/an, a[1]/an, a[2]/an};
    const bool msc = !fv && CleanBoundaryFace(x, 0, m, k, j, i) && rwsc[i != is];
    Real sm = 1.0, ho = 0.0, hm = 0.0, u;
    for (int d = 0; d < n && msc; ++d) {
      const int r = CleanFvMirror(dir_, 0, n, d, t, sg, u);
      if (r < 0) ho += dir_(d,3)*u*inten(m, d, k, j, ic);
      if (r >= 0) hm += dir_(d,3)*u*inten(m, r, k, j, ic);
    }
    if (ho > 0.0 && hm > 0.0) sm = ho/hm;
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      Real iv = gl ? CleanFvFace(gi, gi, false, m, d, 0, k, j, i, f > 0.0)
                   : CleanFvFace(fi_, fq_, o2, m, doff + d, 0, k, j, i, f > 0.0, nl);
      if (CleanBoundaryFace(x, 0, m, k, j, i)) {   // the sweep's x1 boundary
        const Real s0 = fv ? sq_(m,CLQ_S0,k,j,is) : exp(cw_(m,1,k+ok,j+og,is));
        iv = (i == is) ? ((f > 0.0) ? (emi ? iemit : s0) : iv)
                       : ((f > 0.0) ? iv : ihi);
        for (int q = 0; q < bmn && i == is && f > 0.0; ++q) {
          if (d == bmr[q]) iv += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
        }
        // (a reflecting one of the FV rays: their mirrored beams, CleanFvWalls)
        const int up = (i != is);
        if (rwl[up] && (f > 0.0) != up) iv = fw_(m,2+up,k,j,d);
        const int r = msc ? CleanFvMirror(dir_, 0, n, d, t, sg, u) : -1;
        if (r >= 0) iv = sm*inten(m, r, k, j, ic);
      }
      add(e, x.area[0](m,k,j,i), d, f*iv, q);
    }
    for (int c = 0; c < 4; ++c) q[c] *= x.area[0](m,k,j,i);
  };
  if (!x.sph && !x.cub) {   // Cartesian: every x1 face carries the rays' upwind flux
    par_for("cln_gd_x1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;
      Real q[4];
      for (int up = 0; up <= ((i == ie) ? 1 : 0); ++up) {
        rflux(m, k, j, i + up, q);
        for (int c = 0; c < 4; ++c) {
          const Real f = q[c]/x.area[0](m,k,j,i+up);
          if (!up) cs(m,CLN_FX+c,k,j,i) = f;
          if (((!up && i == is) || up) && f1) cs(m,CLN_F1+c,k,j,i) = f;
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
  // rr: the rays' r face fluxes, all faces at once (the columns only read them)
  // (pk: one array over the sweep, freed below)
  if (pk && one && s0 == is) {
    cln_gdrf = DvceArray5D<Real>("cln_gd_rf", nmb1 + 1, 4, ke - ks + 1, nj, ie - is + 2);
  }
  DvceArray5D<Real> rf = pk ? cln_gdrf : DvceArray5D<Real>("cln_gd_rf", rr ? nmb1 + 1 : 1,
                       4, rr ? ke - ks + 1 : 1, rr ? nj : 1, rr ? ie - is + 2 : 1);
  if (pk && !one) cln_gdrf = DvceArray5D<Real>();
  if (rr && (one || !pk)) {
    par_for("cln_gd_rf", DevExeSpace(), 0, nmb1, ks, ke, js, je, one ? s0 : is,
            one ? ((s1 == ie) ? ie + 1 : s1) : ie + 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4];
      rflux(m, k, j, i, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
  }
  if (one) return;
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
    par_for("cln_gd_cel", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
    par_for("cln_gd_col", DevExeSpace(), 0, nmb1, ks, ke, js, je,
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
          if (!fv) rflux(m, k, j, i + 1, q);
          for (int c = 0; c < 4 && fv; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
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
  par_for("cln_rays", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
    CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
    CleanGhostFill(cln_v, 0, 2);
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

}  // namespace radm1
