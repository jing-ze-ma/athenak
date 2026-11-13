//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_lz.cpp
//! \brief radlz-1113 <rad_m1>/implicit_lz_form: the axial angular momentum L_z of gas +
//! radiation conserved by construction on the cubed sphere (AM_CONSERVING_CS_1112.md,
//! sect. 10).
//!
//! Under implicit_gas_src = conservative(2) the gas takes everything of a face except its
//! transport part Zt_f = om_f dt G_f (G_f = gr + off of the face kernel), so gas + radiation
//! (cell F) change by -(1/2) sum_f Zt_f e_f per cell, e_f = e_xi/sin | e_eta/sin the vector
//! a face-normal value adds to the cell F.  The radial faces carry no z-torque (F^1 lies
//! along rhat_c).  With a_cf = (V_c/2) (x_c x e_f)_z and the flux-form face torque T_c of
//! implicit_am_cs (sum_c T_c = the boundary torque), the cell's mismatch is
//!   mu_c = T_c + sum_{lat f of c} a_cf om_f G_f   (implicit_lz_mode = stress)
//! or, by default (full), with the face's whole Z_f/dt in place of om_f G_f: every part of
//! a face that the gas + radiation budget books (a seam face's midpoint theta leaves
//! (1 - om) F0^old, the blend its AP-HLL part) is then closed as well.
//! Every lateral face takes the extra stress h_f = -W_f (q_L + q_R), W_f = a_L,f + a_R,f,
//! q_c = mu_c / sum_{lat f of c} W_f^2, so sum_f W_f h_f = -sum_c mu_c exactly and
//!   dL_z(gas + rad) = -dt sum_c sum_f a_cf (om_f G_f + h_f) = dt sum_c T_c.
//! h enters the face flux as -th0 chat c dt h (th0 = 1/(1 + chat dt kt_f)) and Z_f as
//! +dt h (a stress: the gas pays its drag share only); it is lagged (TRHS) like the skew
//! cross term.  A laterally uniform stress has G = 0 and T_c = 0, so h = 0.
//! A seam neighbour's arm is the mirror image of the cell's own: a_B = (V/2) (R x_A x
//! (-R e_A))_z, R the reflection across the seam plane (normal n_f).

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {
// the arm a = (V/2) (x x e)_z of a cell (centre x = r rhat, e = e_t/sin), or, with refl,
// of its mirror image across the plane with unit normal n (x -> R x, e -> -R e)
KOKKOS_INLINE_FUNCTION
Real LzArm(const Real dv, const Real *x, const Real *e, const bool refl,
           const Real *n) {
  if (!refl) {return 0.5*dv*(x[0]*e[1] - x[1]*e[0]);}
  const Real xn = x[0]*n[0] + x[1]*n[1] + x[2]*n[2];
  const Real en = e[0]*n[0] + e[1]*n[1] + e[2]*n[2];
  Real xr[3], er[3];
  for (int q = 0; q < 3; ++q) {
    xr[q] = x[q] - 2.0*xn*n[q];
    er[q] = -(e[q] - 2.0*en*n[q]);
  }
  return 0.5*dv*(xr[0]*er[1] - xr[1]*er[0]);
}

// the cell-centre position x and the vector e = e_t/sin (dir 0: xi, 1: eta) of cell (k,j)
// of a panel MeshBlock at radius r
KOKKOS_INLINE_FUNCTION
void LzCellVec(const int pnl, const int k, const int j, const int ks, const int js,
               const int nx2, const int nx3, const Real x2l, const Real x2u,
               const Real x3l, const Real x3u, const Real r, const Real sn,
               const int dir, Real *x, Real *e) {
  const Real xc = 0.25*M_PI*CellCenterX(j-js, nx2, x2l, x2u);
  const Real ec = 0.25*M_PI*CellCenterX(k-ks, nx3, x3l, x3u);
  Real rh[3], e1[3], e2[3];
  cubed_sphere::PanelToCart(pnl, xc, ec, rh);
  cubed_sphere::PanelTangents(pnl, xc, ec, e1, e2);
  for (int q = 0; q < 3; ++q) {
    x[q] = r*rh[q];
    e[q] = ((dir == 0) ? e1[q] : e2[q])/sn;
  }
}

// tau_f = (M_f x n_f)_z of a planar lateral face (exact first moment, as m1_am_cs) with
// the edge points (xa, ea) and (xb, eb) and the radial moment mom = (r_r^3 - r_l^3)/3
KOKKOS_INLINE_FUNCTION
Real LzTau(const int pnl, const Real xa, const Real ea, const Real xb, const Real eb,
           const Real mom, const Real *nf) {
  Real qa[3], qb[3];
  cubed_sphere::PanelToCart(pnl, xa, ea, qa);
  cubed_sphere::PanelToCart(pnl, xb, eb, qb);
  const Real cx = qa[1]*qb[2] - qa[2]*qb[1];
  const Real cy = qa[2]*qb[0] - qa[0]*qb[2];
  const Real cz = qa[0]*qb[1] - qa[1]*qb[0];
  const Real dot = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2];
  const Real th = sqrt(cx*cx + cy*cy + cz*cz)/(1.0 + dot);
  const Real m0 = mom*th*(qa[0] + qb[0]), m1 = mom*th*(qa[1] + qb[1]);
  return m0*nf[1] - m1*nf[0];
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::LzGeom
//! \brief the geometry of implicit_lz_form: the cell arms a_c (x2, x3), and per lateral
//! face W_f = a_L + a_R (a seam neighbour by reflection) and tau_f.  Physical faces 0.

void RadiationM1::LzGeom() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2a = indcs.nx2, nx3a = indcs.nx3;
  const int nc1 = indcs.nx1 + 2*indcs.ng, nc2 = indcs.nx2 + 2*indcs.ng;
  const int nc3 = indcs.nx3 + 2*indcs.ng;
  const int nmb = pmy_pack->nmb_thispack, nmb1 = nmb - 1;
  Kokkos::realloc(lz_arm, nmb, 2, nc3, nc2, nc1);
  Kokkos::realloc(lz_w2, nmb, nc3, nc2+1, nc1);
  Kokkos::realloc(lz_t2, nmb, nc3, nc2+1, nc1);
  Kokkos::realloc(lz_g2, nmb, nc3, nc2+1, nc1);
  Kokkos::realloc(lz_o2, nmb, nc3, nc2+1, nc1);
  Kokkos::realloc(lz_w3, nmb, nc3+1, nc2, nc1);
  Kokkos::realloc(lz_t3, nmb, nc3+1, nc2, nc1);
  Kokkos::realloc(lz_g3, nmb, nc3+1, nc2, nc1);
  Kokkos::realloc(lz_o3, nmb, nc3+1, nc2, nc1);
  Kokkos::deep_copy(lz_arm, 0.0);
  Kokkos::deep_copy(lz_w2, 0.0);
  Kokkos::deep_copy(lz_t2, 0.0);
  Kokkos::deep_copy(lz_g2, 0.0);
  Kokkos::deep_copy(lz_o2, 0.0);
  Kokkos::deep_copy(lz_w3, 0.0);
  Kokkos::deep_copy(lz_t3, 0.0);
  Kokkos::deep_copy(lz_g3, 0.0);
  Kokkos::deep_copy(lz_o3, 0.0);
  auto vola = pmy_pack->pcoord->volume;
  auto x1v = pmy_pack->pcoord->x1v;
  auto x1f = pmy_pack->pcoord->xx1f;
  auto sinc = pmy_pack->pcoord->sin_cell;
  auto &msz = pmy_pack->pmb->mb_size;
  auto &mpn = pmy_pack->pmb->mb_panel;
  auto &mbbcs = m1bcs;
  auto fn2 = cs_fn2;
  auto fn3 = cs_fn3;
  auto cse = cs_seam.d_view;
  auto arm_ = lz_arm;
  auto w2_ = lz_w2, t2_ = lz_t2, w3_ = lz_w3, t3_ = lz_t3;
  par_for("m1_lz_arm", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x2l = msz.d_view(m).x2min, x2u = msz.d_view(m).x2max;
    const Real x3l = msz.d_view(m).x3min, x3u = msz.d_view(m).x3max;
    const int pnl = mpn.d_view(m);
    const Real nn[3] = {0.0, 0.0, 1.0};
    for (int d = 0; d < 2; ++d) {
      Real x[3], e[3];
      LzCellVec(pnl, k, j, ks, js, nx2a, nx3a, x2l, x2u, x3l, x3u, x1v(m,i),
                sinc(m,k,j), d, x, e);
      arm_(m,d,k,j,i) = LzArm(vola(m,k,j,i), x, e, false, nn);
    }
  });
  // the x2 faces (xi = const)
  par_for("m1_lz_w2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    const BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    const bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    const bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) return;
    const Real x2l = msz.d_view(m).x2min, x2u = msz.d_view(m).x2max;
    const Real x3l = msz.d_view(m).x3min, x3u = msz.d_view(m).x3max;
    const int pnl = mpn.d_view(m);
    const Real nf[3] = {fn2(m,j,0), fn2(m,j,1), fn2(m,j,2)};
    const bool slo = (j == js) && (cse(m,0) != 0), shi = (j == je+1) && (cse(m,1) != 0);
    const Real r = x1v(m,i);
    Real x[3], e[3], w = 0.0;
    if (slo || shi) {
      const int jin = slo ? j : (j - 1);
      LzCellVec(pnl, k, jin, ks, js, nx2a, nx3a, x2l, x2u, x3l, x3u, r, sinc(m,k,jin), 0,
                x, e);
      const Real dv = vola(m,k,jin,i);
      w = LzArm(dv, x, e, false, nf) + LzArm(dv, x, e, true, nf);
    } else {
      for (int s = 0; s < 2; ++s) {
        const int jq = j - 1 + s;
        LzCellVec(pnl, k, jq, ks, js, nx2a, nx3a, x2l, x2u, x3l, x3u, r, sinc(m,k,jq), 0,
                  x, e);
        w += LzArm(vola(m,k,jq,i), x, e, false, nf);
      }
    }
    w2_(m,k,j,i) = w;
    const Real rl = x1f(m,i), rr = x1f(m,i+1);
    const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
    const Real xa = 0.25*M_PI*LeftEdgeX(j-js, nx2a, x2l, x2u);
    const Real el = 0.25*M_PI*LeftEdgeX(k-ks, nx3a, x3l, x3u);
    const Real er = 0.25*M_PI*LeftEdgeX(k+1-ks, nx3a, x3l, x3u);
    t2_(m,k,j,i) = LzTau(pnl, xa, el, xa, er, mom, nf);
  });
  // the x3 faces (eta = const)
  par_for("m1_lz_w3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    const BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    const bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    const bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) return;
    const Real x2l = msz.d_view(m).x2min, x2u = msz.d_view(m).x2max;
    const Real x3l = msz.d_view(m).x3min, x3u = msz.d_view(m).x3max;
    const int pnl = mpn.d_view(m);
    const Real nf[3] = {fn3(m,k,0), fn3(m,k,1), fn3(m,k,2)};
    const bool slo = (k == ks) && (cse(m,2) != 0), shi = (k == ke+1) && (cse(m,3) != 0);
    const Real r = x1v(m,i);
    Real x[3], e[3], w = 0.0;
    if (slo || shi) {
      const int kin = slo ? k : (k - 1);
      LzCellVec(pnl, kin, j, ks, js, nx2a, nx3a, x2l, x2u, x3l, x3u, r, sinc(m,kin,j), 1,
                x, e);
      const Real dv = vola(m,kin,j,i);
      w = LzArm(dv, x, e, false, nf) + LzArm(dv, x, e, true, nf);
    } else {
      for (int s = 0; s < 2; ++s) {
        const int kq = k - 1 + s;
        LzCellVec(pnl, kq, j, ks, js, nx2a, nx3a, x2l, x2u, x3l, x3u, r, sinc(m,kq,j), 1,
                  x, e);
        w += LzArm(vola(m,kq,j,i), x, e, false, nf);
      }
    }
    w3_(m,k,j,i) = w;
    const Real rl = x1f(m,i), rr = x1f(m,i+1);
    const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
    const Real ea = 0.25*M_PI*LeftEdgeX(k-ks, nx3a, x3l, x3u);
    const Real xl = 0.25*M_PI*LeftEdgeX(j-js, nx2a, x2l, x2u);
    const Real xr = 0.25*M_PI*LeftEdgeX(j+1-js, nx2a, x2l, x2u);
    t3_(m,k,j,i) = LzTau(pnl, xl, ea, xr, ea, mom, nf);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::LzApply
//! \brief one implicit_lz_form pass (called by ImplicitTransverseTerms after its face
//! kernels): mu_c and q_c, the halo of q, then h_f into f0x2 / f0x3 and gcz2 / gcz3.

void RadiationM1::LzApply() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2a = indcs.nx2, nx3a = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  auto arm_ = lz_arm;
  auto w2_ = lz_w2, t2_ = lz_t2, g2_ = lz_g2, o2_ = lz_o2;
  auto w3_ = lz_w3, t3_ = lz_t3, g3_ = lz_g3, o3_ = lz_o3;
  auto &mbbcs = m1bcs;
  const int ilz = iw_lz;
  const Real dt = dt_sub, ch = chat, cl = c_light;
  auto d_ = lz_d;
  const bool dg = (impl_lz_diag > lz_ndiag);
  if (dg) {Kokkos::deep_copy(lz_d, 0.0);}
  // vet_gd D_r,lat: the radial faces' tangential stress (as m1_am_cs)
  auto tt_ = VlatTT();
  const int c0t = M1_TT_LAT0;
  const bool vlc = vlat_now && vlat_ready;
  const bool apr = (impl_flux != M1_IFLUX_CENTRAL);
  auto ifr_ = ifw;
  auto a1a = pmy_pack->pcoord->area.x1f;
  auto x1fa = pmy_pack->pcoord->xx1f;
  auto &msz = pmy_pack->pmb->mb_size;
  auto &mpn = pmy_pack->pmb->mb_panel;
  // implicit_lz_mode = full (default): mu from the face's whole Z_f (everything the gas +
  // radiation budget books: the stress om G, and the non-stress parts of a face, e.g. the
  // (1 - om) F0^old of a seam face's midpoint theta or of the blend); stress: om G only
  const bool lzfull = (impl_lz_mode == 0);
  auto gz2_ = gcz2;
  auto gz3_ = gcz3;
  const Real idt = 1.0/dt;
  par_for("m1_lz_mu", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real tz = 0.0, ag = 0.0, s2 = 0.0;
    const Real a2 = arm_(m,0,k,j,i), a3 = arm_(m,1,k,j,i);
    for (int n = 0; n < 2; ++n) {
      const int jf = j + n;
      const Real w = w2_(m,k,jf,i);
      const Real pf = 0.5*(M1DDiag(iw_,vd_,dfull,m,1,k,jf-1,i)*iw_(m,M1_IW_EP,k,jf-1,i)
                           + M1DDiag(iw_,vd_,dfull,m,1,k,jf,i)*iw_(m,M1_IW_EP,k,jf,i));
      const Real t = o2_(m,k,jf,i)*t2_(m,k,jf,i)*pf;
      tz += n ? -t : t;
      ag += a2*(lzfull ? gz2_(m,k,jf,i)*idt : g2_(m,k,jf,i));
      s2 += w*w;
    }
    for (int n = 0; n < 2; ++n) {
      const int kf = k + n;
      const Real w = w3_(m,kf,j,i);
      const Real pf = 0.5*(M1DDiag(iw_,vd_,dfull,m,2,kf-1,j,i)*iw_(m,M1_IW_EP,kf-1,j,i)
                           + M1DDiag(iw_,vd_,dfull,m,2,kf,j,i)*iw_(m,M1_IW_EP,kf,j,i));
      const Real t = o3_(m,kf,j,i)*t3_(m,kf,j,i)*pf;
      tz += n ? -t : t;
      ag += a3*(lzfull ? gz3_(m,kf,j,i)*idt : g3_(m,kf,j,i));
      s2 += w*w;
    }
    if (vlc) {
      const Real xc = 0.25*M_PI*CellCenterX(j-js, nx2a, msz.d_view(m).x2min,
                                            msz.d_view(m).x2max);
      const Real ec = 0.25*M_PI*CellCenterX(k-ks, nx3a, msz.d_view(m).x3min,
                                            msz.d_view(m).x3max);
      Real rh[3];
      cubed_sphere::PanelToCart(mpn.d_view(m), xc, ec, rh);
      for (int n = 0; n < 2; ++n) {
        const int ifc = i + n;
        const bool bnd = (n == 0) ? (i == is && mbbcs.d_view(m,BoundaryFace::inner_x1) !=
                                     BoundaryFlag::block)
                                  : (i == ie && mbbcs.d_view(m,BoundaryFace::outer_x1) !=
                                     BoundaryFlag::block);
        if (bnd) continue;
        const Real w0 = 0.5*(tt_(m,c0t+1,k,j,ifc-1)*iw_(m,M1_IW_EP,k,j,ifc-1)
                             + tt_(m,c0t+1,k,j,ifc)*iw_(m,M1_IW_EP,k,j,ifc));
        const Real w1 = 0.5*(tt_(m,c0t+2,k,j,ifc-1)*iw_(m,M1_IW_EP,k,j,ifc-1)
                             + tt_(m,c0t+2,k,j,ifc)*iw_(m,M1_IW_EP,k,j,ifc));
        const Real ar = a1a(m,k,j,ifc)*x1fa(m,ifc);
        Real t = ar*(rh[0]*w1 - rh[1]*w0);
        if (apr) {t *= 1.0 - ifr_(m,M1_IFW_AL,k,j,ifc);}
        tz += n ? -t : t;
      }
    }
    const Real mu = tz + ag;
    iw_(m,ilz,k,j,i) = (s2 > 0.0) ? mu/s2 : 0.0;
    if (dg) {
      Kokkos::atomic_add(&d_(0), dt*mu);
      Kokkos::atomic_add(&d_(1), fabs(dt*mu));
      Kokkos::atomic_add(&d_(3), dt*tz);
    }
  });
  ImplicitHaloExchange(1, ilz);
  auto f2_ = f0x2;
  auto f3_ = f0x3;
  par_for("m1_lz_h2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    const BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    const bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    const bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) return;
    const Real h = -w2_(m,k,j,i)*(iw_(m,ilz,k,j-1,i) + iw_(m,ilz,k,j,i));
    const Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
    f2_(m,k,j,i) -= ch*cl*dt*h/(1.0 + ch*dt*ktf);
    gz2_(m,k,j,i) += dt*h;
  });
  par_for("m1_lz_h3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    const BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    const bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    const bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) return;
    const Real h = -w3_(m,k,j,i)*(iw_(m,ilz,k-1,j,i) + iw_(m,ilz,k,j,i));
    const Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
    f3_(m,k,j,i) -= ch*cl*dt*h/(1.0 + ch*dt*ktf);
    gz3_(m,k,j,i) += dt*h;
  });
  if (!dg) return;
  // diagnostic: the realised cell-mean torque of h, sum_c sum_f a_cf dt h_f (= -sum dt mu)
  par_for("m1_lz_dg", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = 0.0;
    for (int n = 0; n < 2; ++n) {
      const int jf = j + n, kf = k + n;
      s += arm_(m,0,k,j,i)*(-w2_(m,k,jf,i)*(iw_(m,ilz,k,jf-1,i) + iw_(m,ilz,k,jf,i)));
      s += arm_(m,1,k,j,i)*(-w3_(m,kf,j,i)*(iw_(m,ilz,kf-1,j,i) + iw_(m,ilz,kf,j,i)));
    }
    Kokkos::atomic_add(&d_(2), dt*s);
    // the write-back's attributed lateral change, -sum_f a_cf Z_f (= dt T_c summed)
    Real z = 0.0;
    for (int n = 0; n < 2; ++n) {
      z -= arm_(m,0,k,j,i)*gz2_(m,k,j+n,i) + arm_(m,1,k,j,i)*gz3_(m,k+n,j,i);
    }
    Kokkos::atomic_add(&d_(4), z);
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), lz_d);
  Real v[5] = {h(0), h(1), h(2), h(3), h(4)};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, v, 5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    std::cout << std::scientific << std::setprecision(6) << "LZ_FORM pass " << lz_ndiag
              << " sum_dt_mu " << v[0] << " sum|dt_mu| " << v[1] << " sum_dt_a.h " << v[2]
              << " sum_dt_T " << v[3] << " (mu + a.h) " << (v[0] + v[2]) << " -sum_aZ " << v[4]
              << std::endl;
  }
  lz_ndiag += 1;
}

}  // namespace radm1
