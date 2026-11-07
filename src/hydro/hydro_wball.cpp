//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_wball.cpp
//! \brief wball-1104: <hydro>/wb_dirs = all -- the x2/x3 background cache with the
//! validity weight (BuildWBCacheAll) and the shared well-balanced gravity source
//! (AddWbGravity).  See hydro_wball.hpp for the scheme.  Cartesian geometry in this
//! milestone (hydro.cpp refuses the spherical grids).

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "hydro.hpp"
#include "hydro_wball.hpp"
#include "hydro_am.hpp"

namespace hydro {

#define WBALL_ALIASES \
  MeshBlockPack *pmy_pack = c.pmbp; (void) pmy_pack; \
  auto &w0 = c.w0; auto &u0 = c.u0; auto &wder = c.wder; (void) wder; \
  auto &wbq0 = c.wbq0; auto &wbq2 = c.wbq2; auto &wbq3 = c.wbq3; auto &wbw = c.wbw; \
  auto &phicc0 = c.phicc; auto &wtemp = c.wtemp; auto &pwb = c.pwb; (void) pwb; \
  struct WbPhiF {DvceArray4D<Real> x1f, x2f, x3f;}; WbPhiF phi0{c.phi1f, c.phi2f, c.phi3f}; AmModule *pam = c.pam; (void) pam; \
  const WBOption wb_option = c.wbo; (void) wb_option; \
  const bool wb_weight = c.wgt, wb_guard_ideal = c.wbgi; \
  const bool use_wellbalance_static = c.wbs; (void) use_wellbalance_static; \
  const Real wbw_r0 = c.r0, wbw_r1 = c.r1, wbw_m0 = c.m0, wbw_m1 = c.m1, wbw_eta = c.eta; \
  (void) wb_weight; (void) wb_guard_ideal; (void) wbw_r0; (void) wbw_r1; (void) wbw_m0; \
  (void) wbw_m1; (void) wbw_eta; (void) wbq0; (void) wbq2; (void) wbq3; (void) wbw; \
  (void) phicc0; (void) wtemp; (void) phi0; (void) u0; (void) w0


//----------------------------------------------------------------------------------------
//! \fn void Hydro::BuildWBCacheAll
//! \brief the x2 and x3 local hydrostatic backgrounds of every cell whose 3-point stencil
//! lies inside the arrays (wbq2 / wbq3, the layout of wbq0) and the weights wbw(0..2)
//! (wbw(0) = 1: x1 stays the global scheme).  Built once per stage, after ConsToPrim.

void WbAllBuild(WbAllCtx &c) {
  WBALL_ALIASES;
  if (pmy_pack->pmesh->use_spherical_polar &&
      (pmy_pack->pcoord->sp_cart_polar_momentum || pmy_pack->pcoord->sp_wellbalanced_src)) {
    std::cout << "### FATAL ERROR in hydro_wball: <hydro>/wb_dirs = all on spherical "
              << "polar needs the default sp momentum/source (no sp_cart_polar_momentum, "
              << "sp_cart_all_momentum, sp_wellbalanced_src)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto eos = c.eos;
  const WBOption wbo = wb_option;
  const bool gen = eos.IsGeneral();
  auto &w0_ = w0;
  auto &phicc = phicc0;
  auto &wt = wtemp;
  auto &ww = wbw;
  const bool wgt = wb_weight;
  const Real phimx = c.phimax;
  auto pw = c.phiwb;
  const Real r0 = wbw_r0, r1 = wbw_r1, m0 = wbw_m0, m1 = wbw_m1, eta = wbw_eta;
  const bool wbgi = wb_guard_ideal;
  const Real gam = eos.gamma, gm1 = eos.gamma - 1.0;
  Kokkos::deep_copy(Kokkos::subview(ww, Kokkos::ALL, 0, Kokkos::ALL, Kokkos::ALL,
                                    Kokkos::ALL), 1.0);
  for (int dir = 2; dir <= 3; ++dir) {
    if (dir == 2 && n2 == 1) continue;
    if (dir == 3 && n3 == 1) continue;
    auto c = (dir == 2) ? wbq2 : wbq3;
    auto phif = (dir == 2) ? phi0.x2f : phi0.x3f;
    const int dk = (dir == 3) ? 1 : 0, dj = (dir == 2) ? 1 : 0;
    const int kl = (n3 > 1) ? dk : 0, ku = (n3 > 1) ? n3 - 1 - dk : 0;
    const int jl = dj, ju = n2 - 1 - dj;
    par_for("wball_cache", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int km = k - dk, kp = k + dk, jm = j - dj, jp = j + dj;
      const Real dm = w0_(m,IDN,km,jm,i), dc = w0_(m,IDN,k,j,i), dp = w0_(m,IDN,kp,jp,i);
      const Real em = w0_(m,IEN,km,jm,i), ec = w0_(m,IEN,k,j,i), ep = w0_(m,IEN,kp,jp,i);
      const Real fm = phicc(m,km,jm,i), fc = phicc(m,k,j,i), fp = phicc(m,kp,jp,i);
      const Real fmh = phif(m,k,j,i), fph = phif(m,kp,jp,i);
      if (gen) {
        WBState s_m, s_mh, s_c, s_ph, s_p;
        WBBackgroundStencil(eos, wbo, dm, dc, dp, em, ec, ep, fm, fmh, fc, fph, fp,
                            s_m, s_mh, s_c, s_ph, s_p, WBT(wt,m,km,jm,i), WBT(wt,m,k,j,i),
                            WBT(wt,m,kp,jp,i));
        c(m,0,k,j,i) = s_m.d; c(m,1,k,j,i) = s_mh.d; c(m,2,k,j,i) = s_c.d;
        c(m,3,k,j,i) = s_ph.d; c(m,4,k,j,i) = s_p.d;
        c(m,5,k,j,i) = s_m.e; c(m,6,k,j,i) = s_mh.e; c(m,7,k,j,i) = s_c.e;
        c(m,8,k,j,i) = s_ph.e; c(m,9,k,j,i) = s_p.e;
        c(m,10,k,j,i) = s_m.p; c(m,11,k,j,i) = s_mh.p; c(m,12,k,j,i) = s_c.p;
        c(m,13,k,j,i) = s_ph.p; c(m,14,k,j,i) = s_p.p;
      } else {
        for (int var=0; var<2; ++var) {
          Real a, b, cc, d, e;
          Hydro::getWBq0(eos, wbo, var, dm, dc, dp, em, ec, ep, fm, fmh, fc, fph, fp,
                         a, b, cc, d, e, -1.0, -1.0, -1.0, false);
          c(m,5*var,k,j,i) = a; c(m,5*var+1,k,j,i) = b; c(m,5*var+2,k,j,i) = cc;
          c(m,5*var+3,k,j,i) = d; c(m,5*var+4,k,j,i) = e;
        }
        for (int q=0; q<5; ++q) c(m,10+q,k,j,i) = gm1*c(m,5+q,k,j,i);
        if (wbgi) {
          WBState s[5];
          for (int q=0; q<5; ++q) {
            s[q].d = c(m,q,k,j,i); s[q].e = c(m,5+q,k,j,i); s[q].p = c(m,10+q,k,j,i);
          }
          WBGuard(s[0], s[1], s[2], s[3], s[4]);
          for (int q=0; q<5; ++q) {
            c(m,q,k,j,i) = s[q].d; c(m,5+q,k,j,i) = s[q].e; c(m,10+q,k,j,i) = s[q].p;
          }
        }
      }
      // the validity weight of this cell in this direction
      Real w = 1.0;
      if (wgt) {
        Real pm, pc, pp, cs2;
        if (gen) {
          pm = eos.Pressure(dm, em, WBT(wt,m,km,jm,i));
          pc = eos.Pressure(dc, ec, WBT(wt,m,k,j,i));
          pp = eos.Pressure(dp, ep, WBT(wt,m,kp,jp,i));
          cs2 = gam*pc/dc;   // the Mach ramp only needs the order of magnitude
        } else {
          pm = gm1*em; pc = gm1*ec; pp = gm1*ep; cs2 = gam*pc/dc;
        }
        const Real dpr = pp - pm, dph = fp - fm;
        const Real rd = fabs(dpr + dc*dph)/(dc*fabs(dph) + fabs(dpr) + eta*pc + 1.0e-300);
        const Real v2 = w0_(m,IVX,k,j,i)*w0_(m,IVX,k,j,i) +
                        w0_(m,IVY,k,j,i)*w0_(m,IVY,k,j,i) +
                        w0_(m,IVZ,k,j,i)*w0_(m,IVZ,k,j,i);
        const Real mach = sqrt(v2/fmax(cs2, 1.0e-300));
        w = WbRamp(rd, r0, r1)*WbRamp(mach, m0, m1);
      }
      // dhjwball-1106: no WB outside wb_phimax (dhj Roche: the L1 equipotential)
      // (the wb_phimax scale = phiwb; x1 neighbours too: the AM background reads their wbq0)
      if (phimx < 1.0e299) {
        const int ia = (i > 0) ? i - 1 : 0, ib = (i < n1 - 1) ? i + 1 : n1 - 1;
        if (fmax(fmax(pw(m,km,jm,i), pw(m,k,j,i)), pw(m,kp,jp,i)) > phimx ||
            fmax(pw(m,k,j,ia), pw(m,k,j,ib)) > phimx) w = 0.0;
      }
      ww(m,dir-1,k,j,i) = w;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::AddWbGravity
//! \brief the shared well-balanced gravitational MOMENTUM source of <hydro>/wb_dirs =
//! all, the face-sum form of hydro_wball.hpp (Cartesian: A/V = 1/dx).  x2, x3 always;
//! x1 (wbq0, w = 1, the pgen's own fg form) when with_x1.  The energy is the caller's
//! (etotgrav, or its own v.(-rho grad Phi)).  Reads the caches of THIS stage: call it
//! from the user source of the stage whose fluxes built them.

void WbAllAmBackground(WbAllCtx &c);
//! dhjwball-1106 <fluid>/wb_poison_test (TEST ONLY): every cached WB background of a cell
//! with Phi_wb > wb_phimax -> a huge value: if any code path still reads it, the run differs
//! from the run without the key.  nv channels of cache(m,n,k,j,i).
void WbPoison(MeshBlockPack *pmbp, DvceArray5D<Real> &cache, const DvceArray4D<Real> &phiwb,
              const Real phimax) {
  auto q = cache; auto pw = phiwb;
  const int nm = pmbp->nmb_thispack - 1, nv = q.extent_int(1) - 1;
  const int n3 = q.extent_int(2) - 1, n2 = q.extent_int(3) - 1, n1 = q.extent_int(4) - 1;
  par_for("wb_poison", DevExeSpace(), 0, nm, 0, nv, 0, n3, 0, n2, 0, n1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    if (pw(m,k,j,i) > phimax) q(m,n,k,j,i) = 1.0e40;
  });
}
void WbAllGravity(WbAllCtx &c, const Real bdt, const bool with_x1) {
  WBALL_ALIASES;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool x2 = (indcs.nx2 > 1), x3 = (indcs.nx3 > 1);
  auto u0_ = u0;
  auto w0_ = w0;
  auto c1 = wbq0, c2 = wbq2, c3 = wbq3;
  auto ww = wbw;
  auto phicc = phicc0;
  auto f2 = phi0.x2f, f3 = phi0.x3f;
  auto &size = pmy_pack->pmb->mb_size;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const bool am = (pam != nullptr) && pam->am_transport &&
                  (pmy_pack->pmesh->use_cubed_sphere || pmy_pack->pmesh->use_spherical_polar);
  if (am) WbAllAmBackground(c);
  auto amu = *c.wbam_u;
  if (pmy_pack->pmesh->use_cubed_sphere) {
    // the cubed sphere: the face-sum form on the covariant (e_xi, e_eta) components, the
    // angular faces carrying their own sin factor exactly as the flux conversion
    // (GnomonicEquiangleFluxX2/X3) and the coordinate source x_ov_rD / y_ov_rC take it:
    // S_a = [A+ s+ (p_bg,+ - p_c) - A- s- (p_bg,- - p_c)]/V (+ the linearised analytic
    // part with the same weights)
    if (pmy_pack->pcoord->cs_wellbalanced_src) {
      std::cout << "### FATAL ERROR in hydro_wball: <hydro>/wb_dirs = all on the cubed "
                << "sphere needs cs_wellbalanced_src = false" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto vol = pmy_pack->pcoord->volume;
    auto a1 = pmy_pack->pcoord->area.x1f, a2 = pmy_pack->pcoord->area.x2f;
    auto a3 = pmy_pack->pcoord->area.x3f;
    auto sx = pmy_pack->pcoord->sin_face_xi, se = pmy_pack->pcoord->sin_face_eta;
    par_for("wball_grav_cs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = w0_(m,IDN,k,j,i);
      const Real phc = phicc(m,k,j,i);
      const Real iv = 1.0/vol(m,k,j,i);
      if (with_x1) {
        const Real pc = c1(m,12,k,j,i);
        u0_(m,IM1,k,j,i) += bdt*iv*(a1(m,k,j,i+1)*(c1(m,13,k,j,i) - pc) +
                                    a1(m,k,j,i)*(pc - c1(m,11,k,j,i)));
      }
      {
        const Real pc = c2(m,12,k,j,i);
        const Real wc = ww(m,1,k,j,i);
        const int jp = (j + 1 > n2 - 1) ? n2 - 1 : j + 1;
        const Real wp = fmin(wc, ww(m,1,k,jp,i)), wm = fmin(wc, ww(m,1,k,j-1,i));
        const Real sp = wp*(c2(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f2(m,k,j+1,i) - phc));
        const Real sm = wm*(c2(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f2(m,k,j,i) - phc));
        u0_(m,IM2,k,j,i) += bdt*iv*(a2(m,k,j+1,i)*sx(m,k,j+1)*sp - a2(m,k,j,i)*sx(m,k,j)*sm);
      }
      {
        const Real pc = c3(m,12,k,j,i);
        const Real wc = ww(m,2,k,j,i);
        const int kp = (k + 1 > n3 - 1) ? n3 - 1 : k + 1;
        const Real wp = fmin(wc, ww(m,2,kp,j,i)), wm = fmin(wc, ww(m,2,k-1,j,i));
        const Real sp = wp*(c3(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f3(m,k+1,j,i) - phc));
        const Real sm = wm*(c3(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f3(m,k,j,i) - phc));
        u0_(m,IM3,k,j,i) += bdt*iv*(a3(m,k+1,j,i)*se(m,k+1,j)*sp - a3(m,k,j,i)*se(m,k,j)*sm);
      }
      if (am) {   // minus the AM correction of the background fluxes, x the cell weight
        const Real wa = fmin(ww(m,1,k,j,i), ww(m,2,k,j,i));
        u0_(m,IM2,k,j,i) -= bdt*wa*amu(m,IM2,k,j,i);
        u0_(m,IM3,k,j,i) -= bdt*wa*amu(m,IM3,k,j,i);
      }
    });
    return;
  }
  if (pmy_pack->pmesh->use_spherical_polar) {
    // the face-sum form with the sp areas: S_d = [A+ s+ - A- s-]/V, s = w (p_bg,f - p_c)
    // + (1 - w)(-rho (Phi_f - Phi_c)); the sp coordinate source adds p_c (A+ - A-)/V
    // (x_ov_rD, z_ov_rE), so a background state balances to round-off
    auto vol = pmy_pack->pcoord->volume;
    auto a1 = pmy_pack->pcoord->area.x1f, a2 = pmy_pack->pcoord->area.x2f;
    auto a3 = pmy_pack->pcoord->area.x3f;
    par_for("wball_grav_sp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = w0_(m,IDN,k,j,i);
      const Real phc = phicc(m,k,j,i);
      const Real iv = 1.0/vol(m,k,j,i);
      if (with_x1) {
        const Real pc = c1(m,12,k,j,i);
        u0_(m,IM1,k,j,i) += bdt*iv*(a1(m,k,j,i+1)*(c1(m,13,k,j,i) - pc) +
                                    a1(m,k,j,i)*(pc - c1(m,11,k,j,i)));
      }
      if (x2) {
        const Real pc = c2(m,12,k,j,i);
        const Real wc = ww(m,1,k,j,i);
        const int jp = (j + 1 > n2 - 1) ? n2 - 1 : j + 1;
        const Real wp = fmin(wc, ww(m,1,k,jp,i)), wm = fmin(wc, ww(m,1,k,j-1,i));
        const Real sp = wp*(c2(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f2(m,k,j+1,i) - phc));
        const Real sm = wm*(c2(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f2(m,k,j,i) - phc));
        u0_(m,IM2,k,j,i) += bdt*iv*(a2(m,k,j+1,i)*sp - a2(m,k,j,i)*sm);
      }
      if (x3) {
        const Real pc = c3(m,12,k,j,i);
        const Real wc = ww(m,2,k,j,i);
        const int kp = (k + 1 > n3 - 1) ? n3 - 1 : k + 1;
        const Real wp = fmin(wc, ww(m,2,kp,j,i)), wm = fmin(wc, ww(m,2,k-1,j,i));
        const Real sp = wp*(c3(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f3(m,k+1,j,i) - phc));
        const Real sm = wm*(c3(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f3(m,k,j,i) - phc));
        u0_(m,IM3,k,j,i) += bdt*iv*(a3(m,k+1,j,i)*sp - a3(m,k,j,i)*sm);
      }
      if (am) {
        const Real wa = x3 ? fmin(ww(m,1,k,j,i), ww(m,2,k,j,i)) : ww(m,1,k,j,i);
        u0_(m,IM2,k,j,i) -= bdt*wa*amu(m,IM2,k,j,i);
        u0_(m,IM3,k,j,i) -= bdt*wa*amu(m,IM3,k,j,i);
      }
    });
    return;
  }
  par_for("wball_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0_(m,IDN,k,j,i);
    const Real phc = phicc(m,k,j,i);
    if (with_x1) {
      // x1: the global scheme, (p_bg,i+1/2 - p_bg,i-1/2)/dx1
      u0_(m,IM1,k,j,i) += bdt*(c1(m,13,k,j,i) - c1(m,11,k,j,i))/size.d_view(m).dx1;
    }
    if (x2) {
      const Real pc = c2(m,12,k,j,i);
      const Real wc = ww(m,1,k,j,i);
      const int jp = (j + 1 > n2 - 1) ? n2 - 1 : j + 1;
      const Real wp = fmin(wc, ww(m,1,k,jp,i)), wm = fmin(wc, ww(m,1,k,j-1,i));
      const Real sp = wp*(c2(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f2(m,k,j+1,i) - phc));
      const Real sm = wm*(c2(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f2(m,k,j,i) - phc));
      u0_(m,IM2,k,j,i) += bdt*(sp - sm)/size.d_view(m).dx2;
    }
    if (x3) {
      const Real pc = c3(m,12,k,j,i);
      const Real wc = ww(m,2,k,j,i);
      const int kp = (k + 1 > n3 - 1) ? n3 - 1 : k + 1;
      const Real wp = fmin(wc, ww(m,2,kp,j,i)), wm = fmin(wc, ww(m,2,k-1,j,i));
      const Real sp = wp*(c3(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f3(m,k+1,j,i) - phc));
      const Real sm = wm*(c3(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f3(m,k,j,i) - phc));
      u0_(m,IM3,k,j,i) += bdt*(sp - sm)/size.d_view(m).dx3;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::WbAmBackground
//! \brief <hydro>/am_transport (or am_vector) corrects the zonal momentum update of the
//! ACTUAL fluxes by T_cons - T_code (hydro_am.hpp).  For a background state the fluxes are
//! the background face pressures, and that correction is a part of the discrete pressure
//! force the WB source must cancel.  Here: the background momentum fluxes of every face
//! (the mean of the two cells' walked face pressures; cs: on the covariant face basis,
//! p sin on the normal component), a static copy of w0 (v = 0), and AmTransportProbe on
//! them (bdt = 1) into wbam_u.  AddWbGravity subtracts it (x the cell weight).

void WbAllAmBackground(WbAllCtx &c) {
  WBALL_ALIASES;
  DvceFaceFld5D<Real> *&wbam_f = *c.wbam_f;
  DvceArray5D<Real> &wbam_u = *c.wbam_u, &wbam_w = *c.wbam_w;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool x3 = (indcs.nx3 > 1);
  if (wbam_f == nullptr) {
    const int n1 = indcs.nx1 + 2*indcs.ng;
    const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
    const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
    const int nv = u0.extent_int(1);
    wbam_f = new DvceFaceFld5D<Real>("wbam_f", nmb1 + 1, nv, n3, n2, n1);
    wbam_u = DvceArray5D<Real>("wbam_u", nmb1 + 1, nv, n3, n2, n1);
    wbam_w = DvceArray5D<Real>("wbam_w", nmb1 + 1, w0.extent_int(1), n3, n2, n1);
  }
  Kokkos::deep_copy(wbam_u, 0.0);
  Kokkos::deep_copy(wbam_f->x1f, 0.0);
  Kokkos::deep_copy(wbam_f->x2f, 0.0);
  if (x3) Kokkos::deep_copy(wbam_f->x3f, 0.0);
  Kokkos::deep_copy(wbam_w, w0);
  auto ws = wbam_w;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  par_for("wbam_static", DevExeSpace(), 0, nmb1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    ws(m,IVX,k,j,i) = 0.0; ws(m,IVY,k,j,i) = 0.0; ws(m,IVZ,k,j,i) = 0.0;
  });
  auto c1 = wbq0, c2 = wbq2, c3 = wbq3;
  auto f1 = wbam_f->x1f, f2 = wbam_f->x2f, f3 = wbam_f->x3f;
  auto sx = pmy_pack->pcoord->sin_face_xi, se = pmy_pack->pcoord->sin_face_eta;
  par_for("wbam_f1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f1(m,IM1,k,j,i) = 0.5*(c1(m,13,k,j,i-1) + c1(m,11,k,j,i));
  });
  par_for("wbam_f2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real p = 0.5*(c2(m,13,k,j-1,i) + c2(m,11,k,j,i));
    f2(m,IM2,k,j,i) = cs ? p*sx(m,k,j) : p;
  });
  if (x3) {
    par_for("wbam_f3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real p = 0.5*(c3(m,13,k-1,j,i) + c3(m,11,k,j,i));
      f3(m,IM3,k,j,i) = cs ? p*se(m,k,j) : p;
    });
  }
  AmFluid fl{wbam_u, wbam_u, wbam_w, wder, DvceArray5D<Real>(), *wbam_f, c.eos,
             false, use_wellbalance_static, pwb};
  pam->AmTransportProbe(fl);
}

//----------------------------------------------------------------------------------------
//! \fn void WbAllSeamFaces
//! \brief the cs shared-seam face states (Hydro/MHD::CsSeamFaceStates) under wb_dirs = all:
//! each edge cell's own face state is the one the sweep forms (WbAllPLM: the WB deviation
//! blended with the plain PLM already in fs by the face weight), density and energy (and
//! the derived pressure of a general EOS, record slot nvh + IDPR).  Rows ia..ib.

void WbAllSeamFaces(WbAllCtx &c, DvceFaceFld5D<Real> &fs, const int nvh, const int nder,
                    const bool wrho, const int ia, const int ib) {
  WBALL_ALIASES;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto f2 = fs.x2f, f3 = fs.x3f;
  auto q_ = w0, qd_ = wder;
  auto c2 = wbq2, c3 = wbq3, ww = wbw;
  auto &size = pmy_pack->pmb->mb_size;
  par_for("cs_seamfs_wb2", DevExeSpace(), 0, nmb1, ks, ke, 0, 1, ia, ib,
  KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x2min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x2max >= 1.0 - 1.0e-12)) return;
    const int jc = lo ? js : je, jf = lo ? js : je+1, jn = lo ? jc - 1 : jc + 1;
    const Real w = fmin(ww(m,1,k,jc,i), ww(m,1,k,jn,i));
    for (int n = 0; n < 2; ++n) {
      const int nv = (n == 0) ? IDN : IEN;
      if (nv == IDN && !wrho) continue;
      f2(m,nv,k,jf,i) = WbSeamBlend(q_(m,nv,k,jc-1,i), q_(m,nv,k,jc,i), q_(m,nv,k,jc+1,i),
                                    c2, (nv == IEN) ? WBVar::wb_eint : WBVar::wb_dens,
                                    m, k, jc, i, lo, w, f2(m,nv,k,jf,i));
    }
    if (nder > 0) {
      f2(m,nvh+IDPR,k,jf,i) = WbSeamBlend(qd_(m,IDPR,k,jc-1,i), qd_(m,IDPR,k,jc,i),
                                          qd_(m,IDPR,k,jc+1,i), c2, WBVar::wb_pres,
                                          m, k, jc, i, lo, w, f2(m,nvh+IDPR,k,jf,i));
    }
  });
  par_for("cs_seamfs_wb3", DevExeSpace(), 0, nmb1, 0, 1, js, je, ia, ib,
  KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
    const bool lo = (sd == 0);
    if (lo && !(size.d_view(m).x3min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x3max >= 1.0 - 1.0e-12)) return;
    const int kc = lo ? ks : ke, kf = lo ? ks : ke+1, kn = lo ? kc - 1 : kc + 1;
    const Real w = fmin(ww(m,2,kc,j,i), ww(m,2,kn,j,i));
    for (int n = 0; n < 2; ++n) {
      const int nv = (n == 0) ? IDN : IEN;
      if (nv == IDN && !wrho) continue;
      f3(m,nv,kf,j,i) = WbSeamBlend(q_(m,nv,kc-1,j,i), q_(m,nv,kc,j,i), q_(m,nv,kc+1,j,i),
                                    c3, (nv == IEN) ? WBVar::wb_eint : WBVar::wb_dens,
                                    m, kc, j, i, lo, w, f3(m,nv,kf,j,i));
    }
    if (nder > 0) {
      f3(m,nvh+IDPR,kf,j,i) = WbSeamBlend(qd_(m,IDPR,kc-1,j,i), qd_(m,IDPR,kc,j,i),
                                          qd_(m,IDPR,kc+1,j,i), c3, WBVar::wb_pres,
                                          m, kc, j, i, lo, w, f3(m,nvh+IDPR,kf,j,i));
    }
  });
}

//----------------------------------------------------------------------------------------
//! Hydro's members: fill the context, call the shared routines
WbAllCtx Hydro::WbCtx() {
  WbAllCtx c;
  c.pmbp = pmy_pack; c.eos = peos->eos_data; c.wbo = wb_option;
  c.wgt = wb_weight; c.wbgi = wb_guard_ideal; c.wbs = use_wellbalance_static;
  c.r0 = wbw_r0; c.r1 = wbw_r1; c.m0 = wbw_m0; c.m1 = wbw_m1; c.eta = wbw_eta;
  c.w0 = w0; c.u0 = u0; c.wder = wder; c.wbq0 = wbq0; c.wbq2 = wbq2; c.wbq3 = wbq3;
  c.wbw = wbw; c.phicc = phicc0; c.wtemp = wtemp; c.pwb = pwb; c.phi1f = phi0.x1f; c.phi2f = phi0.x2f; c.phi3f = phi0.x3f; c.pam = pam;
  c.wbam_u = &wbam_u; c.wbam_w = &wbam_w; c.wbam_f = &wbam_f;
  c.phimax = wb_phimax; c.phiwb = phicc_wb;
  return c;
}
void Hydro::BuildWBCacheAll() {WbAllCtx c = WbCtx(); WbAllBuild(c);}
void Hydro::AddWbGravity(const Real bdt, const bool with_x1) {
  WbAllCtx c = WbCtx(); WbAllGravity(c, bdt, with_x1);
}

}  // namespace hydro
