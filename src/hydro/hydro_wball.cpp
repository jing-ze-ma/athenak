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

//----------------------------------------------------------------------------------------
//! \fn void Hydro::BuildWBCacheAll
//! \brief the x2 and x3 local hydrostatic backgrounds of every cell whose 3-point stencil
//! lies inside the arrays (wbq2 / wbq3, the layout of wbq0) and the weights wbw(0..2)
//! (wbw(0) = 1: x1 stays the global scheme).  Built once per stage, after ConsToPrim.

void Hydro::BuildWBCacheAll() {
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
  auto eos = peos->eos_data;
  const WBOption wbo = wb_option;
  const bool gen = eos.IsGeneral();
  auto &w0_ = w0;
  auto &phicc = phicc0;
  auto &wt = wtemp;
  auto &ww = wbw;
  const bool wgt = wb_weight;
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

void Hydro::AddWbGravity(const Real bdt, const bool with_x1) {
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
  if (pmy_pack->pmesh->use_cubed_sphere) {
    // the cubed sphere: the face-sum form on the covariant (e_xi, e_eta) components, the
    // angular faces carrying their own sin factor exactly as the flux conversion
    // (GnomonicEquiangleFluxX2/X3) and the coordinate source x_ov_rD / y_ov_rC take it:
    // S_a = [A+ s+ (p_bg,+ - p_c) - A- s- (p_bg,- - p_c)]/V (+ the linearised analytic
    // part with the same weights)
    if ((pam != nullptr) && pam->am_transport) {
      std::cout << "### FATAL ERROR in hydro_wball: <hydro>/wb_dirs = all on the cubed "
                << "sphere with am_transport/am_vector: not yet" << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
    // <hydro>/am_transport (the scalar zonal correction, hydro_am.hpp) turns the x3 flux
    // force into -(Rbar_3/R_c)(F+ A+ - F- A-)/V (Rbar_3 the phi face's mean lever arm, R_c
    // = x1v sin(x2v)): the x3 source takes the same factor, so the pair stays balanced
    const bool amx3 = (pam != nullptr) && pam->am_transport;
    if (amx3 && pam->am_vector) {
      std::cout << "### FATAL ERROR in hydro_wball: <hydro>/wb_dirs = all on spherical "
                << "polar with am_vector: not yet (am_transport scalar only)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto x1v = pmy_pack->pcoord->x1v, x1f = pmy_pack->pcoord->xx1f;
    auto x2v = pmy_pack->pcoord->x2v, x2f = pmy_pack->pcoord->xx2f;
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
        Real fam = 1.0;
        if (amx3) {
          const Real rl = x1f(m,i), rr = x1f(m,i+1);
          const Real thl = x2f(m,j), thr = x2f(m,j+1);
          const Real rc = x1v(m,i)*sin(x2v(m,j));
          const Real r3 = (rr*rr*rr - rl*rl*rl)/3.0, r2h = 0.5*(rr*rr - rl*rl);
          const Real rb3 = (thr > thl) ? r3*(cos(thl) - cos(thr))/(r2h*(thr - thl)) : 0.0;
          if (rc > 0.0) fam = rb3/rc;
        }
        // the WB part with the am_transport factor, the analytic part plain
        const Real sp = fam*wp*(c3(m,13,k,j,i) - pc) + (1.0 - wp)*(-d*(f3(m,k+1,j,i) - phc));
        const Real sm = fam*wm*(c3(m,11,k,j,i) - pc) + (1.0 - wm)*(-d*(f3(m,k,j,i) - phc));
        u0_(m,IM3,k,j,i) += bdt*iv*(a3(m,k+1,j,i)*sp - a3(m,k,j,i)*sm);
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

}  // namespace hydro
