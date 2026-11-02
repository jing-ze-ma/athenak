//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter_rt_src.cpp
//! \brief deep_hot_jupiter_rt problem generator (tusplit-1102 piece): the source terms --
//! AmFix, the roche_wb_tan probe, SourceFunc, the double-grey two-stream RT and the
//! operator-split correlated-k RT (DhjCkRtSplit, DhjRtFreezeSplit).  Moved verbatim
//! from deep_hot_jupiter_rt.cpp.

#include <algorithm>
#include <cstdio>

#include "deep_hot_jupiter_rt.hpp"
#include "driver/dt_retry.hpp"
#include "utils/imex_sdirk2.hpp"   // the tableau shared with rad_m1 hesdirk2


//----------------------------------------------------------------------------------------
//! \fn void AmFix
//! \brief problem/coriolis_am and problem/am_transport: add the zonal torque
//!     dT = [coriolis_am] (-Omega Q_c - T_cc) + [am_transport] (T_cons - T_code)
//! of every active cell to its momentum (AmOut::c2/c3; no radial or meridional part).
//! Reads this stage's w0, wder, bcc0 and uflx and adds to u0 with the stage's bdt, like
//! every other term of SourceFunc; T_cc is the torque that SourceFunc's rotating-frame
//! source applied from the same w0.

void AmFix(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const AmV v = AmMakeV(pm, amh::trans);
  const bool fix = amh::fix;
  const bool trans = amh::trans;
  const Real omega = v.omega;
  par_for("dhj_am_fix", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    AmOut o;
    AmCell(v, m, k, j, i, o);
    if (!(o.wc > 0.0)) return;
    Real dt = 0.0;
    if (fix) dt += -omega*o.q - o.tcc;
    if (trans) dt += o.tcons - (o.tc[0] + o.tc[1] + o.tc[2] + o.tc[3]);
    u0(m,IM2,k,j,i) += dt*o.c2*bdt;
    u0(m,IM3,k,j,i) += dt*o.c3*bdt;
  });
}

//----------------------------------------------------------------------------------------
//! \fn WbTanProbe (roche-wbtan-1101, problem/roche_wb_tan): wbtan_common.hpp Probe
//! (hydro or MHD, cubed sphere or spherical polar), with the pgen's own AmFix of W_eq as
//! the extra term -- its problem/am_transport part (T_cons - T_code: a transport
//! correction whose sum is the boundary torque) ONLY.  NOT the problem/coriolis_am part:
//! -Omega Q_c of W_eq is the Coriolis torque on the equilibrium's NUMERICAL mass fluxes
//! (the first stage's radial redistribution, sum = Omega dI/dt); frozen into g_t it was
//! a constant net z torque (Ttide -5.5e36 cs C8 WASP deck, -9.8e34 = -1.2e-3 Iw Omega
//! per rotation on the C32 3x deck; wbtan-torquefix-1101).  The live AmFix applies it to
//! the evolving state every stage, as without the key.

void WbTanAmFixExtra(Mesh *pm, DvceArray5D<Real> &us, DvceArray5D<Real> &u0s) {
  if (!amh::trans) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  Kokkos::deep_copy(u0s, 0.0);
  const bool sv_fix = amh::fix;
  amh::fix = false;
  AmFix(pm, 1.0);   // adds to the fluid's u0, which IS u0s here
  amh::fix = sv_fix;
  auto u0a = u0s;
  auto usa = us;
  par_for("wbtan_amfix", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    usa(m,IM2,k,j,i) += u0a(m,IM2,k,j,i);
    usa(m,IM3,k,j,i) += u0a(m,IM3,k,j,i);
  });
}

void WbTanProbe(Mesh *pm, const Real bdt) {
  wbtan::Probe(pm, dhj_wbt, WbTanAmFixExtra, "deep_hot_jupiter_rt");
  (void) bdt;
}

void SourceFunc(Mesh *pm, Real bdt) {
  if (dhj_wbt.on && !dhj_wbt.ready) WbTanProbe(pm, bdt);
  EBudStage(pm, bdt);
  DhjFluxStage(pm);
  if (amh::on) AmHstStage(pm);
  if (amh::fix || amh::trans) AmFix(pm, bdt);
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0, w0wb;
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    bool use_etotgrav = false;
    bool use_wellbalance_static = false;
    bool use_wellbalance_dynamic = false;
    DvceArray5D<Real> wbq0_;    // the well-balanced background cache (BuildWBCache)
    Real gamma;
    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      wbq0_ = pmbp->phydro->wbq0;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
      w0wb = pmbp->phydro->w0wb;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      wbq0_ = pmbp->pmhd->wbq0;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      w0wb = pmbp->pmhd->w0wb;
    }
    
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    
//    Real ap = 9.44e9;
//    Real omega = 2.06e-5;
//    Real grav_acc = -942.0;
//    Real Rgas = 4.593e7;
    
//    ParameterInput* pin;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    const bool tide = pm->pgen->hot_jupiter_param.stellar_tide;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    // rot_potential: the radial centrifugal force is carried by the potential ONLY
    // through
    // the well-balanced source (the plain gravity source knows gravity alone), so the
    // explicit radial term is dropped only when that scheme is on; with etotgrav the
    // energy flux carries the potential either way
    const bool rotpot = pm->pgen->hot_jupiter_param.rot_potential;
    // problem/sponge_top (default true), problem/sponge_bottom (default false since
    // 09-30): switch the two velocity sponges at the end of this kernel off; captured by
    // value.
    const bool sponge_top_ = pm->pgen->hot_jupiter_param.sponge_top;
    const bool sponge_top_rad_ = dhj_sponge_top_radial;
    const Real sp_plo_ = dhj_sp_top_plo_bar, sp_phi_ = dhj_sp_top_phi_bar;
    const Real sp_tau_ = dhj_sp_top_tau;
    const bool sponge_bottom_ = pm->pgen->hot_jupiter_param.sponge_bottom;
    const bool rotpot_src = rotpot && use_wellbalance_dynamic;
    // roche-dhj-1031 problem/roche: the TANGENTIAL tidal force -grad_t TidePot (the radial
    // part is in phicc0/phi0 -> the well-balanced source; the work in the etotgrav flux,
    // both required with roche).  rgms = 0 when off: not applied, bitwise.
    const Real rgms = pm->pgen->hot_jupiter_param.roche ?
                      pm->pgen->hot_jupiter_param.roche_gms : 0.0;
    const Real rra = pm->pgen->hot_jupiter_param.roche_a;
    // roche-wbtan-1101 problem/roche_wb_tan: the tangential centrifugal + tidal force is
    // rho g_t (WbTanProbe) in place of the analytic projections below; Coriolis and the
    // radial terms unchanged; no energy term (etotgrav books the potential work)
    const bool wbt_ = dhj_wbt.on;
    // relaxdrag-1101: the start-up Rayleigh drag switch and length (R/5, 1 by default)
    const bool relax_ = dhj_relax_drag;
    const Real rdf_ = dhj_relax_drag_rot/5.0;
    DvceArray5D<Real> wbtg_;
    if (wbt_) wbtg_ = *dhj_wbt.gt;

    // Radial grid stretch, copied out of the Mesh so the device lambdas below capture
    // plain values. Radii rebuilt here with CellCenterX/LeftEdgeX are UNSTRETCHED, so
    // each is passed through ApplyRStretch to match the Coordinates arrays.
    const bool str_r_ = pm->use_grid_stretch_r;
    const bool str_rp_ = pm->use_grid_stretch_r_poly;
    const Real fstr_r_ = pm->fStretchR;
    const Real rmin_ = pm->mesh_size.x1min;
    const Real rmax_ = pm->mesh_size.x1max;
    Real cpoly_[NSTRETCH_R_POLY];
    for (int n=0; n<NSTRETCH_R_POLY; ++n) cpoly_[n] = pm->fStretchRPoly[n];
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    
    Real time = pm->time;
    
    // problem/ck_impl_once: the radiation has left the RK stage (see
    // DhjCkRtSplit); nothing radiative is evaluated here.
    // problem/rt_test_freeze: the radiation is applied by DhjRtFreezeSplit instead
    // problem/ck_time_scheme = hesdirk2: the SDIRK stage solve of this Heun stage
    if (two_stream_rt::ck_t2_on) {
      DhjCkT2Stage(pm, bdt);
    } else if (!two_stream_rt::ck_impl_once && !rt_test_freeze) {
      picket_fence_two_stream_RT(pm, bdt);
    }

    // which fluid module is on, as VALUES: dereferencing the host pointer `pmbp`
    // inside a device lambda is illegal on a discrete GPU (see CLAUDE.md).
    const bool is_hydro_ = (pmbp->phydro != nullptr);
    const bool is_mhd_ = (pmbp->pmhd != nullptr);
    // problem/budget_dt: the whole kernel's energy and mass change per shell, weighted
    // by what survives the rk2 combination (0.5 dt/bdt: 0.5 in stage 1, 1 in stage 2)
    const bool eb_on = ebud::on && ebud::started;
    DvceArray2D<Real> ebacc;
    if (eb_on) ebacc = *ebud::acc;
    const Real eb_w = 0.5*pm->dt/bdt;
    par_for("usrsource", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real eb_e0 = u0(m,IEN,k,j,i);
        const Real eb_d0 = u0(m,IDN,k,j,i);
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          // ASSIGN the x1v/x2v/x3v declared above -- do NOT redeclare them.  These three
          // lines used to read `Real x1v = ...`, which SHADOWED the outer variables and
          // threw the values away at the closing brace, leaving the outer ones
          // uninitialised for every non-spherical-polar grid.  Nothing exercised that
          // path until the cubed sphere did, and the result was the whole domain going
          // NaN on the first source-term application.
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
          if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        }
        
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
          z = x1v-ap;
          r = x1v;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          z = x1v-ap;
          r = x1v;
        } else {
          lam = x3v*iap;
          phi = x2v*iap;
          z = x1v;
        }
        Real rho = w0(m,IDN,k,j,i);
        Real p = eos.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,i),p);
        
        Real area_r = area1(m,k,j,i+1);
        Real area_l = area1(m,k,j,i);
        Real vol = volume(m,k,j,i);
        
        // gravity
        const Real grav_r = GravAccAt(grav_acc, ap, r, grav_pmass);
        Real src = bdt*grav_r*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance_static) {
            src = bdt*grav_r*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_dynamic) {
          // The source is the background's own pressure difference across the cell, so it
          // has to come from the same entry point the reconstruction uses: getWBq0, which
          // dispatches to the ideal-gas closed forms or to the general-EOS background and
          // returns pressure directly, not an energy to be multiplied by (gamma-1).
          Real pl,pr,dum1,dum2,dum3;
          if (is_hydro_) {
              WBReadCache(wbq0_, WBVar::wb_pres, m, k, j, i, dum1, pl, dum2, pr, dum3);
          } else if (is_mhd_) {
              WBReadCache(wbq0_, WBVar::wb_pres, m, k, j, i, dum1, pl, dum2, pr, dum3);
          }
          src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol;
        }
        u0(m,IM1,k,j,i) += src;
        
        // Forces in the corotating frame
        Real wbt_a2 = 0.0, wbt_a3 = 0.0;   // problem/roche_wb_tan: Coriolis, per geometry
        if (use_spherical_polar) {
          Real vtheta = w0(m,IVY,k,j,i);
          Real vphi = w0(m,IVZ,k,j,i);
          Real vr = w0(m,IVX,k,j,i);
          Real sine = sin(theta);
          Real cosine = cos(theta);
          Real oor = SQR(omega)*r*sine;
          Real cor = 2.0*omega*vphi;
          if (wbt_) {
            // geometry front-end only: Coriolis on (e_theta, e_phi); ApplyTan below
            wbt_a2 = cor*cosine;
            wbt_a3 = -2.0*omega*(vr*sine+vtheta*cosine);
          } else {
          u0(m,IM2,k,j,i) += rho*(cor+oor)*cosine*bdt;
          u0(m,IM3,k,j,i) += -rho*2.0*omega*(vr*sine+vtheta*cosine)*bdt;
          }
          // rot_potential: the RADIAL centrifugal force lives in the potential (the
          // well-balanced source and, with etotgrav, the energy flux carry it); only the
          // theta part and Coriolis stay explicit here
          u0(m,IM1,k,j,i) += rho*(cor + (rotpot_src ? 0.0 : oor))*sine*bdt;
          // The centrifugal WORK: with rot_potential AND etotgrav the face potentials
          // phi0.x1f/x2f/x3f are TotPotAt (gravity + centrifugal, theta-dependent), so
          // the etotgrav energy flux rho v Phi_tot already does ALL of it, radial and
          // horizontal; adding the theta part here counted it twice (dhj-rotpot-energy,
          // 09-28: a spurious ~1e31 erg/s source/sink pair in the WASP-121b runs).
          // Without etotgrav (or without rot_potential) the source does all of it.
          if (!(rotpot && use_etotgrav)) {
            u0(m,IEN,k,j,i) += rho*oor*(((rotpot && use_etotgrav) ? 0.0 : vr*sine)
                                        + vtheta*cosine)*bdt;
          }
          // The host star's tidal term, on top of the planet's own centrifugal `oor`.
          // Together they make up the Hill acceleration Omega^2 (3x, 0, -z); see
          // TideAccR/TideAccT/TideAccP. mu is the substellar direction cosine, the same
          // sin(theta)*cos(phi) the radiative transfer uses for the stellar beam.
          if (tide) {
            const Real cosph = cos(phi);
            const Real sinph = sin(phi);
            const Real mu = sine*cosph;
            const Real atr = TideAccR(omega, r, mu, tide);
            const Real att = TideAccT(omega, r, sine, cosine, cosph, tide);
            const Real atp = TideAccP(omega, r, sine, sinph, cosph, tide);
            u0(m,IM1,k,j,i) += rho*atr*bdt;
            u0(m,IM2,k,j,i) += rho*att*bdt;
            u0(m,IM3,k,j,i) += rho*atp*bdt;
            // Tidal work. Unconditional, like the centrifugal work above: the tidal
            // potential is deliberately NOT carried in phicc0, so use_etotgrav does not
            // already account for it.
            u0(m,IEN,k,j,i) += rho*(atr*vr + att*vtheta + atp*vphi)*bdt;
          }
          // roche-dhj-1031: tangential -grad TidePot (stellar_tide is refused with it)
          if (rgms != 0.0 && !wbt_) {
            const Real cph = cos(phi), sph = sin(phi);
            Real gs, gy, gz;
            roche::TideGrad(rgms, rra, r*sine*cph, r*sine*sph, r*cosine, gs, gy, gz);
            u0(m,IM2,k,j,i) -= rho*(gs*cosine*cph + gy*cosine*sph - gz*sine)*bdt;
            u0(m,IM3,k,j,i) -= rho*(-gs*sph + gy*cph)*bdt;
          }
        } else if (use_cubed_sphere_) {
          // CUBED SPHERE: the same corotating-frame forces the spherical-polar branch
          // above applies -- full Coriolis AND centrifugal -- instead of the equatorial
          // beta-plane below, which is written for the CARTESIAN BOX and which the cubed
          // sphere used to fall through into.  The coordinate block a few lines up
          // branches three ways; this one branched only two, so the omission was silent.
          //
          // The beta-plane itself is not wrong -- f = 2 omega lam IS the standard
          // equatorial form (f = beta*y with beta = 2 omega/R and y = R lam) -- but it is
          // an EQUATORIAL approximation on a GLOBAL grid: against the correct 2 omega
          // sin(lam) it is 11 % too strong at 45 deg, 21 % at 60 deg and 57 % at the
          // pole, it drops the centrifugal term entirely, and its two terms assume an
          // orthonormal (x,y) pair while IM2/IM3 here are COVARIANT components on the
          // non-orthogonal gnomonic tangent basis.  Every cubed-sphere vs spherical-polar
          // comparison of this problem was confounded by that difference.
          //
          // Done in CARTESIAN and then projected, which avoids writing the curvilinear
          // Coriolis by hand in a non-orthogonal basis:
          //     a = -2 Omega x v + Omega^2 (x, y, 0),   Omega = omega zhat,
          // and since m_i = V.e_i is COVARIANT, the momentum source is rho (a.e_i) with
          // e_i = (rhat, e_xi, e_eta).  w0's two angular velocities are CONTRAVARIANT on
          // that same pair -- see Coordinates::GnomonicEquiangleRaiseVelMHD, which builds
          // them as v2 = (m2 - c m3)/(d det) -- so V is their plain combination with the
          // basis vectors, with no metric factor.  zhat is the rotation axis on both
          // grids: CSCellAngles takes the colatitude from the Cartesian z of the same
          // chart, in the convention the spherical-polar branch uses.
          const Real xi_c = 0.25*M_PI*x2v;
          const Real eta_c = 0.25*M_PI*x3v;
          const int pnl = mbpanel_.d_view(m);
          Real qc[3], e1[3], e2[3];
          cubed_sphere::PanelToCart(pnl, xi_c, eta_c, qc);
          cubed_sphere::PanelTangents(pnl, xi_c, eta_c, e1, e2);
          const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
          const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn, rh2 = qc[2]*qn;
          const Real v1 = w0(m,IVX,k,j,i);
          const Real v2 = w0(m,IVY,k,j,i);
          const Real v3 = w0(m,IVZ,k,j,i);
          const Real vcx = v1*rh0 + v2*e1[0] + v3*e2[0];
          const Real vcy = v1*rh1 + v2*e1[1] + v3*e2[1];
          const Real xx = r*rh0;
          const Real yy = r*rh1;
          const Real acx =  2.0*omega*vcy + SQR(omega)*xx;
          const Real acy = -2.0*omega*vcx + SQR(omega)*yy;
          // rot_potential: drop the radial projection of the centrifugal acceleration,
          // Omega^2 r sin^2(theta) = Omega^2 (x rh0 + y rh1), which the potential carries
          const Real acr_rot = rotpot_src ? SQR(omega)*(xx*rh0 + yy*rh1) : 0.0;
          u0(m,IM1,k,j,i) += rho*(acx*rh0 + acy*rh1 - acr_rot)*bdt;
          if (wbt_) {
            // geometry front-end only: Coriolis on (e_xi, e_eta); ApplyTan below
            const Real ccx =  2.0*omega*vcy, ccy = -2.0*omega*vcx;
            wbt_a2 = ccx*e1[0] + ccy*e1[1];
            wbt_a3 = ccx*e2[0] + ccy*e2[1];
          } else {
          u0(m,IM2,k,j,i) += rho*(acx*e1[0] + acy*e1[1])*bdt;
          u0(m,IM3,k,j,i) += rho*(acx*e2[0] + acy*e2[1])*bdt;
          }
          // Only the centrifugal part does work -- Coriolis is perpendicular to v.  With
          // rot_potential AND etotgrav the etotgrav flux rho v Phi_tot (TotPotAt on all
          // three face sets) already does the whole centrifugal work, so no energy
          // source may remain (the old code subtracted only the radial part and counted
          // the horizontal one twice; see the spherical-polar branch above).
          if (!(rotpot && use_etotgrav)) {
            u0(m,IEN,k,j,i) += rho*(SQR(omega)*(xx*vcx + yy*vcy)
                                    - ((rotpot && use_etotgrav) ?
                                       SQR(omega)*(xx*rh0 + yy*rh1)*v1 : 0.0))*bdt;
          }
          // roche-dhj-1031: tangential -grad TidePot.  The substellar point is along -x
          // (CSCellAngles: phi = atan2(-y, -x)), so the stellar frame is s = -x, y' = -y,
          // z; the Cartesian acceleration is (gs, gy', -gz), projected like the
          // centrifugal one above on e1, e2 (covariant momentum components)
          if (rgms != 0.0 && !wbt_) {
            Real gs, gy, gz;
            roche::TideGrad(rgms, rra, -r*rh0, -r*rh1, r*rh2, gs, gy, gz);
            u0(m,IM2,k,j,i) += rho*(gs*e1[0] + gy*e1[1] - gz*e1[2])*bdt;
            u0(m,IM3,k,j,i) += rho*(gs*e2[0] + gy*e2[1] - gz*e2[2])*bdt;
          }
        } else {
          // corotating beta-plane approximation e.g. Fromang+2016 (roche_wb_tan refused)
          Real omega1 = omega*lam;
          u0(m,IM2,k,j,i) += -2.0*rho*omega1*(-w0(m,IVZ,k,j,i))*bdt;
          u0(m,IM3,k,j,i) += -2.0*rho*omega1*w0(m,IVY,k,j,i)*bdt;
        }
        // problem/roche_wb_tan: ONE force path for sp and cs (wbtan_common.hpp): the
        // tangential centrifugal + tidal force rho g_t plus the front-end's Coriolis
        if (wbt_) wbtan::ApplyTan(u0, wbtg_, m, k, j, i, rho, wbt_a2, wbt_a3, bdt);

//        // Newtonian cooling
//        Real Teq, itrad;
//        get_eq_Tp(lam, phi, p, Teq);
//        get_itrad(p,itrad);
////        Real t0 = 2.16e6;
////        Real ff = (t0-time)/t0;
////        ff = (ff < 0.0) ? 0.0 : ff;
////        Real gg = 2.0*ff;
////        itrad /= pow(10.0,gg);
//        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
//        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;

        // Rayleigh drag (initial relaxation)
        Real dyntime = 2.0*M_PI/omega;
        Real tau1 = dyntime / 10.0;
        Real tau2 = dyntime;
        Real t1   = 2.0 * dyntime * rdf_;
        Real t2   = 5.0 * dyntime * rdf_;
        Real itdrag, fredux;
        if(time < t1) {
          itdrag = 1.0 / tau1;
        } else if(time < t2)
          {
            Real alp = (time - t1) / (t2 - t1);
            itdrag   = 1.0 / tau1 * pow(tau1 / tau2, alp);
          }
        else {
          itdrag = 0.0;
        }
        if (relax_ && time < t2) {
          fredux = itdrag*bdt;///(1.0+itdrag*bdt);
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }

        // Top sponge layer
        Real bar = 1.0e6;
        Real logpl = log(sp_phi_*bar);
        Real logpt = log(sp_plo_*bar);
        Real logp = log(p);
        Real fdrag = 1.0 - (logp-logpt)/(logpl-logpt); // high p = 0, low p = 1
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/sp_tau_;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
        if (sponge_top_) {
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          if (!sponge_top_rad_) {
            u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
            u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
          }
        }
        // Bottom sponge layer
        logpl = log(1.0e2*bar);
        logpt = log(5.0e1*bar);
        logp = log(p);
        fdrag = (logp-logpt)/(logpl-logpt); // high p = 1, low p = 0
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
//        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        if (sponge_bottom_) {
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }
        if (eb_on) {
          Kokkos::atomic_add(&ebacc(4, i-is),
                             (u0(m,IEN,k,j,i) - eb_e0)*volume(m,k,j,i)*eb_w);
          Kokkos::atomic_add(&ebacc(5, i-is),
                             (u0(m,IDN,k,j,i) - eb_d0)*volume(m,k,j,i)*eb_w);
        }
    });

    return;
}

void double_gray_two_stream_RT(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;
        
    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;
        
    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;

    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 100.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 5;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
    int nclip = 0;
    const Real demax = rt_de_max;
    par_reduce_clip3("2stream_rt", 0, nmb1, ks, ke, js, je, nclip,
    KOKKOS_LAMBDA(const int m, const int k, const int j, int &nc) {
        constexpr int NN = 270;
        Real tau_ir_down_f[NN];
        Real F_v_down_f[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real B[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real rtop = x1f_(m,ie+1);
        Real rbot = x1f_(m,is);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = -1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real pm1 = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie),w0(m,IEN,k,j,ie));
        Real pf = exp((log(p)+log(pm1))/2.0);
        Real tau_v_f = 0.0;//pf/(grav/kap_v);
        Real tau_ir_f = 0.0;//pf/(grav/kap_ir);
        
        tau_ir_down_f[ie+1] = tau_ir_f;
        F_v_down_f[ie+1] = Fstar*mu0;//*exp(-tau_v_f/muf);
        F_v_down_f[ie+1] = (mu0 > 0.0) ? F_v_down_f[ie+1] : 0.0;
        
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
          Real kap_ir = 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
            
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);

//          Real muf = sqrt(1.0-SQR(ap/rf)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r1));
          Real muf = (mu0 < mucr) ? mucr : mu0;
          if (test_oned) muf = mu0;
          Real dtau_v = kap_v*rho*dr;
          Real dtau_ir = kap_ir*rho*dr;
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf); // Zhang+2023
          Real trans = exp(-dtau_v/muf);
          F_v_down_f[i] = F_v_down_f[i+1]*trans*fac;
          tau_ir_down_f[i] = tau_ir_down_f[i+1] + dtau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        I_ir_down_f[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf);
          Real dtau = tau_ir_down_f[i]-tau_ir_down_f[i+1];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i]+B[i+1])/2.0;
          I_ir_down_f[i] = (I_ir_down_f[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        I_ir_up_f[is] = Iint;
        // up-sweep
        for (int i=is+1; i<ie+2; ++i) {
          Real rf = x1f_(m,i);
          Real rfm1 = x1f_(m,i-1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rfm1/rf);
          Real dtau = tau_ir_down_f[i-1]-tau_ir_down_f[i];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i-1]+B[i])/2.0;
          I_ir_up_f[i] = (I_ir_up_f[i-1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        // flux divergence
        for (int i=is; i<ie+1; ++i) {
          Real Ft = 2.0*M_PI*mu1*(I_ir_up_f[i+1]-I_ir_down_f[i+1])-F_v_down_f[i+1];
          Real Fb = 2.0*M_PI*mu1*(I_ir_up_f[i]-I_ir_down_f[i])-F_v_down_f[i];
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) src = -(Ft*area_t-Fb*area_b)/vol;
          Real de = src*bdt;
          if (demax > 0.0) {
            const Real dl = LimitRTSource(de, w0(m,IEN,k,j,i), demax);
            if (dl != de) { ++nc; de = dl; }
          }
          u0(m,IEN,k,j,i) += de;
        }
        
    });
    RTSourceLimiterWarn(nclip);
    
    return;
}

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;
    
    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real grav = 942.0;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 500.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 4;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(0), n1);
        constexpr int NN = 270;
        Real tau_ir_down[NN];
        Real I_down[NN];
        Real B[NN];
        Real Q_v[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(50.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real tau_v = 0.0;//p/(grav/kap_v);
        Real tau_ir = 0.0;//p/(grav/kap_ir);
//        if (test_oned) {
//          tau_v = p/(grav/kap_v)/mu0;
//          tau_ir = p/(grav/kap_ir)/mu1;
//        }
        tau_ir_down[ie+1] = tau_ir;
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
            Real kap_ir = 1.0e-2; // 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
            
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
////          Real mu = sqrt(1.0-SQR(r1/r)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r));
//          Real mu = (mu0 < mucr) ? mucr : mu0;
//          if (test_oned) mu = mu0;
//          Real dtau_v = kap_v*rho*dr;
          Real delta = dr/r;
          Real drcor = r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          if (test_oned) drcor = dr/mu0;
          Real dtau_v = kap_v*rho*drcor;
          Real dtau_ir = kap_ir*rho*dr;
            
          tau_v += dtau_v;
          Real fac = 1.0;
//          if (correct_spherical) fac = SQR(r1/r); // Zhang+2023
//          Real Q_v = kap_v*rho*Fstar*fac*exp(-tau_v/mu); // Zhang+2023
//          Q_v = (mu0 > 0.0) ? Q_v : 0.0;
          Real Qv = kap_v*rho*Fstar*fac*exp(-tau_v);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
//          u0(m,IEN,k,j,i) += Q_v*bdt;
        
          tau_ir_down[i] = tau_ir;
          tau_ir += dtau_ir/mu1;
          if (i==is) tau_ir_down[i-1] = tau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        Real rtop = x1v_(m,ie+1);
        I_down[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real r = x1v_(m,i);
          Real rp1 = x1v_(m,i+1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rp1/r);
          Real dtau = tau_ir_down[i]-tau_ir_down[i+1];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i)+B(i+1))/2.0;
          I_down[i] = (I_down[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        Real rbot = x1v_(m,is-1);
        Real I_up = Iint;
        // up-sweep
        for (int i=is; i<ie+1; ++i) {
          Real r = x1v_(m,i);
          Real rm1 = x1v_(m,i-1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rm1/r);
          Real dtau = tau_ir_down[i-1]-tau_ir_down[i];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i-1)+B(i))/2.0;
          I_up = (I_up*trans + Bavg*(1.0-trans))*fac;
          Real J = (I_up+I_down[i])/2.0;
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
            Real kap_ir = 1.0e-2; //2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
//          Real Q_ir = 4.0*M_PI*kap_ir*rho*(J-B(i));
//          u0(m,IEN,k,j,i) += Q_ir*bdt;
            
          Real kk = -4.0*M_PI*kap_ir*rho*boltz_sigma/M_PI*bdt;
          Real e0, e;
          if (eos.IsGeneral()) {
            // Same implicit balance e - kk T^4 = bb, but a general EOS has no e = c_v T
            // with constant c_v, so Newton-Raphson runs on the internal energy directly:
            // F(e) = e - kk T(e)^4 - bb, with dT/de = 1/(d c_v). c_v is per unit mass in
            // CODE units while T here is in Kelvin, hence the temp_cgs factor.
            e0 = w0(m,IEN,k,j,i);
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            e = e0;
            // `tc` is T in CODE temperature, carried alongside the kelvin `T` so that
            // neither EOS call has to solve for it. The two-argument SpecificHeatCv(d,e)
            // and Temperature(d,e) are the COLD-START forms -- the first is documented
            // "setup-time use only" because it inverts e(d,T) itself, and the second
            // brackets from scratch -- so using them here cost two full root finds per
            // Newton step, per cell, per stage. T and e are consistent at the top of
            // every iteration, so c_v can be evaluated at the temperature already in
            // hand, and the refresh below only needs the previous T as a warm start.
            Real tc = T/eos.temp_cgs;
            for (int n=0; n<100; ++n) {
              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e,tc));
              Real de = e - kk*SQR(SQR(T)) - bb;
              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
              if (fabs(de) <= 1.0e-10*e)
                break;
              tc = eos.Temperature(rho,e,tc);
              T = tc*eos.temp_cgs;
            }
          } else {
            Real cv = Rgas*rho*igm1;
            e0 = cv*T;
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            // Newton-Raphson
            for (int n=0; n<100; ++n) {
              e = cv*T;
              Real de = e - kk*SQR(SQR(T)) - bb;
              T -= de / (cv - 4.0*kk*T*T*T);
              if (fabs(de) <= 1.0e-10*e)
                break;
            }
          }
          u0(m,IEN,k,j,i) += e-e0;
        }
        
    });
    
    return;
}




//----------------------------------------------------------------------------------------
//! \fn void DhjCycleDiag
//! \brief PER-CYCLE, SINGLE-MESHBLOCK dump of everything that touches the energy of a
//! cell: the RT two-stream source (rate, applied de, clip flag), the radiative-conduction
//! face fluxes with the kappa / effective diffusivity / dt candidate behind them, the
//! optical-depth blend weights, and the hydro state.
//!
//! WHY.  A dt collapse in a cubed-sphere vertex column is a competition between the
//! explicit RT source and the flux-limited radiative conduction, and neither operator
//! appears in any output: the bin dumps carry only the state, at a cadence of a rotation.
//! One meshblock every cycle is ~10 MB and answers the question directly -- which of the
//! two put the energy in, and whether the source limiter was already saturated.
//!
//! Off unless problem/diag_gid >= 0, and then it costs one host copy of one meshblock.
//! Anything not allocated (the RT arrays before the first RT call, the conduction arrays
//! with no conduction module) is written as zeros, so the record layout is fixed.

//--------------------------------------------------------------------------------
//! \fn void DhjCkRtSplit
//! \brief problem/ck_impl_once.  The correlated-k radiation applied ONCE per
//! hydro step, with the FULL cycle dt, after the last RK stage
//! (ProblemGenerator::user_split_once).  Under rk2 that is one implicit column
//! solve per step instead of two, and the operator split is first order in dt --
//! the same order the in-stage source already had, since an SSP-RK stage average
//! of a per-stage relaxation is not second order in a stiff source either.

void DhjCkRtSplit(Mesh *pm, const Real dt) {
  if (dtloc_every > 0 && pm->ncycle % dtloc_every == 0 &&
      (pm->pmb_pack->phydro != nullptr || pm->pmb_pack->pmhd != nullptr) &&
      two_stream_rt::rt_pb_ptr != nullptr) {
    // the cell that set the fluid dt: hydro, or MHD (same members, mhd_newdt.cpp)
    hydro::Hydro *ph = pm->pmb_pack->phydro;
    mhd::MHD *pq = pm->pmb_pack->pmhd;
    const int dm = (ph != nullptr) ? ph->dtnew_m : pq->dtnew_m;
    const int dk = (ph != nullptr) ? ph->dtnew_k : pq->dtnew_k;
    const int dj = (ph != nullptr) ? ph->dtnew_j : pq->dtnew_j;
    const int di = (ph != nullptr) ? ph->dtnew_i : pq->dtnew_i;
    const Real dtf = (ph != nullptr) ? ph->dtnew : pq->dtnew;
    Real pbar = -1.0;
    if (dm >= 0) {
      auto sv = Kokkos::subview(*two_stream_rt::rt_pb_ptr, dm, dk, dj, di);
      auto hv = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sv);
      pbar = hv();
    }
    std::printf("### dtloc ncycle=%d rank=%d dt=%.5e %s=%.5e gid=%d k=%d j=%d i=%d "
                "ie=%d p_bar=%.4e\n", static_cast<int>(pm->ncycle),
                global_variable::my_rank, dt, (ph != nullptr) ? "dthyd" : "dtmhd", dtf,
                pm->pmb_pack->gids + std::max(dm, 0), dk, dj, di, pm->mb_indcs.ie, pbar);
  }
  // problem/ck_impl_every > 1: the cadence (full call every N cycles, linearised step in
  // between); see utils/two_stream_column_ck.hpp
  EBudCk(pm, dt, true);
  const int64_t nf0 = two_stream_rt::ck_cad_nfull;
  if (two_stream_rt::ck_impl_every > 1) {
    two_stream_rt::CkCadStep(pm, dt);
  } else {
    picket_fence_two_stream_RT(pm, dt);
  }
  EBudCk(pm, dt, false);
  if (dhj_ck_count_every > 0) {
    DhjCkCount((two_stream_rt::ck_impl_every <= 1) || (two_stream_rt::ck_cad_nfull != nf0));
    if (pm->ncycle % dhj_ck_count_every == 0) DhjCkT2Report(pm, false);
  }
}

//----------------------------------------------------------------------------------------
//! ckhesdirk-1102: problem/ck_time_scheme = hesdirk2.
//!
//! DESIGN.  The ck heating S(e) (the backward-Euler correlated-k column solve and the
//! cadence's linearised step CkCadLin, both BE in T) becomes the STIFF IMPLICIT part of
//! the IMEX H-ESDIRK2 step of rad_m1_time2.cpp, whose explicit part is the Heun hydro,
//! unchanged.  Implicit tableau A = [[0,0,0],[1-g,g,0],[1/2,b2,g]], g = 1 - 1/sqrt 2,
//! b2 = g/(2(1-g)), stiffly accurate, FSAL (K1 = the last slope of the previous step).
//! The ck runs IN each Heun stage (SourceFunc, where the in-stage ck call sat, after
//! RKUpdate and before the stage's own exchange; MHD: MHDSrcTerms before CT), no longer
//! in the operator split (user_split_func is not enrolled).  Each stage solve is the
//! EXISTING BE solve with dt -> g dt, its old vector = (stage start) + inc:
//!   stage 1: start = U^n + dt F(U^n),               inc = (1-g) dt K1
//!   stage 2: start = (U^n + Y2)/2 + dt/2 F(Y2),       inc = dt (g/2 K1 + (b2-g/2) K2)
//! (the Heun average supplies the rest of the tableau row, as in M1).  inc is added to
//! u0(IEN) BEFORE the solve, so the solve's old vector and first iterate are rhs; its
//! slope is K = (Y - rhs)/(g dt), i.e. what it actually applied (caps, fallbacks and
//! floors of the solve included): stage 1 -> K2, stage 2 -> the next K1.
//! Admissibility: a negative inc is clamped per cell to -(1 - ck_t2_pos_th) e_int
//! (default th 0.5; counted), no global redo.
//! Cadence (ck_impl_every > 1): stage 1 = CkCadStep as today (full on ncycle % N == 0);
//! stage 2 is ALWAYS the linearised step on the operator stage 1 (or the last full call)
//! stored (ck_cad_half = 2, shared with rkgfuse-1102 ck_strang), so the number of full solves per cycle is unchanged;
//! ck_impl_every = 1: both stages are full solves (2 per cycle).  xstep counts cycles.
//! Without a valid K1 (first step, a restart file without the T2 slabs) the step is
//! backward Euler: no ck in stage 1, one BE solve with the full dt in stage 2 (= the
//! Lie split of today, applied before the rest of stage 2's SourceFunc), K1 stored.
//! Restart: K1 and its validity travel in the ck restart block (slabs T2K1, T2OK), so a
//! continuation is the straight run.  dt_retry: K1 and ok are registered pgen state
//! (dtretry::RegisterView/Value) and come back with a retried cycle (dt_retry is not
//! wired for dhj without rad_m1 today).  MHD: the same hook; the RKG resistive step
//! stays operator split as today.  The key absent: nothing here runs (bitwise).

namespace {
const imex_sdirk2::Tableau kCkT2Tab = imex_sdirk2::Make(0);   // hesdirk2
const Real kCkT2G = kCkT2Tab.g;
DvceArray4D<Real> *t2_k2 = nullptr;
DvceArray4D<Real> *t2_r = nullptr;
Real t2_pos_th = 0.5;
double t2_nclamp = 0.0;
int64_t t2_nbe = 0, t2_nstep = 0;
bool t2_first = true;
int64_t cnt_full = 0, cnt_nfall = 0, cnt_ncap = 0, cnt_nonconv = 0;

void T2Alloc(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto u0 = (pmbp->phydro != nullptr) ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  auto &k1 = *two_stream_rt::ck_t2_k1_ptr;
  if (k1.size() > 0 && k1.extent(0) == u0.extent(0)) return;
  k1 = DvceArray4D<Real>("ck_t2_k1", u0.extent(0), u0.extent(2), u0.extent(3),
                         u0.extent(4));
  *t2_k2 = DvceArray4D<Real>("ck_t2_k2", u0.extent(0), u0.extent(2), u0.extent(3),
                             u0.extent(4));
  *t2_r = DvceArray4D<Real>("ck_t2_r", u0.extent(0), u0.extent(2), u0.extent(3),
                            u0.extent(4));
}
}  // namespace

void DhjCkCount(const bool full) {
  if (!full) return;
  namespace ts = two_stream_rt;
  ++cnt_full;
  cnt_nfall += ts::ck_impl_nfall;
  cnt_ncap += ts::ck_impl_ncap;
  cnt_nonconv += ts::ck_impl_nonconv;
}

void DhjCkT2Report(Mesh *pm, const bool final) {
  if (global_variable::my_rank != 0) return;
  namespace ts = two_stream_rt;
  std::printf("### ck_count%s ncycle=%d scheme=%s full=%lld nfall=%lld ncap=%lld "
              "nonconv=%lld cadfull=%lld cadlin=%lld t2step=%lld t2be=%lld "
              "t2clamp=%.0f\n", final ? "_final" : "", static_cast<int>(pm->ncycle),
              ts::ck_t2_on ? "hesdirk2" : "split", static_cast<long long>(cnt_full),
              static_cast<long long>(cnt_nfall), static_cast<long long>(cnt_ncap),
              static_cast<long long>(cnt_nonconv),
              static_cast<long long>(ts::ck_cad_nfull),
              static_cast<long long>(ts::ck_cad_nlin), static_cast<long long>(t2_nstep),
              static_cast<long long>(t2_nbe), t2_nclamp);
  std::fflush(stdout);
}

void DhjCkT2Init(ParameterInput *pin, Mesh *pm) {
  namespace ts = two_stream_rt;
  t2_pos_th = pin->GetOrAddReal("problem", "ck_t2_pos_th", 0.5);
  ts::ck_t2_k1_ptr = new DvceArray4D<Real>();
  t2_k2 = new DvceArray4D<Real>();
  t2_r = new DvceArray4D<Real>();
  dtretry::RegisterView("ck_t2_k1", ts::ck_t2_k1_ptr);
  dtretry::RegisterValue("ck_t2_ok", &ts::ck_t2_ok);
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: ck_time_scheme = hesdirk2 (IMEX H-ESDIRK2, ck "
              << "solved in each Heun stage with g dt, g = " << kCkT2G
              << "; cadence stage 2 linearised; pos_th " << t2_pos_th << ")" << std::endl;
  }
  (void) pm;
}

void DhjCkT2Stage(Mesh *pm, const Real bdt) {
  namespace ts = two_stream_rt;
  if (t2_first) {
    // the restart's K1 (or none: the first step is backward Euler)
    t2_first = false;
    T2Alloc(pm);
    ts::ck_t2_ok = false;
    if (ts::ck_rst_have && ts::CkRstStaged(ts::kCkSlabT2K1) &&
        ts::CkRstStaged(ts::kCkSlabT2Ok)) {
      const bool okv = (ts::ck_rst_stage[ts::kCkSlabT2Ok][0] > 0.5);
      if (ts::CkRstPut(ts::kCkSlabT2K1, *ts::ck_t2_k1_ptr)) ts::ck_t2_ok = okv;
      ts::CkRstStageFree(ts::kCkSlabT2Ok);
    } else if (global_variable::restart_run && global_variable::my_rank == 0) {
      std::cout << "### ck_time_scheme = hesdirk2: the restart file has no K1; the "
                << "first step is backward Euler" << std::endl;
    }
  }
  T2Alloc(pm);
  const Real dt = pm->dt;
  const int stage = (bdt == dt) ? 1 : 2;
  const Real g = kCkT2G;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto u0 = (pmbp->phydro != nullptr) ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto k1 = *ts::ck_t2_k1_ptr;
  auto k2 = *t2_k2;
  auto rr = *t2_r;
  const bool cad = (ts::ck_impl_every > 1);
  auto solve = [&](const Real sdt, const bool lin) {
    const int64_t nf0 = ts::ck_cad_nfull;
    if (cad) {
      // rkgfuse-1102's ck_cad_half: 1 = may be the cycle's full call, 2 = linearised
      ts::ck_cad_half = lin ? 2 : 1;
      ts::CkCadStep(pm, sdt);
      ts::ck_cad_half = 0;
    } else {
      picket_fence_two_stream_RT(pm, sdt);
    }
    if (dhj_ck_count_every > 0) DhjCkCount(!cad || (ts::ck_cad_nfull != nf0));
  };
  auto report = [&]() {
    if (dhj_ck_count_every > 0 && pm->ncycle % dhj_ck_count_every == 0) {
      DhjCkT2Report(pm, false);
    }
  };
  if (!ts::ck_t2_ok) {
    if (stage == 1) return;
    par_for("ck_t2_be0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      rr(m,k,j,i) = u0(m,IEN,k,j,i);
    });
    solve(dt, false);
    const Real idt = 1.0/dt;
    par_for("ck_t2_be1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      k1(m,k,j,i) = (u0(m,IEN,k,j,i) - rr(m,k,j,i))*idt;
    });
    ts::ck_t2_ok = true;
    ++t2_nbe;
    report();
    return;
  }
  const bool two = (stage == 2);
  Real c1, c2;
  imex_sdirk2::IncCoeffs(kCkT2Tab, stage, dt, c1, c2);
  const ts::CkCadState st = ts::CkCadMakeState(pm);
  const Real lof = 1.0 - t2_pos_th;
  Real ncl = 0.0;
  Kokkos::parallel_reduce("ck_t2_inc", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sum) {
    Real inc = c1*k1(m,k,j,i);
    if (two) inc += c2*k2(m,k,j,i);
    if (inc < 0.0) {
      const Real e = st.Eint(m, k, j, i);
      const Real lo = (e > 0.0) ? -lof*e : 0.0;
      if (inc < lo) {
        inc = lo;
        sum += 1.0;
      }
    }
    if (!(inc == inc)) inc = 0.0;
    u0(m,IEN,k,j,i) += inc;
    rr(m,k,j,i) = u0(m,IEN,k,j,i);
  }, Kokkos::Sum<Real>(ncl));
  t2_nclamp += ncl;
  solve(g*dt, two);
  const Real igdt = 1.0/(g*dt);
  auto kd = two ? k1 : k2;
  par_for("ck_t2_slope", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    kd(m,k,j,i) = (u0(m,IEN,k,j,i) - rr(m,k,j,i))*igdt;
  });
  if (two) {
    ++t2_nstep;
    report();
  }
}

//--------------------------------------------------------------------------------
//! \fn void DhjRtFreezeSplit
//! \brief problem/rt_test_freeze.  Reset the state to the frozen one (throwing away what
//! the RK stages did), apply the radiation once with dt = rt_test_dt (else the hydro dt),
//! write the diagnostic record, and keep the result as the new frozen state.  Test mode
//! only; see the note at rt_test_freeze.

void DhjRtFreezeSplit(Mesh *pm, const Real dt) {
  namespace ts = two_stream_rt;
  MeshBlockPack *pmbp = pm->pmb_pack;
  DvceArray5D<Real> u0 = pmbp->phydro->u0;
  auto &frz = *rt_frz_ptr;
  auto &frzT = *rt_frzT_ptr;
  DvceArray4D<Real> wt = pmbp->phydro->wtemp;
  Kokkos::deep_copy(u0, frz);
  Kokkos::deep_copy(wt, frzT);
  const Real dtr = (rt_test_dt > 0.0) ? rt_test_dt : dt;
  picket_fence_two_stream_RT(pm, dtr);
  ++rt_test_ncall;
  if (!rt_test_out.empty()) {
    auto &indcs = pm->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nmb = pmbp->nmb_thispack;
    auto eos = pmbp->phydro->peos->eos_data;
    const bool etg = pmbp->phydro->use_etotgrav;
    auto phi = pmbp->phydro->phicc0;
    // T [K], e, and c_v [erg/g/K] at the NEW state, on the device
    DvceArray5D<Real> dg("rt_test_dg", nmb, 3, u0.extent(2), u0.extent(3), u0.extent(4));
    par_for("rt_test_diag", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real d = u0(m,IDN,k,j,i);
      const Real e = EintFromCons(u0, m, k, j, i, 0.0, false, etg,
                                  etg ? phi(m,k,j,i) : 0.0);
      const Real t = eos.Temperature(d, e);
      dg(m,0,k,j,i) = t*eos.temp_cgs;
      dg(m,1,k,j,i) = e;
      dg(m,2,k,j,i) = eos.SpecificHeatCv(d, e, t)/eos.temp_cgs;
    });
    auto hdg = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dg);
    auto hu = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), u0);
    auto hf = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), frz);
    auto hx = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                  pmbp->pcoord->xx1f);
    const bool haveS = ts::ck_implicit && ts::ck_src_ptr != nullptr &&
                       ts::rt_icut_ptr != nullptr;
    DvceArray4D<Real>::HostMirror hs;
    DvceArray3D<int>::HostMirror hc;
    if (haveS) {
      hs = Kokkos::create_mirror_view(*ts::ck_src_ptr);
      Kokkos::deep_copy(hs, *ts::ck_src_ptr);
      hc = Kokkos::create_mirror_view(*ts::rt_icut_ptr);
      Kokkos::deep_copy(hc, *ts::rt_icut_ptr);
    }
    const int m0 = rt_test_col[0];
    const int k0 = (rt_test_col[1] >= 0) ? rt_test_col[1] : ks;
    const int j0 = (rt_test_col[2] >= 0) ? rt_test_col[2] : js;
    // E: every column against the reference column, BITWISE, on T and on u0(IEN)
    std::int64_t ndT = 0, ndE = 0;
    Real mxT = 0.0;
    // C: the column budget sum V (e - e*) - dt sum V S, and the solver's own residual
    // norm max |R|/(e + eps emax) re-evaluated on the FINAL state
    Real bud0 = 0.0, dE0 = 0.0, dS0 = 0.0, budmax = 0.0, rnmax = 0.0, ve0 = 0.0;
    for (int m=0; m<nmb; ++m) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          Real sE = 0.0, sS = 0.0, sA = 0.0, emax = 0.0, sV = 0.0;
          const int ic = haveS ? hc(m,k,j) : is;
          for (int i=is; i<=ie; ++i) {
            const Real t = hdg(m,0,k,j,i), tr = hdg(m0,0,k0,j0,i);
            if (t != tr) ++ndT;
            if (hu(m,IEN,k,j,i) != hu(m0,IEN,k0,j0,i)) ++ndE;
            const Real rl = std::fabs(t - tr)/tr;
            if (rl > mxT) mxT = rl;
            if (i >= ic && hdg(m,1,k,j,i) > emax) emax = hdg(m,1,k,j,i);
          }
          if (!haveS) continue;
          for (int i=ic; i<=ie; ++i) {
            const Real r0 = hx(m,i), r1 = hx(m,i+1);
            const Real v = (r1*r1*r1 - r0*r0*r0)/3.0;
            const Real de = hu(m,IEN,k,j,i) - hf(m,IEN,k,j,i);
            const Real sd = dtr*hs(m,k,j,i);
            sE += v*de;
            sS += v*sd;
            sA += v*std::fabs(sd);
            sV += v*hdg(m,1,k,j,i);
            const Real rn = std::fabs(de - sd)/
                            (hdg(m,1,k,j,i) + ts::ck_impl_norm_eps*emax);
            if (rn > rnmax) rnmax = rn;
          }
          const Real b = (sA > 0.0) ? std::fabs(sE - sS)/sA : 0.0;
          if (b > budmax) budmax = b;
          if (m == m0 && k == k0 && j == j0) {
            bud0 = (sA > 0.0) ? (sE - sS)/sA : 0.0;
            dE0 = sE;
            dS0 = sS;
            ve0 = sV;
          }
        }
      }
    }
    if (rt_test_fp == nullptr) {
      rt_test_fp = std::fopen(rt_test_out.c_str(), "w");
      std::fprintf(rt_test_fp, "# S ncall trt dt ndiffT ndiffE maxrelT passes nonconv "
                   "res gap capped fallback bud0 dE0 dtS0 budmax resfinal sumVe0\n");
      std::fprintf(rt_test_fp, "X");
      for (int i=is; i<=ie; ++i) std::fprintf(rt_test_fp, " %.17e",
                                              0.5*(hx(m0,i) + hx(m0,i+1)));
      std::fprintf(rt_test_fp, "\nD");
      for (int i=is; i<=ie; ++i) std::fprintf(rt_test_fp, " %.17e",
                                              hu(m0,IDN,k0,j0,i));
      std::fprintf(rt_test_fp, "\n");
      if (haveS) std::fprintf(rt_test_fp, "I %d\n", hc(m0,k0,j0) - is);
    }
    std::fprintf(rt_test_fp, "S %lld %.17e %.17e %lld %lld %.6e %d %d %.6e %.6e %d %d "
                 "%.9e %.9e %.9e %.6e %.6e %.9e\n", static_cast<long long>(rt_test_ncall),
                 rt_test_ncall*dtr, dtr, static_cast<long long>(ndT),
                 static_cast<long long>(ndE), mxT,
                 ts::ck_implicit ? ts::ck_impl_last_it : 0,
                 ts::ck_implicit ? ts::ck_impl_nonconv : 0, ts::ck_impl_last_res,
                 ts::ck_impl_last_gap, ts::ck_impl_ncap, ts::ck_impl_nfall,
                 bud0, dE0, dS0,
                 budmax, rnmax, ve0);
    if (rt_test_ncall == 1 || rt_test_ncall % rt_test_every == 0) {
      const char tag[3] = {'T', 'E', 'C'};
      for (int q=0; q<3; ++q) {
        std::fprintf(rt_test_fp, "%c %lld", tag[q],
                     static_cast<long long>(rt_test_ncall));
        for (int i=is; i<=ie; ++i) std::fprintf(rt_test_fp, " %.17e",
                                                hdg(m0,q,k0,j0,i));
        std::fprintf(rt_test_fp, "\n");
      }
    }
    std::fflush(rt_test_fp);
  }
  Kokkos::deep_copy(frz, u0);
  // the new state's temperature cache, solved with the old one as the guess
  {
    auto &indcs = pm->mb_indcs;
    auto eos = pmbp->phydro->peos->eos_data;
    const bool etg = pmbp->phydro->use_etotgrav;
    auto phi = pmbp->phydro->phicc0;
    par_for("rt_test_frzT", DevExeSpace(), 0, pmbp->nmb_thispack-1, indcs.ks, indcs.ke,
            indcs.js, indcs.je, indcs.is, indcs.ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real e = EintFromCons(u0, m, k, j, i, 0.0, false, etg,
                                  etg ? phi(m,k,j,i) : 0.0);
      frzT(m,k,j,i) = eos.Temperature(u0(m,IDN,k,j,i), e, frzT(m,k,j,i));
    });
    Kokkos::deep_copy(wt, frzT);
  }
}
