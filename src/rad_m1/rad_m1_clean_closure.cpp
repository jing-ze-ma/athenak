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

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_clean.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosure
//! \brief builds the closure for t^{n+1}.  Both stage solves and a backward-Euler redo end
//! at t^{n+1}, so one closure serves the whole step; it is rebuilt for a new step or a
//! new dt.

void RadiationM1::CleanClosure() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_havec && time == cln_tc && dt == cln_dtb) return;

  CleanPredictState();

  // formal solution of the predicted state, with its own source function (no relaxed
  // source or relaxed scattering J, which would make the closure depend on dt)
  vcol_srelax = false;
  vcol_sbl = false;
  vcol_sesc = false;
  vscat_jrel = false;
  VetShortChar();
  if (cln_opsaved) Kokkos::deep_copy(opac, cln_opb);

  CleanClosureFromRays();
  CleanClosureFill();
  CleanGhostFill(cln_s, 0, CLN_N);

  cln_tc = time;
  cln_dtb = dt;
  cln_havec = true;
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
//!     force on the gas.  A loss is applied at constant relative rate, E^n exp(r/E^n).
//!  3. exchange: the implicit local relaxation of (eint, E) over dt, starting from the
//!     hydro eint and the transported E.  It is stiff and local, so it is solved as such
//!     and cannot overshoot equilibrium.
//! Density is not changed by radiation within a step.  The prediction is within O(dt^2)
//! of U^{n+1}.  Without a valid K1 (first step) the coupling rates are zero.

void RadiationM1::CleanPredictState() {
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
  const bool rates = cln_pred && t2_ok && (t2k1.extent_int(0) > nmb1);

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
      const Real eint = fmax(uh(m,IEN,k,j,i) - ekin, 0.0);
      t = eos.Temperature(d, eint);
      if (use_rates) {
        for (int n = 0; n < 3; ++n) mom[n] += dt*k1(m,M1_T2_M1+n,k,j,i);
        const Real ekin_new = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
        const Real de = dt*(k1(m,M1_T2_E,k,j,i) + k1(m,M1_T2_EN,k,j,i))
                        - (ekin_new - ekin);
        e = (de >= 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
        // local exchange: E (1 + h c kE) = e + h c kP a T^4, eint(T) + E = eint + e
        const Real hce = dt*cl*opac_(m,M1_OP_E,k,j,i);
        const Real hcp = dt*cl*opac_(m,M1_OP_P,k,j,i)*arad_;
        t = CleanExchangeTemperature(eos, d, hcp/(1.0 + hce), eint + e*hce/(1.0 + hce), t);
        e = (e + hcp*SQR(SQR(t)))/(1.0 + hce);
      }
    } else if (use_rates) {
      const Real de = dt*k1(m,M1_T2_E,k,j,i);
      e = (de >= 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
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
//! \fn void RadiationM1::CleanClosureFromRays
//! \brief closure of every cell whose rays carry its radiation (J >= clean_ray_min E):
//! D = K/J, h+- = H+-/J, C1 = <|I - S|>/J.  Other cells are marked for the fill.

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

  par_for("cln_rays", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = vc(m,M1_VET_J,k,j,i);
    if (jr > 0.0 && jr >= ray_min*iw_(m,M1_IW_EN,k,j,i)) {
      for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = vc(m,M1_VET_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,M1_VET_HP1+a,k,j,i)/jr;
        const Real h = vc(m,M1_VET_H1+a,k,j,i)/jr;
        cs(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cs(m,CLN_HM+a,k,j,i) = fmin(h - hp, 0.0);
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
      cs(m,CLN_C1,k,j,i) = 1.0;
      cs(m,CLN_J,k,j,i) = 0.0;
      valid(m,0,k,j,i) = 0.0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanClosureFill
//! \brief cells without ray radiation take the mean closure (D, h+-) of their valid
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
      Real sum[12] = {0.0};
      int nv = 0;
      for (int kk = k-dk; kk <= k+dk; ++kk) {
        for (int jj = j-dj; jj <= j+dj; ++jj) {
          for (int ii = i-1; ii <= i+1; ++ii) {
            if (valid(m,0,kk,jj,ii) <= 0.0) continue;
            ++nv;
            for (int n = 0; n < 12; ++n) sum[n] += cs(m,CLN_D11+n,kk,jj,ii);
          }
        }
      }
      if (nv == 0) return;
      valid(m,1,k,j,i) = 1.0;
      for (int n = 0; n < 12; ++n) cs(m,CLN_D11+n,k,j,i) = sum[n]/nv;
    });
    par_for("cln_fill_mark", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,0,k,j,i) = valid(m,1,k,j,i);
    });
  }
}

}  // namespace radm1
