//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_newdt.cpp
//! \brief function to compute hydro timestep across all MeshBlock(s) in a MeshBlockPack

#include <math.h>

#include <limits>
#include <cstdlib>
#include <iostream>
#include <algorithm> // min
#include <array>
#include <vector>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro.hpp"
#include "hydro_mergelat.hpp"
#include "hydro_ringavg.hpp"
#include "diffusion/conduction.hpp"
#include "diffusion/viscosity.hpp"
#include "srcterms/srcterms.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "rad_vet/rad_vet_closure.hpp"

namespace hydro {
// <hydro>/dt_report (read only when named, in the Hydro constructor; 0 = off).
// Diagnostics only: no state change.
int dt_report_key = 0;
// <hydro>/dt_floor_courant (read only when named, in the Hydro constructor; 0 = off)
Real dt_floor_courant_key = 0.0;

//----------------------------------------------------------------------------------------
// \!fn void Hydro::NewTimeStep()
// \brief calculate the minimum timestep within a MeshBlockPack for hydrodynamic problems

TaskStatus Hydro::NewTimeStep(Driver *pdrive, int stage) {
  if (stage != (pdrive->nexp_stages)) {
    return TaskStatus::complete; // only execute last stage
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();
  // MinLoc alongside the three Mins: when dt collapses the only question that matters is
  // WHICH cell did it, and the location rides along for free.  See the report in
  // Mesh::NewTimeStep.
  Kokkos::ValLocScalar<Real, int> hloc;
  hloc.val = std::numeric_limits<float>::max();
  hloc.loc = -1;
  if (dt_diag.h_view.extent(0) == 0) Kokkos::realloc(dt_diag, ndtdiag);

  // <hydro>/dt_report = N > 0: the report cycles also keep every cell's CFL dt (dt_cell)
  // for the per-region minima (polar rows / rest)
  const int dtrep = dt_report_key;
  const bool dtcell_on = (dtrep > 0) && ((pmy_pack->pmesh->ncycle % dtrep) == 0) &&
                         (pmy_pack->pmesh->use_spherical_polar);
  if (dtcell_on && dt_cell.extent(0) != static_cast<size_t>(pmy_pack->nmb_thispack)) {
    Kokkos::realloc(dt_cell, pmy_pack->nmb_thispack, nx3, nx2, nx1);
  }
  auto dtc_ = dt_cell;
  // capture class variables for kernel
  auto &w0_ = w0;
  auto &eos = pmy_pack->phydro->peos->eos_data;
  // derived thermodynamic variables (general EOS only; empty view for an ideal gas)
  auto &wder_ = pmy_pack->phydro->wder;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &is_special_relativistic_ = pmy_pack->pcoord->is_special_relativistic;
  auto &is_general_relativistic_ = pmy_pack->pcoord->is_general_relativistic;
  auto &is_dynamical_relativistic_ = pmy_pack->pcoord->is_dynamical_relativistic;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  // CUBED SPHERE: the signal speed NORMAL to an angular face is not the stored velocity
  // component.  dx2 is the arc length along e_xi while the face separation is
  // dx2*sin_cell, and w0(IVY) is CONTRAVARIANT on the unit basis so the normal advection
  // speed is sin_cell*v^xi -- the two sin_cell factors CANCEL in dx2/|v2|, which is why
  // the kinematic branch below needs no correction.  The SOUND SPEED does not scale that
  // way: it is a physical speed, so the constraint is dx2/(|v2| + cs/sin_cell).  Leaving
  // it out overstates dt by up to 13.7% at a cube vertex (sin_cell = 0.88 there) and not
  // at all on the panel axes.  sp and Cartesian are untouched.
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  auto &sncell_ = pmy_pack->pcoord->sin_cell;
  const bool multi_d_ = pmy_pack->pmesh->multi_d;
  const bool three_d_ = pmy_pack->pmesh->three_d;
  // the merged cells' lateral CFL width is the GROUP width f dx2, f dx3
  const bool mlat_ = (pmlat != nullptr && pmlat->on);
  DvceArray2D<int> mlfi_;
  if (mlat_) {
    pmlat->Build();
    mlfi_ = pmlat->fi;
  }
  const bool mlw_ = (mlat_ && pmlat->radial);
  // the averaged polar rows' phi CFL width is the GROUP width N_j dx3
  const bool rav_ = (pravg != nullptr);
  DvceArray2D<Real> rnf_;
  if (rav_) rnf_ = pravg->nfac;

  // RADIATION SIGNAL SPEED (<hydro>/rad_signal_speed; default on with an implicit M1
  // transport, see hydro.cpp).  The gas sound speed alone understates the acoustic speed
  // where the M1 radiation is trapped: there P_rad takes part in the compression.  In
  // the strong-coupling, no-diffusion (equilibrium) limit gas + LTE radiation is one
  // fluid with P = P_g + P_r, beta = P_g/P, and the adiabatic exponent of an ideal gas of
  // ratio gamma mixed with radiation is (Chandrasekhar 1939, ch. II; Mihalas & Mihalas
  // 1984, sect. 101)
  //   Gamma_1 = beta + (4 - 3 beta)^2 (gamma - 1) / (beta + 12 (gamma - 1)(1 - beta)),
  // c_eq^2 = Gamma_1 P / rho  (-> gamma P_g/rho at beta = 1, -> 4/3 P_r/rho at beta = 0).
  // For a general EOS (Gamma_1 of the gas from the table) the decoupled sum
  // c_eq^2 = c_g^2 + (4/3) P_r/rho is used, an upper bound of the coupled value for an
  // ideal gas (by <= 3 % in speed at gamma = 5/3).  P_r along direction d is the
  // diagonal of the closure in use, P_dd = E [(1 - chi)/2 + (3 chi - 1)/2 n_d^2]:
  // (chi, n) of tau_ten (closure = tau | vet_col), of vet_cell (vet_sc; D_dd itself
  // with vet_tensor = full), M1Chi(|F|/cE) about F (m1 | minerbo | kershaw), 1/3
  // (eddington); E/3 wherever the closure gives no tensor yet (tau_ten before its first
  // build).  It is E/3 in the optically thick limit.  The radiation term is tapered by
  // the cell optical depth along each direction, tau_d = rho kappa_T dx_d (kappa_T the
  // module's transport opacity, opac(VET_OP_T)):
  //   c_d^2 = c_g^2 + (1 - exp(-tau_d)) (c_eq,d^2 - c_g^2),
  // so optically thin cells keep the gas speed.  The fully coupled speed is the LARGEST
  // acoustic speed of the medium (radiative diffusion at the grid scale only lowers it),
  // so the bound is conservative.  E, F, the tensor and opac are the M1 state of the
  // last radiation step (the module runs after the hydro stages): a lag of one step.
  // Before the first step they come from RadiationVet::RssPrime (Driver::Initialize):
  // the restart file's copy, or opac of the initial state on a fresh start.
  // Off: no change.
  // rad_signal_speed_inertia: the relativistic sound speed of gas + radiation, the
  // radiation enthalpy w_r = (E + P_dd)/c^2 in the inertia (Mihalas & Mihalas 1984, sect.
  // 101): c_eq^2 -> c_eq^2/(1 + w_r/rho), at most c^2/3 + c_g^2 where the radiation
  // dominates (a floored cell beside a bright field)
  // rad_signal_speed_diff: the taper is the grid-scale TRAPPING criterion instead of the
  // optical depth.  A radiation-pressure perturbation of wavenumber k = pi/h diffuses
  // away at the rate c_hat k^2/(3 rho kappa_T) and the wave oscillates at c_eq k: it
  // takes part in the compression only when diffusion is slower, tau_d c_eq/c_hat
  // > 3/pi ~ 1 (Mihalas & Mihalas 1984, sect. 101), so the weight is
  // 1 - exp(-tau_d c_eq,d/c_hat).  Cells with 1 < tau_d < c_hat/c_eq (a radiation-
  // dominated convective zone) then keep about the gas speed, which bounds the speed
  // of the radiatively damped wave from above; deep cells (tau_d c_eq/c_hat >> 1)
  // keep c_eq.
  // rad_signal_speed_phase: the weight is the IN-PHASE part of the radiation pressure
  // perturbation, 1 - exp(-x^2), x = tau_d c_eq,d/c_hat.  For x << 1 the perturbation
  // lags the compression by ~pi/2: its in-phase (restoring) part is O(x^2) and the O(x)
  // part damps the wave (the implicit exchange's job), so the linear root of gas +
  // moments stays at c_g until the mode is overdamped near x ~ 0.3 and reaches c_eq for
  // x >~ 1 (Mihalas & Mihalas 1984, sect. 101).  The weight is >= 1 - exp(-x) for
  // x >= 1, so dt is never larger where the radiation is trapped.  In thin gas c_eq^2
  // ~ 1/rho: the radiation term of csd^2 grows as rho^-1/2 with the linear weight and
  // tends to [(4/3) D_dd E kappa_T h/c_hat]^2 (density independent, << c_g^2) with
  // this one.
  const bool rss_ = rad_signal_speed, rsi_ = rss_inertia, rsd_ = rss_diff;
  const bool rsp_ = rss_phase;
  DvceArray5D<Real> erad_, kopc_, tten_, vcel_;
  radvet::CleanGeo gn_;
  int rss_mode = 0;          // 0 isotropic E/3, 1 M1Chi, 2 tau_ten, 3 vet_sc, 4 vet full
  int rss_kind = 0;
  Real rss_cl = 1.0, rss_ch = 1.0;
  if (rss_) {
    auto *prm = pmy_pack->pradvet;
    if (prm == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/rad_signal_speed = true needs a <rad_vet> block"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    erad_ = prm->u0;
    kopc_ = prm->opac;
    rss_cl = prm->c_light;
    rss_ch = prm->chat;
    rss_kind = prm->chi_kind;
    // the same selection decides what a restart file carries (RadiationVet::RssPrime)
    rss_mode = prm->RssMode();
    if (rss_mode == 2) tten_ = prm->tau_ten;
    if (rss_mode >= 3) vcel_ = (rss_mode == 5) ? prm->cln_s : prm->vet_cell;
    if (rss_mode == 5) gn_ = prm->cln_gn;   // (cubed sphere: D Cartesian)
  }

  // the radiation signal speeds csd of cell (m,k,j,i) at gas speed cs and widths h,
  // with the D_dd, E and kappa_T they use: one formula for the time step and its DT
  // report
  auto rssd = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            const Real cs, const Real h[3], Real csd[3], Real dd[3],
                            Real &er, Real &kt) {
    const Real d = w0_(m,IDN,k,j,i);
    for (int a = 0; a < 3; ++a) dd[a] = 1.0/3.0;
    er = fmax(erad_(m,radvet::VET_E,k,j,i), 0.0);
    // diagonal of the Eddington tensor, D_dd = P_dd/E
    if (rss_mode >= 4) {   // (5: clean's cln_s, 0 before its first closure: 1/3)
      const int d0 = (rss_mode == 5) ? radvet::CLN_D11 : radvet::VETRAY_D11;
      if (rss_mode == 4 || vcel_(m,d0,k,j,i) > 0.0) {
        for (int a = 0; a < 3; ++a) {
          dd[a] = (rss_mode == 5 && use_cubed_sphere) ?
                  radvet::CleanDuu(vcel_, gn_, (a == 0) ? radvet::CLG_R :
                                   radvet::CLG_E + 3*(a - 1), m, k, j, i) :
                  vcel_(m,d0+a,k,j,i);
        }
      }
    } else if (rss_mode >= 1) {
      Real chi, n1, n2, n3;
      if (rss_mode == 2) {
        chi = tten_(m,0,k,j,i); n1 = tten_(m,1,k,j,i);
        n2 = tten_(m,2,k,j,i); n3 = tten_(m,3,k,j,i);
      } else if (rss_mode == 3) {
        chi = vcel_(m,radvet::VETRAY_CHI,k,j,i);
        n1 = vcel_(m,radvet::VETRAY_N1,k,j,i);
        n2 = vcel_(m,radvet::VETRAY_N1+1,k,j,i);
        n3 = vcel_(m,radvet::VETRAY_N1+2,k,j,i);
      } else {
        const Real f1 = erad_(m,radvet::VET_F1,k,j,i);
        const Real f2 = erad_(m,radvet::VET_F2,k,j,i);
        const Real f3 = erad_(m,radvet::VET_F3,k,j,i);
        const Real fn = sqrt(f1*f1 + f2*f2 + f3*f3);
        const Real inv = (fn > 0.0) ? 1.0/fn : 0.0;
        chi = (er > 0.0) ? radvet::M1Chi(fn/(rss_cl*er), rss_kind) : 1.0/3.0;
        n1 = f1*inv; n2 = f2*inv; n3 = f3*inv;
      }
      if (chi > 0.0) {   // tau_ten is zero before its first build: keep 1/3
        const Real dg = 0.5*(1.0 - chi), an = 0.5*(3.0*chi - 1.0);
        dd[0] = dg + an*n1*n1;
        dd[1] = dg + an*n2*n2;
        dd[2] = dg + an*n3*n3;
      }
    }
    const Real cg2 = cs*cs;
    kt = fmax(kopc_(m,radvet::VET_OP_T,k,j,i), 0.0);
    Real pg = 0.0, gm1 = 0.0;
    const bool ideal = (eos.is_ideal && !eos.IsGeneral());
    if (ideal) {
      pg = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
      gm1 = eos.gamma - 1.0;
    }
    for (int a = 0; a < 3; ++a) {
      const Real pr = fmax(dd[a], 0.0)*er;
      Real ceq2;
      if (ideal) {
        const Real ptot = pg + pr;
        const Real beta = pg/ptot;
        const Real gam1 = beta + SQR(4.0 - 3.0*beta)*gm1
                                 /(beta + 12.0*gm1*(1.0 - beta));
        ceq2 = gam1*ptot/d;
      } else {
        ceq2 = cg2 + (4.0/3.0)*pr/d;
      }
      if (rsi_) ceq2 /= 1.0 + (er + pr)/(d*rss_cl*rss_cl);
      const Real dc2 = fmax(ceq2 - cg2, 0.0);
      const Real tw = rsd_ ? kt*h[a]*sqrt(ceq2)/rss_ch : kt*h[a];
      csd[a] = sqrt(cg2 + dc2*(1.0 - exp(-(rsp_ ? tw*tw : tw))));
    }
  };

  // <hydro>/dt_floor_courant = C (read only when named; 0 = off = bitwise): a cell of
  // floor gas, rho < <rad_vet>/floorcap_rho, takes dt <= C/sum_d (|v_d| + c_d)/h_d
  // instead of cfl_number h_d/(|v_d| + c_d), so the floor atmosphere does not set dt
  const Real dfc_ = (dt_floor_courant_key > 0.0) ?
                    dt_floor_courant_key/pmy_pack->pmesh->cfl_no : 0.0;
  Real fcr_ = 0.0;
  if (dfc_ > 0.0) {
    if (pmy_pack->pradvet == nullptr || pmy_pack->pradvet->fcap_rho <= 0.0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/dt_floor_courant needs <rad_vet>/floorcap_rho > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    fcr_ = pmy_pack->pradvet->fcap_rho;
  }

  if (pdrive->time_evolution == TimeEvolution::kinematic) {
    // find smallest (dx/v) in each direction for advection problems
    Kokkos::parallel_reduce("HydroNudt1",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
      // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      if (use_cubed_sphere || use_spherical_polar) {
        const Real mf = mlat_ ? static_cast<Real>(mlfi_(m,i)) : 1.0;
        min_dt1 = fmin((dx1_(m,k,j,i)/fabs(w0_(m,IVX,k,j,i))), min_dt1);
        if (mlat_) {
          min_dt2 = fmin((mf*dx2_(m,k,j,i)/fabs(w0_(m,IVY,k,j,i))), min_dt2);
          min_dt3 = fmin((mf*dx3_(m,k,j,i)/fabs(w0_(m,IVZ,k,j,i))), min_dt3);
        } else {
        const Real h3 = rav_ ? rnf_(m,j)*dx3_(m,k,j,i) : dx3_(m,k,j,i);
        min_dt2 = fmin((dx2_(m,k,j,i)/fabs(w0_(m,IVY,k,j,i))), min_dt2);
        min_dt3 = fmin((h3/fabs(w0_(m,IVZ,k,j,i))), min_dt3);
        }
      } else {
      min_dt1 = fmin((mbsize.d_view(m).dx1/fabs(w0_(m,IVX,k,j,i))), min_dt1);
      min_dt2 = fmin((mbsize.d_view(m).dx2/fabs(w0_(m,IVY,k,j,i))), min_dt2);
      min_dt3 = fmin((mbsize.d_view(m).dx3/fabs(w0_(m,IVZ,k,j,i))), min_dt3);
      }
    }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2),Kokkos::Min<Real>(dt3));
  } else {
    // find smallest dx/(v +/- Cs) in each direction for hydrodynamic problems
    Kokkos::parallel_reduce("HydroNudt2",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3,
                  Kokkos::ValLocScalar<Real, int> &mres) {
      // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      Real max_dv1 = 0.0, max_dv2 = 0.0, max_dv3 = 0.0;

      if (is_general_relativistic_ || is_dynamical_relativistic_) {
        max_dv1 = 1.0;
        max_dv2 = 1.0;
        max_dv3 = 1.0;
      } else if (is_special_relativistic_) {
        Real v2 = SQR(w0_(m,IVX,k,j,i)) + SQR(w0_(m,IVY,k,j,i)) + SQR(w0_(m,IVZ,k,j,i));
        Real lor = sqrt(1.0 + v2);
        // FIXME ERM: Ideal fluid for now
        Real p = eos.IdealGasPressure(w0_(m,IEN,k,j,i));

        Real lm, lp;
        eos.IdealSRHydroSoundSpeeds(w0_(m,IDN,k,j,i), p, w0_(m,IVX,k,j,i), lor, lp, lm);
        max_dv1 = fmax(fabs(lm), lp);

        eos.IdealSRHydroSoundSpeeds(w0_(m,IDN,k,j,i), p, w0_(m,IVY,k,j,i), lor, lp, lm);
        max_dv2 = fmax(fabs(lm), lp);

        eos.IdealSRHydroSoundSpeeds(w0_(m,IDN,k,j,i), p, w0_(m,IVZ,k,j,i), lor, lp, lm);
        max_dv3 = fmax(fabs(lm), lp);
      } else {
        Real cs;
        if (eos.IsGeneral()) {
          // pressure and Gamma_1 were evaluated once per cell in ConsToPrim. Using them
          // here matters: with ionization or radiation pressure the ideal-gas expression
          // below gives the wrong sound speed, and hence the wrong CFL timestep.
          cs = eos.SoundSpeedFromP(w0_(m,IDN,k,j,i), wder_(m,IDPR,k,j,i),
                                   wder_(m,IDG1,k,j,i));
        } else if (eos.is_ideal) {
          Real p = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
          cs = eos.IdealHydroSoundSpeed(w0_(m,IDN,k,j,i), p);
        } else         {
          cs = eos.iso_cs;
        }
        Real cs1 = cs, cs2 = cs, cs3 = cs;
        if (rss_) {
          Real h[3];
          if (use_cubed_sphere || use_spherical_polar) {
            h[0] = dx1_(m,k,j,i); h[1] = dx2_(m,k,j,i); h[2] = dx3_(m,k,j,i);
            if (mlat_) {
              const Real mf = mlw_ ? MlatWidth(w0_, mlfi_(m,i), js, ks, m, k, j, i) :
                              static_cast<Real>(mlfi_(m,i));
              h[1] *= mf;
              h[2] *= mf;
            }
          } else {
            h[0] = mbsize.d_view(m).dx1; h[1] = mbsize.d_view(m).dx2;
            h[2] = mbsize.d_view(m).dx3;
          }
          Real csd[3], dd[3], er, kt;
          rssd(m, k, j, i, cs, h, csd, dd, er, kt);
          cs1 = csd[0]; cs2 = csd[1]; cs3 = csd[2];
        }
        max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cs1;
        max_dv2 = fabs(w0_(m,IVY,k,j,i))
                 + (cs_ ? cs2/sncell_(m,k,j) : cs2);
        max_dv3 = fabs(w0_(m,IVZ,k,j,i))
                 + (cs_ ? cs3/sncell_(m,k,j) : cs3);
      }
      // NON-FINITE GUARD: a NaN/inf state gives a NaN dx/v, which fmin DROPS, so dt
      // would be set by the finite cells alone and could keep doubling (the whole domain
      // NaN, dt doubled per step to tlim, rc 0).  Flag the cell
      // with a negative dt instead; Mesh::NewTimeStep stops the run on dt <= 0.
      if (!(max_dv1 + max_dv2 + max_dv3 < 1.0e300)) {
        min_dt1 = -1.0;
        mres.val = -1.0;
        mres.loc = idx;
        return;
      }
      // the per-direction CFL widths h_d/(|v_d| + c_d)
      Real q1, q2, q3;
      const bool curv = (use_cubed_sphere || use_spherical_polar);
      if (curv) {
        q1 = dx1_(m,k,j,i)/max_dv1;
        if (mlat_) {
          const Real mf = mlw_ ? MlatWidth(w0_, mlfi_(m,i), js, ks, m, k, j, i) :
                          static_cast<Real>(mlfi_(m,i));
          q2 = mf*dx2_(m,k,j,i)/max_dv2;
          q3 = mf*dx3_(m,k,j,i)/max_dv3;
        } else {
          q2 = dx2_(m,k,j,i)/max_dv2;
          q3 = (rav_ ? rnf_(m,j)*dx3_(m,k,j,i) : dx3_(m,k,j,i))/max_dv3;
        }
      } else {
        q1 = mbsize.d_view(m).dx1/max_dv1;
        q2 = mbsize.d_view(m).dx2/max_dv2;
        q3 = mbsize.d_view(m).dx3/max_dv3;
      }
      // <hydro>/dt_floor_courant: floor gas (rho < floorcap_rho) bound on the SUM of the
      // directions' Courant numbers instead of cfl_number per direction
      if (dfc_ > 0.0 && w0_(m,IDN,k,j,i) < fcr_) {
        Real s = 1.0/q1;
        if (multi_d_) s += 1.0/q2;
        if (three_d_) s += 1.0/q3;
        q1 = q2 = q3 = dfc_/s;
      }
      min_dt1 = fmin(q1, min_dt1);
      min_dt2 = fmin(q2, min_dt2);
      min_dt3 = fmin(q3, min_dt3);
      Real cell_dt = q1;
      if (multi_d_) cell_dt = fmin(cell_dt, q2);
      if (three_d_) cell_dt = fmin(cell_dt, q3);
      if (dtcell_on && curv && !mlat_) {dtc_(m,k-ks,j-js,i-is) = cell_dt;}
      if (cell_dt < mres.val) { mres.val = cell_dt; mres.loc = idx; }
    }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2),Kokkos::Min<Real>(dt3),
       Kokkos::MinLoc<Real, int>(hloc));
  }

  // compute minimum of dt1/dt2/dt3 for 1D/2D/3D problems
  dtnew = dt1;
  if (pmy_pack->pmesh->multi_d) { dtnew = std::min(dtnew, dt2); }
  if (pmy_pack->pmesh->three_d) { dtnew = std::min(dtnew, dt3); }

  // decode the winning cell, and -- only when dt has just collapsed -- its state
  if (hloc.loc >= 0 && hloc.loc < nmkji) {
    dtnew_m = (hloc.loc)/nkji;
    dtnew_k = (hloc.loc - dtnew_m*nkji)/nji + ks;
    dtnew_j = (hloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji)/nx1 + js;
    dtnew_i = (hloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji - (dtnew_j-js)*nx1) + is;
  } else {
    dtnew_m = dtnew_k = dtnew_j = dtnew_i = -1;
  }
  dt_diag_valid = false;
  if (dtnew_m >= 0 && (dtnew_prev < 0.0 || dtnew < 0.25*dtnew_prev)) {
    auto dd = dt_diag;
    const int dm = dtnew_m, dk = dtnew_k, dj = dtnew_j, di = dtnew_i;
    auto &x1v_ = pmy_pack->pcoord->x1v;
    auto &wtemp_ = pmy_pack->phydro->wtemp;
    auto eos_ = eos;
    const bool curvilinear_ = (use_cubed_sphere || use_spherical_polar);
    const int is_ = is, nx1_ = nx1;
    par_for("hyd_dtdiag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      const Real d = w0_(dm,IDN,dk,dj,di);
      Real pr, cs;
      if (eos_.IsGeneral()) {
        pr = wder_(dm,IDPR,dk,dj,di);
        cs = eos_.SoundSpeedFromP(d, pr, wder_(dm,IDG1,dk,dj,di));
      } else if (eos_.is_ideal) {
        pr = eos_.IdealGasPressure(w0_(dm,IEN,dk,dj,di));
        cs = eos_.IdealHydroSoundSpeed(d, pr);
      } else {
        pr = d*SQR(eos_.iso_cs);
        cs = eos_.iso_cs;
      }
      // x1v is allocated only on spherical-polar/cubed-sphere grids; on a Cartesian
      // grid it is a (1,1) dummy, so build the cell-centre position from the block size
      dd.d_view(0) = curvilinear_ ? x1v_(dm,di)
                   : CellCenterX(di-is_, nx1_, mbsize.d_view(dm).x1min,
                                 mbsize.d_view(dm).x1max);
      dd.d_view(1) = d;
      dd.d_view(2) = eos_.IsGeneral() ? wtemp_(dm,dk,dj,di) : (pr/d);
      dd.d_view(3) = pr;
      dd.d_view(4) = cs;
      dd.d_view(5) = w0_(dm,IVX,dk,dj,di);
      dd.d_view(6) = w0_(dm,IVY,dk,dj,di);
      dd.d_view(7) = w0_(dm,IVZ,dk,dj,di);
    });
    dt_diag.modify_device();
    dt_diag.sync_host();
    dt_diag_valid = true;
  }
  dtnew_prev = dtnew;
  // <hydro>/merge_lat_diag: the global gas M, E, P, L sums
  if (pmlat != nullptr && pmlat->diag_every > 0 &&
      (pmy_pack->pmesh->ncycle % pmlat->diag_every) == 0) {
    pmlat->Diag(u0, "h");
  }

  // dt report (<hydro>/dt_report = N > 0, every N cycles, no state change): the rank
  // holding the global minimum hydro dt prints the limiting cell, its direction, gas and
  // radiation-modified signal speeds, the taper optical depth rho kappa_T dx, E_rad,
  // c dt rho kappa_P and the thin-freeze flag; the global minimum of the GAS-ONLY dt
  // dx/(|v| + c_gas) is printed alongside.
  if (dtrep > 0 && pdrive->time_evolution != TimeEvolution::kinematic &&
      !is_special_relativistic_ && !is_general_relativistic_ &&
      !is_dynamical_relativistic_ && (pmy_pack->pmesh->ncycle % dtrep) == 0) {
    // gas-only dt (no radiation term)
    Real dtg = std::numeric_limits<float>::max();
    Kokkos::parallel_reduce("HydroNudtGas",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &mn) {
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real cs;
      if (eos.IsGeneral()) {
        cs = eos.SoundSpeedFromP(w0_(m,IDN,k,j,i), wder_(m,IDPR,k,j,i),
                                 wder_(m,IDG1,k,j,i));
      } else if (eos.is_ideal) {
        cs = eos.IdealHydroSoundSpeed(w0_(m,IDN,k,j,i),
                                      eos.IdealGasPressure(w0_(m,IEN,k,j,i)));
      } else {
        cs = eos.iso_cs;
      }
      Real h1, h2, h3;
      if (use_cubed_sphere || use_spherical_polar) {
        h1 = dx1_(m,k,j,i); h2 = dx2_(m,k,j,i); h3 = dx3_(m,k,j,i);
      } else {
        h1 = mbsize.d_view(m).dx1; h2 = mbsize.d_view(m).dx2; h3 = mbsize.d_view(m).dx3;
      }
      const Real sc = cs_ ? 1.0/sncell_(m,k,j) : 1.0;
      Real c = h1/(fabs(w0_(m,IVX,k,j,i)) + cs);
      if (multi_d_) c = fmin(c, h2/(fabs(w0_(m,IVY,k,j,i)) + cs*sc));
      if (three_d_) c = fmin(c, h3/(fabs(w0_(m,IVZ,k,j,i)) + cs*sc));
      mn = fmin(mn, c);
    }, Kokkos::Min<Real>(dtg));
    struct {double v; int r;} loc{static_cast<double>(dtnew), global_variable::my_rank};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &loc, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dtg, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
    // the cell diagnostics for the global minimum and (sp) the minima of
    // the polar rows (7 rows from each pole; lines DTREG polar) and of the rest (DTREG
    // rest), each printed by its owning rank
    auto pmy_x2v_ = pmy_pack->pcoord->x2v;
    auto pmy_x3v_ = pmy_pack->pcoord->x3v;
    // the targets: kind 0 the global minimum (its owner), 1 polar rows, 2 the rest
    std::vector<std::array<int, 7>> tgt;
    if (loc.r == global_variable::my_rank && dtnew_m >= 0) {
      tgt.push_back({0, dtnew_m, dtnew_k, dtnew_j, dtnew_i, 0, 0});
    }
    if (dtcell_on) {
      // the polar rows (global row < npol from either pole) and the rest
      const int npol = 7;
      const int nrow = pmy_pack->pmesh->mesh_indcs.nx2;
      auto &mbl = pmy_pack->pmesh->lloc_eachmb;
      std::vector<int> lx2h(pmy_pack->nmb_thispack);
      for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
        lx2h[m] = mbl[pmy_pack->gids + m].lx2;
      }
      DvceArray1D<int> lx2d("dtrep_lx2", pmy_pack->nmb_thispack);
      {
        auto h = Kokkos::create_mirror_view(lx2d);
        for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {h(m) = lx2h[m];}
        Kokkos::deep_copy(lx2d, h);
      }
      for (int reg = 0; reg < 2; ++reg) {
        Kokkos::ValLocScalar<Real, int> rl;
        rl.val = std::numeric_limits<float>::max();
        rl.loc = -1;
        const int n2b = nx2, npl = npol;
        Kokkos::parallel_reduce("HydroDtReg",
                                Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
        KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &mr) {
          const int m = idx/nkji;
          const int k = (idx - m*nkji)/nji;
          const int j = (idx - m*nkji - k*nji)/nx1;
          const int i = idx - m*nkji - k*nji - j*nx1;
          const int jg = lx2d(m)*n2b + j;
          const int jp = (jg < nrow - 1 - jg) ? jg : (nrow - 1 - jg);
          const bool pol = (jp < npl);
          if (pol != (reg == 0)) {return;}
          const Real v = dtc_(m,k,j,i);
          if (v < mr.val) {mr.val = v; mr.loc = idx;}
        }, Kokkos::MinLoc<Real, int>(rl));
        struct {double v; int r;} rg{static_cast<double>(rl.val),
                                     global_variable::my_rank};
#if MPI_PARALLEL_ENABLED
        MPI_Allreduce(MPI_IN_PLACE, &rg, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
#endif
        if (rg.r == global_variable::my_rank && rl.loc >= 0) {
          const int m = rl.loc/nkji;
          const int k = (rl.loc - m*nkji)/nji + ks;
          const int j = (rl.loc - m*nkji - (k-ks)*nji)/nx1 + js;
          const int i = (rl.loc - m*nkji - (k-ks)*nji - (j-js)*nx1) + is;
          tgt.push_back({1 + reg, m, k, j, i, (lx2h[m]*nx2 + j - js), npol});
        }
      }
    }
    for (const auto &tg : tgt) {
      DualArray1D<Real> rp("dtrep", 32);
      const int dm = tg[1], dk = tg[2], dj = tg[3], di = tg[4];
      auto &x1v_ = pmy_pack->pcoord->x1v;
      auto &wtemp_ = pmy_pack->phydro->wtemp;
      const bool curv = (use_cubed_sphere || use_spherical_polar);
      DvceArray4D<Real> tfz_;
      bool havetf = false;
      if (pmy_pack->pradvet != nullptr) {
        tfz_ = pmy_pack->pradvet->thin_frz;
        havetf = (tfz_.extent(0) > 0);
      }
      const Real cdt = (pmy_pack->pradvet != nullptr ? pmy_pack->pradvet->c_light : 0.0)
                       *pmy_pack->pmesh->dt;
      const int is_ = is, nx1_ = nx1;
      auto rv = rp.d_view;
      par_for("hyd_dtrep", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
        const int m = dm, k = dk, j = dj, i = di;
        const Real d = w0_(m,IDN,k,j,i);
        Real pr, cs, g1 = 0.0;
        if (eos.IsGeneral()) {
          pr = wder_(m,IDPR,k,j,i);
          g1 = wder_(m,IDG1,k,j,i);
          cs = eos.SoundSpeedFromP(d, pr, g1);
        } else if (eos.is_ideal) {
          pr = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
          cs = eos.IdealHydroSoundSpeed(d, pr);
        } else {
          pr = d*SQR(eos.iso_cs);
          cs = eos.iso_cs;
        }
        Real h[3];
        if (curv) {
          h[0] = dx1_(m,k,j,i); h[1] = dx2_(m,k,j,i); h[2] = dx3_(m,k,j,i);
          if (rav_) {h[2] *= rnf_(m,j);}   // the ring-avg group width
        } else {
          h[0] = mbsize.d_view(m).dx1; h[1] = mbsize.d_view(m).dx2;
          h[2] = mbsize.d_view(m).dx3;
        }
        Real csd[3] = {cs, cs, cs};
        Real er = 0.0, kt = 0.0, kp = 0.0;
        Real dd[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
        if (rss_) {
          rssd(m, k, j, i, cs, h, csd, dd, er, kt);
          kp = kopc_(m,radvet::VET_OP_P,k,j,i);
        }
        const Real sc = cs_ ? 1.0/sncell_(m,k,j) : 1.0;
        const Real v[3] = {w0_(m,IVX,k,j,i), w0_(m,IVY,k,j,i), w0_(m,IVZ,k,j,i)};
        int adir = 0;
        Real best = 1.0e300;
        for (int a = 0; a < 3; ++a) {
          if ((a == 1 && !multi_d_) || (a == 2 && !three_d_)) continue;
          const Real t = h[a]/(fabs(v[a]) + csd[a]*(a > 0 ? sc : 1.0));
          if (t < best) {best = t; adir = a;}
        }
        Real fl = 0.0;   // dt_floor_courant: the floor cell's summed bound
        if (dfc_ > 0.0 && d < fcr_) {
          Real sm = 0.0;
          for (int a = 0; a < 3; ++a) {
            if ((a == 1 && !multi_d_) || (a == 2 && !three_d_)) continue;
            sm += (fabs(v[a]) + csd[a]*(a > 0 ? sc : 1.0))/h[a];
          }
          best = dfc_/sm;
          fl = 1.0;
        }
        rv(0) = curv ? x1v_(m,i) : CellCenterX(i-is_, nx1_, mbsize.d_view(m).x1min,
                                               mbsize.d_view(m).x1max);
        rv(1) = d;
        rv(2) = eos.IsGeneral() ? wtemp_(m,k,j,i) : (pr/d);
        rv(3) = pr;
        rv(4) = g1;
        rv(5) = cs;
        rv(6) = adir + 1;
        rv(7) = v[adir];
        rv(8) = csd[adir];
        rv(9) = h[adir];
        rv(10) = kt*h[adir];
        rv(11) = er;
        rv(12) = dd[adir];
        rv(13) = cdt*kp;
        rv(14) = havetf ? tfz_(m,k,j,i) : -1.0;
        rv(15) = best;
        rv(16) = h[adir]/(fabs(v[adir]) + cs*(adir > 0 ? sc : 1.0));
        rv(17) = sc;
        rv(18) = (curv ? pmy_x2v_(m,j) : 0.0);
        rv(19) = (curv ? pmy_x3v_(m,k) : 0.0);
        rv(20) = fl;
      });
      rp.modify_device();
      rp.sync_host();
      auto &q = rp.h_view;
      if (tg[0] == 0) {
      const Real cf = pmy_pack->pmesh->cfl_no;
      std::cout << "DTREP cycle=" << pmy_pack->pmesh->ncycle
                << " time=" << pmy_pack->pmesh->time
                << " cfl*dt_hyd=" << cf*dtnew << " cfl*dt_gasonly=" << cf*dtg
                << " rank=" << global_variable::my_rank
                << " gid=" << (pmy_pack->gids + dm) << " kji=" << dk << "," << dj << ","
                << di << " dir=" << static_cast<int>(q(6)) << " r=" << q(0)
                << " rho=" << q(1) << " T=" << q(2) << " pg=" << q(3) << " G1=" << q(4)
                << " cgas=" << q(5) << " v=" << q(7) << " ceff=" << q(8)
                << " h=" << q(9) << " tauT_h=" << q(10) << " Erad=" << q(11)
                << " Ddd=" << q(12) << " cdt_kP=" << q(13) << " frz=" << q(14)
                << " sc=" << q(17) << " cell_dt=" << q(15) << " cell_dt_gas=" << q(16)
                << " theta_deg=" << q(18)*180.0/M_PI << " phi=" << q(19)
                << " floor=" << static_cast<int>(q(20)) << std::endl;
      }
      if (tg[0] > 0) {
        const Real cf = pmy_pack->pmesh->cfl_no;
        std::cout << (tg[0] == 1 ? "DTREG polar" : "DTREG rest ") << " cycle="
                  << pmy_pack->pmesh->ncycle << " npol=" << tg[6]
                  << " cfl*cell_dt=" << cf*q(15) << " row=" << tg[5]
                  << " theta_deg=" << q(18)*180.0/M_PI << " phi=" << q(19)
                  << " r=" << q(0) << " dir=" << static_cast<int>(q(6))
                  << " rho=" << q(1) << " T=" << q(2) << " v=" << q(7)
                  << " cgas=" << q(5) << " ceff=" << q(8) << " h=" << q(9)
                  << " cfl*cell_dt_gas=" << cf*q(16) << std::endl;
      }
    }
  }

  // compute timestep for diffusion
  if (pcond != nullptr) {
    pcond->NewTimeStep(w0, peos->eos_data);
  }
  if (pvisc != nullptr) {
    pvisc->NewTimeStep(w0, peos->eos_data);
  }
  // compute source terms timestep
  if (psrc != nullptr) {
    psrc->NewTimeStep(w0, peos->eos_data);
  }

  return TaskStatus::complete;
}
} // namespace hydro
