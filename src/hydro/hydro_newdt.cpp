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

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro.hpp"
#include "diffusion/conduction.hpp"
#include "diffusion/viscosity.hpp"
#include "srcterms/srcterms.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

namespace hydro {

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
  // module's transport opacity, opac(M1_OP_T)):
  //   c_d^2 = c_g^2 + (1 - exp(-tau_d)) (c_eq,d^2 - c_g^2),
  // so optically thin cells keep the gas speed.  The fully coupled speed is the LARGEST
  // acoustic speed of the medium (radiative diffusion at the grid scale only lowers it),
  // so the bound is conservative.  E, F, the tensor and opac are the M1 state of the
  // last radiation step (the module runs after the hydro stages): a lag of one step.
  // Off: no change.
  const bool rss_ = rad_signal_speed;
  DvceArray5D<Real> erad_, kopc_, tten_, vcel_;
  int rss_mode = 0;          // 0 isotropic E/3, 1 M1Chi, 2 tau_ten, 3 vet_sc, 4 vet full
  int rss_kind = 0;
  Real rss_cl = 1.0;
  if (rss_) {
    auto *prm = pmy_pack->pradm1;
    if (prm == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/rad_signal_speed = true needs a <rad_m1> block"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    erad_ = prm->u0;
    kopc_ = prm->opac;
    rss_cl = prm->c_light;
    rss_kind = prm->chi_kind;
    if (prm->tau_closure && prm->tau_ten.extent(0) > 0) {
      rss_mode = 2;
      tten_ = prm->tau_ten;
    } else if (prm->vet_sc && prm->vet_cell.extent(0) > 0) {
      rss_mode = prm->vet_full ? 4 : 3;
      vcel_ = prm->vet_cell;
    } else if (!prm->eddington && !prm->tau_closure && !prm->vet_sc) {
      rss_mode = 1;
    }
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
        min_dt1 = fmin((dx1_(m,k,j,i)/fabs(w0_(m,IVX,k,j,i))), min_dt1);
        min_dt2 = fmin((dx2_(m,k,j,i)/fabs(w0_(m,IVY,k,j,i))), min_dt2);
        min_dt3 = fmin((dx3_(m,k,j,i)/fabs(w0_(m,IVZ,k,j,i))), min_dt3);
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
          const Real d = w0_(m,IDN,k,j,i);
          const Real er = fmax(erad_(m,radm1::M1_E,k,j,i), 0.0);
          // diagonal of the Eddington tensor, D_dd = P_dd/E
          Real dd[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
          if (rss_mode == 4) {
            for (int a = 0; a < 3; ++a) dd[a] = vcel_(m,radm1::M1_VET_D11+a,k,j,i);
          } else if (rss_mode >= 1) {
            Real chi, n1, n2, n3;
            if (rss_mode == 2) {
              chi = tten_(m,0,k,j,i); n1 = tten_(m,1,k,j,i);
              n2 = tten_(m,2,k,j,i); n3 = tten_(m,3,k,j,i);
            } else if (rss_mode == 3) {
              chi = vcel_(m,radm1::M1_VET_CHI,k,j,i);
              n1 = vcel_(m,radm1::M1_VET_N1,k,j,i);
              n2 = vcel_(m,radm1::M1_VET_N1+1,k,j,i);
              n3 = vcel_(m,radm1::M1_VET_N1+2,k,j,i);
            } else {
              const Real f1 = erad_(m,radm1::M1_F1,k,j,i);
              const Real f2 = erad_(m,radm1::M1_F2,k,j,i);
              const Real f3 = erad_(m,radm1::M1_F3,k,j,i);
              const Real fn = sqrt(f1*f1 + f2*f2 + f3*f3);
              const Real inv = (fn > 0.0) ? 1.0/fn : 0.0;
              chi = (er > 0.0) ? radm1::M1Chi(fn/(rss_cl*er), rss_kind) : 1.0/3.0;
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
          const Real kt = fmax(kopc_(m,radm1::M1_OP_T,k,j,i), 0.0);
          Real h[3];
          if (use_cubed_sphere || use_spherical_polar) {
            h[0] = dx1_(m,k,j,i); h[1] = dx2_(m,k,j,i); h[2] = dx3_(m,k,j,i);
          } else {
            h[0] = mbsize.d_view(m).dx1; h[1] = mbsize.d_view(m).dx2;
            h[2] = mbsize.d_view(m).dx3;
          }
          Real pg = 0.0, gm1 = 0.0;
          const bool ideal = (eos.is_ideal && !eos.IsGeneral());
          if (ideal) {
            pg = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
            gm1 = eos.gamma - 1.0;
          }
          Real csd[3];
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
            const Real dc2 = fmax(ceq2 - cg2, 0.0);
            csd[a] = sqrt(cg2 + dc2*(1.0 - exp(-kt*h[a])));
          }
          cs1 = csd[0]; cs2 = csd[1]; cs3 = csd[2];
        }
        max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cs1;
        max_dv2 = fabs(w0_(m,IVY,k,j,i))
                 + (cs_ ? cs2/sncell_(m,k,j) : cs2);
        max_dv3 = fabs(w0_(m,IVZ,k,j,i))
                 + (cs_ ? cs3/sncell_(m,k,j) : cs3);
      }
      Real cell_dt;
      if (use_cubed_sphere || use_spherical_polar) {
        min_dt1 = fmin((dx1_(m,k,j,i)/max_dv1), min_dt1);
        min_dt2 = fmin((dx2_(m,k,j,i)/max_dv2), min_dt2);
        min_dt3 = fmin((dx3_(m,k,j,i)/max_dv3), min_dt3);
        cell_dt = dx1_(m,k,j,i)/max_dv1;
        if (multi_d_) cell_dt = fmin(cell_dt, dx2_(m,k,j,i)/max_dv2);
        if (three_d_) cell_dt = fmin(cell_dt, dx3_(m,k,j,i)/max_dv3);
      } else {
      min_dt1 = fmin((mbsize.d_view(m).dx1/max_dv1), min_dt1);
      min_dt2 = fmin((mbsize.d_view(m).dx2/max_dv2), min_dt2);
      min_dt3 = fmin((mbsize.d_view(m).dx3/max_dv3), min_dt3);
      cell_dt = mbsize.d_view(m).dx1/max_dv1;
      if (multi_d_) cell_dt = fmin(cell_dt, mbsize.d_view(m).dx2/max_dv2);
      if (three_d_) cell_dt = fmin(cell_dt, mbsize.d_view(m).dx3/max_dv3);
      }
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
