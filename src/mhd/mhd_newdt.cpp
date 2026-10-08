//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_newdt.cpp
//! \brief function to compute MHD timestep across all MeshBlock(s) in a MeshBlockPack

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
#include "mhd.hpp"
#include "diffusion/conduction.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/resistivity.hpp"
#include "srcterms/srcterms.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
// \!fn void MHD::NewTimeStep()
// \brief calculate the minimum timestep within a MeshBlockPack for MHD problems

TaskStatus MHD::NewTimeStep(Driver *pdriver, int stage) {
  if (stage != (pdriver->nexp_stages)) {
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
  Kokkos::ValLocScalar<Real, int> mloc;
  mloc.val = std::numeric_limits<float>::max();
  mloc.loc = -1;
  if (dt_diag.h_view.extent(0) == 0) Kokkos::realloc(dt_diag, ndtdiag);

  // capture class variables for kernel
  auto &w0_ = w0;
  auto &eos = pmy_pack->pmhd->peos->eos_data;
  // derived thermodynamic variables (general EOS only; empty view for an ideal gas)
  auto &wder_ = pmy_pack->pmhd->wder;
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
  // the kinematic branch below needs no correction.  The FAST SPEED does not scale that
  // way: cf is a physical speed, so the constraint is dx2/(|v2| + cf/sin_cell).  Leaving
  // it out overstates dt by up to 13.7% at a cube vertex (sin_cell = 0.88 there) and not
  // at all on the panel axes.  sp and Cartesian are untouched.
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  auto &sncell_ = pmy_pack->pcoord->sin_cell;
  const bool multi_d_ = pmy_pack->pmesh->multi_d;
  const bool three_d_ = pmy_pack->pmesh->three_d;

  // RADIATION SIGNAL SPEED (<mhd>/rad_signal_speed; default on with an implicit M1
  // transport, see mhd.cpp).  The MHD form of <hydro>/rad_signal_speed (hydro_newdt.cpp,
  // where the closure diagonal, the equilibrium sound speed and the optical-depth taper
  // are derived).  In the strong-coupling (equilibrium-diffusion) limit gas + radiation
  // is one fluid of total pressure P = P_g + P_r and adiabatic exponent Gamma_1 (ideal
  // gas: the Chandrasekhar mixture; general EOS: the decoupled sum).  The pressure enters
  // the ideal-MHD wave speeds only through the sound speed a^2 = Gamma_1 P/rho, so the
  // fast magnetosonic speed along direction d is the usual one with a^2 replaced by the
  // tapered radiation-modified sound speed of hydro_newdt.cpp,
  //   c_d^2  = c_g^2 + (1 - exp(-tau_d)) (c_eq,d^2 - c_g^2),
  //   cf_d^2 = (1/2) [ c_d^2 + B^2/rho + sqrt((c_d^2 + B^2/rho)^2 - 4 c_d^2 B_d^2/rho) ],
  // c_g the gas sound speed of the fast speed below (ideal: gamma p/rho; general EOS:
  // Gamma_1 p/rho; isothermal: iso_cs^2).  cf_d grows monotonically with c_d^2, so this
  // is >= the gas fast speed and -> it in thin cells (tau_d -> 0); at B = 0 it is the
  // hydro signal speed.  Off: no change (the gas fast speed, bitwise).
  const bool rss_ = rad_signal_speed;
  DvceArray5D<Real> erad_, kopc_, tten_, vcel_;
  int rss_mode = 0;          // 0 isotropic E/3, 1 M1Chi, 2 tau_ten, 3 vet_sc, 4 vet full
  int rss_kind = 0;
  Real rss_cl = 1.0;
  if (rss_) {
    auto *prm = pmy_pack->pradm1;
    if (prm == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/rad_signal_speed = true needs a <rad_m1> block"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    erad_ = prm->u0;
    kopc_ = prm->opac;
    rss_cl = prm->c_light;
    rss_kind = prm->chi_kind;
    // the same selection decides what a restart file carries (RadiationM1::RssPrime)
    rss_mode = prm->RssMode();
    if (rss_mode == 2) tten_ = prm->tau_ten;
    if (rss_mode >= 3) vcel_ = prm->vet_cell;
  }

  if (pdriver->time_evolution == TimeEvolution::kinematic) {
    // find smallest (dx/v) in each direction for advection problems
    Kokkos::parallel_reduce("MHDNudt1",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
      // compute m,k,j,i indices of thread and call function
      int m = (idx)/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      // NON-FINITE GUARD (see the MHDNudt2 kernel below): a NaN velocity gives a NaN
      // dx/|v|, which fmin drops; flag it with a negative dt so the run stops
      if (!(fabs(w0_(m,IVX,k,j,i)) + fabs(w0_(m,IVY,k,j,i)) + fabs(w0_(m,IVZ,k,j,i))
            < 1.0e300)) {
        min_dt1 = -1.0;
        return;
      }
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
    // find smallest dx/(v +/- Cf) in each direction for mhd problems
    auto &bcc0_ = bcc0;
    // the slot the Newtonian gas pressure is read from: an isothermal EOS has no energy
    // (IEN is out of range) and never uses that pressure, so it reads IDN instead
    const int ipr_ = eos.is_ideal ? static_cast<int>(IEN) : static_cast<int>(IDN);

    Kokkos::parallel_reduce("MHDNudt2",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
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

      // timestep in GR MHD
      if (is_general_relativistic_ || is_dynamical_relativistic_) {
        max_dv1 = 1.0;
        max_dv2 = 1.0;
        max_dv3 = 1.0;
      // timestep in SR MHD
      } else if (is_special_relativistic_) {
        Real &wd = w0_(m,IDN,k,j,i);
        Real &ux = w0_(m,IVX,k,j,i);
        Real &uy = w0_(m,IVY,k,j,i);
        Real &uz = w0_(m,IVZ,k,j,i);
        Real &bcc1 = bcc0_(m,IBX,k,j,i);
        Real &bcc2 = bcc0_(m,IBY,k,j,i);
        Real &bcc3 = bcc0_(m,IBZ,k,j,i);

        Real v2 = SQR(ux) + SQR(uy) + SQR(uz);
        Real lor = sqrt(1.0 + v2);
        // FIXME ERM: Ideal fluid for now
        Real p = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
        // Calculate 4-magnetic field in left state
        Real b_0 = bcc1*ux + bcc2*uy + bcc3*uz;
        Real b_1 = (bcc1 + b_0 * ux) / lor;
        Real b_2 = (bcc2 + b_0 * uy) / lor;
        Real b_3 = (bcc3 + b_0 * uz) / lor;
        Real b_sq = -SQR(b_0) + SQR(b_1) + SQR(b_2) + SQR(b_3);

        Real lm, lp;
        eos.IdealSRMHDFastSpeeds(wd, p, ux, lor, b_sq, lp, lm);
        max_dv1 = fmax(fabs(lm), lp);

        eos.IdealSRMHDFastSpeeds(wd, p, uy, lor, b_sq, lp, lm);
        max_dv2 = fmax(fabs(lm), lp);

        eos.IdealSRMHDFastSpeeds(wd, p, uz, lor, b_sq, lp, lm);
        max_dv3 = fmax(fabs(lm), lp);
      // timestep in Newtonian MHD
      } else {
        Real &w_d = w0_(m,IDN,k,j,i);
        Real &w_bx = bcc0_(m,IBX,k,j,i);
        Real &w_by = bcc0_(m,IBY,k,j,i);
        Real &w_bz = bcc0_(m,IBZ,k,j,i);
        Real cf;
        Real p = eos.IdealGasPressure(w0_(m,ipr_,k,j,i));
        // For a general EOS the pressure and Gamma_1 were evaluated once per cell in
        // ConsToPrim. Using them here matters: with ionization or radiation pressure the
        // ideal-gas expression gives the wrong fast speed, hence the wrong CFL timestep.
        Real pg = 0.0, g1 = 0.0;
        if (eos.IsGeneral()) {
          pg = wder_(m,IDPR,k,j,i);
          g1 = wder_(m,IDG1,k,j,i);
        }

        if (eos.IsGeneral()) {
          cf = eos.FastSpeedFromP(w_d, pg, g1, w_bx, w_by, w_bz);
        } else if (eos.is_ideal) {
          cf = eos.IdealMHDFastSpeed(w_d, p, w_bx, w_by, w_bz);
        } else {
          cf = eos.IdealMHDFastSpeed(w_d, w_bx, w_by, w_bz);
        }
        max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cf;

        if (eos.IsGeneral()) {
          cf = eos.FastSpeedFromP(w_d, pg, g1, w_by, w_bz, w_bx);
        } else if (eos.is_ideal) {
          cf = eos.IdealMHDFastSpeed(w_d, p, w_by, w_bz, w_bx);
        } else {
          cf = eos.IdealMHDFastSpeed(w_d, w_by, w_bz, w_bx);
        }
        max_dv2 = fabs(w0_(m,IVY,k,j,i))
                 + (cs_ ? cf/sncell_(m,k,j) : cf);

        if (eos.IsGeneral()) {
          cf = eos.FastSpeedFromP(w_d, pg, g1, w_bz, w_bx, w_by);
        } else if (eos.is_ideal) {
          cf = eos.IdealMHDFastSpeed(w_d, p, w_bz, w_bx, w_by);
        } else {
          cf = eos.IdealMHDFastSpeed(w_d, w_bz, w_bx, w_by);
        }
        max_dv3 = fabs(w0_(m,IVZ,k,j,i))
                 + (cs_ ? cf/sncell_(m,k,j) : cf);

        if (rss_) {
          // <mhd>/rad_signal_speed: the same closure diagonal, equilibrium sound speed
          // and taper as hydro_newdt.cpp, then the fast speed with a^2 -> c_d^2
          const Real d = w_d;
          const Real er = fmax(erad_(m,radm1::M1_E,k,j,i), 0.0);
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
          // the gas sound speed squared of the fast speeds above
          Real cg2;
          if (eos.IsGeneral()) {
            cg2 = g1*pg/d;
          } else if (eos.is_ideal) {
            cg2 = eos.gamma*p/d;
          } else {
            cg2 = eos.iso_cs*eos.iso_cs;
          }
          const Real kt = fmax(kopc_(m,radm1::M1_OP_T,k,j,i), 0.0);
          Real h[3];
          if (use_cubed_sphere || use_spherical_polar) {
            h[0] = dx1_(m,k,j,i); h[1] = dx2_(m,k,j,i); h[2] = dx3_(m,k,j,i);
          } else {
            h[0] = mbsize.d_view(m).dx1; h[1] = mbsize.d_view(m).dx2;
            h[2] = mbsize.d_view(m).dx3;
          }
          const bool ideal = (eos.is_ideal && !eos.IsGeneral());
          const Real gm1 = ideal ? (eos.gamma - 1.0) : 0.0;
          const Real bb[3] = {w_bx, w_by, w_bz};
          Real cfd[3];
          for (int a = 0; a < 3; ++a) {
            const Real pr = fmax(dd[a], 0.0)*er;
            Real ceq2;
            if (ideal) {
              const Real ptot = p + pr;
              const Real beta = p/ptot;
              const Real gam1 = beta + SQR(4.0 - 3.0*beta)*gm1
                                       /(beta + 12.0*gm1*(1.0 - beta));
              ceq2 = gam1*ptot/d;
            } else {
              ceq2 = cg2 + (4.0/3.0)*pr/d;
            }
            const Real dc2 = fmax(ceq2 - cg2, 0.0);
            const Real asq = d*(cg2 + dc2*(1.0 - exp(-kt*h[a])));   // rho c_d^2
            const Real bn2 = bb[a]*bb[a];
            const Real ct2 = bb[(a+1)%3]*bb[(a+1)%3] + bb[(a+2)%3]*bb[(a+2)%3];
            const Real qsq = bn2 + ct2 + asq;
            const Real tmp = bn2 + ct2 - asq;
            cfd[a] = sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2))/d);
          }
          max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cfd[0];
          max_dv2 = fabs(w0_(m,IVY,k,j,i))
                   + (cs_ ? cfd[1]/sncell_(m,k,j) : cfd[1]);
          max_dv3 = fabs(w0_(m,IVZ,k,j,i))
                   + (cs_ ? cfd[2]/sncell_(m,k,j) : cfd[2]);
        }
      }

      // NON-FINITE GUARD (the MHD port of hydro_newdt.cpp): a NaN/inf state gives a NaN
      // dx/v, which fmin DROPS, so dt would be set by the finite cells alone and could
      // keep doubling (BSG production 10-01: the whole domain NaN, dt doubled per step
      // to tlim, rc 0).  Flag the cell with a negative dt instead; Mesh::NewTimeStep
      // stops the run on dt <= 0 and reports the cell.
      if (!(max_dv1 + max_dv2 + max_dv3 < 1.0e300)) {
        min_dt1 = -1.0;
        mres.val = -1.0;
        mres.loc = idx;
        return;
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
       Kokkos::MinLoc<Real, int>(mloc));
  }

  // compute minimum of dt1/dt2/dt3 for 1D/2D/3D problems
  dtnew = dt1;
  if (pmy_pack->pmesh->multi_d) { dtnew = std::min(dtnew, dt2); }
  if (pmy_pack->pmesh->three_d) { dtnew = std::min(dtnew, dt3); }

  // decode the winning cell, and -- only when dt has just collapsed -- its state
  if (mloc.loc >= 0 && mloc.loc < nmkji) {
    dtnew_m = (mloc.loc)/nkji;
    dtnew_k = (mloc.loc - dtnew_m*nkji)/nji + ks;
    dtnew_j = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji)/nx1 + js;
    dtnew_i = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji - (dtnew_j-js)*nx1) + is;
  } else {
    dtnew_m = dtnew_k = dtnew_j = dtnew_i = -1;
  }
  dt_diag_valid = false;
  if (dtnew_m >= 0 && (dtnew_prev < 0.0 || dtnew < 0.25*dtnew_prev)) {
    auto dd = dt_diag;
    const int dm = dtnew_m, dk = dtnew_k, dj = dtnew_j, di = dtnew_i;
    auto &x1v_ = pmy_pack->pcoord->x1v;
    auto &wtemp_ = pmy_pack->pmhd->wtemp;
    auto &bccd_ = bcc0;
    auto eos_ = eos;
    const bool curvilinear_ = (use_cubed_sphere || use_spherical_polar);
    const bool relativistic_ = (is_special_relativistic_ || is_general_relativistic_ ||
                                is_dynamical_relativistic_);
    const int is_ = is, nx1_ = nx1;
    par_for("mhd_dtdiag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      const Real d = w0_(dm,IDN,dk,dj,di);
      const Real bx = bccd_(dm,IBX,dk,dj,di);
      const Real by = bccd_(dm,IBY,dk,dj,di);
      const Real bz = bccd_(dm,IBZ,dk,dj,di);
      const Real bsq = SQR(bx) + SQR(by) + SQR(bz);
      Real pr, cf;
      if (eos_.IsGeneral()) {
        pr = wder_(dm,IDPR,dk,dj,di);
        cf = eos_.FastSpeedFromP(d, pr, wder_(dm,IDG1,dk,dj,di), bx, by, bz);
      } else if (eos_.is_ideal) {
        pr = eos_.IdealGasPressure(w0_(dm,IEN,dk,dj,di));
        cf = eos_.IdealMHDFastSpeed(d, pr, bx, by, bz);
      } else {
        pr = d*SQR(eos_.iso_cs);
        cf = eos_.IdealMHDFastSpeed(d, bx, by, bz);
      }
      // the fast speed above is the NEWTONIAN one; in SR/GR the signal speed is set by
      // the light cone in the reduction, so reporting it here would only mislead
      if (relativistic_) cf = 0.0;
      // x1v is allocated only on spherical-polar/cubed-sphere grids; on a Cartesian
      // grid it is a (1,1) dummy, so build the cell-centre position from the block size
      dd.d_view(0) = curvilinear_ ? x1v_(dm,di)
                   : CellCenterX(di-is_, nx1_, mbsize.d_view(dm).x1min,
                                 mbsize.d_view(dm).x1max);
      dd.d_view(1) = d;
      dd.d_view(2) = eos_.IsGeneral() ? wtemp_(dm,dk,dj,di) : (pr/d);
      dd.d_view(3) = pr;
      dd.d_view(4) = bsq;
      // beta is the number that usually explains an MHD collapse; a field-free cell
      // has no beta, so report a sentinel rather than a division by zero
      dd.d_view(5) = (bsq > 0.0) ? (2.0*pr/bsq) : -1.0;
      dd.d_view(6) = w0_(dm,IVX,dk,dj,di);
      dd.d_view(7) = w0_(dm,IVY,dk,dj,di);
      dd.d_view(8) = w0_(dm,IVZ,dk,dj,di);
      dd.d_view(9) = cf;
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
  if (presist != nullptr) {
    presist->NewTimeStep(w0, peos->eos_data);
  }
  // compute source terms timestep
  if (psrc != nullptr) {
    psrc->NewTimeStep(w0, peos->eos_data);
  }

  return TaskStatus::complete;
}
} // namespace mhd
