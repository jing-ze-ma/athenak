//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_xe_tracer.cpp
//! \brief advected electron fraction, <mhd>/xe_tracer (xetracer-1101); see xe_tracer.hpp
//! for the model.  The source runs once per cycle in "before_timeintegrator" (after the
//! RT half step, before SaveMHDState), over EVERY cell including ghosts (pointwise, so a
//! ghost that mirrors an interior cell gets the identical value), then the radial ghosts
//! at physical x1 boundaries are made zero-gradient in z, and eta (Ohmic/ambipolar) is
//! recomputed from the relaxed z.

#include <iostream>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/resistivity.hpp"
#include "pgen/pgen.hpp"
#include "mhd/mhd.hpp"
#include "mhd/xe_tracer.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::XeTracerSource
//! \brief recombination source over the cycle dt, then eta from the new z

TaskStatus MHD::XeTracerSource(Driver *pdrive, int stage) {
  XeTracerUpdate(pmy_pack->pmesh->dt, false);
  if (presist != nullptr) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int ng = indcs.ng;
    const int n1m1 = indcs.nx1 + 2*ng - 1;
    const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    const Real rgas = pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas;
    if (presist->OhmicNeedsUpdate()) {
      presist->SetResistivity(w0, peos->eos_data, rgas, presist->eta_b,
                              0, n1m1, 0, n2m1, 0, n3m1);
    }
    if (presist->AmbipolarNeedsUpdate()) {
      presist->SetAmbipolar(w0, peos->eos_data, rgas, 0, n1m1, 0, n2m1, 0, n3m1);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::XeTracerUpdate
//! \brief init = true: z = local Saha x_eq.  Otherwise one backward-Euler recombination
//! step of n_e = z n_tot over dt (code units), capped at max(x_eq, gas-phase donors).
//! Writes w0 and u0 (= rho z) of the tracer in every cell.

void MHD::XeTracerUpdate(Real dt, bool init) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = peos->eos_data;
  const bool gen = eos.IsGeneral();
  const Real gm1 = eos.gamma - 1.0;
  const Real rgas = pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas;
  if (!gen && !(rgas > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/xe_tracer on an ideal gas needs <problem>/Rgas"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const xetracer::XeTracerParams par = xe_par;
  const Real dt_s = dt*par.time_cgs;
  const Real nh_fac = par.xh/(eos_cgs::a_hyd*eos_cgs::m_u);
  const int iz = ixe;
  auto u0_ = u0;
  auto w0_ = w0;
  auto &wtemp_ = wtemp;
  par_for("xe_tracer", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real kb = 1.380649e-16;
    const Real rho = w0_(m,IDN,k,j,i);
    Real t, ntot, rhoc;
    if (gen) {
      const Real tcode = wtemp_(m,k,j,i);
      t = tcode*eos.temp_cgs;
      const Real mu = eos.MeanMolecularWeight(rho, w0_(m,IEN,k,j,i), tcode);
      ntot = rho*eos.pres_cgs/(mu*eos.temp_cgs*kb);
      rhoc = rho*eos.dens_cgs;
    } else {
      t = gm1*w0_(m,IEN,k,j,i)/rgas/rho;
      ntot = rho*rgas/kb;
      rhoc = rho;
    }
    const Real xeq = eos.ElectronFractionKelvin(rho, t);
    Real z;
    if (init) {
      z = xeq;
    } else {
      const Real neq = xeq*ntot;
      const Real a = par.alpha_fac*xetracer::AlphaRR(t)*dt_s;
      Real n = xetracer::Relax(w0_(m,iz,k,j,i)*ntot, neq, a);
      const Real ncap = fmax(neq, xetracer::DonorCap(par, t, ntot*kb*t, nh_fac*rhoc));
      n = fmin(n, ncap);
      z = n/ntot;
    }
    w0_(m,iz,k,j,i) = z;
    u0_(m,iz,k,j,i) = u0_(m,IDN,k,j,i)*z;
  });

  // zero-gradient z in the radial ghosts at physical x1 boundaries (the dhj inner wall
  // and top): the ghost T there is a frozen or extrapolated state, and its local Saha
  // value would feed a spurious slope into the PLM reconstruction of the top cells
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const int is = indcs.is, ie = indcs.ie;
  par_for("xe_tracer_x1bc", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(int m, int k, int j) {
    const BoundaryFlag fi = mb_bcs.d_view(m, BoundaryFace::inner_x1);
    if (fi != BoundaryFlag::block && fi != BoundaryFlag::periodic &&
        fi != BoundaryFlag::panel) {
      const Real z = w0_(m,iz,k,j,is);
      for (int g=1; g<=ng; ++g) {
        w0_(m,iz,k,j,is-g) = z;
        u0_(m,iz,k,j,is-g) = u0_(m,IDN,k,j,is-g)*z;
      }
    }
    const BoundaryFlag fo = mb_bcs.d_view(m, BoundaryFace::outer_x1);
    if (fo != BoundaryFlag::block && fo != BoundaryFlag::periodic &&
        fo != BoundaryFlag::panel) {
      const Real z = w0_(m,iz,k,j,ie);
      for (int g=1; g<=ng; ++g) {
        w0_(m,iz,k,j,ie+g) = z;
        u0_(m,iz,k,j,ie+g) = u0_(m,IDN,k,j,ie+g)*z;
      }
    }
  });
  return;
}

}  // namespace mhd
