//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file xe_tracer_test.cpp
//! \brief tests of <mhd>/xe_tracer (xetracer-1101) on a uniform MHD box with the general
//! EOS table (cgs).  B = 0, uniform rho and T, uniform velocity <problem>/vx.
//!   <problem>/z_fac  : initial z = z_fac x local Saha x_eq (uniform)
//!   <problem>/z_amp  : plus a Gaussian bump z_amp x x_eq exp(-(x-xc)^2/w^2), xc at the
//!                      box centre, w = <problem>/z_width (advection test: run one box
//!                      crossing with <mhd>/xe_alpha_fac = 0)
//! The finalize hook prints one line per cell of the first MeshBlock row (j=js,k=ks):
//!   XE i x z x_eq n_tot T dt ncycle time
//! which tst/scripts check against the backward-Euler recursion and the exact solution.

#include <cmath>
#include <cstdio>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {
void XeFinal(ParameterInput *pin, Mesh *pm);
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  pgen_final_func = XeFinal;
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "xe_tracer_test needs <mhd>" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool trc = pmbp->pmhd->xe_tracer;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  const Real dens = pin->GetOrAddReal("problem", "rho", 1.0e-8);
  const Real tkel = pin->GetOrAddReal("problem", "T", 2500.0);
  const Real vx = pin->GetOrAddReal("problem", "vx", 0.0);
  const Real zfac = pin->GetOrAddReal("problem", "z_fac", 10.0);
  const Real zamp = pin->GetOrAddReal("problem", "z_amp", 0.0);
  const Real zw = pin->GetOrAddReal("problem", "z_width", 0.1);
  const Real xc = 0.5*(pmy_mesh_->mesh_size.x1min + pmy_mesh_->mesh_size.x1max);
  const Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  const EOS_Data eos = pmbp->pmhd->peos->eos_data;
  const Real tcode = tkel/eos.temp_cgs;
  const Real eint = eos.EnergyFromTemperature(dens, tcode);
  const Real xeq = eos.ElectronFractionKelvin(dens, tkel);
  const int iz = trc ? pmbp->pmhd->ixe : 0;
  auto w0_ = pmbp->pmhd->w0;
  par_for("pgen_xe", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    w0_(m,IDN,k,j,i) = dens;
    w0_(m,IVX,k,j,i) = vx;
    w0_(m,IVY,k,j,i) = 0.0;
    w0_(m,IVZ,k,j,i) = 0.0;
    w0_(m,IEN,k,j,i) = eint;
    Real dx = (x - xc)/(zw*lx);
    if (trc) w0_(m,iz,k,j,i) = xeq*(zfac + zamp*exp(-dx*dx));
  });
  auto &b0 = pmbp->pmhd->b0;
  Kokkos::deep_copy(b0.x1f, 0.0);
  Kokkos::deep_copy(b0.x2f, 0.0);
  Kokkos::deep_copy(b0.x3f, 0.0);
  Kokkos::deep_copy(pmbp->pmhd->bcc0, 0.0);
  // the test sets z itself: no Saha initialisation at the first ConToPrim
  pmbp->pmhd->xe_init_pending = false;
  pmbp->pmhd->peos->PrimToCons(w0_, pmbp->pmhd->bcc0, pmbp->pmhd->u0, is, ie, js, je,
                               ks, ke);
  if (global_variable::my_rank == 0) {
    std::printf("XEINIT rho %.17e T %.17e tcode %.17e x_eq %.17e\n", dens, tkel, tcode,
                xeq);
  }
}

namespace {
void XeFinal(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  if (!pmbp->pmhd->xe_tracer) return;
  const int iz = pmbp->pmhd->ixe;
  const EOS_Data eos = pmbp->pmhd->peos->eos_data;
  auto w0h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pmbp->pmhd->w0);
  auto th = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pmbp->pmhd->wtemp);
  auto &size = pmbp->pmb->mb_size;
  if (global_variable::my_rank != 0) return;
  const Real kb = 1.380649e-16;
  for (int i=is; i<=ie; ++i) {
    Real x = CellCenterX(i-is, indcs.nx1, size.h_view(0).x1min, size.h_view(0).x1max);
    const Real rho = w0h(0,IDN,ks,js,i);
    const Real tc = th(0,ks,js,i);
    const Real tk = tc*eos.temp_cgs;
    const Real mu = eos.MeanMolecularWeight(rho, w0h(0,IEN,ks,js,i), tc);
    const Real ntot = rho*eos.pres_cgs/(mu*eos.temp_cgs*kb);
    std::printf("XE %d %.17e %.17e %.17e %.17e %.17e %.17e %d %.17e\n", i-is, x,
                w0h(0,iz,ks,js,i), eos.ElectronFractionKelvin(rho, tk), ntot, tk, pm->dt,
                pm->ncycle, pm->time);
  }
}
}  // namespace
