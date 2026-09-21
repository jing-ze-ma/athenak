//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_update.cpp
//! \brief explicit stage update of the M1 moments.
//!
//! The stage weights are the module's OWN (PD-ARS with S == 0, i.e. Heun / SSP-RK2;
//! design sect. 5), not the Driver's: the module is sub-cycled and the substep length
//! is dt_sub = dt_mesh/N_sub, never pmesh->dt.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::Update

TaskStatus RadiationM1::Update(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  Real g0 = gam0[stage-1];
  Real g1 = gam1[stage-1];
  Real beta_dt = (beta[stage-1])*dt_sub;

  auto u0_ = u0;
  auto u1_ = u1;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto &mbsize = pmy_pack->pmb->mb_size;

  par_for("m1_update", DevExeSpace(), 0, nmb1, 0, M1_NVAR-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    Real divf = (flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i))/mbsize.d_view(m).dx1;
    if (multi_d) {
      divf += (flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i))/mbsize.d_view(m).dx2;
    }
    if (three_d) {
      divf += (flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i))/mbsize.d_view(m).dx3;
    }
    u0_(m,n,k,j,i) = g0*u0_(m,n,k,j,i) + g1*u1_(m,n,k,j,i) - beta_dt*divf;
  });

  return TaskStatus::complete;
}

} // namespace radm1
