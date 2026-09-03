//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_update.cpp
//! \brief Performs explicit update of MHD conserved variables (u0) for each stage of the
//! SSP RK integrators (e.g. RK1, RK2, RK3) implemented in AthenaK, using weighted average
//! and partial time update of flux divergence. Source terms are added in the
//! MHDSrcTerms() function.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn  void MHD::Update
//  \brief Explicit RK update including flux divergence terms

TaskStatus MHD::RKUpdate(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  Real &gam0 = pdriver->gam0[stage-1];
  Real &gam1 = pdriver->gam1[stage-1];
  Real beta_dt = (pdriver->beta[stage-1])*(pmy_pack->pmesh->dt);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  int nv1 = nmhd + nscalars - 1;
  auto u0_ = u0;
  auto u1_ = u1;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  // cubed-sphere RHS-split diagnostic: drop the flux divergence entirely
  const Real dfac = cs_diag_no_divf ? 0.0 : 1.0;
    
  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &volume = pmy_pack->pcoord->volume;
  auto &area1 = pmy_pack->pcoord->area.x1f;
  auto &area2 = pmy_pack->pcoord->area.x2f;
  auto &area3 = pmy_pack->pcoord->area.x3f;

  // hierarchical parallel loop that updates conserved variables to intermediate step
  // using weights and fractional time step appropriate to stages of time-integrator used
  // Vector inner loop used for good performance on cpus
  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1);

  par_for_outer("mhd_update",DevExeSpace(),scr_size,scr_level,0,nmb1,0,nv1,ks,ke,js,je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int n, const int k, const int j) {
    ScrArray1D<Real> divf(member.team_scratch(scr_level), ncells1);
      
    if (use_cubed_sphere || use_spherical_polar) {
          
      // compute dF1/dx1
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) = (flx1(m,n,k,j,i+1)*area1(m,k,j,i+1) - flx1(m,n,k,j,i)*area1(m,k,j,i))/volume(m,k,j,i);
      });
      member.team_barrier();

      // Add dF2/dx2
      // Fluxes must be summed in pairs to symmetrize round-off error in each dir
      if (multi_d) {
        par_for_inner(member, is, ie, [&](const int i) {
          divf(i) += (flx2(m,n,k,j+1,i)*area2(m,k,j+1,i) - flx2(m,n,k,j,i)*area2(m,k,j,i))/volume(m,k,j,i);
      });
      member.team_barrier();
      }

      // Add dF3/dx3
      // Fluxes must be summed in pairs to symmetrize round-off error in each dir
      if (three_d) {
        par_for_inner(member, is, ie, [&](const int i) {
          divf(i) += (flx3(m,n,k+1,j,i)*area3(m,k+1,j,i) - flx3(m,n,k,j,i)*area3(m,k,j,i))/volume(m,k,j,i);
        });
        member.team_barrier();
      }
            
    } else {

    // compute dF1/dx1
    par_for_inner(member, is, ie, [&](const int i) {
      divf(i) = (flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i))/mbsize.d_view(m).dx1;
    });
    member.team_barrier();

    // Add dF2/dx2
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (multi_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += (flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i))/mbsize.d_view(m).dx2;
      });
      member.team_barrier();
    }

    // Add dF3/dx3
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (three_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += (flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i))/mbsize.d_view(m).dx3;
      });
      member.team_barrier();
    }
        
    }

    par_for_inner(member, is, ie, [&](const int i) {
      u0_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i)
                       - dfac*beta_dt*divf(i);
    });
  });
  return TaskStatus::complete;
}
} // namespace mhd
