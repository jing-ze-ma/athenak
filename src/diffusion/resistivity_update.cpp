//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_update.cpp
//! \brief Update of the MHD energy u0(IEN) for each super-stage of the RKG
//! super-time-stepping of the resistive terms (second-order Runge-Kutta-Gegenbauer).

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "diffusion/resistivity.hpp"

//----------------------------------------------------------------------------------------
//! \fn TaskStatus Resistivity::RKUpdate
//! \brief RKG recurrence for the energy: stage 1 u_1 = u_0 - dt f0 div F_0, stage j
//! u_j = mu u_{j-1} + nu u_{j-2} + (1-mu-nu) u_0 - dt (mu~ div F_{j-1} + gamma~ div F_0).

TaskStatus Resistivity::RKUpdate(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto dt = pmy_pack->pmesh->dt*rkg_dtfac;   // dt/2 for a Strang half
  auto u0_ = pmy_pack->pmhd->u0;
  auto u2_ = u2;
  auto u_ideal_ = u_ideal;
  auto flx1 = pmy_pack->pmhd->uflx.x1f;
  auto flx2 = pmy_pack->pmhd->uflx.x2f;
  auto flx3 = pmy_pack->pmhd->uflx.x3f;
  auto flx1_ideal = uflx_ideal.x1f;
  auto flx2_ideal = uflx_ideal.x2f;
  auto flx3_ideal = uflx_ideal.x3f;
  auto &mbsize = pmy_pack->pmb->mb_size;

  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &volume = pmy_pack->pcoord->volume;
  auto &area1 = pmy_pack->pcoord->area.x1f;
  auto &area2 = pmy_pack->pcoord->area.x2f;
  auto &area3 = pmy_pack->pcoord->area.x3f;

  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1);

  int n = (IEN);
  Real fjm1;
  Real f0;
  if (stage == 1) {
    fjm1 = 0.0;
    f0 = bjm1*w1;
  } else {
    fjm1 = mut;
    f0 = gat;
  }
  // local copies of the RKG weights, so the kernels do not capture `this`
  const Real mu_ = mu;
  const Real nu_ = nu;
  // with rkg_lean the energy slot of u2 (u_{j-2}, read only here, in the active cells)
  // is rotated in place: the update saves the u0 it overwrites for the next stage
  const bool lean_ = rkg_lean;

  // the energy flux through the x1 faces is and ie+1 of every column, followed through
  // the same recurrence: E_1 = f0 F_0, E_j = mu E_{j-1} + nu E_{j-2} + fjm1 F_j + f0 F_0
  if (track_x1_bflx) {
    const int nmb = nmb1 + 1;
    const int nsl = rkg_strang ? 8 : 6;   // slots 6|7 keep the first Strang half
    if (static_cast<int>(rkg_bflx.extent(0)) != nmb ||
        static_cast<int>(rkg_bflx.extent(1)) != nsl) {
      Kokkos::realloc(rkg_bflx, nmb, nsl, indcs.nx3 + 2*indcs.ng, indcs.nx2 + 2*indcs.ng);
    }
    auto bf = rkg_bflx;
    const bool first = (stage == 1);
    const Real fj_ = fjm1, f0_ = f0;
    par_for("rkg_bflx", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      for (int s2=0; s2<2; ++s2) {
        const Real f = flx1(m,n,k,j,(s2 == 0) ? is : ie+1);
        if (first) {
          bf(m,s2,k,j) = f;
          bf(m,2+s2,k,j) = 0.0;
          bf(m,4+s2,k,j) = f0_*f;
        } else {
          const Real e2 = bf(m,2+s2,k,j), e1 = bf(m,4+s2,k,j);
          bf(m,2+s2,k,j) = e1;
          bf(m,4+s2,k,j) = mu_*e1 + nu_*e2 + fj_*f + f0_*bf(m,s2,k,j);
        }
      }
    });
  }

  // curvilinear meshes: one cell per thread, flux differences weighted by face areas
  if (use_cubed_sphere || use_spherical_polar) {
    const bool md_ = multi_d, td_ = three_d;
    par_for("resmhd_update", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real divf = ((fjm1*flx1(m,n,k,j,i+1)+f0*flx1_ideal(m,n,k,j,i+1))*area1(m,k,j,i+1)
                 - (fjm1*flx1(m,n,k,j,i)+f0*flx1_ideal(m,n,k,j,i))*area1(m,k,j,i))
                 /volume(m,k,j,i);
      if (md_) {
        divf += ((fjm1*flx2(m,n,k,j+1,i)+f0*flx2_ideal(m,n,k,j+1,i))*area2(m,k,j+1,i)
               - (fjm1*flx2(m,n,k,j,i)+f0*flx2_ideal(m,n,k,j,i))*area2(m,k,j,i))
               /volume(m,k,j,i);
      }
      if (td_) {
        divf += ((fjm1*flx3(m,n,k+1,j,i)+f0*flx3_ideal(m,n,k+1,j,i))*area3(m,k+1,j,i)
               - (fjm1*flx3(m,n,k,j,i)+f0*flx3_ideal(m,n,k,j,i))*area3(m,k,j,i))
               /volume(m,k,j,i);
      }
      if (stage == 1) {
        u0_(m,n,k,j,i) = u0_(m,n,k,j,i) - dt*divf;
      } else {
        const Real uold_ = u0_(m,n,k,j,i);
        u0_(m,n,k,j,i) = mu_*uold_ + nu_*u2_(m,n,k,j,i)
                         + (1.0-mu_-nu_)*u_ideal_(m,n,k,j,i) - dt*divf;
        if (lean_) u2_(m,n,k,j,i) = uold_;
      }
    });
    return TaskStatus::complete;
  }

  // Cartesian: hierarchical parallel loop, vector inner loop for good performance on cpus
  scr_level = TeamScratchLevel(scr_size, 0);
  par_for_outer("resmhd_update",DevExeSpace(),scr_size,scr_level,0,nmb1,ks,ke,js,je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> divf(member.team_scratch(scr_level), ncells1);

    // compute dF1/dx1
    par_for_inner(member, is, ie, [&](const int i) {
      divf(i) = ((fjm1*flx1(m,n,k,j,i+1)+f0*flx1_ideal(m,n,k,j,i+1))
               - (fjm1*flx1(m,n,k,j,i)+f0*flx1_ideal(m,n,k,j,i)))/mbsize.d_view(m).dx1;
    });
    member.team_barrier();

    // Add dF2/dx2
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (multi_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += ((fjm1*flx2(m,n,k,j+1,i)+f0*flx2_ideal(m,n,k,j+1,i))
                  - (fjm1*flx2(m,n,k,j,i)+f0*flx2_ideal(m,n,k,j,i)))/mbsize.d_view(m).dx2;
      });
      member.team_barrier();
    }

    // Add dF3/dx3
    // Fluxes must be summed in pairs to symmetrize round-off error in each dir
    if (three_d) {
      par_for_inner(member, is, ie, [&](const int i) {
        divf(i) += ((fjm1*flx3(m,n,k+1,j,i)+f0*flx3_ideal(m,n,k+1,j,i))
                  - (fjm1*flx3(m,n,k,j,i)+f0*flx3_ideal(m,n,k,j,i)))/mbsize.d_view(m).dx3;
      });
      member.team_barrier();
    }

    par_for_inner(member, is, ie, [&](const int i) {
      if (stage == 1) {
        u0_(m,n,k,j,i) = u0_(m,n,k,j,i) - dt*divf(i);
      } else {
        const Real uold_ = u0_(m,n,k,j,i);
        u0_(m,n,k,j,i) = mu_*uold_ + nu_*u2_(m,n,k,j,i)
                         + (1.0-mu_-nu_)*u_ideal_(m,n,k,j,i) - dt*divf(i);
        if (lean_) u2_(m,n,k,j,i) = uold_;
      }
    });
  });
  return TaskStatus::complete;
}
