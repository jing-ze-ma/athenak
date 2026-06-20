#ifndef DIFFUSION_CURRENT_DENSITY_HPP_
#define DIFFUSION_CURRENT_DENSITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file current_density.hpp
//  \brief Inlined function to compute current density in 1-D pencils in i-direction

#include "athena.hpp"
#include "mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn CurrentDensity()
//  \brief Calculates the three components of the current density at cell edges
//  Each component of J is centered identically to the edge-electric-field
//               _____________
//               |\           \
//               | \           \
//               |  \___________\
//               |   |           |
//               \   |           |
//              J2*  *J3         |
//   x2 x3         \ |           |
//    \ |           \|_____*_____|
//     \|__x1             J1

KOKKOS_INLINE_FUNCTION
void CurrentDensity(MeshBlockPack *pmy_pack, TeamMember_t const &member, const int m, const int k, const int j,
     const int il, const int iu, const DvceFaceFld4D<Real> &b, const RegionSize &size,
     ScrArray1D<Real> &j1, ScrArray1D<Real> &j2, ScrArray1D<Real> &j3) {
  const bool use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &area1 = pmy_pack->pcoord->areaedge.x1e;
  auto &area2 = pmy_pack->pcoord->areaedge.x2e;
  auto &area3 = pmy_pack->pcoord->areaedge.x3e;
  auto &dx1 = pmy_pack->pcoord->dxface.x1f;
  auto &dx2 = pmy_pack->pcoord->dxface.x2f;
  auto &dx3 = pmy_pack->pcoord->dxface.x3f;
  if (use_spherical_polar) {
      par_for_inner(member, il, iu, [&](const int i) {
        j1(i) = 0.0;
        j2(i) = -(dx3(m,k,j,i)*b.x3f(m,k,j,i) - dx3(m,k,j,i-1)*b.x3f(m,k,j,i-1))/area2(m,k,j,i);
        j3(i) =  (dx2(m,k,j,i)*b.x2f(m,k,j,i) - dx2(m,k,j,i-1)*b.x2f(m,k,j,i-1))/area3(m,k,j,i);
      });
      member.team_barrier();

      if (pmy_pack->pmesh->two_d) {
        par_for_inner(member, il, iu, [&](const int i) {
          j1(i) += (dx3(m,k,j,i)*b.x3f(m,k,j,i) - dx3(m,k,j-1,i)*b.x3f(m,k,j-1,i))/area1(m,k,j,i);
          j3(i) -= (dx1(m,k,j,i)*b.x1f(m,k,j,i) - dx1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/area3(m,k,j,i);
        });
        member.team_barrier();
      }

      if (pmy_pack->pmesh->three_d) {
        par_for_inner(member, il, iu, [&](const int i) {
          j1(i) -= (dx2(m,k,j,i)*b.x2f(m,k,j,i) - dx2(m,k-1,j,i)*b.x2f(m,k-1,j,i))/area1(m,k,j,i);
          j2(i) += (dx1(m,k,j,i)*b.x1f(m,k,j,i) - dx1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/area2(m,k,j,i);
        });
        member.team_barrier();
      }
  } else {
  par_for_inner(member, il, iu, [&](const int i) {
    j1(i) = 0.0;
    j2(i) = -(b.x3f(m,k,j,i) - b.x3f(m,k,j,i-1))/size.dx1;
    j3(i) =  (b.x2f(m,k,j,i) - b.x2f(m,k,j,i-1))/size.dx1;
  });
  member.team_barrier();

  if (pmy_pack->pmesh->two_d) {
    par_for_inner(member, il, iu, [&](const int i) {
      j1(i) += (b.x3f(m,k,j,i) - b.x3f(m,k,j-1,i))/size.dx2;
      j3(i) -= (b.x1f(m,k,j,i) - b.x1f(m,k,j-1,i))/size.dx2;
    });
    member.team_barrier();
  }

  if (pmy_pack->pmesh->three_d) {
    par_for_inner(member, il, iu, [&](const int i) {
      j1(i) -= (b.x2f(m,k,j,i) - b.x2f(m,k-1,j,i))/size.dx3;
      j2(i) += (b.x1f(m,k,j,i) - b.x1f(m,k-1,j,i))/size.dx3;
    });
    member.team_barrier();
  }
  }
  return;
}

#endif // DIFFUSION_CURRENT_DENSITY_HPP_
