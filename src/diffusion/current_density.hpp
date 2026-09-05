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
//! \struct CurrentDensityGeom
//! \brief The mesh flags and coordinate arrays CurrentDensity() reads, gathered by value.
//!
//! CurrentDensity() runs on the device, so it must not be handed a MeshBlockPack*:
//! chasing pmy_pack->pmesh or pmy_pack->pcoord there dereferences a HOST pointer from
//! device code. That is legal only where the two share coherent memory (an MI300A APU)
//! and faults or returns garbage on any discrete GPU. Everything the kernel needs is
//! therefore collected on the host by MakeCurrentDensityGeom() and captured by value:
//! Kokkos Views are device handles and copy correctly into a lambda, plain bools copy as
//! themselves.
//!
//! The six coordinate arrays are only read on the spherical-polar branch, but they are
//! always valid Views -- Coordinates leaves them at their 1x1x1x1 placeholder unless that
//! geometry is in use -- so copying them unconditionally is safe and costs a few hundred
//! bytes of kernel argument space.

struct CurrentDensityGeom {
  bool use_spherical_polar;
  bool multi_d;     // nx2 > 1: the x2-derivative terms of J1 and J3 exist (2D AND 3D)
  bool three_d;     // nx3 > 1: the x3-derivative terms of J1 and J2 exist
  DvceArray4D<Real> area1, area2, area3;  // Coordinates::areaedge.x{1,2,3}e
  DvceArray4D<Real> dx1, dx2, dx3;        // Coordinates::dxface.x{1,2,3}f
};

//----------------------------------------------------------------------------------------
//! \fn MakeCurrentDensityGeom()
//  \brief Host-side gather of everything CurrentDensity() needs. Call this OUTSIDE any
//  par_for, then capture the result by value.

inline CurrentDensityGeom MakeCurrentDensityGeom(MeshBlockPack *pmy_pack) {
  CurrentDensityGeom geom;
  geom.use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  // Mesh::two_d is EXCLUSIVE (nx2 > 1 and nx3 == 1), so a 3D run has two_d == false.
  // Gating the x2-derivative terms on it dropped d(B3)/dx2 from J1 and d(B1)/dx2 from J3
  // in every 3D run between 1123736f (2026-06-20) and this fix: curl B was missing two
  // of its six terms on the Cartesian and spherical-polar grids (the cubed sphere has
  // its own curl in resistivity_gnomonic.cpp).  Found by sp_test iprob = 11, whose
  // force-free field decayed at half the exact rate, slower still with resolution.
  geom.multi_d = pmy_pack->pmesh->multi_d;
  geom.three_d = pmy_pack->pmesh->three_d;
  geom.area1 = pmy_pack->pcoord->areaedge.x1e;
  geom.area2 = pmy_pack->pcoord->areaedge.x2e;
  geom.area3 = pmy_pack->pcoord->areaedge.x3e;
  geom.dx1 = pmy_pack->pcoord->dxface.x1f;
  geom.dx2 = pmy_pack->pcoord->dxface.x2f;
  geom.dx3 = pmy_pack->pcoord->dxface.x3f;
  return geom;
}

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
void CurrentDensity(const CurrentDensityGeom &geom, TeamMember_t const &member,
     const int m, const int k, const int j,
     const int il, const int iu, const DvceFaceFld4D<Real> &b, const RegionSize &size,
     ScrArray1D<Real> &j1, ScrArray1D<Real> &j2, ScrArray1D<Real> &j3) {
  const bool use_spherical_polar = geom.use_spherical_polar;
  const auto &area1 = geom.area1;
  const auto &area2 = geom.area2;
  const auto &area3 = geom.area3;
  const auto &dx1 = geom.dx1;
  const auto &dx2 = geom.dx2;
  const auto &dx3 = geom.dx3;
  if (use_spherical_polar) {
      par_for_inner(member, il, iu, [&](const int i) {
        j1(i) = 0.0;
        j2(i) = -(dx3(m,k,j,i)*b.x3f(m,k,j,i) - dx3(m,k,j,i-1)*b.x3f(m,k,j,i-1))/area2(m,k,j,i);
        j3(i) =  (dx2(m,k,j,i)*b.x2f(m,k,j,i) - dx2(m,k,j,i-1)*b.x2f(m,k,j,i-1))/area3(m,k,j,i);
      });
      member.team_barrier();

      if (geom.multi_d) {
        par_for_inner(member, il, iu, [&](const int i) {
          j1(i) += (dx3(m,k,j,i)*b.x3f(m,k,j,i) - dx3(m,k,j-1,i)*b.x3f(m,k,j-1,i))/area1(m,k,j,i);
          j3(i) -= (dx1(m,k,j,i)*b.x1f(m,k,j,i) - dx1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/area3(m,k,j,i);
        });
        member.team_barrier();
      }

      if (geom.three_d) {
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

  if (geom.multi_d) {
    par_for_inner(member, il, iu, [&](const int i) {
      j1(i) += (b.x3f(m,k,j,i) - b.x3f(m,k,j-1,i))/size.dx2;
      j3(i) -= (b.x1f(m,k,j,i) - b.x1f(m,k,j-1,i))/size.dx2;
    });
    member.team_barrier();
  }

  if (geom.three_d) {
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
