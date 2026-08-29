//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_gnomonic.cpp
//! \brief The resistive EMF on the CUBED SPHERE.
//!
//! WHY THIS NEEDS ITS OWN FILE. `current_density.hpp` branches on `use_spherical_polar`
//! only, so on the cubed sphere it fell through to the Cartesian formula and divided by
//! `size.dx1/dx2/dx3` -- which on this grid are index spacings, not lengths (x2 and x3
//! are equiangular coordinates on [-1,1]). That is wrong by an O(1) factor everywhere.
//! Its curvilinear branch is not usable either, for a reason that is not a matter of
//! filling in different metric coefficients:
//!
//!   THE TANGENT BASIS IS NOT ORTHOGONAL. `b0.x2f` and `b0.x3f` hold B.nhat, the flux
//!   density through the face, and on the cubed sphere nhat_xi is perpendicular to
//!   e_eta, NOT parallel to e_xi. Stokes' theorem needs the component ALONG each edge,
//!   B.e_xi, and the two differ by the metric:
//!       B.e_xi = (B.nhat_xi + c B.nhat_eta)/s,   c = cos_cell, s = sin_cell,
//!   so the correction is O(c) -- O(1) near a panel corner, not O(h). This is the same
//!   trap that made the cell-centred EMF inconsistent in mhd_corner_e.cpp.
//!
//! The curl is therefore taken in two passes. Pass one applies Stokes' theorem on the
//! loop through the four cell CENTRES around each edge (lengths from Coordinates::dxface,
//! enclosed area from Coordinates::areaedge), which yields J in the FACE-NORMAL frame
//! {rhat, nhat_xi, nhat_eta}. Pass two rotates it into the frame the edge EMFs are stored
//! in -- E.that along the edge, which is what mhd_ct.cpp's dxedge*E/area consumes and
//! what GnomonicEquiangleEmfX1 produces for the ideal EMF -- and multiplies by eta. The
//! two passes cannot be fused: the rotation mixes J at x2-edges with J at x3-edges, which
//! are different points.

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "resistivity.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
//! \brief B.e_xi at the centre of an x2 face, from the face-normal components. The
//! partner component lives on the four x3 faces that share this face's corners.

KOKKOS_INLINE_FUNCTION
Real BcovXi(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
            const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
            const int i) {
  const Real b3a = 0.25*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i)
                       + b.x3f(m,k,j-1,i) + b.x3f(m,k+1,j-1,i));
  return (b.x2f(m,k,j,i) + cosf(m,k,j)*b3a)/sinf(m,k,j);
}

//----------------------------------------------------------------------------------------
//! \brief B.e_eta at the centre of an x3 face.

KOKKOS_INLINE_FUNCTION
Real BcovEta(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
             const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
             const int i) {
  const Real b2a = 0.25*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i)
                       + b.x2f(m,k-1,j,i) + b.x2f(m,k-1,j+1,i));
  return (b.x3f(m,k,j,i) + cosf(m,k,j)*b2a)/sinf(m,k,j);
}

//----------------------------------------------------------------------------------------
//! \fn Resistivity::AddEMFGnomonicResist
//! \brief Adds E_resistive = eta J to the edge-centred electric fields, cubed sphere.

void Resistivity::AddEMFGnomonicResist(const DvceFaceFld4D<Real> &b0,
                                       DvceEdgeFld4D<Real> &efld) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  if (!(pmy_pack->pmesh->three_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Resistivity on the cubed sphere requires a 3D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &cxi = pmy_pack->pcoord->cos_face_xi;
  auto &sxi = pmy_pack->pcoord->sin_face_xi;
  auto &cet = pmy_pack->pcoord->cos_face_eta;
  auto &set = pmy_pack->pcoord->sin_face_eta;
  auto &dxf1 = pmy_pack->pcoord->dxface.x1f;
  auto &dxf2 = pmy_pack->pcoord->dxface.x2f;
  auto &dxf3 = pmy_pack->pcoord->dxface.x3f;
  auto &ae1 = pmy_pack->pcoord->areaedge.x1e;
  auto &ae2 = pmy_pack->pcoord->areaedge.x2e;
  auto &ae3 = pmy_pack->pcoord->areaedge.x3e;
  auto b = b0;
  auto jn1 = jnorm.x1e;
  auto jn2 = jnorm.x2e;
  auto jn3 = jnorm.x3e;

  // ---- pass 1: J in the face-normal frame, by Stokes' theorem on the dual loop --------
  // Each loop is bounded by the four cell centres around the edge, so every side passes
  // through a face centre and carries that face's field. Pass 2 needs J one layer beyond
  // its own edge range in the two transverse directions, hence the -1 / +1 below.
  par_for("res_cs_j1", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    jn1(m,k,j,i) = (dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                  - dxf3(m,k,j-1,i)*BcovEta(b,cet,set,m,k,j-1,i)
                  - dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                  + dxf2(m,k-1,j,i)*BcovXi(b,cxi,sxi,m,k-1,j,i))/ae1(m,k,j,i);
  });

  par_for("res_cs_j2", DevExeSpace(), 0,nmb1, ks,ke+1, js-1,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    jn2(m,k,j,i) = (-dxf3(m,k,j,i)*BcovEta(b,cet,set,m,k,j,i)
                    + dxf3(m,k,j,i-1)*BcovEta(b,cet,set,m,k,j,i-1)
                    + dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                    - dxf1(m,k-1,j,i)*b.x1f(m,k-1,j,i))/ae2(m,k,j,i);
  });

  par_for("res_cs_j3", DevExeSpace(), 0,nmb1, ks-1,ke+1, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    jn3(m,k,j,i) = (dxf2(m,k,j,i)*BcovXi(b,cxi,sxi,m,k,j,i)
                  - dxf2(m,k,j,i-1)*BcovXi(b,cxi,sxi,m,k,j,i-1)
                  - dxf1(m,k,j,i)*b.x1f(m,k,j,i)
                  + dxf1(m,k,j-1,i)*b.x1f(m,k,j-1,i))/ae3(m,k,j,i);
  });

  // ---- pass 2: rotate into the edge (covariant) frame and multiply by eta -------------
  // J.e_xi = (J.nhat_xi + c J.nhat_eta)/s and likewise for eta, the same transform
  // GnomonicEquiangleRaiseVelMHD applies to the field. The radial edge needs no rotation:
  // rhat is orthogonal to both tangent directions.
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto eta_b_ = eta_b;

  par_for("res_cs_e1", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j-1,i));
    e1(m,k,j,i) += et*jn1(m,k,j,i);
  });

  par_for("res_cs_e2", DevExeSpace(), 0,nmb1, ks,ke+1, js,je, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j,i-1)
                        + eta_b_(m,k-1,j,i) + eta_b_(m,k-1,j,i-1));
    const Real j3a = 0.25*(jn3(m,k,j,i) + jn3(m,k,j+1,i)
                         + jn3(m,k-1,j,i) + jn3(m,k-1,j+1,i));
    e2(m,k,j,i) += et*(jn2(m,k,j,i) + cet(m,k,j)*j3a)/set(m,k,j);
  });

  par_for("res_cs_e3", DevExeSpace(), 0,nmb1, ks,ke, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real et = 0.25*(eta_b_(m,k,j,i) + eta_b_(m,k,j-1,i)
                        + eta_b_(m,k,j,i-1) + eta_b_(m,k,j-1,i-1));
    const Real j2a = 0.25*(jn2(m,k,j,i) + jn2(m,k+1,j,i)
                         + jn2(m,k,j-1,i) + jn2(m,k+1,j-1,i));
    e3(m,k,j,i) += et*(jn3(m,k,j,i) + cxi(m,k,j)*j2a)/sxi(m,k,j);
  });

  return;
}
