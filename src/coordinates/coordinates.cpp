//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coordinates.cpp
//! \brief
#include <iostream> // cout
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "cartesian_ks.hpp"
#include "coordinates.hpp"
#include "cell_locations.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
// constructor, initializes coordinates data

Coordinates::Coordinates(ParameterInput *pin, MeshBlockPack *ppack) :
    pmy_pack(ppack),
    excision_floor("excision_floor",1,1,1,1),
    excision_flux("excision_flux",1,1,1,1),
    volume("volume",1,1,1,1),
    area("area",1,1,1,1),
    dx1("dx1",1,1,1,1), dx2("dx2",1,1,1,1), dx3("dx3",1,1,1,1),
    sin_cell("sin_cell",1,1,1), cos_cell("cos_cell",1,1,1),
    sin_face1("sin_face1",1,1,1), cos_face1("cos_face1",1,1,1),
    sin_face2("sin_face2",1,1,1), cos_face2("cos_face2",1,1,1),
    x_ov_rD("x_ov_rD",1,1,1,1), y_ov_rC("y_ov_rC",1,1,1,1), z_ov_rE("z_ov_rC",1,1,1,1) {
        
  if (pmy_pack->pmesh->use_cubed_sphere || pmy_pack->pmesh->use_spherical_polar) {
    // Total number of MeshBlocks on this rank to be used in array dimensioning
    int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(volume, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(area.x1f, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(area.x2f, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(area.x3f, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(dx1, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(dx2, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(dx3, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(x_ov_rD, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(y_ov_rC, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(z_ov_rE, nmb, ncells3, ncells2, ncells1);
  if (pmy_pack->pmesh->use_cubed_sphere) {
    Kokkos::realloc(sin_cell, nmb, ncells2, ncells1);
    Kokkos::realloc(cos_cell, nmb, ncells2, ncells1);
    Kokkos::realloc(sin_face1, nmb, ncells2, ncells1);
    Kokkos::realloc(cos_face1, nmb, ncells2, ncells1);
    Kokkos::realloc(sin_face2, nmb, ncells2, ncells1);
    Kokkos::realloc(cos_face2, nmb, ncells2, ncells1);
    CoordGnomonicEquiangle();
  }
  if (pmy_pack->pmesh->use_spherical_polar) {
    CoordSphericalPolar();
  }
  }
        
  // Check for relativistic dynamics
  // WGC: idea for handling new EOS
  is_dynamical_relativistic = (pin->DoesBlockExist("adm") || pin->DoesBlockExist("z4c"))
                         && (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd"));
  if(!is_dynamical_relativistic) {
    is_special_relativistic = pin->GetOrAddBoolean("coord","special_rel",false);
    is_general_relativistic = pin->GetOrAddBoolean("coord","general_rel",false);
  } else {
    is_special_relativistic = is_general_relativistic = false;
  }
  if (is_special_relativistic && is_general_relativistic) {
    std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
              << "Cannot specify both SR and GR at same time" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Read properties of metric and excision from input file for GR.
  if (is_general_relativistic || is_dynamical_relativistic) {
    coord_data.is_minkowski = pin->GetOrAddBoolean("coord","minkowski",false);
    if (!(coord_data.is_minkowski)) {
      coord_data.bh_spin = pin->GetReal("coord","a");
      coord_data.bh_excise = pin->GetOrAddBoolean("coord","excise",true);
    } else {
      coord_data.bh_spin = 0.0;
      coord_data.bh_excise = false;
    }

    if (coord_data.bh_excise) {
      // Set the density and pressure to which cells inside the excision radius will
      // be reset to.  Primitive velocities will be set to zero.
      coord_data.dexcise = pin->GetReal("coord","dexcise");
      coord_data.pexcise = pin->GetReal("coord","pexcise");
      coord_data.flux_excise_r = (pin->DoesBlockExist("radiation")) ?
        1.0+sqrt(1.0-SQR(coord_data.bh_spin)) :
        pin->GetOrAddReal("coord","flux_excise_r",1.0);
      coord_data.rexcise =
        (pin->DoesBlockExist("radiation")) ? 1.0+sqrt(1.0-SQR(coord_data.bh_spin)) : 1.0;

      coord_data.excision_scheme = ExcisionScheme::fixed;
      if (is_dynamical_relativistic) {
        std::string emethod = pin->GetOrAddString("coord","excision_scheme","fixed");
        if (emethod.compare("fixed") == 0) {
          coord_data.excision_scheme = ExcisionScheme::fixed;
        } else if (emethod.compare("lapse") == 0) {
          coord_data.excision_scheme = ExcisionScheme::lapse;
          coord_data.excise_lapse = pin->GetOrAddReal("coord","excise_lapse", 0.25);
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
                    << __LINE__ << std::endl
                    << "Unknown excision method: " << emethod << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      // boolean masks allocation
      int nmb = ppack->nmb_thispack;
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(excision_floor, nmb, ncells3, ncells2, ncells1);
      Kokkos::realloc(excision_flux, nmb, ncells3, ncells2, ncells1);
      if (coord_data.excision_scheme == ExcisionScheme::fixed) {
        SetExcisionMasks(excision_floor, excision_flux);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn
// Coordinate (geometric) source term function for GR hydrodynamics

void Coordinates::CoordSrcTerms(const DvceArray5D<Real> &prim, const EOS_Data &eos,
                                const Real dt, DvceArray5D<Real> &cons) {
  // capture variables for kernel
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is; int ie = indcs.ie;
  int js = indcs.js; int je = indcs.je;
  int ks = indcs.ks; int ke = indcs.ke;
  auto &size = pmy_pack->pmb->mb_size;
  auto &flat = coord_data.is_minkowski;
  auto &spin = coord_data.bh_spin;

  Real gamma_prime = eos.gamma / (eos.gamma - 1.0);

  int nmb1 = pmy_pack->nmb_thispack - 1;
  par_for("coord_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Extract components of metric
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    Real glower[4][4], gupper[4][4];
    ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

    // Extract primitives
    const Real &rho  = prim(m,IDN,k,j,i);
    const Real &uu1  = prim(m,IVX,k,j,i);
    const Real &uu2  = prim(m,IVY,k,j,i);
    const Real &uu3  = prim(m,IVZ,k,j,i);
    Real pgas = eos.IdealGasPressure(prim(m,IEN,k,j,i));

    // Calculate 4-velocity (exploiting symmetry of metric)
    Real uu_sq = glower[1][1]*uu1*uu1 +2.0*glower[1][2]*uu1*uu2 +2.0*glower[1][3]*uu1*uu3
               + glower[2][2]*uu2*uu2 +2.0*glower[2][3]*uu2*uu3
               + glower[3][3]*uu3*uu3;
    Real alpha = sqrt(-1.0/gupper[0][0]);
    Real gamma = sqrt(1.0 + uu_sq);
    Real u0 = gamma / alpha;
    Real u1 = uu1 - alpha * gamma * gupper[0][1];
    Real u2 = uu2 - alpha * gamma * gupper[0][2];
    Real u3 = uu3 - alpha * gamma * gupper[0][3];

    // Calculate stress-energy tensor
    Real wtot = rho + gamma_prime * pgas;
    Real ptot = pgas;
    Real tt[4][4];
    tt[0][0] = wtot * u0 * u0 + ptot * gupper[0][0];
    tt[0][1] = wtot * u0 * u1 + ptot * gupper[0][1];
    tt[0][2] = wtot * u0 * u2 + ptot * gupper[0][2];
    tt[0][3] = wtot * u0 * u3 + ptot * gupper[0][3];
    tt[1][1] = wtot * u1 * u1 + ptot * gupper[1][1];
    tt[1][2] = wtot * u1 * u2 + ptot * gupper[1][2];
    tt[1][3] = wtot * u1 * u3 + ptot * gupper[1][3];
    tt[2][2] = wtot * u2 * u2 + ptot * gupper[2][2];
    tt[2][3] = wtot * u2 * u3 + ptot * gupper[2][3];
    tt[3][3] = wtot * u3 * u3 + ptot * gupper[3][3];

    // compute derivates of metric.
    Real dg_dx1[4][4], dg_dx2[4][4], dg_dx3[4][4];
    ComputeMetricDerivatives(x1v, x2v, x3v, flat, spin, dg_dx1, dg_dx2, dg_dx3);

    // Calculate source terms, exploiting symmetries
    Real s_1 = 0.0, s_2 = 0.0, s_3 = 0.0;
    s_1 += 0.5*dg_dx1[0][0] * tt[0][0];
    s_1 +=     dg_dx1[0][1] * tt[0][1];
    s_1 +=     dg_dx1[0][2] * tt[0][2];
    s_1 +=     dg_dx1[0][3] * tt[0][3];
    s_1 += 0.5*dg_dx1[1][1] * tt[1][1];
    s_1 +=     dg_dx1[1][2] * tt[1][2];
    s_1 +=     dg_dx1[1][3] * tt[1][3];
    s_1 += 0.5*dg_dx1[2][2] * tt[2][2];
    s_1 +=     dg_dx1[2][3] * tt[2][3];
    s_1 += 0.5*dg_dx1[3][3] * tt[3][3];

    s_2 += 0.5*dg_dx2[0][0] * tt[0][0];
    s_2 +=     dg_dx2[0][1] * tt[0][1];
    s_2 +=     dg_dx2[0][2] * tt[0][2];
    s_2 +=     dg_dx2[0][3] * tt[0][3];
    s_2 += 0.5*dg_dx2[1][1] * tt[1][1];
    s_2 +=     dg_dx2[1][2] * tt[1][2];
    s_2 +=     dg_dx2[1][3] * tt[1][3];
    s_2 += 0.5*dg_dx2[2][2] * tt[2][2];
    s_2 +=     dg_dx2[2][3] * tt[2][3];
    s_2 += 0.5*dg_dx2[3][3] * tt[3][3];

    s_3 += 0.5*dg_dx3[0][0] * tt[0][0];
    s_3 +=     dg_dx3[0][1] * tt[0][1];
    s_3 +=     dg_dx3[0][2] * tt[0][2];
    s_3 +=     dg_dx3[0][3] * tt[0][3];
    s_3 += 0.5*dg_dx3[1][1] * tt[1][1];
    s_3 +=     dg_dx3[1][2] * tt[1][2];
    s_3 +=     dg_dx3[1][3] * tt[1][3];
    s_3 += 0.5*dg_dx3[2][2] * tt[2][2];
    s_3 +=     dg_dx3[2][3] * tt[2][3];
    s_3 += 0.5*dg_dx3[3][3] * tt[3][3];

    // Add source terms to conserved quantities
    cons(m,IM1,k,j,i) += dt * s_1;
    cons(m,IM2,k,j,i) += dt * s_2;
    cons(m,IM3,k,j,i) += dt * s_3;
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn
// Coordinate (geometric) source term function for GR MHD
//
// TODO(@user): Most of this function just copies the Hydro version.  Only difference is
// the inclusion of the magnetic field in computing the stress-energy tensor.  There must
// be a smarter way to generalize these two functions and avoid duplicated code.
// Functions distinguished only by argument list.

void Coordinates::CoordSrcTerms(const DvceArray5D<Real> &prim,
                                const DvceArray5D<Real> &bcc, const EOS_Data &eos,
                                const Real dt, DvceArray5D<Real> &cons) {
  // capture variables for kernel
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is; int ie = indcs.ie;
  int js = indcs.js; int je = indcs.je;
  int ks = indcs.ks; int ke = indcs.ke;
  auto &size = pmy_pack->pmb->mb_size;
  auto &flat = coord_data.is_minkowski;
  auto &spin = coord_data.bh_spin;

  Real gamma_prime = eos.gamma / (eos.gamma - 1.0);

  int nmb1 = pmy_pack->nmb_thispack - 1;
  par_for("coord_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Extract components of metric
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    Real glower[4][4], gupper[4][4];
    ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

    // Extract primitives
    const Real &rho  = prim(m,IDN,k,j,i);
    const Real &uu1  = prim(m,IVX,k,j,i);
    const Real &uu2  = prim(m,IVY,k,j,i);
    const Real &uu3  = prim(m,IVZ,k,j,i);
    Real pgas = eos.IdealGasPressure(prim(m,IEN,k,j,i));

    // Calculate 4-velocity
    Real uu_sq = glower[1][1]*uu1*uu1 +2.0*glower[1][2]*uu1*uu2 +2.0*glower[1][3]*uu1*uu3
               + glower[2][2]*uu2*uu2 +2.0*glower[2][3]*uu2*uu3
               + glower[3][3]*uu3*uu3;
    Real alpha = sqrt(-1.0/gupper[0][0]);
    Real gamma = sqrt(1.0 + uu_sq);
    Real u0 = gamma / alpha;
    Real u1 = uu1 - alpha * gamma * gupper[0][1];
    Real u2 = uu2 - alpha * gamma * gupper[0][2];
    Real u3 = uu3 - alpha * gamma * gupper[0][3];

    // lower vector indices
    Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
    Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
    Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;

    // calculate 4-magnetic field
    const Real &bb1 = bcc(m,IBX,k,j,i);
    const Real &bb2 = bcc(m,IBY,k,j,i);
    const Real &bb3 = bcc(m,IBZ,k,j,i);
    Real b0 = u_1*bb1 + u_2*bb2 + u_3*bb3;
    Real b1 = (bb1 + b0 * u1) / u0;
    Real b2 = (bb2 + b0 * u2) / u0;
    Real b3 = (bb3 + b0 * u3) / u0;

    // lower vector indices
    Real b_0 = glower[0][0]*b0 + glower[0][1]*b1 + glower[0][2]*b2 + glower[0][3]*b3;
    Real b_1 = glower[1][0]*b0 + glower[1][1]*b1 + glower[1][2]*b2 + glower[1][3]*b3;
    Real b_2 = glower[2][0]*b0 + glower[2][1]*b1 + glower[2][2]*b2 + glower[2][3]*b3;
    Real b_3 = glower[3][0]*b0 + glower[3][1]*b1 + glower[3][2]*b2 + glower[3][3]*b3;
    Real b_sq = b_0*b0 + b_1*b1 + b_2*b2 + b_3*b3;

    // Calculate stress-energy tensor
    Real wtot = rho + gamma_prime * pgas + b_sq;
    Real ptot = pgas + 0.5*b_sq;
    Real tt[4][4];
    tt[0][0] = wtot * u0 * u0 + ptot * gupper[0][0] - b0 * b0;
    tt[0][1] = wtot * u0 * u1 + ptot * gupper[0][1] - b0 * b1;
    tt[0][2] = wtot * u0 * u2 + ptot * gupper[0][2] - b0 * b2;
    tt[0][3] = wtot * u0 * u3 + ptot * gupper[0][3] - b0 * b3;
    tt[1][1] = wtot * u1 * u1 + ptot * gupper[1][1] - b1 * b1;
    tt[1][2] = wtot * u1 * u2 + ptot * gupper[1][2] - b1 * b2;
    tt[1][3] = wtot * u1 * u3 + ptot * gupper[1][3] - b1 * b3;
    tt[2][2] = wtot * u2 * u2 + ptot * gupper[2][2] - b2 * b2;
    tt[2][3] = wtot * u2 * u3 + ptot * gupper[2][3] - b2 * b3;
    tt[3][3] = wtot * u3 * u3 + ptot * gupper[3][3] - b3 * b3;

    // compute derivates of metric.
    Real dg_dx1[4][4], dg_dx2[4][4], dg_dx3[4][4];
    ComputeMetricDerivatives(x1v, x2v, x3v, flat, spin, dg_dx1, dg_dx2, dg_dx3);

    // Calculate source terms
    Real s_1 = 0.0, s_2 = 0.0, s_3 = 0.0;
    s_1 += 0.5*dg_dx1[0][0] * tt[0][0];
    s_1 +=     dg_dx1[0][1] * tt[0][1];
    s_1 +=     dg_dx1[0][2] * tt[0][2];
    s_1 +=     dg_dx1[0][3] * tt[0][3];
    s_1 += 0.5*dg_dx1[1][1] * tt[1][1];
    s_1 +=     dg_dx1[1][2] * tt[1][2];
    s_1 +=     dg_dx1[1][3] * tt[1][3];
    s_1 += 0.5*dg_dx1[2][2] * tt[2][2];
    s_1 +=     dg_dx1[2][3] * tt[2][3];
    s_1 += 0.5*dg_dx1[3][3] * tt[3][3];

    s_2 += 0.5*dg_dx2[0][0] * tt[0][0];
    s_2 +=     dg_dx2[0][1] * tt[0][1];
    s_2 +=     dg_dx2[0][2] * tt[0][2];
    s_2 +=     dg_dx2[0][3] * tt[0][3];
    s_2 += 0.5*dg_dx2[1][1] * tt[1][1];
    s_2 +=     dg_dx2[1][2] * tt[1][2];
    s_2 +=     dg_dx2[1][3] * tt[1][3];
    s_2 += 0.5*dg_dx2[2][2] * tt[2][2];
    s_2 +=     dg_dx2[2][3] * tt[2][3];
    s_2 += 0.5*dg_dx2[3][3] * tt[3][3];

    s_3 += 0.5*dg_dx3[0][0] * tt[0][0];
    s_3 +=     dg_dx3[0][1] * tt[0][1];
    s_3 +=     dg_dx3[0][2] * tt[0][2];
    s_3 +=     dg_dx3[0][3] * tt[0][3];
    s_3 += 0.5*dg_dx3[1][1] * tt[1][1];
    s_3 +=     dg_dx3[1][2] * tt[1][2];
    s_3 +=     dg_dx3[1][3] * tt[1][3];
    s_3 += 0.5*dg_dx3[2][2] * tt[2][2];
    s_3 +=     dg_dx3[2][3] * tt[2][3];
    s_3 += 0.5*dg_dx3[3][3] * tt[3][3];

    // Add source terms to conserved quantities
    cons(m,IM1,k,j,i) += dt * s_1;
    cons(m,IM2,k,j,i) += dt * s_2;
    cons(m,IM3,k,j,i) += dt * s_3;
  });

  return;
}

void Coordinates::CoordGnomonicEquiangle() {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  par_for("cscoord", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    Real r_c, r_l, r_r, dxr;
    if (pmy_pack->pmesh->three_d) {
        // --- radial ---
        r_c  = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
        r_l  = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
        r_r  = LeftEdgeX(k+1-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
        dxr = size.d_view(m).dx3;
    } else {
        r_c = 1.0;
        r_l = 1.0;
        r_r = 1.0;
        dxr = 1.0;
    }

    // --- angles ---
    Real xi   = M_PI/4.0 * CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real eta  = M_PI/4.0 * CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real xil  = M_PI/4.0 * LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real xir  = M_PI/4.0 * LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real etal = M_PI/4.0 * LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real etar = M_PI/4.0 * LeftEdgeX(j+1-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);

    Real x  = tan(xi);
    Real y  = tan(eta);
    Real xl = tan(xil);
    Real xr = tan(xir);
    Real yl = tan(etal);
    Real yr = tan(etar);

    Real C = sqrt(1.0 + SQR(x));
    Real Cl = sqrt(1.0 + SQR(xl));
    Real Cr = sqrt(1.0 + SQR(xr));
    Real D = sqrt(1.0 + SQR(y));
    Real Dl = sqrt(1.0 + SQR(yl));
    Real Dr = sqrt(1.0 + SQR(yr));
      
    sin_cell(m,j,i) = sqrt(1.0 + SQR(x) + SQR(y)) / (C * D);
    cos_cell(m,j,i) = - x * y / (C * D);
    sin_face1(m,j,i) = sqrt(1.0 + SQR(xl) + SQR(y)) / (Cl * D);
    cos_face1(m,j,i) = - xl * y / (Cl * D);
    sin_face2(m,j,i) = sqrt(1.0 + SQR(x) + SQR(yl)) / (C * Dl);
    cos_face2(m,j,i) = - x * yl / (C * Dl);
    Real sin_face1r = sqrt(1.0 + SQR(xr) + SQR(y)) / (Cr * D);
    Real sin_face2r = sqrt(1.0 + SQR(x) + SQR(yr)) / (C * Dr);
    if (i == ie) {
      sin_face1(m,j,i+1) = sin_face1r;
      cos_face1(m,j,i+1) = - xr * y / (Cr * D);
    }
    if (j == je) {
      sin_face2(m,j+1,i) = sin_face2r;
      cos_face2(m,j+1,i) = - x * yr / (C * Dr);
    }

    // --- angular face widths ---
    Real dth_xi = acos((1.0 + xl*xr + y*y) / sqrt(1.0 + xl*xl + y*y) / sqrt(1.0 + xr*xr + y*y));
    Real dth_xil = acos((1.0 + xl*xr + yl*yl) / sqrt(1.0 + xl*xl + yl*yl) / sqrt(1.0 + xr*xr + yl*yl));
    Real dth_xir = acos((1.0 + xl*xr + yr*yr) / sqrt(1.0 + xl*xl + yr*yr) / sqrt(1.0 + xr*xr + yr*yr));
    Real dth_eta = acos((1.0 + x*x + yl*yr) / sqrt(1.0 + x*x + yl*yl) / sqrt(1.0 + x*x + yr*yr));
    Real dth_etal = acos((1.0 + xl*xl + yl*yr) / sqrt(1.0 + xl*xl + yl*yl) / sqrt(1.0 + xl*xl + yr*yr));
    Real dth_etar = acos((1.0 + xr*xr + yl*yr) / sqrt(1.0 + xr*xr + yl*yl) / sqrt(1.0 + xr*xr + yr*yr));

    // --- radial faces ---
    area.x3f(m,k,j,i) = SQR(r_l) * dth_xi * dth_eta * sin_cell(m,j,i);
    Real area3r = SQR(r_r) * dth_xi * dth_eta * sin_cell(m,j,i);
    if (k == ke) area.x3f(m,k+1,j,i) = area3r;

    // --- xi faces ---
    area.x1f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etal;
    Real area1r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etar;
    if (i == ie) area.x1f(m,k,j,i+1) = area1r;

    // --- eta faces ---
    area.x2f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xil;
    Real area2r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xir;
    if (j == je) area.x2f(m,k,j+1,i) = area2r;

    // --- volume (SNAP trapezoidal) ---

    volume(m,k,j,i) = 1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l) * dth_xi * dth_eta * sin_cell(m,j,i);
    dx1(m,k,j,i) = r_c * dth_xi;
    dx2(m,k,j,i) = r_c * dth_eta;
      
    x_ov_rD(m,k,j,i) = (area1r * sin_face1r - area.x1f(m,k,j,i) * sin_face1(m,j,i)) / volume(m,k,j,i);
    y_ov_rC(m,k,j,i) = (area2r * sin_face2r - area.x2f(m,k,j,i) * sin_face2(m,j,i)) / volume(m,k,j,i);
    z_ov_rE(m,k,j,i) = (area3r - area.x3f(m,k,j,i)) / volume(m,k,j,i);
  });
}

void Coordinates::SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0, const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
    
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  par_for("cssrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      
    Real radius = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real dr = size.d_view(m).dx3;
      
    Real v1 = w0(m,IVX,k,j,i);
    Real v2 = w0(m,IVY,k,j,i);
    Real v3 = w0(m,IVZ,k,j,i);
    Real pr = w0(m,IEN,k,j,i)*gm1;
    Real rho = w0(m,IDN,k,j,i);
      
    Real cosine = cos_cell(m,j,i);
    Real sine2 = SQR(sin_cell(m,j,i));
    Real v_1 = v1 + v2 * cosine;
    Real v_2 = v2 + v1 * cosine;
      
    Real src3 = z_ov_rE(m,k,j,i) * (pr + 0.5*rho*(v1*v_1+v2*v_2));
    Real src1 = x_ov_rD(m,k,j,i) * (pr + rho*SQR(v2)*sine2);// - rho*v3*v_1 / radius;
    Real src2 = y_ov_rC(m,k,j,i) * (pr + rho*SQR(v1)*sine2);// - rho*v3*v_2 / radius;
      
    src1 -= dr/2.0/radius * (uflx.x3f(m,IM1,k,j,i)*area.x3f(m,k,j,i)+uflx.x3f(m,IM1,k+1,j,i)*area.x3f(m,k+1,j,i))/volume(m,k,j,i);
    src2 -= dr/2.0/radius * (uflx.x3f(m,IM2,k,j,i)*area.x3f(m,k,j,i)+uflx.x3f(m,IM2,k+1,j,i)*area.x3f(m,k+1,j,i))/volume(m,k,j,i);
      
    u0(m,IVX,k,j,i) += src1*bdt;
    u0(m,IVY,k,j,i) += src2*bdt;
    u0(m,IVZ,k,j,i) += src3*bdt;
    
  });
}

void Coordinates::CoordSphericalPolar() {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  par_for("spcoord", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    Real r_c, r_l, r_r;
    if (pmy_pack->pmesh->three_d) {
        // --- radial ---
        r_c  = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
        r_l  = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
        r_r  = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
        // stretch the radial grid
        if (pmy_pack->pmesh->use_grid_stretch) {
          Real r0 = pmy_pack->pmesh->mesh_size.x1min;
          Real r1 = pmy_pack->pmesh->mesh_size.x1max;
          StretchR(r0,r1,r_c);
          StretchR(r0,r1,r_l);
          StretchR(r0,r1,r_r);
        }
    } else {
        r_c = 1.0;
        r_l = 0.0;
        r_r = 1.0;
    }

    // --- angles ---
    Real theta  = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real phi    = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real thetal = LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real thetar = LeftEdgeX(j+1-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real phil   = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phir   = LeftEdgeX(k+1-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);

    Real sinl   = fabs(sin(thetal));
    Real sinr   = fabs(sin(thetar));
    Real sinc   = fabs(sin(theta));
    Real cosl   = cos(thetal);
    Real cosr   = cos(thetar);

    // --- radial faces ---
    area.x1f(m,k,j,i) = SQR(r_l) * fabs(cosl-cosr) * (phir-phil);
    Real area1r = SQR(r_r) * fabs(cosl-cosr) * (phir-phil);
    if (i == ie) area.x1f(m,k,j,i+1) = area1r;

    // --- theta faces ---
    area.x2f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * sinl * (phir-phil);
    Real area2r = 0.5 * (SQR(r_r)-SQR(r_l)) * sinr * (phir-phil);
    if (j == je) area.x2f(m,k,j+1,i) = area2r;

    // --- phi faces ---
    area.x3f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * (thetar-thetal);
    Real area3r = 0.5 * (SQR(r_r)-SQR(r_l)) * (thetar-thetal);
    if (k == ke) area.x3f(m,k+1,j,i) = area3r;

    // --- volume ---

    volume(m,k,j,i) = 1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l) * fabs(cosl-cosr) * (phir-phil);
    dx2(m,k,j,i) = r_c * (thetar-thetal);
    dx3(m,k,j,i) = r_c * sinc * (phir-phil);
    dx1(m,k,j,i) = r_r-r_l;
      
    z_ov_rE(m,k,j,i) = (area1r - area.x1f(m,k,j,i)) / volume(m,k,j,i);
    x_ov_rD(m,k,j,i) = (area2r - area.x2f(m,k,j,i)) / volume(m,k,j,i);
    y_ov_rC(m,k,j,i) = (sinr-sinl)/(sinr+sinl);
  });
}

void Coordinates::SrcTermsSphericalPolar(const DvceArray5D<Real> &w0, const DvceArray5D<Real> &w0wb, const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
    
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  par_for("spsrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      
    Real r_l  = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real r_r  = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    // stretch the radial grid
    if (pmy_pack->pmesh->use_grid_stretch) {
      Real r0 = pmy_pack->pmesh->mesh_size.x1min;
      Real r1 = pmy_pack->pmesh->mesh_size.x1max;
      StretchR(r0,r1,r_l);
      StretchR(r0,r1,r_r);
    }
    Real factor = (r_r-r_l)/(r_r+r_l);
      
    Real v1 = w0(m,IVX,k,j,i);
    Real v2 = w0(m,IVY,k,j,i);
    Real v3 = w0(m,IVZ,k,j,i);
    Real pr = w0(m,IEN,k,j,i)*gm1;
    if (pmy_pack->phydro->use_wellbalance) pr -= w0wb(m,IEN,k,j,i)*gm1;
    Real rho = w0(m,IDN,k,j,i);
      
    Real src1 = z_ov_rE(m,k,j,i) * (pr + 0.5*rho*(v2*v2+v3*v3));
    Real src2 = x_ov_rD(m,k,j,i) * (pr + rho*SQR(v3));
    Real src3 = -y_ov_rC(m,k,j,i) * (uflx.x2f(m,IM3,k,j,i)*area.x2f(m,k,j,i)+uflx.x2f(m,IM3,k,j+1,i)*area.x2f(m,k,j+1,i))/volume(m,k,j,i);
      
    src2 -= factor * (uflx.x1f(m,IM2,k,j,i)*area.x1f(m,k,j,i)+uflx.x1f(m,IM2,k,j,i+1)*area.x1f(m,k,j,i+1))/volume(m,k,j,i);
    src3 -= factor * (uflx.x1f(m,IM3,k,j,i)*area.x1f(m,k,j,i)+uflx.x1f(m,IM3,k,j,i+1)*area.x1f(m,k,j,i+1))/volume(m,k,j,i);
      
    u0(m,IM1,k,j,i) += src1*bdt;
    u0(m,IM2,k,j,i) += src2*bdt;
    u0(m,IM3,k,j,i) += src3*bdt;
    
  });
}
