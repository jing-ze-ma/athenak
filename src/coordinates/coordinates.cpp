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
    area("area",1,1,1,1), areaedge("areae",1,1,1,1),
    dx1("dx1",1,1,1,1), dx2("dx2",1,1,1,1), dx3("dx3",1,1,1,1),
    x1v("x1v",1,1), x2v("x2v",1,1), x3v("x3v",1,1),
    xx1f("xx1f",1,1), xx2f("xx2f",1,1), xx3f("xx3f",1,1), dxedge("dxe",1,1,1,1), dxface("dxf",1,1,1,1),
    sin_cell("sin_cell",1,1,1), cos_cell("cos_cell",1,1,1),
    sin_face_xi("sin_face_xi",1,1,1), cos_face_xi("cos_face_xi",1,1,1),
    sin_face_eta("sin_face_eta",1,1,1), cos_face_eta("cos_face_eta",1,1,1),
    x_ov_rD("x_ov_rD",1,1,1,1), y_ov_rC("y_ov_rC",1,1,1,1), z_ov_rE("z_ov_rC",1,1,1,1) {
        
  if (pmy_pack->pmesh->use_cubed_sphere || pmy_pack->pmesh->use_spherical_polar) {
    // Total number of MeshBlocks on this rank to be used in array dimensioning
    int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(volume, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(area.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(area.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(area.x3f, nmb, ncells3+1, ncells2, ncells1);
    Kokkos::realloc(dx1, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(dx2, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(dx3, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(x1v, nmb, ncells1);
    Kokkos::realloc(x2v, nmb, ncells2);
    Kokkos::realloc(x3v, nmb, ncells3);
    Kokkos::realloc(xx1f, nmb, ncells1+1);
    Kokkos::realloc(xx2f, nmb, ncells2+1);
    Kokkos::realloc(xx3f, nmb, ncells3+1);
    Kokkos::realloc(dxedge.x1e, nmb, ncells3+1, ncells2+1, ncells1);
    Kokkos::realloc(dxedge.x2e, nmb, ncells3+1, ncells2, ncells1+1);
    Kokkos::realloc(dxedge.x3e, nmb, ncells3, ncells2+1, ncells1+1);
    Kokkos::realloc(dxface.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(dxface.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(dxface.x3f, nmb, ncells3+1, ncells2, ncells1);
    Kokkos::realloc(areaedge.x1e, nmb, ncells3+1, ncells2+1, ncells1);
    Kokkos::realloc(areaedge.x2e, nmb, ncells3+1, ncells2, ncells1+1);
    Kokkos::realloc(areaedge.x3e, nmb, ncells3, ncells2+1, ncells1+1);
    Kokkos::realloc(x_ov_rD, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(y_ov_rC, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(z_ov_rE, nmb, ncells3, ncells2, ncells1);
  if (pmy_pack->pmesh->use_cubed_sphere) {
    // The gnomonic angles depend on the two PANEL-TANGENTIAL axes only -- xi = x2 and
    // eta = x3 -- so these are 3D arrays indexed (m,k,j), with no i (radial) extent.
    // +1 in the STAGGERED direction: the xi-face arrays live on x2 faces and are written
    // at j+1 for the last cell, the eta-face arrays on x3 faces and written at k+1.
    Kokkos::realloc(sin_cell, nmb, ncells3, ncells2);
    Kokkos::realloc(cos_cell, nmb, ncells3, ncells2);
    Kokkos::realloc(sin_face_xi, nmb, ncells3, ncells2+1);
    Kokkos::realloc(cos_face_xi, nmb, ncells3, ncells2+1);
    Kokkos::realloc(sin_face_eta, nmb, ncells3+1, ncells2);
    Kokkos::realloc(cos_face_eta, nmb, ncells3+1, ncells2);
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
//! \fn GnomonicCornerOmega / GnomonicSolidAngle
//! \brief EXACT solid angle of the gnomonic cell bounded by X = xl,xr and Y = yl,yr,
//! where X = tan(xi) and Y = tan(eta).
//!
//! A cell of the equiangular cubed sphere is a spherical quadrilateral bounded by four
//! GREAT circles (a curve of constant xi spans the plane through the origin containing
//! (1,X,0) and (0,0,1)), so its area is the classical solid angle of a rectangular
//! pyramid: one potential evaluated at the four corners.  Being a corner difference it
//! TELESCOPES -- subdividing a cell gives sub-areas summing to the parent to round-off.
//!
//! That property is the whole point.  The midpoint form this replaces,
//! dth_xi*dth_eta*sin_cell, is a perfectly good O(h^2) quadrature but is NOT additive
//! across a refinement level: four fine faces overshoot their parent by ~3e-4 relative
//! at nx2 = 32.  On a uniform grid that is invisible -- the radial momentum balance
//! cancels identically for ANY area values, because the geometric source term is built
//! as z_ov_rE = (area1r - area.x1f)/volume.  At a coarse/fine boundary it does not:
//! conservative flux correction hands the coarse cell sum(A_fine*F), while its source
//! term still uses A_coarse, and the leftover p*(sum(A_fine) - A_coarse)/V drives a
//! spurious radial acceleration out of a state that ought to be exactly static.
//! With the exact solid angle the two agree, and Omega cancels out of z_ov_rE entirely.

KOKKOS_INLINE_FUNCTION
Real GnomonicCornerOmega(Real X, Real Y) {
  return atan(X*Y/sqrt(1.0 + X*X + Y*Y));
}

KOKKOS_INLINE_FUNCTION
Real GnomonicSolidAngle(Real xl, Real xr, Real yl, Real yr) {
  return GnomonicCornerOmega(xr,yr) - GnomonicCornerOmega(xl,yr)
       - GnomonicCornerOmega(xr,yl) + GnomonicCornerOmega(xl,yl);
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
  (void)ie; (void)je; (void)ke;

  // AXES. x1 is the RADIAL direction, matching CoordSphericalPolar; the two
  // panel-tangential angles are xi = x2 and eta = x3. See the note in mesh.hpp.
  //
  // Cover the GHOST zones too: leaving sin/cos_cell, the areas and the volume at zero in
  // the ghosts breaks reconstruction and the momentum raise, both of which read them.
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;

  // THE RADIAL GRID STRETCH.  x1 is radial on this grid exactly as it is in spherical
  // polar, so the SAME 1-D map applies -- it is a function of the radial coordinate
  // alone and knows nothing about the angular grid, which is why it transplants
  // unchanged.  Only use_grid_stretch_theta stays refused (mesh.cpp): theta is not a
  // coordinate of this grid.
  const bool str_r_ = pmy_pack->pmesh->use_grid_stretch_r;
  const bool str_rp_ = pmy_pack->pmesh->use_grid_stretch_r_poly;
  const Real fstr_r_ = pmy_pack->pmesh->fStretchR;
  const Real rmin_ = pmy_pack->pmesh->mesh_size.x1min;
  const Real rmax_ = pmy_pack->pmesh->mesh_size.x1max;
  Real cpoly_[NSTRETCH_R_POLY];
  for (int n=0; n<NSTRETCH_R_POLY; ++n) { cpoly_[n] = pmy_pack->pmesh->fStretchRPoly[n]; }

  // THE 1-D COORDINATE ARRAYS.  They are ALLOCATED for the cubed sphere (the constructor
  // guard covers both curvilinear grids) but were never FILLED here, so every consumer
  // read ZEROS.  That is not a quiet inaccuracy: the hot-Jupiter pgen's hydrostatic outer
  // boundary forms (x1f[i+1] - x1v[i])/(x1f[i+1] - x1f[i]), which is 0/0, and it turned
  // the entire domain to NaN in ONE cycle while the run reported a normal exit.
  //
  // x1 is the RADIUS, exactly as in CoordSphericalPolar, so it gets the same volume-
  // centroid definition (Mignone 2014) and the two grids agree cell for cell.  x2 and x3
  // are the PANEL coordinates on [-1,1] and are NOT angles -- anything that needs an
  // angle must go through the panel chart, because the direction depends on which panel
  // the cell is on.
  {
    auto x1v_ = x1v;  auto xx1f_ = xx1f;
    auto x2v_ = x2v;  auto xx2f_ = xx2f;
    auto x3v_ = x3v;  auto xx3f_ = xx3f;
    par_for("cs_coord1d_1", DevExeSpace(), 0,nmb1, 0,n1m1,
    KOKKOS_LAMBDA(const int m, const int i) {
      Real r_l = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                                            size.d_view(m).x1max);
      Real r_r = LeftEdgeX(i-is+1, indcs.nx1, size.d_view(m).x1min,
                                              size.d_view(m).x1max);
      if (str_r_) {
        StretchR(fstr_r_, rmin_, rmax_, r_l);
        StretchR(fstr_r_, rmin_, rmax_, r_r);
      } else if (str_rp_) {
        StretchRPoly(cpoly_, rmin_, rmax_, r_l);
        StretchRPoly(cpoly_, rmin_, rmax_, r_r);
      }
      const Real q = r_l/r_r;
      x1v_(m,i) = 0.25*(q*q + 1.0)/((1.0/3.0)*(q*q + q + 1.0))*(r_r + r_l);
      xx1f_(m,i) = r_l;
      if (i == n1m1) { xx1f_(m,i+1) = r_r; }
    });
    par_for("cs_coord1d_2", DevExeSpace(), 0,nmb1, 0,n2m1,
    KOKKOS_LAMBDA(const int m, const int j) {
      x2v_(m,j) = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                               size.d_view(m).x2max);
      xx2f_(m,j) = LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                                              size.d_view(m).x2max);
      if (j == n2m1) {
        xx2f_(m,j+1) = LeftEdgeX(j-js+1, indcs.nx2, size.d_view(m).x2min,
                                                    size.d_view(m).x2max);
      }
    });
    par_for("cs_coord1d_3", DevExeSpace(), 0,nmb1, 0,n3m1,
    KOKKOS_LAMBDA(const int m, const int k) {
      x3v_(m,k) = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                               size.d_view(m).x3max);
      xx3f_(m,k) = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                              size.d_view(m).x3max);
      if (k == n3m1) {
        xx3f_(m,k+1) = LeftEdgeX(k-ks+1, indcs.nx3, size.d_view(m).x3min,
                                                    size.d_view(m).x3max);
      }
    });
  }

  auto volume_ = volume;
  auto area_ = area;
  auto dx1_ = dx1;
  auto dx2_ = dx2;
  auto dx3_ = dx3;
  auto dxedge_ = dxedge;
  auto dxface_ = dxface;
  auto areaedge_ = areaedge;
  auto sin_cell_ = sin_cell;
  auto cos_cell_ = cos_cell;
  auto sin_face_xi_ = sin_face_xi;
  auto cos_face_xi_ = cos_face_xi;
  auto sin_face_eta_ = sin_face_eta;
  auto cos_face_eta_ = cos_face_eta;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  par_for("cscoord", DevExeSpace(), 0,nmb1,0,n3m1,0,n2m1,0,n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    // --- radial (x1) ---
    // With a single radial cell CellCenterX/LeftEdgeX return the midpoint and the two
    // ends
    // of [x1min,x1max], i.e. a shell of finite thickness. Do NOT collapse this to
    // r_l = r_r: every face area below carries a factor (r_r^2 - r_l^2) and the volume a
    // factor (r_r^3 - r_l^3), so a zero-thickness shell makes them identically zero and
    // x_ov_rD / y_ov_rC / z_ov_rE then divide by that zero.
    Real r_c = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real r_l = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real r_r = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    // the radial stretch, applied exactly as CoordSphericalPolar applies it
    if (str_r_) {
      StretchR(fstr_r_, rmin_, rmax_, r_c);
      StretchR(fstr_r_, rmin_, rmax_, r_l);
      StretchR(fstr_r_, rmin_, rmax_, r_r);
    } else if (str_rp_) {
      StretchRPoly(cpoly_, rmin_, rmax_, r_c);
      StretchRPoly(cpoly_, rmin_, rmax_, r_l);
      StretchRPoly(cpoly_, rmin_, rmax_, r_r);
    }

    // --- angles: xi on x2, eta on x3 ---
    Real &x2min = size.d_view(m).x2min; Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min; Real &x3max = size.d_view(m).x3max;
    Real xi   = M_PI/4.0 * CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real eta  = M_PI/4.0 * CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    Real xil  = M_PI/4.0 * LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
    Real xir  = M_PI/4.0 * LeftEdgeX(j+1-js, indcs.nx2, x2min, x2max);
    Real etal = M_PI/4.0 * LeftEdgeX(k-ks, indcs.nx3, x3min, x3max);
    Real etar = M_PI/4.0 * LeftEdgeX(k+1-ks, indcs.nx3, x3min, x3max);

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

    sin_cell_(m,k,j) = sqrt(1.0 + SQR(x) + SQR(y)) / (C * D);
    cos_cell_(m,k,j) = - x * y / (C * D);
    sin_face_xi_(m,k,j) = sqrt(1.0 + SQR(xl) + SQR(y)) / (Cl * D);
    cos_face_xi_(m,k,j) = - xl * y / (Cl * D);
    sin_face_eta_(m,k,j) = sqrt(1.0 + SQR(x) + SQR(yl)) / (C * Dl);
    cos_face_eta_(m,k,j) = - x * yl / (C * Dl);
    Real sin_face_xir = sqrt(1.0 + SQR(xr) + SQR(y)) / (Cr * D);
    Real sin_face_etar = sqrt(1.0 + SQR(x) + SQR(yr)) / (C * Dr);
    if (j == n2m1) {
      sin_face_xi_(m,k,j+1) = sin_face_xir;
      cos_face_xi_(m,k,j+1) = - xr * y / (Cr * D);
    }
    if (k == n3m1) {
      sin_face_eta_(m,k+1,j) = sin_face_etar;
      cos_face_eta_(m,k+1,j) = - x * yr / (C * Dr);
    }

    // --- angular face widths ---
    Real dth_xi = acos((1.0 + xl*xr + y*y) / sqrt(1.0 + xl*xl + y*y) / sqrt(1.0 + xr*xr + y*y));
    Real dth_xil = acos((1.0 + xl*xr + yl*yl) / sqrt(1.0 + xl*xl + yl*yl) / sqrt(1.0 + xr*xr + yl*yl));
    Real dth_xir = acos((1.0 + xl*xr + yr*yr) / sqrt(1.0 + xl*xl + yr*yr) / sqrt(1.0 + xr*xr + yr*yr));
    Real dth_eta = acos((1.0 + x*x + yl*yr) / sqrt(1.0 + x*x + yl*yl) / sqrt(1.0 + x*x + yr*yr));
    Real dth_etal = acos((1.0 + xl*xl + yl*yr) / sqrt(1.0 + xl*xl + yl*yl) / sqrt(1.0 + xl*xl + yr*yr));
    Real dth_etar = acos((1.0 + xr*xr + yl*yr) / sqrt(1.0 + xr*xr + yl*yl) / sqrt(1.0 + xr*xr + yr*yr));

    // --- radial faces (x1) ---
    // EXACT solid angle, not dth_xi*dth_eta*sin_cell: only the exact form is additive
    // across a refinement level.  See GnomonicSolidAngle above.
    const Real domega = GnomonicSolidAngle(xl, xr, yl, yr);
    area_.x1f(m,k,j,i) = SQR(r_l) * domega;
    Real area1r = SQR(r_r) * domega;
    if (i == n1m1) area_.x1f(m,k,j,i+1) = area1r;

    // --- xi faces (x2) ---
    area_.x2f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etal;
    Real area2r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etar;
    if (j == n2m1) area_.x2f(m,k,j+1,i) = area2r;

    // --- eta faces (x3) ---
    area_.x3f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xil;
    Real area3r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xir;
    if (k == n3m1) area_.x3f(m,k+1,j,i) = area3r;

    // --- volume_ (exact solid angle x exact radial shell) ---
    volume_(m,k,j,i) = 1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l) * domega;
    dx1_(m,k,j,i) = r_r - r_l;
    dx2_(m,k,j,i) = r_c * dth_xi;
    dx3_(m,k,j,i) = r_c * dth_eta;

    // --- edge lengths, consumed by the curvilinear CT update in mhd_ct.cpp -------------
    // Without these the CT branch multiplies every EMF by an uninitialised (zero) edge
    // length, so B is left exactly unchanged and the field silently never evolves.
    // Each edge takes its length from the SAME arc the adjoining face area uses, so the
    // two stay consistent by construction:
    //   x1e is radial, so it is just r_r - r_l and depends on neither angle;
    //   x2e runs along xi  and is evaluated at an ETA face, hence dth_xi{l,r} -- the same
    //     arcs that give area_.x3f;
    //   x3e runs along eta and is evaluated at a XI  face, hence dth_eta{l,r} -- the same
    //     arcs that give area_.x2f.
    // The i / j / k == n*m1 clauses fill the far edge of the last cell, exactly as
    // CoordSphericalPolar does.
    dxedge_.x1e(m,k,j,i) = r_r - r_l;
    if (j == n2m1) dxedge_.x1e(m,k,j+1,i) = r_r - r_l;
    if (k == n3m1) dxedge_.x1e(m,k+1,j,i) = r_r - r_l;
    if (j == n2m1 && k == n3m1) dxedge_.x1e(m,k+1,j+1,i) = r_r - r_l;

    dxedge_.x2e(m,k,j,i) = r_l * dth_xil;
    if (i == n1m1) dxedge_.x2e(m,k,j,i+1) = r_r * dth_xil;
    if (k == n3m1) dxedge_.x2e(m,k+1,j,i) = r_l * dth_xir;
    if (i == n1m1 && k == n3m1) dxedge_.x2e(m,k+1,j,i+1) = r_r * dth_xir;

    dxedge_.x3e(m,k,j,i) = r_l * dth_etal;
    if (i == n1m1) dxedge_.x3e(m,k,j,i+1) = r_r * dth_etal;
    if (j == n2m1) dxedge_.x3e(m,k,j+1,i) = r_l * dth_etar;
    if (i == n1m1 && j == n2m1) dxedge_.x3e(m,k,j+1,i+1) = r_r * dth_etar;

    // --- DUAL mesh: edge-centred areas and centre-to-centre lengths ------------------
    // Consumed by the resistive current density (diffusion/current_density.hpp), which
    // applies Stokes' theorem on the loop joining the four cell CENTRES around an edge.
    // Every side of that loop passes through a face centre, so its length is the
    // centre-to-centre distance stored in dxface, and the area it encloses is areaedge_.
    // The x1 loop lies in the tangent surface, where the two basis vectors are NOT
    // orthogonal, so its area carries the sin of the angle between them -- exactly as
    // area_.x1f does. The x2 and x3 loops each contain rhat, which IS orthogonal to both
    // tangent directions, so they carry no such factor.
    Real r_cm = CellCenterX(i-1-is, indcs.nx1, size.d_view(m).x1min,
                                            size.d_view(m).x1max);
    Real xim  = M_PI/4.0 * CellCenterX(j-1-js, indcs.nx2, x2min, x2max);
    Real etam = M_PI/4.0 * CellCenterX(k-1-ks, indcs.nx3, x3min, x3max);
    Real xm = tan(xim);
    Real ym = tan(etam);
    // angle subtended between two directions on the same panel, at fixed other angle
    auto arc_xi = [](Real xa, Real xb, Real yy) {
      return acos((1.0 + xa*xb + yy*yy)/sqrt(1.0 + xa*xa + yy*yy)
                                       /sqrt(1.0 + xb*xb + yy*yy));
    };
    auto arc_eta = [](Real ya, Real yb, Real xx) {
      return acos((1.0 + xx*xx + ya*yb)/sqrt(1.0 + xx*xx + ya*ya)
                                       /sqrt(1.0 + xx*xx + yb*yb));
    };
    if (i > 0) dxface_.x1f(m,k,j,i) = r_c - r_cm;
    if (j > 0) dxface_.x2f(m,k,j,i) = r_c * arc_xi(xm, x, y);
    if (k > 0) dxface_.x3f(m,k,j,i) = r_c * arc_eta(ym, y, x);
    if (j > 0 && k > 0) {
      areaedge_.x1e(m,k,j,i) = SQR(r_c) * arc_xi(xm, x, yl) * arc_eta(ym, y, xl)
                              * sqrt(1.0 + SQR(xl) + SQR(yl)) / (Cl * Dl);
    }
    if (i > 0 && k > 0) {
      areaedge_.x2e(m,k,j,i) = 0.5*(SQR(r_c) - SQR(r_cm)) * arc_eta(ym, y, x);
    }
    if (i > 0 && j > 0) {
      areaedge_.x3e(m,k,j,i) = 0.5*(SQR(r_c) - SQR(r_cm)) * arc_xi(xm, x, y);
    }

    // Geometric coefficients consumed by SrcTermsGnomonicEquiangle. z_ov_rE is the RADIAL
    // one, matching its meaning in CoordSphericalPolar/SrcTermsSphericalPolar.
    x_ov_rD_(m,k,j,i) = (area2r * sin_face_xir
                        - area_.x2f(m,k,j,i) * sin_face_xi_(m,k,j)) / volume_(m,k,j,i);
    y_ov_rC_(m,k,j,i) = (area3r * sin_face_etar
                        - area_.x3f(m,k,j,i) * sin_face_eta_(m,k,j)) / volume_(m,k,j,i);
    z_ov_rE_(m,k,j,i) = (area1r - area_.x1f(m,k,j,i)) / volume_(m,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleRaiseVel
//! \brief Recompute the primitive velocity and internal energy from the conserved state
//! on the cubed sphere, using the metric of the non-orthogonal gnomonic tangent basis.
//!
//! WHY THIS EXISTS. On the cubed sphere the two angular basis vectors are unit vectors
//! separated by an angle theta with cos(theta) = cos_cell = -xy/(CD), so the metric on
//! that basis is g = [[1, c],[c, 1]] and is NOT the identity. The rest of the scheme is
//! built on a consistent split:
//!   * GnomonicEquianglePrimFaceX* rotate w0's velocity into a locally ORTHONORMAL frame
//!     before the Riemann solve, treating (IVX,IVY) as CONTRAVARIANT components v^1, v^2;
//!   * GnomonicEquiangleFluxX* rotate the returned flux back and then LOWER the index, so
//!     what the update accumulates in u0(IM1,IM2) is the COVARIANT momentum rho*v_i;
//!   * SrcTermsGnomonicEquiangle reads w0 as contravariant and forms v_i = v^i + c v^j
//!     itself, i.e. it agrees with both of the above.
//! The generic EOS does not, and cannot: SingleC2P_IdealHyd sets v = m/rho and subtracts
//! a kinetic energy 0.5*(m.m)/rho, both of which assume an orthonormal basis. Applying it
//! unchanged leaves the primitive velocity holding COVARIANT components where every
//! consumer expects contravariant ones, and the internal energy short by the cross term
//! rho c v^1 v^2. The resulting scheme is inconsistent, not merely inaccurate: a rigid
//! rotation, which is an exact steady solution of the Euler equations, drifts by an
//! amount that does NOT decrease under grid refinement.
//!
//! This runs immediately after ConsToPrim and overwrites what it wrote for IVX/IVY/IVZ
//! and IEN. Density and the scalars are left exactly as ConsToPrim left them. The ENERGY
//! FLOORS are re-applied here rather than inherited: ConsToPrim floored an internal
//! energy computed with an orthonormal kinetic energy, and the metric cross term added
//! below moves it, so the floor has to be retested against the corrected value (and u0
//! updated to match). See the note at the floor itself.
//! Non-relativistic ideal/general hydro only -- the relativistic inversions are coupled
//! and would need the metric inside the root find, not after it.

void Coordinates::GnomonicEquiangleRaiseVel(DvceArray5D<Real> &u0,
    DvceArray5D<Real> &w0, const EOS_Data &eos_data, const DvceArray5D<Real> &wder,
    const DvceArray4D<Real> &wtemp, const int il, const int iu, const int jl,
    const int ju, const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  // A GENERAL EOS CACHES p, Gamma_1 AND T, AND ConsToPrim SOLVED FOR THEM BEFORE THIS
  // CORRECTION.  Its inversion subtracts an ORTHONORMAL kinetic energy, so the internal
  // energy it inverted is wrong by the metric cross term -- exactly the amount corrected
  // below -- and wder/wtemp are written NOWHERE ELSE.  Left stale they are read by the
  // fluxes, the timestep, FOFC, prolongation and the geometric source term itself, so the
  // error is O(cos_cell) x tangential KE: zero on the panel midlines and largest at the
  // panel corners.  Re-solve them here from the CORRECTED internal energy, warm-started
  // on the temperature ConsToPrim already found, which makes the root find cheap.
  const bool gen_ = eos_data.IsGeneral();
  auto eos_ = eos_data;
  auto wder_ = wder;
  auto wtemp_ = wtemp;

  par_for("cs_raisev", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real det = 1.0 - c*c;
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    const Real m2 = u0(m,IM2,k,j,i);   // xi
    const Real m3 = u0(m,IM3,k,j,i);   // eta
    // v^i = g^{ij} m_j / rho, with the metric acting on the ANGULAR pair only
    const Real v1 = m1/d;
    const Real v2 = (m2 - c*m3)/(d*det);
    const Real v3 = (m3 - c*m2)/(d*det);
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    // KE = 0.5 rho g_ij v^i v^j = 0.5 m_i v^i, which is the cross-term-correct form.
    Real eint = u0(m,IEN,k,j,i) - 0.5*(m1*v1 + m2*v2 + m3*v3);
    // ---------------------------------------------------------------------------------
    // RE-APPLY THE FLOORS.  ConsToPrim floored the state it inverted, but that state
    // carried an ORTHONORMAL kinetic energy; the metric cross term above MOVES the
    // internal energy, so a cell ConsToPrim left comfortably above the floor can land
    // below it -- or below zero -- here.  Leaving that unfloored hands a non-positive
    // internal energy straight to the tabulated inversion below, which is undefined
    // there and returns NaN in T, p and Gamma_1; the NaN then leaves the cell through
    // the reconstruction stencil and takes the whole grid down within ~100 cycles.
    // The sequence deliberately mirrors SingleC2P_GeneralHyd, including its guard
    // against inverting a non-positive energy.  NOTE the pressure floor alone is NOT
    // sufficient under a tabulated EOS: at upper-atmosphere densities e(d,pfloor) lies
    // far below the table's lowest temperature, so it is the TEMPERATURE floor that
    // actually keeps the lookup in range.  u0 is updated in step, exactly as ConsToPrim
    // updates cons when a floor fires, so the conserved energy cannot keep sinking and
    // re-trip the floor on every cycle.
    if (gen_) {
      Real temp = -1.0, pnew = 0.0, g1new = 0.0;
      const bool e_positive = (eint > 0.0);
      bool stale = !e_positive;
      if (e_positive) {
        eos_.TemperaturePressureGamma1(d, eint, wtemp_(m,k,j,i), temp, pnew, g1new);
      }
      if (!e_positive || pnew < eos_.pfloor) {
        eint = eos_.EnergyFromPressure(d, eos_.pfloor, temp);
        stale = true;
      }
      if (temp < eos_.tfloor) {
        eint = eos_.EnergyFromTemperature(d, eos_.tfloor);
        temp = eos_.tfloor;
        stale = true;
      }
      if (stale) {
        eos_.PressureAndGamma1(d, eint, temp, pnew, g1new);
        u0(m,IEN,k,j,i) = eint + 0.5*(m1*v1 + m2*v2 + m3*v3);
      }
      wder_(m,IDPR,k,j,i) = pnew;
      wder_(m,IDG1,k,j,i) = g1new;
      wtemp_(m,k,j,i) = temp;
    } else {
      const Real eold = eint;
      eos_.ApplyEnergyFloor(d, eint);
      if (eint != eold) { u0(m,IEN,k,j,i) = eint + 0.5*(m1*v1+m2*v2+m3*v3); }
    }
    w0(m,IEN,k,j,i) = eint;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleRaiseVelMHD
//! \brief The MHD counterpart of GnomonicEquiangleRaiseVel: it does everything that
//! routine does, and additionally puts the CELL-CENTRED FIELD into an orthonormal frame
//! and subtracts the correct magnetic energy.
//!
//! WHY THE FIELD NEEDS THIS TOO. b0.x*f holds the flux density B.nhat through its own
//! face -- the convention forced by the curvilinear CT update in mhd_ct.cpp, which uses
//! area/dxedge -- so the plain face average that ConsToPrim forms is
//! (B.rhat, B.nhat_xi, B.nhat_eta). Those three directions are NOT mutually orthogonal:
//! nhat_xi.nhat_eta = -c. Summing their squares therefore does not give |B|^2, and the
//! error is O(1) near a panel corner where c is largest. Every consumer that squares and
//! sums bcc -- the magnetic energy here, PrimToCons, the fast speed in mhd_newdt -- is
//! wrong by that amount.
//!
//! The fix is to store bcc in the ORTHONORMAL frame {rhat, e_xi, (e_eta - c e_xi)/s}
//! instead, which costs one rotation of a single component and makes every one of those
//! sums correct with no further change. It also SIMPLIFIES the flux path: the x1 and x3
//! sweeps then need no field rotation at all, and only the x2 sweep does.
//!
//! Recomputing bcc from the faces rather than transforming what ConsToPrim wrote keeps
//! this idempotent. As in the hydro version, the floors are left exactly as ConsToPrim
//! applied them -- it tested them against an internal energy that lacked both the metric
//! cross term and the corrected magnetic energy.

void Coordinates::GnomonicEquiangleRaiseVelMHD(DvceArray5D<Real> &u0,
    const DvceFaceFld4D<Real> &b0, DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
    const EOS_Data &eos_data, const DvceArray5D<Real> &wder,
    const DvceArray4D<Real> &wtemp, const int il, const int iu, const int jl,
    const int ju, const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  auto &sin_cell_ = sin_cell;
  // See the note in GnomonicEquiangleRaiseVel: a general EOS's cached p, Gamma_1 and T
  // were solved for BEFORE this correction and are written nowhere else, so they have to
  // be re-solved from the corrected internal energy.  Here the magnetic energy is also
  // rebuilt, since ConsToPrim formed it from the NON-ORTHOGONAL face-normal triple.
  const bool gen_ = eos_data.IsGeneral();
  auto eos_ = eos_data;
  auto wder_ = wder;
  auto wtemp_ = wtemp;

  par_for("cs_raisev_mhd", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real sn = sin_cell_(m,k,j);
    const Real det = 1.0 - c*c;

    // cell-centred field, then into the orthonormal frame. The eta slot already IS the
    // third axis of that frame (B.(e_eta - c e_xi)/s = s B^eta), so only the xi slot
    // moves: B.e_xi = B^xi + c B^eta = (b_xi + c b_eta)/s.
    const Real bx = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    const Real by = (by_n + c*bz_n)/sn;
    const Real bz = bz_n;
    bcc0(m,IBX,k,j,i) = bx;
    bcc0(m,IBY,k,j,i) = by;
    bcc0(m,IBZ,k,j,i) = bz;

    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    const Real m2 = u0(m,IM2,k,j,i);   // xi
    const Real m3 = u0(m,IM3,k,j,i);   // eta
    const Real v1 = m1/d;
    const Real v2 = (m2 - c*m3)/(d*det);
    const Real v3 = (m3 - c*m2)/(d*det);
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    // the frame is orthonormal, so the magnetic energy IS the sum of squares.
    // Keep the two subtractions separate and in this order: folding them into a single
    // (kinetic + magnetic) sum re-associates the rounding and perturbs every ideal-EOS
    // cubed-sphere answer in the last bits, for nothing.
    Real eint = u0(m,IEN,k,j,i) - 0.5*(m1*v1 + m2*v2 + m3*v3)
                                - 0.5*(bx*bx + by*by + bz*bz);
    // Re-apply the floors to the CORRECTED internal energy; see the extended note in
    // GnomonicEquiangleRaiseVel above. Without this the tabulated inversion below is
    // handed a non-positive energy and returns NaN.
    if (gen_) {
      Real temp = -1.0, pnew = 0.0, g1new = 0.0;
      const bool e_positive = (eint > 0.0);
      bool stale = !e_positive;
      if (e_positive) {
        eos_.TemperaturePressureGamma1(d, eint, wtemp_(m,k,j,i), temp, pnew, g1new);
      }
      if (!e_positive || pnew < eos_.pfloor) {
        eint = eos_.EnergyFromPressure(d, eos_.pfloor, temp);
        stale = true;
      }
      if (temp < eos_.tfloor) {
        eint = eos_.EnergyFromTemperature(d, eos_.tfloor);
        temp = eos_.tfloor;
        stale = true;
      }
      if (stale) {
        eos_.PressureAndGamma1(d, eint, temp, pnew, g1new);
        u0(m,IEN,k,j,i) = eint + 0.5*(m1*v1 + m2*v2 + m3*v3)
                                 + 0.5*(bx*bx + by*by + bz*bz);
      }
      wder_(m,IDPR,k,j,i) = pnew;
      wder_(m,IDG1,k,j,i) = g1new;
      wtemp_(m,k,j,i) = temp;
    } else {
      const Real eold = eint;
      eos_.ApplyEnergyFloor(d, eint);
      if (eint != eold) {
        u0(m,IEN,k,j,i) = eint + 0.5*(m1*v1 + m2*v2 + m3*v3)
                               + 0.5*(bx*bx + by*by + bz*bz);
      }
    }
    w0(m,IEN,k,j,i) = eint;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleLowerMom
//! \brief The inverse of GnomonicEquiangleRaiseVel: build the conserved state from
//! primitives on the cubed sphere. Use this wherever a problem generator would otherwise
//! call EquationOfState::PrimToCons, which assumes an orthonormal basis.

void Coordinates::GnomonicEquiangleLowerMom(const DvceArray5D<Real> &w0,
    DvceArray5D<Real> &u0, const int il, const int iu, const int jl, const int ju,
    const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;

  par_for("cs_lowerm", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real d = w0(m,IDN,k,j,i);
    const Real v1 = w0(m,IVX,k,j,i);
    const Real v2 = w0(m,IVY,k,j,i);
    const Real v3 = w0(m,IVZ,k,j,i);
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;            // m_i = rho g_ij v^j
    u0(m,IM2,k,j,i) = d*(v2 + c*v3);
    u0(m,IM3,k,j,i) = d*(v3 + c*v2);
    u0(m,IEN,k,j,i) = w0(m,IEN,k,j,i)
                    + 0.5*d*(v1*v1 + v2*v2 + v3*v3 + 2.0*c*v2*v3);
  });
  return;
}

// Hydro and MHD share one implementation; the MHD entry point differs only in adding the
// Maxwell stress, so the two are kept in the same kernel rather than duplicated the way
// SrcTermsSphericalPolar{Hydro,MHD} are.
void Coordinates::SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
    const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {
  SrcTermsGnomonicEquiangleImpl(w0, w0, false, wder, uflx, eos_data, bdt, u0);
  return;
}

void Coordinates::SrcTermsGnomonicEquiangleMHD(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc0, const DvceArray5D<Real> &wder,
    const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
    DvceArray5D<Real> &u0) {
  SrcTermsGnomonicEquiangleImpl(w0, bcc0, true, wder, uflx, eos_data, bdt, u0);
  return;
}

void Coordinates::SrcTermsGnomonicEquiangleImpl(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc0, const bool is_mhd,
    const DvceArray5D<Real> &wder, const DvceFaceFld5D<Real> uflx,
    const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
    
  // by-value copy of the EOS, capturable in the device lambda below
  auto eos_ = eos_data;
  // The live gas pressure is read from the derived-variable cache rather than recomputed:
  // p(d,e) is a root find for a general EOS, and ConsToPrim has already evaluated it once
  // for this very state. An ideal gas allocates no wder array, and (gamma-1)*e is free.
  // This routine serves both Hydro and MHD, so it cannot pick the module's wder itself --
  // a two-fluid run has both -- and takes it as an argument instead.
  const bool gen_ = eos_.IsGeneral();
  auto &wder_ = wder;
  auto &bcc_ = bcc0;
  const bool mhd_ = is_mhd;

  auto volume_ = volume;
  auto area_ = area;
  auto sin_cell_ = sin_cell;
  auto cos_cell_ = cos_cell;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  par_for("cssrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      
    // x1 is RADIAL; xi = x2 and eta = x3 are the panel-tangential angles.
    Real radius = CellCenterX(i-is, indcs.nx1,
                              size.d_view(m).x1min, size.d_view(m).x1max);
    Real dr = size.d_view(m).dx1;

    Real v1 = w0(m,IVX,k,j,i);   // radial
    Real v2 = w0(m,IVY,k,j,i);   // xi
    Real v3 = w0(m,IVZ,k,j,i);   // eta
    Real pr = gen_ ? wder_(m,IDPR,k,j,i)
                            : eos_.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
    Real rho = w0(m,IDN,k,j,i);

    Real sine = sin_cell_(m,k,j);
    Real cosine = cos_cell_(m,k,j);
    Real sine2 = SQR(sine);
    // lower the index on the two ANGULAR components with g = [[1,c],[c,1]]
    Real v_2 = v2 + v3 * cosine;
    Real v_3 = v3 + v2 * cosine;

    // MHD adds the Maxwell stress: T^ij = rho v^i v^j + (p + B^2/2) g^ij - B^i B^j, so
    // every p that multiplies a metric derivative below becomes the TOTAL pressure and
    // every Reynolds stress becomes Reynolds minus Maxwell. This is the same substitution
    // that separates SrcTermsSphericalPolarMHD from ...Hydro.
    Real ptot = pr;
    Real b_tsq = 0.0, bb_2 = 0.0, bb_3 = 0.0;
    if (mhd_) {
      // bcc is the ORTHONORMAL triple (B.rhat, B.e_xi, B.(e_eta - c e_xi)/s) -- see
      // GnomonicEquiangleRaiseVelMHD -- so |B|^2 really is the sum of squares. The
      // Maxwell partners of rho*v3^2*sine2 and rho*v2^2*sine2 are (s*B^eta)^2 and
      // (s*B^xi)^2, and B^eta = b3/s while s*B^xi = s*b2 - c*b3.
      const Real b1 = bcc_(m,IBX,k,j,i);
      const Real b2 = bcc_(m,IBY,k,j,i);
      const Real b3 = bcc_(m,IBZ,k,j,i);
      b_tsq = b2*b2 + b3*b3;
      bb_2 = SQR(sine*b2 - cosine*b3);
      bb_3 = b3*b3;
      ptot += 0.5*(b1*b1 + b_tsq);
    }

    // radial source: (A_out - A_in)/V * (p + angular kinetic energy), the same form as
    // SrcTermsSphericalPolar's src1. For MHD the tangential Maxwell stress cancels all
    // but the RADIAL magnetic pressure, leaving z_ov_rE*(p + B_r^2/2 + rho|v_t|^2/2).
    Real src1 = z_ov_rE_(m,k,j,i) * (ptot + 0.5*(rho*(v2*v_2+v3*v_3) - b_tsq));
    Real src2 = x_ov_rD_(m,k,j,i) * (ptot + rho*SQR(v3)*sine2 - bb_3);
    Real src3 = y_ov_rC_(m,k,j,i) * (ptot + rho*SQR(v2)*sine2 - bb_2);

    // the radial-face flux of each ANGULAR momentum, which the r-dependence of the
    // angular basis converts into a source
    src2 -= dr/2.0/radius * (uflx.x1f(m,IM2,k,j,i)*area_.x1f(m,k,j,i)
             + uflx.x1f(m,IM2,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
    src3 -= dr/2.0/radius * (uflx.x1f(m,IM3,k,j,i)*area_.x1f(m,k,j,i)
             + uflx.x1f(m,IM3,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);

    u0(m,IM1,k,j,i) += src1*bdt;
    u0(m,IM2,k,j,i) += src2*bdt;
    u0(m,IM3,k,j,i) += src3*bdt;
    
  });
}

void Coordinates::CoordSphericalPolar() {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  // Copy everything the lambda needs out of the Mesh before the kernel: dereferencing a
  // host-side pointer on device only happens to work here under unified memory.
  const bool str_r = pmy_pack->pmesh->use_grid_stretch_r;
  const bool str_rp = pmy_pack->pmesh->use_grid_stretch_r_poly;
  const Real fstr_r = pmy_pack->pmesh->fStretchR;
  const Real rmin = pmy_pack->pmesh->mesh_size.x1min;
  const Real rmax = pmy_pack->pmesh->mesh_size.x1max;
  Real cpoly[NSTRETCH_R_POLY];
  for (int n=0; n<NSTRETCH_R_POLY; ++n) cpoly[n] = pmy_pack->pmesh->fStretchRPoly[n];

  auto volume_ = volume;
  auto area_ = area;
  auto dx1_ = dx1;
  auto dx2_ = dx2;
  auto dx3_ = dx3;
  auto x1v_ = x1v;
  auto x2v_ = x2v;
  auto x3v_ = x3v;
  auto xx1f_ = xx1f;
  auto xx2f_ = xx2f;
  auto xx3f_ = xx3f;
  auto dxedge_ = dxedge;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  const bool use_stretch_theta = pmy_pack->pmesh->use_grid_stretch_theta;
  const Real f_stretch_theta = pmy_pack->pmesh->fStretchTheta;
  par_for("spcoord", DevExeSpace(), 0,nmb1,0,n3m1,0,n2m1,0,n1m1,//ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    // --- radial ---
    Real r_c, r_l, r_r;
    r_c  = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    r_l  = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    r_r  = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    // stretch the radial grid
    if (str_r) {
      StretchR(fstr_r,rmin,rmax,r_c);
      StretchR(fstr_r,rmin,rmax,r_l);
      StretchR(fstr_r,rmin,rmax,r_r);
    } else if (str_rp) {
      StretchRPoly(cpoly,rmin,rmax,r_c);
      StretchRPoly(cpoly,rmin,rmax,r_l);
      StretchRPoly(cpoly,rmin,rmax,r_r);
    }

    // --- angles ---
    Real theta  = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real thetal = LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real thetar = LeftEdgeX(j+1-js, indcs.nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    // stretch the radial grid
    if (use_stretch_theta) {
      Real factor_stretch = f_stretch_theta;
      StretchTheta(factor_stretch,theta);
      StretchTheta(factor_stretch,thetal);
      StretchTheta(factor_stretch,thetar);
    }
    Real phi    = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phil   = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phir   = LeftEdgeX(k+1-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);

    Real sintl  = sin(thetal);
    Real sintr  = sin(thetar);
    Real sinl   = fabs(sintl);
    Real sinr   = fabs(sintr);
    Real sinc   = fabs(sin(theta));
    Real cosl   = cos(thetal);
    Real cosr   = cos(thetar);
      
    if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar) {
      if (j == js) {
        thetal = 0.0;
        sintl = 0.0;
        sinl = 0.0;
        cosl = 1.0;
      }
      if (j == js-1) {
        thetar = 0.0;
        sintr = 0.0;
        sinr = 0.0;
        cosr = 1.0;
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar) {
      if (j == je) {
        thetar = M_PI;
        sintr = 0.0;
        sinr = 0.0;
        cosr = -1.0;
      }
      if (j == je+1) {
        thetal = M_PI;
        sintl = 0.0;
        sinl = 0.0;
        cosl = -1.0;
      }
    }
      
    // --- radial edge lengths ---
    dxedge_.x1e(m,k,j,i) = r_r - r_l;
    if (j == n2m1) dxedge_.x1e(m,k,j+1,i) = r_r - r_l;
    if (k == n3m1) dxedge_.x1e(m,k+1,j,i) = r_r - r_l;
    if (j == n2m1 && k == n3m1) dxedge_.x1e(m,k+1,j+1,i) = r_r - r_l;
      
    // --- theta edge lengths ---
    dxedge_.x2e(m,k,j,i) = r_l * (thetar-thetal);
    if (i == n1m1) dxedge_.x2e(m,k,j,i+1) = r_r * (thetar-thetal);
    if (k == n3m1) dxedge_.x2e(m,k+1,j,i) = r_l * (thetar-thetal);
    if (i == n1m1 && k == n3m1) dxedge_.x2e(m,k+1,j,i+1) = r_r * (thetar-thetal);
      
    // --- phi edge lengths ---
    dxedge_.x3e(m,k,j,i) = r_l * sinl * (phir-phil);
    if (i == n1m1) dxedge_.x3e(m,k,j,i+1) = r_r * sinl * (phir-phil);
    if (j == n2m1) dxedge_.x3e(m,k,j+1,i) = r_l * sinr * (phir-phil);
    if (i == n1m1 && j == n2m1) dxedge_.x3e(m,k,j+1,i+1) = r_r * sinr * (phir-phil);

    // --- radial faces ---
    area_.x1f(m,k,j,i) = SQR(r_l) * fabs(cosl-cosr) * (phir-phil);
    Real area1r = SQR(r_r) * fabs(cosl-cosr) * (phir-phil);
    if (i == n1m1) area_.x1f(m,k,j,i+1) = area1r;
    xx1f_(m,i) = r_l;
    if (i == n1m1) xx1f_(m,i+1) = r_r;

    // --- theta faces ---
    area_.x2f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * sinl * (phir-phil);
    Real area2r = 0.5 * (SQR(r_r)-SQR(r_l)) * sinr * (phir-phil);
    if (j == n2m1) area_.x2f(m,k,j+1,i) = area2r;
    xx2f_(m,j) = thetal;
    if (j == n2m1) xx2f_(m,j+1) = thetar;

    // --- phi faces ---
    area_.x3f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * (thetar-thetal);
    Real area3r = 0.5 * (SQR(r_r)-SQR(r_l)) * (thetar-thetal);
    if (k == n3m1) area_.x3f(m,k+1,j,i) = area3r;
    xx3f_(m,k) = phil;
    if (k == n3m1) xx3f_(m,k+1) = phir;

    // --- volume ---

    volume_(m,k,j,i) = 1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l) * fabs(cosl-cosr) * (phir-phil);
      
    // Mignone 2014 correction
    x1v_(m,i) = 1.0/4.0*(SQR(r_l/r_r)+1.0) / (1.0/3.0*(SQR(r_l/r_r)+(r_l/r_r)+1.0)) *(r_r+r_l); // 1.0/4.0*(SQR(SQR(r_r))-SQR(SQR(r_l))) / (1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l)); //0.5*(r_r+r_l); //
    x2v_(m,j) = -((sintr-thetar*cosr) - (sintl-thetal*cosl)) / (cosr-cosl); // 0.5*(thetar+thetal); // 
    x3v_(m,k) = 0.5*(phir+phil);
    dx2_(m,k,j,i) = x1v_(m,i) * (thetar-thetal);
    dx3_(m,k,j,i) = x1v_(m,i) * fabs(sin(x2v_(m,j))) * (phir-phil);
    dx1_(m,k,j,i) = r_r-r_l;
      
    z_ov_rE_(m,k,j,i) = (area1r - area_.x1f(m,k,j,i)) / volume_(m,k,j,i);
    x_ov_rD_(m,k,j,i) = (area2r - area_.x2f(m,k,j,i)) / volume_(m,k,j,i);
    y_ov_rC_(m,k,j,i) = (sinr-sinl)/(sinr+sinl);
  });
  auto dxface_ = dxface;
  auto areaedge_ = areaedge;
  par_for("spcoord_add", DevExeSpace(), 0,nmb1,0,n3m1,0,n2m1,0,n1m1,//ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (i > 0) dxface_.x1f(m,k,j,i) = x1v_(m,i) - x1v_(m,i-1);
    if (j > 0) dxface_.x2f(m,k,j,i) = x1v_(m,i)*(x2v_(m,j)-x2v_(m,j-1));
    if (k > 0) {
      dxface_.x3f(m,k,j,i) = x1v_(m,i)*fabs(sin(x2v_(m,j)))*(x3v_(m,k)-x3v_(m,k-1));
    }
    if (j > 0 && k > 0) areaedge_.x1e(m,k,j,i) = SQR(x1v_(m,i)) * fabs(cos(x2v_(m,j))-cos(x2v_(m,j-1))) * (x3v_(m,k)-x3v_(m,k-1));
    if (k > 0 && i > 0) areaedge_.x2e(m,k,j,i) = 0.5 * (SQR(x1v_(m,i))-SQR(x1v_(m,i-1))) * fabs(sin(x2v_(m,j))) * (x3v_(m,k)-x3v_(m,k-1));
    if (i > 0 && j > 0) areaedge_.x3e(m,k,j,i) = 0.5 * (SQR(x1v_(m,i))-SQR(x1v_(m,i-1))) * (x2v_(m,j)-x2v_(m,j-1));
  });
}

void Coordinates::SrcTermsSphericalPolarHydro(const DvceArray5D<Real> &w0,
    const DvceArray4D<Real> &pwb, const DvceFaceFld5D<Real> uflx,
    const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
 
  // by-value copy of the EOS, capturable in the device lambda below
  auto eos_ = eos_data;
  // The live gas pressure is read from the derived-variable cache rather than recomputed:
  // p(d,e) is a root find for a general EOS, and ConsToPrim has already evaluated it once
  // for this very state. An ideal gas allocates no wder array, and (gamma-1)*e is free.
  const bool gen_ = eos_.IsGeneral();
  auto &wder_ = pmy_pack->phydro->wder;

  auto volume_ = volume;
  auto area_ = area;
  auto xx1f_ = xx1f;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  const bool use_wb_static_ = pmy_pack->phydro->use_wellbalance_static;
  par_for("spsrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      
    Real r_l = xx1f_(m,i);
    Real r_r = xx1f_(m,i+1);
    Real factor = (r_r-r_l)/(r_r+r_l);
      
    Real v1 = w0(m,IVX,k,j,i);
    Real v2 = w0(m,IVY,k,j,i);
    Real v3 = w0(m,IVZ,k,j,i);
    Real rho = w0(m,IDN,k,j,i);
    Real pr = gen_ ? wder_(m,IDPR,k,j,i)
                            : eos_.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
    // subtract the background pressure, precomputed once by SetWbBackgroundPressure()
    // -- for a general EOS it depends on the background DENSITY as well as its internal
    // energy, so it cannot be written as (gamma-1)*e_bg and is not free to recompute
    if (use_wb_static_) {
      pr -= pwb(m,k,j,i);
    }
    Real m_ii_h = pr + 0.5*rho*(v2*v2+v3*v3);
    Real m_pp = pr + rho*SQR(v3);
      
    Real src1 = z_ov_rE_(m,k,j,i) * m_ii_h;
    Real src2 = x_ov_rD_(m,k,j,i) * m_pp;
    Real src3 = -y_ov_rC_(m,k,j,i) * (uflx.x2f(m,IM3,k,j,i)*area_.x2f(m,k,j,i)+uflx.x2f(m,IM3,k,j+1,i)*area_.x2f(m,k,j+1,i))/volume_(m,k,j,i);
      
    src2 -= factor * (uflx.x1f(m,IM2,k,j,i)*area_.x1f(m,k,j,i)+uflx.x1f(m,IM2,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
    src3 -= factor * (uflx.x1f(m,IM3,k,j,i)*area_.x1f(m,k,j,i)+uflx.x1f(m,IM3,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
      
    u0(m,IM1,k,j,i) += src1*bdt;
    u0(m,IM2,k,j,i) += src2*bdt;
    u0(m,IM3,k,j,i) += src3*bdt;
    
  });
}

void Coordinates::SrcTermsSphericalPolarMHD(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc0, const DvceArray4D<Real> &pwb,
    const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
    DvceArray5D<Real> &u0) {

  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
    
  // by-value copy of the EOS, capturable in the device lambda below
  auto eos_ = eos_data;
  // The live gas pressure is read from the derived-variable cache rather than recomputed:
  // p(d,e) is a root find for a general EOS, and ConsToPrim has already evaluated it once
  // for this very state. An ideal gas allocates no wder array, and (gamma-1)*e is free.
  const bool gen_ = eos_.IsGeneral();
  auto &wder_ = pmy_pack->pmhd->wder;

  auto volume_ = volume;
  auto area_ = area;
  auto xx1f_ = xx1f;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  const bool use_wb_static_ = pmy_pack->pmhd->use_wellbalance_static;
  par_for("spsrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      
    Real r_l = xx1f_(m,i);
    Real r_r = xx1f_(m,i+1);
    Real factor = (r_r-r_l)/(r_r+r_l);
      
    Real v1 = w0(m,IVX,k,j,i);
    Real v2 = w0(m,IVY,k,j,i);
    Real v3 = w0(m,IVZ,k,j,i);
    Real rho = w0(m,IDN,k,j,i);
    Real pr = gen_ ? wder_(m,IDPR,k,j,i)
                            : eos_.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
    // subtract the background pressure, which for a general EOS depends on the background
    // DENSITY as well as its internal energy and so cannot be written as (gamma-1)*e_bg
    if (use_wb_static_) {
      pr -= pwb(m,k,j,i);
    }
    Real m_ii_h = pr + 0.5*rho*(v2*v2+v3*v3);
    Real m_pp = pr + rho*SQR(v3);
    m_ii_h += 0.5*SQR(bcc0(m,IBX,k,j,i));
    m_pp += 0.5*( SQR(bcc0(m,IBX,k,j,i)) + SQR(bcc0(m,IBY,k,j,i)) - SQR(bcc0(m,IBZ,k,j,i)) );
      
    Real src1 = z_ov_rE_(m,k,j,i) * m_ii_h;
    Real src2 = x_ov_rD_(m,k,j,i) * m_pp;
    Real src3 = -y_ov_rC_(m,k,j,i) * (uflx.x2f(m,IM3,k,j,i)*area_.x2f(m,k,j,i)+uflx.x2f(m,IM3,k,j+1,i)*area_.x2f(m,k,j+1,i))/volume_(m,k,j,i);
      
    src2 -= factor * (uflx.x1f(m,IM2,k,j,i)*area_.x1f(m,k,j,i)+uflx.x1f(m,IM2,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
    src3 -= factor * (uflx.x1f(m,IM3,k,j,i)*area_.x1f(m,k,j,i)+uflx.x1f(m,IM3,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
      
    u0(m,IM1,k,j,i) += src1*bdt;
    u0(m,IM2,k,j,i) += src2*bdt;
    u0(m,IM3,k,j,i) += src3*bdt;
    
  });
}
