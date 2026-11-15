//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coordinates.cpp
//! \brief Implementation of the Coordinates class: geometry of the curvilinear grids
//! (cubed sphere, spherical polar), the metric velocity conversions on the cubed sphere,
//! and the geometric source terms (GR, cubed sphere, spherical polar).

#include <algorithm>
#include <cmath>
#include <iostream> // cout
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "cartesian_ks.hpp"
#include "coordinates.hpp"
#include "cell_locations.hpp"
#include "cubed_sphere.hpp"
#include "gnomonic_raisevel.hpp"
#include "gnomonic_raisevel_mhd.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "utils/recipe.hpp"

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
    xx1f("xx1f",1,1), xx2f("xx2f",1,1), xx3f("xx3f",1,1), dxedge("dxe",1,1,1,1),
    dxface("dxf",1,1,1,1),
    sin_cell("sin_cell",1,1,1), cos_cell("cos_cell",1,1,1),
    sin_face_xi("sin_face_xi",1,1,1), cos_face_xi("cos_face_xi",1,1,1),
    sin_face_eta("sin_face_eta",1,1,1), cos_face_eta("cos_face_eta",1,1,1),
    x_ov_rD("x_ov_rD",1,1,1,1), y_ov_rC("y_ov_rC",1,1,1,1), z_ov_rE("z_ov_rC",1,1,1,1) {
  // defaults-1009 (user 10-09): sp_x2_periodic_image defaults ON on a FRESH run on a
  // theta-PERIODIC spherical-polar mesh (the only place it acts, BuildX2vRecon); a
  // restart whose file lacks the key and every other mesh keep the old default (off)
  const bool img_def = !global_variable::restart_run &&
      pmy_pack->pmesh->use_spherical_polar && pmy_pack->pmesh->multi_d &&
      pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic &&
      !pmy_pack->pmesh->use_grid_stretch_theta;
  // curvilinear source-term options, from <mhd> or <hydro> (read only from an existing
  // block, since AddPhysics builds modules from the blocks present)
  if (pin->DoesBlockExist("mhd")) {
    sp_x2_periodic_image = pin->GetOrAddBoolean("mhd","sp_x2_periodic_image",img_def);
  } else if (pin->DoesBlockExist("hydro")) {
    sp_x2_periodic_image = pin->GetOrAddBoolean("hydro","sp_x2_periodic_image",
                                                img_def);
  }

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
    Kokkos::realloc(x2v_rec, nmb, ncells2);
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
    // gnomonic trig depends on (xi, eta) only, indexed (m,k,j); face arrays are
    // staggered in j (xi faces) or k (eta faces)
    Kokkos::realloc(sin_cell, nmb, ncells3, ncells2);
    Kokkos::realloc(cos_cell, nmb, ncells3, ncells2);
    Kokkos::realloc(sin_face_xi, nmb, ncells3, ncells2+1);
    Kokkos::realloc(cos_face_xi, nmb, ncells3, ncells2+1);
    Kokkos::realloc(sin_face_eta, nmb, ncells3+1, ncells2);
    Kokkos::realloc(cos_face_eta, nmb, ncells3+1, ncells2);
    CoordGnomonicEquiangle();
  }
  if (pmy_pack->pmesh->use_spherical_polar) {
    // prod1103d default (USER 11-03): sp_phi_global = true on a FRESH spherical-polar
    // FULL-SPHERE start (polar boundaries at both poles, x3 periodic over 2 pi;
    // recipe::NewDefaults) when the key is not named: the phi faces/centres come from the
    // global index, so the IC is MeshBlock-layout exact.  Recorded for restarts; a
    // restart whose input lacks the key keeps false (bitwise); an explicit key wins.
    {
      Mesh *pm = pmy_pack->pmesh;
      const bool fullsph = pm->use_spherical_polar && pm->use_polar_boundary &&
      pin->GetOrAddString("mesh", "ix2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ox2_bc", "") == "polar" &&
      pin->GetOrAddString("mesh", "ix3_bc", "") == "periodic" &&
      std::abs(pm->mesh_size.x3max - pm->mesh_size.x3min - 2.0*M_PI) < 1.0e-8;
      if (fullsph && !pin->DoesParameterExist("mesh", "sp_phi_global")) {
        if (recipe::NewDefaults(pin)) {
          pin->GetOrAddBoolean("mesh", "sp_phi_global", true);
          recipe::Log("<mesh>/sp_phi_global = true (default on a fresh sp full-sphere start)");
        } else {
          recipe::Log("<mesh>/sp_phi_global = false (" + std::string(
              recipe::GetMode(pin) == recipe::Mode::none ? "recipe none" :
              "restart input has no key") + ": the old default)");
        }
      }
    }
    // layinv-1105: with <mesh>/layout_invariant (fresh default) both angles come from the
    // global index (theta: sp_theta_global; phi: sp_phi_global, now also on wedges)
    // unless the key is named; without layout_invariant the old rules, bitwise.
    {
      const bool li = pmy_pack->pmesh->layout_invariant;
      sp_phi_global = pin->DoesParameterExist("mesh", "sp_phi_global") ?
                      pin->GetBoolean("mesh", "sp_phi_global") : li;
      sp_theta_global = pin->DoesParameterExist("mesh", "sp_theta_global") ?
                        pin->GetBoolean("mesh", "sp_theta_global") : li;
    }
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
//! \fn Real GnomonicCornerOmega, GnomonicSolidAngle
//! \brief Exact solid angle of the gnomonic cell bounded by X = xl,xr and Y = yl,yr
//! (X = tan(xi), Y = tan(eta)), as a difference of a corner potential. Unlike the
//! midpoint form it is additive across refinement levels, so flux correction at
//! coarse/fine boundaries stays consistent with the geometric source term.


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

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::CoordGnomonicEquiangle
//! \brief Computes the geometry of the gnomonic equiangular cubed sphere (x1 radial,
//! x2 = xi, x3 = eta): coordinate arrays, face areas, volumes, edge lengths, dual-mesh
//! lengths and areas, the trig of the angle between the tangent basis vectors and the
//! geometric source coefficients.

void Coordinates::CoordGnomonicEquiangle() {
  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  (void)ie; (void)je; (void)ke;

  // x1 = radius, xi = x2, eta = x3 (see mesh.hpp). Geometry covers the ghost zones too.
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;

  // radial grid stretch (same map as on the spherical-polar grid)
  const bool str_r_ = pmy_pack->pmesh->use_grid_stretch_r;
  const bool str_rp_ = pmy_pack->pmesh->use_grid_stretch_r_poly;
  const Real fstr_r_ = pmy_pack->pmesh->fStretchR;
  const Real rmin_ = pmy_pack->pmesh->mesh_size.x1min;
  const Real rmax_ = pmy_pack->pmesh->mesh_size.x1max;
  Real cpoly_[NSTRETCH_R_POLY];
  for (int n=0; n<NSTRETCH_R_POLY; ++n) {
    cpoly_[n] = pmy_pack->pmesh->fStretchRPoly[n];
  }

  // 1-D coordinate arrays: x1 is the radius with the volume-centroid centre (as in
  // CoordSphericalPolar); x2 and x3 are the panel coordinates (not angles)
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
      if (i == n1m1) {
        xx1f_(m,i+1) = r_r;
      }
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
    // radial (x1) faces and centre
    Real r_c = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real r_l = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real r_r = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    // radial stretch
    if (str_r_) {
      StretchR(fstr_r_, rmin_, rmax_, r_c);
      StretchR(fstr_r_, rmin_, rmax_, r_l);
      StretchR(fstr_r_, rmin_, rmax_, r_r);
    } else if (str_rp_) {
      StretchRPoly(cpoly_, rmin_, rmax_, r_c);
      StretchRPoly(cpoly_, rmin_, rmax_, r_l);
      StretchRPoly(cpoly_, rmin_, rmax_, r_r);
    }
    // cell centre at the volume centroid of the shell
    r_c = RadialCentroid(r_l, r_r);

    // angles: xi on x2, eta on x3
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

    // angular face widths
    Real dth_xi = acos((1.0 + xl*xr + y*y) / sqrt(1.0 + xl*xl + y*y)
                    / sqrt(1.0 + xr*xr + y*y));
    Real dth_xil = acos((1.0 + xl*xr + yl*yl) / sqrt(1.0 + xl*xl + yl*yl)
                    / sqrt(1.0 + xr*xr + yl*yl));
    Real dth_xir = acos((1.0 + xl*xr + yr*yr) / sqrt(1.0 + xl*xl + yr*yr)
                    / sqrt(1.0 + xr*xr + yr*yr));
    Real dth_eta = acos((1.0 + x*x + yl*yr) / sqrt(1.0 + x*x + yl*yl)
                    / sqrt(1.0 + x*x + yr*yr));
    Real dth_etal = acos((1.0 + xl*xl + yl*yr) / sqrt(1.0 + xl*xl + yl*yl)
                    / sqrt(1.0 + xl*xl + yr*yr));
    Real dth_etar = acos((1.0 + xr*xr + yl*yr) / sqrt(1.0 + xr*xr + yl*yl)
                    / sqrt(1.0 + xr*xr + yr*yr));

    // radial faces (x1): exact solid angle, additive across refinement levels
    const Real domega = GnomonicSolidAngle(xl, xr, yl, yr);
    area_.x1f(m,k,j,i) = SQR(r_l) * domega;
    Real area1r = SQR(r_r) * domega;
    if (i == n1m1) area_.x1f(m,k,j,i+1) = area1r;

    // xi faces (x2)
    area_.x2f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etal;
    Real area2r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_etar;
    if (j == n2m1) area_.x2f(m,k,j+1,i) = area2r;

    // eta faces (x3)
    area_.x3f(m,k,j,i) = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xil;
    Real area3r = 0.5 * (SQR(r_r)-SQR(r_l)) * dth_xir;
    if (k == n3m1) area_.x3f(m,k+1,j,i) = area3r;

    // volume (exact solid angle x radial shell)
    volume_(m,k,j,i) = 1.0/3.0*(r_r*r_r*r_r-r_l*r_l*r_l) * domega;
    dx1_(m,k,j,i) = r_r - r_l;
    dx2_(m,k,j,i) = r_c * dth_xi;
    dx3_(m,k,j,i) = r_c * dth_eta;

    // edge lengths for the CT update, from the same arcs as the adjoining face areas:
    // x1e radial, x2e along xi at an eta face, x3e along eta at a xi face
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

    // dual mesh: edge-centred areas (areaedge) and centre-to-centre lengths (dxface) of
    // the loop joining the four cell centres around an edge, used by the current density.
    // The x1 loop lies in the tangent surface and carries the sin of the angle between
    // the basis vectors; the x2 and x3 loops contain rhat and do not.
    Real r_cm = CellCenterX(i-1-is, indcs.nx1, size.d_view(m).x1min,
                                            size.d_view(m).x1max);
    // radial centre-to-centre lengths use the stretched centroids
    Real r_lm = LeftEdgeX(i-1-is, indcs.nx1, size.d_view(m).x1min,
                                             size.d_view(m).x1max);
    if (str_r_) {
      StretchR(fstr_r_, rmin_, rmax_, r_lm);
    } else if (str_rp_) {
      StretchRPoly(cpoly_, rmin_, rmax_, r_lm);
    }
    r_cm = RadialCentroid(r_lm, r_l);   // centroid of cell i-1
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

    // geometric source coefficients for SrcTermsGnomonicEquiangle (z_ov_rE is radial)
    x_ov_rD_(m,k,j,i) = (area2r * sin_face_xir
                        - area_.x2f(m,k,j,i) * sin_face_xi_(m,k,j)) / volume_(m,k,j,i);
    y_ov_rC_(m,k,j,i) = (area3r * sin_face_etar
                        - area_.x3f(m,k,j,i) * sin_face_eta_(m,k,j)) / volume_(m,k,j,i);
    z_ov_rE_(m,k,j,i) = (area1r - area_.x1f(m,k,j,i)) / volume_(m,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::GnomonicEquiangleRaiseVel
//! \brief Recomputes the contravariant primitive velocity and the internal energy from
//! the covariant conserved momentum on the cubed sphere, using the metric
//! g = [[1, c],[c, 1]] of the non-orthogonal gnomonic tangent basis (c = cos_cell).
//! Runs after ConsToPrim and overwrites IVX/IVY/IVZ and IEN; the energy floors and the
//! velocity ceiling are re-applied to the corrected internal energy.
//! Non-relativistic ideal/general hydro only.

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

  // isothermal EOS: no energy, only the velocity is raised

  if (!eos_data.is_ideal) {
    par_for("cs_raisev_iso", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = cos_cell_(m,k,j);
      const Real det = 1.0 - c*c;
      const Real d = u0(m,IDN,k,j,i);
      const Real m2 = u0(m,IM2,k,j,i);   // xi
      const Real m3 = u0(m,IM3,k,j,i);   // eta
      w0(m,IVX,k,j,i) = u0(m,IM1,k,j,i)/d;
      w0(m,IVY,k,j,i) = (m2 - c*m3)/(d*det);
      w0(m,IVZ,k,j,i) = (m3 - c*m2)/(d*det);
    });
    return;
  }

  // default kernel when no floor switch is enabled (EOS_Data::floors_legacy): kept
  // separate to preserve its floating-point results
  if (eos_data.floors_legacy) {
    par_for("cs_raisev", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = cos_cell_(m,k,j);
      const Real det = 1.0 - c*c;
      const Real d = u0(m,IDN,k,j,i);
      const Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
      const Real m2 = u0(m,IM2,k,j,i);   // xi
      const Real m3 = u0(m,IM3,k,j,i);   // eta
      // contravariant velocity v^i = g^{ij} m_j / rho
      const Real v1 = m1/d;
      const Real v2 = (m2 - c*m3)/(d*det);
      const Real v3 = (m3 - c*m2)/(d*det);
      w0(m,IVX,k,j,i) = v1;
      w0(m,IVY,k,j,i) = v2;
      w0(m,IVZ,k,j,i) = v3;
      // kinetic energy 0.5 m_i v^i
      Real eint = u0(m,IEN,k,j,i) - 0.5*(m1*v1 + m2*v2 + m3*v3);
      // re-apply the floors to the internal energy corrected by the metric kinetic
      // energy (as in SingleC2P_GeneralHyd), updating u0 consistently
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
        if (eint != eold) {
          u0(m,IEN,k,j,i) = eint + 0.5*(m1*v1+m2*v2+m3*v3);
        }
      }
      w0(m,IEN,k,j,i) = eint;
    });
    return;
  }

  // density-floor corrections and velocity ceiling that need the metric kinetic energy
  const bool keepv_ = eos_data.dfloor_keep_velocity && eos_data.defer_cons_floors;
  const bool keept_ = eos_data.dfloor_keep_temperature && eos_data.defer_cons_floors;
  const Real vceil_ = eos_data.vceil;
  // the ceiling acts in active cells only (ghosts are overwritten by the next exchange)
  auto &aidx = pmy_pack->pmesh->mb_indcs;
  const int ais = aidx.is, aie = aidx.ie;
  const int ajs = aidx.js, aje = aidx.je;
  const int aks = aidx.ks, ake = aidx.ke;
  auto dflfv_ = (pmy_pack->phydro != nullptr) ? pmy_pack->phydro->dfl_fv
                                              : DvceArray4D<Real>("cs_fv_dm",1,1,1,1);

  // energy created by the floors (applied here on the cubed sphere, see
  // EOS_Data::defer_cons_floors)
  Real efloor_de_ = 0.0;
  // energy removed by the velocity ceiling (EventCounters::vceil_de)
  Real vceil_de_ = 0.0;
  int nceilv_ = 0;
  // count of cells where efloor_as_tfloor rebuilt e from T
  int ntset_ = 0;
  const int nkji = (ku - kl + 1)*(ju - jl + 1)*(iu - il + 1);
  const int nji = (ju - jl + 1)*(iu - il + 1);
  const int ni = (iu - il + 1);
  // HIP: LaunchBounds<256, 2> avoids register spills of the default 1024-thread bound
#if defined(KOKKOS_ENABLE_HIP)
  using RvPol = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 2>>;
#else
  using RvPol = Kokkos::RangePolicy<>;
#endif
  Kokkos::parallel_reduce("cs_raisev",
  RvPol(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int &idx, int &sumv, int &sumts, Real &sumde, Real &sumvde) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + kl;
    const int j = (idx - m*nkji - (k - kl)*nji)/ni + jl;
    const int i = (idx - m*nkji - (k - kl)*nji - (j - jl)*ni) + il;
    // per-cell metric inversion, floors and ceiling (gnomonic_raisevel.hpp)
    const Real c = cos_cell_(m,k,j);
    const Real d = u0(m,IDN,k,j,i);
    Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    Real m2 = u0(m,IM2,k,j,i);   // xi
    Real m3 = u0(m,IM3,k,j,i);   // eta
    Real etot = u0(m,IEN,k,j,i);
    // ceiling only in active cells
    const bool act_ = (i >= ais && i <= aie && j >= ajs && j <= aje &&
                       k >= aks && k <= ake);
    Real v1, v2, v3, eint, pnew, g1new, temp, de, vde;
    bool ceil_used, floored, tset_used;
    GnomonicRaiseVelFloors(c, eos_, gen_, keepv_, keept_, vceil_, act_,
                           (keepv_ || keept_) ? dflfv_(m,k,j,i) : 1.0,
                           gen_ ? wtemp_(m,k,j,i) : 0.0, d,
                           m1, m2, m3, etot, v1, v2, v3, eint, pnew, g1new, temp,
                           ceil_used, floored, tset_used, de, vde);
    u0(m,IM1,k,j,i) = m1;
    u0(m,IM2,k,j,i) = m2;
    u0(m,IM3,k,j,i) = m3;
    u0(m,IEN,k,j,i) = etot;
    if (gen_) {
      wder_(m,IDPR,k,j,i) = pnew;
      wder_(m,IDG1,k,j,i) = g1new;
      wtemp_(m,k,j,i) = temp;
    }
    if (ceil_used) {
      sumv++;
    }
    if (tset_used) {
      sumts++;
    }
    sumde += de;
    sumvde += vde;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    w0(m,IEN,k,j,i) = eint;
  }, Kokkos::Sum<int>(nceilv_), Kokkos::Sum<int>(ntset_),
     Kokkos::Sum<Real>(efloor_de_), Kokkos::Sum<Real>(vceil_de_));
  pmy_pack->pmesh->ecounter.efloor_de += efloor_de_;
  pmy_pack->pmesh->ecounter.vceil_de += vceil_de_;
  pmy_pack->pmesh->ecounter.neos_vceil += nceilv_;
  pmy_pack->pmesh->ecounter.neos_tset += ntset_;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::GnomonicEquiangleRaiseVelMHD
//! \brief MHD version of GnomonicEquiangleRaiseVel. In addition, rebuilds the
//! cell-centred field from the faces in the orthonormal frame
//! {rhat, e_xi, (e_eta - c e_xi)/s}, since the face-normal triple
//! (B.rhat, B.nhat_xi, B.nhat_eta) held in b0 is not orthogonal, and subtracts the corresponding magnetic energy.

void Coordinates::GnomonicEquiangleRaiseVelMHD(DvceArray5D<Real> &u0,
    const DvceFaceFld4D<Real> &b0, DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
    const EOS_Data &eos_data, const DvceArray5D<Real> &wder,
    const DvceArray4D<Real> &wtemp, const int il, const int iu, const int jl,
    const int ju, const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  auto &sin_cell_ = sin_cell;
  // a general EOS re-solves p, Gamma_1 and T from the corrected internal energy, which
  // also uses the magnetic energy of the orthonormal cell-centred field
  const bool gen_ = eos_data.IsGeneral();
  auto eos_ = eos_data;
  auto wder_ = wder;
  auto wtemp_ = wtemp;
  // isothermal EOS (no energy): only the orthonormal field and the raised velocity
  if (!eos_data.is_ideal) {
    auto &x1v_i = x1v;
    auto &x1f_i = xx1f;
    par_for("cs_raisev_mhd_iso", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = cos_cell_(m,k,j);
      const Real sn = sin_cell_(m,k,j);
      const Real det = 1.0 - c*c;
      const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      bcc0(m,IBX,k,j,i) = CellCenteredRadialFld(b0.x1f(m,k,j,i), b0.x1f(m,k,j,i+1),
                                                x1f_i(m,i), x1f_i(m,i+1), x1v_i(m,i));
      bcc0(m,IBY,k,j,i) = (by_n + c*bz_n)/sn;
      bcc0(m,IBZ,k,j,i) = bz_n;
      const Real d = u0(m,IDN,k,j,i);
      const Real m2 = u0(m,IM2,k,j,i);   // xi
      const Real m3 = u0(m,IM3,k,j,i);   // eta
      w0(m,IVX,k,j,i) = u0(m,IM1,k,j,i)/d;
      w0(m,IVY,k,j,i) = (m2 - c*m3)/(d*det);
      w0(m,IVZ,k,j,i) = (m3 - c*m2)/(d*det);
    });
    return;
  }
  // radial field at the volume centroid by weighted face interpolation; ghosts behind a
  // reflecting radial wall use mirror-consistent weights (CellCenteredRadialFldWall)
  const bool str_x1_ = true;
  auto &wbc_ = pmy_pack->pmb->mb_bcs;
  const int wis_ = pmy_pack->pmesh->mb_indcs.is, wie_ = pmy_pack->pmesh->mb_indcs.ie;
  auto &x1v_ = x1v;
  auto &x1f_ = xx1f;

  // velocity ceiling on the metric |v| (applied here on the cubed sphere)
  const Real vceil_ = eos_data.vceil;
  // the ceiling acts in active cells only
  auto &aidx = pmy_pack->pmesh->mb_indcs;
  const int ais = aidx.is, aie = aidx.ie;
  const int ajs = aidx.js, aje = aidx.je;
  const int aks = aidx.ks, ake = aidx.ke;

  int nceilv_ = 0;
  const int nkji = (ku - kl + 1)*(ju - jl + 1)*(iu - il + 1);
  const int nji = (ju - jl + 1)*(iu - il + 1);
  const int ni = (iu - il + 1);
  Kokkos::parallel_reduce("cs_raisev_mhd",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int &idx, int &sumv) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + kl;
    const int j = (idx - m*nkji - (k - kl)*nji)/ni + jl;
    const int i = (idx - m*nkji - (k - kl)*nji - (j - jl)*ni) + il;
    const Real c = cos_cell_(m,k,j);
    const Real sn = sin_cell_(m,k,j);

    // cell-centred field in the orthonormal frame {rhat, e_xi, (e_eta - c e_xi)/s}: only
    // the xi slot changes, B.e_xi = (b_xi + c b_eta)/s. x2/x3 are plain face averages
    // (the panel coordinates are uniform).
    const bool rin_ = wbc_.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect;
    const bool rout_ = wbc_.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect;
    const Real bx = str_x1_ ?
        CellCenteredRadialFldWall(b0.x1f, x1f_, x1v_, m, k, j, i, wis_, wie_,
                                  rin_, rout_) :
        0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
    const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    const Real by = (by_n + c*bz_n)/sn;
    const Real bz = bz_n;
    bcc0(m,IBX,k,j,i) = bx;
    bcc0(m,IBY,k,j,i) = by;
    bcc0(m,IBZ,k,j,i) = bz;

    // per-cell metric inversion, floors and ceiling (gnomonic_raisevel_mhd.hpp)
    const Real d = u0(m,IDN,k,j,i);
    Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    Real m2 = u0(m,IM2,k,j,i);   // xi
    Real m3 = u0(m,IM3,k,j,i);   // eta
    Real etot = u0(m,IEN,k,j,i);
    const bool act_ = (i >= ais && i <= aie && j >= ajs && j <= aje &&
                       k >= aks && k <= ake);
    Real v1, v2, v3, eint, pnew, g1new, temp;
    bool ceil_used, floored;
    GnomonicRaiseVelMHDFloors(c, eos_, gen_, vceil_, act_,
                              gen_ ? wtemp_(m,k,j,i) : 0.0, d, bx, by, bz,
                              m1, m2, m3, etot, v1, v2, v3, eint, pnew, g1new, temp,
                              ceil_used, floored);
    if (ceil_used) {
      u0(m,IM1,k,j,i) = m1;
      u0(m,IM2,k,j,i) = m2;
      u0(m,IM3,k,j,i) = m3;
      sumv++;
    }
    if (ceil_used || floored) {
      u0(m,IEN,k,j,i) = etot;
    }
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    if (gen_) {
      wder_(m,IDPR,k,j,i) = pnew;
      wder_(m,IDG1,k,j,i) = g1new;
      wtemp_(m,k,j,i) = temp;
    }
    w0(m,IEN,k,j,i) = eint;
  }, Kokkos::Sum<int>(nceilv_));
  pmy_pack->pmesh->ecounter.neos_vceil += nceilv_;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::GnomonicCellCenteredB
//! \brief Cell-centred field of GnomonicEquiangleRaiseVelMHD only (same expressions), for
//! the intermediate RKG super-stages (<mhd>/rkg_lean), which change only B and the total
//! energy.


void Coordinates::GnomonicCellCenteredB(const DvceFaceFld4D<Real> &b0,
    DvceArray5D<Real> &bcc0, const int il, const int iu, const int jl, const int ju,
    const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  auto &sin_cell_ = sin_cell;
  // ghosts behind a reflecting radial wall use mirror-consistent radial weights
  auto &wbc_ = pmy_pack->pmb->mb_bcs;
  const int wis_ = pmy_pack->pmesh->mb_indcs.is, wie_ = pmy_pack->pmesh->mb_indcs.ie;
  auto &x1v_ = x1v;
  auto &x1f_ = xx1f;
  par_for("cs_bcc_only", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real sn = sin_cell_(m,k,j);
    const bool rin_ = wbc_.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect;
    const bool rout_ = wbc_.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect;
    const Real bx = CellCenteredRadialFldWall(b0.x1f, x1f_, x1v_, m, k, j, i, wis_, wie_,
                                              rin_, rout_);
    const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    bcc0(m,IBX,k,j,i) = bx;
    bcc0(m,IBY,k,j,i) = (by_n + c*bz_n)/sn;
    bcc0(m,IBZ,k,j,i) = bz_n;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::GnomonicEquiangleLowerMom
//! \brief Computes the conserved state from primitives on the cubed sphere (inverse of
//! GnomonicEquiangleRaiseVel), lowering the velocity index with the metric.

void Coordinates::GnomonicEquiangleLowerMom(const DvceArray5D<Real> &w0,
    DvceArray5D<Real> &u0, const int il, const int iu, const int jl, const int ju,
    const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;

  // isothermal EOS: only density and momentum
  bool iso_ = false;
  if (pmy_pack->phydro != nullptr) {
    iso_ = !pmy_pack->phydro->peos->eos_data.is_ideal;
  } else if (pmy_pack->pmhd != nullptr) {
    iso_ = !pmy_pack->pmhd->peos->eos_data.is_ideal;
  }
  if (iso_) {
    par_for("cs_lowerm_iso", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = cos_cell_(m,k,j);
      const Real d = w0(m,IDN,k,j,i);
      const Real v2 = w0(m,IVY,k,j,i);
      const Real v3 = w0(m,IVZ,k,j,i);
      u0(m,IDN,k,j,i) = d;
      u0(m,IM1,k,j,i) = d*w0(m,IVX,k,j,i);
      u0(m,IM2,k,j,i) = d*(v2 + c*v3);
      u0(m,IM3,k,j,i) = d*(v3 + c*v2);
    });
    return;
  }

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

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::SrcTermsGnomonicEquiangle
//! \brief Hydro geometric source terms on the cubed sphere.
void Coordinates::SrcTermsGnomonicEquiangle(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &wder, const DvceArray4D<Real> &pwb,
    const DvceFaceFld5D<Real> uflx,
    const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {
  SrcTermsGnomonicEquiangleImpl(w0, w0, false, wder, pwb,
                                pmy_pack->phydro->use_wellbalance_static,
                                uflx, eos_data, bdt, u0);
  return;
}
//----------------------------------------------------------------------------------------
//! \fn void Coordinates::SrcTermsGnomonicEquiangleMHD
//! \brief MHD geometric source terms on the cubed sphere.
void Coordinates::SrcTermsGnomonicEquiangleMHD(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc0, const DvceArray5D<Real> &wder,
    const DvceArray4D<Real> &pwb,
    const DvceFaceFld5D<Real> uflx, const EOS_Data &eos_data, const Real bdt,
    DvceArray5D<Real> &u0) {
  SrcTermsGnomonicEquiangleImpl(w0, bcc0, true, wder, pwb,
                                pmy_pack->pmhd->use_wellbalance_static,
                                uflx, eos_data, bdt, u0);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Coordinates::SrcTermsGnomonicEquiangleImpl
//! \brief Geometric source terms of the momentum equation on the cubed sphere, shared by
//! hydro and MHD (is_mhd adds the Maxwell stress).
void Coordinates::SrcTermsGnomonicEquiangleImpl(const DvceArray5D<Real> &w0,
    const DvceArray5D<Real> &bcc0, const bool is_mhd,
    const DvceArray5D<Real> &wder, const DvceArray4D<Real> &pwb,
    const bool use_wb_static, const DvceFaceFld5D<Real> uflx,
    const EOS_Data &eos_data, const Real bdt, DvceArray5D<Real> &u0) {
  auto &size = pmy_pack->pmb->mb_size;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  // by-value copy of the EOS, capturable in the device lambda below
  auto eos_ = eos_data;
  // gas pressure from the derived-variable cache (wder) for a general EOS, (gamma-1)*e
  // for an ideal gas
  const bool gen_ = eos_.IsGeneral();
  auto &wder_ = wder;
  auto &bcc_ = bcc0;
  const bool mhd_ = is_mhd;
  auto pwb_ = pwb;
  const bool wbs_ = use_wb_static;

  auto volume_ = volume;
  auto area_ = area;
  auto sin_cell_ = sin_cell;
  auto cos_cell_ = cos_cell;
  auto x_ov_rD_ = x_ov_rD;
  auto y_ov_rC_ = y_ov_rC;
  auto z_ov_rE_ = z_ov_rE;
  auto xx1f_ = xx1f;
  par_for("cssrc", DevExeSpace(), 0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // radial factors from the stretched faces
    const Real r_l = xx1f_(m,i);
    const Real r_r = xx1f_(m,i+1);
    const Real factor = (r_r-r_l)/(r_r+r_l);

    Real v1 = w0(m,IVX,k,j,i);   // radial
    Real v2 = w0(m,IVY,k,j,i);   // xi
    Real v3 = w0(m,IVZ,k,j,i);   // eta
    // isothermal EOS: no energy in w0 (IEN is out of range), p = d c_s^2
    Real pr = (!eos_.is_ideal) ? w0(m,IDN,k,j,i)*SQR(eos_.iso_cs)
              : gen_ ? wder_(m,IDPR,k,j,i)
              : eos_.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
    // subtract the static well-balanced background pressure (also removed from the
    // fluxes by RemoveWbFlux)
    if (wbs_) pr -= pwb_(m,k,j,i);
    Real rho = w0(m,IDN,k,j,i);

    Real sine = sin_cell_(m,k,j);
    Real cosine = cos_cell_(m,k,j);
    Real sine2 = SQR(sine);
    // lower the index of the angular components with g = [[1,c],[c,1]]
    Real v_2 = v2 + v3 * cosine;
    Real v_3 = v3 + v2 * cosine;

    // MHD: T^ij = rho v^i v^j + (p + B^2/2) g^ij - B^i B^j, so p becomes the total
    // pressure and the Reynolds stress becomes Reynolds minus Maxwell
    Real ptot = pr;
    Real b_tsq = 0.0, bb_2 = 0.0, bb_3 = 0.0;
    if (mhd_) {
      // bcc is the orthonormal triple (B.rhat, B.e_xi, B.(e_eta - c e_xi)/s), so |B|^2 is
      // the sum of squares; B^eta = b3/s and s*B^xi = s*b2 - c*b3
      const Real b1 = bcc_(m,IBX,k,j,i);
      const Real b2 = bcc_(m,IBY,k,j,i);
      const Real b3 = bcc_(m,IBZ,k,j,i);
      b_tsq = b2*b2 + b3*b3;
      bb_2 = SQR(sine*b2 - cosine*b3);
      bb_3 = b3*b3;
      ptot += 0.5*(b1*b1 + b_tsq);
    }

    // radial source z_ov_rE*(p + rho|v_t|^2/2), for MHD with p -> p + B_r^2/2
    Real src1 = z_ov_rE_(m,k,j,i) * (ptot + 0.5*(rho*(v2*v_2+v3*v_3) - b_tsq));
    Real src2 = x_ov_rD_(m,k,j,i) * (ptot + rho*SQR(v3)*sine2 - bb_3);
    Real src3 = y_ov_rC_(m,k,j,i) * (ptot + rho*SQR(v2)*sine2 - bb_2);

    // radial-face flux of each angular momentum, converted into a source by the
    // r-dependence of the angular basis
    src2 -= factor * (uflx.x1f(m,IM2,k,j,i)*area_.x1f(m,k,j,i)
             + uflx.x1f(m,IM2,k,j,i+1)*area_.x1f(m,k,j,i+1))/volume_(m,k,j,i);
    src3 -= factor * (uflx.x1f(m,IM3,k,j,i)*area_.x1f(m,k,j,i)
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
  // accpolar-1103 <mesh>/sp_phi_global: the global phi index offset of every block
  const bool phg = sp_phi_global;
  const int nx3t = pmy_pack->pmesh->mesh_indcs.nx3;
  const Real mx3min = pmy_pack->pmesh->mesh_size.x3min;
  const Real mx3max = pmy_pack->pmesh->mesh_size.x3max;
  DvceArray1D<int> k0g("sp_phi_k0", nmb1 + 1);
  {
    auto h = Kokkos::create_mirror_view(k0g);
    for (int m = 0; m <= nmb1; ++m) {
      h(m) = static_cast<int>(pmy_pack->pmesh->lloc_eachmb[pmy_pack->gids + m].lx3)*
             indcs.nx3;
    }
    Kokkos::deep_copy(k0g, h);
  }
  // layinv-1105 <mesh>/sp_theta_global (and phi with refinement): the global index of the
  // block's first cell at its own level and the global cell count at that level
  // (root-level count << (level - root level)); a refinement-free mesh has level = root,
  // i.e. exactly the cell faces/centres a one-block run computes.  The theta stretch is a
  // function of the uniform coordinate, applied afterwards exactly as before.
  const bool thg = sp_theta_global;
  const int nx2t = pmy_pack->pmesh->mesh_indcs.nx2;
  const Real mx2min = pmy_pack->pmesh->mesh_size.x2min;
  const Real mx2max = pmy_pack->pmesh->mesh_size.x2max;
  DvceArray2D<int> lvg("sp_glob_lev", nmb1 + 1, 4);   // j0, nx2 at level, k0, nx3 at level
  {
    auto h = Kokkos::create_mirror_view(lvg);
    const int rl = pmy_pack->pmesh->root_level;
    for (int m = 0; m <= nmb1; ++m) {
      const auto &ll = pmy_pack->pmesh->lloc_eachmb[pmy_pack->gids + m];
      const int sh = ll.level - rl;
      h(m,0) = static_cast<int>(ll.lx2)*indcs.nx2;
      h(m,1) = nx2t << sh;
      h(m,2) = static_cast<int>(ll.lx3)*indcs.nx3;
      h(m,3) = nx3t << sh;
    }
    Kokkos::deep_copy(lvg, h);
  }
  Real tplat[8];   // fullsph-1027 theta_plateau parameters (tplat[2] = 0: sinh map)
  for (int n=0; n<8; ++n) tplat[n] = pmy_pack->pmesh->fStretchThetaPlat[n];
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
    if (thg) {   // layinv-1105: from the global theta index (layout independent)
      const int jg = lvg(m,0) + j - js, n2g = lvg(m,1);
      theta  = CellCenterX(jg, n2g, mx2min, mx2max);
      thetal = LeftEdgeX(jg, n2g, mx2min, mx2max);
      thetar = LeftEdgeX(jg + 1, n2g, mx2min, mx2max);
    }
    // stretch the radial grid
    if (use_stretch_theta) {
      Real factor_stretch = f_stretch_theta;
      ApplyThetaStretch(factor_stretch,tplat,theta);
      ApplyThetaStretch(factor_stretch,tplat,thetal);
      ApplyThetaStretch(factor_stretch,tplat,thetar);
    }
    Real phi    = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phil   = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real phir   = LeftEdgeX(k+1-ks, indcs.nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    if (phg) {   // accpolar-1103: from the global phi index (layout independent)
      // layinv-1105: the cell count at the block's level (= nx3t on a uniform mesh)
      const int kg = k0g(m) + k - ks, n3g = lvg(m,3);
      phi  = CellCenterX(kg, n3g, mx3min, mx3max);
      phil = LeftEdgeX(kg, n3g, mx3min, mx3max);
      phir = LeftEdgeX(kg + 1, n3g, mx3min, mx3max);
    }

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
    if (j > 0 && k > 0) {
      // The dual face around an r-edge spans theta from the centroid of cell j-1 to
      // that of cell j.  At a POLE the two centroids are mirror images through the
      // axis (the ghost centroid is the negative of the active one), so the plain
      // |cos a - cos b| is IDENTICALLY ZERO there and the resistive current
      // J_r = circulation/area at the polar edge divided by zero (NaN in the first
      // cycle of any resistive run with a tangential field at the pole).  The dual
      // face straddling the pole is two caps-slices, r^2 [(1-|cos a|)+(1-|cos b|)] dphi.
      const Real ca = cos(x2v_(m,j)), cb = cos(x2v_(m,j-1));
      const bool pole =
          (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) ||
          (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1);
      const Real dcos = pole ? (2.0 - fabs(ca) - fabs(cb)) : fabs(ca - cb);
      areaedge_.x1e(m,k,j,i) = SQR(x1v_(m,i)) * dcos * (x3v_(m,k)-x3v_(m,k-1));
    }
    if (k > 0 && i > 0) areaedge_.x2e(m,k,j,i) = 0.5 * (SQR(x1v_(m,i))-SQR(x1v_(m,i-1))) * fabs(sin(x2v_(m,j))) * (x3v_(m,k)-x3v_(m,k-1));
    if (i > 0 && j > 0) areaedge_.x3e(m,k,j,i) = 0.5 * (SQR(x1v_(m,i))-SQR(x1v_(m,i-1))) * (x2v_(m,j)-x2v_(m,j-1));
  });
  BuildX2vRecon();
}

//----------------------------------------------------------------------------------------
//! \fn Coordinates::BuildX2vRecon
//! \brief theta positions used by the angular x2 reconstruction (GridPiecewiseLinearX2).
//! Equal to x2v except, with <hydro|mhd>/sp_x2_periodic_image on a theta-PERIODIC mesh,
//! in the ghost rows of blocks at x2min/x2max: there they are the PERIODIC IMAGES of the
//! cells the ghosts hold (centroid of [theta_l + L, theta_r + L] minus L, L = x2max -
//! x2min).  x2v is the sin-weighted centroid, which is not invariant under theta ->
//! theta + L, so with the plain x2v the faces x2min and x2max reconstruct the SAME cell
//! values with different stencil spacings, F2(x2min) != F2(x2max), and mass, momentum
//! and energy leak through the periodic seam (measured +2.6e23 g/cycle on the AG Car A
//! 480x8x8 wedge, mass_1009).  With the images both faces see bitwise-equivalent
//! stencils up to round-off.  Default off: x2v_rec = x2v, results bitwise unchanged.

void Coordinates::BuildX2vRecon() {
  Kokkos::deep_copy(x2v_rec, x2v);
  Mesh *pm = pmy_pack->pmesh;
  if (!sp_x2_periodic_image || !pm->use_spherical_polar || !pm->multi_d) return;
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) return;
  if (pm->use_grid_stretch_theta) return;  // images are not translations then
  auto &indcs = pm->mb_indcs;
  const int js = indcs.js, je = indcs.je, ng = indcs.ng;
  const Real x2min = pm->mesh_size.x2min, x2max = pm->mesh_size.x2max;
  const Real len = x2max - x2min;
  auto h_x2v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x2v_rec);
  auto h_x2f = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), xx2f);
  auto &size = pmy_pack->pmb->mb_size;
  auto centroid = [](Real tl, Real tr) {
    return -((sin(tr)-tr*cos(tr)) - (sin(tl)-tl*cos(tl))) / (cos(tr)-cos(tl));
  };
  const Real tol = 1.0e-12*len;
  for (int m=0; m<pmy_pack->nmb_thispack; ++m) {
    if (fabs(size.h_view(m).x2min - x2min) < tol) {
      for (int j=js-ng; j<js; ++j) {
        h_x2v(m,j) = centroid(h_x2f(m,j) + len, h_x2f(m,j+1) + len) - len;
      }
    }
    if (fabs(size.h_view(m).x2max - x2max) < tol) {
      for (int j=je+1; j<=je+ng; ++j) {
        h_x2v(m,j) = centroid(h_x2f(m,j) - len, h_x2f(m,j+1) - len) + len;
      }
    }
  }
  Kokkos::deep_copy(x2v_rec, h_x2v);
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
    // isothermal EOS: no energy in w0 (IEN is out of range), p = d c_s^2
    Real pr = (!eos_.is_ideal) ? w0(m,IDN,k,j,i)*SQR(eos_.iso_cs)
              : gen_ ? wder_(m,IDPR,k,j,i)
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
    // isothermal EOS: no energy in w0 (IEN is out of range), p = d c_s^2
    Real pr = (!eos_.is_ideal) ? w0(m,IDN,k,j,i)*SQR(eos_.iso_cs)
              : gen_ ? wder_(m,IDPR,k,j,i)
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
