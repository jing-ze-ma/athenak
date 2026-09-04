//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sp_test.cpp
//! \brief SPHERICAL-POLAR counterpart of cs_test iprob=8: a uniform gas AT REST threaded
//! by a UNIFORM Cartesian field B = b0 zhat.  curl B = 0, so the exact evolution is
//! nothing; after one step the momentum per cell is the spurious force of the geometric
//! source term against the flux divergence (the cancellation the cubed sphere's
//! well-balanced source was built for).  Measure it from two consecutive dumps.
//!
//! The face-normal components are FACE AVERAGES, not point values: on an r-face the
//! area-weighted mean of b0 cos(theta) is b0 (cos th_l + cos th_r)/2, and on a theta-face
//! -b0 sin(theta_f) is constant.  With those the discrete divergence is zero EXACTLY
//! (both terms reduce to b0 (r_r^2 - r_l^2) dphi (sin^2 th_r - sin^2 th_l)/2), so no
//! monopole force contaminates the measurement.

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (!pmy_mesh_->use_spherical_polar || pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "sp_test requires mesh/use_spherical_polar = true and an <mhd> block"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  const Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  const Real b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  // econsistent = true: the magnetic energy is built from the SAME cell-centred bcc that
  // ConsToPrim subtracts, so the recovered pressure is exactly uniform and only the
  // geometric-source/flux-divergence mismatch survives.  false (cs_test's choice): the
  // analytic b0^2/2, so the |bcc|^2 - b0^2 defect appears as an O(dtheta^2) pressure
  // ripple, largest in the polar rows -- which then dominates the one-step force.
  const bool econs = pin->GetOrAddBoolean("problem", "econsistent", true);
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0f = pmbp->pmhd->b0;
  auto &bcc = pmbp->pmhd->bcc0;
  auto &x2f_ = pmbp->pcoord->xx2f;   // theta at the theta-faces, (m, j)

  par_for("sp_test_b", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (j <= je && k <= ke) {
      b0f.x1f(m,k,j,i) = b0*0.5*(std::cos(x2f_(m,j)) + std::cos(x2f_(m,j+1)));
    }
    if (i <= ie && k <= ke) {
      b0f.x2f(m,k,j,i) = -b0*std::sin(x2f_(m,j));
    }
    if (i <= ie && j <= je) {
      b0f.x3f(m,k,j,i) = 0.0;
    }
  });

  par_for("sp_test_u", DevExeSpace(), 0,nmb1, ks,ke, js,je, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IDN,k,j,i) = d0;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    bcc(m,IBX,k,j,i) = 0.5*(b0f.x1f(m,k,j,i) + b0f.x1f(m,k,j,i+1));
    bcc(m,IBY,k,j,i) = 0.5*(b0f.x2f(m,k,j,i) + b0f.x2f(m,k,j+1,i));
    bcc(m,IBZ,k,j,i) = 0.5*(b0f.x3f(m,k,j,i) + b0f.x3f(m,k+1,j,i));
    const Real bsq = econs ? (SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                             + SQR(bcc(m,IBZ,k,j,i))) : b0*b0;
    u0(m,IEN,k,j,i) = p0/gm1 + 0.5*bsq;
  });
  return;
}
