//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file kh.cpp
//  \brief Problem generator for Ohmic diffusion

#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "coordinates/adm.hpp"
#include "pgen.hpp"

#include <Kokkos_Random.hpp>

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Problem Generator for Ohmic diffusion tests

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  // read problem parameters from input file

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &panel = pmbp->pmb->mb_panel;
    
    Real dens = 1.0;
    Real pres = 1.0;
    Real bx0 = 1.0;
    Real sigma = 0.1;

  // Select either Hydro or MHD
  Real gm1;
  int nfluid, nscalars;
  if (pmbp->phydro != nullptr) {
    gm1 = (pmbp->phydro->peos->eos_data.gamma) - 1.0;
    nfluid = pmbp->phydro->nhydro;
    nscalars = pmbp->phydro->nscalars;
  } else if (pmbp->pmhd != nullptr) {
    gm1 = (pmbp->pmhd->peos->eos_data.gamma) - 1.0;
    nfluid = pmbp->pmhd->nmhd;
    nscalars = pmbp->pmhd->nscalars;
  }
  if (pmbp->padm != nullptr) {
    gm1 = 1.0;
  }
  auto &w0_ = (pmbp->phydro != nullptr)? pmbp->phydro->w0 : pmbp->pmhd->w0;

  // initialize primitive variables
  Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
  par_for("pgen_ohmic", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    // set primitives in both newtonian and SR hydro
    w0_(m,IDN,k,j,i) = dens;
    w0_(m,IEN,k,j,i) = pres/gm1;
    w0_(m,IVX,k,j,i) = 0.0;
    w0_(m,IVY,k,j,i) = 0.0;
    w0_(m,IVZ,k,j,i) = 0.0;
  });

  // initialize magnetic fields if MHD
  if (pmbp->pmhd != nullptr) {
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real x1fl = LeftEdgeX(i-is, nx1, x1min, x1max);
      Real x1fr = LeftEdgeX(i+1-is, nx1, x1min, x1max);
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = bx0*exp(-SQR(x1v/sigma)/2.0);
      b0.x3f(m,k,j,i) = 0.0;
      if (i==ie) b0.x1f(m,k,j,i+1) = 0.0;
      if (j==je) b0.x2f(m,k,j+1,i) = bx0*exp(-SQR(x1v/sigma)/2.0);
      if (k==ke) b0.x3f(m,k+1,j,i) = 0.0;
      bcc0(m,IBX,k,j,i) = 0.0;
      bcc0(m,IBY,k,j,i) = bx0*exp(-SQR(x1v/sigma)/2.0);
      bcc0(m,IBZ,k,j,i) = 0.0;
    });
  }

  // Convert primitives to conserved
  if (pmbp->padm == nullptr) {
    if (pmbp->phydro != nullptr) {
      auto &u0_ = pmbp->phydro->u0;
      pmbp->phydro->peos->PrimToCons(w0_, u0_, is, ie, js, je, ks, ke);
    } else if (pmbp->pmhd != nullptr) {
      auto &u0_ = pmbp->pmhd->u0;
      auto &bcc0_ = pmbp->pmhd->bcc0;
      pmbp->pmhd->peos->PrimToCons(w0_, bcc0_, u0_, is, ie, js, je, ks, ke);
    }
  }

  return;
}
