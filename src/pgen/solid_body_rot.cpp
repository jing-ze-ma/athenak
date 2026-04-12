//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file solid_body_rot.cpp
//  \brief Problem generator for solid body rotation

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

void UserBoundary(Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Problem Generator for solid body rotation tests

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  // read problem parameters from input file
  int iprob  = pin->GetReal("problem","iprob");
    
    user_bcs_func = UserBoundary;

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &panel = pmbp->pmb->mb_panel;

  // Select either Hydro or MHD
  Real gm1, p0;
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

  bool is_relativistic = false;
  if (pmbp->pcoord->is_special_relativistic ||
      pmbp->pcoord->is_general_relativistic ||
      pmbp->pcoord->is_dynamical_relativistic) {
    is_relativistic = true;
  }

  // initialize primitive variables
  par_for("pgen_cs", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      
      if (pmy_mesh_->use_cubed_sphere) {
      Real xi = M_PI/4.0*x1v;
      Real eta = M_PI/4.0*x2v;
      Real omega = 1.0;
      
      Real x = tan(xi);
      Real y = tan(eta);
      Real delta = sqrt(1.0 + x*x + y*y);
      Real C = sqrt(1.0 + x*x);
      Real D = sqrt(1.0 + y*y);
      
      Real theta = xi+M_PI/2.0;
      Real phi = eta;
      Real xx = x3v*sin(theta)*cos(phi);
      Real yy = x3v*sin(theta)*sin(phi);
      Real vx = -omega*y*x3v;
      Real vy = omega*x*x3v;
      Real vz = 0.0;
      
      
    Real dens = 1.0;
    Real pres = 1.0;

    // set primitives in both newtonian and SR hydro
    w0_(m,IDN,k,j,i) = dens;
    w0_(m,IEN,k,j,i) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
      w0_(m,IVX,k,j,i) = ((1.0+y*y)*vx-x*y*vy)/delta/D;//*C*C/delta/delta/delta;
      w0_(m,IVY,k,j,i) = (-x*y*vx+(1.0+x*x)*vy)/delta/C;//*D*D/delta/delta/delta;
      w0_(m,IVZ,k,j,i) = (x*vx+y*vy)/delta;
      }
      
      if (pmy_mesh_->use_spherical_polar) {
          Real r = x1v;
          Real theta = x2v;
          Real phi = x3v;
          
          Real omega = 1.0e-2;
          Real xx = r*sin(theta)*cos(phi);
          Real yy = r*sin(theta)*sin(phi);
          Real vx = -omega*yy;
          Real vy = omega*xx;
          Real vz = 0.0;
          
          Real dens = 1.0;
          Real pres = 1.0;

          // set primitives in both newtonian and SR hydro
          w0_(m,IDN,k,j,i) = dens;
          w0_(m,IEN,k,j,i) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
          w0_(m,IVX,k,j,i) = 0.0;
          w0_(m,IVY,k,j,i) = 0.0;
          w0_(m,IVZ,k,j,i) = omega*r*sin(theta);
      }
    
  });

  // initialize magnetic fields if MHD
  if (pmbp->pmhd != nullptr) {
    // Read magnetic field strength
    Real bx = pin->GetReal("problem","b0");
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0.x1f(m,k,j,i) = bx;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = 0.0;
      if (i==ie) b0.x1f(m,k,j,i+1) = bx;
      if (j==je) b0.x2f(m,k,j+1,i) = 0.0;
      if (k==ke) b0.x3f(m,k+1,j,i) = 0.0;
      bcc0(m,IBX,k,j,i) = bx;
      bcc0(m,IBY,k,j,i) = 0.0;
      bcc0(m,IBZ,k,j,i) = 0.0;
    });
  }

  // Initialize the ADM variables if enabled
  if (pmbp->padm != nullptr) {
    pmbp->padm->SetADMVariables(pmbp);
    pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
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

void UserBoundary(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  int nmb = pm->pmb_pack->nmb_thispack;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0_;
    DvceArray5D<Real> w0_;
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
    const bool use_wellbalance = pmbp->phydro->use_wellbalance;

    Real gm1;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
        gm1 = (pmbp->phydro->peos->eos_data.gamma) - 1.0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->phydro->w0;
        gm1 = (pmbp->pmhd->peos->eos_data.gamma) - 1.0;
    }
    

  par_for("usrboundary", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),0,(ng-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user || mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
          
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          
        if (pm->use_cubed_sphere) {
          Real xi = M_PI/4.0*x1v;
          Real eta = M_PI/4.0*x2v;
          Real omega = 1.0;
          
          Real x = tan(xi);
          Real y = tan(eta);
          Real delta = sqrt(1.0 + x*x + y*y);
          Real C = sqrt(1.0 + x*x);
          Real D = sqrt(1.0 + y*y);
        
        Real theta = xi+M_PI/2.0;
        Real phi = eta;
        Real xx = x3v*sin(theta)*cos(phi);
        Real yy = x3v*sin(theta)*sin(phi);
          Real vx = -omega*yy;
          Real vy = omega*xx;
          Real vz = 0.0;
          
          
        Real dens = 1.0;
        Real pres = 1.0;

        // set primitives in both newtonian and SR hydro
        w0_(m,IDN,k,j,i) = dens;
        w0_(m,IEN,k,j,i) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
          w0_(m,IVX,k,j,i) = ((1.0+y*y)*vx-x*y*vy)*C*C/delta/delta/delta;
          w0_(m,IVY,k,j,i) = (-x*y*vx+(1.0+x*x)*vy)*D*D/delta/delta/delta;
          w0_(m,IVZ,k,j,i) = 0.0;
        }
        
        if (pm->use_spherical_polar) {
            Real r = x1v;
            Real theta = x2v;
            Real phi = x3v;
            
            Real omega = 1.0e-2;
            Real xx = r*sin(theta)*cos(phi);
            Real yy = r*sin(theta)*sin(phi);
            Real vx = -omega*yy;
            Real vy = omega*xx;
            Real vz = 0.0;
            
            Real dens = 1.0;
            Real pres = 1.0;

            // set primitives in both newtonian and SR hydro
            w0_(m,IDN,k,j,i) = dens;
            w0_(m,IEN,k,j,i) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
            w0_(m,IVX,k,j,i) = 0.0;
            w0_(m,IVY,k,j,i) = 0.0;
            w0_(m,IVZ,k,j,i) = omega*r*sin(theta);
        }
        
            u0_(m,IDN,k,j,i) = w0_(m,IDN,k,j,i);
            u0_(m,IM1,k,j,i) = w0_(m,IDN,k,j,i)*w0_(m,IVX,k,j,i);
            u0_(m,IM2,k,j,i) = w0_(m,IDN,k,j,i)*w0_(m,IVY,k,j,i);
            u0_(m,IM3,k,j,i) = w0_(m,IDN,k,j,i)*w0_(m,IVZ,k,j,i);
            u0_(m,IEN,k,j,i) = w0_(m,IEN,k,j,i) + 0.5*w0_(m,IDN,k,j,i)*(SQR(w0_(m,IVX,k,j,i))+SQR(w0_(m,IVY,k,j,i)));
    }
  });
  return;
}
