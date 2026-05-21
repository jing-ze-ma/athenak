//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sph_mhd_balance.cpp
//  \brief Problem generator for steady state MHD problem in spherical-polar coordinates

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
//  \brief Problem Generator for spherical MHD balance tests

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
    
    auto &use_spherical_polar = pmbp->pmesh->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;

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
//    Real &x1min = size.d_view(m).x1min;
//    Real &x1max = size.d_view(m).x1max;
//    int nx1 = indcs.nx1;
//    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
//
//    Real &x2min = size.d_view(m).x2min;
//    Real &x2max = size.d_view(m).x2max;
//    int nx2 = indcs.nx2;
//    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
//
//      Real &x3min = size.d_view(m).x3min;
//      Real &x3max = size.d_view(m).x3max;
//      int nx3 = indcs.nx3;
//      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      
//      if (use_spherical_polar) {
          Real x1v = x1v_(m,i);
          Real x2v = x2v_(m,j);
          Real x3v = x3v_(m,k);
          
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
//      }
    
  });

  // initialize magnetic fields if MHD
  if (pmbp->pmhd != nullptr) {
    // Read magnetic field strength
//    Real bx = pin->GetReal("problem","b0");
    Real Bini = 0.0e-4;
    auto &b0 = pmbp->pmhd->b0;
    auto &bcc0 = pmbp->pmhd->bcc0;
    par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
//      Real &x1min = size.d_view(m).x1min;
//      Real &x1max = size.d_view(m).x1max;
//      int nx1 = indcs.nx1;
//      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
//
//        Real &x2min = size.d_view(m).x2min;
//        Real &x2max = size.d_view(m).x2max;
//        int nx2 = indcs.nx2;
//        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        
        Real x1v = x1v_(m,i);
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
      bcc0(m,IBX,k,j,i) = 0.0;
      bcc0(m,IBY,k,j,i) = 0.0;
      bcc0(m,IBZ,k,j,i) = Bini/(x1v*sin(x2v));
        
      Real x1f = x1f_(m,i);
      b0.x1f(m,k,j,i) = 0.0;
      b0.x2f(m,k,j,i) = 0.0;
      b0.x3f(m,k,j,i) = Bini/(x1v*sin(x2v));
      Real x1fr = x1f_(m,i+1);
      if (i==ie) b0.x1f(m,k,j,i+1) = 0.0;
      if (j==je) b0.x2f(m,k,j+1,i) = 0.0;
      if (k==ke) b0.x3f(m,k+1,j,i) = Bini/(x1v*sin(x2v));
    });
  }

//  // Initialize the ADM variables if enabled
//  if (pmbp->padm != nullptr) {
//    pmbp->padm->SetADMVariables(pmbp);
//    pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
//  }

  // Convert primitives to conserved
//  if (pmbp->padm == nullptr) {
    if (pmbp->phydro != nullptr) {
      auto &u0_ = pmbp->phydro->u0;
      pmbp->phydro->peos->PrimToCons(w0_, u0_, is, ie, js, je, ks, ke);
    } else if (pmbp->pmhd != nullptr) {
      auto &u0_ = pmbp->pmhd->u0;
      auto &bcc0_ = pmbp->pmhd->bcc0;
      pmbp->pmhd->peos->PrimToCons(w0_, bcc0_, u0_, is, ie, js, je, ks, ke);
    }
//  }

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
    
    auto &use_spherical_polar = pmbp->pmesh->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;

    DvceArray5D<Real> u0_;
    DvceArray5D<Real> w0_;
    bool use_etotgrav;

    Real gm1;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
        gm1 = (pmbp->phydro->peos->eos_data.gamma) - 1.0;
        use_etotgrav = pmbp->phydro->use_etotgrav;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->pmhd->w0;
//      b0_ = pmbp->pmhd->b0;
//      bcc0_ = pmbp->pmhd->bcc0;
        gm1 = (pmbp->pmhd->peos->eos_data.gamma) - 1.0;
        use_etotgrav = pmbp->pmhd->use_etotgrav;
    }
    

  par_for("usrboundaryx1", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
//      Real &x2min = size.d_view(m).x2min;
//      Real &x2max = size.d_view(m).x2max;
//      int nx2 = indcs.nx2;
//      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
//
//      Real &x3min = size.d_view(m).x3min;
//      Real &x3max = size.d_view(m).x3max;
//      int nx3 = indcs.nx3;
//      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      
      Real x2v = x2v_(m,j);
      Real x3v = x3v_(m,k);
      
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int i=0; i<ng; ++i) {
//          Real &x1min = size.d_view(m).x1min;
//          Real &x1max = size.d_view(m).x1max;
//          int nx1 = indcs.nx1;
//          Real x1v = CellCenterX((is-i-1)-is, nx1, x1min, x1max);
            
            Real x1v = x1v_(m,(is-i-1));
        
//        if (pm->use_spherical_polar) {
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
            w0_(m,IDN,k,j,(is-i-1)) = dens;
            w0_(m,IEN,k,j,(is-i-1)) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
            w0_(m,IVX,k,j,(is-i-1)) = 0.0;
            w0_(m,IVY,k,j,(is-i-1)) = 0.0;
            w0_(m,IVZ,k,j,(is-i-1)) = omega*r*sin(theta);
//        }
        
            u0_(m,IDN,k,j,(is-i-1)) = w0_(m,IDN,k,j,(is-i-1));
            u0_(m,IM1,k,j,(is-i-1)) = w0_(m,IDN,k,j,(is-i-1))*w0_(m,IVX,k,j,(is-i-1));
            u0_(m,IM2,k,j,(is-i-1)) = w0_(m,IDN,k,j,(is-i-1))*w0_(m,IVY,k,j,(is-i-1));
            u0_(m,IM3,k,j,(is-i-1)) = w0_(m,IDN,k,j,(is-i-1))*w0_(m,IVZ,k,j,(is-i-1));
            u0_(m,IEN,k,j,(is-i-1)) = w0_(m,IEN,k,j,(is-i-1)) + 0.5*w0_(m,IDN,k,j,(is-i-1))*(SQR(w0_(m,IVX,k,j,(is-i-1)))+SQR(w0_(m,IVY,k,j,(is-i-1)))+SQR(w0_(m,IVZ,k,j,(is-i-1))));
        }
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
        for (int i=0; i<ng; ++i) {
//          Real &x1min = size.d_view(m).x1min;
//          Real &x1max = size.d_view(m).x1max;
//          int nx1 = indcs.nx1;
//          Real x1v = CellCenterX((ie+i+1)-is, nx1, x1min, x1max);
            
            Real x1v = x1v_(m,(ie+i+1));
        
//        if (pm->use_spherical_polar) {
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
            w0_(m,IDN,k,j,(ie+i+1)) = dens;
            w0_(m,IEN,k,j,(ie+i+1)) = (pres + 0.5*dens*(vx*vx+vy*vy))/gm1;
            w0_(m,IVX,k,j,(ie+i+1)) = 0.0;
            w0_(m,IVY,k,j,(ie+i+1)) = 0.0;
            w0_(m,IVZ,k,j,(ie+i+1)) = omega*r*sin(theta);
//        }
        
            u0_(m,IDN,k,j,(ie+i+1)) = w0_(m,IDN,k,j,(ie+i+1));
            u0_(m,IM1,k,j,(ie+i+1)) = w0_(m,IDN,k,j,(ie+i+1))*w0_(m,IVX,k,j,(ie+i+1));
            u0_(m,IM2,k,j,(ie+i+1)) = w0_(m,IDN,k,j,(ie+i+1))*w0_(m,IVY,k,j,(ie+i+1));
            u0_(m,IM3,k,j,(ie+i+1)) = w0_(m,IDN,k,j,(ie+i+1))*w0_(m,IVZ,k,j,(ie+i+1));
            u0_(m,IEN,k,j,(ie+i+1)) = w0_(m,IEN,k,j,(ie+i+1)) + 0.5*w0_(m,IDN,k,j,(ie+i+1))*(SQR(w0_(m,IVX,k,j,(ie+i+1)))+SQR(w0_(m,IVY,k,j,(ie+i+1)))+SQR(w0_(m,IVZ,k,j,(ie+i+1))));
        }
      }
  });
    
  if (pmbp->pmhd != nullptr) {
    auto &b0_ = pmbp->pmhd->b0;
    auto &bcc0_ = pmbp->pmhd->bcc0;
      
      par_for("usrboundarymhd", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
//          Real &x2min = size.d_view(m).x2min;
//          Real &x2max = size.d_view(m).x2max;
//          int nx2 = indcs.nx2;
//          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
//
//            Real &x3min = size.d_view(m).x3min;
//            Real &x3max = size.d_view(m).x3max;
//            int nx3 = indcs.nx3;
//            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          
          Real x2v = x2v_(m,j);
          Real x3v = x3v_(m,k);
          
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int i=0; i<ng; ++i) {
//            Real &x1min = size.d_view(m).x1min;
//            Real &x1max = size.d_view(m).x1max;
//            int nx1 = indcs.nx1;
//            Real x1v = CellCenterX((is-i-1)-is, nx1, x1min, x1max);
              
              Real x1v = x1v_(m,(is-i-1));
            
            Real Bini = 0.0e-4;
                
              bcc0_(m,IBX,k,j,(is-i-1)) = 0.0;
              bcc0_(m,IBY,k,j,(is-i-1)) = 0.0;
              bcc0_(m,IBZ,k,j,(is-i-1)) = Bini/(x1v*sin(x2v));
                
              b0_.x1f(m,k,j,(is-i-1)) = 0.0;
              b0_.x2f(m,k,j,(is-i-1)) = 0.0;
              b0_.x3f(m,k,j,(is-i-1)) = Bini/(x1v*sin(x2v));
              if (j==n2-1) b0_.x2f(m,k,j+1,(is-i-1)) = 0.0;
              if (k==n3-1) b0_.x3f(m,k+1,j,(is-i-1)) = Bini/(x1v*sin(x2v));
            
            u0_(m,IEN,k,j,(is-i-1)) += 0.5*(SQR(bcc0_(m,IBX,k,j,(is-i-1)))+SQR(bcc0_(m,IBY,k,j,(is-i-1)))+SQR(bcc0_(m,IBZ,k,j,(is-i-1))));
          }
        }
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
//              Real &x1min = size.d_view(m).x1min;
//              Real &x1max = size.d_view(m).x1max;
//              int nx1 = indcs.nx1;
//              Real x1v = CellCenterX((ie+i+1)-is, nx1, x1min, x1max);
                
                Real x1v = x1v_(m,(ie+i+1));
              
              Real Bini = 0.0e-4;
                  
                bcc0_(m,IBX,k,j,(ie+i+1)) = 0.0;
                bcc0_(m,IBY,k,j,(ie+i+1)) = 0.0;
                bcc0_(m,IBZ,k,j,(ie+i+1)) = Bini/(x1v*sin(x2v));
                  
                b0_.x1f(m,k,j,(ie+i+1)) = 0.0;
                b0_.x2f(m,k,j,(ie+i+1)) = 0.0;
                b0_.x3f(m,k,j,(ie+i+1)) = Bini/(x1v*sin(x2v));
                if (j==n2-1) b0_.x2f(m,k,j+1,(ie+i+1)) = 0.0;
                if (k==n3-1) b0_.x3f(m,k+1,j,(ie+i+1)) = Bini/(x1v*sin(x2v));
              
              u0_(m,IEN,k,j,(ie+i+1)) += 0.5*(SQR(bcc0_(m,IBX,k,j,(ie+i+1)))+SQR(bcc0_(m,IBY,k,j,(ie+i+1)))+SQR(bcc0_(m,IBZ,k,j,(ie+i+1))));
            }
          }
      });
  }
  return;
}
