//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sph_mhd_balance.cpp
//  \brief Problem generator for MHD blast wave problem in spherical-polar coordinates

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
//  \brief Problem Generator for MHD blast wave tests

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &panel = pmbp->pmb->mb_panel;
    
  auto &use_spherical_polar = pmbp->pmesh->use_spherical_polar;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &area3 = pmbp->pcoord->area.x3f;
  auto &dxe1 = pmbp->pcoord->dxedge.x1e;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;

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
  auto &w0_ = (pmbp->phydro != nullptr)? pmbp->phydro->w0 : pmbp->pmhd->w0;
    
  Real rho0 = pin->GetReal("problem","rho0");
  Real p0 = pin->GetReal("problem","p0");
  Real r0 = pin->GetReal("problem","r0");
  Real E0 = pin->GetReal("problem","E0");
  Real b00 = pin->GetReal("problem","b0");
  Real angle = pin->GetReal("problem","angle");
  angle *= M_PI/2.0/90.0;
  Real x0 = pin->GetReal("problem","x0");
  Real y0 = pin->GetReal("problem","y0");
  Real z0 = pin->GetReal("problem","z0");
    
  Real V0 = 4.0/3.0*M_PI*r0*r0*r0;
  Real p1 = E0*gm1/V0;

  // initialize primitive variables
  par_for("pgen_prob", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
        
    Real x1v, x2v, x3v;
    Real x, y, z;
    Real x00, y00, z00;
    if (use_spherical_polar) {
      x1v = x1v_(m,i);
      x2v = x2v_(m,j);
      x3v = x3v_(m,k);
      x = x1v*sin(x2v)*cos(x3v);
      y = x1v*sin(x2v)*sin(x3v);
      z = x1v*cos(x2v);
      x00 = x0*sin(y0)*cos(z0);
      y00 = x0*sin(y0)*sin(z0);
      z00 = x0*cos(y0);
    } else {
      x1v = CellCenterX(i-is, nx1, x1min, x1max);
      x2v = CellCenterX(j-js, nx2, x2min, x2max);
      x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      x = x1v;
      y = x2v;
      z = x3v;
      x00 = x0;
      y00 = y0;
      z00 = z0;
    }
      
    Real r = sqrt(SQR(x-x00)+SQR(y-y00)+SQR(z-z00));
    Real p = (r < r0) ? p1 : p0;

    // set primitives in hydro
    w0_(m,IDN,k,j,i) = rho0;
    w0_(m,IEN,k,j,i) = p/gm1;
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
        
      if (use_spherical_polar) {
        Real x1v = x1v_(m,i);
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        Real x1fl = x1f_(m,i);
        Real x2fl = x2f_(m,j);
        Real x3fl = x3f_(m,k);
        Real x1fr = x1f_(m,i+1);
        Real x2fr = x2f_(m,j+1);
        Real x3fr = x3f_(m,k+1);
          
        Real A1 = 0.0;
//        Real A2 = -0.5*b00*x1fl*sin(x3fl-angle);
//        Real A2ip = -0.5*b00*x1fr*sin(x3fl-angle);
//        Real A2kp = -0.5*b00*x1fl*sin(x3fr-angle);
//        Real A2ipkp = -0.5*b00*x1fr*sin(x3fr-angle);
//        Real A3 = -0.5*b00*x1fl*cos(x2fl)*cos(x3v-angle);
//        Real A3ip = -0.5*b00*x1fr*cos(x2fl)*cos(x3v-angle);
//        Real A3jp = -0.5*b00*x1fl*cos(x2fr)*cos(x3v-angle);
//        Real A3ipjp = -0.5*b00*x1fr*cos(x2fr)*cos(x3v-angle);
        Real A2 = -0.5*b00*x1fl*cos(angle)*sin(x3fl);
        Real A2ip = -0.5*b00*x1fr*cos(angle)*sin(x3fl);
        Real A2kp = -0.5*b00*x1fl*cos(angle)*sin(x3fr);
        Real A2ipkp = -0.5*b00*x1fr*cos(angle)*sin(x3fr);
        Real A3 = 0.5*b00*x1fl*(sin(angle)*sin(x2fl)-cos(angle)*cos(x2fl)*cos(x3v));
          Real A3ip = 0.5*b00*x1fr*(sin(angle)*sin(x2fl)-cos(angle)*cos(x2fl)*cos(x3v));
          Real A3jp = 0.5*b00*x1fl*(sin(angle)*sin(x2fr)-cos(angle)*cos(x2fr)*cos(x3v));
          Real A3ipjp = 0.5*b00*x1fr*(sin(angle)*sin(x2fr)-cos(angle)*cos(x2fr)*cos(x3v));
          
        b0.x1f(m,k,j,i) = (dxe3(m,k,j+1,i)*A3jp - dxe3(m,k,j,i)*A3)/area1(m,k,j,i) - (dxe2(m,k+1,j,i)*A2kp - dxe2(m,k,j,i)*A2)/area1(m,k,j,i);
        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) {
          b0.x2f(m,k,j,i) = - ((x1fr*(x3fr-x3fl))*A3ip - (x1fl*(x3fr-x3fl))*A3) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
        } else {
          b0.x2f(m,k,j,i) = - (dxe3(m,k,j,i+1)*A3ip - dxe3(m,k,j,i)*A3)/area2(m,k,j,i);
        }
        b0.x3f(m,k,j,i) = (dxe2(m,k,j,i+1)*A2ip - dxe2(m,k,j,i)*A2)/area3(m,k,j,i);
        Real b0x1fip = (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j,i+1)*A3ip)/area1(m,k,j,i+1) - (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k,j,i+1)*A2ip)/area1(m,k,j,i+1);
        Real b0x2fjp;
        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je) {
          b0x2fjp = - ((x1fr*(x3fr-x3fl))*A3ipjp - (x1fl*(x3fr-x3fl))*A3jp) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
        } else {
          b0x2fjp = - (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j+1,i)*A3jp)/area2(m,k,j+1,i);
        }
        Real b0x3fkp = (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k+1,j,i)*A2kp)/area3(m,k+1,j,i);
        if (i==ie) b0.x1f(m,k,j,i+1) = b0x1fip;
        if (j==je) b0.x2f(m,k,j+1,i) = b0x2fjp;
        if (k==ke) b0.x3f(m,k+1,j,i) = b0x3fkp;
          
        Real lw, rw;
        lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
        lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0x2fjp;
        lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0x3fkp;
            
      } else {
        Real bx = b00*cos(angle);
        Real by = b00*sin(angle);
        Real bz = 0.0;
    
        bcc0(m,IBX,k,j,i) = bx;
        bcc0(m,IBY,k,j,i) = by;
        bcc0(m,IBZ,k,j,i) = bz;
      
        b0.x1f(m,k,j,i) = bx;
        b0.x2f(m,k,j,i) = by;
        b0.x3f(m,k,j,i) = bz;
        if (i==ie) b0.x1f(m,k,j,i+1) = bx;
        if (j==je) b0.x2f(m,k,j+1,i) = by;
        if (k==ke) b0.x3f(m,k+1,j,i) = bz;
      }

    });
  }

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
