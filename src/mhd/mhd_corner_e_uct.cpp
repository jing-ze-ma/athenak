//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_corner_e_uct.cpp
//  \brief
//  Reconstructing corner E with the UCT method (Mignone & Del Zanna 2021)

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "diffusion/resistivity.hpp"
#include "mhd.hpp"

#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn  void MHD::CornerE_UCT
//  \brief calculate the corner electric fields.

TaskStatus MHD::CornerE_UCT(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;

  //---- 3-D problem:
  // Use MZ21 algorithm to compute all three of E1, E2, and E3

  if (pmy_pack->pmesh->three_d) {
    // E1=-(v X B)=VzBy-VyBz
    // E2=-(v X B)=VxBz-VzBx
    // E3=-(v X B)=VyBx-VxBy

    // capture class variables for the kernels
    auto b1 = b0.x1f;
    auto b2 = b0.x2f;
    auto b3 = b0.x3f;
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
      auto v31l_ = v31l;
      auto v31r_ = v31r;
      auto v21l_ = v21l;
      auto v21r_ = v21r;
      auto v11l_ = v11l;
      auto v11r_ = v11r;
      auto nus1_ = nus1;
      auto d1l_ = d1l;
      auto d1r_ = d1r;
      auto lam1l_ = lam1l;
      auto lam1r_ = lam1r;
      auto v12l_ = v12l;
      auto v12r_ = v12r;
      auto v32l_ = v32l;
      auto v32r_ = v32r;
      auto v22l_ = v22l;
      auto v22r_ = v22r;
      auto nus2_ = nus2;
      auto d2l_ = d2l;
      auto d2r_ = d2r;
      auto lam2l_ = lam2l;
      auto lam2r_ = lam2r;
      auto v23l_ = v23l;
      auto v23r_ = v23r;
      auto v13l_ = v13l;
      auto v13r_ = v13r;
      auto v33l_ = v33l;
      auto v33r_ = v33r;
      auto nus3_ = nus3;
      auto d3l_ = d3l;
      auto d3r_ = d3r;
      auto lam3l_ = lam3l;
      auto lam3r_ = lam3r;
      auto &x1f_ = pmy_pack->pcoord->xx1f;
      auto &x1v_ = pmy_pack->pcoord->x1v;
      auto &x2f_ = pmy_pack->pcoord->xx2f;
      auto &x2v_ = pmy_pack->pcoord->x2v;
      auto &x3f_ = pmy_pack->pcoord->xx3f;
      auto &x3v_ = pmy_pack->pcoord->x3v;
      
    par_for("emf3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real al,ar,alm,arm,axW,axE,ayS,ayN,dxW,dxE,dyS,dyN;
      Real vyN,vyS,vxE,vxW,ByE,ByW,BxN,BxS;
      Real alpxp,alpxm,alpyp,alpym;
        
      Real nus1c = nus1_(m,k,j,i);
      Real nus2c = nus2_(m,k,j,i);
      Real nus3c = nus3_(m,k,j,i);
      Real d1lc = d1l_(m,k,j,i);
      Real d1rc = d1r_(m,k,j,i);
      Real d2lc = d2l_(m,k,j,i);
      Real d2rc = d2r_(m,k,j,i);
      Real d3lc = d3l_(m,k,j,i);
      Real d3rc = d3r_(m,k,j,i);
      Real b1c = b1(m,k,j,i);
      Real b2c = b2(m,k,j,i);
      Real b3c = b3(m,k,j,i);
        
      // reconstruct E3 to corner using FS18
      Real dum;
      Real b1l2,b1r2,b2l1,b2r1;
      Real v1l1l2,v1r1l2,v1l1r2,v1r1r2,v1l2l1,v1r2l1,v1l2r1,v1r2r1;
      Real v2l1l2,v2r1l2,v2l1r2,v2r1r2,v2l2l1,v2r2l1,v2l2r1,v2r2r1;
      Real v1l1l2a,v1l1r2a,v1r1l2a,v1r1r2a,v2l1l2a,v2l1r2a,v2r1l2a,v2r1r2a;
      Real e3l1l2,e3l1r2,e3r1l2,e3r1r2;
      Real al1l2,ar1l2,al1r2,ar1r2,adenom,d1,d2;
      Real dxl2l,dxr2l,dxlh2l,dxrh2l,dxl2r,dxr2r,dxlh2r,dxrh2r,dxl1l,dxr1l,dxlh1l,dxrh1l,dxl1r,dxr1r,dxlh1r,dxrh1r;
        
      if (use_spherical_polar) {
          dxl2l = x2v_(m,j-1)-x2v_(m,j-2);
          dxr2l = x2v_(m,j)-x2v_(m,j-1);
          dxlh2l = x2v_(m,j-1)-x2f_(m,j-1);
          dxrh2l = x2f_(m,j)-x2v_(m,j-1);
          dxl2r = x2v_(m,j)-x2v_(m,j-1);
          dxr2r = x2v_(m,j+1)-x2v_(m,j);
          dxlh2r = x2v_(m,j)-x2f_(m,j);
          dxrh2r = x2f_(m,j+1)-x2v_(m,j);
          dxl1l = x1v_(m,i-1)-x1v_(m,i-2);
          dxr1l = x1v_(m,i)-x1v_(m,i-1);
          dxlh1l = x1v_(m,i-1)-x1f_(m,i-1);
          dxrh1l = x1f_(m,i)-x1v_(m,i-1);
          dxl1r = x1v_(m,i)-x1v_(m,i-1);
          dxr1r = x1v_(m,i+1)-x1v_(m,i);
          dxlh1r = x1v_(m,i)-x1f_(m,i);
          dxrh1r = x1f_(m,i+1)-x1v_(m,i);
          
          PLM_nonuniform(b1(m,k,j-2,i),b1(m,k,j-1,i),b1(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,b1l2,dum);
          PLM_nonuniform(b1(m,k,j-1,i),b1(m,k,j,i),b1(m,k,j+1,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,b1r2);
          PLM_nonuniform(b2(m,k,j,i-2),b2(m,k,j,i-1),b2(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,b2l1,dum);
          PLM_nonuniform(b2(m,k,j,i-1),b2(m,k,j,i),b2(m,k,j,i+1),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,b2r1);
            
          PLM_nonuniform(v11l_(m,k,j-2,i),v11l_(m,k,j-1,i),v11l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1l1l2,dum);
          PLM_nonuniform(v11r_(m,k,j-2,i),v11r_(m,k,j-1,i),v11r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1r1l2,dum);
          PLM_nonuniform(v21l_(m,k,j-2,i),v21l_(m,k,j-1,i),v21l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2l1l2,dum);
          PLM_nonuniform(v21r_(m,k,j-2,i),v21r_(m,k,j-1,i),v21r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2r1l2,dum);
          PLM_nonuniform(v11l_(m,k,j-1,i),v11l_(m,k,j,i),v11l_(m,k,j+1,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1l1r2);
          PLM_nonuniform(v11r_(m,k,j-1,i),v11r_(m,k,j,i),v11r_(m,k,j+1,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1r1r2);
          PLM_nonuniform(v21l_(m,k,j-1,i),v21l_(m,k,j,i),v21l_(m,k,j+1,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2l1r2);
          PLM_nonuniform(v21r_(m,k,j-1,i),v21r_(m,k,j,i),v21r_(m,k,j+1,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2r1r2);
          PLM_nonuniform(v12l_(m,k,j,i-2),v12l_(m,k,j,i-1),v12l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1l2l1,dum);
          PLM_nonuniform(v12r_(m,k,j,i-2),v12r_(m,k,j,i-1),v12r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1r2l1,dum);
          PLM_nonuniform(v22l_(m,k,j,i-2),v22l_(m,k,j,i-1),v22l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2l2l1,dum);
          PLM_nonuniform(v22r_(m,k,j,i-2),v22r_(m,k,j,i-1),v22r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2r2l1,dum);
          PLM_nonuniform(v12l_(m,k,j,i-1),v12l_(m,k,j,i),v12l_(m,k,j,i+1),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1l2r1);
          PLM_nonuniform(v12r_(m,k,j,i-1),v12r_(m,k,j,i),v12r_(m,k,j,i+1),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1r2r1);
          PLM_nonuniform(v22l_(m,k,j,i-1),v22l_(m,k,j,i),v22l_(m,k,j,i+1),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2l2r1);
          PLM_nonuniform(v22r_(m,k,j,i-1),v22r_(m,k,j,i),v22r_(m,k,j,i+1),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2r2r1);
      } else {
          PLM(b1(m,k,j-2,i),b1(m,k,j-1,i),b1(m,k,j,i),b1l2,dum);
          PLM(b1(m,k,j-1,i),b1(m,k,j,i),b1(m,k,j+1,i),dum,b1r2);
          PLM(b2(m,k,j,i-2),b2(m,k,j,i-1),b2(m,k,j,i),b2l1,dum);
          PLM(b2(m,k,j,i-1),b2(m,k,j,i),b2(m,k,j,i+1),dum,b2r1);
            
          PLM(v11l_(m,k,j-2,i),v11l_(m,k,j-1,i),v11l_(m,k,j,i),v1l1l2,dum);
          PLM(v11r_(m,k,j-2,i),v11r_(m,k,j-1,i),v11r_(m,k,j,i),v1r1l2,dum);
          PLM(v21l_(m,k,j-2,i),v21l_(m,k,j-1,i),v21l_(m,k,j,i),v2l1l2,dum);
          PLM(v21r_(m,k,j-2,i),v21r_(m,k,j-1,i),v21r_(m,k,j,i),v2r1l2,dum);
          PLM(v11l_(m,k,j-1,i),v11l_(m,k,j,i),v11l_(m,k,j+1,i),dum,v1l1r2);
          PLM(v11r_(m,k,j-1,i),v11r_(m,k,j,i),v11r_(m,k,j+1,i),dum,v1r1r2);
          PLM(v21l_(m,k,j-1,i),v21l_(m,k,j,i),v21l_(m,k,j+1,i),dum,v2l1r2);
          PLM(v21r_(m,k,j-1,i),v21r_(m,k,j,i),v21r_(m,k,j+1,i),dum,v2r1r2);
          PLM(v12l_(m,k,j,i-2),v12l_(m,k,j,i-1),v12l_(m,k,j,i),v1l2l1,dum);
          PLM(v12r_(m,k,j,i-2),v12r_(m,k,j,i-1),v12r_(m,k,j,i),v1r2l1,dum);
          PLM(v22l_(m,k,j,i-2),v22l_(m,k,j,i-1),v22l_(m,k,j,i),v2l2l1,dum);
          PLM(v22r_(m,k,j,i-2),v22r_(m,k,j,i-1),v22r_(m,k,j,i),v2r2l1,dum);
          PLM(v12l_(m,k,j,i-1),v12l_(m,k,j,i),v12l_(m,k,j,i+1),dum,v1l2r1);
          PLM(v12r_(m,k,j,i-1),v12r_(m,k,j,i),v12r_(m,k,j,i+1),dum,v1r2r1);
          PLM(v22l_(m,k,j,i-1),v22l_(m,k,j,i),v22l_(m,k,j,i+1),dum,v2l2r1);
          PLM(v22r_(m,k,j,i-1),v22r_(m,k,j,i),v22r_(m,k,j,i+1),dum,v2r2r1);
      }
        
      v1l1l2a = 0.5*(v1l1l2+v1l2l1);
      v1l1r2a = 0.5*(v1l1r2+v1r2l1);
      v1r1l2a = 0.5*(v1r1l2+v1l2r1);
      v1r1r2a = 0.5*(v1r1r2+v1r2r1);
      v2l1l2a = 0.5*(v2l1l2+v2l2l1);
      v2l1r2a = 0.5*(v2l1r2+v2r2l1);
      v2r1l2a = 0.5*(v2r1l2+v2l2r1);
      v2r1r2a = 0.5*(v2r1r2+v2r2r1);
      
      e3l1l2 = v2l1l2a*b1l2 - v1l1l2a*b2l1;
      e3l1r2 = v2l1r2a*b1r2 - v1l1r2a*b2l1;
      e3r1l2 = v2r1l2a*b1l2 - v1r1l2a*b2r1;
      e3r1r2 = v2r1r2a*b1r2 - v1r1r2a*b2r1;
        
      alpxp = fmax(0.0,fmax(lam1r_(m,k,j-1,i),lam1r_(m,k,j,i)));
      alpxm = -fmin(0.0,fmin(lam1l_(m,k,j-1,i),lam1l_(m,k,j,i)));
      alpyp = fmax(0.0,fmax(lam2r_(m,k,j,i-1),lam2r_(m,k,j,i)));
      alpym = -fmin(0.0,fmin(lam2l_(m,k,j,i-1),lam2l_(m,k,j,i)));
    
      al1l2 = alpxp*alpyp;
      ar1l2 = alpxm*alpyp;
      al1r2 = alpxp*alpym;
      ar1r2 = alpxm*alpym;
      adenom = (alpxp+alpxm)*(alpyp+alpym);
      al1l2 /= adenom;
      ar1l2 /= adenom;
      al1r2 /= adenom;
      ar1r2 /= adenom;
      d1 = -alpyp*alpym/(alpyp+alpym);
      d2 = alpxp*alpxm/(alpxp+alpxm);
      
      e3(m,k,j,i) = al1l2*e3l1l2 + al1r2*e3l1r2 + ar1l2*e3r1l2 + ar1r2*e3r1r2 + d1*(b1r2-b1l2) + d2*(b2r1-b2l1);
        
        // reconstruct E1 to corner using FS18
          
        if (use_spherical_polar) {
            dxl2l = x3v_(m,k-1)-x3v_(m,k-2);
            dxr2l = x3v_(m,k)-x3v_(m,k-1);
            dxlh2l = x3v_(m,k-1)-x3f_(m,k-1);
            dxrh2l = x3f_(m,k)-x3v_(m,k-1);
            dxl2r = x3v_(m,k)-x3v_(m,k-1);
            dxr2r = x3v_(m,k+1)-x3v_(m,k);
            dxlh2r = x3v_(m,k)-x3f_(m,k);
            dxrh2r = x3f_(m,k+1)-x3v_(m,k);
            dxl1l = x2v_(m,j-1)-x2v_(m,j-2);
            dxr1l = x2v_(m,j)-x2v_(m,j-1);
            dxlh1l = x2v_(m,j-1)-x2f_(m,j-1);
            dxrh1l = x2f_(m,j)-x2v_(m,j-1);
            dxl1r = x2v_(m,j)-x2v_(m,j-1);
            dxr1r = x2v_(m,j+1)-x2v_(m,j);
            dxlh1r = x2v_(m,j)-x2f_(m,j);
            dxrh1r = x2f_(m,j+1)-x2v_(m,j);
            
            PLM_nonuniform(b2(m,k-2,j,i),b2(m,k-1,j,i),b2(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,b1l2,dum);
            PLM_nonuniform(b2(m,k-1,j,i),b2(m,k,j,i),b2(m,k+1,j,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,b1r2);
            PLM_nonuniform(b3(m,k,j-2,i),b3(m,k,j-1,i),b3(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,b2l1,dum);
            PLM_nonuniform(b3(m,k,j-1,i),b3(m,k,j,i),b3(m,k,j+1,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,b2r1);
              
            PLM_nonuniform(v22l_(m,k-2,j,i),v22l_(m,k-1,j,i),v22l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1l1l2,dum);
            PLM_nonuniform(v22r_(m,k-2,j,i),v22r_(m,k-1,j,i),v22r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1r1l2,dum);
            PLM_nonuniform(v32l_(m,k-2,j,i),v32l_(m,k-1,j,i),v32l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2l1l2,dum);
            PLM_nonuniform(v32r_(m,k-2,j,i),v32r_(m,k-1,j,i),v32r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2r1l2,dum);
            PLM_nonuniform(v22l_(m,k-1,j,i),v22l_(m,k,j,i),v22l_(m,k+1,j,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1l1r2);
            PLM_nonuniform(v22r_(m,k-1,j,i),v22r_(m,k,j,i),v22r_(m,k+1,j,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1r1r2);
            PLM_nonuniform(v32l_(m,k-1,j,i),v32l_(m,k,j,i),v32l_(m,k+1,j,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2l1r2);
            PLM_nonuniform(v32r_(m,k-1,j,i),v32r_(m,k,j,i),v32r_(m,k+1,j,i),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2r1r2);
            PLM_nonuniform(v23l_(m,k,j-2,i),v23l_(m,k,j-1,i),v23l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1l2l1,dum);
            PLM_nonuniform(v23r_(m,k,j-2,i),v23r_(m,k,j-1,i),v23r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1r2l1,dum);
            PLM_nonuniform(v33l_(m,k,j-2,i),v33l_(m,k,j-1,i),v33l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2l2l1,dum);
            PLM_nonuniform(v33r_(m,k,j-2,i),v33r_(m,k,j-1,i),v33r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2r2l1,dum);
            PLM_nonuniform(v23l_(m,k,j-1,i),v23l_(m,k,j,i),v23l_(m,k,j+1,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1l2r1);
            PLM_nonuniform(v23r_(m,k,j-1,i),v23r_(m,k,j,i),v23r_(m,k,j+1,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1r2r1);
            PLM_nonuniform(v33l_(m,k,j-1,i),v33l_(m,k,j,i),v33l_(m,k,j+1,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2l2r1);
            PLM_nonuniform(v33r_(m,k,j-1,i),v33r_(m,k,j,i),v33r_(m,k,j+1,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2r2r1);
        } else {
            PLM(b2(m,k-2,j,i),b2(m,k-1,j,i),b2(m,k,j,i),b1l2,dum);
            PLM(b2(m,k-1,j,i),b2(m,k,j,i),b2(m,k+1,j,i),dum,b1r2);
            PLM(b3(m,k,j-2,i),b3(m,k,j-1,i),b3(m,k,j,i),b2l1,dum);
            PLM(b3(m,k,j-1,i),b3(m,k,j,i),b3(m,k,j+1,i),dum,b2r1);
              
            PLM(v22l_(m,k-2,j,i),v22l_(m,k-1,j,i),v22l_(m,k,j,i),v1l1l2,dum);
            PLM(v22r_(m,k-2,j,i),v22r_(m,k-1,j,i),v22r_(m,k,j,i),v1r1l2,dum);
            PLM(v32l_(m,k-2,j,i),v32l_(m,k-1,j,i),v32l_(m,k,j,i),v2l1l2,dum);
            PLM(v32r_(m,k-2,j,i),v32r_(m,k-1,j,i),v32r_(m,k,j,i),v2r1l2,dum);
            PLM(v22l_(m,k-1,j,i),v22l_(m,k,j,i),v22l_(m,k+1,j,i),dum,v1l1r2);
            PLM(v22r_(m,k-1,j,i),v22r_(m,k,j,i),v22r_(m,k+1,j,i),dum,v1r1r2);
            PLM(v32l_(m,k-1,j,i),v32l_(m,k,j,i),v32l_(m,k+1,j,i),dum,v2l1r2);
            PLM(v32r_(m,k-1,j,i),v32r_(m,k,j,i),v32r_(m,k+1,j,i),dum,v2r1r2);
            PLM(v23l_(m,k,j-2,i),v23l_(m,k,j-1,i),v23l_(m,k,j,i),v1l2l1,dum);
            PLM(v23r_(m,k,j-2,i),v23r_(m,k,j-1,i),v23r_(m,k,j,i),v1r2l1,dum);
            PLM(v33l_(m,k,j-2,i),v33l_(m,k,j-1,i),v33l_(m,k,j,i),v2l2l1,dum);
            PLM(v33r_(m,k,j-2,i),v33r_(m,k,j-1,i),v33r_(m,k,j,i),v2r2l1,dum);
            PLM(v23l_(m,k,j-1,i),v23l_(m,k,j,i),v23l_(m,k,j+1,i),dum,v1l2r1);
            PLM(v23r_(m,k,j-1,i),v23r_(m,k,j,i),v23r_(m,k,j+1,i),dum,v1r2r1);
            PLM(v33l_(m,k,j-1,i),v33l_(m,k,j,i),v33l_(m,k,j+1,i),dum,v2l2r1);
            PLM(v33r_(m,k,j-1,i),v33r_(m,k,j,i),v33r_(m,k,j+1,i),dum,v2r2r1);
        }
        
        v1l1l2a = 0.5*(v1l1l2+v1l2l1);
        v1l1r2a = 0.5*(v1l1r2+v1r2l1);
        v1r1l2a = 0.5*(v1r1l2+v1l2r1);
        v1r1r2a = 0.5*(v1r1r2+v1r2r1);
        v2l1l2a = 0.5*(v2l1l2+v2l2l1);
        v2l1r2a = 0.5*(v2l1r2+v2r2l1);
        v2r1l2a = 0.5*(v2r1l2+v2l2r1);
        v2r1r2a = 0.5*(v2r1r2+v2r2r1);
        
        e3l1l2 = v2l1l2a*b1l2 - v1l1l2a*b2l1;
        e3l1r2 = v2l1r2a*b1r2 - v1l1r2a*b2l1;
        e3r1l2 = v2r1l2a*b1l2 - v1r1l2a*b2r1;
        e3r1r2 = v2r1r2a*b1r2 - v1r1r2a*b2r1;
          
        alpxp = fmax(0.0,fmax(lam2r_(m,k-1,j,i),lam2r_(m,k,j,i)));
        alpxm = -fmin(0.0,fmin(lam2l_(m,k-1,j,i),lam2l_(m,k,j,i)));
        alpyp = fmax(0.0,fmax(lam3r_(m,k,j-1,i),lam3r_(m,k,j,i)));
        alpym = -fmin(0.0,fmin(lam3l_(m,k,j-1,i),lam3l_(m,k,j,i)));
      
        al1l2 = alpxp*alpyp;
        ar1l2 = alpxm*alpyp;
        al1r2 = alpxp*alpym;
        ar1r2 = alpxm*alpym;
        adenom = (alpxp+alpxm)*(alpyp+alpym);
        al1l2 /= adenom;
        ar1l2 /= adenom;
        al1r2 /= adenom;
        ar1r2 /= adenom;
        d1 = -alpyp*alpym/(alpyp+alpym);
        d2 = alpxp*alpxm/(alpxp+alpxm);
        
//        bool do_inner =
//              (mb_bcs.d_view(m, BoundaryFace::inner_x2) ==
//               BoundaryFlag::polar) && j == js;
//        bool do_outer =
//              (mb_bcs.d_view(m, BoundaryFace::outer_x2) ==
//               BoundaryFlag::polar) && j == je+1;
//        if (do_inner || do_outer) {
//            e1(m,k,j,i) = al1l2*e3l1l2 + al1r2*e3l1r2 + ar1l2*e3r1l2 + ar1r2*e3r1r2 + d1*(b1r2-b1l2) + d2*(b2r1-b2l1);
//        } else {
        e1(m,k,j,i) = al1l2*e3l1l2 + al1r2*e3l1r2 + ar1l2*e3r1l2 + ar1r2*e3r1r2 + d1*(b1r2-b1l2) + d2*(b2r1-b2l1);
//        }
        
        // reconstruct E2 to corner using FS18
          
        if (use_spherical_polar) {
            dxl2l = x1v_(m,i-1)-x1v_(m,i-2);
            dxr2l = x1v_(m,i)-x1v_(m,i-1);
            dxlh2l = x1v_(m,i-1)-x1f_(m,i-1);
            dxrh2l = x1f_(m,i)-x1v_(m,i-1);
            dxl2r = x1v_(m,i)-x1v_(m,i-1);
            dxr2r = x1v_(m,i+1)-x1v_(m,i);
            dxlh2r = x1v_(m,i)-x1f_(m,i);
            dxrh2r = x1f_(m,i+1)-x1v_(m,i);
            dxl1l = x3v_(m,k-1)-x3v_(m,k-2);
            dxr1l = x3v_(m,k)-x3v_(m,k-1);
            dxlh1l = x3v_(m,k-1)-x3f_(m,k-1);
            dxrh1l = x3f_(m,k)-x3v_(m,k-1);
            dxl1r = x3v_(m,k)-x3v_(m,k-1);
            dxr1r = x3v_(m,k+1)-x3v_(m,k);
            dxlh1r = x3v_(m,k)-x3f_(m,k);
            dxrh1r = x3f_(m,k+1)-x3v_(m,k);
            
            PLM_nonuniform(b3(m,k,j,i-2),b3(m,k,j,i-1),b3(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,b1l2,dum);
            PLM_nonuniform(b3(m,k,j,i-1),b3(m,k,j,i),b3(m,k,j,i+1),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,b1r2);
            PLM_nonuniform(b1(m,k-2,j,i),b1(m,k-1,j,i),b1(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,b2l1,dum);
            PLM_nonuniform(b1(m,k-1,j,i),b1(m,k,j,i),b1(m,k+1,j,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,b2r1);
              
            PLM_nonuniform(v33l_(m,k,j,i-2),v33l_(m,k,j,i-1),v33l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1l1l2,dum);
            PLM_nonuniform(v33r_(m,k,j,i-2),v33r_(m,k,j,i-1),v33r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v1r1l2,dum);
            PLM_nonuniform(v13l_(m,k,j,i-2),v13l_(m,k,j,i-1),v13l_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2l1l2,dum);
            PLM_nonuniform(v13r_(m,k,j,i-2),v13r_(m,k,j,i-1),v13r_(m,k,j,i),dxl2l,dxr2l,dxlh2l,dxrh2l,v2r1l2,dum);
            PLM_nonuniform(v33l_(m,k,j,i-1),v33l_(m,k,j,i),v33l_(m,k,j,i+1),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1l1r2);
            PLM_nonuniform(v33r_(m,k,j,i-1),v33r_(m,k,j,i),v33r_(m,k,j,i+1),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v1r1r2);
            PLM_nonuniform(v13l_(m,k,j,i-1),v13l_(m,k,j,i),v13l_(m,k,j,i+1),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2l1r2);
            PLM_nonuniform(v13r_(m,k,j,i-1),v13r_(m,k,j,i),v13r_(m,k,j,i+1),dxl2r,dxr2r,dxlh2r,dxrh2r,dum,v2r1r2);
            PLM_nonuniform(v31l_(m,k-2,j,i),v31l_(m,k-1,j,i),v31l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1l2l1,dum);
            PLM_nonuniform(v31r_(m,k-2,j,i),v31r_(m,k-1,j,i),v31r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v1r2l1,dum);
            PLM_nonuniform(v11l_(m,k-2,j,i),v11l_(m,k-1,j,i),v11l_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2l2l1,dum);
            PLM_nonuniform(v11r_(m,k-2,j,i),v11r_(m,k-1,j,i),v11r_(m,k,j,i),dxl1l,dxr1l,dxlh1l,dxrh1l,v2r2l1,dum);
            PLM_nonuniform(v31l_(m,k-1,j,i),v31l_(m,k,j,i),v31l_(m,k+1,j,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1l2r1);
            PLM_nonuniform(v31r_(m,k-1,j,i),v31r_(m,k,j,i),v31r_(m,k+1,j,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v1r2r1);
            PLM_nonuniform(v11l_(m,k-1,j,i),v11l_(m,k,j,i),v11l_(m,k+1,j,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2l2r1);
            PLM_nonuniform(v11r_(m,k-1,j,i),v11r_(m,k,j,i),v11r_(m,k+1,j,i),dxl1r,dxr1r,dxlh1r,dxrh1r,dum,v2r2r1);
        } else {
            PLM(b3(m,k,j,i-2),b3(m,k,j,i-1),b3(m,k,j,i),b1l2,dum);
            PLM(b3(m,k,j,i-1),b3(m,k,j,i),b3(m,k,j,i+1),dum,b1r2);
            PLM(b1(m,k-2,j,i),b1(m,k-1,j,i),b1(m,k,j,i),b2l1,dum);
            PLM(b1(m,k-1,j,i),b1(m,k,j,i),b1(m,k+1,j,i),dum,b2r1);
              
            PLM(v33l_(m,k,j,i-2),v33l_(m,k,j,i-1),v33l_(m,k,j,i),v1l1l2,dum);
            PLM(v33r_(m,k,j,i-2),v33r_(m,k,j,i-1),v33r_(m,k,j,i),v1r1l2,dum);
            PLM(v13l_(m,k,j,i-2),v13l_(m,k,j,i-1),v13l_(m,k,j,i),v2l1l2,dum);
            PLM(v13r_(m,k,j,i-2),v13r_(m,k,j,i-1),v13r_(m,k,j,i),v2r1l2,dum);
            PLM(v33l_(m,k,j,i-1),v33l_(m,k,j,i),v33l_(m,k,j,i+1),dum,v1l1r2);
            PLM(v33r_(m,k,j,i-1),v33r_(m,k,j,i),v33r_(m,k,j,i+1),dum,v1r1r2);
            PLM(v13l_(m,k,j,i-1),v13l_(m,k,j,i),v13l_(m,k,j,i+1),dum,v2l1r2);
            PLM(v13r_(m,k,j,i-1),v13r_(m,k,j,i),v13r_(m,k,j,i+1),dum,v2r1r2);
            PLM(v31l_(m,k-2,j,i),v31l_(m,k-1,j,i),v31l_(m,k,j,i),v1l2l1,dum);
            PLM(v31r_(m,k-2,j,i),v31r_(m,k-1,j,i),v31r_(m,k,j,i),v1r2l1,dum);
            PLM(v11l_(m,k-2,j,i),v11l_(m,k-1,j,i),v11l_(m,k,j,i),v2l2l1,dum);
            PLM(v11r_(m,k-2,j,i),v11r_(m,k-1,j,i),v11r_(m,k,j,i),v2r2l1,dum);
            PLM(v31l_(m,k-1,j,i),v31l_(m,k,j,i),v31l_(m,k+1,j,i),dum,v1l2r1);
            PLM(v31r_(m,k-1,j,i),v31r_(m,k,j,i),v31r_(m,k+1,j,i),dum,v1r2r1);
            PLM(v11l_(m,k-1,j,i),v11l_(m,k,j,i),v11l_(m,k+1,j,i),dum,v2l2r1);
            PLM(v11r_(m,k-1,j,i),v11r_(m,k,j,i),v11r_(m,k+1,j,i),dum,v2r2r1);
        }
        
        v1l1l2a = 0.5*(v1l1l2+v1l2l1);
        v1l1r2a = 0.5*(v1l1r2+v1r2l1);
        v1r1l2a = 0.5*(v1r1l2+v1l2r1);
        v1r1r2a = 0.5*(v1r1r2+v1r2r1);
        v2l1l2a = 0.5*(v2l1l2+v2l2l1);
        v2l1r2a = 0.5*(v2l1r2+v2r2l1);
        v2r1l2a = 0.5*(v2r1l2+v2l2r1);
        v2r1r2a = 0.5*(v2r1r2+v2r2r1);
        
        e3l1l2 = v2l1l2a*b1l2 - v1l1l2a*b2l1;
        e3l1r2 = v2l1r2a*b1r2 - v1l1r2a*b2l1;
        e3r1l2 = v2r1l2a*b1l2 - v1r1l2a*b2r1;
        e3r1r2 = v2r1r2a*b1r2 - v1r1r2a*b2r1;
          
        alpxp = fmax(0.0,fmax(lam3r_(m,k,j,i-1),lam3r_(m,k,j,i)));
        alpxm = -fmin(0.0,fmin(lam3l_(m,k,j,i-1),lam3l_(m,k,j,i)));
        alpyp = fmax(0.0,fmax(lam1r_(m,k-1,j,i),lam1r_(m,k,j,i)));
        alpym = -fmin(0.0,fmin(lam1l_(m,k-1,j,i),lam1l_(m,k,j,i)));
      
        al1l2 = alpxp*alpyp;
        ar1l2 = alpxm*alpyp;
        al1r2 = alpxp*alpym;
        ar1r2 = alpxm*alpym;
        adenom = (alpxp+alpxm)*(alpyp+alpym);
        al1l2 /= adenom;
        ar1l2 /= adenom;
        al1r2 /= adenom;
        ar1r2 /= adenom;
        d1 = -alpyp*alpym/(alpyp+alpym);
        d2 = alpxp*alpxm/(alpxp+alpxm);
        
        e2(m,k,j,i) = al1l2*e3l1l2 + al1r2*e3l1r2 + ar1l2*e3r1l2 + ar1r2*e3r1r2 + d1*(b1r2-b1l2) + d2*(b2r1-b2l1);
    
        
        
//      // integrate E3 to corner using MZ21
////      alm = 0.5*(1.0+nus1_(m,k,j-1,i));
////      arm = 1.0 - alm;
////      al = 0.5*(1.0+nus1c);
////      ar = 1.0 - al;
////      axW = 0.5*(alm+al);
////      axE = 0.5*(arm+ar);
////      alm = 0.5*(1.0+nus2_(m,k,j,i-1));
////      arm = 1.0 - alm;
////      al = 0.5*(1.0+nus2c);
////      ar = 1.0 - al;
////      ayS = 0.5*(alm+al);
////      ayN = 0.5*(arm+ar);
////      dxW = 0.5*(d1l_(m,k,j-1,i)+d1lc);
////      dxE = 0.5*(d1r_(m,k,j-1,i)+d1rc);
////      dyS = 0.5*(d2l_(m,k,j,i-1)+d2lc);
////      dyN = 0.5*(d2r_(m,k,j,i-1)+d2rc);
//      alpxp = fmax(0.0,fmax(lam1r_(m,k,j-1,i),lam1r_(m,k,j,i)));
//      alpxm = -fmin(0.0,fmin(lam1l_(m,k,j-1,i),lam1l_(m,k,j,i)));
//      alpyp = fmax(0.0,fmax(lam2r_(m,k,j,i-1),lam2r_(m,k,j,i)));
//      alpym = -fmin(0.0,fmin(lam2l_(m,k,j,i-1),lam2l_(m,k,j,i)));
//      axW = alpxp/(alpxp+alpxm);
//      axE = 1.0 - axW;
//      ayS = alpyp/(alpyp+alpym);
//      ayN = 1.0 - ayS;
//      dxW = alpxp*alpxm/(alpxp+alpxm);
//      dxE = dxW;
//      dyS = alpyp*alpym/(alpyp+alpym);
//      dyN = dyS;
//
//      vyN = axW*v21l_(m,k,j,i) + axE*v21r_(m,k,j,i);
//      vyS = axW*v21l_(m,k,j,i-1) + axE*v21r_(m,k,j,i-1);
//      vxE = ayS*v12l_(m,k,j,i) + ayN*v12r_(m,k,j,i);
//      vxW = ayS*v12l_(m,k,j-1,i) + ayN*v12r_(m,k,j-1,i);
//      ByE = b2c;
//      ByW = b2(m,k,j,i-1);
//      BxN = b1c;
//      BxS = b1(m,k,j-1,i);
//      e3(m,k,j,i) = - (axW*vxW*ByW + axE*vxE*ByE) + (ayN*vyN*BxN + ayS*vyS*BxS) + (dxE*ByE - dxW*ByW) - (dyN*BxN - dyS*BxS);
//
//        // integrate E1 to corner using MZ21
////        alm = 0.5*(1.0+nus2_(m,k-1,j,i));
////        arm = 1.0 - alm;
////        al = 0.5*(1.0+nus2c);
////        ar = 1.0 - al;
////        axW = 0.5*(alm+al);
////        axE = 0.5*(arm+ar);
////        alm = 0.5*(1.0+nus3_(m,k,j-1,i));
////        arm = 1.0 - alm;
////        al = 0.5*(1.0+nus3c);
////        ar = 1.0 - al;
////        ayS = 0.5*(alm+al);
////        ayN = 0.5*(arm+ar);
////        dxW = 0.5*(d2l_(m,k-1,j,i)+d2lc);
////        dxE = 0.5*(d2r_(m,k-1,j,i)+d2rc);
////        dyS = 0.5*(d3l_(m,k,j-1,i)+d3lc);
////        dyN = 0.5*(d3r_(m,k,j-1,i)+d3rc);
//        alpxp = fmax(0.0,fmax(lam2r_(m,k-1,j,i),lam2r_(m,k,j,i)));
//        alpxm = -fmin(0.0,fmin(lam2l_(m,k-1,j,i),lam2l_(m,k,j,i)));
//        alpyp = fmax(0.0,fmax(lam3r_(m,k,j-1,i),lam3r_(m,k,j,i)));
//        alpym = -fmin(0.0,fmin(lam3l_(m,k,j-1,i),lam3l_(m,k,j,i)));
//        axW = alpxp/(alpxp+alpxm);
//        axE = 1.0 - axW;
//        ayS = alpyp/(alpyp+alpym);
//        ayN = 1.0 - ayS;
//        dxW = alpxp*alpxm/(alpxp+alpxm);
//        dxE = dxW;
//        dyS = alpyp*alpym/(alpyp+alpym);
//        dyN = dyS;
//
//        vyN = axW*v32l_(m,k,j,i) + axE*v32r_(m,k,j,i);
//        vyS = axW*v32l_(m,k,j-1,i) + axE*v32r_(m,k,j-1,i);
//        vxE = ayS*v23l_(m,k,j,i) + ayN*v23r_(m,k,j,i);
//        vxW = ayS*v23l_(m,k-1,j,i) + ayN*v23r_(m,k-1,j,i);
//        ByE = b3c;
//        ByW = b3(m,k,j-1,i);
//        BxN = b2c;
//        BxS = b2(m,k-1,j,i);
//        e1(m,k,j,i) = - (axW*vxW*ByW + axE*vxE*ByE) + (ayN*vyN*BxN + ayS*vyS*BxS) + (dxE*ByE - dxW*ByW) - (dyN*BxN - dyS*BxS);
//
//        // integrate E2 to corner using MZ21
////        alm = 0.5*(1.0+nus3_(m,k,j,i-1));
////        arm = 1.0 - alm;
////        al = 0.5*(1.0+nus3c);
////        ar = 1.0 - al;
////        axW = 0.5*(alm+al);
////        axE = 0.5*(arm+ar);
////        alm = 0.5*(1.0+nus1_(m,k-1,j,i));
////        arm = 1.0 - alm;
////        al = 0.5*(1.0+nus1c);
////        ar = 1.0 - al;
////        ayS = 0.5*(alm+al);
////        ayN = 0.5*(arm+ar);
////        dxW = 0.5*(d3l_(m,k,j,i-1)+d3lc);
////        dxE = 0.5*(d3r_(m,k,j,i-1)+d3rc);
////        dyS = 0.5*(d1l_(m,k-1,j,i)+d1lc);
////        dyN = 0.5*(d1r_(m,k-1,j,i)+d1rc);
//        alpxp = fmax(0.0,fmax(lam3r_(m,k,j,i-1),lam3r_(m,k,j,i)));
//        alpxm = -fmin(0.0,fmin(lam3l_(m,k,j,i-1),lam3l_(m,k,j,i)));
//        alpyp = fmax(0.0,fmax(lam1r_(m,k-1,j,i),lam1r_(m,k,j,i)));
//        alpym = -fmin(0.0,fmin(lam1l_(m,k-1,j,i),lam1l_(m,k,j,i)));
//        axW = alpxp/(alpxp+alpxm);
//        axE = 1.0 - axW;
//        ayS = alpyp/(alpyp+alpym);
//        ayN = 1.0 - ayS;
//        dxW = alpxp*alpxm/(alpxp+alpxm);
//        dxE = dxW;
//        dyS = alpyp*alpym/(alpyp+alpym);
//        dyN = dyS;
//
//        vyN = axW*v13l_(m,k,j,i) + axE*v13r_(m,k,j,i);
//        vyS = axW*v13l_(m,k-1,j,i) + axE*v13r_(m,k-1,j,i);
//        vxE = ayS*v31l_(m,k,j,i) + ayN*v31r_(m,k,j,i);
//        vxW = ayS*v31l_(m,k,j,i-1) + ayN*v31r_(m,k,j,i-1);
//        ByE = b1c;
//        ByW = b1(m,k-1,j,i);
//        BxN = b3c;
//        BxS = b3(m,k,j,i-1);
//        e2(m,k,j,i) = - (axW*vxW*ByW + axE*vxE*ByE) + (ayN*vyN*BxN + ayS*vyS*BxS) + (dxE*ByE - dxW*ByW) - (dyN*BxN - dyS*BxS);

    });
  }
  if (pmy_pack->pmesh->use_polar_boundary) PolarAzimuthalAverageErUCT();
  return TaskStatus::complete;
}

void MHD::PolarAzimuthalAverageErUCT() {
  auto e1 = efld.x1e;

  auto &indcs = pmy_pack->pmesh->mb_indcs;

  int nx1 = indcs.nx1;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  const int ncells1 = nx1 + 2*ng;

  const int mesh_nx3 = pmy_pack->pmesh->mesh_indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;

  auto &mb_bcs = pmy_pack->pmb->mb_bcs;

  // Same treatment as PolarAzimuthalAverageEr in mhd_corner_e.cpp -- see that function for
  // why.  The two share their work arrays: only one of the two corner-E paths runs.
  // The host mirrors must be in the SAME guard condition as the device arrays and are
  // tested separately, because the MHD constructor pre-sizes inner_local/outer_local
  // (mhd.cpp) but not these two.  A guard on inner_local alone therefore never fires and
  // leaves the mirrors at their default extent 0.  Nothing notices in serial -- only the
  // MPI branch below touches them -- so every multi-rank run with polar boundaries died
  // on the first cycle with "deep_copy extents of views don't match: (0) inner(n)".
  if (static_cast<int>(inner_local.extent(0)) != ncells1 ||
      static_cast<int>(polar_inner_h.extent(0)) != ncells1) {
    Kokkos::realloc(inner_local, ncells1);
    Kokkos::realloc(outer_local, ncells1);
    Kokkos::realloc(polar_inner_h, ncells1);
    Kokkos::realloc(polar_outer_h, ncells1);
  }
  if (static_cast<int>(polar_part_in.extent(0)) != nmb ||
      static_cast<int>(polar_part_in.extent(1)) != ncells1) {
    Kokkos::realloc(polar_part_in,  nmb, ncells1);
    Kokkos::realloc(polar_part_out, nmb, ncells1);
  }
  auto part_in  = polar_part_in;
  auto part_out = polar_part_out;
  auto inner_   = inner_local;
  auto outer_   = outer_local;

  // one thread per (MeshBlock, radius); phi stays serial and ascending within a block
  par_for("polar_partial_uct", DevExeSpace(), 0, nmb-1, 0, ncells1-1,
  KOKKOS_LAMBDA(const int m, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    if ((do_inner || do_outer) && (i >= is) && (i <= ie)) {
      for (int k=ks+1; k<=ke+1; ++k) {
        if (do_inner) inner_sum += e1(m,k,js,i);
        if (do_outer) outer_sum += e1(m,k,je+1,i);
      }
    }
    part_in(m,i)  = inner_sum;
    part_out(m,i) = outer_sum;
  });

  // combine per-block partials in ascending block order
  par_for("polar_combine_uct", DevExeSpace(), 0, ncells1-1,
  KOKKOS_LAMBDA(const int i) {
    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    for (int m=0; m<nmb; ++m) {
      inner_sum += part_in(m,i);
      outer_sum += part_out(m,i);
    }
    inner_(i) = inner_sum;
    outer_(i) = outer_sum;
  });

#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    Kokkos::deep_copy(polar_inner_h, inner_local);
    Kokkos::deep_copy(polar_outer_h, outer_local);
    MPI_Allreduce(MPI_IN_PLACE, polar_inner_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, polar_outer_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    Kokkos::deep_copy(inner_local, polar_inner_h);
    Kokkos::deep_copy(outer_local, polar_outer_h);
  }
#endif

  par_for("polar_refill_uct", DevExeSpace(), 0, nmb-1, ks, ke+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    if (do_inner) e1(m,k,js  ,i) = inner_(i)/static_cast<Real>(mesh_nx3);
    if (do_outer) e1(m,k,je+1,i) = outer_(i)/static_cast<Real>(mesh_nx3);
  });
}

} // namespace mhd
