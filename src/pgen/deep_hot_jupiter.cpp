//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter.cpp
//! \brief Problem generator for the deep hot Jupiter.
//!
//! REFERENCE: Heng, Menou, Phillipps, MNRAS, 413, 2380 (2011); Deitrick, Mendonça, Schroffenegger, Grimm, Tsai, Heng, ApJS, 248, 30 (2020)

// C++ headers
#include <cmath>
#include <iostream> // cout

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"

#include <Kokkos_Random.hpp>

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);
KOKKOS_INLINE_FUNCTION
void StretchR(const Real r0, const Real r1, Real &r) {
  Real xi = (r-r0)/(r1-r0);
  Real a = 2.0;
  if (fabs(a) < 1e-12) {
    r = r0 + (r1 - r0)*xi;
  }
  Real denom = 1.0 - exp(-a);
  r = r0 + (r1 - r0)*(1.0 - exp(-a*xi))/denom;
};
KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_itrad(const Real &p, Real &itrad);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_eq_Tp(const Real &lam, const Real &phi, const Real &p, Real &Teq);
template <typename View1D>
void get_wb_eos_arr(const int &N, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
//KOKKOS_INLINE_FUNCTION
//void get_init_eos_arr(const Real &lam, const Real &phi, const int &N, DvceArray1D<Real> zarr, DvceArray1D<Real> logparr);
//KOKKOS_INLINE_FUNCTION
//void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, const Real &lam, const Real &phi, Real &rho, Real &p);


//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  const bool use_etotgrav = pmy_mesh_->pmb_pack->phydro->use_etotgrav;
  const bool use_wellbalance = pmy_mesh_->pmb_pack->phydro->use_wellbalance;
  const bool use_wellbalance_local = pmy_mesh_->pmb_pack->phydro->use_wellbalance_local;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  const bool use_grid_stretch = pmy_mesh_->use_grid_stretch;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  user_bcs_func = HydrostaticEquilibrium;
  if (restart) return;
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
    
    Real r0, r1;
    if (use_grid_stretch) {
        r0 = pmy_mesh_->mesh_size.x1min;
        r1 = pmy_mesh_->mesh_size.x1max;
    }
    
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_;
    
    auto phi0_x1f = pmbp->phydro->phi0.x1f;
    auto phi0_x2f = pmbp->phydro->phi0.x2f;
    auto phi0_x3f = pmbp->phydro->phi0.x3f;

    auto phicc0 = pmbp->phydro->phicc0;
    
    auto u0wb = pmbp->phydro->u0wb;
    auto w0wb = pmbp->phydro->w0wb;
    auto w0facewb_x1f = pmbp->phydro->w0facewb.x1f;
    auto w0facewb_x2f = pmbp->phydro->w0facewb.x2f;
    auto w0facewb_x3f = pmbp->phydro->w0facewb.x3f;

  Real gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
  Real grav_acc = -942.0;
  Real ap = 9.44e9;
  Real iap = 1.0/ap;
  Real omega = 2.06e-5;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_wb_eos_arr(N, zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    
//    DvceArray1D<Real> zinitarr("zinitarr", N);
//    DvceArray1D<Real> logpinitarr("logpinitarr", N);
    DualArray1D<Real> zarr_init("zarrinit", N);
    DualArray1D<Real> logparr_init("logparrinit", N);
    get_init_eos_arr(N, zarr_init.h_view, logparr_init.h_view);
    zarr_init.modify_host();
    zarr_init.sync_device();
    logparr_init.modify_host();
    logparr_init.sync_device();
  
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      if (use_grid_stretch) StretchR(r0,r1,x1v);
      Real r = x1v;
      if (use_spherical_polar) x1v -= ap;

      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        
      Real lam, theta, phi;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -x2v+M_PI/2.0;
        phi = x3v-M_PI;
      } else {
        theta = -(x2v*iap-M_PI/2.0);
        lam = x2v*iap;
        phi = x1v*iap;
      }
        
      Real pwb, denwb;
      get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
      Real p, den;
      get_init_eos(zarr_init.d_view,logparr_init.d_view,x1v,den,p);
//      get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//      get_init_eos(zinitarr,logpinitarr,x3v,lam,phi,den,p);
//      p = pwb;
//      den = denwb;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = p*igm1;
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
        if (use_etotgrav) {
            u0_(m,IEN,k,j,i) += den*phicc;
        }
    });
            
    par_for("probwb", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,x1v);
        Real r = x1v;
        if (use_spherical_polar) x1v -= ap;

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
          
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          
        Real lam, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -x2v+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          theta = -(x3v*iap-M_PI/2.0);
          lam = x3v*iap;
          phi = x2v*iap;
        }
          
        Real pwb, denwb;
        get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
        Real p, den;
        get_init_eos(zarr_init.d_view,logparr_init.d_view,x1v,den,p);
//        get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//        get_init_eos(zinitarr,logpinitarr,x1v,lam,phi,den,p);
//        p = pwb;
//        den = denwb;
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
        
      if (use_etotgrav || use_wellbalance_local) {
          x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              if (use_grid_stretch) StretchR(r0,r1,x1v);
              r = x1v;
              if (use_spherical_polar) x1v -= ap;
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x2f(m,k,j,i) = phicc;
          if (j == je) {
              x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x3f(m,k+1,j,i) = phicc;
          }
      }
        if (use_wellbalance) {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            if (use_grid_stretch) StretchR(r0,r1,x1v);
            if (use_spherical_polar) x1v -= ap;
            Real denwb;
            get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = pwb*igm1;
            if (i == ie) {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                if (use_grid_stretch) StretchR(r0,r1,x1v);
                if (use_spherical_polar) x1v -= ap;
                get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = pwb*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            if (use_grid_stretch) StretchR(r0,r1,x1v);
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = pwb*igm1;
            if (j == je) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                if (use_grid_stretch) StretchR(r0,r1,x1v);
                if (use_spherical_polar) x1v -= ap;
                get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = pwb*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            if (use_grid_stretch) StretchR(r0,r1,x1v);
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = pwb*igm1;
            if (k == ke) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
                if (use_grid_stretch) StretchR(r0,r1,x1v);
                if (use_spherical_polar) x1v -= ap;
                get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = pwb*igm1;
            }
        }
    });
    if (use_etotgrav || use_wellbalance_local) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
            if (use_grid_stretch) StretchR(r0,r1,x1v);
            Real r = x1v;
            if (use_spherical_polar) x1v -= ap;

            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
              
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
              
            Real lam, theta, phi;
            if (use_spherical_polar) {
              theta = x2v;
              lam = -x2v+M_PI/2.0;
              phi = x3v-M_PI;
            } else {
              theta = -(x3v*iap-M_PI/2.0);
              lam = x3v*iap;
              phi = x2v*iap;
            }
            
            Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
            phicc0(m,k,j,i) = phicc;
        });
    }
    if (use_wellbalance) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          int nx1 = indcs.nx1;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
          if (use_spherical_polar) x1v -= ap;

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
            
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            
          Real pwb, denwb;
          get_wb_eos(zarr.d_view,logparr.d_view,x1v,denwb,pwb);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = pwb*igm1;
          if (use_etotgrav) {
              Real phicc = - grav_acc * x1v;
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = pwb*igm1;
        });
    }

    // initialize magnetic fields if MHD
    if (pmbp->pmhd != nullptr) {
      // Read magnetic field strength
      Real bx = pin->GetReal("problem","b0");
      auto &b0 = pmbp->pmhd->b0;
      par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        b0.x1f(m,k,j,i) = bx;
        b0.x2f(m,k,j,i) = 0.0;
        b0.x3f(m,k,j,i) = 0.0;
        if (i==ie) b0.x1f(m,k,j,i+1) = bx;
        if (j==je) b0.x2f(m,k,j+1,i) = 0.0;
        if (k==ke) b0.x3f(m,k+1,j,i) = 0.0;
        u0_(m,IEN,k,j,i) += 0.5*bx*bx;
      });
    }

  return;
}


void HydrostaticEquilibrium(Mesh *pm) {
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
    const bool use_wellbalance_local = pmbp->phydro->use_wellbalance_local;
    
    auto phi0_x1f = pmbp->phydro->phi0.x1f;
    auto phicc0 = pmbp->phydro->phicc0;

    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->pmhd->w0;
    }
    
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real igm1 = 1.0/(gamma-1.0);
    Real gigm1 = gamma*igm1;
    Real gm1ig = (gamma-1.0)/gamma;
    Real ig = 1.0/gamma;
    
    par_for("usrboundaryx1", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
          for (int i=0; i<ng; ++i) {
              if (use_wellbalance_local) {
                Real rho_i = u0_(m,IDN,k,j,ie);
                Real e_i = u0_(m,IEN,k,j,ie) - 0.5*(SQR(u0_(m,IM1,k,j,ie))+SQR(u0_(m,IM2,k,j,ie))+SQR(u0_(m,IM3,k,j,ie)))/rho_i;
                if (use_etotgrav) e_i -= rho_i*phicc0(m,k,j,ie);
                Real phi_i = phicc0(m,k,j,ie);
//                Real q0_i = pow(e_i,gm1ig);
//                Real factor_i = pow(rho_i,gamma)/e_i*ig;
                Real q0_i = log(e_i);
                Real factor_i = rho_i/e_i*igm1;
                Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
                Real q0_ip = q0_i - factor_i * dphi_i;
//                Real e0_ip = pow(q0_ip,gigm1);
//                Real rho0_ip = pow(e0_ip/e_i,ig)*rho_i;
                Real e0_ip = exp(q0_ip);
                Real rho0_ip = e0_ip/e_i*rho_i;
                u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
                u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
                u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
                Real mom = fmax(0.0,u0_(m,IM1,k,j,ie)/rho_i*rho0_ip);
//                Real mom = u0_(m,IM3,ke,j,i)/rho_i*rho0_ip;
                u0_(m,IM1,k,j,(ie+i+1)) = mom;
                u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(ie+i+1)))+SQR(u0_(m,IM2,k,j,(ie+i+1)))+SQR(u0_(m,IM3,k,j,(ie+i+1))))/rho0_ip;
                if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
              } else {
            u0_(m,IDN,k,j,(ie+i+1)) = u0_(m,IDN,k,j,ie);
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie);
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie);
            Real mom = fmax(0.0,u0_(m,IM1,k,j,ie));
            u0_(m,IEN,k,j,(ie+i+1)) = u0_(m,IEN,k,j,ie) - 0.5*(SQR(u0_(m,IM1,k,j,ie))-SQR(mom))/u0_(m,IDN,k,j,ie);
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
              }
          }
        }
    });
  return;
}


void SourceFunc(Mesh *pm, Real bdt) {
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0, w0wb;
    u0 = pmbp->phydro->u0;
    w0 = pmbp->phydro->w0;
    w0wb = pmbp->phydro->w0wb;
    auto phi0_x1f = pmbp->phydro->phi0.x1f;
    auto phicc0 = pmbp->phydro->phicc0;
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
    const bool use_wellbalance = pmbp->phydro->use_wellbalance;
    const bool use_spherical_polar = pm->use_spherical_polar;
    const bool use_grid_stretch = pm->use_grid_stretch;
    const bool use_wellbalance_local = pmbp->phydro->use_wellbalance_local;
    
    Real r0, r1;
    if (use_grid_stretch) {
        r0 = pm->mesh_size.x1min;
        r1 = pm->mesh_size.x1max;
    }
    
    Real ap = 9.44e9;
    Real iap = 1.0/ap;
    Real omega = 2.06e-5;
    Real grav_acc = -942.0;
    Real Rgas = 4.593e7;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    
    Real time = pm->time;

    par_for("usrsource", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
        
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        if (use_grid_stretch) StretchR(r0,r1,x1v);
        
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
          z = x1v-ap;
          r = x1v;
        } else {
          lam = x3v*iap;
          phi = x2v*iap;
          z = x1v;
        }
        Real rho = w0(m,IDN,k,j,i);
        Real p = w0(m,IEN,k,j,i)*gm1;
        Real T = p/Rgas/rho;
        
        Real area_r = area1(m,k,j,i+1);
        Real area_l = area1(m,k,j,i);
        Real vol = volume(m,k,j,i);
        
        // gravity
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_local) {
//          Real dphi_i = phi0_x1f(m,k,j,i+1)-phicc0(m,k,j,i);
//          Real dphi_imh = phicc0(m,k,j,i)-phi0_x1f(m,k,j,i);
//          src = -bdt*(area_r*dphi_i+area_l*dphi_imh)/vol*w0(m,IDN,k,j,i);
          Real e_imh,e_iph,dum1,dum2,dum3;
          pmbp->phydro->getWBerho(IEN, gamma,
            w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
            phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
            dum1,e_imh,dum2,e_iph,dum3);
          Real pl = e_imh*gm1;
          Real pr = e_iph*gm1;
          src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol;
        }
        u0(m,IM1,k,j,i) += src;
        
        // Forces in the corotating frame
        if (use_spherical_polar) {
          Real vtheta = w0(m,IVY,k,j,i);
          Real vphi = w0(m,IVZ,k,j,i);
          Real vr = w0(m,IVX,k,j,i);
          Real sine = sin(theta);
          Real cosine = cos(theta);
          Real oor = SQR(omega)*r*sine;
          Real cor = 2.0*omega*vphi;
          u0(m,IM2,k,j,i) += rho*(cor+oor)*cosine*bdt;
          u0(m,IM3,k,j,i) += -rho*2.0*omega*(vr*sine+vtheta*cosine)*bdt;
          u0(m,IM1,k,j,i) += rho*(cor+oor)*sine*bdt;
//          if (!use_etotgrav)
          u0(m,IEN,k,j,i) += rho*oor*(vr*sine+vtheta*cosine)*bdt;
        } else {
          // corotating beta-plane approximation e.g. Fromang+2016
          Real omega1 = omega*lam;
          u0(m,IM2,k,j,i) += -2.0*rho*omega1*(-w0(m,IVZ,k,j,i))*bdt;
          u0(m,IM3,k,j,i) += -2.0*rho*omega1*w0(m,IVY,k,j,i)*bdt;
        }

        // Newtonian cooling
        Real Teq, itrad;
        get_eq_Tp(lam, phi, p, Teq);
        get_itrad(p,itrad);
//        Real t0 = 2.16e6;
//        Real ff = (t0-time)/t0;
//        ff = (ff < 0.0) ? 0.0 : ff;
//        Real gg = 2.0*ff;
//        itrad /= pow(10.0,gg);
        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;
        
        // Rayleigh drag (initial relaxation)
        Real dyntime = 2.0*M_PI/omega;
        Real tau1 = dyntime / 10.0;
        Real tau2 = dyntime;
        Real t1   = 2.0 * dyntime;
        Real t2   = 5.0 * dyntime;
        Real itdrag;
        if(time < t1) {
          itdrag = 1.0 / tau1;
        } else if(time < t2)
          {
            Real alp = (time - t1) / (t2 - t1);
            itdrag   = 1.0 / tau1 * pow(tau1 / tau2, alp);
          }
        else {
          itdrag = 0.0;
        }
        if (time < t2) {
          Real fredux = itdrag*bdt;///(1.0+itdrag*bdt);
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }
        
////        // Top sponge layer
////        Real bar = 1.0e6;
////        Real logpl = log(1.0e-4*bar);
////        Real logpt = log(1.0e-5*bar);
////        Real logp = log(p);
////        Real fdrag = 1.0 - (logp-logpt)/(logpl-logpt); // high p = 0, low p = 1
////        fdrag = fmax(fdrag,0.0);
////        fdrag = fmin(fdrag,1.0);
////        Real itdrag = fdrag/1.0e3;
////        Real fredux = itdrag*bdt/(1.0+itdrag*bdt);
////        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
////        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
////        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
    });

    return;
}


KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td) {
    
    Real bar = 1.0e6;
    Real pl = 1.0e-3*bar;
    Real pt = log10(p/bar);
    Real ptl = log10(pl/bar);
    
    Real fn[13];
    Real fd[14];
    
    fn[0] = 1388.77348;
    fn[1] = 279.575848;
    fn[2] = -213.835822;
    fn[3] = 21.0010475;
    fn[4] = 100.938036;
    fn[5] = 12.7972336;
    fn[6] = -13.9266925;
    fn[7] = -3.70783272;
    fn[8] = 0.522370269;
    fn[9] = 0.320837882;
    fn[10]= 0.0451831612;
    fn[11]= 2.18195583e-3;
    fn[12]= 3.98938097e-6;

    fd[0] = 2152.06036;
    fd[1] = 29.3485512;
    fd[2] = -183.318696;
    fd[3] = 46.3893130;
    fd[4] = 19.8116485;
    fd[5] = -28.5473177;
    fd[6] = -2.52726545;
    fd[7] = 8.43627538;
    fd[8] = 2.62945375;
    fd[9] = -0.297098168;
    fd[10]= -0.286871487;
    fd[11]= -0.0590629443;
    fd[12]= -5.38679474e-3;
    fd[13]= -1.89972415e-4;
    
    Real Tnstar = 0.0;
    Real Tdstar = 0.0;
    Real Tnstar_pl = 0.0;
    Real Tdstar_pl = 0.0;
    Real ptn = 1.0;
    Real ptln = 1.0;
    for(int ilogp=0; ilogp<13; ilogp++) {
        Tnstar += fn[ilogp]*ptn;
        Tnstar_pl += fn[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }
    ptn = 1.0;
    ptln = 1.0;
    for(int ilogp=0; ilogp<14; ilogp++) {
        Tdstar += fd[ilogp]*ptn;
        Tdstar_pl += fd[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }

    Tn = Tnstar;
    Td = Tdstar;
//    if (p < pl) {
//        Tn = Tnstar_pl * exp(0.1*log10(p/pl));
//        if (Tn < 250) {
//            Tn = 250;
//        }
//        Td = Tdstar_pl * exp(0.015*log10(p/pl));
//        if (Td < 1000) {
//            Td = 1000;
//        }
//    }
//    if (Td < Tn) {
//        Td = Tn;
//    }
    Real x = log10(p / pl);

    Real Tn_new = Tnstar_pl;// * exp(0.1   * x);
    Real Td_new = Tdstar_pl;// * exp(0.015 * x);

    Tn_new = fmax(Tn_new, 250.0);
    Td_new = fmax(Td_new, 1000.0);

    // blend
    Tn = (p < pl) ? Tn_new : Tn;
    Td = (p < pl) ? Td_new : Td;

    // enforce Td >= Tn
    Td = fmax(Td, Tn);

    return;
}

KOKKOS_INLINE_FUNCTION
void get_itrad(const Real &p, Real &itrad) {
    
    Real bar = 1.0e6;
    Real pl = 1.0e-5*bar;
    Real pt = log10(p/bar);
    Real ptl = log10(pl/bar);
    
    Real ft[5];
    
    ft[0] = 5.4659686;
    ft[1] = 1.4940124;
    ft[2] = 0.66079196;
    ft[3] = 0.16475329;
    ft[4] = 0.014241552;
    
    Real logtrad = 0.0;
    Real logtrad_pl = 0.0;
    Real ptn = 1.0;
    Real ptln = 1.0;
    for(int it=0; it<5; it++) {
        logtrad += ft[it]*ptn;
        logtrad_pl += ft[it]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }
//    if (p < pl) {
//        logtrad = logtrad_pl;
//    }
//    itrad = 1.0 / pow(logtrad,10.0);
//
//    if (p >= 10.0*bar) {
//        itrad = 0.0;
//    }
    logtrad = (p < pl) ? logtrad_pl : logtrad;

    Real trad = pow(10.0, logtrad);

    itrad = 1.0 / trad;

    itrad *= (p < 10.0 * bar);

    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
    T = Tn;
    
//    Real Tn2 = Tn*Tn;
//    Real Tn4 = Tn2*Tn2;
//    Real Td2 = Td*Td;
//    Real Td4 = Td2*Td2;
//    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;
//
//    Real Tmid = sqrt(sqrt(Tmid4));
//    T = Tmid;
//
//    Real bar = 1.0e6;
//    Real ph = 1.0e-2*bar;
//    Real pl = 0.1*ph;
//    Real logph = log(ph);
//    Real logpl = log(pl);
//    Real logp = log(p);
//    Real btrop = (logp < logpl)? 0.0 : ((logp-logpl)/(logph-logpl));
//    btrop = (logp < logph)? btrop : 1.0;
//    T = btrop*Tmid + (1.0-btrop)*Tn;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_Tp(const Real &p, Real &T) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
//    T = Tn;
    
    Real Tn2 = Tn*Tn;
    Real Tn4 = Tn2*Tn2;
    Real Td2 = Td*Td;
    Real Td4 = Td2*Td2;
    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;

    Real Tmid = sqrt(sqrt(Tmid4));
    T = Tmid;

//    Real bar = 1.0e6;
//    Real ph = 1.0e-2*bar;
//    Real pl = 0.1*ph;
//    Real logph = log(ph);
//    Real logpl = log(pl);
//    Real logp = log(p);
//    Real btrop = (logp < logpl)? 0.0 : ((logp-logpl)/(logph-logpl));
//    btrop = (logp < logph)? btrop : 1.0;
//    T = btrop*Tmid + (1.0-btrop)*Tn;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_eq_Tp(const Real &lam, const Real &phi, const Real &p, Real &Teq) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
    Real Tn2 = Tn*Tn;
    Real Tn4 = Tn2*Tn2;
    Real Td2 = Td*Td;
    Real Td4 = Td2*Td2;
    Real Teq4 = Tn4 + (Td4-Tn4) * cos(lam) * cos(phi);
    
    Teq = sqrt(sqrt(Teq4));
    Teq = (phi > -M_PI/2.0 && phi < M_PI/2.0) ? Teq : Tn;
    
    return;
}

template <typename View1D>
void get_wb_eos_arr(const int &N, View1D zarr, View1D logparr) {
    
    Real Rgas = 4.593e7;
    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 220.0*bar;

    Real zmin = 0.0;
    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_wb_Tp(p,T);

        logparr(n+1) = logparr(n) + fac/T;
    }

    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 220.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_wb_Tp(p,T);
        rho = p/Rgas/T;
    } else {
        Real T0;
        get_wb_Tp(p0,T0);
        Real rho0 = p0/Rgas/T0;
        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr) {

    Real Rgas = 4.593e7;
    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 220.0*bar;

    Real zmin = 0.0;
    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_init_Tp(p,T);

        logparr(n+1) = logparr(n) + fac/T;
    }


    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {

    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 220.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tp(p,T);
        rho = p/Rgas/T;
    } else {
        Real T0;
        get_init_Tp(p0,T0);
        Real rho0 = p0/Rgas/T0;
        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

//KOKKOS_INLINE_FUNCTION
//void get_init_eos_arr(const Real &lam, const Real &phi, const int &N, DvceArray1D<Real> zarr, DvceArray1D<Real> logparr) {
//
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
//    Real bar = 1.0e6;
//    Real p0 = 220.0*bar;
//
//    Real zmin = 0.0;
//    Real zmax = 1.2e9;
//    Real dz = (zmax-zmin)/N;
//    Real fac = grav_acc/Rgas*dz;
//    logparr(0) = std::log(p0);
//
//    Real ph = 1.0e-2*bar;
//    Real pl = 0.1*ph;
//    Real logph = log(ph);
//    Real logpl = log(pl);
//
//    for(int n=0; n<N; n++) {
//        Real Thse, Teq;
//        Real logp = logparr(n);
//        Real p = exp(logp);
//        zarr(n) = zmin + n*dz;
//
////        get_wb_Tp(p,Thse);
//        get_eq_Tp(lam,phi,p,Teq);
////        Real btrop = (logp < logpl)? 0.0 : ((logp-logpl)/(logph-logpl));
////        btrop = (logp < logph)? btrop : 1.0;
////        Real T = btrop*Thse + (1.0-btrop)*Teq;
//        Real T = Teq;
//
//        logparr(n+1) = logparr(n) + fac/T;
//    }
//
//    return;
//}
//
//KOKKOS_INLINE_FUNCTION
//void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, const Real &lam, const Real &phi, Real &rho, Real &p) {
//
//    Real Rgas = 4.593e7;
//    Real bar = 1.0e6;
//    Real p0 = 220.0*bar;
//    Real dz = zarr(1)-zarr(0);
//    Real T, Thse, Teq;
//
//    Real ph = 1.0e-2*bar;
//    Real pl = 0.1*ph;
//    Real logph = log(ph);
//    Real logpl = log(pl);
//
//    if (z >= 0.0) {
//        int Nt = std::floor(z/dz);
//        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
//        p = std::exp(logp);
////        get_wb_Tp(p,Thse);
//        get_eq_Tp(lam,phi,p,Teq);
////        Real btrop = (logp < logpl)? 0.0 : ((logp-logpl)/(logph-logpl));
////        btrop = (logp < logph)? btrop : 1.0;
////        Real T = btrop*Thse + (1.0-btrop)*Teq;
//        Real T = Teq;
//        rho = p/Rgas/T;
//    } else {
//        Real T0;
//        get_wb_Tp(p0,T0);
//        Real rho0 = p0/Rgas/T0;
//        Real grav_acc = -942.0;
//        Real H0 = -p0/rho0/grav_acc;
//        Real iH0 = 1.0/H0;
//        p = p0 * std::exp(-z*iH0);
//        rho = rho0 * std::exp(-z*iH0);
//    }
//
//    return;
//}
