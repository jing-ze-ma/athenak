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

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt);
void double_gray_two_stream_RT(Mesh *pm, Real bdt);

void picket_fence_two_stream_RT(Mesh *pm, Real bdt);
KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2);
KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr);
KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real &taulim, const Real &A, const Real &B, const Real (&C)[3], const Real (&D)[3], const Real (&E)[3], const Real (&gamv)[3], const Real &tau, Real &T);
template <typename View1D>
void get_picket_fence_pT_arr(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void adjust_ad_pT_arr(const int &N, View1D Tarr, View1D lgparr);

KOKKOS_INLINE_FUNCTION
void StretchR(const Real r0, const Real r1, Real &r) {
  Real xi = (r-r0)/(r1-r0);
  Real a = 2.0;
  if (fabs(a) < 1.0e-12) {
    r = r0 + (r1 - r0)*xi;
  }
  Real denom = 1.0 - exp(-a);
  r = r0 + (r1 - r0)*(1.0 - exp(-a*xi))/denom;
};
KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_wb_eos_arr(const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);


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
//    if (use_grid_stretch) {
    r0 = pmy_mesh_->mesh_size.x1min;
    r1 = pmy_mesh_->mesh_size.x1max;
//    }
    
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
    
//  Real Teq = 1469.0;
//  Real grav_acc = -942.0;
//  Real ap = 9.44e9;
//  Real omega = 2.06e-5;
//  Real Rgas = 4.593e7;
//  Real met = 0.0;
    
    Real Teq = pin->GetReal("problem","Teq");
    Real grav_acc = -pin->GetReal("problem","grav");
    Real ap = pin->GetReal("problem","ap");
    Real omega = pin->GetReal("problem","omega");
    Real Rgas = pin->GetReal("problem","Rgas");
    Real met = pin->GetReal("problem","met");
  
    Real iap = 1.0/ap;
    
    Real grav = -grav_acc;
    Real Tirr = Teq*sqrt(2);
    Real Tint;
    get_Tint(Teq, Tint);
    Real mus = 1.0;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_wb_eos_arr(Rgas, grav_acc, N, (r1-r0)*1.1, zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    
    DualArray1D<Real> Tarr_init("Tarrinit", N);
    DualArray1D<Real> lgparr_init("lgparrinit", N);
    get_picket_fence_pT_arr(Tint, Tirr, met, grav, mus, N, Tarr_init.h_view, lgparr_init.h_view);
    
    Tarr_init.modify_host();
    Tarr_init.sync_device();
    lgparr_init.modify_host();
    lgparr_init.sync_device();
    
//    DvceArray1D<Real> zinitarr("zinitarr", N);
//    DvceArray1D<Real> logpinitarr("logpinitarr", N);
    DualArray1D<Real> zarr_init("zarrinit", N);
    DualArray1D<Real> logparr_init("logparrinit", N);
    get_init_eos_arr(Rgas, grav_acc, Tarr_init.h_view, lgparr_init.h_view, N, (r1-r0)*1.1, zarr_init.h_view, logparr_init.h_view);
    
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
      get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
      Real p, den;
      get_init_eos(Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
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
        get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
        Real p, den;
        get_init_eos(Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
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
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
                get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
                get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
                get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
          get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
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
            Real rho_i = u0_(m,IDN,k,j,ie);
            Real e_i = u0_(m,IEN,k,j,ie) - 0.5*(SQR(u0_(m,IM1,k,j,ie))+SQR(u0_(m,IM2,k,j,ie))+SQR(u0_(m,IM3,k,j,ie)))/rho_i;
            if (use_etotgrav) e_i -= rho_i*phicc0(m,k,j,ie);
            Real phi_i = phicc0(m,k,j,ie);
            Real q0_i = log(e_i);
            Real factor_i = rho_i/e_i*igm1;
            Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
            Real q0_ip = q0_i - factor_i * dphi_i;
            Real e0_ip = exp(q0_ip);
            Real rho0_ip = e0_ip/e_i*rho_i;
            u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
              Real mom = u0_(m,IM1,k,j,ie)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,ie)/rho_i*rho0_ip);
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(ie+i+1)))+SQR(u0_(m,IM2,k,j,(ie+i+1)))+SQR(u0_(m,IM3,k,j,(ie+i+1))))/rho0_ip;
            if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
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
//    if (use_grid_stretch) {
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
//    }
    
//    Real ap = 9.44e9;
//    Real omega = 2.06e-5;
//    Real grav_acc = -942.0;
//    Real Rgas = 4.593e7;
    
//    ParameterInput* pin;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    
    Real iap = 1.0/ap;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    
    Real time = pm->time;
    
    picket_fence_two_stream_RT(pm, bdt);

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

//        // Newtonian cooling
//        Real Teq, itrad;
//        get_eq_Tp(lam, phi, p, Teq);
//        get_itrad(p,itrad);
////        Real t0 = 2.16e6;
////        Real ff = (t0-time)/t0;
////        ff = (ff < 0.0) ? 0.0 : ff;
////        Real gg = 2.0*ff;
////        itrad /= pow(10.0,gg);
//        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
//        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;
        
        // Rayleigh drag (initial relaxation)
        Real dyntime = 2.0*M_PI/omega;
        Real tau1 = dyntime / 10.0;
        Real tau2 = dyntime;
        Real t1   = 2.0 * dyntime;
        Real t2   = 5.0 * dyntime;
        Real itdrag, fredux;
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
          fredux = itdrag*bdt;///(1.0+itdrag*bdt);
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }
        
        // Top sponge layer
        Real bar = 1.0e6;
        Real logpl = log(1.0e-6*bar);
        Real logpt = log(1.0e-7*bar);
        Real logp = log(p);
        Real fdrag = 1.0 - (logp-logpt)/(logpl-logpt); // high p = 0, low p = 1
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
    });

    return;
}

void double_gray_two_stream_RT(Mesh *pm, Real bdt) {
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

    DvceArray5D<Real> u0, w0;
    u0 = pmbp->phydro->u0;
    w0 = pmbp->phydro->w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    const bool use_grid_stretch = pm->use_grid_stretch;
    const bool correct_spherical = false;
    const bool test_oned = true;
    
    Real r0, r1;
//    if (use_grid_stretch) {
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
//    }
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 100.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 5;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
        constexpr int NN = 270;
        Real tau_ir_down_f[NN];
        Real F_v_down_f[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real B[NN];
        
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        
        Real rtop = LeftEdgeX(ie+1-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rtop);
        Real rbot = LeftEdgeX(is-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rbot);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = -1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = w0(m,IEN,k,j,ie+1)*gm1;
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = p/Rgas/rho;
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real pm1 = w0(m,IEN,k,j,ie)*gm1;
        Real pf = exp((log(p)+log(pm1))/2.0);
        Real tau_v_f = 0.0;//pf/(grav/kap_v);
        Real tau_ir_f = 0.0;//pf/(grav/kap_ir);
        
        tau_ir_down_f[ie+1] = tau_ir_f;
        F_v_down_f[ie+1] = Fstar*mu0;//*exp(-tau_v_f/muf);
        F_v_down_f[ie+1] = (mu0 > 0.0) ? F_v_down_f[ie+1] : 0.0;
        
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
          Real kap_ir = 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
            
          Real rf = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rf);
          Real rf1 = LeftEdgeX(i+1-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rf1);

//          Real muf = sqrt(1.0-SQR(ap/rf)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r1));
          Real muf = (mu0 < mucr) ? mucr : mu0;
          if (test_oned) muf = mu0;
          Real dtau_v = kap_v*rho*dr;
          Real dtau_ir = kap_ir*rho*dr;
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf); // Zhang+2023
          Real trans = exp(-dtau_v/muf);
          F_v_down_f[i] = F_v_down_f[i+1]*trans*fac;
          tau_ir_down_f[i] = tau_ir_down_f[i+1] + dtau_ir;
        }
        p = w0(m,IEN,k,j,is-1)*gm1;
        rho = w0(m,IDN,k,j,is-1);
        T = p/Rgas/rho;
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        I_ir_down_f[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rf = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rf);
          Real rf1 = LeftEdgeX(i+1-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rf1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf);
          Real dtau = tau_ir_down_f[i]-tau_ir_down_f[i+1];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i]+B[i+1])/2.0;
          I_ir_down_f[i] = (I_ir_down_f[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        I_ir_up_f[is] = Iint;
        // up-sweep
        for (int i=is+1; i<ie+2; ++i) {
          Real rf = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rf);
          Real rfm1 = LeftEdgeX(i-1-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rfm1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rfm1/rf);
          Real dtau = tau_ir_down_f[i-1]-tau_ir_down_f[i];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i-1]+B[i])/2.0;
          I_ir_up_f[i] = (I_ir_up_f[i-1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        // flux divergence
        for (int i=is; i<ie+1; ++i) {
          Real Ft = 2.0*M_PI*mu1*(I_ir_up_f[i+1]-I_ir_down_f[i+1])-F_v_down_f[i+1];
          Real Fb = 2.0*M_PI*mu1*(I_ir_up_f[i]-I_ir_down_f[i])-F_v_down_f[i];
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) src = -(Ft*area_t-Fb*area_b)/vol;
          u0(m,IEN,k,j,i) += src*bdt;
        }
        
    });
    
    return;
}

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt) {
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

    DvceArray5D<Real> u0, w0;
    u0 = pmbp->phydro->u0;
    w0 = pmbp->phydro->w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    const bool use_grid_stretch = pm->use_grid_stretch;
    const bool correct_spherical = false;
    const bool test_oned = false;
    
    Real r0, r1;
//    if (use_grid_stretch) {
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
//    }
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real grav = 942.0;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 500.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 4;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(0), n1);
        constexpr int NN = 270;
        Real tau_ir_down[NN];
        Real I_down[NN];
        Real B[NN];
        Real Q_v[NN];
        
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(50.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = w0(m,IEN,k,j,ie+1)*gm1;
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = p/Rgas/rho;
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real tau_v = 0.0;//p/(grav/kap_v);
        Real tau_ir = 0.0;//p/(grav/kap_ir);
//        if (test_oned) {
//          tau_v = p/(grav/kap_v)/mu0;
//          tau_ir = p/(grav/kap_ir)/mu1;
//        }
        tau_ir_down[ie+1] = tau_ir;
        for (int i=ie; i>is-1; --i) {
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
            Real kap_ir = 1.0e-2; // 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,x1v);
            
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
////          Real mu = sqrt(1.0-SQR(r1/r)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r));
//          Real mu = (mu0 < mucr) ? mucr : mu0;
//          if (test_oned) mu = mu0;
//          Real dtau_v = kap_v*rho*dr;
          Real delta = dr/r;
          Real drcor = r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          if (test_oned) drcor = dr/mu0;
          Real dtau_v = kap_v*rho*drcor;
          Real dtau_ir = kap_ir*rho*dr;
            
          tau_v += dtau_v;
          Real fac = 1.0;
//          if (correct_spherical) fac = SQR(r1/r); // Zhang+2023
//          Real Q_v = kap_v*rho*Fstar*fac*exp(-tau_v/mu); // Zhang+2023
//          Q_v = (mu0 > 0.0) ? Q_v : 0.0;
          Real Qv = kap_v*rho*Fstar*fac*exp(-tau_v);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
//          u0(m,IEN,k,j,i) += Q_v*bdt;
        
          tau_ir_down[i] = tau_ir;
          tau_ir += dtau_ir/mu1;
          if (i==is) tau_ir_down[i-1] = tau_ir;
        }
        p = w0(m,IEN,k,j,is-1)*gm1;
        rho = w0(m,IDN,k,j,is-1);
        T = p/Rgas/rho;
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        Real rtop = CellCenterX(ie+1-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rtop);
        I_down[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real r = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,r);
          Real rp1 = CellCenterX(i+1-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rp1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rp1/r);
          Real dtau = tau_ir_down[i]-tau_ir_down[i+1];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i)+B(i+1))/2.0;
          I_down[i] = (I_down[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        Real rbot = CellCenterX(is-1-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rbot);
        Real I_up = Iint;
        // up-sweep
        for (int i=is; i<ie+1; ++i) {
          Real r = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,r);
          Real rm1 = CellCenterX(i-1-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,rm1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rm1/r);
          Real dtau = tau_ir_down[i-1]-tau_ir_down[i];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i-1)+B(i))/2.0;
          I_up = (I_up*trans + Bavg*(1.0-trans))*fac;
          Real J = (I_up+I_down[i])/2.0;
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
            Real kap_ir = 1.0e-2; //2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
//          Real Q_ir = 4.0*M_PI*kap_ir*rho*(J-B(i));
//          u0(m,IEN,k,j,i) += Q_ir*bdt;
            
          Real cv = Rgas*rho*igm1;
          Real e0 = cv*T;
          Real kk = -4.0*M_PI*kap_ir*rho*boltz_sigma/M_PI*bdt;
          Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
          Real e;
          // Newton-Raphson
          for (int n=0; n<100; ++n) {
            e = cv*T;
            Real de = e - kk*SQR(SQR(T)) - bb;
            T -= de / (cv - 4.0*kk*T*T*T);
            if (fabs(de) <= 1.0e-10*e)
              break;
          }
          u0(m,IEN,k,j,i) += e-e0;
        }
        
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
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
//    Real Tn, Td;
//    get_daynight_Tp(p, Tn, Td);
//
////    T = Tn;
//
//    Real Tn2 = Tn*Tn;
//    Real Tn4 = Tn2*Tn2;
//    Real Td2 = Td*Td;
//    Real Td4 = Td2*Td2;
//    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;
//
//    Real Tmid = sqrt(sqrt(Tmid4));
//    T = Tmid;
//
////    Real bar = 1.0e6;
////    Real Teq = 1469.0;
//////    Real Tirr = Teq*sqrt(2);
//////    Real Tint = 100.0;
//////    Real g = 942.0;
//////    Real mus = cos(50.0/90.0*M_PI/2.0);
//////    Real fH = 0.5;
//////    Real fK = 1.0/3.0;
//////    Real kap_v = 4.0e-3;
//////    Real kap_ir = 1.0e-2;
//////    Real gam = kap_v/kap_ir;
//////    Real tau = p/(g/kap_ir);
//////    tau = (tau < 0.0) ? 0.0 : tau;
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
////////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
//////    T = sqrt(sqrt(T4)); //2581.5574;
////    T = Teq;
    
    
    return;
}

template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    
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
//    Real Teq = 1469.0;
////    Real Tirr = Teq*sqrt(2);
////    Real Tint = 100.0;
////    Real g = 942.0;
////    Real mus = cos(50.0/90.0*M_PI/2.0);
////    Real fH = 0.5;
////    Real fK = 1.0/3.0;
////    Real kap_v = 4.0e-3;
////    Real kap_ir = 1.0e-2;
////    Real gam = kap_v/kap_ir;
////    Real tau = p/(g/kap_ir);
////    tau = (tau < 0.0) ? 0.0 : tau;
////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
////    T = sqrt(sqrt(T4)); //2581.5574;
//    T = Teq;
    
    return;
}


template <typename View1D>
void get_wb_eos_arr(const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr) {
    
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
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
void get_wb_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
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
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

template <typename View1D>
void get_init_eos_arr(const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr) {

//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_init_Tp_host(N, Tarr, lgparr, p, T);

        logparr(n+1) = logparr(n) + fac/T;
    }


    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {

//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tp(N, Tarr, lgparr, p,T);
        rho = p/Rgas/T;
    } else {
        Real T0;
        get_init_Tp(N, Tarr, lgparr, p0,T0);
        Real rho0 = p0/Rgas/T0;
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A) {
  // Parmentier+2015
  Real X = log10(Teff);
  Real g = gg*0.01;
  Real a,b;
  if (Teff < 250.0) {
    a = -0.335*pow(g,0.07);
    b = 0.0;
  } else if (Teff < 750.0) {
    a = -0.335*pow(g,0.07) + 2.149*pow(g,0.135);
    b = -0.896*pow(g,0.135);
  } else if (Teff < 1250.0) {
    a = -0.335*pow(g,0.07) - 0.428*pow(g,0.135);
    b = 0.0;
  } else {
    a = 16.947 - 3.174*pow(g,0.07) - 4.051*pow(g,0.135);
    b = -5.472 + 0.917*pow(g,0.07) + 1.170*pow(g,0.135);
  }
  Real log10A = a + b*X;
  A = pow(10.0,log10A);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2) {
  // Parmentier & Giollot 2014; Parmentier+2015; Roth+2024
  Real X = log10(Teff);
  Real a3, a2, a1, b3, b2, b1, ab, bb;
  
  Real ap = -2.36;
  Real bp = 13.92;
  Real cp = -19.38;
  if (Teq < 1800.0 && Teff >= 1400.0) {
    ap = -12.45;
    bp = 82.25;
    cp = -134.42;
  }
  
  if (Teff < 2000.0) {
    ab = 0.84;
    bb = 0.0;
  } else {
    ab = 6.21;
    bb = -1.63;
  }
  if (Teq < 1800.0 && Teff >= 1400.0) {
    ab = 3.0;
    bb = -0.69;
  }
  
  if (Teff < 200.0) {
    a3 = -3.03;
    b3 = -0.2;
    a2 = -7.37;
    b2 = 2.53;
    a1 = -5.51;
    b1 = 2.48;
  } else if (Teff < 300.0) {
    a3 = -13.87;
    b3 = 4.51;
    a2 = 13.99;
    b2 = -6.75;
    a1 = 1.23;
    b1 = -0.45;
  } else if (Teff < 600.0) {
    a3 = -11.95;
    b3 = 3.74;
    a2 = -15.18;
    b2 = 5.02;
    a1 = 8.65;
    b1 = -3.45;
  } else if (Teff < 1400.0) {
    a3 = -6.97;
    b3 = 1.94;
    a2 = -10.41;
    b2 = 3.31;
    a1 = -12.96;
    b1 = 4.33;
  } else if (Teff < 2000.0) {
    a3 = -3.65;
    b3 = 0.89;
    a2 = -19.95;
    b2 = 6.34;
    a1 = -23.75;
    b1 = 7.76;
  } else {
    a3 = -6.02;
    b3 = 1.61;
    a2 = 13.56;
    b2 = -3.81;
    a1 = 12.65;
    b1 = -3.27;
  }
  if (Teq < 1800.0  && Teff >= 1400.0) {
    if (Teff < 2000.0) {
      a3 = 0.02;
      b3 = -0.28;
      a2 = 6.96;
      b2 = -2.21;
      a1 = -1.68;
      b1 = 0.75;
    } else {
      a3 = -16.54;
      b3 = 4.74;
      a2 = -2.4;
      b2 = 0.62;
      a1 = 10.37;
      b1 = -2.91;
    }
  }
  Real log10gamv1 = a1 + b1*X;
  Real log10gamv2 = a2 + b2*X;
  Real log10gamv3 = a3 + b3*X;
  Real log10gamp = ap*SQR(X) + bp*X + cp;
  beta = ab + bb*X;
    
  gamv1 = pow(10.0,log10gamv1);
  gamv2 = pow(10.0,log10gamv2);
  gamv3 = pow(10.0,log10gamv3);
  Real gamp = pow(10.0,log10gamp);
  Real dum = (gamp-1.0)/(2.0*beta*(1.0-beta));
  Real R = 1.0 + dum + sqrt(SQR(dum)+dum);
  gamir1 = beta + R - beta*R;
  gamir2 = gamir1/R;
    
  return;
}

KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr) {
  // Freedman+2014
  Real T1 = T;
  Real p1 = p;
  if (T > 4000.0) T1 = 4000.0;
  if (T < 75.0) T1 = 75.0;
  if (p > 3.0e8) p1 = 3.0e8;
  if (p < 1.0) p1 = 1.0;
  Real lgT = log10(T1);
  Real lgp = log10(p1);
  Real c1 = 10.602;
  Real c2 = 2.882;
  Real c3 = 6.09e-15;
  Real c4 = 2.954;
  Real c5 = -2.526;
  Real c6 = 0.843;
  Real c7 = -5.490;
  Real c8, c9, c10, c11, c12, c13;
  if (T1 < 800.0) {
    c8 = -14.051;
    c9 = 3.055;
    c10 = 0.024;
    c11 = 1.877;
    c12 = -0.445;
    c13 = 0.8321;
  } else {
    c8 = 82.241;
    c9 = -55.456;
    c10 = 8.754;
    c11 = 0.7048;
    c12 = -0.0414;
    c13 = 0.8321;
  }
    
  Real lgkl = c1*atan(lgT-c2) - c3/(lgp+c4)*exp(SQR(lgT-c5)) + c6*met + c7;
  Real lgkh = c8 + c9*lgT + c10*SQR(lgT) + lgp*(c11+c12*lgT) + c13*met*(0.5+1.0/M_PI*atan((lgT-2.5)/0.2));
  Real kl = pow(10.0,lgkl);
  Real kh = pow(10.0,lgkh);
  kapr = kl + kh;
  return;
}

KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint) {
  // Thorngren+2019 +Erratum
  Real boltz_sigma = 5.6704e-5;
  Real F = 4.0*boltz_sigma*SQR(SQR(Teq));
  Tint = 0.39*Teq*exp(-SQR(log10(F)-9.0-0.14)/1.095);
  Tint = (Tint < 100.0) ? 100.0 : Tint;
  return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real Teq = Tirr/sqrt(2);
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);
    
    Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
    Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
    get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
    
    Real R = gamir1/gamir2;
    Real gamp = gamir1 + gamir2 - SQR(gamir2)*R;
    taulim = sqrt(R/3.0)*sqrt(beta*SQR(R-1.0)-SQR(beta*(R-1.0))+R)/SQR(gamir1);
    Real At1 = SQR(gamir1)*log(1.0+1.0/(taulim*gamir1));
    Real At2 = SQR(gamir2)*log(1.0+1.0/(taulim*gamir2));
    
    Real a0 = 1.0/gamir1 + 1.0/gamir2;
    Real a1 = -1.0/(3.0*SQR(taulim))*(gamp/(1.0-gamp)*(gamir1+gamir2-2.0)/(gamir1+gamir2) + (gamir1+gamir2)*taulim - (At1+At2)*SQR(taulim));
    Real b0 = 1.0/(gamir1*gamir2/(gamir1-gamir2)*(At1-At2)/3.0 - SQR(gamir1*gamir2)/sqrt(3.0*gamp) - SQR(gamir1*gamir2)*gamir1*gamir2/(1.0-gamir1)/(1.0-gamir2)/(gamir1+gamir2));
    A = 1.0/3.0*(a0+a1*b0);
    B = -1.0/3.0*SQR(gamir1*gamir2)/gamp*b0;
    
//    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    
    for (int iv=0; iv<3; ++iv) {
        Real gamv;
        if (iv==0) gamv = gamv1;
        if (iv==1) gamv = gamv2;
        if (iv==2) gamv = gamv3;
        gamvv[iv] = gamv;
        
        Real longf = (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(gamir1+gamir2) - 3.0*gamv*(6.0*SQR(gamir1*gamir2)-SQR(gamv)*(SQR(gamir1)+SQR(gamir2)));
        Real a2 = SQR(taulim)/(gamp*SQR(gamv)) * longf/(1.0-SQR(gamv*taulim));
        Real Av1 = SQR(gamir1)*log(1.0+gamv/gamir1);
        Real Av2 = SQR(gamir2)*log(1.0+gamv/gamir2);
        Real a3 = -SQR(taulim) * (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(Av1+Av2) / (gamp*SQR(gamv)*gamv*(1.0-SQR(gamv*taulim)));
        Real b1 = gamir1*gamir2*(3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*SQR(taulim) / (gamp*SQR(gamv)*(SQR(gamv*taulim)-1.0));
        Real b2 = 3.0*(gamir1+gamir2)*SQR(gamv)*gamv / ((3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv)));
        Real b3 = (Av2-Av1) / (gamv*(gamir1-gamir2));
        
        C[iv] = -1.0/3.0*(b0*b1*(1.0+b2+b3)*a1 + a2 + a3);
        D[iv] = -B*b1*(1.0+b2+b3);
        E[iv] = (3.0-SQR(gamv/gamir1))*(3.0-SQR(gamv/gamir2)) / (9.0*gamv*(SQR(gamv*taulim)-1.0));
//        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C + D*exp(-tau/taulim) + E*exp(-gamv*tau));
    }
    
//    T = sqrt(sqrt(T4));
    return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real taulim, const Real &A, const Real &B, const Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamv)[3], const Real &tau, Real &T) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    for (int iv=0; iv<3; ++iv) {
        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C[iv] + D[iv]*exp(-tau/taulim) + E[iv]*exp(-gamv[iv]*tau));
    }
    T = sqrt(sqrt(T4));
    return;
}

template <typename View1D>
void get_picket_fence_pT_arr(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr) {
    Real bar = 1.0e6;
    Real tautop = 1.0e-6;
    Real Ttop;
    Real taulim, A, B, C[3], D[3], E[3], gamv[3];
    get_picket_fence_Ttau_coeff(Tint, Tirr, met, grav, mus, taulim, A, B, C, D, E, gamv);
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    Real kapr, ptop;
    ptop = 2.0; // 2e-6 bar
    get_kapr(Ttop, ptop, met, kapr);
    tautop = kapr/grav*ptop;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    
    Real lgptop = log10(ptop);
    Real pbot = 300.0*bar;
    Real lgpbot = log10(pbot);
    Real dlgp = (lgpbot-lgptop)/(N-1);
    Real lgp = lgptop;
    Real p = pow(10.0,lgp);
    Real tau = tautop;
    Real lgtau = log10(tau);
    Real T;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
    get_kapr(Ttop, p, met, kapr);
    for (int ip=0; ip<N; ++ ip) {
        Tarr(ip) = T;
        lgparr(ip) = lgp;
        Real K = p*kapr/grav*log(10.0);
        tau += K*dlgp;
        lgp += dlgp;
        p = pow(10.0,lgp);
//        lgp += dlgp;
//        Real pp = pow(10.0,lgp);
//        Real dp = pp-p;
//        Real K = kapr/grav;
//        tau += K*dp;
//        p = pp;
        get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
        get_kapr(T, p, met, kapr);
    }
    
    adjust_ad_pT_arr(N, Tarr, lgparr);
    
    return;
}

template <typename View1D>
void adjust_ad_pT_arr(const int &N, View1D Tarr, View1D lgparr) {

  // --- find convective boundary (search from top) ---
  int ic = -1;

  for (int ip = 0; ip < N-1; ++ip) {
    int i_inv = N - 1 - ip;
    int i_inv_p1 = N - 1 - (ip + 1);

    Real lgT  = log10(Tarr(i_inv));
    Real lgT1 = log10(Tarr(i_inv_p1));

    Real lgp  = lgparr(i_inv);
    Real lgp1 = lgparr(i_inv_p1);

    Real nabla = (lgT - lgT1) / (lgp - lgp1);

    Real T = Tarr(i_inv);
    Real nabla_ad = 0.32 - 0.1 * (T / 3000.0);

    if (nabla < nabla_ad) {
      ic = i_inv_p1;  // map back to original indexing
      break;
    }
  }

  // fallback if no crossing found
  if (ic > 0) {
    // --- enforce adiabat downward ---
    for (int ip = ic; ip < N-1; ++ip) {
      Real T = Tarr(ip);
      Real nabla_ad = 0.32 - 0.1 * (T / 3000.0);

      Tarr(ip+1) = pow(10.0,
        nabla_ad * (lgparr(ip+1) - lgparr(ip)) + log10(T)
      );
    }
  }
    
  return;
}

void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
    // Noti+2023; Lee+2021
    
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

    DvceArray5D<Real> u0, w0;
    u0 = pmbp->phydro->u0;
    w0 = pmbp->phydro->w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    const bool use_grid_stretch = pm->use_grid_stretch;
    const bool correct_spherical = false;
    const bool test_oned = false;
    
    Real r0, r1;
//    if (use_grid_stretch) {
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
//    }
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real met = 0.0;
    
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    Real met = pm->pgen->hot_jupiter_param.met;
    
    Real iap = 1.0/ap;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Tirr4 = SQR(SQR(Tirr));
    Real Fstar = boltz_sigma*Tirr4;
    Real Tint;
    get_Tint(Teq, Tint);
    Real Tint4 = SQR(SQR(Tint));
    Real Iint = boltz_sigma/M_PI*Tint4;
    
    Real mug[2];
    Real wg[2];
    mug[0] = 0.21132487;
    mug[1] = 0.78867513;
    wg[0] = 0.5;
    wg[1] = 0.5;
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);

//    size_t scr_size = 8 * ScrArray1D<Real>::shmem_size(n1);
//    int scr_level = 0;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, scr_level,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_down_r_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> B(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_ir_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> kapJ_ir(member.team_scratch(scr_level), n1);
        constexpr int NN = 270;
        Real tau_down_r_f[NN];
//        Real F_v_down_f[NN];
        Real B[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real F_ir_f[NN];
        Real Q_v[NN];
//        Real kapJ_ir[NN];
        
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        
        Real rtop = LeftEdgeX(ie+1-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rtop);
        Real rbot = LeftEdgeX(is-is, indcs.nx1, x1min, x1max);
        if (use_grid_stretch) StretchR(r0,r1,rbot);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        Real mus = (mu0 > 0.0) ? mu0 : 0.0;
        Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
        Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
        get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
        
        // 3 V Bands
        // top
        Real p = w0(m,IEN,k,j,ie+1)*gm1;
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = p/Rgas/rho;
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr;
        get_kapr(T, p, met, kapr);
        Real tau_r_f = kapr*p/grav;
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = tau_r_f/(kapr*rho);
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr;
          get_kapr(T, p, met, kapr);
          Real dr = dx1(m,k,j,i);
          tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
          Real r = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
          if (use_grid_stretch) StretchR(r0,r1,r);
////          Real delta = (drtop+(rtop-r))/r;
//          Real delta = dr/r;
//          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
//          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
          Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
          Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
//          Real trans1 = exp(-gamv1*tausl);
//          Real trans2 = exp(-gamv2*tausl);
//          Real trans3 = exp(-gamv3*tausl);
//          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
          Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
          Real Qv = kapr*rho*(1.0-albedo)*Fstar*1.0/3.0*(gamv1*trans1+gamv2*trans2+gamv3*trans3);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
        }
        
        // 2 IR Bands
        for (int i=is; i<ie+2; ++i) {
          F_ir_f[i] = 0.0;
//          kapJ_ir[i] = 0.0;
        }
        // two quadrature
        for (int n=0; n<2; ++n) {
          Real mugg = mug[n];
          Real wgg = wg[n];
          for (int vir=0; vir<2; ++vir) {
            Real gamir, fb;
            if (vir == 0) {
              gamir = gamir1;
              fb = beta;
            } else {
              gamir = gamir2;
              fb = 1.0-beta;
            }

            // top
            Real dtauir = gamir*tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/mugg);
            I_ir_down_f[ie+1] = (1.0-trans)*(fb*B[ie+1]);
            // down-sweep
            for (int i=ie; i>is-1; --i) {
              Real dtauir = gamir*(tau_down_r_f[i]-tau_down_r_f[i+1]);
              Real trans = exp(-dtauir/mugg);
              Real e0 = -expm1(-dtauir/mugg);
              Real e1 = dtauir/mugg - e0;
              Real alpa = 0.5*e0*(fb*B[i+1]+fb*B[i])/(fb*B[i+1]);
              Real alp = (dtauir > 1.0e-3) ? (e0 - e1/(dtauir/mugg)) : alpa;
              Real bet = (dtauir > 1.0e-3) ? (e1/(dtauir/mugg)) : 0.0;
              I_ir_down_f[i] = (1.0-e0)*I_ir_down_f[i+1] + alp*fb*B[i+1] + bet*fb*B[i];
            }
              
            // bottom
            I_ir_up_f[is] = Iint + I_ir_down_f[is];
            // up-sweep
            for (int i=is+1; i<ie+2; ++i) {
              Real dtauir = gamir*(tau_down_r_f[i-1]-tau_down_r_f[i]);
              Real trans = exp(-dtauir/mugg);
              Real e0 = -expm1(-dtauir/mugg);
              Real e1 = dtauir/mugg - e0;
              Real beto = 0.5*e0*(fb*B[i]+fb*B[i-1])/(fb*B[i]);
              Real bet = (dtauir > 1.0e-3) ? (e1/(dtauir/mugg)) : beto;
              Real gam = (dtauir > 1.0e-3) ? (e0 - e1/(dtauir/mugg)) : 0.0;
              I_ir_up_f[i] = (1.0-e0)*I_ir_up_f[i-1] + bet*fb*B[i] + gam*fb*B[i-1];
            }
              
            for (int i=is; i<ie+2; ++i) {
              Real F_ir_down_f = 2.0*M_PI*wgg*mugg*I_ir_down_f[i];
              Real F_ir_up_f = 2.0*M_PI*wgg*mugg*I_ir_up_f[i];
              F_ir_f[i] += F_ir_up_f - F_ir_down_f;
            }
//            for (int i=is; i<ie+1; ++i) {
//              Real J = wgg*I_ir_down_f[i] + wgg*I_ir_up_f[i];
//              Real p = w0(m,IEN,k,j,i)*gm1;
//              Real rho = w0(m,IDN,k,j,i);
//              Real T = p/Rgas/rho;
//              Real kapr;
//              get_kapr(T, p, met, kapr);
//              kapJ_ir[i] += gamir*kapr*J;
//            }
          }
        }
        
//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();
        
//        par_for_inner(member, is, ie, [&](const int i) {
        for (int i=is; i<ie+1; ++i) {
          // source term as flux divergence
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
            Real Ft = F_ir_f[i+1];//-F_v_down_f(i+1);
            Real Fb = F_ir_f[i];//-F_v_down_f(i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) {
            src = -(Ft*area_t-Fb*area_b)/vol;
          }
            src += Q_v[i];
          Real du_flux = src*bdt;
            
          // source term semi-implicit
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
          Real kapr;
          get_kapr(T, p, met, kapr);
          Real cv = Rgas*rho*igm1;
          Real e0 = cv*T;
          Real kk = 0.0;
          Real bb = du_flux + e0;
//          Real bb = Q_v(i)*bdt + e0;
          for (int vir=0; vir<2; ++vir) {
            Real gamir, fb;
            if (vir == 0) {
              gamir = gamir1;
              fb = beta;
            } else {
              gamir = gamir2;
              fb = 1.0-beta;
            }
            kk += -4.0*M_PI*gamir*kapr*rho*fb*boltz_sigma/M_PI*bdt;
            bb += 4.0*M_PI*gamir*kapr*rho*fb*B[i]*bdt;
          }
//          bb += 4.0*M_PI*rho*kapJ_ir(i)*bdt;
          int ierr=0;
          Real e;
          // Newton-Raphson
          for (int n=0; n<100; ++n) {
            e = cv*T;
            Real de = e - kk*SQR(SQR(T)) - bb;
            T -= de / (cv - 4.0*kk*T*T*T);
            if (T < 0.0) {
              e = e0;
              ierr = 1;
              break;
            }
            if (fabs(de) <= 1.0e-10*e)
              break;
          }
          Real du_src = e-e0;

          Real du = (fabs(du_flux) < e0 && ierr == 1) ? du_flux : du_src;
//          Real du = du_src;
//          Real du = du_flux;
          u0(m,IEN,k,j,i) += du;
        }
//        });
        
//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();
        
    });
    
    return;
}

