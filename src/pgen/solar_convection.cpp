//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file solar_convection.cpp
//! \brief Problem generator for the solar convection.
//!
//! REFERENCE:

// C++ headers
#include <cmath>
#include <iostream> // cout

// Athena++ headers
#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"
#include "diffusion/resistivity.hpp"

#include <Kokkos_Random.hpp>

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);

void two_stream_RT(Mesh *pm, Real bdt);
KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &rho, const Real &T, Real &kapr);

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
template <typename View1D>
void get_wb_eos_arr(const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr, View1D Tarr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const DvceArray1D<Real> &Tarr, const Real &z, Real &rho, Real &p);


//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  bool use_etotgrav = false;
  bool use_wellbalance_static = false;
  bool use_wellbalance_dynamic = false;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  user_bcs_func = HydrostaticEquilibrium;
  if (restart) return;
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "solar convection problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pin->GetInteger("mesh", "nx1") != pin->GetInteger("meshblock", "nx1")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "solar convection problem generator only allows one meshblock in r direction for the RT to work properly" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    
  Real r0, r1;
  r0 = pmy_mesh_->mesh_size.x1min;
  r1 = pmy_mesh_->mesh_size.x1max;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_;
  DvceArray5D<Real> w0_;
    
  DvceArray4D<Real> phi0_x1f;
  DvceArray4D<Real> phi0_x2f;
  DvceArray4D<Real> phi0_x3f;
  DvceArray4D<Real> phicc0;
  DvceArray5D<Real> u0wb;
  DvceArray5D<Real> w0wb;
  DvceArray5D<Real> w0facewb_x1f;
  DvceArray5D<Real> w0facewb_x2f;
  DvceArray5D<Real> w0facewb_x3f;
    
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

  Real gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    w0_ = pmbp->phydro->w0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
    use_etotgrav = pmbp->phydro->use_etotgrav;
    use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
    phi0_x1f = pmbp->phydro->phi0.x1f;
    phi0_x2f = pmbp->phydro->phi0.x2f;
    phi0_x3f = pmbp->phydro->phi0.x3f;
    phicc0 = pmbp->phydro->phicc0;
    u0wb = pmbp->phydro->u0wb;
    w0wb = pmbp->phydro->w0wb;
    w0facewb_x1f = pmbp->phydro->w0facewb.x1f;
    w0facewb_x2f = pmbp->phydro->w0facewb.x2f;
    w0facewb_x3f = pmbp->phydro->w0facewb.x3f;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    w0_ = pmbp->pmhd->w0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
    use_etotgrav = pmbp->pmhd->use_etotgrav;
    use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
    phi0_x1f = pmbp->pmhd->phi0.x1f;
    phi0_x2f = pmbp->pmhd->phi0.x2f;
    phi0_x3f = pmbp->pmhd->phi0.x3f;
    phicc0 = pmbp->pmhd->phicc0;
    u0wb = pmbp->pmhd->u0wb;
    w0wb = pmbp->pmhd->w0wb;
    w0facewb_x1f = pmbp->pmhd->w0facewb.x1f;
    w0facewb_x2f = pmbp->pmhd->w0facewb.x2f;
    w0facewb_x3f = pmbp->pmhd->w0facewb.x3f;
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
    
    Real Tbot = 1.0e4; //pin->GetReal("problem","Tbot");
    Real rhobot = 1.0e-5; //pin->GetReal("problem","rhobot");
    Real grav_acc = -2.74e4; //-pin->GetReal("problem","grav");
    Real ap = 6.957e10; //pin->GetReal("problem","ap");
    Real omega = 0.0; //pin->GetReal("problem","omega");
    Real Rgas = 1.38e8; //pin->GetReal("problem","Rgas");
  
    Real iap = 1.0/ap;
    
    Real grav = -grav_acc;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    DualArray1D<Real> Tarr("Tarr", N);
    get_wb_eos_arr(Rgas, grav_acc, N, (r1-r0)*1.1, zarr.h_view, logparr.h_view, Tarr.h_view);
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    Tarr.modify_host();
    Tarr.sync_device();
    Real z_seed = 0.45*(r1-r0);   // seed perturbations in the CZ (below the ~50%-height photosphere)
  
    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
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
      if (use_spherical_polar) {
        x1v = x1v_(m,i);
        x2v = x2v_(m,j);
        x3v = x3v_(m,k);
      } else {
        x1v = CellCenterX(i-is, nx1, x1min, x1max);
        x2v = CellCenterX(j-js, nx2, x2min, x2max);
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }
      Real r = x1v;
      if (use_spherical_polar) x1v -= ap;
        
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
      get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);

      auto rand_gen = rand_pool64.get_state();  // per-thread random state
      Real amp = 0.01;                          // white-noise perturbation amplitude
      Real dden = denwb*amp*2.0*(rand_gen.frand()-0.5);
      rand_pool64.free_state(rand_gen);
      dden = (x1v < z_seed) ? dden : 0.0;     // seed perturbations only in the convection zone
      Real p = pwb;
      Real den = denwb + dden;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = p*igm1;
        
      w0_(m,IDN,k,j,i) = den;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IEN,k,j,i) = p*igm1;
        
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
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        }
        Real r = x1v;
        if (use_spherical_polar) x1v -= ap;
          
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
        
      if (use_etotgrav || use_wellbalance_dynamic) {
          if (use_spherical_polar) {
            x1v = x1f_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
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
              if (use_spherical_polar) {
                x1v = x1f_(m,i+1);
              } else {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              }
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
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2f_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
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
              if (use_spherical_polar) {
                x2v = x2f_(m,j+1);
              } else {
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              }
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
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3f_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          }
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
              if (use_spherical_polar) {
                x3v = x3f_(m,k+1);
              } else {
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              }
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
        if (use_wellbalance_static) {
            if (use_spherical_polar) {
              x1v = x1f_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            Real denwb,pwb;
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = pwb*igm1;
            if (i == ie) {
                if (use_spherical_polar) {
                  x1v = x1f_(m,i+1);
                } else {
                  x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                }
                if (use_spherical_polar) x1v -= ap;
                get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = pwb*igm1;
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2f_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = pwb*igm1;
            if (j == je) {
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = pwb*igm1;
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3f_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = pwb*igm1;
            if (k == ke) {
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = pwb*igm1;
            }
        }
    });
    if (use_etotgrav || use_wellbalance_dynamic) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
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
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            Real r = x1v;
            if (use_spherical_polar) x1v -= ap;
              
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
    if (use_wellbalance_static) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
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
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          Real r = x1v;
          if (use_spherical_polar) x1v -= ap;
            
          Real pwb, denwb;
          get_wb_eos(Rgas, grav_acc, zarr.d_view,logparr.d_view,Tarr.d_view,x1v,denwb,pwb);
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
      Real bbot = pin->GetReal("problem","bbot");
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),0, n3m1, 0, n2m1, 0, n1m1, //ks,ke,js,je,is,ie,
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
            Real A2 = 0.0;
            Real A2ip = 0.0;
            Real A2kp = 0.0;
            Real A2ipkp = 0.0;
            Real A3 = 0.5*bbot*r0*sin(x2fl)/SQR(x1fl/r0);
            Real A3ip = 0.5*bbot*r0*sin(x2fl)/SQR(x1fr/r0);
            Real A3jp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fl/r0);
            Real A3ipjp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fr/r0);
//            Real A3 = 0.5*bbot*x1fl*sin(x2fl);
//            Real A3ip = 0.5*bbot*x1fr*sin(x2fl);
//            Real A3jp = 0.5*bbot*x1fl*sin(x2fr);
//            Real A3ipjp = 0.5*bbot*x1fr*sin(x2fr);
              
            b0.x1f(m,k,j,i) = (dxe3(m,k,j+1,i)*A3jp - dxe3(m,k,j,i)*A3)/area1(m,k,j,i) - (dxe2(m,k+1,j,i)*A2kp - dxe2(m,k,j,i)*A2)/area1(m,k,j,i);
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1)) {
              b0.x2f(m,k,j,i) = - ((x1fr*(x3fr-x3fl))*A3ip - (x1fl*(x3fr-x3fl))*A3) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0.x2f(m,k,j,i) = - (dxe3(m,k,j,i+1)*A3ip - dxe3(m,k,j,i)*A3)/area2(m,k,j,i);
            }
            b0.x3f(m,k,j,i) = (dxe2(m,k,j,i+1)*A2ip - dxe2(m,k,j,i)*A2)/area3(m,k,j,i);
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j,i+1)*A3ip)/area1(m,k,j,i+1) - (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k,j,i+1)*A2ip)/area1(m,k,j,i+1);
            Real b0x2fjp;
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js-1) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je)) {
              b0x2fjp = - ((x1fr*(x3fr-x3fl))*A3ipjp - (x1fl*(x3fr-x3fl))*A3jp) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0x2fjp = - (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j+1,i)*A3jp)/area2(m,k,j+1,i);
            }
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k+1,j,i)*A2kp)/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;
              
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
              
            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))+SQR(bcc0(m,IBY,k,j,i))+SQR(bcc0(m,IBZ,k,j,i)));
          }
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
    Real gamma;
    bool use_etotgrav = false;
    bool use_wellbalance_dynamic = false;
    
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    DvceArray5D<Real> bcc0;
    DvceArray4D<Real> b0_x1f;
    DvceArray4D<Real> b0_x2f;
    DvceArray4D<Real> b0_x3f;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &area2 = pmbp->pcoord->area.x2f;
    auto &area3 = pmbp->pcoord->area.x3f;
    auto &volume = pmbp->pcoord->volume;
    auto &z_ov_rE = pmbp->pcoord->z_ov_rE;
    Real grav_acc = -2.74e4;

    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      bcc0 = pmbp->pmhd->bcc0;
      b0_x1f = pmbp->pmhd->b0.x1f;
      b0_x2f = pmbp->pmhd->b0.x2f;
      b0_x3f = pmbp->pmhd->b0.x3f;
    }
    Real bbot = 0.0;
    
    int nvar = u0_.extent_int(1);
    
    Real igm1 = 1.0/(gamma-1.0);
    Real gm1 = gamma-1.0;
    Real gigm1 = gamma*igm1;
    Real gm1ig = (gamma-1.0)/gamma;
    Real ig = 1.0/gamma;
    // Fixed deep-adiabat reference for the S&N open-bottom INFLOW (this sets the energy input,
    // replacing volumetric heating). Anchored to the equilibrium IC base state (constant in
    // time) so the boundary does NOT feed back on the evolving interior -- an interior-tied
    // extrapolation pressurised the box and ran away (Tbase 13470->95000K).
    Real Rgas_bc = 1.38e8;
    Real p_base  = 4.38e6;      // IC base pressure at z=0 (C=20 equilibrium)
    Real T_base  = 13718.0;     // IC base temperature (deep adiabat, unchanged by opacity)
    Real rho_base = 2.31e-6;    // IC base density (C=20: denser than C=100's 7.92e-7)
    Real K_bot = p_base/pow(rho_base, gamma);   // deep-adiabat constant (unused here; bottom is CO5BOLD)
    Real T_top_fix = 3500.0;    // top: internal energy held constant at the temperature minimum
                                // (lowered from 4860 -> matches the cooler atmospheric top, so the
                                //  hydrostatic ghost is not over-dense and stops leaking mass in)
    (void)p_base; (void)T_base; (void)rho_base; (void)K_bot;
    
    if (pmbp->pmhd != nullptr) {
      par_for("usrboundaryx1_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = is+i;
              int iex = is-i-1;
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin+1)/area1(m,k,j,iex);
              b0_x1f(m,k,j,iex) = b0_x1f(m,k,j,iin+1)*fac1;
            }
          }
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = ie-i;
              int iex = ie+i+1;
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin)/area1(m,k,j,iex+1);
              b0_x1f(m,k,j,iex+1) = b0_x1f(m,k,j,iin)*fac1;
            }
        }
      });
    }
    
    if (pmbp->phydro != nullptr) {
      if (use_etotgrav) {
        pmbp->phydro->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->phydro->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }
    else if (pmbp->pmhd != nullptr) {
      auto b0 = pmbp->pmhd->b0;
      if (use_etotgrav) {
        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->pmhd->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }


    par_for("usrboundaryx1_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real rho_i = u0_(m,IDN,k,j,is);
          for (int i=0; i<ng; ++i) {
            if (pmbp->pmhd != nullptr) {
              u0_(m,IEN,k,j,is-i-1) -= 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
              Real lw, rw;
              lw = (x1f_(m,(is-i-1)+1)-x1v_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              rw = (x1v_(m,(is-i-1))-x1f_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              bcc0(m,IBX,k,j,(is-i-1)) = lw*b0_x1f(m,k,j,(is-i-1)) + rw*b0_x1f(m,k,j,(is-i-1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(is-i-1)) = lw*b0_x2f(m,k,j,(is-i-1)) + rw*b0_x2f(m,k,j+1,(is-i-1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(is-i-1)) = lw*b0_x3f(m,k,j,(is-i-1)) + rw*b0_x3f(m,k+1,j,(is-i-1));
              u0_(m,IEN,k,j,is-i-1) += 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
            }
            // Stein-Nordlund-style OPEN bottom: fixed-entropy inflow / zero-gradient outflow.
            // Upflows (vx>=0) enter carrying the fixed deep adiabat K_bot (the energy source);
            // downflows (vx<0) leave freely. Pressure follows an isothermal hydrostatic
            // extrapolation from the interior so the boundary stays ~hydrostatic (sets the mass).
            int igc = is-i-1;
            Real rho_is = u0_(m,IDN,k,j,is);
            Real vx_is = u0_(m,IM1,k,j,is)/rho_is;
            Real vy_is = u0_(m,IM2,k,j,is)/rho_is;
            Real vz_is = u0_(m,IM3,k,j,is)/rho_is;
            Real eint_is = u0_(m,IEN,k,j,is)
                         - 0.5*rho_is*(vx_is*vx_is+vy_is*vy_is+vz_is*vz_is);
            if (use_etotgrav) eint_is -= rho_is*phicc0(m,k,j,is);
            if (pmbp->pmhd != nullptr) {
              eint_is -= 0.5*(SQR(bcc0(m,IBX,k,j,is))+SQR(bcc0(m,IBY,k,j,is))+SQR(bcc0(m,IBZ,k,j,is)));
            }
            // Bottom ghosts: isothermal-hydrostatic extrapolation from the (CO5BOLD-relaxed)
            // bottom active layer; velocity zero-gradient. The CO5BOLD relaxation on the is layer
            // (in SourceFunc) sets the inflow entropy and enforces zero net mass flux.
            Real dphi_b = phicc0(m,k,j,igc) - phicc0(m,k,j,is);
            Real eint_g = exp(log(eint_is) - (rho_is/eint_is*igm1)*dphi_b);   // isothermal hydrostatic
            Real rho_g = eint_g/eint_is*rho_is;
            Real ke_g = 0.5*rho_g*(vx_is*vx_is + vy_is*vy_is + vz_is*vz_is);
            u0_(m,IDN,k,j,igc) = rho_g;
            u0_(m,IM1,k,j,igc) = rho_g*vx_is;
            u0_(m,IM2,k,j,igc) = rho_g*vy_is;
            u0_(m,IM3,k,j,igc) = rho_g*vz_is;
            u0_(m,IEN,k,j,igc) = eint_g + ke_g;
            if (use_etotgrav) u0_(m,IEN,k,j,igc) += rho_g*phicc0(m,k,j,igc);
            if (pmbp->pmhd != nullptr) {
              u0_(m,IEN,k,j,igc) += 0.5*(SQR(bcc0(m,IBX,k,j,igc))+SQR(bcc0(m,IBY,k,j,igc))+SQR(bcc0(m,IBZ,k,j,igc)));
            }
          }
        }
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
          Real rho_i = w0_(m,IDN,k,j,ie);
          Real e_i = w0_(m,IEN,k,j,ie);
          Real phi_i = phicc0(m,k,j,ie);
          Real q0_i = log(e_i);
          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            Real dM1mag = 0.0;
            if (pmbp->pmhd != nullptr) {
              Real lw, rw;
              lw = (x1f_(m,(ie+i+1)+1)-x1v_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
              rw = (x1v_(m,(ie+i+1))-x1f_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
              bcc0(m,IBX,k,j,(ie+i+1)) = lw*b0_x1f(m,k,j,(ie+i+1)) + rw*b0_x1f(m,k,j,(ie+i+1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(ie+i+1)) = lw*b0_x2f(m,k,j,(ie+i+1)) + rw*b0_x2f(m,k,j+1,(ie+i+1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(ie+i+1)) = lw*b0_x3f(m,k,j,(ie+i+1)) + rw*b0_x3f(m,k+1,j,(ie+i+1));

              Real pb = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
              Real pbp1 = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+2)))+SQR(bcc0(m,IBY,k,j,(ie+i+2)))+SQR(bcc0(m,IBZ,k,j,(ie+i+2))));
              Real M11 = pb - SQR(b0_x1f(m,k,j,(ie+i+1)));
              Real M11p1 = pbp1 - SQR(b0_x1f(m,k,j,(ie+i+2)));
              Real M12 = - b0_x2f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M12p1 = - b0_x2f(m,k,j+1,(ie+i+1)) * bcc0(m,IBX,k,j+1,(ie+i+1));
              Real M13 = - b0_x3f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M13p1 = - b0_x3f(m,k+1,j,(ie+i+1)) * bcc0(m,IBX,k+1,j,(ie+i+1));
              dM1mag = -( (M11p1*area1(m,k,j,(ie+i+2))-M11*area1(m,k,j,(ie+i+1))) + (M12p1*area2(m,k,j+1,(ie+i+1))-M12*area2(m,k,j,(ie+i+1))) + (M13p1*area3(m,k+1,j,(ie+i+1))-M13*area3(m,k,j,(ie+i+1))) )/volume(m,k,j,(ie+i+1));
              dM1mag += z_ov_rE(m,k,j,(ie+i+1)) * 0.5*SQR(bcc0(m,IBX,k,j,(ie+i+1)));
            }
            // S&N/CO5BOLD top (temperature minimum): HYDROSTATIC density gradient at the FIXED
            // T_top_fix, specific internal energy held CONSTANT (= cv*T_top_fix) in time & space,
            // velocity zero-gradient (open). The interior-tied isothermal extrapolation drained
            // and went singular (nan); this fixed-T lid is stable.
            Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
            Real rho0_ip = rho_i*exp(-dphi_i/(Rgas_bc*T_top_fix));
            Real e0_ip = rho0_ip*Rgas_bc*igm1*T_top_fix;   // internal energy density (const specific e)
            (void)q0_i; (void)factor_i; (void)dM1mag;
            u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
            // Stein-Nordlund-style OPEN/transmitting top: copy the vertical velocity so mass
            // flows BOTH ways (outflow of upflows, inflow to refill downflows). The old
            // outflow-only (fmax) form drained the ultra-thin top atmosphere to ~0 density.
            Real mom = u0_(m,IM1,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(ie+i+1)))+SQR(u0_(m,IM2,k,j,(ie+i+1)))+SQR(u0_(m,IM3,k,j,(ie+i+1))))/rho0_ip;
            if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(ie+i+1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
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
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    bool use_etotgrav = false;
    bool use_wellbalance_static = false;
    bool use_wellbalance_dynamic = false;
    Real gamma;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
      w0wb = pmbp->phydro->w0wb;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      w0wb = pmbp->pmhd->w0wb;
    }
    
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    
    Real grav_acc = -2.74e4;
    Real ap = 6.957e10;
    Real omega = 0.0;
    Real Rgas = 1.38e8;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;

    Real time = pm->time;

    // Bottom convective-flux heating layer: inject the stellar energy budget
    // F_base = sigma*Teff^4 as volumetric heating spread over the bottom fraction f_heat
    // of the box. This represents the convective enthalpy flux entering the deep CZ (the
    // box bottom is deep in the convection zone, so this energy is NOT radiative). It
    // balances the surface radiative cooling from two_stream_RT, so total energy should
    // PLATEAU rather than swing. Qheat [erg/cm^3/s] = F_base / heating-layer thickness.
    Real boltz_sigma = 5.6704e-5;
    Real Teff_heat = 5778.0;
    Real Fbase = boltz_sigma*SQR(SQR(Teff_heat));   // ~6.3e10 erg/cm^2/s
    Real f_heat = 0.10;                             // heat the bottom 10% of the box
    Real d_heat = f_heat*(r1-r0);
    Real z_heat_top = r0 + d_heat;
    Real Qheat = Fbase/d_heat;                      // uniform volumetric heating rate

    two_stream_RT(pm, bdt);

    par_for("usrsource", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        }
        
        Real area_r, area_l, vol;
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
          z = x1v-ap;
          r = x1v;
          area_r = area1(m,k,j,i+1);
          area_l = area1(m,k,j,i);
          vol = volume(m,k,j,i);
        } else {
          lam = x3v*iap;
          phi = x2v*iap;
          z = x1v;
        }
        Real rho = w0(m,IDN,k,j,i);
        Real p = w0(m,IEN,k,j,i)*gm1;
        Real T = p/Rgas/rho;
        
        // gravity
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance_static) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_dynamic) {
          Real e_imh,e_iph,dum1,dum2,dum3;
          if (pmbp->phydro != nullptr) {
              pmbp->phydro->getWBerho(IEN, gamma,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,e_imh,dum2,e_iph,dum3);
          } else if (pmbp->pmhd != nullptr) {
              pmbp->pmhd->getWBerho(IEN, gamma,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,e_imh,dum2,e_iph,dum3);
          }
          Real pl = e_imh*gm1;
          Real pr = e_iph*gm1;
          src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol;
        }
        u0(m,IM1,k,j,i) += src;

        // Top WAVE-DAMPING layer: relax velocities toward zero in the upper atmosphere so that
        // convection-driven acoustic/gravity waves are ABSORBED rather than reflected off the
        // fixed-T top (which was building standing modes -> vx_rms rising to ~2.7 km/s at the top).
        // Top ~20%, tau=50s (10x stronger + deeper than before) to absorb waves within ~one
        // crossing; ff^2 ramp avoids an impedance jump. Sits well above the photosphere (~45%).
        {
          Real zt = r1;
          Real zb = r1 - 0.20*(r1-r0);       // sponge covers the top ~20% of the domain
          Real tau_sponge = 5.0e1;           // velocity damping timescale [s] (was 500)
          Real ff = (z-zb)/(zt-zb);
          ff = (ff < 0.0) ? 0.0 : ((ff > 1.0) ? 1.0 : ff);
          Real fredux = ff*ff*bdt/tau_sponge;
          fredux = (fredux > 1.0) ? 1.0 : fredux;
          Real rhoc = u0(m,IDN,k,j,i);
          Real ke_old = 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))
                            +SQR(u0(m,IM3,k,j,i)))/rhoc;
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
          Real ke_new = 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))
                            +SQR(u0(m,IM3,k,j,i)))/rhoc;
          u0(m,IEN,k,j,i) += (ke_new-ke_old);   // drain damped KE from total energy
        }

        // Bottom convective-flux heating: DISABLED for now. The equilibrium IC already starts
        // balanced; the explicit base heating overheated the low-heat-capacity base before
        // convection could develop. Re-enable once startup is robust (ramp/thicker layer).
        if (false && x1v < z_heat_top) {
          u0(m,IEN,k,j,i) += Qheat*bdt;
        }

//        // Forces in the corotating frame
//        if (use_spherical_polar) {
//          Real vtheta = w0(m,IVY,k,j,i);
//          Real vphi = w0(m,IVZ,k,j,i);
//          Real vr = w0(m,IVX,k,j,i);
//          Real sine = sin(theta);
//          Real cosine = cos(theta);
//          Real oor = SQR(omega)*r*sine;
//          Real cor = 2.0*omega*vphi;
//          u0(m,IM2,k,j,i) += rho*(cor+oor)*cosine*bdt;
//          u0(m,IM3,k,j,i) += -rho*2.0*omega*(vr*sine+vtheta*cosine)*bdt;
//          u0(m,IM1,k,j,i) += rho*(cor+oor)*sine*bdt;
////          if (!use_etotgrav)
//          u0(m,IEN,k,j,i) += rho*oor*(vr*sine+vtheta*cosine)*bdt;
//        } else {
//          // corotating beta-plane approximation e.g. Fromang+2016
//          Real omega1 = omega*lam;
//          u0(m,IM2,k,j,i) += -2.0*rho*omega1*(-w0(m,IVZ,k,j,i))*bdt;
//          u0(m,IM3,k,j,i) += -2.0*rho*omega1*w0(m,IVY,k,j,i)*bdt;
//        }

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

//        // Rayleigh drag (initial relaxation)
//        Real dyntime = 2.0*M_PI/omega;
//        Real tau1 = dyntime / 10.0;
//        Real tau2 = dyntime;
//        Real t1   = 2.0 * dyntime;
//        Real t2   = 5.0 * dyntime;
        Real itdrag, fredux;
//        if(time < t1) {
//          itdrag = 1.0 / tau1;
//        } else if(time < t2)
//          {
//            Real alp = (time - t1) / (t2 - t1);
//            itdrag   = 1.0 / tau1 * pow(tau1 / tau2, alp);
//          }
//        else {
//          itdrag = 0.0;
//        }
//        if (time < t2) {
//          fredux = itdrag*bdt;///(1.0+itdrag*bdt);
//          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
//          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
//          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
//        }
//
//        // Top sponge layer
//        Real ztop = 1.25e8;
//        Real zbot = 1.2e8;
//        Real fdrag = (z-zbot)/(ztop-zbot);
//        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
//        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
//        itdrag = fdrag/1.0e4;
//        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
//        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
//        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
//        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
//        // Bottom sponge layer
//        logpl = log(1.0e2*bar);
//        logpt = log(5.0e1*bar);
//        logp = log(p);
//        fdrag = (logp-logpt)/(logpl-logpt); // high p = 1, low p = 0
//        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
//        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
//        itdrag = fdrag/1.0e3;
//        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
////        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
//        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
//        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
    });

    // ---- CO5BOLD (Freytag+ 2012) OPEN-BOTTOM relaxation on the lowest active layer (is) ----
    // Applied as a source term with the stage timestep bdt. Steps: (2) relax the entropy of
    // UPFLOWING gas toward the prescribed deep-adiabat s_inflow (energy input, sets Teff);
    // (3) damp pressure fluctuations toward the global horizontal mean <P>; (4) restore the mean
    // density (mass conservation / closed to plane-parallel waves); (5) subtract <rho*v3>/<rho>
    // from the vertical velocity (zero net mass flux). Horizontal means are GLOBAL over the plane.
    {
      Real igm1 = 1.0/gm1;
      Real cv = Rgas*igm1;
      Real rho_base = 2.31e-6, p_base = 4.38e6;   // C=20 equilibrium deep-adiabat base state
      Real K_bot = p_base/pow(rho_base, gamma);
      Real s_inflow = cv*log(K_bot);
      Real CsChange = 0.1, CPChange = 0.3;
      // Pass A: global <rho0>, <P>, N over the bottom plane
      Real srho0=0.0, sumP=0.0, sN=0.0;
      Kokkos::parallel_reduce("co5A", nmb1+1, KOKKOS_LAMBDA(int m, Real &a, Real &b, Real &c) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int k=ks; k<=ke; ++k) { for (int jj=js; jj<=je; ++jj) {
            Real r = u0(m,IDN,k,jj,is);
            Real ke = 0.5*(SQR(u0(m,IM1,k,jj,is))+SQR(u0(m,IM2,k,jj,is))+SQR(u0(m,IM3,k,jj,is)))/r;
            Real ei = u0(m,IEN,k,jj,is) - ke;
            if (use_etotgrav) ei -= r*phicc0(m,k,jj,is);
            a += r; b += ei*gm1; c += 1.0;
          } }
        }
      }, srho0, sumP, sN);
#if MPI_PARALLEL_ENABLED
      { Real g[3]={srho0,sumP,sN};
        MPI_Allreduce(MPI_IN_PLACE,g,3,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
        srho0=g[0]; sumP=g[1]; sN=g[2]; }
#endif
      Real meanP = sumP/sN;
      // Pass B: steps 2 (inflow entropy) + 3 (pressure damping)
      par_for("co5B", DevExeSpace(), 0, nmb1, ks, ke, js, je, KOKKOS_LAMBDA(int m,int k,int jj) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real dz = size.d_view(m).dx1;
          Real r = u0(m,IDN,k,jj,is);
          Real vx = u0(m,IM1,k,jj,is)/r, vy = u0(m,IM2,k,jj,is)/r, vz = u0(m,IM3,k,jj,is)/r;
          Real ei = u0(m,IEN,k,jj,is) - 0.5*r*(vx*vx+vy*vy+vz*vz);
          if (use_etotgrav) ei -= r*phicc0(m,k,jj,is);
          Real es = ei/r;
          Real P = ei*gm1, T = P/(r*Rgas), cs = sqrt(gamma*P/r);
          Real s = cv*log(P/pow(r,gamma));
          if (vx > 0.0) {                        // step 2: relax upflow entropy toward s_inflow
            Real rlx = CsChange*bdt*cs/dz, ss = s_inflow - s;
            r  += rlx*(-r*r*T*(gamma-1.0)/(P*gamma))*ss;
            es += rlx*T*(1.0/gamma)*ss;
          }
          Real P1 = r*es*gm1, cs2 = gamma*P1/r;   // step 3: damp pressure toward <P>
          Real rlxp = CPChange*bdt*sqrt(cs2)/dz;
          r  += rlxp*(1.0/cs2)*(meanP - P1);
          es += rlxp*(1.0/(gamma*r))*(meanP - P1);
          u0(m,IDN,k,jj,is) = r;
          u0(m,IM1,k,jj,is) = r*vx; u0(m,IM2,k,jj,is) = r*vy; u0(m,IM3,k,jj,is) = r*vz;
          Real E = r*es + 0.5*r*(vx*vx+vy*vy+vz*vz);
          if (use_etotgrav) E += r*phicc0(m,k,jj,is);
          u0(m,IEN,k,jj,is) = E;
        }
      });
      // Pass C: global <rho2>, <rho2*v3>, <v3>
      Real srho2=0.0, srho2v=0.0, sv=0.0;
      Kokkos::parallel_reduce("co5C", nmb1+1, KOKKOS_LAMBDA(int m, Real &a, Real &b, Real &c) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int k=ks; k<=ke; ++k) { for (int jj=js; jj<=je; ++jj) {
            Real r = u0(m,IDN,k,jj,is); Real vx = u0(m,IM1,k,jj,is)/r;
            a += r; b += r*vx; c += vx;
          } }
        }
      }, srho2, srho2v, sv);
#if MPI_PARALLEL_ENABLED
      { Real g[3]={srho2,srho2v,sv};
        MPI_Allreduce(MPI_IN_PLACE,g,3,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
        srho2=g[0]; srho2v=g[1]; sv=g[2]; }
#endif
      Real drho4 = (srho0 - srho2)/sN;                  // step 4: restore mean density (uniform)
      Real cvel  = (srho2v + drho4*sv)/srho0;           // step 5: <rho3 v3>/<rho0>
      // Pass D: steps 4 + 5
      par_for("co5D", DevExeSpace(), 0, nmb1, ks, ke, js, je, KOKKOS_LAMBDA(int m,int k,int jj) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real r0 = u0(m,IDN,k,jj,is);
          Real vx = u0(m,IM1,k,jj,is)/r0, vy = u0(m,IM2,k,jj,is)/r0, vz = u0(m,IM3,k,jj,is)/r0;
          Real ei = u0(m,IEN,k,jj,is) - 0.5*r0*(vx*vx+vy*vy+vz*vz);
          if (use_etotgrav) ei -= r0*phicc0(m,k,jj,is);
          Real es = ei/r0;
          Real r = r0 + drho4;                           // step 4
          Real vxn = vx - cvel;                          // step 5
          u0(m,IDN,k,jj,is) = r;
          u0(m,IM1,k,jj,is) = r*vxn; u0(m,IM2,k,jj,is) = r*vy; u0(m,IM3,k,jj,is) = r*vz;
          Real E = r*es + 0.5*r*(vxn*vxn+vy*vy+vz*vz);
          if (use_etotgrav) E += r*phicc0(m,k,jj,is);
          u0(m,IEN,k,jj,is) = E;
        }
      });
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    // Two-zone hydrostatic reference profile:
    //   (1) mildly superadiabatic convection zone (p >= p_ph): T = Tbot*(p/p0)^nabla_cz
    //   (2) isothermal, convectively stable atmosphere (p < p_ph): T = T_ph
    // For gamma=5/3, nabla_ad = (gamma-1)/gamma = 0.4, so nabla_cz=0.45 is only mildly
    // superadiabatic (excess 0.05); the isothermal atmosphere (nabla=0) is strongly stable.
    Real Tbot = 1.0e4;
    Real rhobot = 1.0e-5;
    Real Rgas = 1.38e8;
    Real p0 = Rgas*Tbot*rhobot;      // pressure at the bottom of the box
    Real nabla_cz = 0.45;            // convection-zone logarithmic T-p slope
    Real p_ph = 0.2*p0;              // photosphere pressure (CZ / atmosphere boundary)
    Real T_ph = pow(p_ph/p0,nabla_cz)*Tbot;
    T = (p >= p_ph) ? pow(p/p0,nabla_cz)*Tbot : T_ph;

    return;
}


// Build the initial reference profile in RADIATIVE-CONVECTIVE EQUILIBRIUM so the box starts
// balanced (surface radiates sigma*Teff^4 = the base heating Fbase) and does NOT heat/restructure.
// Integrated from the TOP (tau=0, skin temperature) downward:
//   * Atmosphere (radiative): gray Eddington relation  T^4 = (3/4) Teff^4 (tau + 2/3),
//     so T(tau=2/3) = Teff = the photosphere temperature (emergent flux = sigma*Teff^4).
//   * Interior (convective): once the radiative gradient exceeds nabla_ad (Schwarzschild),
//     switch to an adiabat T = T_sw (p/p_sw)^nabla_ad. Base T ~= Teff + g*depth/cp.
// Opacity is the code's own get_kapr, so the photosphere self-consistently lands where tau=2/3.
// A shooting loop on the top pressure p_top places the photosphere at ~50% of the box height.
template <typename View1D>
void get_wb_eos_arr(const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr, View1D Tarr) {

    Real g = -grav_acc;                       // positive gravitational acceleration
    Real gamma = 5.0/3.0;
    Real grad_ad = (gamma-1.0)/gamma;         // 0.4
    Real Teff = 5778.0;                       // photosphere effective temperature (matches Fbase)
    Real T_skin = Teff*std::pow(0.5, 0.25);   // Eddington skin temperature at tau=0
    Real dz = zmax/N;
    Real z_ph_target = 0.50*zmax/1.1;         // put photosphere ~50% up the physical box

    // one top-down integration for a given p_top; returns z where tau=2/3 (or -1), fills arrays if store
    auto integrate = [&](Real p_top, bool store) -> Real {
        Real tau = 0.0, T = T_skin, p = p_top;
        bool convective = false;
        Real p_sw = p, T_sw = T;
        Real z_ph = -1.0;
        for (int n=N-1; n>=0; --n) {
            if (store) { zarr(n) = n*dz; logparr(n) = std::log(p); Tarr(n) = T; }
            if (n==0) break;
            Real rho = p/(Rgas*T);
            Real kappa; get_kapr(rho, T, kappa);
            Real dtau = kappa*rho*dz;
            Real tau_new = tau + dtau;
            if (z_ph < 0.0 && tau_new >= 2.0/3.0) z_ph = n*dz;
            Real p_new = p + rho*g*dz;         // hydrostatic, downward
            Real T_new;
            if (!convective) {
                Real T_rad = Teff*std::pow(0.75*(tau_new+2.0/3.0), 0.25);
                Real grad = (std::log(T_rad)-std::log(T))/(std::log(p_new)-std::log(p));
                if (grad > grad_ad) { convective = true; p_sw = p; T_sw = T; T_new = T_sw*std::pow(p_new/p_sw, grad_ad); }
                else { T_new = T_rad; }
            } else {
                T_new = T_sw*std::pow(p_new/p_sw, grad_ad);
            }
            tau = tau_new; p = p_new; T = T_new;
        }
        return z_ph;
    };

    // shoot on p_top: higher p_top -> denser atmosphere -> tau=2/3 reached higher up (larger z)
    Real plo = 1.0e-4, phi = 1.0e8, p_top = 1.0;
    for (int it=0; it<80; ++it) {
        p_top = std::sqrt(plo*phi);
        Real z_ph = integrate(p_top, false);
        if (z_ph < 0.0 || z_ph < z_ph_target) plo = p_top;   // too little opacity -> raise p_top
        else phi = p_top;                                     // photosphere too high -> lower p_top
    }
    integrate(std::sqrt(plo*phi), true);
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_eos(const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const DvceArray1D<Real> &Tarr, const Real &z, Real &rho, Real &p) {

    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real w = (z-zarr(Nt))/(zarr(Nt+1)-zarr(Nt));
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))*w;
        T = Tarr(Nt) + (Tarr(Nt+1)-Tarr(Nt))*w;
        p = std::exp(logp);
        rho = p/Rgas/T;
    } else {
        // isothermal hydrostatic extrapolation below the base using the tabulated base state
        Real p0 = std::exp(logparr(0));
        Real T0 = Tarr(0);
        Real rho0 = p0/Rgas/T0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &rho, const Real &T, Real &kapr) {
  Real T4 = T/1.0e4;
  Real rhom6 = rho/1.0e-6;
  Real Tp3 = T4*T4*T4;
  Real Tp9 = Tp3*Tp3*Tp3;
  kapr = 2.0e1*sqrt(rhom6)*Tp9;   // opacity normalization C: 100 -> 20 (denser photosphere ->
                                  // slower, ~solar convective velocity; rho_ph ~ C^(-2/3))
  Real rhom5 = rho/1.0e-5;
  Real T5 = T/1.0e5;
  Real Tp7 = SQR(SQR(T5))*T5*T5*T5;
  Real kapd = 3.0e3*rhom5/sqrt(Tp7);
  kapr = fmin(kapr,kapd);
  Real kapsc = (T < 1.0e4) ? 0.0 : (0.2*(1.0+0.7));   // electron scattering (ionized gas)
  kapr += kapsc;
  
  kapr = fmax(kapr, 1.0e-2);      // OPACITY FLOOR: keep the cool atmosphere radiatively coupled
                                  // (kappa ~ T^9 otherwise collapses) so the RT DAMPS the
                                  // acoustic/gravity waves instead of letting them run up undamped
  return;
}

void two_stream_RT(Mesh *pm, Real bdt) {
    
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
    const bool use_spherical_polar = pm->use_spherical_polar;
        
    Real gamma;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
    }
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;

    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    
    Real grav = 2.74e4;
    Real ap = 6.957e10;
    Real Rgas = 1.38e8;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    // NO radiative flux injected at the base: the box bottom sits deep in the convection
    // zone, where radiation carries only ~0.2% of the total flux (F_rad ~ 1e8 vs
    // sigma*Teff^4 ~ 6e10). The sigma*Teff^4 energy enters as CONVECTIVE enthalpy flux via
    // the bottom heating layer in SourceFunc, not as radiation here. RT is a net cooler.
    Real Tint4 = 0.0;
    Real Iint = boltz_sigma/M_PI*Tint4;
    
    Real mug[2];
    Real wg[2];
    mug[0] = 0.21132487;
    mug[1] = 0.78867513;
    wg[0] = 0.5;
    wg[1] = 0.5;

    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
        constexpr int NN = 1040;
        Real tau_down_r_f[NN];
        Real B[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real F_ir_f[NN];
        
        // top
        Real p = w0(m,IEN,k,j,ie+1)*gm1;
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = p/Rgas/rho;
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr;
        get_kapr(rho, T, kapr);
        Real tau_r_f = 0.0; //kapr*p/grav;
        tau_down_r_f[ie+1] = tau_r_f;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real p = w0(m,IEN,k,j,i)*gm1;
          Real rho = w0(m,IDN,k,j,i);
          Real T = p/Rgas/rho;
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr;
          get_kapr(rho, T, kapr);
          Real dr;
          if (use_spherical_polar) {
            dr = dx1(m,k,j,i);
          } else {
            dr = size.d_view(m).dx1;
          }
          tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
        }
        
        // two-stream RT
        for (int i=is; i<ie+2; ++i) {
          F_ir_f[i] = 0.0;
        }
        // two quadrature
        for (int n=0; n<2; ++n) {
          Real mugg = mug[n];
          Real wgg = wg[n];

            // top
            Real dtauir = tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/mugg);
            I_ir_down_f[ie+1] = (1.0-trans)*B[ie+1];
            // down-sweep
            for (int i=ie; i>is-1; --i) {
              Real dtauir = tau_down_r_f[i]-tau_down_r_f[i+1];
              Real x = dtauir/mugg;
              Real e0 = -expm1(-x);
              Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              I_ir_down_f[i] = (1.0-e0)*I_ir_down_f[i+1] + alp*B[i+1] + bet*B[i];
            }
              
            // bottom
            I_ir_up_f[is] = Iint + I_ir_down_f[is];
            // up-sweep
            for (int i=is+1; i<ie+2; ++i) {
              Real dtauir = tau_down_r_f[i-1]-tau_down_r_f[i];
              Real x = dtauir/mugg;
              Real e0 = -expm1(-x);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              Real gam = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              I_ir_up_f[i] = (1.0-e0)*I_ir_up_f[i-1] + bet*B[i] + gam*B[i-1];
            }
              
            for (int i=is; i<ie+2; ++i) {
              Real F_ir_down_f = 2.0*M_PI*wgg*mugg*I_ir_down_f[i];
              Real F_ir_up_f = 2.0*M_PI*wgg*mugg*I_ir_up_f[i];
              F_ir_f[i] += F_ir_up_f - F_ir_down_f;
            }
        }
        
        for (int i=is; i<ie+1; ++i) {
          // source term as flux divergence
          Real Ft = F_ir_f[i+1];
          Real Fb = F_ir_f[i];
          Real dr;
          if (use_spherical_polar) {
            dr = dx1(m,k,j,i);
          } else {
            dr = size.d_view(m).dx1;
          }
          Real src = -(Ft-Fb)/dr;
          Real du_flux = src*bdt;
            
          // NOTE: implicit (Newton-Raphson) local-emission update temporarily reverted to
          // explicit to isolate its effect on atmospheric cooling. See git history / the
          // commented block below to re-enable.
          //          Real p = w0(m,IEN,k,j,i)*gm1;
          //          Real rho = w0(m,IDN,k,j,i);
          //          Real T = p/Rgas/rho;
          //          Real kapr;
          //          get_kapr(rho, T, kapr);
          //          Real cv = Rgas*rho*igm1;
          //          Real e0 = cv*T;
          //          Real kk = 0.0;
          //          Real bb = du_flux + e0;
          //          kk += -4.0*M_PI*kapr*rho*boltz_sigma/M_PI*bdt;
          //          bb += 4.0*M_PI*kapr*rho*B[i]*bdt;
          //          int ierr=0;
          //          Real e;
          //          for (int n=0; n<100; ++n) {
          //            e = cv*T;
          //            Real de = e - kk*SQR(SQR(T)) - bb;
          //            T -= de / (cv - 4.0*kk*T*T*T);
          //            if (T < 0.0) { e = e0; ierr = 1; break; }
          //            if (fabs(de) <= 1.0e-10*e) break;
          //          }
          //          Real du_src = e-e0;
          //          Real du = (fabs(du_flux) < e0 && ierr == 1) ? du_flux : du_src;
          Real du = du_flux;
          u0(m,IEN,k,j,i) += du;
        }
    });
    
    return;
}

