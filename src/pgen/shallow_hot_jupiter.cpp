//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file shallow_hot_jupiter.cpp
//! \brief Problem generator for the shallow hot Jupiter.
//!
//! REFERENCE: Heng, Menou, Phillipps, MNRAS, 413, 2380 (2011); Mendonça, Grimm, Grosheintz, Heng, ApJ, 829, 115 (2016); Ge, Li, Zhang, Lee, ApJ, 898, 130 (2020)

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
void get_init_Tz(const Real &z, Real &T);
template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);


//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  const bool use_etotgrav = pmy_mesh_->pmb_pack->phydro->use_etotgrav;
  const bool use_wellbalance = pmy_mesh_->pmb_pack->phydro->use_wellbalance;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  user_bcs_func = HydrostaticEquilibrium;
  if (restart) return;
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "shallow hot Jupiter problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

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
  Real ap = 1.0e10;
  Real omega = 2.1e-5;
  Real grav_acc = -800;
  Real p0 = 1e6;
  Real T0 = 1600;
  Real Rgas = 3.779e7;
  Real rho0 = p0/(Rgas * T0);
  Real H0 = -p0/rho0/grav_acc;
  Real iH0 = 1.0/H0;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_init_eos_arr(N, zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
  
    par_for("probini", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
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
        
      Real p = p0 * std::exp(-x3v*iH0);
      Real den = rho0 * std::exp(-x3v*iH0);
      get_init_eos(zarr.d_view,logparr.d_view,x3v,den,p);
      Real phicc = - grav_acc * x3v;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = p*igm1;
      if (use_etotgrav) {
          u0_(m,IEN,k,j,i) += den*phicc;
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          phi0_x2f(m,k,j,i) = phicc;
          if (j == je) {
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          phi0_x3f(m,k,j,i) = - grav_acc * x3v;
          if (k == ke) {
              x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              phi0_x3f(m,k+1,j,i) = - grav_acc * x3v;
          }
      }
        if (use_wellbalance) {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            p = p0 * std::exp(-x3v*iH0);
            Real denwb = rho0 * std::exp(-x3v*iH0);
            get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = p*igm1;
            if (i == ie) {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                p = p0 * std::exp(-x3v*iH0);
                denwb = rho0 * std::exp(-x3v*iH0);
                get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            p = p0 * std::exp(-x3v*iH0);
            denwb = rho0 * std::exp(-x3v*iH0);
            get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = p*igm1;
            if (j == je) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                p = p0 * std::exp(-x3v*iH0);
                denwb = rho0 * std::exp(-x3v*iH0);
                get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            p = p0 * std::exp(-x3v*iH0);
            denwb = rho0 * std::exp(-x3v*iH0);
            get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = p*igm1;
            if (k == ke) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
                p = p0 * std::exp(-x3v*iH0);
                denwb = rho0 * std::exp(-x3v*iH0);
                get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = p*igm1;
            }
        }
    });
    if (use_etotgrav) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {

            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            
            Real phicc = - grav_acc * x3v;
            phicc0(m,k,j,i) = phicc;
        });
        auto &mb_bcs = pmbp->pmb->mb_bcs;
        par_for("wbgravbc_x3", DevExeSpace(), 0,(pmbp->nmb_thispack-1),0,n2m1,0,n1m1,
        KOKKOS_LAMBDA(int m, int j, int i) {
          // apply physical boundaries to inner_x3
          switch (mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
            case BoundaryFlag::reflect:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ks-k-1,j,i) =  phicc0(m,ks+k,j,i);
              }
              break;
            case BoundaryFlag::outflow:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ks-k-1,j,i) = phicc0(m,ks,j,i);
              }
              break;
            case BoundaryFlag::diode:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ks-k-1,j,i) = phicc0(m,ks,j,i);
              }
              break;
            case BoundaryFlag::vacuum:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ks-k-1,j,i) = 0.0;
              }
              break;
            default:
              break;
          }

          // apply physical boundaries to outer_x3
          switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
            case BoundaryFlag::reflect:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ke+k+1,j,i) =  phicc0(m,ke-k,j,i);
              }
              break;
            case BoundaryFlag::outflow:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ke+k+1,j,i) = phicc0(m,ke,j,i);
              }
              break;
            case BoundaryFlag::diode:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ke+k+1,j,i) = phicc0(m,ke,j,i);
              }
              break;
            case BoundaryFlag::vacuum:
              for (int k=0; k<ng; ++k) {
                  phicc0(m,ke+k+1,j,i) = 0.0;
              }
              break;
            default:
              break;
          }
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

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
            
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            
          Real p = p0 * std::exp(-x3v*iH0);
          Real denwb = rho0 * std::exp(-x3v*iH0);
          get_init_eos(zarr.d_view,logparr.d_view,x3v,denwb,p);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = p*igm1;
          if (use_etotgrav) {
              Real phicc = - grav_acc * x3v;
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = p*igm1;
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
    DvceArray5D<Real> u0wb;
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
    const bool use_wellbalance = pmbp->phydro->use_wellbalance;

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
    Real grav_acc = -800;
    Real p0 = 1e6;
    Real T0 = 1600;
    Real Rgas = 3.779e7;
    Real rho0 = p0/(Rgas * T0);
    Real H0 = -p0/rho0/grav_acc;
    Real iH0 = 1.0/H0;
    
    u0wb = pmbp->phydro->u0wb;
    

  par_for("usrboundary", DevExeSpace(),0,(nmb-1),0,(ng-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
        if (use_wellbalance)
        {
            u0_(m,IDN,k,j,i) = u0wb(m,IDN,k,j,i);
            u0_(m,IM1,k,j,i) = 0.0;
            u0_(m,IM2,k,j,i) = 0.0;
            u0_(m,IM3,k,j,i) = 0.0;
            u0_(m,IEN,k,j,i) = u0wb(m,IEN,k,j,i);
        }
//        else
//        {
//            int nx3 = indcs.nx3;
//            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
//
//            Real p = p0 * std::exp(-x3v*iH0);
//            Real den = rho0 * std::exp(-x3v*iH0);
////            get_init_eos(x3v,den,p);
//            Real phicc = - grav_acc * x3v;
//
//            u0_(m,IDN,k,j,i) = den;
//            u0_(m,IM1,k,j,i) = 0.0;
//            u0_(m,IM2,k,j,i) = 0.0;
//            u0_(m,IM3,k,j,i) = 0.0;
//            u0_(m,IEN,k,j,i) = p*igm1;
//            if (use_etotgrav) {
//                u0_(m,IEN,k,j,i) += den*phicc;
//            }
//        }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
        if (use_wellbalance)
        {
            u0_(m,IDN,(ke+k+1),j,i) = u0wb(m,IDN,(ke+k+1),j,i);
            u0_(m,IM1,(ke+k+1),j,i) = 0.0;
            u0_(m,IM2,(ke+k+1),j,i) = 0.0;
            u0_(m,IM3,(ke+k+1),j,i) = 0.0;
            u0_(m,IEN,(ke+k+1),j,i) = u0wb(m,IEN,(ke+k+1),j,i);
        }
//        else
//        {
//            int nx3 = indcs.nx3;
//            Real x3v = CellCenterX((ke+k+1)-ks, nx3, x3min, x3max);
//
//            Real p = p0 * std::exp(-x3v*iH0);
//            Real den = rho0 * std::exp(-x3v*iH0);
////            get_init_eos(x3v,den,p);
//            Real phicc = - grav_acc * x3v;
//
//            u0_(m,IDN,(ke+k+1),j,i) = den;
//            u0_(m,IM1,(ke+k+1),j,i) = 0.0;
//            u0_(m,IM2,(ke+k+1),j,i) = 0.0;
//            u0_(m,IM3,(ke+k+1),j,i) = 0.0;
//            u0_(m,IEN,(ke+k+1),j,i) = p*igm1;
//            if (use_etotgrav) {
//                u0_(m,IEN,(ke+k+1),j,i) += den*phicc;
//            }
//        }
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
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
    const bool use_wellbalance = pmbp->phydro->use_wellbalance;
    
    Real ap = 1.0e10;
    Real iap = 1.0/ap;
    Real omega = 2.1e-5;
    Real grav_acc = -800;
    Real p0 = 1e6;
    Real T0 = 1600;
    Real Rgas = 3.779e7;
    Real rho0 = p0/(Rgas * T0);
    Real H0 = -p0/rho0/grav_acc;
    Real iH0 = 1.0/H0;
    Real gamma = pmbp->phydro->peos->eos_data.gamma;
    Real gm1 = gamma-1.0;
    
    // coefficients for Teq in Newtonian cooling
    Real Tsurf = 1600;
    Real Gtrop = 2.0e-6;
    Real dT = 10.0;
    Real dTep = 300.0;
    Real zstra = 2.0e8;
    Real Ps = 1e6;
//    Real sstra = p0 * std::exp(-zstra*iH0)/Ps;
    Real sstra = 0.125;
    
    // coefficients for trad in Newtonian cooling
    Real trad = 1.5e5;
    Real itrad = 1.0/trad;

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
        
        Real lam = x2v*iap;
        Real phi = x1v*iap;
        Real z = x3v;
        Real rho = w0(m,IDN,k,j,i);
        Real p = w0(m,IEN,k,j,i)*gm1;
        Real T = p/Rgas/rho;
        
        // gravity
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IM3,k,j,i);
        }
        if (use_wellbalance) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        u0(m,IM3,k,j,i) += src;
        
        // corotating beta-plane approximation e.g. Fromang+2016
        Real omega3 = omega*lam;
        u0(m,IM1,k,j,i) += -2.0*rho*omega3*(-w0(m,IM2,k,j,i))*bdt;
        u0(m,IM2,k,j,i) += -2.0*rho*omega3*w0(m,IM1,k,j,i)*bdt;

        // Newtonian cooling
        // Teq
        Real Tlow = Tsurf - Gtrop*(zstra+(z-zstra)/2.0) + sqrt(SQR(Gtrop*(z-zstra)/2.0) + SQR(dT));
        Real Thigh = Tsurf - Gtrop*zstra + dT;
        Real sig = p0 * std::exp(-z*iH0)/Ps;
        if (use_wellbalance) {
            sig = w0wb(m,IEN,k,j,i)*gm1/Ps;
        }
        Real blow = sin(M_PI/2.0*(sig-sstra)/(1.0-sstra));
        Real bhigh = 0.0;
        Real Tvert = (z > zstra)? Thigh : Tlow;
        Real btrop = (sig < sstra)? bhigh : blow;
        Real Teq = Tvert + btrop*dTep*cos(phi)*cos(lam); //exp(-SQR(lam/0.7)/2.0); //

        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;
    });

    return;
}


KOKKOS_INLINE_FUNCTION
void get_init_Tz(const Real &z, Real &T) {

    // coefficients for Tvert
    Real Tsurf = 1600;
    Real Gtrop = 2.0e-6;
    Real dT = 10.0;
    Real zstra = 2.0e8;
    Real Ps = 1e6;

    Real Tlow = Tsurf - Gtrop*(zstra+(z-zstra)/2.0) + sqrt(SQR(Gtrop*(z-zstra)/2.0) + SQR(dT));
    Real Thigh = Tsurf - Gtrop*zstra + dT;
    T = (z > zstra)? Thigh : Tlow;

    return;
}

template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr) {
    
    Real Rgas = 3.779e7;
    Real grav_acc = -800;
    Real T0 = 1600;
    Real p0 = 1e6;

//    const int N = 100;
//    DualArray1D<Real> zarr("zarr", N);
//    DualArray1D<Real> logparr("logparr", N);
//    DualArray1D<Real> Tarr("Tarr", N);
    Real zmin = 0.0;
    Real zmax = 5e8;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);
//    Tarr[0] = T0;
//    par_for("integrate", DevExeSpace(), 0, N,
//    KOKKOS_LAMBDA(const int n) {
//        Real zr = zmin+(n+0.5)*dz;
//        Real T;
//        zarr[n] = zmin+n*dz;
//        get_init_Tz(zr,T);
////        get_init_Tz(zarr[n],Tarr[n]);
//        if (n < N) {
//            logparr[n+1] = logparr[n] + fac/T;
//        }
//    });
    for(int n=0; n<N; n++) {
        Real zr = zmin + (n+0.5)*dz;
        Real T;
        zarr(n) = zmin + n*dz;

        get_init_Tz(zr, T);

        logparr(n+1) = logparr(n) + fac/T;
    }

    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
    Real Rgas = 3.779e7;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tz(z,T);
        rho = p/Rgas/T;
    } else {
        Real T0;
        get_init_Tz(0.0,T0);
        Real p0 = std::exp(logparr(0));
        Real rho0 = p0/Rgas/T0;
        Real grav_acc = -800;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}
