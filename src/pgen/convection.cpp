//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file convection.cpp
//! \brief Problem generator for the turbulent convection.
//!
//! REFERENCE: Leidi, Andrassy, Barsukow, Higl, Edelmann, Röpke, A&A, 686, A34 (2024)

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
void get_gz(const Real &z, Real &g);
KOKKOS_INLINE_FUNCTION
void get_gamz(const Real &z, Real &gam);
KOKKOS_INLINE_FUNCTION
void get_qdotz(const Real &z, Real &qdot);
template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr, View1D logrhoarr, View1D phiarr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const DvceArray1D<Real> &logrhoarr, const DvceArray1D<Real> &phiarr, const Real &z, Real &rho, Real &p, Real &phi);


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
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    DualArray1D<Real> logrhoarr("logrhoarr", N);
    DualArray1D<Real> phiarr("phiarr", N);
    get_init_eos_arr(N, zarr.h_view, logparr.h_view, logrhoarr.h_view, phiarr.h_view);
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    logrhoarr.modify_host();
    logrhoarr.sync_device();
    phiarr.modify_host();
    phiarr.sync_device();
  
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
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
        
      Real p, den, phicc, dummy;
      get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,den,p,phicc);
      Real qdot;
      get_qdotz(x3v,qdot);
      Real qdot0 = 3.795720e-5;
      Real &dx3 = size.d_view(m).dx3;
      qdot *= sin(4.0*M_PI*dx3)/(4.0*M_PI*dx3);
      Real drho = 1.1e-5*qdot/qdot0*(sin(3.0*M_PI*x1v)+cos(3.0*M_PI*x1v))*(sin(3.0*M_PI*x2v)-cos(M_PI*x2v));
      den += drho;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = p*igm1;
        
    });
    
    par_for("probwb", DevExeSpace(), 0, (pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
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
        
        Real p, den, phicc, dummy;
        get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,den,p,phicc);
        Real qdot;
        get_qdotz(x3v,qdot);
        Real qdot0 = 3.795720e-5;
        Real &dx3 = size.d_view(m).dx3;
        qdot *= sin(4.0*M_PI*dx3)/(4.0*M_PI*dx3);
        Real drho = 1.1e-5*qdot/qdot0*(sin(3.0*M_PI*x1v)+cos(3.0*M_PI*x1v))*(sin(3.0*M_PI*x2v)-cos(M_PI*x2v));
        den += drho;
        
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
          get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,dummy,dummy,phicc);
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,dummy,dummy,phicc);
              phi0_x3f(m,k+1,j,i) = phicc;
          }
      }
        if (use_wellbalance) {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            Real denwb;
            get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = p*igm1;
            if (i == ie) {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = p*igm1;
            if (j == je) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
                x3v = CellCenterX(k-ks, nx3, x3min, x3max);
                get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = p*igm1;
            if (k == ke) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
                get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,dummy);
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
            
            Real phicc, dummy;
            get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,dummy,dummy,phicc);
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
            
          Real p, denwb, phicc;
          get_init_eos(zarr.d_view,logparr.d_view,logrhoarr.d_view,phiarr.d_view,x3v,denwb,p,phicc);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = p*igm1;
          if (use_etotgrav) {
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
//  auto &indcs = pm->mb_indcs;
//  int &ng = indcs.ng;
//  int n1 = indcs.nx1 + 2*ng;
//  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
//  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
//  int &is = indcs.is;  int &ie  = indcs.ie;
//  int &js = indcs.js;  int &je  = indcs.je;
//  int &ks = indcs.ks;  int &ke  = indcs.ke;
//  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
//  int nmb = pm->pmb_pack->nmb_thispack;
//  MeshBlockPack *pmbp = pm->pmb_pack;
//  auto &size = pmbp->pmb->mb_size;
//
//    DvceArray5D<Real> u0_;
//    DvceArray5D<Real> u0wb;
//    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
//    const bool use_wellbalance = pmbp->phydro->use_wellbalance;
//
//    if (pmbp->phydro != nullptr) {
//      u0_ = pmbp->phydro->u0;
//    } else if (pmbp->pmhd != nullptr) {
//      u0_ = pmbp->pmhd->u0;
//    }
//
//    u0wb = pmbp->phydro->u0wb;
//
//
//  par_for("usrboundary", DevExeSpace(),0,(nmb-1),0,(ng-1),0,(n2-1),0,(n1-1),
//  KOKKOS_LAMBDA(int m, int k, int j, int i) {
//    Real &x3min = size.d_view(m).x3min;
//    Real &x3max = size.d_view(m).x3max;
//    if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
//        if (use_wellbalance)
//        {
//            u0_(m,IDN,k,j,i) = u0wb(m,IDN,k,j,i);
//            u0_(m,IM1,k,j,i) = 0.0;
//            u0_(m,IM2,k,j,i) = 0.0;
//            u0_(m,IM3,k,j,i) = 0.0;
//            u0_(m,IEN,k,j,i) = u0wb(m,IEN,k,j,i);
//        }
//    }
//    if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
//        if (use_wellbalance)
//        {
//            u0_(m,IDN,(ke+k+1),j,i) = u0wb(m,IDN,(ke+k+1),j,i);
//            u0_(m,IM1,(ke+k+1),j,i) = 0.0;
//            u0_(m,IM2,(ke+k+1),j,i) = 0.0;
//            u0_(m,IM3,(ke+k+1),j,i) = 0.0;
//            u0_(m,IEN,(ke+k+1),j,i) = u0wb(m,IEN,(ke+k+1),j,i);
//        }
//    }
//  });
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
    DvceArray4D<Real> phicc0;
    auto phi0_x3f = pmbp->phydro->phi0.x3f;
    u0 = pmbp->phydro->u0;
    w0 = pmbp->phydro->w0;
    w0wb = pmbp->phydro->w0wb;
    phicc0 = pmbp->phydro->phicc0;
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;
    const bool use_wellbalance = pmbp->phydro->use_wellbalance;
    const bool use_wellbalance_local = pmbp->phydro->use_wellbalance_local;

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
        
        Real grav_acc;
        get_gz(x3v,grav_acc);
        
        // gravity
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IM3,k,j,i);
        }
        if (use_wellbalance) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
//        if (use_wellbalance_local) {
////            src = bdt*(phicc0(m,k-1,j,i)-phicc0(m,k+1,j,i))/(2.0*size.d_view(m).dx3)*w0(m,IDN,k,j,i);
//            src = bdt*(phi0_x3f(m,k,j,i)-phi0_x3f(m,k+1,j,i))/(size.d_view(m).dx3)*w0(m,IDN,k,j,i);
//        }
        u0(m,IM3,k,j,i) += src;

        // cooling
        Real qdot;
        get_qdotz(x3v,qdot);
        Real &dx3 = size.d_view(m).dx3;
        qdot *= sin(4.0*M_PI*dx3)/(4.0*M_PI*dx3);
        u0(m,IEN,k,j,i) += qdot*bdt;
    });

    return;
}

KOKKOS_INLINE_FUNCTION
void get_gz(const Real &z, Real &g) {

    Real g0 = -1.414870;
    Real fg = 1.0;
    Real hg = 0.5*(1.0+sin(16.0*M_PI*(z-1.0/32.0)));
    fg = (z >= 1.0+1.0/16.0 && z <= 3.0-1.0/16.0) ? fg : hg;
    g = g0*fg*pow(z,-5.0/4.0);
    Real g1 = 0.5*(1.0+sin(16.0*M_PI*(1.0-1.0/32.0)));
    g = (z < 1.0) ? g1 : g;

    return;
}

KOKKOS_INLINE_FUNCTION
void get_gamz(const Real &z, Real &gam) {

    Real gam0 = 1.666667;
    Real gam1 = 1.3;
    Real etal = 0.0;
    Real etam = 0.5*(1.0+sin(8.0*M_PI*z));
    Real etah = 1.0;
    Real eta = (z < 2.0-1.0/16.0) ? etal : etam;
    eta = (z > 2.0+1.0/16.0) ? etah : eta;
    Real gamm = gam0 + eta*(gam1-gam0);
    
    gam = (z < 2.0-1.0/16.0) ? gam0 : gamm;
    gam = (z > 2.0+1.0/16.0) ? gam1 : gam;

    return;
}

KOKKOS_INLINE_FUNCTION
void get_qdotz(const Real &z, Real &qdot) {

    Real qdot0 = 3.795720e-5;
    Real qdotl = qdot0*sin(8*M_PI*z);
    Real qdoth = 0.0;
    qdot = (z < 1.0+1.0/8.0) ? qdotl : qdoth;

    return;
}

template <typename View1D>
void get_init_eos_arr(const int &N, View1D zarr, View1D logparr, View1D logrhoarr, View1D phiarr) {
    
    Real rho0 = 1.0;
    Real p0 = 1.0/1.666667;

    Real zmin = 1.0;
    Real zmax = 3.2;
    Real dz = (zmax-zmin)/N;
    Real gam, g, rho, p;
    logparr(0) = std::log(p0);
    logrhoarr(0) = std::log(rho0);
    phiarr(0) = 0.0;
    rho = rho0;
    p = p0;
    for(int n=0; n<N; n++) {
        Real zr = zmin + (n+0.5)*dz;
        zarr(n) = zmin + n*dz;

        get_gamz(zr, gam);
        get_gz(zr,g);

        Real dlogp = rho/p*g*dz;
        Real dlogrho = dlogp/gam;
        logparr(n+1) = logparr(n) + dlogp;
        logrhoarr(n+1) = logrhoarr(n) + dlogrho;
        phiarr(n+1) = phiarr(n) - g*dz;
        rho = exp(logrhoarr(n+1));
        p = exp(logparr(n+1));
    }
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const DvceArray1D<Real> &logrhoarr, const DvceArray1D<Real> &phiarr, const Real &z, Real &rho, Real &p, Real &phi) {
    
    Real dz = zarr(1)-zarr(0);
    Real zmin = 1.0;

    if (z >= zmin) {
        int Nt = std::floor((z-zmin)/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        Real logrho = logrhoarr(Nt) + (logrhoarr(Nt+1)-logrhoarr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        phi = phiarr(Nt) + (phiarr(Nt+1)-phiarr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        rho = std::exp(logrho);
    } else {
        Real rho0 = 1.0;
        Real p0 = 1.0/1.666667;
        Real g0;
        get_gz(1.0,g0);
        Real H0 = -p0/rho0/g0;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-(z-1.0)*iH0);
        rho = rho0 * std::exp(-(z-1.0)*iH0);
        phi = -g0*(z-1.0);
    }

    return;
}
