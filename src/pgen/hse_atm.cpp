//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hse_atm.cpp
//! \brief Problem generator for an atmosphere in hydrostatic equilibrium.
//!
//! Hydrostatic atmosphere test for well-balanced scheme and low-Mach Riemann solver
//! Sets up different initial conditions selected by flag "iprob"
//!   - iprob=1 : isothermal
//!   - iprob=2 : polytrope
//! REFERENCE: Edelmann et al., A&A, 652, A53 (2021)

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
void GravitySource(Mesh *pm, Real bdt);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the HSE atmosphere test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  const bool use_etotgrav = pmy_mesh_->pmb_pack->phydro->use_etotgrav;
  const bool use_wellbalance = pmy_mesh_->pmb_pack->phydro->use_wellbalance;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  int iprob  = pin->GetReal("problem","iprob");
  if (user_srcs) user_srcs_func = GravitySource;
  user_bcs_func = HydrostaticEquilibrium;
  if (restart) return;
  if (pmy_mesh_->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "hot bubble problem generator only works in 2D/3D" << std::endl;
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

  Real p0, gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
    p0 = 1.0/(pmbp->phydro->peos->eos_data.gamma);
    p0 = pin->GetOrAddReal("problem", "p0", p0);
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
    p0 = 1.0/(pmbp->pmhd->peos->eos_data.gamma);
    p0 = pin->GetOrAddReal("problem", "p0", p0);
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
  Real gdgm1 = gamma*igm1;
  Real grav_acc = pin->GetReal("problem","g0");
  Real rho0 = pin->GetReal("problem","rho0");
  Real H0 = -p0/rho0/grav_acc;
  Real iH0 = 1.0/H0;
  Real nu, inum1;
  if (iprob == 2) {
      nu = pin->GetReal("problem","nu");
      inum1 = 1.0/(nu-1.0);
      H0 *= nu*inum1;
  }

  // 2D PROBLEM ----------------------------------------------------------------

//  if (pmbp->pmesh->two_d) {
//    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
    par_for("hotbubble", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        
      Real p = p0 * std::exp(-x2v*iH0);
      Real den = rho0 * std::exp(-x2v*iH0);
      if (iprob == 2) {
          p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
          den = rho0 * std::pow(1.0-x2v*iH0,inum1);
      }
      Real phicc = - grav_acc * x2v;

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
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              phi0_x3f(m,k+1,j,i) = phicc;
          }
          x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
          phi0_x2f(m,k,j,i) = - grav_acc * x2v;
          if (j == je) {
              x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              phi0_x2f(m,k,j+1,i) = - grav_acc * x2v;
          }
      }
        if (use_wellbalance) {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            p = p0 * std::exp(-x2v*iH0);
            Real denwb = rho0 * std::exp(-x2v*iH0);
            if (iprob == 2) {
                p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
            }
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = p*igm1;
            if (i == ie) {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                p = p0 * std::exp(-x2v*iH0);
                denwb = rho0 * std::exp(-x2v*iH0);
                if (iprob == 2) {
                    p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                    denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
                }
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            p = p0 * std::exp(-x2v*iH0);
            denwb = rho0 * std::exp(-x2v*iH0);
            if (iprob == 2) {
                p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
            }
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = p*igm1;
            if (j == je) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
                p = p0 * std::exp(-x2v*iH0);
                denwb = rho0 * std::exp(-x2v*iH0);
                if (iprob == 2) {
                    p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                    denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
                }
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = p*igm1;
            }
            
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            p = p0 * std::exp(-x2v*iH0);
            denwb = rho0 * std::exp(-x2v*iH0);
            if (iprob == 2) {
                p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
            }
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = p*igm1;
            if (k == ke) {
                x1v = CellCenterX(i-is, nx1, x1min, x1max);
                x2v = CellCenterX(j-js, nx2, x2min, x2max);
                p = p0 * std::exp(-x2v*iH0);
                denwb = rho0 * std::exp(-x2v*iH0);
                if (iprob == 2) {
                    p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                    denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
                }
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
        par_for("hotbubblegrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {

            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
            
            Real phicc = - grav_acc * x2v;
            phicc0(m,k,j,i) = phicc;
        });
    }
    if (use_wellbalance) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("hotbubble", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          int nx1 = indcs.nx1;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
            
          Real p = p0 * std::exp(-x2v*iH0);
          Real denwb = rho0 * std::exp(-x2v*iH0);
            if (iprob == 2) {
                p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
                denwb = rho0 * std::pow(1.0-x2v*iH0,inum1);
            }
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = p*igm1;
          if (use_etotgrav) {
              Real phicc = - grav_acc * x2v;
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

//  // 3D PROBLEM ----------------------------------------------------------------
//
//  } else {
//    Kokkos::Random_XorShift64_Pool<> rand_pool64(pmbp->gids);
//    par_for("rt3d", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
//    KOKKOS_LAMBDA(int m, int k, int j, int i) {
//      Real &x1min = size.d_view(m).x1min;
//      Real &x1max = size.d_view(m).x1max;
//      int nx1 = indcs.nx1;
//      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
//
//      Real &x2min = size.d_view(m).x2min;
//      Real &x2max = size.d_view(m).x2max;
//      int nx2 = indcs.nx2;
//      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
//
//      Real &x3min = size.d_view(m).x3min;
//      Real &x3max = size.d_view(m).x3max;
//      int nx3 = indcs.nx3;
//      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
//
//      Real den=1.0;
//      if (x3v > 0.0) den *= drat;
//
//      if (iprob == 1) {
//        u0_(m,IM3,k,j,i) = (1.0+cos(kx*x1v))*(1.0+cos(ky*x2v))*(1.0+cos(kz*x3v))/8.0;
//      } else {
//        auto rand_gen = rand_pool64.get_state();  // get random number state this thread
//        Real r = 2.0*static_cast<Real>(rand_gen.frand()) - 1.0;
//        u0_(m,IM3,k,j,i) = r * (1.0 + cos(kz*x3v))/2.0;
//        rand_pool64.free_state(rand_gen);  // free state for use by other threads
//      }
//
//      u0_(m,IDN,k,j,i) = den;
//      u0_(m,IM1,k,j,i) = 0.0;
//      u0_(m,IM2,k,j,i) = 0.0;
//      u0_(m,IM3,k,j,i) *= (den*amp);
//      u0_(m,IEN,k,j,i) = (p0 + grav_acc*den*x3v)/gm1 + 0.5*SQR(u0_(m,IM3,k,j,i))/den;
//    });
//  }

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
    
    int iprob = 2;

    DvceArray5D<Real> u0_;
    const bool use_etotgrav = pmbp->phydro->use_etotgrav;

    Real p0, gamma;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      p0 = 1.0/(pmbp->phydro->peos->eos_data.gamma);
      p0 = 1.0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      p0 = 1.0/(pmbp->pmhd->peos->eos_data.gamma);
      p0 = 1.0;
    }
    Real gm1 = gamma - 1.0;
    Real igm1 = 1.0/gm1;
    Real gdgm1 = gamma*igm1;
    Real grav_acc = -1.0;
    Real rho0 = 1.0;
    Real H0 = -p0/rho0/grav_acc;
    Real iH0 = 1.0/H0;
    Real nu,inum1;
    if (iprob == 2) {
        nu = 1.6;
        inum1 = 1.0/(nu-1.0);
        H0 *= nu*inum1;
    }

  par_for("hotbubble_x2", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(ng-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
        int nx2 = indcs.nx2;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
          
        Real p = p0 * std::exp(-x2v*iH0);
        Real den = rho0 * std::exp(-x2v*iH0);
        if (iprob == 2) {
            p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
            den = rho0 * std::pow(1.0-x2v*iH0,inum1);
        }
        Real phicc = - grav_acc * x2v;

        u0_(m,IDN,k,j,i) = den;
        u0_(m,IM1,k,j,i) = 0.0;
        u0_(m,IM2,k,j,i) = 0.0;
        u0_(m,IM3,k,j,i) = 0.0;
        u0_(m,IEN,k,j,i) = p*igm1;
        if (use_etotgrav) {
            u0_(m,IEN,k,j,i) += den*phicc;
        }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
        int nx2 = indcs.nx2;
        Real x2v = CellCenterX((je+j+1)-js, nx2, x2min, x2max);
          
        Real p = p0 * std::exp(-x2v*iH0);
        Real den = rho0 * std::exp(-x2v*iH0);
        if (iprob == 2) {
            p = p0 * std::pow(1.0-x2v*iH0,nu*inum1);
            den = rho0 * std::pow(1.0-x2v*iH0,inum1);
        }
        Real phicc = - grav_acc * x2v;

        u0_(m,IDN,k,(je+j+1),i) = den;
        u0_(m,IM1,k,(je+j+1),i) = 0.0;
        u0_(m,IM2,k,(je+j+1),i) = 0.0;
        u0_(m,IM3,k,(je+j+1),i) = 0.0;
        u0_(m,IEN,k,(je+j+1),i) = p*igm1;
        if (use_etotgrav) {
            u0_(m,IEN,k,(je+j+1),i) += den*phicc;
        }
    }
  });
  return;
}


void GravitySource(Mesh *pm, Real bdt) {
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

    par_for("varying_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
      Real g = -1.0;
      Real src = bdt*g*w0(m,IDN,k,j,i);
      if (!use_etotgrav) {
        u0(m,IEN,k,j,i) += src*w0(m,IM2,k,j,i);
      }
      if (use_wellbalance) {
        src = bdt*g*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
      }
      u0(m,IM2,k,j,i) += src;
    });

    return;
}
