//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_bcs.cpp
//! \brief physical boundary conditions for the <rad_m1> moments (E, F1, F2, F3).
//!
//! The ghost fills read the CONSERVED radiation variables only, never any primitive
//! array (the restart-bitwise rule, design sect. 6).  Supported flags:
//!
//!   outflow        : zero-gradient copy of all four moments
//!   reflect        : mirror copy with the face-NORMAL flux component sign-flipped
//!   vacuum, diode  : free-streaming outflow with zero incoming radiation -- copy E and
//!                    the transverse flux, clamp the normal component so it can only
//!                    point OUT of the domain, then re-apply |F| <= c E
//!   inflow         : the constant state in u_in(n,face), set by the problem generator
//!   user           : left untouched here; ProblemGenerator::user_bcs_func fills it
//!   periodic       : handled by the mesh, never reaches this function

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

//! The per-cell ghost fill lives in rad_m1_closure.hpp as radm1::M1FillGhost, so that a
//! problem generator with `user` x1 walls can run exactly the same fill.

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValues::RadM1BCs()

void MeshBoundaryValues::RadM1BCs(MeshBlockPack *ppack, DualArray2D<Real> u_in,
                                  DvceArray5D<Real> u0, Real cl, Real e_floor) {
  auto &pm = ppack->pmesh;
  auto &indcs = ppack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  auto &mb_bcs = ppack->pmb->mb_bcs;

  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int nmb = ppack->nmb_thispack;
  int nvar = radm1::M1_NVAR;

  //--------------------------------------------------------------------------------- x1
  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    par_for("radm1bc_x1", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x1)) {
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,is-i-1,k,j,is,1,0,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,is-i-1,k,j,is+i,1,1,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::vacuum:
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,is-i-1,k,j,is,1,2,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            for (int n=0; n<nvar; ++n) {
              u0(m,n,k,j,is-i-1) = u_in.d_view(n,BoundaryFace::inner_x1);
            }
          }
          break;
        default:
          break;
      }

      switch (mb_bcs.d_view(m,BoundaryFace::outer_x1)) {
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,ie+i+1,k,j,ie,1,0,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,ie+i+1,k,j,ie-i,1,1,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::vacuum:
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            radm1::M1FillGhost(u0,m,k,j,ie+i+1,k,j,ie,1,2,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            for (int n=0; n<nvar; ++n) {
              u0(m,n,k,j,ie+i+1) = u_in.d_view(n,BoundaryFace::outer_x1);
            }
          }
          break;
        default:
          break;
      }
    });
  }
  if (pm->one_d) return;

  //--------------------------------------------------------------------------------- x2
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    int &js = indcs.js;
    int &je = indcs.je;
    par_for("radm1bc_x2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int i) {
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x2)) {
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,js-j-1,i,k,js,i,2,0,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,js-j-1,i,k,js+j,i,2,1,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::vacuum:
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,js-j-1,i,k,js,i,2,2,-1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            for (int n=0; n<nvar; ++n) {
              u0(m,n,k,js-j-1,i) = u_in.d_view(n,BoundaryFace::inner_x2);
            }
          }
          break;
        default:
          break;
      }

      switch (mb_bcs.d_view(m,BoundaryFace::outer_x2)) {
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,je+j+1,i,k,je,i,2,0,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,je+j+1,i,k,je-j,i,2,1,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::vacuum:
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            radm1::M1FillGhost(u0,m,k,je+j+1,i,k,je,i,2,2,1.0,cl,e_floor);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            for (int n=0; n<nvar; ++n) {
              u0(m,n,k,je+j+1,i) = u_in.d_view(n,BoundaryFace::outer_x2);
            }
          }
          break;
        default:
          break;
      }
    });
  }
  if (pm->two_d) return;

  //--------------------------------------------------------------------------------- x3
  if (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic) return;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  par_for("radm1bc_x3", DevExeSpace(), 0,(nmb-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int j, int i) {
    switch (mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ks-k-1,j,i,ks,j,i,3,0,-1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ks-k-1,j,i,ks+k,j,i,3,1,-1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::vacuum:
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ks-k-1,j,i,ks,j,i,3,2,-1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          for (int n=0; n<nvar; ++n) {
            u0(m,n,ks-k-1,j,i) = u_in.d_view(n,BoundaryFace::inner_x3);
          }
        }
        break;
      default:
        break;
    }

    switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ke+k+1,j,i,ke,j,i,3,0,1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ke+k+1,j,i,ke-k,j,i,3,1,1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::vacuum:
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          radm1::M1FillGhost(u0,m,ke+k+1,j,i,ke,j,i,3,2,1.0,cl,e_floor);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          for (int n=0; n<nvar; ++n) {
            u0(m,n,ke+k+1,j,i) = u_in.d_view(n,BoundaryFace::outer_x3);
          }
        }
        break;
      default:
        break;
    }
  });

  return;
}
