//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bfield_bcs.cpp
//  \brief

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
//! \!fn void BoundaryValues::BFieldBCs()
//! \brief Apply physical boundary conditions for all field variables at faces of MB which
//! are at the edge of the computational domain

void MeshBoundaryValues::BFieldBCs(MeshBlockPack *ppack, DualArray2D<Real> b_in,
                               DvceFaceFld4D<Real> b0) {
  // loop over all MeshBlocks in this MeshBlockPack
  auto &pm = ppack->pmesh;
  auto &indcs = ppack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  auto &mb_bcs = ppack->pmb->mb_bcs;

  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int nmb = ppack->nmb_thispack;
  auto &x2f_ = ppack->pcoord->xx2f;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    par_for("bfield-bc_x1", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
      // apply physical boundaries to inner_x1
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,is-i-1) = -b0.x1f(m,k,j,is+i+1);
            b0.x2f(m,k,j,is-i-1) =  b0.x2f(m,k,j,is+i);
            if (j == n2-1) {b0.x2f(m,k,j+1,is-i-1) = b0.x2f(m,k,j+1,is+i);}
            b0.x3f(m,k,j,is-i-1) =  b0.x3f(m,k,j,is+i);
            if (k == n3-1) {b0.x3f(m,k+1,j,is-i-1) = b0.x3f(m,k+1,j,is+i);}
          }
          break;
        case BoundaryFlag::outflow:
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,is-i-1) = b0.x1f(m,k,j,is);
            b0.x2f(m,k,j,is-i-1) = b0.x2f(m,k,j,is);
            if (j == n2-1) {b0.x2f(m,k,j+1,is-i-1) = b0.x2f(m,k,j+1,is);}
            b0.x3f(m,k,j,is-i-1) = b0.x3f(m,k,j,is);
            if (k == n3-1) {b0.x3f(m,k+1,j,is-i-1) = b0.x3f(m,k+1,j,is);}
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,is-i-1) = b_in.d_view(IBX,BoundaryFace::inner_x1);
            b0.x2f(m,k,j,is-i-1) = b_in.d_view(IBY,BoundaryFace::inner_x1);
            if (j == n2-1) {
              b0.x2f(m,k,j+1,is-i-1) = b_in.d_view(IBY,BoundaryFace::inner_x1);
            }
            b0.x3f(m,k,j,is-i-1) = b_in.d_view(IBZ,BoundaryFace::inner_x1);
            if (k == n3-1) {
              b0.x3f(m,k+1,j,is-i-1) = b_in.d_view(IBZ,BoundaryFace::inner_x1);
            }
          }
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x1
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,ie+i+2) = -b0.x1f(m,k,j,ie-i);
            b0.x2f(m,k,j,ie+i+1) =  b0.x2f(m,k,j,ie-i);
            if (j == n2-1) {b0.x2f(m,k,j+1,ie+i+1) = b0.x2f(m,k,j+1,ie-i);}
            b0.x3f(m,k,j,ie+i+1) =  b0.x3f(m,k,j,ie-i);
            if (k == n3-1) {b0.x3f(m,k+1,j,ie+i+1) = b0.x3f(m,k+1,j,ie-i);}
          }
          break;
        case BoundaryFlag::outflow:
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,ie+i+2) = b0.x1f(m,k,j,ie+1);
            b0.x2f(m,k,j,ie+i+1) = b0.x2f(m,k,j,ie);
            if (j == n2-1) {b0.x2f(m,k,j+1,ie+i+1) = b0.x2f(m,k,j+1,ie);}
            b0.x3f(m,k,j,ie+i+1) = b0.x3f(m,k,j,ie);
            if (k == n3-1) {b0.x3f(m,k+1,j,ie+i+1) = b0.x3f(m,k+1,j,ie);}
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            b0.x1f(m,k,j,ie+i+2) = b_in.d_view(IBX,BoundaryFace::outer_x1);
            b0.x2f(m,k,j,ie+i+1) = b_in.d_view(IBY,BoundaryFace::outer_x1);
            if (j == n2-1) {
              b0.x2f(m,k,j+1,ie+i+1) = b_in.d_view(IBY,BoundaryFace::outer_x1);
            }
            b0.x3f(m,k,j,ie+i+1) = b_in.d_view(IBZ,BoundaryFace::outer_x1);
            if (k == n3-1) {
              b0.x3f(m,k+1,j,ie+i+1) = b_in.d_view(IBZ,BoundaryFace::outer_x1);
            }
          }
          break;
        default:
          break;
      }
    });
  }
  if (pm->one_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    int &js = indcs.js;
    int &je = indcs.je;
    par_for("bfield-bc_x2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int i) {
      // apply physical boundaries to inner_x2
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,js-j-1,i) =  b0.x1f(m,k,js+j,i);
            if (i == n1-1) {b0.x1f(m,k,js-j-1,i+1) = b0.x1f(m,k,js+j,i+1);}
            b0.x2f(m,k,js-j-1,i) = -b0.x2f(m,k,js+j+1,i);
            b0.x3f(m,k,js-j-1,i) =  b0.x3f(m,k,js+j,i);
            if (k == n3-1) {b0.x3f(m,k+1,js-j-1,i) = b0.x3f(m,k+1,js+j,i);}
          }
          break;
        case BoundaryFlag::outflow:
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,js-j-1,i) = b0.x1f(m,k,js,i);
            if (i == n1-1) {b0.x1f(m,k,js-j-1,i+1) = b0.x1f(m,k,js,i+1);}
            b0.x2f(m,k,js-j-1,i) = b0.x2f(m,k,js,i);
            b0.x3f(m,k,js-j-1,i) = b0.x3f(m,k,js,i);
            if (k == n3-1) {b0.x3f(m,k+1,js-j-1,i) = b0.x3f(m,k+1,js,i);}
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,js-j-1,i) = b_in.d_view(IBX,BoundaryFace::inner_x2);
            if (i == n1-1) {
              b0.x1f(m,k,js-j-1,i+1) = b_in.d_view(IBX,BoundaryFace::inner_x2);
            }
            b0.x2f(m,k,js-j-1,i) = b_in.d_view(IBY,BoundaryFace::inner_x2);
            b0.x3f(m,k,js-j-1,i) = b_in.d_view(IBZ,BoundaryFace::inner_x2);
            if (k == n3-1) {
              b0.x3f(m,k+1,js-j-1,i) = b_in.d_view(IBZ,BoundaryFace::inner_x2);
            }
          }
          break;
        case BoundaryFlag::polar:
          b0.x2f(m,k,js,i) = 0.5*(b0.x2f(m,k,js-1,i)+b0.x2f(m,k,js+1,i)); // /cos(x2f_(m,js+1));
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x2
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,je+j+1,i) =  b0.x1f(m,k,je-j,i);
            if (i == n1-1) {b0.x1f(m,k,je+j+1,i+1) = b0.x1f(m,k,je-j,i+1);}
            b0.x2f(m,k,je+j+2,i) = -b0.x2f(m,k,je-j,i);
            b0.x3f(m,k,je+j+1,i) =  b0.x3f(m,k,je-j,i);
            if (k == n3-1) {b0.x3f(m,k+1,je+j+1,i) = b0.x3f(m,k+1,je-j,i);}
          }
          break;
        case BoundaryFlag::outflow:
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,je+j+1,i) = b0.x1f(m,k,je,i);
            if (i == n1-1) {b0.x1f(m,k,je+j+1,i+1) = b0.x1f(m,k,je,i+1);}
            b0.x2f(m,k,je+j+2,i) = b0.x2f(m,k,je+1,i);
            b0.x3f(m,k,je+j+1,i) = b0.x3f(m,k,je,i);
            if (k == n3-1) {b0.x3f(m,k+1,je+j+1,i) = b0.x3f(m,k+1,je,i);}
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            b0.x1f(m,k,je+j+1,i) = b_in.d_view(IBX,BoundaryFace::outer_x2);
            if (i == n1-1) {
              b0.x1f(m,k,je+j+1,i+1) = b_in.d_view(IBX,BoundaryFace::outer_x2);
            }
            b0.x2f(m,k,je+j+2,i) = b_in.d_view(IBY,BoundaryFace::outer_x2);
            b0.x3f(m,k,je+j+1,i) = b_in.d_view(IBZ,BoundaryFace::outer_x2);
            if (k == n3-1) {
              b0.x3f(m,k+1,je+j+1,i) = b_in.d_view(IBZ,BoundaryFace::outer_x2);
            }
          }
          break;
        case BoundaryFlag::polar:
          b0.x2f(m,k,je+1,i) = 0.5*(b0.x2f(m,k,je,i)+b0.x2f(m,k,je+2,i)); // /(-cos(x2f_(m,je)));
          break;
        default:
          break;
      }
    });
  }
//  if (ppack->pmesh->use_polar_boundary) PolarAzimuthalAverageBxBy(ppack,b0);
  if (pm->two_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic) return;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  par_for("bfield-bc_x3", DevExeSpace(), 0,(nmb-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int j, int i) {
    // apply physical boundaries to inner_x3
    switch (mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ks-k-1,j,i) =  b0.x1f(m,ks+k,j,i);
          if (i == n1-1) {b0.x1f(m,ks-k-1,j,i+1) = b0.x1f(m,ks+k,j,i+1);}
          b0.x2f(m,ks-k-1,j,i) =  b0.x2f(m,ks+k,j,i);
          if (j == n2-1) {b0.x2f(m,ks-k-1,j+1,i) = b0.x2f(m,ks+k,j+1,i);}
          b0.x3f(m,ks-k-1,j,i) = -b0.x3f(m,ks+k+1,j,i);
        }
        break;
      case BoundaryFlag::outflow:
      case BoundaryFlag::diode:
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ks-k-1,j,i) = b0.x1f(m,ks,j,i);
          if (i == n1-1) {b0.x1f(m,ks-k-1,j,i+1) = b0.x1f(m,ks,j,i+1);}
          b0.x2f(m,ks-k-1,j,i) = b0.x2f(m,ks,j,i);
          if (j == n2-1) {b0.x2f(m,ks-k-1,j+1,i) = b0.x2f(m,ks,j+1,i);}
          b0.x3f(m,ks-k-1,j,i) = b0.x3f(m,ks,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ks-k-1,j,i) = b_in.d_view(IBX,BoundaryFace::inner_x3);
          if (i == n1-1) {
            b0.x1f(m,ks-k-1,j,i+1) = b_in.d_view(IBX,BoundaryFace::inner_x3);
          }
          b0.x2f(m,ks-k-1,j,i) = b_in.d_view(IBY,BoundaryFace::inner_x3);
          if (j == n2-1) {
            b0.x2f(m,ks-k-1,j+1,i) = b_in.d_view(IBY,BoundaryFace::inner_x3);
          }
          b0.x3f(m,ks-k-1,j,i) = b_in.d_view(IBZ,BoundaryFace::inner_x3);
        }
        break;
      default:
        break;
    }

    // apply physical boundaries to outer_x3
    switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ke+k+1,j,i) =  b0.x1f(m,ke-k,j,i);
          if (i == n1-1) {b0.x1f(m,ke+k+1,j,i+1) = b0.x1f(m,ke-k,j,i+1);}
          b0.x2f(m,ke+k+1,j,i) =  b0.x2f(m,ke-k,j,i);
          if (j == n2-1) {b0.x2f(m,ke+k+1,j+1,i) = b0.x2f(m,ke-k,j+1,i);}
          b0.x3f(m,ke+k+2,j,i) = -b0.x3f(m,ke-k,j,i);
        }
        break;
      case BoundaryFlag::outflow:
      case BoundaryFlag::diode:
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ke+k+1,j,i) = b0.x1f(m,ke,j,i);
          if (i == n1-1) {b0.x1f(m,ke+k+1,j,i+1) = b0.x1f(m,ke,j,i+1);}
          b0.x2f(m,ke+k+1,j,i) = b0.x2f(m,ke,j,i);
          if (j == n2-1) {b0.x2f(m,ke+k+1,j+1,i) = b0.x2f(m,ke,j+1,i);}
          b0.x3f(m,ke+k+2,j,i) = b0.x3f(m,ke+1,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          b0.x1f(m,ke+k+1,j,i) = b_in.d_view(IBX,BoundaryFace::outer_x3);
          if (i == n1-1) {
            b0.x1f(m,ke+k+1,j,i+1) = b_in.d_view(IBX,BoundaryFace::outer_x3);
          }
          b0.x2f(m,ke+k+1,j,i) = b_in.d_view(IBY,BoundaryFace::outer_x3);
          if (j == n2-1) {
            b0.x2f(m,ke+k+1,j+1,i) = b_in.d_view(IBY,BoundaryFace::outer_x3);
          }
          b0.x3f(m,ke+k+2,j,i) = b_in.d_view(IBZ,BoundaryFace::outer_x3);
        }
        break;
      default:
        break;
    }
  });

  return;
}

void MeshBoundaryValues::PolarAzimuthalAverageBxBy(MeshBlockPack *ppack, DvceFaceFld4D<Real> b0) {
  auto bx2f = b0.x2f;
  auto bx3f = b0.x3f;
  auto &x2f_ = ppack->pcoord->xx2f;
  auto &x2v_ = ppack->pcoord->x2v;
  auto &x3f_ = ppack->pcoord->xx3f;
  auto &x3v_ = ppack->pcoord->x3v;

  auto &indcs = ppack->pmesh->mb_indcs;

  int nx1 = indcs.nx1;
  int nx3 = indcs.nx3;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  int n1m1 = nx1 + 2*ng - 1;
  int n3m1 = nx3 + 2*ng - 1;

  const int mesh_nx3 = ppack->pmesh->mesh_indcs.nx3;
  const int nmb = ppack->nmb_thispack;

  auto &mb_bcs = ppack->pmb->mb_bcs;

  // =====================================================
  // local device sums
  // =====================================================

  Kokkos::View<Real*> inner_local_bx("inner_local_bx", n1m1+1);
  Kokkos::View<Real*> outer_local_bx("outer_local_bx", n1m1+1);
  Kokkos::View<Real*> inner_local_by("inner_local_by", n1m1+1);
  Kokkos::View<Real*> outer_local_by("outer_local_by", n1m1+1);
//  Kokkos::View<Real*> inner_local_bxp("inner_local_bxp", n1m1+1);
//  Kokkos::View<Real*> outer_local_bxp("outer_local_bxp", n1m1+1);
//  Kokkos::View<Real*> inner_local_byp("inner_local_byp", n1m1+1);
//  Kokkos::View<Real*> outer_local_byp("outer_local_byp", n1m1+1);

  Kokkos::parallel_for(
      "polar_local_sum_b",
      Kokkos::RangePolicy<DevExeSpace>(0, n1m1),
      KOKKOS_LAMBDA(const int i) {

        Real inner_sum_x = 0.0;
        Real outer_sum_x = 0.0;
        Real inner_sum_y = 0.0;
        Real outer_sum_y = 0.0;
//        Real inner_sum_xp = 0.0;
//        Real outer_sum_xp = 0.0;
//        Real inner_sum_yp = 0.0;
//        Real outer_sum_yp = 0.0;

        for (int m = 0; m < nmb; ++m) {

          bool do_inner =
              (mb_bcs.d_view(m, BoundaryFace::inner_x2) ==
               BoundaryFlag::polar);

          bool do_outer =
              (mb_bcs.d_view(m, BoundaryFace::outer_x2) ==
               BoundaryFlag::polar);
            
          if (!(do_inner || do_outer)) continue;

          for (int k = ks; k <= ke; ++k) {
            Real bav, dum;
            if (do_inner) {
              bav = 0.5*(bx2f(m,k,js-1,i)+bx2f(m,k,js+1,i));
              inner_sum_x += 2.0*bav*cos(x3v_(m,k));
              inner_sum_y += 2.0*bav*sin(x3v_(m,k));
//              Real dxl = x2v_(m,js)-x2v_(m,js-1);
//              Real dxr = x2v_(m,js+1)-x2v_(m,js);
//              Real dxlh = x2v_(m,js)-x2f_(m,js);
//              Real dxrh = x2f_(m,js+1)-x2v_(m,js);
//              ppack->pmhd->PLM_nonuniform(bx3f(m,k,js-1,i),bx3f(m,k,js,i),bx3f(m,k,js+1,i),dxl,dxr,dxlh,dxrh,bav,dum);
//              inner_sum_x += -2.0*bav*sin(x3f_(m,k));
//              inner_sum_y += 2.0*bav*cos(x3f_(m,k));
//              bav = 0.5*(bx3f(m,k,js-1,i)+bx3f(m,k,js+1,i));
//              inner_sum_xp += -2.0*bav*sin(x3v_(m,k));
//              inner_sum_yp += 2.0*bav*cos(x3v_(m,k));
            }
            if (do_outer) {
              bav = 0.5*(bx2f(m,k,je,i)+bx2f(m,k,je+2,i));
              outer_sum_x += 2.0*bav*cos(x3v_(m,k));
              outer_sum_y += 2.0*bav*sin(x3v_(m,k));
//              Real dxl = x2v_(m,je)-x2v_(m,je-1);
//              Real dxr = x2v_(m,je+1)-x2v_(m,je);
//              Real dxlh = x2v_(m,je)-x2f_(m,je);
//              Real dxrh = x2f_(m,je+1)-x2v_(m,je);
//              ppack->pmhd->PLM_nonuniform(bx3f(m,k,je-1,i),bx3f(m,k,je,i),bx3f(m,k,je+1,i),dxl,dxr,dxlh,dxrh,dum,bav);
//              outer_sum_x += -2.0*bav*sin(x3f_(m,k));
//              outer_sum_y += 2.0*bav*cos(x3f_(m,k));
//              bav = 0.5*(bx3f(m,k,je,i)+bx3f(m,k,je+2,i));
//              outer_sum_xp += -2.0*bav*sin(x3v_(m,k));
//              outer_sum_yp += 2.0*bav*cos(x3v_(m,k));
            }
          }
        }

        inner_local_bx(i) = inner_sum_x;
        outer_local_bx(i) = outer_sum_x;
        inner_local_by(i) = inner_sum_y;
        outer_local_by(i) = outer_sum_y;
//        inner_local_bxp(i) = inner_sum_xp;
//        outer_local_bxp(i) = outer_sum_xp;
//        inner_local_byp(i) = inner_sum_yp;
//        outer_local_byp(i) = outer_sum_yp;
      });

  Kokkos::fence();

  // =====================================================
  // copy to host for MPI
  // =====================================================

  auto inner_h_bx = Kokkos::create_mirror_view(inner_local_bx);
  auto outer_h_bx = Kokkos::create_mirror_view(outer_local_bx);
  auto inner_h_by = Kokkos::create_mirror_view(inner_local_by);
  auto outer_h_by = Kokkos::create_mirror_view(outer_local_by);
//  auto inner_h_bxp = Kokkos::create_mirror_view(inner_local_bxp);
//  auto outer_h_bxp = Kokkos::create_mirror_view(outer_local_bxp);
//  auto inner_h_byp = Kokkos::create_mirror_view(inner_local_byp);
//  auto outer_h_byp = Kokkos::create_mirror_view(outer_local_byp);

  Kokkos::deep_copy(inner_h_bx, inner_local_bx);
  Kokkos::deep_copy(outer_h_bx, outer_local_bx);
  Kokkos::deep_copy(inner_h_by, inner_local_by);
  Kokkos::deep_copy(outer_h_by, outer_local_by);
//  Kokkos::deep_copy(inner_h_bxp, inner_local_bxp);
//  Kokkos::deep_copy(outer_h_bxp, outer_local_bxp);
//  Kokkos::deep_copy(inner_h_byp, inner_local_byp);
//  Kokkos::deep_copy(outer_h_byp, outer_local_byp);

  // =====================================================
  // MPI reduction
  // =====================================================

  std::vector<Real> inner_global_bx(n1m1+1);
  std::vector<Real> outer_global_bx(n1m1+1);
  std::vector<Real> inner_global_by(n1m1+1);
  std::vector<Real> outer_global_by(n1m1+1);
//  std::vector<Real> inner_global_bxp(n1m1+1);
//  std::vector<Real> outer_global_bxp(n1m1+1);
//  std::vector<Real> inner_global_byp(n1m1+1);
//  std::vector<Real> outer_global_byp(n1m1+1);

  MPI_Allreduce(inner_h_bx.data(),
                inner_global_bx.data(),
                n1m1+1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(outer_h_bx.data(),
                outer_global_bx.data(),
                n1m1+1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(inner_h_by.data(),
                inner_global_by.data(),
                n1m1+1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(outer_h_by.data(),
                outer_global_by.data(),
                n1m1+1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

//    MPI_Allreduce(inner_h_bxp.data(),
//                  inner_global_bxp.data(),
//                  n1m1+1,
//                  MPI_DOUBLE,
//                  MPI_SUM,
//                  MPI_COMM_WORLD);
//
//    MPI_Allreduce(outer_h_bxp.data(),
//                  outer_global_bxp.data(),
//                  n1m1+1,
//                  MPI_DOUBLE,
//                  MPI_SUM,
//                  MPI_COMM_WORLD);
//
//    MPI_Allreduce(inner_h_byp.data(),
//                  inner_global_byp.data(),
//                  n1m1+1,
//                  MPI_DOUBLE,
//                  MPI_SUM,
//                  MPI_COMM_WORLD);
//
//    MPI_Allreduce(outer_h_byp.data(),
//                  outer_global_byp.data(),
//                  n1m1+1,
//                  MPI_DOUBLE,
//                  MPI_SUM,
//                  MPI_COMM_WORLD);

  // =====================================================
  // move global sums back to device
  // =====================================================

  Kokkos::View<Real*> inner_global_d_bx("inner_global_d_bx", n1m1+1);
  Kokkos::View<Real*> outer_global_d_bx("outer_global_d_bx", n1m1+1);
  Kokkos::View<Real*> inner_global_d_by("inner_global_d_by", n1m1+1);
  Kokkos::View<Real*> outer_global_d_by("outer_global_d_by", n1m1+1);
//    Kokkos::View<Real*> inner_global_d_bxp("inner_global_d_bxp", n1m1+1);
//    Kokkos::View<Real*> outer_global_d_bxp("outer_global_d_bxp", n1m1+1);
//    Kokkos::View<Real*> inner_global_d_byp("inner_global_d_byp", n1m1+1);
//    Kokkos::View<Real*> outer_global_d_byp("outer_global_d_byp", n1m1+1);

  auto inner_global_h_bx = Kokkos::create_mirror_view(inner_global_d_bx);
  auto outer_global_h_bx = Kokkos::create_mirror_view(outer_global_d_bx);
  auto inner_global_h_by = Kokkos::create_mirror_view(inner_global_d_by);
  auto outer_global_h_by = Kokkos::create_mirror_view(outer_global_d_by);
//    auto inner_global_h_bxp = Kokkos::create_mirror_view(inner_global_d_bxp);
//    auto outer_global_h_bxp = Kokkos::create_mirror_view(outer_global_d_bxp);
//    auto inner_global_h_byp = Kokkos::create_mirror_view(inner_global_d_byp);
//    auto outer_global_h_byp = Kokkos::create_mirror_view(outer_global_d_byp);

  for (int i=0; i<=n1m1; ++i) {
    inner_global_h_bx(i) = inner_global_bx[i];
    outer_global_h_bx(i) = outer_global_bx[i];
    inner_global_h_by(i) = inner_global_by[i];
    outer_global_h_by(i) = outer_global_by[i];
//      inner_global_h_bxp(i) = inner_global_bxp[i];
//      outer_global_h_bxp(i) = outer_global_bxp[i];
//      inner_global_h_byp(i) = inner_global_byp[i];
//      outer_global_h_byp(i) = outer_global_byp[i];
  }

  Kokkos::deep_copy(inner_global_d_bx, inner_global_h_bx);
  Kokkos::deep_copy(outer_global_d_bx, outer_global_h_bx);
  Kokkos::deep_copy(inner_global_d_by, inner_global_h_by);
  Kokkos::deep_copy(outer_global_d_by, outer_global_h_by);
//    Kokkos::deep_copy(inner_global_d_bxp, inner_global_h_bxp);
//    Kokkos::deep_copy(outer_global_d_bxp, outer_global_h_bxp);
//    Kokkos::deep_copy(inner_global_d_byp, inner_global_h_byp);
//    Kokkos::deep_copy(outer_global_d_byp, outer_global_h_byp);


  // =====================================================
  // refill polar EMFs
  // =====================================================

  Kokkos::parallel_for(
      "polar_refill_b",
      Kokkos::RangePolicy<DevExeSpace>(0, nmb),
      KOKKOS_LAMBDA(const int m) {

        bool do_inner =
            (mb_bcs.d_view(m, BoundaryFace::inner_x2) ==
             BoundaryFlag::polar);

        bool do_outer =
            (mb_bcs.d_view(m, BoundaryFace::outer_x2) ==
             BoundaryFlag::polar);

        if (!(do_inner || do_outer)) return;
          
        Real costi = cos(x2f_(m,js+1));
        Real costo = -cos(x2f_(m,je));

        for (int k=0; k<=n3m1; ++k) {
          Real cosp = cos(x3v_(m,k));
          Real sinp = sin(x3v_(m,k));
          for (int i=0; i<=n1m1; ++i) {

            if (do_inner) {
//              Real inner_bx_avg =
//                  (-inner_global_d_bx(i)/costi + 2.0*inner_global_d_bxp(i)) /
//                  static_cast<Real>(mesh_nx3);
//              Real inner_by_avg =
//                  (-inner_global_d_by(i)/costi + 2.0*inner_global_d_byp(i)) /
//                  static_cast<Real>(mesh_nx3);
                Real inner_bx_avg =
                    (inner_global_d_bx(i)/costi) /
                    static_cast<Real>(mesh_nx3);
                Real inner_by_avg =
                    (inner_global_d_by(i)/costi) /
                    static_cast<Real>(mesh_nx3);
                bx2f(m,k,js,i) = inner_bx_avg*cosp + inner_by_avg*sinp; // 0.5*(bx2f(m,k,js-1,i)+bx2f(m,k,js+1,i)); //
            }

            if (do_outer) {
              Real outer_bx_avg =
                  (outer_global_d_bx(i)/costo) /
                  static_cast<Real>(mesh_nx3);
              Real outer_by_avg =
                  (outer_global_d_by(i)/costo) /
                  static_cast<Real>(mesh_nx3);
                bx2f(m,k,je+1,i) = outer_bx_avg*cosp + outer_by_avg*sinp; // 0.5*(bx2f(m,k,je,i)+bx2f(m,k,je+2,i)); //
            }
          }
        }
      });

  Kokkos::fence();
}
