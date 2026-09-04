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
#include "globals.hpp"
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
  // Re-enabled behind <mesh>/use_polar_average_b (see mesh.hpp).  Must run after the
  // per-block polar case above has set x2f at the pole face, and needs every block's
  // contribution, hence the mesh-wide sum inside it.
  if (ppack->pmesh->use_polar_average_b) PolarAzimuthalAverageBxBy(ppack,b0);
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

void MeshBoundaryValues::PolarAzimuthalAverageBxBy(MeshBlockPack *ppack,
                                                  DvceFaceFld4D<Real> b0) {
  auto bx2f = b0.x2f;
  auto &x2f_ = ppack->pcoord->xx2f;
  auto &x3v_ = ppack->pcoord->x3v;

  auto &indcs = ppack->pmesh->mb_indcs;

  int nx1 = indcs.nx1;
  int nx3 = indcs.nx3;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  const int ncells1 = nx1 + 2*ng;
  const int ncells3 = nx3 + 2*ng;

  const int mesh_nx3 = ppack->pmesh->mesh_indcs.nx3;
  const int nmb = ppack->nmb_thispack;

  auto &mb_bcs = ppack->pmb->mb_bcs;

  // Same rewrite as MHD::PolarAzimuthalAverageEr (see mhd_corner_e.cpp): the sum used to
  // run on ncells1 threads with the MeshBlock and phi loops serial inside each of them,
  // and the refill on one thread per MeshBlock.  Four quantities are carried at once --
  // 0 = inner bx, 1 = outer bx, 2 = inner by, 3 = outer by -- so this needs two
  // allocations where the old code made eight, plus four host mirrors and four vectors.
  //
  // NOTE this function is currently dead: its only call site, in BFieldBCs above, is
  // commented out.  It is cleaned up rather than left as a trap for whoever re-enables it.
  DvceArray3D<Real> part("polar_b_part", 4, nmb, ncells1);
  DvceArray2D<Real> totals("polar_b_totals", 4, ncells1);

  par_for("polar_local_sum_b", DevExeSpace(), 0, nmb-1, 0, ncells1-1,
  KOKKOS_LAMBDA(const int m, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    Real inner_bx = 0.0, outer_bx = 0.0, inner_by = 0.0, outer_by = 0.0;
    // The old RangePolicy stopped at n1m1 = ncells1-1 exclusive, leaving the last radial
    // index at the zero the View was initialised with, while the refill below reads it.
    // Kept exactly so, rather than silently changing what lands in the outermost ghost.
    if ((do_inner || do_outer) && (i < ncells1-1)) {
      for (int k=ks; k<=ke; ++k) {
        Real cosp = cos(x3v_(m,k));
        Real sinp = sin(x3v_(m,k));
        if (do_inner) {
          Real bav = 0.5*(bx2f(m,k,js-1,i)+bx2f(m,k,js+1,i));
          inner_bx += 2.0*bav*cosp;
          inner_by += 2.0*bav*sinp;
        }
        if (do_outer) {
          Real bav = 0.5*(bx2f(m,k,je,i)+bx2f(m,k,je+2,i));
          outer_bx += 2.0*bav*cosp;
          outer_by += 2.0*bav*sinp;
        }
      }
    }
    part(0,m,i) = inner_bx;
    part(1,m,i) = outer_bx;
    part(2,m,i) = inner_by;
    part(3,m,i) = outer_by;
  });

  par_for("polar_combine_b", DevExeSpace(), 0, 3, 0, ncells1-1,
  KOKKOS_LAMBDA(const int q, const int i) {
    Real sum = 0.0;
    for (int m=0; m<nmb; ++m) {
      sum += part(q,m,i);
    }
    totals(q,i) = sum;
  });

#if MPI_PARALLEL_ENABLED
  // Only a distributed run needs this; the old code paid the round trip unconditionally.
  // The Allreduce was also unguarded (breaking non-MPI builds) and hard-coded MPI_DOUBLE.
  if (global_variable::nranks > 1) {
    auto totals_h = Kokkos::create_mirror_view(totals);
    Kokkos::deep_copy(totals_h, totals);
    MPI_Allreduce(MPI_IN_PLACE, totals_h.data(), 4*ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    Kokkos::deep_copy(totals, totals_h);
  }
#endif

  par_for("polar_refill_b", DevExeSpace(), 0, nmb-1, 0, ncells3-1, 0, ncells1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);
    if (!(do_inner || do_outer)) return;

    Real cosp = cos(x3v_(m,k));
    Real sinp = sin(x3v_(m,k));

    if (do_inner) {
      Real costi = cos(x2f_(m,js+1));
      Real bx_avg = (totals(0,i)/costi)/static_cast<Real>(mesh_nx3);
      Real by_avg = (totals(2,i)/costi)/static_cast<Real>(mesh_nx3);
      bx2f(m,k,js,i) = bx_avg*cosp + by_avg*sinp;
    }
    if (do_outer) {
      Real costo = -cos(x2f_(m,je));
      Real bx_avg = (totals(1,i)/costo)/static_cast<Real>(mesh_nx3);
      Real by_avg = (totals(3,i)/costo)/static_cast<Real>(mesh_nx3);
      bx2f(m,k,je+1,i) = bx_avg*cosp + by_avg*sinp;
    }
  });
}
