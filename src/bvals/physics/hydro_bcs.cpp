//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_bcs.cpp
//  \brief

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GravEtotGhosts()
//! \brief <hydro|mhd>/etotgrav: re-reference the potential energy of the ghost cells a
//! built-in BC just filled along direction `dir`.
//!
//! With etotgrav the conserved energy is E = e + KE (+ emag) + rho*Phi.  The reflect,
//! outflow and diode BCs COPY u0 from an active cell, so the ghost inherits rho*Phi of
//! the SOURCE cell, and ConToPrim then removes rho*Phi of the GHOST: the ghost's internal
//! energy came out wrong by rho*(Phi_src - Phi_ghost), i.e. rho g dz per ghost layer
//! wherever the potential varies across the face (a bottom wall is heated, an outflow top
//! can go negative).  After the copy the ghost's density is the source's, so
//!     E_ghost += rho_ghost*(Phi_ghost - Phi_src)
//! makes it the source's e + KE plus the ghost's own potential energy.  Inflow ghosts
//! are set from u_in, which carries no potential energy: they get + rho_ghost*Phi_ghost.
//! Vacuum ghosts are zero (rho = 0) and need nothing.  Where the potential does not vary
//! across the face the correction is rho*0 = +0.0 and E is unchanged bit for bit.

namespace {
void GravEtotGhosts(MeshBlockPack *ppack, DvceArray5D<Real> u0, DvceArray4D<Real> phi,
                    const int dir) {
  auto &indcs = ppack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  auto &mb_bcs = ppack->pmb->mb_bcs;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  const int nmb = ppack->nmb_thispack;
  int lo, hi, na, nb, fin, fout;
  if (dir == 1) {
    lo = indcs.is; hi = indcs.ie; na = n3; nb = n2;
    fin = BoundaryFace::inner_x1; fout = BoundaryFace::outer_x1;
  } else if (dir == 2) {
    lo = indcs.js; hi = indcs.je; na = n3; nb = n1;
    fin = BoundaryFace::inner_x2; fout = BoundaryFace::outer_x2;
  } else {
    lo = indcs.ks; hi = indcs.ke; na = n2; nb = n1;
    fin = BoundaryFace::inner_x3; fout = BoundaryFace::outer_x3;
  }
  par_for("hydrobc_etg", DevExeSpace(), 0, (nmb-1), 0, (na-1), 0, (nb-1),
  KOKKOS_LAMBDA(int m, int a, int b) {
    // (a, b) = (k, j) for dir 1, (k, i) for dir 2, (j, i) for dir 3; `x` is the index
    // along dir
    auto fix = [&](const int xg, const int xs) {
      const int k = (dir == 3) ? xg : a;
      const int j = (dir == 2) ? xg : ((dir == 3) ? a : b);
      const int i = (dir == 1) ? xg : b;
      Real dphi = phi(m,k,j,i);
      if (xs >= 0) {
        const int ks = (dir == 3) ? xs : k;
        const int js = (dir == 2) ? xs : j;
        const int is = (dir == 1) ? xs : i;
        dphi -= phi(m,ks,js,is);
      }
      u0(m,IEN,k,j,i) += u0(m,IDN,k,j,i)*dphi;
    };
    for (int face = 0; face < 2; ++face) {
      const BoundaryFlag flag = mb_bcs.d_view(m, (face == 0) ? fin : fout);
      for (int g=0; g<ng; ++g) {
        const int xg = (face == 0) ? (lo - g - 1) : (hi + g + 1);
        switch (flag) {
          case BoundaryFlag::reflect:
            fix(xg, (face == 0) ? (lo + g) : (hi - g));
            break;
          case BoundaryFlag::outflow:
          case BoundaryFlag::diode:
            fix(xg, (face == 0) ? lo : hi);
            break;
          case BoundaryFlag::inflow:
            fix(xg, -1);
            break;
          default:
            break;
        }
      }
    }
  });
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \!fn void BoundaryValues::HydroBCs()
//! \brief Apply physical boundary conditions for all Hydro variables at faces of MB which
//! are at the edge of the computational domain

void MeshBoundaryValues::HydroBCs(MeshBlockPack *ppack, DualArray2D<Real> u_in,
                                  DvceArray5D<Real> u0, DvceArray4D<Real> phi) {
  const bool etg = (phi.size() > 0);
  // loop over all MeshBlocks in this MeshBlockPack
  auto &pm = ppack->pmesh;
  auto &indcs = ppack->pmesh->mb_indcs;
  int &ng = indcs.ng;
  auto &mb_bcs = ppack->pmb->mb_bcs;

  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int nvar = u0.extent_int(1);  // TODO(@user): 2nd index from L of in array must be NVAR
  int nmb = ppack->nmb_thispack;

  // only apply BCs if not (periodic) or (shear_periodic)
  if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic &&
      pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::shear_periodic) {
    int &is = indcs.is;
    int &ie = indcs.ie;
    par_for("hydrobc_x1", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int n, int k, int j) {
      // apply physical boundaries to inner_x1
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,is-i-1) = -u0(m,n,k,j,is+i);
            } else {
              u0(m,n,k,j,is-i-1) =  u0(m,n,k,j,is+i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = u0(m,n,k,j,is);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = u_in.d_view(n,BoundaryFace::inner_x1);
          }
          break;
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,is-i-1) = fmin(0.0,u0(m,n,k,j,is));
            } else {
              u0(m,n  ,k,j,is-i-1) = u0(m,n,k,j,is);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,is-i-1) = 0.0;
          }
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x1
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x1)) {
        case BoundaryFlag::reflect:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {  // reflect 1-velocity
              u0(m,n,k,j,ie+i+1) = -u0(m,n,k,j,ie-i);
            } else {
              u0(m,n,k,j,ie+i+1) =  u0(m,n,k,j,ie-i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = u0(m,n,k,j,ie);
          }
          break;
        case BoundaryFlag::inflow:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = u_in.d_view(n,BoundaryFace::outer_x1);
          }
          break;
        case BoundaryFlag::diode:
          for (int i=0; i<ng; ++i) {
            if (n==(IVX)) {
              u0(m,n,k,j,ie+i+1) = fmax(0.0,u0(m,n,k,j,ie));
            } else {
              u0(m,n  ,k,j,ie+i+1) = u0(m,n,k,j,ie);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int i=0; i<ng; ++i) {
            u0(m,n,k,j,ie+i+1) = 0.0;
          }
          break;
        default:
          break;
      }
    });
    if (etg) GravEtotGhosts(ppack, u0, phi, 1);
  }

  if (pm->one_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic) {
    int &js = indcs.js;
    int &je = indcs.je;
    par_for("hydrobc_x2", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int n, int k, int i) {
      // apply physical boundaries to inner_x2
      switch (mb_bcs.d_view(m,BoundaryFace::inner_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {  // reflect 2-velocity
              u0(m,n,k,js-j-1,i) = -u0(m,n,k,js+j,i);
            } else {
              u0(m,n,k,js-j-1,i) =  u0(m,n,k,js+j,i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = u0(m,n,k,js,i);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = u_in.d_view(n,BoundaryFace::inner_x2);
          }
          break;
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {
              u0(m,n,k,js-j-1,i) = fmin(0.0,u0(m,n,k,js,i));
            } else {
              u0(m,n,k,js-j-1,i) = u0(m,n,k,js,i);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,js-j-1,i) = 0.0;
          }
          break;
        default:
          break;
      }

      // apply physical boundaries to outer_x2
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x2)) {
        case BoundaryFlag::reflect:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {  // reflect 2-velocity
              u0(m,n,k,je+j+1,i) = -u0(m,n,k,je-j,i);
            } else {
              u0(m,n,k,je+j+1,i) =  u0(m,n,k,je-j,i);
            }
          }
          break;
        case BoundaryFlag::outflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = u0(m,n,k,je,i);
          }
          break;
        case BoundaryFlag::inflow:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = u_in.d_view(n,BoundaryFace::outer_x2);
          }
          break;
        case BoundaryFlag::diode:
          for (int j=0; j<ng; ++j) {
            if (n==(IVY)) {
              u0(m,n,k,je+j+1,i) = fmax(0.0,u0(m,n,k,je,i));
            } else {
              u0(m,n,k,je+j+1,i) = u0(m,n,k,je,i);
            }
          }
          break;
        case BoundaryFlag::vacuum:
          for (int j=0; j<ng; ++j) {
            u0(m,n,k,je+j+1,i) = 0.0;
          }
          break;
        default:
          break;
      }
    });
    if (etg) GravEtotGhosts(ppack, u0, phi, 2);
  }
  if (pm->two_d) return;

  // only apply BCs if not periodic
  if (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic) return;
  int &ks = indcs.ks;
  int &ke = indcs.ke;
  par_for("hydrobc_x3", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int n, int j, int i) {
    // apply physical boundaries to inner_x3
    switch (mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {  // reflect 3-velocity
            u0(m,n,ks-k-1,j,i) = -u0(m,n,ks+k,j,i);
          } else {
            u0(m,n,ks-k-1,j,i) =  u0(m,n,ks+k,j,i);
          }
        }
        break;
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = u0(m,n,ks,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = u_in.d_view(n,BoundaryFace::inner_x3);
        }
        break;
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {
            u0(m,n,ks-k-1,j,i) = fmin(0.0,u0(m,n,ks,j,i));
          } else {
            u0(m,n,ks-k-1,j,i) = u0(m,n,ks,j,i);
          }
        }
        break;
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ks-k-1,j,i) = 0.0;
        }
        break;
      default:
        break;
    }

    // apply physical boundaries to outer_x3
    switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
      case BoundaryFlag::reflect:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {  // reflect 3-velocity
            u0(m,n,ke+k+1,j,i) = -u0(m,n,ke-k,j,i);
          } else {
            u0(m,n,ke+k+1,j,i) =  u0(m,n,ke-k,j,i);
          }
        }
        break;
      case BoundaryFlag::outflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = u0(m,n,ke,j,i);
        }
        break;
      case BoundaryFlag::inflow:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = u_in.d_view(n,BoundaryFace::outer_x3);
        }
        break;
      case BoundaryFlag::diode:
        for (int k=0; k<ng; ++k) {
          if (n==(IVZ)) {
            u0(m,n,ke+k+1,j,i) = fmax(0.0,u0(m,n,ke,j,i));
          } else {
            u0(m,n,ke+k+1,j,i) = u0(m,n,ke,j,i);
          }
        }
        break;
      case BoundaryFlag::vacuum:
        for (int k=0; k<ng; ++k) {
          u0(m,n,ke+k+1,j,i) = 0.0;
        }
        break;
      default:
        break;
    }
  });
  if (etg) GravEtotGhosts(ppack, u0, phi, 3);

  return;
}
