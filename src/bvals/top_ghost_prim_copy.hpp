#ifndef BVALS_TOP_GHOST_PRIM_COPY_HPP_
#define BVALS_TOP_GHOST_PRIM_COPY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file top_ghost_prim_copy.hpp
//! \brief <hydro|mhd>/top_ghost_prim_copy: with the zero-gradient outer-x1 user boundary
//! (problem/top_bc = outflow1), copies the internal energy, temperature and derived EOS
//! quantities of the last active cell exactly into the ghost cells after ConsToPrim, so
//! the ghost pressure equals the cell pressure (deterministic well-balanced PLM choice).

#include "athena.hpp"
#include "mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn void TopGhostPrimCopy
//! \brief Copies e (or p), T and the derived EOS variables of cell ie into the outer-x1
//! ghost cells.

inline void TopGhostPrimCopy(MeshBlockPack *pp, const DvceArray5D<Real> &w0,
                             const DvceArray5D<Real> &wder,
                             const DvceArray4D<Real> &wtemp, const bool gen) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int ie = indcs.ie, ng = indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb1 = pp->nmb_thispack - 1;
  const int nder = (gen && wder.is_allocated()) ? wder.extent_int(1) : 0;
  const bool hast = gen && wtemp.is_allocated();
  auto &mb_bcs = pp->pmb->mb_bcs;
  par_for("topghostprim", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (mb_bcs.d_view(m, BoundaryFace::outer_x1) != BoundaryFlag::user) return;
    for (int i = 1; i <= ng; ++i) {
      w0(m,IEN,k,j,ie+i) = w0(m,IEN,k,j,ie);
      if (hast) wtemp(m,k,j,ie+i) = wtemp(m,k,j,ie);
      for (int n = 0; n < nder; ++n) wder(m,n,k,j,ie+i) = wder(m,n,k,j,ie);
    }
  });
}
#endif  // BVALS_TOP_GHOST_PRIM_COPY_HPP_
