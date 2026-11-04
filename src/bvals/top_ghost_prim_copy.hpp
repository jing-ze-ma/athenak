#ifndef BVALS_TOP_GHOST_PRIM_COPY_HPP_
#define BVALS_TOP_GHOST_PRIM_COPY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file top_ghost_prim_copy.hpp
//! \brief wbplm-1103 <hydro|mhd>/top_ghost_prim_copy (default false = bitwise; acts only
//! with the one-way zero-gradient top, diode_ox1 = problem/top_bc outflow1): after
//! ConsToPrim, the outer-x1 user-BC ghost cells take the internal energy, temperature and
//! derived EOS quantities (p, Gamma_1, ...) of the last active cell BITWISE.
//!
//! outflow1 means a zero-gradient copy of rho and e, but the ghost conserved energy is
//! assembled as e + KE + rho phi_ghost (+ B_ghost^2/2) and the inversion takes those terms
//! back off, so the ghost pressure came back as the cell's plus a few ulp of either sign.
//! The x1 well-balanced PLM tests the sign of sR = (p_ghost - p_ie)/dx and switched the top
//! face between two reconstructions on that roundoff (wbplm-1103 notes).  With the copy
//! sR is exactly 0 and the choice is deterministic.  rho and the velocities are already
//! bitwise copies (same rho, copied/clipped momenta) and are not touched.
#include "athena.hpp"
#include "mesh/mesh.hpp"

inline void TopGhostPrimCopy(MeshBlockPack *pp, const DvceArray5D<Real> &w0,
                             const DvceArray5D<Real> &wder, const DvceArray4D<Real> &wtemp,
                             const bool gen) {
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
