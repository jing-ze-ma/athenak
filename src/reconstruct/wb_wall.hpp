#ifndef RECONSTRUCT_WB_WALL_HPP_
#define RECONSTRUCT_WB_WALL_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_wall.hpp
//! \brief <hydro|mhd>/wall_wb_onesided: the wall cell of a closed x1 wall (reflect, or a
//! user face under wall_closed_ix1) reconstructs its well-balanced deviation with the
//! one-sided slope towards the interior at both faces.
//!
//! The ghosts of a closed wall carry no consistent WB deviation, so the limited slope of
//! the first active cell is flattened and the cell sees an unbalanced force, held by a
//! stationary odd-even v_x1 pattern.  The one-sided slope reconstructs a linear deviation
//! exactly, so the wall cell is balanced like an interior cell.
//!
//! Shared by Hydro::Fluxes and MHD::Fluxes on every geometry.  Channels: e, rho (wb_rho)
//! and the general-EOS pressure (nder > 0), from the x1 WB cache wbq0; a channel keeps
//! its faces where a one-sided face value is not positive.

#include "athena.hpp"
#include "hydro/hydro.hpp"   // WBReadCache, WBVar, IDPR

namespace wbwall {

//! One channel: background at the wall face bw and at the other face bo, deviation of
//! the cell q1c and of its interior neighbour q1n, centre-to-neighbour distance dn,
//! centre-to-wall-face hw and centre-to-other-face ho (all > 0).  Returns false (faces
//! untouched) when a face value would not be positive.
KOKKOS_INLINE_FUNCTION
bool OneSided(const Real bw, const Real bo, const Real q1c, const Real q1n, const Real dn,
              const Real hw, const Real ho, Real &fw, Real &fo) {
  const Real s = (q1n - q1c)/dn;     // deviation gradient pointing into the domain
  fw = bw + q1c - s*hw;
  fo = bo + q1c + s*ho;
  return (fw > 0.0 && fo > 0.0);
}

//! The wall cell ic of an x1 wall (inner: ic = is, wall face ic, neighbour ic+1; outer:
//! ic = ie, wall face ic+1, neighbour ic-1).  wl(n,i)/wr(n,i) = left/right states of face
//! i (the sweep's convention: cell i writes wr(n,i) and wl(n,i+1)); dl/dr the same for
//! the derived variables.  Called from ONE thread of the team (Kokkos::single).
template <typename ScrT>
KOKKOS_INLINE_FUNCTION
void WallCellX1(const bool inner, const int m, const int k, const int j, const int ic,
                const DvceArray5D<Real> &w, const DvceArray5D<Real> &wd,
                const DvceArray5D<Real> &wbq, const DvceArray2D<Real> &xv,
                const DvceArray2D<Real> &xf, const bool wb_rho, const int nder,
                ScrT &wl, ScrT &wr, ScrT &dl, ScrT &dr) {
  const int in = inner ? ic + 1 : ic - 1;
  const Real hl = xv(m,ic) - xf(m,ic), hr = xf(m,ic+1) - xv(m,ic);
  const Real dn = inner ? (xv(m,in) - xv(m,ic)) : (xv(m,ic) - xv(m,in));
  for (int c=0; c<3; ++c) {
    if (c == 0 && !wb_rho) continue;
    if (c == 2 && nder == 0) continue;
    const WBVar wv = (c == 0) ? WBVar::wb_dens : ((c == 1) ? WBVar::wb_eint : WBVar::wb_pres);
    Real bm, bmh, bc, bph, bp;
    WBReadCache(wbq, wv, m, k, j, ic, bm, bmh, bc, bph, bp);
    const int n = (c == 0) ? IDN : IEN;
    const Real qc = (c == 2) ? wd(m,IDPR,k,j,ic) : w(m,n,k,j,ic);
    const Real qn = (c == 2) ? wd(m,IDPR,k,j,in) : w(m,n,k,j,in);
    Real fw, fo;
    if (inner) {
      if (!OneSided(bmh, bph, qc - bc, qn - bp, dn, hl, hr, fw, fo)) continue;
      if (c == 2) {dr(IDPR,ic) = fw; dl(IDPR,ic+1) = fo;} else {wr(n,ic) = fw; wl(n,ic+1) = fo;}
    } else {
      if (!OneSided(bph, bmh, qc - bc, qn - bm, dn, hr, hl, fw, fo)) continue;
      if (c == 2) {dl(IDPR,ic+1) = fw; dr(IDPR,ic) = fo;} else {wl(n,ic+1) = fw; wr(n,ic) = fo;}
    }
  }
}

}  // namespace wbwall
#endif  // RECONSTRUCT_WB_WALL_HPP_
