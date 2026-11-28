#ifndef RAD_VET_RAYS_SC_CART_SUM_HPP_
#define RAD_VET_RAYS_SC_CART_SUM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_sum.hpp
//! \brief the Cartesian slab's per-cell face sums (VetScFaceSum, one thread per cell and
//! component), every sum added in ray order.
//!
//!  called from: rays/sc_cart_setup.cpp

#include "athena.hpp"
#include "rad_vet/rad_vet.hpp"

namespace radvet {

// kinds of a per-cell face sum, as passed to a writer w(kind, u, n, value) (u = 0 up,
// 1 down): F2 to the cell, FN to its upwind cell
enum VetScSum {VSS_F2, VSS_FN};

//----------------------------------------------------------------------------------------
//! \fn VetScFaceSum
//! \brief the face sums n = c, c + 4, c + 8, c + 12 of component c (1, mu1, mu2, mu3) of
//! one cell over rays [0, nray): ray r with active(r) (up rays m1 > 0 to u = 0) and
//! face-sum scalar g(q, r), added in ray order; written to w (up f2, down f2, up fn,
//! down fn of each n)
template <class A, class FA, class FG, class W>
KOKKOS_INLINE_FUNCTION
void VetScFaceSum(const int c, const A &ang, const int nray, const FA &active,
                  const FG &g, const W &w) {
  Real f2[2][4] = {}, fn[2][4] = {};
  for (int r = 0; r < nray; ++r) {
    if (!active(r)) continue;
    const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2);
    const bool up = (m1 > 0.0);
    const int u = up ? 0 : 1;
    auto cm = [&](const Real q) {
      return (c == 0) ? q : ((c == 1) ? q*m1 : ((c == 2) ? q*m2 : q*m3));
    };
    f2[u][1] += cm(g(0, r));
    fn[u][0] += cm(g(1, r));
    fn[u][3] += cm(g(2, r));
    if (up) {f2[0][2] += cm(g(3, r));} else {fn[1][2] += cm(g(3, r));}
    f2[u][1] += cm(g(4, r));
  }
  for (int b = 0; b < 4; ++b) {
    w(VSS_F2, 0, 4*b + c, f2[0][b]);
    w(VSS_F2, 1, 4*b + c, f2[1][b]);
    w(VSS_FN, 0, 4*b + c, fn[0][b]);
    w(VSS_FN, 1, 4*b + c, fn[1][b]);
  }
}

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_CART_SUM_HPP_
