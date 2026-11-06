#ifndef RAD_M1_RAD_M1_VGDLG_HPP_
#define RAD_M1_RAD_M1_VGDLG_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vgdlg.hpp
//! \brief layinv-1105 <mesh>/layout_invariant: the (theta, phi) positions and spacings
//! the vet_gd sweep uses on a spherical-polar mesh, from the GLOBAL cell index.
//!
//! The sweep places the band rows/columns at x2min_b + (j - js + 0.5) dx2_b of its block
//! b (uniform in index space).  The block bounds and dx2_b = (x2max_b - x2min_b)/nx2 are
//! rounded per block, so the same global cell gets positions differing in the last bit
//! between MeshBlock layouts; the ray geometry (and a tangent cell's branch test) follows.
//! With on: x2a + (J + 0.5) dx2 with the global row J and the MESH spacing (the same
//! numbers for every layout); off: the old per-block expressions, bitwise.

#include "athena.hpp"

namespace radm1 {

struct VgdLG {
  bool on = false;
  Real x2a = 0.0, dx2 = 1.0, x3a = 0.0, dx3 = 1.0;
  Kokkos::View<int**, DevMemSpace> o;            // (m, 0|1): global j|k of the first cell
  Kokkos::View<int**, DevMemSpace>::HostMirror oh;

  template <class TB> KOKKOS_INLINE_FUNCTION
  Real Th(const TB &mb, const int m, const int j, const int js) const {
    return on ? x2a + (static_cast<Real>(o(m,0) + j - js) + 0.5)*dx2
              : mb.d_view(m).x2min + (j - js + 0.5)*mb.d_view(m).dx2;
  }
  template <class TB> KOKKOS_INLINE_FUNCTION
  Real Ph(const TB &mb, const int m, const int k, const int ks) const {
    return on ? x3a + (static_cast<Real>(o(m,1) + k - ks) + 0.5)*dx3
              : mb.d_view(m).x3min + (k - ks + 0.5)*mb.d_view(m).dx3;
  }
  template <class TB> KOKKOS_INLINE_FUNCTION
  Real Dx2(const TB &mb, const int m) const {return on ? dx2 : mb.d_view(m).dx2;}
  template <class TB> KOKKOS_INLINE_FUNCTION
  Real Dx3(const TB &mb, const int m) const {return on ? dx3 : mb.d_view(m).dx3;}
  // js + floor((t - x2min_b)/dx2_b + 0.5): the band row (origin js) nearest the face t
  template <class TB> KOKKOS_INLINE_FUNCTION
  Real Jw(const TB &mb, const int m, const int js, const Real t) const {
    return on ? (js - o(m,0)) + floor((t - x2a)/dx2 + 0.5)
              : js + floor((t - mb.d_view(m).x2min)/mb.d_view(m).dx2 + 0.5);
  }
  // the band-row / column offset to the global index (0 when off)
  KOKKOS_INLINE_FUNCTION int Sj(const int m, const int js) const {return on ? o(m,0) - js : 0;}
  KOKKOS_INLINE_FUNCTION int Sk(const int m, const int ks) const {return on ? o(m,1) - ks : 0;}
  // the fractional index j + d with the integer part at its global magnitude (exact
  // removal): the same rounding of d on every layout; s = 0 gives j + d, bitwise
  KOKKOS_INLINE_FUNCTION static Real Fs(const int j, const int s, const Real d) {
    return (static_cast<Real>(j + s) + d) - s;
  }
  // host versions (mb.h_view)
  template <class TB>
  Real ThH(const TB &mb, const int m, const int j, const int js) const {
    return on ? x2a + (static_cast<Real>(oh(m,0) + j - js) + 0.5)*dx2
              : mb.h_view(m).x2min + (j - js + 0.5)*mb.h_view(m).dx2;
  }
  template <class TB>
  Real PhH(const TB &mb, const int m, const int k, const int ks) const {
    return on ? x3a + (static_cast<Real>(oh(m,1) + k - ks) + 0.5)*dx3
              : mb.h_view(m).x3min + (k - ks + 0.5)*mb.h_view(m).dx3;
  }
  template <class TB>
  Real Dx2H(const TB &mb, const int m) const {return on ? dx2 : mb.h_view(m).dx2;}
  template <class TB>
  Real Dx3H(const TB &mb, const int m) const {return on ? dx3 : mb.h_view(m).dx3;}
};

}  // namespace radm1

#endif  // RAD_M1_RAD_M1_VGDLG_HPP_
