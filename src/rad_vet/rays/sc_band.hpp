#ifndef RAD_VET_RAYS_SC_BAND_HPP_
#define RAD_VET_RAYS_SC_BAND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_band.hpp
//! \brief (theta, phi) cell positions and spacings of the spherical-polar vet_gd sweep.
//! <mesh>/layout_invariant on: x2a + (J + 0.5) dx2 from the global index J and the mesh
//! spacing (identical for every MeshBlock layout); off: per-block x2min_b, dx2_b.

#include "athena.hpp"

namespace radvet {

struct ScLG {
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
  KOKKOS_INLINE_FUNCTION int Sj(const int m, const int js) const {
    return on ? o(m,0) - js : 0;
  }
  KOKKOS_INLINE_FUNCTION int Sk(const int m, const int ks) const {
    return on ? o(m,1) - ks : 0;
  }
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

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_BAND_HPP_
