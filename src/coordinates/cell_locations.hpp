#ifndef COORDINATES_CELL_LOCATIONS_HPP_
#define COORDINATES_CELL_LOCATIONS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cell_locationss.hpp
//  \brief functions to compute locations on a uniform Cartesian grid
// They provide functionality of the Coordinates class in the C++ version of the code.
// Very similar to cc_pos.c function in C version of the code (Athena4.2)
// Not incoporated in Coordinates class so that they can be used anywhere (for exmaple
// to compute locations of MeshBlocks in Mesh).

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn void LeftEdgeX()
// returns x-posn of left edge of i^th cell where index range [0,N] maps to [xmin,xmax]
// returns ghost cell posn if i outside range [0,N] (e.g. i=-1 is x-posn of first ghost
// cell). Averages linear interpolation from each side to symmetrize r.o. error

KOKKOS_INLINE_FUNCTION
static Real LeftEdgeX(int ith, int n, Real xmin, Real xmax) {
  Real x = (static_cast<Real>(ith)) / (static_cast<Real>(n));
  return (x*xmax - x*xmin) - (0.5*xmax - 0.5*xmin) + (0.5*xmin + 0.5*xmax);
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenterX()
// returns cell-center posn of i^th cell where index range [0,N] maps to [xmin,xmax]
// returns ghost cell posn if i outside range [0,N] (e.g. i=-1 is cc-posn of first ghost
// cell). Averages linear interpolation from each side to symmetrize r.o. error

KOKKOS_INLINE_FUNCTION
static Real CellCenterX(int ith, int n, Real xmin, Real xmax) {
  Real x = (static_cast<Real>(ith) + 0.5) / (static_cast<Real>(n));
  return (x*xmax - x*xmin) - (0.5*xmax - 0.5*xmin) + (0.5*xmin + 0.5*xmax);
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenterIndex()
// returns i-index of cell containing x position

// TODO(@user): set trap if out-of-range

KOKKOS_INLINE_FUNCTION
static int CellCenterIndex(Real x, int n, Real xmin, Real xmax) {
  return static_cast<int>(((x-xmin)/(xmax-xmin))*static_cast<Real>(n));
}

//----------------------------------------------------------------------------------------
//! \fn Real CellCenteredRadialFld()
//! \brief Cell-centred value of a radial face-centred field by linear interpolation to
//! x1v, which on a stretched radial grid is not the midpoint of the two faces.

KOKKOS_INLINE_FUNCTION
static Real CellCenteredRadialFld(const Real bl, const Real br, const Real x1f_l,
                                  const Real x1f_r, const Real x1v) {
  const Real idx = 1.0/(x1f_r - x1f_l);
  return ((x1f_r - x1v)*idx)*bl + ((x1v - x1f_l)*idx)*br;
}

//----------------------------------------------------------------------------------------
//! \fn Real CellCenteredRadialFldWall()
//! \brief CellCenteredRadialFld with mirror-consistent weights behind a reflecting
//! radial wall: a ghost cell uses the exchanged weights of its interior image cell, so
//! the ghost value is exactly minus the image value.
template <typename FaceView, typename PosView>
KOKKOS_INLINE_FUNCTION
Real CellCenteredRadialFldWall(const FaceView &b1f, const PosView &x1f,
                               const PosView &x1v, const int m, const int k,
                               const int j, const int i, const int is, const int ie,
                               const bool refl_in, const bool refl_out) {
  int ii = i;                 // cell whose weights are used
  bool mirror = false;
  if (refl_in && i < is) {
    ii = 2*is - 1 - i;  mirror = true;
  } else if (refl_out && i > ie) {
    ii = 2*ie + 1 - i;  mirror = true;
  }
  const Real idx = 1.0/(x1f(m,ii+1) - x1f(m,ii));
  const Real wl = (x1f(m,ii+1) - x1v(m,ii))*idx;    // weight of the image's inner face
  const Real wr = (x1v(m,ii) - x1f(m,ii))*idx;      // weight of the image's outer face
  if (!mirror) return wl*b1f(m,k,j,i) + wr*b1f(m,k,j,i+1);
  // ghost inner face is the mirror of the image outer face, and vice versa
  return wr*b1f(m,k,j,i) + wl*b1f(m,k,j,i+1);
}

#endif // COORDINATES_CELL_LOCATIONS_HPP_
