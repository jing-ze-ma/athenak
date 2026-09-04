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
//! \brief the cell-centred value of a face-centred field in the RADIAL direction, as a
//! linear interpolation to the cell centre rather than the plain average of its faces.
//!
//! WHY THIS IS NOT 0.5*(bl + br).  On a STRETCHED radial grid the cell centre is not the
//! midpoint of its two faces: x1v is S(midpoint of the uniform cell) while the faces are
//! S(left) and S(right), and S is nonlinear, so S(mid) != 0.5*(S(l) + S(r)).  The plain
//! average therefore returns the field at the FACE MIDPOINT, not at x1v.  Spherical polar
//! has always used the weighted form; the cubed sphere used the plain average even though
//! its x1 is radial and stretched in exactly the same way.
//!
//! HOW BIG.  On the production hot-Jupiter grid (128 radial cells, 4.07x spread in cell
//! width, the 4-coefficient polynomial stretch) the weight departs from 0.5 by at most
//! 4.4e-3 and on average 1.9e-3.  This is a second-order accuracy fix, not a bug fix --
//! it will not move a dynamical result.  The reason to make it is CONSISTENCY: everything
//! that rebuilds bcc must agree cell for cell, or the difference reappears as an energy
//! mismatch that the C2P floors then bake into u.e.
//!
//! Reduces to exactly 0.5*(bl + br) only up to round-off on a uniform grid, so it is
//! applied ONLY where x1 really is a stretched radial coordinate; everywhere else keeps
//! the plain average and its bit-for-bit answers.

KOKKOS_INLINE_FUNCTION
static Real CellCenteredRadialFld(const Real bl, const Real br, const Real x1f_l,
                                  const Real x1f_r, const Real x1v) {
  const Real idx = 1.0/(x1f_r - x1f_l);
  return ((x1f_r - x1v)*idx)*bl + ((x1v - x1f_l)*idx)*br;
}

#endif // COORDINATES_CELL_LOCATIONS_HPP_
