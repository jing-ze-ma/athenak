#ifndef BVALS_CS_SEAM_ROUND_HPP_
#define BVALS_CS_SEAM_ROUND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_seam_round.hpp
//! \brief Position-to-index rounding for the cubed-sphere seam stencil tables.
//! Seam-mapped
//! positions can fall exactly on a rounding boundary (half-integers); a small fixed bias
//! breaks such ties the same way for every MeshBlock layout.

#include <cmath>

#include "athena.hpp"

namespace cs_seam {
constexpr Real kSeamTie = 1.0e-9;

//----------------------------------------------------------------------------------------
//! \fn int cs_seam::SeamRound
//! \brief Nearest integer, ties rounded up.

KOKKOS_FORCEINLINE_FUNCTION
int SeamRound(const Real x) { return static_cast<int>(floor(x + 0.5 + kSeamTie)); }

//----------------------------------------------------------------------------------------
//! \fn Real cs_seam::SeamPos
//! \brief Seam position sc + d with the integer part taken at its magnitude in the panel
//! (s = first cell of the block in the panel), so the rounding of d does not depend on
//! the block layout.

KOKKOS_FORCEINLINE_FUNCTION
Real SeamPos(const int sc, const int s, const Real d) {
  return (static_cast<Real>(sc + s) + d) - static_cast<Real>(s);
}

//----------------------------------------------------------------------------------------
//! \fn int cs_seam::SeamFloor
//! \brief Floor, with exact integers kept.

KOKKOS_FORCEINLINE_FUNCTION
int SeamFloor(const Real x) { return static_cast<int>(floor(x + kSeamTie)); }
}  // namespace cs_seam
#endif  // BVALS_CS_SEAM_ROUND_HPP_
