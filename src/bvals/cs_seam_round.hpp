#ifndef BVALS_CS_SEAM_ROUND_HPP_
#define BVALS_CS_SEAM_ROUND_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
//========================================================================================
//! \file cs_seam_round.hpp
//! \brief fc2vtx-1103: TIE-SAFE position -> index rounding for the cubed-sphere seam
//! stencil tables (bvals_cc.cpp CC resample, bvals_fc.cpp order-2 tables and order-3
//! windows).  A seam-mapped position computed from a MeshBlock's own coordinates can sit
//! EXACTLY on a rounding boundary by symmetry -- at a panel-end face the map shifts by
//! the normal offset, so a ghost cell centre maps to a half-integer -- and plain
//! floor(pos + 0.5) then follows the last bit of pos, i.e. the block layout.  The fixed
//! bias makes every layout break such a tie the same way; any argument farther than
//! kSeamTie from a boundary rounds exactly as before (bitwise).
#include <cmath>
#include "athena.hpp"

namespace cs_seam {
constexpr Real kSeamTie = 1.0e-9;
//! nearest integer, ties up (layout-independent)
KOKKOS_FORCEINLINE_FUNCTION
int SeamRound(const Real x) { return static_cast<int>(floor(x + 0.5 + kSeamTie)); }
//! layoutfix-1106 <mesh>/cs_seam_pos_global: sc + d with the integer part at its magnitude in
//! the PANEL (s = the block's first cell in the panel): sc + d rounds d to the ulp of the
//! block-local index, so the weights differ in the last bits between MeshBlock layouts; the
//! removal of s is exact.  s = 0 is sc + d, bitwise
KOKKOS_FORCEINLINE_FUNCTION
Real SeamPos(const int sc, const int s, const Real d) {
  return (static_cast<Real>(sc + s) + d) - static_cast<Real>(s);
}
//! floor, an exact integer stays itself (layout-independent)
KOKKOS_FORCEINLINE_FUNCTION
int SeamFloor(const Real x) { return static_cast<int>(floor(x + kSeamTie)); }
}  // namespace cs_seam
#endif  // BVALS_CS_SEAM_ROUND_HPP_
