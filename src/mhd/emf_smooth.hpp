#ifndef MHD_EMF_SMOOTH_HPP_
#define MHD_EMF_SMOOTH_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emf_smooth.hpp
//! \brief the GS07 upwind weight (EmfSmoothW, lhlld2-1102), shared by MHD::CornerE and
//! the cube-vertex EMF of MHD::CsVertexEmf3 (cornerfix-1103).  Moved verbatim from
//! mhd_corner_e.cpp.

#include "athena.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn Real EmfSmoothW(fd, eps, rhoe)
//! \brief lhlld2-1102, <mhd>/emf_smooth = eps (read only when named; default off =
//! the GS07 sign switch, bitwise).  Weight of the LEFT (lower-index) upwind choice in the
//! GS07 corner EMF, as a smooth function of the face mass flux fd:
//!   w = 0.5*(1 + fd/sqrt(fd^2 + delta^2)),  delta = eps*sqrt(rho*e)  (~ eps*rho*c_s),
//! instead of the step function w = (fd >= 0).  The step makes the semi-discrete
//! right-hand side DISCONTINUOUS wherever a face mass flux changes sign (any flow at
//! rest or turning, and every low-Mach flow, where fd is round-off-sized noise about 0):
//! the two upwind choices differ by O(dx^2 d^2E/dx^2), so a Runge-Kutta stage that lands
//! on the other side of the switch than its predecessor costs an O(dt * jump) local
//! error, and the global temporal error is FIRST order at fixed grid (in-plane-B
//! magnetized Gresho vortex, M = 0.01: B order 0.8-1.3, v 1.2).  eps = 0 is the GS07
//! rule for a zero mass flux (the average of both sides) applied to every face.
//! For |fd| >> delta the weight is 1 or 0 to O(delta^2/fd^2): GS07 is recovered in any
//! flow that is not slow compared with eps*c_s.

KOKKOS_INLINE_FUNCTION
Real EmfSmoothW(const Real fd, const Real eps, const Real rhoe) {
  const Real del2 = eps*eps*rhoe;
  if (del2 <= 0.0) {   // eps = 0: GS07 zero-mass-flux rule everywhere
    return 0.5;
  }
  return 0.5*(1.0 + fd/sqrt(fd*fd + del2));
}
}  // namespace mhd
#endif  // MHD_EMF_SMOOTH_HPP_
