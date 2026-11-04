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

//----------------------------------------------------------------------------------------
//! \fn Real EmfRampW(fd, cdx, rhosum)
//! \brief ctramp-1103, <mhd>/emf_ramp = athenapp (read only when named; default off,
//! bitwise): the Athena++ corner-EMF weight (Hydro::GetWeightForCT, used by
//! Field::ComputeCornerE): v = clamp(1024 dt F_rho / (dx (rho_L + rho_R)), -1/2, 1/2),
//! weights (1/2 + v) on the LEFT (lower-index, F_rho > 0 upwind) choice and (1/2 - v) on
//! the right.  cdx = 1024*dt/dx of the face direction; Athena++ uses the FULL step
//! dt = Mesh::dt in every RK stage (not beta*dt) and dx = CenterWidth of the face's
//! right cell; rho_L, rho_R are there the reconstructed face states, here the two
//! adjacent cell-centred densities (no face states are kept after the flux sweep).
//! A linear ramp of the face Courant number: upwind (GS07) once |v_n| dt/dx > 1/512,
//! the plain average as dt -> 0 -- so the semi-discrete limit is the AVERAGE scheme.
KOKKOS_INLINE_FUNCTION
Real EmfRampW(const Real fd, const Real cdx, const Real rhosum) {
  const Real v = cdx*fd/rhosum;
  return 0.5 + fmax(-0.5, fmin(0.5, v));
}

//! \brief the weight of the LEFT choice for MHD::CornerE: emf_ramp when cdx > 0 (then
//! ra, rb are the cell densities), else emf_smooth (ra, rb = rho*e, see EmfSmoothW)
KOKKOS_INLINE_FUNCTION
Real EmfWeightL(const Real cdx, const Real fd, const Real eps, const Real ra,
                const Real rb) {
  return (cdx > 0.0) ? EmfRampW(fd, cdx, ra + rb) : EmfSmoothW(fd, eps, 0.5*(ra + rb));
}
}  // namespace mhd
#endif  // MHD_EMF_SMOOTH_HPP_
