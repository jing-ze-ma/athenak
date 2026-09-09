#ifndef UTILS_ATM_COLUMN_HPP_
#define UTILS_ATM_COLUMN_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file atm_column.hpp
//! \brief Small pure helpers for a spherical atmosphere on a radial column: the cell's
//! direction on the cubed sphere, the gravity model (constant or point mass, with the
//! host star's tide), and the cached-temperature guess.
//!
//! Extracted VERBATIM from src/pgen/deep_hot_jupiter_rt.cpp so that the two-stream
//! radiation module and more than one problem generator can share them.  Every one takes
//! all of its inputs as arguments and holds no state.

#include <math.h>

#include "athena.hpp"
#include "coordinates/cubed_sphere.hpp"

namespace atm_column {

//----------------------------------------------------------------------------------------
//! \fn Real TGuess
//! \brief the cached temperature of a cell, or "no guess" when there is no cache.
//!
//! Hydro/MHD::wtemp is allocated ONLY under a general EOS -- an ideal gas has nothing to
//! cache, since T is algebraic. Under `eos = ideal` the view is therefore empty, and
//! indexing it is a read through a null pointer: silent in a Release build until it
//! happens to land outside the mapped heap, which is what made the shipped ideal-gas
//! input segfault inside the RT a few cycles in. The consumers all ignore the guess on
//! the ideal branch anyway, so returning a non-positive "no guess" here is exact.
KOKKOS_INLINE_FUNCTION
Real TGuess(const DvceArray4D<Real> &wt, const int m, const int k, const int j,
            const int i) {
  return (wt.extent(0) > 0) ? wt(m,k,j,i) : -1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void CSCellAngles
//! \brief The cell's spherical direction on the CUBED SPHERE, in the same convention the
//! spherical-polar branch uses: theta the COLATITUDE, lam the LATITUDE, phi the LONGITUDE
//! measured from the SUBSTELLAR point.
//!
//! On the cubed sphere x2 and x3 are the two panel-tangential coordinates on [-1,1], not
//! angles, so there is no expression in them alone -- the direction depends on WHICH
//! PANEL the cell is on.  Go through the chart: xi = pi/4 * x2, eta = pi/4 * x3, then
//! PanelToCart gives the unit direction and the angles follow from it.
//!
//! THE LONGITUDE ORIGIN MUST MATCH.  The spherical-polar branch uses phi = x3v - M_PI, so
//! the substellar point sits at coordinate longitude pi, i.e. along -x.  atan2(-cy,-cx)
//! is atan2(cy,cx) rotated by pi and already lands in (-pi, pi], which is the same
//! interval that branch produces -- no wrapping needed.
KOKKOS_INLINE_FUNCTION
void CSCellAngles(const int panel, const Real x2v, const Real x3v,
                  Real &theta, Real &lam, Real &phi) {
  Real q[3];
  cubed_sphere::PanelToCart(panel, 0.25*M_PI*x2v, 0.25*M_PI*x3v, q);
  const Real cx = q[0], cy = q[1], cz = q[2];
  const Real cr = sqrt(cx*cx + cy*cy + cz*cz);
  const Real cth = (cr > 0.0) ? (cz/cr) : 0.0;
  theta = acos(fmin(1.0, fmax(-1.0, cth)));
  lam = 0.5*M_PI - theta;
  phi = atan2(-cy, -cx);
}

KOKKOS_INLINE_FUNCTION
Real GravAccAt(const Real grav_acc, const Real ap, const Real r, const bool pmass) {
  return pmass ? grav_acc*SQR(ap/r) : grav_acc;
}

KOKKOS_INLINE_FUNCTION
Real GravPotAt(const Real grav_acc, const Real ap, const Real r, const Real z,
               const bool pmass) {
  return pmass ? (-grav_acc)*ap*(1.0 - ap/r) : (-grav_acc)*z;
}

KOKKOS_INLINE_FUNCTION
Real TideAccR(const Real omega, const Real r, const Real mu, const bool tide) {
  return tide ? SQR(omega)*r*(3.0*SQR(mu) - 1.0) : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real TideAccT(const Real omega, const Real r, const Real sinth, const Real costh,
              const Real cosph, const bool tide) {
  return tide ? 3.0*SQR(omega)*r*sinth*costh*SQR(cosph) : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real TideAccP(const Real omega, const Real r, const Real sinth, const Real sinph,
              const Real cosph, const bool tide) {
  return tide ? -3.0*SQR(omega)*r*sinth*sinph*cosph : 0.0;
}

KOKKOS_INLINE_FUNCTION
Real EffGravAt(const Real grav, const Real ap, const Real r, const bool pmass,
               const Real omega, const Real mu, const bool tide) {
  const Real g = GravAccAt(grav, ap, r, pmass);
  return tide ? fmax(g - TideAccR(omega, r, mu, tide), 0.1*g) : g;
}

}  // namespace atm_column

#endif  // UTILS_ATM_COLUMN_HPP_
