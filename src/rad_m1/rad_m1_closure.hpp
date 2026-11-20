#ifndef RAD_M1_RAD_M1_CLOSURE_HPP_
#define RAD_M1_RAD_M1_CLOSURE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_closure.hpp
//! \brief Levermore (1984) M1 closure, the flux limiter and the Skinner & Ostriker
//! (2013, eq. 41a) closed-form wave speeds, as device-inline functions.
//!
//! Conventions (docs/dev/rad_m1_design.md sect. 1 and 3): E is the lab-frame radiation
//! energy density, F_i the lab-frame flux, f_i = F_i/(c E) the reduced flux with
//! |f| <= 1, n_i = f_i/|f| and
//!
//!   chi  = (3 + 4 f^2)/(5 + 2 sqrt(4 - 3 f^2))
//!   P_ij = E [ (1-chi)/2 delta_ij + (3 chi - 1)/2 n_i n_j ]
//!
//! All of this is algebraic in flat space: there is no root find anywhere in here.

#include <math.h>
#include "athena.hpp"

namespace radm1 {

// <rad_m1>/thick_flux: which asymptotic-preserving correction the face flux carries
// (design sect. 3).  Device code branches on these, so they are plain ints.
constexpr int M1_THICK_NONE   = 0;   // plain HLL
constexpr int M1_THICK_APHLL  = 1;   // Bloch et al. (2021) alpha on the E-flux only
constexpr int M1_THICK_SCALED = 2;   // Jiang (2021) App. A: both wave speeds x eps(tau)

// <rad_m1>/ap_form: which algebraic form of the thick_flux = ap_hll E-flux is used
// (design sect. 10 and the 1c measurement recorded there).
constexpr int M1_APFORM_UNIFIED = 0;  // alpha F_HLL + (1-alpha) F_diff(1 - dEr/dEc)
constexpr int M1_APFORM_ALPHA2  = 1;  // the 1b pair: Berthon for dc, alpha^2 for plm

//----------------------------------------------------------------------------------------
//! \fn M1Chi
//! \brief closure factor chi as a function of the reduced flux magnitude; chi(0) = 1/3
//! (isotropic), chi(1) = 1 (free streaming) for every kind.
//!   M1_CHI_LEVERMORE  Levermore (1984), the Lorentz-boosted isotropic field (default)
//!   M1_CHI_MINERBO    Minerbo (1978) classical maximum entropy, the usual polynomial fit
//!                     chi = 1/3 + (2 f^2/15)(3 - f + 3 f^2)
//!   M1_CHI_KERSHAW    Kershaw (1976), chi = (1 + 2 f^2)/3
//! The explicit HLL wave speeds are Levermore's, so the other two are accepted with the
//! implicit transports only (docs/dev/rad_m1_closure_survey.md).

constexpr int M1_CHI_LEVERMORE = 0;
constexpr int M1_CHI_MINERBO   = 1;
constexpr int M1_CHI_KERSHAW   = 2;

KOKKOS_INLINE_FUNCTION
Real M1Chi(const Real fnorm, const int kind = M1_CHI_LEVERMORE) {
  Real f2 = fnorm*fnorm;
  if (f2 > 1.0) {f2 = 1.0;}
  if (kind == M1_CHI_MINERBO) {
    Real f = sqrt(f2);
    return 1.0/3.0 + (2.0*f2/15.0)*(3.0 - f + 3.0*f2);
  } else if (kind == M1_CHI_KERSHAW) {
    return (1.0 + 2.0*f2)/3.0;
  }
  Real s = sqrt(fmax(4.0 - 3.0*f2, 0.0));
  return (3.0 + 4.0*f2)/(5.0 + 2.0*s);
}

//----------------------------------------------------------------------------------------
//! \fn M1ApplyLimits
//! \brief Apply the energy floor and the flux limit |F| <= c E to one cell state.  The
//! flux is scaled by (c E/|F|), NOT by the square of that ratio.  Returns true if either
//! limit was active.

KOKKOS_INLINE_FUNCTION
bool M1ApplyLimits(const Real cl, const Real e_floor,
                   Real &e, Real &f1, Real &f2, Real &f3) {
  bool limited = false;
  if (!(e > e_floor)) {
    e = e_floor;
    limited = true;
  }
  Real fmag = sqrt(f1*f1 + f2*f2 + f3*f3);
  Real fmax = cl*e;
  if (fmag > fmax) {
    Real scale = (fmag > 0.0) ? (fmax/fmag) : 0.0;
    f1 *= scale;
    f2 *= scale;
    f3 *= scale;
    limited = true;
  }
  return limited;
}

//----------------------------------------------------------------------------------------
//! \fn M1FluxNormCs
//! \brief STAGE CS2 (cubed sphere): the physical |F| of a flux whose (f2, f3) are
//! COVARIANT on the panel basis at the cell angle c = cos_cell:
//! |F|^2 = f1^2 + g^ab f_a f_b, g^ab = [[1, -c], [-c, 1]]/(1 - c^2).  The one place this
//! metric norm is written; the write-back clip (M1ApplyLimitsCs) and
//! implicit_realisable_coupling both use it.

KOKKOS_INLINE_FUNCTION
Real M1FluxNormCs(const Real c, const Real f1, const Real f2, const Real f3) {
  return sqrt(f1*f1 + (f2*f2 + f3*f3 - 2.0*c*f2*f3)/(1.0 - c*c));
}

//----------------------------------------------------------------------------------------
//! \fn M1ApplyLimitsCs
//! \brief STAGE CS2 (cubed sphere): M1ApplyLimits for a flux whose (f2, f3) are COVARIANT
//! on the panel basis at the cell angle c = cos_cell, with |F| = M1FluxNormCs.

KOKKOS_INLINE_FUNCTION
bool M1ApplyLimitsCs(const Real cl, const Real e_floor, const Real c,
                     Real &e, Real &f1, Real &f2, Real &f3) {
  bool limited = false;
  if (!(e > e_floor)) {
    e = e_floor;
    limited = true;
  }
  Real fmag = M1FluxNormCs(c, f1, f2, f3);
  Real fmax = cl*e;
  if (fmag > fmax) {
    Real scale = (fmag > 0.0) ? (fmax/fmag) : 0.0;
    f1 *= scale;
    f2 *= scale;
    f3 *= scale;
    limited = true;
  }
  return limited;
}

//----------------------------------------------------------------------------------------
//! \fn M1FillGhost
//! \brief ONE ghost cell of the physical x1/x2/x3 boundaries, from the CONSERVED moments
//! of the active cell (ka,ja,ia) it continues.  Moved here verbatim from the anonymous
//! namespace of bvals/physics/rad_m1_bcs.cpp so that a problem generator whose mesh flag
//! is `user` -- and whose own BC routine therefore replaces RadM1BCs on that face -- can
//! reuse the SAME fill instead of re-deriving it (box_convection's M1 top wall).
//!
//! mode 0 = plain copy (outflow), 1 = reflect (flip the normal flux), 2 = vacuum, i.e. a
//! DARK ghost: no incoming radiation, ever.  inorm = 1,2,3 selects the face-normal flux
//! component; sgn = +1 when "outgoing" is the positive coordinate direction.

KOKKOS_INLINE_FUNCTION
void M1FillGhost(const DvceArray5D<Real> &u, const int m,
                 const int kg, const int jg, const int ig,
                 const int ka, const int ja, const int ia,
                 const int inorm, const int mode, const Real sgn,
                 const Real cl, const Real efl) {
  Real e  = u(m,M1_E, ka,ja,ia);
  Real f1 = u(m,M1_F1,ka,ja,ia);
  Real f2 = u(m,M1_F2,ka,ja,ia);
  Real f3 = u(m,M1_F3,ka,ja,ia);
  if (mode == 1) {
    if (inorm == 1) {
      f1 = -f1;
    } else if (inorm == 2) {
      f2 = -f2;
    } else {
      f3 = -f3;
    }
  } else if (mode == 2) {
    // VACUUM: the ghost is ALWAYS dark, E = e_floor and F = 0, i.e. no incoming
    // radiation whatever the interior state does.  A beam at f = 1 still leaves
    // unimpeded, because the HLL left wave speed of a dark state is b_L = 0.
    //
    // The earlier "copy the interior when the interior flux points out" fill was
    // wrong for two reasons.  Zeroing only the normal component (the version before
    // that) leaves the ghost a reservoir of E with no flux and the HLL flux then PUMPS
    // energy in (38 % growth of a grazing beam).  Copying the WHOLE state is the same
    // defect in a milder form: it represents incoming radiation equal to the interior's
    // own inward half, so dE/dn = 0 at the face and NO relation between F and E is
    // imposed.  At a plane-parallel free surface the physical condition is f = 1/2, not
    // 1, and with the copy the surface value of E floats -- measured on the He column:
    // E at the top grew 137x over 1000 s and the emergent flux fell to 0.36 of the
    // imposed one.  A dark ghost supplies the Marshak-like relation through the HLL
    // flux and the same column stays flat to 2.4 %.
    (void) sgn;
    e = efl;
    f1 = 0.0;
    f2 = 0.0;
    f3 = 0.0;
  }
  (void) cl;
  u(m,M1_E, kg,jg,ig) = e;
  u(m,M1_F1,kg,jg,ig) = f1;
  u(m,M1_F2,kg,jg,ig) = f2;
  u(m,M1_F3,kg,jg,ig) = f3;
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_CLOSURE_HPP_
