#ifndef RAD_M1_RAD_M1_VETSEG_HPP_
#define RAD_M1_RAD_M1_VETSEG_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetseg.hpp
//! \brief segfix-1011: <rad_m1>/vet_seg = split, the short-characteristics segment from
//! the upwind point u to the cell centre 0 integrated as TWO sub-segments split at the
//! cell face: the upwind part with chi_u (optical depth hu), the own part with chi_0
//! (h0).
//!
//! The source is piecewise linear in tau, continuous at the face value
//!   S_f = S_u + hu*g,   g = gf + th*(minmod(gb, gf) - gf),
//!   gf = (S_0 - S_u)/(hu + h0)      the segment's mean slope (dS/dtau along the ray),
//!   gb = S_u - I_u                  the upwind cell's own slope: I - S = -dS/dtau of the
//!                                   formal solution inside an opaque cell,
//!   th = max(hu - h0, 0)/(hu + h0)*(1 - exp(-hu))   the opacity contrast, upwind thick.
//! With th = 0 (any segment with hu <= h0, in particular constant chi) S_f is the face
//! value of the single linear-in-tau S of vet_seg = mean and the two sub-segments are
//! that segment exactly (to round-off).  In a smooth field th = O(dx), so the change is
//! second order.  At an opaque -> thin interface (hu >> 1 >> h0) with a uniform hot
//! upwind body (I_u = S_u, gb = 0) or a single hot cell (gb opposite to gf), S_f = S_u
//! and the emergent intensity is ~S_u (vet_seg = mean gives ~S_0 + S_u/hu).  An
//! unresolved diffusion photosphere (S linear in tau inside the opaque cell: gb = gf)
//! keeps the mean result.  S_f lies between S_u and S_0 and all weights are >= 0:
//! positive, monotone.
//! vet_seg = emit (emit = true): the limiter acts only where S_u > S_0, the opaque cell
//! EMITTING into the thinner one; an opaque cell heated from its thin side (S_u <= S_0:
//! the unresolved skin is near the incident field) keeps the mean face value.
//! vet_seg = jump (coreL-1012, opt-in; emit = 2): the face value of a material
//! discontinuity is set by the chi CONTRAST, not by the upwind cell's opacity: with
//! th = (hu - h0)/(hu + h0) (S_u > S_0 as emit) S_f = S_u + hu*(gf + th*((1 - e^-hu) mm
//! - gf)): a thin edge cell of an emitting body keeps its own S up to the face (emit
//! fell back to the mean ramp S_u -> S_0 there, losing half of the outer half-cell's
//! emission: the T6b core luminosity of the formal solution -6 % at tau_cell 0.13 and
//! not first order).  And where the OWN cell is the thicker and S_0 > S_u (a ray
//! entering an emitting body), the own part of the segment starts from
//! S_f0 = S_f + (h0 - hu)/(h0 + hu) (S_0 - S_f) instead of the continuous S_f (S jumps
//! at the face): the inward intensity of the edge cells (H +5 %, J -9 % at tau_cell
//! 0.13 with emit).  In a smooth field both contrasts are O(dx): second order kept.

#include <cmath>

#include "athena.hpp"

// the positive first-order SC weights of one linear-in-tau sub-segment
KOKKOS_INLINE_FUNCTION
void VetSegW(const Real dtau, Real &ex, Real &w0, Real &wu) {
  ex = exp(-dtau);
  if (dtau < 1.0e-3) {
    w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
    wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
  } else {
    const Real g = (1.0 - ex)/dtau;
    w0 = 1.0 - g;
    wu = g - ex;
  }
}

// the face source S_f of the split segment (see above); emit: 0 split, 1 emit, 2 jump
KOKKOS_INLINE_FUNCTION
Real VetSegFace(const Real iu, const Real su, const Real s0, const Real hu,
                const Real h0, const int emit) {
  const Real ht = hu + h0;
  if (!(ht > 0.0)) {return su;}
  const Real gf = (s0 - su)/ht;
  Real g = gf;
  if (hu > h0 && (!emit || su > s0)) {
    const Real gb = su - iu;
    const Real mm = (gb*gf > 0.0) ? ((fabs(gb) < fabs(gf)) ? gb : gf) : 0.0;
    if (emit == 2) {
      // coreL-1012 jump: the contrast alone moves S_f to the upwind cell's own value;
      // only the slope estimate gb (valid in an opaque cell) keeps the factor 1 - e^-hu
      const Real th = (hu - h0)/ht;
      g = gf + th*((1.0 - exp(-hu))*mm - gf);
    } else {
      const Real th = (hu - h0)/ht*(1.0 - exp(-hu));
      g = gf + th*(mm - gf);
    }
  }
  return su + hu*g;
}

// the intensity at the cell centre from I_u at the upwind point, the split segment
KOKKOS_INLINE_FUNCTION
Real VetSegSplit(const Real iu, const Real su, const Real s0, const Real hu,
                 const Real h0, const int emit) {
  const Real sf = VetSegFace(iu, su, s0, hu, h0, emit);
  Real exu, w0u, wuu, ex0, w00, wu0;
  VetSegW(hu, exu, w0u, wuu);
  VetSegW(h0, ex0, w00, wu0);
  const Real ia = iu*exu + wuu*su + w0u*sf;
  if (emit == 2 && h0 > hu && s0 > su) {
    // coreL-1012 jump: the OWN cell the thicker and the hotter (a ray entering an
    // emitting body): its side of the face keeps S_0 by the contrast (S may jump at
    // the face; the upwind part keeps S_u -> S_f)
    const Real sf0 = sf + (h0 - hu)/(h0 + hu)*(s0 - sf);
    return ia*ex0 + wu0*sf0 + w00*s0;
  }
  return ia*ex0 + wu0*sf + w00*s0;
}

// unify-sc-gd-1012 (#6): one SC segment from the upwind point (I_u, chi_u, S_u) to the
// cell centre (chi_0, S_0) over the path length hl = dx1/|mu1|, vet_seg split or mean
KOKKOS_INLINE_FUNCTION
Real VetSegOne(const Real iu, const Real su, const Real s0, const Real cu, const Real c0,
               const Real hl, const bool spl, const int spe) {
  if (spl) {return VetSegSplit(iu, su, s0, 0.5*cu*hl, 0.5*c0*hl, spe);}
  Real ex, w0, wu;
  VetSegW(0.5*(cu + c0)*hl, ex, w0, wu);
  return iu*ex + wu*su + w0*s0;
}

// unify-sc-gd-1012 (#6) vet_interp: the upwind segment from the (bi)linear stencil of
// the upwind plane, weights w[4] (sum 1), corner intensities, extinctions, sources.
//   1 log:    I_u = sum w I_q, chi_u and S_u geometric (sum w ln), one segment (vet_gd's
//             interpolation);
//   2 corner: one segment per corner, I = sum w I(segment from corner q): the OUTCOME is
//             interpolated.  At a sharp opaque edge the plain interpolation (0 lin: chi_u,
//             S_u linear, one segment) puts the opaque chi_u on a ray that passes the
//             edge and replaces its intensity by the interpolated S_u (a ray grazing an
//             opaque disc of S = B loses the limb: L_fs 0.58-0.66 of exact, resolution
//             independent, scemu.py); corner keeps the thin corner's intensity with its
//             weight (0.95-1.0).  Smooth fields: corner - lin = O(dx^2).
//   3 auto:   corner where the stencil's extinction contrast max/min chi exceeds
//             kVetCornerRatio (a material edge), log elsewhere (a smooth field: the two
//             agree to O(dx^2) there, and log costs one segment instead of up to four)
constexpr Real kVetCornerRatio = 4.0;

KOKKOS_INLINE_FUNCTION
Real VetUpCorner(int mode, const Real w[4], const Real iq[4], const Real cq[4],
                 const Real sq[4], const Real s0, const Real c0, const Real hl,
                 const bool spl, const int spe) {
  if (mode == 3) {
    Real cmn = 1.0e300, cmx = 0.0;
    for (int q = 0; q < 4; ++q) {
      if (w[q] > 0.0) {cmn = fmin(cmn, cq[q]); cmx = fmax(cmx, cq[q]);}
    }
    mode = (cmx > kVetCornerRatio*cmn) ? 2 : 1;
  }
  if (mode == 2) {
    Real iv = 0.0;
    for (int q = 0; q < 4; ++q) {
      if (w[q] > 0.0) {iv += w[q]*VetSegOne(iq[q], sq[q], s0, cq[q], c0, hl, spl, spe);}
    }
    return iv;
  }
  Real iu = 0.0, lc = 0.0, ls = 0.0;
  for (int q = 0; q < 4; ++q) {
    iu += w[q]*iq[q];
    lc += w[q]*log(fmax(cq[q], 1.0e-300));
    ls += w[q]*log(fmax(sq[q], 1.0e-300));
  }
  return VetSegOne(iu, exp(ls), s0, exp(lc), c0, hl, spl, spe);
}

#endif  // RAD_M1_RAD_M1_VETSEG_HPP_
