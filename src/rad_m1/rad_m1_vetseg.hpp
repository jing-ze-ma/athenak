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
#include <cstring>

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

// mixprec-1014 (study, build flag -DM1_VGD_FP32T; used only by the vet_gd shell kernel):
// copies of VetSegW / VetSegFace / VetSegSplit whose exponentials have a float core,
// VetExpF: range-reduced in Real (n = rint(x / ln 2), r = x - n ln 2), exp(r) in float,
// scaled by 2^n: relative error ~1e-7 at any argument.  The originals above are unchanged.
KOKKOS_INLINE_FUNCTION
Real VetExpF(const Real x) {
  if (!(x > -708.0)) {return (x == x) ? 0.0 : x;}
  if (x > 709.0) {return exp(x);}
  const Real nn = rint(x*1.4426950408889634);
  const float r = static_cast<float>(x - nn*0.6931471805599453);
  // 2^n built from its exponent bits (n in [-1021, 1023] by the range tests above)
  const long long eb = static_cast<long long>(static_cast<int>(nn) + 1023) << 52;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  const Real sc = __longlong_as_double(eb);
#else
  Real sc;
  std::memcpy(&sc, &eb, sizeof(sc));
#endif
  return static_cast<Real>(expf(r))*sc;
}
KOKKOS_INLINE_FUNCTION
void VetSegWF(const Real dtau, Real &ex, Real &w0, Real &wu) {
  ex = VetExpF(-dtau);
  if (dtau < 1.0e-3) {
    w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
    wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
  } else {
    const Real g = (1.0 - ex)/dtau;
    w0 = 1.0 - g;
    wu = g - ex;
  }
}
KOKKOS_INLINE_FUNCTION
Real VetSegFaceF(const Real iu, const Real su, const Real s0, const Real hu,
                 const Real h0, const int emit) {
  const Real ht = hu + h0;
  if (!(ht > 0.0)) {return su;}
  const Real gf = (s0 - su)/ht;
  Real g = gf;
  if (hu > h0 && (!emit || su > s0)) {
    const Real gb = su - iu;
    const Real mm = (gb*gf > 0.0) ? ((fabs(gb) < fabs(gf)) ? gb : gf) : 0.0;
    if (emit == 2) {
      const Real th = (hu - h0)/ht;
      g = gf + th*((1.0 - VetExpF(-hu))*mm - gf);
    } else {
      const Real th = (hu - h0)/ht*(1.0 - VetExpF(-hu));
      g = gf + th*(mm - gf);
    }
  }
  return su + hu*g;
}
KOKKOS_INLINE_FUNCTION
Real VetSegSplitF(const Real iu, const Real su, const Real s0, const Real hu,
                  const Real h0, const int emit) {
  const Real sf = VetSegFaceF(iu, su, s0, hu, h0, emit);
  Real exu, w0u, wuu, ex0, w00, wu0;
  VetSegWF(hu, exu, w0u, wuu);
  VetSegWF(h0, ex0, w00, wu0);
  const Real ia = iu*exu + wuu*su + w0u*sf;
  if (emit == 2 && h0 > hu && s0 > su) {
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


// open-od-1013: vet_path = cells.  The x1-plane sweep steps ONE x1 column per segment; a
// ray steeper than the cell diagonal (|sh2| = dx1|mu2|/(|mu1| dx2) > 1, or |sh3| > 1)
// jumps |sh| rows in that one segment and sees chi and S only at its two ends, so an
// opaque body narrower than the jump is invisible to it (the single opaque disc: the
// formal solution J near the x2 axis 0.86 of the exact staircase J at 15 cells per
// diameter, first order in dx/diameter; the x1-near rays are exact).  vet_path = cells
// subdivides the step at the NODES where the characteristic crosses a row (x2) or layer
// (x3) centre between the foot and the cell: each sub-segment is the code's own segment
// (VetSegOne: vet_seg split / jump, linear S in tau), the first one from the foot
// stencil as before (VetUpCorner / lin), the node chi and S linear in x1 between the
// upwind and own columns (and in the other transverse index), or, at a material edge
// (stencil chi contrast > kVetCornerRatio, as vet_interp auto), the cell that contains
// the node.  A linear S with constant chi gives the single segment again (linear-in-tau
// sub-segments of a linear source are exact): the diffusion limit is kept.  Rays with
// |sh2|, |sh3| <= 1 have no node: bitwise the single segment.  Only a step whose cells
// (foot stencil, node cells, own cell) have a chi contrast > kVetCornerRatio, or an S
// contrast > kVetCornerRatio among those not thin over the step (chi hl > 0.01), is
// subdivided (the material edge the single segment misses); a smooth step keeps the
// single segment bitwise.
// The accessor acc(col, dk, dj, chi, S) returns the extinction and source of the column
// col (0 upwind, 1 own) at the transverse offsets (dk, dj) from the cell, and false
// where that cell lies outside an OPEN (vacuum) face (the ray there carries I = 0).
template <class A>
KOKKOS_INLINE_FUNCTION
bool VetPathNode(const A &acc, const int vim, const Real s, const int dj0, const Real aj,
                 const int dk0, const Real ak, Real &cn, Real &sn) {
  const Real wc[2] = {1.0 - s, s}, wj[2] = {1.0 - aj, aj}, wk[2] = {1.0 - ak, ak};
  Real wv[8], cv[8], sv[8];
  int nv = 0;
  Real cmn = 1.0e300, cmx = 0.0;
  for (int c = 0; c < 2; ++c) {
    for (int b = 0; b < 2; ++b) {
      for (int a = 0; a < 2; ++a) {
        const Real w = wc[c]*wk[b]*wj[a];
        if (!(w > 0.0)) continue;
        acc(c, dk0 + b, dj0 + a, cv[nv], sv[nv]);
        wv[nv] = w;
        cmn = fmin(cmn, cv[nv]);
        cmx = fmax(cmx, cv[nv]);
        ++nv;
      }
    }
  }
  const int cc = (s < 0.5) ? 0 : 1;
  const int dj = dj0 + ((aj < 0.5) ? 0 : 1), dk = dk0 + ((ak < 0.5) ? 0 : 1);
  Real ce, se;
  const bool in = acc(cc, dk, dj, ce, se);
  if (cmx > kVetCornerRatio*cmn) {   // a material edge: the cell that holds the node
    cn = ce;
    sn = se;
    return in;
  }
  Real ws = 0.0, lc = 0.0, ls = 0.0;
  const bool lg = (vim == 1 || vim == 3);
  for (int q = 0; q < nv; ++q) {
    ws += wv[q];
    lc += wv[q]*(lg ? log(fmax(cv[q], 1.0e-300)) : cv[q]);
    ls += wv[q]*(lg ? log(fmax(sv[q], 1.0e-300)) : sv[q]);
  }
  cn = lg ? exp(lc/ws) : lc/ws;
  sn = lg ? exp(ls/ws) : ls/ws;
  return in;
}

// the subdivided step (see above): returns false (iv untouched) where the step has no
// node; iul, cul, sul: the lin (vet_interp = lin) upwind values; hl = dx1/|mu1|
template <class A>
KOKKOS_INLINE_FUNCTION
bool VetPathStep(const A &acc, const int vim, const Real w[4], const Real iq[4],
                 const Real cq[4], const Real sq[4], const Real iul, const Real cul,
                 const Real sul, const Real s0, const Real c0, const Real hl,
                 const bool spl, const int spe, const Real sh2, const Real sh3,
                 Real &iv) {
  const Real a2 = fabs(sh2), a3 = fabs(sh3);
  const int n2 = (a2 > 1.0) ? (static_cast<int>(ceil(a2)) - 1) : 0;
  const int n3 = (a3 > 1.0) ? (static_cast<int>(ceil(a3)) - 1) : 0;
  if (n2 + n3 == 0) return false;
  const int g2 = (sh2 < 0.0) ? -1 : 1, g3 = (sh3 < 0.0) ? -1 : 1;
  // the gate: subdivide only where the cells the step passes (the cell that holds each
  // node, the foot stencil, the cell itself) have a chi or S contrast > kVetCornerRatio;
  // a smooth path keeps the single segment (bitwise; the diffusion limit untouched)
  {
    // S counts only in cells not thin over the step (chi hl > 1e-2): a skipped emitter
    const Real tmin = 1.0e-2/hl;
    Real cmn = c0, cmx = c0, smn = 1.0e300, smx = 0.0;
    if (c0 > tmin) {smn = s0; smx = s0;}
    for (int q = 0; q < 4; ++q) {
      if (w[q] > 0.0) {
        cmn = fmin(cmn, cq[q]); cmx = fmax(cmx, cq[q]);
        if (cq[q] > tmin) {smn = fmin(smn, sq[q]); smx = fmax(smx, sq[q]);}
      }
    }
    int u2 = 0, u3 = 0;
    while (u2 < n2 || u3 < n3) {
      const Real s2 = (u2 < n2) ? (1.0 - (n2 - u2)/a2) : 2.0;
      const Real s3 = (u3 < n3) ? (1.0 - (n3 - u3)/a3) : 2.0;
      int dj, dk;
      Real sn;
      if (s2 <= s3) {
        sn = s2;
        dj = g2*(n2 - u2);
        dk = static_cast<int>(floor(sh3*(1.0 - sn) + 0.5));
        ++u2;
        if (s3 == s2) ++u3;
      } else {
        sn = s3;
        dk = g3*(n3 - u3);
        dj = static_cast<int>(floor(sh2*(1.0 - sn) + 0.5));
        ++u3;
      }
      Real cv, sv;
      acc((sn < 0.5) ? 0 : 1, dk, dj, cv, sv);
      cmn = fmin(cmn, cv); cmx = fmax(cmx, cv);
      if (cv > tmin) {smn = fmin(smn, sv); smx = fmax(smx, sv);}
    }
    if (!(cmx > kVetCornerRatio*cmn) && !(smx > kVetCornerRatio*smn)) return false;
  }
  int t2 = 0, t3 = 0;
  Real sp = 0.0, cp = 0.0, spv = 0.0;
  bool first = true;
  Real ivc = 0.0;
  while (true) {
    const Real s2 = (t2 < n2) ? (1.0 - (n2 - t2)/a2) : 2.0;
    const Real s3 = (t3 < n3) ? (1.0 - (n3 - t3)/a3) : 2.0;
    Real sn, cn, snv;
    bool in = true;
    if (s2 > 1.0 && s3 > 1.0) {
      sn = 1.0;
      cn = c0;
      snv = s0;
    } else if (s2 <= s3) {
      sn = s2;
      const Real z = sh3*(1.0 - sn);
      const Real fz = floor(z);
      in = VetPathNode(acc, vim, sn, g2*(n2 - t2), 0.0, static_cast<int>(fz), z - fz,
                       cn, snv);
      ++t2;
      if (s3 == s2) ++t3;   // a row-layer corner: one node
    } else {
      sn = s3;
      const Real y = sh2*(1.0 - sn);
      const Real fy = floor(y);
      in = VetPathNode(acc, vim, sn, static_cast<int>(fy), y - fy, g3*(n3 - t3), 0.0,
                       cn, snv);
      ++t3;
    }
    if (first) {
      ivc = (vim != 0) ? VetUpCorner(vim, w, iq, cq, sq, snv, cn, hl*sn, spl, spe)
                       : VetSegOne(iul, sul, snv, cul, cn, hl*sn, spl, spe);
      first = false;
    } else {
      ivc = VetSegOne(ivc, spv, snv, cp, cn, hl*(sn - sp), spl, spe);
    }
    if (!in) ivc = 0.0;   // the node lies beyond an open face: the ray entered there
    if (sn >= 1.0) break;
    sp = sn;
    cp = cn;
    spv = snv;
  }
  iv = ivc;
  return true;
}

#endif  // RAD_M1_RAD_M1_VETSEG_HPP_
