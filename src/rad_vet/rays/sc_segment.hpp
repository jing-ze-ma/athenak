#ifndef RAD_VET_RAYS_SC_SEGMENT_HPP_
#define RAD_VET_RAYS_SC_SEGMENT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_segment.hpp
//! \brief Short-characteristics segment kernels (vet_seg, vet_interp), vet_gd launch
//! policies and the PLM foot of the transverse interpolation.
//!
//! vet_seg = split: the segment u -> 0 is two sub-segments split at the cell face
//! (upwind optical depth hu, own h0), S piecewise linear in tau, continuous at
//!   S_f = S_u + hu*g,  g = gf + th*(minmod(gb, gf) - gf),
//!   gf = (S_0 - S_u)/(hu + h0),  gb = S_u - I_u (slope inside an opaque cell),
//!   th = max(hu - h0, 0)/(hu + h0)*(1 - exp(-hu)).
//! th = 0 for hu <= h0 gives vet_seg = mean; th = O(dx) in smooth fields (2nd order).
//! S_f lies between S_u and S_0, all weights >= 0 (positive, monotone).
//! emit (1): limiter only where S_u > S_0 (the opaque cell emits into the thinner one).
//! jump (2): th = (hu - h0)/(hu + h0)*q with plateau factor
//!   q = 1 - (S_uu - S_u)/(S_u - S_0) in [0, 1] (S_uu one cell further upwind, q = 1
//!   without one), so an isothermal body keeps S_u to the face while a diffusion skin
//!   keeps the mean ramp; where the own cell is thicker and hotter the own part starts
//!   from S_f0 = S_f + (h0 - hu)/(h0 + hu)*q0*(S_0 - S_f) (mirror plateau factor q0).

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

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

// range policy passing the functor as a kernel argument (HintLightWeight) for the
// per-shell vet_gd kernels; avoids the constant-memory launch and its host wait
inline auto ScLight(const DevExeSpace &ex, const int a, const int b) {
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
}

// vet_gd shell kernel policy: ScLight with launch bounds on HIP, so the register
// allocator may use 512/VET_VGD_LB VGPRs instead of spilling (arithmetic unchanged)
#ifndef VET_VGD_LB
#define VET_VGD_LB 1
#endif
inline auto ScShellPol(const DevExeSpace &ex, const int a, const int b) {
#if defined(KOKKOS_ENABLE_HIP) && (VET_VGD_LB > 0)
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<Kokkos::LaunchBounds<256, VET_VGD_LB>>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
#else
  return ScLight(ex, a, b);
#endif
}

// CUDA gx-false shell kernel launch bounds (VET_VGD_CUDA_LB; 0 = plain ScShellPol,
// 2 = <128, 5> on A100)
#if defined(KOKKOS_ENABLE_CUDA)
#ifndef VET_VGD_CUDA_LB
#if defined(KOKKOS_ARCH_AMPERE80)
#define VET_VGD_CUDA_LB 2
#else
#define VET_VGD_CUDA_LB 0
#endif
#endif
template <unsigned T, unsigned B>
inline auto ScLbPol(const DevExeSpace &ex, const int a, const int b) {
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<Kokkos::LaunchBounds<T, B>>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
}
#endif

// copy a host kernel closure byte-wise into the device buffer buf (grown as needed) and
// return its device address: the shell kernels call it through this pointer instead of
// a ~2.7 kB kernel argument (copied to every thread's stack); deep_copy fences
template <class F>
inline const F *ScClosureToDevice(const F &f, DvceArray1D<char> &buf) {
  const size_t n = sizeof(F);
  if (buf.extent(0) < n) {buf = DvceArray1D<char>("vet_vgd_body", n);}
  using HostBytes = Kokkos::View<const char *, Kokkos::HostSpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  HostBytes h(reinterpret_cast<const char *>(&f), n);
  Kokkos::deep_copy(Kokkos::subview(buf, std::make_pair(static_cast<size_t>(0), n)), h);
  return reinterpret_cast<const F *>(buf.data());
}

// fatal-error exit of the vet_gd setup: file, tag, msg
[[noreturn]] inline void ScFatalAt(const char *file, const char *tag,
                                   const std::string &msg) {
  std::cout << "### FATAL ERROR in " << file << std::endl << tag << msg << std::endl;
  // _Exit: skip static destructors behind live Kokkos views
  std::_Exit(EXIT_FAILURE);
}

// face source S_f of the split segment; emit: 0 split, 1 emit, 2 jump, 3 jump with
// q = 1 (no S_uu given)
KOKKOS_INLINE_FUNCTION
Real VetSegFace(const Real iu, const Real su, const Real s0, const Real hu,
                const Real h0, const int emit, const Real suu = -1.0) {
  const Real ht = hu + h0;
  if (!(ht > 0.0)) {return su;}
  const Real gf = (s0 - su)/ht;
  Real g = gf;
  if (hu > h0 && (!emit || su > s0)) {
    const Real gb = su - iu;
    const Real mm = (gb*gf > 0.0) ? ((fabs(gb) < fabs(gf)) ? gb : gf) : 0.0;
    if (emit >= 2) {
      // jump: S_f moves toward S_u as far as the upwind cell is a plateau (q); only
      // the opaque-cell slope gb keeps the factor 1 - e^-hu
      const Real q = (suu > su) ? fmax(1.0 - (suu - su)/(su - s0), 0.0) : 1.0;
      const Real th = (hu - h0)/ht*q;
      g = gf + th*((1.0 - exp(-hu))*mm - gf);
    } else {
      const Real th = (hu - h0)/ht*(1.0 - exp(-hu));
      g = gf + th*(mm - gf);
    }
  }
  return su + hu*g;
}

// intensity at the face of the split segment (its upwind part)
KOKKOS_INLINE_FUNCTION
Real VetSegFaceI(const Real iu, const Real su, const Real s0, const Real hu,
                 const Real h0, const int emit) {
  const Real sf = VetSegFace(iu, su, s0, hu, h0, emit);
  Real ex, w0, wu;
  VetSegW(hu, ex, w0, wu);
  return iu*ex + wu*su + w0*sf;
}
// intensity at the cell centre from I_u, split segment; ifc (if given) gets the face I
KOKKOS_INLINE_FUNCTION
Real VetSegSplit(const Real iu, const Real su, const Real s0, const Real hu,
                 const Real h0, const int emit, Real *ifc = nullptr,
                 const Real suu = -1.0, const Real sdd = -1.0) {
  const Real sf = VetSegFace(iu, su, s0, hu, h0, emit, suu);
  Real exu, w0u, wuu, ex0, w00, wu0;
  VetSegW(hu, exu, w0u, wuu);
  VetSegW(h0, ex0, w00, wu0);
  const Real ia = iu*exu + wuu*su + w0u*sf;
  if (ifc != nullptr) *ifc = ia;
  if (emit >= 2 && h0 > hu && s0 > su) {
    // jump, own cell thicker and hotter: its side of the face keeps S_0 as far as it
    // is a plateau (q0, S_dd one cell downstream; mirror of VetSegFace's q)
    const Real q0 = (sdd > s0) ? fmax(1.0 - (sdd - s0)/(s0 - su), 0.0) : 1.0;
    const Real sf0 = sf + (h0 - hu)/(h0 + hu)*q0*(s0 - sf);
    return ia*ex0 + wu0*sf0 + w00*s0;
  }
  return ia*ex0 + wu0*sf + w00*s0;
}

// one SC segment from the upwind point (I_u, chi_u, S_u) to the cell centre (chi_0, S_0)
// over hl = dx1/|mu1|, split or mean; ifc (if given) gets the face I; suu, sdd as above
KOKKOS_INLINE_FUNCTION
Real VetSegOne(const Real iu, const Real su, const Real s0, const Real cu, const Real c0,
               const Real hl, const bool spl, const int spe, Real *ifc = nullptr,
               const Real suu = -1.0, const Real sdd = -1.0) {
  if (spl) {return VetSegSplit(iu, su, s0, 0.5*cu*hl, 0.5*c0*hl, spe, ifc, suu, sdd);}
  if (ifc != nullptr) *ifc = VetSegFaceI(iu, su, s0, 0.5*cu*hl, 0.5*c0*hl, 0);
  Real ex, w0, wu;
  VetSegW(0.5*(cu + c0)*hl, ex, w0, wu);
  return iu*ex + wu*su + w0*s0;
}

// vet_interp: upwind segment from the (bi)linear stencil of the upwind plane, weights
// w[4] (sum 1) and corner I, chi, S:
//   1 log:    I_u = sum w I_q, chi_u and S_u geometric means, one segment;
//   2 corner: one segment per corner, I = sum w I_q(segment): keeps a thin corner's
//             intensity at a sharp opaque edge (lin would put the opaque chi_u on a
//             grazing ray); smooth fields: corner - lin = O(dx^2);
//   3 auto:   corner where max/min chi > kVetCornerRatio (material edge), else log.
// ifc (if given) gets the face I; suq the corners' S one x1 cell further upwind (jump),
// sdd the own cell's S one x1 cell downstream.
constexpr Real kVetCornerRatio = 4.0;

KOKKOS_INLINE_FUNCTION
Real VetUpCorner(int mode, const Real w[4], const Real iq[4], const Real cq[4],
                 const Real sq[4], const Real s0, const Real c0, const Real hl,
                 const bool spl, const int spe, Real *ifc = nullptr,
                 const Real *suq = nullptr, const Real sdd = -1.0) {
  if (mode == 3) {
    Real cmn = 1.0e300, cmx = 0.0;
    for (int q = 0; q < 4; ++q) {
      if (w[q] > 0.0) {cmn = fmin(cmn, cq[q]); cmx = fmax(cmx, cq[q]);}
    }
    mode = (cmx > kVetCornerRatio*cmn) ? 2 : 1;
  }
  if (mode == 2) {
    Real iv = 0.0, fq = 0.0;
    if (ifc != nullptr) *ifc = 0.0;
    for (int q = 0; q < 4; ++q) {
      if (!(w[q] > 0.0)) continue;
      iv += w[q]*VetSegOne(iq[q], sq[q], s0, cq[q], c0, hl, spl, spe,
                           (ifc != nullptr) ? &fq : nullptr,
                           (suq != nullptr) ? suq[q] : -1.0, sdd);
      if (ifc != nullptr) *ifc += w[q]*fq;
    }
    return iv;
  }
  Real iu = 0.0, lc = 0.0, ls = 0.0, lu = 0.0;
  for (int q = 0; q < 4; ++q) {
    iu += w[q]*iq[q];
    lc += w[q]*log(fmax(cq[q], 1.0e-300));
    ls += w[q]*log(fmax(sq[q], 1.0e-300));
    if (suq != nullptr) lu += w[q]*log(fmax(suq[q], 1.0e-300));
  }
  return VetSegOne(iu, exp(ls), s0, exp(lc), c0, hl, spl, spe, ifc,
                   (suq != nullptr) ? exp(lu) : -1.0, sdd);
}

// second-order foot of the transverse interpolation: PLM remap of the upwind plane with
// van Leer slopes.  Adds a(1 - a)/2 (s_a - s_b) along x2 and x3 to the donor-cell
// (linear) value; the part a(1 - a)/2 s_a crosses the cell's low face (fl[0] x2, fl[1]
// x3), keeping the update in flux form.  ip(dk, dj) = I at (ka + dk, ja + dj).
template <class IP>
KOKKOS_INLINE_FUNCTION
Real VetPlmFoot(const IP &ip, const Real a2, const Real a3, const bool thrd, Real fl[2]) {
  auto sl = [](const Real m, const Real c, const Real p) {
    const Real dl = c - m, dr = p - c;
    return (dl*dr > 0.0) ? 2.0*dl*dr/(dl + dr) : 0.0;
  };
  const Real k2 = 0.5*a2*(1.0 - a2), k3 = 0.5*a3*(1.0 - a3);
  Real dv = 0.0;
  fl[0] = 0.0;
  fl[1] = 0.0;
  for (int b = 0; b < (thrd ? 2 : 1); ++b) {
    const Real wr = (b == 0) ? 1.0 - a3 : a3, wc = (b == 0) ? 1.0 - a2 : a2;
    const Real sa = sl(ip(b,-1), ip(b,0), ip(b,1));
    dv += wr*k2*(sa - sl(ip(b,0), ip(b,1), ip(b,2)));
    fl[0] += wr*k2*sa;
    if (!thrd) continue;
    const Real ta = sl(ip(-1,b), ip(0,b), ip(1,b));
    dv += wc*k3*(ta - sl(ip(0,b), ip(1,b), ip(2,b)));
    fl[1] += wc*k3*ta;
  }
  return dv;
}

#endif  // RAD_VET_RAYS_SC_SEGMENT_HPP_
