#ifndef RAD_VET_RAYS_SC_SEGMENT_HPP_
#define RAD_VET_RAYS_SC_SEGMENT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_segment.hpp
//! \brief Short-characteristics segment kernels (vet_seg) and the vet_gd launch policies.
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

// a material edge of the upwind stencil: max/min chi > kVetCornerRatio (the corner rule
// of the curved core, ScCornerEdge)
constexpr Real kVetCornerRatio = 4.0;

#endif  // RAD_VET_RAYS_SC_SEGMENT_HPP_
