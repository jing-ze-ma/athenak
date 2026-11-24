#ifndef RAD_VET_RAD_VET_VETSEG_HPP_
#define RAD_VET_RAD_VET_VETSEG_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetseg.hpp
//! \brief <rad_vet>/vet_seg = split, the short-characteristics segment from
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
//! th = (hu - h0)/(hu + h0) q (S_u > S_0 as emit) S_f = S_u + hu*(gf + th*((1 - e^-hu)
//! mm - gf)).  q = 1 - (S_uu - S_u)/(S_u - S_0) in [0, 1], S_uu the source of the cell
//! one further upwind (rad_vet_vgdray.hpp; q = 1 where the sweep has none): an
//! isothermal body (S_uu = S_u) or a single hot cell keeps S_u to the face, while a
//! resolved monotone temperature gradient over an unresolved optical depth (an old-prod
//! photosphere: tau per cell ~300 -> 0.005 within 2 cells, S_uu - S_u >= S_u - S_0)
//! is a radiative-diffusion skin, S linear in tau to the face (the mean face value; S_u
//! there made the HG/RSG restarts 16-47 % overluminous, STAGFIX5).  A thin edge cell
//! of an emitting body keeps its own S up to the face (emit
//! fell back to the mean ramp S_u -> S_0 there, losing half of the outer half-cell's
//! emission: the T6b core luminosity of the formal solution -6 % at tau_cell 0.13 and
//! not first order).  And where the OWN cell is the thicker and S_0 > S_u (a ray
//! entering an emitting body), the own part of the segment starts from
//! S_f0 = S_f + (h0 - hu)/(h0 + hu) (S_0 - S_f) instead of the continuous S_f (S jumps
//! at the face): the inward intensity of the edge cells (H +5 %, J -9 % at tau_cell
//! 0.13 with emit).  In a smooth field both contrasts are O(dx): second order kept.

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

// a range policy whose functor rides as a KERNEL ARGUMENT
// (HintLightWeight, as VetFor in rad_vet_vet.cpp) for the per-shell vet_gd sweep kernels
// (sp and cs: shell, pack, scatter/unpack); without it the 0.5-4 kB functors take the
// constant-memory launch with a host wait on the previous such launch (the same
// kernels, the same numbers)
inline auto VgdLight(const DevExeSpace &ex, const int a, const int b) {
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
}

// the vet_gd shell kernel's policy (VgdLight with launch bounds on HIP).
// Kokkos compiles a plain HIP RangePolicy kernel for up to 1024 threads per block, which
// caps it at 128 VGPRs on gfx942 and spills the shell body to scratch; it launches 256
// per block anyway.  LaunchBounds<256, M1_VGD_LB> lets the register allocator use up to
// 512 / M1_VGD_LB VGPRs.  The arithmetic is unchanged (register allocation only).
// 1 (was 2): R136 D6 apudev 13084234, vet_gd sweeps 5.40 -> 5.09 s, main
// loop medians 20.53 -> 20.33 s, output md5 identical
#ifndef M1_VGD_LB
#define M1_VGD_LB 1
#endif
inline auto VgdShellPol(const DevExeSpace &ex, const int a, const int b) {
#if defined(KOKKOS_ENABLE_HIP) && (M1_VGD_LB > 0)
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<Kokkos::LaunchBounds<256, M1_VGD_LB>>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
#else
  return VgdLight(ex, a, b);
#endif
}

// the CUDA gx-false shell kernel, launch bounds <T, B> (M1_VGD_CUDA_LB;
// 0 = the plain VgdShellPol launch).  default 2 = <128, 5> on A100 (96
// registers, 5 blocks per SM instead of 4; R136 r136.00097 4x1 -26 ms/cycle, hst
// identical); 0 elsewhere (not measured)
#if defined(KOKKOS_ENABLE_CUDA)
#ifndef M1_VGD_CUDA_LB
#if defined(KOKKOS_ARCH_AMPERE80)
#define M1_VGD_CUDA_LB 2
#else
#define M1_VGD_CUDA_LB 0
#endif
#endif
template <unsigned T, unsigned B>
inline auto VgdLbPol(const DevExeSpace &ex, const int a, const int b) {
  return Kokkos::Experimental::require(
      Kokkos::RangePolicy<Kokkos::LaunchBounds<T, B>>(ex, a, b),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
}
#endif

// vet_gd_body_global: copy a (host-side) kernel closure byte-wise into the
// device buffer buf (grown when too small) and return its device address.  The per-shell
// kernels then call the closure through this pointer instead of carrying it as a kernel
// argument: a ~2.7 kB argument whose operator() is called is copied by nvcc from the
// parameter space into every thread's local stack (cuobjdump STACK 2680 B per thread of
// the shell kernels).  The bytes are those Kokkos itself copies for a global-memory
// launch; deep_copy fences, so a buffer is never overwritten under a running kernel.
template <class F>
inline const F *VgdClosureToDevice(const F &f, DvceArray1D<char> &buf) {
  const size_t n = sizeof(F);
  if (buf.extent(0) < n) {buf = DvceArray1D<char>("vet_vgd_body", n);}
  using HostBytes = Kokkos::View<const char *, Kokkos::HostSpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  HostBytes h(reinterpret_cast<const char *>(&f), n);
  Kokkos::deep_copy(Kokkos::subview(buf, std::make_pair(static_cast<size_t>(0), n)), h);
  return reinterpret_cast<const F *>(buf.data());
}

// the fatal-error exit of the vet_gd setup (sp VgdFatal, cs VgdCsFatal): file, tag, msg
[[noreturn]] inline void VgdFatalAt(const char *file, const char *tag,
                                    const std::string &msg) {
  std::cout << "### FATAL ERROR in " << file << std::endl << tag << msg << std::endl;
  // (_Exit: no static destructors behind the live Kokkos views; they crashed std::exit)
  std::_Exit(EXIT_FAILURE);
}

// the vet_gd_replace surface q (vet_col_surface_q) of the top cell ie of
// column (m,k,j) (VetGdMoments): the geometry RH gives the projection a = n_d . r_hat
// of direction d (VgdRhatSp: the sp cell centre angles).  rr = r_v(ie)/r_f(ie + 1)
struct VgdRhatSp {
  Real st, ct, sp, cp;
  template <typename VD>
  KOKKOS_INLINE_FUNCTION Real A(const VD &dir_, const int d) const {
    return dir_(d,0)*st*cp + dir_(d,1)*st*sp + dir_(d,2)*ct;
  }
};
template <typename RH, typename VD, typename VI>
KOKKOS_INLINE_FUNCTION
Real VgdSurfaceQ(const RH &rh, const VD &dir_, const VI &vi_, const int n, const int m,
                 const int k, const int j, const int og, const int ie, const Real rr,
                 const bool sqf, const Real jca, const Real jcb, const Real jcap,
                 const Real qlo, const Real qhi, const Real q0) {
  Real jt = 0.0, jb = 0.0, ht = 0.0;
  for (int d = 0; d < n; ++d) {
    const Real a = rh.A(dir_, d);
    const Real w = dir_(d,3);
    jt += w*vi_(m,d,k+og,j+og,ie);
    jb += w*vi_(m,d,k+og,j+og,ie-1);
    ht += w*a*vi_(m,d,k+og,j+og,ie);
  }
  const Real hf = ht*rr*rr;
  const Real jq = sqf ? fmin(fmax(jca*jt - jcb*jb, 0.0), jcap*jt) : jt;
  return (jq > 0.0) ? fmin(fmax(hf/jq, qlo), qhi) : q0;
}

// the face source S_f of the split segment (see above); emit: 0 split, 1 emit, 2 jump,
// 3 jump without the plateau weights (q = q0 = 1: the caller passes no S_uu, S_dd)
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
      // coreL-1012 jump: the contrast moves S_f to the upwind cell's own value, as far
      // as that cell is a plateau: q = 1 - (S_uu - S_u)/(S_u - S_0) in [0, 1] (S_uu the
      // source one cell further upwind; q = 1 without one); only the slope estimate gb
      // (valid in an opaque cell) keeps the factor 1 - e^-hu
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

// the intensity at the face of the split segment (its upwind part)
KOKKOS_INLINE_FUNCTION
Real VetSegFaceI(const Real iu, const Real su, const Real s0, const Real hu,
                 const Real h0, const int emit) {
  const Real sf = VetSegFace(iu, su, s0, hu, h0, emit);
  Real ex, w0, wu;
  VetSegW(hu, ex, w0, wu);
  return iu*ex + wu*su + w0*sf;
}
// the intensity at the cell centre from I_u at the upwind point, the split segment;
// ifc (if given) takes the intensity at the face
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
    // coreL-1012 jump: the OWN cell the thicker and the hotter (a ray entering an
    // emitting body): its side of the face keeps S_0 by the contrast (S may jump at
    // the face; the upwind part keeps S_u -> S_f), as far as the own cell is a plateau:
    // q0 = 1 - (S_dd - S_0)/(S_0 - S_u) in [0, 1] (S_dd the source one cell further
    // downstream; q0 = 1 without one), the mirror of VetSegFace's q
    const Real q0 = (sdd > s0) ? fmax(1.0 - (sdd - s0)/(s0 - su), 0.0) : 1.0;
    const Real sf0 = sf + (h0 - hu)/(h0 + hu)*q0*(s0 - sf);
    return ia*ex0 + wu0*sf0 + w00*s0;
  }
  return ia*ex0 + wu0*sf + w00*s0;
}

// one SC segment from the upwind point (I_u, chi_u, S_u) to the
// cell centre (chi_0, S_0) over the path length hl = dx1/|mu1|, vet_seg split or mean;
// ifc (if given) takes the intensity at the face (half the path); suu, sdd: VetSegFace,
// VetSegSplit
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

// vet_interp: the upwind segment from the (bi)linear stencil of
// the upwind plane, weights w[4] (sum 1), corner intensities, extinctions, sources.
//   1 log:    I_u = sum w I_q, chi_u and S_u geometric (sum w ln), one segment (vet_gd's
//             interpolation);
//   2 corner: one segment per corner, I = sum w I(segment from corner q): the OUTCOME is
//             interpolated.  At a sharp opaque edge the plain interpolation (0 lin:
//             chi_u, S_u linear, one segment) puts the opaque chi_u on a ray that passes
//             the edge and replaces its intensity by the interpolated S_u (a ray grazing
//             an opaque disc of S = B loses the limb: L_fs 0.58-0.66 of exact,
//             resolution independent, scemu.py); corner keeps the thin corner's intensity
//             with its weight (0.95-1.0).  Smooth fields: corner - lin = O(dx^2).
//   3 auto:   corner where the stencil's extinction contrast max/min chi exceeds
//             kVetCornerRatio (a material edge), log elsewhere (a smooth field: the two
//             agree to O(dx^2) there, and log costs one segment instead of up to four)
// ifc (if given) takes the intensity at the face; suq (if given) the corners' S one
// x1 cell further upwind (vet_seg = jump, VetSegFace), interpolated as sq; sdd the own
// cell's S one x1 cell downstream (VetSegSplit).
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

// scheme = clean: the second-order foot of the transverse interpolation, the PLM remap
// of the upwind plane with van Leer slopes s.  The linear interpolation at fraction a
// between nodes a and b is the donor-cell remap, whose diffusion accumulates along a
// ray (first order in an optically thin multi-D field); the remap adds
// a (1 - a)/2 (s_a - s_b) along x2 (rows ka, kb weighted as in the bilinear value) and
// along x3.  The part a (1 - a)/2 s_a crosses the low face of the cell (fl[0] x2, fl[1]
// x3), so the update stays in flux form.  ip(dk, dj) = I at (ka + dk, ja + dj).
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

#endif  // RAD_VET_RAD_VET_VETSEG_HPP_
