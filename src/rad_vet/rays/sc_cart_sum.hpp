#ifndef RAD_VET_RAYS_SC_CART_SUM_HPP_
#define RAD_VET_RAYS_SC_CART_SUM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_sum.hpp
//! \brief the Cartesian slab's per-cell face sums (VetScSumGroup, one thread per cell and
//! group of sums), every sum added in ray order.
//!
//!  called from: rays/sc_cart_setup.cpp

#include "athena.hpp"
#include "rad_vet/rad_vet.hpp"

namespace radvet {

// kinds of a per-cell sum, as passed to a writer w(kind, u, n, value) (u = 0 up, 1 down)
enum VetScSum {VSS_J, VSS_HP, VSS_KP, VSS_NL, VSS_F2, VSS_FN};

// groups of the ray-parallel sum: moments J..9 in two halves, the half-range and K
// sums with the non-locality, and the face sums by their (1, mu1, mu2, mu3) component
constexpr int VETSUM_NG = 7;

// the moment term n (VETRAY_J + n) of a ray of weighted intensity a
KOKKOS_INLINE_FUNCTION
Real VetScMom(const int n, const Real a, const Real m1, const Real m2, const Real m3) {
  switch (n) {
    case 0: return a;
    case 1: return a*m1*m1;
    case 2: return a*m2*m2;
    case 3: return a*m3*m3;
    case 4: return a*m1*m2;
    case 5: return a*m1*m3;
    case 6: return a*m2*m3;
    case 7: return a*m1;
    case 8: return a*m2;
    default: return a*m3;
  }
}

// face sums n of the up / down rays to the writer: f2 then fn of each (mbo, the banded
// sweep), else both f2 then both fn (one writer call order per form: same bits)
template <class W>
KOKKOS_INLINE_FUNCTION
void VetScFaceW(const W &w, const int n, const Real f2u, const Real f2d, const Real fnu,
                const Real fnd, const bool mbo) {
  w(VSS_F2, 0, n, f2u);
  if (mbo) {w(VSS_FN, 0, n, fnu);}
  w(VSS_F2, 1, n, f2d);
  if (!mbo) {w(VSS_FN, 0, n, fnu);}
  w(VSS_FN, 1, n, fnd);
}

//----------------------------------------------------------------------------------------
//! \fn VetScSumGroup
//! \brief group G of the sums of one cell over rays [0, nray): ray r with active(r)
//! (up rays m1 > 0 to u = 0), intensity iv(r), own-cell source s0(r) (non-locality) and
//! face-sum scalar g(q, r); written to w (mbo: up f2 fn, then down f2 fn).  Each
//! group holds at most 26 sums: VETSUM_NG threads per cell, few registers each.
template <int G, class A, class FA, class FI, class FS, class FG, class W>
KOKKOS_INLINE_FUNCTION
void VetScSumGroup(const A &ang, const int nray, const FA &active, const FI &iv,
                   const FS &s0, const FG &g, const W &w, const bool hrs,
                   const bool kps, const bool nls, const bool f2s, const bool mbo) {
  if (G <= 1) {   // moments VETRAY_J + 5 G .. + 4
    Real s[2][5] = {};
    for (int r = 0; r < nray; ++r) {
      if (!active(r)) continue;
      const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2), wr = ang(r,3);
      const int u = (m1 > 0.0) ? 0 : 1;
      const Real a = wr*iv(r);
      for (int n = 0; n < 5; ++n) {s[u][n] += VetScMom(5*G + n, a, m1, m2, m3);}
    }
    for (int n = 0; n < 5; ++n) {
      w(VSS_J, 0, 5*G + n, s[0][n]);
      w(VSS_J, 1, 5*G + n, s[1][n]);
    }
  } else if (G == 2) {   // half-range sums, K sums, non-locality
    if (!hrs && !kps && !nls) return;
    Real hp[2][3] = {}, kp[2][VETRAY_NKP] = {}, nl[2] = {};
    for (int r = 0; r < nray; ++r) {
      if (!active(r)) continue;
      const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2), wr = ang(r,3);
      const int u = (m1 > 0.0) ? 0 : 1;
      const Real ir = iv(r);
      const Real a = wr*ir;
      if (hrs) {
        hp[u][0] += a*fmax(m1, 0.0);
        hp[u][1] += a*fmax(m2, 0.0);
        hp[u][2] += a*fmax(m3, 0.0);
      }
      if (kps) {
        const Real mu[3] = {m1, m2, m3};
        for (int d = 0; d < 3; ++d) {
          if (!(mu[d] > 0.0)) continue;
          for (int b = 0; b < 3; ++b) kp[u][3*d+b] += a*mu[b]*mu[d];
        }
      }
      if (nls) {nl[u] += wr*fabs(ir - s0(r));}
    }
    if (hrs) {
      for (int n = 0; n < 3; ++n) {w(VSS_HP, 0, n, hp[0][n]); w(VSS_HP, 1, n, hp[1][n]);}
    }
    if (kps) {
      for (int n = 0; n < VETRAY_NKP; ++n) {
        w(VSS_KP, 0, n, kp[0][n]);
        w(VSS_KP, 1, n, kp[1][n]);
      }
    }
    if (nls) {w(VSS_NL, 0, 0, nl[0]); w(VSS_NL, 1, 0, nl[1]);}
  } else {   // face sums n = c, c + 4, c + 8, c + 12 of component c = G - 3
    if (!f2s) return;
    constexpr int c = G - 3;
    Real f2[2][4] = {}, fn[2][4] = {};
    for (int r = 0; r < nray; ++r) {
      if (!active(r)) continue;
      const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2);
      const bool up = (m1 > 0.0);
      const int u = up ? 0 : 1;
      auto cm = [&](const Real q) {
        return (c == 0) ? q : ((c == 1) ? q*m1 : ((c == 2) ? q*m2 : q*m3));
      };
      f2[u][1] += cm(g(0, r));
      fn[u][0] += cm(g(1, r));
      fn[u][3] += cm(g(2, r));
      if (up) {f2[0][2] += cm(g(3, r));} else {fn[1][2] += cm(g(3, r));}
      f2[u][1] += cm(g(4, r));
    }
    for (int b = 0; b < 4; ++b) {
      VetScFaceW(w, 4*b + c, f2[0][b], f2[1][b], fn[0][b], fn[1][b], mbo);
    }
  }
}

// all groups of one cell: thread (cell, grp) of a ray-parallel sum kernel
template <class A, class FA, class FI, class FS, class FG, class W>
KOKKOS_INLINE_FUNCTION
void VetScSum(const int grp, const A &ang, const int nray, const FA &active,
              const FI &iv, const FS &s0, const FG &g, const W &w, const bool hrs,
              const bool kps, const bool nls, const bool f2s, const bool mbo = false) {
  switch (grp) {
    case 0: VetScSumGroup<0>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    case 1: VetScSumGroup<1>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    case 2: VetScSumGroup<2>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    case 3: VetScSumGroup<3>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    case 4: VetScSumGroup<4>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    case 5: VetScSumGroup<5>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                             mbo);
      break;
    default: VetScSumGroup<6>(ang, nray, active, iv, s0, g, w, hrs, kps, nls, f2s,
                              mbo);
  }
}

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_CART_SUM_HPP_
