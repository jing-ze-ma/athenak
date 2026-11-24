#ifndef RAD_M1_RAD_M1_ZALESAK_HPP_
#define RAD_M1_RAD_M1_ZALESAK_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_zalesak.hpp
//! \brief <rad_m1>/zalesak_converge: the Zalesak limiter R of time2_face_fct and
//! clean_wb_fct as the greatest fixed point of a monotone map, converged per rank between
//! R exchanges

#include <cfloat>
#include <cstdio>
#include <utility>

#include "globals.hpp"

namespace radm1 {

//! \brief the Zalesak rounds stopped at their cap with R still falling: the limited E is
//! then not bounded (never silent)
inline void ZalesakCapWarn(const char *who, const int cycle, const int t2s) {
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: %s Zalesak rounds at the cap 4096, R not "
                "converged (E above the floor not guaranteed)\n", cycle, t2s, who);
  }
}

//! \brief a cell's limiter bound: the largest R <= 1 with q - R po above the round-off
//! of the sums that later form its E.  q is what the cell holds without its outflow po
//! (its inflow at the donors' R), s the magnitudes of the terms of q and po at R = 1
//! (independent of R, so the bound is monotone in q).  E is a sum of at most 8 such terms
//! (the content, 6 faces, a floor), each addition off by at most eps of a partial sum
//! <= s, and R po by eps po <= eps s: 16 eps s covers both.  R = 0 where q is not above
//! it (no outflow; 'left' if q itself is not > 0); R = 1 without outflow.

KOKKOS_INLINE_FUNCTION
Real ZalesakBound(const Real q, const Real po, const Real s) {
  const Real b = q - 16.0*DBL_EPSILON*s;
  return (po > 0.0) ? ((b > 0.0) ? fmin(1.0, b/po) : 0.0) : 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::ZalesakConverge
//! \brief R = 1, then R_i = min(R_i, rf(m, k, j, i, c)) until no R falls, rf monotone in
//! the face neighbours' R (read from slot c with fabs).  The map is monotone, so from
//! R = 1 every schedule ends at its greatest fixed point: the result does not depend on
//! ranks, layout or round order.  Each rank runs Jacobi rounds on its own cells with the
//! ghosts held (an R that fell carries its sign bit; a round evaluates only the cells
//! next to one) until none falls; then one count over the ranks of the falls within 2 ng
//! of a MeshBlock face.  None: the ghosts are current, so R is the fixed point.  Else one
//! exchange; a face ghost that fell takes the sign bit.  The R end in slot
//! Time2FctSlot(0), face ghosts included, |R|; returns the falls (all ranks); q indexes
//! the counters (zl_halos: exchanges, zl_sub: local rounds).

template <class F>
Real RadiationM1::ZalesakConverge(const F &rf, const int q) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1, zd = 2*id.ng;
  const bool x2 = pmy_pack->pmesh->multi_d, x3 = pmy_pack->pmesh->three_d;
  auto iw_ = iw;
  const int ca = Time2FctSlot(0), cb = Time2FctSlot(1);
  auto slot = [&](const int c) {
    return Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  };
  Kokkos::deep_copy(DevExeSpace(), slot(ca), 1.0);   // (physical ghosts stay 1)
  Kokkos::deep_copy(DevExeSpace(), slot(cb), 1.0);
  const int n1 = ie - is + 1, n2 = je - js + 1, n3 = ke - ks + 1;
  const int nji = n2*n1, nkji = n3*nji;
  const int f1 = n3*n2, f2 = x2 ? n3*n1 : 0, f3 = x3 ? n2*n1 : 0, nf = 2*(f1 + f2 + f3);
  int sa = ca, sb = cb;
  bool all = true;
  Real nfall = 0.0;
  bool done = false;
  for (int r = 0; r < 4096 && !done; ++r) {
    Real nb = 0.0;
    for (int s = 0; s < 4096; ++s) {   // (local rounds, ghosts held)
      Real nv = 0.0, nvb = 0.0;
      const int s0 = sa, s1 = sb;
      const bool a0 = all;
      Kokkos::parallel_reduce("m1_zlc_r", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1 + 1)*nkji),
      KOKKOS_LAMBDA(const int idx, Real &a, Real &b) {
        const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
        const int i = is + idx%n1;
        const Real rc = fabs(iw_(m,s0,k,j,i));
        bool act = a0;
        for (int d = 0; d < 3 && !act; ++d) {
          if ((d == 1 && !x2) || (d == 2 && !x3)) continue;
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          act = Kokkos::signbit(iw_(m,s0,k-dk,j-dj,i-di)) ||
                Kokkos::signbit(iw_(m,s0,k+dk,j+dj,i+di));
        }
        Real rn = rc;
        if (act) {
          const Real f = rf(m, k, j, i, s0);
          if (f < rc) {
            rn = -f;
            a += 1.0;
            if (i < is + zd || i > ie - zd || (x2 && (j < js + zd || j > je - zd)) ||
                (x3 && (k < ks + zd || k > ke - zd))) {b += 1.0;}
          }
        }
        iw_(m,s1,k,j,i) = rn;
      }, Kokkos::Sum<Real>(nv), Kokkos::Sum<Real>(nvb));
      all = false;
      std::swap(sa, sb);
      zl_sub[q] += 1.0;
      nfall += nv;
      nb += nvb;
      if (!(nv > 0.0)) break;
    }
    if (!(CleanReduce(nb, false) > 0.0)) {
      done = true;
      break;
    }
    if (sa != ca) {   // (every rank exchanges slot ca: the ranks' round counts differ)
      Kokkos::deep_copy(DevExeSpace(), slot(ca), slot(sa));
      std::swap(sa, sb);
    }
    Time2FctHalo(sa);
    zl_halos[q] += 1.0;
    // (each face ghost of both slots: the exchanged |R|, its sign bit set if it fell)
    const int s0 = sa, s1 = sb;
    Kokkos::parallel_for("m1_zlc_g", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         (nmb1 + 1)*nf),
    KOKKOS_LAMBDA(const int t) {
      const int m = t/nf;
      int u = t - m*nf, k, j, i;
      if (u < 2*f1) {
        i = (u < f1) ? is - 1 : ie + 1;
        u %= f1;
        k = ks + u/n2;
        j = js + u%n2;
      } else if (u < 2*(f1 + f2)) {
        u -= 2*f1;
        j = (u < f2) ? js - 1 : je + 1;
        u %= f2;
        k = ks + u/n1;
        i = is + u%n1;
      } else {
        u -= 2*(f1 + f2);
        k = (u < f3) ? ks - 1 : ke + 1;
        u %= f3;
        j = js + u/n1;
        i = is + u%n1;
      }
      const Real g = fabs(iw_(m,s0,k,j,i));
      iw_(m,s0,k,j,i) = iw_(m,s1,k,j,i) = (g < fabs(iw_(m,s1,k,j,i))) ? -g : g;
    });
  }
  if (!done) {ZalesakCapWarn(q ? "clean_wb_fct" : "time2_face_fct",
                             pmy_pack->pmesh->ncycle, t2_solve);}
  auto r0 = slot(sa);
  auto r1 = slot(ca);
  par_for("m1_zlc_abs", DevExeSpace(), 0, nmb1, 0, r0.extent_int(1) - 1, 0,
          r0.extent_int(2) - 1, 0, r0.extent_int(3) - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    r1(m,k,j,i) = fabs(r0(m,k,j,i));
  });
  return CleanReduce(nfall, false);
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_ZALESAK_HPP_
