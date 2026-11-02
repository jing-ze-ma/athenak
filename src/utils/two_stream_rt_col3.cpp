//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt_col3.cpp
//! \brief the rt_implicit_column = 3 column-solve launchers (RTCol3Launch,
//! RTCol3TeamLaunch).
//! tusplit-1102: moved verbatim out of the two-stream headers (which keep the
//! declarations and the doc comments) so that they are compiled once, in their
//! own translation unit.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "utils/two_stream_rt.hpp"

namespace two_stream_rt {

void RTCol3Launch(const RTCol3 &c, const int nmb1, const int ks, const int ke,
                  const int js, const int je) {
  RTCol3 cc = c;
  par_for("rt_col3", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    cc.Solve(m, k, j);
  });
}

void RTCol3TeamLaunch(const RTCol3 &c, const int nmb1, const int ks, const int ke,
                      const int js, const int je) {
  RTCol3 cc = c;
  const int nk = ke - ks + 1, nj = je - js + 1;
  const int nteam = (nmb1 + 1)*nk*nj;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nteam, cc.nseg);
  // problem/rt_impl_mixed picks the factor type; the two extra instantiations exist
  // only when the switch is set, and mixed = 0 is the untouched double kernel.
  if (cc.mixed == 1) {
    Kokkos::parallel_for("rt_col3_p", policy, KOKKOS_LAMBDA(TeamMember_t tm) {
      const int q = tm.league_rank();
      const int m = q/(nk*nj);
      const int k = ks + (q - m*nk*nj)/nj;
      const int j = js + (q - m*nk*nj)%nj;
      RTCol3TeamSolve<float, false>(cc, tm, m, k, j);
    });
  } else if (cc.mixed == 2) {
    Kokkos::parallel_for("rt_col3_p", policy, KOKKOS_LAMBDA(TeamMember_t tm) {
      const int q = tm.league_rank();
      const int m = q/(nk*nj);
      const int k = ks + (q - m*nk*nj)/nj;
      const int j = js + (q - m*nk*nj)%nj;
      RTCol3TeamSolve<float, true>(cc, tm, m, k, j);
    });
  } else {
    Kokkos::parallel_for("rt_col3_p", policy, KOKKOS_LAMBDA(TeamMember_t tm) {
      const int q = tm.league_rank();
      const int m = q/(nk*nj);
      const int k = ks + (q - m*nk*nj)/nj;
      const int j = js + (q - m*nk*nj)%nj;
      RTCol3TeamSolve<Real, false>(cc, tm, m, k, j);
    });
  }
}

}  // namespace two_stream_rt
