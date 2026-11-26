//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_oldvector.cpp
//! \brief the old vector of a stage, and the lagged donor form of the off-diagonal stress
//! flux.
//!
//!  called from: rad_vet_tasks.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_flux_seam.cpp, rad_vet_init.cpp,
//!    rad_vet_ledger.cpp

#include <cstdio>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "utils/p_ledger.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::OldVector
//! \brief E and gas state the stage starts from (stage increments added), with the gas
//! internal energy and temperature

void RadiationVet::OldVector() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  const bool feedback = cln_feedback;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fluid = fl.on;
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto u0_ = u0;
  auto inc = t2inc;
  auto w = cln_w;
  SolveArgs x = MakeArgs();

  // counted: gas eint < 0 at the stage start (taken as 0 by the exchange), E < 0, and a
  // non-finite start (a NaN max(e0, 0) would hide: the solve fails, the step is redone).
  // A flat range; the counts are integers, exact in any order.
  Real nn[3] = {0.0, 0.0, 0.0};
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("cln_old", SolveRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &ng, Real &ne, Real &nf) {
    const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
    const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
    Real t = 0.0;
    if (fluid) {
      if (feedback && stage) {
        for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += inc(m,VET_T2_M1+n,k,j,i);
        uh(m,IEN,k,j,i) += inc(m,VET_T2_EN,k,j,i);
      }
      const Real d = uh(m,IDN,k,j,i);
      Real p[3];
      GasMomentum(x, uh, m, k, j, i, p);
      const Real ekin = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/d;
      // a deficit eint < 0 stays in the gas row, paid by its implicit exchange (e(T) > 0)
      const Real eint = GasEint(x, uh, ekin, m, k, j, i);
      if (eint < 0.0) ng += 1.0;
      if (!(Kokkos::isfinite(eint) && Kokkos::isfinite(d))) nf += 1.0;
      t = eos.Temperature(d, fmax(eint, 0.0));
      w(m,CLW_EI,k,j,i) = eint;
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,VET_E,k,j,i) + (stage ? inc(m,VET_T2_E,k,j,i) : 0.0);
    if (w(m,CLW_EO,k,j,i) < 0.0) ne += 1.0;
    if (!Kokkos::isfinite(w(m,CLW_EO,k,j,i))) nf += 1.0;
  }, Kokkos::Sum<Real>(nn[0]), Kokkos::Sum<Real>(nn[1]), Kokkos::Sum<Real>(nn[2]));
  SumOverRanks(nn, 3);   // (counts: one reduction)
  // clean_cs_lz_abs = 2: the increment's L_z, of the slopes' sums (Time2FormStage)
  if (cln_lzbook && fluid && feedback && stage) LzBook(cln_lzc, true, cln_lzr[3]);
  // p_ledger: the increment's gas momentum withheld (cone cut, fwb, reference force)
  if (fluid && feedback && stage && !cln_vzero && !t2_inprobe) p_ledger::Start();
  if (global_variable::my_rank == 0 && nn[2] > 0.0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage start not finite in %.0f cells "
                "(gas or E; the solve fails and the step is redone)\n",
                pmy_pack->pmesh->ncycle, t2_solve, nn[2]);
  }
  if (global_variable::my_rank == 0 && cln_verbose > 0 && (nn[0] > 0.0 || nn[1] > 0.0)) {
    std::printf("CLEAN OLDVEC cycle %d t2s %d eint<0 %.0f (carried) E<0 %.0f\n",
                pmy_pack->pmesh->ncycle, t2_solve, nn[0], nn[1]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FaceDonor
//! \brief clean_stag_xdonor (staggered): the off-diagonal stress flux of each face,
//! C = -(1-w) theta h c^2 d_t(D_nt E) (FaceCross), of the iterate E_k enters the
//! pass as a donor-cell coefficient, C(E) -> q E_donor, q = C(E_k)/E_k,donor, |q| <= c.
//! q adds to al (>= 0) or ar (<= 0) only, so the operator stays an M-matrix for any D, h:
//! a non-negative right-hand side gives a non-negative solution, without a mask.  At the
//! fixed point q E_donor = C(E) where |C| <= c E_donor; beyond, the cross flux is limited
//! at free streaming.  The gas momentum keeps the central stress.
//! Lagged C alternates pass to pass where nothing else holds E (optically thin), so a
//! face whose change turns back, d = l x the last (l < 0), takes d/(1 - l), the fixed
//! point of that mode.  first: a solve's first iterate, no history.

void RadiationVet::FaceDonor(const bool first, const bool lat, const bool hold) {
  if (!pmy_pack->pmesh->multi_d) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  cln_xd_on = false;
  GhostDefer(cln_w, CLW_EK, 1, 1);   // one exchange, depth 1
  GhostFill(cln_w, CLW_EO, 1, 0, 1);   // donors' E_old of the lag split
  if (lat) CsLateralStress(CLW_EK);
  // hold (face_xdonor_freeze): keep C, q and the explicit part of the last pass
  cln_xdh = hold;
  if (hold) {cln_xd_on = true; return;}
  SolveArgs x = MakeArgs();   // x.xq = -1: the central cross flux
  const int xq = cln_xqs;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    // zero every face (one flat range); the pass's split (XSplit) sets the
    // active ones
    const int n4 = fc.extent_int(4), n34 = fc.extent_int(3)*n4;
    const int n234 = fc.extent_int(2)*n34;
    AthFor("cln_xdonor", Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0,
                         (nmb1 + 1)*n234),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/n234, k = (idx - m*n234)/n34;
      const int j = (idx - m*n234 - k*n34)/n4, i = idx%n4;
      for (int q = 0; q < 3; ++q) fc(m,xq+q,k,j,i) = 0.0;
      if (k < ks || k > ke+dk || j < js || j > je+dj || i < is || i > ie+di) return;
      const Real cn = FaceCross(x, CLW_EK, dir, m, k, j, i);
      const Real cp = fc(m,CLF_XC,k,j,i), d = cn - cp, dp = fc(m,CLF_XC+1,k,j,i);
      // l < 0 where C turns back: step d/(1 - l), the fixed point of a mode C -> l C;
      // continuous in l
      const Real l = (d*dp < 0.0) ? d/dp : 0.0;
      const Real cu = first ? cn : cp + d/(1.0 - l);
      fc(m,CLF_XC,k,j,i) = cu;
      fc(m,CLF_XC+1,k,j,i) = first ? 0.0 : cu - cp;
    });
  }
  cln_xd_on = true;
}

}  // namespace radvet
