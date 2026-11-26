//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_fct.cpp
//! \brief limited face energy fluxes of the stages (time2_face_fct, clean_wb_fct) and the
//! energy flux's memory d.
//!
//!  called from: rad_vet_imex.cpp, rad_vet_ledger.cpp, rad_vet_restart.cpp,
//!    rad_vet_update.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_init.cpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxMemoryKeep
//! \brief clean_stag_fkin: the energy flux memory d (CLF_N + 3) at the stage-loop start
//! into CLF_DS (save) or back (a rejected stage's write-backs are discarded)

void RadiationVet::FluxMemoryKeep(const bool save) {
  if (!cln_fsolve || !cln_fkin) return;
  if (!cln_init) InitArrays();
  auto AL = Kokkos::ALL;
  for (auto fc : {cln_fc1, cln_fc2, cln_fc3}) {
    auto d = Kokkos::subview(fc.s, AL, CLF_N + 3, AL, AL, AL);
    auto s = Kokkos::subview(fc.s, AL, CLF_DS, AL, AL, AL);
    if (save) {Kokkos::deep_copy(s, d);} else {Kokkos::deep_copy(d, s);}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxMemoryInc
//! \brief clean_stag_fkin: d of a hesdirk2 stage's old vector: the stage start state
//! (stage 2: w1 d^n + w0 Y1) plus the increment c1 K1 + c2 K2 of d's own slopes (CLF_KD,
//! CLF_KD + 1), as every other component of the integrated state

void RadiationVet::FluxMemoryInc(const bool two, const Real c1, const Real c2,
                                 const Real w1, const Real w0) {
  if (!cln_fsolve || !cln_fkin) return;
  for (auto fc : {cln_fc1, cln_fc2, cln_fc3}) {
    par_for("cln_dinc", DevExeSpace(), 0, fc.extent_int(0) - 1, 0, fc.extent_int(2) - 1,
            0, fc.extent_int(3) - 1, 0, fc.extent_int(4) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real d = fc(m,CLF_N+3,k,j,i) + c1*fc(m,CLF_KD,k,j,i);
      if (two) d += w1*fc(m,CLF_DS,k,j,i) + (w0 - 1.0)*fc(m,CLF_N+3,k,j,i) +
                    c2*fc(m,CLF_KD+1,k,j,i);
      fc(m,CLF_N+3,k,j,i) = d;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxMemoryRetry
//! \brief dt_retry: the memory d and its slope K1 (CLF_N + 3, CLF_KD), the only face
//! values a cycle takes from the last one, to / from cln_rstd (the restart's layout)

void RadiationVet::FluxMemoryRetry(const bool save) {
  const int nmb = cln_fc1.extent_int(0);
  if (!cln_fsolve || !cln_fkin || nmb == 0) return;
  auto AL = Kokkos::ALL;
  if (save && (cln_rstd.extent_int(0) != nmb || cln_rstd.extent_int(1) != 6)) {
    Kokkos::realloc(cln_rstd, nmb, 6, cln_fc1.extent_int(2), cln_fc1.extent_int(3),
                    cln_fc1.extent_int(4));
  }
  for (int e = 0; e < cln_rstd.extent_int(1) && cln_rstd.extent_int(0) == nmb; ++e) {
    auto d = Kokkos::subview((e%3 == 0) ? cln_fc1.s : ((e%3 == 1) ? cln_fc2.s :
                             cln_fc3.s), AL, (e < 3) ? CLF_N + 3 : CLF_KD, AL, AL, AL);
    auto r = Kokkos::subview(cln_rstd, AL, e, AL, AL, AL);
    if (save) {Kokkos::deep_copy(r, d);} else {Kokkos::deep_copy(d, r);}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::WbFct
//! \brief clean_wb_fct: Zalesak limit of the written faces' energy flux Phi (cln_wbdiv's,
//! fp 0..2) so that E > e_floor after the write-back by construction.  A cell's E is
//! Q - P: Q its E without face transport (cln_wbcell probe, fp 3, less e_floor and
//! round-off) plus inflow, P its outflow.  A cell with Q - R P not > 0 takes R = min(1,
//! (1 - 1e-12) Q/P) (0, counted 'left', where Q <= 0); each face is scaled by its donor's
//! R; rounds until no R falls (an R that fell carries its sign bit, so the next round
//! evaluates only its neighbours).  The limited faces' (1 - s) Phi leave qf (CLW_Q), pf
//! and the x1 boundary ledger.  Face F is not touched: momentum and L_z exact.

void RadiationVet::WbFct(const SolveArgs &x, const SolveArgs *px,
                         const DvceArray5D<Real> &fp, const bool pfo,
                         const DvceArray5D<Real> &pf) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1, nf = x.nf;
  if (iw.extent_int(1) <= Time2FctSlot(1)) {
    RadVetFatal("<rad_vet>/clean_wb_fct needs time2_face_fct's iw slots");
  }
  auto iw_ = iw;
  const int ca = Time2FctSlot(0), cb = Time2FctSlot(1);
  for (const int c : {ca, cb}) {   // R = 1, ghosts too: physical ones stay unlimited
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), 1.0);
  }
  // a cell's Q and P at the donors' R in slot c; a face's scale, its donor's R in ca
  auto qp = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, const int c,
                          Real &q, Real &po) {
    const SolveArgs &x = *px;
    q = fp(m,3,k,j,i);
    po = 0.0;
    for (int n = 0; n < nf; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      const Real gl = x.h*DivW(x, n, m, k, j, i, 0)*fp(m,n,k,j,i);   // + = in
      const Real gh = x.h*DivW(x, n, m, k, j, i, 1)*fp(m,n,k+dk,j+dj,i+di);
      if (gl > 0.0) {q += fabs(iw_(m,c,k-dk,j-dj,i-di))*gl;} else {po -= gl;}
      if (gh > 0.0) {po += gh;} else {q -= fabs(iw_(m,c,k+dk,j+dj,i+di))*gh;}
    }
  };
  auto sf = KOKKOS_LAMBDA(const int n, const int m, const int k, const int j,
                          const int i) {
    const Real f = fp(m,n,k,j,i);
    return (f > 0.0) ? iw_(m,ca,k-(n == 2),j-(n == 1),i-(n == 0)) :
                       ((f < 0.0) ? iw_(m,ca,k,j,i) : 1.0);
  };
  const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
  int nr = 0, sa = ca, sb = cb;
  for (; nr < 4096; ++nr) {   // R only falls: until none newly limited
    Real nv = 0.0;
    const bool all = (nr == 0);
    Kokkos::parallel_reduce("cln_wbfr", SolveRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &a) {
      const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
      const int i = is + idx%n1;
      const Real rc = fabs(iw_(m,sa,k,j,i));
      bool act = all;
      for (int n = 0; n < nf && !act; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        act = Kokkos::signbit(iw_(m,sa,k-dk,j-dj,i-di)) ||
              Kokkos::signbit(iw_(m,sa,k+dk,j+dj,i+di));
      }
      Real r = rc;
      if (act) {
        Real q, po;
        qp(m, k, j, i, sa, q, po);
        const Real rq = (q > 0.0 && po > 0.0) ? fmin(1.0, (1.0 - 1.0e-12)*q/po) : 0.0;
        if (!(q - rc*po > 0.0) && rq < rc) {
          r = -rq;
          a += 1.0;
        }
      }
      iw_(m,sb,k,j,i) = r;
    }, Kokkos::Sum<Real>(nv));
    Time2FctHalo(sb);
    std::swap(sa, sb);
    if (!(GlobalReduce(nv, false) > 0.0)) break;
  }
  if (nr == 4096) {
    ZalesakCapWarn("clean_wb_fct", pmy_pack->pmesh->ncycle, t2_solve);
  }
  if (sa != ca) {   // the last round's R, no sign bit set, ghosts too
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw, Kokkos::ALL, sa,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (nr == 0) return;   // nothing limited: bitwise
  Real lv[2] = {0.0, 0.0};   // 'left': not > 0 at R; the outflow the limit kept
  Kokkos::parallel_reduce("cln_wbfl", SolveRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a, Real &b) {
    const SolveArgs &x = *px;
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    Real q, po;
    qp(m, k, j, i, ca, q, po);
    const Real rc = iw_(m,ca,k,j,i);
    if (!(q - rc*po > 0.0)) a += 1.0;
    b += (1.0 - rc)*po*CellVolume(x, m, k, j, i);
    Real dq = 0.0;   // qf less the limited faces' share
    for (int n = 0; n < nf; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      dq += DivF(x, n, m, k, j, i, (1.0 - sf(n, m, k, j, i))*fp(m,n,k,j,i),
                 (1.0 - sf(n, m, k+dk, j+dj, i+di))*fp(m,n,k+dk,j+dj,i+di));
    }
    x.w(m,CLW_Q,k,j,i) -= x.h*dq;
  }, Kokkos::Sum<Real>(lv[0]), Kokkos::Sum<Real>(lv[1]));
  lv[0] = GlobalReduce(lv[0], false);
  lv[1] = GlobalReduce(lv[1], false);
  if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(lv[1]))) {
    std::printf("CLEAN %sWBFCT cycle %d t2s %d rounds %d left %.0f kept %.6e\n",
                std::isfinite(lv[1]) ? "" : "WARNING ", pmy_pack->pmesh->ncycle,
                t2_solve, nr, lv[0], lv[1]);
  }
  for (int n = 0; n < nf && pfo; ++n) {   // the stored Phi: what moved E
    SolveFor("cln_wbfs", nmb1, ks, ke+(n == 2), js, je+(n == 1), is, ie+(n == 0),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pf(m,n,k,j,i) -= (1.0 - sf(n, m, k, j, i))*fp(m,n,k,j,i);
    });
  }
  if (cln_lacc.extent_int(0) > nmb1) {   // the history's boundary energy (cln_wbbnd)
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;
    const Real h = x.h;
    par_for("cln_wbfb", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      if (IsBoundaryFace(*px, 0, m, k, j, i)) {
        la(m,s,k,j) -= h*ar(m,k,j,i)*(1.0 - sf(0, m, k, j, i))*fp(m,0,k,j,i);
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FctHalo
//! \brief exchange of the limiter R (iw slot c): first ghost layer only, the clean direct
//! lists, without the deferred clean exchanges

void RadiationVet::Time2FctHalo(const int c) {
  std::vector<GSeg<Real>> hp;
  hp.swap(cln_hpend);
  GhostFill(iw, c, 1, false, 0, 1);
  cln_hpend.swap(hp);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FaceFct
//! \brief time2_face_fct: Zalesak positivity limiter on the face fluxes of the stage
//! increment t2inc (E and face F channels).  G_f = A_f (c1 t2pf1 + c2 t2pf2) is the
//! energy the increment carries across face f; only cells whose old vector E + dE is not
//! > 0 are limited, by R_i = (1 - 1e-12) Q_i/P_i, and a face takes its donor's R.  Rounds
//! repeat until no cell is newly limited (R only decreases; cap 4096).  Bitwise where no
//! face is limited.

void RadiationVet::Time2FaceFct(int stage, Real c1, Real c2) {
  const bool two = (stage != 1);
  auto rhalo = [&](const int c) {Time2FctHalo(c);};
  // no stored face flux for K1 (first step, restart): K1 is the last stage-2 slope, whose
  // solution is U^n, so its face flux is the stage-start t2f* = f0x*
  if (!two && !t2pf1_ok && t2pf1.extent_int(0) >= pmy_pack->nmb_thispack &&
      t2f1.extent_int(0) >= pmy_pack->nmb_thispack) {
    auto pf_ = t2pf1;
    auto g1_ = t2f1, g2_ = t2f2, g3_ = t2f3;
    const bool q2 = trans_on && t2f2.extent_int(0) > 0;
    const bool q3 = trans_on && trans_x3 && t2f3.extent_int(0) > 0;
    const Real crp = chat/c_light;
    par_for("vet_t2_ffct_fill", DevExeSpace(), 0, pf_.extent_int(0) - 1,
            0, pf_.extent_int(2) - 1, 0, pf_.extent_int(3) - 1, 0, pf_.extent_int(4) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pf_(m,0,k,j,i) = crp*g1_(m,k,j,i);
      pf_(m,1,k,j,i) = q2 ? crp*g2_(m,k,j,i) : 0.0;
      pf_(m,2,k,j,i) = q3 ? crp*g3_(m,k,j,i) : 0.0;
    });
    t2pf1_ok = true;
    t2_ffct_fill += 1.0;
  }
  if (!t2pf1_ok || (two && !t2pf2_ok) || t2pf1.extent_int(0) < pmy_pack->nmb_thispack) {
    t2_ffct_nop += 1.0;  // no stored face flux for this slope
    return;
  }
  if (iw.extent_int(1) <= Time2FctSlot(0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time2_face_fct needs the Krylov work slots"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  t2_ffct_calls += 1.0;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto in_ = t2inc;
  auto u0_ = u0;
  auto p1_ = t2pf1;
  auto p2_ = t2pf2;
  auto iw_ = iw;
  const int ca = Time2FctSlot(0), cb = Time2FctSlot(1);
  const bool tr2 = trans_on, tr3 = trans_on && trans_x3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto vol = pmy_pack->pcoord->volume;
  auto ar1 = pmy_pack->pcoord->area.x1f;
  auto ar2 = pmy_pack->pcoord->area.x2f;
  auto ar3 = pmy_pack->pcoord->area.x3f;
  const bool curv = (cs_geom || sph_geom) && vol.extent_int(0) > nmb1 &&
                    vol.extent_int(3) > 1 && ar1.extent_int(0) > nmb1;
  const bool hoi = impl_hold && (hold_m.extent_int(0) > nmb1);
  auto hm_ = hold_m;
  const Real cc1 = c1, cc2 = two ? c2 : 0.0;
  // energy across the lower face of cell (k,j,i) in direction d (+ = towards +x_d)
  auto gf = KOKKOS_LAMBDA(const int m, const int d, const int k, const int j,
                          const int i) {
    const Real phi = cc1*p1_(m,d,k,j,i) + cc2*p2_(m,d,k,j,i);
    Real a;
    if (curv) {
      a = (d == 0) ? ar1(m,k,j,i) : ((d == 1) ? ar2(m,k,j,i) : ar3(m,k,j,i));
    } else {
      const auto &s = mbsize.d_view(m);
      a = (d == 0) ? s.dx2*s.dx3 : ((d == 1) ? s.dx1*s.dx3 : s.dx1*s.dx2);
    }
    return a*phi;
  };
  auto cvol = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const auto &s = mbsize.d_view(m);
    return curv ? vol(m,k,j,i) : s.dx1*s.dx2*s.dx3;
  };
  auto un_ = u1;   // (E^n: the limiter keeps th V E^n in a limited cell)
  const Real eth = t2_ffct_emin;
  // w = V (E + dE) unlimited, outflow po, inflow pi, the inflow pl the donors' R let
  // through, the change dw of V E_new by R (slot sa); lim: an adjacent face is limited
  auto cell = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            Real &w, Real &po, Real &pi, Real &pl, Real &dw,
                            bool &lim, const int sa) {
    const Real rc = fabs(iw_(m,sa,k,j,i));
    po = 0.0;
    pi = 0.0;
    pl = 0.0;
    dw = 0.0;
    lim = false;
    for (int d = 0; d < 3; ++d) {
      if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
      const int dk = (d == 2) ? 1 : 0, dj = (d == 1) ? 1 : 0, di = (d == 0) ? 1 : 0;
      // lower face: + = inflow, donor the lower neighbour
      Real g = gf(m, d, k, j, i);
      if (g > 0.0) {
        pi += g;
        const Real rn = fabs(iw_(m,sa,k-dk,j-dj,i-di));
        pl += rn*g;
        if (rn < 1.0) {dw -= (1.0 - rn)*g; lim = true;}
      } else if (g < 0.0) {
        po -= g;
        if (rc < 1.0) {dw -= (1.0 - rc)*g; lim = true;}
      }
      // upper face: + = outflow, donor this cell
      g = gf(m, d, k+dk, j+dj, i+di);
      if (g > 0.0) {
        po += g;
        if (rc < 1.0) {dw += (1.0 - rc)*g; lim = true;}
      } else if (g < 0.0) {
        pi -= g;
        const Real rn = fabs(iw_(m,sa,k+dk,j+dj,i+di));
        pl -= rn*g;
        if (rn < 1.0) {dw += (1.0 - rn)*g; lim = true;}
      }
    }
    w = cvol(m,k,j,i)*(u0_(m,VET_E,k,j,i) + in_(m,VET_T2_E,k,j,i));
  };
  // R = 1 everywhere incl. ghosts (physical ghosts stay unlimited); rounds read slot sa
  // and write sb, then swap (both slots keep the same ghosts)
  for (const int c : {ca, cb}) {
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), 1.0);
  }
  int sa = ca, sb = cb;
  Real ncell = 0.0;
  int nround = 0;
  // exchange R only after a round that changed a cell within 2 ng of a
  // MeshBlock face; the first exchange always runs
  const int zd = 2*id.ng;
  const bool zx2 = pmy_pack->pmesh->multi_d, zx3 = pmy_pack->pmesh->three_d;
  bool zfirst = true;
  // an R that fell is written with its sign bit set; the next round
  // evaluates only cells with such a face neighbour, the others keep |R|
  // no exchange: the face ghosts of the slot read next are those of the
  // other slot
  auto fghost = [&](const int c, const int cf) {
    const int n1 = ie - is + 1, n2 = je - js + 1, n3 = ke - ks + 1;
    const int f1 = n3*n2, f2 = tr2 ? n3*n1 : 0, f3 = tr3 ? n2*n1 : 0;
    const int nf = 2*(f1 + f2 + f3);
    Kokkos::parallel_for("vet_t2_ffct_fg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         (nmb1 + 1)*nf),
    KOKKOS_LAMBDA(const int t) {
      const int m = t/nf;
      int u = t - m*nf, k, j, i;
      if (u < 2*f1) {
        i = (u < f1) ? is - 1 : ie + 1;
        u %= f1;
        k = ks + u/n2;
        j = js + u % n2;
      } else if (u < 2*(f1 + f2)) {
        u -= 2*f1;
        j = (u < f2) ? js - 1 : je + 1;
        u %= f2;
        k = ks + u/n1;
        i = is + u % n1;
      } else {
        u -= 2*(f1 + f2);
        k = (u < f3) ? ks - 1 : ke + 1;
        u %= f3;
        j = js + u/n1;
        i = is + u % n1;
      }
      iw_(m,c,k,j,i) = fabs(iw_(m,cf,k,j,i));
    });
  };
  for (int r = 0; r < 4096; ++r) {   // until none newly limited
    Real nv[2] = {0.0, 0.0};
    const bool zsel = r > 0;   // (every cell in round 0)
    const int s0 = sa, s1 = sb;
    auto rnd = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                             Real &a, Real &ab) {
      const Real rc = fabs(iw_(m,s0,k,j,i));
      Real rn = rc;
      bool act = !zsel;   // (a face neighbour's R fell)
      for (int d = 0; d < 3 && !act; ++d) {
        if ((d == 1 && !tr2) || (d == 2 && !tr3)) {continue;}
        const int dk = (d == 2), dj = (d == 1), di = (d == 0);
        act = Kokkos::signbit(iw_(m,s0,k-dk,j-dj,i-di)) ||
              Kokkos::signbit(iw_(m,s0,k+dk,j+dj,i+di));
      }
      if (act && !(hoi && hm_(m,k,j,i) > 0.0)) {
        Real w, po, pi, pl, dw;
        bool lim;
        cell(m, k, j, i, w, po, pi, pl, dw, lim, s0);
        const Real wf = eth*cvol(m,k,j,i)*un_(m,VET_E,k,j,i);
        if (!(w + dw > wf)) {
          // held without face transport + remaining inflow, above th V E^n
          const Real q = w + po - pi + pl - wf;
          const Real rq = (q > 0.0 && po > 0.0) ? fmin(1.0, (1.0 - 1.0e-12)*q/po) : 0.0;
          if (rq < rc) {
            rn = -rq;
            a += 1.0;
            if ((i < is + zd || i > ie - zd ||
                        (zx2 && (j < js + zd || j > je - zd)) ||
                        (zx3 && (k < ks + zd || k > ke - zd)))) {ab += 1.0;}
          }
        }
      }
      iw_(m,s1,k,j,i) = rn;
    };
    par_reduce_flat2("vet_t2_ffct_r", 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a,
                  Real &ab) {
      rnd(m, k, j, i, a, ab);
    }, Kokkos::Sum<Real>(nv[0]), Kokkos::Sum<Real>(nv[1]));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, nv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const Real nn = nv[0];
    if (!(nn > 0.0)) {break;}
    ncell += nn;
    nround = r + 1;
    std::swap(sa, sb);   // (the round's R is read next)
    if (zfirst || nv[1] > 0.0) {
      rhalo(sa);
      zl_halos += 1.0;
      zfirst = false;
    } else {
      fghost(sa, sb);
    }
  }
  {   // (the last round changed nothing: |R| of its read slot, ghosts too)
    auto r0 = Kokkos::subview(iw, Kokkos::ALL, sa, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    auto r1 = Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    par_for("vet_t2_ffct_abs", DevExeSpace(), 0, nmb1, 0, r0.extent_int(1) - 1, 0,
            r0.extent_int(2) - 1, 0, r0.extent_int(3) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      r1(m,k,j,i) = fabs(r0(m,k,j,i));
    });
  }
  if (nround >= 4096) {
    ZalesakCapWarn("time2_face_fct", pmy_pack->pmesh->ncycle, t2_solve);
  }
  zl_rounds += static_cast<Real>(nround);
  t2_ffct_cells += ncell;
  t2_ffct_rmax = std::max(t2_ffct_rmax, static_cast<Real>(nround));
  if (nround == 0) {return;}   // nothing limited: t2inc untouched (bitwise)
  // apply: the E channel of every cell next to a limited face, then the face F channels
  Real lv[3] = {0.0, 0.0, 0.0}, qmin = 1.0e300;
  // a flat range (faster than a 4-D MDRange reduction)
  const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("vet_t2_ffct_e", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &sl, Real &sw, Real &se, Real &sq) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    Real w, po, pi, pl, dw;
    bool lim;
    cell(m, k, j, i, w, po, pi, pl, dw, lim, ca);
    const Real rc = iw_(m,ca,k,j,i);
    if (lim) {in_(m,VET_T2_E,k,j,i) += dw/cvol(m,k,j,i);}
    if (rc < 1.0) {sw += (1.0 - rc)*po;}
    if (!(hoi && hm_(m,k,j,i) > 0.0) && !(w + dw > 0.0)) {sl += 1.0;}
    if (lim) {   // (the limited old vector against E^n)
      const Real q = (u0_(m,VET_E,k,j,i) + in_(m,VET_T2_E,k,j,i))/un_(m,VET_E,k,j,i);
      if (q < eth) {se += 1.0;}
      sq = fmin(sq, q);
    }
  }, Kokkos::Sum<Real>(lv[0]), Kokkos::Sum<Real>(lv[1]), Kokkos::Sum<Real>(lv[2]),
     Kokkos::Min<Real>(qmin));
  // clean staggered: face F is the radiation momentum, left alone by the limiter, so
  // momentum and L_z stay exact
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, lv, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &qmin, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  t2_ffct_qmin = std::min(t2_ffct_qmin, qmin);
  t2_ffct_nemp += lv[2];
  t2_ffct_left += lv[0];
  t2_ffct_de += lv[1];
}

int RadiationVet::Time2FctSlot(const int q) const {
  return q ? VET_IW_S2 : VET_IW_S1;
}

}  // namespace radvet
