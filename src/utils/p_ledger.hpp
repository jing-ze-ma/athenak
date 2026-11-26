#ifndef UTILS_P_LEDGER_HPP_
#define UTILS_P_LEDGER_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file p_ledger.hpp
//! \brief the momentum ledger, <hydro>/p_ledger = true (read only when named): the
//! change of P = gas + radiation momentum over each accepted step, booked by operator.
//! P: Cartesian (p_x, p_y, p_z) on a Cartesian mesh; the radial p_r on spherical polar
//! and cubed sphere (gas m_r; radiation the x1 face momenta).  Gas: sum V m; radiation:
//! each cell's halves of its faces' V F_f / c^2 (f0x1..3; a face two cells share counts
//! (V_L + V_R)/2).  At the eint_ledger marks P is summed exactly (repro::ReproSum: any
//! layout and rank count) and its change since the previous mark booked to the class
//! of the operator in between:
//!   hflx  the hydro fluxes through the physical boundary faces (outward +, Stage),
//!   ntel  the flux divergence less hflx: what does not telescope (cubed-sphere seams),
//!   src   coordinate, physics, am and user sources, the polar ring average (named
//!         physical sources: gravity, the geometric terms, the sponge),
//!   rad   the radiation update and write-back: the radiation's boundary stress and
//!         geometric terms (physical),
//!   cut   the face-flux cone cut the gas does not take (flux_cone_gas = false),
//!   fwb   the fwb residual's momentum where it is not projected out,
//!   flr   floors: dfloor_flux_limit's momentum, ConToPrim, the floor-gas cap,
//!   mrg   the hydro's merged lateral cells (merge_lat) average,
//!   gap   any other change between marks (0 when the marks cover the step).
//! cut, fwb and the reference force the radiation withholds from the gas (force_reference:
//! the hydro's well-balanced background gives it, booked to src) are the write-back's
//! sums (Withheld) moved out of rad, and the stage starts' shares of them (the slopes'
//! increments, as LzBook's).
//! The first cycle's start (gap_cycle) takes the baseline: the start-up is not booked.
//! The RK average u0 = gam0 u0 + gam1 u1 (Stage) scales the gas books by gam0, so they
//! are the change in the step's final state.  Per accepted step (Commit), rank 0 prints
//! CLEAN PLEDGER with the step's books (x 1/dt) and END = P(t + dt) - P(t) - sum of the
//! books (round-off); the history takes them (Hist).  Reads the state only: a run with
//! the ledger is bitwise the run without it.
//! <hydro>/p_ledger_lz (read only when named; 3-D spherical polar with <rad_vet>): also
//! L_z = gas + radiation (LzTotal less the L_z ledger) in slot 1, the same classes (hflx:
//! am_bflx's torque of the x1 faces; cut, fwb and the reference force stay in rad), and
//! the L_z ledger's growth lzb and rbt, the radiation's torque through the physical x1
//! faces (LzBook, a flux outside the ledger: Lz_rflx), so rad + lzb + rbt is the
//! write-back's L_z neither holds; CLEAN LZLEDGER, user.hst LzL_* in place of Pr_*.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "outputs/outputs.hpp"
#include "parameter_input.hpp"
#include "rad_vet/rad_vet.hpp"
#include "utils/repro_sum.hpp"

namespace p_ledger {
constexpr int kNc = 9, kGap = kNc - 1;   // hflx ntel src rad cut fwb flr mrg gap
inline bool on = false, lz = false;   // (p_ledger_lz)
inline Real lb0[2], lrate[2];             // L_z ledger, its x1 torque at the base; rates
inline Real lz0[2];   // gas, radiation L_z at the last step's end (p_ledger_lz0, restart)
inline bool have_base = false, started = false;
inline Real pg[3], pr[3], pn[3];          // gas, radiation at the last mark; gas of u1
inline Real bg[kNc][3], br[kNc][3];       // this step's books: gas part, radiation part
inline Real p0[3], rate[kNc][3], ptot[3], pend[3], gross = 0.0;   // the last step's
inline Real ws[3][2][3], wc[3][3];   // withheld b (cut, fwb, fref): slopes K1, K2; increment
inline ParameterInput *ppin = nullptr;   // (the run's: the slopes K1 in p_ledger_k1)
inline const char *nm[kNc] = {"hflx", "ntel", "src", "rad", "cut", "fwb", "flr", "mrg",
                              "gap"};

inline int Class(const char *n) {
  const char *src[] = {"psrc", "coord", "am", "ringavg", "usrc_rest", "vet_post"};
  const char *flr[] = {"dfl", "c2p", "vet_c2p", "floorcap"};
  if (!std::strcmp(n, "flux")) return 1;
  for (auto s : src) if (!std::strcmp(n, s)) return 2;
  if (!std::strcmp(n, "vet_ravg")) return 3;
  for (auto s : flr) if (!std::strcmp(n, s)) return 6;
  if (!std::strcmp(n, "mrglat")) return 7;
  return kGap;
}

//! the exact sums of gas momentum (u, components 0..2), radiation momentum (3..5) and
//! the gross sum of their sizes (6) over the active cells and ranks
inline void Sums(Mesh *pm, const DvceArray5D<Real> &u, const bool rad, Real out[7]) {
  auto *pp = pm->pmb_pack;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, js = ic.js, ks = ic.ks;
  const int nx1 = ic.nx1, nji = ic.nx2*nx1, nkji = ic.nx3*nji;
  const bool cur = pm->use_spherical_polar || pm->use_cubed_sphere;
  const int nc = cur ? 1 : 3;
  const bool hr = rad && pp->pradvet != nullptr && pp->pradvet->f0x1.extent_int(0) > 0;
  const bool md = pm->multi_d, td = pm->three_d;
  auto vol = pp->pcoord->volume;   // (curvilinear only)
  auto &sz = pp->pmb->mb_size;
  DvceArray4D<Real> f[3];
  Real ic2 = 0.0;
  if (hr) {
    f[0] = pp->pradvet->f0x1;  f[1] = pp->pradvet->f0x2;  f[2] = pp->pradvet->f0x3;
    ic2 = 1.0/SQR(pp->pradvet->c_light);
  }
  repro::ReproVal<7> v;
  Kokkos::parallel_reduce("p_ledger", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<7> &a) {
    const int m = idx/nkji, k = (idx%nkji)/nji + ks, j = (idx%nji)/nx1 + js;
    const int i = idx%nx1 + is;
    const auto &z = sz.d_view(m);
    const Real w = cur ? vol(m,k,j,i) : z.dx1*z.dx2*z.dx3;
    for (int c = 0; c < nc; ++c) {
      const Real g = w*u(m,IM1+c,k,j,i);
      const bool hf = hr && (c == 0 || (c == 1 && md) || (c == 2 && td));
      const Real r = hf ? 0.5*w*ic2*(f[c](m,k,j,i) + f[c](m,k+(c == 2),j+(c == 1),
                                                          i+(c == 0))) : 0.0;
      a.s[c].Add(g);
      a.s[3+c].Add(r);
      a.s[6].Add(fabs(g) + fabs(r));
    }
  }, repro::ReproRed<7, Kokkos::HostSpace>(v));
  repro::AllMerge(v.s, 7);
  for (int q = 0; q < 7; ++q) out[q] = v.s[q].Value();
  if (lz && pp->pradvet != nullptr) {   // (p_ledger_lz: L_z, LzTotal's parts)
    out[1] = pp->pradvet->LzTotal(1, &u);
    out[4] = rad ? pp->pradvet->LzTotal(2) : 0.0;
  }
}

//! p_ledger_lz: the L_z ledger and its x1-face torque part (LzBook), pending included
inline void LzLedger(Mesh *pm, Real l[2]) {
  auto *v = pm->pmb_pack->pradvet;
  l[0] = (v != nullptr) ? v->cln_lzl + v->cln_lzp : 0.0;
  l[1] = (v != nullptr) ? v->cln_lzr[0] + v->cln_lzr[4] : 0.0;
}

//! p_ledger_lz: the step's gas L_z through the physical x1 faces (am_bflx, RK-weighted)
inline Real AmFlux(Mesh *pm) {
  hydro::AmModule *pam = hydro::AmModuleOf(pm->pmb_pack);
  if (pam == nullptr || pam->am_bflx.extent_int(0) == 0) return 0.0;
  auto bf = pam->am_bflx;
  auto &ic = pm->mb_indcs;
  const int ks = ic.ks, js = ic.js, n2 = ic.nx2, nkj = ic.nx3*n2;
  repro::ReproVal<1> v;
  Kokkos::parallel_reduce("p_ledger_lz", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pm->pmb_pack->nmb_thispack*2*nkj),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<1> &a) {
    const int m = idx/(2*nkj), s = (idx/nkj)%2, k = (idx%nkj)/n2 + ks, j = idx%n2 + js;
    a.s[0].Add(bf(m,2,s,k,j));
  }, repro::ReproRed<1, Kokkos::HostSpace>(v));
  repro::AllMerge(v.s, 1);
  return v.s[0].Value();
}

//! <hydro>/p_ledger on: the withheld slopes K1 the restart file keeps (p_ledger_k1, every
//! accepted step), so the ledger restarts bitwise; a file without them: NOTE, 0 (the
//! first step's stage-start shares of cut, fwb and the reference force stay in rad)
inline void Init(ParameterInput *pin) {
  ppin = pin;
  if (pin->DoesParameterExist("hydro", "p_ledger_lz0")) {   // (p_ledger_lz)
    std::istringstream is(pin->GetString("hydro", "p_ledger_lz0"));
    is >> lz0[0] >> lz0[1];
  }
  if (pin->DoesParameterExist("hydro", "p_ledger_k1")) {
    std::istringstream is(pin->GetString("hydro", "p_ledger_k1"));
    for (int b = 0; b < 3; ++b) for (int q = 0; q < 3; ++q) is >> ws[b][0][q];
  } else if (global_variable::restart_run && global_variable::my_rank == 0) {
    std::printf("### NOTE <hydro>/p_ledger: the restart file holds no p_ledger_k1; the "
                "first step's stage-start shares of cut, fwb and the reference force are "
                "booked as rad\n");
  }
}

//! book the change since the previous mark to class c (book = false: baseline only)
inline void Take(Mesh *pm, const int c, const bool book) {
  Real s[7];
  Sums(pm, pm->pmb_pack->phydro->u0, true, s);
  auto *v = pm->pmb_pack->pradvet;
  if (lz && (v == nullptr || !v->cln_init)) {   // (a restart before LzTotal: the file's)
    s[1] = lz0[0];
    s[4] = lz0[1];
  }
  for (int q = 0; q < 3; ++q) {
    if (book && have_base) {
      bg[c][q] += s[q] - pg[q];
      br[c][q] += s[3+q] - pr[q];
    }
    pg[q] = s[q];
    pr[q] = s[3+q];
  }
  gross = s[6];
  for (int q = 0; q < 3 && !have_base; ++q) p0[q] = pg[q] + pr[q];
  if (lz && !have_base) LzLedger(pm, lb0);
  have_base = true;
}

inline void Discard(Mesh *pm);
inline void Mark(Mesh *pm, const char *name) {
  if (!on) return;
  if (!started && !std::strcmp(name, "gap_cycle")) {   // (the start-up not booked)
    started = true;
    Discard(pm);
    return;
  }
  // (no sum where the next mark always books the same class: the source marks before
  // usrc_rest, and gap_pre, which Stage's own sum follows)
  for (auto s : {"gap_pre", "psrc", "coord", "am", "ringavg"}) {
    if (!std::strcmp(name, s)) return;
  }
  Take(pm, Class(name), true);
}

//! drop this step's books, new baseline (dt_retry / backward-Euler redo)
inline void Discard(Mesh *pm) {
  if (!on) return;
  std::memset(bg, 0, sizeof(bg));
  std::memset(br, 0, sizeof(br));
  have_base = false;
  Take(pm, kGap, false);
}

//! a write-back's sum a of the gas momentum the radiation step withholds (b: 0 the cone
//! cut, 1 the fwb residual, 2 the reference force; P components), out of rad to cut, fwb,
//! src; k >= 0 also the stage's slope K_{k+1} = a/h (probe: the slope only)
inline void Withheld(const int b, const Real a[3], const int k, const Real h,
                     const bool book) {
  if (!on) return;
  for (int q = 0; q < (lz ? 1 : 3); ++q) {   // (p_ledger_lz: slot 1 is L_z)
    if (k >= 0) ws[b][k][q] = a[q]/h;
    if (!book) continue;
    bg[3][q] += a[q];
    bg[(b == 2) ? 2 : 4 + b][q] -= a[q];
  }
}

//! the momentum a that dfloor_flux_limit's held cells send out by physical faces
//! (DflMomentum): out of flr, a boundary flux (hflx, outward +)
inline void Outflow(const Real a[3]) {
  if (!on) return;
  for (int q = 0; q < (lz ? 1 : 3); ++q) {
    bg[6][q] += a[q];
    bg[0][q] += a[q];
  }
}

//! the stage start's increment c1 K1 + c2 K2 of the withheld slopes (Time2FormStage),
//! k1 = 1: K1 <- K2 (an accepted probe); Start books it as the stage start adds it
inline void Inc(const Real c1, const Real c2, const int k1 = 0) {
  for (int b = 0; b < 3; ++b) {
    for (int q = 0; q < 3; ++q) {
      if (k1) ws[b][0][q] = ws[b][1][q];
      wc[b][q] = c1*ws[b][0][q] + c2*ws[b][1][q];
    }
  }
}
inline void Start() {
  for (int b = 0; b < 3; ++b) Withheld(b, wc[b], -1, 1.0, true);
}

//! before the RK update of a stage: u0 -> gam0 u0 + gam1 u1 scales the gas books; the
//! boundary fluxes beta dt sum (+-) A F of the physical faces go to hflx, -them to ntel
//! (whose flux-divergence mark follows)
inline void Stage(MeshBlockPack *pp, const Real gam0, const Real gam1, const Real bdt,
                  const int stage, const DvceFaceFld5D<Real> &uflx,
                  const DvceArray5D<Real> &u1) {
  if (!on) return;
  Mesh *pm = pp->pmesh;
  Take(pm, kGap, true);
  if (stage == 1) {
    Real s[7];
    Sums(pm, u1, false, s);
    for (int q = 0; q < 3; ++q) pn[q] = s[q];
  }
  for (int q = 0; q < 3; ++q) {
    for (int c = 0; c < kNc; ++c) bg[c][q] *= gam0;
    pg[q] = gam0*pg[q] + gam1*pn[q];
  }
  auto &ic = pm->mb_indcs;
  const int is = ic.is, ie = ic.ie, js = ic.js, je = ic.je, ks = ic.ks, ke = ic.ke;
  const int nmb = pp->nmb_thispack, n1 = ic.nx1, n2 = ic.nx2, n3 = ic.nx3;
  const bool cur = pm->use_spherical_polar || pm->use_cubed_sphere;
  const int nc = cur ? 1 : 3, nd = pm->three_d ? 3 : (pm->multi_d ? 2 : 1);
  auto bcs = pp->pmb->mb_bcs;
  auto &sz = pp->pmb->mb_size;
  auto a1 = pp->pcoord->area.x1f, a2 = pp->pcoord->area.x2f, a3 = pp->pcoord->area.x3f;
  auto f1 = uflx.x1f, f2 = uflx.x2f, f3 = uflx.x3f;
  repro::ReproSum b[3];
  for (int q = 0; q < 3; ++q) b[q].Init();
  for (int d = 0; d < nd; ++d) {
    repro::ReproVal<3> v;   // faces of direction d: (m, side, the other two)
    const int na = (d == 0) ? n2 : n1, nb = (d == 2) ? n2 : n3;
    Kokkos::parallel_reduce("p_ledger_bf", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            nmb*2*na*nb),
    KOKKOS_LAMBDA(const int idx, repro::ReproVal<3> &a) {
      const int m = idx/(2*na*nb), s = (idx/(na*nb))%2, p = idx%na, q = (idx/na)%nb;
      const BoundaryFace bf = static_cast<BoundaryFace>(2*d + s);
      const BoundaryFlag fl = bcs.d_view(m, bf);
      if (fl == BoundaryFlag::block || fl == BoundaryFlag::periodic ||
          fl == BoundaryFlag::panel || fl == BoundaryFlag::polar) return;
      const int i = (d == 0) ? (s ? ie + 1 : is) : is + p;
      const int j = (d == 1) ? (s ? je + 1 : js) : ((d == 0) ? js + p : js + q);
      const int k = (d == 2) ? (s ? ke + 1 : ks) : ks + q;
      const auto &z = sz.d_view(m);
      const Real ar = cur ? ((d == 0) ? a1(m,k,j,i) : ((d == 1) ? a2(m,k,j,i)
                                                                : a3(m,k,j,i)))
                    : ((d == 0) ? z.dx2*z.dx3 : ((d == 1) ? z.dx1*z.dx3 : z.dx1*z.dx2));
      for (int c = 0; c < nc; ++c) {
        const Real fx = (d == 0) ? f1(m,IM1+c,k,j,i) : ((d == 1) ? f2(m,IM1+c,k,j,i)
                                                                 : f3(m,IM1+c,k,j,i));
        a.s[c].Add((s ? 1.0 : -1.0)*bdt*ar*fx);
      }
    }, repro::ReproRed<3, Kokkos::HostSpace>(v));
    for (int q = 0; q < 3; ++q) b[q].Join(v.s[q]);
  }
  repro::AllMerge(b, 3);
  for (int q = 0; q < nc; ++q) {
    bg[0][q] += b[q].Value();
    bg[1][q] += b[q].Value();
  }
}

//! an accepted step: the books x 1/dt to rate, END, the CLEAN PLEDGER line
inline void Commit(Mesh *pm) {
  if (!on) return;
  Take(pm, kGap, true);
  const Real dt = pm->dt;
  Real lb[2];
  if (lz) {   // (p_ledger_lz: the gas's x1 torque, out of ntel; the ledger's growth)
    const Real hz = AmFlux(pm);
    bg[0][1] += hz;
    bg[1][1] += hz;
    LzLedger(pm, lb);
    for (int q = 0; q < 2; ++q) {
      lrate[q] = (lb[q] - lb0[q])/dt;
      lb0[q] = lb[q];
    }
  }
  const int nc = (pm->use_spherical_polar || pm->use_cubed_sphere) ? 1 : 3;
  for (int q = 0; q < 3; ++q) {   // dP = -hflx + the other books
    ptot[q] = pg[q] + pr[q];
    pend[q] = ptot[q] - p0[q] + bg[0][q];
    rate[0][q] = bg[0][q]/dt;
    for (int c = 1; c < kNc; ++c) {
      pend[q] -= bg[c][q] + br[c][q];
      rate[c][q] = (bg[c][q] + br[c][q])/dt;
    }
    p0[q] = ptot[q];
  }
  std::memset(bg, 0, sizeof(bg));
  std::memset(br, 0, sizeof(br));
  if (ppin != nullptr) {
    char t[64];
    std::string k;
    for (int b = 0; b < 3; ++b) {
      for (int q = 0; q < 3; ++q) {
        std::snprintf(t, sizeof(t), "%s%.17e", k.empty() ? "" : " ", ws[b][0][q]);
        k += t;
      }
    }
    ppin->SetString("hydro", "p_ledger_k1", k);
    if (lz) {
      std::snprintf(t, sizeof(t), "%.17e %.17e", pg[1], pr[1]);
      ppin->SetString("hydro", "p_ledger_lz0", t);
    }
  }
  if (global_variable::my_rank != 0) return;
  std::printf("CLEAN PLEDGER cycle %d dt %.6e", pm->ncycle, dt);
  for (int c = 0; c < kNc; ++c) {
    std::printf(" %s", nm[c]);
    for (int q = 0; q < nc; ++q) std::printf(" %.6e", rate[c][q]);
  }
  std::printf(" P");
  for (int q = 0; q < nc; ++q) std::printf(" %.16e", ptot[q]);
  std::printf(" END");
  for (int q = 0; q < nc; ++q) std::printf(" %.6e", pend[q]);
  std::printf(" gross %.6e\n", gross);
  if (!lz) return;
  std::printf("CLEAN LZLEDGER cycle %d dt %.6e", pm->ncycle, dt);
  for (int c = 0; c < kNc; ++c) std::printf(" %s %.6e", nm[c], rate[c][1]);
  std::printf(" lzb %.6e rbt %.6e Lz %.16e END %.6e\n", lrate[0], lrate[1], ptot[1],
              pend[1]);
}

//! the history columns P*_tot, P*_hflx, ..., P*_gap (rank 0 holds the global values)
inline void Hist(HistoryData *pd, Mesh *pm) {
  if (!on) return;
  const bool cur = pm->use_spherical_polar || pm->use_cubed_sphere;
  const char *cn[3] = {cur ? "Pr" : "Px", lz ? "LzL" : "Py", "Pz"};
  const bool r0 = (global_variable::my_rank == 0);
  for (int q = lz; q < (lz ? 2 : (cur ? 1 : 3)); ++q) {   // (p_ledger_lz: L_z alone)
    for (int c = -1; c < kNc + 2*lz; ++c) {
      const int n = pd->nhist++;
      pd->label[n] = std::string(cn[q]) + "_" + ((c < 0) ? "tot" : ((c < kNc) ? nm[c] :
                     ((c == kNc) ? "lzb" : "rbt")));
      pd->hdata[n] = r0 ? ((c < 0) ? ptot[q] : ((c < kNc) ? rate[c][q] :
                                                  lrate[c - kNc])) : 0.0;
    }
  }
}

}  // namespace p_ledger

#endif  // UTILS_P_LEDGER_HPP_
