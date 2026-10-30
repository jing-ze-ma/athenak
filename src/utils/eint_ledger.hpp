#ifndef UTILS_EINT_LEDGER_HPP_
#define UTILS_EINT_LEDGER_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eint_ledger.hpp
//! \brief floorheat-1029 <problem>/eint_ledger (DIAGNOSTIC, opt-in, default off): a
//! per-operator INTERNAL-ENERGY ledger of the low-density (floor) gas.
//!
//! Mark(pm, "name") is called at the task boundaries (hydro flux update, coordinate
//! sources, am_transport, polar ring averaging, the pgen's source kernel / density floor /
//! reference-force work, ConToPrim, the M1 solve, the M1 ring averaging, the M1 fluid
//! refresh).  At each mark the hydro u0 of every active cell is turned into
//!   rho, KE = |m|^2/(2 rho), e = E - KE - rho Phi (etotgrav)
//! and compared with the previous mark.  The difference is booked to the operator that
//! ran in between ("name"):
//!   - rank-reduced sums over four cell classes (by the state at the previous mark):
//!     floor gas (rho < eint_ledger_rho) off / near the poles (theta < eint_ledger_theta
//!     from either pole), and the HOT floor gas (e/rho > eint_ledger_hot) off / near
//!     the poles; per class: d(e) vol, d(KE) vol, d(rho) vol, H vol with
//!     H = d(e) - (e/rho)_prev d(rho) (heating at fixed specific energy), cell count;
//!   - per cell, the specific heating d(e/rho), accumulated over the run.
//! A cycle's entries are committed only when the cycle is accepted (a dt_retry or a
//! hesdirk2 backward-Euler redo discards the attempt and takes a new baseline).
//! Output (rank 0): per cycle, one ELEDGER line per operator with the cycle's class sums,
//! and one ELHOT line: the hottest floor cell (max e/rho, r > eint_ledger_rmin) and its
//! per-operator d(e/rho) in this cycle.  Every eint_ledger_dump cycles each rank writes
//! eledger.<cycle>.r<rank>.bin (run-accumulated per-cell d(e/rho) per operator + the
//! state) into the run directory.
//! Reads u0 only: a run with the ledger on is bitwise the run without it.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace eint_ledger {

constexpr int kNcls = 5;   // + class 4: the region (polarpert-1029)
constexpr int kNq = 5;
constexpr int kNs = kNcls*kNq;
constexpr int kNopMax = 24;

struct Sums {
  Real v[kNs];
  KOKKOS_INLINE_FUNCTION Sums() {for (int n = 0; n < kNs; ++n) v[n] = 0.0;}
  KOKKOS_INLINE_FUNCTION Sums &operator+=(const Sums &o) {
    for (int n = 0; n < kNs; ++n) v[n] += o.v[n];
    return *this;
  }
};

}  // namespace eint_ledger

namespace Kokkos {
template <>
struct reduction_identity<eint_ledger::Sums> {
  KOKKOS_FORCEINLINE_FUNCTION static eint_ledger::Sums sum() {return eint_ledger::Sums();}
};
}  // namespace Kokkos

namespace eint_ledger {

inline bool on = false;
inline Real rho_cut = 1.0e-7;
inline Real hot_cut = 1.0e4;
inline Real theta_cut = 0.2618;   // rad (15 deg)
inline Real rmin = 0.0;
inline int dump_every = 0;
// polarpert-1029: region class 4 = r in [reg_rmin, reg_rmax], theta within reg_th of a pole
// (any density); q = d(rho v_r^2/2) vol, d(rho v_th^2/2) vol, d(e) vol, d(rho v_ph^2/2) vol, n
inline Real reg_rmin = 0.0, reg_rmax = -1.0, reg_th = 0.0;
inline std::string dump_dir = ".";

inline bool have_base = false;
inline std::vector<std::string> names;
inline std::vector<std::vector<Real>> cyc_sums;   // [op][kNs] this cycle (rank-local)
inline std::vector<std::vector<Real>> tot_sums;   // [op][kNs] accepted cycles
inline DvceArray4D<Real> d_prv, e_prv, k_prv;
inline DvceArray5D<Real> cyc, tot;                // (op, m, k, j, i) active cells
inline Kokkos::View<float*****, DevMemSpace> cyv;  // (3 op + q, m, k, j, i): cycle d(v_q)
inline DvceArray4D<Real> v1_prv, v2_prv, v3_prv;

inline int OpIndex(const char *name) {
  for (int n = 0; n < static_cast<int>(names.size()); ++n) {
    if (names[n] == name) return n;
  }
  if (static_cast<int>(names.size()) >= kNopMax) return kNopMax - 1;
  names.push_back(name);
  cyc_sums.push_back(std::vector<Real>(kNs, 0.0));
  tot_sums.push_back(std::vector<Real>(kNs, 0.0));
  return static_cast<int>(names.size()) - 1;
}

inline void Alloc(Mesh *pm) {
  auto &ic = pm->mb_indcs;
  const int nmb = pm->pmb_pack->nmb_thispack;
  if (cyc.extent_int(1) == nmb) return;
  Kokkos::realloc(d_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(e_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(k_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(cyc, kNopMax, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(tot, kNopMax, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(cyv, 3*kNopMax, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(v1_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(v2_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::realloc(v3_prv, nmb, ic.nx3, ic.nx2, ic.nx1);
  Kokkos::deep_copy(cyv, 0.0f);
  Kokkos::deep_copy(cyc, 0.0);
  Kokkos::deep_copy(tot, 0.0);
  have_base = false;
  if (global_variable::my_rank == 0) {
    std::printf("ELEDGER on: rho_cut %.3e hot_cut %.3e theta_cut %.4f rmin %.4f dump %d"
                " (classes: flo_np flo_p hot_np hot_p; q: dE dKE dM H n)\n", rho_cut,
                hot_cut, theta_cut, rmin, dump_every);
  }
}

//! the core: book the change since the previous mark to op (book = false: baseline only)
inline void Take(Mesh *pm, const int op, const bool book) {
  Alloc(pm);
  auto *pmbp = pm->pmb_pack;
  auto *ph = pmbp->phydro;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, js = ic.js, ks = ic.ks;
  const int ni = ic.nx1, nj = ic.nx2, nk = ic.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = ph->u0;
  const bool grav = ph->use_etotgrav;
  auto phi = ph->phicc0;
  auto vol = pmbp->pcoord->volume;
  auto &x2v = pmbp->pcoord->x2v;
  auto dp = d_prv, ep = e_prv, kp = k_prv;
  auto cy = cyc;
  auto cv = cyv;
  auto a1 = v1_prv, a2 = v2_prv, a3 = v3_prv;
  const Real rc = rho_cut, hc = hot_cut, tc = theta_cut;
  const Real rr0 = reg_rmin, rr1 = reg_rmax, rth = reg_th;
  auto &rx1 = pmbp->pcoord->x1v;
  const bool bk = book && have_base;
  const int nji = nj*ni, nkji = nk*nji;
  Sums s;
  Kokkos::parallel_reduce("eint_ledger", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Sums &sum) {
    const int m = idx/nkji;
    const int kk = (idx - m*nkji)/nji;
    const int jj = (idx - m*nkji - kk*nji)/ni;
    const int ii = idx%ni;
    const int k = kk + ks, j = jj + js, i = ii + is;
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
    const Real ke = (d > 0.0) ? 0.5*(m1*m1 + m2*m2 + m3*m3)/d : 0.0;
    const Real ei = u0(m,IEN,k,j,i) - ke - (grav ? d*phi(m,k,j,i) : 0.0);
    if (bk) {
      const Real d0 = dp(m,kk,jj,ii), e0 = ep(m,kk,jj,ii), k0 = kp(m,kk,jj,ii);
      if (d0 < rc) {
        const Real th = x2v(m,j);
        const bool pol = (th < tc) || (th > M_PI - tc);
        const Real eps0 = (d0 > 0.0) ? e0/d0 : 0.0;
        const Real v = vol(m,k,j,i);
        const Real q[kNq] = {(ei - e0)*v, (ke - k0)*v, (d - d0)*v,
                             ((ei - e0) - eps0*(d - d0))*v, 1.0};
        const int c0 = pol ? 1 : 0;
        for (int n = 0; n < kNq; ++n) sum.v[c0*kNq + n] += q[n];
        if (eps0 > hc) {
          for (int n = 0; n < kNq; ++n) sum.v[(c0 + 2)*kNq + n] += q[n];
        }
      }
      {
        const Real r = rx1(m,i), th = x2v(m,j);
        if (r >= rr0 && r <= rr1 && (th < rth || th > M_PI - rth) && d > 0.0 && d0 > 0.0) {
          const Real v = vol(m,k,j,i);
          const Real b1 = a1(m,kk,jj,ii), b2 = a2(m,kk,jj,ii), b3 = a3(m,kk,jj,ii);
          sum.v[4*kNq + 0] += 0.5*(m1*m1/d - d0*b1*b1)*v;
          sum.v[4*kNq + 1] += 0.5*(m2*m2/d - d0*b2*b2)*v;
          sum.v[4*kNq + 2] += (ei - e0)*v;
          sum.v[4*kNq + 3] += 0.5*(m3*m3/d - d0*b3*b3)*v;
          sum.v[4*kNq + 4] += 1.0;
        }
      }
      if (d > 0.0 && d0 > 0.0) {
        cy(op,m,kk,jj,ii) += ei/d - e0/d0;
        cv(3*op,m,kk,jj,ii) += static_cast<float>(m1/d - a1(m,kk,jj,ii));
        cv(3*op+1,m,kk,jj,ii) += static_cast<float>(m2/d - a2(m,kk,jj,ii));
        cv(3*op+2,m,kk,jj,ii) += static_cast<float>(m3/d - a3(m,kk,jj,ii));
      }
    }
    dp(m,kk,jj,ii) = d;
    ep(m,kk,jj,ii) = ei;
    kp(m,kk,jj,ii) = ke;
    a1(m,kk,jj,ii) = (d > 0.0) ? m1/d : 0.0;
    a2(m,kk,jj,ii) = (d > 0.0) ? m2/d : 0.0;
    a3(m,kk,jj,ii) = (d > 0.0) ? m3/d : 0.0;
  }, Kokkos::Sum<Sums>(s));
  if (bk) {
    for (int n = 0; n < kNs; ++n) cyc_sums[op][n] += s.v[n];
  }
  have_base = true;
}

inline void Mark(Mesh *pm, const char *name) {
  if (!on) return;
  Take(pm, OpIndex(name), true);
}

//! drop this cycle's entries and take a new baseline (dt_retry / backward-Euler redo)
inline void Discard(Mesh *pm) {
  if (!on) return;
  Alloc(pm);
  Kokkos::deep_copy(cyc, 0.0);
  Kokkos::deep_copy(cyv, 0.0f);
  for (auto &v : cyc_sums) for (auto &x : v) x = 0.0;
  Take(pm, 0, false);
}

inline void WriteDump(Mesh *pm) {
  auto *pmbp = pm->pmb_pack;
  auto &ic = pm->mb_indcs;
  const int nmb = pmbp->nmb_thispack;
  const int ni = ic.nx1, nj = ic.nx2, nk = ic.nx3;
  const int nop = static_cast<int>(names.size());
  char fn[512];
  std::snprintf(fn, sizeof(fn), "%s/eledger.%06d.r%04d.bin", dump_dir.c_str(),
                pm->ncycle, global_variable::my_rank);
  FILE *f = std::fopen(fn, "wb");
  if (f == nullptr) {std::printf("ELEDGER cannot open %s\n", fn); return;}
  int hdr[6] = {nop, nmb, nk, nj, ni, pm->ncycle};
  std::fwrite(hdr, sizeof(int), 6, f);
  for (int n = 0; n < nop; ++n) {
    char nm[32] = {0};
    std::snprintf(nm, sizeof(nm), "%s", names[n].c_str());
    std::fwrite(nm, 1, 32, f);
  }
  auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
  auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
  auto ht = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tot);
  auto hd = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d_prv);
  auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(), e_prv);
  std::vector<float> buf;
  for (int m = 0; m < nmb; ++m) {
    int gid = pmbp->gids + m;
    std::fwrite(&gid, sizeof(int), 1, f);
    std::vector<double> x;
    for (int i = 0; i < ni; ++i) x.push_back(hx1(m, i + ic.is));
    for (int j = 0; j < nj; ++j) x.push_back(hx2(m, j + ic.js));
    for (int k = 0; k < nk; ++k) x.push_back(hx3(m, k + ic.ks));
    std::fwrite(x.data(), sizeof(double), x.size(), f);
    buf.resize(static_cast<size_t>(nk)*nj*ni);
    for (int n = 0; n < nop + 2; ++n) {
      size_t p = 0;
      for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i) {
        Real v;
        if (n < nop) {
          v = ht(n,m,k,j,i);
        } else if (n == nop) {
          v = hd(m,k,j,i);
        } else {
          v = (hd(m,k,j,i) > 0.0) ? he(m,k,j,i)/hd(m,k,j,i) : 0.0;
        }
        buf[p++] = static_cast<float>(v);
      }
      std::fwrite(buf.data(), sizeof(float), buf.size(), f);
    }
  }
  std::fclose(f);
}

//! the cycle is accepted: totals, the ELEDGER / ELHOT lines, the dump
inline void Commit(Mesh *pm) {
  if (!on || !have_base) return;
  const int nop = static_cast<int>(names.size());
  std::vector<Real> loc(static_cast<size_t>(nop)*kNs), glb(loc.size());
  for (int n = 0; n < nop; ++n) for (int q = 0; q < kNs; ++q) loc[n*kNs + q] = cyc_sums[n][q];
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(loc.data(), glb.data(), static_cast<int>(loc.size()), MPI_ATHENA_REAL,
                MPI_SUM, MPI_COMM_WORLD);
#else
  glb = loc;
#endif
  if (global_variable::my_rank == 0) {
    for (int n = 0; n < nop; ++n) {
      std::printf("ELEDGER c=%d t=%.7e op=%-11s", pm->ncycle, pm->time, names[n].c_str());
      for (int c = 0; c < kNcls; ++c) {
        std::printf(" |");
        for (int q = 0; q < kNq - 1; ++q) std::printf(" %+.3e", glb[n*kNs + c*kNq + q]);
        std::printf(" %.0f", glb[n*kNs + c*kNq + kNq - 1]);
      }
      std::printf("\n");
    }
  }
  auto *pmbp = pm->pmb_pack;
  if (reg_rmax > reg_rmin) {
    // polarpert-1029: the region's state at the last mark: sum rho v_r^2/2 vol,
    // rho v_th^2/2 vol, rho vol (mass), rho v_r vol (mean flow), and the count
    auto &icr = pm->mb_indcs;
    const int ni = icr.nx1, nj = icr.nx2, nk = icr.nx3, nmb = pmbp->nmb_thispack;
    const int nji = nj*ni, nkji = nk*nji;
    auto dp = d_prv; auto b1 = v1_prv, b2 = v2_prv;
    auto &x1 = pmbp->pcoord->x1v; auto &x2 = pmbp->pcoord->x2v;
    auto vol = pmbp->pcoord->volume;
    const int is = icr.is, js = icr.js, ks = icr.ks;
    const Real rr0 = reg_rmin, rr1 = reg_rmax, rth = reg_th;
    Sums st;
    Kokkos::parallel_reduce("eint_ledger_reg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Sums &sm) {
      const int m = idx/nkji;
      const int kk = (idx - m*nkji)/nji;
      const int jj = (idx - m*nkji - kk*nji)/ni;
      const int ii = idx%ni;
      const Real r = x1(m, ii + is), th = x2(m, jj + js);
      if (!(r >= rr0 && r <= rr1 && (th < rth || th > M_PI - rth))) return;
      const Real d = dp(m,kk,jj,ii), v = vol(m, kk + ks, jj + js, ii + is);
      const Real u1 = b1(m,kk,jj,ii), u2 = b2(m,kk,jj,ii);
      sm.v[0] += 0.5*d*u1*u1*v; sm.v[1] += 0.5*d*u2*u2*v; sm.v[2] += d*v;
      sm.v[3] += d*u1*v; sm.v[4] += 1.0;
    }, Kokkos::Sum<Sums>(st));
    double g[5];
    double l5[5] = {st.v[0], st.v[1], st.v[2], st.v[3], st.v[4]};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(l5, g, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (int q = 0; q < 5; ++q) g[q] = l5[q];
#endif
    if (global_variable::my_rank == 0) {
      std::printf("ELREG c=%d t=%.7e KEr %.6e KEth %.6e M %.6e <v_r> %+.4e vrms_r %.4e "
                  "vrms_th %.4e n %.0f\n", pm->ncycle, pm->time, g[0], g[1], g[2],
                  g[3]/g[2], sqrt(2.0*g[0]/g[2]), sqrt(2.0*g[1]/g[2]), g[4]);
    }
  }
  // the hottest floor cell of this rank (state at the last mark), then the global max
  auto &ic = pm->mb_indcs;
  const int ni = ic.nx1, nj = ic.nx2, nk = ic.nx3, nmb = pmbp->nmb_thispack;
  const int nji = nj*ni, nkji = nk*nji;
  auto dp = d_prv, ep = e_prv;
  auto &x1v = pmbp->pcoord->x1v;
  const Real rc = rho_cut, rm = rmin;
  const int is = ic.is;
  Kokkos::MaxLoc<Real, int>::value_type ml;
  Kokkos::parallel_reduce("eint_ledger_hot", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Kokkos::MaxLoc<Real, int>::value_type &v) {
    const int m = idx/nkji;
    const int ii = idx%ni;
    const int kk = (idx - m*nkji)/nji;
    const int jj = (idx - m*nkji - kk*nji)/ni;
    const Real d = dp(m,kk,jj,ii);
    if (!(d < rc) || d <= 0.0 || x1v(m, ii + is) < rm) return;
    const Real eps = ep(m,kk,jj,ii)/d;
    if (eps > v.val) {v.val = eps; v.loc = idx;}
  }, Kokkos::MaxLoc<Real, int>(ml));
  struct {double v; int r;} in{static_cast<double>(ml.val), global_variable::my_rank}, out;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
  out = in;
#endif
  if (out.r == global_variable::my_rank && ml.loc >= 0 && ml.val > 0.0) {
    const int idx = ml.loc;
    const int m = idx/nkji;
    const int ii = idx%ni;
    const int kk = (idx - m*nkji)/nji;
    const int jj = (idx - m*nkji - kk*nji)/ni;
    std::vector<Real> hc(nop, 0.0);
    for (int n = 0; n < nop; ++n) {
      Kokkos::deep_copy(hc[n], Kokkos::subview(cyc, n, m, kk, jj, ii));
    }
    auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
    auto hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
    auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
    Real hd;
    Kokkos::deep_copy(hd, Kokkos::subview(dp, m, kk, jj, ii));
    std::printf("ELHOT c=%d gid %d (k,j,i) (%d,%d,%d) r %.4f th %.3f deg ph %.4f rho %.3e "
                "e/rho %.4e :", pm->ncycle, pmbp->gids + m, kk, jj, ii, hx1(m, ii + is),
                hx2(m, jj + ic.js)*180.0/M_PI, hx3(m, kk + ic.ks), hd,
                static_cast<double>(ml.val));
    for (int n = 0; n < nop; ++n) {
      if (hc[n] != 0.0) std::printf(" %s %+.3e", names[n].c_str(), hc[n]);
    }
    std::printf("\n");
    for (int dk = -4; dk <= 4; ++dk) {
      const int kq = kk + dk;
      if (kq < 0 || kq >= nk) continue;
      Real dq;
      Kokkos::deep_copy(dq, Kokkos::subview(dp, m, kq, jj, ii));
      std::printf("ELHOTN c=%d dk %+d rho %.3e :", pm->ncycle, dk, dq);
      for (int n = 0; n < nop; ++n) {
        Real de;
        float w1, w2, w3;
        Kokkos::deep_copy(de, Kokkos::subview(cyc, n, m, kq, jj, ii));
        Kokkos::deep_copy(w1, Kokkos::subview(cyv, 3*n, m, kq, jj, ii));
        Kokkos::deep_copy(w2, Kokkos::subview(cyv, 3*n+1, m, kq, jj, ii));
        Kokkos::deep_copy(w3, Kokkos::subview(cyv, 3*n+2, m, kq, jj, ii));
        if (de != 0.0 || w1 != 0.0f || w2 != 0.0f || w3 != 0.0f) {
          std::printf(" %s [%+.2e %+.1f %+.1f %+.1f]", names[n].c_str(), de, w1, w2, w3);
        }
      }
      std::printf("\n");
    }
    std::fflush(stdout);
  }
  // accept: cycle -> totals
  auto cy = cyc, tt = tot;
  Kokkos::parallel_for("eint_ledger_acc", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                       nop*nmb*nkji),
  KOKKOS_LAMBDA(const int idx) {
    const int n = idx/(nmb*nkji);
    const int r = idx - n*nmb*nkji;
    const int m = r/nkji;
    const int kk = (r - m*nkji)/nji;
    const int jj = (r - m*nkji - kk*nji)/ni;
    const int ii = r%ni;
    tt(n,m,kk,jj,ii) += cy(n,m,kk,jj,ii);
    cy(n,m,kk,jj,ii) = 0.0;
  });
  Kokkos::deep_copy(cyv, 0.0f);
  for (int n = 0; n < nop; ++n) {
    for (int q = 0; q < kNs; ++q) {tot_sums[n][q] += cyc_sums[n][q]; cyc_sums[n][q] = 0.0;}
  }
  if (dump_every > 0 && (pm->ncycle % dump_every) == 0) {
    Kokkos::fence();
    WriteDump(pm);
  }
  std::fflush(stdout);
}

//! free the device views before Kokkos::finalize (pgen final hook)
inline void Free() {
  d_prv = DvceArray4D<Real>(); e_prv = DvceArray4D<Real>(); k_prv = DvceArray4D<Real>();
  v1_prv = DvceArray4D<Real>(); v2_prv = DvceArray4D<Real>(); v3_prv = DvceArray4D<Real>();
  cyc = DvceArray5D<Real>(); tot = DvceArray5D<Real>();
  cyv = Kokkos::View<float*****, DevMemSpace>();
}

}  // namespace eint_ledger

#endif  // UTILS_EINT_LEDGER_HPP_
