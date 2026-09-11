#ifndef UTILS_RUNAWAY_SCAN_HPP_
#define UTILS_RUNAWAY_SCAN_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file runaway_scan.hpp
//! \brief <problem>/runaway_scan: an OPERATOR LEDGER for the photospheric runaway.
//!
//! WHAT IT ANSWERS.  In the red-giant open-top runs a photospheric cell reaches
//! 1e16-1e17 K during the t ~ 5.7e5 burst.  Radial conduction is implicit and the
//! angular faces are capped, so conduction cannot lift a cell above its neighbours; the
//! remaining suspects are the hydro flux divergence (a vacuum Riemann state into an
//! evacuated cell), the explicit source terms (gravity / well-balanced), the two-stream
//! apply (its semi-implicit relaxation and the Newton rescue) and ConToPrim's floors.
//! This scan is called at each of those task boundaries and reports, per cycle,
//!   (a) the global max of T_cell/T_neighbour and which cell holds it, and
//!   (b) the internal-energy INCREMENT that boundary gave the watched cell,
//! so the printed sequence is a ledger: which operator added how much, in order.
//!
//! HOW IT READS THE STATE.  Mid-stage w0 is stale, so the internal energy is taken from
//! u0 with exactly the extraction ConToPrim uses -- total energy minus the gravitational
//! term (etotgrav) minus the kinetic energy formed with the cubed-sphere metric (see
//! Coordinates::GnomonicEquiangleRaiseVel) -- and T = eos.Temperature of that.  The scan
//! only READS u0; it writes nothing but its own scratch, so a run with it on is
//! bit-identical to one with it off.
//!
//! COST.  One EOS temperature inversion per cell per boundary, restricted to
//! r > runaway_scan::rmin.  Diagnostic only: default off, and nothing is allocated or
//! launched until it is switched on.

#include <math.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "units/units.hpp"
#include "hydro/hydro.hpp"

namespace runaway_scan {

inline bool on = false;             // <problem>/runaway_scan
inline Real rmin = 3.3e12;          // only cells outside this radius are scanned
inline Real ratio_print = 3.0;      // print a ledger line when max T/T_nb exceeds this
inline Real ratio_state = 10.0;     // dump the 7-cell state when a cell first exceeds
inline int lines = 0;
inline const int maxlines = 2000;
inline const int maxcells = 50;

// device scratch, allocated once on the first call
inline DvceArray4D<Real> *t_cur = nullptr;   // T of every cell, this boundary
inline DvceArray4D<Real> *e_cur = nullptr;   // internal energy density, this boundary
inline DvceArray4D<Real> *e_prv = nullptr;   // ...and at the previous boundary
inline DvceArray1D<Real> *rec = nullptr;     // the readback record
// the watched cell: latched at the end of a cycle, followed through the next one so the
// per-boundary increments form a ledger for ONE cell rather than for whichever cell
// happens to be worst at each boundary
inline int w_m = -1, w_k = -1, w_j = -1, w_i = -1;
inline int w_cycle = -1;
// cells whose 7-cell state has already been dumped (gid<<24 | k<<16 | j<<8 ... packed)
inline std::set<int64_t> dumped;

//----------------------------------------------------------------------------------------
//! \fn void runaway_scan::Scan
//! \brief one ledger entry: T/T_nb over the domain and the watched cell's increment.

inline void Scan(Mesh *pm, const char *opname) {
  if (!on || lines >= maxlines) return;
  // one line, once, so a silent run means "no cell was ever that far above its
  // neighbours" and not "the switch never reached the code"
  {
    static bool announced = false;
    if (!announced && global_variable::my_rank == 0) {
      announced = true;
      std::cout << "### runaway_scan ON: rmin = " << rmin << ", ledger line above T/T_nb"
                << " = " << ratio_print << ", 7-cell dump above " << ratio_state
                << std::endl;
    }
  }
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const int ncells1 = indcs.nx1 + 2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int ncells3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  if (t_cur == nullptr) {
    t_cur = new DvceArray4D<Real>("rsc_t", nmb1+1, ncells3, ncells2, ncells1);
    e_cur = new DvceArray4D<Real>("rsc_e", nmb1+1, ncells3, ncells2, ncells1);
    e_prv = new DvceArray4D<Real>("rsc_ep", nmb1+1, ncells3, ncells2, ncells1);
    rec = new DvceArray1D<Real>("rsc_rec", 32);
    Kokkos::deep_copy(*e_prv, 0.0);
  }
  auto tc = *t_cur;
  auto ec = *e_cur;
  auto ep = *e_prv;
  auto rc = *rec;

  auto u0 = pmbp->phydro->u0;
  auto eos_ = pmbp->phydro->peos->eos_data;
  const bool gen = eos_.IsGeneral();
  auto wtemp_ = pmbp->phydro->wtemp;
  auto phicc_ = pmbp->phydro->phicc0;
  const bool etg = pmbp->phydro->use_etotgrav;
  const bool cs_ = pm->use_cubed_sphere;
  auto cosc_ = pmbp->pcoord->cos_cell;
  auto x1v_ = pmbp->pcoord->x1v;
  const Real rmin_ = rmin;
  const bool three_d = pm->three_d;
  // the ledger prints temperatures in KELVIN; the scan itself works in code units, where
  // the ratio T/T_nb is the same number either way
  const Real tk = (pmbp->punit != nullptr) ? pmbp->punit->temperature_cgs() : 1.0;

  // ---- pass 1: internal energy and temperature of every cell, one ghost layer included
  // so a cell on a block face still has six neighbours (those ghosts carry the last
  // exchanged state, which is what the operator itself saw).
  const int il = is-1, iu = ie+1;
  const int jl = (indcs.nx2 > 1) ? js-1 : js, ju = (indcs.nx2 > 1) ? je+1 : je;
  const int kl = three_d ? ks-1 : ks, ku = three_d ? ke+1 : ke;
  par_for("rsc_T", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    tc(m,k,j,i) = -1.0;
    ec(m,k,j,i) = 0.0;
    if (x1v_(m,i) < rmin_) return;
    const Real d = u0(m,IDN,k,j,i);
    if (!(d > 0.0)) return;
    Real ekin;
    if (cs_) {
      const Real c = cosc_(m,k,j);
      const Real q1 = u0(m,IM1,k,j,i);
      const Real q2 = u0(m,IM2,k,j,i);
      const Real q3 = u0(m,IM3,k,j,i);
      ekin = 0.5*(q1*q1 + (q2*q2 + q3*q3 - 2.0*c*q2*q3)/(1.0 - c*c))/d;
    } else {
      ekin = 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                  + SQR(u0(m,IM3,k,j,i)))/d;
    }
    Real ei = u0(m,IEN,k,j,i) - ekin;
    if (etg) ei -= d*phicc_(m,k,j,i);
    ec(m,k,j,i) = ei;
    if (!(ei > 0.0) || !isfinite(ei)) return;
    const Real t = eos_.Temperature(d, ei, gen ? wtemp_(m,k,j,i) : -1.0);
    if (isfinite(t) && t > 0.0) tc(m,k,j,i) = t;
  });

  // ---- pass 2: the ratio to the hottest face neighbour, and where it is worst
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  Kokkos::ValLocScalar<Real, int> rloc;
  rloc.val = -1.0;
  rloc.loc = -1;
  Kokkos::parallel_reduce("rsc_ratio",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nkji),
  KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &res) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k-ks)*nji)/nx1 + js;
    const int i = (idx - m*nkji - (k-ks)*nji - (j-js)*nx1) + is;
    const Real t = tc(m,k,j,i);
    if (!(t > 0.0)) return;
    Real tnb = 0.0;
    tnb = fmax(tnb, tc(m,k,j,i-1));
    tnb = fmax(tnb, tc(m,k,j,i+1));
    if (nx2 > 1) {
      tnb = fmax(tnb, tc(m,k,j-1,i));
      tnb = fmax(tnb, tc(m,k,j+1,i));
    }
    if (nx3 > 1) {
      tnb = fmax(tnb, tc(m,k-1,j,i));
      tnb = fmax(tnb, tc(m,k+1,j,i));
    }
    if (!(tnb > 0.0)) return;
    const Real r = t/tnb;
    if (r > res.val) { res.val = r; res.loc = idx; }
  }, Kokkos::MaxLoc<Real, int>(rloc));

  // the winning cell of THIS rank, and the watched cell's increment, in one small kernel
  int lm = -1, lk = -1, lj = -1, li = -1;
  if (rloc.loc >= 0) {
    lm = rloc.loc/nkji;
    lk = (rloc.loc - lm*nkji)/nji + ks;
    lj = (rloc.loc - lm*nkji - (lk-ks)*nji)/nx1 + js;
    li = (rloc.loc - lm*nkji - (lk-ks)*nji - (lj-js)*nx1) + is;
  }
  // latch a watched cell for the cycle if there is none yet
  if (w_cycle != pm->ncycle) {
    w_cycle = pm->ncycle;
    w_m = lm; w_k = lk; w_j = lj; w_i = li;
  }
  const int wm = w_m, wk = w_k, wj = w_j, wi = w_i;
  const int dm = lm, dk = lk, dj = lj, di = li;
  const bool have = (dm >= 0);
  const bool havew = (wm >= 0);
  if (have || havew) {
    par_for("rsc_rec", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      for (int q=0; q<32; ++q) rc(q) = 0.0;
      if (have) {
        rc(0) = x1v_(dm,di);
        rc(1) = u0(dm,IDN,dk,dj,di);
        rc(2) = tc(dm,dk,dj,di);
        rc(3) = ec(dm,dk,dj,di);
        rc(4) = ec(dm,dk,dj,di) - ep(dm,dk,dj,di);
        Real tnb = 0.0;
        tnb = fmax(tnb, tc(dm,dk,dj,di-1));
        tnb = fmax(tnb, tc(dm,dk,dj,di+1));
        tnb = fmax(tnb, tc(dm,dk,dj-1,di));
        tnb = fmax(tnb, tc(dm,dk,dj+1,di));
        tnb = fmax(tnb, tc(dm,dk-1,dj,di));
        tnb = fmax(tnb, tc(dm,dk+1,dj,di));
        rc(5) = tnb;
        rc(6) = u0(dm,IM1,dk,dj,di)/u0(dm,IDN,dk,dj,di);
      }
      if (havew) {
        rc(10) = tc(wm,wk,wj,wi);
        rc(11) = ec(wm,wk,wj,wi);
        rc(12) = ec(wm,wk,wj,wi) - ep(wm,wk,wj,wi);
        rc(13) = u0(wm,IDN,wk,wj,wi);
        rc(14) = x1v_(wm,wi);
        Real tnb = 0.0;
        tnb = fmax(tnb, tc(wm,wk,wj,wi-1));
        tnb = fmax(tnb, tc(wm,wk,wj,wi+1));
        tnb = fmax(tnb, tc(wm,wk,wj-1,wi));
        tnb = fmax(tnb, tc(wm,wk,wj+1,wi));
        tnb = fmax(tnb, tc(wm,wk-1,wj,wi));
        tnb = fmax(tnb, tc(wm,wk+1,wj,wi));
        rc(15) = tnb;
      }
      // the 7-cell state of the worst cell, for the first-crossing dump
      if (have) {
        // order: centre, i-1, i+1, j-1, j+1, k-1, k+1
        const int ok[7] = {0, 0, 0, 0, 0, -1, 1};
        const int oj[7] = {0, 0, 0, -1, 1, 0, 0};
        const int oi[7] = {0, -1, 1, 0, 0, 0, 0};
        for (int q=0; q<7; ++q) {
          const int kk = dk + ok[q], jj = dj + oj[q], ii = di + oi[q];
          rc(16+q) = tc(dm,kk,jj,ii);
          rc(23+q) = u0(dm,IDN,kk,jj,ii);
        }
      }
    });
  }
  auto hr = Kokkos::create_mirror_view(rc);
  Kokkos::deep_copy(hr, rc);
  Kokkos::deep_copy(ep, ec);   // this boundary becomes the baseline for the next

  // ---- the global winner across ranks, so only one line per boundary is printed
  Real gval = have ? rloc.val : -1.0;
  int grank = global_variable::my_rank;
#if MPI_PARALLEL_ENABLED
  {
    struct { double v; int r; } lin, out;
    lin.v = static_cast<double>(gval);
    lin.r = global_variable::my_rank;
    MPI_Allreduce(&lin, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    gval = static_cast<Real>(out.v);
    grank = out.r;
  }
#endif
  if (!(gval > ratio_print)) return;
  // THE LINE BUDGET MUST BE THE SAME NUMBER ON EVERY RANK.  It used to be incremented
  // only on the rank that owned the winning cell, so the "lines >= maxlines" early
  // return at the top of Scan fired on that rank alone -- it skipped the MPI_Allreduce
  // above while all the other ranks entered it, and the run deadlocked there (the same
  // shape as the event-log deadlock of 27da6380).  gval is global, so counting here,
  // before the owner filter, makes the budget -- and hence the early return -- identical
  // on all ranks.
  ++lines;
  if (global_variable::my_rank != grank) return;
  const int gid = pmbp->gids + ((dm >= 0) ? dm : 0);
  std::cout << "### runaway [" << opname << "] cycle " << pm->ncycle
            << " t = " << pm->time
            << " | max T/T_nb = " << gval << " at gid " << gid << " ("
            << dk << "," << dj << "," << di << ") r = " << hr(0)
            << " T = " << hr(2)*tk << " T_nb = " << hr(5)*tk << " rho = " << hr(1)
            << " v_r = " << hr(6)
            << " | watch gid " << (pmbp->gids + ((wm >= 0) ? wm : 0)) << " ("
            << wk << "," << wj << "," << wi << ") T = " << hr(10)*tk
            << " T_nb = " << hr(15)*tk << " e = " << hr(11)
            << " de = " << hr(12)
            << " de/e = " << ((hr(11) != 0.0) ? hr(12)/hr(11) : 0.0)
            << std::endl;
  // first crossing of the state threshold: the seven-cell picture, once per cell
  if (gval > ratio_state && dm >= 0 &&
      static_cast<int>(dumped.size()) < maxcells) {
    const int64_t key = ((static_cast<int64_t>(gid)*512 + di)*64 + dk)*64 + dj;
    if (dumped.find(key) == dumped.end()) {
      dumped.insert(key);
      std::cout << "    runaway FIRST CROSSING gid " << gid << " (" << dk << ","
                << dj << "," << di << ") at [" << opname << "] cycle " << pm->ncycle
                << std::endl << "      T   (c,i-1,i+1,j-1,j+1,k-1,k+1) =";
      for (int q=0; q<7; ++q) std::cout << " " << hr(16+q)*tk;
      std::cout << std::endl << "      rho (same order)              =";
      for (int q=0; q<7; ++q) std::cout << " " << hr(23+q);
      std::cout << std::endl;
    }
  }
}

}  // namespace runaway_scan

#endif  // UTILS_RUNAWAY_SCAN_HPP_
