//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter_rt_diag.cpp
//! \brief deep_hot_jupiter_rt problem generator (tusplit-1102 piece): diagnostics -- the
//! energy/mass budget (EBud*), the flux and angular-momentum history columns, the
//! photosphere/OLR dumps and the per-cycle diagnostic.  Moved verbatim from
//! deep_hot_jupiter_rt.cpp.

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "deep_hot_jupiter_rt.hpp"




//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// problem/budget_dt > 0 (DIAGNOSTIC, branch dhj-budget-0929): time-integrated per-shell
// energy and mass ledger.  Every MeshBlock spans the full x1 range, so the shell index
// is i-is.  rk2 only: a change made in stage 1 survives with weight gam0(2) = 0.5, one
// made in stage 2 or in the operator split with weight 1; a stage's flux reaches the
// final state with weight 0.5 dt in both stages.  Rows of acc (index = shell or face):
//   0 FXE  x1 face energy flux fx1(IEN) A, weight 0.5 dt (incl. rho v Phi)      [faces]
//   1 FXM  x1 face mass flux                                                    [faces]
//   2 HDVE x2/x3 energy flux divergence sum (A F) over the cell, 0.5 dt         [shell]
//   3 HDVM x2/x3 mass flux divergence                                           [shell]
//   4 SRCE energy change of the whole user source (SourceFunc), weighted        [shell]
//   5 SRCM mass change of SourceFunc (must be 0)                                [shell]
//   6 FLRE ConToPrim floor energy (stage-weighted + split), from c2p_acc         [shell]
//   7 FLRM ConToPrim floor mass                                                 [shell]
//   8 QCK  energy change across the operator-split ck call                      [shell]
//   9 FIR  ck net thermal flux of the arrays that call leaves, A dt (i >= icut) [faces]
//  10 QSW  ck stellar deposit of those arrays, V dt                             [shell]
//  11 [0..3] top face: int F_E dt, int F_M dt of stage 1, then of stage 2 (unweighted)
//  12 [0..3] top face applied (0.5 dt): E, M over inflow faces (F_M < 0), then outflow
//  13 FLRS the part of FLRE made by the split conversion after the ck call      [shell]
//  14..18 CKX (= QCK - booked A F dt, V Qb dt) per cell class, summed per shell:
//      14 full call, cell ends at e_floor;  15 full call, cell on the ck_impl_demax
//      bound;  16 full call, other;  17 linearised cadence step, |de| capped;  18 lin,
//      other.  Their sum over the rows is QCK - (FIR_i - FIR_i+1 + QSW_i).
//  19 [0..4] cell counts of classes 14..18, [5] full calls, [6] lin steps
// The dump (file ebudget.txt, rank 0) is taken right BEFORE the ck call of a cycle,
// after the floors so far were drained: time, the cumulative rows, then the shell
// sums E = sum u0(IEN) V and M = sum u0(IDN) V.  Ledger per shell over two dumps:
//   dE_i = FXE_i - FXE_i+1 - HDVE_i + SRCE_i + FLRE_i + QCK_i  (and the same for M).

void EBudShellEM(Mesh *pm, DvceArray2D<Real> t) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto u0 = pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  const int nmb1 = pmbp->nmb_thispack - 1;
  Kokkos::deep_copy(t, 0.0);
  par_for("eb_em", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Kokkos::atomic_add(&t(0, i-is), u0(m,IEN,k,j,i)*vol(m,k,j,i));
    Kokkos::atomic_add(&t(1, i-is), u0(m,IDN,k,j,i)*vol(m,k,j,i));
  });
}

// the ConToPrim floor accumulator (c2p_acc, weighted by c2p_track::Weight) into the
// ledger and into the pending sums the flux_hst_floor columns read; then zero it
void EBudDrainFloors(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto c2p = pmbp->phydro->c2p_acc;
  auto acc = *ebud::acc;
  auto fh = *ebud::fhst;
  const bool st = ebud::started;
  const int n = static_cast<int>(c2p.extent(1));
  par_for("eb_flr", DevExeSpace(), 0, n-1,
  KOKKOS_LAMBDA(const int i) {
    for (int r=0; r<4; ++r) fh(r,i) += c2p(r,i);
    if (st) {
      acc(6,i) += c2p(0,i) + c2p(2,i);
      acc(7,i) += c2p(1,i) + c2p(3,i);
      acc(13,i) += c2p(2,i);
    }
    for (int r=0; r<4; ++r) c2p(r,i) = 0.0;
  });
}

// the flux terms of one RK stage (called at the top of SourceFunc: uflx holds the
// stage's fluxes as RKUpdate applied them)
void EBudStage(Mesh *pm, const Real bdt) {
  if (!ebud::on || !ebud::started) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto acc = *ebud::acc;
  auto fx1 = pmbp->phydro->uflx.x1f;
  auto fx2 = pmbp->phydro->uflx.x2f;
  auto fx3 = pmbp->phydro->uflx.x3f;
  auto a1 = pmbp->pcoord->area.x1f;
  auto a2 = pmbp->pcoord->area.x2f;
  auto a3 = pmbp->pcoord->area.x3f;
  const Real w = 0.5*pm->dt;
  const Real dtf = pm->dt;
  const int r11 = (bdt > 0.75*pm->dt) ? 0 : 2;   // rk2: beta = 1 (stage 1), 0.5 (2)
  const int nmb1 = pmbp->nmb_thispack - 1;
  par_for("eb_stage", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real fe = fx1(m,IEN,k,j,i)*a1(m,k,j,i);
    const Real fm = fx1(m,IDN,k,j,i)*a1(m,k,j,i);
    Kokkos::atomic_add(&acc(0, i-is), fe*w);
    Kokkos::atomic_add(&acc(1, i-is), fm*w);
    if (i <= ie) {
      const Real he = fx2(m,IEN,k,j+1,i)*a2(m,k,j+1,i) - fx2(m,IEN,k,j,i)*a2(m,k,j,i)
                    + fx3(m,IEN,k+1,j,i)*a3(m,k+1,j,i) - fx3(m,IEN,k,j,i)*a3(m,k,j,i);
      const Real hm = fx2(m,IDN,k,j+1,i)*a2(m,k,j+1,i) - fx2(m,IDN,k,j,i)*a2(m,k,j,i)
                    + fx3(m,IDN,k+1,j,i)*a3(m,k+1,j,i) - fx3(m,IDN,k,j,i)*a3(m,k,j,i);
      Kokkos::atomic_add(&acc(2, i-is), he*w);
      Kokkos::atomic_add(&acc(3, i-is), hm*w);
    } else {
      Kokkos::atomic_add(&acc(11, r11), fe*dtf);
      Kokkos::atomic_add(&acc(11, r11+1), fm*dtf);
      const int r12 = (fm < 0.0) ? 0 : 2;
      Kokkos::atomic_add(&acc(12, r12), fe*w);
      Kokkos::atomic_add(&acc(12, r12+1), fm*w);
    }
  });
}

void EBudWrite(Mesh *pm, const Real t) {
  const int n = pm->mb_indcs.nx1 + 1;
  auto tm = *ebud::tmp;
  EBudShellEM(pm, tm);
  auto ha = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), *ebud::acc);
  auto ht = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), tm);
  std::vector<Real> buf((ebud::NQ + 2)*n);
  for (int q=0; q<ebud::NQ; ++q)
    for (int i=0; i<n; ++i) buf[q*n+i] = ha(q, i);
  for (int q=0; q<2; ++q)
    for (int i=0; i<n; ++i) buf[(ebud::NQ+q)*n+i] = ht(q, i);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, buf.data(), static_cast<int>(buf.size()),
                MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    FILE *fp = std::fopen("ebudget.txt", "a");
    std::fprintf(fp, "# t= %.17e ncycle= %d\n", t, static_cast<int>(pm->ncycle));
    for (int i=0; i<n; ++i) {
      std::fprintf(fp, "%d", i);
      for (int q=0; q<ebud::NQ+2; ++q) std::fprintf(fp, " %.17e", buf[q*n+i]);
      std::fprintf(fp, "\n");
    }
    std::fclose(fp);
  }
}

// wraps the operator-split ck call (user_split_once: once per cycle, full dt)
void EBudCk(Mesh *pm, const Real dt, const bool before) {
  if (!ebud::on) return;
  namespace ts = two_stream_rt;
  auto tm = *ebud::tmp;
  const Real tstate = pm->time + dt;
  if (before) {
    EBudDrainFloors(pm);
    if (!ebud::started) {
      Kokkos::deep_copy(*ebud::acc, 0.0);
      ebud::started = true;
      ebud::t_next = tstate + ebud::dt_out;
      EBudWrite(pm, tstate);
    } else if (tstate >= ebud::t_next) {
      EBudWrite(pm, tstate);
      ebud::t_next += ebud::dt_out;
    }
    EBudShellEM(pm, tm);   // QCK 'before' (row 0)
    if (ebud::cpre != nullptr) {
      MeshBlockPack *pb = pm->pmb_pack;
      auto &ib = pm->mb_indcs;
      auto cp = *ebud::cpre;
      auto u0b = pb->phydro->u0;
      const ts::CkCadState st = ts::CkCadMakeState(pm);
      par_for("eb_cpre", DevExeSpace(), 0, pb->nmb_thispack - 1, ib.ks, ib.ke, ib.js,
              ib.je, ib.is, ib.ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        cp(m,0,k,j,i) = u0b(m,IEN,k,j,i);
        cp(m,1,k,j,i) = st.Eint(m, k, j, i);
      });
      ebud::nfull0 = ts::ck_cad_nfull;
    }
    return;
  }
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto acc = *ebud::acc;
  auto u0 = pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  auto a1 = pmbp->pcoord->area.x1f;
  const int nmb1 = pmbp->nmb_thispack - 1;
  par_for("eb_ck", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Kokkos::atomic_add(&acc(8, i-is), u0(m,IEN,k,j,i)*vol(m,k,j,i));
  });
  par_for("eb_ck0", DevExeSpace(), 0, ie-is,
  KOKKOS_LAMBDA(const int i) {
    acc(8, i) -= tm(0, i);
  });
  const bool rad = ts::rt_face_flux_ready() && (ts::rt_Qb_ptr != nullptr)
                   && (ts::rt_icut_ptr != nullptr);
  if (rad && ebud::cpre != nullptr) {
    // CKX per cell class (rows 14..19)
    auto fb = ts::rt_face_flux();
    auto qb = *ts::rt_Qb_ptr;
    auto icut = ts::rt_cut_index();
    const int nblk = ts::rt_face_nblk();
    auto cp = *ebud::cpre;
    const ts::CkCadState st = ts::CkCadMakeState(pm);
    const bool full = (ts::ck_cad_nfull != ebud::nfull0) || (ts::ck_impl_every <= 1);
    const Real demax = ts::ck_impl_demax;
    Real capf = ts::ck_impl_dtmax;
    if (demax > 0.0 && demax < capf) capf = demax;
    par_for("eb_ckx", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int ic = icut(m,k,j);
      Real bk = 0.0;
      if (i >= ic) {
        Real f0 = 0.0, f1 = 0.0, q = 0.0;
        for (int b=0; b<nblk; ++b) {
          f0 += fb(m,b,i,k,j);
          f1 += fb(m,b,i+1,k,j);
          q += qb(m,b,i,k,j);
        }
        bk = (f0*a1(m,k,j,i) - f1*a1(m,k,j,i+1) + q*vol(m,k,j,i))*dt;
      }
      const Real x = (u0(m,IEN,k,j,i) - cp(m,0,k,j,i))*vol(m,k,j,i) - bk;
      const Real eb = cp(m,1,k,j,i);
      const Real ea = st.Eint(m, k, j, i);
      const Real d = u0(m,IDN,k,j,i);
      int c = 0;
      if (full) {
        Real efl = st.eos.EnergyFromTemperature(d, st.eos.tfloor);
        const Real ep = st.eos.EnergyFromPressure(d, st.eos.pfloor);
        if (ep > efl) efl = ep;
        if (ea <= efl*(1.0 + 1.0e-8)) {
          c = 0;
        } else if (demax > 0.0 && eb > 0.0 &&
                   (ea <= eb*(1.0 - demax)*(1.0 + 1.0e-8) ||
                    ea >= eb*(1.0 + demax)*(1.0 - 1.0e-8))) {
          c = 1;
        } else {
          c = 2;
        }
      } else {
        c = (eb > 0.0 && fabs(ea - eb) >= capf*eb*(1.0 - 1.0e-8)) ? 3 : 4;
      }
      Kokkos::atomic_add(&acc(14 + c, i-is), x);
      Kokkos::atomic_add(&acc(19, c), 1.0);
    });
    par_for("eb_ckn", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      acc(19, full ? 5 : 6) += 1.0;
    });
  }
  if (rad) {
    auto fb = ts::rt_face_flux();
    auto qb = *ts::rt_Qb_ptr;
    auto icut = ts::rt_cut_index();
    const int nblk = ts::rt_face_nblk();
    par_for("eb_rad", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int ic = icut(m,k,j);
      if (i < ic) return;
      Real f = 0.0, q = 0.0;
      for (int b=0; b<nblk; ++b) {
        f += fb(m,b,i,k,j);
        if (i <= ie) q += qb(m,b,i,k,j);
      }
      Kokkos::atomic_add(&acc(9, i-is), f*a1(m,k,j,i)*dt);
      if (i <= ie) Kokkos::atomic_add(&acc(10, i-is), q*vol(m,k,j,i)*dt);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DhjFluxStage
//! \brief problem/flux_hst_rkavg: add this RK stage's fluid face fluxes at the top face
//! and at the inner wall, times the stage's weight, into dhj_fhst (reset on stage 1)

void DhjFluxStage(Mesh *pm) {
  if (!dhj_fhst_avg) return;
  if (pm->ncycle != dhj_fhst_cyc) {
    dhj_fhst_cyc = pm->ncycle;
    dhj_fhst_stage = 0;
  }
  if (dhj_fhst_stage >= dhj_fhst_nst) return;     // never: one call per stage
  const Real w = dhj_fhst_w[dhj_fhst_stage];
  const bool first = (dhj_fhst_stage == 0);
  ++dhj_fhst_stage;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto fx1 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
  auto acc = *dhj_fhst_ptr;
  par_for("dhj_fhst_stage", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real v[4] = {fx1(m,IEN,k,j,ie+1), fx1(m,IDN,k,j,ie+1), fx1(m,IEN,k,j,is),
                       fx1(m,IDN,k,j,is)};
    for (int n=0; n<4; ++n) acc(m,n,k,j) = (first ? 0.0 : acc(m,n,k,j)) + w*v[n];
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmHstStage
//! \brief problem/am_hst: one RK stage's torque sums (AmCell, from the stage's w0 and
//! uflx: what SourceFunc and the update use), added to amh::acc with the stage's weight
//! in the cycle (w_s dt; rk2 1/2, 1/2), so that acc / (t - t0) is the rate the update
//! applied over the history interval.
//!   0 Tcc        the cell-centre Coriolis torque, sum -2 Omega rho R v_R vol
//!   1 Tcont      -Omega sum_c Q_c (namespace amh) = -Omega dI/dt - Omega (W outflow)
//!   2 Tbnd       the boundary part of Tcont, -Omega sum Phi^out W_f, x1 faces is, ie+1
//!   3 Tdif_deep  Tcc - Tcont over cells with p > problem/am_hst_pdeep (100 bar)
//!   4 Tdif_rad   the radial-velocity part of Tcc minus the radial-face part of Tcont
//!   5-9          (am_hst_split) T_code x1, x2, x3, metric source, and T_cons

void AmHstStage(Mesh *pm) {
  if (pm->ncycle != amh::cyc) {
    amh::cyc = pm->ncycle;
    amh::stage = 0;
  }
  if (amh::stage >= amh::nst) return;
  const Real wdt = amh::w[amh::stage]*pm->dt;
  ++amh::stage;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const AmV v = AmMakeV(pm, amh::split);
  const Real pdeep = amh::pdeep;
  const Real omega = v.omega;
  const int ncell = nmb*nx3*nx2*nx1;
  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("dhj_am_stage", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/(nx3*nx2*nx1);
    int rem = idx - m*(nx3*nx2*nx1);
    const int k = ks + rem/(nx2*nx1);
    rem -= (k - ks)*(nx2*nx1);
    const int j = js + rem/nx1;
    const int i = is + rem - (j - js)*nx1;
    AmOut o;
    AmCell(v, m, k, j, i, o);
    const Real rho = v.w0(m,IDN,k,j,i);
    const Real p = v.gen ? v.wder(m,IDPR,k,j,i) : v.eos.Pressure(rho, v.w0(m,IEN,k,j,i));
    const Real tcont = -omega*o.q;
    array_sum::GlobalSum h;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
    h.the_array[0] = o.tcc;
    h.the_array[1] = tcont;
    h.the_array[2] = o.tbnd;
    h.the_array[3] = (p > pdeep) ? (o.tcc - tcont) : 0.0;
    h.the_array[4] = o.tccr + omega*o.qr;
    for (int n=0; n<4; ++n) h.the_array[5+n] = o.tc[n];
    h.the_array[9] = o.tcons;
    h.the_array[10] = o.twall;
    mb_sum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  Kokkos::fence();
  for (int n=0; n<11; ++n) amh::acc[n] += wdt*sum.the_array[n];
  // relaxdrag-1101 Tdrag: the drag removes itdrag bdt of every momentum component of u0
  // in SourceFunc, so its z torque is -itdrag Lrel(u0) (same time and u0 as the kernel)
  const Real itd = DhjRelaxDragRate(pm->time, omega);
  if (itd > 0.0) amh::acc[11] += wdt*(-itd*DhjLrelU0(pm));
}

//----------------------------------------------------------------------------------------
//! \fn DhjRelaxDragRate (relaxdrag-1101)
//! \brief 1/tau of the start-up Rayleigh drag at time t, as SourceFunc applies it (0 when
//! off or after relax_drag_rot rotations)

Real DhjRelaxDragRate(const Real time, const Real omega) {
  if (!dhj_relax_drag || omega == 0.0) return 0.0;
  const Real rdf = dhj_relax_drag_rot/5.0;
  const Real dyntime = 2.0*M_PI/omega;
  const Real tau1 = dyntime/10.0, tau2 = dyntime;
  const Real t1 = 2.0*dyntime*rdf, t2 = 5.0*dyntime*rdf;
  if (time < t1) return 1.0/tau1;
  if (time < t2) return 1.0/tau1*pow(tau1/tau2, (time - t1)/(t2 - t1));
  return 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn DhjLrelU0 (relaxdrag-1101)
//! \brief Lrel = sum vol (x P_y - y P_x) of this rank's u0, as AmHstColumns forms it

Real DhjLrelU0(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v = pmbp->pcoord->x2v;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;
  const bool cs = pm->use_cubed_sphere;
  Real lsum = 0.0;
  Kokkos::parallel_reduce("dhj_lrel_u0",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nx3*nx2*nx1),
  KOKKOS_LAMBDA(const int idx, Real &ls) {
    const int m = idx/(nx3*nx2*nx1);
    int rem = idx - m*(nx3*nx2*nx1);
    const int k = ks + rem/(nx2*nx1);
    rem -= (k - ks)*(nx2*nx1);
    const int j = js + rem/nx1;
    const int i = is + rem - (j - js)*nx1;
    const Real r = x1v(m,i);
    Real lz;
    if (cs) {
      const Real xi_c = 0.25*M_PI*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                              size.d_view(m).x2max);
      const Real eta_c = 0.25*M_PI*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                               size.d_view(m).x3max);
      const int pnl = mbpanel.d_view(m);
      Real qc[3], e1[3], e2[3];
      cubed_sphere::PanelToCart(pnl, xi_c, eta_c, qc);
      cubed_sphere::PanelTangents(pnl, xi_c, eta_c, e1, e2);
      const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
      const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn;
      const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
      const Real det = 1.0 - c*c;
      const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
      const Real p2 = (m2 - c*m3)/det, p3 = (m3 - c*m2)/det;
      const Real px = m1*rh0 + p2*e1[0] + p3*e2[0];
      const Real py = m1*rh1 + p2*e1[1] + p3*e2[1];
      lz = r*rh0*py - r*rh1*px;
    } else {
      lz = r*sin(x2v(m,j))*u0(m,IM3,k,j,i);
    }
    ls += lz*vol(m,k,j,i);
  }, Kokkos::Sum<Real>(lsum));
  return lsum;
}



//----------------------------------------------------------------------------------------
//! \fn void DhjPhotosphereDump
//! \brief The tau = 2/3 PHOTOSPHERE, per band and g-point, from the run's OWN opacity.
//!
//! WHY THIS EXISTS.  An isobar is not the observed surface.  What a telescope sees is the
//! level where the vertical optical depth reaches ~2/3, and on an ultra-hot Jupiter that
//! level is not one surface: it sits DEEP in a window between water bands and HIGH in a
//! band centre, and it rises on the dayside where H- opacity climbs with temperature.
//! Presenting a single pressure level as "the observed surface" hides all three effects.
//!
//! The opacity is not re-derived here.  It is the identical lookup the longwave sweep
//! uses -- `ck_kappa(...) + rt_kc` -- so the photosphere reported is the one the run's
//! own radiative transfer saw, not a model of it.  Requires the correlated-k path and
//! at least one RT evaluation, because the per-cell tables are allocated lazily: restart
//! and take a single cycle.
//!
//! DEFINITION, stated because it is a choice.  Optical depth is accumulated VERTICALLY
//! (kappa*rho*dr, with no 1/mu slant factor) downwards from the top, starting from the
//! hydrostatic column above the domain that the sweep's own top boundary term adds.  For
//! each (band, g-point) the crossing tau = 2/3 is located by interpolating in log tau
//! against log p.  A band's reported level is the g-WEIGHTED MEAN of log p over its
//! g-points, i.e. the emission-weighted level; the spread between its lowest and
//! g-point is the window-to-line-core range, which is the point of the diagnostic.  The
//! broadband level weights bands by the local Planck fraction f_b(T).
//!
//! A g-point whose tau = 2/3 lies ABOVE the domain top is reported as p = 0: the model
//! does not contain that emitting level, which is a fact about the grid worth seeing
//! rather than clamping silently away.

void DhjPhotosphereDump(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const std::string fname = pin->GetOrAddString("problem", "photosphere_dump", "");
  if (fname.empty()) return;
  if (!rt_ck || ck_lk_ptr == nullptr || rt_xT_ptr == nullptr || rt_kc_ptr == nullptr) {
    if (global_variable::my_rank == 0) {
      std::cout << "### photosphere_dump: needs the correlated-k path AND at least one RT"
                << " evaluation (the per-cell tables are allocated lazily) -- restart and"
                << " run one cycle. Nothing written." << std::endl;
    }
    return;
  }
  const Real TAU_PH = 2.0/3.0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;

  auto xT = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_xT_ptr);
  auto xP = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_xP_ptr);
  auto pbr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_pb_ptr);
  auto Tc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_T_ptr);
  auto kc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_kc_ptr);
  auto cf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_cf_ptr);
  auto lk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_lk_ptr);
  auto gw = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_gw_ptr);
  auto wl = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_wl_ptr);
  auto pf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_pf_ptr);
  auto w0h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                 (pmbp->phydro != nullptr) ? pmbp->phydro->w0 : pmbp->pmhd->w0);
  auto dx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->dx1);
  auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto x2v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
  auto x3v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);

  const HotJupiterParam &hj = pm->pgen->hot_jupiter_param;
  std::ostringstream os;
  os.precision(7);
  os << std::scientific;
  if (global_variable::my_rank == 0) {
    os << "# AthenaK tau=2/3 photosphere, from the run's own correlated-k opacity\n"
       << "# t = " << pm->time << " s   cycle " << pm->ncycle << "\n"
       << "# bands " << CK_NB << "  g-points " << CK_NG << "\n"
       << "# band edges [um], descending:";
    for (int b=0; b<=CK_NB; ++b) os << " " << wl(b);
    os << "\n# lat_deg lon_deg band lam_lo_um lam_hi_um p_eff_bar p_deep_bar p_high_bar"
       << " r_eff_cm\n"
       << "# p = 0 means tau=2/3 lies above the domain top; band index " << CK_NB
       << " is the broadband, Planck-weighted level\n";
  }

  Real tau[CK_NG], plev[CK_NG], rlev[CK_NG];
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      const Real lon = (x3v(m,k) - M_PI)*180.0/M_PI;
      for (int j=js; j<=je; ++j) {
        const Real lat = 90.0 - x2v(m,j)*180.0/M_PI;
        const Real mu0 = cf(m,k,j,3);
        Real pbb[CK_NB+1], rbb[CK_NB+1], pdeep[CK_NB+1], phigh[CK_NB+1];
        for (int b=0; b<CK_NB; ++b) {
          // start each g-point from the hydrostatic column ABOVE the domain, exactly as
          // the sweep's own top boundary term does
          const Real ptop = pbr(m,k,j,ie+1);
          const int iT0 = static_cast<int>(xT(m,k,j,ie+1));
          const int iP0 = static_cast<int>(xP(m,k,j,ie+1));
          const Real fT0 = xT(m,k,j,ie+1) - static_cast<Real>(iT0);
          const Real fP0 = xP(m,k,j,ie+1) - static_cast<Real>(iP0);
          const Real gtop = EffGravAt(hj.grav, hj.ap, x1v(m,ie+1), hj.grav_point_mass,
                                      hj.omega, mu0, hj.stellar_tide)
                            *cf(m,k,j,4);   // roche-dhj-1031: |grad Phi|/g_pm, else 1
          for (int g=0; g<CK_NG; ++g) {
            const Real kap = ck_kappa(lk, iT0, fT0, iP0, fP0, b, g) + kc(m,b,ie+1,k,j);
            tau[g] = kap*ptop*1.0e6/gtop;
            plev[g] = 0.0;
            rlev[g] = 0.0;
          }
          for (int i=ie; i>=is; --i) {
            const int iT = static_cast<int>(xT(m,k,j,i));
            const int iP = static_cast<int>(xP(m,k,j,i));
            const Real fT = xT(m,k,j,i) - static_cast<Real>(iT);
            const Real fP = xP(m,k,j,i) - static_cast<Real>(iP);
            const Real drho = w0h(m,IDN,k,j,i)*dx1(m,k,j,i);
            for (int g=0; g<CK_NG; ++g) {
              if (plev[g] != 0.0) continue;
              const Real t0 = tau[g];
              tau[g] += (ck_kappa(lk, iT, fT, iP, fP, b, g) + kc(m,b,i,k,j))*drho;
              if (tau[g] >= TAU_PH && t0 < TAU_PH) {
                const Real f = (t0 > 0.0)
                    ? (std::log(TAU_PH) - std::log(t0))/(std::log(tau[g]) - std::log(t0))
                    : 1.0;
                const Real lp0 = std::log(pbr(m,k,j,i+1));
                const Real lp1 = std::log(pbr(m,k,j,i));
                plev[g] = std::exp(lp0 + f*(lp1 - lp0));
                rlev[g] = x1v(m,i+1) + f*(x1v(m,i) - x1v(m,i+1));
              }
            }
          }
          Real sw = 0.0, slp = 0.0, sr = 0.0, lo = 0.0, hi = 0.0;
          for (int g=0; g<CK_NG; ++g) {
            if (plev[g] <= 0.0) continue;
            sw += gw(g);
            slp += gw(g)*std::log(plev[g]);
            sr += gw(g)*rlev[g];
            lo = (lo == 0.0) ? plev[g] : fmax(lo, plev[g]);
            hi = (hi == 0.0) ? plev[g] : fmin(hi, plev[g]);
          }
          pbb[b] = (sw > 0.0) ? std::exp(slp/sw) : 0.0;
          rbb[b] = (sw > 0.0) ? sr/sw : 0.0;
          pdeep[b] = lo;
          phigh[b] = hi;
        }
        Real Bb[CK_NB], swb = 0.0, slb = 0.0, srb = 0.0;
        ck_planck_bands(pf, ck_pf_lTmin, ck_pf_idlT, 1.0, Tc(m,k,j,ie), Bb);
        for (int b=0; b<CK_NB; ++b) {
          if (pbb[b] <= 0.0) continue;
          swb += Bb[b];
          slb += Bb[b]*std::log(pbb[b]);
          srb += Bb[b]*rbb[b];
        }
        pbb[CK_NB] = (swb > 0.0) ? std::exp(slb/swb) : 0.0;
        rbb[CK_NB] = (swb > 0.0) ? srb/swb : 0.0;
        pdeep[CK_NB] = pbb[CK_NB];
        phigh[CK_NB] = pbb[CK_NB];
        for (int b=0; b<=CK_NB; ++b) {
          os << lat << " " << lon << " " << b << " "
             << ((b < CK_NB) ? wl(b+1) : 0.0) << " " << ((b < CK_NB) ? wl(b) : 0.0) << " "
             << pbb[b] << " " << pdeep[b] << " " << phigh[b] << " " << rbb[b] << "\n";
        }
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  {
    // gather the per-rank text on rank 0, in rank order
    std::string mine = os.str();
    int len = static_cast<int>(mine.size());
    int nrank = global_variable::nranks;
    std::vector<int> lens(nrank), offs(nrank, 0);
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    int tot = 0;
    if (global_variable::my_rank == 0) {
      for (int r=0; r<nrank; ++r) { offs[r] = tot; tot += lens[r]; }
    }
    std::vector<char> all((global_variable::my_rank == 0) ? tot : 1);
    MPI_Gatherv(mine.data(), len, MPI_CHAR, all.data(), lens.data(), offs.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      std::ofstream f(fname);
      f.write(all.data(), tot);
      f.close();
      std::cout << "photosphere: written to '" << fname << "'" << std::endl;
    }
  }
#else
  {
    std::ofstream f(fname);
    f << os.str();
    f.close();
    std::cout << "photosphere: written to '" << fname << "'" << std::endl;
  }
#endif
  return;
}


void DhjCycleDiag(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (diag_gid < pmbp->gids || diag_gid > pmbp->gide) return;   // not this rank
  if (pm->ncycle < diag_cycle_min || pm->ncycle > diag_cycle_max) return;
  const int m = diag_gid - pmbp->gids;

  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nf = n1 + 1;

  auto &u0 = (pmbp->phydro != nullptr) ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  auto &w0 = (pmbp->phydro != nullptr) ? pmbp->phydro->w0 : pmbp->pmhd->w0;
  Conduction *pc = (pmbp->phydro != nullptr) ? pmbp->phydro->pcond
                                             : ((pmbp->pmhd != nullptr)
                                                ? pmbp->pmhd->pcond : nullptr);
  const int nw = static_cast<int>(w0.extent_int(1));
  const int nu = static_cast<int>(u0.extent_int(1));

  if (!diag_dir_made) {
    mkdir("cyclediag", 0775);
    diag_dir_made = true;
  }
  std::ostringstream fn;
  fn << "cyclediag/" << diag_basename << ".cyclediag." << std::setw(8)
     << std::setfill('0') << pm->ncycle << ".dat";
  // write to <name>.tmp and rename() only after close(): a rank killed part-way through
  // by another rank's MPI_Abort then leaves no truncated file under the final name
  const std::string fname = fn.str();
  const std::string ftmp = fname + ".tmp";
  std::ofstream f(ftmp.c_str(), std::ios::binary);
  if (!f.is_open()) {
    std::cout << "### cyclediag: could not open '" << ftmp << "'" << std::endl;
    return;
  }

  // ---- the variable list, fixed whatever is allocated -------------------------------
  std::vector<std::string> names;
  const char *wn[5] = {"w_dens", "w_velx", "w_vely", "w_velz", "w_eint"};
  const char *un[5] = {"u_dens", "u_m1", "u_m2", "u_m3", "u_ener"};
  for (int n=0; n<nw; ++n) {
    if (n < 5) {
      names.push_back(wn[n]);
    } else {
      names.push_back("w_" + std::to_string(n));
    }
  }
  for (int n=0; n<nu; ++n) {
    if (n < 5) {
      names.push_back(un[n]);
    } else {
      names.push_back("u_" + std::to_string(n));
    }
  }
  names.push_back("rt_src");
  names.push_back("rt_de");
  names.push_back("rt_clip");
  names.push_back("cond_f1");
  names.push_back("cond_f2");
  names.push_back("cond_f3");
  names.push_back("cond_kappa");
  names.push_back("cond_keff");
  names.push_back("cond_dtcell");
  names.push_back("rt_T");
  names.push_back("rt_icut");
  names.push_back("rad_w");
  names.push_back("rad_tauf");
  // appended AFTER the v1 set, so a 23-variable and a 27-variable file share their
  // first 23 variables byte for byte and one name-driven reader handles both
  names.push_back("rt_Ft");
  names.push_back("rt_Fb");
  names.push_back("rt_Qs");
  names.push_back("rt_Em");

  std::ostringstream hd;
  hd.precision(17);
  hd << "cyclediag v1 cycle=" << pm->ncycle << " time=" << pm->time
     << " dt=" << pm->dt << " gid=" << diag_gid
     << " n1=" << n1 << " n2=" << n2 << " n3=" << n3
     << " is=" << indcs.is << " ie=" << indcs.ie
     << " js=" << indcs.js << " je=" << indcs.je
     << " ks=" << indcs.ks << " ke=" << indcs.ke
     << " nvar=" << names.size() << " nface=" << nf << "\n";
  f << hd.str();
  for (size_t n=0; n<names.size(); ++n) {
    f << names[n] << ((n+1 == names.size()) ? "\n" : " ");
  }

  // ---- writers ----------------------------------------------------------------------
  std::vector<double> buf;
  auto put3 = [&](const int nk, const int nj, const int ni, auto hv) {
    buf.assign(static_cast<size_t>(nk)*nj*ni, 0.0);
    size_t p = 0;
    for (int k=0; k<nk; ++k) {
      for (int j=0; j<nj; ++j) {
        for (int i=0; i<ni; ++i) {
          buf[p++] = static_cast<double>(hv(k,j,i));
        }
      }
    }
    f.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()*sizeof(double)));
  };
  auto zeros = [&](const size_t nel) {
    buf.assign(nel, 0.0);
    f.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()*sizeof(double)));
  };
  // one variable out of a (m,n,k,j,i) device array, via a host mirror of the subview
  auto put5 = [&](const DvceArray5D<Real> &a, const int n) {
    auto s = Kokkos::subview(a, m, n, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), s);
    put3(static_cast<int>(h.extent_int(0)), static_cast<int>(h.extent_int(1)),
         static_cast<int>(h.extent_int(2)), h);
  };

  const size_t ncell = static_cast<size_t>(n3)*n2*n1;
  for (int n=0; n<nw; ++n) put5(w0, n);
  for (int n=0; n<nu; ++n) put5(u0, n);

  if (rt_diag_ptr != nullptr) {
    for (int n=0; n<3; ++n) put5(*rt_diag_ptr, n);
  } else {
    zeros(3*ncell);
  }

  if (pc != nullptr && pc->diag && pc->cond_diag.extent_int(0) > m) {
    for (int n=0; n<6; ++n) put5(pc->cond_diag, n);
  } else {
    zeros(6*ncell);
  }

  // temperature [K], from the RT's own per-cell table
  if (rt_T_ptr != nullptr) {
    auto s = Kokkos::subview(*rt_T_ptr, m, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), s);
    put3(static_cast<int>(h.extent_int(0)), static_cast<int>(h.extent_int(1)),
         static_cast<int>(h.extent_int(2)), h);
  } else {
    zeros(ncell);
  }

  // the correlated-k cut index, one per column, promoted to Real
  if (rt_icut_ptr != nullptr) {
    auto s = Kokkos::subview(*rt_icut_ptr, m, Kokkos::ALL(), Kokkos::ALL());
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), s);
    buf.assign(static_cast<size_t>(n3)*n2, 0.0);
    size_t p = 0;
    for (int k=0; k<static_cast<int>(h.extent_int(0)); ++k) {
      for (int j=0; j<static_cast<int>(h.extent_int(1)); ++j) {
        buf[p++] = static_cast<double>(h(k,j));
      }
    }
    f.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()*sizeof(double)));
  } else {
    zeros(static_cast<size_t>(n3)*n2);
  }

  // the blend weight and the face optical depth, both on x1 faces (n3,n2,n1+1)
  const size_t nfc = static_cast<size_t>(n3)*n2*nf;
  if (pc != nullptr && pc->rad_tau_mode && pc->rad_w.extent_int(0) > m) {
    auto sw = Kokkos::subview(pc->rad_w, m, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto hw = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sw);
    put3(static_cast<int>(hw.extent_int(0)), static_cast<int>(hw.extent_int(1)),
         static_cast<int>(hw.extent_int(2)), hw);
    auto st = Kokkos::subview(pc->rad_tauf, m, Kokkos::ALL(), Kokkos::ALL(),
                              Kokkos::ALL());
    auto ht = Kokkos::create_mirror_view_and_copy(HostMemSpace(), st);
    put3(static_cast<int>(ht.extent_int(0)), static_cast<int>(ht.extent_int(1)),
         static_cast<int>(ht.extent_int(2)), ht);
  } else {
    zeros(2*nfc);
  }

  // the flux components behind rt_src: Ft, Fb, Qs, Em (v2, appended at the end)
  if (rt_diag_ptr != nullptr && rt_diag_ptr->extent_int(1) >= 7) {
    for (int n=3; n<7; ++n) put5(*rt_diag_ptr, n);
  } else {
    zeros(4*ncell);
  }

  f.close();
  if (std::rename(ftmp.c_str(), fname.c_str()) != 0) {
    std::cout << "### cyclediag: could not rename '" << ftmp << "'" << std::endl;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void DhjFinalDumps / DhjOlrDump
//! \brief problem/olr_dump (default off): per column, the NET longwave flux through the
//! top face ie+1 per ck band (rt_Fb blocks summed per band: chain
//! c = ((b CK_NG)+g) nq + q, RT_NB chains per block), the cos of the stellar zenith
//! angle and the absorbed stellar flux sum_i Q_sw V / A(ie+1).  Code flux units (cgs).
//! Reads the arrays of the LAST ck call, like DhjFluxHistory; nothing is re-run.
//! Rows: gid k-ks j-js mu0
//! F_sw_abs F_band[0..nb-1] (bands in table order, edges in the header).

void DhjFinalDumps(ParameterInput *pin, Mesh *pm) {
  if (!pin->GetOrAddString("problem", "photosphere_dump", "").empty()) {
    DhjPhotosphereDump(pin, pm);
  }
  DhjOlrDump(pin, pm);
}

void DhjOlrDump(ParameterInput *pin, Mesh *pm) {
  namespace ts = two_stream_rt;
  const std::string fname = pin->GetOrAddString("problem", "olr_dump", "");
  if (fname.empty()) return;
  if (!ts::rt_face_flux_ready() || ts::rt_Qb_ptr == nullptr || ts::rt_cf_ptr == nullptr
      || !rt_ck || ck_wl_ptr == nullptr) {
    if (global_variable::my_rank == 0) {
      std::cout << "### olr_dump: needs the ck split path and one RT call;"
                << " nothing written"
                << std::endl;
    }
    return;
  }
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto fb = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ts::rt_face_flux());
  auto qb = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ts::rt_Qb_ptr);
  auto cf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ts::rt_cf_ptr);
  auto ic = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ts::rt_cut_index());
  auto wl = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_wl_ptr);
  auto area = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                  pmbp->pcoord->area.x1f);
  auto vol = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->volume);
  const int nblk = fb.extent_int(1);
  const int nqb = nblk/CK_NB;          // blocks per band
  std::ostringstream os;
  os.precision(9);
  os << std::scientific;
  if (global_variable::my_rank == 0) {
    os << "# AthenaK olr_dump: net top longwave flux per ck band [erg/s/cm2]\n"
       << "# t = " << pm->time << " s   cycle " << pm->ncycle << "   bands " << CK_NB
       << "   blocks " << nblk << "\n# band edges [um] (table order):";
    for (int b=0; b<=CK_NB; ++b) os << " " << wl(b);
    os << "\n# gid kk jj mu0 F_sw_abs F_band[0.." << CK_NB-1 << "]\n";
  }
  for (int m=0; m<nmb; ++m) {
    const int gid = pmbp->pmb->mb_gid.h_view(m);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        const Real at = area(m,k,j,ie+1);
        Real qsw = 0.0;
        for (int i=ic(m,k,j); i<=ie; ++i) {
          for (int b=0; b<nblk; ++b) qsw += qb(m,b,i,k,j)*vol(m,k,j,i);
        }
        os << gid << " " << k-ks << " " << j-js << " " << cf(m,k,j,3) << " " << qsw/at;
        for (int b=0; b<CK_NB; ++b) {
          Real f = 0.0;
          for (int q=0; q<nqb; ++q) f += fb(m,b*nqb+q,ie+1,k,j);
          os << " " << f;
        }
        os << "\n";
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  {
    std::string mine = os.str();
    int len = static_cast<int>(mine.size());
    int nrank = global_variable::nranks;
    std::vector<int> lens(nrank), offs(nrank, 0);
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    int tot = 0;
    if (global_variable::my_rank == 0) {
      for (int r=0; r<nrank; ++r) { offs[r] = tot; tot += lens[r]; }
    }
    std::vector<char> all((global_variable::my_rank == 0) ? tot : 1);
    MPI_Gatherv(mine.data(), len, MPI_CHAR, all.data(), lens.data(), offs.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      std::ofstream f(fname);
      f.write(all.data(), tot);
      f.close();
      std::cout << "olr_dump: written to '" << fname << "'" << std::endl;
    }
  }
#else
  {
    std::ofstream f(fname);
    f << os.str();
    f.close();
    std::cout << "olr_dump: written to '" << fname << "'" << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void DhjFluxHistory
//! \brief problem/flux_hst: area-integrated flux columns, appended after the hydro/MHD
//! history columns.  All in cgs (erg/s and g/s; the dhj runs use cgs code units), summed
//! over the whole sphere (MPI_SUM over ranks by the history output).  Nothing re-runs the
//! solver: the radiative columns read the correlated-k two-stream's own arrays, i.e. the
//! LAST ck call that touched each column (under ck_impl_every > 1 the last full call,
//! or a later guard call on that column; within the call, its final pass).  The fluid
//! columns read the RK-weighted sum of the stages' Riemann fluxes
//! (problem/flux_hst_rkavg, default since dhj-ck-conserve 09-29; false = the last RK
//! stage, as before; the flux_hst_wall split columns always read the last stage).
//! Radiation is zero until the first ck call.
//!   Lir_top   net thermal (longwave) luminosity out through the top face ie+1:
//!             sum F_net(ie+1) A(ie+1), F_net = up - down, what the update differences
//!   Lsw_refl  stellar power reflected: albedo/(1-albedo) Lsw_abs, per column.  The
//!             deposit is (1-albedo) times the local beam absorption, so the power the
//!             model reflects is albedo times the power it INTERCEPTS, on the same
//!             effective area (pi r_abs^2, r_abs where the slant depth reaches ~1), NOT
//!             albedo F_star pi r_top^2 (the pre-dhj-wall-0927 definition, which counted
//!             light that never reaches the absorbing layer).
//!   Lsw_abs   stellar power absorbed (the net stellar heating): sum Q_sw V over cells
//!   Lrad_bot  net radiative luminosity up through the deepest ck face (icut; the inner
//!             wall when ck_pcut_bar covers the column); under ck_int_at_cut (route B)
//!             this is where sigma T_int^4 enters
//!   Etot_top  total energy flux out through the top: fluid + Lir_top - Lsw_abs
//!   Etot_bot  total energy flux up through the inner wall: fluid + Lrad_bot.  Lrad_bot
//!             is the flux the ck kernel APPLIES at the cut face (the deposit telescopes
//!             to it), so the column means the same for every ck_sph_face.  Under route
//!             B the faces differ in what the wall delivers: face 0 carries pi I_int =
//!             sigma T_int^4 exactly; faces 5 and 6 (Marshak-type wall J + beta H =
//!             (B + I_int)(1/2 + beta/4) with the pair at mu_eff) carry
//!             2 pi mu (1 + beta/2)/(1 + beta mu) I_int, = 4/(2 + sqrt 3) = 1.0718
//!             sigma T_int^4 in the isotropic deep limit (beta = 2, mu = 1/sqrt 3):
//!             a real 7.2 % larger internal heating, not a diagnostic offset
//!             (f6default_1002/RESULTS.md).
//!   Mdot_top  mass flux out through the top face,  Mdot_bot  up through the inner wall
//!   Lsw_in    stellar power incident on the DOMAIN TOP sphere, F_star sum max(mu0,0)
//!             A(ie+1) = F_star pi r_top^2.  Starlight bookkeeping closes as
//!             Lsw_in = Lsw_refl + Lsw_abs + Lsw_trans, where the remainder Lsw_trans is
//!             the light crossing the optically thin shell r_abs < b < r_top (and the
//!             pseudo-spherical twilight) without being absorbed: it leaves the domain.
//!             (The ck beam is a LOCAL absorption rate, by design -- see ck_beam_sph in
//!             two_stream_rt.hpp; its photon budget is exact, tests_ck_sph 916dc953.)
//! With problem/flux_hst_wall = true four more columns follow Lsw_in: Mdot_bot_in,
//! Mdot_bot_out (sums over inner-wall faces with fx1(IDN,is) > 0 / < 0) and Efl_bot_in,
//! Efl_bot_out (the fluid energy flux fx1(IEN,is) A over the same faces).

//----------------------------------------------------------------------------------------
//! \brief problem/am_hst: the eight angular-momentum columns from column n0 on.  Lrel and
//! Iw are the state at the row (u0): Lrel = sum vol (x P_y - y P_x), P = rho V rebuilt
//! from the covariant momenta, Iw = sum vol rho (x^2 + y^2), both with the cell-centre
//! position the Coriolis source uses, so L_abs = Lrel + Omega Iw is the quantity
//! problem/coriolis_am conserves.  Tcor (the torque the source applied: Tcont with
//! coriolis_am, else Tcc) and the AmHstStage sums are rates averaged over the interval
//! since the previous row [g cm^2 s^-2].  Then
//!     d(L_abs)/dt = Tcor + Omega dI/dt + (all other torques),
//!     Omega dI/dt = Tbnd - Tcont  (+ floors),
//! so without the fix Tcc - Tcont + Tbnd is the Coriolis part of d(L_abs)/dt.

void AmHstColumns(HistoryData *pdata, Mesh *pm, const int n0) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v = pmbp->pcoord->x2v;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;
  const bool cs = pm->use_cubed_sphere;
  const int ncell = nmb*nx3*nx2*nx1;
  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("dhj_am_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/(nx3*nx2*nx1);
    int rem = idx - m*(nx3*nx2*nx1);
    const int k = ks + rem/(nx2*nx1);
    rem -= (k - ks)*(nx2*nx1);
    const int j = js + rem/nx1;
    const int i = is + rem - (j - js)*nx1;
    const Real r = x1v(m,i);
    const Real dv = vol(m,k,j,i);
    Real lz, wc;
    if (cs) {
      const Real xi_c = 0.25*M_PI*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                              size.d_view(m).x2max);
      const Real eta_c = 0.25*M_PI*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                               size.d_view(m).x3max);
      const int pnl = mbpanel.d_view(m);
      Real qc[3], e1[3], e2[3];
      cubed_sphere::PanelToCart(pnl, xi_c, eta_c, qc);
      cubed_sphere::PanelTangents(pnl, xi_c, eta_c, e1, e2);
      const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
      const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn;
      const Real xx = r*rh0, yy = r*rh1;
      // covariant m2, m3 on the unit pair e1, e2 (e1.e2 = c): raise, then rebuild P
      const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
      const Real det = 1.0 - c*c;
      const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
      const Real p2 = (m2 - c*m3)/det, p3 = (m3 - c*m2)/det;
      const Real px = m1*rh0 + p2*e1[0] + p3*e2[0];
      const Real py = m1*rh1 + p2*e1[1] + p3*e2[1];
      lz = xx*py - yy*px;
      wc = xx*xx + yy*yy;
    } else {
      const Real rc = r*sin(x2v(m,j));
      lz = rc*u0(m,IM3,k,j,i);
      wc = rc*rc;
    }
    array_sum::GlobalSum h;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
    h.the_array[0] = lz*dv;
    h.the_array[1] = u0(m,IDN,k,j,i)*wc*dv;
    mb_sum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  Kokkos::fence();
  const Real tint = pm->time - amh::t0;
  Real rate[amh::NACC];
  for (int n=0; n<amh::NACC; ++n) {
    rate[n] = (tint > 0.0) ? amh::acc[n]/tint : 0.0;
    amh::acc[n] = 0.0;
  }
  amh::t0 = pm->time;
  amh::tdrag_rate = rate[11];
  pdata->hdata[n0] = sum.the_array[0];
  pdata->hdata[n0+1] = sum.the_array[1];
  pdata->hdata[n0+2] = amh::fix ? rate[1] : rate[0];
  for (int n=0; n<3; ++n) pdata->hdata[n0+3+n] = rate[n];
  if (amh::split) {
    for (int n=0; n<6; ++n) pdata->hdata[n0+6+n] = rate[5+n];
  } else {
    for (int n=0; n<2; ++n) pdata->hdata[n0+6+n] = rate[3+n];
  }
}

//----------------------------------------------------------------------------------------
//! \fn DhjTideTorque (tidetrq-1101)
//! \brief the z torque [g cm^2 s^-2] of the APPLIED tangential tide about the spin axis,
//! sum over this rank's active cells of rho (r x a_t)_z vol = rho R (a_t . phi_hat) vol
//! (the history output sums the ranks), at the current w0 (instantaneous, like Lrel/Iw):
//! the explicit problem/roche tangential -grad TidePot, or with problem/roche_wb_tan the
//! probed g_t (0 before the probe).  The radial part carries no z torque.  History
//! column Ttide, present only with problem/roche, to close L_z in the planet runs.

Real DhjTideTorque(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool mhd = (pmbp->pmhd != nullptr);
  auto w0 = mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v_ = pmbp->pcoord->x2v;
  auto x3v_ = pmbp->pcoord->x3v;
  auto &size = pmbp->pmb->mb_size;
  auto &panel = pmbp->pmb->mb_panel;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, js = ic.js, ks = ic.ks;
  const int nx1 = ic.nx1, nx2 = ic.nx2, nx3 = ic.nx3;
  const int nmb = pmbp->nmb_thispack;
  const bool cs = pm->use_cubed_sphere, sp = pm->use_spherical_polar;
  const Real rgms = pm->pgen->hot_jupiter_param.roche_gms;
  const Real rra = pm->pgen->hot_jupiter_param.roche_a;
  const bool wbt = dhj_wbt.on;
  if ((wbt && !dhj_wbt.ready) || !(cs || sp)) return 0.0;
  DvceArray5D<Real> g;
  if (wbt) g = *dhj_wbt.gt;
  Real tq = 0.0;
  Kokkos::parallel_reduce("dhj_tide_torque",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nx3*nx2*nx1),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/(nx3*nx2*nx1);
    int rem = idx - m*(nx3*nx2*nx1);
    const int k = ks + rem/(nx2*nx1);
    rem -= (k - ks)*(nx2*nx1);
    const int j = js + rem/nx1;
    const int i = is + rem - (j - js)*nx1;
    const Real r = x1v(m,i);
    const Real rv = w0(m,IDN,k,j,i)*vol(m,k,j,i);
    if (sp) {
      const Real theta = x2v_(m,j), phi = x3v_(m,k) - M_PI;
      const Real sine = sin(theta), cosine = cos(theta);
      Real aphi;
      if (wbt) {
        aphi = g(m,1,k,j,i);
      } else {
        const Real cph = cos(phi), sph = sin(phi);
        Real gs, gy, gz;
        roche::TideGrad(rgms, rra, r*sine*cph, r*sine*sph, r*cosine, gs, gy, gz);
        aphi = -(-gs*sph + gy*cph);
      }
      sum += rv*r*sine*aphi;
    } else {
      const Real x2 = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      const Real x3 = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
      const int pnl = panel.d_view(m);
      Real qc[3], e1[3], e2[3];
      cubed_sphere::PanelToCart(pnl, 0.25*M_PI*x2, 0.25*M_PI*x3, qc);
      cubed_sphere::PanelTangents(pnl, 0.25*M_PI*x2, 0.25*M_PI*x3, e1, e2);
      const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
      const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn, rh2 = qc[2]*qn;
      Real ax, ay;
      if (wbt) {
        // covariant (g2, g3) on the unit pair (e1, e2) -> the tangential vector
        const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
        const Real id = 1.0/(1.0 - c*c);
        const Real a1 = (g(m,0,k,j,i) - c*g(m,1,k,j,i))*id;
        const Real a2 = (g(m,1,k,j,i) - c*g(m,0,k,j,i))*id;
        ax = a1*e1[0] + a2*e2[0];
        ay = a1*e1[1] + a2*e2[1];
      } else {
        // as the source: Cartesian acceleration (gs, gy, -gz) in the chart frame
        Real gs, gy, gz;
        roche::TideGrad(rgms, rra, -r*rh0, -r*rh1, r*rh2, gs, gy, gz);
        ax = gs;
        ay = gy;
      }
      sum += rv*r*(rh0*ay - rh1*ax);
    }
  }, Kokkos::Sum<Real>(tq));
  return tq;
}

void DhjFluxHistory(HistoryData *pdata, Mesh *pm) {
  namespace ts = two_stream_rt;
  constexpr int NFH0 = 9;
  const bool fwall = dhj_flux_hst_wall;
  const bool ffloor = dhj_flux_hst_floor;
  // problem/am_hst: eight columns after all the others
  const int NFL = (fwall ? (NFH0 + 4) : NFH0) + (ffloor ? 3 : 0);
  const int NFH = NFL + (amh::on ? (amh::split ? 12 : 8) : 0);
  // tidetrq-1101: Ttide (DhjTideTorque) after all the others, only with problem/roche
  const bool ttide = pm->pgen->hot_jupiter_param.roche;
  // floortrc-1101: with <mhd>/floor_tracer and problem/flux_hst_floor, Mdrain [g/s] and
  // Edrain [erg/s] after all the others: the floor-tracer drain (MHD::FloorTracerDrain)
  // over the same interval as Mfloor/Efloor, which are then NET (added - drained); the
  // gross floor addition is Mfloor + Mdrain
  // relaxdrag-1101: Tdrag (the start-up Rayleigh drag's z torque, a rate over the
  // interval like Tcor) after Ttide, only with problem/am_hst
  const int NTD = NFH + (ttide ? 1 : 0);
  const int NDR = NTD + (amh::on ? 1 : 0);   // Mdrain/Edrain after Ttide and Tdrag
  // (the floor tracer is a default with the Alfven floor: when the two columns do not
  // fit, they are left out with a note instead of stopping the run)
  bool fdrain = ffloor && (pm->pmb_pack->pmhd != nullptr) &&
                pm->pmb_pack->pmhd->floor_tracer;
  if (fdrain && NDR + 2 > NHISTORY_VARIABLES) {
    static bool warned = false;
    if (!warned && global_variable::my_rank == 0) {
      std::cout << "### NOTE: dhj flux history: no room for Mdrain/Edrain ("
                << NDR + 2 << " > NHISTORY_VARIABLES = " << NHISTORY_VARIABLES
                << "; problem/am_hst_split = false frees 4); Mfloor/Efloor are NET of "
                << "the floor-tracer drain" << std::endl;
    }
    warned = true;
    fdrain = false;
  }
  pdata->nhist = NDR + (fdrain ? 2 : 0);
  if (ttide) pdata->label[NFH] = "Ttide";
  if (amh::on) pdata->label[NTD] = "Tdrag";
  if (fdrain) {
    pdata->label[NDR] = "Mdrain";
    pdata->label[NDR+1] = "Edrain";
  }
  pdata->label[0] = "Lir_top";
  pdata->label[1] = "Lsw_refl";
  pdata->label[2] = "Lsw_abs";
  pdata->label[3] = "Lrad_bot";
  pdata->label[4] = "Etot_top";
  pdata->label[5] = "Etot_bot";
  pdata->label[6] = "Mdot_top";
  pdata->label[7] = "Mdot_bot";
  pdata->label[8] = "Lsw_in";
  if (fwall) {
    pdata->label[9] = "Mdot_bot_in";
    pdata->label[10] = "Mdot_bot_out";
    pdata->label[11] = "Efl_bot_in";
    pdata->label[12] = "Efl_bot_out";
  }
  if (ffloor) {
    pdata->label[NFL-3] = "Efloor";
    pdata->label[NFL-2] = "Mfloor";
    pdata->label[NFL-1] = "Efloor_rt";
  }
  if (amh::on) {
    pdata->label[NFL] = "Lrel";
    pdata->label[NFL+1] = "Iw";
    pdata->label[NFL+2] = "Tcor";
    pdata->label[NFL+3] = "Tcc";
    pdata->label[NFL+4] = "Tcont";
    pdata->label[NFL+5] = "Tbnd";
    if (amh::split) {
      pdata->label[NFL+6] = "Tx1";
      pdata->label[NFL+7] = "Tx2";
      pdata->label[NFL+8] = "Tx3";
      pdata->label[NFL+9] = "Tsrc";
      pdata->label[NFL+10] = "Tcons";
      pdata->label[NFL+11] = "Twall";
    } else {
      pdata->label[NFL+6] = "Tdif_deep";
      pdata->label[NFL+7] = "Tdif_rad";
    }
  }
  for (int n=0; n<NFH; ++n) pdata->hdata[n] = 0.0;

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int ncol = nmb*nx3*nx2;
  DvceFaceFld5D<Real> uflx = (pmbp->pmhd != nullptr) ? pmbp->pmhd->uflx
                                                     : pmbp->phydro->uflx;
  auto fx1 = uflx.x1f;
  auto area = pmbp->pcoord->area.x1f;
  auto vol = pmbp->pcoord->volume;
  bool rad = ts::rt_face_flux_ready() && (ts::rt_Qb_ptr != nullptr)
             && (ts::rt_cf_ptr != nullptr) && (ts::rt_icut_ptr != nullptr);
  int nblk = rad ? ts::rt_face_nblk() : 0;
  DvceArray5D<Real> fb, qb;
  DvceArray4D<Real> cf;
  DvceArray3D<int> icut;
  Real fstar = ts::rt_hist_fstar;
  Real alb = ts::rt_hist_albedo;
  if (rad) {
    fb = ts::rt_face_flux();
    qb = *ts::rt_Qb_ptr;
    cf = *ts::rt_cf_ptr;
    icut = ts::rt_cut_index();
  } else if (ts::RtFhStageViews(nmb, fb, qb, cf, icut, fstar, alb)) {
    // problem/flux_hst_rst: a restart, before its first RT call: the cache the
    // restart file carries (chain blocks pre-summed, so nblk = 1 gives the same sums)
    rad = true;
    nblk = 1;
  }
  // problem/flux_hst_rkavg: the stage-weighted fluid fluxes once a whole cycle is in
  // the accumulator (else, e.g. the first row after a restart, the last stage's)
  const bool fav = dhj_fhst_avg && (dhj_fhst_stage == dhj_fhst_nst);
  auto fav_ = fav ? *dhj_fhst_ptr : DvceArray4D<Real>();
  // stswall-1029: under RKG super-stepping the resistive (Poynting) energy flux is not in
  // uflx during the RK stages, so neither column above contains it.  The flux the last
  // super-step applied through ie+1 and is (Resistivity::rkg_bflx, slots 5 and 4) is
  // FOLDED INTO Etot_top and Etot_bot rather than given its own column: the history is
  // at its NREDUCTION_VARIABLES cap.  Explicit resistivity adds it to uflx in the RK
  // stages, so the columns already hold it there.  0 at the conducting wall (E_t = 0).
  const bool rkb = (pmbp->pmhd != nullptr) && (pmbp->pmhd->presist != nullptr) &&
                   pmbp->pmhd->presist->track_x1_bflx &&
                   (pmbp->pmhd->presist->rkg_bflx.extent_int(0) == pmbp->nmb_thispack);
  auto rkb_ = rkb ? pmbp->pmhd->presist->rkg_bflx : DvceArray4D<Real>();
  const Real arat = (alb < 1.0) ? alb/(1.0 - alb) : 0.0;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("dhj_flux_hist",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int k = ks + kj/nx2;
    const int j = js + kj - (kj/nx2)*nx2;
    const Real at = area(m,k,j,ie+1);
    const Real ab = area(m,k,j,is);
    Real lir = 0.0, lrf = 0.0, lsw = 0.0, lbot = 0.0, lin = 0.0;
    if (rad) {
      const int ic = icut(m,k,j);
      if (ic <= ie) {
        Real ft = 0.0, fc = 0.0;
        for (int b=0; b<nblk; ++b) {
          ft += fb(m,b,ie+1,k,j);
          fc += fb(m,b,ic,k,j);
        }
        lir = ft*at;
        lbot = fc*area(m,k,j,ic);
        for (int i=ic; i<=ie; ++i) {
          Real q = 0.0;
          for (int b=0; b<nblk; ++b) q += qb(m,b,i,k,j);
          lsw += q*vol(m,k,j,i);
        }
      }
      const Real mu0 = cf(m,k,j,3);
      lin = (mu0 > 0.0) ? fstar*mu0*at : 0.0;
      lrf = arat*lsw;
    }
    array_sum::GlobalSum hvars;
    hvars.the_array[0] = lir;
    hvars.the_array[1] = lrf;
    hvars.the_array[2] = lsw;
    hvars.the_array[3] = lbot;
    const Real fet = (fav ? fav_(m,0,k,j) : fx1(m,IEN,k,j,ie+1))
                     + (rkb ? rkb_(m,5,k,j) : 0.0);
    const Real fmt = fav ? fav_(m,1,k,j) : fx1(m,IDN,k,j,ie+1);
    const Real feb = (fav ? fav_(m,2,k,j) : fx1(m,IEN,k,j,is))
                     + (rkb ? rkb_(m,4,k,j) : 0.0);
    const Real fmb = fav ? fav_(m,3,k,j) : fx1(m,IDN,k,j,is);
    hvars.the_array[4] = fet*at + lir - lsw;
    hvars.the_array[5] = feb*ab + lbot;
    hvars.the_array[6] = fmt*at;
    hvars.the_array[7] = fmb*ab;
    hvars.the_array[8] = lin;
    for (int n=NFH0; n<NHISTORY_VARIABLES; ++n) hvars.the_array[n] = 0.0;
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();
  for (int n=0; n<NFH0; ++n) pdata->hdata[n] = sum_this_mb.the_array[n];
  if (ffloor) {
    // problem/flux_hst_floor: the ConToPrim floor energy and mass of this rank since the
    // previous row (the history output sums the ranks), as a rate; then restart the sum
    DvceArray2D<Real> acc = (pmbp->pmhd != nullptr) ? pmbp->pmhd->c2p_acc
                                                    : pmbp->phydro->c2p_acc;
    // problem/budget_dt: the ledger drains c2p_acc; the columns read its pending copy
    if (ebud::on) {
      EBudDrainFloors(pm);
      acc = *ebud::fhst;
    }
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), acc);
    Real ef = 0.0, mf = 0.0, efr = 0.0;
    for (int i=0; i<static_cast<int>(h.extent(1)); ++i) {
      ef += h(0,i) + h(2,i);
      mf += h(1,i) + h(3,i);
      efr += h(2,i);
    }
    Kokkos::deep_copy(acc, 0.0);
    const Real tint = pm->time - dhj_floor_t0;
    dhj_floor_t0 = pm->time;
    pdata->hdata[NFL-3] = (tint > 0.0) ? ef/tint : 0.0;
    pdata->hdata[NFL-2] = (tint > 0.0) ? mf/tint : 0.0;
    pdata->hdata[NFL-1] = (tint > 0.0) ? efr/tint : 0.0;
    if (fdrain) {
      auto *pmh = pmbp->pmhd;
      pdata->hdata[NDR] = (tint > 0.0) ? pmh->fltr_dm/tint : 0.0;
      pdata->hdata[NDR+1] = (tint > 0.0) ? pmh->fltr_de/tint : 0.0;
      pmh->fltr_dm = 0.0;
      pmh->fltr_de = 0.0;
    }
  }
  if (amh::on) AmHstColumns(pdata, pm, NFL);
  if (ttide) pdata->hdata[NFH] = DhjTideTorque(pm);
  if (amh::on) pdata->hdata[NTD] = amh::tdrag_rate;
  if (!fwall) return;
  // problem/flux_hst_wall: the inner-wall in/out split, its own small reduction
  array_sum::GlobalSum wsum;
  Kokkos::parallel_reduce("dhj_flux_hist_wall",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncol),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &mb_sum) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int k = ks + kj/nx2;
    const int j = js + kj - (kj/nx2)*nx2;
    const Real ab = area(m,k,j,is);
    const Real mb = fx1(m,IDN,k,j,is)*ab;
    const Real eb = fx1(m,IEN,k,j,is)*ab;
    array_sum::GlobalSum hvars;
    for (int n=0; n<NHISTORY_VARIABLES; ++n) hvars.the_array[n] = 0.0;
    if (mb > 0.0) {
      hvars.the_array[0] = mb;  hvars.the_array[2] = eb;
    } else {
      hvars.the_array[1] = mb;  hvars.the_array[3] = eb;
    }
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(wsum));
  Kokkos::fence();
  for (int n=0; n<4; ++n) pdata->hdata[NFH0+n] = wsum.the_array[n];
}
