#ifndef HYDRO_DUAL_ENERGY_HPP_
#define HYDRO_DUAL_ENERGY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dual_energy.hpp
//! \brief dual-energy selection of the internal energy, the per-cell audit of
//! the energy floors, and the payment of the energy they create from the cell's own
//! radiation (hydro and MHD, every geometry; general EOS only).  Every key but floor_pay
//! is read only when named; absent, nothing here is allocated or called (bitwise).
//!
//! Keys (<hydro> or <mhd>):
//!  dual_energy = true: two scalar slots after every other scalar, the internal energy
//!    e (specific value e/rho) and the volume tracer 1 (specific value 1/rho).  Both are
//!    advected by the mass flux of the same Riemann solver as any passive scalar, so
//!    after a stage the volume slot is 1 - beta dt div(v) in the exact discrete (metric)
//!    form of the grid.  The p dV work is added at the ConToPrim:
//!    e_d = e_s exp(p_n (V_s - 1)/e_s), p_n the pressure of the stage-start state (2nd
//!    order in the RK stage, positive by construction).  e_d replaces E - KE (- EB) only
//!    where E - KE <= 0 or < dual_energy_eta E (KE- or field-dominated).  The TOTAL
//!    energy stays the conserved variable: E is rebuilt as e_d + KE (+ EB) and the change
//!    is counted (and paid from the radiation under floor_pay = rad).  Every ConToPrim
//!    resyncs the slots to the final state (e_s = e, V_s = 1), so the radiation exchange,
//!    source terms and floors are all seen; only a ConToPrim that follows a hydro update
//!    (armed by Fluxes) selects.
//!  dual_energy_eta = 1e-3: the selection threshold e/E.
//!  floor_pay = rad_column (default of a fresh run with eos = general and <rad_vet>,
//!    recorded; none otherwise, and the only value without them) | rad | none.
//!  floor_pay = rad: the energy the floors (and the dual selection) create in an active
//!    cell, e_final - (E - KE), is taken from the cell's E_rad after the radiation stage
//!    (RadiationVet::GasUpdate): gas + (c/chat) E_rad is exact per cell; what E_rad
//!    cannot cover (down to e_floor) is counted as unpaid.  floor_pay = rad_column:
//!    that rest is taken from the E_rad of the cell's radial column (MeshBlock, nearest
//!    first, as implicit_pos_gas_cons); what is still left is created, CLEAN WARNING.
//!  floor_pay_skip = true (needs floor_pay = rad*): also the cells whose density or
//!    momentum the inversion changed (dfloor, vceil): their total energy change
//!    E_after - E_before (kinetic energy removed, the floors) is given to / taken from
//!    E_rad the same way, so gas + radiation is exact there too.  Off (default), it is
//!    only counted: floor_audit skip_de (+ under etotgrav the floor mass's
//!    (d_after - d_before) Phi, as Etot counts it; never paid, it depends on the zero of
//!    Phi) and skip_dm, the run totals at the end.
//!  floor_audit = N: report the cumulative per-cell counters every N cycles (rank 0):
//!    cells with E - KE <= 0, KE-dominated, below e(rho, tfloor), below the EOS table
//!    edge, dual selected, floored (metric-correct on the cubed sphere), and the volume
//!    integrated energies.  Diagnostic only.
//!  floor_audit_worst = true (needs floor_audit): with each report, every rank prints
//!    its active cell of the largest (E_after - E_before) x vol of the inversion, with
//!    the cell's state before and after.  Diagnostic only.
//!  Momentum ledger (hydro with floor_from_ekin, or vceil with the bookkeeping above on):
//!    the momentum the accepted step's inversions took from the active cells, each with
//!    its end-of-step (RK) weight, by cause: ekin (floor_from_ekin), vceil, other (a
//!    density floor, or both with a raised E).  Components (r, theta, phi), Cartesian
//!    (x, y, z); one CLEAN FLOORMOM line per step in which any fired: cells, signed sums
//!    and sums of |dp| x vol.  Exact sums (repro_sum.hpp): any MeshBlock layout and rank
//!    count.  Per step and per job, as floor_pay's counts (not in the restart file).
//!  fofc_retest = true (hydro): FOFC re-tests every flagged trial cell with the
//!    metric-correct internal energy (and the dual selection) and unflags those no floor
//!    would touch (the cubed-sphere ConsToPrim test uses an ORTHONORMAL kinetic energy).

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "eos/eos.hpp"
#include "utils/repro_sum.hpp"

class MeshBlockPack;

namespace fluid {

// counter slots (cumulative; energies volume integrated, code units)
enum DualCnt {
  DC_N = 0, DC_KEDOM = 1, DC_LE0 = 2, DC_SEL = 3, DC_SEL_DE = 4, DC_FLR = 5,
  DC_FLR_DE = 6, DC_SUBT = 7, DC_SUBEDGE = 8, DC_SKIP = 9, DC_SUBT_SEL = 10,
  DC_FLR_KED = 11,
  DC_REFRESH = 12,   // the same 12 slots for a ConToPrim that follows no hydro update
  DC_PAID = 24, DC_UNPAID = 25, DC_NPAY = 26, DC_RETEST = 27, DC_UNFLAG = 28,
  DC_GIVEN = 29, DC_TOPE_U = 30, DC_TOPE_R = 31,
  DC_FLR_R = 32, DC_FLR_RDE = 33,   // (+ DC_REFRESH = 44, 45: refresh) sum r, sum de r
  DC_SKIP_DE = 34, DC_SKIP_DM = 35,  // (+ 46, 47) the skip cells' E and mass change x vol
  DC_NCNT = 48
};

struct DualEnergy {
  bool on = false;      // <blk>/dual_energy
  bool pay = false;     // <blk>/floor_pay = rad (or rad_column)
  bool pcol = false;    // <blk>/floor_pay = rad_column: the rest from the radial column
  bool pretry = false;  // <blk>/floor_pay_retry: unpaid energy fails the cycle (dt_retry)
  bool pskip = false;   // <blk>/floor_pay_skip: the density floor / velocity ceiling too
  bool retest = false;  // <blk>/fofc_retest
  bool worst = false;   // <blk>/floor_audit_worst
  bool any = false;     // Pre/Post bookkeeping active
  bool armed = false;   // a hydro update evolved the slots since the last Post
  bool valid = false;   // the slots hold a resynced state
  bool init = false;    // arrays allocated
  bool led = false;     // the momentum ledger
  bool book = true;     // false: a conversion whose state is put back (Time2GasTime)
  Real eta = 1.0e-3;
  Real wst = 1.0;       // end-of-step weight of a change made after the current stage
  int report = 0;       // <blk>/floor_audit
  int last_report = -1;
  int iden = -1, ivol = -1;
  std::string blk;
  DvceArray5D<Real> pre;    // (nmb, 3, n3, n2, n1): E, d, |m|^2 before the ConToPrim
  DvceArray4D<Real> pst;    // pressure of the last resynced state
  DvceArray4D<Real> payv;   // energy created in active cells, not yet paid (per volume)
  DvceArray4D<Real> phi;    // <hydro>/etotgrav: the potential, rho*Phi of a floor's mass
  Real hc[DC_NCNT] = {};    // this rank's cumulative counters
  Real hck[DC_NCNT] = {};   // hc at the step's start (Mark; Discard puts it back)
  // the ledger [ekin, vceil, other][cells, dp_1-3, |dp_1-3|]: step, its start, the job
  static constexpr int nlq = 7;
  repro::ReproSum lsum[3*nlq], lsk[3*nlq];
  Real ltot[3*nlq] = {};

  DualEnergy() {
    for (int q = 0; q < 3*nlq; ++q) {
      lsum[q].Init();
      lsk[q].Init();
    }
  }

  // read the keys; `nfluid` = nhydro or nmhd, `du_off` = tracers.DuOff()
  void Read(ParameterInput *pin, const std::string &block, const EOS_Data &eos,
            int nfluid, int du_off);
  void Allocate(MeshBlockPack *pp, const DvceArray5D<Real> &u0);
  // before the ConToPrim (after RemoveGravEtot): E, d, |m|^2 of every cell
  void Pre(MeshBlockPack *pp, const DvceArray5D<Real> &u0);
  // after ConToPrim (+ the cubed-sphere RaiseVel), before AddGravEtot
  void Post(MeshBlockPack *pp, DvceArray5D<Real> &u0, DvceArray5D<Real> &w0,
            DvceArray5D<Real> &wder, DvceArray4D<Real> &wtemp, const EOS_Data &eos);
  void Report(MeshBlockPack *pp);
  void Worst(MeshBlockPack *pp, const DvceArray5D<Real> &u0, const DvceArray5D<Real> &w0);
  void Book(MeshBlockPack *pp, const DvceArray5D<Real> &u0, const EOS_Data &eos);
  // the step's line (and zero) / the job's totals
  void LedgerLine(const std::string &when);
  void LedgerTotal();
  // a discarded attempt (dt_retry, the sdirk2 redo): its floors' payment and counts go
  // (Mark, at the step's start before the checkpoint, reports first: a report closes the
  // steps before, and a retried attempt, which does not mark, reports nothing again)
  void Mark(MeshBlockPack *pp) {
    Report(pp);
    for (int n = 0; n < DC_NCNT; ++n) hck[n] = hc[n];
    for (int q = 0; q < 3*nlq; ++q) lsk[q] = lsum[q];
  }
  void Discard() {
    for (int n = 0; n < DC_NCNT; ++n) hc[n] = hck[n];
    for (int q = 0; q < 3*nlq; ++q) lsum[q] = lsk[q];
    if (pay && init) Kokkos::deep_copy(payv, 0.0);
  }
};

//! floor_pay = rad: pay `payv` of the fluid from E_rad (slot ie of urad) of the same
//! cell (every cell, ghosts mirror their neighbour); ch_cl = chat/c; efl = e_floor.
//! Called ONCE PER STEP by the Driver after every radiation stage (so no RK / time2
//! register of the radiation re-weights or drops the payment); resets payv.  Returns
//! the energy left unpaid (all ranks; 0 unless rad_column or floor_pay_retry).
Real DualPayFromRad(MeshBlockPack *pp, DualEnergy &de, DvceArray5D<Real> &urad,
                    const int ie, const Real ch_cl, const Real efl);

}  // namespace fluid

#endif  // HYDRO_DUAL_ENERGY_HPP_
