#ifndef HYDRO_DUAL_ENERGY_HPP_
#define HYDRO_DUAL_ENERGY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dual_energy.hpp
//! \brief duale-1112: dual-energy selection of the internal energy, the per-cell audit of
//! the energy floors, and the payment of the energy they create from the cell's own
//! radiation (hydro and MHD, every geometry; general EOS only).  Every key is read only
//! when named; absent, nothing here is allocated or called (bitwise).
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
//!  floor_pay = rad: the energy the floors (and the dual selection) create in an active
//!    cell, e_final - (E - KE), is taken from the cell's E_rad after the radiation stage
//!    (RadiationM1::HydroConToPrim): gas + (c/chat) E_rad is exact per cell; what E_rad
//!    cannot cover (down to e_floor) is counted as unpaid.
//!  floor_audit = N: report the cumulative per-cell counters every N cycles (rank 0):
//!    cells with E - KE <= 0, KE-dominated, below e(rho, tfloor), below the EOS table
//!    edge, dual selected, floored (metric-correct on the cubed sphere), and the volume
//!    integrated energies.  Diagnostic only.
//!  fofc_retest = true (hydro): FOFC re-tests every flagged trial cell with the
//!    metric-correct internal energy (and the dual selection) and unflags those no floor
//!    would touch (the cubed-sphere ConsToPrim test uses an ORTHONORMAL kinetic energy).

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "eos/eos.hpp"

class MeshBlockPack;

namespace fluid {

// counter slots (cumulative; energies volume integrated, code units)
enum DualCnt {
  DC_N = 0, DC_KEDOM = 1, DC_LE0 = 2, DC_SEL = 3, DC_SEL_DE = 4, DC_FLR = 5,
  DC_FLR_DE = 6, DC_SUBT = 7, DC_SUBEDGE = 8, DC_SKIP = 9, DC_SUBT_SEL = 10,
  DC_FLR_KED = 11,
  DC_REFRESH = 12,   // the same 12 slots for a ConToPrim that follows no hydro update
  DC_PAID = 24, DC_UNPAID = 25, DC_NPAY = 26, DC_RETEST = 27, DC_UNFLAG = 28,
  DC_GIVEN = 29, DC_TOPE_U = 30, DC_TOPE_R = 31, DC_NCNT = 32
};

struct DualEnergy {
  bool on = false;      // <blk>/dual_energy
  bool pay = false;     // <blk>/floor_pay = rad
  bool retest = false;  // <blk>/fofc_retest
  bool any = false;     // Pre/Post bookkeeping active
  bool armed = false;   // a hydro update evolved the slots since the last Post
  bool valid = false;   // the slots hold a resynced state
  bool init = false;    // arrays allocated
  Real eta = 1.0e-3;
  int report = 0;       // <blk>/floor_audit
  int last_report = -1;
  int iden = -1, ivol = -1;
  std::string blk;
  DvceArray5D<Real> pre;    // (nmb, 3, n3, n2, n1): E, d, |m|^2 before the ConToPrim
  DvceArray4D<Real> pst;    // pressure of the last resynced state
  DvceArray4D<Real> payv;   // energy created in active cells, not yet paid (per volume)
  Real hc[DC_NCNT] = {};    // this rank's cumulative counters

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
};

//! floor_pay = rad: pay `payv` of the fluid from E_rad (slot ie of urad) of the same
//! cell (every cell, ghosts mirror their neighbour); ch_cl = chat/c; efl = e_floor.
//! Called ONCE PER STEP by the Driver after every radiation stage (so no RK / time2
//! register of the radiation re-weights or drops the payment); resets payv.
void DualPayFromRad(MeshBlockPack *pp, DualEnergy &de, DvceArray5D<Real> &urad,
                    const int ie, const Real ch_cl, const Real efl);

}  // namespace fluid

#endif  // HYDRO_DUAL_ENERGY_HPP_
