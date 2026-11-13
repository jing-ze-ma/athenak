#ifndef PGEN_HE_DRIVE_HPP_
#define PGEN_HE_DRIVE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file he_drive.hpp
//! \brief hg262-drive-1113: the time-dependent DRIVEN inner boundary of he_star_m1
//! (he_drive.cpp).  Every key is read only when named; with none of them the problem
//! generator never calls into this file beyond the cheap On() queries (bitwise).
//!
//!  problem/he_bc_inner_table = <file>  (needs he_bc_inner = inflow): columns t, rho,
//!      T [K], v_r (code units); the inner ghosts hold the initial column's shape scaled
//!      to (rho, e(rho, T)) of the table at r_in, with v_r = v_tab (r_in/r)^2 rho_tab/rho
//!      (constant mass flux through the ghosts) and v_t = 0.
//!  problem/he_gm_live = true: the LIVE monopole gravity G m(r, t) = G M_int(t) +
//!      G int_{r_in}^{r} 4 pi r'^2 <rho> dr', applied as the INCREMENT over the potential
//!      in force when the reference was taken (the IC's on a fresh start, the restart
//!      state's when the keys are switched on at a restart): Phi = Phi_ref + dPhi, dPhi =
//!      the exact potential of dm(r) = dM_int + sum of the shell-mean mass changes
//!      (uniform density inside each cell).  Unchanged shell masses give dPhi = 0
//!      EXACTLY (the reference balance is kept bitwise).  Refreshed at the END of every
//!      problem/he_gm_live_every-th cycle (user_cycle_func, after the accepted step and
//!      its outputs); the energy keeps e_int: E += rho dPhi (the work of the
//!      time-dependent potential, booked as W_grav).
//!  problem/he_gm_table = <file> (optional, with he_gm_live): columns t, M_int (code mass
//!      inside r_in); G M_int(t) = he_gm + G (M_int(t) - M_int(t_ref)).
//!  problem/he_drive_t0 (default 0): table time = t - he_drive_t0 (both pgen tables).
//!  problem/he_drive_ledger (default: true when any drive key is named): the RK-exact
//!      ledgers of the mass and total energy through the x1 boundary faces, the sponge,
//!      W_grav and int F_in A dt, written with M_tot, E_tot every
//!      problem/he_drive_log_every cycles to <basename>.hsdrive.txt (rank 0).
//! All state travels in the restart file (trailer HSDRV001, the last one).

#include <string>
#include <vector>

#include "athena.hpp"

class Mesh;
class ParameterInput;

namespace hsdrive {

//! everything the problem generator knows that the drive needs
struct SetupArgs {
  Real gcode = 0.0;      // G in code units (0: no <units>)
  Real tk = 1.0;         // kelvin per code temperature unit
  Real gm = 0.0;         // he_gm = G m(r_in) of the column
  Real rho_rin = 0.0;    // the initial column at r_in: rho, eint (scales of the ghosts)
  Real eint_rin = 0.0;
  Real fconst = 0.0;     // <rad_m1>/implicit_flux_x1min (constant)
  bool binf = false;     // he_bc_inner = inflow
  bool pointmass = false;
  bool bc_hse = false;   // he_bc_outer = hse (reads the static G m(r))
};

bool On();               // any drive feature: the ledger and the cycle hook run
bool InflowOn();         // he_bc_inner_table
bool LiveOn();           // he_gm_live
//! read the keys (every start), after the pgen built its potentials
void Setup(ParameterInput *pin, Mesh *pm, bool restart, const SetupArgs &a);
//! strip the HSDRV001 trailer from the pgen restart block (call on every restart before
//! any other trailer is stripped; keeps the numbers for RestartApply)
void RestartStrip(std::vector<char> *blk);
//! at the end of a restart's UserProblem: re-apply the stored dPhi (no energy change:
//! the restart state carries it), take the reference if there is none, and do the
//! refresh the straight run did right after writing this restart
void RestartApply(Mesh *pm);
//! the inner-ghost table state at the BC time (InflowOn only; T in code units)
void InflowAt(Mesh *pm, Real *d, Real *t, Real *v);
//! once per hydro stage, from the user source after every u0 update (On only):
//! dsp = the stage's change of sum IEN V by the sponge
void StageLedger(Mesh *pm, Real bdt, bool binf, bool zflux, Real dsp);
//! user_cycle_func (On only)
void Cycle(Mesh *pm);
//! append the trailer (On only; called by the pgen restart writer, last)
void RstAppend(std::vector<char> *out);
//! the inner flux the closure / history read: the table flux at t when the rad_m1 flux
//! table is on, else fconst (bitwise the old constant)
Real InnerFlux(Mesh *pm, Real fconst);
void Final();

}  // namespace hsdrive

#endif  // PGEN_HE_DRIVE_HPP_
