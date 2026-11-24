//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_init.cpp
//! \brief set-up of the implicit radiation (<rad_m1>/transport = implicit): the key
//! reading and checks (ImplicitInit), the x1 partition (ImplicitPartitionInit), the
//! on-rank halo copies of the work array (ImplicitHalo*), the run-end report
//! (ImplicitReport, EventTotals) and the face-flux restart slabs (FaxRst*).

#include <float.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "utils/eint_ledger.hpp"
#include "rad_m1/rad_m1_t2reset.hpp"
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "utils/recipe.hpp"
#include "utils/repro_sum.hpp"
#include "rad_m1/rad_m1_layreduce.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "reconstruct/plm.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "utils/eint_from_cons.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
// reset-sector exemption of the hesdirk2 admissibility (rad_m1_t2reset.hpp)
bool g_t2_reset_on = false;
DvceArray4D<Real> g_t2_reset;
double g_t2_reset_nskip = 0.0;
static int M1HaloDepthEnv(MeshBlockPack *pp, ParameterInput *pin);   // (below)
// DEBUG <rad_m1>/dbg_eint_cells (file scope so the class
// layout is unchanged): 0 = off
static int g_dbg_eint_n = 0;

namespace {
//----------------------------------------------------------------------------------------
//! \fn ImplBCFromString
//! \brief one x1 boundary type of the implicit solve, from the input string, or from the
//! mesh boundary flag when the input does not name one.

int ImplBCFromString(const std::string &s, const BoundaryFlag mbc) {
  if (s.compare("auto") == 0) {
    if (mbc == BoundaryFlag::periodic) return M1_IBC_PERIODIC;
    if (mbc == BoundaryFlag::reflect) return M1_IBC_REFLECT;
    return M1_IBC_MARSHAK;
  }
  if (s.compare("marshak") == 0) return M1_IBC_MARSHAK;
  if (s.compare("flux") == 0) return M1_IBC_FLUX;
  if (s.compare("reflect") == 0) return M1_IBC_REFLECT;
  if (s.compare("periodic") == 0) return M1_IBC_PERIODIC;
  if (s.compare("efix") == 0) return M1_IBC_EFIX;
  if (s.compare("outflow") == 0) return M1_IBC_OUTFLOW;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    << std::endl << "<rad_m1>/implicit_bc_x1min|max = '" << s << "' is not a valid "
    << "choice (auto | marshak | flux | reflect | periodic | efix | outflow (scheme = "
    << "clean))" << std::endl;
  std::exit(EXIT_FAILURE);
  return M1_IBC_MARSHAK;
}

//----------------------------------------------------------------------------------------
//! \fn ImplFatal
//! \brief one fatal-error exit with a message, used by the option parsers

void ImplFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitInit
//! \brief read the <rad_m1> parameters of the implicit solver, check the restrictions of
//! 3a and allocate the face array and the work array.  Called at the END of the
//! constructor, from the delimited hook there, and a no-op in explicit mode.

void RadiationM1::ImplicitInit(ParameterInput *pin) {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;

  // ---- MILESTONE 3b phase B.  `implicit` = the full 7-point solve; `implicit_x1` keeps
  // every branch below on the 3a/3a2/3c arithmetic, bit for bit.
  const bool full = (transport == M1_TRANSPORT_IMPLICIT);
  // <rad_m1>/time_scheme (read by Time2Init below).  DEFAULT hesdirk2 since m1-defaults2
  // (tests_m1/runs_5g_defaults2; /viper/ptmp2/jinma/h2val_0924: hesdirk2 at cfl 0.9 is
  // 1.57x cheaper than be at cfl 0.3 at equal accuracy) wherever hesdirk2 is accepted:
  // transport = implicit (not implicit_x1) on a Cartesian mesh or (m1-sph2,
  // tests_m1/runs_5h_sph2) the spherical-polar wedge (no cubed sphere or polar
  // boundary), a closure other than tau (vet_col included since m1-sph2), and
  // <time>/integrator = rk2.  Elsewhere the default stays be, silently.  A restart whose
  // file lacks the key (written when it was read only when named) keeps be.  The
  // resolved value is echoed (so later restarts keep it); explicit input overrides.
  {auto *pmh = pmy_pack->pmesh;
  const std::string integ = pin->DoesParameterExist("time","integrator") ?
                            pin->GetString("time","integrator") : "rk2";
  // (USER: second order in time on cs too): the cubed sphere takes the same
  // hesdirk2 default as sp (was excluded; W1 measured hesdirk2 on cs at 1.95-1.99)
  const bool h2def = full && !pmh->use_polar_boundary &&
                     (!tau_closure || vet_col) && (integ.compare("rk2") == 0) &&
                     !global_variable::restart_run;
  (void) pin->GetOrAddString("rad_m1","time_scheme", h2def ? "hesdirk2" : "be");
  // plain be = the hesdirk2 backward-Euler path (t2_beonly, rad_m1.hpp)
  // wherever hesdirk2 is accepted (the h2def conditions without the restart rule);
  // <rad_m1>/be_legacy = true keeps the old step.  time2_be_only = true: the same path
  // under time_scheme = hesdirk2.  defaults-1030 (USER 10-30): this is the be default on
  // a fresh run; a restart whose file lacks be_legacy keeps the old step (be_legacy =
  // true, bitwise); the resolved value is echoed.
  {const std::string tsn = pin->GetString("rad_m1","time_scheme");
  const bool h2ok = full && !pmh->use_polar_boundary && (!tau_closure || vet_col) &&
                    (integ.compare("rk2") == 0);
  if (tsn.compare("be") == 0) {
    const bool leg = pin->GetOrAddBoolean("rad_m1","be_legacy",
                                          global_variable::restart_run);
    t2_beonly = h2ok && !leg;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> time_scheme = be: be_legacy = " << (leg ? "true" : "false")
                << (t2_beonly ? " (the hesdirk2 backward-Euler path)" : " (the old step)")
                << std::endl;
    }
  } else if (tsn.compare("hesdirk2") == 0 &&
             pin->DoesParameterExist("rad_m1","time2_be_only")) {
    t2_beonly = pin->GetBoolean("rad_m1","time2_be_only");
  }}}
  trans_on = false;
  trans_x3 = false;
  impl_cfl = pin->GetOrAddReal("rad_m1","implicit_cfl",-1.0);
  impl_tol = pin->GetOrAddReal("rad_m1","implicit_tol",1.0e-8);
  // implicit_tol_rel = kappa (default 0 = off; opt-in study key).  A solve
  // that starts from the step predictor measures its first change r_0 = max |x_1 -
  // x_pred|/x: the predictor (an extrapolation of the last increments, O(dt^3) off
  // under predictor_order 2) misses the fixed point by about the local truncation
  // error of the step, so an iteration error kappa r_0 is kappa x that error and
  // shrinks with dt as it does (the absolute implicit_tol does not scale).  The
  // solve then stops at tol_eff = max(implicit_tol, kappa r_0) (both the plain and
  // the contraction-estimate tests); cold starts keep implicit_tol.
  impl_tol_rel = pin->GetOrAddReal("rad_m1","implicit_tol_rel",0.0);
  impl_tol_rel_max = pin->GetOrAddReal("rad_m1","implicit_tol_rel_max",1.0e-4);
  if (impl_tol_rel < 0.0 || impl_tol_rel >= 1.0) {
    ImplFatal("<rad_m1>/implicit_tol_rel must lie in [0,1)");
  }
  // implicit_final_resid_log (DIAGNOSTIC, read only when named so that a
  // run without it has the parameter dump of the run before the key)
  if (pin->DoesParameterExist("rad_m1","implicit_final_resid_log")) {
    impl_frlog = pin->GetBoolean("rad_m1","implicit_final_resid_log");
  }
  // (rad_m1.hpp): the guarded contraction estimate and the residual stop,
  // read only when named (a run without them has the parameter dump of)
  if (pin->DoesParameterExist("rad_m1","implicit_conv_est_guard")) {
    impl_ce_guard = pin->GetReal("rad_m1","implicit_conv_est_guard");
    if (impl_ce_guard < 0.0) {
      ImplFatal("<rad_m1>/implicit_conv_est_guard must be >= 0");
    }
  }
  if (pin->DoesParameterExist("rad_m1","implicit_conv_est_nofloor")) {
    impl_ce_nofl = pin->GetBoolean("rad_m1","implicit_conv_est_nofloor");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_final_resid_stop")) {
    impl_frs = pin->GetReal("rad_m1","implicit_final_resid_stop");
    if (impl_frs < 0.0) {ImplFatal("<rad_m1>/implicit_final_resid_stop must be >= 0");}
  }
  if (pin->DoesParameterExist("rad_m1","implicit_final_resid_stop_rel")) {
    impl_frs_rel = pin->GetReal("rad_m1","implicit_final_resid_stop_rel");
    if (impl_frs_rel < 0.0) {
      ImplFatal("<rad_m1>/implicit_final_resid_stop_rel must be >= 0");
    }
  }
  impl_frs_on = (impl_frs > 0.0) || (impl_frs_rel > 0.0);
  // a lagged transverse coupling needs more outer passes than a pure column solve, so
  // the ceiling is raised (not the count: the loop still stops at convergence).
  impl_maxit = pin->GetOrAddInteger("rad_m1","implicit_maxit", full ? 200 : 30);
  impl_lin_tol = pin->GetOrAddReal("rad_m1","implicit_lin_tol",1.0e-10);
  impl_lin_maxit = pin->GetOrAddInteger("rad_m1","implicit_lin_maxit",200);
  // (USER 10-19): default bicgstab under the new default set
  // (utils/recipe.hpp); a restart whose input lacks the key / recipe none: line_jacobi
  if (!pin->DoesParameterExist("rad_m1","implicit_solver")) {
    const bool nd = recipe::NewDefaults(pin);
    recipe::Log(std::string("<rad_m1>/implicit_solver = ") +
                (nd ? "bicgstab (default)" : "line_jacobi (the old default)"));
    pin->GetOrAddString("rad_m1","implicit_solver", nd ? "bicgstab" : "line_jacobi");
  }
  {std::string sv = pin->GetOrAddString("rad_m1","implicit_solver","line_jacobi");
  if (sv.compare("line_jacobi") == 0) {
    impl_solver = M1_ISOLV_LINE_JACOBI;
  } else if (sv.compare("bicgstab") == 0) {
    impl_solver = M1_ISOLV_BICGSTAB;
  } else {
    ImplFatal("<rad_m1>/implicit_solver = '" + sv
              + "' is not a choice (line_jacobi | bicgstab)");
  }
  }
  // ---- MILESTONE 3b phase D: what happens to the off-diagonal Eddington terms, and how
  // fast the closure is allowed to move between Picard passes.  The defaults reproduce
  // phase C on every 1-D / implicit_x1 configuration (there are no off-diagonal terms and
  // no transverse closure there), and `operator` is the default in multi-D, where lagging
  // them has no fixed point in an optically thin cell.
  std::string sod = pin->GetOrAddString("rad_m1","implicit_offdiag","auto");
  bool od_auto = (sod.compare("auto") == 0);
  if (od_auto || sod.compare("lagged") == 0) {
    impl_offdiag = M1_OD_LAGGED;
  } else if (sod.compare("operator") == 0) {
    impl_offdiag = M1_OD_OPERATOR;
  } else if (sod.compare("none") == 0) {
    impl_offdiag = M1_OD_NONE;
  } else {
    ImplFatal("<rad_m1>/implicit_offdiag = '" + sod
              + "' is not a choice (auto | lagged | operator | none)");
  }
  impl_crelax = pin->GetOrAddReal("rad_m1","implicit_closure_relax",1.0);
  if (!(impl_crelax > 0.0) || impl_crelax > 1.0) {
    ImplFatal("<rad_m1>/implicit_closure_relax must lie in (0,1]");
  }
  impl_crelax_thin = pin->GetOrAddBoolean("rad_m1","implicit_closure_relax_thin",false);
  // implicit_closure_lag: DEFAULT step since defaults-0927 (2026-09-27) for transport =
  // implicit on a multi-D mesh (runs_3b5 recommends step in 2-D/3-D, every modern input
  // names it, pass is ~4x slower per simulated second; switch_inventory_2026-09-24.md);
  // pass elsewhere and on a restart whose file lacks the key.  The value is recorded.
  const bool md_full = full && pmy_pack->pmesh->multi_d;
  {std::string sc = pin->GetOrAddString("rad_m1","implicit_closure_lag",
      (md_full && !global_variable::restart_run) ? "step" : "pass");
  if (sc.compare("pass") == 0) {
    impl_clag_step = false;
  } else if (sc.compare("step") == 0) {
    impl_clag_step = true;
  } else {
    ImplFatal("<rad_m1>/implicit_closure_lag = '" + sc
              + "' is not a choice (pass | step)");
  }
  }
  // runs_5c_thinstab: the step-to-step closure relaxation of optically thin cells.
  // DEFAULT 1.5 since defaults-0927 (2026-09-27; memory m1-thin-cell-cure, merged
  // b7e48012: the cure of the multi-D thin-cell instability of the chi(f) closures, He
  // slab 332 fallbacks -> 0 at 0 extra passes, be and hesdirk2) for a chi(f) closure (m1,
  // minerbo, kershaw) under transport = implicit on a multi-D mesh with
  // implicit_closure_lag = step -- its own requirements; the fixed-tensor closures
  // (eddington, vet_sc, tau, vet_col) do not need it.  0 (off) elsewhere and on a
  // restart whose file lacks the key; the value is recorded where the default applies.
  impl_ctrelax = 0.0;
  ctr_init = false;
  {const bool chif = !(eddington || vet_sc || tau_closure);
  if (chif && md_full && impl_clag_step) {
    impl_ctrelax = pin->GetOrAddReal("rad_m1","implicit_closure_thin_relax",
                                     global_variable::restart_run ? 0.0 : 1.5);
  } else if (pin->DoesParameterExist("rad_m1","implicit_closure_thin_relax")) {
    impl_ctrelax = pin->GetReal("rad_m1","implicit_closure_thin_relax");
    // (memory): a fixed-tensor closure (eddington, vet_sc, tau, vet_col) never
    // relaxes (impl_lag_kernel: ctr && !edd && !vetsc && !tkeep && !tauc), so the 4-slab
    // closure memory stayed all-zero and was still written to every restart: off
    if (!chif) {impl_ctrelax = 0.0;}
  }
  if (impl_ctrelax < 0.0) {
    ImplFatal("<rad_m1>/implicit_closure_thin_relax must be >= 0");
  }
  }
  // ---- the TRANSVERSE realizability limiter.  `none` is not merely the default: it
  // takes the OLD expression for theta everywhere (see ImplicitTransTheta), so an input
  // file that does not name it is bitwise unchanged.
  {std::string st = pin->GetOrAddString("rad_m1","implicit_trans_limit","none");
  if (st.compare("none") == 0) {
    impl_tlim = M1_TLIM_NONE;
  } else if (st.compare("lp") == 0) {
    impl_tlim = M1_TLIM_LP;
  } else if (st.compare("lpr") == 0) {
    impl_tlim = M1_TLIM_LPR;
  } else {
    ImplFatal("<rad_m1>/implicit_trans_limit = '" + st
              + "' is not a choice (none | lp | lpr)");
  }
  }
  impl_tfmax = pin->GetOrAddReal("rad_m1","implicit_trans_fmax",1.0);
  if (!(impl_tfmax > 0.0)) {
    ImplFatal("<rad_m1>/implicit_trans_fmax must be positive");
  }
  // ---- MILESTONE 3e: ANDERSON acceleration of the outer (Picard) iteration.  `none` is
  // not merely the default: nothing below is allocated and ImplicitAccelSave/Apply are
  // never called, so an input file that does not name it is bitwise unchanged.
  {std::string sa = pin->GetOrAddString("rad_m1","implicit_accel","none");
  if (sa.compare("none") == 0) {
    impl_accel = M1_IACC_NONE;
  } else if (sa.compare("anderson") == 0) {
    impl_accel = M1_IACC_ANDERSON;
  } else {
    ImplFatal("<rad_m1>/implicit_accel = '" + sa
              + "' is not a choice (none | anderson)");
  }
  }
  // ---- MILESTONE 3g: the GAS-RADIATION energy coupling.  implicit_eos_cache defaults
  // false; neither key allocates anything nor is referenced when false.
  // implicit_gas_newton: DEFAULT true since defaults-0927 (2026-09-27;
  // runs_3g_newton_T VERDICT PASS, and every GPU timing and gate since 3g ran it;
  // switch_inventory_2026-09-24.md decision 8) for transport = implicit; false under
  // implicit_x1 (the 3a arithmetic) and on a restart whose file lacks the key.  The
  // value is recorded.
  impl_gas_newton = pin->GetOrAddBoolean("rad_m1","implicit_gas_newton",
                                         full && !global_variable::restart_run);
  impl_eos_cache = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache",false);
  dbg_etally = pin->GetOrAddBoolean("rad_m1","dbg_energy_tally",false);
  // read only when named (an input without it is unchanged, also in the
  // parameters a restart file records)
  impl_ecnt = pin->GetOrAddInteger("rad_m1","implicit_eos_cache_nt",2);
  impl_eccheck = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache_check",true);
  // the check is a MEASUREMENT only (nothing reads igm or ec_emax/ec_tmax but the final
  // report): implicit_eos_cache_check_every = N takes it on the cycles N divides.  Read
  // only when named; the default 1 is the old every-step check.
  impl_eccheck_every = 1;
  if (pin->DoesParameterExist("rad_m1","implicit_eos_cache_check_every")) {
    impl_eccheck_every = pin->GetInteger("rad_m1","implicit_eos_cache_check_every");
    if (impl_eccheck_every < 1) {
      ImplFatal("<rad_m1>/implicit_eos_cache_check_every must be >= 1");
    }
  }
  // the x1 LINE SOLVE (preconditioner and line-Jacobi pass): thomas = one thread per
  // column, serial recurrence (the original); pcr = one team per column, parallel
  // cyclic reduction in team scratch (GPU).  Same system; the answers agree to
  // round-off, not bitwise.  The gathered stack sweep (part_nblk > 1) is Thomas always.
  // DEFAULT pcr since bench/m1_defaults_0923 (15x faster line solve on the GPU, +1 %
  // wall on the host; thomas reproduces the earlier default bitwise).
  {std::string ls = pin->GetOrAddString("rad_m1","implicit_line_solver","pcr");
  if (ls.compare("thomas") == 0) {
    impl_line_solver = 0;
  } else if (ls.compare("pcr") == 0) {
    impl_line_solver = 1;
  } else {
    ImplFatal("<rad_m1>/implicit_line_solver = '" + ls
              + "' is not a choice (thomas | pcr)");
  }
  }
  impl_pcr_team = pin->GetOrAddInteger("rad_m1","implicit_pcr_team",0);
  impl_pcr_check = pin->GetOrAddBoolean("rad_m1","implicit_pcr_check",false);
  if (impl_pcr_team < 0 || impl_pcr_team > 1024) {
    ImplFatal("<rad_m1>/implicit_pcr_team must lie in [0,1024]");
  }
  // SPEED-UP 3: the host synchronisations of the BiCGStab loop.  0 = the original loop;
  // 1 = the same recurrence with fused reductions and vector updates (3 blocking
  // reductions per iteration instead of 5); 2 = 1 with alpha kept on the device, so
  // rhat.v does not block either (1 rank only; with more ranks 2 acts as 1).  Levels
  // 1 and 2 sum in a different order than 0: same answer to round-off, not bitwise.
  // DEFAULT 1 since bench/m1_defaults_0923 (0 reproduces the earlier default).
  impl_bcg_sync = pin->GetOrAddInteger("rad_m1","implicit_bcg_sync",1);
  if (impl_bcg_sync < 0 || impl_bcg_sync > 2) {
    ImplFatal("<rad_m1>/implicit_bcg_sync must be 0, 1 or 2");
  }
  if (impl_bcg_sync == 2) {
    bcg_rvd = Kokkos::View<Real, DevMemSpace>("m1_bcg_rvd");
  }
  // GPU COST OF THE KRYLOV ITERATION (bench/m1_fast_0923).  Off, nothing is allocated
  // or called and the path is bitwise the pre-0923 code.  Since
  // tests_m1/runs_3p_fastdefault they DEFAULT ON (the `FAST` set of
  // tests_m1/runs_3k_gpu3d/README_FAST.md: halo_direct, od_cache, krylov_fuse = 3,
  // op_stencil, precond = rbgs_fwd; 2.0-2.4x on the GPU, the converged answer moved at
  // the +-1-ulp level) for the closures whose tensor is fixed within a step (eddington,
  // vet_sc, tau) whenever the configuration admits them: transport = implicit with
  // bicgstab, implicit_bcg_sync = 1, implicit_line_solver = pcr, ONE MeshBlock along x1,
  // and (op_stencil) no periodic x1 wrap; with implicit_lin_cnorm > 0 krylov_fuse is 1.
  // m1 / minerbo / kershaw keep them off (not gated there).  A key the input names keeps
  // its value.
  //  implicit_halo_direct = true   the implicit exchanges as ONE on-rank copy kernel
  //    when every neighbour of every block is on its own rank at the same level
  //    (otherwise the ordinary exchange, as before).  Bitwise.
  //  implicit_od_cache = true      the off-diagonal Eddington operator from a per-cell
  //    cache of sum_e d_e P_de, fused with the 7-point row into one kernel.  Bitwise.
  //  implicit_krylov_fuse = 1      the pcr preconditioner reads and writes the Krylov
  //    vectors itself and carries the p / s updates (bitwise); = 2 also puts the
  //    (rhat,v) and (t,s),(t,t) reductions in the operator kernel (round-off: the sums
  //    are ordered differently).  Needs implicit_bcg_sync = 1, implicit_line_solver =
  //    pcr, one block along x1; 2 needs implicit_od_cache.  3 = 2 with TWO blocking
  //    reductions per iteration (ImplicitBiCGStabTwo; round-off).
  //  implicit_op_stencil = true    the frozen operator of each Picard pass written out
  //    once as a 19-point stencil (7-point row + off-diagonal Eddington terms), applied
  //    by one kernel per Krylov product (round-off: the terms are grouped differently).
  //    Needs implicit_bcg_sync = 1; not with a periodic x1 wrap.
  //  implicit_precond_float = true  the line solves of the fused path's preconditioner
  //    in float (the operator, the vectors and every reduction stay double: a
  //    preconditioner only sets the convergence rate).  Needs implicit_krylov_fuse >= 1.
  //  implicit_precond = line | rbgs | rbgs_fwd   the preconditioner of the fused path:
  //    x1 line Jacobi (the original), symmetric / forward red-black transverse line
  //    Gauss-Seidel, block-local (ImplicitPrecondX).  Needs implicit_krylov_fuse >= 1.
  const bool fixcl = eddington || vet_sc || tau_closure;
  bool fdef = fixcl && full && (impl_solver == M1_ISOLV_BICGSTAB) &&
              (impl_bcg_sync == 1) && (impl_line_solver == 1) &&
              (pmy_pack->pmesh->mesh_indcs.nx1 == pmy_pack->pmesh->mb_indcs.nx1);
  // implicit_lin_cnorm > 0 (m1-cnorm-fastpath, bsg_1001/CNORM_FAST_1003.md): the fast
  // set stays on, with implicit_krylov_fuse = 1 instead of 3.  The per-cell norm is
  // taken in the reductions of ImplicitBiCGStab (fuse 0/1/2), which ImplicitBiCGStabTwo
  // (fuse 3) does not have; fuse 1 is the p-update in the preconditioner (bitwise vs 0)
  // and costs the same as 3 (He wedge 0.392-0.422 vs 0.410-0.422 s/cycle).  Before, any
  // cnorm > 0 dropped the whole fast set (fdef = false: 2.2x slower on the He wedge).
  const bool cn_def = pin->DoesParameterExist("rad_m1","implicit_lin_cnorm") &&
                      pin->GetReal("rad_m1","implicit_lin_cnorm") > 0.0;
  // the x1 wrap as the boundary parsing below will read it (same keys, same defaults)
  const bool x1per = (ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      pmy_pack->pmesh->mesh_bcs[static_cast<int>(BoundaryFace::inner_x1)])
      == M1_IBC_PERIODIC);
  impl_halo_direct = pin->GetOrAddBoolean("rad_m1","implicit_halo_direct",fdef);
  impl_odc = pin->GetOrAddBoolean("rad_m1","implicit_od_cache",fdef);
  impl_kfuse = pin->GetOrAddInteger("rad_m1","implicit_krylov_fuse",
                                    fdef ? (cn_def ? 1 : 3) : 0);
  // a restart file written without cnorm echoes the defaulted fuse 3: with cnorm named
  // now it runs fuse 1 (round-off vs 3), echoed; a fuse 3 the input names still FATALs
  if (impl_kfuse == 3 && cn_def &&
      pin->IsParameterDefaulted("rad_m1","implicit_krylov_fuse")) {
    impl_kfuse = 1;
    pin->SetInteger("rad_m1","implicit_krylov_fuse",1);
    std::cout << "rad_m1: implicit_krylov_fuse 3 (defaulted in the restart file) -> 1, "
              << "since implicit_lin_cnorm > 0" << std::endl;
  }
  impl_stencil = pin->GetOrAddBoolean("rad_m1","implicit_op_stencil",fdef && !x1per);
  impl_prec_float = pin->GetOrAddBoolean("rad_m1","implicit_precond_float",false);
  // tests_m1/runs_4a_accel levers (implicit_vimp_fold, implicit_one_pass,
  // implicit_predictor_order, implicit_fast_kernels): resolved further down, once
  // implicit_predictor and implicit_vimp are known (their defaults depend on them)
  impl_vfold = false;
  impl_onep = 0;
  impl_onep_s = 3.0;
  for (int t = 0; t < 3; ++t) {
    onep_qa[t] = -1.0;
    onep_qb[t] = -1.0;
    onep_cnt[t] = 0.0;
    onep_off[t] = 0.0;
    onep_actr[t] = 0.0;
    onep_prb[t] = 0.0;
  }
  impl_onep_n = 0.0;
  impl_onep_nchk = 0.0;
  impl_onep_auto = false;
  impl_onep_awin = 2;
  impl_onep_arep = 64;
  impl_onep_ndis = 0.0;
  impl_onep_nprb = 0.0;
  impl_onep_nren = 0.0;
  impl_onep_nfpr = 0.0;
  ew_tight = false;
  impl_pord = 1;
  impl_opsplit = false;
  if (pin->DoesParameterExist("rad_m1","implicit_op_split_red")) {
    impl_opsplit = pin->GetBoolean("rad_m1","implicit_op_split_red");
  }
  // default OFF since audit-m1-0928 (was on since 2026-09-25, m1-fast3: -12 % on the
  // static T-S4 wedge).  On moving gas on MI300A it gains 2.2-2.6 % (1 GPU) and 0.2-2.2 %
  // (2 GPUs) per simulated second (He box, sph_wedge; 5 interleaved reps), under the ~3 %
  // bar for a result-changing default (round-off).  Opt-in by name.
  impl_opteam = pin->GetOrAddBoolean("rad_m1","implicit_op_team_red",false);
  impl_fastk = false;
  impl_odskip = false;
  if (impl_stencil && impl_bcg_sync != 1) {
    ImplFatal("<rad_m1>/implicit_op_stencil needs implicit_bcg_sync = 1");
  }
  // Spherical-polar fast-path default: rbgs_fwd.  m1-sp-order2b had made mg_gc the sp
  // default (runs_5p_coarse2: radiation-only T-S4 wedge, -18 %); on the rad-hydro
  // sph_wedge with MOVING gas rbgs_fwd is -8.7 % against mg_gc (m1-wedge,
  // docs/dev/m1_wedge_0926.md 7.2; user decision 09-26), so the sp default is rbgs_fwd
  // again.  mg_gc stays available by name; the resolved value is echoed.
  // defaults-0927 (2026-09-27): on the Cartesian fast path the default is mg (with
  // implicit_mg_levels = 3; tests_m1/runs_5m_precond: He box -7 %/cycle on 1 GPU, -6 % on
  // 2; H200 2 GPUs: plain mg fastest, handover 09-26 sect. 6) on a multi-D, single-level,
  // non-cs/sp/polar mesh, unless implicit_krylov_dev > 0 is named (it refuses mg).  A
  // restart whose file lacks the key keeps rbgs_fwd.
  auto *pmq = pmy_pack->pmesh;
  const bool kdev_named = pin->DoesParameterExist("rad_m1","implicit_krylov_dev") &&
                          pin->GetInteger("rad_m1","implicit_krylov_dev") > 0;
  // (user 10-09): mg also on the spherical-polar wedge (the AG Car A and He
  // giant productions name it there); the cubed sphere, polar boundaries, SMR/AMR, 1-D,
  // krylov_dev and the non-fast path keep their old default
  const bool mgdef = !pmq->use_cubed_sphere && !pmq->use_polar_boundary &&
                     !pmq->multilevel && pmq->multi_d && !kdev_named &&
                     !global_variable::restart_run;
  // (USER 10-30): mg_gc on a spherical-polar mesh WITH polar boundaries
  // (the full-sphere accretor; fullsph_1027/mgsm_check) on a fresh single-level 3-D run;
  // every other mesh keeps its default above, a named key wins, and a restart whose file
  // lacks the key keeps rbgs_fwd (bitwise)
  const bool gcdef = sph_geom && !cs_geom && pmq->use_polar_boundary && pmq->three_d &&
                     !pmq->multilevel && !kdev_named && !global_variable::restart_run;
  const char *pcdef = fdef ? (mgdef ? "mg" : (gcdef ? "mg_gc" : "rbgs_fwd")) : "line";
  {std::string pc = pin->GetOrAddString("rad_m1","implicit_precond", pcdef);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_precond = " << pc << std::endl;
  }
  if (pc.compare("line") == 0) {
    impl_prec = 0;
  } else if (pc.compare("rbgs") == 0) {
    impl_prec = 1;
  } else if (pc.compare("rbgs_fwd") == 0) {
    impl_prec = 2;
  } else if (pc.compare("mg") == 0) {
    impl_prec = 3;
  } else if (pc.compare("mg_gc") == 0) {
    impl_prec = 4;
  } else if (pc.compare("mg_gf") == 0) {
    impl_prec = 5;
  } else {
    ImplFatal("<rad_m1>/implicit_precond = '" + pc
              + "' is not a choice (line | rbgs | rbgs_fwd | mg | mg_gc | mg_gf)");
  }
  }
  // implicit_precond = mg (rad_m1_precond.cpp, tests_m1/runs_5m_precond): keys read only
  // under mg, so every other configuration is untouched
  mg_nlev = 0;
  mg_halo = true;
  if (impl_prec == 3) {
    // default 3 since defaults-0927 (the recommended value, runs_5m_precond; was 2)
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",3);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 2) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 2");}
  }
  // implicit_precond = mg_gc (rad_m1_precond.cpp, tests_m1/runs_5p_coarse2): mg with
  // implicit_mg_levels >= 1 (1 = rbgs_fwd + the global coarse space only) plus the
  // global band coarse space; its keys are read only under mg_gc
  if (impl_prec == 4) {
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",1);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 1) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 1 (mg_gc)");}
    gc_b2 = pin->GetOrAddInteger("rad_m1","implicit_gc_bands2",1);
    gc_b3 = pin->GetOrAddInteger("rad_m1","implicit_gc_bands3",1);
    gc_on = true;
  }
  // implicit_precond = mg_gf (rad_m1_precond.cpp, tests_m1/runs_5t_fast5box): mg plus
  // the global Fourier coarse space; its keys are read only under mg_gf
  if (impl_prec == 5) {
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",3);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 1) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 1 (mg_gf)");}
    gf_k = pin->GetOrAddInteger("rad_m1","implicit_gf_modes",2);
    if (gf_k < 0 || gf_k > 3) {ImplFatal("<rad_m1>/implicit_gf_modes must be 0..3");}
    gf_on = true;
  }
  // (rad_m1_precond.cpp): the smoother keys of mg / mg_gc / mg_gf, each read
  // only when named (default false = the old map, bitwise); mesh-dependent checks in
  // ImplicitMGBuild
  mg_x3line = false; mg_sym = false; mg_frows = false;
  // (USER 10-30): implicit_mg_x3line DEFAULT true under mg / mg_gc / mg_gf
  // on a spherical-polar 3-D mesh WITH polar boundaries on a fresh run (the full-sphere
  // accretor: the coarse polar phi rings); elsewhere false; a named key wins; a restart
  // whose file lacks the key keeps false (bitwise)
  {const bool x3def = (impl_prec >= 3) && sph_geom && !cs_geom &&
                      pmq->use_polar_boundary && pmq->three_d &&
                      !global_variable::restart_run;
  // (added to the input only where the default is on, so the restart file of such a run
  // carries it and other meshes' parameter dumps are unchanged)
  if (x3def) {
    mg_x3line = pin->GetOrAddBoolean("rad_m1","implicit_mg_x3line", true);
  } else if (pin->DoesParameterExist("rad_m1","implicit_mg_x3line")) {
    mg_x3line = pin->GetBoolean("rad_m1","implicit_mg_x3line");
  }
  if (impl_prec >= 3 && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_mg_x3line = " << (mg_x3line ? "true" : "false")
              << std::endl;
  }}
  if (pin->DoesParameterExist("rad_m1","implicit_mg_sym")) {
    mg_sym = pin->GetBoolean("rad_m1","implicit_mg_sym");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_mg_full_rows")) {
    mg_frows = pin->GetBoolean("rad_m1","implicit_mg_full_rows");
  }
  tmr_x3 = pin->DoesParameterExist("rad_m1","implicit_timers_mg_x3") &&
           pin->GetBoolean("rad_m1","implicit_timers_mg_x3");
  mg_x3halo = 1;
  if (pin->DoesParameterExist("rad_m1","implicit_mg_x3line_halo")) {
    mg_x3halo = pin->GetInteger("rad_m1","implicit_mg_x3line_halo");
    if (mg_x3halo != 1 && mg_x3halo != 2) {
      ImplFatal("<rad_m1>/implicit_mg_x3line_halo must be 1 or 2");
    }
  }
  if ((mg_x3line || mg_sym || mg_frows) && impl_prec < 3) {
    ImplFatal("<rad_m1>/implicit_mg_x3line | implicit_mg_sym | implicit_mg_full_rows "
              "need implicit_precond = mg | mg_gc | mg_gf");
  }
  if (mg_frows && !impl_stencil) {
    ImplFatal("<rad_m1>/implicit_mg_full_rows needs implicit_op_stencil = true");
  }
  if (mg_frows && mg_nlev >= 2 && !mg_halo) {
    ImplFatal("<rad_m1>/implicit_mg_full_rows needs implicit_mg_halo = true (the fine "
              "residual from the full operator reads the halo)");
  }
  if ((mg_x3line || mg_sym || mg_frows) && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> smoother: x3line " << mg_x3line << ", sym " << mg_sym
              << ", full_rows " << mg_frows << std::endl;
  }
  // implicit_bcg_rho_direct (runs_5p_coarse2; read only when named, default false =
  // the recurrence, bitwise): implicit_krylov_fuse = 3 sums (rhat, r) directly
  impl_rho_direct = false;
  if (pin->DoesParameterExist("rad_m1","implicit_bcg_rho_direct")) {
    impl_rho_direct = pin->GetBoolean("rad_m1","implicit_bcg_rho_direct");
  }
  if (impl_kfuse < 0 || impl_kfuse > 3) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse must be 0, 1, 2 or 3");
  }
  if (impl_kfuse == 3 && pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0) > 0.0) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse = 3 does not take implicit_lin_cnorm");
  }
  if (impl_kfuse > 0 && (impl_bcg_sync != 1 || impl_line_solver != 1)) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse needs implicit_bcg_sync = 1 and "
              "implicit_line_solver = pcr");
  }
  if (impl_kfuse >= 2 && !impl_odc) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse >= 2 needs implicit_od_cache = true");
  }
  if (impl_prec_float && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond_float needs implicit_krylov_fuse >= 1");
  }
  if (impl_prec > 0 && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond = rbgs needs implicit_krylov_fuse >= 1");
  }
  // multi-rank Krylov (rad_m1_krylov.cpp, tests_m1/runs_3w_krylov).  krylov_pipe
  // defaults OFF.  halo_mpi (bitwise the ordinary exchange) defaults ON since
  // m1-defaults wherever its preconditions hold: implicit_halo_direct, a same-level mesh
  // (no SMR/AMR), no cubed-sphere seams and no polar boundary; otherwise the ordinary
  // exchange, silently.  A key the input (or restart echo) names keeps its value.
  // m1-sphhalo (tests_m1/runs_5i_sphhalo): also on the spherical-polar WEDGE (no pole:
  // every ghost is the plain copy, no vector slot flips), except on a restart whose file
  // lacks the key (written before m1-sphhalo), which keeps the ordinary exchange.
  impl_kpipe = pin->GetOrAddBoolean("rad_m1","implicit_krylov_pipe",false);
  {auto *pmh = pmy_pack->pmesh;
  const bool hmdef = impl_halo_direct && !pmh->multilevel && !pmh->use_cubed_sphere &&
                     !pmh->use_polar_boundary &&
                     !(pmh->use_spherical_polar && global_variable::restart_run);
  impl_halo_mpi = pin->GetOrAddBoolean("rad_m1","implicit_halo_mpi",hmdef);}
  hm_state = 0;
  hm_comm = nullptr;
  // implicit_halo_ipc (m1-perf-0928, rad_m1_krylov.cpp): default false
  impl_halo_ipc = pin->GetOrAddBoolean("rad_m1","implicit_halo_ipc",false);
  if (impl_halo_ipc && !impl_halo_mpi) {
    ImplFatal("<rad_m1>/implicit_halo_ipc needs implicit_halo_mpi = true");
  }
  // implicit_halo_overlap (rad_m1_krylov.cpp, tests_m1/runs_3y_halo_overlap: round-off
  // vs off, restarts bitwise).  DEFAULT true since m1-accmerge wherever it is valid:
  // implicit_halo_mpi on and more than one rank; otherwise false, silently.  A restart
  // whose file lacks the key (written before m1-accmerge, which read it only when
  // named) keeps false; the resolved value is echoed; explicit input overrides.
  impl_halo_ovl = pin->GetOrAddBoolean("rad_m1","implicit_halo_overlap",
                      impl_halo_mpi && (global_variable::nranks > 1) &&
                      !global_variable::restart_run);
  // implicit_halo_ovl_faces (tests_m1/runs_4l_sync): default OFF since audit-m1-0928
  // (was true wherever the overlap is on, from 09-24 H200 weak scaling).  On MI300A with
  // moving gas it gains -0.1..+1.0 % (2 GPUs, one node) and -1.4..+0.8 % (2 nodes x 2
  // GPUs) per simulated second (He box, sph_wedge): no gain for a result-changing
  // (round-off) default.  Opt-in by name; the resolved value is echoed.
  impl_ovl_faces = pin->GetOrAddBoolean("rad_m1","implicit_halo_ovl_faces",false);
  for (int f = 0; f < 6; ++f) {hm_face[f] = 1;}
  if (impl_halo_ovl && !impl_halo_mpi) {
    ImplFatal("<rad_m1>/implicit_halo_overlap needs implicit_halo_mpi = true");
  }
  if (impl_kpipe && impl_kfuse != 3) {
    ImplFatal("<rad_m1>/implicit_krylov_pipe needs implicit_krylov_fuse = 3");
  }
  if (impl_halo_mpi && !impl_halo_direct) {
    ImplFatal("<rad_m1>/implicit_halo_mpi needs implicit_halo_direct = true");
  }
  // implicit_krylov_dev (tests_m1/runs_4k_launch): default 0, read only when named
  impl_kdev = 0;
  kdev_last[0] = kdev_last[1] = kdev_last[2] = 0;
  kdev_slot = 0;
  kdev_x1p = -1;
  kdev_nchk = 0.0;
  kdev_nq = 0.0;
  impl_kdev_halo = 0;
  if (pin->DoesParameterExist("rad_m1","implicit_krylov_dev")) {
    impl_kdev = pin->GetInteger("rad_m1","implicit_krylov_dev");
    if (impl_kdev < 0) {ImplFatal("<rad_m1>/implicit_krylov_dev must be >= 0");}
    if (impl_kdev > 0 && impl_prec >= 3) {
      ImplFatal("<rad_m1>/implicit_krylov_dev does not take implicit_precond = mg(_gc)");
    }
    if (impl_kdev > 0 && (impl_kfuse != 3 || impl_kpipe || !impl_halo_direct)) {
      ImplFatal("<rad_m1>/implicit_krylov_dev needs implicit_krylov_fuse = 3, "
                "implicit_halo_direct = true and implicit_krylov_pipe = false");
    }
    impl_kdev_halo = 0;
    if (impl_kdev > 0 && pin->DoesParameterExist("rad_m1","implicit_krylov_dev_halo")) {
      impl_kdev_halo = pin->GetInteger("rad_m1","implicit_krylov_dev_halo");
    }
    if (impl_kdev > 0) {
      kdv = DvceArray1D<Real>("m1_kdv", M1_KD_SIZE);
      kdh = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_kdh", M1_KD_SIZE);
    }
  }
  // implicit_op_check (rad_m1_opcheck.cpp, tests_m1/gates): read only when named
  impl_opchk = 0;
  opchk_n = 0;
  impl_opchk_tol = 1.0e-12;
  if (pin->DoesParameterExist("rad_m1","implicit_op_check")) {
    impl_opchk = pin->GetInteger("rad_m1","implicit_op_check");
    if (pin->DoesParameterExist("rad_m1","implicit_op_check_tol")) {
      impl_opchk_tol = pin->GetReal("rad_m1","implicit_op_check_tol");
    }
  }
  impl_dump_cyc = -1;
  impl_dump_done = false;
  if (pin->DoesParameterExist("rad_m1","implicit_dump_op")) {
    impl_dump_cyc = pin->GetInteger("rad_m1","implicit_dump_op");
  }
  // the Picard pass count (bench/m1_picard_0923): a per-pass log, off by default
  impl_plog = pin->GetOrAddInteger("rad_m1","implicit_picard_log",0);
  // implicit_tsolve_opac (nc_cure_1002): the Newton-fallback gas-temperature root find
  // with kappa_P(T), kappa_E(T) evaluated at the trial T (M1ImplTemperatureOpac) instead
  // of at the lagged iterate; it removes the |dT|/T = 0.365 limit cycles of shock-heated
  // near-void cells (He presn wedge, /viper/ptmp2/jinma/he_nc_cure_1002).  ...only from
  // Picard pass implicit_tsolve_opac_start on: a solve that converges before it is
  // bitwise the frozen-opacity one.  DEFAULT true with start 20 on fresh runs (user
  // 10-02; part of the R2 default with mode 1 and switch 10 below); a restart whose
  // file lacks the key keeps the old behaviour (off) and says so.  The resolved values
  // are recorded.
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac")) {
    impl_tsolve_opac = pin->GetBoolean("rad_m1","implicit_tsolve_opac");
  } else if (global_variable::restart_run) {
    impl_tsolve_opac = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_tsolve_opac: keeping the old "
                << "default (off); set it to switch" << std::endl;
    }
  } else {
    impl_tsolve_opac = pin->GetOrAddBoolean("rad_m1","implicit_tsolve_opac",true);
  }
  impl_tsolve_opac_start = pin->GetOrAddInteger("rad_m1","implicit_tsolve_opac_start",20);
  // implicit_tsolve_opac_mode (tsolve_root_fix_1003): 0 = the root find above; bit 1
  // (value 1) takes the root of the kappa(T) equation NEAREST the lagged T
  // (M1ImplTemperatureOpacNear) instead of the halve/double bracket's; bit 2 (value 2)
  // is a slope guard: a cell whose d ln kappa_P / d ln T at the lagged T is below
  // implicit_tsolve_opac_slope_min (default -4, where kappa_P T^4 falls with T) takes
  // the frozen-opacity root find instead.  3 = both.  Bit 4 (value 8, OPT-IN,
  // r2b_shift_1003) takes the stable root nearest the end point of the exact exchange
  // ODE from T_gn (M1ImplTemperatureOpacOde; wins over bit 1).  BSG R2E (mode 8, rst
  // 00008 on, bsg_1001/r2b_shift_1003): robust (no divergence) but ~10 % slower than
  // mode 1 (2.12 vs 1.93 s/cycle) with more NON-CONV (54 vs 31, worst resid 7.9e-4 vs
  // 3.8e-5), so the default stays 1.  DEFAULT 1 on fresh runs (R2
  // default, user 10-03): mode 0 picks the wrong root of the 3-root gas-T equation at H
  // recombination (kappa_P ~ T^30) and 2-cycles; with mode 1 the BSG envelope (rst
  // 00008, switch 10) runs 0 FATAL / 4 NON-CONV where mode 0 diverged
  // (bsg_1001/TSOLVE_FIX_1003.md, AB_SOLVER_1003.md), BSG R2 vs A0 lies within 2-3x the
  // noise spread (bsg_1001/r2b_shift_1003), He wedge/box: 0 NON-CONV, box bitwise
  // (he_mltpp_1002/GATE_MODE1.md).  A restart whose file lacks the key keeps the old
  // behaviour (0) and says so; a named value always wins.
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac_mode")) {
    impl_tsolve_opac_mode = pin->GetInteger("rad_m1","implicit_tsolve_opac_mode");
  } else if (global_variable::restart_run) {
    impl_tsolve_opac_mode = 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_tsolve_opac_mode: keeping "
                << "the old default (0); set it to switch" << std::endl;
    }
  } else {
    impl_tsolve_opac_mode = pin->GetOrAddInteger("rad_m1","implicit_tsolve_opac_mode",1);
  }
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac_slope_min")) {
    impl_tsolve_opac_smin = pin->GetReal("rad_m1","implicit_tsolve_opac_slope_min");
  }
  // DEBUG (nc_cure_1002): implicit_nc_dump = N > 0 prints, in the last N passes of a
  // solve that runs to implicit_maxit, the local state of each rank's worst-residual
  // cell (optical depths, coupling stiffness, energy ratio, beta, f, the T-Newton
  // flags).  Read only when named; prints only; no number of the solve changes.
  impl_ncdump = pin->DoesParameterExist("rad_m1","implicit_nc_dump") ?
                pin->GetInteger("rad_m1","implicit_nc_dump") : 0;
  // DEBUG: dbg_eint_cells = N > 0 prints, in every solve
  // whose gas write-back gives eint <= 0 somewhere, up to N of those cells per rank with
  // their state.  Read only when named; prints only; no number of the solve changes.
  g_dbg_eint_n = pin->DoesParameterExist("rad_m1","dbg_eint_cells") ?
               pin->GetInteger("rad_m1","dbg_eint_cells") : 0;
  // ...and with time_scheme = be, where Time2Init returns before reading it, the
  // DEBUG (he_estall_1003): implicit_stall_trace = N > 0: in a solve still unconverged at
  // pass implicit_stall_trace_p0 (default 40), fix each rank's worst cell of that pass
  // and print its NCD line (tag 3) before step (c) and (tag 2) after the residual of
  // every pass p0 .. p0+N-1, for at most implicit_stall_trace_max (default 30) solves.
  // Prints only; read only when named.
  if (pin->DoesParameterExist("rad_m1","implicit_stall_trace")) {
    impl_strace = pin->GetInteger("rad_m1","implicit_stall_trace");
    impl_strace_p0 = pin->GetOrAddInteger("rad_m1","implicit_stall_trace_p0",40);
    impl_strace_max = pin->GetOrAddInteger("rad_m1","implicit_stall_trace_max",30);
  }
  tmr_c0 = pin->DoesParameterExist("rad_m1","implicit_timers") ?
           pin->GetInteger("rad_m1","implicit_timers") : 0;   // read only when named
  // ...and the options that cut it.  Since bench/m1_defaults_0923 they DEFAULT ON for
  // the closures that do not read the iterate (eddington, vet_sc, tau: 1.6x on the GPU,
  // statistics moved at the solver-tolerance level) and stay OFF for m1 / minerbo /
  // kershaw, whose closure moves with the iterate: there the lresid test is NOT
  // redundant and the inexact (Eisenstat-Walker) pass 0 raised NON-CONVERGED steps
  // (tests_m1/runs_3k_gpu3d/README_DEFAULTS.md).  The earlier default path is
  // implicit_lres_test = true, implicit_conv_est = false, implicit_lin_ew_max = 0,
  // implicit_predictor = none (with implicit_bcg_sync = 0 and implicit_line_solver =
  // thomas it is bitwise the pre-flip code):
  //  implicit_lres_test = false  drop the pass-to-pass transverse-change test under
  //    bicgstab, where the transverse coupling is IN the operator and is solved to
  //    implicit_lin_tol on every pass (the test lags the Picard test by one pass);
  //  implicit_conv_est = true    stop once q/(1-q) times the last change (q = the
  //    measured contraction of the last two passes, required < 1/2) is below
  //    implicit_tol, i.e. without the confirming pass;
  //  implicit_lin_ew_max > 0     Eisenstat-Walker (choice 2) inner tolerance
  //    max(lin_tol max|b|, eta_k max|r0|), eta_0 = ew_max,
  //    eta_k = min(ew_max, gamma (|r0_k|/|r0_{k-1}|)^2) with the gamma eta_{k-1}^2
  //    safeguard (bcg_sync >= 1 only);
  //  implicit_predictor = step   start the loop from the previous step's implicit
  //    increment scaled by dt/dt_prev (closures that do not read the iterate only;
  //    the default for them since 0923, restart-safe).
  //  The Eisenstat-Walker default is 1e-2 with bcg_sync >= 1 and 0 (off) with
  //  bcg_sync = 0, so an input that asks for the original BiCGStab loop still runs.
  //  implicit_res_dmin, implicit_res_rmax (read only when named; default 0 = off): the
  //    Picard convergence norm (resid and lresid) leaves out the cells with gas density
  //    rho < res_dmin (needs <hydro>) or radius x1v > res_rmax (spherical polar only),
  //    e.g. the density-floor top whose residual stalls at 1e-7..1e-6.  Those cells are
  //    still solved on every pass; only the stopping test ignores them.  Off: bitwise.
  if (pin->DoesParameterExist("rad_m1","implicit_res_dmin")) {
    impl_res_dmin = pin->GetReal("rad_m1","implicit_res_dmin");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_res_rmax")) {
    impl_res_rmax = pin->GetReal("rad_m1","implicit_res_rmax");
  }
  //  implicit_realisable_coupling (gd_physfix_1005, read only when named; default false =
  //    bitwise the old code): the gas sees the REALISABLE radiation flux.  In optically
  //    thin cells with c dt/dx >> 1 the Picard iterate of the cell flux is unrealisable
  //    (BSG far-thin top: mean |F| = 10-20 c E, ~200 c E in the cells that stalled);
  //    the write-back clips it to |F| = c E (M1ApplyLimits), but the comoving correction
  //    E0 - E = -2 beta.F/c + O(beta^2) and the momentum deposit dt (rho k_t)_f F0_f/c
  //    took the UNCLIPPED flux: E0 - E came out -0.85 E (|E0 - E| <= 2 |beta| E for any
  //    intensity field), which moved the gas-T equation into the kappa_P(T) valley of the
  //    low-density table and made it 3-rooted (the Picard T cycle), and the floor gas got
  //    10-200x the largest force any radiation field of energy E can exert (rho k E).
  //    With the key, per cell, s = min(1, c E/|F_iter|) scales the flux in beta.F of E0,
  //    and the momentum deposit (and so its work, which the radiation loses: total
  //    energy stays exact).  s = 1 wherever the iterate is realisable.
  if (pin->DoesParameterExist("rad_m1","implicit_realisable_coupling")) {
    impl_real_couple = pin->GetBoolean("rad_m1","implicit_realisable_coupling");
  }
  // implicit_cs_e0_metric (cubed sphere; read only when named): the work
  //    array's transverse v and F are FACE-NORMAL components (v.n_xi, v.n_eta) on the
  //    skew panel basis, n_xi.n_eta = -cos_cell.  Their Euclidean products miss the
  //    metric G^-1 = [[1, c], [c, 1]]/(1 - c^2) in beta.F, beta^2 and beta.P.beta of
  //    E0 - E (EQ_AUDIT_1112 #4).  With the key the three are built with it, i.e. with
  //    the same metric as the work v^a dm_a and the KE of the write-back.
  //    prod-1114 (user 11-12): DEFAULT true (recorded); name false for the old E0.
  impl_cs_e0m = pin->GetOrAddBoolean("rad_m1","implicit_cs_e0_metric",true);
  // (read only when named):
  //  implicit_cs_closure_metric (cs): |F|, f and the unit n of a chi(f) closure with the
  //    panel metric, and A = v + D.v, beta.D.beta with D the MIXED tensor on face-normal
  //    components (the n = F/|F| closures; vet_col's n = r is unchanged by it)
  //  implicit_comoving_lat (sp, cs): the lateral D_r,lat and tangential (cs Q, sp a and
  //    D_tp) parts of the VET tensor in A = v + D.v of F0 and in beta.D.beta of E0,
  //    wherever the face equations carry them (M1LatDx)
  if (pin->DoesParameterExist("rad_m1","implicit_cs_closure_metric")) {
    impl_cs_cnm = pin->GetBoolean("rad_m1","implicit_cs_closure_metric");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_comoving_lat")) {
    impl_comov_lat = pin->GetBoolean("rad_m1","implicit_comoving_lat");
  }
  //  implicit_thin_freeze (m1-picard-aa, 095ca6ee): cells with c dt rho kappa_P below it
  //    at pass 0 keep their start-of-solve opacities for every Picard pass and take the
  //    frozen-opacity gas-T find.  DEFAULT 1e-2 on fresh runs (user 10-04): on the He
  //    presn 128^2 wedge (he_mltpp_1002/FREEZE_AB.md, rst 00055, 1.5 ks) it is neutral
  //    within the cfl-0.2999 noise member below the photosphere, -32 % s/cycle, Picard
  //    mean 39.6 -> 16.3, 0 NON-CONV; dt -0.9 % (as on the BSG, PICARD_ROBUST_1003.md).
  //    A restart whose file lacks the key keeps the old behaviour (0 = off) and says so;
  //    a named value always wins.
  if (pin->DoesParameterExist("rad_m1","implicit_thin_freeze")) {
    impl_thin_frz = pin->GetReal("rad_m1","implicit_thin_freeze");
  } else if (global_variable::restart_run) {
    impl_thin_frz = 0.0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_thin_freeze: keeping the old "
                << "default (0 = off); set it to switch" << std::endl;
    }
  } else {
    impl_thin_frz = pin->GetOrAddReal("rad_m1","implicit_thin_freeze",1.0e-2);
  }
  if (impl_thin_frz < 0.0) {ImplFatal("<rad_m1>/implicit_thin_freeze must be >= 0");}
  //  implicit_resid_fatal: a Picard solve that ends NON-CONVERGED with resid > this
  //    value (or a non-finite resid or lin_resid) is a diverged solve: FATAL before its
  //    state spreads (BSG production 10-01: one 200-pass solve with resid 1.2e6 at
  //    t = 2.91e5 s, the whole domain NaN ~15 cycles later).  0 = off.  DEFAULT 1e2
  // (/viper/ptmp2/jinma/defaults_gate_1002/item5): the largest
  //    NON-CONVERGED resid of a run that went on cleanly was 0.997 (He presn onv128;
  //    scout64 0.37); every resid > 1e3 seen came from a run that was failing.  A restart
  //    whose file lacks the key keeps the old default (off) and says so; the resolved
  //    value is recorded.
  if (pin->DoesParameterExist("rad_m1","implicit_resid_fatal")) {
    impl_res_fatal = pin->GetReal("rad_m1","implicit_resid_fatal");
    if (!(impl_res_fatal >= 0.0)) {
      ImplFatal("<rad_m1>/implicit_resid_fatal must be >= 0");
    }
  } else if (recipe::NewDefaults(pin)) {
    // (USER 10-19): 1e-2 under the new default set (utils/recipe.hpp)
    impl_res_fatal = pin->GetOrAddReal("rad_m1","implicit_resid_fatal",1.0e-2);
    recipe::Log("<rad_m1>/implicit_resid_fatal = 1e-2 (default)");
  } else if (global_variable::restart_run) {
    recipe::Log("<rad_m1>/implicit_resid_fatal = 0 (off: restart input has no key, the "
                "old default kept; set it to switch)");
  } else {
    impl_res_fatal = pin->GetOrAddReal("rad_m1","implicit_resid_fatal",1.0e2);
    recipe::Log("<rad_m1>/implicit_resid_fatal = 1e2 (recipe none: the old default)");
  }
  if (impl_res_dmin < 0.0 || impl_res_rmax < 0.0) {
    ImplFatal("<rad_m1>/implicit_res_dmin and implicit_res_rmax must be >= 0");
  }
  if (impl_res_rmax > 0.0 && !pmy_pack->pmesh->use_spherical_polar) {
    ImplFatal("<rad_m1>/implicit_res_rmax needs a spherical-polar mesh");
  }
  if (impl_res_dmin > 0.0 && pmy_pack->phydro == nullptr) {
    ImplFatal("<rad_m1>/implicit_res_dmin needs <hydro> (the gas density)");
  }
  impl_res_mask = (impl_res_dmin > 0.0) || (impl_res_rmax > 0.0);
  //  implicit_resid_fatal_masked (resmask_1003; read only when the mask is on, so a
  //    mask-off run is unchanged): FATAL when the largest residual of the excluded cells
  //    at the last pass of ANY solve (converged or not) is above this value or
  //    non-finite.  Default 1.0: the worst excluded residual of the BSG res_dmin = 1e-16
  //    A/B arm over 881 steps was 9.6e-2 (/viper/ptmp2/jinma/bsg_1001/SEAM_AB_1003.md).
  //    0 = off.
  if (impl_res_mask) {
    impl_res_fatal_masked = pin->GetOrAddReal("rad_m1","implicit_resid_fatal_masked",1.0);
    if (!(impl_res_fatal_masked >= 0.0)) {
      ImplFatal("<rad_m1>/implicit_resid_fatal_masked must be >= 0");
    }
  }
  impl_lres_test = pin->GetOrAddBoolean("rad_m1","implicit_lres_test",!fixcl);
  impl_conv_est = pin->GetOrAddBoolean("rad_m1","implicit_conv_est",fixcl);
  // implicit_src_stable (stall_1002): the cancellation-free form of the gas-eliminated
  // source row.  DEFAULT true (user 10-02) for a fresh run; false (the plain form,
  // bitwise) on a restart whose file lacks the key, with one line on rank 0.  The value
  // is recorded, so later restarts keep it; a named key always wins.
  {const bool named = pin->DoesParameterExist("rad_m1","implicit_src_stable");
  impl_src_stable = pin->GetOrAddBoolean("rad_m1","implicit_src_stable",
                                         !global_variable::restart_run);
  if (!named && global_variable::restart_run && global_variable::my_rank == 0) {
    std::cout << "rad_m1: implicit_src_stable = false (restart file lacks the key; "
              << "name it to switch the stable source row on)" << std::endl;
  }
  }
  impl_ew_max = pin->GetOrAddReal("rad_m1","implicit_lin_ew_max",
                                  (fixcl && impl_bcg_sync >= 1) ? 1.0e-2 : 0.0);
  impl_ew_gam = pin->GetOrAddReal("rad_m1","implicit_lin_ew_gamma",0.9);
  //  implicit_lin_cnorm > 0     the inner test on max_i |r_i|/(s_i E^k_i) < cnorm,
  //    s_i = 1 + SRCB_i the row excess: a bound on the per-cell relative error of E
  //    the residual implies (bcg_sync >= 1 only).
  impl_lin_cnorm = pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0);
  if (impl_lin_cnorm < 0.0 || (impl_lin_cnorm > 0.0 && impl_bcg_sync == 0)) {
    ImplFatal("<rad_m1>/implicit_lin_cnorm must be >= 0 and needs "
              "implicit_bcg_sync >= 1");
  }
  // implicit_bcg_max_restarts / implicit_bcg_fallback (defaults 2 / line
  // = the original, bitwise).  After more than max_restarts breakdowns (a vanishing
  // rho / omega / rhat.v, or a recurrence residual that met the tolerance while the true
  // residual did not) the fused BiCGStab gives up.  `line` then DISCARDS the Krylov
  // iterate and makes one line-Jacobi update from E^k (ImplicitBiCGStabEnd); `best`
  // first forms the true residual of the iterate and keeps the iterate whenever that
  // residual is smaller than the initial one (the same norm the tolerance uses), so the
  // Picard pass gets the better of the two instead of always the cruder.  Fused path
  // (bcg_sync >= 1, krylov_fuse <= 2) only.
  impl_bcg_maxrst = pin->GetOrAddInteger("rad_m1","implicit_bcg_max_restarts",2);
  {std::string fb = pin->GetOrAddString("rad_m1","implicit_bcg_fallback","line");
  if (fb.compare("line") == 0) {
    impl_bcg_keep = false;
  } else if (fb.compare("best") == 0) {
    impl_bcg_keep = true;
  } else {
    ImplFatal("<rad_m1>/implicit_bcg_fallback = '" + fb
              + "' is not a choice (line | best)");
  }
  }
  if (impl_bcg_maxrst < 0) {
    ImplFatal("<rad_m1>/implicit_bcg_max_restarts must be >= 0");
  }
  // implicit_precond_fallback (default none = bitwise the old path).
  // A linear solve whose BiCGStab stops at implicit_lin_maxit without meeting its test
  // (not a breakdown: that is implicit_bcg_fallback) is redone ONCE with this
  // preconditioner, from the same initial guess (the Picard iterate) and with the same
  // stopping test (the Eisenstat-Walker state is put back first).  The redo's iterate
  // replaces the first one whatever its outcome.  Motivation: cbh vacuum beams
  // (blend_1022/CBH_RESULTS.md), where mg's coarse corrections spread the error left at
  // the cap into near-vacuum cells and the Picard test never passes; rbgs / line keep it
  // local.  A fallback equal to implicit_precond is a no-op.
  {std::string pf = pin->GetOrAddString("rad_m1","implicit_precond_fallback","none");
  if (pf.compare("none") == 0) {
    impl_prec_fb = -1;
  } else if (pf.compare("line") == 0) {
    impl_prec_fb = 0;
  } else if (pf.compare("rbgs") == 0) {
    impl_prec_fb = 1;
  } else {
    ImplFatal("<rad_m1>/implicit_precond_fallback = '" + pf
              + "' is not a choice (none | rbgs | line)");
  }
  }
  if (impl_prec_fb == impl_prec) {impl_prec_fb = -1;}
  if (impl_prec_fb > 0 && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond_fallback = rbgs needs "
              "implicit_krylov_fuse >= 1");
  }
  // implicit_face_opac_n (F2): rad_m1.hpp.  defaults-1009 (user 10-09):
  // default true on a fresh run; a restart whose file lacks the key keeps false
  impl_face_ktn = pin->GetOrAddBoolean("rad_m1","implicit_face_opac_n",
                                       !global_variable::restart_run);
  if (impl_ew_max < 0.0 || impl_ew_max >= 1.0 || !(impl_ew_gam > 0.0)) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max must lie in [0,1), ew_gamma > 0");
  }
  if (impl_ew_max > 0.0 && impl_bcg_sync == 0) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max needs implicit_bcg_sync >= 1");
  }
  // implicit_det_reduce (bsg_1001/DETERMINISM.md): the Krylov sums of the fused
  // BiCGStab in a fixed order (M1DetReduce) and the cross-rank sums in rank order, so
  // that two runs of one binary on one decomposition are bitwise identical on the GPU.
  // DEFAULT: on for a device build (GPU: measured 2-7 % CHEAPER per cycle than the
  // Kokkos reductions on the BSG 3-D and He ext wedges, DETERMINISM.md) wherever it is
  // supported, off on a host build (the strided slots cost 1.7x on the CPU column).
  // Both settings are run-to-run bitwise on the GPU for a fixed decomposition.
  const bool det_ok = (impl_solver == M1_ISOLV_BICGSTAB) && (impl_bcg_sync == 1) &&
                      (impl_kfuse != 2) && !impl_kpipe && (impl_kdev == 0) &&
                      (impl_accel == M1_IACC_NONE);
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA) \
    || defined(KOKKOS_ENABLE_SYCL)
  const bool det_def = det_ok;
#else
  const bool det_def = false;
#endif
  impl_det = pin->GetOrAddBoolean("rad_m1","implicit_det_reduce",det_def);
  if (impl_det && !det_ok) {
    ImplFatal("<rad_m1>/implicit_det_reduce supports implicit_solver = bicgstab with "
              "implicit_bcg_sync = 1, implicit_krylov_fuse = 0, 1 or 3 (not krylov_pipe "
              "or krylov_dev) and implicit_accel = none "
              "(solver " + std::to_string(impl_solver) + ", bcg_sync " +
              std::to_string(impl_bcg_sync) + ", krylov_fuse " +
              std::to_string(impl_kfuse) + ", accel " + std::to_string(impl_accel) + ")");
  }
  if (impl_det) {det_part = DvceArray1D<Real>("m1_det_part", 4*1024 + 4);}
  {std::string pr = pin->GetOrAddString("rad_m1","implicit_predictor",
                                        fixcl ? "step" : "none");
  // predictor = step is the default for the closures that do not read the iterate: its
  // state (ipred, pred_ok, pred_dt) travels in the restart file (radm1::kM1PredRstMagic,
  // tests_m1/runs_3k_gpu3d/README_PREDRST.md), so a restart stays bitwise.
  if (pr.compare("none") == 0) {
    impl_pred = false;
  } else if (pr.compare("step") == 0) {
    impl_pred = true;
  } else {
    ImplFatal("<rad_m1>/implicit_predictor = '" + pr
              + "' is not a choice (none | step)");
  }
  // implicit_predictor_frozen: with a fixed-tensor closure the
  // step predictor was switched off under berthon / blend + implicit_recon_lag = step
  // except on the sp wedge (spfx: the face reduced flux is f(D_rr) of the lagged
  // closure there, not the iterate's).  Under implicit_flux_beam = halfrange with an
  // implicit_thin_switch every face coefficient is built from the formal solution and
  // E^n alone (M1ThinFace), on every mesh, so the predictor moves only the starting
  // point there as on the wedge.  true: that rule on every mesh (opt-in); false: the old
  // sp-only gate (DEFAULT, fresh and restart; echoed).  the fresh default true
  // left a Picard-tolerance error at implicit_tol 1e-12 that does not shrink with dt (the
  // Cartesian rw_p10t1 / rw_p10t100 T ladders lost their order: E 1.92 -> 0.65 at l 2-4;
  // clean at tol 1e-14 or with false), so FIN3 keeps the old gate.
  if (impl_pred) {
    impl_pred_frz = pin->GetOrAddBoolean("rad_m1","implicit_predictor_frozen",false);
  }
  }
  if (impl_ecnt < 0 || impl_ecnt > M1_EC_NTMAX) {
    ImplFatal("<rad_m1>/implicit_eos_cache_nt must lie in [0,8]");
  }
  impl_and_m = pin->GetOrAddInteger("rad_m1","implicit_anderson_m",5);
  impl_and_beta = pin->GetOrAddReal("rad_m1","implicit_anderson_beta",1.0);
  impl_and_start = pin->GetOrAddInteger("rad_m1","implicit_anderson_start",1);
  if (impl_accel == M1_IACC_ANDERSON) {
    if (impl_and_m < 1 || impl_and_m > M1_AND_MMAX) {
      ImplFatal("<rad_m1>/implicit_anderson_m must lie in [1,10]");
    }
    if (!(impl_and_beta > 0.0) || impl_and_beta > 1.0) {
      ImplFatal("<rad_m1>/implicit_anderson_beta must lie in (0,1]");
    }
    if (impl_and_start < 1) {
      ImplFatal("<rad_m1>/implicit_anderson_start must be >= 1 (pass 0 has no history)");
    }
  }
  // implicit_opac_update: default true with force_reference = wb_arad
  // (m1-keydefault-0927, docs/dev/ke_dt_0926.md 8), false otherwise
  impl_opac_update = pin->GetOrAddBoolean("rad_m1","implicit_opac_update",
                                          force_ref == M1_FREF_WB_ARAD);
  // implicit_opac_newton (m1-perf-0928): the Picard loop under implicit_opac_update is a
  // fixed-point iteration on the FLUX opacity rho kappa_T(T) of the x1 face terms, and
  // it converges only linearly (~0.17 per pass on the He box: 6.8 passes per solve at
  // tol 1e-8 against 2.0 with the opacities frozen).  This key adds the Newton term of
  // that dependence to the x1 rows (the chain dG/dkappa_face dkappa/dT dT/dE' through
  // the local gas Newton step), so the pass converges ~quadratically.  The added terms
  // are proportional to the predicted temperature change of the gas Newton step, which
  // vanishes at the fixed point: the converged state is the same.
  // DEFAULT true (m1-perf-0928, user) wherever it applies: implicit_opac_update and
  // implicit_gas_newton on; hydro or MHD (FluidRef), Cartesian or sp rows.  False on a
  // restart whose file lacks the key (the implicit_gas_newton convention).  Named true
  // without its two prerequisites it is fatal (below).
  // (USER 11-08 "make some defaults", r8 key set): DEFAULT false (fresh and
  // restart); recipe none keeps the old rule above
  const bool rcp_none = recipe::GetMode(pin) == recipe::Mode::none;
  impl_opac_newton = pin->GetOrAddBoolean("rad_m1","implicit_opac_newton",
                                          !rcp_none ? false :
                                          (impl_opac_update && impl_gas_newton &&
                                           !global_variable::restart_run));
  // implicit_opac_newton_guard (m1-opn-guard, 09-29): the Newton term above has the sign
  // of the face flux times d kappa/dT and no bound; in fast optically thin cells (He
  // presn wedge, across the Fe bump) it took a row diagonal of +802 to -9200, the row
  // lost the M-matrix property and the solved E went negative (floor clips, stage
  // failures, NaN).  A face's term is kept in a row only if the diagonal stays >= guard
  // x its value without it; otherwise that face is left Picard in that row (the fixed
  // point is the same: the term vanishes there).  <= 0 turns the guard off (the
  // unguarded rows of m1-perf-0928, bitwise).
  impl_opn_guard = pin->GetOrAddReal("rad_m1","implicit_opac_newton_guard",0.5);
  // implicit_gas_newton_u4 (default false, opt-in study key): the gas Newton update of
  // step (f) in a cell whose linearised row is EMISSION-dominated (4 c dt rho kappa_P a
  // T^3 > rho c_v) moves u = T^4 by the Newton step, T' = (T^4 + 4 T^3 dT)^(1/4), instead
  // of T by dT: then a T'^4 is exactly the emission the radiation row was built with, and
  // a hot start (T >> T_eq in a tenuous cell: plain Newton on T^4 shrinks T by only 1/4
  // per pass) lands in one step.  Same fixed point (dT = 0 there).
  impl_gn_u4 = pin->GetOrAddBoolean("rad_m1","implicit_gas_newton_u4",false);
  // implicit_opac_lag = none (default) | pred: pred evaluates the opacities
  // in passes 0 .. N-1 only (pass 0 at the predicted T of the step predictor, or T^n on a
  // cold start; N = implicit_opac_lag_passes, default 1) and freezes them (and drops the
  // opac_newton face terms) for the remaining passes -- a change of the discretisation,
  // not only of the iteration (N = 2 keeps order 2 also without a predictor).
  {const std::string ol = pin->GetOrAddString("rad_m1","implicit_opac_lag","none");
  if (ol.compare("none") == 0) {
    impl_opac_lagp = false;
  } else if (ol.compare("pred") == 0) {
    impl_opac_lagp = true;
  } else {
    ImplFatal("<rad_m1>/implicit_opac_lag = '" + ol + "' is not a choice (none | pred)");
  }
  // (r8 key set): DEFAULT 3 (was 1; recipe none keeps 1)
  impl_opac_lagn = pin->GetOrAddInteger("rad_m1","implicit_opac_lag_passes",
                     (recipe::GetMode(pin) == recipe::Mode::none) ? 1 : 3);
  if (impl_opac_lagn < 1) {ImplFatal("<rad_m1>/implicit_opac_lag_passes must be >= 1");}
  }
  // m1-positivity set: DEFAULT ON (implicit_g0_exchange = true,
  // implicit_g0_limit = 1, implicit_pos_gas = true, implicit_pos_floor = true and
  // implicit_opac_newton_guard_mode = 6 below; every he_star_m1 production ran this set;
  // He box / wedge IDENT, BSG column within the 1-ulp twin spread:
  // /viper/ptmp2/jinma/defaults_gate_1002/item2).  A restart whose file lacks a key
  // keeps its old default (off; guard mode 2 is recorded in every file) and says so once.
  // The resolved values are recorded.  implicit_pos_gas_frac: read only when named.
  {const bool rst = global_variable::restart_run;
  bool old_kept = false;
  auto pbool = [&](const char *key, bool &v) {
    if (pin->DoesParameterExist("rad_m1",key)) {
      v = pin->GetBoolean("rad_m1",key);
    } else if (rst) {
      old_kept = true;
    } else {
      v = pin->GetOrAddBoolean("rad_m1",key,true);
    }
  };
  if (pin->DoesParameterExist("rad_m1","implicit_g0_limit")) {
    impl_g0_lim = pin->GetReal("rad_m1","implicit_g0_limit");
  } else if (rst) {
    old_kept = true;
  } else {
    impl_g0_lim = pin->GetOrAddReal("rad_m1","implicit_g0_limit",1.0);
  }
  pbool("implicit_g0_exchange", impl_g0_exch);
  pbool("implicit_pos_gas", impl_pos_gas);
  pbool("implicit_pos_floor", impl_pos_floor);
  // (user 10-09): default = implicit_pos_floor on a fresh run (it needs
  // it); a restart whose file lacks the key keeps false
  // (c): a restart whose file lacks the key is no longer silently false: a
  // WARN line says so, and <rad_m1>/implicit_pos_floor_solve_rst_default = true (read
  // only when named) gives such a restart the fresh-run default instead
  {const bool pfsn = pin->DoesParameterExist("rad_m1","implicit_pos_floor_solve");
  bool pfs_rd = false;
  if (!pfsn && rst && impl_pos_floor) {
    // prod-1112-cand: DEFAULT true (name implicit_pos_floor_solve_rst_default = false for
    // the prod-1111 behaviour: such a restart keeps false)
    pfs_rd = pin->DoesParameterExist("rad_m1","implicit_pos_floor_solve_rst_default") ?
             pin->GetBoolean("rad_m1","implicit_pos_floor_solve_rst_default") : true;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> WARN: restart input lacks implicit_pos_floor_solve: "
                << (pfs_rd ? "true (implicit_pos_floor_solve_rst_default)"
                           : "false (the e_floor raise of the accept is NOT "
                             "charged; name the key, or "
                             "implicit_pos_floor_solve_rst_default = true)") << std::endl;
    }
  }
  impl_pos_floor_s2 = pin->GetOrAddBoolean("rad_m1","implicit_pos_floor_solve",
                                           impl_pos_floor && (!rst || pfs_rd));
  }
  // (a) implicit_cons_final, dbg_cons_ledger (rad_m1.hpp): read only when
  // named; the solver checks are made at the first solve
  // prod-1112-cand: implicit_cons_final DEFAULT true (name false for the
  // behaviour, bitwise).  Not recorded when absent: a deck that cannot take it (no
  // implicit_pos_floor here; line_jacobi multi-D or polar_ring_avg at the first solve)
  // falls back to off with a one-line notice; only a NAMED true is fatal there.
  cf_named = pin->DoesParameterExist("rad_m1","implicit_cons_final");
  impl_cons_final = cf_named ? pin->GetBoolean("rad_m1","implicit_cons_final") : true;
  impl_lin_tolcap = pin->DoesParameterExist("rad_m1","implicit_lin_tolcap") ?
                    pin->GetReal("rad_m1","implicit_lin_tolcap") : 0.0;
  if (impl_lin_tolcap < 0.0) {ImplFatal("<rad_m1>/implicit_lin_tolcap must be >= 0");}
  impl_lin_tcstage = pin->DoesParameterExist("rad_m1","implicit_lin_tolcap_stage") ?
                     pin->GetInteger("rad_m1","implicit_lin_tolcap_stage") : 0;
  impl_lin_sktrue = pin->DoesParameterExist("rad_m1","implicit_lin_skip_true") ?
                    pin->GetReal("rad_m1","implicit_lin_skip_true") : 0.0;
  impl_lin_tclast = pin->DoesParameterExist("rad_m1","implicit_lin_tolcap_last") ?
                    pin->GetReal("rad_m1","implicit_lin_tolcap_last") : 0.0;
  dbg_cons_ledger = pin->DoesParameterExist("rad_m1","dbg_cons_ledger") &&
                    pin->GetBoolean("rad_m1","dbg_cons_ledger");
  if (impl_cons_final && !impl_pos_floor && !cf_named) {
    impl_cons_final = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_cons_final (default) off: implicit_pos_floor = "
                << "false" << std::endl;
    }
  }
  if (impl_cons_final && !impl_pos_floor) {
    ImplFatal("<rad_m1>/implicit_cons_final needs implicit_pos_floor = true (an E below "
              "e_floor is raised with the energy taken from the gas)");
  }
  if (impl_pos_floor_s2 && !impl_pos_floor) {
    ImplFatal("<rad_m1>/implicit_pos_floor_solve needs implicit_pos_floor = true");
  }
  if (old_kept && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> restart input lacks m1-positivity keys (implicit_g0_exchange,"
              << " g0_limit, pos_gas, pos_floor): keeping the old default (off) for "
              << "those; name them to switch" << std::endl;
  }
  }  // end of the restart-default scope
  if (pin->DoesParameterExist("rad_m1","implicit_pos_gas_frac")) {
    impl_pos_gas_frac = pin->GetReal("rad_m1","implicit_pos_gas_frac");
  }
  Kokkos::realloc(pos_cnt_d, M1_POS_N);
  Kokkos::deep_copy(pos_cnt_d, 0.0);
  // implicit_wb_vkick_reject (rad_m1.hpp; read only when named: the old scheme's, which
  // scheme = clean never reads)
  if (pin->DoesParameterExist("rad_m1","implicit_wb_vkick_reject")) {
    impl_wb_vkrej = pin->GetReal("rad_m1","implicit_wb_vkick_reject");
  }
  if (!(impl_wb_vkrej >= 0.0)) {
    ImplFatal("<rad_m1>/implicit_wb_vkick_reject must be >= 0 (0 = off)");
  }
  if (impl_wb_vkrej > 0.0) {
    if (pin->DoesParameterExist("rad_m1","implicit_wb_vkick_be_retry")) {
      impl_wb_vkbr = pin->GetBoolean("rad_m1","implicit_wb_vkick_be_retry");
    }
    Kokkos::realloc(vkick_cnt_d, 2);
    Kokkos::deep_copy(vkick_cnt_d, 0.0);
  }
  // implicit_am_cs_kecap (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_am_cs_kecap")) {
    impl_am_kecap = pin->GetReal("rad_m1","implicit_am_cs_kecap");
    if (!(impl_am_kecap >= 0.0 && impl_am_kecap <= 1.0)) {
      ImplFatal("<rad_m1>/implicit_am_cs_kecap must be in [0, 1] (0 = off)");
    }
    if (impl_am_kecap > 0.0) {
      Kokkos::realloc(amkc_d, 4);
      Kokkos::deep_copy(amkc_d, 0.0);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_am_cs_kecap=" << impl_am_kecap << ": the am_cs "
                  << "deposit of a cell gives the gas at most this fraction of its e_int "
                  << "as kinetic energy (state and stored slope)" << std::endl;
      }
    }
  }
  // amcsE-1111 implicit_am_cs_energy / _frac / implicit_am_cs_cons_diag (rad_m1.hpp)
  if (pin->DoesParameterExist("rad_m1","implicit_am_cs_energy")) {
    const std::string ae = pin->GetString("rad_m1","implicit_am_cs_energy");
    if (ae == "gas") {
      impl_am_en = 0;
    } else if (ae == "rad") {
      impl_am_en = 1;
    } else {
      ImplFatal("<rad_m1>/implicit_am_cs_energy must be gas or rad");
    }
    if (pin->DoesParameterExist("rad_m1","implicit_am_cs_energy_frac")) {
      impl_am_en_frac = pin->GetReal("rad_m1","implicit_am_cs_energy_frac");
      if (!(impl_am_en_frac >= 0.0 && impl_am_en_frac < 1.0)) {
        ImplFatal("<rad_m1>/implicit_am_cs_energy_frac must be in [0, 1)");
      }
    }
    if (impl_am_en == 1) {
      Kokkos::realloc(amen_d, 4);
      Kokkos::deep_copy(amen_d, 0.0);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_am_cs_energy=rad: the kinetic energy of the "
                  << "am_cs "
                  << "deposit is paid by the cell's radiation energy (at most "
                  << impl_am_en_frac << " of E; the rest by e_int)" << std::endl;
      }
    }
  }
  if (pin->DoesParameterExist("rad_m1","implicit_am_cs_cons_diag")) {
    impl_am_cons_diag = pin->GetInteger("rad_m1","implicit_am_cs_cons_diag");
    if (impl_am_cons_diag > 0) {Kokkos::realloc(amcd_d, 24);}
  }
  // implicit_am_cs_fold (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_am_cs_fold")) {
    const std::string af = pin->GetString("rad_m1","implicit_am_cs_fold");
    if (af == "applied") {
      impl_am_fold_ap = true;
    } else if (af != "stored") {
      ImplFatal("<rad_m1>/implicit_am_cs_fold must be stored or applied");
    }
    if (impl_am_fold_ap && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_am_cs_fold=applied: am_cs measures the gas+rad L "
                << "change from the write-back's F^old (cmem: face means of F_c + A E^n)"
                << std::endl;
    }
  }
  // implicit_am_cs_deposit (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_am_cs_deposit")) {
    const std::string ad = pin->GetString("rad_m1","implicit_am_cs_deposit");
    if (ad == "local") {
      impl_am_dep = 0;
    } else if (ad == "patch") {
      impl_am_dep = 1;
    } else {
      ImplFatal("<rad_m1>/implicit_am_cs_deposit must be local or patch");
    }
    if (impl_am_dep == 1) {
      if (!pmy_pack->pmesh->use_cubed_sphere || !pmy_pack->pmesh->three_d) {
        ImplFatal("<rad_m1>/implicit_am_cs_deposit = patch needs a 3-D cubed sphere");
      }
      Kokkos::realloc(amvp_d, 8);
      Kokkos::deep_copy(amvp_d, 0.0);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_am_cs_deposit=patch: the am_cs L vector is "
                  << "deposited as the minimum-KE rigid rotation of the lateral 3x3 "
                  << "patch; capped remainders go to the radial column" << std::endl;
      }
    }
  }
  // implicit_clip_kecap (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_clip_kecap")) {
    impl_clip_kecap = pin->GetReal("rad_m1","implicit_clip_kecap");
    if (!(impl_clip_kecap >= 0.0 && impl_clip_kecap <= 1.0)) {
      ImplFatal("<rad_m1>/implicit_clip_kecap must be in [0, 1] (0 = off)");
    }
    if (impl_clip_kecap > 0.0) {
      Kokkos::realloc(clkc_d, 3);
      Kokkos::deep_copy(clkc_d, 0.0);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_clip_kecap=" << impl_clip_kecap << ": the "
                  << "conservative write-back's F-clip momentum gives the gas at most "
                  << "this "
                  << "fraction of its e_int as kinetic energy (cs / sp volumes only)"
                  << std::endl;
      }
    }
  }
  // implicit_wb_vcap (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_wb_vcap")) {
    impl_wb_vcap = pin->GetReal("rad_m1","implicit_wb_vcap");
    if (!(impl_wb_vcap >= 0.0)) {
      ImplFatal("<rad_m1>/implicit_wb_vcap must be >= 0 (0 = off)");
    }
    if (impl_wb_vcap > 0.0) {
      Kokkos::realloc(vcap_cnt_d, 4);
      Kokkos::deep_copy(vcap_cnt_d, 0.0);
    }
    if (pin->DoesParameterExist("rad_m1","implicit_wb_vcap_measure")) {
      impl_wb_vcap_meas = pin->GetBoolean("rad_m1","implicit_wb_vcap_measure");
    }
    if (pin->DoesParameterExist("rad_m1","implicit_wb_vcap_print")) {
      impl_wb_vcap_np = pin->GetInteger("rad_m1","implicit_wb_vcap_print");
    }
  }
  // implicit_opac_newton_guard_mode: 0 = the diagonal test alone; + 1 = also the
  // right-hand side, + 2 = also the sign of the neighbour entry (M1OpnGuardFace; He presn
  // wedge gate: 7 floor clips with the diagonal alone, 0 with both), + 4 (m1-positivity)
  // = also the row diagonal dominance (m1_impl_asm).  DEFAULT 6 (the
  // positivity set above); 2 on a restart whose file lacks the key (files written since
  // 09-29 record it).
  impl_opn_guard_mode = pin->GetOrAddInteger("rad_m1","implicit_opac_newton_guard_mode",
                                             global_variable::restart_run ? 2 : 6);
  // implicit_opac_newton_slope_off / _slope_max (read only when
  // named; 0 = off).  At an opacity cliff (He recombination: |d ln kappa/d ln T| ~ 5-14
  // in the TOPS tables) the d(rho kappa_T)/dT term of a strongly coupled cell makes the E
  // iterate cycle (He giant N897 photosphere, HEGIANT_PICARD.md: 77 passes, NON-CONV);
  // these keys drop or damp the term in cliff-steep cells only.
  if (pin->DoesParameterExist("rad_m1","implicit_opac_newton_slope_off")) {
    impl_opn_soff = pin->GetReal("rad_m1","implicit_opac_newton_slope_off");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_opac_newton_slope_max")) {
    impl_opn_smax = pin->GetReal("rad_m1","implicit_opac_newton_slope_max");
  }
  impl_allow_multid = pin->GetOrAddBoolean("rad_m1","implicit_allow_multid",false);
  marshak_q = pin->GetOrAddReal("rad_m1","marshak_q",0.5);
  // implicit_marshak_face (m1-sp-order2, tests_m1/runs_5o_sporder2).  DEFAULT linear on
  // the spherical-polar wedge since m1-sp-order2b (tests_m1/runs_5q_sporder2b) with the
  // fixed-tensor closures (eddington, vet_sc, tau, vet_col); cell with the lagged
  // closures (m1, minerbo, kershaw: the runs_5h pp_np atmosphere, closure_lag = step at
  // implicit_cfl 1e6, went from 30 round-off NON-CONVERGED solves to 54 with Picard
  // blow-ups).  A restart whose file lacks the key keeps cell; the resolved value is
  // echoed.  Elsewhere read only when named (parameter dump and the cell form kept).
  if (sph_geom || pin->DoesParameterExist("rad_m1","implicit_marshak_face")) {
    const bool mfdef = !global_variable::restart_run &&
                       (eddington || vet_sc || tau_closure);
    const std::string smf = sph_geom ?
        pin->GetOrAddString("rad_m1","implicit_marshak_face", mfdef ? "linear" : "cell") :
        pin->GetString("rad_m1","implicit_marshak_face");
    if (smf.compare("cell") == 0) {
      impl_mface_lin = false;
    } else if (smf.compare("linear") == 0) {
      impl_mface_lin = true;
    } else {
      ImplFatal("<rad_m1>/implicit_marshak_face = '" + smf
                + "' is not a choice (cell | linear)");
    }
    if (impl_mface_lin && !sph_geom) {
      ImplFatal("<rad_m1>/implicit_marshak_face = linear is implemented on the "
                "spherical-polar wedge only");
    }
  }
  // implicit_face_weight (see rad_m1.hpp): distance (default on the spherical-polar
  // wedge) | equal (default elsewhere; the only choice on a Cartesian grid, which is
  // uniform).  On the wedge the key is recorded (GetOrAdd), so a restart keeps what the
  // run started with.  A RESTART whose embedded input predates the key keeps `equal`,
  // the behaviour it was run with, unless the key is given on the command line or in
  // an -i overlay (ModifyFromCmdline adds absent keys on a restart).  On a Cartesian grid
  // it is read only when named, so the parameter dump there is unchanged.
  if (sph_geom || pin->DoesParameterExist("rad_m1","implicit_face_weight")) {
    const bool named = pin->DoesParameterExist("rad_m1","implicit_face_weight");
    const std::string sfw = pin->GetOrAddString("rad_m1","implicit_face_weight",
        global_variable::restart_run ? "equal" : "distance");
    if (sfw.compare("equal") == 0) {
      impl_face_wdist = false;
    } else if (sfw.compare("distance") == 0) {
      impl_face_wdist = true;
    } else {
      ImplFatal("<rad_m1>/implicit_face_weight = '" + sfw
                + "' is not a choice (equal | distance)");
    }
    if (impl_face_wdist && !sph_geom) {
      ImplFatal("<rad_m1>/implicit_face_weight = distance is implemented on the "
                "spherical-polar wedge only (a Cartesian grid is uniform)");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1>/implicit_face_weight = " << sfw
                << (named ? " (from the input)"
                    : (global_variable::restart_run
                       ? " (restart file predates the key: kept the old equal)"
                       : " (default)")) << std::endl;
    }
  }
  // (scheme = clean: not read; the clean solve has its own face, and CleanPreflight
  // lists the keys an input names)
  std::string srn = pin->GetOrAddString("rad_m1","implicit_recon","dc");
  if (srn.compare("dc") == 0) {
    impl_recon = M1_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = M1_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_m1>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // implicit_enthalpy (see rad_m1_implicit.hpp).  DEFAULT plm since m1-defaults for the
  // closures fixed within a step (eddington, vet_sc, tau; tests_m1/runs_3s_space2: second
  // order, fewer Picard passes, NON-CONVERGED 0); m1 / minerbo / kershaw keep upwind
  // (not gated; tests_m1/runs_4b_defaults).  Up to d0c59f7c the key was read only
  // when named, so a restart file written then carries it only if the input named it;
  // such a restart that does NOT carry it keeps the old default, upwind.  The resolved
  // value is echoed, so later restarts keep it.
  const bool enth_up = global_variable::restart_run || !fixcl;
  const std::string sen = pin->GetOrAddString("rad_m1","implicit_enthalpy",
                                              enth_up ? "upwind" : "plm");
  if (sen.compare("upwind") == 0) {
    impl_enth = M1_IENTH_UPWIND;
  } else if (sen.compare("central") == 0) {
    impl_enth = M1_IENTH_CENTRAL;
  } else if (sen.compare("plm") == 0) {
    impl_enth = M1_IENTH_PLM;
  } else if (sen.compare("plmc") == 0) {   // diagnostic (rad_m1_implicit.hpp)
    impl_enth = M1_IENTH_PLMC;
  } else if (sen.compare("lin") == 0) {    // diagnostic (rad_m1_implicit.hpp)
    impl_enth = M1_IENTH_LIN;
  } else {
    ImplFatal("<rad_m1>/implicit_enthalpy = '" + sen
              + "' is not a choice (upwind | central | plm)");
  }
  // implicit_vimp (rad_m1_implicit.hpp).  DEFAULT on since m1-defaults2
  // (tests_m1/runs_5g_defaults2) only where the resolved time_scheme is hesdirk2, with
  // which it was validated (under be it over-damps P = 100 waves), and where it is valid:
  // a multi-D mesh with implicit_solver = bicgstab, nghost >= 2, hydro with coupling,
  // gas_feedback and dbg_gas_force.  Otherwise off, silently.  A restart whose file lacks
  // the key (read only when named before) keeps off; the resolved value is echoed;
  // explicit input overrides.  only where the gas is LIVE, <time>/
  // evolution = dynamic (a kinematic or static medium never takes the solve's radiative
  // kick, so the implicit velocity has nothing to predict; the validated test decks name
  // false there).  Before defaults-1012 such a fresh deck took true.
  {auto *pmh = pmy_pack->pmesh;
  const bool gas_live = !pin->DoesParameterExist("time","evolution") ||
                        (pin->GetString("time","evolution").compare("dynamic") == 0);
  const bool vdef = ((pin->GetString("rad_m1","time_scheme").compare("hesdirk2") == 0) ||
                     t2_beonly) &&
                    full && pmh->multi_d && (impl_solver == M1_ISOLV_BICGSTAB) &&
                    (pmh->mb_indcs.ng >= 2) && fl_on && gas_live &&
                    coupling && gas_feedback && dbg_gas_force &&
                    !global_variable::restart_run;
  impl_vimp = pin->GetOrAddBoolean("rad_m1","implicit_vimp",vdef);}
  // <rad_m1>/doppler_vel = stage | vstar (rad_m1.hpp dop_vstar)
  if (pin->DoesParameterExist("rad_m1","doppler_vel")) {
    const std::string dv = pin->GetString("rad_m1","doppler_vel");
    if (dv == "vstar") {
      dop_vstar = true;
    } else if (dv == "vstar_n") {
      dop_vstar = true;
      dop_vstar_n = true;
    } else if (dv == "vstar_nbe") {
      dop_vstar = true;
      dop_vstar_n = true;
      dop_vstar_nbe = true;
    } else if (dv != "stage") {
      ImplFatal("<rad_m1>/doppler_vel must be stage | vstar | vstar_n | vstar_nbe");
    }
    if (dop_vstar) {
      if (impl_vimp && pin->DoesParameterExist("rad_m1","implicit_vimp") &&
          global_variable::my_rank == 0) {
        std::cout << "<rad_m1> doppler_vel = vstar: implicit_vimp = true overridden to "
                     "false"
                  << std::endl;
      }
      impl_vimp = false;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> doppler_vel = vstar: v-terms at the "
                     "Heun "
                  << "stage velocity before t2inc, implicit_vimp off" << std::endl;
      }
    }
  }
  // <rad_m1>/doppler_face = donor | central (rad_m1.hpp dop_face; read only
  // when named, default donor = bitwise)
  if (pin->DoesParameterExist("rad_m1","doppler_face")) {
    const std::string df = pin->GetString("rad_m1","doppler_face");
    if (df == "central") {
      dop_face = 1;
    } else if (df != "donor") {
      ImplFatal("<rad_m1>/doppler_face must be donor | central");
    }
    if (dop_face == 1) {
      if (impl_vimp && global_variable::my_rank == 0) {
        std::cout << "<rad_m1> doppler_face = central: implicit_vimp overridden to false "
                  << "(its Newton rows linearise the donor flux)" << std::endl;
      }
      impl_vimp = false;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> doppler_face = central: x1 enthalpy flux "
                  << "a E = v (1 + chi) E as the face mean 0.5 (a_L E_L + a_R E_R) "
                     "(Athena++ "
                  << "DO near-central face), no donor / plm correction" << std::endl;
      }
    }
  }
  {  // (USER defaults group 1): DEFAULT true; inert unless a pgen marks hold_m
    impl_hold = pin->GetOrAddBoolean("rad_m1","implicit_hold",true);
    if (impl_hold) {
      auto &idh = pmy_pack->pmesh->mb_indcs;
      const int nh1 = idh.nx1 + 2*idh.ng;
      const int nh2 = (idh.nx2 > 1) ? (idh.nx2 + 2*idh.ng) : 1;
      const int nh3 = (idh.nx3 > 1) ? (idh.nx3 + 2*idh.ng) : 1;
      const int nmh = std::max(pmy_pack->nmb_thispack, pmy_pack->pmesh->nmb_maxperrank);
      Kokkos::realloc(hold_m, nmh, nh3, nh2, nh1);
      Kokkos::deep_copy(hold_m, 0.0);
    }
  }
  if (pin->DoesParameterExist("rad_m1","implicit_vimp_kill")) {
    vimp_kill = pin->GetBoolean("rad_m1","implicit_vimp_kill");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_vimp_lag_local")) {
    vimp_lag_local = pin->GetBoolean("rad_m1","implicit_vimp_lag_local");
  }
  if (pin->DoesParameterExist("rad_m1","time2_dtr_poslim")) {
    t2_dtr_pl = pin->GetBoolean("rad_m1","time2_dtr_poslim");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_pos_transport")) {
    const std::string pt = pin->GetString("rad_m1","implicit_pos_transport");
    if (pt == "none") {
      impl_pos_trans = 0;
    } else if (pt == "upwind") {
      impl_pos_trans = 1;
    } else {
      ImplFatal("<rad_m1>/implicit_pos_transport must be none | upwind");
    }
    if (pin->DoesParameterExist("rad_m1","implicit_pos_transport_mem")) {
      pt_gam = pin->GetReal("rad_m1","implicit_pos_transport_mem");
    }
    if (impl_pos_trans != 0) {
      pt_cnt = DvceArray1D<Real>("m1_pt_cnt", 2);
      Kokkos::deep_copy(pt_cnt, 0.0);
    }
  }
  if (pin->DoesParameterExist("rad_m1","implicit_vimp_ewf")) {
    impl_vimp_ewf = pin->GetBoolean("rad_m1","implicit_vimp_ewf");
  }
  // DIAGNOSTIC: a scale of the Jacobian P (1 = Newton).  The converged state does not
  // depend on it; only the Picard contraction does.
  impl_vimp_jscale = 1.0;
  if (pin->DoesParameterExist("rad_m1","implicit_vimp_jscale")) {
    impl_vimp_jscale = pin->GetReal("rad_m1","implicit_vimp_jscale");
  }
  // ---- THE runs_4a_accel LEVERS (tests_m1/runs_4a_accel, tests_m1/runs_4j_accmerge).
  // DEFAULT ON since m1-accmerge for time_scheme = be (and since m1-h2fast for hesdirk2,
  // tests_m1/runs_5f_h2fast, whose stage solves also take time2_one_pass_safety and
  // time2_lin_tol_fac, rad_m1_time2.cpp) with transport = implicit and a
  // closure whose tensor is fixed within a step (eddington, vet_sc, tau), wherever each
  // is valid:
  //   implicit_fast_kernels = true   (bitwise-exact kernel shortcuts)
  //   implicit_vimp_fold = true      where implicit_vimp and implicit_op_stencil are on
  //   implicit_one_pass = 8          where the predictor is on
  //   implicit_predictor_order = 2   where the predictor is on
  // (1 GPU, 3-D He box: be 51.6 -> 43.5 ms/cycle Eddington, 57.5 -> 50.7 vet_sc with
  // one_pass = 4; one_pass = 8 is the faster of the two in runs_4a_accel).  Up to the
  // merge of m1-accel the keys were read only when named, so a restart file written
  // before this default carries them only if its input named them; such a restart that
  // does NOT carry a key keeps the old value (off / 0 / 1).  The resolved values are
  // always echoed, so later restarts keep them.  Explicit input values override.
  //  implicit_vimp_fold: fold the implicit_vimp operator part into the stored stencil
  //    (x2/x3 +-1 into slots 3-6, +-2 neighbours into slots 19-24) instead of
  //    M1VimpRow per apply (round-off).  Needs implicit_op_stencil.
  //  implicit_one_pass = N: see ImplicitSolve; its state travels in M1ONEP01.
  //  implicit_predictor_order = 2: the predictor extrapolates the increment rate g =
  //    dE/dt linearly in time, g* = g1 + h dt1 with h = (g1 - g2)/dt2, h stored per cell
  //    in ipred channels 3 (E) and 4 (T).  Under hesdirk2 the backward-Euler steps (first
  //    step, fallbacks) then leave ipred alone, so that it holds stage-1 increments only.
  //  implicit_fast_kernels: closure = eddington: D_ab = 0 off the diagonal, so the
  //    stencil build and the right-hand side skip the off-diagonal (od) terms, which are
  //    exactly zero; hesdirk2 + vet_sc without extrapolation (time2_vet_extrap = false):
  //    the tensor save of Time2VetExtrapolate is a plain copy.
  {
  // m1-h2fast (tests_m1/runs_5f_h2fast): the same defaults for time_scheme = hesdirk2,
  // whose stage solves already track the one_pass q and the predictor per solve kind
  bool ts_ok = true;
  if (pin->DoesParameterExist("rad_m1","time_scheme")) {
    const std::string tsn = pin->GetString("rad_m1","time_scheme");
    ts_ok = (tsn.compare("be") == 0) || (tsn.compare("hesdirk2") == 0);
  }
  const bool ldef = full && ts_ok && fixcl && !global_variable::restart_run;
  impl_fastk = pin->GetOrAddBoolean("rad_m1","implicit_fast_kernels",ldef);
  // implicit_vimp_fold: default OFF since audit-m1-0928 (gain 0.9-2.3 % per simulated
  // second on the moving He box and sph_wedge, MI300A 1 and 2 GPUs; under the ~3 % bar
  // for a result-changing default).  Opt-in by name.
  impl_vfold = pin->GetOrAddBoolean("rad_m1","implicit_vimp_fold",false);
  // implicit_one_pass: default 0 (off) everywhere since audit-m1-0928 (8 + auto gained
  // -2.8..+0.0 % on the moving He box and sph_wedge on MI300A, i.e. nothing, and changes
  // results 5-22x the round-off spread; before, 0 only with force_reference = wb_arad,
  // m1-keydefault-0927).  Opt-in by name (e.g. 8).
  impl_onep = pin->GetOrAddInteger("rad_m1","implicit_one_pass",0);
  if (impl_onep != 0 && impl_onep < 2) {
    ImplFatal("<rad_m1>/implicit_one_pass (the check period) must be 0 (off) or >= 2");
  }
  if (impl_onep > 0) {
    impl_onep_s = pin->GetOrAddReal("rad_m1","implicit_one_pass_safety",3.0);
    if (!(impl_onep_s >= 1.0)) {
      ImplFatal("<rad_m1>/implicit_one_pass_safety must be >= 1");
    }
    // implicit_one_pass_auto (m1-onepass-auto): see ImplicitSolve.  A kind of solve that
    // goes auto_window check periods (auto_window * implicit_one_pass eligible solves)
    // without one solve accepted after one pass is switched to one_pass = 0; after
    // auto_reprobe such solves it is switched on again for one check period.  Default
    // true; false on a restart whose file does not carry the key (written before it
    // existed), which then continues as before (see above).
    impl_onep_auto = pin->GetOrAddBoolean("rad_m1","implicit_one_pass_auto",
                                          !global_variable::restart_run);
    impl_onep_awin = pin->GetOrAddInteger("rad_m1","implicit_one_pass_auto_window",2);
    impl_onep_arep = pin->GetOrAddInteger("rad_m1","implicit_one_pass_auto_reprobe",64);
    if (impl_onep_awin < 1 || impl_onep_arep < 1) {
      ImplFatal("<rad_m1>/implicit_one_pass_auto_window and _auto_reprobe must be >= 1");
    }
  }
  impl_pord = pin->GetOrAddInteger("rad_m1","implicit_predictor_order",
                                   (ldef && impl_pred) ? 2 : 1);
  if (impl_pord != 1 && impl_pord != 2) {
    ImplFatal("<rad_m1>/implicit_predictor_order must be 1 or 2");
  }
  impl_odskip = impl_fastk && eddington;
  if (impl_vfold && !impl_stencil) {
    ImplFatal("<rad_m1>/implicit_vimp_fold needs implicit_op_stencil = true");
  }
  }
  // LIMIT 4 of the 3a findings is NOT implemented in 3a2: a column still has to live
  // inside one MeshBlock along x1 (the fatal below).  The option is parsed so that the
  // input files and the gate scripts can already name it, and `gather` fatals rather
  // than silently doing something else.
  std::string spt = pin->GetOrAddString("rad_m1","implicit_partition","none");
  if (spt.compare("none") == 0) {
    impl_part = M1_IPART_NONE;
  } else if (spt.compare("gather") == 0) {
    impl_part = M1_IPART_GATHER;
  } else {
    ImplFatal("<rad_m1>/implicit_partition = '" + spt
              + "' is not a choice (none | gather)");
  }
  impl_recon_w = pin->GetOrAddReal("rad_m1","implicit_recon_w",-1.0);
  impl_res_floor = pin->GetOrAddReal("rad_m1","implicit_res_floor",0.0);
  // implicit_gas_newton_switch = W (he_estall_1003).  DEFAULT 10 on fresh runs (R2
  // default, user 10-03; it was off by default in 69e8f149 only because the switch
  // with the mode-0 root find diverged the BSG envelope, bsg_1001/AB_SOLVER_1003.md,
  // rst 00008: W = 10 died at cycle 14575, resid 0.31; with implicit_tsolve_opac_mode
  // 1 it is robust there, TSOLVE_FIX_1003.md).  A restart whose file lacks the key
  // keeps the old behaviour (0 = off, bitwise the pre-key code) and says so; a named
  // value always wins.  Under implicit_gas_newton, a solve
  // whose Picard residual is detected STALLED drops the gas Newton update (and the
  // opac_newton face terms, which need it) for the rest of that solve: the gas T is the
  // bracketed root of the exact backward-Euler gas equation with kappa_P, kappa_E at the
  // trial T (the implicit_tsolve_opac root find, forced on) in every remaining pass.
  // Stall test at pass it (resid = the global MPI_MAX value, so every rank agrees):
  //   it + 1 >= implicit_gas_newton_switch_min (default 20), it + 1 > W,
  //   min(resid over the last W passes) >= implicit_stall_fac (default 0.5) x the
  //   minimum before that window, and resid rose at least twice inside the window
  //   (a slowly but monotonically contracting solve keeps Newton).
  // Solves that converge never switch.  Diagnosis (he_estall_1003): in hot (2-3 MK),
  // low-density void cells the Newton row omits d kappa_P/dT, the lagged-opacity
  // iteration has slope < -1 and the gas T flips between two values every pass (a
  // resid 1e-8..2e-6 2-cycle that ran ~29 % of the He presn wedge solves to
  // implicit_maxit); the root find with kappa(T) has no such cycle.  He wedge -30..-36 %
  // wall, accuracy within the noise spread; He box bitwise (never switches).
  if (pin->DoesParameterExist("rad_m1","implicit_gas_newton_switch")) {
    impl_gn_sw = pin->GetInteger("rad_m1","implicit_gas_newton_switch");
  } else if (global_variable::restart_run) {
    impl_gn_sw = 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_gas_newton_switch: keeping "
                << "the old default (off); set it to switch" << std::endl;
    }
  } else {
    impl_gn_sw = pin->GetOrAddInteger("rad_m1","implicit_gas_newton_switch",10);
  }
  if (impl_gn_sw > 0) {
    impl_gn_sw_min = pin->GetOrAddInteger("rad_m1","implicit_gas_newton_switch_min",20);
    impl_stall_fac = pin->GetOrAddReal("rad_m1","implicit_stall_fac",0.5);
  }
  if (impl_gn_sw < 0 || impl_gn_sw > 64 || (impl_gn_sw > 0 && impl_gn_sw < 3) ||
      impl_gn_sw_min < 1 || !(impl_stall_fac > 0.0 && impl_stall_fac <= 1.0)) {
    ImplFatal("<rad_m1>/implicit_gas_newton_switch must be 0 or lie in [3,64], "
              "implicit_gas_newton_switch_min >= 1 and implicit_stall_fac in (0,1]");
  }
  // on the multi-D wedge berthon | blend default to `step`.  The lagged
  // reduced flux of a cell is the MEAN of its two faces' upwind ratios, and in free
  // streaming the upwind berthon flux reproduces any lagged f (F_f = c f E_up), so
  // recomputing the coefficients every Picard pass only AVERAGES f between neighbours,
  // pass after pass: a neutral iteration in which the zero-flux inner face and the
  // Marshak face leak inward (measured, gate b at CFL 0.4: Picard never meets 1e-12
  // in 200 passes and f falls from 1 to 0.69 in one step).  With `step` the face
  // coefficients come from the start-of-step state and the row is linear in E.
  const bool splag = full && (sph_geom || impl_flux_all) &&
                     (impl_flux != M1_IFLUX_CENTRAL);
  std::string slg = pin->GetOrAddString("rad_m1","implicit_recon_lag",
                                        splag ? "step" : "picard");
  impl_recon_freeze = (slg.compare("step") == 0);
  if (pin->DoesParameterExist("rad_m1","implicit_lag_restage")) {
    impl_lag_restage = pin->GetBoolean("rad_m1","implicit_lag_restage");
  } else if (impl_aphll_mem) {
    impl_lag_restage = true;
  }
  // MILESTONE 3c: freeze the deferred correction, and with it the plm limiter's choice,
  // after this many Picard passes.  The 3a2 finding is that the limiter keeps switching
  // on a handful of cells and the iteration is a small limit cycle that never meets
  // implicit_tol, so plm_dc costs implicit_maxit passes per step; freezing the
  // correction after a few passes leaves an ordinary linear system to converge.
  // <= 0 (the default) never freezes, i.e. reproduces 3a2.
  impl_recon_npass = pin->GetOrAddInteger("rad_m1","implicit_recon_npass",-1);
  if (pin->DoesParameterExist("rad_m1","implicit_recon_dgpass")) {
    impl_recon_dgpass = pin->GetBoolean("rad_m1","implicit_recon_dgpass");
  }
  if (!impl_recon_freeze && slg.compare("picard") != 0) {
    ImplFatal("<rad_m1>/implicit_recon_lag = '" + slg
              + "' is not a choice (step | picard)");
  }
  // LIMIT 3 of the 3a findings: the DEFAULT is now `true`.  It is a bug fix, not a
  // tuning knob (the boundary face used to hand its whole momentum to one interior cell,
  // which gave that cell 1.5 face-shares of radiative force: He column bottom-cell |v1|
  // 10.3 -> 0.47 v_MLT).  The key is kept so that `false` reproduces runs_3a/RESULTS.txt.
  impl_bmom_half = pin->GetOrAddBoolean("rad_m1","implicit_bmom_half",true);
  // the moving-medium fixes (rad_m1.hpp); every one off reproduces a61c1e75
  // the same default on the cubed sphere (the transverse deposit and the
  // fixed point are the cs blocks of the write-back; cs-s4-1017 had false there)
  impl_rad_inertia = pin->GetOrAddBoolean("rad_m1","implicit_rad_inertia",true);
  // implicit_rad_inertia_solve (read only when named; rad_m1.hpp)
  if (pin->DoesParameterExist("rad_m1","implicit_rad_inertia_solve")) {
    impl_rinr_solve = pin->GetBoolean("rad_m1","implicit_rad_inertia_solve");
  }
  impl_nofb_work = pin->GetOrAddBoolean("rad_m1","implicit_nofb_work",true);
  // implicit_work_fix (read only when named; rad_m1.hpp; DEFAULT true since
  // USER 10-25; false = the old write-back)
  if (pin->DoesParameterExist("rad_m1","implicit_work_fix")) {
    impl_wk_fix = pin->GetBoolean("rad_m1","implicit_work_fix");
  }
  // cs_vertex_mono / cs_vertex_mono_n (read only when named; rad_m1.hpp;
  // DEFAULT true on the cubed sphere, USER 10-29)
  cs_vtx_mono = cs_geom;
  if (pin->DoesParameterExist("rad_m1","cs_vertex_mono")) {
    cs_vtx_mono = pin->GetBoolean("rad_m1","cs_vertex_mono") && cs_geom;
  }
  if (pin->DoesParameterExist("rad_m1","implicit_tstart_old")) {
    impl_tstart_old = pin->GetBoolean("rad_m1","implicit_tstart_old");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_predictor_ecap")) {
    impl_pred_ecap = pin->GetReal("rad_m1","implicit_predictor_ecap");
  }
  // (USER 11-08 "make some defaults", r8 key set): implicit_predictor_tcap
  // DEFAULT 1.1 (was 2) and cs_cross_mono DEFAULT true (was false; cubed sphere only, see
  // xmono); read only when named as before; recipe none keeps the old values
  const bool rcp_none8 = recipe::GetMode(pin) == recipe::Mode::none;
  impl_pred_tcap = rcp_none8 ? 2.0 : 1.1;
  cs_cross_mono = !rcp_none8;
  if (pin->DoesParameterExist("rad_m1","implicit_predictor_tcap")) {   // seamE-1106
    impl_pred_tcap = pin->GetReal("rad_m1","implicit_predictor_tcap");
  }
  if (pin->DoesParameterExist("rad_m1","cs_cross_mono")) {   // seamE-1106
    cs_cross_mono = pin->GetBoolean("rad_m1","cs_cross_mono");
  }
  if (cs_cross_mono && pin->DoesParameterExist("rad_m1","cs_cross_mono_k")) {
    cs_cross_mono_k = pin->GetReal("rad_m1","cs_cross_mono_k");
  }
  if (pin->DoesParameterExist("rad_m1","cs_vertex_mono_n")) {
    cs_vtx_n = pin->GetInteger("rad_m1","cs_vertex_mono_n");
    if (cs_vtx_n < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1>/cs_vertex_mono_n must be >= 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // implicit_pos_floor_cons (read only when named; rad_m1.hpp; DEFAULT true
  // (USER 10-29); false = the old floor that creates energy)
  if (pin->DoesParameterExist("rad_m1","implicit_pos_floor_cons")) {
    impl_pos_flr_cons = pin->GetBoolean("rad_m1","implicit_pos_floor_cons");
  }
  if (impl_pos_flr_cons && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_pos_floor_cons: the E floor's uncovered raise is "
                 "taken "
              << "from the radiation of the same column (no energy created)" << std::endl;
  }
  // implicit_pos_gas_cons = true and implicit_pos_floor_gas_keep = 1: the defaults of a
  // fresh clean run (v20.4, recorded; a restart file without them keeps the old off)
  const bool cnw = clean && !global_variable::restart_run;
  if (cnw || pin->DoesParameterExist("rad_m1","implicit_pos_gas_cons")) {
    impl_pos_gas_cons = pin->GetOrAddBoolean("rad_m1","implicit_pos_gas_cons",true);
  }
  if (cnw || pin->DoesParameterExist("rad_m1","implicit_pos_floor_gas_keep")) {
    impl_pos_flr_gk = pin->GetOrAddReal("rad_m1","implicit_pos_floor_gas_keep",1.0);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_pos_floor_gas_keep=" << impl_pos_flr_gk
                << ": the E floor's raise leaves the gas at least that share of its eint"
                << std::endl;
    }
  }
  // implicit_fclip_diag (rad_m1.hpp; read only when named)
  if (pin->DoesParameterExist("rad_m1","implicit_fclip_diag") &&
      pin->GetBoolean("rad_m1","implicit_fclip_diag")) {
    impl_fclip_diag = true;
    Kokkos::realloc(fcd_cnt_d, 20);
    Kokkos::deep_copy(fcd_cnt_d, 0.0);
  }
  if (cs_vtx_mono && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> cs_vertex_mono: two-point lateral faces within " << cs_vtx_n
              << " cell(s) of the 8 cube vertices (no skew cross term, seam foot values, "
              << "D_r,lat = 0 within " << cs_vtx_n + 1 << ")" << std::endl;
  }
  // the gas source of the coupled solve (rad_m1.hpp, NOTE-consrc.md).
  // (user sign-off of): DEFAULT conservative2 (with its
  // implicit_cs2_hll = drag, implicit_cs2_shift = false) on a FRESH run with gas coupling
  // (hydro, coupling, gas_feedback); matched without it and on a restart whose file
  // lacks the key (files echo it).  Explicit input overrides.
  {
    // the same default on the cubed sphere (the conservative sources are
    // ported: cs blocks of the write-back, face-normal A^old / dA, m1bcs, the paired dA
    // exchange; cs-s4-1017 had matched there)
    const bool gsdef = !global_variable::restart_run && fl_on && coupling &&
                       gas_feedback;
    const std::string gs = pin->GetOrAddString("rad_m1","implicit_gas_src",
                                               gsdef ? "conservative2" : "matched");
    if (gs.compare("matched") == 0) {
      impl_gas_src_cons = false;
    } else if (gs.compare("conservative") == 0) {
      impl_gas_src_cons = true;
    } else if (gs.compare("conservative_full") == 0) {
      // the AP-HLL / blend part of a face too (no drag-only exception): exact total
      // momentum on every face, at the price of the HLL flux's own (non-momentum-
      // equation) changes going to the gas
      impl_gas_src_cons = true;
      impl_gas_src_full = true;
    } else if (gs.compare("conservative2") == 0) {
      // conservative + the theta split of the whole dA + the vimp
      // prediction of the conservative velocity (rad_m1.hpp)
      impl_gas_src_cons = true;
      impl_gas_src_cs2 = true;
      cs2_dbg = pin->GetOrAddInteger("rad_m1","implicit_cs2_dbg",0);
      // implicit_cs2_vimp: the vimp prediction of the conservative
      // write-back velocity (ImplicitCs2Vimp, (b) in rad_m1.hpp).  DEFAULT true (FIN3);
      // false (diagnostic, bitwise implicit_cs2_dbg = 2) drops it: rw_p10t1 gas order 2
      // but rw_p10t100 E 1.75 at l 16-32, so not adopted.
      cs2_vimp = pin->GetOrAddBoolean("rad_m1","implicit_cs2_vimp",true);
      // implicit_cs2_vimp_npass: the prediction lags the neighbours' dA
      // and the cross terms by one Picard pass, and the solves stop at 2 passes, so the
      // E rows' dv misses the write-back's by O(dt) relative (implicit_cs2_dbg bit 32:
      // 2.7e-3 / 1.4e-3 at rw_p10t1 l 1 / 2): a PER-STEP O(dt) term (rw_p10t1 hesdirk2 T
      // rho 1.58 at l 16-32, rw_p10t100 rho 1.29 at l 32-64).  At least npass Picard
      // passes while the prediction is active close it (npass 3: rw_p10t1 rho 1.98-2.0
      // to l 64, rw_p10t100 rho 2.03; AG Car E1 +9 %, BSG +4 % s/cycle on MI300A).
      // DEFAULT 3 on a fresh run; a restart whose file lacks the key keeps 0 (FIN3,
      // bitwise).  0 = off.
      cs2_vnp = pin->GetOrAddInteger("rad_m1","implicit_cs2_vimp_npass",
                                     global_variable::restart_run ? 0 : 3);
      if (cs2_vnp < 0) {
        ImplFatal("<rad_m1>/implicit_cs2_vimp_npass must be >= 0");
      }
      // consfix2: the hydro's dA as F0's explicit tendency (ImplicitCs2Shift) is opt-in;
      // off, it stays in the solve's dA (clean order 2 on rw_p10t100 / rsw / rlat)
      cs2_shift = pin->GetOrAddBoolean("rad_m1","implicit_cs2_shift",false);
      // implicit_cs2_clip (read only when named): gas (default: the gas is charged with
      // the realisability-limited cell flux, so the limiter's change of F -- up to |F| -
      // cE in free-streaming cells -- is a momentum kick to the gas: ~100 km/s in 1e-8
      // floor gas, v12lz accretor) | rad (the limit's change stays with the radiation, as
      // on the matched path) (user 11-12): implicit_cs2_wbsolve DEFAULT true (recorded);
      // name false for the 3 fixed-point passes of before
      cs2_wbs = pin->GetOrAddBoolean("rad_m1","implicit_cs2_wbsolve",true);
      if (pin->DoesParameterExist("rad_m1","implicit_cs2_clip")) {
        const std::string cm = pin->GetString("rad_m1","implicit_cs2_clip");
        if (cm.compare("rad") == 0) {
          cs2_clip_rad = true;
        } else if (cm.compare("gas") != 0) {
          ImplFatal("<rad_m1>/implicit_cs2_clip = '" + cm + "' is not a choice "
                    "(gas | rad)");
        }
      }
      // implicit_cs2_dacap (read only when named; rad_m1.hpp)
      if (pin->DoesParameterExist("rad_m1","implicit_cs2_dacap")) {
        cs2_cap = pin->GetReal("rad_m1","implicit_cs2_dacap");
        if (!(cs2_cap >= 0.0)) {
          ImplFatal("<rad_m1>/implicit_cs2_dacap must be >= 0 (0 = off)");
        }
        // the speed of the cap is at least vmin c (default 1e-4): a (near-)static cell
        // keeps the truncation-level dA of its enthalpy flux (read only with the cap on)
        if (cs2_cap > 0.0) {
          cs2_cap_vmin = pin->GetOrAddReal("rad_m1","implicit_cs2_dacap_vmin",1.0e-4);
          if (!(cs2_cap_vmin >= 0.0)) {
            ImplFatal("<rad_m1>/implicit_cs2_dacap_vmin must be >= 0");
          }
        }
      }
      // consfix3: the AP-HLL (thin) part of a face has no face equation, so the theta
      // split cannot act there; `drag` (default): that part charges the gas only its
      // drag, as matched (the radiation keeps its lab flux), `inertia`: all of dA (the
      // pre-consfix3 conservative2; floor gas with E/(rho c^2) >~ 1 then takes kicks ~c)
      {
        const std::string hm = pin->GetOrAddString("rad_m1","implicit_cs2_hll","drag");
        if (hm.compare("drag") == 0) {
          cs2_hll_drag = true;
        } else if (hm.compare("inertia") == 0) {
          cs2_hll_drag = false;
        } else {
          ImplFatal("<rad_m1>/implicit_cs2_hll = '" + hm + "' is not a choice (drag | "
                    "inertia)");
        }
      }
    } else {
      ImplFatal("<rad_m1>/implicit_gas_src = '" + gs + "' is not a choice (matched | "
                "conservative | conservative2 | conservative_full)");
    }
  }
  {
    const std::string fi = pin->GetOrAddString("rad_m1","implicit_f0_init","auto");
    if (fi.compare("zero") == 0) {
      impl_f0_init = 0;
    } else if (fi.compare("lab") == 0) {
      impl_f0_init = 1;
    } else if (fi.compare("auto") == 0) {
      impl_f0_init = 2;
    } else {
      ImplFatal("<rad_m1>/implicit_f0_init = '" + fi + "' is not a choice (auto | zero | "
                "lab)");
    }
  }
  impl_f0_init_done = false;
  if (impl_lin_maxit < 1) {
    ImplFatal("<rad_m1>/implicit_lin_maxit must be >= 1");
  }
  if (!(impl_tol > 0.0) || impl_maxit < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/implicit_tol must be positive and implicit_maxit >= 1"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &mbcs = pmy_pack->pmesh->mesh_bcs;
  ibc_x1min = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      mbcs[static_cast<int>(BoundaryFace::inner_x1)]);
  ibc_x1max = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1max","auto"),
      mbcs[static_cast<int>(BoundaryFace::outer_x1)]);
  if (!clean && (ibc_x1min == M1_IBC_OUTFLOW || ibc_x1max == M1_IBC_OUTFLOW)) {
    ImplFatal("<rad_m1>/implicit_bc_x1min|max = outflow needs scheme = clean");
  }
  // the transverse faces (read only when named; auto = the old rule)
  {
    auto tbc = [&](const char *key, const BoundaryFace f) -> int {
      const bool per = (mbcs[static_cast<int>(f)] == BoundaryFlag::periodic);
      std::string v = "auto";
      if (pin->DoesParameterExist("rad_m1", key)) {v = pin->GetString("rad_m1", key);}
      if (v.compare("auto") == 0) {return per ? M1_IBC_PERIODIC : M1_IBC_REFLECT;}
      int r = -1;
      if (v.compare("reflect") == 0) {r = M1_IBC_REFLECT;}
      if (v.compare("marshak") == 0) {r = M1_IBC_MARSHAK;}
      // zero-gradient outgoing flux, x2 faces only
      if (v.compare("outflow") == 0 && (f == BoundaryFace::inner_x2 ||
                                        f == BoundaryFace::outer_x2)) {
        r = M1_IBC_OUTFLOW;
      }
      if (v.compare("flux") == 0 && (f == BoundaryFace::inner_x2 ||
                                     f == BoundaryFace::outer_x2)) {r = M1_IBC_X2FLUX;}
      if (r < 0 || per) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "<rad_m1>/" << key << " = '" << v << "': auto | reflect | "
          << "marshak (| outflow | flux on x2), the last ones on a NON-periodic mesh face"
          << std::endl;
        std::exit(EXIT_FAILURE);
      }
      return r;
    };
    ibc_x2min = tbc("implicit_bc_x2min", BoundaryFace::inner_x2);
    ibc_x2max = tbc("implicit_bc_x2max", BoundaryFace::outer_x2);
    ibc_x3min = tbc("implicit_bc_x3min", BoundaryFace::inner_x3);
    ibc_x3max = tbc("implicit_bc_x3max", BoundaryFace::outer_x3);
  }
  iflux_x1min = pin->GetOrAddReal("rad_m1","implicit_flux_x1min",0.0);
  // (read only when named; test only): implicit_flux_x1min_mod_amp = A,
  // implicit_flux_x1min_mod_period = P > 0: the inner face flux F0 (1 + A sin(2 pi t/P)),
  // F0 = implicit_flux_x1min, set before each hesdirk2 stage solve at its time
  // (Time2FormStage: t_n + dt), so the forcing of a time-order test is time-consistent
  if (pin->DoesParameterExist("rad_m1","implicit_flux_x1min_mod_period")) {
    fmod_per = pin->GetReal("rad_m1","implicit_flux_x1min_mod_period");
    fmod_amp = pin->GetReal("rad_m1","implicit_flux_x1min_mod_amp");
    fmod_f0 = iflux_x1min;
  }
  // (read only when named; test only): implicit_flux_x1min_ramp_t0 = t0,
  // implicit_flux_x1min_ramp_w = w > 0: the inner flux times 0.5 (1 + tanh((t - t0)/w)),
  // a smooth switch-on (no start-up discontinuity), also at the stage times; the value at
  // the start time is set here
  if (pin->DoesParameterExist("rad_m1","implicit_flux_x1min_ramp_w")) {
    fmod_rw = pin->GetReal("rad_m1","implicit_flux_x1min_ramp_w");
    fmod_rt0 = pin->GetReal("rad_m1","implicit_flux_x1min_ramp_t0");
    if (fmod_per <= 0.0) {fmod_f0 = iflux_x1min;}
    if (fmod_rw > 0.0) {
      const Real t = pmy_pack->pmesh->time;
      iflux_x1min = fmod_f0*0.5*(1.0 + std::tanh((t - fmod_rt0)/fmod_rw));
    }
  }
  // (read only when named; absent = bitwise): the inner face flux from
  // a (t, F or L) table (rad_m1.hpp ftab_on), set here at the start time and at the
  // stage times by the hooks above
  if (pin->DoesParameterExist("rad_m1","implicit_flux_x1min_table")) {
    const std::string fn = pin->GetString("rad_m1","implicit_flux_x1min_table");
    const std::string err = ftab.Read(fn, 1);
    if (!err.empty()) {ImplFatal("<rad_m1>/implicit_flux_x1min_table: " + err);}
    if (ibc_x1min != M1_IBC_FLUX) {
      ImplFatal("<rad_m1>/implicit_flux_x1min_table needs implicit_bc_x1min = flux");
    }
    ftab_t0 = pin->GetOrAddReal("rad_m1","implicit_flux_x1min_table_t0",0.0);
    const std::string un = pin->GetOrAddString("rad_m1","implicit_flux_x1min_table_units",
                                               "flux");
    if (un == "flux") {
      ftab_scale = 1.0;
    } else if (un == "luminosity") {
      const Real r0 = pmy_pack->pmesh->mesh_size.x1min;
      ftab_scale = 1.0/(4.0*M_PI*r0*r0);
    } else {
      ImplFatal("<rad_m1>/implicit_flux_x1min_table_units must be flux|luminosity");
    }
    ftab_on = true;
    if (fmod_per <= 0.0 && fmod_rw <= 0.0) {fmod_f0 = iflux_x1min;}
    iflux_x1min = InnerFluxAt(pmy_pack->pmesh->time);
    if (global_variable::my_rank == 0) {
      std::cout << "rad_m1: implicit_flux_x1min_table " << fn << ": " << ftab.t.size()
                << " rows, t0 " << ftab_t0 << ", units " << un << ", F(t = "
                << pmy_pack->pmesh->time << ") = " << iflux_x1min << std::endl;
    }
  }
  iflux_x1max = pin->GetOrAddReal("rad_m1","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_m1","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_m1","implicit_ebath_x1max",0.0);
  if (pin->DoesParameterExist("rad_m1","athena_face_diss")) {   // col1d A7
    afd_on = pin->GetBoolean("rad_m1","athena_face_diss");
    if (afd_on) {
      afd_force = pin->GetOrAddBoolean("rad_m1","athena_face_diss_force",false);
      afd_tc = pin->GetOrAddReal("rad_m1","athena_face_taucell",4.0);
      afd_mu1 = pin->GetOrAddReal("rad_m1","athena_face_mu1",0.5063576644955016);
      afd_taumin = pin->GetOrAddReal("rad_m1","athena_face_taumin",0.0);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> athena_face_diss: Athena++ HLL tau-factor dissipation at "
                  << "every interior x1 face, taucell " << afd_tc << ", <mu>_+ "
                  << afd_mu1 << ", tau_cell min " << afd_taumin
                  << (afd_force ? ", in the face flux (force consistent)"
                                : ", energy row only (no matching force)") << std::endl;
      }
    }
  }
  if (pin->DoesParameterExist("rad_m1","implicit_ghost_x1min")) {
    impl_ghx = pin->GetBoolean("rad_m1","implicit_ghost_x1min");
  }
  if (impl_ghx) {
    impl_ghx_e = pin->GetReal("rad_m1","implicit_ghost_e_x1min");
    if (ibc_x1min != M1_IBC_FLUX || !sph_geom || !(impl_ghx_e > 0.0)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/implicit_ghost_x1min = true needs implicit_bc_x1min = "
        << "flux, spherical-polar geometry and implicit_ghost_e_x1min > 0" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // runs_4f_drift: without it the Marshak end face carries only the comoving flux
  // c q (E - E_bath), so the lab-frame radiation energy advected with the gas,
  // A E = v (1 + chi) E, cannot leave (or enter) through the end: the end cell's E is
  // raised by ~(A/(c q)) E at an outflow end (lowered at an inflow end), and the
  // Lowrie-Edwards shocks drift / carry an N-independent T error.  DEFAULT true since
  // m1-bcadv (tests_m1/runs_4h_bcadv).  On a restart whose file lacks the key (written
  // before 5c9e4432) the old value false is kept; files since then echo their value.
  impl_bc_advect = pin->GetOrAddBoolean("rad_m1","implicit_bc_advect",
                                        !global_variable::restart_run);
  // implicit_vimp_bc (ImplicitVimpBuild).  DEFAULT true since
  // where implicit_vimp is on (the only place it acts): a fresh run takes
  // true, a restart whose file lacks the key keeps the old value false; the resolved
  // value is echoed (later restarts keep it); explicit input overrides.  Without
  // implicit_vimp the key is read only when named (default false, not echoed).
  if (impl_vimp) {
    impl_vimp_bc = pin->GetOrAddBoolean("rad_m1","implicit_vimp_bc",
                                        !global_variable::restart_run);
  } else if (pin->DoesParameterExist("rad_m1","implicit_vimp_bc")) {
    impl_vimp_bc = pin->GetBoolean("rad_m1","implicit_vimp_bc");
  }
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- the restrictions of 3a
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mindcs = pmy_pack->pmesh->mesh_indcs;
  part_nblk = 1;
  if (!clean && mindcs.nx1 != indcs.nx1) {   // the clean solve has its own exchange
    if (impl_part != M1_IPART_GATHER) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/transport = implicit_x1 with implicit_partition = none "
        << "needs exactly ONE MeshBlock along x1: <meshblock>/nx1 must equal <mesh>/nx1 ("
        << indcs.nx1 << " vs " << mindcs.nx1 << ").  Set <rad_m1>/implicit_partition = "
        << "gather for the line solve partitioned over MeshBlocks and ranks" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((mindcs.nx1 % indcs.nx1) != 0) {
      ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nx1 to be an exact "
                "multiple of <meshblock>/nx1 (uniform mesh only)");
    }
    part_nblk = mindcs.nx1/indcs.nx1;
  }
  if (part_nblk > 1 && (ibc_x1min == M1_IBC_PERIODIC)) {
    ImplFatal("<rad_m1>/implicit_partition = gather does not support PERIODIC x1 across "
              "more than one MeshBlock (the cyclic Thomas sweep of 3a wraps inside one "
              "block).  Use one MeshBlock along x1, or a non-periodic x1 boundary pair");
  }
  if (part_nblk > 1 && indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nghost >= 2");
  }
  if (full) {
    trans_on = pmy_pack->pmesh->multi_d;
    trans_x3 = pmy_pack->pmesh->three_d;
    if (!trans_on) {
      // a 1-D mesh has no transverse direction: the solve IS the x1 column solve, and
      // gate G3 is exactly this statement.
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1>: transport=implicit on a 1-D mesh is the x1 column solve"
                  << std::endl;
      }
    }
  }
  // which faces carry the per-face beam-vector coefficients
  blat_on = full && trans_on && impl_flux_all;
  bvec_x1 = full && impl_flux_all && !sph_geom;
  if (impl_flux_all && full && (impl_recon == M1_IRECON_PLMDC)) {
    ImplFatal("<rad_m1>/implicit_flux_faces = all takes implicit_recon = dc only");
  }
  // (stage S4): halfrange on the cubed sphere too, with vet_gd (VetGdHalfRange
  // about the cs axes, VetGdCsHrPrep, ImplicitLatFaceCoefCs)
  if (impl_beam_hr && (part_nblk > 1 || impl_recon == M1_IRECON_PLMDC)) {
    ImplFatal("<rad_m1>/implicit_flux_beam = halfrange needs one MeshBlock along x1 and "
              "implicit_recon = dc");
  }
  if (impl_beam_hr && cs_geom && !vgd_on && !impl_hr_model) {
    ImplFatal("<rad_m1>/implicit_flux_beam = halfrange on the cubed sphere needs vet_gd "
              "(the ray data of the cs sweep) or implicit_hr_model");
  }
  if (cs_geom && (impl_beam_fs || impl_hr_sfgd != 0 || (impl_hr_sf && !vgd_on))) {
    ImplFatal("<rad_m1> cubed sphere: implicit_flux_beam = fs and implicit_hr_face = sf "
              "(implicit_hr_face_gd x1 | all) are not ported (the vet_sc / sp-wedge "
              "data); use implicit_flux_beam = halfrange with implicit_hr_face_gd = "
              "none");
  }
  if (impl_beam_fs && !(blat_on && (vet_sc || (vgd_on && sph_geom)))) {
    ImplFatal("<rad_m1>/implicit_flux_beam = fs needs implicit_flux_faces = all on a "
              "multi-D mesh and closure = vet_sc, or vet_gd on the sp wedge");
  }
  if (blat_on && part_nblk > 1) {
    ImplFatal("<rad_m1>/implicit_flux_faces = all needs one MeshBlock along x1");
  }
  if ((indcs.nx2 > 1 || indcs.nx3 > 1) && !impl_allow_multid && !full) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does NOT transport along x2/x3. "
      << "Set <rad_m1>/implicit_allow_multid = true to run a set of INDEPENDENT x1 "
      << "columns (F_2 = F_3 = 0 everywhere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does not support SMR/AMR"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // one solve per hydro step: no sub-cycling, one stage
  subcycle = false;
  nstage = 1;
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (impl_cfl > 0.0) || (!other_sets_dt);

  // ---- arrays
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(f0x1, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1, 0.0);
  // (hgmem: scheme = clean: f0x*n, gcz*, cs2f*, thin_frz and ktn are ImplicitSolve's,
  // not made; the keys that name them have no effect there)
  if (!clean) {Kokkos::realloc(f0x1n, nmb, ncells3, ncells2, ncells1+1);}
  Kokkos::deep_copy(f0x1n, 0.0);
  // without the half-range beam the stored flux has no blend: inert (echoed)
  // also the berthon | blend AP-HLL faces (the lagged reduced flux of
  // the face coefficients is then F_c/(c E_up), evolved by the face equation, instead of
  // the applied flux's own ratio, which has no F^n memory at w = 1: a per-STEP map)
  if (impl_hr_cmem && !impl_beam_hr && impl_flux == M1_IFLUX_CENTRAL) {
    impl_hr_cmem = false;
    if (impl_physflux) {
      ImplFatal("<rad_m1>/implicit_physflux needs implicit_flux_beam = halfrange or an "
                "AP-HLL face (it is cmem's two-flux split made physical)");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_hr_cmem: inert without implicit_flux_beam = "
                << "halfrange" << std::endl;
    }
  }
  if (impl_hr_cmem) {   // the applied face flux (F_c lives in f0x1)
    // implicit_vimp reads the APPLIED fluxes fax* (ImplicitVimpBuild)
    // the conservative gas source (Z_f) and implicit_rad_inertia take f0x* = F_c as the
    // face's old flux: exact without gas momentum exchange (static / no feedback
    // tests); with feedback the face momentum budget is off by (1 - w)(F_c^n - F_a^n)
    if (impl_physflux && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_physflux = true: F_c is the momentum flux of the "
                << "cell F, the drag and the gas source; fax* = E transport only"
                << std::endl;
    }
    if ((impl_gas_src_cons || impl_rad_inertia) && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_hr_cmem: WARNING the conservative gas source / "
                << "rad inertia read F_c as the old face flux (TEST key with feedback)"
                << std::endl;
    }
    Kokkos::realloc(fax1, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::deep_copy(fax1, 0.0);
    impl_cm_fainit = pin->DoesParameterExist("rad_m1","implicit_hr_cmem_fainit") &&
                     pin->GetBoolean("rad_m1","implicit_hr_cmem_fainit");
    impl_cm_rst = pin->DoesParameterExist("rad_m1","implicit_hr_cmem_rst") &&
                  pin->GetBoolean("rad_m1","implicit_hr_cmem_rst");
    // implicit_hr_cmem_gas_fa
    impl_cmgfa = pin->DoesParameterExist("rad_m1","implicit_hr_cmem_gas_fa") &&
                 pin->GetBoolean("rad_m1","implicit_hr_cmem_gas_fa");
    if (impl_pf_fa && !impl_pf_zfc) {impl_cmgfa = true;}   // mom = fa
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_hr_cmem = true: separate central-branch face "
                << "memory F_c (applied flux in fax*)" << std::endl;
    }
    // implicit_hr_cmem_implicit (rad_m1.hpp)
    impl_cmi = pin->DoesParameterExist("rad_m1","implicit_hr_cmem_implicit") &&
               pin->GetBoolean("rad_m1","implicit_hr_cmem_implicit");
    if (impl_cmi) {
      if (impl_pos_trans != 0) {
        ImplFatal("<rad_m1>/implicit_hr_cmem_implicit: not with implicit_pos_transport = "
                  "upwind (both rewrite the x1 memory term)");
      }
      cmi_cnt = DvceArray1D<Real>("m1_cmi_cnt", 3);
      Kokkos::deep_copy(cmi_cnt, 0.0);
      auto cch = Kokkos::create_mirror_view(cmi_cnt);
      cch(2) = 1.0;
      Kokkos::deep_copy(cmi_cnt, cch);
      if (trans_on) {
        Kokkos::realloc(cmm2, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(cmq2, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(cmd2, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(cmtr, nmb, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(cmm2, 0.0);
        Kokkos::deep_copy(cmq2, 0.0);
        Kokkos::deep_copy(cmd2, 0.0);
        Kokkos::deep_copy(cmtr, 0.0);
        if (trans_x3) {
          Kokkos::realloc(cmm3, nmb, ncells3+1, ncells2, ncells1);
          Kokkos::realloc(cmq3, nmb, ncells3+1, ncells2, ncells1);
          Kokkos::realloc(cmd3, nmb, ncells3+1, ncells2, ncells1);
          Kokkos::deep_copy(cmm3, 0.0);
          Kokkos::deep_copy(cmq3, 0.0);
          Kokkos::deep_copy(cmd3, 0.0);
        }
      }
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_hr_cmem_implicit = true: the F_c memory part of "
                  << "every face flux implicit upwind (q E'_d), its remainder Zalesak-"
                  << "limited by the donor E^k" << std::endl;
      }
    }
  }
  if (impl_gas_src_cons && !clean) {
    Kokkos::realloc(gcz1, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::deep_copy(gcz1, 0.0);
  }
  if ((impl_gas_src_cs2 || impl_rinr_solve) && !clean) {   // dA face means in F0^n
    Kokkos::realloc(cs2f1, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::deep_copy(cs2f1, 0.0);
    if (trans_on) {
      Kokkos::realloc(cs2f2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(cs2f2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(cs2f3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(cs2f3, 0.0);
      }
    }
  }
  // MILESTONE 3b phase C: the BiCGStab wrapper needs the four/six transverse off-diagonal
  // coefficients and ten Krylov vectors on top of the phase-B work array.  They are
  // allocated only when it is selected AND the mesh is multi-D, so line_jacobi keeps the
  // footprint (and, on a 1-D mesh, the solve IS the column solve and there is no system
  // to wrap).
  bicg_on = trans_on && (impl_solver == M1_ISOLV_BICGSTAB);
  // MILESTONE 3b phase D.  `auto` = the off-diagonal Eddington terms go INTO the operator
  // wherever there is a Krylov solver to carry them, and stay LAGGED otherwise, which is
  // what every line_jacobi and implicit_x1 configuration did before phase D.
  if (od_auto) {
    impl_offdiag = bicg_on ? M1_OD_OPERATOR : M1_OD_LAGGED;
  }
  if (impl_offdiag == M1_OD_OPERATOR && trans_on && !bicg_on) {
    ImplFatal("<rad_m1>/implicit_offdiag = operator needs implicit_solver = bicgstab: "
              "the off-diagonal Eddington terms are a 9-/19-point coupling and the "
              "x1 line solve of line_jacobi cannot carry them");
  }
  od_now = impl_offdiag;
  if (impl_vimp && !bicg_on) {
    ImplFatal("<rad_m1>/implicit_vimp needs transport = implicit on a multi-D mesh with "
              "implicit_solver = bicgstab: the implicit velocity couples cells two "
              "apart");
  }
  if (impl_vimp && pmy_pack->pmesh->mb_indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_vimp needs <mesh>/nghost >= 2");
  }
  int niw = full ? (bicg_on ? M1_NIW_K : M1_NIW) : M1_NIW_X1;
  // scheme = clean reads none of the components appended below, nor the Krylov
  // vectors (its iw ends at M1_NIW; time2_face_fct's scratch, Time2FctSlot)
  const int niw_clean = full ? M1_NIW : niw;
  // MILESTONE 3g: the five per-cell components of the gas coupling are APPENDED, so
  // every index above keeps the value it had and the array grows only when asked for.
  iw_gas = -1;
  if (impl_gas_newton || impl_eos_cache) {
    iw_gas = niw;
    niw += M1_NIW_GAS;
  }
  iw_vimp = -1;
  if (impl_vimp) {
    iw_vimp = niw;
    niw += M1_NIW_VIMP;
    // time_scheme = hesdirk2 (read in Time2Init, later): the DA components
    if (t2_beonly || (pin->DoesParameterExist("rad_m1","time_scheme") &&
        pin->GetString("rad_m1","time_scheme").compare("hesdirk2") == 0)) {
      niw += M1_NIW_VIMP_T2;
    }
  }
  iw_muscl = -1;
  if (impl_muscl) {
    iw_muscl = niw;
    niw += M1_NIW_MUSCL;
  }
  if (impl_hr_sfgd != 0) {vgd_thx = true;}   // J in vgd_hr
  iw_hsf = -1;
  if (impl_hr_sf) {
    iw_hsf = niw;
    niw += M1_NIW_HSF;
  }
  // <rad_m1>/implicit_face_opac = arith (default: the cell mean of
  // M1FaceAvgX1, bitwise) | kirchhoff: every interior x1 face takes the Kirchhoff mean
  //   1/(rho kappa_T)_f = int_{T_L}^{T_R} T^3/(rho_f kappa_T(rho_f,T)) dT / int T^3 dT,
  // rho_f = sqrt(rho_L rho_R) (implicit_face_opac_ng Gauss points, default 8), so the
  // diffusive face flux ~ (Phi(T_R) - Phi(T_L))/dr, Phi = int T^3/(rho kappa) dT, is
  // monotone in each cell's own T across an unresolved opacity cliff (table opacity only)
  {const std::string sfo = pin->GetOrAddString("rad_m1","implicit_face_opac","arith");
  iw_kir = -1;
  if (sfo.compare("kirchhoff") == 0) {
    if (opacity_type != M1_OPAC_TABLE) {
      ImplFatal("<rad_m1>/implicit_face_opac = kirchhoff needs opacity = table (or "
                "harmonic)");
    }
    iw_kir = niw;
    niw += 1;
    impl_kir_ng = pin->GetOrAddInteger("rad_m1","implicit_face_opac_ng",8);
    if (impl_kir_ng < 2 || impl_kir_ng > 8) {
      ImplFatal("<rad_m1>/implicit_face_opac_ng must be in [2, 8]");
    }
  } else if (sfo.compare("harmonic") == 0) {
    // distance-weighted harmonic mean of the cell rho kappa_T (= the
    // distance-weighted arithmetic mean of the diffusivity 1/(rho kappa_T)); no lookups
    iw_kir = niw;
    niw += 1;
    impl_kir_harm = true;
  } else if (sfo.compare("arith") != 0) {
    ImplFatal("<rad_m1>/implicit_face_opac = '" + sfo
              + "' is not a choice (arith | kirchhoff | harmonic)");
  }
  }
  iw_cs2 = -1;
  // the conservative gas sources run on the cubed sphere (were refused
  // before): the write-back's conservative loop re-forms the cell F inside
  // the loop (cs block), A^old / dA are face-normal (ImplicitCs2DeltaA cs correction),
  // the theta split reads m1bcs (a seam is an interior face) and the transverse dA
  // pair crosses a seam on the destination's face normals (ImplicitCsPairHalo).  Only
  // the opt-in implicit_cs2_shift (stage-start velocity shift) is not ported.
  // <rad_m1>/implicit_am_cs (read only when named; default true on the cubed
  // sphere when <hydro>/am_transport is on): the gas's zonal momentum takes the torque
  // that makes the radiation-gas exchange conserve L_z = sum V r x (m + F/(c chat)) to
  // round-off (ImplicitSolve, after the write-back)
  if (cs_geom) {
    // (<hydro> or <mhd>: the fluid that owns the angular-momentum module)
    auto named = [&](const char *b, const char *k) {
      return pin->DoesParameterExist(b, k) && pin->GetBoolean(b, k);
    };
    impl_am_cs = named("hydro","am_transport") || named("hydro","am_vector") ||
                 named("mhd","am_transport") || named("mhd","am_vector");
    if (pin->DoesParameterExist("rad_m1","implicit_am_cs")) {
      impl_am_cs = pin->GetBoolean("rad_m1","implicit_am_cs");
    } else {
      recipe::Log(std::string("<rad_m1>/implicit_am_cs = ") +
                  (impl_am_cs ? "true" : "false") +
                  " (cubed sphere: follows <hydro|mhd>/am_transport)");
    }
    am_vec = impl_am_cs && (named("hydro","am_vector") || named("mhd","am_vector"));
    if (pin->DoesParameterExist("rad_m1","implicit_am_cs_diag")) {
      impl_am_cs_diag = pin->GetInteger("rad_m1","implicit_am_cs_diag");
    }
  }
  if (cs_geom && impl_gas_src_cs2 && cs2_shift) {
    ImplFatal("<rad_m1>/implicit_cs2_shift = true is not ported to the cubed sphere; "
              "use the default false");
  }
  if (impl_gas_src_cs2 || impl_rinr_solve) {   // A^old (3), dA (3)
    iw_cs2 = niw;
    niw += 9;   // + the per-axis mean face-equation weight om (implicit_cs2_hll = drag)
    if (impl_gas_src_cs2 && cs2_cap > 0.0) {niw += 1;}   // s (iw_cs2 + 9)
  }
  // implicit_face_mmat (rad_m1.hpp; read only when named: the old scheme's face
  // M-matrix, which scheme = clean never reads)
  if (pin->DoesParameterExist("rad_m1","implicit_face_mmat")) {
    impl_fmm = pin->GetBoolean("rad_m1","implicit_face_mmat");
  }
  // implicit_lz_form (rad_m1.hpp, AM_CONSERVING_CS_1112.md sect. 10)
  if (pin->DoesParameterExist("rad_m1","implicit_lz_form") &&
      pin->GetBoolean("rad_m1","implicit_lz_form")) {
    impl_lz = true;
    if (!cs_geom || !trans_on || !trans_x3) {
      ImplFatal("<rad_m1>/implicit_lz_form needs the cubed sphere in 3-D with implicit "
                "transport");
    }
    if (!impl_gas_src_cons || impl_gas_src_full) {
      ImplFatal("<rad_m1>/implicit_lz_form needs implicit_gas_src = conservative or "
                "conservative2 (Z_f books the face transport)");
    }
    if (impl_hr_cmem || impl_physflux) {
      ImplFatal("<rad_m1>/implicit_lz_form is not combined (yet) with implicit_hr_cmem "
                "or implicit_physflux (their momentum face F_c is not the applied flux)");
    }
    iw_lz = niw;
    niw += 1;
    if (pin->DoesParameterExist("rad_m1","implicit_lz_mode")) {
      const std::string lm = pin->GetString("rad_m1","implicit_lz_mode");
      if (lm.compare("stress") == 0) {
        impl_lz_mode = 1;
      } else if (lm.compare("full") != 0) {
        ImplFatal("<rad_m1>/implicit_lz_mode must be full or stress");
      }
    }
    if (pin->DoesParameterExist("rad_m1","implicit_lz_diag")) {
      impl_lz_diag = pin->GetInteger("rad_m1","implicit_lz_diag");
    }
    Kokkos::realloc(lz_d, 6);
    Kokkos::deep_copy(lz_d, 0.0);
  }
  // implicit_lzx (rad_m1.hpp, LZPHYS_1113.md); it closes the solve's whole L_z budget, so
  // the implicit_am_cs deposit (a partial closure, after the write-back) is off inert
  // (message, no slot) on a Cartesian mesh and with the gas fixed: no gas + radiation L_z
  // budget to close there.  Read only when named (the old scheme's; scheme = clean never
  // reads it)
  const bool lzxq = pin->DoesParameterExist("rad_m1","implicit_lzx") &&
                    pin->GetBoolean("rad_m1","implicit_lzx");
  if (lzxq && (!(cs_geom || sph_geom) || !(coupling && gas_feedback))) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_lzx: inert (needs the cubed sphere or spherical "
                << "polar and coupling with gas_feedback)" << std::endl;
    }
  } else if (lzxq) {
    impl_lzx = true;
    if (!trans_on || !trans_x3) {
      ImplFatal("<rad_m1>/implicit_lzx needs a 3-D mesh with implicit transport");
    }
    if (!impl_gas_src_cons || impl_gas_src_full) {
      ImplFatal("<rad_m1>/implicit_lzx needs implicit_gas_src = conservative or "
                "conservative2 (Z_f books the face)");
    }
    if (impl_hr_cmem || impl_physflux || impl_lz) {
      ImplFatal("<rad_m1>/implicit_lzx is not combined with implicit_hr_cmem, "
                "implicit_physflux or implicit_lz_form");
    }
    if (impl_am_cs) {
      impl_am_cs = false;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> implicit_lzx: implicit_am_cs deposit off (lzx closes the "
                  << "solve's L_z)" << std::endl;
      }
    }
    if (sph_geom && !cs_geom && !(pin->GetString("mesh","ix3_bc") == "periodic" &&
                      pin->GetString("mesh","ox3_bc") == "periodic")) {
      ImplFatal("<rad_m1>/implicit_lzx on spherical polar needs periodic x3 (phi walls "
                "would carry a booked wall torque, not implemented)");
    }
    if (pin->DoesParameterExist("rad_m1","implicit_lzx_weight")) {
      const std::string lw = pin->GetString("rad_m1","implicit_lzx_weight");
      if (lw.compare("energy") == 0) {
        impl_lzx_w = 1;
      } else if (lw.compare("none") == 0) {   // diagnostics only (LZX lines)
        impl_lzx_w = 2;
      } else if (lw.compare("blend") == 0) {
        impl_lzx_w = 3;
      } else if (lw.compare("impulse") != 0) {
        ImplFatal("<rad_m1>/implicit_lzx_weight must be impulse, blend, energy or none");
      }
    }
    if (pin->DoesParameterExist("rad_m1","implicit_lzx_diag")) {
      impl_lzx_diag = pin->GetInteger("rad_m1","implicit_lzx_diag");
    }
    iw_lzx = niw;
    niw += 1;
  }
  // (rad_m1.hpp): halo exchange of the Zalesak R only when it can change
  if (pin->DoesParameterExist("rad_m1","zalesak_halo_skip")) {
    zl_hskip = pin->GetBoolean("rad_m1","zalesak_halo_skip");
  }
  if (pin->DoesParameterExist("rad_m1","zalesak_local_rounds")) {
    zl_local = pin->GetInteger("rad_m1","zalesak_local_rounds");
  }
  // zalesak_cap_rounds (inert without implicit_face_mmat; read only when named)
  if (pin->DoesParameterExist("rad_m1","zalesak_cap_rounds")) {
    zl_cap = pin->GetInteger("rad_m1","zalesak_cap_rounds");
  }
  if (pin->DoesParameterExist("rad_m1","zalesak_colsweep")) {
    zl_colsweep = pin->GetBoolean("rad_m1","zalesak_colsweep");
  }
  // (zalesak_active: default on for a fresh clean run, v20.4, recorded)
  if ((clean && !global_variable::restart_run) ||
      pin->DoesParameterExist("rad_m1","zalesak_active")) {
    zl_active = pin->GetOrAddBoolean("rad_m1","zalesak_active",true) && (zl_local == 0) &&
                !zl_colsweep;
  }
  if (pin->DoesParameterExist("rad_m1","zalesak_converge")) {
    zl_conv = pin->GetBoolean("rad_m1","zalesak_converge");
  }
  zl_prof = (std::getenv("M1_FCT_PROF") != nullptr);
  if (impl_fmm && pin->DoesParameterExist("rad_m1","implicit_face_mmat_rowsel")) {
    fmm_rowsel = pin->GetBoolean("rad_m1","implicit_face_mmat_rowsel");
  }
  if (impl_fmm && fmm_rowsel) {   // (rad_m1.hpp)
    Kokkos::realloc(fmm_rf, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(fmm_ns, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(fmm_ns, 0);
  }
  if (impl_fmm) {
    if (impl_pos_trans == 1 || afd_on || impl_muscl_damp > 0.0) {
      ImplFatal("<rad_m1>/implicit_face_mmat is not combined with implicit_pos_transport"
                " = upwind, athena_face_diss or implicit_hr_damp");
    }
    iw_fmm = niw;
    niw += 11;   // R, R new, b, 4 transverse dt A/V, 4 applied dC (stage 2b)
    impl_fmm_real = pin->DoesParameterExist("rad_m1","implicit_face_mmat_real") &&
                    pin->GetBoolean("rad_m1","implicit_face_mmat_real");
    if (impl_fmm_real) {
      ImplFatal("<rad_m1>/implicit_face_mmat_real (stage 2a, measured negative) is not "
                "combined with the persistent lam of stage 2b");
    }
    if (!clean) {   // (scheme = clean: ImplicitSolve's face M-matrix; never read)
      Kokkos::realloc(fmm_lam, nmb, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(fmm_g, nmb, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(fmm_lo, nmb, 2, ncells3, ncells2, ncells1+1);
      Kokkos::realloc(fmm_fh, nmb, ncells3, ncells2, ncells1+1);
    }
    Kokkos::deep_copy(fmm_fh, 0.0);
    // DEBUG/ablation implicit_face_mmat_trans = false (read only when named): x1 only
    const bool fmtr = !pin->DoesParameterExist("rad_m1","implicit_face_mmat_trans") ||
                      pin->GetBoolean("rad_m1","implicit_face_mmat_trans");
    if (trans_on && fmtr && !clean) {
      Kokkos::realloc(fmm_t, nmb, 14, ncells3+1, ncells2+1, ncells1);
      Kokkos::deep_copy(fmm_t, 0.0);
    }
    Kokkos::deep_copy(fmm_lam, 1.0);
    Kokkos::deep_copy(fmm_g, 0.0);
    Kokkos::deep_copy(fmm_lo, 0.0);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_face_mmat = true: x1 face flux lam F_H + (1-lam) "
                << "F_L, lam = Zalesak factor of the explicit face flux"
                << std::endl;
    }
  }
  if (impl_eos_cache) {
    impl_nec = M1EosCacheNComp(impl_ecnt);
    Kokkos::realloc(ecache, nmb, impl_nec, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ecache, -1.0);
  }
  if (impl_cmi) {iw_cmi = niw; niw += 1;}   // the donor factor lam
  MemArenaInit(pin, nmb, clean ? niw_clean : niw, ncells3, ncells2, ncells1);
  if (impl_thin_frz > 0.0 && !clean) {
    Kokkos::realloc(thin_frz, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thin_frz, 0.0);
  }
  if (impl_opac_newton) {
    if (!impl_opac_update || !impl_gas_newton) {
      ImplFatal("<rad_m1>/implicit_opac_newton needs implicit_opac_update = true and "
                "implicit_gas_newton = true");
    }
    Kokkos::realloc(ktd, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ktd, 0.0);
    Kokkos::realloc(opn_nskip_d, 1);
    Kokkos::deep_copy(opn_nskip_d, 0.0);
  }
  if (impl_face_ktn) {
    if (part_nblk > 1) {
      ImplFatal("<rad_m1>/implicit_face_opac_n needs ONE MeshBlock along x1 (the T^n "
                "snapshot has no x1 halo)");
    }
    if (!clean) {Kokkos::realloc(ktn, nmb, ncells3, ncells2, ncells1);}
    Kokkos::deep_copy(ktn, 0.0);
  }
  if (impl_ctrelax > 0.0) {
    if (!trans_on || !impl_clag_step) {
      ImplFatal("<rad_m1>/implicit_closure_thin_relax needs transport = implicit on a "
                "multi-D mesh and implicit_closure_lag = step");
    }
    Kokkos::realloc(ctr_mem, nmb, 4, ncells3, ncells2, ncells1);
  }
  Kokkos::deep_copy(iw, 0.0);
  if (impl_pred && !(clean && cln_not2h)) {
    Kokkos::realloc(ipred, nmb, (impl_pord == 2) ? 5 : 3, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ipred, 0.0);
    pred_ok = false;
  }
  // (scheme = clean: the face work arrays ifw, ifw2, ifw3 are ImplicitSolve's; not made)
  if (!clean) {Kokkos::realloc(ifw, nmb, M1_NIFW, ncells3, ncells2, ncells1+1);}
  Kokkos::deep_copy(ifw, 0.0);
  if (impl_pos_trans != 0 && sph_geom && pt_gam > 0.0) {
    Kokkos::realloc(ptq1, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::deep_copy(ptq1, 0.0);
  }
  if (trans_on) {
    Kokkos::realloc(f0x2, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2, 0.0);
    if (!clean) {Kokkos::realloc(f0x2n, nmb, ncells3, ncells2+1, ncells1);}
    Kokkos::deep_copy(f0x2n, 0.0);
    if (impl_hr_cmem) {
      Kokkos::realloc(fax2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(fax2, 0.0);
    }
    if (trans_x3) {
      Kokkos::realloc(f0x3, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3, 0.0);
      if (!clean) {Kokkos::realloc(f0x3n, nmb, ncells3+1, ncells2, ncells1);}
      Kokkos::deep_copy(f0x3n, 0.0);
      if (impl_hr_cmem) {
        Kokkos::realloc(fax3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(fax3, 0.0);
      }
      if (impl_gas_src_cons && !clean) {
        Kokkos::realloc(gcz3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(gcz3, 0.0);
      }
    }
    if (impl_gas_src_cons && !clean) {
      Kokkos::realloc(gcz2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(gcz2, 0.0);
    }
    if (impl_tlim != M1_TLIM_NONE || blat_on) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
      }
    }
    // the outer Marshak q of a half-range solve with ray data = h+_1 of the top cell (the
    // vacuum sends nothing back), per column, in vcol_q (read only when named; default
    // cell = bitwise): implicit_hr_q = col keeps the outer Marshak q of the closure's own
    // surface estimate (vet_col_surface_q / vet_gd_replace: H at the FACE over J of the
    // top cell) instead of the top cell's CENTRE ratio h+_1 = H+_N/J_N.  At tau_cell ~ 1
    // the centre value lies an optical depth tau_cell/2 below the vacuum face, where the
    // outgoing intensity is far less one-sided (gda s1: h+ = 0.38 at tau 0.63), so c q
    // E_N under-estimates the escaping flux by O(tau_cell): the plateau of ord1020 (PROD
    // retains +1.2 % energy at s1/s2). Only with closure vet_col (+ vet_gd), whose
    // surface q exists; vet_sc keeps h+_1
    bool hrq_col = false;
    if (pin->DoesParameterExist("rad_m1","implicit_hr_q")) {
      const std::string sq = pin->GetString("rad_m1","implicit_hr_q");
      if (sq.compare("col") == 0) {
        hrq_col = true;
      } else if (sq.compare("cell") != 0) {
        ImplFatal("<rad_m1>/implicit_hr_q = '" + sq + "' is not a choice (cell | col)");
      }
    }
    hr_q = impl_beam_hr && (vet_sc || (vgd_on && sph_geom)) &&
           (ibc_x1max == M1_IBC_MARSHAK) && !(hrq_col && vet_col);
    if (hr_q && vcol_q.extent_int(0) < nmb) {
      Kokkos::realloc(vcol_q, nmb, ncells3, ncells2);
      Kokkos::deep_copy(vcol_q, marshak_q);
    }
    if (blat_on && !clean) {
      // implicit_thin_switch also keeps the lagged face flux DG there
      const int nfw = (impl_thsw != M1_THSW_NONE) ? M1_NIFW : 3;
      Kokkos::realloc(ifw2, nmb, nfw, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(ifw2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(ifw3, nmb, nfw, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(ifw3, 0.0);
      }
    }
    if (impl_tlim == M1_TLIM_LPR && !blat_on && !clean) {
      // the sp rows (M1SphTransRow, the face/cell kernels' blend terms) read
      // the limited face theta thx2/thx3 through their blend flag; with no blend its
      // weights stay identically zero (AL = HCL = HCR = 0), so only theta changes
      Kokkos::realloc(ifw2, nmb, 3, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(ifw2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(ifw3, nmb, 3, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(ifw3, 0.0);
      }
    }
    if (impl_tlim != M1_TLIM_NONE) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      Kokkos::realloc(klx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(klx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
        Kokkos::realloc(klx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(klx3, 0.0);
      }
    }
    // the transverse halo rides the module's ORDINARY cell-centred boundary machinery on
    // a scratch array, which is what gives it periodic wrap, corner/edge neighbours and
    // MPI for free; the hand-rolled x1 halo of 3b is then not used at all under
    // transport = implicit (this exchange carries its six quantities and seven more).
    // (scheme = clean: made at its first exchange, ImplicitHaloExchange; the clean solve
    // has its own halo, CleanGhostFill)
    if (!clean) {Kokkos::realloc(thw, nmb, M1_NHALO_T, ncells3, ncells2, ncells1);}
    Kokkos::realloc(thw_c, nmb, M1_NHALO_T, 1, 1, 1);
    pbval_th = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_th->exch_depth = M1HaloDepthEnv(pmy_pack, pin);
    impl_hdepth = pbval_th->exch_depth;   // G5: also the halo_mpi cell lists
    if (!clean) {pbval_th->InitializeBuffers(M1_NHALO_T);}   // (clean: with thw)
    // COMPONENT ROLES of the scratch halos (SetVectorPairs).  Slots 2,3 hold F1 and KT,
    // SCALARS under the polar flip and the seam transform; the real tangential pairs
    // are (N2,N3), (A2,A3), (V2,V3).  The default (IVY, IVZ) rule would get all of that
    // wrong (a latent trap while rad_m1 is a fatal on sp/cs, see rad_m1.cpp).
    {
      auto slot = [](const int iw) {
        int s = -1;
        for (int n=0; n<M1_NHALO_T; ++n) {
          if (M1HaloCompT(n) == iw) {s = n;}
        }
        return s;
      };
      pbval_th->SetVectorPairs(M1_NHALO_T, {{slot(M1_IW_N2), slot(M1_IW_N3)},
                                            {slot(M1_IW_A2), slot(M1_IW_A3)},
                                            {slot(M1_IW_V2), slot(M1_IW_V3)}});
    }
    // the NARROW exchange of the same list: the M1_NHALO_Q components a Picard pass can
    // move once the closure is frozen.  A separate array because the exchange takes the
    // variable count from the array's second extent, and a prefix subview of a
    // LayoutRight array is not one.
    // (scheme = clean: thq too is made at its first exchange, as thw)
    if (!clean) {Kokkos::realloc(thq, nmb, M1_NHALO_Q, ncells3, ncells2, ncells1);}
    Kokkos::deep_copy(thq, 0.0);
    Kokkos::realloc(thq_c, nmb, M1_NHALO_Q, 1, 1, 1);
    pbval_tq = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_tq->exch_depth = M1HaloDepthEnv(pmy_pack, pin);
    if (!clean) {pbval_tq->InitializeBuffers(M1_NHALO_Q);}   // (clean: with thq)
    pbval_tq->SetVectorPairs(M1_NHALO_Q, {});   // E', G0, F1, KT: all scalars
    // the deep interior of the scratch arrays is neither read by a send nor written by a
    // receive when every neighbour is at the SAME level (a same-level buffer reaches ng
    // cells in from the active boundary, buffs_cc.cpp) and neither the cubed-sphere
    // resample nor the polar transform is in play; then the copies to and from iw can
    // skip it.  SMR/AMR is already a fatal above.
    halo_shell = !(pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere ||
                   pmy_pack->pmesh->use_polar_boundary);
    // (test switch, env M1_CS_HALO_SHELL=1): the shell copy on the cubed
    // sphere too (the seam resample and the corner fill read the outermost ng active
    // cells only, which the shell holds; gate bitwise before use)
    {
      const char *e = std::getenv("M1_CS_HALO_SHELL");
      if (e != nullptr && std::atoi(e) != 0 && pmy_pack->pmesh->use_cubed_sphere &&
          !pmy_pack->pmesh->multilevel && !pmy_pack->pmesh->use_polar_boundary) {
        halo_shell = true;
      }
    }
    {
      // ONE more exchange object, for the single Krylov vector the operator application
      // needs in its ghost zones -- and for the one-component exchange of E alone, which
      // is why it is allocated under every solver.  It is used strictly SEQUENTIALLY
      // with pbval_th and pbval_tq (each exchange runs its
      // InitRecv/Send/Recv/Clear chain to completion before the next
      // starts) and every rank issues the identical SEQUENCE of exchanges -- the Picard
      // count, the BiCGStab count and every breakdown decision are taken from GLOBAL
      // reductions -- so MPI's non-overtaking guarantee keeps the two streams apart even
      // though they share the tag space.
      Kokkos::realloc(krw, nmb, 1, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(krw, 0.0);
      Kokkos::realloc(krw_c, nmb, 1, 1, 1, 1);
      pbval_kr = new MeshBoundaryValuesCC(pmy_pack, pin, false);
      pbval_kr->exch_depth = M1HaloDepthEnv(pmy_pack, pin);
    pbval_kr->InitializeBuffers(1);
      pbval_kr->SetVectorPairs(1, {});   // one scalar Krylov vector
      // implicit_vimp: the per-pass exchange of the Jacobian rows and dv^k (same
      // sequential-use argument as above).  The scratch arrays of the solve from here on
      // are not allocated with scheme = clean, which never runs it.
      if (impl_vimp && !clean) {
        Kokkos::realloc(vmw, nmb, M1_NVIMP_X, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(vmw, 0.0);
        Kokkos::realloc(vmw_c, nmb, M1_NVIMP_X, 1, 1, 1);
        pbval_vm = new MeshBoundaryValuesCC(pmy_pack, pin, false);
        pbval_vm->exch_depth = M1HaloDepthEnv(pmy_pack, pin);
    pbval_vm->InitializeBuffers(M1_NVIMP_X);
        // (DV2, DV3) is the tangential pair; the nine per-direction Jacobian rows P are
        // treated as scalars (their proper curvilinear treatment is a later stage)
        pbval_vm->SetVectorPairs(M1_NVIMP_X, {{M1_IV_DV + 1, M1_IV_DV + 2}});
      }
      // the (scalar, face-normal pair) exchange of the conservative
      // sources' dA on the cubed sphere (ImplicitCsPairHalo)
      if (cs_geom && !clean) {
        Kokkos::realloc(d3w, nmb, 3, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(d3w, 0.0);
        Kokkos::realloc(d3w_c, nmb, 3, 1, 1, 1);
        pbval_d3 = new MeshBoundaryValuesCC(pmy_pack, pin, false);
        pbval_d3->exch_depth = M1HaloDepthEnv(pmy_pack, pin);
        pbval_d3->InitializeBuffers(3);
        pbval_d3->SetVectorPairs(3, {{1, 2}});
        pbval_d3->cs_noresample = true;
        pbval_d3->cs_perm_pairs = true;
      }
    }
    if (impl_halo_direct) {ImplicitHaloDirectInit();}
    if (impl_odc && !clean) {
      Kokkos::realloc(odc, nmb, 3, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(odc, 0.0);
    }
    if (impl_stencil && !clean) {
      if (ibc_x1min == M1_IBC_PERIODIC) {
        ImplFatal("<rad_m1>/implicit_op_stencil does not take a periodic x1 wrap");
      }
      Kokkos::realloc(ost, nmb, impl_vfold ? 25 : 19, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(ost, 0.0);
    }
  }
  // MILESTONE 3e: the Anderson histories.  Allocated ONLY when the acceleration is on.
  if (impl_accel != M1_IACC_NONE && !clean) {
    aa_nc = trans_on ? 4 : 2;
    Kokkos::realloc(aa_xc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xc, 0.0);
    Kokkos::realloc(aa_fc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fc, 0.0);
    Kokkos::realloc(aa_xp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xp, 0.0);
    Kokkos::realloc(aa_fp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fp, 0.0);
    Kokkos::realloc(aa_dx, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_dx, 0.0);
    Kokkos::realloc(aa_df, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_df, 0.0);
    Kokkos::realloc(aa_sc, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_sc, 1.0);
  }
  ImplicitPartitionInit();

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: transport=" << (full ? "implicit" : "implicit_x1")
              << (full ? (trans_x3 ? " (3-D)" : (trans_on ? " (2-D)" : " (1-D)")) : "")
              << " solver=" << (full ? (bicg_on ? "bicgstab" : "line_jacobi") : "thomas")
              << " lin_tol=" << impl_lin_tol
              << " lin_maxit=" << impl_lin_maxit << std::endl;
    std::cout << "         implicit_cfl="
              << impl_cfl << " tol=" << impl_tol << " maxit=" << impl_maxit
              << " opac_update=" << (impl_opac_update ? "true" : "false")
              << " marshak_q=" << marshak_q << std::endl;
    std::cout << "         implicit_flux="
              << ((impl_flux == M1_IFLUX_APHLL) ? "ap_hll" :
                  ((impl_flux == M1_IFLUX_BERTHON) ? "berthon" :
                   ((impl_flux == M1_IFLUX_BLEND) ? "blend" : "central")))
              << ((impl_flux == M1_IFLUX_BLEND) ?
                  (":" + pin->GetString("rad_m1","implicit_blend") + "/" +
                   pin->GetString("rad_m1","implicit_blend_fmode") + "/" +
                   pin->GetString("rad_m1","implicit_blend_mode")) : std::string(""))
              << " implicit_recon="
              << ((impl_recon == M1_IRECON_PLMDC) ? "plm_dc" : "dc")
              << " implicit_partition="
              << ((impl_part == M1_IPART_GATHER) ? "gather" : "none") << std::endl;
    if (impl_vimp) {
      std::cout << "         implicit_vimp=true (v' of the enthalpy flux implicit)"
                << std::endl;
    }
    if (impl_enth != M1_IENTH_UPWIND) {
      std::cout << "         implicit_enthalpy="
                << ((impl_enth == M1_IENTH_PLM) ? "plm"
                    : (impl_enth == M1_IENTH_PLMC) ? "plmc"
                    : (impl_enth == M1_IENTH_LIN) ? "lin" : "central")
                << " (deferred correction)" << std::endl;
    }
    if (ibc_x1min == M1_IBC_MARSHAK || ibc_x1max == M1_IBC_MARSHAK) {
      std::cout << "         implicit_bc_advect=" << (impl_bc_advect ? "true" : "false")
                << " (Marshak x1 end carries A E)" << std::endl;
    }
    if (trans_on) {
      std::cout << "         implicit_offdiag="
                << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                    ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
                << (od_auto ? " (auto)" : "")
                << " closure_relax=" << impl_crelax
                << (impl_crelax_thin ? " (thin faces only)" : "")
                << " closure_lag=" << (impl_clag_step ? "step" : "pass")
                << " trans_limit="
                << ((impl_tlim == M1_TLIM_LP) ? "lp"
                    : ((impl_tlim == M1_TLIM_LPR) ? "lpr" : "none"))
                << ((impl_tlim == M1_TLIM_LP)
                    ? (" fmax=" + std::to_string(impl_tfmax)) : std::string(""))
                << std::endl;
    }
    std::cout << "         implicit_accel="
              << ((impl_accel == M1_IACC_ANDERSON) ? "anderson" : "none")
              << ((impl_accel == M1_IACC_ANDERSON)
                  ? (" m=" + std::to_string(impl_and_m)
                     + " beta=" + std::to_string(impl_and_beta)
                     + " start=" + std::to_string(impl_and_start)
                     + " ncomp=" + std::to_string(aa_nc))
                  : std::string("")) << std::endl;
    if (impl_gas_newton || impl_eos_cache) {
      std::cout << "         implicit_gas_newton="
                << (impl_gas_newton ? "true" : "false")
                << " implicit_eos_cache=" << (impl_eos_cache ? "true" : "false")
                << (impl_eos_cache ? (" nt=" + std::to_string(impl_ecnt)
                                      + " ncomp=" + std::to_string(impl_nec)
                                      + " check="
                                      + (impl_eccheck ? "true" : "false"))
                                   : std::string("")) << std::endl;
    }
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {VetInit(pin);}
  ImplBeamInit(pin);   // (rad_m1_vet.cpp)
  Time2Init(pin);
  MRInit(pin);         // implicit_mr_every (rad_m1_mr.cpp)
  // (memory): u1 is used only by the explicit stage copy, hesdirk2 and the
  // multi-rate step; a backward-Euler implicit run without multi-rate frees it
  if (transport >= M1_TRANSPORT_IMPLICIT_X1 && time_scheme == M1_TIME_BE &&
      mr_every <= 1) {
    Kokkos::realloc(u1, u1.extent(0), M1_NVAR, 1, 1, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPartitionInit
//! \brief milestone 3b, LIMIT 4: build the topology of the x1 stacks and allocate the
//! gather/scatter and halo buffers.  A "stack" is the set of MeshBlocks that share the
//! same (x2,x3) footprint, ordered by their x1 logical location; its ROOT is the block
//! with the lowest one.  Uniform mesh only (the caller has already fatalled on SMR/AMR
//! and on a non-integer block count along x1), so the stack of a block is found by a
//! scan of the global LogicalLocation list.
//!
//! Nothing is allocated and nothing is communicated when part_nblk == 1: the 3a path is
//! then bitwise what it was.

void RadiationM1::ImplicitPartitionInit() {
  part_nroot = 0;
  part_nx1g = 0;
  part_nlay = 0;
  part_nqa = M1_NHALO_A;
  part_any_mpi = false;

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int g0 = pmy_pack->gids;
  // the three per-block tables are ALWAYS allocated: ImplicitSolve reads part_pos in
  // every kernel, and with one MeshBlock per column the defaults (position 0 of a stack
  // of one) select exactly the 3a branches.
  Kokkos::realloc(part_pos, nmb);
  Kokkos::realloc(part_slot, nmb);
  Kokkos::realloc(part_nbr, 2*nmb);
  for (int m=0; m<nmb; ++m) {
    part_pos.h_view(m) = 0;
    part_slot.h_view(m) = -1;
    part_nbr.h_view(2*m) = -1;
    part_nbr.h_view(2*m+1) = -1;
  }
  if (part_nblk <= 1) {
    part_pos.modify_host();
    part_pos.sync_device();
    part_slot.modify_host();
    part_slot.sync_device();
    part_nbr.modify_host();
    part_nbr.sync_device();
    return;
  }
  part_nx1g = part_nblk*indcs.nx1;
  part_nlay = std::min(indcs.ng, 2);

  // (1) for every LOCAL block: its position in the stack, the global ids of the stack
  // members and of its two x1 neighbours.
  part_rootgid.assign(nmb, -1);
  part_rootrank.assign(nmb, -1);
  part_nbrrank.assign(2*nmb, -1);
  part_nbrgid.assign(2*nmb, -1);
  std::vector<int> stack(part_nblk);
  for (int m=0; m<nmb; ++m) {
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {stack[p] = -1;}
    for (int g=0; g<pm->nmb_total; ++g) {
      LogicalLocation &lg = pm->lloc_eachmb[g];
      if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level) {
        if (lg.lx1 >= 0 && lg.lx1 < part_nblk) {stack[lg.lx1] = g;}
      }
    }
    for (int p=0; p<part_nblk; ++p) {
      if (stack[p] < 0) {
        ImplFatal("<rad_m1>/implicit_partition = gather: the x1 stack of a MeshBlock is "
                  "incomplete (a non-uniform mesh?)");
      }
    }
    int pos = static_cast<int>(lm.lx1);
    part_pos.h_view(m) = pos;
    part_rootgid[m] = stack[0];
    part_rootrank[m] = pm->rank_eachmb[stack[0]];
    if (pos > 0) {
      part_nbrgid[2*m] = stack[pos-1];
      part_nbrrank[2*m] = pm->rank_eachmb[stack[pos-1]];
    }
    if (pos < part_nblk-1) {
      part_nbrgid[2*m+1] = stack[pos+1];
      part_nbrrank[2*m+1] = pm->rank_eachmb[stack[pos+1]];
    }
    for (int s=0; s<2; ++s) {
      int gn = part_nbrgid[2*m+s];
      bool loc = (gn >= 0) && (part_nbrrank[2*m+s] == global_variable::my_rank);
      part_nbr.h_view(2*m+s) = loc ? (gn - g0) : -1;
      if (gn >= 0 && part_nbrrank[2*m+s] != global_variable::my_rank) {
        part_any_mpi = true;
      }
    }
    if (part_rootrank[m] != global_variable::my_rank) {part_any_mpi = true;}
  }

  // (2) the local ROOTS and their member lists
  part_mrank.clear();
  part_mgid.clear();
  std::vector<int> rootmb;
  for (int m=0; m<nmb; ++m) {
    if (part_pos.h_view(m) != 0) continue;
    rootmb.push_back(m);
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {
      int gfound = -1;
      for (int g=0; g<pm->nmb_total; ++g) {
        LogicalLocation &lg = pm->lloc_eachmb[g];
        if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level && lg.lx1 == p) {
          gfound = g;
        }
      }
      part_mgid.push_back(gfound);
      part_mrank.push_back(pm->rank_eachmb[gfound]);
      if (pm->rank_eachmb[gfound] != global_variable::my_rank) {part_any_mpi = true;}
    }
  }
  part_nroot = static_cast<int>(rootmb.size());
  for (int m=0; m<nmb; ++m) {
    int sl = -1;
    if (part_rootrank[m] == global_variable::my_rank) {
      for (int s=0; s<part_nroot; ++s) {
        if (rootmb[s] + g0 == part_rootgid[m]) {sl = s;}
      }
    }
    part_slot.h_view(m) = sl;
  }
  part_pos.modify_host();
  part_pos.sync_device();
  part_slot.modify_host();
  part_slot.sync_device();
  part_nbr.modify_host();
  part_nbr.sync_device();

  // (3) buffers
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(part_sys, std::max(part_nroot,1), 6, ncells3, ncells2, part_nx1g);
  Kokkos::deep_copy(part_sys, 0.0);
  int nrow = 4*ncells3*ncells2*indcs.nx1;
  Kokkos::realloc(part_sbuf, nmb, nrow);
  Kokkos::realloc(part_rbuf, std::max(part_nroot*part_nblk,1), nrow);
  Kokkos::realloc(part_sbuf_h, nmb, nrow);
  Kokkos::realloc(part_rbuf_h, std::max(part_nroot*part_nblk,1), nrow);
  int nhal = M1_NHALO_A*part_nlay*ncells3*ncells2;
  Kokkos::realloc(part_hbuf, nmb, 4, nhal);
  Kokkos::realloc(part_hbuf_h, nmb, 4, nhal);

  if (global_variable::my_rank == 0) {
    std::cout << "         implicit_partition=gather: " << part_nblk
              << " MeshBlocks per x1 column, " << part_nx1g << " rows per gathered line, "
              << part_nlay << " halo layers" << std::endl;
  }
}

// the depth of the implicit scratch exchanges on the cubed sphere,
// <rad_m1>/implicit_halo_depth (default 2; 0 or >= ng = all ng layers, the old
// behaviour), env M1_HALO_DEPTH overrides.  The implicit kernels read 2 ghost layers
// (the vimp enthalpy stencil +-2; poisoned-ghost gate M1_POISON_GHOST: NaN beyond depth
// 2 leaves the outputs identical, beyond depth 1 does not), so depth 2 sends the same
// values that are read with 2/3 of the bytes (RSG 2-node loop -2.5 %, bitwise).
// G5 (USER 11-04): also on the spherical-polar mesh (wedge and full sphere:
// the pole partner's rows are reversed within the strip, so a depth-2 strip fills the
// 2 ghost rows next to the pole; poisoned-ghost gate incl. the polar rows) when x1 is
// one MeshBlock (the gate poisons x2/x3 ghosts only; every sp production is x1-unsplit);
// not on Cartesian meshes (the old behaviour, all ng layers).
static int M1HaloDepthEnv(MeshBlockPack *pp, ParameterInput *pin) {
  Mesh *pm = pp->pmesh;
  const bool sp = pm->use_spherical_polar && (pm->mesh_indcs.nx1 == pm->mb_indcs.nx1);
  if (!pm->use_cubed_sphere && !sp) {return 0;}
  int d = 2;
  if (pin != nullptr && pin->DoesParameterExist("rad_m1", "implicit_halo_depth")) {
    d = pin->GetInteger("rad_m1", "implicit_halo_depth");
  }
  const char *e = std::getenv("M1_HALO_DEPTH");
  if (e != nullptr) {d = std::atoi(e);}
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloCopy
//! \brief copy `nq` components between iw and the scratch halo array `sc` (`topack`
//! picks the direction), over the HALO SHELL only.
//!
//! What the exchange reads out of `sc` is the outermost ng ACTIVE cells, and what it
//! writes back into it is the ghost zones (buffs_cc.cpp, isame); with same-level
//! neighbours only -- which is all the implicit solve allows -- nothing else in `sc` is
//! ever touched.  The box [is+ng,ie-ng] x [js+ng,je-ng] x [ks+ng,ke-ng] is therefore
//! neither sent nor received and need not be copied in either direction: on the way in
//! its value is never read, and on the way back it is a copy of what the pack put there.
//! At a 2-D 84 x 32 block with ng = 2 that is 480 cells per component instead of 3168.
//! `halo_shell` falls back to the full copy for the exchanges that reach deeper (the
//! cubed-sphere resample and the polar transform read the whole strip).
//!
//! The loop carries (m,n,k,j) and runs the contiguous i direction inside, so the index
//! arithmetic of the flattened range policy (an integer division per index) is paid once
//! per ROW instead of once per cell.

void RadiationM1::ImplicitHaloCopy(DvceArray5D<Real> &sc, int nq, int c0, bool topack) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // the deep interior, as an index box.  A direction with no neighbour (a degenerate
  // dimension) is interior everywhere; a direction too thin to have one (nx <= 2 ng), or
  // a run whose exchange reaches deeper, makes the box empty and the copy full.
  int ng = indcs.ng;
  int ilo = halo_shell ? (indcs.is + ng) : n1;
  int ihi = halo_shell ? (indcs.ie - ng) : (n1 - 1);
  int jlo = (indcs.nx2 > 1) ? (indcs.js + ng) : 0;
  int jhi = (indcs.nx2 > 1) ? (indcs.je - ng) : 0;
  int klo = (indcs.nx3 > 1) ? (indcs.ks + ng) : 0;
  int khi = (indcs.nx3 > 1) ? (indcs.ke - ng) : 0;
  if (ilo > ihi) {
    ilo = n1;
    ihi = n1 - 1;
  }
  const int il_ = ilo, iu_ = ihi, jl_ = jlo, ju_ = jhi, kl_ = klo, ku_ = khi;
  const int nc0 = c0, nn1 = n1;
  const bool pack_ = topack;
  auto iw_ = iw;
  auto sc_ = sc;
  // test switches (env): M1_HCOPY_CELL=1 the one-thread-per-cell copy for
  // every exchange (the same copies); M1_POISON_GHOST=1 (gate only, NOT a run mode): on
  // the way back, every ghost cell more than ONE cell outside the active range in x2 or
  // x3 is set to NaN (if the outputs stay identical, no implicit kernel reads those
  // layers of the exchanged components)
  static int s_hcell = -1, s_poison = -1;
  if (s_hcell < 0) {
    const char *e = std::getenv("M1_HCOPY_CELL");
    s_hcell = (e != nullptr) ? std::atoi(e) : 0;
    const char *f = std::getenv("M1_POISON_GHOST");
    s_poison = (f != nullptr) ? std::atoi(f) : 0;
  }
  const bool poi = (s_poison != 0) && !topack;
  const int pd = s_poison;   // M1_POISON_GHOST=D: the ghost layers beyond depth D
  const int pjs = indcs.js, pje = indcs.je, pks = indcs.ks, pke = indcs.ke;
  const bool p2 = (indcs.nx2 > 1), p3 = (indcs.nx3 > 1);
  if (impl_halo_direct || s_hcell != 0 || poi) {
    // implicit_halo_direct on a mesh whose neighbours are not all on this rank: the
    // same shell copy with ONE thread per cell (coalesced along i) instead of one per
    // row -- the row loop runs ~13k threads with strided rows on the 3-D box.
#if defined(KOKKOS_ENABLE_HIP)
    // (HIP): one 64-lane vector per (m,n,k,j) row, 4 rows per team, i.e.
    // the same cells, each copied by the same assignment, with the index division paid
    // per row instead of per cell (rocprofv3, HG small cs deck: m1_impl_hcpyf 0.40 s,
    // 7466 launches, ~2.6 TB/s).  Bitwise (plain copies).  Env M1_HCOPY_TEAM=0 restores
    // the flat launch.
    static int s_hteam = -1;
    if (s_hteam < 0) {
      const char *e = std::getenv("M1_HCOPY_TEAM");
      s_hteam = (e != nullptr) ? std::atoi(e) : 1;
    }
    if (s_hteam != 0 && !poi) {
      const int nr2 = n2, nr3 = n3, nqq = nq;
      const int nrow = (nmb1 + 1)*nq*n3*n2;
      constexpr int kTr = 4;
      Kokkos::TeamPolicy<DevExeSpace, Kokkos::LaunchBounds<256, 2>>
          pol(DevExeSpace(), (nrow + kTr - 1)/kTr, kTr, 64);
      Kokkos::parallel_for("m1_impl_hcpyt", pol, KOKKOS_LAMBDA(TeamMember_t tm) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, kTr), [&](const int r) {
          const int row = tm.league_rank()*kTr + r;
          if (row >= nrow) return;
          const int m = row/(nqq*nr3*nr2);
          const int n = (row - m*(nqq*nr3*nr2))/(nr3*nr2);
          const int k = (row - (m*nqq + n)*(nr3*nr2))/nr2;
          const int j = row - ((m*nqq + n)*nr3 + k)*nr2;
          const bool inr = (k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_);
          const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tm, nn1), [&](const int i) {
            if (inr && (i >= il_) && (i <= iu_)) return;
            if (pack_) {
              sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
            } else {
              iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
            }
          });
        });
      });
      return;
    }
#endif
    par_for("m1_impl_hcpyf", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_) && (i >= il_) &&
          (i <= iu_)) {
        return;
      }
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      if (pack_) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      } else if (poi && ((p2 && (j < pjs - pd || j > pje + pd)) ||
                         (p3 && (k < pks - pd || k > pke + pd)))) {
        iw_(m,nc,k,j,i) = std::numeric_limits<Real>::quiet_NaN();
      } else {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    });
    return;
  }
  par_for("m1_impl_hcpy", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j) {
    const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
    // the two i runs this row copies: the whole row unless the row is interior
    int ia = 0, ib = nn1 - 1, ic = nn1, id = nn1 - 1;
    if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_)) {
      ib = il_ - 1;
      ic = iu_ + 1;
    }
    if (pack_) {
      for (int i=ia; i<=ib; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
    } else {
      for (int i=ia; i<=ib; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloExchange
//! \brief the common core: pack, exchange through the ordinary cell-centred boundary
//! machinery on the PERSISTENT scratch array whose width matches `nq`, unpack.  The
//! three widths (M1_NHALO_T, M1_NHALO_Q, 1) have one array and one boundary object each;
//! they are used strictly sequentially, so they may share the MPI tag space.

void RadiationM1::ImplicitHaloExchange(int nq, int c0) {
  if (halo_direct_on) {   // implicit_halo_direct, every neighbour on this rank
    ImplicitHaloDirect(nq, c0);
    return;
  }
  if (impl_halo_mpi) {     // implicit_halo_mpi: on-rank copy + one message per rank
    if (hm_state == 0) {ImplicitHaloMPIInit();}
    if (hm_state == 1) {
      ImplicitHaloMPI(nq, c0);
      return;
    }
  }
  DvceArray5D<Real> *pa, *pc;
  MeshBoundaryValuesCC *pb;
  if (nq == M1_NHALO_T) {
    pa = &thw; pc = &thw_c; pb = pbval_th;
  } else if (nq == M1_NHALO_Q) {
    pa = &thq; pc = &thq_c; pb = pbval_tq;
  } else if (nq == M1_NVIMP_X) {
    pa = &vmw; pc = &vmw_c; pb = pbval_vm;
  } else {
    pa = &krw; pc = &krw_c; pb = pbval_kr;
  }
  if (pa->extent_int(0) == 0) {   // (thw under scheme = clean: made at its first use)
    auto &id = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(*pa, pmy_pack->nmb_thispack, nq, (id.nx3 > 1) ? id.nx3 + 2*id.ng : 1,
                    (id.nx2 > 1) ? id.nx2 + 2*id.ng : 1, id.nx1 + 2*id.ng);
    pb->InitializeBuffers(nq);
  }
  ImplicitHaloCopy(*pa, nq, c0, true);
  while (pb->InitRecv(nq) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
  ImplicitHaloCopy(*pa, nq, c0, false);
  // STAGE CS2 (C3): (N2, N3) of the work array are FACE-NORMAL components (as A and V)
  // and cross a seam re-expressed on the destination's face normals (cs_perm_pairs,
  // CubedS1Init; was the signed axis permutation); it is renormalised in the
  // ghost cells all the same (a no-op for the cs closures: n = 0 or rhat),
  // so that a later change of the seam map cannot hand the closure a non-unit n.  A zero
  // n (eddington; vet_col about r_hat has n = r_hat) is left alone.
  if (cs_geom && nq == M1_NHALO_T) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int n1 = indcs.nx1 + 2*indcs.ng, n2 = indcs.nx2 + 2*indcs.ng;
    const int n3 = indcs.nx3 + 2*indcs.ng;
    auto iw_ = iw;
    par_for("m1_cs_nrenorm", DevExeSpace(), 0, pmy_pack->nmb_thispack-1, 0, n3-1, 0, n2-1,
            0, n1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) return;
      const Real a = iw_(m,M1_IW_N1,k,j,i), b = iw_(m,M1_IW_N2,k,j,i);
      const Real d = iw_(m,M1_IW_N3,k,j,i);
      const Real nn = a*a + b*b + d*d;
      if (!(nn > 0.0)) return;
      const Real f = 1.0/sqrt(nn);
      iw_(m,M1_IW_N1,k,j,i) = a*f;
      iw_(m,M1_IW_N2,k,j,i) = b*f;
      iw_(m,M1_IW_N3,k,j,i) = d*f;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirectInit
//! \brief implicit_halo_direct: tabulate, for every MeshBlock of the pack and each of the
//! 26 directions (ox1,ox2,ox3), the LOCAL index of the same-level neighbour that fills
//! that ghost region, or -1 where there is none (a physical boundary: its ghost zones
//! are not filled by the ordinary exchange either).  The direct copy is used only when
//! EVERY rank finds every neighbour on its own rank at the same level, outside the
//! cubed-sphere and polar transforms (so all ranks take the same branch and no MPI
//! exchange is left half-posted); otherwise the ordinary exchange runs, as before.

void RadiationM1::ImplicitHaloDirectInit() {
  auto *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack;
  hd_src = DualArray2D<int>("m1_hd_src", nmb, 27);
  int ok = (pm->multilevel || pm->use_cubed_sphere || pm->use_polar_boundary) ? 0 : 1;
  auto &nb = pmy_pack->pmb->nghbr;
  auto &lev = pmy_pack->pmb->mb_lev;
  const int e2 = pm->multi_d ? 1 : 0;
  const int e3 = pm->three_d ? 1 : 0;
  for (int m = 0; m < nmb; ++m) {
    for (int d = 0; d < 27; ++d) {hd_src.h_view(m,d) = -1;}
    for (int o3 = -e3; o3 <= e3; ++o3) {
      for (int o2 = -e2; o2 <= e2; ++o2) {
        for (int o1 = -1; o1 <= 1; ++o1) {
          if (o1 == 0 && o2 == 0 && o3 == 0) continue;
          int n = NeighborIndex(o1, o2, o3, 0, 0);
          if (n < 0 || n >= pmy_pack->pmb->nnghbr) {ok = 0; continue;}
          const NeighborBlock &q = nb.h_view(m,n);
          if (q.gid < 0) continue;
          if (q.rank != global_variable::my_rank || q.lev != lev.h_view(m)) {
            ok = 0;
            continue;
          }
          hd_src.h_view(m, (o1+1) + 3*(o2+1) + 9*(o3+1)) = q.gid - pmy_pack->gids;
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  {int g = ok;
  MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  ok = g;}
#endif
  halo_direct_on = (ok == 1);
  hd_src.modify_host();
  hd_src.sync_device();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirect
//! \brief implicit_halo_direct: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list) filled by ONE kernel that copies each ghost
//! cell from the active cell of the neighbour that owns it (all neighbours are on this
//! rank at the same level, ImplicitHaloDirectInit).  The same numbers the pack /
//! exchange / unpack chain delivers: 1 launch where there were 4, and no host work.

void RadiationM1::ImplicitHaloDirect(int nq, int c0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int n1 = nx1 + 2*indcs.ng;
  const int n2 = (nx2 > 1) ? (nx2 + 2*indcs.ng) : 1;
  const int n3 = (nx3 > 1) ? (nx3 + 2*indcs.ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool md = (nx2 > 1), td = (nx3 > 1);
  auto iw_ = iw;
  auto tab = hd_src.d_view;
  const int nc0 = c0;
  par_for("m1_impl_hdir", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    const int o1 = (i < is) ? -1 : ((i > ie) ? 1 : 0);
    const int o2 = md ? ((j < js) ? -1 : ((j > je) ? 1 : 0)) : 0;
    const int o3 = td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    if (o1 == 0 && o2 == 0 && o3 == 0) return;
    const int src = tab(m, (o1+1) + 3*(o2+1) + 9*(o3+1));
    if (src < 0) return;
    const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
    iw_(m,nc,k,j,i) = iw_(src,nc,k - o3*nx3,j - o2*nx2,i - o1*nx1);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitKrylovHalo
//! \brief milestone 3b phase C: put ONE component of the work array into the scratch
//! array `krw`, exchange it with all six neighbours through the module's ordinary
//! cell-centred boundary machinery (which is what supplies periodic wrap, edge/corner
//! neighbours and MPI), and copy the ghost zones back.  PHYSICAL boundaries are not
//! filled: the row's coefficient towards such a neighbour is identically zero, so the
//! ghost value is multiplied by zero and never read in anger.

void RadiationM1::ImplicitKrylovHalo(int comp) {
  const int kir_ = iw_kir;   // (-1 = off)
  auto cx1f_k = pmy_pack->pcoord->xx1f; (void) cx1f_k;
  ImplicitHaloExchange(1, comp);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::TmrMark
//! \brief implicit_timers: fence, charge the time since the last mark to category c
//! (c < 0 only sets the mark)

void RadiationM1::TmrMark(int c) {
  if (!tmr_on) return;
  Kokkos::fence();
  const double t = tmr_t.seconds();
  if (c >= 0) {tmr_acc[c] += t - tmr_last;}
  tmr_last = t;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::EventTotals
//! \brief fixbundle-1009 F4: the solver / limiter totals of this run (or restart) for the
//! event log, all ranks (collective).  v[0..6] are decided on global reductions and are
//! the same on every rank; v[7..12] are per-rank device counters, summed here.
//!   0 BiCGStab breakdowns, 1 line-Jacobi fallbacks, 2 fallbacks that kept the Krylov
//!   iterate (implicit_bcg_fallback = best), 3 APPLIED gas-Newton fallbacks,
//!   4 solves switched to the gas root find, 5 od positivity drops, 6 vet_col_lat
//!   positivity drops, 7 opac-Newton face terms dropped by the guard, 8 E floor raises
//!   (write-back, cell-solves), 9 energy the floor CREATED (x volume), 10 of which the
//!   solve undershoot (implicit_pos_floor_solve), 11 ApplyClosureLimits E raises,
//!   12 ApplyClosureLimits |F| clips (active cells, every stage)

void RadiationM1::EventTotals(Real *v) {
  v[0] = bcg_nbreak;
  v[1] = bcg_nfall;
  v[2] = bcg_nkeep;
  v[3] = newt_nfb;
  v[4] = impl_gn_nsw;
  v[5] = od_nfall;
  v[6] = vlat_nfall;
  for (int q = 7; q < 13; ++q) {v[q] = 0.0;}
  if (opn_nskip_d.extent_int(0) > 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    v[7] = h(0);
  }
  if (pos_cnt_d.extent_int(0) >= M1_POS_N) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    v[8] = h(M1_POS_FLR);
    v[9] = h(M1_POS_FLR_UN);
    v[10] = h(M1_POS_FLRS);
  }
  if (acl_cnt_d.extent_int(0) == 2) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acl_cnt_d);
    v[11] = h(0);
    v[12] = h(1);
  }
#if MPI_PARALLEL_ENABLED
  Real g[6];
  MPI_Allreduce(&v[7], g, 6, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < 6; ++q) {v[7+q] = g[q];}
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitReport
//! \brief one line at the end of the run with the solver statistics

void RadiationM1::ImplicitReport() {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;
  if ((impl_lin_tolcap > 0.0 || impl_lin_sktrue > 0.0) &&
      global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_lin_tolcap=" << impl_lin_tolcap
              << " tolcap_last=" << impl_lin_tclast << " capped passes (rank 0)="
              << impl_tclast_n << " stage=" << impl_lin_tcstage << " skip_true="
              << impl_lin_sktrue << " skipped checks (rank 0)=" << impl_sktrue_n
              << std::endl;
  }
  if ((impl_cons_final || dbg_cons_ledger) && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_cons_final=" << impl_cons_final
              << " solves=" << cf_tot[2] << " extra passes after a mode fallback="
              << cf_tot[3];
    if (dbg_cons_ledger) {
      std::cout << "; run sum V (c/chat) r=" << cf_tot[0] << " sum |.|=" << cf_tot[1]
                << (impl_cons_final ? " (applied to E: removed from the leak)"
                                    : " (NOT applied: the solves' energy error)");
    }
    std::cout << std::endl;
  }
  if (cs_geom && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> cubed sphere: seam face averages " << cs_seam_avg_n
              << (cs_seam_avg_on ? "" : " (disabled by rad_m1/cs_seam_avg)")
              << ", largest change of an x2 face value (rank 0) " << cs_seam_dmax
              << " of max|F0|" << std::endl;
  }
  // implicit_opac_newton_guard: the dropped (row, face) pairs of all ranks (collective:
  // every rank reaches this line with the same impl_opac_newton and guard)
  const bool opgr = impl_opac_newton && (impl_opn_guard > 0.0);
  if (opgr) {
    auto hsk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    opn_nskip = hsk(0);
#if MPI_PARALLEL_ENABLED
    Real g = 0.0;
    MPI_Allreduce(&opn_nskip, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    opn_nskip = g;
#endif
  }
  if (impl_g0_lim > 0.0 || impl_g0_exch || impl_pos_gas || impl_pos_floor) {
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = hpc(q);}
#if MPI_PARALLEL_ENABLED
    Real g[M1_POS_N];
    MPI_Allreduce(pos_cnt, g, M1_POS_N, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = g[q];}
    Real gd = 0.0;
    MPI_Allreduce(&pos_un_drop, &gd, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    pos_un_drop = gd;
#endif
  }
  // s3c: the hydro inversions' floor energy (EventCounters flr_de_v), all ranks
  Real hfl = pmy_pack->pmesh->ecounter.flr_de_v;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &hfl, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  // the Picard iteration count is MPI_MAX-reduced every step (see ImplicitSolve), so
  // every rank holds the same three numbers and no reduction is needed here
  // implicit_am_cs_kecap: the capped am_cs deposits
  if (impl_am_kecap > 0.0 && amkc_d.extent_int(0) == 4) {
    auto ha = Kokkos::create_mirror_view_and_copy(HostMemSpace(), amkc_d);
    Real av[4] = {ha(0), ha(1), ha(2), ha(3)};
#if MPI_PARALLEL_ENABLED
    Real ga[4];
    MPI_Allreduce(av, ga, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < 4; ++q) {av[q] = ga[q];}
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_am_cs_kecap=" << impl_am_kecap
                << ": capped cell-solves=" << av[0] << " kinetic energy not given="
                << av[1] << " (code units x volume) |L_z| not deposited=" << av[2]
                << " of all deposits " << av[3] << " (all attempts)" << std::endl;
    }
  }
  // amcsE-1111 implicit_am_cs_energy = rad: who paid the deposits' kinetic energy
  if (impl_am_en == 1 && amen_d.extent_int(0) == 4) {
    auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(), amen_d);
    Real ev[4] = {he(0), he(1), he(2), he(3)};
#if MPI_PARALLEL_ENABLED
    Real ge[4];
    MPI_Allreduce(ev, ge, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < 4; ++q) {ev[q] = ge[q];}
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_am_cs_energy=rad: cell-solves beyond the E share="
                << ev[0] << " KE paid by E=" << ev[1] << " by e_int=" << ev[2]
                << " not given=" << ev[3] << " (code units x volume, all attempts)"
                << std::endl;
    }
  }
  // implicit_am_cs_deposit = patch: the patch / column bookkeeping
  if (impl_am_dep == 1 && amvp_d.extent_int(0) == 8) {
    auto hp = Kokkos::create_mirror_view_and_copy(HostMemSpace(), amvp_d);
    Real pv[8];
    for (int q = 0; q < 8; ++q) {pv[q] = hp(q);}
#if MPI_PARALLEL_ENABLED
    Real gp[8];
    MPI_Allreduce(pv, gp, 8, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < 8; ++q) {pv[q] = gp[q];}
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_am_cs_deposit=patch: sum |del|=" << pv[0]
                << " capped cell-solves=" << pv[1] << " |L| moved to columns=" << pv[2]
                << " |L| not deposited=" << pv[3] << " singular patches=" << pv[4]
                << " column KE paid by E=" << pv[5] << " by e_int=" << pv[6]
                << " column cells beyond the e_int cap=" << pv[7]
                << " (code units x volume, all attempts)" << std::endl;
    }
  }
  // implicit_clip_kecap: the capped clip cell-solves
  if (impl_clip_kecap > 0.0 && clkc_d.extent_int(0) == 3) {
    auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), clkc_d);
    Real cv[3] = {hc(0), hc(1), hc(2)};
#if MPI_PARALLEL_ENABLED
    Real gc[3];
    MPI_Allreduce(cv, gc, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < 3; ++q) {cv[q] = gc[q];}
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_clip_kecap=" << impl_clip_kecap
                << ": capped clip cell-solves=" << cv[0] << " kinetic energy not given="
                << cv[1] << " |dm| not given=" << cv[2] << " (x volume, all attempts)"
                << std::endl;
    }
  }
  // implicit_wb_vcap: the capped cell-solves and the energy returned to E
  if (impl_wb_vcap > 0.0 && vcap_cnt_d.extent_int(0) == 4) {
    auto hv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vcap_cnt_d);
    Real vv[3] = {hv(0), hv(1), hv(3)};
#if MPI_PARALLEL_ENABLED
    Real gv[3];
    MPI_Allreduce(vv, gv, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    vv[0] = gv[0]; vv[1] = gv[1]; vv[2] = gv[2];
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> implicit_wb_vcap=" << impl_wb_vcap
                << (impl_wb_vcap_meas ? " (measure only)" : "")
                << ": cell-solves beyond it=" << vv[0]
                << " kinetic energy "
                << (impl_wb_vcap_meas ? "beyond it=" : "taken from the gas=") << vv[1]
                << " of which not returned to E (bound E)=" << vv[2]
                << " (code units x volume, all attempts)" << std::endl;
    }
  }
  if (impl_wb_vkrej > 0.0 && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> implicit_wb_vkick_reject=" << impl_wb_vkrej
              << ": stages rejected=" << vkick_nrej
              << " backward-Euler solves with kicks=" << vkick_nbe << " (cells "
              << vkick_nbecell << ")";
    if (impl_wb_vkbr) {
      std::cout << " of which dt_retry=" << vkick_nberetry
                << " accepted after the last cut=" << vkick_nbefall;
    }
    std::cout << std::endl;
  }
  if (impl_fclip_diag && fcd_cnt_d.extent_int(0) == 20) {
    auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fcd_cnt_d);
    Real fv[20];
    for (int q = 0; q < 20; ++q) {fv[q] = hf(q);}
#if MPI_PARALLEL_ENABLED
    Real gf[20];
    MPI_Allreduce(fv, gf, 20, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < 20; ++q) {fv[q] = gf[q];}
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> FCLIPDIAG in-loop (cs2_clip=rad, unpaired): n=" << fv[0]
                << " mean(|u|/|q|-1)=" << fv[1]/fmax(fv[0], 1.0)
                << " face-mean>cE n=" << fv[2]
                << " of which face-mean<=cE^old n=" << fv[3]
                << " mean(cE/cE^old)=" << fv[4]/fmax(fv[0], 1.0)
                << " unpaired |dF|/(c chat) x vol=" << fv[5] << std::endl;
      std::cout << "<rad_m1> FCLIPDIAG final: n(metric)=" << fv[6]
                << " mean(|F|/cE-1)=" << fv[7]/fmax(fv[6], 1.0)
                << " n(plain norm)=" << fv[8] << " made by pos moves=" << fv[9]
                << " made by cons_final=" << fv[10]
                << " E lowered by pos moves=" << fv[11]
                << " removed (|F|-cE)/(c chat) x vol=" << fv[12]
                << " paired (consrc kick) n=" << fv[13]
                << " face-mean>cE(pre-move) n=" << fv[14] << " E<E^old/2 n=" << fv[15]
                << std::endl;
      std::cout << "<rad_m1> FCLIPDIAG implicit_realisable_coupling rescale: lag (beta.F"
                << " of E0) cell-passes=" << fv[16] << " mean(|F|/cE-1)="
                << fv[17]/fmax(fv[16], 1.0) << " | write-back drag cell-solves=" << fv[18]
                << " mean(|F|/cE-1)=" << fv[19]/fmax(fv[18], 1.0) << std::endl;
    }
  }
  if (global_variable::my_rank != 0) return;
  if (tmr_c0 > 0) {
    // implicit_timers: seconds per category, then the counters (see TmrMark)
    std::cout << "<rad_m1> timers from cycle " << tmr_c0 << ": stages=" << tmr_cnt[7]
              << " solves=" << tmr_cnt[0] << " passes=" << tmr_cnt[1]
              << " kry_it=" << tmr_cnt[2] << " passes_s1=" << tmr_cnt[3]
              << " passes_s2=" << tmr_cnt[4] << " it_s1=" << tmr_cnt[5]
              << " it_s2=" << tmr_cnt[6] << std::endl;
    const char *nm[12] = {"closure", "opacity", "solve_pre", "tensor", "pred",
                          "pass_setup", "krylov", "pass_post", "solve_end", "hyd_c2p",
                          "m1_bvals", "mg_x3"};
    std::cout << "<rad_m1> timers(s):";
    // mg_x3 = the x3 line sweeps (out of krylov) under implicit_timers_mg_x3
    const int ntm = (mg_x3line && tmr_x3) ? 12 : 11;
    for (int q = 0; q < ntm; ++q) {std::cout << " " << nm[q] << "=" << tmr_acc[q];}
    std::cout << std::endl;
  }
  MRReport();
  Real mean = (impl_nstep > 0.0) ? (impl_itsum/impl_nstep) : 0.0;
  std::cout << "<rad_m1> implicit transport: solves=" << impl_nstep
            << " Picard iterations mean=" << mean << " max=" << impl_itmax
            << " NON-CONVERGED=" << impl_nfail << std::endl;
  if (impl_opac_lagp) {
    std::cout << "<rad_m1> implicit_opac_lag=pred passes=" << impl_opac_lagn
              << ": lagged solves=" << impl_ol_nlag << std::endl;
  }
  if (impl_tol_rel > 0.0) {
    std::cout << "<rad_m1> implicit_tol_rel=" << impl_tol_rel << ": predicted solves "
              << impl_tolr_n << ", mean tol_eff "
              << ((impl_tolr_n > 0.0) ? impl_tolr_sum/impl_tolr_n : 0.0) << std::endl;
  }
  if (impl_frlog) {
    // implicit_final_resid_log: the per-run sums of ImplicitFinalResid
    const Real *S = impl_fr_sum;
    std::cout << "<rad_m1> implicit_final_resid_log: solves=" << S[0]
              << " max rc=" << S[14] << " max rT=" << S[15]
              << " | solves with max rc > tol " << S[1] << ", > 10 tol " << S[2]
              << ", > 1e-6 " << S[3] << ", > 1e-4 " << S[4] << ", > 1e-2 " << S[5]
              << " | max rT > tol " << S[6] << ", > 10 tol " << S[7]
              << " | max(rc,rT) > tol " << S[19] << ", > 10 tol " << S[8]
              << " | conv_est stops " << S[16] << " (of them > 10 tol " << S[17] << ")"
              << " | cells > 10 tol " << S[12] << " of " << S[13]
              << " | no gas row " << S[18]
              << " | time loop " << S[11] << " s, shadow eval " << S[9]
              << " s, snapshot " << S[10] << " s" << std::endl;
  }
  if (impl_frs_on || impl_ce_guard > 0.0 || impl_ce_nofl) {
    // the guarded contraction estimate and the residual stop
    const Real *P = impl_frs_sum;
    std::cout << "<rad_m1> picstop: implicit_final_resid_stop=" << impl_frs << " (rel "
              << impl_frs_rel << " x tol_eff) checks=" << P[0] << " failed+continued="
              << P[1] << " failed-no-pass-left(NC)=" << P[2]
              << " solves accepted after a failed check=" << P[3]
              << " passes after a failed check=" << P[4]
              << " | implicit_conv_est_guard=" << impl_ce_guard << " refused=" << P[5]
              << " | implicit_conv_est_nofloor=" << (impl_ce_nofl ? 1 : 0)
              << " refused=" << P[6] << " extra stops=" << P[7]
              << " | log worst-cell time " << P[8] << " s" << std::endl;
  }
  if (impl_gn_sw > 0) {
    std::cout << "<rad_m1> implicit_gas_newton_switch=" << impl_gn_sw << " (min pass "
              << impl_gn_sw_min << "): solves switched to the root find=" << impl_gn_nsw
              << std::endl;
  }
  if (impl_res_mask) {
    std::cout << "<rad_m1> implicit_res mask: rho < " << impl_res_dmin << " or r > "
              << impl_res_rmax << " (0 = unused) left out of the stopping test; max"
              << " excluded residual at the last pass=" << impl_res_exmax
              << ", steps whose excluded residual >= tol=" << impl_res_nex << std::endl;
  }
  if (opgr) {
    std::cout << "<rad_m1> implicit_opac_newton_guard=" << impl_opn_guard
              << ": Newton face terms dropped (all ranks, all passes)=" << opn_nskip
              << " per solve=" << ((impl_nstep > 0.0) ? (opn_nskip/impl_nstep) : 0.0)
              << std::endl;
  }
  if (impl_g0_lim > 0.0 || impl_g0_exch || impl_pos_gas || impl_pos_floor) {
    std::cout << "<rad_m1> m1-positivity: implicit_g0_exchange=" << impl_g0_exch
              << " implicit_g0_limit=" << impl_g0_lim
              << " clipped cell-passes=" << pos_cnt[M1_POS_G0]
              << " | implicit_pos_gas=" << impl_pos_gas << " (frac "
              << impl_pos_gas_frac << ") cell-solves=" << pos_cnt[M1_POS_GAS]
              << " energy moved rad->gas=" << pos_cnt[M1_POS_GAS_DE]
              << " | E floor raises cell-solves=" << pos_cnt[M1_POS_FLR]
              << " energy from gas=" << pos_cnt[M1_POS_FLR_DE]
              << " energy created=" << pos_cnt[M1_POS_FLR_UN] << " (code units x volume)"
              << " [in discarded attempts, rolled back: " << pos_un_drop << "]"
              << " of which solve undershoot (pos_floor_solve)=" << pos_cnt[M1_POS_FLRS]
              << " taken from the column's radiation (implicit_pos_floor_cons)="
              << pos_cnt[M1_POS_FLR_MV]
              << " | |F| > c E scaled back cell-solves=" << pos_cnt[M1_POS_FCLIP]
              << " sum(|F|/(cE) - 1)=" << pos_cnt[M1_POS_FCLIPM]
              << " | hydro EOS floors energy created=" << hfl << std::endl;
  }
  if (impl_fmm) {   // (all ranks, all passes)
    std::cout << "<rad_m1> implicit_face_mmat: checked passes=" << fmm_pass
              << " limited passes=" << fmm_lpass << " limited cells=" << fmm_cells
              << " limited x1 faces=" << fmm_faces << " max rounds=" << fmm_rmax
              << " left (b<=0 after the limit)=" << fmm_left
              << " realisability-bound faces (|F_H| > c min E)=" << fmm_rfaces
              << " (rank 0: rounds=" << zl_rounds[1] << " R exchanges=" << zl_halos[1]
              << " zalesak_halo_skip=" << (zl_hskip ? "on" : "off")
              << (zl_local > 0 ? " zalesak_local_rounds=" + std::to_string(zl_local)
                  : std::string(""))
              << " rowsel=" << (fmm_rowsel ? "on" : "off") << " re-assembled rows="
              << fmm_selrows << " of " << fmm_selall;
    if (zl_prof) {std::cout << " second-assembly wall=" << fmm_asm2_wall << " s";}
    if (zl_prof) {
      std::cout << " round-loop wall=" << zl_wall[1] << " s (exchanges " << zl_xwall[1]
                << " s)";
    }
    if (zl_active) {std::cout << " zalesak_active=on";}
    if (zl_colsweep) {std::cout << " zalesak_colsweep=on";}
    if (zl_cap > 0) {
      std::cout << " zalesak_cap_rounds=" << zl_cap << " capped passes=" << zl_capped;
    }
    if (zl_local > 0) {std::cout << " local sub-rounds=" << zl_sub[1];}
    std::cout << ")" << std::endl;
  }
  if (impl_onep > 0) {
    std::cout << "<rad_m1> implicit_one_pass: period=" << impl_onep
              << " safety=" << impl_onep_s
              << " solves accepted after one pass=" << impl_onep_n
              << " contraction measurements=" << impl_onep_nchk
              << " last q (be, stage 1, stage 2)=" << onep_qa[0] << " " << onep_qa[1]
              << " " << onep_qa[2] << std::endl;
    if (impl_onep_auto) {
      // the counters are those of this run (from its start or restart); the state is
      // per kind: on, off, or on as a re-probe not yet confirmed
      std::cout << "<rad_m1> implicit_one_pass_auto: window=" << impl_onep_awin
                << " reprobe=" << impl_onep_arep << "; this run: switch-offs="
                << impl_onep_ndis << " re-probes=" << impl_onep_nprb
                << " (failed=" << impl_onep_nfpr << " confirmed=" << impl_onep_nren
                << "); state at end (be, stage 1, stage 2)=";
      for (int t = 0; t < 3; ++t) {
        const char *st = (onep_off[t] > 0.5) ? "off"
                         : ((onep_prb[t] > 0.5) ? "probe" : "on");
        std::cout << st << ((t < 2) ? "," : "");
      }
      std::cout << std::endl;
    }
  }
  if (impl_accel == M1_IACC_ANDERSON) {
    Real apst = (impl_nstep > 0.0) ? (aa_nacc/impl_nstep) : 0.0;
    std::cout << "<rad_m1> anderson: m=" << impl_and_m
              << " beta=" << impl_and_beta
              << " start=" << impl_and_start
              << " accelerated passes=" << aa_nacc
              << " (" << apst << " per solve)"
              << " history restarts=" << aa_nrst << std::endl;
  }
  if (trans_on) {
    Real lmean = (impl_nstep > 0.0) ? (impl_linsum/impl_nstep) : 0.0;
    std::cout << "<rad_m1> implicit transverse ("
              << (bicg_on ? "bicgstab" : "line_jacobi")
              << "): final 7-point linear "
              << "residual (max|r|/max|b|) mean=" << lmean << " max=" << impl_linmax
              << std::endl;
    // say which test stopped the solves.  The max|r|/max|b| above is
    // compared with implicit_lin_tol only when implicit_lres_test is on (default off for
    // the fixed closures: eddington, vet_sc, tau / vet_col); the inner BiCGStab stops on
    // implicit_lin_cnorm (per-cell |r|/(s E)) when > 0, else on |r|/max|b| < lin_tol,
    // loosened by the Eisenstat-Walker term when implicit_lin_ew_max > 0
    std::cout << "<rad_m1>   (diagnostic; Picard test on it: "
              << (impl_lres_test ? "yes, lin_tol=" : "no; lin_tol=") << impl_lin_tol
              << "; inner stop: "
              << ((impl_lin_cnorm > 0.0) ? "cnorm=" : "lin_tol relative to max|b|");
    if (impl_lin_cnorm > 0.0) {std::cout << impl_lin_cnorm;}
    if (impl_ew_max > 0.0) {
      std::cout << ", loosened by Eisenstat-Walker ew_max=" << impl_ew_max;
    }
    std::cout << ")" << std::endl;
  }
  if (bicg_on) {
    Real imean = (bcg_nsolve > 0.0) ? (bcg_itsum/bcg_nsolve) : 0.0;
    Real rper = (bcg_itsum > 0.0) ? (bcg_nred/bcg_itsum) : 0.0;
    std::cout << "<rad_m1> bicgstab: outer passes=" << impl_itsum
              << " linear solves=" << bcg_nsolve
              << " inner iterations mean=" << imean << " max=" << bcg_itmax
              << " total=" << bcg_itsum << std::endl;
    std::cout << "<rad_m1> bicgstab: breakdowns=" << bcg_nbreak
              << " line_jacobi fallbacks=" << bcg_nfall
              << " global reductions=" << bcg_nred
              << " (" << rper << " per inner iteration)" << std::endl;
    if (impl_bcg_keep || impl_bcg_maxrst != 2) {
      std::cout << "<rad_m1> bicgstab: implicit_bcg_max_restarts=" << impl_bcg_maxrst
                << " fallback=" << (impl_bcg_keep ? "best" : "line")
                << ": Krylov iterate kept instead of line-Jacobi=" << bcg_nkeep
                << std::endl;
    }
    if (impl_prec_fb >= 0) {
      std::cout << "<rad_m1> bicgstab: implicit_precond_fallback="
                << ((impl_prec_fb == 1) ? "rbgs" : "line")
                << ": solves capped at lin_maxit and redone=" << bcg_npfb
                << " their inner iterations=" << bcg_pfbit << std::endl;
    }
    if (impl_kdev > 0) {
      std::cout << "<rad_m1> implicit_krylov_dev: period=" << impl_kdev
                << " host status reads=" << kdev_nchk << " queued iterations="
                << kdev_nq << std::endl;
    }
  }
  if (impl_gas_newton || impl_eos_cache) {
    Real fpc = (gas_ncell > 0.0) ? (newt_nfb/gas_ncell) : 0.0;
    Real mpc = (gas_ncell > 0.0) ? (ec_nmiss/gas_ncell) : 0.0;
    std::cout << "<rad_m1> gas coupling: newton="
              << (impl_gas_newton ? "true" : "false")
              << " eos_cache=" << (impl_eos_cache ? "true" : "false")
              << " cell-passes=" << gas_ncell
              << " newton fallbacks=" << newt_nfb << " (" << fpc << " per cell-pass)"
              << std::endl;
    if (impl_eos_cache) {
      std::cout << "<rad_m1> eos_cache: nt=" << impl_ecnt
                << " misses=" << ec_nmiss << " (" << mpc << " per cell-pass)"
                << " max |de|/e vs the table=" << ec_emax
                << " max |dq|/q of the exchanged energy=" << ec_tmax << std::endl;
    }
  }
  if (impl_pcr_check) {
    std::cout << "<rad_m1> line solver=" << ((impl_line_solver == 1) ? "pcr" : "thomas")
              << " pcr_check: calls=" << pcr_chk_n
              << " max over calls of max|x_pcr - x_thomas|/max|x|=" << pcr_chk_max
              << " (rank 0)" << std::endl;
  }
  if (trans_on) {
    std::cout << "<rad_m1> offdiag="
              << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                  ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
              << " closure_relax=" << impl_crelax
              << " closure_lag=" << (impl_clag_step ? "step" : "pass")
              << " positivity fallbacks=" << od_nfall
              << " (vet_col_lat " << vlat_nfall << ")"
              << " min E from the solve=" << od_emin << std::endl;
  }
  if (vgd_on && vgd_tan && cs_geom) {
    std::cout << "<rad_m1> vet_gd_tangential (cubed sphere): Q active in " << cstan_na
              << " of " << cstan_np << " transverse assemblies (rank 0); positivity "
              << "fallbacks (vet_col_lat) " << vlat_nfall << std::endl;
  }
  if (vlat_posfix_lag) {
    std::cout << "<rad_m1> vet_col_lat_posfix = lag: triggers moved to the lagged "
              << "lateral form " << vlat_nlag << ", dropped after a second trigger "
              << vlat_nfall
              << std::endl;
  }
  if (impl_vimp) {
    std::cout << "<rad_m1> implicit_vimp positivity fallbacks=" << vimp_nfall
              << " min E from the solve=" << vimp_emin;
    if (vimp_kill) {
      std::cout << " local kills=" << vkill_ev << " cells(rank 0)=" << vkill_nc;
    }
    std::cout << std::endl;
  }
  if (vimp_lag_local) {
    std::cout << "<rad_m1> implicit_vimp_lag_local: lag-path triggers with the vimp kick"
              << " dropped locally (not for the whole solve)=" << vimp_lag_local_n
              << std::endl;
  }
  if (t2_dtr_pl) {
    std::cout << "<rad_m1> time2_dtr_poslim: BE-type solves with written-back eint <= 0"
              << " all covered by implicit_pos_gas (no dt_retry cut)=" << t2_dtr_pl_n
              << std::endl;
  }
  if (impl_pos_trans != 0 && pt_cnt.extent_int(0) >= 2) {
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pt_cnt);
    std::cout << "<rad_m1> implicit_pos_transport = upwind (rank 0):"
              << " lagged positive x1 off-diagonals=" << h(1)
              << " memory-upwinded x1 face passes=" << h(0) << std::endl;
  }
  if (impl_cmi && cmi_cnt.extent_int(0) >= 3) {
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cmi_cnt);
    std::cout << "<rad_m1> implicit_hr_cmem_implicit (rank 0): cell-passes with an "
              << "explicit memory outflow=" << h(1) << " limited (lam < 1)=" << h(0)
              << " min lam=" << h(2) << std::endl;
  }
  if (impl_muscl && impl_muscl_pred) {
    std::cout << "<rad_m1> implicit_hr_recon_pred: dc-phase passes=" << muscl_pred_dc
              << " plm-phase passes=" << muscl_pred_hi << std::endl;
  }
  if (impl_muscl && impl_muscl_bound) {
    std::cout << "<rad_m1> implicit_hr_pos = bound: cell-axes set to dc=" << muscl_nbound
              << std::endl;
  }
  if (impl_muscl) {
    std::cout << "<rad_m1> implicit_hr_recon = plm positivity fallbacks=" << muscl_nfall
              << " min E from the solve=" << muscl_emin << std::endl;
  }
  std::cout << "<rad_m1> floor clips (all ranks, every solve incl. BE): solved E <= "
            << "e_floor cell-solves=" << flr_ne
            << " written-back gas eint <= 0 cell-solves="
            << flr_ng << std::endl;
}
//----------------------------------------------------------------------------------------
//! implicit_hr_cmem_rst: fax1..3 to / from cell-sized restart slabs
//! (rad_m1.hpp).  Slab 2d: face d indices 0..n-1; slab 2d+1: face n in the first plane.

void RadiationM1::FaxRstPack(DvceArray5D<Real> out, int nmb) {
  const int n3 = out.extent_int(2), n2 = out.extent_int(3), n1 = out.extent_int(4);
  const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
  const int s0 = FaxRstNslab() - ((t2_pfst && t2pf1_ok) ? nd : 0);
  // (time2_face_fct: t2pf1 behind the fax slabs, one per direction, files before two,
  // the second 0; the slabs include the ghosts, so every face a solve set)
  if (s0 < FaxRstNslab()) {
    auto pf = t2pf1;
    par_for("pf_rst_pack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < nd; ++d) out(m,s0+d,k,j,i) = pf(m,d,k,j,i);
    });
  }
  if (s0 == 0) {return;}
  const bool h2 = fax2.extent_int(0) > 0, h3 = fax3.extent_int(0) > 0;
  auto f1 = fax1, f2 = fax2, f3 = fax3;
  par_for("fax_rst_pack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    out(m,0,k,j,i) = f1(m,k,j,i);
    out(m,1,k,j,i) = (i == 0) ? f1(m,k,j,n1) : 0.0;
    if (h2) {
      out(m,2,k,j,i) = f2(m,k,j,i);
      out(m,3,k,j,i) = (j == 0) ? f2(m,k,n2,i) : 0.0;
    }
    if (h3) {
      const int s = h2 ? 4 : 2;
      out(m,s,k,j,i) = f3(m,k,j,i);
      out(m,s+1,k,j,i) = (k == 0) ? f3(m,n3,j,i) : 0.0;
    }
  });
}

bool RadiationM1::FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl) {
  const int n3 = in.extent_int(2), n2 = in.extent_int(3), n1 = in.extent_int(4);
  // time2_face_fct: t2pf1 behind the fax slabs (a file without them: the first stage
  // takes the stage-start face F, Time2FaceFct)
  const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
  const bool cm = impl_hr_cmem && fax1.extent_int(0) >= nmb;
  const int na = cm ? 2 + ((fax2.extent_int(0) > 0) ? 2 : 0) +
                      ((fax3.extent_int(0) > 0) ? 2 : 0) : 0;
  if (t2_pfst && t2pf1.extent_int(0) >= nmb && (nsl == na + nd || nsl == na + 2*nd)) {
    auto pf = t2pf1;
    const int st = (nsl == na + nd) ? 1 : 2;
    par_for("pf_rst_unpack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < nd; ++d) pf(m,d,k,j,i) = in(m,na+st*d,k,j,i);
    });
    t2pf1_ok = true;
    nsl = na;
  } else if (t2_pfst && global_variable::my_rank == 0) {
    std::cout << "NOTE <rad_m1> restart: no time2_face_fct face fluxes in the file; the "
              << "first stage takes the stage-start face F" << std::endl;
  }
  if (!cm || nsl == 0) {return nsl == 0;}
  const bool h2 = fax2.extent_int(0) > 0, h3 = fax3.extent_int(0) > 0;
  if (nsl != 2 + (h2 ? 2 : 0) + (h3 ? 2 : 0)) {return false;}
  auto f1 = fax1, f2 = fax2, f3 = fax3;
  par_for("fax_rst_unpack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f1(m,k,j,i) = in(m,0,k,j,i);
    if (i == 0) {f1(m,k,j,n1) = in(m,1,k,j,0);}
    if (h2) {
      f2(m,k,j,i) = in(m,2,k,j,i);
      if (j == 0) {f2(m,k,n2,i) = in(m,3,k,0,i);}
    }
    if (h3) {
      const int s = h2 ? 4 : 2;
      f3(m,k,j,i) = in(m,s,k,j,i);
      if (k == 0) {f3(m,n3,j,i) = in(m,s+1,0,j,i);}
    }
  });
  impl_cm_fainit_done = true;   // restored exactly: no fainit
  return true;
}

} // namespace radm1
