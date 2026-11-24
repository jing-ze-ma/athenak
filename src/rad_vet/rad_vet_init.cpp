//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_init.cpp
//! \brief set-up of the implicit radiation (<rad_vet>/transport = implicit): the keys the
//! clean steps read (ImplicitInit), the x1 partition (ImplicitPartitionInit), the ghost
//! exchanges of the work array (ImplicitHalo*), the run-end report (ImplicitReport,
//! EventTotals) and the face-flux restart slabs (FaxRst*).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <string>
#include <vector>

#include "rad_vet/rad_vet_t2reset.hpp"
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "coordinates/coordinates.hpp"
#include "utils/recipe.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {
// reset-sector exemption of the hesdirk2 admissibility (rad_vet_t2reset.hpp)
bool g_t2_reset_on = false;
DvceArray4D<Real> g_t2_reset;
double g_t2_reset_nskip = 0.0;
static int M1HaloDepth(MeshBlockPack *pp, ParameterInput *pin);   // (below)

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
    << std::endl << "<rad_vet>/implicit_bc_x1min|max = '" << s << "' is not a valid "
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
//! \fn void RadiationVet::ImplicitInit
//! \brief read the <rad_vet> parameters of the implicit transport that the clean steps
//! use (time scheme, boundaries, predictor, gas source, positivity, halos), check them
//! and allocate the face and work arrays.  Called at the END of the constructor.

void RadiationVet::ImplicitInit(ParameterInput *pin) {
  // <rad_vet>/time_scheme (read by Time2Init below): hesdirk2 on a fresh run wherever it
  // is accepted (no polar boundary, a closure other than tau (vet_col included),
  // <time>/integrator = rk2), be elsewhere and on a restart whose file lacks the key.
  // The resolved value is echoed (so later restarts keep it).
  {auto *pmh = pmy_pack->pmesh;
  const std::string integ = pin->DoesParameterExist("time","integrator") ?
                            pin->GetString("time","integrator") : "rk2";
  const bool h2ok = !pmh->use_polar_boundary && (!tau_closure || vet_col) &&
                    (integ.compare("rk2") == 0);
  (void) pin->GetOrAddString("rad_vet","time_scheme",
                             (h2ok && !global_variable::restart_run) ? "hesdirk2" : "be");
  // plain be = the hesdirk2 backward-Euler path (t2_beonly, rad_vet.hpp) wherever
  // hesdirk2 is accepted; <rad_vet>/be_legacy = true keeps the old step (the default on a
  // restart whose file lacks it).  time2_be_only = true: the same path under hesdirk2.
  const std::string tsn = pin->GetString("rad_vet","time_scheme");
  if (tsn.compare("be") == 0) {
    const bool leg = pin->GetOrAddBoolean("rad_vet","be_legacy",
                                          global_variable::restart_run);
    t2_beonly = h2ok && !leg;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> time_scheme = be: be_legacy = " << (leg ? "true" : "false")
                << (t2_beonly ? " (the hesdirk2 backward-Euler path)" : " (the old step)")
                << std::endl;
    }
  } else if (tsn.compare("hesdirk2") == 0 &&
             pin->DoesParameterExist("rad_vet","time2_be_only")) {
    t2_beonly = pin->GetBoolean("rad_vet","time2_be_only");
  }}
  trans_on = false;
  trans_x3 = false;
  impl_cfl = pin->GetOrAddReal("rad_vet","implicit_cfl",-1.0);
  // implicit_solver: bicgstab under the new default set (utils/recipe.hpp); a restart
  // whose input lacks the key / recipe none: line_jacobi (polar_ring_avg needs bicgstab)
  if (!pin->DoesParameterExist("rad_vet","implicit_solver")) {
    const bool nd = recipe::NewDefaults(pin);
    recipe::Log(std::string("<rad_vet>/implicit_solver = ") +
                (nd ? "bicgstab (default)" : "line_jacobi (the old default)"));
    pin->GetOrAddString("rad_vet","implicit_solver", nd ? "bicgstab" : "line_jacobi");
  }
  {std::string sv = pin->GetOrAddString("rad_vet","implicit_solver","line_jacobi");
  if (sv.compare("line_jacobi") == 0) {
    impl_solver = M1_ISOLV_LINE_JACOBI;
  } else if (sv.compare("bicgstab") == 0) {
    impl_solver = M1_ISOLV_BICGSTAB;
  } else {
    ImplFatal("<rad_vet>/implicit_solver = '" + sv
              + "' is not a choice (line_jacobi | bicgstab)");
  }
  }
  // implicit_offdiag (the spherical-polar checks read it too): auto = operator with
  // bicgstab on a multi-D mesh, lagged otherwise
  std::string sod = pin->GetOrAddString("rad_vet","implicit_offdiag","auto");
  bool od_auto = (sod.compare("auto") == 0);
  if (od_auto || sod.compare("lagged") == 0) {
    impl_offdiag = M1_OD_LAGGED;
  } else if (sod.compare("operator") == 0) {
    impl_offdiag = M1_OD_OPERATOR;
  } else if (sod.compare("none") == 0) {
    impl_offdiag = M1_OD_NONE;
  } else {
    ImplFatal("<rad_vet>/implicit_offdiag = '" + sod
              + "' is not a choice (auto | lagged | operator | none)");
  }
  // implicit_closure_lag: step on a multi-D mesh, pass elsewhere and on a restart whose
  // file lacks the key (recorded).  With step, a chi(f) closure (m1, minerbo, kershaw)
  // relaxes optically thin cells step to step by implicit_closure_thin_relax (default
  // 1.5; 0 on a restart whose file lacks it): its memory ctr_mem travels in the restart
  // file.  The fixed-tensor closures (eddington, vet_sc, tau, vet_col) never relax.
  const bool md_full = pmy_pack->pmesh->multi_d;
  bool clag_step = false;
  {std::string sc = pin->GetOrAddString("rad_vet","implicit_closure_lag",
      (md_full && !global_variable::restart_run) ? "step" : "pass");
  if (sc.compare("step") == 0) {
    clag_step = true;
  } else if (sc.compare("pass") != 0) {
    ImplFatal("<rad_vet>/implicit_closure_lag = '" + sc
              + "' is not a choice (pass | step)");
  }
  }
  impl_ctrelax = 0.0;
  ctr_init = false;
  {const bool chif = !(eddington || vet_sc || tau_closure);
  if (chif && md_full && clag_step) {
    impl_ctrelax = pin->GetOrAddReal("rad_vet","implicit_closure_thin_relax",
                                     global_variable::restart_run ? 0.0 : 1.5);
  } else if (pin->DoesParameterExist("rad_vet","implicit_closure_thin_relax")) {
    impl_ctrelax = pin->GetReal("rad_vet","implicit_closure_thin_relax");
    if (!chif) {impl_ctrelax = 0.0;}
  }
  if (impl_ctrelax < 0.0) {
    ImplFatal("<rad_vet>/implicit_closure_thin_relax must be >= 0");
  }
  }
  // the transverse realizability limiter (none: the plain theta)
  {std::string st = pin->GetOrAddString("rad_vet","implicit_trans_limit","none");
  if (st.compare("none") == 0) {
    impl_tlim = M1_TLIM_NONE;
  } else if (st.compare("lp") == 0) {
    impl_tlim = M1_TLIM_LP;
  } else if (st.compare("lpr") == 0) {
    impl_tlim = M1_TLIM_LPR;
  } else {
    ImplFatal("<rad_vet>/implicit_trans_limit = '" + st
              + "' is not a choice (none | lp | lpr)");
  }
  }
  dbg_etally = pin->GetOrAddBoolean("rad_vet","dbg_energy_tally",false);
  // the fast set of the fixed-tensor closures (eddington, vet_sc, tau) with bicgstab and
  // ONE MeshBlock along x1: implicit_halo_direct (the on-rank exchanges as one copy
  // kernel when every neighbour is on its rank at the same level; bitwise),
  // implicit_krylov_fuse (3; 1 with implicit_lin_cnorm > 0) and implicit_op_stencil (not
  // with a periodic x1 wrap).  A key the input names keeps its value.
  const bool fixcl = eddington || vet_sc || tau_closure;
  const bool fdef = fixcl && (impl_solver == M1_ISOLV_BICGSTAB) &&
                    (pmy_pack->pmesh->mesh_indcs.nx1 == pmy_pack->pmesh->mb_indcs.nx1);
  const bool cn_def = pin->DoesParameterExist("rad_vet","implicit_lin_cnorm") &&
                      pin->GetReal("rad_vet","implicit_lin_cnorm") > 0.0;
  // the x1 wrap as the boundary parsing below will read it (same keys, same defaults)
  const bool x1per = (ImplBCFromString(
      pin->GetOrAddString("rad_vet","implicit_bc_x1min","auto"),
      pmy_pack->pmesh->mesh_bcs[static_cast<int>(BoundaryFace::inner_x1)])
      == M1_IBC_PERIODIC);
  impl_halo_direct = pin->GetOrAddBoolean("rad_vet","implicit_halo_direct",fdef);
  impl_kfuse = pin->GetOrAddInteger("rad_vet","implicit_krylov_fuse",
                                    fdef ? (cn_def ? 1 : 3) : 0);
  // a restart file written without cnorm echoes the defaulted fuse 3: with cnorm named
  // now it runs fuse 1, echoed
  if (impl_kfuse == 3 && cn_def &&
      pin->IsParameterDefaulted("rad_vet","implicit_krylov_fuse")) {
    impl_kfuse = 1;
    pin->SetInteger("rad_vet","implicit_krylov_fuse",1);
    std::cout << "rad_vet: implicit_krylov_fuse 3 (defaulted in the restart file) -> 1, "
              << "since implicit_lin_cnorm > 0" << std::endl;
  }
  if (impl_kfuse < 0 || impl_kfuse > 3) {
    ImplFatal("<rad_vet>/implicit_krylov_fuse must be 0, 1, 2 or 3");
  }
  impl_stencil = pin->GetOrAddBoolean("rad_vet","implicit_op_stencil",fdef && !x1per);
  impl_onep = 0;
  for (int t = 0; t < 3; ++t) {
    onep_qa[t] = -1.0;
    onep_qb[t] = -1.0;
    onep_cnt[t] = 0.0;
    onep_off[t] = 0.0;
    onep_actr[t] = 0.0;
    onep_prb[t] = 0.0;
  }
  impl_onep_auto = false;
  impl_pord = 1;
  // multi-rank exchanges: implicit_halo_mpi (bitwise the ordinary exchange) wherever
  // implicit_halo_direct holds on a same-level mesh without cubed-sphere seams or polar
  // boundaries (on the spherical-polar wedge not on a restart whose file lacks the key);
  // implicit_halo_overlap on with more than one rank on a fresh run.  A key the input
  // (or restart echo) names keeps its value.
  {auto *pmh = pmy_pack->pmesh;
  const bool hmdef = impl_halo_direct && !pmh->multilevel && !pmh->use_cubed_sphere &&
                     !pmh->use_polar_boundary &&
                     !(pmh->use_spherical_polar && global_variable::restart_run);
  impl_halo_mpi = pin->GetOrAddBoolean("rad_vet","implicit_halo_mpi",hmdef);}
  hm_state = 0;
  hm_comm = nullptr;
  impl_halo_ipc = pin->GetOrAddBoolean("rad_vet","implicit_halo_ipc",false);
  if (impl_halo_ipc && !impl_halo_mpi) {
    ImplFatal("<rad_vet>/implicit_halo_ipc needs implicit_halo_mpi = true");
  }
  impl_halo_ovl = pin->GetOrAddBoolean("rad_vet","implicit_halo_overlap",
                      impl_halo_mpi && (global_variable::nranks > 1) &&
                      !global_variable::restart_run);
  impl_ovl_faces = pin->GetOrAddBoolean("rad_vet","implicit_halo_ovl_faces",false);
  for (int f = 0; f < 6; ++f) {hm_face[f] = 1;}
  if (impl_halo_ovl && !impl_halo_mpi) {
    ImplFatal("<rad_vet>/implicit_halo_overlap needs implicit_halo_mpi = true");
  }
  if (impl_halo_mpi && !impl_halo_direct) {
    ImplFatal("<rad_vet>/implicit_halo_mpi needs implicit_halo_direct = true");
  }
  tmr_c0 = pin->DoesParameterExist("rad_vet","implicit_timers") ?
           pin->GetInteger("rad_vet","implicit_timers") : 0;   // read only when named
  // implicit_resid_fatal: the residual above which a stage stops the run (default 1e-2
  // with the new defaults; 0 = off on a restart whose file lacks the key; 1e2 under
  // recipe none)
  if (pin->DoesParameterExist("rad_vet","implicit_resid_fatal")) {
    impl_res_fatal = pin->GetReal("rad_vet","implicit_resid_fatal");
    if (!(impl_res_fatal >= 0.0)) {
      ImplFatal("<rad_vet>/implicit_resid_fatal must be >= 0");
    }
  } else if (recipe::NewDefaults(pin)) {
    impl_res_fatal = pin->GetOrAddReal("rad_vet","implicit_resid_fatal",1.0e-2);
    recipe::Log("<rad_vet>/implicit_resid_fatal = 1e-2 (default)");
  } else if (global_variable::restart_run) {
    recipe::Log("<rad_vet>/implicit_resid_fatal = 0 (off: restart input has no key, the "
                "old default kept; set it to switch)");
  } else {
    impl_res_fatal = pin->GetOrAddReal("rad_vet","implicit_resid_fatal",1.0e2);
    recipe::Log("<rad_vet>/implicit_resid_fatal = 1e2 (recipe none: the old default)");
  }
  // implicit_predictor = step: the step predictor (ipred, pred_ok, pred_dt travel in the
  // restart file); none elsewhere
  {std::string pr = pin->GetOrAddString("rad_vet","implicit_predictor",
                                        fixcl ? "step" : "none");
  if (pr.compare("none") == 0) {
    impl_pred = false;
  } else if (pr.compare("step") == 0) {
    impl_pred = true;
  } else {
    ImplFatal("<rad_vet>/implicit_predictor = '" + pr
              + "' is not a choice (none | step)");
  }
  }
  // the positivity keys implicit_pos_gas and implicit_pos_floor: default true on a fresh
  // run; a restart whose file lacks them keeps the old default (off)
  {const bool rst = global_variable::restart_run;
  bool old_kept = false;
  auto pbool = [&](const char *key, bool &v) {
    if (pin->DoesParameterExist("rad_vet",key)) {
      v = pin->GetBoolean("rad_vet",key);
    } else if (rst) {
      old_kept = true;
    } else {
      v = pin->GetOrAddBoolean("rad_vet",key,true);
    }
  };
  pbool("implicit_pos_gas", impl_pos_gas);
  pbool("implicit_pos_floor", impl_pos_floor);
  if (old_kept && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> restart input lacks implicit_pos_gas / implicit_pos_floor: "
              << "keeping the old default (off) for those; name them to switch"
              << std::endl;
  }
  }
  if (pin->DoesParameterExist("rad_vet","implicit_pos_gas_frac")) {
    impl_pos_gas_frac = pin->GetReal("rad_vet","implicit_pos_gas_frac");
  }
  Kokkos::realloc(pos_cnt_d, M1_POS_N);
  Kokkos::deep_copy(pos_cnt_d, 0.0);
  marshak_q = pin->GetOrAddReal("rad_vet","marshak_q",0.5);
  // implicit_marshak_face = linear (the spherical-polar wedge only): the Marshak face
  // value extrapolated linearly; default linear on a fresh run with a fixed-tensor
  // closure, cell elsewhere
  if (sph_geom || pin->DoesParameterExist("rad_vet","implicit_marshak_face")) {
    const bool mfdef = !global_variable::restart_run &&
                       (eddington || vet_sc || tau_closure);
    const std::string smf = sph_geom ?
        pin->GetOrAddString("rad_vet","implicit_marshak_face",
                            mfdef ? "linear" : "cell") :
        pin->GetString("rad_vet","implicit_marshak_face");
    if (smf.compare("cell") == 0) {
      impl_mface_lin = false;
    } else if (smf.compare("linear") == 0) {
      impl_mface_lin = true;
    } else {
      ImplFatal("<rad_vet>/implicit_marshak_face = '" + smf
                + "' is not a choice (cell | linear)");
    }
    if (impl_mface_lin && !sph_geom) {
      ImplFatal("<rad_vet>/implicit_marshak_face = linear is implemented on the "
                "spherical-polar wedge only");
    }
  }
  std::string srn = pin->GetOrAddString("rad_vet","implicit_recon","dc");
  if (srn.compare("dc") == 0) {
    impl_recon = M1_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = M1_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_vet>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // implicit_enthalpy: the face enthalpy flux a_f E_f (plm on a fresh run with a
  // fixed-tensor closure; upwind elsewhere and on a restart whose file lacks the key)
  const bool enth_up = global_variable::restart_run || !fixcl;
  const std::string sen = pin->GetOrAddString("rad_vet","implicit_enthalpy",
                                              enth_up ? "upwind" : "plm");
  if (sen.compare("upwind") == 0) {
    impl_enth = M1_IENTH_UPWIND;
  } else if (sen.compare("central") == 0) {
    impl_enth = M1_IENTH_CENTRAL;
  } else if (sen.compare("plm") == 0) {
    impl_enth = M1_IENTH_PLM;
  } else if (sen.compare("plmc") == 0) {
    impl_enth = M1_IENTH_PLMC;
  } else if (sen.compare("lin") == 0) {
    impl_enth = M1_IENTH_LIN;
  } else {
    ImplFatal("<rad_vet>/implicit_enthalpy = '" + sen
              + "' is not a choice (upwind | central | plm)");
  }
  // implicit_vimp: the gas velocity implicit in the enthalpy flux, default on a fresh
  // multi-D hesdirk2 / be-path run with live gas coupling and bicgstab
  {auto *pmh = pmy_pack->pmesh;
  const bool gas_live = !pin->DoesParameterExist("time","evolution") ||
                        (pin->GetString("time","evolution").compare("dynamic") == 0);
  const bool vdef = ((pin->GetString("rad_vet","time_scheme").compare("hesdirk2") == 0) ||
                     t2_beonly) &&
                    pmh->multi_d && (impl_solver == M1_ISOLV_BICGSTAB) &&
                    (pmh->mb_indcs.ng >= 2) && fl_on && gas_live &&
                    coupling && gas_feedback && dbg_gas_force &&
                    !global_variable::restart_run;
  impl_vimp = pin->GetOrAddBoolean("rad_vet","implicit_vimp",vdef);}
  // <rad_vet>/doppler_face = donor | central (rad_vet.hpp dop_face; read only when named)
  if (pin->DoesParameterExist("rad_vet","doppler_face")) {
    const std::string df = pin->GetString("rad_vet","doppler_face");
    if (df == "central") {
      dop_face = 1;
    } else if (df != "donor") {
      ImplFatal("<rad_vet>/doppler_face must be donor | central");
    }
    if (dop_face == 1) {
      if (impl_vimp && global_variable::my_rank == 0) {
        std::cout << "<rad_vet> doppler_face = central: implicit_vimp overridden to "
                     "false "
                  << "(its Newton rows linearise the donor flux)" << std::endl;
      }
      impl_vimp = false;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> doppler_face = central: x1 enthalpy flux "
                  << "a E = v (1 + chi) E as the face mean 0.5 (a_L E_L + a_R E_R) "
                     "(Athena++ "
                  << "DO near-central face), no donor / plm correction" << std::endl;
      }
    }
  }
  {  // implicit_hold (default true): inert unless a pgen marks hold_m
    impl_hold = pin->GetOrAddBoolean("rad_vet","implicit_hold",true);
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
  // implicit_predictor_order = 2 (the step predictor from the last two increments) on a
  // fresh be / hesdirk2 run with a fixed-tensor closure; 1 elsewhere
  {
  bool ts_ok = true;
  if (pin->DoesParameterExist("rad_vet","time_scheme")) {
    const std::string tsn = pin->GetString("rad_vet","time_scheme");
    ts_ok = (tsn.compare("be") == 0) || (tsn.compare("hesdirk2") == 0);
  }
  const bool ldef = ts_ok && fixcl && !global_variable::restart_run;
  impl_pord = pin->GetOrAddInteger("rad_vet","implicit_predictor_order",
                                   (ldef && impl_pred) ? 2 : 1);
  if (impl_pord != 1 && impl_pord != 2) {
    ImplFatal("<rad_vet>/implicit_predictor_order must be 1 or 2");
  }
  }
  // cs_vertex_mono (default on the cubed sphere): two-point lateral faces within
  // cs_vertex_mono_n cells of the 8 cube vertices
  cs_vtx_mono = cs_geom;
  if (pin->DoesParameterExist("rad_vet","cs_vertex_mono")) {
    cs_vtx_mono = pin->GetBoolean("rad_vet","cs_vertex_mono") && cs_geom;
  }
  impl_pred_tcap = (recipe::GetMode(pin) == recipe::Mode::none) ? 2.0 : 1.1;
  if (pin->DoesParameterExist("rad_vet","implicit_predictor_tcap")) {
    impl_pred_tcap = pin->GetReal("rad_vet","implicit_predictor_tcap");
  }
  if (pin->DoesParameterExist("rad_vet","cs_vertex_mono_n")) {
    cs_vtx_n = pin->GetInteger("rad_vet","cs_vertex_mono_n");
    if (cs_vtx_n < 0) {
      ImplFatal("<rad_vet>/cs_vertex_mono_n must be >= 0");
    }
  }
  // the E floor's raise (implicit_pos_floor): implicit_pos_floor_cons takes its uncovered
  // part from the radiation of the same column, implicit_pos_gas_cons and
  // implicit_pos_floor_gas_keep (defaults true and 1 on a fresh run) bound what it takes
  // from the gas
  if (pin->DoesParameterExist("rad_vet","implicit_pos_floor_cons")) {
    impl_pos_flr_cons = pin->GetBoolean("rad_vet","implicit_pos_floor_cons");
  }
  if (impl_pos_flr_cons && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> implicit_pos_floor_cons: the E floor's uncovered raise is "
                 "taken "
              << "from the radiation of the same column (no energy created)" << std::endl;
  }
  const bool cnw = !global_variable::restart_run;
  if (cnw || pin->DoesParameterExist("rad_vet","implicit_pos_gas_cons")) {
    impl_pos_gas_cons = pin->GetOrAddBoolean("rad_vet","implicit_pos_gas_cons",true);
  }
  if (cnw || pin->DoesParameterExist("rad_vet","implicit_pos_floor_gas_keep")) {
    impl_pos_flr_gk = pin->GetOrAddReal("rad_vet","implicit_pos_floor_gas_keep",1.0);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> implicit_pos_floor_gas_keep=" << impl_pos_flr_gk
                << ": the E floor's raise leaves the gas at least that share of its eint"
                << std::endl;
    }
  }
  if (cs_vtx_mono && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> cs_vertex_mono: two-point lateral faces within " << cs_vtx_n
              << " cell(s) of the 8 cube vertices (no skew cross term, seam foot values, "
              << "D_r,lat = 0 within " << cs_vtx_n + 1 << ")" << std::endl;
  }
  // the gas source of the coupled solve: conservative2 (the cell-F slope csk1 travels in
  // the restart file) on a FRESH run with gas coupling (hydro, coupling, gas_feedback);
  // matched without it and on a restart whose file lacks the key (files echo it)
  {
    const bool gsdef = !global_variable::restart_run && fl_on && coupling &&
                       gas_feedback;
    const std::string gs = pin->GetOrAddString("rad_vet","implicit_gas_src",
                                               gsdef ? "conservative2" : "matched");
    if (gs.compare("matched") == 0) {
      impl_gas_src_cons = false;
    } else if (gs.compare("conservative") == 0 ||
               gs.compare("conservative_full") == 0) {
      impl_gas_src_cons = true;
    } else if (gs.compare("conservative2") == 0) {
      impl_gas_src_cons = true;
      cs2_vnp = pin->GetOrAddInteger("rad_vet","implicit_cs2_vimp_npass",
                                     global_variable::restart_run ? 0 : 3);
      if (cs2_vnp < 0) {
        ImplFatal("<rad_vet>/implicit_cs2_vimp_npass must be >= 0");
      }
    } else {
      ImplFatal("<rad_vet>/implicit_gas_src = '" + gs + "' is not a choice (matched | "
                "conservative | conservative2 | conservative_full)");
    }
  }

  auto &mbcs = pmy_pack->pmesh->mesh_bcs;
  ibc_x1min = ImplBCFromString(
      pin->GetOrAddString("rad_vet","implicit_bc_x1min","auto"),
      mbcs[static_cast<int>(BoundaryFace::inner_x1)]);
  ibc_x1max = ImplBCFromString(
      pin->GetOrAddString("rad_vet","implicit_bc_x1max","auto"),
      mbcs[static_cast<int>(BoundaryFace::outer_x1)]);
  // the transverse faces (read only when named; auto = the old rule)
  {
    auto tbc = [&](const char *key, const BoundaryFace f) -> int {
      const bool per = (mbcs[static_cast<int>(f)] == BoundaryFlag::periodic);
      std::string v = "auto";
      if (pin->DoesParameterExist("rad_vet", key)) {v = pin->GetString("rad_vet", key);}
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
          << std::endl << "<rad_vet>/" << key << " = '" << v << "': auto | reflect | "
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
  iflux_x1min = pin->GetOrAddReal("rad_vet","implicit_flux_x1min",0.0);
  // (read only when named; test only): implicit_flux_x1min_mod_amp = A,
  // implicit_flux_x1min_mod_period = P > 0: the inner face flux F0 (1 + A sin(2 pi t/P)),
  // F0 = implicit_flux_x1min, set before each hesdirk2 stage solve at its time
  // (Time2FormStage: t_n + dt), so the forcing of a time-order test is time-consistent
  if (pin->DoesParameterExist("rad_vet","implicit_flux_x1min_mod_period")) {
    fmod_per = pin->GetReal("rad_vet","implicit_flux_x1min_mod_period");
    fmod_amp = pin->GetReal("rad_vet","implicit_flux_x1min_mod_amp");
    fmod_f0 = iflux_x1min;
  }
  // (read only when named; test only): implicit_flux_x1min_ramp_t0 = t0,
  // implicit_flux_x1min_ramp_w = w > 0: the inner flux times 0.5 (1 + tanh((t - t0)/w)),
  // a smooth switch-on (no start-up discontinuity), also at the stage times; the value at
  // the start time is set here
  if (pin->DoesParameterExist("rad_vet","implicit_flux_x1min_ramp_w")) {
    fmod_rw = pin->GetReal("rad_vet","implicit_flux_x1min_ramp_w");
    fmod_rt0 = pin->GetReal("rad_vet","implicit_flux_x1min_ramp_t0");
    if (fmod_per <= 0.0) {fmod_f0 = iflux_x1min;}
    if (fmod_rw > 0.0) {
      const Real t = pmy_pack->pmesh->time;
      iflux_x1min = fmod_f0*0.5*(1.0 + std::tanh((t - fmod_rt0)/fmod_rw));
    }
  }
  // (read only when named; absent = bitwise): the inner face flux from
  // a (t, F or L) table (rad_vet.hpp ftab_on), set here at the start time and at the
  // stage times by the hooks above
  if (pin->DoesParameterExist("rad_vet","implicit_flux_x1min_table")) {
    const std::string fn = pin->GetString("rad_vet","implicit_flux_x1min_table");
    const std::string err = ftab.Read(fn, 1);
    if (!err.empty()) {ImplFatal("<rad_vet>/implicit_flux_x1min_table: " + err);}
    if (ibc_x1min != M1_IBC_FLUX) {
      ImplFatal("<rad_vet>/implicit_flux_x1min_table needs implicit_bc_x1min = flux");
    }
    ftab_t0 = pin->GetOrAddReal("rad_vet","implicit_flux_x1min_table_t0",0.0);
    const std::string un =
        pin->GetOrAddString("rad_vet","implicit_flux_x1min_table_units","flux");
    if (un == "flux") {
      ftab_scale = 1.0;
    } else if (un == "luminosity") {
      const Real r0 = pmy_pack->pmesh->mesh_size.x1min;
      ftab_scale = 1.0/(4.0*M_PI*r0*r0);
    } else {
      ImplFatal("<rad_vet>/implicit_flux_x1min_table_units must be flux|luminosity");
    }
    ftab_on = true;
    if (fmod_per <= 0.0 && fmod_rw <= 0.0) {fmod_f0 = iflux_x1min;}
    iflux_x1min = InnerFluxAt(pmy_pack->pmesh->time);
    if (global_variable::my_rank == 0) {
      std::cout << "rad_vet: implicit_flux_x1min_table " << fn << ": " << ftab.t.size()
                << " rows, t0 " << ftab_t0 << ", units " << un << ", F(t = "
                << pmy_pack->pmesh->time << ") = " << iflux_x1min << std::endl;
    }
  }
  iflux_x1max = pin->GetOrAddReal("rad_vet","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_vet","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_vet","implicit_ebath_x1max",0.0);
  impl_bc_advect = pin->GetOrAddBoolean("rad_vet","implicit_bc_advect",
                                        !global_variable::restart_run);
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  part_nblk = 1;
  trans_on = pmy_pack->pmesh->multi_d;
  trans_x3 = pmy_pack->pmesh->three_d;
  if (!trans_on && global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: transport=implicit on a 1-D mesh is the x1 column solve"
              << std::endl;
  }
  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/transport = implicit does not support SMR/AMR"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  nstage = 1;
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (impl_cfl > 0.0) || (!other_sets_dt);

  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(f0x1, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1, 0.0);
  bicg_on = trans_on && (impl_solver == M1_ISOLV_BICGSTAB);
  if (od_auto) {
    impl_offdiag = bicg_on ? M1_OD_OPERATOR : M1_OD_LAGGED;
  }
  if (impl_offdiag == M1_OD_OPERATOR && trans_on && !bicg_on) {
    ImplFatal("<rad_vet>/implicit_offdiag = operator needs implicit_solver = bicgstab");
  }
  od_now = impl_offdiag;
  if (impl_vimp && !bicg_on) {
    ImplFatal("<rad_vet>/implicit_vimp needs transport = implicit on a multi-D mesh with "
              "implicit_solver = bicgstab");
  }
  if (impl_vimp && indcs.ng < 2) {
    ImplFatal("<rad_vet>/implicit_vimp needs <mesh>/nghost >= 2");
  }
  // the cubed sphere's angular-momentum transport: the hesdirk2 spin slope (am_vec)
  // follows <hydro|mhd>/am_vector unless implicit_am_cs = false
  if (cs_geom) {
    auto named = [&](const char *b, const char *k) {
      return pin->DoesParameterExist(b, k) && pin->GetBoolean(b, k);
    };
    bool am_cs = named("hydro","am_transport") || named("hydro","am_vector") ||
                 named("mhd","am_transport") || named("mhd","am_vector");
    if (pin->DoesParameterExist("rad_vet","implicit_am_cs")) {
      am_cs = pin->GetBoolean("rad_vet","implicit_am_cs");
    } else {
      recipe::Log(std::string("<rad_vet>/implicit_am_cs = ") +
                  (am_cs ? "true" : "false") +
                  " (cubed sphere: follows <hydro|mhd>/am_transport)");
    }
    am_vec = am_cs && (named("hydro","am_vector") || named("mhd","am_vector"));
  }
  // the Zalesak limiter of the hesdirk2 face fluxes: zalesak_active (default true on a
  // fresh run) unless zalesak_local_rounds > 0 or zalesak_colsweep; zalesak_halo_skip
  if (pin->DoesParameterExist("rad_vet","zalesak_halo_skip")) {
    zl_hskip = pin->GetBoolean("rad_vet","zalesak_halo_skip");
  }
  if (pin->DoesParameterExist("rad_vet","zalesak_local_rounds")) {
    zl_local = pin->GetInteger("rad_vet","zalesak_local_rounds");
  }
  if (!global_variable::restart_run ||
      pin->DoesParameterExist("rad_vet","zalesak_active")) {
    const bool zcs = pin->DoesParameterExist("rad_vet","zalesak_colsweep") &&
                     pin->GetBoolean("rad_vet","zalesak_colsweep");
    zl_active = pin->GetOrAddBoolean("rad_vet","zalesak_active",true) &&
                (zl_local == 0) && !zcs;
  }
  if (pin->DoesParameterExist("rad_vet","zalesak_converge")) {
    zl_conv = pin->GetBoolean("rad_vet","zalesak_converge");
  }
  MemArenaInit(pin, nmb, M1_NIW, ncells3, ncells2, ncells1);
  if (impl_ctrelax > 0.0) {
    if (!trans_on || !clag_step) {
      ImplFatal("<rad_vet>/implicit_closure_thin_relax needs transport = implicit on a "
                "multi-D mesh and implicit_closure_lag = step");
    }
    Kokkos::realloc(ctr_mem, nmb, 4, ncells3, ncells2, ncells1);
  }
  Kokkos::deep_copy(iw, 0.0);
  if (impl_pred && !cln_not2h) {
    Kokkos::realloc(ipred, nmb, (impl_pord == 2) ? 5 : 3, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ipred, 0.0);
    pred_ok = false;
  }
  if (trans_on) {
    Kokkos::realloc(f0x2, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2, 0.0);
    if (trans_x3) {
      Kokkos::realloc(f0x3, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3, 0.0);
    }
    if (impl_tlim != M1_TLIM_NONE) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
      }
    }
    // implicit_hr_q = col (read only when named; default cell = bitwise): the outer
    // Marshak q of the half-range beam from the closure's own surface estimate (closure
    // vet_col + vet_gd) instead of the top cell's centre ratio H+_N/J_N
    bool hrq_col = false;
    if (pin->DoesParameterExist("rad_vet","implicit_hr_q")) {
      const std::string sq = pin->GetString("rad_vet","implicit_hr_q");
      if (sq.compare("col") == 0) {
        hrq_col = true;
      } else if (sq.compare("cell") != 0) {
        ImplFatal("<rad_vet>/implicit_hr_q = '" + sq + "' is not a choice (cell | col)");
      }
    }
    hr_q = impl_beam_hr && (vet_sc || (vgd_on && sph_geom)) &&
           (ibc_x1max == M1_IBC_MARSHAK) && !(hrq_col && vet_col);
    if (hr_q && vcol_q.extent_int(0) < nmb) {
      Kokkos::realloc(vcol_q, nmb, ncells3, ncells2);
      Kokkos::deep_copy(vcol_q, marshak_q);
    }
    // the transverse ghost exchanges (the clean solve makes thw / thq at their first use)
    Kokkos::realloc(thw_c, nmb, M1_NHALO_T, 1, 1, 1);
    pbval_th = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_th->exch_depth = M1HaloDepth(pmy_pack, pin);
    impl_hdepth = pbval_th->exch_depth;   // also the halo_mpi cell lists
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
    Kokkos::realloc(thq_c, nmb, M1_NHALO_Q, 1, 1, 1);
    pbval_tq = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_tq->exch_depth = M1HaloDepth(pmy_pack, pin);
    pbval_tq->SetVectorPairs(M1_NHALO_Q, {});   // E', G0, F1, KT: all scalars
    halo_shell = !(pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere ||
                   pmy_pack->pmesh->use_polar_boundary);
    Kokkos::realloc(krw, nmb, 1, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(krw, 0.0);
    Kokkos::realloc(krw_c, nmb, 1, 1, 1, 1);
    pbval_kr = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_kr->exch_depth = M1HaloDepth(pmy_pack, pin);
    pbval_kr->InitializeBuffers(1);
    pbval_kr->SetVectorPairs(1, {});   // one scalar Krylov vector
    if (impl_halo_direct) {ImplicitHaloDirectInit();}
  }
  ImplicitPartitionInit();

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: transport=implicit"
              << (trans_x3 ? " (3-D)" : (trans_on ? " (2-D)" : " (1-D)"))
              << " implicit_cfl=" << impl_cfl << " marshak_q=" << marshak_q
              << " implicit_enthalpy="
              << ((impl_enth == M1_IENTH_PLM) ? "plm"
                  : (impl_enth == M1_IENTH_PLMC) ? "plmc"
                  : (impl_enth == M1_IENTH_LIN) ? "lin"
                  : (impl_enth == M1_IENTH_UPWIND) ? "upwind" : "central")
              << " implicit_vimp=" << (impl_vimp ? "true" : "false") << std::endl;
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {VetInit(pin);}
  ImplBeamInit(pin);   // (rad_vet_vet.cpp)
  Time2Init(pin);
  MRInit(pin);         // implicit_mr_every (rad_vet_mr.cpp)
  if (time_scheme == M1_TIME_BE && mr_every <= 1) {
    Kokkos::realloc(u1, u1.extent(0), M1_NVAR, 1, 1, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitPartitionInit
//! \brief milestone 3b, LIMIT 4: build the topology of the x1 stacks and allocate the
//! gather/scatter and halo buffers.  A "stack" is the set of MeshBlocks that share the
//! same (x2,x3) footprint, ordered by their x1 logical location; its ROOT is the block
//! with the lowest one.  Uniform mesh only (the caller has already fatalled on SMR/AMR
//! and on a non-integer block count along x1), so the stack of a block is found by a
//! scan of the global LogicalLocation list.
//!
//! Nothing is allocated and nothing is communicated when part_nblk == 1: the 3a path is
//! then bitwise what it was.

void RadiationVet::ImplicitPartitionInit() {
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
        ImplFatal("<rad_vet>/implicit_partition = gather: the x1 stack of a MeshBlock is "
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
// <rad_vet>/implicit_halo_depth (default 2; 0 or >= ng = all ng layers, the old
// behaviour).  The implicit kernels read 2 ghost layers
// (the vimp enthalpy stencil +-2; poisoned-ghost gate M1_POISON_GHOST: NaN beyond depth
// 2 leaves the outputs identical, beyond depth 1 does not), so depth 2 sends the same
// values that are read with 2/3 of the bytes (RSG 2-node loop -2.5 %, bitwise).
// G5 (USER 11-04): also on the spherical-polar mesh (wedge and full sphere:
// the pole partner's rows are reversed within the strip, so a depth-2 strip fills the
// 2 ghost rows next to the pole; poisoned-ghost gate incl. the polar rows) when x1 is
// one MeshBlock (the gate poisons x2/x3 ghosts only; every sp production is x1-unsplit);
// not on Cartesian meshes (the old behaviour, all ng layers).
static int M1HaloDepth(MeshBlockPack *pp, ParameterInput *pin) {
  Mesh *pm = pp->pmesh;
  const bool sp = pm->use_spherical_polar && (pm->mesh_indcs.nx1 == pm->mb_indcs.nx1);
  if (!pm->use_cubed_sphere && !sp) {return 0;}
  int d = 2;
  if (pin != nullptr && pin->DoesParameterExist("rad_vet", "implicit_halo_depth")) {
    d = pin->GetInteger("rad_vet", "implicit_halo_depth");
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitHaloCopy
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

void RadiationVet::ImplicitHaloCopy(DvceArray5D<Real> &sc, int nq, int c0, bool topack) {
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
  // test switch (env M1_POISON_GHOST=1; gate only, NOT a run mode): on the way back,
  // every ghost cell more than ONE cell outside the active range in x2 or x3 is set to
  // NaN (if the outputs stay identical, no implicit kernel reads those layers of the
  // exchanged components)
  static int s_poison = -1;
  if (s_poison < 0) {
    s_poison = 0;
#if VET_DEBUG
    const char *f = std::getenv("M1_POISON_GHOST");
    s_poison = (f != nullptr) ? std::atoi(f) : 0;
#endif
  }
  const bool poi = (s_poison != 0) && !topack;
  const int pd = s_poison;   // M1_POISON_GHOST=D: the ghost layers beyond depth D
  const int pjs = indcs.js, pje = indcs.je, pks = indcs.ks, pke = indcs.ke;
  const bool p2 = (indcs.nx2 > 1), p3 = (indcs.nx3 > 1);
  if (impl_halo_direct || poi) {
    // implicit_halo_direct on a mesh whose neighbours are not all on this rank: the
    // same shell copy with ONE thread per cell (coalesced along i) instead of one per
    // row -- the row loop runs ~13k threads with strided rows on the 3-D box.
#if defined(KOKKOS_ENABLE_HIP)
    // (HIP): one 64-lane vector per (m,n,k,j) row, 4 rows per team, i.e.
    // the same cells, each copied by the same assignment, with the index division paid
    // per row instead of per cell (rocprofv3, HG small cs deck: m1_impl_hcpyf 0.40 s,
    // 7466 launches, ~2.6 TB/s).  Bitwise (plain copies).
    if (!poi) {
      const int nr2 = n2, nr3 = n3, nqq = nq;
      const int nrow = (nmb1 + 1)*nq*n3*n2;
      constexpr int kTr = 4;
      Kokkos::TeamPolicy<DevExeSpace, Kokkos::LaunchBounds<256, 2>>
          pol(DevExeSpace(), (nrow + kTr - 1)/kTr, kTr, 64);
      Kokkos::parallel_for("vet_impl_hcpyt", pol, KOKKOS_LAMBDA(TeamMember_t tm) {
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
    par_for("vet_impl_hcpyf", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
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
  par_for("vet_impl_hcpy", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1,
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
//! \fn void RadiationVet::ImplicitHaloExchange
//! \brief the common core: pack, exchange through the ordinary cell-centred boundary
//! machinery on the PERSISTENT scratch array whose width matches `nq`, unpack.  The
//! three widths (M1_NHALO_T, M1_NHALO_Q, 1) have one array and one boundary object each;
//! they are used strictly sequentially, so they may share the MPI tag space.

void RadiationVet::ImplicitHaloExchange(int nq, int c0) {
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
    par_for("vet_cs_nrenorm", DevExeSpace(), 0, pmy_pack->nmb_thispack-1, 0, n3-1,
            0, n2-1, 0, n1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
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
//! \fn void RadiationVet::ImplicitHaloDirectInit
//! \brief implicit_halo_direct: tabulate, for every MeshBlock of the pack and each of the
//! 26 directions (ox1,ox2,ox3), the LOCAL index of the same-level neighbour that fills
//! that ghost region, or -1 where there is none (a physical boundary: its ghost zones
//! are not filled by the ordinary exchange either).  The direct copy is used only when
//! EVERY rank finds every neighbour on its own rank at the same level, outside the
//! cubed-sphere and polar transforms (so all ranks take the same branch and no MPI
//! exchange is left half-posted); otherwise the ordinary exchange runs, as before.

void RadiationVet::ImplicitHaloDirectInit() {
  auto *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack;
  hd_src = DualArray2D<int>("vet_hd_src", nmb, 27);
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
//! \fn void RadiationVet::ImplicitHaloDirect
//! \brief implicit_halo_direct: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list) filled by ONE kernel that copies each ghost
//! cell from the active cell of the neighbour that owns it (all neighbours are on this
//! rank at the same level, ImplicitHaloDirectInit).  The same numbers the pack /
//! exchange / unpack chain delivers: 1 launch where there were 4, and no host work.

void RadiationVet::ImplicitHaloDirect(int nq, int c0) {
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
  par_for("vet_impl_hdir", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
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
//! \fn void RadiationVet::ImplicitKrylovHalo
//! \brief milestone 3b phase C: put ONE component of the work array into the scratch
//! array `krw`, exchange it with all six neighbours through the module's ordinary
//! cell-centred boundary machinery (which is what supplies periodic wrap, edge/corner
//! neighbours and MPI), and copy the ghost zones back.  PHYSICAL boundaries are not
//! filled: the row's coefficient towards such a neighbour is identically zero, so the
//! ghost value is multiplied by zero and never read in anger.

void RadiationVet::ImplicitKrylovHalo(int comp) {
  ImplicitHaloExchange(1, comp);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::TmrMark
//! \brief implicit_timers: fence, charge the time since the last mark to category c
//! (c < 0 only sets the mark)

void RadiationVet::TmrMark(int c) {
  if (!tmr_on) return;
  Kokkos::fence();
  const double t = tmr_t.seconds();
  if (c >= 0) {tmr_acc[c] += t - tmr_last;}
  tmr_last = t;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::EventTotals
//! \brief the limiter totals of this run (or restart) for the event log, all ranks
//! (collective): per-rank device counters, summed here.  0..6 (the old solver's
//! BiCGStab / Newton / positivity fallbacks) stay 0; 7 opac-Newton face terms dropped,
//! 8 E floor raises (write-back, cell-solves), 9 energy the floor CREATED (x volume),
//! 10 of which the solve undershoot, 11 ApplyClosureLimits E raises, 12
//! ApplyClosureLimits |F| clips (active cells, every stage)

void RadiationVet::EventTotals(Real *v) {
  for (int q = 0; q < 13; ++q) {v[q] = 0.0;}
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
//! \fn void RadiationVet::ImplicitReport
//! \brief the end-of-run lines of the implicit transport: the timers (implicit_timers),
//! the multi-rate report and the positivity bookings of the clean steps

void RadiationVet::ImplicitReport() {
  // the positivity counters and the hydro inversions' floor energy, all ranks
  // (collective: every rank reaches this line with the same pos keys)
  const bool pos = impl_pos_gas || impl_pos_floor;
  Real hfl = pmy_pack->pmesh->ecounter.flr_de_v;
  if (pos) {
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = hpc(q);}
#if MPI_PARALLEL_ENABLED
    Real g[M1_POS_N];
    MPI_Allreduce(pos_cnt, g, M1_POS_N, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = g[q];}
    Real gd = 0.0;
    MPI_Allreduce(&pos_un_drop, &gd, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    pos_un_drop = gd;
    MPI_Allreduce(MPI_IN_PLACE, &hfl, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  }
  if (global_variable::my_rank != 0) return;
  if (tmr_c0 > 0) {
    // implicit_timers: seconds per category (see TmrMark)
    std::cout << "<rad_vet> timers from cycle " << tmr_c0 << ": stages=" << tmr_cnt[7]
              << std::endl;
    const char *nm[11] = {"closure", "opacity", "solve_pre", "tensor", "pred",
                          "pass_setup", "krylov", "pass_post", "solve_end", "hyd_c2p",
                          "vet_bvals"};
    std::cout << "<rad_vet> timers(s):";
    for (int q = 0; q < 11; ++q) {std::cout << " " << nm[q] << "=" << tmr_acc[q];}
    std::cout << std::endl;
  }
  MRReport();
  if (pos) {
    std::cout << "<rad_vet> vet-positivity: implicit_pos_gas=" << impl_pos_gas << " (frac "
              << impl_pos_gas_frac << ") cell-solves=" << pos_cnt[M1_POS_GAS]
              << " energy moved rad->gas=" << pos_cnt[M1_POS_GAS_DE]
              << " | E floor raises cell-solves=" << pos_cnt[M1_POS_FLR]
              << " energy from gas=" << pos_cnt[M1_POS_FLR_DE]
              << " energy created=" << pos_cnt[M1_POS_FLR_UN] << " (code units x volume)"
              << " [in discarded attempts, rolled back: " << pos_un_drop << "]"
              << " taken from the column's radiation (implicit_pos_floor_cons)="
              << pos_cnt[M1_POS_FLR_MV]
              << " | hydro EOS floors energy created=" << hfl << std::endl;
  }
}
//----------------------------------------------------------------------------------------
//! implicit_hr_cmem_rst: fax1..3 to / from cell-sized restart slabs
//! (rad_vet.hpp).  Slab 2d: face d indices 0..n-1; slab 2d+1: face n in the first plane.

void RadiationVet::FaxRstPack(DvceArray5D<Real> out, int nmb) {
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

bool RadiationVet::FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl) {
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
    std::cout << "NOTE <rad_vet> restart: no time2_face_fct face fluxes in the file; the "
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

} // namespace radvet
