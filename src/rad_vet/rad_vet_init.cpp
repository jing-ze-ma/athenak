//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_init.cpp
//! \brief set-up of the implicit radiation (<rad_vet>/transport = implicit): the keys the
//! clean steps read (ImplicitInit), the run-end report (ImplicitReport, EventTotals) and
//! the face-flux restart slabs (FaxRst*).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "rad_vet/rad_vet_t2reset.hpp"
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
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
  const bool fixcl = eddington || vet_sc || tau_closure;
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
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
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
  // the Zalesak limiter of the hesdirk2 face fluxes: zalesak_halo_skip, zalesak_active
  // (default true on a fresh run)
  if (pin->DoesParameterExist("rad_vet","zalesak_halo_skip")) {
    zl_hskip = pin->GetBoolean("rad_vet","zalesak_halo_skip");
  }
  if (!global_variable::restart_run ||
      pin->DoesParameterExist("rad_vet","zalesak_active")) {
    zl_active = pin->GetOrAddBoolean("rad_vet","zalesak_active",true);
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
  }

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
//! time2_face_fct: t2pf1 to / from cell-sized restart slabs (rad_vet.hpp FaxRstNslab),
//! one per direction (files before: two, the second 0); the slabs include the ghosts, so
//! every face a solve set

void RadiationVet::FaxRstPack(DvceArray5D<Real> out, int nmb) {
  const int n3 = out.extent_int(2), n2 = out.extent_int(3), n1 = out.extent_int(4);
  const int nd = FaxRstNslab();
  auto pf = t2pf1;
  par_for("pf_rst_pack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int d = 0; d < nd; ++d) out(m,d,k,j,i) = pf(m,d,k,j,i);
  });
}

bool RadiationVet::FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl) {
  const int n3 = in.extent_int(2), n2 = in.extent_int(3), n1 = in.extent_int(4);
  // (a file without them: the first stage takes the stage-start face F, Time2FaceFct)
  const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
  if (t2_pfst && t2pf1.extent_int(0) >= nmb && (nsl == nd || nsl == 2*nd)) {
    auto pf = t2pf1;
    const int st = (nsl == nd) ? 1 : 2;
    par_for("pf_rst_unpack", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < nd; ++d) pf(m,d,k,j,i) = in(m,st*d,k,j,i);
    });
    t2pf1_ok = true;
    return true;
  }
  if (t2_pfst && global_variable::my_rank == 0) {
    std::cout << "NOTE <rad_vet> restart: no time2_face_fct face fluxes in the file; the "
              << "first stage takes the stage-start face F" << std::endl;
  }
  return false;
}

} // namespace radvet
