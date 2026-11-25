//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_init.cpp
//! \brief implicit-transport set-up (ImplicitInit), run-end report, face-flux restart
//! slabs, and the device memory arena shared by vet_gd and the clean work vectors.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro.hpp"
#include "utils/recipe.hpp"
#include "utils/scratch_view.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_clean.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {
//----------------------------------------------------------------------------------------
//! \fn ImplBCFromString
//! \brief x1 boundary type of the implicit solve from its key; auto: from the mesh flag

int ImplBCFromString(const std::string &s, const BoundaryFlag mbc) {
  if (s.compare("auto") == 0) {
    if (mbc == BoundaryFlag::periodic) return VET_IBC_PERIODIC;
    if (mbc == BoundaryFlag::reflect) return VET_IBC_REFLECT;
    return VET_IBC_MARSHAK;
  }
  if (s.compare("marshak") == 0) return VET_IBC_MARSHAK;
  if (s.compare("flux") == 0) return VET_IBC_FLUX;
  if (s.compare("reflect") == 0) return VET_IBC_REFLECT;
  if (s.compare("periodic") == 0) return VET_IBC_PERIODIC;
  if (s.compare("efix") == 0) return VET_IBC_EFIX;
  if (s.compare("outflow") == 0) return VET_IBC_OUTFLOW;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    << std::endl << "<rad_vet>/implicit_bc_x1min|max = '" << s << "' is not a valid "
    << "choice (auto | marshak | flux | reflect | periodic | efix | outflow (scheme = "
    << "clean))" << std::endl;
  std::exit(EXIT_FAILURE);
  return VET_IBC_MARSHAK;
}

//----------------------------------------------------------------------------------------
//! \fn ImplFatal
//! \brief fatal-error exit with a message

void ImplFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitInit
//! \brief reads and checks the implicit-transport keys, allocates the face and work
//! arrays.  Called at the end of the constructor.

void RadiationVet::ImplicitInit(ParameterInput *pin) {
  // time_scheme (read by Time2Init): hesdirk2 on a fresh run where accepted (no polar
  // boundary, not tau unless vet_col, integrator rk2), else be
  {auto *pmh = pmy_pack->pmesh;
  const std::string integ = pin->DoesParameterExist("time","integrator") ?
                            pin->GetString("time","integrator") : "rk2";
  const bool h2ok = !pmh->use_polar_boundary && (!tau_closure || vet_col) &&
                    (integ.compare("rk2") == 0);
  (void) pin->GetOrAddString("rad_vet","time_scheme",
                             (h2ok && !global_variable::restart_run) ? "hesdirk2" : "be");
  // be runs the hesdirk2 backward-Euler path (t2_beonly) where hesdirk2 is accepted;
  // be_legacy: plain BE step (default on restart); time2_be_only: BE path under hesdirk2
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
  // implicit_solver: bicgstab under the recipe defaults, else line_jacobi
  if (!pin->DoesParameterExist("rad_vet","implicit_solver")) {
    const bool nd = recipe::NewDefaults(pin);
    recipe::Log(std::string("<rad_vet>/implicit_solver = ") +
                (nd ? "bicgstab (default)" : "line_jacobi (the old default)"));
    pin->GetOrAddString("rad_vet","implicit_solver", nd ? "bicgstab" : "line_jacobi");
  }
  {std::string sv = pin->GetOrAddString("rad_vet","implicit_solver","line_jacobi");
  if (sv.compare("line_jacobi") == 0) {
    impl_solver = VET_ISOLV_LINE_JACOBI;
  } else if (sv.compare("bicgstab") == 0) {
    impl_solver = VET_ISOLV_BICGSTAB;
  } else {
    ImplFatal("<rad_vet>/implicit_solver = '" + sv
              + "' is not a choice (line_jacobi | bicgstab)");
  }
  }
  // implicit_offdiag: auto = operator with multi-D bicgstab, else lagged
  std::string sod = pin->GetOrAddString("rad_vet","implicit_offdiag","auto");
  bool od_auto = (sod.compare("auto") == 0);
  if (od_auto || sod.compare("lagged") == 0) {
    impl_offdiag = VET_OD_LAGGED;
  } else if (sod.compare("operator") == 0) {
    impl_offdiag = VET_OD_OPERATOR;
  } else if (sod.compare("none") == 0) {
    impl_offdiag = VET_OD_NONE;
  } else {
    ImplFatal("<rad_vet>/implicit_offdiag = '" + sod
              + "' is not a choice (auto | lagged | operator | none)");
  }
  // implicit_closure_lag: step (multi-D fresh run) | pass.  With step, chi(f) closures
  // relax thin cells by implicit_closure_thin_relax (memory ctr_mem is in the restart);
  // fixed-tensor closures never relax
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
  // implicit_trans_limit: transverse realizability limiter
  {std::string st = pin->GetOrAddString("rad_vet","implicit_trans_limit","none");
  if (st.compare("none") == 0) {
    impl_tlim = VET_TLIM_NONE;
  } else if (st.compare("lp") == 0) {
    impl_tlim = VET_TLIM_LP;
  } else if (st.compare("lpr") == 0) {
    impl_tlim = VET_TLIM_LPR;
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
  // implicit_resid_fatal: stage residual that stops the run (0 = off)
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
  }
  // implicit_predictor: step | none (the predictor state is in the restart)
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
  // implicit_pos_gas, implicit_pos_floor (default true)
  impl_pos_gas = pin->GetOrAddBoolean("rad_vet","implicit_pos_gas",
      recipe::Fix(pin, "rad_vet", "implicit_pos_gas", "true", "false"));
  impl_pos_floor = pin->GetOrAddBoolean("rad_vet","implicit_pos_floor",
      recipe::Fix(pin, "rad_vet", "implicit_pos_floor", "true", "false"));
  if (pin->DoesParameterExist("rad_vet","implicit_pos_gas_frac")) {
    impl_pos_gas_frac = pin->GetReal("rad_vet","implicit_pos_gas_frac");
  }
  Kokkos::realloc(pos_cnt_d, VET_POS_N);
  Kokkos::deep_copy(pos_cnt_d, 0.0);
  marshak_q = pin->GetOrAddReal("rad_vet","marshak_q",0.5);
  // implicit_marshak_face: cell | linear (Marshak face value; spherical wedge only)
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
    impl_recon = VET_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = VET_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_vet>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // implicit_enthalpy: reconstruction of the face enthalpy flux a_f E_f
  const bool enth_up = global_variable::restart_run || !fixcl;
  const std::string sen = pin->GetOrAddString("rad_vet","implicit_enthalpy",
                                              enth_up ? "upwind" : "plm");
  if (sen.compare("upwind") == 0) {
    impl_enth = VET_IENTH_UPWIND;
  } else if (sen.compare("central") == 0) {
    impl_enth = VET_IENTH_CENTRAL;
  } else if (sen.compare("plm") == 0) {
    impl_enth = VET_IENTH_PLM;
  } else if (sen.compare("plmc") == 0) {
    impl_enth = VET_IENTH_PLMC;
  } else if (sen.compare("lin") == 0) {
    impl_enth = VET_IENTH_LIN;
  } else {
    ImplFatal("<rad_vet>/implicit_enthalpy = '" + sen
              + "' is not a choice (upwind | central | plm)");
  }
  // implicit_vimp: gas velocity implicit in the enthalpy flux
  {auto *pmh = pmy_pack->pmesh;
  const bool gas_live = !pin->DoesParameterExist("time","evolution") ||
                        (pin->GetString("time","evolution").compare("dynamic") == 0);
  const bool vdef = ((pin->GetString("rad_vet","time_scheme").compare("hesdirk2") == 0) ||
                     t2_beonly) &&
                    pmh->multi_d && (impl_solver == VET_ISOLV_BICGSTAB) &&
                    (pmh->mb_indcs.ng >= 2) && fl_on && gas_live &&
                    coupling && gas_feedback && dbg_gas_force &&
                    !global_variable::restart_run;
  impl_vimp = pin->GetOrAddBoolean("rad_vet","implicit_vimp",vdef);}
  // doppler_face: donor | central x1 enthalpy flux (read only when named)
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
  {  // implicit_hold: inert unless a pgen marks hold_m
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
  // implicit_predictor_order: 1 | 2 (2 = from the last two increments)
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
  // cs_vertex_mono: two-point lateral faces near the 8 cube vertices (cs_vertex_mono_n)
  cs_vtx_mono = pin->GetOrAddBoolean("rad_vet","cs_vertex_mono",true) && cs_geom;
  impl_pred_tcap = 1.1;
  impl_pred_tcap = pin->GetOrAddReal("rad_vet","implicit_predictor_tcap",impl_pred_tcap);
  if (pin->DoesParameterExist("rad_vet","cs_vertex_mono_n")) {
    cs_vtx_n = pin->GetInteger("rad_vet","cs_vertex_mono_n");
    if (cs_vtx_n < 0) {
      ImplFatal("<rad_vet>/cs_vertex_mono_n must be >= 0");
    }
  }
  // E floor raise: implicit_pos_floor_cons takes the uncovered part from the column's
  // radiation; implicit_pos_gas_cons / implicit_pos_floor_gas_keep bound the gas share
  impl_pos_flr_cons =
      pin->GetOrAddBoolean("rad_vet","implicit_pos_floor_cons",impl_pos_flr_cons);
  if (impl_pos_flr_cons && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> implicit_pos_floor_cons: the E floor's uncovered raise is "
                 "taken "
              << "from the radiation of the same column (no energy created)" << std::endl;
  }
  impl_pos_gas_cons = pin->GetOrAddBoolean("rad_vet","implicit_pos_gas_cons",
      recipe::Fix(pin, "rad_vet", "implicit_pos_gas_cons", "true", "false"));
  impl_pos_flr_gk = pin->GetOrAddReal("rad_vet","implicit_pos_floor_gas_keep",
      recipe::Fix(pin, "rad_vet", "implicit_pos_floor_gas_keep", "1", "-1") ? 1.0 : -1.0);
  if (global_variable::my_rank == 0 && impl_pos_flr_gk > 0.0) {
    std::cout << "<rad_vet> implicit_pos_floor_gas_keep=" << impl_pos_flr_gk
              << ": the E floor's raise leaves the gas at least that share of its eint"
              << std::endl;
  }
  if (cs_vtx_mono && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> cs_vertex_mono: two-point lateral faces within " << cs_vtx_n
              << " cell(s) of the 8 cube vertices (no skew cross term, seam foot values, "
              << "D_r,lat = 0 within " << cs_vtx_n + 1 << ")" << std::endl;
  }
  // implicit_gas_src: gas source of the coupled solve (conservative2: csk1 in restart)
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
  // transverse boundaries (read only when named; auto: periodic or reflect)
  {
    auto tbc = [&](const char *key, const BoundaryFace f) -> int {
      const bool per = (mbcs[static_cast<int>(f)] == BoundaryFlag::periodic);
      std::string v = "auto";
      if (pin->DoesParameterExist("rad_vet", key)) {v = pin->GetString("rad_vet", key);}
      if (v.compare("auto") == 0) {return per ? VET_IBC_PERIODIC : VET_IBC_REFLECT;}
      int r = -1;
      if (v.compare("reflect") == 0) {r = VET_IBC_REFLECT;}
      if (v.compare("marshak") == 0) {r = VET_IBC_MARSHAK;}
      // zero-gradient outgoing flux, x2 faces only
      if (v.compare("outflow") == 0 && (f == BoundaryFace::inner_x2 ||
                                        f == BoundaryFace::outer_x2)) {
        r = VET_IBC_OUTFLOW;
      }
      if (v.compare("flux") == 0 && (f == BoundaryFace::inner_x2 ||
                                     f == BoundaryFace::outer_x2)) {r = VET_IBC_X2FLUX;}
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
  iflux_x1max = pin->GetOrAddReal("rad_vet","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_vet","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_vet","implicit_ebath_x1max",0.0);
  if ((ibc_x1min == VET_IBC_PERIODIC) != (ibc_x1max == VET_IBC_PERIODIC)) {
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
  bicg_on = trans_on && (impl_solver == VET_ISOLV_BICGSTAB);
  if (od_auto) {
    impl_offdiag = bicg_on ? VET_OD_OPERATOR : VET_OD_LAGGED;
  }
  if (impl_offdiag == VET_OD_OPERATOR && trans_on && !bicg_on) {
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
  // cubed sphere: the spin slope am_vec follows <hydro|mhd>/am_vector unless
  // implicit_am_cs = false
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
  // Zalesak limiter of the hesdirk2 face fluxes: zalesak_halo_skip, zalesak_active
  if (pin->DoesParameterExist("rad_vet","zalesak_halo_skip")) {
    zl_hskip = pin->GetBoolean("rad_vet","zalesak_halo_skip");
  }
  zl_active = pin->GetOrAddBoolean("rad_vet","zalesak_active",
      recipe::Fix(pin, "rad_vet", "zalesak_active", "true", "false"));
  MemArenaInit(pin, nmb, VET_NIW, ncells3, ncells2, ncells1);
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
    if (impl_tlim != VET_TLIM_NONE) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
      }
    }
    // implicit_hr_q: cell | col (outer Marshak q of the half-range beam: top-cell
    // H+_N/J_N, or the vet_col surface estimate)
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
           (ibc_x1max == VET_IBC_MARSHAK) && !(hrq_col && vet_col);
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
              << ((impl_enth == VET_IENTH_PLM) ? "plm"
                  : (impl_enth == VET_IENTH_PLMC) ? "plmc"
                  : (impl_enth == VET_IENTH_LIN) ? "lin"
                  : (impl_enth == VET_IENTH_UPWIND) ? "upwind" : "central")
              << " implicit_vimp=" << (impl_vimp ? "true" : "false") << std::endl;
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {VetInit(pin);}
  ImplBeamInit(pin);
  Time2Init(pin);
  if (time_scheme == VET_TIME_BE) {
    Kokkos::realloc(u1, u1.extent(0), VET_NVAR, 1, 1, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::TmrMark
//! \brief implicit_timers: fence, charge the time since the last mark to category c
//! (c < 0: only set the mark)

void RadiationVet::TmrMark(int c) {
  if (!tmr_on) return;
  Kokkos::fence();
  const double t = tmr_t.seconds();
  if (c >= 0) {tmr_acc[c] += t - tmr_last;}
  tmr_last = t;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::EventTotals
//! \brief limiter totals for the event log, summed over ranks (collective).
//! Slots: 0..6 unused (0); 7 opac-Newton face terms dropped; 8 E floor raises;
//! 9 floor energy created (x volume); 10 of which solve undershoot; 11-12 closure clips

void RadiationVet::EventTotals(Real *v) {
  for (int q = 0; q < 13; ++q) {v[q] = 0.0;}
  if (opn_nskip_d.extent_int(0) > 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    v[7] = h(0);
  }
  if (pos_cnt_d.extent_int(0) >= VET_POS_N) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    v[8] = h(VET_POS_FLR);
    v[9] = h(VET_POS_FLR_UN);
    v[10] = h(VET_POS_FLRS);
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
//! \brief end-of-run report: timers and positivity bookings

void RadiationVet::ImplicitReport() {
  // positivity counters and hydro floor energy (collective: same pos keys on all ranks)
  const bool pos = impl_pos_gas || impl_pos_floor;
  Real hfl = pmy_pack->pmesh->ecounter.flr_de_v;
  if (pos) {
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    for (int q = 0; q < VET_POS_N; ++q) {pos_cnt[q] = hpc(q);}
#if MPI_PARALLEL_ENABLED
    Real g[VET_POS_N];
    MPI_Allreduce(pos_cnt, g, VET_POS_N, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < VET_POS_N; ++q) {pos_cnt[q] = g[q];}
    Real gd = 0.0;
    MPI_Allreduce(&pos_un_drop, &gd, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    pos_un_drop = gd;
    MPI_Allreduce(MPI_IN_PLACE, &hfl, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  }
  if (global_variable::my_rank != 0) return;
  if (tmr_c0 > 0) {
    // seconds per category (see TmrMark)
    std::cout << "<rad_vet> timers from cycle " << tmr_c0 << ": stages=" << tmr_cnt[7]
              << std::endl;
    const char *nm[11] = {"closure", "opacity", "solve_pre", "tensor", "pred",
                          "pass_setup", "krylov", "pass_post", "solve_end", "hyd_c2p",
                          "vet_bvals"};
    std::cout << "<rad_vet> timers(s):";
    for (int q = 0; q < 11; ++q) {std::cout << " " << nm[q] << "=" << tmr_acc[q];}
    std::cout << std::endl;
  }
  if (pos) {
    std::cout << "<rad_vet> vet-positivity: implicit_pos_gas=" << impl_pos_gas
              << " (frac " << impl_pos_gas_frac << ") cell-solves=" <<
                  pos_cnt[VET_POS_GAS]
              << " energy moved rad->gas=" << pos_cnt[VET_POS_GAS_DE]
              << " | E floor raises cell-solves=" << pos_cnt[VET_POS_FLR]
              << " energy from gas=" << pos_cnt[VET_POS_FLR_DE]
              << " energy created=" << pos_cnt[VET_POS_FLR_UN] << " (code units x volume)"
              << " [in discarded attempts, rolled back: " << pos_un_drop << "]"
              << " taken from the column's radiation (implicit_pos_floor_cons)="
              << pos_cnt[VET_POS_FLR_MV]
              << " | hydro EOS floors energy created=" << hfl << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::FaxRstPack / FaxRstUnpack
//! \brief time2_face_fct: t2pf1 to / from cell-sized restart slabs (with ghosts), one
//! per direction; files with two slabs per direction: the first one is read

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
  // file without them: the first stage takes the stage-start face F
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

//========================================================================================
// device memory arena (vet_gd intensities, clean work vectors over them)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaInit
//! \brief enables the arena with vet_gd; allocates iw (not in the arena: some of its
//! channels live across a VetGdBuild)

void RadiationVet::MemArenaInit(ParameterInput *pin, const int nmb, const int niw,
                                const int n3, const int n2, const int n1) {
  marena_on = vgd_on;
  Kokkos::realloc(iw, nmb, niw, n3, n2, n1);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::MemArenaVgd
//! \brief binds a.d (n Reals) into the arena after the previous intensity array;
//! false: not bound, the caller allocates a.d

bool RadiationVet::MemArenaVgd(VgdRag &a, const int64_t n) {
  if (!marena_on || marena_sealed) return false;
  const int64_t off = marena_sw;
  const int64_t need = off + n;
  MemArenaGrow(need, false);   // (the earlier array holds zeros, as a new arena)
  a.d = DvceArray1D<Real>(marena.data() + off, n);
  marena_sw = need + ((64 - (need % 64)) % 64);   // 64-Real alignment
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaGrow
//! \brief grows the arena to need Reals and rebinds the views in it (same offsets);
//! keep: copy the old contents, else free the old allocation first, nothing copied

void RadiationVet::MemArenaGrow(const int64_t need, const bool keep) {
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  if (need <= have) return;
  Kokkos::fence();
  const Real *o0 = marena.data(), *o1 = marena.data() + have;
  if (!keep) marena = DvceArray1D<Real>();
  DvceArray1D<Real> nb("vet_arena", need);
  if (have > 0 && keep) {
    Kokkos::deep_copy(Kokkos::subview(nb, std::make_pair(int64_t(0), have)), marena);
  }
  VgdRag *bound[2] = {&vgd_i, &vgd_i_xb};
  for (VgdRag *b : bound) {
    if (have > 0 && b->d.data() != nullptr && b->d.data() >= o0 && b->d.data() < o1) {
      b->d = DvceArray1D<Real>(nb.data() + (b->d.data() - o0), b->d.extent(0));
    }
  }
  for (DvceArray5D<Real> *w : {&cln_w.o, &cln_w.s}) {   // bound by CleanWoBind
    if (have > 0 && w->data() != nullptr && w->data() >= o0 && w->data() < o1) {
      *w = DvceArray5D<Real>(nb.data() + (w->data() - o0), w->extent(0), w->extent(1),
                             w->extent(2), w->extent(3), w->extent(4));
    }
  }
  marena = nb;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaSeal
//! \brief called by every DtRetry::Save before it copies; the first call after the
//! intensities exist moves cln_w.o into the arena (copied, rebound), then seals it

void RadiationVet::MemArenaSeal() {
  // each checkpoint takes the handles anew, so a retry after the move restores them
  if (marena_sealed || !marena_on || marena.data() == nullptr) return;
  marena_sealed = true;
  DvceArray5D<Real> *cand[] = {&cln_w.o, &ost, &fmm_t, &odc, &d3w, &aa_xc, &aa_fc,
                               &aa_xp, &aa_fp};
  const char *nm[] = {"cln_wo", "ost", "fmm_t", "odc", "d3w", "aa_xc", "aa_fc", "aa_xp",
                      "aa_fp"};
  // grow to hold cln_w.o where it outgrows the intensities
  MemArenaGrow(static_cast<int64_t>(cln_w.o.span()));
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  int64_t off = 0;
  std::string moved;
  Kokkos::fence();
  for (int q = 0; q < static_cast<int>(sizeof(cand)/sizeof(cand[0])); ++q) {
    DvceArray5D<Real> &a = *cand[q];
    const int64_t n = static_cast<int64_t>(a.span());
    // only cln_w.o (q = 0) moves
    if (q > 0 || a.data() == nullptr || n == 0 ||
        !a.span_is_contiguous() || off + n > have || a.data() == marena.data()) {
      // already in the arena (CleanWoBind)
      if (a.data() == marena.data() && n > 0) off += n + ((64 - (n % 64)) % 64);
      continue;
    }
    DvceArray5D<Real> b(marena.data() + off, a.extent(0), a.extent(1), a.extent(2),
                        a.extent(3), a.extent(4));
    Kokkos::deep_copy(b, a);
    a = b;
    off += n + ((64 - (n % 64)) % 64);
    moved += std::string(" ") + nm[q];
  }
  MemArenaBind();
  Kokkos::fence();
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> mem_arena: intensities " << marena_sw*sizeof(Real)/1.0e9
              << " GB in one allocation; solve scratch moved over them (" << off*
              sizeof(Real)/1.0e9 << " GB on rank 0):" << moved << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaBind
//! \brief once cln_w.o lies in the arena: moves the scratch arrays that are dead across a
//! closure build, written before they are read in each use, over the intensities too:
//! after cln_w the solver's (matrix rows, lateral sums, multigrid levels, write-back
//! scratch), and over both, from the start, the hydro's flux-stage arrays (dead across
//! the radiation stage: x2, x3 fluxes, FOFC's trial state, the well-balanced cache, the
//! seam, block-end and dfloor_flux_limit records).  The history reads uflx.x1f: kept.
//! Behind both: DualEnergy's pre-conversion state, used within Hydro::ConToPrim, and K2
//! (t2k2, t2pf2)

void RadiationVet::MemArenaBind() {
  const Real *a0 = marena.data();
  const int64_t have = static_cast<int64_t>(marena.extent(0));
  if (marena_bound || a0 == nullptr || cln_w.o.data() != a0) return;
  marena_bound = true;
  Kokkos::fence();
  auto align = [](const int64_t n) {return n + ((64 - (n % 64)) % 64);};
  int64_t off = 0;
  for (auto *w : {&cln_w.o, &cln_w.s}) {
    const int64_t o = w->data() - a0;
    if (w->data() != nullptr && o >= 0 && o < have) {
      off = std::max(off, align(o + static_cast<int64_t>(w->span())));
    }
  }
  int64_t moved = 0;
  auto bind = [&](DvceArray5D<Real> &a) {
    const int64_t n = static_cast<int64_t>(a.span());
    if (a.data() == nullptr || n <= 1 || off + n > have) return;
    a = DvceArray5D<Real>(marena.data() + off, a.extent(0), a.extent(1), a.extent(2),
                          a.extent(3), a.extent(4));
    off = align(off + n);
    moved += n;
  };
  // (cubed sphere: cln_st holds the build's face-layer chunks, CleanGdFaceSums)
  if (!pmy_pack->pmesh->use_cubed_sphere) bind(cln_st);
  bind(cln_yl);
  bind(cln_vf);
  for (auto &a : cln_mg) bind(a);
  // (the closure built once per stage: no pass keeps the faces' coefficients over one)
  if (!cln_closit && !(cln_fwbst && t2_tab == 0)) {
    for (CleanF *f : {&cln_fc1, &cln_fc2, &cln_fc3}) bind(f->o);
  }
  // (the write-back's: their first use makes them, here already over the arena)
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int n1 = id.nx1 + 2*id.ng, n2 = (id.nx2 > 1) ? id.nx2 + 2*id.ng : 1;
  const int n3 = (id.nx3 > 1) ? id.nx3 + 2*id.ng : 1, nmb = pmy_pack->nmb_thispack;
  for (const auto &s : {std::make_pair("cln_wbfp", 4), std::make_pair("cln_xf", 3)}) {
    const int64_t n = static_cast<int64_t>(nmb)*s.second*n3*n2*n1;
    if (off + n > have) continue;
    ScratchBind(s.first, DvceArray5D<Real>(marena.data() + off, nmb, s.second, n3, n2,
                                           n1));
    off = align(off + n);
    moved += n;
  }
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph != nullptr) {
    const int64_t os = off;
    off = 0;
    bind(ph->wbq0);
    if (ph->rsolver_method != Hydro_RSolver::advect) {   // (it leaves IEN's fluxes zero)
      bind(ph->uflx.x2f);
      bind(ph->uflx.x3f);
    }
    bind(ph->utest);
    auto record = [&](DvceFaceFld5D<Real> *r) {
      if (r == nullptr) return;
      bind(r->x2f);
      bind(r->x3f);
    };
    record(ph->pseamfs);
    if (ph->pmlfs != ph->pseamfs) record(ph->pmlfs);
    if (ph->use_dfl) ph->DfloorFluxInit();
    record(ph->pdfl_fs);
    marena_hyd = off;
    // (behind both: the conserved state before each conversion to primitives)
    off = std::max(off, os);
    bind(ph->dual.pre);
  }
  // K2 and its face fluxes: from the write-back of stage 1 to the end of the step; with
  // extrap1 no closure is built between them (sdirk2: the stage-2 solve writes K1)
  if (cln_extrap == 3 && t2_tab == 2) {
    bind(t2k2);
    bind(t2pf2);
    marena_k2 = true;
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> mem_arena: " << moved*sizeof(Real)/1.0e9
              << " GB of solver and hydro scratch over the intensities" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::MemArenaHydro
//! \brief the hydro's arrays over the intensities, zero at each flux stage as their own
//! allocations are where nothing writes them (FOFC's trial reads faces beyond the sweeps)

void RadiationVet::MemArenaHydro() {
  if (marena_hyd == 0) return;
  auto h = Kokkos::subview(marena, std::make_pair(int64_t(0), marena_hyd));
  Kokkos::deep_copy(DevExeSpace(), h, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanWoBind
//! \brief allocates cln_w.o at the first solve: over the intensities in the arena when
//! it exists (dead outside a build), else its own allocation.  cln_w.s, if CleanInit
//! left it unmade, follows cln_w.o (after it in the arena, or its own allocation)

void RadiationVet::CleanWoBind() {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int n1 = id.nx1 + 2*id.ng, n2 = (id.nx2 > 1) ? id.nx2 + 2*id.ng : 1;
  const int n3 = (id.nx3 > 1) ? id.nx3 + 2*id.ng : 1, nmb = pmy_pack->nmb_thispack;
  if (cln_w.o.data() == nullptr) {
    if (marena.data() == nullptr) {
      Kokkos::realloc(cln_w.o, nmb, cln_wo_nc, n3, n2, n1);
    } else {
      MemArenaGrow(static_cast<int64_t>(nmb)*cln_wo_nc*n3*n2*n1, false);
      cln_w.o = DvceArray5D<Real>(marena.data(), nmb, cln_wo_nc, n3, n2, n1);
    }
  }
  if (cln_w.s.data() != nullptr) {
    MemArenaBind();
    return;
  }
  if (cln_w.o.data() != marena.data()) {
    Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
    return;
  }
  const int64_t no = static_cast<int64_t>(cln_w.o.span());
  const int64_t off = no + ((64 - (no % 64)) % 64);   // 64-Real alignment
  MemArenaGrow(off + static_cast<int64_t>(nmb)*CLW_NS*n3*n2*n1, false);
  cln_w.s = DvceArray5D<Real>(marena.data() + off, nmb, CLW_NS, n3, n2, n1);
  MemArenaBind();
}

} // namespace radvet
