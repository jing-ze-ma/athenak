//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_params_implicit.cpp
//! \brief keys of the implicit transport (ImplicitInit) and of its IMEX stages
//! (Time2Init).
//!
//!  called from: rad_vet.cpp
//!  calls into: rad_vet_init.cpp, rays/sc_cart_setup.cpp

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/imex_sdirk2.hpp"
#include "utils/recipe.hpp"

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

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplicitInit
//! \brief reads and checks the implicit-transport keys, allocates the face and work
//! arrays.  Called at the end of the constructor.

void RadiationVet::ImplicitInit(ParameterInput *pin) {
  // time_scheme (read by Time2Init): sdirk2 on a fresh run where accepted (no polar
  // boundary, integrator rk2), else be
  {auto *pmh = pmy_pack->pmesh;
  const std::string integ = pin->DoesParameterExist("time","integrator") ?
                            pin->GetString("time","integrator") : "rk2";
  const bool h2ok = !pmh->use_polar_boundary &&
                    (integ.compare("rk2") == 0);
  recipe::Old(pin, "rad_vet", "time_scheme", "be");
  // the old name hesdirk2 of sdirk2 (removed): stored -> NOTE, named -> FATAL
  if (pin->DoesParameterExist("rad_vet","time_scheme") &&
      pin->GetString("rad_vet","time_scheme").compare("hesdirk2") == 0) {
    recipe::Rule(pin, "rad_vet", "time_scheme", "sdirk2");
    pin->SetString("rad_vet", "time_scheme", "sdirk2");
  }
  (void) pin->GetOrAddString("rad_vet","time_scheme",
                             (h2ok && !global_variable::restart_run) ? "sdirk2" : "be");
  // be runs the sdirk2 backward-Euler path (t2_beonly) where sdirk2 is accepted
  t2_beonly = h2ok && pin->GetString("rad_vet","time_scheme").compare("be") == 0;}
  trans_on = false;
  trans_x3 = false;
  impl_cfl = pin->GetOrAddReal("rad_vet","implicit_cfl",-1.0);
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
  // implicit_pos_floor (default true)
  impl_pos_floor = pin->GetOrAddBoolean("rad_vet","implicit_pos_floor",
      recipe::Fix(pin, "rad_vet", "implicit_pos_floor", "true", "false"));
  Kokkos::realloc(pos_cnt_d, VET_POS_N);
  Kokkos::deep_copy(pos_cnt_d, 0.0);
  marshak_q = pin->GetOrAddReal("rad_vet","marshak_q",0.5);
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
  // the old scheme's predictor order: the channels of its restart files
  if (global_variable::restart_run) {
    impl_pord = pin->DoesParameterExist("rad_vet","implicit_predictor_order") ?
                pin->GetInteger("rad_vet","implicit_predictor_order") : 1;
  } else {
    recipe::Rule(pin, "rad_vet", "implicit_predictor_order", "2");
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
  // cubed sphere: the spin slope follows <hydro|mhd>/am_vector
  for (const char *b : {"hydro", "mhd"}) {
    am_vec = am_vec || (cs_geom && pin->DoesParameterExist(b, "am_vector") &&
                        pin->GetBoolean(b, "am_vector"));
  }
  MemArenaInit(pin, nmb, VET_NIW, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(iw, 0.0);
  if (trans_on) {
    Kokkos::realloc(f0x2, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2, 0.0);
    if (trans_x3) {
      Kokkos::realloc(f0x3, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3, 0.0);
    }
  }

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: transport=implicit"
              << (trans_x3 ? " (3-D)" : (trans_on ? " (2-D)" : " (1-D)"))
              << " implicit_cfl=" << impl_cfl << " marshak_q=" << marshak_q
              << std::endl;
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {ScCartInit(pin);}
  ImplBeamInit(pin);
  Time2Init(pin);
  if (time_scheme == VET_TIME_BE) {
    Kokkos::realloc(u1, u1.extent(0), VET_NVAR, 1, 1, 1);
  }
}

namespace {

const Real kT2G = imex_sdirk2::Gamma();   // shared tableau

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Init
//! \brief reads time_scheme (its default set by ImplicitInit) and the time2_* keys, and
//! allocates the stage state

void RadiationVet::Time2Init(ParameterInput *pin) {
  time_scheme = VET_TIME_BE;
  if (pin->DoesParameterExist("rad_vet", "time_scheme")) {
    std::string s = pin->GetString("rad_vet", "time_scheme");
    if (s.compare("be") == 0) {
      time_scheme = VET_TIME_BE;
    } else if (s.compare("sdirk2") == 0) {
      time_scheme = VET_TIME_SDIRK2;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/time_scheme = '" << s
                << "' is not a choice (be | sdirk2)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // be with t2_beonly runs as sdirk2 with every step on the backward-Euler path
  if (time_scheme == VET_TIME_BE && t2_beonly) {
    time_scheme = VET_TIME_SDIRK2;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> time_scheme = be: the sdirk2 backward-Euler path on "
                   "every step" << std::endl;
    }
  }
  if (time_scheme == VET_TIME_BE) {t2_beonly = false; return;}

  std::string integ = pin->GetString("time", "integrator");
  if (integ.compare("rk2") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time_scheme = sdirk2 needs <time>/integrator "
              << "= rk2 (its explicit part IS the Heun hydro)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // stage coefficients (utils/imex_sdirk2.hpp): SDIRK2, no explicit K1 (Time2GasTime).
  // time2_tableau is removed (hesdirk2, trbdf2 retired): a stored other value runs sdirk2
  // with a NOTE, a named one is FATAL
  recipe::Rule(pin, "rad_vet", "time2_tableau", "sdirk2");
  {
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make();
    t2_a21 = tb.a21;
    t2_w0 = tb.w0;
    t2_w1 = tb.w1;
    t2_c31 = tb.c31;
    t2_c32 = tb.c32;
    {
      // time2_gtime_slope: Time2GasTime's slope-limited exchange.  Default true
      const bool nm = pin->DoesParameterExist("rad_vet", "time2_gtime_slope");
      t2_gts = pin->GetOrAddBoolean("rad_vet", "time2_gtime_slope",
          recipe::Fix(pin, "rad_vet", "time2_gtime_slope", "true", "false"));
      recipe::OldStored(pin, "rad_vet", "time2_gtime_slope", "false", "true");
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_gtime_slope = " << (t2_gts ? "true" : "false")
                  << (nm ? "" : " (default)") << std::endl;
      }
      // time2_gtime_pos: the positive form of the flux state (always with
      // time2_gtime_slope).  Default true
      const bool np = pin->DoesParameterExist("rad_vet", "time2_gtime_pos");
      t2_gtp = pin->GetOrAddBoolean("rad_vet", "time2_gtime_pos",
          recipe::Fix(pin, "rad_vet", "time2_gtime_pos", "true", "false"));
      recipe::OldStored(pin, "rad_vet", "time2_gtime_pos", "false", "true");
      if (t2_gts && !t2_gtp) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_vet>/time2_gtime_slope = true needs "
                  << "time2_gtime_pos = true" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> time2_gtime_pos = " << (t2_gtp ? "true" : "false")
                  << (np ? "" : " (default)") << std::endl;
      }
    }
  }
  // a restart's time2_vet (alias time2_vet_col) = extrap2: the extrap2 channels of a file
  // of version < 4 (Time2RstNchWant)
  if (sph_geom && global_variable::restart_run) {
    const std::string kvc = VetKeyName(pin, "time2_vet", {"time2_vet_col"});
    t2_x2old = pin->DoesParameterExist("rad_vet", kvc) &&
               pin->GetString("rad_vet", kvc).compare("extrap2") == 0;
  }
  t2_dbg_fail = -1;
  if (pin->DoesParameterExist("rad_vet", "time2_dbg_fail")) {
    t2_dbg_fail = pin->GetInteger("rad_vet", "time2_dbg_fail");
  }

  // the stage increment's Zalesak limiter (Time2FaceFct): the clean solve stores the
  // energy flux Phi that moved E (WriteBack)
  t2_pfst = true;
  t2_ffct_emin = 0.1;   // (rad_vet.hpp)
  t2_abe = pin->GetOrAddInteger("rad_vet", "time2_auto_be", t2_abe);
  if (t2_abe < 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_vet>/time2_auto_be must be >= 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (global_variable::my_rank == 0 && t2_abe > 0) {
    std::cout << "<rad_vet> time2_auto_be=" << t2_abe << ": after a stage fallback the "
              << "next 1, 2, 4, ... " << t2_abe << " steps take the backward-Euler redo "
              << "path at once" << std::endl;
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(t2k1, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2k2, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::realloc(t2inc, nmb, VET_T2_NK, n3, n2, n1);
  Kokkos::deep_copy(t2k1, 0.0);
  Kokkos::deep_copy(t2k2, 0.0);
  Kokkos::deep_copy(t2inc, 0.0);
  if (t2_pfst) {  // zero: no face flux known (start, restart)
    Kokkos::realloc(t2pf1, nmb, 3, n3, n2, n1);
    Kokkos::realloc(t2pf2, nmb, 3, n3, n2, n1);
    Kokkos::deep_copy(t2pf1, 0.0);
    Kokkos::deep_copy(t2pf2, 0.0);
    t2pf1_ok = false;
    t2pf2_ok = false;
  }
  Kokkos::realloc(t2f1, nmb, n3, n2, n1+1);
  if (trans_on) {
    Kokkos::realloc(t2f2, nmb, n3, n2+1, n1);
    if (trans_x3) {Kokkos::realloc(t2f3, nmb, n3+1, n2, n1);}
  }
  if (vet_sc) {
    Kokkos::realloc(vet_prev, nmb, VET_T2_NVET, n3, n2, n1);
    Kokkos::deep_copy(vet_prev, 0.0);
  }
  t2_ok = false;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> time_scheme=sdirk2 (g=" << kT2G << "; first step, restarts "
              << "without the slope and failed stages take backward Euler)" << std::endl;
  }
}

}  // namespace radvet
