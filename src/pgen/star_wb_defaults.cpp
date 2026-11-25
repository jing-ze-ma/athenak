//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file star_wb_defaults.cpp
//! \brief recipe::PgenEarlyDefaults of the
//! star_vet and ry_per_accretor builds (src/CMakeLists.txt links this file for those
//! two problems; deep_hot_jupiter_rt has its own, every other pgen the no-op
//! early_defaults_none.cpp).  The x1 well-balanced reconstruction of every star_vet
//! and accretor production deck:
//!   <fluid>/wellbalance_dynamic = true, <fluid>/wb_x1 = true, <fluid>/wb_option =
//! polytropic on the cubed sphere or the spherical-polar mesh.  Rules as every recipe
//! default: only when recipe::NewDefaults(pin) (fresh start);
//! a named key always wins (wellbalance_dynamic named false: nothing is added); the value
//! is recorded (GetOrAdd) so the run's restarts carry it; each decision is logged in the
//! start-up defaults block.  A restart whose input lacks a key keeps the old default
//! (off). <fluid>/wb_rmax is NOT defaulted (per-deck radius; 0 = WB everywhere).

#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "utils/recipe.hpp"

namespace recipe {

void PgenEarlyDefaults(ParameterInput *pin) {
  if (!NewDefaults(pin)) return;
  auto named = [pin](const std::string &b, const std::string &k) {
    return pin->DoesParameterExist(b, k);
  };
  const std::string fb = pin->DoesBlockExist("mhd") ? "mhd" : "hydro";
  if (!pin->DoesBlockExist(fb)) return;
  const bool cs = named("mesh", "use_cubed_sphere") &&
                  pin->GetBoolean("mesh", "use_cubed_sphere");
  const bool sp = named("mesh", "use_spherical_polar") &&
                  pin->GetBoolean("mesh", "use_spherical_polar");
  if (!(cs || sp)) {
    if (!named(fb, "wellbalance_dynamic")) {
      Log("<" + fb + ">/wellbalance_dynamic = false (default skipped: not on the cubed "
          "sphere / spherical polar grid)");
    }
    return;
  }
  const std::string geo = cs ? "cs" : "sp";
  if (named(fb, "wellbalance_dynamic") && !pin->GetBoolean(fb, "wellbalance_dynamic")) {
    return;   // explicitly off: no WB defaults at all
  }
  for (const char *k : {"wellbalance_dynamic", "wb_x1"}) {
    if (named(fb, k)) continue;
    pin->GetOrAddBoolean(fb, k, true);
    Log("<" + fb + ">/" + k + " = true (default on a fresh " + geo + " star start)");
  }
  if (!named(fb, "wb_option")) {
    pin->GetOrAddString(fb, "wb_option", "polytropic");
    Log("<" + fb + ">/wb_option = polytropic (default with wellbalance_dynamic)");
  }
  // star_vet with <rad_vet> (a fresh start): the well-balanced radiative envelope, the
  // steady start, the imposed inner flux, floored cells keep their velocity
  if (!star_pgen || !pin->DoesBlockExist("rad_vet") || global_variable::restart_run) {
    return;
  }
  auto add = [&](const std::string &b, const char *k, const char *v) {
    if (named(b, k)) return;
    pin->GetOrAddString(b, k, v);
    Log("<" + b + ">/" + k + " = " + v + " (default on a fresh star start)");
  };
  add("rad_vet", "force_reference", "wb_arad");
  add("rad_vet", "clean_init", "steady");
  if (named("rad_vet", "implicit_flux_x1min")) {
    add("rad_vet", "implicit_bc_x1min", "flux");
  }
  if (fb == "hydro") add("hydro", "dfloor_keep_velocity", "true");
}

}  // namespace recipe
