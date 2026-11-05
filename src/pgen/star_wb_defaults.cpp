//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file star_wb_defaults.cpp
//! \brief defaults-1105 (USER 11-05, share-audit P12): recipe::PgenEarlyDefaults of the
//! he_star_m1 and ry_per_accretor builds (src/CMakeLists.txt links this file for those two
//! problems; deep_hot_jupiter_rt has its own, every other pgen the no-op
//! early_defaults_none.cpp).  The x1 well-balanced reconstruction of every he_star (HG,
//! RSG, R136a1) and accretor production deck, as the dhj G4 default:
//!   <fluid>/wellbalance_dynamic = true, <fluid>/wb_x1 = true, <fluid>/wb_option = polytropic
//! on the cubed sphere or the spherical-polar mesh.  Rules as every recipe default: only
//! when recipe::NewDefaults(pin) (fresh start, or <job>/recipe = prod1019); a named key
//! always wins (wellbalance_dynamic named false: nothing is added); the value is recorded
//! (GetOrAdd) so the run's restarts carry it; each decision is logged in the start-up
//! defaults block.  A restart whose input lacks a key keeps the old default (off).
//! <fluid>/wb_rmax is NOT defaulted (per-deck radius; 0 = WB everywhere).
//! (<rad_m1>/vet_rotate_advance, the other P12 key, is a rad_m1 default:
//! RadiationM1::VetRotInit.)

#include <string>

#include "athena.hpp"
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
}

}  // namespace recipe
