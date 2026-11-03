//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter_rt_defaults.cpp
//! \brief prod1103c (USER 11-03 "we should also make some other settings default"): the
//! deep_hot_jupiter_rt FRESH-start defaults of keys that physics modules read before the
//! problem generator exists (time, fluid, resistivity, EOS).  main.cpp calls
//! recipe::PgenEarlyDefaults(pin) once, after the input and the command line are read and
//! before the Mesh is built; only this pgen's build links this file (other builds link
//! pgen/early_defaults_none.cpp).  Rules as every recipe default: applied only when
//! recipe::NewDefaults(pin) (fresh start, or <job>/recipe = prod1019); a key that is named
//! (deck, command line, restart input) always wins; the value is recorded (GetOrAdd) so
//! the run's restarts carry it; each decision is logged in the start-up defaults block.
//! A restart whose input lacks a key keeps the old default (nothing is added here).
//! The values are those of every Caltech dhj production deck (prod_1103/DEFAULTS_1103.txt
//! section A2).

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
  // the fluid block (MHD decks have <mhd>, hydro decks <hydro>), the mesh kind
  const std::string fb = pin->DoesBlockExist("mhd") ? "mhd" : "hydro";
  const bool fluid = pin->DoesBlockExist(fb);
  const bool cs = named("mesh", "use_cubed_sphere") &&
                  pin->GetBoolean("mesh", "use_cubed_sphere");
  // <time>/dt_min = 1e-3 s (the abort floor of every dhj production deck; old default 0)
  if (!named("time", "dt_min")) {
    pin->GetOrAddReal("time", "dt_min", 1.0e-3);
    Log("<time>/dt_min = 1e-3 (default on a fresh dhj start)");
  }
  // <fluid>/wellbalance_dynamic = true + wb_x1 = true + wb_option = polytropic (the x1
  // well-balanced reconstruction of every dhj production deck; old default off) on the
  // cubed sphere; each key only when not named, wb_option only when the scheme ends on
  if (fluid && cs) {
    for (const char *k : {"wellbalance_dynamic", "wb_x1"}) {
      if (named(fb, k)) continue;
      pin->GetOrAddBoolean(fb, k, true);
      Log("<" + fb + ">/" + k + " = true (default on a fresh cs dhj start)");
    }
    if (pin->GetBoolean(fb, "wellbalance_dynamic") && !named(fb, "wb_option")) {
      pin->GetOrAddString(fb, "wb_option", "polytropic");
      Log("<" + fb + ">/wb_option = polytropic (default with wellbalance_dynamic)");
    }
  }
}

}  // namespace recipe
