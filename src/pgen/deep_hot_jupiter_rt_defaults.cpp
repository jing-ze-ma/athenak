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
  const bool sp = named("mesh", "use_spherical_polar") &&
                  pin->GetBoolean("mesh", "use_spherical_polar");
  const std::string geo = cs ? "cs" : (sp ? "sp" : "Cartesian");
  // <time>/dt_min = 1e-3 s (the abort floor of every dhj production deck; old default 0)
  if (!named("time", "dt_min")) {
    pin->GetOrAddReal("time", "dt_min", 1.0e-3);
    Log("<time>/dt_min = 1e-3 (default on a fresh dhj start)");
  }
  // <fluid>/wellbalance_dynamic = true + wb_x1 = true + wb_option = polytropic (the x1
  // well-balanced reconstruction of every dhj production deck; old default off) on every
  // dhj mesh (share-1104 G4: the K&M x1 scheme has no geometry restriction -- cubed
  // sphere, spherical polar, and the Cartesian plane-parallel box with constant g, where
  // x1 is the height); each key only when not named, wb_option only when the scheme ends on
  if (fluid) {
    for (const char *k : {"wellbalance_dynamic", "wb_x1"}) {
      if (named(fb, k)) continue;
      pin->GetOrAddBoolean(fb, k, true);
      Log("<" + fb + ">/" + k + " = true (default on a fresh " + geo + " dhj start)");
    }
    if (pin->GetBoolean(fb, "wellbalance_dynamic") && !named(fb, "wb_option")) {
      pin->GetOrAddString(fb, "wb_option", "polytropic");
      Log("<" + fb + ">/wb_option = polytropic (default with wellbalance_dynamic)");
    }
  }
  // <fluid>/am_vector = true with problem/am_transport = false (the angular-momentum
  // vector form + spin reservoir, as every dhj production deck; the pgen correction must
  // be off with it, YieldToPgen FATAL otherwise).  Only as the PAIR, on the cubed sphere
  // or the spherical-polar mesh (share-1104 G4; AmModule's own exclusions: 1-D sp and
  // 2-D sp with sp_cart_polar_momentum; no angular momentum on the Cartesian box), when
  // neither the fluid's AM keys nor problem/am_transport are named.
  if (fluid && !named(fb, "am_vector") && !named(fb, "am_transport") &&
      !named(fb, "frame_omega") && !named("problem", "am_transport")) {
    std::string why;
    if (!(cs || sp)) {
      why = "Cartesian mesh";
    } else if (sp) {
      const int nx2 = named("mesh", "nx2") ? pin->GetInteger("mesh", "nx2") : 1;
      const int nx3 = named("mesh", "nx3") ? pin->GetInteger("mesh", "nx3") : 1;
      auto on = [pin, &named, &fb](const char *k) {
        return named(fb, k) && pin->GetBoolean(fb, k);
      };
      if (nx2 < 2) {
        why = "1-D spherical polar";
      } else if (nx3 < 2 && (on("sp_cart_polar_momentum") || on("sp_cart_all_momentum"))) {
        why = "2-D spherical polar with sp_cart_polar_momentum";
      }
    }
    if (why.empty()) {
      pin->GetOrAddBoolean(fb, "am_vector", true);
      pin->GetOrAddBoolean("problem", "am_transport", false);
      Log("<" + fb + ">/am_vector = true + problem/am_transport = false (default pair on a "
          "fresh " + geo + " dhj start)");
    } else {
      Log("<" + fb + ">/am_vector = false (default pair skipped: " + why + ")");
    }
  }
  // <mhd>/max_eta = 1e15 cm^2/s with ohmic_resistivity = eos (a required key until now;
  // user 10-03 test 1e13/1e14/1e15, every dhj MHD production deck)
  if (fb == "mhd" && named("mhd", "ohmic_resistivity") &&
      pin->GetString("mhd", "ohmic_resistivity") == "eos" && !named("mhd", "max_eta")) {
    pin->GetOrAddReal("mhd", "max_eta", 1.0e15);
    Log("<mhd>/max_eta = 1e15 (default with ohmic_resistivity = eos on a fresh dhj start)");
  }
  // <mhd>/valfven_max_cgs = 2e6 cm/s (the Alfven-speed density floor of every dhj MHD
  // production deck; old: off).  Also switches on its own defaults in mhd.cpp
  // (valfven_noinfall_fac 4, floor_tracer), as a named key does.
  if (fb == "mhd" && fluid && !named("mhd", "valfven_max_cgs")) {
    pin->GetOrAddReal("mhd", "valfven_max_cgs", 2.0e6);
    Log("<mhd>/valfven_max_cgs = 2e6 (default on a fresh dhj MHD start)");
  }
  // <mhd>/use_rkg_sts = true (RKG2 super-time-stepping of the resistive terms, every dhj
  // MHD production deck; USER 11-03 overrides the 11-01 "keep false" for this scope) with
  // resistive dhj MHD: ohmic_resistivity named, or ambipolar named and not none.  Its
  // own fresh defaults (rkg_lean, rkg_dt_safety, rkg_strang, ...) follow in resistivity.cpp.
  if (fb == "mhd" && fluid && !named("mhd", "use_rkg_sts") &&
      (named("mhd", "ohmic_resistivity") ||
       (named("mhd", "ambipolar") && pin->GetString("mhd", "ambipolar") != "none"))) {
    pin->GetOrAddBoolean("mhd", "use_rkg_sts", true);
    Log("<mhd>/use_rkg_sts = true (default with resistive dhj MHD on a fresh start)");
  }
  // <fluid>/wb_plm_range = true, wb_plm_smooth = true, top_ghost_prim_copy = true
  // (wbplm-1103; USER 11-03) for the dhj well-balanced x1 PLM on the non-uniform x1 grid
  // (cubed sphere or spherical polar, wellbalance_dynamic + wb_x1, reconstruct plm) with the
  // outflow1 top (problem/top_bc outflow1, the fresh default; top_ghost_prim_copy acts
  // only there).  Not supported elsewhere (the keys FATAL off the non-uniform WB PLM):
  // then nothing is added and the reason is logged.  Each key only when not named.
  {
    const bool wb = fluid && named(fb, "wellbalance_dynamic") &&
                    pin->GetBoolean(fb, "wellbalance_dynamic") && named(fb, "wb_x1") &&
                    pin->GetBoolean(fb, "wb_x1");
    const bool plm = !named(fb, "reconstruct") || pin->GetString(fb, "reconstruct") == "plm";
    const bool top1 = !named("problem", "top_bc") ||
                      pin->GetString("problem", "top_bc") == "outflow1";
    const char *ks[3] = {"wb_plm_range", "wb_plm_smooth", "top_ghost_prim_copy"};
    std::string why;
    if (!(cs || sp)) why = "not on the cubed sphere / spherical polar grid";
    else if (!wb) why = "no x1 well-balanced reconstruction";
    else if (!plm) why = "reconstruct is not plm";
    else if (!top1) why = "problem/top_bc is not outflow1";
    for (int q = 0; q < 3; ++q) {
      if (named(fb, ks[q])) continue;
      if (why.empty()) {
        pin->GetOrAddBoolean(fb, ks[q], true);
        Log("<" + fb + ">/" + ks[q] + " = true (default: dhj well-balanced PLM, outflow1 top)");
      } else if (fluid) {
        Log("<" + fb + ">/" + ks[q] + " = false (default skipped: " + why + ")");
      }
    }
  }
}

}  // namespace recipe
