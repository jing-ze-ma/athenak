//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.cpp
//! \brief <job>/recipe and the start-up defaults block (recipe.hpp)

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "utils/recipe.hpp"

namespace recipe {
namespace {
std::vector<std::string> &Lines() {
  static std::vector<std::string> v;
  return v;
}
}  // namespace

Mode GetMode(ParameterInput *pin) {
  if (!pin->DoesParameterExist("job", "recipe")) return Mode::standard;
  const std::string r = pin->GetString("job", "recipe");
  if (r == "prod1019") return Mode::prod1019;
  if (r == "none") return Mode::none;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<job>/recipe = '" << r
            << "' is not a choice (prod1019 | none)" << std::endl;
  std::exit(EXIT_FAILURE);
}

bool NewDefaults(ParameterInput *pin) {
  const Mode m = GetMode(pin);
  if (m == Mode::none) return false;
  if (m == Mode::prod1019) return true;
  return !global_variable::restart_run;
}

void Log(const std::string &line) {Lines().push_back(line);}

void Old(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &value) {
  if (!global_variable::restart_run || pin->DoesParameterExist(block, key)) return;
  Log("<" + block + ">/" + key + " = " + value +
      " (a restart file without the key: the old default)");
}

void CleanDefaults(ParameterInput *pin) {
  if (global_variable::restart_run || !pin->DoesBlockExist("rad_vet")) return;
  auto named = [&](const char *b, const char *k) {return pin->DoesParameterExist(b, k);};
  if (!named("rad_vet", "transport")) {
    pin->SetString("rad_vet", "transport", "implicit");
    Log("<rad_vet>/transport = implicit (the default of a fresh run)");
  }
  const std::string sc = named("rad_vet", "scheme") ? pin->GetString("rad_vet",
                                                                     "scheme") :
      (pin->GetString("rad_vet", "transport") == "implicit" ? "clean" : "default");
  if (sc != "clean") return;
  // spherical polar and the cubed sphere: the vet_gd rays, unless the input names a
  // closure (or the analytic clean_closure)
  auto on = [&](const char *k) {return named("mesh", k) && pin->GetBoolean("mesh", k);};
  if ((on("use_spherical_polar") || on("use_cubed_sphere")) &&
      !named("rad_vet", "closure") && !named("rad_vet", "clean_closure")) {
    pin->SetString("rad_vet", "closure", "vet_col");
    if (!named("rad_vet", "vet_gd")) pin->SetBoolean("rad_vet", "vet_gd", true);
    Log("<rad_vet>/closure = vet_col, vet_gd = true (the default of a fresh clean run on "
        "spherical polar or the cubed sphere)");
  }
  if (!named("time", "cfl_number")) {
    pin->GetOrAddReal("time", "cfl_number", 0.3);
    Log("<time>/cfl_number = 0.3 (the default of a fresh clean VET run)");
  }
  auto add = [&](const char *b, const char *k, const char *v) {
    if (named(b, k)) return;
    pin->GetOrAddString(b, k, v);
    Log(std::string("<") + b + ">/" + k + " = " + v +
        " (the default of a fresh clean VET run)");
  };
  add("time", "evolution", "dynamic");
  add("mesh", "nghost", "3");
  if (pin->DoesBlockExist("hydro") && pin->GetString("time", "evolution") == "dynamic") {
    add("hydro", "rsolver", "hllc");
  }
}

void PrintLog() {
  auto &v = Lines();
  if (global_variable::my_rank == 0 && !v.empty()) {
    std::cout
        << "--- defaults (<job>/recipe sets the default set: prod1019 | none) ---"
        << std::endl;
    for (const auto &s : v) {std::cout << "  " << s << std::endl;}
    std::cout << "--- end defaults ---" << std::endl;
  }
  v.clear();
}
}  // namespace recipe
