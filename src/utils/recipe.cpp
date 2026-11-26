//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.cpp
//! \brief the start-up defaults block (recipe.hpp)

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

bool NewDefaults(ParameterInput *pin) {return !global_variable::restart_run;}

void CheckRemoved(ParameterInput *pin) {
  if (!pin->DoesParameterExist("job", "recipe")) return;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<job>/recipe was removed (it is named with value '"
            << pin->GetString("job", "recipe") << "'): delete it and name every key the "
            << "run needs" << std::endl;
  std::exit(EXIT_FAILURE);
}

void Log(const std::string &line) {Lines().push_back(line);}

bool Fix(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &value, const std::string &old) {
  if (!global_variable::restart_run) return true;
  if (!pin->DoesParameterExist(block, key) && global_variable::my_rank == 0) {
    std::cout << "### NOTE: <" << block << ">/" << key << " = " << value << " (restart "
              << "file without it: the fix applies; name " << key << "=" << old
              << " to keep the old behaviour)" << std::endl;
  }
  return true;
}

void Old(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &old) {
  if (!global_variable::restart_run || pin->DoesParameterExist(block, key)) return;
  Log("<" + block + ">/" + key + " = " + old + " (restart file without it: the old "
      "default kept)");
}

void OldStored(ParameterInput *pin, const std::string &block, const std::string &key,
               const std::string &old, const std::string &value) {
  if (!global_variable::restart_run || global_variable::my_rank != 0 ||
      !pin->DoesParameterExist(block, key) || pin->GetString(block, key) != old) return;
  std::cout << "### NOTE: <" << block << ">/" << key << " = " << old << " (stored in the "
            << "restart file or named) is the old behaviour; name " << key << "="
            << value << " to switch" << std::endl;
}

void Rule(ParameterInput *pin, const std::string &block, const std::string &key,
          const std::string &value) {
  if (!pin->DoesParameterExist(block, key) ||
      ParameterInput::SameValue(pin->GetString(block, key), value)) return;
  for (auto &b : pin->block) {   // only the restart file holds it (not named): NOTE
    if (b.block_name != block || b.GetPtrToLine(key)->origin !=
        ParameterInput::kOriginRestart) continue;
    if (global_variable::my_rank == 0) {
      std::cout << "### NOTE: <" << block << "/" << key << "> = " << pin->GetString(block,
                key) << " is stored in the restart file (not named): the key was removed, "
                << "the code runs " << value << std::endl;
    }
    return;
  }
  std::cout << "### FATAL ERROR <" << block << "/" << key << "> = "
            << pin->GetString(block, key) << ": the key was removed, the code runs "
            << value << " here (name " << key << "=" << value << " to continue so)"
            << std::endl;
  std::exit(EXIT_FAILURE);
}

void CleanDefaults(ParameterInput *pin) {
  if (global_variable::restart_run || !pin->DoesBlockExist("rad_vet")) return;
  auto named = [&](const char *b, const char *k) {return pin->DoesParameterExist(b, k);};
  if (!named("rad_vet", "transport")) {
    pin->SetString("rad_vet", "transport", "implicit");
    Log("<rad_vet>/transport = implicit (the default of a fresh run)");
  }
  if (pin->GetString("rad_vet", "transport") != "implicit") return;
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
  add("time", "integrator", "rk2");
  add("mesh", "nghost", "3");
  if (pin->DoesBlockExist("hydro") && pin->GetString("time", "evolution") == "dynamic") {
    add("hydro", "rsolver", "hllc");
  }
}

void PrintLog() {
  auto &v = Lines();
  if (global_variable::my_rank == 0 && !v.empty()) {
    std::cout
        << "--- defaults ---"
        << std::endl;
    for (const auto &s : v) {std::cout << "  " << s << std::endl;}
    std::cout << "--- end defaults ---" << std::endl;
  }
  v.clear();
}
}  // namespace recipe
