//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.cpp
//! \brief amdef-1019: <job>/recipe and the start-up defaults block (recipe.hpp)

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
            << "<job>/recipe = '" << r << "' is not a choice (prod1019 | none)" << std::endl;
  std::exit(EXIT_FAILURE);
}

bool NewDefaults(ParameterInput *pin) {
  const Mode m = GetMode(pin);
  if (m == Mode::none) return false;
  if (m == Mode::prod1019) return true;
  return !global_variable::restart_run;
}

void Log(const std::string &line) {Lines().push_back(line);}

void PrintLog() {
  auto &v = Lines();
  if (global_variable::my_rank == 0 && !v.empty()) {
    std::cout << "--- defaults (amdef-1019; <job>/recipe sets the default set: prod1019 | "
              << "none) ---" << std::endl;
    for (const auto &s : v) {std::cout << "  " << s << std::endl;}
    std::cout << "--- end defaults ---" << std::endl;
  }
  v.clear();
}
}  // namespace recipe
