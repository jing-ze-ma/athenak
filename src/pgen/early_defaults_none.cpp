//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file early_defaults_none.cpp
//! \brief prod1103c: recipe::PgenEarlyDefaults for every problem generator that has no
//! early (pre-physics) fresh-start defaults of its own: a no-op.  deep_hot_jupiter_rt
//! builds pgen/deep_hot_jupiter_rt_defaults.cpp instead (src/CMakeLists.txt).

#include "parameter_input.hpp"
#include "utils/recipe.hpp"

namespace recipe {
void PgenEarlyDefaults(ParameterInput *pin) { (void) pin; }
}  // namespace recipe
