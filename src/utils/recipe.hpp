#ifndef UTILS_RECIPE_HPP_
#define UTILS_RECIPE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.hpp
//! \brief the default set of a run and <job>/recipe.
//!
//! The recipe defaults (<hydro|mhd>/am_transport where supported,
//! <hydro|mhd>/etotgrav where supported, <rad_vet>/implicit_solver = bicgstab,
//! <rad_vet>/implicit_resid_fatal = 1e-2, problem/theta_wall_closed with theta_bc = hse;
//! <time>/dt_retry = true was already the fresh default) apply when NewDefaults():
//!   <job>/recipe not named    : a FRESH run (a restart whose input lacks a key keeps
//! that key's OLD default and says so: running chains unchanged)
//!   <job>/recipe = prod1019   : fresh AND restart (a preset: the same values; explicit
//!                               keys in the deck / on the command line still win)
//!   <job>/recipe = none       : never: every key takes its plain default exactly
//! Every key a default decides is logged (Log) and printed in ONE block at start-up
//! (PrintLog, rank 0, after the problem generator).

#include <string>

class ParameterInput;

namespace recipe {
enum class Mode {standard, prod1019, none};
Mode GetMode(ParameterInput *pin);        // validates <job>/recipe (read without adding)
bool NewDefaults(ParameterInput *pin);    // see above
void Log(const std::string &line);        // one defaulted key (or why it stayed off)
void Old(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &value);       // Log: a restart file without key runs value
void PrintLog();                          // the start-up block (rank 0), then clears
// a problem generator's fresh-start defaults of keys read before the pgen
// exists (main.cpp, before the Mesh); pgen/<problem>_defaults.cpp or the no-op
// pgen/early_defaults_none.cpp
void PgenEarlyDefaults(ParameterInput *pin);
// a fresh M1 run's defaults of keys read before <rad_vet> exists (main.cpp, after
// PgenEarlyDefaults): <rad_vet>/transport = implicit; with scheme = clean (named or
// the default with transport = implicit) closure = vet_col + vet_gd on spherical polar
// and the cubed sphere and <time>/cfl_number = 0.3 (rad_vet.cpp holds the rest)
void CleanDefaults(ParameterInput *pin);
}  // namespace recipe

#endif  // UTILS_RECIPE_HPP_
