#ifndef UTILS_RECIPE_HPP_
#define UTILS_RECIPE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.hpp
//! \brief amdef-1019 (USER 10-19): the default set of a run and <job>/recipe.
//!
//! The defaults introduced 10-19 (<hydro|mhd>/am_transport where supported,
//! <hydro|mhd>/etotgrav where supported, <rad_m1>/implicit_solver = bicgstab,
//! <rad_m1>/implicit_resid_fatal = 1e-2, problem/theta_wall_closed with theta_bc = hse;
//! <time>/dt_retry = true was already the fresh default) apply when NewDefaults():
//!   <job>/recipe not named    : a FRESH run (a restart whose input lacks a key keeps that
//!                               key's OLD default and says so: running chains unchanged)
//!   <job>/recipe = prod1019   : fresh AND restart (a preset: the same values; explicit
//!                               keys in the deck / on the command line still win)
//!   <job>/recipe = none       : never: every key takes its pre-10-19 default exactly
//! Every key a default decides is logged (Log) and printed in ONE block at start-up
//! (PrintLog, rank 0, after the problem generator).

#include <string>

class ParameterInput;

namespace recipe {
enum class Mode {standard, prod1019, none};
Mode GetMode(ParameterInput *pin);        // validates <job>/recipe (read without adding)
bool NewDefaults(ParameterInput *pin);    // see above
void Log(const std::string &line);        // one defaulted key (or why it stayed off)
void PrintLog();                          // the start-up block (rank 0), then clears
}  // namespace recipe

#endif  // UTILS_RECIPE_HPP_
