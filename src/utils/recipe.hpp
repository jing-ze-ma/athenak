#ifndef UTILS_RECIPE_HPP_
#define UTILS_RECIPE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file recipe.hpp
//! \brief the default set of a run.
//!
//! A FRESH run takes the new defaults; a restart whose input lacks a key keeps that key's
//! OLD default and says so (running chains unchanged).  An explicit key always wins.
//! Every key a default decides is logged (Log) and printed in ONE block at start-up
//! (PrintLog, rank 0, after the problem generator).

#include <string>

class ParameterInput;

namespace recipe {
bool NewDefaults(ParameterInput *pin);    // a fresh run (see above)
void CheckRemoved(ParameterInput *pin);   // FATAL if <job>/recipe is named
void Log(const std::string &line);        // one defaulted key (or why it stayed off)
bool Fix(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &value, const std::string &old);
// true: the fresh default of a key that corrects the scheme applies; a restart whose file
// lacks it takes it too, with a NOTE (naming the old value keeps the old behaviour)
void Old(ParameterInput *pin, const std::string &block, const std::string &key,
         const std::string &old);
// a restart whose file lacks the key keeps its old default (a state the file holds or the
// run's set-up): logged
void OldStored(ParameterInput *pin, const std::string &block, const std::string &key,
               const std::string &old, const std::string &value);
// NOTE: a restart runs the stored (or named) old value of such a key; how to switch
void PrintLog();                          // the start-up block (rank 0), then clears
// a problem generator's fresh-start defaults of keys read before the pgen
// exists (main.cpp, before the Mesh); pgen/<problem>_defaults.cpp or the no-op
// pgen/early_defaults_none.cpp
void PgenEarlyDefaults(ParameterInput *pin);
// a fresh VET run's defaults of keys read before <rad_vet> exists (main.cpp, after
// PgenEarlyDefaults): <rad_vet>/transport = implicit; with scheme = clean (named or
// the default with transport = implicit) <time>/cfl_number = 0.3, evolution = dynamic,
// <mesh>/nghost = 3, <hydro>/rsolver = hllc (rad_vet.cpp holds the rest)
void CleanDefaults(ParameterInput *pin);
// star_vet sets it: the star defaults of PgenEarlyDefaults (star_wb_defaults.cpp)
inline bool star_pgen = false;
}  // namespace recipe

#endif  // UTILS_RECIPE_HPP_
