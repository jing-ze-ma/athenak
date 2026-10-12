#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file globals.hpp
//  \brief namespace containing external global variables

namespace global_variable {
extern int my_rank, nranks;
extern bool restart_run;   // true when started with -r <restart_file>
extern bool final_outputs; // Driver::Finalize is writing its (unscheduled) outputs
}

#endif // GLOBALS_HPP_
