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
extern int rst_mhd_nscalars;   // tracer-1101: the restart file's <mhd>/nscalars (user
                               // scalars only; -1: no restart or no <mhd>)
extern int rst_mhd_tracers;    // tracer-1101: its MHD tracers, bit 0 x_e, bit 1 floor
extern int rst_hydro_nscalars; // acctrc-1019: <hydro>/nscalars of the restart file's OWN
                               // input (-1: no restart or no <hydro>); see pgen.cpp
extern int rst_hydro_tracers;  // share-1104 P4: its hydro tracers (bit 1 floor), as
                               // rst_mhd_tracers
}

#endif // GLOBALS_HPP_
