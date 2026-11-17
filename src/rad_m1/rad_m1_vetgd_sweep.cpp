//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_sweep.cpp
//! \brief vgdtu-1117: RadiationM1::VetGdSweep (the vet_gd shell sweep and its shell-body
//! kernels) in its own translation unit: rad_m1_vetgd.cpp compiled with M1_VGD_SWEEP_TU
//! keeps its includes, its file-local helpers and this one function.  The sweep's kernels
//! are the A100 build's longest nvcc/ptxas job; split off, the two halves build in parallel.

#define M1_VGD_SWEEP_TU
#include "rad_m1_vetgd.cpp"
