//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetgd_sweep.cpp
//! \brief RadiationVet::VetGdSweep from rad_vet_vetgd.cpp (VET_VGD_SWEEP_TU); split off
//! because its kernels are the longest nvcc/ptxas job, so the halves build in parallel.

#define VET_VGD_SWEEP_TU
#include "rad_vet_vetgd.cpp"
