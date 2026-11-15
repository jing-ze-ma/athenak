//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_wellbalance.cpp
//! \brief Implements functions for deviation-based well-balanced scheme.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro.hpp"

// the body is shared with MHD (hydro/wellbalance_impl.inc)
#define WB_CLASS Hydro
#define WB_NVAR nhydro
namespace hydro {
#include "hydro/wellbalance_impl.inc"
} // namespace hydro
#undef WB_CLASS
#undef WB_NVAR
