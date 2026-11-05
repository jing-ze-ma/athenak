//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_etotgrav.cpp
//! \brief Implements functions for including time-independent gravitational potential energy in the energy conservation.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"
#include "hydro/etotgrav_shared.hpp"

namespace mhd {
// share-1104 P6: thin wrappers of the shared fluid::EtotGrav* kernels
void MHD::AddGravFlux(const DvceFaceFld4D<Real> &phi0, DvceFaceFld5D<Real> &flx) {
  fluid::EtotGravFlux(pmy_pack, phi0, flx);
}
void MHD::AddGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons,
                         const int il, const int iu, const int jl, const int ju,
                         const int kl, const int ku) {
  fluid::EtotGravCons(pmy_pack, phicc0, cons, +1, il, iu, jl, ju, kl, ku);
}
void MHD::RemoveGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons,
                            const int il, const int iu, const int jl, const int ju,
                            const int kl, const int ku) {
  fluid::EtotGravCons(pmy_pack, phicc0, cons, -1, il, iu, jl, ju, kl, ku);
}
}  // namespace mhd
