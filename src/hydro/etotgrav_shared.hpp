#ifndef HYDRO_ETOTGRAV_SHARED_HPP_
#define HYDRO_ETOTGRAV_SHARED_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file etotgrav_shared.hpp
//! \brief share-1104 P6: the <hydro|mhd>/etotgrav kernels, ONE implementation for both
//! fluids (defined in hydro_etotgrav.cpp).  Hydro::AddGravFlux/AddGravEtot/RemoveGravEtot
//! and MHD::AddGravFlux/AddGravEtot/RemoveGravEtot are thin wrappers.  The MHD energy
//! carries no extra term here: the gravitational energy is rho*Phi in both.

#include "athena.hpp"

class MeshBlockPack;

namespace fluid {
void EtotGravFlux(MeshBlockPack *pp, const DvceFaceFld4D<Real> &phi0,
                  DvceFaceFld5D<Real> &flx);
// sgn = +1 adds rho*Phi to cons(IEN), -1 removes it
void EtotGravCons(MeshBlockPack *pp, const DvceArray4D<Real> &phicc0,
                  DvceArray5D<Real> &cons, const int sgn, const int il, const int iu,
                  const int jl, const int ju, const int kl, const int ku);
}  // namespace fluid

#endif  // HYDRO_ETOTGRAV_SHARED_HPP_
