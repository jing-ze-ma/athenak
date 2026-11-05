#ifndef HYDRO_RK_UPDATE_SHARED_HPP_
#define HYDRO_RK_UPDATE_SHARED_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rk_update_shared.hpp
//! \brief share-1104 P6: the explicit RK stage update shared by Hydro and MHD (defined in
//! hydro/hydro_update.cpp).

#include "athena.hpp"

class MeshBlockPack;
class Driver;

namespace fluid {
void RKUpdateCore(MeshBlockPack *pmy_pack, Driver *pdriver, const int stage,
                  const int nvar, const DvceFaceFld5D<Real> &uflx,
                  DvceArray5D<Real> &u0, const DvceArray5D<Real> &u1);
}  // namespace fluid

#endif  // HYDRO_RK_UPDATE_SHARED_HPP_
