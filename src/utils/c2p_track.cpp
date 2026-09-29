//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file c2p_track.cpp
//! \brief the device kernels of the ConToPrim floor bookkeeping (see c2p_track.hpp)

#include "utils/c2p_track.hpp"

namespace c2p_track {

void Before(MeshBlockPack *pmbp, const DvceArray5D<Real> &u0,
            DvceArray4D<Real> re, DvceArray4D<Real> rd) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  par_for("c2p_track_b", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    re(m,k-ks,j-js,i-is) = u0(m,IEN,k,j,i);
    rd(m,k-ks,j-js,i-is) = u0(m,IDN,k,j,i);
  });
}

void After(MeshBlockPack *pmbp, const DvceArray5D<Real> &u0,
           DvceArray4D<Real> re, DvceArray4D<Real> rd, DvceArray2D<Real> acc,
           const Real w, const bool split) {
  const int r0 = split ? 2 : 0;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto vol = pmbp->pcoord->volume;
  par_for("c2p_track_a", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real de = u0(m,IEN,k,j,i) - re(m,k-ks,j-js,i-is);
    const Real dm = u0(m,IDN,k,j,i) - rd(m,k-ks,j-js,i-is);
    if (de != 0.0) Kokkos::atomic_add(&acc(r0,i-is), w*de*vol(m,k,j,i));
    if (dm != 0.0) Kokkos::atomic_add(&acc(r0+1,i-is), w*dm*vol(m,k,j,i));
  });
}

}  // namespace c2p_track
