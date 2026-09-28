#ifndef UTILS_C2P_TRACK_HPP_
#define UTILS_C2P_TRACK_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file c2p_track.hpp
//! \brief floor bookkeeping shared by Hydro::C2PTrack and MHD::C2PTrack: the energy and
//! mass a ConToPrim call adds to the ACTIVE cells, volume integrated per x1 index.
//!
//! ConToPrim resets the conserved density and energy of a cell whose inversion hit a
//! floor.  Nothing books that change: it is neither a flux nor a source term, so an
//! energy budget built from the fluxes and sources misses it.  The owner (Hydro or MHD)
//! calls C2PTrackBefore just before the conversion and C2PTrackAfter just after it;
//! only the IDN and IEN differences of the active cells are summed.  Rows 0/1 (energy,
//! mass) take the conversions inside the RK stages, i.e. what the HYDRO update left
//! below a floor; rows 2/3 the operator-split conversion that follows the split
//! radiation (RTOpSplitBvals / Hydro::ConToPrimSplit), i.e. what the RADIATION step
//! left below a floor.  Off by default and
//! diagnostic only (the state is never written), so a run that does not enable it is
//! bitwise unchanged.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"

namespace c2p_track {

// allocate the accumulator (4, nx1) and the pre-conversion copies of the active cells
inline void Enable(MeshBlockPack *pmbp, DvceArray2D<Real> &acc, DvceArray4D<Real> &re,
                   DvceArray4D<Real> &rd) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int nmb = pmbp->nmb_thispack;
  const int n3 = indcs.nx3, n2 = indcs.nx2, n1 = indcs.nx1;
  if (acc.extent(1) != static_cast<size_t>(n1)) {
    Kokkos::realloc(acc, 4, n1);
    Kokkos::deep_copy(acc, 0.0);
  }
  Kokkos::realloc(re, nmb, n3, n2, n1);
  Kokkos::realloc(rd, nmb, n3, n2, n1);
}

// weight of a change made now in the value the step ends with: inside RK stage s the
// later stages keep gam0 of u0 each; outside the stages (pdrive null, stage 0, or the
// operator-split conversion) the change is carried in full
inline Real Weight(Driver *pdrive, const int stage, const bool split) {
  if (pdrive == nullptr || split || stage < 1 || stage > pdrive->nexp_stages) return 1.0;
  Real w = 1.0;
  for (int s = stage + 1; s <= pdrive->nexp_stages; ++s) w *= pdrive->gam0[s-1];
  return w;
}

inline void Before(MeshBlockPack *pmbp, const DvceArray5D<Real> &u0,
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

inline void After(MeshBlockPack *pmbp, const DvceArray5D<Real> &u0,
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

#endif  // UTILS_C2P_TRACK_HPP_
