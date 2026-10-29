#ifndef HYDRO_HYDRO_RINGAVG_HPP_
#define HYDRO_HYDRO_RINGAVG_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_ringavg.hpp
//! \brief ringavg-1029: conservative POLAR RING AVERAGING (PLUTO RING_AVERAGE style),
//! opt-in, for full-sphere spherical-polar meshes with polar boundaries.
//!
//! The theta rows next to a pole have a phi width r sin(theta) dphi that is a small
//! fraction of the equator's, and set the global hydro dt.  In each such row the phi
//! cells are taken in consecutive groups of N_j cells and, after every hydro stage
//! update (flux divergence + coordinate sources + am_transport, BEFORE the user sources,
//! floors, ConToPrim and the boundary exchange), the conserved state of each group is
//! replaced by its group average: rho, IM1, IM2, IM3 (componentwise in the local
//! (r, theta, phi) basis, so L_z = sum V R m_phi is kept exactly; all cells of a row have
//! the same volume), the passive scalars, and the total energy (with <hydro>/etotgrav the
//! rho Phi part is redistributed with the cell's own Phi: E_i = <E> + rho_bar (Phi_i -
//! <Phi>), which keeps sum E).  Hydro::NewTimeStep then uses N_j dphi as the phi width
//! in those rows.  Groups never straddle a MeshBlock (N_j divides the block nx3) and the
//! group boundaries are aligned across the pole (polar ghosts land on whole groups).
//! The radiation (rad_m1) is NOT averaged: the implicit transport has no phi CFL.
//!
//! Keys (nothing is constructed unless polar_ring_avg = true: default bitwise):
//!   <mesh>/polar_ring_avg     = true    needs use_spherical_polar + use_polar_boundary
//!   <mesh>/polar_ring_avg_n   = N_0,N_1,...  group size per row counted from EACH pole
//!                                 (override; rows beyond the list: 1).  Default: the
//!                                 largest power of 2 with N_j sin(theta_j) <=
//!                                 polar_ring_avg_tol, capped to divide the block nx3
//!   <mesh>/polar_ring_avg_tol = 1.05    (default rule only)

#include "athena.hpp"

class MeshBlockPack;
class ParameterInput;

namespace hydro {

class RingAvg {
 public:
  RingAvg(MeshBlockPack *ppack, ParameterInput *pin);

  // per (m, j) group size as a Real (1 outside the averaged rows), incl. x2 ghosts;
  // Hydro::NewTimeStep multiplies dx3 by it
  DvceArray2D<Real> nfac;

  //! u: hydro conserved (IDN, IM1..3, IEN, then nvar - nhyd scalars)
  void Average(DvceArray5D<Real> &u, const int nvar, const DvceArray4D<Real> &phi,
               const bool grav);

 private:
  MeshBlockPack *pmy_pack;
  int nent = 0;                 // number of averaged (block, row) entries on this rank
  DvceArray1D<int> em, ej, en;  // block, row j, group size of each entry
};

}  // namespace hydro
#endif  // HYDRO_HYDRO_RINGAVG_HPP_
