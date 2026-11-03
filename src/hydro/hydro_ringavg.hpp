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
//! The radiation (rad_m1) is averaged only with <rad_m1>/polar_ring_avg = true
//! (rad_m1_ringavg.cpp: the implicit transport of each group as ONE coarse cell).
//!
//! Keys (nothing is constructed unless polar_ring_avg = true: default bitwise):
//!   <mesh>/polar_ring_avg     = true    needs use_spherical_polar + use_polar_boundary
//!   <mesh>/polar_ring_avg_n   = N_0,N_1,...  group size per row counted from EACH pole
//!                                 (override; rows beyond the list: 1).  Default: the
//!                                 largest power of 2 with N_j sin(theta_j) <=
//!                                 polar_ring_avg_tol, capped to divide the block nx3
//!   <mesh>/polar_ring_avg_tol = 1.05    (default rule only)
//!
//! ringsmooth-1102 (opt-in, read only when named; absent = the group average, bitwise):
//!   <mesh>/polar_ring_avg_filter = fourier   instead of the piecewise-constant group
//!       average, a conservative FOURIER TRUNCATION of each whole ring (all phi blocks, all
//!       ranks: partial DFT sums + MPI_Allreduce): rho, IM1..3, the scalars and the internal
//!       energy e = E - KE (- rho Phi with etotgrav) keep their modes m <= m_c(N_j), m_c =
//!       floor(mfac (nphi/pi) asin(1/N_j)) (the largest m whose fine-grid phase speed
//!       sin(pi m/nphi)/dphi fits the group-width CFL: the same stability margin as the
//!       unaveraged rows), all higher m are removed.  The state stays SMOOTH in phi (no
//!       staircase), the low-m content (tidal m = 1, 2 equilibrium, cross-pole flow) is
//!       untouched.  Energy: the ring-sum defect dE = sum E - sum E_f (the kinetic and rho
//!       Phi energy of the removed modes) is put back UNIFORMLY on the ring as heat (or, after
//!       the M1 solve with polar_ring_avg_ke = rad, into the radiation E), never as velocity.
//!       Mass, momentum components, L_z and total energy are conserved per ring.  A ring
//!       (pole, row, i) whose filtered rho or e falls below thr x the old value somewhere
//!       falls back to the group average (counted, reported at the end).
//!   <mesh>/polar_ring_avg_mfac = 1.0, <mesh>/polar_ring_avg_thr = 0.5  (fourier only)

#include <vector>

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

  // the group table, read by the radiation (<rad_m1>/polar_ring_avg, rad_m1_ringavg.cpp)
  int nent = 0;                 // number of averaged (block, row) entries on this rank
  DvceArray1D<int> em, ej, en;  // block, row j, group size of each entry

  // ringsmooth-1102 <mesh>/polar_ring_avg_filter = fourier (hydro_ringfourier.cpp)
  bool fourier = false;
  //! u as Average; urad != nullptr: the energy defect goes to urad(m, irad) x hc (M1 E)
  //! instead of gas heat
  void Filter(DvceArray5D<Real> &u, const int nvar, const DvceArray4D<Real> &phi,
              const bool grav, DvceArray5D<Real> *urad, const int irad, const Real hc);
  ~RingAvg();

 private:
  MeshBlockPack *pmy_pack;
  // ringsmooth-1102 fourier tables (FourierInit)
  void FourierInit(ParameterInput *pin, const std::vector<int> &vjp);
  void FourierAlloc(const int nvar);
  Real fthr = 0.5;
  int nphi = 0, nx1g = 0, nring = 0, nslot = 0, fnvar = 0, nfoff = 0, fmcmax = 0;
  DvceArray1D<int> ering, ek0, eig0;          // per entry: ring, global k0, global i0
  DvceArray1D<int> rmc, roff;                 // per ring: cutoff m_c, coef offset (nvar 0)
  DvceArray1D<int> sring, sig0, sbeg, slist;  // slots (ring, i0): CSR of local entries
  DvceArray1D<Real> tcos, tsin;               // cos, sin(2 pi q/nphi)
  DualArray1D<Real> coef, rsum, rmin;         // ring coefficients; ring (sum E, sum E_f); mins
  DvceArray4D<Real> scr;                      // (entry*nvar, k, i) filtered state
  DvceArray3D<Real> part;                     // (entry, i, 4) partial sums / mins
  double nfall = 0.0, ncall = 0.0;
};

}  // namespace hydro
#endif  // HYDRO_HYDRO_RINGAVG_HPP_
