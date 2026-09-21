//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel_frozen.cpp
//! \brief Coordinates::GnomonicEquiangleRaiseVel{,MHD}Frozen(), the CUBED-SPHERE half of
//! the conserved-to-primitive conversion used ONCE, by the first ConToPrim of a RESTARTED
//! run under a general EOS.
//!
//! WHY THEY EXIST.  general_hyd_frozen.cpp / general_mhd_frozen.cpp froze the thermo-
//! dynamic cache inside ConsToPrim, but on the cubed sphere ConsToPrim is only half of
//! the inversion: Coordinates::GnomonicEquiangleRaiseVel{,MHD} runs immediately after it
//! and is the routine that actually OWNS the general-EOS cache there.  It re-solves
//! wtemp and wder from the metric-corrected internal energy (the `if (gen_)` block of
//! gnomonic_raisevel.hpp / gnomonic_raisevel_mhd.hpp, reached from
//! src/hydro/hydro_tasks.cpp:966 and src/mhd/mhd_tasks.cpp:735), and it re-applies the
//! DEFERRED floors and the velocity ceiling.  So on that grid the frozen ConsToPrim
//! bought nothing: the cache the restart file restored was overwritten a moment later by
//! a fresh table inversion, and two of the deferred floors read `Hydro::dfl_fv`, which
//! the frozen ConsToPrim does not write at all, so they ran against the PREVIOUS cycle's
//! density-floor record.
//!
//! WHAT THEY DO, AND THE ORDER ARGUMENT.  The state the restart file holds was written
//! AFTER the straight run's own raise-velocity step, floors, ceiling and cache solve.
//! Reproducing that state therefore means redoing the ALGEBRA of the raise -- raising the
//! covariant momentum with the metric, forming the metric kinetic energy (and, for MHD,
//! the cell-centred field in the orthonormal frame and its magnetic energy) and taking
//! the internal energy off the conserved total -- and NOTHING else: no floor, no ceiling,
//! no inversion, no write to u0, wtemp or wder.  Every floor test would be false anyway,
//! because the conserved state in the file is the one those floors already produced, and
//! re-running them could only move the answer; the cache in the file is already the
//! answer the straight run's inversion gave, and re-solving it from its own output moves
//! it by an ULP (the table root find stops on a step size, so its result depends on the
//! guess it started from).  The next ConToPrim, at the end of the first stage, is the
//! ordinary one again -- Hydro::c2p_freeze_derived / MHD::c2p_freeze_derived are cleared
//! by the call that honours them.
//!
//! WHY THIS IS A FILE OF ITS OWN.  The same reason general_hyd_frozen.cpp is: adding a
//! kernel to a translation unit moves how the device compiler contracts the arithmetic of
//! the kernels already in it, and coordinates.cpp holds the DEFAULT cubed-sphere path,
//! whose answer must not move by a unit in the last place.
//!
//! The arithmetic below is copied from GnomonicRaiseVelFloors / GnomonicRaiseVelMHDFloors
//! and from the cell-centred field block of Coordinates::GnomonicEquiangleRaiseVelMHD, in
//! the same expressions and the same order.  Keep them in step.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleRaiseVelFrozen
//! \brief the algebraic half of GnomonicEquiangleRaiseVel, with the cache frozen.

void Coordinates::GnomonicEquiangleRaiseVelFrozen(DvceArray5D<Real> &u0,
    DvceArray5D<Real> &w0, const int il, const int iu, const int jl, const int ju,
    const int kl, const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  par_for("cs_raisev_frozen", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real det = 1.0 - c*c;
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    const Real m2 = u0(m,IM2,k,j,i);   // xi
    const Real m3 = u0(m,IM3,k,j,i);   // eta
    // v^i = g^{ij} m_j / rho, with the metric acting on the ANGULAR pair only
    const Real v1 = m1/d;
    const Real v2 = (m2 - c*m3)/(d*det);
    const Real v3 = (m3 - c*m2)/(d*det);
    // KE = 0.5 rho g_ij v^i v^j = 0.5 m_i v^i, which is the cross-term-correct form.
    const Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    w0(m,IEN,k,j,i) = u0(m,IEN,k,j,i) - ekin;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleRaiseVelMHDFrozen
//! \brief the algebraic half of GnomonicEquiangleRaiseVelMHD, with the cache frozen.

void Coordinates::GnomonicEquiangleRaiseVelMHDFrozen(DvceArray5D<Real> &u0,
    const DvceFaceFld4D<Real> &b0, DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
    const int il, const int iu, const int jl, const int ju, const int kl,
    const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  auto &sin_cell_ = sin_cell;
  // the cell centre is the volume CENTROID on every radial grid, so the radial field at
  // the centre is the weighted interpolation, unconditionally -- exactly the branch
  // GnomonicEquiangleRaiseVelMHD takes (its str_x1_ is `true`)
  auto &x1v_ = x1v;
  auto &x1f_ = xx1f;
  par_for("cs_raisev_mhd_frozen", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real sn = sin_cell_(m,k,j);
    // cell-centred field, then into the orthonormal frame; see the note in
    // Coordinates::GnomonicEquiangleRaiseVelMHD, whose expressions these are
    const Real bx = CellCenteredRadialFld(b0.x1f(m,k,j,i), b0.x1f(m,k,j,i+1),
                                          x1f_(m,i), x1f_(m,i+1), x1v_(m,i));
    const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
    const Real by = (by_n + c*bz_n)/sn;
    const Real bz = bz_n;
    bcc0(m,IBX,k,j,i) = bx;
    bcc0(m,IBY,k,j,i) = by;
    bcc0(m,IBZ,k,j,i) = bz;

    const Real det = 1.0 - c*c;
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i);   // radial: orthogonal to both angles
    const Real m2 = u0(m,IM2,k,j,i);   // xi
    const Real m3 = u0(m,IM3,k,j,i);   // eta
    const Real v1 = m1/d;
    const Real v2 = (m2 - c*m3)/(d*det);
    const Real v3 = (m3 - c*m2)/(d*det);
    // the two energies are kept separate and subtracted in THIS order, as
    // GnomonicRaiseVelMHDFloors does: folding them into one sum re-associates the
    // rounding.
    const Real ekin = 0.5*(m1*v1 + m2*v2 + m3*v3);
    const Real emag = 0.5*(bx*bx + by*by + bz*bz);
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    w0(m,IEN,k,j,i) = u0(m,IEN,k,j,i) - ekin - emag;
  });
  return;
}
