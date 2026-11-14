//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_raisevel_frozen.cpp
//! \brief Coordinates::GnomonicEquiangleRaiseVel{,MHD}Frozen(): the metric part of the
//! cubed-sphere cons-to-prim conversion used by the first ConToPrim of a restarted run
//! under a general EOS. Only the algebra is redone (raise the covariant momentum, form
//! the metric kinetic and magnetic energies, subtract them from the total energy);
//! floors, ceiling and EOS inversion are skipped so the restored EOS cache is kept.
//! The expressions match GnomonicRaiseVelFloors / GnomonicRaiseVelMHDFloors. Kept in a
//! separate translation unit so it does not affect code generation of coordinates.cpp.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"

//----------------------------------------------------------------------------------------
//! \fn Coordinates::GnomonicEquiangleRaiseVelFrozen
//! \brief Algebraic part of GnomonicEquiangleRaiseVel with the EOS cache frozen.

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
    // contravariant velocity v^i = g^{ij} m_j / rho
    const Real v1 = m1/d;
    const Real v2 = (m2 - c*m3)/(d*det);
    const Real v3 = (m3 - c*m2)/(d*det);
    // kinetic energy 0.5 m_i v^i
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
//! \brief Algebraic part of GnomonicEquiangleRaiseVelMHD with the EOS cache frozen.

void Coordinates::GnomonicEquiangleRaiseVelMHDFrozen(DvceArray5D<Real> &u0,
    const DvceFaceFld4D<Real> &b0, DvceArray5D<Real> &bcc0, DvceArray5D<Real> &w0,
    const int il, const int iu, const int jl, const int ju, const int kl,
    const int ku) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cos_cell_ = cos_cell;
  auto &sin_cell_ = sin_cell;
  // radial field at the cell centroid by weighted interpolation of the faces
  // the cell centre is the volume CENTROID on every radial grid, so the radial field at
  // the centre is the weighted interpolation, unconditionally -- exactly the branch
  // GnomonicEquiangleRaiseVelMHD takes (its str_x1_ is `true`)
  auto &x1v_ = x1v;
  auto &x1f_ = xx1f;
  par_for("cs_raisev_mhd_frozen", DevExeSpace(), 0,nmb1, kl,ku, jl,ju, il,iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = cos_cell_(m,k,j);
    const Real sn = sin_cell_(m,k,j);
    // cell-centred field in the orthonormal frame (as GnomonicEquiangleRaiseVelMHD)
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
    // kinetic and magnetic energies subtracted separately, as in the floor function
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
