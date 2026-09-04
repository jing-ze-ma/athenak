//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_seam_econsist.cpp
//! \brief CUBED SPHERE: make a seam ghost cell's TOTAL ENERGY consistent with its own
//! ghost magnetic field.
//!
//! THE DEFECT, measured by <mhd>/cs_seam_diag (mhd_seam_diag.cpp).  The cell-centred halo
//! carries the CONSERVED variables, so IEN arrives as the TOTAL energy with
//! 0.5|B_src|^2 already inside it, resampled across the seam AS A SCALAR.  The face field
//! arrives through a SEPARATE exchange and is resampled COMPONENTWISE, then squared by
//! ConsToPrim.  Resample-then-square is not square-then-resample, so
//!
//!     e_int = E_ghost - KE - 0.5|B_ghost|^2
//!
//! keeps the whole residual, divided by beta.  Measured on the stratified test: at a seam
//! below beta = 1 the internal-energy anomaly is 285x the density anomaly in the same
//! cells, and it grows 30x across beta = 1 while density's does not.
//!
//! THE FIX.  Exchange the FIELD-FREE part of the energy, E - 0.5|B_cc|^2, across the
//! seam, and rebuild the ghost's total energy from it and the ghost's OWN field:
//!
//!     E_ghost  <-  R[E_src - 0.5|B_src|^2]  +  0.5|B_cc,ghost|^2
//!
//! RESAMPLE THE SMALL QUANTITY, DO NOT SUBTRACT TWO LARGE ONES.  The first version of
//! this exchanged 0.5|B|^2 and corrected additively,
//! E_ghost <- E_ghost - R[0.5|B_src|^2] + 0.5|B_ghost|^2, which is algebraically the same
//! thing ONLY IF R IS LINEAR.  It is not: the along-seam resample carries a MONOTONE
//! CLAMP, and at low beta 0.5|B|^2 dominates E, so the clamp acts differently on the two
//! large numbers and most of the residual survives.  Measured: the additive form removed
//! 14 % of the seam low-beta anomaly, this one removes far more, for the same traffic.
//!
//! WHY IT NEEDS ITS OWN EXCHANGE, AFTER CT.  The magnetic energy that belongs inside
//! E^{n+1} is the one built from B^{n+1}, and CT does not run until AFTER SendU: at
//! SendU the face field is still B^n.  Correcting with B^n would trade an O(h^2)
//! inconsistency for an O(dt) one.  So the scalar rides its own cell-centred exchange
//! placed after RecvB.  That is safe against the existing halo traffic because every
//! MeshBoundaryValues object MPI_Comm_dup's its own communicator, so its tags -- which
//! encode only receiver lid and bufid -- cannot collide with the u0 or b0 exchange.
//!
//! WHAT IT DOES NOT FIX.  The kinetic term has the same structure: momentum is
//! transformed exactly and then resampled, so KE = |m|^2/2rho in the ghost is again
//! square-of-resample against the resample-of-square inside E.  That residual grows
//! with the Mach number rather than with 1/beta, and in a subsonic atmosphere it is far
//! smaller, so it is deliberately left alone here.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "mhd.hpp"
#include "coordinates/cell_locations.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn void MHD::FillSeamME
//! \brief fill me0 with E - 0.5|B_cc|^2 over the ACTIVE cells, using the post-CT field.
//!
//! The triple must be built EXACTLY as Coordinates::GnomonicEquiangleRaiseVelMHD builds
//! it -- orthonormal frame, and the stretched-grid weighting in the radial slot -- or the
//! value subtracted here is not the value ConsToPrim will add back, and the correction
//! would inject the difference instead of removing it.

void MHD::FillSeamME() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto me_ = me0;
  auto u0_ = u0;
  auto b0_ = b0;
  auto &ccell = pmy_pack->pcoord->cos_cell;
  auto &scell = pmy_pack->pcoord->sin_cell;
  const bool str_x1_ = (pmy_pack->pmesh->use_grid_stretch_r ||
                        pmy_pack->pmesh->use_grid_stretch_r_poly);
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;

  par_for("seam_me_fill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c = ccell(m,k,j);
    const Real sn = scell(m,k,j);
    const Real bx = str_x1_ ?
        CellCenteredRadialFld(b0_.x1f(m,k,j,i), b0_.x1f(m,k,j,i+1),
                              x1f_(m,i), x1f_(m,i+1), x1v_(m,i)) :
        0.5*(b0_.x1f(m,k,j,i) + b0_.x1f(m,k,j,i+1));
    const Real by_n = 0.5*(b0_.x2f(m,k,j,i) + b0_.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j,i));
    const Real by = (by_n + c*bz_n)/sn;
    const Real bz = bz_n;
    me_(m,0,k,j,i) = u0_(m,IEN,k,j,i) - 0.5*(bx*bx + by*by + bz*bz);
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::SeamEnergyFix
//! \brief swap the resampled source magnetic energy for the ghost cell's own.
//!
//! Applied over the TANGENTIAL ghost layers only -- x2 and x3, where the seams are.  The
//! radial ghosts are physical boundaries whose me0 was never exchanged, so touching them
//! would write whatever the buffer happened to hold.
//!
//! AND ONLY ACROSS A GENUINE SEAM.  Rebuilding as (E - 0.5B^2) + 0.5B^2 is not bitwise
//! identity in floating point, so applying it at a same-panel block face -- where the
//! halo is a plain copy and there is nothing to correct -- would perturb every ghost cell
//! by round-off for no reason.  A block whose x2min sits on the mesh x2min has its x2-lo
//! face ON a panel edge; that comparison needs no neighbour lookup.  A corner ghost
//! counts as a seam if EITHER of its two faces is one.

void MHD::SeamEnergyFix() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u0_ = u0;
  auto me_ = me0;
  auto b0_ = b0;
  auto &ccell = pmy_pack->pcoord->cos_cell;
  auto &scell = pmy_pack->pcoord->sin_cell;
  const bool str_x1_ = (pmy_pack->pmesh->use_grid_stretch_r ||
                        pmy_pack->pmesh->use_grid_stretch_r_poly);
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;

  const Real mx2min = pmy_pack->pmesh->mesh_size.x2min;
  const Real mx2max = pmy_pack->pmesh->mesh_size.x2max;
  const Real mx3min = pmy_pack->pmesh->mesh_size.x3min;
  const Real mx3max = pmy_pack->pmesh->mesh_size.x3max;
  const Real tol = 1.0e-10*(mx2max - mx2min);
  auto &size_ = pmy_pack->pmb->mb_size;

  par_for("seam_e_fix", DevExeSpace(), 0, nmb1, ks-ng, ke+ng, js-ng, je+ng, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // tangential ghost layers only
    if (j >= js && j <= je && k >= ks && k <= ke) { return; }
    bool seam = false;
    if (j < js) { seam = seam || (fabs(size_.d_view(m).x2min - mx2min) < tol); }
    if (j > je) { seam = seam || (fabs(size_.d_view(m).x2max - mx2max) < tol); }
    if (k < ks) { seam = seam || (fabs(size_.d_view(m).x3min - mx3min) < tol); }
    if (k > ke) { seam = seam || (fabs(size_.d_view(m).x3max - mx3max) < tol); }
    if (!seam) { return; }
    const Real c = ccell(m,k,j);
    const Real sn = scell(m,k,j);
    const Real bx = str_x1_ ?
        CellCenteredRadialFld(b0_.x1f(m,k,j,i), b0_.x1f(m,k,j,i+1),
                              x1f_(m,i), x1f_(m,i+1), x1v_(m,i)) :
        0.5*(b0_.x1f(m,k,j,i) + b0_.x1f(m,k,j,i+1));
    const Real by_n = 0.5*(b0_.x2f(m,k,j,i) + b0_.x2f(m,k,j+1,i));
    const Real bz_n = 0.5*(b0_.x3f(m,k,j,i) + b0_.x3f(m,k+1,j,i));
    const Real by = (by_n + c*bz_n)/sn;
    const Real bz = bz_n;
    u0_(m,IEN,k,j,i) = me_(m,0,k,j,i) + 0.5*(bx*bx + by*by + bz*bz);
  });
  return;
}

}  // namespace mhd
