//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_boris.cpp
//! \brief boris-1102 <mhd>/boris_c: the cell part of the Alfven-speed (simplified Boris)
//! limiter.  See the note on EOS_Data::boris_c2inv (eos.hpp) for the formulation and
//! rsolvers/boris_hlle_mhd.hpp for the face part.
//!
//! Each RK stage writes u0 = gam0 u0 + gam1 u1 - beta dt div F + beta dt S.  RKUpdate
//! first saves the momentum reference M_ref = gam0 M0 + gam1 M1 (BorisSaveRef); after
//! ALL momentum sources of the stage (end of MHDSrcTerms) the stage's momentum change
//! dM = M - M_ref is split along the unit field b of the stage-start state and its
//! perpendicular part multiplied by f:
//!   dM' = f dM + (1 - f) (b.dM) b,     f = 1/sqrt(1 + (v_A/c_eff)^4).
//! Because the change of EVERYTHING (pressure gradient, Lorentz force, gravity, Coriolis,
//! tide, ...) is scaled, a state in force balance (dM = 0) stays in balance: the steady
//! state is unchanged.  The total energy is corrected by the kinetic-energy difference,
//!   E += (|M_ref + dM'|^2 - |M_ref + dM|^2)/(2 rho),
//! so the internal energy (and T) evolves exactly as it would from the same momentum:
//! the energy the Lorentz force puts into perpendicular motion but that the increased
//! inertia does not let become kinetic energy is the electric-field energy
//! rho v_perp^2 (1/f - 1)/2 of semi-relativistic MHD, which the code's E does not carry.
//! Non-conservation (both measured, printed every ndiag cycles): momentum by the
//! perpendicular (1 - f) dM_perp (zero for steady states), total energy by minus the
//! change of that electric energy (bounded, oscillatory in waves).
//! On the cubed sphere the momentum is covariant on the non-orthogonal gnomonic basis;
//! the split is done in the orthonormal frame {r, e_xi, (e_eta - c e_xi)/s} in which
//! bcc0 is stored (Coordinates::GnomonicEquiangleRaiseVelMHD).  Cartesian and spherical
//! polar momenta are orthonormal already.

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "coordinates/coordinates.hpp"
#include "mhd.hpp"

namespace mhd {

void MHD::BorisInit(ParameterInput *pin) {
  auto &eos = peos->eos_data;
  if (!eos.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/boris_c needs an ideal or general EOS (not isothermal)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (use_fofc || rad_signal_speed) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/boris_c is not implemented with fofc or rad_signal_speed"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!(rsolver_method == MHD_RSolver::hlle || rsolver_method == MHD_RSolver::hlld ||
        rsolver_method == MHD_RSolver::lhlld || rsolver_method == MHD_RSolver::llf)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/boris_c needs rsolver = llf | hlle | hlld | lhlld" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  boris_on = true;
  boris_ndiag = pin->GetOrAddInteger("time","ndiag",1);
  if (boris_ndiag < 1) boris_ndiag = 1;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(boris_mref, pmy_pack->nmb_thispack, 3, ncells3, ncells2, ncells1);
}

//----------------------------------------------------------------------------------------
//! \fn MHD::BorisSaveRef: M_ref = gam0 M0 + gam1 M1 (the stage's zero-RHS momentum)

void MHD::BorisSaveRef(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const Real gam0 = pdriver->gam0[stage-1];
  const Real gam1 = pdriver->gam1[stage-1];
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u0_ = u0;
  auto u1_ = u1;
  auto mref = boris_mref;
  par_for("boris_ref", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int n = 0; n < 3; ++n) {
      mref(m,n,k,j,i) = gam0*u0_(m,IM1+n,k,j,i) + gam1*u1_(m,IM1+n,k,j,i);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn MHD::BorisCorrect: dM_perp -> f dM_perp, E corrected by the KE difference

void MHD::BorisCorrect(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto u0_ = u0;
  auto w0_ = w0;
  auto bcc_ = bcc0;
  auto wder_ = wder;
  auto mref = boris_mref;
  auto eos = peos->eos_data;
  const bool gen = eos.IsGeneral();
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool curv = cs || pmy_pack->pmesh->use_spherical_polar;
  auto cosc = pmy_pack->pcoord->cos_cell;
  auto sinc = pmy_pack->pcoord->sin_cell;
  auto vol = pmy_pack->pcoord->volume;
  auto &mbsize = pmy_pack->pmb->mb_size;

  Real de_sum = 0.0, de_abs = 0.0, dm_abs = 0.0, ncell = 0.0;
  Real fmin_ = 1.0;
  Kokkos::parallel_reduce("boris_corr", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &s_de, Real &s_dea, Real &s_dma, Real &s_n,
                Real &s_f) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k-ks)*nji)/nx1 + js;
    const int i = (idx - m*nkji - (k-ks)*nji - (j-js)*nx1) + is;
    const Real bx = bcc_(m,IBX,k,j,i), by = bcc_(m,IBY,k,j,i), bz = bcc_(m,IBZ,k,j,i);
    const Real b2 = bx*bx + by*by + bz*bz;
    if (!(b2 > 0.0)) return;
    const Real p = gen ? wder_(m,IDPR,k,j,i) : eos.IdealGasPressure(w0_(m,IEN,k,j,i));
    const Real f = eos.BorisF(eos.BorisX(w0_(m,IDN,k,j,i), b2, p));
    if (f < s_f) s_f = f;
    if (!(f < 1.0)) return;
    if (f < 0.99) s_n += 1.0;
    Real c = 0.0, s = 1.0;
    if (cs) { c = cosc(m,k,j); s = sinc(m,k,j); }
    // full momentum and its reference, orthonormal frame
    const Real m1 = u0_(m,IM1,k,j,i), m2 = u0_(m,IM2,k,j,i), m3 = u0_(m,IM3,k,j,i);
    const Real r1 = mref(m,0,k,j,i), r2 = mref(m,1,k,j,i), r3 = mref(m,2,k,j,i);
    const Real mo[3] = {m1, m2, (m3 - c*m2)/s};
    const Real ro[3] = {r1, r2, (r3 - c*r2)/s};
    const Real bi = 1.0/sqrt(b2);
    const Real bh[3] = {bx*bi, by*bi, bz*bi};
    Real dm[3], mn[3];
    for (int n = 0; n < 3; ++n) dm[n] = mo[n] - ro[n];
    const Real bd = bh[0]*dm[0] + bh[1]*dm[1] + bh[2]*dm[2];
    for (int n = 0; n < 3; ++n) mn[n] = ro[n] + f*dm[n] + (1.0 - f)*bd*bh[n];
    const Real d = u0_(m,IDN,k,j,i);
    const Real dke = 0.5*((mn[0]*mn[0] + mn[1]*mn[1] + mn[2]*mn[2])
                        - (mo[0]*mo[0] + mo[1]*mo[1] + mo[2]*mo[2]))/d;
    u0_(m,IM1,k,j,i) = mn[0];
    u0_(m,IM2,k,j,i) = mn[1];
    u0_(m,IM3,k,j,i) = s*mn[2] + c*mn[1];
    u0_(m,IEN,k,j,i) += dke;
    const Real dv = curv ? vol(m,k,j,i)
                         : mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*mbsize.d_view(m).dx3;
    s_de  += dke*dv;
    s_dea += fabs(dke)*dv;
    s_dma += sqrt(SQR(mn[0]-mo[0]) + SQR(mn[1]-mo[1]) + SQR(mn[2]-mo[2]))*dv;
  }, Kokkos::Sum<Real>(de_sum), Kokkos::Sum<Real>(de_abs), Kokkos::Sum<Real>(dm_abs),
     Kokkos::Sum<Real>(ncell), Kokkos::Min<Real>(fmin_));
  boris_de_sum += de_sum;
  boris_de_abs += de_abs;
  boris_dm_abs += dm_abs;
  if (stage == pdriver->nexp_stages) {
    boris_ncell = ncell;
    boris_fmin = fmin_;
    const int ncyc = pmy_pack->pmesh->ncycle;
    if (global_variable::my_rank == 0 && pmy_pack->gids == 0 &&
        (ncyc % boris_ndiag) == 0) {
      std::cout << "boris: cycle=" << ncyc << " cells(f<0.99)=" << ncell
                << " f_min=" << fmin_ << " cum dE=" << boris_de_sum
                << " cum|dE|=" << boris_de_abs << " cum|dM|=" << boris_dm_abs
                << " (rank 0 pack)" << std::endl;
    }
  }
}

}  // namespace mhd
