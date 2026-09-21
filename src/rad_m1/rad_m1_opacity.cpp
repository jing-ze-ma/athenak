//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_opacity.cpp
//! \brief fill RadiationM1::opac with rho*kappa_{P,E,transport} once per substep stage.
//!
//! This runs BEFORE the fluxes (design sect. 5 task order), because the face value of
//! rho*(kappa_F + kappa_s) is what the thick-limit flux needs, and it covers the ghost
//! cells as well so that the first and last active face have both of their cells.
//!
//! The state is read from hydro's CONSERVED u0 -- density, the kinetic energy that is
//! subtracted to get the internal energy, and from it the temperature -- and never from
//! w0: the coupling writes u0 inside the substep, so w0 is one ConToPrim behind.

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_opacity.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::Opacity

TaskStatus RadiationM1::Opacity(Driver *pdrive, int stage) {
  // with no opacity at all the array stays at the zero it was allocated with, and the
  // thick-limit flux and the coupling are both switched off: nothing to do
  if (opac_zero) return TaskStatus::complete;
  if (pmy_pack->phydro == nullptr) return TaskStatus::complete;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto opac_ = opac;
  auto uh = pmy_pack->phydro->u0;
  auto eos = pmy_pack->phydro->peos->eos_data;
  int otype = opacity_type;
  Real kp = kappa_p, ke = kappa_e, kf = kappa_f, ks = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  bool need_t = (otype != M1_OPAC_CONST);

  par_for("m1_opacity", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real d = uh(m,IDN,k,j,i);
    Real t = 0.0;
    if (need_t) {
      Real ke_dens = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                          SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      Real eint = uh(m,IEN,k,j,i) - ke_dens;
      t = eos.Temperature(d, fmax(eint, 0.0));
    }
    Real op, oe, of, os;
    M1Opacities(otype, d, t, kp, ke, kf, ks, rref, tref, aa, bb, op, oe, of, os);
    opac_(m,M1_OP_P,k,j,i) = d*op;
    opac_(m,M1_OP_E,k,j,i) = d*oe;
    opac_(m,M1_OP_T,k,j,i) = d*(of + os);
  });

  return TaskStatus::complete;
}

} // namespace radm1
