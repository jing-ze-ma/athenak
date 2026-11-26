//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_opacity.cpp
//! \brief the opacity task: rho kappa of the stage's gas (and |B|^2/2 under MHD).
//!
//!  called from: rad_vet_imex.cpp, rad_vet_newdt.cpp
//!  calls into: rad_vet_newdt.cpp

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "utils/eint_from_cons.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::EmagBuild
//! \brief |B|^2/2 per cell (incl. ghosts) of the field b1 + w (b0 - b1) -> emag0, with the
//! MHD C2P's face average: w = 1 the current field b0, w = g sdirk2's stage-1 field
//! (Time2GasTime). No-op for hydro.

void RadiationVet::EmagBuild(Real w) {
  mhd::MHD *pm = pmy_pack->pmhd;
  if (pm == nullptr) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*(indcs.ng);
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const FaceMix b1f{pm->b0.x1f, pm->b1.x1f, w}, b2f{pm->b0.x2f, pm->b1.x2f, w};
  const FaceMix b3f{pm->b0.x3f, pm->b1.x3f, w};
  auto em = emag0;
  const bool sph = pmy_pack->pmesh->use_spherical_polar;
  const bool csr = pmy_pack->pmesh->use_cubed_sphere;
  auto &x1v = pmy_pack->pcoord->x1v;
  auto &x1f = pmy_pack->pcoord->xx1f;
  auto &x2v = pmy_pack->pcoord->x2v;
  auto &x2f = pmy_pack->pcoord->xx2f;
  auto &x3v = pmy_pack->pcoord->x3v;
  auto &x3f = pmy_pack->pcoord->xx3f;
  if (csr) {
    // cubed sphere: |B|^2 of the orthonormal triple of the MHD C2P (faces hold B.nhat,
    // nhat_xi.nhat_eta = -c): bx = B.rhat, by = (b_xi + c b_eta)/s, bz = b_eta
    auto ccl = pmy_pack->pcoord->cos_cell;
    auto csn = pmy_pack->pcoord->sin_cell;
    par_for("vet_emag_cs", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = ccl(m,k,j), sn = csn(m,k,j);
      const Real bx = CellCenteredRadialFld(b1f(m,k,j,i), b1f(m,k,j,i+1), x1f(m,i),
                                            x1f(m,i+1), x1v(m,i));
      const Real byn = 0.5*(b2f(m,k,j,i) + b2f(m,k,j+1,i));
      const Real bz = 0.5*(b3f(m,k,j,i) + b3f(m,k+1,j,i));
      const Real by = (byn + c*bz)/sn;
      em(m,k,j,i) = 0.5*(bx*bx + by*by + bz*bz);
    });
    return;
  }
  par_for("vet_emag", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    em(m,k,j,i) = VetEmagCell(b1f, b2f, b3f, sph, csr, x1v, x1f, x2v, x2f, x3v, x3f,
                             m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
// opacities

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::Opacity
//! \brief fill opac = rho*kappa (four channels) from the gas state of the stage

TaskStatus RadiationVet::Opacity(Driver *pdrive, int stage) {
  // first gas reader of both VET chains: refresh |B|^2/2 of the stage's field
  if (fl_mhd && !emag_hold) EmagBuild(1.0);
  if (opac_zero) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  if (opacity_type == VET_OPAC_TABLE && otab.nT <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/opacity = table, but no tables were handed to the "
      << "module: the problem generator must call RadiationVet::SetOpacityTables"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto opac_ = opac;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  // with <hydro>/etotgrav the conserved energy includes rho*Phi: subtract it as well
  const bool etg = fl.use_etotgrav;
  auto phicc = fl.phicc0;
  const bool mhd = fl_mhd;      // and |B|^2/2
  auto emag_ = emag0;
  // cubed sphere: the momentum is covariant, KE = 0.5 m_a v^a (gnomonic_raisevel.hpp)
  const bool csk = pmy_pack->pmesh->use_cubed_sphere;
  auto ccl = pmy_pack->pcoord->cos_cell;
  int otype = opacity_type;
  Real kp = kappa_p, ke = kappa_e, kf = kappa_f, ks = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  bool need_t = (otype != VET_OPAC_CONST);
  // copy the tables (held by value) into the closure instead of capturing this
  VetOpacTab ot = otab;
  const bool osc = vscat;   // vet_scatter: rho kappa_e -> opac(VET_OP_S)
  if (ot.scat == 2 && !eos.tbl.active) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/vet_scatter_kappa_e = eos needs a general (tabulated) "
      << "EOS: the electron fraction comes from its Saha composition" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  par_for_lb("vet_opacity", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real d = uh(m,IDN,k,j,i);
    Real t = 0.0;
    if (need_t) {
      Real ke_dens = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                          SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        ke_dens = CsKeMom(uh(m,IM1,k,j,i), m2, m3, c)
                  /fmax(d, 1.0e-300);
      }
      Real eint = uh(m,IEN,k,j,i) - ke_dens;
      if (etg) eint -= d*phicc(m,k,j,i);
      if (mhd) eint -= emag_(m,k,j,i);
      t = eos.Temperature(d, fmax(eint, 0.0));
    }
    Real op, oe, of, os;
    VetCellOpacities(otype, ot, eos, d, t, kp, ke, kf, ks, rref, tref, aa, bb, op, oe, of,
                    os);
    opac_(m,VET_OP_P,k,j,i) = d*op;
    opac_(m,VET_OP_E,k,j,i) = d*oe;
    opac_(m,VET_OP_T,k,j,i) = d*(of + os);
    if (osc) opac_(m,VET_OP_S,k,j,i) = d*os;
  });

  return TaskStatus::complete;
}

}  // namespace radvet
