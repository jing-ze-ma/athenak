//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_atm.cpp
//! \brief 2D polytropic atmosphere in constant gravity, for the STATIC well-balanced
//! scheme.
//!
//! The stratification is the isentropic (polytropic) hydrostatic solution
//!     p(y) = p0 (1 - y/H0)^(gamma/(gamma-1)),   rho(y) = (p/A0)^(1/gamma),
//!     A0   = p0/rho0^gamma,   H0 = -[gamma/(gamma-1)] p0/(rho0 g),
//! with g = <problem>/g0 < 0 pointing along -x2.  The same profile is handed to the
//! hydro module as the static well-balanced background (u0wb, w0wb, w0facewb), so with
//!   <hydro>/wellbalance_static = true, wellbalance_static_reconst = true
//! the deviation reconstruction sees an identically zero perturbation and the initial
//! state is a machine-exact fixed point.  With those flags false the very same problem
//! runs as an ordinary stratified atmosphere, which is what makes a WB-on / WB-off
//! comparison meaningful.
//!
//! <problem>/vamp adds a strong radially diverging velocity kick inside a circle of
//! radius <problem>/r0 about (x0, y0).  At Mach >~ 10 the kick opens a near-vacuum, so
//! the trial conserved state genuinely needs the density/pressure floors -- which is
//! what first-order flux correction (FOFC) exists for.  vamp = 0 leaves the balanced
//! state untouched.
//!
//! The x2 boundaries must be set to "user": the ghost zones are filled with the same
//! analytic background, so that the balance holds right up to the physical boundary.
//! Gravity is applied through <problem>/user_srcs = true (WbAtmGravity below), which is
//! the only place that knows the source has to act on rho - rho_background when the
//! static well-balanced scheme is on.
//!
//! The same problem also drives the DYNAMIC well-balanced scheme (Kappeli & Mishra,
//!   <hydro>/wellbalance_dynamic = true, wb_x2 = true, wb_rho = true,
//!   wb_option = isentropic),
//! which does NOT use u0wb/w0wb at all: it rebuilds a local hydrostatic background from
//! the current state and the gravitational potential.  What it needs from the pgen is
//! that potential -- phicc0 and phi0.x{1,2,3}f, filled over the full arrays below -- and
//! a gravitational source written as the background's own pressure difference (see
//! WbAtmGravity).  The atmosphere here is isentropic, which is exactly the closed-form
//! background wb_option = isentropic assumes, so the balance is held to round-off.

// C++ headers
#include <cmath>
#include <iostream>

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen/pgen.hpp"

namespace {
// parameters of the background, shared by the pgen, the BCs and the source term
struct WbAtmVars {
  Real g0, p0, rho0, gamma, a0, h0;
};
WbAtmVars wbatm;

//! \fn WbAtmPres
//! \brief gas pressure of the hydrostatic background at height y
KOKKOS_INLINE_FUNCTION
Real WbAtmPres(const Real y, const Real p0, const Real h0, const Real gdgm1) {
  return p0*std::pow((1.0 - y/h0), gdgm1);
}
}  // namespace

void WbAtmBCs(Mesh *pm);
void WbAtmGravity(Mesh *pm, const Real bdt);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::WbAtmosphere()
//! \brief sets up the well-balanced atmosphere and its background arrays

void ProblemGenerator::WbAtmosphere(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "wb_atm problem generator requires <hydro>" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(pmy_mesh_->two_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "wb_atm problem generator only works in 2D" << std::endl;
    exit(EXIT_FAILURE);
  }

  Real g0 = pin->GetOrAddReal("problem", "g0", -1.0);
  Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real gamma = pmbp->phydro->peos->eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real gdgm1 = gamma/gm1;
  wbatm.g0 = g0;
  wbatm.p0 = p0;
  wbatm.rho0 = rho0;
  wbatm.gamma = gamma;
  wbatm.a0 = p0/std::pow(rho0, gamma);
  wbatm.h0 = -gdgm1*p0/rho0/g0;
  Real a0 = wbatm.a0;
  Real h0 = wbatm.h0;

  // the source term must know about the background, and the BCs about both
  user_srcs_func = WbAtmGravity;
  user_bcs_func = WbAtmBCs;
  if (restart) return;

  Real vamp = pin->GetOrAddReal("problem", "vamp", 0.0);
  Real x0 = pin->GetOrAddReal("problem", "x0", 0.0);
  Real y0 = pin->GetOrAddReal("problem", "y0", 1.0);
  Real r0 = pin->GetOrAddReal("problem", "r0", 0.2);

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = indcs.nx2 + 2*ng - 1;
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->phydro->u0;
  const bool wbstatic = pmbp->phydro->use_wellbalance_static;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;

  // initial conditions: background + (optional) velocity kick
  par_for("wb_atm_ic", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    Real p = WbAtmPres(x2v, p0, h0, gdgm1);
    Real den = std::pow(p/a0, 1.0/gamma);

    Real v1 = 0.0, v2 = 0.0;
    Real rad = std::sqrt(SQR(x1v - x0) + SQR(x2v - y0));
    if (rad < r0) {
      v1 = vamp*(x1v - x0)/r0;
      v2 = vamp*(x2v - y0)/r0;
    }
    u0(m,IDN,k,j,i) = den;
    u0(m,IM1,k,j,i) = den*v1;
    u0(m,IM2,k,j,i) = den*v2;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = p/gm1 + 0.5*den*(v1*v1 + v2*v2);
  });

  // ---- DYNAMIC well-balanced scheme (Kappeli & Mishra) -------------------------------
  // It builds its own local hydrostatic background out of the CURRENT state and the
  // gravitational potential, so it needs no u0wb/w0wb; what it does need is the
  // potential, cell-centred (phicc0) and on the faces (phi0), over the FULL arrays
  // including the ghost zones -- the x2 reconstruction stencil reaches one cell beyond
  // the active range, and the source term below differences the x2 faces.  With
  // g = g0 = const along x2 the potential is simply Phi = -g0*y.
  if (wbdyn) {
    auto phicc0 = pmbp->phydro->phicc0;
    auto phi1 = pmbp->phydro->phi0.x1f;
    auto phi2 = pmbp->phydro->phi0.x2f;
    auto phi3 = pmbp->phydro->phi0.x3f;
    par_for("wb_atm_phi", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, 0, 0, n2m1,
    0, n1m1, KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
      Real phicc = -g0*x2v;
      phicc0(m,k,j,i) = phicc;
      // the x1- and x3-faces sit at the cell-centred height
      phi1(m,k,j,i) = phicc;
      if (i == n1m1) { phi1(m,k,j,i+1) = phicc; }
      phi3(m,k,j,i) = phicc;
      phi3(m,k+1,j,i) = phicc;
      // the x2-face carries its own (lower-edge) height
      phi2(m,k,j,i) = -g0*LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
      if (j == n2m1) {
        phi2(m,k,j+1,i) = -g0*LeftEdgeX(j+1-js, indcs.nx2, x2min, x2max);
      }
    });
  }

  if (!wbstatic) return;

  // static well-balanced background, over the FULL arrays including ghost zones
  auto u0wb = pmbp->phydro->u0wb;
  auto w0wb = pmbp->phydro->w0wb;
  auto wf1 = pmbp->phydro->w0facewb.x1f;
  auto wf2 = pmbp->phydro->w0facewb.x2f;
  auto wf3 = pmbp->phydro->w0facewb.x3f;
  par_for("wb_atm_bg", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, 0, 0, n2m1,
  0, n1m1, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real p = WbAtmPres(x2v, p0, h0, gdgm1);
    Real den = std::pow(p/a0, 1.0/gamma);

    u0wb(m,IDN,k,j,i) = den;
    u0wb(m,IM1,k,j,i) = 0.0;
    u0wb(m,IM2,k,j,i) = 0.0;
    u0wb(m,IM3,k,j,i) = 0.0;
    u0wb(m,IEN,k,j,i) = p/gm1;
    w0wb(m,IDN,k,j,i) = den;
    w0wb(m,IM1,k,j,i) = 0.0;
    w0wb(m,IM2,k,j,i) = 0.0;
    w0wb(m,IM3,k,j,i) = 0.0;
    w0wb(m,IEN,k,j,i) = p/gm1;

    // x1- and x3-faces carry the CELL-CENTRED height, the x2-face its own left edge
    wf1(m,IDN,k,j,i) = den;
    wf1(m,IM1,k,j,i) = 0.0;
    wf1(m,IM2,k,j,i) = 0.0;
    wf1(m,IM3,k,j,i) = 0.0;
    wf1(m,IEN,k,j,i) = p/gm1;
    wf3(m,IDN,k,j,i) = den;
    wf3(m,IM1,k,j,i) = 0.0;
    wf3(m,IM2,k,j,i) = 0.0;
    wf3(m,IM3,k,j,i) = 0.0;
    wf3(m,IEN,k,j,i) = p/gm1;
    wf3(m,IDN,k+1,j,i) = den;
    wf3(m,IM1,k+1,j,i) = 0.0;
    wf3(m,IM2,k+1,j,i) = 0.0;
    wf3(m,IM3,k+1,j,i) = 0.0;
    wf3(m,IEN,k+1,j,i) = p/gm1;
    if (i == n1m1) {
      wf1(m,IDN,k,j,i+1) = den;
      wf1(m,IM1,k,j,i+1) = 0.0;
      wf1(m,IM2,k,j,i+1) = 0.0;
      wf1(m,IM3,k,j,i+1) = 0.0;
      wf1(m,IEN,k,j,i+1) = p/gm1;
    }

    Real x2f = LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
    Real pf = WbAtmPres(x2f, p0, h0, gdgm1);
    Real denf = std::pow(pf/a0, 1.0/gamma);
    wf2(m,IDN,k,j,i) = denf;
    wf2(m,IM1,k,j,i) = 0.0;
    wf2(m,IM2,k,j,i) = 0.0;
    wf2(m,IM3,k,j,i) = 0.0;
    wf2(m,IEN,k,j,i) = pf/gm1;
    if (j == n2m1) {
      x2f = LeftEdgeX(j+1-js, indcs.nx2, x2min, x2max);
      pf = WbAtmPres(x2f, p0, h0, gdgm1);
      denf = std::pow(pf/a0, 1.0/gamma);
      wf2(m,IDN,k,j+1,i) = denf;
      wf2(m,IM1,k,j+1,i) = 0.0;
      wf2(m,IM2,k,j+1,i) = 0.0;
      wf2(m,IM3,k,j+1,i) = 0.0;
      wf2(m,IEN,k,j+1,i) = pf/gm1;
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WbAtmBCs
//! \brief x2 ghost zones held at the analytic hydrostatic background, so that the
//! balance is exact right up to the physical boundary.

void WbAtmBCs(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int &js = indcs.js; int &je = indcs.je;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto u0 = pmbp->phydro->u0;

  Real p0 = wbatm.p0, h0 = wbatm.h0, a0 = wbatm.a0, gamma = wbatm.gamma;
  Real gm1 = gamma - 1.0;
  Real gdgm1 = gamma/gm1;

  par_for("wb_atm_bc", DevExeSpace(), 0, nmb1, 0, 0, 0, (ng-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
      Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
      Real p = WbAtmPres(x2v, p0, h0, gdgm1);
      Real den = std::pow(p/a0, 1.0/gamma);
      u0(m,IDN,k,j,i) = den;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = p/gm1;
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
      Real x2v = CellCenterX((je+j+1)-js, indcs.nx2, x2min, x2max);
      Real p = WbAtmPres(x2v, p0, h0, gdgm1);
      Real den = std::pow(p/a0, 1.0/gamma);
      u0(m,IDN,k,(je+j+1),i) = den;
      u0(m,IM1,k,(je+j+1),i) = 0.0;
      u0(m,IM2,k,(je+j+1),i) = 0.0;
      u0(m,IM3,k,(je+j+1),i) = 0.0;
      u0(m,IEN,k,(je+j+1),i) = p/gm1;
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WbAtmGravity
//! \brief constant gravity along -x2.  With the static well-balanced scheme on, the
//! background's own weight is already carried by the flux (RemoveWbFlux), so the
//! momentum source must act on the DEVIATION rho - rho_background only.  With the
//! DYNAMIC scheme on there is no stored background to subtract: the momentum source is
//! instead the local hydrostatic background's OWN pressure difference across the cell,
//! (p_{j+1/2} - p_{j-1/2})/dx2, obtained through the same entry point the
//! reconstruction uses (Hydro::getWBq0) so that the two cancel identically for a
//! balanced state.  (Handing the dynamic scheme the plain -rho g source would NOT be
//! well balanced, since the flux it assembles already carries the full background
//! pressure.)  The energy source is physical and always uses the full density.

void WbAtmGravity(Mesh *pm, const Real bdt) {
  auto &indcs = pm->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pm->pmb_pack;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto u0 = pmbp->phydro->u0;
  auto w0 = pmbp->phydro->w0;
  auto w0wb = pmbp->phydro->w0wb;
  auto phicc0 = pmbp->phydro->phicc0;
  auto phi2 = pmbp->phydro->phi0.x2f;
  auto &size = pmbp->pmb->mb_size;
  const bool wbstatic = pmbp->phydro->use_wellbalance_static;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;
  // by-value copies, capturable in the device lambda (never dereference pmbp there)
  const auto eos = pmbp->phydro->peos->eos_data;
  const WBOption wbo = pmbp->phydro->wb_option;
  Real g0 = wbatm.g0;

  par_for("wb_atm_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real src = bdt*g0*w0(m,IDN,k,j,i);
    u0(m,IEN,k,j,i) += src*w0(m,IVY,k,j,i);
    if (wbstatic) {
      src = bdt*g0*(w0(m,IDN,k,j,i) - w0wb(m,IDN,k,j,i));
    }
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
          w0(m,IDN,k,j-1,i), w0(m,IDN,k,j,i), w0(m,IDN,k,j+1,i),
          w0(m,IEN,k,j-1,i), w0(m,IEN,k,j,i), w0(m,IEN,k,j+1,i),
          phicc0(m,k,j-1,i), phi2(m,k,j,i), phicc0(m,k,j,i), phi2(m,k,j+1,i),
          phicc0(m,k,j+1,i), d1, pl, d2, pr, d3);
      src = bdt*(pr - pl)/size.d_view(m).dx2;
    }
    u0(m,IM2,k,j,i) += src;
  });
  return;
}
