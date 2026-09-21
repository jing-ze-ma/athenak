//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_tests2.cpp
//! \brief the milestone-1d test problems of the grey M1 module, dispatched from
//! rad_m1_tests.cpp when <problem>/m1_test names none of the 1a-1c tests.  Kept in a
//! separate file so that the 1a-1c generators are untouched.
//!
//!   shadow    (design sect. 8, T2) Hayes & Norman (2003) / HERACLES shadow test: a
//!             dense elliptical clump with a smooth (Fermi) edge in a thin ambient
//!             medium, illuminated from the inner-x1 face by a hot source streaming in
//!             +x.  The half domain is modelled (the clump is centred on y = 0 and the
//!             inner-x2 face is reflecting).
//!
//!             UNITS: cgs, except that the gas temperature is carried in KELVIN --
//!             the ideal-gas EOS returns T = (gamma-1) e/rho, so the internal energy is
//!             initialised as rho T/(gamma-1) with T in kelvin and `arad` is then the
//!             cgs radiation constant.  The gas is a PRESCRIBED static background: a
//!             user source term re-imposes (rho, rho v, E_gas) at the end of every
//!             hydro stage and <rad_m1>/gas_feedback = false keeps the coupling from
//!             writing it back, so T = T0 exactly for the whole run and the opacity
//!             sigma = 0.1 (T/T0)^-3.5 (rho/rho0)^2 per cm is a pure function of rho.
//!             (Gas heating on a light-crossing time is negligible in any case: the gas
//!             heat capacity here exceeds a T0^4 by ten orders of magnitude.)

#include <math.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "globals.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "pgen/pgen.hpp"

namespace {
// shadow: the clump geometry and the two background states, needed by the user BC and
// by the source term that re-imposes the gas
Real m1_sh2_rho0 = 1.0, m1_sh2_ratio = 1.0e3, m1_sh2_t0 = 290.0;
Real m1_sh2_x0 = 0.5, m1_sh2_y0 = 0.0, m1_sh2_ax = 0.1, m1_sh2_ay = 0.06;
Real m1_sh2_delta = 10.0, m1_sh2_gm1 = 2.0/3.0;
Real m1_sh2_ein = 1.0, m1_sh2_fin = 1.0 - 1.0e-6, m1_sh2_c = 1.0;

//----------------------------------------------------------------------------------------
//! \fn M1ShadowRho
//! \brief the clump density profile: a Fermi-like smoothing of the ellipse boundary,
//! rho = rho0 [1 + (ratio - 1)/(1 + exp(delta (D - 1)))], D^2 = (dx/ax)^2 + (dy/ay)^2.

KOKKOS_INLINE_FUNCTION
Real M1ShadowRho(const Real x, const Real y, const Real rho0, const Real ratio,
                 const Real x0, const Real y0, const Real ax, const Real ay,
                 const Real delta) {
  Real dd = sqrt(SQR((x - x0)/ax) + SQR((y - y0)/ay));
  Real ex = delta*(dd - 1.0);
  // exp() of a large positive argument overflows; the limit is rho0
  Real w = (ex > 60.0) ? 0.0 : (1.0/(1.0 + exp(ex)));
  return rho0*(1.0 + (ratio - 1.0)*w);
}
} // namespace

// prototypes for the hooks enrolled below
void RadM1ShadowBC(Mesh *pm);
void RadM1ShadowGas(Mesh *pm, const Real bdt);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationM1Tests2()
//! \brief the milestone-1d branches, reached from RadiationM1Tests() when
//! <problem>/m1_test names none of the 1a-1c tests.

void ProblemGenerator::RadiationM1Tests2(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  std::string test = pin->GetOrAddString("problem","m1_test","beam");
  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;
  Real ar = pmbp->pradm1->arad;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->pradm1->u0;

  if (test.compare("shadow") == 0) {
    Real rho0 = pin->GetOrAddReal("problem","shadow_rho0",1.0);
    Real ratio = pin->GetOrAddReal("problem","shadow_ratio",1.0e3);
    Real t0 = pin->GetOrAddReal("problem","shadow_t0",290.0);
    Real tr = pin->GetOrAddReal("problem","shadow_tr",1740.0);
    Real x0 = pin->GetOrAddReal("problem","shadow_x0",0.5);
    Real y0 = pin->GetOrAddReal("problem","shadow_y0",0.0);
    Real ax = pin->GetOrAddReal("problem","shadow_ax",0.1);
    Real ay = pin->GetOrAddReal("problem","shadow_ay",0.06);
    Real delta = pin->GetOrAddReal("problem","shadow_delta",10.0);
    Real fin = pin->GetOrAddReal("problem","shadow_f",1.0 - 1.0e-6);
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    m1_sh2_rho0 = rho0;
    m1_sh2_ratio = ratio;
    m1_sh2_t0 = t0;
    m1_sh2_x0 = x0;
    m1_sh2_y0 = y0;
    m1_sh2_ax = ax;
    m1_sh2_ay = ay;
    m1_sh2_delta = delta;
    m1_sh2_gm1 = gm1;
    m1_sh2_c = cl;
    m1_sh2_fin = fin;
    m1_sh2_ein = ar*tr*tr*tr*tr;
    user_bcs_func = RadM1ShadowBC;
    user_srcs_func = RadM1ShadowGas;
    if (global_variable::my_rank == 0) {
      std::cout << "  m1_test = shadow: E_ambient=" << ar*t0*t0*t0*t0
                << " E_source=" << m1_sh2_ein << " ratio="
                << (m1_sh2_ein/(ar*t0*t0*t0*t0)) << std::endl;
    }
    if (restart) return;
    Real e_amb = fmax(ar*t0*t0*t0*t0, efl);
    auto uh = pmbp->phydro->u0;
    par_for("m1_shadow_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      Real d = M1ShadowRho(x1v, x2v, rho0, ratio, x0, y0, ax, ay, delta);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*t0/gm1;
      u0(m,radm1::M1_E,k,j,i) = e_amb;
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = '" << test << "' not implemented "
      << "(beam | pulse1d | thick_pulse | tophat | jump | equil | advect_pulse "
      << "| advect_uniform | advect_shear | marshak | shadow)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1ShadowGas()
//! \brief re-impose the PRESCRIBED static gas of the shadow test (ghost zones included)
//! at the end of every hydro stage, so that rho and T are exactly the analytic ones and
//! the opacity law is a pure function of position.  bdt is unused.

void RadM1ShadowGas(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto uh = pmbp->phydro->u0;
  Real rho0 = m1_sh2_rho0, ratio = m1_sh2_ratio, t0 = m1_sh2_t0;
  Real x0 = m1_sh2_x0, y0 = m1_sh2_y0, ax = m1_sh2_ax, ay = m1_sh2_ay;
  Real delta = m1_sh2_delta, gm1 = m1_sh2_gm1;

  par_for("m1_shadow_reset", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real d = M1ShadowRho(x1v, x2v, rho0, ratio, x0, y0, ax, ay, delta);
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = d*t0/gm1;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1ShadowBC()
//! \brief inner-x1 illumination of the shadow test: the whole face is a source at
//! T_r = shadow_tr streaming in +x with |f| just below 1.  Every other face is left to
//! the <rad_m1> physical BCs of the input file (reflect at y = 0, vacuum elsewhere).

void RadM1ShadowBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradm1 == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradm1->u0;
  Real cl = m1_sh2_c, ein = m1_sh2_ein, fin = m1_sh2_fin;

  par_for("m1_shadow_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        u0(m,radm1::M1_E,k,j,is-i-1) = ein;
        u0(m,radm1::M1_F1,k,j,is-i-1) = cl*ein*fin;
        u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
        u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
      }
    }
  });
  return;
}
