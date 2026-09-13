//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_diff2d.cpp
//! \brief TRANSVERSE radiative diffusion of a 2D Gaussian: the test problem of the
//! implicit transverse operator (<hydro>/rad_implicit_ang, conduction_transverse.cpp).
//!
//! A Cartesian box at rest, uniform in x1, with a small Gaussian temperature
//! perturbation in (x2,x3) on a uniform density.  For an amplitude of 1e-4 the radiative
//! conductivity kappa_rad = 16 sigma T^3/(3 kappa_R rho) is constant to a few parts in
//! 1e4, so the perturbation obeys the linear heat equation dT/dt = D grad^2 T with
//! D = kappa_rad/(rho c_v) and the exact solution is the spreading Gaussian
//!     dT(r,t) = amp T0 (s0^2/s^2) exp(-r^2/(2 s^2)),   s^2 = s0^2 + 2 D t.
//! Run in kinematic mode (rsolver = advect, v = 0) so that conduction is the only
//! operator that moves anything and the comparison is with the analytic solution rather
//! than with a converged run.  The domain is periodic in x2/x3; s0 is small enough that
//! the periodic images are below 1e-7 of the amplitude at the final time.
//!
//! D is computed here from the SAME Freedman Rosseland mean and the same unit
//! conversion the conduction module uses, so the reference solution carries no
//! independently-tuned constant.  It is printed at start-up together with the time
//! 2 D t = s0^2 at which the variance has doubled, which is what <time>/tlim is set to.
//!
//! The errors against that solution are written to <basename>-errs.dat by the standard
//! OutputErrors path, exactly as the linear-wave convergence tests do.

// C++ headers
#include <cmath>
#include <iostream>
#include <string>

// AthenaK headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "utils/rosseland.hpp"

void RadDiff2DErrors(ParameterInput *pin, Mesh *pm);

namespace {
// true while setting the initial conditions, false while filling the reference register
bool set_initial_conditions = true;
struct RadDiff2DVars {
  Real d0, p0, amp, sig0, x20, x30, dcoeff;
};
RadDiff2DVars rdv;
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::RadDiff2D()
//! \brief sets the initial conditions, or (on the second call, from RadDiff2DErrors) the
//! analytic solution at the current time into the u1 register

void ProblemGenerator::RadDiff2D(ParameterInput *pin, const bool restart) {
  pgen_final_func = RadDiff2DErrors;
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr || pmbp->phydro->pcond == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "rad_diff2d needs a <hydro> block with isotropic_conduction"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  if (!(eos.is_ideal)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "rad_diff2d needs an ideal EOS" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real gm1 = eos.gamma - 1.0;

  if (set_initial_conditions) {
    rdv.d0 = pin->GetOrAddReal("problem", "d0", 1.0);
    rdv.p0 = pin->GetOrAddReal("problem", "p0", 1.0);
    rdv.amp = pin->GetOrAddReal("problem", "amp", 1.0e-4);
    rdv.sig0 = pin->GetOrAddReal("problem", "sig0", 0.06);
    rdv.x20 = pin->GetOrAddReal("problem", "x20", 0.0);
    rdv.x30 = pin->GetOrAddReal("problem", "x30", 0.0);

    // D = kappa_rad/(rho c_v) in code units, from the background state.  The two lines
    // below are RadiativeKappa() and the code-unit conversion of conduction.cpp,
    // written out here because both live inside that translation unit.
    Conduction *pc = pmbp->phydro->pcond;
    const Real temp_unit = pmbp->punit->temperature_cgs();
    const Real pres_unit = pmbp->punit->pressure_cgs();
    const Real dens_unit = pmbp->punit->density_cgs();
    const Real len_unit = pmbp->punit->length_cgs();
    const Real eflx_unit = pres_unit*pmbp->punit->velocity_cgs();
    const Real t0 = rdv.p0/rdv.d0;
    const Real tk = t0*temp_unit;
    const Real sigma_sb = 5.670374419e-5;
    const Real kr = pc->rad_kappa_fac
                    *RosselandFreedman2014(tk, rdv.p0*pres_unit, pc->rad_met);
    const Real kap = 16.0*sigma_sb*tk*tk*tk/(3.0*kr*rdv.d0*dens_unit);
    const Real kcode = kap*temp_unit/len_unit/eflx_unit;
    rdv.dcoeff = kcode*gm1/rdv.d0;
    if (global_variable::my_rank == 0) {
      std::cout << "rad_diff2d: kappa_R = " << kr << " cm^2/g, kappa_rad = " << kap
                << " erg/cm/s/K, D = " << rdv.dcoeff << " (code), s0^2/(2 D) = "
                << (rdv.sig0*rdv.sig0/(2.0*rdv.dcoeff)) << " (code time)" << std::endl;
    }
  }

  // capture for the kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  const Real d0 = rdv.d0, t0 = rdv.p0/rdv.d0, amp = rdv.amp, gmm1 = gm1;
  const Real x20 = rdv.x20, x30 = rdv.x30;
  // the domain is periodic in x2/x3, and the periodic solution of the heat equation is
  // the sum over the images of the infinite-domain Gaussian.  Two images each way is
  // exact to well below round-off for every width this test uses.
  const Real l2 = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  const Real l3 = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  const Real tnow = set_initial_conditions ? 0.0 : pmbp->pmesh->time;
  const Real s2 = rdv.sig0*rdv.sig0 + 2.0*rdv.dcoeff*tnow;
  const Real fac = amp*t0*(rdv.sig0*rdv.sig0)/s2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;

  auto &u1 = (set_initial_conditions) ? pmbp->phydro->u0 : pmbp->phydro->u1;
  par_for("pgen_raddiff2d", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je,
  is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    Real g = 0.0;
    for (int n2=-2; n2<=2; ++n2) {
      for (int n3=-2; n3<=2; ++n3) {
        const Real r2 = SQR(x2v - x20 + n2*l2) + SQR(x3v - x30 + n3*l3);
        g += exp(-0.5*r2/s2);
      }
    }
    const Real temp = t0 + fac*g;
    u1(m,IDN,k,j,i) = d0;
    u1(m,IM1,k,j,i) = 0.0;
    u1(m,IM2,k,j,i) = 0.0;
    u1(m,IM3,k,j,i) = 0.0;
    u1(m,IEN,k,j,i) = d0*temp/gmm1;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadDiff2DErrors()
//! \brief fills u1 with the analytic solution at the final time and writes the errors

void RadDiff2DErrors(ParameterInput *pin, Mesh *pm) {
  set_initial_conditions = false;
  pm->pgen->RadDiff2D(pin, false);
  pm->pgen->OutputErrors(pin, pm);
  return;
}
