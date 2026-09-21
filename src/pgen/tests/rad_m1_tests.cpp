//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_tests.cpp
//! \brief test problems for the grey photon M1 module, selected by <problem>/m1_test.
//!
//!   beam    (design sect. 8, T1): 2-D unit box, no opacity, a beam entering at 45
//!           degrees through a patch of the inner-x1 face.  The rest of the domain
//!           starts (and the other faces stay) at the vacuum floor.
//!   pulse1d (convergence check): 1-D periodic Gaussian in E with f = |F|/(cE) = 1
//!           EXACTLY.  At f = 1 the M1 closure gives chi = 1, P = E, the two
//!           eigenvalues merge at +c, and the system reduces to dE/dt + c dE/dx = 0:
//!           the pulse translates rigidly at c, so the exact solution after one box
//!           crossing is the initial condition.  That is the only free-streaming case
//!           that is exactly solvable under M1 (with f < 1 the profile spreads), which
//!           is why f = 1 is used here rather than 0.9.
//!
//! MILESTONE 1b (design sect. 8, T3/T3b/T3c/T5).  All four carry a <hydro> block, which
//! is where rho -- and so rho*kappa -- comes from; the three static ones run with
//! <rad_m1>/gas_feedback = false, which is what "static gas" means here: the radiation
//! feels the opacity, the gas is not pushed or heated by it.
//!
//!   thick_pulse (T3)  1-D Gaussian in E in a scattering-dominated static medium.  The
//!                     variance must grow at 2D, D = c/(3 rho kappa), at every
//!                     tau_cell = rho kappa dx.  `nyquist_amp` seeds the odd-even mode
//!                     on top, whose decay rate is the second half of the gate.
//!   jump        (T3b) constant-flux steady state across a jump in rho*kappa_F (Bloch
//!                     et al. 2021 sect. 5.2).  The gas carries the jump as a density
//!                     step at CONSTANT pressure, so hydro itself stays in equilibrium.
//!                     The analytic steady state is the initial condition and the
//!                     boundaries hold E fixed (user BC below), so the run only has to
//!                     relax the discrete solution.
//!   tophat      (T3c) top-hat in E at tau_cell = 1e4: the limiter-clipped extremum at
//!                     which the uncorrected HLL flux is ~tau_cell times too diffusive.
//!   equil       (T5)  uniform single-zone equilibration, HERACLES/Turner & Stone
//!                     numbers, code units = cgs.
//!
//! MILESTONE 1c (design sect. 8, T4/T4b/T6).  These three carry a MOVING or a HEATING
//! medium, so hydro is evolved for real and gas_feedback is on.
//!
//!   advect_pulse   (T4)  Krumholz et al. (2007) / QUOKKA advecting radiation pulse: a
//!                        Gaussian temperature pulse in radiative equilibrium whose
//!                        density is set by UNIFORM TOTAL (gas + radiation) pressure,
//!                        advected at a uniform v.  Tests the enthalpy-flux split.
//!   advect_uniform (T4b) uniform medium at v, initialised at the EXACT fixed point of
//!                        the sect. 1 sources (F0 = 0, E0 = arad T^4): the gas
//!                        temperature must not move.  Tests the SOURCE form.
//!   marshak        (T6)  cold slab heated by a radiation bath held in the inner-x1
//!                        ghosts; CONSTANT c_v, not Su-Olson (see the input file).

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
// beam parameters, needed by the user boundary function
Real m1_beam_e = 1.0;
Real m1_beam_n1 = 1.0/M_SQRT2;
Real m1_beam_n2 = 1.0/M_SQRT2;
Real m1_beam_f = 1.0 - 1.0e-6;
Real m1_beam_y0 = 0.0;
Real m1_beam_y1 = 0.125;
int  m1_test_id = 0;   // 0 = beam, 1 = pulse1d, 2 = jump
// jump test: the two boundary values of E held fixed in the ghost zones
Real m1_jump_el = 1.0;
Real m1_jump_er = 0.0;
} // namespace

// prototypes for the user BCs
void RadM1BeamBC(Mesh *pm);
void RadM1FixedEBC(Mesh *pm);

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1SetUniformGas
//! \brief fill hydro's conserved u0 with a gas at rest of density d and internal energy
//! density e.  The M1 tests never evolve a velocity, so the conserved state is just
//! (d, 0, 0, 0, e); it is written directly rather than through PrimToCons so that the
//! ghost zones are filled too (the opacity kernel reads them).

void M1SetUniformGas(MeshBlockPack *pmbp, Real d, Real e) {
  if (pmbp->phydro == nullptr) return;
  auto &indcs = pmbp->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto uh = pmbp->phydro->u0;
  par_for("m1_gas_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = e;
  });
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationM1Tests()

void ProblemGenerator::RadiationM1Tests(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pradm1 == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "pgen_name = rad_m1_beam requires a <rad_m1> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string test = pin->GetOrAddString("problem","m1_test","beam");
  if (test.compare("beam") != 0 && test.compare("pulse1d") != 0 &&
      pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = '" << test << "' needs a <hydro> block"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int nx1 = indcs.nx1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->pradm1->u0;

  if (test.compare("beam") == 0) {
    m1_test_id = 0;
    // beam direction and patch, all optional
    Real theta = pin->GetOrAddReal("problem","beam_angle",0.25*M_PI);
    m1_beam_n1 = cos(theta);
    m1_beam_n2 = sin(theta);
    m1_beam_e = pin->GetOrAddReal("problem","beam_e",1.0);
    m1_beam_f = pin->GetOrAddReal("problem","beam_f",1.0-1.0e-6);
    m1_beam_y0 = pin->GetOrAddReal("problem","beam_y0",0.0);
    m1_beam_y1 = pin->GetOrAddReal("problem","beam_y1",0.125);
    user_bcs_func = RadM1BeamBC;

    // vacuum everywhere
    par_for("m1_beam_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radm1::M1_E,k,j,i) = efl;
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("pulse1d") == 0) {
    m1_test_id = 1;
    Real amp = pin->GetOrAddReal("problem","pulse_amp",1.0);
    Real wid = pin->GetOrAddReal("problem","pulse_width",0.1);
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.5);
    Real bg = pin->GetOrAddReal("problem","pulse_bg",1.0e-4);

    par_for("m1_pulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real e = bg + amp*exp(-SQR((x1v - x0)/wid));
      u0(m,radm1::M1_E,k,j,i) = e;
      u0(m,radm1::M1_F1,k,j,i) = cl*e;   // f = 1 exactly: rigid translation at +c
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("thick_pulse") == 0 || test.compare("tophat") == 0) {
    // T3 / T3c.  Static uniform gas; the opacity is whatever <rad_m1> says, and with
    // kappa_s only the gas neither emits nor absorbs.
    bool tophat = (test.compare("tophat") == 0);
    m1_test_id = tophat ? 3 : 2;
    Real amp = pin->GetOrAddReal("problem","pulse_amp",1.0);
    Real wid = pin->GetOrAddReal("problem","pulse_width",0.05);
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.5);
    Real bg = pin->GetOrAddReal("problem","pulse_bg",1.0e-3);
    Real nyq = pin->GetOrAddReal("problem","nyquist_amp",0.0);
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real egas = pin->GetOrAddReal("problem","gas_eint",1.0);
    M1SetUniformGas(pmbp, dgas, egas);

    // the Nyquist seed is (-1)^(global cell index), so it must not restart at each
    // MeshBlock: the index is rebuilt from the block's own x1min
    Real mx1min = pmy_mesh_->mesh_size.x1min;
    Real dxm = (pmy_mesh_->mesh_size.x1max - mx1min)/
               static_cast<Real>(pmy_mesh_->mesh_indcs.nx1);
    par_for("m1_thickpulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real e = bg;
      if (tophat) {
        e += (fabs(x1v - x0) < wid) ? amp : 0.0;
      } else {
        e += amp*exp(-SQR((x1v - x0)/wid));
      }
      if (nyq != 0.0) {
        int ig = static_cast<int>(floor((x1v - mx1min)/dxm));
        e += ((ig % 2) == 0) ? nyq : -nyq;
      }
      u0(m,radm1::M1_E,k,j,i) = fmax(e, efl);
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("jump") == 0) {
    // T3b.  rho jumps by jump_ratio at x = x_jump at CONSTANT gas pressure; with
    // opacity = const that is a jump of the same ratio in rho*kappa_F.  E is initialised
    // with the analytic constant-flux steady state, piecewise linear with slopes in the
    // ratio of the two opacities, and its two end values are held in the ghost zones.
    m1_test_id = 4;
    Real xj = pin->GetOrAddReal("problem","x_jump",0.5);
    Real ratio = pin->GetOrAddReal("problem","jump_ratio",1.0e3);
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real pgas = pin->GetOrAddReal("problem","gas_pres",1.0);
    Real e_l = pin->GetOrAddReal("problem","e_left",1.0);
    Real flux = pin->GetOrAddReal("problem","jump_flux",1.0e-3);
    Real kf = pmbp->pradm1->kappa_f + pmbp->pradm1->kappa_s;
    Real x1l = pmy_mesh_->mesh_size.x1min;
    Real x1r = pmy_mesh_->mesh_size.x1max;
    // dE/dx = -3 rho kappa F/c on each side of the jump
    Real s_l = -3.0*dgas*kf*flux/cl;
    Real s_r = ratio*s_l;
    Real ej = e_l + s_l*(xj - x1l);
    m1_jump_el = e_l;
    m1_jump_er = ej + s_r*(x1r - xj);
    user_bcs_func = RadM1FixedEBC;

    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real eint = pgas/gm1;
    auto uh = pmbp->phydro->u0;
    par_for("m1_jump_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real d = (x1v < xj) ? dgas : (ratio*dgas);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = eint;
      Real e = (x1v < xj) ? (e_l + s_l*(x1v - x1l)) : (ej + s_r*(x1v - xj));
      u0(m,radm1::M1_E,k,j,i) = fmax(e, efl);
      u0(m,radm1::M1_F1,k,j,i) = flux;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("advect_pulse") == 0) {
    // T4.  Krumholz, Klein & McKee (2007) Sect. 4.3 / QUOKKA AP paper (2404.08247)
    // Sect. 5 advecting radiation pulse, in code units with T = P/rho (so the cgs
    // temperature is T_K = T_code mu m_H/k and arad must be given in the same units).
    //
    //   T(x)   = T0 + (T1 - T0) exp(-(x - x0)^2/(2 w^2))
    //   rho(x) = [P_tot - arad T^4/3]/T,   P_tot = rho0 T0 + arad T0^4/3
    //
    // i.e. the density is adjusted so that the TOTAL (gas + radiation) pressure is
    // uniform, which is what Krumholz et al. prescribe; a temperature pulse at uniform
    // GAS pressure alone would be pushed apart by the radiation pressure gradient.
    // The gas is in radiative equilibrium with the pulse, E = arad T^4, and the flux is
    // the dynamic-diffusion one, F = (4/3) v E - (c/(3 rho kappa)) dE/dx, evaluated
    // analytically.  v is UNIFORM and the gas really is advected by hydro.
    m1_test_id = 6;
    Real t0 = pin->GetReal("problem","pulse_t0");
    Real t1 = pin->GetReal("problem","pulse_t1");
    Real d0 = pin->GetReal("problem","gas_rho");
    Real wid = pin->GetReal("problem","pulse_width");
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.0);
    Real vx = pin->GetOrAddReal("problem","pulse_v",0.0);
    Real ar = pmbp->pradm1->arad;
    Real kap = pmbp->pradm1->kappa_f + pmbp->pradm1->kappa_s;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real ptot = d0*t0 + ar*t0*t0*t0*t0/3.0;
    auto uh = pmbp->phydro->u0;
    par_for("m1_advpulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real s = (x1v - x0)/wid;
      Real g = exp(-0.5*s*s);
      Real tt = t0 + (t1 - t0)*g;
      Real t4 = tt*tt*tt*tt;
      Real d = (ptot - ar*t4/3.0)/tt;
      Real e = ar*t4;
      // dE/dx = 4 arad T^3 dT/dx, dT/dx = (T1 - T0) g (-(x-x0)/w^2)
      Real dtdx = -(t1 - t0)*g*(x1v - x0)/(wid*wid);
      Real dedx = 4.0*ar*tt*tt*tt*dtdx;
      Real ff = (4.0/3.0)*vx*e - cl*dedx/(3.0*d*kap);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = d*vx;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*vx*vx;
      u0(m,radm1::M1_E,k,j,i) = fmax(e, efl);
      u0(m,radm1::M1_F1,k,j,i) = ff;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("advect_uniform") == 0) {
    // T4b (QUOKKA AP paper Sect. 5.5).  Uniform medium moving at v, initialised at the
    // EXACT fixed point of the design sect. 1 source terms, so that a correct scheme
    // moves nothing at all: F0_i = 0 and E0 = arad T^4.  In 1-D, with f = F/(cE),
    //   F0/c = 0  <=>  f = beta (1 + chi(f))            (f -> (4/3) beta as chi -> 1/3)
    //   E0     = E [(1 + beta^2) - 2 beta f + beta^2 chi] = arad T^4
    // so E = arad T^4/[...] = arad T^4 (1 + (4/3) beta^2 + ...) and F = c E f.  Setting
    // E = arad T^4 instead (the naive reading of "F = (4/3) v E") leaves a REAL O(beta^2)
    // relaxation that has nothing to do with the scheme, which is why the fixed point is
    // solved for here.  The gas temperature must then not move at all.
    m1_test_id = 7;
    Real tg = pin->GetReal("problem","gas_temp");
    Real d0 = pin->GetReal("problem","gas_rho");
    Real vx = pin->GetOrAddReal("problem","pulse_v",0.0);
    Real ar = pmbp->pradm1->arad;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real bb = vx/cl;
    // fixed point of f = beta (1 + chi(f)), a handful of Picard steps (contraction
    // factor ~ beta)
    Real ff = (4.0/3.0)*bb;
    Real chi = 1.0/3.0;
    for (int it=0; it<50; ++it) {
      chi = radm1::M1Chi(fabs(ff));
      ff = bb*(1.0 + chi);
    }
    Real erad = ar*tg*tg*tg*tg/((1.0 + bb*bb) - 2.0*bb*ff + bb*bb*chi);
    Real frad = cl*erad*ff;
    if (global_variable::my_rank == 0) {
      std::cout << "  m1_test = advect_uniform: beta=" << bb << " f=" << ff
                << " E=" << erad << " F=" << frad << " F/((4/3)vE)="
                << frad/((4.0/3.0)*vx*erad) << std::endl;
    }
    auto uh = pmbp->phydro->u0;
    par_for("m1_advunif_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      uh(m,IDN,k,j,i) = d0;
      uh(m,IM1,k,j,i) = d0*vx;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d0*tg/gm1 + 0.5*d0*vx*vx;
      u0(m,radm1::M1_E,k,j,i) = erad;
      u0(m,radm1::M1_F1,k,j,i) = frad;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("marshak") == 0) {
    // T6.  Non-equilibrium Marshak wave: a cold, static, uniform slab of constant
    // opacity heated from x1min by an incident radiation bath.  The material heat
    // capacity is the code's CONSTANT c_v, not the Su-Olson alpha T^3 (which the EOS
    // cannot represent); tests_m1/t6_marshak.py's own S_N reference is run with the same
    // constant c_v (--cv), so the comparison is like for like.  The gas does not move
    // (gas_feedback still writes the energy exchange; the momentum exchange is left in,
    // and stays ~1e-10 of the pressure over the run).
    m1_test_id = 8;
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real tini = pin->GetReal("problem","gas_temp");
    Real ebath = pin->GetReal("problem","e_bath");
    Real ar = pmbp->pradm1->arad;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    m1_jump_el = ebath;              // held in the inner-x1 ghost zones
    m1_jump_er = ar*tini*tini*tini*tini;
    user_bcs_func = RadM1FixedEBC;
    Real e_ini = fmax(ar*tini*tini*tini*tini, efl);
    M1SetUniformGas(pmbp, dgas, dgas*tini/gm1);
    par_for("m1_marshak_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radm1::M1_E,k,j,i) = e_ini;
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("equil") == 0) {
    // T5.  Uniform single zone: nothing has a gradient, so the transport is exactly
    // zero and what is exercised is the implicit solve alone.
    m1_test_id = 5;
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0e-7);
    Real egas = pin->GetReal("problem","gas_eint");
    Real erad = pin->GetReal("problem","e_rad");
    M1SetUniformGas(pmbp, dgas, egas);
    par_for("m1_equil_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radm1::M1_E,k,j,i) = erad;
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = '" << test << "' not implemented "
      << "(beam | pulse1d | thick_pulse | tophat | jump | equil | advect_pulse "
      << "| advect_uniform | marshak)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1FixedEBC()
//! \brief Dirichlet in E, zero-gradient in F, on both x1 faces: the steady state of T3b
//! is set by the two end values of E, while F has to be free for the solution to find
//! its own constant flux.

void RadM1FixedEBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradm1 == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradm1->u0;
  Real el = m1_jump_el, er = m1_jump_er;

  par_for("m1_fixede_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        u0(m,radm1::M1_E,k,j,is-i-1) = el;
        u0(m,radm1::M1_F1,k,j,is-i-1) = u0(m,radm1::M1_F1,k,j,is);
        u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
        u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        u0(m,radm1::M1_E,k,j,ie+i+1) = er;
        u0(m,radm1::M1_F1,k,j,ie+i+1) = u0(m,radm1::M1_F1,k,j,ie);
        u0(m,radm1::M1_F2,k,j,ie+i+1) = 0.0;
        u0(m,radm1::M1_F3,k,j,ie+i+1) = 0.0;
      }
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1BeamBC()
//! \brief inner-x1 inflow PATCH for the beam test: inside [beam_y0, beam_y1] the ghost
//! cells hold E = beam_e and F = c E f n (f just below 1, so the closure is at chi ~ 1
//! without ever dividing by a vanishing sqrt(4-3f^2)); outside the patch, and on the
//! inner-x2 face, they hold the vacuum floor so nothing enters.  Outer faces are left to
//! the <rad_m1> physical BCs (vacuum) set in the input file.

void RadM1BeamBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradm1 == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int nx2 = indcs.nx2;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradm1->u0;

  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;
  Real be = m1_beam_e, bf = m1_beam_f, bn1 = m1_beam_n1, bn2 = m1_beam_n2;
  Real by0 = m1_beam_y0, by1 = m1_beam_y1;

  par_for("m1_beam_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      bool inbeam = (x2v >= by0 && x2v <= by1);
      for (int i=0; i<ng; ++i) {
        if (inbeam) {
          u0(m,radm1::M1_E,k,j,is-i-1) = be;
          u0(m,radm1::M1_F1,k,j,is-i-1) = cl*be*bf*bn1;
          u0(m,radm1::M1_F2,k,j,is-i-1) = cl*be*bf*bn2;
          u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
        } else {
          u0(m,radm1::M1_E,k,j,is-i-1) = efl;
          u0(m,radm1::M1_F1,k,j,is-i-1) = 0.0;
          u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
          u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
        }
      }
    }
  });
  return;
}
