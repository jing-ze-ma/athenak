//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tests.cpp
//! \brief test problems for the grey photon VET module, selected by <problem>/vet_test.
//!
//!   beam    (design sect. 8, T1): 2-D unit box, no opacity, at the vacuum floor; with
//!           beam_packet a free-streaming Gaussian packet travelling at beam_angle.
//!   pulse1d (convergence check): 1-D periodic Gaussian in E with f = |F|/(cE) = 1
//!           EXACTLY.  At f = 1 the M1 closure gives chi = 1, P = E, the two
//!           eigenvalues merge at +c, and the system reduces to dE/dt + c dE/dx = 0:
//!           the pulse translates rigidly at c, so the exact solution after one box
//!           crossing is the initial condition.  That is the only free-streaming case
//!           that is exactly solvable under M1 (with f < 1 the profile spreads), which
//!           is why f = 1 is used here rather than 0.9.
//!
//! MILESTONE 1b (design sect. 8, T3/T3c/T5).  All three carry a <hydro> block, which
//! is where rho -- and so rho*kappa -- comes from; the two static ones run with
//! <rad_vet>/gas_feedback = false, which is what "static gas" means here: the radiation
//! feels the opacity, the gas is not pushed or heated by it.
//!
//!   thick_pulse (T3)  1-D Gaussian in E in a scattering-dominated static medium.  The
//!                     variance must grow at 2D, D = c/(3 rho kappa), at every
//!                     tau_cell = rho kappa dx.  `nyquist_amp` seeds the odd-even mode
//!                     on top, whose decay rate is the second half of the gate.
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
//!   marshak        (T6)  cold slab heated by a radiation bath at x1min (the <rad_vet>
//!                        x1 boundary); CONSTANT c_v, not Su-Olson (see the input file).

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
#include "mhd/mhd.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "globals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "pgen/pgen.hpp"

namespace {
// beam packet parameters
Real vet_beam_e = 1.0;
Real vet_beam_n1 = 1.0/M_SQRT2;
Real vet_beam_n2 = 1.0/M_SQRT2;
Real vet_beam_f = 1.0 - 1.0e-6;
int  vet_test_id = 0;   // 0 = beam, 1 = pulse1d, 2 = thick_pulse
// m1-mhd2: a UNIFORM field problem/vet_b0_1..3 under <mhd> (force free, div B = 0), so
// that the uniform-gas tests (equil, marshak, thick_pulse/tophat, advect_uniform) run the
// |B|^2/2 bookkeeping of every VET temperature; 0 under <hydro> (never read)
Real vet_emag = 0.0;
} // namespace

namespace {
//----------------------------------------------------------------------------------------
//! \fn VetSetUniformGas
//! \brief fill hydro's conserved u0 with a gas at rest of density d and internal energy
//! density e.  The VET tests never evolve a velocity, so the conserved state is just
//! (d, 0, 0, 0, e); it is written directly rather than through PrimToCons so that the
//! ghost zones are filled too (the opacity kernel reads them).

void VetSetUniformGas(MeshBlockPack *pmbp, Real d, Real e) {
  // m1-mhd2: <hydro> or <mhd>; under <mhd> IEN also carries the uniform |B|^2/2
  radvet::FluidRef flr = radvet::FluidRef::Get(pmbp);
  if (!flr.on) return;
  const Real eb = flr.mhd ? vet_emag : 0.0;
  auto &indcs = pmbp->pmesh->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto uh = flr.u0;
  const bool mhdg = flr.mhd;
  par_for("vet_gas_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = mhdg ? (e + eb) : e;
  });
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationVetTests()

void ProblemGenerator::RadiationVetTests(ParameterInput *pin, const bool restart) {
  // NOTE on RESTARTS (milestone 1c-B).  The evolved fields come back from the restart
  // file, but the function POINTERS do not: a test with user boundary conditions or a
  // user source term has to re-enrol them, and the namespace constants those hooks read
  // have to be recomputed, or the run dies in pgen.cpp's "not enrolled during restart"
  // check.  So this function is NOT skipped on a restart; each branch runs its
  // parameter reads and its enrolment and then returns just before writing the arrays.
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pradvet == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "pgen_name = rad_vet_beam requires a <rad_vet> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string test = pin->GetOrAddString("problem","vet_test","beam");
  // m1-mhd2: the uniform-gas tests below and radwave have an <mhd> mode (FluidRef,
  // optional uniform field problem/vet_b0_1..3).  The others (advect_pulse and the rest
  // of rad_vet_tests2.cpp) read phydro directly.
  const bool mhd_ok = (test.compare("beam") == 0 || test.compare("pulse1d") == 0 ||
                       test.compare("equil") == 0 ||
                       test.compare("marshak") == 0 || test.compare("thick_pulse") == 0 ||
                       test.compare("tophat") == 0 || test.compare("radwave") == 0 ||
                       test.compare("advect_uniform") == 0 ||
                       test.compare("sph_shell") == 0);
  if (pmbp->pmhd != nullptr && !mhd_ok) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/vet_test = '" << test << "' is hydro only; with <mhd> "
      << "only beam, pulse1d, equil, marshak, thick_pulse, tophat, "
      << "radwave and advect_uniform (and the box_convection pgen) are supported"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // m1-mhd2: the uniform field, read only under <mhd> (a hydro run's parameter dump is
  // unchanged); set on every face incl. the ghosts, so it is exactly force free
  if (pmbp->pmhd == nullptr && (pin->DoesParameterExist("problem", "vet_b0_1") ||
                                pin->DoesParameterExist("problem", "vet_b0_2") ||
                                pin->DoesParameterExist("problem", "vet_b0_3"))) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "problem/vet_b0_* needs an <mhd> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmbp->pmhd != nullptr) {
    const Real bx = pin->GetOrAddReal("problem", "vet_b0_1", 0.0);
    const Real by = pin->GetOrAddReal("problem", "vet_b0_2", 0.0);
    const Real bz = pin->GetOrAddReal("problem", "vet_b0_3", 0.0);
    vet_emag = 0.5*(bx*bx + by*by + bz*bz);
    if (!restart) {
      auto b = pmbp->pmhd->b0;
      Kokkos::deep_copy(b.x1f, bx);
      Kokkos::deep_copy(b.x2f, by);
      Kokkos::deep_copy(b.x3f, bz);
    }
  }
  if (test.compare("beam") != 0 && test.compare("pulse1d") != 0 &&
      pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/vet_test = '" << test << "' needs a <hydro> block"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real cl = pmbp->pradvet->c_light;
  Real efl = pmbp->pradvet->e_floor;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->pradvet->u0;

  if (test.compare("beam") == 0) {
    vet_test_id = 0;
    // beam direction, amplitude and reduced flux, all optional
    Real theta = pin->GetOrAddReal("problem","beam_angle",0.25*M_PI);
    vet_beam_n1 = cos(theta);
    vet_beam_n2 = sin(theta);
    vet_beam_e = pin->GetOrAddReal("problem","beam_e",1.0);
    vet_beam_f = pin->GetOrAddReal("problem","beam_f",1.0-1.0e-6);

    if (restart) return;
    // problem/beam_packet = true (default off, read only when named): a free-streaming
    // Gaussian PACKET travelling along (n1, n2) with reduced flux beam_f.  Transverse
    // 1/e half-width packet_w, longitudinal packet_l, centre (packet_x0, packet_y0).
    // Use periodic or outflow boundaries.
    if (pin->DoesParameterExist("problem","beam_packet") &&
        pin->GetBoolean("problem","beam_packet")) {
      Real px0 = pin->GetReal("problem","packet_x0");
      Real py0 = pin->GetReal("problem","packet_y0");
      Real pw = pin->GetReal("problem","packet_w");
      Real pl = pin->GetReal("problem","packet_l");
      Real be = vet_beam_e, bf = vet_beam_f, bn1 = vet_beam_n1, bn2 = vet_beam_n2;
      int nx2 = indcs.nx2;
      par_for("vet_packet_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
        Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
        Real spar = (x1v - px0)*bn1 + (x2v - py0)*bn2;
        Real sper = -(x1v - px0)*bn2 + (x2v - py0)*bn1;
        Real e = be*exp(-SQR(sper/pw) - SQR(spar/pl));
        u0(m,radvet::VET_E,k,j,i) = efl + e;
        u0(m,radvet::VET_F1,k,j,i) = bf*cl*e*bn1;
        u0(m,radvet::VET_F2,k,j,i) = bf*cl*e*bn2;
        u0(m,radvet::VET_F3,k,j,i) = 0.0;
      });
      return;
    }
    // vacuum everywhere
    par_for("vet_beam_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radvet::VET_E,k,j,i) = efl;
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("pulse1d") == 0) {
    vet_test_id = 1;
    Real amp = pin->GetOrAddReal("problem","pulse_amp",1.0);
    Real wid = pin->GetOrAddReal("problem","pulse_width",0.1);
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.5);
    Real bg = pin->GetOrAddReal("problem","pulse_bg",1.0e-4);

    if (restart) return;
    par_for("vet_pulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real e = bg + amp*exp(-SQR((x1v - x0)/wid));
      u0(m,radvet::VET_E,k,j,i) = e;
      u0(m,radvet::VET_F1,k,j,i) = cl*e;   // f = 1 exactly: rigid translation at +c
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("thick_pulse") == 0 || test.compare("tophat") == 0) {
    // T3 / T3c.  Static uniform gas; the opacity is whatever <rad_vet> says, and with
    // kappa_s only the gas neither emits nor absorbs.
    bool tophat = (test.compare("tophat") == 0);
    vet_test_id = tophat ? 3 : 2;
    Real amp = pin->GetOrAddReal("problem","pulse_amp",1.0);
    Real wid = pin->GetOrAddReal("problem","pulse_width",0.05);
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.5);
    Real bg = pin->GetOrAddReal("problem","pulse_bg",1.0e-3);
    Real nyq = pin->GetOrAddReal("problem","nyquist_amp",0.0);
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real egas = pin->GetOrAddReal("problem","gas_eint",1.0);
    // MILESTONE 3b phase B: the same pulse in 2-D / 3-D.  The Gaussian becomes
    // ISOTROPIC in whichever directions the mesh has, exp(-r^2/wid^2) with r the
    // distance from (pulse_x0, pulse_y0, pulse_z0), so that the variance of E must grow
    // at 2 D in EVERY direction (gate G1).  A 1-D mesh is bitwise what it was.
    Real y0 = pin->GetOrAddReal("problem","pulse_y0",0.0);
    Real z0 = pin->GetOrAddReal("problem","pulse_z0",0.0);
    // problem/pulse_1d keeps the Gaussian a function of x1 ALONE on a multi-D mesh:
    // gate G3 of 3b phase B, where transport = implicit must reproduce implicit_x1.
    bool p1d = pin->GetOrAddBoolean("problem","pulse_1d",false);
    bool md = pmy_mesh_->multi_d && !p1d;
    bool td = pmy_mesh_->three_d && !p1d;
    int nx2 = indcs.nx2, nx3 = indcs.nx3;
    // gas_vx: a uniform x1 velocity of the (evolved, still uniform) gas, for
    // the dynamic-diffusion test (Jiang 2021 Sect. 4.4 / Ma et al. 2025 Sect. 4.2.1)
    Real gvx = pin->GetOrAddReal("problem","gas_vx",0.0);
    if (restart) return;
    VetSetUniformGas(pmbp, dgas, egas);
    if (gvx != 0.0) {
      auto uhv = pmbp->phydro->u0;
      par_for("vet_thickpulse_v", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        uhv(m,IM1,k,j,i) = uhv(m,IDN,k,j,i)*gvx;
        uhv(m,IEN,k,j,i) += 0.5*uhv(m,IDN,k,j,i)*gvx*gvx;
      });
    }

    // the Nyquist seed is (-1)^(global cell index), so it must not restart at each
    // MeshBlock: the index is rebuilt from the block's own x1min
    Real mx1min = pmy_mesh_->mesh_size.x1min;
    Real dxm = (pmy_mesh_->mesh_size.x1max - mx1min)/
               static_cast<Real>(pmy_mesh_->mesh_indcs.nx1);
    par_for("vet_thickpulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real r2 = SQR(x1v - x0);
      Real rmx = fabs(x1v - x0);
      if (md) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        r2 += SQR(x2v - y0);
        rmx = fmax(rmx, fabs(x2v - y0));
      }
      if (td) {
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        r2 += SQR(x3v - z0);
        rmx = fmax(rmx, fabs(x3v - z0));
      }
      Real e = bg;
      if (tophat) {
        e += (rmx < wid) ? amp : 0.0;
      } else {
        e += amp*exp(-r2/SQR(wid));
      }
      if (nyq != 0.0) {
        int ig = static_cast<int>(floor((x1v - mx1min)/dxm));
        e += ((ig % 2) == 0) ? nyq : -nyq;
      }
      u0(m,radvet::VET_E,k,j,i) = fmax(e, efl);
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
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
    vet_test_id = 6;
    Real t0 = pin->GetReal("problem","pulse_t0");
    Real t1 = pin->GetReal("problem","pulse_t1");
    Real d0 = pin->GetReal("problem","gas_rho");
    Real wid = pin->GetReal("problem","pulse_width");
    Real x0 = pin->GetOrAddReal("problem","pulse_x0",0.0);
    Real vx = pin->GetOrAddReal("problem","pulse_v",0.0);
    Real ar = pmbp->pradvet->arad;
    Real kap = pmbp->pradvet->kappa_f + pmbp->pradvet->kappa_s;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real ptot = d0*t0 + ar*t0*t0*t0*t0/3.0;
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    par_for("vet_advpulse_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
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
      u0(m,radvet::VET_E,k,j,i) = fmax(e, efl);
      u0(m,radvet::VET_F1,k,j,i) = ff;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
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
    vet_test_id = 7;
    Real tg = pin->GetReal("problem","gas_temp");
    Real d0 = pin->GetReal("problem","gas_rho");
    Real vx = pin->GetOrAddReal("problem","pulse_v",0.0);
    Real ar = pmbp->pradvet->arad;
    radvet::FluidRef flr = radvet::FluidRef::Get(pmbp);   // m1-mhd2
    Real gm1 = flr.eos.gamma - 1.0;
    const Real eb = flr.mhd ? vet_emag : 0.0;
    Real bb = vx/cl;
    // fixed point of f = beta (1 + chi(f)), a handful of Picard steps (contraction
    // factor ~ beta)
    Real ff = (4.0/3.0)*bb;
    Real chi = 1.0/3.0;
    for (int it=0; it<50; ++it) {
      chi = radvet::VetChi(fabs(ff));
      ff = bb*(1.0 + chi);
    }
    Real erad = ar*tg*tg*tg*tg/((1.0 + bb*bb) - 2.0*bb*ff + bb*bb*chi);
    Real frad = cl*erad*ff;
    if (global_variable::my_rank == 0) {
      std::cout << "  vet_test = advect_uniform: beta=" << bb << " f=" << ff
                << " E=" << erad << " F=" << frad << " F/((4/3)vE)="
                << frad/((4.0/3.0)*vx*erad) << std::endl;
    }
    if (restart) return;
    auto uh = flr.u0;
    const bool mhdu = flr.mhd;
    par_for("vet_advunif_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      uh(m,IDN,k,j,i) = d0;
      uh(m,IM1,k,j,i) = d0*vx;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d0*tg/gm1 + 0.5*d0*vx*vx;
      if (mhdu) uh(m,IEN,k,j,i) += eb;
      u0(m,radvet::VET_E,k,j,i) = erad;
      u0(m,radvet::VET_F1,k,j,i) = frad;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("marshak") == 0) {
    // T6.  Non-equilibrium Marshak wave: a cold, static, uniform slab of constant
    // opacity heated from x1min by an incident radiation bath.  The material heat
    // capacity is the code's CONSTANT c_v, not the Su-Olson alpha T^3 (which the EOS
    // cannot represent); tests_m1/t6_marshak.py's own S_N reference is run with the same
    // constant c_v (--cv), so the comparison is like for like.  The gas does not move
    // (gas_feedback still writes the energy exchange; the momentum exchange is left in,
    // and stays ~1e-10 of the pressure over the run).
    vet_test_id = 8;
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real tini = pin->GetReal("problem","gas_temp");
    Real ar = pmbp->pradvet->arad;
    Real gm1 = radvet::FluidRef::Get(pmbp).eos.gamma - 1.0;   // m1-mhd2
    Real e_ini = fmax(ar*tini*tini*tini*tini, efl);
    if (restart) return;
    VetSetUniformGas(pmbp, dgas, dgas*tini/gm1);
    par_for("vet_marshak_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radvet::VET_E,k,j,i) = e_ini;
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("hot_blob") == 0) {
    // (TEST ONLY, preconditioner orientation study): an optically thick hot
    // region in a thinner, cooler medium, gas at rest, E = arad T^4, F = 0, evolved as
    // rad-hydro.  The region is an ellipsoid with semi-axes blob_a1..a3 about
    // (blob_x1, blob_x2, blob_x3); a semi-axis <= 0 drops that direction, so (a, 0, 0) is
    // a slab normal to x1 and (a, a, a) a sphere.  Smooth edge: s = (1 -
    // tanh((r-1)/w))/2.
    vet_test_id = 10;
    Real rc = pin->GetReal("problem","blob_rho");
    Real tc = pin->GetReal("problem","blob_temp");
    Real rb = pin->GetReal("problem","bg_rho");
    Real tb = pin->GetReal("problem","bg_temp");
    Real ew = pin->GetOrAddReal("problem","blob_edge",0.1);
    Real ca[3], aa[3];
    for (int d = 0; d < 3; ++d) {
      ca[d] = pin->GetOrAddReal("problem","blob_x" + std::to_string(d+1),0.0);
      aa[d] = pin->GetOrAddReal("problem","blob_a" + std::to_string(d+1),0.0);
    }
    const Real c1 = ca[0], c2 = ca[1], c3 = ca[2];
    const Real ia1 = (aa[0] > 0.0) ? 1.0/aa[0] : 0.0;
    const Real ia2 = (aa[1] > 0.0) ? 1.0/aa[1] : 0.0;
    const Real ia3 = (aa[2] > 0.0) ? 1.0/aa[2] : 0.0;
    Real ar = pmbp->pradvet->arad;
    Real gm1 = radvet::FluidRef::Get(pmbp).eos.gamma - 1.0;
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    int nx2 = indcs.nx2, nx3 = indcs.nx3;
    par_for("vet_hotblob_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      Real x3v = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
      Real r = sqrt(SQR((x1v - c1)*ia1) + SQR((x2v - c2)*ia2) + SQR((x3v - c3)*ia3));
      Real s = 0.5*(1.0 - tanh((r - 1.0)/ew));
      Real d = rb + (rc - rb)*s;
      Real t = tb + (tc - tb)*s;
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*t/gm1;
      u0(m,radvet::VET_E,k,j,i) = fmax(ar*t*t*t*t, efl);
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("equil") == 0) {
    // T5.  Uniform single zone: nothing has a gradient, so the transport is exactly
    // zero and what is exercised is the implicit solve alone.
    vet_test_id = 5;
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0e-7);
    Real egas = pin->GetReal("problem","gas_eint");
    Real erad = pin->GetReal("problem","e_rad");
    Real gvx = pin->GetOrAddReal("problem","gas_vx",0.0);   // (J21 4.1)
    if (restart) return;
    VetSetUniformGas(pmbp, dgas, egas);
    if (gvx != 0.0) {
      auto uhv = pmbp->phydro->u0;
      par_for("vet_equil_v", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        uhv(m,IM1,k,j,i) = uhv(m,IDN,k,j,i)*gvx;
        uhv(m,IEN,k,j,i) += 0.5*uhv(m,IDN,k,j,i)*gvx*gvx;
      });
    }
    par_for("vet_equil_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radvet::VET_E,k,j,i) = erad;
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else {
    // milestone 1d (shadow, radshock): src/pgen/tests/rad_vet_tests2.cpp, which also
    // owns the "not implemented" fatal
    RadiationVetTests2(pin, restart);
  }
  return;
}
