//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tests2.cpp
//! \brief the milestone-1d test problems of the grey VET module, dispatched from
//! rad_vet_tests.cpp when <problem>/vet_test names none of the 1a-1c tests.  Kept in a
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
//!             hydro stage and <rad_vet>/gas_feedback = false keeps the coupling from
//!             writing it back, so T = T0 exactly for the whole run and the opacity
//!             sigma = 0.1 (T/T0)^-3.5 (rho/rho0)^2 per cm is a pure function of rho.
//!             (Gas heating on a light-crossing time is negligible in any case: the gas
//!             heat capacity here exceeds a T0^4 by ten orders of magnitude.)
//!
//!   radshock  (design sect. 8, T7) Lowrie & Edwards (2008) steady radiative shock,
//!             initialised from the semi-analytic profile written by
//!             tests_m1/t7_radshock.py --write-ref, so that only a local relaxation is
//!             needed rather than the many flow-through times a step start would cost.
//!             With <problem>/vet_shock_ref unset the two far-field states are laid down
//!             with a tanh jump of <problem>/vet_shock_w cells instead.  Both x1 faces
//!             are Dirichlet at the far-field states (user BC), in the shock frame.
//!
//!   atmosphere (milestone 3a2, gate I6b) STATIC plane-parallel grey column with a
//!             constant imposed flux at the bottom and a free surface at the top.  The
//!             gas is PRESCRIBED (a user source term re-imposes it every stage and
//!             <rad_vet>/gas_feedback = false), the opacity is pure scattering
//!             (kappa_P = kappa_E = 0) and rho falls exponentially with height,
//!
//!               rho(z) = atm_rho_top exp[(x1max - z)/atm_scale_h],
//!               tau(z) = kappa rho_top H (exp[(x1max - z)/H] - 1)   (downward from the
//!                                                                    TOP MESH FACE),
//!
//!             so the steady state is the closed-form M1 moment solution of
//!             bench/vet_stage2/ic/build_ic.py: with F = const,
//!             dP_rad/dtau = F/c and P_rad = chi(f) E = (F/c)(chi/f) give
//!
//!               chi(f)/f = q0 + tau,      E(tau) = F/(c f(tau)),
//!
//!             on the DECREASING branch of chi/f (which has a minimum 0.89806 at
//!             f = 0.8250; q0 > that, so the solution never leaves it).  The constant q0
//!             is fixed by the surface: the Marshak condition F = c q E_top means
//!             f(0) = q = <rad_vet>/marshak_q exactly, hence the CONSISTENT PAIR is
//!
//!               q0 = chi(marshak_q)/marshak_q,   and for marshak_q = 1/2,
//!               chi(1/2) = 0.46481586, q0 = 0.92963172.
//!
//!             The initial condition is the EDDINGTON solution E = 3(F/c)(tau + 2/3),
//!             which is 40 % off at the surface, so the gate measures the convergence to
//!             the steady state and not the initial condition.
//!
//!   radwave   (milestone 3b phase E, gate T10) LINEAR RADIATION-MODIFIED ACOUSTIC WAVE
//!             in a uniform, optically thick, radiation-pressure-significant medium,
//!             laid down along x1, x2, x3 or the x1-x2 diagonal (<problem>/radwave_dir).
//!             The SAME wave rotated onto a transverse axis is what isolates the
//!             transverse gas-radiation coupling of transport = implicit from the x1
//!             path that milestones 3a-3b already validate.  See the branch below for
//!             the eigenmode and the mixture exponents.

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "hydro/hydro_am.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "globals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "pgen/pgen.hpp"

namespace {
// shadow: the clump geometry and the background state, needed by the source term that
// re-imposes the gas
Real vet_sh2_rho0 = 1.0, vet_sh2_ratio = 1.0e3, vet_sh2_t0 = 290.0;
Real vet_sh2_x0 = 0.5, vet_sh2_y0 = 0.0, vet_sh2_ax = 0.1, vet_sh2_ay = 0.06;
Real vet_sh2_delta = 10.0, vet_sh2_gm1 = 2.0/3.0;
bool vet_sh2_dsq = false, vet_sh2_peq = false;
// vet_test = cyl: up to 3 cylinders (along z) in a thin background, each
// an EMITTER (opaque, T fixed) or an ABSORBER (opaque, T ~ 0); the gas is prescribed
// and re-imposed every stage (RadVetCylGas).  (x, y, radius, rho, T) per cylinder.
int vet_cy_n = 0;
// up to 8 objects; (x, y, r, rho, T, rect, half-length, half-width,
// angle):
// rect = 1 is an oriented rectangle about (x, y), long axis at `angle` (radians) to x1
Real vet_cy_p[24][9] = {};
// the object table captured by value in the reset kernel (namespace scope: CUDA does
// not take a function-local type in an extended lambda)
struct VetCylP {Real v[24][9];};
Real vet_cy_rhob = 1.0e-6, vet_cy_tb = 1.0e-3, vet_cy_gm1 = 2.0/3.0;

// radshock: the two Dirichlet end states (rho, v, E_gas, E_rad, F_rad)
Real vet_rs_l[5] = {1.0, 0.0, 1.0, 0.0, 0.0};
Real vet_rs_r[5] = {1.0, 0.0, 1.0, 0.0, 0.0};

// atmosphere: the prescribed gas
Real vet_at_rho = 0.128, vet_at_h = 0.113, vet_at_ztop = 1.0, vet_at_egas = 1.0;

// sph_atm (STAGE S2, tests_m1/runs_5b_sp_s2): the held gas rho = rho0 (r/r_in)^-n
Real vet_sa_rho0 = 1.0, vet_sa_rin = 1.0, vet_sa_n = 0.0;
// atm_scale_h = H > 0: rho = rho0 exp(-(r - r_in)/H) instead
Real vet_sa_h = 0.0;

// STAGE CS1 (cubed sphere, sph_shell): the lateral mode axis of P_l(mu), mu = rhat . axis
Real vet_cs_ax[3] = {0.0, 0.0, 1.0};

//----------------------------------------------------------------------------------------
//! \fn VetLegendre
//! \brief P_l(x) by the upward recurrence

KOKKOS_INLINE_FUNCTION
Real VetLegendre(const int l, const Real x) {
  Real p0 = 1.0, p1 = x;
  if (l == 0) return p0;
  for (int n = 1; n < l; ++n) {
    const Real p2 = ((2.0*n + 1.0)*x*p1 - n*p0)/(n + 1.0);
    p0 = p1;
    p1 = p2;
  }
  return p1;
}

//----------------------------------------------------------------------------------------
//! \fn VetCsMu
//! \brief rhat . axis of the centre of cell (k,j) of a cubed-sphere MeshBlock

KOKKOS_INLINE_FUNCTION
Real VetCsMu(const int panel, const int k, const int j, const int ks, const int js,
            const int nx2, const int nx3, const Real x2mn, const Real x2mx,
            const Real x3mn, const Real x3mx, const Real ax, const Real ay,
            const Real az) {
  const Real xi = 0.25*M_PI*CellCenterX(j - js, nx2, x2mn, x2mx);
  const Real eta = 0.25*M_PI*CellCenterX(k - ks, nx3, x3mn, x3mx);
  Real q[3];
  cubed_sphere::PanelToCart(panel, xi, eta, q);
  return q[0]*ax + q[1]*ay + q[2]*az;
}

//----------------------------------------------------------------------------------------
//! \fn VetAtmTau
//! \brief the Rosseland optical depth measured DOWNWARD from the top mesh face z = ztop
//! for rho(z) = rho_top exp[(ztop - z)/H] and a constant kappa.

KOKKOS_INLINE_FUNCTION
Real VetAtmTau(const Real z, const Real rho_top, const Real kap, const Real hh,
              const Real ztop) {
  Real s = (ztop - z)/hh;
  // expm1 keeps the top cells (s << 1) accurate: tau -> kappa rho_top (ztop - z) there
  return kap*rho_top*hh*expm1(s);
}

//----------------------------------------------------------------------------------------
//! \fn VetShadowRho
//! \brief the clump density profile: a Fermi-like smoothing of the ellipse boundary,
//! rho = rho0 [1 + (ratio - 1)/(1 + exp(delta (D - 1)))], D^2 = (dx/ax)^2 + (dy/ay)^2.

KOKKOS_INLINE_FUNCTION
Real VetShadowRho(const Real x, const Real y, const Real rho0, const Real ratio,
                 const Real x0, const Real y0, const Real ax, const Real ay,
                 const Real delta, const bool dsq = false) {
  // shadow_dsq: D^2 - 1 in the exponent (Jiang 2021 / Ma et al. 2025)
  Real dd = SQR((x - x0)/ax) + SQR((y - y0)/ay);
  if (!dsq) dd = sqrt(dd);
  Real ex = delta*(dd - 1.0);
  // exp() of a large positive argument overflows; the limit is rho0
  Real w = (ex > 60.0) ? 0.0 : (1.0/(1.0 + exp(ex)));
  return rho0*(1.0 + (ratio - 1.0)*w);
}

} // namespace

// prototypes for the hooks enrolled below
void RadVetShadowGas(Mesh *pm, const Real bdt);
void RadVetCylGas(Mesh *pm, const Real bdt);
void RadVetShockBC(Mesh *pm);
void RadVetAtmGas(Mesh *pm, const Real bdt);
void RadVetAmHistory(HistoryData *pdata, Mesh *pm);
// problem/am_hist_ebnd (read only when named): RadVetAmHistory adds column 15
// Eb = the accumulated radiation energy the scheme put in through the physical x1 faces
// (Ebnd, the hesdirk2 stage weights); absent = the 14 columns as before
static bool vet_am_ebnd = false;

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationVetTests2()
//! \brief the milestone-1d branches, reached from RadiationVetTests() when
//! <problem>/vet_test names none of the 1a-1c tests.

void ProblemGenerator::RadiationVetTests2(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  std::string test = pin->GetOrAddString("problem","vet_test","beam");
  // removed test options (the held / structured sph_shell gas, the run-end shell report
  // and the radwave / cubed-sphere user histories) stop the run when named
  for (const char *k : {"gas_hold", "gas_hold_v", "gas_hold_rho", "core_r", "clump_a",
                        "screen_r0", "clump_tmod_period", "flux_mod_period",
                        "shell_report", "cell_dump", "rw_hist"}) {
    if (pin->DoesParameterExist("problem", k)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<problem>/" << k << " was removed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  Real cl = pmbp->pradvet->c_light;
  Real efl = pmbp->pradvet->e_floor;
  Real ar = pmbp->pradvet->arad;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto u0 = pmbp->pradvet->u0;

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
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    vet_sh2_rho0 = rho0;
    vet_sh2_ratio = ratio;
    vet_sh2_t0 = t0;
    vet_sh2_x0 = x0;
    vet_sh2_y0 = y0;
    vet_sh2_ax = ax;
    vet_sh2_ay = ay;
    vet_sh2_delta = delta;
    // shadow_dsq (exponent in D^2) and shadow_peq (pressure equilibrium,
    // T = t0 rho0/rho: Jiang et al. 2012 / Jiang 2021 / Ma et al. 2025)
    vet_sh2_dsq = pin->GetOrAddBoolean("problem","shadow_dsq",false);
    vet_sh2_peq = pin->GetOrAddBoolean("problem","shadow_peq",false);
    const bool dsq = vet_sh2_dsq, peq = vet_sh2_peq;
    vet_sh2_gm1 = gm1;
    const Real ein = ar*tr*tr*tr*tr;
    user_srcs_func = RadVetShadowGas;
    if (global_variable::my_rank == 0) {
      std::cout << "  vet_test = shadow: E_ambient=" << ar*t0*t0*t0*t0
                << " E_source=" << ein << " ratio="
                << (ein/(ar*t0*t0*t0*t0)) << std::endl;
    }
    if (restart) return;
    Real e_amb = fmax(ar*t0*t0*t0*t0, efl);
    auto uh = pmbp->phydro->u0;
    par_for("vet_shadow_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
      Real d = VetShadowRho(x1v, x2v, rho0, ratio, x0, y0, ax, ay, delta, dsq);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = (peq ? rho0 : d)*t0/gm1;
      u0(m,radvet::VET_E,k,j,i) = e_amb;
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("cyl") == 0) {
    // crossing beams / shadow by opaque cylinders (RadVetCylGas)
    vet_cy_n = pin->GetOrAddInteger("problem","cyl_n",2);
    if (vet_cy_n < 1 || vet_cy_n > 24) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "problem/cyl_n must be 1..24" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int c = 0; c < vet_cy_n; ++c) {
      const std::string b = "cyl" + std::to_string(c) + "_";
      vet_cy_p[c][0] = pin->GetReal("problem", b + "x");
      vet_cy_p[c][1] = pin->GetReal("problem", b + "y");
      vet_cy_p[c][2] = pin->GetReal("problem", b + "r");
      vet_cy_p[c][3] = pin->GetReal("problem", b + "rho");
      vet_cy_p[c][4] = pin->GetReal("problem", b + "t");
      // cyl{c}_shape = rect, with _hl, _hw and _ang (degrees)
      if (pin->GetOrAddString("problem", b + "shape", "cyl").compare("rect") == 0) {
        vet_cy_p[c][5] = 1.0;
        vet_cy_p[c][6] = pin->GetReal("problem", b + "hl");
        vet_cy_p[c][7] = pin->GetReal("problem", b + "hw");
        vet_cy_p[c][8] = pin->GetReal("problem", b + "ang")*M_PI/180.0;
      }
    }
    vet_cy_rhob = pin->GetOrAddReal("problem","cyl_rho_bg",1.0e-6);
    vet_cy_tb = pin->GetOrAddReal("problem","cyl_t_bg",1.0e-3);
    vet_cy_gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    user_srcs_func = RadVetCylGas;
    if (restart) return;
    RadVetCylGas(pmy_mesh_, 0.0);
    const Real eb = fmax(ar*SQR(SQR(vet_cy_tb)), efl);
    par_for("vet_cyl_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,radvet::VET_E,k,j,i) = eb;
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("radshock") == 0) {
    // the two far-field states, in the cgs scaling printed by t7_radshock.py --athinput
    Real dl = pin->GetReal("problem","vet_shock_rho_l");
    Real vl = pin->GetReal("problem","vet_shock_v_l");
    Real tl = pin->GetReal("problem","vet_shock_t_l");
    Real dr = pin->GetReal("problem","vet_shock_rho_r");
    Real vr = pin->GetReal("problem","vet_shock_v_r");
    Real trr = pin->GetReal("problem","vet_shock_t_r");
    Real xs = pin->GetReal("problem","vet_shock_xs");
    Real wid = pin->GetOrAddReal("problem","vet_shock_w",2.0);
    std::string ref = pin->GetOrAddString("problem","vet_shock_ref","none");
    // the code temperature per kelvin.  In honest cgs with the ideal-gas EOS,
    // T_code = (gamma-1) e/rho = [k/(mu m_H)] T_kelvin, so tunit = k/(mu m_H) and
    // <rad_vet>/arad must be a_cgs/tunit^4 for E = arad T_code^4 to be erg/cm^3.
    Real tun = pin->GetOrAddReal("problem","vet_shock_tunit",1.0);
    tl *= tun;
    trr *= tun;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    vet_rs_l[0] = dl;
    vet_rs_l[1] = dl*vl;
    vet_rs_l[2] = dl*tl/gm1 + 0.5*dl*vl*vl;
    vet_rs_l[3] = ar*tl*tl*tl*tl;
    vet_rs_l[4] = 0.0;
    vet_rs_r[0] = dr;
    vet_rs_r[1] = dr*vr;
    vet_rs_r[2] = dr*trr/gm1 + 0.5*dr*vr*vr;
    vet_rs_r[3] = ar*trr*trr*trr*trr;
    vet_rs_r[4] = 0.0;
    user_bcs_func = RadVetShockBC;
    if (restart) return;

    // the semi-analytic profile, if one was given: columns x, rho, v, T_gas, T_rad
    int nref = 0;
    std::vector<Real> vx, vd, vv, vtg, vtr;
    if (ref.compare("none") != 0) {
      std::ifstream fp(ref.c_str());
      if (!fp.is_open()) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cannot open <problem>/vet_shock_ref = '" << ref << "'"
          << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string line;
      while (std::getline(fp, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        Real a, b, c, d, e;
        if (!(ss >> a >> b >> c >> d >> e)) continue;
        vx.push_back(a);
        vd.push_back(b);
        vv.push_back(c);
        vtg.push_back(d);
        vtr.push_back(e);
      }
      nref = static_cast<int>(vx.size());
      if (nref < 2) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "vet_shock_ref '" << ref << "' has fewer than 2 rows"
          << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    DualArray2D<Real> rf("rs_ref", 5, (nref > 0) ? nref : 1);
    for (int n=0; n<nref; ++n) {
      rf.h_view(0,n) = vx[n];
      rf.h_view(1,n) = vd[n];
      rf.h_view(2,n) = vv[n];
      rf.h_view(3,n) = vtg[n];
      rf.h_view(4,n) = vtr[n];
    }
    rf.modify_host();
    rf.sync_device();
    auto rf_ = rf;
    int nr = nref;
    auto uh = pmbp->phydro->u0;
    // the reference table's x is measured FROM THE SHOCK, so it is shifted by xs here
    Real dxr = (nref > 1) ? (rf.h_view(0,1) - rf.h_view(0,0)) : 1.0;
    Real x0r = ((nref > 0) ? rf.h_view(0,0) : 0.0) + xs;
    Real mdx = (pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min)/
               static_cast<Real>(pmy_mesh_->mesh_indcs.nx1);
    par_for("vet_radshock_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real d, v, tg, trd;
      if (nr > 1) {
        // the reference is on a uniform grid: linear interpolation, clamped at the ends
        Real s = (x1v - x0r)/dxr;
        int n = static_cast<int>(floor(s));
        if (n < 0) {n = 0;}
        if (n > nr-2) {n = nr-2;}
        Real w = fmin(fmax(s - static_cast<Real>(n), 0.0), 1.0);
        d = (1.0-w)*rf_.d_view(1,n) + w*rf_.d_view(1,n+1);
        v = (1.0-w)*rf_.d_view(2,n) + w*rf_.d_view(2,n+1);
        tg = tun*((1.0-w)*rf_.d_view(3,n) + w*rf_.d_view(3,n+1));
        trd = tun*((1.0-w)*rf_.d_view(4,n) + w*rf_.d_view(4,n+1));
      } else {
        Real w = 0.5*(1.0 + tanh((x1v - xs)/(fmax(wid,1.0e-30)*mdx)));
        d = (1.0-w)*dl + w*dr;
        v = (1.0-w)*vl + w*vr;
        tg = (1.0-w)*tl + w*trr;
        trd = tg;
      }
      Real erad = ar*trd*trd*trd*trd;
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = d*v;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*tg/gm1 + 0.5*d*v*v;
      u0(m,radvet::VET_E,k,j,i) = fmax(erad, efl);
      // the comoving flux is small in these thick shocks; the lab flux is dominated by
      // the enthalpy term, which the stiff source relaxes to within a few substeps
      u0(m,radvet::VET_F1,k,j,i) = (4.0/3.0)*v*erad;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("atmosphere") == 0) {
    // milestone 3a2, gate I6b: the static grey plane-parallel column.  See the file
    // header for the closed-form steady state and the (marshak_q, q0) pair.
    Real rht = pin->GetOrAddReal("problem","atm_rho_top",0.128);
    Real hh = pin->GetOrAddReal("problem","atm_scale_h",0.113);
    Real fin = pin->GetOrAddReal("problem","atm_flux",1.0);
    Real tg = pin->GetOrAddReal("problem","atm_temp",1.0);
    Real kap = pmbp->pradvet->kappa_s + pmbp->pradvet->kappa_f;
    Real ztop = pmy_mesh_->mesh_size.x1max;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real qq = pmbp->pradvet->marshak_q;
    Real q0 = radvet::VetChi(qq)/qq;
    vet_at_rho = rht;
    vet_at_h = hh;
    vet_at_ztop = ztop;
    vet_at_egas = tg/gm1;
    // the closed-form E at the bottom MESH FACE, which is where tau_bot sits (reported)
    Real taub = VetAtmTau(pmy_mesh_->mesh_size.x1min, rht, kap, hh, ztop);
    Real ebot;
    {
      // invert chi(f)/f = q0 + tau on the decreasing branch by bisection on f
      Real target = q0 + taub;
      Real flo = 1.0e-14, fhi = 0.8250;
      for (int n=0; n<200; ++n) {
        Real fm = 0.5*(flo + fhi);
        if (radvet::VetChi(fm)/fm > target) {flo = fm;} else {fhi = fm;}
      }
      ebot = fin/(cl*0.5*(flo + fhi));
    }
    user_srcs_func = RadVetAtmGas;
    if (global_variable::my_rank == 0) {
      std::cout << "  vet_test = atmosphere: kappa=" << kap << " rho_top=" << rht
                << " H=" << hh << " tau_bot=" << taub << " tau_topcell="
                << VetAtmTau(ztop - (ztop - pmy_mesh_->mesh_size.x1min)/
                            static_cast<Real>(pmy_mesh_->mesh_indcs.nx1), rht, kap, hh,
                            ztop)
                << " F=" << fin << " marshak_q=" << qq << " q0=" << q0
                << " E(tau_bot)=" << ebot << std::endl;
    }
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    Real egas = vet_at_egas;
    // atm_seed (runs_5c_thinstab): a TRANSVERSE perturbation of the initial E,
    // E *= 1 + atm_seed cos(2 pi atm_seed_k (x2 - x2min)/L2), or (-1)^j for
    // atm_seed_k <= 0.  Default 0: the IC is untouched.
    const Real aseed = pin->GetOrAddReal("problem","atm_seed",0.0);
    const int aseedk = pin->GetOrAddInteger("problem","atm_seed_k",1);
    const Real ax2min = pmy_mesh_->mesh_size.x2min;
    const Real ax2len = pmy_mesh_->mesh_size.x2max - ax2min;
    par_for("vet_atm_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real d = rht*exp((ztop - x1v)/hh);
      Real tau = VetAtmTau(x1v, rht, kap, hh, ztop);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*egas;
      // the EDDINGTON guess, deliberately not the M1 answer
      u0(m,radvet::VET_E,k,j,i) = fmax(3.0*(fin/cl)*(tau + 2.0/3.0), efl);
      if (aseed != 0.0) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        Real sm = (aseedk > 0) ? cos(2.0*M_PI*aseedk*(x2v - ax2min)/ax2len)
                               : ((((j - js) & 1) == 0) ? 1.0 : -1.0);
        u0(m,radvet::VET_E,k,j,i) *= 1.0 + aseed*sm;
      }
      u0(m,radvet::VET_F1,k,j,i) = fin;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("radwave") == 0) {
    // MILESTONE 3b phase E (design sect. 14).  LINEAR RADIATION-MODIFIED ACOUSTIC WAVE
    // in a uniform, optically thick, radiation-pressure-significant medium.  The point
    // of the test is that the SAME wave can be laid down along x1, along x2, along x3 or
    // on the diagonal of the x1-x2 plane, so that the transverse gas-radiation coupling
    // (the x2/x3 momentum deposit, the work term, the enthalpy flux A_2/A_3 and the
    // velocity-dependent E0 corrections) is exercised against a path -- x1 -- that every
    // earlier gate has already validated.
    //
    // With tau per wavelength >> 1 and c >> c_s the gas and the radiation are in
    // EQUILIBRIUM DIFFUSION: T_rad = T_gas, and the mixture behaves as a single fluid
    // with Chandrasekhar's generalised adiabatic exponents.  With beta = P_gas/P_tot,
    //
    //   Gamma_1   = beta + (4 - 3 beta)^2 (gamma-1)/[beta + 12 (gamma-1)(1 - beta)]
    //   Gamma_3-1 = (Gamma_1 - beta)/(4 - 3 beta)          (= dlnT/dlnrho on the adiabat)
    //   c_s^2     = Gamma_1 (P_gas + P_rad)/rho
    //
    // (beta = 1 gives Gamma_1 = gamma and Gamma_3-1 = gamma-1; beta = 0 gives 4/3 and
    // 1/3).  The eigenmode of a RIGHT-travelling wave of amplitude A is then
    //
    //   drho/rho = A cos(phi),  dv = n c_s A cos(phi),  dT/T = (Gamma_3-1) A cos(phi),
    //   dE_rad   = 4 E_rad (Gamma_3-1) A cos(phi),
    //   F_rad    = n [ (4/3) dv E_rad + (c/(3 kappa rho)) |k| dE_rad sin(phi) ],
    //
    // the second piece of F being the diffusive flux -c grad E/(3 kappa rho), which is
    // O(1/tau) and is what damps the wave.  phi = k.x with k one wavelength across the
    // box in each direction the wave runs.
    std::string wdir = pin->GetOrAddString("problem","radwave_dir","x1");
    Real amp = pin->GetOrAddReal("problem","radwave_amp",1.0e-6);
    Real d0 = pin->GetOrAddReal("problem","radwave_rho",1.0);
    Real t0 = pin->GetOrAddReal("problem","radwave_t",1.0);
    // m1-mhd2: <hydro> or <mhd>; under <mhd> the uniform field problem/vet_b0_1..3 (set
    // on b0 by RadiationVetTests) adds its |B|^2/2 to IEN.  B parallel to k leaves the
    // wave exactly the hydro one.
    radvet::FluidRef flr = radvet::FluidRef::Get(pmbp);
    Real gm1 = flr.eos.gamma - 1.0;
    const bool mhdw = flr.mhd;
    Real ebw = 0.0;
    if (mhdw) {
      const Real bx = pin->GetOrAddReal("problem", "vet_b0_1", 0.0);
      const Real by = pin->GetOrAddReal("problem", "vet_b0_2", 0.0);
      const Real bz = pin->GetOrAddReal("problem", "vet_b0_3", 0.0);
      ebw = 0.5*(bx*bx + by*by + bz*bz);
    }
    Real kapf = pmbp->pradvet->kappa_f;
    Real pgas = d0*t0;
    Real er0 = ar*t0*t0*t0*t0;
    Real ptot = pgas + er0/3.0;
    Real bet = pgas/ptot;
    Real gam1 = bet + SQR(4.0 - 3.0*bet)*gm1/(bet + 12.0*gm1*(1.0 - bet));
    Real g3m1 = (gam1 - bet)/(4.0 - 3.0*bet);
    Real cs = sqrt(gam1*ptot/d0);
    int mx = 0, my = 0, mz = 0;
    if (wdir.compare("x1") == 0) {
      mx = 1;
    } else if (wdir.compare("x2") == 0) {
      my = 1;
    } else if (wdir.compare("x3") == 0) {
      mz = 1;
    } else if (wdir.compare("xy") == 0) {
      mx = 1;
      my = 1;
    } else if (wdir.compare("xyz") == 0) {
      // the body diagonal of a 3-D box (tests_m1/runs_3r_radwave)
      mx = 1;
      my = 1;
      mz = 1;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<problem>/radwave_dir = '" << wdir
        << "' is not a choice (x1 | x2 | x3 | xy | xyz)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto &msz = pmy_mesh_->mesh_size;
    Real twopi = 2.0*M_PI;
    Real kx = twopi*mx/(msz.x1max - msz.x1min);
    Real ky = twopi*my/(msz.x2max - msz.x2min);
    Real kz = twopi*mz/(msz.x3max - msz.x3min);
    Real kmag = sqrt(kx*kx + ky*ky + kz*kz);
    Real nx = kx/kmag, ny = ky/kmag, nz = kz/kmag;
    Real lam = twopi/kmag;
    Real dvv = amp*cs;
    Real derad = 4.0*er0*g3m1*amp;
    Real fdif = (cl/(3.0*kapf*d0))*kmag*derad;
    if (global_variable::my_rank == 0) {
      std::cout << "  vet_test = radwave: dir=" << wdir << " beta=" << bet
                << " Gamma1=" << gam1 << " Gamma3-1=" << g3m1 << std::endl
                << "    c_s=" << cs << " c/c_s=" << (cl/cs) << " lambda=" << lam
                << " period=" << (lam/cs) << " tau_lambda=" << (kapf*d0*lam)
                << std::endl;
    }
    // STAGE CS2 (gate (d)), read only when named: on the cubed sphere
    // problem/radwave_cs_l = l lays down a TRAVELLING SECTORAL wave instead of the
    // plane wave: S = (x' + i y')^l = sin^l(theta') e^{i l phi'} about the axis
    // radwave_cs_axis = z (its equator crosses four panel seams) | vertex ((1,-1,0)/
    // sqrt2: the great circle through four cube vertices), with the equilibrium-diffusion
    // eigenmode of the plane wave written through the surface gradient,
    //   drho/rho = A Re S, dT/T = (G3-1) A Re S, dE = derad Re S,
    //   v = (A c_s/k) grad_s Im S,  F = (4/3) E0 v - (fdif/k) grad_s Re S,
    // k^2 = l(l+1) <1/r^2> (the shell's volume mean), v and F in COVARIANT panel
    // components (v.e_xi, v.e_eta), F0 on the faces zero (as the plane wave without
    // radwave_eig).
    if (pmy_mesh_->use_cubed_sphere &&
        pin->DoesParameterExist("problem","radwave_cs_l")) {
      const int cl_l = pin->GetInteger("problem","radwave_cs_l");
      const std::string cax = pin->GetOrAddString("problem","radwave_cs_axis","z");
      Real ax[3] = {0.0, 0.0, 1.0}, e1[3] = {1.0, 0.0, 0.0};
      if (cax == "vertex") {
        const Real h = 1.0/std::sqrt(2.0);
        ax[0] = h; ax[1] = -h; ax[2] = 0.0;
        e1[0] = 0.0; e1[1] = 0.0; e1[2] = 1.0;
      }
      Real e2[3] = {ax[1]*e1[2] - ax[2]*e1[1], ax[2]*e1[0] - ax[0]*e1[2],
                    ax[0]*e1[1] - ax[1]*e1[0]};
      const Real r1 = msz.x1min, r2 = msz.x1max;
      const Real ir2 = (r2 - r1)/((r2*r2*r2 - r1*r1*r1)/3.0);
      const Real kk = std::sqrt(cl_l*(cl_l + 1.0)*ir2);
      const Real derad_c = 4.0*er0*g3m1*amp;
      const Real fdif_c = (cl/(3.0*kapf*d0))*kk*derad_c;
      if (global_variable::my_rank == 0) {
        std::cout << "  radwave_cs: sectoral l=" << cl_l << " axis " << cax << " k=" << kk
                  << " lambda=" << (2.0*M_PI/kk) << " period=" << (2.0*M_PI/(kk*cs))
                  << std::endl;
      }
      if (restart) return;
      auto uh = flr.u0;
      auto x1v_ = pmbp->pcoord->x1v;
      auto mbpan = pmbp->pmb->mb_panel;
      const Real e10 = e1[0], e11 = e1[1], e12 = e1[2];
      const Real e20 = e2[0], e21 = e2[1], e22 = e2[2];
      const int ll = cl_l;
      const Real ampv = amp*cs/kk, er0_ = er0;
      par_for("vet_radwave_cs", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        auto &sz = size.d_view(m);
        const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, sz.x2min, sz.x2max);
        const Real eta = 0.25*M_PI*CellCenterX(k-ks, nx3, sz.x3min, sz.x3max);
        Real q[3], ta[3], tb[3];
        cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, q);
        cubed_sphere::PanelTangents(mbpan.d_view(m), xi, eta, ta, tb);
        const Real xp = q[0]*e10 + q[1]*e11 + q[2]*e12;
        const Real yp = q[0]*e20 + q[1]*e21 + q[2]*e22;
        // S = (xp + i yp)^l and (xp + i yp)^(l-1)
        Real sr = 1.0, si = 0.0, pr = 1.0, pi_ = 0.0;
        for (int n = 0; n < ll; ++n) {
          if (n == ll - 1) {pr = sr; pi_ = si;}
          const Real tr = sr*xp - si*yp;
          si = sr*yp + si*xp;
          sr = tr;
        }
        const Real r = x1v_(m,i);
        // grad_s S on the sphere of radius r: (l/r)[(xp+iyp)^(l-1) (e1 + i e2) - S q]
        Real gr[3], gi[3];
        const Real e1v[3] = {e10, e11, e12}, e2v[3] = {e20, e21, e22};
        for (int c = 0; c < 3; ++c) {
          gr[c] = (ll/r)*(pr*e1v[c] - pi_*e2v[c] - sr*q[c]);
          gi[c] = (ll/r)*(pr*e2v[c] + pi_*e1v[c] - si*q[c]);
        }
        Real v[3], f[3];
        for (int c = 0; c < 3; ++c) {
          v[c] = ampv*gi[c];
          f[c] = (4.0/3.0)*er0_*v[c] - (fdif_c/kk)*gr[c];
        }
        const Real vxi = v[0]*ta[0] + v[1]*ta[1] + v[2]*ta[2];
        const Real veta = v[0]*tb[0] + v[1]*tb[1] + v[2]*tb[2];
        const Real fxi = f[0]*ta[0] + f[1]*ta[1] + f[2]*ta[2];
        const Real feta = f[0]*tb[0] + f[1]*tb[1] + f[2]*tb[2];
        const Real d = d0*(1.0 + amp*sr);
        const Real tt = t0*(1.0 + g3m1*amp*sr);
        uh(m,IDN,k,j,i) = d;
        uh(m,IM1,k,j,i) = 0.0;
        uh(m,IM2,k,j,i) = d*vxi;
        uh(m,IM3,k,j,i) = d*veta;
        uh(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        u0(m,radvet::VET_E,k,j,i) = fmax(er0_ + derad_c*sr, efl);
        u0(m,radvet::VET_F1,k,j,i) = 0.0;
        u0(m,radvet::VET_F2,k,j,i) = fxi;
        u0(m,radvet::VET_F3,k,j,i) = feta;
      });
      return;
    }
    if (restart) return;
    auto uh = flr.u0;
    // <problem>/radwave_eig = true (default false; tests_m1/runs_3r_radwave): lay down
    // the EXACT linear eigenmode of the gas + moment system instead of the
    // equilibrium-diffusion one.  The complex amplitudes per unit drho/rho are computed
    // offline (tests_m1/runs_3r_radwave/radwave_disp.py) and passed as
    //   dv = A Re[v^ e^{i phi}],  dT/T0 = A Re[t^ e^{i phi}],
    //   dE/E0 = A Re[e^ e^{i phi}],
    //   F/(c E0) = A Re[f^ e^{i phi}]   (radwave_eig_{v,t,e,f}_{re,im}),
    // and the implicit face fluxes F0 = F - (4/3) v E0 are set on every face, so that the
    // first implicit step does not start from F0 = 0.  Read only when named, so that an
    // input without it is untouched.
    const bool eig = pin->DoesParameterExist("problem","radwave_eig") &&
                     pin->GetBoolean("problem","radwave_eig");
    if (eig) {
      const Real vre = pin->GetReal("problem","radwave_eig_v_re");
      const Real vim = pin->GetReal("problem","radwave_eig_v_im");
      const Real tre = pin->GetReal("problem","radwave_eig_t_re");
      const Real tim = pin->GetReal("problem","radwave_eig_t_im");
      const Real ere = pin->GetReal("problem","radwave_eig_e_re");
      const Real eim = pin->GetReal("problem","radwave_eig_e_im");
      const Real fre = pin->GetReal("problem","radwave_eig_f_re");
      const Real fim = pin->GetReal("problem","radwave_eig_f_im");
      const Real ce0 = cl*er0;
      // (<mhd>, read only when named) the magnetosonic and Alfven modes: the transverse
      // amplitudes along t (zhat x n, the x1-x2 plane) and s = n x t, dv_q = A
      // Re[vq^ e^{i phi}], dB_q = A Re[bq^ e^{i phi}], F_q/(c E0) = A Re[fq^ e^{i phi}]
      // (radwave_eig_{vt,bt,ft,vs,bs,fs}_{re,im}), dB the discrete curl of A_s = i bt/k,
      // A_t = -i bs/k on the edges; radwave_eig_d (default 1) scales drho (0: Alfven)
      auto ekey = [&](const char *q, const char *c) {
        const std::string nm = std::string("radwave_eig_") + q + c;
        return pin->DoesParameterExist("problem",nm) ? pin->GetReal("problem",nm) : 0.0;
      };
      const Real tv[2][4] = {{ekey("vt","_re"), ekey("vt","_im"), ekey("ft","_re"),
                              ekey("ft","_im")},
                             {ekey("vs","_re"), ekey("vs","_im"), ekey("fs","_re"),
                              ekey("fs","_im")}};
      const Real btr = ekey("bt","_re"), bti = ekey("bt","_im");
      const Real bsr = ekey("bs","_re"), bsi = ekey("bs","_im");
      const Real ad = pin->DoesParameterExist("problem","radwave_eig_d") ?
                      pin->GetReal("problem","radwave_eig_d") : 1.0;
      const Real hp = sqrt(nx*nx + ny*ny);
      const Real tq[2][3] = {{(hp > 0.0) ? -ny/hp : 1.0, (hp > 0.0) ? nx/hp : 0.0, 0.0},
                             {(hp > 0.0) ? -nz*nx/hp : 0.0, (hp > 0.0) ? -nz*ny/hp : 1.0,
                              hp}};
      const bool mwv = mhdw && (btr != 0.0 || bti != 0.0 || bsr != 0.0 || bsi != 0.0);
      // <problem>/radwave_v0 (tests_m1/runs_3s_space2): a uniform background drift along
      // the wave vector, with the matching lab flux (4/3) v0 E.  Read only when named.
      const bool drift = pin->DoesParameterExist("problem","radwave_v0");
      const Real v0 = drift ? pin->GetReal("problem","radwave_v0") : 0.0;
      par_for("vet_radwave_eig", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        auto &sz = size.d_view(m);
        Real x1v = CellCenterX(i-is, nx1, sz.x1min, sz.x1max);
        Real x2v = (nx2 > 1) ? CellCenterX(j-js, nx2, sz.x2min, sz.x2max) : 0.0;
        Real x3v = (nx3 > 1) ? CellCenterX(k-ks, nx3, sz.x3min, sz.x3max) : 0.0;
        Real ph = kx*x1v + ky*x2v + kz*x3v;
        Real cp = cos(ph), sp = sin(ph);
        Real d = d0*(1.0 + amp*ad*cp);
        Real tt = t0*(1.0 + amp*(tre*cp - tim*sp));
        Real vv = amp*(vre*cp - vim*sp);
        if (drift) {vv += v0;}
        uh(m,IDN,k,j,i) = d;
        uh(m,IM1,k,j,i) = d*vv*nx;
        uh(m,IM2,k,j,i) = d*vv*ny;
        uh(m,IM3,k,j,i) = d*vv*nz;
        uh(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*vv*vv;
        if (mhdw) uh(m,IEN,k,j,i) += ebw;
        Real fq[2];
        for (int q = 0; q < 2; ++q) {
          const Real vq = amp*(tv[q][0]*cp - tv[q][1]*sp);
          fq[q] = amp*ce0*(tv[q][2]*cp - tv[q][3]*sp);
          uh(m,IM1,k,j,i) += d*vq*tq[q][0];
          uh(m,IM2,k,j,i) += d*vq*tq[q][1];
          uh(m,IM3,k,j,i) += d*vq*tq[q][2];
          uh(m,IEN,k,j,i) += 0.5*d*vq*vq;
        }
        Real ff = amp*ce0*(fre*cp - fim*sp);
        u0(m,radvet::VET_E,k,j,i) = fmax(er0*(1.0 + amp*(ere*cp - eim*sp)), efl);
        if (drift) {
          ff += (4.0/3.0)*v0*(u0(m,radvet::VET_E,k,j,i) - er0)
                + (4.0/3.0)*v0*er0;
        }
        u0(m,radvet::VET_F1,k,j,i) = ff*nx + fq[0]*tq[0][0] + fq[1]*tq[1][0];
        u0(m,radvet::VET_F2,k,j,i) = ff*ny + fq[0]*tq[0][1] + fq[1]*tq[1][1];
        u0(m,radvet::VET_F3,k,j,i) = ff*nz + fq[0]*tq[0][2] + fq[1]*tq[1][2];
      });
      if (mwv) {
        auto b0 = pmbp->pmhd->b0;
        auto b1f = b0.x1f, b2f = b0.x2f, b3f = b0.x3f;
        par_for("vet_radwave_b", DevExeSpace(), 0,nmb1,0,n3,0,n2,0,n1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          auto &sz = size.d_view(m);
          const Real xl = LeftEdgeX(i-is, nx1, sz.x1min, sz.x1max);
          const Real xc = CellCenterX(i-is, nx1, sz.x1min, sz.x1max);
          const Real yl = LeftEdgeX(j-js, nx2, sz.x2min, sz.x2max);
          const Real yc = CellCenterX(j-js, nx2, sz.x2min, sz.x2max);
          const Real zl = LeftEdgeX(k-ks, nx3, sz.x3min, sz.x3max);
          const Real zc = CellCenterX(k-ks, nx3, sz.x3min, sz.x3max);
          const Real d1 = sz.dx1, d2 = sz.dx2, d3 = sz.dx3;
          // component c of A at (x1, x2, x3)
          auto av = [&](const int c, const Real x1, const Real x2, const Real x3) {
            const Real ph = kx*x1 + ky*x2 + kz*x3, cp = cos(ph), sp = sin(ph);
            return amp/kmag*(-(btr*sp + bti*cp)*tq[1][c] + (bsr*sp + bsi*cp)*tq[0][c]);
          };
          if (j < n2 && k < n3) {
            b1f(m,k,j,i) += (av(2, xl, yl + d2, zc) - av(2, xl, yl, zc))/d2 -
                            (av(1, xl, yc, zl + d3) - av(1, xl, yc, zl))/d3;
          }
          if (i < n1 && k < n3) {
            b2f(m,k,j,i) += (av(0, xc, yl, zl + d3) - av(0, xc, yl, zl))/d3 -
                            (av(2, xl + d1, yl, zc) - av(2, xl, yl, zc))/d1;
          }
          if (i < n1 && j < n2) {
            b3f(m,k,j,i) += (av(1, xl + d1, yc, zl) - av(1, xl, yc, zl))/d1 -
                            (av(0, xc, yl + d2, zl) - av(0, xc, yl, zl))/d2;
          }
        });
        par_for("vet_radwave_be", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          uh(m,IEN,k,j,i) += 0.5*(SQR(0.5*(b1f(m,k,j,i) + b1f(m,k,j,i+1))) +
                                  SQR(0.5*(b2f(m,k,j,i) + b2f(m,k,j+1,i))) +
                                  SQR(0.5*(b3f(m,k,j,i) + b3f(m,k+1,j,i)))) - ebw;
        });
      }
      // the comoving face fluxes of the implicit transport (allocated only there)
      const Real f0re = ce0*fre - (4.0/3.0)*er0*vre;
      const Real f0im = ce0*fim - (4.0/3.0)*er0*vim;
      for (int dir = 0; dir < 3; ++dir) {
        DvceArray4D<Real> ff0 = (dir == 0) ? pmbp->pradvet->f0x1 :
                                ((dir == 1) ? pmbp->pradvet->f0x2 : pmbp->pradvet->f0x3);
        if (ff0.extent_int(0) < nmb1 + 1) continue;
        const Real nd = (dir == 0) ? nx : ((dir == 1) ? ny : nz);
        const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
        const int e1 = ff0.extent_int(3) - 1;
        par_for("vet_radwave_eigf", DevExeSpace(), 0,nmb1,0,e3,0,e2,0,e1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          auto &sz = size.d_view(m);
          Real x1v = (dir == 0) ? LeftEdgeX(i-is, nx1, sz.x1min, sz.x1max)
                                : CellCenterX(i-is, nx1, sz.x1min, sz.x1max);
          Real x2v = 0.0, x3v = 0.0;
          if (nx2 > 1) {
            x2v = (dir == 1) ? LeftEdgeX(j-js, nx2, sz.x2min, sz.x2max)
                             : CellCenterX(j-js, nx2, sz.x2min, sz.x2max);
          }
          if (nx3 > 1) {
            x3v = (dir == 2) ? LeftEdgeX(k-ks, nx3, sz.x3min, sz.x3max)
                             : CellCenterX(k-ks, nx3, sz.x3min, sz.x3max);
          }
          Real ph = kx*x1v + ky*x2v + kz*x3v;
          ff0(m,k,j,i) = nd*amp*(f0re*cos(ph) - f0im*sin(ph));
          for (int q = 0; q < 2; ++q) {
            ff0(m,k,j,i) += tq[q][dir]*amp*(
                (ce0*tv[q][2] - (4.0/3.0)*er0*tv[q][0])*cos(ph) -
                (ce0*tv[q][3] - (4.0/3.0)*er0*tv[q][1])*sin(ph));
          }
        });
      }
      return;
    }
    par_for("vet_radwave_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real x2v = (nx2 > 1) ? CellCenterX(j-js, nx2, x2min, x2max) : 0.0;
      Real x3v = (nx3 > 1) ? CellCenterX(k-ks, nx3, x3min, x3max) : 0.0;
      Real ph = kx*x1v + ky*x2v + kz*x3v;
      Real cp = cos(ph), sp = sin(ph);
      Real d = d0*(1.0 + amp*cp);
      Real tt = t0*(1.0 + g3m1*amp*cp);
      Real vv = dvv*cp;
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = d*vv*nx;
      uh(m,IM2,k,j,i) = d*vv*ny;
      uh(m,IM3,k,j,i) = d*vv*nz;
      uh(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*vv*vv;
      if (mhdw) uh(m,IEN,k,j,i) += ebw;
      Real ee = er0 + derad*cp;
      Real ff = (4.0/3.0)*vv*er0 + fdif*sp;
      u0(m,radvet::VET_E,k,j,i) = fmax(ee, efl);
      u0(m,radvet::VET_F1,k,j,i) = ff*nx;
      u0(m,radvet::VET_F2,k,j,i) = ff*ny;
      u0(m,radvet::VET_F3,k,j,i) = ff*nz;
    });
  } else if (test.compare("sph_shell") == 0) {
    // STAGE S1 (tests_m1/runs_5a_sp_s1): a static uniform gas on a spherical-polar
    // wedge and a spherically symmetric radiation field, E = e_out + e_amp
    // exp(-((r - r0)/w)^2) at the cell centroid r = x1v, F = 0.  e_amp = 0 (default)
    // is the flat start of the steady-diffusion gate T-S1 (imposed flux at r_in, the
    // outer end cell held at e_out by implicit_bc_x1max = efix); e_amp > 0 is the
    // symmetric transient of T-sym.  The gas does not move (run it with
    // <rad_vet>/gas_feedback = false and a scattering opacity).
    // (On a Cartesian mesh r is x1: the planar twin used to compare solver behaviour.)
    // <hydro> or <mhd> (the MHD mode: an optional azimuthal field, below)
    if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<problem>/vet_test = sph_shell needs a <hydro> or <mhd> block"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real tgas = pin->GetOrAddReal("problem","gas_temp",1.0e-8);
    Real eout = pin->GetReal("problem","e_out");
    Real eamp = pin->GetOrAddReal("problem","e_amp",0.0);
    Real r0 = pin->GetOrAddReal("problem","r0",0.0);
    Real wid = pin->GetOrAddReal("problem","width",1.0);
    // STAGE S2 (T-S2): e_in > 0 sets the first active cell and the inner ghosts to e_in
    // (the Dirichlet source of the free-streaming point source, implicit_bc_x1min =
    // efix), with F_r = c e_in there
    Real ein = pin->GetOrAddReal("problem","e_in",0.0);
    // m1-sp-order2 (tests_m1/runs_5o_sporder2), each read only when named so that an
    // input without them keeps its parameter dump: lat_amp > 0 multiplies E by
    // 1 + lat_amp cos(lat_l2 pi (x2 - x2min)/(x2max - x2min)) cos(2 pi lat_l3 (x3 -
    // x3min)/(x3max - x3min)) (a smooth lateral mode, even at reflecting x2 walls; sp
    // only);
    // gas_vr sets the gas velocity v_r = gas_vr sin(pi (x1 - x1min)/(x1max - x1min)).
    const Real lamp = pin->DoesParameterExist("problem","lat_amp") ?
                      pin->GetReal("problem","lat_amp") : 0.0;
    const Real ll2 = pin->DoesParameterExist("problem","lat_l2") ?
                     pin->GetReal("problem","lat_l2") : 1.0;
    const Real ll3 = pin->DoesParameterExist("problem","lat_l3") ?
                     pin->GetReal("problem","lat_l3") : 0.0;
    const Real gvr = pin->DoesParameterExist("problem","gas_vr") ?
                     pin->GetReal("problem","gas_vr") : 0.0;
    // m1-sp-order2b (tests_m1/runs_5q_sporder2b), read only when named: gas_vth sets
    // v_theta = gas_vth sin(pi (x2 - x2min)/(x2max - x2min)) (sp: a lateral flow, zero at
    // reflecting theta walls); gas_teq = true puts the gas at T = (E/arad)^(1/4) of the
    // initial E (radiative equilibrium, the start of the stiff coupled tests)
    const Real gvt = pin->DoesParameterExist("problem","gas_vth") ?
                     pin->GetReal("problem","gas_vth") : 0.0;
    // STAGE CS3 (read only when named): gas_rot = Omega, a solid-body rotation of the
    // gas about z, v = Omega z x r, on sp (v_phi = Omega r sin theta) and on the cubed
    // sphere (covariant v.e_xi, v.e_eta); the moving-gas O(v/c) terms of the order gate
    const Real grot = pin->DoesParameterExist("problem","gas_rot") ?
                      pin->GetReal("problem","gas_rot") : 0.0;
    // (read only when named): gas_rot_axis = "ax,ay,az", the rotation axis of
    // gas_rot (normalised; default z)
    bool gax = false;
    Real gax_x = 0.0, gay = 0.0, gaz = 1.0;
    if (pin->DoesParameterExist("problem","gas_rot_axis")) {
      gax = true;
      std::string sa = pin->GetString("problem","gas_rot_axis");
      for (auto &ch : sa) {if (ch == ',') ch = ' ';}
      std::istringstream iss(sa);
      iss >> gax_x >> gay >> gaz;
      const Real an = std::sqrt(gax_x*gax_x + gay*gay + gaz*gaz);
      gax_x /= an; gay /= an; gaz /= an;
    }
    const bool gteq = pin->DoesParameterExist("problem","gas_teq") ?
                      pin->GetBoolean("problem","gas_teq") : false;
    // (read only when named): rad_eq_boost = true puts the radiation at the
    // exact fixed point of the implicit source in the gas_rot flow, with D = I/3:
    // F0 = 0 -> F = (4/3) v E (covariant on cs), E0 = E (1 - (4/3) beta^2) = arad T^4
    const bool reqb = pin->DoesParameterExist("problem","rad_eq_boost") &&
                      pin->GetBoolean("problem","rad_eq_boost");
    auto ccl_ = pmbp->pcoord->cos_cell;
    const Real ara = (pmbp->pradvet != nullptr) ? pmbp->pradvet->arad : 1.0;
    Real gm1 = ((pmbp->phydro != nullptr) ? pmbp->phydro->peos : pmbp->pmhd->peos)
               ->eos_data.gamma - 1.0;
    // (the history / final functions are enrolled before the restart return)
    // STAGE CS0 (m1-cs-implicit): on the cubed sphere x1 is r as well
    const bool spr = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
    // STAGE CS1 (cubed sphere), read only when named: lat_l = l multiplies E by
    // 1 + lat_amp P_l(rhat . axis), lat_axis = z | edge (1,1,0) | vertex (1,1,1)
    const bool csm = pmy_mesh_->use_cubed_sphere;
    int csl = 0;
    if (csm && pin->DoesParameterExist("problem","lat_l")) {
      csl = pin->GetInteger("problem","lat_l");
      std::string ax = pin->GetOrAddString("problem","lat_axis","z");
      Real a0 = 0.0, a1 = 0.0, a2 = 1.0;
      if (ax == "edge") {a0 = 1.0; a1 = 1.0; a2 = 0.0;}
      if (ax == "vertex") {a0 = 1.0; a1 = 1.0; a2 = 1.0;}
      const Real an = std::sqrt(a0*a0 + a1*a1 + a2*a2);
      vet_cs_ax[0] = a0/an; vet_cs_ax[1] = a1/an; vet_cs_ax[2] = a2/an;
    }
    // (read only when named): problem/am_hist = true takes the user history to the
    // double-precision energy / angular-momentum budget RadVetAmHistory
    // (cubed sphere and spherical polar)
    if (spr && pin->DoesParameterExist("problem","am_hist") &&
        pin->GetBoolean("problem","am_hist")) {
      user_hist_func = RadVetAmHistory;
      if (pin->DoesParameterExist("problem","am_hist_ebnd") &&
          pin->GetBoolean("problem","am_hist_ebnd") && pmbp->pradvet != nullptr) {
        vet_am_ebnd = true;
        pmbp->pradvet->ebnd_on = true;
      }
    }
    if (restart) return;
    auto uh = (pmbp->phydro != nullptr) ? pmbp->phydro->u0 : pmbp->pmhd->u0;
    auto x1v = pmbp->pcoord->x1v;
    auto x2v = pmbp->pcoord->x2v;
    auto x3v = pmbp->pcoord->x3v;
    const bool sp = pmy_mesh_->use_spherical_polar;
    // (MHD, cubed sphere, read only when named): b_azi = b, b_azi_axis = "ax,ay,az": the
    // field B = b a x r (azimuthal about the inclined axis a: B.rhat = 0 on every face,
    // so the reflecting shell is a CLOSED MHD system).  Faces: the radial ones carry 0, a
    // planar xi / eta face the exact flux b (a x int_f x dA).nhat_f divided by the code's
    // face area (b0.x2f/x3f are face-normal PROJECTIONS B.nhat), so the CT divergence is
    // round-off; the cell energy takes |B(rc)|^2/2.
    Real bazi = 0.0, bax = 0.0, bay = 0.0, baz = 1.0;
    if (pmbp->pmhd != nullptr && pin->DoesParameterExist("problem","b_azi")) {
      bazi = pin->GetReal("problem","b_azi");
      if (pin->DoesParameterExist("problem","b_azi_axis")) {
        std::string sa = pin->GetString("problem","b_azi_axis");
        for (auto &ch : sa) {if (ch == ',') ch = ' ';}
        std::istringstream iss(sa);
        iss >> bax >> bay >> baz;
      }
      const Real an = std::sqrt(bax*bax + bay*bay + baz*baz);
      bax /= an; bay /= an; baz /= an;
      // spherical polar (the sp reference of the cs MHD tests): the axis must be z there,
      // B = b r sin(theta) phihat; the phi faces carry the exact flux
      // b (r_r^3 - r_l^3)/3 (cos theta_l - cos theta_r) / area, the others 0
      const bool spz = pmy_mesh_->use_spherical_polar && std::fabs(bax) < 1.0e-14 &&
                       std::fabs(bay) < 1.0e-14 && baz > 0.0;
      if (!pmy_mesh_->use_cubed_sphere && !spz) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "problem/b_azi is implemented on the cubed sphere and (axis z) "
          << "on spherical polar only" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (spz) {
        auto b0 = pmbp->pmhd->b0;
        auto a3 = pmbp->pcoord->area.x3f;
        auto xf1 = pmbp->pcoord->xx1f;
        auto xf2 = pmbp->pcoord->xx2f;
        const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
        const Real bb = bazi;
        par_for("vet_sh_bazi_sp", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
          if (j <= je && k <= ke) {b0.x1f(m,k,j,i) = 0.0;}
          if (i > ie) return;
          if (k <= ke) {b0.x2f(m,k,j,i) = 0.0;}
          if (j <= je) {
            const Real rl = xf1(m,i), rr = xf1(m,i+1);
            const Real fl = bb*(rr*rr*rr - rl*rl*rl)/3.0*
                            (cos(xf2(m,j)) - cos(xf2(m,j+1)));
            b0.x3f(m,k,j,i) = fl/a3(m,k,j,i);
          }
        });
      }
    }
    if (pmbp->pmhd != nullptr && bazi != 0.0 && pmy_mesh_->use_cubed_sphere) {
      auto b0 = pmbp->pmhd->b0;
      auto a2 = pmbp->pcoord->area.x2f;
      auto a3 = pmbp->pcoord->area.x3f;
      auto xf1 = pmbp->pcoord->xx1f;
      auto pnl = pmbp->pmb->mb_panel;
      const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
      const int nx2 = indcs.nx2, nx3 = indcs.nx3;
      const Real bb = bazi, ax0 = bax, ax1 = bay, ax2 = baz;
      par_for("vet_sh_bazi", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const Real x2l = size.d_view(m).x2min, x2u = size.d_view(m).x2max;
        const Real x3l = size.d_view(m).x3min, x3u = size.d_view(m).x3max;
        const int p = pnl.d_view(m);
        if (j <= je && k <= ke) {b0.x1f(m,k,j,i) = 0.0;}
        if (i > ie) return;
        const Real rl = xf1(m,i), rr = xf1(m,i+1);
        const Real mom = (rr*rr*rr - rl*rl*rl)/3.0;
        Real fa[3], fb[3], fn[3];
        cubed_sphere::PanelFrame(p, fa, fb, fn);
        auto flux = [&](const Real xa, const Real ea, const Real xb, const Real eb,
                        const Real *nf) {
          Real qa[3], qb[3];
          cubed_sphere::PanelToCart(p, xa, ea, qa);
          cubed_sphere::PanelToCart(p, xb, eb, qb);
          const Real cx = qa[1]*qb[2] - qa[2]*qb[1], cy = qa[2]*qb[0] - qa[0]*qb[2];
          const Real cz = qa[0]*qb[1] - qa[1]*qb[0];
          const Real dt = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2];
          const Real th = std::sqrt(cx*cx + cy*cy + cz*cz)/(1.0 + dt);
          const Real mv[3] = {mom*th*(qa[0] + qb[0]), mom*th*(qa[1] + qb[1]),
                              mom*th*(qa[2] + qb[2])};
          // b (a x M).n
          const Real w[3] = {ax1*mv[2] - ax2*mv[1], ax2*mv[0] - ax0*mv[2],
                             ax0*mv[1] - ax1*mv[0]};
          return bb*(w[0]*nf[0] + w[1]*nf[1] + w[2]*nf[2]);
        };
        if (k <= ke) {   // x2 face j: xi = xf
          const Real xf = 0.25*M_PI*LeftEdgeX(j-js, nx2, x2l, x2u);
          const Real el = 0.25*M_PI*LeftEdgeX(k-ks, nx3, x3l, x3u);
          const Real er = 0.25*M_PI*LeftEdgeX(k+1-ks, nx3, x3l, x3u);
          const Real t = std::tan(xf);
          Real nf[3] = {fa[0] - t*fn[0], fa[1] - t*fn[1], fa[2] - t*fn[2]};
          const Real nn = 1.0/std::sqrt(nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2]);
          for (int q = 0; q < 3; ++q) {nf[q] *= nn;}
          b0.x2f(m,k,j,i) = flux(xf, el, xf, er, nf)/a2(m,k,j,i);
        }
        if (j <= je) {   // x3 face k: eta = ef
          const Real ef = 0.25*M_PI*LeftEdgeX(k-ks, nx3, x3l, x3u);
          const Real xl = 0.25*M_PI*LeftEdgeX(j-js, nx2, x2l, x2u);
          const Real xr = 0.25*M_PI*LeftEdgeX(j+1-js, nx2, x2l, x2u);
          const Real t = std::tan(ef);
          Real nf[3] = {fb[0] - t*fn[0], fb[1] - t*fn[1], fb[2] - t*fn[2]};
          const Real nn = 1.0/std::sqrt(nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2]);
          for (int q = 0; q < 3; ++q) {nf[q] *= nn;}
          b0.x3f(m,k,j,i) = flux(xl, ef, xr, ef, nf)/a3(m,k,j,i);
        }
      });
    }
    const Real cax = vet_cs_ax[0], cay = vet_cs_ax[1], caz = vet_cs_ax[2];
    // (cubed sphere, read only when named): rho_lat_amp = a multiplies the gas
    // density by 1 + a P_1(rhat . axis) with the axis along a cube edge, (1,1,0)/sqrt(2):
    // a smooth lateral opacity structure for the seam / vertex order gate (a static cold
    // gas: the pressure stays negligible)
    Real rla = 0.0;
    const Real rlx = 1.0/std::sqrt(2.0), rly = rlx, rlz = 0.0;
    if (csm && pin->DoesParameterExist("problem","rho_lat_amp")) {
      rla = pin->GetReal("problem","rho_lat_amp");
    }
    // (spherical polar, read only when named): rho_lat_sp = a multiplies the
    // gas density by 1 + a cos(2 pi (theta - x2min)/L2) cos(2 pi (phi - x3min)/L3): the
    // periodic-wedge counterpart of rho_lat_amp (S7 steady thin test, cs vs sp)
    Real rls = 0.0;
    if (sp && pin->DoesParameterExist("problem","rho_lat_sp")) {
      rls = pin->GetReal("problem","rho_lat_sp");
    }
    // (cs and sp, read only when named): rho_rad_h = H > 0 multiplies the gas
    // density by exp(-(r - x1min)/H): a smooth radial atmosphere (thick at the base,
    // thin at the top, the photosphere resolved) for the steady S7 thin-regime order test
    Real rrh = 0.0;
    if (pin->DoesParameterExist("problem","rho_rad_h")) {
      rrh = pin->GetReal("problem","rho_rad_h");
    }
    auto mbpan = pmbp->pmb->mb_panel;
    auto &msz = pmy_mesh_->mesh_size;
    const Real x1a = msz.x1min, x1l = msz.x1max - msz.x1min;
    const Real x2a = msz.x2min, x2l = msz.x2max - msz.x2min;
    const Real x3a = msz.x3min, x3l = msz.x3max - msz.x3min;
    par_for("vet_sph_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real r = spr ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                            size.d_view(m).x1max);
      uh(m,IDN,k,j,i) = dgas;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = dgas*tgas/gm1;
      Real ebz = 0.0;   // |B(rc)|^2/2 of b_azi, re-added after gas_teq below
      if (bazi != 0.0 && sp) {   // sp: B = b r sin(theta) phihat at the cell centre
        const Real st = sin(x2v(m,j));
        ebz = 0.5*bazi*bazi*r*r*st*st;
        uh(m,IEN,k,j,i) += ebz;
      }
      if (bazi != 0.0 && !sp) {   // |B(rc)|^2/2, B = b a x rc
        const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                               size.d_view(m).x2max);
        const Real eta = 0.25*M_PI*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                                size.d_view(m).x3max);
        Real q[3];
        cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, q);
        const Real w[3] = {bay*q[2] - baz*q[1], baz*q[0] - bax*q[2], bax*q[1] - bay*q[0]};
        ebz = 0.5*bazi*bazi*r*r*(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        uh(m,IEN,k,j,i) += ebz;
      }
      Real x = (r - r0)/wid;
      u0(m,radvet::VET_E,k,j,i) = fmax(eout + eamp*exp(-x*x), efl);
      if (lamp != 0.0 && sp) {
        u0(m,radvet::VET_E,k,j,i) *= 1.0 + lamp*cos(ll2*M_PI*(x2v(m,j) - x2a)/x2l)
                                   *cos(2.0*M_PI*ll3*(x3v(m,k) - x3a)/x3l);
      }
      if (lamp != 0.0 && csm && csl > 0) {
        const Real mu = VetCsMu(mbpan.d_view(m), k, j, ks, js, nx2, nx3,
                               size.d_view(m).x2min, size.d_view(m).x2max,
                               size.d_view(m).x3min, size.d_view(m).x3max, cax, cay, caz);
        u0(m,radvet::VET_E,k,j,i) *= 1.0 + lamp*VetLegendre(csl, mu);
      }
      if (rla != 0.0) {
        const Real mr = VetCsMu(mbpan.d_view(m), k, j, ks, js, nx2, nx3,
                               size.d_view(m).x2min, size.d_view(m).x2max,
                               size.d_view(m).x3min, size.d_view(m).x3max, rlx, rly, rlz);
        const Real fr = 1.0 + rla*VetLegendre(1, mr);
        uh(m,IDN,k,j,i) *= fr;
        uh(m,IEN,k,j,i) *= fr;
      }
      if (rls != 0.0) {
        const Real fr = 1.0 + rls*cos(2.0*M_PI*(x2v(m,j) - x2a)/x2l)
                        *cos(2.0*M_PI*(x3v(m,k) - x3a)/x3l);
        uh(m,IDN,k,j,i) *= fr;
        uh(m,IEN,k,j,i) *= fr;
      }
      if (rrh > 0.0) {
        const Real fr = exp(-(r - x1a)/rrh);
        uh(m,IDN,k,j,i) *= fr;
        uh(m,IEN,k,j,i) *= fr;
      }
      if (gteq) {
        // CSMHD: the gas in equilibrium with E keeps the field's |B|^2/2 (b_azi) in IEN
        uh(m,IEN,k,j,i) = dgas*sqrt(sqrt(u0(m,radvet::VET_E,k,j,i)/ara))/gm1 + ebz;
      }
      if (gvr != 0.0) {
        Real vr = gvr*sin(M_PI*(r - x1a)/x1l);
        uh(m,IM1,k,j,i) = dgas*vr;
        uh(m,IEN,k,j,i) += 0.5*dgas*vr*vr;
      }
      if (gvt != 0.0 && sp) {
        Real vt = gvt*sin(M_PI*(x2v(m,j) - x2a)/x2l);
        uh(m,IM2,k,j,i) = dgas*vt;
        uh(m,IEN,k,j,i) += 0.5*dgas*vt*vt;
      }
      if (grot != 0.0 && sp && !gax) {
        const Real vp = grot*r*sin(x2v(m,j));
        uh(m,IM3,k,j,i) = dgas*vp;
        uh(m,IEN,k,j,i) += 0.5*dgas*vp*vp;
      }
      if (grot != 0.0 && sp && gax) {
        // gas_rot_axis: v = Omega a x r on the orthonormal (theta, phi) basis
        const Real th = x2v(m,j), ph = x3v(m,k);
        const Real q[3] = {sin(th)*cos(ph), sin(th)*sin(ph), cos(th)};
        const Real et[3] = {cos(th)*cos(ph), cos(th)*sin(ph), -sin(th)};
        const Real ep[3] = {-sin(ph), cos(ph), 0.0};
        const Real w[3] = {grot*r*(gay*q[2] - gaz*q[1]), grot*r*(gaz*q[0] - gax_x*q[2]),
                           grot*r*(gax_x*q[1] - gay*q[0])};
        const Real vt = w[0]*et[0] + w[1]*et[1] + w[2]*et[2];
        const Real vp = w[0]*ep[0] + w[1]*ep[1] + w[2]*ep[2];
        uh(m,IM2,k,j,i) = dgas*vt;
        uh(m,IM3,k,j,i) = dgas*vp;
        uh(m,IEN,k,j,i) += 0.5*dgas*(vt*vt + vp*vp);
      }
      if (grot != 0.0 && csm) {
        const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                               size.d_view(m).x2max);
        const Real eta = 0.25*M_PI*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                                size.d_view(m).x3max);
        Real q[3], ta[3], tb[3];
        cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, q);
        cubed_sphere::PanelTangents(mbpan.d_view(m), xi, eta, ta, tb);
        Real w[3] = {-grot*r*q[1], grot*r*q[0], 0.0};
        if (gax) {   // gas_rot_axis: v = Omega a x r
          w[0] = grot*r*(gay*q[2] - gaz*q[1]);
          w[1] = grot*r*(gaz*q[0] - gax_x*q[2]);
          w[2] = grot*r*(gax_x*q[1] - gay*q[0]);
        }
        uh(m,IM2,k,j,i) = dgas*(w[0]*ta[0] + w[1]*ta[1] + w[2]*ta[2]);
        uh(m,IM3,k,j,i) = dgas*(w[0]*tb[0] + w[1]*tb[1] + w[2]*tb[2]);
        uh(m,IEN,k,j,i) += 0.5*dgas*(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
      }
      u0(m,radvet::VET_F1,k,j,i) = 0.0;
      if (ein > 0.0 && i <= is) {
        u0(m,radvet::VET_E,k,j,i) = ein;
        u0(m,radvet::VET_F1,k,j,i) = cl*ein;
      }
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
      if (reqb) {
        // the moving radiative equilibrium of the cell velocity (gas_rot,
        // gas_vr, gas_vth), the exact fixed point F0 = 0, E0 = arad T^4 for D = I/3
        // (chi = 1/3): f = beta (1 + chi), E = arad T^4 /
        // ((1 + beta^2) - 2 beta f + beta^2 chi), F = c E f vhat (cs: covariant)
        const Real dd = uh(m,IDN,k,j,i);
        const Real w1 = uh(m,IM1,k,j,i)/dd, w2 = uh(m,IM2,k,j,i)/dd;
        const Real w3 = uh(m,IM3,k,j,i)/dd;
        Real vq = w1*w1 + w2*w2 + w3*w3;
        if (csm) {
          const Real c = ccl_(m,k,j);
          vq = w1*w1 + (w2*w2 + w3*w3 - 2.0*c*w2*w3)/(1.0 - c*c);
        }
        const Real bb = sqrt(vq)/cl;
        const Real ff = (4.0/3.0)*bb, chq = 1.0/3.0;
        const Real er = ara*SQR(SQR(tgas))/((1.0 + bb*bb) - 2.0*bb*ff + bb*bb*chq);
        const Real fs = (bb > 0.0) ? cl*er*ff/(bb*cl) : 0.0;   // F = fs v
        u0(m,radvet::VET_E,k,j,i) = er;
        u0(m,radvet::VET_F1,k,j,i) = fs*w1;
        u0(m,radvet::VET_F2,k,j,i) = fs*w2;
        u0(m,radvet::VET_F3,k,j,i) = fs*w3;
      }
    });
    // (read only when named): beam_f = f > 0 makes the shell a radially
    // BEAMED pulse, F_r = f c E on the cells and the face-normal comoving flux of the
    // implicit solve f0x1 = f c E(r_f) on every x1 face (the gas is static), so that
    // the face-eliminated central form starts from the same beam its memory term
    // carries.  A shell with F = 0 instead splits into an ingoing and an outgoing half.
    if (pin->DoesParameterExist("problem","beam_f") && pmbp->pradvet != nullptr &&
        pmbp->pradvet->f0x1.extent_int(0) > 0) {
      const Real bf = pin->GetReal("problem","beam_f");
      auto f0 = pmbp->pradvet->f0x1;
      auto xf = pmbp->pcoord->xx1f;
      par_for("vet_sph_beam", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        u0(m,radvet::VET_F1,k,j,i) = bf*cl*u0(m,radvet::VET_E,k,j,i);
        if (i >= is && i <= is+nx1) {
          Real rf = spr ? xf(m,i) : LeftEdgeX(i-is, nx1, size.d_view(m).x1min,
                                              size.d_view(m).x1max);
          Real x = (rf - r0)/wid;
          f0(m,k,j,i) = bf*cl*fmax(eout + eamp*exp(-x*x), efl);
        }
      });
    }
  } else if (test.compare("sph_atm") == 0) {
    // STAGE S2 (tests_m1/runs_5b_sp_s2, gate T-S4): an extended grey atmosphere on the
    // spherical-polar wedge.  Gas rho = atm_rho0 (r/r_in)^(-atm_rho_n) at T = atm_temp,
    // E = e_out, F = 0; the luminosity enters through implicit_bc_x1min = flux and leaves
    // through a Marshak outer face.  The gas evolves (risk 6).  force_reference = wb_arad
    // takes the reference acceleration of the exact steady state, F = F_in (r_in/r)^2
    // (what a well-balanced hydro source would carry), see atm_aref below.
    if (pmbp->phydro == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<problem>/vet_test = sph_atm needs a <hydro> block" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    vet_sa_rho0 = pin->GetOrAddReal("problem","atm_rho0",1.0);
    vet_sa_n = pin->GetOrAddReal("problem","atm_rho_n",0.0);
    vet_sa_rin = pmy_mesh_->mesh_size.x1min;
    if (pin->DoesParameterExist("problem","atm_scale_h")) {
      vet_sa_h = pin->GetReal("problem","atm_scale_h");
    }
    const Real sah = vet_sa_h;
    Real tgas = pin->GetOrAddReal("problem","atm_temp",1.0);
    Real eout = pin->GetOrAddReal("problem","e_out",1.0e-10);
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    auto x1v = pmbp->pcoord->x1v;
    // STAGE CS1: on the cubed sphere x1 is r as well
    const bool sp = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
    const bool csm = pmy_mesh_->use_cubed_sphere;
    if (pmbp->pradvet->force_ref == radvet::VET_FREF_WB_ARAD) {
      Real fin = pin->GetOrAddReal("rad_vet","implicit_flux_x1min",0.0);
      Real kt = pmbp->pradvet->kappa_f + pmbp->pradvet->kappa_s;
      Real rin = vet_sa_rin;
      // a local View: the module's copy (arad_ref) keeps it alive and is destroyed with
      // the module, before Kokkos::finalize
      DvceArray4D<Real> aref("vet_sa_aref", nmb1+1, n3, n2, n1);
      // atm_aref = face (default): the reference in the form the implicit coupling
      // deposits the force, the mean of kappa_t F_f/c over the two x1 faces with the
      // exact steady F_f = F_in (r_in/r_f)^2; = cell: kappa_t F(r_c)/c at the centroid
      const bool afc = (pin->GetOrAddString("problem","atm_aref","face")
                        .compare("cell") == 0);
      auto xf = pmbp->pcoord->xx1f;
      par_for("vet_sa_aref", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real r = sp ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                             size.d_view(m).x1max);
        Real rl = sp ? xf(m,i) : LeftEdgeX(i-is, nx1, size.d_view(m).x1min,
                                           size.d_view(m).x1max);
        Real rr = sp ? xf(m,i+1) : LeftEdgeX(i+1-is, nx1, size.d_view(m).x1min,
                                             size.d_view(m).x1max);
        aref(m,k,j,i) = afc ? (kt*fin*rin*rin/(cl*r*r))
                            : (0.5*kt*fin*rin*rin*(1.0/(rl*rl) + 1.0/(rr*rr))/cl);
      });
      pmbp->pradvet->SetForceReference(aref);
    }
    // atm_init = eddington: start from the spherical EDDINGTON solution of the same
    // atmosphere, E = 3 int_r^rout rho kappa_t F/c dr' + F_out/(c q), F = F_in
    // (r_in/r)^2 on the cells and the x1 faces, the gas at T = (E/a)^(1/4) when it
    // absorbs (a start close to the answer, as the Cartesian atmosphere gate of runs_3a2
    // does)
    const bool edd0 = (pin->GetOrAddString("problem","atm_init","vacuum")
                       .compare("eddington") == 0);
    // atm_seed: a small deterministic cell-to-cell perturbation of the initial E (a seed
    // for transverse modes, which a symmetric start on a Cartesian mesh never has)
    const Real seed = pin->GetOrAddReal("problem","atm_seed",0.0);
    // m1-sp-order2: the smooth lateral mode of sph_shell (lat_amp, lat_l2, lat_l3; read
    // only when named, sp only)
    const Real lamp = pin->DoesParameterExist("problem","lat_amp") ?
                      pin->GetReal("problem","lat_amp") : 0.0;
    const Real ll2 = pin->DoesParameterExist("problem","lat_l2") ?
                     pin->GetReal("problem","lat_l2") : 1.0;
    const Real ll3 = pin->DoesParameterExist("problem","lat_l3") ?
                     pin->GetReal("problem","lat_l3") : 0.0;
    auto x2v = pmbp->pcoord->x2v;
    auto x3v = pmbp->pcoord->x3v;
    // STAGE CS1 (cubed sphere), read only when named: the P_l(rhat . axis) mode of
    // sph_shell (lat_l, lat_axis, lat_amp)
    int csl = 0;
    if (csm && pin->DoesParameterExist("problem","lat_l")) {
      csl = pin->GetInteger("problem","lat_l");
      std::string ax = pin->GetOrAddString("problem","lat_axis","z");
      Real a0 = 0.0, a1 = 0.0, a2 = 1.0;
      if (ax == "edge") {a0 = 1.0; a1 = 1.0; a2 = 0.0;}
      if (ax == "vertex") {a0 = 1.0; a1 = 1.0; a2 = 1.0;}
      const Real an = std::sqrt(a0*a0 + a1*a1 + a2*a2);
      vet_cs_ax[0] = a0/an; vet_cs_ax[1] = a1/an; vet_cs_ax[2] = a2/an;
    }
    const Real cax = vet_cs_ax[0], cay = vet_cs_ax[1], caz = vet_cs_ax[2];
    auto mbpan = pmbp->pmb->mb_panel;
    const Real x2a = pmy_mesh_->mesh_size.x2min;
    const Real x2l = pmy_mesh_->mesh_size.x2max - x2a;
    const Real x3a = pmy_mesh_->mesh_size.x3min;
    const Real x3l = pmy_mesh_->mesh_size.x3max - x3a;
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    Real rho0 = vet_sa_rho0, rin = vet_sa_rin, en = vet_sa_n;
    Real fin = pin->GetOrAddReal("rad_vet","implicit_flux_x1min",0.0);
    Real kt = pmbp->pradvet->kappa_f + pmbp->pradvet->kappa_s;
    Real rout = pmy_mesh_->mesh_size.x1max;
    Real qm = pmbp->pradvet->marshak_q;
    // with absorption the gas starts at the radiation temperature; pure scattering keeps
    // atm_temp (uniform: the static state of risk 6)
    const bool teq = (pmbp->pradvet->kappa_p > 0.0);
    par_for("vet_sa_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real r = sp ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                           size.d_view(m).x1max);
      Real d = (sah > 0.0) ? rho0*exp(-(r - rin)/sah) : rho0*pow(r/rin, -en);
      Real e = fmax(eout, efl), f = 0.0, tt = tgas;
      if (edd0) {
        f = fin*rin*rin/(r*r);
        e = 3.0*rho0*kt*fin*pow(rin, en + 2.0)/(cl*(en + 1.0))
            *(pow(r, -(en + 1.0)) - pow(rout, -(en + 1.0)))
            + fin*rin*rin/(rout*rout*cl*qm);
        if (teq) {tt = pow(e/ar, 0.25);}
      }
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*tt/gm1;
      u0(m,radvet::VET_E,k,j,i) = e*(1.0 + seed*cos(2.3*j + 1.7*k + 0.9*i));
      if (lamp != 0.0 && sp && !csm) {
        u0(m,radvet::VET_E,k,j,i) *= 1.0 + lamp*cos(ll2*M_PI*(x2v(m,j) - x2a)/x2l)
                                   *cos(2.0*M_PI*ll3*(x3v(m,k) - x3a)/x3l);
      }
      if (lamp != 0.0 && csm && csl > 0) {
        const Real mu = VetCsMu(mbpan.d_view(m), k, j, ks, js, nx2, nx3,
                               size.d_view(m).x2min, size.d_view(m).x2max,
                               size.d_view(m).x3min, size.d_view(m).x3max, cax, cay, caz);
        u0(m,radvet::VET_E,k,j,i) *= 1.0 + lamp*VetLegendre(csl, mu);
      }
      u0(m,radvet::VET_F1,k,j,i) = f;
      u0(m,radvet::VET_F2,k,j,i) = 0.0;
      u0(m,radvet::VET_F3,k,j,i) = 0.0;
    });
    // the comoving x1 face fluxes of the implicit transport (allocated only there)
    auto ff0 = pmbp->pradvet->f0x1;
    if (edd0 && ff0.extent_int(0) >= nmb1 + 1) {
      auto xf = pmbp->pcoord->xx1f;
      const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
      const int e1 = ff0.extent_int(3) - 1;
      par_for("vet_sa_icf", DevExeSpace(), 0,nmb1,0,e3,0,e2,0,e1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real r = sp ? xf(m,i) : LeftEdgeX(i-is, nx1, size.d_view(m).x1min,
                                          size.d_view(m).x1max);
        ff0(m,k,j,i) = fin*rin*rin/(r*r);
      });
    }
  } else if (test.compare("rhse") == 0) {
    RadiationVetRhse(pin, restart);   // rad_vet_rhse.cpp
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/vet_test = '" << test << "' not implemented "
      << "(beam | pulse1d | thick_pulse | tophat | equil | advect_pulse "
      << "| advect_uniform | marshak | shadow | radshock "
      << "| atmosphere | radwave | sph_shell | sph_atm | rhse)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadVetShadowGas()
//! \brief re-impose the PRESCRIBED static gas of the shadow test (ghost zones included)
//! at the end of every hydro stage, so that rho and T are exactly the analytic ones and
//! the opacity law is a pure function of position.  bdt is unused.

void RadVetShadowGas(Mesh *pm, const Real bdt) {
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
  Real rho0 = vet_sh2_rho0, ratio = vet_sh2_ratio, t0 = vet_sh2_t0;
  Real x0 = vet_sh2_x0, y0 = vet_sh2_y0, ax = vet_sh2_ax, ay = vet_sh2_ay;
  Real delta = vet_sh2_delta, gm1 = vet_sh2_gm1;
  const bool dsq = vet_sh2_dsq, peq = vet_sh2_peq;

  par_for("vet_shadow_reset", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real d = VetShadowRho(x1v, x2v, rho0, ratio, x0, y0, ax, ay, delta, dsq);
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = (peq ? rho0 : d)*t0/gm1;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadVetCylGas()
//! \brief vet_test = cyl: the prescribed static gas (ghosts included),
//! re-imposed at the end of every hydro stage: background (cyl_rho_bg, cyl_t_bg), and
//! inside cylinder c (cell centre within cylc_r of (cylc_x, cylc_y)) (cylc_rho, cylc_t).

void RadVetCylGas(Mesh *pm, const Real bdt) {
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
  const int nc = vet_cy_n;
  VetCylP p;
  for (int c = 0; c < 24; ++c) {
    for (int q = 0; q < 9; ++q) {p.v[c][q] = vet_cy_p[c][q];}
  }
  const Real rb = vet_cy_rhob, tb = vet_cy_tb, gm1 = vet_cy_gm1;
  par_for("vet_cyl_reset", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real xv = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    Real yv = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    Real d = rb, t = tb;
    for (int c = 0; c < nc; ++c) {
      bool in;
      if (p.v[c][5] > 0.5) {
        const Real ca = cos(p.v[c][8]), sa = sin(p.v[c][8]);
        const Real dx = xv - p.v[c][0], dy = yv - p.v[c][1];
        in = (fabs(dx*ca + dy*sa) <= p.v[c][6]) && (fabs(-dx*sa + dy*ca) <= p.v[c][7]);
      } else {
        in = (SQR(xv - p.v[c][0]) + SQR(yv - p.v[c][1]) <= p.v[c][2]*p.v[c][2]);
      }
      if (in) {
        d = p.v[c][3];
        t = p.v[c][4];
      }
    }
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = d*t/gm1;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadVetShockBC()
//! \brief Dirichlet at the two far-field states of the radiative shock, for BOTH the
//! radiation moments and the gas: the shock is held in its own frame.

void RadVetShockBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradvet == nullptr || pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradvet->u0;
  auto uh = pmbp->phydro->u0;
  Real cl = pmbp->pradvet->c_light;
  Real dl = vet_rs_l[0], ml = vet_rs_l[1], el = vet_rs_l[2], rl = vet_rs_l[3];
  Real dr = vet_rs_r[0], mr = vet_rs_r[1], er = vet_rs_r[2], rr = vet_rs_r[3];

  par_for("vet_shock_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        uh(m,IDN,k,j,is-i-1) = dl;
        uh(m,IM1,k,j,is-i-1) = ml;
        uh(m,IM2,k,j,is-i-1) = 0.0;
        uh(m,IM3,k,j,is-i-1) = 0.0;
        uh(m,IEN,k,j,is-i-1) = el;
        u0(m,radvet::VET_E,k,j,is-i-1) = rl;
        u0(m,radvet::VET_F1,k,j,is-i-1) = (4.0/3.0)*(ml/dl)*rl;
        u0(m,radvet::VET_F2,k,j,is-i-1) = 0.0;
        u0(m,radvet::VET_F3,k,j,is-i-1) = 0.0;
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        uh(m,IDN,k,j,ie+i+1) = dr;
        uh(m,IM1,k,j,ie+i+1) = mr;
        uh(m,IM2,k,j,ie+i+1) = 0.0;
        uh(m,IM3,k,j,ie+i+1) = 0.0;
        uh(m,IEN,k,j,ie+i+1) = er;
        u0(m,radvet::VET_E,k,j,ie+i+1) = rr;
        u0(m,radvet::VET_F1,k,j,ie+i+1) = (4.0/3.0)*(mr/dr)*rr;
        u0(m,radvet::VET_F2,k,j,ie+i+1) = 0.0;
        u0(m,radvet::VET_F3,k,j,ie+i+1) = 0.0;
      }
    }
  });
  (void) cl;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadVetAtmGas()
//! \brief re-impose the PRESCRIBED exponential atmosphere (ghost zones included) at the
//! end of every hydro stage, so that rho*kappa is an exact function of height for the
//! whole run and the steady radiation state is the closed-form M1 one.  bdt is unused.

void RadVetAtmGas(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int nx1 = indcs.nx1;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &size = pmbp->pmb->mb_size;
  auto uh = pmbp->phydro->u0;
  Real rht = vet_at_rho, hh = vet_at_h, ztop = vet_at_ztop, egas = vet_at_egas;
  (void) bdt;

  par_for("vet_atm_reset", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real d = rht*exp((ztop - x1v)/hh);
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = d*egas;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadVetAmHistory
//! \brief (sph_shell, problem/am_hist): double-precision budget sums over the
//! active cells, V = Coordinates volume, arm r_c rhat (cell-centre radius, as the
//! <hydro>/am_transport correction): Eg = sum V u(IEN), Er = sum V E, Lg = sum V r x m
//! (x, y, z), Lr = sum V r x F/c^2 (x, y, z).  cs: covariant m, F raised with the cell's
//! e_xi.e_eta (the clean scheme's F in its frame of y, CsFrame); sp: orthonormal.

void RadVetAmHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = vet_am_ebnd ? 15 : 14;
  const char *lab[15] = {"Eg", "Er", "Lgx", "Lgy", "Lgz", "Lrx", "Lry", "Lrz", "Sx", "Sy",
                         "Sz", "Lfx", "Lfy", "Lfz", "Eb"};
  for (int n = 0; n < pdata->nhist; ++n) {pdata->label[n] = lab[n];}
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradvet == nullptr || (pmbp->phydro == nullptr && pmbp->pmhd == nullptr)) {
    return;
  }
  auto u0 = pmbp->pradvet->u0;
  auto uh = (pmbp->phydro != nullptr) ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto ccell = pmbp->pcoord->cos_cell;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpan = pmbp->pmb->mb_panel;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const bool csm = pm->use_cubed_sphere;
  const Real ic2 = 1.0/(pmbp->pradvet->c_light*pmbp->pradvet->c_light);
  // (clean_stag_split: the radiation momentum is the faces' F^H, not the cell flux; on
  // the cubed sphere the cell flux is the face momenta's vector.  The clean scheme's cs
  // cell flux is in the frame of y, (r, theta, phi) about rhat, not covariant)
  radvet::RadiationVet *pr = pmbp->pradvet;
  const bool sfh = pr->cln_split && !csm;
  const bool fry = csm;
  const bool md = pm->multi_d, td = pm->three_d;
  auto g1 = pr->f0x1, g2 = pr->f0x2, g3 = pr->f0x3;
  // <hydro>/am_vector: the spin reservoir, L = am_spin rhat_c (total)
  hydro::AmModule *pam = hydro::AmModuleOf(pmbp);
  const bool spn = (pam != nullptr) && pam->am_vector;
  auto asp = spn ? pam->am_spin : DvceArray4D<Real>();
  // frame_omega: the frame part of the absolute L, sum V rho rc x (Omega x rc)
  const Real fo0 = (pam != nullptr) ? pam->frame_om[0] : 0.0;
  const Real fo1 = (pam != nullptr) ? pam->frame_om[1] : 0.0;
  const Real fo2 = (pam != nullptr) ? pam->frame_om[2] : 0.0;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("VetAmHist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks; j += js; i += is;
    array_sum::GlobalSum hv;
    for (int n = 0; n < NHISTORY_VARIABLES; ++n) {hv.the_array[n] = 0.0;}
    Real rh[3], ta[3], tb[3];
    Real c = 0.0;
    if (csm) {
      const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                             size.d_view(m).x2max);
      const Real eta = 0.25*M_PI*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                              size.d_view(m).x3max);
      cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, rh);
      cubed_sphere::PanelTangents(mbpan.d_view(m), xi, eta, ta, tb);
      c = ccell(m,k,j);
    } else {
      const Real th = x2v(m,j), ph = x3v(m,k);
      rh[0] = sin(th)*cos(ph); rh[1] = sin(th)*sin(ph); rh[2] = cos(th);
      ta[0] = cos(th)*cos(ph); ta[1] = cos(th)*sin(ph); ta[2] = -sin(th);
      tb[0] = -sin(ph); tb[1] = cos(ph); tb[2] = 0.0;
    }
    const Real r = x1v(m,i), v = vol(m,k,j,i), s2 = 1.0 - c*c;
    auto lv = [&](const Real a1, const Real a2c, const Real a3c, Real *l) {
      const Real a2 = (a2c - c*a3c)/s2, a3 = (a3c - c*a2c)/s2;   // contravariant
      Real p[3];
      for (int q = 0; q < 3; ++q) {p[q] = a1*rh[q] + a2*ta[q] + a3*tb[q];}
      l[0] = r*(rh[1]*p[2] - rh[2]*p[1]);
      l[1] = r*(rh[2]*p[0] - rh[0]*p[2]);
      l[2] = r*(rh[0]*p[1] - rh[1]*p[0]);
    };
    Real lg[3], lr[3];
    lv(uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i), lg);
    if (fry) {   // (F = F_r rhat + F_theta theta_hat + F_phi phi_hat)
      const Real rho = std::sqrt(rh[0]*rh[0] + rh[1]*rh[1]);
      const Real ph[3] = {-rh[1]/rho, rh[0]/rho, 0.0};
      const Real th[3] = {ph[1]*rh[2], -ph[0]*rh[2], ph[0]*rh[1] - ph[1]*rh[0]};
      Real p[3];
      for (int q = 0; q < 3; ++q) {
        p[q] = u0(m,radvet::VET_F1,k,j,i)*rh[q] + u0(m,radvet::VET_F2,k,j,i)*th[q] +
               u0(m,radvet::VET_F3,k,j,i)*ph[q];
      }
      lr[0] = r*(rh[1]*p[2] - rh[2]*p[1]);
      lr[1] = r*(rh[2]*p[0] - rh[0]*p[2]);
      lr[2] = r*(rh[0]*p[1] - rh[1]*p[0]);
    } else if (sfh) {
      lv(0.5*(g1(m,k,j,i) + g1(m,k,j,i+1)), md ? 0.5*(g2(m,k,j,i) + g2(m,k,j+1,i)) : 0.0,
         td ? 0.5*(g3(m,k,j,i) + g3(m,k+1,j,i)) : 0.0, lr);
    } else {
      lv(u0(m,radvet::VET_F1,k,j,i), u0(m,radvet::VET_F2,k,j,i),
         u0(m,radvet::VET_F3,k,j,i),
         lr);
    }
    hv.the_array[0] = v*uh(m,IEN,k,j,i);
    hv.the_array[1] = v*u0(m,radvet::VET_E,k,j,i);
    for (int q = 0; q < 3; ++q) {
      hv.the_array[2+q] = v*lg[q];
      hv.the_array[5+q] = v*lr[q]*ic2;
      hv.the_array[8+q] = spn ? asp(m,k,j,i)*rh[q] : 0.0;
    }
    {
      const Real x[3] = {r*rh[0], r*rh[1], r*rh[2]};
      const Real ox[3] = {fo1*x[2] - fo2*x[1], fo2*x[0] - fo0*x[2], fo0*x[1] - fo1*x[0]};
      const Real d = v*uh(m,IDN,k,j,i);
      hv.the_array[11] = d*(x[1]*ox[2] - x[2]*ox[1]);
      hv.the_array[12] = d*(x[2]*ox[0] - x[0]*ox[2]);
      hv.the_array[13] = d*(x[0]*ox[1] - x[1]*ox[0]);
    }
    mb_sum += hv;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  for (int n = 0; n < 14; ++n) {pdata->hdata[n] = sum_this_mb.the_array[n];}
  // ebnd_acc is global on every rank; the history output sums over ranks
  if (vet_am_ebnd) {
    pdata->hdata[14] = (global_variable::my_rank == 0) ? pmbp->pradvet->ebnd_acc : 0.0;
  }
  // Lz_tot: gas + radiation + ledger L_z of the clean scheme (LzTotal; with
  // clean_cs_lz_tel the radiation's is that of its face momenta), last column
  const int nl = pdata->nhist;
  const Real lzt = pmbp->pradvet->LzTotal();   // (collective)
  pdata->label[nl] = "Lz_tot";
  pdata->hdata[nl] = (global_variable::my_rank == 0) ? lzt : 0.0;
  pdata->nhist = nl + 1;
}
