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
//!
//!   radshock  (design sect. 8, T7) Lowrie & Edwards (2008) steady radiative shock,
//!             initialised from the semi-analytic profile written by
//!             tests_m1/t7_radshock.py --write-ref, so that only a local relaxation is
//!             needed rather than the many flow-through times a step start would cost.
//!             With <problem>/m1_shock_ref unset the two far-field states are laid down
//!             with a tanh jump of <problem>/m1_shock_w cells instead.  Both x1 faces
//!             are Dirichlet at the far-field states (user BC), in the shock frame.
//!
//!   atmosphere (milestone 3a2, gate I6b) STATIC plane-parallel grey column with a
//!             constant imposed flux at the bottom and a free surface at the top.  The
//!             gas is PRESCRIBED (a user source term re-imposes it every stage and
//!             <rad_m1>/gas_feedback = false), the opacity is pure scattering
//!             (kappa_P = kappa_E = 0) and rho falls exponentially with height,
//!
//!               rho(z) = atm_rho_top exp[(x1max - z)/atm_scale_h],
//!               tau(z) = kappa rho_top H (exp[(x1max - z)/H] - 1)   (downward from the
//!                                                                    TOP MESH FACE),
//!
//!             so the steady state is the closed-form M1 moment solution of
//!             bench/m1_stage2/ic/build_ic.py: with F = const,
//!             dP_rad/dtau = F/c and P_rad = chi(f) E = (F/c)(chi/f) give
//!
//!               chi(f)/f = q0 + tau,      E(tau) = F/(c f(tau)),
//!
//!             on the DECREASING branch of chi/f (which has a minimum 0.89806 at
//!             f = 0.8250; q0 > that, so the solution never leaves it).  The constant q0
//!             is fixed by the surface: the Marshak condition F = c q E_top means
//!             f(0) = q = <rad_m1>/marshak_q exactly, hence the CONSISTENT PAIR is
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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
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

// radshock: the two Dirichlet end states (rho, v, E_gas, E_rad, F_rad)
Real m1_rs_l[5] = {1.0, 0.0, 1.0, 0.0, 0.0};
Real m1_rs_r[5] = {1.0, 0.0, 1.0, 0.0, 0.0};

// atmosphere: the prescribed gas and the two end states of the reference BC
Real m1_at_rho = 0.128, m1_at_h = 0.113, m1_at_ztop = 1.0, m1_at_egas = 1.0;
Real m1_at_flux = 1.0, m1_at_c = 1.0, m1_at_ebot = 1.0, m1_at_kap = 1.0;

//----------------------------------------------------------------------------------------
//! \fn M1AtmTau
//! \brief the Rosseland optical depth measured DOWNWARD from the top mesh face z = ztop
//! for rho(z) = rho_top exp[(ztop - z)/H] and a constant kappa.

KOKKOS_INLINE_FUNCTION
Real M1AtmTau(const Real z, const Real rho_top, const Real kap, const Real hh,
              const Real ztop) {
  Real s = (ztop - z)/hh;
  // expm1 keeps the top cells (s << 1) accurate: tau -> kappa rho_top (ztop - z) there
  return kap*rho_top*hh*expm1(s);
}

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
void RadM1ShockBC(Mesh *pm);
void RadM1AtmBC(Mesh *pm);
void RadM1AtmGas(Mesh *pm, const Real bdt);

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
  int &ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
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
  } else if (test.compare("radshock") == 0) {
    // the two far-field states, in the cgs scaling printed by t7_radshock.py --athinput
    Real dl = pin->GetReal("problem","m1_shock_rho_l");
    Real vl = pin->GetReal("problem","m1_shock_v_l");
    Real tl = pin->GetReal("problem","m1_shock_t_l");
    Real dr = pin->GetReal("problem","m1_shock_rho_r");
    Real vr = pin->GetReal("problem","m1_shock_v_r");
    Real trr = pin->GetReal("problem","m1_shock_t_r");
    Real xs = pin->GetReal("problem","m1_shock_xs");
    Real wid = pin->GetOrAddReal("problem","m1_shock_w",2.0);
    std::string ref = pin->GetOrAddString("problem","m1_shock_ref","none");
    // the code temperature per kelvin.  In honest cgs with the ideal-gas EOS,
    // T_code = (gamma-1) e/rho = [k/(mu m_H)] T_kelvin, so tunit = k/(mu m_H) and
    // <rad_m1>/arad must be a_cgs/tunit^4 for E = arad T_code^4 to be erg/cm^3.
    Real tun = pin->GetOrAddReal("problem","m1_shock_tunit",1.0);
    tl *= tun;
    trr *= tun;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    m1_rs_l[0] = dl;
    m1_rs_l[1] = dl*vl;
    m1_rs_l[2] = dl*tl/gm1 + 0.5*dl*vl*vl;
    m1_rs_l[3] = ar*tl*tl*tl*tl;
    m1_rs_l[4] = 0.0;
    m1_rs_r[0] = dr;
    m1_rs_r[1] = dr*vr;
    m1_rs_r[2] = dr*trr/gm1 + 0.5*dr*vr*vr;
    m1_rs_r[3] = ar*trr*trr*trr*trr;
    m1_rs_r[4] = 0.0;
    user_bcs_func = RadM1ShockBC;
    if (restart) return;

    // the semi-analytic profile, if one was given: columns x, rho, v, T_gas, T_rad
    int nref = 0;
    std::vector<Real> vx, vd, vv, vtg, vtr;
    if (ref.compare("none") != 0) {
      std::ifstream fp(ref.c_str());
      if (!fp.is_open()) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cannot open <problem>/m1_shock_ref = '" << ref << "'"
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
          << std::endl << "m1_shock_ref '" << ref << "' has fewer than 2 rows"
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
    par_for("m1_radshock_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
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
      u0(m,radm1::M1_E,k,j,i) = fmax(erad, efl);
      // the comoving flux is small in these thick shocks; the lab flux is dominated by
      // the enthalpy term, which the stiff source relaxes to within a few substeps
      u0(m,radm1::M1_F1,k,j,i) = (4.0/3.0)*v*erad;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else if (test.compare("atmosphere") == 0) {
    // milestone 3a2, gate I6b: the static grey plane-parallel column.  See the file
    // header for the closed-form steady state and the (marshak_q, q0) pair.
    Real rht = pin->GetOrAddReal("problem","atm_rho_top",0.128);
    Real hh = pin->GetOrAddReal("problem","atm_scale_h",0.113);
    Real fin = pin->GetOrAddReal("problem","atm_flux",1.0);
    Real tg = pin->GetOrAddReal("problem","atm_temp",1.0);
    Real kap = pmbp->pradm1->kappa_s + pmbp->pradm1->kappa_f;
    Real ztop = pmy_mesh_->mesh_size.x1max;
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real qq = pmbp->pradm1->marshak_q;
    Real q0 = radm1::M1Chi(qq)/qq;
    m1_at_rho = rht;
    m1_at_h = hh;
    m1_at_ztop = ztop;
    m1_at_egas = tg/gm1;
    m1_at_flux = fin;
    m1_at_c = cl;
    m1_at_kap = kap;
    // the reference (explicit-scheme) ghost state at the bottom: the closed-form E at
    // the bottom MESH FACE, which is where tau_bot sits
    Real taub = M1AtmTau(pmy_mesh_->mesh_size.x1min, rht, kap, hh, ztop);
    {
      // invert chi(f)/f = q0 + tau on the decreasing branch by bisection on f
      Real target = q0 + taub;
      Real flo = 1.0e-14, fhi = 0.8250;
      for (int n=0; n<200; ++n) {
        Real fm = 0.5*(flo + fhi);
        if (radm1::M1Chi(fm)/fm > target) {flo = fm;} else {fhi = fm;}
      }
      m1_at_ebot = fin/(cl*0.5*(flo + fhi));
    }
    user_bcs_func = RadM1AtmBC;
    user_srcs_func = RadM1AtmGas;
    if (global_variable::my_rank == 0) {
      std::cout << "  m1_test = atmosphere: kappa=" << kap << " rho_top=" << rht
                << " H=" << hh << " tau_bot=" << taub << " tau_topcell="
                << M1AtmTau(ztop - (ztop - pmy_mesh_->mesh_size.x1min)/
                            static_cast<Real>(pmy_mesh_->mesh_indcs.nx1), rht, kap, hh,
                            ztop)
                << " F=" << fin << " marshak_q=" << qq << " q0=" << q0
                << " E(tau_bot)=" << m1_at_ebot << std::endl;
    }
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    Real egas = m1_at_egas;
    // atm_seed (runs_5c_thinstab): a TRANSVERSE perturbation of the initial E,
    // E *= 1 + atm_seed cos(2 pi atm_seed_k (x2 - x2min)/L2), or (-1)^j for
    // atm_seed_k <= 0.  Default 0: the IC is untouched.
    const Real aseed = pin->GetOrAddReal("problem","atm_seed",0.0);
    const int aseedk = pin->GetOrAddInteger("problem","atm_seed_k",1);
    const Real ax2min = pmy_mesh_->mesh_size.x2min;
    const Real ax2len = pmy_mesh_->mesh_size.x2max - ax2min;
    par_for("m1_atm_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
      Real d = rht*exp((ztop - x1v)/hh);
      Real tau = M1AtmTau(x1v, rht, kap, hh, ztop);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*egas;
      // the EDDINGTON guess, deliberately not the M1 answer
      u0(m,radm1::M1_E,k,j,i) = fmax(3.0*(fin/cl)*(tau + 2.0/3.0), efl);
      if (aseed != 0.0) {
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        Real sm = (aseedk > 0) ? cos(2.0*M_PI*aseedk*(x2v - ax2min)/ax2len)
                               : ((((j - js) & 1) == 0) ? 1.0 : -1.0);
        u0(m,radm1::M1_E,k,j,i) *= 1.0 + aseed*sm;
      }
      u0(m,radm1::M1_F1,k,j,i) = fin;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
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
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Real kapf = pmbp->pradm1->kappa_f;
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
      std::cout << "  m1_test = radwave: dir=" << wdir << " beta=" << bet
                << " Gamma1=" << gam1 << " Gamma3-1=" << g3m1 << std::endl
                << "    c_s=" << cs << " c/c_s=" << (cl/cs) << " lambda=" << lam
                << " period=" << (lam/cs) << " tau_lambda=" << (kapf*d0*lam)
                << std::endl;
    }
    if (restart) return;
    auto uh = pmbp->phydro->u0;
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
      // <problem>/radwave_v0 (tests_m1/runs_3s_space2): a uniform background drift along
      // the wave vector, with the matching lab flux (4/3) v0 E.  Read only when named.
      const bool drift = pin->DoesParameterExist("problem","radwave_v0");
      const Real v0 = drift ? pin->GetReal("problem","radwave_v0") : 0.0;
      par_for("m1_radwave_eig", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        auto &sz = size.d_view(m);
        Real x1v = CellCenterX(i-is, nx1, sz.x1min, sz.x1max);
        Real x2v = (nx2 > 1) ? CellCenterX(j-js, nx2, sz.x2min, sz.x2max) : 0.0;
        Real x3v = (nx3 > 1) ? CellCenterX(k-ks, nx3, sz.x3min, sz.x3max) : 0.0;
        Real ph = kx*x1v + ky*x2v + kz*x3v;
        Real cp = cos(ph), sp = sin(ph);
        Real d = d0*(1.0 + amp*cp);
        Real tt = t0*(1.0 + amp*(tre*cp - tim*sp));
        Real vv = amp*(vre*cp - vim*sp);
        if (drift) {vv += v0;}
        uh(m,IDN,k,j,i) = d;
        uh(m,IM1,k,j,i) = d*vv*nx;
        uh(m,IM2,k,j,i) = d*vv*ny;
        uh(m,IM3,k,j,i) = d*vv*nz;
        uh(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*vv*vv;
        Real ff = amp*ce0*(fre*cp - fim*sp);
        u0(m,radm1::M1_E,k,j,i) = fmax(er0*(1.0 + amp*(ere*cp - eim*sp)), efl);
        if (drift) {
          ff += (4.0/3.0)*v0*(u0(m,radm1::M1_E,k,j,i) - er0)
                + (4.0/3.0)*v0*er0;
        }
        u0(m,radm1::M1_F1,k,j,i) = ff*nx;
        u0(m,radm1::M1_F2,k,j,i) = ff*ny;
        u0(m,radm1::M1_F3,k,j,i) = ff*nz;
      });
      // the comoving face fluxes of the implicit transport (allocated only there)
      const Real f0re = ce0*fre - (4.0/3.0)*er0*vre;
      const Real f0im = ce0*fim - (4.0/3.0)*er0*vim;
      for (int dir = 0; dir < 3; ++dir) {
        DvceArray4D<Real> ff0 = (dir == 0) ? pmbp->pradm1->f0x1 :
                                ((dir == 1) ? pmbp->pradm1->f0x2 : pmbp->pradm1->f0x3);
        if (ff0.extent_int(0) < nmb1 + 1) continue;
        const Real nd = (dir == 0) ? nx : ((dir == 1) ? ny : nz);
        const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
        const int e1 = ff0.extent_int(3) - 1;
        par_for("m1_radwave_eigf", DevExeSpace(), 0,nmb1,0,e3,0,e2,0,e1,
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
        });
      }
      return;
    }
    par_for("m1_radwave_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
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
      Real ee = er0 + derad*cp;
      Real ff = (4.0/3.0)*vv*er0 + fdif*sp;
      u0(m,radm1::M1_E,k,j,i) = fmax(ee, efl);
      u0(m,radm1::M1_F1,k,j,i) = ff*nx;
      u0(m,radm1::M1_F2,k,j,i) = ff*ny;
      u0(m,radm1::M1_F3,k,j,i) = ff*nz;
    });
  } else if (test.compare("sph_shell") == 0) {
    // STAGE S1 (tests_m1/runs_5a_sp_s1): a static uniform gas on a spherical-polar
    // wedge and a spherically symmetric radiation field, E = e_out + e_amp
    // exp(-((r - r0)/w)^2) at the cell centroid r = x1v, F = 0.  e_amp = 0 (default)
    // is the flat start of the steady-diffusion gate T-S1 (imposed flux at r_in, the
    // outer end cell held at e_out by implicit_bc_x1max = efix); e_amp > 0 is the
    // symmetric transient of T-sym.  The gas does not move (run it with
    // <rad_m1>/gas_feedback = false and a scattering opacity).
    // (On a Cartesian mesh r is x1: the planar twin used to compare solver behaviour.)
    if (pmbp->phydro == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<problem>/m1_test = sph_shell needs a <hydro> block"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Real dgas = pin->GetOrAddReal("problem","gas_rho",1.0);
    Real tgas = pin->GetOrAddReal("problem","gas_temp",1.0e-8);
    Real eout = pin->GetReal("problem","e_out");
    Real eamp = pin->GetOrAddReal("problem","e_amp",0.0);
    Real r0 = pin->GetOrAddReal("problem","r0",0.0);
    Real wid = pin->GetOrAddReal("problem","width",1.0);
    Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    if (restart) return;
    auto uh = pmbp->phydro->u0;
    auto x1v = pmbp->pcoord->x1v;
    const bool sp = pmy_mesh_->use_spherical_polar;
    par_for("m1_sph_ic", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real r = sp ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                           size.d_view(m).x1max);
      uh(m,IDN,k,j,i) = dgas;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = dgas*tgas/gm1;
      Real x = (r - r0)/wid;
      u0(m,radm1::M1_E,k,j,i) = fmax(eout + eamp*exp(-x*x), efl);
      u0(m,radm1::M1_F1,k,j,i) = 0.0;
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = '" << test << "' not implemented "
      << "(beam | pulse1d | thick_pulse | tophat | jump | equil | advect_pulse "
      << "| advect_uniform | advect_shear | marshak | shadow | radshock "
      << "| atmosphere | radwave | sph_shell)" << std::endl;
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

//----------------------------------------------------------------------------------------
//! \fn void RadM1ShockBC()
//! \brief Dirichlet at the two far-field states of the radiative shock, for BOTH the
//! radiation moments and the gas: the shock is held in its own frame.

void RadM1ShockBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pradm1 == nullptr || pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int nmb1 = (pmbp->nmb_thispack - 1);
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->pradm1->u0;
  auto uh = pmbp->phydro->u0;
  Real cl = pmbp->pradm1->c_light;
  Real dl = m1_rs_l[0], ml = m1_rs_l[1], el = m1_rs_l[2], rl = m1_rs_l[3];
  Real dr = m1_rs_r[0], mr = m1_rs_r[1], er = m1_rs_r[2], rr = m1_rs_r[3];

  par_for("m1_shock_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        uh(m,IDN,k,j,is-i-1) = dl;
        uh(m,IM1,k,j,is-i-1) = ml;
        uh(m,IM2,k,j,is-i-1) = 0.0;
        uh(m,IM3,k,j,is-i-1) = 0.0;
        uh(m,IEN,k,j,is-i-1) = el;
        u0(m,radm1::M1_E,k,j,is-i-1) = rl;
        u0(m,radm1::M1_F1,k,j,is-i-1) = (4.0/3.0)*(ml/dl)*rl;
        u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
        u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        uh(m,IDN,k,j,ie+i+1) = dr;
        uh(m,IM1,k,j,ie+i+1) = mr;
        uh(m,IM2,k,j,ie+i+1) = 0.0;
        uh(m,IM3,k,j,ie+i+1) = 0.0;
        uh(m,IEN,k,j,ie+i+1) = er;
        u0(m,radm1::M1_E,k,j,ie+i+1) = rr;
        u0(m,radm1::M1_F1,k,j,ie+i+1) = (4.0/3.0)*(mr/dr)*rr;
        u0(m,radm1::M1_F2,k,j,ie+i+1) = 0.0;
        u0(m,radm1::M1_F3,k,j,ie+i+1) = 0.0;
      }
    }
  });
  (void) cl;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1AtmGas()
//! \brief re-impose the PRESCRIBED exponential atmosphere (ghost zones included) at the
//! end of every hydro stage, so that rho*kappa is an exact function of height for the
//! whole run and the steady radiation state is the closed-form M1 one.  bdt is unused.

void RadM1AtmGas(Mesh *pm, const Real bdt) {
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
  Real rht = m1_at_rho, hh = m1_at_h, ztop = m1_at_ztop, egas = m1_at_egas;
  (void) bdt;

  par_for("m1_atm_reset", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),0,(n1-1),
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
//! \fn void RadM1AtmBC()
//! \brief the radiation ghost cells of the atmosphere test, for the EXPLICIT reference
//! run only (the implicit scheme imposes its boundaries on the FACES and never reads
//! these).  Bottom: the closed-form (E, F) at tau_bot, i.e. the imposed flux with the
//! E that carries it.  Top: a DARK ghost (M1FillGhost mode 2), which is what supplies
//! the free-surface relation through the HLL flux.

void RadM1AtmBC(Mesh *pm) {
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
  Real cl = pmbp->pradm1->c_light;
  Real efl = pmbp->pradm1->e_floor;
  Real fin = m1_at_flux, ebot = m1_at_ebot;

  par_for("m1_atm_bc", DevExeSpace(), 0,nmb1,0,(n3-1),0,(n2-1),
  KOKKOS_LAMBDA(int m, int k, int j) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        u0(m,radm1::M1_E,k,j,is-i-1) = ebot;
        u0(m,radm1::M1_F1,k,j,is-i-1) = fin;
        u0(m,radm1::M1_F2,k,j,is-i-1) = 0.0;
        u0(m,radm1::M1_F3,k,j,is-i-1) = 0.0;
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      for (int i=0; i<ng; ++i) {
        radm1::M1FillGhost(u0, m, k, j, ie+i+1, k, j, ie, 1, 2, 1.0, cl, efl);
      }
    }
  });
  return;
}
