//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_rhse.cpp
//! \brief <problem>/m1_test = rhse: a static grey atmosphere in hydrostatic and radiative
//! equilibrium, plane parallel (Cartesian, x1 = height, gravity rh_g along -x1) or
//! spherical (spherical polar or cubed sphere, x1 = r, gravity rh_gm/r^2).  The test of
//! the well-balanced clean scheme (<rad_m1>/force_reference = clean with
//! <hydro>/wellbalance_static).
//!
//! Gas: ideal, p = rho T.  The flux F (F0 at the bottom, F0 (r_in/r)^2 on a sphere) with
//! the constant opacity kappa = kappa_f + kappa_s gives the radiative acceleration
//! a = rh_eddf g, a fixed fraction of gravity.  T is rh_temp, or with absorption
//! (kappa_p > 0) the Eddington radiation temperature (E/arad)^(1/4); rho follows from
//! dp/dx1 = -rho g (1 - rh_eddf), integrated by the midpoint rule.  E is the Eddington
//! solution 3 (F/c) (tau + 2/3) of the isothermal column, F the flux above, the bath at
//! x1min E(x1min).  The radiation relaxes to the scheme's own steady state with
//! <rad_m1>/clean_init = steady (which with force_reference = clean also relaxes the gas
//! internal energy); after it the gas becomes the static background (all cells, ghost
//! zones included), so the discrete state of both is the reference equilibrium.
//!
//! Gravity (user_srcs): rho g, on rho - rho_background with wellbalance_static, plus its
//! work; the background is in the restart file (pgen_rst_slabs).  Mesh x1 boundaries:
//! reflect.  <problem>/rh_seed adds a cell-to-cell x1 velocity of that amplitude after
//! the background is set (a lateral seed for stability tests; in 3-D the unsplit hydro
//! step needs cfl_number <= 1/3, else its grid-scale lateral mode grows).  History: gas
//! x1 momentum, kinetic energy, max |v|/c_s (one rank), radiation E and F1/c (volume
//! sums).
//!
//! The tidally distorted envelope (spherical polar or cubed sphere, <hydro>/
//! wellbalance_dynamic, a test of <rad_m1>/force_reference = wb_arad_md with <hydro>/
//! wb_dirs = all): Phi = -rh_gm/r - rh_tide h(r) P2(X/r) (RhPhi: a quadrupole tide along
//! (cos t, sin t, 0), t = rh_tide_phi, whose x1 boundaries are equipotentials; the IC's
//! at rh_tide_phi0) and its radiation source (c rh_eddf/kappa) lap Phi (esrc),
//! rho = rh_rho0 exp(-(1 - rh_eddf) (Phi - Phi(x1min))/T), T = rh_temp (kappa_p = 0), the
//! diffusion E(Phi) = 3 rh_eddf T (rho - rho_b)/(1 - rh_eddf) + 2 F_b/c (b: the
//! unperturbed x1max; at least 2 F_b/c rho/rho_b) and F_r of the sphere; the radiation
//! relaxes with clean_init = steady (x1min: implicit_bc_x1min = flux).  Gravity: the
//! walked background (AddWbGravity, x1 in the face-sum form of the x1 cache) and its work
//! v.(-rho grad Phi); with
//! force_reference = wb_arad the radial part rh_eddf dPhi/dr of the radiation force
//! (effective x1 potential (1 - rh_eddf) Phi), with wb_arad_md the module's own, and
//! their work (force_reference_work = split).  problem/wb_lobe_cut (default true with a
//! tide): no WB and no reference above Phi(L1) = -3 rh_gm/(2 r_L1), r_L1 = (rh_gm/(2
//! rh_tide))^(1/3).  History: + max |v_r|/c_s and max |v_lat|/c_s, all cells and r <
//! (x1min + x1max)/2 (one rank).

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "outputs/outputs.hpp"
#include "rad_m1/rad_m1.hpp"
#include "pgen/pgen.hpp"

void RadM1AmHistory(HistoryData *pdata, Mesh *pm);   // rad_m1_tests2.cpp

namespace {
Real rh_g = 1.0, rh_gm = 0.0, rh_eddf = 0.0, rh_seed = 0.0;
bool rh_sph = false;
bool rh_wb = false;     // the tidally distorted envelope (<hydro>/wellbalance_dynamic)
Real rh_tide = 0.0, rh_tph = 0.0, rh_tph0 = 0.0, rh_t0 = 1.0, rh_a = 1.0, rh_b = 2.0;
bool rh_am = false;     // problem/am_hist: + the torque budget (RhseHist)
Real rh_tq[3] = {0.0, 0.0, 0.0};   // accumulated torque_z: WB gravity, withheld, analytic
Real rh_dtp = 0.0, rh_sig[2] = {0.0, 0.0};   // (the step's rk2 combination so far)
int rh_cyc = -1;
Real RhseLz(Mesh *pm, const int w);

//! the position (Cartesian) of (r, x2, x3): sp theta, phi; cs xi, eta over pi/4
KOKKOS_INLINE_FUNCTION
void RhPos(const bool cs, const int pan, const Real r, const Real a2, const Real a3,
           Real x[3]) {
  if (cs) {
    cubed_sphere::PanelToCart(pan, 0.25*M_PI*a2, 0.25*M_PI*a3, x);
  } else {
    x[0] = sin(a2)*cos(a3); x[1] = sin(a2)*sin(a3); x[2] = cos(a2);
  }
  for (int n = 0; n < 3; ++n) x[n] *= r;
}

//! Phi and -grad Phi (Cartesian) of the point mass and the tide along (cos t, sin t, 0),
//! -q h(r) (3 X^2/r^2 - 1)/2, h = 16 (r - ra)^2 (rb - r)^2/(rb - ra)^4: the x1min and
//! x1max spheres are equipotentials with a uniform dPhi/dr (the imposed inner flux and
//! the outer vacuum of the diffusion equilibrium); lap the Laplacian of the tide
KOKKOS_INLINE_FUNCTION
Real RhPhi(const Real gm, const Real q, const Real t, const Real ra, const Real rb,
           const Real x[3], Real g[3], Real *lap = nullptr) {
  const Real r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2], r = sqrt(r2);
  const Real ct = cos(t), st = sin(t), xt = ct*x[0] + st*x[1], p2 = 1.5*xt*xt/r2 - 0.5;
  const Real k = 16.0/SQR(SQR(rb - ra)), u = r - ra, w = rb - r;
  const Real f = -q*k*u*u*w*w, f1 = -2.0*q*k*u*w*(w - u);
  const Real f2 = -2.0*q*k*(w*w - 4.0*u*w + u*u);
  // -grad (f P2): radial -f' P2, angular -f grad P2 = -f 3 xt (t - xt x/r^2)/r^2
  const Real gp = 3.0*f*xt/r2;
  for (int n = 0; n < 3; ++n) g[n] = -gm*x[n]/(r2*r) - f1*p2*x[n]/r + gp*xt*x[n]/r2;
  g[0] -= gp*ct;
  g[1] -= gp*st;
  if (lap != nullptr) *lap = (f2 + 2.0*f1/r - 6.0*f/r2)*p2;
  return -gm/r + f*p2;
}

//! the unit vectors of the hydro's components at the cell centre: rhat and (sp) e_theta,
//! e_phi or (cs) the panel tangents e_xi, e_eta (covariant components a.e)
KOKKOS_INLINE_FUNCTION
void RhBasis(const bool cs, const int pan, const Real a2, const Real a3, Real e[3][3]) {
  RhPos(cs, pan, 1.0, a2, a3, e[0]);
  if (cs) {
    cubed_sphere::PanelTangents(pan, 0.25*M_PI*a2, 0.25*M_PI*a3, e[1], e[2]);
  } else {
    e[1][0] = cos(a2)*cos(a3); e[1][1] = cos(a2)*sin(a3); e[1][2] = -sin(a2);
    e[2][0] = -sin(a3); e[2][1] = cos(a3); e[2][2] = 0.0;
  }
}

//! gravity at x1 (height or radius)
KOKKOS_INLINE_FUNCTION
Real RhGravity(const bool sph, const Real g, const Real gm, const Real x) {
  return sph ? gm/(x*x) : g;
}

//! the isothermal column at x1 = x on [a, b]: rho, tau from the top face and the flux;
//! H = t0/(g (1 - ef)) on a plane, the tau integral by the midpoint rule on a sphere
KOKKOS_INLINE_FUNCTION
void RhseColumn(const bool sph, const Real x, const Real a, const Real b, const Real rho0,
                const Real t0, const Real g, const Real gm, const Real ef, const Real f0,
                const Real kap, Real &d, Real &tau, Real &f) {
  if (sph) {
    const Real q = gm*(1.0 - ef)/t0;
    d = rho0*exp(q*(1.0/x - 1.0/a));
    f = f0*a*a/(x*x);
    tau = 0.0;
    const Real dx = (b - x)/64.0;
    for (int n = 0; n < 64; ++n) {
      tau += kap*rho0*exp(q*(1.0/(x + (n + 0.5)*dx) - 1.0/a))*dx;
    }
  } else {
    const Real hh = t0/(g*(1.0 - ef));
    d = rho0*exp(-(x - a)/hh);
    f = f0;
    tau = kap*rho0*hh*(exp(-(x - a)/hh) - exp(-(b - a)/hh));
  }
}

//! the Eddington E = 3 (F/c) (tau + 2/3) of that column
KOKKOS_INLINE_FUNCTION
Real RhseEdd(const bool sph, const Real x, const Real a, const Real b, const Real rho0,
             const Real t0, const Real g, const Real gm, const Real ef, const Real f0,
             const Real kap, const Real cl) {
  Real d, tau, f;
  RhseColumn(sph, x, a, b, rho0, t0, g, gm, ef, f0, kap, d, tau, f);
  return 3.0*(f/cl)*(tau + 2.0/3.0);
}

void RhseGravity(Mesh *pm, const Real bdt);
void RhseWbGravity(Mesh *pm, const Real bdt);
void RhsePotential(Mesh *pm, const Real ef);
void RhseRelax(Mesh *pm);
void RhseBackground(Mesh *pm);
void RhseBgDerive(Mesh *pm);
void RhseHist(HistoryData *pdata, Mesh *pm);
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationM1Rhse
//! \brief m1_test = rhse, see the file header

void ProblemGenerator::RadiationM1Rhse(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  radm1::RadiationM1 *pm1 = pmbp->pradm1;
  hydro::Hydro *ph = pmbp->phydro;
  if (ph == nullptr || !ph->peos->eos_data.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<problem>/m1_test = rhse needs <hydro> with an ideal gas"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rh_sph = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
  rh_a = pmy_mesh_->mesh_size.x1min;
  rh_b = pmy_mesh_->mesh_size.x1max;
  rh_g = pin->GetOrAddReal("problem","rh_g",1.0);
  rh_gm = pin->GetOrAddReal("problem","rh_gm",1.0);
  rh_eddf = pin->GetOrAddReal("problem","rh_eddf",0.3);
  const Real rho0 = pin->GetOrAddReal("problem","rh_rho0",1.0);
  const Real t0 = pin->GetOrAddReal("problem","rh_temp",1.0);
  rh_t0 = t0;
  rh_seed = pin->GetOrAddReal("problem","rh_seed",0.0);
  rh_wb = ph->use_wellbalance_dynamic;
  if (rh_wb && (!rh_sph || pm1->kappa_p > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "m1_test = rhse with <hydro>/wellbalance_dynamic: "
              << "spherical polar or cubed sphere, kappa_p = 0" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (rh_wb) {
    rh_tide = pin->GetOrAddReal("problem","rh_tide",0.0);
    rh_tph = pin->GetOrAddReal("problem","rh_tide_phi",0.0);
    rh_tph0 = pin->GetOrAddReal("problem","rh_tide_phi0",rh_tph);
  }
  user_srcs_func = rh_wb ? RhseWbGravity : RhseGravity;
  user_hist_func = RhseHist;
  // (read only when named) problem/am_hist: the energy / angular-momentum budget instead
  if (pin->DoesParameterExist("problem","am_hist") &&
      pin->GetBoolean("problem","am_hist")) {
    user_hist_func = RadM1AmHistory;
    rh_am = rh_wb;
    if (rh_am) user_hist_func = RhseHist;
  }
  if (ph->use_wellbalance_static) {
    // the background travels in the restart file
    user_icpost_func = RhseBackground;
    pgen_rst_slabs = &ph->u0wb;
    if (restart && pgen_rst_slabs_in.extent_int(1) == 5) {
      Kokkos::deep_copy(ph->u0wb, pgen_rst_slabs_in);
      RhseBgDerive(pmy_mesh_);
    }
  }

  const Real cl = pm1->c_light, ar = pm1->arad, kap = pm1->kappa_f + pm1->kappa_s;
  const Real gm1 = ph->peos->eos_data.gamma - 1.0;
  const bool teq = (pm1->kappa_p > 0.0);
  const Real a = pmy_mesh_->mesh_size.x1min, b = pmy_mesh_->mesh_size.x1max;
  const bool sph = rh_sph;
  const Real g = rh_g, gm = rh_gm, ef = rh_eddf;
  // flux at the bottom: kappa F0/c = rh_eddf g(x1min)
  const Real f0 = ef*RhGravity(sph, g, gm, a)*cl/kap;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng, nmb1 = pmbp->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto x1v = pmbp->pcoord->x1v;
  auto &size = pmbp->pmb->mb_size;
  const int is = indcs.is, nx1 = indcs.nx1;
  auto uh = ph->u0;
  auto u0 = pm1->u0;
  const Real efl = pm1->e_floor;
  pm1->iebath_x1min = RhseEdd(sph, a, a, b, rho0, t0, g, gm, ef, f0, kap, cl);
  // (the tidal envelope: E(Phi) of the diffusion limit, 2 F_b/c at rho_b of the
  // unperturbed x1max, at least 2 F_b/c rho/rho_b; the bath of the unperturbed x1min)
  const Real qt = rh_tide;
  const Real ce = 3.0*ef*t0/(1.0 - ef), e0 = 2.0*ef*gm/(b*b*kap);
  const Real rb = rho0*exp(-(1.0 - ef)*gm*(1.0/a - 1.0/b)/t0);
  if (rh_wb) {
    pm1->iebath_x1min = ce*(rho0 - rb) + e0;
    // problem/wb_lobe_cut: no WB above Phi(L1), on the x1 WB potential (wb_arad: its
    // (1 - rh_eddf) Phi)
    const bool rfa = (pm1->force_ref == radm1::M1_FREF_WB_ARAD);
    if (qt > 0.0 && pin->GetOrAddBoolean("problem","wb_lobe_cut",true)) {
      // (L1: the first radius on the tide's axis where g_r turns outward; none inside
      // the domain: no cut)
      auto gax = [&](const Real r) {
        Real x[3] = {r*cos(rh_tph), r*sin(rh_tph), 0.0}, gg[3];
        (void) RhPhi(gm, qt, rh_tph, a, b, x, gg);
        return gg[0]*cos(rh_tph) + gg[1]*sin(rh_tph);
      };
      Real rl = -1.0;
      for (int n = 1; n <= 1000 && rl < 0.0; ++n) {
        Real lo = a + (n - 1)*(b - a)/1000.0, hi = a + n*(b - a)/1000.0;
        if (!(gax(hi) >= 0.0)) continue;
        for (int it = 0; it < 60; ++it) {
          const Real mid = 0.5*(lo + hi);
          if (gax(mid) >= 0.0) {hi = mid;} else {lo = mid;}
        }
        rl = 0.5*(lo + hi);
      }
      if (rl > 0.0) {
        Real x[3] = {rl*cos(rh_tph), rl*sin(rh_tph), 0.0}, gg[3];
        ph->wb_phimax = RhPhi(gm, qt, rh_tph, a, b, x, gg)*(rfa ? 1.0 - ef : 1.0);
      }
      if (global_variable::my_rank == 0) {
        std::printf("rhse: problem/wb_lobe_cut: r_L1 %.6f, Phi_wb(L1) %.6e\n", rl,
                    ph->wb_phimax);
      }
    }
    if (rfa) ph->EnableWBEffectivePotential();
    // the tide's radiation source (c rh_eddf/kappa) lap Phi = div F of the equilibrium
    // (esrc, every start)
    if (qt > 0.0) {
      Kokkos::realloc(pm1->esrc, nmb1+1, n3, n2, n1);
      pm1->esrc_on = true;
      auto es = pm1->esrc;
      auto x2v = pmbp->pcoord->x2v, x3v = pmbp->pcoord->x3v;
      auto &pan = pmbp->pmb->mb_panel;
      const bool cs = pmy_mesh_->use_cubed_sphere;
      const Real tp = rh_tph, sc = cl*ef/kap;
      par_for("rhse_tide_esrc", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real x[3], gg[3], lp;
        RhPos(cs, cs ? pan.d_view(m) : 0, x1v(m,i), x2v(m,j), x3v(m,k), x);
        (void) RhPhi(gm, qt, tp, a, b, x, gg, &lp);
        es(m,k,j,i) = sc*lp;
      });
    }
    // problem/rh_relax (read only when named; wb_arad_md, one MeshBlock in x1): after
    // the steady start the gas of each column in the x1 walk's equilibrium with the
    // reference (the steady radiation's force)
    if (pm1->force_ref == radm1::M1_FREF_WB_MD &&
        pin->DoesParameterExist("problem","rh_relax") &&
        pin->GetBoolean("problem","rh_relax")) {
      user_icpost_func = RhseRelax;
    }
    RhsePotential(pmy_mesh_, rfa ? ef : 0.0);
    if (pm1->force_ref != radm1::M1_FREF_NONE) {
      pm1->fref_wsplit_ok = true;   // (RhseWbGravity gives the reference work)
      pm1->SetForceReference(pm1->arad_ref);
    }
  }
  if (restart) return;
  if (rh_wb) {
    auto x2v = pmbp->pcoord->x2v, x3v = pmbp->pcoord->x3v;
    auto &pan = pmbp->pmb->mb_panel;
    const bool cs = pmy_mesh_->use_cubed_sphere;
    const Real t0p = rh_tph0;
    par_for("rhse_tide_ic", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real x[3], gg[3];
      RhPos(cs, cs ? pan.d_view(m) : 0, x1v(m,i), x2v(m,j), x3v(m,k), x);
      const Real d = rho0*exp(-(1.0 - ef)*(RhPhi(gm, qt, t0p, a, b, x, gg) + gm/a)/t0);
      const Real r = x1v(m,i);
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      uh(m,IEN,k,j,i) = d*t0/gm1;
      u0(m,radm1::M1_E,k,j,i) = fmax(fmax(ce*(d - rb) + e0, e0*d/rb), efl);
      u0(m,radm1::M1_F1,k,j,i) = f0*a*a/(r*r);
      u0(m,radm1::M1_F2,k,j,i) = 0.0;
      u0(m,radm1::M1_F3,k,j,i) = 0.0;
    });
    return;
  }
  par_for("rhse_ic", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x = sph ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                                size.d_view(m).x1max);
    Real d, tau, f;
    RhseColumn(sph, x, a, b, rho0, t0, g, gm, ef, f0, kap, d, tau, f);
    const Real e = 3.0*(f/cl)*(tau + 2.0/3.0);
    // p from the bottom by the midpoint rule on 64 steps, T(x) of the Eddington E
    Real p = rho0*t0, tt = t0;
    if (teq) {
      p = rho0*pow(RhseEdd(sph, a, a, b, rho0, t0, g, gm, ef, f0, kap, cl)/ar, 0.25);
      const Real dx = (x - a)/64.0;
      for (int n = 0; n < 64; ++n) {
        const Real xm = a + (n + 0.5)*dx;
        tt = pow(RhseEdd(sph, xm, a, b, rho0, t0, g, gm, ef, f0, kap, cl)/ar, 0.25);
        p *= exp(-RhGravity(sph, g, gm, xm)*(1.0 - ef)*dx/tt);
      }
      tt = pow(e/ar, 0.25);
      d = p/tt;
    }
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = d*tt/gm1;
    u0(m,radm1::M1_E,k,j,i) = fmax(e, efl);
    u0(m,radm1::M1_F1,k,j,i) = f;
    u0(m,radm1::M1_F2,k,j,i) = 0.0;
    u0(m,radm1::M1_F3,k,j,i) = 0.0;
  });
  if (ph->use_wellbalance_static) RhseBackground(pmy_mesh_);
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void RhseBackground
//! \brief the current gas, everywhere, as the static well-balanced background; the face
//! background is the mean of the two cells.  Then a cell-to-cell x1 velocity seed of
//! amplitude rh_seed on top (a lateral seed for stability tests; default 0)

void RhseBackground(Mesh *pm) {
  hydro::Hydro *ph = pm->pmb_pack->phydro;
  auto uh = ph->u0;
  Kokkos::deep_copy(ph->u0wb, Kokkos::subview(uh, Kokkos::ALL, std::make_pair(0, 5),
                                              Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  RhseBgDerive(pm);
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  const int n3 = uh.extent_int(2), n2 = uh.extent_int(3), n1 = uh.extent_int(4);
  const Real sd = rh_seed;
  if (sd == 0.0) return;
  par_for("rhse_seed", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dv = sd*sin(6.2831853*(0.618*i + 1.414*j + 2.236*k + 0.577*m));
    const Real m1 = uh(m,IM1,k,j,i);
    uh(m,IM1,k,j,i) += uh(m,IDN,k,j,i)*dv;
    uh(m,IEN,k,j,i) += 0.5*(SQR(uh(m,IM1,k,j,i)) - m1*m1)/uh(m,IDN,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RhseBgDerive
//! \brief the primitive and face background of the conserved one (u0wb)

void RhseBgDerive(Mesh *pm) {
  hydro::Hydro *ph = pm->pmb_pack->phydro;
  auto u0wb = ph->u0wb, w0wb = ph->w0wb;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  const int n3 = u0wb.extent_int(2), n2 = u0wb.extent_int(3), n1 = u0wb.extent_int(4);
  par_for("rhse_bg", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = u0wb(m,IDN,k,j,i);
    w0wb(m,IDN,k,j,i) = d;
    for (int n = IVX; n <= IVZ; ++n) w0wb(m,n,k,j,i) = u0wb(m,n,k,j,i)/d;
    w0wb(m,IEN,k,j,i) = u0wb(m,IEN,k,j,i) - 0.5*(SQR(u0wb(m,IM1,k,j,i)) +
                        SQR(u0wb(m,IM2,k,j,i)) + SQR(u0wb(m,IM3,k,j,i)))/d;
  });
  for (int dir = 0; dir < 3; ++dir) {
    if ((dir == 1 && n2 == 1) || (dir == 2 && n3 == 1)) continue;
    auto wf = (dir == 0) ? ph->w0facewb.x1f : ((dir == 1) ? ph->w0facewb.x2f
                                                           : ph->w0facewb.x3f);
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
    par_for("rhse_bgf", DevExeSpace(), 0, nmb1, 0, n3-1+dk, 0, n2-1+dj, 0, n1-1+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = (k > 0) ? k - dk : k, jl = (j > 0) ? j - dj : j;
      const int il = (i > 0) ? i - di : i;
      const int kr = (k < n3) ? k : k - dk, jr = (j < n2) ? j : j - dj;
      const int ir = (i < n1) ? i : i - di;
      for (int n = 0; n < 5; ++n) {
        wf(m,n,k,j,i) = 0.5*(w0wb(m,n,kl,jl,il) + w0wb(m,n,kr,jr,ir));
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RhseGravity
//! \brief gravity along -x1 and its work; on rho - rho_background with wellbalance_static

void RhseGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nx1 = indcs.nx1;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto w0wb = ph->w0wb;
  auto x1v = pmbp->pcoord->x1v;
  auto &size = pmbp->pmb->mb_size;
  const bool wbs = ph->use_wellbalance_static, sph = rh_sph;
  const Real g = rh_g, gm = rh_gm;
  par_for("rhse_grav", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x = sph ? x1v(m,i) : CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                                size.d_view(m).x1max);
    const Real d = w0(m,IDN,k,j,i) - (wbs ? w0wb(m,IDN,k,j,i) : 0.0);
    const Real s = -bdt*RhGravity(sph, g, gm, x)*d;
    u0(m,IM1,k,j,i) += s;
    u0(m,IEN,k,j,i) += s*w0(m,IVX,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RhsePotential
//! \brief the tidal envelope: Phi in every cell and on every face (ghosts included,
//! every start); force_reference = wb_arad (ef > 0): the x1 WB potential (1 - ef) Phi
//! and the reference ef dPhi/dr of the x1 faces (none outside the lobe)

void RhsePotential(Mesh *pm, const Real ef) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng, n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto x1v = pmbp->pcoord->x1v, x2v = pmbp->pcoord->x2v, x3v = pmbp->pcoord->x3v;
  auto f1 = pmbp->pcoord->xx1f, f2 = pmbp->pcoord->xx2f, f3 = pmbp->pcoord->xx3f;
  auto pc = ph->phicc0, p1 = ph->phi0.x1f, p2 = ph->phi0.x2f, p3 = ph->phi0.x3f;
  auto wc = ph->phicc_wb, w1 = ph->phi_wb_x1f;
  auto aref = pmbp->pradm1->arad_ref;
  const bool cs = pm->use_cubed_sphere, fa = (ef > 0.0);
  auto &pan = pmbp->pmb->mb_panel;
  const Real gm = rh_gm, q = rh_tide, t = rh_tph, ra = rh_a, rb = rh_b;
  const Real pmx = ph->wb_phimax;
  par_for("rhse_phi", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int pn = cs ? pan.d_view(m) : 0;
    Real x[3], g[3];
    auto phi = [&](const Real r, const Real b2, const Real b3) {
      RhPos(cs, pn, r, b2, b3, x);
      return RhPhi(gm, q, t, ra, rb, x, g);
    };
    pc(m,k,j,i) = phi(x1v(m,i), x2v(m,j), x3v(m,k));
    p1(m,k,j,i) = phi(f1(m,i), x2v(m,j), x3v(m,k));
    const Real pr = phi(f1(m,i+1), x2v(m,j), x3v(m,k));
    if (i == n1 - 1) p1(m,k,j,i+1) = pr;
    p2(m,k,j,i) = phi(x1v(m,i), f2(m,j), x3v(m,k));
    if (j == n2 - 1) p2(m,k,j+1,i) = phi(x1v(m,i), f2(m,j+1), x3v(m,k));
    p3(m,k,j,i) = phi(x1v(m,i), x2v(m,j), f3(m,k));
    if (k == n3 - 1) p3(m,k+1,j,i) = phi(x1v(m,i), x2v(m,j), f3(m,k+1));
    if (!fa) return;
    const bool in = !(pc(m,k,j,i) > pmx/(1.0 - ef));
    wc(m,k,j,i) = (1.0 - ef)*pc(m,k,j,i);
    w1(m,k,j,i) = (1.0 - ef)*p1(m,k,j,i);
    if (i == n1 - 1) w1(m,k,j,i+1) = (1.0 - ef)*pr;
    aref(m,k,j,i) = in ? ef*(pr - p1(m,k,j,i))/(f1(m,i+1) - f1(m,i)) : 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real RhseLz
//! \brief the torque budget of the tidal envelope (spherical polar, this rank): w = 0
//! sum V r sin(theta) m_phi of the gas, 1 the same of rho a_phi of wb_arad_md, 2 of the
//! analytic tidal force rho (-grad Phi)_phi

Real RhseLz(Mesh *pm, const int w) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  radm1::RadiationM1 *pm1 = pmbp->pradm1;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto u0 = pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v, x2v = pmbp->pcoord->x2v, x3v = pmbp->pcoord->x3v;
  auto fr = pm1->cln_fref;
  const bool fm = (pm1->force_ref == radm1::M1_FREF_WB_MD) && pm1->cln_fref_set;
  if (w == 1 && !fm) return 0.0;
  const Real gm = rh_gm, q = rh_tide, t = rh_tph, ra = rh_a, rb = rh_b;
  Real a = 0.0;
  Kokkos::parallel_reduce("rhse_lz", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int n, Real &sl) {
    const int m = n/nkji, k = ks + (n%nkji)/nji, j = js + (n%nji)/nx1, i = is + n%nx1;
    const Real rs = vol(m,k,j,i)*x1v(m,i)*sin(x2v(m,j));
    if (w == 0) {
      sl += rs*u0(m,IM3,k,j,i);
    } else if (w == 1) {
      sl += rs*u0(m,IDN,k,j,i)*fr(m,2,k,j,i);
    } else {
      Real x[3], g[3], e[3][3];
      RhPos(false, 0, x1v(m,i), x2v(m,j), x3v(m,k), x);
      (void) RhPhi(gm, q, t, ra, rb, x, g);
      RhBasis(false, 0, x2v(m,j), x3v(m,k), e);
      sl += rs*u0(m,IDN,k,j,i)*(g[0]*e[2][0] + g[1]*e[2][1] + g[2]*e[2][2]);
    }
  }, Kokkos::Sum<Real>(a));
  return a;
}

//----------------------------------------------------------------------------------------
//! \fn void RhseRelax
//! \brief problem/rh_relax: rho (isothermal, T = rh_temp) of each column from its x1min
//! cell along the x1 walk of wb_arad_md, rho_{i+1} = rho_i exp(-(Phi_{i+1} - Phi_i - P_1
//! (i) - P_0(i+1))/T), the x1 background of the steady radiation's force

void RhseRelax(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 != pm->mesh_indcs.nx1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "problem/rh_relax needs one MeshBlock in x1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto uh = ph->u0;
  auto pc = ph->phicc0;
  auto fr = pmbp->pradm1->cln_fref;
  const Real t0 = rh_t0;
  const Real gm1 = ph->peos->eos_data.gamma - 1.0;
  par_for("rhse_relax", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int i = is; i < ie; ++i) {
      const Real d = uh(m,IDN,k,j,i)*exp(-(pc(m,k,j,i+1) - pc(m,k,j,i) -
                     fr(m,4,k,j,i) - fr(m,3,k,j,i+1))/t0);
      uh(m,IDN,k,j,i+1) = d;
      uh(m,IEN,k,j,i+1) = d*t0/gm1;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RhseWbGravity
//! \brief the tidal envelope's gravity: the walked background (x2, x3 AddWbGravity; x1
//! the face-sum form of the x1 cache, the plain -rho dPhi/dr above wb_phimax), its work
//! v.(-rho grad Phi) and with force_reference_work = split the reference work rho a.v of
//! wb_arad (radial) or wb_arad_md (the module's acceleration, cs: on the hydro's basis)

void RhseWbGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  radm1::RadiationM1 *pm1 = pmbp->pradm1;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto u0 = ph->u0, w0 = ph->w0, c1 = ph->wbq0;
  auto pc = ph->phicc_wb, p1 = ph->phi0.x1f;
  auto vol = pmbp->pcoord->volume, a1 = pmbp->pcoord->area.x1f;
  auto dx1 = pmbp->pcoord->dx1;
  auto x1v = pmbp->pcoord->x1v, x2v = pmbp->pcoord->x2v, x3v = pmbp->pcoord->x3v;
  const bool cs = pm->use_cubed_sphere, all = ph->wb_all;
  auto &pan = pmbp->pmb->mb_panel;
  const Real gm = rh_gm, q = rh_tide, t = rh_tph, ra = rh_a, rb = rh_b;
  const Real pmx = ph->wb_phimax;
  const int fref = pm1->force_ref;
  const bool fws = pm1->fref_wsplit && (fref != radm1::M1_FREF_NONE);
  auto aref = pm1->arad_ref;
  auto fr = pm1->cln_fref;
  const bool mdo = (fref == radm1::M1_FREF_WB_MD) && pm1->cln_fref_set;
  const bool acc = fws && pm1->FrefWaccOn();
  auto wacc = acc ? pm1->FrefWacc() : DvceArray4D<Real>();
  const Real g0 = acc ? pm1->FrefWaccGam0(bdt, pm->dt) : 0.0;
  // (problem/am_hist, spherical polar: the torque budget, RhseLz; first the force the
  // radiation withheld in the step before)
  const bool tq = rh_am && !cs;
  // (integrator = rk2: a stage's increment enters the step as sig = 0.5 sig + inc, the
  // gas's combination)
  const bool ns = tq && pm->ncycle != rh_cyc;
  if (ns) {
    rh_tq[1] -= rh_dtp*RhseLz(pm, 1);
    for (int n = 0; n < 2; ++n) {rh_tq[2*n] += rh_sig[n]; rh_sig[n] = 0.0;}
    rh_cyc = pm->ncycle;
    rh_dtp = pm->dt;
  }
  const Real l0 = tq ? RhseLz(pm, 0) : 0.0;
  if (tq) rh_sig[1] = (ns ? 0.0 : 0.5*rh_sig[1]) + bdt*RhseLz(pm, 2);
  if (all) ph->AddWbGravity(bdt, false);
  par_for("rhse_wbgrav", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int pn = cs ? pan.d_view(m) : 0;
    const Real d = w0(m,IDN,k,j,i);
    Real x[3], g[3], e[3][3], gh[3];
    RhPos(cs, pn, x1v(m,i), x2v(m,j), x3v(m,k), x);
    (void) RhPhi(gm, q, t, ra, rb, x, g);
    RhBasis(cs, pn, x2v(m,j), x3v(m,k), e);
    for (int n = 0; n < 3; ++n) gh[n] = g[0]*e[n][0] + g[1]*e[n][1] + g[2]*e[n][2];
    const Real pcc = c1(m,12,k,j,i);
    Real s1 = bdt*(a1(m,k,j,i+1)*(c1(m,13,k,j,i) - pcc) +
                   a1(m,k,j,i)*(pcc - c1(m,11,k,j,i)))/vol(m,k,j,i);
    if (pc(m,k,j,i) > pmx) s1 = -bdt*d*(p1(m,k,j,i+1) - p1(m,k,j,i))/dx1(m,k,j,i);
    u0(m,IM1,k,j,i) += s1;
    if (!all) {   // (spherical polar only: the plain lateral force)
      u0(m,IM2,k,j,i) += bdt*d*gh[1];
      u0(m,IM3,k,j,i) += bdt*d*gh[2];
    }
    Real w = 0.0;
    for (int n = 0; n < 3; ++n) w += w0(m,IVX+n,k,j,i)*gh[n];
    u0(m,IEN,k,j,i) += bdt*d*w;
    if (!fws) return;
    Real ar[3] = {(fref == radm1::M1_FREF_WB_ARAD) ? aref(m,k,j,i) : 0.0, 0.0, 0.0};
    for (int n = 0; n < 3 && mdo; ++n) {
      ar[n] = cs ? fr(m,0,k,j,i)*e[n][0] + fr(m,1,k,j,i)*e[n][1] + fr(m,2,k,j,i)*e[n][2]
                 : fr(m,n,k,j,i);
    }
    Real wr = 0.0;
    for (int n = 0; n < 3; ++n) wr += bdt*d*ar[n]*w0(m,IVX+n,k,j,i);
    u0(m,IEN,k,j,i) += wr;
    if (acc) wacc(m,k,j,i) = g0*wacc(m,k,j,i) + wr;
  });
  if (tq) rh_sig[0] = (ns ? 0.0 : 0.5*rh_sig[0]) + RhseLz(pm, 0) - l0;
}

//----------------------------------------------------------------------------------------
//! \fn void RhseHist
//! \brief gas x1 momentum and kinetic energy, max |v|/c_s (one rank), radiation E, F1/c

void RhseHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (rh_am) {   // (the budget of RadM1AmHistory + the accumulated torques, RhseLz)
    RadM1AmHistory(pdata, pm);
    const int n0 = pdata->nhist;
    const char *lab[4] = {"Tgrav", "Twith", "Tanal", "Twpend"};
    for (int n = 0; n < 4; ++n) pdata->label[n0+n] = lab[n];
    for (int n = 0; n < 3; ++n) {
      pdata->hdata[n0+n] = rh_tq[n] + ((n == 1) ? 0.0 : rh_sig[n/2]);
    }
    pdata->hdata[n0+3] = (rh_cyc >= 0) ? -rh_dtp*RhseLz(pm, 1) : 0.0;
    pdata->nhist = n0 + 4;
    return;
  }
  pdata->nhist = rh_wb ? 9 : 5;
  pdata->label[0] = "M1_gas";
  pdata->label[1] = "KE_gas";
  pdata->label[2] = "vcs_max";
  pdata->label[3] = "E_rad";
  pdata->label[4] = "Y1_rad";
  pdata->label[5] = "vr_max";
  pdata->label[6] = "vlat_max";
  pdata->label[7] = "vr_deep";
  pdata->label[8] = "vlat_deep";
  const Real rdp = 0.5*(pm->mesh_size.x1min + pm->mesh_size.x1max);
  auto x1v = pmbp->pcoord->x1v;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w0 = pmbp->phydro->w0;
  auto ur = pmbp->pradm1->u0;
  auto vol = pmbp->pcoord->volume;
  auto &size = pmbp->pmb->mb_size;
  const bool sph = rh_sph;
  const Real gam = pmbp->phydro->peos->eos_data.gamma, ic = 1.0/pmbp->pradm1->c_light;
  Real s[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (rh_wb) {   // (the tidal envelope: radial and lateral |v|/c_s, all and r < rdp)
    Real v[4] = {0.0, 0.0, 0.0, 0.0};
    Kokkos::parallel_reduce("rhse_histv", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            pmbp->nmb_thispack*nkji),
    KOKKOS_LAMBDA(const int q, Real &v0, Real &v1, Real &v2, Real &v3) {
      const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/nx1, i = is + q%nx1;
      const Real cs = std::sqrt(gam*(gam - 1.0)*w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i));
      const Real vr = fabs(w0(m,IVX,k,j,i))/cs;
      const Real vl = std::sqrt(SQR(w0(m,IVY,k,j,i)) + SQR(w0(m,IVZ,k,j,i)))/cs;
      v0 = fmax(v0, vr); v1 = fmax(v1, vl);
      if (x1v(m,i) < rdp) {v2 = fmax(v2, vr); v3 = fmax(v3, vl);}
    }, Kokkos::Max<Real>(v[0]), Kokkos::Max<Real>(v[1]), Kokkos::Max<Real>(v[2]),
    Kokkos::Max<Real>(v[3]));
    for (int n = 0; n < 4; ++n) pdata->hdata[5+n] = v[n];
  }
  Kokkos::parallel_reduce("rhse_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int q, Real &s0, Real &s1, Real &s2, Real &s3, Real &s4) {
    const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/nx1, i = is + q%nx1;
    const Real dv = sph ? vol(m,k,j,i) : size.d_view(m).dx1*size.d_view(m).dx2*
                                          size.d_view(m).dx3;
    const Real d = w0(m,IDN,k,j,i);
    const Real v2 = SQR(w0(m,IVX,k,j,i)) + SQR(w0(m,IVY,k,j,i)) + SQR(w0(m,IVZ,k,j,i));
    s0 += dv*d*w0(m,IVX,k,j,i);
    s1 += dv*0.5*d*v2;
    s2 = fmax(s2, std::sqrt(v2*d/(gam*(gam - 1.0)*w0(m,IEN,k,j,i))));
    s3 += dv*ur(m,radm1::M1_E,k,j,i);
    s4 += dv*ic*ur(m,radm1::M1_F1,k,j,i);
  }, Kokkos::Sum<Real>(s[0]), Kokkos::Sum<Real>(s[1]), Kokkos::Max<Real>(s[2]),
  Kokkos::Sum<Real>(s[3]), Kokkos::Sum<Real>(s[4]));
  for (int n = 0; n < 5; ++n) pdata->hdata[n] = s[n];
}
}  // namespace
