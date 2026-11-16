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

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
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
  rh_g = pin->GetOrAddReal("problem","rh_g",1.0);
  rh_gm = pin->GetOrAddReal("problem","rh_gm",1.0);
  rh_eddf = pin->GetOrAddReal("problem","rh_eddf",0.3);
  const Real rho0 = pin->GetOrAddReal("problem","rh_rho0",1.0);
  const Real t0 = pin->GetOrAddReal("problem","rh_temp",1.0);
  rh_seed = pin->GetOrAddReal("problem","rh_seed",0.0);
  user_srcs_func = RhseGravity;
  user_hist_func = RhseHist;
  // (read only when named) problem/am_hist: the energy / angular-momentum budget instead
  if (pin->DoesParameterExist("problem","am_hist") && pin->GetBoolean("problem","am_hist")) {
    user_hist_func = RadM1AmHistory;
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
  if (restart) return;
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
//! \fn void RhseHist
//! \brief gas x1 momentum and kinetic energy, max |v|/c_s (one rank), radiation E, F1/c

void RhseHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  pdata->nhist = 5;
  pdata->label[0] = "M1_gas";
  pdata->label[1] = "KE_gas";
  pdata->label[2] = "vcs_max";
  pdata->label[3] = "E_rad";
  pdata->label[4] = "Y1_rad";
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
