//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_rhse.cpp
//! \brief <problem>/vet_test = rhse: a static grey atmosphere in hydrostatic and
//! radiative equilibrium, plane parallel (Cartesian, x1 = height, gravity rh_g along -x1)
//! or
//! spherical (spherical polar or cubed sphere, x1 = r, gravity rh_gm/r^2).  The test of
//! the well-balanced clean scheme (<rad_vet>/force_reference = clean with
//! <hydro|mhd>/wellbalance_static, or = wb_arad_md with wellbalance_dynamic and wb_dirs =
//! all).  <hydro> or <mhd>, eos = ideal or general.
//!
//! Gas: ideal, p = rho T.  The flux F (F0 at the bottom, F0 (r_in/r)^2 on a sphere) with
//! the constant opacity kappa = kappa_f + kappa_s gives the radiative acceleration
//! a = rh_eddf g, a fixed fraction of gravity.  T is rh_temp, or with absorption
//! (kappa_p > 0) the Eddington radiation temperature (E/arad)^(1/4); rho follows from
//! dp/dx1 = -rho g (1 - rh_eddf), integrated by the midpoint rule.  E is the Eddington
//! solution 3 (F/c) (tau + 2/3) of the isothermal column, F the flux above, the bath at
//! x1min E(x1min).  The radiation relaxes to the scheme's own steady state with
//! <rad_vet>/clean_init = steady (which with force_reference = clean also relaxes the gas
//! internal energy); after it the gas becomes the static background (all cells, ghost
//! zones included), so the discrete state of both is the reference equilibrium.
//!
//! Gravity (user_srcs): rho g, on rho - rho_background with wellbalance_static, plus its
//! work; the background is in the restart file (pgen_rst_slabs).  Mesh x1 boundaries:
//! reflect.  <problem>/rh_seed adds a cell-to-cell x1 velocity of that amplitude after
//! the background is set (a lateral seed for stability tests; in 3-D the unsplit hydro
//! step needs cfl_number <= 1/3, else its grid-scale lateral mode grows).  History: gas
//! x1 momentum, kinetic energy, max |v|/c_s (all ranks), radiation E and F1/c (volume
//! sums).
//! <mhd> (read only when named): rh_b1, rh_b2, rh_b3, uniform on a Cartesian mesh;
//! on a sphere rh_b1 is the monopole B_r = rh_b1 (x1min/r)^2 (both force free), rh_b3 the
//! azimuthal field B_phi = rh_b3 r sin(theta)/x1min (spherical polar: not force free).
//! rh_rot = Omega: the velocity Omega z x x on top, after the background (L_z budgets).
//! wellbalance_dynamic (isothermal column, kappa_p = 0): Phi and the walked background
//! (AddWbGravity) in place of rho - rho_background; the gas of each column is relaxed to
//! the x1 walk of the steady radiation's force (wb_arad_md, one MeshBlock in x1).

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "pgen/pgen.hpp"
#include "utils/repro_sum.hpp"

void RadVetAmHistory(HistoryData *pdata, Mesh *pm);   // rad_vet_tests2.cpp

namespace {
Real rh_g = 1.0, rh_gm = 0.0, rh_eddf = 0.0, rh_seed = 0.0, rh_rot = 0.0, rh_t0 = 1.0;
int rh_round = 0;
bool rh_sph = false;

//! the fluid (<hydro> or <mhd>) and its well-balancing members
struct RhFluid : public radvet::FluidRef {
  DvceArray5D<Real> u0wb, w0wb;
  DvceFaceFld5D<Real> *wf = nullptr;
  bool wbs = false, wball = false;
  template <class F>
  void Wb(F *p) {
    u0wb = p->u0wb; w0wb = p->w0wb; wf = &(p->w0facewb);
    wbs = p->use_wellbalance_static; wball = p->wb_all;
  }
  explicit RhFluid(MeshBlockPack *pp) : radvet::FluidRef(radvet::FluidRef::Get(pp)) {
    if (pm != nullptr) {
      Wb(pm);
    } else if (ph != nullptr) {
      Wb(ph);
    }
  }
};

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
void RhseRelax(Mesh *pm);
void RhseBackground(Mesh *pm);
void RhseBgDerive(Mesh *pm);
DvceArray4D<Real> RhseEmag(Mesh *pm);
void RhseHist(HistoryData *pdata, Mesh *pm);
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationVetRhse
//! \brief vet_test = rhse, see the file header

void ProblemGenerator::RadiationVetRhse(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  radvet::RadiationVet *pvet = pmbp->pradvet;
  RhFluid fl(pmbp);
  if (!fl.on || !fl.eos.is_ideal) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<problem>/vet_test = rhse needs <hydro> or <mhd> with an "
              << "ideal or general EOS" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rh_sph = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
  rh_g = pin->GetOrAddReal("problem","rh_g",1.0);
  rh_gm = pin->GetOrAddReal("problem","rh_gm",1.0);
  rh_eddf = pin->GetOrAddReal("problem","rh_eddf",0.3);
  const Real rho0 = pin->GetOrAddReal("problem","rh_rho0",1.0);
  const Real t0 = pin->GetOrAddReal("problem","rh_temp",1.0);
  rh_seed = pin->GetOrAddReal("problem","rh_seed",0.0);
  rh_t0 = t0;
  auto key = [&](const char *name) {
    return pin->DoesParameterExist("problem", name) ? pin->GetReal("problem", name) : 0.0;
  };
  rh_rot = key("rh_rot");
  const bool wbd = fl.use_wellbalance_dynamic;
  const bool fmd = (pvet->force_ref == radvet::VET_FREF_WB_MD);
  if (wbd && (!fl.wball || !fmd || pvet->kappa_p > 0.0 || fl.use_etotgrav)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "vet_test = rhse with wellbalance_dynamic needs wb_dirs = "
              << "all, <rad_vet>/force_reference = wb_arad_md, kappa_p = 0, no etotgrav"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  user_srcs_func = wbd ? RhseWbGravity : RhseGravity;
  user_hist_func = RhseHist;
  // (read only when named) problem/am_hist: the energy / angular-momentum budget instead
  if (pin->DoesParameterExist("problem","am_hist") &&
      pin->GetBoolean("problem","am_hist")) {
    user_hist_func = RadVetAmHistory;
  }
  if (fl.wbs) {
    // the background travels in the restart file
    user_icpost_func = RhseBackground;
    pgen_rst_slabs = (fl.pm != nullptr) ? &fl.pm->u0wb : &fl.ph->u0wb;
    if (restart && pgen_rst_slabs_in.extent_int(1) == 5) {
      Kokkos::deep_copy(fl.u0wb, pgen_rst_slabs_in);
      RhseBgDerive(pmy_mesh_);
    }
  }

  const Real cl = pvet->c_light, ar = pvet->arad, kap = pvet->kappa_f + pvet->kappa_s;
  const Real gm1 = fl.eos.gamma - 1.0;
  auto eos = fl.eos;
  const bool gen = eos.IsGeneral();
  const bool teq = (pvet->kappa_p > 0.0);
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
  auto uh = fl.u0;
  auto u0 = pvet->u0;
  const Real efl = pvet->e_floor;
  pvet->iebath_x1min = RhseEdd(sph, a, a, b, rho0, t0, g, gm, ef, f0, kap, cl);
  if (wbd) {
    user_icpost_func = RhseRelax;
    pvet->fref_wsplit_ok = true;   // (RhseWbGravity gives the reference work)
    pvet->SetForceReference(pvet->arad_ref);
    auto pc = fl.phicc0, p1 = fl.phi0->x1f, p2 = fl.phi0->x2f, p3 = fl.phi0->x3f;
    auto xf = pmbp->pcoord->xx1f;
    par_for("rhse_phi", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real xl = size.d_view(m).x1min, xr = size.d_view(m).x1max;
      auto phi = [&](const Real x) {return sph ? -gm/x : g*x;};
      const Real c = phi(sph ? x1v(m,i) : CellCenterX(i-is, nx1, xl, xr));
      pc(m,k,j,i) = c;
      p1(m,k,j,i) = phi(sph ? xf(m,i) : LeftEdgeX(i-is, nx1, xl, xr));
      p1(m,k,j,i+1) = phi(sph ? xf(m,i+1) : LeftEdgeX(i+1-is, nx1, xl, xr));
      p2(m,k,j,i) = c;
      if (j == n2 - 1) p2(m,k,j+1,i) = c;
      p3(m,k,j,i) = c;
      if (k == n3 - 1) p3(m,k+1,j,i) = c;
    });
  }
  if (restart) return;
  if (fl.mhd) {
    const Real b1 = key("rh_b1"), b2 = key("rh_b2"), b3 = key("rh_b3");
    auto b0 = fl.pm->b0;
    auto xf = pmbp->pcoord->xx1f, x2f = pmbp->pcoord->xx2f;
    const bool sp = pmy_mesh_->use_spherical_polar;
    Kokkos::deep_copy(b0.x1f, b1);
    Kokkos::deep_copy(b0.x2f, sph ? 0.0 : b2);
    Kokkos::deep_copy(b0.x3f, sph ? 0.0 : b3);
    auto bx1 = b0.x1f, bx3 = b0.x3f;
    if (sph) {
      par_for("rhse_b", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        bx1(m,k,j,i) = b1*SQR(a/xf(m,i));
        bx1(m,k,j,i+1) = b1*SQR(a/xf(m,i+1));
        if (!sp) return;   // (axisymmetric: the mean of sin(theta) over the face)
        bx3(m,k,j,i) = b3*x1v(m,i)/a*sin(0.5*(x2f(m,j) + x2f(m,j+1)));
        if (n3 > 1) bx3(m,k+1,j,i) = bx3(m,k,j,i);
      });
    }
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
    uh(m,IEN,k,j,i) = gen ? eos.EnergyFromPressure(d, d*tt) : d*tt/gm1;
    u0(m,radvet::VET_E,k,j,i) = fmax(e, efl);
    u0(m,radvet::VET_F1,k,j,i) = f;
    u0(m,radvet::VET_F2,k,j,i) = 0.0;
    u0(m,radvet::VET_F3,k,j,i) = 0.0;
  });
  if (fl.mhd) {
    auto em = RhseEmag(pmy_mesh_);
    par_for("rhse_em", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      uh(m,IEN,k,j,i) += em(m,k,j,i);
    });
  }
  if (fl.wbs) RhseBackground(pmy_mesh_);
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn DvceArray4D<Real> RhseEmag
//! \brief |B|^2/2 of every cell as the MHD ConsToPrim forms it

DvceArray4D<Real> RhseEmag(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto uh = pmbp->pmhd->u0;
  const int n3 = uh.extent_int(2), n2 = uh.extent_int(3), n1 = uh.extent_int(4);
  DvceArray4D<Real> em("rhse_emag", pmbp->nmb_thispack, n3, n2, n1);
  auto b1f = pmbp->pmhd->b0.x1f, b2f = pmbp->pmhd->b0.x2f, b3f = pmbp->pmhd->b0.x3f;
  const bool sp = pm->use_spherical_polar, cs = pm->use_cubed_sphere;
  auto x1v = pmbp->pcoord->x1v, x1f = pmbp->pcoord->xx1f, x2v = pmbp->pcoord->x2v;
  auto x2f = pmbp->pcoord->xx2f, x3v = pmbp->pcoord->x3v, x3f = pmbp->pcoord->xx3f;
  par_for("rhse_emag", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    em(m,k,j,i) = radvet::VetEmagCell(b1f, b2f, b3f, sp, cs, x1v, x1f, x2v, x2f, x3v, x3f,
                                      m, k, j, i);
  });
  return em;
}

//----------------------------------------------------------------------------------------
//! \fn void RhseBackground
//! \brief the current gas, everywhere, as the static well-balanced background; the face
//! background is the mean of the two cells.  Then a cell-to-cell x1 velocity seed of
//! amplitude rh_seed on top (a lateral seed for stability tests; default 0)

void RhseBackground(Mesh *pm) {
  RhFluid fl(pm->pmb_pack);
  auto uh = fl.u0;
  if (fl.wbs) {
    Kokkos::deep_copy(fl.u0wb, Kokkos::subview(uh, Kokkos::ALL, std::make_pair(0, 5),
                                               Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  const int n3 = uh.extent_int(2), n2 = uh.extent_int(3), n1 = uh.extent_int(4);
  if (fl.wbs && fl.mhd) {   // (the background is the gas: no magnetic energy)
    auto em = RhseEmag(pm);
    auto ub = fl.u0wb;
    par_for("rhse_bgem", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      ub(m,IEN,k,j,i) -= em(m,k,j,i);
    });
  }
  if (fl.wbs) RhseBgDerive(pm);
  const Real sd = rh_seed;
  if (sd != 0.0) {
    par_for("rhse_seed", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dv = sd*sin(6.2831853*(0.618*i + 1.414*j + 2.236*k + 0.577*m));
      const Real m1 = uh(m,IM1,k,j,i);
      uh(m,IM1,k,j,i) += uh(m,IDN,k,j,i)*dv;
      uh(m,IEN,k,j,i) += 0.5*(SQR(uh(m,IM1,k,j,i)) - m1*m1)/uh(m,IDN,k,j,i);
    });
  }
  if (rh_rot == 0.0 || !rh_sph) return;
  // v = Omega z x x, covariant on the cell's tangents
  const Real om = rh_rot;
  const bool cs = pm->use_cubed_sphere;
  auto &id = pm->mb_indcs;
  const int js = id.js, ks = id.ks, nx2 = id.nx2, nx3 = id.nx3;
  auto &size = pm->pmb_pack->pmb->mb_size;
  auto &pan = pm->pmb_pack->pmb->mb_panel;
  auto x1v = pm->pmb_pack->pcoord->x1v, x2v = pm->pmb_pack->pcoord->x2v;
  par_for("rhse_rot", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i), r = x1v(m,i);
    Real v2 = 0.0, v3 = 0.0, st = 1.0;
    if (cs) {
      const auto &z = size.d_view(m);
      const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, z.x2min, z.x2max);
      const Real et = 0.25*M_PI*CellCenterX(k-ks, nx3, z.x3min, z.x3max);
      Real q[3], ta[3], tb[3];
      cubed_sphere::PanelToCart(pan.d_view(m), xi, et, q);
      cubed_sphere::PanelTangents(pan.d_view(m), xi, et, ta, tb);
      v2 = om*r*(q[0]*ta[1] - q[1]*ta[0]);
      v3 = om*r*(q[0]*tb[1] - q[1]*tb[0]);
      st = sqrt(q[0]*q[0] + q[1]*q[1]);
    } else {
      st = sin(x2v(m,j));
      v3 = om*r*st;
    }
    uh(m,IM2,k,j,i) += d*v2;
    uh(m,IM3,k,j,i) += d*v3;
    uh(m,IEN,k,j,i) += 0.5*d*SQR(om*r*st);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RhseBgDerive
//! \brief the primitive and face background of the conserved one (u0wb)

void RhseBgDerive(Mesh *pm) {
  RhFluid fl(pm->pmb_pack);
  auto u0wb = fl.u0wb, w0wb = fl.w0wb;
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
    auto wf = (dir == 0) ? fl.wf->x1f : ((dir == 1) ? fl.wf->x2f : fl.wf->x3f);
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
  RhFluid fl(pmbp);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nx1 = indcs.nx1;
  auto u0 = fl.u0;
  auto w0 = fl.w0;
  auto w0wb = fl.w0wb;
  auto x1v = pmbp->pcoord->x1v;
  auto &size = pmbp->pmb->mb_size;
  const bool wbs = fl.wbs, sph = rh_sph;
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
//! \fn void RhseRelax
//! \brief wellbalance_dynamic, after the steady start: rho (isothermal, T = rh_temp) of
//! each column from its x1min cell along the x1 walk of wb_arad_md, rho_{i+1} = rho_i
//! exp(-(Phi_{i+1} - Phi_i - P_1(i) - P_0(i+1))/T), eight rounds with the steady start;
//! then the rotation (RhseBackground)

void RhseRelax(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  RhFluid fl(pmbp);
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 != pm->mesh_indcs.nx1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "vet_test = rhse with wellbalance_dynamic needs one "
              << "MeshBlock in x1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int is = indcs.is, ie = indcs.ie;
  auto uh = fl.u0;
  auto pc = fl.phicc0;
  auto fr = pmbp->pradvet->cln_fref;
  const Real t0 = rh_t0, gm1 = fl.eos.gamma - 1.0;
  auto eos = fl.eos;
  const bool gen = eos.IsGeneral(), mhd = fl.mhd;
  auto em = mhd ? RhseEmag(pm) : DvceArray4D<Real>("rhse_e0", 1, 1, 1, 1);
  // (with the x2, x3 ghosts, which the next steady start reads)
  par_for("rhse_relax", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, uh.extent_int(2)-1,
          0, uh.extent_int(3)-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int i = is; i < ie; ++i) {
      const Real d = uh(m,IDN,k,j,i)*exp(-(pc(m,k,j,i+1) - pc(m,k,j,i) -
                     fr(m,4,k,j,i) - fr(m,3,k,j,i+1))/t0);
      uh(m,IDN,k,j,i+1) = d;
      uh(m,IEN,k,j,i+1) = (gen ? eos.EnergyFromPressure(d, d*t0) : d*t0/gm1) +
                          (mhd ? em(m,k,j,i+1) : 0.0);
    }
  });
  pm->pgen->user_icpost_again = (++rh_round < 8);
  if (!pm->pgen->user_icpost_again) RhseBackground(pm);
}

//----------------------------------------------------------------------------------------
//! \fn void RhseWbGravity
//! \brief wellbalance_dynamic: the walked background's gravity (AddWbGravity, all
//! directions), its work -rho v.grad Phi and the reference work rho a.v of wb_arad_md
//! (the module's acceleration; cubed sphere: on the hydro's basis)

void RhseWbGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  RhFluid fl(pmbp);
  radvet::RadiationVet *pvet = pmbp->pradvet;
  if (fl.pm != nullptr) {
    fl.pm->AddWbGravity(bdt, true);
  } else {
    fl.ph->AddWbGravity(bdt, true);
  }
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto u0 = fl.u0, w0 = fl.w0;
  auto x1v = pmbp->pcoord->x1v;
  auto &size = pmbp->pmb->mb_size;
  auto &pan = pmbp->pmb->mb_panel;
  const bool sph = rh_sph, cs = pm->use_cubed_sphere;
  const Real g = rh_g, gm = rh_gm;
  auto fr = pvet->cln_fref;
  const bool mdo = pvet->cln_fref_set && pvet->fref_wsplit;
  const bool acc = pvet->fref_wsplit && pvet->FrefWaccOn();
  auto wacc = acc ? pvet->FrefWacc() : DvceArray4D<Real>();
  const Real g0 = acc ? pvet->FrefWaccGam0(bdt, pm->dt) : 0.0;
  par_for("rhse_wbgrav", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ks+nx3-1, js,
          js+nx2-1, is, is+nx1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const auto &z = size.d_view(m);
    const Real x = sph ? x1v(m,i) : CellCenterX(i-is, nx1, z.x1min, z.x1max);
    const Real d = w0(m,IDN,k,j,i);
    u0(m,IEN,k,j,i) -= bdt*RhGravity(sph, g, gm, x)*d*w0(m,IVX,k,j,i);
    Real ar[3] = {0.0, 0.0, 0.0};
    for (int n = 0; n < 3 && mdo; ++n) ar[n] = fr(m,n,k,j,i);
    if (cs && mdo) {
      const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, z.x2min, z.x2max);
      const Real et = 0.25*M_PI*CellCenterX(k-ks, nx3, z.x3min, z.x3max);
      Real q[3], ta[3], tb[3];
      cubed_sphere::PanelToCart(pan.d_view(m), xi, et, q);
      cubed_sphere::PanelTangents(pan.d_view(m), xi, et, ta, tb);
      const Real fc[3] = {ar[0], ar[1], ar[2]};
      for (int n = 0; n < 3; ++n) ar[n] = 0.0;
      for (int n = 0; n < 3; ++n) {
        ar[0] += fc[n]*q[n]; ar[1] += fc[n]*ta[n]; ar[2] += fc[n]*tb[n];
      }
    }
    Real wr = 0.0;
    for (int n = 0; n < 3; ++n) wr += bdt*d*ar[n]*w0(m,IVX+n,k,j,i);
    u0(m,IEN,k,j,i) += wr;
    if (acc) wacc(m,k,j,i) = g0*wacc(m,k,j,i) + wr;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RhseHist
//! \brief gas x1 momentum and kinetic energy, max |v|/c_s (all ranks), radiation E, F1/c

void RhseHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  pdata->nhist = 5;
  pdata->label[0] = "vet_gas";
  pdata->label[1] = "KE_gas";
  pdata->label[2] = "vcs_max";
  pdata->label[3] = "E_rad";
  pdata->label[4] = "Y1_rad";
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  RhFluid fl(pmbp);
  auto w0 = fl.w0;
  auto ur = pmbp->pradvet->u0;
  auto vol = pmbp->pcoord->volume;
  auto &size = pmbp->pmb->mb_size;
  const bool sph = rh_sph;
  const Real gam = fl.eos.gamma, ic = 1.0/pmbp->pradvet->c_light;
  // (the sums in a fixed order, repro::HistSums: rank-count and layout invariant)
  Real s[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  array_sum::GlobalSum hs;
  repro::HistSums("rhse_hist", pmbp->nmb_thispack, nkji, nx1,
  KOKKOS_LAMBDA(const int q, array_sum::GlobalSum &a) {
    const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/nx1, i = is + q%nx1;
    const Real dv = sph ? vol(m,k,j,i) : size.d_view(m).dx1*size.d_view(m).dx2*
                                          size.d_view(m).dx3;
    const Real d = w0(m,IDN,k,j,i);
    const Real v2 = SQR(w0(m,IVX,k,j,i)) + SQR(w0(m,IVY,k,j,i)) + SQR(w0(m,IVZ,k,j,i));
    a.the_array[0] += dv*d*w0(m,IVX,k,j,i);
    a.the_array[1] += dv*0.5*d*v2;
    a.the_array[3] += dv*ur(m,radvet::VET_E,k,j,i);
    a.the_array[4] += dv*ic*ur(m,radvet::VET_F1,k,j,i);
  }, hs);
  Kokkos::parallel_reduce("rhse_vmax", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int q, Real &s2) {
    const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/nx1, i = is + q%nx1;
    const Real d = w0(m,IDN,k,j,i);
    const Real v2 = SQR(w0(m,IVX,k,j,i)) + SQR(w0(m,IVY,k,j,i)) + SQR(w0(m,IVZ,k,j,i));
    s2 = fmax(s2, std::sqrt(v2*d/(gam*(gam - 1.0)*w0(m,IEN,k,j,i))));
  }, Kokkos::Max<Real>(s[2]));
  for (int n : {0, 1, 3, 4}) s[n] = hs.the_array[n];
#if MPI_PARALLEL_ENABLED   // (the maximum on the master only: the history writer sums)
  MPI_Allreduce(MPI_IN_PLACE, &s[2], 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  if (global_variable::my_rank != 0) s[2] = 0.0;
#endif
  for (int n = 0; n < 5; ++n) pdata->hdata[n] = s[n];
}
}  // namespace
