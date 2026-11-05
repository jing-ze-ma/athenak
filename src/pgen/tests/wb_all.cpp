//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_all.cpp
//! \brief wball-1104 test problem (<problem>/pgen_name = wb_all): a hydrostatic
//! atmosphere in an arbitrary, grid-OBLIQUE potential on a Cartesian 2-D/3-D mesh, for
//! the all-direction dynamic well-balanced scheme (<hydro>/wb_dirs = all).
//!
//! <problem>/phi      = uniform  : Phi = -(g . x), g = (gx, gy, gz)
//!                    = pointtide: Phi = -GM/sqrt(|x - x0|^2 + soft^2)
//!                                       - kt/2 (2 X^2 - Y^2 - Z^2),  X = x - x0 etc.
//! <problem>/strat    = isothermal : p = p0 exp(-Phi/T0), rho = p/T0       (T0 = p0/rho0)
//!                    = polytropic : rho = rho0 (1 - Phi/h0)^n, p = K rho^(1+1/n),
//!                                   K = p0/rho0^(1+1/n), h0 = (n+1) p0/rho0
//!                    = expT       : T = T0 exp(beta Phi) (NOT a polytrope), p = p0 exp(
//!                                   -(1 - exp(-beta Phi))/(beta T0)), rho = p/T
//! (Phi = 0 at the reference point; every state is the exact continuous hydrostatic
//! solution at the cell centre.)
//! <problem>/flow     = none | shear | jet (phi = uniform only): a velocity along the
//!   equipotentials t = n x z-hat-ish (perpendicular to g in the g/x1 plane), a function
//!   of s = n . x only (n = g/|g|), so the flow is an EXACT steady solution:
//!   shear: v = U tanh((s - s0)/delta);  jet: v = U exp(-((s - s0)/delta)^2)
//! Gravity (user source): <hydro>/wb_dirs = all -> Hydro::AddWbGravity (x1 included);
//! otherwise the plain central -rho grad Phi from the face potentials (x1: the dynamic
//! x1 scheme's own source when wellbalance_dynamic + wb_x1).  Energy: v . (-rho grad
//! Phi), central.  All ghost zones: the exact solution (user BCs on every face).
//! Finalize prints max |v|/c_s, the L1 and max relative density error and the L1 error
//! of v against the exact solution: "WBALL_DIAG ...".

#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include <functional>
#include "pgen/pgen.hpp"

namespace {
struct WbAllV {
  int phi = 0;              // 0 uniform, 1 pointtide
  Real gx = 0.0, gy = -1.0, gz = 0.0;
  Real gm = 1.0, x0 = 0.0, y0 = -2.0, z0 = 0.0, soft = 0.0, kt = 0.0;
  int strat = 0;            // 0 isothermal, 1 polytropic, 2 expT
  Real p0 = 1.0, rho0 = 1.0, pn = 1.5, beta = 0.3, gm1 = 0.4;
  Real phiref = 0.0;        // Phi at the reference point (0, 0, 0)
  int flow = 0;             // 0 none, 1 shear, 2 jet
  Real U = 0.0, s0 = 0.0, delta = 0.1;
  Real om = 0.0;            // <problem>/omega: + the centrifugal -om^2 (x^2 + y^2)/2
  bool sp = false;          // spherical polar: (x1, x2, x3) = (r, theta, phi)
};
WbAllV wv_;

//! the fluid the problem runs with (<hydro> or <mhd>): the pieces wb_all needs
struct WbFl {
  DvceArray5D<Real> u0, w0, wbq0;
  bool all = false, dyn = false, wbx1 = false, mhd = false;
  DvceArray4D<Real> phicc;
  DvceFaceFld4D<Real> phi0;
  EOS_Data eos;
  std::function<void(Real, bool)> addgrav;
};
WbFl WbFluid(MeshBlockPack *pmbp) {
  WbFl f;
  if (pmbp->pmhd != nullptr) {
    auto *p = pmbp->pmhd;
    f.u0 = p->u0; f.w0 = p->w0; f.wbq0 = p->wbq0; f.all = p->wb_all;
    f.dyn = p->use_wellbalance_dynamic; f.wbx1 = p->use_wb_x1; f.mhd = true;
    f.phicc = p->phicc0; f.phi0 = p->phi0; f.eos = p->peos->eos_data;
    f.addgrav = [p](Real b, bool x1) {p->AddWbGravity(b, x1);};
  } else {
    auto *p = pmbp->phydro;
    f.u0 = p->u0; f.w0 = p->w0; f.wbq0 = p->wbq0; f.all = p->wb_all;
    f.dyn = p->use_wellbalance_dynamic; f.wbx1 = p->use_wb_x1;
    f.phicc = p->phicc0; f.phi0 = p->phi0; f.eos = p->peos->eos_data;
    f.addgrav = [p](Real b, bool x1) {p->AddWbGravity(b, x1);};
  }
  return f;
}
Real wb_bx_ = 0.0;   // <problem>/bx: a uniform x1 field (MHD, Cartesian)

KOKKOS_INLINE_FUNCTION
Real WbAllPhiRaw(const WbAllV &v, const Real x1, const Real x2, const Real x3) {
  Real x = x1, y = x2, z = x3;
  if (v.sp) {
    x = x1*sin(x2)*cos(x3); y = x1*sin(x2)*sin(x3); z = x1*cos(x2);
  }
  const Real cf = -0.5*v.om*v.om*(x*x + y*y);
  if (v.phi == 0) return -(v.gx*x + v.gy*y + v.gz*z) + cf;
  const Real X = x - v.x0, Y = y - v.y0, Z = z - v.z0;
  return -v.gm/sqrt(X*X + Y*Y + Z*Z + v.soft*v.soft) - 0.5*v.kt*(2.0*X*X - Y*Y - Z*Z) + cf;
}
KOKKOS_INLINE_FUNCTION
Real WbAllPhi(const WbAllV &v, const Real x, const Real y, const Real z) {
  return WbAllPhiRaw(v, x, y, z) - v.phiref;
}
//! exact hydrostatic (p, rho) at potential phi (relative to the reference point)
KOKKOS_INLINE_FUNCTION
void WbAllState(const WbAllV &v, const Real phi, Real &p, Real &d) {
  const Real t0 = v.p0/v.rho0;
  if (v.strat == 0) {
    p = v.p0*exp(-phi/t0); d = p/t0;
  } else if (v.strat == 1) {
    const Real h0 = (v.pn + 1.0)*t0;
    d = v.rho0*pow(fmax(1.0 - phi/h0, 1.0e-30), v.pn);
    p = v.p0*pow(d/v.rho0, 1.0 + 1.0/v.pn);
  } else {
    const Real b = v.beta;
    p = v.p0*exp(-(1.0 - exp(-b*phi))/(b*t0));
    d = p/(t0*exp(b*phi));
  }
}
//! the exact velocity (uniform phi only): along t = (-ny, nx, 0) for a g in the x-y plane
KOKKOS_INLINE_FUNCTION
void WbAllVel(const WbAllV &v, const Real x, const Real y, const Real z, Real &v1,
              Real &v2, Real &v3) {
  v1 = 0.0; v2 = 0.0; v3 = 0.0;
  if (v.flow == 0 || v.phi != 0 || v.sp) return;
  const Real gn = sqrt(v.gx*v.gx + v.gy*v.gy + v.gz*v.gz);
  const Real nx = v.gx/gn, ny = v.gy/gn, nz = v.gz/gn;
  const Real s = nx*x + ny*y + nz*z;
  const Real a = (s - v.s0)/v.delta;
  const Real u = (v.flow == 1) ? v.U*tanh(a) : v.U*exp(-a*a);
  // a unit vector perpendicular to n: (-ny, nx, 0) when n is not along z, else x-hat
  Real tx = -ny, ty = nx, tz = 0.0;
  const Real tn = sqrt(tx*tx + ty*ty);
  if (tn < 1.0e-12) {tx = 1.0; ty = 0.0;} else {tx /= tn; ty /= tn;}
  (void) nz;
  v1 = u*tx; v2 = u*ty; v3 = u*tz;
}


//! the cell-centre and face coordinates of one cell, Cartesian or spherical polar
struct WbAllG {
  bool sp = false, threed = false, cs = false;
  int is = 0, js = 0, ks = 0, nx1 = 1, nx2 = 1, nx3 = 1;
  DvceArray2D<Real> x1v, x2v, x3v, x1f, x2f, x3f;
  DualArray1D<RegionSize> size;
  DualArray1D<int> panel;
};
WbAllG WbAllGeo(Mesh *pm) {
  WbAllG g;
  auto *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  g.sp = pm->use_spherical_polar; g.threed = pm->three_d;
  g.is = indcs.is; g.js = indcs.js; g.ks = indcs.ks;
  g.nx1 = indcs.nx1; g.nx2 = indcs.nx2; g.nx3 = indcs.nx3;
  g.size = pmbp->pmb->mb_size;
  g.cs = pm->use_cubed_sphere;
  if (g.cs) g.panel = pmbp->pmb->mb_panel;
  if (g.sp || g.cs) {
    g.x1v = pmbp->pcoord->x1v; g.x2v = pmbp->pcoord->x2v; g.x3v = pmbp->pcoord->x3v;
    g.x1f = pmbp->pcoord->xx1f; g.x2f = pmbp->pcoord->xx2f; g.x3f = pmbp->pcoord->xx3f;
  }
  return g;
}
//! centre (c = 0) or left face (c = 1, 2, 3: of direction c) coordinates
KOKKOS_INLINE_FUNCTION
void WbAllX(const WbAllG &g, const int c, const int m, const int k, const int j,
            const int i, Real &x, Real &y, Real &z) {
  if (g.cs) {   // the Cartesian position (WbAllPhi then takes it as Cartesian)
    const Real r = (c == 1) ? g.x1f(m,i) : g.x1v(m,i);
    const Real xi = (c == 2) ? g.x2f(m,j) : g.x2v(m,j);
    const Real et = (c == 3) ? g.x3f(m,k) : g.x3v(m,k);
    Real q[3];
    cubed_sphere::PanelToCart(g.panel.d_view(m), 0.25*M_PI*xi, 0.25*M_PI*et, q);
    x = r*q[0]; y = r*q[1]; z = r*q[2];
    return;
  }
  if (g.sp) {
    x = (c == 1) ? g.x1f(m,i) : g.x1v(m,i);
    y = (c == 2) ? g.x2f(m,j) : g.x2v(m,j);
    z = (c == 3) ? g.x3f(m,k) : g.x3v(m,k);
    return;
  }
  const RegionSize &sz = g.size.d_view(m);
  x = (c == 1) ? LeftEdgeX(i - g.is, g.nx1, sz.x1min, sz.x1max)
               : CellCenterX(i - g.is, g.nx1, sz.x1min, sz.x1max);
  y = (c == 2) ? LeftEdgeX(j - g.js, g.nx2, sz.x2min, sz.x2max)
               : CellCenterX(j - g.js, g.nx2, sz.x2min, sz.x2max);
  if (!g.threed) {z = 0.0; return;}
  z = (c == 3) ? LeftEdgeX(k - g.ks, g.nx3, sz.x3min, sz.x3max)
               : CellCenterX(k - g.ks, g.nx3, sz.x3min, sz.x3max);
}
//! the plain central acceleration (orthonormal components) at a cell centre from the
//! potential at the face centres
KOKKOS_INLINE_FUNCTION
void WbAllAcc(const WbAllV &v, const WbAllG &g, const int m, const int k, const int j,
              const int i, Real &g1, Real &g2, Real &g3) {
  Real x, y, z, xl, yl, zl, xr, yr, zr;
  if (g.cs) {
    // -grad Phi by central differences in Cartesian, on (rhat, e_xi, e_eta) (covariant)
    WbAllX(g, 0, m, k, j, i, x, y, z);
    const Real r = sqrt(x*x + y*y + z*z), h = 1.0e-5*r;
    const Real gx = -(WbAllPhi(v, x + h, y, z) - WbAllPhi(v, x - h, y, z))/(2.0*h);
    const Real gy = -(WbAllPhi(v, x, y + h, z) - WbAllPhi(v, x, y - h, z))/(2.0*h);
    const Real gz = -(WbAllPhi(v, x, y, z + h) - WbAllPhi(v, x, y, z - h))/(2.0*h);
    Real e1[3], e2[3];
    cubed_sphere::PanelTangents(g.panel.d_view(m), 0.25*M_PI*g.x2v(m,j),
                                0.25*M_PI*g.x3v(m,k), e1, e2);
    g1 = (gx*x + gy*y + gz*z)/r;
    g2 = gx*e1[0] + gy*e1[1] + gz*e1[2];
    g3 = gx*e2[0] + gy*e2[1] + gz*e2[2];
    return;
  }
  WbAllX(g, 0, m, k, j, i, x, y, z);
  WbAllX(g, 1, m, k, j, i, xl, yl, zl);
  WbAllX(g, 1, m, k, j, i+1, xr, yr, zr);
  g1 = -(WbAllPhi(v, xr, y, z) - WbAllPhi(v, xl, y, z))/(xr - xl);
  WbAllX(g, 2, m, k, j, i, xl, yl, zl);
  WbAllX(g, 2, m, k, j+1, i, xr, yr, zr);
  g2 = -(WbAllPhi(v, x, yr, z) - WbAllPhi(v, x, yl, z))/((g.sp ? x : 1.0)*(yr - yl));
  g3 = 0.0;
  if (g.threed) {
    WbAllX(g, 3, m, k, j, i, xl, yl, zl);
    WbAllX(g, 3, m, k+1, j, i, xr, yr, zr);
    g3 = -(WbAllPhi(v, x, y, zr) - WbAllPhi(v, x, y, zl))/
         ((g.sp ? x*sin(y) : 1.0)*(zr - zl));
  }
}

void WbAllBCs(Mesh *pm);
void WbAllSrc(Mesh *pm, const Real bdt);
void WbAllFinal(ParameterInput *pin, Mesh *pm);

//! fill cells [kl..ku]x[jl..ju]x[il..iu] with the exact state
void WbAllFillRange(Mesh *pm, const int il, const int iu, const int jl, const int ju,
                    const int kl, const int ku, const bool only_ghosts) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto &bcs = pmbp->pmb->mb_bcs;
  auto u0 = WbFluid(pmbp).u0;
  const Real bx = wb_bx_;
  const WbAllV v = wv_;
  const bool threed = pm->three_d;
  const WbAllG gg = WbAllGeo(pm);
  par_for("wball_fill", DevExeSpace(), 0, pmbp->nmb_thispack - 1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (only_ghosts) {
      // a ghost cell of a PHYSICAL (user) boundary only: periodic/internal ghosts are
      // the exchange's
      bool g = false;
      if (i < is && bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) g = true;
      if (i > ie && bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) g = true;
      if (j < js && bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) g = true;
      if (j > je && bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) g = true;
      if (threed && k < ks &&
          bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) g = true;
      if (threed && k > ke &&
          bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) g = true;
      if (!g) return;
    }
    Real x, y, z;
    WbAllX(gg, 0, m, k, j, i, x, y, z);
    Real p, d, v1, v2, v3;
    WbAllState(v, WbAllPhi(v, x, y, z), p, d);
    WbAllVel(v, x, y, z, v1, v2, v3);
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*v2;
    u0(m,IM3,k,j,i) = d*v3;
    u0(m,IEN,k,j,i) = p/v.gm1 + 0.5*d*(v1*v1 + v2*v2 + v3*v3) + 0.5*bx*bx;
  });
}
}  // namespace

void ProblemGenerator::WbAll(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if ((pmbp->phydro == nullptr && pmbp->pmhd == nullptr) ||
      (pmy_mesh_->one_d && !pmy_mesh_->use_cubed_sphere)) {
    std::cout << "### FATAL ERROR in wb_all: needs <hydro> on a 2-D/3-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  WbFl fl = WbFluid(pmbp);
  if (!fl.eos.is_ideal || fl.eos.IsGeneral()) {
    std::cout << "### FATAL ERROR in wb_all: ideal gas only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  WbAllV v;
  const std::string ph_ = pin->GetOrAddString("problem", "phi", "uniform");
  v.phi = (ph_.compare("pointtide") == 0) ? 1 : 0;
  v.gx = pin->GetOrAddReal("problem", "gx", 0.0);
  v.gy = pin->GetOrAddReal("problem", "gy", -1.0);
  v.gz = pin->GetOrAddReal("problem", "gz", 0.0);
  v.gm = pin->GetOrAddReal("problem", "gm", 1.0);
  v.x0 = pin->GetOrAddReal("problem", "x0", 0.0);
  v.y0 = pin->GetOrAddReal("problem", "y0", -2.0);
  v.z0 = pin->GetOrAddReal("problem", "z0", 0.0);
  v.soft = pin->GetOrAddReal("problem", "soft", 0.0);
  v.kt = pin->GetOrAddReal("problem", "kt", 0.0);
  const std::string st = pin->GetOrAddString("problem", "strat", "isothermal");
  v.strat = (st.compare("polytropic") == 0) ? 1 : ((st.compare("expT") == 0) ? 2 : 0);
  v.p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  v.rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  v.pn = pin->GetOrAddReal("problem", "poly_n", 1.5);
  v.beta = pin->GetOrAddReal("problem", "beta", 0.3);
  v.gm1 = fl.eos.gamma - 1.0;
  v.om = pin->GetOrAddReal("problem", "omega", 0.0);
  v.sp = pmy_mesh_->use_spherical_polar;
  v.phiref = 0.0;
  // the reference point: the origin (Cartesian) or (r_ref, pi/2, 0) (spherical polar)
  v.phiref = v.sp ? WbAllPhiRaw(v, pin->GetOrAddReal("problem", "r_ref", 1.0), 0.5*M_PI,
                                0.0)
                  : (pmy_mesh_->use_cubed_sphere
                     ? WbAllPhiRaw(v, pin->GetOrAddReal("problem", "r_ref", 1.0), 0.0, 0.0)
                     : WbAllPhiRaw(v, 0.0, 0.0, 0.0));
  const std::string fl = pin->GetOrAddString("problem", "flow", "none");
  v.flow = (fl.compare("shear") == 0) ? 1 : ((fl.compare("jet") == 0) ? 2 : 0);
  v.U = pin->GetOrAddReal("problem", "U", 0.0);
  v.s0 = pin->GetOrAddReal("problem", "s0", 0.0);
  v.delta = pin->GetOrAddReal("problem", "delta", 0.1);
  wv_ = v;
  wb_bx_ = (pmbp->pmhd != nullptr) ? pin->GetOrAddReal("problem", "bx", 0.0) : 0.0;
  user_srcs_func = WbAllSrc;
  user_bcs_func = WbAllBCs;
  pgen_final_func = WbAllFinal;

  // the potential over the FULL arrays (the WB stencils reach into the ghost zones)
  if (fl.dyn) {
    auto &indcs = pmy_mesh_->mb_indcs;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3, ng = indcs.ng;
    const int n1 = nx1 + 2*ng, n2 = nx2 + 2*ng, n3 = (nx3 > 1) ? nx3 + 2*ng : 1;
    auto &size = pmbp->pmb->mb_size;
    auto pcc = fl.phicc;
    auto p1 = fl.phi0.x1f, p2 = fl.phi0.x2f, p3 = fl.phi0.x3f;
    const bool threed = pmy_mesh_->three_d;
    const WbAllG gg = WbAllGeo(pmy_mesh_);
    par_for("wball_phi", DevExeSpace(), 0, pmbp->nmb_thispack - 1, 0, n3 - 1, 0, n2 - 1,
            0, n1 - 1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real x, y, z, a, b, c;
      WbAllX(gg, 0, m, k, j, i, x, y, z);
      pcc(m,k,j,i) = WbAllPhi(v, x, y, z);
      WbAllX(gg, 1, m, k, j, i, a, b, c);
      p1(m,k,j,i) = WbAllPhi(v, a, b, c);
      WbAllX(gg, 2, m, k, j, i, a, b, c);
      p2(m,k,j,i) = WbAllPhi(v, a, b, c);
      if (threed) {
        WbAllX(gg, 3, m, k, j, i, a, b, c);
        p3(m,k,j,i) = WbAllPhi(v, a, b, c);
      }
      if (i == n1 - 1) {
        WbAllX(gg, 1, m, k, j, i+1, a, b, c);
        p1(m,k,j,i+1) = WbAllPhi(v, a, b, c);
      }
      if (j == n2 - 1) {
        WbAllX(gg, 2, m, k, j+1, i, a, b, c);
        p2(m,k,j+1,i) = WbAllPhi(v, a, b, c);
      }
      if (threed && k == n3 - 1) {
        WbAllX(gg, 3, m, k+1, j, i, a, b, c);
        p3(m,k+1,j,i) = WbAllPhi(v, a, b, c);
      }
    });
  }
  if (restart) return;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng, n2 = indcs.nx2 + 2*ng;
  const int n3 = pmy_mesh_->three_d ? indcs.nx3 + 2*ng : 1;
  WbAllFillRange(pmy_mesh_, 0, n1 - 1, 0, n2 - 1, 0, n3 - 1, false);
  if (pmbp->pmhd != nullptr) {   // the field: uniform B_x1 = bx (0: field-free MHD)
    const Real bx = wb_bx_;
    Kokkos::deep_copy(pmbp->pmhd->b0.x1f, bx);
    Kokkos::deep_copy(pmbp->pmhd->b0.x2f, 0.0);
    Kokkos::deep_copy(pmbp->pmhd->b0.x3f, 0.0);
    auto bcc = pmbp->pmhd->bcc0;
    par_for("wball_bcc", DevExeSpace(), 0, pmbp->nmb_thispack - 1, 0, n3 - 1, 0, n2 - 1,
            0, n1 - 1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      bcc(m,IBX,k,j,i) = bx; bcc(m,IBY,k,j,i) = 0.0; bcc(m,IBZ,k,j,i) = 0.0;
    });
  }
}

namespace {
void WbAllBCs(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng, n2 = indcs.nx2 + 2*ng;
  const int n3 = pm->three_d ? indcs.nx3 + 2*ng : 1;
  WbAllFillRange(pm, 0, n1 - 1, 0, n2 - 1, 0, n3 - 1, true);
}

void WbAllSrc(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  WbFl fl = WbFluid(pmbp);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto u0 = fl.u0;
  auto w0 = fl.w0;
  auto c1 = fl.wbq0;
  const bool all = fl.all;
  const bool dyn1 = fl.dyn && fl.wbx1 && !all;
  const bool threed = pm->three_d;
  const WbAllV v = wv_;
  const WbAllG gg = WbAllGeo(pm);
  const bool sp = gg.sp;
  auto vol = pmbp->pcoord->volume;
  auto a1 = pmbp->pcoord->area.x1f;
  if (all) fl.addgrav(bdt, true);
  par_for("wball_src", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    Real g1, g2, g3;
    WbAllAcc(v, gg, m, k, j, i, g1, g2, g3);
    if (!all) {
      if (dyn1 && sp) {
        const Real pc = c1(m,12,k,j,i);
        u0(m,IM1,k,j,i) += bdt*(a1(m,k,j,i+1)*(c1(m,13,k,j,i) - pc) +
                                a1(m,k,j,i)*(pc - c1(m,11,k,j,i)))/vol(m,k,j,i);
      } else if (dyn1) {
        u0(m,IM1,k,j,i) += bdt*(c1(m,13,k,j,i) - c1(m,11,k,j,i))/size.d_view(m).dx1;
      } else {
        u0(m,IM1,k,j,i) += bdt*d*g1;
      }
      u0(m,IM2,k,j,i) += bdt*d*g2;
      if (threed) u0(m,IM3,k,j,i) += bdt*d*g3;
    }
    u0(m,IEN,k,j,i) += bdt*d*(w0(m,IVX,k,j,i)*g1 + w0(m,IVY,k,j,i)*g2 +
                              w0(m,IVZ,k,j,i)*g3);
  });
}

void WbAllFinal(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  WbFl fl = WbFluid(pmbp);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  auto u0 = fl.u0;
  const bool threed = pm->three_d;
  const WbAllV v = wv_;
  const Real gam = fl.eos.gamma;
  const Real bx = wb_bx_;
  const WbAllG gg = WbAllGeo(pm);
  Real mach = 0.0, l1d = 0.0, mxd = 0.0, l1v = 0.0, ncell = 0.0, l1vw = 0.0, nw = 0.0;
  Kokkos::parallel_reduce("wball_diag",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a, Real &b, Real &c, Real &e, Real &f, Real &g,
                Real &h) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx%nx1 + is;
    Real x, y, z;
    WbAllX(gg, 0, m, k, j, i, x, y, z);
    Real p, d, v1, v2, v3;
    WbAllState(v, WbAllPhi(v, x, y, z), p, d);
    WbAllVel(v, x, y, z, v1, v2, v3);
    const Real dn = u0(m,IDN,k,j,i);
    const Real u1 = u0(m,IM1,k,j,i)/dn, u2 = u0(m,IM2,k,j,i)/dn, u3 = u0(m,IM3,k,j,i)/dn;
    const Real ek = 0.5*dn*(u1*u1 + u2*u2 + u3*u3);
    const Real pn = (gam - 1.0)*(u0(m,IEN,k,j,i) - ek - 0.5*bx*bx);
    const Real cs = sqrt(gam*fmax(pn, 1.0e-300)/dn);
    const Real dv = sqrt((u1-v1)*(u1-v1) + (u2-v2)*(u2-v2) + (u3-v3)*(u3-v3));
    a = fmax(a, dv/cs);
    const Real rd = fabs(dn - d)/d;
    b += rd; c = fmax(c, rd); e += dv/cs; f += 1.0;
    // the static part (exact |v| < 0.01 c_s): the spurious velocity where the flow is not
    const Real vv = sqrt(v1*v1 + v2*v2 + v3*v3);
    if (vv < 0.01*cs) {g += dv/cs; h += 1.0;}
  }, Kokkos::Max<Real>(mach), Kokkos::Sum<Real>(l1d), Kokkos::Max<Real>(mxd),
     Kokkos::Sum<Real>(l1v), Kokkos::Sum<Real>(ncell), Kokkos::Sum<Real>(l1vw),
     Kokkos::Sum<Real>(nw));
#if MPI_PARALLEL_ENABLED
  Real s[5] = {l1d, l1v, ncell, l1vw, nw}, x[2] = {mach, mxd};
  MPI_Allreduce(MPI_IN_PLACE, s, 5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, x, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  l1d = s[0]; l1v = s[1]; ncell = s[2]; l1vw = s[3]; nw = s[4]; mach = x[0]; mxd = x[1];
#endif
  if (global_variable::my_rank == 0) {
    std::printf("WBALL_DIAG t %.6e cycle %d max|dv|/cs %.6e L1|dv|/cs %.6e "
                "L1|drho|/rho %.6e max|drho|/rho %.6e static-part L1|dv|/cs %.6e "
                "(%.0f cells)\n", pm->time, pm->ncycle, mach, l1v/ncell, l1d/ncell, mxd,
                (nw > 0.0) ? l1vw/nw : 0.0, nw);
  }
}
}  // namespace
