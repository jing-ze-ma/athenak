//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file am_wb.cpp
//! \brief test problem (pgen_name = am_wb): a spherical shell of gas in a point-mass
//! potential, on the cubed sphere or spherical polar (3-D or 2-D sp), for the
//! angular-momentum transport correction (<hydro>/am_transport) together with the
//! geometric-source variants wellbalance_static[_reconst].
//!
//! Isothermal-profile hydrostatic atmosphere (ideal gas, p = rho cs2):
//!     rho(r) = rho0 exp[(gm/cs2) (1/r - 1/r_in)],   g = -gm/r^2 rhat (problem/gm, cs2,
//!     rho0; gravity applied as a user source when problem/user_srcs = true),
//! optionally rotating rigidly about z (problem/rot = Omega) with a density perturbation
//! rho (1 + pert sin^2(theta) sin(3 phi) [1 + cos(theta)]) (problem/pert): a closed shell
//! (reflecting radial walls; cs whole sphere, sp polar or periodic wedge) whose axial
//! angular momentum may change only through the transport (L_z budget), and a static
//! hydrostatic state (rot = pert = 0) for the well-balanced property.
//! Cell noise and flows of the merge_lat tests (read only when named): problem/pert_rand
//! = D, every cell's density times 10^(D (h - 1/2)), h a hash of its position; v_rand, a
//! hashed Cartesian velocity of that size per component; vr, v_r = vr sin(pi (r -
//! x1min)/(x1max - x1min)); vx, a uniform velocity along x; vortex = beta, the isentropic
//! vortex (axis z, centre vortex_x, vortex_y, radius vortex_w; needs gm = 0).
//! With <mhd> the field is the discrete curl of A = bnoise s(x) along x1 (rhat) on the x1
//! edges (problem/bnoise): B_1 = 0 and div B = 0 to round-off. <mhd> also runs on a
//! Cartesian mesh (gm = 0; periodic, unit box: s has its period; the basis is x, y, z).
//! With <hydro>/wellbalance_static the same profile is the static background (u0wb, w0wb,
//! w0facewb; radial faces at x1f, lateral faces at the cell-centre radius) and gravity
//! acts on rho - rho_background (as wb_atm).  ix1_bc / ox1_bc = user hold the ghosts at
//! the analytic profile (v = 0).
//! History (problem/user_hist = true, hst data_format %.17e): M, E (total incl. the
//! potential energy sum V rho Phi), Lx, Ly, Lz (sum V r_c x m, arm r_c = x1v rhat_c as
//! am_transport), KE, and KEin = KE over 1.1 x1min < r < 0.9 x1max (away from the walls).

#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

namespace {
struct AmWbVars {
  Real gm, cs2, rho0, rin, gamma;
};
AmWbVars amwb;

KOKKOS_INLINE_FUNCTION
Real AmWbRho(const AmWbVars &a, const Real r) {
  return a.rho0*exp((a.gm/a.cs2)*(1.0/r - 1.0/a.rin));
}

//! the cell direction (rhat) and its tangential basis (cs: unit panel tangents; sp:
//! theta-hat, phi-hat) of cell (m,k,j)
KOKKOS_INLINE_FUNCTION
void AmWbDir(const bool csm, const RegionSize &sz, const int pnl, const int js,
             const int ks, const int nx2, const int nx3, const Real th, const Real ph,
             const int k, const int j, Real *rh, Real *ta, Real *tb, Real &c) {
  if (csm) {
    const Real xi = 0.25*M_PI*CellCenterX(j-js, nx2, sz.x2min, sz.x2max);
    const Real eta = 0.25*M_PI*CellCenterX(k-ks, nx3, sz.x3min, sz.x3max);
    Real q[3];
    cubed_sphere::PanelToCart(pnl, xi, eta, q);
    const Real qn = 1.0/sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    for (int n = 0; n < 3; ++n) rh[n] = q[n]*qn;
    cubed_sphere::PanelTangents(pnl, xi, eta, ta, tb);
    c = ta[0]*tb[0] + ta[1]*tb[1] + ta[2]*tb[2];
  } else {
    rh[0] = sin(th)*cos(ph); rh[1] = sin(th)*sin(ph); rh[2] = cos(th);
    ta[0] = cos(th)*cos(ph); ta[1] = cos(th)*sin(ph); ta[2] = -sin(th);
    tb[0] = -sin(ph); tb[1] = cos(ph); tb[2] = 0.0;
    c = 0.0;
  }
}
}  // namespace

void AmWbGravity(Mesh *pm, const Real bdt);
void AmWbBCs(Mesh *pm);
void AmWbHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::AmWb()

void ProblemGenerator::AmWb(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  const bool mhd = pmbp->pmhd != nullptr;
  const bool cart = !(pmy_mesh_->use_cubed_sphere || pmy_mesh_->use_spherical_polar);
  if ((pmbp->phydro == nullptr && !mhd) || (cart && !mhd)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "am_wb needs <hydro> on the cubed sphere or spherical polar"
              << ", or <mhd>" << std::endl;
    exit(EXIT_FAILURE);
  }
  amwb.gm = pin->GetOrAddReal("problem", "gm", 1.0);
  amwb.cs2 = pin->GetOrAddReal("problem", "cs2", 0.5);
  amwb.rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  amwb.rin = cart ? 1.0 : pin->GetReal("mesh", "x1min");
  amwb.gamma = (mhd ? pmbp->pmhd->peos : pmbp->phydro->peos)->eos_data.gamma;
  user_srcs_func = AmWbGravity;
  user_bcs_func = AmWbBCs;
  user_hist_func = AmWbHistory;
  if (restart) return;

  const Real rot = pin->GetOrAddReal("problem", "rot", 0.0);
  const Real pert = pin->GetOrAddReal("problem", "pert", 0.0);
  auto key = [&](const char *name, const Real def) {
    return pin->DoesParameterExist("problem", name) ? pin->GetReal("problem", name) : def;
  };
  const Real drand = key("pert_rand", 0.0), vrand = key("v_rand", 0.0);
  const Real vr = key("vr", 0.0), vx0 = key("vx", 0.0), beta = key("vortex", 0.0);
  const Real xv = key("vortex_x", 0.0), yv = key("vortex_y", 0.0);
  const Real wv = key("vortex_w", 1.0);
  const bool noise = drand != 0.0 || vrand != 0.0 || vr != 0.0 || vx0 != 0.0 ||
                     beta != 0.0;
  const Real rout = pmy_mesh_->mesh_size.x1max;
  const AmWbVars a = amwb;
  const Real gm1 = a.gamma - 1.0;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const bool csm = pmy_mesh_->use_cubed_sphere;
  auto &size = pmbp->pmb->mb_size;
  auto &pan = pmbp->pmb->mb_panel;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto u0 = mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  par_for("am_wb_ic", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real rh[3], ta[3], tb[3], c;
    Real r;
    if (cart) {   // (the x, y, z basis; r = 1, rh = the position in the noise)
      const auto &z = size.d_view(m);
      r = 1.0;
      rh[0] = CellCenterX(i-is, nx1, z.x1min, z.x1max);
      rh[1] = CellCenterX(j-js, nx2, z.x2min, z.x2max);
      rh[2] = CellCenterX(k-ks, nx3, z.x3min, z.x3max);
      for (int n = 0; n < 3; ++n) {ta[n] = (n == 1); tb[n] = (n == 2);}
    } else {
      AmWbDir(csm, size.d_view(m), csm ? pan.d_view(m) : 0, js, ks, nx2, nx3,
              csm ? 0.0 : x2v(m,j), csm ? 0.0 : x3v(m,k), k, j, rh, ta, tb, c);
      r = x1v(m,i);
    }
    const Real e1[3] = {cart ? 1.0 : rh[0], cart ? 0.0 : rh[1], cart ? 0.0 : rh[2]};
    const Real st2 = rh[0]*rh[0] + rh[1]*rh[1];
    const Real ph = atan2(rh[1], rh[0]);
    Real rho = AmWbRho(a, r);
    rho *= 1.0 + pert*st2*sin(3.0*ph)*(1.0 + rh[2]);
    const Real pr = AmWbRho(a, r)*a.cs2;
    // v = Omega z x x (Cartesian), covariant on the cell basis
    const Real vx = -rot*r*rh[1], vy = rot*r*rh[0];
    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = cart ? rho*vx : 0.0;
    u0(m,IM2,k,j,i) = rho*(vx*ta[0] + vy*ta[1]);
    u0(m,IM3,k,j,i) = rho*(vx*tb[0] + vy*tb[1]);
    u0(m,IEN,k,j,i) = pr/gm1 + 0.5*rho*(vx*vx + vy*vy);
    if (!noise) return;
    Real p = pr, v[3] = {vx0 + vx, vy, 0.0};
    if (beta != 0.0) {
      const Real dx = (r*rh[0] - xv)/wv, dy = (r*rh[1] - yv)/wv;
      const Real g = sqrt(a.cs2)*beta*exp(0.5*(1.0 - dx*dx - dy*dy))/(2.0*M_PI);
      const Real t = 1.0 - gm1*g*g/(2.0*a.gamma*a.cs2);
      rho *= pow(t, 1.0/gm1);
      p = rho*a.cs2*t;
      v[0] -= g*dy;
      v[1] += g*dx;
    }
    auto hash = [&](const Real s) {
      const Real h = sin(r*(12.9898*rh[0] + 78.233*rh[1] + 37.719*rh[2]) + s)*43758.5453;
      return h - floor(h);
    };
    const Real vrad = cart ? vr : vr*sin(M_PI*(r - a.rin)/(rout - a.rin));
    for (int n = 0; n < 3; ++n) v[n] += vrad*e1[n] + vrand*(2.0*hash(1.0 + n) - 1.0);
    rho *= pow(10.0, drand*(hash(0.0) - 0.5));
    u0(m,IDN,k,j,i) = rho;
    u0(m,IM1,k,j,i) = rho*(v[0]*e1[0] + v[1]*e1[1] + v[2]*e1[2]);
    u0(m,IM2,k,j,i) = rho*(v[0]*ta[0] + v[1]*ta[1] + v[2]*ta[2]);
    u0(m,IM3,k,j,i) = rho*(v[0]*tb[0] + v[1]*tb[1] + v[2]*tb[2]);
    u0(m,IEN,k,j,i) = p/gm1 + 0.5*rho*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  });

  if (mhd) {
    const Real bn = key("bnoise", 0.0);
    auto b0 = pmbp->pmhd->b0;
    auto dl = pmbp->pcoord->dxedge.x1e;
    auto a2 = pmbp->pcoord->area.x2f, a3 = pmbp->pcoord->area.x3f;
    auto cc = pmbp->pcoord->cos_cell;
    auto x2f = pmbp->pcoord->xx2f, x3f = pmbp->pcoord->xx3f;
    const Real t0 = pmy_mesh_->mesh_size.x2min, t1 = pmy_mesh_->mesh_size.x2max;
    const Real nphi = 4.0*M_PI/(pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min);
    par_for("am_wb_b", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke+1, js, je+1, is,
            ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const auto &sz = size.d_view(m);
      auto circ = [&](const int kk, const int jj) {   // A dl along the x1 edge
        const Real y2 = LeftEdgeX(jj-js, nx2, sz.x2min, sz.x2max);
        const Real y3 = LeftEdgeX(kk-ks, nx3, sz.x3min, sz.x3max);
        if (cart) {
          const Real y1 = CellCenterX(i-is, nx1, sz.x1min, sz.x1max);
          return bn*sz.dx1*sin(4.0*M_PI*y2 + 2.0*M_PI*y1)*sin(6.0*M_PI*y3 + 1.0);
        }
        Real q[3] = {sin(x2f(m,jj))*cos(x3f(m,kk)), sin(x2f(m,jj))*sin(x3f(m,kk)),
                     cos(x2f(m,jj))};
        if (csm) cubed_sphere::PanelToCart(pan.d_view(m), 0.25*M_PI*y2, 0.25*M_PI*y3, q);
        Real w = 1.0;   // (a theta wedge with walls: no field through them)
        if (!csm) w = SQR(sin(M_PI*(x2f(m,jj) - t0)/(t1 - t0)));
        const Real p3 = csm ? sin(9.0*q[0] + 3.0*x1v(m,i))*sin(8.0*q[1] + 1.0)*
                              sin(7.0*q[2] + 2.0)
                            : sin(nphi*x3f(m,kk) + 3.0*x1v(m,i))*sin(7.0*q[2] + 2.0);
        return bn*w*dl(m,kk,jj,i)*p3;
      };
      const Real s2 = cart ? sz.dx1*sz.dx3 : a2(m,k,j,i);
      const Real s3 = cart ? sz.dx1*sz.dx2 : a3(m,k,j,i);
      if (j <= je && k <= ke) b0.x1f(m,k,j,i) = 0.0;
      if (i > ie) return;
      if (k <= ke) b0.x2f(m,k,j,i) = (s2 > 0.0) ? (circ(k+1,j) - circ(k,j))/s2 : 0.0;
      if (j <= je) b0.x3f(m,k,j,i) = (circ(k,j) - circ(k,j+1))/s3;
    });
    par_for("am_wb_em", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real b2 = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
      const Real b3 = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      if (csm) {
        const Real c = cc(m,k,j);
        b2 = (b2 + c*b3)/sqrt(1.0 - c*c);
      }
      u0(m,IEN,k,j,i) += 0.5*(b2*b2 + b3*b3);
    });
    return;
  }
  if (!pmbp->phydro->use_wellbalance_static) return;
  // static background over the full arrays (ghosts included)
  const int n1m1 = indcs.nx1 + 2*indcs.ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng - 1) : 0;
  auto u0wb = pmbp->phydro->u0wb;
  auto w0wb = pmbp->phydro->w0wb;
  auto wf1 = pmbp->phydro->w0facewb.x1f;
  auto wf2 = pmbp->phydro->w0facewb.x2f;
  auto wf3 = pmbp->phydro->w0facewb.x3f;
  auto x1f = pmbp->pcoord->xx1f;
  par_for("am_wb_bg", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    auto set = [&](const DvceArray5D<Real> &w, const int kk, const int jj, const int ii,
                   const Real r) {
      const Real d = AmWbRho(a, r);
      w(m,IDN,kk,jj,ii) = d;
      w(m,IM1,kk,jj,ii) = 0.0;
      w(m,IM2,kk,jj,ii) = 0.0;
      w(m,IM3,kk,jj,ii) = 0.0;
      w(m,IEN,kk,jj,ii) = d*a.cs2/gm1;
    };
    const Real r = x1v(m,i);
    set(u0wb, k, j, i, r);
    set(w0wb, k, j, i, r);
    set(wf1, k, j, i, x1f(m,i));
    if (i == n1m1) set(wf1, k, j, i+1, x1f(m,i+1));
    set(wf2, k, j, i, r);
    if (j == n2m1) set(wf2, k, j+1, i, r);
    set(wf3, k, j, i, r);
    if (k == n3m1) set(wf3, k+1, j, i, r);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmWbGravity
//! \brief g = -gm/r^2 rhat on the radial momentum (on rho - rho_bg with
//! wellbalance_static) and its work on the total energy (full rho)

void AmWbGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  if (pmbp->pmhd != nullptr) return;   // (gm = 0)
  auto u0 = pmbp->phydro->u0;
  auto w0 = pmbp->phydro->w0;
  auto w0wb = pmbp->phydro->w0wb;
  const bool wbs = pmbp->phydro->use_wellbalance_static;
  auto x1v = pmbp->pcoord->x1v;
  const Real gm = amwb.gm;
  par_for("am_wb_grav", DevExeSpace(), 0, pmbp->nmb_thispack-1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real g = -gm/(r*r);
    const Real d = w0(m,IDN,k,j,i);
    u0(m,IEN,k,j,i) += bdt*g*d*w0(m,IVX,k,j,i);
    u0(m,IM1,k,j,i) += bdt*g*(wbs ? (d - w0wb(m,IDN,k,j,i)) : d);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmWbBCs
//! \brief user x1 boundaries: the analytic hydrostatic profile at rest in the ghosts

void AmWbBCs(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng, is = indcs.is, ie = indcs.ie;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  auto &bcs = pmbp->pmb->mb_bcs;
  auto u0 = pmbp->phydro->u0;
  auto x1v = pmbp->pcoord->x1v;
  const AmWbVars a = amwb;
  const Real gm1 = a.gamma - 1.0;
  par_for("am_wb_bc", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, n3-1, 0, n2-1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int g) {
    for (int s = 0; s < 2; ++s) {
      const BoundaryFace f = s ? BoundaryFace::outer_x1 : BoundaryFace::inner_x1;
      if (bcs.d_view(m,f) != BoundaryFlag::user) continue;
      const int i = s ? (ie + 1 + g) : (is - 1 - g);
      const Real d = AmWbRho(a, x1v(m,i));
      u0(m,IDN,k,j,i) = d;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = d*a.cs2/gm1;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmWbHistory
//! \brief M, E (+ potential), Lx, Ly, Lz, KE, KEin (double precision sums)

void AmWbHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 7;
  const char *lab[7] = {"M", "Etot", "Lx", "Ly", "Lz", "KE", "KEin"};
  for (int n = 0; n < 7; ++n) {pdata->label[n] = lab[n];}
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const bool csm = pm->use_cubed_sphere;
  auto u0 = pmbp->pmhd != nullptr ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto vol = pmbp->pcoord->volume;
  auto x1v = pmbp->pcoord->x1v;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto &size = pmbp->pmb->mb_size;
  auto &pan = pmbp->pmb->mb_panel;
  const Real gm = amwb.gm;
  const bool cart = !(csm || pm->use_spherical_polar);
  const Real rlo = 1.1*pm->mesh_size.x1min, rhi = 0.9*pm->mesh_size.x1max;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const int nmkji = pmbp->nmb_thispack*nkji;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("AmWbHist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks; j += js; i += is;
    array_sum::GlobalSum hv;
    for (int n = 0; n < NHISTORY_VARIABLES; ++n) {hv.the_array[n] = 0.0;}
    if (cart) {   // M, E and the momentum in place of L
      const auto &z = size.d_view(m);
      const Real dv = z.dx1*z.dx2*z.dx3;
      for (int q = 0; q < 2; ++q) hv.the_array[q] = dv*u0(m,q == 0 ? IDN : IEN,k,j,i);
      for (int q = 0; q < 3; ++q) hv.the_array[2+q] = dv*u0(m,IM1+q,k,j,i);
      mb_sum += hv;
      return;
    }
    Real rh[3], ta[3], tb[3], c;
    AmWbDir(csm, size.d_view(m), csm ? pan.d_view(m) : 0, js, ks, nx2, nx3,
            csm ? 0.0 : x2v(m,j), csm ? 0.0 : x3v(m,k), k, j, rh, ta, tb, c);
    const Real r = x1v(m,i), v = vol(m,k,j,i), s2 = 1.0 - c*c;
    const Real d = u0(m,IDN,k,j,i);
    const Real m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
    const Real a2 = (m2 - c*m3)/s2, a3 = (m3 - c*m2)/s2;   // contravariant
    Real p[3];
    for (int q = 0; q < 3; ++q) {p[q] = u0(m,IM1,k,j,i)*rh[q] + a2*ta[q] + a3*tb[q];}
    const Real ke = 0.5*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/d;
    hv.the_array[0] = v*d;
    hv.the_array[1] = v*(u0(m,IEN,k,j,i) - d*gm/r);
    hv.the_array[2] = v*r*(rh[1]*p[2] - rh[2]*p[1]);
    hv.the_array[3] = v*r*(rh[2]*p[0] - rh[0]*p[2]);
    hv.the_array[4] = v*r*(rh[0]*p[1] - rh[1]*p[0]);
    hv.the_array[5] = v*ke;
    hv.the_array[6] = (r > rlo && r < rhi) ? v*ke : 0.0;
    mb_sum += hv;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  for (int n = 0; n < pdata->nhist; ++n) {pdata->hdata[n] = sum_this_mb.the_array[n];}
}
