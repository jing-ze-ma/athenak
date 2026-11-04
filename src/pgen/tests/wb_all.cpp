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
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
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
};
WbAllV wv_;

KOKKOS_INLINE_FUNCTION
Real WbAllPhiRaw(const WbAllV &v, const Real x, const Real y, const Real z) {
  if (v.phi == 0) return -(v.gx*x + v.gy*y + v.gz*z);
  const Real X = x - v.x0, Y = y - v.y0, Z = z - v.z0;
  return -v.gm/sqrt(X*X + Y*Y + Z*Z + v.soft*v.soft) - 0.5*v.kt*(2.0*X*X - Y*Y - Z*Z);
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
  if (v.flow == 0 || v.phi != 0) return;
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
  auto u0 = pmbp->phydro->u0;
  const WbAllV v = wv_;
  const bool threed = pm->three_d;
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
    const Real x = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real z = threed ? CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                        size.d_view(m).x3max) : 0.0;
    Real p, d, v1, v2, v3;
    WbAllState(v, WbAllPhi(v, x, y, z), p, d);
    WbAllVel(v, x, y, z, v1, v2, v3);
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*v2;
    u0(m,IM3,k,j,i) = d*v3;
    u0(m,IEN,k,j,i) = p/v.gm1 + 0.5*d*(v1*v1 + v2*v2 + v3*v3);
  });
}
}  // namespace

void ProblemGenerator::WbAll(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr || pmy_mesh_->one_d) {
    std::cout << "### FATAL ERROR in wb_all: needs <hydro> on a 2-D/3-D mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto *ph = pmbp->phydro;
  if (!ph->peos->eos_data.is_ideal || ph->peos->eos_data.IsGeneral()) {
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
  v.gm1 = ph->peos->eos_data.gamma - 1.0;
  v.phiref = 0.0;
  v.phiref = WbAllPhiRaw(v, 0.0, 0.0, 0.0);
  const std::string fl = pin->GetOrAddString("problem", "flow", "none");
  v.flow = (fl.compare("shear") == 0) ? 1 : ((fl.compare("jet") == 0) ? 2 : 0);
  v.U = pin->GetOrAddReal("problem", "U", 0.0);
  v.s0 = pin->GetOrAddReal("problem", "s0", 0.0);
  v.delta = pin->GetOrAddReal("problem", "delta", 0.1);
  wv_ = v;
  user_srcs_func = WbAllSrc;
  user_bcs_func = WbAllBCs;
  pgen_final_func = WbAllFinal;

  // the potential over the FULL arrays (the WB stencils reach into the ghost zones)
  if (ph->use_wellbalance_dynamic) {
    auto &indcs = pmy_mesh_->mb_indcs;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3, ng = indcs.ng;
    const int n1 = nx1 + 2*ng, n2 = nx2 + 2*ng, n3 = (nx3 > 1) ? nx3 + 2*ng : 1;
    auto &size = pmbp->pmb->mb_size;
    auto pcc = ph->phicc0;
    auto p1 = ph->phi0.x1f, p2 = ph->phi0.x2f, p3 = ph->phi0.x3f;
    const bool threed = pmy_mesh_->three_d;
    par_for("wball_phi", DevExeSpace(), 0, pmbp->nmb_thispack - 1, 0, n3 - 1, 0, n2 - 1,
            0, n1 - 1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      const Real y = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      const Real z = threed ? CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                          size.d_view(m).x3max) : 0.0;
      const Real xf = LeftEdgeX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      const Real yf = LeftEdgeX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      const Real zf = threed ? LeftEdgeX(k - ks, nx3, size.d_view(m).x3min,
                                         size.d_view(m).x3max) : 0.0;
      pcc(m,k,j,i) = WbAllPhi(v, x, y, z);
      p1(m,k,j,i) = WbAllPhi(v, xf, y, z);
      p2(m,k,j,i) = WbAllPhi(v, x, yf, z);
      p3(m,k,j,i) = WbAllPhi(v, x, y, zf);
      if (i == n1 - 1) {
        const Real xr = LeftEdgeX(i + 1 - is, nx1, size.d_view(m).x1min,
                                  size.d_view(m).x1max);
        p1(m,k,j,i+1) = WbAllPhi(v, xr, y, z);
      }
      if (j == n2 - 1) {
        const Real yr = LeftEdgeX(j + 1 - js, nx2, size.d_view(m).x2min,
                                  size.d_view(m).x2max);
        p2(m,k,j+1,i) = WbAllPhi(v, x, yr, z);
      }
      if (threed && k == n3 - 1) {
        const Real zr = LeftEdgeX(k + 1 - ks, nx3, size.d_view(m).x3min,
                                  size.d_view(m).x3max);
        p3(m,k+1,j,i) = WbAllPhi(v, x, y, zr);
      }
    });
  }
  if (restart) return;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng, n2 = indcs.nx2 + 2*ng;
  const int n3 = pmy_mesh_->three_d ? indcs.nx3 + 2*ng : 1;
  WbAllFillRange(pmy_mesh_, 0, n1 - 1, 0, n2 - 1, 0, n3 - 1, false);
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
  auto *ph = pmbp->phydro;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto c1 = ph->wbq0;
  const bool all = ph->wb_all;
  const bool dyn1 = ph->use_wellbalance_dynamic && ph->use_wb_x1 && !all;
  const bool threed = pm->three_d;
  const WbAllV v = wv_;
  if (all) ph->AddWbGravity(bdt, true);
  par_for("wball_src", DevExeSpace(), 0, pmbp->nmb_thispack - 1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    const Real dx1 = size.d_view(m).dx1, dx2 = size.d_view(m).dx2;
    const Real dx3 = threed ? size.d_view(m).dx3 : 1.0;
    const Real x = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real z = threed ? CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                        size.d_view(m).x3max) : 0.0;
    // the central analytic acceleration from the face potentials
    const Real g1 = -(WbAllPhi(v, x + 0.5*dx1, y, z) - WbAllPhi(v, x - 0.5*dx1, y, z))/dx1;
    const Real g2 = -(WbAllPhi(v, x, y + 0.5*dx2, z) - WbAllPhi(v, x, y - 0.5*dx2, z))/dx2;
    const Real g3 = threed ? -(WbAllPhi(v, x, y, z + 0.5*dx3) -
                               WbAllPhi(v, x, y, z - 0.5*dx3))/dx3 : 0.0;
    if (!all) {
      if (dyn1) {
        u0(m,IM1,k,j,i) += bdt*(c1(m,13,k,j,i) - c1(m,11,k,j,i))/dx1;
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
  auto *ph = pmbp->phydro;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  auto u0 = ph->u0;
  const bool threed = pm->three_d;
  const WbAllV v = wv_;
  const Real gam = ph->peos->eos_data.gamma;
  Real mach = 0.0, l1d = 0.0, mxd = 0.0, l1v = 0.0, ncell = 0.0, l1vw = 0.0, nw = 0.0;
  Kokkos::parallel_reduce("wball_diag",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a, Real &b, Real &c, Real &e, Real &f, Real &g,
                Real &h) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx%nx1 + is;
    const Real x = CellCenterX(i - is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(j - js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real z = threed ? CellCenterX(k - ks, nx3, size.d_view(m).x3min,
                                        size.d_view(m).x3max) : 0.0;
    Real p, d, v1, v2, v3;
    WbAllState(v, WbAllPhi(v, x, y, z), p, d);
    WbAllVel(v, x, y, z, v1, v2, v3);
    const Real dn = u0(m,IDN,k,j,i);
    const Real u1 = u0(m,IM1,k,j,i)/dn, u2 = u0(m,IM2,k,j,i)/dn, u3 = u0(m,IM3,k,j,i)/dn;
    const Real ek = 0.5*dn*(u1*u1 + u2*u2 + u3*u3);
    const Real pn = (gam - 1.0)*(u0(m,IEN,k,j,i) - ek);
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
