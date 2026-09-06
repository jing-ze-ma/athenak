//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wb_column.cpp
//! \brief A plane-parallel hydrostatic column for the dynamic well-balanced scheme, with
//! gravity along x1 or x3, for hydro, MHD and resistive MHD (ideal gas).
//!
//!   problem/axis     1 or 3: the vertical axis (gravity points to -axis)
//!   problem/g0       gravity at the bottom (z = xmin of that axis)
//!   problem/ap       > 0: point-mass form g(z) = g0 (ap/(ap + z))^2, else constant g0
//!   problem/rho0     density at the bottom
//!   problem/t0       T = p/rho at the bottom (code units)
//!   problem/tgrad    dT/dz: a LINEAR temperature profile, so neither the isothermal nor
//!                    the isentropic background is exact and the polytropic one is
//!   problem/b0       MHD only: uniform field along x2 (code units), 0 = none
//!   problem/user_srcs must be true (the gravity source lives here)
//!
//! The initial column is integrated on a fine host grid, dln p/dz = -g/T, and evaluated
//! in every cell INCLUDING the ghosts; the potential Phi(z) = int g dz (increasing
//! upward, the code's convention) is stored at every centroid and every face, ghosts
//! included.  The gravity source is -rho g along the axis, or under wellbalance_dynamic
//! the background's own pressure difference across the cell, from the same stencil the
//! reconstruction uses (the x1 cache when axis = 1 and wb_x1 is on; a direct walk
//! otherwise).  Use ix/ox_bc = user on the vertical axis (the column continued into the
//! ghosts, normal velocity mirrored) and periodic sides.

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "utils/wb_background.hpp"
#include "pgen.hpp"

void WBColumnGravity(Mesh *pm, Real bdt);
void WBColumnBC(Mesh *pm);
void WBColumnFinal(ParameterInput *pin, Mesh *pm);

namespace {
int axis_ = 1;
Real g0_ = 1.0, ap_ = 0.0;
// the column profile and everything the user boundary needs (it gets only a Mesh*)
DvceArray1D<Real> lnp_d_;
Real zlo_ = 0.0, dzf_ = 1.0, zmin_ = 0.0, t0_ = 1.0, tgrad_ = 0.0, gm1_ = 0.4, b0_ = 0.0;
int nfine_ = 0;
bool etotgrav_ = false;
KOKKOS_INLINE_FUNCTION Real GravAt(const Real g0, const Real ap, const Real z) {
  return (ap > 0.0) ? g0*SQR(ap/(ap + z)) : g0;
}
KOKKOS_INLINE_FUNCTION Real PotAt(const Real g0, const Real ap, const Real z) {
  return (ap > 0.0) ? g0*ap*(1.0 - ap/(ap + z)) : g0*z;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_srcs_func = WBColumnGravity;
  user_bcs_func = WBColumnBC;
  pgen_final_func = WBColumnFinal;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;

  axis_ = pin->GetOrAddInteger("problem", "axis", 1);
  g0_ = pin->GetOrAddReal("problem", "g0", 1.0);
  ap_ = pin->GetOrAddReal("problem", "ap", 0.0);
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real t0 = pin->GetOrAddReal("problem", "t0", 1.0);
  const Real tgrad = pin->GetOrAddReal("problem", "tgrad", 0.0);
  const Real b0 = pin->GetOrAddReal("problem", "b0", 0.0);
  const int axis = axis_;
  const Real g0 = g0_, ap = ap_;
  if (axis != 1 && axis != 3) {
    std::cout << "### FATAL ERROR in wb_column: problem/axis must be 1 or 3" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const Real gamma = is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                            : pmbp->phydro->peos->eos_data.gamma;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;

  // the column on a fine host grid over the mesh's vertical extent plus its ghosts
  const Real zmin = (axis == 1) ? pmy_mesh_->mesh_size.x1min : pmy_mesh_->mesh_size.x3min;
  const Real zmax = (axis == 1) ? pmy_mesh_->mesh_size.x1max : pmy_mesh_->mesh_size.x3max;
  const int nz = (axis == 1) ? pmy_mesh_->mesh_indcs.nx1 : pmy_mesh_->mesh_indcs.nx3;
  const Real dz = (zmax - zmin)/nz;
  const int nfine = 20000;
  const Real zlo = zmin - (ng + 1)*dz, zhi = zmax + (ng + 1)*dz;
  const Real dzf = (zhi - zlo)/(nfine - 1);
  DualArray1D<Real> lnp("lnp", nfine);
  // integrate from the bottom of the ACTIVE column, both ways
  int i0 = static_cast<int>((zmin - zlo)/dzf + 0.5);
  auto tprof = [&](const Real z) { return t0 + tgrad*(z - zmin); };
  lnp.h_view(i0) = std::log(rho0*t0);
  for (int i = i0+1; i < nfine; ++i) {
    const Real za = zlo + (i-1)*dzf, zb = zlo + i*dzf, zm = 0.5*(za + zb);
    lnp.h_view(i) = lnp.h_view(i-1) - dzf*GravAt(g0, ap, zm - zmin)/tprof(zm);
  }
  for (int i = i0-1; i >= 0; --i) {
    const Real za = zlo + i*dzf, zb = zlo + (i+1)*dzf, zm = 0.5*(za + zb);
    lnp.h_view(i) = lnp.h_view(i+1) + dzf*GravAt(g0, ap, zm - zmin)/tprof(zm);
  }
  lnp.template modify<HostMemSpace>();
  lnp.template sync<DevExeSpace>();
  auto lnp_d = lnp.d_view;
  const Real gm1 = gamma - 1.0;
  lnp_d_ = lnp.d_view; zlo_ = zlo; dzf_ = dzf; zmin_ = zmin; t0_ = t0; tgrad_ = tgrad;
  gm1_ = gm1; b0_ = b0; nfine_ = nfine; etotgrav_ = etotgrav;
  if (restart) return;

  // the state and the potentials, every cell and every face, ghosts included
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray4D<Real> ph2 = is_mhd ? pmbp->pmhd->phi0.x2f : pmbp->phydro->phi0.x2f;
  DvceArray4D<Real> ph3 = is_mhd ? pmbp->pmhd->phi0.x3f : pmbp->phydro->phi0.x3f;
  const bool have_phi = (etotgrav || wbdyn);
  par_for("wbcol_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
    const Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    const Real x1l = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
    const Real x3l = LeftEdgeX(k-ks, indcs.nx3, x3min, x3max);
    const Real x1r = LeftEdgeX(i+1-is, indcs.nx1, x1min, x1max);
    const Real x3r = LeftEdgeX(k+1-ks, indcs.nx3, x3min, x3max);
    const Real z = (axis == 1) ? x1v : x3v;
    // linear interpolation of ln p on the fine grid
    Real s = (z - zlo)/dzf;
    int ii = static_cast<int>(s);
    ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
    const Real f = s - ii;
    const Real p = exp(lnp_d(ii)*(1.0 - f) + lnp_d(ii+1)*f);
    const Real t = t0 + tgrad*(z - zmin);
    const Real d = p/t;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = p/gm1 + 0.5*b0*b0;
    if (have_phi) {
      const Real phi_c = PotAt(g0, ap, z - zmin);
      phicc(m,k,j,i) = phi_c;
      if (etotgrav) u0(m,IEN,k,j,i) += d*phi_c;
      if (axis == 1) {
        ph1(m,k,j,i) = PotAt(g0, ap, x1l - zmin);
        if (i == n1m1) ph1(m,k,j,i+1) = PotAt(g0, ap, x1r - zmin);
        ph2(m,k,j,i) = phi_c; ph3(m,k,j,i) = phi_c;
        if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
        if (k == n3m1) ph3(m,k+1,j,i) = phi_c;
      } else {
        ph3(m,k,j,i) = PotAt(g0, ap, x3l - zmin);
        if (k == n3m1) ph3(m,k+1,j,i) = PotAt(g0, ap, x3r - zmin);
        ph1(m,k,j,i) = phi_c; ph2(m,k,j,i) = phi_c;
        if (i == n1m1) ph1(m,k,j,i+1) = phi_c;
        if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
      }
    }
  });
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    par_for("wbcol_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      b.x1f(m,k,j,i) = 0.0;
      b.x2f(m,k,j,i) = b0;
      b.x3f(m,k,j,i) = 0.0;
      if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
      if (j == n2m1) b.x2f(m,k,j+1,i) = b0;
      if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WBColumnGravity
//! \brief the gravity source along the chosen axis; the well-balanced form under
//! wellbalance_dynamic.

void WBColumnGravity(Mesh *pm, Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;
  const bool wbx1 = is_mhd ? pmbp->pmhd->use_wb_x1 : pmbp->phydro->use_wb_x1;
  const WBOption wbo = is_mhd ? pmbp->pmhd->wb_option : pmbp->phydro->wb_option;
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray4D<Real> ph3 = is_mhd ? pmbp->pmhd->phi0.x3f : pmbp->phydro->phi0.x3f;
  DvceArray5D<Real> wbq0 = is_mhd ? pmbp->pmhd->wbq0 : pmbp->phydro->wbq0;
  const int axis = axis_;
  const Real g0 = g0_, ap = ap_;
  const bool use_cache = (axis == 1) && wbx1;

  par_for("wbcol_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
    const Real zmin = (axis == 1) ? x1min : x3min;   // NOTE: mesh xmin == block xmin here
    const Real z = (axis == 1) ? CellCenterX(i-is, indcs.nx1, x1min, x1max)
                               : CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    const Real dzc = (axis == 1) ? (x1max - x1min)/indcs.nx1 : (x3max - x3min)/indcs.nx3;
    const int iv = (axis == 1) ? IVX : IVZ;
    const int im = (axis == 1) ? IM1 : IM3;
    const Real d = w0(m,IDN,k,j,i);
    const Real gz = GravAt(g0, ap, z - zmin);
    Real src = -bdt*gz*d;
    if (!etotgrav) u0(m,IEN,k,j,i) += src*w0(m,iv,k,j,i);
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      if (use_cache) {
        WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      } else if (axis == 1) {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i), w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i), w0(m,IEN,k,j,i+1),
            phicc(m,k,j,i-1), ph1(m,k,j,i), phicc(m,k,j,i), ph1(m,k,j,i+1),
            phicc(m,k,j,i+1), d1, pl, d2, pr, d3);
      } else {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k-1,j,i), w0(m,IDN,k,j,i), w0(m,IDN,k+1,j,i),
            w0(m,IEN,k-1,j,i), w0(m,IEN,k,j,i), w0(m,IEN,k+1,j,i),
            phicc(m,k-1,j,i), ph3(m,k,j,i), phicc(m,k,j,i), ph3(m,k+1,j,i),
            phicc(m,k+1,j,i), d1, pl, d2, pr, d3);
      }
      // Cartesian: equal face areas, so the background's pressure difference over dz
      src = bdt*(pr - pl)/dzc;
    }
    u0(m,im,k,j,i) += src;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WBColumnBC
//! \brief the hydrostatic column continued into the ghost cells of the vertical axis,
//! with the normal velocity mirrored (a reflecting wall whose ghost state is the column
//! itself, so the well-balanced stencil of the wall cell sees a consistent profile).

void WBColumnBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie, ks = indcs.ks, ke = indcs.ke;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const int axis = axis_;
  const Real g0 = g0_, ap = ap_, zlo = zlo_, dzf = dzf_, zmin = zmin_, t0 = t0_;
  const Real tgrad = tgrad_, gm1 = gm1_, b0 = b0_;
  const int nfine = nfine_;
  const bool etotgrav = etotgrav_;
  auto lnp_d = lnp_d_;
  auto fill = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                          const int i, const int km, const int jm,
                                          const int im) {
    // (k,j,i) the ghost cell, (km,jm,im) its mirror active cell
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
    const Real z = (axis == 1) ? CellCenterX(i-is, indcs.nx1, x1min, x1max)
                               : CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    Real s = (z - zlo)/dzf;
    int ii = static_cast<int>(s);
    ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
    const Real f = s - ii;
    const Real p = exp(lnp_d(ii)*(1.0 - f) + lnp_d(ii+1)*f);
    const Real t = t0 + tgrad*(z - zmin);
    const Real d = p/t;
    const int iv = (axis == 1) ? IVX : IVZ;
    Real v1 = w0(m,IVX,km,jm,im), v2 = w0(m,IVY,km,jm,im), v3 = w0(m,IVZ,km,jm,im);
    if (iv == IVX) v1 = -v1; else v3 = -v3;
    w0(m,IDN,k,j,i) = d; w0(m,IEN,k,j,i) = p/gm1;
    w0(m,IVX,k,j,i) = v1; w0(m,IVY,k,j,i) = v2; w0(m,IVZ,k,j,i) = v3;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1; u0(m,IM2,k,j,i) = d*v2; u0(m,IM3,k,j,i) = d*v3;
    Real e = p/gm1 + 0.5*d*(v1*v1 + v2*v2 + v3*v3) + 0.5*b0*b0;
    if (etotgrav) e += d*PotAt(g0, ap, z - zmin);
    u0(m,IEN,k,j,i) = e;
  };
  if (axis == 1) {
    par_for("wbcol_bc_x1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, ng-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int n) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        fill(m, k, j, is-1-n, k, j, is+n);
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
        fill(m, k, j, ie+1+n, k, j, ie-n);
      }
    });
  } else {
    par_for("wbcol_bc_x3", DevExeSpace(), 0, nmb1, 0, n2m1, 0, n1m1, 0, ng-1,
    KOKKOS_LAMBDA(const int m, const int j, const int i, const int n) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == BoundaryFlag::user) {
        fill(m, ks-1-n, j, i, ks+n, j, i);
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == BoundaryFlag::user) {
        fill(m, ke+1+n, j, i, ke-n, j, i);
      }
    });
  }
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    par_for("wbcol_bc_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool gx1 = (i < is || i > ie) && axis == 1;
      const bool gx3 = (k < ks || k > ke) && axis == 3;
      if (gx1 || gx3) {
        b.x1f(m,k,j,i) = 0.0; b.x2f(m,k,j,i) = b0; b.x3f(m,k,j,i) = 0.0;
        if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
        if (j == n2m1) b.x2f(m,k,j+1,i) = b0;
        if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
      }
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WBColumnFinal
//! \brief releases the file-scope device view before Kokkos::finalize()

void WBColumnFinal(ParameterInput *pin, Mesh *pm) {
  lnp_d_ = DvceArray1D<Real>();
  return;
}
