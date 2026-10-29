//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ohm_drag.cpp
//! \brief ohmdrag-1029 (TEST): Ohmic (Lorentz) drag with a steep radial eta(r) on the
//! cubed sphere and on spherical polar, against a 1-D reference (ohmdrag_1029/ref1d.py).
//!
//! Shell r in [1,2] between two perfectly conducting walls at rest: either the
//! production wall (mesh x1 bc = reflect, default <mhd>/wall_conducting = full) or
//! (x1 bc = user) OdRadialBC, the dipole held exactly in the ghosts plus the induced
//! tangential field continued with (r B_t)' = 0.  An axial dipole
//! B_p = m(3 cos rhat - zhat)
//! /r^3, set as the discrete curl of A = m (-y, x, 0)/r^3, threads both walls.  Shellular
//! rotation v = Omega(r) zhat x r with Omega = om0 sin^4(pi (r - 1)) (Omega, Omega' and
//! Omega'' vanish on both walls).  eta(r) = eta_bot 10^(-ndec (1 + tanh((r - rt)/w))/2)
//! is written into Resistivity::eta_b once (<mhd>/ohmic_resistivity = constant never
//! touches eta_b again).
//!
//! KINEMATIC (problem/kinematic = true): a user source puts the momentum back to the
//! prescribed rotation after every stage, keeping the internal energy (E += dKE), so v is
//! held and B_p stays the dipole exactly (v_p = 0, B_p potential).  The induced field is
//! then exactly B_phi = b(r,t) sin cos (one l = 2 mode), u = r b obeying
//!     u_t = (eta u')' - 6 eta u / r^2 + 2 m Omega' / r,   u' = 0 on both walls.
//! The L_z the reset removes per stage, divided by beta*dt, is the Lorentz torque
//! density.
//! DYNAMIC (kinematic = false): the same state left free; the user history carries the
//! L_z / energy budget and an independent eta J^2 and s (J x B)_phi from the code's own
//! resistive edge EMF.
//!
//! problem/ref_file: the ref1d.py table at t = tlim; the final function writes the L1
//! errors of B_phi (faces), J and eta J (edges), the accumulated and instantaneous Joule
//! heating and the torque density, plus totals, to problem/err_file.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"
#include "outputs/outputs.hpp"
#include "pgen.hpp"

namespace {

Real od_m = 1.0, od_om0 = 0.05, od_etab = 0.1, od_ndec = 4.0, od_rt = 1.6, od_w = 0.06;
Real od_d0 = 1.0, od_p0 = 1.0;
bool od_kin = true;
bool od_snap = false;      // e_int snapshot taken (start of cycle 1)
char od_ref[1024] = "none", od_err[1024] = "ohm_drag-errs.dat";

// aux slots
enum {AVC1 = 0, AVC2, AVC3, AVT1, AVT2, AVT3, ALZ2, ALZ3, AEI0, ATRQ, ACOS, NAUX};

struct OdData {
  DvceFaceFld4D<Real> binit;
  DvceArray5D<Real> aux;
  OdData() : binit("od_binit", 1, 1, 1, 1), aux("od_aux", 1, 1, 1, 1, 1) {}
};
OdData *od = nullptr;

KOKKOS_INLINE_FUNCTION
Real OdOmega(const Real r, const Real om0) {
  const Real s = sin(M_PI*(r - 1.0));
  return om0*s*s*s*s;
}

KOKKOS_INLINE_FUNCTION
Real OdEta(const Real r, const Real etab, const Real ndec, const Real rt, const Real w) {
  return etab*pow(10.0, -ndec*0.5*(1.0 + tanh((r - rt)/w)));
}

//! position unit vector q and the unit tangent pair at native angular coordinates
//! (cs: xi, eta in radians on panel p; sp: theta, phi).  sp: e1 = thetahat, e2 = phihat.
KOKKOS_INLINE_FUNCTION
void OdGeo(const bool cs, const int p, const Real a2, const Real a3,
           Real q[3], Real e1[3], Real e2[3]) {
  if (cs) {
    cubed_sphere::PanelToCart(p, a2, a3, q);
    cubed_sphere::PanelTangents(p, a2, a3, e1, e2);
  } else {
    const Real st = sin(a2), ct = cos(a2), sp = sin(a3), cp = cos(a3);
    q[0] = st*cp; q[1] = st*sp; q[2] = ct;
    e1[0] = ct*cp; e1[1] = ct*sp; e1[2] = -st;
    e2[0] = -sp; e2[1] = cp; e2[2] = 0.0;
  }
}

//! face normals: cs nhat_xi = (e1 - c e2)/s, nhat_eta = (e2 - c e1)/s; sp the tangents
KOKKOS_INLINE_FUNCTION
void OdNormals(const Real e1[3], const Real e2[3], Real n1[3], Real n2[3]) {
  const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
  const Real s = sqrt(1.0 - c*c);
  for (int d=0; d<3; ++d) {
    n1[d] = (e1[d] - c*e2[d])/s;
    n2[d] = (e2[d] - c*e1[d])/s;
  }
}

//! reference fields at (r, cos theta, sin theta) from b, u' = (r b)'
KOKKOS_INLINE_FUNCTION
void OdRefJ(const Real r, const Real b, const Real up, const Real q[3],
            Real jc[3]) {
  // J = J_r rhat + J_th thetahat, J_r = b (3 cos^2 - 1)/r, J_th = -u' sin cos / r
  const Real ct = q[2];
  const Real st = sqrt(fmax(0.0, 1.0 - ct*ct));
  const Real jr = b*(3.0*ct*ct - 1.0)/r;
  const Real jt = -up*st*ct/r;
  Real th[3];
  if (st > 1.0e-14) {
    th[0] = ct*q[0]/st; th[1] = ct*q[1]/st; th[2] = -st;
  } else {
    th[0] = 1.0; th[1] = 0.0; th[2] = 0.0;
  }
  for (int d=0; d<3; ++d) jc[d] = jr*q[d] + jt*th[d];
}

// reference table on a uniform cell-centred grid in r
struct OdRef {
  DvceArray2D<Real> tab;   // (5, N): r b up Ib Id
  int n = 0;
  Real h = 0.0;
  Real qtot = 0.0, torque = 0.0;
};

KOKKOS_INLINE_FUNCTION
Real OdInterp(const DvceArray2D<Real> &tab, const int n, const Real h, const int col,
              const Real r) {
  Real x = (r - 1.0)/h - 0.5;
  if (x < 0.0) x = 0.0;
  if (x > n - 1.0) x = n - 1.0;
  int i0 = static_cast<int>(x);
  if (i0 > n - 2) i0 = n - 2;
  const Real f = x - i0;
  return (1.0 - f)*tab(col, i0) + f*tab(col, i0 + 1);
}

void OdSource(Mesh *pm, const Real bdt);
void OdRadialBC(Mesh *pm);
void OdHistory(HistoryData *pdata, Mesh *pm);
void OdFinal(ParameterInput *pin, Mesh *pm);
void OdResistEMF(MeshBlockPack *pmbp);

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->pmhd->presist == nullptr) {
    std::cout << "### FATAL ERROR ohm_drag: needs <mhd> with ohmic_resistivity"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const bool cs = pmy_mesh_->use_cubed_sphere;
  if (!cs && !pmy_mesh_->use_spherical_polar) {
    std::cout << "### FATAL ERROR ohm_drag: cubed sphere or spherical polar only"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  od_m = pin->GetOrAddReal("problem", "m_dip", 1.0);
  od_om0 = pin->GetOrAddReal("problem", "om0", 0.05);
  od_etab = pin->GetOrAddReal("problem", "eta_bot", 0.1);
  od_ndec = pin->GetOrAddReal("problem", "ndec", 4.0);
  od_rt = pin->GetOrAddReal("problem", "rt", 1.6);
  od_w = pin->GetOrAddReal("problem", "w", 0.06);
  od_d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  od_p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  od_kin = pin->GetOrAddBoolean("problem", "kinematic", true);
  std::snprintf(od_ref, sizeof(od_ref), "%s",
                pin->GetOrAddString("problem", "ref_file", "none").c_str());
  std::snprintf(od_err, sizeof(od_err), "%s",
                pin->GetOrAddString("problem", "err_file", "ohm_drag-errs.dat").c_str());
  user_srcs_func = OdSource;
  user_bcs_func = OdRadialBC;
  user_hist_func = OdHistory;
  pgen_final_func = OdFinal;
  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng, n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto &size = pmbp->pmb->mb_size;
  auto &mbp = pmbp->pmb->mb_panel;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  const Real sc = cs ? 0.25*M_PI : 1.0;
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  const Real m_ = od_m, om0 = od_om0, d0 = od_d0, p0 = od_p0;
  const Real etab = od_etab, ndec = od_ndec, rt = od_rt, w = od_w;

  od = new OdData();
  Kokkos::realloc(od->aux, nmb, static_cast<int>(NAUX), n3, n2, n1);
  auto aux = od->aux;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &eta_b = pmbp->pmhd->presist->eta_b;

  // cell state, aux geometry and eta_b, ghosts included
  par_for("od_cells", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real a2 = sc*CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real a3 = sc*CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    const int p = cs ? mbp.d_view(m) : 0;
    Real q[3], e1[3], e2[3];
    OdGeo(cs, p, a2, a3, q, e1, e2);
    const Real r = x1v(m, i);
    const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
    const Real om = OdOmega(r, om0);
    const Real v[3] = {-om*r*q[1], om*r*q[0], 0.0};
    const Real d1 = v[0]*e1[0] + v[1]*e1[1] + v[2]*e1[2];
    const Real d2 = v[0]*e2[0] + v[1]*e2[1] + v[2]*e2[2];
    const Real det = 1.0 - c*c;
    const Real va = (d1 - c*d2)/det, vb = (d2 - c*d1)/det;
    // lever arms: Lz = a z1 + b z2, (a,b) = G^{-1} (m2, m3)
    const Real z1 = r*(q[0]*e1[1] - q[1]*e1[0]);
    const Real z2 = r*(q[0]*e2[1] - q[1]*e2[0]);
    aux(m, AVC1, k, j, i) = 0.0;
    aux(m, AVC2, k, j, i) = d1;
    aux(m, AVC3, k, j, i) = d2;
    aux(m, AVT1, k, j, i) = 0.0;
    aux(m, AVT2, k, j, i) = va;
    aux(m, AVT3, k, j, i) = vb;
    aux(m, ALZ2, k, j, i) = (z1 - c*z2)/det;
    aux(m, ALZ3, k, j, i) = (z2 - c*z1)/det;
    aux(m, AEI0, k, j, i) = p0/gm1;
    aux(m, ATRQ, k, j, i) = 0.0;
    aux(m, ACOS, k, j, i) = c;
    eta_b(m, k, j, i) = OdEta(r, etab, ndec, rt, w);
    const Real ek = 0.5*d0*(d1*va + d2*vb);
    const Real b2 = m_*m_*(3.0*q[2]*q[2] + 1.0)/(r*r*r*r*r*r);
    w0(m, IDN, k, j, i) = d0;
    w0(m, IVX, k, j, i) = 0.0;
    w0(m, IVY, k, j, i) = va;
    w0(m, IVZ, k, j, i) = vb;
    w0(m, IEN, k, j, i) = p0/gm1;
    u0(m, IDN, k, j, i) = d0;
    u0(m, IM1, k, j, i) = 0.0;
    u0(m, IM2, k, j, i) = d0*d1;
    u0(m, IM3, k, j, i) = d0*d2;
    u0(m, IEN, k, j, i) = p0/gm1 + ek + 0.5*b2;
  });

  // faces: discrete curl of A = m (-y, x, 0)/r^3 on the CT edges (A.rhat = 0)
  auto &b0f = pmbp->pmhd->b0;
  auto &ar1 = pmbp->pcoord->area.x1f;
  auto &ar2 = pmbp->pcoord->area.x2f;
  auto &ar3 = pmbp->pcoord->area.x3f;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;
  par_for("od_faces", DevExeSpace(), 0, nmb-1, 0, n3, 0, n2, 0, n1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const int p = cs ? mbp.d_view(m) : 0;
    const Real x2mn = size.d_view(m).x2min, x2mx = size.d_view(m).x2max;
    const Real x3mn = size.d_view(m).x3min, x3mx = size.d_view(m).x3max;
    // A on the edge's own unit tangent; along2: the x2 edge (x2 centre, x3 face)
    auto Aedge = [&](const int ii, const int jj, const int kk, const bool along2) {
      const Real rf = x1f(m, ii);
      const Real a2 = sc*(along2 ? CellCenterX(jj-js, nx2, x2mn, x2mx)
                                 : LeftEdgeX(jj-js, nx2, x2mn, x2mx));
      const Real a3 = sc*(along2 ? LeftEdgeX(kk-ks, nx3, x3mn, x3mx)
                                 : CellCenterX(kk-ks, nx3, x3mn, x3mx));
      Real q[3], e1[3], e2[3];
      OdGeo(cs, p, a2, a3, q, e1, e2);
      const Real ax = -m_*q[1]/(rf*rf), ay = m_*q[0]/(rf*rf);
      const Real *t = along2 ? e1 : e2;
      return ax*t[0] + ay*t[1];
    };
    if (j < n2 && k < n3) {
      const Real a = ar1(m, k, j, i);
      b0f.x1f(m, k, j, i) = (a > 0.0) ?
          (dxe3(m, k, j+1, i)*Aedge(i, j+1, k, false)
         - dxe3(m, k, j, i)*Aedge(i, j, k, false)
         - dxe2(m, k+1, j, i)*Aedge(i, j, k+1, true)
         + dxe2(m, k, j, i)*Aedge(i, j, k, true))
          /a : 0.0;
    }
    if (i < n1 && k < n3) {
      const Real a = ar2(m, k, j, i);
      b0f.x2f(m, k, j, i) = (a > 1.0e-14) ?
          -(dxe3(m, k, j, i+1)*Aedge(i+1, j, k, false)
          - dxe3(m, k, j, i)*Aedge(i, j, k, false))
          /a : 0.0;
    }
    if (i < n1 && j < n2) {
      const Real a = ar3(m, k, j, i);
      b0f.x3f(m, k, j, i) = (a > 0.0) ?
          (dxe2(m, k, j, i+1)*Aedge(i+1, j, k, true)
         - dxe2(m, k, j, i)*Aedge(i, j, k, true))
          /a : 0.0;
    }
  });
  od->binit.x1f = DvceArray4D<Real>("od_b1", nmb, n3, n2, n1+1);
  od->binit.x2f = DvceArray4D<Real>("od_b2", nmb, n3, n2+1, n1);
  od->binit.x3f = DvceArray4D<Real>("od_b3", nmb, n3+1, n2, n1);
  Kokkos::deep_copy(od->binit.x1f, b0f.x1f);
  Kokkos::deep_copy(od->binit.x2f, b0f.x2f);
  Kokkos::deep_copy(od->binit.x3f, b0f.x3f);
  od_snap = false;
  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! kinematic hold (after every stage) + e_int snapshot at the start of cycle 1

void OdSource(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto aux = od->aux;
  auto &u0 = pmbp->pmhd->u0;
  auto &w0 = pmbp->pmhd->w0;
  if (!od_snap && pm->ncycle >= 1) {
    par_for("od_snap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      aux(m, AEI0, k, j, i) = w0(m, IEN, k, j, i);
    });
    od_snap = true;
  }
  if (!od_kin) return;
  const Real ibdt = 1.0/bdt;
  par_for("od_hold", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real d = u0(m, IDN, k, j, i);
    const Real c = aux(m, ACOS, k, j, i);
    const Real m1 = u0(m, IM1, k, j, i), m2 = u0(m, IM2, k, j, i);
    const Real m3 = u0(m, IM3, k, j, i);
    const Real ekold = 0.5*(m1*m1 + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))/d;
    const Real n2 = d*aux(m, AVC2, k, j, i), n3 = d*aux(m, AVC3, k, j, i);
    const Real eknew = 0.5*d*(aux(m, AVC2, k, j, i)*aux(m, AVT2, k, j, i)
                            + aux(m, AVC3, k, j, i)*aux(m, AVT3, k, j, i));
    aux(m, ATRQ, k, j, i) = (aux(m, ALZ2, k, j, i)*(m2 - n2)
                           + aux(m, ALZ3, k, j, i)*(m3 - n3))*ibdt;
    u0(m, IM1, k, j, i) = 0.0;
    u0(m, IM2, k, j, i) = n2;
    u0(m, IM3, k, j, i) = n3;
    u0(m, IEN, k, j, i) += eknew - ekold;
  });
}

//----------------------------------------------------------------------------------------
//! problem/wall_bc = exact (mesh ix1_bc = ox1_bc = user): the radial ghosts carry the
//! dipole exactly (x1 ghost faces = the initial faces) plus the induced tangential field
//! continued with (r B_t)' = 0 (conducting wall at rest, E_theta = eta J_theta = 0),
//! b_g = binit_g + (r_m/r_g)(b_m - binit_m).  Hydro ghosts: rho = d0, the held rotation
//! (kinematic) or the mirror velocity (v_r odd), and e_int = the mirror cell's; the ghost
//! u0(IEN) reuses the magnetic energy the last C2P subtracted there (u0 - w0 - e_k), so
//! the next C2P returns exactly the mirror e_int.

void OdRadialBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  auto &u0 = pmbp->pmhd->u0;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto bi = od->binit;
  auto aux = od->aux;
  const bool kin = od_kin;
  const Real d0 = od_d0;
  par_for("od_rbc", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, n3, 0, n2, 0, ng-1,
  KOKKOS_LAMBDA(int m, int k, int j, int n) {
    for (int side=0; side<2; ++side) {
      if (side == 0 &&
          mb_bcs.d_view(m, BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
      if (side == 1 &&
          mb_bcs.d_view(m, BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
      const int ig = (side == 0) ? is - 1 - n : ie + 1 + n;
      const int im = (side == 0) ? is + n : ie - n;
      const Real rat = x1v(m, im)/x1v(m, ig);
      if (k < n3) {
        b0.x2f(m, k, j, ig) = bi.x2f(m, k, j, ig)
                            + rat*(b0.x2f(m, k, j, im) - bi.x2f(m, k, j, im));
      }
      if (j < n2) {
        b0.x3f(m, k, j, ig) = bi.x3f(m, k, j, ig)
                            + rat*(b0.x3f(m, k, j, im) - bi.x3f(m, k, j, im));
      }
      if (j < n2 && k < n3) {
        const int ifc = (side == 0) ? ig : ig + 1;
        b0.x1f(m, k, j, ifc) = bi.x1f(m, k, j, ifc);
        // hydro ghost
        const Real ekold = 0.5*(u0(m, IM1, k, j, ig)*w0(m, IVX, k, j, ig)
                              + u0(m, IM2, k, j, ig)*w0(m, IVY, k, j, ig)
                              + u0(m, IM3, k, j, ig)*w0(m, IVZ, k, j, ig));
        const Real emag = u0(m, IEN, k, j, ig) - w0(m, IEN, k, j, ig) - ekold;
        const Real e = w0(m, IEN, k, j, im);
        Real v1, v2, v3;
        if (kin) {
          v1 = 0.0; v2 = aux(m, AVT2, k, j, ig); v3 = aux(m, AVT3, k, j, ig);
        } else {
          v1 = -w0(m, IVX, k, j, im);
          v2 = w0(m, IVY, k, j, im);
          v3 = w0(m, IVZ, k, j, im);
        }
        const Real c = aux(m, ACOS, k, j, ig);
        const Real c2 = v2 + c*v3, c3 = v3 + c*v2;
        w0(m, IDN, k, j, ig) = d0;
        w0(m, IVX, k, j, ig) = v1;
        w0(m, IVY, k, j, ig) = v2;
        w0(m, IVZ, k, j, ig) = v3;
        w0(m, IEN, k, j, ig) = e;
        u0(m, IDN, k, j, ig) = d0;
        u0(m, IM1, k, j, ig) = d0*v1;
        u0(m, IM2, k, j, ig) = d0*c2;
        u0(m, IM3, k, j, ig) = d0*c3;
        u0(m, IEN, k, j, ig) = e + emag + 0.5*d0*(v1*v1 + v2*c2 + v3*c3);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! the code's own resistive edge EMF eta J . that from the current b0 (walls zeroed as in
//! the stage), left in presist->efld_resist

void OdResistEMF(MeshBlockPack *pmbp) {
  auto *pr = pmbp->pmhd->presist;
  pr->ClearResistiveEMFs(pr->efld_resist);
  pr->AddEMFGeneralResist(pmbp->pmhd->b0, pr->efld_resist);
  if (pr->wall_cond_reflect || pr->wall_cond_ix1) {
    pmbp->pmhd->WallEmfZero(pr->efld_resist, pr->wall_cond_reflect, pr->wall_cond_ix1);
  }
}

//----------------------------------------------------------------------------------------
//! history: budgets

void OdHistory(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  pdata->nhist = 12;
  pdata->label[0] = "Lz";
  pdata->label[1] = "KE";
  pdata->label[2] = "ME";
  pdata->label[3] = "Eint";
  pdata->label[4] = "Etot";
  pdata->label[5] = "Jdot";       // sum eta |J|^2 dV (code edge EMF)
  pdata->label[6] = "Tvol";       // sum s (J x B)_phi dV
  pdata->label[7] = "Twall_ex";   // wall Maxwell torque, Neumann-extrapolated B_phi
  pdata->label[8] = "Twall_0";    // wall Maxwell torque, first-cell B_phi
  pdata->label[9] = "Treset";     // kinematic: torque removed by the hold
  pdata->label[10] = "dEint";     // sum (e - e_snapshot) dV
  pdata->label[11] = "mass";
  OdResistEMF(pmbp);
  const bool cs = pm->use_cubed_sphere;
  const Real sc = cs ? 0.25*M_PI : 1.0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto &mbp = pmbp->pmb->mb_panel;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  auto &vol = pmbp->pcoord->volume;
  auto &ar1 = pmbp->pcoord->area.x1f;
  auto &u0 = pmbp->pmhd->u0;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto &eta_b = pmbp->pmhd->presist->eta_b;
  auto e1f = pmbp->pmhd->presist->efld_resist.x1e;
  auto e2f = pmbp->pmhd->presist->efld_resist.x2e;
  auto e3f = pmbp->pmhd->presist->efld_resist.x3e;
  auto aux = od->aux;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("od_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &s) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    const int p = cs ? mbp.d_view(m) : 0;
    const Real a2 = sc*CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real a3 = sc*CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    Real q[3], t1[3], t2[3];
    OdGeo(cs, p, a2, a3, q, t1, t2);
    const Real c = aux(m, ACOS, k, j, i), det = 1.0 - c*c, sn = sqrt(det);
    const Real r = x1v(m, i), dv = vol(m, k, j, i);
    // J at the centre from the 4 surrounding edges of each family
    const Real eta1 = 0.25*(eta_b(m, k, j, i) + eta_b(m, k, j-1, i)
                              + eta_b(m, k-1, j, i)
                            + eta_b(m, k-1, j-1, i));   // x1 edges share i
    const Real jr = 0.25*(e1f(m, k, j, i) + e1f(m, k, j+1, i) + e1f(m, k+1, j, i)
                          + e1f(m, k+1, j+1, i))/eta1;
    const Real etl = 0.5*(eta_b(m, k, j, i) + eta_b(m, k, j, i-1));
    const Real etr = 0.5*(eta_b(m, k, j, i) + eta_b(m, k, j, i+1));
    const Real ja = 0.25*((e2f(m, k, j, i) + e2f(m, k+1, j, i))/etl
                        + (e2f(m, k, j, i+1) + e2f(m, k+1, j, i+1))/etr);
    const Real jb = 0.25*((e3f(m, k, j, i) + e3f(m, k, j+1, i))/etl
                        + (e3f(m, k, j, i+1) + e3f(m, k, j+1, i+1))/etr);
    const Real alj = (ja - c*jb)/det, bej = (jb - c*ja)/det;
    const Real eta = eta_b(m, k, j, i);
    const Real j2 = jr*jr + ja*alj + jb*bej;
    // B at the centre
    const Real br = 0.5*(b0.x1f(m, k, j, i) + b0.x1f(m, k, j, i+1));
    const Real bn1 = 0.5*(b0.x2f(m, k, j, i) + b0.x2f(m, k, j+1, i));
    const Real bn2 = 0.5*(b0.x3f(m, k, j, i) + b0.x3f(m, k+1, j, i));
    Real jc[3], bc[3];
    for (int d=0; d<3; ++d) {
      jc[d] = jr*q[d] + alj*t1[d] + bej*t2[d];
      bc[d] = br*q[d] + (bn1/sn)*t1[d] + (bn2/sn)*t2[d];
    }
    const Real f[3] = {jc[1]*bc[2] - jc[2]*bc[1], jc[2]*bc[0] - jc[0]*bc[2],
                       jc[0]*bc[1] - jc[1]*bc[0]};
    const Real tq = r*(q[0]*f[1] - q[1]*f[0]);
    const Real ke = 0.5*(u0(m, IM1, k, j, i)*w0(m, IVX, k, j, i)
                       + u0(m, IM2, k, j, i)*w0(m, IVY, k, j, i)
                       + u0(m, IM3, k, j, i)*w0(m, IVZ, k, j, i));
    const Real lz = aux(m, ALZ2, k, j, i)*u0(m, IM2, k, j, i)
                  + aux(m, ALZ3, k, j, i)*u0(m, IM3, k, j, i);
    // wall Maxwell torque: psi = r (q x B)_z (= r sin B_phi), one face per wall column
    Real twx = 0.0, tw0 = 0.0;
    const BoundaryFlag fi = mb_bcs.d_view(m, BoundaryFace::inner_x1);
    const BoundaryFlag fo = mb_bcs.d_view(m, BoundaryFace::outer_x1);
    const BoundaryFlag rfl = BoundaryFlag::reflect, usr = BoundaryFlag::user;
    const bool iw = (i == is) && (fi == rfl || fi == usr);
    const bool ow = (i == ie) && (fo == rfl || fo == usr);
    if (iw || ow) {
      const int i1 = iw ? i + 1 : i - 1;
      auto psi = [&](const int ii) {
        const Real a1 = 0.5*(b0.x2f(m, k, j, ii) + b0.x2f(m, k, j+1, ii))/sn;
        const Real b1 = 0.5*(b0.x3f(m, k, j, ii) + b0.x3f(m, k+1, j, ii))/sn;
        Real bt[3];
        for (int d=0; d<3; ++d) bt[d] = a1*t1[d] + b1*t2[d];
        return x1v(m, ii)*(q[0]*bt[1] - q[1]*bt[0]);
      };
      const Real p0v = psi(i), p1v = psi(i1);
      const int iface = iw ? i : i + 1;
      const Real rw = iw ? 1.0 : 1.0;   // psi already carries its own radius
      const Real sgn = iw ? -1.0 : 1.0;
      const Real brw = b0.x1f(m, k, j, iface), aw = ar1(m, k, j, iface);
      // psi is r sin B_phi; the lever arm on the wall is r_w sin, and (r B_phi)' = 0
      twx = sgn*rw*(9.0*p0v - p1v)/8.0*brw*aw;
      tw0 = sgn*p0v*brw*aw;
    }
    array_sum::GlobalSum h;
    h.the_array[0] = lz*dv;
    h.the_array[1] = ke*dv;
    h.the_array[2] = (u0(m, IEN, k, j, i) - w0(m, IEN, k, j, i) - ke)*dv;
    h.the_array[3] = w0(m, IEN, k, j, i)*dv;
    h.the_array[4] = u0(m, IEN, k, j, i)*dv;
    h.the_array[5] = eta*j2*dv;
    h.the_array[6] = tq*dv;
    h.the_array[7] = twx;
    h.the_array[8] = tw0;
    h.the_array[9] = aux(m, ATRQ, k, j, i)*dv;
    h.the_array[10] = (w0(m, IEN, k, j, i) - aux(m, AEI0, k, j, i))*dv;
    h.the_array[11] = u0(m, IDN, k, j, i)*dv;
    for (int n=12; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
    s += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  for (int n=0; n<pdata->nhist; ++n) pdata->hdata[n] = sum.the_array[n];
}

//----------------------------------------------------------------------------------------
//! final: errors against the ref1d.py table at t = tlim

void OdFinal(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (std::string(od_ref) == "none") {
    delete od; od = nullptr;
    return;
  }
  // read the table (every rank)
  std::ifstream fin(od_ref);
  if (!fin) {
    std::cout << "### ohm_drag: cannot open ref_file " << od_ref << std::endl;
    delete od; od = nullptr;
    return;
  }
  std::vector<Real> col[5];
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    Real v[5];
    if (!(ss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4])) continue;
    for (int c=0; c<5; ++c) col[c].push_back(v[c]);
  }
  const int nt = static_cast<int>(col[0].size());
  DvceArray2D<Real> tab("od_tab", 5, nt);
  auto htab = Kokkos::create_mirror_view(tab);
  for (int c=0; c<5; ++c) for (int n=0; n<nt; ++n) htab(c, n) = col[c][n];
  Kokkos::deep_copy(tab, htab);
  const Real hr = (col[0][nt-1] - col[0][0])/(nt - 1);

  OdResistEMF(pmbp);
  const bool cs = pm->use_cubed_sphere;
  const Real sc = cs ? 0.25*M_PI : 1.0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto &mbp = pmbp->pmb->mb_panel;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto &vol = pmbp->pcoord->volume;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto &eta_b = pmbp->pmhd->presist->eta_b;
  auto e1f = pmbp->pmhd->presist->efld_resist.x1e;
  auto e2f = pmbp->pmhd->presist->efld_resist.x2e;
  auto e3f = pmbp->pmhd->presist->efld_resist.x3e;
  auto aux = od->aux;
  auto bi = od->binit;
  const Real m_ = od_m, etab = od_etab, ndec = od_ndec, rt = od_rt, w = od_w;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  // slots: 0/1 B_phi faces |err|, |ref| (vol-weighted)
  // 4/5 eta J tangential edges; 6/7 eta J radial edges; 8/9 J tang; 10/11 J radial
  // 12/13 accumulated heating (all cells); 14/15 heating interior (no wall row)
  // 16/17 instantaneous eta J^2; 18/19 torque density; 20 sum dE; 21 sum torque density;
  // 22 sum eta J^2 (code)
  array_sum::GlobalSum sum;
  Kokkos::parallel_reduce("od_err", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &s) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    array_sum::GlobalSum h;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
    const int p = cs ? mbp.d_view(m) : 0;
    const Real x2mn = size.d_view(m).x2min, x2mx = size.d_view(m).x2max;
    const Real x3mn = size.d_view(m).x3min, x3mx = size.d_view(m).x3max;
    const Real a2c = sc*CellCenterX(j-js, nx2, x2mn, x2mx);
    const Real a3c = sc*CellCenterX(k-ks, nx3, x3mn, x3mx);
    const Real a2f = sc*LeftEdgeX(j-js, nx2, x2mn, x2mx);
    const Real a3f = sc*LeftEdgeX(k-ks, nx3, x3mn, x3mx);
    const Real r = x1v(m, i), dv = vol(m, k, j, i);
    Real q[3], t1[3], t2[3], n1[3], n2[3];
    // (a) induced tangential field on the x2 face (j) and x3 face (k) of this cell
    {
      const Real b = OdInterp(tab, nt, hr, 1, r);
      // x2 face
      OdGeo(cs, p, a2f, a3c, q, t1, t2);
      OdNormals(t1, t2, n1, n2);
      Real st = sqrt(fmax(0.0, 1.0 - q[2]*q[2]));
      Real bphi = b*st*q[2];
      Real ph[3] = {-q[1], q[0], 0.0};   // sin * phihat
      Real ref = (st > 1.0e-14) ? bphi*(ph[0]*n1[0] + ph[1]*n1[1])/st : 0.0;
      Real err = fabs(b0.x2f(m, k, j, i) - bi.x2f(m, k, j, i) - ref);
      h.the_array[0] += err*dv;
      h.the_array[1] += fabs(ref)*dv;
      // x3 face
      OdGeo(cs, p, a2c, a3f, q, t1, t2);
      OdNormals(t1, t2, n1, n2);
      st = sqrt(fmax(0.0, 1.0 - q[2]*q[2]));
      bphi = b*st*q[2];
      ph[0] = -q[1]; ph[1] = q[0];
      ref = (st > 1.0e-14) ? bphi*(ph[0]*n2[0] + ph[1]*n2[1])/st : 0.0;
      err = fabs(b0.x3f(m, k, j, i) - bi.x3f(m, k, j, i) - ref);
      h.the_array[0] += err*dv;
      h.the_array[1] += fabs(ref)*dv;
    }
    // (b) edges: x1 edge (r centre, x2 face, x3 face) and, away from the walls, the x2
    // edge (r face, x2 centre, x3 face) and x3 edge (r face, x2 face, x3 centre)
    {
      const Real b = OdInterp(tab, nt, hr, 1, r), up = OdInterp(tab, nt, hr, 2, r);
      OdGeo(cs, p, a2f, a3f, q, t1, t2);
      Real jc[3];
      OdRefJ(r, b, up, q, jc);
      const Real jref = jc[0]*q[0] + jc[1]*q[1] + jc[2]*q[2];
      const Real eta1 = 0.25*(eta_b(m, k, j, i) + eta_b(m, k, j-1, i)
                              + eta_b(m, k-1, j, i)
                              + eta_b(m, k-1, j-1, i));
      const Real etr = OdEta(r, etab, ndec, rt, w);
      const Real e = e1f(m, k, j, i);
      h.the_array[6] += fabs(e - etr*jref);
      h.the_array[7] += fabs(etr*jref);
      h.the_array[10] += fabs(e/eta1 - jref);
      h.the_array[11] += fabs(jref);
      if (i > is) {
        const Real rf = x1f(m, i);
        const Real bf = OdInterp(tab, nt, hr, 1, rf), upf = OdInterp(tab, nt, hr, 2, rf);
        const Real etf = OdEta(rf, etab, ndec, rt, w);
        const Real etc = 0.5*(eta_b(m, k, j, i) + eta_b(m, k, j, i-1));
        // x2 edge: tangent t1
        OdGeo(cs, p, a2c, a3f, q, t1, t2);
        OdRefJ(rf, bf, upf, q, jc);
        Real jt = jc[0]*t1[0] + jc[1]*t1[1] + jc[2]*t1[2];
        Real ec = e2f(m, k, j, i);
        h.the_array[4] += fabs(ec - etf*jt);
        h.the_array[5] += fabs(etf*jt);
        h.the_array[8] += fabs(ec/etc - jt);
        h.the_array[9] += fabs(jt);
        // x3 edge: tangent t2
        OdGeo(cs, p, a2f, a3c, q, t1, t2);
        OdRefJ(rf, bf, upf, q, jc);
        jt = jc[0]*t2[0] + jc[1]*t2[1] + jc[2]*t2[2];
        ec = e3f(m, k, j, i);
        h.the_array[4] += fabs(ec - etf*jt);
        h.the_array[5] += fabs(etf*jt);
        h.the_array[8] += fabs(ec/etc - jt);
        h.the_array[9] += fabs(jt);
      }
    }
    // (c)-(e) cell centre: heating, eta J^2, torque density
    {
      OdGeo(cs, p, a2c, a3c, q, t1, t2);
      const Real ct = q[2], st2 = fmax(0.0, 1.0 - ct*ct);
      const Real b = OdInterp(tab, nt, hr, 1, r), up = OdInterp(tab, nt, hr, 2, r);
      const Real ib = OdInterp(tab, nt, hr, 3, r), id = OdInterp(tab, nt, hr, 4, r);
      const Real etr = OdEta(r, etab, ndec, rt, w);
      const Real p2 = 3.0*ct*ct - 1.0;
      const Real qref = etr/(r*r)*(p2*p2*ib + st2*ct*ct*id);
      const Real de = w0(m, IEN, k, j, i) - aux(m, AEI0, k, j, i);
      h.the_array[12] += fabs(de - qref)*dv;
      h.the_array[13] += fabs(qref)*dv;
      if (i > is && i < is + nx1 - 1) {
        h.the_array[14] += fabs(de - qref)*dv;
        h.the_array[15] += fabs(qref)*dv;
      }
      h.the_array[20] += de*dv;
      // code J at the centre (as in the history)
      const Real c = aux(m, ACOS, k, j, i), det = 1.0 - c*c;
      const Real eta1 = 0.25*(eta_b(m, k, j, i) + eta_b(m, k, j-1, i)
                              + eta_b(m, k-1, j, i)
                              + eta_b(m, k-1, j-1, i));
      const Real jr = 0.25*(e1f(m, k, j, i) + e1f(m, k, j+1, i) + e1f(m, k+1, j, i)
                            + e1f(m, k+1, j+1, i))/eta1;
      const Real etl = 0.5*(eta_b(m, k, j, i) + eta_b(m, k, j, i-1));
      const Real etrr = 0.5*(eta_b(m, k, j, i) + eta_b(m, k, j, i+1));
      const Real ja = 0.25*((e2f(m, k, j, i) + e2f(m, k+1, j, i))/etl
                          + (e2f(m, k, j, i+1) + e2f(m, k+1, j, i+1))/etrr);
      const Real jb = 0.25*((e3f(m, k, j, i) + e3f(m, k, j+1, i))/etl
                          + (e3f(m, k, j, i+1) + e3f(m, k, j+1, i+1))/etrr);
      const Real j2 = jr*jr + (ja*ja + jb*jb - 2.0*c*ja*jb)/det;
      const Real jr0 = b*p2/r, jt0 = -up*sqrt(st2)*ct/r;
      const Real jq = eta_b(m, k, j, i)*j2, jqr = etr*(jr0*jr0 + jt0*jt0);
      h.the_array[16] += fabs(jq - jqr)*dv;
      h.the_array[17] += fabs(jqr)*dv;
      h.the_array[22] += jq*dv;
      const Real tref = m_*st2/(r*r*r)*(b*(2.0*ct*ct - st2) + 2.0*up*ct*ct);
      const Real tq = aux(m, ATRQ, k, j, i);
      h.the_array[18] += fabs(tq - tref)*dv;
      h.the_array[19] += fabs(tref)*dv;
      h.the_array[21] += tq*dv;
    }
    s += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum));
  // radial profiles (meshblock nx1 = mesh nx1 required): sum dE dV, Q_ref dV, |dE-Q| dV,
  // torque density dV, reference torque density dV, |B_phi err| dV, |B_phi ref| dV
  const bool rprof = pin->GetOrAddBoolean("problem", "radial_profile", false) &&
                     (nx1 == pm->mesh_indcs.nx1);
  constexpr int NPR = 7;
  DvceArray2D<Real> prof("od_prof", nx1, NPR);
  if (rprof) {
    par_for("od_prof", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ks+nx3-1,
            js, js+nx2-1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const int p = cs ? mbp.d_view(m) : 0;
      const Real a2c = sc*CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                      size.d_view(m).x2max);
      const Real a3c = sc*CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                      size.d_view(m).x3max);
      Real q[3], t1[3], t2[3];
      OdGeo(cs, p, a2c, a3c, q, t1, t2);
      const Real r = x1v(m, i), dv = vol(m, k, j, i);
      const Real ct = q[2], st2 = fmax(0.0, 1.0 - ct*ct);
      const Real b = OdInterp(tab, nt, hr, 1, r), up = OdInterp(tab, nt, hr, 2, r);
      const Real ib = OdInterp(tab, nt, hr, 3, r), id = OdInterp(tab, nt, hr, 4, r);
      const Real etr = OdEta(r, etab, ndec, rt, w);
      const Real p2 = 3.0*ct*ct - 1.0;
      const Real qref = etr/(r*r)*(p2*p2*ib + st2*ct*ct*id);
      const Real de = w0(m, IEN, k, j, i) - aux(m, AEI0, k, j, i);
      const Real tref = m_*st2/(r*r*r)*(b*(2.0*ct*ct - st2) + 2.0*up*ct*ct);
      Real n1[3], n2[3];
      OdNormals(t1, t2, n1, n2);
      // cell-centre B_phi from the two x3 faces (projection on nhat_eta)
      const Real bn = 0.5*(b0.x3f(m, k, j, i) - bi.x3f(m, k, j, i)
                         + b0.x3f(m, k+1, j, i) - bi.x3f(m, k+1, j, i));
      const Real st = sqrt(st2);
      const Real ref = (st > 1.0e-14) ? b*ct*(-q[1]*n2[0] + q[0]*n2[1]) : 0.0;
      const int ir = i - is;
      Kokkos::atomic_add(&prof(ir, 0), de*dv);
      Kokkos::atomic_add(&prof(ir, 1), qref*dv);
      Kokkos::atomic_add(&prof(ir, 2), fabs(de - qref)*dv);
      Kokkos::atomic_add(&prof(ir, 3), aux(m, ATRQ, k, j, i)*dv);
      Kokkos::atomic_add(&prof(ir, 4), tref*dv);
      Kokkos::atomic_add(&prof(ir, 5), fabs(bn - ref)*dv);
      Kokkos::atomic_add(&prof(ir, 6), fabs(ref)*dv);
    });
    auto hp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), prof);
    std::vector<Real> pv(nx1*NPR);
    for (int i=0; i<nx1; ++i) for (int n=0; n<NPR; ++n) pv[i*NPR + n] = hp(i, n);
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, pv.data(), nx1*NPR, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::string pf = std::string(od_err) + ".prof";
      FILE *fp = std::fopen(pf.c_str(), "a");
      if (fp != nullptr) {
        std::fprintf(fp, "# %s nx1 %d nx2 %d t %.6e: i dE Qref |dE-Q| trq trqref "
                     "|dBphi| |Bphi|\n", cs ? "cs" : "sp", pm->mesh_indcs.nx1,
                     pm->mesh_indcs.nx2, pm->time);
        for (int i=0; i<nx1; ++i) {
          std::fprintf(fp, "%3d", i);
          for (int n=0; n<NPR; ++n) std::fprintf(fp, " %13.5e", pv[i*NPR + n]);
          std::fprintf(fp, "\n");
        }
        std::fclose(fp);
      }
    }
  }
  Real v[NREDUCTION_VARIABLES];
  for (int n=0; n<NREDUCTION_VARIABLES; ++n) v[n] = sum.the_array[n];
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, v, NREDUCTION_VARIABLES, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    // reference totals from the table header
    Real qtot = 0.0, trq = 0.0;
    std::ifstream f2(od_ref);
    std::string l2;
    while (std::getline(f2, l2)) {
      const auto pq = l2.find("Q_total");
      if (pq != std::string::npos) {
        std::istringstream ss(l2.substr(pq + 7));
        std::string tk;
        ss >> qtot >> tk >> trq;
        break;
      }
    }
    char buf[2048];
    std::snprintf(buf, sizeof(buf),
        "%s nx1 %d nx2 %d t %.6e eta_bot %g | Bphi_L1 %.4e | etaJ_tan %.4e etaJ_rad %.4e"
        " | J_tan %.4e J_rad %.4e | heat_L1 %.4e heat_L1_int %.4e heat_tot %.6e"
        " heat_ref %.6e ratio %.6f | jq_L1 %.4e jq_tot %.6e | trq_L1 %.4e trq_tot %.6e"
        " trq_ref %.6e ratio %.6f",
        cs ? "cs" : "sp", pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->time, od_etab,
        v[0]/v[1], v[4]/v[5], v[6]/v[7], v[8]/v[9], v[10]/v[11], v[12]/v[13],
        v[14]/v[15], v[20], qtot, v[20]/qtot, v[16]/v[17], v[22], v[18]/v[19], v[21],
        trq, v[21]/trq);
    std::cout << "OHMDRAG " << buf << std::endl;
    FILE *fo = std::fopen(od_err, "a");
    if (fo != nullptr) {
      std::fprintf(fo, "%s\n", buf);
      std::fclose(fo);
    }
  }
  delete od;
  od = nullptr;
}

}  // namespace
