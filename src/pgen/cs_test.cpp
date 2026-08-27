//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_test.cpp
//! \brief Clean validation problems for the gnomonic-equiangular CUBED SPHERE grid.
//!
//! The pre-existing cubed_sphere.cpp / solid_body_rot.cpp pgens are scratch files: they
//! do not compile against the current well-balancing API, they write the KINETIC energy
//! into the primitive energy slot (w0(IEN) is the INTERNAL energy density p/(gamma-1),
//! not the total), and their velocity transforms are only valid on one panel. Diagnosing
//! the grid through them means debugging two things at once, so this file provides a
//! deliberately minimal instrument instead.
//!
//! problem/iprob selects:
//!   1 = STATIC UNIFORM STATE.  rho = d0, p = p0, v = 0 everywhere on all six panels.
//!       The exact solution is "unchanged, forever". This isolates the GEOMETRY and the
//!       inter-panel communication from anything to do with the velocity representation:
//!       with v = 0 every velocity rotation in coordinates.hpp is identically zero, so a
//!       failure here can only come from the areas/volumes, the geometric pressure source
//!       terms, or the panel boundary exchange. Run it before anything else.
//!
//!   3 = RIGID-ROTATION EQUILIBRIUM.  The sharpest available test of the VELOCITY
//!       representation. A barotropic fluid in rigid rotation v = Omega zhat x r is an
//!       EXACT steady solution of the compressible Euler equations with no gravity,
//!       provided the pressure supplies the centripetal force:
//!           (v.grad)v = -Omega^2 R Rhat   =>   grad p = rho Omega^2 R Rhat
//!           =>  p = p0 + 0.5 rho Omega^2 R^2,   R = cylindrical radius.
//!       Note v.rhat = 0 everywhere, so nothing crosses the radial boundaries; the exact
//!       solution is imposed in the radial ghosts by the user BC below.
//!       Unlike iprob=1 this has a nonzero velocity with nontrivial components in every
//!       panel basis, so it exercises the contravariant/covariant handling and the
//!       velocity transforms across panel seams. Any drift is a bug, and WHERE it appears
//!       is diagnostic: the non-orthogonality angle has cos(theta) = -xy/(CD), which
//!       vanishes at each panel CENTRE and is largest at the panel CORNERS.
//!
//!   2 = SOLID-BODY ROTATION.  Uniform rho and p, plus a rigid rotation about the
//!       Cartesian z axis, and a Gaussian passive scalar blob to advect. The velocity
//!       field is a Killing vector of the sphere, so the flow is a steady state and the
//!       blob should return to its starting point undistorted after one period. This is
//!       the test that exercises the contravariant velocity components and their
//!       transformation across panel seams.
//!
//! Both are run in 2D (nx3 = 1), where Coordinates::CoordGnomonicEquiangle sets the
//! radius to 1 and the mesh is the unit sphere.

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "pgen.hpp"

namespace {

// Parameters shared between the initial condition and the radial boundary condition.
// ProblemGenerator carries no ParameterInput pointer, so stash them at file scope when
// UserProblem reads them.
Real cs_d0 = 1.0, cs_p0 = 1.0, cs_omega = 1.0;
int  cs_iprob = 1;
Real cs_amp = 0.5, cs_r0 = 1.0;
int  cs_exact_panel_ghosts = 0;

//----------------------------------------------------------------------------------------
//! \brief iprob = 5: density a function of the SPHERICAL RADIUS alone, uniform pressure,
//! zero velocity. With grad p = 0 and v = 0 this is an exact steady solution of Euler for
//! any rho(r), so the active cells must not move at all.
//!
//! Its purpose is to isolate the RADIAL index handling in the panel halo exchange from
//! everything angular. A function of r alone is invariant under any angular remapping, so
//! the ~21%-of-a-cell angular misalignment between a panel's gnomonic ghost extension and
//! its neighbour's cells is INVISIBLE here. What is not invisible is a reversal or a
//! swap of the x3 (radial) index across a seam: that flips the radial profile in the
//! ghosts and shows up as an O(1) error, not an O(dx) one.

KOKKOS_INLINE_FUNCTION
Real RadialProfile(const Real r, const Real d0, const Real amp, const Real r0) {
  return d0*(1.0 + amp*(r - r0));
}

//----------------------------------------------------------------------------------------
//! \brief iprob = 6 density. A smooth function of the CARTESIAN position with no
//! symmetry a panel map could hide behind: every direction and both parities differ.

//! The angular part uses the DIRECTION cosines, so it is bounded by |sum of coeffs| and
//! the whole profile stays positive for amp <= 0.3 -- a negative rho here is silently
//! floored and then shows up as a fake O(1) "error" in the active cells.

KOKKOS_INLINE_FUNCTION
Real CartProfile(const Real cx, const Real cy, const Real cz, const Real r,
                 const Real d0, const Real r0, const Real amp) {
  return d0*(1.0 + amp*(0.31*cx + 0.57*cy + 0.79*cz + 0.23*cx*cy - 0.17*cy*cz)
                 + 0.3*(r - r0));
}

//----------------------------------------------------------------------------------------
//! \brief Thin wrappers over coordinates/cubed_sphere.hpp. The panel frames used to be
//! written out a second time here, and that copy had panels 3 and 4 interchanged -- a
//! labelling that is self-consistent within any single panel and therefore invisible to
//! every single-panel test, while making a globally defined initial condition (iprob = 3)
//! a different field on two of the six panels. There is now exactly one copy of the
//! frames; do not reintroduce another.

KOKKOS_INLINE_FUNCTION
void PanelToCart(const int p, const Real xi, const Real eta,
                 Real &cx, Real &cy, Real &cz) {
  Real q[3];
  cubed_sphere::PanelToCart(p, xi, eta, q);
  cx = q[0];
  cy = q[1];
  cz = q[2];
  return;
}

//----------------------------------------------------------------------------------------
//! \brief Project a Cartesian vector onto the panel's UNIT tangent basis, returning the
//! contravariant components (a,b) such that V = a*e1hat + b*e2hat. The basis is not
//! orthogonal, so this inverts the Gram matrix [[1,cs],[cs,1]]. `cs` is returned for
//! cross-checking against Coordinates::cos_cell.

KOKKOS_INLINE_FUNCTION
void CartToPanelVec(const int p, const Real xi, const Real eta,
                    const Real vx, const Real vy, const Real vz,
                    Real &a, Real &b, Real &cs) {
  Real e1[3], e2[3];
  cubed_sphere::PanelTangents(p, xi, eta, e1, e2);
  cs = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
  const Real d1 = vx*e1[0] + vy*e1[1] + vz*e1[2];
  const Real d2 = vx*e2[0] + vy*e2[1] + vz*e2[2];
  const Real det = 1.0 - cs*cs;
  a = (d1 - cs*d2)/det;
  b = (d2 - cs*d1)/det;
  return;
}

//----------------------------------------------------------------------------------------
//! \brief The exact rigid-rotation equilibrium (iprob = 3) at one point.
//! r is the spherical radius, (xi,eta) the equiangular coordinates on panel p.

KOKKOS_INLINE_FUNCTION
void RigidRotState(const int p, const Real xi, const Real eta, const Real r,
                   const Real d0, const Real p0, const Real omega, const Real gm1,
                   Real &dn, Real &ie, Real &v1, Real &v2, Real &v3) {
  Real cx, cy, cz;
  PanelToCart(p, xi, eta, cx, cy, cz);
  // Cartesian position and the rigid-rotation velocity there.
  const Real px = r*cx, py = r*cy;
  const Real vx = -omega*py;
  const Real vy =  omega*px;
  const Real R2 = px*px + py*py;      // cylindrical radius squared
  dn = d0;
  ie = (p0 + 0.5*d0*omega*omega*R2)/gm1;
  Real a, b, cs;
  CartToPanelVec(p, xi, eta, vx, vy, 0.0, a, b, cs);
  v1 = 0.0;                            // v.rhat = 0 for rotation about the z axis
  v2 = a;                              // xi   (x2)
  v3 = b;                              // eta  (x3)
}

} // namespace

void CSTestRadialBC(Mesh *pm);
void CSTestGhostCheck(ParameterInput *pin, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem
//! \brief Sets the initial conditions for the cubed-sphere validation problems.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "cs_test requires a <hydro> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!pmy_mesh_->use_cubed_sphere) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "cs_test requires mesh/use_cubed_sphere = true" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const int iprob = pin->GetInteger("problem", "iprob");
  const Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  const Real omega = pin->GetOrAddReal("problem", "omega", 1.0);
  cs_d0 = d0; cs_p0 = p0; cs_omega = omega;
  cs_iprob = iprob;
  cs_amp = pin->GetOrAddReal("problem", "amp", 0.5);
  cs_exact_panel_ghosts = pin->GetOrAddInteger("problem",
                                              "exact_panel_ghosts", 0);
  cs_r0 = pmy_mesh_->mesh_size.x1min;
  if (iprob >= 3 && iprob <= 7 && iprob != 4) user_bcs_func = CSTestRadialBC;
  if (iprob >= 3 && iprob <= 7 && iprob != 4) pgen_final_func = CSTestGhostCheck;
  const Real blob_w = pin->GetOrAddReal("problem", "blob_width", 0.3);
  const Real amp_ = cs_amp, r0_ = cs_r0;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &js = indcs.js; int &ks = indcs.ks;
  int &ie = indcs.ie; int &je = indcs.je; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;

  auto &w0 = pmbp->phydro->w0;
  const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  const int nhyd = pmbp->phydro->nhydro;
  const int nscal = pmbp->phydro->nscalars;

  // The equiangular coordinates are xi = (pi/4) x1, eta = (pi/4) x2, matching
  // Coordinates::CoordGnomonicEquiangle, so x1 and x2 both run over [-1,1].
  par_for("pgen_cs_test", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // x1 is RADIAL; xi = x2 and eta = x3, matching CoordGnomonicEquiangle.
    const Real x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                  size.d_view(m).x2max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                  size.d_view(m).x3max);
    const Real xi  = 0.25*M_PI*x2v;
    const Real eta = 0.25*M_PI*x3v;
    const Real rad = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                  size.d_view(m).x1max);

    if (iprob == 7) {
      // TAG FIELD. v = 0 and uniform p make ANY rho an exact steady state, so rho can
      // carry a label instead of a physical profile: 1000*panel + (j-js) + (k-ks)/1000.
      // Decoding a ghost cell then names the exact (panel, j, k) the halo fetched it
      // from -- the only diagnostic that distinguishes a wrong neighbour panel from a
      // wrong orientation from a wrong offset.
      w0(m,IDN,k,j,i) = 1000.0*mbpanel.d_view(m) + (j-js) + 0.001*(k-ks);
      w0(m,IEN,k,j,i) = p0/gm1;
      w0(m,IVX,k,j,i) = 0.0;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      return;
    }

    if (iprob == 6) {
      Real cx, cy, cz;
      PanelToCart(mbpanel.d_view(m), xi, eta, cx, cy, cz);
      w0(m,IDN,k,j,i) = CartProfile(cx, cy, cz, rad, d0, r0_, amp_);
      w0(m,IEN,k,j,i) = p0/gm1;
      w0(m,IVX,k,j,i) = 0.0;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      return;
    }

    if (iprob == 5) {
      w0(m,IDN,k,j,i) = RadialProfile(rad, d0, amp_, r0_);
      w0(m,IEN,k,j,i) = p0/gm1;
      w0(m,IVX,k,j,i) = 0.0;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      return;
    }

    if (iprob == 3) {
      Real dn, ie_, v1, v2, v3;
      RigidRotState(mbpanel.d_view(m), xi, eta, rad, d0, p0, omega, gm1,
                    dn, ie_, v1, v2, v3);
      w0(m,IDN,k,j,i) = dn;
      w0(m,IEN,k,j,i) = ie_;
      w0(m,IVX,k,j,i) = v1;
      w0(m,IVY,k,j,i) = v2;
      w0(m,IVZ,k,j,i) = v3;
      if (nscal > 0) w0(m,nhyd,k,j,i) = 1.0;
      return;
    }

    w0(m,IDN,k,j,i) = d0;
    // w0(IEN) is the INTERNAL energy density, i.e. p/(gamma-1). See ideal_hyd.cpp.
    w0(m,IEN,k,j,i) = p0/gm1;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;

    if (iprob == 2) {
      // Rigid rotation about the Cartesian z axis: v_cart = omega * (zhat x r).
      // Decompose onto the panel's own tangent basis. The gnomonic coordinate basis
      // vectors are e_xi = d(rhat)/d(xi) and e_eta = d(rhat)/d(eta); the code stores
      // CONTRAVARIANT components on the UNIT-normalised versions of those two vectors
      // (see GnomonicEquianglePrimFaceX1 in coordinates.hpp, which multiplies by
      // sin_cell/cos_cell to reach a locally orthonormal frame).
      const Real x = std::tan(xi);
      const Real y = std::tan(eta);
      const Real delta = std::sqrt(1.0 + x*x + y*y);
      const Real C = std::sqrt(1.0 + x*x);
      const Real D = std::sqrt(1.0 + y*y);

      const int p = mbpanel.d_view(m);
      Real cx, cy, cz;
      PanelToCart(p, xi, eta, cx, cy, cz);
      // Cartesian velocity of the rigid rotation at this point (unit sphere).
      const Real vx = -omega*cy;
      const Real vy =  omega*cx;
      const Real vz =  0.0;

      // Tangent basis: differentiate the panel map wrt xi and eta. Rather than hard-code
      // six sets of derivatives, difference the map -- this is an initial condition, so a
      // centred difference at machine-epsilon-safe width is accurate enough and cannot
      // silently disagree with PanelToCart above.
      const Real h = 1.0e-6;
      Real xp, yp, zp, xm, ym, zm;
      PanelToCart(p, xi+h, eta, xp, yp, zp);
      PanelToCart(p, xi-h, eta, xm, ym, zm);
      Real e1x = (xp-xm)/(2.0*h), e1y = (yp-ym)/(2.0*h), e1z = (zp-zm)/(2.0*h);
      PanelToCart(p, xi, eta+h, xp, yp, zp);
      PanelToCart(p, xi, eta-h, xm, ym, zm);
      Real e2x = (xp-xm)/(2.0*h), e2y = (yp-ym)/(2.0*h), e2z = (zp-zm)/(2.0*h);
      // Normalise to unit tangent vectors.
      const Real n1 = std::sqrt(e1x*e1x + e1y*e1y + e1z*e1z);
      const Real n2 = std::sqrt(e2x*e2x + e2y*e2y + e2z*e2z);
      e1x /= n1; e1y /= n1; e1z /= n1;
      e2x /= n2; e2y /= n2; e2z /= n2;

      // Solve v = a e1 + b e2 for (a,b) with the non-orthogonal Gram matrix
      // [[1, cs],[cs, 1]], cs = e1.e2. These a,b are the contravariant components on the
      // unit-normalised basis, which is what the solver expects in IVX/IVY.
      const Real cs = e1x*e2x + e1y*e2y + e1z*e2z;
      const Real d1 = vx*e1x + vy*e1y + vz*e1z;
      const Real d2 = vx*e2x + vy*e2y + vz*e2z;
      const Real det = 1.0 - cs*cs;
      w0(m,IVX,k,j,i) = (d1 - cs*d2)/det;
      w0(m,IVY,k,j,i) = (d2 - cs*d1)/det;
      w0(m,IVZ,k,j,i) = 0.0;

      // A cross-check that the stored geometry agrees with this file's panel map: cs here
      // must equal cos_cell(m,j,i) = -xy/(CD) from CoordGnomonicEquiangle.
      (void)delta; (void)C; (void)D;
    }

    // Passive scalar: a Gaussian blob centred on the +x axis, used to see advection.
    if (nscal > 0) {
      Real cx, cy, cz;
      PanelToCart(mbpanel.d_view(m), xi, eta, cx, cy, cz);
      // angular distance from (1,0,0)
      const Real cosd = fmin(1.0, fmax(-1.0, cx));
      const Real dist = std::acos(cosd);
      w0(m,nhyd,k,j,i) = std::exp(-SQR(dist/blob_w));
    }
  });

  // Convert primitives to conserved.
  auto &u0 = pmbp->phydro->u0;
  // NOT peos->PrimToCons: that assumes an orthonormal basis. On the cubed sphere the
  // conserved momentum is the COVARIANT one. See Coordinates::GnomonicEquiangleLowerMom.
  pmbp->pcoord->GnomonicEquiangleLowerMom(w0, u0, is, ie, js, je, ks, ke);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn CSTestRadialBC
//! \brief Impose the exact iprob=3 solution in the RADIAL (x3) ghost zones.
//!
//! The lateral (x1/x2) ghosts are filled by the panel exchange; only the two radial
//! boundaries are physical here. Because v.rhat = 0 in the exact solution nothing flows
//! through them, so this BC is not doing any work beyond holding the exact state -- which
//! is what makes any drift in the interior attributable to the scheme rather than the BC.

void CSTestRadialBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &ks = indcs.ks;
  int &ng = indcs.ng;
  int n2 = indcs.nx2 + 2*ng;
  int n3 = indcs.nx3 + 2*ng;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &w0 = pmbp->phydro->w0;
  auto &u0 = pmbp->phydro->u0;
  const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  const int nhyd = pmbp->phydro->nhydro;
  const int nscal = pmbp->phydro->nscalars;

  const Real d0 = cs_d0, p0 = cs_p0, omega = cs_omega;
  const int iprob = cs_iprob;
  const Real amp_ = cs_amp, r0_ = cs_r0;

  // The RADIAL direction is x1. The lateral (x2/x3) ghosts are filled by the panel
  // exchange; only the two radial boundaries are physical here.
  par_for("cs_test_rbc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
          0,(n2-1), 0,(ng-1),
  KOKKOS_LAMBDA(int m, int k, int j, int ig) {
    const Real x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                  size.d_view(m).x2max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                  size.d_view(m).x3max);
    const Real xi  = 0.25*M_PI*x2v;
    const Real eta = 0.25*M_PI*x3v;
    const int p = mbpanel.d_view(m);

    for (int side=0; side<2; ++side) {
      if (side == 0 &&
          mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
      if (side == 1 &&
          mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
      const int i = (side == 0) ? (is-ng+ig) : (ie+1+ig);
      const Real rad = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                    size.d_view(m).x1max);
      Real dn, ie_, v1, v2, v3;
      if (iprob == 7) {
        dn = 1000.0*p + (j-js) + 0.001*(k-ks);
        ie_ = p0/gm1;
        v1 = 0.0; v2 = 0.0; v3 = 0.0;
      } else if (iprob == 6) {
        Real cx, cy, cz;
        PanelToCart(p, xi, eta, cx, cy, cz);
        dn = CartProfile(cx, cy, cz, rad, d0, r0_, amp_);
        ie_ = p0/gm1;
        v1 = 0.0; v2 = 0.0; v3 = 0.0;
      } else if (iprob == 5) {
        dn = RadialProfile(rad, d0, amp_, r0_);
        ie_ = p0/gm1;
        v1 = 0.0; v2 = 0.0; v3 = 0.0;
      } else {
        RigidRotState(p, xi, eta, rad, d0, p0, omega, gm1, dn, ie_, v1, v2, v3);
      }
      w0(m,IDN,k,j,i) = dn;
      w0(m,IEN,k,j,i) = ie_;
      w0(m,IVX,k,j,i) = v1;
      w0(m,IVY,k,j,i) = v2;
      w0(m,IVZ,k,j,i) = v3;
      // Covariant momentum, matching Coordinates::GnomonicEquiangleLowerMom: the metric
      // acts on the ANGULAR pair (IM2,IM3) only. ccos is the cosine of the angle between
      // the two tangent basis vectors at this point.
      Real ca, cb, ccos;
      CartToPanelVec(p, xi, eta, 1.0, 0.0, 0.0, ca, cb, ccos);
      u0(m,IDN,k,j,i) = dn;
      u0(m,IM1,k,j,i) = dn*v1;
      u0(m,IM2,k,j,i) = dn*(v2 + ccos*v3);
      u0(m,IM3,k,j,i) = dn*(v3 + ccos*v2);
      u0(m,IEN,k,j,i) = ie_
                      + 0.5*dn*(v1*v1 + v2*v2 + v3*v3 + 2.0*ccos*v2*v3);
      for (int n=nhyd; n<nhyd+nscal; ++n) {
        w0(m,n,k,j,i) = 1.0;
        u0(m,n,k,j,i) = dn;
      }
    }
  });

  // ---------------------------------------------------------------------------------
  // DIAGNOSTIC MODE problem/exact_panel_ghosts. Overwrite the x2/x3 (panel) ghost zones
  // with the exact solution, evaluated on THIS panel's own extended gnomonic map and
  // projected onto THIS panel's own tangent basis. That is what a perfect seam halo
  // would produce, so it isolates the halo exchange from the rest of the scheme: any
  // residual left over is NOT the seam.
  const int pgmode = cs_exact_panel_ghosts;
  if (pgmode > 0 && iprob == 3) {
    int n1 = indcs.nx1 + 2*ng;
    int &je = indcs.je; int &ke = indcs.ke;
    par_for("cs_test_pghost", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
            0,(n2-1), 0,(n1-1),
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const bool joff = (j < js) || (j > je);
      const bool koff = (k < ks) || (k > ke);
      if (!joff && !koff) return;
      // mode 2: FACE ghosts only, leave the corner blocks to the exchange
      if (pgmode == 2 && joff && koff) return;
      const Real x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                    size.d_view(m).x2max);
      const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                    size.d_view(m).x3max);
      const Real xi  = 0.25*M_PI*x2v;
      const Real eta = 0.25*M_PI*x3v;
      const int p = mbpanel.d_view(m);
      const Real rad = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                    size.d_view(m).x1max);
      Real dn, ie_, v1, v2, v3;
      RigidRotState(p, xi, eta, rad, d0, p0, omega, gm1, dn, ie_, v1, v2, v3);
      Real ca, cb, ccos;
      CartToPanelVec(p, xi, eta, 1.0, 0.0, 0.0, ca, cb, ccos);
      w0(m,IDN,k,j,i) = dn;
      w0(m,IEN,k,j,i) = ie_;
      w0(m,IVX,k,j,i) = v1;
      w0(m,IVY,k,j,i) = v2;
      w0(m,IVZ,k,j,i) = v3;
      u0(m,IDN,k,j,i) = dn;
      u0(m,IM1,k,j,i) = dn*v1;
      u0(m,IM2,k,j,i) = dn*(v2 + ccos*v3);
      u0(m,IM3,k,j,i) = dn*(v3 + ccos*v2);
      u0(m,IEN,k,j,i) = ie_
                      + 0.5*dn*(v1*v1 + v2*v2 + v3*v3 + 2.0*ccos*v2*v3);
      for (int n=nhyd; n<nhyd+nscal; ++n) {
        w0(m,n,k,j,i) = 1.0;
        u0(m,n,k,j,i) = dn;
      }
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn CSTestGhostCheck
//! \brief Diagnostic for iprob = 5, called from Driver::Finalize().
//!
//! Reports, per (panel, face), the max error in the density relative to the exact rho(r):
//!   * over the ACTIVE cells -- this must be zero to round-off, since rho(r) with uniform
//!     pressure and v = 0 is an exact steady solution;
//!   * over the x1 and x2 FACE ghost zones, which are filled entirely by the panel halo
//!     exchange. A function of r alone cannot be disturbed by any ANGULAR remapping, so a
//!     nonzero value here is direct evidence that the exchange is mishandling the RADIAL
//!     index. Corner ghosts are excluded (j,k held in the active range) because they are
//!     filled from edge neighbours and are not the thing under test.

void CSTestGhostCheck(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;

  auto &w0 = pmbp->phydro->w0;
  auto w0_h = Kokkos::create_mirror_view(w0);
  Kokkos::deep_copy(w0_h, w0);
  auto &size = pmbp->pmb->mb_size;
  size.template sync<HostMemSpace>();
  auto &mbpanel = pmbp->pmb->mb_panel;
  mbpanel.template sync<HostMemSpace>();

  const Real d0 = cs_d0, amp = cs_amp, r0 = cs_r0;

  // --- iprob = 3: LOCALISE the spurious radial velocity ------------------------------
  // The exact rigid-rotation solution has v1 = 0 identically (x1 is radial), so v1 is a
  // pure error field. Where it lives says which mechanism produces it: an error from the
  // panel halo (whose ghost cells sit ~21% of a cell from the neighbour cells that fill
  // them) is confined to a band a few cells wide along each seam, whereas an error in the
  // geometric source terms is a bulk effect. Measure this after ONE step: by the end of a
  // run the seam error has propagated and filled the domain.
  // How wrong is the panel halo itself, component by component? Compare the FACE ghosts
  // against the exact solution evaluated on THIS panel's own extended gnomonic map and
  // projected onto THIS panel's own tangent basis -- i.e. against what a perfect seam
  // halo would have delivered. rho isolates the ALONG-SEAM OFFSET (a scalar cannot see a
  // basis error); IVY/IVZ additionally see the tangent-basis transform across the seam.
  // O(dx) here (halving per refinement) means offset only; O(1) means the vector
  // transform is wrong too.
  if (cs_iprob == 3) {
    std::cout << "### CS PANEL-HALO ERROR (iprob=3, face ghosts vs exact)" << std::endl;
    std::cout << "  panel     d(rho)        d(v_r)        d(v_xi)       d(v_eta)"
              << std::endl;
    const Real p0_ = cs_p0, om_ = cs_omega;
    const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      const int p = mbpanel.h_view(m);
      Real e[4] = {0.0, 0.0, 0.0, 0.0};
      for (int i=is; i<=ie; ++i) {
        const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                      size.h_view(m).x1max);
        auto probe = [&](int j, int k) {
          const Real xi  = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                                  size.h_view(m).x2max);
          const Real eta = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                                  size.h_view(m).x3max);
          Real dn, ie_, v1, v2, v3;
          RigidRotState(p, xi, eta, rad, d0, p0_, om_, gm1, dn, ie_, v1, v2, v3);
          e[0] = fmax(e[0], fabs(w0_h(m,IDN,k,j,i) - dn));
          e[1] = fmax(e[1], fabs(w0_h(m,IVX,k,j,i) - v1));
          e[2] = fmax(e[2], fabs(w0_h(m,IVY,k,j,i) - v2));
          e[3] = fmax(e[3], fabs(w0_h(m,IVZ,k,j,i) - v3));
        };
        for (int k=ks; k<=ke; ++k) {
          for (int g=0; g<ng; ++g) { probe(js-1-g,k); probe(je+1+g,k); }
        }
        for (int j=js; j<=je; ++j) {
          for (int g=0; g<ng; ++g) { probe(j,ks-1-g); probe(j,ke+1+g); }
        }
      }
      std::printf("  %4d   %12.5e  %12.5e  %12.5e  %12.5e\n", p, e[0],e[1],e[2],e[3]);
    }
  }

  if (cs_iprob == 3) {
    const int nseam = ng;
    std::cout << "### CS v_r LOCALISATION (iprob=3, exact v_r = 0)" << std::endl;
    std::cout << "  panel   max|v_r|     at(i,j,k)      sum v_r^2 seam  interior"
              << "     ncell seam/int" << std::endl;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      Real vmax = 0.0; int im=0, jm=0, km=0;
      Real s_seam = 0.0, s_int = 0.0;
      int n_seam = 0, n_int = 0;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          // the panel-tangential axes are x2 and x3, so a seam band is defined by j,k
          const bool near = (j-js < nseam) || (je-j < nseam) ||
                            (k-ks < nseam) || (ke-k < nseam);
          for (int i=is; i<=ie; ++i) {
            const Real vr = w0_h(m,IVX,k,j,i);
            if (fabs(vr) > vmax) { vmax = fabs(vr); im=i-is; jm=j-js; km=k-ks; }
            if (near) { s_seam += vr*vr; ++n_seam; } else { s_int += vr*vr; ++n_int; }
          }
        }
      }
      std::printf("  %4d   %10.3e   %3d %3d %3d   %12.5e   %12.5e   %6d %6d\n",
                  mbpanel.h_view(m), vmax, im, jm, km, s_seam, s_int, n_seam, n_int);
    }
    return;
  }

  // --- iprob = 7: WHERE did each ghost cell come from? -------------------------------
  // Decode the tag rho = 1000*panel + j + k/1000 out of the face ghosts and print, for
  // each (panel, face), the source panel and the source (j,k) of the two extreme
  // along-seam cells and of the seam midline. Compare against the expected source: the
  // seam is exact at the shared edge, so ghost layer g of panel P must come from active
  // layer g of the neighbour, at the along-seam index given by the panel table.
  if (cs_iprob == 7) {
    const int n2a = indcs.nx2, n3a = indcs.nx3;
    std::cout << "### CS GHOST SOURCE MAP (iprob=7)" << std::endl;
    // First: did anything MOVE? v = 0 with uniform p is exactly steady for any rho, so
    // every active cell must still carry its own tag to round-off.
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      Real e = 0.0; int ei=-1, ej=-1, ek=-1;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            const Real ex = 1000.0*mbpanel.h_view(m) + (j-js) + 0.001*(k-ks);
            const Real d = fabs(w0_h(m,IDN,k,j,i) - ex);
            if (d > e) { e = d; ei=i-is; ej=j-js; ek=k-ks; }
          }
        }
      }
      std::printf("  panel %d  max active drift %12.5e at (%d,%d,%d)\n",
                  mbpanel.h_view(m), e, ei, ej, ek);
    }
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      const int p = mbpanel.h_view(m);
      const int i = is;
      for (int f=0; f<4; ++f) {
        std::printf("  panel %d  face %s :", p,
                    (f==0?"-x2":(f==1?"+x2":(f==2?"-x3":"+x3"))));
        for (int g=0; g<ng; ++g) {
          // three samples along the seam: first active index, middle, last
          for (int t=0; t<3; ++t) {
            const int a = (t==0) ? 0 : ((t==1) ? (f<2?n3a/2:n2a/2)
                                              : ((f<2?n3a:n2a)-1));
            int j = js+a, k = ks+a;
            if (f == 0) {
              j = js-1-g;
            } else if (f == 1) {
              j = je+1+g;
            } else if (f == 2) {
              k = ks-1-g;
            } else {
              k = ke+1+g;
            }
            const Real tag = w0_h(m,IDN,k,j,i);
            const int sp = static_cast<int>(std::floor(tag/1000.0));
            const Real rem = tag - 1000.0*sp;
            const int sj = static_cast<int>(std::floor(rem + 1.0e-6));
            const int sk = static_cast<int>(std::lround((rem - sj)*1000.0));
            std::printf("  g%d[%2d]=P%d(%2d,%2d)", g, a, sp, sj, sk);
          }
        }
        std::printf("\n");
      }
    }
    return;
  }

  // --- iprob = 6: what does the panel halo do to an ANGULARLY varying scalar? --------
  // rho = f(cartesian position) with v = 0 and uniform p is an exact steady state, so the
  // active cells are a null test, and every ghost cell has a KNOWN exact value: f
  // evaluated at the physical point the ghost cell occupies under this panel's own
  // gnomonic map, extended past the seam. Compare per FACE (corner ghosts excluded --
  // they are measured to be irrelevant to the flux stencil at ng = 2).
  //   * O(1) error, flat in resolution -> the exchange is reaching the wrong panel or the
  //     wrong along-seam orientation.
  //   * O(dx) error, halving with resolution -> only the along-seam offset is left, i.e.
  //     the exchange is topologically right and the halo needs interpolation.
  if (cs_iprob == 6) {
    const Real amp6 = cs_amp;  const Real r0 = cs_r0;
    std::cout << "### CS GHOST CHECK (iprob=6, rho = f(x,y,z), v=0)" << std::endl;
    std::cout << "  panel   active-err       -x2          +x2          -x3"
              << "          +x3       (max |rho-exact| in the FACE ghosts)" << std::endl;
    Real gmax = 0.0;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      const int p = mbpanel.h_view(m);
      Real eact = 0.0, ef[4] = {0.0, 0.0, 0.0, 0.0};
      int ai=-1, aj=-1, ak=-1; Real eahave=0.0, eawant=0.0;
      for (int i=is; i<=ie; ++i) {
        const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                      size.h_view(m).x1max);
        // exact value at the cell (j,k) of THIS panel, ghost indices included
        auto exact = [&](int j, int k) {
          const Real xi  = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                                  size.h_view(m).x2max);
          const Real eta = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                                  size.h_view(m).x3max);
          Real cx, cy, cz;
          PanelToCart(p, xi, eta, cx, cy, cz);
          return CartProfile(cx, cy, cz, rad, d0, r0, amp6);
        };
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real ee = fabs(w0_h(m,IDN,k,j,i) - exact(j,k));
            if (ee > eact) { eact = ee; ai=i-is; aj=j-js; ak=k-ks;
                             eahave = w0_h(m,IDN,k,j,i); eawant = exact(j,k); }
          }
          for (int g=0; g<ng; ++g) {
            ef[0] = fmax(ef[0], fabs(w0_h(m,IDN,k,js-1-g,i) - exact(js-1-g,k)));
            ef[1] = fmax(ef[1], fabs(w0_h(m,IDN,k,je+1+g,i) - exact(je+1+g,k)));
          }
        }
        for (int j=js; j<=je; ++j) {
          for (int g=0; g<ng; ++g) {
            ef[2] = fmax(ef[2], fabs(w0_h(m,IDN,ks-1-g,j,i) - exact(j,ks-1-g)));
            ef[3] = fmax(ef[3], fabs(w0_h(m,IDN,ke+1+g,j,i) - exact(j,ke+1+g)));
          }
        }
      }
      for (int f=0; f<4; ++f) gmax = fmax(gmax, ef[f]);
      std::printf("  %4d   %12.5e  %11.4e  %11.4e  %11.4e  %11.4e  act@(%d,%d,%d)"
                  " have %.6f want %.6f\n",
                  p, eact, ef[0], ef[1], ef[2], ef[3], ai, aj, ak, eahave, eawant);
    }
    std::printf("  GHOST-MAX  %12.5e\n", gmax);
    return;
  }

  // --- iprob = 5: is the panel halo corrupting the RADIAL index? ----------------------
  // rho = rho(r) alone, so any ANGULAR remapping leaves it invariant; a nonzero error in
  // the x2/x3 face ghosts is direct evidence the exchange is mishandling x1. Corner
  // ghosts are excluded (the other tangential index is held in its active range).
  std::cout << "### CS GHOST CHECK (iprob=5, rho = rho(r) only)" << std::endl;
  std::cout << "  panel   active-err     x2-ghost-err   x3-ghost-err" << std::endl;

  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    Real eact = 0.0, eg2 = 0.0, eg3 = 0.0;
    for (int i=is; i<=ie; ++i) {
      const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                    size.h_view(m).x1max);
      const Real ex = RadialProfile(rad, d0, amp, r0);
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          eact = fmax(eact, fabs(w0_h(m,IDN,k,j,i) - ex));
        }
        for (int g=0; g<ng; ++g) {
          eg2 = fmax(eg2, fabs(w0_h(m,IDN,k,js-1-g,i) - ex));
          eg2 = fmax(eg2, fabs(w0_h(m,IDN,k,je+1+g,i) - ex));
        }
      }
      for (int j=js; j<=je; ++j) {
        for (int g=0; g<ng; ++g) {
          eg3 = fmax(eg3, fabs(w0_h(m,IDN,ks-1-g,j,i) - ex));
          eg3 = fmax(eg3, fabs(w0_h(m,IDN,ke+1+g,j,i) - ex));
        }
      }
    }
    std::printf("  %4d   %12.5e   %12.5e   %12.5e\n",
                mbpanel.h_view(m), eact, eg2, eg3);
  }
  return;
}
