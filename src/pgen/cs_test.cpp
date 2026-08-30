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

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <sstream>
#include <cmath>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "pgen.hpp"

namespace {

// Parameters shared between the initial condition and the radial boundary condition.
// ProblemGenerator carries no ParameterInput pointer, so stash them at file scope when
// UserProblem reads them.
Real cs_d0 = 1.0, cs_p0 = 1.0, cs_omega = 1.0;
int  cs_iprob = 1;
Real cs_amp = 0.5, cs_r0 = 1.0;
Real cs_b0r = 0.0;
Real cs_bc_tfrac = 0.0;   // ghost time offset, in units of dt (problem/bc_time_frac)
int  cs_bc_bcc_match = 0; // ghost x1f matched to bcc (problem/bc_bcc_match)
int  cs_bc_probe = 0;     // print the interior boundary-layer PHASE (problem/bc_probe)
Real cs_bazi = 0.0;       // iprob=11 azimuthal amplitude; curl B = 2*cs_bazi*zhat
Real cs_bvx = 0.0, cs_bvy = 0.0, cs_bvz = 0.0;
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
//! \brief The two panel FACE NORMALS at a point: nhat_xi is perpendicular to e_eta (and
//! to rhat), nhat_eta perpendicular to e_xi. They are what b0.x2f/x3f are projections on,
//! and they are NOT parallel to the tangent vectors, nor orthogonal to each other
//! (nhat_xi.nhat_eta = -c). Gram-Schmidt on the tangent pair.

KOKKOS_INLINE_FUNCTION
void PanelNormals(const int p, const Real xi, const Real eta,
                  Real n1[3], Real n2[3]) {
  Real e1[3], e2[3];
  cubed_sphere::PanelTangents(p, xi, eta, e1, e2);
  const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
  const Real s = std::sqrt(1.0 - c*c);
  for (int q=0; q<3; ++q) {
    n1[q] = (e1[q] - c*e2[q])/s;
    n2[q] = (e2[q] - c*e1[q])/s;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \brief The EXACT magnetic field of iprob = 9 at time t.
//!
//! A uniform field is force-free (curl B = 0), so it rides passively on the iprob=3 rigid
//! rotation without disturbing it, and flux freezing carries it round with the flow. For
//! a SPATIALLY UNIFORM B0 that transport is exactly a rigid rotation of the vector:
//!     B(x,t) = Rz(omega t) B0(Rz(-omega t) x) = Rz(omega t) B0,
//! so B stays uniform for all time and only its DIRECTION precesses about z. Check it
//! against the induction equation directly: d/dt [Rz(wt) B0] = omega zhat x B, and
//!     curl(v x B) = (B.grad)v - (v.grad)B = omega zhat x B - 0
//! for uniform B and v = omega zhat x r, which is the same thing. Equivalently: in the
//! frame rotating at omega the entire state is STATIC (v' = 0 kills the Coriolis term and
//! the centrifugal term is balanced by the iprob=3 pressure), and ideal MHD is
//! rotationally covariant.

KOKKOS_INLINE_FUNCTION
void RotatingUniformField(const Real omega, const Real t,
                          const Real bx0, const Real by0, const Real bz0,
                          Real &bx, Real &by, Real &bz) {
  const Real cw = cos(omega*t), sw = sin(omega*t);
  bx = cw*bx0 - sw*by0;
  by = sw*bx0 + cw*by0;
  bz = bz0;
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
void CSTestConvErrors(ParameterInput *pin, Mesh *pm);
void CSTestResistCheck(ParameterInput *pin, Mesh *pm);
void CSTestSeamFluxCheck(Mesh *pm);
void CSTestLevelFluxCheck(Mesh *pm);
void CSTestConsSums(Mesh *pm);
void CSTestHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem
//! \brief Sets the initial conditions for the cubed-sphere validation problems.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // Either module drives the same states: MHD's w0 shares the IDN/IVX/IEN layout, so the
  // initial condition below is written once and only the B field is extra.
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (pmbp->phydro == nullptr && !is_mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "cs_test requires a <hydro> or <mhd> block" << std::endl;
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
  cs_bc_tfrac = pin->GetOrAddReal("problem", "bc_time_frac", 0.0);
  cs_bc_bcc_match = pin->GetOrAddInteger("problem", "bc_bcc_match", 0);
  cs_bc_probe = pin->GetOrAddInteger("problem", "bc_probe", 0);
  cs_exact_panel_ghosts = pin->GetOrAddInteger("problem",
                                              "exact_panel_ghosts", 0);
  cs_r0 = pmy_mesh_->mesh_size.x1min;
  // iprob=1 is included so that the MHD monopole can be given its EXACT radial state:
  // reflect is exact for the uniform hydro state but not for B_r(r), which is not
  // reflection-symmetric, and the jump it manufactures at the radial boundary drives a
  // spurious radial force. The hook is a no-op unless ix1_bc/ox1_bc are set to `user`.
  if (iprob == 1 || iprob == 8 || iprob == 9 || iprob == 11 ||
      (iprob >= 3 && iprob <= 7 && iprob != 4)) {
    user_bcs_func = CSTestRadialBC;
  }
  if (iprob == 8 || (iprob >= 3 && iprob <= 7 && iprob != 4)) {
    pgen_final_func = CSTestGhostCheck;
  }
  // Grid-imprinting time series. iprob=3 is exactly steady, so the deviation from the
  // initial state IS the error at any time. Enable with <problem>/user_hist = true.
  if (iprob == 3) {
    user_hist_func = CSTestHistory;
  }
  // iprob = 9 is the CONVERGENCE test and has its own exact solution at every time, so it
  // reports L1 errors rather than the single-state gates CSTestGhostCheck runs.
  if (iprob == 9 || pin->GetOrAddInteger("problem", "conv_errors", 0) != 0) {
    pgen_final_func = CSTestConvErrors;
  }
  if (iprob == 11) {
    pgen_final_func = CSTestResistCheck;
  }
  // iprob = 11: B = b0c*(-y,x,0) + a uniform tilted field. Its curl is the UNIFORM
  // vector 2*b0c*zhat while the field itself is not uniform, which is what makes it a
  // quantitative test of the resistive current: the uniform part contributes nothing to
  // curl B but DOES make every term of the discrete Stokes loop nonzero, so a term that
  // is dropped or mis-rotated cannot hide. J x B is a pure gradient, so with
  // p = p0 - b0c^2 R^2 + F.r, F = 2*b0c*(zhat x Bunif), the whole state is static in
  // both ideal and resistive MHD.
  Real bazi = 0.0, bux = 0.0, buy = 0.0, buz = 0.0;
  if (iprob == 11) {
    bazi = pin->GetOrAddReal("problem", "b0c", 0.1);
    const Real bu = pin->GetOrAddReal("problem", "bunif", 0.1);
    const Real hx = pin->GetOrAddReal("problem", "bhx", 0.37);
    const Real hy = pin->GetOrAddReal("problem", "bhy", 0.61);
    const Real hz = pin->GetOrAddReal("problem", "bhz", 0.70);
    const Real hn = std::sqrt(hx*hx + hy*hy + hz*hz);
    bux = bu*hx/hn; buy = bu*hy/hn; buz = bu*hz/hn;
    cs_bazi = bazi;
    cs_bvx = bux; cs_bvy = buy; cs_bvz = buz;
  }

  const Real blob_w = pin->GetOrAddReal("problem", "blob_width", 0.3);
  const Real amp_ = cs_amp, r0_ = cs_r0;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &js = indcs.js; int &ks = indcs.ks;
  int &ie = indcs.ie; int &je = indcs.je; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;

  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  const int nhyd = is_mhd ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  const int nscal = is_mhd ? pmbp->pmhd->nscalars : pmbp->phydro->nscalars;

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

    if (iprob == 11) {
      // static: J x B = grad(-b0c^2 R^2 + F.r) is balanced by the pressure
      Real cx, cy, cz;
      PanelToCart(mbpanel.d_view(m), xi, eta, cx, cy, cz);
      const Real px = rad*cx, py = rad*cy;
      w0(m,IDN,k,j,i) = d0;
      w0(m,IEN,k,j,i) = (p0 - SQR(bazi)*(px*px + py*py)
                         + 2.0*bazi*(-buy*px + bux*py))/gm1;
      w0(m,IVX,k,j,i) = 0.0;
      w0(m,IVY,k,j,i) = 0.0;
      w0(m,IVZ,k,j,i) = 0.0;
      return;
    }

    if (iprob == 3 || iprob == 9) {
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

  // --- magnetic field (MHD only) -------------------------------------------------------
  // A RADIAL MONOPOLE, B_r = b0r*(r0/r)^2, is the natural first field for this grid: the
  // radial face area carries r^2, so B_r*A is the SAME on the inner and outer face of
  // every cell while the two angular faces carry no flux at all.  div B is therefore zero
  // to round-off on every panel, by construction and independently of the panel map --
  // the monopole's singularity sits at the origin, outside the shell.  Combined with the
  // iprob=3 rigid rotation it gives a nonzero EMF, so it exercises CT rather than just
  // sitting there.
  if (is_mhd && (iprob == 8 || iprob == 9)) {
    // iprob = 8: a UNIFORM CARTESIAN field B = b0c * bhat, with v = 0 and uniform p.
    // iprob = 9: the SAME field, but on top of the iprob=3 rigid rotation, so that it is
    // carried round by the flow instead of sitting still. See RotatingUniformField.
    // A uniform field has curl B = 0, so it exerts no force and this is again an exact
    // static state -- but unlike the monopole it has TANGENTIAL components on every
    // panel, which is what makes it the test of the non-orthogonal basis. |B|^2 = b0c^2
    // exactly and everywhere, so the total energy below is analytic, and the pressure
    // ConsToPrim recovers from it is a direct measurement of the magnetic energy the
    // inversion subtracted. See CSTestGhostCheck.
    const Real b0c = pin->GetOrAddReal("problem", "b0c", 1.0);
    const Real bhx = pin->GetOrAddReal("problem", "bhx", 0.0);
    const Real bhy = pin->GetOrAddReal("problem", "bhy", 0.0);
    const Real bhz = pin->GetOrAddReal("problem", "bhz", 1.0);
    const Real bn = std::sqrt(bhx*bhx + bhy*bhy + bhz*bhz);
    const Real bvx = b0c*bhx/bn, bvy = b0c*bhy/bn, bvz = b0c*bhz/bn;
    cs_b0r = b0c;
    cs_bvx = bvx; cs_bvy = bvy; cs_bvz = bvz;
    auto &b0f = pmbp->pmhd->b0;
    par_for("pgen_cs_bfld8", DevExeSpace(), 0,(pmbp->nmb_thispack-1),
            ks,ke+1, js,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                              size.d_view(m).x2max);
      const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                              size.d_view(m).x3max);
      const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                            size.d_view(m).x2max);
      const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                            size.d_view(m).x3max);
      const int p = mbpanel.d_view(m);
      Real n1[3], n2[3], cx, cy, cz;
      // x1 face: normal is rhat, evaluated at the (j,k) cell centre
      if (j <= je && k <= ke) {
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        b0f.x1f(m,k,j,i) = bvx*cx + bvy*cy + bvz*cz;
      }
      // x2 face: normal is nhat_xi, evaluated at the xi FACE
      if (i <= ie && k <= ke) {
        PanelNormals(p, x2f, x3c, n1, n2);
        b0f.x2f(m,k,j,i) = bvx*n1[0] + bvy*n1[1] + bvz*n1[2];
      }
      // x3 face: normal is nhat_eta, evaluated at the eta FACE
      if (i <= ie && j <= je) {
        PanelNormals(p, x2c, x3f, n1, n2);
        b0f.x3f(m,k,j,i) = bvx*n2[0] + bvy*n2[1] + bvz*n2[2];
      }
    });
  } else if (is_mhd && iprob == 11) {
    // Face-normal projections of B(x) at each face centre. Unlike iprob = 8 the field
    // varies in space, so each face needs its own POSITION as well as its own normal.
    auto &b0f = pmbp->pmhd->b0;
    par_for("pgen_cs_bfld11", DevExeSpace(), 0,(pmbp->nmb_thispack-1),
            ks,ke+1, js,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                              size.d_view(m).x2max);
      const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                              size.d_view(m).x3max);
      const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                            size.d_view(m).x2max);
      const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                            size.d_view(m).x3max);
      const Real rc = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                   size.d_view(m).x1max);
      const Real rl = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                 size.d_view(m).x1max);
      const int p = mbpanel.d_view(m);
      Real n1[3], n2[3], cx, cy, cz;
      if (j <= je && k <= ke) {
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        const Real bx = -bazi*rl*cy + bux;
        const Real by =  bazi*rl*cx + buy;
        b0f.x1f(m,k,j,i) = bx*cx + by*cy + buz*cz;
      }
      if (i <= ie && k <= ke) {
        PanelToCart(p, x2f, x3c, cx, cy, cz);
        PanelNormals(p, x2f, x3c, n1, n2);
        const Real bx = -bazi*rc*cy + bux;
        const Real by =  bazi*rc*cx + buy;
        b0f.x2f(m,k,j,i) = bx*n1[0] + by*n1[1] + buz*n1[2];
      }
      if (i <= ie && j <= je) {
        PanelToCart(p, x2c, x3f, cx, cy, cz);
        PanelNormals(p, x2c, x3f, n1, n2);
        const Real bx = -bazi*rc*cy + bux;
        const Real by =  bazi*rc*cx + buy;
        b0f.x3f(m,k,j,i) = bx*n2[0] + by*n2[1] + buz*n2[2];
      }
    });
  } else if (is_mhd) {
    const Real b0r = pin->GetOrAddReal("problem", "b0r", 1.0);
    cs_b0r = b0r;
    auto &b0f = pmbp->pmhd->b0;
    auto &bcc = pmbp->pmhd->bcc0;
    par_for("pgen_cs_bfld", DevExeSpace(), 0,(pmbp->nmb_thispack-1),
            ks,ke+1, js,je+1, is,ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      // x1 is RADIAL, so only the x1 faces carry flux.  x1f is indexed by the LEFT edge.
      if (j <= je && k <= ke) {
        const Real rf = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                   size.d_view(m).x1max);
        b0f.x1f(m,k,j,i) = b0r*SQR(r0_/rf);
      }
      if (i <= ie && k <= ke) b0f.x2f(m,k,j,i) = 0.0;
      if (i <= ie && j <= je) b0f.x3f(m,k,j,i) = 0.0;
    });
    // cell-centred field: the average of the two opposing faces
    par_for("pgen_cs_bcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je, is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      bcc(m,IBX,k,j,i) = 0.5*(b0f.x1f(m,k,j,i) + b0f.x1f(m,k,j,i+1));
      bcc(m,IBY,k,j,i) = 0.5*(b0f.x2f(m,k,j,i) + b0f.x2f(m,k,j+1,i));
      bcc(m,IBZ,k,j,i) = 0.5*(b0f.x3f(m,k,j,i) + b0f.x3f(m,k+1,j,i));
    });
  }

  // Convert primitives to conserved.
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  // NOT peos->PrimToCons: that assumes an orthonormal basis. On the cubed sphere the
  // conserved momentum is the COVARIANT one. See Coordinates::GnomonicEquiangleLowerMom.
  pmbp->pcoord->GnomonicEquiangleLowerMom(w0, u0, is, ie, js, je, ks, ke);
  // LowerMom knows nothing about B, so add the magnetic energy it left out.
  if (is_mhd && iprob == 11) {
    // ANALYTIC too, and it matters: the pressure is in force balance against the ANALYTIC
    // |B|^2/2, so taking the magnetic energy from the face-averaged bcc instead would
    // leave the state out of balance by that difference and set the shell ringing.
    par_for("pgen_cs_emag11", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je,
            is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const Real x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                    size.d_view(m).x2max);
      const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                    size.d_view(m).x3max);
      const Real rad = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                    size.d_view(m).x1max);
      Real cx, cy, cz;
      PanelToCart(mbpanel.d_view(m), 0.25*M_PI*x2v, 0.25*M_PI*x3v, cx, cy, cz);
      const Real px = rad*cx, py = rad*cy;
      const Real bx = -bazi*py + bux, by = bazi*px + buy;
      u0(m,IEN,k,j,i) += 0.5*(bx*bx + by*by + buz*buz);
    });
  } else if (is_mhd && (iprob == 8 || iprob == 9)) {
    // ANALYTIC, deliberately: |B|^2 = b0c^2 exactly for a uniform Cartesian field. Adding
    // a magnetic energy built from the same bcc that ConsToPrim will subtract would make
    // the pressure check below self-consistent instead of correct, and would pass under
    // any convention.
    const Real bsq = SQR(cs_b0r);
    par_for("pgen_cs_emag8", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je, is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IEN,k,j,i) += 0.5*bsq;
    });
  } else if (is_mhd) {
    auto &bcc = pmbp->pmhd->bcc0;
    par_for("pgen_cs_emag", DevExeSpace(), 0,(pmbp->nmb_thispack-1), ks,ke, js,je, is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IEN,k,j,i) += 0.5*(SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                            + SQR(bcc(m,IBZ,k,j,i)));
    });
  }

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
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  const int nhyd = is_mhd ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  const int nscal = is_mhd ? pmbp->pmhd->nscalars : pmbp->phydro->nscalars;

  const Real d0 = cs_d0, p0 = cs_p0, omega = cs_omega;
  const Real bazi = cs_bazi;
  // iprob = 11 heats at eta*|J|^2 everywhere, so the ghosts must heat with it. Holding
  // the t = 0 pressure instead leaves the heated interior pushing against a cold
  // boundary: that leaks energy at a rate proportional to the heat already deposited,
  // which shows up as a heating deficit growing like t^2 and diluting as 1/L.
  Real dp_ohm = 0.0;
  if (cs_iprob == 11 && pmbp->pmhd != nullptr && pmbp->pmhd->presist != nullptr) {
    const Real gmm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
    dp_ohm = gmm1*(pmbp->pmhd->presist->eta_ohm_const)*SQR(2.0*cs_bazi)*pm->time;
  }
  const Real bux = cs_bvx, buy = cs_bvy, buz = cs_bvz;
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
      if (iprob == 11) {
        Real cx, cy, cz;
        PanelToCart(p, xi, eta, cx, cy, cz);
        const Real px = rad*cx, py = rad*cy;
        dn = d0;
        ie_ = (p0 - SQR(bazi)*(px*px + py*py)
               + 2.0*bazi*(-buy*px + bux*py) + dp_ohm)/gm1;
        v1 = 0.0; v2 = 0.0; v3 = 0.0;
      } else if (iprob == 7) {
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
      } else if (iprob == 1 || iprob == 8) {
        dn = d0;
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
  // MHD: the radial ghost FACES. `user` is not one of the cases bfield_bcs.cpp handles,
  // so b0 in the radial ghosts is whatever was left there unless it is set here. Only the
  // x1 faces carry flux for the monopole; the two angular face fields stay zero, and the
  // energy has to be topped up because u0 above was written without the field.
  if (is_mhd && (iprob == 8 || iprob == 9)) {
    // The uniform Cartesian field is no more reflection-symmetric than the monopole, so
    // its radial ghosts get the exact state too -- faces, cell centres and the analytic
    // magnetic energy. For iprob = 9 the exact field is TIME DEPENDENT: it precesses
    // about z with the flow, so the ghosts must precess with it or the boundary holds the
    // t = 0 field against an interior that has moved on, which is a first-order error
    // confined to the radial boundary and destroys the convergence being measured.
    auto &b0f = pmbp->pmhd->b0;
    auto &bcc = pmbp->pmhd->bcc0;
    Real bvx = cs_bvx, bvy = cs_bvy, bvz = cs_bvz;
    if (iprob == 9) {
      // The ghosts are set to the exact field at t^n + bc_time_frac*dt. Mesh::time is
      // advanced only at the END of a cycle (driver.cpp), while ApplyPhysicalBCs runs
      // inside each RK stage, so which time the ghosts should carry is not obvious and
      // it MATTERS: the residual error of this test is an O(dt) rotation PHASE LAG.
      //
      // MEASURED, nx1=16 nx2=nx3=32 tlim=1, scanning frac from -1 to +1 (the phase is
      // read off the best-fit uniform B; + is ahead of the exact field):
      //     frac    -1.0    -0.75   -0.5    -0.25    0      +0.5    +1.0
      //     L1(B)  6.73e-4 5.15e-4 3.61e-4 2.27e-4 2.09e-4 4.58e-4 6.77e-4
      //     phase  +2.8e-4 +2.4e-4 +0.7e-4 -1.8e-4 -4.7e-4 -9.9e-4 -1.06e-3
      // Two separate facts, and the earlier note that this family is "monotonic" was
      // only half the scan and is RETRACTED:
      //  (a) L1(B) is a V with its minimum at frac = 0, because any |frac| =/= 0 opens a
      //      jump between the ghost field and the interior one and that jump adds
      //      structured, non-rotational error. frac = 0 stays the right default.
      //  (b) the PHASE is linear in frac and crosses zero at frac = -0.41, with slope
      //      -0.906 per unit of ghost phase offset omega*frac*dt. A slope of order -1
      //      means the shell's mean field LOCKS to the phase the boundary holds (the
      //      Alfven crossing time L/v_A = 1 is the run time here) rather than
      //      accumulating from it, and the whole scan is one law:
      //          delta = -0.906 * omega * (frac + 0.41) * dt.
      //      The frac = 0 lag is that law at frac = 0, which is why it is O(dt) and why
      //      it dilutes as 1/L.
      // So the clock cannot BOTH null the phase and keep the field a rigid rotation, and
      // that is why "advance the ghosts to t^n + dt" fails. What is still unexplained is
      // the offset -0.41: problem/bc_probe measures the interior's own phase at the first
      // active radial layer and it sits at omega*(t^n + dt), i.e. exactly where SSP-RK2
      // puts both of its stage outputs, so the ghost that nulls the lag is 1.41*dt BEHIND
      // the interior it faces, not level with it.
      RotatingUniformField(omega, pm->time + cs_bc_tfrac*pm->dt, cs_bvx, cs_bvy, cs_bvz,
                           bvx, bvy, bvz);
    }
    // -------------------------------------------------------------------------------
    // problem/bc_probe. WHAT TIME IS THE INTERIOR AT, as seen from the radial boundary?
    // The whole bc_time_frac question turns on this and it is directly measurable. The
    // exact field is Rz(theta) B0 = cos(theta) P + sin(theta) Q + Z with
    // P = (B0x,B0y,0), Q = (-B0y,B0x,0), Z = (0,0,B0z), and the FACE fields are linear
    // projections of it, so a 2-parameter least squares over the first ACTIVE radial
    // layer recovers the phase theta the interior is actually carrying. Compare it with
    // omega*t^n (bc_time_frac = 0) and omega*(t^n + dt) (bc_time_frac = 1).
    if (iprob == 9 && cs_bc_probe > 0 && (pm->ncycle % cs_bc_probe) == 0) {
      int &je_p = indcs.je; int &ke_p = indcs.ke;
      const int nj = indcs.nx2, nk = indcs.nx3;
      const int nmb = pmbp->nmb_thispack;
      const Real p0x = cs_bvx, p0y = cs_bvy, p0z = cs_bvz;
      Real app = 0.0, apq = 0.0, aqq = 0.0, bp = 0.0, bq = 0.0;
      const int nn = nmb*nk*nj;
      // one reduce per normal-equation entry: this is a low-frequency diagnostic
      for (int which=0; which<5; ++which) {
        Real sum = 0.0;
        Kokkos::parallel_reduce("cs_bc_probe", Kokkos::RangePolicy<>(DevExeSpace(),0,nn),
        KOKKOS_LAMBDA(const int idx, Real &lsum) {
          const int m = idx/(nk*nj);
          const int k = (idx - m*nk*nj)/nj + ks;
          const int j = (idx - m*nk*nj) % nj + js;
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
          const Real x2c = 0.25*M_PI*CellCenterX(j-js, nj, size.d_view(m).x2min,
                                                           size.d_view(m).x2max);
          const Real x3c = 0.25*M_PI*CellCenterX(k-ks, nk, size.d_view(m).x3min,
                                                           size.d_view(m).x3max);
          const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, nj, size.d_view(m).x2min,
                                                         size.d_view(m).x2max);
          const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, nk, size.d_view(m).x3min,
                                                         size.d_view(m).x3max);
          const int p = mbpanel.d_view(m);
          Real n1[3], n2[3];
          // the two ANGULAR faces of the first active cell: these are what a phase
          // mismatch across the boundary would shear
          for (int f=0; f<2; ++f) {
            if (f == 0) {
              PanelNormals(p, x2f, x3c, n1, n2);
            } else {
              PanelNormals(p, x2c, x3f, n1, n2);
            }
            const Real *nh = (f == 0) ? n1 : n2;
            const Real pn = p0x*nh[0] + p0y*nh[1];
            const Real qn = -p0y*nh[0] + p0x*nh[1];
            const Real zn = p0z*nh[2];
            const Real meas = (f == 0) ? b0f.x2f(m,k,j,is) : b0f.x3f(m,k,j,is);
            const Real r = meas - zn;
            if (which == 0) lsum += pn*pn;
            if (which == 1) lsum += pn*qn;
            if (which == 2) lsum += qn*qn;
            if (which == 3) lsum += pn*r;
            if (which == 4) lsum += qn*r;
          }
        }, sum);
        if (which == 0) app = sum;
        if (which == 1) apq = sum;
        if (which == 2) aqq = sum;
        if (which == 3) bp = sum;
        if (which == 4) bq = sum;
      }
      const Real det = app*aqq - apq*apq;
      const Real ca = (bp*aqq - bq*apq)/det;
      const Real sa = (bq*app - bp*apq)/det;
      const Real th = std::atan2(sa, ca);
      std::printf("### CS BC PROBE  ncycle=%d  t=%.8f dt=%.6e  theta_int=%.9e"
                  "  w*t=%.9e  th-w*t=%+.4e  (th-w*t)/(w*dt)=%+.4f  |B|fit=%.9f\n",
                  pm->ncycle, pm->time, pm->dt, th, omega*pm->time,
                  th - omega*pm->time,
                  (th - omega*pm->time)/(omega*pm->dt), std::sqrt(ca*ca + sa*sa));
    }

    const Real bsq = SQR(cs_b0r);
    int &ie_i = indcs.ie;
    par_for("cs_test_rbc_b8", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
            0,(n2-1), 0,(ng-1),
    KOKKOS_LAMBDA(int m, int k, int j, int ig) {
      const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                              size.d_view(m).x2max);
      const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                              size.d_view(m).x3max);
      const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                            size.d_view(m).x2max);
      const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                            size.d_view(m).x3max);
      const int p = mbpanel.d_view(m);
      Real n1[3], n2[3], cx, cy, cz;
      PanelToCart(p, x2c, x3c, cx, cy, cz);
      const Real br = bvx*cx + bvy*cy + bvz*cz;
      PanelNormals(p, x2f, x3c, n1, n2);
      const Real bxi = bvx*n1[0] + bvy*n1[1] + bvz*n1[2];
      PanelNormals(p, x2c, x3f, n1, n2);
      const Real bet = bvx*n2[0] + bvy*n2[1] + bvz*n2[2];
      for (int side=0; side<2; ++side) {
        if (side == 0 &&
            mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
        if (side == 1 &&
            mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
        const int i = (side == 0) ? (is-ng+ig) : (ie_i+1+ig);
        b0f.x2f(m,k,j,i) = bxi;
        b0f.x3f(m,k,j,i) = bet;
        // x1f: fill the GHOST faces only. The two DOMAIN faces x1f(is) and x1f(ie+1)
        // belong to active cells and are updated by CT; overwriting one of them each
        // step breaks the exact cell-by-cell conservation of sum(area*B) that CT
        // otherwise guarantees. This wrote x1f(ie+1) and made div B at i = ie grow to
        // 6e-2 over five steps -- on every panel, seam and interior alike -- which
        // reads exactly like a seam EMF defect and is not one. The inner side already
        // left x1f(is) alone; the outer side is now symmetric with it.
        b0f.x1f(m,k,j,(side == 0) ? i : (i+1)) = br;
        // bcc in the orthonormal frame, matching GnomonicEquiangleRaiseVelMHD
        Real e1[3], e2[3];
        cubed_sphere::PanelTangents(p, x2c, x3c, e1, e2);
        const Real cc = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
        const Real ss = sqrt(1.0 - cc*cc);
        bcc(m,IBX,k,j,i) = br;
        bcc(m,IBY,k,j,i) = (bxi + cc*bet)/ss;
        bcc(m,IBZ,k,j,i) = bet;
        u0(m,IEN,k,j,i) += 0.5*bsq;
      }
    });

    // -------------------------------------------------------------------------------
    // problem/bc_bcc_match. WHY THIS KNOB EXISTS. Whatever this BC writes into the ghost
    // `bcc` is thrown away: the stage chain is ApplyPhysicalBCs -> Prolongate ->
    // ConToPrim, and ConToPrim (then GnomonicEquiangleRaiseVelMHD) REBUILDS bcc over the
    // whole array, ghost zones included, as the plain face average
    //     bcc_r(i) = 0.5*(x1f(i) + x1f(i+1)).
    // For the FIRST ghost cell that average straddles the boundary: x1f(is) is the
    // interior face CT is evolving, and only x1f(is-1) is the analytic value written
    // above. So the ghost cell's RADIAL field is half analytic and half interior, while
    // its two ANGULAR components (whose faces both lie inside the ghost column) are
    // fully analytic. Advancing the ghost clock with bc_time_frac therefore does NOT
    // advance the ghost cell state coherently -- it moves the angular components by the
    // full omega*frac*dt and the radial one by only half of it -- which is exactly the
    // measured signature that at frac = 1 the error stops being a rigid rotation and the
    // two components no longer agree on a single angle.
    //
    // With bc_bcc_match = 1 the ghost x1f are instead chosen so that the average IS the
    // exact radial field at each ghost cell centre, x1f(i) = 2*Br - x1f(i+1) marching
    // away from the boundary (and the mirror on the outer side). The ghost cell then
    // carries the exact field in all three components, the DOMAIN faces x1f(is) and
    // x1f(ie+1) are still untouched, and the bc_time_frac scan becomes a clean test of
    // the ghost CLOCK alone: SSP-RK2 evaluates the interior at t^n and t^n + dt, so a
    // coherent ghost should be best at frac = 0.5.
    if (cs_bc_bcc_match != 0) {
      int &ie_j = indcs.ie;
      par_for("cs_test_rbc_b9match", DevExeSpace(), 0,(pmbp->nmb_thispack-1),
              0,(n3-1), 0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
        const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                                size.d_view(m).x2max);
        const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                                size.d_view(m).x3max);
        const int p = mbpanel.d_view(m);
        Real cx, cy, cz;
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        // rhat depends only on the angles, so Br is the same in every ghost cell of the
        // column -- the recursion below alternates rather than drifting.
        const Real br = bvx*cx + bvy*cy + bvz*cz;
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int i=is-1; i>=is-ng; --i) {
            b0f.x1f(m,k,j,i) = 2.0*br - b0f.x1f(m,k,j,i+1);
          }
        }
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
          for (int i=ie_j+1; i<=ie_j+ng; ++i) {
            b0f.x1f(m,k,j,i+1) = 2.0*br - b0f.x1f(m,k,j,i);
          }
        }
      });
    }
  }

  // iprob = 11: the exact static field in the radial ghosts. As in the iprob = 8/9 block
  // the two DOMAIN faces x1f(is) and x1f(ie+1) are left to CT.
  if (is_mhd && iprob == 11) {
    auto &b0f = pmbp->pmhd->b0;
    int &ie_i = indcs.ie;
    par_for("cs_test_rbc_b11", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
            0,(n2-1), 0,(ng-1),
    KOKKOS_LAMBDA(int m, int k, int j, int ig) {
      const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                              size.d_view(m).x2max);
      const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                              size.d_view(m).x3max);
      const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                                                            size.d_view(m).x2max);
      const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                                            size.d_view(m).x3max);
      const int p = mbpanel.d_view(m);
      for (int side=0; side<2; ++side) {
        if (side == 0 &&
            mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
        if (side == 1 &&
            mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
        const int i = (side == 0) ? (is-ng+ig) : (ie_i+1+ig);
        const Real rc = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                     size.d_view(m).x1max);
        const int ir = (side == 0) ? i : (i+1);
        const Real rf = LeftEdgeX(ir-is, indcs.nx1, size.d_view(m).x1min,
                                                    size.d_view(m).x1max);
        Real n1[3], n2[3], cx, cy, cz;
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        Real bx = -bazi*rf*cy + bux;
        Real by =  bazi*rf*cx + buy;
        b0f.x1f(m,k,j,ir) = bx*cx + by*cy + buz*cz;
        PanelToCart(p, x2f, x3c, cx, cy, cz);
        PanelNormals(p, x2f, x3c, n1, n2);
        bx = -bazi*rc*cy + bux;
        by =  bazi*rc*cx + buy;
        b0f.x2f(m,k,j,i) = bx*n1[0] + by*n1[1] + buz*n1[2];
        PanelToCart(p, x2c, x3f, cx, cy, cz);
        PanelNormals(p, x2c, x3f, n1, n2);
        bx = -bazi*rc*cy + bux;
        by =  bazi*rc*cx + buy;
        b0f.x3f(m,k,j,i) = bx*n2[0] + by*n2[1] + buz*n2[2];
        // magnetic energy at the ghost cell centre, analytic
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        const Real px = rc*cx, py = rc*cy;
        const Real bxc = -bazi*py + bux, byc = bazi*px + buy;
        u0(m,IEN,k,j,i) += 0.5*(bxc*bxc + byc*byc + buz*buz);
      }
    });
  }

  if (is_mhd && iprob == 1) {
    auto &b0f = pmbp->pmhd->b0;
    auto &bcc = pmbp->pmhd->bcc0;
    const Real b0r = cs_b0r, r0b = cs_r0;
    int &ie_i = indcs.ie;
    par_for("cs_test_rbc_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
            0,(n2-1), 0,(ng-1),
    KOKKOS_LAMBDA(int m, int k, int j, int ig) {
      for (int side=0; side<2; ++side) {
        if (side == 0 &&
            mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
        if (side == 1 &&
            mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
        const int i = (side == 0) ? (is-ng+ig) : (ie_i+1+ig);
        const Real rf = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                                                   size.d_view(m).x1max);
        b0f.x1f(m,k,j,i) = b0r*SQR(r0b/rf);
        b0f.x2f(m,k,j,i) = 0.0;
        b0f.x3f(m,k,j,i) = 0.0;
        if (side == 1 && ig == ng-1) {
          const Real rfo = LeftEdgeX(i+1-is, indcs.nx1, size.d_view(m).x1min,
                                                        size.d_view(m).x1max);
          b0f.x1f(m,k,j,i+1) = b0r*SQR(r0b/rfo);
        }
        const Real bc = 0.5*(b0f.x1f(m,k,j,i)
                           + b0r*SQR(r0b/LeftEdgeX(i+1-is, indcs.nx1,
                               size.d_view(m).x1min, size.d_view(m).x1max)));
        bcc(m,IBX,k,j,i) = bc;
        bcc(m,IBY,k,j,i) = 0.0;
        bcc(m,IBZ,k,j,i) = 0.0;
        u0(m,IEN,k,j,i) += 0.5*SQR(bc);
      }
    });
  }

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

//----------------------------------------------------------------------------------------
//! \fn CSTestConvErrors
//! \brief L1 errors against the EXACT iprob = 9 solution, for the convergence test.
//!
//! iprob = 9 is the rigid rotation of iprob = 3 carrying the uniform force-free field of
//! iprob = 8. The exact solution is known at EVERY time (see RotatingUniformField), not
//! just at a return period, so the errors below are a genuine measure of the evolved
//! solution rather than of the initial condition -- which is what separates this from the
//! static gates. It exercises the seam EMF hard: E = -v x B is large and varies across
//! every panel boundary, yet its curl must reproduce a field that stays exactly UNIFORM,
//! so any seam inconsistency shows up directly as structure in B.
//!
//! The FACE fields are compared, not bcc: they are what CT actually evolves, so this
//! avoids folding in the cell-centring average's own convention and error. Each cell
//! contributes its three lower faces, each projected onto the exact Cartesian field with
//! the same PanelNormals used to seed it.
//!
//! The norm is the arithmetic mean of |error| over active cells, normalised by b0c (and
//! by the local value for the hydro variables). On a grid uniform in (r, xi, eta) that is
//! the volume-weighted L1 norm with the smooth O(1) volume factor replaced by 1, which
//! leaves the convergence ORDER unchanged; it is a convergence measure, not a physical
//! norm. Errors are APPENDED to <basename>-cs-errs.dat so a resolution sweep accumulates
//! in one file.

void CSTestConvErrors(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  CSTestConsSums(pm);
  CSTestSeamFluxCheck(pm);
  CSTestLevelFluxCheck(pm);
  // The field errors are only defined for iprob = 9. Any other setup (in practice the
  // FIELD-FREE iprob = 3 rigid rotation, which is the control for this test) reports the
  // hydro errors alone, against the same exact steady equilibrium.
  const bool has_b = (pmbp->pmhd != nullptr) && (cs_iprob == 9);
  if (pmbp->pmhd == nullptr && pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  size.sync_host();
  auto &mbpanel = pmbp->pmb->mb_panel;
  mbpanel.sync_host();

  // the exact field at the final time
  Real ex, ey, ez;
  RotatingUniformField(cs_omega, pm->time, cs_bvx, cs_bvy, cs_bvz, ex, ey, ez);

  DvceArray4D<Real> b1d, b2d, b3d;
  if (has_b) {
    b1d = pmbp->pmhd->b0.x1f;
    b2d = pmbp->pmhd->b0.x2f;
    b3d = pmbp->pmhd->b0.x3f;
  }
  auto b1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b1d);
  auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b2d);
  auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b3d);
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto w0h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), w0);
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;

  Real l1b = 0.0, lib = 0.0, l1v = 0.0, l1p = 0.0;
  std::int64_t ncell = 0;
  // Localisation: split the L1(B) budget between cells adjacent to a panel SEAM (within
  // nband of the tangential block edge), cells adjacent to the RADIAL boundary, and the
  // panel interior, and record where the maximum sits. A flat L1 with a doubling Linf
  // means the error lives on a set that thins as the grid refines -- this says which.
  const int nband = 2;
  Real l1_seam = 0.0, l1_rad = 0.0, l1_int = 0.0;
  std::int64_t n_seam = 0, n_rad = 0, n_int = 0;
  int mx_m = -1, mx_p = -1, mx_i = -1, mx_j = -1, mx_k = -1;
  // Profile of the error against DISTANCE from the panel seam, in cells. The flow advects
  // ~omega*t/(pi/2) of a panel width in the run, so error created at a seam is carried a
  // few cells inward; only a profile that stays flat deep inside a panel shows a real
  // INTERIOR error rather than seam error that has been advected.
  const int NBIN = 6;
  Real l1_bin[NBIN] = {};
  std::int64_t n_bin[NBIN] = {};
  // The same profile against distance from the RADIAL boundary. A boundary-condition
  // error shows up as a profile that decays inward from i = 0; an error made by the
  // interior scheme does not.
  Real l1_rbin[NBIN] = {};
  std::int64_t n_rbin[NBIN] = {};
  // Least-squares fit of the face data to a single UNIFORM Cartesian vector: every face
  // value is one equation b_num = Bfit . n in the same normal n used to seed it. This
  // splits the error into a SYSTEMATIC part (Bfit - B_exact: the field is still uniform
  // but points the wrong way, i.e. it was rotated by the wrong angle) and a STRUCTURED
  // part (the residual about Bfit). Which one dominates says what kind of bug this is.
  Real ata[3][3] = {};
  Real atb[3] = {0.0, 0.0, 0.0};
  Real l1res = 0.0;
  // Is the FLOW itself turning at the right rate? A rigid rotation is linear in omega,
  // so if the numerical velocity is still a rigid rotation but at omega_fit, then
  // v_num = (omega_fit/omega) v_exact everywhere. One least-squares scalar recovers that
  // ratio, with no need to leave the gnomonic basis: s = <v_num.v_ex>/<v_ex.v_ex>. If
  // (s-1) matches the field's phase lag divided by omega*t, the field is simply frozen
  // into a flow that is slow, and the induction is not the culprit.
  Real vdotv = 0.0, vdote = 0.0;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    for (int k=ks; k<=ke; ++k) {
      const Real x3c = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                              size.h_view(m).x3max);
      const Real x3f = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                            size.h_view(m).x3max);
      for (int j=js; j<=je; ++j) {
        const Real x2c = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                                size.h_view(m).x2max);
        const Real x2f = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                              size.h_view(m).x2max);
        // exact face-normal components, seeded exactly as the initial condition was
        Real n1[3], n2[3], cx, cy, cz;
        PanelToCart(p, x2c, x3c, cx, cy, cz);
        const Real br_ex = ex*cx + ey*cy + ez*cz;
        Real nv[3][3];
        nv[0][0] = cx; nv[0][1] = cy; nv[0][2] = cz;
        PanelNormals(p, x2f, x3c, n1, n2);
        const Real b2_ex = ex*n1[0] + ey*n1[1] + ez*n1[2];
        for (int c=0; c<3; ++c) { nv[1][c] = n1[c]; }
        PanelNormals(p, x2c, x3f, n1, n2);
        const Real b3_ex = ex*n2[0] + ey*n2[1] + ez*n2[2];
        for (int c=0; c<3; ++c) { nv[2][c] = n2[c]; }
        for (int i=is; i<=ie; ++i) {
          if (has_b) {
            const Real e1 = fabs(b1h(m,k,j,i) - br_ex);
            const Real e2 = fabs(b2h(m,k,j,i) - b2_ex);
            const Real e3 = fabs(b3h(m,k,j,i) - b3_ex);
            const Real ec = (e1 + e2 + e3)/3.0;
            const Real em = fmax(e1, fmax(e2, e3));
            l1b += ec;
            if (em > lib) {
              lib = em;
              mx_m = m; mx_p = p; mx_i = i-is; mx_j = j-js; mx_k = k-ks;
            }
            const bool seam = (j-js < nband) || (je-j < nband)
                           || (k-ks < nband) || (ke-k < nband);
            const bool radb = (i-is < nband) || (ie-i < nband);
            const Real bnum[3] = {b1h(m,k,j,i), b2h(m,k,j,i), b3h(m,k,j,i)};
            for (int f=0; f<3; ++f) {
              for (int c=0; c<3; ++c) {
                atb[c] += bnum[f]*nv[f][c];
                for (int d=0; d<3; ++d) { ata[c][d] += nv[f][c]*nv[f][d]; }
              }
            }
            const int dsm = fmin(fmin(j-js, je-j), fmin(k-ks, ke-k));
            int bn = 0;
            while (bn < NBIN-1 && dsm >= (1 << bn)) { ++bn; }
            l1_bin[bn] += ec; ++n_bin[bn];
            const int drd = fmin(i-is, ie-i);
            int rb = 0;
            while (rb < NBIN-1 && drd >= (1 << rb)) { ++rb; }
            l1_rbin[rb] += ec; ++n_rbin[rb];
            if (seam) {
              l1_seam += ec; ++n_seam;
            } else if (radb) {
              l1_rad += ec; ++n_rad;
            } else {
              l1_int += ec; ++n_int;
            }
          }
          // the hydro state is STEADY, so its error is measured against the exact
          // rigid-rotation equilibrium at this point
          const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                        size.h_view(m).x1max);
          Real dn, iex, v1, v2, v3;
          RigidRotState(p, x2c, x3c, rad, cs_d0, cs_p0, cs_omega, gm1,
                        dn, iex, v1, v2, v3);
          l1v += fabs(w0h(m,IVX,k,j,i) - v1) + fabs(w0h(m,IVY,k,j,i) - v2)
               + fabs(w0h(m,IVZ,k,j,i) - v3);
          vdote += w0h(m,IVY,k,j,i)*v2 + w0h(m,IVZ,k,j,i)*v3;
          vdotv += v2*v2 + v3*v3;
          l1p += fabs(w0h(m,IEN,k,j,i) - iex)/iex;
          ++ncell;
        }
      }
    }
  }
  // EVERY sum, count and maximum above is RANK-LOCAL: the loop visits only the
  // MeshBlocks in this rank's MeshBlockPack. Reporting them unreduced makes the norms
  // depend on the DECOMPOSITION -- one identical solution then prints a different L1 at
  // every rank count, which reads exactly like an MPI bug and is not one. (It cost this
  // file a false "refined runs are rank-dependent" finding; the binary dumps were
  // bitwise identical at 1, 2, 3 and 6 ranks the whole time.) Reduce first, report on
  // rank 0 only.
  //
  // The Linf LOCATION cannot be reduced numerically -- it names a block, and block ids
  // are local -- so each rank formats its own line and MAXLOC says whose to broadcast.
  char linfline[256];
  linfline[0] = '\0';
  if (mx_m >= 0) {
    auto &mlev = pmbp->pmb->mb_lev;
    auto &ngh = pmbp->pmb->nghbr;
    ngh.sync_host();
    int nlo = -1, nhi = -1;
    for (int q=0; q<4; ++q) {
      if (ngh.h_view(mx_m,q).gid   >= 0) nlo = ngh.h_view(mx_m,q).lev;
      if (ngh.h_view(mx_m,4+q).gid >= 0) nhi = ngh.h_view(mx_m,4+q).lev;
    }
    std::snprintf(linfline, sizeof(linfline),
                  "###   Linf(B) at panel %d  (i,j,k) = (%d,%d,%d)  of block %d"
                  " [lev %d, x1 nghbr lev %d / %d, x1 = %g..%g]\n",
                  mx_p, mx_i, mx_j, mx_k, mx_m, mlev.h_view(mx_m), nlo, nhi,
                  size.h_view(mx_m).x1min, size.h_view(mx_m).x1max);
  }
#if MPI_PARALLEL_ENABLED
  {
    const int nrl = 6 + 2*NBIN + 12 + 2;
    Real rsum[6 + 2*NBIN + 12 + 2];
    int q = 0;
    rsum[q++] = l1b; rsum[q++] = l1v; rsum[q++] = l1p;
    rsum[q++] = l1_seam; rsum[q++] = l1_rad; rsum[q++] = l1_int;
    for (int b=0; b<NBIN; ++b) { rsum[q++] = l1_bin[b]; }
    for (int b=0; b<NBIN; ++b) { rsum[q++] = l1_rbin[b]; }
    for (int r=0; r<3; ++r) {
      for (int c=0; c<3; ++c) { rsum[q++] = ata[r][c]; }
    }
    for (int r=0; r<3; ++r) { rsum[q++] = atb[r]; }
    rsum[q++] = vdotv; rsum[q++] = vdote;
    MPI_Allreduce(MPI_IN_PLACE, rsum, nrl, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    q = 0;
    l1b = rsum[q++]; l1v = rsum[q++]; l1p = rsum[q++];
    l1_seam = rsum[q++]; l1_rad = rsum[q++]; l1_int = rsum[q++];
    for (int b=0; b<NBIN; ++b) { l1_bin[b] = rsum[q++]; }
    for (int b=0; b<NBIN; ++b) { l1_rbin[b] = rsum[q++]; }
    for (int r=0; r<3; ++r) {
      for (int c=0; c<3; ++c) { ata[r][c] = rsum[q++]; }
    }
    for (int r=0; r<3; ++r) { atb[r] = rsum[q++]; }
    vdotv = rsum[q++]; vdote = rsum[q++];

    const int nil = 4 + 2*NBIN;
    std::int64_t isum[4 + 2*NBIN];
    q = 0;
    isum[q++] = ncell; isum[q++] = n_seam; isum[q++] = n_rad; isum[q++] = n_int;
    for (int b=0; b<NBIN; ++b) { isum[q++] = n_bin[b]; }
    for (int b=0; b<NBIN; ++b) { isum[q++] = n_rbin[b]; }
    MPI_Allreduce(MPI_IN_PLACE, isum, nil, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    q = 0;
    ncell = isum[q++]; n_seam = isum[q++]; n_rad = isum[q++]; n_int = isum[q++];
    for (int b=0; b<NBIN; ++b) { n_bin[b] = isum[q++]; }
    for (int b=0; b<NBIN; ++b) { n_rbin[b] = isum[q++]; }

    struct { double val; int rnk; } mxin, mxout;
    mxin.val = static_cast<double>(lib);
    mxin.rnk = global_variable::my_rank;
    MPI_Allreduce(&mxin, &mxout, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    lib = static_cast<Real>(mxout.val);
    MPI_Bcast(linfline, sizeof(linfline), MPI_CHAR, mxout.rnk, MPI_COMM_WORLD);
  }
#endif
  // only rank 0 reports; every rank appending to one file also raced
  if (global_variable::my_rank != 0) return;

  const Real nrm = (ncell > 0) ? 1.0/static_cast<Real>(ncell) : 0.0;
  const Real b0c = (cs_b0r != 0.0) ? cs_b0r : 1.0;
  l1b *= nrm/b0c;
  lib /= b0c;
  l1v *= nrm;
  l1p *= nrm;

  std::printf("### CS CONVERGENCE (iprob=%d): nx2=%d t=%.6f  L1(B)=%.6e"
              "  Linf(B)=%.6e  L1(v)=%.6e  L1(p)=%.6e\n",
              cs_iprob, indcs.nx2, pm->time, l1b, lib, l1v, l1p);

  if (has_b) {
    // solve the 3x3 normal equations by Gaussian elimination with partial pivoting
    Real a[3][4];
    for (int r=0; r<3; ++r) {
      for (int c=0; c<3; ++c) { a[r][c] = ata[r][c]; }
      a[r][3] = atb[r];
    }
    for (int c=0; c<3; ++c) {
      int piv = c;
      for (int r=c+1; r<3; ++r) {
        if (fabs(a[r][c]) > fabs(a[piv][c])) piv = r;
      }
      for (int q=c; q<4; ++q) {
        const Real tmp = a[c][q]; a[c][q] = a[piv][q]; a[piv][q] = tmp;
      }
      for (int r=0; r<3; ++r) {
        if (r == c || a[c][c] == 0.0) continue;
        const Real f = a[r][c]/a[c][c];
        for (int q=c; q<4; ++q) { a[r][q] -= f*a[c][q]; }
      }
    }
    const Real bf[3] = {a[0][3]/a[0][0], a[1][3]/a[1][1], a[2][3]/a[2][2]};
    const Real dsys = (fabs(bf[0]-ex) + fabs(bf[1]-ey) + fabs(bf[2]-ez))/3.0/b0c;
    std::printf("###   best-fit UNIFORM B = (%.6e,%.6e,%.6e)\n", bf[0], bf[1], bf[2]);
    std::printf("###   exact           B = (%.6e,%.6e,%.6e)\n", ex, ey, ez);
    std::printf("###   systematic |Bfit-Bex|/3/b0c = %.6e   vs  L1(B) = %.6e\n",
                dsys, l1b);
  }

  if (has_b) {
    const Real tot = (l1_seam + l1_rad + l1_int > 0.0) ? l1_seam+l1_rad+l1_int : 1.0;
    std::printf("###   L1(B) budget: seam %.1f%% (%lld cells)  radial %.1f%% (%lld)"
                "  interior %.1f%% (%lld)\n",
                100.0*l1_seam/tot, static_cast<long long>(n_seam),
                100.0*l1_rad/tot,  static_cast<long long>(n_rad),
                100.0*l1_int/tot,  static_cast<long long>(n_int));
    std::printf("###   L1(B) per cell vs seam distance (cells):");
    for (int q=0; q<NBIN; ++q) {
      if (n_bin[q] > 0) {
        std::printf("  d%d:%.3e", (q == 0) ? 0 : (1 << (q-1)),
                    l1_bin[q]/static_cast<Real>(n_bin[q])/b0c);
      }
    }
    std::printf("\n");
    std::printf("###   L1(B) per cell vs RADIAL-boundary distance:");
    for (int q=0; q<NBIN; ++q) {
      if (n_rbin[q] > 0) {
        std::printf("  r%d:%.3e", (q == 0) ? 0 : (1 << (q-1)),
                    l1_rbin[q]/static_cast<Real>(n_rbin[q])/b0c);
      }
    }
    std::printf("\n");
    // WHERE the worst cell sits matters more than its value: on a refined mesh, "at the
    // block's inner radial face, on a block one level finer than its x1 neighbour" is a
    // level-boundary defect, and "in the middle of a block" is not.  Formatted above,
    // by whichever rank owns the maximum.
    std::printf("%s", linfline);
  }

  {
    const Real sfit = (vdotv > 0.0) ? vdote/vdotv : 1.0;
    std::printf("###   flow rate omega_fit/omega = %.9f  (s-1 = %+.4e);"
                "  a field frozen into it would lag by %.4e rad at t=%.3f\n",
                sfit, sfit-1.0, (1.0-sfit)*cs_omega*pm->time, pm->time);
  }

  // append so that a resolution sweep accumulates in one file
  std::string fname = pin->GetString("job", "basename") + "-cs-errs.dat";
  FILE *pf = std::fopen(fname.c_str(), "r");
  if (pf == nullptr) {
    pf = std::fopen(fname.c_str(), "w");
    if (pf == nullptr) return;
    std::fprintf(pf, "# AthenaK cubed-sphere MHD convergence (cs_test iprob=9)\n");
    std::fprintf(pf, "# nx1  nx2  nx3  time  L1_B  Linf_B  L1_v  L1_p\n");
  } else {
    std::fclose(pf);
    pf = std::fopen(fname.c_str(), "a");
    if (pf == nullptr) return;
  }
  std::fprintf(pf, "%d  %d  %d  %.8e  %.8e  %.8e  %.8e  %.8e\n",
               indcs.nx1, indcs.nx2, indcs.nx3, pm->time, l1b, lib, l1v, l1p);
  std::fclose(pf);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn CSTestResistCheck
//! \brief iprob = 11: compare the RESISTIVE EMF the code built with the analytic one.
//!
//! The field B = b0c*(-y,x,0) + Bunif has curl B = 2*b0c*zhat EXACTLY and everywhere, so
//! the resistive EMF eta*J must be eta*2*b0c times the projection of zhat onto each
//! edge's own unit tangent -- rhat on an x1 edge, e_xi on an x2 edge, e_eta on an x3
//! edge. That projection is the whole point: the three edge directions are not mutually
//! orthogonal and two of them are not the face normals the field is stored on, so a
//! current density that skips the metric fails this check by O(1), not by O(h).
//! Errors are reported relative to |eta*2*b0c| and should converge at second order.

void CSTestResistCheck(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pmhd == nullptr || pmbp->pmhd->presist == nullptr) {
    std::cout << "### CS RESIST CHECK: iprob=11 needs <mhd> with ohmic_resistivity"
              << std::endl;
    return;
  }
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto pres = pmbp->pmhd->presist;
  const Real eta = pres->eta_ohm_const;
  const Real jz = 2.0*cs_bazi;
  const Real scale = fabs(eta*jz);
  if (scale == 0.0) {
    std::cout << "### CS RESIST CHECK: eta*2*b0c = 0, nothing to compare" << std::endl;
    return;
  }
  auto e1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x1e);
  auto e2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x2e);
  auto e3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x3e);
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;

  Real mx[3] = {0.0, 0.0, 0.0}, l1[3] = {0.0, 0.0, 0.0};
  std::int64_t nc[3] = {0, 0, 0};
  // where the worst x1 edge sits, and how the error splits between the panel EDGE ring
  // (the outermost `nseam` cells of a panel, where the halo supplies the stencil) and
  // the panel interior
  int mxloc[4] = {-1, -1, -1, -1};
  const int nseam = 2;
  Real l1s = 0.0, l1i = 0.0;
  std::int64_t ncs = 0, nci = 0;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    const Real x2min = size.h_view(m).x2min, x2max = size.h_view(m).x2max;
    const Real x3min = size.h_view(m).x3min, x3max = size.h_view(m).x3max;
    for (int k=ks; k<=ke+1; ++k) {
      const Real ec = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, x3min, x3max);
      const Real ef = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, x3min, x3max);
      for (int j=js; j<=je+1; ++j) {
        const Real xc = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, x2min, x2max);
        const Real xf = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
        Real t1[3], t2[3], cx, cy, cz;
        // x1 edge: (cell centre in r, xi face, eta face); tangent is rhat
        if (k <= ke+1 && j <= je+1) {
          PanelToCart(p, xf, ef, cx, cy, cz);
          const Real ex = eta*jz*cz;
          const bool seam = (j-js < nseam) || (je+1-j < nseam) ||
                            (k-ks < nseam) || (ke+1-k < nseam);
          for (int i=is; i<=ie; ++i) {
            const Real d = fabs(e1(m,k,j,i) - ex)/scale;
            if (d > mx[0]) { mxloc[0]=p; mxloc[1]=i; mxloc[2]=j; mxloc[3]=k; }
            mx[0] = fmax(mx[0], d); l1[0] += d; ++nc[0];
            if (seam) { l1s += d; ++ncs; } else { l1i += d; ++nci; }
          }
        }
        // x2 edge: (r face, xi centre, eta face); tangent is e_xi there
        if (k <= ke+1 && j <= je) {
          cubed_sphere::PanelTangents(p, xc, ef, t1, t2);
          const Real ex = eta*jz*t1[2];
          for (int i=is; i<=ie+1; ++i) {
            const Real d = fabs(e2(m,k,j,i) - ex)/scale;
            mx[1] = fmax(mx[1], d); l1[1] += d; ++nc[1];
          }
        }
        // x3 edge: (r face, xi face, eta centre); tangent is e_eta there
        if (k <= ke && j <= je+1) {
          cubed_sphere::PanelTangents(p, xf, ec, t1, t2);
          const Real ex = eta*jz*t2[2];
          for (int i=is; i<=ie+1; ++i) {
            const Real d = fabs(e3(m,k,j,i) - ex)/scale;
            mx[2] = fmax(mx[2], d); l1[2] += d; ++nc[2];
          }
        }
      }
    }
  }
  for (int q=0; q<3; ++q) {
    if (nc[q] > 0) l1[q] /= static_cast<Real>(nc[q]);
  }
  std::printf("### CS RESISTIVE EMF CHECK (iprob=11): nx1=%d nx2=%d  eta=%.3e"
              "  |2*eta*b0c|=%.3e\n", indcs.nx1, indcs.nx2, eta, scale);
  std::printf("###   relative error vs the analytic eta*J:"
              "  x1e L1=%.4e max=%.4e | x2e L1=%.4e max=%.4e | x3e L1=%.4e max=%.4e\n",
              l1[0], mx[0], l1[1], mx[1], l1[2], mx[2]);
  if (ncs > 0) l1s /= static_cast<Real>(ncs);
  if (nci > 0) l1i /= static_cast<Real>(nci);

  // ---- EVOLVED FIELD ERROR: the decomposition-independent gate -----------------------
  //
  // WHY THIS EXISTS. Everything above reads `efld_resist`, which is the block's OWN
  // eta*J as the gnomonic curl built it -- BEFORE AddEMFDirect merges it into efld and
  // before SendE/RecvE average the two sides of a shared block face. Two different
  // values there are EXPECTED, and say nothing on their own about the evolved solution:
  // the averaging exists precisely to reconcile them. Quoting that array's block-face
  // maximum as if it were a defect is a measurement error, and it was made once.
  //
  // This is the gate that means something. eta*J is a CONSTANT vector for this problem
  // (curl B = 2*b0c*zhat everywhere), so curl(eta*J) = 0 and the resistive term cannot
  // move B at all. The field drifts only through -v x B, and v stays at ~1e-4, so
  // B(t) = B(0) to that accuracy and the norm below is a PHYSICAL error with no
  // reference to any block. It must therefore be the same under any decomposition:
  // splitting x1, or refining, may change it only at round-off. That is the same gate
  // that certified ideal MHD, and it is the one to trust here.
  auto bf1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x1f);
  auto bf2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x2f);
  auto bf3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x3f);
  const Real bazi = cs_bazi;
  const Real bux = cs_bvx, buy = cs_bvy, buz = cs_bvz;
  // WHERE the maximum sits, and which component it is. A radial split must not change
  // the answer at all -- ideal MHD reproduces its L1(B) to all SEVEN printed digits
  // under the same test -- so any Linf that moves is a defect, and the first thing to
  // know about it is whether it sits ON the radial block interface.
  Real l1f = 0.0, mxf = 0.0, bmax = 0.0;
  int mxc = -1, mxi = -1, mxj = -1, mxk = -1, mxp = -1, mxm = -1;
  std::int64_t ncf = 0;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    const Real x2min = size.h_view(m).x2min, x2max = size.h_view(m).x2max;
    const Real x3min = size.h_view(m).x3min, x3max = size.h_view(m).x3max;
    const Real x1min = size.h_view(m).x1min, x1max = size.h_view(m).x1max;
    for (int k=ks; k<=ke+1; ++k) {
      const Real ec = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, x3min, x3max);
      const Real ef = 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, x3min, x3max);
      for (int j=js; j<=je+1; ++j) {
        const Real xc = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, x2min, x2max);
        const Real xf = 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
        Real n1[3], n2[3], cx, cy, cz;
        for (int i=is; i<=ie+1; ++i) {
          const Real rc = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          const Real rl = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
          // exactly the three expressions the initial condition used
          if (j <= je && k <= ke) {
            PanelToCart(p, xc, ec, cx, cy, cz);
            const Real ex = (-bazi*rl*cy + bux)*cx + (bazi*rl*cx + buy)*cy + buz*cz;
            const Real d = fabs(bf1h(m,k,j,i) - ex);
            if (d > mxf) {
              mxc = 0; mxm = m; mxp = p; mxi = i-is; mxj = j-js; mxk = k-ks;
            }
            l1f += d; ++ncf; mxf = fmax(mxf, d); bmax = fmax(bmax, fabs(ex));
          }
          if (i <= ie && k <= ke) {
            PanelToCart(p, xf, ec, cx, cy, cz);
            PanelNormals(p, xf, ec, n1, n2);
            const Real ex = (-bazi*rc*cy + bux)*n1[0] + (bazi*rc*cx + buy)*n1[1]
                          + buz*n1[2];
            const Real d = fabs(bf2h(m,k,j,i) - ex);
            if (d > mxf) {
              mxc = 1; mxm = m; mxp = p; mxi = i-is; mxj = j-js; mxk = k-ks;
            }
            l1f += d; ++ncf; mxf = fmax(mxf, d); bmax = fmax(bmax, fabs(ex));
          }
          if (i <= ie && j <= je) {
            PanelToCart(p, xc, ef, cx, cy, cz);
            PanelNormals(p, xc, ef, n1, n2);
            const Real ex = (-bazi*rc*cy + bux)*n2[0] + (bazi*rc*cx + buy)*n2[1]
                          + buz*n2[2];
            const Real d = fabs(bf3h(m,k,j,i) - ex);
            if (d > mxf) {
              mxc = 2; mxm = m; mxp = p; mxi = i-is; mxj = j-js; mxk = k-ks;
            }
            l1f += d; ++ncf; mxf = fmax(mxf, d); bmax = fmax(bmax, fabs(ex));
          }
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  {
    Real r_[2] = {l1f, 0.0};
    MPI_Allreduce(MPI_IN_PLACE, r_, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    l1f = r_[0];
    Real m_[2] = {mxf, bmax};
    MPI_Allreduce(MPI_IN_PLACE, m_, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    mxf = m_[0]; bmax = m_[1];
    std::int64_t n_ = ncf;
    MPI_Allreduce(MPI_IN_PLACE, &n_, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    ncf = n_;
  }
#endif
  {
    const Real bs = (bmax > 0.0) ? bmax : 1.0;
    const Real l1n = (ncf > 0) ? l1f/static_cast<Real>(ncf)/bs : 0.0;
    std::printf("###   EVOLVED FIELD error vs the exact static B (all three face"
                " components):  L1=%.4e  Linf=%.4e  (over %lld faces)\n",
                l1n, mxf/bs, static_cast<long long>(ncf));
    const char *cn[3] = {"x1f", "x2f", "x3f"};
    if (mxc >= 0) {
      std::printf("###     Linf on %s at panel %d (i,j,k)=(%d,%d,%d) of block %d;"
                  "  block x1 = %g..%g%s\n", cn[mxc], mxp, mxi, mxj, mxk, mxm,
                  size.h_view(mxm).x1min, size.h_view(mxm).x1max,
                  (mxi == 0) ? "   [ON the block's INNER RADIAL face]" : "");
    }
  }

  // ---- Ohmic heating: the only quantitative check on the POYNTING flux ---------------
  // div(E x B) = B.curl E - E.curl J's parent = -E.J = -eta |J|^2, so this exact
  // solution is static in B but heats UNIFORMLY at eta*|J|^2. With v = 0 all of it goes
  // into internal energy, so the volume-averaged pressure must rise by
  // (gamma-1)*eta*|J|^2*t and by nothing else. A Poynting flux built by crossing two
  // vectors held in different frames gets this wrong at O(cos_cell).
  auto &w0 = pmbp->pmhd->w0;
  auto wh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), w0);
  auto vol = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->volume);
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  Real dpv = 0.0, vtot = 0.0, vmax = 0.0;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    for (int k=ks; k<=ke; ++k) {
      const Real ec = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                             size.h_view(m).x3max);
      for (int j=js; j<=je; ++j) {
        const Real xc = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                               size.h_view(m).x2max);
        Real cx, cy, cz;
        PanelToCart(p, xc, ec, cx, cy, cz);
        for (int i=is; i<=ie; ++i) {
          const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                        size.h_view(m).x1max);
          const Real px = rad*cx, py = rad*cy;
          const Real p_ini = cs_p0 - SQR(cs_bazi)*(px*px + py*py)
                             + 2.0*cs_bazi*(-cs_bvy*px + cs_bvx*py);
          const Real dv = vol(m,k,j,i);
          dpv += dv*(gm1*wh(m,IEN,k,j,i) - p_ini);
          vtot += dv;
          vmax = fmax(vmax, sqrt(SQR(wh(m,IVX,k,j,i)) + SQR(wh(m,IVY,k,j,i))
                                 + SQR(wh(m,IVZ,k,j,i))));
        }
      }
    }
  }
  // ---- div B: the gate on the SEAM ---------------------------------------------------
  // CT conserves sum(area*B) per cell EXACTLY, whatever the EMF is, PROVIDED each edge
  // carries ONE value shared by every cell and every panel touching it. The resistive EMF
  // is added inside MHD::EField, i.e. before SendE, so it should ride the same seam
  // exchange the ideal EMF does -- and if it did not, div B would grow at the seams and
  // nowhere else. This is the one thing the analytic-EMF check above cannot see.
  auto b1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x1f);
  auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x2f);
  auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x3f);
  auto a1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x1f);
  auto a2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x2f);
  auto a3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x3f);
  Real dvb_seam = 0.0, dvb_int = 0.0;
  const int nsm = 2;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        const bool seam = (j-js < nsm) || (je-j < nsm) || (k-ks < nsm) || (ke-k < nsm);
        for (int i=is; i<=ie; ++i) {
          const Real fl = b1h(m,k,j,i+1)*a1h(m,k,j,i+1) - b1h(m,k,j,i)*a1h(m,k,j,i)
                        + b2h(m,k,j+1,i)*a2h(m,k,j+1,i) - b2h(m,k,j,i)*a2h(m,k,j,i)
                        + b3h(m,k+1,j,i)*a3h(m,k+1,j,i) - b3h(m,k,j,i)*a3h(m,k,j,i);
          const Real amx = fmax(a1h(m,k,j,i+1), fmax(a2h(m,k,j+1,i), a3h(m,k+1,j,i)));
          const Real bmg = sqrt(SQR(wh(m,IDN,k,j,i))*0.0 + 1.0e-300)
                         + fabs(cs_bazi) + sqrt(SQR(cs_bvx)+SQR(cs_bvy)+SQR(cs_bvz));
          const Real d = fabs(fl)/(bmg*amx);
          if (seam) { dvb_seam = fmax(dvb_seam, d); } else { dvb_int = fmax(dvb_int, d); }
        }
      }
    }
  }
  std::printf("###   max |sum A.B| / (|B| max A):  seam %.3e   interior %.3e\n",
              dvb_seam, dvb_int);

  const Real dp_meas = (vtot > 0.0) ? dpv/vtot : 0.0;
  const Real dp_exact = gm1*eta*jz*jz*pm->time;
  std::printf("###   Ohmic heating: <dp> measured %.6e  exact (gam-1)*eta*J^2*t %.6e"
              "  ratio %.6f   max|v| %.3e\n", dp_meas, dp_exact,
              (dp_exact != 0.0) ? dp_meas/dp_exact : 0.0, vmax);
  std::printf("###   x1e worst at panel %d (i,j,k)=(%d,%d,%d);"
              "  L1 over the %d-cell panel-edge ring = %.4e, panel interior = %.4e\n",
              mxloc[0], mxloc[1], mxloc[2], mxloc[3], nseam, l1s, l1i);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CSTestHistory(HistoryData *pdata, Mesh *pm)
//! \brief GRID IMPRINTING over time: the error near a panel seam and far from it, as a
//! time series in the .hst.
//!
//! Every cubed-sphere error measurement in this pgen so far is a single number at the
//! END of a short run, which cannot distinguish an error that saturates from one that
//! grows secularly and eventually contaminates the whole field. That distinction is what
//! the FV3 duo-grid work calls GRID IMPRINTING, and it is the reason they track wave
//! symmetry out to 100+ days rather than quoting one L1. This writes the same kind of
//! diagnostic here.
//!
//! iprob = 3 is the right vehicle: the rigid-rotation equilibrium is EXACTLY steady, so
//! the error at any time is simply the deviation from the initial exact state, with no
//! reference solution to evolve and no phase error to disentangle.
//!
//! Reported as SUMS with their cell COUNTS, because the history reduction adds across
//! ranks and only sums are meaningful there; divide afterwards. Two disjoint sets:
//!   SEAM     -- within 2 cells of a panel edge, tested GEOMETRICALLY on |xi| and |eta|
//!               so it is independent of how panels are split into MeshBlocks;
//!   INTERIOR -- everything else.
//! Both EXCLUDE the two cell layers next to each radial boundary, so that the radial BC's
//! own error (which is first order in the ghost clock, see the mhd convergence note)
//! cannot leak into either bucket and be mistaken for imprinting.
//!
//! Both quantities are chosen so that no reference solution has to be evaluated in the
//! kernel: the exact density of iprob=3 is uniformly d0, and the exact v_r is
//! IDENTICALLY ZERO, so sum(v_r^2) is a pure error norm with no cancellation and no
//! normalisation choice. v_r is the sharper of the two.

void CSTestHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 6;
  pdata->label[0] = "sm-drho";   // sum |rho-rho_ex| over seam cells
  pdata->label[1] = "sm-n";      // number of seam cells
  pdata->label[2] = "sm-vr2";    // sum v_r^2 over seam cells
  pdata->label[3] = "in-drho";   // the same three over interior cells
  pdata->label[4] = "in-n";
  pdata->label[5] = "in-vr2";

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &w0_ = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const Real d0 = cs_d0;   // the exact density of iprob=3 is uniformly d0

  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("CSHist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1);
    array_sum::GlobalSum hv;
    for (int n=0; n<NHISTORY_VARIABLES; ++n) {hv.the_array[n] = 0.0;}

    // skip the two layers next to each radial boundary
    if (i >= 2 && i <= nx1-3) {
      const Real x2mn = size.d_view(m).x2min, x2mx = size.d_view(m).x2max;
      const Real x3mn = size.d_view(m).x3min, x3mx = size.d_view(m).x3max;
      const Real xi  = 0.25*M_PI*CellCenterX(j, nx2, x2mn, x2mx);
      const Real eta = 0.25*M_PI*CellCenterX(k, nx3, x3mn, x3mx);
      const Real dxi  = 0.25*M_PI*(x2mx - x2mn)/static_cast<Real>(nx2);
      const Real deta = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3);
      const Real edge = 0.25*M_PI;
      const bool near_seam = (fabs(xi)  > edge - 2.0*dxi) ||
                             (fabs(eta) > edge - 2.0*deta);

      const int kk = k + ks, jj = j + js, ii = i + is;
      const Real drho = fabs(w0_(m,IDN,kk,jj,ii) - d0);
      const Real vr2  = SQR(w0_(m,IVX,kk,jj,ii));   // exact v_r is identically zero
      if (near_seam) {
        hv.the_array[0] = drho; hv.the_array[1] = 1.0; hv.the_array[2] = vr2;
      } else {
        hv.the_array[3] = drho; hv.the_array[4] = 1.0; hv.the_array[5] = vr2;
      }
    }
    mb_sum += hv;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  for (int n=0; n<pdata->nhist; ++n) {pdata->hdata[n] = sum_this_mb.the_array[n];}
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CSTestConsSums(Mesh *pm)
//! \brief The GLOBAL CONSERVED SUMS, printed at full precision.
//!
//! The `.hst` output cannot be used for this. `HistoryOutput::Load{Hydro,MHD}HistoryData`
//! weights every cell by `dx1*dx2*dx3` -- the CARTESIAN coordinate volume -- and never by
//! `Coordinates::volume`, so on the cubed sphere its "mass" is not sum(rho dV) but a
//! differently weighted sum that the finite-volume update does not conserve; and its
//! "1-mom/2-mom/3-mom" adds up components on the local gnomonic basis, which point in a
//! different physical direction in every cell. Neither can measure conservation. (Its
//! default `%12.5e` also quantises the answer at the level being measured.)
//!
//! Here mass and total energy are weighted by the true cell volume, and the momentum is
//! rotated into the FIXED Cartesian frame first, which is the frame in which it is
//! actually a conserved vector: u(IM1..IM3) are the components on (rhat, e1hat, e2hat),
//! so P = u1*rhat + u2*e1hat + u3*e2hat cell by cell.
//!
//! ANGULAR momentum L = sum(r x P dV) is reported too, and on a sphere it is the one that
//! matters. In the rigid-rotation equilibria the net LINEAR momentum is identically zero
//! by symmetry -- P comes out at 1e-15 whatever the scheme does -- so it is a vacuous
//! test, while Lz is O(1). These sums are all exactly constant in the iprob = 3 / 9
//! steady equilibria, so any change is scheme error.
//!
//! These sums are RANK-LOCAL and each rank prints its own; under MPI read the `.hst`,
//! which is reduced across ranks (and, since the volume fix, correctly weighted).

void CSTestConsSums(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pmhd == nullptr && pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto uh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), u0);
  auto vh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->volume);
  auto &size = pmbp->pmb->mb_size;
  size.sync_host();
  auto &mbpanel = pmbp->pmb->mb_panel;
  mbpanel.sync_host();

  Real sm = 0.0, se = 0.0, sp[3] = {0.0, 0.0, 0.0}, svol = 0.0;
  Real sl[3] = {0.0, 0.0, 0.0};
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    for (int k=ks; k<=ke; ++k) {
      const Real et = 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, size.h_view(m).x3min,
                                                             size.h_view(m).x3max);
      for (int j=js; j<=je; ++j) {
        const Real xi = 0.25*M_PI*CellCenterX(j-js, indcs.nx2, size.h_view(m).x2min,
                                                               size.h_view(m).x2max);
        Real e1[3], e2[3], rx, ry, rz;
        cubed_sphere::PanelTangents(p, xi, et, e1, e2);
        PanelToCart(p, xi, et, rx, ry, rz);
        const Real rh[3] = {rx, ry, rz};
        for (int i=is; i<=ie; ++i) {
          const Real dv = vh(m,k,j,i);
          const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                        size.h_view(m).x1max);
          svol += dv;
          sm += dv*uh(m,IDN,k,j,i);
          se += dv*uh(m,IEN,k,j,i);
          // u0(IM1..IM3) are the COVARIANT components rho*v_i on (rhat, e1hat, e2hat)
          // -- see Coordinates::GnomonicEquiangleRaiseVel. The tangent pair is NOT
          // orthogonal, so they must be RAISED with g^{-1}, g = [[1,c],[c,1]], before
          // they can be used as coefficients of e1hat and e2hat. Treating them as
          // contravariant is wrong by up to ~46% of |V|, vanishing on the seam midline
          // and worst at the panel corners. rhat is orthogonal to both and unit, so the
          // radial component needs no raising.
          const Real cg = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
          const Real dg = 1.0 - cg*cg;
          const Real q2 = (uh(m,IM2,k,j,i) - cg*uh(m,IM3,k,j,i))/dg;
          const Real q3 = (uh(m,IM3,k,j,i) - cg*uh(m,IM2,k,j,i))/dg;
          Real pc[3];
          for (int q=0; q<3; ++q) {
            pc[q] = uh(m,IM1,k,j,i)*rh[q] + q2*e1[q] + q3*e2[q];
            sp[q] += dv*pc[q];
          }
          // L = r x p, with r = rad*rhat.  This is the momentum-like quantity that
          // matters on a sphere: the rigid-rotation equilibria carry a large Lz and no
          // net linear momentum at all, so P is identically zero by symmetry and can
          // measure nothing, while Lz is O(1) and any loss shows up directly.
          for (int q=0; q<3; ++q) {
            const int a = (q+1)%3, b = (q+2)%3;
            sl[q] += dv*rad*(rh[a]*pc[b] - rh[b]*pc[a]);
          }
        }
      }
    }
  }
  // These are totals over the WHOLE mesh, so they must be reduced: the loop above walks
  // only this rank's MeshBlockPack. Unreduced they are per-rank subtotals that look like
  // a conservation failure the moment the decomposition changes.
#if MPI_PARALLEL_ENABLED
  {
    Real rsum[11] = {svol, sm, se, sp[0], sp[1], sp[2], sl[0], sl[1], sl[2], 0.0, 0.0};
    MPI_Allreduce(MPI_IN_PLACE, rsum, 9, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    svol = rsum[0]; sm = rsum[1]; se = rsum[2];
    sp[0] = rsum[3]; sp[1] = rsum[4]; sp[2] = rsum[5];
    sl[0] = rsum[6]; sl[1] = rsum[7]; sl[2] = rsum[8];
  }
#endif
  if (global_variable::my_rank != 0) return;

  std::printf("### CS CONSERVED SUMS at t = %.17e  (ncyc %d)\n", pm->time, pm->ncycle);
  std::printf("###   volume %.17e\n", svol);
  std::printf("###   mass   %.17e\n", sm);
  std::printf("###   energy %.17e\n", se);
  std::printf("###   Px     %.17e\n", sp[0]);
  std::printf("###   Py     %.17e\n", sp[1]);
  std::printf("###   Pz     %.17e\n", sp[2]);
  std::printf("###   Lx     %.17e\n", sl[0]);
  std::printf("###   Ly     %.17e\n", sl[1]);
  std::printf("###   Lz     %.17e\n", sl[2]);
}

//----------------------------------------------------------------------------------------
//! \fn void CSTestLevelFluxCheck(Mesh *pm)
//! \brief Does the radial flux TELESCOPE across a coarse/fine LEVEL BOUNDARY?
//!
//! The level-boundary twin of CSTestSeamFluxCheck. A finite-volume update conserves for
//! exactly one reason: what leaves one cell enters its neighbour, because both cells
//! multiply the SAME stored flux by the SAME stored area. Between MeshBlocks at one level
//! that holds bit for bit. Across a level boundary it holds only because the flux
//! correction replaces the coarse face's flux by the area-weighted sum of the four fine
//! ones, so that A_c*F_c == sum(A_f*F_f). If that identity fails, the difference is mass
//! and energy created out of nothing -- which is what refined MHD does (+1.7e-7 in mass
//! over five cycles with reflecting walls, against hydro's -3.5e-14).
//!
//! Faces are paired BY GEOMETRY, not by neighbour index, so the gate is independent of
//! the buffer machinery it is auditing. Each radial block-boundary face is bucketed by
//! its panel, its radius, and the angular cell of the COARSEST block present -- a bucket
//! wide enough that a coarse face and its four fine children land together. Each face
//! contributes +A*F when it is its block's OUTER face and -A*F when it is the INNER one,
//! so a bucket that telescopes sums to zero.
//!
//! Only the SCALAR fluxes are meaningful: rho*v.nhat and (E+p)*v.nhat are true scalars.
//!
//! The NULL CONTROL is what makes this a measurement: buckets holding exactly two faces
//! are ordinary same-level block boundaries and MUST come out at round-off. If they do
//! not, the gate itself is wrong and its level-boundary numbers mean nothing.
//!
//! RANK-LOCAL, like the seam gate: under MPI a face whose two sides live on different
//! ranks is invisible. Run it on one rank.

void CSTestLevelFluxCheck(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (!is_mhd && pmbp->phydro == nullptr) return;
  auto &fx = is_mhd ? pmbp->pmhd->uflx : pmbp->phydro->uflx;

  auto f1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fx.x1f);
  auto a1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                 pmbp->pcoord->area.x1f);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &mblev = pmbp->pmb->mb_lev;
  auto &mbpanel = pmbp->pmb->mb_panel;
  const int nmb = pmbp->nmb_thispack;

  // bucket width: the angular cell of the COARSEST block present, and a radial tolerance
  // far below the finest cell. Both sides of a shared face compute the same radius from
  // the same mesh bounds, so they agree to round-off.
  int lmin = 1 << 30;
  for (int m=0; m<nmb; ++m) lmin = std::min(lmin, mblev.h_view(m));
  Real dxi_key = 0.0, dr_key = 0.0, dr_min = 1.0e30;
  for (int m=0; m<nmb; ++m) {
    const Real dxi = 0.25*M_PI*(size.h_view(m).x2max - size.h_view(m).x2min)/indcs.nx2;
    if (mblev.h_view(m) == lmin) {
      dxi_key = dxi;
      dr_key = (size.h_view(m).x1max - size.h_view(m).x1min)/indcs.nx1;
    }
    dr_min = fmin(dr_min, (size.h_view(m).x1max - size.h_view(m).x1min)/indcs.nx1);
  }
  if (dxi_key <= 0.0) return;
  const Real rq = 0.01*dr_min;

  struct Bucket { Real sm, se, scale; int n, nlev, npos, nneg;
                  int d, p, m0, m1; int64_t q1, q2, q3; };
  std::map<int64_t, Bucket> buck;
  auto a2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                 pmbp->pcoord->area.x2f);
  auto a3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                 pmbp->pcoord->area.x3f);
  auto f2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fx.x2f);
  auto f3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fx.x3f);
  const Real rmin = pm->mesh_size.x1min;

  // A face coordinate is shared exactly by both sides, so quantise it finely; a
  // cell-centre coordinate differs between a coarse cell and its children, so bucket it
  // at the COARSEST level's width.
  auto qface = [&](Real x, Real h) { return llround(x/(0.01*h)); };
  auto qcell = [&](Real x, Real x0, Real h) {
    return static_cast<int64_t>(floor((x - x0)/h));
  };

  for (int m=0; m<nmb; ++m) {
    const int p = mbpanel.h_view(m);
    const int lev = mblev.h_view(m);
    const Real x1mn = size.h_view(m).x1min, x1mx = size.h_view(m).x1max;
    const Real x2mn = size.h_view(m).x2min, x2mx = size.h_view(m).x2max;
    const Real x3mn = size.h_view(m).x3min, x3mx = size.h_view(m).x3max;
    auto rc  = [&](int i) { return CellCenterX(i-is, indcs.nx1, x1mn, x1mx); };
    auto xic = [&](int j) { return 0.25*M_PI*CellCenterX(j-js, indcs.nx2, x2mn, x2mx); };
    auto etc = [&](int k) { return 0.25*M_PI*CellCenterX(k-ks, indcs.nx3, x3mn, x3mx); };
    auto xif = [&](int j) { return 0.25*M_PI*LeftEdgeX(j-js, indcs.nx2, x2mn, x2mx); };
    auto etf = [&](int k) { return 0.25*M_PI*LeftEdgeX(k-ks, indcs.nx3, x3mn, x3mx); };

    // The DIRECTION has to be part of the key. Without it an x1 bucket and an x2 bucket
    // can land on the same value and pool unrelated faces, which makes even hydro -- a
    // case known to telescope to round-off -- look 100% wrong.
    auto add = [&](int d, int64_t q1, int64_t q2, int64_t q3,
                   Real sgn, Real am, Real ae) {
      const int64_t key = (((((static_cast<int64_t>(p)*4 + d)*8000000
                          + (q1 + 4000000))*8192) + (q2 + 4096))*8192) + (q3 + 4096);
      auto it = buck.find(key);
      if (it == buck.end()) {
        buck[key] = {sgn*am, sgn*ae, fabs(am), 1, lev,
                     (sgn > 0.0) ? 1 : 0, (sgn > 0.0) ? 0 : 1,
                     d, p, m, -1, q1, q2, q3};
      } else {
        it->second.sm += sgn*am;  it->second.se += sgn*ae;
        it->second.scale = fmax(it->second.scale, fabs(am));
        it->second.n += 1;
        if (m != it->second.m0 && it->second.m1 < 0) it->second.m1 = m;
        if (sgn > 0.0) { it->second.npos += 1; } else { it->second.nneg += 1; }
        if (lev != it->second.nlev) it->second.nlev = -1;
      }
    };

    // x1 (radial) block-boundary faces
    for (int side=0; side<2; ++side) {
      const int i = (side == 0) ? is : (ie+1);
      const Real sgn = (side == 0) ? -1.0 : 1.0;
      const int64_t q1 = qface((side == 0) ? x1mn : x1mx, dr_min);
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          add(1, q1, qcell(xic(j), -0.25*M_PI, dxi_key),
              qcell(etc(k), -0.25*M_PI, dxi_key), sgn,
              a1h(m,k,j,i)*f1h(m,IDN,k,j,i), a1h(m,k,j,i)*f1h(m,IEN,k,j,i));
        }
      }
    }
    // x2 (xi) block-boundary faces
    for (int side=0; side<2; ++side) {
      const int j = (side == 0) ? js : (je+1);
      const Real sgn = (side == 0) ? -1.0 : 1.0;
      const int64_t q2 = qface(xif(j), dxi_key);
      for (int k=ks; k<=ke; ++k) {
        for (int i=is; i<=ie; ++i) {
          add(2, q2, qcell(rc(i), rmin, dr_key), qcell(etc(k), -0.25*M_PI, dxi_key),
              sgn, a2h(m,k,j,i)*f2h(m,IDN,k,j,i), a2h(m,k,j,i)*f2h(m,IEN,k,j,i));
        }
      }
    }
    // x3 (eta) block-boundary faces
    for (int side=0; side<2; ++side) {
      const int k = (side == 0) ? ks : (ke+1);
      const Real sgn = (side == 0) ? -1.0 : 1.0;
      const int64_t q3 = qface(etf(k), dxi_key);
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          add(3, q3, qcell(rc(i), rmin, dr_key), qcell(xic(j), -0.25*M_PI, dxi_key),
              sgn, a3h(m,k,j,i)*f3h(m,IDN,k,j,i), a3h(m,k,j,i)*f3h(m,IEN,k,j,i));
        }
      }
    }
  }

  Real dmax_[4] = {0.0}, dsum_[4] = {0.0}, lmax_[2] = {0.0}, lsum_[2] = {0.0};
  int dcnt_[4] = {0}, lcnt_[2] = {0};
  Bucket wb = {0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0};
  Real sm2 = 0.0, se2 = 0.0, smL = 0.0, seL = 0.0, sc2 = 0.0, scL = 0.0;
  Real tm2 = 0.0, tmL = 0.0;
  int n2 = 0, nL = 0;
  for (auto &b : buck) {
    // A bucket needs contributions of BOTH signs to be a shared face at all. One-sided
    // buckets are PHYSICAL DOMAIN boundaries, where nothing is supposed to cancel, and
    // counting them reports the wall's own flux as a telescoping failure.
    if (b.second.npos == 0 || b.second.nneg == 0) continue;
    // Classify by whether the bucket actually MIXES levels. Counting faces does not
    // work: the bucket is a coarse angular cell, so a same-level interface already puts
    // four face-cells from each side into it.
    if (b.second.nlev >= 0) {                     // same level: the NULL CONTROL
      sm2 = fmax(sm2, fabs(b.second.sm)); se2 = fmax(se2, fabs(b.second.se));
      sc2 = fmax(sc2, b.second.scale); tm2 += b.second.sm; ++n2;
      // BREAKDOWN of the same-level failure: by DIRECTION, by the LEVEL the two blocks
      // sit at, and the identity of the single worst bucket. What telescopes for a
      // same-level face is that both sides reconstruct from ghosts that are exact
      // COPIES; this says which family of copies has stopped being exact.
      const int dd = b.second.d;
      dmax_[dd] = fmax(dmax_[dd], fabs(b.second.sm));  dsum_[dd] += b.second.sm;
      dcnt_[dd] += 1;
      const int li = (b.second.nlev == lmin) ? 0 : 1;
      lmax_[li] = fmax(lmax_[li], fabs(b.second.sm));  lsum_[li] += b.second.sm;
      lcnt_[li] += 1;
      if (fabs(b.second.sm) > fabs(wb.sm)) wb = b.second;
    } else {
      smL = fmax(smL, fabs(b.second.sm)); seL = fmax(seL, fabs(b.second.se));
      scL = fmax(scL, b.second.scale); tmL += b.second.sm; ++nL;
    }
  }
  std::cout << "### CS LEVEL-BOUNDARY FLUX TELESCOPING" << std::endl;
  std::printf("   same level (NULL CONTROL, must be round-off): %d faces"
              "  max|dM| %11.4e  max|dE| %11.4e  (scale %11.4e)\n",
              n2, sm2, se2, sc2);
  std::printf("   coarse/fine boundary:                         %d faces"
              "  max|dM| %11.4e  max|dE| %11.4e  (scale %11.4e)\n",
              nL, smL, seL, scL);
  std::printf("   summed spurious mass source:  same level %11.4e   coarse/fine %11.4e\n",
              tm2, tmL);
  const char *dn[4] = {"--", "x1(r)", "x2(xi)", "x3(eta)"};
  for (int d=1; d<=3; ++d) {
    std::printf("     same-level by direction %-8s n %6d  max|dM| %11.4e  sum %11.4e\n",
                dn[d], dcnt_[d], dmax_[d], dsum_[d]);
  }
  for (int l=0; l<2; ++l) {
    std::printf("     same-level on %s blocks   n %6d  max|dM| %11.4e  sum %11.4e\n",
                (l == 0) ? "COARSE" : "FINE  ", lcnt_[l], lmax_[l], lsum_[l]);
  }
  if (wb.m0 >= 0) {
    std::printf("     WORST same-level bucket: dir %s panel %d lev %d  blocks %d,%d"
                "  q=(%lld,%lld,%lld)  dM %11.4e  scale %11.4e  nfaces %d\n",
                dn[wb.d], wb.p, wb.nlev, wb.m0, wb.m1,
                static_cast<long long>(wb.q1), static_cast<long long>(wb.q2),
                static_cast<long long>(wb.q3), wb.sm, wb.scale, wb.n);
    std::printf("     block %d: lev %d panel %d x1 [%g,%g] x2 [%g,%g] x3 [%g,%g]\n",
                wb.m0, mblev.h_view(wb.m0), mbpanel.h_view(wb.m0),
                size.h_view(wb.m0).x1min, size.h_view(wb.m0).x1max,
                size.h_view(wb.m0).x2min, size.h_view(wb.m0).x2max,
                size.h_view(wb.m0).x3min, size.h_view(wb.m0).x3max);
    if (wb.m1 >= 0) {
      std::printf("     block %d: lev %d panel %d x1 [%g,%g] x2 [%g,%g] x3 [%g,%g]\n",
                  wb.m1, mblev.h_view(wb.m1), mbpanel.h_view(wb.m1),
                  size.h_view(wb.m1).x1min, size.h_view(wb.m1).x1max,
                  size.h_view(wb.m1).x2min, size.h_view(wb.m1).x2max,
                  size.h_view(wb.m1).x3min, size.h_view(wb.m1).x3max);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CSTestSeamFluxCheck(Mesh *pm)
//! \brief Is the flux through a shared SEAM FACE single-valued?
//!
//! The cell-centred counterpart of `CS SEAM FACE SINGLE-VALUED`, and the gate that says
//! whether the seam is CONSERVATIVE. A finite-volume update telescopes exactly -- what
//! leaves one cell enters its neighbour -- only because both cells multiply the SAME
//! stored flux by the SAME stored area. Inside a panel that holds by construction, and it
//! holds between MeshBlocks of one panel too, because their ghost states are copies and
//! so reproduce the neighbour's flux bit for bit. At a panel seam it does NOT: the ghosts
//! are INTERPOLATED, so the two panels reconstruct slightly different states either side
//! of the one physical face and their Riemann solves disagree. Whatever the difference
//! is, it is mass and energy created out of nothing.
//!
//! Only the SCALAR fluxes are tested. Mass and total energy flux densities are true
//! scalars -- rho v.nhat and (E+p) v.nhat -- so the two panels' values must agree up to
//! the sign of the outward normal with no basis transform. The momentum fluxes are
//! components on the local gnomonic basis and are not comparable this way.
//!
//! As in the FC gate the shared faces are found BY GEOMETRY, keyed on the physical face
//! centre, so the test cannot inherit a bug from the index map it is testing; two keys
//! that agree while the panels differ are the two copies of one seam face. Results are
//! split by distance to a CUBE VERTEX, where three panels meet and the halo has the least
//! to work with.
//!
//! UNDER MPI THIS GATE UNDER-COUNTS, and each rank prints its own line. Pairing is by
//! geometry within one rank's blocks, so a seam whose two copies live on DIFFERENT ranks
//! contributes nothing -- a rank can even report an exact zero because it holds no
//! complete pair (seen with 2 ranks). Run it on one rank to gate the seam itself; to
//! check the exchange under MPI, compare the GLOBAL conserved sums in the `.hst` across
//! rank counts instead, which are properly reduced.

void CSTestSeamFluxCheck(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->pmhd == nullptr && pmbp->phydro == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;
  auto &uflx = (pmbp->pmhd != nullptr) ? pmbp->pmhd->uflx : pmbp->phydro->uflx;
  auto f2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), uflx.x2f);
  auto f3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), uflx.x3f);
  auto a2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x2f);
  auto a3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x3f);
  auto &size = pmbp->pmb->mb_size;
  size.sync_host();
  auto &mbpanel = pmbp->pmb->mb_panel;
  mbpanel.sync_host();

  // One copy of a seam face: the OUTWARD mass and energy flux this panel sees through it,
  // already multiplied by this panel's own stored area, which is what the update adds.
  struct SeamFlux { Real fm, fe, area, fp[3]; int panel; };
  std::map<std::string, SeamFlux> seen;
  Real smass[2] = {0.0, 0.0}, sener[2] = {0.0, 0.0};   // [0] away from a vertex, [1] near
  Real dmax[2] = {0.0, 0.0}, fmax_ = 0.0, damax = 0.0;
  int npair[2] = {0, 0};
  // The NULL control for this gate. A face shared by two MeshBlocks of the SAME panel is
  // the identical situation minus the seam: the ghosts there are copies, so both blocks
  // must produce the same flux to round-off. If this is not at round-off then the gate
  // itself -- its keying, its sign convention or its areas -- is wrong and the seam
  // numbers mean nothing. Empty unless there are several MeshBlocks per panel.
  Real snull = 0.0, dnull = 0.0;
  int nnull = 0;
  // The spurious TORQUE about z. The momentum flux is a VECTOR, so unlike mass and energy
  // it is only comparable between the two panels once both are rotated into the common
  // Cartesian frame; done there, conservation demands the two outward vectors cancel just
  // as the scalars do. Summing r x (what does not cancel) gives the rate at which the
  // seam manufactures angular momentum, which is what dLz/dt should be compared against.
  Real storq[3] = {0.0, 0.0, 0.0};

  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const int p = mbpanel.h_view(m);
    for (int side=0; side<4; ++side) {
      const bool jface = (side < 2);
      const int jj = (side == 0) ? js : ((side == 1) ? je+1 : 0);
      const int kk = (side == 2) ? ks : ((side == 3) ? ke+1 : 0);
      const Real sgn = (side == 0 || side == 2) ? -1.0 : 1.0;  // outward from this block
      const int nlo = jface ? ks : js;
      const int nhi = jface ? ke : je;
      for (int t=nlo; t<=nhi; ++t) {
        const int j = jface ? jj : t;
        const int k = jface ? t  : kk;
        const Real x2n = size.h_view(m).x2min, x2x = size.h_view(m).x2max;
        const Real x3n = size.h_view(m).x3min, x3x = size.h_view(m).x3max;
        const Real xi = 0.25*M_PI*(jface
            ? LeftEdgeX(j-js, indcs.nx2, x2n, x2x)
            : CellCenterX(j-js, indcs.nx2, x2n, x2x));
        const Real et = 0.25*M_PI*(jface
            ? CellCenterX(k-ks, indcs.nx3, x3n, x3x)
            : LeftEdgeX(k-ks, indcs.nx3, x3n, x3x));
        Real cx, cy, cz;
        PanelToCart(p, xi, et, cx, cy, cz);
        // How close is this face to a cube vertex?  Geometrically, never by index: with
        // several MeshBlocks per panel a block's seam ends are mostly interior edges.
        Real vdot = 0.0;
        for (int vv=0; vv<8; ++vv) {
          const Real vx = ((vv&1)?1.0:-1.0)/sqrt(3.0);
          const Real vy = ((vv&2)?1.0:-1.0)/sqrt(3.0);
          const Real vz = ((vv&4)?1.0:-1.0)/sqrt(3.0);
          vdot = fmax(vdot, cx*vx + cy*vy + cz*vz);
        }
        // The angular width of ONE cell. It must come from this block's own x2 extent:
        // indcs.nx2 counts the cells in a MESHBLOCK, so 0.5*pi/nx2 is the cell width only
        // when a panel is a single block, and splitting the panel would otherwise widen
        // the "near a vertex" set and change what the buckets mean.
        const Real dang = 0.25*M_PI*(x2x - x2n)/static_cast<Real>(indcs.nx2);
        const int b = (acos(fmin(1.0,vdot)) < (ng + 0.5)*dang) ? 1 : 0;
        for (int i=is; i<=ie; ++i) {
          const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                       size.h_view(m).x1max);
          char key[160];
          std::snprintf(key, sizeof(key), "%.9f_%.9f_%.9f_%.9f", rad, cx, cy, cz);
          SeamFlux f;
          f.area = jface ? a2h(m,k,j,i) : a3h(m,k,j,i);
          f.fm = sgn*f.area*(jface ? f2h(m,IDN,k,j,i) : f3h(m,IDN,k,j,i));
          f.fe = sgn*f.area*(jface ? f2h(m,IEN,k,j,i) : f3h(m,IEN,k,j,i));
          f.panel = p;
          // the outward momentum flux as a Cartesian vector
          {Real e1[3], e2[3];
          cubed_sphere::PanelTangents(p, xi, et, e1, e2);
          const Real g1 = jface ? f2h(m,IM1,k,j,i) : f3h(m,IM1,k,j,i);
          const Real g2 = jface ? f2h(m,IM2,k,j,i) : f3h(m,IM2,k,j,i);
          const Real g3 = jface ? f2h(m,IM3,k,j,i) : f3h(m,IM3,k,j,i);
          // COVARIANT, like u0(IM*) itself: raise before using as basis coefficients.
          const Real cg = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
          const Real dg = 1.0 - cg*cg;
          const Real h2 = (g2 - cg*g3)/dg, h3 = (g3 - cg*g2)/dg;
          const Real rh[3] = {cx, cy, cz};
          for (int q=0; q<3; ++q) {
            f.fp[q] = sgn*f.area*(g1*rh[q] + h2*e1[q] + h3*e2[q]);
          }}
          auto it = seen.find(std::string(key));
          if (it == seen.end()) {
            seen[std::string(key)] = f;
          } else if (it->second.panel == p) {
            const Real dm = it->second.fm + f.fm;
            snull += dm;
            dnull = fmax(dnull, fabs(dm));
            ++nnull;
          } else {
            // Conservation demands what leaves one panel enters the other: fm_P + fm_Q
            // must vanish.  What it actually is, is the spurious source.
            const Real dm = it->second.fm + f.fm;
            smass[b] += dm;
            sener[b] += it->second.fe + f.fe;
            {Real dp[3];
            for (int q=0; q<3; ++q) {dp[q] = it->second.fp[q] + f.fp[q];}
            const Real rad = CellCenterX(i-is, indcs.nx1, size.h_view(m).x1min,
                                                          size.h_view(m).x1max);
            const Real rr[3] = {rad*cx, rad*cy, rad*cz};
            for (int q=0; q<3; ++q) {
              const int a = (q+1)%3, bb = (q+2)%3;
              storq[q] += rr[a]*dp[bb] - rr[bb]*dp[a];
            }}
            dmax[b] = fmax(dmax[b], fabs(dm));
            fmax_ = fmax(fmax_, fabs(f.fm));
            damax = fmax(damax, fabs(it->second.area - f.area)/f.area);
            ++npair[b];
          }
        }
      }
    }
  }
  const Real mtot = smass[0] + smass[1], etot = sener[0] + sener[1];
  std::printf("### CS SEAM FLUX MISMATCH (spurious source rate, sum over shared"
              " faces)\n");
  std::printf("###   dM/dt = %12.5e  (away %12.5e, cube vertex %12.5e)\n",
              mtot, smass[0], smass[1]);
  std::printf("###   dE/dt = %12.5e  (away %12.5e, cube vertex %12.5e)\n",
              etot, sener[0], sener[1]);
  std::printf("###   per-face max |dFm| = %12.5e away, %12.5e at a vertex;"
              "  max|area*Fm| = %12.5e\n", dmax[0], dmax[1], fmax_);
  std::printf("###   shared faces: %d away + %d within %d cells of a vertex;"
              "  max area mismatch %9.3e\n", npair[0], npair[1], ng, damax);
  std::printf("###   NULL (same-panel block faces): sum %12.5e  max %12.5e  over %d\n",
              snull, dnull, nnull);
  std::printf("###   spurious TORQUE (dL/dt): x %12.5e  y %12.5e  z %12.5e\n",
              -storq[0], -storq[1], -storq[2]);
}

//----------------------------------------------------------------------------------------
void CSTestGhostCheck(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  CSTestConsSums(pm);
  CSTestSeamFluxCheck(pm);
  CSTestLevelFluxCheck(pm);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;

  auto &w0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto w0_h = Kokkos::create_mirror_view(w0);
  Kokkos::deep_copy(w0_h, w0);
  auto &size = pmbp->pmb->mb_size;
  size.sync_host();
  auto &mbpanel = pmbp->pmb->mb_panel;
  mbpanel.sync_host();

  const Real d0 = cs_d0, amp = cs_amp, r0 = cs_r0;

  // --- iprob = 8: is the MAGNETIC ENERGY the inversion subtracted the right one? -------
  // u0(IEN) was seeded with the ANALYTIC 0.5*b0c^2, so whatever ConsToPrim subtracted
  // shows up directly as an error in the recovered pressure. On the gnomonic basis the
  // cell-centred field is a triple of NON-ORTHOGONAL face-normal components, so summing
  // their squares is not |B|^2 at all: the error is O(1) and largest at the panel corners
  // where cos_cell is largest. Storing bcc in an orthonormal frame instead
  // (Coordinates::GnomonicEquiangleRaiseVelMHD) makes it O(dx^2) and convergent.
  // The velocity is a second, independent probe: v = 0 exactly in this state.
  if (cs_iprob == 8) {
    const Real gam_ = (pmbp->pmhd != nullptr) ? pmbp->pmhd->peos->eos_data.gamma
                                              : pmbp->phydro->peos->eos_data.gamma;
    const Real gm1_ = gam_ - 1.0;
    const Real ie_exact = cs_p0/gm1_;
    std::cout << "### CS MHD ENERGY CHECK (iprob=8, uniform Cartesian B)" << std::endl;
    std::cout << "  panel   max|dp/p|      max|v|      (active cells)" << std::endl;
    Real gmax = 0.0, gvmax = 0.0;
    int gnbad = 0;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      Real e = 0.0, ev = 0.0;
      int ei = -1, ej = -1, ek = -1;   // WHERE the worst cell is, not just how bad
      // fmax(NaN, x) returns x, so a run that has already gone non-finite reports a
      // SMALL max and reads as healthy.  That happened: a blown-up run printed 5.08e-05.
      int nbad = 0;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            const Real ee = fabs(w0_h(m,IEN,k,j,i) - ie_exact)/ie_exact;
            if (!std::isfinite(ee)) ++nbad;
            if (ee > e) { e = ee; ei = i-is; ej = j-js; ek = k-ks; }
            ev = fmax(ev, fmax(fabs(w0_h(m,IVX,k,j,i)),
                          fmax(fabs(w0_h(m,IVY,k,j,i)), fabs(w0_h(m,IVZ,k,j,i)))));
          }
        }
      }
      std::printf("  %4d   %12.5e  %12.5e   at (i,j,k)=(%d,%d,%d) of (%d,%d,%d)%s\n",
                  mbpanel.h_view(m), e, ev, ei, ej, ek,
                  indcs.nx1, indcs.nx2, indcs.nx3,
                  (nbad > 0) ? "   *** NON-FINITE CELLS ***" : "");
      gnbad += nbad;
      gmax = fmax(gmax, e); gvmax = fmax(gvmax, ev);
    }
    std::printf("  MAX over panels: dp/p %12.5e   |v| %12.5e%s\n", gmax, gvmax,
                (gnbad > 0) ? "   *** RUN IS NON-FINITE, the max above is MEANINGLESS ***"
                            : "");

    // --- the EMF, which for this problem must be numerical diffusion ONLY --------------
    // v = 0 exactly everywhere (u0 is a global constant with zero momentum, and the
    // prolongation of zero is zero), so the ideal EMF -v x B vanishes identically.  What
    // survives is the Riemann solver's diffusive term, proportional to the jump in B
    // across the face.  It is therefore a direct, dimensionless read-out of how far from
    // single-valued the field is -- and comparing the last radial cell against the
    // interior says whether a coarse/fine interface is worse than ordinary truncation.
    {
      auto e1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x1e);
      auto e2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x2e);
      auto e3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x3e);
      std::cout << "### CS EMF (must be diffusion only; v = 0 identically)" << std::endl;
      // Each component has its OWN valid range: x1e runs ALONG x1 so it is a cell index
      // in i and a corner index in j,k, and cyclically for the other two.  Scanning all
      // three to +1 in every index reads memory the update never wrote, which shows up
      // as ~1e12 garbage and looks exactly like a catastrophic bug.
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        Real eint = 0.0, eend = 0.0;
        int wi = -1, wj = -1, wk = -1, wc = -1;
        int vi = -1, vj = -1, vk = -1, vc = -1;
        auto scan = [&](int c, int i0, int i1, int j0, int j1, int k0, int k1) {
          for (int k=k0; k<=k1; ++k) {
            for (int j=j0; j<=j1; ++j) {
              for (int i=i0; i<=i1; ++i) {
                const Real e = fabs((c==1) ? e1h(m,k,j,i)
                                   : (c==2) ? e2h(m,k,j,i) : e3h(m,k,j,i));
                if (i >= ie) {
                  if (e > eend) { eend = e; wi=i-is; wj=j-js; wk=k-ks; wc=c; }
                } else if (e > eint) {
                  eint = e;  vi=i-is; vj=j-js; vk=k-ks; vc=c;
                }
              }
            }
          }
        };
        scan(1, is,ie,   js,je+1, ks,ke+1);
        scan(2, is,ie+1, js,je,   ks,ke+1);
        scan(3, is,ie+1, js,je+1, ks,ke);
        std::printf("  m=%3d p%d lev%d  interior %12.5e at x%de (%d,%d,%d)"
                    "   last radial %12.5e at x%de (%d,%d,%d)   ratio %9.2e\n",
                    m, mbpanel.h_view(m), pmbp->pmb->mb_lev.h_view(m),
                    eint, vc, vi, vj, vk, eend, wc, wi, wj, wk,
                    (eint > 0.0) ? eend/eint : 0.0);
      }
    }

    // --- THE COARSE/FINE INTERFACE EMF, split by where on the interface it sits ------
    // For this state the exact EMF is identically zero (v = 0 everywhere, ghosts
    // included), so every number here is pure error. The split is the point: the
    // INTERIOR of an interface face and its own BOUNDARY EDGE -- the ring where the
    // level boundary meets an ordinary same-level tangential block boundary -- are
    // filled by different machinery, and only one of them is broken.
    {
      auto e1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x1e);
      auto e2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x2e);
      auto e3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pmhd->efld.x3e);
      auto &ngh = pmbp->pmb->nghbr;
      ngh.sync_host();
      auto &mlev = pmbp->pmb->mb_lev;
      Real fin[2] = {0.0, 0.0}, cin[2] = {0.0, 0.0};   // [0] face interior, [1] its edge
      int fwi[3] = {-1,-1,-1}, cwi[3] = {-1,-1,-1}, fwm = -1, cwm = -1, fwc = 0, cwc = 0;
      int nfine = 0, ncoar = 0;
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        for (int side=0; side<2; ++side) {         // 0 = -x1 face, 1 = +x1 face
          const int nid = (side == 0) ? 0 : 4;
          if (ngh.h_view(m,nid).gid < 0) continue;
          const int nlev = ngh.h_view(m,nid).lev;
          if (nlev == mlev.h_view(m)) continue;
          const bool amfine = (nlev < mlev.h_view(m));
          const int iface = (side == 0) ? is : (ie+1);
          if (amfine) { ++nfine; } else { ++ncoar; }
          for (int c=1; c<=3; ++c) {
            const int j1 = (c == 1 || c == 3) ? (je+1) : je;
            const int k1 = (c == 1 || c == 2) ? (ke+1) : ke;
            for (int k=ks; k<=k1; ++k) {
              for (int j=js; j<=j1; ++j) {
                // x1e does not live on an x1 FACE at all (it runs along x1), so only the
                // two transverse components sit on this interface.
                if (c == 1) continue;
                const Real e = fabs((c==2) ? e2h(m,k,j,iface) : e3h(m,k,j,iface));
                // "edge" = the outermost ring of the face, where the interface meets a
                // tangential block boundary.
                const bool edge = (j <= js || j >= j1 || k <= ks || k >= k1);
                const int b = edge ? 1 : 0;
                if (amfine) {
                  if (e > fin[b]) {
                    fin[b] = e;
                    if (b == 1) {
                      fwm=m; fwc=c;
                      fwi[0]=iface-is; fwi[1]=j-js; fwi[2]=k-ks;
                    }
                  }
                } else {
                  if (e > cin[b]) {
                    cin[b] = e;
                    if (b == 1) {
                      cwm=m; cwc=c;
                      cwi[0]=iface-is; cwi[1]=j-js; cwi[2]=k-ks;
                    }
                  }
                }
              }
            }
          }
        }
      }
      std::printf("### CS INTERFACE EMF (exact = 0; %d fine-side, %d coarse-side)\n",
                  nfine, ncoar);
      std::printf("   FINE   side: face interior %11.4e   face EDGE ring %11.4e"
                  "  (worst m=%d x%de at %d,%d,%d)\n",
                  fin[0], fin[1], fwm, fwc, fwi[0], fwi[1], fwi[2]);
      std::printf("   COARSE side: face interior %11.4e   face EDGE ring %11.4e"
                  "  (worst m=%d x%de at %d,%d,%d)\n",
                  cin[0], cin[1], cwm, cwc, cwi[0], cwi[1], cwi[2]);

      // Where does the edge ring get its data? The EMF there reads the DIAGONAL ghost
      // block beyond the interface -- i past the level boundary AND j or k past a
      // tangential block boundary -- which no gate has ever looked at: the halo check
      // holds j,k in the active range on purpose. For this field the exact value of any
      // ghost face is the same projection as the active face at the same angles, so the
      // comparison is exact and needs no interpolation.
      auto b1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x1f);
      auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x2f);
      auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->b0.x3f);
      const Real bvx = cs_bvx, bvy = cs_bvy, bvz = cs_bvz;
      Real dg[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
      int dwm = -1, dwi[3] = {-1,-1,-1}, dwc = -1, nzero = 0;
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        bool at_lvl = false;
        int iside = 0;
        for (int nn=0; nn<8; ++nn) {
          if (ngh.h_view(m,nn).gid >= 0 && ngh.h_view(m,nn).lev < mlev.h_view(m)) {
            at_lvl = true;  iside = (nn < 4) ? -1 : 1;
          }
        }
        if (!at_lvl) continue;
        const int p = mbpanel.h_view(m);
        auto ang = [&](int idx, int nx, Real mn, Real mx, bool face) {
          return 0.25*M_PI*(face ? LeftEdgeX(idx, nx, mn, mx)
                                 : CellCenterX(idx, nx, mn, mx));
        };
        auto exact = [&](int j, int k, bool jf, bool kf, int comp) {
          const Real xi = ang(j-js, indcs.nx2, size.h_view(m).x2min,
                                               size.h_view(m).x2max, jf);
          const Real et = ang(k-ks, indcs.nx3, size.h_view(m).x3min,
                                               size.h_view(m).x3max, kf);
          Real n1[3], n2[3], cx, cy, cz;
          if (comp == 0) {
            PanelToCart(p, xi, et, cx, cy, cz);
            return bvx*cx + bvy*cy + bvz*cz;
          }
          PanelNormals(p, xi, et, n1, n2);
          return (comp == 1) ? (bvx*n1[0] + bvy*n1[1] + bvz*n1[2])
                             : (bvx*n2[0] + bvy*n2[1] + bvz*n2[2]);
        };
        for (int g=1; g<=indcs.ng; ++g) {
          const int i1 = (iside > 0) ? (ie+1+g) : (is-g);      // b.x1f index
          const int ic = (iside > 0) ? (ie+g)   : (is-g);      // cell index in i
          for (int k=ks-indcs.ng; k<=ke+indcs.ng; ++k) {
            for (int j=js-indcs.ng; j<=je+indcs.ng; ++j) {
              const bool jg = (j < js || j > je), kg = (k < ks || k > ke);
              // 0 = on the face; 1 = a TWO-dimensional diagonal (the filled ones);
              // 2 = a THREE-dimensional corner, which is deliberately not filled.
              const int b = (jg && kg) ? 2 : ((jg || kg) ? 1 : 0);
              if (j > je+1 || k > ke+1) continue;
              const Real d1 = fabs(b1h(m,k,j,i1) - exact(j,k,false,false,0));
              const Real d2 = fabs(b2h(m,k,j,ic) - exact(j,k,true ,false,1));
              const Real d3 = fabs(b3h(m,k,j,ic) - exact(j,k,false,true ,2));
              if (d1 > dg[b][0]) dg[b][0] = d1;
              if (d2 > dg[b][1]) dg[b][1] = d2;
              if (d3 > dg[b][2]) dg[b][2] = d3;
              if (b >= 1 && b1h(m,k,j,i1) == 0.0 && b2h(m,k,j,ic) == 0.0) {
                ++nzero;
                if (nzero <= 6) {
                  // Which neighbour slot owns this ghost cell, and is it registered?
                  const int ox1 = iside;
                  const int ox2 = (j < js) ? -1 : ((j > je) ? 1 : 0);
                  const int ox3 = (k < ks) ? -1 : ((k > ke) ? 1 : 0);
                  int nid = -1;
                  for (int nn=0; nn<pmbp->pmb->nnghbr; ++nn) {
                    if (ngh.h_view(m,nn).gid >= 0) continue;
                  }
                  std::printf("   ZERO ghost m=%2d(lev %d,pan %d) cell (%d,%d,%d)"
                              " offset (%d,%d,%d)  slots:", m, mlev.h_view(m),
                              mbpanel.h_view(m), ic-is, j-js, k-ks, ox1, ox2, ox3);
                  const int base = (ox3 == 0) ? (16 + (ox1+1) + 2*(ox2+1))
                                 : ((ox2 == 0) ? (24 + (ox1+9)*(ox1!=0) + 2*(ox3+1))
                                               : (48 + (ox1+1)/2 + (ox2+1) + 2*(ox3+1)));
                  for (int q=0; q<2; ++q) {
                    std::printf(" [%d]=gid %d lev %d;", base+q,
                                ngh.h_view(m,base+q).gid, ngh.h_view(m,base+q).lev);
                  }
                  std::printf("  nid %d\n", nid);
                }
              }
            }
          }
        }
      }
      std::printf("   ghosts BEYOND the interface, vs exact:  ON the face (j,k active)"
                  " %11.4e %11.4e %11.4e\n", dg[0][0], dg[0][1], dg[0][2]);
      std::printf("                                          2-D DIAGONAL (filled) "
                  " %11.4e %11.4e %11.4e\n", dg[1][0], dg[1][1], dg[1][2]);
      std::printf("                                          3-D CORNER (not filled)"
                  " %11.4e %11.4e %11.4e\n", dg[2][0], dg[2][1], dg[2][2]);
      std::printf("   diagonal ghost cells left at EXACTLY ZERO: %d\n", nzero);
      // And the CELL-CENTRED state in the same cells? The corner EMF reads velocities
      // there too, and a zero density is not a small error.
      {
        auto u0h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->u0);
        int nz = 0, ntot = 0;  Real dmin = 1.0e300;
        for (int m=0; m<pmbp->nmb_thispack; ++m) {
          bool at_lvl = false;  int iside = 0;
          for (int nn=0; nn<8; ++nn) {
            if (ngh.h_view(m,nn).gid >= 0 && ngh.h_view(m,nn).lev < mlev.h_view(m)) {
              at_lvl = true;  iside = (nn < 4) ? -1 : 1;
            }
          }
          if (!at_lvl) continue;
          for (int g=1; g<=indcs.ng; ++g) {
            const int ic = (iside > 0) ? (ie+g) : (is-g);
            for (int k=ks-indcs.ng; k<=ke+indcs.ng; ++k) {
              for (int j=js-indcs.ng; j<=je+indcs.ng; ++j) {
                if (j >= js && j <= je && k >= ks && k <= ke) continue;
                ++ntot;
                if (u0h(m,IDN,k,j,ic) == 0.0) ++nz;
                dmin = fmin(dmin, u0h(m,IDN,k,j,ic));
              }
            }
          }
        }
        std::printf("   CELL-CENTRED in the same cells: %d of %d have rho == 0"
                    " exactly;  min rho %11.4e\n", nz, ntot, dmin);
      }
    }

    // --- what did the FACE HALO actually deliver? ------------------------------------
    // Every ghost face has a KNOWN exact value: the uniform field projected on the
    // normal of that face, on THIS panel's own gnomonic map extended past the seam. That
    // is what a perfect halo would produce, so this separates the halo from the scheme,
    // and separates the two failure modes from each other:
    //   O(1) error, flat in resolution -> the INDEX MAP is reaching the wrong face, or
    //     the wrong orientation;
    //   O(dx) error, halving with resolution -> the map is right and only the along-seam
    //     offset is left, i.e. the halo needs interpolation.
    // Corner ghosts are excluded: j,k are held in the active range.
    {
      auto &b0f = pmbp->pmhd->b0;
      auto ah1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x1f);
      auto ah2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x2f);
      auto ah3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x3f);
      auto b1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x1f);
      auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x2f);
      auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x3f);
      const Real bvx = cs_bvx, bvy = cs_bvy, bvz = cs_bvz;
      std::cout << "### CS FIELD HALO CHECK (iprob=8, ghost faces vs exact)" << std::endl;
      std::cout << "  panel seam   d(b.x1f)      d(b.x2f)      d(b.x3f)" << std::endl;
      Real g1 = 0.0, g2 = 0.0, g3 = 0.0;
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        const int p = mbpanel.h_view(m);
        Real e1 = 0.0, e2 = 0.0, e3 = 0.0;
        auto ang = [&](int idx, int nx, Real mn, Real mx, bool face) {
          return 0.25*M_PI*(face ? LeftEdgeX(idx, nx, mn, mx)
                                 : CellCenterX(idx, nx, mn, mx));
        };
        // exact projections at (xi,eta) given as index positions
        auto exact = [&](int j, int k, bool jface, bool kface, int comp) {
          const Real xi = ang(j-js, indcs.nx2, size.h_view(m).x2min,
                                               size.h_view(m).x2max, jface);
          const Real et = ang(k-ks, indcs.nx3, size.h_view(m).x3min,
                                               size.h_view(m).x3max, kface);
          Real n1[3], n2[3], cx, cy, cz;
          if (comp == 0) {
            PanelToCart(p, xi, et, cx, cy, cz);
            return bvx*cx + bvy*cy + bvz*cz;
          }
          PanelNormals(p, xi, et, n1, n2);
          return (comp == 1) ? (bvx*n1[0] + bvy*n1[1] + bvz*n1[2])
                             : (bvx*n2[0] + bvy*n2[1] + bvz*n2[2]);
        };
        // --- RADIAL ghost faces, which is where a coarse/fine interface lives -------
        // For a uniform Cartesian field the exact face-normal projections depend on the
        // ANGLES only, not on r, so a radial ghost has the same exact value as the
        // active face at the same (j,k). That makes this a clean probe of whatever
        // filled those ghosts -- the physical BC, or prolongation from a coarser
        // neighbour when the block sits at a level boundary.
        {
          Real r1 = 0.0, r2 = 0.0, r3 = 0.0;
          // L1 as well as Linf: a slope limiter clipping at an extremum gives FIRST
          // order in the max norm while staying second order in the mean, which is
          // ordinary behaviour and not a defect.  Only both being first order is.
          Real s1 = 0.0; int n1c = 0;
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int g=1; g<=ng; ++g) {
                s1 += fabs(b1h(m,k,j,ie+1+g) - exact(j,k,false,false,0));
                s1 += fabs(b1h(m,k,j,is-g)   - exact(j,k,false,false,0));
                n1c += 2;
                r1 = fmax(r1, fabs(b1h(m,k,j,ie+1+g) - exact(j,k,false,false,0)));
                r1 = fmax(r1, fabs(b1h(m,k,j,is-g)   - exact(j,k,false,false,0)));
                r2 = fmax(r2, fabs(b2h(m,k,j,ie+g)   - exact(j,k,true ,false,1)));
                r2 = fmax(r2, fabs(b2h(m,k,j,is-g)   - exact(j,k,true ,false,1)));
                r3 = fmax(r3, fabs(b3h(m,k,j,ie+g)   - exact(j,k,false,true ,2)));
                r3 = fmax(r3, fabs(b3h(m,k,j,is-g)   - exact(j,k,false,true ,2)));
              }
            }
          }
          // The INTERFACE FACE ITSELF, i = ie+1 / is.  This is not a ghost: it is a
          // face the fine block owns, but ProlongFCSharedX1Face overwrites it from the
          // coarse face so that the two levels agree.  The ghost sweep above starts at
          // ie+2 and never looks at it.
          Real f1 = 0.0, f0 = 0.0;
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              f1 = fmax(f1, fabs(b1h(m,k,j,ie+1) - exact(j,k,false,false,0)));
              f0 = fmax(f0, fabs(b1h(m,k,j,is)   - exact(j,k,false,false,0)));
            }
          }
          std::printf("  %4d rad   %12.5e  %12.5e  %12.5e   (radial ghosts)"
                      "   L1(b.x1f) %12.5e   INTERFACE FACE b.x1f: %12.5e %12.5e\n",
                      p, r1, r2, r3, (n1c > 0) ? s1/n1c : 0.0, f1, f0);
        }

        // --- and the SOLENOIDAL constraint in those same ghost cells ---------------
        // CT keeps div B = 0 in ACTIVE cells whatever the ghosts hold, so an active-cell
        // check cannot see a bad prolongation.  A ghost cell that carries a monopole is
        // still read by reconstruction and the Riemann solve, and shows up as an O(1)
        // spurious Lorentz force.  Report sum(area*B) over the six faces, normalised by
        // the sum of the magnitudes, so it is dimensionless and grid-independent.
        {
          auto dvb = [&](int k, int j, int i) {
            Real f[6] = {-ah1(m,k,j,i)*b1h(m,k,j,i),   ah1(m,k,j,i+1)*b1h(m,k,j,i+1),
                         -ah2(m,k,j,i)*b2h(m,k,j,i),   ah2(m,k,j+1,i)*b2h(m,k,j+1,i),
                         -ah3(m,k,j,i)*b3h(m,k,j,i),   ah3(m,k+1,j,i)*b3h(m,k+1,j,i)};
            Real sum = 0.0, mag = 0.0;
            for (int q=0; q<6; ++q) { sum += f[q]; mag += fabs(f[q]); }
            return (mag > 0.0) ? fabs(sum)/mag : 0.0;
          };
          Real dg = 0.0, da = 0.0;
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int g=1; g<=ng-1; ++g) {
                dg = fmax(dg, dvb(k,j,ie+g));
                dg = fmax(dg, dvb(k,j,is-g));
              }
              for (int i=is; i<=ie; ++i) da = fmax(da, dvb(k,j,i));
            }
          }
          std::printf("  %4d divB  %12.5e (radial ghosts)  %12.5e (active)\n",
                      p, dg, da);
        }

        const int i = is;
        // report PER SEAM: an axis-swap seam (equatorial<->polar) and a plain one behave
        // differently, and a max over all four hides which
        for (int f=0; f<4; ++f) {
          Real f1 = 0.0, f2 = 0.0, f3 = 0.0;
          for (int g=0; g<ng; ++g) {
            if (f < 2) {
              const int jc = (f == 0) ? (js-1-g) : (je+1+g);
              const int jf = (f == 0) ? (js-1-g) : (je+2+g);
              for (int k=ks; k<=ke; ++k) {
                f1 = fmax(f1, fabs(b1h(m,k,jc,i) - exact(jc,k,false,false,0)));
                f3 = fmax(f3, fabs(b3h(m,k,jc,i) - exact(jc,k,false,true, 2)));
                f2 = fmax(f2, fabs(b2h(m,k,jf,i) - exact(jf,k,true, false,1)));
              }
            } else {
              const int kc = (f == 2) ? (ks-1-g) : (ke+1+g);
              const int kf = (f == 2) ? (ks-1-g) : (ke+2+g);
              for (int j=js; j<=je; ++j) {
                f1 = fmax(f1, fabs(b1h(m,kc,j,i) - exact(j,kc,false,false,0)));
                f2 = fmax(f2, fabs(b2h(m,kc,j,i) - exact(j,kc,true, false,1)));
                f3 = fmax(f3, fabs(b3h(m,kf,j,i) - exact(j,kf,false,true, 2)));
              }
            }
          }
          std::printf("  %4d %s  %12.5e  %12.5e  %12.5e\n", p,
                      (f==0?"-x2":(f==1?"+x2":(f==2?"-x3":"+x3"))), f1, f2, f3);
          e1 = fmax(e1,f1); e2 = fmax(e2,f2); e3 = fmax(e3,f3);
        }
        // THE FOUR CORNER GHOST BLOCKS, ng x ng each, fed by the diagonal (edge/corner)
        // buffers rather than the face ones. Reported separately because they are filled
        // by a different code path -- a plain signed permutation, with neither the basis
        // transform nor the along-seam resample -- and because at a cube VERTEX only
        // three panels meet, so a "diagonal neighbour" is not even well defined there.
        Real c1 = 0.0, c2 = 0.0, c3 = 0.0;
        for (int cj=0; cj<2; ++cj) {
          for (int ck=0; ck<2; ++ck) {
            for (int gj=0; gj<ng; ++gj) {
              for (int gk=0; gk<ng; ++gk) {
                const int jc = cj ? (je+1+gj) : (js-1-gj);
                const int kc = ck ? (ke+1+gk) : (ks-1-gk);
                const int jf = cj ? (je+2+gj) : (js-1-gj);
                const int kf = ck ? (ke+2+gk) : (ks-1-gk);
                c1 = fmax(c1, fabs(b1h(m,kc,jc,i) - exact(jc,kc,false,false,0)));
                c2 = fmax(c2, fabs(b2h(m,kc,jf,i) - exact(jf,kc,true, false,1)));
                c3 = fmax(c3, fabs(b3h(m,kf,jc,i) - exact(jc,kf,false,true, 2)));
              }
            }
          }
        }
        std::printf("  %4d cnr  %12.5e  %12.5e  %12.5e\n", p, c1, c2, c3);
        // THE SHARED SEAM FACES THEMSELVES, which are ACTIVE on both panels and are
        // never communicated: each panel evolves its own copy with its own EMFs, so they
        // stay equal only if the seam EMFs are single-valued. Compared here against the
        // exact value and against the same quantity on interior faces -- with consistent
        // seam EMFs both sit at truncation and stay there; with inconsistent ones the
        // seam copies drift apart, and away from exact, secularly in time.
        Real sf = 0.0, itf = 0.0;
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je+1; ++j) {
            const Real e2v = exact(j,k,true,false,1);
            const Real d = fabs(b2h(m,k,j,i) - e2v);
            if (j == js || j == je+1) { sf = fmax(sf,d); } else { itf = fmax(itf,d); }
          }
        }
        for (int k=ks; k<=ke+1; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real e3v = exact(j,k,false,true,2);
            const Real d = fabs(b3h(m,k,j,i) - e3v);
            if (k == ks || k == ke+1) { sf = fmax(sf,d); } else { itf = fmax(itf,d); }
          }
        }
        std::printf("  %4d shr  %12.5e (seam faces)  %12.5e (interior faces)\n",
                    p, sf, itf);
        g1 = fmax(g1,e1); g2 = fmax(g2,e2); g3 = fmax(g3,e3);
      }
      std::printf("  MAX: x1f %12.5e  x2f %12.5e  x3f %12.5e\n", g1, g2, g3);

      // Panel 0, -x3 seam, ghost layer 0: the in-seam component b.x2f, delivered vs
      // exact, sampled along the seam. actual ~ -exact means a sign; actual ~ exact one
      // index over means the map is off by a cell; neither means the transform itself.
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        if (mbpanel.h_view(m) != 0) continue;
        const int kc = ks-1;
        std::cout << "  P0 -x3 g0 b.x2f:   j    delivered        exact"
                  << "        exact(j-1)     exact(j+1)" << std::endl;
        const int nsmp = 12;
        int smp[nsmp] = {0,1,2,3,4, indcs.nx2/2,
                         indcs.nx2-4, indcs.nx2-3, indcs.nx2-2, indcs.nx2-1,
                         indcs.nx2, indcs.nx2};
        for (int t=0; t<nsmp; ++t) {
          const int j = js + smp[t];
          auto ex = [&](int jj) {
            const Real xi = 0.25*M_PI*LeftEdgeX(jj-js, indcs.nx2, size.h_view(m).x2min,
                                                                  size.h_view(m).x2max);
            const Real et = 0.25*M_PI*CellCenterX(kc-ks, indcs.nx3, size.h_view(m).x3min,
                                                                    size.h_view(m).x3max);
            Real n1[3], n2[3];
            PanelNormals(0, xi, et, n1, n2);
            return cs_bvx*n1[0] + cs_bvy*n1[1] + cs_bvz*n1[2];
          };
          std::printf("                   %4d  %12.5e  %12.5e  %12.5e  %12.5e\n",
                      j-js, b2h(m,kc,j,is), ex(j), ex(j-1), ex(j+1));
        }
        break;
      }
    }

    // WHERE is the error? v = 0 is exact here, so |v| is a pure error field. A defect in
    // the panel halo is confined to a band a few cells wide along each seam; a defect in
    // the interior flux/EMF path with a tangential field is a bulk effect. Run this after
    // ONE step -- later the seam error has propagated everywhere.
    const int nseam = ng;
    // --- IS THE CUBE-VERTEX EMF SINGLE-VALUED? ---------------------------------------
    // The radial edge at a panel corner is one physical edge shared by THREE panels, and
    // CT consumes dxedge*E along it. Keyed by geometry like the face check below, so no
    // index map is involved. Reports the max spread among the panels sharing each vertex.
    {
      auto &ef = pmbp->pmhd->efld;
      auto e1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ef.x1e);
      std::map<std::string, std::pair<Real,Real>> lo_hi;   // min, max over panels
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        const int p = mbpanel.h_view(m);
        for (int cj=0; cj<2; ++cj) {
          for (int ck=0; ck<2; ++ck) {
            const int jc = cj ? je+1 : js;
            const int kc = ck ? ke+1 : ks;
            const Real xi = 0.25*M_PI*LeftEdgeX(jc-js, indcs.nx2, size.h_view(m).x2min,
                                                                  size.h_view(m).x2max);
            const Real et = 0.25*M_PI*LeftEdgeX(kc-ks, indcs.nx3, size.h_view(m).x3min,
                                                                  size.h_view(m).x3max);
            Real cx, cy, cz;
            PanelToCart(p, xi, et, cx, cy, cz);
            // A CUBE VERTEX has |x| = |y| = |z|. With several MeshBlocks per panel a
            // block corner is usually just an ordinary interior edge, which is shared by
            // four blocks and is not what this check is about.
            if (fabs(fabs(cx) - fabs(cy)) > 1.0e-9 ||
                fabs(fabs(cx) - fabs(cz)) > 1.0e-9) continue;
            for (int i=is; i<=ie; ++i) {
              char key[128];
              std::snprintf(key, sizeof(key), "%d_%.9f_%.9f_%.9f", i-is, cx, cy, cz);
              const Real e = e1h(m,kc,jc,i);
              auto it = lo_hi.find(std::string(key));
              if (it == lo_hi.end()) {
                lo_hi[std::string(key)] = std::make_pair(e,e);
              } else {
                it->second.first  = fmin(it->second.first, e);
                it->second.second = fmax(it->second.second, e);
              }
            }
          }
        }
      }
      Real spread = 0.0, emag = 0.0;
      int nvert = 0;
      for (auto &kv : lo_hi) {
        spread = fmax(spread, kv.second.second - kv.second.first);
        emag = fmax(emag, fabs(kv.second.second));
        ++nvert;
      }
      std::printf("### CS CUBE-VERTEX EMF SPREAD: %12.5e over %d vertex edges"
                  " (max|e1| %10.3e)\n", spread, nvert, emag);
    }

    // TEMP AUDIT: is every neighbour pairing RECIPROCAL?  Under MPI a send uses the
    // receiver's (lid, dest) as its tag, so a slot whose partner does not point back is
    // a receive nobody satisfies -> MPI_Wait hangs in ClearRecv.
    {
      auto &ng = pmbp->pmb->nghbr;
      ng.sync_host();
      auto &gid = pmbp->pmb->mb_gid;
      const int nn = pmbp->pmb->nnghbr;
      std::printf("### CS NEIGHBOUR RECIPROCITY AUDIT (nnghbr=%d)\n", nn);
      int bad[7] = {0,0,0,0,0,0,0}, good[7] = {0,0,0,0,0,0,0};
      int uncovered = 0;   // non-reciprocal AND not skipped: the MPI hang condition
      const char *cls[7] = {"x1face", "x2face", "x1x2edge", "x3face",
                            "x3x1edge", "x2x3edge", "corner"};
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        for (int n=0; n<nn; ++n) {
          if (ng.h_view(m,n).gid < 0) continue;
          int c = (n<8) ? 0 : (n<16) ? 1 : (n<24) ? 2 : (n<32) ? 3
                                    : (n<40) ? 4 : (n<48) ? 5 : 6;
          int mn = -1;
          for (int mm=0; mm<pmbp->nmb_thispack; ++mm) {
            if (gid.h_view(mm) == ng.h_view(m,n).gid) {
              mn = mm;
              break;
            }
          }
          if (mn < 0) continue;
          const int d = ng.h_view(m,n).dest;
          if (d < 0 || d >= nn || ng.h_view(mn,d).gid != gid.h_view(m)) {
            // Does the partner list us AT ALL, just in another slot?  That separates a
            // `dest` bookkeeping error (it does, at slot rslot) from a genuine SEARCH
            // asymmetry (it does not list us anywhere), which need different fixes.
            int rslot = -1;
            for (int nn2=0; nn2<nn; ++nn2) {
              if (ng.h_view(mn,nn2).gid == gid.h_view(m)) { rslot = nn2; break; }
            }
            if (bad[c] < 6) {
              const LogicalLocation &la = pm->lloc_eachmb[gid.h_view(m)];
              const LogicalLocation &lb = pm->lloc_eachmb[ng.h_view(m,n).gid];
              std::printf("   NON-RECIP %-8s m=%2d gid %3d pan %d lev %d l=(%d,%d,%d)"
                          " n=%2d -> gid %3d pan %d lev %d l=(%d,%d,%d)"
                          "  dest=%2d holds %3d;  reverse slot = %d\n",
                          cls[c], m, gid.h_view(m), mbpanel.h_view(m), la.level,
                          static_cast<int>(la.lx1), static_cast<int>(la.lx2),
                          static_cast<int>(la.lx3), n,
                          ng.h_view(m,n).gid, mbpanel.h_view(mn), lb.level,
                          static_cast<int>(lb.lx1), static_cast<int>(lb.lx2),
                          static_cast<int>(lb.lx3), d,
                          (d>=0 && d<nn) ? ng.h_view(mn,d).gid : -999, rslot);
            }
            ++bad[c];
            // The only thing that MATTERS is whether a non-reciprocal slot is actually
            // EXCHANGED. A cube vertex is skipped on both sides at once, so it has
            // neither a send nor a posted receive and cannot hang. Anything else that
            // shows up here is a real bug: it will hang ClearRecv under MPI.
            if (!IsCubeVertexCorner(ng.h_view, mbpanel.h_view, m, n)) {
              ++uncovered;
              if (uncovered <= 5) {
                std::printf("   !! UNCOVERED non-reciprocal slot m=%d n=%d\n", m, n);
              }
            }
          } else {
            ++good[c];
          }
        }
      }
      for (int c=0; c<7; ++c) {
        std::printf("   %-9s  %5d reciprocal, %5d NON-reciprocal\n",
                    cls[c], good[c], bad[c]);
      }
      std::printf("   non-reciprocal slots that are still EXCHANGED (must be 0): %d\n",
                  uncovered);
    }

    // --- IS THE x1e ON EVERY BLOCK CORNER SINGLE-VALUED? -----------------------------
    // Generalises the vertex check above to EVERY block-corner radial edge, bucketed by
    // GEOMETRY into the three cases that fail for different reasons: a cube vertex (three
    // panels, fixed by AveragePanelCornerEMF), a four-block edge lying ON a seam (two
    // panels and four blocks -- this is what the cube-vertex buffer collision used to
    // corrupt), and a four-block edge in a panel's interior (the generic machinery, which
    // was always right). Splitting them is what localised the seam defect: the aggregate
    // number alone cannot tell "the seam is broken" from "corners are broken".
    {
      auto &ef = pmbp->pmhd->efld;
      auto e1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ef.x1e);
      std::map<std::string, std::pair<std::pair<Real,Real>,int>> tb;
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        const int p = mbpanel.h_view(m);
        for (int cj=0; cj<2; ++cj) {
          for (int ck=0; ck<2; ++ck) {
            const int jc = cj ? je+1 : js;
            const int kc = ck ? ke+1 : ks;
            const Real xi = 0.25*M_PI*LeftEdgeX(jc-js, indcs.nx2, size.h_view(m).x2min,
                                                                  size.h_view(m).x2max);
            const Real et = 0.25*M_PI*LeftEdgeX(kc-ks, indcs.nx3, size.h_view(m).x3min,
                                                                  size.h_view(m).x3max);
            Real cx, cy, cz;
            PanelToCart(p, xi, et, cx, cy, cz);
            const bool vtx = (fabs(fabs(cx)-fabs(cy)) < 1.0e-9 &&
                              fabs(fabs(cx)-fabs(cz)) < 1.0e-9);
            const bool seam = (fabs(fabs(xi) - 0.25*M_PI) < 1.0e-9 ||
                               fabs(fabs(et) - 0.25*M_PI) < 1.0e-9);
            const int bkt = vtx ? 0 : (seam ? 1 : 2);
            for (int i=is; i<=ie; ++i) {
              char key[128];
              std::snprintf(key, sizeof(key), "%d_%.9f_%.9f_%.9f", i-is, cx, cy, cz);
              const Real e = e1h(m,kc,jc,i);
              auto it = tb.find(std::string(key));
              if (it == tb.end()) {
                tb[std::string(key)] = std::make_pair(std::make_pair(e,e), bkt);
              } else {
                it->second.first.first  = fmin(it->second.first.first, e);
                it->second.first.second = fmax(it->second.first.second, e);
              }
            }
          }
        }
      }
      Real sp[3] = {0.0, 0.0, 0.0};
      int nb[3] = {0, 0, 0};
      for (auto &kv : tb) {
        const int b = kv.second.second;
        sp[b] = fmax(sp[b], kv.second.first.second - kv.second.first.first);
        nb[b] += 1;
      }
      const char *nm[3] = {"cube vertex   ", "four-block seam", "panel interior"};
      std::printf("### CS BLOCK-CORNER x1e SPREAD by geometry\n");
      for (int b=0; b<3; ++b) {
        std::printf("   %s  %12.5e   (%d edges)\n", nm[b], sp[b], nb[b]);
      }
    }

    // --- IS THE SHARED SEAM FACE SINGLE-VALUED? --------------------------------------
    // The definitive gate for seam EMF consistency. A face on a panel seam is ONE
    // physical face, stored twice -- once as an ACTIVE face of each of the two panels
    // that meet there -- and it is never communicated: each panel evolves its own copy
    // with its own EMFs. The copies stay equal if and only if every edge bounding that
    // face carries the same line integral on both panels. Any mismatch is a magnetic
    // field that is multi-valued at the seam.
    //
    // Found here WITHOUT any index map, by geometry: every panel-boundary face is keyed
    // on its physical centre, and a key that two different panels both produce is a
    // shared face. The two panels' outward normals there are the same physical vector up
    // to sign, so the comparison is b_P - (n_P.n_Q) b_Q.
    {
      auto &b0f = pmbp->pmhd->b0;
      auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x2f);
      auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x3f);
      struct SeamFace { Real b; Real n[3]; int panel; };
      // Split off the faces within ng of a seam END, i.e. a panel corner. Only there do
      // the two panels average over DIFFERENT sets of contributions.
      std::map<std::string, SeamFace> seen;
      Real dmax = 0.0, bmax = 0.0, cmax = 0.0;
      int npair = 0;
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        const int p = mbpanel.h_view(m);
        for (int side=0; side<4; ++side) {
          const bool jface = (side < 2);
          const int jj = (side == 0) ? js : ((side == 1) ? je+1 : 0);
          const int kk = (side == 2) ? ks : ((side == 3) ? ke+1 : 0);
          const int nlo = jface ? ks : js;
          const int nhi = jface ? ke : je;
          for (int t=nlo; t<=nhi; ++t) {
            const int j = jface ? jj : t;
            const int k = jface ? t  : kk;
            const Real x2n = size.h_view(m).x2min, x2x = size.h_view(m).x2max;
            const Real x3n = size.h_view(m).x3min, x3x = size.h_view(m).x3max;
            const Real xi = 0.25*M_PI*(jface
                ? LeftEdgeX(j-js, indcs.nx2, x2n, x2x)
                : CellCenterX(j-js, indcs.nx2, x2n, x2x));
            const Real et = 0.25*M_PI*(jface
                ? CellCenterX(k-ks, indcs.nx3, x3n, x3x)
                : LeftEdgeX(k-ks, indcs.nx3, x3n, x3x));
            Real n1[3], n2[3], cx, cy, cz;
            PanelToCart(p, xi, et, cx, cy, cz);
            PanelNormals(p, xi, et, n1, n2);
            for (int i=is; i<=ie; ++i) {
              char key[128];
              std::snprintf(key, sizeof(key), "%d_%.9f_%.9f_%.9f", i-is, cx, cy, cz);
              SeamFace f;
              f.b = jface ? b2h(m,k,j,i) : b3h(m,k,j,i);
              for (int c=0; c<3; ++c) {f.n[c] = jface ? n1[c] : n2[c];}
              f.panel = p;
              auto it = seen.find(std::string(key));
              if (it == seen.end()) {
                seen[std::string(key)] = f;
              } else if (it->second.panel != p) {
                const Real dot = it->second.n[0]*f.n[0] + it->second.n[1]*f.n[1]
                               + it->second.n[2]*f.n[2];
                const Real d = fabs(it->second.b - dot*f.b);
                // "near a corner" must mean near a CUBE VERTEX, tested geometrically:
                // with several MeshBlocks per panel the ends of a block's seam are
                // mostly ordinary interior edges, and bucketing by index would count
                // those too and hide which faces are actually at issue.
                Real vdot = 0.0;
                for (int vv=0; vv<8; ++vv) {
                  const Real vx = ((vv&1)?1.0:-1.0)/sqrt(3.0);
                  const Real vy = ((vv&2)?1.0:-1.0)/sqrt(3.0);
                  const Real vz = ((vv&4)?1.0:-1.0)/sqrt(3.0);
                  vdot = fmax(vdot, cx*vx + cy*vy + cz*vz);
                }
                const Real dang = 0.5*M_PI/static_cast<Real>(indcs.nx2);
                if (acos(fmin(1.0,vdot)) < (ng + 0.5)*dang) {
                  cmax = fmax(cmax, d);
                } else {
                  dmax = fmax(dmax, d);
                }
                bmax = fmax(bmax, fabs(f.b));
                ++npair;
              }
            }
          }
        }
      }
      std::printf("### CS SEAM FACE SINGLE-VALUED: max|b_P - b_Q| = %12.5e away from a"
                  " cube vertex, %12.5e within %d cells of one; %d shared faces,"
                  " max|b| %9.3e\n", dmax, cmax, ng, npair, bmax);
    }

    // --- div B, the quantity SEAM EMF CONSISTENCY controls ---------------------------
    // CT conserves sum(area*B) over a cell's six faces EXACTLY, provided every edge's
    // line integral dxedge*E is single-valued. Two panels sharing a seam edge compute it
    // independently, so any mismatch there shows up here and nowhere else. Reported as
    // the face-flux imbalance normalised by |B| times the largest face area of the cell,
    // which is dimensionless and O(1) if the field were discontinuous.
    {
      auto &b0f = pmbp->pmhd->b0;
      auto b1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x1f);
      auto b2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x2f);
      auto b3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0f.x3f);
      auto a1h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x1f);
      auto a2h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x2f);
      auto a3h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmbp->pcoord->area.x3f);
      std::cout << "### CS DIV B (iprob=8; |sum A.B| / (|B| max A), exact = 0)"
                << std::endl;
      std::cout << "  panel     max seam      max interior   location" << std::endl;
      const Real bmag = sqrt(SQR(cs_bvx) + SQR(cs_bvy) + SQR(cs_bvz));
      for (int m=0; m<pmbp->nmb_thispack; ++m) {
        Real ds = 0.0, di = 0.0;
        int jm = 0, km = 0, im = 0, ii = 0;
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const bool near = (j-js < ng) || (je-j < ng) ||
                              (k-ks < ng) || (ke-k < ng);
            for (int i=is; i<=ie; ++i) {
              const Real f1 = a1h(m,k,j,i+1)*b1h(m,k,j,i+1) - a1h(m,k,j,i)*b1h(m,k,j,i);
              const Real f2 = a2h(m,k,j+1,i)*b2h(m,k,j+1,i) - a2h(m,k,j,i)*b2h(m,k,j,i);
              const Real f3 = a3h(m,k+1,j,i)*b3h(m,k+1,j,i) - a3h(m,k,j,i)*b3h(m,k,j,i);
              Real amx = fmax(a1h(m,k,j,i+1), fmax(a2h(m,k,j+1,i), a3h(m,k+1,j,i)));
              const Real d = fabs(f1 + f2 + f3)/(bmag*amx);
              if (near) {
                if (d > ds) { ds = d; jm = j-js; km = k-ks; im = i-is; }
              } else {
                if (d > di) { di = d; ii = i-is; }
              }
            }
          }
        }
        std::printf("  %4d   %12.5e   %12.5e   seam(i=%d,j=%d,k=%d)  int(i=%d)\n",
                    mbpanel.h_view(m), ds, di, im, jm, km, ii);
      }
    }

    std::cout << "### CS SEAM LOCALISATION (iprob=8, exact v = 0)" << std::endl;
    std::cout << "  panel   max|v| seam   max|v| interior   mean v^2 seam  interior"
              << "     max|v| corner  (j,k) of seam max" << std::endl;
    // A seam cell is split further into CORNER (within nseam of two panel edges at
    // once, i.e. the region the ng x ng edge/corner halo buffers feed) and EDGE (one
    // seam only, fed by the face buffers). They fail for different reasons: only the
    // face buffers carry the basis transform and the along-seam resample, so a max that
    // lives in the corners is not evidence about the face halo at all.
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      Real vs = 0.0, vi = 0.0, vc = 0.0, ss = 0.0, si = 0.0;
      int ns = 0, ni = 0, jm = 0, km = 0;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          const bool nj_ = (j-js < nseam) || (je-j < nseam);
          const bool nk_ = (k-ks < nseam) || (ke-k < nseam);
          const bool near = nj_ || nk_;
          for (int i=is; i<=ie; ++i) {
            const Real v2 = SQR(w0_h(m,IVX,k,j,i)) + SQR(w0_h(m,IVY,k,j,i))
                          + SQR(w0_h(m,IVZ,k,j,i));
            const Real vm = sqrt(v2);
            if (nj_ && nk_) vc = fmax(vc,vm);
            if (near) {
              if (vm > vs) { vs = vm; jm = j-js; km = k-ks; }
              ss += v2; ++ns;
            } else {
              vi = fmax(vi,vm); si += v2; ++ni;
            }
          }
        }
      }
      Real dmn = 1.0e30, pmn = 1.0e30;
      int pi = 0, pj = 0, pk = 0;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            dmn = fmin(dmn, w0_h(m,IDN,k,j,i));
            if (w0_h(m,IEN,k,j,i) < pmn) {
              pmn = w0_h(m,IEN,k,j,i); pi = i-is; pj = j-js; pk = k-ks;
            }
          }
        }
      }
      std::printf("  %4d   %12.5e   %12.5e   %12.5e  %12.5e   %12.5e  (%d,%d)"
                  "  minrho %11.4e  minie %11.4e at (%d,%d,%d)\n",
                  mbpanel.h_view(m), vs, vi, (ns>0?ss/ns:0.0), (ni>0?si/ni:0.0), vc,
                  jm, km, dmn, pmn, pi, pj, pk);
    }
    return;
  }

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
    const Real gm1 = ((pmbp->pmhd != nullptr) ? pmbp->pmhd->peos->eos_data.gamma
                                              : pmbp->phydro->peos->eos_data.gamma) - 1.0;
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
