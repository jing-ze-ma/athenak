//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sp_test.cpp
//! \brief Exact-solution tests for the SPHERICAL-POLAR grid, the counterparts of the
//! cubed sphere's cs_test.  problem/iprob selects:
//!
//!   8 (default) = UNIFORM FIELD AT REST (cs_test iprob=8).  A uniform gas AT REST
//!       threaded by a UNIFORM Cartesian field B = b0 zhat (or xhat, problem/bdir = 0).
//!       curl B = 0, so the exact evolution is nothing; after one step the momentum per
//!       cell is the spurious force of the geometric source term against the flux
//!       divergence.  Measure it from two consecutive dumps.  <mhd> required.
//!
//!   3 = RIGID ROTATION (cs_test iprob=3).  v = omega zhat x r, i.e. v_phi = omega r
//!       sin(theta), rho = d0, and p = p0 + d0 omega^2 R^2 / 2 (R = r sin theta) supplies
//!       the centripetal force, so this is an EXACT steady state of the Euler equations.
//!       With an <mhd> block a uniform axial field B = b0 zhat is added: v x B =
//!       omega b0 R Rhat is curl-free and J = 0, so the MHD state is exactly steady too.
//!       v.rhat = 0, so nothing crosses the radial walls; ix1_bc/ox1_bc = user hold the
//!       exact state in the radial ghosts.  The deviation from the initial state at any
//!       time IS the error.  <hydro> or <mhd>.
//!
//!   11 = RESISTIVE DECAY.  A force-free Cartesian field B = b0 (sin(alpha z),
//!       cos(alpha z), 0) has curl B = alpha B, so J x B = 0 and |B|^2 = b0^2 is UNIFORM.
//!       With constant eta and v = 0 the exact resistive solution is
//!           B(t) = B(0) exp(-eta alpha^2 t),
//!           p(t) = p0 + (gamma-1) b0^2 (1 - exp(-2 eta alpha^2 t)) / 2,
//!       uniform in space and at rest for all time: the Ohmic heating eta J^2 =
//!       eta alpha^2 |B|^2 is uniform and the Poynting flux eta J x B vanishes, so the
//!       total energy density is exactly constant.  The face fields are built from the
//!       vector potential A = B/alpha by exact Stokes loops (Gauss-Legendre edge
//!       integrals), so div B is zero to round-off and the same construction, scaled by
//!       the decay factor, gives the exact solution at any time -- in the radial ghosts
//!       and in the final error measurement.  A gate that a curl on the polar row, the
//!       resistive dual mesh, or the resistive dt gets wrong by O(1), and a second-order
//!       convergence test for all three face fields and the pressure.  <mhd> with
//!       ohmic_resistivity = constant required.
//!
//!   12 = TOROIDAL FIELD (cs_test iprob=11's azimuthal part).  B = b0c R phihat with
//!       R = r sin(theta), i.e. b0c (-y, x, 0): curl B = 2 b0c zhat is a UNIFORM current,
//!       J x B = -2 b0c^2 R Rhat is balanced by p = p0 - b0c^2 R^2, and B_r = B_theta = 0
//!       everywhere, so nothing crosses the radial walls in the exact solution: the
//!       CLOSED
//!       conservation gate for MHD (mass, energy, L_z), which the force-free field cannot
//!       give because it threads the walls.  With constant eta the field is STILL static
//!       (curl of the uniform eta J vanishes) while the uniform Ohmic heating raises p by
//!       (gamma-1) 4 eta b0c^2 t everywhere; the total energy must rise by exactly
//!       eta J^2 V t, which the finalizer checks against the code's own energy sum, and
//!       the resistive EMF on every edge is checked against eta J.  <mhd>, resistivity
//!       optional.
//!
//! iprob 3, 11 and 12 report L1 / Linf errors against the exact state at the end of
//!       the run
//! (pgen_final_func), split into the POLAR rows and the interior, and append them to
//! <basename>-errs.dat; with problem/user_hist = true the same L1 errors are written to
//! the history file every dt.
//!
//! The face-normal components of the uniform field (iprob 3, 8) are FACE AVERAGES, not
//! point values: on an r-face the area-weighted mean of b0 cos(theta) is
//! b0 (cos th_l + cos th_r)/2, and on a theta-face -b0 sin(theta_f) is constant.  With
//! those the discrete divergence is zero EXACTLY (both terms reduce to
//! b0 (r_r^2 - r_l^2) dphi (sin^2 th_r - sin^2 th_l)/2), so no monopole force
//! contaminates the measurement.

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <string>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"
#include "outputs/outputs.hpp"
#include "pgen.hpp"
#include "hydro/hydro.hpp"

// FLUX PROBE (problem/flux_probe = 1, needs problem/user_srcs = true): on the first stage
// of the first cycle, print the code's momentum flux at one polar cell's outer theta-face
// and two radial faces next to the EXACT stress projection T.n of the analytic field,
// with each face's area/volume weight.  Tells which face carries the polar-row residual.
namespace {
Real sp_b0 = 1.0, sp_p0 = 1.0, sp_d0 = 1.0, sp_omega = 0.2, sp_alpha = 1.0;
Real sp_eta = 0.0, sp_bc_tfrac = 0.0, sp_svol = 1.0, sp_b0c = 0.3, sp_e0 = 0.0;
int sp_bdir = 2, sp_probe_done = 0, sp_iprob = 8, sp_axis = 2;

//----------------------------------------------------------------------------------------
// The exact states.

//! rigid rotation v = omega ahat x r about ahat = zhat (axis 2) or xhat (axis 0), at
//! (r, theta, phi): (v_theta, v_phi) and the pressure p0 + d0 omega^2 R^2 / 2 with R the
//! distance from the axis.  About xhat the flow runs THROUGH the poles (v_theta =
//! -omega r sin(phi), v_phi = -omega r cos(theta) cos(phi)) and p is not axisymmetric:
//! the polar-row test for a scalar gradient across the axis.
KOKKOS_INLINE_FUNCTION
void RigidRotState(Real r, Real th, Real ph, int axis, Real d0, Real p0, Real omega,
                   Real &vth, Real &vphi, Real &p) {
  Real rcyl2;
  if (axis == 0) {
    vth = -omega*r*sin(ph);
    vphi = -omega*r*cos(th)*cos(ph);
    rcyl2 = r*r*(SQR(sin(th)*sin(ph)) + SQR(cos(th)));
  } else {
    vth = 0.0;
    vphi = omega*r*sin(th);
    rcyl2 = SQR(r*sin(th));
  }
  p = p0 + 0.5*d0*omega*omega*rcyl2;
}

//! face averages of the uniform axial field b0 zhat: B_r over [thl, thr], B_th at thl
KOKKOS_INLINE_FUNCTION
void UniformZFaces(Real b0, Real thl, Real thr, Real &b1, Real &b2) {
  b1 = b0*0.5*(cos(thl) + cos(thr));
  b2 = -b0*sin(thl);
}

//! 4-point Gauss-Legendre nodes and weights on [-1, 1]
KOKKOS_INLINE_FUNCTION
void GL4(int n, Real &x, Real &w) {
  const Real xa = 0.3399810435848563, xb = 0.8611363115940526;
  const Real wa = 0.6521451548625461, wb = 0.3478548451374538;
  if (n == 0) { x = -xb; w = wb; } else if (n == 1) { x = -xa; w = wa; }
  else if (n == 2) { x = xa; w = wa; } else { x = xb; w = wb; }
}

//! spherical components of the force-free vector potential A = B/alpha,
//! B = b0 (sin(alpha z), cos(alpha z), 0), at (r, theta, phi)
KOKKOS_INLINE_FUNCTION
void FFVecPot(Real r, Real th, Real ph, Real b0, Real alpha,
              Real &ar, Real &at, Real &aph) {
  const Real z = r*cos(th);
  const Real ax = (b0/alpha)*sin(alpha*z);
  const Real ay = (b0/alpha)*cos(alpha*z);
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  ar  = ax*st*cp + ay*st*sp;
  at  = ax*ct*cp + ay*ct*sp;
  aph = -ax*sp + ay*cp;
}

//! spherical components of the force-free FIELD itself at a point (for the current check)
KOKKOS_INLINE_FUNCTION
void FFField(Real r, Real th, Real ph, Real b0, Real alpha,
             Real &br, Real &bt, Real &bph) {
  const Real z = r*cos(th);
  const Real bx = b0*sin(alpha*z), by = b0*cos(alpha*z);
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  br  = bx*st*cp + by*st*sp;
  bt  = bx*ct*cp + by*ct*sp;
  bph = -bx*sp + by*cp;
}

//! line integrals of A along the three edge types (exact to Gauss-Legendre 4)
KOKKOS_INLINE_FUNCTION
Real FFEdge1(Real rl, Real rr, Real th, Real ph, Real b0, Real alpha) {
  Real s = 0.0;
  for (int n=0; n<4; ++n) {
    Real x, w; GL4(n, x, w);
    const Real r = 0.5*(rl + rr) + 0.5*(rr - rl)*x;
    Real ar, at, aph; FFVecPot(r, th, ph, b0, alpha, ar, at, aph);
    s += w*ar;
  }
  return 0.5*(rr - rl)*s;
}
KOKKOS_INLINE_FUNCTION
Real FFEdge2(Real r, Real thl, Real thr, Real ph, Real b0, Real alpha) {
  Real s = 0.0;
  for (int n=0; n<4; ++n) {
    Real x, w; GL4(n, x, w);
    const Real th = 0.5*(thl + thr) + 0.5*(thr - thl)*x;
    Real ar, at, aph; FFVecPot(r, th, ph, b0, alpha, ar, at, aph);
    s += w*at*r;
  }
  return 0.5*(thr - thl)*s;
}
KOKKOS_INLINE_FUNCTION
Real FFEdge3(Real r, Real th, Real phl, Real phr, Real b0, Real alpha) {
  Real s = 0.0;
  for (int n=0; n<4; ++n) {
    Real x, w; GL4(n, x, w);
    const Real ph = 0.5*(phl + phr) + 0.5*(phr - phl)*x;
    Real ar, at, aph; FFVecPot(r, th, ph, b0, alpha, ar, at, aph);
    s += w*aph*r*sin(th);
  }
  return 0.5*(phr - phl)*s;
}

//! the three face-averaged field components from Stokes loops of A around each face.
//! Orientation: rhat = thhat x phhat, thhat = phhat x rhat, phhat = rhat x thhat.
KOKKOS_INLINE_FUNCTION
Real FFFaceB1(Real r, Real thl, Real thr, Real phl, Real phr, Real b0, Real alpha) {
  const Real flux = FFEdge2(r, thl, thr, phl, b0, alpha)
                  + FFEdge3(r, thr, phl, phr, b0, alpha)
                  - FFEdge2(r, thl, thr, phr, b0, alpha)
                  - FFEdge3(r, thl, phl, phr, b0, alpha);
  const Real area = r*r*fabs(cos(thl) - cos(thr))*(phr - phl);
  return (area > 0.0) ? flux/area : 0.0;
}
KOKKOS_INLINE_FUNCTION
Real FFFaceB2(Real rl, Real rr, Real th, Real phl, Real phr, Real b0, Real alpha) {
  const Real flux = FFEdge3(rl, th, phl, phr, b0, alpha)
                  + FFEdge1(rl, rr, th, phr, b0, alpha)
                  - FFEdge3(rr, th, phl, phr, b0, alpha)
                  - FFEdge1(rl, rr, th, phl, b0, alpha);
  const Real area = 0.5*(rr*rr - rl*rl)*fabs(sin(th))*(phr - phl);
  return (area > 0.0) ? flux/area : 0.0;
}
KOKKOS_INLINE_FUNCTION
Real FFFaceB3(Real rl, Real rr, Real thl, Real thr, Real ph, Real b0, Real alpha) {
  const Real flux = FFEdge1(rl, rr, thl, ph, b0, alpha)
                  + FFEdge2(rr, thl, thr, ph, b0, alpha)
                  - FFEdge1(rl, rr, thr, ph, b0, alpha)
                  - FFEdge2(rl, thl, thr, ph, b0, alpha);
  const Real area = 0.5*(rr*rr - rl*rl)*(thr - thl);
  return (area > 0.0) ? flux/area : 0.0;
}

//! the face ON the pole has zero area; bfield_bcs.cpp sets it to the mean of the two
//! adjacent theta-faces, the one across the pole entering with the sign the polar
//! exchange gives it (theta-hat reverses through the axis).  The exact counterpart is
//! half the difference of the adjacent face at phi and at phi + pi, with thn the angle
//! of that adjacent face (dtheta at the north pole, pi - dtheta at the south).
KOKKOS_INLINE_FUNCTION
Real FFPoleFaceB2(Real rl, Real rr, Real thn, Real phl, Real phr, Real b0, Real alpha) {
  return 0.5*(FFFaceB2(rl, rr, thn, phl, phr, b0, alpha)
              - FFFaceB2(rl, rr, thn, phl + M_PI, phr + M_PI, b0, alpha));
}

//! iprob 12: the phi-face average of B_phi = b0c r sin(theta) (area element r dr dtheta)
KOKKOS_INLINE_FUNCTION
Real TorFaceB3(Real rl, Real rr, Real tl, Real tr, Real b0c) {
  return b0c*((rr*rr*rr - rl*rl*rl)/3.0)*(cos(tl) - cos(tr))
         /(0.5*(rr*rr - rl*rl)*(tr - tl));
}
//! iprob 12: the exact pressure p0 - b0c^2 R^2 + (gamma-1) eta J^2 t, J = 2 b0c
KOKKOS_INLINE_FUNCTION
Real TorPressure(Real r, Real th, Real t, Real p0, Real b0c, Real eta, Real gm1) {
  return p0 - b0c*b0c*SQR(r*sin(th)) + gm1*eta*4.0*b0c*b0c*t;
}

//! exact decay factor and pressure of iprob = 11 at time t
KOKKOS_INLINE_FUNCTION
Real FFDecay(Real eta, Real alpha, Real t) { return exp(-eta*alpha*alpha*t); }
KOKKOS_INLINE_FUNCTION
Real FFPressure(Real p0, Real b0, Real gm1, Real decay) {
  return p0 + 0.5*gm1*b0*b0*(1.0 - decay*decay);
}
}  // namespace

void SPTestFluxProbe(Mesh *pm, const Real bdt);
void SPTestRadialBC(Mesh *pm);
void SPTestHistory(HistoryData *pdata, Mesh *pm);
void SPTestErrors(ParameterInput *pin, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  const int iprob = pin->GetOrAddInteger("problem", "iprob", 8);
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (!pmy_mesh_->use_spherical_polar
      || (pmbp->pmhd == nullptr && pmbp->phydro == nullptr) || (iprob != 3 && !is_mhd)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "sp_test requires mesh/use_spherical_polar = true and an <mhd> block"
              << " (iprob = 3 also accepts <hydro>)" << std::endl;
    exit(EXIT_FAILURE);
  }
  const Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  const Real b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  const Real omega = pin->GetOrAddReal("problem", "omega", 0.2);
  const Real alpha = pin->GetOrAddReal("problem", "alpha", 0.5*M_PI);
  const Real b0c = pin->GetOrAddReal("problem", "b0c", 0.3);   // iprob 12 amplitude
  sp_b0c = b0c;
  // field direction: bdir = 2 (default) is b0 zhat, bdir = 0 is b0 xhat.  For xhat the
  // face averages are B_r = b0 <sin th> cos-average, B_th = b0 cos(th_f) <cos ph>,
  // B_ph = -b0 sin(ph_f); each is the exact area-weighted mean over its face.
  const int bdir = pin->GetOrAddInteger("problem", "bdir", 2);
  // iprob 3: rotation axis, 2 = zhat (default), 0 = xhat (flow through the poles; HYDRO
  // only: with a field it would need B along xhat, which the MHD branch does not build)
  const int axis = pin->GetOrAddInteger("problem", "rot_axis", 2);
  if (iprob == 3 && axis != 2 && is_mhd) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "sp_test iprob = 3 with rot_axis != 2 is hydro only" << std::endl;
    exit(EXIT_FAILURE);
  }
  sp_axis = axis;
  sp_b0 = b0; sp_p0 = p0; sp_d0 = d0; sp_bdir = bdir; sp_iprob = iprob;
  sp_omega = omega; sp_alpha = alpha;
  // the time the radial ghosts of iprob = 11 are evaluated at is t^n + bc_time_frac*dt;
  // Mesh::time advances only at the end of a cycle while the BCs run inside each stage
  sp_bc_tfrac = pin->GetOrAddReal("problem", "bc_time_frac", 0.0);
  if (pin->GetOrAddInteger("problem", "flux_probe", 0) != 0) {
    user_srcs_func = SPTestFluxProbe;
  }
  // econsistent = true: the magnetic energy is built from the SAME cell-centred bcc that
  // ConsToPrim subtracts, so the recovered pressure is exactly uniform and only the
  // geometric-source/flux-divergence mismatch survives.  false (cs_test's choice): the
  // analytic b0^2/2, so the |bcc|^2 - b0^2 defect appears as an O(dtheta^2) pressure
  // ripple, largest in the polar rows -- which then dominates the one-step force.
  const bool econs = pin->GetOrAddBoolean("problem", "econsistent", true);
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  if (iprob == 11) {
    if (pmbp->pmhd->presist == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "sp_test iprob = 11 needs <mhd>/ohmic_resistivity = constant"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    sp_eta = pmbp->pmhd->presist->eta_ohm_const;
  }
  if (iprob == 12) {
    sp_eta = (pmbp->pmhd->presist != nullptr) ? pmbp->pmhd->presist->eta_ohm_const : 0.0;
  }
  if (iprob == 3 || iprob == 11 || iprob == 12) {
    user_bcs_func = SPTestRadialBC;
    user_hist_func = SPTestHistory;
    pgen_final_func = SPTestErrors;
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  const int ng = indcs.ng;
  const int nmb1 = pmbp->nmb_thispack - 1;
  // whole extent INCLUDING ghosts, so the state is consistent before the first BC call
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto &x2f_ = pmbp->pcoord->xx2f;   // theta at the theta-faces, (m, j)
  auto &x3f_ = pmbp->pcoord->xx3f;   // phi at the phi-faces, (m, k)
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x1v_ = pmbp->pcoord->x1v;    // volume centroids
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &mb_bcs = pmbp->pmb->mb_bcs;

  // total active volume, for the volume-weighted L1 norms of the history and the errors
  {
    auto &vol_ = pmbp->pcoord->volume;
    Real svol = 0.0;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
    const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
    Kokkos::parallel_reduce("sp_test_vol", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1+1)*nkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji;
      const int j = (idx - m*nkji - k*nji)/nx1;
      const int i = idx - m*nkji - k*nji - j*nx1;
      sum += vol_(m, k+ks, j+js, i+is);
    }, Kokkos::Sum<Real>(svol));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &svol, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    sp_svol = svol;
  }

  if (is_mhd) {
    auto &b0f = pmbp->pmhd->b0;
    par_for("sp_test_b", DevExeSpace(), 0,nmb1, 0,n3, 0,n2, 0,n1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const int jr = (j < n2) ? j+1 : j, kr = (k < n3) ? k+1 : k;
      const Real tl = x2f_(m,j), tr = x2f_(m,jr);
      const Real pl = x3f_(m,k), pr = x3f_(m,kr);
      if (iprob == 11) {
        const Real rl = x1f_(m,i), rr = x1f_(m,(i < n1) ? i+1 : i);
        if (j < n2 && k < n3) {
          b0f.x1f(m,k,j,i) = FFFaceB1(rl, tl, tr, pl, pr, b0, alpha);
        }
        if (i < n1 && k < n3) {
          const bool npole = (j == js) &&
              (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar);
          const bool spole = (j == je+1) &&
              (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
          if (npole) {
            b0f.x2f(m,k,j,i) = FFPoleFaceB2(rl, rr, x2f_(m,j+1), pl, pr, b0, alpha);
          } else if (spole) {
            b0f.x2f(m,k,j,i) = FFPoleFaceB2(rl, rr, x2f_(m,j-1), pl, pr, b0, alpha);
          } else {
            b0f.x2f(m,k,j,i) = FFFaceB2(rl, rr, tl, pl, pr, b0, alpha);
          }
        }
        if (i < n1 && j < n2) {
          b0f.x3f(m,k,j,i) = FFFaceB3(rl, rr, tl, tr, pl, b0, alpha);
        }
      } else if (iprob == 12) {
        const Real rl = x1f_(m,i), rr = x1f_(m,(i < n1) ? i+1 : i);
        if (j < n2 && k < n3) { b0f.x1f(m,k,j,i) = 0.0; }
        if (i < n1 && k < n3) { b0f.x2f(m,k,j,i) = 0.0; }
        if (i < n1 && j < n2) { b0f.x3f(m,k,j,i) = TorFaceB3(rl, rr, tl, tr, b0c); }
      } else if (bdir == 2 || iprob == 3) {
        Real b1, b2;
        UniformZFaces(b0, tl, tr, b1, b2);
        if (j < n2 && k < n3) { b0f.x1f(m,k,j,i) = b1; }
        if (i < n1 && k < n3) { b0f.x2f(m,k,j,i) = b2; }
        if (i < n1 && j < n2) { b0f.x3f(m,k,j,i) = 0.0; }
      } else {
        // B = b0 xhat: B_r = b0 sin(th) cos(ph), B_th = b0 cos(th) cos(ph),
        // B_ph = -b0 sin(ph)
        if (j < n2 && k < n3) {
          // <sin th>_area = Int sin^2 / Int sin ; <cos ph> = (sin pr - sin pl)/(pr - pl)
          const Real isin2 = 0.5*(tr - tl) - 0.25*(sin(2.0*tr) - sin(2.0*tl));
          const Real isin  = cos(tl) - cos(tr);
          b0f.x1f(m,k,j,i) = b0*(isin2/isin)*(sin(pr) - sin(pl))/(pr - pl);
        }
        if (i < n1 && k < n3) {
          b0f.x2f(m,k,j,i) = b0*cos(tl)*(sin(pr) - sin(pl))/(pr - pl);
        }
        if (i < n1 && j < n2) {
          b0f.x3f(m,k,j,i) = -b0*sin(pl);
        }
      }
    });
  }

  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  DvceArray5D<Real> bcc;
  DvceArray4D<Real> bx1f, bx2f, bx3f;
  if (is_mhd) {
    bcc = pmbp->pmhd->bcc0;
    bx1f = pmbp->pmhd->b0.x1f; bx2f = pmbp->pmhd->b0.x2f; bx3f = pmbp->pmhd->b0.x3f;
  }
  par_for("sp_test_u", DevExeSpace(), 0,nmb1, 0,n3-1, 0,n2-1, 0,n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real r = x1v_(m,i), th = x2v_(m,j), ph = x3v_(m,k);
    Real dn = d0, pgas = p0, vth = 0.0, vphi = 0.0;
    if (iprob == 3) RigidRotState(r, th, ph, axis, d0, p0, omega, vth, vphi, pgas);
    if (iprob == 12) pgas = TorPressure(r, th, 0.0, p0, b0c, 0.0, gm1);
    u0(m,IDN,k,j,i) = dn;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = dn*vth;
    u0(m,IM3,k,j,i) = dn*vphi;
    Real bsq = 0.0;
    if (is_mhd) {
      // the SAME position-weighted interpolation to the volume centroid that ConsToPrim
      // uses on this grid (ideal_mhd.cpp), so econsistent really is consistent
      Real lw, rw;
      lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
      rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
      bcc(m,IBX,k,j,i) = lw*bx1f(m,k,j,i) + rw*bx1f(m,k,j,i+1);
      lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
      rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
      bcc(m,IBY,k,j,i) = lw*bx2f(m,k,j,i) + rw*bx2f(m,k,j+1,i);
      lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
      rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
      bcc(m,IBZ,k,j,i) = lw*bx3f(m,k,j,i) + rw*bx3f(m,k+1,j,i);
      bsq = (econs || iprob != 8) ? (SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                                    + SQR(bcc(m,IBZ,k,j,i))) : b0*b0;
    }
    u0(m,IEN,k,j,i) = pgas/gm1 + 0.5*dn*(vth*vth + vphi*vphi) + 0.5*bsq;
  });
  // the initial total energy, for the iprob 12 Ohmic budget (same sum as history.cpp)
  {
    auto &vol_ = pmbp->pcoord->volume;
    Real e0 = 0.0;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
    const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
    Kokkos::parallel_reduce("sp_test_e0", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1+1)*nkji),
    KOKKOS_LAMBDA(const int &idx, Real &sum) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji;
      const int j = (idx - m*nkji - k*nji)/nx1;
      const int i = idx - m*nkji - k*nji - j*nx1;
      sum += vol_(m, k+ks, j+js, i+is)*u0(m, IEN, k+ks, j+js, i+is);
    }, Kokkos::Sum<Real>(e0));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &e0, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    sp_e0 = e0;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn SPTestRadialBC
//! \brief the exact state of iprob = 3 / 11 in the radial ghosts (ix1_bc/ox1_bc = user).
//! `user` is not one of the cases bfield_bcs.cpp handles, so the ghost FACES must be set
//! here too, and the cell-centred field and energy rebuilt from them.

void SPTestRadialBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  const int nhyd = is_mhd ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  const int nscal = is_mhd ? pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;
  const int iprob = sp_iprob;
  const Real d0 = sp_d0, p0 = sp_p0, omega = sp_omega, alpha = sp_alpha;
  const int axis = sp_axis;
  const Real b0c = sp_b0c, eta_ = sp_eta;
  const Real tbc = pm->time + sp_bc_tfrac*pm->dt;
  const Real decay = (iprob == 11) ? FFDecay(sp_eta, alpha, tbc) : 1.0;
  const Real b0 = sp_b0*decay;
  const Real p11 = FFPressure(p0, sp_b0, gm1, decay);

  // the ghost FACES first (MHD): face i on the inner side, face i+1 on the outer, and
  // every angular face of the ghost cell including the last j / k face
  if (is_mhd) {
    auto &b0f = pmbp->pmhd->b0;
    par_for("sp_test_rbc_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
            0,(n2-1), 0,(ng-1),
    KOKKOS_LAMBDA(int m, int k, int j, int ig) {
      for (int side=0; side<2; ++side) {
        if (side == 0 &&
            mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
        if (side == 1 &&
            mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
        const int i = (side == 0) ? (is-ng+ig) : (ie+1+ig);
        const int ifc = (side == 0) ? i : i+1;   // the face NOT shared with the interior
        const Real rl = x1f_(m,i), rr = x1f_(m,i+1);
        const Real tl = x2f_(m,j), tr = x2f_(m,j+1);
        const Real pl = x3f_(m,k), pr = x3f_(m,k+1);
        if (iprob == 11) {
          const bool npole = (j == js) &&
              (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar);
          const bool spole = (j == je+1) &&
              (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
          b0f.x1f(m,k,j,ifc) = FFFaceB1(x1f_(m,ifc), tl, tr, pl, pr, b0, alpha);
          if (npole) {
            b0f.x2f(m,k,j,i) = FFPoleFaceB2(rl, rr, tr, pl, pr, b0, alpha);
          } else if (spole) {
            b0f.x2f(m,k,j,i) = FFPoleFaceB2(rl, rr, x2f_(m,j-1), pl, pr, b0, alpha);
          } else {
            b0f.x2f(m,k,j,i) = FFFaceB2(rl, rr, tl, pl, pr, b0, alpha);
          }
          b0f.x3f(m,k,j,i) = FFFaceB3(rl, rr, tl, tr, pl, b0, alpha);
          if (j == n2-1) {
            b0f.x2f(m,k,j+1,i) = FFFaceB2(rl, rr, tr, pl, pr, b0, alpha);
          }
          if (k == n3-1) {
            b0f.x3f(m,k+1,j,i) = FFFaceB3(rl, rr, tl, tr, pr, b0, alpha);
          }
        } else if (iprob == 12) {
          b0f.x1f(m,k,j,ifc) = 0.0;
          b0f.x2f(m,k,j,i) = 0.0;
          b0f.x3f(m,k,j,i) = TorFaceB3(rl, rr, tl, tr, b0c);
          if (j == n2-1) { b0f.x2f(m,k,j+1,i) = 0.0; }
          if (k == n3-1) { b0f.x3f(m,k+1,j,i) = TorFaceB3(rl, rr, tl, tr, b0c); }
        } else {
          Real b1, b2;
          UniformZFaces(b0, tl, tr, b1, b2);
          b0f.x1f(m,k,j,ifc) = b1;
          b0f.x2f(m,k,j,i) = b2;
          b0f.x3f(m,k,j,i) = 0.0;
          if (j == n2-1) { b0f.x2f(m,k,j+1,i) = -b0*sin(tr); }
          if (k == n3-1) { b0f.x3f(m,k+1,j,i) = 0.0; }
        }
      }
    });
  }

  DvceArray5D<Real> bcc;
  DvceArray4D<Real> bx1f, bx2f, bx3f;
  if (is_mhd) {
    bcc = pmbp->pmhd->bcc0;
    bx1f = pmbp->pmhd->b0.x1f; bx2f = pmbp->pmhd->b0.x2f; bx3f = pmbp->pmhd->b0.x3f;
  }
  par_for("sp_test_rbc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0,(n3-1),
          0,(n2-1), 0,(ng-1),
  KOKKOS_LAMBDA(int m, int k, int j, int ig) {
    for (int side=0; side<2; ++side) {
      if (side == 0 &&
          mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) continue;
      if (side == 1 &&
          mb_bcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::user) continue;
      const int i = (side == 0) ? (is-ng+ig) : (ie+1+ig);
      const Real r = x1v_(m,i), th = x2v_(m,j), ph = x3v_(m,k);
      Real dn = d0, pgas = p0, vth = 0.0, vphi = 0.0;
      if (iprob == 3) {
        RigidRotState(r, th, ph, axis, d0, p0, omega, vth, vphi, pgas);
      } else if (iprob == 12) {
        pgas = TorPressure(r, th, tbc, p0, b0c, eta_, gm1);
      } else {
        pgas = p11;
      }
      Real bsq = 0.0;
      if (is_mhd) {
        Real lw, rw;
        lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        bcc(m,IBX,k,j,i) = lw*bx1f(m,k,j,i) + rw*bx1f(m,k,j,i+1);
        lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        bcc(m,IBY,k,j,i) = lw*bx2f(m,k,j,i) + rw*bx2f(m,k,j+1,i);
        lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        bcc(m,IBZ,k,j,i) = lw*bx3f(m,k,j,i) + rw*bx3f(m,k+1,j,i);
        bsq = SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) + SQR(bcc(m,IBZ,k,j,i));
      }
      w0(m,IDN,k,j,i) = dn;
      w0(m,IVX,k,j,i) = 0.0;
      w0(m,IVY,k,j,i) = vth;
      w0(m,IVZ,k,j,i) = vphi;
      w0(m,IEN,k,j,i) = pgas/gm1;
      u0(m,IDN,k,j,i) = dn;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = dn*vth;
      u0(m,IM3,k,j,i) = dn*vphi;
      u0(m,IEN,k,j,i) = pgas/gm1 + 0.5*dn*(vth*vth + vphi*vphi) + 0.5*bsq;
      for (int n=nhyd; n<nhyd+nscal; ++n) {
        w0(m,n,k,j,i) = 1.0;
        u0(m,n,k,j,i) = dn;
      }
    }
  });
  (void)n1;
}

//----------------------------------------------------------------------------------------
//! \fn SPTestHistory
//! \brief volume-weighted L1 errors against the exact state, per unit total volume:
//! |v - v_ex|, |p - p_ex|, |rho - d0| and, for MHD, |bcc - bcc_ex| (bcc_ex from the exact
//! faces at the cell centroid, the SAME interpolation ConsToPrim uses).

void SPTestHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 5;
  pdata->label[0] = "L1-v";
  pdata->label[1] = "L1-p";
  pdata->label[2] = "L1-rho";
  pdata->label[3] = "L1-b";
  // z angular momentum: the cell-average momentum times the EXACT cell integral of the
  // lever arm, Int r sin(theta) dV = (r_r^4 - r_l^4)/4 * Int sin^2 * dphi
  pdata->label[4] = "Lz";

  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &w0_ = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  DvceArray5D<Real> bcc;
  if (is_mhd) bcc = pmbp->pmhd->bcc0;
  auto &vol_ = pmbp->pcoord->volume;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, je = indcs.je;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  const int iprob = sp_iprob;
  const Real d0 = sp_d0, p0 = sp_p0, omega = sp_omega, alpha = sp_alpha;
  const int axis = sp_axis;
  const Real b0c = sp_b0c, eta_ = sp_eta;
  const Real decay = (iprob == 11) ? FFDecay(sp_eta, alpha, pm->time) : 1.0;
  const Real b0 = sp_b0*decay;
  const Real p11 = FFPressure(p0, sp_b0, gm1, decay);
  const Real isvol = 1.0/sp_svol;

  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("SPHist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    const Real r = x1v_(m,i), th = x2v_(m,j), ph = x3v_(m,k);
    Real pex = p0, vth = 0.0, vphi = 0.0;
    if (iprob == 3) {
      RigidRotState(r, th, ph, axis, d0, p0, omega, vth, vphi, pex);
    } else if (iprob == 12) {
      pex = TorPressure(r, th, pm->time, p0, b0c, eta_, gm1);
    } else {
      pex = p11;
    }
    const Real dv = sqrt(SQR(w0_(m,IVX,k,j,i)) + SQR(w0_(m,IVY,k,j,i) - vth)
                         + SQR(w0_(m,IVZ,k,j,i) - vphi));
    const Real dp = fabs(gm1*w0_(m,IEN,k,j,i) - pex);
    const Real dd = fabs(w0_(m,IDN,k,j,i) - d0);
    Real db = 0.0;
    if (is_mhd) {
      const Real rl = x1f_(m,i), rr = x1f_(m,i+1);
      const Real tl = x2f_(m,j), tr = x2f_(m,j+1);
      const Real pl = x3f_(m,k), pr = x3f_(m,k+1);
      Real f1l, f1r, f2l, f2r, f3l, f3r;
      if (iprob == 11) {
        const bool npole = (j == js) &&
            (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar);
        const bool spole = (j == je) &&
            (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
        f1l = FFFaceB1(rl, tl, tr, pl, pr, b0, alpha);
        f1r = FFFaceB1(rr, tl, tr, pl, pr, b0, alpha);
        f2l = npole ? FFPoleFaceB2(rl, rr, tr, pl, pr, b0, alpha)
                    : FFFaceB2(rl, rr, tl, pl, pr, b0, alpha);
        f2r = spole ? FFPoleFaceB2(rl, rr, tl, pl, pr, b0, alpha)
                    : FFFaceB2(rl, rr, tr, pl, pr, b0, alpha);
        f3l = FFFaceB3(rl, rr, tl, tr, pl, b0, alpha);
        f3r = FFFaceB3(rl, rr, tl, tr, pr, b0, alpha);
      } else if (iprob == 12) {
        f1l = 0.0; f1r = 0.0; f2l = 0.0; f2r = 0.0;
        f3l = TorFaceB3(rl, rr, tl, tr, b0c); f3r = f3l;
      } else {
        UniformZFaces(b0, tl, tr, f1l, f2l);
        f1r = f1l; f2r = -b0*sin(tr); f3l = 0.0; f3r = 0.0;
      }
      Real lw, rw;
      lw = (rr-r)/(rr-rl); rw = (r-rl)/(rr-rl);
      const Real bx = lw*f1l + rw*f1r;
      lw = (tr-th)/(tr-tl); rw = (th-tl)/(tr-tl);
      const Real by = lw*f2l + rw*f2r;
      const Real ph = x3v_(m,k);
      lw = (pr-ph)/(pr-pl); rw = (ph-pl)/(pr-pl);
      const Real bz = lw*f3l + rw*f3r;
      db = sqrt(SQR(bcc(m,IBX,k,j,i) - bx) + SQR(bcc(m,IBY,k,j,i) - by)
                + SQR(bcc(m,IBZ,k,j,i) - bz));
    }
    const Real w = vol_(m,k,j,i)*isvol;
    array_sum::GlobalSum hv;
    hv.the_array[0] = w*dv;
    hv.the_array[1] = w*dp;
    hv.the_array[2] = w*dd;
    hv.the_array[3] = w*db;
    {
      const Real rl = x1f_(m,i), rr = x1f_(m,i+1);
      const Real tl = x2f_(m,j), tr = x2f_(m,j+1);
      const Real isin2 = 0.5*(tr - tl) - 0.25*(sin(2.0*tr) - sin(2.0*tl));
      const Real lever = 0.25*(rr*rr*rr*rr - rl*rl*rl*rl)*isin2*(x3f_(m,k+1) - x3f_(m,k));
      hv.the_array[4] = w0_(m,IDN,k,j,i)*w0_(m,IVZ,k,j,i)*lever;
    }
    for (int n=5; n<NREDUCTION_VARIABLES; ++n) { hv.the_array[n] = 0.0; }
    mb_sum += hv;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  for (int n=0; n<pdata->nhist; ++n) { pdata->hdata[n] = sum_this_mb.the_array[n]; }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn SPTestErrors
//! \brief L1 (volume-weighted, per unit volume) and Linf errors against the exact state
//! at the final time, for the whole shell and split into the POLAR rows (cells within
//! problem/err_nband cells of either pole) and the interior; the face fields are compared
//! face by face against the exact face averages.  Printed and appended to
//! <basename>-errs.dat.

void SPTestErrors(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  auto &w0d = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto w0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), w0d);
  auto vol = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->volume);
  auto x1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto x2f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx2f);
  auto x3f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx3f);
  auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto x2v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
  auto x3v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
  DvceArray4D<Real> b1d, b2d, b3d;
  if (is_mhd) {
    b1d = pmbp->pmhd->b0.x1f; b2d = pmbp->pmhd->b0.x2f; b3d = pmbp->pmhd->b0.x3f;
  }
  auto b1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b1d);
  auto b2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b2d);
  auto b3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b3d);

  const int iprob = sp_iprob;
  const Real d0 = sp_d0, p0 = sp_p0, omega = sp_omega, alpha = sp_alpha;
  const int axis = sp_axis;
  const Real b0c = sp_b0c, eta_ = sp_eta;
  const Real decay = (iprob == 11) ? FFDecay(sp_eta, alpha, pm->time) : 1.0;
  const Real b0 = sp_b0*decay;
  const Real p11 = FFPressure(p0, sp_b0, gm1, decay);
  // polar rows: within nband cells of either pole, by ANGLE (independent of blocks)
  const int nband = pin->GetOrAddInteger("problem", "err_nband", 2);
  const Real dth_pole = nband*M_PI/static_cast<Real>(pm->mesh_indcs.nx2);

  // sums: [0..3] L1 v,p,rho,b (volume / face weighted), [4] volume, [5] face area
  // for polar (index P) and interior (index I); maxima: v, p, b
  enum { NS = 6 };
  Real sP[NS] = {}, sI[NS] = {};
  Real mxv = 0.0, mxp = 0.0, mxb = 0.0, mxv_pol = 0.0, mxb_pol = 0.0;
  int mxb_f = -1, mxb_j = -1, mxb_i = -1;   // face type (0-5), local j, i of max B error
  Real mxdivb = 0.0;
  const Real lshell = pm->mesh_size.x1max - pm->mesh_size.x1min;
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        const Real th = x2v(m,j);
        const bool polar = (th < dth_pole) || (th > M_PI - dth_pole);
        Real *s = polar ? sP : sI;
        for (int i=is; i<=ie; ++i) {
          const Real r = x1v(m,i), ph = x3v(m,k);
          Real pex = p0, vth = 0.0, vphi = 0.0;
          if (iprob == 3) {
            RigidRotState(r, th, ph, axis, d0, p0, omega, vth, vphi, pex);
          } else if (iprob == 12) {
            pex = TorPressure(r, th, pm->time, p0, b0c, eta_, gm1);
          } else {
            pex = p11;
          }
          const Real dv = std::sqrt(SQR(w0(m,IVX,k,j,i)) + SQR(w0(m,IVY,k,j,i) - vth)
                                    + SQR(w0(m,IVZ,k,j,i) - vphi));
          const Real dp = std::fabs(gm1*w0(m,IEN,k,j,i) - pex);
          const Real dd = std::fabs(w0(m,IDN,k,j,i) - d0);
          const Real w = vol(m,k,j,i);
          s[0] += w*dv; s[1] += w*dp; s[2] += w*dd; s[4] += w;
          mxv = std::fmax(mxv, dv); mxp = std::fmax(mxp, dp);
          if (polar) mxv_pol = std::fmax(mxv_pol, dv);
          if (is_mhd) {
            // the discrete divergence, in units of b0/L with L the shell thickness: the
            // construction is exact to round-off; anything else here says a face was
            // built with the wrong sign or a wrong area, which would masquerade as a
            // monopole force in the velocity error below
            {
              const Real rl = x1f(m,i), rr = x1f(m,i+1);
              const Real tl = x2f(m,j), tr = x2f(m,j+1);
              const Real pl = x3f(m,k), pr = x3f(m,k+1);
              const Real a1l = rl*rl*std::fabs(std::cos(tl) - std::cos(tr))*(pr - pl);
              const Real a1r = rr*rr*std::fabs(std::cos(tl) - std::cos(tr))*(pr - pl);
              const Real a2l = 0.5*(rr*rr - rl*rl)*std::fabs(std::sin(tl))*(pr - pl);
              const Real a2r = 0.5*(rr*rr - rl*rl)*std::fabs(std::sin(tr))*(pr - pl);
              const Real a3 = 0.5*(rr*rr - rl*rl)*(tr - tl);
              const Real divb = (a1r*b1(m,k,j,i+1) - a1l*b1(m,k,j,i)
                                 + a2r*b2(m,k,j+1,i) - a2l*b2(m,k,j,i)
                                 + a3*(b3(m,k+1,j,i) - b3(m,k,j,i)))/vol(m,k,j,i);
              mxdivb = std::fmax(mxdivb, std::fabs(divb)*lshell/sp_b0);
            }
            // the three faces on the low side of every cell, plus the high faces on the
            // last cell of each row: every face exactly once
            const Real rl = x1f(m,i), rr = x1f(m,i+1);
            const Real tl = x2f(m,j), tr = x2f(m,j+1);
            const Real pl = x3f(m,k), pr = x3f(m,k+1);
            for (int f=0; f<6; ++f) {
              if (f == 3 && i != ie) continue;
              if (f == 4 && j != je) continue;
              if (f == 5 && k != ke) continue;
              Real ex, code, area;
              if (f == 0 || f == 3) {
                const Real rf = (f == 0) ? rl : rr;
                if (iprob == 11) { ex = FFFaceB1(rf, tl, tr, pl, pr, b0, alpha); }
                else if (iprob == 12) { ex = 0.0; }
                else { Real bb; UniformZFaces(b0, tl, tr, ex, bb); }
                code = b1(m,k,j,(f == 0) ? i : i+1);
                area = rf*rf*std::fabs(std::cos(tl) - std::cos(tr))*(pr - pl);
              } else if (f == 1 || f == 4) {
                const Real tf = (f == 1) ? tl : tr;
                if (iprob == 11) { ex = FFFaceB2(rl, rr, tf, pl, pr, b0, alpha); }
                else if (iprob == 12) { ex = 0.0; }
                else { ex = -b0*std::sin(tf); }
                code = b2(m,k,(f == 1) ? j : j+1,i);
                area = 0.5*(rr*rr - rl*rl)*std::fabs(std::sin(tf))*(pr - pl);
              } else {
                const Real pf = (f == 2) ? pl : pr;
                ex = (iprob == 11) ? FFFaceB3(rl, rr, tl, tr, pf, b0, alpha)
                   : ((iprob == 12) ? TorFaceB3(rl, rr, tl, tr, b0c) : 0.0);
                code = b3(m,(f == 2) ? k : k+1,j,i);
                area = 0.5*(rr*rr - rl*rl)*(tr - tl);
              }
              // skip the degenerate theta-face ON the pole (sin(pi) is 1e-16, not 0)
              if (area <= 1.0e-12*rl*rl) continue;
              const Real db = std::fabs(code - ex);
              s[3] += area*db; s[5] += area;
              if (db > mxb) { mxb = db; mxb_f = f; mxb_j = j - js; mxb_i = i - is; }
              if (polar) mxb_pol = std::fmax(mxb_pol, db);
            }
          }
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, sP, NS, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, sI, NS, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  Real mx[6] = {mxv, mxp, mxb, mxv_pol, mxb_pol, mxdivb};
  MPI_Allreduce(MPI_IN_PLACE, mx, 6, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  mxv = mx[0]; mxp = mx[1]; mxb = mx[2]; mxv_pol = mx[3]; mxb_pol = mx[4]; mxdivb = mx[5];
#endif
  if (global_variable::my_rank != 0) return;

  Real sT[NS];
  for (int n=0; n<NS; ++n) sT[n] = sP[n] + sI[n];
  auto l1 = [](const Real *s, int n) {
    const Real w = (n == 3) ? s[5] : s[4];
    return (w > 0.0) ? s[n]/w : 0.0;
  };
  std::printf("### SP TEST ERRORS iprob=%d  t=%.6e  nx1=%d nx2=%d nx3=%d  polar band %d"
              " rows\n", iprob, pm->time, pm->mesh_indcs.nx1, pm->mesh_indcs.nx2,
              pm->mesh_indcs.nx3, nband);
  if (iprob == 11) {
    std::printf("###   exact decay factor %.6e  (eta=%g alpha=%g)  exact p=%.9e\n",
                decay, sp_eta, alpha, p11);
  }
  std::printf("###   L1 (per unit volume):  v %.6e  p %.6e  rho %.6e  B(faces) %.6e\n",
              l1(sT,0), l1(sT,1), l1(sT,2), l1(sT,3));
  std::printf("###   Linf:                  v %.6e  p %.6e  B %.6e  (B max on this rank"
              " at face type %d, j-js=%d, i-is=%d)\n", mxv, mxp, mxb, mxb_f, mxb_j, mxb_i);
  std::printf("###   POLAR rows:   L1 v %.6e  p %.6e  B %.6e   Linf v %.6e  B %.6e"
              "  (%.2f%% of the volume)\n", l1(sP,0), l1(sP,1), l1(sP,3), mxv_pol, mxb_pol,
              100.0*sP[4]/sT[4]);
  std::printf("###   INTERIOR:     L1 v %.6e  p %.6e  B %.6e\n", l1(sI,0), l1(sI,1),
              l1(sI,3));
  if (is_mhd) {
    std::printf("###   max |div B| L/b0 = %.3e  (round-off expected: the faces come from"
                " Stokes loops)\n", mxdivb);
  }

  // iprob 11: the code's resistive EMF eta*J on every edge against the exact
  // eta*alpha*B (curl B = alpha B), in units of eta*alpha*b0*decay, polar rows vs
  // interior, with the north and south maxima
  if ((iprob == 11 || iprob == 12) && is_mhd && pmbp->pmhd->presist != nullptr) {
    auto pres = pmbp->pmhd->presist;
    auto e1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x1e);
    auto e2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x2e);
    auto e3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pres->efld_resist.x3e);
    const Real eta = pres->eta_ohm_const;
    const Real scale = (iprob == 12) ? eta*2.0*b0c : eta*alpha*b0;
    // exact J on an edge: alpha*B (iprob 11) or 2 b0c zhat (iprob 12), then eta*J.t
    auto exactJ = [&](Real r, Real th, Real ph, Real &jr, Real &jt, Real &jp) {
      if (iprob == 12) {
        jr = 2.0*b0c*std::cos(th); jt = -2.0*b0c*std::sin(th); jp = 0.0;
      } else {
        FFField(r, th, ph, b0, alpha, jr, jt, jp);
        jr *= alpha; jt *= alpha; jp *= alpha;
      }
    };
    Real l1[3][2] = {}, mx[3][2] = {};
    std::int64_t nc[3][2] = {};
    int mxj[3] = {-1, -1, -1}, mxi[3] = {-1, -1, -1};
    Real mxn[3] = {0.0, 0.0, 0.0}, mxs[3] = {0.0, 0.0, 0.0};   // north / south pole rows
    const int nbp = nband;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            const Real rl = x1f(m,i), rr = x1f(m,(i <= ie) ? i+1 : i);
            const Real tl = x2f(m,j), tr = x2f(m,(j <= je) ? j+1 : j);
            const Real pl = x3f(m,k), pr = x3f(m,(k <= ke) ? k+1 : k);
            Real br, bt, bph;
            // x1 edge: along r at (theta_f(j), phi_f(k)), cells i <= ie
            if (i <= ie) {
              exactJ(0.5*(rl + rr), tl, pl, br, bt, bph);
              const bool pol = (j - js < nbp) || (je + 1 - j < nbp);
              const Real d = std::fabs(e1(m,k,j,i) - eta*br)/scale;
              l1[0][pol] += d; nc[0][pol]++;
              if (d > mx[0][pol]) {
                mx[0][pol] = d;
                if (pol) { mxj[0] = j - js; mxi[0] = i - is; }
              }
              if ((j - js < nbp)) mxn[0] = std::fmax(mxn[0], d);
              if ((je + 1 - j < nbp)) mxs[0] = std::fmax(mxs[0], d);
            }
            // x2 edge: along theta at (r_f(i), phi_f(k)), cells j <= je
            if (j <= je) {
              exactJ(rl, 0.5*(tl + tr), pl, br, bt, bph);
              const bool pol = (j - js < nbp) || (je - j < nbp);
              const Real d = std::fabs(e2(m,k,j,i) - eta*bt)/scale;
              l1[1][pol] += d; nc[1][pol]++;
              if (d > mx[1][pol]) {
                mx[1][pol] = d;
                if (pol) { mxj[1] = j - js; mxi[1] = i - is; }
              }
              if ((j - js < nbp)) mxn[1] = std::fmax(mxn[1], d);
              if ((je - j < nbp)) mxs[1] = std::fmax(mxs[1], d);
            }
            // x3 edge: along phi at (r_f(i), theta_f(j)), cells k <= ke; skip the pole
            if (k <= ke && std::fabs(std::sin(tl)) > 1.0e-12) {
              exactJ(rl, tl, 0.5*(pl + pr), br, bt, bph);
              const bool pol = (j - js < nbp) || (je + 1 - j < nbp);
              const Real d = std::fabs(e3(m,k,j,i) - eta*bph)/scale;
              l1[2][pol] += d; nc[2][pol]++;
              if (d > mx[2][pol]) {
                mx[2][pol] = d;
                if (pol) { mxj[2] = j - js; mxi[2] = i - is; }
              }
              if ((j - js < nbp)) mxn[2] = std::fmax(mxn[2], d);
              if ((je + 1 - j < nbp)) mxs[2] = std::fmax(mxs[2], d);
            }
          }
        }
      }
    }
    if (global_variable::my_rank == 0) {
      const char *nm[3] = {"E1 (r edges)   ", "E2 (theta edges)", "E3 (phi edges) "};
      std::printf("###   RESISTIVE EMF vs exact eta*J, in units eta*|J| (this rank):\n");
      for (int c=0; c<3; ++c) {
        std::printf("###     %s  interior L1 %.3e Linf %.3e | polar L1 %.3e Linf %.3e"
                    " (max at j-js=%d i-is=%d; north %.3e south %.3e)\n", nm[c],
                    (nc[c][0] > 0) ? l1[c][0]/nc[c][0] : 0.0, mx[c][0],
                    (nc[c][1] > 0) ? l1[c][1]/nc[c][1] : 0.0, mx[c][1], mxj[c], mxi[c],
                    mxn[c], mxs[c]);
      }
    }
  }

  // iprob 12: the Ohmic energy budget.  The exact solution gains eta J^2 V t; compare
  // with the code's own total-energy sum (the same sum history.cpp forms).
  if (iprob == 12 && is_mhd) {
    auto u0h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmhd->u0);
    Real etot = 0.0;
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) { etot += vol(m,k,j,i)*u0h(m,IEN,k,j,i); }
        }
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &etot, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const Real gain = eta_*4.0*b0c*b0c*sp_svol*pm->time;
    if (global_variable::my_rank == 0) {
      std::printf("###   OHMIC BUDGET: E(0) %.10e  E(t) %.10e  gained %.6e  exact"
                  " eta J^2 V t %.6e  ratio %.6f  (eta=0: dE/E = %.3e)\n", sp_e0, etot,
                  etot - sp_e0, gain, (gain != 0.0) ? (etot - sp_e0)/gain : 0.0,
                  (etot - sp_e0)/sp_e0);
    }
  }

  std::string fname = pin->GetString("job", "basename") + "-errs.dat";
  FILE *pfile = std::fopen(fname.c_str(), "r");
  const bool exists = (pfile != nullptr);
  if (exists) std::fclose(pfile);
  pfile = std::fopen(fname.c_str(), "a");
  if (pfile == nullptr) {
    std::cout << "### WARNING: could not open " << fname << std::endl;
    return;
  }
  if (!exists) {
    std::fprintf(pfile, "# iprob nx1 nx2 nx3 time  L1v L1p L1rho L1b  Linfv Linfp Linfb"
                        "  L1v_polar L1b_polar Linfv_polar Linfb_polar  L1v_int L1b_int\n");
  }
  std::fprintf(pfile, "%d %d %d %d %.6e  %.6e %.6e %.6e %.6e  %.6e %.6e %.6e"
                      "  %.6e %.6e %.6e %.6e  %.6e %.6e\n",
               iprob, pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3, pm->time,
               l1(sT,0), l1(sT,1), l1(sT,2), l1(sT,3), mxv, mxp, mxb,
               l1(sP,0), l1(sP,3), mxv_pol, mxb_pol, l1(sI,0), l1(sI,3));
  std::fclose(pfile);
}

//----------------------------------------------------------------------------------------
//! \fn SPTestFluxProbe

void SPTestFluxProbe(Mesh *pm, const Real bdt) {
  if (sp_probe_done++ != 0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  auto x1f = Kokkos::create_mirror_view(pmbp->pcoord->xx1f); Kokkos::deep_copy(x1f, pmbp->pcoord->xx1f);
  auto x2f = Kokkos::create_mirror_view(pmbp->pcoord->xx2f); Kokkos::deep_copy(x2f, pmbp->pcoord->xx2f);
  auto x2v = Kokkos::create_mirror_view(pmbp->pcoord->x2v);  Kokkos::deep_copy(x2v, pmbp->pcoord->x2v);
  auto a1 = Kokkos::create_mirror_view(pmbp->pcoord->area.x1f); Kokkos::deep_copy(a1, pmbp->pcoord->area.x1f);
  auto a2 = Kokkos::create_mirror_view(pmbp->pcoord->area.x2f); Kokkos::deep_copy(a2, pmbp->pcoord->area.x2f);
  auto vol = Kokkos::create_mirror_view(pmbp->pcoord->volume); Kokkos::deep_copy(vol, pmbp->pcoord->volume);
  auto f1 = Kokkos::create_mirror_view(pmbp->pmhd->uflx.x1f); Kokkos::deep_copy(f1, pmbp->pmhd->uflx.x1f);
  auto f2 = Kokkos::create_mirror_view(pmbp->pmhd->uflx.x2f); Kokkos::deep_copy(f2, pmbp->pmhd->uflx.x2f);
  auto bcc = Kokkos::create_mirror_view(pmbp->pmhd->bcc0); Kokkos::deep_copy(bcc, pmbp->pmhd->bcc0);
  auto b2f = Kokkos::create_mirror_view(pmbp->pmhd->b0.x2f); Kokkos::deep_copy(b2f, pmbp->pmhd->b0.x2f);
  auto w = Kokkos::create_mirror_view(pmbp->pmhd->w0); Kokkos::deep_copy(w, pmbp->pmhd->w0);
  const int m = 0, k = ks, i = is + 3;
  const Real b0 = sp_b0, p0 = sp_p0;
  std::printf("### SP FLUX PROBE (bdir=%d, b0=%g, p0=%g) block %d k=%d i=%d, cell j=js\n",
              sp_bdir, b0, p0, m, k, i);
  std::printf("###   V=%.6e  theta faces: %.6e %.6e  x2v=%.6e\n", vol(m,k,js,i),
              x2f(m,js), x2f(m,js+1), x2v(m,js));
  std::printf("###   cell bcc: Br=%.8e Bth=%.8e Bph=%.8e  p(w)=%.8e\n", bcc(m,IBX,k,js,i),
              bcc(m,IBY,k,js,i), bcc(m,IBZ,k,js,i), w(m,IEN,k,js,i)*(5.0/3.0-1.0));
  // exact stress projections for B = b0 zhat, v = 0:
  //   theta-face (n = thhat):  F_r = -B_th B_r = b0^2 sin cos,  F_th = p + b0^2 cos(2th)/2, F_ph = 0
  //   r-face     (n = rhat):   F_r = p + b0^2 (sin^2 - cos^2)/2 = p - b0^2 cos(2th)/2,  F_th = -B_r B_th = b0^2 sin cos
  for (int jf = js; jf <= js + 2; ++jf) {
    const Real th = x2f(m,jf);
    const Real exr = b0*b0*std::sin(th)*std::cos(th), exth = p0 + 0.5*b0*b0*std::cos(2.0*th);
    std::printf("###   THETA face j=%d (th=%.6f) A/V=%.6e  B_n(face)=%.8e\n"
                "###       code F_r=%.8e F_th=%.8e F_ph=%.3e | exact F_r=%.8e F_th=%.8e\n"
                "###       (code-exact)*A/V: r %.3e  th %.3e   in units (b0^2/2)/r: r %.3e th %.3e\n",
                jf, th, a2(m,k,jf,i)/vol(m,k,js,i), b2f(m,k,jf,i),
                f2(m,IM1,k,jf,i), f2(m,IM2,k,jf,i), f2(m,IM3,k,jf,i), exr, exth,
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,js,i),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,js,i),
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)));
  }
  // r-faces of the polar cell: the exact flux is the AREA MEAN over the face's theta range
  const Real tl = x2f(m,js), tr = x2f(m,js+1);
  const Real isin = std::cos(tl) - std::cos(tr);
  // <cos 2th> over the r-face with weight sin th: Int (2cos^2-1) sin dth / Int sin dth
  const Real icos2sin = (std::cos(tl)*std::cos(tl)*std::cos(tl) - std::cos(tr)*std::cos(tr)*std::cos(tr))/3.0;
  const Real avg_cos2 = 2.0*icos2sin/isin - 1.0;
  // <sin cos> = Int sin^2 cos / Int sin = (sin^3 tr - sin^3 tl)/3 / isin
  const Real avg_sc = (std::pow(std::sin(tr),3) - std::pow(std::sin(tl),3))/3.0/isin;
  for (int ii = i; ii <= i+1; ++ii) {
    const Real exr = p0 - 0.5*b0*b0*avg_cos2, exth = b0*b0*avg_sc;
    std::printf("###   R face i=%d (r=%.4f) A/V=%.6e\n"
                "###       code F_r=%.8e F_th=%.8e | exact(area-mean) F_r=%.8e F_th=%.8e\n"
                "###       (code-exact)*A/V in (b0^2/2)/r: r %.3e th %.3e\n",
                ii, x1f(m,ii), a1(m,k,js,ii)/vol(m,k,js,i),
                f1(m,IM1,k,js,ii), f1(m,IM2,k,js,ii), exr, exth,
                (f1(m,IM1,k,js,ii)-exr)*a1(m,k,js,ii)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)),
                (f1(m,IM2,k,js,ii)-exth)*a1(m,k,js,ii)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)));
  }
  // and the interior cell j = js+4 for reference (its outer theta-face)
  {
    const int jc = js + 4, jf = jc + 1; const Real th = x2f(m,jf);
    const Real exr = b0*b0*std::sin(th)*std::cos(th), exth = p0 + 0.5*b0*b0*std::cos(2.0*th);
    std::printf("###   INTERIOR ref: theta face j=%d (th=%.6f) A/V=%.6e  (code-exact)*A/V in (b0^2/2)/r: r %.3e th %.3e\n",
                jf, th, a2(m,k,jf,i)/vol(m,k,jc,i),
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,jc,i)/(0.5*b0*b0/x1f(m,i)),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,jc,i)/(0.5*b0*b0/x1f(m,i)));
  }
}
