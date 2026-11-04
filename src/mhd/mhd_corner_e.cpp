//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_corner_e.cpp
//  \brief
//  Also includes contributions to electric field from "source terms" such as the
//  shearing box.

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "diffusion/resistivity.hpp"
#include "mhd.hpp"
#include "reconstruct/plm.hpp"

#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn Real EmfSmoothW(fd, eps, rhoe)
//! \brief lhlld2-1102, <mhd>/emf_smooth = eps (read only when named; default off =
//! the GS07 sign switch, bitwise).  Weight of the LEFT (lower-index) upwind choice in the
//! GS07 corner EMF, as a smooth function of the face mass flux fd:
//!   w = 0.5*(1 + fd/sqrt(fd^2 + delta^2)),  delta = eps*sqrt(rho*e)  (~ eps*rho*c_s),
//! instead of the step function w = (fd >= 0).  The step makes the semi-discrete
//! right-hand side DISCONTINUOUS wherever a face mass flux changes sign (any flow at
//! rest or turning, and every low-Mach flow, where fd is round-off-sized noise about 0):
//! the two upwind choices differ by O(dx^2 d^2E/dx^2), so a Runge-Kutta stage that lands
//! on the other side of the switch than its predecessor costs an O(dt * jump) local
//! error, and the global temporal error is FIRST order at fixed grid (in-plane-B
//! magnetized Gresho vortex, M = 0.01: B order 0.8-1.3, v 1.2).  eps = 0 is the GS07
//! rule for a zero mass flux (the average of both sides) applied to every face.
//! For |fd| >> delta the weight is 1 or 0 to O(delta^2/fd^2): GS07 is recovered in any
//! flow that is not slow compared with eps*c_s.

KOKKOS_INLINE_FUNCTION
Real EmfSmoothW(const Real fd, const Real eps, const Real rhoe) {
  const Real del2 = eps*eps*rhoe;
  if (del2 <= 0.0) {   // eps = 0: GS07 zero-mass-flux rule everywhere
    return 0.5;
  }
  return 0.5*(1.0 + fd/sqrt(fd*fd + del2));
}

//! \brief lhlld2-1102: the GS07 upwind choice with smooth weights: a is the choice for
//! fd >= 0 (its reference EMF in cell A), b the one for fd < 0 (cell B); rA, rB = rho*e
//! (rho*c_iso^2 isothermal) of those two cells.
KOKKOS_INLINE_FUNCTION
Real EmfPick(const Real fd, const Real eps, const Real ra, const Real rb,
             const Real a, const Real b) {
  return b + EmfSmoothW(fd, eps, 0.5*(ra + rb))*(a - b);
}

//----------------------------------------------------------------------------------------
//! <mhd>/emf_scheme = uct | uct_hlld (lhlld2-1102, ucthlld-1102): the UCT corner EMF of
//! Mignone & Del Zanna 2021 (JCP 424, 109748, "MD21", eq. 35 with eq. 36 face averages),
//! evaluated per corner inside the GS07 corner kernels of MHD::CornerE (same loops, same
//! wall handling; the resistive/RKG EMFs, the cubed-sphere seam sum/average and the
//! vertex average of flux_correct_fc act on its result exactly as on GS07's).  For an
//! edge along t with face families p, q (E3: p=x1,q=x2; E1: p=x2,q=x3; E2: p=x3,q=x1):
//!   E.t = { -[a_p^W (v_p B_q)^W + a_p^E (v_p B_q)^E] + [a_q^S (v_q B_p)^S + a_q^N (v_q B_p)^N]
//!           + (d_p^E B_q^E - d_p^W B_q^W) - (d_q^N B_p^N - d_q^S B_p^S) } / sin(p,q)
//! W/E = the lower/upper side along p, S/N along q; a, d = the means over the two p- (q-)
//! faces at the edge of the coefficients UctFaceStore saved (a^W = a^L, a^E = 1 - a^L,
//! d^W = d^L, d^E = d^R); B_q = the staggered field of the q faces (B.n_q), PLM along p
//! to the edge; v_p = v.n_p of the q faces (built from their a-weighted face-frame
//! velocity) carried to the edge with half the PLM slope of the cell-centred v.n_p of the
//! two cells sharing the face.  In 1-D this is exactly the face flux.  Writing E.t with
//! the NORMAL components v.n and B.n of the two face families makes it frame-free:
//!   (v x B).t = [(v.n_p)(B.n_q) - (v.n_q)(B.n_p)] / sin(p,q)  for n_p, n_q ⊥ t,
//! and a jump of B.n_q across a p-face Riemann fan is sin(p,q) times the jump of the
//! face-parallel B_t that the a/d coefficients act on, so the d terms carry the same
//! 1/sin.  Cartesian: sin = 1 and v.n = the velocity component.  CUBED SPHERE (gnomonic
//! frames, gnomonic_kernels.hpp; c, s = cos, sin of the xi-eta angle; w0 holds the
//! contravariant (v^xi, v^eta) on the unit tangents, so v.n_xi = s v^xi, v.n_eta = s v^eta):
//!   x1 face {rhat, e_xi, f2 = n_eta}: v.n_eta = slot 5, v.n_xi = s slot 4 - c slot 5;
//!   x2 face {n_xi, e_eta, rhat}: v.rhat = slot 5, v.n_eta = s slot 4 - c slot 3;
//!   x3 face {n_eta, rhat, e_xi}: v.rhat = slot 4, v.n_xi = s slot 5 - c slot 3;
//! E1 = E.rhat has sin(xi,eta) = s at the corner; E2 = E.e_xi and E3 = E.e_eta have
//! sin = 1 (rhat is orthogonal to n_xi, n_eta), and those are the edge projections the
//! face EMFs and the GS07 corner EMF carry (GnomonicEquiangleEmfX1).
struct UctCorner {
  DvceArray5D<Real> u1, u2, u3, w;            // face stores; primitives
  DvceArray4D<Real> b1, b2, b3;               // staggered field
  DvceArray3D<Real> sc, cc, sx, cx, se, ce;   // cs trig: cell, xi face, eta face
  bool cs;
  int js, ks, nx2, nx3;
  DualArray1D<RegionSize> size;
};

KOKKOS_INLINE_FUNCTION
Real UctHalf(const Real qm, const Real q0, const Real qp, const bool upper) {
  Real ql, qr;
  PLM(qm, q0, qp, ql, qr);
  return upper ? (ql - q0) : (qr - q0);
}
KOKKOS_INLINE_FUNCTION
Real UctEdge(const Real aw, const Real dw, const Real de, const Real vpw, const Real vpe,
             const Real bqw, const Real bqe, const Real as, const Real ds, const Real dn,
             const Real vqs, const Real vqn, const Real bps, const Real bpn) {
  return -(aw*vpw*bqw + (1.0 - aw)*vpe*bqe) + (as*vqs*bps + (1.0 - as)*vqn*bpn)
         + (de*bqe - dw*bqw) - (dn*bpn - ds*bps);
}
// cell-centred v.n_xi and v.n_eta (Cartesian: v2, v3)
KOKKOS_INLINE_FUNCTION
Real UctVnXi(const UctCorner &d, const int m, const int k, const int j, const int i) {
  return d.cs ? d.sc(m,k,j)*d.w(m,IVY,k,j,i) : d.w(m,IVY,k,j,i);
}
KOKKOS_INLINE_FUNCTION
Real UctVnEta(const UctCorner &d, const int m, const int k, const int j, const int i) {
  return d.cs ? d.sc(m,k,j)*d.w(m,IVZ,k,j,i) : d.w(m,IVZ,k,j,i);
}

// E3 (x3 edge at x1 face i, x2 face j): p = x1, q = x2
KOKKOS_INLINE_FUNCTION
Real UctE3(const UctCorner &d, const int m, const int k, const int j, const int i) {
  const Real aw = 0.5*(d.u1(m,0,k,j-1,i) + d.u1(m,0,k,j,i));
  const Real dw = 0.5*(d.u1(m,1,k,j-1,i) + d.u1(m,1,k,j,i));
  const Real de = 0.5*(d.u1(m,2,k,j-1,i) + d.u1(m,2,k,j,i));
  const Real as = 0.5*(d.u2(m,0,k,j,i-1) + d.u2(m,0,k,j,i));
  const Real ds = 0.5*(d.u2(m,1,k,j,i-1) + d.u2(m,1,k,j,i));
  const Real dn = 0.5*(d.u2(m,2,k,j,i-1) + d.u2(m,2,k,j,i));
  // v.n_1 (v1) of the x2 faces (slot 5) along x1
  const Real vpw = d.u2(m,5,k,j,i-1) + 0.5*(
    UctHalf(d.w(m,IVX,k,j-1,i-2), d.w(m,IVX,k,j-1,i-1), d.w(m,IVX,k,j-1,i), true) +
    UctHalf(d.w(m,IVX,k,j  ,i-2), d.w(m,IVX,k,j  ,i-1), d.w(m,IVX,k,j  ,i), true));
  const Real vpe = d.u2(m,5,k,j,i) + 0.5*(
    UctHalf(d.w(m,IVX,k,j-1,i-1), d.w(m,IVX,k,j-1,i), d.w(m,IVX,k,j-1,i+1), false) +
    UctHalf(d.w(m,IVX,k,j  ,i-1), d.w(m,IVX,k,j  ,i), d.w(m,IVX,k,j  ,i+1), false));
  // v.n_2 of the x1 faces along x2
  Real fs = d.u1(m,4,k,j-1,i), fn = d.u1(m,4,k,j,i);
  if (d.cs) {
    fs = d.sc(m,k,j-1)*fs - d.cc(m,k,j-1)*d.u1(m,5,k,j-1,i);
    fn = d.sc(m,k,j  )*fn - d.cc(m,k,j  )*d.u1(m,5,k,j  ,i);
  }
  const Real vqs = fs + 0.5*(
    UctHalf(UctVnXi(d,m,k,j-2,i-1), UctVnXi(d,m,k,j-1,i-1), UctVnXi(d,m,k,j,i-1), true) +
    UctHalf(UctVnXi(d,m,k,j-2,i  ), UctVnXi(d,m,k,j-1,i  ), UctVnXi(d,m,k,j,i  ), true));
  const Real vqn = fn + 0.5*(
    UctHalf(UctVnXi(d,m,k,j-1,i-1), UctVnXi(d,m,k,j,i-1), UctVnXi(d,m,k,j+1,i-1), false) +
    UctHalf(UctVnXi(d,m,k,j-1,i  ), UctVnXi(d,m,k,j,i  ), UctVnXi(d,m,k,j+1,i  ), false));
  Real dum, bqw, bqe, bps, bpn;
  PLM(d.b2(m,k,j,i-2), d.b2(m,k,j,i-1), d.b2(m,k,j,i  ), bqw, dum);
  PLM(d.b2(m,k,j,i-1), d.b2(m,k,j,i  ), d.b2(m,k,j,i+1), dum, bqe);
  PLM(d.b1(m,k,j-2,i), d.b1(m,k,j-1,i), d.b1(m,k,j  ,i), bps, dum);
  PLM(d.b1(m,k,j-1,i), d.b1(m,k,j  ,i), d.b1(m,k,j+1,i), dum, bpn);
  return UctEdge(aw, dw, de, vpw, vpe, bqw, bqe, as, ds, dn, vqs, vqn, bps, bpn);
}

// E1 (x1 edge at x2 face j, x3 face k): p = x2, q = x3; / sin(xi,eta) at the corner
KOKKOS_INLINE_FUNCTION
Real UctE1(const UctCorner &d, const int m, const int k, const int j, const int i) {
  const Real aw = 0.5*(d.u2(m,0,k-1,j,i) + d.u2(m,0,k,j,i));
  const Real dw = 0.5*(d.u2(m,1,k-1,j,i) + d.u2(m,1,k,j,i));
  const Real de = 0.5*(d.u2(m,2,k-1,j,i) + d.u2(m,2,k,j,i));
  const Real as = 0.5*(d.u3(m,0,k,j-1,i) + d.u3(m,0,k,j,i));
  const Real ds = 0.5*(d.u3(m,1,k,j-1,i) + d.u3(m,1,k,j,i));
  const Real dn = 0.5*(d.u3(m,2,k,j-1,i) + d.u3(m,2,k,j,i));
  // v.n_2 of the x3 faces along x2
  Real fw = d.u3(m,5,k,j-1,i), fe = d.u3(m,5,k,j,i);
  if (d.cs) {
    fw = d.se(m,k,j-1)*fw - d.ce(m,k,j-1)*d.u3(m,3,k,j-1,i);
    fe = d.se(m,k,j  )*fe - d.ce(m,k,j  )*d.u3(m,3,k,j  ,i);
  }
  const Real vpw = fw + 0.5*(
    UctHalf(UctVnXi(d,m,k-1,j-2,i), UctVnXi(d,m,k-1,j-1,i), UctVnXi(d,m,k-1,j,i), true) +
    UctHalf(UctVnXi(d,m,k  ,j-2,i), UctVnXi(d,m,k  ,j-1,i), UctVnXi(d,m,k  ,j,i), true));
  const Real vpe = fe + 0.5*(
    UctHalf(UctVnXi(d,m,k-1,j-1,i), UctVnXi(d,m,k-1,j,i), UctVnXi(d,m,k-1,j+1,i), false) +
    UctHalf(UctVnXi(d,m,k  ,j-1,i), UctVnXi(d,m,k  ,j,i), UctVnXi(d,m,k  ,j+1,i), false));
  // v.n_3 of the x2 faces along x3
  Real fs = d.u2(m,4,k-1,j,i), fn = d.u2(m,4,k,j,i);
  if (d.cs) {
    fs = d.sx(m,k-1,j)*fs - d.cx(m,k-1,j)*d.u2(m,3,k-1,j,i);
    fn = d.sx(m,k  ,j)*fn - d.cx(m,k  ,j)*d.u2(m,3,k  ,j,i);
  }
  const Real vqs = fs + 0.5*(
    UctHalf(UctVnEta(d,m,k-2,j-1,i), UctVnEta(d,m,k-1,j-1,i), UctVnEta(d,m,k,j-1,i), true) +
    UctHalf(UctVnEta(d,m,k-2,j  ,i), UctVnEta(d,m,k-1,j  ,i), UctVnEta(d,m,k,j  ,i), true));
  const Real vqn = fn + 0.5*(
    UctHalf(UctVnEta(d,m,k-1,j-1,i), UctVnEta(d,m,k,j-1,i), UctVnEta(d,m,k+1,j-1,i), false) +
    UctHalf(UctVnEta(d,m,k-1,j  ,i), UctVnEta(d,m,k,j  ,i), UctVnEta(d,m,k+1,j  ,i), false));
  Real dum, bqw, bqe, bps, bpn;
  PLM(d.b3(m,k,j-2,i), d.b3(m,k,j-1,i), d.b3(m,k,j  ,i), bqw, dum);
  PLM(d.b3(m,k,j-1,i), d.b3(m,k,j  ,i), d.b3(m,k,j+1,i), dum, bqe);
  PLM(d.b2(m,k-2,j,i), d.b2(m,k-1,j,i), d.b2(m,k  ,j,i), bps, dum);
  PLM(d.b2(m,k-1,j,i), d.b2(m,k  ,j,i), d.b2(m,k+1,j,i), dum, bpn);
  const Real e = UctEdge(aw, dw, de, vpw, vpe, bqw, bqe, as, ds, dn, vqs, vqn, bps, bpn);
  if (!d.cs) return e;
  // sin of the xi-eta angle at the corner (xi_{j-1/2}, eta_{k-1/2}), as Coordinates
  // computes the face values (coordinates.cpp, "cscoord")
  const RegionSize &sz = d.size.d_view(m);
  const Real x = tan(0.25*M_PI*LeftEdgeX(j-d.js, d.nx2, sz.x2min, sz.x2max));
  const Real y = tan(0.25*M_PI*LeftEdgeX(k-d.ks, d.nx3, sz.x3min, sz.x3max));
  const Real sn = sqrt(1.0 + SQR(x) + SQR(y))/(sqrt(1.0 + SQR(x))*sqrt(1.0 + SQR(y)));
  return e/sn;
}

// E2 (x2 edge at x3 face k, x1 face i): p = x3, q = x1
KOKKOS_INLINE_FUNCTION
Real UctE2(const UctCorner &d, const int m, const int k, const int j, const int i) {
  const Real aw = 0.5*(d.u3(m,0,k,j,i-1) + d.u3(m,0,k,j,i));
  const Real dw = 0.5*(d.u3(m,1,k,j,i-1) + d.u3(m,1,k,j,i));
  const Real de = 0.5*(d.u3(m,2,k,j,i-1) + d.u3(m,2,k,j,i));
  const Real as = 0.5*(d.u1(m,0,k-1,j,i) + d.u1(m,0,k,j,i));
  const Real ds = 0.5*(d.u1(m,1,k-1,j,i) + d.u1(m,1,k,j,i));
  const Real dn = 0.5*(d.u1(m,2,k-1,j,i) + d.u1(m,2,k,j,i));
  // v.n_3 of the x1 faces (slot 5) along x3
  const Real vpw = d.u1(m,5,k-1,j,i) + 0.5*(
    UctHalf(UctVnEta(d,m,k-2,j,i-1), UctVnEta(d,m,k-1,j,i-1), UctVnEta(d,m,k,j,i-1), true) +
    UctHalf(UctVnEta(d,m,k-2,j,i  ), UctVnEta(d,m,k-1,j,i  ), UctVnEta(d,m,k,j,i  ), true));
  const Real vpe = d.u1(m,5,k,j,i) + 0.5*(
    UctHalf(UctVnEta(d,m,k-1,j,i-1), UctVnEta(d,m,k,j,i-1), UctVnEta(d,m,k+1,j,i-1), false) +
    UctHalf(UctVnEta(d,m,k-1,j,i  ), UctVnEta(d,m,k,j,i  ), UctVnEta(d,m,k+1,j,i  ), false));
  // v.n_1 (v1) of the x3 faces (slot 4) along x1
  const Real vqs = d.u3(m,4,k,j,i-1) + 0.5*(
    UctHalf(d.w(m,IVX,k-1,j,i-2), d.w(m,IVX,k-1,j,i-1), d.w(m,IVX,k-1,j,i), true) +
    UctHalf(d.w(m,IVX,k  ,j,i-2), d.w(m,IVX,k  ,j,i-1), d.w(m,IVX,k  ,j,i), true));
  const Real vqn = d.u3(m,4,k,j,i) + 0.5*(
    UctHalf(d.w(m,IVX,k-1,j,i-1), d.w(m,IVX,k-1,j,i), d.w(m,IVX,k-1,j,i+1), false) +
    UctHalf(d.w(m,IVX,k  ,j,i-1), d.w(m,IVX,k  ,j,i), d.w(m,IVX,k  ,j,i+1), false));
  Real dum, bqw, bqe, bps, bpn;
  PLM(d.b1(m,k-2,j,i), d.b1(m,k-1,j,i), d.b1(m,k  ,j,i), bqw, dum);
  PLM(d.b1(m,k-1,j,i), d.b1(m,k  ,j,i), d.b1(m,k+1,j,i), dum, bqe);
  PLM(d.b3(m,k,j,i-2), d.b3(m,k,j,i-1), d.b3(m,k,j,i  ), bps, dum);
  PLM(d.b3(m,k,j,i-1), d.b3(m,k,j,i  ), d.b3(m,k,j,i+1), dum, bpn);
  return UctEdge(aw, dw, de, vpw, vpe, bqw, bqe, as, ds, dn, vqs, vqn, bps, bpn);
}

//----------------------------------------------------------------------------------------
//! \fn  void MHD::CornerE
//  \brief calculate the corner electric fields.

TaskStatus MHD::CornerE(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  auto &flat = pmy_pack->pcoord->coord_data.is_minkowski;
  auto &spin = pmy_pack->pcoord->coord_data.bh_spin;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  const bool use_bs = use_bs_emf;
  const Real emf_sm_ = emf_smooth;   // lhlld2-1102: < 0 = off (GS07 sign switch)
  // rho*e of a cell for EmfSmoothW (isothermal: rho*c_iso^2); unused when emf_sm_ < 0
  const bool emf_has_e_ = (nmhd == 5);
  const Real emf_cs2_ = (emf_has_e_ ? 0.0 : SQR(peos->eos_data.iso_cs));
  // A REFLECTING radial wall (spherical polar or cubed sphere): the x1 mass flux there
  // is zero (the two wall states are mirror images, mhd_fluxes.cpp), so GS07 has no
  // upwind side.  ">= 0" then always took the GHOST cell's EMF gradient, which behind a
  // mirror wall is minus the interior one: the wall-edge EMF picked up an O(v_r B_t)
  // term, and CT turned it into a growing B.rhat at a closed wall (closed MHD vertex
  // blast, cs_test iprob=12: ~3e-3 of |B|, rising with resolution).  At such a wall the
  // x1 choices take the average of the two sides -- the GS07 rule for a zero contact
  // speed -- which cancels exactly for mirrored states.  The transverse (x2/x3) choices
  // at a wall edge read the mass flux of the INTERIOR row on both sides: the ghost-row
  // flux is the mirror of it everywhere except on a cubed-sphere seam face, where the
  // shared seam flux is applied to the active rows only, and there the two sides picked
  // different upwind cells.  Each piece is needed: with the x1 average, the interior
  // row, or the mirrored ghost bcc1 (CellCenteredRadialFldWall) left out, the blast wall
  // B.rhat/|B| at nx2 = 32 is 1.6e-3, 2.9e-3, 4.1e-4 (6.7e-5 with all three).
  const bool sph_wall_ = use_cubed_sphere || pmy_pack->pmesh->use_spherical_polar;
  const BoundaryFlag rflag_ = BoundaryFlag::reflect;
  // STRESS-FREE closed USER inner wall (dhj wall_mhd_mirror = stress_free; mhd_fluxes):
  // its mirror state has v_t and B_n even, v_n and B_t odd, so every tangential EMF of
  // the mirror row equals the interior one.  The ghost row itself holds the problem's
  // own state (the dhj ghosts are the initial column at rest), so at a wall edge the
  // GS07 terms read the INTERIOR row in place of the ghost row (im_ = i), and the x1
  // choices are those of a zero-mass-flux wall (wall_in, as for reflect).
  const bool wsf_ = wall_closed_ix1 && wall_stress_free;
  const BoundaryFlag uflag_ = BoundaryFlag::user;
  // <mhd>/emf_scheme = uct | uct_hlld: the MD21 corner formula (UctE1/E2/E3) replaces
  // the GS07 weights per corner inside the kernels below; the cell-centred EMFs are not
  // needed then.  At a reflecting / stress-free wall edge it takes the mean of the two
  // wall-face EMFs (the solver's mirror-state flux), reading no ghost row.
  const bool uct_ = emf_uct || emf_uct_hlld;
  UctCorner ucd;
  if (uct_) {
    ucd.u1 = uct1; ucd.u2 = uct2; ucd.u3 = uct3; ucd.w = w0;
    ucd.b1 = b0.x1f; ucd.b2 = b0.x2f; ucd.b3 = b0.x3f;
    ucd.cs = use_cubed_sphere;
    if (use_cubed_sphere) {
      ucd.sc = pmy_pack->pcoord->sin_cell;     ucd.cc = pmy_pack->pcoord->cos_cell;
      ucd.sx = pmy_pack->pcoord->sin_face_xi;  ucd.cx = pmy_pack->pcoord->cos_face_xi;
      ucd.se = pmy_pack->pcoord->sin_face_eta; ucd.ce = pmy_pack->pcoord->cos_face_eta;
    }
    ucd.js = js; ucd.ks = ks; ucd.nx2 = indcs.nx2; ucd.nx3 = indcs.nx3;
    ucd.size = size;
  }

  //---- 1-D problem:
  //  copy face-centered E-fields to edges and return.
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    par_for("emf1", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(int m, int i) {
      e2(m,ks  ,js  ,i) = e2x1_(m,ks,js,i);
      e2(m,ke+1,js  ,i) = e2x1_(m,ks,js,i);
      e3(m,ks  ,js  ,i) = e3x1_(m,ks,js,i);
      e3(m,ks  ,je+1,i) = e3x1_(m,ks,js,i);
    });
  }

  //---- 2-D problem:
  // Copy face-centered E1 and E2 to edges, use GS07 algorithm to compute E3

  if (pmy_pack->pmesh->two_d) {
    // Compute cell-centered E3 = -(v X B) = VyBx-VxBy
    auto w0_ = w0;
    auto bcc_ = bcc0;
    auto e3cc_ = e3_cc;

    // compute cell-centered EMF in dynamical GRMHD
    if (uct_) {
      // emf_scheme = uct | uct_hlld: no cell-centred EMF needed
    } else if (pmy_pack->padm != nullptr) {
      auto &adm = pmy_pack->padm->adm;
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        // Calculate the spatial components of the three-velocity
        const Real &ux = w0_(m,IVX,ks,j,i);
        const Real &uy = w0_(m,IVY,ks,j,i);
        const Real &uz = w0_(m,IVZ,ks,j,i);
        Real iW = 1.0/sqrt(1.0
                    + adm.g_dd(m,0,0,ks,j,i)*ux*ux + 2.0*adm.g_dd(m,0,1,ks,j,i)*ux*uy
                    + 2.0*adm.g_dd(m,0,2,ks,j,i)*ux*uz + adm.g_dd(m,1,1,ks,j,i)*uy*uy
                    + 2.0*adm.g_dd(m,1,2,ks,j,i)*uy*uz + adm.g_dd(m,2,2,ks,j,i)*uz*uz);
        Real v1 = ux*iW;
        Real v2 = uy*iW;
        //Real v3 = uz*iW;

        const Real &alpha = adm.alpha(m,ks,j,i);
        e3cc_(m,ks,j,i) = bcc_(m,IBX,ks,j,i)*(alpha*v2 - adm.beta_u(m, 1, ks, j, i))
                        - bcc_(m,IBY,ks,j,i)*(alpha*v1 - adm.beta_u(m, 0, ks, j, i));
      });
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      // compute cell-centered EMF in GR MHD
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        // Extract components of metric
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(0, indcs.nx3, x3min, x3max);

        Real glower[4][4], gupper[4][4];
        ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

        const Real &ux = w0_(m,IVX,ks,j,i);
        const Real &uy = w0_(m,IVY,ks,j,i);
        const Real &uz = w0_(m,IVZ,ks,j,i);
        const Real &bx = bcc_(m,IBX,ks,j,i);
        const Real &by = bcc_(m,IBY,ks,j,i);
        const Real &bz = bcc_(m,IBZ,ks,j,i);
        // Calculate 4-velocity
        Real tmp = glower[1][1]*ux*ux + 2.0*glower[1][2]*ux*uy + 2.0*glower[1][3]*ux*uz
                 + glower[2][2]*uy*uy + 2.0*glower[2][3]*uy*uz
                 + glower[3][3]*uz*uz;
        Real alpha = sqrt(-1.0/gupper[0][0]);
        Real gamma = sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = ux - alpha * gamma * gupper[0][1];
        Real u2 = uy - alpha * gamma * gupper[0][2];
        Real u3 = uz - alpha * gamma * gupper[0][3];
        // lower vector indices
        Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
        Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
        Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;
        // calculate 4-magnetic field
        Real b0 = u_1*bx + u_2*by + u_3*bz;
        Real b1 = (bx + b0 * u1) / u0;
        Real b2 = (by + b0 * u2) / u0;
        //Real b3 = (bz + b0 * u3) / u0;

        e3cc_(m,ks,j,i) = b1 * u2 - b2 * u1;
      });

    // compute cell-centered EMF in SR MHD
    } else if (pmy_pack->pcoord->is_special_relativistic) {
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        const Real &u1 = w0_(m,IVX,ks,j,i);
        const Real &u2 = w0_(m,IVY,ks,j,i);
        const Real &u3 = w0_(m,IVZ,ks,j,i);
        Real u0 = sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
        e3cc_(m,ks,j,i) = (u2 * bcc_(m,IBX,ks,j,i) - u1 * bcc_(m,IBY,ks,j,i)) / u0;
      });

    // compute cell-centered EMF in Newtonian MHD
    } else {
      par_for("e_cc_2d", DevExeSpace(), 0, nmb1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int j, int i) {
        e3cc_(m,ks,j,i) = w0_(m,IVY,ks,j,i)*bcc_(m,IBX,ks,j,i) -
                          w0_(m,IVX,ks,j,i)*bcc_(m,IBY,ks,j,i);
      });
    }

    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    auto e1x2_ = e1x2;
    auto e3x2_ = e3x2;
    auto flx1 = uflx.x1f;
    auto flx2 = uflx.x2f;

    // integrate E3 to corner using SG07
    //  Note e1[is:ie,  js:je+1,ks:ke+1]
    //       e2[is:ie+1,js:je,  ks:ke+1]
    //       e3[is:ie+1,js:je+1,ks:ke  ]
    par_for("emf2", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      // reflecting radial wall: GS07 x1 upwind choices averaged (see sph_wall_)
      const bool wall_sf = wsf_ && i == is &&
                           mb_bcs.d_view(m,BoundaryFace::inner_x1) == uflag_;
      const bool wall_in = (sph_wall_ && i == is &&
                            mb_bcs.d_view(m,BoundaryFace::inner_x1) == rflag_) || wall_sf;
      const int im_ = wall_sf ? i : i-1;   // the row standing in for the ghost row
      const bool wall_out = sph_wall_ && i == ie+1 &&
                            mb_bcs.d_view(m,BoundaryFace::outer_x1) == rflag_;
      const bool wall1 = wall_in || wall_out;
      // the row whose transverse mass flux picks the upwind side: the interior one on
      // BOTH sides of a wall (the ghost row is its mirror image except where a seam
      // flux was shared over the active rows only)
      const int il_ = wall_in ? i : i-1;
      const int ir_ = wall_out ? i-1 : i;
      e2(m,ks  ,j,i) = e2x1_(m,ks,j,i);
      e2(m,ke+1,j,i) = e2x1_(m,ks,j,i);
      e1(m,ks  ,j,i) = e1x2_(m,ks,j,i);
      e1(m,ke+1,j,i) = e1x2_(m,ks,j,i);
      if (uct_) {   // emf_scheme = uct | uct_hlld (UctE3)
        e3(m,ks,j,i) = wall1 ? 0.5*(e3x1_(m,ks,j-1,i) + e3x1_(m,ks,j,i))
                             : UctE3(ucd, m, ks, j, i);
        return;
      }

      Real e3_l2, e3_r2, e3_l1, e3_r1;
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e3_l2 = 0.5*((e3x2_(m,ks,j,im_) - e3cc_(m,ks,j-1,im_)) +
              (e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j-1,i  )));
      } else if (flx1(m,IDN,ks,j-1,i) >= 0.0) {
        e3_l2 = e3x2_(m,ks,j,im_) - e3cc_(m,ks,j-1,im_);
      } else {
        e3_l2 = e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j-1,i  );
      }
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e3_r2 = 0.5*((e3x2_(m,ks,j,im_) - e3cc_(m,ks,j  ,im_)) +
              (e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j  ,i  )));
      } else if (flx1(m,IDN,ks,j,i) >= 0.0) {
        e3_r2 = e3x2_(m,ks,j,im_) - e3cc_(m,ks,j  ,im_);
      } else {
        e3_r2 = e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j  ,i  );
      }
      if (flx2(m,IDN,ks,j,il_) >= 0.0) {
        e3_l1 = e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,im_);
      } else {
        e3_l1 = e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,im_);
      }
      if (flx2(m,IDN,ks,j,ir_) >= 0.0) {
        e3_r1 = e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,i  );
      } else {
        e3_r1 = e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,i  );
      }
      if (emf_sm_ >= 0.0) {
        // lhlld2-1102 <mhd>/emf_smooth: smooth GS07 upwind weights (see EmfSmoothW)
        auto re = [&](const int jj, const int ii) {
          return w0_(m,IDN,ks,jj,ii)*(emf_has_e_ ? fabs(w0_(m,IEN,ks,jj,ii)) : emf_cs2_);
        };
        if (!wall1) {
          e3_l2 = EmfPick(flx1(m,IDN,ks,j-1,i), emf_sm_, re(j-1,im_), re(j-1,i),
                          e3x2_(m,ks,j,im_) - e3cc_(m,ks,j-1,im_),
                          e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j-1,i  ));
          e3_r2 = EmfPick(flx1(m,IDN,ks,j,i), emf_sm_, re(j,im_), re(j,i),
                          e3x2_(m,ks,j,im_) - e3cc_(m,ks,j  ,im_),
                          e3x2_(m,ks,j,i  ) - e3cc_(m,ks,j  ,i  ));
        }
        e3_l1 = EmfPick(flx2(m,IDN,ks,j,il_), emf_sm_, re(j-1,im_), re(j,im_),
                        e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,im_),
                        e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,im_));
        e3_r1 = EmfPick(flx2(m,IDN,ks,j,ir_), emf_sm_, re(j-1,i), re(j,i),
                        e3x1_(m,ks,j-1,i) - e3cc_(m,ks,j-1,i  ),
                        e3x1_(m,ks,j  ,i) - e3cc_(m,ks,j  ,i  ));
      }
      e3(m,ks,j,i) = 0.25*(e3_l1 + e3_r1 + e3_l2 + e3_r2 +
             e3x2_(m,ks,j,im_) + e3x2_(m,ks,j,i) + e3x1_(m,ks,j-1,i) + e3x1_(m,ks,j,i));
    });
  }

  //---- 3-D problem:
  // Use GS07 algorithm to compute all three of E1, E2, and E3

  if (pmy_pack->pmesh->three_d) {
    // Compute cell-centered electric fields
    // E1=-(v X B)=VzBy-VyBz
    // E2=-(v X B)=VxBz-VzBx
    // E3=-(v X B)=VyBx-VxBy
    auto w0_ = w0;
    auto bcc_ = bcc0;
    auto e1cc_ = e1_cc;
    auto e2cc_ = e2_cc;
    auto e3cc_ = e3_cc;

    // compute cell-centered EMFs in dynamical GRMHD
    if (uct_) {
      // emf_scheme = uct | uct_hlld: no cell-centred EMF needed
    } else if (pmy_pack->padm != nullptr) {
      auto &adm = pmy_pack->padm->adm;
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        // Calculate something that resembles the spatial components of the four-velocity
        // normalized by W.
        const Real &ux = w0_(m,IVX,k,j,i);
        const Real &uy = w0_(m,IVY,k,j,i);
        const Real &uz = w0_(m,IVZ,k,j,i);
        const Real &bx = bcc_(m,IBX,k,j,i);
        const Real &by = bcc_(m,IBY,k,j,i);
        const Real &bz = bcc_(m,IBZ,k,j,i);
        Real iW = 1.0/sqrt(1.0
                         + adm.g_dd(m,0,0,k,j,i)*ux*ux + 2.0*adm.g_dd(m,0,1,k,j,i)*ux*uy
                         + 2.0*adm.g_dd(m,0,2,k,j,i)*ux*uz + adm.g_dd(m,1,1,k,j,i)*uy*uy
                         + 2.0*adm.g_dd(m,1,2,k,j,i)*uy*uz + adm.g_dd(m,2,2,k,j,i)*uz*uz);
        const Real &alpha = adm.alpha(m, k, j, i);
        Real v1c = alpha*ux*iW - adm.beta_u(m, 0, k, j, i);
        Real v2c = alpha*uy*iW - adm.beta_u(m, 1, k, j, i);
        Real v3c = alpha*uz*iW - adm.beta_u(m, 2, k, j, i);

        e1cc_(m,k,j,i) = by * v3c - bz * v2c;
        e2cc_(m,k,j,i) = bz * v1c - bx * v3c;
        e3cc_(m,k,j,i) = bx * v2c - by * v1c;
      });
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      // compute cell-centered EMFs in GR MHD
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        // Extract components of metric
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);

        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

        Real glower[4][4], gupper[4][4];
        ComputeMetricAndInverse(x1v, x2v, x3v, flat, spin, glower, gupper);

        const Real &ux = w0_(m,IVX,k,j,i);
        const Real &uy = w0_(m,IVY,k,j,i);
        const Real &uz = w0_(m,IVZ,k,j,i);
        const Real &bx = bcc_(m,IBX,k,j,i);
        const Real &by = bcc_(m,IBY,k,j,i);
        const Real &bz = bcc_(m,IBZ,k,j,i);
        // Calculate 4-velocity
        Real tmp = glower[1][1]*ux*ux + 2.0*glower[1][2]*ux*uy + 2.0*glower[1][3]*ux*uz
                 + glower[2][2]*uy*uy + 2.0*glower[2][3]*uy*uz
                 + glower[3][3]*uz*uz;
        Real alpha = sqrt(-1.0/gupper[0][0]);
        Real gamma = sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = ux - alpha * gamma * gupper[0][1];
        Real u2 = uy - alpha * gamma * gupper[0][2];
        Real u3 = uz - alpha * gamma * gupper[0][3];
        // lower vector indices
        Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
        Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
        Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;
        // calculate 4-magnetic field
        Real b0 = u_1*bx + u_2*by + u_3*bz;
        Real b1 = (bx + b0 * u1) / u0;
        Real b2 = (by + b0 * u2) / u0;
        Real b3 = (bz + b0 * u3) / u0;

        e1cc_(m,k,j,i) = b2 * u3 - b3 * u2;
        e2cc_(m,k,j,i) = b3 * u1 - b1 * u3;
        e3cc_(m,k,j,i) = b1 * u2 - b2 * u1;
      });

    // compute cell-centered EMFs in SR MHD
    } else if (pmy_pack->pcoord->is_special_relativistic) {
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const Real &u1 = w0_(m,IVX,k,j,i);
        const Real &u2 = w0_(m,IVY,k,j,i);
        const Real &u3 = w0_(m,IVZ,k,j,i);
        Real u0 = sqrt(1.0 + SQR(u1) + SQR(u2) + SQR(u3));
        e1cc_(m,k,j,i) = (u3 * bcc_(m,IBY,k,j,i) - u2 * bcc_(m,IBZ,k,j,i)) / u0;
        e2cc_(m,k,j,i) = (u1 * bcc_(m,IBZ,k,j,i) - u3 * bcc_(m,IBX,k,j,i)) / u0;
        e3cc_(m,k,j,i) = (u2 * bcc_(m,IBX,k,j,i) - u1 * bcc_(m,IBY,k,j,i)) / u0;
      });

    // CUBED SPHERE, the gnomonic form of the same thing.  v x B is only a cross product
    // in an ORTHONORMAL frame, and w0's velocity is CONTRAVARIANT on the non-orthogonal
    // unit tangent basis while bcc is already in the orthonormal triple
    // (B.rhat, B.e_xi, B.f2), f2 = (e_eta - c e_xi)/s.  So put the velocity in that same
    // triple first -- u = (v1, v2 + c v3, s v3), exactly what GnomonicEquianglePrimFaceX1
    // does -- take the cross product there, and project the result onto the EDGE
    // directions, which is what the face EMFs carry and what the corner EMF must be: the
    // r and xi axes ARE edges, and E.e_eta = c*E.e_xi + s*E.f2, the same rotation
    // GnomonicEquiangleEmfX1 applies to the x1-sweep EMF.
    //
    // (rhat, e_xi, f2) is right-handed, so the cross product needs no extra sign:
    // e_xi x f2 = e_xi x (e_eta - c e_xi)/s = (e_xi x e_eta)/s = rhat.
    } else if (use_cubed_sphere) {
      auto sin_ = pmy_pack->pcoord->sin_cell;
      auto cos_ = pmy_pack->pcoord->cos_cell;
      par_for("e_cc_3d_cs", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        const Real sn = sin_(m,k,j);
        const Real cs = cos_(m,k,j);
        const Real u1 = w0_(m,IVX,k,j,i);
        const Real u2 = w0_(m,IVY,k,j,i) + cs*w0_(m,IVZ,k,j,i);
        const Real u3 = sn*w0_(m,IVZ,k,j,i);
        const Real b1 = bcc_(m,IBX,k,j,i);
        const Real b2 = bcc_(m,IBY,k,j,i);
        const Real b3 = bcc_(m,IBZ,k,j,i);
        const Real er  = u3*b2 - u2*b3;
        const Real exi = u1*b3 - u3*b1;
        const Real ef2 = u2*b1 - u1*b2;
        e1cc_(m,k,j,i) = er;
        e2cc_(m,k,j,i) = exi;
        e3cc_(m,k,j,i) = cs*exi + sn*ef2;
      });

    // compute cell-centered EMFs in Newtonian MHD
    } else {
      par_for("e_cc_3d", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        e1cc_(m,k,j,i) = w0_(m,IVZ,k,j,i)*bcc_(m,IBY,k,j,i) -
                         w0_(m,IVY,k,j,i)*bcc_(m,IBZ,k,j,i);
        e2cc_(m,k,j,i) = w0_(m,IVX,k,j,i)*bcc_(m,IBZ,k,j,i) -
                         w0_(m,IVZ,k,j,i)*bcc_(m,IBX,k,j,i);
        e3cc_(m,k,j,i) = w0_(m,IVY,k,j,i)*bcc_(m,IBX,k,j,i) -
                         w0_(m,IVX,k,j,i)*bcc_(m,IBY,k,j,i);
      });
    }

    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2x1_ = e2x1;
    auto e3x1_ = e3x1;
    auto e1x2_ = e1x2;
    auto e3x2_ = e3x2;
    auto e1x3_ = e1x3;
    auto e2x3_ = e2x3;
    auto flx1 = uflx.x1f;
    auto flx2 = uflx.x2f;
    auto flx3 = uflx.x3f;
    // polar-row face-field dissipation (see <mhd>/polar_emf_diss in mhd.cpp)
    const bool polar_diss = use_polar_emf_diss && pmy_pack->pmesh->use_polar_boundary;
    auto b1_ = b0.x1f;
    auto b3_ = b0.x3f;
    auto &eos = peos->eos_data;
    auto wder_ = wder;

    // Integrate E1, E2, E3 to corners
    //  Note e1[is:ie,  js:je+1,ks:ke+1]
    //       e2[is:ie+1,js:je,  ks:ke+1]
    //       e3[is:ie+1,js:je+1,ks:ke  ]
    // CUBED SPHERE: the GS07 upwind correction below needs a CELL-CENTRED reference EMF,
    // and the flat-space v x B is not valid there (contravariant velocity, orthonormal
    // bcc); e_cc_3d_cs above supplies the gnomonic form, so the cubed sphere takes the
    // same GS07 branch as every other grid (a user decision, 2026-09-05; on the cs_test
    // gate: same order, 1.8x the smooth-field L1(B) of the plain average).
    // The plain four-face average (Balsara-Spicer) survives only as the <mhd>/bs_emf
    // diagnostic: it does not reduce to the upwind flux in 1D and leaves the odd-even
    // field mode undamped (it collapsed dt 5x on the spherical-polar dhj run).
    if (use_bs) {
      par_for("emf3_cs", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        e1(m,k,j,i) = 0.25*(e1x2_(m,k-1,j,i) + e1x2_(m,k,j,i) +
                            e1x3_(m,k,j-1,i) + e1x3_(m,k,j,i));
        e2(m,k,j,i) = 0.25*(e2x3_(m,k,j,i-1) + e2x3_(m,k,j,i) +
                            e2x1_(m,k-1,j,i) + e2x1_(m,k,j,i));
        e3(m,k,j,i) = 0.25*(e3x2_(m,k,j,i-1) + e3x2_(m,k,j,i) +
                            e3x1_(m,k,j-1,i) + e3x1_(m,k,j,i));
      });
      if (pmy_pack->pmesh->use_polar_boundary) PolarAzimuthalAverageEr();
      return TaskStatus::complete;
    }

    par_for("emf3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // lhlld2-1102 emf_smooth: rho*e of a cell (EmfSmoothW); unused when emf_sm_ < 0
      auto re3 = [&](const int kk, const int jj, const int ii) {
        return w0_(m,IDN,kk,jj,ii)*(emf_has_e_ ? fabs(w0_(m,IEN,kk,jj,ii)) : emf_cs2_);
      };
      // reflecting radial wall: GS07 x1 upwind choices averaged (see sph_wall_)
      const bool wall_sf = wsf_ && i == is &&
                           mb_bcs.d_view(m,BoundaryFace::inner_x1) == uflag_;
      const bool wall_in = (sph_wall_ && i == is &&
                            mb_bcs.d_view(m,BoundaryFace::inner_x1) == rflag_) || wall_sf;
      const int im_ = wall_sf ? i : i-1;   // the row standing in for the ghost row
      const bool wall_out = sph_wall_ && i == ie+1 &&
                            mb_bcs.d_view(m,BoundaryFace::outer_x1) == rflag_;
      const bool wall1 = wall_in || wall_out;
      // the row whose transverse mass flux picks the upwind side: the interior one on
      // BOTH sides of a wall (the ghost row is its mirror image except where a seam
      // flux was shared over the active rows only)
      const int il_ = wall_in ? i : i-1;
      const int ir_ = wall_out ? i-1 : i;
      // In a polar cell row, add Rusanov dissipation built from the FACE fields on either
      // side of this edge: -a_r/2 * [B3(i) - B3(i-1)] and +a_phi/2 * [B1(k) - B1(k-1)],
      // signs fixed by the CT update (bx3f -= d(e2)/dx1, bx1f += d(e2)/dx3). The GS05
      // terms above see the (-1)^(i+k) face checkerboard only through the phi-face EMF
      // and bias the corner value DOWNWIND for it; this term restores an upwind bias.
      auto polar_diss_e2 = [&](Real &e2v) {
      if (polar_diss) {
        const bool do_pole =
            (j == js && mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar) ||
            (j == je && mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
        if (do_pole) {
          Real ar = 0.0, ap = 0.0;
          for (int dk=-1; dk<=0; ++dk) {
            for (int di=-1; di<=0; ++di) {
              // radial: ACTIVE cells only (the physical radial ghosts); phi: the ghost
              // cell itself -- a copy of the neighbour's active cell, w0 / bcc / wder
              // included (ConsToPrim runs over all cells after the halo exchange).
              // layoutaudit-1103: the phi index used to be clamped to THIS MeshBlock's
              // ks..ke, so at an internal phi block edge the fast speed came from one
              // side only: the polar-row e2 (and B) depended on the phi decomposition
              // (sp_test resistive/MHD decks, 2 vs 4 phi blocks: B 6e-8 after 1 cycle).
              int ic = (i+di < is) ? is : ((i+di > ie) ? ie : i+di);
              int kc = k+dk;
              Real d = w0_(m,IDN,kc,j,ic);
              Real bx = bcc_(m,IBX,kc,j,ic);
              Real by = bcc_(m,IBY,kc,j,ic);
              Real bz = bcc_(m,IBZ,kc,j,ic);
              Real cf;
              if (eos.IsGeneral()) {
                cf = eos.FastSpeedFromP(d, wder_(m,IDPR,kc,j,ic), wder_(m,IDG1,kc,j,ic),
                                        bx, by, bz);
              } else if (eos.is_ideal) {
                Real p = eos.IdealGasPressure(w0_(m,IEN,kc,j,ic));
                cf = eos.IdealMHDFastSpeed(d, p, bx, by, bz);
              } else {
                cf = eos.IdealMHDFastSpeed(d, bx, by, bz);
              }
              ar = fmax(ar, fabs(w0_(m,IVX,kc,j,ic)) + cf);
              ap = fmax(ap, fabs(w0_(m,IVZ,kc,j,ic)) + cf);
            }
          }
          // THIRD differences, not first.  The first-difference (Rusanov) form damps a
          // mode of wavenumber kappa at the rate a*kappa^2*dx/2, and in the polar row
          // dx_phi = r sin(theta) dphi ~ dtheta*dphi while the smooth m = 1 pattern of
          // B_r has kappa ~ 1/(r sin theta): the rate is a*dphi/(2 r sin theta) = O(1),
          // INDEPENDENT of resolution.  That distorted every smooth field threading the
          // pole (sp_test iprob=11: polar-row v_r error 2.5e-2 at nx2 = 16 AND 32, and
          // 100x smaller with the term off).  -D3/4 equals D1 exactly on the (-1)^k
          // checkerboard this term exists to damp and is O(kappa^3 dx^3) on smooth modes.
          e2v -= 0.5*ar*0.25*(-b3_(m,k,j,i+1) + 3.0*b3_(m,k,j,i)
                                      - 3.0*b3_(m,k,j,i-1) + b3_(m,k,j,i-2));
          e2v += 0.5*ap*0.25*(-b1_(m,k+1,j,i) + 3.0*b1_(m,k,j,i)
                                      - 3.0*b1_(m,k-1,j,i) + b1_(m,k-2,j,i));
        }
      }
      };
      if (uct_) {   // emf_scheme = uct | uct_hlld (UctE1/E2/E3)
        e1(m,k,j,i) = UctE1(ucd, m, k, j, i);
        e2(m,k,j,i) = wall1 ? 0.5*(e2x1_(m,k-1,j,i) + e2x1_(m,k,j,i))
                            : UctE2(ucd, m, k, j, i);
        e3(m,k,j,i) = wall1 ? 0.5*(e3x1_(m,k,j-1,i) + e3x1_(m,k,j,i))
                            : UctE3(ucd, m, k, j, i);
        // spdef-1103: the polar-row face-field dissipation acts on the UCT corner e2
        // too (UCT's d-terms do not see the (-1)^(i+k) mode that closes through the
        // polar ghost row; same term, same <mhd>/polar_emf_diss key as GS07)
        polar_diss_e2(e2(m,k,j,i));
        return;
      }
      // integrate E1 to corner using SG07
      Real e1_l3, e1_r3, e1_l2, e1_r2;
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e1_l3 = EmfPick(flx2(m,IDN,k-1,j,i), emf_sm_, re3(k-1,j-1,i), re3(k-1,j  ,i),
                        e1x3_(m,k,j-1,i) - e1cc_(m,k-1,j-1,i),
                        e1x3_(m,k,j  ,i) - e1cc_(m,k-1,j  ,i));
      } else if (flx2(m,IDN,k-1,j,i) >= 0.0) {
        e1_l3 = e1x3_(m,k,j-1,i) - e1cc_(m,k-1,j-1,i);
      } else {
        e1_l3 = e1x3_(m,k,j  ,i) - e1cc_(m,k-1,j  ,i);
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e1_r3 = EmfPick(flx2(m,IDN,k,j,i), emf_sm_, re3(k  ,j-1,i), re3(k  ,j  ,i),
                        e1x3_(m,k,j-1,i) - e1cc_(m,k  ,j-1,i),
                        e1x3_(m,k,j  ,i) - e1cc_(m,k  ,j  ,i));
      } else if (flx2(m,IDN,k,j,i) >= 0.0) {
        e1_r3 = e1x3_(m,k,j-1,i) - e1cc_(m,k  ,j-1,i);
      } else {
        e1_r3 = e1x3_(m,k,j  ,i) - e1cc_(m,k  ,j  ,i);
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e1_l2 = EmfPick(flx3(m,IDN,k,j-1,i), emf_sm_, re3(k-1,j-1,i), re3(k  ,j-1,i),
                        e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j-1,i),
                        e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j-1,i));
      } else if (flx3(m,IDN,k,j-1,i) >= 0.0) {
        e1_l2 = e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j-1,i);
      } else {
        e1_l2 = e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j-1,i);
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e1_r2 = EmfPick(flx3(m,IDN,k,j,i), emf_sm_, re3(k-1,j  ,i), re3(k  ,j  ,i),
                        e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j  ,i),
                        e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j  ,i));
      } else if (flx3(m,IDN,k,j,i) >= 0.0) {
        e1_r2 = e1x2_(m,k-1,j,i) - e1cc_(m,k-1,j  ,i);
      } else {
        e1_r2 = e1x2_(m,k  ,j,i) - e1cc_(m,k  ,j  ,i);
      }
      e1(m,k,j,i) = 0.25*(e1_l3 + e1_r3 + e1_l2 + e1_r2 +
                e1x2_(m,k-1,j,i) + e1x2_(m,k,j,i) + e1x3_(m,k,j-1,i) + e1x3_(m,k,j,i));

      // integrate E2 to corner using SG07
      Real e2_l3, e2_r3, e2_l1, e2_r1;
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e2_l3 = 0.5*((e2x3_(m,k,j,im_) - e2cc_(m,k-1,j,im_)) +
              (e2x3_(m,k,j,i  ) - e2cc_(m,k-1,j,i  )));
      } else if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e2_l3 = EmfPick(flx1(m,IDN,k-1,j,i), emf_sm_, re3(k-1,j,im_), re3(k-1,j,i  ),
                        e2x3_(m,k,j,im_) - e2cc_(m,k-1,j,im_),
                        e2x3_(m,k,j,i  ) - e2cc_(m,k-1,j,i  ));
      } else if (flx1(m,IDN,k-1,j,i) >= 0.0) {
        e2_l3 = e2x3_(m,k,j,im_) - e2cc_(m,k-1,j,im_);
      } else {
        e2_l3 = e2x3_(m,k,j,i  ) - e2cc_(m,k-1,j,i  );
      }
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e2_r3 = 0.5*((e2x3_(m,k,j,im_) - e2cc_(m,k  ,j,im_)) +
              (e2x3_(m,k,j,i  ) - e2cc_(m,k  ,j,i  )));
      } else if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e2_r3 = EmfPick(flx1(m,IDN,k,j,i), emf_sm_, re3(k  ,j,im_), re3(k  ,j,i  ),
                        e2x3_(m,k,j,im_) - e2cc_(m,k  ,j,im_),
                        e2x3_(m,k,j,i  ) - e2cc_(m,k  ,j,i  ));
      } else if (flx1(m,IDN,k,j,i) >= 0.0) {
        e2_r3 = e2x3_(m,k,j,im_) - e2cc_(m,k  ,j,im_);
      } else {
        e2_r3 = e2x3_(m,k,j,i  ) - e2cc_(m,k  ,j,i  );
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e2_l1 = EmfPick(flx3(m,IDN,k,j,il_), emf_sm_, re3(k-1,j,im_), re3(k  ,j,im_),
                        e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,im_),
                        e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,im_));
      } else if (flx3(m,IDN,k,j,il_) >= 0.0) {
        e2_l1 = e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,im_);
      } else {
        e2_l1 = e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,im_);
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e2_r1 = EmfPick(flx3(m,IDN,k,j,ir_), emf_sm_, re3(k-1,j,i  ), re3(k  ,j,i  ),
                        e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,i  ),
                        e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,i  ));
      } else if (flx3(m,IDN,k,j,ir_) >= 0.0) {
        e2_r1 = e2x1_(m,k-1,j,i) - e2cc_(m,k-1,j,i  );
      } else {
        e2_r1 = e2x1_(m,k  ,j,i) - e2cc_(m,k  ,j,i  );
      }
      e2(m,k,j,i) = 0.25*(e2_l3 + e2_r3 + e2_l1 + e2_r1 +
                e2x3_(m,k,j,im_) + e2x3_(m,k,j,i) + e2x1_(m,k-1,j,i) + e2x1_(m,k,j,i));

      polar_diss_e2(e2(m,k,j,i));   // see polar_diss_e2 above (GS07 path)

      // integrate E3 to corner using SG07
      Real e3_l2, e3_r2, e3_l1, e3_r1;
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e3_l2 = 0.5*((e3x2_(m,k,j,im_) - e3cc_(m,k,j-1,im_)) +
              (e3x2_(m,k,j,i  ) - e3cc_(m,k,j-1,i  )));
      } else if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e3_l2 = EmfPick(flx1(m,IDN,k,j-1,i), emf_sm_, re3(k,j-1,im_), re3(k,j-1,i  ),
                        e3x2_(m,k,j,im_) - e3cc_(m,k,j-1,im_),
                        e3x2_(m,k,j,i  ) - e3cc_(m,k,j-1,i  ));
      } else if (flx1(m,IDN,k,j-1,i) >= 0.0) {
        e3_l2 = e3x2_(m,k,j,im_) - e3cc_(m,k,j-1,im_);
      } else {
        e3_l2 = e3x2_(m,k,j,i  ) - e3cc_(m,k,j-1,i  );
      }
      if (wall1) {  // reflecting radial wall: no upwind side, see above
        e3_r2 = 0.5*((e3x2_(m,k,j,im_) - e3cc_(m,k,j  ,im_)) +
              (e3x2_(m,k,j,i  ) - e3cc_(m,k,j  ,i  )));
      } else if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e3_r2 = EmfPick(flx1(m,IDN,k,j,i), emf_sm_, re3(k,j  ,im_), re3(k,j  ,i  ),
                        e3x2_(m,k,j,im_) - e3cc_(m,k,j  ,im_),
                        e3x2_(m,k,j,i  ) - e3cc_(m,k,j  ,i  ));
      } else if (flx1(m,IDN,k,j,i) >= 0.0) {
        e3_r2 = e3x2_(m,k,j,im_) - e3cc_(m,k,j  ,im_);
      } else {
        e3_r2 = e3x2_(m,k,j,i  ) - e3cc_(m,k,j  ,i  );
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e3_l1 = EmfPick(flx2(m,IDN,k,j,il_), emf_sm_, re3(k,j-1,im_), re3(k,j  ,im_),
                        e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,im_),
                        e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,im_));
      } else if (flx2(m,IDN,k,j,il_) >= 0.0) {
        e3_l1 = e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,im_);
      } else {
        e3_l1 = e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,im_);
      }
      if (emf_sm_ >= 0.0) {  // lhlld2-1102 emf_smooth
        e3_r1 = EmfPick(flx2(m,IDN,k,j,ir_), emf_sm_, re3(k,j-1,i  ), re3(k,j  ,i  ),
                        e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,i  ),
                        e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,i  ));
      } else if (flx2(m,IDN,k,j,ir_) >= 0.0) {
        e3_r1 = e3x1_(m,k,j-1,i) - e3cc_(m,k,j-1,i  );
      } else {
        e3_r1 = e3x1_(m,k,j  ,i) - e3cc_(m,k,j  ,i  );
      }
      e3(m,k,j,i) = 0.25*(e3_l1 + e3_r1 + e3_l2 + e3_r2 +
                e3x2_(m,k,j,im_) + e3x2_(m,k,j,i) + e3x1_(m,k,j-1,i) + e3x1_(m,k,j,i));
    });
  }
  // ucthshare-1102 <mhd>/cs_seam_mhd = shared: the seam-edge / cube-vertex UCT EMFs from
  // shared inputs (mhd_cs_seam.cpp)
  if (cs_seam_mhd >= 2 && pmy_pack->pmesh->three_d) {UctSeamEdges();}
  if (pmy_pack->pmesh->use_polar_boundary) PolarAzimuthalAverageEr();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::WallEmfZero
//! \brief condwall-1028: perfectly conducting wall at rest.  On every REFLECT mesh face
//! of a MeshBlock (refl) and/or on the is face of a USER inner x1 boundary (usr_ix1, the
//! dhj closed wall under wall_mhd_mirror = conducting) the edge EMFs of `e` tangential to
//! that face are set to exactly 0, so CT leaves the face-normal field of the wall faces
//! unchanged.  Generic: Cartesian, spherical polar, cubed sphere (x1 = r there; the panel
//! faces are never reflect).  Edge extents: e1(k: ks..ke+1, j: js..je+1, i: is..ie),
//! e2(k: ks..ke+1, j: js..je, i: is..ie+1), e3(k: ks..ke, j: js..je+1, i: is..ie+1); in
//! 1-D/2-D the degenerate ke+1 / je+1 copies are zeroed as well.  Called on the ideal
//! corner EMF (wall_conducting = full, or the dhj conducting wall) and on the resistive
//! EMF (Resistivity::AddFluxGeneralResist): the explicit and the RKG paths both see it.

void MHD::WallEmfZero(DvceEdgeFld4D<Real> &e, const bool refl, const bool usr_ix1) {
  auto e1 = e.x1e;
  auto e2 = e.x2e;
  auto e3 = e.x3e;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const BoundaryFlag rf = BoundaryFlag::reflect, uf = BoundaryFlag::user;
  // x1 faces: e2 and e3 at i = is / ie+1
  par_for("wall_emf0_x1", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const BoundaryFlag bi = mb_bcs.d_view(m,BoundaryFace::inner_x1);
    const BoundaryFlag bo = mb_bcs.d_view(m,BoundaryFace::outer_x1);
    if ((refl && bi == rf) || (usr_ix1 && bi == uf)) {
      if (j <= je) {e2(m,k,j,is) = 0.0;}
      if (k <= ke) {e3(m,k,j,is) = 0.0;}
    }
    if (refl && bo == rf) {
      if (j <= je) {e2(m,k,j,ie+1) = 0.0;}
      if (k <= ke) {e3(m,k,j,ie+1) = 0.0;}
    }
  });
  if (!refl || !md) return;
  // x2 faces: e1 and e3 at j = js / je+1
  par_for("wall_emf0_x2", DevExeSpace(), 0, nmb1, ks, ke+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == rf) {
      if (i <= ie) {e1(m,k,js,i) = 0.0;}
      if (k <= ke) {e3(m,k,js,i) = 0.0;}
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == rf) {
      if (i <= ie) {e1(m,k,je+1,i) = 0.0;}
      if (k <= ke) {e3(m,k,je+1,i) = 0.0;}
    }
  });
  if (!td) return;
  // x3 faces: e1 and e2 at k = ks / ke+1
  par_for("wall_emf0_x3", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int j, const int i) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x3) == rf) {
      if (i <= ie) {e1(m,ks,j,i) = 0.0;}
      if (j <= je) {e2(m,ks,j,i) = 0.0;}
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x3) == rf) {
      if (i <= ie) {e1(m,ke+1,j,i) = 0.0;}
      if (j <= je) {e2(m,ke+1,j,i) = 0.0;}
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::WallFluxZero
//! \brief condwall-1028: the Riemann MASS, ENERGY and passive-scalar fluxes through a
//! conducting wall are set to 0: every REFLECT mesh face (refl, <mhd>/wall_conducting =
//! full) and/or the is face of a USER inner x1 boundary (usr_ix1, the dhj conducting
//! wall).  At a closed wall at rest v_n = 0, and E_t = 0 (WallEmfZero) removes the
//! Poynting flux E x B . n, so the exact mass and energy fluxes are 0.  The Riemann flux
//! is not: the mirror state (v_n odd, B_t even) is a symmetry of MHD only when B_n = 0.
//! With B_n != 0 the exact mirror would need B_t ODD (a stress-free wall), so the solver
//! returned S_M != 0, a mass flux O(B_n B_t) (cw_test threaded box: 2.6e-2 of M through
//! the walls in t = 1 on Cartesian, 2e-3 on sp) and the work -B_n (v_t . B_t).  Called on
//! the Riemann (+FOFC) fluxes only, before any diffusive flux (conduction / viscosity
//! through the wall are kept).  The momentum fluxes (pressure and Maxwell stress) stay.

void MHD::WallFluxZero(const bool refl, const bool usr_ix1) {
  auto f1 = uflx.x1f;
  auto f2 = uflx.x2f;
  auto f3 = uflx.x3f;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  // the zeroed slots: IDN, IEN (ideal EOS only) and the passive scalars
  const int ien = (nmhd > IEN) ? IEN : -1;
  const int ns0 = nmhd, ns1 = nmhd + nscalars;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  const BoundaryFlag rf = BoundaryFlag::reflect, uf = BoundaryFlag::user;
  par_for("wall_flux0_x1", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const BoundaryFlag bi = mb_bcs.d_view(m,BoundaryFace::inner_x1);
    const bool zi = (refl && bi == rf) || (usr_ix1 && bi == uf);
    const bool zo = refl && mb_bcs.d_view(m,BoundaryFace::outer_x1) == rf;
    for (int n=0; n<ns1; ++n) {
      if (n != IDN && n != ien && n < ns0) continue;
      if (zi) {f1(m,n,k,j,is) = 0.0;}
      if (zo) {f1(m,n,k,j,ie+1) = 0.0;}
    }
  });
  if (!refl || !md) return;
  par_for("wall_flux0_x2", DevExeSpace(), 0, nmb1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    const bool zi = mb_bcs.d_view(m,BoundaryFace::inner_x2) == rf;
    const bool zo = mb_bcs.d_view(m,BoundaryFace::outer_x2) == rf;
    for (int n=0; n<ns1; ++n) {
      if (n != IDN && n != ien && n < ns0) continue;
      if (zi) {f2(m,n,k,js,i) = 0.0;}
      if (zo) {f2(m,n,k,je+1,i) = 0.0;}
    }
  });
  if (!td) return;
  par_for("wall_flux0_x3", DevExeSpace(), 0, nmb1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int j, const int i) {
    const bool zi = mb_bcs.d_view(m,BoundaryFace::inner_x3) == rf;
    const bool zo = mb_bcs.d_view(m,BoundaryFace::outer_x3) == rf;
    for (int n=0; n<ns1; ++n) {
      if (n != IDN && n != ien && n < ns0) continue;
      if (zi) {f3(m,n,ks,j,i) = 0.0;}
      if (zo) {f3(m,n,ke+1,j,i) = 0.0;}
    }
  });
}

void MHD::PolarAzimuthalAverageEr() {
  auto e1 = efld.x1e;

  auto &indcs = pmy_pack->pmesh->mb_indcs;

  int nx1 = indcs.nx1;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  int &ng = indcs.ng;
  const int ncells1 = nx1 + 2*ng;

  const int mesh_nx3 = pmy_pack->pmesh->mesh_indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;

  auto &mb_bcs = pmy_pack->pmb->mb_bcs;

  // Work arrays are members, sized once.  The previous version allocated four Views, a
  // pair of host mirrors and two std::vectors on EVERY call -- twice per cycle.
  // The host mirrors must be in the SAME guard condition as the device arrays and are
  // tested separately, because the MHD constructor pre-sizes inner_local/outer_local
  // (mhd.cpp) but not these two.  A guard on inner_local alone therefore never fires and
  // leaves the mirrors at their default extent 0.  Nothing notices in serial -- only the
  // MPI branch below touches them -- so every multi-rank run with polar boundaries died
  // on the first cycle with "deep_copy extents of views don't match: (0) inner(n)".
  if (static_cast<int>(inner_local.extent(0)) != ncells1 ||
      static_cast<int>(polar_inner_h.extent(0)) != ncells1) {
    Kokkos::realloc(inner_local, ncells1);
    Kokkos::realloc(outer_local, ncells1);
    Kokkos::realloc(polar_inner_h, ncells1);
    Kokkos::realloc(polar_outer_h, ncells1);
  }
  if (static_cast<int>(polar_part_in.extent(0)) != nmb ||
      static_cast<int>(polar_part_in.extent(1)) != ncells1) {
    Kokkos::realloc(polar_part_in,  nmb, ncells1);
    Kokkos::realloc(polar_part_out, nmb, ncells1);
  }
  auto part_in  = polar_part_in;
  auto part_out = polar_part_out;
  auto inner_   = inner_local;
  auto outer_   = outer_local;

  // =====================================================
  // stage 1: per-(MeshBlock, radius) partial sums over phi
  // =====================================================
  // One thread per (m,i) rather than one per i.  The old kernel launched ncells1 threads
  // -- 68 on this problem, well under a single MI300A workgroup -- and made each of them
  // walk every MeshBlock and every phi cell serially.  The sum over phi within one block
  // is left serial and in ascending k, exactly as before.
  par_for("polar_partial", DevExeSpace(), 0, nmb-1, 0, ncells1-1,
  KOKKOS_LAMBDA(const int m, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    if ((do_inner || do_outer) && (i >= is) && (i <= ie)) {
      for (int k=ks+1; k<=ke+1; ++k) {
        if (do_inner) inner_sum += e1(m,k,js,i);
        if (do_outer) outer_sum += e1(m,k,je+1,i);
      }
    }
    part_in(m,i)  = inner_sum;
    part_out(m,i) = outer_sum;
  });

  // =====================================================
  // stage 2: combine the per-block partials
  // =====================================================
  // Ascending block order, so the result does not depend on how the runtime scheduled
  // stage 1.  (It is NOT bit-identical to the old single-loop accumulation, which carried
  // one running total across block boundaries -- see the branch notes.)
  par_for("polar_combine", DevExeSpace(), 0, ncells1-1,
  KOKKOS_LAMBDA(const int i) {
    Real inner_sum = 0.0;
    Real outer_sum = 0.0;
    for (int m=0; m<nmb; ++m) {
      inner_sum += part_in(m,i);
      outer_sum += part_out(m,i);
    }
    inner_(i) = inner_sum;
    outer_(i) = outer_sum;
  });

#if MPI_PARALLEL_ENABLED
  // Only a genuinely distributed run needs the device->host->device round trip.  On one
  // rank the old code paid it anyway, twice per cycle, for an Allreduce over itself.
  if (global_variable::nranks > 1) {
    Kokkos::deep_copy(polar_inner_h, inner_local);
    Kokkos::deep_copy(polar_outer_h, outer_local);
    MPI_Allreduce(MPI_IN_PLACE, polar_inner_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, polar_outer_h.data(), ncells1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    Kokkos::deep_copy(inner_local, polar_inner_h);
    Kokkos::deep_copy(outer_local, polar_outer_h);
  }
#endif

  // =====================================================
  // refill polar EMFs with the azimuthal average
  // =====================================================
  // Parallel over (m,k,i); the old kernel ran one thread per MeshBlock and looped k and i
  // inside it.  Pure writes, so this is bit-identical to before.
  par_for("polar_refill", DevExeSpace(), 0, nmb-1, ks, ke+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    bool do_inner = (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::polar);
    bool do_outer = (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::polar);

    if (do_inner) e1(m,k,js  ,i) = inner_(i)/static_cast<Real>(mesh_nx3);
    if (do_outer) e1(m,k,je+1,i) = outer_(i)/static_cast<Real>(mesh_nx3);
  });
}

} // namespace mhd
