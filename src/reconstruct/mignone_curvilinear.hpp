#ifndef RECONSTRUCT_MIGNONE_CURVILINEAR_HPP_
#define RECONSTRUCT_MIGNONE_CURVILINEAR_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mignone_curvilinear.hpp
//! \brief Higher-order reconstruction on a genuinely NON-UNIFORM radial (x1) grid on
//! the cubed sphere and on spherical polar (x1 = spherical radius r on both grids, per
//! the radial unification -- see coordinates/grid_stretch.hpp:RadialCentroid, which is
//! exactly what x1v stores on both). This follows
//!
//!   A. Mignone, "High-order conservative reconstruction schemes for finite volume
//!   methods in cylindrical and spherical coordinates", J. Comput. Phys. 270, 784
//!   (2014). arXiv:1404.0537.
//!
//! Equation numbers in the comments below refer to that paper. The key idea (Mignone
//! Sec. 2.2) is to reconstruct directly from VOLUME averages using the exact
//! geometry-dependent moments (his Eq. 23) of a power-law Jacobian dV/dxi = xi^m
//! (m=2 for a spherical radial coordinate), rather than assuming uniform spacing or
//! reconstructing blindly in a grid-index coordinate. Eq. (23) is closed-form and
//! EXACT for an arbitrary (non-uniform) face layout, so no assumption of uniform
//! Delta-r is made anywhere in this file; Mignone's own remark 5 (Sec. 2.2) endorses
//! solving his Eq. (21) by direct numerical inversion whenever the grid is
//! non-uniform, which is what MignoneWeights3/4 below do (Cramer's rule, since the
//! systems are only 3x3 / 4x4).
//!
//! PLM is NOT reimplemented here: the existing GridPiecewiseLinearX1 (hydro.hpp /
//! mhd.hpp), built around PLM_nonuniform's modified van-Leer limiter, is ALREADY
//! Mignone's curvilinear PLM (his Eq. 30 with the modified van-Leer limiter, Eq. 37,
//! using cF_i/cB_i = his Eq. 33) evaluated with x1v = his volume centroid Eq. (17):
//! PLM_nonuniform's cL = dxL/dxLh and cR = dxR/dxRh reduce exactly to his cB_i, cF_i
//! once x1v is the volume centroid, which RadialCentroid() confirms it is. That is
//! why the default (plm) x1 reconstruction on cs/sp is left bitwise unchanged.
//!
//! What IS implemented here, for <mhd|hydro>/reconstruct_x1:
//!  * wenoz -> Mignone's WENO3 (Sec. 3.2, third order). This is a REDUCTION from the
//!    5th-order uniform-grid WENO-Z used elsewhere in this code (reconstruct/wenoz.hpp)
//!    -- Mignone (2014) gives no curvilinear generalization of a higher-order WENO,
//!    so third order is what the published, non-invented scheme supports. Verify the
//!    measured order in tests_cs_radial_recon/.
//!  * ppm4 -> Mignone's curvilinear PPM (Sec. 3.3): the unlimited 4th-order interface
//!    value from his Eq. (20)/(21) (the face-centered 4-point stencil, matching the
//!    stencil already used by the uniform-grid PPM4 in reconstruct/ppm.hpp), followed
//!    by his curvilinear monotonicity limiter (Eqs. 45-48) with the closed-form h_i^pm
//!    for a spherical radial coordinate (Eq. 48, "spherical, xi = r" case).
//!  * ppmx (the Colella & Sekora extremum-preserving limiter) is REFUSED for x1 on
//!    cs/sp: Mignone (2014) only curvilinearizes the CW84-style limiter (Eqs. 45-48),
//!    not a Colella-Sekora extremum-preserving one, and inventing that generalization
//!    is out of scope here. See the fatal in hydro.cpp / mhd.cpp.

#include <math.h>
#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn MignoneBeta()
//! \brief Eq. (23): beta_{i+s,n} for a power-law Jacobian dV/dxi = xi^m, with xi_c = 0
//! (remark 2: the final interface values do not depend on this choice). EXACT for any
//! (non-uniform) face pair [xm, xp]; m=2 is the spherical radial case used throughout
//! this file (cs and sp both use x1 = spherical r).

KOKKOS_INLINE_FUNCTION
Real MignoneBeta(const int n, const Real xm, const Real xp, const int m) {
  const Real pw = static_cast<Real>(n + m + 1);
  const Real pm = static_cast<Real>(m + 1);
  return pm/pw * (pow(xp, pw) - pow(xm, pw)) / (pow(xp, pm) - pow(xm, pm));
}

//----------------------------------------------------------------------------------------
//! \fn MignoneWeights3()
//! \brief Direct (Cramer's rule) solve of Eq. (21) for p=3 (iL=iR=1, stencil
//! {i-1,i,i+1}), evaluated at a single target abscissa xt. xf0..xf3 are the FOUR face
//! positions bounding the three cells (xf0,xf1 bound cell i-1; xf1,xf2 bound cell i;
//! xf2,xf3 bound cell i+1). Returns the three weights w[-1],w[0],w[+1] such that
//! Q(xt) = w_{-1}<Q>_{i-1} + w_0<Q>_i + w_{+1}<Q>_{i+1}.

KOKKOS_INLINE_FUNCTION
void MignoneWeights3(const Real xf0, const Real xf1, const Real xf2, const Real xf3,
                      const Real xt, const int m,
                      Real &wm1, Real &w0, Real &wp1) {
  // beta matrix rows (s=-1,0,1), columns n=0,1,2 -- Eq. (21) is B^T w = rhs
  const Real b00 = 1.0;
  const Real b01 = MignoneBeta(1, xf0, xf1, m), b02 = MignoneBeta(2, xf0, xf1, m);
  const Real b10 = 1.0;
  const Real b11 = MignoneBeta(1, xf1, xf2, m), b12 = MignoneBeta(2, xf1, xf2, m);
  const Real b20 = 1.0;
  const Real b21 = MignoneBeta(1, xf2, xf3, m), b22 = MignoneBeta(2, xf2, xf3, m);
  // B^T
  const Real a00=b00, a01=b10, a02=b20;
  const Real a10=b01, a11=b11, a12=b21;
  const Real a20=b02, a21=b12, a22=b22;
  const Real r0 = 1.0, r1 = xt, r2 = xt*xt;
  auto det3 = [](Real m00, Real m01, Real m02, Real m10, Real m11, Real m12,
                 Real m20, Real m21, Real m22) {
    return m00*(m11*m22-m12*m21) - m01*(m10*m22-m12*m20) + m02*(m10*m21-m11*m20);
  };
  const Real D  = det3(a00,a01,a02, a10,a11,a12, a20,a21,a22);
  const Real D0 = det3(r0,a01,a02,  r1,a11,a12,  r2,a21,a22);
  const Real D1 = det3(a00,r0,a02,  a10,r1,a12,  a20,r2,a22);
  const Real D2 = det3(a00,a01,r0,  a10,a11,r1,  a20,a21,r2);
  wm1 = D0/D; w0 = D1/D; wp1 = D2/D;
}

//----------------------------------------------------------------------------------------
//! \fn MignoneWeights4()
//! \brief Direct (Cramer's rule) solve of Eq. (21) for p=4, the FACE-centered stencil
//! {i-1,i,i+1,i+2} used by PPM (remark 1: the interface value at x_{i+1/2} is computed
//! ONCE and shared as Q+_i and Q-_{i+1}). xf0..xf4 are the five face positions bounding
//! the four cells. The target is normally one of the interior faces, xf2 (=x_{i+1/2}).

KOKKOS_INLINE_FUNCTION
void MignoneWeights4(const Real xf0, const Real xf1, const Real xf2, const Real xf3,
                      const Real xf4, const Real xt, const int m,
                      Real &wm1, Real &w0, Real &wp1, Real &wp2) {
  Real b[4][4];
  const Real xl[4] = {xf0, xf1, xf2, xf3};
  const Real xr[4] = {xf1, xf2, xf3, xf4};
  for (int s=0; s<4; ++s) {
    b[s][0] = 1.0;
    b[s][1] = MignoneBeta(1, xl[s], xr[s], m);
    b[s][2] = MignoneBeta(2, xl[s], xr[s], m);
    b[s][3] = MignoneBeta(3, xl[s], xr[s], m);
  }
  // A = B^T (4x4), rhs = [1, xt, xt^2, xt^3]
  Real a[4][4];
  for (int r=0; r<4; ++r) { for (int c=0; c<4; ++c) { a[r][c] = b[c][r]; } }
  const Real rhs[4] = {1.0, xt, xt*xt, xt*xt*xt};

  // 4x4 determinant via cofactor expansion along the first row, reused for Cramer.
  auto det3 = [](Real m00, Real m01, Real m02, Real m10, Real m11, Real m12,
                 Real m20, Real m21, Real m22) {
    return m00*(m11*m22-m12*m21) - m01*(m10*m22-m12*m20) + m02*(m10*m21-m11*m20);
  };
  auto det4 = [&](Real M[4][4]) {
    Real d = 0.0;
    for (int c=0; c<4; ++c) {
      Real minor[3][3]; int mc=0;
      for (int cc=0; cc<4; ++cc) {
        if (cc==c) continue;
        for (int rr=1; rr<4; ++rr) { minor[rr-1][mc] = M[rr][cc]; }
        mc++;
      }
      Real cof = det3(minor[0][0],minor[0][1],minor[0][2],
                       minor[1][0],minor[1][1],minor[1][2],
                       minor[2][0],minor[2][1],minor[2][2]);
      d += ((c%2==0) ? 1.0 : -1.0) * M[0][c] * cof;
    }
    return d;
  };
  const Real D = det4(a);
  Real w[4];
  for (int col=0; col<4; ++col) {
    Real Ac[4][4];
    for (int r=0; r<4; ++r) { for (int c=0; c<4; ++c) { Ac[r][c] = a[r][c]; } }
    for (int r=0; r<4; ++r) { Ac[r][col] = rhs[r]; }
    w[col] = det4(Ac)/D;
  }
  wm1 = w[0]; w0 = w[1]; wp1 = w[2]; wp2 = w[3];
}

//----------------------------------------------------------------------------------------
//! \fn MignoneWENO3()
//! \brief Third-order curvilinear WENO reconstruction, Mignone Sec. 3.2 (Eqs. 28-29,
//! 39-42). q_{im1,i,ip1} are the cell averages; x_{im1,i,ip1} the volume CENTROIDS
//! (x1v); xfm,xfp the two faces of cell i (x_{i-1/2}, x_{i+1/2}); dxi = xfp-xfm.
//! dp0, dm0 are the geometry-only linear weights d+_{i,0}, d-_{i,0}: derived by
//! matching the coefficient of <Q>_{i+1} (resp. <Q>_{i-1}) between the exact 3-point
//! (p=3) reconstruction (MignoneWeights3) and the d*Q[f]+(1-d)*Q[b] decomposition of
//! Eq. (39); this reproduces Eq. (42) without depending on its exact typeset form.

KOKKOS_INLINE_FUNCTION
void MignoneWENO3(const Real q_im1, const Real q_i, const Real q_ip1,
                   const Real x_im1, const Real x_i, const Real x_ip1,
                   const Real xfm, const Real xfp, const Real dxi,
                   const Real dp0, const Real dm0,
                   Real &ql_ip1, Real &qr_i) {
  // Eq. (29): forward/backward undivided differences (exact for a linear function on
  // any grid spacing).
  const Real dQF = dxi*(q_ip1 - q_i)/(x_ip1 - x_i);
  const Real dQB = dxi*(q_i - q_im1)/(x_i - x_im1);
  // Eq. (28): the four 2-point extrapolated interface states
  const Real Qpf = q_i + dQF*(xfp - x_i)/dxi;   // Q+,[f]_i
  const Real Qpb = q_i + dQB*(xfp - x_i)/dxi;   // Q+,[b]_i
  const Real Qmf = q_i + dQF*(xfm - x_i)/dxi;   // Q-,[f]_i
  const Real Qmb = q_i + dQB*(xfm - x_i)/dxi;   // Q-,[b]_i
  // Eq. (40)-(41): Yamaleev-Carpenter-style nonlinear weights, curvilinear version.
  // Cref = 20, as used throughout Mignone (2014).
  const Real tau = fabs(dQF - dQB);
  const Real beta0 = dQF*dQF, beta1 = dQB*dQB;
  const Real Qref = 20.0*fmax(fabs(q_im1), fmax(fabs(q_i), fabs(q_ip1)));
  const Real refsq = Qref*Qref;
  const Real tausq = tau*tau;
  // Regularization epsilon: for a field component that is identically uniform (often
  // exactly zero, e.g. one gnomonic component of an otherwise nonzero uniform field),
  // tau, beta0, beta1 and Qref all vanish TOGETHER, and the ratio tausq/(beta+refsq)
  // becomes 0/0. The same regularization is already used in the uniform-grid WENOZ()
  // (reconstruct/wenoz.hpp, epsL = 1e-42); the correct limit here is omega -> d (the
  // linear weights), which this epsilon reproduces without perturbing a resolved,
  // nonzero profile.
  const Real weno_eps = 1.0e-40;
  const Real ap0 = dp0*(1.0 + tausq/(beta0 + refsq + weno_eps));
  const Real ap1 = (1.0-dp0)*(1.0 + tausq/(beta1 + refsq + weno_eps));
  const Real am0 = dm0*(1.0 + tausq/(beta0 + refsq + weno_eps));
  const Real am1 = (1.0-dm0)*(1.0 + tausq/(beta1 + refsq + weno_eps));
  ql_ip1 = (ap0*Qpf + ap1*Qpb)/(ap0+ap1);
  qr_i   = (am0*Qmf + am1*Qmb)/(am0+am1);
}

//----------------------------------------------------------------------------------------
//! \fn MignonePPM4()
//! \brief Curvilinear 4th-order PPM (Mignone Sec. 3.3): unlimited face values from the
//! face-centered 4-point stencil (MignoneWeights4), clipped to neighbouring bounds
//! (Eq. 45), then the CW84-style monotonicity limiter generalized with the closed-form
//! h_i^pm for a spherical radial coordinate (Eqs. 46, 48; ratio_p/ratio_m below reduce
//! to the classic factor of 2 -- see the existing PPM4()'s "2.0*qd"/"2.0*qc" -- exactly
//! when h+ = h- = 3, its Cartesian value). fm2..fp3 are the SIX face positions bounding
//! cells {i-2,i-1,i,i+1,i+2}, i.e. fm2=x_{i-5/2} .. fp3=x_{i+5/2}, matching the call
//! convention of the existing uniform PPM4/PPMX (q_im2..q_ip2 around cell i). r_i, dr_i
//! in Eq. (48) are taken to be the volume centroid and face-to-face width of cell i,
//! consistent with how a "cell radius" is used everywhere else in this scheme.

KOKKOS_INLINE_FUNCTION
void MignonePPM4(const Real q_im2, const Real q_im1, const Real q_i, const Real q_ip1,
                  const Real q_ip2,
                  const Real fm2, const Real fm1, const Real f0, const Real fp1,
                  const Real fp2, const Real fp3, const Real xv_i, const int m,
                  Real &ql_ip1, Real &qr_i) {
  Real wa, wb, wc, wd;
  // Q-_i: face x_{i-1/2} = f0, face-centered stencil {i-2,i-1,i,i+1}
  MignoneWeights4(fm2, fm1, f0, fp1, fp2, f0, m, wa, wb, wc, wd);
  Real Qm_i = wa*q_im2 + wb*q_im1 + wc*q_i + wd*q_ip1;
  // Q+_i: face x_{i+1/2} = fp1, face-centered stencil {i-1,i,i+1,i+2}
  MignoneWeights4(fm1, f0, fp1, fp2, fp3, fp1, m, wa, wb, wc, wd);
  Real Qp_i = wa*q_im1 + wb*q_i + wc*q_ip1 + wd*q_ip2;

  // Eq. (45): clip to the bounds set by the immediate neighbours.
  Qp_i = fmin(Qp_i, fmax(q_i, q_ip1));  Qp_i = fmax(Qp_i, fmin(q_i, q_ip1));
  Qm_i = fmin(Qm_i, fmax(q_i, q_im1));  Qm_i = fmax(Qm_i, fmin(q_i, q_im1));

  // Eq. (48), spherical radial case: h_i^pm = 3 + 2*dr_i*(pm 10*r_i + dr_i)/
  // (20*r_i^2 + dr_i^2); "spherical, xi=r" in the paper. m selects the case so the
  // same routine could later be extended to m=0,1; only m=2 is exercised here.
  const Real dr_i = fp1 - f0;
  Real hp, hm;
  if (m == 2) {
    const Real denom = 20.0*xv_i*xv_i + dr_i*dr_i;
    hp = 3.0 + 2.0*dr_i*( 10.0*xv_i + dr_i)/denom;
    hm = 3.0 + 2.0*dr_i*(-10.0*xv_i + dr_i)/denom;
  } else {
    hp = 3.0; hm = 3.0;  // Cartesian (m=0); unused in this file but kept for safety
  }

  const Real dQp = Qp_i - q_i;
  const Real dQm = Qm_i - q_i;
  if (dQp*dQm >= 0.0) {
    Qp_i = q_i;
    Qm_i = q_i;
  } else {
    const Real ratio_p = (hm + 1.0)/(hp - 1.0);
    const Real ratio_m = (hp + 1.0)/(hm - 1.0);
    if (fabs(dQp) >= ratio_p*fabs(dQm)) { Qp_i = q_i - ratio_p*dQm; }
    if (fabs(dQm) >= ratio_m*fabs(dQp)) { Qm_i = q_i - ratio_m*dQp; }
  }

  ql_ip1 = Qp_i;
  qr_i   = Qm_i;
}

//----------------------------------------------------------------------------------------
//! \fn MignoneWENO3X1() / MignonePPM4X1()
//! \brief Wrapper functions for the radial (x1) sweep on cs/sp, called over
//! [is-1,ie+1] to get BOTH L/R states over [is,ie] (same convention as PLM/PPM/WENOZ
//! X1 wrappers in reconstruct/). m_geom = 2 (spherical r) on both grids.

KOKKOS_INLINE_FUNCTION
void MignoneWENO3X1(TeamMember_t const &member, const EOS_Data &eos,
     const bool apply_floors,
     const int m, const int k, const int j, const int il, const int iu,
     const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv,
     const DvceArray2D<Real> &xf, ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
  int nvar = q.extent_int(1);
  const Real &dfloor_ = eos.dfloor;
  const int m_geom = 2;
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      const Real xim1 = xv(m,i-1), xi = xv(m,i), xip1 = xv(m,i+1);
      const Real xfm = xf(m,i), xfp = xf(m,i+1);
      const Real dxi = xfp - xfm;
      Real wm1, w0, wp1;
      // d+_{i,0}: match the coefficient of <Q>_{i+1} between the exact 3-point
      // reconstruction at x_{i+1/2} and the Eq.(39) decomposition (see file header).
      MignoneWeights3(xf(m,i-1), xf(m,i), xf(m,i+1), xf(m,i+2), xfp, m_geom,
                       wm1, w0, wp1);
      const Real cF = (xip1 - xi)/(xfp - xi);
      const Real dp0 = wp1*cF;
      // d-_{i,0}: match the coefficient of <Q>_{i-1} at x_{i-1/2}.
      MignoneWeights3(xf(m,i-1), xf(m,i), xf(m,i+1), xf(m,i+2), xfm, m_geom,
                       wm1, w0, wp1);
      const Real cB = (xi - xim1)/(xi - xfm);
      const Real dm0 = 1.0 - wm1*cB;
      MignoneWENO3(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), xim1, xi, xip1,
                   xfm, xfp, dxi, dp0, dm0, ql(n,i+1), qr(n,i));
      if (apply_floors) {
        if (n==IDN) {
          ql(IDN,i+1) = fmax(ql(IDN,i+1), dfloor_);
          qr(IDN,i  ) = fmax(qr(IDN,i  ), dfloor_);
        }
        if (n==IEN) {
          eos.ApplyEnergyFloor(ql(IDN,i+1), ql(IEN,i+1));
          eos.ApplyEnergyFloor(qr(IDN,i  ), qr(IEN,i  ));
        }
      }
    });
  }
  return;
}

KOKKOS_INLINE_FUNCTION
void MignonePPM4X1(TeamMember_t const &member, const EOS_Data &eos,
     const bool apply_floors,
     const int m, const int k, const int j, const int il, const int iu,
     const DvceArray5D<Real> &q, const DvceArray2D<Real> &xv,
     const DvceArray2D<Real> &xf, ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
  int nvar = q.extent_int(1);
  const Real &dfloor_ = eos.dfloor;
  const int m_geom = 2;
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      MignonePPM4(q(m,n,k,j,i-2), q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1),
                  q(m,n,k,j,i+2),
                  xf(m,i-2), xf(m,i-1), xf(m,i), xf(m,i+1), xf(m,i+2), xf(m,i+3),
                  xv(m,i), m_geom, ql(n,i+1), qr(n,i));
      if (apply_floors) {
        if (n==IDN) {
          ql(IDN,i+1) = fmax(ql(IDN,i+1), dfloor_);
          qr(IDN,i  ) = fmax(qr(IDN,i  ), dfloor_);
        }
        if (n==IEN) {
          eos.ApplyEnergyFloor(ql(IDN,i+1), ql(IEN,i+1));
          eos.ApplyEnergyFloor(qr(IDN,i  ), qr(IEN,i  ));
        }
      }
    });
  }
  return;
}
#endif  // RECONSTRUCT_MIGNONE_CURVILINEAR_HPP_
