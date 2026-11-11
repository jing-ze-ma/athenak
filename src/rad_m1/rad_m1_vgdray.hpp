#ifndef RAD_M1_RAD_M1_VGDRAY_HPP_
#define RAD_M1_RAD_M1_VGDRAY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vgdray.hpp
//! \brief share3-1031: the ONE upwind core of the vet_gd ray solve, shared by the sp sweep
//! (rad_m1_vetgd_body.inc, VetGdSweep) and the cubed-sphere sweep (VetGdCsSweep): the
//! branch / chord geometry of a (cell, direction), the top-face inflow, the inner-face
//! value (vet_gd_inner = emit | diffusion), the upwind intensity over the stencil
//! (renormalised over the cells of the right branch) and the closing segment (vet_seg =
//! jump below the inner face, vet_interp = corner, vet_seg = split, the plain segment).
//!
//! The geometry differs only in the upwind STENCIL, passed as an accessor A:
//!   A.N()            the number of stencil cells
//!   A.K(e), A.J(e)   band position of cell e
//!   A.W(e)           its weight
//!   A.Interp(f)      the stencil interpolation of f(kk, jj) (sp: the nested bilinear
//!                    VgdLerp; cs: the weighted sum) -- each in its own former order
//! The field values come in as callables: fI(kk, jj) the upwind-shell intensity (sp: the
//! band / wall read rd; cs: the owner read rdi), fC(c, kk, jj, ii) ln chi (c = 0) or ln S
//! (c = 1) of a cell.  What stays separate: the stencil construction itself (sp: the
//! (theta, phi) bilinear foot with the band clamp; cs: VgdCsStencil on the panel map),
//! the sp-only branches (vet_gd_band_exit, vet_gd_theta_wall = vacuum, vet_path =
//! cells, the inner-face beam mask), and the band / halo layouts.  Every expression is
//! the former code's, term for term and in its order (bitwise).

#include <cmath>

#include "athena.hpp"
#include "rad_m1/rad_m1_vetseg.hpp"

//----------------------------------------------------------------------------------------
// the upwind shell iu, the chord ds back to it and the own-cell part dso (vet_seg =
// split: the path from the cell centre back to the radial face it shares with the upwind
// shell; typ 3: half the chord).  typ: 1 upwind shell, 2 inner face, 3 own shell (turned)
template <typename VX>
KOKKOS_INLINE_FUNCTION
void VgdRayChord(const bool inw, const int l, const int lcut, const int m, const int i,
                 const Real zr, const Real p2, const VX &cx1v, const VX &cx1f,
                 const bool spl, int &typ, int &iu, Real &ds, Real &dso) {
  typ = 1;
  iu = i + 1;
  if (inw) {
    const Real ru = cx1v(m,i+1);
    ds = sqrt(ru*ru - p2) - zr;
  } else {
    const Real rd = (l > lcut) ? cx1v(m,i-1) : cx1f(m,i);
    if (l > lcut && p2 <= rd*rd*(1.0 + 1.0e-13)) {
      iu = i - 1;
      ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
    } else if (l == lcut && p2 <= rd*rd) {
      typ = 2;
      iu = i;
      ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
    } else {
      typ = 3;
      iu = i;
      ds = 2.0*zr;
    }
  }
  dso = 0.5*ds;
  if (spl && typ == 1) {
    const Real rf = inw ? cx1f(m,i+1) : cx1f(m,i);
    dso = inw ? (sqrt(fmax(rf*rf - p2, 0.0)) - zr)
              : (zr - sqrt(fmax(rf*rf - p2, 0.0)));
    dso = fmin(fmax(dso, 0.0), ds);
  }
}

//----------------------------------------------------------------------------------------
// from the vacuum top face (the end cell's S extrapolated, as vet_col_lat), or with the
// top bath (vgdbath-1011) its incident intensity attenuated; sbl = S of the shell below
KOKKOS_INLINE_FUNCTION
Real VgdTopInflow(const Real ch0, const Real ds, const Real s0, const Real sbl,
                  const bool tbh, const Real itop) {
  const Real stp = fmax(1.5*s0 - 0.5*sbl, 1.0e-300);
  Real ex, w0, wu;
  VetSegW(ch0*ds, ex, w0, wu);
  if (tbh) {
    return itop*ex + wu*stp + w0*s0;
  } else {
    return wu*stp + w0*s0;
  }
}

//----------------------------------------------------------------------------------------
// typ 2, the inner face: emi (vet_gd_inner = emit): the boundary emission ie0 carried to
// the cell with the face source; else the diffusion intensity at the inner face of the
// first shell, S + mu (-dS/dtau_r) (bfc, vet_bc_diff = face: S of the FACE, 1.5 S_a -
// 0.5 S_b).  rf = the inner face radius, ra/rb = the first two shell radii
template <typename A, typename FC>
KOKKOS_INLINE_FUNCTION
Real VgdInnerFace(const A &acc, const FC &fC, const int i, const bool emi,
                  const Real ie0, const bool bfc, const Real rf, const Real ra,
                  const Real rb, const Real p2, const Real ch0, const Real ds,
                  const Real s0, const bool zf = false) {
  const int i1 = i + 1;
  auto sbt = [&](const int kk, const int jj) {
    return fmax(1.5*exp(fC(1,kk,jj,i)) - 0.5*exp(fC(1,kk,jj,i1)), 1.0e-300);
  };
  if (emi) {
    Real ex, w0, wu;
    VetSegW(ch0*ds, ex, w0, wu);
    const Real qv = acc.Interp(sbt);
    return ie0*ex + wu*qv + w0*s0;
  }
  const Real muf = sqrt(fmax(1.0 - p2/(rf*rf), 0.0));
  auto ib = [&](const int kk, const int jj) {
    const Real sa = exp(fC(1,kk,jj,i)), sb = exp(fC(1,kk,jj,i1));
    // hrsplit-1110 vet_gd_inner_reflect = zeroflux (zf): at a reflecting inner x1 the
    // diffusion start without the gradient term (H = 0 at the zero-flux face)
    const Real g0 = zf ? 0.0 : -(sb - sa)/(rb - ra)/exp(fC(0,kk,jj,i));
    return (bfc ? (1.5*sa - 0.5*sb) : sa) + g0*muf;
  };
  Real ex, w0, wu;
  VetSegW(ch0*ds, ex, w0, wu);
  const Real bv = fmax(acc.Interp(ib), 0.0);
  const Real qv = acc.Interp(sbt);
  return bv*ex + wu*qv + w0*s0;
}

//----------------------------------------------------------------------------------------
// the upwind intensity: each (cell, direction) holds ONE intensity, of the pass its
// n_d . r_hat selects; stencil cells of the other pass (ok false) are dropped and the
// weights renormalised (all valid: the plain interpolation); none: -1 (the closing
// segment then starts from the upwind source)
template <typename A, typename FI>
KOKKOS_INLINE_FUNCTION
Real VgdUpwindI(const A &acc, const FI &fI, const bool *ok, const bool all) {
  if (all) {return acc.Interp(fI);}
  Real ws = 0.0, vs = 0.0;
  for (int e = 0; e < acc.N(); ++e) {
    const Real we = acc.W(e);
    if (ok[e] && we > 0.0) {ws += we; vs += we*fI(acc.K(e), acc.J(e));}
  }
  return (ws > 0.0) ? vs/ws : -1.0;
}

//----------------------------------------------------------------------------------------
// the closing segment from the upwind point (I = ivu, chi = chu, S = su; ivu < 0: none)
// to the cell (ch0, s0) over the chord ds (own part dso):
//  * vet_seg = jump (spe 2) on a turned ray (typ 3, l > lcut) that dips below its own
//    inner face rfi: piecewise constant per cell (the chord below the face takes chi, S
//    of the cell below, log-mean of the own column and the upwind point);
//  * vet_interp = corner (vim 2, or 3 at a material edge, edgef()): one segment per
//    stencil cell from its own (I, chi, S), the outcome interpolated over the cells of
//    the right branch;
//  * vet_seg = split / the plain linear-in-tau segment.
// (k, j): the own column (its shell below, i - 1, for the jump branch)
template <typename A, typename FI, typename FC, typename FE>
KOKKOS_INLINE_FUNCTION
Real VgdSegFinish(const A &acc, const FI &fI, const FC &fC, const FE &edgef,
                  const bool *ok, const int iu, const int k, const int j,
                  const int i, const Real ivu, const Real chu, const Real su,
                  const Real ch0, const Real s0, const Real ds, const Real dso,
                  const bool spl, const int spe, const int typ, const int l,
                  const int lcut, const Real p2, const Real rfi, const Real zr,
                  const int vim) {
  if (spl && spe == 2 && typ == 3 && l > lcut && p2 < rfi*rfi) {
    auto lind = [&](const int c) {
      const Real a = acc.Interp([&](const int kk, const int jj) {
        return fC(c, kk, jj, i - 1);
      });
      return exp(0.5*(a + fC(c, k, j, i - 1)));
    };
    const Real zf = sqrt(fmax(rfi*rfi - p2, 0.0));
    const Real e1 = exp(-chu*(zr - zf)), e2 = exp(-lind(0)*2.0*zf);
    const Real e3 = exp(-ch0*(zr - zf));
    const Real ia = ((ivu < 0.0) ? su : ivu)*e1 + su*(1.0 - e1);
    return fmax((ia*e2 + lind(1)*(1.0 - e2))*e3 + s0*(1.0 - e3), 0.0);
  } else if (vim == 2 || (vim == 3 && edgef())) {
    Real ws = 0.0, vs = 0.0;
    for (int e = 0; e < acc.N(); ++e) {
      const Real we = acc.W(e);
      if (!(we > 0.0) || !(ok[e] || ivu < 0.0)) {continue;}
      const int kk = acc.K(e), jj = acc.J(e);
      const Real cq = exp(fC(0,kk,jj,iu)), sq = exp(fC(1,kk,jj,iu));
      const Real iq = (ivu < 0.0) ? sq : fI(kk, jj);
      Real vq;
      if (spl) {
        vq = VetSegSplit(iq, sq, s0, cq*(ds - dso), ch0*dso, spe);
      } else {
        Real ex, w0, wu;
        VetSegW(0.5*(cq + ch0)*ds, ex, w0, wu);
        vq = iq*ex + wu*sq + w0*s0;
      }
      ws += we;
      vs += we*vq;
    }
    return fmax(vs/ws, 0.0);
  } else if (spl) {
    return fmax(VetSegSplit((ivu < 0.0) ? su : ivu, su, s0, chu*(ds - dso), ch0*dso, spe),
                0.0);
  } else {
    Real ex, w0, wu;
    VetSegW(0.5*(chu + ch0)*ds, ex, w0, wu);
    return fmax(((ivu < 0.0) ? su : ivu)*ex + wu*su + w0*s0, 0.0);
  }
}

//----------------------------------------------------------------------------------------
// vetgdup-1110 vet_gd_upwind: the shell's mean angular profile Ibar(mu) of shell iu on
// the half range of the branch (lo: mu < 0; 2 kh bins over [-1, 1], width 1/kh; bin b
// holds a0 + a1 t + a2 t^2, t = (mu - c_b) kh, flag pf(iu,b,3) > 0 when fitted); between
// bin centres the linear blend of the two bins' quadratics (smooth: its second derivative
// carries the fits' own curvature, not the O(1) jump of a piecewise-linear profile).
// Returns false when a bin it needs holds no fit.  f1, f2: d/dmu, d^2/dmu^2
template <typename P>
KOKKOS_INLINE_FUNCTION
bool VgdUpProf(const P &pf, const int iu, const int kh, const bool lo, const Real mu,
               const Real tol, Real &f, Real &f1, Real &f2) {
  const Real mc = lo ? fmin(fmax(mu, -1.0), 0.0) : fmin(fmax(mu, 0.0), 1.0);
  const Real p = (mc + 1.0)*kh - 0.5;
  const int h0 = lo ? 0 : kh, h1 = lo ? (kh - 1) : (2*kh - 1);
  int b0 = static_cast<int>(floor(p));
  Real w = p - b0;
  int b1 = b0 + 1;
  if (b0 < h0) {b0 = h0; w = 1.0;}
  if (b1 > h1) {b1 = h1; w = 0.0;}
  if (b0 > h1) {b0 = h1;}
  if (b1 < h0) {b1 = h0;}
  Real q[2], q1[2], q2[2];
  for (int s = 0; s < 2; ++s) {
    const int b = (s == 0) ? b0 : b1;
    if (!(pf(iu,b,3) > 0.0)) {return false;}
    const Real t = (mc + 1.0)*kh - (b + 0.5);
    const Real a0 = pf(iu,b,0), a1 = pf(iu,b,1), a2 = pf(iu,b,2);
    q[s] = a0 + t*(a1 + t*a2);
    q1[s] = (a1 + 2.0*a2*t)*kh;
    q2[s] = 2.0*a2*kh*kh;
  }
  // the two bins' fits must agree at mu (else the profile does not resolve the
  // angular structure there, e.g. a narrow beam: no shift)
  if (b1 > b0 && fabs(q[1] - q[0]) > tol*fmax(fmax(fabs(q[0]), fabs(q[1])), 1.0e-300)) {
    return false;
  }
  const Real dw = (b1 > b0) ? static_cast<Real>(kh) : 0.0;
  f = (1.0 - w)*q[0] + w*q[1];
  f1 = (1.0 - w)*q1[0] + w*q1[1] + (q[1] - q[0])*dw;
  f2 = (1.0 - w)*q2[0] + w*q2[1] + 2.0*(q1[1] - q1[0])*dw;
  return true;
}

// the additive corrections cq[e] of the stencil cells (cosines mq(e), intensities iq(e),
// used where ok[e] or all) for the upwind point of cosine mu: upf 1 Ibar(mu) - Ibar(mq),
// 2 the Taylor terms of Ibar about mq to second order.  false (no correction anywhere)
// when a bin it needs holds no fit, or when the profile does not reproduce a stencil
// cell's own intensity to tol (|iq - Ibar(mq)| > tol max(|iq|, |Ibar|): the angular
// structure is not resolved there, e.g. the edge of a beam; the read stays bilinear)
template <typename P, typename A, typename MQ, typename IQ>
KOKKOS_INLINE_FUNCTION
bool VgdUpCorr(const P &pf, const int iu, const int kh, const bool lo, const int upf,
               const Real tol, const Real mu, const A &acc, const MQ &mq, const IQ &iq,
               const bool *ok, const bool all, Real *cq) {
  Real fu = 0.0, f1, f2;
  if (upf == 1 && !VgdUpProf(pf, iu, kh, lo, mu, tol, fu, f1, f2)) {return false;}
  for (int e = 0; e < acc.N(); ++e) {
    cq[e] = 0.0;
    // only the cells the upwind read itself reads (all: every stencil cell; else the
    // valid ones of nonzero weight): no read outside the halo read masks
    if (!(all || (ok[e] && acc.W(e) > 0.0))) {continue;}
    const Real m = mq(e);
    Real f;
    if (!VgdUpProf(pf, iu, kh, lo, m, tol, f, f1, f2)) {return false;}
    const Real v = iq(e);
    if (fabs(v - f) > tol*fmax(fmax(fabs(v), fabs(f)), 1.0e-300)) {return false;}
    const Real dm = mu - m;
    cq[e] = (upf == 1) ? (fu - f) : dm*(f1 + 0.5*dm*f2);
  }
  return true;
}

#endif  // RAD_M1_RAD_M1_VGDRAY_HPP_
