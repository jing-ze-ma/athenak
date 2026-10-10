//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_implicit.cpp
//! \brief IMPLICIT transport along x1 columns, <rad_m1>/transport = implicit_x1
//! (milestone 3a of docs/dev/rad_m1_implicit_design.md).
//!
//! ONE backward-Euler solve per hydro step, at the TRUE speed of light, no sub-cycling
//! and no PD-ARS.  The normal flux lives on the x1 FACES and is eliminated there, which
//! turns the coupled (E, F) system into one tridiagonal M-matrix system for E per column:
//!
//!   F0'_f = theta_f [ F0^n_f - chat c dt (w_i E'_i - w_{i-1} E'_{i-1})/dx
//!                     - chat dt v_f g0_f ],    theta_f = 1/(1 + chat dt (rho k_t)_f)
//!   E'_i + (dt/dx)(chat/c)[ (F0' + A')_{i+1/2} - (F0' + A')_{i-1/2} ]
//!        = E^n_i + dt chat (rho kappa_P a T'^4 - rho kappa_E E0'_i)
//!
//! with w = P_11/E from the LAGGED closure (= chi in 1-D), A = a E the enthalpy flux
//! (a = v1 (1 + w)) upwinded with the face velocity, and the emission term linearised in
//! T and eliminated into the diagonal (design sect. 2).  Both the diffusion and the
//! upwind-advection off-diagonals are non-positive and the source adds
//! dt chat rho kappa_E rho c_v/B >= 0 to the diagonal, so the matrix is an M-matrix and
//! E' > 0 at any dt.
//!
//! OUTER LOOP: Picard on (w, theta, a, de0, g0, the upwind directions, optionally the
//! opacities) plus the safeguarded scalar root find for T' at fixed E'.  At convergence
//! the linearised emission correction vanishes identically, so the converged state
//! satisfies the NONLINEAR backward-Euler equations (see the note above
//! M1ImplTemperature).
//!
//! RESTRICTIONS of 3a, all fatal:
//!   * exactly ONE MeshBlock along x1 (each column is solved by a plain Thomas sweep
//!     inside one block; the partitioned line solve of two_stream_column_partition.hpp
//!     is NOT used);
//!   * nx2 = nx3 = 1 unless <rad_m1>/implicit_allow_multid = true, in which case the
//!     columns are INDEPENDENT: there is no transport along x2/x3 in this mode and
//!     F_2 = F_3 = 0 always;
//!   * no SMR/AMR.

#include <float.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1FLev
//! \brief blendall-1009: the reduced flux |f| the M1 (Levermore) relation assigns to an
//! Eddington factor chi, chi = (3 + 4 f^2)/(5 + 2 sqrt(4 - 3 f^2)) inverted:
//! sqrt(4 - 3 f^2) = (5 - 3 chi)/2, chi clamped to [1/3, 1] (the sp-blend-1008 f(D_rr)).

KOKKOS_INLINE_FUNCTION
Real M1FLev(const Real chi) {
  const Real cq = fmin(fmax(chi, 1.0/3.0), 1.0);
  const Real sq = 0.5*(5.0 - 3.0*cq);
  return sqrt(fmax((4.0 - sq*sq)/3.0, 0.0));
}

//----------------------------------------------------------------------------------------
//! \fn M1BeamFace
//! \brief blendall-1009 (implicit_flux_faces = all): the berthon / blend coefficients of
//! ONE face of the multi-D solve, from the BEAM VECTOR f = |f| n of its two cells: n the
//! axis of the LAGGED uniaxial closure (closure = m1: the lagged flux direction; vet_sc:
//! the formal solution's H/|H|), |f| = f_Lev(chi) (implicit_flux_beam = closure) or the
//! formal solution's |H|/J (= fs, vet_sc only).
//!   fl0, fr0: |f| of the two cells;
//!   nl, nr: the face-normal component n_d of the two cells (mu of the wave speeds),
//!   tauf:   the face optical depth (rho kappa_T)_f dx_f.
//! Wave speeds: the multi-D M1 eigenvalues lam(|f|, mu = n_d) (Skinner & Ostriker); the
//! advective part carries the NORMAL component f_d = |f| n_d; alpha (Bloch eq. 25) takes
//! the mean |f_d|; the blend weight w_f takes |f| (fmode max/mean), i.e. the beam measure
//! of the closure, not of its normal projection, so an oblique beam is a beam on every
//! face it crosses.  Exactly the 1-D formulas otherwise (x1 m1_impl_aphll): berthon =
//! both
//! HLL parts weighted alpha, the fmax/fmin M-matrix guards, AL = w_f.  A beam parallel to
//! the face (|f| = 1, mu = 0) has lam = 0 on both sides: zero face flux, as it should.

KOKKOS_INLINE_FUNCTION
void M1BeamFace(const Real fl0, const Real nl, const Real fr0, const Real nr,
                const Real tauf, const bool edd, const bool blend, const int bkind,
                const int bfm, const Real btau0, const Real bflo, const Real bfhi,
                const Real ch, Real &alw, Real &ccl, Real &ccr) {
  const Real fl = edd ? 0.0 : fmin(fmax(fl0, 0.0), 1.0);
  const Real fr = edd ? 0.0 : fmin(fmax(fr0, 0.0), 1.0);
  const Real fdl = fl*nl, fdr = fr*nr;
  Real bl, br;
  if (edd) {
    br = ch/sqrt(3.0);
    bl = -br;
  } else {
    Real lml, lpl, lmr, lpr;
    M1WaveSpeeds(fl, fmin(fmax(nl, -1.0), 1.0), lml, lpl);
    M1WaveSpeeds(fr, fmin(fmax(nr, -1.0), 1.0), lmr, lpr);
    bl = ch*fmin(fmin(lml, lmr), 0.0);
    br = ch*fmax(fmax(lpl, lpr), 0.0);
  }
  Real al = 1.0;
  if (tauf > 0.0) {
    const Real fbar = 0.5*(fabs(fdl) + fabs(fdr));
    const Real guard = fmax(1.0 - fbar*fbar, 0.0);
    const Real lp = br/ch, lm = bl/ch;
    const Real den = 1.0 - 3.0*tauf*guard*lp*lm/(lp - lm + 1.0e-300);
    al = 1.0/fmax(den, 1.0);
  }
  const Real invb = 1.0/(br - bl + 1.0e-300);
  const Real adl = br*ch*fdl*invb;
  const Real adr = -bl*ch*fdr*invb;
  const Real dk = -br*bl*invb;
  Real wf = 1.0;
  if (blend) {
    wf = M1BlendWeight(bkind, bfm, tauf, btau0, fl, fr, bflo, bfhi);
  }
  ccl = wf*fmax(al*adl + al*dk, 0.0);
  ccr = wf*fmin(al*adr - al*dk, 0.0);
  alw = wf;
}

//----------------------------------------------------------------------------------------
//! \fn M1HrCell
//! \brief hrup-1009 (implicit_flux_beam = halfrange): the HALF-RANGE ratios hp = H+/J >=
//! 0
//! and hm = H-/J <= 0 of cell (m,k,j,i) along axis d (0, 1, 2), hp + hm = H_d/J.
//!  src 1 (vet_sc): the ray sums of the sweep, sum w I max(mu_d, 0) (M1_VET_HP1+d) and
//!    H_d (M1_VET_H1+d), over J;
//!  src 2 (vet_gd, sp): vgd_hr (twin-corrected); hp < 0 there flags "no formal solution"
//!    (below the vet_gd cut), which the caller turns into a central face;
//!  src 0 (no ray data): the ISOTROPIC + BEAM model of the lagged closure, I = (1 - f) J
//!    isotropic + f J along n: hp = (1 - f)/4 + f max(n_d, 0),
//!    hm = -(1 - f)/4 + f min(n_d, 0), with f = f_Lev(chi), n the closure axis; exact in
//! both limits (f = 0: +-1/4;
//!    f = 1: the beam) and hp + hm = f n_d (the flux of the closure).
//!  src 3 (1-D x1 solve): the same model with |f| and the sign of n_1 from fx.

template <class V, class VH>
KOKKOS_INLINE_FUNCTION
void M1HrCell(const int src, const V &iw, const V &vc, const VH &gh, const int m,
              const int k, const int j, const int i, const int d, const Real fx,
              Real &hp, Real &hm, Real &fb) {
  if (src == 1) {
    const Real jj = vc(m,M1_VET_J,k,j,i);
    if (jj > 0.0) {
      const Real hq = vc(m,M1_VET_HP1+d,k,j,i);
      const Real hh = vc(m,M1_VET_H1+d,k,j,i);
      hp = fmin(fmax(hq/jj, 0.0), 1.0);
      hm = fmin(fmax((hh - hq)/jj, -1.0), 0.0);
      const Real h1 = vc(m,M1_VET_H1,k,j,i), h2 = vc(m,M1_VET_H1+1,k,j,i);
      const Real h3 = vc(m,M1_VET_H1+2,k,j,i);
      fb = fmin(sqrt(h1*h1 + h2*h2 + h3*h3)/jj, 1.0);
      return;
    }
    hp = 0.25;
    hm = -0.25;
    fb = 0.0;
    return;
  }
  if (src == 2) {
    hp = gh(m,2*d,k,j,i);
    hm = gh(m,2*d+1,k,j,i);
    Real f2 = 0.0;
    for (int a = 0; a < 3; ++a) {f2 += SQR(gh(m,2*a,k,j,i) + gh(m,2*a+1,k,j,i));}
    fb = fmin(sqrt(f2), 1.0);
    return;
  }
  if (src < 0) {   // no ray data and no model: the face stays central
    hp = -1.0;
    hm = 0.0;
    fb = 0.0;
    return;
  }
  Real f, nd;
  if (src == 3) {
    f = fmin(fabs(fx), 1.0);
    nd = (fx >= 0.0) ? 1.0 : -1.0;
  } else {
    f = M1FLev(iw(m,M1_IW_WCHI,k,j,i));
    nd = iw(m,M1_IW_N1+d,k,j,i);
  }
  hp = 0.25*(1.0 - f) + f*fmax(nd, 0.0);
  hm = -0.25*(1.0 - f) + f*fmin(nd, 0.0);
  fb = f;
}

//----------------------------------------------------------------------------------------
//! \fn M1HrJFlag
//! \brief xthinfix-1009 (implicit_blend_xthin_mode = beam_kn): 1 when the formal-solution
//! J of the cell (vet_sc, src 1) is out of equilibrium with the iterate E by more than a
//! factor 3 (a light front crossing the cell this step: the lagged E has no gradient to
//! test yet), else 0.  No J for the other sources: 0.

template <class V, class VC>
KOKKOS_INLINE_FUNCTION
Real M1HrJFlag(const int src, const V &iw, const VC &vc, const int m, const int k,
               const int j, const int i) {
  if (src != 1) {return 0.0;}
  const Real jj = vc(m,M1_VET_J,k,j,i);
  const Real ee = iw(m,M1_IW_EP,k,j,i);
  return (jj > 3.0*ee || 3.0*jj < ee) ? 1.0 : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn M1HrFace
//! \brief hrup-1009: AL, HCL, HCR of one face from hp of the LEFT cell and hm of the
//! RIGHT cell: the upwind half-range flux c (hp_L E_L + hm_R E_R) weighted by the AP
//! weight w (implicit_blend = tau: exp(-(tau_f/tau0)^2), exactly 0.0 deep in the thick
//! interior -> the face-eliminated central flux to round-off; berthon: 1).  HCL >= 0 and
//! HCR <= 0 by construction: the M-matrix column argument of sect. 7 holds unchanged.

//! Weights (implicit_blend): tau  w = exp(-(tau_f/tau0)^2);  tau_f  that times the
//! smoothstep of the beam measure |H|/J of the two cells (fmode max/mean, f_lo..f_hi);
//! idort (the user's IDORT discrete-ordinates weight)  w = 1/(1 + x + x^2/tau0),
//! x = alpha tau_f: w -> 1 thin, w -> tau0/x^2 thick, so the excess (upwind) diffusion
//! falls like 1/tau;  knudsen  w = R^4/(R^4 + R0^4) exp(-(tau_f/tau0)^2), R = |E_R -
//! E_L| / (tau_f max(E_L, E_R)) the face KNUDSEN number of the lagged E (the photon mean
//! free path over the gradient length; R = 3 f in a diffusion regime, R -> inf in
//! vacuum): central wherever the field is DIFFUSIVE however thin the cell (the grey
//! atmosphere's semi-thin cells, G1), upwind where the mean free path exceeds the
//! gradient length (beams, shadows in vacuum).  The face flux is the CONVEX blend
//! (1 - w) F_central + w F_hr for
//! every kind: w = 1 is the pure half-range upwind flux, w = 0 the face-eliminated
//! central flux, and both column-M-matrix structures survive any w in [0, 1].

KOKKOS_INLINE_FUNCTION
void M1HrFace(const Real hpl, const Real hmr, const Real fbl, const Real fbr,
              const Real tauf, const bool blend, const int bkind, const int bfm,
              const Real btau0, const Real bflo, const Real bfhi, const Real balpha,
              const Real ch, const Real el, const Real er, const Real br0,
              const Real hml, const Real hpr, const Real cdx, const Real bx0,
              const int bxmode, const Real bxwmin, const Real bxr0, const Real bxj,
              Real &alw, Real &ccl, Real &ccr) {
  if (hpl < 0.0 || hmr > 0.0) {   // no formal solution in one of the two cells
    alw = 0.0;
    ccl = 0.0;
    ccr = 0.0;
    return;
  }
  Real w = 1.0;
  if (blend) {
    if (bkind == M1_IBLEND_IDORT || bkind == M1_IBLEND_IDF || bkind == M1_IBLEND_IDA) {
      const Real x = balpha*tauf;
      w = 1.0/(1.0 + x + x*x/btau0);
      if (bkind != M1_IBLEND_IDORT) {
        // the beam-aware factor: smoothstep(flo..fhi) of the FS flux factor |H|/J
        // (idort_f)
        // or of the half-range ASYMMETRY along the face normal, |h+ + h-|/(h+ - h-) (the
        // normal flux over the normal |mu|-moment; 0 isotropic, 1 one-sided) (idort_a);
        // max or mean over the two cells (implicit_blend_fmode)
        Real gl = fbl, gr = fbr;
        if (bkind == M1_IBLEND_IDA) {
          gl = fabs(hpl + hml)/fmax(hpl - hml, 1.0e-300);
          gr = fabs(hpr + hmr)/fmax(hpr - hmr, 1.0e-300);
        }
        w *= M1BlendWeight(M1_IBLEND_F, bfm, 0.0, 1.0, gl, gr, bflo, bfhi);
      }
    } else if (bkind == M1_IBLEND_KN) {
      const Real em = fmax(fmax(el, er), 1.0e-300);
      const Real r = fabs(er - el)/fmax(tauf*em, 1.0e-300*em);
      const Real r4 = SQR(SQR(fmin(r, 1.0e30)));
      const Real x = tauf/btau0;
      w = r4/(r4 + SQR(SQR(br0)))*exp(-x*x);
    } else {
      w = M1BlendWeight(bkind, bfm, tauf, btau0, fbl, fbr, bflo, bfhi);
    }
  }
  // implicit_blend_xthin = X0 > 0: a face TRANSPARENT over the step goes upwind whatever
  // the weight, w -> 1 - (1 - w) (1 - X^2/(X^2 + X0^2)), X = (c dt/dx)/(1 + c dt chi_f)
  // the
  // ratio of the central face coefficient to the upwind one (X -> c dt/dx in vacuum, ->
  // 1/tau_cell in a thick cell).  Mixed faces with X >> 1 make the 7-point system
  // ill-conditioned (BiCGStab 200-550 its or stagnation on xb20) and the central flux has
  // no diffusion limit there anyway.
  // xthinfix-1009: where the override acts (implicit_blend_xthin_mode)
  //   all (0):     every face (hrup-1009);
  //   beam (1):    only faces whose weight is already > implicit_blend_xthin_wmin (the
  //                mixed beam/central faces); diffusive faces keep the plain weight, so
  //                smooth fields stay central and converge (order_1009);
  //   beam_kn (2): beam, plus the TRANSPARENT faces with w <= wmin whose field is NOT
  //                diffusive: weight max(R^4/(R^4 + R0^4), bxj) X^8/(X^8 + X0^8),
  //                R = |E_R - E_L|/(tau_f max(E_L, E_R)) the face Knudsen ratio (= 3|F|/
  //                cE in the diffusion limit, -> infinity in vacuum next to a beam),
  //                R0 = implicit_blend_xthin_r0, and bxj = 1 where the formal-solution J
  //                and E differ by > 3x (a light front, M1HrJFlag).  The steep X^8 keeps
  //                the override off for X < X0 (a thin time-dependent field at c dt/dx
  //                < X0 stays central and consistent).  As in the always-dissipative HLLE
  //                of Jiang+ 2012 / Menon+ 2022, transparent non-diffusive faces stay
  //                upwind (realisable, well conditioned); a diffusive transition layer
  //                stays central (xthinfix_1009).
  if (bx0 > 0.0) {
    const Real x = cdx/(1.0 + cdx*tauf);
    Real s2 = x*x/(x*x + bx0*bx0);   // the hrup-1009 expression, kept for bitwise GPU
    if (bxmode == 3) {
      // steep (xthinfix-1009): every face, independent of w, with the steep transparency
      // gate X^8/(X^8 + X0^8): ~0 for X < X0 (a field the step resolves in time stays
      // central and time-consistent), ~1 for X >> X0 (the light crosses many cells per
      // step: the quasi-static upwind limit, well conditioned and realisable)
      const Real x8 = SQR(SQR(SQR(fmin(x/bx0, 1.0e30))));
      s2 = x8/(x8 + 1.0);
    } else if (bxmode != 0 && !(w > bxwmin)) {
      Real g = 0.0;
      if (bxmode == 2) {
        const Real em = fmax(fmax(el, er), 1.0e-300);
        const Real r = fabs(er - el)/fmax(tauf*em, 1.0e-300*em);
        const Real r4 = SQR(SQR(fmin(r, 1.0e30)));
        const Real x8 = SQR(SQR(SQR(fmin(x/bx0, 1.0e30))));
        g = fmax(r4/(r4 + SQR(SQR(bxr0))), bxj)*x8/(x8 + 1.0);
      }
      s2 = g;
    }
    w = 1.0 - (1.0 - w)*(1.0 - s2);
  }
  alw = w;
  ccl = w*ch*hpl;
  ccr = w*ch*hmr;
}

//----------------------------------------------------------------------------------------
//! \fn M1EnthIdx
//! \brief implicit_enthalpy: the cell index a face stencil may read along one direction,
//! for the raw (unwrapped) index ii.  With a periodic wrap inside the block (cyc) the
//! index is wrapped; otherwise it must lie in [lo - hlo, hi + hhi], hlo/hhi being the
//! number of ghost layers filled on that side (0 at a physical boundary).

KOKKOS_INLINE_FUNCTION
int M1EnthIdx(const int ii, const int lo, const int hi, const bool cyc, const int hlo,
              const int hhi, bool &ok) {
  if (cyc) {
    const int n = hi - lo + 1;
    int r = ii;
    if (r < lo) {r += n;}
    if (r > hi) {r -= n;}
    ok = true;
    return r;
  }
  ok = (ii >= lo - hlo) && (ii <= hi + hhi);
  return ok ? ii : lo;
}

//----------------------------------------------------------------------------------------
//! \fn M1EnthCorr
//! \brief implicit_enthalpy: (high-order face enthalpy flux) - (donor-cell face flux the
//! matrix carries), both at the lagged iterate.  el2, el, er, er2 are E and al2, al, ar,
//! ar2 the advective coefficients a at the cells L-1, L, R, R+1 of the face, and ok says
//! whether L-1 AND R+1 exist; vf is the face velocity whose sign picks the donor cell of
//! the matrix.  Requiring BOTH outer cells for the plm form makes the choice the same for
//! the two MeshBlocks that share a face (each of them misses one of the two when the halo
//! is too thin), so the face flux stays single valued.
//!   central: a_f = (a_L + a_R)/2, E_f = (E_L + E_R)/2.
//!   plm:     a_f = the mean of the two van Leer plm face values of a (it reduces to the
//!            central mean at an extremum of a and is a 4-point interpolation where a is
//!            smooth and monotone), E_f = the van Leer plm value of E from the side
//!            upwind of a_f (monotone, E_f <= 2 E_donor).  Central where !ok.

KOKKOS_INLINE_FUNCTION
Real M1EnthCorr(const int mode, const Real el2, const Real el, const Real er,
                const Real er2, const bool ok, const Real al2, const Real al,
                const Real ar, const Real ar2, const Real vf) {
  const Real alow = (vf > 0.0) ? (al*el) : (ar*er);
  Real af = 0.5*(al + ar);
  Real ef = 0.5*(el + er);
  if (mode == M1_IENTH_PLM && ok) {
    Real dum, afl, afr;
    PLM(al2, al, ar, afl, dum);
    PLM(al, ar, ar2, dum, afr);
    af = 0.5*(afl + afr);
    if (af > 0.0) {
      PLM(el2, el, er, ef, dum);
    } else if (af < 0.0) {
      PLM(el, er, er2, dum, ef);
    }
  }
  return af*ef - alow;
}

//----------------------------------------------------------------------------------------
//! \fn M1EnthEf
//! \brief implicit_vimp: the face value of E that the enthalpy flux of implicit_enthalpy
//! = mode carries at the lagged iterate (the donor cell for upwind, else the E_f of
//! M1EnthCorr), which multiplies the implicit velocity increment of the face.

KOKKOS_INLINE_FUNCTION
Real M1EnthEf(const int mode, const Real el2, const Real el, const Real er,
              const Real er2, const bool ok, const Real al2, const Real al,
              const Real ar, const Real ar2, const Real vf) {
  if (mode == M1_IENTH_UPWIND) {return (vf > 0.0) ? el : er;}
  Real ef = 0.5*(el + er);
  if (mode == M1_IENTH_PLM && ok) {
    Real dum, afl, afr;
    PLM(al2, al, ar, afl, dum);
    PLM(al, ar, ar2, dum, afr);
    Real af = 0.5*(afl + afr);
    if (af > 0.0) {
      PLM(el2, el, er, ef, dum);
    } else if (af < 0.0) {
      PLM(el, er, er2, dum, ef);
    }
  }
  return ef;
}

//----------------------------------------------------------------------------------------
//! \fn M1EnthCorrT, M1EnthEfT
//! \brief time_scheme = hesdirk2 stage solves (time2_enth_vel = start): M1EnthCorr and
//! M1EnthEf with the plm face value of a reconstructed from the STAGE-START a = a - d
//! (d = DA, the part of a carried by the old vector's FSAL velocity increment), and d
//! added back as its face MEAN.  With d = 0 they are M1EnthCorr / M1EnthEf.

KOKKOS_INLINE_FUNCTION
Real M1EnthAfT(const int mode, const bool ok, const Real al2, const Real al,
               const Real ar, const Real ar2, const Real dl2, const Real dl,
               const Real dr, const Real dr2, const bool afc) {
  if (afc) {return 0.5*(al + ar);}   // time2_enth_vel = central: a_f central
  const Real sl2 = al2 - dl2, sl = al - dl, sr = ar - dr, sr2 = ar2 - dr2;
  Real af = 0.5*(sl + sr);
  if (mode == M1_IENTH_PLM && ok) {
    Real dum, afl, afr;
    PLM(sl2, sl, sr, afl, dum);
    PLM(sl, sr, sr2, dum, afr);
    af = 0.5*(afl + afr);
  }
  return af + 0.5*(dl + dr);
}

KOKKOS_INLINE_FUNCTION
Real M1EnthCorrT(const int mode, const Real el2, const Real el, const Real er,
                 const Real er2, const bool ok, const Real al2, const Real al,
                 const Real ar, const Real ar2, const Real dl2, const Real dl,
                 const Real dr, const Real dr2, const Real vf, const bool afc) {
  const Real alow = (vf > 0.0) ? (al*el) : (ar*er);
  const Real af = M1EnthAfT(mode, ok, al2, al, ar, ar2, dl2, dl, dr, dr2, afc);
  Real ef = 0.5*(el + er);
  if (mode == M1_IENTH_PLM && ok) {
    Real dum;
    if (af > 0.0) {
      PLM(el2, el, er, ef, dum);
    } else if (af < 0.0) {
      PLM(el, er, er2, dum, ef);
    }
  }
  return af*ef - alow;
}

KOKKOS_INLINE_FUNCTION
Real M1EnthEfT(const int mode, const Real el2, const Real el, const Real er,
               const Real er2, const bool ok, const Real al2, const Real al,
               const Real ar, const Real ar2, const Real dl2, const Real dl,
               const Real dr, const Real dr2, const Real vf, const bool afc) {
  if (mode == M1_IENTH_UPWIND) {return (vf > 0.0) ? el : er;}
  Real ef = 0.5*(el + er);
  if (mode == M1_IENTH_PLM && ok) {
    const Real af = M1EnthAfT(mode, ok, al2, al, ar, ar2, dl2, dl, dr, dr2, afc);
    Real dum;
    if (af > 0.0) {
      PLM(el2, el, er, ef, dum);
    } else if (af < 0.0) {
      PLM(el, er, er2, dum, ef);
    }
  }
  return ef;
}

//----------------------------------------------------------------------------------------
//! \fn M1VimpRow
//! \brief implicit_vimp: the part of the operator row outside the 7-point row (x1 +-2,
//! x2 and x3 -2..+2 without 0), applied to the component cx.  b = RadiationM1::iw_vimp.

KOKKOS_INLINE_FUNCTION
Real M1VimpRow(const DvceArray5D<Real> &iw, const int b, const int cx, const int m,
               const int k, const int j, const int i, const int is, const int ie,
               const bool cyclic, const bool thrd) {
  int im2 = i - 2, ip2 = i + 2;
  if (cyclic) {
    const int n = ie - is + 1;
    while (im2 < is) {im2 += n;}
    while (ip2 > ie) {ip2 -= n;}
  }
  Real y = iw(m,b+M1_IV_X1M2,k,j,i)*iw(m,cx,k,j,im2)
           + iw(m,b+M1_IV_X1P2,k,j,i)*iw(m,cx,k,j,ip2)
           + iw(m,b+M1_IV_X2M2,k,j,i)*iw(m,cx,k,j-2,i)
           + iw(m,b+M1_IV_X2M2+1,k,j,i)*iw(m,cx,k,j-1,i)
           + iw(m,b+M1_IV_X2M2+2,k,j,i)*iw(m,cx,k,j+1,i)
           + iw(m,b+M1_IV_X2M2+3,k,j,i)*iw(m,cx,k,j+2,i);
  if (thrd) {
    y += iw(m,b+M1_IV_X3M2,k,j,i)*iw(m,cx,k-2,j,i)
         + iw(m,b+M1_IV_X3M2+1,k,j,i)*iw(m,cx,k-1,j,i)
         + iw(m,b+M1_IV_X3M2+2,k,j,i)*iw(m,cx,k+1,j,i)
         + iw(m,b+M1_IV_X3M2+3,k,j,i)*iw(m,cx,k+2,j,i);
  }
  return y;
}

//----------------------------------------------------------------------------------------
//! \fn M1MusclRow
//! \brief implicit_hr_recon = plm (xthinfix-1009 Fix B): the plm part of the half-range
//! face fluxes in the row of cell (k,j,i), applied to the component cx: diag, x1 -2..+2,
//! x2 -2..+2, x3 -2..+2 (no 0).  b = RadiationM1::iw_muscl (ImplicitMusclBuild).

KOKKOS_INLINE_FUNCTION
Real M1MusclRow(const DvceArray5D<Real> &iw, const int b, const int cx, const int m,
                const int k, const int j, const int i, const int is, const int ie,
                const bool cyclic, const bool thrd) {
  int im2 = i - 2, im1 = i - 1, ip1 = i + 1, ip2 = i + 2;
  if (cyclic) {
    const int n = ie - is + 1;
    while (im2 < is) {im2 += n;}
    while (im1 < is) {im1 += n;}
    while (ip1 > ie) {ip1 -= n;}
    while (ip2 > ie) {ip2 -= n;}
  }
  Real y = iw(m,b+M1_IM_D,k,j,i)*iw(m,cx,k,j,i)
           + iw(m,b+M1_IM_X1,k,j,i)*iw(m,cx,k,j,im2)
           + iw(m,b+M1_IM_X1+1,k,j,i)*iw(m,cx,k,j,im1)
           + iw(m,b+M1_IM_X1+2,k,j,i)*iw(m,cx,k,j,ip1)
           + iw(m,b+M1_IM_X1+3,k,j,i)*iw(m,cx,k,j,ip2)
           + iw(m,b+M1_IM_X2,k,j,i)*iw(m,cx,k,j-2,i)
           + iw(m,b+M1_IM_X2+1,k,j,i)*iw(m,cx,k,j-1,i)
           + iw(m,b+M1_IM_X2+2,k,j,i)*iw(m,cx,k,j+1,i)
           + iw(m,b+M1_IM_X2+3,k,j,i)*iw(m,cx,k,j+2,i);
  if (thrd) {
    y += iw(m,b+M1_IM_X3,k,j,i)*iw(m,cx,k-2,j,i)
         + iw(m,b+M1_IM_X3+1,k,j,i)*iw(m,cx,k-1,j,i)
         + iw(m,b+M1_IM_X3+2,k,j,i)*iw(m,cx,k+1,j,i)
         + iw(m,b+M1_IM_X3+3,k,j,i)*iw(m,cx,k+2,j,i);
  }
  return y;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn ImplBCFromString
//! \brief one x1 boundary type of the implicit solve, from the input string, or from the
//! mesh boundary flag when the input does not name one.

int ImplBCFromString(const std::string &s, const BoundaryFlag mbc) {
  if (s.compare("auto") == 0) {
    if (mbc == BoundaryFlag::periodic) return M1_IBC_PERIODIC;
    if (mbc == BoundaryFlag::reflect) return M1_IBC_REFLECT;
    return M1_IBC_MARSHAK;
  }
  if (s.compare("marshak") == 0) return M1_IBC_MARSHAK;
  if (s.compare("flux") == 0) return M1_IBC_FLUX;
  if (s.compare("reflect") == 0) return M1_IBC_REFLECT;
  if (s.compare("periodic") == 0) return M1_IBC_PERIODIC;
  if (s.compare("efix") == 0) return M1_IBC_EFIX;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    << std::endl << "<rad_m1>/implicit_bc_x1min|max = '" << s << "' is not a valid "
    << "choice (auto | marshak | flux | reflect | periodic | efix)" << std::endl;
  std::exit(EXIT_FAILURE);
  return M1_IBC_MARSHAK;
}

//----------------------------------------------------------------------------------------
//! \fn ImplFatal
//! \brief one fatal-error exit with a message, used by the 3a2 option parsers

void ImplFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitInit
//! \brief read the <rad_m1> parameters of the implicit solver, check the restrictions of
//! 3a and allocate the face array and the work array.  Called at the END of the
//! constructor, from the delimited hook there, and a no-op in explicit mode.

void RadiationM1::ImplicitInit(ParameterInput *pin) {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;

  // ---- MILESTONE 3b phase B.  `implicit` = the full 7-point solve; `implicit_x1` keeps
  // every branch below on the 3a/3a2/3c arithmetic, bit for bit.
  const bool full = (transport == M1_TRANSPORT_IMPLICIT);
  // <rad_m1>/time_scheme (read by Time2Init below).  DEFAULT hesdirk2 since m1-defaults2
  // (tests_m1/runs_5g_defaults2; /viper/ptmp2/jinma/h2val_0924: hesdirk2 at cfl 0.9 is
  // 1.57x cheaper than be at cfl 0.3 at equal accuracy) wherever hesdirk2 is accepted:
  // transport = implicit (not implicit_x1) on a Cartesian mesh or (m1-sph2,
  // tests_m1/runs_5h_sph2) the spherical-polar wedge (no cubed sphere or polar
  // boundary), a closure other than tau (vet_col included since m1-sph2), and
  // <time>/integrator = rk2.  Elsewhere the default stays be, silently.  A restart whose
  // file lacks the key (written when it was read only when named) keeps be.  The
  // resolved value is echoed (so later restarts keep it); explicit input overrides.
  {auto *pmh = pmy_pack->pmesh;
  const std::string integ = pin->DoesParameterExist("time","integrator") ?
                            pin->GetString("time","integrator") : "rk2";
  const bool h2def = full && !pmh->use_cubed_sphere && !pmh->use_polar_boundary &&
                     (!tau_closure || vet_col) && (integ.compare("rk2") == 0) &&
                     !global_variable::restart_run;
  (void) pin->GetOrAddString("rad_m1","time_scheme", h2def ? "hesdirk2" : "be");}
  trans_on = false;
  trans_x3 = false;
  impl_cfl = pin->GetOrAddReal("rad_m1","implicit_cfl",-1.0);
  impl_tol = pin->GetOrAddReal("rad_m1","implicit_tol",1.0e-8);
  // a lagged transverse coupling needs more outer passes than a pure column solve, so
  // the ceiling is raised (not the count: the loop still stops at convergence).
  impl_maxit = pin->GetOrAddInteger("rad_m1","implicit_maxit", full ? 200 : 30);
  impl_lin_tol = pin->GetOrAddReal("rad_m1","implicit_lin_tol",1.0e-10);
  impl_lin_maxit = pin->GetOrAddInteger("rad_m1","implicit_lin_maxit",200);
  {std::string sv = pin->GetOrAddString("rad_m1","implicit_solver","line_jacobi");
  if (sv.compare("line_jacobi") == 0) {
    impl_solver = M1_ISOLV_LINE_JACOBI;
  } else if (sv.compare("bicgstab") == 0) {
    impl_solver = M1_ISOLV_BICGSTAB;
  } else {
    ImplFatal("<rad_m1>/implicit_solver = '" + sv
              + "' is not a choice (line_jacobi | bicgstab)");
  }
  }
  // ---- MILESTONE 3b phase D: what happens to the off-diagonal Eddington terms, and how
  // fast the closure is allowed to move between Picard passes.  The defaults reproduce
  // phase C on every 1-D / implicit_x1 configuration (there are no off-diagonal terms and
  // no transverse closure there), and `operator` is the default in multi-D, where lagging
  // them has no fixed point in an optically thin cell.
  std::string sod = pin->GetOrAddString("rad_m1","implicit_offdiag","auto");
  bool od_auto = (sod.compare("auto") == 0);
  if (od_auto || sod.compare("lagged") == 0) {
    impl_offdiag = M1_OD_LAGGED;
  } else if (sod.compare("operator") == 0) {
    impl_offdiag = M1_OD_OPERATOR;
  } else if (sod.compare("none") == 0) {
    impl_offdiag = M1_OD_NONE;
  } else {
    ImplFatal("<rad_m1>/implicit_offdiag = '" + sod
              + "' is not a choice (auto | lagged | operator | none)");
  }
  impl_crelax = pin->GetOrAddReal("rad_m1","implicit_closure_relax",1.0);
  if (!(impl_crelax > 0.0) || impl_crelax > 1.0) {
    ImplFatal("<rad_m1>/implicit_closure_relax must lie in (0,1]");
  }
  impl_crelax_thin = pin->GetOrAddBoolean("rad_m1","implicit_closure_relax_thin",false);
  // implicit_closure_lag: DEFAULT step since defaults-0927 (2026-09-27) for transport =
  // implicit on a multi-D mesh (runs_3b5 recommends step in 2-D/3-D, every modern input
  // names it, pass is ~4x slower per simulated second; switch_inventory_2026-09-24.md);
  // pass elsewhere and on a restart whose file lacks the key.  The value is recorded.
  const bool md_full = full && pmy_pack->pmesh->multi_d;
  {std::string sc = pin->GetOrAddString("rad_m1","implicit_closure_lag",
      (md_full && !global_variable::restart_run) ? "step" : "pass");
  if (sc.compare("pass") == 0) {
    impl_clag_step = false;
  } else if (sc.compare("step") == 0) {
    impl_clag_step = true;
  } else {
    ImplFatal("<rad_m1>/implicit_closure_lag = '" + sc
              + "' is not a choice (pass | step)");
  }
  }
  // runs_5c_thinstab: the step-to-step closure relaxation of optically thin cells.
  // DEFAULT 1.5 since defaults-0927 (2026-09-27; memory m1-thin-cell-cure, merged
  // b7e48012: the cure of the multi-D thin-cell instability of the chi(f) closures, He
  // slab 332 fallbacks -> 0 at 0 extra passes, be and hesdirk2) for a chi(f) closure (m1,
  // minerbo, kershaw) under transport = implicit on a multi-D mesh with
  // implicit_closure_lag = step -- its own requirements; the fixed-tensor closures
  // (eddington, vet_sc, tau, vet_col) do not need it.  0 (off) elsewhere and on a
  // restart whose file lacks the key; the value is recorded where the default applies.
  impl_ctrelax = 0.0;
  ctr_init = false;
  {const bool chif = !(eddington || vet_sc || tau_closure);
  if (chif && md_full && impl_clag_step) {
    impl_ctrelax = pin->GetOrAddReal("rad_m1","implicit_closure_thin_relax",
                                     global_variable::restart_run ? 0.0 : 1.5);
  } else if (pin->DoesParameterExist("rad_m1","implicit_closure_thin_relax")) {
    impl_ctrelax = pin->GetReal("rad_m1","implicit_closure_thin_relax");
    // accel-1009 (memory): a fixed-tensor closure (eddington, vet_sc, tau, vet_col) never
    // relaxes (impl_lag_kernel: ctr && !edd && !vetsc && !tkeep && !tauc), so the 4-slab
    // closure memory stayed all-zero and was still written to every restart: off
    if (!chif) {impl_ctrelax = 0.0;}
  }
  if (impl_ctrelax < 0.0) {
    ImplFatal("<rad_m1>/implicit_closure_thin_relax must be >= 0");
  }
  }
  // ---- the TRANSVERSE realizability limiter.  `none` is not merely the default: it
  // takes the OLD expression for theta everywhere (see ImplicitTransTheta), so an input
  // file that does not name it is bitwise unchanged.
  {std::string st = pin->GetOrAddString("rad_m1","implicit_trans_limit","none");
  if (st.compare("none") == 0) {
    impl_tlim = M1_TLIM_NONE;
  } else if (st.compare("lp") == 0) {
    impl_tlim = M1_TLIM_LP;
  } else {
    ImplFatal("<rad_m1>/implicit_trans_limit = '" + st
              + "' is not a choice (none | lp)");
  }
  }
  impl_tfmax = pin->GetOrAddReal("rad_m1","implicit_trans_fmax",1.0);
  if (!(impl_tfmax > 0.0)) {
    ImplFatal("<rad_m1>/implicit_trans_fmax must be positive");
  }
  // ---- MILESTONE 3e: ANDERSON acceleration of the outer (Picard) iteration.  `none` is
  // not merely the default: nothing below is allocated and ImplicitAccelSave/Apply are
  // never called, so an input file that does not name it is bitwise unchanged.
  {std::string sa = pin->GetOrAddString("rad_m1","implicit_accel","none");
  if (sa.compare("none") == 0) {
    impl_accel = M1_IACC_NONE;
  } else if (sa.compare("anderson") == 0) {
    impl_accel = M1_IACC_ANDERSON;
  } else {
    ImplFatal("<rad_m1>/implicit_accel = '" + sa
              + "' is not a choice (none | anderson)");
  }
  }
  // ---- MILESTONE 3g: the GAS-RADIATION energy coupling.  implicit_eos_cache defaults
  // false; neither key allocates anything nor is referenced when false.
  // implicit_gas_newton: DEFAULT true since defaults-0927 (2026-09-27;
  // runs_3g_newton_T VERDICT PASS, and every GPU timing and gate since 3g ran it;
  // switch_inventory_2026-09-24.md decision 8) for transport = implicit; false under
  // implicit_x1 (the 3a arithmetic) and on a restart whose file lacks the key.  The
  // value is recorded.
  impl_gas_newton = pin->GetOrAddBoolean("rad_m1","implicit_gas_newton",
                                         full && !global_variable::restart_run);
  impl_eos_cache = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache",false);
  pin_report_newton_fb = pin->GetOrAddBoolean("rad_m1","report_newton_fb",false);
  dbg_etally = pin->GetOrAddBoolean("rad_m1","dbg_energy_tally",false);
  impl_ecnt = pin->GetOrAddInteger("rad_m1","implicit_eos_cache_nt",2);
  impl_eccheck = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache_check",true);
  // the check is a MEASUREMENT only (nothing reads igm or ec_emax/ec_tmax but the final
  // report): implicit_eos_cache_check_every = N takes it on the cycles N divides.  Read
  // only when named; the default 1 is the old every-step check.
  impl_eccheck_every = 1;
  if (pin->DoesParameterExist("rad_m1","implicit_eos_cache_check_every")) {
    impl_eccheck_every = pin->GetInteger("rad_m1","implicit_eos_cache_check_every");
    if (impl_eccheck_every < 1) {
      ImplFatal("<rad_m1>/implicit_eos_cache_check_every must be >= 1");
    }
  }
  // the x1 LINE SOLVE (preconditioner and line-Jacobi pass): thomas = one thread per
  // column, serial recurrence (the original); pcr = one team per column, parallel
  // cyclic reduction in team scratch (GPU).  Same system; the answers agree to
  // round-off, not bitwise.  The gathered stack sweep (part_nblk > 1) is Thomas always.
  // DEFAULT pcr since bench/m1_defaults_0923 (15x faster line solve on the GPU, +1 %
  // wall on the host; thomas reproduces the earlier default bitwise).
  {std::string ls = pin->GetOrAddString("rad_m1","implicit_line_solver","pcr");
  if (ls.compare("thomas") == 0) {
    impl_line_solver = 0;
  } else if (ls.compare("pcr") == 0) {
    impl_line_solver = 1;
  } else {
    ImplFatal("<rad_m1>/implicit_line_solver = '" + ls
              + "' is not a choice (thomas | pcr)");
  }
  }
  impl_pcr_team = pin->GetOrAddInteger("rad_m1","implicit_pcr_team",0);
  impl_pcr_check = pin->GetOrAddBoolean("rad_m1","implicit_pcr_check",false);
  if (impl_pcr_team < 0 || impl_pcr_team > 1024) {
    ImplFatal("<rad_m1>/implicit_pcr_team must lie in [0,1024]");
  }
  // SPEED-UP 3: the host synchronisations of the BiCGStab loop.  0 = the original loop;
  // 1 = the same recurrence with fused reductions and vector updates (3 blocking
  // reductions per iteration instead of 5); 2 = 1 with alpha kept on the device, so
  // rhat.v does not block either (1 rank only; with more ranks 2 acts as 1).  Levels
  // 1 and 2 sum in a different order than 0: same answer to round-off, not bitwise.
  // DEFAULT 1 since bench/m1_defaults_0923 (0 reproduces the earlier default).
  impl_bcg_sync = pin->GetOrAddInteger("rad_m1","implicit_bcg_sync",1);
  if (impl_bcg_sync < 0 || impl_bcg_sync > 2) {
    ImplFatal("<rad_m1>/implicit_bcg_sync must be 0, 1 or 2");
  }
  if (impl_bcg_sync == 2) {
    bcg_rvd = Kokkos::View<Real, DevMemSpace>("m1_bcg_rvd");
  }
  // GPU COST OF THE KRYLOV ITERATION (bench/m1_fast_0923).  Off, nothing is allocated
  // or called and the path is bitwise the pre-0923 code.  Since
  // tests_m1/runs_3p_fastdefault they DEFAULT ON (the `FAST` set of
  // tests_m1/runs_3k_gpu3d/README_FAST.md: halo_direct, od_cache, krylov_fuse = 3,
  // op_stencil, precond = rbgs_fwd; 2.0-2.4x on the GPU, the converged answer moved at
  // the +-1-ulp level) for the closures whose tensor is fixed within a step (eddington,
  // vet_sc, tau) whenever the configuration admits them: transport = implicit with
  // bicgstab, implicit_bcg_sync = 1, implicit_line_solver = pcr, ONE MeshBlock along x1,
  // and (op_stencil) no periodic x1 wrap; with implicit_lin_cnorm > 0 krylov_fuse is 1.
  // m1 / minerbo / kershaw keep them off (not gated there).  A key the input names keeps
  // its value.
  //  implicit_halo_direct = true   the implicit exchanges as ONE on-rank copy kernel
  //    when every neighbour of every block is on its own rank at the same level
  //    (otherwise the ordinary exchange, as before).  Bitwise.
  //  implicit_od_cache = true      the off-diagonal Eddington operator from a per-cell
  //    cache of sum_e d_e P_de, fused with the 7-point row into one kernel.  Bitwise.
  //  implicit_krylov_fuse = 1      the pcr preconditioner reads and writes the Krylov
  //    vectors itself and carries the p / s updates (bitwise); = 2 also puts the
  //    (rhat,v) and (t,s),(t,t) reductions in the operator kernel (round-off: the sums
  //    are ordered differently).  Needs implicit_bcg_sync = 1, implicit_line_solver =
  //    pcr, one block along x1; 2 needs implicit_od_cache.  3 = 2 with TWO blocking
  //    reductions per iteration (ImplicitBiCGStabTwo; round-off).
  //  implicit_op_stencil = true    the frozen operator of each Picard pass written out
  //    once as a 19-point stencil (7-point row + off-diagonal Eddington terms), applied
  //    by one kernel per Krylov product (round-off: the terms are grouped differently).
  //    Needs implicit_bcg_sync = 1; not with a periodic x1 wrap.
  //  implicit_precond_float = true  the line solves of the fused path's preconditioner
  //    in float (the operator, the vectors and every reduction stay double: a
  //    preconditioner only sets the convergence rate).  Needs implicit_krylov_fuse >= 1.
  //  implicit_precond = line | rbgs | rbgs_fwd   the preconditioner of the fused path:
  //    x1 line Jacobi (the original), symmetric / forward red-black transverse line
  //    Gauss-Seidel, block-local (ImplicitPrecondX).  Needs implicit_krylov_fuse >= 1.
  const bool fixcl = eddington || vet_sc || tau_closure;
  bool fdef = fixcl && full && (impl_solver == M1_ISOLV_BICGSTAB) &&
              (impl_bcg_sync == 1) && (impl_line_solver == 1) &&
              (pmy_pack->pmesh->mesh_indcs.nx1 == pmy_pack->pmesh->mb_indcs.nx1);
  // implicit_lin_cnorm > 0 (m1-cnorm-fastpath, bsg_1001/CNORM_FAST_1003.md): the fast
  // set stays on, with implicit_krylov_fuse = 1 instead of 3.  The per-cell norm is
  // taken in the reductions of ImplicitBiCGStab (fuse 0/1/2), which ImplicitBiCGStabTwo
  // (fuse 3) does not have; fuse 1 is the p-update in the preconditioner (bitwise vs 0)
  // and costs the same as 3 (He wedge 0.392-0.422 vs 0.410-0.422 s/cycle).  Before, any
  // cnorm > 0 dropped the whole fast set (fdef = false: 2.2x slower on the He wedge).
  const bool cn_def = pin->DoesParameterExist("rad_m1","implicit_lin_cnorm") &&
                      pin->GetReal("rad_m1","implicit_lin_cnorm") > 0.0;
  // the x1 wrap as the boundary parsing below will read it (same keys, same defaults)
  const bool x1per = (ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      pmy_pack->pmesh->mesh_bcs[static_cast<int>(BoundaryFace::inner_x1)])
      == M1_IBC_PERIODIC);
  impl_halo_direct = pin->GetOrAddBoolean("rad_m1","implicit_halo_direct",fdef);
  impl_odc = pin->GetOrAddBoolean("rad_m1","implicit_od_cache",fdef);
  impl_kfuse = pin->GetOrAddInteger("rad_m1","implicit_krylov_fuse",
                                    fdef ? (cn_def ? 1 : 3) : 0);
  // a restart file written without cnorm echoes the defaulted fuse 3: with cnorm named
  // now it runs fuse 1 (round-off vs 3), echoed; a fuse 3 the input names still FATALs
  if (impl_kfuse == 3 && cn_def &&
      pin->IsParameterDefaulted("rad_m1","implicit_krylov_fuse")) {
    impl_kfuse = 1;
    pin->SetInteger("rad_m1","implicit_krylov_fuse",1);
    std::cout << "rad_m1: implicit_krylov_fuse 3 (defaulted in the restart file) -> 1, "
              << "since implicit_lin_cnorm > 0" << std::endl;
  }
  impl_stencil = pin->GetOrAddBoolean("rad_m1","implicit_op_stencil",fdef && !x1per);
  impl_prec_float = pin->GetOrAddBoolean("rad_m1","implicit_precond_float",false);
  // tests_m1/runs_4a_accel levers (implicit_vimp_fold, implicit_one_pass,
  // implicit_predictor_order, implicit_fast_kernels): resolved further down, once
  // implicit_predictor and implicit_vimp are known (their defaults depend on them)
  impl_vfold = false;
  impl_onep = 0;
  impl_onep_s = 3.0;
  for (int t = 0; t < 3; ++t) {
    onep_qa[t] = -1.0;
    onep_qb[t] = -1.0;
    onep_cnt[t] = 0.0;
    onep_off[t] = 0.0;
    onep_actr[t] = 0.0;
    onep_prb[t] = 0.0;
  }
  impl_onep_n = 0.0;
  impl_onep_nchk = 0.0;
  impl_onep_auto = false;
  impl_onep_awin = 2;
  impl_onep_arep = 64;
  impl_onep_ndis = 0.0;
  impl_onep_nprb = 0.0;
  impl_onep_nren = 0.0;
  impl_onep_nfpr = 0.0;
  ew_tight = false;
  impl_pord = 1;
  impl_opsplit = false;
  if (pin->DoesParameterExist("rad_m1","implicit_op_split_red")) {
    impl_opsplit = pin->GetBoolean("rad_m1","implicit_op_split_red");
  }
  // default OFF since audit-m1-0928 (was on since 2026-09-25, m1-fast3: -12 % on the
  // static T-S4 wedge).  On moving gas on MI300A it gains 2.2-2.6 % (1 GPU) and 0.2-2.2 %
  // (2 GPUs) per simulated second (He box, sph_wedge; 5 interleaved reps), under the ~3 %
  // bar for a result-changing default (round-off).  Opt-in by name.
  impl_opteam = pin->GetOrAddBoolean("rad_m1","implicit_op_team_red",false);
  impl_fastk = false;
  impl_odskip = false;
  if (impl_stencil && impl_bcg_sync != 1) {
    ImplFatal("<rad_m1>/implicit_op_stencil needs implicit_bcg_sync = 1");
  }
  // Spherical-polar fast-path default: rbgs_fwd.  m1-sp-order2b had made mg_gc the sp
  // default (runs_5p_coarse2: radiation-only T-S4 wedge, -18 %); on the rad-hydro
  // sph_wedge with MOVING gas rbgs_fwd is -8.7 % against mg_gc (m1-wedge,
  // docs/dev/m1_wedge_0926.md 7.2; user decision 09-26), so the sp default is rbgs_fwd
  // again.  mg_gc stays available by name; the resolved value is echoed.
  // defaults-0927 (2026-09-27): on the Cartesian fast path the default is mg (with
  // implicit_mg_levels = 3; tests_m1/runs_5m_precond: He box -7 %/cycle on 1 GPU, -6 % on
  // 2; H200 2 GPUs: plain mg fastest, handover 09-26 sect. 6) on a multi-D, single-level,
  // non-cs/sp/polar mesh, unless implicit_krylov_dev > 0 is named (it refuses mg).  A
  // restart whose file lacks the key keeps rbgs_fwd.
  auto *pmq = pmy_pack->pmesh;
  const bool kdev_named = pin->DoesParameterExist("rad_m1","implicit_krylov_dev") &&
                          pin->GetInteger("rad_m1","implicit_krylov_dev") > 0;
  // defaults-1009 (user 10-09): mg also on the spherical-polar wedge (the AG Car A and He
  // giant productions name it there); the cubed sphere, polar boundaries, SMR/AMR, 1-D,
  // krylov_dev and the non-fast path keep their old default
  const bool mgdef = !pmq->use_cubed_sphere && !pmq->use_polar_boundary &&
                     !pmq->multilevel && pmq->multi_d && !kdev_named &&
                     !global_variable::restart_run;
  const char *pcdef = fdef ? (mgdef ? "mg" : "rbgs_fwd") : "line";
  {std::string pc = pin->GetOrAddString("rad_m1","implicit_precond", pcdef);
  if (pc.compare("line") == 0) {
    impl_prec = 0;
  } else if (pc.compare("rbgs") == 0) {
    impl_prec = 1;
  } else if (pc.compare("rbgs_fwd") == 0) {
    impl_prec = 2;
  } else if (pc.compare("mg") == 0) {
    impl_prec = 3;
  } else if (pc.compare("mg_gc") == 0) {
    impl_prec = 4;
  } else if (pc.compare("mg_gf") == 0) {
    impl_prec = 5;
  } else {
    ImplFatal("<rad_m1>/implicit_precond = '" + pc
              + "' is not a choice (line | rbgs | rbgs_fwd | mg | mg_gc | mg_gf)");
  }
  }
  // implicit_precond = mg (rad_m1_precond.cpp, tests_m1/runs_5m_precond): keys read only
  // under mg, so every other configuration is untouched
  mg_nlev = 0;
  mg_halo = true;
  if (impl_prec == 3) {
    // default 3 since defaults-0927 (the recommended value, runs_5m_precond; was 2)
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",3);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 2) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 2");}
  }
  // implicit_precond = mg_gc (rad_m1_precond.cpp, tests_m1/runs_5p_coarse2): mg with
  // implicit_mg_levels >= 1 (1 = rbgs_fwd + the global coarse space only) plus the
  // global band coarse space; its keys are read only under mg_gc
  if (impl_prec == 4) {
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",1);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 1) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 1 (mg_gc)");}
    gc_b2 = pin->GetOrAddInteger("rad_m1","implicit_gc_bands2",1);
    gc_b3 = pin->GetOrAddInteger("rad_m1","implicit_gc_bands3",1);
    gc_on = true;
  }
  // implicit_precond = mg_gf (rad_m1_precond.cpp, tests_m1/runs_5t_fast5box): mg plus
  // the global Fourier coarse space; its keys are read only under mg_gf
  if (impl_prec == 5) {
    mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",3);
    mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
    if (mg_nlev < 1) {ImplFatal("<rad_m1>/implicit_mg_levels must be >= 1 (mg_gf)");}
    gf_k = pin->GetOrAddInteger("rad_m1","implicit_gf_modes",2);
    if (gf_k < 0 || gf_k > 3) {ImplFatal("<rad_m1>/implicit_gf_modes must be 0..3");}
    gf_on = true;
  }
  // implicit_bcg_rho_direct (runs_5p_coarse2; read only when named, default false =
  // the recurrence, bitwise): implicit_krylov_fuse = 3 sums (rhat, r) directly
  impl_rho_direct = false;
  if (pin->DoesParameterExist("rad_m1","implicit_bcg_rho_direct")) {
    impl_rho_direct = pin->GetBoolean("rad_m1","implicit_bcg_rho_direct");
  }
  if (impl_kfuse < 0 || impl_kfuse > 3) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse must be 0, 1, 2 or 3");
  }
  if (impl_kfuse == 3 && pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0) > 0.0) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse = 3 does not take implicit_lin_cnorm");
  }
  if (impl_kfuse > 0 && (impl_bcg_sync != 1 || impl_line_solver != 1)) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse needs implicit_bcg_sync = 1 and "
              "implicit_line_solver = pcr");
  }
  if (impl_kfuse >= 2 && !impl_odc) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse >= 2 needs implicit_od_cache = true");
  }
  if (impl_prec_float && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond_float needs implicit_krylov_fuse >= 1");
  }
  if (impl_prec > 0 && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond = rbgs needs implicit_krylov_fuse >= 1");
  }
  // multi-rank Krylov (rad_m1_krylov.cpp, tests_m1/runs_3w_krylov).  krylov_pipe
  // defaults OFF.  halo_mpi (bitwise the ordinary exchange) defaults ON since
  // m1-defaults wherever its preconditions hold: implicit_halo_direct, a same-level mesh
  // (no SMR/AMR), no cubed-sphere seams and no polar boundary; otherwise the ordinary
  // exchange, silently.  A key the input (or restart echo) names keeps its value.
  // m1-sphhalo (tests_m1/runs_5i_sphhalo): also on the spherical-polar WEDGE (no pole:
  // every ghost is the plain copy, no vector slot flips), except on a restart whose file
  // lacks the key (written before m1-sphhalo), which keeps the ordinary exchange.
  impl_kpipe = pin->GetOrAddBoolean("rad_m1","implicit_krylov_pipe",false);
  {auto *pmh = pmy_pack->pmesh;
  const bool hmdef = impl_halo_direct && !pmh->multilevel && !pmh->use_cubed_sphere &&
                     !pmh->use_polar_boundary &&
                     !(pmh->use_spherical_polar && global_variable::restart_run);
  impl_halo_mpi = pin->GetOrAddBoolean("rad_m1","implicit_halo_mpi",hmdef);}
  hm_state = 0;
  hm_comm = nullptr;
  // implicit_halo_ipc (m1-perf-0928, rad_m1_krylov.cpp): default false
  impl_halo_ipc = pin->GetOrAddBoolean("rad_m1","implicit_halo_ipc",false);
  if (impl_halo_ipc && !impl_halo_mpi) {
    ImplFatal("<rad_m1>/implicit_halo_ipc needs implicit_halo_mpi = true");
  }
  // implicit_halo_overlap (rad_m1_krylov.cpp, tests_m1/runs_3y_halo_overlap: round-off
  // vs off, restarts bitwise).  DEFAULT true since m1-accmerge wherever it is valid:
  // implicit_halo_mpi on and more than one rank; otherwise false, silently.  A restart
  // whose file lacks the key (written before m1-accmerge, which read it only when
  // named) keeps false; the resolved value is echoed; explicit input overrides.
  impl_halo_ovl = pin->GetOrAddBoolean("rad_m1","implicit_halo_overlap",
                      impl_halo_mpi && (global_variable::nranks > 1) &&
                      !global_variable::restart_run);
  // implicit_halo_ovl_faces (tests_m1/runs_4l_sync): default OFF since audit-m1-0928
  // (was true wherever the overlap is on, from 09-24 H200 weak scaling).  On MI300A with
  // moving gas it gains -0.1..+1.0 % (2 GPUs, one node) and -1.4..+0.8 % (2 nodes x 2
  // GPUs) per simulated second (He box, sph_wedge): no gain for a result-changing
  // (round-off) default.  Opt-in by name; the resolved value is echoed.
  impl_ovl_faces = pin->GetOrAddBoolean("rad_m1","implicit_halo_ovl_faces",false);
  for (int f = 0; f < 6; ++f) {hm_face[f] = 1;}
  if (impl_halo_ovl && !impl_halo_mpi) {
    ImplFatal("<rad_m1>/implicit_halo_overlap needs implicit_halo_mpi = true");
  }
  if (impl_kpipe && impl_kfuse != 3) {
    ImplFatal("<rad_m1>/implicit_krylov_pipe needs implicit_krylov_fuse = 3");
  }
  if (impl_halo_mpi && !impl_halo_direct) {
    ImplFatal("<rad_m1>/implicit_halo_mpi needs implicit_halo_direct = true");
  }
  // implicit_krylov_dev (tests_m1/runs_4k_launch): default 0, read only when named
  impl_kdev = 0;
  kdev_last[0] = kdev_last[1] = kdev_last[2] = 0;
  kdev_slot = 0;
  kdev_x1p = -1;
  kdev_nchk = 0.0;
  kdev_nq = 0.0;
  impl_kdev_halo = 0;
  if (pin->DoesParameterExist("rad_m1","implicit_krylov_dev")) {
    impl_kdev = pin->GetInteger("rad_m1","implicit_krylov_dev");
    if (impl_kdev < 0) {ImplFatal("<rad_m1>/implicit_krylov_dev must be >= 0");}
    if (impl_kdev > 0 && impl_prec >= 3) {
      ImplFatal("<rad_m1>/implicit_krylov_dev does not take implicit_precond = mg(_gc)");
    }
    if (impl_kdev > 0 && (impl_kfuse != 3 || impl_kpipe || !impl_halo_direct)) {
      ImplFatal("<rad_m1>/implicit_krylov_dev needs implicit_krylov_fuse = 3, "
                "implicit_halo_direct = true and implicit_krylov_pipe = false");
    }
    impl_kdev_halo = 0;
    if (impl_kdev > 0 && pin->DoesParameterExist("rad_m1","implicit_krylov_dev_halo")) {
      impl_kdev_halo = pin->GetInteger("rad_m1","implicit_krylov_dev_halo");
    }
    if (impl_kdev > 0) {
      kdv = DvceArray1D<Real>("m1_kdv", M1_KD_SIZE);
      kdh = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_kdh", M1_KD_SIZE);
    }
  }
  // implicit_op_check (rad_m1_opcheck.cpp, tests_m1/gates): read only when named
  impl_opchk = 0;
  opchk_n = 0;
  impl_opchk_tol = 1.0e-12;
  if (pin->DoesParameterExist("rad_m1","implicit_op_check")) {
    impl_opchk = pin->GetInteger("rad_m1","implicit_op_check");
    if (pin->DoesParameterExist("rad_m1","implicit_op_check_tol")) {
      impl_opchk_tol = pin->GetReal("rad_m1","implicit_op_check_tol");
    }
  }
  impl_dump_cyc = -1;
  impl_dump_done = false;
  if (pin->DoesParameterExist("rad_m1","implicit_dump_op")) {
    impl_dump_cyc = pin->GetInteger("rad_m1","implicit_dump_op");
  }
  // the Picard pass count (bench/m1_picard_0923): a per-pass log, off by default
  impl_plog = pin->GetOrAddInteger("rad_m1","implicit_picard_log",0);
  // implicit_tsolve_opac (nc_cure_1002): the Newton-fallback gas-temperature root find
  // with kappa_P(T), kappa_E(T) evaluated at the trial T (M1ImplTemperatureOpac) instead
  // of at the lagged iterate; it removes the |dT|/T = 0.365 limit cycles of shock-heated
  // near-void cells (He presn wedge, /viper/ptmp2/jinma/he_nc_cure_1002).  ...only from
  // Picard pass implicit_tsolve_opac_start on: a solve that converges before it is
  // bitwise the frozen-opacity one.  DEFAULT true with start 20 on fresh runs (user
  // 10-02; part of the R2 default with mode 1 and switch 10 below); a restart whose
  // file lacks the key keeps the old behaviour (off) and says so.  The resolved values
  // are recorded.
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac")) {
    impl_tsolve_opac = pin->GetBoolean("rad_m1","implicit_tsolve_opac");
  } else if (global_variable::restart_run) {
    impl_tsolve_opac = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_tsolve_opac: keeping the old "
                << "default (off); set it to switch" << std::endl;
    }
  } else {
    impl_tsolve_opac = pin->GetOrAddBoolean("rad_m1","implicit_tsolve_opac",true);
  }
  impl_tsolve_opac_start = pin->GetOrAddInteger("rad_m1","implicit_tsolve_opac_start",20);
  // implicit_tsolve_opac_mode (tsolve_root_fix_1003): 0 = the root find above; bit 1
  // (value 1) takes the root of the kappa(T) equation NEAREST the lagged T
  // (M1ImplTemperatureOpacNear) instead of the halve/double bracket's; bit 2 (value 2)
  // is a slope guard: a cell whose d ln kappa_P / d ln T at the lagged T is below
  // implicit_tsolve_opac_slope_min (default -4, where kappa_P T^4 falls with T) takes
  // the frozen-opacity root find instead.  3 = both.  Bit 4 (value 8, OPT-IN,
  // r2b_shift_1003) takes the stable root nearest the end point of the exact exchange
  // ODE from T_gn (M1ImplTemperatureOpacOde; wins over bit 1).  BSG R2E (mode 8, rst
  // 00008 on, bsg_1001/r2b_shift_1003): robust (no divergence) but ~10 % slower than
  // mode 1 (2.12 vs 1.93 s/cycle) with more NON-CONV (54 vs 31, worst resid 7.9e-4 vs
  // 3.8e-5), so the default stays 1.  DEFAULT 1 on fresh runs (R2
  // default, user 10-03): mode 0 picks the wrong root of the 3-root gas-T equation at H
  // recombination (kappa_P ~ T^30) and 2-cycles; with mode 1 the BSG envelope (rst
  // 00008, switch 10) runs 0 FATAL / 4 NON-CONV where mode 0 diverged
  // (bsg_1001/TSOLVE_FIX_1003.md, AB_SOLVER_1003.md), BSG R2 vs A0 lies within 2-3x the
  // noise spread (bsg_1001/r2b_shift_1003), He wedge/box: 0 NON-CONV, box bitwise
  // (he_mltpp_1002/GATE_MODE1.md).  A restart whose file lacks the key keeps the old
  // behaviour (0) and says so; a named value always wins.
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac_mode")) {
    impl_tsolve_opac_mode = pin->GetInteger("rad_m1","implicit_tsolve_opac_mode");
  } else if (global_variable::restart_run) {
    impl_tsolve_opac_mode = 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_tsolve_opac_mode: keeping "
                << "the old default (0); set it to switch" << std::endl;
    }
  } else {
    impl_tsolve_opac_mode = pin->GetOrAddInteger("rad_m1","implicit_tsolve_opac_mode",1);
  }
  if (pin->DoesParameterExist("rad_m1","implicit_tsolve_opac_slope_min")) {
    impl_tsolve_opac_smin = pin->GetReal("rad_m1","implicit_tsolve_opac_slope_min");
  }
  // DEBUG (nc_cure_1002): implicit_nc_dump = N > 0 prints, in the last N passes of a
  // solve that runs to implicit_maxit, the local state of each rank's worst-residual
  // cell (optical depths, coupling stiffness, energy ratio, beta, f, the T-Newton
  // flags).  Read only when named; prints only; no number of the solve changes.
  impl_ncdump = pin->DoesParameterExist("rad_m1","implicit_nc_dump") ?
                pin->GetInteger("rad_m1","implicit_nc_dump") : 0;
  // DEBUG (he_estall_1003): implicit_stall_trace = N > 0: in a solve still unconverged at
  // pass implicit_stall_trace_p0 (default 40), fix each rank's worst cell of that pass
  // and print its NCD line (tag 3) before step (c) and (tag 2) after the residual of
  // every pass p0 .. p0+N-1, for at most implicit_stall_trace_max (default 30) solves.
  // Prints only; read only when named.
  if (pin->DoesParameterExist("rad_m1","implicit_stall_trace")) {
    impl_strace = pin->GetInteger("rad_m1","implicit_stall_trace");
    impl_strace_p0 = pin->GetOrAddInteger("rad_m1","implicit_stall_trace_p0",40);
    impl_strace_max = pin->GetOrAddInteger("rad_m1","implicit_stall_trace_max",30);
  }
  impl_dtrace = pin->DoesParameterExist("rad_m1","implicit_det_trace") ?
                pin->GetInteger("rad_m1","implicit_det_trace") : 0;
  tmr_c0 = pin->DoesParameterExist("rad_m1","implicit_timers") ?
           pin->GetInteger("rad_m1","implicit_timers") : 0;   // read only when named
  // ...and the options that cut it.  Since bench/m1_defaults_0923 they DEFAULT ON for
  // the closures that do not read the iterate (eddington, vet_sc, tau: 1.6x on the GPU,
  // statistics moved at the solver-tolerance level) and stay OFF for m1 / minerbo /
  // kershaw, whose closure moves with the iterate: there the lresid test is NOT
  // redundant and the inexact (Eisenstat-Walker) pass 0 raised NON-CONVERGED steps
  // (tests_m1/runs_3k_gpu3d/README_DEFAULTS.md).  The earlier default path is
  // implicit_lres_test = true, implicit_conv_est = false, implicit_lin_ew_max = 0,
  // implicit_predictor = none (with implicit_bcg_sync = 0 and implicit_line_solver =
  // thomas it is bitwise the pre-flip code):
  //  implicit_lres_test = false  drop the pass-to-pass transverse-change test under
  //    bicgstab, where the transverse coupling is IN the operator and is solved to
  //    implicit_lin_tol on every pass (the test lags the Picard test by one pass);
  //  implicit_conv_est = true    stop once q/(1-q) times the last change (q = the
  //    measured contraction of the last two passes, required < 1/2) is below
  //    implicit_tol, i.e. without the confirming pass;
  //  implicit_lin_ew_max > 0     Eisenstat-Walker (choice 2) inner tolerance
  //    max(lin_tol max|b|, eta_k max|r0|), eta_0 = ew_max,
  //    eta_k = min(ew_max, gamma (|r0_k|/|r0_{k-1}|)^2) with the gamma eta_{k-1}^2
  //    safeguard (bcg_sync >= 1 only);
  //  implicit_predictor = step   start the loop from the previous step's implicit
  //    increment scaled by dt/dt_prev (closures that do not read the iterate only;
  //    the default for them since 0923, restart-safe).
  //  The Eisenstat-Walker default is 1e-2 with bcg_sync >= 1 and 0 (off) with
  //  bcg_sync = 0, so an input that asks for the original BiCGStab loop still runs.
  //  implicit_res_dmin, implicit_res_rmax (read only when named; default 0 = off): the
  //    Picard convergence norm (resid and lresid) leaves out the cells with gas density
  //    rho < res_dmin (needs <hydro>) or radius x1v > res_rmax (spherical polar only),
  //    e.g. the density-floor top whose residual stalls at 1e-7..1e-6.  Those cells are
  //    still solved on every pass; only the stopping test ignores them.  Off: bitwise.
  if (pin->DoesParameterExist("rad_m1","implicit_res_dmin")) {
    impl_res_dmin = pin->GetReal("rad_m1","implicit_res_dmin");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_res_rmax")) {
    impl_res_rmax = pin->GetReal("rad_m1","implicit_res_rmax");
  }
  //  implicit_realisable_coupling (gd_physfix_1005, read only when named; default false =
  //    bitwise the old code): the gas sees the REALISABLE radiation flux.  In optically
  //    thin cells with c dt/dx >> 1 the Picard iterate of the cell flux is unrealisable
  //    (BSG far-thin top: mean |F| = 10-20 c E, ~200 c E in the cells that stalled);
  //    the write-back clips it to |F| = c E (M1ApplyLimits), but the comoving correction
  //    E0 - E = -2 beta.F/c + O(beta^2) and the momentum deposit dt (rho k_t)_f F0_f/c
  //    took the UNCLIPPED flux: E0 - E came out -0.85 E (|E0 - E| <= 2 |beta| E for any
  //    intensity field), which moved the gas-T equation into the kappa_P(T) valley of the
  //    low-density table and made it 3-rooted (the Picard T cycle), and the floor gas got
  //    10-200x the largest force any radiation field of energy E can exert (rho k E).
  //    With the key, per cell, s = min(1, c E/|F_iter|) scales the flux in beta.F of E0,
  //    and the momentum deposit (and so its work, which the radiation loses: total
  //    energy stays exact).  s = 1 wherever the iterate is realisable.
  if (pin->DoesParameterExist("rad_m1","implicit_realisable_coupling")) {
    impl_real_couple = pin->GetBoolean("rad_m1","implicit_realisable_coupling");
  }
  //  implicit_thin_freeze (m1-picard-aa, 095ca6ee): cells with c dt rho kappa_P below it
  //    at pass 0 keep their start-of-solve opacities for every Picard pass and take the
  //    frozen-opacity gas-T find.  DEFAULT 1e-2 on fresh runs (user 10-04): on the He
  //    presn 128^2 wedge (he_mltpp_1002/FREEZE_AB.md, rst 00055, 1.5 ks) it is neutral
  //    within the cfl-0.2999 noise member below the photosphere, -32 % s/cycle, Picard
  //    mean 39.6 -> 16.3, 0 NON-CONV; dt -0.9 % (as on the BSG, PICARD_ROBUST_1003.md).
  //    A restart whose file lacks the key keeps the old behaviour (0 = off) and says so;
  //    a named value always wins.
  if (pin->DoesParameterExist("rad_m1","implicit_thin_freeze")) {
    impl_thin_frz = pin->GetReal("rad_m1","implicit_thin_freeze");
  } else if (global_variable::restart_run) {
    impl_thin_frz = 0.0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_thin_freeze: keeping the old "
                << "default (0 = off); set it to switch" << std::endl;
    }
  } else {
    impl_thin_frz = pin->GetOrAddReal("rad_m1","implicit_thin_freeze",1.0e-2);
  }
  if (impl_thin_frz < 0.0) {ImplFatal("<rad_m1>/implicit_thin_freeze must be >= 0");}
  //  implicit_resid_fatal: a Picard solve that ends NON-CONVERGED with resid > this
  //    value (or a non-finite resid or lin_resid) is a diverged solve: FATAL before its
  //    state spreads (BSG production 10-01: one 200-pass solve with resid 1.2e6 at
  //    t = 2.91e5 s, the whole domain NaN ~15 cycles later).  0 = off.  DEFAULT 1e2
  //    since defaults-1002 (/viper/ptmp2/jinma/defaults_gate_1002/item5): the largest
  //    NON-CONVERGED resid of a run that went on cleanly was 0.997 (He presn onv128;
  //    scout64 0.37); every resid > 1e3 seen came from a run that was failing.  A restart
  //    whose file lacks the key keeps the old default (off) and says so; the resolved
  //    value is recorded.
  if (pin->DoesParameterExist("rad_m1","implicit_resid_fatal")) {
    impl_res_fatal = pin->GetReal("rad_m1","implicit_resid_fatal");
    if (!(impl_res_fatal >= 0.0)) {
      ImplFatal("<rad_m1>/implicit_resid_fatal must be >= 0");
    }
  } else if (global_variable::restart_run) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_resid_fatal: keeping the old "
                << "default (off); set it to switch" << std::endl;
    }
  } else {
    impl_res_fatal = pin->GetOrAddReal("rad_m1","implicit_resid_fatal",1.0e2);
  }
  if (impl_res_dmin < 0.0 || impl_res_rmax < 0.0) {
    ImplFatal("<rad_m1>/implicit_res_dmin and implicit_res_rmax must be >= 0");
  }
  if (impl_res_rmax > 0.0 && !pmy_pack->pmesh->use_spherical_polar) {
    ImplFatal("<rad_m1>/implicit_res_rmax needs a spherical-polar mesh");
  }
  if (impl_res_dmin > 0.0 && pmy_pack->phydro == nullptr) {
    ImplFatal("<rad_m1>/implicit_res_dmin needs <hydro> (the gas density)");
  }
  impl_res_mask = (impl_res_dmin > 0.0) || (impl_res_rmax > 0.0);
  //  implicit_resid_fatal_masked (resmask_1003; read only when the mask is on, so a
  //    mask-off run is unchanged): FATAL when the largest residual of the excluded cells
  //    at the last pass of ANY solve (converged or not) is above this value or
  //    non-finite.  Default 1.0: the worst excluded residual of the BSG res_dmin = 1e-16
  //    A/B arm over 881 steps was 9.6e-2 (/viper/ptmp2/jinma/bsg_1001/SEAM_AB_1003.md).
  //    0 = off.
  if (impl_res_mask) {
    impl_res_fatal_masked = pin->GetOrAddReal("rad_m1","implicit_resid_fatal_masked",1.0);
    if (!(impl_res_fatal_masked >= 0.0)) {
      ImplFatal("<rad_m1>/implicit_resid_fatal_masked must be >= 0");
    }
  }
  impl_lres_test = pin->GetOrAddBoolean("rad_m1","implicit_lres_test",!fixcl);
  impl_conv_est = pin->GetOrAddBoolean("rad_m1","implicit_conv_est",fixcl);
  // implicit_src_stable (stall_1002): the cancellation-free form of the gas-eliminated
  // source row.  DEFAULT true (user 10-02) for a fresh run; false (the plain form,
  // bitwise) on a restart whose file lacks the key, with one line on rank 0.  The value
  // is recorded, so later restarts keep it; a named key always wins.
  {const bool named = pin->DoesParameterExist("rad_m1","implicit_src_stable");
  impl_src_stable = pin->GetOrAddBoolean("rad_m1","implicit_src_stable",
                                         !global_variable::restart_run);
  if (!named && global_variable::restart_run && global_variable::my_rank == 0) {
    std::cout << "rad_m1: implicit_src_stable = false (restart file lacks the key; "
              << "name it to switch the stable source row on)" << std::endl;
  }
  }
  impl_ew_max = pin->GetOrAddReal("rad_m1","implicit_lin_ew_max",
                                  (fixcl && impl_bcg_sync >= 1) ? 1.0e-2 : 0.0);
  impl_ew_gam = pin->GetOrAddReal("rad_m1","implicit_lin_ew_gamma",0.9);
  //  implicit_lin_cnorm > 0     the inner test on max_i |r_i|/(s_i E^k_i) < cnorm,
  //    s_i = 1 + SRCB_i the row excess: a bound on the per-cell relative error of E
  //    the residual implies (bcg_sync >= 1 only).
  impl_lin_cnorm = pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0);
  if (impl_lin_cnorm < 0.0 || (impl_lin_cnorm > 0.0 && impl_bcg_sync == 0)) {
    ImplFatal("<rad_m1>/implicit_lin_cnorm must be >= 0 and needs "
              "implicit_bcg_sync >= 1");
  }
  // implicit_bcg_max_restarts / implicit_bcg_fallback (sp-blend2-1008; defaults 2 / line
  // = the original, bitwise).  After more than max_restarts breakdowns (a vanishing
  // rho / omega / rhat.v, or a recurrence residual that met the tolerance while the true
  // residual did not) the fused BiCGStab gives up.  `line` then DISCARDS the Krylov
  // iterate and makes one line-Jacobi update from E^k (ImplicitBiCGStabEnd); `best`
  // first forms the true residual of the iterate and keeps the iterate whenever that
  // residual is smaller than the initial one (the same norm the tolerance uses), so the
  // Picard pass gets the better of the two instead of always the cruder.  Fused path
  // (bcg_sync >= 1, krylov_fuse <= 2) only.
  impl_bcg_maxrst = pin->GetOrAddInteger("rad_m1","implicit_bcg_max_restarts",2);
  {std::string fb = pin->GetOrAddString("rad_m1","implicit_bcg_fallback","line");
  if (fb.compare("line") == 0) {
    impl_bcg_keep = false;
  } else if (fb.compare("best") == 0) {
    impl_bcg_keep = true;
  } else {
    ImplFatal("<rad_m1>/implicit_bcg_fallback = '" + fb
              + "' is not a choice (line | best)");
  }
  }
  if (impl_bcg_maxrst < 0) {
    ImplFatal("<rad_m1>/implicit_bcg_max_restarts must be >= 0");
  }
  // implicit_face_opac_n (fixbundle-1009 F2): rad_m1.hpp.  defaults-1009 (user 10-09):
  // default true on a fresh run; a restart whose file lacks the key keeps false
  impl_face_ktn = pin->GetOrAddBoolean("rad_m1","implicit_face_opac_n",
                                       !global_variable::restart_run);
  if (impl_ew_max < 0.0 || impl_ew_max >= 1.0 || !(impl_ew_gam > 0.0)) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max must lie in [0,1), ew_gamma > 0");
  }
  if (impl_ew_max > 0.0 && impl_bcg_sync == 0) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max needs implicit_bcg_sync >= 1");
  }
  // implicit_det_reduce (bsg_1001/DETERMINISM.md): the Krylov sums of the fused
  // BiCGStab in a fixed order (M1DetReduce) and the cross-rank sums in rank order, so
  // that two runs of one binary on one decomposition are bitwise identical on the GPU.
  // DEFAULT: on for a device build (GPU: measured 2-7 % CHEAPER per cycle than the
  // Kokkos reductions on the BSG 3-D and He ext wedges, DETERMINISM.md) wherever it is
  // supported, off on a host build (the strided slots cost 1.7x on the CPU column).
  // Both settings are run-to-run bitwise on the GPU for a fixed decomposition.
  const bool det_ok = (impl_solver == M1_ISOLV_BICGSTAB) && (impl_bcg_sync == 1) &&
                      (impl_kfuse != 2) && !impl_kpipe && (impl_kdev == 0) &&
                      (impl_accel == M1_IACC_NONE);
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA) \
    || defined(KOKKOS_ENABLE_SYCL)
  const bool det_def = det_ok;
#else
  const bool det_def = false;
#endif
  impl_det = pin->GetOrAddBoolean("rad_m1","implicit_det_reduce",det_def);
  if (impl_det && !det_ok) {
    ImplFatal("<rad_m1>/implicit_det_reduce supports implicit_solver = bicgstab with "
              "implicit_bcg_sync = 1, implicit_krylov_fuse = 0, 1 or 3 (not krylov_pipe "
              "or krylov_dev) and implicit_accel = none "
              "(solver " + std::to_string(impl_solver) + ", bcg_sync " +
              std::to_string(impl_bcg_sync) + ", krylov_fuse " +
              std::to_string(impl_kfuse) + ", accel " + std::to_string(impl_accel) + ")");
  }
  if (impl_det) {det_part = DvceArray1D<Real>("m1_det_part", 4*1024 + 4);}
  {std::string pr = pin->GetOrAddString("rad_m1","implicit_predictor",
                                        fixcl ? "step" : "none");
  // predictor = step is the default for the closures that do not read the iterate: its
  // state (ipred, pred_ok, pred_dt) travels in the restart file (radm1::kM1PredRstMagic,
  // tests_m1/runs_3k_gpu3d/README_PREDRST.md), so a restart stays bitwise.
  if (pr.compare("none") == 0) {
    impl_pred = false;
  } else if (pr.compare("step") == 0) {
    impl_pred = true;
  } else {
    ImplFatal("<rad_m1>/implicit_predictor = '" + pr
              + "' is not a choice (none | step)");
  }
  }
  if (impl_ecnt < 0 || impl_ecnt > M1_EC_NTMAX) {
    ImplFatal("<rad_m1>/implicit_eos_cache_nt must lie in [0,8]");
  }
  impl_and_m = pin->GetOrAddInteger("rad_m1","implicit_anderson_m",5);
  impl_and_beta = pin->GetOrAddReal("rad_m1","implicit_anderson_beta",1.0);
  impl_and_start = pin->GetOrAddInteger("rad_m1","implicit_anderson_start",1);
  if (impl_accel == M1_IACC_ANDERSON) {
    if (impl_and_m < 1 || impl_and_m > M1_AND_MMAX) {
      ImplFatal("<rad_m1>/implicit_anderson_m must lie in [1,10]");
    }
    if (!(impl_and_beta > 0.0) || impl_and_beta > 1.0) {
      ImplFatal("<rad_m1>/implicit_anderson_beta must lie in (0,1]");
    }
    if (impl_and_start < 1) {
      ImplFatal("<rad_m1>/implicit_anderson_start must be >= 1 (pass 0 has no history)");
    }
  }
  // implicit_opac_update: default true with force_reference = wb_arad
  // (m1-keydefault-0927, docs/dev/ke_dt_0926.md 8), false otherwise
  impl_opac_update = pin->GetOrAddBoolean("rad_m1","implicit_opac_update",
                                          force_ref == M1_FREF_WB_ARAD);
  // DIAGNOSTIC (m1-perf-0928): which opacities the Picard-loop update moves
  dbg_opac_part = pin->GetOrAddInteger("rad_m1","dbg_opac_part",0);
  // implicit_opac_newton (m1-perf-0928): the Picard loop under implicit_opac_update is a
  // fixed-point iteration on the FLUX opacity rho kappa_T(T) of the x1 face terms, and
  // it converges only linearly (~0.17 per pass on the He box: 6.8 passes per solve at
  // tol 1e-8 against 2.0 with the opacities frozen).  This key adds the Newton term of
  // that dependence to the x1 rows (the chain dG/dkappa_face dkappa/dT dT/dE' through
  // the local gas Newton step), so the pass converges ~quadratically.  The added terms
  // are proportional to the predicted temperature change of the gas Newton step, which
  // vanishes at the fixed point: the converged state is the same.
  // DEFAULT true (m1-perf-0928, user) wherever it applies: implicit_opac_update and
  // implicit_gas_newton on; hydro or MHD (FluidRef), Cartesian or sp rows.  False on a
  // restart whose file lacks the key (the implicit_gas_newton convention).  Named true
  // without its two prerequisites it is fatal (below).
  impl_opac_newton = pin->GetOrAddBoolean("rad_m1","implicit_opac_newton",
                                          impl_opac_update && impl_gas_newton &&
                                          !global_variable::restart_run);
  // implicit_opac_newton_guard (m1-opn-guard, 09-29): the Newton term above has the sign
  // of the face flux times d kappa/dT and no bound; in fast optically thin cells (He
  // presn wedge, across the Fe bump) it took a row diagonal of +802 to -9200, the row
  // lost the M-matrix property and the solved E went negative (floor clips, stage
  // failures, NaN).  A face's term is kept in a row only if the diagonal stays >= guard
  // x its value without it; otherwise that face is left Picard in that row (the fixed
  // point is the same: the term vanishes there).  <= 0 turns the guard off (the
  // unguarded rows of m1-perf-0928, bitwise).
  impl_opn_guard = pin->GetOrAddReal("rad_m1","implicit_opac_newton_guard",0.5);
  // m1-positivity set: DEFAULT ON since defaults-1002 (implicit_g0_exchange = true,
  // implicit_g0_limit = 1, implicit_pos_gas = true, implicit_pos_floor = true and
  // implicit_opac_newton_guard_mode = 6 below; every he_star_m1 production ran this set;
  // He box / wedge IDENT, BSG column within the 1-ulp twin spread:
  // /viper/ptmp2/jinma/defaults_gate_1002/item2).  A restart whose file lacks a key
  // keeps its old default (off; guard mode 2 is recorded in every file) and says so once.
  // The resolved values are recorded.  implicit_pos_gas_frac: read only when named.
  {const bool rst = global_variable::restart_run;
  bool old_kept = false;
  auto pbool = [&](const char *key, bool &v) {
    if (pin->DoesParameterExist("rad_m1",key)) {
      v = pin->GetBoolean("rad_m1",key);
    } else if (rst) {
      old_kept = true;
    } else {
      v = pin->GetOrAddBoolean("rad_m1",key,true);
    }
  };
  if (pin->DoesParameterExist("rad_m1","implicit_g0_limit")) {
    impl_g0_lim = pin->GetReal("rad_m1","implicit_g0_limit");
  } else if (rst) {
    old_kept = true;
  } else {
    impl_g0_lim = pin->GetOrAddReal("rad_m1","implicit_g0_limit",1.0);
  }
  pbool("implicit_g0_exchange", impl_g0_exch);
  pbool("implicit_pos_gas", impl_pos_gas);
  pbool("implicit_pos_floor", impl_pos_floor);
  // defaults-1009 (user 10-09): default = implicit_pos_floor on a fresh run (it needs
  // it); a restart whose file lacks the key keeps false
  impl_pos_floor_s2 = pin->GetOrAddBoolean("rad_m1","implicit_pos_floor_solve",
                                           impl_pos_floor && !rst);
  if (impl_pos_floor_s2 && !impl_pos_floor) {
    ImplFatal("<rad_m1>/implicit_pos_floor_solve needs implicit_pos_floor = true");
  }
  if (old_kept && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> restart input lacks m1-positivity keys (implicit_g0_exchange,"
              << " g0_limit, pos_gas, pos_floor): keeping the old default (off) for "
              << "those; name them to switch" << std::endl;
  }
  }  // end of the restart-default scope
  if (pin->DoesParameterExist("rad_m1","implicit_pos_gas_frac")) {
    impl_pos_gas_frac = pin->GetReal("rad_m1","implicit_pos_gas_frac");
  }
  Kokkos::realloc(pos_cnt_d, M1_POS_N);
  Kokkos::deep_copy(pos_cnt_d, 0.0);
  // implicit_opac_newton_guard_mode: 0 = the diagonal test alone; + 1 = also the
  // right-hand side, + 2 = also the sign of the neighbour entry (M1OpnGuardFace; He presn
  // wedge gate: 7 floor clips with the diagonal alone, 0 with both), + 4 (m1-positivity)
  // = also the row diagonal dominance (m1_impl_asm).  DEFAULT 6 since defaults-1002 (the
  // positivity set above); 2 on a restart whose file lacks the key (files written since
  // 09-29 record it).
  impl_opn_guard_mode = pin->GetOrAddInteger("rad_m1","implicit_opac_newton_guard_mode",
                                             global_variable::restart_run ? 2 : 6);
  // opacnewt-cliff-1007: implicit_opac_newton_slope_off / _slope_max (read only when
  // named; 0 = off).  At an opacity cliff (He recombination: |d ln kappa/d ln T| ~ 5-14
  // in the TOPS tables) the d(rho kappa_T)/dT term of a strongly coupled cell makes the E
  // iterate cycle (He giant N897 photosphere, HEGIANT_PICARD.md: 77 passes, NON-CONV);
  // these keys drop or damp the term in cliff-steep cells only.
  if (pin->DoesParameterExist("rad_m1","implicit_opac_newton_slope_off")) {
    impl_opn_soff = pin->GetReal("rad_m1","implicit_opac_newton_slope_off");
  }
  if (pin->DoesParameterExist("rad_m1","implicit_opac_newton_slope_max")) {
    impl_opn_smax = pin->GetReal("rad_m1","implicit_opac_newton_slope_max");
  }
  impl_allow_multid = pin->GetOrAddBoolean("rad_m1","implicit_allow_multid",false);
  marshak_q = pin->GetOrAddReal("rad_m1","marshak_q",0.5);
  // implicit_marshak_face (m1-sp-order2, tests_m1/runs_5o_sporder2).  DEFAULT linear on
  // the spherical-polar wedge since m1-sp-order2b (tests_m1/runs_5q_sporder2b) with the
  // fixed-tensor closures (eddington, vet_sc, tau, vet_col); cell with the lagged
  // closures (m1, minerbo, kershaw: the runs_5h pp_np atmosphere, closure_lag = step at
  // implicit_cfl 1e6, went from 30 round-off NON-CONVERGED solves to 54 with Picard
  // blow-ups).  A restart whose file lacks the key keeps cell; the resolved value is
  // echoed.  Elsewhere read only when named (parameter dump and the cell form kept).
  if (sph_geom || pin->DoesParameterExist("rad_m1","implicit_marshak_face")) {
    const bool mfdef = !global_variable::restart_run &&
                       (eddington || vet_sc || tau_closure);
    const std::string smf = sph_geom ?
        pin->GetOrAddString("rad_m1","implicit_marshak_face", mfdef ? "linear" : "cell") :
        pin->GetString("rad_m1","implicit_marshak_face");
    if (smf.compare("cell") == 0) {
      impl_mface_lin = false;
    } else if (smf.compare("linear") == 0) {
      impl_mface_lin = true;
    } else {
      ImplFatal("<rad_m1>/implicit_marshak_face = '" + smf
                + "' is not a choice (cell | linear)");
    }
    if (impl_mface_lin && !sph_geom) {
      ImplFatal("<rad_m1>/implicit_marshak_face = linear is implemented on the "
                "spherical-polar wedge only");
    }
  }
  // implicit_face_weight (see rad_m1.hpp): distance (default on the spherical-polar
  // wedge) | equal (default elsewhere; the only choice on a Cartesian grid, which is
  // uniform).  On the wedge the key is recorded (GetOrAdd), so a restart keeps what the
  // run started with.  A RESTART whose embedded input predates the key keeps `equal`,
  // the behaviour it was run with, unless the key is given on the command line or in
  // an -i overlay (ModifyFromCmdline adds absent keys on a restart).  On a Cartesian grid
  // it is read only when named, so the parameter dump there is unchanged.
  if (sph_geom || pin->DoesParameterExist("rad_m1","implicit_face_weight")) {
    const bool named = pin->DoesParameterExist("rad_m1","implicit_face_weight");
    const std::string sfw = pin->GetOrAddString("rad_m1","implicit_face_weight",
        global_variable::restart_run ? "equal" : "distance");
    if (sfw.compare("equal") == 0) {
      impl_face_wdist = false;
    } else if (sfw.compare("distance") == 0) {
      impl_face_wdist = true;
    } else {
      ImplFatal("<rad_m1>/implicit_face_weight = '" + sfw
                + "' is not a choice (equal | distance)");
    }
    if (impl_face_wdist && !sph_geom) {
      ImplFatal("<rad_m1>/implicit_face_weight = distance is implemented on the "
                "spherical-polar wedge only (a Cartesian grid is uniform)");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1>/implicit_face_weight = " << sfw
                << (named ? " (from the input)"
                    : (global_variable::restart_run
                       ? " (restart file predates the key: kept the old equal)"
                       : " (default)")) << std::endl;
    }
  }
  // ---- milestone 3a2 options.  All three default to the 3a behaviour, so an input file
  // that does not name them reproduces RESULTS.txt of runs_3a exactly.
  std::string sfx = pin->GetOrAddString("rad_m1","implicit_flux","central");
  if (sfx.compare("central") == 0) {
    impl_flux = M1_IFLUX_CENTRAL;
  } else if (sfx.compare("ap_hll") == 0) {
    impl_flux = M1_IFLUX_APHLL;
  } else if (sfx.compare("berthon") == 0) {
    impl_flux = M1_IFLUX_BERTHON;
  } else if (sfx.compare("blend") == 0) {
    impl_flux = M1_IFLUX_BLEND;
  } else {
    ImplFatal("<rad_m1>/implicit_flux = '" + sfx
              + "' is not a choice (central | ap_hll | berthon | blend)");
  }
  // 3b phase B: the transverse operator is built for the face-eliminated (central) form
  // only.  The HLL/berthon/blend face fluxes carry per-face coefficients (ifw) that exist
  // for the x1 faces alone, so anything but `central` would silently be central in x2/x3.
  // sp-blend-1008: on the spherical-polar WEDGE (not the cubed sphere) berthon | blend
  // act on the x1 (radial) faces, with the A_f/V_i areas and volumes and the centroid
  // distance of the sp row; the x2/x3 faces stay central (a beam in an outflow is
  // radial).  Every other multi-D mesh keeps the refusal.
  // blendall-1009: implicit_flux_faces = all extends berthon | blend to EVERY face of the
  // multi-D solve -- the Cartesian x1, x2, x3 faces and the sp LATERAL faces -- with the
  // per-face coefficients of the lagged closure's beam vector (M1BeamFace).  Default
  // `x1`: the behaviour above, bitwise.
  {
    std::string sff = pin->GetOrAddString("rad_m1","implicit_flux_faces","x1");
    if (sff.compare("all") == 0) {
      impl_flux_all = true;
    } else if (sff.compare("x1") != 0) {
      ImplFatal("<rad_m1>/implicit_flux_faces = '" + sff
                + "' is not a choice (x1 | all)");
    }
  }
  {
    std::string sbb = pin->GetOrAddString("rad_m1","implicit_flux_beam","closure");
    if (sbb.compare("fs") == 0) {
      impl_beam_fs = true;
    } else if (sbb.compare("halfrange") == 0) {
      impl_beam_hr = true;
      // without ray data (closure m1 / eddington / vet_col without vet_gd) the faces
      // stay central unless implicit_hr_model = true (the isotropic + beam model of the
      // lagged closure: it is NEUTRAL in a thick layer -- the lagged flux carries the
      // whole flux with no gradient -- and fails G1, tests: blendall_1009/gates/g1)
      impl_hr_model = pin->GetOrAddBoolean("rad_m1","implicit_hr_model",false);
    } else if (sbb.compare("closure") != 0) {
      ImplFatal("<rad_m1>/implicit_flux_beam = '" + sbb + "' is not a choice "
                "(closure | fs | halfrange)");
    }
  }
  if (full && impl_flux_all) {
    if (!(impl_flux == M1_IFLUX_BERTHON || impl_flux == M1_IFLUX_BLEND) || cs_geom) {
      ImplFatal("<rad_m1>/implicit_flux_faces = all needs implicit_flux = berthon | "
                "blend "
                "on a Cartesian or spherical-polar (not cubed-sphere) mesh");
    }
  }
  if (full && impl_flux != M1_IFLUX_CENTRAL) {
    const bool spok = (sph_geom || impl_flux_all) && !cs_geom &&
                      (impl_flux == M1_IFLUX_BERTHON || impl_flux == M1_IFLUX_BLEND);
    if (!spok) {
      ImplFatal("<rad_m1>/transport = implicit supports implicit_flux = central only "
                "(berthon | blend on the x1 faces of the spherical-polar wedge); use "
                "transport = implicit_x1 for ap_hll | berthon | blend");
    }
  }
  // ---- milestone 3c.  The weight of implicit_flux = blend.  Inert for every other
  // flux, and the two ends of the blend are BITWISE central and berthon.
  std::string sbl = pin->GetOrAddString("rad_m1","implicit_blend","tau_f");
  if (sbl.compare("tau") == 0) {
    impl_blend = M1_IBLEND_TAU;
  } else if (sbl.compare("f") == 0) {
    impl_blend = M1_IBLEND_F;
  } else if (sbl.compare("tau_f") == 0) {
    impl_blend = M1_IBLEND_TAUF;
  } else if (sbl.compare("idort") == 0) {
    impl_blend = M1_IBLEND_IDORT;   // hrup-1009, half-range only
  } else if (sbl.compare("knudsen") == 0) {
    impl_blend = M1_IBLEND_KN;      // hrup-1009, half-range only
  } else if (sbl.compare("idort_f") == 0) {
    impl_blend = M1_IBLEND_IDF;     // hrup-1009, half-range only
  } else if (sbl.compare("idort_a") == 0) {
    impl_blend = M1_IBLEND_IDA;     // hrup-1009, half-range only
  } else {
    ImplFatal("<rad_m1>/implicit_blend = '" + sbl
              + "' is not a choice (tau | f | tau_f | idort | knudsen | idort_f | "
              "idort_a)");
  }
  if (impl_blend >= M1_IBLEND_IDORT &&
      pin->GetOrAddString("rad_m1","implicit_flux_beam",
                          "closure").compare("halfrange") != 0) {
    ImplFatal("<rad_m1>/implicit_blend = idort | knudsen | idort_f | idort_a needs "
              "implicit_flux_beam = halfrange");
  }
  impl_blend_alpha = pin->GetOrAddReal("rad_m1","implicit_blend_alpha",1.0);
  impl_blend_r0 = pin->GetOrAddReal("rad_m1","implicit_blend_r0",1.5);
  impl_blend_xthin = pin->GetOrAddReal("rad_m1","implicit_blend_xthin",0.0);
  {
    // xthinfix-1009: where the xthin override acts (all faces = hrup-1009 behaviour)
    std::string sxm = pin->GetOrAddString("rad_m1","implicit_blend_xthin_mode","all");
    if (sxm.compare("all") == 0) {
      impl_blend_xthin_mode = 0;
    } else if (sxm.compare("beam") == 0) {
      impl_blend_xthin_mode = 1;
    } else if (sxm.compare("beam_kn") == 0) {
      impl_blend_xthin_mode = 2;
    } else if (sxm.compare("steep") == 0) {
      impl_blend_xthin_mode = 3;
    } else {
      ImplFatal("<rad_m1>/implicit_blend_xthin_mode = '" + sxm
                + "' is not a choice (all | beam | beam_kn | steep)");
    }
    if (impl_blend_xthin_mode == 2) {
      impl_blend_xthin_r0 = pin->GetOrAddReal("rad_m1","implicit_blend_xthin_r0",1.5);
    }
    impl_blend_xthin_wmin = pin->GetOrAddReal("rad_m1","implicit_blend_xthin_wmin",0.0);
  }
  if (pin->DoesParameterExist("rad_m1","implicit_hr_recon")) {
    // xthinfix-1009 Fix B: read only when given (keys-off parameter dump unchanged)
    std::string shr = pin->GetString("rad_m1","implicit_hr_recon");
    if (shr.compare("plm") == 0) {
      impl_muscl = true;
    } else if (shr.compare("dc") != 0) {
      ImplFatal("<rad_m1>/implicit_hr_recon = '" + shr + "' is not a choice (dc | plm)");
    }
    if (impl_muscl) {
      impl_muscl_nfresh = pin->GetOrAddInteger("rad_m1","implicit_hr_recon_fresh",2);
      impl_muscl_damp = pin->GetOrAddReal("rad_m1","implicit_hr_damp",0.0);
    }
    // the plm face states multiply the face coefficients HCL/HCR of ANY upwind part:
    // the half-range flux, or the berthon / blend AP-HLL part (implicit_flux_beam =
    // closure | fs); a central face has none
    if (impl_muscl && pin->GetOrAddString("rad_m1","implicit_flux",
                                          "central").compare("central") == 0) {
      ImplFatal("<rad_m1>/implicit_hr_recon = plm needs implicit_flux = berthon | blend");
    }
    if (impl_muscl && cs_geom) {
      ImplFatal("<rad_m1>/implicit_hr_recon = plm is not implemented on the cubed "
                "sphere");
    }
  }
  std::string sbm = pin->GetOrAddString("rad_m1","implicit_blend_fmode","max");
  if (sbm.compare("max") == 0) {
    impl_blend_fmode = M1_IBFM_MAX;
  } else if (sbm.compare("mean") == 0) {
    impl_blend_fmode = M1_IBFM_MEAN;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_fmode = '" + sbm
              + "' is not a choice (max | mean)");
  }
  std::string sbw = pin->GetOrAddString("rad_m1","implicit_blend_mode","flux");
  if (sbw.compare("flux") == 0) {
    impl_blend_mode = M1_IBMODE_FLUX;
  } else if (sbw.compare("dissipation") == 0) {
    impl_blend_mode = M1_IBMODE_DISSIP;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_mode = '" + sbw
              + "' is not a choice (flux | dissipation)");
  }
  // sp-blend-1008: the tau-only weight fails the grey atmosphere, the opacity jump and
  // the He column (design sect. 9), and `dissipation` is rejected there too; neither is
  // offered on the multi-D wedge
  if (impl_beam_hr) {
    // hrup-1009: the half-range face flux takes the AP weight of the face optical depth
    // alone (implicit_blend = tau), or w = 1 (berthon); no f-gate: the half-range ratios
    // ARE the angular information
    if (!((impl_flux == M1_IFLUX_BLEND && impl_blend != M1_IBLEND_F &&
           impl_blend_mode == M1_IBMODE_FLUX) || impl_flux == M1_IFLUX_BERTHON) ||
        (full && pmy_pack->pmesh->multi_d && !impl_flux_all)) {
      ImplFatal("<rad_m1>/implicit_flux_beam = halfrange needs implicit_flux = "
                "berthon, or "
                "blend with implicit_blend = tau | tau_f | idort and "
                "implicit_blend_mode = "
                "flux; on a multi-D mesh also implicit_flux_faces = all");
    }
  }
  if (full && impl_flux == M1_IFLUX_BLEND && !impl_beam_hr &&
      (impl_blend == M1_IBLEND_TAU || impl_blend_mode == M1_IBMODE_DISSIP)) {
    ImplFatal("<rad_m1>/transport = implicit with implicit_flux = blend takes "
              "implicit_blend = f | tau_f and implicit_blend_mode = flux only");
  }
  impl_blend_tau0 = pin->GetOrAddReal("rad_m1","implicit_blend_tau0",1.0);
  impl_blend_flo = pin->GetOrAddReal("rad_m1","implicit_blend_flo",0.6);
  impl_blend_fhi = pin->GetOrAddReal("rad_m1","implicit_blend_fhi",0.9);
  if (!(impl_blend_tau0 > 0.0) || !(impl_blend_fhi > impl_blend_flo)) {
    ImplFatal("<rad_m1>: implicit_blend_tau0 must be positive and implicit_blend_fhi "
              "must exceed implicit_blend_flo");
  }
  std::string srn = pin->GetOrAddString("rad_m1","implicit_recon","dc");
  if (srn.compare("dc") == 0) {
    impl_recon = M1_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = M1_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_m1>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // implicit_enthalpy (see rad_m1_implicit.hpp).  DEFAULT plm since m1-defaults for the
  // closures fixed within a step (eddington, vet_sc, tau; tests_m1/runs_3s_space2: second
  // order, fewer Picard passes, NON-CONVERGED 0); m1 / minerbo / kershaw keep upwind
  // (not gated; tests_m1/runs_4b_defaults).  Up to d0c59f7c the key was read only
  // when named, so a restart file written then carries it only if the input named it;
  // such a restart that does NOT carry it keeps the old default, upwind.  The resolved
  // value is echoed, so later restarts keep it.
  const bool enth_up = global_variable::restart_run || !fixcl;
  const std::string sen = pin->GetOrAddString("rad_m1","implicit_enthalpy",
                                              enth_up ? "upwind" : "plm");
  if (sen.compare("upwind") == 0) {
    impl_enth = M1_IENTH_UPWIND;
  } else if (sen.compare("central") == 0) {
    impl_enth = M1_IENTH_CENTRAL;
  } else if (sen.compare("plm") == 0) {
    impl_enth = M1_IENTH_PLM;
  } else {
    ImplFatal("<rad_m1>/implicit_enthalpy = '" + sen
              + "' is not a choice (upwind | central | plm)");
  }
  // implicit_vimp (rad_m1_implicit.hpp).  DEFAULT on since m1-defaults2
  // (tests_m1/runs_5g_defaults2) only where the resolved time_scheme is hesdirk2, with
  // which it was validated (under be it over-damps P = 100 waves), and where it is valid:
  // a multi-D mesh with implicit_solver = bicgstab, nghost >= 2, hydro with coupling,
  // gas_feedback and dbg_gas_force.  Otherwise off, silently.  A restart whose file lacks
  // the key (read only when named before) keeps off; the resolved value is echoed;
  // explicit input overrides.
  {auto *pmh = pmy_pack->pmesh;
  const bool vdef = (pin->GetString("rad_m1","time_scheme").compare("hesdirk2") == 0) &&
                    full && pmh->multi_d && (impl_solver == M1_ISOLV_BICGSTAB) &&
                    (pmh->mb_indcs.ng >= 2) && fl_on &&
                    coupling && gas_feedback && dbg_gas_force &&
                    !global_variable::restart_run;
  impl_vimp = pin->GetOrAddBoolean("rad_m1","implicit_vimp",vdef);}
  // DIAGNOSTIC: a scale of the Jacobian P (1 = Newton).  The converged state does not
  // depend on it; only the Picard contraction does.
  impl_vimp_jscale = 1.0;
  if (pin->DoesParameterExist("rad_m1","implicit_vimp_jscale")) {
    impl_vimp_jscale = pin->GetReal("rad_m1","implicit_vimp_jscale");
  }
  // ---- THE runs_4a_accel LEVERS (tests_m1/runs_4a_accel, tests_m1/runs_4j_accmerge).
  // DEFAULT ON since m1-accmerge for time_scheme = be (and since m1-h2fast for hesdirk2,
  // tests_m1/runs_5f_h2fast, whose stage solves also take time2_one_pass_safety and
  // time2_lin_tol_fac, rad_m1_time2.cpp) with transport = implicit and a
  // closure whose tensor is fixed within a step (eddington, vet_sc, tau), wherever each
  // is valid:
  //   implicit_fast_kernels = true   (bitwise-exact kernel shortcuts)
  //   implicit_vimp_fold = true      where implicit_vimp and implicit_op_stencil are on
  //   implicit_one_pass = 8          where the predictor is on
  //   implicit_predictor_order = 2   where the predictor is on
  // (1 GPU, 3-D He box: be 51.6 -> 43.5 ms/cycle Eddington, 57.5 -> 50.7 vet_sc with
  // one_pass = 4; one_pass = 8 is the faster of the two in runs_4a_accel).  Up to the
  // merge of m1-accel the keys were read only when named, so a restart file written
  // before this default carries them only if its input named them; such a restart that
  // does NOT carry a key keeps the old value (off / 0 / 1).  The resolved values are
  // always echoed, so later restarts keep them.  Explicit input values override.
  //  implicit_vimp_fold: fold the implicit_vimp operator part into the stored stencil
  //    (x2/x3 +-1 into slots 3-6, +-2 neighbours into slots 19-24) instead of
  //    M1VimpRow per apply (round-off).  Needs implicit_op_stencil.
  //  implicit_one_pass = N: see ImplicitSolve; its state travels in M1ONEP01.
  //  implicit_predictor_order = 2: the predictor extrapolates the increment rate g =
  //    dE/dt linearly in time, g* = g1 + h dt1 with h = (g1 - g2)/dt2, h stored per cell
  //    in ipred channels 3 (E) and 4 (T).  Under hesdirk2 the backward-Euler steps (first
  //    step, fallbacks) then leave ipred alone, so that it holds stage-1 increments only.
  //  implicit_fast_kernels: closure = eddington: D_ab = 0 off the diagonal, so the
  //    stencil build and the right-hand side skip the off-diagonal (od) terms, which are
  //    exactly zero; hesdirk2 + vet_sc without extrapolation (time2_vet_extrap = false):
  //    the tensor save of Time2VetExtrapolate is a plain copy.
  {
  // m1-h2fast (tests_m1/runs_5f_h2fast): the same defaults for time_scheme = hesdirk2,
  // whose stage solves already track the one_pass q and the predictor per solve kind
  bool ts_ok = true;
  if (pin->DoesParameterExist("rad_m1","time_scheme")) {
    const std::string tsn = pin->GetString("rad_m1","time_scheme");
    ts_ok = (tsn.compare("be") == 0) || (tsn.compare("hesdirk2") == 0);
  }
  const bool ldef = full && ts_ok && fixcl && !global_variable::restart_run;
  impl_fastk = pin->GetOrAddBoolean("rad_m1","implicit_fast_kernels",ldef);
  // implicit_vimp_fold: default OFF since audit-m1-0928 (gain 0.9-2.3 % per simulated
  // second on the moving He box and sph_wedge, MI300A 1 and 2 GPUs; under the ~3 % bar
  // for a result-changing default).  Opt-in by name.
  impl_vfold = pin->GetOrAddBoolean("rad_m1","implicit_vimp_fold",false);
  // implicit_one_pass: default 0 (off) everywhere since audit-m1-0928 (8 + auto gained
  // -2.8..+0.0 % on the moving He box and sph_wedge on MI300A, i.e. nothing, and changes
  // results 5-22x the round-off spread; before, 0 only with force_reference = wb_arad,
  // m1-keydefault-0927).  Opt-in by name (e.g. 8).
  impl_onep = pin->GetOrAddInteger("rad_m1","implicit_one_pass",0);
  if (impl_onep != 0 && impl_onep < 2) {
    ImplFatal("<rad_m1>/implicit_one_pass (the check period) must be 0 (off) or >= 2");
  }
  if (impl_onep > 0) {
    impl_onep_s = pin->GetOrAddReal("rad_m1","implicit_one_pass_safety",3.0);
    if (!(impl_onep_s >= 1.0)) {
      ImplFatal("<rad_m1>/implicit_one_pass_safety must be >= 1");
    }
    // implicit_one_pass_auto (m1-onepass-auto): see ImplicitSolve.  A kind of solve that
    // goes auto_window check periods (auto_window * implicit_one_pass eligible solves)
    // without one solve accepted after one pass is switched to one_pass = 0; after
    // auto_reprobe such solves it is switched on again for one check period.  Default
    // true; false on a restart whose file does not carry the key (written before it
    // existed), which then continues as before (see above).
    impl_onep_auto = pin->GetOrAddBoolean("rad_m1","implicit_one_pass_auto",
                                          !global_variable::restart_run);
    impl_onep_awin = pin->GetOrAddInteger("rad_m1","implicit_one_pass_auto_window",2);
    impl_onep_arep = pin->GetOrAddInteger("rad_m1","implicit_one_pass_auto_reprobe",64);
    if (impl_onep_awin < 1 || impl_onep_arep < 1) {
      ImplFatal("<rad_m1>/implicit_one_pass_auto_window and _auto_reprobe must be >= 1");
    }
  }
  impl_pord = pin->GetOrAddInteger("rad_m1","implicit_predictor_order",
                                   (ldef && impl_pred) ? 2 : 1);
  if (impl_pord != 1 && impl_pord != 2) {
    ImplFatal("<rad_m1>/implicit_predictor_order must be 1 or 2");
  }
  impl_odskip = impl_fastk && eddington;
  if (impl_vfold && !impl_stencil) {
    ImplFatal("<rad_m1>/implicit_vimp_fold needs implicit_op_stencil = true");
  }
  }
  // LIMIT 4 of the 3a findings is NOT implemented in 3a2: a column still has to live
  // inside one MeshBlock along x1 (the fatal below).  The option is parsed so that the
  // input files and the gate scripts can already name it, and `gather` fatals rather
  // than silently doing something else.
  std::string spt = pin->GetOrAddString("rad_m1","implicit_partition","none");
  if (spt.compare("none") == 0) {
    impl_part = M1_IPART_NONE;
  } else if (spt.compare("gather") == 0) {
    impl_part = M1_IPART_GATHER;
  } else {
    ImplFatal("<rad_m1>/implicit_partition = '" + spt
              + "' is not a choice (none | gather)");
  }
  impl_recon_w = pin->GetOrAddReal("rad_m1","implicit_recon_w",-1.0);
  impl_res_floor = pin->GetOrAddReal("rad_m1","implicit_res_floor",0.0);
  // implicit_gas_newton_switch = W (he_estall_1003).  DEFAULT 10 on fresh runs (R2
  // default, user 10-03; it was off by default in 69e8f149 only because the switch
  // with the mode-0 root find diverged the BSG envelope, bsg_1001/AB_SOLVER_1003.md,
  // rst 00008: W = 10 died at cycle 14575, resid 0.31; with implicit_tsolve_opac_mode
  // 1 it is robust there, TSOLVE_FIX_1003.md).  A restart whose file lacks the key
  // keeps the old behaviour (0 = off, bitwise the pre-key code) and says so; a named
  // value always wins.  Under implicit_gas_newton, a solve
  // whose Picard residual is detected STALLED drops the gas Newton update (and the
  // opac_newton face terms, which need it) for the rest of that solve: the gas T is the
  // bracketed root of the exact backward-Euler gas equation with kappa_P, kappa_E at the
  // trial T (the implicit_tsolve_opac root find, forced on) in every remaining pass.
  // Stall test at pass it (resid = the global MPI_MAX value, so every rank agrees):
  //   it + 1 >= implicit_gas_newton_switch_min (default 20), it + 1 > W,
  //   min(resid over the last W passes) >= implicit_stall_fac (default 0.5) x the
  //   minimum before that window, and resid rose at least twice inside the window
  //   (a slowly but monotonically contracting solve keeps Newton).
  // Solves that converge never switch.  Diagnosis (he_estall_1003): in hot (2-3 MK),
  // low-density void cells the Newton row omits d kappa_P/dT, the lagged-opacity
  // iteration has slope < -1 and the gas T flips between two values every pass (a
  // resid 1e-8..2e-6 2-cycle that ran ~29 % of the He presn wedge solves to
  // implicit_maxit); the root find with kappa(T) has no such cycle.  He wedge -30..-36 %
  // wall, accuracy within the noise spread; He box bitwise (never switches).
  if (pin->DoesParameterExist("rad_m1","implicit_gas_newton_switch")) {
    impl_gn_sw = pin->GetInteger("rad_m1","implicit_gas_newton_switch");
  } else if (global_variable::restart_run) {
    impl_gn_sw = 0;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> restart input has no implicit_gas_newton_switch: keeping "
                << "the old default (off); set it to switch" << std::endl;
    }
  } else {
    impl_gn_sw = pin->GetOrAddInteger("rad_m1","implicit_gas_newton_switch",10);
  }
  if (impl_gn_sw > 0) {
    impl_gn_sw_min = pin->GetOrAddInteger("rad_m1","implicit_gas_newton_switch_min",20);
    impl_stall_fac = pin->GetOrAddReal("rad_m1","implicit_stall_fac",0.5);
  }
  if (impl_gn_sw < 0 || impl_gn_sw > 64 || (impl_gn_sw > 0 && impl_gn_sw < 3) ||
      impl_gn_sw_min < 1 || !(impl_stall_fac > 0.0 && impl_stall_fac <= 1.0)) {
    ImplFatal("<rad_m1>/implicit_gas_newton_switch must be 0 or lie in [3,64], "
              "implicit_gas_newton_switch_min >= 1 and implicit_stall_fac in (0,1]");
  }
  // sp-blend-1008: on the multi-D wedge berthon | blend default to `step`.  The lagged
  // reduced flux of a cell is the MEAN of its two faces' upwind ratios, and in free
  // streaming the upwind berthon flux reproduces any lagged f (F_f = c f E_up), so
  // recomputing the coefficients every Picard pass only AVERAGES f between neighbours,
  // pass after pass: a neutral iteration in which the zero-flux inner face and the
  // Marshak face leak inward (measured, gate b at CFL 0.4: Picard never meets 1e-12
  // in 200 passes and f falls from 1 to 0.69 in one step).  With `step` the face
  // coefficients come from the start-of-step state and the row is linear in E.
  const bool splag = full && (sph_geom || impl_flux_all) &&
                     (impl_flux != M1_IFLUX_CENTRAL);
  std::string slg = pin->GetOrAddString("rad_m1","implicit_recon_lag",
                                        splag ? "step" : "picard");
  impl_recon_freeze = (slg.compare("step") == 0);
  // MILESTONE 3c: freeze the deferred correction, and with it the plm limiter's choice,
  // after this many Picard passes.  The 3a2 finding is that the limiter keeps switching
  // on a handful of cells and the iteration is a small limit cycle that never meets
  // implicit_tol, so plm_dc costs implicit_maxit passes per step; freezing the
  // correction after a few passes leaves an ordinary linear system to converge.
  // <= 0 (the default) never freezes, i.e. reproduces 3a2.
  impl_recon_npass = pin->GetOrAddInteger("rad_m1","implicit_recon_npass",-1);
  if (pin->DoesParameterExist("rad_m1","implicit_recon_dgpass")) {
    impl_recon_dgpass = pin->GetBoolean("rad_m1","implicit_recon_dgpass");
  }
  if (!impl_recon_freeze && slg.compare("picard") != 0) {
    ImplFatal("<rad_m1>/implicit_recon_lag = '" + slg
              + "' is not a choice (step | picard)");
  }
  // LIMIT 3 of the 3a findings: the DEFAULT is now `true`.  It is a bug fix, not a
  // tuning knob (the boundary face used to hand its whole momentum to one interior cell,
  // which gave that cell 1.5 face-shares of radiative force: He column bottom-cell |v1|
  // 10.3 -> 0.47 v_MLT).  The key is kept so that `false` reproduces runs_3a/RESULTS.txt.
  impl_bmom_half = pin->GetOrAddBoolean("rad_m1","implicit_bmom_half",true);
  if (impl_lin_maxit < 1) {
    ImplFatal("<rad_m1>/implicit_lin_maxit must be >= 1");
  }
  if (!(impl_tol > 0.0) || impl_maxit < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/implicit_tol must be positive and implicit_maxit >= 1"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &mbcs = pmy_pack->pmesh->mesh_bcs;
  ibc_x1min = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      mbcs[static_cast<int>(BoundaryFace::inner_x1)]);
  ibc_x1max = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1max","auto"),
      mbcs[static_cast<int>(BoundaryFace::outer_x1)]);
  iflux_x1min = pin->GetOrAddReal("rad_m1","implicit_flux_x1min",0.0);
  iflux_x1max = pin->GetOrAddReal("rad_m1","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_m1","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_m1","implicit_ebath_x1max",0.0);
  // runs_4f_drift: without it the Marshak end face carries only the comoving flux
  // c q (E - E_bath), so the lab-frame radiation energy advected with the gas,
  // A E = v (1 + chi) E, cannot leave (or enter) through the end: the end cell's E is
  // raised by ~(A/(c q)) E at an outflow end (lowered at an inflow end), and the
  // Lowrie-Edwards shocks drift / carry an N-independent T error.  DEFAULT true since
  // m1-bcadv (tests_m1/runs_4h_bcadv).  On a restart whose file lacks the key (written
  // before 5c9e4432) the old value false is kept; files since then echo their value.
  impl_bc_advect = pin->GetOrAddBoolean("rad_m1","implicit_bc_advect",
                                        !global_variable::restart_run);
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- the restrictions of 3a
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mindcs = pmy_pack->pmesh->mesh_indcs;
  part_nblk = 1;
  if (mindcs.nx1 != indcs.nx1) {
    if (impl_part != M1_IPART_GATHER) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/transport = implicit_x1 with implicit_partition = none "
        << "needs exactly ONE MeshBlock along x1: <meshblock>/nx1 must equal <mesh>/nx1 ("
        << indcs.nx1 << " vs " << mindcs.nx1 << ").  Set <rad_m1>/implicit_partition = "
        << "gather for the line solve partitioned over MeshBlocks and ranks" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((mindcs.nx1 % indcs.nx1) != 0) {
      ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nx1 to be an exact "
                "multiple of <meshblock>/nx1 (uniform mesh only)");
    }
    part_nblk = mindcs.nx1/indcs.nx1;
  }
  if (part_nblk > 1 && (ibc_x1min == M1_IBC_PERIODIC)) {
    ImplFatal("<rad_m1>/implicit_partition = gather does not support PERIODIC x1 across "
              "more than one MeshBlock (the cyclic Thomas sweep of 3a wraps inside one "
              "block).  Use one MeshBlock along x1, or a non-periodic x1 boundary pair");
  }
  if (part_nblk > 1 && indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nghost >= 2");
  }
  if (full) {
    trans_on = pmy_pack->pmesh->multi_d;
    trans_x3 = pmy_pack->pmesh->three_d;
    if (!trans_on) {
      // a 1-D mesh has no transverse direction: the solve IS the x1 column solve, and
      // gate G3 is exactly this statement.
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1>: transport=implicit on a 1-D mesh is the x1 column solve"
                  << std::endl;
      }
    }
  }
  // blendall-1009: which faces carry the per-face beam-vector coefficients
  blat_on = full && trans_on && impl_flux_all;
  bvec_x1 = full && impl_flux_all && !sph_geom;
  if (impl_flux_all && full && (impl_recon == M1_IRECON_PLMDC)) {
    ImplFatal("<rad_m1>/implicit_flux_faces = all takes implicit_recon = dc only");
  }
  if (impl_beam_hr && (part_nblk > 1 || cs_geom || impl_recon == M1_IRECON_PLMDC)) {
    ImplFatal("<rad_m1>/implicit_flux_beam = halfrange needs one MeshBlock along x1, no "
              "cubed sphere and implicit_recon = dc");
  }
  if (impl_beam_fs && !(blat_on && (vet_sc || (vgd_on && sph_geom)))) {
    ImplFatal("<rad_m1>/implicit_flux_beam = fs needs implicit_flux_faces = all on a "
              "multi-D mesh and closure = vet_sc, or vet_gd on the sp wedge");
  }
  if (blat_on && part_nblk > 1) {
    ImplFatal("<rad_m1>/implicit_flux_faces = all needs one MeshBlock along x1");
  }
  if ((indcs.nx2 > 1 || indcs.nx3 > 1) && !impl_allow_multid && !full) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does NOT transport along x2/x3. "
      << "Set <rad_m1>/implicit_allow_multid = true to run a set of INDEPENDENT x1 "
      << "columns (F_2 = F_3 = 0 everywhere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does not support SMR/AMR"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // one solve per hydro step: no sub-cycling, one stage
  subcycle = false;
  nstage = 1;
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (impl_cfl > 0.0) || (!other_sets_dt);

  // ---- arrays
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(f0x1, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1, 0.0);
  Kokkos::realloc(f0x1n, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1n, 0.0);
  // MILESTONE 3b phase C: the BiCGStab wrapper needs the four/six transverse off-diagonal
  // coefficients and ten Krylov vectors on top of the phase-B work array.  They are
  // allocated only when it is selected AND the mesh is multi-D, so line_jacobi keeps the
  // footprint (and, on a 1-D mesh, the solve IS the column solve and there is no system
  // to wrap).
  bicg_on = trans_on && (impl_solver == M1_ISOLV_BICGSTAB);
  // MILESTONE 3b phase D.  `auto` = the off-diagonal Eddington terms go INTO the operator
  // wherever there is a Krylov solver to carry them, and stay LAGGED otherwise, which is
  // what every line_jacobi and implicit_x1 configuration did before phase D.
  if (od_auto) {
    impl_offdiag = bicg_on ? M1_OD_OPERATOR : M1_OD_LAGGED;
  }
  if (impl_offdiag == M1_OD_OPERATOR && trans_on && !bicg_on) {
    ImplFatal("<rad_m1>/implicit_offdiag = operator needs implicit_solver = bicgstab: "
              "the off-diagonal Eddington terms are a 9-/19-point coupling and the "
              "x1 line solve of line_jacobi cannot carry them");
  }
  od_now = impl_offdiag;
  if (impl_vimp && !bicg_on) {
    ImplFatal("<rad_m1>/implicit_vimp needs transport = implicit on a multi-D mesh with "
              "implicit_solver = bicgstab: the implicit velocity couples cells two "
              "apart");
  }
  if (impl_vimp && pmy_pack->pmesh->mb_indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_vimp needs <mesh>/nghost >= 2");
  }
  int niw = full ? (bicg_on ? M1_NIW_K : M1_NIW) : M1_NIW_X1;
  // MILESTONE 3g: the five per-cell components of the gas coupling are APPENDED, so
  // every index above keeps the value it had and the array grows only when asked for.
  iw_gas = -1;
  if (impl_gas_newton || impl_eos_cache) {
    iw_gas = niw;
    niw += M1_NIW_GAS;
  }
  iw_vimp = -1;
  if (impl_vimp) {
    iw_vimp = niw;
    niw += M1_NIW_VIMP;
    // time_scheme = hesdirk2 (read in Time2Init, later): the DA components
    if (pin->DoesParameterExist("rad_m1","time_scheme") &&
        pin->GetString("rad_m1","time_scheme").compare("hesdirk2") == 0) {
      niw += M1_NIW_VIMP_T2;
    }
  }
  iw_muscl = -1;
  if (impl_muscl) {
    iw_muscl = niw;
    niw += M1_NIW_MUSCL;
  }
  if (impl_eos_cache) {
    impl_nec = M1EosCacheNComp(impl_ecnt);
    Kokkos::realloc(ecache, nmb, impl_nec, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ecache, -1.0);
  }
  Kokkos::realloc(iw, nmb, niw, ncells3, ncells2, ncells1);
  if (impl_thin_frz > 0.0) {
    Kokkos::realloc(thin_frz, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thin_frz, 0.0);
  }
  if (impl_opac_newton) {
    if (!impl_opac_update || !impl_gas_newton) {
      ImplFatal("<rad_m1>/implicit_opac_newton needs implicit_opac_update = true and "
                "implicit_gas_newton = true");
    }
    Kokkos::realloc(ktd, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ktd, 0.0);
    Kokkos::realloc(opn_nskip_d, 1);
    Kokkos::deep_copy(opn_nskip_d, 0.0);
  }
  if (impl_face_ktn) {
    if (part_nblk > 1) {
      ImplFatal("<rad_m1>/implicit_face_opac_n needs ONE MeshBlock along x1 (the T^n "
                "snapshot has no x1 halo)");
    }
    Kokkos::realloc(ktn, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ktn, 0.0);
  }
  if (impl_ctrelax > 0.0) {
    if (!trans_on || !impl_clag_step) {
      ImplFatal("<rad_m1>/implicit_closure_thin_relax needs transport = implicit on a "
                "multi-D mesh and implicit_closure_lag = step");
    }
    Kokkos::realloc(ctr_mem, nmb, 4, ncells3, ncells2, ncells1);
  }
  Kokkos::deep_copy(iw, 0.0);
  if (impl_pred) {
    Kokkos::realloc(ipred, nmb, (impl_pord == 2) ? 5 : 3, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ipred, 0.0);
    pred_ok = false;
  }
  Kokkos::realloc(ifw, nmb, M1_NIFW, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(ifw, 0.0);
  if (trans_on) {
    Kokkos::realloc(f0x2, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2, 0.0);
    Kokkos::realloc(f0x2n, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2n, 0.0);
    if (trans_x3) {
      Kokkos::realloc(f0x3, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3, 0.0);
      Kokkos::realloc(f0x3n, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3n, 0.0);
    }
    if (impl_tlim != M1_TLIM_NONE || blat_on) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
      }
    }
    // hrup-1009: the outer Marshak q of a half-range solve with ray data = h+_1 of the
    // top cell (the vacuum sends nothing back), per column, in vcol_q
    hr_q = impl_beam_hr && (vet_sc || (vgd_on && sph_geom)) &&
           (ibc_x1max == M1_IBC_MARSHAK);
    if (hr_q && vcol_q.extent_int(0) < nmb) {
      Kokkos::realloc(vcol_q, nmb, ncells3, ncells2);
      Kokkos::deep_copy(vcol_q, marshak_q);
    }
    if (blat_on) {
      Kokkos::realloc(ifw2, nmb, 3, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(ifw2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(ifw3, nmb, 3, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(ifw3, 0.0);
      }
    }
    if (impl_tlim != M1_TLIM_NONE) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      Kokkos::realloc(klx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(klx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
        Kokkos::realloc(klx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(klx3, 0.0);
      }
    }
    // the transverse halo rides the module's ORDINARY cell-centred boundary machinery on
    // a scratch array, which is what gives it periodic wrap, corner/edge neighbours and
    // MPI for free; the hand-rolled x1 halo of 3b is then not used at all under
    // transport = implicit (this exchange carries its six quantities and seven more).
    Kokkos::realloc(thw, nmb, M1_NHALO_T, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thw, 0.0);
    Kokkos::realloc(thw_c, nmb, M1_NHALO_T, 1, 1, 1);
    pbval_th = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_th->InitializeBuffers(M1_NHALO_T);
    // COMPONENT ROLES of the scratch halos (SetVectorPairs).  Slots 2,3 hold F1 and KT,
    // SCALARS under the polar flip and the seam transform; the real tangential pairs
    // are (N2,N3), (A2,A3), (V2,V3).  The default (IVY, IVZ) rule would get all of that
    // wrong (a latent trap while rad_m1 is a fatal on sp/cs, see rad_m1.cpp).
    {
      auto slot = [](const int iw) {
        int s = -1;
        for (int n=0; n<M1_NHALO_T; ++n) {
          if (M1HaloCompT(n) == iw) {s = n;}
        }
        return s;
      };
      pbval_th->SetVectorPairs(M1_NHALO_T, {{slot(M1_IW_N2), slot(M1_IW_N3)},
                                            {slot(M1_IW_A2), slot(M1_IW_A3)},
                                            {slot(M1_IW_V2), slot(M1_IW_V3)}});
    }
    // the NARROW exchange of the same list: the M1_NHALO_Q components a Picard pass can
    // move once the closure is frozen.  A separate array because the exchange takes the
    // variable count from the array's second extent, and a prefix subview of a
    // LayoutRight array is not one.
    Kokkos::realloc(thq, nmb, M1_NHALO_Q, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thq, 0.0);
    Kokkos::realloc(thq_c, nmb, M1_NHALO_Q, 1, 1, 1);
    pbval_tq = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_tq->InitializeBuffers(M1_NHALO_Q);
    pbval_tq->SetVectorPairs(M1_NHALO_Q, {});   // E', G0, F1, KT: all scalars
    // the deep interior of the scratch arrays is neither read by a send nor written by a
    // receive when every neighbour is at the SAME level (a same-level buffer reaches ng
    // cells in from the active boundary, buffs_cc.cpp) and neither the cubed-sphere
    // resample nor the polar transform is in play; then the copies to and from iw can
    // skip it.  SMR/AMR is already a fatal above.
    halo_shell = !(pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere ||
                   pmy_pack->pmesh->use_polar_boundary);
    {
      // ONE more exchange object, for the single Krylov vector the operator application
      // needs in its ghost zones -- and for the one-component exchange of E alone, which
      // is why it is allocated under every solver.  It is used strictly SEQUENTIALLY
      // with pbval_th and pbval_tq (each exchange runs its
      // InitRecv/Send/Recv/Clear chain to completion before the next
      // starts) and every rank issues the identical SEQUENCE of exchanges -- the Picard
      // count, the BiCGStab count and every breakdown decision are taken from GLOBAL
      // reductions -- so MPI's non-overtaking guarantee keeps the two streams apart even
      // though they share the tag space.
      Kokkos::realloc(krw, nmb, 1, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(krw, 0.0);
      Kokkos::realloc(krw_c, nmb, 1, 1, 1, 1);
      pbval_kr = new MeshBoundaryValuesCC(pmy_pack, pin, false);
      pbval_kr->InitializeBuffers(1);
      pbval_kr->SetVectorPairs(1, {});   // one scalar Krylov vector
      // implicit_vimp: the per-pass exchange of the Jacobian rows and dv^k (same
      // sequential-use argument as above)
      if (impl_vimp) {
        Kokkos::realloc(vmw, nmb, M1_NVIMP_X, ncells3, ncells2, ncells1);
        Kokkos::deep_copy(vmw, 0.0);
        Kokkos::realloc(vmw_c, nmb, M1_NVIMP_X, 1, 1, 1);
        pbval_vm = new MeshBoundaryValuesCC(pmy_pack, pin, false);
        pbval_vm->InitializeBuffers(M1_NVIMP_X);
        // (DV2, DV3) is the tangential pair; the nine per-direction Jacobian rows P are
        // treated as scalars (their proper curvilinear treatment is a later stage)
        pbval_vm->SetVectorPairs(M1_NVIMP_X, {{M1_IV_DV + 1, M1_IV_DV + 2}});
      }
    }
    if (impl_halo_direct) {ImplicitHaloDirectInit();}
    if (impl_odc) {
      Kokkos::realloc(odc, nmb, 3, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(odc, 0.0);
    }
    if (impl_stencil) {
      if (ibc_x1min == M1_IBC_PERIODIC) {
        ImplFatal("<rad_m1>/implicit_op_stencil does not take a periodic x1 wrap");
      }
      Kokkos::realloc(ost, nmb, impl_vfold ? 25 : 19, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(ost, 0.0);
    }
  }
  // MILESTONE 3e: the Anderson histories.  Allocated ONLY when the acceleration is on.
  if (impl_accel != M1_IACC_NONE) {
    aa_nc = trans_on ? 4 : 2;
    Kokkos::realloc(aa_xc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xc, 0.0);
    Kokkos::realloc(aa_fc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fc, 0.0);
    Kokkos::realloc(aa_xp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xp, 0.0);
    Kokkos::realloc(aa_fp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fp, 0.0);
    Kokkos::realloc(aa_dx, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_dx, 0.0);
    Kokkos::realloc(aa_df, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_df, 0.0);
    Kokkos::realloc(aa_sc, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_sc, 1.0);
  }
  ImplicitPartitionInit();

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: transport=" << (full ? "implicit" : "implicit_x1")
              << (full ? (trans_x3 ? " (3-D)" : (trans_on ? " (2-D)" : " (1-D)")) : "")
              << " solver=" << (full ? (bicg_on ? "bicgstab" : "line_jacobi") : "thomas")
              << " lin_tol=" << impl_lin_tol
              << " lin_maxit=" << impl_lin_maxit << std::endl;
    std::cout << "         implicit_cfl="
              << impl_cfl << " tol=" << impl_tol << " maxit=" << impl_maxit
              << " opac_update=" << (impl_opac_update ? "true" : "false")
              << " marshak_q=" << marshak_q << std::endl;
    std::cout << "         implicit_flux="
              << ((impl_flux == M1_IFLUX_APHLL) ? "ap_hll" :
                  ((impl_flux == M1_IFLUX_BERTHON) ? "berthon" :
                   ((impl_flux == M1_IFLUX_BLEND) ? "blend" : "central")))
              << ((impl_flux == M1_IFLUX_BLEND) ? (":" + sbl + "/" + sbm + "/" + sbw)
                                                : std::string(""))
              << " implicit_recon="
              << ((impl_recon == M1_IRECON_PLMDC) ? "plm_dc" : "dc")
              << " implicit_partition="
              << ((impl_part == M1_IPART_GATHER) ? "gather" : "none") << std::endl;
    if (impl_vimp) {
      std::cout << "         implicit_vimp=true (v' of the enthalpy flux implicit)"
                << std::endl;
    }
    if (impl_enth != M1_IENTH_UPWIND) {
      std::cout << "         implicit_enthalpy="
                << ((impl_enth == M1_IENTH_PLM) ? "plm" : "central")
                << " (deferred correction)" << std::endl;
    }
    if (ibc_x1min == M1_IBC_MARSHAK || ibc_x1max == M1_IBC_MARSHAK) {
      std::cout << "         implicit_bc_advect=" << (impl_bc_advect ? "true" : "false")
                << " (Marshak x1 end carries A E)" << std::endl;
    }
    if (trans_on) {
      std::cout << "         implicit_offdiag="
                << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                    ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
                << (od_auto ? " (auto)" : "")
                << " closure_relax=" << impl_crelax
                << (impl_crelax_thin ? " (thin faces only)" : "")
                << " closure_lag=" << (impl_clag_step ? "step" : "pass")
                << " trans_limit="
                << ((impl_tlim == M1_TLIM_LP) ? "lp" : "none")
                << ((impl_tlim == M1_TLIM_LP)
                    ? (" fmax=" + std::to_string(impl_tfmax)) : std::string(""))
                << std::endl;
    }
    std::cout << "         implicit_accel="
              << ((impl_accel == M1_IACC_ANDERSON) ? "anderson" : "none")
              << ((impl_accel == M1_IACC_ANDERSON)
                  ? (" m=" + std::to_string(impl_and_m)
                     + " beta=" + std::to_string(impl_and_beta)
                     + " start=" + std::to_string(impl_and_start)
                     + " ncomp=" + std::to_string(aa_nc))
                  : std::string("")) << std::endl;
    if (impl_gas_newton || impl_eos_cache) {
      std::cout << "         implicit_gas_newton="
                << (impl_gas_newton ? "true" : "false")
                << " implicit_eos_cache=" << (impl_eos_cache ? "true" : "false")
                << (impl_eos_cache ? (" nt=" + std::to_string(impl_ecnt)
                                      + " ncomp=" + std::to_string(impl_nec)
                                      + " check="
                                      + (impl_eccheck ? "true" : "false"))
                                   : std::string("")) << std::endl;
    }
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {VetInit(pin);}
  Time2Init(pin);
  MRInit(pin);         // implicit_mr_every (rad_m1_mr.cpp)
  // accel-1009 (memory): u1 is used only by the explicit stage copy, hesdirk2 and the
  // multi-rate step; a backward-Euler implicit run without multi-rate frees it
  if (transport >= M1_TRANSPORT_IMPLICIT_X1 && time_scheme == M1_TIME_BE &&
      mr_every <= 1) {
    Kokkos::realloc(u1, u1.extent(0), M1_NVAR, 1, 1, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPartitionInit
//! \brief milestone 3b, LIMIT 4: build the topology of the x1 stacks and allocate the
//! gather/scatter and halo buffers.  A "stack" is the set of MeshBlocks that share the
//! same (x2,x3) footprint, ordered by their x1 logical location; its ROOT is the block
//! with the lowest one.  Uniform mesh only (the caller has already fatalled on SMR/AMR
//! and on a non-integer block count along x1), so the stack of a block is found by a
//! scan of the global LogicalLocation list.
//!
//! Nothing is allocated and nothing is communicated when part_nblk == 1: the 3a path is
//! then bitwise what it was.

void RadiationM1::ImplicitPartitionInit() {
  part_nroot = 0;
  part_nx1g = 0;
  part_nlay = 0;
  part_nqa = M1_NHALO_A;
  part_any_mpi = false;

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int g0 = pmy_pack->gids;
  // the three per-block tables are ALWAYS allocated: ImplicitSolve reads part_pos in
  // every kernel, and with one MeshBlock per column the defaults (position 0 of a stack
  // of one) select exactly the 3a branches.
  Kokkos::realloc(part_pos, nmb);
  Kokkos::realloc(part_slot, nmb);
  Kokkos::realloc(part_nbr, 2*nmb);
  for (int m=0; m<nmb; ++m) {
    part_pos.h_view(m) = 0;
    part_slot.h_view(m) = -1;
    part_nbr.h_view(2*m) = -1;
    part_nbr.h_view(2*m+1) = -1;
  }
  if (part_nblk <= 1) {
    part_pos.modify_host();
    part_pos.sync_device();
    part_slot.modify_host();
    part_slot.sync_device();
    part_nbr.modify_host();
    part_nbr.sync_device();
    return;
  }
  part_nx1g = part_nblk*indcs.nx1;
  part_nlay = std::min(indcs.ng, 2);

  // (1) for every LOCAL block: its position in the stack, the global ids of the stack
  // members and of its two x1 neighbours.
  part_rootgid.assign(nmb, -1);
  part_rootrank.assign(nmb, -1);
  part_nbrrank.assign(2*nmb, -1);
  part_nbrgid.assign(2*nmb, -1);
  std::vector<int> stack(part_nblk);
  for (int m=0; m<nmb; ++m) {
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {stack[p] = -1;}
    for (int g=0; g<pm->nmb_total; ++g) {
      LogicalLocation &lg = pm->lloc_eachmb[g];
      if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level) {
        if (lg.lx1 >= 0 && lg.lx1 < part_nblk) {stack[lg.lx1] = g;}
      }
    }
    for (int p=0; p<part_nblk; ++p) {
      if (stack[p] < 0) {
        ImplFatal("<rad_m1>/implicit_partition = gather: the x1 stack of a MeshBlock is "
                  "incomplete (a non-uniform mesh?)");
      }
    }
    int pos = static_cast<int>(lm.lx1);
    part_pos.h_view(m) = pos;
    part_rootgid[m] = stack[0];
    part_rootrank[m] = pm->rank_eachmb[stack[0]];
    if (pos > 0) {
      part_nbrgid[2*m] = stack[pos-1];
      part_nbrrank[2*m] = pm->rank_eachmb[stack[pos-1]];
    }
    if (pos < part_nblk-1) {
      part_nbrgid[2*m+1] = stack[pos+1];
      part_nbrrank[2*m+1] = pm->rank_eachmb[stack[pos+1]];
    }
    for (int s=0; s<2; ++s) {
      int gn = part_nbrgid[2*m+s];
      bool loc = (gn >= 0) && (part_nbrrank[2*m+s] == global_variable::my_rank);
      part_nbr.h_view(2*m+s) = loc ? (gn - g0) : -1;
      if (gn >= 0 && part_nbrrank[2*m+s] != global_variable::my_rank) {
        part_any_mpi = true;
      }
    }
    if (part_rootrank[m] != global_variable::my_rank) {part_any_mpi = true;}
  }

  // (2) the local ROOTS and their member lists
  part_mrank.clear();
  part_mgid.clear();
  std::vector<int> rootmb;
  for (int m=0; m<nmb; ++m) {
    if (part_pos.h_view(m) != 0) continue;
    rootmb.push_back(m);
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {
      int gfound = -1;
      for (int g=0; g<pm->nmb_total; ++g) {
        LogicalLocation &lg = pm->lloc_eachmb[g];
        if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level && lg.lx1 == p) {
          gfound = g;
        }
      }
      part_mgid.push_back(gfound);
      part_mrank.push_back(pm->rank_eachmb[gfound]);
      if (pm->rank_eachmb[gfound] != global_variable::my_rank) {part_any_mpi = true;}
    }
  }
  part_nroot = static_cast<int>(rootmb.size());
  for (int m=0; m<nmb; ++m) {
    int sl = -1;
    if (part_rootrank[m] == global_variable::my_rank) {
      for (int s=0; s<part_nroot; ++s) {
        if (rootmb[s] + g0 == part_rootgid[m]) {sl = s;}
      }
    }
    part_slot.h_view(m) = sl;
  }
  part_pos.modify_host();
  part_pos.sync_device();
  part_slot.modify_host();
  part_slot.sync_device();
  part_nbr.modify_host();
  part_nbr.sync_device();

  // (3) buffers
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(part_sys, std::max(part_nroot,1), 6, ncells3, ncells2, part_nx1g);
  Kokkos::deep_copy(part_sys, 0.0);
  int nrow = 4*ncells3*ncells2*indcs.nx1;
  Kokkos::realloc(part_sbuf, nmb, nrow);
  Kokkos::realloc(part_rbuf, std::max(part_nroot*part_nblk,1), nrow);
  Kokkos::realloc(part_sbuf_h, nmb, nrow);
  Kokkos::realloc(part_rbuf_h, std::max(part_nroot*part_nblk,1), nrow);
  int nhal = M1_NHALO_A*part_nlay*ncells3*ncells2;
  Kokkos::realloc(part_hbuf, nmb, 4, nhal);
  Kokkos::realloc(part_hbuf_h, nmb, 4, nhal);

  if (global_variable::my_rank == 0) {
    std::cout << "         implicit_partition=gather: " << part_nblk
              << " MeshBlocks per x1 column, " << part_nx1g << " rows per gathered line, "
              << part_nlay << " halo layers" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitX1Halo
//! \brief exchange the x1 ghost layers of the work array `iw` between the MeshBlocks of
//! one x1 stack.  `eponly` picks the one-quantity set (M1_IW_EP, the new iterate, needed
//! by the face update and by the next pass) instead of the six LAGGED quantities
//! (M1HaloCompA).
//!
//! What makes the partitioned solve BITWISE identical to a single-block solve is that
//! the ghost value a block reads is the VERY NUMBER its neighbour computed, copied, and
//! never a quantity recomputed from a hydro ghost: this routine moves nothing else.

void RadiationM1::ImplicitX1Halo(bool eponly) {
  if (part_nblk <= 1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nl = part_nlay;
  int nq = eponly ? 1 : M1_NHALO_A;
  const bool ep1 = eponly;
  auto iw_ = iw;
  auto nbr_ = part_nbr;
  int nj = je - js + 1, nk = ke - ks + 1;

  // (1) same-rank neighbours: a plain device copy, neighbour ACTIVE -> my GHOST
  par_for("m1_impl_halo_loc", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int mlo = nbr_.d_view(2*m);
    int mhi = nbr_.d_view(2*m+1);
    if (mlo >= 0) {iw_(m,nc,k,j,is-1-l) = iw_(mlo,nc,k,j,ie-l);}
    if (mhi >= 0) {iw_(m,nc,k,j,ie+1+l) = iw_(mhi,nc,k,j,is+l);}
  });
  if (!part_any_mpi) return;

#if MPI_PARALLEL_ENABLED
  // (2) remote neighbours.  Buffer slots: 0 = send to lo, 1 = send to hi,
  // 2 = recv from lo, 3 = recv from hi.
  auto hb_ = part_hbuf;
  par_for("m1_impl_halo_pack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    hb_(m,0,idx) = iw_(m,nc,k,j,is+l);
    hb_(m,1,idx) = iw_(m,nc,k,j,ie-l);
  });
  int nbuf = nq*nl*nk*nj;
  Kokkos::deep_copy(part_hbuf_h, part_hbuf);
  std::vector<MPI_Request> req;
  // the tag identifies (receiving local block, receiving side, which halo set); the
  // source rank is named in the receive, so it need only be unique per rank pair.
  int te = eponly ? 1 : 0;
  int *gr = pmy_pack->pmesh->gids_eachrank;
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_hbuf_h(m,2+s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*m + 2*s + te, MPI_COMM_WORLD, &req.back());
    }
  }
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      // the neighbour receives this message into ITS slot 2+(1-s)
      int lidn = part_nbrgid[2*m+s] - gr[rk];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_hbuf_h(m,s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*lidn + 2*(1-s) + te, MPI_COMM_WORLD, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  Kokkos::deep_copy(part_hbuf, part_hbuf_h);
  auto nbrk = part_nbr;
  par_for("m1_impl_halo_unpack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    if (nbrk.d_view(2*m) < 0) {iw_(m,nc,k,j,is-1-l) = hb_(m,2,idx);}
    if (nbrk.d_view(2*m+1) < 0) {iw_(m,nc,k,j,ie+1+l) = hb_(m,3,idx);}
  });
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransverseHalo
//! \brief milestone 3b phase B: put the M1_NHALO_T LAGGED quantities the x2/x3 faces (and
//! the lagged off-diagonal Eddington terms of every face) read at a neighbouring cell
//! into the scratch array `thw`, exchange it with ALL six neighbours through the module's
//! ORDINARY cell-centred boundary machinery -- which is what supplies periodic wrap,
//! edge/corner neighbours and MPI without another hand-rolled protocol -- and copy the
//! ghost zones back into `iw`.
//!
//! The values a block reads in its ghost zones are therefore the VERY NUMBERS its
//! neighbour computed, so the two blocks that share a face assemble that face's flux from
//! bit-identical inputs: the face flux is single-valued and the scheme is conservative
//! across a MeshBlock boundary (gate G4).
//!
//! PHYSICAL (non-periodic) boundaries are NOT filled here: every kernel below branches on
//! the boundary flag instead and imposes F = 0 on a physical x2/x3 face (reflecting).

void RadiationM1::ImplicitTransverseHalo(int nq) {
  if (!trans_on) return;
  ImplicitHaloExchange(nq, -1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloCopy
//! \brief copy `nq` components between iw and the scratch halo array `sc` (`topack`
//! picks the direction), over the HALO SHELL only.
//!
//! What the exchange reads out of `sc` is the outermost ng ACTIVE cells, and what it
//! writes back into it is the ghost zones (buffs_cc.cpp, isame); with same-level
//! neighbours only -- which is all the implicit solve allows -- nothing else in `sc` is
//! ever touched.  The box [is+ng,ie-ng] x [js+ng,je-ng] x [ks+ng,ke-ng] is therefore
//! neither sent nor received and need not be copied in either direction: on the way in
//! its value is never read, and on the way back it is a copy of what the pack put there.
//! At a 2-D 84 x 32 block with ng = 2 that is 480 cells per component instead of 3168.
//! `halo_shell` falls back to the full copy for the exchanges that reach deeper (the
//! cubed-sphere resample and the polar transform read the whole strip).
//!
//! The loop carries (m,n,k,j) and runs the contiguous i direction inside, so the index
//! arithmetic of the flattened range policy (an integer division per index) is paid once
//! per ROW instead of once per cell.

void RadiationM1::ImplicitHaloCopy(DvceArray5D<Real> &sc, int nq, int c0, bool topack) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // the deep interior, as an index box.  A direction with no neighbour (a degenerate
  // dimension) is interior everywhere; a direction too thin to have one (nx <= 2 ng), or
  // a run whose exchange reaches deeper, makes the box empty and the copy full.
  int ng = indcs.ng;
  int ilo = halo_shell ? (indcs.is + ng) : n1;
  int ihi = halo_shell ? (indcs.ie - ng) : (n1 - 1);
  int jlo = (indcs.nx2 > 1) ? (indcs.js + ng) : 0;
  int jhi = (indcs.nx2 > 1) ? (indcs.je - ng) : 0;
  int klo = (indcs.nx3 > 1) ? (indcs.ks + ng) : 0;
  int khi = (indcs.nx3 > 1) ? (indcs.ke - ng) : 0;
  if (ilo > ihi) {
    ilo = n1;
    ihi = n1 - 1;
  }
  const int il_ = ilo, iu_ = ihi, jl_ = jlo, ju_ = jhi, kl_ = klo, ku_ = khi;
  const int nc0 = c0, nn1 = n1;
  const bool pack_ = topack;
  auto iw_ = iw;
  auto sc_ = sc;
  if (impl_halo_direct) {
    // implicit_halo_direct on a mesh whose neighbours are not all on this rank: the
    // same shell copy with ONE thread per cell (coalesced along i) instead of one per
    // row -- the row loop runs ~13k threads with strided rows on the 3-D box.
    par_for("m1_impl_hcpyf", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_) && (i >= il_) &&
          (i <= iu_)) {
        return;
      }
      const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
      if (pack_) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      } else {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    });
    return;
  }
  par_for("m1_impl_hcpy", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j) {
    const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
    // the two i runs this row copies: the whole row unless the row is interior
    int ia = 0, ib = nn1 - 1, ic = nn1, id = nn1 - 1;
    if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_)) {
      ib = il_ - 1;
      ic = iu_ + 1;
    }
    if (pack_) {
      for (int i=ia; i<=ib; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
    } else {
      for (int i=ia; i<=ib; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloExchange
//! \brief the common core: pack, exchange through the ordinary cell-centred boundary
//! machinery on the PERSISTENT scratch array whose width matches `nq`, unpack.  The
//! three widths (M1_NHALO_T, M1_NHALO_Q, 1) have one array and one boundary object each;
//! they are used strictly sequentially, so they may share the MPI tag space.

void RadiationM1::ImplicitHaloExchange(int nq, int c0) {
  if (halo_direct_on) {   // implicit_halo_direct, every neighbour on this rank
    ImplicitHaloDirect(nq, c0);
    return;
  }
  if (impl_halo_mpi) {     // implicit_halo_mpi: on-rank copy + one message per rank
    if (hm_state == 0) {ImplicitHaloMPIInit();}
    if (hm_state == 1) {
      ImplicitHaloMPI(nq, c0);
      return;
    }
  }
  DvceArray5D<Real> *pa, *pc;
  MeshBoundaryValuesCC *pb;
  if (nq == M1_NHALO_T) {
    pa = &thw; pc = &thw_c; pb = pbval_th;
  } else if (nq == M1_NHALO_Q) {
    pa = &thq; pc = &thq_c; pb = pbval_tq;
  } else if (nq == M1_NVIMP_X) {
    pa = &vmw; pc = &vmw_c; pb = pbval_vm;
  } else {
    pa = &krw; pc = &krw_c; pb = pbval_kr;
  }
  ImplicitHaloCopy(*pa, nq, c0, true);
  while (pb->InitRecv(nq) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
  ImplicitHaloCopy(*pa, nq, c0, false);
  // STAGE CS2 (C3): (N2, N3) of the work array are FACE-NORMAL components (as A and V)
  // and cross a seam by the signed axis permutation (cs_perm_pairs, CubedS1Init), which
  // keeps |n|; it is renormalised in the ghost cells all the same (to round-off a no-op),
  // so that a later change of the seam map cannot hand the closure a non-unit n.  A zero
  // n (eddington; vet_col about r_hat has n = r_hat) is left alone.
  if (cs_geom && nq == M1_NHALO_T) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int n1 = indcs.nx1 + 2*indcs.ng, n2 = indcs.nx2 + 2*indcs.ng;
    const int n3 = indcs.nx3 + 2*indcs.ng;
    auto iw_ = iw;
    par_for("m1_cs_nrenorm", DevExeSpace(), 0, pmy_pack->nmb_thispack-1, 0, n3-1, 0, n2-1,
            0, n1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) return;
      const Real a = iw_(m,M1_IW_N1,k,j,i), b = iw_(m,M1_IW_N2,k,j,i);
      const Real d = iw_(m,M1_IW_N3,k,j,i);
      const Real nn = a*a + b*b + d*d;
      if (!(nn > 0.0)) return;
      const Real f = 1.0/sqrt(nn);
      iw_(m,M1_IW_N1,k,j,i) = a*f;
      iw_(m,M1_IW_N2,k,j,i) = b*f;
      iw_(m,M1_IW_N3,k,j,i) = d*f;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitLatFaceCoef
//! \brief blendall-1009 (implicit_flux_faces = all): the berthon / blend coefficients
//! AL, HCL, HCR of every x2 (and x3) face, from the lagged closure of the two cells
//! (M1BeamFace with mu = n_2 or n_3, the face optical depth (rho kappa_T)_f dx_f with the
//! arithmetic face mean the transverse theta uses and the centre distance: dx2/dx3 on
//! Cartesian, dxface.x2f/x3f on sp).  Every input (WCHI, N1..N3, KT) is in the
//! transverse halo.  A physical face is reflecting (F = 0): all three are zero there.
//! Called from ImplicitTransverseTerms before ImplicitTransTheta, which multiplies the
//! face theta by 1 - AL.

void RadiationM1::ImplicitLatFaceCoef() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto fw2 = ifw2;
  auto fw3 = ifw3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  auto cdxf = pmy_pack->pcoord->dxface;
  const bool sph = sph_geom;
  const bool edd = eddington;
  const bool blend = (impl_flux == M1_IFLUX_BLEND);
  const int bkind = impl_blend, bfm = impl_blend_fmode;
  const Real btau0 = impl_blend_tau0;
  const Real bflo = impl_blend_flo, bfhi = impl_blend_fhi;
  const Real ch = chat;
  const Real balph = impl_blend_alpha, br0 = impl_blend_r0;
  const Real bx0 = impl_blend_xthin, dtl = dt_sub;
  const int bxmode = impl_blend_xthin_mode;
  const Real bxwmin = impl_blend_xthin_wmin, bxr0 = impl_blend_xthin_r0;
  if (impl_beam_hr) {
    // hrup-1009: per direction d, h+_d -> S1 and h-_d -> S3 of every active cell
    // (M1HrCell), both through the transverse halo, then the faces (M1HrFace)
    const bool hgd = vgd_on && sph && (vgd_hr.extent_int(0) > 0);
    const int hsrc = vet_sc ? 1 : (hgd ? 2 : (impl_hr_model ? 0 : -1));
    auto hvc_ = vet_cell;
    auto hgh_ = vgd_hr;
    for (int d = 1; d <= (trans_x3 ? 2 : 1); ++d) {
      par_for("m1_impl_lfchr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real hp, hm, fb;
        M1HrCell(hsrc, iw_, hvc_, hgh_, m, k, j, i, d, 0.0, hp, hm, fb);
        iw_(m,M1_IW_S1,k,j,i) = hp;
        iw_(m,M1_IW_S3,k,j,i) = hm;
        // beam_kn: the light-front flag (M1HrJFlag) rides in RES as fb + 2 (fb <= 1)
        if (bxmode == 2) {fb += 2.0*M1HrJFlag(hsrc, iw_, hvc_, m, k, j, i);}
        iw_(m,M1_IW_RES,k,j,i) = fb;
      });
      ImplicitHaloExchange(1, M1_IW_S1);
      ImplicitHaloExchange(1, M1_IW_S3);
      ImplicitHaloExchange(1, M1_IW_RES);
      auto fw = (d == 1) ? fw2 : fw3;
      const int kup = (d == 2) ? ke + 1 : ke, jup = (d == 1) ? je + 1 : je;
      const int inx = (d == 1) ? BoundaryFace::inner_x2 : BoundaryFace::inner_x3;
      const int onx = (d == 1) ? BoundaryFace::outer_x2 : BoundaryFace::outer_x3;
      par_for("m1_impl_lfchf", DevExeSpace(), 0, nmb1, ks, kup, js, jup, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        BoundaryFlag blo = mbbcs.d_view(m,inx);
        BoundaryFlag bhi = mbbcs.d_view(m,onx);
        bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
        bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
        const int c = (d == 1) ? j : k;
        const int cs = (d == 1) ? js : ks, ce = (d == 1) ? je : ke;
        if ((c == cs && plo) || (c == ce + 1 && phi)) {
          fw(m,M1_IFW_AL,k,j,i) = 0.0;
          fw(m,M1_IFW_HCL,k,j,i) = 0.0;
          fw(m,M1_IFW_HCR,k,j,i) = 0.0;
          return;
        }
        const int km = (d == 2) ? k - 1 : k, jm = (d == 1) ? j - 1 : j;
        Real dx = (d == 1) ? mbsize.d_view(m).dx2 : mbsize.d_view(m).dx3;
        if (sph) {dx = (d == 1) ? cdxf.x2f(m,k,j,i) : cdxf.x3f(m,k,j,i);}
        const Real tauf = 0.5*(iw_(m,M1_IW_KT,km,jm,i) + iw_(m,M1_IW_KT,k,j,i))*dx;
        Real al, hl, hr;
        Real fbl = iw_(m,M1_IW_RES,km,jm,i), fbr = iw_(m,M1_IW_RES,k,j,i), jfl = 0.0;
        if (fbl > 1.5) {fbl -= 2.0; jfl = 1.0;}
        if (fbr > 1.5) {fbr -= 2.0; jfl = 1.0;}
        M1HrFace(iw_(m,M1_IW_S1,km,jm,i), iw_(m,M1_IW_S3,k,j,i), fbl, fbr, tauf, blend,
                 bkind, bfm, btau0, bflo, bfhi, balph,
                 ch, iw_(m,M1_IW_EP,km,jm,i), iw_(m,M1_IW_EP,k,j,i), br0,
                 iw_(m,M1_IW_S3,km,jm,i), iw_(m,M1_IW_S1,k,j,i), ch*dtl/dx, bx0, bxmode,
                 bxwmin, bxr0, jfl, al, hl, hr);
        fw(m,M1_IFW_AL,k,j,i) = al;
        fw(m,M1_IFW_HCL,k,j,i) = hl;
        fw(m,M1_IFW_HCR,k,j,i) = hr;
      });
    }
    return;
  }
  // implicit_flux_beam = fs: |f| = |H|/J of the vet_sc formal solution, into the Thomas
  // scratch M1_IW_S1 (free at this point of the pass) and through the transverse halo
  const bool bfs = impl_beam_fs;
  // vet_gd (sp): the beam vector of the twin-free vet_gd intensities, J and H in the
  // local (r, theta, phi) basis: |f| = |H|/J -> S1, n_theta -> S3, n_phi -> RES (all
  // three free at this point of the pass), below the vet_gd cut (vlat_icut) or where
  // J <= 0 the lagged closure's f_Lev(chi) and n.  Exchanged through the halo.
  const bool bgd = bfs && vgd_on && sph;
  if (bgd) {
    auto vi_ = vgd_i;
    auto dir_ = vgd_dir;
    const int nd = vgd_n;
    const int og = vgd_w - indcs.ng;
    const int ilo = is + vlat_icut;
    par_for("m1_impl_lfcgd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real fm = M1FLev(iw_(m,M1_IW_WCHI,k,j,i));
      Real n2 = iw_(m,M1_IW_N2,k,j,i), n3 = iw_(m,M1_IW_N3,k,j,i);
      if (i >= ilo) {
        const Real th = mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
        const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        Real jm = 0.0, hr = 0.0, ht = 0.0, hp = 0.0;
        for (int d = 0; d < nd; ++d) {
          const Real wi = dir_(d,3)*vi_(m,d,k+og,j+og,i);
          jm += wi;
          hr += wi*(dir_(d,0)*st*cp + dir_(d,1)*st*sp + dir_(d,2)*ct);
          ht += wi*(dir_(d,0)*ct*cp + dir_(d,1)*ct*sp - dir_(d,2)*st);
          hp += wi*(-dir_(d,0)*sp + dir_(d,1)*cp);
        }
        const Real hm = sqrt(hr*hr + ht*ht + hp*hp);
        if (jm > 0.0 && hm > 0.0) {
          fm = fmin(hm/jm, 1.0);
          n2 = ht/hm;
          n3 = hp/hm;
        }
      }
      iw_(m,M1_IW_S1,k,j,i) = fm;
      iw_(m,M1_IW_S3,k,j,i) = n2;
      iw_(m,M1_IW_RES,k,j,i) = n3;
    });
    ImplicitHaloExchange(1, M1_IW_S1);
    ImplicitHaloExchange(1, M1_IW_S3);
    ImplicitHaloExchange(1, M1_IW_RES);
  } else if (bfs) {
    auto vc_ = vet_cell;
    par_for("m1_impl_lfcfs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jj = vc_(m,M1_VET_J,k,j,i);
      const Real h1 = vc_(m,M1_VET_H1,k,j,i), h2 = vc_(m,M1_VET_H1+1,k,j,i);
      const Real h3 = vc_(m,M1_VET_H1+2,k,j,i);
      iw_(m,M1_IW_S1,k,j,i) = (jj > 0.0) ? fmin(sqrt(h1*h1 + h2*h2 + h3*h3)/jj, 1.0)
                                         : 0.0;
    });
    ImplicitHaloExchange(1, M1_IW_S1);
  }
  const int cn2 = bgd ? M1_IW_S3 : M1_IW_N2;
  const int cn3 = bgd ? M1_IW_RES : M1_IW_N3;
  auto fmag = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i) {
    return bfs ? iw_(m,M1_IW_S1,k,j,i) : M1FLev(iw_(m,M1_IW_WCHI,k,j,i));
  };
  par_for("m1_impl_lfc2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      fw2(m,M1_IFW_AL,k,j,i) = 0.0;
      fw2(m,M1_IFW_HCL,k,j,i) = 0.0;
      fw2(m,M1_IFW_HCR,k,j,i) = 0.0;
      return;
    }
    const int jm = j - 1;
    Real dx = mbsize.d_view(m).dx2;
    if (sph) {dx = cdxf.x2f(m,k,j,i);}
    const Real tauf = 0.5*(iw_(m,M1_IW_KT,k,jm,i) + iw_(m,M1_IW_KT,k,j,i))*dx;
    Real al, hl, hr;
    M1BeamFace(fmag(m,k,jm,i), iw_(m,cn2,k,jm,i),
               fmag(m,k,j,i), iw_(m,cn2,k,j,i), tauf, edd, blend, bkind,
               bfm, btau0, bflo, bfhi, ch, al, hl, hr);
    fw2(m,M1_IFW_AL,k,j,i) = al;
    fw2(m,M1_IFW_HCL,k,j,i) = hl;
    fw2(m,M1_IFW_HCR,k,j,i) = hr;
  });
  if (std::getenv("M1_BLDBG") != nullptr) {   // DEBUG (blendall-1009): face statistics
    auto h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ifw2);
    auto h1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ifw);
    auto hw = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
    Real s[3][4] = {{0.0}};
    int n = 0;
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        for (int c = 0; c < 3; ++c) {
          s[0][c] += h2(0,c,ks,j,i);
          s[1][c] += h1(0,c,ks,j,i);
        }
        s[2][0] += hw(0,M1_IW_WCHI,ks,j,i);
        s[2][1] += fabs(hw(0,M1_IW_N1,ks,j,i));
        s[2][2] += fabs(hw(0,M1_IW_N2,ks,j,i));
        ++n;
      }
    }
    std::cout << "BLDBG x2 <AL,HCL,HCR> " << s[0][0]/n << " " << s[0][1]/n << " "
              << s[0][2]/n << "  x1 " << s[1][0]/n << " " << s[1][1]/n << " "
              << s[1][2]/n << "  <chi,|n1|,|n2|> " << s[2][0]/n << " " << s[2][1]/n
              << " " << s[2][2]/n << std::endl;
  }
  if (!trans_x3) return;
  par_for("m1_impl_lfc3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) {
      fw3(m,M1_IFW_AL,k,j,i) = 0.0;
      fw3(m,M1_IFW_HCL,k,j,i) = 0.0;
      fw3(m,M1_IFW_HCR,k,j,i) = 0.0;
      return;
    }
    const int km = k - 1;
    Real dx = mbsize.d_view(m).dx3;
    if (sph) {dx = cdxf.x3f(m,k,j,i);}
    const Real tauf = 0.5*(iw_(m,M1_IW_KT,km,j,i) + iw_(m,M1_IW_KT,k,j,i))*dx;
    Real al, hl, hr;
    M1BeamFace(fmag(m,km,j,i), iw_(m,cn3,km,j,i),
               fmag(m,k,j,i), iw_(m,cn3,k,j,i), tauf, edd, blend, bkind,
               bfm, btau0, bflo, bfhi, ch, al, hl, hr);
    fw3(m,M1_IFW_AL,k,j,i) = al;
    fw3(m,M1_IFW_HCL,k,j,i) = hl;
    fw3(m,M1_IFW_HCR,k,j,i) = hr;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransTheta
//! \brief the TRANSVERSE realizability limiter: fill thx2/thx3, the ONE face theta every
//! use of the transverse face elimination reads (the face-flux kernels, the cell-terms
//! kernel, ImplicitOffDiagOp and the post-solve reconstruction).
//!
//! The face-eliminated transverse flux
//!   F_f' = theta_f [F_f^n - c^2 dt G_f - c dt v_f g0_f],  G_f = gr_f + off_f,
//!   theta_f = 1/(1 + c dt kt_f)
//! has NO free-streaming bound.  Its steady state is F_f = -c G_f/kt_f, the unlimited
//! diffusive flux, which at c dt/dx ~ 7e3 is super-luminal wherever rho kappa is tiny:
//! in the optically thin top of the 2-D He slab |F_2|/(c E) saturates at the post-solve
//! clip of 1 and E goes horizontally unphysical within 2 s.
//!
//! Under implicit_trans_limit = lp each face gets a LAGGED limiter opacity
//!   klim_f = |G_f| / (phi_f E_f),   E_f = (E_L + E_R)/2 of the lagged iterate,
//!   phi_f  = fmax sqrt(max(1 - f1_f^2, 0.01)),  f1_f the face mean of the lagged cell
//!            x1 reduced flux F1/(c E), clipped to [-1,1],
//! and theta_f = 1/(1 + c dt (kt_f + klim_f)).  Then in steady state
//!   |F_f| = c |G_f| / (kt_f + |G_f|/(phi_f E_f)) <= c phi_f E_f,
//! i.e. the transverse flux can never exceed the room the x1 flux leaves in the
//! realizability cone, while where R = |G_f|/(kt_f E_f) << 1 -- every optically thick
//! face -- the change is O(R) and the diffusion limit is untouched.  This is a flux
//! limiter of exactly the Levermore-Pomraning form, written as an opacity so that it
//! enters the ALREADY linear face elimination.
//!
//! klim only ever makes theta_f SMALLER, and theta_f multiplies the transverse
//! off-diagonals of the 7-point row (and its own diagonal contribution by the same
//! factor), so the row stays diagonally dominant with the same signs: the M-matrix
//! property that guarantees E' > 0 at any dt is preserved.
//!
//! klim follows implicit_closure_lag: with `step` it is computed ONCE per step, from the
//! entry state (newk is then true only on the first Picard pass), and frozen for every
//! later pass and every Krylov application of that step; with `pass` it is recomputed at
//! each pass, from that pass' lagged iterate.  Either way it is a lagged coefficient and
//! the linear system stays linear.

void RadiationM1::ImplicitTransTheta(bool newk) {
  if (impl_tlim == M1_TLIM_NONE && !blat_on) return;
  if (impl_tlim == M1_TLIM_NONE) {
    // blendall-1009 without the limiter: theta_f = (1 - AL_f)/(1 + c dt kt_f), the plain
    // face theta times the central weight of the blend; read wherever lm is true
    auto &idc = pmy_pack->pmesh->mb_indcs;
    const int is0 = idc.is, ie0 = idc.ie, js0 = idc.js, je0 = idc.je;
    const int ks0 = idc.ks, ke0 = idc.ke;
    const int nmb0 = pmy_pack->nmb_thispack - 1;
    auto iw0 = iw;
    auto t2 = thx2;
    auto t3 = thx3;
    auto w2 = ifw2;
    auto w3 = ifw3;
    const Real ch0 = chat, dt0 = dt_sub;
    par_for("m1_impl_tbl2", DevExeSpace(), 0, nmb0, ks0, ke0, js0, je0+1, is0, ie0,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real ktf = 0.5*(iw0(m,M1_IW_KT,k,j-1,i) + iw0(m,M1_IW_KT,k,j,i));
      Real th = 1.0/(1.0 + ch0*dt0*ktf);
      t2(m,k,j,i) = (1.0 - w2(m,M1_IFW_AL,k,j,i))*th;
    });
    if (trans_x3) {
      par_for("m1_impl_tbl3", DevExeSpace(), 0, nmb0, ks0, ke0+1, js0, je0, is0, ie0,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real ktf = 0.5*(iw0(m,M1_IW_KT,k-1,j,i) + iw0(m,M1_IW_KT,k,j,i));
        Real th = 1.0/(1.0 + ch0*dt0*ktf);
        t3(m,k,j,i) = (1.0 - w3(m,M1_IFW_AL,k,j,i))*th;
      });
    }
    return;
  }
  const bool blt = blat_on;
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto th2_ = thx2;
  auto kl2_ = klx2;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  auto th3_ = thx3;
  auto kl3_ = klx3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs;   // CS1: seams open
  Real cl = c_light, ch = chat, dt = dt_sub;
  Real efl = e_floor;
  Real fmx = impl_tfmax;
  const bool thrd = trans_x3;
  const bool nk = newk;
  const int odm = od_now;

  par_for("m1_impl_tlim2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      // a physical x2 face is reflecting: its flux is zero and its theta is never read
      kl2_(m,k,j,i) = 0.0;
      th2_(m,k,j,i) = 1.0;
      return;
    }
    int jm = j - 1;
    if (nk) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag a6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int kl = (!thrd || ((a5 != BoundaryFlag::block) &&
                          (a5 != BoundaryFlag::periodic))) ? ks : ks-1;
      int ku = (!thrd || ((a6 != BoundaryFlag::block) &&
                          (a6 != BoundaryFlag::periodic))) ? ke : ke+1;
      int jl = plo ? js : js-1;
      int ju = phi ? je : je+1;
      Real el = fmax(iw_(m,M1_IW_EP,k,jm,i), efl);
      Real er = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real dl = M1DDiag(iw_,vd_,dfull,m,1,k,jm,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
      Real gf = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/dx2;
      if (odm != M1_OD_NONE) {
        gf += 0.5*(M1OffDiv(iw_,m,1,k,jm,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      Real rl = iw_(m,M1_IW_F1,k,jm,i)/(cl*el);
      Real rr = iw_(m,M1_IW_F1,k,j,i)/(cl*er);
      rl = fmin(fmax(rl, -1.0), 1.0);
      rr = fmin(fmax(rr, -1.0), 1.0);
      Real f1f = fmin(fmax(0.5*(rl + rr), -1.0), 1.0);
      Real phif = fmx*sqrt(fmax(1.0 - f1f*f1f, 0.01));
      Real ef = fmax(0.5*(el + er), 1.0e-300);
      kl2_(m,k,j,i) = fabs(gf)/(phif*ef);
    }
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,jm,i) + iw_(m,M1_IW_KT,k,j,i));
    th2_(m,k,j,i) = 1.0/(1.0 + ch*dt*(ktf + kl2_(m,k,j,i)));
    if (blt) {th2_(m,k,j,i) *= 1.0 - bw2_(m,M1_IFW_AL,k,j,i);}
  });

  if (!thrd) return;
  par_for("m1_impl_tlim3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) {
      kl3_(m,k,j,i) = 0.0;
      th3_(m,k,j,i) = 1.0;
      return;
    }
    int km = k - 1;
    if (nk) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
      BoundaryFlag a4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      bool pj1 = (a3 != BoundaryFlag::block) && (a3 != BoundaryFlag::periodic);
      bool pj2 = (a4 != BoundaryFlag::block) && (a4 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int jl = pj1 ? js : js-1;
      int ju = pj2 ? je : je+1;
      int kl = plo ? ks : ks-1;
      int ku = phi ? ke : ke+1;
      Real el = fmax(iw_(m,M1_IW_EP,km,j,i), efl);
      Real er = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real dl = M1DDiag(iw_,vd_,dfull,m,2,km,j,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      Real gf = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/dx3;
      if (odm != M1_OD_NONE) {
        gf += 0.5*(M1OffDiv(iw_,m,2,km,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      Real rl = iw_(m,M1_IW_F1,km,j,i)/(cl*el);
      Real rr = iw_(m,M1_IW_F1,k,j,i)/(cl*er);
      rl = fmin(fmax(rl, -1.0), 1.0);
      rr = fmin(fmax(rr, -1.0), 1.0);
      Real f1f = fmin(fmax(0.5*(rl + rr), -1.0), 1.0);
      Real phif = fmx*sqrt(fmax(1.0 - f1f*f1f, 0.01));
      Real ef = fmax(0.5*(el + er), 1.0e-300);
      kl3_(m,k,j,i) = fabs(gf)/(phif*ef);
    }
    Real ktf = 0.5*(iw_(m,M1_IW_KT,km,j,i) + iw_(m,M1_IW_KT,k,j,i));
    th3_(m,k,j,i) = 1.0/(1.0 + ch*dt*(ktf + kl3_(m,k,j,i)));
    if (blt) {th3_(m,k,j,i) *= 1.0 - bw3_(m,M1_IFW_AL,k,j,i);}
  });
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1SphTransRow
//! \brief STAGE S1 (spherical-polar wedge): the transverse part of the E row, i.e. the
//! diagonal `dia`, the lagged right-hand side `tt` and (bicgstab) the neighbour
//! coefficients CJM/CJP/CKM/CKP, with dt A_f/V_i per face and the gradient over the
//! centre-to-centre arc length.  Same terms as the Cartesian m1_impl_tcell; the face
//! fluxes fp/fm (x2) and gp/gm (x3) come from it unchanged.

template <typename V, typename VD, typename VV, typename VF, typename VT>
KOKKOS_INLINE_FUNCTION
void M1SphTransRow(const V &iw_, const VD &vd_, const bool dfull, const VV &cvol,
                   const VF &carea, const VF &cdxf, const int m, const int k,
                   const int j, const int i, const int js, const int je, const int ks,
                   const int ke, const bool p2lo, const bool p2hi, const bool thrd,
                   const BoundaryFlag b5, const BoundaryFlag b6, const bool bcg,
                   const Real ch, const Real cl, const Real dt, const Real fp,
                   const Real fm, const Real gp, const Real gm, Real &dia, Real &tt,
                   const bool blt, const VT &th2, const VT &th3, const V &bw2,
                   const V &bw3) {
  const Real cr = ch/cl;
  const Real iv = dt/cvol(m,k,j,i);
  dia = 0.0;
  Real cjp = 0.0, cjm = 0.0, ckp = 0.0, ckm = 0.0;
  const Real d2c = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
  const Real a2c = iw_(m,M1_IW_A2,k,j,i);
  const Real n2p = carea.x2f(m,k,j+1,i)*iv, n2m = carea.x2f(m,k,j,i)*iv;
  if (!(j == je && p2hi)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    if (blt) {   // blendall-1009: (1 - AL) th, and the berthon part below
      th = th2(m,k,j+1,i);
      dia += n2p*bw2(m,M1_IFW_HCL,k,j+1,i);
      cjp += n2p*bw2(m,M1_IFW_HCR,k,j+1,i);
    }
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,j,i) + iw_(m,M1_IW_V2,k,j+1,i));
    if (vf > 0.0) {
      dia += n2p*cr*a2c;
    } else {
      cjp += n2p*cr*iw_(m,M1_IW_A2,k,j+1,i);
    }
    Real g = th*ch*ch*dt/cdxf.x2f(m,k,j+1,i);
    dia += n2p*g*d2c;
    cjp -= n2p*g*M1DDiag(iw_,vd_,dfull,m,1,k,j+1,i);
  }
  if (!(j == js && p2lo)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    if (blt) {
      th = th2(m,k,j,i);
      dia -= n2m*bw2(m,M1_IFW_HCR,k,j,i);
      cjm -= n2m*bw2(m,M1_IFW_HCL,k,j,i);
    }
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,j-1,i) + iw_(m,M1_IW_V2,k,j,i));
    if (vf > 0.0) {
      cjm -= n2m*cr*iw_(m,M1_IW_A2,k,j-1,i);
    } else {
      dia -= n2m*cr*a2c;
    }
    Real g = th*ch*ch*dt/cdxf.x2f(m,k,j,i);
    dia += n2m*g*d2c;
    cjm -= n2m*g*M1DDiag(iw_,vd_,dfull,m,1,k,j-1,i);
  }
  tt = cr*(n2p*fp - n2m*fm);
  if (thrd) {
    bool p3lo = (b5 != BoundaryFlag::block) && (b5 != BoundaryFlag::periodic);
    bool p3hi = (b6 != BoundaryFlag::block) && (b6 != BoundaryFlag::periodic);
    const Real d3c = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
    const Real a3c = iw_(m,M1_IW_A3,k,j,i);
    const Real n3p = carea.x3f(m,k+1,j,i)*iv, n3m = carea.x3f(m,k,j,i)*iv;
    if (!(k == ke && p3hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (blt) {
        th = th3(m,k+1,j,i);
        dia += n3p*bw3(m,M1_IFW_HCL,k+1,j,i);
        ckp += n3p*bw3(m,M1_IFW_HCR,k+1,j,i);
      }
      Real vf = 0.5*(iw_(m,M1_IW_V3,k,j,i) + iw_(m,M1_IW_V3,k+1,j,i));
      if (vf > 0.0) {
        dia += n3p*cr*a3c;
      } else {
        ckp += n3p*cr*iw_(m,M1_IW_A3,k+1,j,i);
      }
      Real g = th*ch*ch*dt/cdxf.x3f(m,k+1,j,i);
      dia += n3p*g*d3c;
      ckp -= n3p*g*M1DDiag(iw_,vd_,dfull,m,2,k+1,j,i);
    }
    if (!(k == ks && p3lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (blt) {
        th = th3(m,k,j,i);
        dia -= n3m*bw3(m,M1_IFW_HCR,k,j,i);
        ckm -= n3m*bw3(m,M1_IFW_HCL,k,j,i);
      }
      Real vf = 0.5*(iw_(m,M1_IW_V3,k-1,j,i) + iw_(m,M1_IW_V3,k,j,i));
      if (vf > 0.0) {
        ckm -= n3m*cr*iw_(m,M1_IW_A3,k-1,j,i);
      } else {
        dia -= n3m*cr*a3c;
      }
      Real g = th*ch*ch*dt/cdxf.x3f(m,k,j,i);
      dia += n3m*g*d3c;
      ckm -= n3m*g*M1DDiag(iw_,vd_,dfull,m,2,k-1,j,i);
    }
    tt += cr*(n3p*gp - n3m*gm);
  }
  if (bcg) {
    iw_(m,M1_IW_CJM,k,j,i) = cjm;
    iw_(m,M1_IW_CJP,k,j,i) = cjp;
    iw_(m,M1_IW_CKM,k,j,i) = ckm;
    iw_(m,M1_IW_CKP,k,j,i) = ckp;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1CsTransRow
//! \brief STAGE CS1 (cubed sphere): M1SphTransRow with the skewed-grid two-point
//! coefficients.  An interior face takes the centre arc times sin of the face angle
//! (dP/dn = (dP/dl_a - cos dP/dl_b)/sin, the cross part lagged in the face flux); a
//! panel-seam face takes the canonical seam area 0.5 (r_r^2 - r_l^2) dth and the mirror
//! pair distance r angm (M1CsSeamGeom), the same numbers on both panels.

template <typename V, typename VD, typename VV, typename VF, typename VS, typename VG,
          typename VB, typename VX>
KOKKOS_INLINE_FUNCTION
void M1CsTransRow(const V &iw_, const VD &vd_, const bool dfull, const VV &cvol,
                  const VF &carea, const VF &cdxf, const VS &sn2, const VS &sn3,
                  const VG &g2, const VG &g3, const VB &cseam, const VX &x1v,
                  const VX &x1f, const int m, const int k, const int j, const int i,
                  const int js, const int je, const int ks, const int ke,
                  const bool p2lo, const bool p2hi, const bool thrd,
                  const bool p3lo, const bool p3hi, const bool bcg,
                  const Real ch, const Real cl, const Real dt, const Real fp,
                  const Real fm, const Real gp, const Real gm, Real &dia, Real &tt) {
  const Real cr = ch/cl;
  const Real iv = dt/cvol(m,k,j,i);
  const Real dsh = 0.5*(x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i));
  const bool s2lo = (cseam(m,0) != 0), s2hi = (cseam(m,1) != 0);
  const bool s3lo = (cseam(m,2) != 0), s3hi = (cseam(m,3) != 0);
  dia = 0.0;
  Real cjp = 0.0, cjm = 0.0, ckp = 0.0, ckm = 0.0;
  const Real d2c = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
  const Real a2c = iw_(m,M1_IW_A2,k,j,i);
  // the x2 faces j+1 (p) and j (m): area and two-point distance
  const bool sp2 = (j == je) && s2hi, sm2 = (j == js) && s2lo;
  const Real a2p = sp2 ? dsh*g2(m,0,k,1) : carea.x2f(m,k,j+1,i);
  const Real a2m = sm2 ? dsh*g2(m,0,k,0) : carea.x2f(m,k,j,i);
  const Real l2p = sp2 ? x1v(m,i)*g2(m,1,k,1) : cdxf.x2f(m,k,j+1,i)*sn2(m,k,j+1);
  const Real l2m = sm2 ? x1v(m,i)*g2(m,1,k,0) : cdxf.x2f(m,k,j,i)*sn2(m,k,j);
  const Real n2p = a2p*iv, n2m = a2m*iv;
  if (!(j == je && p2hi)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,j,i) + iw_(m,M1_IW_V2,k,j+1,i));
    if (vf > 0.0) {
      dia += n2p*cr*a2c;
    } else {
      cjp += n2p*cr*iw_(m,M1_IW_A2,k,j+1,i);
    }
    Real g = th*ch*ch*dt/l2p;
    dia += n2p*g*d2c;
    cjp -= n2p*g*M1DDiag(iw_,vd_,dfull,m,1,k,j+1,i);
  }
  if (!(j == js && p2lo)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,j-1,i) + iw_(m,M1_IW_V2,k,j,i));
    if (vf > 0.0) {
      cjm -= n2m*cr*iw_(m,M1_IW_A2,k,j-1,i);
    } else {
      dia -= n2m*cr*a2c;
    }
    Real g = th*ch*ch*dt/l2m;
    dia += n2m*g*d2c;
    cjm -= n2m*g*M1DDiag(iw_,vd_,dfull,m,1,k,j-1,i);
  }
  tt = cr*(n2p*fp - n2m*fm);
  if (thrd) {
    const Real d3c = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
    const Real a3c = iw_(m,M1_IW_A3,k,j,i);
    const bool sp3 = (k == ke) && s3hi, sm3 = (k == ks) && s3lo;
    const Real a3p = sp3 ? dsh*g3(m,0,j,1) : carea.x3f(m,k+1,j,i);
    const Real a3m = sm3 ? dsh*g3(m,0,j,0) : carea.x3f(m,k,j,i);
    const Real l3p = sp3 ? x1v(m,i)*g3(m,1,j,1) : cdxf.x3f(m,k+1,j,i)*sn3(m,k+1,j);
    const Real l3m = sm3 ? x1v(m,i)*g3(m,1,j,0) : cdxf.x3f(m,k,j,i)*sn3(m,k,j);
    const Real n3p = a3p*iv, n3m = a3m*iv;
    if (!(k == ke && p3hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V3,k,j,i) + iw_(m,M1_IW_V3,k+1,j,i));
      if (vf > 0.0) {
        dia += n3p*cr*a3c;
      } else {
        ckp += n3p*cr*iw_(m,M1_IW_A3,k+1,j,i);
      }
      Real g = th*ch*ch*dt/l3p;
      dia += n3p*g*d3c;
      ckp -= n3p*g*M1DDiag(iw_,vd_,dfull,m,2,k+1,j,i);
    }
    if (!(k == ks && p3lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V3,k-1,j,i) + iw_(m,M1_IW_V3,k,j,i));
      if (vf > 0.0) {
        ckm -= n3m*cr*iw_(m,M1_IW_A3,k-1,j,i);
      } else {
        dia -= n3m*cr*a3c;
      }
      Real g = th*ch*ch*dt/l3m;
      dia += n3m*g*d3c;
      ckm -= n3m*g*M1DDiag(iw_,vd_,dfull,m,2,k-1,j,i);
    }
    tt += cr*(n3p*gp - n3m*gm);
  }
  if (bcg) {
    iw_(m,M1_IW_CJM,k,j,i) = cjm;
    iw_(m,M1_IW_CJP,k,j,i) = cjp;
    iw_(m,M1_IW_CKM,k,j,i) = ckm;
    iw_(m,M1_IW_CKP,k,j,i) = ckp;
  }
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransverseTerms
//! \brief milestone 3b phase B: one line-Jacobi evaluation of the TRANSVERSE (x2, x3)
//! part of the operator at the current Picard iterate.
//!
//! Per direction d and face f, exactly the face-eliminated form the x1 faces use:
//!
//!   F0_f' = theta_f [ F0_f^n - c^2 dt (P_dd,R - P_dd,L)/dx_d - c dt v_f g0_f
//!                     - c^2 dt (sum_{e != d} d_e P_de)_f ]
//!   theta_f = 1/(1 + c dt (rho kappa_t)_f),  (rho kappa_t)_f the arithmetic face mean
//!
//! with P_dd = D_dd E of the LAGGED closure, D_dd = (1-chi)/2 + (3 chi - 1)/2 n_d^2, and
//! the off-diagonal divergence fully lagged (centred differences of the previous pass' E
//! and closure), so it is a pure right-hand-side term.  The advective enthalpy flux
//! A_d = v_d E + (v.P)_d is upwinded with the face velocity, exactly as in x1.
//!
//! The result is split into
//!   M1_IW_TDIA = dT/dE_c >= 0, which goes on the matrix DIAGONAL, and
//!   M1_IW_TRHS = -(T(E^k) - TDIA E^k_c), the neighbours' lagged contribution.
//! Keeping the diagonal part on the matrix is what preserves the full 7-point M-matrix:
//! the column sums stay >= 1 and E' > 0 at any dt.
//!
//! It also measures the TRUE residual of the full 7-point linear system.  After the line
//! solve E^{k+1} satisfies  D1(E^{k+1}) + TDIA E^{k+1} + U(E^k) = b  with
//! U(E) = T(E) - TDIA E_c, so the residual of the FULL system at E^{k+1} is exactly
//! U(E^k) - U(E^{k+1}): the change of the lagged off-diagonal term between two passes.
//! That is what M1_IW_LRES holds, and it is tested separately from the Picard residual --
//! a lagged coupling can stall and look converged on |dE|/E alone.

void RadiationM1::ImplicitTransverseTerms(bool first) {
  if (!trans_on) return;
  // blendall-1009: the x2/x3 face coefficients, frozen for the step under
  // implicit_recon_lag = step (the default with implicit_flux_faces = all)
  if (blat_on && (first || !impl_recon_freeze)) {ImplicitLatFaceCoef();}
  // the transverse realizability limiter: one evaluation of theta for the whole pass,
  // which everything below and the Krylov operator then READ.  klim itself follows the
  // closure lag (frozen for the step under implicit_closure_lag = step).
  ImplicitTransTheta(!impl_clag_step || first);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto f2_ = f0x2;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  auto f2n_ = f0x2n;
  auto f3_ = f0x3;
  auto f3n_ = f0x3n;
  // the transverse realizability limiter.  Under `none` the OLD expression for theta is
  // taken, term for term, so the arithmetic of every earlier configuration is bitwise
  // unchanged; under `lp` every use below reads the SAME thx2/thx3 the Krylov operator
  // reads.
  const bool lm = (impl_tlim != M1_TLIM_NONE) || blat_on;
  auto th2_ = thx2;
  auto th3_ = thx3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs;   // CS1: seams open
  Real cl = c_light, ch = chat, dt = dt_sub;
  // blendall-1009: the berthon part HCL E_L + HCR E_R of the x2/x3 faces (ifw2/ifw3)
  const bool blt = blat_on;
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  const bool thrd = trans_x3;
  const bool fst = first;
  const Real wmem = dbg_trans_memory;
  // MILESTONE 3b phase C: also store the transverse OFF-DIAGONAL coefficients of the
  // frozen row, which is what the BiCGStab operator applies and what turns TRHS back into
  // the right-hand side of the full system.
  const bool bcg = bicg_on;
  // MILESTONE 3b phase D: implicit_offdiag.  The face flux stored here is the PHYSICAL
  // one and keeps the off-diagonal term at the current iterate under `lagged` and
  // `operator` alike; what changes between the two is the row the solve is given (see
  // ImplicitSolve step (e), where the term is subtracted from the right-hand side and
  // handed to the operator instead).  Under `none` it is dropped everywhere.
  const int odm = od_now;
  // STAGE S1 (rad_m1_sph.cpp): on the spherical-polar wedge the face-flux equations take
  // the centre-to-centre arc length dxface (r dtheta, r sin(theta) dphi) and the E row
  // takes dt A_f/V_i per face.  Every Cartesian expression below is left textually
  // untouched; the sp forms are separate `if (sph)` blocks that OVERWRITE the results
  // (hipcc contracts FMAs differently once a factor becomes a run-time select, which
  // broke the GPU Cartesian bitwise gate of the first version).
  const bool sph = sph_geom;
  auto cvol = pmy_pack->pcoord->volume;
  auto carea = pmy_pack->pcoord->area;
  auto cdxf = pmy_pack->pcoord->dxface;
  // STAGE S2 (sph_q: a chi(f) closure on the wedge): the lagged part of the face
  // equation is M1SphCurv (curvature always, off-diagonal terms under `lagged`), which
  // OVERWRITES the Cartesian `off` below
  const bool sphq = sph_q;
  const bool odl = (odm != M1_OD_NONE);
  // vet_col_lat (m1-vetcol-lat): the lagged lateral off-diagonal term M1SphLat
  const bool vlat = vlat_now && vlat_ready && sph;
  auto vlt_ = tau_ten;
  const int c0l = M1_TT_LAT0;
  // vet_gd (m1-vet-gd): the lagged tangential anisotropy M1SphTan
  const bool vtan = vgd_on && vgd_tan && vlat_ready && sph;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx2v = pmy_pack->pcoord->x2v;
  auto cx3v = pmy_pack->pcoord->x3v;
  // STAGE CS1 (cubed sphere): the face-normal gradient on the skewed panel grid, as
  // `if (csg)` overwrites of the sp result.  An interior face takes
  //   dP/dn = (dP/dl_a - cos dP/dl_b)/sin,
  // the two-point part over the centre arc dxface (implicit: M1CsTransRow) and the cross
  // term from M1CsTanDer of the iterate (lagged into TRHS through the face flux, the
  // Picard contraction <= |cos|/... <= 1/2).  A PANEL-SEAM face pairs the cell with its
  // mirror image (the no-resample ghost): the pair chord is normal to the seam, so the
  // two-point difference over the canonical pair distance r angm IS the normal gradient,
  // at the foot of the chord; the lagged correction del * d/ds moves it along the seam to
  // the face midpoint (del = smid - sfoot, d/ds over the neighbouring seam faces, one-
  // sided at a cube vertex).  Both panels evaluate the same canonical numbers.
  const bool csg = cs_geom;
  auto cseam = cs_seam.d_view;
  auto csg2_ = csg2;
  auto csg3_ = csg3;
  auto csn2 = pmy_pack->pcoord->sin_face_xi;
  auto ccs2 = pmy_pack->pcoord->cos_face_xi;
  auto csn3 = pmy_pack->pcoord->sin_face_eta;
  auto ccs3 = pmy_pack->pcoord->cos_face_eta;
  auto cx1f = pmy_pack->pcoord->xx1f;

  // (1) the x2 face fluxes
  const bool must = muscl_now;   // xthinfix-1009 Fix B (sig of the latest build)
  const int imbt = iw_muscl;
  par_for_lb("m1_impl_f2face", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      f2_(m,k,j,i) = 0.0;
      return;
    }
    Real dx1 = mbsize.d_view(m).dx1;
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
    BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
    BoundaryFlag a5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag a6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
    bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
    int il = pi1 ? is : is-1;
    int iu = pi2 ? ie : ie+1;
    int kl = (!thrd || ((a5 != BoundaryFlag::block) &&
                        (a5 != BoundaryFlag::periodic))) ? ks : ks-1;
    int ku = (!thrd || ((a6 != BoundaryFlag::block) &&
                        (a6 != BoundaryFlag::periodic))) ? ke : ke+1;
    int jl = plo ? js : js-1;
    int ju = phi ? je : je+1;
    int jm = j - 1;
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,jm,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
    Real dl = M1DDiag(iw_,vd_,dfull,m,1,k,jm,i);
    Real dr = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
    if (vtan) {   // vet_gd: D_tt = (1 - D_rr)/2 + a in the compact face gradient
      dl += vlt_(m,c0l+3,k,jm,i);
      dr += vlt_(m,c0l+3,k,j,i);
    }
    Real gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/dx2;
    if (sph) {
      gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/cdxf.x2f(m,k,j,i);
    }
    if (csg) {
      const bool s2lo = (cseam(m,0) != 0), s2hi = (cseam(m,1) != 0);
      const bool s3lo = (cseam(m,2) != 0), s3hi = (cseam(m,3) != 0);
      if ((j == js && s2lo) || (j == je+1 && s2hi)) {
        const int sd = (j == js) ? 0 : 1;
        auto dpr = [&](const int kq) {
          return (M1DDiag(iw_,vd_,dfull,m,1,kq,j,i)*iw_(m,M1_IW_EP,kq,j,i)
                  - M1DDiag(iw_,vd_,dfull,m,1,kq,jm,i)*iw_(m,M1_IW_EP,kq,jm,i))
                 /(cx1v(m,i)*csg2_(m,1,kq,sd));
        };
        const Real d0 = dpr(k);
        const bool up = !(k == ke && s3hi), dn = !(k == ks && s3lo);
        Real dd = 0.0;
        if (up && dn) {
          dd = (dpr(k+1) - dpr(k-1))/(csg2_(m,2,k+1,sd) - csg2_(m,2,k-1,sd));
        } else if (up) {
          dd = (dpr(k+1) - d0)/(csg2_(m,2,k+1,sd) - csg2_(m,2,k,sd));
        } else if (dn) {
          dd = (d0 - dpr(k-1))/(csg2_(m,2,k,sd) - csg2_(m,2,k-1,sd));
        }
        gr = d0 + (csg2_(m,3,k,sd) - csg2_(m,2,k,sd))*dd;
      } else {
        const Real ge = 0.5*(M1CsTanDer(iw_,vd_,dfull,cdxf,m,2,k,jm,i,ks,ke,s3lo,s3hi)
                             + M1CsTanDer(iw_,vd_,dfull,cdxf,m,2,k,j,i,ks,ke,s3lo,s3hi));
        gr = ((dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/cdxf.x2f(m,k,j,i)
              - ccs2(m,k,j)*ge)/csn2(m,k,j);
      }
    }
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,jm,i) + iw_(m,M1_IW_V2,k,j,i));
    Real g0f = 0.5*(iw_(m,M1_IW_G0,k,jm,i) + iw_(m,M1_IW_G0,k,j,i));
    Real off = 0.0;
    if (odm != M1_OD_NONE) {
      off = 0.5*(M1OffDiv(iw_,m,1,k,jm,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                          dfull)
                 + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                            dfull));
    }
    if (sphq) {
      off = 0.5*(M1SphCurv(iw_,cx1v,cx2v,cx3v,m,1,k,jm,i,odl,thrd,il,iu,jl,ju,kl,ku,
                           M1_IW_EP)
                 + M1SphCurv(iw_,cx1v,cx2v,cx3v,m,1,k,j,i,odl,thrd,il,iu,jl,ju,kl,ku,
                             M1_IW_EP));
    }
    if (vlat) {
      off += 0.5*(M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,1,k,jm,i,thrd,il,iu,jl,
                             ju,kl,ku,M1_IW_EP)
                 + M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,1,k,j,i,thrd,il,iu,jl,
                            ju,kl,ku,M1_IW_EP));
    }
    if (vtan) {
      off += 0.5*(M1SphTan(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,1,k,jm,i,thrd,jl,ju,kl,ku,
                           M1_IW_EP)
                 + M1SphTan(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,1,k,j,i,thrd,jl,ju,kl,ku,
                            M1_IW_EP));
    }
    if (csg) {
      // STAGE CS2: THE LATERAL-CLOSURE HOOK of the cs transverse face (CS2_RESULTS.md,
      // "interface").  With the closures cs runs (eddington, vet_col about r_hat) the
      // tensor is uniaxial about r_hat, P = p_t (I - r r) + p_r r r, and its tangential
      // divergence is grad_t p_t: the two-point + cross-term gradient above carries all
      // of it, the sp curvature M1SphCurv (theta-based) does not apply, and the lagged
      // remainder is zero.  A lateral/off-diagonal correction of the tensor (tau_ten
      // slots M1_TT_LAT, the planned SC correction of vet_col) enters HERE as a lagged
      // term of the face equation, (d_r P_r,lat + curvature + d_t of the non-isotropic
      // tangential part), exactly as `off` does on Cartesian and sp meshes.
      off = 0.0;
    }
    f2_(m,k,j,i) = th*(wmem*f2n_(m,k,j,i) - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*off);
    if (blt) {
      f2_(m,k,j,i) += (cl/ch)*(bw2_(m,M1_IFW_HCL,k,j,i)*iw_(m,M1_IW_EP,k,jm,i)
                               + bw2_(m,M1_IFW_HCR,k,j,i)*iw_(m,M1_IW_EP,k,j,i));
      if (must) {
        const Real hcl = bw2_(m,M1_IFW_HCL,k,j,i), hcr = bw2_(m,M1_IFW_HCR,k,j,i);
        Real gm = 0.0;
        if (hcl != 0.0) {
          gm += 0.25*hcl*iw_(m,imbt+M1_IM_SIG+1,k,jm,i)
                *(iw_(m,M1_IW_EP,k,j,i) - iw_(m,M1_IW_EP,k,jm-1,i));
        }
        if (hcr != 0.0) {
          gm -= 0.25*hcr*iw_(m,imbt+M1_IM_SIG+1,k,j,i)
                *(iw_(m,M1_IW_EP,k,j+1,i) - iw_(m,M1_IW_EP,k,jm,i));
        }
        f2_(m,k,j,i) += (cl/ch)*gm;
      }
    }
  });

  // (2) the x3 face fluxes
  if (thrd) {
    par_for_lb("m1_impl_f3face", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
      bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
      if ((k == ks && plo) || (k == ke+1 && phi)) {
        f3_(m,k,j,i) = 0.0;
        return;
      }
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
      BoundaryFlag a4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      bool pj1 = (a3 != BoundaryFlag::block) && (a3 != BoundaryFlag::periodic);
      bool pj2 = (a4 != BoundaryFlag::block) && (a4 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int jl = pj1 ? js : js-1;
      int ju = pj2 ? je : je+1;
      int kl = plo ? ks : ks-1;
      int ku = phi ? ke : ke+1;
      int km = k - 1;
      Real ktf = 0.5*(iw_(m,M1_IW_KT,km,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real dl = M1DDiag(iw_,vd_,dfull,m,2,km,j,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      if (vtan) {   // vet_gd: D_pp = (1 - D_rr)/2 - a in the compact face gradient
        dl += vlt_(m,c0l+5,km,j,i);
        dr += vlt_(m,c0l+5,k,j,i);
      }
      Real gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/dx3;
      if (sph) {
        gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/cdxf.x3f(m,k,j,i);
      }
      if (csg) {
        const bool s2lo = (cseam(m,0) != 0), s2hi = (cseam(m,1) != 0);
        const bool s3lo = (cseam(m,2) != 0), s3hi = (cseam(m,3) != 0);
        if ((k == ks && s3lo) || (k == ke+1 && s3hi)) {
          const int sd = (k == ks) ? 0 : 1;
          auto dpr = [&](const int jq) {
            return (M1DDiag(iw_,vd_,dfull,m,2,k,jq,i)*iw_(m,M1_IW_EP,k,jq,i)
                    - M1DDiag(iw_,vd_,dfull,m,2,km,jq,i)*iw_(m,M1_IW_EP,km,jq,i))
                   /(cx1v(m,i)*csg3_(m,1,jq,sd));
          };
          const Real d0 = dpr(j);
          const bool up = !(j == je && s2hi), dn = !(j == js && s2lo);
          Real dd = 0.0;
          if (up && dn) {
            dd = (dpr(j+1) - dpr(j-1))/(csg3_(m,2,j+1,sd) - csg3_(m,2,j-1,sd));
          } else if (up) {
            dd = (dpr(j+1) - d0)/(csg3_(m,2,j+1,sd) - csg3_(m,2,j,sd));
          } else if (dn) {
            dd = (d0 - dpr(j-1))/(csg3_(m,2,j,sd) - csg3_(m,2,j-1,sd));
          }
          gr = d0 + (csg3_(m,3,j,sd) - csg3_(m,2,j,sd))*dd;
        } else {
          const Real ge = 0.5*(M1CsTanDer(iw_,vd_,dfull,cdxf,m,1,km,j,i,js,je,s2lo,s2hi)
                               + M1CsTanDer(iw_,vd_,dfull,cdxf,m,1,k,j,i,js,je,s2lo,
                                            s2hi));
          gr = ((dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/cdxf.x3f(m,k,j,i)
                - ccs3(m,k,j)*ge)/csn3(m,k,j);
        }
      }
      Real vf = 0.5*(iw_(m,M1_IW_V3,km,j,i) + iw_(m,M1_IW_V3,k,j,i));
      Real g0f = 0.5*(iw_(m,M1_IW_G0,km,j,i) + iw_(m,M1_IW_G0,k,j,i));
      Real off = 0.0;
      if (odm != M1_OD_NONE) {
        off = 0.5*(M1OffDiv(iw_,m,2,km,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      if (sphq) {
        off = 0.5*(M1SphCurv(iw_,cx1v,cx2v,cx3v,m,2,km,j,i,odl,thrd,il,iu,jl,ju,kl,ku,
                             M1_IW_EP)
                   + M1SphCurv(iw_,cx1v,cx2v,cx3v,m,2,k,j,i,odl,thrd,il,iu,jl,ju,kl,ku,
                               M1_IW_EP));
      }
      if (vlat) {
        off += 0.5*(M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,2,km,j,i,thrd,il,iu,jl,
                               ju,kl,ku,M1_IW_EP)
                   + M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,2,k,j,i,thrd,il,iu,jl,
                              ju,kl,ku,M1_IW_EP));
      }
      if (vtan) {
        off += 0.5*(M1SphTan(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,2,km,j,i,thrd,jl,ju,kl,ku,
                             M1_IW_EP)
                   + M1SphTan(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,2,k,j,i,thrd,jl,ju,kl,ku,
                              M1_IW_EP));
      }
      if (csg) {off = 0.0;}   // STAGE CS2: the lateral-closure hook (see the x2 face)
      f3_(m,k,j,i) = th*(wmem*f3n_(m,k,j,i)
                         - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*off);
      if (blt) {
        f3_(m,k,j,i) += (cl/ch)*(bw3_(m,M1_IFW_HCL,k,j,i)*iw_(m,M1_IW_EP,km,j,i)
                                 + bw3_(m,M1_IFW_HCR,k,j,i)*iw_(m,M1_IW_EP,k,j,i));
        if (must) {
          const Real hcl = bw3_(m,M1_IFW_HCL,k,j,i), hcr = bw3_(m,M1_IFW_HCR,k,j,i);
          Real gm = 0.0;
          if (hcl != 0.0) {
            gm += 0.25*hcl*iw_(m,imbt+M1_IM_SIG+2,km,j,i)
                  *(iw_(m,M1_IW_EP,k,j,i) - iw_(m,M1_IW_EP,km-1,j,i));
          }
          if (hcr != 0.0) {
            gm -= 0.25*hcr*iw_(m,imbt+M1_IM_SIG+2,k,j,i)
                  *(iw_(m,M1_IW_EP,k+1,j,i) - iw_(m,M1_IW_EP,km,j,i));
          }
          f3_(m,k,j,i) += (cl/ch)*gm;
        }
      }
    });
  }

  // (3) the cell terms: the diagonal part, the lagged right-hand side and the residual
  // of the full 7-point system.  implicit_enthalpy adds its deferred correction to the
  // face fluxes fp/fm/gp/gm, which reach only the right-hand side (TRHS), never the row.
  const int enm = impl_enth;
  const bool enth2 = (enm != M1_IENTH_UPWIND);
  const int ngh = indcs.ng;
  // hesdirk2 stage solve, time2_enth_vel = start: see M1EnthCorrT
  const bool t2vs = (t2_afmode != 0) && impl_vimp &&
                    (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  const bool t2afc = (t2_afmode == 2);
  const int t2da = impl_vimp ? (iw_vimp + M1_IV_DA) : 0;
  par_for_lb("m1_impl_tcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    Real cr = ch/cl;
    Real ec = iw_(m,M1_IW_EP,k,j,i);
    Real dia = 0.0, tt = 0.0;
    BoundaryFlag b3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag b4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool p2lo = (b3 != BoundaryFlag::block) && (b3 != BoundaryFlag::periodic);
    bool p2hi = (b4 != BoundaryFlag::block) && (b4 != BoundaryFlag::periodic);
    Real nu2 = dt/dx2;
    Real d2c = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
    Real a2c = iw_(m,M1_IW_A2,k,j,i);
    Real fp = 0.0, fm = 0.0;
    Real cjp = 0.0, cjm = 0.0;
    if (!(j == je && p2hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
      Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V2,k,j,i) + iw_(m,M1_IW_V2,k,j+1,i));
      fp = f2_(m,k,j+1,i);
      if (vf > 0.0) {
        fp += a2c*ec;
        dia += nu2*cr*a2c;
      } else {
        fp += iw_(m,M1_IW_A2,k,j+1,i)*iw_(m,M1_IW_EP,k,j+1,i);
        if (bcg) {cjp += nu2*cr*iw_(m,M1_IW_A2,k,j+1,i);}
      }
      if (enth2) {
        bool o0, o3;
        int j0 = M1EnthIdx(j-1, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o0);
        int j3 = M1EnthIdx(j+2, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o3);
        if (t2vs) {
          fp += M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j0,i), ec, iw_(m,M1_IW_EP,k,j+1,i),
                            iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                            a2c, iw_(m,M1_IW_A2,k,j+1,i), iw_(m,M1_IW_A2,k,j3,i),
                            iw_(m,t2da+1,k,j0,i), iw_(m,t2da+1,k,j,i),
                            iw_(m,t2da+1,k,j+1,i), iw_(m,t2da+1,k,j3,i), vf, t2afc);
        } else {
        fp += M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j0,i), ec, iw_(m,M1_IW_EP,k,j+1,i),
                         iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                         a2c, iw_(m,M1_IW_A2,k,j+1,i), iw_(m,M1_IW_A2,k,j3,i), vf);
        }
      }
      dia += nu2*th*ch*ch*dt*d2c/dx2;
      if (bcg) {
        Real d2p = M1DDiag(iw_,vd_,dfull,m,1,k,j+1,i);
        cjp -= nu2*th*ch*ch*dt*d2p/dx2;
      }
      if (blt) {   // blendall-1009: HCL >= 0 on the diagonal, HCR <= 0 off it
        dia += nu2*bw2_(m,M1_IFW_HCL,k,j+1,i);
        if (bcg) {cjp += nu2*bw2_(m,M1_IFW_HCR,k,j+1,i);}
      }
    }
    if (!(j == js && p2lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V2,k,j-1,i) + iw_(m,M1_IW_V2,k,j,i));
      fm = f2_(m,k,j,i);
      if (vf > 0.0) {
        fm += iw_(m,M1_IW_A2,k,j-1,i)*iw_(m,M1_IW_EP,k,j-1,i);
        if (bcg) {cjm -= nu2*cr*iw_(m,M1_IW_A2,k,j-1,i);}
      } else {
        fm += a2c*ec;
        dia -= nu2*cr*a2c;
      }
      if (enth2) {
        bool o0, o3;
        int j0 = M1EnthIdx(j-2, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o0);
        int j3 = M1EnthIdx(j+1, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o3);
        if (t2vs) {
          fm += M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j0,i), iw_(m,M1_IW_EP,k,j-1,i), ec,
                            iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                            iw_(m,M1_IW_A2,k,j-1,i), a2c, iw_(m,M1_IW_A2,k,j3,i),
                            iw_(m,t2da+1,k,j0,i), iw_(m,t2da+1,k,j-1,i),
                            iw_(m,t2da+1,k,j,i), iw_(m,t2da+1,k,j3,i), vf, t2afc);
        } else {
        fm += M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j0,i), iw_(m,M1_IW_EP,k,j-1,i), ec,
                         iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                         iw_(m,M1_IW_A2,k,j-1,i), a2c, iw_(m,M1_IW_A2,k,j3,i), vf);
        }
      }
      dia += nu2*th*ch*ch*dt*d2c/dx2;
      if (bcg) {
        Real d2m = M1DDiag(iw_,vd_,dfull,m,1,k,j-1,i);
        cjm -= nu2*th*ch*ch*dt*d2m/dx2;
      }
      if (blt) {
        dia -= nu2*bw2_(m,M1_IFW_HCR,k,j,i);
        if (bcg) {cjm -= nu2*bw2_(m,M1_IFW_HCL,k,j,i);}
      }
    }
    tt += nu2*cr*(fp - fm);
    Real gps = 0.0, gms = 0.0;   // (sp) the x3 face fluxes, copied out below
    if (bcg) {
      iw_(m,M1_IW_CJM,k,j,i) = cjm;
      iw_(m,M1_IW_CJP,k,j,i) = cjp;
      iw_(m,M1_IW_CKM,k,j,i) = 0.0;
      iw_(m,M1_IW_CKP,k,j,i) = 0.0;
    }
    if (thrd) {
      BoundaryFlag b5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag b6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool p3lo = (b5 != BoundaryFlag::block) && (b5 != BoundaryFlag::periodic);
      bool p3hi = (b6 != BoundaryFlag::block) && (b6 != BoundaryFlag::periodic);
      Real nu3 = dt/dx3;
      Real d3c = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      Real a3c = iw_(m,M1_IW_A3,k,j,i);
      Real gp = 0.0, gm = 0.0;
      Real ckp = 0.0, ckm = 0.0;
      if (!(k == ke && p3hi)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
        Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V3,k,j,i) + iw_(m,M1_IW_V3,k+1,j,i));
        gp = f3_(m,k+1,j,i);
        if (vf > 0.0) {
          gp += a3c*ec;
          dia += nu3*cr*a3c;
        } else {
          gp += iw_(m,M1_IW_A3,k+1,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
          if (bcg) {ckp += nu3*cr*iw_(m,M1_IW_A3,k+1,j,i);}
        }
        if (enth2) {
          bool o0, o3;
          int k0 = M1EnthIdx(k-1, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o0);
          int k3 = M1EnthIdx(k+2, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o3);
          if (t2vs) {
            gp += M1EnthCorrT(enm, iw_(m,M1_IW_EP,k0,j,i), ec, iw_(m,M1_IW_EP,k+1,j,i),
                              iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                              iw_(m,M1_IW_A3,k0,j,i), a3c, iw_(m,M1_IW_A3,k+1,j,i),
                              iw_(m,M1_IW_A3,k3,j,i), iw_(m,t2da+2,k0,j,i),
                              iw_(m,t2da+2,k,j,i), iw_(m,t2da+2,k+1,j,i),
                              iw_(m,t2da+2,k3,j,i), vf, t2afc);
          } else {
          gp += M1EnthCorr(enm, iw_(m,M1_IW_EP,k0,j,i), ec, iw_(m,M1_IW_EP,k+1,j,i),
                           iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                           iw_(m,M1_IW_A3,k0,j,i), a3c, iw_(m,M1_IW_A3,k+1,j,i),
                           iw_(m,M1_IW_A3,k3,j,i), vf);
          }
        }
        dia += nu3*th*ch*ch*dt*d3c/dx3;
        if (bcg) {
          Real d3p = M1DDiag(iw_,vd_,dfull,m,2,k+1,j,i);
          ckp -= nu3*th*ch*ch*dt*d3p/dx3;
        }
        if (blt) {
          dia += nu3*bw3_(m,M1_IFW_HCL,k+1,j,i);
          if (bcg) {ckp += nu3*bw3_(m,M1_IFW_HCR,k+1,j,i);}
        }
      }
      if (!(k == ks && p3lo)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
        Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V3,k-1,j,i) + iw_(m,M1_IW_V3,k,j,i));
        gm = f3_(m,k,j,i);
        if (vf > 0.0) {
          gm += iw_(m,M1_IW_A3,k-1,j,i)*iw_(m,M1_IW_EP,k-1,j,i);
          if (bcg) {ckm -= nu3*cr*iw_(m,M1_IW_A3,k-1,j,i);}
        } else {
          gm += a3c*ec;
          dia -= nu3*cr*a3c;
        }
        if (enth2) {
          bool o0, o3;
          int k0 = M1EnthIdx(k-2, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o0);
          int k3 = M1EnthIdx(k+1, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o3);
          if (t2vs) {
            gm += M1EnthCorrT(enm, iw_(m,M1_IW_EP,k0,j,i), iw_(m,M1_IW_EP,k-1,j,i), ec,
                              iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                              iw_(m,M1_IW_A3,k0,j,i), iw_(m,M1_IW_A3,k-1,j,i), a3c,
                              iw_(m,M1_IW_A3,k3,j,i), iw_(m,t2da+2,k0,j,i),
                              iw_(m,t2da+2,k-1,j,i), iw_(m,t2da+2,k,j,i),
                              iw_(m,t2da+2,k3,j,i), vf, t2afc);
          } else {
          gm += M1EnthCorr(enm, iw_(m,M1_IW_EP,k0,j,i), iw_(m,M1_IW_EP,k-1,j,i), ec,
                           iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                           iw_(m,M1_IW_A3,k0,j,i), iw_(m,M1_IW_A3,k-1,j,i), a3c,
                           iw_(m,M1_IW_A3,k3,j,i), vf);
          }
        }
        dia += nu3*th*ch*ch*dt*d3c/dx3;
        if (bcg) {
          Real d3m = M1DDiag(iw_,vd_,dfull,m,2,k-1,j,i);
          ckm -= nu3*th*ch*ch*dt*d3m/dx3;
        }
        if (blt) {
          dia -= nu3*bw3_(m,M1_IFW_HCR,k,j,i);
          if (bcg) {ckm -= nu3*bw3_(m,M1_IFW_HCL,k,j,i);}
        }
      }
      tt += nu3*cr*(gp - gm);
      gps = gp;
      gms = gm;
      if (bcg) {
        iw_(m,M1_IW_CKM,k,j,i) = ckm;
        iw_(m,M1_IW_CKP,k,j,i) = ckp;
      }
    }
    if (sph) {
      // STAGE S1: the same row with dt A_f/V_i per face and the centre distances.  The
      // face fluxes fp, fm, gps, gms (face-normal F0 + enthalpy flux, corrections
      // included) carry no geometry and are reused; everything the geometry multiplies
      // is rebuilt.  (implicit_trans_limit = none on sp, so theta is the plain form.)
      M1SphTransRow(iw_, vd_, dfull, cvol, carea, cdxf, m, k, j, i, js, je, ks, ke,
                    p2lo, p2hi, thrd, mbbcs.d_view(m,BoundaryFace::inner_x3),
                    mbbcs.d_view(m,BoundaryFace::outer_x3), bcg, ch, cl, dt, fp, fm,
                    gps, gms, dia, tt, blt, th2_, th3_, bw2_, bw3_);
    }
    if (csg) {
      // STAGE CS1: the same row on the skewed panel grid and the seam pairs
      BoundaryFlag c5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag c6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool q3lo = (c5 != BoundaryFlag::block) && (c5 != BoundaryFlag::periodic);
      bool q3hi = (c6 != BoundaryFlag::block) && (c6 != BoundaryFlag::periodic);
      M1CsTransRow(iw_, vd_, dfull, cvol, carea, cdxf, csn2, csn3, csg2_, csg3_, cseam,
                   cx1v, cx1f, m, k, j, i, js, je, ks, ke, p2lo, p2hi, thrd, q3lo, q3hi,
                   bcg, ch, cl, dt, fp, fm, gps, gms, dia, tt);
    }
    Real unew = tt - dia*ec;
    Real uold = -iw_(m,M1_IW_TRHS,k,j,i);
    iw_(m,M1_IW_LRES,k,j,i) = fst ? 1.0e30 : fabs(unew - uold);
    iw_(m,M1_IW_TDIA,k,j,i) = dia;
    iw_(m,M1_IW_TRHS,k,j,i) = -unew;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGatherSolve
//! \brief gather the assembled rows (a,b,c,r) of every x1 stack onto its root block, run
//! there the IDENTICAL serial Thomas sweep the single-block solve runs, and scatter the
//! solution back into iw(M1_IW_S2).  Two messages per Picard iteration per non-root
//! block (a gather of 4*nx1 reals per column and a scatter of nx1).
//!
//! Cost: the root sweeps part_nblk*nx1 rows serially per column.  That is acceptable
//! here (the columns are independent and the kernel is parallel over (root,k,j), and
//! these columns are 84-512 cells long), but it is the serial bottleneck of the scheme;
//! the scalable successor is Schur condensation to the block-interface unknowns (one
//! reduced tridiagonal system of part_nblk rows per column, solved after two local
//! sweeps) or cyclic reduction, neither of which can be bitwise.

void RadiationM1::ImplicitGatherSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nx1 = indcs.nx1, nx1g = part_nx1g;
  int nj = je - js + 1, nk = ke - ks + 1;
  auto iw_ = iw;
  auto sys_ = part_sys;
  auto pos_ = part_pos;
  auto slot_ = part_slot;

  // (1) same-rank members: copy their rows straight into the root's gathered system
  par_for("m1_impl_gth_loc", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    sys_(sl,n,k,j,pos_.d_view(m)*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
  });

#if MPI_PARALLEL_ENABLED
  int nrow = 4*nk*nj*nx1;
  std::vector<MPI_Request> req;
  if (part_any_mpi) {
    auto sb_ = part_sbuf;
    par_for("m1_impl_gth_pack", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      sb_(m,((n*nk + (k-ks))*nj + (j-js))*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
    });
    Kokkos::deep_copy(part_sbuf_h, part_sbuf);
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(&part_rbuf_h(s*part_nblk+p,0), nrow, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]), MPI_COMM_WORLD, &req.back());
      }
    }
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_sbuf_h(m,0), nrow, MPI_ATHENA_REAL, part_rootrank[m],
                2*m, MPI_COMM_WORLD, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    req.clear();
    if (part_nroot > 0) {
      Kokkos::deep_copy(part_rbuf, part_rbuf_h);
      auto rb_ = part_rbuf;
      // which (slot,position) pairs are remote: encoded as a host loop over kernels
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_gth_unp", DevExeSpace(), 0, 3, ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int n, const int k, const int j, const int i) {
            sys_(ss,n,k,j,pp*nx1 + i) =
                rb_(row,((n*nk + (k-ks))*nj + (j-js))*nx1 + i);
          });
        }
      }
    }
  }
#endif

  // (2) the gathered Thomas sweep: the same arithmetic, over nx1g rows
  if (part_nroot > 0) {
    par_for("m1_impl_gth_thomas", DevExeSpace(), 0, part_nroot-1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int s, const int k, const int j) {
      Real bet = sys_(s,1,k,j,0);
      sys_(s,5,k,j,0) = sys_(s,3,k,j,0)/bet;
      for (int i=1; i<nx1g; ++i) {
        sys_(s,4,k,j,i) = sys_(s,2,k,j,i-1)/bet;
        bet = sys_(s,1,k,j,i) - sys_(s,0,k,j,i)*sys_(s,4,k,j,i);
        sys_(s,5,k,j,i) = (sys_(s,3,k,j,i) - sys_(s,0,k,j,i)*sys_(s,5,k,j,i-1))/bet;
      }
      for (int i=nx1g-2; i>=0; --i) {
        sys_(s,5,k,j,i) -= sys_(s,4,k,j,i+1)*sys_(s,5,k,j,i+1);
      }
    });
  }

  // (3) scatter
#if MPI_PARALLEL_ENABLED
  if (part_any_mpi) {
    int nsol = nk*nj*nx1;
    if (part_nroot > 0) {
      auto rb_ = part_rbuf;
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_sct_pack", DevExeSpace(), ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int k, const int j, const int i) {
            rb_(row,((k-ks)*nj + (j-js))*nx1 + i) = sys_(ss,5,k,j,pp*nx1 + i);
          });
        }
      }
      Kokkos::deep_copy(part_rbuf_h, part_rbuf);
    }
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_sbuf_h(m,0), nsol, MPI_ATHENA_REAL, part_rootrank[m],
                2*m + 1, MPI_COMM_WORLD, &req.back());
    }
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&part_rbuf_h(s*part_nblk+p,0), nsol, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]) + 1, MPI_COMM_WORLD,
                  &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    Kokkos::deep_copy(part_sbuf, part_sbuf_h);
    auto sb_ = part_sbuf;
    par_for("m1_impl_sct_unp", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      iw_(m,M1_IW_S2,k,j,i) = sb_(m,((k-ks)*nj + (j-js))*(ie-is+1) + (i-is));
    });
  }
#endif
  par_for("m1_impl_sct_loc", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    iw_(m,M1_IW_S2,k,j,i) = sys_(sl,5,k,j,pos_.d_view(m)*(ie-is+1) + (i-is));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTridiagSolve
//! \brief the x1 LINE SOLVE of the assembled rows: read (M1_IW_TA, TB, TC, TR) and leave
//! the solution in M1_IW_S2.  Plain Thomas, cyclic Thomas (Sherman-Morrison) when the x1
//! boundaries are periodic inside one MeshBlock, or the gathered stack sweep when the
//! column spans several blocks.
//!
//! It is BOTH the line-Jacobi pass (called once per Picard iteration with the lagged
//! right-hand side) and the PRECONDITIONER of the BiCGStab wrapper (called twice per
//! inner iteration with a Krylov vector in M1_IW_TR).  Extracting it changed no
//! arithmetic: the two sweeps below are verbatim what the Picard loop used to run inline.

void RadiationM1::ImplicitTridiagSolve() {
  if (part_nblk > 1) {
    ImplicitGatherSolve();
    return;
  }
  if (!impl_pcr_check) {
    if (impl_line_solver == 1) {
      ImplicitPCRSolve();
    } else {
      ImplicitThomasSolve();
    }
    return;
  }
  // implicit_pcr_check: run the OTHER solver first, keep its answer, then the selected
  // one (whose answer stays in M1_IW_S2), and record max|dx|/max|x| over the pack.
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  if (pcr_chk.extent(0) == 0) {
    Kokkos::realloc(pcr_chk, nmb1+1, indcs.nx3 + 2*indcs.ng*(indcs.nx3 > 1 ? 1 : 0),
                    indcs.nx2 + 2*indcs.ng*(indcs.nx2 > 1 ? 1 : 0), indcs.nx1+2*indcs.ng);
  }
  auto ck_ = pcr_chk;
  if (impl_line_solver == 1) {
    ImplicitThomasSolve();
  } else {
    ImplicitPCRSolve();
  }
  par_for("m1_impl_pcrck_cp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    ck_(m,k,j,i) = iw_(m,M1_IW_S2,k,j,i);
  });
  if (impl_line_solver == 1) {
    ImplicitPCRSolve();
  } else {
    ImplicitThomasSolve();
  }
  const int nkji = (ke-ks+1)*(je-js+1)*(ie-is+1);
  const int nji = (je-js+1)*(ie-is+1), ni = ie-is+1;
  Real dmax = 0.0, xmax = 0.0;
  Kokkos::parallel_reduce("m1_impl_pcrck_red",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dm, Real &xm) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + is;
    k += ks;
    j += js;
    dm = fmax(dm, fabs(iw_(m,M1_IW_S2,k,j,i) - ck_(m,k,j,i)));
    xm = fmax(xm, fabs(ck_(m,k,j,i)));
  }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(xmax));
  Real rel = (xmax > 0.0) ? dmax/xmax : dmax;
  pcr_chk_max = fmax(pcr_chk_max, rel);
  pcr_chk_n += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitThomasSolve
//! \brief the Thomas / cyclic-Thomas line solve, one thread per (m,k,j) column (the
//! original ImplicitTridiagSolve body, unchanged)

void RadiationM1::ImplicitThomasSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  par_for("m1_impl_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (!cyclic) {
      Real bet = iw_(m,M1_IW_TB,k,j,is);
      iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
      for (int i=is+1; i<=ie; ++i) {
        iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
        bet = iw_(m,M1_IW_TB,k,j,i) - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
        iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
      }
      for (int i=ie-1; i>=is; --i) {
        iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
      }
    } else {
      // cyclic tridiagonal, Sherman-Morrison (Press et al. `cyclic`).  alpha is the
      // BOTTOM-LEFT corner, c(ie) (row ie coupling to cell is), and beta the TOP-RIGHT
      // one, a(is) (row is coupling to cell ie) -- that is the convention u and v below
      // are built for, and swapping them is invisible on a symmetric matrix but wrong
      // as soon as upwind advection makes the two off-diagonals differ: it produced a
      // spurious dipole across the seam (-14 % in the first cell, +11 % in the last)
      // on the very first step of T4b, where the exact answer is "nothing moves".
      Real alpha = iw_(m,M1_IW_TC,k,j,ie);
      Real beta = iw_(m,M1_IW_TA,k,j,is);
      Real gam = -iw_(m,M1_IW_TB,k,j,is);
      Real bb0 = iw_(m,M1_IW_TB,k,j,is) - gam;
      Real bbn = iw_(m,M1_IW_TB,k,j,ie) - alpha*beta/gam;
      // solve A' y = r and A' z = u with u = (gam,0,...,0,alpha)
      Real bet = bb0;
      iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
      iw_(m,M1_IW_S3,k,j,is) = gam/bet;
      for (int i=is+1; i<=ie; ++i) {
        Real bd = (i == ie) ? bbn : iw_(m,M1_IW_TB,k,j,i);
        iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
        bet = bd - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
        iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
        Real uu = (i == ie) ? alpha : 0.0;
        iw_(m,M1_IW_S3,k,j,i) = (uu
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S3,k,j,i-1))/bet;
      }
      for (int i=ie-1; i>=is; --i) {
        iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
        iw_(m,M1_IW_S3,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S3,k,j,i+1);
      }
      // x = y - z (v.y)/(1 + v.z),  v = (1,0,...,0,beta/gam)
      Real vy = iw_(m,M1_IW_S2,k,j,is) + (beta/gam)*iw_(m,M1_IW_S2,k,j,ie);
      Real vz = iw_(m,M1_IW_S3,k,j,is) + (beta/gam)*iw_(m,M1_IW_S3,k,j,ie);
      Real fac = vy/(1.0 + vz);
      for (int i=is; i<=ie; ++i) {
        iw_(m,M1_IW_S2,k,j,i) -= fac*iw_(m,M1_IW_S3,k,j,i);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPCRSolve
//! \brief implicit_line_solver = pcr: the same x1 line solve as the Thomas sweep of
//! ImplicitTridiagSolve (non-cyclic, or cyclic by the same Sherman-Morrison split), by
//! PARALLEL CYCLIC REDUCTION with one Kokkos team per (m,k,j) column.  The rows are
//! loaded into team scratch with the team's threads running along i (coalesced), then
//! ceil(log2 nx1) PCR rounds (double-buffered, one team barrier each) decouple every
//! row, and x_i = r_i/b_i.  Out-of-range neighbours of a round are the identity row.
//! Under the cyclic split the second right-hand side u = (gam,0,...,0,alpha) rides the
//! same elimination.  Writes only M1_IW_S2 (the Thomas scratch S1/S3 is not touched and
//! nothing else reads it).  O(n log n) work, so on a CPU (team size 1) it is slower than
//! Thomas; it is meant for the GPU, where the Thomas sweep has one thread per column.

void RadiationM1::ImplicitPCRSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = ie - is + 1;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const int nkj = nk*nj;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int nv = cyclic ? 5 : 4;   // a, b, c, r (+ u) per buffer
  size_t scr_size = ScrArray1D<Real>::shmem_size(2*nv*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  Kokkos::TeamPolicy<DevExeSpace> policy;
  if (!std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    // a GPU: an explicit team size (the threads run along i)
    int ts = impl_pcr_team;
    if (ts == 0) {
      ts = 1;
      while (ts < nx && ts < 256) ts *= 2;
    }
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, ts);
  } else {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, Kokkos::AUTO);
  }
  // level-0 team scratch unless the device cannot fit it (A100: 48 KB per block)
  const int lv = TeamScratchLevel(scr_size, 0, policy.team_size());
  Kokkos::parallel_for("m1_impl_pcr",
                       policy.set_scratch_size(lv, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/nkj;
    const int k = (tm.league_rank() - m*nkj)/nj + ks;
    const int j = (tm.league_rank() - m*nkj)%nj + js;
    ScrArray1D<Real> sw(tm.team_scratch(lv), 2*nv*nx);
    // buffer q (0/1), variable v: sw(q*nv*nx + v*nx + i); v = 0 a, 1 b, 2 c, 3 r, 4 u
    Real alpha = 0.0, beta = 0.0, gam = 1.0;
    if (cyclic) {
      alpha = iw_(m,M1_IW_TC,k,j,ie);
      beta = iw_(m,M1_IW_TA,k,j,is);
      gam = -iw_(m,M1_IW_TB,k,j,is);
    }
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      const int ii = i + is;
      Real bd = iw_(m,M1_IW_TB,k,j,ii);
      if (cyclic) {
        if (i == 0) bd -= gam;
        if (i == nx-1) bd -= alpha*beta/gam;
        sw(4*nx + i) = (i == 0) ? gam : ((i == nx-1) ? alpha : 0.0);
      }
      sw(i) = (i == 0) ? 0.0 : iw_(m,M1_IW_TA,k,j,ii);
      sw(nx + i) = bd;
      sw(2*nx + i) = (i == nx-1) ? 0.0 : iw_(m,M1_IW_TC,k,j,ii);
      sw(3*nx + i) = iw_(m,M1_IW_TR,k,j,ii);
    });
    tm.team_barrier();
    int src = 0;
    for (int rd=0, s=1; rd<nround; ++rd, s*=2) {
      const int o = src*nv*nx, d = (1-src)*nv*nx;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        const int im = i - s, ip = i + s;
        Real ai = sw(o + i), bi = sw(o + nx + i), ci = sw(o + 2*nx + i);
        Real ri = sw(o + 3*nx + i);
        Real ui = cyclic ? sw(o + 4*nx + i) : 0.0;
        Real an = 0.0, cn = 0.0;
        if (im >= 0) {
          Real f = -ai/sw(o + nx + im);
          an = f*sw(o + im);
          bi += f*sw(o + 2*nx + im);
          ri += f*sw(o + 3*nx + im);
          if (cyclic) ui += f*sw(o + 4*nx + im);
        }
        if (ip < nx) {
          Real g = -ci/sw(o + nx + ip);
          cn = g*sw(o + 2*nx + ip);
          bi += g*sw(o + ip);
          ri += g*sw(o + 3*nx + ip);
          if (cyclic) ui += g*sw(o + 4*nx + ip);
        }
        sw(d + i) = an;
        sw(d + nx + i) = bi;
        sw(d + 2*nx + i) = cn;
        sw(d + 3*nx + i) = ri;
        if (cyclic) sw(d + 4*nx + i) = ui;
      });
      tm.team_barrier();
      src = 1 - src;
    }
    const int o = src*nv*nx;
    if (!cyclic) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        iw_(m,M1_IW_S2,k,j,i+is) = sw(o + 3*nx + i)/sw(o + nx + i);
      });
    } else {
      // x = y - z (v.y)/(1 + v.z),  v = (1,0,...,0,beta/gam)
      Real y0 = sw(o + 3*nx)/sw(o + nx);
      Real yn = sw(o + 4*nx - 1)/sw(o + 2*nx - 1);
      Real z0 = sw(o + 4*nx)/sw(o + nx);
      Real zn = sw(o + 5*nx - 1)/sw(o + 2*nx - 1);
      Real fac = (y0 + (beta/gam)*yn)/(1.0 + z0 + (beta/gam)*zn);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        Real bi = sw(o + nx + i);
        iw_(m,M1_IW_S2,k,j,i+is) = sw(o + 3*nx + i)/bi - fac*(sw(o + 4*nx + i)/bi);
      });
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPCRSolveX
//! \brief implicit_krylov_fuse >= 1: the pcr line solve of ImplicitPCRSolve, term for
//! term, with the right-hand side taken from component `rc` of iw (upd = 0) or MADE in
//! the load phase as the BiCGStab update p = r + c1 (p - c2 v) (upd = 1, written to
//! M1_IW_KP) or s = r - c1 v (upd = 2, written to M1_IW_KS), and the answer written
//! straight into component `zc`.  That removes the staging copy into M1_IW_TR, the copy
//! out of M1_IW_S2 and the separate p / s kernels: 1 launch where there were 3.
//!
//! implicit_precond = rbgs: `col` = 0 / 1 solves only the columns (k,j) whose parity
//! (k-ks)+(j-js) is `col` (one launch per colour, half the teams), and `sub` >= 0
//! subtracts the transverse 5-point coupling sum_nb C_nb z_nb of component `sub` from
//! the right-hand side, over the neighbours INSIDE the MeshBlock only (the other colour,
//! which the previous half-sweep has just written).  col < 0: every column, sub ignored.

namespace {
//! the kernel of ImplicitPCRSolveX, with the elimination carried in T (Real, or float
//! under implicit_precond_float: a preconditioner only has to be a FIXED linear map,
//! so its precision sets the convergence rate, not the converged answer)
template <typename T>
void M1PCRX(const DvceArray5D<Real> &iw_, Kokkos::TeamPolicy<DevExeSpace> policy,
            const int is, const int ie, const int js, const int je, const int ks,
            const int ke, const int nkj, const int njl, const bool colr, const int cl_,
            const int cs_, const bool thrd, const bool cyclic, const int cr,
            const int cz, const int up, const Real a1, const Real a2) {
  const int nx = ie - is + 1;
  const int nv = cyclic ? 5 : 4;   // a, b, c, r (+ u) per buffer
  size_t scr_size = ScrArray1D<T>::shmem_size(2*nv*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  const int lv = TeamScratchLevel(scr_size, 0, policy.team_size());
  Kokkos::parallel_for("m1_impl_pcrx",
                       policy.set_scratch_size(lv, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/nkj;
    const int kk = (tm.league_rank() - m*nkj)/njl;
    const int jj = (tm.league_rank() - m*nkj)%njl;
    const int k = kk + ks;
    const int j = colr ? (js + 2*jj + ((cl_ + kk) & 1)) : (jj + js);
    if (j > je) return;   // team-uniform
    ScrArray1D<T> sw(tm.team_scratch(lv), 2*nv*nx);
    Real alpha = 0.0, beta = 0.0, gam = 1.0;
    if (cyclic) {
      alpha = iw_(m,M1_IW_TC,k,j,ie);
      beta = iw_(m,M1_IW_TA,k,j,is);
      gam = -iw_(m,M1_IW_TB,k,j,is);
    }
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      const int ii = i + is;
      Real bd = iw_(m,M1_IW_TB,k,j,ii);
      if (cyclic) {
        if (i == 0) bd -= gam;
        if (i == nx-1) bd -= alpha*beta/gam;
        sw(4*nx + i) = static_cast<T>((i == 0) ? gam : ((i == nx-1) ? alpha : 0.0));
      }
      sw(i) = static_cast<T>((i == 0) ? 0.0 : iw_(m,M1_IW_TA,k,j,ii));
      sw(nx + i) = static_cast<T>(bd);
      sw(2*nx + i) = static_cast<T>((i == nx-1) ? 0.0 : iw_(m,M1_IW_TC,k,j,ii));
      Real rr;
      if (up == 1) {
        rr = iw_(m,M1_IW_KR,k,j,ii)
             + a1*(iw_(m,M1_IW_KP,k,j,ii) - a2*iw_(m,M1_IW_KV,k,j,ii));
        iw_(m,M1_IW_KP,k,j,ii) = rr;
      } else if (up == 2) {
        rr = iw_(m,M1_IW_KR,k,j,ii) - a1*iw_(m,M1_IW_KV,k,j,ii);
        iw_(m,M1_IW_KS,k,j,ii) = rr;
      } else {
        rr = iw_(m,cr,k,j,ii);
      }
      if (cs_ >= 0) {
        if (j > js) rr -= iw_(m,M1_IW_CJM,k,j,ii)*iw_(m,cs_,k,j-1,ii);
        if (j < je) rr -= iw_(m,M1_IW_CJP,k,j,ii)*iw_(m,cs_,k,j+1,ii);
        if (thrd) {
          if (k > ks) rr -= iw_(m,M1_IW_CKM,k,j,ii)*iw_(m,cs_,k-1,j,ii);
          if (k < ke) rr -= iw_(m,M1_IW_CKP,k,j,ii)*iw_(m,cs_,k+1,j,ii);
        }
      }
      sw(3*nx + i) = static_cast<T>(rr);
    });
    tm.team_barrier();
    int src = 0;
    for (int rd=0, s=1; rd<nround; ++rd, s*=2) {
      const int o = src*nv*nx, d = (1-src)*nv*nx;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        const int im = i - s, ip = i + s;
        T ai = sw(o + i), bi = sw(o + nx + i), ci = sw(o + 2*nx + i);
        T ri = sw(o + 3*nx + i);
        T ui = cyclic ? sw(o + 4*nx + i) : static_cast<T>(0.0);
        T an = 0.0, cn = 0.0;
        if (im >= 0) {
          T f = -ai/sw(o + nx + im);
          an = f*sw(o + im);
          bi += f*sw(o + 2*nx + im);
          ri += f*sw(o + 3*nx + im);
          if (cyclic) ui += f*sw(o + 4*nx + im);
        }
        if (ip < nx) {
          T g = -ci/sw(o + nx + ip);
          cn = g*sw(o + 2*nx + ip);
          bi += g*sw(o + ip);
          ri += g*sw(o + 3*nx + ip);
          if (cyclic) ui += g*sw(o + 4*nx + ip);
        }
        sw(d + i) = an;
        sw(d + nx + i) = bi;
        sw(d + 2*nx + i) = cn;
        sw(d + 3*nx + i) = ri;
        if (cyclic) sw(d + 4*nx + i) = ui;
      });
      tm.team_barrier();
      src = 1 - src;
    }
    const int o = src*nv*nx;
    if (!cyclic) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        iw_(m,cz,k,j,i+is) = static_cast<Real>(sw(o + 3*nx + i)/sw(o + nx + i));
      });
    } else {
      T y0 = sw(o + 3*nx)/sw(o + nx);
      T yn = sw(o + 4*nx - 1)/sw(o + 2*nx - 1);
      T z0 = sw(o + 4*nx)/sw(o + nx);
      T zn = sw(o + 5*nx - 1)/sw(o + 2*nx - 1);
      T bg = static_cast<T>(beta/gam);
      T fac = (y0 + bg*yn)/(static_cast<T>(1.0) + z0 + bg*zn);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        T bi = sw(o + nx + i);
        iw_(m,cz,k,j,i+is) = static_cast<Real>(sw(o + 3*nx + i)/bi
                                               - fac*(sw(o + 4*nx + i)/bi));
      });
    }
  });
}
} // namespace

void RadiationM1::ImplicitPCRSolveX(int rc, int zc, int upd, Real c1, Real c2, int col,
                                    int sub) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = ie - is + 1;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const bool colr = (col >= 0);
  const int njl = colr ? (nj + 1)/2 : nj;   // league columns per (m,k)
  const int nkj = nk*njl;
  const int cl_ = colr ? col : 0;
  const int cs_ = (colr && sub >= 0) ? sub : -1;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  Kokkos::TeamPolicy<DevExeSpace> policy;
  if (!std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    int ts = impl_pcr_team;
    if (ts == 0) {
      ts = 1;
      while (ts < nx && ts < 256) ts *= 2;
    }
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, ts);
  } else {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, Kokkos::AUTO);
  }
  if (impl_prec_float) {
    M1PCRX<float>(iw, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                  trans_x3, cyclic, rc, zc, upd, c1, c2);
  } else {
    M1PCRX<Real>(iw, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                 trans_x3, cyclic, rc, zc, upd, c1, c2);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPrecondX
//! \brief the preconditioner of the implicit_krylov_fuse path, z = M^{-1} r with r from
//! `rc` or made by the p / s update `upd` (ImplicitPCRSolveX):
//!  implicit_precond = line: M = the x1 line part of the row (as ImplicitPrecond);
//!  implicit_precond = rbgs: ONE symmetric red-black transverse line Gauss-Seidel sweep
//!    from z = 0 (red, black, red; each colour one line solve of its columns against
//!    the just-updated other colour), block-local: the x2/x3 couplings across a
//!    MeshBlock face are left out of M (they stay in the operator), so it needs no
//!    communication.  M is symmetric positive definite whenever the 5-point part is, and
//!    stronger than line Jacobi on the transverse coupling (an approximate inverse of
//!    the whole 5-point-per-line system instead of its x1 part alone).
//!  implicit_precond = rbgs_fwd: the forward half only (red, black).

void RadiationM1::ImplicitPrecondX(int rc, int zc, int upd, Real c1, Real c2) {
  if (impl_prec >= 3) {
    ImplicitMGApply(rc, zc, upd, c1, c2);
    return;
  }
  if (impl_prec == 0) {
    ImplicitPCRSolveX(rc, zc, upd, c1, c2, -1, -1);
    return;
  }
  // the right-hand side the third half-sweep re-reads: the update the first one wrote
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 0, -1);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 1, zc);
  if (impl_prec == 1) {
    ImplicitPCRSolveX(rr, zc, 0, 0.0, 0.0, 0, zc);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitODCache
//! \brief implicit_od_cache: odc(m,d,k,j,i) = M1OffDiv(x, d) at every cell a face of
//! ImplicitOffDiagOp reads it from -- the active box plus one layer in x1, x2 (and x3)
//! -- with exactly the index limits that routine passes, so every face value
//! 0.5*(odc_L + odc_R) is the very number ImplicitOffDiagOp forms (bitwise), from 3
//! evaluations per cell instead of 12.  Cells the faces never read are computed and
//! ignored.

void RadiationM1::ImplicitODCache(int xc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const bool thrd = trans_x3;
  const int e3 = thrd ? 1 : 0;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto od_ = odc;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = (cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs).d_view;   // CS1: seams open
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  const int cx = xc;
  par_for("m1_impl_odc", DevExeSpace(), 0, nmb1, ks-e3, ke+e3, js-1, je+1, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dx1 = mbsize(m).dx1;
    Real dx2 = mbsize(m).dx2;
    Real dx3 = mbsize(m).dx3;
    BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
    BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
    if (!p2lo) {jl = js-1;}
    if (!p2hi) {ju = je+1;}
    if (thrd && !p3lo) {kl = ks-1;}
    if (thrd && !p3hi) {ku = ke+1;}
    od_(m,0,k,j,i) = M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull);
    od_(m,1,k,j,i) = M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull);
    if (thrd) {
      od_(m,2,k,j,i) = M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirectInit
//! \brief implicit_halo_direct: tabulate, for every MeshBlock of the pack and each of the
//! 26 directions (ox1,ox2,ox3), the LOCAL index of the same-level neighbour that fills
//! that ghost region, or -1 where there is none (a physical boundary: its ghost zones
//! are not filled by the ordinary exchange either).  The direct copy is used only when
//! EVERY rank finds every neighbour on its own rank at the same level, outside the
//! cubed-sphere and polar transforms (so all ranks take the same branch and no MPI
//! exchange is left half-posted); otherwise the ordinary exchange runs, as before.

void RadiationM1::ImplicitHaloDirectInit() {
  auto *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack;
  hd_src = DualArray2D<int>("m1_hd_src", nmb, 27);
  int ok = (pm->multilevel || pm->use_cubed_sphere || pm->use_polar_boundary) ? 0 : 1;
  auto &nb = pmy_pack->pmb->nghbr;
  auto &lev = pmy_pack->pmb->mb_lev;
  const int e2 = pm->multi_d ? 1 : 0;
  const int e3 = pm->three_d ? 1 : 0;
  for (int m = 0; m < nmb; ++m) {
    for (int d = 0; d < 27; ++d) {hd_src.h_view(m,d) = -1;}
    for (int o3 = -e3; o3 <= e3; ++o3) {
      for (int o2 = -e2; o2 <= e2; ++o2) {
        for (int o1 = -1; o1 <= 1; ++o1) {
          if (o1 == 0 && o2 == 0 && o3 == 0) continue;
          int n = NeighborIndex(o1, o2, o3, 0, 0);
          if (n < 0 || n >= pmy_pack->pmb->nnghbr) {ok = 0; continue;}
          const NeighborBlock &q = nb.h_view(m,n);
          if (q.gid < 0) continue;
          if (q.rank != global_variable::my_rank || q.lev != lev.h_view(m)) {
            ok = 0;
            continue;
          }
          hd_src.h_view(m, (o1+1) + 3*(o2+1) + 9*(o3+1)) = q.gid - pmy_pack->gids;
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  {int g = ok;
  MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  ok = g;}
#endif
  halo_direct_on = (ok == 1);
  hd_src.modify_host();
  hd_src.sync_device();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirect
//! \brief implicit_halo_direct: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list) filled by ONE kernel that copies each ghost
//! cell from the active cell of the neighbour that owns it (all neighbours are on this
//! rank at the same level, ImplicitHaloDirectInit).  The same numbers the pack /
//! exchange / unpack chain delivers: 1 launch where there were 4, and no host work.

void RadiationM1::ImplicitHaloDirect(int nq, int c0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int n1 = nx1 + 2*indcs.ng;
  const int n2 = (nx2 > 1) ? (nx2 + 2*indcs.ng) : 1;
  const int n3 = (nx3 > 1) ? (nx3 + 2*indcs.ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool md = (nx2 > 1), td = (nx3 > 1);
  auto iw_ = iw;
  auto tab = hd_src.d_view;
  const int nc0 = c0;
  par_for("m1_impl_hdir", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    const int o1 = (i < is) ? -1 : ((i > ie) ? 1 : 0);
    const int o2 = md ? ((j < js) ? -1 : ((j > je) ? 1 : 0)) : 0;
    const int o3 = td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    if (o1 == 0 && o2 == 0 && o3 == 0) return;
    const int src = tab(m, (o1+1) + 3*(o2+1) + 9*(o3+1));
    if (src < 0) return;
    const int nc = (nc0 >= 0) ? (nc0 + n) : M1HaloCompT(n);
    iw_(m,nc,k,j,i) = iw_(src,nc,k - o3*nx3,j - o2*nx2,i - o1*nx1);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitKrylovHalo
//! \brief milestone 3b phase C: put ONE component of the work array into the scratch
//! array `krw`, exchange it with all six neighbours through the module's ordinary
//! cell-centred boundary machinery (which is what supplies periodic wrap, edge/corner
//! neighbours and MPI), and copy the ghost zones back.  PHYSICAL boundaries are not
//! filled: the row's coefficient towards such a neighbour is identically zero, so the
//! ghost value is multiplied by zero and never read in anger.

void RadiationM1::ImplicitKrylovHalo(int comp) {
  ImplicitHaloExchange(1, comp);
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1OdFaces
//! \brief implicit_od_cache: the six face terms of ImplicitOffDiagOp at cell (m,k,j,i),
//! in the same order and with the same arithmetic, from the per-cell cache od(m,d,...)
//! = M1OffDiv(x, d) instead of re-evaluating M1OffDiv at both cells of every face.
//! The caller has already returned on an M1_IBC_EFIX row.

KOKKOS_INLINE_FUNCTION
Real M1OdFaces(const DvceArray5D<Real> &iw_, const DvceArray5D<Real> &od_,
               const DvceArray4D<Real> &th2_, const DvceArray4D<Real> &th3_,
               const bool lm, const int m, const int k, const int j, const int i,
               const int is, const int ie, const int js, const int je, const int ks,
               const int ke, const bool cyclic, const bool botb, const bool topb,
               const bool p2lo, const bool p2hi, const bool p3lo, const bool p3hi,
               const bool thrd, const Real dx1, const Real dx2, const Real dx3,
               const Real ch, const Real cl, const Real dt, const bool bx1,
               const DvceArray5D<Real> &fw_) {
  Real cr = ch/cl;
  Real kk = ch*cl*dt;
  Real y = 0.0;
  if (i < ie || cyclic || !topb) {
    int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    if (bx1) {th *= 1.0 - fw_(m,M1_IFW_AL,k,j,i+1);}   // blendall-1009
    Real od = 0.5*(od_(m,0,k,j,i) + od_(m,0,k,j,ip));
    y -= (dt/dx1)*cr*th*kk*od;
  }
  if (i > is || cyclic || !botb) {
    int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    if (bx1) {th *= 1.0 - fw_(m,M1_IFW_AL,k,j,i);}
    Real od = 0.5*(od_(m,0,k,j,im) + od_(m,0,k,j,i));
    y += (dt/dx1)*cr*th*kk*od;
  }
  Real nu2 = dt/dx2;
  if (!(j == je && p2hi)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
    Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,1,k,j,i) + od_(m,1,k,j+1,i));
    y -= nu2*cr*th*kk*od;
  }
  if (!(j == js && p2lo)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,1,k,j-1,i) + od_(m,1,k,j,i));
    y += nu2*cr*th*kk*od;
  }
  if (thrd) {
    Real nu3 = dt/dx3;
    if (!(k == ke && p3hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
      Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(od_(m,2,k,j,i) + od_(m,2,k+1,j,i));
      y -= nu3*cr*th*kk*od;
    }
    if (!(k == ks && p3lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(od_(m,2,k-1,j,i) + od_(m,2,k,j,i));
      y += nu3*cr*th*kk*od;
    }
  }
  return y;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOffDiagOpC
//! \brief implicit_od_cache: ImplicitOffDiagOp (y += sgn L_off(x)) from the cache, and,
//! with `with7`, the whole operator y = A x in ONE kernel (the 7-point row of
//! ImplicitApplyOp, then + L_off(x) exactly as ImplicitApplyOp + ImplicitOffDiagOp
//! form it).  The reduction rides in the operator kernel (implicit_krylov_fuse >= 2):
//! red = 1: out[0] = (rhat,y);  red = 2: out[0] = (y,s), out[1] = (y,y);
//! red = 3: as 2 plus out[2] = (rhat,y);  red = 4: as 1 plus out[3] = max|r|.

void RadiationM1::ImplicitOffDiagOpC(int xc, int yc, Real sgn, bool with7, int red,
                                     Real *out) {
  ImplicitODCache(xc);
  const bool bx1 = bvec_x1 && trans_on;   // blendall-1009: x1 faces keep 1 - AL
  auto fw_ = ifw;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto od_ = odc;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = (cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs).d_view;   // CS1: seams open
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  const bool lm = (impl_tlim != M1_TLIM_NONE) || blat_on;
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const int cx = xc, cy = yc;
  const Real sg = sgn;
  const bool w7 = with7;
  const int rm = red;
  // the standalone call (with7 = false) is made only under the operator form
  const bool odon = w7 ? (od_now == M1_OD_OPERATOR) : true;
  const bool vim = w7 && vimp_now;
  const int ivb = iw_vimp;
  const bool mus = w7 && muscl_now;   // xthinfix-1009 Fix B
  const int imb = iw_muscl;
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  // one cell's result; returns y (the value stored)
  auto row = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) -> Real {
    Real y7 = 0.0;
    if (w7) {
      int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
      int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
      y7 = iw_(m,M1_IW_TB,k,j,i)*iw_(m,cx,k,j,i)
           + iw_(m,M1_IW_TA,k,j,i)*iw_(m,cx,k,j,im)
           + iw_(m,M1_IW_TC,k,j,i)*iw_(m,cx,k,j,ip)
           + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,cx,k,j-1,i)
           + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,cx,k,j+1,i);
      if (thrd) {
        y7 += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,cx,k-1,j,i)
              + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,cx,k+1,j,i);
      }
      if (vim) {y7 += M1VimpRow(iw_, ivb, cx, m, k, j, i, is, ie, cyclic, thrd);}
      if (mus) {y7 += M1MusclRow(iw_, imb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    } else {
      y7 = iw_(m,cy,k,j,i);
    }
    if (!odon) {
      iw_(m,cy,k,j,i) = y7;
      return y7;
    }
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                    (i == ie && topb && bchi == M1_IBC_EFIX))) {
      iw_(m,cy,k,j,i) = y7;
      return y7;
    }
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    Real y = M1OdFaces(iw_, od_, th2_, th3_, lm, m, k, j, i, is, ie, js, je, ks, ke,
                       cyclic, botb, topb, p2lo, p2hi, p3lo, p3hi, thrd,
                       mbsize(m).dx1, mbsize(m).dx2, mbsize(m).dx3, ch, cl, dt,
                       bx1, fw_);
    Real out = y7 + sg*y;
    iw_(m,cy,k,j,i) = out;
    return out;
  };
  if (rm == 0) {
    par_for("m1_impl_opc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      row(m, k, j, i);
    });
    return;
  }
  // 256-thread blocks: the default 1024-thread block of a reduction with a 4-Real
  // value takes 33 kB of LDS, i.e. ONE block per CU (measured 97 us vs 47 us for the
  // same stencil as a par_for)
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  Real a0 = 0.0, a1 = 0.0, a2 = 0.0, amx = 0.0;
  Kokkos::parallel_reduce("m1_impl_opcr", pol,
  KOKKOS_LAMBDA(const int idx, Real &l0, Real &l1, Real &l2, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real y = row(m, k, j, i);
    if (rm == 1 || rm == 4) {
      l0 += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        lmx = (a > lmx) ? a : lmx;
      }
    } else {
      l0 += y*iw_(m,M1_IW_KS,k,j,i);
      l1 += y*y;
      if (rm == 3) {l2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, a0, a1, a2, Kokkos::Max<Real>(amx));
  out[0] = a0;
  out[1] = a1;
  out[2] = a2;
  out[3] = amx;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1StIdx
//! \brief implicit_op_stencil: the slot of the 19-point stencil for the offset
//! (di,dj,dk) in {-1,0,1}^3 with at most two non-zero entries: 0 centre; 1..6 the faces
//! i-1,i+1,j-1,j+1,k-1,k+1; 7..10 the (i,j) edges, 11..14 the (i,k) edges, 15..18 the
//! (j,k) edges, each ordered (-,-),(+,-),(-,+),(+,+) in (first, second) axis.

KOKKOS_INLINE_FUNCTION
constexpr int M1StIdx(const int di, const int dj, const int dk) {
  if (dk == 0) {
    if (dj == 0) {return (di == 0) ? 0 : ((di < 0) ? 1 : 2);}
    if (di == 0) {return (dj < 0) ? 3 : 4;}
    return 7 + ((di > 0) ? 1 : 0) + ((dj > 0) ? 2 : 0);
  }
  if (dj == 0) {
    if (di == 0) {return (dk < 0) ? 5 : 6;}
    return 11 + ((di > 0) ? 1 : 0) + ((dk > 0) ? 2 : 0);
  }
  return 15 + ((dj > 0) ? 1 : 0) + ((dk > 0) ? 2 : 0);
}

//! D_de of the frozen closure at one cell (the coefficient M1POff multiplies x by)
KOKKOS_INLINE_FUNCTION
Real M1DOffC(const DvceArray5D<Real> &iw, const DvceArray5D<Real> &vd, const bool full,
             const int m, const int a, const int b, const int k, const int j,
             const int i) {
  if (full) {return vd(m,M1_VET_D11+2+a+b,k,j,i);}
  return M1EddOff(iw(m,M1_IW_WCHI,k,j,i), iw(m,M1_IW_N1+a,k,j,i),
                  iw(m,M1_IW_N1+b,k,j,i));
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitStencilBuild
//! \brief implicit_op_stencil: the frozen operator of this Picard pass -- the 7-point
//! row (TA,TB,TC,CJM..CKP) plus, under implicit_offdiag = operator, the off-diagonal
//! Eddington terms of ImplicitOffDiagOp -- written out ONCE as a 19-point stencil
//! st(m,0..18,k,j,i).  Every face term of ImplicitOffDiagOp,
//!   -/+ (dt/dx_d) (chat/c) theta_f chat c dt * 0.5 [OD_d(L) + OD_d(R)],
//!   OD_d(q) = sum_{e!=d} [D_de x](q_a) - [D_de x](q_b)) / ((a-b) dx_e)
//! with the same one-sided clamps at physical faces, is linear in x with coefficients
//! frozen over the pass, so it distributes onto the centre, the 6 face and the 12 edge
//! neighbours.  Applying the stencil (ImplicitStencilOp) is then one read of 19
//! coefficients per cell instead of re-deriving D_de and theta at every face in every
//! Krylov iteration.  Same operator; the sums are grouped differently (round-off).
//! Not for a periodic x1 wrap (cyclic), which keeps the od_cache path.

void RadiationM1::ImplicitStencilBuild() {
  if (ibc_x1min == M1_IBC_PERIODIC) {   // a problem generator may set it late
    ImplFatal("<rad_m1>/implicit_op_stencil does not take a periodic x1 wrap");
  }
  const bool bx1 = bvec_x1 && trans_on;   // blendall-1009: x1 faces keep 1 - AL
  auto fw_ = ifw;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto st_ = ost;
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = (cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs).d_view;   // CS1: seams open
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool thrd = trans_x3;
  const bool lm = (impl_tlim != M1_TLIM_NONE) || blat_on;
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const bool odon = (od_now == M1_OD_OPERATOR) && !impl_odskip;
  // implicit_vimp_fold: the vimp part of the row goes into the stencil (slots 3-6 and
  // 19-24), and ImplicitStencilOp does not call M1VimpRow
  const bool vfold = impl_vfold && vimp_now;
  const int ivb = iw_vimp;
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Real emax = 0.0;
  Kokkos::parallel_reduce("m1_impl_stb",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>(DevExeSpace(), 0,
                                                                (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &lmx) M1_INL {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real c[19];
    for (int o = 0; o < 19; ++o) {c[o] = 0.0;}
    c[0] = iw_(m,M1_IW_TB,k,j,i);
    c[1] = iw_(m,M1_IW_TA,k,j,i);
    c[2] = iw_(m,M1_IW_TC,k,j,i);
    c[3] = iw_(m,M1_IW_CJM,k,j,i);
    c[4] = iw_(m,M1_IW_CJP,k,j,i);
    if (thrd) {
      c[5] = iw_(m,M1_IW_CKM,k,j,i);
      c[6] = iw_(m,M1_IW_CKP,k,j,i);
    }
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    bool efix = (i == is && botb && bclo == M1_IBC_EFIX) ||
                (i == ie && topb && bchi == M1_IBC_EFIX);
    if (odon && !efix) {
      Real dxv[3] = {mbsize(m).dx1, mbsize(m).dx2, mbsize(m).dx3};
      BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
      BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
      BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
      BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
      BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
      BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
      bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
      bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
      bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
      bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
      // the index limits M1OffDiv is given, per axis
      int lo[3] = {is, js, ks}, hi[3] = {ie, je, ke};
      if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {lo[0] = is-1;}
      if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {hi[0] = ie+1;}
      if (!p2lo) {lo[1] = js-1;}
      if (!p2hi) {hi[1] = je+1;}
      if (thrd && !p3lo) {lo[2] = ks-1;}
      if (thrd && !p3hi) {hi[2] = ke+1;}
      const Real cr = ch/cl, kk = ch*cl*dt;
      // m1-fast5-sp: the loops over (face axis d, side sd, cell qs, cross axis e) are
      // unrolled at compile time, so every slot of c is a constant (two candidates per
      // term, chosen by the one-sided clamp): c stays in registers instead of a 176 B
      // per-lane scratch array.  The same terms are added in the same order (bitwise).
      auto qe = [&](auto dc, auto sdc, auto qsc, auto ec, const Real w) M1_INL {
        constexpr int d = decltype(dc)::value, sd = decltype(sdc)::value;
        constexpr int qs = decltype(qsc)::value, e = decltype(ec)::value;
        int q[3] = {i, j, k};
        if (qs == 1) {q[d] += sd;}
        const bool ua = (q[e] + 1 <= hi[e]);
        const bool ub = (q[e] - 1 >= lo[e]);
        if (!ua && !ub) return;
        int qa[3] = {q[0], q[1], q[2]}, qb[3] = {q[0], q[1], q[2]};
        if (ua) {qa[e] = q[e] + 1;}
        if (ub) {qb[e] = q[e] - 1;}
        const Real f = w/((qa[e] - qb[e])*dxv[e]);
        const Real da = M1DOffC(iw_, vd_, dfull, m, d, e, qa[2], qa[1], qa[0]);
        const Real db = M1DOffC(iw_, vd_, dfull, m, d, e, qb[2], qb[1], qb[0]);
        constexpr int o0 = (d == 0) ? qs*sd : 0, o1 = (d == 1) ? qs*sd : 0;
        constexpr int o2 = (d == 2) ? qs*sd : 0;
        constexpr int s0 = M1StIdx(o0, o1, o2);
        constexpr int sa = M1StIdx(o0 + (e == 0), o1 + (e == 1), o2 + (e == 2));
        constexpr int sb = M1StIdx(o0 - (e == 0), o1 - (e == 1), o2 - (e == 2));
        // value selects, not a store to a selected slot (that would index c at run time)
        const Real va = f*da, vb = f*db;
        c[sa] = ua ? (c[sa] + va) : c[sa];
        c[s0] = ua ? c[s0] : (c[s0] + va);
        c[sb] = ub ? (c[sb] - vb) : c[sb];
        c[s0] = ub ? c[s0] : (c[s0] - vb);
      };
      auto face = [&](auto dc, auto sdc) M1_INL {
        constexpr int d = decltype(dc)::value, sd = decltype(sdc)::value;
        // does this face carry the term (the conditions of ImplicitOffDiagOp)?
        bool has;
        if (d == 0) {
          has = (sd > 0) ? (i < ie || !topb) : (i > is || !botb);
        } else if (d == 1) {
          has = (sd > 0) ? !(j == je && p2hi) : !(j == js && p2lo);
        } else {
          has = (sd > 0) ? !(k == ke && p3hi) : !(k == ks && p3lo);
        }
        if (!has) return;
        // the face theta, from the two cells' transport opacities (or the limiter)
        int nb[3] = {i, j, k};
        nb[d] += sd;
        Real th;
        if (d > 0 && lm) {
          th = (d == 1) ? th2_(m,k,(sd > 0) ? j+1 : j,i)
                        : th3_(m,(sd > 0) ? k+1 : k,j,i);
        } else {
          Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,nb[2],nb[1],nb[0]));
          th = 1.0/(1.0 + ch*dt*ktf);
          if (d == 0 && bx1) {   // blendall-1009
            th *= 1.0 - fw_(m,M1_IFW_AL,k,j,(sd > 0) ? i+1 : i);
          }
        }
        // upper face: y -= w od; lower face: y += w od; od = 0.5 (OD(c) + OD(nb))
        const Real w = -static_cast<Real>(sd)*(dt/dxv[d])*cr*th*kk*0.5;
        auto qsall = [&](auto qsc) M1_INL {
          if constexpr (d != 0) {qe(dc, sdc, qsc, std::integral_constant<int, 0>{}, w);}
          if constexpr (d != 1) {qe(dc, sdc, qsc, std::integral_constant<int, 1>{}, w);}
          if constexpr (d != 2) {
            if (thrd) {qe(dc, sdc, qsc, std::integral_constant<int, 2>{}, w);}
          }
        };
        qsall(std::integral_constant<int, 0>{});
        qsall(std::integral_constant<int, 1>{});
      };
      using I0 = std::integral_constant<int, 0>;
      using I1 = std::integral_constant<int, 1>;
      using I2 = std::integral_constant<int, 2>;
      using IM = std::integral_constant<int, -1>;
      face(I0{}, IM{});
      face(I0{}, I1{});
      face(I1{}, IM{});
      face(I1{}, I1{});
      if (thrd) {
        face(I2{}, IM{});
        face(I2{}, I1{});
      }
    }
    for (int o = 0; o < 19; ++o) {st_(m,o,k,j,i) = c[o];}
    for (int o = 7; o < 19; ++o) {lmx = fmax(lmx, fabs(c[o]));}
    if (vfold) {
      st_(m,3,k,j,i) = c[3] + iw_(m,ivb+M1_IV_X2M2+1,k,j,i);
      st_(m,4,k,j,i) = c[4] + iw_(m,ivb+M1_IV_X2M2+2,k,j,i);
      st_(m,19,k,j,i) = iw_(m,ivb+M1_IV_X1M2,k,j,i);
      st_(m,20,k,j,i) = iw_(m,ivb+M1_IV_X1P2,k,j,i);
      st_(m,21,k,j,i) = iw_(m,ivb+M1_IV_X2M2,k,j,i);
      st_(m,22,k,j,i) = iw_(m,ivb+M1_IV_X2M2+3,k,j,i);
      if (thrd) {
        st_(m,5,k,j,i) = c[5] + iw_(m,ivb+M1_IV_X3M2+1,k,j,i);
        st_(m,6,k,j,i) = c[6] + iw_(m,ivb+M1_IV_X3M2+2,k,j,i);
        st_(m,23,k,j,i) = iw_(m,ivb+M1_IV_X3M2,k,j,i);
        st_(m,24,k,j,i) = iw_(m,ivb+M1_IV_X3M2+3,k,j,i);
      }
    }
  }, Kokkos::Max<Real>(emax));
  // no edge coefficient anywhere (the Eddington closure: D_ab = 0 off the diagonal, or
  // implicit_offdiag not operator): ImplicitStencilOp reads the 7 face/centre slots only,
  // which adds the same non-zero terms in the same order
#if MPI_PARALLEL_ENABLED
  {Real g = emax;
  MPI_Allreduce(&emax, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  emax = g;}
#endif
  st_edges = (emax > 0.0);
  // vet_col_lat_offdiag = operator: the lateral off-diagonal terms (rad_m1_vetlat.cpp)
  if (VlatOp() && vlat_ready && sph_geom) {
    VetLatStencilAdd();
    st_edges = true;
  }
}

namespace {
//! implicit_op_team_red: three sums and one max (|r| >= 0: 0 is the max identity), for
//! the team-level reduction of ImplicitStencilOp (the reference is a per-thread value)
struct M1R4Val {
  Real s[3];
  Real mx;
};
struct M1R4Red {
 public:
  using reducer = M1R4Red;
  using value_type = M1R4Val;
  using result_view_type = Kokkos::View<value_type, Kokkos::HostSpace,
                                        Kokkos::MemoryUnmanaged>;

 private:
  result_view_type value;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit M1R4Red(value_type &v) : value(&v) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    for (int q = 0; q < 3; ++q) {d.s[q] += s.s[q];}
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    for (int q = 0; q < 3; ++q) {v.s[q] = 0.0;}
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return true;}
};
// the team size of implicit_op_team_red: 256 cells on a GPU; a host backend allows only
// 1-thread teams (Serial; 256 exceeds OpenMP's limit), so there each cell is its own
// partial (round-off)
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA) \
    || defined(KOKKOS_ENABLE_SYCL)
constexpr int kM1TeamRed = 256;
#else
constexpr int kM1TeamRed = 1;
#endif
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitStencilOp
//! \brief implicit_op_stencil: y = A x from the 19-point stencil of ImplicitStencilBuild
//! (the caller has filled the ghost zones of x).  `red` and `out` as ImplicitOffDiagOpC.

void RadiationM1::ImplicitStencilOp(int xc, int yc, int red, Real *out,
                                    bool red_only) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto st_ = ost;
  const bool thrd = trans_x3;
  const int cx = xc, cy = yc;
  const int rm = red;
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  const bool edg = st_edges;
  const bool vfold = impl_vfold && vimp_now;
  const bool vim = vimp_now && !vfold;
  const int ivb = iw_vimp;
  const bool mus = muscl_now;   // xthinfix-1009 Fix B
  const int imb = iw_muscl;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  auto row = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) -> Real {
    Real y = st_(m,0,k,j,i)*iw_(m,cx,k,j,i)
             + st_(m,1,k,j,i)*iw_(m,cx,k,j,i-1) + st_(m,2,k,j,i)*iw_(m,cx,k,j,i+1)
             + st_(m,3,k,j,i)*iw_(m,cx,k,j-1,i) + st_(m,4,k,j,i)*iw_(m,cx,k,j+1,i);
    if (edg) {
      y += st_(m,7,k,j,i)*iw_(m,cx,k,j-1,i-1) + st_(m,8,k,j,i)*iw_(m,cx,k,j-1,i+1)
           + st_(m,9,k,j,i)*iw_(m,cx,k,j+1,i-1) + st_(m,10,k,j,i)*iw_(m,cx,k,j+1,i+1);
    }
    if (thrd) {
      y += st_(m,5,k,j,i)*iw_(m,cx,k-1,j,i) + st_(m,6,k,j,i)*iw_(m,cx,k+1,j,i);
      if (edg) {
        y += st_(m,11,k,j,i)*iw_(m,cx,k-1,j,i-1) + st_(m,12,k,j,i)*iw_(m,cx,k-1,j,i+1)
             + st_(m,13,k,j,i)*iw_(m,cx,k+1,j,i-1) + st_(m,14,k,j,i)*iw_(m,cx,k+1,j,i+1)
             + st_(m,15,k,j,i)*iw_(m,cx,k-1,j-1,i) + st_(m,16,k,j,i)*iw_(m,cx,k-1,j+1,i)
             + st_(m,17,k,j,i)*iw_(m,cx,k+1,j-1,i) + st_(m,18,k,j,i)*iw_(m,cx,k+1,j+1,i);
      }
    }
    if (vim) {y += M1VimpRow(iw_, ivb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    if (mus) {y += M1MusclRow(iw_, imb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    if (vfold) {
      y += st_(m,19,k,j,i)*iw_(m,cx,k,j,i-2) + st_(m,20,k,j,i)*iw_(m,cx,k,j,i+2)
           + st_(m,21,k,j,i)*iw_(m,cx,k,j-2,i) + st_(m,22,k,j,i)*iw_(m,cx,k,j+2,i);
      if (thrd) {
        y += st_(m,23,k,j,i)*iw_(m,cx,k-2,j,i) + st_(m,24,k,j,i)*iw_(m,cx,k+2,j,i);
      }
    }
    iw_(m,cy,k,j,i) = y;
    return y;
  };
  if ((rm == 0 || impl_opsplit) && !red_only) {
    par_for("m1_impl_sto", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      row(m, k, j, i);
    });
    if (rm == 0) {return;}
  }
  const bool spl = impl_opsplit || red_only;   // y is already in cy
  // 256-thread blocks: the default 1024-thread block of a reduction with a 4-Real
  // value takes 33 kB of LDS, i.e. ONE block per CU (measured 97 us vs 47 us for the
  // same stencil as a par_for)
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  Real a0 = 0.0, a1 = 0.0, a2 = 0.0, amx = 0.0;
  if (impl_opteam && !spl) {
    // implicit_op_team_red (m1-fast3): one cell per thread, as the plain par_for, and
    // the sums of each team of 256 cells into opt_part; then one small reduction over
    // the teams.  y is bitwise; the sums are summed in another order (round-off).
    constexpr int ts = kM1TeamRed;
    const int nw = (nmb1 + 1)*nkji;
    const int nl = (nw + ts - 1)/ts;
    if (opt_part.extent_int(0) < nl) {
      opt_part = DvceArray2D<Real>("m1_opt_part", nl, 4);
    }
    auto pt_ = opt_part;
    Kokkos::parallel_for("m1_impl_stot", Kokkos::TeamPolicy<DevExeSpace>(nl, ts),
    KOKKOS_LAMBDA(const TeamMember_t &tm) {
      const int idx = tm.league_rank()*ts + tm.team_rank();
      M1R4Val v;
      v.s[0] = 0.0; v.s[1] = 0.0; v.s[2] = 0.0; v.mx = 0.0;
      if (idx < nw) {
        int m = idx/nkji;
        int r = idx - m*nkji;
        int k = r/nji;
        r -= k*nji;
        int j = r/ni;
        int i = r - j*ni;
        k += ks; j += js; i += is;
        const Real y = row(m, k, j, i);
        if (rm == 1 || rm == 4) {
          v.s[0] = iw_(m,M1_IW_KRH,k,j,i)*y;
          if (rm == 4) {v.mx = fabs(iw_(m,M1_IW_KR,k,j,i));}
        } else {
          v.s[0] = y*iw_(m,M1_IW_KS,k,j,i);
          v.s[1] = y*y;
          if (rm == 3) {v.s[2] = iw_(m,M1_IW_KRH,k,j,i)*y;}
        }
      }
      tm.team_reduce(M1R4Red(v));
      if (tm.team_rank() == 0) {
        const int l = tm.league_rank();
        pt_(l,0) = v.s[0];
        pt_(l,1) = v.s[1];
        pt_(l,2) = v.s[2];
        pt_(l,3) = v.mx;
      }
    });
    Kokkos::parallel_reduce("m1_impl_stot2", Kokkos::RangePolicy<DevExeSpace>(0, nl),
    KOKKOS_LAMBDA(const int l, Real &l0, Real &l1, Real &l2, Real &lmx) {
      l0 += pt_(l,0);
      l1 += pt_(l,1);
      l2 += pt_(l,2);
      lmx = (pt_(l,3) > lmx) ? pt_(l,3) : lmx;
    }, a0, a1, a2, Kokkos::Max<Real>(amx));
    out[0] = a0;
    out[1] = a1;
    out[2] = a2;
    out[3] = amx;
    return;
  }
  if (spl) {
    // implicit_op_split_red (m1-fast3): the read-only reduction as its own lambda, so
    // the kernel carries no stencil code.  Same policy, same reducer, same per-cell
    // terms in the same order as the fused kernel below: bitwise the same sums.
    Kokkos::parallel_reduce("m1_impl_stos", pol,
    KOKKOS_LAMBDA(const int idx, Real &l0, Real &l1, Real &l2, Real &lmx) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/ni;
      int i = r - j*ni;
      k += ks; j += js; i += is;
      const Real y = iw_(m,cy,k,j,i);
      if (rm == 1 || rm == 4) {
        l0 += iw_(m,M1_IW_KRH,k,j,i)*y;
        if (rm == 4) {
          Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
          lmx = (a > lmx) ? a : lmx;
        }
      } else {
        l0 += y*iw_(m,M1_IW_KS,k,j,i);
        l1 += y*y;
        if (rm == 3) {l2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
      }
    }, a0, a1, a2, Kokkos::Max<Real>(amx));
    out[0] = a0;
    out[1] = a1;
    out[2] = a2;
    out[3] = amx;
    return;
  }
  Kokkos::parallel_reduce("m1_impl_stor", pol,
  KOKKOS_LAMBDA(const int idx, Real &l0, Real &l1, Real &l2, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real y = row(m, k, j, i);
    if (rm == 1 || rm == 4) {
      l0 += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        lmx = (a > lmx) ? a : lmx;
      }
    } else {
      l0 += y*iw_(m,M1_IW_KS,k,j,i);
      l1 += y*y;
      if (rm == 3) {l2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, a0, a1, a2, Kokkos::Max<Real>(amx));
  out[0] = a0;
  out[1] = a1;
  out[2] = a2;
  out[3] = amx;
}

namespace {
//! implicit_halo_overlap: three sums and one max (|r| >= 0, so 0 is the identity of
//! the max slot), reduced into pinned host memory so the launch does not block the host
struct M1HoVal {
  Real s[3];
  Real mx;
};
struct M1HoRed {
 public:
  using reducer = M1HoRed;
  using value_type = M1HoVal;
  using result_view_type = Kokkos::View<value_type, Kokkos::SharedHostPinnedSpace,
                                        Kokkos::MemoryUnmanaged>;

 private:
  result_view_type value;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit M1HoRed(const result_view_type &v) : value(v) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    for (int q = 0; q < 3; ++q) {d.s[q] += s.s[q];}
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    for (int q = 0; q < 3; ++q) {v.s[q] = 0.0;}
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return false;}
};
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitStencilOpPart
//! \brief implicit_halo_overlap: ImplicitStencilOp on PART of the active cells, with the
//! same per-cell arithmetic.  part = 1: the interior box, every cell at least `w` cells
//! from each face of a non-degenerate direction (it reads no ghost zone, so it may run
//! while the halo exchange is in flight); under implicit_halo_ovl_faces only from the
//! faces whose ghosts arrive by MPI (hm_face); part = 2: the rest (the shell), after
//! the halo.  Under implicit_vimp_fold the row is the folded one of ImplicitStencilOp
//! (before m1-sync the unfolded M1VimpRow was ADDED to the folded stencil, counting
//! the x2/x3 +-1 vimp terms twice: tests_m1/runs_4l_sync).  `red` as
//! ImplicitOffDiagOpC; the partial sums go (asynchronously) into the 4 pinned Reals
//! at `hs` (s0, s1, s2, max), which the caller combines after a fence.

void RadiationM1::ImplicitStencilOpPart(int xc, int yc, int red, int part, int w,
                                        Real *hs) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto st_ = ost;
  const bool thrd = trans_x3;
  const int cx = xc, cy = yc;
  const int rm = red;
  const bool edg = st_edges;
  // implicit_vimp_fold: the vimp part of the row is in the stencil (slots 3-6 and
  // 19-24), exactly as in ImplicitStencilOp; M1VimpRow only when it is not folded
  const bool vfold = impl_vfold && vimp_now;
  const bool vim = vimp_now && !vfold;
  const int ivb = iw_vimp;
  const bool mus = muscl_now;   // xthinfix-1009 Fix B
  const int imb = iw_muscl;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  // the depth of the shell at each face (x1-, x1+, x2-, x2+, x3-, x3+): w at every face
  // of a non-degenerate direction; under implicit_halo_ovl_faces only at the faces
  // whose ghosts arrive by MPI
  int fw[6];
  for (int f = 0; f < 6; ++f) {fw[f] = (!impl_ovl_faces || hm_face[f]) ? w : 0;}
  if (indcs.nx2 == 1) {fw[2] = fw[3] = 0;}
  if (indcs.nx3 == 1) {fw[4] = fw[5] = 0;}
  // the interior box
  const int il = is + fw[0], iu = ie - fw[1];
  const int jl = js + fw[2], ju = je - fw[3];
  const int kl = ks + fw[4], ku = ke - fw[5];
  auto row = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) -> Real {
    Real y = st_(m,0,k,j,i)*iw_(m,cx,k,j,i)
             + st_(m,1,k,j,i)*iw_(m,cx,k,j,i-1) + st_(m,2,k,j,i)*iw_(m,cx,k,j,i+1)
             + st_(m,3,k,j,i)*iw_(m,cx,k,j-1,i) + st_(m,4,k,j,i)*iw_(m,cx,k,j+1,i);
    if (edg) {
      y += st_(m,7,k,j,i)*iw_(m,cx,k,j-1,i-1) + st_(m,8,k,j,i)*iw_(m,cx,k,j-1,i+1)
           + st_(m,9,k,j,i)*iw_(m,cx,k,j+1,i-1) + st_(m,10,k,j,i)*iw_(m,cx,k,j+1,i+1);
    }
    if (thrd) {
      y += st_(m,5,k,j,i)*iw_(m,cx,k-1,j,i) + st_(m,6,k,j,i)*iw_(m,cx,k+1,j,i);
      if (edg) {
        y += st_(m,11,k,j,i)*iw_(m,cx,k-1,j,i-1) + st_(m,12,k,j,i)*iw_(m,cx,k-1,j,i+1)
             + st_(m,13,k,j,i)*iw_(m,cx,k+1,j,i-1) + st_(m,14,k,j,i)*iw_(m,cx,k+1,j,i+1)
             + st_(m,15,k,j,i)*iw_(m,cx,k-1,j-1,i) + st_(m,16,k,j,i)*iw_(m,cx,k-1,j+1,i)
             + st_(m,17,k,j,i)*iw_(m,cx,k+1,j-1,i) + st_(m,18,k,j,i)*iw_(m,cx,k+1,j+1,i);
      }
    }
    if (vim) {y += M1VimpRow(iw_, ivb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    if (mus) {y += M1MusclRow(iw_, imb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    if (vfold) {
      y += st_(m,19,k,j,i)*iw_(m,cx,k,j,i-2) + st_(m,20,k,j,i)*iw_(m,cx,k,j,i+2)
           + st_(m,21,k,j,i)*iw_(m,cx,k,j-2,i) + st_(m,22,k,j,i)*iw_(m,cx,k,j+2,i);
      if (thrd) {
        y += st_(m,23,k,j,i)*iw_(m,cx,k-2,j,i) + st_(m,24,k,j,i)*iw_(m,cx,k+2,j,i);
      }
    }
    iw_(m,cy,k,j,i) = y;
    return y;
  };
  // the cells of the part, flattened per block: part 1 the interior box (k,j,i order),
  // part 2 the shell only -- the end planes at both x3 ends, then per interior k the
  // end rows at both x2 ends and per interior (k,j) the end cells at both x1 ends
  const bool inner = (part == 1);
  const int n1 = ie - is + 1, n2 = je - js + 1;
  const int w1l = fw[0], w1h = fw[1], w2l = fw[2], w2h = fw[3];
  const int w3l = fw[4], w3h = fw[5];
  const int w1 = w1l + w1h, w2 = w2l + w2h, w3 = w3l + w3h;
  const int m1 = iu - il + 1, m2 = ju - jl + 1, m3 = ku - kl + 1;
  const int npl = w3*n2*n1;                 // shell: the x3 end planes
  const int nkc = w2*n1 + m2*w1;            // shell: per interior k
  const int ncell = inner ? (m3*m2*m1) : (npl + m3*nkc);
  if (ncell == 0) {   // no shell (every face in place before the interior operator)
    if (rm != 0) {
      for (int q = 0; q < 4; ++q) {hs[q] = 0.0;}
    }
    return;
  }
  auto cell = KOKKOS_LAMBDA(const int idx, int &m, int &k, int &j, int &i) {
    m = idx/ncell;
    int r = idx - m*ncell;
    if (inner) {
      k = r/(m2*m1);
      r -= k*m2*m1;
      j = r/m1;
      i = il + r - j*m1;
      j += jl;
      k += kl;
      return;
    }
    if (r < npl) {
      const int q = r/(n2*n1);
      r -= q*n2*n1;
      k = (q < w3l) ? (ks + q) : (ku + 1 + (q - w3l));
      j = js + r/n1;
      i = is + r%n1;
      return;
    }
    r -= npl;
    const int kk = r/nkc;
    r -= kk*nkc;
    k = kl + kk;
    if (r < w2*n1) {
      const int q = r/n1;
      j = (q < w2l) ? (js + q) : (ju + 1 + (q - w2l));
      i = is + r%n1;
      return;
    }
    r -= w2*n1;
    j = jl + r/w1;
    const int q = r%w1;
    i = (q < w1l) ? (is + q) : (iu + 1 + (q - w1l));
  };
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*ncell);
  if (rm == 0) {
    Kokkos::parallel_for("m1_impl_stop", pol, KOKKOS_LAMBDA(const int idx) {
      int m, k, j, i;
      cell(idx, m, k, j, i);
      row(m, k, j, i);
    });
    return;
  }
  M1HoRed::result_view_type res(reinterpret_cast<M1HoVal *>(hs));
  if (impl_opteam) {
    // implicit_op_team_red (m1-fast3): as in ImplicitStencilOp, one cell per thread in
    // teams, per-team partials in opt_part, then a small reduction into the pinned hs
    // (still asynchronous).  Same stream, so part 2 cannot overwrite opt_part before
    // part 1's small reduction has read it.  Round-off vs the fused reduction.
    constexpr int ts = kM1TeamRed;
    const int nw = (nmb1 + 1)*ncell;
    const int nl = (nw + ts - 1)/ts;
    if (opt_part.extent_int(0) < nl) {
      // sized for the whole pack at once, so parts 1 and 2 never reallocate
      const int nall = (nmb1 + 1)*(ke - ks + 1)*(je - js + 1)*(ie - is + 1);
      opt_part = DvceArray2D<Real>("m1_opt_part", std::max(nl, (nall + ts - 1)/ts), 4);
    }
    auto pt_ = opt_part;
    Kokkos::parallel_for("m1_impl_stopt", Kokkos::TeamPolicy<DevExeSpace>(nl, ts),
    KOKKOS_LAMBDA(const TeamMember_t &tm) {
      const int idx = tm.league_rank()*ts + tm.team_rank();
      M1R4Val v;
      v.s[0] = 0.0; v.s[1] = 0.0; v.s[2] = 0.0; v.mx = 0.0;
      if (idx < nw) {
        int m, k, j, i;
        cell(idx, m, k, j, i);
        const Real y = row(m, k, j, i);
        if (rm == 1 || rm == 4) {
          v.s[0] = iw_(m,M1_IW_KRH,k,j,i)*y;
          if (rm == 4) {v.mx = fabs(iw_(m,M1_IW_KR,k,j,i));}
        } else {
          v.s[0] = y*iw_(m,M1_IW_KS,k,j,i);
          v.s[1] = y*y;
          if (rm == 3) {v.s[2] = iw_(m,M1_IW_KRH,k,j,i)*y;}
        }
      }
      tm.team_reduce(M1R4Red(v));
      if (tm.team_rank() == 0) {
        const int l = tm.league_rank();
        pt_(l,0) = v.s[0];
        pt_(l,1) = v.s[1];
        pt_(l,2) = v.s[2];
        pt_(l,3) = v.mx;
      }
    });
    Kokkos::parallel_reduce("m1_impl_stopt2", Kokkos::RangePolicy<DevExeSpace>(0, nl),
    KOKKOS_LAMBDA(const int l, M1HoVal &v) {
      v.s[0] += pt_(l,0);
      v.s[1] += pt_(l,1);
      v.s[2] += pt_(l,2);
      v.mx = (pt_(l,3) > v.mx) ? pt_(l,3) : v.mx;
    }, M1HoRed(res));
    return;
  }
  Kokkos::parallel_reduce("m1_impl_stopr", pol,
  KOKKOS_LAMBDA(const int idx, M1HoVal &v) {
    int m, k, j, i;
    cell(idx, m, k, j, i);
    Real y = row(m, k, j, i);
    if (rm == 1 || rm == 4) {
      v.s[0] += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        v.mx = (a > v.mx) ? a : v.mx;
      }
    } else {
      v.s[0] += y*iw_(m,M1_IW_KS,k,j,i);
      v.s[1] += y*y;
      if (rm == 3) {v.s[2] += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, M1HoRed(res));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOpX
//! \brief y = A x (+ the reductions `red` of ImplicitOffDiagOpC) on the fast path: from
//! the 19-point stencil under implicit_op_stencil, else from the od cache.  The caller
//! has filled the ghost zones of x.

void RadiationM1::ImplicitOpX(int xc, int yc, int red, Real *out) {
  if (impl_stencil) {
    ImplicitStencilOp(xc, yc, red, out);
  } else {
    ImplicitOffDiagOpC(xc, yc, 1.0, true, red, out);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOffDiagOp
//! \brief MILESTONE 3b phase D: accumulate sgn * L_off(x) into the component yc, where
//! L_off is the contribution the OFF-DIAGONAL Eddington terms of the face equations make
//! to the cell row:
//!
//!   L_off(x)_c = - sum_d (dt/dx_d) (chat/c) chat c dt
//!                  [ theta_{d,+} od_{d,+}(x) - theta_{d,-} od_{d,-}(x) ],
//!   od_{d,f}(x) = (1/2) [ (sum_{e!=d} d_e D_de x)_L + (sum_{e!=d} d_e D_de x)_R ],
//!
//! i.e. exactly the term the face fluxes of ImplicitTransverseTerms and of the x1
//! assembly carry, with the LAGGED energy replaced by the argument x.  The closure
//! (chi, n, hence D_ab) is frozen inside a Picard pass, so this IS a linear operator: a
//! 9-point stencil in 2-D and a 19-point one in 3-D, built from the same one-layer halo
//! the 7-point operator already exchanges (a cell-centred exchange fills the edge and
//! corner neighbours, which is what the cross derivatives read).
//!
//! It is used twice per Picard pass under implicit_offdiag = operator: once with x = E^k
//! to REMOVE the lagged term from the right-hand side the assembly produced, and inside
//! every operator application to put it back on the left.  Both use the same routine, so
//! the two cannot drift apart.
//!
//! The boundary logic mirrors the assembly face by face: a physical x1 face carries an
//! imposed flux and no off-diagonal term, a physical x2/x3 face is reflecting (F = 0),
//! and an M1_IBC_EFIX Dirichlet row is replaced whole and gets nothing at all.

void RadiationM1::ImplicitOffDiagOp(int xc, int yc, Real sgn) {
  if (impl_odc) {   // the same numbers from the per-cell cache
    ImplicitOffDiagOpC(xc, yc, sgn, false, 0, nullptr);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  // capture the DEVICE Views only: whole DualViews pushed this functor past 512 bytes,
  // i.e. onto Kokkos HIP's constant-memory launch, which waits on the previous such
  // launch (hip_event_synchronize) -- one hidden host stall per call, ~290 per cycle.
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = (cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs).d_view;   // CS1: seams open
  auto pos_ = part_pos.d_view;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  const int nblkx1 = part_nblk;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  // the TRANSVERSE face theta of the operator must be the very number the face fluxes
  // and the cell terms used, limiter or not (ImplicitTransTheta).  The x1 faces are not
  // limited.
  const bool lm = (impl_tlim != M1_TLIM_NONE) || blat_on;
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  Real cl = c_light, ch = chat, dt = dt_sub;
  const int cx = xc, cy = yc;
  const Real sg = sgn;
  const bool bx1 = bvec_x1 && trans_on;   // blendall-1009: x1 faces keep 1 - AL
  auto fw_ = ifw;
  par_for("m1_impl_odop", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                    (i == ie && topb && bchi == M1_IBC_EFIX))) {
      return;
    }
    Real dx1 = mbsize(m).dx1;
    Real dx2 = mbsize(m).dx2;
    Real dx3 = mbsize(m).dx3;
    BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
    BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
    if (!p2lo) {jl = js-1;}
    if (!p2hi) {ju = je+1;}
    if (thrd && !p3lo) {kl = ks-1;}
    if (thrd && !p3hi) {ku = ke+1;}
    Real cr = ch/cl;
    Real kk = ch*cl*dt;
    Real y = 0.0;
    // ---- the two x1 faces
    if (i < ie || cyclic || !topb) {
      int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (bx1) {th *= 1.0 - fw_(m,M1_IFW_AL,k,j,i+1);}
      Real od = 0.5*(M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,0,k,j,ip,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y -= (dt/dx1)*cr*th*kk*od;
    }
    if (i > is || cyclic || !botb) {
      int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (bx1) {th *= 1.0 - fw_(m,M1_IFW_AL,k,j,i);}
      Real od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y += (dt/dx1)*cr*th*kk*od;
    }
    // ---- the two x2 faces
    Real nu2 = dt/dx2;
    if (!(j == je && p2hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
      Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,1,k,j+1,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y -= nu2*cr*th*kk*od;
    }
    if (!(j == js && p2lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,1,k,j-1,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y += nu2*cr*th*kk*od;
    }
    // ---- the two x3 faces
    if (thrd) {
      Real nu3 = dt/dx3;
      if (!(k == ke && p3hi)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
        Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real od = 0.5*(M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull)
                       + M1OffDiv(iw_,m,2,k+1,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,
                                  vd_,dfull));
        y -= nu3*cr*th*kk*od;
      }
      if (!(k == ks && p3lo)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
        Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real od = 0.5*(M1OffDiv(iw_,m,2,k-1,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull)
                       + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                  dfull));
        y += nu3*cr*th*kk*od;
      }
    }
    iw_(m,cy,k,j,i) += sg*y;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitApplyOp
//! \brief milestone 3b phase C: y = A x for the FROZEN 7-point operator of the current
//! Picard pass -- the assembled tridiagonal row (M1_IW_TA, TB, TC, with TB already
//! carrying the transverse diagonal M1_IW_TDIA) plus the four/six transverse
//! off-diagonals (M1_IW_CJM..CKP).  One halo exchange of x, then one stencil kernel.
//!
//! The x1 wrap of a periodic single-block column is handled exactly as the assembly
//! handles it; with several blocks stacked along x1 the coupling to the neighbouring
//! block is through the ghost cell the halo just filled, which is the same cell the
//! gathered line solve treats as an interior row.

void RadiationM1::ImplicitApplyOp(int xc, int yc) {
  if (impl_stencil) {   // the 19-point stencil of this pass
    ImplicitHaloOp(xc, yc, 0, nullptr);
    return;
  }
  ImplicitKrylovHalo(xc);
  if (impl_odc) {   // 7-point row + off-diagonal Eddington terms in one kernel
    ImplicitOffDiagOpC(xc, yc, 1.0, true, 0, nullptr);
    if (VlatOp()) {VetLatOp(xc, yc, 1.0);}
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  const int cx = xc, cy = yc;
  const bool vim = vimp_now;
  const int ivb = iw_vimp;
  const bool mus = muscl_now;   // xthinfix-1009 Fix B
  const int imb = iw_muscl;
  par_for("m1_impl_op", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
    int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
    Real y = iw_(m,M1_IW_TB,k,j,i)*iw_(m,cx,k,j,i)
             + iw_(m,M1_IW_TA,k,j,i)*iw_(m,cx,k,j,im)
             + iw_(m,M1_IW_TC,k,j,i)*iw_(m,cx,k,j,ip)
             + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,cx,k,j-1,i)
             + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,cx,k,j+1,i);
    if (thrd) {
      y += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,cx,k-1,j,i)
           + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,cx,k+1,j,i);
    }
    if (vim) {y += M1VimpRow(iw_, ivb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    if (mus) {y += M1MusclRow(iw_, imb, cx, m, k, j, i, is, ie, cyclic, thrd);}
    iw_(m,cy,k,j,i) = y;
  });
  // MILESTONE 3b phase D: the off-diagonal Eddington coupling, when it is part of the
  // operator.  The halo of x above is a full cell-centred exchange, so the edge and
  // corner neighbours the cross derivatives read are already filled: the 9-/19-point
  // stencil costs one extra kernel and no communication at all.
  if (od_now == M1_OD_OPERATOR) {
    ImplicitOffDiagOp(xc, yc, 1.0);
  }
  if (VlatOp()) {VetLatOp(xc, yc, 1.0);}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPrecond
//! \brief milestone 3b phase C: z = M^{-1} r, with M the x1 line part of the row (the
//! exact tridiagonal solve including the full diagonal).  It stages r through M1_IW_TR,
//! which the assembly has already been read out of by the time the Krylov loop runs.

void RadiationM1::ImplicitPrecond(int rc, int zc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const int cr = rc, cz = zc;
  if (cr >= 0) {
    par_for("m1_impl_prein", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_TR,k,j,i) = iw_(m,cr,k,j,i);
    });
  }
  ImplicitTridiagSolve();
  par_for("m1_impl_preout", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,cz,k,j,i) = iw_(m,M1_IW_S2,k,j,i);
  });
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1GlobalSum2
//! \brief two global SUMs in one MPI_Allreduce.  The per-rank part is a Kokkos reduction
//! over the same fixed index range every time, so its summation order is reproducible on
//! a given rank count, and the cross-rank part is one collective.

void M1GlobalSum2(Real &a, Real &b, const bool det = false) {
#if MPI_PARALLEL_ENABLED
  if (det && global_variable::nranks > 1) {
    // implicit_det_reduce: gather the per-rank values, sum them in rank order
    std::vector<Real> all(2*global_variable::nranks);
    Real loc[2] = {a, b};
    MPI_Allgather(loc, 2, MPI_ATHENA_REAL, all.data(), 2, MPI_ATHENA_REAL,
                  MPI_COMM_WORLD);
    a = 0.0;
    b = 0.0;
    for (int r = 0; r < global_variable::nranks; ++r) {
      a += all[2*r];
      b += all[2*r + 1];
    }
    return;
  }
  Real loc[2] = {a, b}, glb[2];
  MPI_Allreduce(loc, glb, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  a = glb[0];
  b = glb[1];
#else
  (void)a;
  (void)b;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1GlobalSumArr
//! \brief n global SUMs in one MPI_Allreduce (milestone 3e: the row of the Anderson
//! normal equations).  Same reproducibility argument as M1GlobalSum2.

void M1GlobalSumArr(Real *a, const int n) {
#if MPI_PARALLEL_ENABLED
  Real glb[NREDUCTION_VARIABLES];
  MPI_Allreduce(a, glb, n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < n; ++q) {a[q] = glb[q];}
#else
  (void)a;
  (void)n;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1SolveSPD
//! \brief solve the small SYMMETRIC system A g = b in place by Gaussian elimination with
//! partial pivoting (n <= M1_AND_MMAX).  Returns false if the matrix is numerically
//! singular, in which case the caller falls back to a plain Picard pass.

bool M1SolveSPD(Real *a, Real *b, const int n) {
  for (int p = 0; p < n; ++p) {
    int piv = p;
    Real amx = fabs(a[p*n + p]);
    for (int r = p+1; r < n; ++r) {
      Real v = fabs(a[r*n + p]);
      if (v > amx) {amx = v; piv = r;}
    }
    if (!(amx > 0.0)) {return false;}
    if (piv != p) {
      for (int q = 0; q < n; ++q) {
        Real t = a[p*n + q];
        a[p*n + q] = a[piv*n + q];
        a[piv*n + q] = t;
      }
      Real t = b[p];
      b[p] = b[piv];
      b[piv] = t;
    }
    Real ip = 1.0/a[p*n + p];
    for (int r = p+1; r < n; ++r) {
      Real fct = a[r*n + p]*ip;
      if (fct == 0.0) {continue;}
      for (int q = p; q < n; ++q) {a[r*n + q] -= fct*a[p*n + q];}
      b[r] -= fct*b[p];
    }
  }
  for (int r = n-1; r >= 0; --r) {
    Real s = b[r];
    for (int q = r+1; q < n; ++q) {s -= a[r*n + q]*b[q];}
    b[r] = s/a[r*n + r];
  }
  for (int r = 0; r < n; ++r) {
    if (!std::isfinite(b[r])) {return false;}
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn M1GlobalMax
void M1GlobalMax(Real &a) {
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&a, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  a = g;
#else
  (void)a;
#endif
}

//----------------------------------------------------------------------------------------
//! \struct M1BcgVal, M1BcgRed
//! \brief implicit_bcg_sync > 0: ONE reduction kernel returns up to three sums and one
//! max (e.g. rhat.r and max|r| straight from the kernel that updates r).  M1BcgRed is a
//! Kokkos reducer of the same shape as Kokkos::Sum; |r| >= 0, so 0 is the identity of
//! the max slot.

struct M1BcgVal {
  Real s0, s1, s2, mx;
};

template <class Space>
struct M1BcgRed {
 public:
  using reducer = M1BcgRed<Space>;
  using value_type = M1BcgVal;
  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool refs;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit M1BcgRed(value_type &v) : value(&v), refs(true) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    d.s0 += s.s0;
    d.s1 += s.s1;
    d.s2 += s.s2;
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    v.s0 = 0.0;
    v.s1 = 0.0;
    v.s2 = 0.0;
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return refs;}
};

//! the flattened (m,k,j,i) index of the fused reductions, as par_for flattens it
KOKKOS_INLINE_FUNCTION
void M1BcgIdx(const int idx, const int nkji, const int nji, const int ni,
              int &m, int &k, int &j, int &i) {
  m = idx/nkji;
  int r = idx - m*nkji;
  k = r/nji;
  r -= k*nji;
  j = r/ni;
  i = r - j*ni;
}

#if MPI_PARALLEL_ENABLED
//! the MPI operation of M1GlobalBcg: sum the first three slots, max the fourth
void M1BcgOpFn(void *in, void *inout, int *len, MPI_Datatype *) {
  Real *a = static_cast<Real *>(in);
  Real *b = static_cast<Real *>(inout);
  for (int q = 0; q < *len; ++q) {
    b[4*q] += a[4*q];
    b[4*q + 1] += a[4*q + 1];
    b[4*q + 2] += a[4*q + 2];
    b[4*q + 3] = (a[4*q + 3] > b[4*q + 3]) ? a[4*q + 3] : b[4*q + 3];
  }
}
#endif

//----------------------------------------------------------------------------------------
//! \fn M1GlobalBcg
//! \brief the three sums and the max of an M1BcgVal over all ranks in ONE MPI_Allreduce
//! (a 4-Real contiguous type with a user operation).  Nothing to do on one rank.

void M1GlobalBcg(M1BcgVal &v, const bool det = false) {
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks == 1) {return;}
  if (det) {
    // implicit_det_reduce: gather the per-rank values, combine them in rank order
    std::vector<Real> all(4*global_variable::nranks);
    Real loc[4] = {v.s0, v.s1, v.s2, v.mx};
    MPI_Allgather(loc, 4, MPI_ATHENA_REAL, all.data(), 4, MPI_ATHENA_REAL,
                  MPI_COMM_WORLD);
    v.s0 = v.s1 = v.s2 = v.mx = 0.0;
    for (int r = 0; r < global_variable::nranks; ++r) {
      v.s0 += all[4*r];
      v.s1 += all[4*r + 1];
      v.s2 += all[4*r + 2];
      v.mx = (all[4*r + 3] > v.mx) ? all[4*r + 3] : v.mx;
    }
    return;
  }
  static MPI_Datatype typ = MPI_DATATYPE_NULL;
  static MPI_Op op = MPI_OP_NULL;
  if (op == MPI_OP_NULL) {
    MPI_Type_contiguous(4, MPI_ATHENA_REAL, &typ);
    MPI_Type_commit(&typ);
    MPI_Op_create(&M1BcgOpFn, 1, &op);
  }
  Real loc[4] = {v.s0, v.s1, v.s2, v.mx}, glb[4];
  MPI_Allreduce(loc, glb, 1, typ, op, MPI_COMM_WORLD);
  v.s0 = glb[0];
  v.s1 = glb[1];
  v.s2 = glb[2];
  v.mx = glb[3];
#else
  (void)v;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1DetReduce
//! \brief <rad_m1>/implicit_det_reduce: the M1BcgVal reduction of fn(idx, v) over
//! idx in [0, n) in a FIXED order that does not depend on the backend, the launch
//! configuration or the thread scheduling.  Level 1: M1_DR_NL teams; slot t of team l
//! sums idx = l*M1_DR_NT + t + q*M1_DR_NL*M1_DR_NT in increasing q (coalesced), then the
//! M1_DR_NT slots are combined by a fixed pairwise tree in team scratch.  Level 2: one
//! team does the same over the M1_DR_NL team partials.  The 4 results are copied to the
//! host.  The tree is written out explicitly (TeamThreadRange + barriers), so the result
//! is the same for any team size Kokkos picks (AUTO): run-to-run bitwise on the GPU.

constexpr int M1_DR_NT = 256;    // slots per team (a power of 2)
constexpr int M1_DR_NL = 1024;   // level-1 teams (a multiple of M1_DR_NT)

template <class F>
void M1DetReduce(const char *name, const int n, const F &fn,
                 const DvceArray1D<Real> &part, M1BcgVal &out) {
  using TP = Kokkos::TeamPolicy<DevExeSpace>;
  using Mem = Kokkos::View<Real*, ScratchMemSpace,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  const size_t sb = Mem::shmem_size(4*M1_DR_NT);
  // part: [0, 4*NL) the level-1 partials, [4*NL, 4*NL + 4) the result
  Kokkos::parallel_for(name, TP(DevExeSpace(), M1_DR_NL, Kokkos::AUTO)
                       .set_scratch_size(0, Kokkos::PerTeam(sb)),
  KOKKOS_LAMBDA(const TP::member_type &tm) {
    Mem sh(tm.team_scratch(0), 4*M1_DR_NT);
    const int l = tm.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, M1_DR_NT), [&](const int t) {
      M1BcgVal v;
      v.s0 = v.s1 = v.s2 = v.mx = 0.0;
      for (int q = l*M1_DR_NT + t; q < n; q += M1_DR_NL*M1_DR_NT) {fn(q, v);}
      sh(4*t) = v.s0;
      sh(4*t + 1) = v.s1;
      sh(4*t + 2) = v.s2;
      sh(4*t + 3) = v.mx;
    });
    tm.team_barrier();
    for (int w = M1_DR_NT/2; w > 0; w /= 2) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, w), [&](const int t) {
        sh(4*t) += sh(4*(t + w));
        sh(4*t + 1) += sh(4*(t + w) + 1);
        sh(4*t + 2) += sh(4*(t + w) + 2);
        sh(4*t + 3) = (sh(4*(t + w) + 3) > sh(4*t + 3)) ? sh(4*(t + w) + 3) : sh(4*t + 3);
      });
      tm.team_barrier();
    }
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      for (int c = 0; c < 4; ++c) {part(4*l + c) = sh(c);}
    });
  });
  Kokkos::parallel_for("m1_det_red2", TP(DevExeSpace(), 1, Kokkos::AUTO)
                       .set_scratch_size(0, Kokkos::PerTeam(sb)),
  KOKKOS_LAMBDA(const TP::member_type &tm) {
    Mem sh(tm.team_scratch(0), 4*M1_DR_NT);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, M1_DR_NT), [&](const int t) {
      Real a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
      for (int q = t; q < M1_DR_NL; q += M1_DR_NT) {
        a0 += part(4*q);
        a1 += part(4*q + 1);
        a2 += part(4*q + 2);
        a3 = (part(4*q + 3) > a3) ? part(4*q + 3) : a3;
      }
      sh(4*t) = a0;
      sh(4*t + 1) = a1;
      sh(4*t + 2) = a2;
      sh(4*t + 3) = a3;
    });
    tm.team_barrier();
    for (int w = M1_DR_NT/2; w > 0; w /= 2) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, w), [&](const int t) {
        sh(4*t) += sh(4*(t + w));
        sh(4*t + 1) += sh(4*(t + w) + 1);
        sh(4*t + 2) += sh(4*(t + w) + 2);
        sh(4*t + 3) = (sh(4*(t + w) + 3) > sh(4*t + 3)) ? sh(4*(t + w) + 3) : sh(4*t + 3);
      });
      tm.team_barrier();
    }
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      for (int c = 0; c < 4; ++c) {part(4*M1_DR_NL + c) = sh(c);}
    });
  });
  Real h[4];
  Kokkos::View<Real*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      hv(h, 4);
  Kokkos::deep_copy(hv, Kokkos::subview(part, Kokkos::make_pair(4*M1_DR_NL,
                                                                4*M1_DR_NL + 4)));
  out.s0 = h[0];
  out.s1 = h[1];
  out.s2 = h[2];
  out.mx = h[3];
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitAccelSave
//! \brief MILESTONE 3e: snapshot the SCALED state x_k entering a Picard pass.
//!
//! The fixed-point vector is the COMPLETE lagged state one pass maps to the next: the
//! solve unknown E and the three cell-centred fluxes F_1, F_2, F_3, which are what the
//! top of the next pass rebuilds the closure (chi, n, hence the whole Eddington tensor)
//! from.  Nothing else a pass reads is independent state: the face fluxes f0x1/f0x2/f0x3
//! are recomputed from E at the end of every pass, the velocities and opacities are
//! frozen over the step, and the temperature is slaved to E by the scalar root find.
//!
//! SCALING.  E and F are not commensurate (F ~ c E), and an unscaled 2-norm would be
//! dominated by the optically thick base, where E is ~ 10 orders above the thin top --
//! which is exactly where the iteration does NOT need help.  Each cell is therefore
//! divided by its OWN energy scale S = max(E^n, e_floor), fixed over the whole step (so
//! the least-squares problem stays linear across passes), and the flux components by c S:
//!   x = ( E/S , F_1/(c S) , F_2/(c S) , F_3/(c S) ) ,
//! i.e. the Anderson residual is the RELATIVE fixed-point residual, cell by cell, and the
//! thin top carries the same weight as the base.

void RadiationM1::ImplicitAccelSave() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto xc_ = aa_xc;
  auto sc_ = aa_sc;
  const Real cl = c_light;
  const int nc = aa_nc;
  par_for("m1_acc_save", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
    Real w = 1.0/(sc_(m,k,j,i)*((c == 0) ? 1.0 : cl));
    xc_(m,c,k,j,i) = iw_(m,M1AccComp(c),k,j,i)*w;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitAccelApply
//! \brief MILESTONE 3e: ANDERSON ACCELERATION (Walker & Ni 2011, SIAM J. Numer.
//! Anal. 49, 1715) of the Picard map, applied at the END of pass `it`, where the iw
//! state is G(x_k).
//!
//! With g_k = G(x_k) - x_k and the histories dX_t = x_{t+1} - x_t and dG_t = g_{t+1} -
//! g_t of the last m passes, the accelerated iterate is
//!
//!   gamma = argmin_gamma || g_k - dG gamma ||_2 ,
//!   x_{k+1} = x_k + beta g_k - (dX + beta dG) gamma ,
//!
//! i.e. the beta-mixed map evaluated at the point of the affine span of the last m
//! iterates whose linear residual model is smallest.  beta = 1 is the plain (undamped)
//! form.  The least-squares problem is solved on the HOST from the NORMAL equations
//! (dG^T dG + lambda I) gamma = dG^T g_k, with lambda = M1_AND_REG tr(dG^T dG)/m a
//! Tikhonov term that makes a nearly dependent history harmless; every inner product is a
//! GLOBAL sum (one MPI_Allreduce of m+1 doubles per history row), so the answer does not
//! depend on the decomposition.
//!
//! REALIZABILITY is re-applied to the accelerated iterate (E >= e_floor, |F| <= c E)
//! before it is written back: the extrapolation is a linear combination of realizable
//! states and the M1 admissible set is convex in (E,F), but the E floor and the per-cell
//! scaling break exact convexity, and an inadmissible lagged state would give the next
//! pass a closure with chi outside [1/3,1].
//!
//! SAFEGUARD: if ||g_k|| exceeds 10x the previous pass', the history is dropped and the
//! pass falls back to plain Picard (x_{k+1} = G(x_k)); the event is counted.  The
//! convergence test of the outer loop is untouched -- it still measures |dE|/E and the
//! true linear residual of the UNACCELERATED map, which is the fixed-point residual.

void RadiationM1::ImplicitAccelApply(int it) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto xc_ = aa_xc;
  auto fc_ = aa_fc;
  auto xp_ = aa_xp;
  auto fp_ = aa_fp;
  auto dx_ = aa_dx;
  auto df_ = aa_df;
  auto sc_ = aa_sc;
  const Real cl = c_light;
  const Real efl = e_floor;
  const int nc = aa_nc;
  const int mm = impl_and_m;

  // (1) the fixed-point residual of this pass, in the scaled variables
  par_for("m1_acc_res", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
    Real w = 1.0/(sc_(m,k,j,i)*((c == 0) ? 1.0 : cl));
    fc_(m,c,k,j,i) = iw_(m,M1AccComp(c),k,j,i)*w - xc_(m,c,k,j,i);
  });

  // (2) its GLOBAL 2-norm
  Real fn2 = 0.0;
  Kokkos::parallel_reduce("m1_acc_nrm",
  Kokkos::MDRangePolicy<Kokkos::Rank<5>>(DevExeSpace(), {0,0,ks,js,is},
                                         {nmb1+1,nc,ke+1,je+1,ie+1}),
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i,
                Real &lsum) {
    lsum += SQR(fc_(m,c,k,j,i));
  }, Kokkos::Sum<Real>(fn2));
  {Real dum = 0.0;
  M1GlobalSum2(fn2, dum);}
  Real fnorm = sqrt(fmax(fn2, 0.0));

  // (3) the divergence safeguard
  bool restart = false;
  if (aa_fnp > 0.0 && fnorm > 10.0*aa_fnp) {
    aa_nh = 0;
    aa_head = 0;
    aa_hasp = false;
    aa_nrst += 1.0;
    restart = true;
  }

  // (4) push (dX, dG) onto the ring, then remember this pass
  if (aa_hasp) {
    const int q = aa_head;
    par_for("m1_acc_push", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
      dx_(m,q,c,k,j,i) = xc_(m,c,k,j,i) - xp_(m,c,k,j,i);
      df_(m,q,c,k,j,i) = fc_(m,c,k,j,i) - fp_(m,c,k,j,i);
    });
    aa_head = (aa_head + 1) % mm;
    aa_nh = std::min(aa_nh + 1, mm);
  }
  Kokkos::deep_copy(DevExeSpace(), aa_xp, aa_xc);
  Kokkos::deep_copy(DevExeSpace(), aa_fp, aa_fc);
  aa_hasp = true;
  aa_fnp = fnorm;

  // nothing to extrapolate from (or the safeguard fired, or the pass is before
  // implicit_anderson_start): leave the iw state as G(x_k), which IS the Picard iterate
  if (restart || it < impl_and_start || aa_nh < 1) {return;}

  // (5) the m x m normal equations, one GLOBAL reduction of m+1 numbers per row.  The
  // history columns in use are the last aa_nh pushed onto the ring.
  const int nh = aa_nh;
  int cid[M1_AND_MMAX];
  for (int t = 0; t < nh; ++t) {cid[t] = (aa_head - nh + t + mm) % mm;}
  Real amat[M1_AND_MMAX*M1_AND_MMAX], bvec[M1_AND_MMAX];
  Real trc = 0.0;
  for (int t = 0; t < nh; ++t) {
    const int qr = cid[t];
    array_sum::GlobalSum row;
    Kokkos::parallel_reduce("m1_acc_dot",
    Kokkos::MDRangePolicy<Kokkos::Rank<5>>(DevExeSpace(), {0,0,ks,js,is},
                                           {nmb1+1,nc,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i,
                  array_sum::GlobalSum &lsum) {
      Real dq = df_(m,qr,c,k,j,i);
      for (int r = 0; r < nh; ++r) {
        lsum.the_array[r] += dq*df_(m,cid[r],c,k,j,i);
      }
      lsum.the_array[nh] += dq*fc_(m,c,k,j,i);
    }, Kokkos::Sum<array_sum::GlobalSum>(row));
    M1GlobalSumArr(row.the_array, nh+1);
    for (int r = 0; r < nh; ++r) {amat[t*nh + r] = row.the_array[r];}
    bvec[t] = row.the_array[nh];
    trc += row.the_array[t];
  }
  const Real lam = M1_AND_REG*fmax(trc, 1.0e-300)/static_cast<Real>(nh);
  for (int t = 0; t < nh; ++t) {amat[t*nh + t] += lam;}
  if (!M1SolveSPD(amat, bvec, nh)) {return;}

  // (6) the accelerated iterate, with the realizability limits re-applied
  Real gam[M1_AND_MMAX];
  for (int t = 0; t < nh; ++t) {gam[t] = bvec[t];}
  const Real bta = impl_and_beta;
  par_for("m1_acc_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real xn[4] = {0.0, 0.0, 0.0, 0.0};
    for (int c = 0; c < nc; ++c) {
      Real v = xc_(m,c,k,j,i) + bta*fc_(m,c,k,j,i);
      for (int t = 0; t < nh; ++t) {
        const int q = cid[t];
        v -= gam[t]*(dx_(m,q,c,k,j,i) + bta*df_(m,q,c,k,j,i));
      }
      xn[c] = v;
    }
    Real s = sc_(m,k,j,i);
    Real e = xn[0]*s;
    if (!(e > efl)) {e = efl;}
    Real f1 = xn[1]*cl*s;
    Real f2 = (nc > 2) ? (xn[2]*cl*s) : 0.0;
    Real f3 = (nc > 3) ? (xn[3]*cl*s) : 0.0;
    Real fm = sqrt(f1*f1 + f2*f2 + f3*f3);
    Real fmx = cl*e;
    if (fm > fmx) {
      Real sf = fmx/fm;
      f1 *= sf;
      f2 *= sf;
      f3 *= sf;
    }
    iw_(m,M1_IW_EP,k,j,i) = e;
    iw_(m,M1_IW_F1,k,j,i) = f1;
    if (nc > 2) {
      iw_(m,M1_IW_F2,k,j,i) = f2;
      iw_(m,M1_IW_F3,k,j,i) = f3;
    }
  });
  aa_nacc += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStab
//! \brief milestone 3b phase C: solve the FROZEN 7-point linear system of one Picard pass
//! by matrix-free BiCGStab, RIGHT-preconditioned by the exact x1 line solve, and leave
//! the answer in M1_IW_S2 -- exactly where the line-Jacobi pass leaves it, so nothing
//! downstream of step (e) of ImplicitSolve knows which solver ran.
//!
//! The system is A x = b with
//!   A x = tridiag(TA,TB,TC) x + CJM x_{j-1} + CJP x_{j+1} + CKM x_{k-1} + CKP x_{k+1},
//!   b   = TR + sum_nb C_nb E^k_nb,
//! i.e. the right-hand side line Jacobi uses PLUS the lagged off-diagonal term it had
//! moved there.  One line-Jacobi pass is x <- M^{-1}(b - sum C x) with the same M and the
//! same A, so the two solvers have the SAME fixed point and converge to the same answer;
//! only the number of passes differs.  The initial guess is the Picard iterate itself.
//!
//! Stopping is on the TRUE residual max|b - A x| relative to max|b|, below
//! implicit_lin_tol: the recursive residual is monitored every iteration (it is free)
//! and, when it passes, ONE extra operator application checks the true one.  If the true
//! residual disagrees the recurrence is RESTARTED from it, which is the standard cure for
//! the drift between the two.
//!
//! BREAKDOWN (rho or rhat.v underflowing, omega vanishing) restarts the recurrence once
//! from the current iterate with a fresh shadow residual; a second breakdown in the same
//! pass falls back to ONE line-Jacobi update, which always exists, and is counted.
//!
//! Three global reductions per iteration (rhat.r with r.r, rhat.v, and t.s with t.t),
//! five scalars in all.

int RadiationM1::ImplicitBiCGStab(Real rhsmax) {
  if (impl_opchk != 0) {ImplicitOpCheck();}   // implicit_op_check (debug, default off)
  if (impl_bcg_sync > 0) {return ImplicitBiCGStabFused(rhsmax);}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  Kokkos::MDRangePolicy<Kokkos::Rank<4>> rng(DevExeSpace(), {0,ks,js,is},
                                             {nmb1+1,ke+1,je+1,ie+1});

  // x0 = the Picard iterate; r0 = b - A x0
  par_for("m1_impl_bcg_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  par_for("m1_impl_bcg_r0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
  });
  Real rnorm = 0.0;
  Kokkos::parallel_reduce("m1_impl_bcg_rn", rng,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
    Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
    lmax = (r > lmax) ? r : lmax;
  }, Kokkos::Max<Real>(rnorm));
  M1GlobalMax(rnorm);
  bcg_nred += 1.0;
  bcg_r0rel = rnorm/bscale;

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = (rnorm/bscale < tol);
  bool fell_back = false;
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    // rho_new = (rhat, r)
    Real rhon = 0.0;
    Kokkos::parallel_reduce("m1_impl_bcg_rho", rng,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
      ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KR,k,j,i);
    }, rhon);
    Real dummy = 0.0;
    M1GlobalSum2(rhon, dummy);
    bcg_nred += 1.0;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      const Real bt = beta, om = omega;
      par_for("m1_impl_bcg_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        iw_(m,M1_IW_KP,k,j,i) = iw_(m,M1_IW_KR,k,j,i)
            + bt*(iw_(m,M1_IW_KP,k,j,i) - om*iw_(m,M1_IW_KV,k,j,i));
      });
      ImplicitPrecond(M1_IW_KP, M1_IW_KY);
      ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
      Real rv = 0.0;
      Kokkos::parallel_reduce("m1_impl_bcg_rv", rng,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
        ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KV,k,j,i);
      }, rv);
      Real d2 = 0.0;
      M1GlobalSum2(rv, d2);
      bcg_nred += 1.0;
      if (!(fabs(rv) > M1_BCG_EPS)) {
        breakdown = true;
      } else {
        alpha = rhon/rv;
        const Real al = alpha;
        par_for("m1_impl_bcg_s", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          iw_(m,M1_IW_KS,k,j,i) = iw_(m,M1_IW_KR,k,j,i) - al*iw_(m,M1_IW_KV,k,j,i);
        });
        ImplicitPrecond(M1_IW_KS, M1_IW_KZ);
        ImplicitApplyOp(M1_IW_KZ, M1_IW_KTT);
        Real ts = 0.0, tt2 = 0.0;
        Kokkos::parallel_reduce("m1_impl_bcg_ts", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
          ls += iw_(m,M1_IW_KTT,k,j,i)*iw_(m,M1_IW_KS,k,j,i);
        }, ts);
        Kokkos::parallel_reduce("m1_impl_bcg_tt", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
          ls += iw_(m,M1_IW_KTT,k,j,i)*iw_(m,M1_IW_KTT,k,j,i);
        }, tt2);
        M1GlobalSum2(ts, tt2);
        bcg_nred += 1.0;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real ow = omega;
        par_for("m1_impl_bcg_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
          iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
        });
        rho = rhon;
        rnorm = 0.0;
        Kokkos::parallel_reduce("m1_impl_bcg_rn2", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
          Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
          lmax = (r > lmax) ? r : lmax;
        }, Kokkos::Max<Real>(rnorm));
        M1GlobalMax(rnorm);
        bcg_nred += 1.0;
        if (rnorm/bscale < tol) {
          // the TRUE residual, which is what the tolerance is about
          ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
          par_for("m1_impl_bcg_true", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
          });
          Real rt = 0.0;
          Kokkos::parallel_reduce("m1_impl_bcg_rnt", rng,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
            Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
            lmax = (r > lmax) ? r : lmax;
          }, Kokkos::Max<Real>(rt));
          M1GlobalMax(rt);
          bcg_nred += 1.0;
          if (rt/bscale < tol) {
            done = true;
          } else {
            breakdown = true;   // restart the recurrence from the true residual
          }
        }
        if (!done && !(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        // FALL BACK: one line-Jacobi update, which is always available -- its right-hand
        // side is b minus the lagged off-diagonal term, i.e. the M1_IW_TR the assembly
        // produced, rebuilt here because the preconditioner has overwritten it.
        fell_back = true;
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      par_for("m1_impl_bcg_rs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
      });
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }

  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitBiCGStabEnd
//! \brief the common end of both BiCGStab loops: the line-Jacobi fallback or the copy of
//! the iterate into M1_IW_S2, and the inner-iteration statistics.

void RadiationM1::ImplicitBiCGStabEnd(int nit, bool fell_back) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  if (fell_back) {
    bcg_nfall += 1.0;
    const bool thrd = trans_x3;
    // implicit_vimp: the operator part outside the row (M1VimpRow) is lagged at E^k
    // exactly like the transverse row coefficients, or the fallback solves a different
    // system from the Krylov solve (measured: a Picard 2-cycle, runs_3v_vimplicit)
    const bool vim = vimp_now;
    const int ivb = iw_vimp;
    const bool mus = muscl_now;   // xthinfix-1009 Fix B (as M1VimpRow, lagged at E^k)
    const int imb = iw_muscl;
    const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
    par_for("m1_impl_bcg_lj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real r = iw_(m,M1_IW_KB,k,j,i)
               - iw_(m,M1_IW_CJM,k,j,i)*iw_(m,M1_IW_EP,k,j-1,i)
               - iw_(m,M1_IW_CJP,k,j,i)*iw_(m,M1_IW_EP,k,j+1,i);
      if (thrd) {
        r -= iw_(m,M1_IW_CKM,k,j,i)*iw_(m,M1_IW_EP,k-1,j,i)
             + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
      }
      if (vim) {r -= M1VimpRow(iw_, ivb, M1_IW_EP, m, k, j, i, is, ie, cyclic, thrd);}
      if (mus) {r -= M1MusclRow(iw_, imb, M1_IW_EP, m, k, j, i, is, ie, cyclic, thrd);}
      iw_(m,M1_IW_TR,k,j,i) = r;
    });
    // MILESTONE 3b phase D: the right-hand side of a LINE-JACOBI update is the assembled
    // TR, which under implicit_offdiag = operator still carries -L_off(E^k) -- and the
    // caller has just added L_off(E^k) to b.  Take it out again, so that the fallback is
    // the same update it was before phase D.
    if (od_now == M1_OD_OPERATOR) {
      ImplicitOffDiagOp(M1_IW_EP, M1_IW_TR, -1.0);
    }
    if (VlatOp()) {VetLatOp(M1_IW_EP, M1_IW_TR, -1.0);}
    ImplicitTridiagSolve();
  } else {
    par_for("m1_impl_bcg_out", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_S2,k,j,i) = iw_(m,M1_IW_KX,k,j,i);
    });
  }
  bcg_nsolve += 1.0;
  bcg_itsum += static_cast<Real>(nit);
  bcg_itmax = std::max(bcg_itmax, static_cast<Real>(nit));
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabFused
//! \brief implicit_bcg_sync = 1 | 2: the SAME right-preconditioned BiCGStab recurrence as
//! ImplicitBiCGStab (same breakdown, restart, true-residual and fallback logic), with its
//! host synchronisations and kernel launches cut.  On the GPU every reduction into a
//! host scalar blocks the host until the queue is empty, and the next launch then starts
//! on an idle device (measured: ~40 us per blocking reduction).  Per iteration:
//!  * rho_{k+1} = (rhat, r_{k+1}) and max|r_{k+1}| come out of the kernel that updates x
//!    and r, i.e. ONE reduction where the original has three kernels and two syncs
//!    (upd, max|r|, and (rhat,r) at the top of the next iteration);
//!  * (t,s) and (t,t) come out of one kernel;
//!  * the vector updates p and s also stage the preconditioner input M1_IW_TR, which
//!    removes the copy kernel of ImplicitPrecond;
//!  * sync level 2 on ONE rank: rhat.v is reduced into a device scalar and s = r - alpha
//!    v reads alpha from there, so rhat.v does not block; the host learns it from the
//!    (t,s),(t,t) reduction and only then applies the |rhat.v| breakdown test.  If that
//!    test fails, s/z/t are discarded and the recurrence restarts from x, which the
//!    iteration had not touched -- exactly the original branch.
//! Blocking reductions per iteration: 5 kernels / 4 MPI_Allreduce (level 0), 3 / 3
//! (level 1), 2 / 0 (level 2, one rank).  The recurrence is identical in exact
//! arithmetic; the sums are ordered differently, so the iterates agree to round-off only.

int RadiationM1::ImplicitBiCGStabFused(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace> pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  using HRed = M1BcgRed<Kokkos::HostSpace>;
  const bool devrv = (impl_bcg_sync == 2) && (global_variable::nranks == 1);
  const int kf = impl_kfuse;   // implicit_krylov_fuse (needs bcg_sync = 1, pcr)
  if (impl_stencil) {ImplicitStencilBuild();}   // the operator of this pass, once
  if (impl_stencil && impl_dump_cyc >= 0) {ImplicitDumpOp();}
  if (impl_prec >= 3) {ImplicitMGBuild();}      // the coarse rows of this pass
  if (kf == 3 && impl_kpipe) {return ImplicitBiCGStabPipe(rhsmax);}
  if (kf == 3 && ImplicitKrylovDevOK()) {return ImplicitBiCGStabDev(rhsmax);}
  if (kf == 3) {return ImplicitBiCGStabTwo(rhsmax);}
  auto rvd_ = bcg_rvd;
  // implicit_det_reduce: every sum of the loop in a fixed order (M1DetReduce) and the
  // cross-rank combination in rank order: bitwise run-to-run for a fixed decomposition
  const bool det = impl_det;
  const int ntot = (nmb1 + 1)*nkji;
  auto dpart = det_part;
  M1BcgVal red;
  auto lred = [&](const char *nm, const auto &fn) {   // this rank's part
    if (det) {
      M1DetReduce(nm, ntot, fn, dpart, red);
    } else {
      Kokkos::parallel_reduce(nm, pol, fn, HRed(red));
    }
  };
  auto bred = [&](const char *nm, const auto &fn) {   // ... and over all ranks
    lred(nm, fn);
    M1GlobalBcg(red, det);
  };
  // implicit_lin_cnorm > 0: the max norm is taken of r_i/(s_i E^k_i), s_i = 1 + SRCB_i
  // the row EXCESS of the M-matrix (diagonal minus the off-diagonal moduli: the
  // transport rows sum to zero, the absorption does not), instead of r_i/max|b|.  For
  // an M-matrix A s >= s componentwise, so |dE_i| <= max_j |r_j|/s_j: the norm bounds
  // the error of E in every cell relative to the local E; the test is r < lin_cnorm.
  const bool cn = (impl_lin_cnorm > 0.0);
  const Real efl = e_floor;

  // x0 = the Picard iterate; r0 = b - A x0, with max|r0| and (r0,r0) in the same kernel
  par_for("m1_impl_bcgf_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  bred("m1_impl_bcgf_r0",
  KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
    int m, k, j, i;
    M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
    k += ks; j += js; i += is;
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
    v.s0 += r*r;
    Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                            *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
    v.mx = (a > v.mx) ? a : v.mx;
  });
  bcg_nred += 1.0;
  Real rnorm = red.mx;
  Real rhon = red.s0;   // (rhat, r) of the NEXT iteration, always known on entry
  const Real rnorm0 = red.mx;   // implicit_bcg_fallback = best: the bar to beat
  bcg_r0rel = cn ? rnorm : (rnorm/bscale);
  // Eisenstat-Walker (implicit_lin_ew_max > 0): max|r0| IS the nonlinear residual of
  // the Picard iterate in the max norm (the system was re-linearised about it), so the
  // forcing term needs nothing that is not already here.  Off: the fixed test, as is.
  const bool ew = (impl_ew_max > 0.0) && !ew_tight;   // ew_tight: runs_4a_accel
  Real tabs = cn ? impl_lin_cnorm : (tol*bscale);
  if (ew) {
    Real eta = impl_ew_max;
    if (ew_fprev > 0.0) {
      eta = impl_ew_gam*SQR(rnorm/ew_fprev);
      const Real sg = impl_ew_gam*SQR(ew_etaprev);
      if (sg > 0.1) {eta = fmax(eta, sg);}
      eta = fmin(eta, impl_ew_max);
    }
    ew_fprev = rnorm;
    ew_etaprev = eta;
    tabs = fmax(tabs, eta*rnorm);
  }
  auto lin_done = [=](const Real r) {
    return (ew || cn) ? (r < tabs) : (r/bscale < tol);
  };

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      const Real bt = beta, om = omega;
      Real rv = 0.0;
      bool rvdone = false;
      if (kf > 0) {
        // implicit_krylov_fuse: the p update rides in the preconditioner's load phase
        ImplicitPrecondX(-1, M1_IW_KY, 1, bt, om);
        if (kf == 2) {
          Real o4[4];
          ImplicitHaloOp(M1_IW_KY, M1_IW_KV, 1, o4);
          rv = o4[0];
          Real d2 = 0.0;
          M1GlobalSum2(rv, d2);
          bcg_nred += 1.0;
          rvdone = true;
        } else {
          ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
        }
      } else {
      par_for("m1_impl_bcgf_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real p = iw_(m,M1_IW_KR,k,j,i)
                 + bt*(iw_(m,M1_IW_KP,k,j,i) - om*iw_(m,M1_IW_KV,k,j,i));
        iw_(m,M1_IW_KP,k,j,i) = p;
        iw_(m,M1_IW_TR,k,j,i) = p;
      });
      ImplicitPrecond(-1, M1_IW_KY);
      ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
      }
      auto rvf = KOKKOS_LAMBDA(const int idx, Real &ls) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KV,k,j,i);
      };
      if (devrv) {
        // no host sync: alpha stays on the device until the (t,s) reduction
        Kokkos::parallel_reduce("m1_impl_bcgf_rv", pol, rvf,
                                Kokkos::Sum<Real, DevMemSpace>(rvd_));
      } else {
        if (!rvdone) {
        if (det) {
          M1BcgVal rr;
          M1DetReduce("m1_impl_bcgf_rv", ntot,
                      KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {rvf(idx, v.s0);},
                      dpart, rr);
          rv = rr.s0;
        } else {
          Kokkos::parallel_reduce("m1_impl_bcgf_rv", pol, rvf, rv);
        }
        Real d2 = 0.0;
        M1GlobalSum2(rv, d2, det);
        bcg_nred += 1.0;
        }
        if (!(fabs(rv) > M1_BCG_EPS)) {
          breakdown = true;
        } else {
          alpha = rhon/rv;
        }
      }
      if (!breakdown) {
        const Real al = alpha, rh = rhon;
        const bool dv = devrv;
        if (kf > 0) {
          ImplicitPrecondX(-1, M1_IW_KZ, 2, al, 0.0);
        } else {
        par_for("m1_impl_bcgf_s", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          Real a = dv ? (rh/rvd_()) : al;
          Real s = iw_(m,M1_IW_KR,k,j,i) - a*iw_(m,M1_IW_KV,k,j,i);
          iw_(m,M1_IW_KS,k,j,i) = s;
          iw_(m,M1_IW_TR,k,j,i) = s;
        });
        ImplicitPrecond(-1, M1_IW_KZ);
        }
        if (kf == 2) {
          red.s2 = 0.0;
          red.mx = 0.0;
          Real o4[4];
          ImplicitHaloOp(M1_IW_KZ, M1_IW_KTT, 2, o4);
          red.s0 = o4[0];
          red.s1 = o4[1];
        } else {
        ImplicitApplyOp(M1_IW_KZ, M1_IW_KTT);
        lred("m1_impl_bcgf_ts",
        KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
          int m, k, j, i;
          M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
          k += ks; j += js; i += is;
          Real t = iw_(m,M1_IW_KTT,k,j,i);
          v.s0 += t*iw_(m,M1_IW_KS,k,j,i);
          v.s1 += t*t;
          if (dv && idx == 0) {v.s2 += rvd_();}   // carries rhat.v to the host, exactly
        });
        }
        M1GlobalBcg(red, det);
        bcg_nred += 1.0;
        if (devrv) {
          rv = red.s2;
          if (!(fabs(rv) > M1_BCG_EPS)) {
            breakdown = true;
          } else {
            alpha = rhon/rv;
          }
        }
      }
      if (!breakdown) {
        const Real ts = red.s0, tt2 = red.s1;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real al = alpha, ow = omega;
        bred("m1_impl_bcgf_upd",
        KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
          int m, k, j, i;
          M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
          k += ks; j += js; i += is;
          iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
          Real r = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
          iw_(m,M1_IW_KR,k,j,i) = r;
          v.s0 += iw_(m,M1_IW_KRH,k,j,i)*r;
          Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                                  *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
          v.mx = (a > v.mx) ? a : v.mx;
        });
        bcg_nred += 1.0;
        if (impl_dtrace > 0) {
          const Real sv[5] = {rv, alpha, omega, red.s0, red.mx};
          DetTraceScalars("bcgf", 5, sv);
        }
        rho = rhon;
        rhon = red.s0;
        rnorm = red.mx;
        if (lin_done(rnorm)) {
          // the TRUE residual, which is what the tolerance is about
          ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
          bred("m1_impl_bcgf_true",
          KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
            int m, k, j, i;
            M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
            k += ks; j += js; i += is;
            Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
            iw_(m,M1_IW_KR,k,j,i) = r;
            Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                                    *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
            v.mx = (a > v.mx) ? a : v.mx;
          });
          bcg_nred += 1.0;
          if (lin_done(red.mx)) {
            done = true;
          } else {
            breakdown = true;   // restart the recurrence from the true residual
          }
        }
        if (!done && !(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > impl_bcg_maxrst) {
        fell_back = true;   // one line-Jacobi update (ImplicitBiCGStabEnd)
        if (impl_bcg_keep) {
          // implicit_bcg_fallback = best: keep the Krylov iterate if its TRUE residual
          // beats the initial one (b - A x0, x0 = the Picard iterate)
          ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
          bred("m1_impl_bcgf_keep",
          KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
            int m, k, j, i;
            M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
            k += ks; j += js; i += is;
            Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
            Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                                    *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
            v.mx = (a > v.mx) ? a : v.mx;
          });
          bcg_nred += 1.0;
          if (red.mx < rnorm0) {
            fell_back = false;
            bcg_nkeep += 1.0;
          }
        }
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      bred("m1_impl_bcgf_rs",
      KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
        v.s0 += r*r;
      });
      bcg_nred += 1.0;
      rhon = red.s0;
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitDetOpRed
//! \brief implicit_det_reduce: y = A x (ImplicitHaloOp with no reduction), then the sums
//! the fused kernel would have returned for `red`, in a FIXED order (M1DetReduce):
//! red 1|4: out[0] = (rhat, y), and 4 also out[3] = max|r|; red 2|3: out[0] = (y, s),
//! out[1] = (y, y), and 3 also out[2] = (rhat, y).  This rank only (the caller sums).

void RadiationM1::ImplicitDetOpRed(int xc, int yc, int red, Real *out) {
  ImplicitHaloOp(xc, yc, 0, nullptr);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.ie - is + 1;
  const int nji = (indcs.je - js + 1)*ni;
  const int nkji = (indcs.ke - ks + 1)*nji;
  const int ntot = pmy_pack->nmb_thispack*nkji;
  auto iw_ = iw;
  const int cy = yc, rm = red;
  M1BcgVal v0;
  M1DetReduce("m1_det_opred", ntot, KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
    int m, k, j, i;
    M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
    k += ks; j += js; i += is;
    const Real y = iw_(m,cy,k,j,i);
    if (rm == 1 || rm == 4) {
      v.s0 += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        const Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        v.mx = (a > v.mx) ? a : v.mx;
      }
    } else {
      v.s0 += y*iw_(m,M1_IW_KS,k,j,i);
      v.s1 += y*y;
      if (rm == 3) {v.s2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, det_part, v0);
  out[0] = v0.s0;
  out[1] = v0.s1;
  out[2] = v0.s2;
  out[3] = v0.mx;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabTwo
//! \brief implicit_krylov_fuse = 3: the right-preconditioned BiCGStab of
//! ImplicitBiCGStabFused with TWO blocking reductions per iteration instead of three.
//!  * (rhat, v) comes out of the operator kernel v = A y, together with max|r| of the
//!    CURRENT r (the one the previous iteration made): the convergence test is taken
//!    one half-iteration late, and on success that half-iteration (p, y, v; x and r are
//!    untouched by it) is discarded;
//!  * (t,s), (t,t) and (rhat,t) come out of the operator kernel t = A z;
//!  * the update x += alpha y + omega z, r = s - omega t is then a plain kernel, and
//!    rho_{k+1} = (rhat, r_{k+1}) = rho_k - alpha (rhat,v) - omega (rhat,t) follows by
//!    recurrence (exact in exact arithmetic, since (rhat,s) = rho_k - alpha (rhat,v)).
//! Breakdown, restart, true-residual and fallback rules are those of the fused loop.
//! Same iterates in exact arithmetic; round-off differs (the recurrence for rho).

int RadiationM1::ImplicitBiCGStabTwo(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  using HRed = M1BcgRed<Kokkos::HostSpace>;
  // implicit_det_reduce: fixed-order sums (M1DetReduce), rank-ordered MPI, and the sums
  // the fused operator kernels return taken by ImplicitDetOpRed instead
  const bool det = impl_det;
  const int ntot = (nmb1 + 1)*nkji;
  auto dpart = det_part;
  auto lred = [&](const char *nm, const auto &fn, M1BcgVal &rv) {
    if (det) {
      M1DetReduce(nm, ntot, fn, dpart, rv);
    } else {
      Kokkos::parallel_reduce(nm, pol, fn, HRed(rv));
    }
    M1GlobalBcg(rv, det);
  };
  auto hop = [&](const int xc, const int yc, const int rd, Real *o) {
    if (det) {
      ImplicitDetOpRed(xc, yc, rd, o);
    } else {
      ImplicitHaloOp(xc, yc, rd, o);
    }
  };

  // x0 = the Picard iterate; r0 = b - A x0, with max|r0| and (r0,r0) in the same kernel
  par_for("m1_impl_bcg2_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  M1BcgVal red;
  lred("m1_impl_bcg2_r0",
  KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
    int m, k, j, i;
    M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
    k += ks; j += js; i += is;
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
    v.s0 += r*r;
    Real a = fabs(r);
    v.mx = (a > v.mx) ? a : v.mx;
  }, red);
  bcg_nred += 1.0;
  Real rnorm = red.mx;
  Real rhon = red.s0;
  bcg_r0rel = rnorm/bscale;
  const bool ew = (impl_ew_max > 0.0) && !ew_tight;   // ew_tight: runs_4a_accel
  Real tabs = tol*bscale;
  if (ew) {
    Real eta = impl_ew_max;
    if (ew_fprev > 0.0) {
      eta = impl_ew_gam*SQR(rnorm/ew_fprev);
      const Real sg = impl_ew_gam*SQR(ew_etaprev);
      if (sg > 0.1) {eta = fmax(eta, sg);}
      eta = fmin(eta, impl_ew_max);
    }
    ew_fprev = rnorm;
    ew_etaprev = eta;
    tabs = fmax(tabs, eta*rnorm);
  }
  auto lin_done = [=](const Real r) {
    return ew ? (r < tabs) : (r/bscale < tol);
  };
  // the true residual of x, which is what the tolerance is about
  auto true_ok = [&]() -> bool {
    ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
    M1BcgVal tr;
    lred("m1_impl_bcg2_true",
    KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
      int m, k, j, i;
      M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
      k += ks; j += js; i += is;
      Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = r;
      Real a = fabs(r);
      v.mx = (a > v.mx) ? a : v.mx;
    }, tr);
    bcg_nred += 1.0;
    return lin_done(tr.mx);
  };

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  bool pend = false;   // r was updated by the last iteration; its max is not known yet
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      ImplicitPrecondX(-1, M1_IW_KY, 1, beta, omega);
      Real o4[4];
      hop(M1_IW_KY, M1_IW_KV, 4, o4);
      M1BcgVal a;
      a.s0 = o4[0]; a.s1 = 0.0; a.s2 = 0.0; a.mx = o4[3];
      M1GlobalBcg(a, det);
      bcg_nred += 1.0;
      if (pend) {
        pend = false;
        if (lin_done(a.mx)) {
          --nit;   // this half-iteration is discarded: x and r are those it started from
          if (true_ok()) {
            done = true;
            break;
          }
          breakdown = true;   // restart the recurrence from the true residual
        }
      }
      const Real rv = a.s0;
      if (!breakdown && !(fabs(rv) > M1_BCG_EPS)) {breakdown = true;}
      if (!breakdown) {
        alpha = rhon/rv;
        ImplicitPrecondX(-1, M1_IW_KZ, 2, alpha, 0.0);
        hop(M1_IW_KZ, M1_IW_KTT, 3, o4);
        red.s0 = o4[0]; red.s1 = o4[1]; red.s2 = o4[2]; red.mx = 0.0;
        M1GlobalBcg(red, det);
        bcg_nred += 1.0;
        const Real ts = red.s0, tt2 = red.s1, rt = red.s2;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real al = alpha, ow = omega;
        rho = rhon;
        if (impl_rho_direct) {
          // implicit_bcg_rho_direct (runs_5p_coarse2): (rhat, r) summed directly in the
          // update kernel (one more blocking reduction) instead of the recurrence below,
          // whose cancellation stalls the solve at ~1e-10 on hard systems
          M1BcgVal ur;
          lred("m1_impl_bcg2_updr",
          KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
            int m, k, j, i;
            M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
            k += ks; j += js; i += is;
            iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
            const Real r = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
            iw_(m,M1_IW_KR,k,j,i) = r;
            v.s0 += iw_(m,M1_IW_KRH,k,j,i)*r;
          }, ur);
          bcg_nred += 1.0;
          rhon = ur.s0;
        } else {
          par_for("m1_impl_bcg2_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
            iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
          });
          rhon = rhon - alpha*rv - omega*rt;
        }
        pend = true;
        if (!(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      pend = false;
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        fell_back = true;
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      lred("m1_impl_bcg2_rs",
      KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
        v.s0 += r*r;
      }, red);
      bcg_nred += 1.0;
      rhon = red.s0;
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }
  // the cap reached with an update whose max is not known yet: test it once
  if (!done && !fell_back && pend) {
    if (true_ok()) {done = true;}
  }
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::SetImplicitX1BC
//! \brief let a problem generator name the x1 boundary types (and the imposed fluxes) of
//! the implicit solve, for a mesh whose x1 flags are `user` and whose own BC routine the
//! module cannot interpret.  Alternative to <rad_m1>/implicit_bc_x1min|max.

void RadiationM1::SetImplicitX1BC(int lo_type, Real lo_flux, int hi_type, Real hi_flux) {
  ibc_x1min = lo_type;
  iflux_x1min = lo_flux;
  ibc_x1max = hi_type;
  iflux_x1max = hi_flux;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::TmrMark
//! \brief implicit_timers: fence, charge the time since the last mark to category c
//! (c < 0 only sets the mark)

void RadiationM1::TmrMark(int c) {
  if (!tmr_on) return;
  Kokkos::fence();
  const double t = tmr_t.seconds();
  if (c >= 0) {tmr_acc[c] += t - tmr_last;}
  tmr_last = t;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::OnePassAuto
//! \brief implicit_one_pass_auto: the counters of one kind of solve (0 be, 1 stage 1,
//! 2 stage 2) after an eligible solve; on = one_pass was on for it, one = it was accepted
//! after one pass.  Host-only bookkeeping, identical on every rank (see ImplicitSolve).

void RadiationM1::OnePassAuto(const int t, const bool on, const bool one) {
  const char *knm[3] = {"be", "stage 1", "stage 2"};
  const bool r0 = (global_variable::my_rank == 0);
  if (!on) {
    onep_actr[t] += 1.0;
    return;
  }
  if (one) {
    onep_actr[t] = 0.0;
    if (onep_prb[t] > 0.5) {
      onep_prb[t] = 0.0;
      impl_onep_nren += 1.0;
      if (r0) {
        std::cout << "<rad_m1> implicit_one_pass_auto: " << knm[t] << " solves: re-probe "
                  << "accepted, one_pass ON again (cycle " << pmy_pack->pmesh->ncycle
                  << ")" << std::endl;
      }
    }
    return;
  }
  onep_actr[t] += 1.0;
  if (onep_actr[t] >= static_cast<Real>(impl_onep_awin*impl_onep)) {
    onep_off[t] = 1.0;
    onep_actr[t] = 0.0;
    if (onep_prb[t] > 0.5) {
      // a failed re-probe: counted in the summary only (it recurs every auto_reprobe)
      onep_prb[t] = 0.0;
      impl_onep_nfpr += 1.0;
    } else {
      impl_onep_ndis += 1.0;
      if (r0) {
        std::cout << "<rad_m1> implicit_one_pass_auto: " << knm[t] << " solves: none "
                  << "accepted after one pass in " << impl_onep_awin*impl_onep
                  << " solves, one_pass OFF (re-probe every " << impl_onep_arep
                  << " solves; cycle " << pmy_pack->pmesh->ncycle << ")" << std::endl;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitReport
//! \brief one line at the end of the run with the Picard statistics

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::EventTotals
//! \brief fixbundle-1009 F4: the solver / limiter totals of this run (or restart) for the
//! event log, all ranks (collective).  v[0..6] are decided on global reductions and are
//! the same on every rank; v[7..12] are per-rank device counters, summed here.
//!   0 BiCGStab breakdowns, 1 line-Jacobi fallbacks, 2 fallbacks that kept the Krylov
//!   iterate (implicit_bcg_fallback = best), 3 APPLIED gas-Newton fallbacks,
//!   4 solves switched to the gas root find, 5 od positivity drops, 6 vet_col_lat
//!   positivity drops, 7 opac-Newton face terms dropped by the guard, 8 E floor raises
//!   (write-back, cell-solves), 9 energy the floor CREATED (x volume), 10 of which the
//!   solve undershoot (implicit_pos_floor_solve), 11 ApplyClosureLimits E raises,
//!   12 ApplyClosureLimits |F| clips (active cells, every stage)

void RadiationM1::EventTotals(Real *v) {
  v[0] = bcg_nbreak;
  v[1] = bcg_nfall;
  v[2] = bcg_nkeep;
  v[3] = newt_nfb;
  v[4] = impl_gn_nsw;
  v[5] = od_nfall;
  v[6] = vlat_nfall;
  for (int q = 7; q < 13; ++q) {v[q] = 0.0;}
  if (opn_nskip_d.extent_int(0) > 0) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    v[7] = h(0);
  }
  if (pos_cnt_d.extent_int(0) >= M1_POS_N) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    v[8] = h(M1_POS_FLR);
    v[9] = h(M1_POS_FLR_UN);
    v[10] = h(M1_POS_FLRS);
  }
  if (acl_cnt_d.extent_int(0) == 2) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acl_cnt_d);
    v[11] = h(0);
    v[12] = h(1);
  }
#if MPI_PARALLEL_ENABLED
  Real g[6];
  MPI_Allreduce(&v[7], g, 6, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < 6; ++q) {v[7+q] = g[q];}
#endif
}

void RadiationM1::ImplicitReport() {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;
  if (cs_geom && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> cubed sphere: seam face averages " << cs_seam_avg_n
              << (cs_seam_avg_on ? "" : " (disabled by rad_m1/cs_seam_avg)")
              << ", largest change of an x2 face value (rank 0) " << cs_seam_dmax
              << " of max|F0|" << std::endl;
  }
  // implicit_opac_newton_guard: the dropped (row, face) pairs of all ranks (collective:
  // every rank reaches this line with the same impl_opac_newton and guard)
  const bool opgr = impl_opac_newton && (impl_opn_guard > 0.0);
  if (opgr) {
    auto hsk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), opn_nskip_d);
    opn_nskip = hsk(0);
#if MPI_PARALLEL_ENABLED
    Real g = 0.0;
    MPI_Allreduce(&opn_nskip, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    opn_nskip = g;
#endif
  }
  if (impl_g0_lim > 0.0 || impl_g0_exch || impl_pos_gas || impl_pos_floor) {
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = hpc(q);}
#if MPI_PARALLEL_ENABLED
    Real g[M1_POS_N];
    MPI_Allreduce(pos_cnt, g, M1_POS_N, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int q = 0; q < M1_POS_N; ++q) {pos_cnt[q] = g[q];}
#endif
  }
  // the Picard iteration count is MPI_MAX-reduced every step (see ImplicitSolve), so
  // every rank holds the same three numbers and no reduction is needed here
  if (global_variable::my_rank != 0) return;
  if (tmr_c0 > 0) {
    // implicit_timers: seconds per category, then the counters (see TmrMark)
    std::cout << "<rad_m1> timers from cycle " << tmr_c0 << ": stages=" << tmr_cnt[7]
              << " solves=" << tmr_cnt[0] << " passes=" << tmr_cnt[1]
              << " kry_it=" << tmr_cnt[2] << " passes_s1=" << tmr_cnt[3]
              << " passes_s2=" << tmr_cnt[4] << " it_s1=" << tmr_cnt[5]
              << " it_s2=" << tmr_cnt[6] << std::endl;
    const char *nm[12] = {"closure", "opacity", "solve_pre", "tensor", "pred",
                          "pass_setup", "krylov", "pass_post", "solve_end", "hyd_c2p",
                          "m1_bvals", "unused"};
    std::cout << "<rad_m1> timers(s):";
    for (int q = 0; q < 11; ++q) {std::cout << " " << nm[q] << "=" << tmr_acc[q];}
    std::cout << std::endl;
  }
  MRReport();
  Real mean = (impl_nstep > 0.0) ? (impl_itsum/impl_nstep) : 0.0;
  std::cout << "<rad_m1> implicit transport: solves=" << impl_nstep
            << " Picard iterations mean=" << mean << " max=" << impl_itmax
            << " NON-CONVERGED=" << impl_nfail << std::endl;
  if (impl_gn_sw > 0) {
    std::cout << "<rad_m1> implicit_gas_newton_switch=" << impl_gn_sw << " (min pass "
              << impl_gn_sw_min << "): solves switched to the root find=" << impl_gn_nsw
              << std::endl;
  }
  if (impl_res_mask) {
    std::cout << "<rad_m1> implicit_res mask: rho < " << impl_res_dmin << " or r > "
              << impl_res_rmax << " (0 = unused) left out of the stopping test; max"
              << " excluded residual at the last pass=" << impl_res_exmax
              << ", steps whose excluded residual >= tol=" << impl_res_nex << std::endl;
  }
  if (opgr) {
    std::cout << "<rad_m1> implicit_opac_newton_guard=" << impl_opn_guard
              << ": Newton face terms dropped (all ranks, all passes)=" << opn_nskip
              << " per solve=" << ((impl_nstep > 0.0) ? (opn_nskip/impl_nstep) : 0.0)
              << std::endl;
  }
  if (impl_g0_lim > 0.0 || impl_g0_exch || impl_pos_gas || impl_pos_floor) {
    std::cout << "<rad_m1> m1-positivity: implicit_g0_exchange=" << impl_g0_exch
              << " implicit_g0_limit=" << impl_g0_lim
              << " clipped cell-passes=" << pos_cnt[M1_POS_G0]
              << " | implicit_pos_gas=" << impl_pos_gas << " (frac "
              << impl_pos_gas_frac << ") cell-solves=" << pos_cnt[M1_POS_GAS]
              << " energy moved rad->gas=" << pos_cnt[M1_POS_GAS_DE]
              << " | E floor raises cell-solves=" << pos_cnt[M1_POS_FLR]
              << " energy from gas=" << pos_cnt[M1_POS_FLR_DE]
              << " energy created=" << pos_cnt[M1_POS_FLR_UN] << " (code units x volume)"
              << " of which solve undershoot (pos_floor_solve)=" << pos_cnt[M1_POS_FLRS]
              << " | |F| > c E scaled back cell-solves=" << pos_cnt[M1_POS_FCLIP]
              << " sum(|F|/(cE) - 1)=" << pos_cnt[M1_POS_FCLIPM] << std::endl;
  }
  if (impl_onep > 0) {
    std::cout << "<rad_m1> implicit_one_pass: period=" << impl_onep
              << " safety=" << impl_onep_s
              << " solves accepted after one pass=" << impl_onep_n
              << " contraction measurements=" << impl_onep_nchk
              << " last q (be, stage 1, stage 2)=" << onep_qa[0] << " " << onep_qa[1]
              << " " << onep_qa[2] << std::endl;
    if (impl_onep_auto) {
      // the counters are those of this run (from its start or restart); the state is
      // per kind: on, off, or on as a re-probe not yet confirmed
      std::cout << "<rad_m1> implicit_one_pass_auto: window=" << impl_onep_awin
                << " reprobe=" << impl_onep_arep << "; this run: switch-offs="
                << impl_onep_ndis << " re-probes=" << impl_onep_nprb
                << " (failed=" << impl_onep_nfpr << " confirmed=" << impl_onep_nren
                << "); state at end (be, stage 1, stage 2)=";
      for (int t = 0; t < 3; ++t) {
        const char *st = (onep_off[t] > 0.5) ? "off"
                         : ((onep_prb[t] > 0.5) ? "probe" : "on");
        std::cout << st << ((t < 2) ? "," : "");
      }
      std::cout << std::endl;
    }
  }
  if (impl_accel == M1_IACC_ANDERSON) {
    Real apst = (impl_nstep > 0.0) ? (aa_nacc/impl_nstep) : 0.0;
    std::cout << "<rad_m1> anderson: m=" << impl_and_m
              << " beta=" << impl_and_beta
              << " start=" << impl_and_start
              << " accelerated passes=" << aa_nacc
              << " (" << apst << " per solve)"
              << " history restarts=" << aa_nrst << std::endl;
  }
  if (trans_on) {
    Real lmean = (impl_nstep > 0.0) ? (impl_linsum/impl_nstep) : 0.0;
    std::cout << "<rad_m1> implicit transverse ("
              << (bicg_on ? "bicgstab" : "line_jacobi")
              << "): final 7-point linear "
              << "residual (max|r|/max|b|) mean=" << lmean << " max=" << impl_linmax
              << std::endl;
    // vgdspeed-1009: say which test stopped the solves.  The max|r|/max|b| above is
    // compared with implicit_lin_tol only when implicit_lres_test is on (default off for
    // the fixed closures: eddington, vet_sc, tau / vet_col); the inner BiCGStab stops on
    // implicit_lin_cnorm (per-cell |r|/(s E)) when > 0, else on |r|/max|b| < lin_tol,
    // loosened by the Eisenstat-Walker term when implicit_lin_ew_max > 0
    std::cout << "<rad_m1>   (diagnostic; Picard test on it: "
              << (impl_lres_test ? "yes, lin_tol=" : "no; lin_tol=") << impl_lin_tol
              << "; inner stop: "
              << ((impl_lin_cnorm > 0.0) ? "cnorm=" : "lin_tol relative to max|b|");
    if (impl_lin_cnorm > 0.0) {std::cout << impl_lin_cnorm;}
    if (impl_ew_max > 0.0) {
      std::cout << ", loosened by Eisenstat-Walker ew_max=" << impl_ew_max;
    }
    std::cout << ")" << std::endl;
  }
  if (bicg_on) {
    Real imean = (bcg_nsolve > 0.0) ? (bcg_itsum/bcg_nsolve) : 0.0;
    Real rper = (bcg_itsum > 0.0) ? (bcg_nred/bcg_itsum) : 0.0;
    std::cout << "<rad_m1> bicgstab: outer passes=" << impl_itsum
              << " linear solves=" << bcg_nsolve
              << " inner iterations mean=" << imean << " max=" << bcg_itmax
              << " total=" << bcg_itsum << std::endl;
    std::cout << "<rad_m1> bicgstab: breakdowns=" << bcg_nbreak
              << " line_jacobi fallbacks=" << bcg_nfall
              << " global reductions=" << bcg_nred
              << " (" << rper << " per inner iteration)" << std::endl;
    if (impl_bcg_keep || impl_bcg_maxrst != 2) {
      std::cout << "<rad_m1> bicgstab: implicit_bcg_max_restarts=" << impl_bcg_maxrst
                << " fallback=" << (impl_bcg_keep ? "best" : "line")
                << ": Krylov iterate kept instead of line-Jacobi=" << bcg_nkeep
                << std::endl;
    }
    if (impl_kdev > 0) {
      std::cout << "<rad_m1> implicit_krylov_dev: period=" << impl_kdev
                << " host status reads=" << kdev_nchk << " queued iterations="
                << kdev_nq << std::endl;
    }
  }
  if (impl_gas_newton || impl_eos_cache) {
    Real fpc = (gas_ncell > 0.0) ? (newt_nfb/gas_ncell) : 0.0;
    Real mpc = (gas_ncell > 0.0) ? (ec_nmiss/gas_ncell) : 0.0;
    std::cout << "<rad_m1> gas coupling: newton="
              << (impl_gas_newton ? "true" : "false")
              << " eos_cache=" << (impl_eos_cache ? "true" : "false")
              << " cell-passes=" << gas_ncell
              << " newton fallbacks=" << newt_nfb << " (" << fpc << " per cell-pass)"
              << std::endl;
    if (impl_eos_cache) {
      std::cout << "<rad_m1> eos_cache: nt=" << impl_ecnt
                << " misses=" << ec_nmiss << " (" << mpc << " per cell-pass)"
                << " max |de|/e vs the table=" << ec_emax
                << " max |dq|/q of the exchanged energy=" << ec_tmax << std::endl;
    }
  }
  if (impl_pcr_check) {
    std::cout << "<rad_m1> line solver=" << ((impl_line_solver == 1) ? "pcr" : "thomas")
              << " pcr_check: calls=" << pcr_chk_n
              << " max over calls of max|x_pcr - x_thomas|/max|x|=" << pcr_chk_max
              << " (rank 0)" << std::endl;
  }
  if (trans_on) {
    std::cout << "<rad_m1> offdiag="
              << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                  ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
              << " closure_relax=" << impl_crelax
              << " closure_lag=" << (impl_clag_step ? "step" : "pass")
              << " positivity fallbacks=" << od_nfall
              << " (vet_col_lat " << vlat_nfall << ")"
              << " min E from the solve=" << od_emin << std::endl;
  }
  if (impl_vimp) {
    std::cout << "<rad_m1> implicit_vimp positivity fallbacks=" << vimp_nfall
              << " min E from the solve=" << vimp_emin << std::endl;
  }
  if (impl_muscl) {
    std::cout << "<rad_m1> implicit_hr_recon = plm positivity fallbacks=" << muscl_nfall
              << " min E from the solve=" << muscl_emin << std::endl;
  }
  std::cout << "<rad_m1> floor clips (all ranks, every solve incl. BE): solved E <= "
            << "e_floor cell-solves=" << flr_ne
            << " written-back gas eint <= 0 cell-solves="
            << flr_ng << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPicardLog
//! \brief DIAGNOSTIC (<rad_m1>/implicit_picard_log): one line per Picard pass with the
//! E and T parts of the Picard residual (and the x1 index of the cell that owns each),
//! the transverse change lresid, the inner iterations and the inner starting residual.

void RadiationM1::ImplicitPicardLog(int it, int nin, Real resid, Real lresid, bool srct) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  using MaxLoc = Kokkos::MaxLoc<Real,int>;
  Real vals[3] = {0.0, 0.0, 0.0};
  int locs[3] = {-1, -1, -1};
  const int comp[3] = {M1_IW_S1, M1_IW_S3, M1_IW_LRES};
  for (int q = 0; q < 3; ++q) {
    if (q == 1 && !srct) {continue;}
    if (q == 2 && !trans_on) {continue;}
    if (q == 0 && !srct) {
      // without the gas coupling RES is the E part alone
      locs[0] = 0;
    }
    const int c = (q == 0 && !srct) ? M1_IW_RES : comp[q];
    MaxLoc::value_type mloc;
    Kokkos::parallel_reduce("m1_impl_plog",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                  MaxLoc::value_type &lmx) {
      Real r = iw_(m,c,k,j,i);
      if (r > lmx.val) {
        lmx.val = r;
        // stall_1002: the full (m,k,j,i) of the cell, i in the low 10 bits
        lmx.loc = (((m*1024 + k)*1024 + j)*1024) + i;
      }
    }, MaxLoc(mloc));
    vals[q] = mloc.val;
    locs[q] = mloc.loc;
  }
  // stall_1002: every rank prints its own maxima, (m,k,j,i) decoded
  {
    std::string sl[3];
    for (int q = 0; q < 3; ++q) {
      int l = locs[q];
      if (l < 0) {sl[q] = "-"; continue;}
      sl[q] = std::to_string(l/(1024*1024*1024)) + ","
              + std::to_string((l/(1024*1024))%1024) + ","
              + std::to_string((l/1024)%1024) + "," + std::to_string(l%1024);
    }
    std::cout << "<rad_m1> plogR rank=" << global_variable::my_rank << " step="
              << static_cast<int>(impl_nstep) << " pass=" << it << " res=" << resid
              << " resE=" << vals[0] << " @" << sl[0] << " resT=" << vals[1] << " @"
              << sl[1] << " lres=" << vals[2] << " @" << sl[2] << std::endl;
  }
  for (int q = 0; q < 3; ++q) {if (locs[q] >= 0) {locs[q] = locs[q] % 1024;}}
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> plog step=" << static_cast<int>(impl_nstep) << " pass=" << it
              << " res=" << resid << " resE=" << vals[0] << " iE=" << locs[0]
              << " resT=" << vals[1] << " iT=" << locs[1]
              << " lres=" << lresid << " iL=" << locs[2]
              << " nin=" << nin << " r0=" << bcg_r0rel << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitNCDump
//! \brief DEBUG (<rad_m1>/implicit_nc_dump, nc_cure_1002): one line per rank with the
//! local state of the rank's worst-residual cell: its geometry and optical depths, the
//! coupling stiffness c dt rho kappa_P, the gas/radiation energy ratio, v/c, the reduced
//! flux, the lagged closure, the T-Newton row and flags, and E, T of its six neighbours.
//! Prints only.

int RadiationM1::ImplicitNCDump(int it, int tag, bool gasx, int igb, int igr, int igf,
                                int igy, int floc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nk = ke - ks + 1, nj = je - js + 1, ni = ie - is + 1;
  auto iw_ = iw;
  using MaxLoc = Kokkos::MaxLoc<Real,int>;
  MaxLoc::value_type mloc;
  mloc.loc = floc;
  if (floc < 0) {
  Kokkos::parallel_reduce("m1_impl_ncd_loc",
  Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                         {nmb1+1,ke+1,je+1,ie+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                MaxLoc::value_type &lmx) {
    // tag 1: the cell with the most negative solved E (floor clips)
    Real r = (tag == 1) ? -iw_(m,M1_IW_S2,k,j,i) : iw_(m,M1_IW_RES,k,j,i);
    if (r > lmx.val) {
      lmx.val = r;
      lmx.loc = ((m*nk + (k-ks))*nj + (j-js))*ni + (i-is);
    }
  }, MaxLoc(mloc));
  }
  if (mloc.loc < 0) {return -1;}
  const int lm = mloc.loc/(nk*nj*ni);
  int rem = mloc.loc - lm*nk*nj*ni;
  const int lk = rem/(nj*ni) + ks;
  rem -= (lk-ks)*nj*ni;
  const int lj = rem/ni + js;
  const int li = rem - (lj-js)*ni + is;
  constexpr int NV = 48;
  DvceArray1D<Real> dv("ncd", NV);
  auto opac_ = opac;
  auto uh = FluidRef::Get(pmy_pack).u0;
  auto x1f = pmy_pack->pcoord->xx1f;
  auto x1v = pmy_pack->pcoord->x1v;
  auto mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light, dt = dt_sub;
  Kokkos::parallel_for("m1_impl_ncd", Kokkos::RangePolicy<>(DevExeSpace(), 0, 1),
  KOKKOS_LAMBDA(const int) {
    const int m = lm, k = lk, j = lj, i = li;
    const Real r = x1v(m,i);
    const Real dr = x1f(m,i+1) - x1f(m,i);
    const Real rdth = r*mbsize.d_view(m).dx2;
    const Real e = iw_(m,M1_IW_EP,k,j,i);
    const Real f1 = iw_(m,M1_IW_F1,k,j,i), f2 = iw_(m,M1_IW_F2,k,j,i);
    const Real f3 = iw_(m,M1_IW_F3,k,j,i);
    const Real kt = iw_(m,M1_IW_KT,k,j,i);
    dv(0) = r; dv(1) = dr; dv(2) = rdth;
    dv(3) = uh(m,IDN,k,j,i);
    dv(4) = iw_(m,M1_IW_EN,k,j,i); dv(5) = e; dv(6) = iw_(m,M1_IW_TP,k,j,i);
    dv(7) = iw_(m,M1_IW_EGN,k,j,i);
    dv(8) = iw_(m,M1_IW_RES,k,j,i); dv(9) = iw_(m,M1_IW_LRES,k,j,i);
    dv(10) = kt*dr; dv(11) = kt*rdth;                       // tau_r, tau_lat
    dv(12) = cl*dt*opac_(m,M1_OP_P,k,j,i);                   // c dt rho kP
    dv(13) = cl*dt*opac_(m,M1_OP_E,k,j,i);                   // c dt rho kE
    dv(14) = iw_(m,M1_IW_EGN,k,j,i)/fmax(e, 1.0e-300);       // e_gas/E
    dv(15) = iw_(m,M1_IW_V1,k,j,i)/cl; dv(16) = iw_(m,M1_IW_V2,k,j,i)/cl;
    dv(17) = iw_(m,M1_IW_V3,k,j,i)/cl;
    dv(18) = sqrt(f1*f1 + f2*f2 + f3*f3)/(cl*fmax(e, 1.0e-300));
    dv(19) = f1/(cl*fmax(e, 1.0e-300));
    dv(20) = iw_(m,M1_IW_WCHI,k,j,i); dv(21) = iw_(m,M1_IW_RF0,k,j,i);
    dv(22) = iw_(m,M1_IW_G0,k,j,i);
    dv(23) = iw_(m,M1_IW_SRCB,k,j,i); dv(24) = iw_(m,M1_IW_SRCR,k,j,i);
    dv(25) = iw_(m,M1_IW_DE0,k,j,i);
    dv(26) = gasx ? iw_(m,igb,k,j,i) : 0.0; dv(27) = gasx ? iw_(m,igr,k,j,i) : 0.0;
    dv(28) = gasx ? iw_(m,igf,k,j,i) : 0.0; dv(29) = gasx ? iw_(m,igy,k,j,i) : 0.0;
    dv(30) = iw_(m,M1_IW_EP,k,j,i-1); dv(31) = iw_(m,M1_IW_EP,k,j,i+1);
    dv(32) = iw_(m,M1_IW_EP,k,j-1,i); dv(33) = iw_(m,M1_IW_EP,k,j+1,i);
    dv(34) = iw_(m,M1_IW_EP,k-1,j,i); dv(35) = iw_(m,M1_IW_EP,k+1,j,i);
    dv(36) = iw_(m,M1_IW_TP,k,j,i-1); dv(37) = iw_(m,M1_IW_TP,k,j,i+1);
    dv(38) = iw_(m,M1_IW_TP,k,j-1,i); dv(39) = iw_(m,M1_IW_TP,k,j+1,i);
    dv(40) = iw_(m,M1_IW_TP,k-1,j,i); dv(41) = iw_(m,M1_IW_TP,k+1,j,i);
    dv(42) = iw_(m,M1_IW_TDIA,k,j,i); dv(43) = iw_(m,M1_IW_TRHS,k,j,i);
    dv(44) = opac_(m,M1_OP_P,k,j,i); dv(45) = opac_(m,M1_OP_E,k,j,i);
    dv(46) = kt; dv(47) = iw_(m,M1_IW_S2,k,j,i);
  });
  auto hv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dv);
  // implicit_stall_trace: the hottest of the six neighbours inside this MeshBlock
  ncd_hot_loc = -1;
  {
    Real tb = -1.0;
    const int di[6] = {-1, 1, 0, 0, 0, 0}, dj[6] = {0, 0, -1, 1, 0, 0};
    const int dk[6] = {0, 0, 0, 0, -1, 1};
    for (int q = 0; q < 6; ++q) {
      const int ni2 = li + di[q], nj2 = lj + dj[q], nk2 = lk + dk[q];
      if (ni2 < is || ni2 > ie || nj2 < js || nj2 > je || nk2 < ks || nk2 > ke) {
        continue;
      }
      if (hv(36 + q) > tb) {
        tb = hv(36 + q);
        ncd_hot_loc = ((lm*nk + (nk2-ks))*nj + (nj2-js))*ni + (ni2-is);
      }
    }
  }
  static const char *nm[NV] = {"r", "dr", "rdth", "rho", "En", "E", "T", "egn", "res",
    "lres", "tau_r", "tau_th", "cdtkP", "cdtkE", "eg/E", "b1", "b2", "b3", "f", "f1",
    "w", "rf0", "g0", "srcb", "srcr", "de0", "Bk", "Rk", "nfb", "Yk", "Eim", "Eip",
    "Ejm", "Ejp", "Ekm", "Ekp", "Tim", "Tip", "Tjm", "Tjp", "Tkm", "Tkp", "tdia", "trhs",
    "rkP", "rkE", "rkT", "S2"};
  std::ostringstream os;
  os.precision(7);
  os << "<rad_m1> NCD rank=" << global_variable::my_rank << " cycle="
     << pmy_pack->pmesh->ncycle << " step=" << static_cast<int>(impl_nstep)
     << " pass=" << it << " tag=" << tag << " cell=" << lm << "," << lk << "," << lj
     << "," << li;
  for (int q = 0; q < NV; ++q) {os << " " << nm[q] << "=" << hv(q);}
  std::cout << os.str() << std::endl;
  return mloc.loc;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::DetTrace
//! \brief <rad_m1>/implicit_det_trace = N (read only when named): for the first N solves,
//! print an order-independent bitwise hash (exact uint64 sum of position-mixed bit
//! patterns) of every component of iw, u0, opac and the fluid u0/w0, per rank.  Two runs
//! of the same binary are bitwise identical iff these lines are.  Diagnostic only.

namespace {
KOKKOS_INLINE_FUNCTION
uint64_t M1DtMix(const Real x, const uint64_t q) {
  union {Real d; uint64_t u;} c;
  c.u = 0;
  c.d = x;
  uint64_t b = c.u ^ (c.u >> 29);
  return b*(2*q + 0x9E3779B97F4A7C15ULL);
}

uint64_t M1DtHash(const DvceArray5D<Real> &a, const int n, const int k0, const int k1,
                  const int j0, const int j1, const int i0, const int i1) {
  // hash over m and k in [k0,k1), j in [j0,j1), i in [i0,i1)
  const int64_t nk = k1 - k0, nj = j1 - j0, ni = i1 - i0;
  const int64_t tot = static_cast<int64_t>(a.extent(0))*nk*nj*ni;
  uint64_t h = 0;
  if (tot <= 0) {return h;}
  Kokkos::parallel_reduce("m1_dt_hash",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>(DevExeSpace(), 0, tot),
  KOKKOS_LAMBDA(const int64_t q, uint64_t &s) {
    int64_t r = q;
    const int i = static_cast<int>(r % ni) + i0;
    r /= ni;
    const int j = static_cast<int>(r % nj) + j0;
    r /= nj;
    const int k = static_cast<int>(r % nk) + k0;
    const int m = static_cast<int>(r/nk);
    s += M1DtMix(a(m,n,k,j,i), static_cast<uint64_t>(q));
  }, Kokkos::Sum<uint64_t>(h));
  return h;
}
} // namespace

void RadiationM1::DetTrace(const char *tag) {
  if (impl_dtrace <= 0 || impl_nstep >= static_cast<Real>(impl_dtrace)) {return;}
  std::ostringstream os;
  os << "<dtr> r" << global_variable::my_rank << " s" << static_cast<int>(impl_nstep)
     << " " << tag << std::hex;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto put = [&](const char *nm, const DvceArray5D<Real> &a) {
    if (a.size() == 0) {return;}
    // full array (ghosts included) and the active cells alone, per component
    const bool cc = (a.extent_int(2) > indcs.ke) && (a.extent_int(3) > indcs.je) &&
                    (a.extent_int(4) > indcs.ie);
    uint64_t all = 0, act = 0;
    std::ostringstream oc;
    oc << std::hex;
    for (int n = 0; n < a.extent_int(1); ++n) {
      const uint64_t h = M1DtHash(a, n, 0, a.extent_int(2), 0, a.extent_int(3), 0,
                                  a.extent_int(4));
      const uint64_t ha = cc ? M1DtHash(a, n, indcs.ks, indcs.ke + 1, indcs.js,
                                        indcs.je + 1, indcs.is, indcs.ie + 1) : h;
      all += h*(2*static_cast<uint64_t>(n) + 1);
      act += ha*(2*static_cast<uint64_t>(n) + 1);
      oc << " " << (ha & 0xffffffULL);
    }
    os << " " << nm << "=" << all << "/a" << act;
    if (a.extent_int(1) > 1) {os << " [" << oc.str() << " ]";}
  };
  radm1::FluidRef fl = radm1::FluidRef::Get(pmy_pack);
  put("u0", u0);
  put("opac", opac);
  put("fu0", fl.u0);
  put("fw0", fl.w0);
  put("iw", iw);
  put("ifw", ifw);
  put("tau", tau_ten);
  if (fl.uflx != nullptr) {
    put("fx1", fl.uflx->x1f);
    put("fx2", fl.uflx->x2f);
    put("fx3", fl.uflx->x3f);
  }
  std::cout << os.str() << std::endl;
  // M1_DTRACE_DUMP=<tag>: raw dumps of the hydro face fluxes, trial state, FOFC counts
  // and w0 at the first call with that tag (diagnostic)
  static std::string dumped = ",";
  const char *dt = std::getenv("M1_DTRACE_DUMP");
  const std::string tg = "," + std::string(tag) + ",";
  if (dt != nullptr && ("," + std::string(dt) + ",").find(tg) != std::string::npos &&
      dumped.find(tg) == std::string::npos && pmy_pack->phydro != nullptr) {
    dumped += std::string(tag) + ",";
    auto dump = [&](const std::string &nm, auto v) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
      std::string fn = "dtd_" + std::string(tag) + "_" + nm + "_r" +
                       std::to_string(global_variable::my_rank) + ".bin";
      FILE *f = std::fopen(fn.c_str(), "wb");
      if (f != nullptr) {
        std::fwrite(h.data(), sizeof(typename decltype(h)::value_type), h.size(), f);
        std::fclose(f);
      }
    };
    auto ph = pmy_pack->phydro;
    dump("fx2", ph->uflx.x2f);
    dump("fx3", ph->uflx.x3f);
    dump("fofc", ph->fofc);
    dump("utest", ph->utest);
    dump("w0", ph->w0);
  }
}

void RadiationM1::DetTraceScalars(const char *tag, int n, const Real *v) {
  if (impl_dtrace <= 0 || impl_nstep >= static_cast<Real>(impl_dtrace)) {return;}
  if (global_variable::my_rank != 0) {return;}
  char buf[64];
  std::ostringstream os;
  os << "<dts> s" << static_cast<int>(impl_nstep) << " " << tag;
  for (int q = 0; q < n; ++q) {
    std::snprintf(buf, sizeof(buf), " %a", v[q]);
    os << buf;
  }
  std::cout << os.str() << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitWorkRow
//! \brief time2_vstage (m1-sp-order2b): the gas work w^k = rho sum_d dv_d^k (v_d +
//! dv_d^k/2) of the iterate's radiative kick (implicit_vimp's dv^k, this pass) per cell.
//! row = true: TR -= (chat/c) w^k (the E row carries the work, so the stage value
//! satisfies its own equation; after the solve alone it is an O(dt) splitting).
//! row = false (after the write-back, which subtracted the full work): E += (chat/c) w^k
//! and the stored slope with it, so E loses work - w^k there (0 at convergence).
//! Dirichlet (efix) end cells are skipped, as in the write-back.

void RadiationM1::ImplicitWorkRow(bool row) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto u0_ = u0;
  auto uh = FluidRef::Get(pmy_pack).u0;
  const int ivb = iw_vimp + M1_IV_DV;
  const bool trans = trans_on, thrd = trans_x3;
  const bool dbgft = dbg_gas_force_trans;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool elo = (ibc_x1min == M1_IBC_EFIX), ehi = (ibc_x1max == M1_IBC_EFIX);
  auto pos_ = part_pos;
  const int nblkx1 = part_nblk;
  const Real cr = chat/c_light;
  const Real fk = 1.0/dt_sub;
  auto kk_ = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;
  const bool slope = (t2_solve != M1_T2S_NONE);
  // force_reference_work = split: the row also carries the reference work the radiation
  // pays (v dt arad_ref per unit mass; dv above is the residual kick)
  const bool fws = (force_ref == M1_FREF_WB_ARAD) && fref_wsplit;
  auto aref_ = arad_ref;
  const Real dtw = dt_sub;
  if (cs_geom) {
    // STAGE CS3: V and dv are FACE-NORMAL components n_a = s v^a on the cubed sphere, so
    // the kinetic energy is 0.5 rho (n_2^2 + n_3^2 + 2 c n_2 n_3)/s^2 (+ the radial
    // part); a kernel of its own, the one below is the pre-CS3 one
    auto cclw = pmy_pack->pcoord->cos_cell;
    par_for("m1_impl_wk_cs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int ipw = pos_.d_view(m);
      if (!cyclic && ((i == is && ipw == 0 && elo) ||
                      (i == ie && ipw == nblkx1-1 && ehi))) {
        return;
      }
      const Real d1 = iw_(m,ivb,k,j,i);
      Real wk = d1*(iw_(m,M1_IW_V1,k,j,i) + 0.5*d1);
      if (trans && dbgft && thrd) {
        const Real c = cclw(m,k,j), is2 = 1.0/(1.0 - c*c);
        const Real a2 = iw_(m,ivb+1,k,j,i), a3 = iw_(m,ivb+2,k,j,i);
        const Real v2 = iw_(m,M1_IW_V2,k,j,i), v3 = iw_(m,M1_IW_V3,k,j,i);
        wk += is2*(a2*(v2 + 0.5*a2) + a3*(v3 + 0.5*a3)
                   + c*(a2*(v3 + 0.5*a3) + a3*(v2 + 0.5*a2)));
      }
      const Real w = cr*uh(m,IDN,k,j,i)*wk;
      if (row) {
        iw_(m,M1_IW_TR,k,j,i) -= w;
      } else {
        u0_(m,M1_E,k,j,i) += w;
        if (slope) {
          kk_(m,M1_T2_E,k,j,i) += w*fk;
        }
      }
    });
  } else {
  par_for("m1_impl_wk", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int ipw = pos_.d_view(m);
    if (!cyclic && ((i == is && ipw == 0 && elo) ||
                    (i == ie && ipw == nblkx1-1 && ehi))) {
      return;
    }
    Real dv = iw_(m,ivb,k,j,i);
    Real wk = dv*(iw_(m,M1_IW_V1,k,j,i) + 0.5*dv);
    if (trans && dbgft) {
      dv = iw_(m,ivb+1,k,j,i);
      wk += dv*(iw_(m,M1_IW_V2,k,j,i) + 0.5*dv);
      if (thrd) {
        dv = iw_(m,ivb+2,k,j,i);
        wk += dv*(iw_(m,M1_IW_V3,k,j,i) + 0.5*dv);
      }
    }
    const Real w = cr*uh(m,IDN,k,j,i)*wk;
    if (row) {
      iw_(m,M1_IW_TR,k,j,i) -= w;
    } else {
      u0_(m,M1_E,k,j,i) += w;
      if (slope) {
        kk_(m,M1_T2_E,k,j,i) += w*fk;
      }
    }
  });
  }   // cs_geom
  if (fws) {
    // force_reference_work = split: the reference work (v dt rho arad_ref), which the
    // radiation pays, rides the row too (a separate kernel: the default one is untouched)
    par_for("m1_impl_wkref", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int ipw = pos_.d_view(m);
      if (!cyclic && ((i == is && ipw == 0 && elo) ||
                      (i == ie && ipw == nblkx1-1 && ehi))) {
        return;
      }
      // at the iterate's velocity (v_old + the residual kick dv^k): the stage value's
      // own velocity, which the implicit stage needs (the old one lags by O(dt))
      const Real vk = iw_(m,M1_IW_V1,k,j,i) + iw_(m,ivb,k,j,i);
      const Real w = cr*uh(m,IDN,k,j,i)*vk*dtw*aref_(m,k,j,i);
      if (row) {
        iw_(m,M1_IW_TR,k,j,i) -= w;
      } else {
        u0_(m,M1_E,k,j,i) += w;
        if (slope) {
          kk_(m,M1_T2_E,k,j,i) += w*fk;
        }
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitVimpBuild
//! \brief implicit_vimp (rad_m1_implicit.hpp, tests_m1/runs_3v_vimplicit): the Newton
//! form of the enthalpy flux with the gas velocity implicit, for this Picard pass.
//!
//! (1) per cell and axis d: the Jacobian row of the velocity change the write-back will
//!     apply, dv_d(E') = P_m E'_{-1} + P_0 E'_0 + P_p E'_{+1} (+ const), built from the
//!     face-normal eliminated face fluxes exactly as step (g) forms them (theta, D_dd,
//!     the AP-HLL blend, Marshak faces; the od and g0 terms are left out of the
//!     Jacobian),
//!     and dv_d^k, the velocity change of the ITERATE's face fluxes (write-back rule).
//!     Both are exchanged (M1_NVIMP_X components).
//! (2) per cell: the operator coefficients of
//!       sum_faces sigma (dt/dx_d)(chat/c) E_f^k [(1 + D_dd)_f (P dv_d)(E')_f]
//!     (P dv)_f the mean of the two cells' rows, cells up to two apart), and the
//!     right-hand side
//!       -sigma (dt/dx_d)(chat/c) E_f^k [(1 + D_dd) (dv_d^k - (P dv_d)(E^k))
//!     + sum_{e != d} D_de dv_e^k]_f.  E_f^k is the face E the implicit_enthalpy mode
//!     carries (M1EnthEf), so the increment rides the same face value as the plm/central
//!     correction; the velocity increment is the face MEAN of the two cells.

void RadiationM1::ImplicitVimpBuild() {
  if (!fl_on || !coupling || !gas_feedback || !dbg_gas_force) {
    ImplFatal("<rad_m1>/implicit_vimp needs hydro with coupling, gas_feedback and "
              "dbg_gas_force on (the implicit velocity IS the solve's radiative force)");
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ngh = indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto ifw_ = ifw;
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  auto f0_ = f0x1;
  auto f2_ = f0x2;
  auto f3_ = f0x3;
  auto th2_ = thx2;
  auto th3_ = thx3;
  const bool lm = (impl_tlim != M1_TLIM_NONE) || blat_on;
  const bool blt = blat_on;   // vimphr-1009: the x2/x3 berthon / half-range part
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = (cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs).d_view;   // CS1: seams open
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool thrd = trans_x3;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real mq = marshak_q;
  // vet_col_surface_q (rad_m1_vetcol.cpp): the OUTER x1 Marshak q of each column from
  // its formal solution; the branches below shadow mq with it (off: mq itself)
  const bool vqs = (vet_col && vcol_sq) || hr_q;   // hrup-1009: + half-range
  auto vq_ = vcol_q;
  const Real mqo = marshak_q;
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);
  const bool bmhalf = impl_bmom_half;
  const bool fref = (force_ref == M1_FREF_WB_ARAD);
  auto aref_ = arad_ref;
  const bool ftr = dbg_gas_force_trans;
  const int enm = impl_enth;
  const int b = iw_vimp;
  auto uh = FluidRef::Get(pmy_pack).u0;
  const Real jsc = impl_vimp_jscale;
  // hesdirk2, time2_enth_vel = start: the plm face a is built from the stage-start a
  // (a - DA) and DA enters as its face mean (M1EnthEfT)
  const bool t2vs = (t2_afmode != 0) && (t2_solve == M1_T2S_STAGE1 ||
                                          t2_solve == M1_T2S_STAGE2);
  const bool t2afc = (t2_afmode == 2);
  const int bda = b + M1_IV_DA;
  // m1-sph2 (tests_m1/runs_5h_sph2): on the spherical-polar wedge the face-flux rows
  // take the centre distance dxface (and, for a chi(f) closure, the S2 radial
  // integrating factor M1SphDrr), and the enthalpy divergence takes dt A_f/V_i per face,
  // exactly as the rows of ImplicitSolve.  The Cartesian expressions are untouched; the
  // sp forms are `if (sph)` overwrites.
  const bool sph = sph_geom;
  const bool sphq = sph_q;
  auto cvol = pmy_pack->pcoord->volume;
  auto carea = pmy_pack->pcoord->area;
  auto cdxf = pmy_pack->pcoord->dxface;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const bool fwd = sph && impl_face_wdist;   // implicit_face_weight = distance
  // STAGE CS3 (C9): implicit_vimp on the cubed sphere.  The work array's transverse
  // velocities and face fluxes are FACE-NORMAL components, in which the velocity change
  // of a face-normal radiative kick is exactly (normal force)/rho, as on a Cartesian
  // mesh (the metric cancels: dm_a covariant = (a + c b)/s gives s v^a = a/rho).  So
  // only the geometry changes: the two-point distance of the transverse face rows (centre
  // arc times sin, the mirror-pair distance at a panel seam), the canonical seam areas,
  // and the Jacobian row of a SEAM ghost (its P rows cross the seam as scalars, which a
  // swapped or reversed axis scrambles): there the cell's own row, reflected, stands in
  // (exact for a mirror-symmetric state; P only sets the Newton direction, the
  // right-hand side carries the true dv^k, so the converged answer does not depend on it)
  // (the effective two-point distances and face areas, cs_dxf_eff / cs_area_eff, are
  // built once in CubedS1Init; the sp branches below read them in place of Coordinates,
  // so the sp arithmetic is textually the one of the wedge)
  if (cs_geom) {
    carea = cs_area_eff;
    cdxf = cs_dxf_eff;
  }

  // (1) the Jacobian rows and dv^k
  par_for("m1_vimp_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dd = fmax(uh(m,IDN,k,j,i), 1.0e-300);
    const Real sc = dt/(cl*dd);
    const Real sj = sc*jsc;
    const Real ktc = iw_(m,M1_IW_KT,k,j,i);
    const int ipos = pos_(m);
    const bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    Real cL[2], cR[2], kf[2], wf[2];
    // ---- x1
    for (int s = 0; s < 2; ++s) {
      const int fi = i + s;
      const bool lo = (fi == is) && botb && !cyclic;
      const bool hi = (fi == ie+1) && topb && !cyclic;
      cL[s] = 0.0;
      cR[s] = 0.0;
      if (lo || hi) {
        kf[s] = ktc;
        wf[s] = bmhalf ? 0.5 : 1.0;
        if ((lo ? bclo : bchi) == M1_IBC_MARSHAK) {
          if (lo) {
            cR[s] = -cl*mq;
          } else {
            const Real mq = vqs ? vq_(m,k,j) : mqo;   // vet_col_surface_q
            cL[s] = cl*mq;
          }
        }
      } else {
        const int im = (cyclic && fi == is) ? ie : (fi-1);
        const int ip = (cyclic && fi == ie+1) ? is : fi;
        kf[s] = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,im), iw_(m,M1_IW_KT,k,j,ip), cx1f, m, im,
                             ip, fwd);
        wf[s] = 0.5;
        const Real th = 1.0/(1.0 + ch*dt*kf[s]);
        const Real bs = th*ch*cl*dt/mbsize(m).dx1;
        cL[s] = bs*M1DDiag(iw_,vd_,dfull,m,0,k,j,im);
        cR[s] = -bs*M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);
        if (aphll) {
          const Real al = ifw_(m,M1_IFW_AL,k,j,fi);
          cL[s] = (1.0 - al)*cL[s] + (cl/ch)*ifw_(m,M1_IFW_HCL,k,j,fi);
          cR[s] = (1.0 - al)*cR[s] + (cl/ch)*ifw_(m,M1_IFW_HCR,k,j,fi);
        }
        if (sph) {
          // sp: the face distance, and the S2 integrating factor
          const Real bsp = th*ch*cl*dt/cdxf.x1f(m,k,j,fi);
          Real wl = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);
          Real wr = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);
          if (sphq) {
            const Real rf = cx1f(m,fi);
            wl = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,im), iw_(m,M1_IW_N1,k,j,im),
                          SQR(cx1v(m,im)/rf));
            wr = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,ip), iw_(m,M1_IW_N1,k,j,ip),
                          SQR(cx1v(m,ip)/rf));
          }
          cL[s] = bsp*wl;
          cR[s] = -bsp*wr;
          // sp-blend-1008: the blend of the face flux the row applies (berthon | blend)
          if (aphll) {
            const Real al = ifw_(m,M1_IFW_AL,k,j,fi);
            cL[s] = (1.0 - al)*cL[s] + (cl/ch)*ifw_(m,M1_IFW_HCL,k,j,fi);
            cR[s] = (1.0 - al)*cR[s] + (cl/ch)*ifw_(m,M1_IFW_HCR,k,j,fi);
          }
        }
      }
    }
    iw_(m,b+0,k,j,i) = sj*wf[0]*kf[0]*cL[0];
    iw_(m,b+1,k,j,i) = sj*(wf[0]*kf[0]*cR[0] + wf[1]*kf[1]*cL[1]);
    iw_(m,b+2,k,j,i) = sj*wf[1]*kf[1]*cR[1];
    Real dv = sc*(wf[0]*kf[0]*f0_(m,k,j,i) + wf[1]*kf[1]*f0_(m,k,j,i+1));
    if (fref) {dv -= dt*aref_(m,k,j,i);}
    iw_(m,b+M1_IV_DV,k,j,i) = dv;
    // ---- x2 and x3: a physical face carries F = 0
    for (int d = 1; d < 3; ++d) {
      const bool on = ftr && ((d == 1) || thrd);
      Real p0 = 0.0, p1 = 0.0, p2 = 0.0, dvd = 0.0;
      if (on) {
        BoundaryFlag qlo = mbbcs(m,(d == 1) ? BoundaryFace::inner_x2
                                             : BoundaryFace::inner_x3);
        BoundaryFlag qhi = mbbcs(m,(d == 1) ? BoundaryFace::outer_x2
                                             : BoundaryFace::outer_x3);
        const bool plo = (qlo != BoundaryFlag::block) && (qlo != BoundaryFlag::periodic);
        const bool phi = (qhi != BoundaryFlag::block) && (qhi != BoundaryFlag::periodic);
        const int c = (d == 1) ? j : k;
        const int cs = (d == 1) ? js : ks, ce = (d == 1) ? je : ke;
        const Real dxd = (d == 1) ? mbsize(m).dx2 : mbsize(m).dx3;
        Real fl[2];
        for (int s = 0; s < 2; ++s) {
          const bool phys = (s == 0) ? (c == cs && plo) : (c == ce && phi);
          const int cf = c + s;                     // the face index along d
          const int kq = (d == 2) ? cf : k, jq = (d == 1) ? cf : j;
          const int kl = (d == 2) ? (cf-1) : k, jl = (d == 1) ? (cf-1) : j;
          fl[s] = (d == 1) ? f2_(m,k,cf,i) : f3_(m,cf,j,i);
          cL[s] = 0.0;
          cR[s] = 0.0;
          if (phys) {
            kf[s] = ktc;
            wf[s] = bmhalf ? 0.5 : 1.0;
          } else {
            kf[s] = 0.5*(iw_(m,M1_IW_KT,kl,jl,i) + iw_(m,M1_IW_KT,kq,jq,i));
            wf[s] = 0.5;
            Real th = 1.0/(1.0 + ch*dt*kf[s]);
            if (lm) {th = (d == 1) ? th2_(m,k,cf,i) : th3_(m,cf,j,i);}
            const Real bs = th*ch*cl*dt/dxd;
            cL[s] = bs*M1DDiag(iw_,vd_,dfull,m,d,kl,jl,i);
            cR[s] = -bs*M1DDiag(iw_,vd_,dfull,m,d,kq,jq,i);
            if (sph) {
              const Real dxs = (d == 1) ? cdxf.x2f(m,k,cf,i) : cdxf.x3f(m,cf,j,i);
              const Real bsp = th*ch*cl*dt/dxs;
              cL[s] = bsp*M1DDiag(iw_,vd_,dfull,m,d,kl,jl,i);
              cR[s] = -bsp*M1DDiag(iw_,vd_,dfull,m,d,kq,jq,i);
            }
            // vimphr-1009: implicit_flux_faces = all.  The face flux the write-back
            // stores (and dv^k is built from) is th (1 - AL) F_central + (c/chat)
            // (HCL E_L + HCR E_R), th already carrying (1 - AL) (thx2/thx3, lm true);
            // the Jacobian row takes the same berthon / half-range part
            if (blt) {
              const Real hcl = (d == 1) ? bw2_(m,M1_IFW_HCL,k,cf,i)
                                        : bw3_(m,M1_IFW_HCL,cf,j,i);
              const Real hcr = (d == 1) ? bw2_(m,M1_IFW_HCR,k,cf,i)
                                        : bw3_(m,M1_IFW_HCR,cf,j,i);
              cL[s] += (cl/ch)*hcl;
              cR[s] += (cl/ch)*hcr;
            }
          }
        }
        p0 = sj*wf[0]*kf[0]*cL[0];
        p1 = sj*(wf[0]*kf[0]*cR[0] + wf[1]*kf[1]*cL[1]);
        p2 = sj*wf[1]*kf[1]*cR[1];
        dvd = sc*(wf[0]*kf[0]*fl[0] + wf[1]*kf[1]*fl[1]);
      }
      iw_(m,b+3*d,k,j,i) = p0;
      iw_(m,b+3*d+1,k,j,i) = p1;
      iw_(m,b+3*d+2,k,j,i) = p2;
      iw_(m,b+M1_IV_DV+d,k,j,i) = dvd;
    }
  });
  ImplicitHaloExchange(M1_NVIMP_X, b);
  if (cs_geom) {
    // STAGE CS3: the Jacobian rows of a SEAM ghost (first layer) for the seam-normal
    // direction: the adjacent cell's own row, reflected (see above)
    auto cseam = cs_seam.d_view;
    par_for("m1_vimp_csrow", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (j == js && cseam(m,0) != 0) {
        for (int r = 0; r < 3; ++r) {iw_(m,b+3+r,k,j-1,i) = -iw_(m,b+3+2-r,k,j,i);}
      }
      if (j == je && cseam(m,1) != 0) {
        for (int r = 0; r < 3; ++r) {iw_(m,b+3+r,k,j+1,i) = -iw_(m,b+3+2-r,k,j,i);}
      }
      if (k == ks && cseam(m,2) != 0) {
        for (int r = 0; r < 3; ++r) {iw_(m,b+6+r,k-1,j,i) = -iw_(m,b+6+2-r,k,j,i);}
      }
      if (k == ke && cseam(m,3) != 0) {
        for (int r = 0; r < 3; ++r) {iw_(m,b+6+r,k+1,j,i) = -iw_(m,b+6+2-r,k,j,i);}
      }
    });
  }

  // (2) the operator coefficients and the right-hand side
  par_for_lb("m1_vimp_j", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int ipos = pos_(m);
    const bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    for (int n = M1_IV_X1M2; n <= M1_IV_JRHS; ++n) {iw_(m,b+n,k,j,i) = 0.0;}
    if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                    (i == ie && topb && bchi == M1_IBC_EFIX))) {
      return;   // the Dirichlet row is replaced whole
    }
    const Real cr = ch/cl;
    BoundaryFlag q[6];
    q[0] = mbbcs(m,BoundaryFace::inner_x1);
    q[1] = mbbcs(m,BoundaryFace::outer_x1);
    q[2] = mbbcs(m,BoundaryFace::inner_x2);
    q[3] = mbbcs(m,BoundaryFace::outer_x2);
    q[4] = mbbcs(m,BoundaryFace::inner_x3);
    q[5] = mbbcs(m,BoundaryFace::outer_x3);
    Real jd = 0.0, jr = 0.0;
    const int nd = thrd ? 3 : 2;
    for (int d = 0; d < nd; ++d) {
      const bool plo = (q[2*d] != BoundaryFlag::block) &&
                       (q[2*d] != BoundaryFlag::periodic);
      const bool phi = (q[2*d+1] != BoundaryFlag::block) &&
                       (q[2*d+1] != BoundaryFlag::periodic);
      const int c = (d == 0) ? i : ((d == 1) ? j : k);
      const int cs = (d == 0) ? is : ((d == 1) ? js : ks);
      const int ce = (d == 0) ? ie : ((d == 1) ? je : ke);
      const bool cyc = (d == 0) && cyclic;
      // physical faces carry no enthalpy flux
      const bool flo = (d == 0) ? (i == is && botb && !cyclic) : (c == cs && plo);
      const bool fhi = (d == 0) ? (i == ie && topb && !cyclic) : (c == ce && phi);
      const int hl = plo ? 0 : ngh;
      const int hh = phi ? 0 : ngh;
      const Real dxd = (d == 0) ? mbsize(m).dx1 : ((d == 1) ? mbsize(m).dx2
                                                            : mbsize(m).dx3);
      const Real nu = dt/dxd;
      const int ac = (d == 0) ? M1_IW_ADV : ((d == 1) ? M1_IW_A2 : M1_IW_A3);
      const int vc = (d == 0) ? M1_IW_V1 : ((d == 1) ? M1_IW_V2 : M1_IW_V3);
      // the cell at offset o along d (x1 wraps inside the block when cyclic)
      Real E[5], A[5], D[5], V[3], P[3][3], DV[3][3], DO[3][3], DA[5];
      bool av[5];
      for (int o = -2; o <= 2; ++o) {
        int kk = k, jj = j, ii = i;
        bool okk;
        int cc = M1EnthIdx(c + o, cs, ce, cyc, hl, hh, okk);
        if (!okk) {cc = c + o;}   // a physical ghost: read, multiplied by zero
        if (d == 0) {ii = cc;} else if (d == 1) {jj = cc;} else {kk = cc;}
        av[o+2] = okk;
        E[o+2] = iw_(m,M1_IW_EP,kk,jj,ii);
        A[o+2] = iw_(m,ac,kk,jj,ii);
        DA[o+2] = t2vs ? iw_(m,bda+d,kk,jj,ii) : 0.0;
        D[o+2] = M1DDiag(iw_,vd_,dfull,m,d,kk,jj,ii);
        if (o >= -1 && o <= 1) {
          V[o+1] = iw_(m,vc,kk,jj,ii);
          for (int r = 0; r < 3; ++r) {P[o+1][r] = iw_(m,b+3*d+r,kk,jj,ii);}
          for (int e = 0; e < 3; ++e) {
            DV[o+1][e] = iw_(m,b+M1_IV_DV+e,kk,jj,ii);
            DO[o+1][e] = (e == d) ? 0.0 : M1DOffC(iw_,vd_,dfull,m,(d < e) ? d : e,
                                                  (d < e) ? e : d,kk,jj,ii);
          }
        }
      }
      Real phf[2] = {0.0, 0.0};
      const Real jr0 = jr;
      for (int s = 0; s < 2; ++s) {
        if ((s == 0) ? flo : fhi) continue;
        const Real sg = (s == 0) ? -1.0 : 1.0;
        const int ol = s - 1;        // offset of the face's L cell
        const int xl = ol + 2;       // its index into E/A/D
        const bool ok = av[xl-1] && av[xl+2];
        const Real vf = 0.5*(V[ol+1] + V[ol+2]);
        const Real ef = t2vs
            ? M1EnthEfT(enm, E[xl-1], E[xl], E[xl+1], E[xl+2], ok, A[xl-1], A[xl],
                        A[xl+1], A[xl+2], DA[xl-1], DA[xl], DA[xl+1], DA[xl+2], vf, t2afc)
            : M1EnthEf(enm, E[xl-1], E[xl], E[xl+1], E[xl+2], ok, A[xl-1],
                       A[xl], A[xl+1], A[xl+2], vf);
        const Real ph = nu*cr*ef*(1.0 + 0.5*(D[xl] + D[xl+1]));
        phf[s] = ph;
        const Real dvk = 0.5*(DV[ol+1][d] + DV[ol+2][d]);
        Real ps = 0.0;
        for (int e = 0; e < 3; ++e) {
          if (e == d) continue;
          ps += 0.5*(DO[ol+1][e] + DO[ol+2][e])*0.5*(DV[ol+1][e] + DV[ol+2][e]);
        }
        ps *= nu*cr*ef;
        // (P dv)(E^k) at the face: mean of the two cells' rows
        Real pe = 0.0;
        for (int h = 0; h < 2; ++h) {
          const int x = xl + h;      // cell index into E
          const int r = ol + 1 + h;  // its index into P
          pe += 0.5*(P[r][0]*E[x-1] + P[r][1]*E[x] + P[r][2]*E[x+1]);
        }
        jr -= sg*(ph*dvk + ps - ph*pe);
      }
      if (sph) {
        // sp: the same two faces with dt A_f/V_i each instead of dt/dx_d (the Cartesian
        // loop above is left as it is and its contribution taken out again)
        jr = jr0;
        const Real iv = dt/cvol(m,k,j,i);
        for (int s = 0; s < 2; ++s) {
          phf[s] = 0.0;
          if ((s == 0) ? flo : fhi) continue;
          Real af;
          if (d == 0) {
            af = carea.x1f(m,k,j,i+s);
          } else if (d == 1) {
            af = carea.x2f(m,k,j+s,i);
          } else {
            af = carea.x3f(m,k+s,j,i);
          }
          const Real nus = af*iv;
          const Real sg = (s == 0) ? -1.0 : 1.0;
          const int ol = s - 1;
          const int xl = ol + 2;
          const bool ok = av[xl-1] && av[xl+2];
          const Real vf = 0.5*(V[ol+1] + V[ol+2]);
          const Real ef = t2vs
              ? M1EnthEfT(enm, E[xl-1], E[xl], E[xl+1], E[xl+2], ok, A[xl-1], A[xl],
                          A[xl+1], A[xl+2], DA[xl-1], DA[xl], DA[xl+1], DA[xl+2], vf,
                          t2afc)
              : M1EnthEf(enm, E[xl-1], E[xl], E[xl+1], E[xl+2], ok, A[xl-1],
                         A[xl], A[xl+1], A[xl+2], vf);
          const Real ph = nus*cr*ef*(1.0 + 0.5*(D[xl] + D[xl+1]));
          phf[s] = ph;
          const Real dvk = 0.5*(DV[ol+1][d] + DV[ol+2][d]);
          Real ps = 0.0;
          for (int e = 0; e < 3; ++e) {
            if (e == d) continue;
            ps += 0.5*(DO[ol+1][e] + DO[ol+2][e])*0.5*(DV[ol+1][e] + DV[ol+2][e]);
          }
          ps *= nus*cr*ef;
          Real pe = 0.0;
          for (int h = 0; h < 2; ++h) {
            const int x = xl + h;
            const int r = ol + 1 + h;
            pe += 0.5*(P[r][0]*E[x-1] + P[r][1]*E[x] + P[r][2]*E[x+1]);
          }
          jr -= sg*(ph*dvk + ps - ph*pe);
        }
      }
      const Real pl = 0.5*phf[0], ph = 0.5*phf[1];
      const Real cm2 = -pl*P[0][0];
      const Real cm1 = ph*P[1][0] - pl*(P[0][1] + P[1][0]);
      const Real c00 = ph*(P[1][1] + P[2][0]) - pl*(P[0][2] + P[1][1]);
      const Real cp1 = ph*(P[1][2] + P[2][1]) - pl*P[1][2];
      const Real cp2 = ph*P[2][2];
      jd += c00;
      if (d == 0) {
        iw_(m,b+M1_IV_X1M2,k,j,i) = cm2;
        iw_(m,b+M1_IV_J1M,k,j,i) = cm1;
        iw_(m,b+M1_IV_J1P,k,j,i) = cp1;
        iw_(m,b+M1_IV_X1P2,k,j,i) = cp2;
      } else {
        const int o0 = (d == 1) ? M1_IV_X2M2 : M1_IV_X3M2;
        iw_(m,b+o0,k,j,i) = cm2;
        iw_(m,b+o0+1,k,j,i) = cm1;
        iw_(m,b+o0+2,k,j,i) = cp1;
        iw_(m,b+o0+3,k,j,i) = cp2;
      }
    }
    iw_(m,b+M1_IV_JD,k,j,i) = jd;
    iw_(m,b+M1_IV_JRHS,k,j,i) = jr;
  });
}

//----------------------------------------------------------------------------------------
//! \fn M1ImplSrcLaunch
//! \brief step (c) of RadiationM1::ImplicitSolve: the linearised emission/absorption
//! source kernel (m1_impl_src), instantiated on the ideal-gas tag.
template <typename Ctx, typename Idl>
void M1ImplSrcLaunch(const Ctx &ctx_, Idl) {
  auto ar = std::get<0>(ctx_);
  auto ch = std::get<1>(ctx_);
  auto cl = std::get<2>(ctx_);
  auto dt = std::get<3>(ctx_);
  auto ec_ = std::get<4>(ctx_);
  auto ecnt = std::get<5>(ctx_);
  auto eos = std::get<6>(ctx_);
  auto gam = std::get<7>(ctx_);
  auto gasx = std::get<8>(ctx_);
  auto gnewt = std::get<9>(ctx_);
  auto igb = std::get<10>(ctx_);
  auto igf = std::get<11>(ctx_);
  auto igm = std::get<12>(ctx_);
  auto igr = std::get<13>(ctx_);
  auto igy = std::get<14>(ctx_);
  auto iw_ = std::get<15>(ctx_);
  auto opac_ = std::get<16>(ctx_);
  auto uh = std::get<17>(ctx_);
  auto usec = std::get<18>(ctx_);
  auto nmb1 = std::get<19>(ctx_);
  auto ks = std::get<20>(ctx_);
  auto ke = std::get<21>(ctx_);
  auto js = std::get<22>(ctx_);
  auto je = std::get<23>(ctx_);
  auto is = std::get<24>(ctx_);
  auto ie = std::get<25>(ctx_);
  auto sstab = std::get<26>(ctx_);
  auto srb = [=] KOKKOS_FUNCTION (Idl idl, const int m, const int k, const int j,
                                  const int i) {
#if defined(KOKKOS_ENABLE_CUDA)
    // nvcc: an extended lambda may not capture a variable for the first time inside
    // an `if constexpr` branch, so name every capture up front (no code is
    // generated; other backends capture exactly what the instantiation reads).
    (void)ar; (void)ch; (void)cl; (void)dt; (void)ec_; (void)ecnt; (void)eos; (void)gam;
    (void)gasx; (void)gnewt; (void)igb; (void)igf; (void)igm; (void)igr; (void)igy;
    (void)iw_; (void)opac_; (void)uh; (void)usec; (void)sstab;
#endif
    Real rkpv = opac_(m,M1_OP_P,k,j,i);
    Real rkev = opac_(m,M1_OP_E,k,j,i);
    if (rkpv == 0.0 && rkev == 0.0) {
      iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
      iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
      return;
    }
    Real dd = uh(m,IDN,k,j,i);
    Real tk = iw_(m,M1_IW_TP,k,j,i);
    Real nmiss = 0.0;
    auto thc = [&]() {
      if constexpr (decltype(idl)::value) {
        return M1EosIdeal{gam};   // never called: usec is false here
      } else {
        return M1EosCached<decltype(eos), decltype(ec_)>{eos, ec_, m, k, j, i, ecnt,
                                                         &nmiss};
      }
    }();
    auto thd = [&]() {
      if constexpr (decltype(idl)::value) {
        return M1EosIdeal{gam};
      } else {
        return M1EosDirect<decltype(eos)>{eos};
      }
    }();
    Real ee, cv;
    if (usec) {
      thc(dd, tk, ee, cv);
    } else {
      thd(dd, tk, ee, cv);
    }
    Real t3 = tk*tk*tk;
    Real t4 = t3*tk;
    Real de0 = iw_(m,M1_IW_DE0,k,j,i);
    Real bk = dd*cv + 4.0*cl*dt*rkpv*ar*t3;
    Real rk = iw_(m,M1_IW_EGN,k,j,i) - ee - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
    // MILESTONE 3g, the Newton SAFEGUARD.  The temperature this pass starts from was
    // produced by the Newton step of the previous pass, whose linearisation dropped
    // the curvature of e(T) and of T^4.  Here -- where e(T_k) has just been
    // evaluated anyway, so the test is FREE -- the exact nonlinear gas residual
    //   y(T) = rho e(T) + c dt rho kappa_P a T^4 - rho e^n - c dt rho kappa_E E0'
    //        = -(R_k + c dt rho kappa_E E')
    // is compared with the value it had BEFORE that step, at the SAME E' (the solve
    // has not moved E since).  If it did not decrease, the Newton step is discarded
    // and the bracketed root find of the pre-3g scheme is run for this cell.
    if (gnewt) {
      Real yprev = iw_(m,igy,k,j,i);
      if (yprev > 0.0) {
        Real ep = iw_(m,M1_IW_EP,k,j,i);
        Real ynow = fabs(rk + cl*dt*rkev*ep);
        // ...but only where the residual still MEANS something.  Once the cell has
        // converged, y is a difference of numbers that cancel to round-off and it
        // stops decreasing monotonically; without this floor every converged cell
        // buys a bracketed root find in every remaining pass (measured: 13 % of all
        // cell-passes fell back, against 0.6 % with it).
        Real ysc = fmax(fabs(iw_(m,M1_IW_EGN,k,j,i)), cl*dt*rkev*fmax(ep, 0.0));
        if (!(ynow < yprev) && ynow > M1_IMPL_TRTOL*ysc) {
          Real tn = tk;
          bool ok = true;
          if (usec) {
            (void) M1ImplTemperatureT(thc, dd, tk, iw_(m,M1_IW_EGN,k,j,i),
                                      cl*dt*rkpv*ar, cl*dt*rkev*(ep + de0), tn, ok);
          } else {
            (void) M1ImplTemperatureT(thd, dd, tk, iw_(m,M1_IW_EGN,k,j,i),
                                      cl*dt*rkpv*ar, cl*dt*rkev*(ep + de0), tn, ok);
          }
          if (ok && tn > 0.0) {
            tk = tn;
            iw_(m,M1_IW_TP,k,j,i) = tk;
            if (usec) {thc(dd, tk, ee, cv);} else {thd(dd, tk, ee, cv);}
            t3 = tk*tk*tk;
            t4 = t3*tk;
            bk = dd*cv + 4.0*cl*dt*rkpv*ar*t3;
            rk = iw_(m,M1_IW_EGN,k,j,i) - ee - cl*dt*rkpv*ar*t4
                 + cl*dt*rkev*de0;
            // fixbundle-1009 F4: count the fallback only when it is APPLIED (output only)
            iw_(m,igf,k,j,i) += 1.0;
          }
        }
      }
      iw_(m,igb,k,j,i) = bk;
      iw_(m,igr,k,j,i) = rk;
    }
    if (gasx) {
      iw_(m,igm,k,j,i) += nmiss;
    }
    Real emis = dt*ch*rkpv*ar;
    if (sstab && bk > 0.0) {
      // implicit_src_stable (stall_1002): the SAME Schur-eliminated row in a form
      // without cancellation.  With K = c dt rho kappa (2e8 at the BSG wall) the plain
      // form takes SRCB = c dt k_E - c dt k_E (bk - rho c_v)/bk, a difference of two
      // numbers ~K that leaves O(1), and SRCR likewise subtracts ~K E terms; both carry
      // an absolute round-off ~eps K that moved the solved E by ~1e-7 relative in every
      // Picard pass (a floor above implicit_tol 1e-8).  Exactly (bk - rho c_v =
      // 4 c dt k_P a T^3, emis 4 T^3 = (ch/cl)(bk - rho c_v)):
      //   SRCB = ch dt k_E w,
      //   SRCR = w (emis T^4 - ch dt k_E de0) + (ch/cl) q (e^n - e_k)
      // with w = rho c_v/bk and q = 1 - w = 4 cl dt k_P a T^3/bk.
      const Real w = dd*cv/bk;
      const Real q = 4.0*cl*dt*rkpv*ar*t3/bk;
      iw_(m,M1_IW_SRCB,k,j,i) = dt*ch*rkev*w;
      iw_(m,M1_IW_SRCR,k,j,i) = w*(emis*t4 - dt*ch*rkev*de0)
                                + (ch/cl)*q*(iw_(m,M1_IW_EGN,k,j,i) - ee);
      return;
    }
    Real kk = (bk > 0.0) ? (emis*4.0*t3*cl*dt*rkev/bk) : 0.0;
    iw_(m,M1_IW_SRCB,k,j,i) = dt*ch*rkev - kk;
    iw_(m,M1_IW_SRCR,k,j,i) = emis*t4 - dt*ch*rkev*de0
                              + ((bk > 0.0) ? (emis*4.0*t3*rk/bk) : 0.0);
  };
  par_for("m1_impl_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    srb(Idl{}, m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \struct M1OpacFn
//! \brief implicit_tsolve_opac: rho kappa_P and rho kappa_E at (d, T), code units, by the
//! same calls as step (a) of the Picard loop (table or analytic).

struct M1OpacFn {
  int otype;
  M1OpacTab ot;
  Real kp, kev, kf, ks, rref, tref, aa, bb;
  KOKKOS_INLINE_FUNCTION
  void operator()(const Real d, const Real t, Real &rkp, Real &rke) const {
    Real op, oe, of, os;
    if (otype == M1_OPAC_TABLE) {
      M1TableOpacities(ot, d, t, op, oe, of, os);
    } else {
      M1Opacities(otype, d, t, kp, kev, kf, ks, rref, tref, aa, bb, op, oe, of, os);
    }
    rkp = d*op;
    rke = d*oe;
  }
};

//----------------------------------------------------------------------------------------
//! \fn M1ImplTsolveLaunch
//! \brief step (f) of RadiationM1::ImplicitSolve: accept E', solve for T'
//! (m1_impl_tsolve), instantiated on the ideal-gas tag.
template <typename Ctx, typename Idl>
void M1ImplTsolveLaunch(const Ctx &ctx_, Idl) {
  auto ar = std::get<0>(ctx_);
  auto cl = std::get<1>(ctx_);
  auto dt = std::get<2>(ctx_);
  auto ec_ = std::get<3>(ctx_);
  auto ecnt = std::get<4>(ctx_);
  auto efl = std::get<5>(ctx_);
  auto eos = std::get<6>(ctx_);
  auto escale = std::get<7>(ctx_);
  auto gasx = std::get<8>(ctx_);
  auto gnewt = std::get<9>(ctx_);
  auto igb = std::get<10>(ctx_);
  auto igf = std::get<11>(ctx_);
  auto igm = std::get<12>(ctx_);
  auto igr = std::get<13>(ctx_);
  auto igy = std::get<14>(ctx_);
  auto iw_ = std::get<15>(ctx_);
  auto opac_ = std::get<16>(ctx_);
  auto plog = std::get<17>(ctx_);
  auto uh = std::get<18>(ctx_);
  auto usec = std::get<19>(ctx_);
  auto nmb1 = std::get<20>(ctx_);
  auto ks = std::get<21>(ctx_);
  auto ke = std::get<22>(ctx_);
  auto js = std::get<23>(ctx_);
  auto je = std::get<24>(ctx_);
  auto is = std::get<25>(ctx_);
  auto ie = std::get<26>(ctx_);
  auto tso = std::get<27>(ctx_);
  auto opf = std::get<28>(ctx_);
  auto tsm = std::get<29>(ctx_);
  auto tsmin = std::get<30>(ctx_);
  auto tf_ = std::get<31>(ctx_);
  auto tfz = std::get<32>(ctx_);
  auto tsb = [=] KOKKOS_FUNCTION (Idl idl, const int m, const int k, const int j,
                                  const int i) M1_INL {
#if defined(KOKKOS_ENABLE_CUDA)
    // nvcc: an extended lambda may not capture a variable for the first time inside
    // an `if constexpr` branch, so name every capture up front (no code is
    // generated; other backends capture exactly what the instantiation reads).
    (void)ar; (void)cl; (void)dt; (void)ec_; (void)ecnt; (void)efl; (void)eos;
    (void)escale; (void)gasx; (void)gnewt; (void)igb; (void)igf; (void)igm; (void)igr;
    (void)igy; (void)iw_; (void)opac_; (void)plog; (void)uh; (void)usec;
    (void)tso; (void)opf; (void)tsm; (void)tsmin; (void)tf_; (void)tfz;
#endif
    Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
    Real eold = iw_(m,M1_IW_EP,k,j,i);
    Real rkpv = opac_(m,M1_OP_P,k,j,i);
    Real rkev = opac_(m,M1_OP_E,k,j,i);
    Real told = iw_(m,M1_IW_TP,k,j,i);
    Real tnew = told;
    if (rkpv > 0.0 || rkev > 0.0) {
      Real dd = uh(m,IDN,k,j,i);
      Real de0 = iw_(m,M1_IW_DE0,k,j,i);
      bool ok = true;
      // MILESTONE 3g.  The gas has ALREADY been eliminated locally to build the row
      // (step (c)): the linearised energy equation is B_k dT = R_k + c dt rho
      // kappa_E E', whose dT is exactly what put -4 a T_k^3 c dt rho kappa_E/B_k on
      // the diagonal and the rest on the right-hand side.  With
      // implicit_gas_newton the SAME relation supplies T', at no table evaluation at
      // all, instead of re-solving the nonlinear equation from scratch in every
      // pass.  It is one Newton step of that equation, so the Picard loop is now a
      // Newton iteration on the coupled (E,T) system, and its fixed point -- where
      // the loop stops, |dT|/T < implicit_tol -- satisfies
      // R_k + c dt rho kappa_E E' = 0, i.e. the EXACT nonlinear backward-Euler gas
      // equation with e(T) and T^4 evaluated (not linearised) at the final T.
      // The elimination only ADDS to the diagonal of the radiation row (the
      // coefficient 4 a T^3 c dt rho kappa_E emis/B_k is >= 0 whenever B_k > 0), so
      // the M-matrix property of sect. 7 is untouched by it.
      bool done = false;
      if (gnewt) {
        Real bk = iw_(m,igb,k,j,i);
        Real rk = iw_(m,igr,k,j,i);
        Real yk = rk + cl*dt*rkev*enew;
        Real dtk = (bk > 0.0) ? (yk/bk) : 0.0;
        if (bk > 0.0 && fabs(dtk) <= M1_NEWT_TRUST*told && (told + dtk) > 0.0) {
          tnew = told + dtk;
          iw_(m,igy,k,j,i) = fabs(yk);
          done = true;
        }
      }
      if (!done) {
        // the pre-3g bracketed root find: also the per-cell FALLBACK of the Newton
        // update (c_v <= 0, a step outside the trust region, a non-positive T).
        // implicit_tsolve_opac: the fallback root find takes kappa_P, kappa_E at the
        // trial T (M1ImplTemperatureOpac) instead of at the lagged iterate.
        // implicit_thin_freeze: a frozen cell never takes the kappa(T) find
        // implicit_tsolve_opac_mode bit 2: the slope guard sends a cell whose kappa_P
        // falls faster than T^smin at the lagged T to the frozen-opacity find below
        bool tsu = tso && !(tfz && tf_(m,k,j,i) > 0.5);
        if (tsu && (tsm & 2)) {tsu = (M1TsoSlope(opf, dd, told) >= tsmin);}
        const bool tnear = ((tsm & 1) != 0);
        // implicit_tsolve_opac_mode bit 4 (value 8): the stable root nearest the end
        // point of the exact exchange ODE (M1ImplTemperatureOpacOde); wins over bit 1
        const bool trd = ((tsm & 8) != 0);
        if (tsu && (tnear || trd)) {
          if constexpr (decltype(idl)::value) {
            M1EosIdeal th{eos.gamma};
            if (trd) {
              (void) M1ImplTemperatureOpacOde(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                  cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            } else {
              (void) M1ImplTemperatureOpacNear(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                   cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            }
          } else if (usec) {
            Real nmiss = 0.0;
            M1EosCached<decltype(eos), decltype(ec_)> thc{eos, ec_, m, k, j, i, ecnt,
                                                          &nmiss};
            if (trd) {
              (void) M1ImplTemperatureOpacOde(thc, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                  cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            } else {
              (void) M1ImplTemperatureOpacNear(thc, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                   cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            }
            if (gasx) {iw_(m,igm,k,j,i) += nmiss;}
          } else {
            M1EosDirect<decltype(eos)> th{eos};
            if (trd) {
              (void) M1ImplTemperatureOpacOde(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                  cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            } else {
              (void) M1ImplTemperatureOpacNear(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                                   cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            }
          }
        } else if (tsu) {
          if constexpr (decltype(idl)::value) {
            M1EosIdeal th{eos.gamma};
            (void) M1ImplTemperatureOpac(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                         cl*dt*ar, cl*dt, enew + de0, tnew, ok);
          } else if (usec) {
            Real nmiss = 0.0;
            M1EosCached<decltype(eos), decltype(ec_)> thc{eos, ec_, m, k, j, i, ecnt,
                                                          &nmiss};
            (void) M1ImplTemperatureOpac(thc, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                         cl*dt*ar, cl*dt, enew + de0, tnew, ok);
            if (gasx) {iw_(m,igm,k,j,i) += nmiss;}
          } else {
            M1EosDirect<decltype(eos)> th{eos};
            (void) M1ImplTemperatureOpac(th, opf, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                         cl*dt*ar, cl*dt, enew + de0, tnew, ok);
          }
        } else if constexpr (decltype(idl)::value) {
          M1EosIdeal th{eos.gamma};
          (void) M1ImplTemperatureT(th, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                    cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew,
                                    ok);
        } else if (usec) {
          Real nmiss = 0.0;
          M1EosCached<decltype(eos), decltype(ec_)> thc{eos, ec_, m, k, j, i, ecnt,
                                                        &nmiss};
          (void) M1ImplTemperatureT(thc, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                    cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew,
                                    ok);
          if (gasx) {iw_(m,igm,k,j,i) += nmiss;}
        } else {
          (void) M1ImplTemperature(eos, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                   cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew, ok);
        }
        if (!ok) {tnew = told;}
        if (gnewt) {
          iw_(m,igy,k,j,i) = 0.0;
          iw_(m,igf,k,j,i) += 1.0;
        }
      }
    }
    iw_(m,M1_IW_EP,k,j,i) = enew;
    iw_(m,M1_IW_TP,k,j,i) = tnew;
    Real re = fabs(enew - eold)/fmax(fmax(fabs(enew), escale), 1.0e-300);
    Real rt = fabs(tnew - told)/fmax(fabs(tnew), 1.0e-300);
    iw_(m,M1_IW_RES,k,j,i) = fmax(re, rt);
    if (plog) {
      iw_(m,M1_IW_S1,k,j,i) = re;
      iw_(m,M1_IW_S3,k,j,i) = rt;
    }
  };
  par_for_lb("m1_impl_tsolve", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
    tsb(Idl{}, m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ImplicitSolve
//! \brief the whole backward-Euler step: the Picard loop, the tridiagonal column solves,
//! the write-back into u0 and into the gas.

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::DbgEnergyTally
//! \brief dbg_energy_tally: all-rank sums over active cells of (gas IEN, E, fref_wacc,
//! esrc) x cell volume into o[0..3]

void RadiationM1::DbgEnergyTally(Real *o) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.nx1, nj = indcs.nx2, nk = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto ur = u0;
  const bool wa = fref_wsplit && FrefWaccOn();
  DvceArray4D<Real> wc;
  if (wa) {wc = FrefWacc();}
  const bool eso = esrc_on;
  auto es_ = esrc;
  const bool sph = sph_geom;
  auto vol = pmy_pack->pcoord->volume;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const bool hy = fl.on;
  for (int q = 0; q < 4; ++q) {
    Real s = 0.0;
    Kokkos::parallel_reduce("m1_dbg_etally", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            nmb*nk*nj*ni),
    KOKKOS_LAMBDA(const int n, Real &acc) {
      const int m = n/(nk*nj*ni);
      const int k = (n/(nj*ni))%nk + ks, j = (n/ni)%nj + js, i = n%ni + is;
      const Real v = sph ? vol(m,k,j,i) : (mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*
                                           mbsize.d_view(m).dx3);
      Real x = 0.0;
      if (q == 0) {x = hy ? uh(m,IEN,k,j,i) : 0.0;}
      if (q == 1) {x = ur(m,M1_E,k,j,i);}
      if (q == 2) {x = wa ? wc(m,k,j,i) : 0.0;}
      if (q == 3) {x = eso ? es_(m,k,j,i) : 0.0;}
      acc += x*v;
    }, Kokkos::Sum<Real>(s));
    o[q] = s;
  }
#if MPI_PARALLEL_ENABLED
  Real g[4];
  MPI_Allreduce(o, g, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < 4; ++q) {o[q] = g[q];}
#endif
}

TaskStatus RadiationM1::ImplicitSolve(Driver *pdrive, int stage) {
  TmrMark(1);   // implicit_timers: Opacity (and anything since the closure limits)
  if (dbg_etally) {
    Real o[4];
    DbgEnergyTally(o);
    if (global_variable::my_rank == 0) {
      std::printf("ETALLY0 cycle=%d stage=%d t=%.16e dt=%.16e ien=%.16e e=%.16e "
                  "wacc=%.16e esrc=%.16e\n", pmy_pack->pmesh->ncycle, stage,
                  pmy_pack->pmesh->time, dt_sub, o[0], o[1], o[2], o[3]);
    }
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  // closure = vet_sc: the cost of the whole solve is timed against the formal solution
  const bool vetsc = vet_sc;
  auto vc_ = vet_cell;
  // vet_tensor = full: every D_ab of the solve is read from vet_cell (M1DDiag, M1POff)
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  Kokkos::Timer vtimer;
  if (vetsc) {Kokkos::fence(); vtimer.reset();}
  DetTrace("solve_in");

  auto u0_ = u0;
  auto iw_ = iw;
  auto ifw_ = ifw;
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);
  const bool berth = (impl_flux == M1_IFLUX_BERTHON) || (impl_flux == M1_IFLUX_BLEND);
  // MILESTONE 3c: the smooth per-face convex blend.  `bkind` etc. are read only when
  // `blend` is true, and `blend` is false for every 3a/3a2 flux, so those paths keep the
  // arithmetic they had (the weight enters as an exact multiplication by 1.0).
  const bool blend = (impl_flux == M1_IFLUX_BLEND);
  const int bkind = impl_blend, bfm = impl_blend_fmode;
  const bool bdis = blend && (impl_blend_mode == M1_IBMODE_DISSIP);
  const Real btau0 = impl_blend_tau0;
  const Real bflo = impl_blend_flo, bfhi = impl_blend_fhi;
  const bool plmdc = (impl_recon == M1_IRECON_PLMDC);
  const Real rwin = impl_recon_w;
  const Real rfl_ = impl_res_floor;
  const bool rfreeze = impl_recon_freeze;
  const int rnpass = impl_recon_npass;
  auto f0_ = f0x1;
  // F0^n, the face flux at the START of the step: the Picard loop overwrites f0x1 with
  // each new iterate, so the backward-Euler right-hand side needs its own copy
  auto f0n_ = f0x1n;
  Kokkos::deep_copy(DevExeSpace(), f0x1n, f0x1);
  // MILESTONE 3b phase B: the transverse couplings.  `trans` is false for every 3a/3a2/3c
  // configuration and for a 1-D mesh, so those paths keep their arithmetic bit for bit.
  const bool trans = trans_on;
  const bool thrd = trans_x3;
  // MILESTONE 3b phase C: implicit_solver = bicgstab.  False for line_jacobi and for
  // every 3a/3a2/3c configuration, so those paths keep their arithmetic bit for bit.
  const bool bicg = bicg_on;
  // MILESTONE 3b phase D.  Every step STARTS in the configured off-diagonal mode; the
  // positivity fallback below may drop this step to `none` (the operator is not an
  // M-matrix, so E' > 0 is no longer guaranteed by construction).
  od_now = impl_offdiag;
  vimp_now = impl_vimp;
  muscl_now = impl_muscl;   // xthinfix-1009 Fix B
  muscl_nbuild = 0;
  // vet_col_lat: the D_r,lat term is on at the start of every step (the operator
  // form may drop it for the rest of the step, positivity below)
  vlat_now = vlat_on && (vlat_odm > 0);
  // the closure under-relaxation and the start-of-step closure freeze.  Both are inert
  // at their defaults (w = 1, lag = pass), so phase C arithmetic is untouched.
  const Real crw = impl_crelax;
  const bool crthin = impl_crelax_thin;
  const bool clagst = impl_clag_step;
  // DIAGNOSTIC dbg_tensor = frozen | tilt: the stored tensor axis survives the step reset
  const bool tpers = (dbg_tensor == 1 || dbg_tensor == 2) && dbg_tensor_init;
  // closure = tau (rad_m1_tau.cpp): (chi, n) from the column optical depth, built on the
  // first pass of the step and read by step (b) on every pass
  const bool tauc = tau_closure;
  if (tauc && !tau_ready) {TauClosureInit();}
  auto tt_ = tau_ten;
  auto f2_ = f0x2;
  auto f3_ = f0x3;
  if (trans) {
    if (cs_geom && thrd) {CubedSeamFaceAverage();}   // STAGE CS1 (C5)
    Kokkos::deep_copy(DevExeSpace(), f0x2n, f0x2);
    if (thrd) {Kokkos::deep_copy(DevExeSpace(), f0x3n, f0x3);}
  }
  // <rad_m1>/time_scheme = hesdirk2 (rad_m1_time2.cpp).  A STAGE solve takes the stage
  // start state (u0, f0x*, the hydro u0) as its first iterate and as the only state the
  // EOS and the opacities see; its OLD vector is that state plus t2inc, which can be
  // non-physical.  Under time_scheme = be t2st is false and nothing below moves.
  const int t2s = t2_solve;
  const bool t2st = (t2s == M1_T2S_STAGE1) || (t2s == M1_T2S_STAGE2);
  // time2_lin_tol: the stage solves' linear tolerance (put back at the end of the solve)
  t2_lin_save = impl_lin_tol;
  if (t2st) {impl_lin_tol = (t2_lin_tol > 0.0) ? t2_lin_tol : (t2_lin_fac*impl_lin_tol);}
  auto t2i_ = t2inc;
  const bool t2k = (t2s != M1_T2S_NONE);
  const bool t2vs = t2st && (t2_afmode != 0) && impl_vimp && trans;
  const bool t2afc = (t2_afmode == 2);
  const bool t2dav = t2vs && (t2_afmode == 1);   // the DA split is needed
  const int t2da = impl_vimp ? (iw_vimp + M1_IV_DA) : 0;
  auto kk_ = (t2s == M1_T2S_STAGE1) ? t2k2 : t2k1;
  if (t2st) {
    par_for("m1_t2_f1n", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      f0n_(m,k,j,i) += t2i_(m,M1_T2_F1,k,j,i);
    });
    if (trans) {
      auto f2n_ = f0x2n;
      par_for("m1_t2_f2n", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        f2n_(m,k,j,i) += t2i_(m,M1_T2_F1+1,k,j,i);
      });
      if (thrd) {
        auto f3n_ = f0x3n;
        par_for("m1_t2_f3n", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          f3n_(m,k,j,i) += t2i_(m,M1_T2_F1+2,k,j,i);
        });
      }
    }
  }
  auto &mbbcs = cs_geom ? m1bcs : pmy_pack->pmb->mb_bcs;   // CS1: seams open
  auto opac_ = opac;
  auto &mbsize = pmy_pack->pmb->mb_size;
  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  Real ar = arad;
  Real dt = dt_sub;
  bool edd = eddington;
  const int chk = chi_kind;
  bool ovc = source_ovc;
  bool feedback = gas_feedback;
  // MILESTONE 3b phase E: the DEBUG switches of the explicit coupling are honoured here
  // too (they used to be hard-wired to `true` on this branch).  Both default to `true`,
  // so every earlier configuration is bitwise unchanged; what they buy is the ability to
  // ask a multi-D run which HALF of the gas coupling drives a flow -- the radiative FORCE
  // (dbg_gas_force) or the energy exchange (dbg_gas_heat) -- exactly as the 1-D
  // pulsation diagnosis did.  `opac_freeze` is still not on this branch.
  const bool dbgf = dbg_gas_force;
  const bool dbgh = dbg_gas_heat;
  const bool dbgft = dbg_gas_force_trans;
  // LIMIT 3 of the 3a findings.  A physical boundary face hands its whole
  // dt (rho k_t)_f F0_f/c to its ONE interior cell in 3a, so that cell receives 1.5
  // face-shares of radiative force where every interior cell receives 1.0; the residual
  // is a steady force the well-balanced reference a_rad_ref does not carry, and it drives
  // the 10.3 v_MLT bottom-cell flow of I8.  With implicit_bmom_half the boundary face
  // gives HALF, like any other face: the cell-averaged radiative force is then
  // (rho k_t F/c) with F the mean of the cell's two faces, everywhere.  The other half
  // leaves the domain with the radiation, which is where it physically goes.
  const bool bmhalf = impl_bmom_half;
  bool fref = (force_ref == M1_FREF_WB_ARAD);
  // force_reference_work = split (ke-dt-0926): the gas gets only the residual's work here
  const bool fws = fref && fref_wsplit;
  auto aref_ = arad_ref;
  Real mq = marshak_q;
  // vet_col_surface_q (rad_m1_vetcol.cpp): the OUTER x1 Marshak q of each column from
  // its formal solution; the branches below shadow mq with it (off: mq itself)
  const bool vqs = (vet_col && vcol_sq) || hr_q;   // hrup-1009: + half-range
  auto vq_ = vcol_q;
  const Real mqo = marshak_q;
  int bclo = ibc_x1min, bchi = ibc_x1max;
  Real fxlo = iflux_x1min, fxhi = iflux_x1max;
  Real eblo = iebath_x1min, ebhi = iebath_x1max;
  const bool badv = impl_bc_advect;
  bool cyclic = (bclo == M1_IBC_PERIODIC);
  // MILESTONE 3b, LIMIT 4.  With more than one MeshBlock along x1 a block is at a
  // PHYSICAL x1 boundary only when it sits at the corresponding end of its stack; the
  // faces it shares with a stack neighbour are ordinary interior faces whose other cell
  // is a GHOST cell, filled by ImplicitX1Halo with the very numbers the neighbour
  // computed.  With part_nblk == 1 every block is both ends and nothing below changes.
  const int nblkx1 = part_nblk;
  auto pos_ = part_pos;
  const int nlay_ = (part_nblk > 1) ? part_nlay : 0;

  // m1-mhd: <hydro> or <mhd> (FluidRef); the name is kept, it means "a fluid exists"
  FluidRef flr = FluidRef::Get(pmy_pack);
  const bool have_hydro = flr.on;
  const bool src_on = have_hydro && coupling && dbgh && !opac_zero;
  // MILESTONE 3g: the gas-radiation energy coupling.  Both are false by default and
  // every branch they guard is then dead, so the pre-3g arithmetic is untouched.
  const bool gnewt = impl_gas_newton && src_on;
  const bool usec = impl_eos_cache && have_hydro;
  const bool gasx = (iw_gas >= 0);
  const int igb = gasx ? (iw_gas + M1_IWG_BK) : 0;
  const int igr = gasx ? (iw_gas + M1_IWG_RK) : 0;
  const int igy = gasx ? (iw_gas + M1_IWG_YR) : 0;
  const int igf = gasx ? (iw_gas + M1_IWG_FB) : 0;
  const int igm = gasx ? (iw_gas + M1_IWG_MS) : 0;
  const int ecnt = impl_ecnt;
  auto ec_ = ecache;
  auto uh = have_hydro ? flr.u0 : u0;
  const bool etg = have_hydro ? flr.use_etotgrav : false;
  auto phicc = have_hydro ? flr.phicc0 : arad_ref;
  // m1-mhd: |B|^2/2 of the current field (refreshed by Opacity, before this task)
  const bool mhd = fl_mhd;
  auto emag_ = emag0;

  // the optional explicit energy source (esrc, default off): dt (chat/c) esrc joins the
  // OLD vector of E (and, below, the hesdirk2 slope)
  const bool eso = esrc_on;
  auto es_ = esrc;
  const Real dtes = dt*chat/c_light;
  //-------------------------------------------------------------------------- start state
  // STAGE CS2 (C3): the cell F of the work array holds FACE-NORMAL components (what the
  // per-pass rebuild from the face fluxes gives); u0 holds them covariant on cs
  const bool csf0 = cs_geom;
  auto cclf0 = pmy_pack->pcoord->cos_cell;
  par_for("m1_impl_i0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = fmax(u0_(m,M1_E,k,j,i), efl);
    iw_(m,M1_IW_EN,k,j,i) = t2st ? (u0_(m,M1_E,k,j,i) + t2i_(m,M1_T2_E,k,j,i)) : e;
    if (eso) {iw_(m,M1_IW_EN,k,j,i) += dtes*es_(m,k,j,i);}
    iw_(m,M1_IW_EP,k,j,i) = e;
    iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
    iw_(m,M1_IW_V1,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
    iw_(m,M1_IW_DE0,k,j,i) = 0.0;
    iw_(m,M1_IW_G0,k,j,i) = 0.0;
    iw_(m,M1_IW_TP,k,j,i) = 0.0;
    iw_(m,M1_IW_EGN,k,j,i) = 0.0;
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
    if (gasx) {
      iw_(m,igy,k,j,i) = 0.0;
      iw_(m,igf,k,j,i) = 0.0;
      iw_(m,igm,k,j,i) = 0.0;
    }
    if (trans) {
      iw_(m,M1_IW_V2,k,j,i) = 0.0;
      iw_(m,M1_IW_V3,k,j,i) = 0.0;
      iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F2,k,j,i);
      iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F3,k,j,i);
      if (csf0) {
        const Real c = cclf0(m,k,j), sn = sqrt(1.0 - c*c);
        const Real a = u0_(m,M1_F2,k,j,i), b = u0_(m,M1_F3,k,j,i);
        iw_(m,M1_IW_F2,k,j,i) = (a - c*b)/sn;
        iw_(m,M1_IW_F3,k,j,i) = (b - c*a)/sn;
      }
      if (!tpers) {
        iw_(m,M1_IW_N1,k,j,i) = 0.0;
        iw_(m,M1_IW_N2,k,j,i) = 0.0;
        iw_(m,M1_IW_N3,k,j,i) = 0.0;
      }
      iw_(m,M1_IW_A2,k,j,i) = 0.0;
      iw_(m,M1_IW_A3,k,j,i) = 0.0;
      iw_(m,M1_IW_TDIA,k,j,i) = 0.0;
      iw_(m,M1_IW_TRHS,k,j,i) = 0.0;
      iw_(m,M1_IW_LRES,k,j,i) = 0.0;
    }
  });

  if (have_hydro) {
    auto eos = flr.eos;
    // STAGE CS2 (C3): on the cubed sphere the hydro momentum is COVARIANT on the panel
    // basis, so the kinetic energy is 0.5 m_a v^a with the metric (gnomonic_raisevel.hpp)
    const bool csk = cs_geom;
    auto ccl = pmy_pack->pcoord->cos_cell;
    par_for_lb("m1_impl_i1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      Real ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                       SQR(uh(m,IM3,k,j,i)))*idd;
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))
               *idd;
      }
      Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      if (mhd) egrv += emag_(m,k,j,i);
      Real eg = uh(m,IEN,k,j,i) - ekin - egrv;
      iw_(m,M1_IW_TP,k,j,i) = eos.Temperature(dd, fmax(eg, 1.0e-300));
      if (t2st) {
        // the start T is set; now the gas becomes the OLD vector of the stage solve
        uh(m,IM1,k,j,i) += t2i_(m,M1_T2_M1,k,j,i);
        uh(m,IM2,k,j,i) += t2i_(m,M1_T2_M1+1,k,j,i);
        uh(m,IM3,k,j,i) += t2i_(m,M1_T2_M1+2,k,j,i);
        uh(m,IEN,k,j,i) += t2i_(m,M1_T2_EN,k,j,i);
        ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                    SQR(uh(m,IM3,k,j,i)))*idd;
        if (csk) {
          // cs-hydro-energy (10-05): the METRIC kinetic energy here too, as above and as
          // the write-back adds back (m1_impl_wb).  The orthonormal form left EGN short
          // by the cross term, and the write-back then created that difference as total
          // energy at every hesdirk2 stage solve (a per-step leak, not truncation).
          const Real c = ccl(m,k,j);
          const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
          ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))
                 *idd;
        }
        eg = uh(m,IEN,k,j,i) - ekin - egrv;
      }
      iw_(m,M1_IW_EGN,k,j,i) = eg;
      iw_(m,M1_IW_V1,k,j,i) = uh(m,IM1,k,j,i)*idd;
      if (trans) {
        iw_(m,M1_IW_V2,k,j,i) = uh(m,IM2,k,j,i)*idd;
        iw_(m,M1_IW_V3,k,j,i) = uh(m,IM3,k,j,i)*idd;
      }
      if (csk && trans) {
        // STAGE CS2 (C3): the transverse V of the work array are the FACE-NORMAL
        // velocities the face fluxes advect with, v.n_xi = s v^xi, v.n_eta = s v^eta
        // (v^a = g^ab m_b/rho contravariant, s = sin_cell)
        const Real c = ccl(m,k,j), s2 = 1.0 - c*c, sn = sqrt(s2);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        iw_(m,M1_IW_V2,k,j,i) = sn*(m2 - c*m3)*idd/s2;
        iw_(m,M1_IW_V3,k,j,i) = sn*(m3 - c*m2)*idd/s2;
      }
    });
    // MILESTONE 3g: the FROZEN-DENSITY e(T) cache.  rho does not move over the step, so
    // the density direction of the tabulated energy surface is collapsed ONCE here and
    // every pass of the Picard loop reads a 1-D cubic in ln T instead of the 2-D table.
    if (usec) {
      auto eos = flr.eos;
      par_for("m1_impl_ecb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        M1EosCacheBuild(eos, ec_, m, k, j, i, ecnt, uh(m,IDN,k,j,i),
                        iw_(m,M1_IW_TP,k,j,i));
      });
    }
  }

  // The SCALE of the Picard convergence test.  With implicit_res_floor = 0 (the 3a
  // default) the residual is the pure relative change |dE|/E, which in a run with a large
  // dynamic range is dominated by cells many orders below the peak: on gate I6 with
  // implicit_recon = plm_dc the tail cells sit 9 orders under the maximum and keep the
  // reported residual above the tolerance for ever, although the SOLUTION is converged
  // (maxit 30 and maxit 100 give the same amplitude to five digits).  A positive
  // implicit_res_floor scales those cells by the column peak instead,
  // res = |dE|/max(E, implicit_res_floor*max(E)).
  Real emax0 = 0.0;
  if (rfl_ > 0.0) {
    Kokkos::parallel_reduce("m1_impl_emax",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = iw_(m,M1_IW_EN,k,j,i);
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(emax0));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&emax0, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    emax0 = g;}
#endif
  }
  const Real escale = rfl_*emax0;

  // the SCALE the TRUE linear residual is measured against: the max norm of the
  // right-hand side of the 7-point system, which is E^n plus the (small) source terms.
  Real rhsmax = 1.0;
  if (trans) {
    Real rm = 0.0;
    Kokkos::parallel_reduce("m1_impl_rhsmax",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = fabs(iw_(m,M1_IW_EN,k,j,i));
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(rm));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&rm, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    rm = g;}
#endif
    rhsmax = fmax(rm, 1.0e-300);
  }
  Real lresid = 0.0;

  // the x1 ghost layers the partitioned solve reads.  E of the start-of-step state is
  // sent once here; the six lagged quantities and the new E are sent inside the loop.
  // Under transport = implicit the SIX-neighbour exchange of ImplicitTransverseHalo
  // carries all of that and seven quantities more, so the hand-rolled x1 halo is not run.
  if (trans) {
    ImplicitTransverseHalo(M1_NHALO_T);
  } else {
    ImplicitX1Halo(true);
  }

  // closure = vet_sc: the formal solution of the start-of-step state.  (chi, n) are
  // then read by step (b) on every Picard pass: the tensor is lagged by one hydro step.
  TmrMark(2);
  if (vetsc) {
    if (t2s == M1_T2S_STAGE1) {
      Time2VetStart();          // at U^n, then D* extrapolated to the stage time
    } else if (t2s != M1_T2S_STAGE2) {
      VetShortChar();           // stage 2 keeps D* of stage 1
      if (t2s == M1_T2S_BESTORE) {t2_vprev = false;}
    }
  }
  // closure = vet_col (rad_m1_vetcol.cpp): the per-column formal solution of the
  // start-of-step state, here for the same reason (T^n, before the predictor moves it)
  // m1-sph2: under time_scheme = hesdirk2 the tensor (and the surface q) is built ONCE
  // per step, at U^n, by the stage-1 solve (Time2VetStart, as vet_sc), and the stage-2
  // solve keeps it; a backward-Euler step builds it here as before
  // m1-sp-order2b: time2_vet_col (rad_m1_time2.cpp) moves the stage builds to t^{n+1}
  if (vet_col) {
    if (t2s == M1_T2S_STAGE1) {
      if (t2_vcmode == 1 || t2_vcmode == 2) {
        Time2VetColAt(1);
      } else {
        Time2VetStart();
        if (t2_vcmode == 3) {
          Time2VetColExtrap();
        }
      }
    } else if (t2s == M1_T2S_STAGE2) {
      if (t2_vcmode == 2) {
        Time2VetColAt(2);
      }
    } else {
      VetColBuild();
      if (t2s == M1_T2S_BESTORE) {
        t2_vcprev = false;
      }
    }
  }
  TmrMark(3);

  // implicit_predictor = step: start the Picard loop from the previous step's implicit
  // increment, scaled by dt/dt_prev.  Only the STARTING POINT moves: E^n (M1_IW_EN),
  // e^n, the EOS cache window and every scale below are those of the step, and the
  // closures it is allowed with do not read the iterate (the Eddington D_ab does not
  // depend on n, vet_sc and tau read their own arrays), so the fixed point is unchanged.
  // It runs AFTER the vet_sc formal solution, which reads T^n (M1_IW_TP) as its source,
  // and the moved E is then sent to the ghost cells.
  // sp-blend-1008: on the multi-D wedge with a fixed-tensor closure the frozen face
  // coefficients read f(D_rr) of the step's lagged closure, not the iterate, so the
  // predicted start does not enter them and the predictor stays on
  const bool spfx = sph_geom && trans_on && !edd && (vetsc || tauc);
  const bool pred = impl_pred && (edd || vetsc || tauc) && !(aphll && rfreeze && !spfx);
  if (impl_pred && (static_cast<int>(ipred.extent(0)) != nmb1 + 1)) {
    pred_ok = false;   // the pack changed size (AMR): start cold
  }
  // hesdirk2: the stage-2 solve keeps its own increment (ipred2); increments are
  // measured from the stage START state (EP here), which under be is E^n = EN
  const bool p2 = (t2s == M1_T2S_STAGE2);
  auto pd_ = p2 ? ipred2 : ipred;
  bool pred_started = false;
  if (pred) {
    const bool pok = p2 ? (pred2_ok && (pred2_dt > 0.0)) : (pred_ok && (pred_dt > 0.0));
    pred_started = pok;
    const Real rat = pok ? (dt/(p2 ? pred2_dt : pred_dt)) : 0.0;
    const bool hh = have_hydro;
    // implicit_predictor_order = 2: + dt dt1 h (h = 0 until two increments are known)
    const bool po2 = (impl_pord == 2);
    const Real r2 = po2 ? dt*(p2 ? pred2_dt : pred_dt) : 0.0;
    par_for("m1_impl_pred", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real tn = iw_(m,M1_IW_TP,k,j,i);
      pd_(m,2,k,j,i) = tn;
      if (!pok) {return;}
      Real ep = (t2st ? iw_(m,M1_IW_EP,k,j,i) : iw_(m,M1_IW_EN,k,j,i))
                + rat*pd_(m,0,k,j,i);
      if (po2) {ep += r2*pd_(m,3,k,j,i);}
      if (ep > efl) {iw_(m,M1_IW_EP,k,j,i) = ep;}
      if (hh) {
        Real tp = tn + rat*pd_(m,1,k,j,i);
        if (po2) {tp += r2*pd_(m,4,k,j,i);}
        if (tp > 0.5*tn && tp < 2.0*tn) {iw_(m,M1_IW_TP,k,j,i) = tp;}
      }
    });
    if (pok) {
      if (trans) {ImplicitTransverseHalo(1);} else {ImplicitX1Halo(true);}
    }
  }

  // implicit_mr_every (rad_m1_mr.cpp), implicit_mr_peq: the stage-A solve of a
  // multi-rate step starts from the LOCAL equilibrium of every cell -- the gas-radiation
  // exchange of the stage alone, backward Euler, no transport:
  //   E' = (E0 + ap a T'^4)/(1 + ae),
  //   rho e(T') + ap/(1+ae) a T'^4 = rho e0 + ae/(1+ae) E0
  // (ap, ae = c g Delta rho kappa_P, kappa_E).  The hydro steps of the window leave the
  // gas out of equilibrium with E, which the extrapolated increment of the last window
  // does not know.  Only the starting point moves; the fixed point is unchanged.
  if (mr_on && mr_peq && t2s == M1_T2S_STAGE1 && src_on) {
    auto eos = flr.eos;
    const Real cdt = cl*dt;
    par_for_lb("m1_mr_peq", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
      const Real ap = cdt*opac_(m,M1_OP_P,k,j,i);
      const Real ae = cdt*opac_(m,M1_OP_E,k,j,i);
      if (!(ap > 0.0) || !(ae > 0.0)) return;
      const Real e0 = fmax(u0_(m,M1_E,k,j,i), efl);
      const Real t0 = iw_(m,M1_IW_TP,k,j,i);
      const Real iae = 1.0/(1.0 + ae);
      Real tq = t0;
      bool ok = true;
      (void) M1ImplTemperature(eos, uh(m,IDN,k,j,i), t0, iw_(m,M1_IW_EGN,k,j,i),
                               ap*iae*ar, ae*iae*e0, tq, ok);
      if (!ok || !(tq > 0.5*t0 && tq < 2.0*t0)) return;
      const Real t2 = tq*tq;
      iw_(m,M1_IW_EP,k,j,i) = fmax((e0 + ap*ar*t2*t2)*iae, efl);
      iw_(m,M1_IW_TP,k,j,i) = tq;
    });
    if (trans) {
      ImplicitTransverseHalo(1);
    } else {
      ImplicitX1Halo(true);
    }
  }

  // MILESTONE 3e: the Anderson histories start empty at every step, and the per-cell
  // scale of the fixed-point vector is frozen at the start-of-step energy (see
  // ImplicitAccelSave).  Nothing here runs under implicit_accel = none.
  const bool accel = (impl_accel != M1_IACC_NONE);
  if (accel) {
    aa_nh = 0;
    aa_head = 0;
    aa_hasp = false;
    aa_fnp = -1.0;
    auto sc_ = aa_sc;
    par_for("m1_acc_scale", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      sc_(m,k,j,i) = t2st ? fmax(u0_(m,M1_E,k,j,i), efl)
                          : fmax(iw_(m,M1_IW_EN,k,j,i), efl);
    });
  }

  //--------------------------------------------------------------------- the Picard loop
  int it = 0;
  Real resid = 0.0;
  bool converged = false;
  // the per-pass log (implicit_picard_log): E and T parts of the Picard residual, the
  // transverse change, the inner iterations and the inner starting residual
  const bool plog = (impl_plog > 0) && (impl_nstep < static_cast<Real>(impl_plog));
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  Real rprev = -1.0;
  const bool rmask = impl_res_mask;   // implicit_res_dmin / _rmax
  Real rexcl_last = 0.0;
  // implicit_one_pass = N (tests_m1/runs_4a_accel).  A solve that starts from the
  // predictor takes its FIRST pass to the full linear tolerance (no Eisenstat-Walker
  // loosening), so that the change of its second pass is the nonlinear (lagged-term)
  // part alone: res_1 = q res_0 with q the Picard contraction.  q is MEASURED, per kind
  // of solve (be / stage 1 / stage 2), on every solve that takes a second pass, and at
  // least every N-th solve of a kind is made to take one.  The other solves are accepted
  // after the first pass when res_0 qe/(1-qe) < implicit_tol, qe = safety * (the larger
  // of the last two measured q), i.e. when the change a second pass would make is
  // bounded below the tolerance.  The state (q, counters) travels in the restart file.
  //   implicit_one_pass_auto (m1-onepass-auto, default true): where no solve is ever
  // accepted after one pass (He box at cfl 0.9, the wedge: +4.6 % for the tighter first
  // passes alone) a kind of solve switches itself to one_pass = 0.  Its eligible solves
  // (onep_el) are counted since the last one accepted after ONE pass (by either test);
  // at auto_window * N of them, i.e. only after a window with no acceptance at all, the
  // kind is switched off: its solves are then exactly those of implicit_one_pass = 0
  // (EW-loosened first pass, no one-pass test, q and the check counter left alone).
  // q is NOT measured while off (the loosened first pass would pollute res_1/res_0);
  // instead, after auto_reprobe switched-off solves the kind is switched on again as a
  // re-probe with NO q (q history cleared), so the first probe solve is a check that
  // measures q afresh, and the counter is set to (auto_window - 1) * N, so that it is
  // off again after one check period unless a solve is accepted, which confirms it.
  // Until a switch-off the logic only counts, so where acceptances keep coming the
  // solves are bitwise those without it.  The state travels in the restart file.
  const int otyp = (t2s == M1_T2S_STAGE1) ? 1 : ((t2s == M1_T2S_STAGE2) ? 2 : 0);
  const bool onep_el = (impl_onep > 0) && pred && pred_started;
  if (onep_el && impl_onep_auto && onep_off[otyp] > 0.5 &&
      onep_actr[otyp] >= static_cast<Real>(impl_onep_arep)) {
    // re-probe
    onep_off[otyp] = 0.0;
    onep_prb[otyp] = 1.0;
    onep_actr[otyp] = static_cast<Real>((impl_onep_awin - 1)*impl_onep);
    onep_qa[otyp] = -1.0;
    onep_qb[otyp] = -1.0;
    onep_cnt[otyp] = 0.0;
    impl_onep_nprb += 1.0;
  }
  const bool onep = onep_el && !(impl_onep_auto && onep_off[otyp] > 0.5);
  const bool ocheck = onep && ((onep_qa[otyp] < 0.0) ||
                               (onep_cnt[otyp] >= static_cast<Real>(impl_onep - 1)));
  Real ores0 = -1.0, ores1 = -1.0;
  // implicit_gas_newton_switch: the residual history of this solve (host side)
  std::vector<Real> rhist;
  const bool strc = (impl_strace > 0) && (impl_strace_n < impl_strace_max);
  int sloc = -1, sloc2 = -1;
  // implicit_gas_newton_switch: set once this solve's residual is detected stalled
  bool gnsw = false;
  const bool gsd = (impl_gn_sw > 0) && gnewt;
  if (gsd) {rhist.reserve(impl_maxit);}
  if (hr_q) {
    // hrup-1009: q = h+_1 of the top cell of each column (M1HrCell), clamped to the
    // vet_col_surface_q range [vet_col_q_min, vet_col_q_max]
    const bool hgd = vgd_on && sph_geom && (vgd_hr.extent_int(0) > 0);
    const int hsrc = vet_sc ? 1 : (hgd ? 2 : 0);
    auto hvc_ = vet_cell;
    auto hgh_ = vgd_hr;
    auto vqw_ = vcol_q;
    auto iwq_ = iw;
    const int ieq = ie;
    const Real qlo = vet_col ? vcol_qmin : 1.0e-3, qhi = vet_col ? vcol_qmax : 1.0;
    const Real q0 = marshak_q;
    par_for("m1_impl_hrq", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      Real hp, hm, fb;
      M1HrCell(hsrc, iwq_, hvc_, hgh_, m, k, j, ieq, 0, 0.0, hp, hm, fb);
      vqw_(m,k,j) = (hp >= 0.0) ? fmin(fmax(hp, qlo), qhi) : q0;
    });
  }
  for (it = 0; it < impl_maxit && !converged; ++it) {
    // vet_gd_rebuild_every = k: the gd tensor from the Picard iterate at passes k, 2k
    if (vgd_on && vlat_ready && vgd_rbe > 0 && it > 0 && (it % vgd_rbe) == 0) {
      VetGdIterRebuild();
    }
    // the gas Newton update of THIS pass (off after a detected stall, see above)
    const bool gnw = gnewt && !gnsw;
    int nin = -1;
    TmrMark((it == 0) ? 4 : 7);
    ew_tight = onep && (it == 0);
    // MILESTONE 3e: x_k, the state this pass maps
    if (accel) {ImplicitAccelSave();}
    // the off-diagonal mode of THIS pass (the positivity fallback can change it)
    const int odm = od_now;
    // the closure moves only after the first pass, and not at all under
    // implicit_closure_lag = step
    const bool dorel = (crw < 1.0) && (it > 0);
    const bool dofreeze = clagst && (it > 0);
    // runs_5c_thinstab: implicit_closure_thin_relax, first pass of the step only
    const bool ctr = (impl_ctrelax > 0.0) && (it == 0);
    const bool ctri = ctr_init;
    const Real ctc = impl_ctrelax;
    const bool ctsph = sph_geom;
    // ...and on a STRETCHED spherical-polar x1 the logical mb_size dx1 is the mean dr:
    // the cell optical depth takes the cell's own width xx1f(i+1) - xx1f(i) there
    // (uniform grids keep the logical dx1, bitwise)
    const bool ctstr = sph_geom && (pmy_pack->pmesh->use_grid_stretch_r ||
                                    pmy_pack->pmesh->use_grid_stretch_r_poly);
    auto ctx1f = pmy_pack->pcoord->xx1f;
    auto cm_ = ctr_mem;
    // DIAGNOSTIC dbg_tensor (VET scaffolding): tkeep = read the stored tensor; ttau =
    // rebuild it from the optical depth (first pass of every step); ttilt = rotate the
    // axis of the tensor computed on the very first pass of the run
    const int tmode = dbg_tensor;
    const bool tkeep = (tmode != 0) && ((it > 0) || (tmode != 3 && dbg_tensor_init));
    const bool ttau = (tmode == 3) && (it == 0);
    const bool ttilt = (tmode == 2) && (it == 0) && !dbg_tensor_init;
    const Real talp = dbg_tensor_tilt;
    const Real tx2min = pmy_pack->pmesh->mesh_size.x2min;
    const Real tx2len = pmy_pack->pmesh->mesh_size.x2max - tx2min;
    if (tmode == 3 && pmy_pack->pmesh->mesh_indcs.nx1 != indcs.nx1) {
      ImplFatal("<rad_m1>/dbg_tensor = tau needs ONE MeshBlock along x1");
    }
    // (a) optional opacity re-evaluation at the current temperature iterate
    // implicit_opac_newton: also at pass 0 (at the iterate T, not T^n), so that every
    // pass's rows carry the derivative at the temperature the opacity was taken at
    const bool opn = impl_opac_newton && have_hydro && !opac_zero;
    // implicit_thin_freeze: mark, at pass 0 from the start-of-solve opacities, the cells
    // with c dt rho kappa_P below the threshold; they skip every opacity update below
    const bool tfz = (impl_thin_frz > 0.0) && have_hydro && !opac_zero;
    if (tfz && it == 0) {
      auto tf_ = thin_frz;
      const Real thr = impl_thin_frz, cdt = c_light*dt;
      par_for("m1_impl_thinfrz", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        tf_(m,k,j,i) = (cdt*opac_(m,M1_OP_P,k,j,i) < thr) ? 1.0 : 0.0;
      });
    }
    // fixbundle-1009 F2: rho kappa_T at T^n (the Opacity task's, copied into iw by the
    // start-state kernel) before pass 0 moves it to the predicted T
    if (impl_face_ktn && it == 0) {
      auto ktn_ = ktn;
      par_for("m1_impl_ktn", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        ktn_(m,k,j,i) = iw_(m,M1_IW_KT,k,j,i);
      });
    }
    if (impl_opac_update && (it > 0 || opn) && have_hydro && !opac_zero) {
      int otype = opacity_type;
      Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kscat = kappa_s;
      Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
      M1OpacTab ot = otab;
      const int opart = dbg_opac_part;
      auto ktd_ = ktd;
      auto tf_ = thin_frz;
      const bool osc = vscat;   // vet_scatter: rho kappa_e -> opac(M1_OP_S)
      const Real soff = impl_opn_soff, smax = impl_opn_smax;
      const bool oscl = (soff > 0.0 || smax > 0.0);
      auto eosv = flr.eos;
      par_for("m1_impl_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (tfz && tf_(m,k,j,i) > 0.5) {
          if (opn) {ktd_(m,k,j,i) = 0.0;}
          return;
        }
        Real d = uh(m,IDN,k,j,i);
        Real t = iw_(m,M1_IW_TP,k,j,i);
        Real op, oe, of, os;
        if (otype == M1_OPAC_TABLE) {
          M1TableOpacities(ot, d, t, op, oe, of, os);
        } else {
          M1Opacities(otype, d, t, kp, kev, kf, kscat, rref, tref, aa, bb, op, oe,
                      of, os);
        }
        if (opart != 1) {
          opac_(m,M1_OP_P,k,j,i) = d*op;
          opac_(m,M1_OP_E,k,j,i) = d*oe;
        }
        if (osc && otype == M1_OPAC_TABLE) {M1ScatterEos(ot, eosv, d, t, of, os);}
        if (opart != 2) {
          opac_(m,M1_OP_T,k,j,i) = d*(of + os);
          iw_(m,M1_IW_KT,k,j,i) = d*(of + os);
          if (osc) {opac_(m,M1_OP_S,k,j,i) = d*os;}
        }
        if (opn) {
          // one-sided difference; the table is bilinear in (log T, log rho), so this is
          // the slope of the cell the iterate sits in (the Jacobian only sets the rate)
          const Real th = t*M1_OPN_H;
          Real op2, oe2, of2, os2;
          if (otype == M1_OPAC_TABLE) {
            M1TableOpacities(ot, d, t + th, op2, oe2, of2, os2);
          } else {
            M1Opacities(otype, d, t + th, kp, kev, kf, kscat, rref, tref, aa, bb, op2,
                        oe2, of2, os2);
          }
          ktd_(m,k,j,i) = d*((of2 + os2) - (of + os))/th;
          if (oscl) {
            // logarithmic slope of kappa_T over the same difference
            const Real kt0 = of + os;
            const Real sl = (kt0 > 0.0) ? fabs((of2 + os2) - kt0)*t/(kt0*th) : 0.0;
            if (soff > 0.0 && sl > soff) {
              ktd_(m,k,j,i) = 0.0;
            } else if (smax > 0.0 && sl > smax) {
              ktd_(m,k,j,i) *= smax/sl;
            }
          }
        }
      });
    }

    if (tauc && it == 0) {TauClosureBuild();}
    // (b) the lagged closure, the enthalpy-flux coefficient, de0 and g0
    const bool vdv = t2st && impl_vimp && t2_fvnew && (it > 0);
    const int ivd = impl_vimp ? (iw_vimp + M1_IV_DV) : 0;
    const bool rcp = impl_real_couple;
    // rad-beam-1008 implicit_recon_dgpass: the lagged reduced flux of a cell is made from
    // the LOW-ORDER part of its face fluxes (the stored face flux minus the plm deferred
    // correction): the correction's face E is not the donor's, and its ratio to the donor
    // E (< 1 on a falling front) slows the front and steepens it into a spike
    const bool rdgx = aphll && plmdc && impl_recon_dgpass;
    const Real clch = c_light/chat;
    // implicit_realisable_coupling on the cubed sphere: the work array's transverse F are
    // FACE-NORMAL; the fix-A kernel forms the covariant pair and the clip's metric norm
    const bool rcs = cs_geom;
    auto rccl = pmy_pack->pcoord->cos_cell;
    auto rcsn = pmy_pack->pcoord->sin_cell;
    if (!rcp) {
      par_for("m1_impl_lag", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
#define M1_RCP 0
#include "rad_m1/rad_m1_impl_lag_kernel.hpp"
#undef M1_RCP
      });
    } else {
      // implicit_realisable_coupling: the same kernel with the key on
      par_for("m1_impl_lag_rc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
#define M1_RCP 1
#include "rad_m1/rad_m1_impl_lag_kernel.hpp"  // NOLINT(build/include)
#undef M1_RCP
      });
    }
    // m1-positivity, implicit_g0_limit = w: the lagged
    // g0 = rho (kappa_E E0 - kappa_P a T^4) of the iterate enters every face-flux
    // equation as - c dt v_f g0_f, a lagged, explicit, centred term.  Where the
    // exchange is stiff and the iterate is out of equilibrium (Fe-bump plume:
    // kappa_P >> kappa_F, c dt rho kappa_P ~ 1e5) chat dt g0 is orders of magnitude
    // above the energy the cell can exchange in the step, and that term alone drives the
    // solved E negative (He presn wedge, dbg_t2_admiss).  At a converged backward-Euler
    // state chat dt g0 = -q, the energy actually exchanged.  Own kernel, so the lag
    // kernel above is untouched.
    // implicit_g0_exchange: from the second pass on, g0 is taken from the EXCHANGE the
    // previous pass's linearised source row gives at the iterate, g0 = -(SRCR - SRCB E^k)
    // /(chat dt), instead of rho (kappa_E E^k - kappa_P a T_k^4).  The two are the same
    // at the fixed point (the write-back's q = SRCR - SRCB E' is the exchanged energy),
    // but
    // in a radiation-dominated, stiffly coupled cell the pointwise form is the difference
    // of two numbers ~1e5 x its converged value and a relative error of T_k of 1e-6 moves
    // it by more than the cell's energy per step.  The first pass keeps the pointwise
    // form (with the clip, if on).
    // The clip (implicit_g0_limit) acts on the FIRST pass only: it is a guard for the
    // pass that has no exchange yet, and is not applied where the fixed point is decided
    // (measured on the He presn wedge: clipping every pass left 1.7e7 cell-passes clipped
    // and 5 of 6 solves NON-CONVERGED; the exchange form alone converged in 17 passes).
    const bool gex = impl_g0_exch && src_on && (it > 0);
    const bool gcl = (impl_g0_lim > 0.0) && src_on && (it == 0);
    if (gcl || gex) {
      const Real w = impl_g0_lim;
      const bool hh = have_hydro;
      auto pc_ = pos_cnt_d;
      par_for("m1_impl_g0lim", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (gex) {
          iw_(m,M1_IW_G0,k,j,i) = -(iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*
                                    fmax(iw_(m,M1_IW_EP,k,j,i), efl))/(ch*dt);
        }
        if (!gcl) return;
        const Real ea = fmax(fmax(iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EN,k,j,i)), 0.0)
                        + (hh ? (ch/cl)*fmax(iw_(m,M1_IW_EGN,k,j,i), 0.0) : 0.0);
        const Real gm = w*ea/(ch*dt);
        const Real g = iw_(m,M1_IW_G0,k,j,i);
        if (fabs(g) > gm) {
          iw_(m,M1_IW_G0,k,j,i) = (g > 0.0) ? gm : -gm;
          Kokkos::atomic_add(&pc_(M1_POS_G0), 1.0);
        }
      });
    }
    if ((tmode == 1 || tmode == 2) && it == 0) {dbg_tensor_init = true;}
    if (ctr) {ctr_init = true;}

    // the x1 halo of the LAGGED quantities (w, a, g0, v1, the comoving reduced flux and
    // the transport opacity).  Every face of the stack is then assembled by both of its
    // blocks from bit-identical numbers.
    if (trans) {
      // ...and only the M1_NHALO_Q of them a pass can still move when the
      // closure is frozen for the step (see M1HaloCompT).
      ImplicitTransverseHalo(dofreeze ? M1_NHALO_Q : M1_NHALO_T);
      // (b1) the lagged transverse operator: the x2/x3 face fluxes of this iterate, their
      // diagonal contribution to the matrix and their lagged right-hand side.
      if (t2dav && (it == 0 || !(edd || vetsc || dofreeze))) {
        for (int c = 0; c < 3; ++c) {ImplicitHaloExchange(1, t2da + c);}
      }
      ImplicitTransverseTerms(it == 0);
    } else {
      ImplicitX1Halo(false);
    }

    // (b2) the FACE coefficients of the asymptotic-preserving HLL blend.  Nothing here
    // runs under implicit_flux = central, where ifw stays identically zero and the row
    // assembled below is bitwise the 3a one.
    // The deferred correction is by default evaluated ONCE per step, at the start-of-step
    // state (implicit_recon_lag = step).  Recomputing it every Picard pass
    // (= picard) makes the loop a limit cycle: the plm limiter keeps switching on a few
    // cells and the strict tolerance is never reached, at 30 iterations per step against
    // 2, although the answer is the same to five digits.  The correction is a lagged,
    // explicit term in any case, so evaluating it at E^n costs nothing in order.
    const int iter = it;
    const bool doface = aphll && (it == 0 || !rfreeze);
    if (doface) {
      const bool dodg = plmdc && (it == 0
                                  || (!rfreeze && (rnpass <= 0 || it < rnpass)));
      // sp-blend-1008 (spherical-polar wedge): the face optical depth is the SAME face
      // kappa_T the sp row's theta uses (M1FaceAvgX1, distance-weighted under
      // implicit_face_weight = distance) times the centroid distance dxface.x1f, and that
      // distance replaces dx1 wherever a face length enters.  The HLL coefficients
      // themselves are face FLUXES per unit E (no length); the row multiplies them by
      // dt A_f/V_i.  The Cartesian expressions are untouched (sphf false).
      const bool sphf = sph_geom;
      const bool fwdf = sph_geom && impl_face_wdist;
      auto cdxff = pmy_pack->pcoord->dxface;
      auto cx1ff = pmy_pack->pcoord->xx1f;
      // sp-blend-1008: with a FIXED-TENSOR closure (vet_col, vet_sc, tau) on the multi-D
      // wedge the reduced flux the HLL part and the weight use is the one the LAGGED
      // closure implies, f(D_rr) of the M1 (Levermore) relation chi = (3 + 4 f^2)/(5 + 2
      // sqrt(4 - 3 f^2)) inverted: sqrt(4 - 3 f^2) = (5 - 3 chi)/2, with the sign of the
      // lagged face flux.  The face-flux ratio alone is not relaxed by anything in the
      // thin limit (berthon keeps whatever f it is given), so a beam from a photosphere
      // kept f ~ 0.89 where the formal solution says 0.97 (AG Car A column, 3 R_ph).
      const bool vfix = sphf && trans && !edd && (vetsc || tauc);
      const bool fktn = impl_face_ktn && (it == 0);   // F2: kappa_T at T^n
      auto ktn_ = ktn;
      const bool bvx = bvec_x1 && !impl_beam_hr;   // blendall-1009
      const bool bfsx = impl_beam_fs;
      // hrup-1009: the half-range face flux on the x1 faces too (every geometry)
      const bool hrx = impl_beam_hr;
      const Real balph = impl_blend_alpha, br0 = impl_blend_r0;
      const Real bx0 = impl_blend_xthin;
      const int bxmode = impl_blend_xthin_mode;
      const Real bxwmin = impl_blend_xthin_wmin, bxr0 = impl_blend_xthin_r0;
      const bool hgd = hrx && vgd_on && sph_geom && (vgd_hr.extent_int(0) > 0);
      const int hsrc = !hrx ? 0 : (vetsc ? 1 : (hgd ? 2 : (!impl_hr_model ? -1 :
                                                            (trans ? 0 : 3))));
      auto hvc_ = vet_cell;
      auto hgh_ = vgd_hr;
      par_for("m1_impl_aphll", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        int ipos = pos_.d_view(m);
        bool phys = (((i == is) && (ipos == 0)) ||
                     ((i == ie+1) && (ipos == nblkx1-1))) && !cyclic;
        if (phys) {
          // a physical boundary face: the flux is IMPOSED there (flux / Marshak /
          // reflect / efix), so there is no Riemann problem and no blend.
          ifw_(m,M1_IFW_AL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCR,k,j,i) = 0.0;
          ifw_(m,M1_IFW_DG,k,j,i) = 0.0;
          return;
        }
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real dx = mbsize.d_view(m).dx1;
        if (sphf) {dx = cdxff.x1f(m,k,j,i);}
        if (bvx) {
          // blendall-1009 (Cartesian, implicit_flux_faces = all): the x1 faces take the
          // beam vector of the lagged closure, as the x2/x3 faces do (M1BeamFace)
          const Real kl0 = fktn ? ktn_(m,k,j,im) : iw_(m,M1_IW_KT,k,j,im);
          const Real kr0 = fktn ? ktn_(m,k,j,ip) : iw_(m,M1_IW_KT,k,j,ip);
          Real alb, hlb, hrb;
          const Real fml = bfsx ? iw_(m,M1_IW_S1,k,j,im)
                                : M1FLev(iw_(m,M1_IW_WCHI,k,j,im));
          const Real fmr = bfsx ? iw_(m,M1_IW_S1,k,j,ip)
                                : M1FLev(iw_(m,M1_IW_WCHI,k,j,ip));
          M1BeamFace(fml, iw_(m,M1_IW_N1,k,j,im), fmr, iw_(m,M1_IW_N1,k,j,ip),
                     0.5*(kl0 + kr0)*dx, edd, blend, bkind, bfm, btau0, bflo, bfhi, ch,
                     alb, hlb, hrb);
          ifw_(m,M1_IFW_AL,k,j,i) = alb;
          ifw_(m,M1_IFW_HCL,k,j,i) = hlb;
          ifw_(m,M1_IFW_HCR,k,j,i) = hrb;
          ifw_(m,M1_IFW_DG,k,j,i) = 0.0;
          return;
        }
        Real rfl = iw_(m,M1_IW_RF0,k,j,im);
        Real rfr = iw_(m,M1_IW_RF0,k,j,ip);
        if (vfix) {
          Real cq = fmin(fmax(M1DDiag(iw_,vd_,dfull,m,0,k,j,im), 1.0/3.0), 1.0);
          Real sq = 0.5*(5.0 - 3.0*cq);
          rfl = copysign(sqrt(fmax((4.0 - sq*sq)/3.0, 0.0)), rfl);
          cq = fmin(fmax(M1DDiag(iw_,vd_,dfull,m,0,k,j,ip), 1.0/3.0), 1.0);
          sq = 0.5*(5.0 - 3.0*cq);
          rfr = copysign(sqrt(fmax((4.0 - sq*sq)/3.0, 0.0)), rfr);
        }
        // closed-form M1 wave speeds of the two LAGGED states (1-D: mu = sign f)
        Real bl, br;
        if (edd) {
          br = ch/sqrt(3.0);
          bl = -br;
        } else {
          Real lml, lpl, lmr, lpr;
          M1WaveSpeeds(fabs(rfl), (rfl >= 0.0) ? 1.0 : -1.0, lml, lpl);
          M1WaveSpeeds(fabs(rfr), (rfr >= 0.0) ? 1.0 : -1.0, lmr, lpr);
          bl = ch*fmin(fmin(lml, lmr), 0.0);
          br = ch*fmax(fmax(lpl, lpr), 0.0);
        }
        // alpha: Bloch et al. (2021) eq. 25 with the (1-f^2) guard and the arithmetic
        // face mean of the CELL optical depth.  lp*lm <= 0, so den >= 1 and alpha <= 1.
        const Real ktl = fktn ? ktn_(m,k,j,im) : iw_(m,M1_IW_KT,k,j,im);
        const Real ktr = fktn ? ktn_(m,k,j,ip) : iw_(m,M1_IW_KT,k,j,ip);
        Real tauf = 0.5*(ktl + ktr)*dx;
        if (sphf) {
          tauf = M1FaceAvgX1(ktl, ktr, cx1ff, m, im, ip, fwdf)*dx;
        }
        if (hrx) {
          // hrup-1009: c (h+_L E_L + h-_R E_R), AP-weighted (M1HrCell, M1HrFace)
          Real hpl, hml, hpr, hmr, fbl, fbr;
          M1HrCell(hsrc, iw_, hvc_, hgh_, m, k, j, im, 0, iw_(m,M1_IW_RF0,k,j,im), hpl,
                   hml, fbl);
          M1HrCell(hsrc, iw_, hvc_, hgh_, m, k, j, ip, 0, iw_(m,M1_IW_RF0,k,j,ip), hpr,
                   hmr, fbr);
          Real alb, hlb, hrb;
          M1HrFace(hpl, hmr, fbl, fbr, tauf, blend, bkind, bfm, btau0, bflo, bfhi, balph,
                   ch, iw_(m,M1_IW_EP,k,j,im), iw_(m,M1_IW_EP,k,j,ip), br0, hml, hpr,
                   ch*dt/dx, bx0, bxmode, bxwmin, bxr0,
                   (bxmode == 2) ? fmax(M1HrJFlag(hsrc, iw_, hvc_, m, k, j, im),
                                        M1HrJFlag(hsrc, iw_, hvc_, m, k, j, ip)) : 0.0,
                   alb, hlb, hrb);
          ifw_(m,M1_IFW_AL,k,j,i) = alb;
          ifw_(m,M1_IFW_HCL,k,j,i) = hlb;
          ifw_(m,M1_IFW_HCR,k,j,i) = hrb;
          ifw_(m,M1_IFW_DG,k,j,i) = 0.0;
          return;
        }
        Real al = 1.0;
        if (tauf > 0.0) {
          Real fbar = 0.5*(fabs(rfl) + fabs(rfr));
          Real guard = fmax(1.0 - fbar*fbar, 0.0);
          Real lp = br/ch, lm = bl/ch;
          Real den = 1.0 - 3.0*tauf*guard*lp*lm/(lp - lm + 1.0e-300);
          al = 1.0/fmax(den, 1.0);
        }
        // F_HLL = [b_R c_h f_L E'_L - b_L c_h f_R E'_R + b_R b_L (E'_R - E'_L)]/(b_R-b_L)
        // is linear in E'.  Split into the ADVECTIVE part (the two physical fluxes) and
        // the DISSIPATION (the jump term), because the two carry different weights:
        //
        //   F = alpha F_adv + alpha^2 F_dis + (1 - alpha) F_diff.
        //
        // The dissipation must carry alpha^2 and not alpha.  At piecewise-constant states
        // -- which is what the MATRIX is built from, whatever implicit_recon says -- the
        // HLL dissipation IS the physical diffusion (b_R b_L dE/(b_R-b_L) -> -c dE/3),
        // so weighting it alpha and adding (1-alpha) F_diff on top counts the diffusion
        // TWICE: measured on gate I1 at tau_cell = 1e3, d(sigma^2)/dt came out 2.0022 x
        // the analytic 2D at every CFL.  alpha^2 ~ 1/tau^2 kills it against F_diff
        // ~ 1/tau and leaves the thin limit (alpha -> 1) exactly the plain HLL flux.
        // This is the `alpha2` form of the explicit scheme (rad_m1_closure.hpp), reached
        // here for the same reason.
        //
        // The two fmax()/fmin() are the M-MATRIX GUARDS.  At alpha = 1 they are provably
        // inactive (the HLL consistency condition b_L <= c_h f <= b_R holds on the M1
        // admissible set), but alpha < 1 rescales the two parts differently and the
        // E'_R coefficient can turn positive; the clamp then drops it to zero, which
        // only makes the face flux more upwind and leaves conservation exact (it is one
        // number per face, used with opposite signs by the two cells).
        Real invb = 1.0/(br - bl + 1.0e-300);
        Real adl = br*ch*rfl*invb;        // E'_L coefficient of F_adv
        Real adr = -bl*ch*rfr*invb;       // E'_R coefficient of F_adv
        Real dk = -br*bl*invb;            // >= 0, the dissipation coefficient
        //
        // implicit_flux = berthon drops F_diff altogether and weights BOTH parts of the
        // HLL flux by alpha, which is what alpha was constructed for: alpha (F_adv +
        // F_dis) is the physical diffusion to first order in 1/tau, the advective part
        // supplying the 1/(0.866 tau) that the dissipation alone is short of.  The
        // F_diff weight is then zero, which is what storing AL = 1 below means.
        Real wdis = berth ? al : (al*al);
        // MILESTONE 3c.  w_f in [0,1] from the LAGGED face quantities: the whole face
        // flux is (1 - w_f) F_central + w_f F_berthon (implicit_blend_mode = flux), or
        // the FULL central flux plus w_f times the HLL DISSIPATION alone
        // (= dissipation).  w_f = 1 with mode = flux reproduces `berthon` bitwise and
        // w_f = 0 reproduces `central` bitwise: the multiplications below are by exactly
        // 1.0 or exactly 0.0.  Every input of the weight lives in iw, which the x1 halo
        // of the partitioned solve already carries.
        Real wf = 1.0;
        if (blend) {
          wf = M1BlendWeight(bkind, bfm, tauf, btau0, rfl, rfr, bflo, bfhi);
        }
        // In `dissipation` mode the advective part of the HLL flux is NOT added (the
        // central flux already carries the transport); only the jump term is.
        Real wadv = bdis ? 0.0 : al;
        Real wdsq = bdis ? al : wdis;
        Real ccl = wf*fmax(wadv*adl + wdsq*dk, 0.0);
        Real ccr = wf*fmin(wadv*adr - wdsq*dk, 0.0);
        // how much of the CENTRAL (face-eliminated) flux the row keeps is 1 - AL.
        ifw_(m,M1_IFW_AL,k,j,i) = bdis ? 0.0 : (berth ? wf : al);
        ifw_(m,M1_IFW_HCL,k,j,i) = ccl;
        ifw_(m,M1_IFW_HCR,k,j,i) = ccr;
        // the plm DEFERRED CORRECTION: the difference between the plm and the dc HLL
        // flux at the PREVIOUS iterate.  It goes to the right-hand side, so the matrix
        // stays the low-order M-matrix.  Both E and the comoving reduced flux are
        // reconstructed, with the same limiter the explicit scheme uses; a face whose
        // 4-cell stencil leaves the block falls back to dc (zero correction).
        Real dg = 0.0;
        if (dodg && al > 0.0) {
          int ilo = is - ((ipos > 0) ? nlay_ : 0);
          int ihi = ie + ((ipos < nblkx1-1) ? nlay_ : 0);
          int imm = (im > ilo) ? (im-1) : (cyclic ? ie : -1);
          int ipp = (ip < ihi) ? (ip+1) : (cyclic ? is : -1);
          if (imm >= 0 && ipp >= 0) {
            Real dum;
            Real elp, erp, flp, frp;
            PLM(iw_(m,M1_IW_EP,k,j,imm), iw_(m,M1_IW_EP,k,j,im),
                iw_(m,M1_IW_EP,k,j,ip), elp, dum);
            PLM(iw_(m,M1_IW_EP,k,j,im), iw_(m,M1_IW_EP,k,j,ip),
                iw_(m,M1_IW_EP,k,j,ipp), dum, erp);
            // rad-beam-1008: with a fixed-tensor closure on the wedge the reduced flux
            // is f(D_rr) of the lagged closure (as rfl, rfr above), in all four cells
            auto fcell = [&](const int ii) {
              Real fv = iw_(m,M1_IW_RF0,k,j,ii);
              if (vfix) {
                const Real cq = fmin(fmax(M1DDiag(iw_,vd_,dfull,m,0,k,j,ii), 1.0/3.0),
                                     1.0);
                const Real sq = 0.5*(5.0 - 3.0*cq);
                fv = copysign(sqrt(fmax((4.0 - sq*sq)/3.0, 0.0)), fv);
              }
              return fv;
            };
            PLM(fcell(imm), rfl, rfr, flp, dum);
            PLM(rfl, rfr, fcell(ipp), dum, frp);
            Real ecl = iw_(m,M1_IW_EP,k,j,im), ecr = iw_(m,M1_IW_EP,k,j,ip);
            // the deferred correction applies to the UPWIND part only, so it carries the
            // same weight w_f the upwind part carries (3c).
            Real gp = wf*(wadv*(br*ch*flp*elp - bl*ch*frp*erp)*invb
                          + wdsq*(-dk)*(erp - elp));
            Real gc = ccl*ecl + ccr*ecr;
            // ADMISSIBILITY of the corrected face flux against the DONOR cell.  The
            // reconstructed face energy may exceed the donor cell's own (plm puts
            // E_i (r-1)/(r+1) on top of E_i for a geometric ratio r), and c times that is
            // then faster than the donor can physically emit: the cell drains below what
            // it receives and, on an exponentially falling background, the drain
            // cascades.  Measured on I6 before this clamp: the 1e-4 background of the
            // free-streaming pulse collapsed onto the floor and the peak grew 7x
            // (amplitude ratio 6.98, Picard never converging).  The low-order flux gc
            // already satisfies this bound, so the clamp never removes the whole
            // correction, only the inadmissible part of it.
            Real gmax = wf*ch*ecl, gmin = -wf*ch*ecr;
            gp = fmin(fmax(gp, gmin), gmax);
            // The DEFERRED-CORRECTION WEIGHT.  A deferred correction is a fixed-point
            // iteration x <- A_low^-1 (b + (A_low - A_high) x), and for advection its
            // contraction factor is ~ 2 nu/(1 + nu) with nu = chat dt/dx: it converges
            // only below CFL ~ 1 and DIVERGES above it (measured: Picard never converges
            // at implicit_cfl = 10 and the I6 pulse amplitude comes out 3.85).  The
            // correction is therefore weighted by w = 1/(1 + nu) unless
            // <rad_m1>/implicit_recon_w names a fixed value.  That makes the contraction
            // factor 2 nu/(1 + nu)^2 <= 1/2 at EVERY CFL, and the fixed point a convex
            // blend of the dc and plm fluxes -- still a monotone flux, second-order where
            // w -> 1 (nu << 1, which is where a propagating front is resolved in time at
            // all) and dc where the step is so long that the front is not resolved.
            Real wdc = (rwin > 0.0) ? rwin : (1.0/(1.0 + ch*dt/dx));
            dg = wdc*(gp - gc);
            // UNDER-RELAXATION across the Picard passes.  The plm limiter keeps switching
            // on a handful of cells and the un-relaxed iteration is a small-amplitude
            // limit cycle that never meets the tolerance (30 passes per step against 2,
            // with the answer already right to five digits).  Averaging with the previous
            // pass leaves the fixed point untouched and breaks the cycle.
            if (iter > 0) {dg = 0.5*(dg + ifw_(m,M1_IFW_DG,k,j,i));}
          }
        }
        ifw_(m,M1_IFW_DG,k,j,i) = dg;
      });
    }
    // rad-beam-1008 implicit_recon_dgpass: the plm deferred correction of the berthon
    // part re-made on EVERY Picard pass from the iterate's E with the step-frozen face
    // coefficients: G = HCL E_L^plm + HCR E_R^plm (the frozen reduced flux, plm in E
    // only) minus the dc flux HCL E_L + HCR E_R, clamped to the donor bound, weighted
    // 1/(1 + nu) and averaged with the previous pass (as the picard-lag correction).
    // The face then converges to a convex dc/plm blend of the END-of-step state, not the
    // explicit correction of E^n.
    if (aphll && plmdc && rfreeze && impl_recon_dgpass) {
      const Real rwin2 = impl_recon_w;
      const bool sphf2 = sph_geom;
      auto cdxf2 = pmy_pack->pcoord->dxface;
      const int nlay2 = nlay_;
      par_for("m1_impl_dgpass", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        int ipos = pos_.d_view(m);
        bool phys = (((i == is) && (ipos == 0)) ||
                     ((i == ie+1) && (ipos == nblkx1-1))) && !cyclic;
        if (phys) {return;}
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        int ilo = is - ((ipos > 0) ? nlay2 : 0);
        int ihi = ie + ((ipos < nblkx1-1) ? nlay2 : 0);
        int imm = (im > ilo) ? (im-1) : (cyclic ? ie : -1);
        int ipp = (ip < ihi) ? (ip+1) : (cyclic ? is : -1);
        Real dg = 0.0;
        const Real ccl = ifw_(m,M1_IFW_HCL,k,j,i), ccr = ifw_(m,M1_IFW_HCR,k,j,i);
        const Real wfr = ifw_(m,M1_IFW_AL,k,j,i);
        if (imm >= 0 && ipp >= 0 && (ccl != 0.0 || ccr != 0.0)) {
          Real dum, elp, erp;
          const Real ecl = iw_(m,M1_IW_EP,k,j,im), ecr = iw_(m,M1_IW_EP,k,j,ip);
          PLM(iw_(m,M1_IW_EP,k,j,imm), ecl, ecr, elp, dum);
          PLM(ecl, ecr, iw_(m,M1_IW_EP,k,j,ipp), dum, erp);
          Real gp = ccl*elp + ccr*erp;
          const Real gc = ccl*ecl + ccr*ecr;
          const Real wq = (wfr > 0.0) ? wfr : 1.0;
          gp = fmin(fmax(gp, -wq*ch*ecr), wq*ch*ecl);
          Real dx = mbsize.d_view(m).dx1;
          if (sphf2) {dx = cdxf2.x1f(m,k,j,i);}
          const Real wdc = (rwin2 > 0.0) ? rwin2 : (1.0/(1.0 + ch*dt/dx));
          dg = wdc*(gp - gc);
          if (iter > 0) {dg = 0.5*(dg + ifw_(m,M1_IFW_DG,k,j,i));}
        }
        ifw_(m,M1_IFW_DG,k,j,i) = dg;
      });
    }

    if (strc && it > impl_strace_p0 && it < impl_strace_p0 + impl_strace) {
      (void) ImplicitNCDump(it, 3, gasx, igb, igr, igf, igy, sloc);
      if (sloc2 >= 0) {(void) ImplicitNCDump(it, 5, gasx, igb, igr, igf, igy, sloc2);}
    }
    // (c) the emission/absorption source, linearised in T about the iterate
    if (src_on) {
      auto eos = flr.eos;
      // m1-fast5-sp: an ideal-gas EOS without the cache gets its own kernel (e and c_v by
      // M1EosIdeal, the ideal branch of ThermoAt: bitwise); the generic one carries the
      // table code (828 B call frame per lane) even when the table is off
      const bool srid = !eos.tbl.active && !usec;
      const Real gam = eos.gamma;
      // nvcc forbids generic (auto) extended lambdas, so the tag-dependent
      // helper and its kernel are the function template M1ImplSrcLaunch
      // (above); srb_ctx is what the helper captured, by value.
      auto srb_ctx = std::make_tuple(ar, ch, cl, dt, ec_, ecnt, eos, gam, gasx, gnw,
                                     igb, igf, igm, igr, igy, iw_, opac_, uh, usec, nmb1,
                                     ks, ke, js, je, is, ie, impl_src_stable);
      if (srid) {
        M1ImplSrcLaunch(srb_ctx, std::true_type{});
      } else {
        M1ImplSrcLaunch(srb_ctx, std::false_type{});
      }
    }

    // implicit_vimp: the Jacobian of the implicit enthalpy velocity for this pass
    if (vimp_now) {ImplicitVimpBuild();}
    if (muscl_now) {ImplicitMusclBuild();}   // xthinfix-1009 Fix B
    const bool vim = vimp_now;
    const int ivb = iw_vimp;

    // STAGE S1 (rad_m1_sph.cpp): the spherical-polar geometry of the x1 faces
    const bool sph = sph_geom;
    auto cvol = pmy_pack->pcoord->volume;
    auto carea = pmy_pack->pcoord->area;
    auto cdxf = pmy_pack->pcoord->dxface;
    // STAGE S2: a chi(f) closure on the wedge -- the radial integrating factor
    // (M1SphDrr) and the lagged curvature (M1SphCurv), inside the sp overwrites only
    const bool sphq = sph_q;
    const bool odl = (odm != M1_OD_NONE);
    // vet_col_lat (m1-vetcol-lat): the lagged lateral off-diagonal term M1SphLat
    const bool vlat = vlat_now && vlat_ready && sph;
    auto vlt_ = tau_ten;
    const int c0l = M1_TT_LAT0;
    auto cx1v = pmy_pack->pcoord->x1v;
    auto cx2v = pmy_pack->pcoord->x2v;
    auto cx3v = pmy_pack->pcoord->x3v;
    auto cx1f = pmy_pack->pcoord->xx1f;
    // implicit_face_weight = distance (sp): x1 face kappa_T, v1, G0 by M1FaceAvgX1
    const bool fwd = sph && impl_face_wdist;
    // implicit_marshak_face = linear (sp): the Marshak faces take the face E
    const bool mfl = impl_mface_lin;
    const bool sqf = vcol_sqf;   // vet_col_surface_face (rad_m1.hpp)
    // m1-sp-order2b (time2_vstage): the gas WORK v.dm of the radiative kick in the E row
    // (lagged at the iterate's kick dv^k of this pass) instead of only after the solve.
    // Subtracted after the solve it is an O(dt) splitting in every stage (passive E in a
    // moving scattering gas: E order 1.1 in time); in the row the stage value satisfies
    // its own equation.  The write-back removes only work - work^k (0 at convergence).
    // force_reference_work = split: the work goes into the row as well (fws)
    const bool wimp = t2st && vim && have_hydro && feedback &&
                      ((t2_fvnew && !fref) || fws);
    t2_wk = wimp;
    // (d) assemble the tridiagonal system of every column
    // implicit_enthalpy: the deferred correction of the x1 enthalpy flux (header)
    const int enm = impl_enth;
    const bool enth2 = (enm != M1_IENTH_UPWIND);
    const int ngh = indcs.ng;
    // implicit_opac_newton: the Newton term of the flux opacity in the Cartesian x1 rows.
    // A face term of the row, G = nu cr om th (f0n - ...) + nu df (w_c E_c - w_n E_n),
    // is proportional to th = 1/(1 + chat dt ktf), so dG/dktf = -chat dt th G.  The
    // face opacity is the mean of the two cells', and a cell's T moves in the local gas
    // Newton step by dT = (R_k + c dt rho kappa_E E')/B_k (m1_impl_src/tsolve), so
    //   G(E', kappa(T_new)) ~ G(E', kappa_k) + q sum_c ktd_c (R_c + c dt kE_c E'_c)/B_c,
    // q = -chat dt th G(E^k)/2.  The E'_c part goes on the row, the rest to the RHS.
    // Only faces whose two cells are in this block (x1 neighbours of the block are
    // not exchanged for ktd): there the term is omitted and the face stays Picard.
    // On sp the rows are rebuilt below with the sp areas; the same term goes there
    // (opns).
    const bool opnr = opn && gnw && !sph_geom;
    const bool opns = opn && gnw && sph_geom;
    // DEBUG dbg_t2_admiss: the sp E row by term, kept for T2AdmissDebug
    const bool dbgr = (t2_dbg_adm > 0) && sph_geom;
    if (dbgr && dbrow.extent_int(0) < nmb1 + 1) {
      Kokkos::realloc(dbrow, nmb1 + 1, 16, iw.extent_int(2), iw.extent_int(3),
                      iw.extent_int(4));
      Kokkos::deep_copy(dbrow, 0.0);
    }
    auto dbrow_ = dbrow;
    auto ktdv = ktd;
    // implicit_opac_newton_guard: record each face's Newton term (gp* face i+1/2, gm*
    // face i-1/2: diagonal, neighbour entry, rr part) and take it out at the end of the
    // row if it leaves the diagonal below opg x its value without it
    const Real opg = impl_opn_guard;
    const int ogm = impl_opn_guard_mode;
    const bool opgd = (opg > 0.0) && (opnr || opns);
    auto nsk_ = opn_nskip_d;
    par_for_lb("m1_impl_asm", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) M1_INL {
      Real dx = mbsize.d_view(m).dx1;
      Real nu = dt/dx;
      Real cr = ch/cl;
      int ipos = pos_.d_view(m);
      bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
      Real wi = iw_(m,M1_IW_WCHI,k,j,i);
      // in 1-D the stored closure IS the diagonal Eddington component (P_11/E = chi);
      // in multi-D it is chi and D_11 has to be built from the lagged flux direction.
      if (trans) {wi = M1DDiag(iw_,vd_,dfull,m,0,k,j,i);}
      Real ai = iw_(m,M1_IW_ADV,k,j,i);
      Real vi = iw_(m,M1_IW_V1,k,j,i);
      Real aa = 0.0, bb = 1.0, cc = 0.0;
      Real rr = iw_(m,M1_IW_EN,k,j,i) + iw_(m,M1_IW_SRCR,k,j,i);
      bb += iw_(m,M1_IW_SRCB,k,j,i);
      Real gpd = 0.0, gpo = 0.0, gpr = 0.0, gmd = 0.0, gmo = 0.0, gmr = 0.0;
      // MILESTONE 3b phase B: the LINE-JACOBI transverse couplings.  Their diagonal part
      // stays on the diagonal (the 7-point M-matrix), the neighbours' lagged part goes to
      // the right-hand side.
      int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
      if (trans) {
        bb += iw_(m,M1_IW_TDIA,k,j,i);
        rr += iw_(m,M1_IW_TRHS,k,j,i);
        BoundaryFlag q1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
        BoundaryFlag q2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
        BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
        BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
        BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
        BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
        if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
        if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
        if ((q3 == BoundaryFlag::block) || (q3 == BoundaryFlag::periodic)) {jl = js-1;}
        if ((q4 == BoundaryFlag::block) || (q4 == BoundaryFlag::periodic)) {ju = je+1;}
        if (thrd && ((q5 == BoundaryFlag::block) ||
                     (q5 == BoundaryFlag::periodic))) {kl = ks-1;}
        if (thrd && ((q6 == BoundaryFlag::block) ||
                     (q6 == BoundaryFlag::periodic))) {ku = ke+1;}
      }
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      // implicit_enthalpy: the x1 ghost layers the plm stencil may read on each side
      int hxl = 0, hxh = 0;
      if (enth2) {
        if (trans) {
          hxl = (il < is) ? ngh : 0;
          hxh = (iu > ie) ? ngh : 0;
        } else {
          hxl = botb ? 0 : nlay_;
          hxh = topb ? 0 : nlay_;
        }
      }

      // ---- face i+1/2
      if (i < ie || cyclic || !topb) {
        int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i+1);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        Real wp = iw_(m,M1_IW_WCHI,k,j,ip);
        if (trans) {wp = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);}
        cc -= nu*df*wp;
        Real vf = 0.5*(vi + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,i) + iw_(m,M1_IW_G0,k,j,ip));
        // the OFF-diagonal Eddington terms of the x1 flux equation, d_2 P_12 + d_3 P_13,
        // fully lagged and therefore a right-hand-side term
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,i,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                             dfull)
                    + M1OffDiv(iw_,m,0,k,j,ip,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,
                               M1_IW_EP,vd_,dfull));
        }
        rr -= nu*cr*om*th*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f - ch*cl*dt*od);
        if (opnr && ip >= is && ip <= ie) {
          const Real gf = nu*cr*om*th*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f - ch*cl*dt*od)
                          + nu*df*(wi*iw_(m,M1_IW_EP,k,j,i) - wp*iw_(m,M1_IW_EP,k,j,ip));
          const Real q = -0.5*ch*dt*th*gf;
          if (opgd) {
            M1OpnCellRec(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                         opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr, gpd, gpr);
            M1OpnCellRec(q, ktdv(m,k,j,ip), iw_(m,igb,k,j,ip), iw_(m,igr,k,j,ip),
                         opac_(m,M1_OP_E,k,j,ip), cl, dt, cc, rr, gpo, gpr);
          } else {
            M1OpnCell(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                      opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr);
            M1OpnCell(q, ktdv(m,k,j,ip), iw_(m,igb,k,j,ip), iw_(m,igr,k,j,ip),
                      opac_(m,M1_OP_E,k,j,ip), cl, dt, cc, rr);
          }
        }
        // the HLL part: its E'_L coefficient is >= 0 (diagonal) and its E'_R coefficient
        // <= 0 (upper off-diagonal), so the blend keeps the M-matrix.
        bb += nu*ifw_(m,M1_IFW_HCL,k,j,i+1);
        cc += nu*ifw_(m,M1_IFW_HCR,k,j,i+1);
        rr -= nu*ifw_(m,M1_IFW_DG,k,j,i+1);
        if (vf > 0.0) {
          bb += nu*cr*ai;
        } else {
          cc += nu*cr*iw_(m,M1_IW_ADV,k,j,ip);
        }
        if (enth2) {
          bool o0, o3;
          int i0 = M1EnthIdx(i-1, is, ie, cyclic, hxl, hxh, o0);
          int i3 = M1EnthIdx(i+2, is, ie, cyclic, hxl, hxh, o3);
          if (t2vs) {
            rr -= nu*cr*M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,i),
                                    iw_(m,M1_IW_EP,k,j,ip), iw_(m,M1_IW_EP,k,j,i3),
                                    o0 && o3, iw_(m,M1_IW_ADV,k,j,i0), ai,
                                    iw_(m,M1_IW_ADV,k,j,ip), iw_(m,M1_IW_ADV,k,j,i3),
                                    iw_(m,t2da,k,j,i0), iw_(m,t2da,k,j,i),
                                    iw_(m,t2da,k,j,ip), iw_(m,t2da,k,j,i3), vf, t2afc);
          } else {
          rr -= nu*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,i),
                                 iw_(m,M1_IW_EP,k,j,ip), iw_(m,M1_IW_EP,k,j,i3),
                                 o0 && o3, iw_(m,M1_IW_ADV,k,j,i0), ai,
                                 iw_(m,M1_IW_ADV,k,j,ip), iw_(m,M1_IW_ADV,k,j,i3), vf);
          }
        }
      } else if (bchi == M1_IBC_MARSHAK) {
        const Real mq = vqs ? vq_(m,k,j) : mqo;   // vet_col_surface_q
        bb += nu*ch*mq;
        rr += nu*ch*mq*ebhi;
        // implicit_bc_advect: the enthalpy flux A E through the end face, upwinded with
        // the end cell's velocity (outflow: this cell's E; inflow: the bath)
        if (badv) {
          if (vi > 0.0) {
            bb += nu*cr*ai;
          } else {
            rr -= nu*cr*ai*ebhi;
          }
        }
      } else if (bchi == M1_IBC_FLUX) {
        rr -= nu*cr*fxhi;
      }

      // ---- face i-1/2
      if (i > is || cyclic || !botb) {
        int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        Real wm = iw_(m,M1_IW_WCHI,k,j,im);
        if (trans) {wm = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);}
        aa -= nu*df*wm;
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + vi);
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,i));
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                             vd_,dfull)
                    + M1OffDiv(iw_,m,0,k,j,i,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                               vd_,dfull));
        }
        rr += nu*cr*om*th*(f0n_(m,k,j,i) - ch*dt*vf*g0f - ch*cl*dt*od);
        if (opnr && im >= is && im <= ie) {
          const Real gf = -nu*cr*om*th*(f0n_(m,k,j,i) - ch*dt*vf*g0f - ch*cl*dt*od)
                          + nu*df*(wi*iw_(m,M1_IW_EP,k,j,i) - wm*iw_(m,M1_IW_EP,k,j,im));
          const Real q = -0.5*ch*dt*th*gf;
          if (opgd) {
            M1OpnCellRec(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                         opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr, gmd, gmr);
            M1OpnCellRec(q, ktdv(m,k,j,im), iw_(m,igb,k,j,im), iw_(m,igr,k,j,im),
                         opac_(m,M1_OP_E,k,j,im), cl, dt, aa, rr, gmo, gmr);
          } else {
            M1OpnCell(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                      opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr);
            M1OpnCell(q, ktdv(m,k,j,im), iw_(m,igb,k,j,im), iw_(m,igr,k,j,im),
                      opac_(m,M1_OP_E,k,j,im), cl, dt, aa, rr);
          }
        }
        aa -= nu*ifw_(m,M1_IFW_HCL,k,j,i);
        bb -= nu*ifw_(m,M1_IFW_HCR,k,j,i);
        rr += nu*ifw_(m,M1_IFW_DG,k,j,i);
        if (vf > 0.0) {
          aa -= nu*cr*iw_(m,M1_IW_ADV,k,j,im);
        } else {
          bb -= nu*cr*ai;
        }
        if (enth2) {
          bool o0, o3;
          int i0 = M1EnthIdx(i-2, is, ie, cyclic, hxl, hxh, o0);
          int i3 = M1EnthIdx(i+1, is, ie, cyclic, hxl, hxh, o3);
          if (t2vs) {
            rr += nu*cr*M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,im),
                                    iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EP,k,j,i3),
                                    o0 && o3, iw_(m,M1_IW_ADV,k,j,i0),
                                    iw_(m,M1_IW_ADV,k,j,im), ai,
                                    iw_(m,M1_IW_ADV,k,j,i3), iw_(m,t2da,k,j,i0),
                                    iw_(m,t2da,k,j,im), iw_(m,t2da,k,j,i),
                                    iw_(m,t2da,k,j,i3), vf, t2afc);
          } else {
          rr += nu*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,im),
                                 iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EP,k,j,i3),
                                 o0 && o3, iw_(m,M1_IW_ADV,k,j,i0),
                                 iw_(m,M1_IW_ADV,k,j,im), ai, iw_(m,M1_IW_ADV,k,j,i3),
                                 vf);
          }
        }
      } else if (bclo == M1_IBC_MARSHAK) {
        bb += nu*ch*mq;
        rr += nu*ch*mq*eblo;
        if (badv) {
          if (vi > 0.0) {
            rr += nu*cr*ai*eblo;
          } else {
            bb -= nu*cr*ai;
          }
        }
      } else if (bclo == M1_IBC_FLUX) {
        rr += nu*cr*fxlo;
      }

      // implicit_vimp: the diagonal and x1 +-1 part of the implicit enthalpy velocity
      // goes into the row (and so into the line preconditioner), its lagged part and
      // J E^k into the right-hand side; the rest is applied by M1VimpRow
      if (vim) {
        aa += iw_(m,ivb+M1_IV_J1M,k,j,i);
        bb += iw_(m,ivb+M1_IV_JD,k,j,i);
        cc += iw_(m,ivb+M1_IV_J1P,k,j,i);
        rr += iw_(m,ivb+M1_IV_JRHS,k,j,i);
      }

      // STAGE S1 (sp): the row rebuilt with dt A_f/V_i per face and the centroid
      // distance dxface.x1f in the face-flux gradient.  The Cartesian row above is left
      // textually untouched and overwritten here.  On sp (SphericalS1Check) the flux is
      // central (om = 1, no HLL/DG part; sp-blend-1008: or berthon | blend, the x1
      // HLL part with dt A_f/V_i) and offdiag none or lagged; m1-sph2 adds the
      // hesdirk2 stage solves (the old vector is generic) and implicit_vimp.
      if (sph) {
        Real iv = dt/cvol(m,k,j,i);
        Real nup = carea.x1f(m,k,j,i+1)*iv, num = carea.x1f(m,k,j,i)*iv;
        aa = 0.0;
        bb = 1.0 + iw_(m,M1_IW_SRCB,k,j,i);
        cc = 0.0;
        rr = iw_(m,M1_IW_EN,k,j,i) + iw_(m,M1_IW_SRCR,k,j,i);
        if (trans) {
          bb += iw_(m,M1_IW_TDIA,k,j,i);
          rr += iw_(m,M1_IW_TRHS,k,j,i);
        }
        // DEBUG dbg_t2_admiss: the row by term (d_[0..6] of bb, d_[7..15] of rr)
        Real d_[16] = {0.0};
        Real pb_ = bb, pr_ = rr;
        if (i < ie || !topb) {
          int ip = (i < ie) ? (i+1) : (ie+1);
          Real ktf = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,i), iw_(m,M1_IW_KT,k,j,ip), cx1f, m,
                                 i, ip, fwd);
          Real th = 1.0/(1.0 + ch*dt*ktf);
          // sp-blend-1008: the central (face-eliminated) part of the face flux keeps the
          // weight 1 - AL of the blend (exactly th when implicit_flux = central)
          Real tw = th;
          if (aphll) {tw = (1.0 - ifw_(m,M1_IFW_AL,k,j,i+1))*th;}
          Real df = tw*ch*ch*dt/cdxf.x1f(m,k,j,i+1);
          Real wp = iw_(m,M1_IW_WCHI,k,j,ip);
          if (trans) {wp = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);}
          Real wiu = wi;
          Real ods = 0.0;   // the S2 curvature term (implicit_opac_newton needs it)
          if (sphq) {
            // S2: the integrating factor, (r_c/r_f)^2 on the q n_r^2 part of each cell
            Real rf = cx1f(m,i+1);
            Real si = SQR(cx1v(m,i)/rf), sp = SQR(cx1v(m,ip)/rf);
            Real n1i = trans ? iw_(m,M1_IW_N1,k,j,i) : 1.0;
            Real n1p = trans ? iw_(m,M1_IW_N1,k,j,ip) : 1.0;
            wiu = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,i), n1i, si);
            wp = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,ip), n1p, sp);
            if (trans) {
              Real od = 0.5*(M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,i,odl,thrd,il,iu,jl,ju,
                                       kl,ku,M1_IW_EP)
                             + M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,ip,odl,thrd,il,iu,jl,
                                         ju,kl,ku,M1_IW_EP));
              if (vlat) {
                od += 0.5*(M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,i,thrd,il,iu,jl,
                                       ju,kl,ku,M1_IW_EP)
                           + M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,ip,thrd,
                                      il,iu,jl,ju,kl,ku,M1_IW_EP));
              }
              rr += nup*cr*tw*ch*cl*dt*od;
              d_[10] += nup*cr*tw*ch*cl*dt*od;
              ods = od;
            }
          }
          bb += nup*df*wiu;
          d_[0] += nup*df*wiu;
          cc -= nup*df*wp;
          Real vf = M1FaceAvgX1(vi, iw_(m,M1_IW_V1,k,j,ip), cx1f, m, i, ip, fwd);
          Real g0f = M1FaceAvgX1(iw_(m,M1_IW_G0,k,j,i), iw_(m,M1_IW_G0,k,j,ip), cx1f, m,
                                 i, ip, fwd);
          rr -= nup*cr*tw*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f);
          d_[11] -= nup*cr*tw*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f);
          // implicit_opac_newton on sp: the same face term G (every part of it is
          // proportional to th), with the sp area factor and the S2 row coefficients
          if (opns && ip <= ie) {
            const Real gf = nup*cr*tw*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f - ch*cl*dt*ods)
                            + nup*df*(wiu*iw_(m,M1_IW_EP,k,j,i)
                                      - wp*iw_(m,M1_IW_EP,k,j,ip));
            const Real q = -0.5*ch*dt*th*gf;
            if (opgd) {
              const Real b0 = bb, r0 = rr;
              M1OpnCellRec(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                           opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr, gpd, gpr);
              M1OpnCellRec(q, ktdv(m,k,j,ip), iw_(m,igb,k,j,ip), iw_(m,igr,k,j,ip),
                           opac_(m,M1_OP_E,k,j,ip), cl, dt, cc, rr, gpo, gpr);
              d_[1] += bb - b0; d_[12] += rr - r0;
            } else {
              {const Real b0 = bb, r0 = rr;
              M1OpnCell(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                        opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr);
              d_[1] += bb - b0; d_[12] += rr - r0;}
              {const Real r0 = rr;
              M1OpnCell(q, ktdv(m,k,j,ip), iw_(m,igb,k,j,ip), iw_(m,igr,k,j,ip),
                        opac_(m,M1_OP_E,k,j,ip), cl, dt, cc, rr);
              d_[12] += rr - r0;}
            }
          }
          // sp-blend-1008: the berthon part of the face flux, ifw HCL E_i + HCR E_ip + DG
          // (HCL >= 0 on the diagonal, HCR <= 0 off it), times dt A_f/V_i
          if (aphll) {
            bb += nup*ifw_(m,M1_IFW_HCL,k,j,i+1);
            cc += nup*ifw_(m,M1_IFW_HCR,k,j,i+1);
            rr -= nup*ifw_(m,M1_IFW_DG,k,j,i+1);
          }
          if (vf > 0.0) {
            bb += nup*cr*ai;
            d_[2] += nup*cr*ai;
          } else {
            cc += nup*cr*iw_(m,M1_IW_ADV,k,j,ip);
          }
          const Real re0_ = rr;
          if (enth2) {
            bool o0, o3;
            int i0 = M1EnthIdx(i-1, is, ie, false, hxl, hxh, o0);
            int i3 = M1EnthIdx(i+2, is, ie, false, hxl, hxh, o3);
            if (t2vs) {
              // m1-sph2: a hesdirk2 stage solve under implicit_vimp, as the Cartesian row
              rr -= nup*cr*M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j,i0),
                                       iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EP,k,j,ip),
                                       iw_(m,M1_IW_EP,k,j,i3), o0 && o3,
                                       iw_(m,M1_IW_ADV,k,j,i0), ai,
                                       iw_(m,M1_IW_ADV,k,j,ip), iw_(m,M1_IW_ADV,k,j,i3),
                                       iw_(m,t2da,k,j,i0), iw_(m,t2da,k,j,i),
                                       iw_(m,t2da,k,j,ip), iw_(m,t2da,k,j,i3), vf, t2afc);
            } else {
            rr -= nup*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,i),
                                    iw_(m,M1_IW_EP,k,j,ip), iw_(m,M1_IW_EP,k,j,i3),
                                    o0 && o3, iw_(m,M1_IW_ADV,k,j,i0), ai,
                                    iw_(m,M1_IW_ADV,k,j,ip), iw_(m,M1_IW_ADV,k,j,i3), vf);
            }
          }
          d_[13] += rr - re0_;
        } else if (bchi == M1_IBC_MARSHAK) {
          const Real mq = vqs ? vq_(m,k,j) : mqo;   // vet_col_surface_q
          bb += nup*ch*mq;
          rr += nup*ch*mq*ebhi;
          if (mfl && (!vqs || sqf) && ie > is) {
            // implicit_marshak_face = linear: c q (E_f - E_bath) with E_f the face E.
            // Not with vet_col_surface_q: its q = H(face)/J(top cell) already makes the
            // face flux the formal solution's H with the CELL E (runs_5e); the face J it
            // would need instead is first order (the grazing layer at the top face)
            // The linear face E is implicit (ca, -cb in the row), the limiter remainder a
            // deferred correction.
            Real ca, cb;
            M1SphMarshakCoef(cx1v(m,ie), cx1v(m,ie-1), cx1f(m,ie+1), ca, cb);
            const Real e0 = iw_(m,M1_IW_EP,k,j,ie), e1 = iw_(m,M1_IW_EP,k,j,ie-1);
            const Real dl = M1SphMarshakFaceE(e0, e1, cx1v(m,ie), cx1v(m,ie-1),
                                              cx1f(m,ie+1)) - (ca*e0 - cb*e1);
            // the face flux c q E_f, and the outflowing enthalpy flux of
            // implicit_bc_advect at the same face E
            Real wq = nup*ch*mq;
            if (badv && vi > 0.0) {wq += nup*cr*ai;}
            bb += wq*(ca - 1.0);
            aa -= wq*cb;
            rr -= wq*dl;
          }
          if (badv) {
            if (vi > 0.0) {
              bb += nup*cr*ai;
            } else {
              rr -= nup*cr*ai*ebhi;
            }
          }
        } else if (bchi == M1_IBC_FLUX) {
          rr -= nup*cr*fxhi;
        }
        if (i > is || !botb) {
          int im = (i > is) ? (i-1) : (is-1);
          Real ktf = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,im), iw_(m,M1_IW_KT,k,j,i), cx1f, m,
                                 im, i, fwd);
          Real th = 1.0/(1.0 + ch*dt*ktf);
          Real tw = th;   // sp-blend-1008, as at face i+1/2
          if (aphll) {tw = (1.0 - ifw_(m,M1_IFW_AL,k,j,i))*th;}
          Real df = tw*ch*ch*dt/cdxf.x1f(m,k,j,i);
          Real wm = iw_(m,M1_IW_WCHI,k,j,im);
          if (trans) {wm = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);}
          Real wil = wi;
          Real ods = 0.0;
          if (sphq) {
            Real rf = cx1f(m,i);
            Real si = SQR(cx1v(m,i)/rf), sm = SQR(cx1v(m,im)/rf);
            Real n1i = trans ? iw_(m,M1_IW_N1,k,j,i) : 1.0;
            Real n1m = trans ? iw_(m,M1_IW_N1,k,j,im) : 1.0;
            wil = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,i), n1i, si);
            wm = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,im), n1m, sm);
            if (trans) {
              Real od = 0.5*(M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,im,odl,thrd,il,iu,jl,ju,
                                       kl,ku,M1_IW_EP)
                             + M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,i,odl,thrd,il,iu,jl,
                                         ju,kl,ku,M1_IW_EP));
              if (vlat) {
                od += 0.5*(M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,im,thrd,il,iu,jl,
                                       ju,kl,ku,M1_IW_EP)
                           + M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,i,thrd,il,iu,jl,
                                      ju,kl,ku,M1_IW_EP));
              }
              rr -= num*cr*tw*ch*cl*dt*od;
              d_[10] -= num*cr*tw*ch*cl*dt*od;
              ods = od;
            }
          }
          bb += num*df*wil;
          d_[0] += num*df*wil;
          aa -= num*df*wm;
          Real vf = M1FaceAvgX1(iw_(m,M1_IW_V1,k,j,im), vi, cx1f, m, im, i, fwd);
          Real g0f = M1FaceAvgX1(iw_(m,M1_IW_G0,k,j,im), iw_(m,M1_IW_G0,k,j,i), cx1f, m,
                                 im, i, fwd);
          rr += num*cr*tw*(f0n_(m,k,j,i) - ch*dt*vf*g0f);
          d_[11] += num*cr*tw*(f0n_(m,k,j,i) - ch*dt*vf*g0f);
          if (opns && im >= is) {
            const Real gf = -num*cr*tw*(f0n_(m,k,j,i) - ch*dt*vf*g0f - ch*cl*dt*ods)
                            + num*df*(wil*iw_(m,M1_IW_EP,k,j,i)
                                      - wm*iw_(m,M1_IW_EP,k,j,im));
            const Real q = -0.5*ch*dt*th*gf;
            if (opgd) {
              const Real b0 = bb, r0 = rr;
              M1OpnCellRec(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                           opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr, gmd, gmr);
              M1OpnCellRec(q, ktdv(m,k,j,im), iw_(m,igb,k,j,im), iw_(m,igr,k,j,im),
                           opac_(m,M1_OP_E,k,j,im), cl, dt, aa, rr, gmo, gmr);
              d_[1] += bb - b0; d_[12] += rr - r0;
            } else {
              {const Real b0 = bb, r0 = rr;
              M1OpnCell(q, ktdv(m,k,j,i), iw_(m,igb,k,j,i), iw_(m,igr,k,j,i),
                        opac_(m,M1_OP_E,k,j,i), cl, dt, bb, rr);
              d_[1] += bb - b0; d_[12] += rr - r0;}
              {const Real r0 = rr;
              M1OpnCell(q, ktdv(m,k,j,im), iw_(m,igb,k,j,im), iw_(m,igr,k,j,im),
                        opac_(m,M1_OP_E,k,j,im), cl, dt, aa, rr);
              d_[12] += rr - r0;}
            }
          }
          if (aphll) {
            aa -= num*ifw_(m,M1_IFW_HCL,k,j,i);
            bb -= num*ifw_(m,M1_IFW_HCR,k,j,i);
            rr += num*ifw_(m,M1_IFW_DG,k,j,i);
          }
          if (vf > 0.0) {
            aa -= num*cr*iw_(m,M1_IW_ADV,k,j,im);
          } else {
            bb -= num*cr*ai;
            d_[2] -= num*cr*ai;
          }
          const Real re1_ = rr;
          if (enth2) {
            bool o0, o3;
            int i0 = M1EnthIdx(i-2, is, ie, false, hxl, hxh, o0);
            int i3 = M1EnthIdx(i+1, is, ie, false, hxl, hxh, o3);
            if (t2vs) {
              rr += num*cr*M1EnthCorrT(enm, iw_(m,M1_IW_EP,k,j,i0),
                                       iw_(m,M1_IW_EP,k,j,im), iw_(m,M1_IW_EP,k,j,i),
                                       iw_(m,M1_IW_EP,k,j,i3), o0 && o3,
                                       iw_(m,M1_IW_ADV,k,j,i0), iw_(m,M1_IW_ADV,k,j,im),
                                       ai, iw_(m,M1_IW_ADV,k,j,i3), iw_(m,t2da,k,j,i0),
                                       iw_(m,t2da,k,j,im), iw_(m,t2da,k,j,i),
                                       iw_(m,t2da,k,j,i3), vf, t2afc);
            } else {
            rr += num*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,im),
                                    iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EP,k,j,i3),
                                    o0 && o3, iw_(m,M1_IW_ADV,k,j,i0),
                                    iw_(m,M1_IW_ADV,k,j,im), ai,
                                    iw_(m,M1_IW_ADV,k,j,i3), vf);
            }
          }
          d_[13] += rr - re1_;
        } else if (bclo == M1_IBC_MARSHAK) {
          bb += num*ch*mq;
          rr += num*ch*mq*eblo;
          if (mfl && ie > is) {
            Real ca, cb;
            M1SphMarshakCoef(cx1v(m,is), cx1v(m,is+1), cx1f(m,is), ca, cb);
            const Real e0 = iw_(m,M1_IW_EP,k,j,is), e1 = iw_(m,M1_IW_EP,k,j,is+1);
            const Real dl = M1SphMarshakFaceE(e0, e1, cx1v(m,is), cx1v(m,is+1),
                                              cx1f(m,is)) - (ca*e0 - cb*e1);
            Real wq = num*ch*mq;
            if (badv && vi < 0.0) {wq -= num*cr*ai;}
            bb += wq*(ca - 1.0);
            cc -= wq*cb;
            rr -= wq*dl;
          }
          if (badv) {
            if (vi > 0.0) {
              rr += num*cr*ai*eblo;
            } else {
              bb -= num*cr*ai;
            }
          }
        } else if (bclo == M1_IBC_FLUX) {
          rr += num*cr*fxlo;
        }
        // m1-sph2: implicit_vimp, as in the Cartesian row above (ImplicitVimpBuild
        // builds its coefficients with the sp areas, volumes and face distances)
        if (vim) {
          aa += iw_(m,ivb+M1_IV_J1M,k,j,i);
          bb += iw_(m,ivb+M1_IV_JD,k,j,i);
          cc += iw_(m,ivb+M1_IV_J1P,k,j,i);
          rr += iw_(m,ivb+M1_IV_JRHS,k,j,i);
        }
        if (dbgr) {
          if (vim) {
            d_[4] = iw_(m,ivb+M1_IV_JD,k,j,i);
            d_[14] = iw_(m,ivb+M1_IV_JRHS,k,j,i);
          }
          d_[5] = pb_;
          d_[6] = bb - pb_ - d_[0] - d_[1] - d_[2] - d_[4];
          d_[7] = pr_;
          d_[15] = rr - pr_ - d_[10] - d_[11] - d_[12] - d_[13] - d_[14];
          for (int q = 0; q < 16; ++q) {dbrow_(m,q,k,j,i) = d_[q];}
        }
      }

      // implicit_opac_newton_guard: face i+1/2 first, then i-1/2 against the diagonal
      // left by the first (the two can together take it no lower than opg^2 x b0)
      if (opgd) {
        Real bc = bb - gpd - gmd, rc = rr + gpr + gmr;
        int ns = M1OpnGuardFace(opg, ogm, gpd, gpo, gpr, bc, rc, bb, cc, rr);
        const int nsp = ns;
        ns += M1OpnGuardFace(opg, ogm, gmd, gmo, gmr, bc, rc, bb, aa, rr);
        if ((ogm & 4) && ns < 2) {
          // m1-positivity, guard_mode bit 4: the ROW DIAGONAL DOMINANCE.  A face's
          // Newton term may leave every entry with the right sign and still make the
          // neighbour entry larger than the diagonal can carry (He presn wedge plume and
          // photosphere: sum|off|/diag = 2.4 with no positive entry, solved E < 0).  A
          // Z-matrix whose rows are strictly dominant is an M-matrix, so the Newton terms
          // still in the row are taken out (the faces stay Picard: same fixed point) when
          // they grow sum|off|/diag well beyond what the row had without them.
          const bool kp = (nsp == 0), km = (ns - nsp == 0);
          Real toff = 0.0;
          if (bicg && trans) {
            toff = fabs(iw_(m,M1_IW_CJM,k,j,i)) + fabs(iw_(m,M1_IW_CJP,k,j,i))
                   + fabs(iw_(m,M1_IW_CKM,k,j,i)) + fabs(iw_(m,M1_IW_CKP,k,j,i));
          }
          const Real b0 = bb - (kp ? gpd : 0.0) - (km ? gmd : 0.0);
          const Real a0 = aa - (km ? gmo : 0.0), c0 = cc - (kp ? gpo : 0.0);
          // the ratio sum|off|/diag with and without the terms: the rows of a smooth
          // diffusion operator sit at ~1 (column sums, not rows, carry conservation), so
          // the test is the growth of that ratio beyond max(1, its value without the
          // terms) x (1 + guard)
          const Real rt0 = (fabs(a0) + fabs(c0) + toff)/fmax(b0, 1.0e-300);
          const Real rt1 = (fabs(aa) + fabs(cc) + toff)/fmax(bb, 1.0e-300);
          if (!(bb > 0.0) || rt1 > fmax(rt0, 1.0)*(1.0 + opg)) {
            bb = b0;
            aa = a0;
            cc = c0;
            rr += (kp ? gpr : 0.0) + (km ? gmr : 0.0);
            ns += (kp ? 1 : 0) + (km ? 1 : 0);
          }
        }
        if (ns > 0) {Kokkos::atomic_add(&nsk_(0), static_cast<Real>(ns));}
      }
      // a Dirichlet end cell: the whole row is replaced, which keeps the matrix an
      // M-matrix and anchors the level of E (see M1_IBC_EFIX)
      if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                      (i == ie && topb && bchi == M1_IBC_EFIX))) {
        aa = 0.0;
        bb = 1.0;
        cc = 0.0;
        rr = iw_(m,M1_IW_EN,k,j,i);
        // the WHOLE row is replaced, transverse couplings included, or the BiCGStab
        // operator would carry an off-diagonal the preconditioner's row does not have.
        if (bicg) {
          iw_(m,M1_IW_CJM,k,j,i) = 0.0;
          iw_(m,M1_IW_CJP,k,j,i) = 0.0;
          iw_(m,M1_IW_CKM,k,j,i) = 0.0;
          iw_(m,M1_IW_CKP,k,j,i) = 0.0;
        }
      }
      iw_(m,M1_IW_TA,k,j,i) = aa;
      iw_(m,M1_IW_TB,k,j,i) = bb;
      iw_(m,M1_IW_TC,k,j,i) = cc;
      iw_(m,M1_IW_TR,k,j,i) = rr;
    });
    // time2_vstage: the gas work of the iterate's kick on the right-hand side (its own
    // kernel, so the assembly kernel of every other configuration is untouched)
    if (wimp) {
      ImplicitWorkRow(true);
    }

    // (e) SOLVE the linear system of this pass.  With implicit_solver = line_jacobi
    // that is ONE x1 line solve with the lagged transverse term already on the
    // right-hand side (the outer Picard loop then IS the Jacobi iteration); with
    // bicgstab the very same system -- the same matrix and the same right-hand side --
    // is solved to implicit_lin_tol by a Krylov iteration preconditioned by that line
    // solve, and the outer loop is left with the nonlinearity alone.
    if (bicg) {
      // b = TR + sum_nb C_nb E^k_nb: undo the move of the lagged off-diagonal term to
      // the right-hand side, so that the operator and the right-hand side describe the
      // same system.  EP still holds E^k here, ghosts included.
      par_for("m1_impl_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real b = iw_(m,M1_IW_TR,k,j,i)
                 + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,M1_IW_EP,k,j-1,i)
                 + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,M1_IW_EP,k,j+1,i);
        if (thrd) {
          b += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,M1_IW_EP,k-1,j,i)
               + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
        }
        iw_(m,M1_IW_KB,k,j,i) = b;
      });
      // MILESTONE 3b phase D.  Under implicit_offdiag = operator the assembly has put
      // -L_off(E^k) on the right-hand side (it is inside TR, through the face fluxes);
      // adding L_off(E^k) back takes it out again, and the operator application adds
      // L_off(x) on the LEFT.  The two changes cancel at x = E^k by construction, so the
      // residual the linear solver measures is the residual of the same system the
      // lagged form measures -- what changes is where the term is solved.
      if (odm == M1_OD_OPERATOR && !impl_odskip) {
        ImplicitOffDiagOp(M1_IW_EP, M1_IW_KB, 1.0);
      }
      // vet_col_lat_offdiag = operator: the same for the lateral off-diagonal term (the
      // assembly put -L_lat(E^k) into TR through M1SphLat; the operator adds L_lat(x))
      if (VlatOp()) {VetLatOp(M1_IW_EP, M1_IW_KB, 1.0);}
      // implicit_hr_damp = beta > 0 (xthinfix-1009, odCMFD-like CONSISTENT damping of the
      // Picard loop under implicit_hr_recon = plm): beta TB (E' - E^k) on the left, i.e.
      // + beta TB on the diagonal of the plm row and + beta TB E^k on the right-hand
      // side.  Zero at the converged fixed point; damps the Picard 2-cycle at mixed plm
      // faces (cyl, XTHINFIX.md) like an under-relaxation 1/(1 + beta).
      if (muscl_now && impl_muscl_damp > 0.0) {
        const Real bdm = impl_muscl_damp;
        const int imd = iw_muscl + M1_IM_D;
        par_for("m1_muscl_damp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const Real d = bdm*iw_(m,M1_IW_TB,k,j,i);
          iw_(m,imd,k,j,i) += d;
          iw_(m,M1_IW_KB,k,j,i) += d*iw_(m,M1_IW_EP,k,j,i);
        });
      }
      kdev_slot = std::min(it, 2);
      TmrMark(5);
      DetTrace("pre_bcg");
      nin = ImplicitBiCGStab(rhsmax);
      DetTrace("post_bcg");
      TmrMark(6);
      if (tmr_on) {
        tmr_cnt[1] += 1.0;
        tmr_cnt[2] += nin;
        if (t2s == M1_T2S_STAGE1) {tmr_cnt[3] += 1.0; tmr_cnt[5] += nin;}
        if (t2s == M1_T2S_STAGE2) {tmr_cnt[4] += 1.0; tmr_cnt[6] += nin;}
      }
      if ((vimp_now || muscl_now) && odm != M1_OD_OPERATOR) {
        // implicit_vimp POSITIVITY: the Newton coupling is not an M-matrix either; a
        // non-positive E drops it for the rest of the step (counted), as for od below.
        // xthinfix-1009 Fix B: the plm half-range row has positive off-diagonals too
        // (the E_c-2 / E_c+2 coefficients of the inflow faces): the same fallback to dc
        Real emin = 1.0e300;
        Kokkos::parallel_reduce("m1_impl_vmmin",
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                               {nmb1+1,ke+1,je+1,ie+1}),
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmin) {
          Real r = iw_(m,M1_IW_S2,k,j,i);
          lmin = (r < lmin) ? r : lmin;
        }, Kokkos::Min<Real>(emin));
#if MPI_PARALLEL_ENABLED
        {Real g;
        MPI_Allreduce(&emin, &g, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
        emin = g;}
#endif
        if (vimp_now) {vimp_emin = std::min(vimp_emin, emin);}
        if (muscl_now) {muscl_emin = std::min(muscl_emin, emin);}
        if (!(emin > 0.0) && muscl_now) {
          // Fix B: dc (sig = 0) in the cells within 2 of a cell solved to E <= 0, for
          // the rest of the step; the next pass rebuilds the row from it (counted)
          ImplicitMusclKill();
          muscl_nfall += 1.0;
        }
        if (!(emin > 0.0) && vimp_now) {
          // DEBUG dbg_t2_admiss (m1-positivity): the rows of the cells this pass solved
          // to E <= 0, WITH the implicit_vimp block, before the block is dropped
          if (t2_dbg_adm_n < t2_dbg_adm) {
            t2_dbg_adm_n += 1;
            if (global_variable::my_rank == 0) {
              std::cout << "<rad_m1> POS_TRIGGER implicit_vimp fallback cycle "
                        << pmy_pack->pmesh->ncycle << " pass " << it << " solve "
                        << (t2st ? "hesdirk2-stage" : "be") << " min E " << emin
                        << std::endl;
            }
            T2AdmissDebug(uh, u0_, t2i_, cl, ch, have_hydro,
                          have_hydro && coupling && dbgh, t2s, it);
          }
          vimp_now = false;
          vimp_nfall += 1.0;
        }
      }
      const bool vlop = VlatOp();
      if (odm == M1_OD_OPERATOR || vlop) {
        // POSITIVITY.  The cross-derivative coefficients have mixed signs, so the
        // 9-/19-point operator is not an M-matrix and E' > 0 is no longer guaranteed.
        // Measure the smallest E the solve produced and, if any cell is non-positive,
        // drop the REST of this step to implicit_offdiag = NONE and count the event.
        // `none` is the M-matrix form that survives: `lagged` is not an alternative
        // here, because lagging these terms in an optically thin cell is exactly what
        // has no fixed point (measured: the seeded He slab blows up in 14 steps with
        // lagged + a frozen closure, and runs with none + a frozen closure).
        Real emin = 1.0e300;
        // m1-perf-0928: a flat range (the MDRange reduction is ~5x slower on CUDA); a
        // min does not depend on the order, so the value is bitwise the same
        const int fni = ie - is + 1, fnj = je - js + 1, fnk = ke - ks + 1;
        Kokkos::parallel_reduce("m1_impl_odmin",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*fnk*fnj*fni),
        KOKKOS_LAMBDA(const int n, Real &lmin) {
          int t = n/fni;
          const int i = is + (n - t*fni);
          const int j = js + (t % fnj);
          t /= fnj;
          const int k = ks + (t % fnk);
          const int m = t/fnk;
          Real r = iw_(m,M1_IW_S2,k,j,i);
          lmin = (r < lmin) ? r : lmin;
        }, Kokkos::Min<Real>(emin));
#if MPI_PARALLEL_ENABLED
        {Real g;
        MPI_Allreduce(&emin, &g, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
        emin = g;}
#endif
        od_emin = std::min(od_emin, emin);
        if (vimp_now) {vimp_emin = std::min(vimp_emin, emin);}
        if (!(emin > 0.0)) {
          if (odm == M1_OD_OPERATOR) {
            od_now = M1_OD_NONE;
            od_nfall += 1.0;
          }
          // vet_col_lat_offdiag = operator: drop D_r,lat for the rest of the step (the
          // assembly's term too: `none`, the M-matrix form), counted
          if (vlop) {
            vlat_now = false;
            vlat_nfall += 1.0;
          }
          if (vimp_now) {
            vimp_now = false;
            vimp_nfall += 1.0;
          }
          if (muscl_now) {
            ImplicitMusclKill();
            muscl_nfall += 1.0;
          }
        }
      }
    } else {
      ImplicitTridiagSolve();
    }

    // (f) accept E', solve for T' and measure the Picard residual
    // F3 implicit_pos_floor_solve: the iterate keeps the solved S2 down to 1e-6 e_floor
    // (no hidden raise to e_floor here), so the face fluxes rebuilt from it are the
    // solve's own and the write-back floor (implicit_pos_floor) raises E to e_floor with
    // the energy charged to the gas or counted
    const Real efls = impl_pos_floor_s2 ? 1.0e-6*efl : efl;
    if (src_on) {
      auto eos = flr.eos;
      // m1-fast5-sp: an ideal-gas EOS without the per-cell cache gets its own kernel,
      // whose root find evaluates e(T) and c_v by the ideal branch of
      // EOS_Data::ThermoAt (the same expressions: bitwise).  The generic kernel carries
      // the table evaluation, 254 VGPRs and 324 B of scratch per lane, even when the
      // table is off (measured, tests_m1/runs_5s_fast5sp).
      const bool tsid = !eos.tbl.active && !usec;
      // nvcc forbids generic (auto) extended lambdas, so the tag-dependent
      // helper and its kernel are the function template M1ImplTsolveLaunch
      // (above); tsb_ctx is what the helper captured, by value.
      auto tsb_ctx = std::make_tuple(ar, cl, dt, ec_, ecnt, efls, eos, escale, gasx, gnw,
                                     igb, igf, igm, igr, igy, iw_, opac_, plog, uh, usec,
                                     nmb1, ks, ke, js, je, is, ie,
                                     (impl_tsolve_opac && (it >= impl_tsolve_opac_start))
                                     || gnsw,
                                     M1OpacFn{opacity_type, otab, kappa_p, kappa_e,
                                              kappa_f, kappa_s, opac_rho_ref,
                                              opac_t_ref, opac_a, opac_b},
                                     impl_tsolve_opac_mode, impl_tsolve_opac_smin,
                                     thin_frz, (impl_thin_frz > 0.0));
      if (tsid) {
        M1ImplTsolveLaunch(tsb_ctx, std::true_type{});
      } else {
        M1ImplTsolveLaunch(tsb_ctx, std::false_type{});
      }
    } else {
      par_for("m1_impl_accept", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efls);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_RES,k,j,i) = fabs(enew - eold)
                                 /fmax(fmax(fabs(enew), escale), 1.0e-300);
      });
    }

    // the new iterate's E has to reach the ghost cells before the FACE update, or the
    // two blocks that share a face would build it from different states.
    if (trans) {
      // E is the ONLY halo quantity the pass changed since the exchange
      // above, so one component goes, not the whole list.
      ImplicitTransverseHalo(1);
    } else {
      ImplicitX1Halo(true);
    }

    // (g) the face fluxes of the new iterate, and the derived cell flux
    const bool musf = muscl_now;   // xthinfix-1009 Fix B: the plm half-range face values
    const int imbf = iw_muscl;
    par_for("m1_impl_face", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      int ipos = pos_.d_view(m);
      bool lo = (i == is) && (ipos == 0);
      bool hi = (i == ie+1) && (ipos == nblkx1-1);
      if ((lo || hi) && !cyclic) {
        int bc = lo ? bclo : bchi;
        int ic = lo ? is : ie;
        Real sgn = lo ? -1.0 : 1.0;
        Real fb = 0.0;
        if (bc == M1_IBC_MARSHAK) {
          const Real mq = (hi && vqs) ? vq_(m,k,j) : mqo;   // vet_col_surface_q
          fb = sgn*cl*mq*(iw_(m,M1_IW_EP,k,j,ic) - (lo ? eblo : ebhi));
          if (mfl && !(hi && vqs && !sqf) && ie > is) {
            // implicit_marshak_face = linear: the face E of the solved iterate
            const int in = lo ? (is+1) : (ie-1);
            const Real ef = M1SphMarshakFaceE(iw_(m,M1_IW_EP,k,j,ic),
                                              iw_(m,M1_IW_EP,k,j,in), cx1v(m,ic),
                                              cx1v(m,in), cx1f(m,i));
            fb = sgn*cl*mq*(ef - (lo ? eblo : ebhi));
          }
        } else if (bc == M1_IBC_FLUX) {
          fb = lo ? fxlo : fxhi;
        } else if (bc == M1_IBC_EFIX) {
          // the Dirichlet row does not define a face flux; the adjacent interior face is
          // copied so that the cell flux and the momentum deposit stay finite.  The
          // boundary cell is held fixed anyway, so nothing downstream depends on it.
          fb = f0_(m,k,j,lo ? (is+1) : ie);
        }
        f0_(m,k,j,i) = fb;
      } else {
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real ktf = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,im), iw_(m,M1_IW_KT,k,j,ip), cx1f, m,
                               im, ip, fwd);
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real vf = M1FaceAvgX1(iw_(m,M1_IW_V1,k,j,im), iw_(m,M1_IW_V1,k,j,ip), cx1f, m,
                              im, ip, fwd);
        Real g0f = M1FaceAvgX1(iw_(m,M1_IW_G0,k,j,im), iw_(m,M1_IW_G0,k,j,ip), cx1f, m,
                               im, ip, fwd);
        Real wp = iw_(m,M1_IW_WCHI,k,j,ip);
        Real wm = iw_(m,M1_IW_WCHI,k,j,im);
        if (trans) {
          wp = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);
          wm = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);
        }
        Real gr = (wp*iw_(m,M1_IW_EP,k,j,ip) - wm*iw_(m,M1_IW_EP,k,j,im))/dx;
        if (sph) {
          gr = (wp*iw_(m,M1_IW_EP,k,j,ip) - wm*iw_(m,M1_IW_EP,k,j,im))/cdxf.x1f(m,k,j,i);
        }
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          Real dx2 = mbsize.d_view(m).dx2;
          Real dx3 = mbsize.d_view(m).dx3;
          int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
          BoundaryFlag q1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
          BoundaryFlag q2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
          BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
          BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
          BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
          BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
          if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
          if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
          if ((q3 == BoundaryFlag::block) || (q3 == BoundaryFlag::periodic)) {jl = js-1;}
          if ((q4 == BoundaryFlag::block) || (q4 == BoundaryFlag::periodic)) {ju = je+1;}
          if (thrd && ((q5 == BoundaryFlag::block) ||
                       (q5 == BoundaryFlag::periodic))) {kl = ks-1;}
          if (thrd && ((q6 == BoundaryFlag::block) ||
                       (q6 == BoundaryFlag::periodic))) {ku = ke+1;}
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                             vd_,dfull)
                    + M1OffDiv(iw_,m,0,k,j,ip,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,
                               M1_IW_EP,vd_,dfull));
        }
        if (sphq) {
          // STAGE S2: the integrating-factor gradient and the lagged curvature, the same
          // expressions the row was assembled with (m1_impl_asm)
          Real rf = cx1f(m,i);
          Real n1p = trans ? iw_(m,M1_IW_N1,k,j,ip) : 1.0;
          Real n1m = trans ? iw_(m,M1_IW_N1,k,j,im) : 1.0;
          Real wps = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,ip), n1p, SQR(cx1v(m,ip)/rf));
          Real wms = M1SphDrr(iw_(m,M1_IW_WCHI,k,j,im), n1m, SQR(cx1v(m,im)/rf));
          gr = (wps*iw_(m,M1_IW_EP,k,j,ip) - wms*iw_(m,M1_IW_EP,k,j,im))
               /cdxf.x1f(m,k,j,i);
          od = 0.0;
          if (trans) {
            int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
            BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
            BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
            BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
            BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
            if ((q3 == BoundaryFlag::block) || (q3 == BoundaryFlag::periodic)) {
              jl = js-1;
            }
            if ((q4 == BoundaryFlag::block) || (q4 == BoundaryFlag::periodic)) {
              ju = je+1;
            }
            if (thrd && ((q5 == BoundaryFlag::block) ||
                         (q5 == BoundaryFlag::periodic))) {kl = ks-1;}
            if (thrd && ((q6 == BoundaryFlag::block) ||
                         (q6 == BoundaryFlag::periodic))) {ku = ke+1;}
            od = 0.5*(M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,im,odl,thrd,il,iu,jl,ju,kl,ku,
                                M1_IW_EP)
                      + M1SphCurv(iw_,cx1v,cx2v,cx3v,m,0,k,j,ip,odl,thrd,il,iu,jl,ju,kl,
                                  ku,M1_IW_EP));
            if (vlat) {
              od += 0.5*(M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,im,thrd,il,iu,jl,
                                     ju,kl,ku,M1_IW_EP)
                         + M1SphLat(iw_,vlt_,c0l,cx1v,cx2v,cx3v,m,0,k,j,ip,thrd,il,iu,jl,
                                    ju,kl,ku,M1_IW_EP));
            }
          }
        }
        Real fn = th*(f0n_(m,k,j,(i == ie+1 && cyclic) ? is : i)
                      - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*od);
        // the SAME blend the row was assembled with, so that the stored comoving face
        // flux (which the restart file carries, which the momentum deposit uses and which
        // the next iterate's reduced flux is built from) is the flux the solve applied.
        // The F0^n MEMORY term sits entirely in the F_diff branch: in the thin limit
        // alpha -> 1 and it must NOT survive, or the lagged flux would fight the upwind
        // HLL flux and the front would be damped exactly as it is under `central`.
        // (3c) the test is on the FLUX FORM, not on al: with implicit_blend_mode =
        // dissipation the central weight is 1 (al = 0) and the HLL part is still there.
        Real al = ifw_(m,M1_IFW_AL,k,j,i);
        if (aphll) {
          Real gh = ifw_(m,M1_IFW_HCL,k,j,i)*iw_(m,M1_IW_EP,k,j,im)
                    + ifw_(m,M1_IFW_HCR,k,j,i)*iw_(m,M1_IW_EP,k,j,ip)
                    + ifw_(m,M1_IFW_DG,k,j,i);
          if (musf) {
            const Real hcl = ifw_(m,M1_IFW_HCL,k,j,i), hcr = ifw_(m,M1_IFW_HCR,k,j,i);
            const int imm = (cyclic && im == is) ? ie : (im-1);
            const int ipp = (cyclic && ip == ie) ? is : (ip+1);
            if (hcl != 0.0) {
              gh += 0.25*hcl*iw_(m,imbf+M1_IM_SIG,k,j,im)
                    *(iw_(m,M1_IW_EP,k,j,ip) - iw_(m,M1_IW_EP,k,j,imm));
            }
            if (hcr != 0.0) {
              gh -= 0.25*hcr*iw_(m,imbf+M1_IM_SIG,k,j,ip)
                    *(iw_(m,M1_IW_EP,k,j,ipp) - iw_(m,M1_IW_EP,k,j,im));
            }
          }
          fn = (1.0 - al)*fn + (cl/ch)*gh;
        }
        f0_(m,k,j,i) = fn;
      }
    });
    if (cyclic) {
      // the wrapped face is stored twice; keep the two copies identical
      par_for("m1_impl_facewrap", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        f0_(m,k,j,ie+1) = f0_(m,k,j,is);
      });
    }
    par_for("m1_impl_f1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_F1,k,j,i) = 0.5*(f0_(m,k,j,i) + f0_(m,k,j,i+1))
                              + iw_(m,M1_IW_ADV,k,j,i)*iw_(m,M1_IW_EP,k,j,i);
      if (trans) {
        Real ep = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_F2,k,j,i) = 0.5*(f2_(m,k,j,i) + f2_(m,k,j+1,i))
                                + iw_(m,M1_IW_A2,k,j,i)*ep;
        if (thrd) {
          iw_(m,M1_IW_F3,k,j,i) = 0.5*(f3_(m,k,j,i) + f3_(m,k+1,j,i))
                                  + iw_(m,M1_IW_A3,k,j,i)*ep;
        }
      }
    });

    // (h) convergence
    // MILESTONE 3b phase B.  The Picard test alone is NOT enough once the transverse
    // couplings are lagged: |dE|/E can stall while the off-diagonal terms are still
    // moving.  The TRUE residual of the full 7-point system is measured separately (see
    // ImplicitTransverseTerms) and both have to be met.
    // m1-perf-0928: both maxima in ONE flat reduction and one MPI_MAX of two doubles
    // (was two MDRange reductions, two host syncs and two MPI_Allreduce); a max does not
    // depend on the order, so the values are bitwise the same.
    resid = 0.0;
    lresid = 0.0;
    Real rexcl = 0.0;
    if (!rmask) {
      const int fni = ie - is + 1, fnj = je - js + 1, fnk = ke - ks + 1;
      const bool ltr = trans;
      Kokkos::parallel_reduce("m1_impl_res",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*fnk*fnj*fni),
      KOKKOS_LAMBDA(const int n, Real &lmax, Real &lmx2) {
        int t = n/fni;
        const int i = is + (n - t*fni);
        const int j = js + (t % fnj);
        t /= fnj;
        const int k = ks + (t % fnk);
        const int m = t/fnk;
        Real r = iw_(m,M1_IW_RES,k,j,i);
        lmax = (r > lmax) ? r : lmax;
        if (ltr) {
          Real q = iw_(m,M1_IW_LRES,k,j,i);
          lmx2 = (q > lmx2) ? q : lmx2;
        }
      }, Kokkos::Max<Real>(resid), Kokkos::Max<Real>(lresid));
    } else {
      // implicit_res_dmin / _rmax: the excluded cells go to their own maximum (rexcl,
      // a diagnostic) and do not enter the stopping test
      const int fni = ie - is + 1, fnj = je - js + 1, fnk = ke - ks + 1;
      const bool ltr = trans;
      const Real rdmin = impl_res_dmin, rrmax = impl_res_rmax;
      auto rx1v = pmy_pack->pcoord->x1v;
      Kokkos::parallel_reduce("m1_impl_res_mask",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*fnk*fnj*fni),
      KOKKOS_LAMBDA(const int n, Real &lmax, Real &lmx2, Real &lmx3) {
        int t = n/fni;
        const int i = is + (n - t*fni);
        const int j = js + (t % fnj);
        t /= fnj;
        const int k = ks + (t % fnk);
        const int m = t/fnk;
        const bool ex = (rdmin > 0.0 && uh(m,IDN,k,j,i) < rdmin) ||
                        (rrmax > 0.0 && rx1v(m,i) > rrmax);
        Real r = iw_(m,M1_IW_RES,k,j,i);
        if (ex) {
          // a NaN would drop out of the max: map it to 1e300 (resid_fatal_masked)
          if (Kokkos::isnan(r)) {r = 1.0e300;}
          lmx3 = (r > lmx3) ? r : lmx3;
        } else {
          lmax = (r > lmax) ? r : lmax;
          if (ltr) {
            Real q = iw_(m,M1_IW_LRES,k,j,i);
            lmx2 = (q > lmx2) ? q : lmx2;
          }
        }
      }, Kokkos::Max<Real>(resid), Kokkos::Max<Real>(lresid), Kokkos::Max<Real>(rexcl));
    }
#if MPI_PARALLEL_ENABLED
    // the convergence test must be GLOBAL: with a partitioned column the ranks would
    // otherwise take different numbers of Picard passes and the gather would deadlock,
    // and even with rank-local columns a per-rank test makes the answer depend on the
    // decomposition.
    {
    Real rl[3] = {resid, lresid, rexcl}, rg[3];
    MPI_Allreduce(rl, rg, rmask ? 3 : 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    resid = rg[0];
    lresid = rg[1];
    if (rmask) {rexcl = rg[2];}
    }
#endif
    rexcl_last = rexcl;
    if (trans) {
      lresid /= rhsmax;
    } else {
      lresid = 0.0;
    }
    converged = (resid < impl_tol) && (!trans || (lresid < impl_lin_tol));
    if (impl_conv_est || !impl_lres_test) {
      // the Picard test, optionally without the confirming pass: q is the contraction
      // of the last two passes, and q/(1-q) times this change bounds what is left
      bool pc = (resid < impl_tol);
      if (impl_conv_est && !pc && rprev > 0.0) {
        Real q = resid/rprev;
        pc = (q < 0.5) && (resid*q/(1.0 - q) < impl_tol);
      }
      if (onep && it == 0 && !pc && !ocheck) {
        const Real qm = std::max(onep_qa[otyp], onep_qb[otyp]);
        // time2_one_pass_safety for the stage solves (otyp 1, 2)
        const Real sf = (otyp != 0 && t2_onep_s > 0.0) ? t2_onep_s : impl_onep_s;
        const Real qe = std::max(sf*qm, 1.0e-6);
        pc = (qe < 0.5) && (resid*qe/(1.0 - qe) < impl_tol);
        if (pc) {impl_onep_n += 1.0;}
      }
      bool lc = !trans || (lresid < impl_lin_tol) || (!impl_lres_test && bicg);
      converged = pc && lc;
    }
    if (gsd) {
      rhist.push_back(resid);
      // the stall test: from pass min on, the minimum of the last nw passes is not below
      // fac x the minimum before them, and resid rose at least twice inside the window
      auto stall_test = [&](int nw, int pmin, Real fac) {
        if (converged || (it + 1 < pmin) || (it + 1 <= nw)) {return false;}
        Real mwin = 1.0e300, mpre = 1.0e300;
        int nup = 0;
        for (int q = 0; q <= it; ++q) {
          if (q > it - nw) {
            mwin = std::min(mwin, rhist[q]);
            if (rhist[q] > rhist[q-1]) {++nup;}
          } else {
            mpre = std::min(mpre, rhist[q]);
          }
        }
        return (mwin >= fac*mpre) && (nup >= 2);
      };
      if (gsd && !gnsw && stall_test(impl_gn_sw, impl_gn_sw_min, impl_stall_fac)) {
        gnsw = true;
        impl_gn_nsw += 1.0;
        // the history restarts: the switched iteration is judged on its own passes
        rhist.clear();
        for (int q = 0; q <= it; ++q) {rhist.push_back(resid);}
      }
    }
    rprev = resid;
    if (it == 0) {ores0 = resid;}
    if (it == 1) {ores1 = resid;}
    if (plog) {ImplicitPicardLog(it, nin, resid, lresid, src_on);}
    if (impl_ncdump > 0 && it >= impl_maxit - impl_ncdump) {
      ImplicitNCDump(it, 0, gasx, igb, igr, igf, igy);
    }
    if (strc && it == impl_strace_p0 && !converged) {
      sloc = ImplicitNCDump(it, 2, gasx, igb, igr, igf, igy);
      sloc2 = ncd_hot_loc;
      if (sloc2 >= 0) {(void) ImplicitNCDump(it, 4, gasx, igb, igr, igf, igy, sloc2);}
      impl_strace_n += 1;
    } else if (strc && sloc >= 0 && it > impl_strace_p0 &&
               it < impl_strace_p0 + impl_strace) {
      (void) ImplicitNCDump(it, 2, gasx, igb, igr, igf, igy, sloc);
      if (sloc2 >= 0) {(void) ImplicitNCDump(it, 4, gasx, igb, igr, igf, igy, sloc2);}
    }
    // MILESTONE 3e: ACCELERATE.  Only on a pass that is followed by another one: the
    // state the step ENDS on must be the one the face fluxes of step (g) were built
    // from, so a converged pass -- and the last pass of a non-converged step -- keeps
    // the plain Picard iterate.  The accelerated E then has to reach the ghost cells
    // before the next pass assembles a row from it.
    if (accel && !converged && (it + 1 < impl_maxit)) {
      ImplicitAccelApply(it);
      if (trans) {
        // the acceleration rewrites the cell flux F1 as well as E
        ImplicitTransverseHalo(M1_NHALO_T);
      } else {
        ImplicitX1Halo(true);
      }
    }
  }
  TmrMark(7);
  if (tmr_on) {tmr_cnt[0] += 1.0;}
  ew_tight = false;
  if (onep) {
    if (ores1 >= 0.0 && ores0 > 0.0) {
      // a second pass was taken: its change measures the contraction
      onep_qb[otyp] = onep_qa[otyp];
      onep_qa[otyp] = ores1/ores0;
      onep_cnt[otyp] = 0.0;
      impl_onep_nchk += 1.0;
    } else {
      onep_cnt[otyp] += 1.0;
    }
  }
  if (onep_el && impl_onep_auto) {OnePassAuto(otyp, onep, converged && it == 1);}
  // implicit_predictor_order = 2 under hesdirk2: a backward-Euler step keeps ipred
  const bool pskip = (impl_pord == 2) && (time_scheme == M1_TIME_HESDIRK2) && !t2st;
  if (pred && !pskip) {
    const bool hh = have_hydro;
    const bool po2 = (impl_pord == 2);
    // the previous increment and its dt, for h = (g_new - g_old)/dt_old
    const bool hv = po2 && (p2 ? (pred2_ok && pred2_dt > 0.0)
                               : (pred_ok && pred_dt > 0.0));
    const Real dto = hv ? (p2 ? pred2_dt : pred_dt) : 1.0;
    const Real dtn = dt;
    par_for("m1_impl_pstore", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real de = iw_(m,M1_IW_EP,k,j,i) - (t2st ? fmax(u0_(m,M1_E,k,j,i), efl)
                                                    : iw_(m,M1_IW_EN,k,j,i));
      const Real dtp = hh ? (iw_(m,M1_IW_TP,k,j,i) - pd_(m,2,k,j,i)) : 0.0;
      if (po2) {
        pd_(m,3,k,j,i) = hv ? (de/dtn - pd_(m,0,k,j,i)/dto)/dto : 0.0;
        pd_(m,4,k,j,i) = hv ? (dtp/dtn - pd_(m,1,k,j,i)/dto)/dto : 0.0;
      }
      pd_(m,0,k,j,i) = de;
      pd_(m,1,k,j,i) = dtp;
    });
    if (p2) {
      pred2_ok = true;
      pred2_dt = dt;
    } else {
      pred_ok = true;
      pred_dt = dt;
    }
  }
  if (trans) {
    // refresh the stored transverse face fluxes with the CONVERGED E, so that the
    // persistent state the restart file carries, the momentum deposit and the derived
    // F_2/F_3 are the fluxes the solve ended on (exactly what step (g) does for x1).
    ImplicitTransverseTerms(false);
    impl_linmax = std::max(impl_linmax, lresid);
    impl_linsum += lresid;
  }
#if MPI_PARALLEL_ENABLED
  {
    // a column never crosses a rank here, but the ITERATION COUNT must be global or the
    // ranks would take different numbers of Picard passes and diverge
    int ilocal = it, iglob;
    MPI_Allreduce(&ilocal, &iglob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    it = iglob;
  }
#endif
  // hesdirk2: is the stage admissible?  Not when the Picard loop did not converge, when
  // a positivity fallback of the operator fired, or when the solved E or T is not
  // positive; the Driver then redoes the step from U^n with backward Euler.
  if (t2st) {
    const bool hh = have_hydro;
    const bool gq = have_hydro && coupling && dbgh;
    Real vmin = 1.0e300;
    // a flat 1-D range (m1-fast4): the rank-4 MDRange reduction ran at 360 us on the
    // wedge; a min is order-independent, so the result is the same bit for bit
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("m1_t2_adm",
    Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &lmin) {
      const int m = idx/nkji;
      int q = idx - m*nkji;
      const int k = q/nji + ks;
      q -= (k - ks)*nji;
      const int j = q/ni + js;
      const int i = q - (j - js)*ni + is;
      Real r = iw_(m,M1_IW_S2,k,j,i);
      if (hh) {r = fmin(r, iw_(m,M1_IW_TP,k,j,i));}
      if (gq) {
        // the gas internal energy the write-back will set (its step (a))
        Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*iw_(m,M1_IW_EP,k,j,i);
        r = fmin(r, iw_(m,M1_IW_EGN,k,j,i) - (cl/ch)*qq);
      }
      if (!(r > 0.0)) {r = -1.0;}
      lmin = (r < lmin) ? r : lmin;
    }, Kokkos::Min<Real>(vmin));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&vmin, &g, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
    vmin = g;}
#endif
    t2_fail = !converged || !(vmin > 0.0) || (od_now != impl_offdiag) ||
              (vimp_now != impl_vimp) ||
              ((t2s == M1_T2S_STAGE1) && (pmy_pack->pmesh->ncycle == t2_dbg_fail));
    if (t2_fail && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> hesdirk2 stage " << ((t2s == M1_T2S_STAGE1) ? 1 : 2)
                << " NOT ADMISSIBLE at cycle " << pmy_pack->pmesh->ncycle
                << " (converged=" << converged << " min(E,T)=" << vmin
                << "): the step is redone with backward Euler" << std::endl;
    }
    // DEBUG <rad_m1>/dbg_t2_admiss = N (default 0 = off): for the first N non-admissible
    // stages, which quantity went <= 0 and where.  Read-only diagnostics after the
    // decision (nothing written back), so results are bitwise with the key on or off.
    if (t2_fail && (t2_dbg_adm_n < t2_dbg_adm)) {
      t2_dbg_adm_n += 1;
      T2AdmissDebug(uh, u0_, t2i_, cl, ch, hh, gq, t2s, it);
    }
  }
  // floor clips (counting only, nothing written): cells whose solved E is at or below
  // e_floor (the iterate is floored to it) and cells whose written-back gas internal
  // energy (the step (a) of the write-back below) is <= 0, in every solve, BE included
  {
    const bool hh = have_hydro;
    const bool gq = have_hydro && coupling && dbgh;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Real nfe = 0.0, nfg = 0.0;
    Kokkos::parallel_reduce("m1_flr_cnt",
    Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &se, Real &sg) {
      const int m = idx/nkji;
      int q = idx - m*nkji;
      const int k = q/nji + ks;
      q -= (k - ks)*nji;
      const int j = q/ni + js;
      const int i = q - (j - js)*ni + is;
      if (!(iw_(m,M1_IW_S2,k,j,i) > efl)) se += 1.0;
      if (hh && gq) {
        Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*iw_(m,M1_IW_EP,k,j,i);
        if (!(iw_(m,M1_IW_EGN,k,j,i) - (cl/ch)*qq > 0.0)) sg += 1.0;
      }
    }, Kokkos::Sum<Real>(nfe), Kokkos::Sum<Real>(nfg));
#if MPI_PARALLEL_ENABLED
    {Real a[2] = {nfe, nfg}, g[2];
    MPI_Allreduce(a, g, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    nfe = g[0]; nfg = g[1];}
#endif
    flr_ne += nfe;
    if (impl_ncdump > 0 && nfe > 0.0) {
      ImplicitNCDump(it, 1, gasx, igb, igr, igf, igy);
    }
    flr_ng += nfg;
    // DEBUG dbg_t2_admiss (m1-positivity): a solve that ENDS with E <= e_floor or a
    // written-back eint <= 0 anywhere (the BE redo included) lists its cells and rows
    if ((nfe > 0.0 || nfg > 0.0) && (t2_dbg_adm_n < t2_dbg_adm)) {
      t2_dbg_adm_n += 1;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> POS_TRIGGER floor cycle " << pmy_pack->pmesh->ncycle
                  << " solve " << (t2st ? "hesdirk2-stage" : "be") << " E<=e_floor "
                  << nfe << " eint_wb<=0 " << nfg << " vimp_now " << vimp_now
                  << std::endl;
      }
      T2AdmissDebug(uh, u0_, t2i_, cl, ch, hh, gq, t2s, it);
    }
  }
  if (rmask) {
    impl_res_exmax = std::max(impl_res_exmax, rexcl_last);
    if (rexcl_last >= impl_tol) {impl_res_nex += 1.0;}
  }
  impl_nstep += 1.0;
  impl_itsum += static_cast<Real>(it);
  impl_itmax = std::max(impl_itmax, static_cast<Real>(it));
  if (!converged) {
    impl_nfail += 1.0;
    // MILESTONE 3e: WHERE the outer iteration stalled.  One line per non-converged step
    // with the cell that owns the largest Picard residual -- the depth index i is what
    // says whether the stall is in the optically thin top or at the base.  Nothing here
    // changes any number; a step that converges (every step of every earlier gate)
    // prints nothing.
    using MaxLoc = Kokkos::MaxLoc<Real,int>;
    MaxLoc::value_type mloc;
    const int nj = je - js + 1, ni = ie - is + 1;
    const Real rdmin = impl_res_dmin, rrmax = impl_res_rmax;
    auto rx1v = pmy_pack->pcoord->x1v;
    Kokkos::parallel_reduce("m1_impl_resloc",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                  MaxLoc::value_type &lmx) {
      Real r = iw_(m,M1_IW_RES,k,j,i);
      if (rmask && ((rdmin > 0.0 && uh(m,IDN,k,j,i) < rdmin) ||
                    (rrmax > 0.0 && rx1v(m,i) > rrmax))) {r = 0.0;}
      if (r > lmx.val) {
        lmx.val = r;
        lmx.loc = ((m*(ke-ks+1) + (k-ks))*nj + (j-js))*ni + (i-is);
      }
    }, MaxLoc(mloc));
    int lmb = mloc.loc/((ke-ks+1)*nj*ni);
    int lrem = mloc.loc - lmb*(ke-ks+1)*nj*ni;
    int lk = lrem/(nj*ni);
    int lj = (lrem - lk*nj*ni)/ni;
    int li = lrem - lk*nj*ni - lj*ni;
    // stall_1002: every rank's own worst cell and its value
    std::cout << "<rad_m1> NC-local rank=" << global_variable::my_rank << " cycle="
              << pmy_pack->pmesh->ncycle << " max=" << mloc.val << " (m,k,j,i)=(" << lmb
              << "," << (lk+ks) << "," << (lj+js) << "," << (li+is) << ")" << std::endl;
    if (vgd_on && std::getenv("VGD_NCDIAG") != nullptr) {
      // vet_gd diagnosis (env VGD_NCDIAG): the 3x3 lateral neighbourhood of the worst
      // cell: D_rr (tau_ten slot 0), E, F1, Picard residual
      const int kk = lk+ks, jj = lj+js, ii = li+is;
      auto pr = std::make_pair(ii, ii+1);
      auto pk = std::make_pair(kk-1, kk+2);
      auto pj = std::make_pair(jj-1, jj+2);
      auto dh = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                  Kokkos::subview(tau_ten, lmb, 0, pk, pj, pr));
      auto eh = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                  Kokkos::subview(iw, lmb, static_cast<int>(M1_IW_EP), pk, pj, pr));
      auto fh = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                  Kokkos::subview(iw, lmb, static_cast<int>(M1_IW_F1), pk, pj, pr));
      auto rh = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                  Kokkos::subview(iw, lmb, static_cast<int>(M1_IW_RES), pk, pj, pr));
      std::cout << "<rad_m1> NC-diag rank=" << global_variable::my_rank << " cycle="
                << pmy_pack->pmesh->ncycle << " edge(k,j)=(" << (kk == ks || kk == ke)
                << "," << (jj == js || jj == je) << ")";
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          std::cout << " [" << (a-1) << (b-1) << " D " << dh(a,b,0) << " E " << eh(a,b,0)
                    << " F " << fh(a,b,0) << " R " << rh(a,b,0) << "]";
        }
      }
      std::cout << std::endl;
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> Picard NON-CONVERGED after " << it << " passes:"
                << " resid=" << resid << " (tol " << impl_tol << ")"
                << " lin_resid=" << lresid;
      if (rmask) {std::cout << " masked_resid=" << rexcl_last;}
      std::cout << " worst cell of rank 0 (m,k,j,i)=(" << lmb << "," << (lk+ks) << ","
                << (lj+js) << "," << (li+is) << ")" << std::endl;
    }
    if (impl_res_fatal > 0.0 && (!(resid <= impl_res_fatal) || !std::isfinite(lresid))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_m1> Picard solve DIVERGED: resid=" << resid
                << " lin_resid=" << lresid << " > implicit_resid_fatal=" << impl_res_fatal
                << " at time " << pmy_pack->pmesh->time << " cycle "
                << pmy_pack->pmesh->ncycle << " (rank " << global_variable::my_rank
                << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // implicit_resid_fatal_masked: the excluded cells are still solved every pass; a
  // diverging one must stop the run before its state spreads (NaN maps to 1e300)
  if (rmask && impl_res_fatal_masked > 0.0 &&
      !(rexcl_last <= impl_res_fatal_masked && rexcl_last < 1.0e300)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<rad_m1> Picard solve DIVERGED in the masked cells: "
              << "masked_resid=" << rexcl_last << " (1e300 = NaN) > "
              << "implicit_resid_fatal_masked=" << impl_res_fatal_masked
              << " (resid=" << resid << ", converged=" << converged << ") at time "
              << pmy_pack->pmesh->time << " cycle " << pmy_pack->pmesh->ncycle
              << " (rank " << global_variable::my_rank << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // MILESTONE 3g: the per-cell counters of the gas coupling, reduced ONCE per step.
  if (gasx) {
    Real sfb = 0.0, sms = 0.0;
    Kokkos::parallel_reduce("m1_impl_gcnt",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lf) {
      lf += iw_(m,igf,k,j,i);
    }, Kokkos::Sum<Real>(sfb));
    Kokkos::parallel_reduce("m1_impl_gms",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
      ls += iw_(m,igm,k,j,i);
    }, Kokkos::Sum<Real>(sms));
    Real ncell = static_cast<Real>(nmb1+1)*static_cast<Real>(ke-ks+1)
                 *static_cast<Real>(je-js+1)*static_cast<Real>(ie-is+1);
    // diagnostic (<rad_m1>/report_newton_fb, default off): WHERE the Newton fallbacks of
    // this step happened (first 8 cells of this rank: cycle, time, m, k, j, i, r)
    if (sfb > 0.0 && pmy_pack->pmesh->ncycle >= 0 &&
        pin_report_newton_fb) {
      // accel-1009: copy only the fallback-count component (the whole work array was
      // copied before: ~0.26 s per step on A100 for AG Car A); the printed lines are
      // the same
      if (nfb_buf.extent(0) != static_cast<size_t>(nmb1 + 1) ||
          nfb_buf.extent(1) != iw_.extent(2) || nfb_buf.extent(2) != iw_.extent(3) ||
          nfb_buf.extent(3) != iw_.extent(4)) {
        Kokkos::realloc(nfb_buf, nmb1 + 1, iw_.extent(2), iw_.extent(3), iw_.extent(4));
      }
      auto nb_ = nfb_buf;
      par_for("m1_impl_nfb_cp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        nb_(m,k,j,i) = iw_(m,igf,k,j,i);
      });
      auto hfb4 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), nfb_buf);
      auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmy_pack->pcoord->x1v);
      auto hfb = [&](const int m, const int, const int k, const int j, const int i) {
        return hfb4(m,k,j,i);
      };
      int nrep = 0;
      for (int m = 0; m <= nmb1 && nrep < 8; ++m) {
        for (int k = ks; k <= ke && nrep < 8; ++k) {
          for (int j = js; j <= je && nrep < 8; ++j) {
            for (int i = is; i <= ie && nrep < 8; ++i) {
              if (hfb(m,igf,k,j,i) > 0.0) {
                std::cout << "<rad_m1> NEWTON-FALLBACK cycle=" << pmy_pack->pmesh->ncycle
                          << " time=" << pmy_pack->pmesh->time << " (m,k,j,i)=(" << m
                          << "," << k << "," << j << "," << i << ") r=" << hx1v(m,i)
                          << " count=" << hfb(m,igf,k,j,i) << std::endl;
                ++nrep;
              }
            }
          }
        }
      }
    }
#if MPI_PARALLEL_ENABLED
    {Real lo[3] = {sfb, sms, ncell}, gl[3];
    MPI_Allreduce(lo, gl, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    sfb = gl[0];
    sms = gl[1];
    ncell = gl[2];}
#endif
    newt_nfb += sfb;
    ec_nmiss += sms;
    gas_ncell += ncell*static_cast<Real>(it);
  }

  // MILESTONE 3g: ONE TRUE-TABLE evaluation per cell at the state the step ends on, to
  // MEASURE what the cache cost.  Two numbers are taken: the relative error of e(T')
  // itself, and the relative error of the energy q = SRCR - SRCB E' the gas actually
  // exchanges with the radiation -- the only channel through which the cache can reach
  // the evolved state.
  //
  // Nothing is CORRECTED here, deliberately.  The gas energy is set from the ASSEMBLED
  // row (see the write-back below), which is what makes e_gas + (c/chat) E change by the
  // fluxes and the work term alone to round-off; overwriting q after the solve with a
  // re-evaluated one would break exactly that algebraic balance.  Consistency with the
  // real table is instead structural: the evolved gas state is (rho, e_gas), T' is not
  // persistent, and the next step re-inverts e_gas through the real table.
  if (usec && impl_eccheck && src_on &&
      (impl_eccheck_every == 1 || pmy_pack->pmesh->ncycle % impl_eccheck_every == 0)) {
    auto eos = flr.eos;
    Real emx = 0.0, qmx = 0.0;
    // m1-fast4: a flat range with LaunchBounds<256,1> (the MDRange form spilled 55
    // VGPRs); a max is order-independent, so the result is bitwise the same
    const int eni = ie - is + 1, enji = (je - js + 1)*eni, enkji = (ke - ks + 1)*enji;
    Kokkos::parallel_reduce("m1_impl_eck",
    Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>(DevExeSpace(), 0,
                                                                  (nmb1 + 1)*enkji),
    KOKKOS_LAMBDA(const int idx, Real &lmax) {
      const int m = idx/enkji;
      int q = idx - m*enkji;
      const int k = q/enji + ks;
      q -= (k - ks)*enji;
      const int j = q/eni + js;
      const int i = q - (j - js)*eni + is;
      iw_(m,igm,k,j,i) = 0.0;
      Real rkpv = opac_(m,M1_OP_P,k,j,i);
      Real rkev = opac_(m,M1_OP_E,k,j,i);
      if (rkpv == 0.0 && rkev == 0.0) {return;}
      Real dd = uh(m,IDN,k,j,i);
      Real tk = iw_(m,M1_IW_TP,k,j,i);
      Real ce, ccv;
      if (!M1EosCacheEval(eos, ec_, m, k, j, i, ecnt, dd, tk, ce, ccv)) {return;}
      Real te, pp, cr, ct, tcv;
      eos.ThermoAt(dd, tk, te, pp, cr, ct, tcv);
      lmax = fmax(lmax, fabs(ce - te)/fmax(fabs(te), 1.0e-300));
      // the same row the pass assembled, re-made with the TRUE table
      Real t3 = tk*tk*tk, t4 = t3*tk;
      Real de0 = iw_(m,M1_IW_DE0,k,j,i);
      Real ep = iw_(m,M1_IW_EP,k,j,i);
      Real emis = dt*ch*rkpv*ar;
      Real bkc = dd*ccv + 4.0*cl*dt*rkpv*ar*t3;
      Real bkt = dd*tcv + 4.0*cl*dt*rkpv*ar*t3;
      Real rkc = iw_(m,M1_IW_EGN,k,j,i) - ce - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
      Real rkt = iw_(m,M1_IW_EGN,k,j,i) - te - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
      Real qc = emis*t4 - dt*ch*rkev*de0 - (dt*ch*rkev)*ep;
      Real qt = qc;
      if (bkc > 0.0) {qc += emis*4.0*t3*(rkc + cl*dt*rkev*ep)/bkc;}
      if (bkt > 0.0) {qt += emis*4.0*t3*(rkt + cl*dt*rkev*ep)/bkt;}
      iw_(m,igm,k,j,i) = fabs(qt - qc)/fmax(fabs(qt), 1.0e-300);
    }, Kokkos::Max<Real>(emx));
    Kokkos::parallel_reduce("m1_impl_eckq",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      lmax = fmax(lmax, iw_(m,igm,k,j,i));
    }, Kokkos::Max<Real>(qmx));
#if MPI_PARALLEL_ENABLED
    {Real lo[2] = {emx, qmx}, gl[2];
    MPI_Allreduce(lo, gl, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    emx = gl[0];
    qmx = gl[1];}
#endif
    ec_emax = std::max(ec_emax, emx);
    ec_tmax = std::max(ec_tmax, qmx);
  }

  //------------------------------------------------------------------------- write back
  // m1-sp-order2b: the cell flux at the stage's own velocity (hesdirk2 + implicit_vimp;
  // time2_vstage, default true on sp, false = the old stage-start form)
  const bool vfx = t2st && impl_vimp && t2_fvnew;
  const bool coupling_ = coupling;  // local: a member read in a kernel captures this
  // m1-positivity (rad_m1.hpp): the energy-conserving gas eint limiter and E floor
  const bool pgas = impl_pos_gas && have_hydro && coupling && dbgh;
  const bool pflr = impl_pos_floor;
  const bool pfs = impl_pos_floor_s2;   // F3
  const bool pany = pgas || pflr;
  const Real pgf = impl_pos_gas_frac;
  auto peos = flr.eos;
  auto pc_ = pos_cnt_d;
  const bool psph = sph_geom;
  auto pvol = pmy_pack->pcoord->volume;
  // implicit_face_weight = distance (sp): the face kappa_T of the momentum deposit is
  // the one the face equation used
  auto cx1f = pmy_pack->pcoord->xx1f;
  const bool fwd = psph && impl_face_wdist;
  // STAGE CS2 (C3), cubed sphere: the transverse face fluxes F0 are FACE-NORMAL, and
  // the hydro momentum and the stored cell F are COVARIANT on the panel basis (the seam
  // transform of pbval_u / pbval_th takes them as such).  With a, b the face-normal x2,
  // x3 values at the cell and (c, s) = (cos, sin) of the cell's angle,
  //   F.e_xi = (a + c b)/s,  F.e_eta = (b + c a)/s
  // (n_xi = (e_xi - c e_eta)/s); the same map turns the face-normal momentum deposit
  // into the covariant dm2, dm3.  The kinetic energy is 0.5 m_a v^a and the work
  // v^a dm_a with the contravariant v^a = g^ab m_b/rho.  `if (csw)` overwrites below;
  // the Cartesian and sp arithmetic is untouched.
  const bool csw = cs_geom;
  auto cclw = pmy_pack->pcoord->cos_cell;
  auto csnw = pmy_pack->pcoord->sin_cell;
  const bool rcpw = impl_real_couple;
  if (!rcpw) {
    par_for("m1_impl_wb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
#define M1_RCP 0
#include "rad_m1/rad_m1_impl_wb_kernel.hpp"
#undef M1_RCP
    });
  } else {
    // implicit_realisable_coupling: the same kernel with the key on
    par_for("m1_impl_wb_rc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
#define M1_RCP 1
#include "rad_m1/rad_m1_impl_wb_kernel.hpp"  // NOLINT(build/include)
#undef M1_RCP
    });
  }
  if (vfx && have_hydro && feedback && csw && trans && thrd) {
    // STAGE CS3: the vimp re-forming of the cell F on the cubed sphere, as a kernel of
    // its own (an overwrite inside the write-back above changed the GPU code of the
    // sp/Cartesian path).  The work array is face-normal, so the transverse velocities
    // are s v^a of the kicked covariant momentum, which the write-back has just stored
    // in uh; E is the final one (the limits are idempotent in E).
    par_for("m1_impl_wb_csvfx", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real idg = 1.0/fmax(uh(m,IDN,k,j,i), 1.0e-300);
      const Real w1 = uh(m,IM1,k,j,i)*idg;
      const Real c = cclw(m,k,j), sn = csnw(m,k,j), s2 = 1.0 - c*c;
      const Real q2 = uh(m,IM2,k,j,i), q3 = uh(m,IM3,k,j,i);
      const Real u2 = sn*(q2 - c*q3)*idg/s2, u3 = sn*(q3 - c*q2)*idg/s2;
      Real chi = iw_(m,M1_IW_WCHI,k,j,i);
      Real n1 = iw_(m,M1_IW_N1,k,j,i), n2 = iw_(m,M1_IW_N2,k,j,i);
      Real n3 = iw_(m,M1_IW_N3,k,j,i);
      Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2), d33 = M1EddDiag(chi,n3);
      Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
      Real d23 = M1EddOff(chi,n2,n3);
      if (dfull) {
        d11 = vd_(m,M1_VET_D11,k,j,i);
        d22 = vd_(m,M1_VET_D11+1,k,j,i);
        d33 = vd_(m,M1_VET_D11+2,k,j,i);
        d12 = vd_(m,M1_VET_D11+3,k,j,i);
        d13 = vd_(m,M1_VET_D11+4,k,j,i);
        d23 = vd_(m,M1_VET_D11+5,k,j,i);
      }
      const Real es = iw_(m,M1_IW_EP,k,j,i);
      const Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
      const Real g2l = f2_(m,k,j,i), g2r = f2_(m,k,j+1,i);
      const Real g3l = f3_(m,k,j,i), g3r = f3_(m,k+1,j,i);
      Real fp1 = 0.5*(fl + fr) + (w1 + (w1*d11 + u2*d12 + u3*d13))*es;
      const Real a2 = 0.5*(g2l + g2r) + (u2 + (w1*d12 + u2*d22 + u3*d23))*es;
      const Real a3 = 0.5*(g3l + g3r) + (u3 + (w1*d13 + u2*d23 + u3*d33))*es;
      Real fp2 = (a2 + c*a3)/sn;   // covariant, as the write-back's transform
      Real fp3 = (a3 + c*a2)/sn;
      Real ep = u0_(m,M1_E,k,j,i);
      M1ApplyLimitsCs(cl, efl, c, ep, fp1, fp2, fp3);
      u0_(m,M1_E,k,j,i) = ep;
      u0_(m,M1_F1,k,j,i) = fp1;
      u0_(m,M1_F2,k,j,i) = fp2;
      u0_(m,M1_F3,k,j,i) = fp3;
    });
  }

  // force_reference_work = split (ke-dt-0926): the write-back above handed the gas the
  // work of the FULL force and took it from E.  Here the gas keeps only the residual
  // kick's work (its exact kinetic-energy change; the WB source gives it the reference
  // work at its own stage), and E pays residual + reference (v' dt rho arad_ref, at the
  // stage value's velocity v').  A separate kernel: the default write-back is untouched.
  if (fws && !fref_wsplit_ok) {
    ImplFatal("<rad_m1>/force_reference_work = split needs a problem generator whose WB "
              "source gives the gas the reference work (box_convection, sph_wedge)");
  }
  if (fws && have_hydro && gas_feedback && coupling && dbg_gas_force) {
    const Real fkw = 1.0/dt;
    const Real chw = chat, clw = c_light;
    auto kkw_ = (t2s == M1_T2S_STAGE1) ? t2k2 : t2k1;
    const bool t2kw = t2k;
    const bool cycw = cyclic;
    const int bclw = ibc_x1min, bchw = ibc_x1max;
    const int nbw = part_nblk;
    auto posw_ = part_pos;
    // time_scheme = be (fref-split-cons-1006): E pays exactly the reference work the
    // gas got in the hydro stages of this step (fref_wacc, the share dt/dt_mesh of it
    // per solve), not v' dt rho a_ref: a separate kernel, the other paths are untouched
    if (FrefWaccOn()) {
      auto wacc_ = FrefWacc();
      const Real fra = dt/pmy_pack->pmesh->dt;
      par_for("m1_impl_fws_x", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const int ipw = posw_.d_view(m);
        const bool efw = !cycw && ((i == is && ipw == 0 && bclw == M1_IBC_EFIX) ||
                                   (i == ie && ipw == nbw-1 && bchw == M1_IBC_EFIX));
        const Real dd = uh(m,IDN,k,j,i);
        const Real idg = 1.0/fmax(dd, 1.0e-300);
        const Real v1 = iw_(m,M1_IW_V1,k,j,i);
        const Real m0 = dd*v1;
        const Real dmref = dt*dd*aref_(m,k,j,i);
        const Real dmr = uh(m,IM1,k,j,i) - m0;
        const Real dm1 = dmr + dmref;
        const Real wf = 0.5*(v1 + (m0 + dm1)*idg)*dm1;
        const Real wr = 0.5*(v1 + (m0 + dmr)*idg)*dmr;
        uh(m,IEN,k,j,i) -= (wf - wr);
        if (!efw) {
          u0_(m,M1_E,k,j,i) += (chw/clw)*(wf - wr - fra*wacc_(m,k,j,i));
        }
      });
    } else {
    par_for("m1_impl_fws", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int ipw = posw_.d_view(m);
      const bool efw = !cycw && ((i == is && ipw == 0 && bclw == M1_IBC_EFIX) ||
                                 (i == ie && ipw == nbw-1 && bchw == M1_IBC_EFIX));
      const Real dd = uh(m,IDN,k,j,i);
      const Real idg = 1.0/fmax(dd, 1.0e-300);
      const Real v1 = iw_(m,M1_IW_V1,k,j,i);
      const Real m0 = dd*v1;
      const Real dmref = dt*dd*aref_(m,k,j,i);
      const Real dmr = uh(m,IM1,k,j,i) - m0;
      const Real dm1 = dmr + dmref;
      const Real wf = 0.5*(v1 + (m0 + dm1)*idg)*dm1;
      const Real wr = 0.5*(v1 + (m0 + dmr)*idg)*dmr;
      uh(m,IEN,k,j,i) -= (wf - wr);
      if (t2kw) {kkw_(m,M1_T2_EN,k,j,i) -= (wf - wr)*fkw;}
      if (!efw) {
        const Real de = (chw/clw)*(wf - wr - (m0 + dmr)*idg*dmref);
        u0_(m,M1_E,k,j,i) += de;
        if (t2kw) {kkw_(m,M1_T2_E,k,j,i) += de*fkw;}
      }
    });
    }
  }

  // time2_vstage: the rows took the work of the last pass's kick, which the write-back
  // subtracted again: give it back (E and the slope; the gas keeps the full work)
  if (t2_wk) {
    ImplicitWorkRow(false);
  }
  // hesdirk2: the face part of the slope, over the active faces
  if (t2k) {
    const Real fk = 1.0/dt;
    par_for("m1_t2_kf1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      kk_(m,M1_T2_F1,k,j,i) = (f0_(m,k,j,i) - f0n_(m,k,j,i))*fk;
    });
    if (trans) {
      auto f2n_ = f0x2n;
      par_for("m1_t2_kf2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        kk_(m,M1_T2_F1+1,k,j,i) = (f2_(m,k,j,i) - f2n_(m,k,j,i))*fk;
      });
      if (thrd) {
        auto f3n_ = f0x3n;
        par_for("m1_t2_kf3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          kk_(m,M1_T2_F1+2,k,j,i) = (f3_(m,k,j,i) - f3n_(m,k,j,i))*fk;
        });
      }
    }
  }

  // DIAGNOSTIC (env VGD_RDIAG = N > 0, vet_gd): realizability of the tensor against the
  // new (E, F) and the sign structure of the last pass's stencil, every N solves
  if (vet_col && tau_ready) {VetGdRealDiag();}
  // time2_vet_col = rebuild: E and T of the stage-1 solution for the stage-2 build
  if (vet_col && t2_vcmode == 2 && t2s == M1_T2S_STAGE1) {
    Time2VetColSaveY1();
  }
  // DIAGNOSTIC <rad_m1>/dbg_cell_lo..hi (default off): per-cell budget after the solve
  if (dbg_cell_lo >= 0 && have_hydro) {
    const int dlo = is + dbg_cell_lo, dhi = is + dbg_cell_hi;
    auto eosd = flr.eos;
    const Real ard = arad;
    const int ncy = pmy_pack->pmesh->ncycle;
    Kokkos::fence();
    auto ttd = tau_ten;
    const bool ttok = tau_ready;
    par_for("m1_dbg_cell", DevExeSpace(), 0, 0, ks, ks, js, js+1, dlo, dhi,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dd = uh(m,IDN,k,j,i);
      const Real ek = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                           SQR(uh(m,IM3,k,j,i)))/dd;
      Real eg = uh(m,IEN,k,j,i) - ek - (etg ? dd*phicc(m,k,j,i) : 0.0);
      const Real tg = eosd.Temperature(dd, fmax(eg, 1.0e-300));
      const Real es = iw_(m,M1_IW_EP,k,j,i), ef = u0_(m,M1_E,k,j,i);
      const Real tps = iw_(m,M1_IW_TP,k,j,i);
      const Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*es;
      const Real fkd = ttok ? ttd(m,0,k,j,i) : -1.0;
      printf("DBGC c=%d j=%d fK-1/3=%.4e i=%d EN=%.10e Esolve=%.10e Efin=%.10e Tp=%.10e "
             "Trs/Tp-1=%.3e Trf/Tg-1=%.3e Tg/Tp-1=%.3e q=%.4e EGN=%.10e eg=%.10e "
             "v0=%.4e v1=%.4e\n", ncy, j - js, fkd - 1.0/3.0, i - is,
             iw_(m,M1_IW_EN,k,j,i), es, ef, tps,
             pow(es/ard, 0.25)/tps - 1.0, pow(ef/ard, 0.25)/tg - 1.0, tg/tps - 1.0, qq,
             iw_(m,M1_IW_EGN,k,j,i), eg, iw_(m,M1_IW_V1,k,j,i), uh(m,IM1,k,j,i)/dd);
    });
    Kokkos::fence();
  }
  if (vetsc) {Kokkos::fence(); vet_itime += vtimer.seconds();}
  impl_lin_tol = t2_lin_save;
  TmrMark(8);
  if (dbg_etally) {
    Real o[4];
    DbgEnergyTally(o);
    if (global_variable::my_rank == 0) {
      std::printf("ETALLY1 cycle=%d stage=%d ien=%.16e e=%.16e wacc=%.16e\n",
                  pmy_pack->pmesh->ncycle, stage, o[0], o[1], o[2]);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::T2AdmissDebug
//! \brief DEBUG (<rad_m1>/dbg_t2_admiss): after a non-admissible hesdirk2 stage, per
//! quantity q = solved E (M1_IW_S2), solved T (M1_IW_TP), written-back gas eint,
//! stage OLD vector E (M1_IW_EN) and gas eint (M1_IW_EGN), stage START E (u0): the
//! number of active cells with q <= 0 and, where there are any, the cell of the minimum
//! of q with its state.  rho/<rho> uses the rank-local unweighted mean over the shell.
//! Every rank prints its own cells.  Nothing is written.

void RadiationM1::T2AdmissDebug(DvceArray5D<Real> uh, DvceArray5D<Real> u0_,
                                DvceArray5D<Real> t2i_, const Real cl, const Real ch,
                                const bool hh, const bool gq, const int t2s,
                                const int it) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto x1v = pmy_pack->pcoord->x1v;
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  const int ntot = (nmb1 + 1)*nkji;
  const char *nm[6] = {"E_solved", "T_solved", "eint_writeback", "E_old_vector",
                       "eint_old_vector", "E_stage_start"};
  const int rank = global_variable::my_rank;
  std::cout << "<rad_m1> DBG_T2_ADMISS rank " << rank << " cycle "
            << pmy_pack->pmesh->ncycle << " time " << pmy_pack->pmesh->time
            << " stage " << ((t2s == M1_T2S_STAGE1) ? 1 : 2) << " picard_passes " << it
            << std::endl;
  for (int q = 0; q < 6; ++q) {
    if ((q == 1 || q == 4) && !hh) continue;
    if (q == 2 && !gq) continue;
    auto qval = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v;
      switch (q) {
        case 0: v = iw_(m,M1_IW_S2,k,j,i); break;
        case 1: v = iw_(m,M1_IW_TP,k,j,i); break;
        case 2: v = iw_(m,M1_IW_EGN,k,j,i) - (cl/ch)*(iw_(m,M1_IW_SRCR,k,j,i) -
                    iw_(m,M1_IW_SRCB,k,j,i)*iw_(m,M1_IW_EP,k,j,i)); break;
        case 3: v = iw_(m,M1_IW_EN,k,j,i); break;
        case 4: v = iw_(m,M1_IW_EGN,k,j,i); break;
        default: v = u0_(m,M1_E,k,j,i); break;
      }
      return v;
    };
    int nbad = 0;
    Kokkos::parallel_reduce("m1_t2dbg_n", Kokkos::RangePolicy<DevExeSpace>(0, ntot),
    KOKKOS_LAMBDA(const int idx, int &n) {
      const int m = idx/nkji;
      int r = idx - m*nkji;
      const int k = r/nji + ks;
      r -= (k - ks)*nji;
      const int j = r/ni + js;
      const int i = r - (j - js)*ni + is;
      if (!(qval(m,k,j,i) > 0.0)) n += 1;
    }, Kokkos::Sum<int>(nbad));
    if (nbad == 0) {
      std::cout << "  " << nm[q] << ": no cell <= 0" << std::endl;
      continue;
    }
    using ML = Kokkos::MinLoc<Real,int>;
    typename ML::value_type mloc;
    Kokkos::parallel_reduce("m1_t2dbg_min", Kokkos::RangePolicy<DevExeSpace>(0, ntot),
    KOKKOS_LAMBDA(const int idx, typename ML::value_type &lm) {
      const int m = idx/nkji;
      int r = idx - m*nkji;
      const int k = r/nji + ks;
      r -= (k - ks)*nji;
      const int j = r/ni + js;
      const int i = r - (j - js)*ni + is;
      Real v = qval(m,k,j,i);
      if (v < lm.val) {lm.val = v; lm.loc = idx;}
    }, ML(mloc));
    const int idx = mloc.loc;
    const int m = idx/nkji;
    int r = idx - m*nkji;
    const int k = r/nji + ks;
    r -= (k - ks)*nji;
    const int j = r/ni + js;
    const int i = r - (j - js)*ni + is;
    // the cell's state and the rank-local shell mean of rho at its i
    DvceArray1D<Real> cv("m1_t2dbg_cv", 16);
    par_for("m1_t2dbg_cell", DevExeSpace(), 0, 0,
    KOKKOS_LAMBDA(const int n) {
      cv(0) = x1v(m,i);
      cv(1) = hh ? uh(m,IDN,k,j,i) : 0.0;
      cv(2) = iw_(m,M1_IW_S2,k,j,i);
      cv(3) = iw_(m,M1_IW_EN,k,j,i);
      cv(4) = u0_(m,M1_E,k,j,i);
      cv(5) = t2i_(m,M1_T2_E,k,j,i);
      cv(6) = iw_(m,M1_IW_EGN,k,j,i);
      cv(7) = iw_(m,M1_IW_EGN,k,j,i) - (cl/ch)*(iw_(m,M1_IW_SRCR,k,j,i) -
              iw_(m,M1_IW_SRCB,k,j,i)*iw_(m,M1_IW_EP,k,j,i));
      cv(8) = hh ? t2i_(m,M1_T2_EN,k,j,i) : 0.0;
      cv(9) = iw_(m,M1_IW_TP,k,j,i);
      cv(10) = iw_(m,M1_IW_V1,k,j,i);
      cv(11) = iw_(m,M1_IW_V2,k,j,i);
      cv(12) = iw_(m,M1_IW_V3,k,j,i);
      cv(13) = iw_(m,M1_IW_EP,k,j,i);
      cv(14) = iw_(m,M1_IW_RES,k,j,i);
    });
    Real rsum = 0.0;
    if (hh) {
      const int nk = ke - ks + 1, nj = je - js + 1, nkj = nk*nj;
      Kokkos::parallel_reduce("m1_t2dbg_sh",
      Kokkos::RangePolicy<DevExeSpace>(0, (nmb1 + 1)*nkj),
      KOKKOS_LAMBDA(const int id, Real &sum) {
        const int mm = id/nkj;
        const int kk = (id - mm*nkj)/nj + ks;
        const int jj = id - mm*nkj - (kk - ks)*nj + js;
        sum += uh(mm,IDN,kk,jj,i);
      }, Kokkos::Sum<Real>(rsum));
      rsum /= static_cast<Real>((nmb1 + 1)*nkj);
    }
    auto hv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cv);
    std::cout << "  " << nm[q] << ": " << nbad << " cells <= 0; min " << mloc.val
              << " at (m,k,j,i)=(" << m << "," << k << "," << j << "," << i << ") r="
              << hv(0) << " rho=" << hv(1) << " rho/<rho>_shell="
              << ((rsum > 0.0) ? hv(1)/rsum : 0.0) << std::endl
              << "    E: solved " << hv(2) << " old_vector " << hv(3) << " stage_start "
              << hv(4) << " t2inc " << hv(5) << " last_iterate " << hv(13)
              << " cell_resid " << hv(14) << std::endl
              << "    gas eint: old_vector " << hv(6) << " writeback " << hv(7)
              << " t2inc(total E) " << hv(8) << "  T_solved " << hv(9)
              << "  v(old vector) " << hv(10) << " " << hv(11) << " " << hv(12)
              << std::endl;
  }
  // every cell with solved E <= 0, written-back eint <= 0 or old-vector eint <= 0 (up to
  // ncap per rank), with histograms by the GLOBAL theta / phi cell index
  const int ncap = 4096;
  const int nfld = 17;
  DvceArray2D<Real> lst("m1_t2dbg_lst", ncap, nfld);
  DvceArray1D<int> lcnt("m1_t2dbg_lcnt", 1);
  // the frozen 7-point row of each listed cell (M1_IW_TA..CKP, KB) with the neighbour
  // solutions, TR, TRHS, TDIA, SRCB, SRCR, KT and WCHI (= f_rr on the sp wedge)
  // + (m1-positivity) the implicit_vimp block of the row (x1 -2,-1,+1,+2 / x2 -2..+2 /
  // x3 -2..+2 couplings, JD, JRHS) when the block exists: 14 more
  const int nrow = 56;
  const int ivb = iw_vimp;
  auto dbr = dbrow;
  const bool hasd = (dbrow.extent_int(0) > 0);
  const bool bcg = bicg_on;
  DvceArray2D<Real> row("m1_t2dbg_row", ncap, nrow);
  DvceArray1D<Real> rmn("m1_t2dbg_rmn", indcs.nx1 + 2*indcs.ng);
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  if (hh) {
    par_for("m1_t2dbg_rmn", DevExeSpace(), is, ie, KOKKOS_LAMBDA(const int i) {
      Real sm = 0.0;
      for (int mm = 0; mm <= nmb1; ++mm) {
        for (int kk = ks; kk <= ke; ++kk) {
          for (int jj = js; jj <= je; ++jj) {sm += uh(mm,IDN,kk,jj,i);}
        }
      }
      rmn(i) = sm/static_cast<Real>((nmb1 + 1)*(ke - ks + 1)*(je - js + 1));
    });
  }
  Kokkos::parallel_for("m1_t2dbg_all", Kokkos::RangePolicy<DevExeSpace>(0, ntot),
  KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji;
    int r = idx - m*nkji;
    const int k = r/nji + ks;
    r -= (k - ks)*nji;
    const int j = r/ni + js;
    const int i = r - (j - js)*ni + is;
    const Real es = iw_(m,M1_IW_S2,k,j,i);
    const Real ego = hh ? iw_(m,M1_IW_EGN,k,j,i) : 1.0;
    const Real egw = gq ? (iw_(m,M1_IW_EGN,k,j,i) - (cl/ch)*(iw_(m,M1_IW_SRCR,k,j,i) -
                     iw_(m,M1_IW_SRCB,k,j,i)*iw_(m,M1_IW_EP,k,j,i))) : 1.0;
    const int flg = ((es > 0.0) ? 0 : 1) + ((egw > 0.0) ? 0 : 2) + ((ego > 0.0) ? 0 : 4);
    if (flg == 0) return;
    const int n = Kokkos::atomic_fetch_add(&lcnt(0), 1);
    if (n >= ncap) return;
    lst(n,0) = m; lst(n,1) = k; lst(n,2) = j; lst(n,3) = i; lst(n,4) = flg;
    lst(n,5) = x1v(m,i); lst(n,6) = x2v(m,j); lst(n,7) = x3v(m,k);
    lst(n,8) = hh ? uh(m,IDN,k,j,i) : 0.0;
    lst(n,9) = (hh && rmn(i) > 0.0) ? uh(m,IDN,k,j,i)/rmn(i) : 0.0;
    lst(n,10) = iw_(m,M1_IW_V1,k,j,i);
    lst(n,11) = iw_(m,M1_IW_V2,k,j,i);
    lst(n,12) = iw_(m,M1_IW_V3,k,j,i);
    lst(n,13) = es;
    lst(n,14) = iw_(m,M1_IW_EN,k,j,i);
    lst(n,15) = ego;
    lst(n,16) = egw;
    row(n,0) = iw_(m,M1_IW_TB,k,j,i);
    row(n,1) = iw_(m,M1_IW_TA,k,j,i);
    row(n,2) = iw_(m,M1_IW_S2,k,j,i-1);
    row(n,3) = iw_(m,M1_IW_TC,k,j,i);
    row(n,4) = iw_(m,M1_IW_S2,k,j,i+1);
    row(n,5) = bcg ? iw_(m,M1_IW_CJM,k,j,i) : 0.0;
    row(n,6) = iw_(m,M1_IW_S2,k,j-1,i);
    row(n,7) = bcg ? iw_(m,M1_IW_CJP,k,j,i) : 0.0;
    row(n,8) = iw_(m,M1_IW_S2,k,j+1,i);
    row(n,9) = bcg ? iw_(m,M1_IW_CKM,k,j,i) : 0.0;
    row(n,10) = iw_(m,M1_IW_S2,k-1,j,i);
    row(n,11) = bcg ? iw_(m,M1_IW_CKP,k,j,i) : 0.0;
    row(n,12) = iw_(m,M1_IW_S2,k+1,j,i);
    row(n,13) = bcg ? iw_(m,M1_IW_KB,k,j,i) : 0.0;
    row(n,14) = iw_(m,M1_IW_TR,k,j,i);
    row(n,15) = iw_(m,M1_IW_TRHS,k,j,i);
    row(n,16) = iw_(m,M1_IW_TDIA,k,j,i);
    row(n,17) = iw_(m,M1_IW_SRCB,k,j,i);
    row(n,18) = iw_(m,M1_IW_SRCR,k,j,i);
    row(n,19) = iw_(m,M1_IW_KT,k,j,i);
    row(n,20) = iw_(m,M1_IW_WCHI,k,j,i);
    row(n,21) = iw_(m,M1_IW_EP,k,j,i);
    row(n,22) = bcg ? iw_(m,M1_IW_KX,k,j,i) : 0.0;
    row(n,23) = iw_(m,M1_IW_G0,k,j,i);
    row(n,24) = iw_(m,M1_IW_DE0,k,j,i);
    row(n,25) = iw_(m,M1_IW_LRES,k,j,i);
    for (int q = 0; q < 16; ++q) {row(n,26+q) = hasd ? dbr(m,q,k,j,i) : 0.0;}
    for (int q = 0; q < 14; ++q) {
      row(n,42+q) = (ivb >= 0) ? iw_(m,ivb+M1_IV_X1M2+q,k,j,i) : 0.0;
    }
  });
  auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), lcnt);
  auto hl = Kokkos::create_mirror_view_and_copy(HostMemSpace(), lst);
  auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), row);
  const int nall = std::min(hc(0), ncap);
  auto &msz = pmy_pack->pmesh->mesh_size;
  const int gnx2 = pmy_pack->pmesh->mesh_indcs.nx2;
  const int gnx3 = pmy_pack->pmesh->mesh_indcs.nx3;
  const Real d2 = (msz.x2max - msz.x2min)/gnx2, d3 = (msz.x3max - msz.x3min)/gnx3;
  std::vector<int> hj(gnx2, 0), hk(gnx3, 0);
  std::cout << "  ALL rank " << rank << ": " << hc(0) << " cells (flag 1 E_solved<=0, "
            << "2 eint_writeback<=0, 4 eint_old<=0); listed " << nall << std::endl
            << "  # m k j i flag r/R_cm theta phi rho rho/<rho> v_r v_th v_ph E_solved "
            << "E_old eint_old eint_wb jg kg" << std::endl;
  for (int n = 0; n < nall; ++n) {
    const int jg = static_cast<int>(std::floor((hl(n,6) - msz.x2min)/d2));
    const int kg = static_cast<int>(std::floor((hl(n,7) - msz.x3min)/d3));
    if (jg >= 0 && jg < gnx2) hj[jg] += 1;
    if (kg >= 0 && kg < gnx3) hk[kg] += 1;
    std::cout << "  C";
    for (int f = 0; f < 5; ++f) std::cout << " " << static_cast<int>(hl(n,f));
    for (int f = 5; f < nfld; ++f) std::cout << " " << hl(n,f);
    std::cout << " " << jg << " " << kg << std::endl;
    if (static_cast<int>(hl(n,4)) & 1) {
      // the row: diag*E, each neighbour coefficient x neighbour solution, b; M-matrix
      // check: sum |off| / |diag| and the number of positive off-diagonal coefficients
      Real soff = 0.0, lhs = hr(n,0)*hl(n,13);
      int npos = 0;
      for (int c = 1; c <= 11; c += 2) {
        soff += std::fabs(hr(n,c));
        if (hr(n,c) > 0.0) npos += 1;
        lhs += hr(n,c)*hr(n,c+1);
      }
      std::cout << "  ROW TB " << hr(n,0) << " TB*E " << hr(n,0)*hl(n,13)
                << " | r- " << hr(n,1) << " x " << hr(n,2) << " r+ " << hr(n,3) << " x "
                << hr(n,4) << " | th- " << hr(n,5) << " x " << hr(n,6) << " th+ "
                << hr(n,7) << " x " << hr(n,8) << " | ph- " << hr(n,9) << " x "
                << hr(n,10) << " ph+ " << hr(n,11) << " x " << hr(n,12) << std::endl
                << "      b(KB) " << hr(n,13) << " lhs(solution) " << lhs
                << " TR " << hr(n,14) << " TRHS " << hr(n,15) << " TDIA " << hr(n,16)
                << " SRCB " << hr(n,17) << " SRCR " << hr(n,18) << " KT " << hr(n,19)
                << " f_rr " << hr(n,20) << " EP " << hr(n,21) << " KX " << hr(n,22)
                << " G0 " << hr(n,23) << " DE0 " << hr(n,24) << " LRES " << hr(n,25)
                << std::endl << "      M-matrix: sum|off|/|diag| "
                << soff/std::fabs(hr(n,0))
                << " positive off-diagonals " << npos << std::endl;
      if (hasd) {
        std::cout << "      TB terms: base(1+SRCB+TDIA) " << hr(n,31) << " diffusion "
                  << hr(n,26) << " opac_newton " << hr(n,27) << " advection " << hr(n,28)
                  << " vimp_JD " << hr(n,30) << " other " << hr(n,32) << std::endl
                  << "      TR terms: base(EN+SRCR+TRHS) " << hr(n,33) << " curvature "
                  << hr(n,36) << " face_F0 " << hr(n,37) << " opac_newton " << hr(n,38)
                  << " enthalpy_corr " << hr(n,39) << " vimp_JRHS " << hr(n,40)
                  << " other " << hr(n,41) << std::endl;
      }
      if (ivb >= 0) {
        std::cout << "      VIMP x1(-2,+2) " << hr(n,42) << " " << hr(n,43)
                  << " x2(-2,-1,+1,+2) " << hr(n,44) << " " << hr(n,45) << " "
                  << hr(n,46) << " " << hr(n,47) << " x3(-2,-1,+1,+2) " << hr(n,48)
                  << " " << hr(n,49) << " " << hr(n,50) << " " << hr(n,51)
                  << " JD " << hr(n,52) << " J1M " << hr(n,53) << " J1P " << hr(n,54)
                  << " JRHS " << hr(n,55) << std::endl;
      }
    }
  }
  std::cout << "  HIST_J rank " << rank << " (global theta index: count)";
  for (int q = 0; q < gnx2; ++q) {
    if (hj[q] > 0) {std::cout << " " << q << ":" << hj[q];}
  }
  std::cout << std::endl << "  HIST_K rank " << rank << " (global phi index: count)";
  for (int q = 0; q < gnx3; ++q) {
    if (hk[q] > 0) {std::cout << " " << q << ":" << hk[q];}
  }
  std::cout << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitMusclBuild
//! \brief implicit_hr_recon = plm (xthinfix-1009 Fix B), once per Picard pass after the
//! face coefficients (x1: ifw, x2/x3: ifw2/ifw3) of the pass are set.
//!  (1) on the first implicit_hr_recon_fresh (2) passes of the step, then FROZEN for the
//!      rest of it (re-made every pass, the limiter switches on a few cells and the
//!      Picard loop limit-cycles: cyl, 200 passes per step; made once from the step-start
//!      E, it is not TVD for a front that crosses many cells in the step: E <= 0 in
//!      vacuum, xb20), sig_d per cell from the Picard iterate E (M1_IW_EP), zero where
//!      the positivity fallback has set KILL this step: van Leer, sig = 4ab/(a+b)^2 for
//!      ab > 0 (a, b the two one-sided differences), else 0, so the face value
//!      E_c +- s_c/2, s_c = sig_c (E_c+1 - E_c-1)/2, lies between E_c and the neighbour
//!      at the lagged state (TVD).  sig = 0 next to a physical boundary (its ghost is not
//!      filled) and on an efix row.  Halo-exchanged (3 components).
//!  (2) the row: per direction and face, dG = (HCL s_L - HCR s_R)/2 with the frozen sig,
//!      row = f+ dG+ - f- dG- (f = dt/dx_d, or dt A_f/V on sp, as the 7-point row's
//!      HCL/HCR terms), as coefficients on E_c-2..E_c+2.  Linear in E for the pass, so
//!      the Krylov operator is the 7-point row + M1MusclRow and the line/mg
//!      preconditioners keep the 7-point (dc) part.  HCL = HCR = 0 at physical faces and
//!      w = 0 faces, so a central face is untouched.

void RadiationM1::ImplicitMusclBuild() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto ifw_ = ifw;
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  auto pos_ = part_pos;
  const int nblkx1 = part_nblk;
  const bool thrd = trans_x3;
  const bool twod = pmy_pack->pmesh->multi_d;
  const Real dt = dt_sub;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const bool sph = sph_geom;
  auto cvol = pmy_pack->pcoord->volume;
  auto carea = pmy_pack->pcoord->area;
  const int b = iw_muscl;
  const bool l2 = twod && (bw2_.extent_int(0) > 0);
  const bool l3 = thrd && (bw3_.extent_int(0) > 0);
  // (1) the frozen limiter (once per step)
  if (muscl_nbuild < impl_muscl_nfresh) {
  const bool first = (muscl_nbuild == 0);
  muscl_nbuild += 1;
  par_for("m1_muscl_sig", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int ipos = pos_.d_view(m);
    const bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    const bool efx = !cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                                 (i == ie && topb && bchi == M1_IBC_EFIX));
    for (int d = 0; d < 3; ++d) {
      Real sg = 0.0;
      bool ok = !efx && ((d == 0) || (d == 1 && l2) || (d == 2 && l3));
      if (ok && d == 0 && !cyclic) {
        ok = !((i == is && botb) || (i == ie && topb));
      } else if (ok && d > 0) {
        const BoundaryFlag blo = mbbcs(m, (d == 1) ? BoundaryFace::inner_x2
                                                   : BoundaryFace::inner_x3);
        const BoundaryFlag bhi = mbbcs(m, (d == 1) ? BoundaryFace::outer_x2
                                                   : BoundaryFace::outer_x3);
        const bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
        const bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
        const int c = (d == 1) ? j : k;
        const int cs = (d == 1) ? js : ks, ce = (d == 1) ? je : ke;
        ok = !((c == cs && plo) || (c == ce && phi));
      }
      if (ok) {
        int km = k, kp = k, jm = j, jp = j, im = i, ip = i;
        if (d == 0) {
          im = (i > is) ? (i-1) : ie;
          ip = (i < ie) ? (i+1) : is;
        } else if (d == 1) {
          jm = j - 1; jp = j + 1;
        } else {
          km = k - 1; kp = k + 1;
        }
        const Real e0 = iw_(m,M1_IW_EP,k,j,i);
        const Real ea = e0 - iw_(m,M1_IW_EP,km,jm,im);
        const Real eb = iw_(m,M1_IW_EP,kp,jp,ip) - e0;
        if (ea*eb > 0.0) {sg = 4.0*ea*eb/SQR(ea + eb);}
      }
      if (first) {
        if (d == 0) {iw_(m,b+M1_IM_KILL,k,j,i) = 0.0;}
      } else if (iw_(m,b+M1_IM_KILL,k,j,i) > 0.5) {
        sg = 0.0;
      }
      iw_(m,b+M1_IM_SIG+d,k,j,i) = sg;
    }
  });
  ImplicitHaloExchange(3, b + M1_IM_SIG);
  }
  // (2) the row
  par_for("m1_muscl_row", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dg = 0.0;
    Real cf[3][4];
    for (int d = 0; d < 3; ++d) {
      for (int q = 0; q < 4; ++q) {cf[d][q] = 0.0;}
      if ((d == 1 && !l2) || (d == 2 && !l3)) {continue;}
      int km = k, kp = k, jm = j, jp = j, im = i, ip = i;
      if (d == 0) {
        im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
        ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
      } else if (d == 1) {
        jm = j - 1; jp = j + 1;
      } else {
        km = k - 1; kp = k + 1;
      }
      // the two faces of the cell along d: HCL/HCR, and the row factors f+ f-
      Real hlp, hrp, hlm, hrm, fp, fm;
      if (d == 0) {
        hlp = ifw_(m,M1_IFW_HCL,k,j,i+1); hrp = ifw_(m,M1_IFW_HCR,k,j,i+1);
        hlm = ifw_(m,M1_IFW_HCL,k,j,i);   hrm = ifw_(m,M1_IFW_HCR,k,j,i);
        fp = dt/mbsize(m).dx1; fm = fp;
        if (sph) {
          const Real iv = dt/cvol(m,k,j,i);
          fp = carea.x1f(m,k,j,i+1)*iv; fm = carea.x1f(m,k,j,i)*iv;
        }
      } else if (d == 1) {
        hlp = bw2_(m,M1_IFW_HCL,k,j+1,i); hrp = bw2_(m,M1_IFW_HCR,k,j+1,i);
        hlm = bw2_(m,M1_IFW_HCL,k,j,i);   hrm = bw2_(m,M1_IFW_HCR,k,j,i);
        fp = dt/mbsize(m).dx2; fm = fp;
        if (sph) {
          const Real iv = dt/cvol(m,k,j,i);
          fp = carea.x2f(m,k,j+1,i)*iv; fm = carea.x2f(m,k,j,i)*iv;
        }
      } else {
        hlp = bw3_(m,M1_IFW_HCL,k+1,j,i); hrp = bw3_(m,M1_IFW_HCR,k+1,j,i);
        hlm = bw3_(m,M1_IFW_HCL,k,j,i);   hrm = bw3_(m,M1_IFW_HCR,k,j,i);
        fp = dt/mbsize(m).dx3; fm = fp;
        if (sph) {
          const Real iv = dt/cvol(m,k,j,i);
          fp = carea.x3f(m,k+1,j,i)*iv; fm = carea.x3f(m,k,j,i)*iv;
        }
      }
      const Real sc = iw_(m,b+M1_IM_SIG+d,k,j,i);
      const Real sp = iw_(m,b+M1_IM_SIG+d,kp,jp,ip);
      const Real sm = iw_(m,b+M1_IM_SIG+d,km,jm,im);
      // dG+ = 1/4 [HCL+ sc (x+1 - x-1) - HCR+ sp (x+2 - x0)]
      // dG- = 1/4 [HCL- sm (x0 - x-2) - HCR- sc (x+1 - x-1)]
      // a zero face coefficient (physical or central face) never multiplies a ghost sig
      const Real ca = (hlp != 0.0) ? 0.25*fp*hlp*sc : 0.0;
      const Real cb = (hrp != 0.0) ? 0.25*fp*hrp*sp : 0.0;
      const Real cc = (hlm != 0.0) ? 0.25*fm*hlm*sm : 0.0;
      const Real cd = (hrm != 0.0) ? 0.25*fm*hrm*sc : 0.0;
      // row = dG+ - dG- : x+2: -cb; x+1: ca + cd; x0: cb - cc; x-1: -(ca + cd); x-2: cc
      cf[d][0] = cc;
      cf[d][1] = -(ca + cd);
      cf[d][2] = ca + cd;
      cf[d][3] = -cb;
      dg += cb - cc;
    }
    iw_(m,b+M1_IM_D,k,j,i) = dg;
    for (int q = 0; q < 4; ++q) {
      iw_(m,b+M1_IM_X1+q,k,j,i) = cf[0][q];
      iw_(m,b+M1_IM_X2+q,k,j,i) = cf[1][q];
      iw_(m,b+M1_IM_X3+q,k,j,i) = cf[2][q];
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitMusclKill
//! \brief implicit_hr_recon = plm POSITIVITY (xthinfix-1009 Fix B): the plm row is not an
//! M-matrix (the E_c-2 / E_c+2 coefficients of the inflow faces are positive), so a
//! solve may give E <= 0.  Then sig = 0 (donor cell) in every cell within 2 cells, along
//! each axis, of a cell the pass solved to E <= 0 (M1_IW_S2), for the rest of the step;
//! the next pass rebuilds the row.  Elsewhere the faces keep plm.

void RadiationM1::ImplicitMusclKill() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const bool twod = pmy_pack->pmesh->multi_d;
  const bool thrd = trans_x3;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int b = iw_muscl;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  ImplicitHaloExchange(1, M1_IW_S2);
  par_for("m1_muscl_kill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // physical ghosts are not filled: read only block/periodic ones
    auto open = [&](const BoundaryFace f) {
      const BoundaryFlag q = mbbcs(m,f);
      return (q == BoundaryFlag::block) || (q == BoundaryFlag::periodic);
    };
    const bool o2l = open(BoundaryFace::inner_x2), o2h = open(BoundaryFace::outer_x2);
    const bool o3l = open(BoundaryFace::inner_x3), o3h = open(BoundaryFace::outer_x3);
    bool neg = false;
    for (int o = -2; o <= 2; ++o) {
      int ii = i + o;
      if (cyclic) {
        const int n = ie - is + 1;
        while (ii < is) {ii += n;}
        while (ii > ie) {ii -= n;}
      } else {
        ii = (ii < is) ? is : ((ii > ie) ? ie : ii);
      }
      neg = neg || !(iw_(m,M1_IW_S2,k,j,ii) > 0.0);
      const int jj = j + o, kk = k + o;
      if (twod && (jj >= js || o2l) && (jj <= je || o2h)) {
        neg = neg || !(iw_(m,M1_IW_S2,k,jj,i) > 0.0);
      }
      if (thrd && (kk >= ks || o3l) && (kk <= ke || o3h)) {
        neg = neg || !(iw_(m,M1_IW_S2,kk,j,i) > 0.0);
      }
    }
    if (neg) {
      for (int d = 0; d < 3; ++d) {iw_(m,b+M1_IM_SIG+d,k,j,i) = 0.0;}
      iw_(m,b+M1_IM_KILL,k,j,i) = 1.0;
    }
  });
  ImplicitHaloExchange(3, b + M1_IM_SIG);
}

} // namespace radm1
