#ifndef UTILS_TWO_STREAM_COLUMN_CK_HPP_
#define UTILS_TWO_STREAM_COLUMN_CK_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_column_ck.hpp
//! \brief problem/ck_implicit: a BACKWARD-EULER column solve for the CORRELATED-K
//! two-stream (deep_hot_jupiter_rt).  Default OFF, in which case nothing in this file is
//! allocated, launched or read and the answer is bitwise today's.
//!
//! WHY.  What the gas receives on the correlated-k path today is NOT the flux divergence
//! the sweep computed.  The apply block relaxes each cell separately toward the
//! equilibrium of ONE sweep -- absorbed field frozen -- and then clips the step with
//! rt_de_max and the floors.  The existing rt_desum diagnostic measures the gap:
//! tests_ck_sph/README.md 0.3 reports -2.9 .. -5.9 % with ck_spherical on and
//! -36.8 .. -40.5 % with it off, on the production radial grid.  A flux-form scheme whose
//! face fluxes telescope exactly (that README 0.1) still loses that fraction of the
//! energy in transients, because the applied de is not the divergence of any flux.
//!
//! WHAT THIS DOES.  Per column, unknown = the cell internal energy e_i, the balance is
//!     R_i(e) = e_i - e*_i - bdt * S_i(e) = 0,
//! S_i the FULL correlated-k source the apply block already forms: the (1-w)-weighted
//! two-stream flux divergence, the tau-blend handover, and the direct beam.  R is
//! evaluated by RE-RUNNING THE ORDINARY SWEEP at the current iterate -- one ck sweep per
//! Newton pass -- so the fixed point is the exact backward-Euler balance whatever the
//! Jacobian gets wrong.  At convergence the gas is handed EXACTLY bdt*S(e^{n+1}) with no
//! relaxation factor, no rt_de_max clip and no per-cell damping: flux-form conservation
//! to the Newton tolerance.
//!
//! THE JACOBIAN.  At frozen opacity the thermal exchange is LINEAR in the band Planck
//! functions B_b(T), so dS_i/dT_j = sum_{bands,g} w_g (dSrc_i/dB_b,j)(dB_b/dT)_j.  The
//! NEAREST-NEIGHBOUR part of dSrc_i/dB is accumulated inside the ck sweep itself out of
//! the layer coefficients already in registers (see "ck_implicit" in the rt_chain_ck
//! kernel); everything two cells away or further, and the direct beam -- which at frozen
//! opacity is temperature-INDEPENDENT and therefore has no Jacobian term at all -- stay
//! in the residual.  Rows:
//!     J_{i,i-1} = -bdt jac0_i/cv_{i-1},  J_{ii} = 1 - bdt jac1_i/cv_i,
//!     J_{i,i+1} = -bdt jac2_i/cv_{i+1}
//! with cv = de/dT.  jac1 < 0 (the cell's own emission) and jac0, jac2 >= 0 (what the
//! neighbours send it), so J is an M-matrix for any bdt and the Thomas sweep is stable
//! without pivoting.
//!
//! WHAT IS NOT IN THE JACOBIAN, deliberately, and why it costs nothing at the fixed
//! point: (i) the direct beam, exactly temperature-independent here; (ii) the
//! ck_spherical
//! PROBE passes -- the two extra passes that build the face-mixing constants Cmx are
//! re-run every Newton pass, so they are exact in R, but their B-dependence is dropped
//! from J (the face-mixing FACTORS (1+beta), (1-beta) and the up ray's area frame change
//! ARE carried, since they are free); (iii) the tau-blend handover d[wF]/dr; (iv) the
//! opacity's own temperature dependence, which is frozen over the step by construction.
//! All four make this a quasi-Newton, not a Newton: the iteration converges linearly
//! rather than quadratically where they matter, and lands on the same root.
//!
//! WHAT IS NOT FOLDED IN.  The radial radiative CONDUCTION below/around ck_pcut_bar
//! (rad_tau_lo/hi) is NOT folded into this tridiagonal -- unlike the grey
//! rt_implicit_column = 3, which merges it into Conduction::ImplicitRadialUpdate.  It
//! stays the separate operator it is today, so the 10 bar hand-over remains split: the
//! blend's handover flux enters S_i explicitly (lagged by one operator), exactly as it
//! does with ck_implicit off.  Folding it in needs the conduction tridiagonal to accept a
//! correlated-k row, which is a second deliverable.
//!
//! HEAT CAPACITY.  cv = e/T, i.e. e ~ T -- the SAME approximation the semi-implicit
//! apply makes (see rt_semi_lin and the e_eq = e (A/E)^(1/4) form).  It is poor where H2
//! is partly dissociated, but it appears ONLY in the Jacobian: it changes how many passes
//! the Newton takes and not what it converges to.
//!
//! CADENCE.  This solves the balance over the interval bdt it is handed, which on this
//! pgen is one RK STAGE (the source term is called per stage; with rk2 that is two ck
//! sweeps per hydro step today).  ck_implicit does not change that: it makes each stage's
//! source implicit, at ck_impl_maxit sweeps per stage instead of one.

#include <math.h>
#include <cstdio>

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"

namespace two_stream_rt {

//----------------------------------------------------------------------------------------
// problem/ck_implicit: the switch.  false = bitwise the semi-implicit per-cell apply.
inline bool ck_implicit = false;
// problem/ck_impl_tol: the RESIDUAL tolerance, max_i |R_i|/(e_i + eps e_max) over the
// column set, with eps = ck_impl_norm_eps.  Measuring a thin top cell against its own
// energy demands an absolute precision nothing else in the scheme has (the same reason
// rt_impl_norm exists on the grey path), so it is measured against its share instead.
inline Real ck_impl_tol = 1.0e-8;
inline Real ck_impl_norm_eps = 1.0e-3;
// problem/ck_impl_dtol: also stop when the Newton STEP is small, max_i |de_i|/e_i.  Both
// tests must pass.
inline Real ck_impl_dtol = 1.0e-8;
// problem/ck_impl_maxit: the pass cap.  Each pass is one full ck sweep.
inline int ck_impl_maxit = 8;
// problem/ck_impl_dtmax: cap on |de|/e per Newton pass.  A belt, not a limiter: it is
// applied to the STEP, so a capped pass is simply a shorter one and the fixed point is
// unchanged.  A column that needed it is counted.
inline Real ck_impl_dtmax = 0.5;
// problem/ck_impl_demax: bound on the TOTAL excursion from e^n over the step,
// |e - e*| <= ck_impl_demax e*.  This is the same bound rt_de_max (default 0.5) imposes
// on the semi-implicit apply today, and it exists for one measured reason: with
// problem/ck_spherical = FALSE the plane-parallel correlated-k source has NO
// positive-energy backward-Euler root in the thin top cells of the production grid (the
// undiluted A F manufactures a cooling rate ~13x the cell's own internal energy per
// stage that is almost independent of T, so the Newton runs the cell to e < 0).  The
// Newton then does not converge, and this bound keeps the answer no worse than today's
// while the fallback counter says so.  With ck_spherical = true it is never reached --
// measured, capped = 0 on every gate run.  0 = no bound.
inline Real ck_impl_demax = 0.5;
// problem/ck_impl_verbose: print the per-call pass count and residual.
inline bool ck_impl_verbose = false;
// problem/ck_impl_debug: print at most this many CAPPED cells per pass, with the row
// that produced them.  Diagnostic only.
inline int ck_impl_debug = 0;
// problem/ck_impl_refresh_kappa: re-look-up the correlated-k opacity (and hence the beam
// transmission tau_ray) at every Newton pass instead of freezing it over the step.  OFF
// by default.  Cost: the pre-opacity kernel is a ck_continuum call plus two table index
// look-ups per cell and the beam's tau_ray is rebuilt inside the sweep anyway, so a
// refresh pass costs one extra rt_pre_opac launch, ~5 % of a sweep -- but it makes the
// Jacobian inconsistent with the residual (dkappa/dT is nowhere in J), so it converges
// more slowly.  Provided because "frozen opacity" is an assumption, not a law.
inline bool ck_impl_refresh_kappa = false;

// the Newton pass index inside one RT call; read by the pass function to decide whether
// the opacity is rebuilt.  -1 = ck_implicit is off.
inline int ck_impl_pass = -1;

// dSrc_i/dT_{i-1,i,i+1}, (m, 3, k, j, i), summed over bands and g-points
inline DvceArray5D<Real> *ck_jac_ptr = nullptr;
// dB_b/dT at the cell, (m, CK_NB, i, k, j) -- the Bb_g layout exactly
inline DvceArray5D<Real> *ck_dbdt_ptr = nullptr;
// the full source S_i the apply block formed, and the internal energy it formed it at
inline DvceArray4D<Real> *ck_src_ptr = nullptr;
inline DvceArray4D<Real> *ck_ei_ptr = nullptr;
// e*_i, the internal energy the step started from
inline DvceArray4D<Real> *ck_estar_ptr = nullptr;
// the Thomas sweep's two work rows, (m,k,j,i)
inline DvceArray4D<Real> *ck_cp_ptr = nullptr;
inline DvceArray4D<Real> *ck_dp_ptr = nullptr;
// 0 = max |R|/(e + eps e_max), 1 = max |de|/e, 2 = capped-step count, 3 = fallback
// count, 4 = sum (e - e*) dx, 5 = sum bdt S dx.  Slots 4 and 5 are the DIRECT analogue
// of the rt_desum line -- what the gas received against what the sweep asked for -- and
// their ratio is the number tests_ck_sph/README.md 0.3 reports for the semi-implicit
// apply.  Here it is 1 to the Newton tolerance by construction.
inline DvceArray1D<Real> *ck_conv_ptr = nullptr;
inline Real ck_impl_last_gap = 0.0;
// host-side record of the last call, for the report line
inline int ck_impl_last_it = 0;
inline Real ck_impl_last_res = 0.0;
inline Real ck_impl_last_dst = 0.0;
inline int ck_impl_ncap = 0;
inline int ck_impl_nfall = 0;
inline int ck_impl_nonconv = 0;

//----------------------------------------------------------------------------------------
//! \fn void CkImplAlloc
//! \brief allocate the column-solve scratch on the first call.  Never reached with
//! ck_implicit off.

inline void CkImplAlloc(const int nmb, const int nb, const int n1, const int n2,
                        const int n3) {
  if (ck_jac_ptr != nullptr &&
      ck_jac_ptr->extent(0) == static_cast<size_t>(nmb) &&
      ck_jac_ptr->extent(4) == static_cast<size_t>(n1)) {
    return;
  }
  if (ck_jac_ptr != nullptr) {
    delete ck_jac_ptr;
    delete ck_dbdt_ptr;
    delete ck_src_ptr;
    delete ck_ei_ptr;
    delete ck_estar_ptr;
    delete ck_cp_ptr;
    delete ck_dp_ptr;
  }
  ck_jac_ptr = new DvceArray5D<Real>("ck_jac", nmb, 3, n3, n2, n1);
  ck_dbdt_ptr = new DvceArray5D<Real>("ck_dbdt", nmb, nb, n1, n3, n2);
  ck_src_ptr = new DvceArray4D<Real>("ck_src", nmb, n3, n2, n1);
  ck_ei_ptr = new DvceArray4D<Real>("ck_ei", nmb, n3, n2, n1);
  ck_estar_ptr = new DvceArray4D<Real>("ck_estar", nmb, n3, n2, n1);
  ck_cp_ptr = new DvceArray4D<Real>("ck_cp", nmb, n3, n2, n1);
  ck_dp_ptr = new DvceArray4D<Real>("ck_dp", nmb, n3, n2, n1);
  if (ck_conv_ptr == nullptr) {
    ck_conv_ptr = new DvceArray1D<Real>("ck_conv", 8);
  }
}

//----------------------------------------------------------------------------------------
//! \fn int CkImplStep
//! \brief form the residual of the CURRENT iterate, test it, and -- if it is not yet
//! converged -- solve the tridiagonal and apply the step to u0.  Returns 1 if a step was
//! taken, 0 if the iterate was already converged (in which case NOTHING is applied and
//! the state handed to the gas satisfies |R| <= tol).
//!
//! The Thomas sweep runs over the cells the correlated-k domain owns, icut .. ie.  Cells
//! below the cut receive no radiative source at all and are not in the system.
//!
//! FALLBACK.  A row whose diagonal is not positive, or a solve that produces a
//! non-finite step, drops that CELL back to the bounded relaxation the semi-implicit
//! apply would have taken, de = deq (1 - exp(-x)) with deq = e((A/E)^(1/4) - 1) -- here
//! with A and E built from the same source, so it degenerates to the explicit step where
//! the emission is unavailable.  It is counted in slot 3.

inline int CkImplStep(Mesh *pm, DvceArray5D<Real> u0, DvceArray3D<int> icut_,
                      DvceArray4D<Real> T_, DvceArray4D<Real> dx1_, const Real bdt) {
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto jac_ = *ck_jac_ptr;
  auto src_ = *ck_src_ptr;
  auto ei_ = *ck_ei_ptr;
  auto est_ = *ck_estar_ptr;
  auto cnv_ = *ck_conv_ptr;
  auto cp_ = *ck_cp_ptr;
  auto dp_ = *ck_dp_ptr;
  Kokkos::deep_copy(cnv_, 0.0);
  const Real tol = ck_impl_tol;
  const Real dtol = ck_impl_dtol;
  const Real eps = ck_impl_norm_eps;
  const Real dcap = ck_impl_dtmax;
  const Real detot = ck_impl_demax;
  const int ndbg = ck_impl_debug;

  // ---- pass 1: the residual norm of the CURRENT iterate --------------------------
  par_for("ck_impl_res", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    Real emax = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real e = ei_(m,k,j,i);
      if (e > emax) emax = e;
    }
    if (!(emax > 0.0)) return;
    Real rn = 0.0;
    Real sg = 0.0, ss = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
      const Real s = fabs(r)/(ei_(m,k,j,i) + eps*emax);
      if (s > rn) rn = s;
      const Real dx = dx1_(m,k,j,i);
      sg += (ei_(m,k,j,i) - est_(m,k,j,i))*dx;
      ss += bdt*src_(m,k,j,i)*dx;
    }
    Kokkos::atomic_max(&cnv_(0), rn);
    Kokkos::atomic_add(&cnv_(4), sg);
    Kokkos::atomic_add(&cnv_(5), ss);
  });
  auto hc = Kokkos::create_mirror_view(cnv_);
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_res = hc(0);
  ck_impl_last_gap = (hc(5) != 0.0) ? (hc(4)/hc(5) - 1.0) : 0.0;
  if (hc(0) <= tol) return 0;

  // ---- pass 2: assemble, solve, apply --------------------------------------------
  par_for("ck_impl_tri", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    // the Thomas sweep.  One thread owns the column, exactly as the sweeps do; the two
    // work rows live in device arrays rather than on the stack so the kernel carries no
    // per-thread footprint that would keep it off an accelerator.
    const int n = ie - ic + 1;
    if (n < 1) return;
    bool bad = false;
    for (int q=0; q<n; ++q) {
      const int i = ic + q;
      const Real ei = ei_(m,k,j,i);
      const Real Ti = T_(m,k,j,i);
      if (!(ei > 0.0) || !(Ti > 0.0)) {
        bad = true;
        break;
      }
      const Real cvi = ei/Ti;
      const Real eim = (q > 0) ? ei_(m,k,j,i-1) : 0.0;
      const Real Tim = (q > 0) ? T_(m,k,j,i-1) : 1.0;
      const Real eip = (q < n-1) ? ei_(m,k,j,i+1) : 0.0;
      const Real Tip = (q < n-1) ? T_(m,k,j,i+1) : 1.0;
      const Real cvm = (q > 0 && Tim > 0.0 && eim > 0.0) ? (eim/Tim) : 1.0;
      const Real cvp = (q < n-1 && Tip > 0.0 && eip > 0.0) ? (eip/Tip) : 1.0;
      const Real a = (q > 0) ? (-bdt*jac_(m,0,k,j,i)/cvm) : 0.0;
      const Real b = 1.0 - bdt*jac_(m,1,k,j,i)/cvi;
      const Real c = (q < n-1) ? (-bdt*jac_(m,2,k,j,i)/cvp) : 0.0;
      const Real d = -(ei - est_(m,k,j,i) - bdt*src_(m,k,j,i));
      if (!(b > 0.0)) {
        bad = true;
        break;
      }
      const Real den = b - a*((q > 0) ? cp_(m,k,j,i-1) : 0.0);
      if (!(fabs(den) > 0.0)) {
        bad = true;
        break;
      }
      cp_(m,k,j,i) = c/den;
      dp_(m,k,j,i) = (d - a*((q > 0) ? dp_(m,k,j,i-1) : 0.0))/den;
    }
    if (bad) {
      Kokkos::atomic_add(&cnv_(3), 1.0);
      // the bounded explicit fallback for the whole column: this is what the
      // semi-implicit apply would have deposited with its relaxation factor set to 1,
      // capped at ck_impl_dtmax.  It keeps the column finite; it does not converge.
      for (int i=ic; i<ie+1; ++i) {
        const Real ei = ei_(m,k,j,i);
        if (!(ei > 0.0)) continue;
        Real de = est_(m,k,j,i) + bdt*src_(m,k,j,i) - ei;
        const Real lim = dcap*ei;
        if (de > lim) de = lim;
        if (de < -lim) de = -lim;
        u0(m,IEN,k,j,i) += de;
      }
      return;
    }
    Real dn = 0.0;
    Real prev = 0.0;
    bool cap = false;
    for (int q=n-1; q>=0; --q) {
      const int i = ic + q;
      Real de = dp_(m,k,j,i) - cp_(m,k,j,i)*prev;
      // NaN guard; <cmath> is not available on the device side
      if (!(de == de)) de = 0.0;
      prev = de;
      const Real ei = ei_(m,k,j,i);
      const Real lim = dcap*ei;
      if (ndbg > 0 && fabs(de) > 0.999*lim &&
          Kokkos::atomic_fetch_add(&cnv_(6), 1.0) < static_cast<Real>(ndbg)) {
        Kokkos::printf("### ck_impl_cap m=%d k=%d j=%d i=%d e=%.6e T=%.6e "
                       "est=%.6e src=%.6e R=%.6e j0=%.6e j1=%.6e j2=%.6e de=%.6e\n",
                       m, k, j, i, ei, T_(m,k,j,i), est_(m,k,j,i), src_(m,k,j,i),
                       ei - est_(m,k,j,i) - bdt*src_(m,k,j,i), jac_(m,0,k,j,i),
                       jac_(m,1,k,j,i), jac_(m,2,k,j,i), de);
      }
      if (de > lim) {
        de = lim;
        cap = true;
      }
      if (de < -lim) {
        de = -lim;
        cap = true;
      }
      // the total-excursion bound, see ck_impl_demax
      const Real es = est_(m,k,j,i);
      if (detot > 0.0 && es > 0.0) {
        Real en = ei + de;
        const Real ehi = es*(1.0 + detot);
        const Real elo = es*(1.0 - detot);
        if (en > ehi) {
          en = ehi;
          cap = true;
        }
        if (en < elo) {
          en = elo;
          cap = true;
        }
        de = en - ei;
      }
      u0(m,IEN,k,j,i) += de;
      const Real s = (ei > 0.0) ? fabs(de)/ei : 0.0;
      if (s > dn) dn = s;
    }
    Kokkos::atomic_max(&cnv_(1), dn);
    if (cap) Kokkos::atomic_add(&cnv_(2), 1.0);
  });
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_dst = hc(1);
  ck_impl_ncap = static_cast<int>(hc(2));
  ck_impl_nfall = static_cast<int>(hc(3));
  // the step test.  The residual of the NEW iterate is measured at the top of the next
  // pass, which is the only honest place to measure it.
  if (hc(1) <= dtol && hc(0) <= tol) return 0;
  return 1;
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_COLUMN_CK_HPP_
