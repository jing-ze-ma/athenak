//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_rhs.cpp
//! \brief the right-hand side of the E solve: exchange, work, face sources.
//!
//!  called from: rad_vet_tasks.cpp
//!  calls into: rad_vet.cpp, rad_vet_bvals.cpp, rad_vet_flux_seam.cpp, rad_vet_init.cpp

#include <cmath>
#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

// 4-D par_for (same index flattening) at MB blocks of 256 per SM
template <int MB, typename F>
static void CleanParForLb(const char *name, const int nmb1, const int ks, const int ke,
                          const int js, const int je, const int is, const int ie,
                          const F &f) {
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int mb = 2;   // HIP: cln_rhs at 2 waves per SIMD (faster despite spills)
#else
  constexpr int mb = CleanLbMb(MB);
#endif
  Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,
                       mb>>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/ni;
    f(m, ks + k, js + j, is + idx - m*nkji - k*nji - j*ni);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagAdvFace
//! \brief the advective flux a E at face (k,j,i) normal to n as the face fluxes carry it:
//! clean_stag_upwind takes the share phi as donor-cell a_f E_f plus the lagged rest
//! (CleanFaceWeights), else the central mean; xh the face momentum state

KOKKOS_INLINE_FUNCTION
static Real CleanStagAdvFace(const CleanArgs &x, const CleanArgs &xh, const int n,
                             const int m, const int kk, const int jj, const int ii) {
  const int kl = kk - (n == 2), jl = jj - (n == 1), il = ii - (n == 0);
  const Real ul = x.w(m,CLW_A1+n,kl,jl,il), ur = x.w(m,CLW_A1+n,kk,jj,ii);
  const Real el = x.w(m,CLW_EK,kl,jl,il), er = x.w(m,CLW_EK,kk,jj,ii);
  const Real ac = 0.5*(ul*el + ur*er);
  if (!x.upw || CleanBoundaryFace(x, n, m, kk, jj, ii)) return ac;
  Real al, ar, phi = 0.0;
  CleanFaceWeights(xh, n, m, kk, jj, ii, al, ar, &phi);
  if (!(phi > 0.0)) return ac;
  const Real ah = CleanStagAdvectCorr(x, n, m, kk, jj, ii) +
                  ((ul + ur > 0.0) ? ul*el : ur*er);
  return ac + phi*(ah - ac);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRightHandSide
//! \brief B = E_old + h c kP a T^4 - r h c kE E_k - W - h div F_source + h esrc, the
//! exchange linearised about the iterate (r = s/(d cv + s), s = 4 h c kP a T^3), and the
//! initial guess X = E_k.  W: the work of the iterate's radiation force on the gas (with
//! wb_arad also the reference work h d arad_ref.v); esrc: the explicit energy source.
//! With gas coupling, the beta terms of Ma et al. (2025, eq. 56): emission takes
//! G = 1 - beta^2/2 (CLW_GM) and B adds the lagged hc kE [(beta^2/2 - beta.a/c) E_k +
//! 2 beta.y_k] (y_k = F/c of the iterate; the gas eint takes its opposite, CLW_EP), with
//! beta = v/c (CLW_BT) and hc kP a T^4 (CLW_RE) for the face fluxes (eq. 57).

void RadiationVet::CleanRightHandSide(const bool beta, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const Real arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  // wb_arad: the reference work rides the row (force_reference_work = split)
  const int fref = force_ref;
  const bool wsp = (fref == VET_FREF_WB_ARAD || fref == VET_FREF_WB_MD) && !x.vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  // backward-Euler step, FrefWaccOn: E pays exactly the reference work the gas got in the
  // hydro stages, its share h/dt of fref_wacc
  const bool wpa = wsp && FrefWaccOn();
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : dt_sub/pmy_pack->pmesh->dt;
  // explicit energy source esrc: h esrc joins the old E (the slope of E carries it)
  const bool eso = esrc_on;
  auto es = esrc;
  const bool bud = cln_bon;
  // CleanSteadyStart's equilibrium relax (icr_fmeas): h -> infinity
  const bool eqr = icr_fmeas && x.vzero && x.feedback;

  if (x.feedback) {
    const bool rx = x.cub;
    CLEAN_HIPX(const CleanArgs *pxb = CleanDevArgs(x, 0);)
    CleanFor<0, 3>("cln_beta", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *pxb;)   // in device memory: no per-thread copy
      Real p[3];
      CleanStagGasMomentum(x, uh, m, k, j, i, p);
      for (int n = 0; n < 3; ++n) {
        x.w(m,CLW_BT+n,k,j,i) = (p[n] + x.w(m,CLW_DM+n,k,j,i))/(uh(m,IDN,k,j,i)*x.c);
      }
      x.w(m,CLW_RE,k,j,i) = x.h*x.c*op(m,VET_OP_P,k,j,i)*arad_*
                            SQR(SQR(x.w(m,CLW_T,k,j,i)));
      // r of the exchange, for the face beta terms (CleanFaceWeights, CleanFaceSource)
      const Real d = uh(m,IDN,k,j,i), t = x.w(m,CLW_T,k,j,i);
      Real b2 = 0.0, e, pg, chi_rho, chi_t, cv;
      for (int n = 0; n < 3; ++n) b2 += SQR(x.w(m,CLW_BT+n,k,j,i));
      eos.ThermoAt(d, t, e, pg, chi_rho, chi_t, cv);
      const Real s = 4.0*x.h*x.c*op(m,VET_OP_P,k,j,i)*(1.0 - 0.5*b2)*arad_*t*t*t;
      x.w(m,CLW_GX,k,j,i) = (s > 0.0) ? s/(d*cv + s) : 0.0;
      if (rx) x.w(m,CLW_JV+2,k,j,i) = d;   // CleanStagDrag's rho: a seam ghost's own
    });
    // E_k's ghosts are those of its last change; BT, RE, GX in one exchange before the
    // lateral sums
    CleanGhostDefer(cln_w, CLW_BT, 4);
    CleanGhostFill(cln_w, rx ? CLW_JV+2 : CLW_GX, rx ? 2 : 1);
    if (lat) CleanStagCsLateral(CLW_EK);
  }
  if (beta) return;
  CleanArgs xh = x;   // face momentum state (clean_stag_split)
  xh.hi = cln_split;
  // a pass without the F solve's F: its F^H, each face once in cln_rhsaf into CLF_N
  const bool fh = x.feedback && !cln_fsk;
  xh.fs = true;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const auto *peos = CleanDev(eos, 2);
  CleanArgs xr = x;   // clean_stag_xdonor: the face sources as they are
  xr.xq = -1;
  const CleanArgs *pxr = CleanDevArgs(xr, 3);
  // each face's source once, for its two cells (CLF_SR), and with gflux its a E (CLF_PT)
  const int nf = x.nf;
  const bool af = x.feedback && x.gf && !x.cub, xsp = (x.xq >= 0) && !cln_xdh;
  // one launch per direction, its face code for a compile-time direction (CleanDirCase)
  for (int dir = 0; dir < nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    CleanFor<4>("cln_rhsaf", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
      CleanDirCase(dir, [&](auto d) {
        constexpr int n = decltype(d)::value;
        if (xsp) CleanStagXSplit(x, *pxr, n, m, k, j, i);
        x.fc[n](m,CLF_SR,k,j,i) = CleanFaceSource(x, n, m, k, j, i);
        if (fh) x.fc[n](m,CLF_N,k,j,i) = CleanFaceFlux<n>(xh, CLW_EK, n, m, k, j, i,
                                                          true);
        if (af) x.fc[n](m,CLF_PT,k,j,i) = CleanStagAdvFace(x, xh, n, m, k, j, i);
      });
    });
  }
  const bool xte = x.feedback;   // CLW_XA, CLW_XW
  const bool gpic = cln_gpic;   // clean_gas_newton_switch
  // sdirk2's stage 2: the chord, CLW_RR holds the last pass's T
  const bool gsc = t2_solve == VET_T2S_STAGE2;
  const Real eps_ = std::numeric_limits<Real>::epsilon();
  CleanParForLb<3>("cln_rhs", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
    const auto &eos = *peos;
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hcke = x.h*x.c*op(m,VET_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0, gm = 1.0, sb = 0.0, dop = 0.0, ek = 0.0, gres = 0.0;
    Real eg = 0.0, wq = 0.0;
    bool gr = false;
    if (x.feedback) {
      Real b2 = 0.0, ba = 0.0, by = 0.0, bp = 0.0, y2 = 0.0, fs[3] = {0.0, 0.0, 0.0};
      // clean_stag_fwork: the F solve's F, which the gas takes
      auto ff = [&](const int d, const int kk, const int jj, const int ii) {
        return xh.fc[d](m,CLF_N,kk,jj,ii);
      };
      // gflux: an open boundary face's flux, absorbed in an opaque half cell, enters the
      // cell mean as its half-cell mean, (1 - e^-q)/q, q = chi dx/2 (CleanMomentumChange)
      Real fw[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      for (int n = 0; n < 3; ++n) {   // trip 3: fw stays in registers
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        for (int sd = 0; sd < 2 && x.gf && n < x.nf; ++sd) {
          if (!CleanBoundaryFace(x, n, m, k+sd*dk, j+sd*dj, i+sd*di) ||
              !CleanOpen(x, m, 2*n + sd)) continue;
          const Real q = 0.5*op(m,VET_OP_T,k,j,i)/CleanDivW(x, n, m, k, j, i, sd);
          fw[2*n+sd] = (q > 1.0e-8) ? -expm1(-q)/q : 1.0;
        }
      }
      for (int n = 0; n < 3; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        if (n < x.nf) fs[n] = fw[2*n]*ff(n, k, j, i) + fw[2*n+1]*ff(n, k+dk, j+dj, i+di);
      }
      if (x.cub) {   // the cell's flux vector (CleanStagCellVector)
        Real c[6];
        for (int d = 0; d < 3; ++d) {
          c[2*d] = fw[2*d]*ff(d, k, j, i);
          c[2*d+1] = fw[2*d+1]*ff(d, k+(d == 2), j+(d == 1), i+(d == 0));
        }
        CleanStagCellVector(x, m, k, j, i, c, fs);
        for (int d = 0; d < 3; ++d) fs[d] *= 2.0;
      }
      for (int n = 0; n < 3; ++n) {
        const Real bt = x.w(m,CLW_BT+n,k,j,i);
        if (n < x.nf) by += bt*0.5*fs[n]/x.c;
        y2 += SQR(0.5*fs[n]/x.c);
        // gflux: a E face by face (cubed sphere: the cell's a E, its face families are
        // not the vector's components)
        if (n < x.nf && x.gf && !x.cub) {
          const int dk = (n == 2), dj = (n == 1), di = (n == 0);
          bp += bt*0.5*(x.fc[n](m,CLF_PT,k,j,i) + x.fc[n](m,CLF_PT,k+dk,j+dj,i+di))/x.c;
        }
        b2 += SQR(bt);
        ba += bt*x.w(m,CLW_A1+n,k,j,i)/x.c;
      }
      if (x.gf && x.cub) bp = ba*x.w(m,CLW_EK,k,j,i);
      // y realisable, |y| <= E; the gas eint takes the same sb (CLW_EP): energy exact
      ek = fmax(x.w(m,CLW_EK,k,j,i), 0.0);
      if (y2 > SQR(ek)) by *= ek/std::sqrt(y2);
      gm = 1.0 - 0.5*b2;
      sb = hcke*((0.5*b2 - ba)*x.w(m,CLW_EK,k,j,i) + 2.0*by);
      // gflux: sb = hcke [dop E + 2 beta.F'/c] with the comoving flux F' = F - a E taken
      // face by face, not as the faces' mean F less the cell's a E (at a front these
      // differ by the neighbour's a E, a lagged source the iteration could not settle)
      dop = 0.5*b2 + ba;
      if (x.gf) sb = hcke*(dop*x.w(m,CLW_EK,k,j,i) + 2.0*(by - bp));
      // clean_stag_freal: the comoving part realisable, |E0 - E| <= 2 |beta| E (a face F'
      // of a hot neighbour through a cold opaque cell is no flux of its E); gas: same sb
      if (x.fre) {
        const Real s0 = hcke*dop*x.w(m,CLW_EK,k,j,i), sl = 2.0*hcke*std::sqrt(b2)*ek;
        sb = s0 + fmin(fmax(sb - s0, -sl), sl);
      }
      if (x.gf) x.w(m,CLW_JV+1,k,j,i) = hcke*dop;   // d sb/dE (CleanNewtonUpdate)
      gr = sb != -x.w(m,CLW_EP,k,j,i);   // T took the last pass's sb
    }
    const Real hckp = x.h*x.c*op(m,VET_OP_P,k,j,i)*gm;
    x.w(m,CLW_GM,k,j,i) = gm;
    x.w(m,CLW_EP,k,j,i) = -sb;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, eg, p, chi_rho, chi_t, cv);
      // exchange slopes: the tangent at the iterate's T (Newton, quadratic) where the
      // exchange in B, EO + (1 - r) hc kP a T^4 - r (e(T) - eint), is >= 0; else the
      // secant through T_o (the old vector's gas T) and T, exact at T, so the solution is
      // the same and the exchange in B is EO + (1 - r) hc kP a T_o^4 >= 0 by construction
      // (the tangent overshoots the root of the concave E -> e map when radiation heats
      // cold gas, giving E < 0).
      Real s = 4.0*hckp*arad_*t*t*t, cg = d*cv;
      // clean_stag_eicarry: a carried deficit eint < 0 takes the chord's T_o at 0
      const Real ei = fmax(x.w(m,CLW_EI,k,j,i), 0.0), eo = x.w(m,CLW_EO,k,j,i);
      const bool sec = eg > ei && eo*(cg + s) + cg*hckp*arad_*SQR(SQR(t)) < s*(eg - ei);
      const Real to = sec ? eos.Temperature(d, ei, t) : t;
      if (sec && t - to > 1.0e-6*t) {
        cg = (eg - ei)/(t - to);
        s = hckp*arad_*(t + to)*(t*t + to*to);
      }
      // stage 2: elsewhere the chord through the last pass's T (the tangent at
      // a c_v jump, e.g. H ionisation, cycles across the jump); the root is unchanged
      const Real tp = x.w(m,CLW_RR,k,j,i);
      // The chord only where it exceeds the tangent d cv by more than the round-off
      // of e(T) - e(T_p) (Newton where e(T) is smooth; bitwise tangent for an ideal gas)
      if (gsc && !sec && tp > 0.0 && tp != t) {
        Real ep, pp, crp, ctp, cvp;
        eos.ThermoAt(d, tp, ep, pp, crp, ctp, cvp);
        const Real dx = (eg - ep - cg*(t - tp))*(t - tp);
        if (dx > 16.0*eps_*(fabs(eg) + fabs(ep) + fabs(cg*(t - tp)))*fabs(t - tp)) {
          cg = (eg - ep)/(t - tp);
          s = hckp*arad_*(t + tp)*(t*t + tp*tp);
        }
      }
      rr = (s > 0.0) ? s/(cg + s) : 0.0;   // no exchange where T = 0
      // 1 - r = cg/(cg + s) without cancellation
      x.w(m,CLW_RY,k,j,i) = (s > 0.0) ? cg/(cg + s) : 1.0;
      // clean_gas_newton_switch: a stalled solve's remaining passes drop the gas Newton
      // slope (Picard on the bracketed root T of the gas row)
      if (gpic) {
        rr = 0.0;
        x.w(m,CLW_RY,k,j,i) = 1.0;
      }
      // the gas row's residual at the iterate with these sb, share r: a change of sb
      // between passes enters E with 1 - r (taken whole it contracts only ~0.5 per pass
      // in radiation-dominated opaque gas)
      if (gr) gres = rr*(x.w(m,CLW_EI,k,j,i) + hcke*ek - sb - eg -
                         hckp*arad_*SQR(SQR(t)));
      // gflux: the Doppler part hcke dop E of sb is taken in the Newton step as the
      // exchange (r -> r + dop (1 - r), sb kept at E_k), not lagged: where hcke dop >> 1
      // the lagged term diverges
      if (x.gf) {
        const Real dc = fmin(fmax(dop, 0.0), 1.0);
        rr += dc*(1.0 - rr);
        x.w(m,CLW_RY,k,j,i) *= 1.0 - dc;
      }
      // work of the radiation force on the gas, (v_old + dm/(2 rho)).dm: exactly the
      // kinetic energy the gas gains; with wb_arad also the reference force's work
      // h rho arad_ref.v at the iterate's velocity
      Real pg[3], rf[3] = {0.0, 0.0, 0.0};
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      if (wsp) CleanStagRefForce(x, fref, fr, aref, d, m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (pg[n] + 0.5*dm)*dm/d;
        if (wsp && !wpa) work += x.h*rf[n]*(pg[n] + dm)/d;
      }
      if (wpa) work += fra*wc(m,k,j,i);
      // a drain W > 0 goes on the diagonal as q X, q = W/E_k (same
      // fixed point, b >= 0, X > 0; the write-back books q X); a source stays in B
      if (work > 0.0 && ek > 0.0) {
        wq = work/ek;
        work = 0.0;
      }
    }
    x.w(m,CLW_RR,k,j,i) = gsc ? t : rr;
    // eqr: the exchange via the gas row of the root T, hckp a T^4 - r hcke E =
    // EI - e(T) - sb + (1 - r) hcke E, avoiding the cancellation of two ~ hcke E terms
    // (round-off at h -> infinity exceeds E in opaque cells)
    const Real sdv = CleanSourceDivergence(x, m, k, j, i);
    const Real hes = eso ? x.h*es(m,k,j,i) : 0.0;
    x.w(m,CLW_B,k,j,i) = eqr ?
        x.w(m,CLW_EO,k,j,i) + x.w(m,CLW_EI,k,j,i) - eg - work
        + x.w(m,CLW_RY,k,j,i)*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - sdv + hes :
        x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
        - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work - sdv + sb + gres + hes;
    x.w(m,CLW_PD,k,j,i) = wq;
    // the row's exchange as linearised (B less its rest) and its work
    if (xte) {
      x.w(m,CLW_XA,k,j,i) = eqr ? 0.0 : x.w(m,CLW_B,k,j,i) - (x.w(m,CLW_EO,k,j,i) - work
          - sdv + hes);
      x.w(m,CLW_XW,k,j,i) = work;
    }
    x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
    if (!bud) return;   // clean_conv = budget: the gross terms of the E and gas rows
    const Real xc = hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) + hckp*arad_*SQR(SQR(t)) +
                    fabs(work) + fabs(sb);
    // face flux scale: the old vector's, at least what the cell's E can stream through
    // it in h, c theta E, theta = 1/(1 + h c chi); a physical boundary face takes the
    // cell's streaming alone
    Real fl = eso ? fabs(es(m,k,j,i)) : 0.0;
    const Real eo = fabs(x.w(m,CLW_EO,k,j,i));
    const Real cs = x.c*eo/(1.0 + x.h*x.c*op(m,VET_OP_T,k,j,i));
    for (int n = 0; n < x.nf; ++n) {
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*(n == 2), jf = j + sd*(n == 1), iff = i + sd*(n == 0);
        const auto &f = x.fc[n];
        fl += CleanDivW(x, n, m, k, j, i, sd)*(CleanBoundaryFace(x, n, m, kf, jf, iff) ?
              cs : fmax(fabs(f(m,CLF_G,kf,jf,iff)), x.c*f(m,CLF_TH,kf,jf,iff)*eo));
      }
    }
    x.w(m,CLW_SE,k,j,i) = eo + xc + x.h*fl;
    x.w(m,CLW_SG,k,j,i) = fabs(x.w(m,CLW_EI,k,j,i)) + xc;
  });
  if (!bud) return;
  // an empty cell (no energy, no flux): its rows at 1e-12 of the largest scale
  auto w = cln_w;
  Real smax[2] = {0.0, 0.0};
  // flat range reduction
  const int n1 = ie - is + 1, n2 = je - js + 1, nji = n2*n1, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("cln_smax", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &se, Real &sg) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    se = fmax(se, w(m,CLW_SE,k,j,i));
    sg = fmax(sg, w(m,CLW_SG,k,j,i));
  }, Kokkos::Max<Real>(smax[0]), Kokkos::Max<Real>(smax[1]));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, smax, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  const Real fse = 1.0e-12*smax[0];
  const Real fsg = fmax(1.0e-12*smax[1], 1.0e-300);
  CleanFor("cln_sfloor", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_SE,k,j,i) = fmax(w(m,CLW_SE,k,j,i), fmax(fse, 1.0e-300));
    w(m,CLW_SG,k,j,i) = fmax(w(m,CLW_SG,k,j,i), fsg);
  });
}

}  // namespace radvet
