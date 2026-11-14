//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean.cpp
//! \brief <rad_m1>/scheme = clean: the clean implicit stage solve
//! (docs/dev/rad_clean_scheme_design.md, appendix A).
//!
//! One call solves one backward-Euler stage (time2 STAGE1/STAGE2 with dt_sub = g dt, or a
//! BE step) for E (cells) with the face fluxes eliminated:
//!   face  : F_f = (1-w) theta [Ft_f - h c^2 (d(D E)/dn + cross)] + w c (h+_L E_L + h-_R E_R)
//!   cell  : E (1 + h c kE) + h div F = Et + h c kP a T^4
//! w = smoothstep in C1 = <|I-S|>/J of the formal solution (no optical depth enters),
//! theta = 1/(1 + h c chi_f).  The diagonal-D and half-range parts give an M-matrix.
//! The closure (D, h+-, C1) is built once per step, for t^{n+1}, by the formal solution
//! of the state predicted by this step's physics (no extrapolation).  Gate-1 scope: Cartesian, one MeshBlock, one
//! rank, chat = c.  Gas feedback: Newton on E with the exact local T per cell.

#include <cmath>
#include <cstdio>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_opacity.hpp"

namespace radm1 {

namespace {
// closure slots (cln_c / cln_s)
constexpr int CLN_D11 = 0;   // D11 D22 D33 D12 D13 D23
constexpr int CLN_HP = 6;    // h+_a = HP_a/J, a = 1 2 3
constexpr int CLN_HM = 9;    // h-_a = H_a/J - h+_a (<= 0)
constexpr int CLN_C1 = 12;   // C1 = NL/J
constexpr int CLN_J = 13;    // J of the formal solution (E units)
constexpr int CLN_H = 14;    // H_a of the formal solution (F_a/c, E units), a = 1 2 3
constexpr int CLN_N = 17;
// work slots (cln_w)
constexpr int CLW_EO = 0;    // old vector E
constexpr int CLW_X = 1;     // solution
constexpr int CLW_B = 2;     // right-hand side
constexpr int CLW_AM = 3;    // x1 line: coefficient of E_{i-1}
constexpr int CLW_A0 = 4;    // diagonal (x1 line + x2/x3 diagonal parts)
constexpr int CLW_AP = 5;    // x1 line: coefficient of E_{i+1}
constexpr int CLW_T = 6;     // gas temperature (fixed)
constexpr int CLW_R = 7, CLW_R0 = 8, CLW_P = 9, CLW_V = 10, CLW_S = 11, CLW_TT = 12;
constexpr int CLW_Y = 13, CLW_Z = 14, CLW_CP = 15;
constexpr int CLW_RR = 16;   // gas feedback: r = 4 h c kP a T^3/(rho cv + 4 h c kP a T^3)
constexpr int CLW_EK = 17;   // gas feedback: the Newton iterate E_k
constexpr int CLW_EI = 18;   // gas feedback: eint of the old vector
constexpr int CLW_A1 = 19;   // gas feedback: enthalpy coefficients a_d = v_d + D_db v_b
constexpr int CLW_DM = 22;   // gas feedback: momentum change of the iterate (3)
constexpr int CLW_N = 25;

struct CleanCtx {
  DvceArray5D<Real> cs;          // stage closure
  DvceArray5D<Real> w;           // work
  DvceArray5D<Real> fc[3];       // face coeff. (m,5,..): w, theta, G, q | corr, raw corr
  DvceArray3D<Real> bmf;         // beam face flux at x1min (m, k-ks, j-js)
  Real h, c, idx[3];
  int bc[6];                     // M1_IBC_* x1min x1max x2min x2max x3min x3max
  int is, ie, js, je, ks, ke;
  bool d2, d3;
  int bform;                     // 0 half-range rays, 1 Marshak q
  int bmn;
  Real q, eb[6];
  bool fb;                       // gas feedback: r in the diagonal, enthalpy flux
  bool hmus;                     // ray face flux: limited upwind reconstruction (else mean)
  bool cfct = true;              // correction limited for positivity (else the q gate)
  bool addc = true;              // add the consistency correction (G part)
  bool hrr = true;               // half-range factors from face ray fluxes (fc 3, 4)
};

// the cell's exact backward-Euler temperature for a given E (constant opacities):
//   d T/gm1 + hckpa T^4 = eold + hckeE,  monotone in T, one root in [0, gm1 rhs/d]
KOKKOS_INLINE_FUNCTION Real CleanTLocal(const Real d, const Real gm1, const Real eold,
                                        const Real hckpa, const Real hckee, const Real tg) {
  const Real rhs = eold + hckee;
  if (!(rhs > 0.0)) return 0.0;
  Real lo = 0.0, hi = gm1*rhs/d;
  Real t = fmin(fmax(tg, 0.0), hi);
  for (int n = 0; n < 200; ++n) {
    const Real t3 = t*t*t;
    const Real g = d*t/gm1 + hckpa*t3*t - rhs;
    if (g > 0.0) {hi = t;} else {lo = t;}
    const Real dg = d/gm1 + 4.0*hckpa*t3;
    Real tn = t - g/dg;
    if (!(tn > lo && tn < hi)) tn = 0.5*(lo + hi);
    if (fabs(tn - t) <= 1.0e-15*fmax(t, 1.0e-300)) {t = tn; break;}
    t = tn;
  }
  return t;
}

KOKKOS_INLINE_FUNCTION Real CleanSmooth(Real x) {
  x = fmin(fmax(x, 0.0), 1.0);
  return x*x*(3.0 - 2.0*x);
}

// the ray-moment consistency correction of the half-range part of face (dir, k, j, i):
// w c (H_face - (h+_L J_L + h-_R J_R)), H_face the rays' flux at the face (limited upwind
// reconstruction of the half-range fluxes, or the face mean), clamped to [-J_R, J_L];
// with it the upwind part acts on E - J only (E = J in a steady streaming field)
KOKKOS_INLINE_FUNCTION Real CleanCorr(const CleanCtx &x, const int dir, const int m,
                                      const int k, const int j, const int i) {
  const int kl = (dir == 2) ? k-1 : k;
  const int jl = (dir == 1) ? j-1 : j;
  const int il = (dir == 0) ? i-1 : i;
  const Real wf = x.fc[dir](m,0,k,j,i);
  if (!(wf > 0.0)) return 0.0;
    const Real jlv = x.cs(m,CLN_J,kl,jl,il), jrv = x.cs(m,CLN_J,k,j,i);
    const Real hup = x.cs(m,CLN_HP+dir,kl,jl,il)*jlv + x.cs(m,CLN_HM+dir,k,j,i)*jrv;
    Real hf;
    if (x.hmus) {
      // upwind-biased limited reconstruction of the rays' half-range fluxes to the
      // face: H+ from the left cell (HP = h+ J), H- from the right (HM = h- J);
      // minmod slopes: 2nd order where smooth, the donor value at a jump
      const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
      auto hpj = [&](const int kk, const int jj, const int ii) {
        return x.cs(m,CLN_HP+dir,kk,jj,ii)*x.cs(m,CLN_J,kk,jj,ii);};
      auto hmj = [&](const int kk, const int jj, const int ii) {
        return x.cs(m,CLN_HM+dir,kk,jj,ii)*x.cs(m,CLN_J,kk,jj,ii);};
      auto mm = [](const Real a, const Real b) {
        return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);};
      const Real pl = hpj(kl,jl,il);
      const Real sp = mm(pl - hpj(kl-dk,jl-dj,il-di), hpj(k,j,i) - pl);
      const Real mr = hmj(k,j,i);
      const Real sm = mm(mr - hmj(kl,jl,il), hmj(k+dk,j+dj,i+di) - mr);
      hf = fmin(fmax((pl + 0.5*sp) + (mr - 0.5*sm), -jrv), jlv);
    } else {
      hf = fmin(fmax(0.5*(x.cs(m,CLN_H+dir,kl,jl,il) + x.cs(m,CLN_H+dir,k,j,i)),
                     -jrv), jlv);
    }
  return wf*x.c*(hf - hup);
}

// face flux of direction dir (0,1,2) at face (k,j,i) (face between cell-1 and cell along
// dir) for the E vector in work slot e.  withg: include the E-independent part G.
KOKKOS_INLINE_FUNCTION Real CleanFlux(const CleanCtx &x, const int e, const int dir,
                                      const int m, const int k, const int j, const int i,
                                      const bool withg, const bool adv = true) {
  const int kl = (dir == 2) ? k-1 : k;
  const int jl = (dir == 1) ? j-1 : j;
  const int il = (dir == 0) ? i-1 : i;
  const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
  const int hi = (dir == 0) ? x.ie+1 : ((dir == 1) ? x.je+1 : x.ke+1);
  const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
  const Real c = x.c;
  // physical boundary faces (non-periodic)
  if (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) {
    const int b = x.bc[2*dir];
    if (b == M1_IBC_REFLECT) return 0.0;
    const Real er = x.w(m,e,k,j,i);
    Real fin = 0.0;
    bool beam = false;
    if (dir == 0 && x.bmn > 0) {
      const Real fb = x.bmf(m,k-x.ks,j-x.js);
      if (fb > 0.0) {fin = fb; beam = true;}
    }
    if (x.bform == 1) {
      if (beam) return withg ? fin : 0.0;
      return c*x.q*((withg ? x.eb[2*dir] : 0.0) - er);
    }
    if (!beam) fin = 0.25*c*x.eb[2*dir];
    return (withg ? fin : 0.0) + c*x.cs(m,CLN_HM+dir,k,j,i)*er;
  }
  if (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC) {
    const int b = x.bc[2*dir+1];
    if (b == M1_IBC_REFLECT) return 0.0;
    const Real el = x.w(m,e,kl,jl,il);
    if (x.bform == 1) return c*x.q*(el - (withg ? x.eb[2*dir+1] : 0.0));
    return c*x.cs(m,CLN_HP+dir,kl,jl,il)*el - (withg ? 0.25*c*x.eb[2*dir+1] : 0.0);
  }
  const Real wf = x.fc[dir](m,0,k,j,i);
  const Real th = x.fc[dir](m,1,k,j,i);
  const Real el = x.w(m,e,kl,jl,il), er = x.w(m,e,k,j,i);
  Real f = 0.0;
  if (wf < 1.0) {
    const int dd = CLN_D11 + dir;   // D11 / D22 / D33
    Real grad = (x.cs(m,dd,k,j,i)*er - x.cs(m,dd,kl,jl,il)*el)*x.idx[dir];
    // cross derivatives d_b(D_{dir,b} E) at the face, b != dir (4-cell average)
    for (int b = 0; b < 3; ++b) {
      if (b == dir) continue;
      if (b == 1 && !x.d2) continue;
      if (b == 2 && !x.d3) continue;
      const int ab = (dir + b == 1) ? 3 : ((dir + b == 2) ? 4 : 5);   // D12 D13 D23
      const int dk = (b == 2), dj = (b == 1), di = (b == 0);
      Real sp = x.cs(m,ab,k+dk,j+dj,i+di)*x.w(m,e,k+dk,j+dj,i+di) +
                x.cs(m,ab,kl+dk,jl+dj,il+di)*x.w(m,e,kl+dk,jl+dj,il+di);
      Real sm = x.cs(m,ab,k-dk,j-dj,i-di)*x.w(m,e,k-dk,j-dj,i-di) +
                x.cs(m,ab,kl-dk,jl-dj,il-di)*x.w(m,e,kl-dk,jl-dj,il-di);
      grad += 0.25*(sp - sm)*x.idx[b];
    }
    f = (1.0 - wf)*th*(-x.h*c*c*grad + (withg ? x.fc[dir](m,2,k,j,i) : 0.0));
    // lab-frame flux: the source -c chi (F - A), A = v E + v.P at the face (central),
    // gives F = theta (Ft - h c^2 grad) + (1 - theta) A
    if (x.fb) {
      f += (1.0 - wf)*(1.0 - th)*0.5*(x.w(m,CLW_A1+dir,kl,jl,il)*el +
                                      x.w(m,CLW_A1+dir,k,j,i)*er);
    }
  }
  if (wf > 0.0 && x.hrr) {
    // consistent half-range face: c (h+_f E_L + h-_f E_R) with h+_f = H+_f/J_L,
    // h-_f = H-_f/J_R from the rays' half-range fluxes at the face (= the rays' face flux
    // when E = J); h+_f >= 0 >= h-_f keeps the M-matrix (positive by construction)
    f += wf*c*(x.fc[dir](m,3,k,j,i)*el + x.fc[dir](m,4,k,j,i)*er);
  } else if (wf > 0.0) {
    f += wf*c*(x.cs(m,CLN_HP+dir,kl,jl,il)*el + x.cs(m,CLN_HM+dir,k,j,i)*er);
    // ray-moment consistency (high-order/low-order), see CleanCorr: gate mode q*corr,
    // fct mode the limited correction stored per face
    if (withg && x.addc) {
      f += x.cfct ? x.fc[dir](m,3,k,j,i) : x.fc[dir](m,3,k,j,i)*CleanCorr(x, dir, m, k, j, i);
    }
  }
  (void)adv;
  return f;
}


// gas feedback: the momentum change dm of the iterate in work slot e (face-flux
// relaxation h c chi (F - (1 - w) A_f) of momentum h chi (F - (1-w) A_f)/c, half of each
// face to each cell) into CLW_DM, and the enthalpy coefficients a_d = v_d + D_db v_b of
// v = (m_old + dm)/rho into CLW_A1 (implicit velocity).  Returns max |dv| change.
KOKKOS_INLINE_FUNCTION Real CleanCellMom(const CleanCtx &x, const int e, const int m,
                                         const int k, const int j, const int i,
                                         const int dir) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  Real s = 0.0;
  for (int side = 0; side < 2; ++side) {
    const int kf = k + side*dk, jf = j + side*dj, if_ = i + side*di;
    const Real th = x.fc[dir](m,1,kf,jf,if_), wf = x.fc[dir](m,0,kf,jf,if_);
    const Real f = CleanFlux(x, e, dir, m, kf, jf, if_, true);
    const Real af = 0.5*(x.w(m,CLW_A1+dir,kf-dk,jf-dj,if_-di)*x.w(m,e,kf-dk,jf-dj,if_-di) +
                         x.w(m,CLW_A1+dir,kf,jf,if_)*x.w(m,e,kf,jf,if_));
    s += 0.5*(1.0/th - 1.0)*(f - (1.0 - wf)*af)/(x.c*x.c);
  }
  return s;
}

// ghost fill of components [c0, c0+nc) of a 5D view: periodic wrap or edge copy
void CleanGhost(DvceArray5D<Real> a, const int c0, const int nc, const int nmb1,
                const int *bc, const int is, const int ie, const int js, const int je,
                const int ks, const int ke, const int ng, const bool d2, const bool d3) {
  const int nx1 = ie - is + 1, nx2 = je - js + 1, nx3 = ke - ks + 1;
  const bool p1 = (bc[0] == M1_IBC_PERIODIC), p2 = (bc[2] == M1_IBC_PERIODIC);
  const bool p3 = (bc[4] == M1_IBC_PERIODIC);
  par_for("cln_gh1", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
    const int c = c0 + n;
    a(m,c,k,j,is-1-g) = p1 ? a(m,c,k,j,ie-g) : a(m,c,k,j,is);
    a(m,c,k,j,ie+1+g) = p1 ? a(m,c,k,j,is+g) : a(m,c,k,j,ie);
  });
  if (d2) {
    par_for("cln_gh2", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, 0, ng-1, is-ng, ie+ng,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int g, const int i) {
      const int c = c0 + n;
      a(m,c,k,js-1-g,i) = p2 ? a(m,c,k,je-g,i) : a(m,c,k,js,i);
      a(m,c,k,je+1+g,i) = p2 ? a(m,c,k,js+g,i) : a(m,c,k,je,i);
    });
  }
  if (d3) {
    const int jlo = d2 ? js-ng : js, jhi = d2 ? je+ng : je;
    par_for("cln_gh3", DevExeSpace(), 0, nmb1, 0, nc-1, 0, ng-1, jlo, jhi, is-ng, ie+ng,
    KOKKOS_LAMBDA(const int m, const int n, const int g, const int j, const int i) {
      const int c = c0 + n;
      a(m,c,ks-1-g,j,i) = p3 ? a(m,c,ke-g,j,i) : a(m,c,ks,j,i);
      a(m,c,ke+1+g,j,i) = p3 ? a(m,c,ks+g,j,i) : a(m,c,ke,j,i);
    });
  }
  (void)nx1; (void)nx2; (void)nx3;
}

// sum over active cells of a(c1)*a(c2)
Real CleanDot(DvceArray5D<Real> a, const int c1, const int c2, const int nmb,
              const int is, const int ie, const int js, const int je,
              const int ks, const int ke) {
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nkji = nk*nj*ni, nji = nj*ni;
  Real s = 0.0;
  Kokkos::parallel_reduce("cln_dot", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &acc) {
    int m = idx/nkji, r = idx - m*nkji;
    int k = r/nji; r -= k*nji;
    int j = r/ni, i = r - j*ni;
    acc += a(m,c1,k+ks,j+js,i+is)*a(m,c2,k+ks,j+js,i+is);
  }, s);
  return s;
}
}  // namespace

//----------------------------------------------------------------------------------------
void RadiationM1::CleanAlloc() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;
  const int nmb = pmy_pack->nmb_thispack;
  Kokkos::realloc(cln_c, nmb, CLN_N, n3, n2, n1);
  Kokkos::realloc(cln_s, nmb, CLN_N, n3, n2, n1);
  Kokkos::realloc(cln_v, nmb, 2, n3, n2, n1);
  Kokkos::realloc(cln_w, nmb, CLW_N, n3, n2, n1);
  Kokkos::realloc(cln_fc1, nmb, 5, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, 5, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, 5, n3, n2, n1);
  Kokkos::deep_copy(cln_w, 0.0);
  Kokkos::deep_copy(cln_fc1, 0.0);
  Kokkos::deep_copy(cln_fc2, 0.0);
  Kokkos::deep_copy(cln_fc3, 0.0);
  int nr = 1;
#if MPI_PARALLEL_ENABLED
  nr = global_variable::nranks;
#endif
  const bool cart = !pmy_pack->pmesh->use_cubed_sphere &&
                    !pmy_pack->pmesh->use_spherical_polar;
  if (nmb != 1 || nr != 1 || !cart || !vet_sc || !vet_full ||
      chat != c_light) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "<rad_m1>/scheme = clean (gate 1) needs: one MeshBlock, one rank, a Cartesian "
      << "mesh, closure = vet_sc with vet_tensor = full, chat = c"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the formal solution must accumulate the half-range sums and C1 for this scheme
  impl_beam_hr = true;
  thsw_c1 = true;
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! closure: formal solution of the state predicted at t^{n+1}, once per step (STAGE1, or
//! a BE solve with a new dt), into cln_c and the stage closure cln_s.  STAGE2 and a BE
//! redo with the same dt keep it.
void RadiationM1::CleanClosure(int t2s) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, ng = indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool d2 = pmy_pack->pmesh->multi_d, d3 = pmy_pack->pmesh->three_d;
  if (t2s == M1_T2S_STAGE2) return;
  const Real tnow = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  // One closure per step, for t^{n+1} (both hesdirk2 stage solves, and a backward-Euler
  // redo, end at t^{n+1}).  A redo with the same dt keeps it; a new step or a new dt
  // rebuilds it.
  if (cln_havec && tnow == cln_tc && dt == cln_dtb) return;
  // The closure belongs to a radiation field, so it is never extrapolated: the formal
  // solution runs on the state PREDICTED at t^{n+1} by the physics of this step, no
  // differences of past states:
  //  - gas rho, momentum, energy: the hydro's explicit (Heun stage-1) state, which is
  //    already at t^{n+1}, plus dt K1 of the radiation coupling (K1 = f(U^n), the FSAL
  //    rate of the implicit part: exchange, radiation force, their work);
  //  - T from that gas; the exchange can only relax T toward the radiation temperature
  //    (E/a)^{1/4}, so the predicted T never crosses it (the physical bound, active
  //    only where dt K1 would carry T past equilibrium);
  //  - E: dt K1_E; a decrease is applied as exp(dt K1_E/E^n) (loss at the rate of U^n:
  //    exact for constant-rate absorption/escape, E > 0); an increase linearly;
  //  - opacities at the predicted (rho, T).
  // Second order: the predicted state is within O(dt^2) of U^{n+1}.  Without a valid
  // K1 (first step, after a restart without slope) the rates are zero.  The flux enters
  // only the diffusion boundary rays: F^n.
  const bool k1 = cln_pred && t2_ok && (t2k1.extent_int(0) > nmb1);
  auto &indcs_ = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs_.nx1 + 2*ng;
  const int n2 = (indcs_.nx2 > 1) ? indcs_.nx2 + 2*ng : 1;
  const int n3 = (indcs_.nx3 > 1) ? indcs_.nx3 + 2*ng : 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool hh = fl.on;
  const bool need_op = hh && !opac_zero && (opacity_type != M1_OPAC_CONST);
  if (need_op && (cln_opb.extent_int(0) != opac.extent_int(0) ||
                  cln_opb.extent_int(1) != opac.extent_int(1))) {
    Kokkos::realloc(cln_opb, opac.extent_int(0), opac.extent_int(1), opac.extent_int(2),
                    opac.extent_int(3), opac.extent_int(4));
  }
  if (need_op) Kokkos::deep_copy(cln_opb, opac);
  {
    auto iw_ = iw;
    auto u0_ = u0;
    auto uh = fl.u0;
    auto eos = fl.eos;
    auto kr = t2k1;
    auto opac_ = opac;
    const int otype = opacity_type;
    const Real kp = kappa_p, ke = kappa_e, kf = kappa_f, ks_ = kappa_s;
    const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
    M1OpacTab ot = otab;
    const bool osc = vscat, nop = need_op;
    const Real ar = arad;
    const int is_ = is, ie_ = ie, js_ = js, je_ = je, ks2 = ks, ke2 = ke;
    par_for("cln_vin", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // K1 lives on active cells only; ghosts take the explicit state alone
      const bool act = (k >= ks2 && k <= ke2 && j >= js_ && j <= je_ &&
                        i >= is_ && i <= ie_);
      const bool rk = k1 && act;
      const Real e0 = u0_(m,M1_E,k,j,i);
      Real e = e0;
      if (rk) {
        const Real de = dt*kr(m,M1_T2_E,k,j,i);
        e = (de >= 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : e0);
      }
      Real d = 1.0, t = 0.0;
      if (hh) {
        d = uh(m,IDN,k,j,i);
        const Real id = 1.0/fmax(d, 1.0e-300);
        Real mx[3], eg = uh(m,IEN,k,j,i);
        for (int n = 0; n < 3; ++n) mx[n] = uh(m,IM1+n,k,j,i);
        const Real eint_h = fmax(eg - 0.5*(mx[0]*mx[0] + mx[1]*mx[1] + mx[2]*mx[2])*id,
                                 0.0);
        const Real th = eos.Temperature(d, eint_h);
        t = th;
        if (rk) {
          for (int n = 0; n < 3; ++n) mx[n] += dt*kr(m,M1_T2_M1+n,k,j,i);
          eg += dt*kr(m,M1_T2_EN,k,j,i);
          const Real eint = eg - 0.5*(mx[0]*mx[0] + mx[1]*mx[1] + mx[2]*mx[2])*id;
          const Real tp = eos.Temperature(d, fmax(eint, 0.0));
          const Real tr = sqrt(sqrt(fmax(e, 0.0)/ar));
          t = tp;
          if (th <= tr) {
            t = fmin(tp, tr);
          } else {
            t = fmax(tp, tr);
          }
          if (!(eint > 0.0)) t = th;   // the rate would empty the gas: keep the explicit T
        }
      }
      if (act) {
        iw_(m,M1_IW_EN,k,j,i) = e;
        iw_(m,M1_IW_TP,k,j,i) = t;
        iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
        iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F1+1,k,j,i);
        iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F1+2,k,j,i);
      }
      if (nop) {
        Real op, oe, of, os;
        if (otype == M1_OPAC_TABLE) {
          M1TableOpacities(ot, d, t, op, oe, of, os);
          M1ScatterEos(ot, eos, d, t, of, os);
        } else {
          M1Opacities(otype, d, t, kp, ke, kf, ks_, rref, tref, aa, bb, op, oe, of, os);
        }
        opac_(m,M1_OP_P,k,j,i) = d*op;
        opac_(m,M1_OP_E,k,j,i) = d*oe;
        opac_(m,M1_OP_T,k,j,i) = d*(of + os);
        if (osc) {opac_(m,M1_OP_S,k,j,i) = d*os;}
      }
    });
  }
  // the formal solution of the (predicted) state itself: no dt-dependent relaxed source
  // or relaxed scattering J
  vcol_srelax = false;
  vcol_sbl = false;
  vcol_sesc = false;
  vscat_jrel = false;
  VetShortChar();
  if (need_op) Kokkos::deep_copy(opac, cln_opb);   // the solve uses the stage opacity
  cln_tc = tnow;
  cln_dtb = dt;
  cln_havec = true;
  {
    // closure of each cell from the rays where they carry the cell's radiation
    // (J >= cln_raymin E); else (the moment field holds radiation the rays do not, e.g.
    // next to a beam the first-order moment transport has widened) from the moment
    // state's own flux: the beam + isotropic closure with al = |F|/(c E), C1 = 1.  No optical depth enters.
    auto vc = vet_cell;
    auto cc = cln_c;
    auto uu = u0;
    auto iwc = iw;
    auto fv = cln_v;
    const Real rmin = cln_raymin, cl = c_light;
    par_for("cln_vcopy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jr = vc(m,M1_VET_J,k,j,i);
      const Real em = iwc(m,M1_IW_EN,k,j,i);   // E of the predicted state
      if (jr > 0.0 && jr >= rmin*em) {
        for (int a = 0; a < 6; ++a) cc(m,CLN_D11+a,k,j,i) = vc(m,M1_VET_D11+a,k,j,i);
        for (int a = 0; a < 3; ++a) {
          const Real hp = vc(m,M1_VET_HP1+a,k,j,i)/jr;
          const Real hh = vc(m,M1_VET_H1+a,k,j,i)/jr;
          cc(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
          cc(m,CLN_HM+a,k,j,i) = fmin(hh - hp, 0.0);
        }
        cc(m,CLN_C1,k,j,i) = fmax(vc(m,M1_VET_NL,k,j,i)/jr, 0.0);
        cc(m,CLN_J,k,j,i) = jr;
        for (int a = 0; a < 3; ++a) cc(m,CLN_H+a,k,j,i) = vc(m,M1_VET_H1+a,k,j,i);
        fv(m,0,k,j,i) = 1.0;
      } else {
        fv(m,0,k,j,i) = 0.0;
        cc(m,CLN_J,k,j,i) = 0.0;
        for (int a = 0; a < 3; ++a) cc(m,CLN_H+a,k,j,i) = 0.0;
        // beam + isotropic closure of the moment state: a fraction al = |f| streams
        // along n = f/|f|, the rest is isotropic (exact for both limits; a beam never
        // flows backward)
        Real f[3], fm = 0.0;
        for (int a = 0; a < 3; ++a) {
          f[a] = (em > 0.0) ? uu(m,M1_F1+a,k,j,i)/(cl*em) : 0.0;
          fm += f[a]*f[a];
        }
        fm = std::sqrt(fm);
        // isotropic (al = 0): a closure built from the cell's own flux feeds back on
        // itself (a sideways stream along a boundary row sustained itself, cbg 11-14)
        const Real al = 0.0*fmin(fm, 1.0);
        Real n[3];
        for (int a = 0; a < 3; ++a) n[a] = (fm > 0.0) ? f[a]/fm : 0.0;
        for (int a = 0; a < 3; ++a) {
          cc(m,CLN_HP+a,k,j,i) = al*fmax(n[a], 0.0) + 0.25*(1.0 - al);
          cc(m,CLN_HM+a,k,j,i) = al*fmin(n[a], 0.0) - 0.25*(1.0 - al);
          cc(m,CLN_D11+a,k,j,i) = al*n[a]*n[a] + (1.0 - al)/3.0;
        }
        cc(m,CLN_D11+3,k,j,i) = al*n[0]*n[1];
        cc(m,CLN_D11+4,k,j,i) = al*n[0]*n[2];
        cc(m,CLN_D11+5,k,j,i) = al*n[1]*n[2];
        cc(m,CLN_C1,k,j,i) = 1.0;
      }
    });
  }
  // cells the rays leave empty take the closure (D, h+-) of their ray-described
  // neighbours, spreading one cell per pass (cln_fillpass passes); the moment transport
  // widens a beam by a few cells, and its radiation there keeps the rays' direction.
  // Cells still empty keep the isotropic closure.  J = H = 0 there (no consistency term).
  {
    int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
    auto cc = cln_c;
    auto fv = cln_v;
    for (int pass = 0; pass < cln_fillpass; ++pass) {
      CleanGhost(cln_c, 0, CLN_N, nmb1, bc, is, ie, js, je, ks, ke, ng, d2, d3);
      CleanGhost(cln_v, 0, 2, nmb1, bc, is, ie, js, je, ks, ke, ng, d2, d3);
      par_for("cln_fill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        fv(m,1,k,j,i) = fv(m,0,k,j,i);
        if (fv(m,0,k,j,i) > 0.0) return;
        Real acc[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        int nv = 0;
        const int dk = d3 ? 1 : 0, dj = d2 ? 1 : 0;
        for (int kk = -dk; kk <= dk; ++kk) {
          for (int jj = -dj; jj <= dj; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
              if (fv(m,0,k+kk,j+jj,i+ii) <= 0.0) continue;
              ++nv;
              for (int a = 0; a < 6; ++a) acc[a] += cc(m,CLN_D11+a,k+kk,j+jj,i+ii);
              for (int a = 0; a < 3; ++a) acc[6+a] += cc(m,CLN_HP+a,k+kk,j+jj,i+ii);
            }
          }
        }
        if (nv == 0) return;
        fv(m,1,k,j,i) = 1.0;
        // only cells empty before this pass are written, only valid ones are read
        for (int a = 0; a < 6; ++a) cc(m,CLN_D11+a,k,j,i) = acc[a]/nv;
        for (int a = 0; a < 3; ++a) cc(m,CLN_HP+a,k,j,i) = acc[6+a]/nv;
      });
      // h- of the filled cells from the neighbours' h- (second sweep, same neighbours)
      par_for("cln_fill2", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (fv(m,0,k,j,i) > 0.0 || fv(m,1,k,j,i) <= 0.0) return;
        Real acc[3] = {0, 0, 0};
        int nv = 0;
        const int dk = d3 ? 1 : 0, dj = d2 ? 1 : 0;
        for (int kk = -dk; kk <= dk; ++kk) {
          for (int jj = -dj; jj <= dj; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
              if (fv(m,0,k+kk,j+jj,i+ii) <= 0.0) continue;
              ++nv;
              for (int a = 0; a < 3; ++a) acc[a] += cc(m,CLN_HM+a,k+kk,j+jj,i+ii);
            }
          }
        }
        for (int a = 0; a < 3; ++a) cc(m,CLN_HM+a,k,j,i) = acc[a]/fmax(nv, 1);
      });
      par_for("cln_fill3", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        fv(m,0,k,j,i) = fv(m,1,k,j,i);
      });
    }
  }
  {
    // the stage closure is the closure of the predicted state (no extrapolation)
    auto cc = cln_c;
    auto cs = cln_s;
    const bool dbge = cln_dbg_edd;
    par_for("cln_cset", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = 0; n < CLN_N; ++n) cs(m,n,k,j,i) = cc(m,n,k,j,i);
      if (dbge) {   // diagnostic only: frozen isotropic closure
        for (int a = 0; a < 3; ++a) {
          cs(m,CLN_D11+a,k,j,i) = 1.0/3.0; cs(m,CLN_D11+3+a,k,j,i) = 0.0;
          cs(m,CLN_HP+a,k,j,i) = 0.25; cs(m,CLN_HM+a,k,j,i) = -0.25;
        }
        cs(m,CLN_C1,k,j,i) = 0.0;
      }
    });
  }
  int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  CleanGhost(cln_s, 0, CLN_N, nmb1, bc, is, ie, js, je, ks, ke, ng, d2, d3);
}

//----------------------------------------------------------------------------------------
//! face coefficients w, theta, G = Ft (the old-vector face flux) for this stage, and the
//! x1-line + diagonal coefficients of the operator for the preconditioner
void RadiationM1::CleanFaceCoef() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool d2 = pmy_pack->pmesh->multi_d, d3 = pmy_pack->pmesh->three_d;
  const Real h = dt_sub, c = c_light;
  const Real c0 = cln_c0, c1 = cln_c1;
  const bool t2 = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;
  auto cw = cln_w;
  const bool gate = cln_gate;
  const bool hrr = cln_hrr, mus = cln_hmus;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !d2) continue;
    if (dir == 2 && !d3) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int ihi = (dir == 0) ? ie+1 : ie, jhi = (dir == 1) ? je+1 : je;
    const int khi = (dir == 2) ? ke+1 : ke;
    const int tf = M1_T2_F1 + dir;
    par_for("cln_fcoef", DevExeSpace(), 0, nmb1, ks, khi, js, jhi, is, ihi,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = (dir == 2) ? k-1 : k, jl = (dir == 1) ? j-1 : j;
      const int il = (dir == 0) ? i-1 : i;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real cmax = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      fc(m,0,k,j,i) = CleanSmooth((cmax - c0)/(c1 - c0));
      fc(m,1,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,2,k,j,i) = f0(m,k,j,i) + (t2 ? inc(m,tf,k,j,i) : 0.0);
      // ray-moment consistency gate q = r^2, r = the worse of min(E,J)/max(E,J) of the
      // two cells (E of the old vector, J of the stage closure): the quasi-static ray
      // flux describes the moment field only where the two agree
      auto rat = [&](const Real e, const Real jj) {
        const Real a = fmax(e, 0.0), b = fmax(jj, 0.0);
        return (fmax(a, b) > 0.0) ? fmin(a, b)/fmax(a, b) : 1.0;};
      const Real r = fmin(rat(cw(m,CLW_EO,k,j,i), cs(m,CLN_J,k,j,i)),
                          rat(cw(m,CLW_EO,kl,jl,il), cs(m,CLN_J,kl,jl,il)));
      fc(m,3,k,j,i) = gate ? r*r : 1.0;
      if (hrr) {
        // the rays' half-range fluxes HP = h+ J (rightward) and HM = h- J at the face:
        // limited upwind reconstruction (HP from the left, HM from the right) or mean
        const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
        auto hpj = [&](const int kk, const int jj, const int ii) {
          return cs(m,CLN_HP+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);};
        auto hmj = [&](const int kk, const int jj, const int ii) {
          return cs(m,CLN_HM+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);};
        auto mm = [](const Real a, const Real b) {
          return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);};
        Real hpf, hmf;
        if (mus) {
          const Real pl = hpj(kl,jl,il);
          hpf = pl + 0.5*mm(pl - hpj(kl-dk,jl-dj,il-di), hpj(k,j,i) - pl);
          const Real mr = hmj(k,j,i);
          hmf = mr - 0.5*mm(mr - hmj(kl,jl,il), hmj(k+dk,j+dj,i+di) - mr);
        } else {
          hpf = 0.5*(hpj(kl,jl,il) + hpj(k,j,i));
          hmf = 0.5*(hmj(kl,jl,il) + hmj(k,j,i));
        }
        const Real jlv = cs(m,CLN_J,kl,jl,il), jrv = cs(m,CLN_J,k,j,i);
        // donor ratios; a donor without ray radiation keeps its own (filled) h
        fc(m,3,k,j,i) = (jlv > 0.0) ? fmin(fmax(hpf/jlv, 0.0), 1.0) : cs(m,CLN_HP+dir,kl,jl,il);
        fc(m,4,k,j,i) = (jrv > 0.0) ? fmax(fmin(hmf/jrv, 0.0), -1.0) : cs(m,CLN_HM+dir,k,j,i);
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! dst = A src (linear part, no G), over the active cells
void RadiationM1::CleanApply(int src, int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanFillGhost(src);
  CleanCtx x;
  x.cs = cln_s; x.w = cln_w; x.fc[0] = cln_fc1; x.fc[1] = cln_fc2; x.fc[2] = cln_fc3;
  x.bmf = ibm_f; x.h = dt_sub; x.c = c_light;
  auto &size = pmy_pack->pmb->mb_size;
  x.idx[0] = 1.0/size.h_view(0).dx1; x.idx[1] = 1.0/size.h_view(0).dx2;
  x.idx[2] = 1.0/size.h_view(0).dx3;
  x.bc[0] = ibc_x1min; x.bc[1] = ibc_x1max; x.bc[2] = ibc_x2min; x.bc[3] = ibc_x2max;
  x.bc[4] = ibc_x3min; x.bc[5] = ibc_x3max;
  x.is = is; x.ie = ie; x.js = js; x.je = je; x.ks = ks; x.ke = ke;
  x.d2 = pmy_pack->pmesh->multi_d; x.d3 = pmy_pack->pmesh->three_d;
  x.bform = cln_bdry_form; x.bmn = ibm_n; x.q = marshak_q;
  for (int n = 0; n < 6; ++n) x.eb[n] = 0.0;
  x.eb[0] = iebath_x1min; x.eb[1] = iebath_x1max;
  x.fb = gas_feedback;
  x.hmus = cln_hmus;
  x.cfct = !cln_gate;
  x.hrr = cln_hrr;
  x.addc = !cln_hrr;
  auto op = opac;
  par_for("cln_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real rr = x.fb ? x.w(m,CLW_RR,k,j,i) : 0.0;
    Real a = x.w(m,src,k,j,i)*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr));
    a += x.h*x.idx[0]*(CleanFlux(x,src,0,m,k,j,i+1,false) - CleanFlux(x,src,0,m,k,j,i,false));
    if (x.d2) a += x.h*x.idx[1]*(CleanFlux(x,src,1,m,k,j+1,i,false) -
                                 CleanFlux(x,src,1,m,k,j,i,false));
    if (x.d3) a += x.h*x.idx[2]*(CleanFlux(x,src,2,m,k+1,j,i,false) -
                                 CleanFlux(x,src,2,m,k,j,i,false));
    x.w(m,dst,k,j,i) = a;
  });
}

//----------------------------------------------------------------------------------------
void RadiationM1::CleanFillGhost(int comp) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  CleanGhost(cln_w, comp, 1, pmy_pack->nmb_thispack - 1, bc, indcs.is, indcs.ie,
             indcs.js, indcs.je, indcs.ks, indcs.ke, indcs.ng, pmy_pack->pmesh->multi_d,
             pmy_pack->pmesh->three_d);
}

//----------------------------------------------------------------------------------------
//! dst = M^{-1} src: x1-line tridiagonal solve (Thomas) with the x2/x3 diagonal parts
void RadiationM1::CleanPrecond(int src, int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  par_for("cln_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = (i > is) ? w(m,CLW_AM,k,j,i) : 0.0;
      const Real ap = (i < ie) ? w(m,CLW_AP,k,j,i) : 0.0;
      const Real den = w(m,CLW_A0,k,j,i) - am*cprev;
      const Real cp = ap/den;
      const Real dp = (w(m,src,k,j,i) - am*dprev)/den;
      w(m,CLW_CP,k,j,i) = cp;
      w(m,dst,k,j,i) = dp;
      cprev = cp; dprev = dp;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,CLW_CP,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! right-preconditioned BiCGStab on A X = B; returns iterations (< 0: not converged)
int RadiationM1::CleanBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nmb1 = nmb - 1;
  auto w = cln_w;
  auto dot = [&](int a, int b) {return CleanDot(w, a, b, nmb, is, ie, js, je, ks, ke);};
  auto axpy = [&](int dst, int a, Real sa, int b, Real sb) {   // dst = sa*a + sb*b
    par_for("cln_axpy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,dst,k,j,i) = sa*w(m,a,k,j,i) + sb*w(m,b,k,j,i);
    });
  };
  CleanApply(CLW_X, CLW_R);
  axpy(CLW_R, CLW_B, 1.0, CLW_R, -1.0);
  axpy(CLW_R0, CLW_R, 1.0, CLW_R, 0.0);
  Kokkos::deep_copy(Kokkos::subview(w, Kokkos::ALL, CLW_P, Kokkos::ALL, Kokkos::ALL,
                    Kokkos::ALL), 0.0);
  Kokkos::deep_copy(Kokkos::subview(w, Kokkos::ALL, CLW_V, Kokkos::ALL, Kokkos::ALL,
                    Kokkos::ALL), 0.0);
  const Real bn = std::sqrt(dot(CLW_B, CLW_B));
  const Real tol = cln_tol*fmax(bn, 1.0e-300);
  Real rn = std::sqrt(dot(CLW_R, CLW_R));
  if (rn <= tol) return 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  for (int it = 1; it <= cln_maxit; ++it) {
    const Real rhon = dot(CLW_R0, CLW_R);
    if (rhon == 0.0) return -it;
    const Real beta = (rhon/rho)*(alpha/omega);
    rho = rhon;
    // p = r + beta (p - omega v)
    par_for("cln_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_P,k,j,i) = w(m,CLW_R,k,j,i) +
                         beta*(w(m,CLW_P,k,j,i) - omega*w(m,CLW_V,k,j,i));
    });
    CleanPrecond(CLW_P, CLW_Y);
    CleanApply(CLW_Y, CLW_V);
    const Real r0v = dot(CLW_R0, CLW_V);
    if (r0v == 0.0) return -it;
    alpha = rho/r0v;
    axpy(CLW_S, CLW_R, 1.0, CLW_V, -alpha);
    const Real sn = std::sqrt(dot(CLW_S, CLW_S));
    if (sn <= tol) {
      axpy(CLW_X, CLW_X, 1.0, CLW_Y, alpha);
      return it;
    }
    CleanPrecond(CLW_S, CLW_Z);
    CleanApply(CLW_Z, CLW_TT);
    const Real tt = dot(CLW_TT, CLW_TT);
    omega = (tt > 0.0) ? dot(CLW_TT, CLW_S)/tt : 0.0;
    par_for("cln_x", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_X,k,j,i) += alpha*w(m,CLW_Y,k,j,i) + omega*w(m,CLW_Z,k,j,i);
      w(m,CLW_R,k,j,i) = w(m,CLW_S,k,j,i) - omega*w(m,CLW_TT,k,j,i);
    });
    rn = std::sqrt(dot(CLW_R, CLW_R));
    if (rn <= tol) return it;
    if (omega == 0.0) return -it;
  }
  return -cln_maxit;
}

//----------------------------------------------------------------------------------------
//! the clean stage solve (replaces ImplicitSolve when <rad_m1>/scheme = clean)
TaskStatus RadiationM1::CleanSolve(Driver *pdrive, int stage) {
  (void)pdrive; (void)stage;
  if (!cln_init) CleanAlloc();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool d2 = pmy_pack->pmesh->multi_d, d3 = pmy_pack->pmesh->three_d;
  const int t2s = t2_solve;
  const bool t2 = (t2s == M1_T2S_STAGE1 || t2s == M1_T2S_STAGE2);
  const Real h = dt_sub, c = c_light, ar = arad;
  auto u0_ = u0;
  auto w = cln_w;
  auto op = opac;
  auto inc = t2inc;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  const bool hh = fl.on;

  // first call: face fluxes from the cell fluxes (restarts keep f0x)
  static bool f0done = false;
  if (!f0done) {
    f0done = true;
    if (pmy_pack->pmesh->ncycle == 0) {
      for (int dir = 0; dir < 3; ++dir) {
        if (dir == 1 && !d2) continue;
        if (dir == 2 && !d3) continue;
        auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
        const int ihi = (dir == 0) ? ie+1 : ie, jhi = (dir == 1) ? je+1 : je;
        const int khi = (dir == 2) ? ke+1 : ke;
        par_for("cln_f0init", DevExeSpace(), 0, nmb1, ks, khi, js, jhi, is, ihi,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const int kl = (dir == 2) ? k-1 : k, jl = (dir == 1) ? j-1 : j;
          const int il = (dir == 0) ? i-1 : i;
          f0(m,k,j,i) = 0.5*(u0_(m,M1_F1+dir,k,j,i) + u0_(m,M1_F1+dir,kl,jl,il));
        });
      }
    }
  }

  const bool fb = gas_feedback && hh;
  const Real gm1 = hh ? eos.gamma - 1.0 : 1.0;
  // gas: the old vector (stages: + t2inc of the gas channels), T, eint; the old E
  par_for("cln_old", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real t = 0.0;
    if (hh) {
      const Real d = uh(m,IDN,k,j,i);
      if (fb && t2) {
        for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += inc(m,M1_T2_M1+n,k,j,i);
        uh(m,IEN,k,j,i) += inc(m,M1_T2_EN,k,j,i);
      }
      const Real kin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                            SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      const Real ei = fmax(uh(m,IEN,k,j,i) - kin, 0.0);
      t = eos.Temperature(d, ei);
      w(m,CLW_EI,k,j,i) = ei;
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,M1_E,k,j,i) + (t2 ? inc(m,M1_T2_E,k,j,i) : 0.0);
  });

  CleanClosure(t2s);
  CleanFaceCoef();

  CleanCtx x;
  x.cs = cln_s; x.w = cln_w; x.fc[0] = cln_fc1; x.fc[1] = cln_fc2; x.fc[2] = cln_fc3;
  x.bmf = ibm_f; x.h = h; x.c = c;
  auto &size = pmy_pack->pmb->mb_size;
  x.idx[0] = 1.0/size.h_view(0).dx1; x.idx[1] = 1.0/size.h_view(0).dx2;
  x.idx[2] = 1.0/size.h_view(0).dx3;
  x.bc[0] = ibc_x1min; x.bc[1] = ibc_x1max; x.bc[2] = ibc_x2min; x.bc[3] = ibc_x2max;
  x.bc[4] = ibc_x3min; x.bc[5] = ibc_x3max;
  x.is = is; x.ie = ie; x.js = js; x.je = je; x.ks = ks; x.ke = ke;
  x.d2 = d2; x.d3 = d3; x.bform = cln_bdry_form; x.bmn = ibm_n; x.q = marshak_q;
  for (int n = 0; n < 6; ++n) x.eb[n] = 0.0;
  x.eb[0] = iebath_x1min; x.eb[1] = iebath_x1max;
  x.fb = fb;
  x.hmus = cln_hmus;

  // gas feedback: enthalpy coefficients a_d = v_d + sum_b D_db v_b (v of the old vector)
  if (fb) {
    par_for("cln_aenth", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
      Real v[3];
      for (int n = 0; n < 3; ++n) v[n] = uh(m,IM1+n,k,j,i)/d;
      const int dmap[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
      for (int a = 0; a < 3; ++a) {
        Real s2 = v[a];
        for (int b = 0; b < 3; ++b) s2 += x.cs(m,CLN_D11+dmap[a][b],k,j,i)*v[b];
        x.w(m,CLW_A1+a,k,j,i) = s2;
      }
    });
    int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
    CleanGhost(cln_w, CLW_A1, 3, nmb1, bc, is, ie, js, je, ks, ke, indcs.ng, d2, d3);
  }
  // the first iterate: E_k = old E, T_k = the exact local T of E_k (fixed T without
  // feedback)
  par_for("cln_it0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real ek = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
    x.w(m,CLW_EK,k,j,i) = ek;
    x.w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) x.w(m,CLW_DM+n,k,j,i) = 0.0;
    if (fb) {
      const Real d = uh(m,IDN,k,j,i);
      x.w(m,CLW_T,k,j,i) = CleanTLocal(d, gm1, x.w(m,CLW_EI,k,j,i),
                                       x.h*x.c*op(m,M1_OP_P,k,j,i)*ar,
                                       x.h*x.c*op(m,M1_OP_E,k,j,i)*ek, x.w(m,CLW_T,k,j,i));
    }
  });
  // G part of each face flux via the zero vector: F(E=0, withg) (CLW_CP is scratch = 0)
  // the raw consistency correction of every interior / periodic face (fct mode)
  x.hrr = cln_hrr;
  x.cfct = !cln_gate && !cln_hrr;
  x.addc = !cln_hrr;
  if (x.cfct) {
    for (int dir = 0; dir < 3; ++dir) {
      if (dir == 1 && !d2) continue;
      if (dir == 2 && !d3) continue;
      const int ihi = (dir == 0) ? ie+1 : ie, jhi = (dir == 1) ? je+1 : je;
      const int khi = (dir == 2) ? ke+1 : ke;
      par_for("cln_corr", DevExeSpace(), 0, nmb1, ks, khi, js, jhi, is, ihi,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
        const int hi = (dir == 0) ? x.ie+1 : ((dir == 1) ? x.je+1 : x.ke+1);
        const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
        const bool bnd = (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) ||
                         (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC);
        x.fc[dir](m,4,k,j,i) = bnd ? 0.0 : CleanCorr(x, dir, m, k, j, i);
        x.fc[dir](m,3,k,j,i) = x.fc[dir](m,4,k,j,i);
      });
    }
  }
  int nit = 0, nnl = 0;
  bool nlok = !fb;
  const int maxnl = fb ? cln_maxnl : 1;
  for (int inl = 0; inl < maxnl; ++inl) {
    ++nnl;
    Kokkos::deep_copy(Kokkos::subview(w, Kokkos::ALL, CLW_CP, Kokkos::ALL, Kokkos::ALL,
                      Kokkos::ALL), 0.0);
    // r, B = EO + h c kP a T_k^4 - r h c kE E_k - h div G.  With the positivity limit
    // of the consistency correction: B0 without it, then each face's correction scaled so
    // that no cell gives away more than B0 (the M-matrix keeps E >= 0 for B >= 0)
    for (int prh = 0; prh < (x.cfct ? 2 : 1); ++prh) {
    x.addc = !cln_hrr && !(x.cfct && prh == 0);
    par_for("cln_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real t = x.w(m,CLW_T,k,j,i);
      const Real kp = op(m,M1_OP_P,k,j,i), ke2 = op(m,M1_OP_E,k,j,i);
      Real rr = 0.0;
      if (fb) {
        const Real s4 = 4.0*x.h*x.c*kp*ar*t*t*t;
        rr = s4/(uh(m,IDN,k,j,i)/gm1 + s4);
      }
      x.w(m,CLW_RR,k,j,i) = rr;
      Real b = x.w(m,CLW_EO,k,j,i) + x.h*x.c*kp*ar*SQR(SQR(t)) -
               rr*x.h*x.c*ke2*x.w(m,CLW_EK,k,j,i);
      const int z = CLW_CP;
      b -= x.h*x.idx[0]*(CleanFlux(x,z,0,m,k,j,i+1,true) - CleanFlux(x,z,0,m,k,j,i,true));
      if (x.d2) b -= x.h*x.idx[1]*(CleanFlux(x,z,1,m,k,j+1,i,true) -
                                   CleanFlux(x,z,1,m,k,j,i,true));
      if (x.d3) b -= x.h*x.idx[2]*(CleanFlux(x,z,2,m,k+1,j,i,true) -
                                   CleanFlux(x,z,2,m,k,j,i,true));
      x.w(m,CLW_B,k,j,i) = b;
    });
    if (x.cfct && prh == 0) {
      par_for("cln_fctr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real pout = x.idx[0]*(fmax(x.fc[0](m,4,k,j,i+1), 0.0) + fmax(-x.fc[0](m,4,k,j,i), 0.0));
        if (x.d2) pout += x.idx[1]*(fmax(x.fc[1](m,4,k,j+1,i), 0.0) +
                                    fmax(-x.fc[1](m,4,k,j,i), 0.0));
        if (x.d3) pout += x.idx[2]*(fmax(x.fc[2](m,4,k+1,j,i), 0.0) +
                                    fmax(-x.fc[2](m,4,k,j,i), 0.0));
        pout *= x.h;
        x.w(m,CLW_TT,k,j,i) = (pout > 0.0) ? fmin(1.0, fmax(x.w(m,CLW_B,k,j,i), 0.0)/pout)
                                           : 1.0;
      });
      CleanFillGhost(CLW_TT);
      for (int dir = 0; dir < 3; ++dir) {
        if (dir == 1 && !d2) continue;
        if (dir == 2 && !d3) continue;
        const int ihi = (dir == 0) ? ie+1 : ie, jhi = (dir == 1) ? je+1 : je;
        const int khi = (dir == 2) ? ke+1 : ke;
        par_for("cln_fctf", DevExeSpace(), 0, nmb1, ks, khi, js, jhi, is, ihi,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const int kl = (dir == 2) ? k-1 : k, jl = (dir == 1) ? j-1 : j;
          const int il = (dir == 0) ? i-1 : i;
          const Real d = x.fc[dir](m,4,k,j,i);
          x.fc[dir](m,3,k,j,i) = d*((d > 0.0) ? x.w(m,CLW_TT,kl,jl,il)
                                              : x.w(m,CLW_TT,k,j,i));
        });
      }
    }
    }
    // preconditioner: the operator's x1 tridiagonal part and the x2/x3 diagonal parts
    // (exact for the diagonal-D, half-range and enthalpy parts; the cross terms are
    // left to the Krylov iteration)
    par_for("cln_pc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real hx = x.h*x.idx[0];
      auto coef = [&](const int dir, const int kk2, const int jj, const int ii, Real &al,
                      Real &arr) {
        const int kl = (dir == 2) ? kk2-1 : kk2, jl = (dir == 1) ? jj-1 : jj;
        const int il = (dir == 0) ? ii-1 : ii;
        const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
        const int hi = (dir == 0) ? x.ie+1 : ((dir == 1) ? x.je+1 : x.ke+1);
        const int pos = (dir == 0) ? ii : ((dir == 1) ? jj : kk2);
        al = 0.0; arr = 0.0;
        if (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) {
          if (x.bc[2*dir] == M1_IBC_REFLECT) return;
          bool beam = false;
          if (dir == 0 && x.bmn > 0) beam = (x.bmf(m,kk2-x.ks,jj-x.js) > 0.0);
          if (x.bform == 1) {arr = beam ? 0.0 : -x.c*x.q; return;}
          arr = x.c*x.cs(m,CLN_HM+dir,kk2,jj,ii);
          return;
        }
        if (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC) {
          if (x.bc[2*dir+1] == M1_IBC_REFLECT) return;
          al = (x.bform == 1) ? x.c*x.q : x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
          return;
        }
        const Real wf = x.fc[dir](m,0,kk2,jj,ii), th = x.fc[dir](m,1,kk2,jj,ii);
        const Real g = (1.0 - wf)*th*x.h*x.c*x.c*x.idx[dir];
        const Real hpf = x.hrr ? x.fc[dir](m,3,kk2,jj,ii) : x.cs(m,CLN_HP+dir,kl,jl,il);
        const Real hmf = x.hrr ? x.fc[dir](m,4,kk2,jj,ii) : x.cs(m,CLN_HM+dir,kk2,jj,ii);
        al = g*x.cs(m,CLN_D11+dir,kl,jl,il) + wf*x.c*hpf;
        arr = -g*x.cs(m,CLN_D11+dir,kk2,jj,ii) + wf*x.c*hmf;
        if (x.fb) {
          al += (1.0 - wf)*(1.0 - th)*0.5*x.w(m,CLW_A1+dir,kl,jl,il);
          arr += (1.0 - wf)*(1.0 - th)*0.5*x.w(m,CLW_A1+dir,kk2,jj,ii);
        }
      };
      Real alp, arp, alm, arm;
      coef(0, k, j, i+1, alp, arp);
      coef(0, k, j, i, alm, arm);
      Real a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - x.w(m,CLW_RR,k,j,i)) +
                hx*(alp - arm);
      x.w(m,CLW_AM,k,j,i) = -hx*alm;
      x.w(m,CLW_AP,k,j,i) = hx*arp;
      if (x.d2) {
        coef(1, k, j+1, i, alp, arp);
        coef(1, k, j, i, alm, arm);
        a0 += x.h*x.idx[1]*(alp - arm);
      }
      if (x.d3) {
        coef(2, k+1, j, i, alp, arp);
        coef(2, k, j, i, alm, arm);
        a0 += x.h*x.idx[2]*(alp - arm);
      }
      x.w(m,CLW_A0,k,j,i) = a0;
      x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
    });
    nit = CleanBiCGStab();
    cln_nsolve += 1.0;
    cln_nit += std::abs(nit);
    cln_maxitseen = fmax(cln_maxitseen, static_cast<Real>(std::abs(nit)));
    if (!fb || nit < 0) break;
    // Newton update: E_{k+1} = X, damped per cell so that E stays positive; T_{k+1} the
    // exact local T of E_{k+1}.  Converged when max |dE|/E and max |dT|/T < cln_tol_nl.
    Real dmax = 0.0;
    {
      const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
      const int nkji = nk*nj*ni, nji = nj*ni;
      Kokkos::parallel_reduce("cln_nlupd", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1+1)*nkji),
      KOKKOS_LAMBDA(const int idx, Real &acc) {
        int m = idx/nkji, r = idx - m*nkji;
        int k = r/nji; r -= k*nji;
        int j = r/ni, i = r - j*ni;
        k += ks; j += js; i += is;
        const Real ek = x.w(m,CLW_EK,k,j,i);
        Real en = x.w(m,CLW_X,k,j,i);
        if (!(en > 0.0)) en = 0.1*ek;
        const Real d = uh(m,IDN,k,j,i);
        const Real told = x.w(m,CLW_T,k,j,i);
        const Real tn = CleanTLocal(d, gm1, x.w(m,CLW_EI,k,j,i),
                                    x.h*x.c*op(m,M1_OP_P,k,j,i)*ar,
                                    x.h*x.c*op(m,M1_OP_E,k,j,i)*en, told);
        x.w(m,CLW_EK,k,j,i) = en;
        x.w(m,CLW_T,k,j,i) = tn;
        const Real de = fabs(en - ek)/fmax(fmax(en, ek), 1.0e-300);
        const Real dtt = fabs(tn - told)/fmax(fmax(tn, told), 1.0e-300);
        acc = fmax(acc, fmax(de, dtt));
      }, Kokkos::Max<Real>(dmax));
    }
    // implicit velocity: dm of the new iterate, v = (m_old + dm)/rho, new a_d
    CleanFillGhost(CLW_EK);
    Real dvmax = 0.0, vmax = 0.0;
    {
      const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
      const int nkji = nk*nj*ni, nji = nj*ni;
      Kokkos::parallel_reduce("cln_dm", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1+1)*nkji),
      KOKKOS_LAMBDA(const int idx, Real &acc) {
        int m = idx/nkji, r = idx - m*nkji;
        int k = r/nji; r -= k*nji;
        int j = r/ni, i = r - j*ni;
        k += ks; j += js; i += is;
        Real dmn[3] = {0.0, 0.0, 0.0};
        dmn[0] = CleanCellMom(x, CLW_EK, m, k, j, i, 0);
        if (x.d2) dmn[1] = CleanCellMom(x, CLW_EK, m, k, j, i, 1);
        if (x.d3) dmn[2] = CleanCellMom(x, CLW_EK, m, k, j, i, 2);
        const Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
        Real dv = 0.0;
        for (int n = 0; n < 3; ++n) {
          dv = fmax(dv, fabs(dmn[n] - x.w(m,CLW_DM+n,k,j,i))/d);
          x.w(m,CLW_DM+n,k,j,i) = dmn[n];
        }
        acc = fmax(acc, dv);
      }, Kokkos::Max<Real>(dvmax));
      Kokkos::parallel_reduce("cln_vmax", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1+1)*nkji),
      KOKKOS_LAMBDA(const int idx, Real &acc) {
        int m = idx/nkji, r = idx - m*nkji;
        int k = r/nji; r -= k*nji;
        int j = r/ni, i = r - j*ni;
        k += ks; j += js; i += is;
        const Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
        Real v2 = 0.0;
        for (int n = 0; n < 3; ++n) v2 += SQR((uh(m,IM1+n,k,j,i) + x.w(m,CLW_DM+n,k,j,i))/d);
        acc = fmax(acc, std::sqrt(v2));
      }, Kokkos::Max<Real>(vmax));
    }
    par_for("cln_anew", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
      Real v[3];
      for (int n = 0; n < 3; ++n) v[n] = (uh(m,IM1+n,k,j,i) + x.w(m,CLW_DM+n,k,j,i))/d;
      const int dmap[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
      for (int a = 0; a < 3; ++a) {
        Real s2 = v[a];
        for (int b = 0; b < 3; ++b) s2 += x.cs(m,CLN_D11+dmap[a][b],k,j,i)*v[b];
        x.w(m,CLW_A1+a,k,j,i) = s2;
      }
    });
    {
      int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
      CleanGhost(cln_w, CLW_A1, 3, nmb1, bc, is, ie, js, je, ks, ke, indcs.ng, d2, d3);
    }
    if (dmax < cln_tol_nl && dvmax <= cln_tol_nl*fmax(vmax, 1.0e-300)) {nlok = true; break;}
  }
  if (fb) {
    // the solution is the converged iterate
    par_for("cln_xek", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
    });
    if (!nlok) nit = -std::abs(nit) - 1;
  }

  // positivity check of the solution.  The linear solve is accurate to cln_tol |B|, so
  // values within 10 cln_tol max|B| of zero are zero to solver accuracy (near-vacuum
  // cells next to a beam): they are set to 0.  Anything below that is a negative.
  Real bmax = 0.0;
  {
    const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
    const int nkji = nk*nj*ni, nji = nj*ni;
    Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1+1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &acc) {
      int m = idx/nkji, r = idx - m*nkji;
      int k = r/nji; r -= k*nji;
      int j = r/ni, i = r - j*ni;
      acc = fmax(acc, fabs(w(m,CLW_B,k+ks,j+js,i+is)));
    }, Kokkos::Max<Real>(bmax));
  }
  const Real thr = 10.0*cln_tol*bmax;
  par_for("cln_zero", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = w(m,CLW_X,k,j,i);
    if (e < 0.0 && e >= -thr) w(m,CLW_X,k,j,i) = 0.0;
  });
  Real nneg = 0.0;
  {
    const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
    const int nkji = nk*nj*ni, nji = nj*ni;
    Kokkos::parallel_reduce("cln_neg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1+1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &acc) {
      int m = idx/nkji, r = idx - m*nkji;
      int k = r/nji; r -= k*nji;
      int j = r/ni, i = r - j*ni;
      if (!(w(m,CLW_X,k+ks,j+js,i+is) >= 0.0)) acc += 1.0;
    }, nneg);
  }
  if (nneg > 0.0 && global_variable::my_rank == 0) {
    // diagnostic: the first negative cells
    auto hw = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_w);
    auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_s);
    auto hf = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_fc1);
    auto hg = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_fc2);
    int np = 0;
    for (int k = ks; k <= ke && np < 4; ++k) {
      for (int j = js; j <= je && np < 4; ++j) {
        for (int i = is; i <= ie && np < 4; ++i) {
          if (hw(0,CLW_X,k,j,i) >= 0.0) continue;
          ++np;
          std::printf("CLEANNEG j %d i %d EO %.3e B %.3e X %.3e C1 %.3f w1 %.2f/%.2f "
                      "w2 %.2f/%.2f G1 %.3e/%.3e G2 %.3e/%.3e hp1 %.3f hm1 %.3f D11 %.3f\n",
                      j-js, i-is, hw(0,CLW_EO,k,j,i), hw(0,CLW_B,k,j,i), hw(0,CLW_X,k,j,i),
                      hs(0,CLN_C1,k,j,i), hf(0,0,k,j,i), hf(0,0,k,j,i+1),
                      d2 ? hg(0,0,k,j,i) : 0.0, d2 ? hg(0,0,k,j+1,i) : 0.0,
                      hf(0,2,k,j,i), hf(0,2,k,j,i+1), d2 ? hg(0,2,k,j,i) : 0.0,
                      d2 ? hg(0,2,k,j+1,i) : 0.0, hs(0,CLN_HP,k,j,i), hs(0,CLN_HM,k,j,i),
                      hs(0,CLN_D11,k,j,i));
        }
      }
    }
  }
  if (nit < 0) cln_nfail += 1.0;
  if (nneg > 0.0) cln_nneg += 1.0;
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN cycle %d t2s %d it %d neg %.0f\n", pmy_pack->pmesh->ncycle, t2s,
                nit, nneg);
  }

  // write-back: faces (full flux of the solution), E, cell F, slopes
  CleanFillGhost(CLW_X);
  const bool sl = (t2s != M1_T2S_NONE);
  auto kk = (t2s == M1_T2S_STAGE1) ? t2k2 : t2k1;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !d2) continue;
    if (dir == 2 && !d3) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    auto fcd = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int ihi = (dir == 0) ? ie+1 : ie, jhi = (dir == 1) ? je+1 : je;
    const int khi = (dir == 2) ? ke+1 : ke;
    const int tf = M1_T2_F1 + dir;
    par_for("cln_fwb", DevExeSpace(), 0, nmb1, ks, khi, js, jhi, is, ihi,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real fn = CleanFlux(x, CLW_X, dir, m, k, j, i, true, false);
      if (sl) kk(m,tf,k,j,i) = (fn - fcd(m,2,k,j,i))/x.h;
      f0(m,k,j,i) = fn;
    });
  }
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  auto fc1 = cln_fc1, fc2 = cln_fc2, fc3 = cln_fc3;
  par_for("cln_wb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = w(m,CLW_X,k,j,i);
    Real dm[3] = {0.0, 0.0, 0.0}, den = 0.0;
    if (fb) {
      // momentum: the face-flux relaxation h c chi F (h chi F/c of momentum), half of
      // each face to each cell; chi_f from theta = 1/(1 + h c chi_f)
      // comoving flux F - (1 - w) A_f, A_f = (a_L E_L + a_R E_R)/2 (lab-frame form)
      auto mom = [&](const Real th, const Real wf, const Real f, const Real af) {
        return (1.0/th - 1.0)*(f - (1.0 - wf)*af)/(c*c);};
      auto afc = [&](const int a, const int kl, const int jl, const int il) {
        return 0.5*(w(m,CLW_A1+a,kl,jl,il)*w(m,CLW_X,kl,jl,il) +
                    w(m,CLW_A1+a,k,j,i)*w(m,CLW_X,k,j,i));};
      auto afh = [&](const int a, const int kh, const int jh, const int ih) {
        return 0.5*(w(m,CLW_A1+a,kh,jh,ih)*w(m,CLW_X,kh,jh,ih) +
                    w(m,CLW_A1+a,k,j,i)*w(m,CLW_X,k,j,i));};
      dm[0] = 0.5*(mom(fc1(m,1,k,j,i), fc1(m,0,k,j,i), f1(m,k,j,i), afc(0,k,j,i-1)) +
                   mom(fc1(m,1,k,j,i+1), fc1(m,0,k,j,i+1), f1(m,k,j,i+1), afh(0,k,j,i+1)));
      if (d2) dm[1] = 0.5*(mom(fc2(m,1,k,j,i), fc2(m,0,k,j,i), f2(m,k,j,i), afc(1,k,j-1,i)) +
                           mom(fc2(m,1,k,j+1,i), fc2(m,0,k,j+1,i), f2(m,k,j+1,i),
                               afh(1,k,j+1,i)));
      if (d3) dm[2] = 0.5*(mom(fc3(m,1,k,j,i), fc3(m,0,k,j,i), f3(m,k,j,i), afc(2,k-1,j,i)) +
                           mom(fc3(m,1,k+1,j,i), fc3(m,0,k+1,j,i), f3(m,k+1,j,i),
                               afh(2,k+1,j,i)));
      const Real d = uh(m,IDN,k,j,i);
      const Real ien0 = uh(m,IEN,k,j,i);
      for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += dm[n];
      const Real kin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                            SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      uh(m,IEN,k,j,i) = d*w(m,CLW_T,k,j,i)/gm1 + kin;
      den = uh(m,IEN,k,j,i) - ien0;
    }
    if (sl) {
      kk(m,M1_T2_E,k,j,i) = (e - w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,M1_T2_M1+n,k,j,i) = dm[n]/h;
      kk(m,M1_T2_EN,k,j,i) = den/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    u0_(m,M1_F1,k,j,i) = 0.5*(f1(m,k,j,i) + f1(m,k,j,i+1));
    u0_(m,M1_F1+1,k,j,i) = d2 ? 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i)) : 0.0;
    u0_(m,M1_F1+2,k,j,i) = d3 ? 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i)) : 0.0;
  });
  if (t2) t2_fail = (nit < 0) || (nneg > 0.0);
  (void)c; (void)op;
  return TaskStatus::complete;
}

void RadiationM1::CleanReport() {
  if (global_variable::my_rank != 0) return;
  std::printf("CLEAN summary: solves %.0f iterations %.0f (max %.0f) failed %.0f "
              "negative %.0f\n", cln_nsolve, cln_nit, cln_maxitseen, cln_nfail, cln_nneg);
}

}  // namespace radm1
