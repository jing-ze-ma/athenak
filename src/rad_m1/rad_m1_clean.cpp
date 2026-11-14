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
//! The closure (D, h+-, C1) is built once per step at U^n and extrapolated linearly to
//! the stage time from the previous build.  Gate-1 scope: Cartesian, one MeshBlock, one
//! rank, chat = c, gas_feedback = false (static gas, fixed T).

#include <cmath>
#include <cstdio>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/m1_fluid.hpp"

namespace radm1 {

namespace {
// closure slots (cln_c / cln_p / cln_s)
constexpr int CLN_D11 = 0;   // D11 D22 D33 D12 D13 D23
constexpr int CLN_HP = 6;    // h+_a = HP_a/J, a = 1 2 3
constexpr int CLN_HM = 9;    // h-_a = H_a/J - h+_a (<= 0)
constexpr int CLN_C1 = 12;   // C1 = NL/J
constexpr int CLN_N = 13;
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
constexpr int CLW_N = 16;

struct CleanCtx {
  DvceArray5D<Real> cs;          // stage closure
  DvceArray5D<Real> w;           // work
  DvceArray5D<Real> fc[3];       // face coefficients (m,3,k,j,i): w, theta, G
  DvceArray3D<Real> bmf;         // beam face flux at x1min (m, k-ks, j-js)
  Real h, c, idx[3];
  int bc[6];                     // M1_IBC_* x1min x1max x2min x2max x3min x3max
  int is, ie, js, je, ks, ke;
  bool d2, d3;
  int bform;                     // 0 half-range rays, 1 Marshak q
  int bmn;
  Real q, eb[6];
};

KOKKOS_INLINE_FUNCTION Real CleanSmooth(Real x) {
  x = fmin(fmax(x, 0.0), 1.0);
  return x*x*(3.0 - 2.0*x);
}

// face flux of direction dir (0,1,2) at face (k,j,i) (face between cell-1 and cell along
// dir) for the E vector in work slot e.  withg: include the E-independent part G.
KOKKOS_INLINE_FUNCTION Real CleanFlux(const CleanCtx &x, const int e, const int dir,
                                      const int m, const int k, const int j, const int i,
                                      const bool withg) {
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
  }
  if (wf > 0.0) {
    f += wf*c*(x.cs(m,CLN_HP+dir,kl,jl,il)*el + x.cs(m,CLN_HM+dir,k,j,i)*er);
  }
  return f;
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
  Kokkos::realloc(cln_p, nmb, CLN_N, n3, n2, n1);
  Kokkos::realloc(cln_s, nmb, CLN_N, n3, n2, n1);
  Kokkos::realloc(cln_w, nmb, CLW_N, n3, n2, n1);
  Kokkos::realloc(cln_fc1, nmb, 3, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, 3, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, 3, n3, n2, n1);
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
      chat != c_light || gas_feedback) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
      << "<rad_m1>/scheme = clean (gate 1) needs: one MeshBlock, one rank, a Cartesian "
      << "mesh, closure = vet_sc with vet_tensor = full, chat = c, gas_feedback = false"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the formal solution must accumulate the half-range sums and C1 for this scheme
  impl_beam_hr = true;
  thsw_c1 = true;
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! closure: formal solution at U^n once per step (STAGE1 or a BE solve), copied to
//! cln_c; the stage closure cln_s = c + r (c - p), r = dt/(t_c - t_p), per cell kept
//! realisable (else c).  STAGE2 keeps the stage-1 cln_s.
void RadiationM1::CleanClosure(int t2s) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, ng = indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool d2 = pmy_pack->pmesh->multi_d, d3 = pmy_pack->pmesh->three_d;
  if (t2s == M1_T2S_STAGE2) return;
  const Real tnow = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  // inputs of the formal solution: E^n, T^n, F^n of the step start
  {
    auto iw_ = iw;
    auto u0_ = u0;
    auto cw_ = cln_w;
    par_for("cln_vin", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_EN,k,j,i) = u0_(m,M1_E,k,j,i);
      iw_(m,M1_IW_TP,k,j,i) = cw_(m,CLW_T,k,j,i);
      iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
      iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F1+1,k,j,i);
      iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F1+2,k,j,i);
    });
  }
  VetShortChar();
  const bool newstep = (tnow > cln_tc);
  if (newstep && cln_havec) {
    Kokkos::deep_copy(cln_p, cln_c);
    cln_havep = true;
    cln_dtp = tnow - cln_tc;
  }
  cln_tc = tnow;
  cln_havec = true;
  {
    auto vc = vet_cell;
    auto cc = cln_c;
    par_for("cln_vcopy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jj = fmax(vc(m,M1_VET_J,k,j,i), 1.0e-300);
      for (int a = 0; a < 6; ++a) cc(m,CLN_D11+a,k,j,i) = vc(m,M1_VET_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,M1_VET_HP1+a,k,j,i)/jj;
        const Real hh = vc(m,M1_VET_H1+a,k,j,i)/jj;
        cc(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cc(m,CLN_HM+a,k,j,i) = fmin(hh - hp, 0.0);
      }
      // no radiation reaches the cell (J = 0): no angular information, take C1 = 1 (the
      // half-range, upwind face), never the central one
      const Real jr = vc(m,M1_VET_J,k,j,i);
      cc(m,CLN_C1,k,j,i) = (jr > 0.0) ? fmax(vc(m,M1_VET_NL,k,j,i)/jr, 0.0) : 1.0;
    });
  }
  const bool ext = cln_extrap && cln_havep && (t2s == M1_T2S_STAGE1) && cln_dtp > 0.0;
  const Real r = ext ? dt/cln_dtp : 0.0;
  {
    auto cc = cln_c;
    auto cp = cln_p;
    auto cs = cln_s;
    par_for("cln_vext", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real v[CLN_N];
      bool ok = true;
      for (int n = 0; n < CLN_N; ++n) {
        v[n] = cc(m,n,k,j,i) + r*(cc(m,n,k,j,i) - cp(m,n,k,j,i));
      }
      if (r > 0.0) {
        // realisability: positive diagonal, |D_ab| <= sqrt(D_aa D_bb), trace <= 1,
        // h+ >= 0 >= h-, C1 >= 0
        for (int a = 0; a < 3; ++a) {
          if (!(v[CLN_D11+a] > 0.0)) ok = false;
          if (v[CLN_HP+a] < 0.0 || v[CLN_HM+a] > 0.0) ok = false;
        }
        if (v[0] + v[1] + v[2] > 1.0 + 1.0e-12) ok = false;
        if (ok) {
          if (v[3]*v[3] > v[0]*v[1] || v[4]*v[4] > v[0]*v[2] || v[5]*v[5] > v[1]*v[2]) {
            ok = false;
          }
        }
        if (v[CLN_C1] < 0.0) v[CLN_C1] = 0.0;
      }
      for (int n = 0; n < CLN_N; ++n) cs(m,n,k,j,i) = ok ? v[n] : cc(m,n,k,j,i);
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
  auto op = opac;
  par_for("cln_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real a = x.w(m,src,k,j,i)*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i));
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

  // gas temperature (fixed: gas_feedback = false) and the old vector of E
  par_for("cln_old", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real t = 0.0;
    if (hh) {
      const Real d = uh(m,IDN,k,j,i);
      const Real kin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                            SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      t = eos.Temperature(d, fmax(uh(m,IEN,k,j,i) - kin, 0.0));
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,M1_E,k,j,i) + (t2 ? inc(m,M1_T2_E,k,j,i) : 0.0);
  });

  CleanClosure(t2s);
  CleanFaceCoef();

  // right-hand side B = EO + h c kP a T^4 - h div G, and the preconditioner coefficients
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
  // G part of each face flux via the zero vector: F(E=0, withg) (CLW_CP is scratch = 0)
  Kokkos::deep_copy(Kokkos::subview(w, Kokkos::ALL, CLW_CP, Kokkos::ALL, Kokkos::ALL,
                    Kokkos::ALL), 0.0);
  par_for("cln_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = x.w(m,CLW_T,k,j,i);
    Real b = x.w(m,CLW_EO,k,j,i) + x.h*x.c*op(m,M1_OP_P,k,j,i)*ar*SQR(SQR(t));
    const int z = CLW_CP;
    b -= x.h*x.idx[0]*(CleanFlux(x,z,0,m,k,j,i+1,true) - CleanFlux(x,z,0,m,k,j,i,true));
    if (x.d2) b -= x.h*x.idx[1]*(CleanFlux(x,z,1,m,k,j+1,i,true) -
                                 CleanFlux(x,z,1,m,k,j,i,true));
    if (x.d3) b -= x.h*x.idx[2]*(CleanFlux(x,z,2,m,k+1,j,i,true) -
                                 CleanFlux(x,z,2,m,k,j,i,true));
    x.w(m,CLW_B,k,j,i) = b;
  });
  // preconditioner: the operator's x1 tridiagonal part and the x2/x3 diagonal parts,
  // measured by applying the linear face flux to unit vectors (exact for the diagonal-D
  // and half-range parts; the cross terms are left to the Krylov iteration)
  par_for("cln_pc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real hx = x.h*x.idx[0];
    // alpha coefficients of face f: F_f = aL E_L + aR E_R (diagonal D + half range)
    auto coef = [&](const int dir, const int kk, const int jj, const int ii, Real &al,
                    Real &arr) {
      const int kl = (dir == 2) ? kk-1 : kk, jl = (dir == 1) ? jj-1 : jj;
      const int il = (dir == 0) ? ii-1 : ii;
      const int lo = (dir == 0) ? x.is : ((dir == 1) ? x.js : x.ks);
      const int hi = (dir == 0) ? x.ie+1 : ((dir == 1) ? x.je+1 : x.ke+1);
      const int pos = (dir == 0) ? ii : ((dir == 1) ? jj : kk);
      al = 0.0; arr = 0.0;
      if (pos == lo && x.bc[2*dir] != M1_IBC_PERIODIC) {
        if (x.bc[2*dir] == M1_IBC_REFLECT) return;
        bool beam = false;
        if (dir == 0 && x.bmn > 0) beam = (x.bmf(m,kk-x.ks,jj-x.js) > 0.0);
        if (x.bform == 1) {arr = beam ? 0.0 : -x.c*x.q; return;}
        arr = x.c*x.cs(m,CLN_HM+dir,kk,jj,ii);
        return;
      }
      if (pos == hi && x.bc[2*dir+1] != M1_IBC_PERIODIC) {
        if (x.bc[2*dir+1] == M1_IBC_REFLECT) return;
        al = (x.bform == 1) ? x.c*x.q : x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
        return;
      }
      const Real wf = x.fc[dir](m,0,kk,jj,ii), th = x.fc[dir](m,1,kk,jj,ii);
      const Real g = (1.0 - wf)*th*x.h*x.c*x.c*x.idx[dir];
      al = g*x.cs(m,CLN_D11+dir,kl,jl,il) + wf*x.c*x.cs(m,CLN_HP+dir,kl,jl,il);
      arr = -g*x.cs(m,CLN_D11+dir,kk,jj,ii) + wf*x.c*x.cs(m,CLN_HM+dir,kk,jj,ii);
    };
    Real alp, arp, alm, arm;
    coef(0, k, j, i+1, alp, arp);
    coef(0, k, j, i, alm, arm);
    Real a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i) + hx*(alp - arm);
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
    x.w(m,CLW_X,k,j,i) = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
  });

  const int nit = CleanBiCGStab();
  cln_nsolve += 1.0;
  cln_nit += std::abs(nit);
  cln_maxitseen = fmax(cln_maxitseen, static_cast<Real>(std::abs(nit)));

  // positivity check of the solution
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
      if (!(w(m,CLW_X,k+ks,j+js,i+is) > 0.0)) acc += 1.0;
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
          if (hw(0,CLW_X,k,j,i) > 0.0) continue;
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
      const Real fn = CleanFlux(x, CLW_X, dir, m, k, j, i, true);
      if (sl) kk(m,tf,k,j,i) = (fn - fcd(m,2,k,j,i))/x.h;
      f0(m,k,j,i) = fn;
    });
  }
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  par_for("cln_wb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = w(m,CLW_X,k,j,i);
    if (sl) {
      kk(m,M1_T2_E,k,j,i) = (e - w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,M1_T2_M1+n,k,j,i) = 0.0;
      kk(m,M1_T2_EN,k,j,i) = 0.0;
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
