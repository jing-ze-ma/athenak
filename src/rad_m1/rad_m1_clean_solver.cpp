//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_solver.cpp
//! \brief Linear algebra of the clean implicit radiation solve: face coefficients, the
//! operator A E = E (1 + h c kE (1-r)) + h div F_lin(E), an x1-line tridiagonal
//! preconditioner and right-preconditioned BiCGStab.

#include <cmath>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_clean.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn CleanArgs RadiationM1::CleanMakeArgs
//! \brief collects the arrays and constants the device functions need

CleanArgs RadiationM1::CleanMakeArgs() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  CleanArgs x;
  x.cs = cln_s;
  x.w = cln_w;
  x.fc[0] = cln_fc1;
  x.fc[1] = cln_fc2;
  x.fc[2] = cln_fc3;
  x.beam = ibm_f;
  x.h = dt_sub;
  x.c = c_light;
  x.idx[0] = 1.0/size.h_view(0).dx1;
  x.idx[1] = 1.0/size.h_view(0).dx2;
  x.idx[2] = 1.0/size.h_view(0).dx3;
  x.bc[0] = ibc_x1min;
  x.bc[1] = ibc_x1max;
  x.bc[2] = ibc_x2min;
  x.bc[3] = ibc_x2max;
  x.bc[4] = ibc_x3min;
  x.bc[5] = ibc_x3max;
  for (int n = 0; n < 6; ++n) x.ebath[n] = 0.0;
  x.ebath[0] = iebath_x1min;
  x.ebath[1] = iebath_x1max;
  x.marshak_q = marshak_q;
  x.marshak = (cln_bdry_form == 1);
  x.use_beam = (ibm_n > 0);
  x.feedback = cln_feedback;
  x.is = indcs.is;
  x.ie = indcs.ie;
  x.js = indcs.js;
  x.je = indcs.je;
  x.ks = indcs.ks;
  x.ke = indcs.ke;
  x.multi_d = pmy_pack->pmesh->multi_d;
  x.three_d = pmy_pack->pmesh->three_d;
  return x;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGhostFill
//! \brief ghost cells of components [c0, c0+nc) of a work array: periodic wrap, else a
//! copy of the edge cell

void RadiationM1::CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool per1 = (ibc_x1min == M1_IBC_PERIODIC);
  const bool per2 = (ibc_x2min == M1_IBC_PERIODIC);
  const bool per3 = (ibc_x3min == M1_IBC_PERIODIC);

  par_for("cln_ghost1", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
    a(m,c0+n,k,j,is-1-g) = per1 ? a(m,c0+n,k,j,ie-g) : a(m,c0+n,k,j,is);
    a(m,c0+n,k,j,ie+1+g) = per1 ? a(m,c0+n,k,j,is+g) : a(m,c0+n,k,j,ie);
  });
  if (pmy_pack->pmesh->multi_d) {
    par_for("cln_ghost2", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, 0, ng-1, is-ng, ie+ng,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int g, const int i) {
      a(m,c0+n,k,js-1-g,i) = per2 ? a(m,c0+n,k,je-g,i) : a(m,c0+n,k,js,i);
      a(m,c0+n,k,je+1+g,i) = per2 ? a(m,c0+n,k,js+g,i) : a(m,c0+n,k,je,i);
    });
  }
  if (pmy_pack->pmesh->three_d) {
    par_for("cln_ghost3", DevExeSpace(), 0, nmb1, 0, nc-1, 0, ng-1, js-ng, je+ng,
            is-ng, ie+ng,
    KOKKOS_LAMBDA(const int m, const int n, const int g, const int j, const int i) {
      a(m,c0+n,ks-1-g,j,i) = per3 ? a(m,c0+n,ke-g,j,i) : a(m,c0+n,ks,j,i);
      a(m,c0+n,ke+1+g,j,i) = per3 ? a(m,c0+n,ks+g,j,i) : a(m,c0+n,ke,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanDot
//! \brief sum over active cells of w(a) w(b)

Real RadiationM1::CleanDot(const int a, const int b) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w = cln_w;
  Real sum = 0.0;
  Kokkos::parallel_reduce("cln_dot", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &s) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    s += w(m,a,k+ks,j+js,i+is)*w(m,b,k+ks,j+js,i+is);
  }, Kokkos::Sum<Real>(sum));
  return sum;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFaceCoefficients
//! \brief per face: the blend weight w(C1), theta, the old-vector flux G, and the
//! half-range factors h+-_f.  h+_f = H+_f/J_L and h-_f = H-_f/J_R, with H+-_f the rays'
//! half-range fluxes reconstructed to the face from the upwind side (minmod-limited);
//! |h+-_f| <= 1 keeps the operator an M-matrix.

void RadiationM1::CleanFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = M1_T2_F1 + dir;

    par_for("cln_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real c1 = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      fc(m,CLF_W,k,j,i) = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_G,k,j,i) = f0(m,k,j,i) + (stage ? inc(m,tf,k,j,i) : 0.0);

      // rays' half-range fluxes H+ = h+ J (from the left) and H- = h- J (from the right)
      auto hp = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HP+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto hm = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HM+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto minmod = [](const Real a, const Real b) {
        return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);
      };
      const Real hpl = hp(kl,jl,il);
      const Real hpf = hpl + 0.5*minmod(hpl - hp(kl-dk,jl-dj,il-di), hp(k,j,i) - hpl);
      const Real hmr = hm(k,j,i);
      const Real hmf = hmr - 0.5*minmod(hmr - hm(kl,jl,il), hm(k+dk,j+dj,i+di) - hmr);
      const Real jl_ = cs(m,CLN_J,kl,jl,il), jr_ = cs(m,CLN_J,k,j,i);
      fc(m,CLF_HP,k,j,i) = (jl_ > 0.0) ? fmin(fmax(hpf/jl_, 0.0), 1.0)
                                       : cs(m,CLN_HP+dir,kl,jl,il);
      fc(m,CLF_HM,k,j,i) = (jr_ > 0.0) ? fmax(fmin(hmf/jr_, 0.0), -1.0)
                                       : cs(m,CLN_HM+dir,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanApply
//! \brief w(dst) = A w(src) over the active cells (the linear part of the operator)

void RadiationM1::CleanApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanGhostFill(cln_w, src, 1);
  CleanArgs x = CleanMakeArgs();
  auto op = opac;

  par_for("cln_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real rr = x.feedback ? x.w(m,CLW_RR,k,j,i) : 0.0;
    x.w(m,dst,k,j,i) = x.w(m,src,k,j,i)*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr))
                       + CleanFluxDivergence(x, src, m, k, j, i, false);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPreconditionerSetup
//! \brief coefficients of the preconditioner: the x1-line tridiagonal part of A plus the
//! diagonal of its x2 and x3 parts (the off-diagonal Eddington terms are left to the
//! Krylov iteration)

void RadiationM1::CleanPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto op = opac;

  par_for("cln_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real al_p, ar_p, al_m, ar_m;
    CleanFaceWeights(x, 0, m, k, j, i+1, al_p, ar_p);
    CleanFaceWeights(x, 0, m, k, j, i, al_m, ar_m);
    const Real rr = x.feedback ? x.w(m,CLW_RR,k,j,i) : 0.0;
    Real a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr) + x.h*x.idx[0]*(al_p - ar_m);
    x.w(m,CLW_AM,k,j,i) = -x.h*x.idx[0]*al_m;
    x.w(m,CLW_AP,k,j,i) = x.h*x.idx[0]*ar_p;
    if (x.multi_d) {
      CleanFaceWeights(x, 1, m, k, j+1, i, al_p, ar_p);
      CleanFaceWeights(x, 1, m, k, j, i, al_m, ar_m);
      a0 += x.h*x.idx[1]*(al_p - ar_m);
    }
    if (x.three_d) {
      CleanFaceWeights(x, 2, m, k+1, j, i, al_p, ar_p);
      CleanFaceWeights(x, 2, m, k, j, i, al_m, ar_m);
      a0 += x.h*x.idx[2]*(al_p - ar_m);
    }
    x.w(m,CLW_A0,k,j,i) = a0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPrecondition
//! \brief w(dst) = M^{-1} w(src): tridiagonal (Thomas) solve along every x1 line

void RadiationM1::CleanPrecondition(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;

  par_for("cln_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = (i > is) ? w(m,CLW_AM,k,j,i) : 0.0;
      const Real ap = (i < ie) ? w(m,CLW_AP,k,j,i) : 0.0;
      const Real den = w(m,CLW_A0,k,j,i) - am*cprev;
      cprev = ap/den;
      dprev = (w(m,src,k,j,i) - am*dprev)/den;
      w(m,CLW_CP,k,j,i) = cprev;
      w(m,dst,k,j,i) = dprev;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,CLW_CP,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanBiCGStab
//! \brief right-preconditioned BiCGStab for A X = B (work vectors CLW_X, CLW_B), starting
//! from the X it is given.  Returns the number of iterations, negative if not converged.

int RadiationM1::CleanBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;

  // w(dst) = sa w(a) + sb w(b)
  auto axpby = [&](const int dst, const int a, const Real sa, const int b, const Real sb) {
    par_for("cln_axpby", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,dst,k,j,i) = sa*w(m,a,k,j,i) + sb*w(m,b,k,j,i);
    });
  };

  CleanApply(CLW_X, CLW_R);
  axpby(CLW_R, CLW_B, 1.0, CLW_R, -1.0);
  axpby(CLW_R0, CLW_R, 1.0, CLW_R, 0.0);
  axpby(CLW_P, CLW_R, 0.0, CLW_R, 0.0);
  axpby(CLW_V, CLW_R, 0.0, CLW_R, 0.0);
  const Real tol = cln_tol*fmax(std::sqrt(CleanDot(CLW_B, CLW_B)), 1.0e-300);
  if (std::sqrt(CleanDot(CLW_R, CLW_R)) <= tol) return 0;

  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  for (int it = 1; it <= cln_maxit; ++it) {
    const Real rho_new = CleanDot(CLW_R0, CLW_R);
    if (rho_new == 0.0) return -it;
    const Real beta = (rho_new/rho)*(alpha/omega);
    rho = rho_new;
    par_for("cln_pupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_P,k,j,i) = w(m,CLW_R,k,j,i)
                         + beta*(w(m,CLW_P,k,j,i) - omega*w(m,CLW_V,k,j,i));
    });
    CleanPrecondition(CLW_P, CLW_Y);
    CleanApply(CLW_Y, CLW_V);
    const Real r0v = CleanDot(CLW_R0, CLW_V);
    if (r0v == 0.0) return -it;
    alpha = rho/r0v;
    axpby(CLW_S, CLW_R, 1.0, CLW_V, -alpha);
    if (std::sqrt(CleanDot(CLW_S, CLW_S)) <= tol) {
      axpby(CLW_X, CLW_X, 1.0, CLW_Y, alpha);
      return it;
    }
    CleanPrecondition(CLW_S, CLW_Z);
    CleanApply(CLW_Z, CLW_Q);
    const Real qq = CleanDot(CLW_Q, CLW_Q);
    omega = (qq > 0.0) ? CleanDot(CLW_Q, CLW_S)/qq : 0.0;
    par_for("cln_xupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_X,k,j,i) += alpha*w(m,CLW_Y,k,j,i) + omega*w(m,CLW_Z,k,j,i);
      w(m,CLW_R,k,j,i) = w(m,CLW_S,k,j,i) - omega*w(m,CLW_Q,k,j,i);
    });
    if (std::sqrt(CleanDot(CLW_R, CLW_R)) <= tol) return it;
    if (omega == 0.0) return -it;
  }
  return -cln_maxit;
}

}  // namespace radm1
