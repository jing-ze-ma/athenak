#ifndef RAD_VET_RAD_VET_MOMENTUM_HPP_
#define RAD_VET_RAD_VET_MOMENTUM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_momentum.hpp
//! \brief device functions of the momentum and L_z exchange and of the gas state.

#include <cstdio>
#include <string>
#include "rad_vet/rad_vet_rows.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsArm
//! \brief cubed sphere: L_z arm lambda = (M x n)_z/A of lateral face (kf,jf,iff), M =
//! int x dA (CLG_M): momentum p n carries L_z lambda p

KOKKOS_INLINE_FUNCTION
Real CleanCsArm(const CleanArgs &x, const int q, const int m, const int kf, const int jf,
                const int iff) {
  return (x.g.At<CLG_M>(m,3*(q-1),kf,jf,iff)*x.g.At<CLG_N>(m,3*q+1,kf,jf,iff) -
          x.g.At<CLG_M>(m,3*(q-1)+1,kf,jf,iff)*x.g.At<CLG_N>(m,3*q,kf,jf,iff))/
         x.area[q](m,kf,jf,iff);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagCsVector
//! \brief cubed sphere: the cell vector u from its six face projections p (c = w p,
//! face 2 q + s) by weighted least squares, (sum w n n^T) u = sum c n; exact for uniform
//! u.  clean_stag_csvec_radial: lateral normals less their r_hat part (radial u exact)

KOKKOS_INLINE_FUNCTION
void CleanStagCsVector(const CleanArgs &x, const int m, const int k, const int j,
                       const int i, const Real c[6], const Real w[6], Real u[3]) {
  Real a[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Real b[3] = {0.0, 0.0, 0.0};
  for (int q = 0; q < 3; ++q) {
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1), iff = i + s*(q == 0);
      Real n[3], nr = 0.0;
      for (int r = 0; r < 3; ++r) {
        n[r] = x.g.At<CLG_N>(m,3*q+r,kf,jf,iff);
        nr += n[r]*x.g(m,CLG_R+r,k,j,i);
      }
      for (int r = 0; r < 3; ++r)
        n[r] -= (q > 0) ? nr*x.g(m,CLG_R+r,k,j,i) : 0.0;
      for (int r = 0; r < 3; ++r) {
        b[r] += c[2*q+s]*n[r];
        for (int t = 0; t < 3; ++t) a[r][t] += w[2*q+s]*n[r]*n[t];
      }
    }
  }
  // (Cramer, written out: no runtime-indexed local)
  Real ca[3][3];
  ca[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ca[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  ca[0][2] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  ca[1][0] = a[2][1]*a[0][2] - a[2][2]*a[0][1];
  ca[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
  ca[1][2] = a[2][0]*a[0][1] - a[2][1]*a[0][0];
  ca[2][0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ca[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  ca[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];
  const Real det = a[0][0]*ca[0][0] + a[0][1]*ca[0][1] + a[0][2]*ca[0][2];
  for (int t = 0; t < 3; ++t) u[t] = (ca[0][t]*b[0] + ca[1][t]*b[1] + ca[2][t]*b[2])/det;
}

//! \brief cell vector u from its six face values c: the mean of each family's two
//! faces; on the cubed sphere CleanStagCsVector with unit weights

KOKKOS_INLINE_FUNCTION
void CleanStagCellVector(const CleanArgs &x, const int m, const int k, const int j,
                         const int i, const Real c[6], Real u[3]) {
  if (x.cub) {
    const Real w[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    CleanStagCsVector(x, m, k, j, i, c, w, u);
    return;
  }
  for (int q = 0; q < 3; ++q) u[q] = 0.5*(c[2*q] + c[2*q+1]);
}

//----------------------------------------------------------------------------------------
//! \brief the face families with an L_z arm, (z x r).e_d != 0: spherical polar x3, the
//! cubed sphere x2, x3 (lateral), Cartesian x1, x2 (L_z about z); the others, without
//! one, are those a cap or cut may take momentum from (clean_stag_fcone_lz)
KOKKOS_INLINE_FUNCTION
bool CleanLzFamily(const bool sph, const bool cub, const int d) {
  return sph ? (d == 2) : (cub ? (d > 0) : (d < 2));
}

//! \brief the arm (z x r).e_dir of momentum along e_dir in cell (k,j,i), or at its face
//! normal to dir: spherical polar R along x3, Cartesian -y along x1 and x along x2 (the
//! cell centres of clean_lz_spin), else 0 (the cubed sphere's faces: CleanCsArm)
KOKKOS_INLINE_FUNCTION
Real CleanDirLzArm(const CleanArgs &x, const int dir, const int m, const int k,
                   const int j, const int i) {
  if (x.sph) return (dir == 2) ? x.x1v(m,i)*x.tg(m,2,j) : 0.0;
  if (x.cub || x.x2v.extent_int(0) == 0) return 0.0;
  return (dir == 0) ? -x.x2v(m,j) : ((dir == 1) ? x.x1v(m,i) : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanMomentumChange
//! \brief momentum the radiation gives cell (k,j,i) along dir over the stage, in
//! conservative form: per face -[(F - G)/c^2 + h div(D E)], shared half and half
//! (clean_stag_split: by chi).  cf, wf: each face's term and weight; tq: its whole term

template <bool kFs = false>
KOKKOS_INLINE_FUNCTION
Real CleanMomentumChange(const CleanArgs &x, const int e, const int dir, const int m,
                         const int k, const int j, const int i, Real *mag = nullptr,
                         Real *cf = nullptr, Real *wf = nullptr, Real *tq = nullptr) {
  // (kFs: x.fs and x.fa known set: stored face fluxes and terms read)
  const bool fs = kFs || x.fs, fa = kFs || x.fa;
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  Real dm = 0.0;
  for (int side = 0; side < 2; ++side) {
    const int kf = k + side*dk, jf = j + side*dj, iff = i + side*di;
    const Real f = fs ? x.fc[dir](m,CLF_N,kf,jf,iff) :
                   CleanFaceFluxCall(x, e, dir, m, kf, jf, iff, true);
    const Real g = x.fc[dir](m,CLF_G,kf,jf,iff);
    // (clean_stag_split: share chi_c/(chi_L + chi_R))
    const int o = 2*side - 1;
    const Real xc = x.op(m,VET_OP_T,k,j,i), xs = xc +
        x.op(m,VET_OP_T,k+o*dk,j+o*dj,i+o*di);
    Real wt = (x.hi && xs > 0.0) ? xc/xs : 0.5;
    if (x.hi && (x.sph || x.cub)) {   // (each half cell its absorption of V_f F)
      const int ko = k + o*dk, jo = j + o*dj, io = i + o*di;
      const Real vc = x.cub ? x.g(m,CLG_F+3,k,j,i) : x.vol(m,k,j,i);
      const Real vo = x.cub ? x.g(m,CLG_F+3,ko,jo,io) : x.vol(m,ko,jo,io);
      const Real xv = xc*vc + x.op(m,VET_OP_T,k+o*dk,j+o*dj,i+o*di)*vo;
      wt = (xv > 0.0) ? 0.5*xc*(vc + vo)/xv : 0.5;
    }
    if (x.hi && CleanBoundaryFace(x, dir, m, kf, jf, iff)) {
      // (boundary flux: the inner half cell's absorption of it, chi dx/2; gflux on an
      // open face: 1 - exp(-chi dx/2))
      const Real dc = 1.0/CleanDivW(x, dir, m, k, j, i, side), q = 0.5*xc*dc;
      const bool gq = x.gf && CleanOpen(x, m, 2*dir + side);
      const Real db = gq ? -expm1(-q)*x.h*f/(x.c*dc) : 0.5*x.h*xc*f/x.c;
      dm += db;
      if (cf) {
        cf[side] = db;
        wf[side] = 0.5;
      }
      if (tq) tq[side] = 2.0*db;
      continue;
    }
    const Real tu = fa ? x.fc[dir](m,CLF_PT,kf,jf,iff) : ((f - g)/(x.c*x.c) +
                    x.h*CleanFaceStressCall(x, e, dir, m, kf, jf, iff, true) -
                    (fs ? x.h*x.fc[dir](m,CLF_N+1,kf,jf,iff) : 0.0));
    if (tq) tq[side] = tu;
    Real t = wt*tu;
    // (mcl, the Newton passes: the conservative term equals the face's absorption
    // h chi_f F/c only where the iterate's F satisfies its equation; capped there at the
    // force of a realisable flux, |F| <= c E, as the production deposit, so that a lagged
    // F does not kick floor gas to ~c; the write-back keeps the conservative term, with
    // clean_wb_mcap the capped one, whose work the E row paid)
    if (x.mcl) {
      const Real cap = wt*(1.0/x.fc[dir](m,CLF_TH,kf,jf,iff) - 1.0)/x.c*
                       fmax(fmax(x.w(m,e,kf,jf,iff), x.w(m,e,kf-dk,jf-dj,iff-di)), 0.0);
      const Real tc = fmin(fmax(t, -cap), cap);
      if (tc != t) {   // (counted; a write-back's also the gas momentum withheld)
        Kokkos::atomic_add(&x.mclc(0), 1.0);
        if (x.mclc.extent_int(0) > 1) {
          Kokkos::atomic_add(&x.mclc(1+dir), (tc - t)*CleanVol(x, m, k, j, i));
        }
      }
      t = tc;
    }
    dm -= t;
    if (cf) {
      cf[side] = -t;
      wf[side] = wt;
    }
    // (mag: the size of its terms, for the round-off level of dm)
    auto de = [&](const int kk, const int jj, const int ii) {
      return fabs(x.cs(m,CLN_D11+dir,kk,jj,ii)*x.w(m,e,kk,jj,ii));
    };
    const Real id = (x.sph || x.cub) ? 1.0/x.dxf[dir](m,kf,jf,iff) : x.idx[dir];
    if (mag) *mag += 0.5*((fabs(f) + fabs(g))/(x.c*x.c) + x.h*id*
                          (de(kf,jf,iff) + de(kf-dk,jf-dj,iff-di)));
  }
  return dm;
}

// dissipation weight candidates (<rad_vet>/clean_alpha)
constexpr int CLA_TAU = 0;          // A(tau)^2, A = (1 - exp(-tau))/tau
constexpr int CLA_TAU_SENSOR = 1;   // max(A(tau)^2 psi, s_ray)
constexpr int CLA_TAU_RAY = 2;      // max(A(tau)^2, s_ray)
constexpr int CLA_ONE = 3;          // 1
constexpr int CLA_ZERO = 4;         // 0

//----------------------------------------------------------------------------------------
//! \fn Real CleanDissipationWeight
//! \brief face dissipation weight al from tau = chi_f dx_d, sensor psi and the rays'
//! streaming weight sray; rc = (1 - al)/tau (analytic to tau = 0 for clean_alpha = tau).
//! A^2 makes the dissipation O(1/tau^2) against the O(1/tau) diffusion flux

KOKKOS_INLINE_FUNCTION
Real CleanDissipationWeight(const int mode, const Real tau, const Real psi,
                            const Real sray, Real &rc, const int p = 2) {
  // A = 1 - b tau
  const Real b = (tau > 1.0e-4) ? (1.0 + expm1(-tau)/tau)/tau
                                : 0.5 - tau*(1.0 - 0.25*tau)/6.0;
  const Real a = (tau > 1.0e-4) ? -expm1(-tau)/tau : 1.0 - b*tau;
  // A^p, p = <rad_vet>/clean_dissip_power
  const Real ap = (p == 3) ? a*a*a : a*a;
  Real al = 0.0;
  switch (mode) {
    case CLA_TAU: al = ap; break;
    case CLA_TAU_SENSOR: al = fmax(ap*psi, sray); break;
    case CLA_TAU_RAY: al = fmax(ap, sray); break;
    case CLA_ONE: al = 1.0; break;
    default: break;
  }
  rc = (mode == CLA_TAU) ? b*(1.0 + a + ((p == 3) ? a*a : 0.0))
                         : (1.0 - al)/fmax(tau, 2.0);
  return al;
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanDiag
//! \brief component D_ab of the Eddington tensor of cell (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanDiag(const CleanArgs &x, const int a, const int b, const int m, const int k,
               const int j, const int i) {
  return x.cs(m,CLN_D11+((a == b) ? a : 2 + a + b),k,j,i);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanNf
//! \brief component b of the unit normal of face (kf,jf,iff) (e_dir off the cubed sphere)

KOKKOS_INLINE_FUNCTION
Real CleanNf(const CleanArgs &x, const int b, const int dir, const int m, const int kf,
             const int jf, const int iff) {
  return x.cub ? x.g.At<CLG_N>(m,3*dir+b,kf,jf,iff) : ((b == dir) ? 1.0 : 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsFrame
//! \brief cubed sphere: frame of y, rows r_hat, theta_hat = phi_hat x r_hat, phi_hat =
//! z x r_hat/|.| (y_cart = Q^T y).  Returns the L_z arm R = x1v |z x r_hat| (> 0: odd
//! panel resolutions are refused)

KOKKOS_INLINE_FUNCTION
Real CleanCsFrame(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  Real q[3][3]) {
  for (int c = 0; c < 3; ++c) q[0][c] = x.g.At<CLG_R>(m,c,k,j,i);
  const Real rho = std::sqrt(q[0][0]*q[0][0] + q[0][1]*q[0][1]);
  q[2][0] = -q[0][1]/rho;
  q[2][1] = q[0][0]/rho;
  q[2][2] = 0.0;
  q[1][0] = q[2][1]*q[0][2];
  q[1][1] = -q[2][0]*q[0][2];
  q[1][2] = q[2][0]*q[0][1] - q[2][1]*q[0][0];
  return x.x1v(m,i)*rho;
}

//! \brief the arm R of L_z = R p_phi of a cell vector p in the frame of y: spherical
//! polar x1v sin(theta), the cubed sphere CleanCsFrame's; Cartesian |(x, y)| with the
//! cell centres of clean_lz_spin, else 0
KOKKOS_INLINE_FUNCTION
Real CleanLzArm(const CleanArgs &x, const int m, const int k, const int j, const int i) {
  Real q[3][3];
  if (!x.sph && !x.cub) {
    return (x.x2v.extent_int(0) > 0) ? sqrt(SQR(x.x1v(m,i)) + SQR(x.x2v(m,j))) : 0.0;
  }
  return x.cub ? CleanCsFrame(x, m, k, j, i, q) : x.x1v(m,i)*x.tg(m,2,j);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanCsFamily
//! \brief cubed sphere: face family normals of the half-range closure: r_hat and the
//! mean normals of the x2 and x3 face pairs

KOKKOS_INLINE_FUNCTION
void CleanCsFamily(const CleanArgs &x, const int m, const int k, const int j, const int i,
                   Real nd[3][3]) {
  for (int d = 0; d < 3; ++d) {
    const int dk = (d == 2), dj = (d == 1);
    Real s = 0.0;
    for (int c = 0; c < 3; ++c) {
      nd[d][c] = (d == 0) ? x.g.At<CLG_R>(m,c,k,j,i) :
                 x.g.At<CLG_N>(m,3*d+c,k,j,i) + x.g.At<CLG_N>(m,3*d+c,k+dk,j+dj,i);
      s += nd[d][c]*nd[d][c];
    }
    for (int c = 0; c < 3; ++c) nd[d][c] /= std::sqrt(s);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanRayFlux
//! \brief rays' flux factor h = H/J from s_b = h+_b + h-_b: s; cubed sphere: solves
//! n_b . h = s_b (CleanCsFamily normals)

KOKKOS_INLINE_FUNCTION
void CleanRayFlux(const CleanArgs &x, const int m, const int k, const int j, const int i,
                  const Real s[3], Real h[3]) {
  for (int b = 0; b < 3; ++b) h[b] = s[b];
  if (!x.cub) return;
  Real nd[3][3], c[3][3];   // c_b = n_{b+1} x n_{b+2}
  CleanCsFamily(x, m, k, j, i, nd);
  for (int b = 0; b < 3; ++b) {
    const Real *p = nd[(b+1)%3], *q = nd[(b+2)%3];
    for (int e = 0; e < 3; ++e) c[b][e] = p[(e+1)%3]*q[(e+2)%3] - p[(e+2)%3]*q[(e+1)%3];
  }
  const Real det = nd[0][0]*c[0][0] + nd[0][1]*c[0][1] + nd[0][2]*c[0][2];
  for (int e = 0; e < 3; ++e) h[e] = (s[0]*c[0][e] + s[1]*c[1][e] + s[2]*c[2][e])/det;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanGasVector, CleanGasMomentum
//! \brief gas momentum in the frame of y.  Cubed sphere: the hydro stores m_r = p.r_hat
//! and covariant m_a = p.e_a on non-orthogonal panel tangents: p_t = g^ab m_b e_a

KOKKOS_INLINE_FUNCTION
void CleanGasVector(const CleanArgs &x, const int m, const int k, const int j,
                    const int i, const Real u[3], Real p[3]) {
  for (int n = 0; n < 3; ++n) p[n] = u[n];
  if (!x.cub) return;
  Real q[3][3], pt[3], c = 0.0;
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) c += x.g.At<CLG_E>(m,b,k,j,i)*x.g.At<CLG_E>(m,3+b,k,j,i);
  const Real u1 = (p[1] - c*p[2])/(1.0 - c*c), u2 = (p[2] - c*p[1])/(1.0 - c*c);
  for (int b = 0; b < 3; ++b) {
    pt[b] = u1*x.g.At<CLG_E>(m,b,k,j,i) + u2*x.g.At<CLG_E>(m,3+b,k,j,i);
  }
  for (int a = 1; a < 3; ++a) p[a] = q[a][0]*pt[0] + q[a][1]*pt[1] + q[a][2]*pt[2];
}

//! (CleanGasMomentum: CleanGasVector of the hydro's momentum)
KOKKOS_INLINE_FUNCTION
void CleanGasMomentum(const CleanArgs &x, const DvceArray5D<Real> &uh, const int m,
                      const int k, const int j, const int i, Real p[3]) {
  const Real u[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
  CleanGasVector(x, m, k, j, i, u, p);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanAdvect
//! \brief a = v + D v of cell (k,j,i), v in the frame of y

KOKKOS_INLINE_FUNCTION
void CleanAdvect(const CleanArgs &x, const int m, const int k, const int j, const int i,
                 const Real v[3], Real a[3]) {
  if (!x.cub) {
    for (int b = 0; b < 3; ++b) {
      a[b] = v[b];
      for (int c = 0; c < 3; ++c) a[b] += CleanDiag(x, b, c, m, k, j, i)*v[c];
    }
    return;
  }
  Real q[3][3], vc[3], ac[3];
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) vc[b] = v[0]*q[0][b] + v[1]*q[1][b] + v[2]*q[2][b];
  for (int b = 0; b < 3; ++b) {
    ac[b] = vc[b];
    for (int c = 0; c < 3; ++c) ac[b] += CleanDiag(x, b, c, m, k, j, i)*vc[c];
  }
  for (int b = 0; b < 3; ++b) a[b] = q[b][0]*ac[0] + q[b][1]*ac[1] + q[b][2]*ac[2];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanGasChange
//! \brief hydro momentum change dh for a change dm in the frame of y (inverse of
//! CleanGasMomentum; cubed sphere: dh_r = dm_r, covariant dh_a = dm_t.e_a)

KOKKOS_INLINE_FUNCTION
void CleanGasChange(const CleanArgs &x, const int m, const int k, const int j,
                    const int i, const Real dm[3], Real dh[3]) {
  for (int n = 0; n < 3; ++n) dh[n] = dm[n];
  if (!x.cub) return;
  Real q[3][3], dt[3];
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) dt[b] = dm[1]*q[1][b] + dm[2]*q[2][b];
  for (int a = 1; a < 3; ++a) {
    dh[a] = dt[0]*x.g.At<CLG_E>(m,3*(a-1),k,j,i) +
            dt[1]*x.g.At<CLG_E>(m,3*(a-1)+1,k,j,i) +
            dt[2]*x.g.At<CLG_E>(m,3*(a-1)+2,k,j,i);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCellFaceDistance
//! \brief distance between the two cell centroids across face (k,j,i)

KOKKOS_INLINE_FUNCTION
Real CleanCellFaceDistance(const CleanArgs &x, const int dir, const int m, const int k,
                           const int j, const int i) {
  return (x.sph || x.cub) ? x.dxf[dir](m,k,j,i) : 1.0/x.idx[dir];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanHcTauW
//! \brief clean_face_hctau / clean_face_jse: weights of a cell value interpolated
//! linearly in optical depth to the face, wl = tau_R/(tau_L + tau_R) (half-cell depths),
//! times (r/r_f)^2 on radial faces; tg: tau_L, tau_R and the two (r/r_f)^2

KOKKOS_INLINE_FUNCTION
void CleanHcTauW(const CleanArgs &x, const int dir, const int m, const int k, const int j,
                 const int i, Real &wl, Real &wr, Real *tg = nullptr) {
  const int il = i - (dir == 0), jl = j - (dir == 1), kl = k - (dir == 2);
  const bool rad = (dir == 0 && (x.sph || x.cub));
  const Real dl = rad ? x.x1f(m,i) - x.x1v(m,il) : 1.0;
  const Real dr = rad ? x.x1v(m,i) - x.x1f(m,i) : 1.0;
  const Real tl = fmax(x.op(m,VET_OP_T,kl,jl,il), 0.0)*dl;
  const Real tr = fmax(x.op(m,VET_OP_T,k,j,i), 0.0)*dr;
  wl = (tl + tr > 0.0) ? tr/(tl + tr) : 0.5;
  wr = 1.0 - wl;
  const Real gl = rad ? SQR(x.x1v(m,il)/x.x1f(m,i)) : 1.0;
  const Real gr = rad ? SQR(x.x1v(m,i)/x.x1f(m,i)) : 1.0;
  wl *= gl;
  wr *= gr;
  const Real hd = rad ? 1.0 : 0.5/x.idx[dir];   // (Cartesian: chi -> the half-cell depth)
  if (tg) tg[0] = tl*hd, tg[1] = tr*hd, tg[2] = gl, tg[3] = gr;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanOverlapFor
//! \brief f over the inner cells (inset by ng) while a ghost exchange is in flight, then
//! fin(), then the shell; ov = false: fin() then all (same values).  fd >= 0: over faces
//! normal to fd.  post: shell first, post(), inner, fin() (the sending side)

// (fin a std::function: nvcc takes no host lambda type in an extended lambda's parent)
template <typename F>
void CleanOverlapFor(const std::string &nm, MeshBlockPack *pp, const bool ov, const F &f,
                     const std::function<void()> &fin, const int fd = -1,
                     const std::function<void()> &post = nullptr) {
  auto &id = pp->pmesh->mb_indcs;
  const int g = id.ng, nmb1 = pp->nmb_thispack - 1, is = id.is, js = id.js, ks = id.ks;
  const bool md = pp->pmesh->multi_d, td = pp->pmesh->three_d;
  const int nx[3] = {id.nx1, id.nx2, id.nx3};
  int n[3], lo[3], hi[3];
  for (int d = 0; d < 3; ++d) {
    const bool ac = (d == 0) || (d == 1 && md) || (d == 2 && td);
    n[d] = nx[d] + (d == fd);
    lo[d] = (fd < 0) ? (ac ? g : 0) : (d == fd);
    hi[d] = (fd < 0) ? (ac ? nx[d] - 1 - g : nx[d] - 1) : nx[d] - 1 - (d == fd);
  }
  const CleanShell sh(n[2], n[1], n[0], 0, lo[0], hi[0], lo[1], hi[1], lo[2], hi[2],
                      true, true);
  if (!ov || sh.Empty()) {
    if (!post) fin();
    par_for(nm, DevExeSpace(), 0, nmb1, ks, id.ke + (fd == 2), js, id.je + (fd == 1), is,
            id.ie + (fd == 0), f);
    if (post) post(), fin();
    return;
  }
  auto inner = [&]() {
    par_for(nm, DevExeSpace(), 0, nmb1, ks + sh.k0, ks + sh.k1, js + sh.j0, js + sh.j1,
            is + sh.i0, is + sh.i1, f);
  };
  if (!post) inner(), fin();
  par_for(nm + "_sh", DevExeSpace(), 0, nmb1, 0, sh.ncell-1,
  KOKKOS_LAMBDA(const int m, const int q) {
    int k, j, i;
    CleanShellCell(sh, q, k, j, i);
    f(m, ks + k, js + j, is + i);
  });
  if (post) post(), inner(), fin();
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanPredictE
//! \brief predicted E = e0 + de, positive: e0 exp(de/e0) where e0 + de <= 0

KOKKOS_INLINE_FUNCTION
Real CleanPredictE(const Real e0, const Real de) {
  return (e0 + de > 0.0) ? e0 + de : ((e0 > 0.0) ? e0*exp(de/e0) : 0.0);
}

//! \brief warn when the Zalesak rounds hit their cap unconverged (E not bounded)
inline void ZalesakCapWarn(const char *who, const int cycle, const int t2s) {
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: %s Zalesak rounds at the cap 4096, R not "
                "converged (E above the floor not guaranteed)\n", cycle, t2s, who);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagGasMomentum
//! \brief staggered: the gas momentum of cell (k,j,i) in the mesh's face-flux components
//! (Cartesian on the cubed sphere)

KOKKOS_INLINE_FUNCTION
static void CleanStagGasMomentum(const CleanArgs &x, const DvceArray5D<Real> &uh,
                                 const int m, const int k, const int j, const int i,
                                 Real p[3]) {
  CleanGasMomentum(x, uh, m, k, j, i, p);
  if (!x.cub) return;
  Real q[3][3];
  const Real y[3] = {p[0], p[1], p[2]};
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) p[b] = y[0]*q[0][b] + y[1]*q[1][b] + y[2]*q[2][b];
}

//! \brief staggered: v from CleanStagGasMomentum's components to the hydro's (changed on
//! the cubed sphere only)

KOKKOS_INLINE_FUNCTION
static void CleanStagHydro(const CleanArgs &x, const int m, const int k, const int j,
                           const int i, Real v[3]) {
  if (!x.cub) return;
  Real q[3][3];
  const Real y[3] = {v[0], v[1], v[2]};
  CleanCsFrame(x, m, k, j, i, q);
  for (int d = 0; d < 3; ++d) v[d] = q[d][0]*y[0] + q[d][1]*y[1] + q[d][2]*y[2];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagRefForce
//! \brief staggered: the reference force density withheld from the gas of cell (k,j,i),
//! in CleanStagGasMomentum's components.  clean: the stored force; wb_arad: d arad_ref
//! along x1 (radial on the cubed sphere: no torque); wb_arad_md: d times the stored
//! acceleration; else 0.  The gas then takes the residual of the radiation force.

KOKKOS_INLINE_FUNCTION
static void CleanStagRefForce(const CleanArgs &x, const int fref,
                              const DvceArray5D<Real> &fr, const DvceArray4D<Real> &aref,
                              const Real d, const int m, const int k, const int j,
                              const int i, Real rf[3]) {
  Real q[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  if (fref == VET_FREF_WB_ARAD && x.cub) CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) {
    rf[b] = (fref == VET_FREF_CLEAN) ? fr(m,b,k,j,i) :
            ((fref == VET_FREF_WB_ARAD) ? d*aref(m,k,j,i)*q[0][b] :
            ((fref == VET_FREF_WB_MD) ? d*fr(m,b,k,j,i) : 0.0));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagMomentum
//! \brief staggered: the radiation's momentum change of cell (k,j,i) over the stage, the
//! cell's shares of its faces' absorption (CleanMomentumChange; mag: the terms' size), in
//! CleanStagGasMomentum's components.
//! Cubed sphere: the face changes combine by CleanStagCsVector (weights 1/2).  Lateral
//! faces are no phi faces, so L_z does not telescope: dlz = -V R (dF_phi/c^2 + dp_phi) is
//! booked as the radiation's (CleanLzBook).  3-D spherical polar: dlz of the phi faces,
//! which telescopes to the L_z crossing the boundaries.

template <bool kFs = false>
KOKKOS_INLINE_FUNCTION
static void CleanStagMomentum(const CleanArgs &x, const int e, const int m, const int k,
                              const int j, const int i, Real dm[3], Real mag[3],
                              Real *dlz = nullptr) {
  Real cf[6], wf[6];
  for (int d = 0; d < 3; ++d) {
    mag[d] = 0.0;
    dm[d] = (d < x.nf) ? CleanMomentumChange<kFs>(x, e, d, m, k, j, i, &mag[d],
                                                  &cf[2*d], &wf[2*d]) : 0.0;
  }
  if (!x.cub) {
    if (dlz && x.sph && x.nf > 2) {
      Real df[2];
      for (int s = 0; s < 2; ++s) {
        const Real f = (kFs || x.fs) ? x.fc[2](m,CLF_N,k+s,j,i)
                                     : CleanFaceFlux(x, e, 2, m, k+s, j, i, true);
        df[s] = f - x.fc[2](m,CLF_G,k+s,j,i);
      }
      *dlz = (-0.5*(df[0] + df[1])/(x.c*x.c) - dm[2])*x.vol(m,k,j,i)*x.x1v(m,i)*
             x.tg(m,2,j);
    }
    return;
  }
  Real dv[3], q[3][3];
  const Real hw[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  CleanStagCsVector(x, m, k, j, i, cf, hw, dv);
  const Real rc = CleanCsFrame(x, m, k, j, i, q);
  const Real dr = q[0][0]*dv[0] + q[0][1]*dv[1] + q[0][2]*dv[2];
  const Real dt = q[1][0]*dv[0] + q[1][1]*dv[1] + q[1][2]*dv[2];
  const Real dp = q[2][0]*dv[0] + q[2][1]*dv[1] + q[2][2]*dv[2];
  if (dlz) {
    const Real one[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Real df[6], fd[3];
    for (int d = 0; d < 6; ++d) {
      const int n = d/2, kf = k + (d%2)*(n == 2), jf = j + (d%2)*(n == 1);
      const int iff = i + (d%2)*(n == 0);
      const Real f = (kFs || x.fs) ? x.fc[n](m,CLF_N,kf,jf,iff)
                                   : CleanFaceFlux(x, e, n, m, kf, jf, iff, true);
      df[d] = f - x.fc[n](m,CLF_G,kf,jf,iff);
    }
    CleanStagCsVector(x, m, k, j, i, df, one, fd);
    *dlz = (-(q[2][0]*fd[0] + q[2][1]*fd[1] + q[2][2]*fd[2])/(x.c*x.c) - dp)*
           x.vol(m,k,j,i)*rc;
  }
  for (int b = 0; b < 3; ++b) dm[b] = dr*q[0][b] + dt*q[1][b] + dp*q[2][b];
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_MOMENTUM_HPP_
