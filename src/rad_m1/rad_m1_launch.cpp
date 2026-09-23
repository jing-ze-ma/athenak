//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_launch.cpp
//! \brief implicit M1: the BiCGStab of implicit_krylov_fuse = 3 with its scalars and its
//! convergence test kept on the DEVICE (tests_m1/runs_4k_launch).
//!
//! <rad_m1>/implicit_krylov_dev = K (default 0 = off, read only when named).  One rank,
//! every neighbour on the rank (implicit_halo_direct), the stored stencil operator.
//!  * alpha, beta, omega, rho and the recursive max|r| live in a device array `kd`.  The
//!    two reductions of an iteration go into a device View, and the reduction's `final`
//!    (run once, on the device, by the last block) updates the scalars and the flags, so
//!    no reduction returns to the host.
//!  * one iteration = 6 launches: [precond red (+ the x, r update of the previous
//!    iteration + the p update)] [precond black] [A y with the halo read straight from
//!    the neighbour block, + (rhat,v), max|r|] [precond red (+ s)] [precond black]
//!    [A z, + (t,s), (t,t), (rhat,t)].  Two's upd and two halo kernels are gone.
//!  * the host queues iterations in batches (the first one as long as the previous
//!    solve, then K) and reads the status once per batch.  Once the device has stopped
//!    (converged candidate, breakdown), every queued kernel returns at once.
//!  * the true-residual confirmation, restarts and the line-Jacobi fallback are the host
//!    code of ImplicitBiCGStabTwo.
//! Same recurrence as ImplicitBiCGStabTwo (the same test one half-iteration late, the
//! same recurrence for rho); round-off may differ (reduction and FMA order).

#include <algorithm>
#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

namespace {
// slots of kd
constexpr int KD_RHO = 0, KD_RHON = 1, KD_AL = 2, KD_OM = 3, KD_RV = 4, KD_TABS = 5,
              KD_BSC = 6, KD_TOL = 7, KD_EW = 8, KD_STOP = 9, KD_NIT = 10, KD_PENDU = 11,
              KD_PEND = 12, KD_AMX = 13, KD_N = 16;
static_assert(KD_N + 4 <= M1_KD_SIZE, "kd too small");

struct M1DVal {
  Real s0, s1, s2, mx;
};

//! the ghost value the direct halo would have put at (m,c,k,j,i): the active cell of the
//! owning neighbour block (hd_src), or the cell itself where there is none
struct M1HDLoad {
  DvceArray5D<Real> a;
  DvceArray2D<int> tab;
  int is, ie, js, je, ks, ke, nx1, nx2, nx3;
  bool md, td, x1p;   // x1p: no block has an x1 neighbour (x1 ghosts are local)
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int c, const int k, const int j, const int i) const {
    if (x1p && (i < is || i > ie)) return a(m,c,k,j,i);
    const int o1 = (i < is) ? -1 : ((i > ie) ? 1 : 0);
    const int o2 = md ? ((j < js) ? -1 : ((j > je) ? 1 : 0)) : 0;
    const int o3 = td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    if (o1 == 0 && o2 == 0 && o3 == 0) return a(m,c,k,j,i);
    const int src = tab(m, (o1+1) + 3*(o2+1) + 9*(o3+1));
    if (src < 0) return a(m,c,k,j,i);
    return a(src,c,k - o3*nx3,j - o2*nx2,i - o1*nx1);
  }
};

//! the stored-stencil row of ImplicitStencilOp (same terms, same order), x read by `ld`
template <class LD>
KOKKOS_INLINE_FUNCTION
Real M1StRowL(const DvceArray5D<Real> &st_, const LD &ld, const int m, const int k,
              const int j, const int i, const bool edg, const bool thrd,
              const bool vfold) {
  Real y = st_(m,0,k,j,i)*ld(k,j,i)
           + st_(m,1,k,j,i)*ld(k,j,i-1) + st_(m,2,k,j,i)*ld(k,j,i+1)
           + st_(m,3,k,j,i)*ld(k,j-1,i) + st_(m,4,k,j,i)*ld(k,j+1,i);
  if (edg) {
    y += st_(m,7,k,j,i)*ld(k,j-1,i-1) + st_(m,8,k,j,i)*ld(k,j-1,i+1)
         + st_(m,9,k,j,i)*ld(k,j+1,i-1) + st_(m,10,k,j,i)*ld(k,j+1,i+1);
  }
  if (thrd) {
    y += st_(m,5,k,j,i)*ld(k-1,j,i) + st_(m,6,k,j,i)*ld(k+1,j,i);
    if (edg) {
      y += st_(m,11,k,j,i)*ld(k-1,j,i-1) + st_(m,12,k,j,i)*ld(k-1,j,i+1)
           + st_(m,13,k,j,i)*ld(k+1,j,i-1) + st_(m,14,k,j,i)*ld(k+1,j,i+1)
           + st_(m,15,k,j,i)*ld(k-1,j-1,i) + st_(m,16,k,j,i)*ld(k-1,j+1,i)
           + st_(m,17,k,j,i)*ld(k+1,j-1,i) + st_(m,18,k,j,i)*ld(k+1,j+1,i);
    }
  }
  if (vfold) {
    y += st_(m,19,k,j,i)*ld(k,j,i-2) + st_(m,20,k,j,i)*ld(k,j,i+2)
         + st_(m,21,k,j,i)*ld(k,j-2,i) + st_(m,22,k,j,i)*ld(k,j+2,i);
    if (thrd) {
      y += st_(m,23,k,j,i)*ld(k-2,j,i) + st_(m,24,k,j,i)*ld(k+2,j,i);
    }
  }
  return y;
}

//! y = A x with the halo read in place, plus the reductions of ImplicitStencilOp red = 4
//! ((rhat,y), max|r|) or red = 3 ((y,s), (y,y), (rhat,y)); `final` advances the scalars
struct M1DevOp {
  using value_type = M1DVal;
  DvceArray5D<Real> iw, st;
  DvceArray1D<Real> kd;
  M1HDLoad hl;
  int is, js, ks, ie, je, ke, ni, nji, nkji, cx, cy, rm, reach;
  bool edg, thrd, vfold, allin;   // allin: the ghosts are filled (halo kernel mode)

  KOKKOS_INLINE_FUNCTION void init(value_type &v) const {
    v.s0 = 0.0; v.s1 = 0.0; v.s2 = 0.0; v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION void join(value_type &d, const value_type &s) const {
    d.s0 += s.s0;
    d.s1 += s.s1;
    d.s2 += s.s2;
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION void operator()(const int idx, value_type &v) const {
    if (kd(KD_STOP) != 0.0) return;
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real y;
    // with x1p every x1 ghost is read in place, so only j and k decide (whole x1 rows
    // take one branch: no divergence inside a wavefront)
    const bool inner = allin || (hl.x1p || ((i - reach >= is) && (i + reach <= ie)))
                       && (j - reach >= js)
                       && (j + reach <= je) && (!thrd || ((k - reach >= ks)
                                                          && (k + reach <= ke)));
    if (inner) {
      auto ld = [&](const int kk, const int jj, const int ii) -> Real {
        return iw(m,cx,kk,jj,ii);
      };
      y = M1StRowL(st, ld, m, k, j, i, edg, thrd, vfold);
    } else {
      auto ld = [&](const int kk, const int jj, const int ii) -> Real {
        return hl(m,cx,kk,jj,ii);
      };
      y = M1StRowL(st, ld, m, k, j, i, edg, thrd, vfold);
    }
    iw(m,cy,k,j,i) = y;
    if (rm == 4) {
      v.s0 += iw(m,M1_IW_KRH,k,j,i)*y;
      Real a = fabs(iw(m,M1_IW_KR,k,j,i));
      v.mx = (a > v.mx) ? a : v.mx;
    } else {
      v.s0 += y*iw(m,M1_IW_KS,k,j,i);
      v.s1 += y*y;
      v.s2 += iw(m,M1_IW_KRH,k,j,i)*y;
    }
  }
  KOKKOS_INLINE_FUNCTION void final(value_type &v) const {
    if (kd(KD_STOP) != 0.0) return;
    if (rm == 4) {
      // the update the preconditioner kernels just made is in x and r now
      kd(KD_PENDU) = 0.0;
      kd(KD_AMX) = v.mx;
      if (kd(KD_PEND) != 0.0) {
        kd(KD_PEND) = 0.0;
        const bool dn = (kd(KD_EW) != 0.0) ? (v.mx < kd(KD_TABS))
                                           : (v.mx/kd(KD_BSC) < kd(KD_TOL));
        if (dn) {   // this half-iteration is discarded; the host confirms on the TRUE r
          kd(KD_STOP) = 1.0;
          return;
        }
      }
      kd(KD_NIT) += 1.0;
      const Real rv = v.s0;
      kd(KD_RV) = rv;
      if (!(fabs(rv) > M1_BCG_EPS)) {
        kd(KD_STOP) = 2.0;
        return;
      }
      kd(KD_AL) = kd(KD_RHON)/rv;
    } else {
      const Real ts = v.s0, tt2 = v.s1, rt = v.s2;
      const Real omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
      const Real alpha = kd(KD_AL);
      const Real rhon = kd(KD_RHON);
      kd(KD_OM) = omega;
      kd(KD_RHO) = rhon;
      const Real rn = rhon - alpha*kd(KD_RV) - omega*rt;
      kd(KD_RHON) = rn;
      kd(KD_PENDU) = 1.0;
      kd(KD_PEND) = 1.0;
      if (!(fabs(omega) > M1_BCG_EPS) || !(fabs(rn) > M1_BCG_EPS)) {kd(KD_STOP) = 2.0;}
    }
  }
};

//! M1PCRX (rad_m1_implicit.cpp) with its scalars read from kd: up = 1 first applies the
//! pending x, r update of the previous iteration (when kd PENDU is set), then p = r +
//! beta (p - omega v); up = 2 makes s = r - alpha v; up = 0 reads cr.  Every team returns
//! at once when the device has stopped.
template <typename T>
void M1PCRXD(const DvceArray5D<Real> &iw_, const DvceArray1D<Real> &kd_,
             Kokkos::TeamPolicy<DevExeSpace> policy,
             const int is, const int ie, const int js, const int je, const int ks,
             const int ke, const int nkj, const int njl, const bool colr, const int cl_,
             const int cs_, const bool thrd, const bool cyclic, const int cr,
             const int cz, const int up) {
  const int nx = ie - is + 1;
  const int nv = cyclic ? 5 : 4;   // a, b, c, r (+ u) per buffer
  size_t scr_size = ScrArray1D<T>::shmem_size(2*nv*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  Kokkos::parallel_for("m1_impl_pcrxd",
                       policy.set_scratch_size(0, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    if (kd_(KD_STOP) != 0.0) return;   // team-uniform
    const int m = tm.league_rank()/nkj;
    const int kk = (tm.league_rank() - m*nkj)/njl;
    const int jj = (tm.league_rank() - m*nkj)%njl;
    const int k = kk + ks;
    const int j = colr ? (js + 2*jj + ((cl_ + kk) & 1)) : (jj + js);
    if (j > je) return;   // team-uniform
    ScrArray1D<T> sw(tm.team_scratch(0), 2*nv*nx);
    Real alpha = 0.0, beta = 0.0, gam = 1.0;
    if (cyclic) {
      alpha = iw_(m,M1_IW_TC,k,j,ie);
      beta = iw_(m,M1_IW_TA,k,j,is);
      gam = -iw_(m,M1_IW_TB,k,j,is);
    }
    const bool pu = (up == 1) && (kd_(KD_PENDU) != 0.0);
    const Real al = kd_(KD_AL), om = kd_(KD_OM);
    const Real bt = (up == 1) ? ((kd_(KD_RHON)/kd_(KD_RHO))*(al/om)) : 0.0;
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
        Real rc = iw_(m,M1_IW_KR,k,j,ii);
        if (pu) {
          iw_(m,M1_IW_KX,k,j,ii) += al*iw_(m,M1_IW_KY,k,j,ii) + om*iw_(m,M1_IW_KZ,k,j,ii);
          rc = iw_(m,M1_IW_KS,k,j,ii) - om*iw_(m,M1_IW_KTT,k,j,ii);
          iw_(m,M1_IW_KR,k,j,ii) = rc;
        }
        rr = rc + bt*(iw_(m,M1_IW_KP,k,j,ii) - om*iw_(m,M1_IW_KV,k,j,ii));
        iw_(m,M1_IW_KP,k,j,ii) = rr;
      } else if (up == 2) {
        rr = iw_(m,M1_IW_KR,k,j,ii) - al*iw_(m,M1_IW_KV,k,j,ii);
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

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::ImplicitKrylovDevOK
//! \brief whether implicit_krylov_dev can run this solve (else ImplicitBiCGStabTwo runs)

bool RadiationM1::ImplicitKrylovDevOK() {
  if (!((impl_kdev > 0) && (global_variable::nranks == 1) && halo_direct_on
        && impl_stencil && !(vimp_now && !impl_vfold))) {
    return false;
  }
  if (kdev_x1p < 0) {   // once: does any block have an x1 neighbour (periodic x1)?
    hd_src.sync_host();
    int p = 1;
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      for (int d = 0; d < 27; ++d) {
        if ((d % 3) != 1 && hd_src.h_view(m,d) >= 0) {p = 0;}
      }
    }
    kdev_x1p = p;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPrecondXD
//! \brief ImplicitPrecondX with the scalars on the device (M1PCRXD)

void RadiationM1::ImplicitPrecondXD(int zc, int upd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = ie - is + 1;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int rr = (upd == 1) ? M1_IW_KP : M1_IW_KS;
  auto sweep = [&](int rc, int up, int col, int sub) {
    const bool colr = (col >= 0);
    const int njl = colr ? (nj + 1)/2 : nj;
    const int nkj = nk*njl;
    const int cl_ = colr ? col : 0;
    const int cs_ = (colr && sub >= 0) ? sub : -1;
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
      M1PCRXD<float>(iw, kdv, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                     trans_x3, cyclic, rc, zc, up);
    } else {
      M1PCRXD<Real>(iw, kdv, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                    trans_x3, cyclic, rc, zc, up);
    }
  };
  if (impl_prec == 0) {
    sweep(-1, upd, -1, -1);
    return;
  }
  sweep(-1, upd, 0, -1);
  sweep(-1, upd, 1, zc);
  if (impl_prec == 1) {sweep(rr, 0, 0, zc);}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOpXD
//! \brief y = A x with the halo read in place and the reductions of red = 4 or 3 into a
//! device View; the reduction's final advances the scalars in kd

void RadiationM1::ImplicitOpXD(int xc, int yc, int red) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  M1DevOp f;
  f.iw = iw;
  f.st = ost;
  f.kd = kdv;
  f.hl.a = iw;
  f.hl.tab = hd_src.d_view;
  f.hl.is = is; f.hl.ie = ie; f.hl.js = js; f.hl.je = je; f.hl.ks = ks; f.hl.ke = ke;
  f.hl.nx1 = indcs.nx1; f.hl.nx2 = indcs.nx2; f.hl.nx3 = indcs.nx3;
  f.hl.md = (indcs.nx2 > 1);
  f.hl.td = (indcs.nx3 > 1);
  f.hl.x1p = kdev_x1p;
  f.is = is; f.js = js; f.ks = ks; f.ie = ie; f.je = je; f.ke = ke;
  f.ni = ie - is + 1;
  f.nji = (je - js + 1)*f.ni;
  f.nkji = (ke - ks + 1)*f.nji;
  f.cx = xc;
  f.cy = yc;
  f.rm = red;
  f.edg = st_edges;
  f.thrd = trans_x3;
  f.vfold = impl_vfold && vimp_now;
  f.reach = f.vfold ? 2 : 1;
  // implicit_krylov_dev_halo = 1: the ghosts by the direct-halo kernel, the operator
  // reads them in place (one launch more, no neighbour lookup in the operator)
  f.allin = (impl_kdev_halo == 1);
  if (f.allin) {ImplicitHaloDirect(1, xc);}
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*f.nkji);
  // the reduction's result goes into the last 4 slots of kd (never read back)
  Kokkos::View<M1DVal, DevMemSpace> res(reinterpret_cast<M1DVal *>(kdv.data() + KD_N));
  Kokkos::parallel_reduce("m1_impl_stord", pol, f, res);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabDev
//! \brief implicit_krylov_dev = K: ImplicitBiCGStabTwo with the scalars on the device
//! (see the file header).  Host syncs per solve: r0, one per batch, the true residual.

int RadiationM1::ImplicitBiCGStabDev(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto kd_ = kdv;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);

  // x0 = the Picard iterate; r0 = b - A x0, with max|r0| and (r0,r0): as Two
  par_for("m1_impl_bcgd_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  Real rs0 = 0.0, rmx = 0.0;
  auto rsf = KOKKOS_LAMBDA(const int idx, Real &l0, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real rr = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = rr;
    iw_(m,M1_IW_KRH,k,j,i) = rr;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
    l0 += rr*rr;
    Real a = fabs(rr);
    lmx = (a > lmx) ? a : lmx;
  };
  Kokkos::parallel_reduce("m1_impl_bcgd_r0", pol, rsf, rs0, Kokkos::Max<Real>(rmx));
  bcg_nred += 1.0;
  Real rnorm = rmx;
  Real rhon = rs0;
  bcg_r0rel = rnorm/bscale;
  const bool ew = (impl_ew_max > 0.0) && !ew_tight;
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
  auto true_ok = [&]() -> bool {
    ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
    Real tmx = 0.0;
    Kokkos::parallel_reduce("m1_impl_bcgd_true", pol,
    KOKKOS_LAMBDA(const int idx, Real &lmx) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/ni;
      int i = r - j*ni;
      k += ks; j += js; i += is;
      Real rr = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = rr;
      Real a = fabs(rr);
      lmx = (a > lmx) ? a : lmx;
    }, Kokkos::Max<Real>(tmx));
    bcg_nred += 1.0;
    return lin_done(tmx);
  };
  // the scalars of a (re)start, set on the device by one tiny kernel
  auto kinit = [&](const Real rn, const Real nitv) {
    const Real ta = tabs, bs = bscale, to = tol, ewv = ew ? 1.0 : 0.0;
    Kokkos::parallel_for("m1_impl_bcgd_init",
                         Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, 1),
    KOKKOS_LAMBDA(const int) {
      kd_(KD_RHO) = 1.0; kd_(KD_RHON) = rn; kd_(KD_AL) = 1.0; kd_(KD_OM) = 1.0;
      kd_(KD_RV) = 0.0; kd_(KD_TABS) = ta; kd_(KD_BSC) = bs; kd_(KD_TOL) = to;
      kd_(KD_EW) = ewv; kd_(KD_STOP) = 0.0; kd_(KD_NIT) = nitv; kd_(KD_PENDU) = 0.0;
      kd_(KD_PEND) = 0.0; kd_(KD_AMX) = 0.0;
    });
  };
  // x += alpha y + omega z, r = s - omega t: the update a stop left pending
  auto flush = [&]() {
    par_for("m1_impl_bcgd_flush", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real al = kd_(KD_AL), om = kd_(KD_OM);
      iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + om*iw_(m,M1_IW_KZ,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KS,k,j,i) - om*iw_(m,M1_IW_KTT,k,j,i);
    });
  };
  auto restart = [&]() {
    ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
    Real r2 = 0.0;
    Kokkos::parallel_reduce("m1_impl_bcgd_rs", pol,
    KOKKOS_LAMBDA(const int idx, Real &l0) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/ni;
      int i = r - j*ni;
      k += ks; j += js; i += is;
      Real rr = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = rr;
      iw_(m,M1_IW_KRH,k,j,i) = rr;
      iw_(m,M1_IW_KP,k,j,i) = 0.0;
      iw_(m,M1_IW_KV,k,j,i) = 0.0;
      l0 += rr*rr;
    }, r2);
    bcg_nred += 1.0;
    return r2;
  };

  int nit = 0;
  int nrestart = 0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  bool first = true;
  if (!done) {kinit(rhon, 0.0);}
  // a restart whose rho is already zero breaks down at once (Two: at the loop top)
  bool brk0 = !(fabs(rhon) > M1_BCG_EPS);
  while (!done) {
    int stop = 0;
    bool pendu = false;
    if (brk0) {
      stop = 2;
      brk0 = false;
    } else {
      if (nit >= impl_lin_maxit) {break;}
      // the first batch: as many as the last solve of this slot took, plus the one
      // half-iteration the test lags by; then K at a time
      int nq = first ? std::max(1, kdev_last[kdev_slot] + 1) : impl_kdev;
      first = false;
      nq = std::min(nq, impl_lin_maxit - nit);
      kdev_nq += static_cast<Real>(nq);
      for (int q = 0; q < nq; ++q) {
        ImplicitPrecondXD(M1_IW_KY, 1);
        ImplicitOpXD(M1_IW_KY, M1_IW_KV, 4);
        ImplicitPrecondXD(M1_IW_KZ, 2);
        ImplicitOpXD(M1_IW_KZ, M1_IW_KTT, 3);
      }
      Kokkos::deep_copy(DevExeSpace(), kdh, kdv);
      DevExeSpace().fence();
      kdev_nchk += 1.0;
      nit = static_cast<int>(kdh(KD_NIT));
      stop = static_cast<int>(kdh(KD_STOP));
      pendu = (kdh(KD_PENDU) != 0.0);
      bcg_nred += 2.0*nq;
    }
    if (stop == 1) {
      if (true_ok()) {
        done = true;
        break;
      }
      stop = 2;   // restart the recurrence from the true residual
    }
    if (stop == 2) {
      if (pendu) {flush();}
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        fell_back = true;
        break;
      }
      rhon = restart();
      kinit(rhon, static_cast<Real>(nit));
      brk0 = !(fabs(rhon) > M1_BCG_EPS);
      first = true;
      continue;
    }
    // stop == 0: the cap reached with an update pending (Two: test it once)
    if (nit >= impl_lin_maxit) {
      if (pendu) {flush();}
      if (true_ok()) {done = true;}
      break;
    }
  }
  kdev_last[kdev_slot] = nit;
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

} // namespace radm1
