//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_solver.cpp
//! \brief the linear solve for E: face coefficients, the operator A E = E (1 + h c kE
//! (1-r)) + h div F_lin(E), x1-line / multigrid preconditioner, right-preconditioned
//! BiCGStab.
//!
//!  called from: rad_vet_newton.cpp, rad_vet_tasks.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_flux_seam.cpp, rad_vet_init.cpp

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_reduce.hpp"
#include "utils/scratch_view.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FaceCoefficients
//! \brief per face: blend weight w(C1), theta, old-vector flux G and h+_f = H+_f/J_L,
//! h-_f = H-_f/J_R (H+- minmod-reconstructed from the upwind side; cubed sphere: about
//! the face normal).  |h+-_f| <= 1 keeps the operator an M-matrix

void RadiationVet::FaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;
  const bool fsk = cln_fsolve;
  // (cubed sphere: seam ghost opacities from the neighbour's cells: single-valued faces)
  if (pmy_pack->pmesh->use_cubed_sphere) {
    GhostFill(opac, 0, opac.extent_int(1), false);
  }
  auto gn = cln_gn;
  auto fi = cln_fvi;
  auto dv = cln_fvd;
  auto fh = cln_fh;
  auto x1f = pmy_pack->pcoord->xx1f;
  auto mbc = pmy_pack->pmb->mb_bcs.d_view;
  const bool fsc = cln_cssc;   // (the sc rays: their folded sums, FaceHalf)
  const bool fho = cln_rfh;    // (vet_rotate_avg_ema_wb: the cells' own)
  const SolveArgs *px = DevArgs(MakeArgs(), 0);
  // (closure_form = eddington | sphere sweeps no rays: the cells' half ranges)
  const int nd = (pmy_pack->pmesh->use_cubed_sphere && !cln_edd && !cln_sphere)
                 ? cln_fvn : 0;

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = VET_T2_F1 + dir;

    SolveFor<0, ATHENA_LB_FACECOEF>("cln_facecoef", nmb1, ks, ke+dk, js, je+dj, is,
                                    ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = FaceChi(op, x1f, mbc, dir, m, k, j, i, is, ie);
      // (g = max(|H|/J, eps C1/1.5) of the cell)
      auto gw = [&](const int kk, const int jj, const int ii) {
        Real h2 = 0.0;
        for (int a = 0; a < 3; ++a) {
          h2 += SQR(cs(m,CLN_HP+a,kk,jj,ii) + cs(m,CLN_HM+a,kk,jj,ii));
        }
        const Real ct = op(m,VET_OP_T,kk,jj,ii);
        Real ep = (ct > 0.0) ? fmin(op(m,VET_OP_P,kk,jj,ii)/ct, 1.0) : 1.0;
        if (ep < 1.0) {   // (tau_d = chi V/A_f, the faces' mean)
          Real td = 2.0*ct/(DivW(*px, dir, m, kk, jj, ii, 0) +
                            DivW(*px, dir, m, kk, jj, ii, 1));
          if (px->cub && dir > 0) {   // (the cell's own family: one w at a seam)
            Real s;
            const int q = CsSlot(*px, dir, m, k, j, i, kk, jj, ii, false, s)/2;
            td = ct*px->g(m,CLG_F+3,kk,jj,ii)/px->g.At<CLG_A>(m,q,kk,jj,ii);
          }
          ep = 1.0 - (1.0 - ep)*(1.0 - exp(-td));
        }
        return fmax(sqrt(h2), ep*cs(m,CLN_C1,kk,jj,ii)/1.5);
      };
      const Real c1 = fmax(gw(k, j, i), gw(kl, jl, il));
      fc(m,CLF_W,k,j,i) = SmoothStep((c1 - c1lo)/(c1hi - c1lo));
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_G,k,j,i) = f0(m,k,j,i) + (stage ? inc(m,tf,k,j,i) : 0.0);
      if (fsk) fc(m,CLF_N+2,k,j,i) = fc(m,CLF_G,k,j,i);   // (fsolve: start of its F)

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
      // (scale h+ by H+_f/H+_L, H+_f = wl H+_L + wr H+_R (HcTauW),
      // factor in [1/2, 2], blended by eps = min(kappa_P/chi, 1) of the face (h- alike);
      // radial faces on curved meshes, all on Cartesian)
      auto jse = [&]() {
        if (!(dir == 0 || !(px->sph || px->cub))) return;
        auto ep = [&](const int kk, const int jj, const int ii) {
          const Real ct = op(m,VET_OP_T,kk,jj,ii);
          return (ct > 0.0) ? fmin(op(m,VET_OP_P,kk,jj,ii)/ct, 1.0) : 1.0;
        };
        const Real wt = fmin(ep(kl,jl,il), ep(k,j,i));
        Real wl, wr;
        HcTauW(*px, dir, m, k, j, i, wl, wr);
        const Real ppl = hp(kl,jl,il), ppr = fmax(hp(k,j,i), 0.0);
        const Real mml = fmin(hm(kl,jl,il), 0.0), mmr = hm(k,j,i);
        auto &fp = fc(m,CLF_HP,k,j,i), &fm = fc(m,CLF_HM,k,j,i);
        if (ppl > 0.0 && fp > 0.0) {
          const Real ph = fmin(fmax(wl + wr*ppr/ppl, 0.5), 2.0);
          fp += wt*(fmin(fp*ph, 1.0) - fp);
        }
        if (mmr < 0.0 && fm < 0.0) {
          const Real ph = fmin(fmax(wr + wl*mml/mmr, 0.5), 2.0);
          fm += wt*(fmax(fm*ph, -1.0) - fm);
        }
      };
      if (nd == 0) {
        jse();
        return;
      }
      // (cubed sphere: about the face normal from the FV intensities of L and R, so a
      // seam face takes the same h+- on both sides)
      // (J_L, J_R, H+_L, H-_R as scalars: nvcc -O3 miscompiles a runtime-indexed local)
      if (fho) {   // (each cell's value out through this face, in its own axes)
        Real s;
        const int nl = CsSlot(*px, dir, m, k, j, i, kl, jl, il, false, s);
        const int nr = CsSlot(*px, dir, m, k, j, i, k, j, i, false, s);
        fc(m,CLF_HP,k,j,i) = (nl & 1) ? cs(m,CLN_HO+nl/2,kl,jl,il)
                                      : -cs(m,CLN_HO+3+nl/2,kl,jl,il);
        fc(m,CLF_HM,k,j,i) = (nr & 1) ? -cs(m,CLN_HO+nr/2,k,j,i)
                                      : cs(m,CLN_HO+3+nr/2,k,j,i);
        jse();
        return;
      }
      Real jl0 = 0.0, jr0 = 0.0, hp0 = 0.0, hm0 = 0.0;
      if (fsc) {
        jl0 = fh(m,4*dir,k,j,i);
        hp0 = fh(m,4*dir+2,k,j,i);
        hm0 = fh(m,4*dir+3,k,j,i);
        jr0 = fh(m,4*dir+1,k,j,i);
      }
      for (int d = 0; d < nd && !fsc; ++d) {
        Real mu = 0.0;
        for (int b = 0; b < 3; ++b) mu += dv(d,b)*gn.At<CLG_N>(m,3*dir+b,k,j,i);
        const Real wl = dv(d,3)*fi(m,d,kl,jl,il), wr = dv(d,3)*fi(m,d,k,j,i);
        jl0 += wl;
        jr0 += wr;
        if (mu > 0.0) hp0 += wl*mu;
        if (!(mu > 0.0)) hm0 += wr*mu;
      }
      fc(m,CLF_HP,k,j,i) = (jl0 > 0.0) ? hp0/jl0 : 0.25;
      fc(m,CLF_HM,k,j,i) = (jr0 > 0.0) ? hm0/jr0 : -0.25;
      jse();
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ApplyOperator
//! \brief w(dst) = A w(src) over the active cells; src's exchange overlaps the inner
//! cells (cubed sphere: only with the xdonor, else the lateral sums need ghosts first)

void RadiationVet::ApplyOperator(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const bool ov = !pmy_pack->pmesh->use_cubed_sphere || cln_xd_on;
  auto wa = (src < cln_w.ns) ? cln_w.s : cln_w.o;   // (src's array, CLW_NS)
  const int ws = (src < cln_w.ns) ? src : src - cln_w.ns;
  // (the 19-point stencil, the cubed sphere's cross terms and the lateral sums read the
  // first ghost layer only: an exchange of depth 1)
  GhostFill(wa, ws, 1, true, ov ? 1 : 0, 1);
  const std::function<void()> fin = [&]() {if (ov) GhostFill(wa, ws, 1, true, 2);};
  // (cubed sphere: the exchange of src's lateral sums Y overlaps the inner cross terms)
  const bool yov = pmy_pack->pmesh->use_cubed_sphere && !cln_xd_on;
  if (!cln_xd_on) CsLateralStress(src, true, yov ? 1 : 0);   // (xdonor: none)
  SolveArgs x = MakeArgs();
  x.fwc = true;
  // (the operator from its rows, PreconditionerSetup)
  auto w = cln_w;
  auto st = cln_st;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  // (cubed sphere: lateral stress from the sums Y; none with the xdonor)
  const bool cub = x.cub, cxf = cub && !(x.xq >= 0 && !x.hi);
  const SolveArgs *px = cxf ? DevArgs(x, 0) : nullptr;
  // (cubed sphere: each face's cross term once, not per adjacent cell)
  DvceArray5D<Real> xf;
  if (cxf) {
    const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
    xf = ScratchView<DvceArray5D<Real>>("cln_xf", pmy_pack->nmb_thispack, 3,
                                        w.s.extent_int(2), w.s.extent_int(3),
                                        w.s.extent_int(4));
    // (pass 0: faces between two active cells; pass 1, after the exchange: end planes)
    for (int ps = 0; ps < 2; ++ps) {
      if (ps == 1 && yov) GhostFill(cln_w, CLW_L, 9, 2, 1);
      const bool xon0 = x.multi_d && !(x.xq >= 0 && !x.hi);
      if (ps == 1 || !xon0) {
        // all directions' faces in one launch: pass 1 the end planes, pass 0 (no inner
        // kernel) the inner faces
        const int ne[3] = {ie - is + 1, je - js + 1, ke - ks + 1};
        Kokkos::Array<int, 7> c0;   // (segment s: faces c0[s] .. c0[s+1]-1 per block)
        Kokkos::Array<int, 6> sd, sk0, sj0, si0, snj, sni;
        int ns = 0;
        c0[0] = 0;
        for (int d = 0; d < x.nf; ++d) {
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          for (int pl = 0; pl < (ps ? 2 : 1); ++pl) {
            const int o = pl*ne[d];
            const int k0 = ps ? ks + o*dk : ks + dk, k1 = ps ? (dk ? ks + o : ke) : ke;
            const int j0 = ps ? js + o*dj : js + dj, j1 = ps ? (dj ? js + o : je) : je;
            const int i0 = ps ? is + o*di : is + di, i1 = ps ? (di ? is + o : ie) : ie;
            sd[ns] = d;
            sk0[ns] = k0;
            sj0[ns] = j0;
            si0[ns] = i0;
            snj[ns] = j1 - j0 + 1;
            sni[ns] = i1 - i0 + 1;
            c0[ns+1] = c0[ns] + (k1 - k0 + 1)*snj[ns]*sni[ns];
            ++ns;
          }
        }
        const int nper = c0[ns];
        Kokkos::parallel_for("cln_xfe", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             pmy_pack->nmb_thispack*nper),
        KOKKOS_LAMBDA(const int t) {
          const int m = t/nper;
          int c = t - m*nper, q = 0;
          while (c >= c0[q+1]) {++q;}
          c -= c0[q];
          const int i = si0[q] + c % sni[q];
          c /= sni[q];
          const int j = sj0[q] + c % snj[q];
          const int k = sk0[q] + c/snj[q];
          xf(m,sd[q],k,j,i) = FaceCross(*px, src, sd[q], m, k, j, i);
        });
        continue;
      }
      for (int d = 0; d < x.nf; ++d) {
        const int dk = (d == 2), dj = (d == 1), di = (d == 0);
        auto gn = x.g;
        auto fcd = x.fc[d];
        // (inside: FaceCross with no boundary or seam, its own kernel)
        SolveFor("cln_xf", pmy_pack->nmb_thispack - 1, ks + dk, ke, js + dj, je, is + di,
                 ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const int kl = k - dk, jl = j - dj, il = i - di;
          auto lat = [&](const int kc, const int jc, const int ic) {
            Real s;   // (CsLat)
            const int q = d ? CsSlot(*px, d, m, k, j, i, kc, jc, ic, false, s)/2
                            : 0;
            const int p = (q == 0) ? 1 : 0, pp = (q == 2) ? 1 : 2;
            Real r = 0.0;
            for (int c = 0; c < 3; ++c) {
              r += gn.At<CLG_N>(m,3*d+c,k,j,i)*(w(m,CLW_L+3*p+c,kc,jc,ic) +
                                           w(m,CLW_L+3*pp+c,kc,jc,ic));
            }
            return r;
          };
          Real dv = 0.0;
          dv += lat(kl, jl, il) + lat(k, j, i);
          xf(m,d,k,j,i) = fcd(m,CLF_QL+2,k,j,i)*
              (dv/(0.5*(gn(m,CLG_F+3,k,j,i) + gn(m,CLG_F+3,kl,jl,il))));
        });
      }
    }
  }
  auto f = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real a = 0.0;
#pragma unroll
    for (int s = 0; s < CLN_ST; ++s) {
      int dk, dj, di;
      StOff(s, dk, dj, di);
      if ((dj != 0 && !md) || (dk != 0 && !td) || (cub && dk*dk + dj*dj + di*di > 1)) {
        continue;   // (no such neighbour)
      }
      a += st(m,s,k-ks,j-js,i-is)*w(m,src,k+dk,j+dj,i+di);
    }
    if (cxf) {
      const SolveArgs &x = *px;   // (in device memory: no per-thread copy)
      for (int d = 0; d < 3; ++d) {
        for (int sd = 0; sd < 2 && d < x.nf; ++sd) {
          const int kf = k + sd*(d == 2), jf = j + sd*(d == 1), iff = i + sd*(d == 0);
          a += (sd ? x.h : -x.h)*DivW(x, d, m, k, j, i, sd)*xf(m,d,kf,jf,iff);
        }
      }
    }
    w(m,dst,k,j,i) = a;
  };
  OverlapFor("cln_apply_st", pmy_pack, ov, f, fin);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::PreconditionerSetup
//! \brief preconditioner coefficients: the x1-line tridiagonal part of A plus the
//! diagonal of its x2, x3 parts (cross terms left to Krylov); for the mg also the x2, x3
//! couplings and the coarse levels.  Also the face weights and the operator rows

void RadiationVet::PreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  SolveArgs x = MakeArgs();
  auto op = opac;
  // the face weights of the pass, once per face for ApplyOperator (x.fwc)
  CLEAN_HIPX(const SolveArgs *pxw = DevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    SolveFor<0, 3>("cln_fwc", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const SolveArgs &x = *pxw;)   // (in device memory: no per-thread copy)
      FaceWeights(x, dir, m, k, j, i, fc(m,CLF_QL,k,j,i), fc(m,CLF_QL+1,k,j,i));
      const Real t0 = fc(m,CLF_TH,k,j,i);   // (clean_stag_vimp: theta_d, FaceCross)
      const Real th = t0 + (x.drag ? fc(m,CLF_DB,k,j,i) : 0.0)*(1.0 - t0);
      fc(m,CLF_QL+2,k,j,i) = (x.xq >= 0) ? 0.0
                                         : -(1.0 - fc(m,CLF_W,k,j,i))*th*x.h*x.c*x.c;
    });
  }
  x.fwc = true;
  const SolveArgs *px = DevArgs(x, 0);
  // the operator's rows (ApplyOperator)
  auto st = cln_st;
  SolveFor<0, 3>("cln_stencil", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const SolveArgs &x = *px;   // (in device memory: no per-thread copy)
    Real c[CLN_ST];
    StencilRow(x, m, k, j, i, c);
    for (int s = 0; s < CLN_ST; ++s) st(m,s,k-ks,j-js,i-is) = c[s];
  });
  // (lateral stress sums' coefficients: geometry and closure only, built once per
  // closure, cln_csyok)
  if (x.cub && !cln_csyok) {
    cln_csyok = true;
    auto cy = cln_csy;
    par_for("cln_csyb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const SolveArgs &x = *px;   // (in device memory: no per-thread copy)
      for (int q = 0; q < 3; ++q) {
        const int dk = (q == 2), dj = (q == 1), di = (q == 0);
        Real y0[3] = {0.0, 0.0, 0.0}, yn[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        const int vl = TransSide(x, q, m, k, j, i, 0);
        const int vh = TransSide(x, q, m, k, j, i, 1);
        // (seam side: ghost extrapolated in panel, CsSeamCoef)
        const bool o2 = (vl + vh == 1) && x.mbc(m,2*q+vl) == BoundaryFlag::panel;
        const Real wq = (vl + vh == 2 || o2) ? 0.25 : 0.5;
        for (int s = -1; s <= 1; s += 2) {
          const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
          const bool ex = o2 && !TransSide(x, q, m, k, j, i, s > 0);
          if (!ex && !TransSide(x, q, m, k, j, i, s > 0)) continue;
          const Real a = wq*s*x.area[q](m,kf,jf,iff);
          for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) {
              const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
              const Real ag = a*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff);
              if (ex) {   // (c0 P0 + c1 P1 + c2 P2 of the other side's cells)
                Real cf[3];
                CsSeamCoef(x, n, q, m, k, j, i, -s, cf);
                y0[b] += cf[0]*ag*x.cs(m,n,k,j,i);
                yn[s < 0][b] += cf[1]*ag*x.cs(m,n,k-s*dk,j-s*dj,i-s*di);
                yn[s > 0][b] += cf[2]*ag*x.cs(m,n,k-2*s*dk,j-2*s*dj,i-2*s*di);
                continue;
              }
              yn[s > 0][b] += ag*x.cs(m,n,k+s*dk,j+s*dj,i+s*di);
              y0[b] -= ag*x.cs(m,n,k,j,i);
            }
          }
        }
        for (int b = 0; b < 3; ++b) {
          cy(m,6*q+b,k-ks,j-js,i-is) = yn[0][b];
          cy(m,6*q+3+b,k-ks,j-js,i-is) = yn[1][b];
        }
        for (int b = 0; b < 3; ++b) cy(m,18+3*q+b,k-ks,j-js,i-is) = y0[b];
      }
    });
  }
  par_for("cln_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    CLEAN_HIPX(const SolveArgs &x = *px;)   // (in device memory: no per-thread copy)
    Real al_p, ar_p, al_m, ar_m;
    FaceWeights(x, 0, m, k, j, i+1, al_p, ar_p);
    FaceWeights(x, 0, m, k, j, i, al_m, ar_m);
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    const Real pd = x.w(m,CLW_PD,k,j,i);   // (the work, RightHandSide)
    Real a0 = 1.0 + x.h*x.c*op(m,VET_OP_E,k,j,i)*qr + pd + x.h*x.idx[0]*(al_p - ar_m);
    x.w(m,CLW_AM,k,j,i) = -x.h*x.idx[0]*al_m;
    x.w(m,CLW_AP,k,j,i) = x.h*x.idx[0]*ar_p;
    if (x.sph || x.cub) {   // (A_f/V weights, DivW)
      const Real wp = DivW(x, 0, m, k, j, i, 1), wm = DivW(x, 0, m, k, j, i, 0);
      a0 = 1.0 + x.h*x.c*op(m,VET_OP_E,k,j,i)*qr + pd + x.h*(wp*al_p - wm*ar_m);
      x.w(m,CLW_AM,k,j,i) = -x.h*wm*al_m;
      x.w(m,CLW_AP,k,j,i) = x.h*wp*ar_p;
    }
    if (x.multi_d) {
      FaceWeights(x, 1, m, k, j+1, i, al_p, ar_p);
      FaceWeights(x, 1, m, k, j, i, al_m, ar_m);
      a0 += (x.sph || x.cub) ? x.h*(DivW(x, 1, m, k, j, i, 1)*al_p -
                                    DivW(x, 1, m, k, j, i, 0)*ar_m)
                             : x.h*x.idx[1]*(al_p - ar_m);
      const bool cv = x.sph || x.cub;
      x.w(m,CLW_TM,k,j,i) = -x.h*al_m*(cv ? DivW(x, 1, m, k, j, i, 0) : x.idx[1]);
      x.w(m,CLW_TM+1,k,j,i) = x.h*ar_p*(cv ? DivW(x, 1, m, k, j, i, 1) : x.idx[1]);
    }
    if (x.three_d) {
      FaceWeights(x, 2, m, k+1, j, i, al_p, ar_p);
      FaceWeights(x, 2, m, k, j, i, al_m, ar_m);
      const bool cv = x.sph || x.cub;
      a0 += cv ? x.h*(DivW(x, 2, m, k, j, i, 1)*al_p -
                      DivW(x, 2, m, k, j, i, 0)*ar_m)
               : x.h*x.idx[2]*(al_p - ar_m);
      x.w(m,CLW_TM+2,k,j,i) = -x.h*al_m*(cv ? DivW(x, 2, m, k, j, i, 0) : x.idx[2]);
      x.w(m,CLW_TM+3,k,j,i) = x.h*ar_p*(cv ? DivW(x, 2, m, k, j, i, 1) : x.idx[2]);
    }
    x.w(m,CLW_A0,k,j,i) = a0;
  });
  CoarseSetup();
}

//----------------------------------------------------------------------------------------
//! \fn MgLevel RadiationVet::GetMgLevel
//! \brief level l of the staggered mg: 0 the mesh (cln_w), l > 0 cln_mg[l-1]; scalar
//! rows in cln_w's slot order (AM, A0, AP, scratch) on every level

MgLevel RadiationVet::GetMgLevel(const int l) {
  MgLevel lv;
  if (l == 0) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    lv.u = cln_w;
    lv.o[0] = indcs.is;
    lv.o[1] = indcs.js;
    lv.o[2] = indcs.ks;
    lv.n[0] = indcs.nx1;
    lv.n[1] = indcs.nx2;
    lv.n[2] = indcs.nx3;
    lv.a0 = CLW_A0;
  } else {
    lv.u = cln_mg[l-1];
    for (int d = 0; d < 3; ++d) {
      lv.o[d] = 0;
      lv.n[d] = cln_mgn[l-1][d];
    }
    lv.a0 = CLS_A0;
  }
  const int nb[3] = {pmy_pack->pmesh->mb_indcs.nx1, pmy_pack->pmesh->mb_indcs.nx2,
                     pmy_pack->pmesh->mb_indcs.nx3};
  for (int d = 0; d < 3; ++d) {   // (layout_tile; default the MeshBlock)
    lv.t[d] = cln_mgt[0] ? cln_mgt[d]*lv.n[d]/nb[d] : lv.n[d];
  }
  const int tm = (l == 0) ? CLW_TM : CLS_TM;
  lv.sm[0] = lv.a0 - 1;
  lv.sp[0] = lv.a0 + 1;
  for (int d = 1; d < 3; ++d) {
    lv.sm[d] = tm + 2*(d-1);
    lv.sp[d] = tm + 2*(d-1) + 1;
  }
  return lv;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CoarseSetup
//! \brief coarse rows of the staggered mg, Galerkin R A P: P constant over 2 x 2 (x2,x3)
//! aggregates at fixed x1 within a MeshBlock, R = P^T; in-aggregate couplings go onto
//! the diagonal.  Curved meshes: the mesh rows weighted by cell volume

void RadiationVet::CoarseSetup() {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nlev = static_cast<int>(cln_mg.size()) + 1;
  const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                    pmy_pack->pmesh->use_cubed_sphere;
  auto vol = pmy_pack->pcoord->volume;
  for (int lc = 1; lc < nlev; ++lc) {
    const MgLevel f = GetMgLevel(lc-1), c = GetMgLevel(lc);
    const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
    const bool vw = (lc == 1) && curv;

    par_for("cls_galerkin", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1,
            0, c.n[0]-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
      Real a0 = 0.0, am[3] = {0.0, 0.0, 0.0}, ap[3] = {0.0, 0.0, 0.0};
      for (int a2 = 0; a2 < r[2]; ++a2) {
        for (int a1 = 0; a1 < r[1]; ++a1) {
          const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
          const int i = f.o[0] + ic;
          const int a[3] = {0, a1, a2};
          const Real v = vw ? vol(m,k,j,i) : 1.0;
          a0 += v*f.u(m,f.a0,k,j,i);
          for (int d = 0; d < 3; ++d) {
            if (f.n[d] == 1) continue;
            const Real lo = v*f.u(m,f.sm[d],k,j,i), hi = v*f.u(m,f.sp[d],k,j,i);
            if (a[d] > 0) {
              a0 += lo;
            } else {
              am[d] += lo;
            }
            if (a[d] < r[d] - 1) {
              a0 += hi;
            } else {
              ap[d] += hi;
            }
          }
        }
      }
      c.u(m,c.a0,kc,jc,ic) = a0;
      for (int d = 0; d < 3; ++d) {
        c.u(m,c.sm[d],kc,jc,ic) = am[d];
        c.u(m,c.sp[d],kc,jc,ic) = ap[d];
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::LineSolve
//! \brief u(dst) = L^{-1} u(src), L the x1-line tridiagonal part, on the lines of colour
//! col = (j + k) mod 2 (-1: all), with the neighbours of colour != col and <= lim on the
//! RHS.  j, k counted from the nearer edge, so the sweep is
//! mirror-symmetric about the x2, x3 midplanes.  GPU: parallel cyclic reduction, one team
//! per line; host: Thomas

void RadiationVet::LineSolve(const MgLevel &lv, const int src, const int dst,
                             const int col, const int lim) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const MgLevel l = lv;
  const int is = l.o[0], ie = l.o[0] + l.n[0] - 1, cp = l.a0 + 2;
  auto w = l.u;
  const int n1 = l.t[1], n2 = l.t[2];   // (positions in the line's tile)
  const int t0 = l.t[0];   // (layout_tile_x1: the line's x1 segments, uncoupled)
  const int c0 = col;
  auto cl = KOKKOS_LAMBDA(int pj, int pk) {   // (the line's colour)
    pj %= n1;
    pk %= n2;
    return ((2*pj >= n1 ? n1 - 1 - pj : pj) + (2*pk >= n2 ? n2 - 1 - pk : pk)) % 2;
  };
  auto on = KOKKOS_LAMBDA(const int pj, const int pk) {   // (a neighbour on the rhs)
    const int c = cl(pj, pk);
    return c != c0 && c <= lim;
  };

  if (!std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>) {
    // (n2 the lines along x2; t2 the tile's along x3, as the host's Thomas)
    const int nx = l.n[0], n2 = l.n[1], n32 = l.n[2]*l.n[1], t2 = l.t[2];
    int nround = 0;
    while ((1 << nround) < nx) ++nround;
    // (one row per lane: a team of nx rounded to a warp, double-buffered in scratch, one
    // barrier per round; the rows' arithmetic is that of the host's PCR, so bitwise)
    int ts = std::min(1024, (nx + 31)/32*32);
    // (a colour: one team per line of that colour, from a list built once per level)
    int nt = (nmb1 + 1)*n32;
    DvceArray1D<int> ll;
    if (col >= 0) {
      const std::array<int, 7> key = {nt, n2, n1, l.t[2], 2, c0, l.o[1]};
      auto it = cln_pcrl.find(key);
      if (it == cln_pcrl.end()) {
        std::vector<int> h;
        for (int q = 0; q < nt; ++q) {
          const int pk = (q % n32)/n2;
          if (cl(q % n32 - pk*n2, pk) == c0) h.push_back(q);
        }
        DvceArray1D<int> d("cln_pcrl", h.size());
        auto hm = Kokkos::create_mirror_view(d);
        for (size_t q = 0; q < h.size(); ++q) hm(q) = h[q];
        Kokkos::deep_copy(d, hm);
        it = cln_pcrl.emplace(key, d).first;
      }
      ll = it->second;
      nt = ll.extent_int(0);
    }
    const bool lst = (col >= 0);
    const size_t scr = ScrArray1D<Real>::shmem_size(8*nx);
    // (the line of team r: m, j, k, and its rows loaded into the team's scratch sw)
    auto ld = KOKKOS_LAMBDA(const TeamMember_t &tm, const ScrArray1D<Real> &sw, int &m,
                            int &k, int &j) {
      const int q = lst ? ll(tm.league_rank()) : tm.league_rank();
      m = q/n32;
      const int pk = (q - m*n32)/n2, pj = q - m*n32 - pk*n2;
      k = l.o[2] + pk;
      j = l.o[1] + pj;
      if (col >= 0 && cl(pj, pk) != c0) return false;   // (team-uniform)
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
        const int i = is + q;
        Real r = w(m,src,k,j,i);
        if (lim >= 0) {
          if (pj%n1 > 0 && on(pj-1, pk)) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
          if (pj%n1 < n1 - 1 && on(pj+1, pk)) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
          if (pk%t2 > 0 && on(pj, pk-1)) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
          if (pk%t2 < t2 - 1 && on(pj, pk+1)) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
        }
        sw(q) = (q%t0 > 0) ? w(m,l.sm[0],k,j,i) : 0.0;
        sw(nx + q) = w(m,l.a0,k,j,i);
        sw(2*nx + q) = (q%t0 < t0 - 1) ? w(m,l.sp[0],k,j,i) : 0.0;
        sw(3*nx + q) = r;
      });
      tm.team_barrier();
      return true;
    };
    // (a team no larger than this kernel launches with: sm_80 cannot take 1024 here,
    // its registers; each row's arithmetic is the same at any team size)
    static std::map<std::pair<size_t, int>, int> tcap;
    for (;;) {
      Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), nt, ts);
      const int lv0 = TeamScratchLevel(scr, 0, pol.team_size());
      pol.set_scratch_size(lv0, Kokkos::PerTeam(scr));
      auto body = KOKKOS_LAMBDA(TeamMember_t tm) {
        ScrArray1D<Real> sw(tm.team_scratch(lv0), 8*nx);
        int m, k, j;
        if (!ld(tm, sw, m, k, j)) return;
        int o = 0;
        for (int rd = 0, s = 1; rd < nround; ++rd, s *= 2) {
          const int d = 4*nx - o;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
            Real a = 0.0, b = sw(o + nx + q), c = 0.0, r = sw(o + 3*nx + q);
            if (q >= s) {
              const Real f = -sw(o + q)/sw(o + nx + q - s);
              a = f*sw(o + q - s);
              b += f*sw(o + 2*nx + q - s);
              r += f*sw(o + 3*nx + q - s);
            }
            if (q + s < nx) {
              const Real g = -sw(o + 2*nx + q)/sw(o + nx + q + s);
              c = g*sw(o + 2*nx + q + s);
              b += g*sw(o + q + s);
              r += g*sw(o + 3*nx + q + s);
            }
            sw(d + q) = a;
            sw(d + nx + q) = b;
            sw(d + 2*nx + q) = c;
            sw(d + 3*nx + q) = r;
          });
          tm.team_barrier();
          o = d;
        }
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
          w(m,dst,k,j,is+q) = sw(o + 3*nx + q)/sw(o + nx + q);
        });
      };
      auto it = tcap.find(std::make_pair(scr, ts));
      if (it == tcap.end()) {
        it = tcap.emplace(std::make_pair(scr, ts),
                          pol.team_size_max(body, Kokkos::ParallelForTag())).first;
      }
      if (ts > 32 && it->second < ts) {
        ts = std::max(32, it->second/32*32);
        continue;
      }
      AthFor("cln_pcr2", pol, body);
      break;
    }
    return;
  }
  par_for("cln_thomas", DevExeSpace(), 0, nmb1, l.o[2], l.o[2]+l.n[2]-1,
          l.o[1], l.o[1]+l.n[1]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int pj = j - l.o[1], pk = k - l.o[2];
    if (col >= 0 && cl(pj, pk) != c0) return;
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = ((i - is)%t0 > 0) ? w(m,l.sm[0],k,j,i) : 0.0;
      const Real ap = ((i - is)%t0 < t0 - 1) ? w(m,l.sp[0],k,j,i) : 0.0;
      Real r = w(m,src,k,j,i);
      if (lim >= 0) {
        if (pj%n1 > 0 && on(pj-1, pk)) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
        if (pj%n1 < n1 - 1 && on(pj+1, pk)) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
        if (pk%n2 > 0 && on(pj, pk-1)) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
        if (pk%n2 < n2 - 1 && on(pj, pk+1)) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
      }
      const Real den = w(m,l.a0,k,j,i) - am*cprev;
      cprev = ap/den;
      dprev = (r - am*dprev)/den;
      w(m,cp,k,j,i) = cprev;
      w(m,dst,k,j,i) = dprev;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,cp,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VCycle
//! \brief u(dst) = M^{-1} u(src) on level l of the staggered mg: z0 = red-black x1-line
//! sweep, r1 = R (r - A z0) (volume-weighted on curved meshes), z = z0 + P z1 with z1
//! from level l+1; the coarsest level only sweeps.  solver_sweep: &1 a second even
//! pass, &2 a post-sweep (odd, then even)

void RadiationVet::VCycle(const int l, const int src, const int dst) {
  const MgLevel f = GetMgLevel(l);
  LineSolve(f, src, dst, 0, -1);
  LineSolve(f, src, dst, 1, 1);
  if (cln_pcsweep & 1) LineSolve(f, src, dst, 0, 1);
  if (l == static_cast<int>(cln_mg.size())) return;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const MgLevel c = GetMgLevel(l+1);
  const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
  const bool vw = (l == 0) && (pmy_pack->pmesh->use_spherical_polar ||
                                pmy_pack->pmesh->use_cubed_sphere);
  auto vol = pmy_pack->pcoord->volume;

  par_for("cls_restrict", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1, 0, c.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
    Real s = 0.0;
    for (int a2 = 0; a2 < r[2]; ++a2) {
      for (int a1 = 0; a1 < r[1]; ++a1) {
        const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
        const int i = f.o[0] + ic;
        const int pos[3] = {ic%f.t[0], (j - f.o[1])%f.t[1], (k - f.o[2])%f.t[2]};
        Real y = f.u(m,f.a0,k,j,i)*f.u(m,dst,k,j,i);
        for (int d = 0; d < 3; ++d) {
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          if (pos[d] > 0) y += f.u(m,f.sm[d],k,j,i)*f.u(m,dst,k-dk,j-dj,i-di);
          if (pos[d] < f.t[d] - 1) y += f.u(m,f.sp[d],k,j,i)*f.u(m,dst,k+dk,j+dj,i+di);
        }
        s += (vw ? vol(m,k,j,i) : 1.0)*(f.u(m,src,k,j,i) - y);
      }
    }
    c.u(m,CLS_B,kc,jc,ic) = s;
  });
  VCycle(l+1, CLS_B, CLS_X);
  par_for("cls_prolong", DevExeSpace(), 0, nmb1, f.o[2], f.o[2]+f.n[2]-1,
          f.o[1], f.o[1]+f.n[1]-1, f.o[0], f.o[0]+f.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f.u(m,dst,k,j,i) += c.u(m,CLS_X,(k - f.o[2])/r[2],(j - f.o[1])/r[1],i - f.o[0]);
  });
  if (cln_pcsweep & 2) {
    LineSolve(f, src, dst, 1, 1);
    LineSolve(f, src, dst, 0, 1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Precondition
//! \brief w(dst) = M^{-1} w(src) (VCycle)

void RadiationVet::Precondition(const int src, const int dst) {
  VCycle(0, src, dst);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::SolveBiCGStab
//! \brief right-preconditioned BiCGStab for A X = B from the given X (cln_w); returns
//! the iteration count, negative if not converged

int RadiationVet::SolveBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  WorkVec w = cln_w;
  const int nc = 1;
  // newton_conv = budget, in the Newton passes (a single solve keeps solver_tol in its own
  // norm): the rows over their gross terms S (RightHandSide), W A X = W B with
  // W = 1/S, right-preconditioned by M^-1 W^-1, so that every norm below is the budget's
  // (cln_ewr -3: newton_closure_budget's pass, in that norm to 0.1 newton_conv_max)
  const bool bud = cln_bon && (cln_ewr > -1.5 || cln_ewr < -2.5);
  const int sx = CLW_X, sb = bud ? CLW_BW : CLW_B;
  const int sr = CLW_R, sr0 = CLW_R0;
  const int sp = CLW_P, sv = CLW_V;
  const int ss = CLW_S, sq = CLW_Q;
  const int sy = CLW_Y, sz = CLW_Z;

  // w(dst) = sa w(a) + sb w(b)
  auto axpby = [&](const int dst, const int a, const Real sa, const int b,
                   const Real sb_) {
    par_for("cln_axpby", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,dst+c,k,j,i) = sa*w(m,a+c,k,j,i) + sb_*w(m,b+c,k,j,i);
      }
    });
  };
  // w(dst) = w(a) * w(s)^p
  auto scale = [&](const int dst, const int a, const int s, const Real p) {
    par_for("cln_scale", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,dst,k,j,i) = (p > 0.0) ? w(m,a,k,j,i)*w(m,s,k,j,i) : w(m,a,k,j,i)/w(m,s,k,j,i);
    });
  };
  auto apply = [&](const int src, const int dst) {
    ApplyOperator(src, dst);
    if (bud) scale(dst, dst, CLW_SE, -1.0);
  };
  auto precondition = [&](const int src, const int dst) {
    if (bud) {
      scale(CLW_PW, src, CLW_SE, 1.0);
      Precondition(CLW_PW, dst);
    } else {
      Precondition(src, dst);
    }
  };
  if (bud) scale(CLW_BW, CLW_B, CLW_SE, -1.0);
  auto dots = [&](const int n, const std::array<int, 4> ab) {
    return GlobalDots(w, nc, n, ab);
  };

  // (re)start from the true residual with r0 = r.  A breakdown (r0.r at rounding level,
  // r0.v or omega zero) restarts; a second one right after a restart fails.  Updates and
  // their dots share a pass (GlobalSums); four reductions per iteration.  Convergence of
  // the recursive residual is confirmed on the true one (<= 10 tol), else restart
  Real rho = 1.0, alpha = 1.0, omega = 1.0, rr = 0.0, r0r0 = 0.0, rho_new = 0.0, bb = 0.0;
  int it0 = 0;
  // clean_bud_maxnorm: max_i |r_i| of the vector a sum pass just formed, in the same
  // reduction (GlobalSumsMax); rmx of r, smx of s
  const bool mn = bud && cln_bmn;
  Real rmx = 0.0, smx = 0.0;
  auto gr = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    return fabs(w(m,sr,k,j,i));
  };
  auto gs = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    return fabs(w(m,ss,k,j,i));
  };
  // (g: device lambdas made above; none may be defined inside this generic lambda)
  auto sums = [&](const char *nm, const int n, const auto &f, const auto &g) {
    if (!mn) {
      const auto d = GlobalSums(pmy_pack, nm, n, f);
      return std::array<Real, 3>{d[0], d[1], 0.0};
    }
    return GlobalSumsMax(pmy_pack, nm, n, f, g);
  };
  auto restart = [&](const int it, const bool force = false) {
    if (it == it0 + 1 && !force) return false;
    apply(sx, sr);
    const int n = (it < 0) ? 2 : 1;
    const auto d = sums("cln_dot_rst", n,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {   // (r = b - A x, r0 = r, p = v = 0)
        w(m,sr+c,k,j,i) = 1.0*w(m,sb+c,k,j,i) + -1.0*w(m,sr+c,k,j,i);
        w(m,sr0+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sp+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sv+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        if (n > 1) s1 += w(m,sb+c,k,j,i)*w(m,sb+c,k,j,i);
      }
    }, gr);
    rmx = d[2];
    rho = alpha = omega = 1.0;
    rr = r0r0 = rho_new = d[0];
    if (it < 0) bb = d[1];
    it0 = it;
    return true;
  };
  restart(-1);
  Real tol = cln_tol*fmax(std::sqrt(bb), 1.0e-300);
  // newton_maxnorm: the budget's tests on max_i |r_i|, the norm of the pass's own test
  // (an L2 norm above sqrt(N) t has its max above t: not reduced)
  const Real sqn = std::sqrt(static_cast<Real>(pmy_pack->pmesh->nmb_total)*
                             (ie - is + 1)*(je - js + 1)*(ke - ks + 1));

  // Eisenstat-Walker forcing (choice 2, safeguarded) on |r0|, the Newton residual of the
  // pass in the budget's norm
  if (bud && cln_ewr < -2.5) {
    tol = 0.1*cln_bmax;
  } else if (bud) {
    const Real r0 = mn ? rmx : std::sqrt(rr);
    Real eta = 0.1;
    if (cln_ewr > 0.0) {
      eta = 0.9*SQR(r0/cln_ewr);
      if (0.9*SQR(cln_eta) > 0.1) eta = fmax(eta, 0.9*SQR(cln_eta));
      eta = fmin(eta, 0.1);
    }
    cln_ewr = r0;
    cln_eta = eta;
    tol = fmax(eta*r0, 0.1*cln_tolb);
  }
  // not with the budget: also
  // require max_i |r_i|/(a_i E_i) <= cnorm, a_i the diagonal, E_i the max of E^k and |X|
  // over the stencil (an M-matrix row's round-off bound, Oettli-Prager), since the |b|
  // test alone leaves faint regions unconverged.  If this stalls with the global test
  // met, it is at round-off: converged
  constexpr Real cnt = 1.0e-10;
  const Real efl = e_floor;
  auto conv = [&](const int s, const Real r2, const Real vm, const Real t) {
    if (!(std::sqrt(r2) <= t) && !(mn && std::sqrt(r2) <= sqn*t && vm <= t)) {
      return false;
    }
    if (bud) return true;
    GhostFill(w, sx, 1, 0, 1);   // (X's ghosts are those of the last apply of X)
    Real a = 0.0;
    const int dj = pmy_pack->pmesh->multi_d, dk = pmy_pack->pmesh->three_d;
    // (a flat range; a max is order-independent)
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("cln_cnorm", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &mx) {
      const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
      const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
      Real e = fmax(w(m,CLW_EK,k,j,i), efl);
      for (int o = -1; o <= 1; o += 2) {
        e = fmax(e, fmax(fabs(w(m,sx,k,j,i+o)), fabs(w(m,sx,k,j+o*dj,i))));
        e = fmax(e, fmax(fabs(w(m,sx,k+o*dk,j,i)), fabs(w(m,sx,k,j,i))));
      }
      mx = fmax(mx, fabs(w(m,s,k,j,i))/fmax(cnt*fabs(w(m,CLW_A0,k,j,i))*e, 1.0e-300));
    }, Kokkos::Max<Real>(a));
    return GlobalReduce(a, true) <= 1.0;
  };
  if (conv(sr, rr, rmx, tol)) return 0;
  auto done = [&](const int it) {
    return restart(it, true) && conv(sr, rr, rmx, 10.0*tol);
  };

  const Real eps = std::numeric_limits<Real>::epsilon();
  for (int it = 1; it <= cln_maxit; ++it) {
    if (!(fabs(rho_new) > eps*std::sqrt(r0r0*rr))) {
      if (restart(it)) continue;
      return -it;
    }
    const Real beta = (rho_new/rho)*(alpha/omega);
    rho = rho_new;
    par_for("cln_pupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,sp+c,k,j,i) = w(m,sr+c,k,j,i)
                          + beta*(w(m,sp+c,k,j,i) - omega*w(m,sv+c,k,j,i));
      }
    });
    precondition(sp, sy);
    apply(sy, sv);
    const Real r0v = dots(1, {sr0, sv, 0, 0})[0];
    if (!std::isfinite(r0v) || r0v == 0.0) {
      if (restart(it)) continue;
      return -it;
    }
    alpha = rho/r0v;
    const Real ma = -alpha;
    const auto ds = sums("cln_dot_s", 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &) {
      for (int c = 0; c < nc; ++c) {   // (s = r - alpha v)
        w(m,ss+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + ma*w(m,sv+c,k,j,i);
        s0 += w(m,ss+c,k,j,i)*w(m,ss+c,k,j,i);
      }
    }, gs);
    const Real s2 = ds[0];
    smx = ds[2];
    if (conv(ss, s2, smx, tol)) {
      axpby(sx, sx, 1.0, sy, alpha);
      if (done(it)) return it;
      continue;
    }
    precondition(ss, sz);
    apply(sz, sq);
    const auto d = dots(2, {sq, sq, sq, ss});
    omega = (d[0] > 0.0) ? d[1]/d[0] : 0.0;
    const Real al_ = alpha, om = omega;
    const auto e = sums("cln_dot_x", 2,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {
        w(m,sx+c,k,j,i) += al_*w(m,sy+c,k,j,i) + om*w(m,sz+c,k,j,i);
        w(m,sr+c,k,j,i) = w(m,ss+c,k,j,i) - om*w(m,sq+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        s1 += w(m,sr0+c,k,j,i)*w(m,sr+c,k,j,i);
      }
    }, gr);
    rr = e[0];
    rho_new = e[1];
    rmx = e[2];
    if (conv(sr, rr, rmx, tol) && done(it)) return it;
    if (omega == 0.0 && !restart(it)) return -it;
  }
  return -cln_maxit;
}

}  // namespace radvet
