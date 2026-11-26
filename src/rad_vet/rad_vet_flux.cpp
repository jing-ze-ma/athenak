//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_flux.cpp
//! \brief the face fluxes F of a given E: one line system per x1 line.
//!
//!  called from: rad_vet_newton.cpp, rad_vet_update.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_flux_bounds.cpp, rad_vet_flux_schur.cpp,
//!    rad_vet_flux_seam.cpp, rad_vet_init.cpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_flux.hpp"
#include "utils/p_ledger.hpp"
#include "utils/repro_sum.hpp"
#include "utils/scratch_view.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxSolve
//! \brief clean_stag_fsolve: with E of the solution, the face momentum F alone, one
//! unknown per face:
//!   (F - G)/h + c^2 [P_R - P_L]/dx_f = -c chi_f (F - A) (+ beta terms),
//!   P_i = (D E)_i - al_i [F]_i/(2c),  [F]_i the jump of cell i's two faces,
//! al_i = max(A(tau_i)^2, w of its faces) (1 where |F^H| > c E_donor, not with fkin); al
//! = 0 gives F^H.  Lines are tridiagonal (Thomas per MeshBlock line, ends lagged and
//! iterated with the ghost exchange).  CLF_N holds F, CLF_N + 1 the dissipation term
//! (al_R [F]_R - al_L [F]_L)/(2 c dx_f) the gas takes.
//! clean_stag_fkin: the streaming share w_i takes the rays' kinetic upwind flux,
//!   P_i = (1 - w_i) [(D E)_i - al_i [F]_i/(2c)] + w_i (h+_i F_lo + h-_i F_hi)/c,
//! in area form on spherical polar (F r^2 of a beam steady); telescoping, M-matrix
//! (h+ >= 0 >= h-).  Cubed sphere: a line's cell terms are the face family's (FluxLine,
//! CsFamily), area form (A_R P_R - A_L P_L)/V_f; seams iterate (up to 64 sweeps,
//! 1e-14), a shared face takes the lower gid's F (bitwise); lateral stress stays D E.

void RadiationVet::FluxSolve(const int e, const bool tie_now, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // lat false: the caller's lateral sums are those of e
  if (lat) CsLateralStress(e);
  SolveArgs x = MakeArgs(), xh = x;
  xh.hi = true;
  auto op = opac;
  const FluxLine ln{cln_vf, x.cub};
  const bool fk = cln_fkin, dok = cln_dok, cub = x.cub;
  const bool rms = cln_rstold && !cln_rmdone && fk && e == CLW_X && dok;
  DvceArray1D<Real> rmc("cln_rmc", 1);
  const bool fwb = cln_fwb;
  if (fwb && cln_fwbc.extent_int(0) != 6) Kokkos::realloc(cln_fwbc, 6);
  if (fwb) Kokkos::deep_copy(cln_fwbc, 0.0);
  auto fwc = cln_fwbc;
  const bool wbc = e == CLW_X;
  const int kd = (t2_solve == VET_T2S_NONE) ? 0 :
                 ((t2_solve == VET_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD);
  const Real c1lo = cln_c0, c1hi = cln_c1;
  // the x1 lines by PCR on a GPU; on the cubed sphere with fschur, in panel
  const bool fpcr = (!cub || cln_fschur) &&
      !std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>;
  // the F of every direction (CLW_Y, CLW_Z, CLW_S or cln_vf) in one exchange, started
  // after the last direction.  Cubed sphere, one MeshBlock along x1: the x1 family's
  // slots (0..2, 9, 10, 15) have no line across a MeshBlock face: not exchanged
  const int f0 = (cub && pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel)
                 ? 2 : 0;
  auto yfill = [&](const int part) {
    if (cub) {   // FluxLine reads the faces' own cells, of the face's family (fm)
      GhostFill(cln_vf, 9 + f0, 6 - f0, true, 0, 1, FluxMask(f0));
      return;
    }
    if (x.nf == 3) GhostDefer(cln_w, CLW_S, 1);
    GhostFill(cln_w, CLW_Y, std::min(x.nf, 2), part);
  };
  // three phases (F^H and al; line solves; F), each direction with its own al and F, so
  // each exchange serves all.  clean_stag_fsweeps N: up to N line solves, each with the
  // line ends' F of the last (converging across MeshBlock faces and the periodic wrap);
  // the cubed sphere without clean_stag_fschur iterates its seams (64 by default)
  const int nsw = std::max(cln_fsweeps, 1);
  const bool tie = !cub || cln_fschur;
  const bool fov = !cub && nsw == 1;
  Real fsd = 0.0, fsf = 0.0;
  for (int q = 0; q < (2 + nsw)*x.nf; ++q) {
    const int dir = q%x.nf, ph = std::min(q/x.nf, 1) + (q >= (1 + nsw)*x.nf);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    if (ph == 0) {
      // F^H of every face (CLF_N) and the cell terms of every cell (FluxLine)
      const SolveArgs *px = DevArgs(x, 0), *pxh = DevArgs(xh, 1);
      SolveFor<4>("cln_fs0", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const SolveArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
        DirCase(dir, [&](auto d) {
          constexpr int n = decltype(d)::value;
          fc(m,CLF_N,k,j,i) = FaceFlux<n>(xh, e, n, m, k, j, i, true);
        });
        fc(m,CLF_PT,k,j,i) = fc(m,CLF_N,k,j,i);   // F^H for FaceDrag (QL is free)
        fc(m,CLF_N+1,k,j,i) = 0.0;
        if (!fk || !x.sph || dir == 2 || IsBoundaryFace(x, dir, m, k, j, i)) return;
        // spherical polar: the w share's transverse pressures (but not the normal
        // stress's, q = 4) are the rays' own (CLN_G) where the cell holds their J
        const Real q1 = SphFaceFactor(x, 1, dir, m, k, j, i);
        const Real qn = 2.0*SphFaceFactor(x, 4, dir, m, k, j, i)/q1;
        auto wq = [&](const int q, const int kk, const int jj, const int ii) {
          const Real ec = x.w(m,e,kk,jj,ii);
          const Real rs = SmoothStep(2.5*fmax(ec, 0.0)/
                                     fmax(x.cs(m,CLN_J,kk,jj,ii), 1.0e-300) - 0.25);
          return SmoothStep((x.cs(m,CLN_C1,kk,jj,ii) - c1lo)/(c1hi - c1lo))*rs*
                 (q1*(SphFaceFactor(x, q, dir, m, k, j, i) -
                      qn*x.cs(m,CLN_D11+dir,kk,jj,ii))*ec -
                  x.cs(m,CLN_G+dir,kk,jj,ii)/x.vol(m,kk,jj,ii));
        };
        fc(m,CLF_N+1,k,j,i) += 0.5*(wq(2, k-dk, j-dj, i-di) + wq(3, k, j, i));
      });
      auto v = cln_vf;
      auto f1 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rc;
        const Real tau = op(m,VET_OP_T,k,j,i)*CellFaceDistance(x, dir, m, k, j, i);
        // fkin: the streaming share is kinetic, without dissipation
        Real a = fk ? 0.0 : DissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
        for (int s = 0; s < 2; ++s) {
          const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
          if (!fk) a = fmax(a, fc(m,CLF_W,kf,jf,iff));
          const Real f = fc(m,CLF_N,kf,jf,iff);
          const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
          if (!fk && !IsBoundaryFace(x, dir, m, kf, jf, iff) &&
              !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
          // cubed sphere faces: the last F solve's
          if (cub) v(m,9+2*dir+s,k,j,i) = fc(m,CLF_N+2,kf,jf,iff);
        }
        // fkin: w of the cell; the part of P - D E without F, -w D E, and the stress the
        // kinetic form misses for the rays' own state, (D - s^2) J (0 for a beam, J/3
        // isotropic), lagged; cubed sphere: D about the family's mean normal, area form
        const Real wc = fk ? SmoothStep((x.cs(m,CLN_C1,k,j,i) - c1lo)/(c1hi - c1lo))
                           : 0.0;
        const Real ec = x.w(m,e,k,j,i);
        const Real sr = x.cs(m,CLN_HP+dir,k,j,i) + x.cs(m,CLN_HM+dir,k,j,i);
        const Real jc = fmin(x.cs(m,CLN_J,k,j,i), fmax(ec, 0.0));
        if (!cub) {
          const Real dd = x.cs(m,CLN_D11+dir,k,j,i);
          x.w(m,CLW_L+3*dir,k,j,i) = (1.0 - wc)*a;
          x.w(m,CLW_L+3*dir+1,k,j,i) = wc*((dd - sr*sr)*jc - dd*ec);
          x.w(m,CLW_L+3*dir+2,k,j,i) = wc;
          // F of the low face, for the exchange (multi-D: the last F solve's)
          x.w(m,FluxLine::Cy(dir),k,j,i) = fc(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,
                                              j,i);
          return;
        }
        Real nd[3][3], dd = 0.0;
        CsFamily(x, m, k, j, i, nd);
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            dd += nd[dir][b]*nd[dir][c]*x.cs(m,CLN_D11+((b == c) ? b : 2 + b + c),k,j,i);
          }
        }
        const Real ac = 0.5*(x.area[dir](m,k,j,i) + x.area[dir](m,k+dk,j+dj,i+di));
        v(m,3*dir,k,j,i) = (1.0 - wc)*a;
        v(m,3*dir+1,k,j,i) = wc*dd*(jc - ec);
        v(m,3*dir+2,k,j,i) = ac*wc;
        v(m,15+dir,k,j,i) = -ac*wc*sr*sr*jc;
      };
      // cubed sphere, one MeshBlock along x1: one depth-1 exchange of the family terms
      // (fm), started once the last direction's sent cells are done
      const bool pv = cub && f0 > 0 && dir == x.nf - 1;
      const std::function<void()> post = [&]() {
        GhostDefer(cln_vf, 3, 6, 1);
        cln_hpend.back().fm = {0, 7, 070, 0};
        GhostDefer(cln_vf, 11, 4, 1);
        cln_hpend.back().fm = {0, 3, 014, 0};
        GhostFill(cln_vf, 16, 2, true, 1, 1, {0, 1, 2, 0});
      };
      const std::function<void()> fin = [&]() {GhostFill(cln_vf, 16, 2, true, 2, 1);};
      if (pv) {
        OverlapFor("cln_fs1", pmy_pack, true, f1, fin, -1, post);
      } else {
        SolveFor("cln_fs1", nmb1, ks, ke, js, je, is, ie, f1);
      }
      if (dir == x.nf - 1 && !pv) {
        if (cub) {
          GhostFill(cln_vf, 0, 18, true, 0, 1);
        } else {
          GhostDefer(cln_w, CLW_L, 3*x.nf);
          yfill(0);
        }
      }
    }
    // clean_stag_fschur: the faces beyond the line ends exact, not lagged.
    // clean_stag_fwb: each row less its streaming share w of its residual at the rays'
    // own state F_ray, so F_ray is the thin rows' exact steady state.  F_ray (lab frame):
    // the rays' flux factor times the donor's old E (no advective a E, as Phi's thin
    // share); a physical face its F^H; also the memory G.  The thick share keeps its
    // row, s takes 1/(1 + tau_f), tau_f = chi_f dx_f.  The residual s enters CLF_N + 1:
    // the gas takes the absorption of F, not s, whose s V/c^2 is booked as the
    // radiation's (L_z: dlz).
    // E and the energy flux memory G + d are unchanged.
    if (fwb && q == x.nf) {
      const SolveArgs *pw = DevArgs(x, 0);
      // clean_stag_fwb_ring, _line: s and its weight of the projected faces
      // (FluxBalanceProject)
      DvceArray4D<Real> rs, rw;
      if (cln_fwbring || cln_fwbpj[0]) {
        rs = ScratchView<DvceArray4D<Real>>("cln_fwrs", nmb1 + 1, u0.extent_int(2),
                                            u0.extent_int(3), u0.extent_int(4));
        rw = ScratchView<DvceArray4D<Real>>("cln_fwrw", nmb1 + 1, u0.extent_int(2),
                                            u0.extent_int(3), u0.extent_int(4));
      }
      auto sl = ScratchView<DvceArray1D<Real>>("cln_csl", 6*CLN_NSL);
      if (wbc) Kokkos::deep_copy(sl, 0.0);
      // the F slots of the lines (FluxLine) hold F_ray, then the solve's start again
      for (int pass = 0; pass < 2; ++pass) {
        for (int d = 0; d < x.nf; ++d) {
          const int ek = (d == 2), ej = (d == 1), ei = (d == 0);
          auto fd = x.fc[d];
          SolveFor("cln_fwbs", nmb1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            const SolveArgs &x = *pw;
            for (int s = 0; s < (cub ? 2 : 1); ++s) {
              const int kf = k + s*ek, jf = j + s*ej, iff = i + s*ei;
              Real f = cub ? fd(m,CLF_N+2,kf,jf,iff)
                           : fd(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,j,i);
              if (pass == 0 && !IsBoundaryFace(x, d, m, kf, jf, iff)) {
                f = FluxRay(x, e, d, m, kf, jf, iff);
              }
              if (cub) {
                ln.v(m,9+2*d+s,k,j,i) = f;
              } else {
                x.w(m,FluxLine::Cy(d),k,j,i) = f;
              }
            }
          });
        }
        yfill(0);
        if (pass == 1) break;
        for (int d = 0; d < x.nf; ++d) {
          const int ek = (d == 2), ej = (d == 1), ei = (d == 0);
          auto fd = x.fc[d];
          const bool pj = (cln_fwbring && d == 2) || cln_fwbpj[d];
          if (pj) {
            Kokkos::deep_copy(rs, 0.0);
            Kokkos::deep_copy(rw, 0.0);
          }
          SolveFor("cln_fwb", nmb1, ks, ke+ek, js, je+ej, is, ie+ei,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            const SolveArgs &x = *pw;
            Real lo, hi, dg, r;
            if (IsBoundaryFace(x, d, m, k, j, i)) return;
            if (!FluxRow(x, ln, fd, d, m, k, j, i, lo, hi, dg, r)) return;
            const Real f0 = FluxRay(x, e, d, m, k, j, i);
            const Real fm = ln.Far(x, d, m, k, j, i, k-ek, j-ej, i-ei);
            const Real fp = ln.Far(x, d, m, k, j, i, k, j, i);
            const Real qd = x.drag ? fd(m,CLF_DB,k,j,i) : 0.0;
            const Real tau = 0.5*(x.op(m,VET_OP_T,k,j,i) +
                                  x.op(m,VET_OP_T,k-ek,j-ej,i-ei))*
                             CellFaceDistance(x, d, m, k, j, i);
            const Real s = fd(m,CLF_W,k,j,i)/(1.0 + tau)*
                (lo*fm + dg*f0 + hi*fp - r - (1.0 - qd)*(f0 - fd(m,CLF_G,k,j,i)));
            if (pj) {
              rs(m,k,j,i) = s;
              rw(m,k,j,i) = fd(m,CLF_W,k,j,i)/(1.0 + tau);
              return;
            }
            fd(m,CLF_N+1,k,j,i) += s/(x.h*x.c*x.c);
            const Real p = fd(m,CLF_TH,k,j,i)*s*0.5*(CellVolume(x, m, k, j, i) +
                           CellVolume(x, m, k-ek, j-ej, i-ei))/(x.c*x.c);
            if (!wbc) return;   // counters: the printed s, the write-back's
            SlotAdd(sl, 0, m, k, j, i, 1.0);
            SlotAdd(sl, 1, m, k, j, i, fabs(p));
            if (d == 0) SlotAdd(sl, 2, m, k, j, i, p);
            // (p_ledger: the gas's s V/c^2 withheld, its own cells' halves)
            SlotAdd(sl, 3 + d, m, k, j, i, -0.5*s/(x.c*x.c)*
                    ((((d == 0) ? i <= ie : ((d == 1) ? j <= je : k <= ke)) ?
                      CellVolume(x, m, k, j, i) : 0.0) +
                     (((d == 0) ? i > is : ((d == 1) ? j > js : k > ks)) ?
                      CellVolume(x, m, k-ek, j-ej, i-ei) : 0.0)));
          });
          if (pj) FluxBalanceProject(d, CLF_N + 1, rs, rw, wbc);
        }
      }
      if (wbc) SlotSum(sl, fwc, 6);
    }
    if (cln_fschur && q == x.nf) FluxSchur();
    const int n = (dir == 0) ? ie - is + 1 : ((dir == 1) ? je - js + 1 : ke - ks + 1);
    const int na = (dir == 0) ? je - js + 1 : ie - is + 1;
    const int nb = (dir == 2) ? je - js + 1 : ke - ks + 1;
    Real dmax = 0.0, fmx = 0.0;
    if (ph == 1) {
      // per face (coalesced): lo (CLW_AM), hi (CLW_AP), dg (CLW_A0, then the Thomas den),
      // r (CLW_R0), faces with an equation (CLW_P), and in the sweep c' (CLW_CP); E-solve
      // scratch slots shared by the directions (the x1 PCR reads them too)
      SolveFor("cln_fs2f", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const int t = (dir == 0) ? i - is : ((dir == 1) ? j - js : k - ks);
        // face tie: a line's high face shared with another MeshBlock or the periodic
        // wrap is its owner's (whose low face it is): an identity row (FluxTie)
        const bool own = ln.Own(x, dir, m, k, j, i, t, n, tie);
        Real lo, hi, dg, r;
        const bool eq = FluxRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, own);
        x.w(m,CLW_AM,k,j,i) = lo;
        x.w(m,CLW_AP,k,j,i) = hi;
        x.w(m,CLW_A0,k,j,i) = dg;
        x.w(m,CLW_R0,k,j,i) = r;
        x.w(m,CLW_P,k,j,i) = !own ? -1.0 : (eq ? 1.0 : 0.0);
      });
    }
    if (ph == 1 && dir == 0 && fpcr) {
      // clean_stag_fs_pcr: on a GPU the x1 lines by parallel cyclic reduction, one team
      // per line (as LineSolve); F to round-off of the Thomas solve
      const int nx = n + 1;
      int nround = 0;
      while ((1 << nround) < nx) ++nround;
      const size_t scr = ScrArray1D<Real>::shmem_size(8*nx);
      Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), (nmb1 + 1)*na*nb, 256);
      const int lv0 = TeamScratchLevel(scr, 0, pol.team_size());
      Kokkos::parallel_for("cln_fs2p", pol.set_scratch_size(lv0, Kokkos::PerTeam(scr)),
      KOKKOS_LAMBDA(TeamMember_t tm) {
        const int m = tm.league_rank()/(na*nb), a = tm.league_rank()%na;
        const int k = ks + (tm.league_rank()/na)%nb, j = js + a;
        ScrArray1D<Real> sw(tm.team_scratch(lv0), 8*nx);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int t) {
          const int i = is + t;
          // rows of cln_fs2f; faces beyond a MeshBlock line end: the neighbours', lagged;
          // ftie: a shared high face its owner's, the identity row
          const Real lo = x.w(m,CLW_AM,k,j,i), hi = x.w(m,CLW_AP,k,j,i);
          const Real dg = x.w(m,CLW_A0,k,j,i), p = x.w(m,CLW_P,k,j,i);
          Real r = x.w(m,CLW_R0,k,j,i);
          if (p > 0.0) {
            if (t == 0) r -= lo*ln.Far(x, dir, m, k, j, i, k, j, i-1);
            if (t == n) r -= hi*ln.Far(x, dir, m, k, j, i, k, j, i);
          }
          if (p < 0.0) r = ln.Near(x, dir, m, k, j, i, k, j, i-(t == 0));
          sw(t) = (t > 0) ? lo : 0.0;
          sw(nx + t) = dg;
          sw(2*nx + t) = (t < n) ? hi : 0.0;
          sw(3*nx + t) = r;
        });
        tm.team_barrier();
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
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int t) {
          ln.Put(x, dir, m, k, j, is+t, t, n, sw(o + 3*nx + t)/sw(o + nx + t));
        });
      });
    } else if (ph == 1) {
      // one exact Thomas solve per line, the faces beyond a MeshBlock line end or
      // periodic wrap the last F solve's (lagged); cubed sphere: also the largest face
      // change and largest F.  layout_invariant_mode = tile: each tile of a lateral line
      // solved as a MeshBlock line, its inner ends from cln_fsjt, joined
      // (FluxSchurTile) or the last sweep's; a face two tiles share is its
      // non-owner's identity row
      const int nk = (dir > 0 && repro::tile_n > 0) ? n/repro::tile_n : 1, nt = n/nk;
      if (nk > 1) FluxTileArray(cln_fsjt[dir], nmb1 + 1, na*nb, nk);
      auto jt = cln_fsjt[dir];
      const bool fsc = cln_fschur, olo = !(cub && tie), ohi = cub || !tie;
      auto line = KOKKOS_LAMBDA(const int idx, Real &d_max, Real &f_max) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        for (int u = 1; u < nk && !fsc; ++u) {
          int k, j, i;
          face(u*nt, k, j, i);
          jt(m,l,25*nk+u-1) = ohi ? ln.Far(x, dir, m, k, j, i, k, j, i)
                                  : ln.Near(x, dir, m, k, j, i, k, j, i);
          jt(m,l,24*nk+u) = olo ? ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di)
                                : ln.Near(x, dir, m, k, j, i, k-dk, j-dj, i-di);
        }
        for (int u = 0; u < nk; ++u) {
          const int o = u*nt;
          const bool il = u > 0, ih = u < nk - 1;
          Real dp = 0.0, cp = 0.0;
          for (int t = 0; t <= nt; ++t) {
            int k, j, i;
            face(o + t, k, j, i);
            const bool ie = (il && t == 0) || (ih && t == nt);
            const bool io = ie && !((t == 0) ? olo : ohi);
            const Real ee = ie ? jt(m,l,((t == 0) ? 24*nk : 25*nk)+u) : 0.0;
            const Real am = io ? 0.0 : x.w(m,CLW_AM,k,j,i);
            const Real ap = io ? 0.0 : x.w(m,CLW_AP,k,j,i);
            const Real pp = io ? -1.0 : x.w(m,CLW_P,k,j,i);
            Real den = io ? 1.0 : x.w(m,CLW_A0,k,j,i);   // dg
            // faces beyond a MeshBlock line end: the neighbours', lagged
            const bool eq = pp > 0.0;
            den -= ((eq && t == 0) ? 0.0 : am)*cp;
            cp = ((eq && t == nt) ? 0.0 : ap)/den;
            x.w(m,CLW_CP,k,j,i) = cp;
            const int t0 = (t == 0) ? 1 : 0;
            Real r = (pp < 0.0) ? (ie ? ee : ln.Near(x, dir, m, k, j, i, k-t0*dk, j-t0*dj,
                                                     i-t0*di)) :
                x.w(m,CLW_R0,k,j,i);
            if (pp > 0.0) {
              if (t == 0) {
                r -= am*(il ? ee : ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di));
              }
              if (t == nt) r -= ap*(ih ? ee : ln.Far(x, dir, m, k, j, i, k, j, i));
            }
            dp = (r - ((t == 0) ? 0.0 : am)*dp)/den;
            x.w(m,CLW_R,k,j,i) = dp;
          }
          Real fn = 0.0;
          // the last face's slot is a ghost: the neighbour's or a physical face
          for (int t = nt; t >= 0; --t) {
            int k, j, i;
            face(o + t, k, j, i);
            fn = x.w(m,CLW_R,k,j,i) - ((t < nt) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
            if (cub) {
              d_max = fmax(d_max, fabs(fn - ln.Get(x, dir, m, k, j, i, t)));
              f_max = fmax(f_max, fmax(fabs(fn), 1.0e-8*x.c*fabs(x.w(m,e,k,j,i))));
            }
            ln.Put(x, dir, m, k, j, i, t, nt, fn);
          }
        }
      };
      // one warp or wavefront per block: few x1 lines, one thread each, latency bound
#if defined(KOKKOS_ENABLE_HIP)
      using LineBounds = Kokkos::LaunchBounds<64,1>;
#else
      using LineBounds = Kokkos::LaunchBounds<32,1>;
#endif
      const auto pol = Kokkos::RangePolicy<DevExeSpace, LineBounds>(DevExeSpace(), 0,
                                                                    (nmb1 + 1)*na*nb);
      if (cub && nsw > 1) {
        AthReduce("cln_fs2", pol, line, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(fmx));
      } else {
        AthFor("cln_fs2", pol, KOKKOS_LAMBDA(const int idx) {
          Real d_max = 0.0, f_max = 0.0;
          line(idx, d_max, f_max);
        });
      }
    }
    if (ph == 1) {
      fsd = (dir == 0) ? dmax : fmax(fsd, dmax);
      fsf = (dir == 0) ? fmx : fmax(fsf, fmx);
    }
    if (ph == 1 && dir == x.nf - 1) {
      yfill(fov ? 1 : 0);
      // sweeps: the cubed sphere until no face changes by 1e-14 of the largest F, else
      // until the owners' high-face F changes by < 1e-12 of itself
      const bool done = nsw > 1 &&
          (cub ? GlobalReduce(fsd, true) <= 1.0e-14*GlobalReduce(fsf, true)
               : (FluxSweepChange() < 1.0e-12 && q >= 2*x.nf));
      if (done) q = (1 + nsw)*x.nf - 1;
    }
    // cubed sphere: a face shared with another MeshBlock takes the lower gid's F, and the
    // cells' slots of it
    if (cub && q == (1 + nsw)*x.nf) FluxCsShared();
    if (ph == 2) {
      // F and the dissipation term of every face
      const SolveArgs *px = DevArgs(x, 0);
      auto f3 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const SolveArgs &x = *px;   // in device memory: no per-thread copy
        const Real f = cub ? fc(m,CLF_N+2,k,j,i) : x.w(m,FluxLine::Cy(dir),k,j,i);
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);   // F^H: the boundary's, kept
        if (IsBoundaryFace(x, dir, m, k, j, i)) return;
        // fkin, the solution's F solve: d = Phi^H - F, Phi^H = F^H + theta d the energy
        // flux's high-order part, memory G + d (FaceSource): its own history, so its
        // steady state does not depend on h
        if (fk && e == CLW_X && dok) {   // and its slope (FluxMemoryInc)
          const Real th = fc(m,CLF_TH,k,j,i);
          const Real dd = fc(m,CLF_N,k,j,i) - f - (1.0 - th)*fc(m,CLF_N+3,k,j,i);
          fc(m,CLF_N+3,k,j,i) += dd;
          if (kd > 0) fc(m,kd,k,j,i) = dd/x.h;
          // clean_rst_mshare, the first solution after a restart from an old scheme's
          // file: the memory G + d at the fixed point of Phi^H = theta (G + d) + B of
          // this E, the VET share's steady flux B/(1 - theta) (above c E of its donor:
          // counted), so the restart has no relaxation transient
          if (rms && th < 1.0) {
            const Real ms = (fc(m,CLF_N,k,j,i) - th*fc(m,CLF_G,k,j,i))/(1.0 - th);
            if (fabs(ms) > x.c*fmax(x.w(m,e,k,j,i), x.w(m,e,k-dk,j-dj,i-di))) {
              Kokkos::atomic_add(&rmc(0), 1.0);
            }
            fc(m,CLF_N+3,k,j,i) = ms - f;
          }
        }
        fc(m,CLF_N,k,j,i) = f;
        fc(m,CLF_N+2,k,j,i) = f;
        Real tl[6], tr[6];
        ln.Cell(x, dir, m, k, j, i, k-dk, j-dj, i-di, tl);
        ln.Cell(x, dir, m, k, j, i, k, j, i, tr);
        if (!(tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0)) return;
        // a physical high face keeps F^H; the ghost fill replaced its slot
        const Real fm = ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di);
        const Real fp = ln.Far(x, dir, m, k, j, i, k, j, i);
        auto pk = [&](const Real t[6], const Real l, const Real h) {
          return t[2]*(t[3]*l + t[4]*h)/x.c + (cub ? t[5] : t[1]);
        };
        const Real pl = pk(tl, fm, f), pr = pk(tr, f, fp);
        // the kinetic stress: the cubed sphere's area form over V_f, else the F row's
        Real dp, ds = 0.0;
        if (cub) {
          dp = (pr - pl + CsNormA(x, dir, m, k, j, i)*(tr[1] - tl[1]))/
               (0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,k-dk,j-dj,i-di)));
        } else {
          dp = (x.sph ? SphFaceFactor(x, 0, dir, m, k, j, i) : x.idx[dir])*(pr - pl);
          ds = x.sph ? SphFaceFactor(x, 4, dir, m, k, j, i)*(pr + pl) : 0.0;
        }
        const Real dc = 2.0*x.c*CellFaceDistance(x, dir, m, k, j, i);
        fc(m,CLF_N+1,k,j,i) += (tr[0]*(fp - f) - tl[0]*(f - fm))/dc - dp - ds;
      };
      // x1: the faces that read no ghost while the exchange of F finishes
      const bool ov = fov && dir == 0;
      const std::function<void()> fin = [&]() {
        if (ov) GhostFill(cln_w, CLW_Y, 1, 2);
      };
      OverlapFor("cln_fs3", pmy_pack, ov, f3, fin, dir);
    }
  }
  if (cln_fcone) FluxCone(e, kd);
  if (tie_now) FluxTie(e, kd);
  if (fwb && e == CLW_X) {   // clean_stag_fwb: radiation momentum of the residuals
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_fwbc);
    const Real n = GlobalReduce(h(0), false), p = GlobalReduce(h(1), false);
    const Real p1 = GlobalReduce(h(2), false);
    if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(p + p1))) {
      std::printf("CLEAN %sFWB cycle %d t2s %d faces %.0f radiation momentum sum |dp| "
                  "%.6e x1 faces sum dp %.6e\n", std::isfinite(p + p1) ? "" : "WARNING ",
                  pmy_pack->pmesh->ncycle, t2_solve, n, p, p1);
    }
    if (p_ledger::on && cln_feedback && !cln_vzero) {   // (the gas's s, unprojected)
      const Real a[3] = {GlobalReduce(h(3), false), GlobalReduce(h(4), false),
                         GlobalReduce(h(5), false)};
      p_ledger::Withheld(1, a, PSlope(t2_solve), dt_sub, !t2_inprobe);
    }
  }
  if (rms) {   // clean_rst_mshare: once, at the old file's first solution
    cln_rmdone = true;
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rmc);
    const Real n = GlobalReduce(h(0), false);
    if (global_variable::my_rank == 0) {
      std::printf(
          "<rad_vet> clean_rst_mshare RESTART CONVERSION (one time): the VET share's "
          "memory at its steady flux (faces above c E of the donor: %.0f)\n", n);
    }
  }
}

}  // namespace radvet
