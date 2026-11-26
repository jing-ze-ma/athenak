//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_flux_seam.cpp
//! \brief faces shared by MeshBlocks and panels: one value each; the cubed sphere's
//! lateral stress sums; the drag terms of the faces.
//!
//!  called from: rad_vet_flux.cpp, rad_vet_flux_bounds.cpp, rad_vet_newton.cpp,
//!    rad_vet_oldvector.cpp, rad_vet_rhs.cpp, rad_vet_solver.cpp, rad_vet_tasks.cpp,
//!    rad_vet_update.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_init.cpp

#include <functional>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_flux.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFaceTie
//! \brief the face tie (spherical polar, Cartesian): a face shared by two MeshBlocks (or
//! one across a periodic wrap) takes its owner's values (the block whose low face it is),
//! as its F does: the momentum exchange T = (F - G)/c^2 + h div P - h N1
//! (CleanMomentumChange), the memory d and its slope.  Both sides then give the gas the
//! halves of one T.  Each block computes its low faces' into CLW_L + 3 dir (+ 1 d, + 2
//! slope); the exchange hands them to the other side.  Cubed sphere (and
//! clean_stag_tie_gid): T in CLF_PT, the lower gid's (CleanStagCsTie).

void RadiationVet::CleanStagFaceTie(const int e, const int kq, const bool pt) {
  // d's slope: clean_stag_fkin only (otherwise the face arrays end before CLF_KD)
  const int kd = cln_fkin ? kq : 0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.fs = true;
  const int nf = x.nf;
  for (int d = 0; d < nf && !pt; ++d) {   // as CleanStagDrag; none at physical faces
    auto fc = x.fc[d];
    CleanFor("cln_tiet", nmb1, ks, ke+(d == 2), js, je+(d == 1), is, ie+(d == 0),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, d, m, k, j, i)) return;
      fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) +
          x.h*CleanFaceStress(x, e, d, m, k, j, i, true) - x.h*fc(m,CLF_N+1,k,j,i);
    });
  }
  CleanStagCsTie((kd > 0) ? 3 : 2, CLF_PT, CLF_N + 3, kd);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFCsShared
//! \brief cubed sphere, the F solve's last step: a face shared with another MeshBlock
//! takes the lower gid's solution (in our axes, CleanCsSlot), bitwise the same on both
//! sides (CLF_N + 2), and the cells' slots of cln_vf take it

void RadiationVet::CleanStagFCsShared() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto v = cln_vf;
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor("cln_fcs3", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) {
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);
        return;
      }
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      Real f = lin ? v(m,9+2*dir+1,kl,jl,il) : v(m,9+2*dir,k,j,i);
      if (!(lin && rin)) {   // the other MeshBlock's solution, in our axes
        const int kg = lin ? k : kl, jg = lin ? j : jl, ig = lin ? i : il;
        const int ko = lin ? kl : k, jo = lin ? jl : j, io = lin ? il : i;
        Real s;
        const int n = CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s);
        if (x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io)) f = s*v(m,9+n,kg,jg,ig);
      }
      fc(m,CLF_N+2,k,j,i) = f;
    });
  }
  CleanFor("cln_fcs4", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int q = 0; q < 3; ++q) {
      v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
      v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
    }
  });
  // one MeshBlock along x1: not the x1 family's F (CleanStagFSolve)
  const int f0 =
      (pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel) ? 2 : 0;
  CleanGhostFill(cln_vf, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));   // (as the F solve)
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagCsTie
//! \brief cubed sphere, the face tie: the nc face slots c0, c1, c2 of a face shared with
//! another MeshBlock take the lower gid's values in our axes (CleanCsSlot), bitwise one
//! value on both sides.  Each cell puts its boundary faces' in cln_vf slot 3 n + q, which
//! the ghosts take.  clean_stag_tie_gid on spherical polar / Cartesian: the same in
//! CLW_L + 3 n + q; a periodic self-neighbour takes its low face's (the owner's).

void RadiationVet::CleanStagCsTie(const int nc, const int c0, const int c1,
                                  const int c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const bool cub = x.cub;
  const int v0 = cub ? 0 : CLW_L, nf = x.nf;
  const CleanW v = cub ? CleanW(cln_vf) : cln_w;
  auto gid = pmy_pack->pmb->mb_gid.d_view;
  auto nb = pmy_pack->pmb->nghbr.d_view;
  CleanFor("cln_cstie", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int c[3] = {c0, c1, c2};
    for (int q = 0; q < nf; ++q) {
      const int u = (q == 0) ? (i == ie) : ((q == 1) ? (j == je) : (k == ke));
      for (int n = 0; n < nc; ++n) {
        v(m,v0+3*n+q,k,j,i) = x.fc[q](m,c[n],k+u*(q == 2),j+u*(q == 1),i+u*(q == 0));
      }
    }
  });
  if (cub) {
    CleanGhostFill(cln_vf, 0, 3*nc, true, 0, 1);
  } else {
    CleanGhostFill(cln_w, CLW_L, 3*nc, 0, 1);
  }
  for (int dir = 0; dir < nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor("cln_cstie2", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      if ((lin && rin) || CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const int kg = lin ? k : k - dk, jg = lin ? j : j - dj, ig = lin ? i : i - di;
      const int ko = lin ? k - dk : k, jo = lin ? j - dj : j, io = lin ? i - di : i;
      // sp, C: NeighborIndex of the face, 0/4 x1, 8/12 x2, 24/28 x3
      const int gn = cub ? 0 : nb(m,((dir == 2) ? 24 : 8*dir) + 4*lin).gid;
      if (cub ? !(x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io))
              : !(gn < gid(m) || (gn == gid(m) && lin))) return;
      const int c[3] = {c0, c1, c2};
      Real s = 1.0;
      const int q = cub ? CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s)/2 : dir;
      for (int n = 0; n < nc; ++n) fc(m,c[n],k,j,i) = s*v(m,v0+3*n+q,kg,jg,ig);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagCsLateral
//! \brief staggered, cubed sphere: for E in work vector e, the lateral stress sums of
//! every cell and face family q, Y_q = 1/4 sum over its two faces of A n.(P' - P), P = D
//! E (Cartesian, n outward, P' of the cell across; one-sided where a side's ghosts do not
//! carry the field, CleanTransSide).  A face's F_n changes by n.div P of its half cells:
//! A n.(P_R - P_L) plus the other families' Y (CleanFaceStress), exact for uniform P.
//! The owner computes Y and the ghosts take it (a seam face sees the true geometry).

void RadiationVet::CleanStagCsLateral(const int e, const bool st, const int part) {
  if (!pmy_pack->pmesh->multi_d) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x = CleanMakeArgs();
  const int ks = indcs.ks, js = indcs.js, is = indcs.is;
  auto cy = cln_csy;
  // the sums' array and first slot (off the cubed sphere cln_yl, x.yl)
  const bool cub = pmy_pack->pmesh->use_cubed_sphere, wl = (CLW_L < cln_w.ns);
  auto ya = !cub ? cln_yl : (wl ? cln_w.s : cln_w.o);
  const int y0 = !cub ? 0 : (wl ? CLW_L : CLW_L - cln_w.ns);
  // st: from the pass's coefficients (cln_csyb, CleanPreconditionerSetup)
  if (st) {
    par_for("cln_cslst", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
            indcs.js, indcs.je, indcs.is, indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int q = 0; q < 3; ++q) {
        const int dk = (q == 2), dj = (q == 1), di = (q == 0);
        // a seam side's slot holds the cell 2 away in panel (cln_csyb)
        const int ol = CleanTransSide(x, q, m, k, j, i, 0) ? -1 : 2;
        const int oh = CleanTransSide(x, q, m, k, j, i, 1) ? 1 : -2;
        for (int b = 0; b < 3; ++b) {
          x.w(m,CLW_L+3*q+b,k,j,i) = cy(m,18+3*q+b,k-ks,j-js,i-is)*x.w(m,e,k,j,i) +
              ((cy(m,6*q+b,k-ks,j-js,i-is) != 0.0) ?   // none past a boundary
               cy(m,6*q+b,k-ks,j-js,i-is)*x.w(m,e,k+ol*dk,j+ol*dj,i+ol*di) : 0.0) +
              ((cy(m,6*q+3+b,k-ks,j-js,i-is) != 0.0) ?
               cy(m,6*q+3+b,k-ks,j-js,i-is)*x.w(m,e,k+oh*dk,j+oh*dj,i+oh*di) : 0.0);
        }
      }
    });
    // read by the faces' own cells: depth 1, the two families other than the face's
    // (CleanStagCsLat; part 1: CleanApply finishes it)
    CleanGhostFill(cln_w, CLW_L, 9, part, 1, {0770, 0707, 0077, 0});
    return;
  }
  auto f = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real p0[6];   // the cell's own P, once
    for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
    for (int q = 0; q < 3; ++q) {
      const int dk = (q == 2), dj = (q == 1), di = (q == 0);
      Real y[3] = {0.0, 0.0, 0.0};
      // a side whose ghosts do not carry the field drops out, the other side's weight
      // doubled (one-sided, CleanTransSide)
      const int vl = CleanTransSide(x, q, m, k, j, i, 0);
      const int vh = CleanTransSide(x, q, m, k, j, i, 1);
      // a seam side keeps its face, the ghost's P extrapolated in panel, clamped
      // (CleanCsSeamCoef); an open boundary drops its side, the other's weight doubled
      const bool o2 = (vl + vh == 1) && x.mbc(m,2*q+vl) == BoundaryFlag::panel;
      const Real wq = (vl + vh == 2 || o2) ? 0.25 : 0.5;
      // clean_stag_lat_sums off the cubed sphere: A/V of a Cartesian face 1/dx, V = 1;
      // no family past the mesh's dimension
      for (int s = -1; s <= 1 && q < x.nf; s += 2) {
        const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
        const bool ex = o2 && !CleanTransSide(x, q, m, k, j, i, s > 0);
        if (!ex && !CleanTransSide(x, q, m, k, j, i, s > 0)) continue;
        const Real a = wq*s*((x.cub || x.sph) ? x.area[q](m,kf,jf,iff) : x.idx[q]);
        if (!ex) {   // each component's P(s) - P(0) and the normal once
          Real dp[6], nf[3];
          for (int n = 0; n < 6; ++n) {
            dp[n] = x.cs(m,CLN_D11+n,k+s*dk,j+s*dj,i+s*di)*x.w(m,e,k+s*dk,j+s*dj,i+s*di) -
                    p0[n];
          }
          for (int c = 0; c < 3; ++c) nf[c] = CleanNf(x, c, q, m, kf, jf, iff);
          for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) y[b] += a*nf[c]*dp[(b == c) ? b : 2 + b + c];
          }
          continue;
        }
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
            auto p = [&](const int o) {
              return x.cs(m,n,k+o*dk,j+o*dj,i+o*di)*x.w(m,e,k+o*dk,j+o*dj,i+o*di);
            };
            Real cf[3] = {0.0, 0.0, 0.0};
            if (ex) CleanCsSeamCoef(x, n, q, m, k, j, i, -s, cf);
            y[b] += a*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff)*(ex ? cf[0]*p(0) + cf[1]*p(-s) +
                                                         cf[2]*p(-2*s) : p(s) - p(0));
          }
        }
      }
      for (int b = 0; b < 3; ++b) ya(m,y0+3*q+b,k,j,i) = y[b];
    }
  };
  // the cells it sends first, its exchange in flight behind the inner ones
  const bool ov = (cln_hmax >= 9);
  const std::function<void()> post = [&]() {
    CleanGhostFill(ya, y0, 9, true, ov ? 1 : 0, 1, {0770, 0707, 0077, 0});   // as above
  };
  const std::function<void()> fin = [&]() {
    if (ov) CleanGhostFill(ya, y0, 9, true, 2, 1);
  };
  CleanOverlapFor("cln_cslat", pmy_pack, ov, f, fin, -1, post);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagDrag
//! \brief clean_stag_vimp (staggered, gas coupling): the gas velocity implicit in the
//! face momentum, eliminated face by face so the E solve stays scalar.  A change dm_f of
//! face f's momentum moves its two half cells (mass rho_c/2, share wt_c = chi_c/(chi_L +
//! chi_R)) by dv_c = 2 wt_c dm_f/rho_c, and the advective flux by (1 - theta) dA_f =
//! beta c^2 dm_f, beta = (1 - theta) sum_c wt_c (1 + D_nn,c) E_c/(rho_c c^2).  So
//!   F = (1 - b) F^H(E) + b [Q - h c^2 d_n(D E)],  b = beta/(1 + beta),
//!   Q = F^k + h c^2 d_n(D E^k) = G - c^2 dm_f^k (CLF_DB, CLF_DQ),
//! and the energy flux takes the same increment; at convergence F = F^H (the lagged
//! velocity's solution).  first: dm_f^k = 0, Q = G (with a reference force: Q of E^k, no
//! drag kick); b = 0 while the gas momentum is frozen.  vel (CleanVelocityUpdate with the
//! F solve's F): first the face's CleanMomentumChange term less its share wt (CLF_PT).

void RadiationVet::CleanStagDrag(const bool first, const bool vel, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // cubed sphere: the lateral stress of the iterate (lat false: the caller's)
  if (lat) CleanStagCsLateral(CLW_EK);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  xh.drag = x.drag && !first;   // the flux the iterate's velocity came from
  auto uh = FluidRef::Get(pmy_pack).u0;
  const bool off = cln_vzero || (first && force_ref == VET_FREF_CLEAN && !cln_fref_set);
  const bool lag = !first || (force_ref != VET_FREF_NONE);
  const int nf = x.nf;   // one launch for the faces of all directions
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  CleanFor("cln_drag", nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
    for (int dir = 0; dir < nf; ++dir) {
      const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      if (i > ie+di || j > je+dj || k > ke+dk) continue;
      auto &fc = x.fc[dir];
      if (CleanBoundaryFace(x, dir, m, k, j, i)) continue;   // no F equation
      const Real s = (vel || lag) ? CleanFaceStress(x, CLW_EK, dir, m, k, j, i, true)
                                  : 0.0;
      // vel: the F solve's fs0 just took F^H of E_k with these weights (CLF_PT)
      const Real fq = (vel && lag) ? fc(m,CLF_PT,k,j,i) : 0.0;
      if (vel) {
        fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) + x.h*s -
                             x.h*fc(m,CLF_N+1,k,j,i);
      }
      const Real q = !lag ? fc(m,CLF_G,k,j,i) :
                     (vel ? fq : CleanFaceFlux(xh, CLW_EK, dir, m, k, j, i, true)) +
                     x.h*x.c*x.c*s;
      // cubed sphere: n.D.n, the owners' half-cell V (as CleanMomentumChange) and the
      // exchange's rho (CLW_JV + 2: a seam ghost's own cell), so both sides build one b
      auto vo = [&](const int kk, const int jj, const int ii) {
        return x.cub ? x.g(m,CLG_F+3,kk,jj,ii) : 1.0;
      };
      const Real vl = vo(k-dk, j-dj, i-di), vr = vo(k, j, i);
      const Real xl = x.op(m,VET_OP_T,k-dk,j-dj,i-di)*vl;
      const Real xs = xl + x.op(m,VET_OP_T,k,j,i)*vr, wl = (xs > 0.0) ? xl/xs : 0.5;
      auto ce = [&](const int kk, const int jj, const int ii) {
        const Real d = CleanStagDnn(x, dir, m, k, j, i, kk, jj, ii);
        return (1.0 + d)*fmax(x.w(m,CLW_EK,kk,jj,ii), 0.0)*0.5*(vl + vr)/
               (vo(kk, jj, ii)*(x.cub ? x.w(m,CLW_JV+2,kk,jj,ii) : uh(m,IDN,kk,jj,ii)));
      };
      const Real bt = off ? 0.0 : (1.0 - fc(m,CLF_TH,k,j,i))*(wl*ce(k-dk, j-dj, i-di) +
                      (1.0 - wl)*ce(k, j, i))/(x.c*x.c);
      fc(m,CLF_DB,k,j,i) = bt/(1.0 + bt);
      fc(m,CLF_DQ,k,j,i) = q;
    }
  });
}

}  // namespace radvet
