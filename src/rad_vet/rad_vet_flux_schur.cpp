//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_flux_schur.cpp
//! \brief clean_stag_fschur: the face flux lines exact across MeshBlock ends.
//!
//!  called from: rad_vet_flux.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_init.cpp

#include <algorithm>
#include <array>
#include <cstdint>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_faces.hpp"
#include "rad_vet/rad_vet_flux.hpp"
#include "utils/repro_sum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::FluxSweepChange
//! \brief clean_stag_fsweeps: the largest relative change of the F the lines' high faces
//! took from their neighbours (CLW_Y, CLW_Z, CLW_S ghosts) since the last sweep (kept in
//! CLW_TM + dir; the first sweep only stores it)

Real RadiationVet::FluxSweepChange() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nx3 = indcs.nx3, n12 = nx1*nx2, n13 = nx1*nx3, n23 = nx2*nx3;
  const int nfc = n23 + n13 + n12;
  SolveArgs x = MakeArgs();
  const int nf = x.nf;
  Real dmax = 0.0;
  Kokkos::parallel_reduce("cln_fsweep", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nfc),
  KOKKOS_LAMBDA(const int idx, Real &mx) {
    const int m = idx/nfc, q = idx%nfc;
    const int d = (q < n23) ? 0 : ((q < n23 + n13) ? 1 : 2);
    if (d >= nf) return;
    int k, j, i;
    if (d == 0) {i = ie + 1; j = js + q%nx2; k = ks + q/nx2;}
    if (d == 1) {j = je + 1; i = is + (q - n23)%nx1; k = ks + (q - n23)/nx1;}
    if (d == 2) {k = ke + 1; i = is + (q - n23 - n13)%nx1; j = js + (q - n23 - n13)/nx1;}
    if (IsBoundaryFace(x, d, m, k, j, i)) return;
    const int cy = (d == 2) ? CLW_S : CLW_Y + d;
    const Real f = x.w(m,cy,k,j,i), f0 = x.w(m,CLW_TM+d,k,j,i);
    mx = fmax(mx, fabs(f - f0)/fmax(fabs(f), 1.0e-300));
    x.w(m,CLW_TM+d,k,j,i) = f;
  }, Kokkos::Max<Real>(dmax));
  return GlobalReduce(dmax, true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxSchur
//! \brief clean_stag_fschur (uniform mesh; with the face tie or the cubed sphere's owner
//! rule): the F solve's lines exact across MeshBlocks and periodic wraps (cubed sphere:
//! round the four-panel rings).  Each MeshBlock line is a segment, its end faces affine
//! in the faces beyond (FluxRowOwn).  Exchange rounds join segments (FluxRowJoin):
//! a MeshBlock sends right its own joined to the one from the left (and mirrored) until
//! they span the line, or once round a periodic one; the faces beyond its ends then
//! replace the lagged ghosts in the line solves.  A segment goes out in three slots of
//! three cells at each end; the six free slots serve two directions at a time.  Across
//! a cubed-sphere seam a reversed neighbour line (CsSlot) is read reversed.

void RadiationVet::FluxSchur() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  SolveArgs x = MakeArgs();
  const FluxLine ln{cln_vf, x.cub};
  const int nn[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nm[3] = {pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3};
  const int fs[6] = {CLW_TM, CLW_TM+1, CLW_TM+2, CLW_TM+3, CLW_BW, CLW_PW};
  // rounds to the left and right: a line of nb MeshBlocks nb - 1 each, a periodic one
  // nb/2 and (nb - 1)/2; the group of directions exchanged together, and their slots
  // (cubed sphere: x2 and x3 lines in one group, as a seam may turn one into the other)
  const bool cub = x.cub;
  int nr[3][2], gr[3], gs[3] = {0, 0, 0}, ngd = 0;
  bool per[3];
  for (int dd = 0; dd < x.nf; ++dd) {
    const int d = cub ? (dd + 1)%3 : dd;
    const int nb = (cub && d > 0) ? 4*nm[d]/nn[d] : nm[d]/nn[d];
    per[d] = (cub && d > 0) || (pm->mesh_bcs[2*d] == BoundaryFlag::periodic);
    nr[d][0] = per[d] ? nb/2 : nb - 1;
    nr[d][1] = per[d] ? (nb - 1)/2 : nb - 1;
    gr[d] = (nr[d][0] > 0) ? ngd/2 : ((per[d]) ? 0 : -1);
    gs[d] = (nr[d][0] > 0) ? 3*(ngd++%2) : 0;
  }
  // layout_invariant_mode = tile, a lateral line of more than one tile
  const int tn = repro::tile_n;
  if (tn > 0 && (nn[1] > tn || (x.nf == 3 && nn[2] > tn))) {
    FluxSchurTile(gr, gs, per, ngd);
    return;
  }
  for (int g = 0; g == 0 || 2*g < ngd; ++g) {
    int nrt = 0, nd = 0;
    for (int d = 0; d < x.nf; ++d) {
      if (gr[d] == g && nr[d][0] > 0) nrt = std::max(nrt, nr[d][0]), ++nd;
    }
    for (int r = 1; r <= nrt + 1; ++r) {
      const bool fin = (r == nrt + 1);
      for (int dir = 0; dir < x.nf; ++dir) {
        const int rl = nr[dir][0], rr = nr[dir][1];
        if (gr[dir] != g || (!fin && (rl == 0 || (r > rl && r != rr + 1)))) continue;
        const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
        const int nb = (dir == 2) ? nn[1] : nn[2];
        auto fc = x.fc[dir];
        const int c0 = fs[gs[dir]], c1 = fs[gs[dir]+1], c2 = fs[gs[dir]+2];
        const int g1 = gs[1], g2 = gs[2];
        const bool pd = per[dir];
        const int rs = cub ? rl : rr;   // rounds the low end sends
        AthForLw("cln_fsch", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             (nmb1 + 1)*na*nb),
        KOKKOS_LAMBDA(const int idx) {
          const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
          auto cell = [&](const int t, int &k, int &j, int &i) {
            k = (dir == 2) ? x.ks + t : x.ks + b;
            j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
            i = (dir == 0) ? x.is + t : x.is + a;
          };
          // number q of a segment: slot q/3, cell t0 + sg (q%3)
          auto slot = [&](const int t0, const int sg, const int q) -> Real & {
            int k, j, i;
            cell(t0 + sg*(q%3), k, j, i);
            return x.w(m,(q < 3) ? c0 : ((q < 6) ? c1 : c2),k,j,i);
          };
          auto keepo = [&](const int q) -> Real & {
            int k, j, i;
            cell(q/4, k, j, i);
            return x.w(m,((gs[dir] == 0) ? CLW_AM : CLW_R) + q%4,k,j,i);
          };
          // a segment in the ghosts beyond the low or high end; cubed sphere: in the
          // neighbour family's slots, reversed (F -> -F) where it runs the other way
          auto ghost = [&](const bool he, Real g[8]) {
            int k, j, i, kc, jc, ic;
            cell(he ? n : -1, k, j, i);
            cell(he ? n : 0, kc, jc, ic);
            Real sg = 1.0;
            const int fq = cub ? CsSlot(x, dir, m, kc, jc, ic, k, j, i, false, sg)/2
                               : dir;
            const int o = (fq == dir) ? gs[dir] : ((fq == 1) ? g1 : g2);
            for (int q = 0; q < 8; ++q) {
              cell(he ? n + q%3 : -1 - q%3, k, j, i);
              g[q] = x.w(m,fs[o + q/3],k,j,i);
            }
            if (sg > 0.0) return;
            const Real h[8] = {g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]};
            g[0] = -h[4], g[1] = h[6], g[2] = h[5], g[3] = h[7];
            g[4] = -h[0], g[5] = h[2], g[6] = h[1], g[7] = h[3];
          };
          int k, j, i;
          cell(0, k, j, i);
          const bool bl = IsBoundaryFace(x, dir, m, k, j, i);
          cell(n, k, j, i);
          const bool bh = IsBoundaryFace(x, dir, m, k, j, i);
          const bool hl = r > 1 && rl > 0 && !bl, hr = r > 1 && rr > 0 && !bh;
          Real o[8], sl[8], sr[8], c[8], cr[8], ya[3], xb[3], yr[3], xr[3];
          // the own segment once, kept for the rounds in cells 0, 1 of four E-solve
          // scratch slots, a set per direction of the group (a one-MeshBlock periodic
          // line sends nothing: only in the last round)
          const bool r1 = r == 1 || rl == 0;
          for (int q = 0; q < 8 && !r1; ++q) o[q] = keepo(q);
          if (r1) FluxRowOwn(x, ln, fc, dir, m, a, b, n, o);
          for (int q = 0; q < 8 && r1 && !fin; ++q) keepo(q) = o[q];
          // the segment to the left in the ghosts; to the right kept after rr (cubed
          // sphere: in the faces' free CLF_N + 2, both ends send for rl rounds)
          auto keep = [&](const int q) -> Real & {
            int k, j, i;
            cell(q, k, j, i);
            return fc(m,CLF_N+2,k,j,i);
          };
          if (hl) ghost(false, sl);
          if (hr && (r - 1 <= rr || cub)) ghost(true, sr);
          for (int q = 0; q < 8; ++q) {
            if (!hl) sl[q] = 0.0;
            if (hr && r - 1 > rr && !cub) sr[q] = slot(0, 1, q);
            if (!hr) sr[q] = 0.0;
            c[q] = o[q], cr[q] = o[q];
          }
          if (hl) FluxRowJoin(sl, o, c, ya, xb);
          if (!fin) {
            if (hr && r <= rs) FluxRowJoin(o, sr, cr, yr, xr);
            for (int q = 0; q < 8; ++q) {
              if (r <= rl) slot(n - 1, -1, q) = c[q];
              if (r <= rs + ((hr && !cub) ? 1 : 0)) {
                slot(0, 1, q) = (r <= rs) ? cr[q] : sr[q];
              }
              if (cub && hr && r == rr + 1) keep(q) = sr[q];
            }
            return;
          }
          for (int q = 0; q < 8 && cub && hr && r - 1 > rr; ++q) sr[q] = keep(q);
          // the whole line (periodic: determinant (1 - c2)(1 - c5) - c1 c6 as a sum)
          for (int q = 0; q < 8; ++q) cr[q] = c[q];
          if (hr) FluxRowJoin(c, sr, cr, yr, xr);
          Real xo = 0.0, yo = 0.0;
          if (pd) {
            const Real det = cr[3]*(cr[6] + cr[7]) + cr[1]*cr[7];
            xo = (cr[0]*(cr[6] + cr[7]) + cr[1]*cr[4])/det;
            yo = (cr[4]*(cr[1] + cr[3]) + cr[6]*cr[0])/det;
          }
          const Real gr = hr ? xr[0] + xr[1]*yo + xr[2]*xo : xo;
          const Real gl = hl ? ya[0] + ya[1]*yo + ya[2]*gr : yo;
          // the faces beyond the own ones (next block's first xR, last block's last yL),
          // in the slots the line solve reads
          Real sg;
          int kf, jf, iff;
          if (!bh) {   // cell n
            const bool on = ln.Own(x, dir, m, k, j, i, n, n, true);
            Real &f = ln.Slot(x, dir, m, k, j, i, k, j, i, on, sg);
            f = sg*gr;
          }
          cell(0, kf, jf, iff);
          cell(-1, k, j, i);
          if (!bl) {
            const bool on = ln.Own(x, dir, m, kf, jf, iff, 0, n, true);
            Real &f = ln.Slot(x, dir, m, kf, jf, iff, k, j, i, on, sg);
            f = sg*gl;
          }
        });
      }
      // the ghosts of a face take its line family's slots (GSeg::fm)
      if (r <= nrt) {
        std::array<uint32_t, 4> fa = {0, 0, 0, 0}, fb = fa;
        for (int d = 0; d < x.nf; ++d) {
          if (gr[d] == g) fa[d] = (gs[d] > 0) ? 3 : 0, fb[d] = (gs[d] > 0) ? 8 : 7;
        }
        if (nd == 2) GhostDefer(cln_w, CLW_BW, 2), cln_hpend.back().fm = fa;
        GhostDefer(cln_w, CLW_TM, (nd == 2) ? 4 : 3);
        cln_hpend.back().fm = fb;
        GhostFlush();
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluxSchurTile
//! \brief layout_invariant_mode = tile: FluxSchur with the nk tiles of a lateral
//! MeshBlock line as MeshBlocks, so the F solve is the same for every whole-tile layout.
//! Tile g's whole line is join(join(L, o_g), R), L the left fold of the nt/2 segments
//! before it (periodic; else all), R the right fold of the (nt - 1)/2 after.  Joins do
//! not associate in round-off, so the folds are built in that order: in pass p each
//! MeshBlock sends right the left fold from its tile p and left the right fold to its
//! tile nk - 1 - p, extended by each block they pass.  Inner tile ends go to cln_fsjt,
//! the MeshBlock's to the ghosts.

void RadiationVet::FluxSchurTile(const int gr[3], const int gs[3],
                                 const bool per[3], const int ngd) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  SolveArgs x = MakeArgs();
  const FluxLine ln{cln_vf, x.cub};
  const bool cub = x.cub;
  const int nn[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nm[3] = {pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3};
  const int fs[6] = {CLW_TM, CLW_TM+1, CLW_TM+2, CLW_TM+3, CLW_BW, CLW_PW};
  const int g1 = gs[1], g2 = gs[2];
  int nk[3] = {1, 1, 1}, nbk[3] = {1, 1, 1};
  for (int d = 0; d < x.nf; ++d) {
    if (d > 0) nk[d] = nn[d]/repro::tile_n;
    nbk[d] = ((cub && d > 0) ? 4 : 1)*nm[d]/nn[d];
  }
  // the tiles' segments
  for (int dir = 0; dir < x.nf; ++dir) {
    if (nk[dir] == 1 && gr[dir] < 0) continue;
    const int nkd = nk[dir], nt = nn[dir]/nkd, na = (dir == 0) ? nn[1] : nn[0];
    const int nb = (dir == 2) ? nn[1] : nn[2];
    FluxTileArray(cln_fsjt[dir], nmb1 + 1, na*nb, nkd);
    auto jt = cln_fsjt[dir];
    auto fc = x.fc[dir];
    AthFor("cln_fsjo", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
      for (int u = 0; u < nkd; ++u) {
        Real s[8];
        FluxRowOwn(x, ln, fc, dir, m, a, b, nt, s, u*nt, u > 0, u < nkd - 1);
        for (int q = 0; q < 8; ++q) jt(m,l,8*u+q) = s[q];
      }
    });
  }
  // the passes of rounds of each group of directions (FluxSchur)
  for (int g = 0; g == 0 || 2*g < ngd; ++g) {
    int np = 0, nrt = 0, nd = 0;
    for (int d = 0; d < x.nf; ++d) {
      if (gr[d] != g || nbk[d] < 2) continue;
      np = std::max(np, per[d] ? nk[d] : 1);
      nrt = std::max(nrt, per[d] ? nbk[d]/2 : nbk[d] - 1);
      ++nd;
    }
    for (int p = 0; p < np; ++p) {
      for (int r = 1; r <= nrt + 1; ++r) {
        for (int dir = 0; dir < x.nf; ++dir) {
          const bool pd = per[dir];
          const int nkd = nk[dir], ns = pd ? nbk[dir]/2 : nbk[dir] - 1;
          if (gr[dir] != g || nbk[dir] < 2 || p >= (pd ? nkd : 1) || r > ns + 1) continue;
          const int ul = pd ? p : 0, ur = pd ? (nkd - p)%nkd : 0;
          // the round the right fold is kept (tile 0 of a periodic line has (nt - 1)/2
          // tiles after it)
          const int rk = (pd && p == 0) ? ((nbk[dir]*nkd - 1)/2)/nkd + 1 : ns + 1;
          const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
          const int nb = (dir == 2) ? nn[1] : nn[2];
          const int c0 = fs[gs[dir]], c1 = fs[gs[dir]+1], c2 = fs[gs[dir]+2];
          const int gd = gs[dir];
          auto jt = cln_fsjt[dir];
          AthFor("cln_fsjr", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
          KOKKOS_LAMBDA(const int idx) {
            const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
            auto cell = [&](const int t, int &k, int &j, int &i) {
              k = (dir == 2) ? x.ks + t : x.ks + b;
              j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
              i = (dir == 0) ? x.is + t : x.is + a;
            };
            auto slot = [&](const int t0, const int sg, const int q) -> Real & {
              int k, j, i;
              cell(t0 + sg*(q%3), k, j, i);
              return x.w(m,(q < 3) ? c0 : ((q < 6) ? c1 : c2),k,j,i);
            };
            auto ghost = [&](const bool he, Real g[8]) {   // as FluxSchur
              int k, j, i, kc, jc, ic;
              cell(he ? n : -1, k, j, i);
              cell(he ? n : 0, kc, jc, ic);
              Real sg = 1.0;
              const int fq = cub ? CsSlot(x, dir, m, kc, jc, ic, k, j, i, false,
                                          sg)/2 : dir;
              const int o = (fq == dir) ? gd : ((fq == 1) ? g1 : g2);
              for (int q = 0; q < 8; ++q) {
                cell(he ? n + q%3 : -1 - q%3, k, j, i);
                g[q] = x.w(m,fs[o + q/3],k,j,i);
              }
              if (sg > 0.0) return;
              const Real h[8] = {g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]};
              g[0] = -h[4], g[1] = h[6], g[2] = h[5], g[3] = h[7];
              g[4] = -h[0], g[5] = h[2], g[6] = h[1], g[7] = h[3];
            };
            int k, j, i;
            cell(0, k, j, i);
            const bool xl = r > 1 && !IsBoundaryFace(x, dir, m, k, j, i);
            cell(n, k, j, i);
            const bool xr = r > 1 && !IsBoundaryFace(x, dir, m, k, j, i);
            Real sl[8], sr[8], f[8], h[8], s[8], ya[3], xb[3];
            if (xl) ghost(false, sl);
            if (xr) ghost(true, sr);
            for (int q = 0; q < 8; ++q) {
              if (xl && r == ns + 1) jt(m,l,8*nkd+8*ul+q) = sl[q];
              if (xr && r == rk) jt(m,l,16*nkd+8*ur+q) = sr[q];
            }
            if (r > ns) return;
            // to the right: the left fold from tile ul (or the one from the left
            // extended), to the left: the right fold to tile nk - 1 - ul (or extended)
            const int u0 = xl ? 0 : ul, u1 = xr ? nkd - 1 : nkd - 1 - ul;
            for (int q = 0; q < 8; ++q) f[q] = xl ? sl[q] : jt(m,l,8*u0+q);
            for (int u = xl ? 0 : u0 + 1; u < nkd; ++u) {
              for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q);
              FluxRowJoin(f, s, h, ya, xb);
              for (int q = 0; q < 8; ++q) f[q] = h[q];
            }
            for (int q = 0; q < 8; ++q) slot(n - 1, -1, q) = f[q];
            for (int q = 0; q < 8; ++q) f[q] = xr ? sr[q] : jt(m,l,8*u1+q);
            for (int u = xr ? nkd - 1 : u1 - 1; u >= 0; --u) {
              for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q);
              FluxRowJoin(s, f, h, ya, xb);
              for (int q = 0; q < 8; ++q) f[q] = h[q];
            }
            for (int q = 0; q < 8; ++q) slot(0, 1, q) = f[q];
          });
        }
        if (r <= nrt) {   // as FluxSchur
          std::array<uint32_t, 4> fa = {0, 0, 0, 0}, fb = fa;
          for (int d = 0; d < x.nf; ++d) {
            if (gr[d] == g) fa[d] = (gs[d] > 0) ? 3 : 0, fb[d] = (gs[d] > 0) ? 8 : 7;
          }
          if (nd == 2) GhostDefer(cln_w, CLW_BW, 2), cln_hpend.back().fm = fa;
          GhostDefer(cln_w, CLW_TM, (nd == 2) ? 4 : 3);
          cln_hpend.back().fm = fb;
          GhostFlush();
        }
      }
    }
  }
  // each tile's whole line, its ends
  for (int dir = 0; dir < x.nf; ++dir) {
    if (nk[dir] == 1 && gr[dir] < 0) continue;
    const bool pd = per[dir];
    const int nkd = nk[dir], n = nn[dir], nbd = nbk[dir], na = (dir == 0) ? nn[1] : nn[0];
    const int nb = (dir == 2) ? nn[1] : nn[2];
    auto jt = cln_fsjt[dir];
    AthFor("cln_fsjf", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
      auto cell = [&](const int t, int &k, int &j, int &i) {
        k = (dir == 2) ? x.ks + t : x.ks + b;
        j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
        i = (dir == 0) ? x.is + t : x.is + a;
      };
      int k, j, i, kf, jf, iff;
      cell(0, k, j, i);
      const bool bl = IsBoundaryFace(x, dir, m, k, j, i);
      cell(n, k, j, i);
      const bool bh = IsBoundaryFace(x, dir, m, k, j, i);
      for (int u = 0; u < nkd; ++u) {
        Real lf[8], rf[8], s[8], h[8], c[8], cr[8], ya[3], xb[3], yr[3], xr[3];
        Real yd[3], xd[3];
        bool hl = false, hr = false;
        auto addl = [&](const int o) {   // L = join(L, the segment at o)
          for (int q = 0; q < 8; ++q) s[q] = jt(m,l,o+q);
          if (hl) FluxRowJoin(lf, s, h, yd, xd);
          for (int q = 0; q < 8; ++q) lf[q] = hl ? h[q] : s[q];
          hl = true;
        };
        auto addr = [&](const int o) {   // R = join(the segment at o, R)
          for (int q = 0; q < 8; ++q) s[q] = jt(m,l,o+q);
          if (hr) FluxRowJoin(s, rf, h, yd, xd);
          for (int q = 0; q < 8; ++q) rf[q] = hr ? h[q] : s[q];
          hr = true;
        };
        if (pd && nbd == 1) {   // round the wrap of one MeshBlock
          for (int t = u - nkd/2; t < u; ++t) addl(8*((t + nkd)%nkd));
          for (int t = u + (nkd - 1)/2; t > u; --t) addr(8*(t%nkd));
        } else {
          if (pd || !bl) addl(8*nkd + 8*(pd ? u : 0));
          for (int t = 0; t < u; ++t) addl(8*t);
          if (pd ? (u > 0 || (nbd*nkd - 1)/2 >= nkd) : !bh) addr(16*nkd + 8*(pd ? u : 0));
          for (int t = nkd - 1; t > u; --t) addr(8*t);
        }
        for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q), c[q] = s[q];
        if (hl) FluxRowJoin(lf, s, c, ya, xb);
        for (int q = 0; q < 8; ++q) cr[q] = c[q];
        if (hr) FluxRowJoin(c, rf, cr, yr, xr);
        Real xo = 0.0, yo = 0.0;
        if (pd) {
          const Real det = cr[3]*(cr[6] + cr[7]) + cr[1]*cr[7];
          xo = (cr[0]*(cr[6] + cr[7]) + cr[1]*cr[4])/det;
          yo = (cr[4]*(cr[1] + cr[3]) + cr[6]*cr[0])/det;
        }
        const Real gR = hr ? xr[0] + xr[1]*yo + xr[2]*xo : xo;
        const Real gL = hl ? ya[0] + ya[1]*yo + ya[2]*gR : yo;
        if (u > 0) jt(m,l,24*nkd+u) = gL;
        if (u < nkd - 1) jt(m,l,25*nkd+u) = gR;
        Real sg;
        cell(n, k, j, i);
        if (u == nkd - 1 && !bh) {   // cell n, as FluxSchur
          const bool on = ln.Own(x, dir, m, k, j, i, n, n, true);
          Real &f = ln.Slot(x, dir, m, k, j, i, k, j, i, on, sg);
          f = sg*gR;
        }
        cell(0, kf, jf, iff);
        cell(-1, k, j, i);
        if (u == 0 && !bl) {
          const bool on = ln.Own(x, dir, m, kf, jf, iff, 0, n, true);
          Real &f = ln.Slot(x, dir, m, kf, jf, iff, k, j, i, on, sg);
          f = sg*gL;
        }
      }
    });
  }
}

}  // namespace radvet
