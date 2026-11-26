//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_closure_faces.cpp
//! \brief the rays' sums and half ranges on faces.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_ring.cpp, rad_vet_tasks.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_init.cpp, rad_vet_params_closure.cpp,
//!    rays/sc_cart_setup.cpp, rays/sc_sph_walls.cpp

#include <algorithm>
#include <cstdint>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_moments.hpp"
#include "utils/scratch_view.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn int FvMirror
//! \brief leaving direction of d0 .. d0+nd-1 nearest the mirror image n - 2 (n . t) t of
//! incoming d at a wall (sg t outward; -1: d not incoming); u = |n . t|.  Ties within
//! 1e-12 take the first (CPU/GPU independent)

KOKKOS_INLINE_FUNCTION
int FvMirror(const DvceArray2D<Real> &dir_, const int d0, const int nd, const int d,
             const Real t[3], const Real sg, Real &u) {
  const Real nt = dir_(d,0)*t[0] + dir_(d,1)*t[1] + dir_(d,2)*t[2];
  u = fabs(nt);
  int r = -1;
  if (!(sg*nt < 0.0)) return r;
  const Real v[3] = {dir_(d,0) - 2.0*nt*t[0], dir_(d,1) - 2.0*nt*t[1],
                     dir_(d,2) - 2.0*nt*t[2]};
  Real bd = -2.0;
  for (int e = d0; e < d0 + nd; ++e) {
    const Real et = dir_(e,0)*t[0] + dir_(e,1)*t[1] + dir_(e,2)*t[2];
    const Real dt = v[0]*dir_(e,0) + v[1]*dir_(e,1) + v[2]*dir_(e,2);
    if (sg*et > 0.0 && dt > bd + 1.0e-12) {bd = dt; r = e;}
  }
  return r;
}

//----------------------------------------------------------------------------------------
//! \struct RayIntensity
//! \brief sc ray intensities in the mesh layout without a copy: vgd_i, the first ghost
//! layer from gl (side 0: j = js - 1, 1: k = ks - 1; six: also 2, 3 high j, k and 4, 5
//! low, high i)

struct RayIntensity {
  RadiationVet::ScRag vi;
  DvceArray5D<Real> gl;
  int ok, og, ks, js;
  bool six = false;
  int ke = 0, je = 0, is = 0, ie = 0;
  // (the ghost side of (k, j, i) or -1, and its two indices)
  KOKKOS_INLINE_FUNCTION
  int Side(const int k, const int j, const int i, int &r, int &t) const {
    if (j < js || (six && j > je)) {r = k - ks; t = i; return (j < js) ? 0 : 2;}
    if (k < ks || (six && k > ke)) {r = j - js; t = i; return (k < ks) ? 1 : 3;}
    if (six && (i < is || i > ie)) {r = k - ks; t = j - js; return (i < is) ? 4 : 5;}
    return -1;
  }
  KOKKOS_INLINE_FUNCTION
  Real operator()(const int m, const int d, const int k, const int j, const int i) const {
    int r, t;
    const int s = Side(k, j, i, r, t);
    return (s >= 0) ? gl(m,d,s,r,t) : vi(m,d,k+ok,j+og,i);
  }
  KOKKOS_INLINE_FUNCTION
  DirRow At(const int m, const int k, const int j, const int i) const {
    int r, t;
    const int s = Side(k, j, i, r, t);
    if (s >= 0) return DirRow{&gl(m,0,s,r,t), static_cast<int64_t>(gl.stride(1))};
    return ScDirRow(vi, m, k+ok, j+og, i);
  }
};

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FaceHalf
//! \brief cubed-sphere sc rays: each face's half-range sums J_L, J_R, H+_L, H-_R about
//! its normal (cln_fh), accumulated set by set after each sweep from vgd_i (or, for a
//! version-5 restart, cln_fvi from row fo); union rows o of dv; zero: first set

void RadiationVet::FaceHalf(const DvceArray2D<Real> &dv, const int o,
                            const bool zero, const int fo) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, n = vgd_n, nmb1 = pmy_pack->nmb_thispack - 1;
  if (cln_fh.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(cln_fh, nmb1 + 1, 12, cln_s.extent_int(2), cln_s.extent_int(3),
                    cln_s.extent_int(4));
  }
  auto fh = cln_fh;
  auto fi = cln_fvi;
  auto gn = cln_gn;
  // (the swept set in vgd_i plus ghost layers; a version-5 restart's in cln_fvi)
  const bool fv = cln_fvi.extent_int(0) > 0;
  const RayIntensity gi{vgd_i, cln_gdgl, vgd_w - ks, vgd_w - indcs.ng, ks, js, true, ke,
                        je, is, ie};
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    SolveForK("cln_facehalf", nmb1, ks, ke+dk, js, je+dj, is, ie+di,   // (HIP: k fastest)
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      // (J_L, J_R, H+_L, H-_R as scalars: nvcc -O3 miscompiles a runtime-indexed local)
      Real jl0 = zero ? 0.0 : fh(m,4*dir,k,j,i), jr0 = zero ? 0.0 : fh(m,4*dir+1,k,j,i);
      Real hp0 = zero ? 0.0 : fh(m,4*dir+2,k,j,i), hm0 = zero ? 0.0 : fh(m,4*dir+3,k,j,i);
      // (the two cells' intensity rows once)
      const DirRow gl_ = fv ? DirRow{nullptr, 0} : gi.At(m, kl, jl, il);
      const DirRow gr_ = fv ? DirRow{nullptr, 0} : gi.At(m, k, j, i);
#pragma unroll 4   // (several rays in flight on hipcc)
      for (int e = 0; e < n; ++e) {
        const int d = o + e;
        Real mu = 0.0;
        for (int b = 0; b < 3; ++b) mu += dv(d,b)*gn.At<CLG_N>(m,3*dir+b,k,j,i);
        const Real il_ = fv ? fi(m,fo+e,kl,jl,il) : gl_[e];
        const Real ir_ = fv ? fi(m,fo+e,k,j,i) : gr_[e];
        const Real wl = dv(d,3)*il_, wr = dv(d,3)*ir_;
        jl0 += wl;
        jr0 += wr;
        if (mu > 0.0) hp0 += wl*mu;
        if (!(mu > 0.0)) hm0 += wr*mu;
      }
      fh(m,4*dir,k,j,i) = jl0;
      fh(m,4*dir+1,k,j,i) = jr0;
      fh(m,4*dir+2,k,j,i) = hp0;
      fh(m,4*dir+3,k,j,i) = hm0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RayFaceSums
//! \brief vet_gd face fluxes (CLN_FX, F2, F3, F1): sum_d w (n_d . A_f) I_f / A_f and the
//! stress columns in the face basis, A_f the exact closed vector area.  Lateral faces
//! (and staggered r faces) carry the upwind intensity.  Else r faces follow the cell
//! balances from the inner face (outer = inner + V chi (S - J, -H) - lateral outflow),
//! so each cell's energy balance is the vet_gd exchange to round-off; CLN_G is the stress
//! curvature left over.  Several MeshBlocks along x1: one level per pass from the one
//! below; a reflecting top spreads the remainder over the outermost level.  Cubed
//! sphere: Cartesian momentum, upwind r faces

void RadiationVet::RayFaceSums() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto fi_ = cln_fvi;
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  auto x3f = pmy_pack->pcoord->xx3f;
  auto x3v = pmy_pack->pcoord->x3v;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  const bool wint = vgd_wint;
  auto cs = cln_s;
  const bool f1 = cln_s.extent_int(1) > CLN_F1;   // (the cell layout's CLN_F1)
  SolveArgs x = MakeArgs();
  // (args in device memory: light launches)
  cln_gdx = DevArgs(x, 3);
  const SolveArgs *px = cln_gdx;
  // nx3 = 1 (axisymmetric): no phi boundary; the low phi face's upwind cell is the cell
  // itself through the wall map, the high face carries the low one's sums in its basis
  const bool axi = !pmy_pack->pmesh->three_d && x.sph;
  auto bnd = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i) {
    const SolveArgs &x = *px;
    return !(axi && dir == 2) && IsBoundaryFace(x, dir, m, k, j, i);
  };
  // inner face: no flux if reflecting, else outward beams from the emission (vet_bc_inner
  // emit) or the first cell, inward beams from the first cell
  const bool open = (ibc_x1min != VET_IBC_REFLECT);
  const bool emi = VetBcEmit();
  const Real iemit = VetBathLo();
  const Real ihi = (ibc_x1max == VET_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool rwsc[2] = {!open, ibc_x1max == VET_IBC_REFLECT};
  auto bmi = ibm_i;
  const int bmn = ibm_nr, bmr[4] = {ibm_ray[0], ibm_ray[1], ibm_ray[2], ibm_ray[3]};
  // intensities in the mesh layout with ghosts (re-indexed by the wall map at a periodic
  // seam); clean_rot_ring: two slots (sets blk mod 2) with cln_scmem
  const int ns2 = (cln_ring && cln_scmem) ? 2 : 1;
  const int doff = (ns2 == 1) ? 0 : cln_rqm*n;
  // gl: one radial MeshBlock (sp: no reflecting top): vgd_i itself plus exchanged first
  // ghost layers, low theta, phi only (cubed sphere, six: all faces), in direction chunks
  const bool six = x.cub && vgd_nb1 == 1 && !cln_scmem;
  const bool gl = (x.sph && !axi && pmy_pack->pmesh->three_d &&
                   vgd_nb1 == 1 && !cln_scmem) || six;
  RayIntensity gi{vi_, DvceArray5D<Real>(), ok, og, ks, js, six, ke, je, is, ie};
  if (gl) {
    const int ng = indcs.ng;
    // (cubed sphere: chunks of clean_halo_vars/4, the build's largest scratch)
    const int nc = std::min(n, six ? std::max(cln_hmax/4, 1) : cln_hmax);
    const int nt = std::max(ke - ks + 1, je - js + 1);
    // (persistent scratch; gl's unwritten rows stay zero, ch is written before read)
    cln_gdgl = ScratchView<DvceArray5D<Real>>("cln_gd_gl", nmb1 + 1, n, six ? 6 : 2, nt,
                                              std::max(cs.extent_int(4), je - js + 1));
    gi.gl = cln_gdgl;
    // (over the ring's Delta, or (six) cln_st off the arena, where it fits: dead here)
    const std::size_t nch = static_cast<std::size_t>(nmb1 + 1)*nc*cs.extent_int(2)*
                            cs.extent_int(3)*cs.extent_int(4);
    const bool orp = cln_rpfree && cln_rpd.span() >= nch;
    const bool ost = six && cln_st.span() >= nch && !MemArenaIn(cln_st.data());
    auto ch = (orp || ost) ? DvceArray5D<Real>(orp ? cln_rpd.data() : cln_st.data(),
                             nmb1 + 1, nc, cs.extent_int(2), cs.extent_int(3),
                             cs.extent_int(4))
            : (six || cln_ring == 1) ?
              DvceArray5D<Real>("cln_gd_ch", nmb1 + 1, nc, cs.extent_int(2),
                                cs.extent_int(3), cs.extent_int(4))
                  : ScratchView<DvceArray5D<Real>>("cln_gd_ch", nmb1 + 1, nc,
                                                   cs.extent_int(2), cs.extent_int(3),
                                                   cs.extent_int(4));
    auto gg = gi.gl;
    for (int d0 = 0; d0 < n; d0 += nc) {
      const int nd = std::min(nc, n - d0);
      // (only the ng layers the exchange reads)
      SolveFor("cln_gd_chc", nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (j >= js + ng && j <= je - ng && k >= ks + ng && k <= ke - ng &&
            !(six && (i < is + ng || i > ie - ng))) {
          return;
        }
        const DirRow r = ScDirRow(vi_, m, k+ok, j+og, i);
        for (int c = 0; c < nd; ++c) ch(m,c,k,j,i) = r[d0+c];
      });
      GhostFill(ch, 0, nd, true, 0, 1);   // (glc, glr read the first layer)
      par_for("cln_gd_glc", DevExeSpace(), 0, nmb1, 0, nd - 1, 0, nt - 1, is, ie,
      KOKKOS_LAMBDA(const int m, const int c, const int t, const int i) {
        if (t <= ke - ks) gg(m,d0+c,0,t,i) = ch(m,c,ks+t,js-1,i);
        if (t <= je - js) gg(m,d0+c,1,t,i) = ch(m,c,ks-1,js+t,i);
        if (six && t <= ke - ks) gg(m,d0+c,2,t,i) = ch(m,c,ks+t,je+1,i);
        if (six && t <= je - js) gg(m,d0+c,3,t,i) = ch(m,c,ke+1,js+t,i);
      });
      if (!six) continue;
      par_for("cln_gd_glr", DevExeSpace(), 0, nmb1, 0, nd - 1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
        gg(m,d0+c,4,k-ks,j-js) = ch(m,c,k,j,is-1);
        gg(m,d0+c,5,k-ks,j-js) = ch(m,c,k,j,ie+1);
      });
    }
  } else {
    if (cln_fvi.extent_int(1) != ns2*n) {
      Kokkos::realloc(cln_fvi, nmb1 + 1, ns2*n, cs.extent_int(2),
                      cs.extent_int(3), cs.extent_int(4));
    }
    fi_ = cln_fvi;
    SolveFor("cln_gd_copy", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int d = 0; d < n; ++d) fi_(m,doff+d,k,j,i) = vi_(m,d,k+ok,j+og,i);
    });
    GhostFill(cln_fvi, doff, n);
  }
  auto inten = [=] KOKKOS_FUNCTION (const int m, const int d, const int k, const int j,
                                    const int i) {
    return gl ? gi(m,d,k,j,i) : fi_(m,doff+d,k,j,i);
  };
  // vector area of the low face normal to dir and the local basis there
  auto face = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                   const int i, Real a[3], Real e[3][3]) {
    const SolveArgs &x = *px;
    FaceVector(x, x3f, dir, m, k, j, i, a);
    if (x.sph) {
      LocalBasis((dir == 1) ? x.x2f(m,j) : x.x2v(m,j),
                 (dir == 2) ? x3f(m,k) : x3v(m,k), e);
      return;
    }
    for (int b = 0; b < 3; ++b) {
      for (int c = 0; c < 3; ++c) e[b][c] = (b == c) ? 1.0 : 0.0;
    }
  };
  // beam of direction d through the low lateral face of vector area a (not a boundary);
  // gl: ro, ru the rows of the cell and its lower neighbour
  auto lat = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                  const int i, const int d, const Real a[3],
                                  const DirRow &ro, const DirRow &ru) {
    const SolveArgs &x = *px;
    const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
    const int kl = k - (dir == 2), jl = j - (dir == 1);
    if (!(f > 0.0)) return f*(gl ? ro[d] : inten(m, d, k, j, i));
    if (axi && dir == 2) return f*inten(m, mp_(m,kl+ok,j+og,d), k, j, i);
    const bool seam = (dir == 1) ? (j == js && x.mbc(m,2) == BoundaryFlag::periodic)
                                 : (k == ks && x.mbc(m,4) == BoundaryFlag::periodic);
    if (!seam || wl_(m,kl+ok,jl+og) == 0) return f*(gl ? ru[d] : inten(m, d, kl, jl, i));
    if (!wint) {
      const int e = mp_(m,kl+ok,jl+og,d);
      return f*(gl ? ru[e] : inten(m, e, kl, jl, i));
    }
    Real v = 0.0;
    for (int t = 0; t < 3; ++t) {
      v += w3_(m,kl+ok,jl+og,d,t)*inten(m, m3_(m,kl+ok,jl+og,d,t), kl, jl, i);
    }
    return f*v;
  };
  // vet_gd_wall_sym / vet_gd_theta_seam_sym: a seam face's sums q averaged with those
  // seen from the other side (own beams exact, incoming via the mirrored wall map)
  const bool wsx = vgd_wsx, wst = vgd_wst;
  const Real wsc = cos(vgd_wsw), wss = sin(vgd_wsw);
  auto mi_ = vgd_mir;
  auto mz_ = vgd_mirz;
  auto seam2 = [=] KOKKOS_FUNCTION (const int dir, const int m, const int k, const int j,
                                    const int i, const Real a[3], const Real e[3][3],
                                    const Real ar, Real q[4]) {
    const SolveArgs &x = *px;
    const int p = (dir == 2) ? 1 : 0;
    if (!(p ? wsx : wst) || (p ? (k != ks) : (j != js)) ||
        x.mbc(m, p ? 4 : 2) != BoundaryFlag::periodic ||
        wl_(m,k-p+ok,j-1+p+og) == 0) {return;}
    // the other side's face vector b and basis eb: phi, rotated about z by W; theta,
    // Z = diag(1, 1, -1): b = -Z a, eb = (Z e_r, -Z e_theta, Z e_phi)
    Real b[3], eb[3][3], s[4] = {0.0}, qb[4];
    for (int c = 0; c < 3; ++c) {
      const Real zc = (c == 2) ? -1.0 : 1.0;
      b[c] = p ? ((c == 0) ? wsc*a[0] - wss*a[1] : (c == 1) ? wss*a[0] + wsc*a[1] : a[2])
               : -zc*a[c];
      for (int t = 0; t < 3; ++t) {
        eb[t][c] = p ? ((c == 0) ? wsc*e[t][0] - wss*e[t][1]
                                 : (c == 1) ? wss*e[t][0] + wsc*e[t][1] : e[t][2])
                     : ((t == 1) ? -zc : zc)*e[t][c];
      }
    }
    const int kl = k - p, jl = j - 1 + p;
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*b[0] + dir_(d,1)*b[1] + dir_(d,2)*b[2];
      const int dm = p ? mi_(d) : mz_(d);
      Real v = (f > 0.0) ? inten(m, d, kl, jl, i) : 0.0;
      for (int t = 0; t < 3 && !(f > 0.0); ++t) {
        const int r = m3_(m,kl+ok,jl+og,dm,t);
        v += w3_(m,kl+ok,jl+og,dm,t)*inten(m, p ? mi_(r) : mz_(r), k, j, i);
      }
      const Real w = dir_(d,3)*f*v;
      s[0] += w;
      for (int t = 0; t < 3; ++t) s[1+t] += w*dir_(d,t);
    }
    for (int t = 0; t < 4; ++t) {
      qb[t] = (t == 0) ? s[0]/ar
                       : (eb[t-1][0]*s[1] + eb[t-1][1]*s[2] + eb[t-1][2]*s[3])/ar;
      q[t] = 0.5*(q[t] + qb[t]);
    }
  };
  // beam sums through the face per unit area (basis e, area ar)
  auto add = [=] KOKKOS_FUNCTION (const Real e[3][3], const Real ar, const int d,
                                  const Real f, Real q[4]) {
    const Real w = dir_(d,3)*f/ar;
    q[0] += w;
    for (int b = 0; b < 3; ++b) {
      q[1+b] += w*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
    }
  };
  // the same with Cartesian momentum, projected and divided by ar once per face (proj)
  auto addc = [=] KOKKOS_FUNCTION (const int d, const Real f, Real s[4]) {
    const Real w = dir_(d,3)*f;
    s[0] += w;
    for (int c = 0; c < 3; ++c) s[1+c] += w*dir_(d,c);
  };
  auto proj = [=] KOKKOS_FUNCTION (const Real e[3][3], const Real ar, const Real s[4],
                                   Real q[4]) {
    q[0] = s[0]/ar;
    for (int b = 0; b < 3; ++b) q[1+b] = (e[b][0]*s[1] + e[b][1]*s[2] + e[b][2]*s[3])/ar;
  };
  // a mirror r face, per unit area: no energy, twice the normal momentum of cell ic's
  // beams that hit it (sign sg), no tangential part
  auto mirror = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                     const int ic, const Real sg, Real q[4]) {
    const SolveArgs &x = *px;
    Real a[3], e[3][3], w[4] = {0.0}, nf[3], wn = 0.0;
    face(0, m, k, j, i, a, e);
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      const Real tr = inten(m, d, k, j, ic);
      if (sg*f > 0.0) add(e, x.area[0](m,k,j,i), d, 2.0*f*tr, w);
    }
    for (int b = 0; b < 3; ++b) nf[b] = NumFaceDirs(x, b, 0, m, k, j, i);
    for (int b = 0; b < 3; ++b) wn += w[1+b]*nf[b];
    q[0] = 0.0;
    for (int b = 0; b < 3; ++b) q[1+b] = wn*nf[b];
  };
  // with gl: a cell's low lateral faces and low r face in one pass over the directions
  // (cln_gd_fs), each sum in the same order as the separate kernels; what few faces add
  // (r boundaries, cubed-sphere high faces, seams) in kernels of their own
  const bool fs = gl && !axi;
  // a lateral face of cell (k,j,i): the low (up 0), or cubed sphere the MeshBlock's high
  // one (CLN_FH of its last cell)
  auto latf = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                   const int dir, const int up) {
    const SolveArgs &x = *px;
    const int kf = k + up*(dir == 2), jf = j + up*(dir == 1);
    Real q[4], sc[4] = {0.0}, a[3], e[3][3];
    face(dir, m, kf, jf, i, a, e);
    const bool bd = bnd(dir, m, kf, jf, i);
    const Real ar = x.area[dir](m,kf,jf,i);
    const DirRow ro = gl ? gi.At(m, kf, jf, i) : DirRow{nullptr, 0};
    const DirRow ru = (gl && !bd) ? gi.At(m, kf - (dir == 2), jf - (dir == 1), i)
                                    : DirRow{nullptr, 0};
    for (int d = 0; d < n; ++d) {
      addc(d, bd ? 0.0 : lat(dir, m, kf, jf, i, d, a, ro, ru), sc);
    }
    proj(e, ar, sc, q);
    if (!bd && !up) {seam2(dir, m, kf, jf, i, a, e, ar, q);}
    const int nf = up ? CLN_FH + 4*(dir - 1) : ((dir == 1) ? CLN_F2 : CLN_F3);
    for (int c = 0; c < 4; ++c) cs(m,nf+c,k,j,i) = q[c];
  };
  for (int dir = 1; dir < 3 && !fs; ++dir) {
    if (!x.sph && dir >= x.nf) continue;
    SolveForK("cln_gd_lat", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      latf(m, k, j, i, dir, 0);
      if (px->cub && ((dir == 1) ? j == je : k == ke)) latf(m, k, j, i, dir, 1);
    });
  }
  if (!fs) GhostFill(cln_s, CLN_F2, 4);
  if (!fs) GhostFill(cln_s, CLN_F3, 4);
  // rr (spherical polar): r faces carry the rays' own upwind fluxes, not the balance
  // march, which loses the diffusion flux of thick cells in V chi (S - J)
  const bool rr = x.sph;
  auto rflux = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                    Real q[4]) {
    const SolveArgs &x = *px;
    Real a[3], e[3][3], sc[4] = {0.0};
    face(0, m, k, j, i, a, e);
    // (reflecting x1 face: mirror images of the leaving beams, FvMirror, scaled so
    // the face returns the energy it receives)
    const int ic = (i == is) ? is : ie;
    const Real sg = (i == is) ? -1.0 : 1.0, an = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    const Real t[3] = {a[0]/an, a[1]/an, a[2]/an};
    const bool msc = IsBoundaryFace(x, 0, m, k, j, i) && rwsc[i != is];
    Real sm = 1.0, ho = 0.0, hm = 0.0, u;
    for (int d = 0; d < n && msc; ++d) {
      const int r = FvMirror(dir_, 0, n, d, t, sg, u);
      if (r < 0) ho += dir_(d,3)*u*inten(m, d, k, j, ic);
      if (r >= 0) hm += dir_(d,3)*u*inten(m, r, k, j, ic);
    }
    if (ho > 0.0 && hm > 0.0) sm = ho/hm;
    const DirRow rl = gl ? gi.At(m, k, j, i - 1) : DirRow{nullptr, 0};
    const DirRow rh = gl ? gi.At(m, k, j, i) : DirRow{nullptr, 0};
    for (int d = 0; d < n; ++d) {
      const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
      Real iv = gl ? ((f > 0.0) ? rl[d] : rh[d])
                   : inten(m, d, k, j, (f > 0.0) ? i - 1 : i);
      if (IsBoundaryFace(x, 0, m, k, j, i)) {   // the sweep's x1 boundary
        const Real s0 = exp(cw_(m,1,k+ok,j+og,is));
        iv = (i == is) ? ((f > 0.0) ? (emi ? iemit : s0) : iv)
                       : ((f > 0.0) ? iv : ihi);
        for (int q = 0; q < bmn && i == is && f > 0.0; ++q) {
          if (d == bmr[q]) iv += bmi(m,q,(bmi.extent_int(2) > 1) ? k - ks : 0,j - js);
        }
        const int r = msc ? FvMirror(dir_, 0, n, d, t, sg, u) : -1;
        if (r >= 0) iv = sm*inten(m, r, k, j, ic);
      }
      addc(d, f*iv, sc);
    }
    proj(e, 1.0, sc, q);   // (the face's sums, not per area)
  };
  auto xf_ = vgd_x1f;
  auto qb_ = vgd_x1b;
  auto qr_ = vgd_x1r;
  const int nj = je - js + 1, nkj = (ke - ks + 1)*nj;
  // the rays' r face fluxes rf, all faces at once (the columns only read them); cubed
  // sphere: above the lowest face, no march there
  const bool rf_on = rr || x.cub;
  const int nk = ke - ks + 1, i1 = x.cub ? 1 : 0;
  DvceArray5D<Real> rf = ScratchView<DvceArray5D<Real>>("cln_gd_rf",
                       rf_on ? nmb1 + 1 : 1, 4, rf_on ? nk : 1, rf_on ? nj : 1,
                       rf_on ? ie - is + 2 : 1);
  Kokkos::deep_copy(DevExeSpace(), rf, 0.0);
  if (fs) {
    SolveForK("cln_gd_fs", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const SolveArgs &x = *px;
      const bool rb = (i == is);
      Real a1[3], a2[3], a0[3] = {0.0}, e[3][3], q[4];
      Real s1[4] = {0.0}, s2[4] = {0.0}, s0[4] = {0.0};
      const bool b1 = bnd(1, m, k, j, i), b2 = bnd(2, m, k, j, i);
      const DirRow ro = gi.At(m, k, j, i);
      const DirRow u1 = b1 ? DirRow{nullptr, 0} : gi.At(m, k, j - 1, i);
      const DirRow u2 = b2 ? DirRow{nullptr, 0} : gi.At(m, k - 1, j, i);
      const DirRow rl = rb ? DirRow{nullptr, 0} : gi.At(m, k, j, i - 1);
      FaceVector(x, x3f, 1, m, k, j, i, a1);
      FaceVector(x, x3f, 2, m, k, j, i, a2);
      if (!rb) FaceVector(x, x3f, 0, m, k, j, i, a0);
      for (int d = 0; d < n; ++d) {
        addc(d, b1 ? 0.0 : lat(1, m, k, j, i, d, a1, ro, u1), s1);
        addc(d, b2 ? 0.0 : lat(2, m, k, j, i, d, a2, ro, u2), s2);
        if (!rb) {
          const Real f = dir_(d,0)*a0[0] + dir_(d,1)*a0[1] + dir_(d,2)*a0[2];
          addc(d, f*((f > 0.0) ? rl[d] : ro[d]), s0);
        }
      }
      face(1, m, k, j, i, a1, e);
      proj(e, x.area[1](m,k,j,i), s1, q);
      for (int c = 0; c < 4; ++c) cs(m,CLN_F2+c,k,j,i) = q[c];
      face(2, m, k, j, i, a2, e);
      proj(e, x.area[2](m,k,j,i), s2, q);
      for (int c = 0; c < 4; ++c) cs(m,CLN_F3+c,k,j,i) = q[c];
      if (rb) return;
      face(0, m, k, j, i, a0, e);
      proj(e, 1.0, s0, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
    for (int dir = 1; dir < 3; ++dir) {
      const bool kp = (dir == 2), sm = kp ? wsx : wst;   // (a seam's sums, both sides)
      for (int up = 0; up < 2; ++up) {
        if (up ? !x.cub : !sm) continue;
        par_for("cln_gd_lats", DevExeSpace(), 0, nmb1, 0, (kp ? nj : nk) - 1, is, ie,
        KOKKOS_LAMBDA(const int m, const int t, const int i) {
          const int k = kp ? (up ? ke : ks) : ks + t, j = kp ? js + t : (up ? je : js);
          if (up) {
            latf(m, k, j, i, dir, 1);
          } else if (!bnd(dir, m, k, j, i)) {
            const SolveArgs &x = *px;
            const int nf = kp ? CLN_F3 : CLN_F2;
            Real a[3], e[3][3], q[4];
            for (int c = 0; c < 4; ++c) q[c] = cs(m,nf+c,k,j,i);
            face(dir, m, k, j, i, a, e);
            seam2(dir, m, k, j, i, a, e, x.area[dir](m,k,j,i), q);
            for (int c = 0; c < 4; ++c) cs(m,nf+c,k,j,i) = q[c];
          }
        });
      }
    }
    GhostFill(cln_s, CLN_F2, 4);
    GhostFill(cln_s, CLN_F3, 4);
    par_for("cln_gd_rb", DevExeSpace(), 0, nmb1, i1, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int b, const int k, const int j) {
      const int i = b ? ie + 1 : is;
      Real q[4];
      rflux(m, k, j, i, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
  } else if (rf_on) {
    SolveForK("cln_gd_rf", nmb1, ks, ke, js, je, is + i1, ie + 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real q[4];
      rflux(m, k, j, i, q);
      for (int c = 0; c < 4; ++c) rf(m,c,k-ks,j-js,i-is) = q[c];
    });
  }
  // the lateral faces' part of cell (k,j,i): their sums into q, the curvature into g
  auto latg = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i,
                                   const Real e[3][3], Real q[4], Real g[2]) {
    const SolveArgs &x = *px;
    Real a[3], ef[3][3];
    for (int dir = 1; dir < 3; ++dir) {
      const int nf = (dir == 1) ? CLN_F2 : CLN_F3;
      for (int up = 0; up < 2; ++up) {
        const int kk = k + (up && dir == 2), jj = j + (up && dir == 1);
        if (bnd(dir, m, kk, jj, i)) continue;
        face(dir, m, kk, jj, i, a, ef);
        const Real sg = (up ? -1.0 : 1.0)*x.area[dir](m,kk,jj,i);
        const int kc = (axi && dir == 2) ? k : kk;
        q[0] += sg*cs(m,nf,kc,jj,i);
        for (int b = 0; b < 3; ++b) {
          if (b < 2) g[b] += ((dir == 2 && axi) ? 0.0 : sg*cs(m,nf+1+b,kc,jj,i));
          for (int c = 0; c < 3; ++c) {
            const Real r = sg*cs(m,nf+1+c,kc,jj,i)*(ef[c][0]*e[b][0] +
                           ef[c][1]*e[b][1] + ef[c][2]*e[b][2]);
            q[1+b] += r;
            if (b < 2) g[b] -= r;
          }
        }
      }
    }
  };
  // rr without a reflecting top: no march, one thread per cell; the columns then only
  // close their ends
  const bool cel = rr && !x.cub && ibc_x1max != VET_IBC_REFLECT;
  if (cel) {
    SolveFor("cln_gd_cel", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const SolveArgs &x = *px;
      Real a[3], e[3][3], q[4], g[2] = {0.0, 0.0};
      face(0, m, k, j, is, a, e);
      for (int c = 0; c < 4; ++c) q[c] = rf(m,c,k-ks,j-js,i-is);
      for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
      const Real qi = q[2];
      latg(m, k, j, i, e, q, g);
      const Real q2 = rf(m,2,k-ks,j-js,i+1-is);
      g[1] += q2*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) - qi*(1.0 - x.x1f(m,i)/x.x1v(m,i));
      for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
    });
  }
  for (int lev = 0; lev < vgd_nb1; ++lev) {
    if (lev > 0) ScX1Pass(true);
    SolveForCol("cln_gd_col", nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const SolveArgs &x = *px;
      Real a[3], e[3][3], q[4] = {0.0};
      const bool in = (xf_(m) & 1);
      const int kj = (k - ks)*nj + j - js;
      face(0, m, k, j, is, a, e);
      for (int d = 0; open && in && !rr && d < n; ++d) {   // (rr: the rays' rf below)
        const Real f = dir_(d,0)*a[0] + dir_(d,1)*a[1] + dir_(d,2)*a[2];
        const Real iv = (emi && f > 0.0) ? iemit : inten(m, d, k, j, is);
        add(e, x.area[0](m,k,j,is), d, f*iv, q);
      }
      if (in && !open && !rr) mirror(m, k, j, is, is, -1.0, q);
      for (int c = 0; c < 4; ++c) {
        q[c] = in ? q[c]*x.area[0](m,k,j,is) : qr_(m,c*nkj + kj);
      }
      for (int c = 0; c < 4 && rr; ++c) q[c] = rf(m,c,k-ks,j-js,0);
      const Real vb = in ? 0.0 : qr_(m,4*nkj + kj);   // the column's volume below
      for (int i = is; i <= ie && !cel; ++i) {
        for (int c = 0; c < 4; ++c) cs(m,CLN_FX+c,k,j,i) = q[c]/x.area[0](m,k,j,i);
        if (x.cub) {   // (no curvature terms: Cartesian momentum)
          for (int c = 0; c < 4; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
          continue;
        }
        Real g[2] = {0.0, 0.0}, qi = q[2];
        // the exchange V chi (S - J, -H_b) and the lateral outflow of the cell
        const Real chi = x.vol(m,k,j,i)*exp(cw_(m,0,k+ok,j+og,i));
        const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
        for (int d = 0; d < n && !rr; ++d) {   // (rr: the face fluxes are the rays')
          const Real di = dir_(d,3)*inten(m, d, k, j, i);
          q[0] += chi*(dir_(d,3)*s0 - di);
          for (int b = 0; b < 3; ++b) {
            q[1+b] -= chi*di*(dir_(d,0)*e[b][0] + dir_(d,1)*e[b][1] + dir_(d,2)*e[b][2]);
          }
        }
        latg(m, k, j, i, e, q, g);
        for (int c = 0; c < 4 && rr; ++c) q[c] = rf(m,c,k-ks,j-js,i+1-is);
        // the angular-momentum weight r_f/r of the y_theta row on the r faces
        g[1] += q[2]*(1.0 - x.x1f(m,i+1)/x.x1v(m,i)) - qi*(1.0 - x.x1f(m,i)/x.x1v(m,i));
        for (int b = 0; b < 2; ++b) cs(m,CLN_G+b,k,j,i) = g[b];
      }
      for (int c = 0; c < 4 && cel; ++c) q[c] = rf(m,c,k-ks,j-js,ie+1-is);
      if (xf_(m) & 4) {
        // reflecting top: the remainder at the wall (all but the mirror's momentum)
        // spread over the column by volume, so no energy leaves (lower levels: below)
        Real qw[4];
        mirror(m, k, j, ie + 1, ie, 1.0, qw);
        for (int c = 0; c < 4; ++c) q[c] -= qw[c]*x.area[0](m,k,j,ie+1);
        Real vt = vb, vc = vb;
        for (int i = is; i <= ie; ++i) vt += x.vol(m,k,j,i);
        for (int i = is + (in ? 1 : 0); i <= ie; ++i) {
          if (i > is) vc += x.vol(m,k,j,i-1);
          for (int c = 0; c < 4; ++c) {
            cs(m,CLN_FX+c,k,j,i) -= q[c]*vc/(vt*x.area[0](m,k,j,i));
          }
        }
        for (int c = 0; c < 4; ++c) qb_(m,c*nkj + kj) = q[c];
        qb_(m,4*nkj + kj) = vt;
        qb_(m,5*nkj + kj) = vb;
        for (int c = 0; c < 4; ++c) q[c] = qw[c]*x.area[0](m,k,j,ie+1);
      }
      for (int c = 0; c < 4; ++c) {
        if (f1) cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
        if (f1) cs(m,CLN_F1+c,k,j,ie) = q[c]/x.area[0](m,k,j,ie+1);
        if (!(xf_(m) & 4)) qb_(m,c*nkj + kj) = q[c];
      }
      if (!(xf_(m) & 4)) {
        Real v = vb;
        for (int i = is; i <= ie; ++i) v += x.vol(m,k,j,i);
        qb_(m,4*nkj + kj) = v;
      }
    });
  }
  // reflecting top over several levels: remainder and volumes handed down, each level
  // takes its share
  for (int lev = vgd_nb1 - 2; lev >= 0 && ibc_x1max == VET_IBC_REFLECT; --lev) {
    ScX1Pass(false);
    par_for("cln_gd_top", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const SolveArgs &x = *px;
      if ((xf_(m) >> 3) != lev) return;
      const int kj = (k - ks)*nj + j - js;
      const Real vt = qr_(m,4*nkj + kj);
      Real vc = qr_(m,5*nkj + kj);
      for (int i = is; i <= ie; ++i) vc -= x.vol(m,k,j,i);
      for (int c = 0; c < 6; ++c) qb_(m,c*nkj + kj) = qr_(m,c*nkj + kj);
      qb_(m,5*nkj + kj) = vc;
      for (int i = is; i <= ie + 1; ++i) {
        if (i > is) vc += x.vol(m,k,j,i-1);
        for (int c = 0; c < 4; ++c) {
          const Real dq = qr_(m,c*nkj + kj)*vc/vt;
          if (i <= ie) cs(m,CLN_FX+c,k,j,i) -= dq/x.area[0](m,k,j,i);
          if (i == is && f1) cs(m,CLN_F1+c,k,j,is) = cs(m,CLN_FX+c,k,j,is);
          if (i == ie + 1 && f1) cs(m,CLN_F1+c,k,j,ie) -= dq/x.area[0](m,k,j,ie+1);
        }
      }
    });
  }
}

}  // namespace radvet
