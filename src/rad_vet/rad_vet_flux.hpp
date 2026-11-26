#ifndef RAD_VET_RAD_VET_FLUX_HPP_
#define RAD_VET_RAD_VET_FLUX_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_flux.hpp
//! \brief rows of the face flux solve, shared by the rad_vet_flux*.cpp files.

#include <array>
#include <cstdint>

#include "athena.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \struct FluxLine
//! \brief the F solve's lines on every mesh (FluxSolve): the cell terms t of cell
//! (k,j,i) at its face (kf,jf,iff) normal to dir (al, lagged pressure p, w, h+, h-, area
//! pressure pa) and the F of the line's faces.  Spherical polar, Cartesian: t in
//! CLW_L + 3 dir, F of each cell's low face in CLW_Y + dir (CLW_S).  Cubed sphere: cln_vf
//! in the owner's axes (an x2 line may continue as an x3 line across a seam, so a cell
//! keeps both faces of each family; family q: terms 3 q + 0..2, 15 + q, F 9 + 2 q + hi)

struct FluxLine {
  DvceArray5D<Real> v;
  bool cub;
  KOKKOS_INLINE_FUNCTION
  static int Cy(const int dir) {return (dir == 2) ? CLW_S : CLW_Y + dir;}
  KOKKOS_INLINE_FUNCTION
  void Cell(const SolveArgs &x, const int dir, const int m, const int kf, const int jf,
            const int iff, const int k, const int j, const int i, Real t[6]) const {
    if (!cub) {
      for (int n = 0; n < 3; ++n) t[n] = x.w(m,CLW_L+3*dir+n,k,j,i);
      t[3] = x.cs(m,CLN_HP+dir,k,j,i);
      t[4] = x.cs(m,CLN_HM+dir,k,j,i);
      t[5] = 0.0;
      return;
    }
    Real s;
    const int q = CsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2;
    for (int n = 0; n < 3; ++n) t[n] = v(m,3*q+n,k,j,i);
    const Real hp = x.cs(m,CLN_HP+q,k,j,i), hm = x.cs(m,CLN_HM+q,k,j,i);
    t[3] = (s > 0.0) ? hp : -hm;
    t[4] = (s > 0.0) ? hm : -hp;
    t[5] = v(m,15+q,k,j,i);
  }
  // slot of the F of face (kf,jf,iff) (far: of the family's other face of cell (k,j,i))
  // in that cell's storage, s its orientation along ours
  KOKKOS_INLINE_FUNCTION
  Real &Slot(const SolveArgs &x, const int dir, const int m, const int kf, const int jf,
             const int iff, const int k, const int j, const int i, const bool far,
             Real &s) const {
    s = 1.0;
    if (cub) return v(m,9+CsSlot(x, dir, m, kf, jf, iff, k, j, i, far, s),k,j,i);
    if (!far) return x.w(m,Cy(dir),kf,jf,iff);
    const int up = (k == kf && j == jf && i == iff) ? 1 : 0;
    return x.w(m,Cy(dir),k+up*(dir == 2),j+up*(dir == 1),i+up*(dir == 0));
  }
  // F of the family's other face of cell (k,j,i) along face (kf,jf,iff), in our axes
  // (sp, Cartesian: a physical face's F^H); Near: of that face
  KOKKOS_INLINE_FUNCTION
  Real Far(const SolveArgs &x, const int dir, const int m, const int kf, const int jf,
           const int iff, const int k, const int j, const int i) const {
    if (!cub) {
      const int up = (k == kf && j == jf && i == iff) ? 1 : 0;
      const int ko = k + up*(dir == 2), jo = j + up*(dir == 1), io = i + up*(dir == 0);
      if (IsBoundaryFace(x, dir, m, ko, jo, io)) return x.fc[dir](m,CLF_N,ko,jo,io);
    }
    Real s;
    const Real f = Slot(x, dir, m, kf, jf, iff, k, j, i, true, s);
    return s*f;
  }
  KOKKOS_INLINE_FUNCTION
  Real Near(const SolveArgs &x, const int dir, const int m, const int kf, const int jf,
            const int iff, const int k, const int j, const int i) const {
    Real s;
    const Real f = Slot(x, dir, m, kf, jf, iff, k, j, i, false, s);
    return s*f;
  }
  // whether the MeshBlock solves face (k,j,i), number t of the n + 1 of its line (tie:
  // the block whose low face it is; cubed sphere with clean_stag_fschur: the lower
  // gid's), else its identity row
  KOKKOS_INLINE_FUNCTION
  bool Own(const SolveArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t, const int n, const bool tie) const {
    if (!cub) return !(tie && t == n && !IsBoundaryFace(x, dir, m, k, j, i));
    if (!tie || (t > 0 && t < n) || IsBoundaryFace(x, dir, m, k, j, i)) return true;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0), o = (t == 0) ? 1 : 0;
    const Real gc = x.g(m,CLG_F+7,k-(1-o)*dk,j-(1-o)*dj,i-(1-o)*di);
    const Real go = x.g(m,CLG_F+7,k-o*dk,j-o*dj,i-o*di);
    return gc < go || (gc == go && t == 0);
  }
  // F of face (k,j,i), number t of the n + 1 faces of its line along dir: get, put
  KOKKOS_INLINE_FUNCTION
  Real Get(const SolveArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t) const {
    if (!cub) return x.w(m,Cy(dir),k,j,i);
    return (t > 0) ? v(m,9+2*dir+1,k-(dir == 2),j-(dir == 1),i-(dir == 0))
                 : v(m,9+2*dir,k,j,i);
  }
  KOKKOS_INLINE_FUNCTION
  void Put(const SolveArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t, const int n, const Real f) const {
    if (!cub) {
      x.w(m,Cy(dir),k,j,i) = f;
      return;
    }
    if (t > 0) v(m,9+2*dir+1,k-(dir == 2),j-(dir == 1),i-(dir == 0)) = f;
    if (t < n) v(m,9+2*dir,k,j,i) = f;
  }
};

//----------------------------------------------------------------------------------------
//! \fn bool FluxRow
//! \brief the F solve's row lo F_{f-1} + dg F_f + hi F_{f+1} = r of face (k,j,i) normal
//! to dir; false where the face has no equation (F = F^H) or another block owns it (own
//! false): the identity row.  sm (clean_stag_fschur): the row sum, without cancellation.
//! The kinetic stress over h c is q (P_R - P_L) + qs (P_R + P_L): spherical polar's area
//! form (1/dx +- s), Cartesian 1/dx, cubed sphere (A_R P_R - A_L P_L)/V_f.

template <class FC>
KOKKOS_INLINE_FUNCTION
bool FluxRow(const SolveArgs &x, const FluxLine &ln, const FC &fc, const int dir,
             const int m, const int k, const int j, const int i, Real &lo, Real &hi,
             Real &dg, Real &r, const bool own = true, Real *sm = nullptr) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  lo = 0.0, hi = 0.0, dg = 1.0, r = fc(m,CLF_N,k,j,i);
  if (sm) *sm = 1.0;
  if (!own) return false;
  const bool bd = IsBoundaryFace(x, dir, m, k, j, i);
  Real tl[6], tr[6];
  ln.Cell(x, dir, m, k, j, i, k-dk, j-dj, i-di, tl);
  ln.Cell(x, dir, m, k, j, i, k, j, i, tr);
  const bool eq = !bd && (tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0);
  if (eq) {
    const Real g = 0.5*x.h*x.c/CellFaceDistance(x, dir, m, k, j, i);
    const Real it_ = 1.0/fc(m,CLF_TH,k,j,i);
    // fkin: h c^2 [w P_kin]_L^R/dx, P_kin = (h+ F_lo + h- F_hi)/c of a cell; the cubed
    // sphere's cells carry their area pressures A P (FaceMet)
    const FaceMet fm = FaceMetric(x, dir, m, k, j, i);
    const Real q = x.h*x.c*(x.cub ? 1.0 : fm.a)/fm.v, qs = x.h*x.c*fm.s;
    const Real ag = x.cub ? fm.a : 1.0;
    lo = -g*tl[0] - q*tl[2]*tl[3] + qs*tl[2]*tl[3];
    hi = -g*tr[0] + q*tr[2]*tr[4] + qs*tr[2]*tr[4];
    const Real sq = q*(tr[2]*(tr[3] + tr[4]) - tl[2]*(tl[3] + tl[4]));
    const Real ss = qs*(tr[2]*(tr[3] + tr[4]) + tl[2]*(tl[3] + tl[4]));
    dg = it_ - lo - hi + sq + ss;
    if (sm) *sm = it_ + sq + ss;
    r = fc(m,CLF_N,k,j,i)*it_ + q*x.c*(ag*(tl[1] - tr[1]) + tl[5] - tr[5]) -
        qs*x.c*(tl[1] + tr[1]) + x.h*x.c*x.c*fc(m,CLF_N+1,k,j,i);
  }
  return eq;
}

//----------------------------------------------------------------------------------------
//! \fn Real FluxRay
//! \brief clean_stag_fwb: the rays' own lab-frame momentum of face (kf,jf,iff): their
//! flux factor (clamped to |g/J| <= 1) times the donor's E (work vector e), plus the
//! upwind advective a E.  Where the donor holds no rays (CLN_J = 0) the thin share's own
//! energy flux c (h+ E_L + h- E_R): a function of E on every face, never the last F, so
//! that the well-balanced rows of a line hold no lagged neighbour

KOKKOS_INLINE_FUNCTION
Real FluxRay(const SolveArgs &x, const int e, const int dir, const int m,
             const int kf, const int jf, const int iff) {
  const int ek = (dir == 2), ej = (dir == 1), ei = (dir == 0);
  const Real g = RayFace(x, dir, m, kf, jf, iff);
  auto ec = [&](const int l) {return fmax(x.w(m,e,kf-l*ek,jf-l*ej,iff-l*ei), 0.0);};
  const int u = (g > 0.0) ? 1 : 0;
  const Real jd = x.cs(m,CLN_J,kf-u*ek,jf-u*ej,iff-u*ei);
  const Real a = 0.5*(FaceNormalPart(x, CLW_A1, dir, m, kf, jf, iff, kf-ek, jf-ej,
                                     iff-ei) +
                      FaceNormalPart(x, CLW_A1, dir, m, kf, jf, iff, kf, jf, iff));
  const Real f = (jd > 0.0) ? x.c*fmin(fmax(g/jd, -1.0), 1.0)*ec(u) :
                 x.c*(x.fc[dir](m,CLF_HP,kf,jf,iff)*ec(1) +
                      x.fc[dir](m,CLF_HM,kf,jf,iff)*ec(0));
  return f + a*ec(a > 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void FluxRowOwn
//! \brief clean_stag_fschur: a MeshBlock's line (m, a, b) along dir (n cells) as segment
//! s of the F line: first own face x = s0 + s1 yL + s2 xR, last own face y = s4 + s5 yL +
//! s6 xR (yL the face before, xR the face after), s3 = 1 - s1 - s2 and s7 = 1 - s5 - s6
//! from the row sums (no cancellation on nearly singular lines); Thomas sweeps up for y,
//! down for x.  layout_invariant_mode = tile: the tile of n cells from cell o, il / ih
//! its ends inside the MeshBlock, owned as between MeshBlocks.

template <class FC>
KOKKOS_INLINE_FUNCTION
void FluxRowOwn(const SolveArgs &x, const FluxLine &ln, const FC &fc, const int dir,
                const int m, const int a, const int b, const int n, Real s[8],
                const int o = 0, const bool il = false, const bool ih = false) {
  auto face = [&](const int t, int &k, int &j, int &i) {
    k = (dir == 2) ? x.ks + o + t : x.ks + b;
    j = (dir == 1) ? x.js + o + t : ((dir == 2) ? x.js + b : x.js + a);
    i = (dir == 0) ? x.is + o + t : x.is + a;
  };
  int k, j, i;
  Real lo, hi, dg, r, sm;
  face(0, k, j, i);
  const int t0 = (il ? !x.cub : ln.Own(x, dir, m, k, j, i, 0, n, true)) ? 0 : 1;
  face(n, k, j, i);
  const bool bh = !ih && IsBoundaryFace(x, dir, m, k, j, i);   // F^H
  const int t1 = (x.cub && !bh && (ih || ln.Own(x, dir, m, k, j, i, n, n, true))) ? n
                                                                                : n - 1;
  const Real fh = bh ? fc(m,CLF_N,k,j,i) : 0.0, wh = bh ? 1.0 : 0.0;
  Real cp = 0.0, dc = 0.0, dl = 1.0, dw = 0.0;   // F_{-1} = yL
  for (int t = t0; t <= t1; ++t) {
    face(t, k, j, i);
    FluxRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, true, &sm);
    const Real den = dg - lo*cp;
    dc = (r - lo*dc)/den;
    dl = -lo*dl/den;
    dw = (sm - lo*dw)/den;
    cp = hi/den;
  }
  s[4] = dc - cp*fh, s[5] = dl, s[6] = bh ? 0.0 : -cp, s[7] = dw - cp*wh;
  Real cq = 0.0, ec = fh, er = 1.0 - wh, ew = wh;   // F_{t1+1} = xR or F^H
  for (int t = t1; t >= t0; --t) {
    face(t, k, j, i);
    FluxRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, true, &sm);
    const Real den = dg - hi*cq;
    ec = (r - hi*ec)/den;
    er = -hi*er/den;
    ew = (sm - hi*ew)/den;
    cq = lo/den;
  }
  s[0] = ec, s[1] = -cq, s[2] = er, s[3] = ew;
}

//----------------------------------------------------------------------------------------
//! \fn void FluxRowJoin
//! \brief clean_stag_fschur: the segment c of two adjacent ones p, q (FluxRowOwn), and
//! p's last face ya and q's first xb in terms of c's outer faces; s3, s7 as sums (no
//! cancellation)

KOKKOS_INLINE_FUNCTION
void FluxRowJoin(const Real p[8], const Real q[8], Real c[8], Real ya[3],
                 Real xb[3]) {
  const Real d = p[5] + p[7] + p[6]*(q[2] + q[3]);
  ya[0] = (p[4] + p[6]*q[0])/d, ya[1] = p[5]/d, ya[2] = p[6]*q[2]/d;
  xb[0] = (q[0] + q[1]*p[4])/d, xb[1] = q[1]*p[5]/d, xb[2] = q[2]/d;
  c[0] = p[0] + p[2]*xb[0], c[1] = p[1] + p[2]*xb[1], c[2] = p[2]*xb[2];
  c[3] = p[3] + p[2]*(q[3] + q[1]*p[7])/d;
  c[4] = q[4] + q[5]*ya[0], c[5] = q[5]*ya[1], c[6] = q[6] + q[5]*ya[2];
  c[7] = q[7] + q[5]*(p[7] + p[6]*q[3])/d;
}

// layout_invariant_mode = tile: per MeshBlock line of nk tiles, the tiles' segments, the
// folds from the left and right (8 nk each) and the tiles' ends (2 nk)
static void FluxTileArray(DvceArray3D<Real> &v, const int nmb, const int nl,
                          const int nk) {
  if (v.extent_int(0) == nmb && v.extent_int(1) == nl && v.extent_int(2) == 26*nk) return;
  v = DvceArray3D<Real>("cln_fsjt", nmb, nl, 26*nk);
}

// the F slots 9 + 2 q + hi of cln_vf the ghosts of a face of family q take
// (GSeg::fm): its own family's; f0 = 2: without the x1 family's
static std::array<uint32_t, 4> FluxMask(const int f0) {
  return (f0 > 0) ? std::array<uint32_t, 4>{0, 3, 014, 0}
                  : std::array<uint32_t, 4>{3, 014, 060, 0};
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_FLUX_HPP_
