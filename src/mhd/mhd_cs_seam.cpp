//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_cs_seam.cpp
//! \brief Shared cubed-sphere seam treatment for MHD, <mhd>/cs_seam_mhd = shared_o
//! (emf_scheme = uct | uct_hlld, 3-D, uniform mesh):
//! - seam faces: every panel-seam face is one Riemann problem; the ghost-side state is
//!   the neighbour panel's own PLM state of its edge cell (bvals/cs_seam_states.hpp),
//!   including the radial ghost rows is-1 / ie+1 (CsSeamFaceStates).
//! - seam edges: every seam-edge and cube-vertex EMF is the per-face UCT sum
//!     E.t = -sum_f [ a_f^L vn^L X^L + a_f^R vn^R X^R - (d_f^R X^R - d_f^L X^L) ],
//!     X^s = B^s.(t x g_f), g_f = S^-1 n_f, S = sum_f n_f n_f^T (plane normal to t),
//!   (the 4-face UCT rule on a regular edge, the 3-face rule g_f = (2/3) n_f at a cube
//!   vertex), each panel adding only the parts it forms on its own grid; the partial
//!   sums are exchanged and added on both panels (UctSeamOwnSide).
//! - without the shared seam, <mesh>/cs_vertex_3cell: the cube-vertex radial EMF from
//!   the three real cells (CsVertexEmf3).
//! - with a resistivity: the two-chart mean of the resistive seam-edge EMF
//!   (CsSeamResist).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "bvals/bvals.hpp"
#include "bvals/cs_seam_states.hpp"
#include "reconstruct/plm.hpp"
#include "mhd.hpp"
#include "hydro/hydro_wball.hpp"
#include "mhd/emf_smooth.hpp"

namespace mhd {

// speed: every seam kernel captures a SeamCtx (13 Views + 3 DualViews,
// ~1 kB), over Kokkos' 512 B constant-memory threshold, so the default launch went
// through the constant-memory path, which waits on the previous constant-memory launch
// (a host stall per kernel; nsys: +18 ms of cudaEventSynchronize per r1m cycle). The
// same loops launched with HintLightWeight pass the functor as a kernel argument
// (< 4 kB) with no host synchronisation; the arithmetic is unchanged (as BvalsTeamFor).
template <typename Function>
inline void seam_for(const char *name, DevExeSpace, const int kl, const int ku,
                     const int jl, const int ju, const int il, const int iu,
                     const Function &function) {
  static_assert(sizeof(Function) <= 3072, "seam kernel functor too large");
  const int nj = ju - jl + 1, ni = iu - il + 1, nji = nj*ni, nkji = (ku - kl + 1)*nji;
  Kokkos::parallel_for(name, Kokkos::Experimental::require(
      Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nkji),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
  KOKKOS_LAMBDA(const int &idx) {
    int k = idx/nji;
    int j = (idx - k*nji)/ni;
    int i = (idx - k*nji - j*ni) + il;
    function(k + kl, j + jl, i);
  });
}
template <typename Function>
inline void seam_for(const char *name, DevExeSpace, const int nl, const int nu,
                     const int kl, const int ku, const int jl, const int ju,
                     const int il, const int iu, const Function &function) {
  static_assert(sizeof(Function) <= 3072, "seam kernel functor too large");
  const int nk = ku - kl + 1, nj = ju - jl + 1, ni = iu - il + 1;
  const int nji = nj*ni, nkji = nk*nji, nnkji = (nu - nl + 1)*nkji;
  Kokkos::parallel_for(name, Kokkos::Experimental::require(
      Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nnkji),
      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
  KOKKOS_LAMBDA(const int &idx) {
    int n = idx/nkji;
    int k = (idx - n*nkji)/nji;
    int j = (idx - n*nkji - k*nji)/ni;
    int i = (idx - n*nkji - k*nji - j*ni) + il;
    function(n + nl, k + kl, j + jl, i);
  });
}

namespace cs_uct {
// FOUR-BLOCK SEAM EDGES: the per-cell seam-row record X (a regular cell-centred halo
// exchange, one layer; only the within-panel along-seam ghosts are read): per family
// (0 = x2 seams, 1 = x3 seams) at offset 38*fam: 0-2 the seam face's a, dL, dR, 3-5
// its face velocity (Cartesian), 6-11 the cell's along-seam half slopes (HL, HU) set 0,
// 12-17 set 1 (ghost row), 18-27 the cell's lower along-seam face record set 0, 28-37
// set 1. So the block on the other side of a four-block seam edge sees this block's
// shared seam-face store and end records exactly (copies), and both evaluate the edge
// from the same inputs.
// Two arrays (a cell-centred exchange moves at most 32 variables, BvVarList): XA = 0-17,
// XB = the two face records (18-37 above, stored from 0). One set of slots serves both
// families: a cell in a seam row of both is a cube-vertex corner cell, never read here.
constexpr int kNXA = 6;

struct SeamCtx {
  DvceArray5D<Real> u1, u2, u3, w;
  DvceArray4D<Real> b1, b2, b3, e21, e31;
  DvceArray5D<Real> o2, o3, n2, n3;            // own / neighbour packets (x2f, x3f)
  DualArray1D<int> panel;
  DualArray1D<RegionSize> size;
  DualArray2D<NeighborBlock> nghbr;
  int is, ie, js, je, ks, ke, nx2, nx3;
};

KOKKOS_INLINE_FUNCTION Real Dot3(const Real *a, const Real *b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
KOKKOS_INLINE_FUNCTION Real Half(const Real qm, const Real q0, const Real qp,
                                 const bool upper) {
  Real ql, qr;
  PLM(qm, q0, qp, ql, qr);
  return upper ? (ql - q0) : (qr - q0);
}
KOKKOS_INLINE_FUNCTION Real Edge(const Real qm, const Real q0, const Real qp,
                                 const bool upper) {
  Real ql, qr;
  PLM(qm, q0, qp, ql, qr);
  return upper ? ql : qr;
}

// local frame at (xi, eta) of panel p: rhat, unit tangents e1 (xi), e2 (eta), and the
// unit face normals nxi = (e1 - c e2)/s, neta = (e2 - c e1)/s
struct Frame {
  Real r[3], e1[3], e2[3], nxi[3], neta[3];
};
KOKKOS_INLINE_FUNCTION void MakeFrame(const int p, const Real xi, const Real eta,
                                      Frame &f) {
  cubed_sphere::PanelTangents(p, xi, eta, f.e1, f.e2);
  cubed_sphere::PanelToCart(p, xi, eta, f.r);
  const Real c = Dot3(f.e1, f.e2);
  const Real s = sqrt(1.0 - c*c);
  for (int q=0; q<3; ++q) {
    f.nxi[q] = (f.e1[q] - c*f.e2[q])/s;
    f.neta[q] = (f.e2[q] - c*f.e1[q])/s;
  }
}
// angles of cell / face positions of block m (x2 = xi, x3 = eta)
KOKKOS_INLINE_FUNCTION Real XiC(const SeamCtx &d, const int m, const int j) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*CellCenterX(j - d.js, d.nx2, z.x2min, z.x2max);
}
KOKKOS_INLINE_FUNCTION Real XiF(const SeamCtx &d, const int m, const int j) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*LeftEdgeX(j - d.js, d.nx2, z.x2min, z.x2max);
}
KOKKOS_INLINE_FUNCTION Real EtaC(const SeamCtx &d, const int m, const int k) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*CellCenterX(k - d.ks, d.nx3, z.x3min, z.x3max);
}
KOKKOS_INLINE_FUNCTION Real EtaF(const SeamCtx &d, const int m, const int k) {
  const RegionSize &z = d.size.d_view(m);
  return 0.25*M_PI*LeftEdgeX(k - d.ks, d.nx3, z.x3min, z.x3max);
}
// Cartesian velocity of cell (k,j,i): w0 holds v.rhat and the CONTRAVARIANT (v^xi,
// v^eta) on the unit tangents
KOKKOS_INLINE_FUNCTION void CellV(const SeamCtx &d, const Frame &f, const int m,
                                  const int k, const int j, const int i, Real v[3]) {
  const Real a = d.w(m,IVX,k,j,i), b = d.w(m,IVY,k,j,i), c = d.w(m,IVZ,k,j,i);
  for (int q=0; q<3; ++q) v[q] = a*f.r[q] + b*f.e1[q] + c*f.e2[q];
}
// Cartesian a-weighted face velocity of a UctFaceStore slot triple (3,4,5)
// dir 1: x1 face {rhat, e_xi, f2 = neta}; 2: x2 face {nxi, e_eta, rhat};
// 3: x3 face {neta, rhat, e_xi}
KOKKOS_INLINE_FUNCTION void FaceV(const DvceArray5D<Real> &u, const Frame &f,
                                  const int dir, const int m, const int k, const int j,
                                  const int i, Real v[3]) {
  const Real s3 = u(m,3,k,j,i), s4 = u(m,4,k,j,i), s5 = u(m,5,k,j,i);
  for (int q=0; q<3; ++q) {
    if (dir == 1) {
      v[q] = s3*f.r[q] + s4*f.e1[q] + s5*f.neta[q];
    } else if (dir == 2) {
      v[q] = s3*f.nxi[q] + s4*f.e2[q] + s5*f.r[q];
    } else {
      v[q] = s3*f.neta[q] + s4*f.r[q] + s5*f.e1[q];
    }
  }
}

// whether the along-seam index runs reversed between block m's seam sd and its
// neighbour panel (the seam_raw index map reverses it; the packets' lower/upper and the
// along-seam EMF sign must follow)
KOKKOS_INLINE_FUNCTION bool SeamRev(const SeamCtx &d, const int m, const int sd) {
  const int nid = (sd == 0) ? 8 : (sd == 1) ? 12 : (sd == 2) ? 24 : 28;
  const int pp = d.panel.d_view(m);
  const int qp = d.nghbr.d_view(m,nid).panel;
  const RegionSize &z = d.size.d_view(m);
  Real xi, eta;
  if (sd < 2) {
    xi = (sd == 0) ? -0.25*M_PI : 0.25*M_PI;
    eta = 0.25*M_PI*0.5*(z.x3min + z.x3max);
  } else {
    eta = (sd == 2) ? -0.25*M_PI : 0.25*M_PI;
    xi = 0.25*M_PI*0.5*(z.x2min + z.x2max);
  }
  Real e1[3], e2[3], r[3];
  cubed_sphere::PanelTangents(pp, xi, eta, e1, e2);
  cubed_sphere::PanelToCart(pp, xi, eta, r);
  const Real *sp = (sd < 2) ? e2 : e1;
  Real sv[3] = {sp[0], sp[1], sp[2]};
  Real xq, eq;
  cubed_sphere::CartToPanel(qp, r, xq, eq);
  Real f1[3], f2[3];
  cubed_sphere::PanelTangents(qp, xq, eq, f1, f2);
  const Real *sq = (fabs(xq) > fabs(eq)) ? f2 : f1;
  return (sv[0]*sq[0] + sv[1]*sq[1] + sv[2]*sq[2]) < 0.0;
}
KOKKOS_INLINE_FUNCTION bool OnPanelEdge(const SeamCtx &d, const int m, const int sd) {
  const RegionSize &z = d.size.d_view(m);
  if (sd == 0) return z.x2min <= -1.0 + 1.0e-12;
  if (sd == 1) return z.x2max >= 1.0 - 1.0e-12;
  if (sd == 2) return z.x3min <= -1.0 + 1.0e-12;
  return z.x3max >= 1.0 - 1.0e-12;
}


// ---- shared_o helpers ------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION Real MinMod(const Real a, const Real b) {
  return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);
}
// own-grid carry to the seam: value at the seam side of row r0 from rows r0, r1, r2
// (r1, r2 = one and two rows away from the seam). shared_o uses mode 3, the one-sided
// QUADRATIC through the three own rows where they are monotone, else the row value (the
// minmod-linear mode 0 was measured 3x worse at the seam, plain linear 2x: cs_test iprob
// 9 per-step n64 seam L1(B) 1.86e-6 / 1.06e-6 / 5.77e-7 for modes 0 / 1 / 3)
KOKKOS_INLINE_FUNCTION Real OwnCarry(const Real v0, const Real v1, const Real v2,
                                     const int mode = 0) {
  if (mode == 1) return v0 + 0.5*(v0 - v1);                          // linear
  if (mode == 2) return v0 + 0.5*(v0 - v1) + 0.125*(v0 - 2.0*v1 + v2)*3.0;  // quadratic
  if (mode == 3) {   // quadratic where the three own rows are monotone, else flat
    return ((v0 - v1)*(v1 - v2) > 0.0) ? (1.875*v0 - 1.25*v1 + 0.375*v2) : v0;
  }
  return v0 + 0.5*MinMod(v0 - v1, v1 - v2);
}
// plane normal to t: (u1, u2) with u1 along n's projection
KOKKOS_INLINE_FUNCTION void PlaneBasis(const Real t[3], const Real n[3], Real u1[3],
                                       Real u2[3]) {
  const Real tn = Dot3(n, t);
  for (int q=0; q<3; ++q) u1[q] = n[q] - tn*t[q];
  const Real un = sqrt(Dot3(u1, u1));
  for (int q=0; q<3; ++q) u1[q] /= un;
  u2[0] = t[1]*u1[2] - t[2]*u1[1];
  u2[1] = t[2]*u1[0] - t[0]*u1[2];
  u2[2] = t[0]*u1[1] - t[1]*u1[0];
}
// t x g_f (3-D), g_f = S^-1 n_f, S = sum over the nn normals (in the plane normal to t)
KOKKOS_INLINE_FUNCTION void TxG(const Real u1[3], const Real u2[3], const Real (*nl)[3],
                                const int nn, const Real nf[3], Real tg[3]) {
  Real s00 = 0.0, s01 = 0.0, s11 = 0.0;
  for (int a=0; a<nn; ++a) {
    const Real p0 = Dot3(nl[a], u1), p1 = Dot3(nl[a], u2);
    s00 += p0*p0; s01 += p0*p1; s11 += p1*p1;
  }
  const Real det = s00*s11 - s01*s01;
  const Real f0 = Dot3(nf, u1), f1 = Dot3(nf, u2);
  const Real g0 = ( s11*f0 - s01*f1)/det, g1 = (-s01*f0 + s00*f1)/det;
  for (int q=0; q<3; ++q) tg[q] = -g1*u1[q] + g0*u2[q];
}
// in-plane B (3-D) from its normal components on two faces
KOKKOS_INLINE_FUNCTION void SolveB(const Real u1[3], const Real u2[3], const Real na[3],
                                   const Real ba, const Real nb[3], const Real bb,
                                   Real bv[3]) {
  const Real a0 = Dot3(na, u1), a1 = Dot3(na, u2), b0 = Dot3(nb, u1), b1 = Dot3(nb, u2);
  const Real dd = a0*b1 - a1*b0;
  const Real x0 = (ba*b1 - a1*bb)/dd, x1 = (a0*bb - b0*ba)/dd;
  for (int q=0; q<3; ++q) bv[q] = x0*u1[q] + x1*u2[q];
}
// one side's share of a face's UCT term: (L side) a vn X + dL X, (R) (1-a) vn X - dR X
KOKKOS_INLINE_FUNCTION Real HalfC(const bool left, const Real a, const Real dl,
                                  const Real dr, const Real vn, const Real x) {
  return left ? (a*vn*x + dl*x) : ((1.0 - a)*vn*x - dr*x);
}

// the OWN partial of the 3-face UCT vertex rule (shared_o, MD21 at three
// faces at 120 deg, g = (2/3) n): this panel's half-terms of its two seam faces at the
// cube-vertex radial edge (kc, jc, i) of block m, corner (jlo, klo). E.r = -(own sum of
// the three panels' partials). Used by UctSeamOwnSide (shared_o) and MHD::CsVertexEmf3.
KOKKOS_INLINE_FUNCTION Real VertexOwnUct(const SeamCtx &d, const int m, const int p,
                                         const int i, const bool jlo, const bool klo) {
  const int jc = jlo ? d.js : d.je+1, kc = klo ? d.ks : d.ke+1;
  const int jv = jlo ? d.js : d.je, kv = klo ? d.ks : d.ke;
  Frame fv, fa, fb;
  MakeFrame(p, XiF(d, m, jc), EtaF(d, m, kc), fv);
  MakeFrame(p, XiF(d, m, jc), EtaC(d, m, kv), fa);
  MakeFrame(p, XiC(d, m, jv), EtaF(d, m, kc), fb);
  Real u1[3], u2[3]; PlaneBasis(fv.r, fa.nxi, u1, u2);
  Real ba_ = d.b2(m,kv,jc,i), bb_ = d.b3(m,kc,jv,i);
  if (d.nx2 > 1) {
    ba_ += 0.5*(d.b2(m,kv,jc,i) - d.b2(m,klo ? kv+1 : kv-1,jc,i));
    bb_ += 0.5*(d.b3(m,kc,jv,i) - d.b3(m,kc,jlo ? jv+1 : jv-1,i));
  }
  Real vA[3], vB[3];
  FaceV(d.u2, fa, 2, m, kv, jc, i, vA);
  FaceV(d.u3, fb, 3, m, kc, jv, i, vB);
  Real Bv[3]; SolveB(u1, u2, fa.nxi, ba_, fb.neta, bb_, Bv);
  Real tga[3], tgb[3];
  const Real *t = fv.r;
  tga[0] = t[1]*fa.nxi[2] - t[2]*fa.nxi[1]; tga[1] = t[2]*fa.nxi[0] - t[0]*fa.nxi[2];
  tga[2] = t[0]*fa.nxi[1] - t[1]*fa.nxi[0];
  tgb[0] = t[1]*fb.neta[2] - t[2]*fb.neta[1]; tgb[1] = t[2]*fb.neta[0] - t[0]*fb.neta[2];
  tgb[2] = t[0]*fb.neta[1] - t[1]*fb.neta[0];
  for (int q=0; q<3; ++q) {tga[q] *= 2.0/3.0; tgb[q] *= 2.0/3.0;}
  return HalfC(!jlo, d.u2(m,0,kv,jc,i), d.u2(m,1,kv,jc,i), d.u2(m,2,kv,jc,i),
               Dot3(vB, fa.nxi), Dot3(Bv, tga))
       + HalfC(!klo, d.u3(m,0,kc,jv,i), d.u3(m,1,kc,jv,i), d.u3(m,2,kc,jv,i),
               Dot3(vA, fb.neta), Dot3(Bv, tgb));
}
DvceArray5D<Real> xrec_;
// the shared RESISTIVE seam record (CsSeamResist): per seam position 0-2 along-seam edge
// at radial face i (H, A, eta2), 3-5 at ie+1, 6-8 / 9-11 radial edge at the lower /
// upper along-seam face of the cell, 12-14 / 15-17 cube-vertex corner at the lower /
// upper end (H, A, eta1)
DvceFaceFld5D<Real> *rsown_ = nullptr, *rsnbr_ = nullptr;
MeshBoundaryValuesCC *pbrs_ = nullptr;
MeshBoundaryValuesCC *pbx_ = nullptr, *pbxr_ = nullptr;
DvceArray5D<Real> xrr_;   // 3-variable four-block patch record (own-side / resistive)
// LOCAL fast path for the four-block records: the along-seam same-panel neighbour block of
// each block as a LOCAL index (columns: x2 inner, x2 outer, x3 inner, x3 outer; -1 = none
// or off-rank). When every needed neighbour is on its own rank (all ranks agree), the
// one-cell record copy is a direct device copy instead of a halo exchange (bitwise the
// same values: the halo copy is exact within a panel).
DvceArray2D<int> nbloc_;
bool xlocal_ = false;   // owned here (one MHD/process)
}  // namespace cs_uct

// release the file-scope Views before Kokkos::finalize (MHD::~MHD)
void CsSeamRelease() {
  cs_uct::xrec_ = DvceArray5D<Real>();
  if (cs_uct::pbx_ != nullptr) {delete cs_uct::pbx_; cs_uct::pbx_ = nullptr;}
  if (cs_uct::pbxr_ != nullptr) {delete cs_uct::pbxr_; cs_uct::pbxr_ = nullptr;}
  cs_uct::xrr_ = DvceArray5D<Real>();
  cs_uct::nbloc_ = DvceArray2D<int>();
  if (cs_uct::rsown_ != nullptr) {delete cs_uct::rsown_; cs_uct::rsown_ = nullptr;}
  if (cs_uct::rsnbr_ != nullptr) {delete cs_uct::rsnbr_; cs_uct::rsnbr_ = nullptr;}
  if (cs_uct::pbrs_ != nullptr) {delete cs_uct::pbrs_; cs_uct::pbrs_ = nullptr;}
}
using namespace cs_uct;  // NOLINT

// the radial ghost rows is-1 / ie+1 of the record <-> its second set at is / ie (the seam
// exchange covers the active radial rows only)
void CsSeamMoveRadialGhosts(MeshBlockPack *pp, DvceFaceFld5D<Real> &fs, const int nfs1,
                            const bool fwd) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pp->nmb_thispack - 1;
  auto f2 = fs.x2f;
  auto f3 = fs.x3f;
  auto &size = pp->pmb->mb_size;
  seam_for("cs_seamfs_rg2", DevExeSpace(), 0, nmb1, ks, ke, 0, 3, 0, nfs1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int c, const int n) {
    const bool lo = (c & 1) == 0;
    const bool up = (c & 2) != 0;    // 0: row is-1 <-> set B at is; 1: ie+1 <-> ie
    if (lo && !(size.d_view(m).x2min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x2max >= 1.0 - 1.0e-12)) return;
    const int jf = lo ? js : je+1;
    const int ig = up ? ie+1 : is-1, ib = up ? ie : is;
    if (fwd) {f2(m,nfs1+n,k,jf,ib) = f2(m,n,k,jf,ig);}
    else {f2(m,n,k,jf,ig) = f2(m,nfs1+n,k,jf,ib);}
  });
  seam_for("cs_seamfs_rg3", DevExeSpace(), 0, nmb1, js, je, 0, 3, 0, nfs1-1,
  KOKKOS_LAMBDA(const int m, const int j, const int c, const int n) {
    const bool lo = (c & 1) == 0;
    const bool up = (c & 2) != 0;
    if (lo && !(size.d_view(m).x3min <= -1.0 + 1.0e-12)) return;
    if (!lo && !(size.d_view(m).x3max >= 1.0 - 1.0e-12)) return;
    const int kf = lo ? ks : ke+1;
    const int ig = up ? ie+1 : is-1, ib = up ? ie : is;
    if (fwd) {f3(m,nfs1+n,kf,j,ib) = f3(m,n,kf,j,ig);}
    else {f3(m,n,kf,j,ig) = f3(m,nfs1+n,kf,j,ib);}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsSeamInit
//! \brief <mhd>/cs_seam_mhd = none | faces | shared, <mhd>/cs_seam_emf_avg (diagnostic)

void MHD::CsSeamInit(ParameterInput *pin) {
  auto *pp = pmy_pack;
  const std::string v = pin->GetString("mhd","cs_seam_mhd");
  if (v == "none") return;
  if (v != "shared_o") {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/cs_seam_mhd = '" << v << "' (none | shared_o)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  cs_seam_mhd = true;
  if (!pp->pmesh->use_cubed_sphere || !pp->pmesh->three_d ||
      recon_method != ReconstructionMethod::plm ||
      (use_wellbalance_dynamic && use_wb_x2) || use_wellbalance_static_reconst_perturb ||
      pp->pmesh->multilevel || pp->pcoord->is_special_relativistic ||
      pp->pcoord->is_general_relativistic || use_fofc) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cs_seam_mhd needs a 3-D cubed sphere, reconstruct = plm, no x2/x3 "
              << "well-balanced or static-perturbation reconstruction, a uniform mesh, "
              << "Newtonian MHD, no fofc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!(emf_uct || emf_uct_hlld)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/cs_seam_mhd = shared_o needs <mhd>/emf_scheme = uct or uct_hlld"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nmb = std::max((pp->nmb_thispack), (pp->pmesh->nmb_maxperrank));
  auto &ind = pp->pmesh->mb_indcs;
  const int n1 = ind.nx1 + 2*ind.ng, n2 = ind.nx2 + 2*ind.ng, n3 = ind.nx3 + 2*ind.ng;
  const int nvars = nmhd + nscalars;
  const int nfs1 = nvars + NDERIVED + 3;
  cs_seam_nfs = 2*nfs1;
  pbval_fs = new MeshBoundaryValuesCC(pp, pin, false);
  pbval_fs->seam_raw = true;
  pbval_fs->seam_xform_b = nvars + NDERIVED + 1;
  pbval_fs->seam_xform_stride = nfs1;
  pbval_fs->InitializeBuffers(cs_seam_nfs);
  pseamfs = new DvceFaceFld5D<Real>("cs_seamfs_mhd", nmb, cs_seam_nfs, n3, n2, n1);
  Kokkos::realloc(pseamfs->x1f, nmb, cs_seam_nfs, 1, 1, 1);
  // the four-block seam-edge record X (only the first ghost layer is read)
  Kokkos::realloc(xrec_, nmb, kNXA, n3, n2, n1);
  pbx_ = new MeshBoundaryValuesCC(pp, pin, false);
  pbx_->exch_depth = 1;
  pbx_->InitializeBuffers(kNXA);
  pbx_->SetVectorPairs(kNXA, {});    // all scalars (only within-panel ghosts are read)
  pbx_->cs_noresample = true;
  Kokkos::realloc(xrr_, nmb, 3, n3, n2, n1);
  pbxr_ = new MeshBoundaryValuesCC(pp, pin, false);
  pbxr_->exch_depth = 1;
  pbxr_->InitializeBuffers(3);
  pbxr_->SetVectorPairs(3, {});
  pbxr_->cs_noresample = true;
  {
    auto &ng_ = pp->pmb->nghbr;
    auto &gid_ = pp->pmb->mb_gid;
    auto &pan_ = pp->pmb->mb_panel;
    auto &lev_ = pp->pmb->mb_lev;
    auto &sz_ = pp->pmb->mb_size;
    ng_.sync_host(); gid_.sync_host(); pan_.sync_host(); lev_.sync_host(); sz_.sync_host();
    Kokkos::realloc(nbloc_, nmb, 4);
    auto hnb = Kokkos::create_mirror_view(nbloc_);
    const int nid[4] = {8, 12, 24, 28};
    int ok = 1;
    for (int m=0; m<nmb; ++m) {
      const RegionSize &z = sz_.h_view(m);
      const bool pe[4] = {z.x2min <= -1.0 + 1.0e-12, z.x2max >= 1.0 - 1.0e-12,
                          z.x3min <= -1.0 + 1.0e-12, z.x3max >= 1.0 - 1.0e-12};
      for (int c=0; c<4; ++c) {
        hnb(m,c) = -1;
        const NeighborBlock &nb = ng_.h_view(m,nid[c]);
        if (nb.gid >= 0 && nb.rank == global_variable::my_rank &&
            nb.panel == pan_.h_view(m) && nb.lev == lev_.h_view(m)) {
          for (int mm=0; mm<nmb; ++mm) {if (gid_.h_view(mm) == nb.gid) hnb(m,c) = mm;}
        }
      }
      // needed: the along-seam neighbours of a seam-side block whose end is not a
      // panel edge (x2 seams run along x3 -> columns 2/3; x3 seams -> columns 0/1)
      for (int c=0; c<4; ++c) {
        const bool need = !pe[c] && ((c < 2) ? (pe[2] || pe[3]) : (pe[0] || pe[1]));
        if (need && hnb(m,c) < 0) ok = 0;
      }
    }
    Kokkos::deep_copy(nbloc_, hnb);
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    xlocal_ = (ok == 1);
  }
  // the shared resistive seam EMF (CsSeamResist; used only with a resistivity)
  pbrs_ = new MeshBoundaryValuesCC(pp, pin, false);
  pbrs_->seam_raw = true;
  pbrs_->seam_xform_off = true;
  const int nr = 4;
  pbrs_->InitializeBuffers(nr);
  rsown_ = new DvceFaceFld5D<Real>("cs_resseam_own", nmb, nr, n3, n2, n1);
  rsnbr_ = new DvceFaceFld5D<Real>("cs_resseam_nbr", nmb, nr, n3, n2, n1);
  Kokkos::realloc(rsown_->x1f, nmb, nr, 1, 1, 1);
  Kokkos::realloc(rsnbr_->x1f, nmb, nr, 1, 1, 1);
  // the own-side seam EMF partials exchanged across the seam (6 per seam face)
  const int nu = 6;
  pbval_us = new MeshBoundaryValuesCC(pp, pin, false);
  pbval_us->seam_raw = true;
  pbval_us->seam_xform_off = true;
  pbval_us->InitializeBuffers(nu);
  pusown = new DvceFaceFld5D<Real>("cs_uctseam_own", nmb, nu, n3, n2, n1);
  pusnbr = new DvceFaceFld5D<Real>("cs_uctseam_nbr", nmb, nu, n3, n2, n1);
  Kokkos::realloc(pusown->x1f, nmb, nu, 1, 1, 1);
  Kokkos::realloc(pusnbr->x1f, nmb, nu, 1, 1, 1);
  if (global_variable::my_rank == 0) {
    std::cout << "MHD: cubed-sphere seam = "
              << "shared_o (faces + seam-edge/vertex UCT EMF as the sum of each "
              << "panel's own-grid contributions)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsSeamFaceStates
//! \brief (1): fill pseamfs with this panel's own seam face states (rows is-1..ie+1), move
//! the radial ghost rows into the second record set, exchange, move them back. After it
//! pseamfs holds, at every seam face of this rank, the NEIGHBOUR's own state.

// copy the record X (variables v0 .. v0+nv-1) of the along-seam same-panel neighbour's end
// cell (seam row) into this block's ghost cell -1 / nc of that row (local fast path)
void SeamXLocalCopy(MeshBlockPack *pp, DvceArray5D<Real> X, const int v0, const int nv) {
  auto &indcs = pp->pmesh->mb_indcs;
  SeamCtx d;
  d.panel = pp->pmb->mb_panel; d.size = pp->pmb->mb_size; d.nghbr = pp->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nc = indcs.nx2;
  auto nb = cs_uct::nbloc_;
  seam_for("seam_xlocal", DevExeSpace(), 0, pp->nmb_thispack-1, 0, 3, 0, 1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int up, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int mn = nb(m, x2t ? (up ? 3 : 2) : (up ? 1 : 0));
    if (mn < 0) return;
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int cb = x2t ? d.ks : d.js;
    const int g = up ? nc : -1, sc = up ? 0 : nc-1;
    for (int v=v0; v<v0+nv; ++v) {
      if (x2t) {X(m,v,cb+g,r0,i) = X(mn,v,cb+sc,r0,i);}
      else     {X(m,v,r0,cb+g,i) = X(mn,v,r0,cb+sc,i);}
    }
  });
}

void MHD::CsSeamFaceStates() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nvars = nmhd + nscalars;
  auto &eos_ = peos->eos_data;
  const int nder = eos_.IsGeneral() ? NDERIVED : 0;
  const int ob = nvars + NDERIVED;
  const int nfs1 = cs_seam_nfs/2;
  DvceArray2D<int> nomlat;
  cs_seam::CsSeamPLMFill(pmy_pack, w0, nvars, 0, -1, eos_.pfloor, eos_.plm_c1, *pseamfs,
                         false, nomlat, is-1, ie+1);
  cs_seam::CsSeamPLMFill(pmy_pack, wder, nder, nvars, IDPR, eos_.pfloor, eos_.plm_c1,
                         *pseamfs, false, nomlat, is-1, ie+1);
  cs_seam::CsSeamPLMFill(pmy_pack, bcc0, 3, ob, -1, eos_.pfloor, eos_.plm_c1, *pseamfs,
                         false, nomlat, is-1, ie+1);
  if (wb_all) {   // the edge cells' own WB seam face states (shared code)
    hydro::WbAllCtx c = WbCtx();
    hydro::WbAllSeamFaces(c, *pseamfs, nvars, nder, use_wb_rho, is-1, ie+1);
  }
  cs_seam::CsSeamCovariant(pmy_pack, *pseamfs, IVY, 0, is-1, ie+1);
  cs_seam::CsSeamCovariant(pmy_pack, *pseamfs, ob+1, 1, is-1, ie+1);
  CsSeamMoveRadialGhosts(pmy_pack, *pseamfs, nfs1, true);
  pbval_fs->InitFluxSeamRecv(cs_seam_nfs);
  pbval_fs->PackAndSendFluxSeamCC(*pseamfs);
  while (pbval_fs->RecvAndUnpackFluxSeamCC(*pseamfs) != TaskStatus::complete) {}
  pbval_fs->ClearFluxSend();
  pbval_fs->ClearFluxRecv();
  CsSeamMoveRadialGhosts(pmy_pack, *pseamfs, nfs1, false);
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::SeamXRecord
//! \brief Fills the per-cell seam-row record X (seam-face a, dL, dR and the face
//! velocity) and exchanges it through a one-layer halo, so that a block at a
//! four-block seam edge sees the seam-face store of its along-seam neighbour.

void MHD::SeamXRecord() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  SeamCtx d;
  d.u1 = uct1; d.u2 = uct2; d.u3 = uct3; d.w = w0;
  d.b1 = b0.x1f; d.b2 = b0.x2f; d.b3 = b0.x3f; d.e21 = e2x1; d.e31 = e3x1;
  d.o2 = pusown->x2f; d.o3 = pusown->x3f; d.n2 = pusnbr->x2f; d.n3 = pusnbr->x3f;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nc = indcs.nx2;
  auto X = xrec_;
  seam_for("uct_seam_xfill", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const int p = d.panel.d_view(m);
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int cb = x2t ? d.ks : d.js;
    const int sf = x2t ? (lo ? d.js : d.je+1) : (lo ? d.ks : d.ke+1);
    const int k0 = x2t ? (cb + c) : r0, j0 = x2t ? r0 : (cb + c);
    const int ksf = x2t ? k0 : sf, jsf = x2t ? sf : j0;
    const int f = 0;
    const DvceArray5D<Real> &us = x2t ? d.u2 : d.u3;
    X(m,f+0,k0,j0,i) = us(m,0,ksf,jsf,i);
    X(m,f+1,k0,j0,i) = us(m,1,ksf,jsf,i);
    X(m,f+2,k0,j0,i) = us(m,2,ksf,jsf,i);
    Frame fs;
    if (x2t) {MakeFrame(p, XiF(d, m, sf), EtaC(d, m, k0), fs);}
    else     {MakeFrame(p, XiC(d, m, j0), EtaF(d, m, sf), fs);}
    Real v[3];
    FaceV(us, fs, x2t ? 2 : 3, m, ksf, jsf, i, v);
    for (int q=0; q<3; ++q) X(m,f+3+q,k0,j0,i) = v[q];
  });
  if (xlocal_) {
    SeamXLocalCopy(pmy_pack, X, 0, kNXA);
  } else {
    while (pbx_->InitRecv(kNXA) == TaskStatus::incomplete) {}
    while (pbx_->PackAndSendCC(X, X) == TaskStatus::incomplete) {}
    while (pbx_->RecvAndUnpackCC(X, X) == TaskStatus::incomplete) {}
    while (pbx_->ClearSend() == TaskStatus::incomplete) {}
    while (pbx_->ClearRecv() == TaskStatus::incomplete) {}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsVertexEmf3
//! \brief <mesh>/cs_vertex_3cell (without the shared seam, cs_seam_mhd <
//! shared): the cube-vertex radial EMF from the THREE real cells and the three seam faces
//! that meet there, for every emf_scheme. Each panel P writes its OWN estimate E_P into
//! x1e at the vertex edge, formed from its corner cell and its two seam faces only, and
//! the existing three-way vertex mean (flux_correct_fc AveragePanelCornerEMF) of the
//! three panels' estimates IS the 3-cell rule:
//! UCT / UCT-HLLD: E_P = -3 x (P's partial of the shared_o 3-face rule, VertexOwnUct),
//! so the mean is exactly shared_o's E = -(own + q + r).
//! GS07: E = (1/3) [ sum_f E_f + sum_f (E_g(X_f) - E_cc(X_f)) ] over the three seam faces
//! f, X_f = the upwind cell of f (mass flux; emf_smooth weights; 1/2 at zero flux)
//! and g(X_f) its other vertex face -- the GS07 corner rule with three faces in
//! place of four (constant E exact; each pick owned by its upwind cell), so
//! E_P = (E_a + E_b)/2 + w_a (E_b - E_cc) + w_b (E_a - E_cc), w_f = P's upwind
//! weight on its seam face f (a = x2 seam face, b = x3 seam face).
//! The old 4-cell stencil read the corner block (e1_cc and the ghost-ghost face EMFs).

void MHD::CsVertexEmf3() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  SeamCtx d;
  d.u2 = uct2; d.u3 = uct3; d.w = w0;
  d.b2 = b0.x2f; d.b3 = b0.x3f;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool uct_ = emf_uct || emf_uct_hlld;
  const Real sm_ = emf_smooth;
  const bool has_e_ = (nmhd == 5);
  const Real cs2_ = has_e_ ? 0.0 : SQR(peos->eos_data.iso_cs);
  auto e1 = efld.x1e;
  auto e1x2_ = e1x2, e1x3_ = e1x3, e1cc_ = e1_cc;
  auto fx2 = uflx.x2f, fx3 = uflx.x3f;
  auto w0_ = w0;
  seam_for("cs_vertex_emf3", DevExeSpace(), 0, nmb1, 0, 3, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int cn, const int i) {
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    if (!OnPanelEdge(d, m, jlo ? 0 : 1) || !OnPanelEdge(d, m, klo ? 2 : 3)) return;
    const int jc = jlo ? d.js : d.je+1, kc = klo ? d.ks : d.ke+1;
    const int jv = jlo ? d.js : d.je, kv = klo ? d.ks : d.ke;
    if (uct_) {
      e1(m,kc,jc,i) = -3.0*VertexOwnUct(d, m, d.panel.d_view(m), i, jlo, klo);
      return;
    }
    auto re = [&](const int k, const int j) {
      return w0_(m,IDN,k,j,i)*(has_e_ ? fabs(w0_(m,IEN,k,j,i)) : cs2_);
    };
    // weight of the LOWER-index side of a face (GS07: the step; emf_smooth: smooth)
    auto wlo = [&](const Real fd, const Real ra, const Real rb) {
      if (sm_ >= 0.0) return EmfSmoothW(fd, sm_, 0.5*(ra + rb));
      return (fd > 0.0) ? 1.0 : ((fd < 0.0) ? 0.0 : 0.5);
    };
    const Real ea = e1x2_(m,kv,jc,i), eb = e1x3_(m,kc,jv,i), ec = e1cc_(m,kv,jv,i);
    const int jn = jlo ? jv-1 : jv+1, kn = klo ? kv-1 : kv+1;   // the seam ghost cells
    const Real r0 = re(kv,jv);
    const Real wa0 = wlo(fx2(m,IDN,kv,jc,i), jlo ? re(kv,jn) : r0, jlo ? r0 : re(kv,jn));
    const Real wb0 = wlo(fx3(m,IDN,kc,jv,i), klo ? re(kn,jv) : r0, klo ? r0 : re(kn,jv));
    const Real wa = jlo ? (1.0 - wa0) : wa0, wb = klo ? (1.0 - wb0) : wb0;
    e1(m,kc,jc,i) = 0.5*(ea + eb) + wa*(eb - ec) + wb*(ea - ec);
  });
  if (uct_) return;
  // GS07: the seam radial edges within ng-1 cells of a vertex. Their 4-cell stencil
  // takes the transverse face EMF of the seam GHOST row (column), whose reconstruction
  // reaches into the corner block. Each panel evaluates the same GS07 formula with the
  // ghost row (column) replaced by its own seam row (column) -- a mirror image across
  // the seam -- and the seam-edge average of the two panels' values combines their own
  // rows. Edges farther along the seam are unchanged.
  const int ng = indcs.ng;
  seam_for("cs_vertex_emf3_adj", DevExeSpace(), 0, nmb1, 0, 7, 1, ng-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int cs, const int s, const int i) {
    const int cn = cs & 3;
    const bool alongj = (cs < 4);   // true: the x3 seam (edges along j), false: x2 seam
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    if (!OnPanelEdge(d, m, jlo ? 0 : 1) || !OnPanelEdge(d, m, klo ? 2 : 3)) return;
    const int jc = jlo ? d.js : d.je+1, kc = klo ? d.ks : d.ke+1;
    const int k = alongj ? kc : (kc + (klo ? s : -s));
    const int j = alongj ? (jc + (jlo ? s : -s)) : jc;
    if (alongj ? !(j > d.js && j <= d.je) : !(k > d.ks && k <= d.ke)) return;
    // the ghost row (k-1 or k) / column (j-1 or j) of this seam -> the own one
    auto mk = [&](const int kk) {
      return (alongj && kk < d.ks) ? d.ks : ((alongj && kk > d.ke) ? d.ke : kk);
    };
    auto mj = [&](const int jj) {
      return (!alongj && jj < d.js) ? d.js : ((!alongj && jj > d.je) ? d.je : jj);
    };
    auto re = [&](const int kk, const int jj) {
      return w0_(m,IDN,mk(kk),mj(jj),i)*(has_e_ ? fabs(w0_(m,IEN,mk(kk),mj(jj),i)) : cs2_);
    };
    auto pick = [&](const Real fd, const Real ra, const Real rb, const Real a, const Real b) {
      if (sm_ >= 0.0) return b + EmfSmoothW(fd, sm_, 0.5*(ra + rb))*(a - b);
      return (fd >= 0.0) ? a : b;
    };
    auto x2 = [&](const int kk, const int jj) {return e1x2_(m,mk(kk),jj,i);};
    auto x3 = [&](const int kk, const int jj) {return e1x3_(m,kk,mj(jj),i);};
    auto cc = [&](const int kk, const int jj) {return e1cc_(m,mk(kk),mj(jj),i);};
    const Real l3 = pick(fx2(m,IDN,mk(k-1),j,i), re(k-1,j-1), re(k-1,j),
                         x3(k,j-1) - cc(k-1,j-1), x3(k,j) - cc(k-1,j));
    const Real r3 = pick(fx2(m,IDN,mk(k),j,i), re(k,j-1), re(k,j),
                         x3(k,j-1) - cc(k,j-1), x3(k,j) - cc(k,j));
    const Real l2 = pick(fx3(m,IDN,k,mj(j-1),i), re(k-1,j-1), re(k,j-1),
                         x2(k-1,j) - cc(k-1,j-1), x2(k,j) - cc(k,j-1));
    const Real r2 = pick(fx3(m,IDN,k,mj(j),i), re(k-1,j), re(k,j),
                         x2(k-1,j) - cc(k-1,j), x2(k,j) - cc(k,j));
    e1(m,k,j,i) = 0.25*(l3 + r3 + l2 + r2 + x2(k-1,j) + x2(k,j) + x3(k,j-1) + x3(k,j));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::UctSeamOwnSide
//! \brief <mhd>/cs_seam_mhd = shared_o (USER 11-02 option 2): every seam / cube-vertex
//! edge EMF as the per-face sum E.t = -sum_f C_f, each panel adding ONLY the parts it can
//! form on its own grid: the full term of its own non-seam face, and the own-side half of
//! each shared seam face's term (a^s vn^s X^s +- d^s X^s), every side state built from the
//! two faces of ONE of its own cells (no stencil across the kink: own-face carries to the
//! seam are one-sided, minmod of the two own-row slopes). The partial sums are exchanged
//! as scalars (E.t is a line integral; only the along-seam edge orientation can flip) and
//! added on both panels: single-valued by construction. Normals: S = sum of the actual
//! face normals (the neighbour's own face = the mirror of this panel's about the seam
//! plane, exact by the cube's symmetry); at a vertex the three seam faces at 120 deg,
//! g = (2/3) n.

void MHD::UctSeamOwnSide() {
  SeamXRecord();   // the four-block seam edges: the neighbour block's seam stores
  const int ocm = 3;   // own-grid carry: monotone one-sided quadratic (OwnCarry)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  SeamCtx d;
  d.u1 = uct1; d.u2 = uct2; d.u3 = uct3; d.w = w0;
  d.b1 = b0.x1f; d.b2 = b0.x2f; d.b3 = b0.x3f; d.e21 = e2x1; d.e31 = e3x1;
  d.o2 = pusown->x2f; d.o3 = pusown->x3f; d.n2 = pusnbr->x2f; d.n3 = pusnbr->x3f;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = indcs.is; d.ie = indcs.ie; d.js = indcs.js; d.je = indcs.je;
  d.ks = indcs.ks; d.ke = indcs.ke; d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nc = indcs.nx2;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool wsf_ = wall_closed_ix1 && wall_stress_free;
  const BoundaryFlag rflag_ = BoundaryFlag::reflect;
  const BoundaryFlag uflag_ = BoundaryFlag::user;
  auto X = xrec_;
  seam_for("uct_seam_own", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const int p = d.panel.d_view(m);
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int iw = lo ? 1 : -1;                 // one row away from the seam
    const int r1 = r0 + iw, r2 = r0 + 2*iw;
    const int cb = x2t ? d.ks : d.js;
    const int sf = x2t ? (lo ? d.js : d.je+1) : (lo ? d.ks : d.ke+1);   // seam face index
    const bool pleft = !lo;                     // this panel on the L side of a seam face
    auto K = [&](const int a, const int r) { return x2t ? (cb + a) : r; };
    auto J = [&](const int a, const int r) { return x2t ? r : (cb + a); };
    auto SK = [&](const int a) { return x2t ? (cb + a) : sf; };
    auto SJ = [&](const int a) { return x2t ? sf : (cb + a); };
    auto FrC = [&](const int a, const int r, Frame &f) {
      MakeFrame(p, XiC(d, m, J(a,r)), EtaC(d, m, K(a,r)), f);
    };
    auto out = [&](const int v) -> Real & {
      return x2t ? d.o2(m, v, cb + c, sf, i) : d.o3(m, v, sf, cb + c, i);
    };
    const DvceArray5D<Real> &us = x2t ? d.u2 : d.u3;     // seam faces
    const DvceArray4D<Real> &bs = x2t ? d.b2 : d.b3;
    const int dirs = x2t ? 2 : 3;
    // one-sided (own rows) half slope toward the seam of the Cartesian velocity of cell a
    auto Hc = [&](const int a, const int ii, Real h[3]) {
      Frame f0, f1, f2;
      FrC(a, r0, f0); FrC(a, r1, f1); FrC(a, r2, f2);
      Real v0[3], v1[3], v2[3];
      CellV(d, f0, m, K(a,r0), J(a,r0), ii, v0);
      CellV(d, f1, m, K(a,r1), J(a,r1), ii, v1);
      CellV(d, f2, m, K(a,r2), J(a,r2), ii, v2);
      for (int q=0; q<3; ++q) h[q] = OwnCarry(v0[q], v1[q], v2[q], ocm) - v0[q];
    };
    auto Bown = [&](const DvceArray4D<Real> &b, const int a, const int ii, const bool fk) {
      // fk: the face index runs along the seam (an along-seam face a), else cell a
      const int kk0 = K(a, r0), jj0 = J(a, r0);
      if (x2t) return OwnCarry(b(m,kk0,r0,ii), b(m,kk0,r1,ii), b(m,kk0,r2,ii), ocm);
      return OwnCarry(b(m,r0,jj0,ii), b(m,r1,jj0,ii), b(m,r2,jj0,ii), ocm);
    };
    const bool vlo = OnPanelEdge(d, m, x2t ? 2 : 0), vhi = OnPanelEdge(d, m, x2t ? 3 : 1);
    // along-seam half slopes of own seam-row cell a (one-sided at a cube-vertex corner)
    auto Hs = [&](const int a, const bool upper, Real h[3]) {
      Frame f0; FrC(a, r0, f0);
      Real v0[3]; CellV(d, f0, m, K(a,r0), J(a,r0), i, v0);
      const bool clo = (a == 0 && vlo), chi = (a == nc-1 && vhi);
      if ((clo || chi) && nc > 1) {
        const int an = clo ? a+1 : a-1;
        Frame fn; FrC(an, r0, fn);
        Real vn[3]; CellV(d, fn, m, K(an,r0), J(an,r0), i, vn);
        for (int q=0; q<3; ++q) {
          const Real sl = clo ? (vn[q] - v0[q]) : (v0[q] - vn[q]);
          h[q] = upper ? 0.5*sl : -0.5*sl;
        }
      } else if (clo || chi) {
        for (int q=0; q<3; ++q) h[q] = 0.0;
      } else {
        Frame fm, fp; FrC(a-1, r0, fm); FrC(a+1, r0, fp);
        Real vm[3], vp[3];
        CellV(d, fm, m, K(a-1,r0), J(a-1,r0), i, vm);
        CellV(d, fp, m, K(a+1,r0), J(a+1,r0), i, vp);
        for (int q=0; q<3; ++q) h[q] = Half(vm[q], v0[q], vp[q], upper);
      }
    };
    auto Bseam = [&](const int a, const bool upper) {   // seam-face field of cell a, along
      auto B = [&](const int x) { return x2t ? bs(m,cb+x,sf,i) : bs(m,sf,cb+x,i); };
      const bool clo = (a == 0 && vlo), chi = (a == nc-1 && vhi);
      if (!(clo || chi)) return Edge(B(a-1), B(a), B(a+1), upper);
      if (nc == 1) return B(a);
      const Real sl = clo ? (B(a+1) - B(a)) : (B(a) - B(a-1));
      return B(a) + (upper ? 0.5 : -0.5)*sl;
    };

    auto XO = [&](const int cc, const int v) {   // X record of the neighbour block's cell
      return x2t ? X(m, v, cb+cc, r0, i) : X(m, v, r0, cb+cc, i);
    };
    // --- along-seam edges at radial faces i (slot 0) and ie+1 (slot 1)
    for (int s=0; s<2; ++s) {
      if (s == 1 && i != d.ie) {out(1) = 0.0; break;}
      const int ifc = (s == 0) ? i : d.ie+1;
      const bool wall_in = ifc == d.is &&
          ((mb_bcs.d_view(m,BoundaryFace::inner_x1) == rflag_) ||
           (wsf_ && mb_bcs.d_view(m,BoundaryFace::inner_x1) == uflag_));
      const bool wall_out = ifc == d.ie+1 &&
                            mb_bcs.d_view(m,BoundaryFace::outer_x1) == rflag_;
      const int k0 = K(c, r0), j0 = J(c, r0);
      if (wall_in || wall_out) {
        out(s) = -0.5*(x2t ? d.e31(m,k0,j0,ifc) : d.e21(m,k0,j0,ifc));
        continue;
      }
      Frame fe;
      if (x2t) {MakeFrame(p, XiF(d, m, sf), EtaC(d, m, cb+c), fe);}
      else     {MakeFrame(p, XiC(d, m, cb+c), EtaF(d, m, sf), fe);}
      const Real *t = x2t ? fe.e2 : fe.e1;
      const Real *nS = x2t ? fe.nxi : fe.neta;
      const Real *nE = fe.r;
      Real u1[3], u2[3]; PlaneBasis(t, nS, u1, u2);
      Real nl[4][3];
      for (int q=0; q<3; ++q) {nl[0][q] = nS[q]; nl[1][q] = nS[q]; nl[2][q] = nE[q];
                               nl[3][q] = nE[q];}
      Real tgE[3], tgS[3];
      TxG(u1, u2, nl, 4, nE, tgE); TxG(u1, u2, nl, 4, nS, tgS);
      // seam faces at radial cells ifc-1 (lo) and ifc (hi): own-cell radial carries
      Real vS[2][3], bS[2];
      Frame fc; FrC(c, r0, fc);
      for (int h=0; h<2; ++h) {
        const int ii = (h == 0) ? ifc-1 : ifc;
        bS[h] = Edge(bs(m,SK(c),SJ(c),ii-1), bs(m,SK(c),SJ(c),ii), bs(m,SK(c),SJ(c),ii+1),
                     h == 0);
        FaceV(us, fe, dirs, m, SK(c), SJ(c), ii, vS[h]);
        Real va[3], vb[3], vc[3];
        CellV(d, fc, m, k0, j0, ii-1, va); CellV(d, fc, m, k0, j0, ii, vb);
        CellV(d, fc, m, k0, j0, ii+1, vc);
        for (int q=0; q<3; ++q) vS[h][q] += Half(va[q], vb[q], vc[q], h == 0);
      }
      // own radial face
      const Real aE = d.u1(m,0,k0,j0,ifc), dlE = d.u1(m,1,k0,j0,ifc);
      const Real drE = d.u1(m,2,k0,j0,ifc);
      Real bE = Bown(d.b1, c, ifc, false);
      Real vE[3], h0[3], h1[3];
      FaceV(d.u1, fc, 1, m, k0, j0, ifc, vE);
      Hc(c, ifc-1, h0); Hc(c, ifc, h1);
      for (int q=0; q<3; ++q) vE[q] += 0.5*(h0[q] + h1[q]);
      Real Blo[3], Bhi[3];
      SolveB(u1, u2, nS, bS[0], nE, bE, Blo);
      SolveB(u1, u2, nS, bS[1], nE, bE, Bhi);
      const Real xl = Dot3(Blo, tgE), xr = Dot3(Bhi, tgE);
      Real part = aE*Dot3(vS[0], nE)*xl + (1.0 - aE)*Dot3(vS[1], nE)*xr - (drE*xr - dlE*xl);
      for (int h=0; h<2; ++h) {
        const int ii = (h == 0) ? ifc-1 : ifc;
        part += HalfC(pleft, us(m,0,SK(c),SJ(c),ii), us(m,1,SK(c),SJ(c),ii),
                      us(m,2,SK(c),SJ(c),ii), Dot3(vE, nS),
                      Dot3((h == 0) ? Blo : Bhi, tgS));
      }
      out(s) = part;
    }
    // --- radial edges at the along-seam faces c (slot 2) and c+1 (slot 3), in-block only
    for (int s=0; s<2; ++s) {
      const int fa = c + s;
      // fa = 0 / nc: a FOUR-BLOCK seam edge (a block end that is not a cube vertex), the
      // neighbour block's seam store from the X record (SeamXRecord)
      if (fa < 0 || fa > nc || (fa == 0 && vlo) || (fa == nc && vhi)) {
        out(2+s) = 0.0; continue;
      }
      Frame fe; MakeFrame(p, XiF(d, m, x2t ? sf : cb+fa), EtaF(d, m, x2t ? cb+fa : sf), fe);
      const Real *t = fe.r;
      const Real *sp = x2t ? fe.e2 : fe.e1;
      const Real *nS = x2t ? fe.nxi : fe.neta;
      // own along-seam face fa
      const int kf = x2t ? (cb + fa) : r0, jf = x2t ? r0 : (cb + fa);
      Frame ff;
      if (x2t) {MakeFrame(p, XiC(d, m, jf), EtaF(d, m, kf), ff);}
      else     {MakeFrame(p, XiF(d, m, jf), EtaC(d, m, kf), ff);}
      const DvceArray5D<Real> &ua = x2t ? d.u3 : d.u2;
      const DvceArray4D<Real> &ba = x2t ? d.b3 : d.b2;
      Real nE[3], nW[3];
      for (int q=0; q<3; ++q) nE[q] = x2t ? ff.neta[q] : ff.nxi[q];
      const Real ens = Dot3(nE, nS);
      for (int q=0; q<3; ++q) nW[q] = nE[q] - 2.0*ens*nS[q];
      Real u1[3], u2[3]; PlaneBasis(t, nS, u1, u2);
      Real nl[4][3];
      for (int q=0; q<3; ++q) {nl[0][q] = nS[q]; nl[1][q] = nS[q]; nl[2][q] = nE[q];
                               nl[3][q] = nW[q];}
      Real tgE[3], tgS[3];
      TxG(u1, u2, nl, 4, nE, tgE); TxG(u1, u2, nl, 4, nS, tgS);
      const Real aE = ua(m,0,kf,jf,i), dlE = ua(m,1,kf,jf,i), drE = ua(m,2,kf,jf,i);
      Real bE = x2t ? OwnCarry(ba(m,kf,r0,i), ba(m,kf,r1,i), ba(m,kf,r2,i), ocm)
                    : OwnCarry(ba(m,r0,jf,i), ba(m,r1,jf,i), ba(m,r2,jf,i), ocm);
      Real vE[3], h0[3], h1[3];
      FaceV(ua, ff, x2t ? 3 : 2, m, kf, jf, i, vE);
      Hc(fa-1, i, h0); Hc(fa, i, h1);
      for (int q=0; q<3; ++q) vE[q] += 0.5*(h0[q] + h1[q]);
      // seam faces of cells fa-1 (S1) and fa (S2), own-cell along-seam carries
      Real vS[2][3], bSv[2];
      for (int h=0; h<2; ++h) {
        const int cc = (h == 0) ? fa-1 : fa;
        Frame fs;
        if (x2t) {MakeFrame(p, XiF(d, m, sf), EtaC(d, m, cb+cc), fs);}
        else     {MakeFrame(p, XiC(d, m, cb+cc), EtaF(d, m, sf), fs);}
        const bool inb = (cc >= 0 && cc < nc);
        if (inb) {
          FaceV(us, fs, dirs, m, SK(cc), SJ(cc), i, vS[h]);
        } else {
          for (int q=0; q<3; ++q) vS[h][q] = XO(cc, 3+q);
        }
        Real hh[3]; Hs(cc, h == 0, hh);
        for (int q=0; q<3; ++q) vS[h][q] += hh[q];
        bSv[h] = Bseam(cc, h == 0);
      }
      Real B1[3], B2[3];
      SolveB(u1, u2, nS, bSv[0], nE, bE, B1);
      SolveB(u1, u2, nS, bSv[1], nE, bE, B2);
      const bool el = Dot3(nE, sp) > 0.0;      // L side of the own face = cell fa-1
      const Real x1 = Dot3(B1, tgE), x2 = Dot3(B2, tgE);
      const Real v1 = Dot3(vS[0], nE), v2 = Dot3(vS[1], nE);
      Real part = el ? (aE*v1*x1 + (1.0 - aE)*v2*x2 - (drE*x2 - dlE*x1))
                     : (aE*v2*x2 + (1.0 - aE)*v1*x1 - (drE*x1 - dlE*x2));
      for (int h=0; h<2; ++h) {
        const int cc = (h == 0) ? fa-1 : fa;
        const bool inb = (cc >= 0 && cc < nc);
        part += HalfC(pleft, inb ? us(m,0,SK(cc),SJ(cc),i) : XO(cc, 0),
                      inb ? us(m,1,SK(cc),SJ(cc),i) : XO(cc, 1),
                      inb ? us(m,2,SK(cc),SJ(cc),i) : XO(cc, 2), Dot3(vE, nS),
                      Dot3((h == 0) ? B1 : B2, tgS));
      }
      out(2+s) = part;
    }
    // --- cube-vertex radial edges at the lower (slot 4) / upper (slot 5) end of the seam
    for (int s=0; s<2; ++s) {
      const bool isv = (s == 0) ? (c == 0 && vlo) : (c == nc-1 && vhi);
      if (!isv) {out(4+s) = 0.0; continue;}
      // corner: x2-seam index jc/jv, x3-seam index kc/kv ( shared helper)
      const bool jlo = x2t ? lo : (s == 0), klo = x2t ? (s == 0) : lo;
      out(4+s) = VertexOwnUct(d, m, p, i, jlo, klo);
    }
  });
  // FOUR-BLOCK seam ends: the CANONICAL partial of such an edge is the one of the block
  // above it (its cell-0 lower-face partial); the block below takes it from the X record
  // (slot 6, exchanged once more), so all four blocks combine the same two numbers even
  // where their copies of the block-boundary face differ (no same-panel averaging)
  auto XR = xrr_;
  seam_for("uct_seam_own_xp", DevExeSpace(), 0, nmb1, 0, 3, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int cb = x2t ? d.ks : d.js;
    const int sf = x2t ? (lo ? d.js : d.je+1) : (lo ? d.ks : d.ke+1);
    const int k0 = x2t ? cb : r0, j0 = x2t ? r0 : cb;
    XR(m,0,k0,j0,i) = x2t ? d.o2(m,2,cb,sf,i) : d.o3(m,2,sf,cb,i);
  });
  if (xlocal_) {
    SeamXLocalCopy(pmy_pack, XR, 0, 1);
  } else {
    while (pbxr_->InitRecv(3) == TaskStatus::incomplete) {}
    while (pbxr_->PackAndSendCC(XR, XR) == TaskStatus::incomplete) {}
    while (pbxr_->RecvAndUnpackCC(XR, XR) == TaskStatus::incomplete) {}
    while (pbxr_->ClearSend() == TaskStatus::incomplete) {}
    while (pbxr_->ClearRecv() == TaskStatus::incomplete) {}
  }
  seam_for("uct_seam_own_xq", DevExeSpace(), 0, nmb1, 0, 3, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    if (OnPanelEdge(d, m, x2t ? 3 : 1)) return;      // upper end is a cube vertex
    const int r0 = x2t ? (lo ? d.js : d.je) : (lo ? d.ks : d.ke);
    const int cb = x2t ? d.ks : d.js;
    const int sf = x2t ? (lo ? d.js : d.je+1) : (lo ? d.ks : d.ke+1);
    const int kg = x2t ? (cb + nc) : r0, jg = x2t ? r0 : (cb + nc);
    if (x2t) {d.o2(m,3,cb+nc-1,sf,i) = XR(m,0,kg,jg,i);}
    else     {d.o3(m,3,sf,cb+nc-1,i) = XR(m,0,kg,jg,i);}
  });
  pbval_us->InitFluxSeamRecv(6);   // (the X record was exchanged at the top)
  pbval_us->PackAndSendFluxSeamCC(*pusown);
  while (pbval_us->RecvAndUnpackFluxSeamCC(*pusnbr) != TaskStatus::complete) {}
  pbval_us->ClearFluxSend();
  pbval_us->ClearFluxRecv();
  // assemble: E.t = -(own + neighbour) (+ the third panel at a vertex)
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  seam_for("uct_seam_own_sum", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const bool rev = SeamRev(d, m, sd);
    const int cb = x2t ? d.ks : d.js;
    const int sf = x2t ? (lo ? d.js : d.je+1) : (lo ? d.ks : d.ke+1);
    auto O = [&](const int v) {
      return x2t ? d.o2(m, v, cb + c, sf, i) : d.o3(m, v, sf, cb + c, i);
    };
    auto N = [&](const int v) {
      return x2t ? d.n2(m, v, cb + c, sf, i) : d.n3(m, v, sf, cb + c, i);
    };
    const Real sg = rev ? -1.0 : 1.0;
    for (int s=0; s<2; ++s) {
      if (s == 1 && i != d.ie) break;
      const int ifc = (s == 0) ? i : d.ie+1;
      const Real e = -(O(s) + sg*N(s));
      if (x2t) {e3(m,cb+c,sf,ifc) = e;} else {e2(m,sf,cb+c,ifc) = e;}
    }
    const bool vlo = OnPanelEdge(d, m, x2t ? 2 : 0), vhi = OnPanelEdge(d, m, x2t ? 3 : 1);
    if (c >= 1 || !vlo) {   // radial edge at the lower face of cell c (c = 0: four-block)
      const Real e = -(O(2) + N(rev ? 3 : 2));
      if (x2t) {e1(m,cb+c,sf,i) = e;} else {e1(m,sf,cb+c,i) = e;}
    }
    if (c == nc-1 && !vhi) {   // the four-block edge at the upper end
      const Real e = -(O(3) + N(rev ? 2 : 3));
      if (x2t) {e1(m,cb+nc,sf,i) = e;} else {e1(m,sf,cb+nc,i) = e;}
    }
  });
  seam_for("uct_seam_own_vtx", DevExeSpace(), 0, nmb1, 0, 3, d.is, d.ie,
  KOKKOS_LAMBDA(const int m, const int cn, const int i) {
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    const int sj = jlo ? 0 : 1, sk = klo ? 2 : 3;
    if (!OnPanelEdge(d, m, sj) || !OnPanelEdge(d, m, sk)) return;
    const int jc = jlo ? d.js : d.je+1, kc = klo ? d.ks : d.ke+1;
    const int jv = jlo ? d.js : d.je, kv = klo ? d.ks : d.ke;
    // own partial (from the x2-seam record, end klo), Q's (x2 seam), R's (x3 seam)
    const Real own = d.o2(m, klo ? 4 : 5, kv, jc, i);
    const bool r2 = SeamRev(d, m, sj), r3 = SeamRev(d, m, sk);
    const Real q = d.n2(m, (klo != r2) ? 4 : 5, kv, jc, i);
    const Real r = d.n3(m, (jlo != r3) ? 4 : 5, kc, jv, i);
    e1(m,kc,jc,i) = -(own + q + r);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CsSeamResist
//! \brief shared_g / shared_o with a resistivity: the RESISTIVE EMF eta J.t on every seam
//! edge (along-seam, radial, cube vertex), single-valued by construction. J.t is the
//! circulation of B around the edge through the cell centres around it (Stokes, as in
//! AddEMFGnomonicResist) split by panel: each panel adds the loop segments on ITS side --
//! the full segment across its own face and HALF of each segment across a shared seam
//! face, B along it from its own grid only (the seam-face normal field, shared, plus the
//! partner component extrapolated from its own two rows: no resampled ghost) -- and half
//! of the enclosed area (a quarter at a cube vertex; along-seam edges: times the sine of
//! the gnomonic angle, the projection on the plane normal to the edge). The partial
//! circulations, areas and the eta of its own cells next to the edge are exchanged as
//! scalars (the along-seam edge orientation sign only) and combined identically:
//! E_res.t = (eta_P + eta_Q)/2 * (H_P + H_Q)/(A_P + A_Q),
//! three panels at a vertex. AddEMFGnomonicResist skips the seam edges in this mode.
//! The radial edges at four-block seam ends use the same-panel ghost copies (exact).

void MHD::CsSeamResist(const DvceFaceFld4D<Real> &b0, const DvceArray4D<Real> &eta,
                       const DvceEdgeFld4D<Real> &jnorm, DvceEdgeFld4D<Real> &efld_) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nc = indcs.nx2;
  SeamCtx d;
  d.panel = pmy_pack->pmb->mb_panel; d.size = pmy_pack->pmb->mb_size;
  d.nghbr = pmy_pack->pmb->nghbr;
  d.is = is; d.ie = ie; d.js = js; d.je = je; d.ks = ks; d.ke = ke;
  d.nx2 = indcs.nx2; d.nx3 = indcs.nx3;
  auto o2 = rsown_->x2f, o3 = rsown_->x3f, n2 = rsnbr_->x2f, n3 = rsnbr_->x3f;
  // ---- cs_seam_res = mean (default): the TWO-CHART MEAN of the standard resistive
  // EMF (pass 2 computed every seam edge from each panel's own + ghost cells), the
  // shared_g idea applied to eta J: each panel exchanges its TOTAL seam-edge EMF (the
  // ideal part is already single-valued, so the mean leaves it unchanged) and both take
  // E = (E_P + s E_Q)/2, three panels at a cube vertex (sorted, then summed),
  // the same numbers in the same order -> single-valued without the post-hoc averaging.
  // Four-block seam ends: the upper block's value is canonical (one-layer record).
  auto e1 = efld_.x1e;
  auto e2 = efld_.x2e;
  auto e3 = efld_.x3e;
  auto XR = xrr_;
  seam_for("res_mean_xp", DevExeSpace(), 0, nmb1, 0, 3, is, ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? js : je) : (lo ? ks : ke);
    const int cb = x2t ? ks : js;
    const int sf = x2t ? (lo ? js : je+1) : (lo ? ks : ke+1);
    if (x2t) {XR(m,0,cb,r0,i) = e1(m,cb,sf,i);} else {XR(m,0,r0,cb,i) = e1(m,sf,cb,i);}
  });
  if (xlocal_) {
    SeamXLocalCopy(pmy_pack, XR, 0, 1);
  } else {
    while (pbxr_->InitRecv(3) == TaskStatus::incomplete) {}
    while (pbxr_->PackAndSendCC(XR, XR) == TaskStatus::incomplete) {}
    while (pbxr_->RecvAndUnpackCC(XR, XR) == TaskStatus::incomplete) {}
    while (pbxr_->ClearSend() == TaskStatus::incomplete) {}
    while (pbxr_->ClearRecv() == TaskStatus::incomplete) {}
  }
  seam_for("res_mean_own", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, is, ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const int r0 = x2t ? (lo ? js : je) : (lo ? ks : ke);
    const int cb = x2t ? ks : js;
    const int sf = x2t ? (lo ? js : je+1) : (lo ? ks : ke+1);
    const bool vhi = OnPanelEdge(d, m, x2t ? 3 : 1);
    auto out = [&](const int v) -> Real & {
      return x2t ? o2(m, v, cb + c, sf, i) : o3(m, v, sf, cb + c, i);
    };
    out(0) = x2t ? e3(m,cb+c,sf,i) : e2(m,sf,cb+c,i);
    out(1) = (i == ie) ? (x2t ? e3(m,cb+c,sf,ie+1) : e2(m,sf,cb+c,ie+1)) : 0.0;
    out(2) = x2t ? e1(m,cb+c,sf,i) : e1(m,sf,cb+c,i);
    if (c == nc-1 && !vhi) {   // four-block upper end: the upper block's value
      const int kg = x2t ? (cb + nc) : r0, jg = x2t ? r0 : (cb + nc);
      out(3) = XR(m,0,kg,jg,i);
    } else {
      out(3) = x2t ? e1(m,cb+c+1,sf,i) : e1(m,sf,cb+c+1,i);
    }
  });
  pbrs_->InitFluxSeamRecv(4);
  pbrs_->PackAndSendFluxSeamCC(*rsown_);
  while (pbrs_->RecvAndUnpackFluxSeamCC(*rsnbr_) != TaskStatus::complete) {}
  pbrs_->ClearFluxSend();
  pbrs_->ClearFluxRecv();
  seam_for("res_mean_sum", DevExeSpace(), 0, nmb1, 0, 3, 0, nc-1, is, ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int c, const int i) {
    if (!OnPanelEdge(d, m, sd)) return;
    const bool x2t = (sd < 2), lo = ((sd & 1) == 0);
    const bool rev = SeamRev(d, m, sd);
    const int cb = x2t ? ks : js;
    const int sf = x2t ? (lo ? js : je+1) : (lo ? ks : ke+1);
    auto O = [&](const int v) { return x2t ? o2(m, v, cb + c, sf, i) : o3(m, v, sf, cb + c, i); };
    auto N = [&](const int v) { return x2t ? n2(m, v, cb + c, sf, i) : n3(m, v, sf, cb + c, i); };
    const Real sg = rev ? -1.0 : 1.0;
    const Real ea = 0.5*(O(0) + sg*N(0));
    if (x2t) {e3(m,cb+c,sf,i) = ea;} else {e2(m,sf,cb+c,i) = ea;}
    if (i == ie) {
      const Real eb = 0.5*(O(1) + sg*N(1));
      if (x2t) {e3(m,cb+c,sf,ie+1) = eb;} else {e2(m,sf,cb+c,ie+1) = eb;}
    }
    const bool vlo = OnPanelEdge(d, m, x2t ? 2 : 0), vhi = OnPanelEdge(d, m, x2t ? 3 : 1);
    if (!(c == 0 && vlo)) {
      const Real er = 0.5*(O(2) + N(rev ? 3 : 2));
      if (x2t) {e1(m,cb+c,sf,i) = er;} else {e1(m,sf,cb+c,i) = er;}
    }
    if (c == nc-1 && !vhi) {
      const Real er = 0.5*(O(3) + N(rev ? 2 : 3));
      if (x2t) {e1(m,cb+nc,sf,i) = er;} else {e1(m,sf,cb+nc,i) = er;}
    }
  });
  seam_for("res_mean_vtx", DevExeSpace(), 0, nmb1, 0, 3, is, ie,
  KOKKOS_LAMBDA(const int m, const int cn, const int i) {
    const bool jlo = (cn & 1) == 0, klo = (cn & 2) == 0;
    const int sj = jlo ? 0 : 1, sk = klo ? 2 : 3;
    if (!OnPanelEdge(d, m, sj) || !OnPanelEdge(d, m, sk)) return;
    const int jc = jlo ? js : je+1, kc = klo ? ks : ke+1;
    const int jv = jlo ? js : je, kv = klo ? ks : ke;
    const bool r2 = SeamRev(d, m, sj), r3 = SeamRev(d, m, sk);
    Real a = o2(m, klo ? 2 : 3, kv, jc, i);
    Real b_ = n2(m, (klo != r2) ? 2 : 3, kv, jc, i);
    Real c_ = n3(m, (jlo != r3) ? 2 : 3, kc, jv, i);
    // sort, so that all three panels add the same numbers in the same order
    if (a > b_) {const Real t = a; a = b_; b_ = t;}
    if (b_ > c_) {const Real t = b_; b_ = c_; c_ = t;}
    if (a > b_) {const Real t = a; a = b_; b_ = t;}
    e1(m,kc,jc,i) = ((a + b_) + c_)/3.0;
  });
}

}  // namespace mhd

// (cs_edge_jb.hpp): the stage field for the shared seam current
