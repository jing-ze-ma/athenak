//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_cs.cpp
//! \brief cs-vetgd-1017 (CS_VETGD.md sect. 4A, stage S3): <rad_m1>/vet_gd on the
//! CUBED SPHERE.  The formal solution is the sp sweep of rad_m1_vetgd.cpp (same shells,
//! passes, segment rules, branch masks, boundary rules, twin, rotation); only the
//! lateral geometry differs, and it lives here in its own kernels (the sp kernels are
//! untouched: bitwise).
//!
//! BAND.  A block's band is its own panel's index space extended by vgd_w cells, as on
//! the sp wedge.  Beyond a panel edge the neighbour panel's cells are stored UNFOLDED:
//! the neighbour cell at depth d from the seam sits at the position where the own
//! panel's index line continues (the mirror image, across the seam's great circle, of
//! the own cell at depth d; the equiangular grid with equal N per panel maps onto
//! itself under that reflection, the CS1 mirror pairing).  Positions beyond TWO edges
//! (the corner quadrant: only three panels meet at a cube vertex) and depths >= N hold
//! no cell.  Every band position is filled per shell by exact copies of its owner's
//! interior value (on-rank lists + one MPI message per partner rank), keyed on the
//! global cell (panel, a, t): decomposition-invariant.
//!
//! READ.  The upwind point x (Cartesian, as on sp) is located on the panel c whose
//! outward normal it is most aligned with, at the equiangular chart coordinates (xi,
//! eta) of c, i.e. at the fractional NATIVE indices (fa, ft) of c (cell centres at
//! integers):
//!   * inside c's centre lattice: bilinear in (fa, ft) on c's 4 cells (= bilinear in
//!     index space, as the sp sweep);
//!   * in the SEAM STRIP (within half a cell of an edge of c, beyond c's edge-column
//!     centres): linear across the seam, weight 1 - dist on c's edge column (dist the
//!     distance beyond its centre in c's cells, <= 1/2) and dist on the neighbour's edge
//!     column, each column interpolated linearly along the seam in its OWN panel's
//!     tangential coordinate of x (continuous across the seam: the two charts agree on
//!     the shared edge);
//!   * in the VERTEX corner (beyond both edges' centres): the three cells around the
//!     cube vertex with spherical barycentric weights (clamped >= 0, renormalised).
//! All weights >= 0, sum 1 (exact for constants).  The branch mask of the sp sweep
//! (drop stencil cells holding the other pass's branch, renormalise) and vet_interp
//! auto / corner, vet_seg split / jump, the inner face (diffusion or emit), the top
//! bath and the tangent rule are the sp sweep's, term for term, over this stencil.
//!
//! HAND-OVER: D_rr - f_K into LAT0 (folded into slot 0 by VetLatBuild), the
//! vet_gd_replace surface q; cs-s4-1017 (stage S4): D_r,lat as the CARTESIAN tangential
//! vector v = D.rhat - D_rr rhat in LAT1..LAT3 (the operator / lagged form on the skew
//! cs faces: M1CsLat0, M1CsLatF, VetLatOpCs, VetLatStencilAddCs); the tangential
//! anisotropy is not handed over on cs (vet_gd_tangential = false).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "utils/repro_sum.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_vetseg.hpp"
#include "rad_m1/rad_m1_vgdray.hpp"
// the compact per-shell exchange over NCCL (build -DM1_VGD_NCCL, env
// ACC_NCCL=1; as the sp list halo, VgdNcclOn in rad_m1_vetgd.cpp)
#if defined(M1_VGD_NCCL) && MPI_PARALLEL_ENABLED && defined(KOKKOS_ENABLE_CUDA)
#include <nccl.h>
#define M1_VGDCS_NCCL 1
#else
#define M1_VGDCS_NCCL 0
#endif

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

// build generation of each compact-list cache entry (NOT a class member, so
// the dt-retry checkpoint does not restore it): DtRetryAfterRestore drops an entry that
// the failed attempt rebuilt in place (its restored key would name the old set)
int64_t g_vgd_csgen[16] = {0};        // (M1_VGD_CSCL_MAX)
int64_t g_vgd_csgen_save[16] = {0};

namespace {

// a stream event of the default device instance (the sweep's stream): record
// after a kernel, query from the host without a full fence.  No device backend: fence
struct VgdCsEv {
  VgdCsEv() = default;
  VgdCsEv(const VgdCsEv &) = delete;
  VgdCsEv &operator=(const VgdCsEv &) = delete;
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEvent_t e = nullptr;
  void Record(const DevExeSpace &ex = DevExeSpace()) {
    if (e == nullptr) {cudaEventCreateWithFlags(&e, cudaEventDisableTiming);}
    cudaEventRecord(e, ex.cuda_stream());
  }
  void StreamWait(const DevExeSpace &ex) const {
    cudaStreamWaitEvent(ex.cuda_stream(), e, 0);
  }
  bool Done() const {return e == nullptr || cudaEventQuery(e) != cudaErrorNotReady;}
  void Wait() const {
    if (e != nullptr) {cudaEventSynchronize(e);}
  }
  ~VgdCsEv() {
    if (e != nullptr) {cudaEventDestroy(e);}
  }
#elif defined(KOKKOS_ENABLE_HIP)
  hipEvent_t e = nullptr;
  void Record(const DevExeSpace &ex = DevExeSpace()) {
    if (e == nullptr) {(void) hipEventCreateWithFlags(&e, hipEventDisableTiming);}
    (void) hipEventRecord(e, ex.hip_stream());
  }
  void StreamWait(const DevExeSpace &ex) const {
    (void) hipStreamWaitEvent(ex.hip_stream(), e, 0);
  }
  bool Done() const {return e == nullptr || hipEventQuery(e) != hipErrorNotReady;}
  void Wait() const {
    if (e != nullptr) {(void) hipEventSynchronize(e);}
  }
  ~VgdCsEv() {
    if (e != nullptr) {(void) hipEventDestroy(e);}
  }
#else
  void Record(const DevExeSpace & = DevExeSpace()) {Kokkos::fence();}
  void StreamWait(const DevExeSpace &) const {}
  bool Done() const {return true;}
  void Wait() const {}
#endif
};
void VgdCsFatal(const std::string &msg) {
  VgdFatalAt(__FILE__, "<rad_m1>/vet_gd (cubed sphere): ", msg);
}

KOKKOS_INLINE_FUNCTION
Real VgdCsDet3(const Real *a, const Real *b, const Real *c) {
  return a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0])
         + a[2]*(b[0]*c[1] - b[1]*c[0]);
}

// the upwind stencil: up to 4 band positions (k, j) of the block and their weights
struct VgdCsSt {
  int n;
  int k[4], j[4];
  Real w[4];
};

// the stencil accessor of the shared upwind core (rad_m1_vgdray.hpp): the
// weighted sum over the stencil cells, in their order
struct VgdCsA {
  const VgdCsSt &st;
  KOKKOS_INLINE_FUNCTION int N() const {return st.n;}
  KOKKOS_INLINE_FUNCTION int K(const int e) const {return st.k[e];}
  KOKKOS_INLINE_FUNCTION int J(const int e) const {return st.j[e];}
  KOKKOS_INLINE_FUNCTION Real W(const int e) const {return st.w[e];}
  template <typename F>
  KOKKOS_INLINE_FUNCTION Real Interp(const F &f) const {
    Real s = 0.0;
    for (int e = 0; e < st.n; ++e) {s += st.w[e]*f(st.k[e], st.j[e]);}
    return s;
  }
};

// the stencil of the unit vector x for a block on panel p whose first interior cell is
// the global (j0, k0) of p; band index = global unfolded index - (j0, k0) + wb, band
// extents c2 x c3.  Returns false when a stencil cell is not in the band (dropped, the
// rest renormalised: counted as a clamped read by the caller).
template <class CM, class NB>
KOKKOS_INLINE_FUNCTION
bool VgdCsStencil(const Real x[3], const int p, const int j0, const int k0, const int wb,
                  const int c2, const int c3, const int nc, const Real h, const CM &cm,
                  const NB &nb, VgdCsSt &st) {
  int qc[4] = {0, 0, 0, 0}, qa[4] = {0, 0, 0, 0}, qt[4] = {0, 0, 0, 0};
  Real qw[4] = {0.0, 0.0, 0.0, 0.0};
  int nq = 0;
  const Real q4 = 0.25*M_PI;
  const Real nm1 = static_cast<Real>(nc - 1);
  const int c = cubed_sphere::FindPanel(x);
  Real xi, eta;
  cubed_sphere::CartToPanel(c, x, xi, eta);
  const Real fa = (xi + q4)/h - 0.5, ft = (eta + q4)/h - 0.5;
  const bool ina = (fa >= 0.0 && fa <= nm1), intt = (ft >= 0.0 && ft <= nm1);
  if (ina && intt) {
    // inside the centre lattice of c: bilinear in its native indices
    int a0 = static_cast<int>(floor(fa));
    int t0 = static_cast<int>(floor(ft));
    a0 = (a0 > nc - 2) ? (nc - 2) : a0;
    t0 = (t0 > nc - 2) ? (nc - 2) : t0;
    const Real ua = fa - a0, ut = ft - t0;
    for (int q = 0; q < 4; ++q) {
      const int da = q % 2, dt = q/2;
      qc[q] = c;
      qa[q] = a0 + da;
      qt[q] = t0 + dt;
      qw[q] = (da ? ua : (1.0 - ua))*(dt ? ut : (1.0 - ut));
    }
    nq = 4;
  } else if (ina != intt) {
    // the seam strip: c's edge column and the neighbour's, linear across the seam
    const bool sa = !ina;                    // the seam normal is c's a axis
    const Real fn = sa ? fa : ft;
    const Real fg = sa ? ft : fa;            // in [0, nc - 1]
    const bool hi = (fn > nm1);
    const int s = sa ? (hi ? 1 : 0) : (hi ? 3 : 2);
    const Real dist = hi ? (fn - nm1) : -fn;
    const int ne = hi ? (nc - 1) : 0;
    int g0 = static_cast<int>(floor(fg));
    g0 = (g0 > nc - 2) ? (nc - 2) : g0;
    const Real ug = fg - g0;
    const Real wc = 1.0 - dist;
    for (int q = 0; q < 2; ++q) {
      qc[q] = c;
      qa[q] = sa ? ne : (g0 + q);
      qt[q] = sa ? (g0 + q) : ne;
      qw[q] = wc*(q ? ug : (1.0 - ug));
    }
    const int c2p = nb(c, s);
    Real xi2, eta2;
    cubed_sphere::CartToPanel(c2p, x, xi2, eta2);
    const Real fa2 = (xi2 + q4)/h - 0.5, ft2 = (eta2 + q4)/h - 0.5;
    const Real oa = fmax(fmax(-fa2, fa2 - nm1), 0.0);
    const Real ot = fmax(fmax(-ft2, ft2 - nm1), 0.0);
    const bool na = (oa >= ot);             // the seam normal is c2p's a axis
    const Real f2n = na ? fa2 : ft2;
    Real f2g = na ? ft2 : fa2;
    const int ne2 = (f2n > 0.5*nm1) ? (nc - 1) : 0;
    f2g = fmin(fmax(f2g, 0.0), nm1);
    int g2 = static_cast<int>(floor(f2g));
    g2 = (g2 > nc - 2) ? (nc - 2) : g2;
    const Real u2 = f2g - g2;
    for (int q = 0; q < 2; ++q) {
      qc[2+q] = c2p;
      qa[2+q] = na ? ne2 : (g2 + q);
      qt[2+q] = na ? (g2 + q) : ne2;
      qw[2+q] = dist*(q ? u2 : (1.0 - u2));
    }
    nq = 4;
  } else {
    // the vertex corner: the three cells around the cube vertex, barycentric
    const int ae = (fa > nm1) ? (nc - 1) : 0, te = (ft > nm1) ? (nc - 1) : 0;
    const int ca = nb(c, (fa > nm1) ? 1 : 0), ct = nb(c, (ft > nm1) ? 3 : 2);
    Real fa_[3], fb_[3], n0[3], n1[3], n2[3], vv[3];
    cubed_sphere::PanelFrame(c, fa_, fb_, n0);
    cubed_sphere::PanelFrame(ca, fa_, fb_, n1);
    cubed_sphere::PanelFrame(ct, fa_, fb_, n2);
    for (int e = 0; e < 3; ++e) {vv[e] = n0[e] + n1[e] + n2[e];}
    Real xv, ev;
    cubed_sphere::CartToPanel(ca, vv, xv, ev);
    const int aa = (xv > 0.0) ? (nc - 1) : 0, ta = (ev > 0.0) ? (nc - 1) : 0;
    cubed_sphere::CartToPanel(ct, vv, xv, ev);
    const int ab = (xv > 0.0) ? (nc - 1) : 0, tb = (ev > 0.0) ? (nc - 1) : 0;
    Real c1[3], c2c[3], c3c[3];
    cubed_sphere::PanelToCart(c, -q4 + (ae + 0.5)*h, -q4 + (te + 0.5)*h, c1);
    cubed_sphere::PanelToCart(ca, -q4 + (aa + 0.5)*h, -q4 + (ta + 0.5)*h, c2c);
    cubed_sphere::PanelToCart(ct, -q4 + (ab + 0.5)*h, -q4 + (tb + 0.5)*h, c3c);
    const Real dd = VgdCsDet3(c1, c2c, c3c);
    Real l1 = VgdCsDet3(x, c2c, c3c)/dd;
    Real l2 = VgdCsDet3(c1, x, c3c)/dd;
    Real l3 = VgdCsDet3(c1, c2c, x)/dd;
    l1 = fmax(l1, 0.0);
    l2 = fmax(l2, 0.0);
    l3 = fmax(l3, 0.0);
    const Real ls = l1 + l2 + l3;
    if (ls > 0.0) {
      l1 /= ls;
      l2 /= ls;
      l3 /= ls;
    } else {
      l1 = 1.0;
    }
    qc[0] = c; qa[0] = ae; qt[0] = te; qw[0] = l1;
    qc[1] = ca; qa[1] = aa; qt[1] = ta; qw[1] = l2;
    qc[2] = ct; qa[2] = ab; qt[2] = tb; qw[2] = l3;
    nq = 3;
  }
  // the global cells -> band positions of this block
  bool ok = true;
  Real ws = 0.0;
  st.n = 0;
  for (int q = 0; q < nq; ++q) {
    int jg, kg;
    const int cq = qc[q];
    if (cq == p) {
      jg = qa[q];
      kg = qt[q];
    } else if (cm(p,cq,0) == 0) {
      ok = false;
      continue;
    } else {
      jg = cm(p,cq,1)*qa[q] + cm(p,cq,2)*qt[q] + cm(p,cq,3);
      kg = cm(p,cq,4)*qa[q] + cm(p,cq,5)*qt[q] + cm(p,cq,6);
    }
    const int jj = jg - j0 + wb, kk = kg - k0 + wb;
    if (jj < 0 || jj >= c2 || kk < 0 || kk >= c3) {
      ok = false;
      continue;
    }
    st.k[st.n] = kk;
    st.j[st.n] = jj;
    st.w[st.n] = qw[q];
    ws += qw[q];
    ++st.n;
  }
  if (!ok && ws > 0.0) {
    for (int q = 0; q < st.n; ++q) {st.w[q] /= ws;}
  }
  return ok;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsInit
//! \brief the cs band: depth vgd_w (from the narrowest cell), the panel tables, the
//! unit vector of every band cell, the exact exchange lists (from VetGdInit)

void RadiationM1::VetGdCsInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  const int nc = pm->mesh_indcs.nx2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  if (nc < 2 || pm->mesh_indcs.nx3 != nc || (nc % nx2) != 0 || (nc % nx3) != 0) {
    VgdCsFatal("needs mesh nx2 = nx3 = N >= 2 per panel, divided by the MeshBlock sizes");
  }
  if (pm->multilevel) {VgdCsFatal("SMR/AMR not supported");}
  vgd_csn = nc;
  const double h = 0.5*M_PI/nc;
  const double q4 = 0.25*M_PI;
  const int nbj = nc/nx2, nbk = nc/nx3;
  const int ngb = pm->nmb_total;
  // the blocks: panel and global first cell; the band of a block comes from the blocks
  // of its own x1 level
  std::vector<int> bp(ngb), bj(ngb), bk(ngb), gof(vgd_nb1*6*nbj*nbk, -1);
  auto gofi = [&](const int g, const int p, const int b3, const int b2) -> int & {
    return gof[((pm->lloc_eachmb[g].lx1*6 + p)*nbk + b3)*nbj + b2];
  };
  for (int g = 0; g < ngb; ++g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    if (l.panel < 0 || l.panel > 5 || l.lx2 < 0 || l.lx2 >= nbj || l.lx3 < 0 ||
        l.lx3 >= nbk) {
      VgdCsFatal("unexpected MeshBlock location (panel / lx2 / lx3)");
    }
    bp[g] = l.panel;
    bj[g] = l.lx2*nx2;
    bk[g] = l.lx3*nx3;
    gofi(g, l.panel, l.lx3, l.lx2) = g;
  }
  for (size_t t = 0; t < gof.size(); ++t) {
    if (gof[t] < 0) {VgdCsFatal("incomplete MeshBlock grid");}
  }
  {
    // the local blocks' sizes agree with the logical locations
    auto &mbs = pmy_pack->pmb->mb_size;
    for (int m = 0; m < nmb; ++m) {
      const double x2 = -1.0 + 2.0*bj[g0 + m]/nc, x3 = -1.0 + 2.0*bk[g0 + m]/nc;
      if (std::fabs(mbs.h_view(m).x2min - x2) > 1.0e-12 ||
          std::fabs(mbs.h_view(m).x3min - x3) > 1.0e-12) {
        VgdCsFatal("MeshBlock size does not match its logical location");
      }
    }
  }
  // panel frames and the neighbour across each side (0 -x2, 1 +x2, 2 -x3, 3 +x3)
  double fa[6][3], fb[6][3], fn[6][3];
  for (int p = 0; p < 6; ++p) {
    Real a[3], b[3], n[3];
    cubed_sphere::PanelFrame(p, a, b, n);
    for (int e = 0; e < 3; ++e) {fa[p][e] = a[e]; fb[p][e] = b[e]; fn[p][e] = n[e];}
  }
  int nbr[6][4];
  for (int p = 0; p < 6; ++p) {
    for (int s = 0; s < 4; ++s) {
      double t[3];
      for (int e = 0; e < 3; ++e) {
        const double v = (s < 2) ? fa[p][e] : fb[p][e];
        t[e] = (s % 2 == 0) ? -v : v;
      }
      nbr[p][s] = -1;
      for (int q = 0; q < 6; ++q) {
        if (t[0]*fn[q][0] + t[1]*fn[q][1] + t[2]*fn[q][2] > 0.5) {nbr[p][s] = q;}
      }
      if (nbr[p][s] < 0) {VgdCsFatal("panel neighbour table");}
    }
  }
  auto cen = [&](const int p, const int a, const int t, double c[3]) {
    Real q[3];
    cubed_sphere::PanelToCart(p, -q4 + (a + 0.5)*h, -q4 + (t + 0.5)*h, q);
    for (int e = 0; e < 3; ++e) {c[e] = q[e];}
  };
  // the unfolded position (J, K) of panel p -> the global cell (q, a, t) it holds
  auto unfold = [&](const int p, const int J, const int K, int &q, int &a,
                    int &t) -> bool {
    const bool ij = (J >= 0 && J < nc), ik = (K >= 0 && K < nc);
    if (ij && ik) {q = p; a = J; t = K; return true;}
    if (!ij && !ik) {return false;}
    int s, jm = J, km = K;
    if (!ij) {
      s = (J < 0) ? 0 : 1;
      jm = (J < 0) ? (-1 - J) : (2*nc - 1 - J);
    } else {
      s = (K < 0) ? 2 : 3;
      km = (K < 0) ? (-1 - K) : (2*nc - 1 - K);
    }
    if (jm < 0 || jm >= nc || km < 0 || km >= nc) {return false;}
    q = nbr[p][s];
    double c[3], e[3];
    cen(p, jm, km, c);
    for (int u = 0; u < 3; ++u) {e[u] = (fn[p][u] - fn[q][u])/std::sqrt(2.0);}
    const double dc = c[0]*e[0] + c[1]*e[1] + c[2]*e[2];
    Real cr[3];
    for (int u = 0; u < 3; ++u) {cr[u] = c[u] - 2.0*dc*e[u];}
    Real xi, eta;
    cubed_sphere::CartToPanel(q, cr, xi, eta);
    const double ga = (xi + q4)/h - 0.5, gt = (eta + q4)/h - 0.5;
    a = static_cast<int>(std::lround(ga));
    t = static_cast<int>(std::lround(gt));
    if (std::fabs(ga - a) > 1.0e-6 || std::fabs(gt - t) > 1.0e-6 || a < 0 || a >= nc ||
        t < 0 || t >= nc) {
      VgdCsFatal("unfold: the mirror cell is off the neighbour's lattice");
    }
    return true;
  };
  // the inverse maps (q's (a, t) -> p's unfolded (J, K)), a signed permutation + shift
  std::vector<int> cmh(6*6*7, 0);
  auto CM = [&](const int p, const int q, const int c) -> int & {
    return cmh[(p*6 + q)*7 + c];
  };
  for (int p = 0; p < 6; ++p) {
    CM(p,p,0) = 1; CM(p,p,1) = 1; CM(p,p,5) = 1;
    for (int s = 0; s < 4; ++s) {
      const int jb = (s == 0) ? -1 : ((s == 1) ? nc : 0);
      const int kb = (s == 2) ? -1 : ((s == 3) ? nc : 0);
      const int sj = (s == 0) ? -1 : 1, sk = (s == 2) ? -1 : 1;
      int q0, a0, t0, q1, a1, t1, q2, a2, t2;
      const bool o0 = unfold(p, jb, kb, q0, a0, t0);
      const bool o1 = unfold(p, jb + ((s < 2) ? sj : 1), kb, q1, a1, t1);
      const bool o2 = unfold(p, jb, kb + ((s < 2) ? 1 : sk), q2, a2, t2);
      if (!o0 || !o1 || !o2 || q1 != q0 || q2 != q0) {VgdCsFatal("unfold strip");}
      const int stj = (s < 2) ? sj : 1, stk = (s < 2) ? 1 : sk;
      const int uj0 = (a1 - a0)*stj, uj1 = (t1 - t0)*stj;   // (a, t) per unit J
      const int uk0 = (a2 - a0)*stk, uk1 = (t2 - t0)*stk;   // (a, t) per unit K
      if (std::abs(uj0) + std::abs(uj1) != 1 || std::abs(uk0) + std::abs(uk1) != 1 ||
          uj0*uk0 + uj1*uk1 != 0) {
        VgdCsFatal("unfold: not a signed permutation");
      }
      CM(p,q0,0) = 1;
      CM(p,q0,1) = uj0;
      CM(p,q0,2) = uj1;
      CM(p,q0,3) = jb - uj0*a0 - uj1*t0;
      CM(p,q0,4) = uk0;
      CM(p,q0,5) = uk1;
      CM(p,q0,6) = kb - uk0*a0 - uk1*t0;
    }
  }
  // check: the forward unfold and the inverse map agree on every strip cell
  for (int p = 0; p < 6; ++p) {
    for (int s = 0; s < 4; ++s) {
      for (int d = 0; d < nc; ++d) {
        for (int g = 0; g < nc; ++g) {
          const int J = (s == 0) ? (-1 - d) : ((s == 1) ? (nc + d) : g);
          const int K = (s == 2) ? (-1 - d) : ((s == 3) ? (nc + d) : g);
          int q, a, t;
          if (!unfold(p, J, K, q, a, t)) {VgdCsFatal("unfold check");}
          const int jj = CM(p,q,1)*a + CM(p,q,2)*t + CM(p,q,3);
          const int kk = CM(p,q,4)*a + CM(p,q,5)*t + CM(p,q,6);
          if (jj != J || kk != K) {VgdCsFatal("unfold: inverse map mismatch");}
        }
      }
    }
  }
  // the band depth: the sp formula with the narrowest angular cell spacing of the panel
  // (+3: one more than sp for the unfolded index kink across a seam); capped at N
  double dmin = 1.0e300;
  for (int a = 0; a < nc; ++a) {
    for (int t = 0; t < nc; ++t) {
      double c0[3], c1[3];
      cen(0, a, t, c0);
      if (a + 1 < nc) {
        cen(0, a + 1, t, c1);
        dmin = std::min(dmin, std::acos(std::min(1.0, c0[0]*c1[0] + c0[1]*c1[1] +
                                                         c0[2]*c1[2])));
      }
      if (t + 1 < nc) {
        cen(0, a, t + 1, c1);
        dmin = std::min(dmin, std::acos(std::min(1.0, c0[0]*c1[0] + c0[1]*c1[1] +
                                                         c0[2]*c1[2])));
      }
    }
  }
  {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    vgd_wsh.assign(indcs.nx1 + 2*indcs.ng, 0);
    vgd_wsi.assign(indcs.nx1 + 2*indcs.ng, 0);
    vgd_wso.assign(indcs.nx1 + 2*indcs.ng, 0);
    int wneed = 0;
    auto xf_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_x1f);
    for (int m = 0; m < nmb; ++m) {
      // (the shells next to another level: its radii from the x1 ghosts)
      const bool lo = (xf_h(m) & 1), hi = (xf_h(m) & 2);
      for (int l = indcs.is; l <= indcs.ie; ++l) {
        const double r = x1v_h(m,l);
        const double rlo = (l > indcs.is || !lo) ? x1v_h(m,l-1) : x1f_h(m,l);
        double a = 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
        if (l < indcs.ie || !hi) {
          a = std::max(a, std::sqrt(x1v_h(m,l+1)*x1v_h(m,l+1) - r*r));
        }
        wneed = std::max(wneed, static_cast<int>(std::ceil(a/(rlo*dmin))) + 3);
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &wneed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    vgd_w = std::min(wneed, nc);
    vgd_capped = (vgd_w < wneed);
    for (int l = indcs.is; l <= indcs.ie; ++l) {
      vgd_wsh[l] = vgd_w;
      vgd_wsi[l] = vgd_w;
      vgd_wso[l] = vgd_w;
    }
    if (me == 0) {
      std::cout << "<rad_m1> vet_gd (cubed sphere): N = " << nc << " per panel, "
                << "narrowest "
                << "cell spacing " << dmin << " rad, dense lateral band " << vgd_w
                << " cells (needed " << wneed << ")"
                << (vgd_capped ? ": capped at N, near-tangent reads beyond are clamped "
                                 "(counted)" : "") << std::endl;
    }
  }
  const int w = vgd_w;
  const int c2w = nx2 + 2*w, c3w = nx3 + 2*w;
  // the device tables
  Kokkos::realloc(vgd_cscm, 6, 6, 7);
  Kokkos::realloc(vgd_csnb, 6, 4);
  Kokkos::realloc(vgd_csblk, nmb, 3);
  {
    auto cm_h = Kokkos::create_mirror_view(vgd_cscm);
    auto nb_h = Kokkos::create_mirror_view(vgd_csnb);
    auto bl_h = Kokkos::create_mirror_view(vgd_csblk);
    for (int p = 0; p < 6; ++p) {
      for (int q = 0; q < 6; ++q) {
        for (int c = 0; c < 7; ++c) {cm_h(p,q,c) = CM(p,q,c);}
      }
      for (int s = 0; s < 4; ++s) {nb_h(p,s) = nbr[p][s];}
    }
    for (int m = 0; m < nmb; ++m) {
      bl_h(m,0) = bp[g0 + m];
      bl_h(m,1) = bj[g0 + m];
      bl_h(m,2) = bk[g0 + m];
    }
    Kokkos::deep_copy(vgd_cscm, cm_h);
    Kokkos::deep_copy(vgd_csnb, nb_h);
    Kokkos::deep_copy(vgd_csblk, bl_h);
  }
  // the unit vector of the cell at every band position (0: none)
  Kokkos::realloc(vgd_csr, nmb, c3w, c2w, 3);
  vgd_csr_h.assign(static_cast<size_t>(nmb)*c3w*c2w*3, 0.0);
  {
    auto cr_h = Kokkos::create_mirror_view(vgd_csr);
    for (int m = 0; m < nmb; ++m) {
      const int g = g0 + m;
      for (int kk = 0; kk < c3w; ++kk) {
        for (int jj = 0; jj < c2w; ++jj) {
          int q, a, t;
          double c[3] = {0.0, 0.0, 0.0};
          if (unfold(bp[g], bj[g] + jj - w, bk[g] + kk - w, q, a, t)) {cen(q, a, t, c);}
          for (int e = 0; e < 3; ++e) {
            cr_h(m,kk,jj,e) = c[e];
            vgd_csr_h[((static_cast<size_t>(m)*c3w + kk)*c2w + jj)*3 + e] = c[e];
          }
        }
      }
    }
    Kokkos::deep_copy(vgd_csr, cr_h);
  }
  // the exchange lists: every receiver block in global order, its band positions in
  // (k, j) order; the same iteration on every rank, so a sender's list to rank R and
  // R's receive list from it are in the same order
  std::vector<int> loc;
  std::map<int, std::vector<int>> snd, rcv;
  for (int gr = 0; gr < ngb; ++gr) {
    const int rr = pm->rank_eachmb[gr];
    for (int kk = 0; kk < c3w; ++kk) {
      for (int jj = 0; jj < c2w; ++jj) {
        if (jj >= w && jj < w + nx2 && kk >= w && kk < w + nx3) {continue;}
        int q, a, t;
        if (!unfold(bp[gr], bj[gr] + jj - w, bk[gr] + kk - w, q, a, t)) {continue;}
        const int go = gofi(gr, q, t/nx3, a/nx2);
        const int ro = pm->rank_eachmb[go];
        if (rr != me && ro != me) {continue;}
        const int ko = t - bk[go] + w, jo = a - bj[go] + w;
        if (rr == me && ro == me) {
          loc.insert(loc.end(), {gr - g0, kk, jj, go - g0, ko, jo});
        } else if (rr == me) {
          rcv[ro].insert(rcv[ro].end(), {gr - g0, kk, jj});
        } else {
          snd[rr].insert(snd[rr].end(), {go - pm->gids_eachrank[me], ko, jo});
        }
      }
    }
  }
  vgd_csnloc = static_cast<int>(loc.size()/6);
  // the owner of every band position for the sweep's direct on-rank reads: the own
  // interior and on-rank owners' interiors, else (remote, none) the position itself
  Kokkos::realloc(vgd_csown, nmb, c3w, c2w, 3);
  {
    auto o_h = Kokkos::create_mirror_view(vgd_csown);
    for (int m = 0; m < nmb; ++m) {
      for (int kk = 0; kk < c3w; ++kk) {
        for (int jj = 0; jj < c2w; ++jj) {
          o_h(m,kk,jj,0) = m;
          o_h(m,kk,jj,1) = kk;
          o_h(m,kk,jj,2) = jj;
        }
      }
    }
    for (int e = 0; e < vgd_csnloc; ++e) {
      const int *x = &loc[6*e];
      o_h(x[0],x[1],x[2],0) = x[3];
      o_h(x[0],x[1],x[2],1) = x[4];
      o_h(x[0],x[1],x[2],2) = x[5];
    }
    Kokkos::deep_copy(vgd_csown, o_h);
  }
  Kokkos::realloc(vgd_csloc, std::max(vgd_csnloc, 1), 6);
  {
    auto l_h = Kokkos::create_mirror_view(vgd_csloc);
    for (int e = 0; e < vgd_csnloc; ++e) {
      for (int c = 0; c < 6; ++c) {l_h(e,c) = loc[6*e + c];}
    }
    Kokkos::deep_copy(vgd_csloc, l_h);
  }
  vgd_csprk.clear();
  for (const auto &x : snd) {vgd_csprk.push_back(x.first);}
  for (const auto &x : rcv) {
    if (snd.find(x.first) == snd.end()) {vgd_csprk.push_back(x.first);}
  }
  std::sort(vgd_csprk.begin(), vgd_csprk.end());
  const int np = static_cast<int>(vgd_csprk.size());
  vgd_csps.assign(np, 0);
  vgd_cspn.assign(np, 0);
  vgd_cspr.assign(np, 0);
  vgd_cspm.assign(np, 0);
  std::vector<int> sl, rl;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    vgd_csps[p] = static_cast<int>(sl.size()/3);
    vgd_cspr[p] = static_cast<int>(rl.size()/3);
    if (snd.count(rk)) {sl.insert(sl.end(), snd[rk].begin(), snd[rk].end());}
    if (rcv.count(rk)) {rl.insert(rl.end(), rcv[rk].begin(), rcv[rk].end());}
    vgd_cspn[p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
    vgd_cspm[p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
  }
  vgd_csnsnd = static_cast<int>(sl.size()/3);
  vgd_csnrcv = static_cast<int>(rl.size()/3);
  Kokkos::realloc(vgd_cssl, std::max(vgd_csnsnd, 1), 3);
  Kokkos::realloc(vgd_csrl, std::max(vgd_csnrcv, 1), 3);
  {
    auto s_h = Kokkos::create_mirror_view(vgd_cssl);
    auto r_h = Kokkos::create_mirror_view(vgd_csrl);
    for (int e = 0; e < vgd_csnsnd; ++e) {
      for (int c = 0; c < 3; ++c) {s_h(e,c) = sl[3*e + c];}
    }
    for (int e = 0; e < vgd_csnrcv; ++e) {
      for (int c = 0; c < 3; ++c) {r_h(e,c) = rl[3*e + c];}
    }
    Kokkos::deep_copy(vgd_cssl, s_h);
    Kokkos::deep_copy(vgd_csrl, r_h);
  }
  // message buffers: one shell of every direction (the source exchange is chunked)
  const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
  const int nvb = std::max(vgd_n, 2*8);
  Kokkos::realloc(vgd_cssb, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_csrb, static_cast<size_t>(nent)*nvb);
  vgd_hmpi = (np > 0);
  // every side of every block carries the band (VgdRagAlloc): no direct neighbour reads
  Kokkos::realloc(vgd_hloc, 8*nmb);
  Kokkos::deep_copy(vgd_hloc, -1);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  // the compact lists (vet_gd_halo_compact > 0 with an MPI band): the
  // receive entry of every band position (-1: own interior, on-rank owner or no cell)
  vgd_cscomp = (vgd_hcomp > 0) && vgd_hmpi;
  Kokkos::realloc(vgd_csrix, nmb, c3w, c2w);
  {
    auto x_h = Kokkos::create_mirror_view(vgd_csrix);
    for (size_t q = 0; q < x_h.size(); ++q) {x_h.data()[q] = -1;}
    for (int e = 0; e < vgd_csnrcv; ++e) {x_h(rl[3*e], rl[3*e+1], rl[3*e+2]) = e;}
    Kokkos::deep_copy(vgd_csrix, x_h);
  }
  Kokkos::realloc(vgd_csmiss, 1);
  // one cache entry per set of vet_rotate_avg (the K sets cycle; 2 before:
  // every rotation rebuilt its lists); the entries beyond vgd_csncl stay unallocated
  vgd_csncl = 2;
  if (vgd_cscomp) {
    vgd_csncl = std::min(std::max(2, vra_k), static_cast<int>(M1_VGD_CSCL_MAX));
  }
  vgd_csper = static_cast<int64_t>(vgd_n)*indcs.nx1;
  if (vgd_cspad) {vgd_csper = 16*((vgd_csper + 15)/16);}
  for (int c = 0; c < vgd_csncl; ++c) {
    vgd_cscl[c].alpha = -9.0;
    vgd_cscl[c].use = 0;
    vgd_cscl[c].msk = decltype(vgd_cscl[c].msk)();   // (hgmem: VetGdCsCompact, vgd_csmb)
  }
  if (me == 0) {
    std::cout << "<rad_m1> vet_gd (cubed sphere): rank 0 band copies on-rank "
              << vgd_csnloc << ", sent " << vgd_csnsnd << ", received " << vgd_csnrcv
              << " (cells per shell), partners " << np << "; MPI band "
              << (vgd_cscomp ? "COMPACT (vet_gd_halo_compact > 0: per shell and pass "
                               "only the read (cell, direction) values)"
                             : (vgd_hmpi ? "dense (every remote band cell, every "
                                           "direction, per shell and pass)"
                                         : "none"))
              << "; list cache " << vgd_csncl << " sets" << std::endl;
  }
#if M1_VGDCS_NCCL
  // the NCCL communicator now (collective: every rank), not at the first
  // exchange (a rank without partners never exchanges)
  if (global_variable::nranks > 1) {
    if (VgdNcclOn() && me == 0) {
      std::cout << "<rad_m1> csnccl-1031: vet_gd cs compact exchange over NCCL "
                << "(stream-ordered, no host fence per shell)" << std::endl;
    }
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsCompact
//! \brief the compact lists of the current direction set (vgd_alpha), from
//! an LRU cache of max(2, vet_rotate_avg) entries (was 2: the current and
//! the next set).  (1) the read mask:
//! a kernel replays the sweep's upwind geometry (VetGdCsSweep, term for term) for every
//! interior (cell, direction) and shell, and marks every REMOTE band value (receive
//! entry, direction, upwind shell) that the sweep reads -- for every first shell at once
//! (the reads of shell l depend on the cut only when l is the cut: both variants are
//! marked), so the lists depend on the set alone; value 1 = an inward value (sent after
//! pass 0 of that shell), 2 = outward (after pass 1); (2) every receiver sends its marks
//! to the owners (one message per partner, once per set); (3) both sides build their
//! lists in the same (pass, shell, partner, entry, direction) order.

void RadiationM1::VetGdCsCompact() {
  int c = -1;
  for (int q = 0; q < vgd_csncl; ++q) {
    if (vgd_cscl[q].alpha == vgd_alpha && vgd_cscl[q].ng == vgd_csxg) {c = q;}
  }
  if (c >= 0 && vgd_cscl[c].msk.data() != nullptr) {
    vgd_cscur = c;
    vgd_cscl[c].use = ++vgd_csuse;
    return;
  }
  // (hgmem: a cached set whose read mask went to another set gets its mask again, (1)
  // only: its lists stand, the mask is the sweep's check of them)
  const bool remask = (c >= 0);
  Kokkos::Timer tb;
  if (!remask) {
    c = 0;   // the least recently used entry (the first of equal ones)
    for (int q = 1; q < vgd_csncl; ++q) {
      if (vgd_cscl[q].use < vgd_cscl[c].use) {c = q;}
    }
  }
  VgdCsCl &E = vgd_cscl[c];
  // hgmem: two read masks only (vgd_csmb), the one of the entry in use (xfade's current
  // set while the next set's is made) is kept
  if (E.msk.data() == nullptr) {
    const auto *k = vgd_cscl[vgd_cscur].msk.data();
    const int b = (k != nullptr && k == vgd_csmb[0].data()) ? 1 : 0;
    for (auto &e : vgd_cscl) {
      if (e.msk.data() == vgd_csmb[b].data()) {e.msk = decltype(e.msk)();}
    }
    const int64_t nmv = static_cast<int64_t>(std::max(vgd_csnrcv, 1))*vgd_csper;
    if (vgd_csmb[b].extent(0) != static_cast<size_t>((nmv + 15)/16)) {
      vgd_csmb[b] = decltype(E.msk)();
      Kokkos::realloc(vgd_csmb[b], (nmv + 15)/16);
    }
    E.msk = vgd_csmb[b];
  }
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int c2b = indcs.nx2 + 2*wb, c3b = indcs.nx3 + 2*wb;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  auto msk_ = E.msk;
  Kokkos::deep_copy(msk_, static_cast<uint32_t>(0));
  auto dir_ = vgd_dir;
  auto mr_ = vgd_mr;
  auto cr_ = vgd_csr;
  auto bl_ = vgd_csblk;
  auto cm_ = vgd_cscm;
  auto nb_ = vgd_csnb;
  auto rx_ = vgd_csrix;
  const int nc = vgd_csn;
  const Real h = 0.5*M_PI/nc;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const int64_t per_ = vgd_csper;   // n*n1, or padded to 16
  Kokkos::Timer tsub;
  // (1) the read mask (VetGdCsSweep's geometry; the values never matter)
  par_for("m1_vgdcs_rmask", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, n - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
    const Real mr = mr_(m,k,j,d);
    const bool inw = !(mr >= 0.0);
    const Real rx = cr_(m,k,j,0), ry = cr_(m,k,j,1), rz = cr_(m,k,j,2);
    // mark the reads of one upwind segment (ds, upwind shell iu, branch want_in)
    auto mark = [&](const Real r, const Real ds, const int iu, const bool want_in) {
      Real xu[3] = {r*rx - ds*nx, r*ry - ds*ny, r*rz - ds*nz};
      const Real ru = sqrt(xu[0]*xu[0] + xu[1]*xu[1] + xu[2]*xu[2]);
      xu[0] /= ru;
      xu[1] /= ru;
      xu[2] /= ru;
      VgdCsSt st;
      VgdCsStencil(xu, bl_(m,0), bl_(m,1), bl_(m,2), wb, c2b, c3b, nc, h, cm_, nb_, st);
      if (st.n == 0) {return;}   // the own column (interior)
      bool ok[4];
      bool all = true;
      for (int e = 0; e < st.n; ++e) {
        const Real mq = mr_(m,st.k[e],st.j[e],d);
        ok[e] = want_in ? (mq < 0.0) : (mq >= 0.0);
        all = all && ok[e];
      }
      for (int e = 0; e < st.n; ++e) {
        if (!all && !(ok[e] && st.w[e] > 0.0)) {continue;}
        const int x = rx_(m,st.k[e],st.j[e]);
        if (x >= 0) {
          const int64_t t = static_cast<int64_t>(x)*per_ + d*n1 + (iu - is);
          Kokkos::atomic_or(&msk_(t >> 4), (want_in ? 1u : 2u) << (2*(t & 15)));
        }
      }
    };
    for (int l = 0; l < n1; ++l) {
      const int i = is + l;
      const Real r = cx1v(m,i);
      const Real zr = r*fabs(mr);
      const Real p2 = fmax(r*r - zr*zr, 0.0);
      if (inw) {
        if (l == n1 - 1) {continue;}
        const Real ru = cx1v(m,i+1);
        mark(r, sqrt(ru*ru - p2) - zr, i + 1, true);
        continue;
      }
      if (zr < 1.0e-12*r) {continue;}
      // as the shell above the cut (l > lcut)
      if (l > 0) {
        const Real rd = cx1v(m,i-1);
        if (p2 <= rd*rd*(1.0 + 1.0e-13)) {
          mark(r, zr - sqrt(fmax(rd*rd - p2, 0.0)), i - 1, false);
        } else {
          mark(r, 2.0*zr, i, true);
        }
      }
      // as the cut shell (l == lcut): the inner face (no intensity read) or turned
      const Real rf = cx1f(m,i);
      if (!(p2 <= rf*rf)) {mark(r, 2.0*zr, i, true);}
    }
  });
  Kokkos::fence();
  vgd_cstb_k += tsub.seconds();
  tsub.reset();
  if (remask) {
    vgd_cscur = c;
    E.use = ++vgd_csuse;
    return;
  }
  const int np = static_cast<int>(vgd_csprk.size());
  const size_t per = static_cast<size_t>(per_);
  // (3) the lists: counts per (pass, shell, group, partner), offsets, fill in (entry, d)
  // order (G = vgd_csxg contiguous direction groups; G = 1 as before)
  const int ng = vgd_csxg;
  std::vector<int> dg(n);
  for (int d = 0; d < n; ++d) {dg[d] = (d*ng)/n;}   // group of d: d in [g n/G, (g+1) n/G)
  const size_t nsl = static_cast<size_t>(2)*n1*ng;
  std::vector<int> sl_h, rl_h;
  int64_t ts = 0, tr = 0;
  if (per % 16 == 0) {
    // the lists straight from the packed 2-bit words (partner blocks of
    // entries are word aligned when n*n1 % 16 == 0): zero words skipped, the marks sent
    // as packed words (8x fewer bytes); the same lists in the same (t-increasing) order
    // as the byte path below
    const size_t wpe = per/16;   // words per entry
    // the host copies of the marks in PERSISTENT buffers (grown, never
    // freed), not a fresh mirror / vector per build.  RSG E 10-31: under an LD_PRELOAD
    // profiler (nsys osrt / cpu sampling) UCX's registration cache missed the munmap of
    // the previous build's buffer, the new buffer reused its addresses and the rendezvous
    // sent stale pages (sender / receiver hashes differ): inconsistent lists, then
    // MPI_ERR_TRUNCATE or a hang in the first sweep.  Same values, same lists (bitwise)
    const size_t nmw = msk_.extent(0);
    if (vgd_csmw_h.size() < nmw) {vgd_csmw_h.resize(nmw);}
    {
      Kokkos::View<uint32_t *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          mh(vgd_csmw_h.data(), nmw);
      Kokkos::deep_copy(mh, msk_);
    }
    const uint32_t *mp = vgd_csmw_h.data();
    const size_t nwr = static_cast<size_t>(vgd_csnrcv)*wpe;
    bool both = false;
    for (size_t w = 0; w < nwr; ++w) {
      const uint32_t x = mp[w];
      if (x != 0u && ((x & (x >> 1) & 0x55555555u) != 0u)) {both = true;}
    }
    if (both) {VgdCsFatal("read mask: a value marked both inward and outward");}
    vgd_cstb_d += tsub.seconds();
    tsub.reset();
    const size_t nsw = std::max<size_t>(static_cast<size_t>(vgd_csnsnd)*wpe, 1);
    if (vgd_csmw_r.size() < nsw) {vgd_csmw_r.resize(nsw);}
    std::fill(vgd_csmw_r.begin(), vgd_csmw_r.begin() + nsw, 0u);
    uint32_t *smw = vgd_csmw_r.data();
#if MPI_PARALLEL_ENABLED
    {
      // own tags (7107 marks, 7108 their checksums; 7105 is the fused
      // twin's sweep exchange) and a checksum per message: a block that arrives other
      // than it was sent is FATAL here, never a silent wrong list
      auto hsum = [](const uint32_t *a, const size_t nw) {
        uint64_t h = 1469598103934665603ull;
        for (size_t w = 0; w < nw; ++w) {h = (h ^ a[w])*1099511628211ull;}
        return h;
      };
      std::vector<uint64_t> hs(np, 0), hr(np, 0);
      std::vector<MPI_Request> req;
      for (int p = 0; p < np; ++p) {
        const int rk = vgd_csprk[p];
        if (vgd_cspn[p] > 0) {
          req.emplace_back();
          MPI_Irecv(smw + static_cast<size_t>(vgd_csps[p])*wpe,
                    static_cast<int>(vgd_cspn[p]*wpe), MPI_UINT32_T, rk, 7107,
                    MPI_COMM_WORLD, &req.back());
          req.emplace_back();
          MPI_Irecv(&hr[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
        }
        if (vgd_cspm[p] > 0) {
          const uint32_t *b = mp + static_cast<size_t>(vgd_cspr[p])*wpe;
          const size_t nw = static_cast<size_t>(vgd_cspm[p])*wpe;
          hs[p] = hsum(b, nw);
          req.emplace_back();
          MPI_Isend(b, static_cast<int>(nw), MPI_UINT32_T, rk, 7107, MPI_COMM_WORLD,
                    &req.back());
          req.emplace_back();
          MPI_Isend(&hs[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
      for (int p = 0; p < np; ++p) {
        if (vgd_cspn[p] <= 0) {continue;}
        const uint64_t h = hsum(smw + static_cast<size_t>(vgd_csps[p])*wpe,
                                static_cast<size_t>(vgd_cspn[p])*wpe);
        if (h != hr[p]) {
          VgdCsFatal("compact lists: the read marks from rank " +
                     std::to_string(vgd_csprk[p]) + " arrived corrupted (rank " +
                     std::to_string(global_variable::my_rank) + ", checksum mismatch): "
                     "the MPI transport delivered other bytes than were sent (seen with "
                     "UCX's registration cache under an LD_PRELOAD profiler such as nsys "
                     "-t osrt / -s cpu: export UCX_RCACHE_ENABLE=n for such runs)");
        }
      }
    }
#endif
    vgd_cstb_m += tsub.seconds();
    tsub.reset();
    // (each partner's counts and fills on its own host thread: disjoint slots and list
    // ranges, the same lists; h200spd: the host lists were ~1.1 s of each ~2 s build)
    auto parts = [np](const auto &f) {
      std::vector<std::thread> th;
      for (int p = 1; p < np; ++p) th.emplace_back(f, p);
      if (np > 0) f(0);
      for (auto &t : th) t.join();
    };
    auto buildw = [&](const uint32_t *mk, const std::vector<int> &st0,
                      const std::vector<int> &cnt, std::vector<int64_t> &off,
                      std::vector<int> &lst) {
      std::vector<int64_t> nn(nsl*np, 0);
      std::vector<std::vector<size_t>> nzw(np);   // the nonzero words of each partner
      parts([&](const int p) {
        const size_t w0 = static_cast<size_t>(st0[p])*wpe;
        const size_t w1 = w0 + static_cast<size_t>(cnt[p])*wpe;
        for (size_t w = w0; w < w1; ++w) {
          uint32_t x = mk[w];
          if (x == 0u) {continue;}
          nzw[p].push_back(w);
          const size_t r0 = (w - static_cast<size_t>(st0[p])*wpe)*16;   // in-block t
          while (x != 0u) {
            const int b = __builtin_ctz(x) >> 1;
            const uint32_t v = (x >> (2*b)) & 3u;
            x &= ~(3u << (2*b));
            const size_t r = (r0 + b) % per;
            const int d = static_cast<int>(r/n1), sh = static_cast<int>(r % n1);
            nn[((static_cast<size_t>(v - 1)*n1 + sh)*ng + dg[d])*np + p] += 1;
          }
        }
      });
      off.assign(nsl*(np + 1), 0);
      int64_t t = 0;
      for (size_t sl = 0; sl < nsl; ++sl) {
        for (int p = 0; p < np; ++p) {
          off[sl*(np + 1) + p] = t;
          t += nn[sl*np + p];
        }
        off[sl*(np + 1) + np] = t;
      }
      lst.assign(std::max<int64_t>(t, 1), 0);
      std::vector<int64_t> cur(nsl*np);
      for (size_t sl = 0; sl < nsl; ++sl) {
        for (int p = 0; p < np; ++p) {cur[sl*np + p] = off[sl*(np + 1) + p];}
      }
      parts([&](const int p) {
        const size_t wb0 = static_cast<size_t>(st0[p])*wpe;
        for (const size_t w : nzw[p]) {
          uint32_t x = mk[w];
          const size_t r0 = (w - wb0)*16;
          while (x != 0u) {   // bits low to high = t increasing (as the byte path)
            const int b = __builtin_ctz(x) >> 1;
            const uint32_t v = (x >> (2*b)) & 3u;
            x &= ~(3u << (2*b));
            const size_t rr = r0 + b;
            const int e = st0[p] + static_cast<int>(rr/per);
            const size_t r = rr % per;
            const int d = static_cast<int>(r/n1), sh = static_cast<int>(r % n1);
            lst[cur[((static_cast<size_t>(v - 1)*n1 + sh)*ng + dg[d])*np + p]++] =
                e*n + d;
          }
        }
      });
      return t;
    };
    // (the send and receive lists at once)
    std::thread tl([&]() {tr = buildw(mp, vgd_cspr, vgd_cspm, E.ro, rl_h);});
    ts = buildw(smw, vgd_csps, vgd_cspn, E.so, sl_h);
    tl.join();
  } else {
    // the packed mask expanded on the host (1 byte per value, as before)
    std::vector<uint8_t> mh(static_cast<size_t>(std::max(vgd_csnrcv, 1))*per, 0);
    {
      auto mp = Kokkos::create_mirror_view_and_copy(HostMemSpace(), msk_);
      bool both = false;
      for (size_t t = 0; t < mh.size(); ++t) {
        const uint8_t v = static_cast<uint8_t>((mp(t >> 4) >> (2*(t & 15))) & 3u);
        both = both || (v == 3);
        mh[t] = v;
      }
      if (both) {VgdCsFatal("read mask: a value marked both inward and outward");}
    }
    // (2) the receivers' marks -> the owners
    std::vector<uint8_t> sm(static_cast<size_t>(std::max(vgd_csnsnd, 1))*per, 0);
  #if MPI_PARALLEL_ENABLED
    {
      std::vector<MPI_Request> req;
      for (int p = 0; p < np; ++p) {
        const int rk = vgd_csprk[p];
        if (vgd_cspn[p] > 0) {
          req.emplace_back();
          MPI_Irecv(sm.data() + static_cast<size_t>(vgd_csps[p])*per,
                    static_cast<int>(vgd_cspn[p]*per), MPI_UNSIGNED_CHAR, rk, 7102,
                    MPI_COMM_WORLD, &req.back());
        }
        if (vgd_cspm[p] > 0) {
          req.emplace_back();
          MPI_Isend(mh.data() + static_cast<size_t>(vgd_cspr[p])*per,
                    static_cast<int>(vgd_cspm[p]*per), MPI_UNSIGNED_CHAR, rk, 7102,
                    MPI_COMM_WORLD, &req.back());
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    }
  #endif
    auto build = [&](const uint8_t *mk, const std::vector<int> &st0,
                     const std::vector<int> &cnt, std::vector<int64_t> &off,
                     std::vector<int> &lst) {
      std::vector<int64_t> nn(nsl*np, 0);
      for (int p = 0; p < np; ++p) {
        for (int e = st0[p]; e < st0[p] + cnt[p]; ++e) {
          const uint8_t *q = mk + static_cast<size_t>(e)*per;
          for (int d = 0; d < n; ++d) {
            for (int sh = 0; sh < n1; ++sh) {
              const uint8_t v = q[static_cast<size_t>(d)*n1 + sh];
              if (v) {nn[((static_cast<size_t>(v - 1)*n1 + sh)*ng + dg[d])*np + p] += 1;}
            }
          }
        }
      }
      off.assign(nsl*(np + 1), 0);
      int64_t t = 0;
      for (size_t sl = 0; sl < nsl; ++sl) {
        for (int p = 0; p < np; ++p) {
          off[sl*(np + 1) + p] = t;
          t += nn[sl*np + p];
        }
        off[sl*(np + 1) + np] = t;
      }
      lst.assign(std::max<int64_t>(t, 1), 0);
      std::vector<int64_t> cur(nsl*np);
      for (size_t sl = 0; sl < nsl; ++sl) {
        for (int p = 0; p < np; ++p) {cur[sl*np + p] = off[sl*(np + 1) + p];}
      }
      for (int p = 0; p < np; ++p) {
        for (int e = st0[p]; e < st0[p] + cnt[p]; ++e) {
          const uint8_t *q = mk + static_cast<size_t>(e)*per;
          for (int d = 0; d < n; ++d) {
            for (int sh = 0; sh < n1; ++sh) {
              const uint8_t v = q[static_cast<size_t>(d)*n1 + sh];
              if (v) {
                lst[cur[((static_cast<size_t>(v - 1)*n1 + sh)*ng + dg[d])*np + p]++] =
                    e*n + d;
              }
            }
          }
        }
      }
      return t;
    };
    ts = build(sm.data(), vgd_csps, vgd_cspn, E.so, sl_h);
    tr = build(mh.data(), vgd_cspr, vgd_cspm, E.ro, rl_h);
  }
  vgd_cstb_h += tsub.seconds();
  tsub.reset();
  // the largest (pass, shell, group) slice (the pipelined slot size)
  for (size_t q = 0; q < nsl; ++q) {
    vgd_csxmax = std::max(vgd_csxmax, E.so[q*(np + 1) + np] - E.so[q*(np + 1)]);
    vgd_csxmax = std::max(vgd_csxmax, E.ro[q*(np + 1) + np] - E.ro[q*(np + 1)]);
  }
  // the order of (entry, direction) within a list is the same on both sides; the
  // counts must agree per (pass, shell, partner): checked by the message sizes below
  Kokkos::realloc(E.si, sl_h.size());
  Kokkos::realloc(E.ri, rl_h.size());
  {
    auto a = Kokkos::create_mirror_view(E.si);
    for (size_t q = 0; q < sl_h.size(); ++q) {a(q) = sl_h[q];}
    Kokkos::deep_copy(E.si, a);
    auto b = Kokkos::create_mirror_view(E.ri);
    for (size_t q = 0; q < rl_h.size(); ++q) {b(q) = rl_h[q];}
    Kokkos::deep_copy(E.ri, b);
  }
  vgd_cstb_u += tsub.seconds();
  E.alpha = vgd_alpha;
  E.ng = ng;
  E.use = ++vgd_csuse;
  vgd_cscur = c;
  g_vgd_csgen[c] += 1;
  vgd_csnbuild += 1.0;
  vgd_cstbuild += tb.seconds();
  if (global_variable::my_rank == 0 && vgd_csnbuild <= 1.0) {
    const double dn = static_cast<double>(vgd_csnsnd)*n*n1*2.0;
    std::cout << "<rad_m1> vet_gd (cubed sphere) compact lists: rank 0 sends " << ts
              << " values per sweep (dense " << dn << ": "
              << ((dn > 0.0) ? ts/dn : 0.0) << "), receives " << tr << "; built in "
              << tb.seconds() << " s (once per direction set)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsHaloC
//! \brief the compact exchange of shell i after pass `pass` of the sweep:
//! the listed (entry, direction) values of every partner (VetGdCsCompact)

void RadiationM1::VetGdCsHaloC(const int pass, const int i) {
#if MPI_PARALLEL_ENABLED
  const VgdCsCl &E = vgd_cscl[vgd_cscur];
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1;
  const int n = vgd_n;
  const int np = static_cast<int>(vgd_csprk.size());
  if (E.ng != 1) {VgdCsFatal("VetGdCsHaloC with direction-grouped lists");}
  const size_t sl = static_cast<size_t>(pass)*n1 + (i - indcs.is);
  const int64_t *so = &E.so[sl*(np + 1)];
  const int64_t *ro = &E.ro[sl*(np + 1)];
  const int64_t s0 = so[0], ns = so[np] - so[0];
  const int64_t r0 = ro[0], nr = ro[np] - ro[0];
  auto a_ = vgd_i;
  auto sb_ = vgd_cssb;
  auto rb_ = vgd_csrb;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto si_ = E.si;
  auto ri_ = E.ri;
  if (static_cast<int64_t>(sb_.extent(0)) < ns ||
      static_cast<int64_t>(rb_.extent(0)) < nr) {
    VgdCsFatal("VetGdCsHaloC: message buffer too small");
  }
  // with NCCL the device buffers (vgd_csdsb / vgd_csdrb) and no host
  // fence: pack, the NCCL group and the unpack all stream-ordered on the sweep's stream
  bool ncl = false;
#if M1_VGDCS_NCCL
  ncl = VgdNcclOn();
  if (ncl && static_cast<int64_t>(vgd_csdsb.extent(0)) <
             static_cast<int64_t>(sb_.extent(0))) {
    Kokkos::realloc(vgd_csdsb, sb_.extent(0));
    Kokkos::realloc(vgd_csdrb, rb_.extent(0));
  }
#endif
  Real *sbp = ncl ? vgd_csdsb.data() : sb_.data();
  Real *rbp = ncl ? vgd_csdrb.data() : rb_.data();
  if (ns > 0) {
    Kokkos::parallel_for("m1_vgdcs_cpack",
                         VgdLight(DevExeSpace(), 0, static_cast<int>(ns)),
    KOKKOS_LAMBDA(const int o) {
      const int x = si_(s0 + o);
      const int e = x/n, d = x - e*n;
      sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
    });
  }
#if M1_VGDCS_NCCL
  if (ncl) {
    static_assert(sizeof(Real) == sizeof(double), "csnccl-1031: ncclDouble");
    ncclComm_t nc = static_cast<ncclComm_t>(acc_ncomm);
    cudaStream_t st = DevExeSpace().cuda_stream();
    Kokkos::Timer tq;
    ncclGroupStart();
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_csprk[p];
      const int64_t cr = ro[p+1] - ro[p], cs = so[p+1] - so[p];
      if (cr > 0) {
        ncclRecv(rbp + (ro[p] - r0), static_cast<size_t>(cr), ncclDouble, rk, nc, st);
      }
      if (cs > 0) {
        ncclSend(sbp + (so[p] - s0), static_cast<size_t>(cs), ncclDouble, rk, nc, st);
        vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
      }
    }
    if (ncclGroupEnd() != ncclSuccess) {VgdCsFatal("csnccl-1031: ncclGroupEnd failed");}
    vgd_tmpi += tq.seconds();
  } else {
#endif
  Kokkos::fence();
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    const int64_t cr = ro[p+1] - ro[p], cs = so[p+1] - so[p];
    if (cr > 0) {
      req.emplace_back();
      MPI_Irecv(rb_.data() + (ro[p] - r0), static_cast<int>(cr), MPI_ATHENA_REAL, rk,
                7103, MPI_COMM_WORLD, &req.back());
    }
    if (cs > 0) {
      req.emplace_back();
      MPI_Isend(sb_.data() + (so[p] - s0), static_cast<int>(cs), MPI_ATHENA_REAL, rk,
                7103, MPI_COMM_WORLD, &req.back());
      vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  vgd_tmpi += tq.seconds();
#if M1_VGDCS_NCCL
  }   // !ncl
#endif
  vgd_nexch += 1.0;
  vgd_csvsent += static_cast<double>(ns);
  vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;
  if (nr > 0) {
    Kokkos::parallel_for("m1_vgdcs_cunpack",
                         VgdLight(DevExeSpace(), 0, static_cast<int>(nr)),
    KOKKOS_LAMBDA(const int o) {
      const int x = ri_(r0 + o);
      const int e = x/n, d = x - e*n;
      a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rbp[o];
    });
  }
#else
  (void) pass;
  (void) i;
#endif
}

// the per-partner offsets of the merged two-set exchange (VetGdCsHaloCX)
namespace {
constexpr int VGD_CSX_NS = 32;
struct VgdCsXIdx {int v[VGD_CSX_NS + 1];};
}   // namespace
//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::VetGdCsHaloCX
//! \brief the compact exchanges of shell i after pass `pass` of a fused
//! two-set sweep in ONE round: one pack kernel, one host fence, ONE message
//! per partner (piece p = [current set's values of p | next set's values of p], tag
//! 7107), one unpack kernel.  The values are those of the two VetGdCsHaloC calls (the
//! same copies into the same band cells: bitwise).  false (nothing done) with NCCL,
//! direction groups or env VGD_CSXM=0: the caller runs the two exchanges

bool RadiationM1::VetGdCsHaloCX(const int pass, const int i, const int eb) {
#if MPI_PARALLEL_ENABLED
  static int s_on = -1;
  if (s_on < 0) {
    const char *ev = std::getenv("VGD_CSXM");
    s_on = (ev != nullptr) ? std::atoi(ev) : 1;
  }
  if (s_on == 0 || eb < 0) {return false;}
#if M1_VGDCS_NCCL
  if (VgdNcclOn()) {return false;}
#endif
  const VgdCsCl &EA = vgd_cscl[vgd_cscur];
  const VgdCsCl &EB = vgd_cscl[eb];
  if (EA.ng != 1 || EB.ng != 1) {return false;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1;
  const int n = vgd_n;
  const int np = static_cast<int>(vgd_csprk.size());
  if (np > VGD_CSX_NS) {return false;}
  const size_t sl = static_cast<size_t>(pass)*n1 + (i - indcs.is);
  const int64_t *soA = &EA.so[sl*(np + 1)];
  const int64_t *roA = &EA.ro[sl*(np + 1)];
  const int64_t *soB = &EB.so[sl*(np + 1)];
  const int64_t *roB = &EB.ro[sl*(np + 1)];
  const int64_t s0A = soA[0], nsA = soA[np] - soA[0];
  const int64_t r0A = roA[0], nrA = roA[np] - roA[0];
  const int64_t s0B = soB[0], nsB = soB[np] - soB[0];
  const int64_t r0B = roB[0], nrB = roB[np] - roB[0];
  if (nsA + nsB >= static_cast<int64_t>(std::numeric_limits<int>::max()) ||
      nrA + nrB >= static_cast<int64_t>(std::numeric_limits<int>::max())) {
    return false;
  }
  // the buffers hold both sets (grown once when the two lists exceed one dense band)
  if (static_cast<int64_t>(vgd_cssb.extent(0)) < nsA + nsB) {
    Kokkos::realloc(vgd_cssb, static_cast<size_t>(nsA + nsB));
  }
  if (static_cast<int64_t>(vgd_csrb.extent(0)) < nrA + nrB) {
    Kokkos::realloc(vgd_csrb, static_cast<size_t>(nrA + nrB));
  }
  // per-set list offsets relative to the shell's first value, and the merged offsets
  VgdCsXIdx sa{}, sb{}, ra{}, rb{}, sm{}, rm{};
  for (int p = 0; p <= np; ++p) {
    sa.v[p] = static_cast<int>(soA[p] - s0A);
    sb.v[p] = static_cast<int>(soB[p] - s0B);
    ra.v[p] = static_cast<int>(roA[p] - r0A);
    rb.v[p] = static_cast<int>(roB[p] - r0B);
    sm.v[p] = sa.v[p] + sb.v[p];
    rm.v[p] = ra.v[p] + rb.v[p];
  }
  auto a_ = vgd_i;
  auto b_ = vgd_i_xb;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto siA_ = EA.si;
  auto riA_ = EA.ri;
  auto siB_ = EB.si;
  auto riB_ = EB.ri;
  Real *sbp = vgd_cssb.data();
  Real *rbp = vgd_csrb.data();
  const int nsa = static_cast<int>(nsA), nsb = static_cast<int>(nsB);
  const int nra = static_cast<int>(nrA), nrb = static_cast<int>(nrB);
  const int npk = np;
  if (nsa + nsb > 0) {
    Kokkos::parallel_for("m1_vgdcs_cpackx", VgdLight(DevExeSpace(), 0, nsa + nsb),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nsa);
      const int o = b ? (c - nsa) : c;
      const VgdCsXIdx &sx = b ? sb : sa;
      int p = 0;
      while (p + 1 < npk && sx.v[p+1] <= o) {++p;}
      const int x = b ? siB_(s0B + o) : siA_(s0A + o);
      const int e = x/n, d = x - e*n;
      const int pos = sm.v[p] + (b ? (sa.v[p+1] - sa.v[p]) : 0) + (o - sx.v[p]);
      sbp[pos] = b ? b_(S_(e,0),d,S_(e,1),S_(e,2),i) : a_(S_(e,0),d,S_(e,1),S_(e,2),i);
    });
  }
  Kokkos::fence();
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  req.reserve(2*np);
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    const int cr = rm.v[p+1] - rm.v[p], cs = sm.v[p+1] - sm.v[p];
    if (cr > 0) {
      req.emplace_back();
      MPI_Irecv(rbp + rm.v[p], cr, MPI_ATHENA_REAL, rk, 7107, MPI_COMM_WORLD,
                &req.back());
    }
    if (cs > 0) {
      req.emplace_back();
      MPI_Isend(sbp + sm.v[p], cs, MPI_ATHENA_REAL, rk, 7107, MPI_COMM_WORLD,
                &req.back());
      vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  vgd_tmpi += tq.seconds();
  vgd_nexch += 2.0;
  vgd_csvsent += static_cast<double>(nsA + nsB);
  vgd_csvdense += static_cast<double>(vgd_csnsnd)*n*2;
  if (nra + nrb > 0) {
    Kokkos::parallel_for("m1_vgdcs_cunpackx", VgdLight(DevExeSpace(), 0, nra + nrb),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nra);
      const int o = b ? (c - nra) : c;
      const VgdCsXIdx &rx = b ? rb : ra;
      int p = 0;
      while (p + 1 < npk && rx.v[p+1] <= o) {++p;}
      const int x = b ? riB_(r0B + o) : riA_(r0A + o);
      const int e = x/n, d = x - e*n;
      const Real v = rbp[rm.v[p] + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (o - rx.v[p])];
      if (b) {
        b_(R_(e,0),d,R_(e,1),R_(e,2),i) = v;
      } else {
        a_(R_(e,0),d,R_(e,1),R_(e,2),i) = v;
      }
    });
  }
  return true;
#else
  (void) pass;
  (void) i;
  (void) eb;
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsHalo
//! \brief fill every band position of a(m, v, k, j, i), v in [0, nv), shells [i0, i1],
//! with its owner's interior value (exact copies; the shell range chunked to the buffer).
//! ws >= 0 with the depth-sorted lists of vet_gd_band_global (the full
//! sphere, VetGdGlobalBandInit): only the positions within depth ws (a prefix of the
//! on-rank list and of every partner's segment); ws < 0 (the cubed sphere): all

template <class V>
void RadiationM1::VetGdCsHalo(V &a, const int nv, const int i0, const int i1,
                              const bool loc, const int ws) {
  auto a_ = a;
  const bool dp = (ws >= 0) && (vgd_gband || vgd_raw) && !vgd_gdnl.empty();
  const int wq = dp ? std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1) : 0;
  const int nl = loc ? (dp ? vgd_gdnl[wq] : vgd_csnloc) : 0;
  {
    const int ni = i1 - i0 + 1;
    if (nl > 0 && ni > 0) {
      auto L_ = vgd_csloc;
      par_for("m1_vgdcs_loc", DevExeSpace(), 0, nl - 1, 0, nv - 1, 0, ni - 1,
      KOKKOS_LAMBDA(const int e, const int v, const int ii) {
        a_(L_(e,0),v,L_(e,1),L_(e,2),i0+ii) = a_(L_(e,3),v,L_(e,4),L_(e,5),i0+ii);
      });
    }
  }
#if MPI_PARALLEL_ENABLED
  if (vgd_hmpi) {
    const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
    // the full-sphere global band packs into DEVICE buffers (vgd_sbuf /
    // vgd_rbuf, CUDA-aware MPI as the sp slot halo); the cubed sphere keeps pinned host
    const int cap = static_cast<int>((dp ? vgd_sbuf.extent(0) : vgd_cssb.extent(0))/nent);
    const int nich = std::max(1, cap/nv);
    for (int c0 = i0; c0 <= i1; c0 += nich) {
      const int c1 = std::min(c0 + nich - 1, i1);
      const int ni = c1 - c0 + 1;
      const int nvi = nv*ni;
      auto sb_ = vgd_cssb;
      auto rb_ = vgd_csrb;
      auto S_ = vgd_cssl;
      auto R_ = vgd_csrl;
      auto DS_ = vgd_gdps;
      auto DR_ = vgd_gdpr;
      auto gsb_ = vgd_sbuf;
      auto grb_ = vgd_rbuf;
      Real *sbp = dp ? gsb_.data() : sb_.data();
      Real *rbp = dp ? grb_.data() : rb_.data();
      if (vgd_csnsnd > 0 && dp) {
        // the depth-ws prefix of every partner segment (same offsets)
        par_for("m1_vgdgb_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DS_(e) > wq) {return;}
          gsb_(static_cast<size_t>(e)*nvi + v*ni + ii) =
              a_(S_(e,0),v,S_(e,1),S_(e,2),c0+ii);
        });
      } else if (vgd_csnsnd > 0) {
        par_for("m1_vgdcs_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          sb_(static_cast<size_t>(e)*nvi + v*ni + ii) =
              a_(S_(e,0),v,S_(e,1),S_(e,2),c0+ii);
        });
      }
      Kokkos::fence();
      Kokkos::Timer tq;
      std::vector<MPI_Request> req;
      for (size_t p = 0; p < vgd_csprk.size(); ++p) {
        const int rk = vgd_csprk[p];
        const int cm = dp ? vgd_gdrn[wq][p] : vgd_cspm[p];
        const int cn = dp ? vgd_gdsn[wq][p] : vgd_cspn[p];
        if (cm > 0) {
          req.emplace_back();
          MPI_Irecv(rbp + static_cast<size_t>(vgd_cspr[p])*nvi, cm*nvi,
                    MPI_ATHENA_REAL, rk, 7101, MPI_COMM_WORLD, &req.back());
        }
        if (cn > 0) {
          req.emplace_back();
          MPI_Isend(sbp + static_cast<size_t>(vgd_csps[p])*nvi, cn*nvi,
                    MPI_ATHENA_REAL, rk, 7101, MPI_COMM_WORLD, &req.back());
          vgd_nbyte += static_cast<Real>(cn)*nvi*sizeof(Real);
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
      vgd_tmpi += tq.seconds();
      vgd_nexch += 1.0;
      if (!loc) {
        vgd_csvsent += static_cast<double>(vgd_csnsnd)*nvi;
        vgd_csvdense += static_cast<double>(vgd_csnsnd)*nvi;
      }
      if (vgd_csnrcv > 0 && dp) {
        par_for("m1_vgdgb_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DR_(e) > wq) {return;}
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              grb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      } else if (vgd_csnrcv > 0) {
        par_for("m1_vgdcs_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              rb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      }
    }
  }
#endif
}

template void RadiationM1::VetGdCsHalo<DvceArray5D<Real>>(DvceArray5D<Real> &a,
                                                         const int nv, const int i0,
                                                         const int i1, const bool loc,
                                                         const int ws);
template void RadiationM1::VetGdCsHalo<RadiationM1::VgdRag>(RadiationM1::VgdRag &a,
                                                            const int nv, const int i0,
                                                            const int i1, const bool loc,
                                                            const int ws);

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsSource
//! \brief the band of ln chi, ln S (vgd_cs) for every shell of the sweep

void RadiationM1::VetGdCsSource() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ilo = indcs.is + (clean ? 0 : vlat_icut), ie = indcs.ie;
  VetGdCsHalo(vgd_cs, 2, ilo, ie, true);
  Kokkos::fence();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsTab
//! \brief the host tables of a direction set d on cs (VgdTabCompute): no walls, the
//! identity map, the branch n . r_hat of every band cell

void RadiationM1::VetGdCsTab(const std::vector<double> &d) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c2 = indcs.nx2 + 2*vgd_w, c3 = indcs.nx3 + 2*vgd_w;
  const int n = vgd_n;
  VgdTabHost &H = vgd_tabh;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        H.w(m,k,j) = 0;
        const double *rh = &vgd_csr_h[((static_cast<size_t>(m)*c3 + k)*c2 + j)*3];
        for (int q = 0; q < n; ++q) {
          H.mp(m,k,j,q) = q;
          H.mr(m,k,j,q) = d[4*q]*rh[0] + d[4*q+1]*rh[1] + d[4*q+2]*rh[2];
        }
      }
    }
  }
  H.nwall = 0;
  H.nexact = 0;
  H.nmap = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsSweep
//! \brief one sweep (inward then outward) of every direction on the cubed sphere (the
//! VetGdSweep body with the cs upwind stencil; every band value from VetGdCsHalo)

void RadiationM1::VetGdCsSweep() {
  if (vgd_xff_prep) {return;}   // nothing to prepare (lists made below)
  // vet_gd_rotate_xfade_fuse: the next set swept in the same shell kernels
  // (its intensities vgd_i_xb, tables vgd_tabd, compact lists and pass lists its own)
  const bool xffs = vgd_xff_now;
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int c2b = indcs.nx2 + 2*wb, c3b = indcs.nx3 + 2*wb;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vgd_scut;
  auto cs_ = vgd_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vlat_cnt;
  auto mr_ = vgd_mr;
  auto cr_ = vgd_csr;
  auto bl_ = vgd_csblk;
  auto cm_ = vgd_cscm;
  auto nb_ = vgd_csnb;
  auto ow_ = vgd_csown;
  const int nc = vgd_csn;
  const Real h = 0.5*M_PI/nc;
  const bool bfc = vet_bc_face;
  const bool izf = vgd_izf && (ibc_x1min == M1_IBC_REFLECT) && lcut == 0;   // hrsplit
  const bool spl = vet_seg_split;
  const int spe = VetSegMode();
  const int vim = vet_interp;
  const bool emi = VetBcEmit() && (lcut == 0) &&
                   (ibc_x1min != M1_IBC_REFLECT);
  const Real iemit = VetBathLo();
  const bool tbh = vgd_tbath && (ibc_x1max == M1_IBC_MARSHAK) && (iebath_x1max != 0.0);
  const Real itop = tbh ? iebath_x1max : 0.0;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const Real lcr = log(kVetCornerRatio);
  vgd_hsweep += xffs ? 2 : 1;
  // the compact lists of this set (built on its first sweep); every remote
  // read is checked against the read mask (a miss would be a value the lists did not
  // send: FATAL below, never silent)
  const bool ccl = vgd_cscomp;
  int xcurb = -1;   // the next set's compact-list entry
  if (ccl) {
    VetGdCsCompact();
    if (xffs) {
      // (the current set's entry first, then the next set's: the order of the two sweeps,
      // so the LRU entry is the same)
      const int ca = vgd_cscur;
      std::swap(vgd_dir, vgd_tabd.dir);
      std::swap(vgd_mr, vgd_tabd.mr);
      std::swap(vgd_alpha, vgd_xf_alt);
      VetGdCsCompact();
      xcurb = vgd_cscur;
      std::swap(vgd_alpha, vgd_xf_alt);
      std::swap(vgd_mr, vgd_tabd.mr);
      std::swap(vgd_dir, vgd_tabd.dir);
      vgd_cscur = ca;
    }
    Kokkos::deep_copy(vgd_csmiss, 0);
  }
  auto msk_ = vgd_cscl[vgd_cscur].msk;
  auto rx_ = vgd_csrix;
  auto miss_ = vgd_csmiss;
  const int64_t per_ = vgd_csper;   // the mask stride per entry
  // vet_gd_band_ring on the cubed sphere (VgdRagAlloc / VetGdInit: with the
  // compact lists and an MPI band): the sp ring of rad_m1_vetgd.cpp (VetGdRingFlush, the
  // same slab operations in the same places of the shell loop), its saved inward band
  // values = the receive list of each shell's inward compact exchange.  The record of
  // shell i (vgd_hxe[i], pass 0) holds the packed band address of every value received
  // by that exchange; the addresses of all shells live in one array kept in the unused
  // entry vgd_hxe[n1c] (pass 1 of the ghost shell 0).  Pure copies: the sweep reads the
  // values of the dense-band layout (bitwise)
  const bool rga = vgd_ring && ccl;
  if (vgd_ring && !ccl) {VgdCsFatal("vet_gd_band_ring without the compact lists");}
  if (xffs && ccl && vgd_csxg > 1) {
    VgdCsFatal("vet_gd_rotate_xfade_fuse: not with direction groups (VetGdXffOk)");
  }
  // (the next set's records in vgd_hxe_alt, unused on the cubed sphere
  // otherwise, from its own compact lists; the ring's second array is its intensities)
  auto ringrec = [&](std::vector<VgdHlEntry> &hxe, const VgdCsCl &E, const Real alr) {
    const int n1c = indcs.nx1 + 2*indcs.ng;
    if (static_cast<int>(hxe.size()) != 2*n1c) {hxe.resize(2*n1c);}
    const int np = static_cast<int>(vgd_csprk.size());
    const int64_t a0 = E.ro[0];
    const int G = E.ng;   // the shell's groups are contiguous in the list
    const int64_t a1 = E.ro[(static_cast<size_t>(n1 - 1)*G + G - 1)*(np + 1) + np];
    const int64_t na = a1 - a0;
    if (static_cast<int64_t>(pmy_pack->nmb_thispack)*n*c3b*c2b >= INT32_MAX ||
        na >= INT32_MAX) {
      VgdCsFatal("vet_gd_band_ring: band addresses >= 2^31 (32-bit records)");
    }
    VgdHlEntry &B = hxe[n1c];
    B.gen = -1;
    B.exact = false;
    if (B.dl.extent_int(0) < na + 1) {
      B.dl = DvceArray1D<int>();
      Kokkos::realloc(B.dl, na + na/8 + 1024);
    }
    auto bd_ = B.dl;
    auto R_ = vgd_csrl;
    auto ri_ = E.ri;
    if (na > 0) {
      par_for("m1_vgdcs_ringrec", DevExeSpace(), 0, static_cast<int>(na) - 1,
      KOKKOS_LAMBDA(const int o) {
        const int x = ri_(a0 + o);
        const int e = x/n, d = x - e*n;
        bd_(o) = ((R_(e,0)*n + d)*c3b + R_(e,1))*c2b + R_(e,2);
      });
    }
    for (int l = 0; l < n1; ++l) {
      VgdHlEntry &X = hxe[is + l];
      const int64_t r0 = E.ro[static_cast<size_t>(l)*G*(np + 1)];
      const int64_t r1 = E.ro[(static_cast<size_t>(l)*G + G - 1)*(np + 1) + np];
      X.dl = DvceArray1D<int>(bd_.data() + (r0 - a0), r1 - r0);
      X.dw = DvceArray1D<int>();
      X.nr = static_cast<int>(r1 - r0);
      X.nw = 0;
      X.exact = true;
      X.gen = 1;
      X.alpha = alr;
      X.scut = vgd_scut;
    }
  };
  if (rga) {
    ringrec(vgd_hxe, vgd_cscl[vgd_cscur], vgd_alpha);
    if (xffs) {ringrec(vgd_hxe_alt, vgd_cscl[xcurb], vgd_xf_alt);}
  }
  if (rga) {vgd_cur = DevExeSpace();}   // the ring's copy kernel: this sweep's stream
  // vet_gd_cs_xgroups = G > 1 (compact lists, MPI band): per (shell, group)
  // one kernel and one exchange, pipelined; the shell kernel itself is the same code
  // over a direction range (every (cell, d) value from the same inputs: bitwise G = 1)
#if MPI_PARALLEL_ENABLED
  const bool pipe = ccl && (vgd_csxg > 1);
#else
  const bool pipe = false;
#endif
  const int G = pipe ? vgd_csxg : 1;
#if MPI_PARALLEL_ENABLED
  struct Xi {
    int pass, l, i, g, slot;
    bool sent;
    std::vector<MPI_Request> rq;
  };
  std::vector<Xi> xfl;          // in flight, oldest first
  size_t xh = 0;                // first not finished
  int64_t xcnt = 0;
  const int xS = 2*G + 2;       // slots: a slot is reused only after its unpack ran
  std::vector<VgdCsEv> xevp(pipe ? xS : 0), xevu(pipe ? xS : 0);
  const int64_t xmx = std::max<int64_t>(vgd_csxmax, 1);
  if (pipe && static_cast<int64_t>(vgd_csxsb.extent(0)) < xS*xmx) {
    Kokkos::fence();
    vgd_csxsb = decltype(vgd_csxsb)();
    vgd_csxrb = decltype(vgd_csxrb)();
    Kokkos::realloc(vgd_csxsb, xS*xmx);
    Kokkos::realloc(vgd_csxrb, xS*xmx);
  }
  const int xnp = static_cast<int>(vgd_csprk.size());
  auto xsl = [&](const Xi &x) {
    return (static_cast<size_t>(x.pass)*n1 + x.l)*G + x.g;
  };
  // post the sends of every item whose pack has run (in item order)
  auto xprog = [&]() {
    const VgdCsCl &E = vgd_cscl[vgd_cscur];
    for (size_t t = xh; t < xfl.size(); ++t) {
      Xi &x = xfl[t];
      if (x.sent) {continue;}
      if (!xevp[x.slot].Done()) {break;}
      const int64_t *so = &E.so[xsl(x)*(xnp + 1)];
      const int64_t s0 = so[0];
      Real *b = vgd_csxsb.data() + x.slot*xmx;
      for (int p = 0; p < xnp; ++p) {
        const int64_t cs = so[p+1] - so[p];
        if (cs > 0) {
          x.rq.emplace_back();
          MPI_Isend(b + (so[p] - s0), static_cast<int>(cs), MPI_ATHENA_REAL,
                    vgd_csprk[p], 7104, MPI_COMM_WORLD, &x.rq.back());
          vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
        }
      }
      vgd_csvsent += static_cast<double>(so[xnp] - s0);
      x.sent = true;
    }
  };
  // launch the pack of (pass, l, g) and post its receives
  auto xstart = [&](const int pass, const int l, const int g) {
    const VgdCsCl &E = vgd_cscl[vgd_cscur];
    Xi x;
    x.pass = pass;
    x.l = l;
    x.i = is + l;
    x.g = g;
    x.slot = static_cast<int>(xcnt++ % xS);
    x.sent = false;
    const size_t sl = xsl(x);
    const int64_t *so = &E.so[sl*(xnp + 1)];
    const int64_t *ro = &E.ro[sl*(xnp + 1)];
    const int64_t s0 = so[0], ns = so[xnp] - so[0];
    const int64_t r0 = ro[0];
    xevu[x.slot].Wait();   // the slot's previous unpack has read its receive buffer
    Real *rb = vgd_csxrb.data() + x.slot*xmx;
    for (int p = 0; p < xnp; ++p) {
      const int64_t cr = ro[p+1] - ro[p];
      if (cr > 0) {
        x.rq.emplace_back();
        MPI_Irecv(rb + (ro[p] - r0), static_cast<int>(cr), MPI_ATHENA_REAL,
                  vgd_csprk[p], 7104, MPI_COMM_WORLD, &x.rq.back());
      }
    }
    if (ns > 0) {
      auto a_ = vgd_i;
      auto sb_ = vgd_csxsb;
      auto S_ = vgd_cssl;
      auto si_ = E.si;
      const int64_t bo = x.slot*xmx;
      const int i = x.i;
      Kokkos::parallel_for("m1_vgdcs_xpack", VgdLight(DevExeSpace(), 0,
                                                     static_cast<int>(ns)),
      KOKKOS_LAMBDA(const int o) {
        const int xx = si_(s0 + o);
        const int e = xx/n, d = xx - e*n;
        sb_(bo + o) = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
      });
    }
    xevp[x.slot].Record();
    if (g == 0) {vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;}
    xfl.push_back(std::move(x));
    xprog();
  };
  // complete the oldest item: its messages, then its unpack (stream order)
  auto xfinish = [&]() {
    const VgdCsCl &E = vgd_cscl[vgd_cscur];
    Xi &x = xfl[xh];
    Kokkos::Timer tq;
    while (true) {
      xprog();
      if (x.sent) {
        int flag = 0;
        MPI_Testall(static_cast<int>(x.rq.size()), x.rq.data(), &flag,
                    MPI_STATUSES_IGNORE);
        if (flag) {break;}
      }
    }
    vgd_tmpi += tq.seconds();
    vgd_nexch += 1.0;
    const size_t sl = xsl(x);
    const int64_t *ro = &E.ro[sl*(xnp + 1)];
    const int64_t r0 = ro[0], nr = ro[xnp] - ro[0];
    if (nr > 0) {
      auto a_ = vgd_i;
      auto rb_ = vgd_csxrb;
      auto R_ = vgd_csrl;
      auto ri_ = E.ri;
      const int64_t bo = x.slot*xmx;
      const int i = x.i;
      Kokkos::parallel_for("m1_vgdcs_xunpack", VgdLight(DevExeSpace(), 0,
                                                     static_cast<int>(nr)),
      KOKKOS_LAMBDA(const int o) {
        const int xx = ri_(r0 + o);
        const int e = xx/n, d = xx - e*n;
        a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rb_(bo + o);
      });
    }
    xevu[x.slot].Record();
    x.rq.clear();
    ++xh;
    if (xh == xfl.size()) {
      xfl.clear();
      xh = 0;
    }
  };
#endif
  // per pass the codes ((m k) j) d of the (cell, direction) pairs of that
  // pass's branch (inward: n.r_hat < 0), in the flat range's order; the shell kernel runs
  // over the list (half the range: the other branch's threads only returned, and a warp
  // of mixed branches cost the whole warp).  Made each sweep from vgd_mr (one scan per
  // pass); the same threads compute the same values.  G > 1 keeps the full range
  const bool usel = (G == 1);
  const int nkl = ke - ks + 1, njl = je - js + 1;
  const int nlt = (nmb1 + 1)*nkl*njl*n;
  int nlist[2] = {0, 0};
  if (usel) {
    if (vgd_cslst.extent_int(0) < 2*nlt) {
      vgd_cslst = DvceArray1D<int>();
      Kokkos::realloc(vgd_cslst, 2*nlt);
    }
    auto L_ = vgd_cslst;
    for (int p = 0; p < 2; ++p) {
      const bool inp = (p == 0);
      const int off = p*nlt;
      int cnt = 0;
      Kokkos::parallel_scan("m1_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0, nlt),
      KOKKOS_LAMBDA(const int c, int &u, const bool fin) {
        int t = c;
        const int d = t % n;
        t /= n;
        const int j = js + t % njl;
        t /= njl;
        const int k = ks + t % nkl;
        const int m = t/nkl;
        if (inp != (mr_(m,k,j,d) >= 0.0)) {
          if (fin) {L_(off + u) = c;}
          u += 1;
        }
      }, cnt);
      nlist[p] = cnt;
    }
  }
  auto lst_ = vgd_cslst;
  // the next set's pass lists (the same scan over its own n . r_hat)
  int nlistb[2] = {0, 0};
  if (xffs) {
    if (!usel) {VgdCsFatal("vet_gd_rotate_xfade_fuse needs the pass lists (G = 1)");}
    if (vgd_cslst_xb.extent_int(0) < 2*nlt) {
      vgd_cslst_xb = DvceArray1D<int>();
      Kokkos::realloc(vgd_cslst_xb, 2*nlt);
    }
    auto L_ = vgd_cslst_xb;
    auto mrb_ = vgd_tabd.mr;
    for (int p = 0; p < 2; ++p) {
      const bool inp = (p == 0);
      const int off = p*nlt;
      int cnt = 0;
      Kokkos::parallel_scan("m1_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0, nlt),
      KOKKOS_LAMBDA(const int c, int &u, const bool fin) {
        int t = c;
        const int d = t % n;
        t /= n;
        const int j = js + t % njl;
        t /= njl;
        const int k = ks + t % nkl;
        const int m = t/nkl;
        if (inp != (mrb_(m,k,j,d) >= 0.0)) {
          if (fin) {L_(off + u) = c;}
          u += 1;
        }
      }, cnt);
      nlistb[p] = cnt;
    }
  }
  auto lstb_ = vgd_cslst_xb;
  // cs early send (see the shell loop)
  static int s_early = -1;
  if (s_early < 0) {
    const char *ee = std::getenv("VGD_CS_EARLY");
    s_early = (ee != nullptr) ? std::atoi(ee) : 1;
  }
  bool early = (s_early != 0) && ccl && !pipe && (G == 1);
#if !MPI_PARALLEL_ENABLED
  early = false;
#endif
#if M1_VGDCS_NCCL
  if (VgdNcclOn()) {early = false;}
#endif
  if (early && vgd_cscl[vgd_cscur].ng != 1) {early = false;}
  if (xffs) {early = false;}   // (the early send: one set per shell kernel)
#if defined(KOKKOS_ENABLE_HIP) && MPI_PARALLEL_ENABLED
  // rsthip-1121: MI300A, ROCm 7.2: more than 2 ranks on one GPU with the default
  // GPU_MAX_HW_QUEUES (4: 5/5, 3: 1/1) hang in the early send's event wait (hgs 6 ranks,
  // 2 MI300A; GPU_MAX_HW_QUEUES=2 or VGD_CS_EARLY=0 run, with the same bits)
  const char *hq = std::getenv("GPU_MAX_HW_QUEUES");
  static bool s_hqchk = false;
  if (early && !s_hqchk && (hq == nullptr || std::atoi(hq) > 2)) {
    s_hqchk = true;
    char bus[32] = {0};
    (void) hipDeviceGetPCIBusId(bus, sizeof(bus), Kokkos::device_id());
    MPI_Comm nc;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nc);
    int nn = 1, sh = 0;
    MPI_Comm_size(nc, &nn);
    std::vector<char> all(32*static_cast<size_t>(nn));
    MPI_Allgather(bus, 32, MPI_CHAR, all.data(), 32, MPI_CHAR, nc);
    MPI_Comm_free(&nc);
    for (int q = 0; q < nn; ++q) {sh += (std::strncmp(bus, &all[32*q], 32) == 0);}
    MPI_Allreduce(MPI_IN_PLACE, &sh, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (sh > 2) {
      VgdCsFatal(std::to_string(sh) + " ranks share one GPU: the early send hangs "
                 "there (HIP hardware queues); set GPU_MAX_HW_QUEUES=2 (same bits) or "
                 "VGD_CS_EARLY=0, or use at most 2 ranks per GPU");
    }
  }
#endif
  // the early send's cells are not computed twice: m1_vgdcs_estamp stamps
  // each (m, k, j, d) m1_vgdcs_eshell computes (vgd_est[0], stamp vgd_esid of the shell;
  // eshell counts the clamps itself) and the shell kernel skips the stamped ones: the
  // same values from the same code, written once (bitwise)
  const int nkb = ke - ks + 1, njb = je - js + 1;
  if (early && (vgd_est[0].extent_int(0) < nmb1 + 1 || vgd_est[0].extent_int(1) < nkb ||
                vgd_est[0].extent_int(2) < njb || vgd_est[0].extent_int(3) < n)) {
    vgd_est[0] = DvceArray4D<int>("vgd_est", nmb1 + 1, nkb, njb, n);   // (zero)
    vgd_esid = 0;
  }
  auto est_ = vgd_est[0];
  int esid = 0;   // (0: no early cells; never a stamp)
  VgdCsEv eev, efork;
  static int s_lb = -1;   // cslb (see the shell kernel)
  if (s_lb < 0) {
    const char *el = std::getenv("VGD_CS_LB");
    // 6 (VGD_CS_LB=0: old launch); 2 on HIP, where the bound is min waves per EU and 6
    // spilled the body (HG small cs deck 13083259: vet_gd sweeps 5.08 -> 2.73 s, outputs
    // bitwise); 8 on CUDA (H200 RSG -1.1 %, HG -0.8 %, bitwise)
#if defined(KOKKOS_ENABLE_HIP)
    s_lb = (el != nullptr) ? std::atoi(el) : 2;
#else
    s_lb = (el != nullptr) ? std::atoi(el) : 8;
#endif
  }
  VetGdRingSweepBegin(xffs);   // (the other set's array second)
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    const int nsh = n1 - lcut;
    for (int q = 0; q < nsh; ++q) {
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      if (rga) {
        // this shell's wide slab: inward from scratch, outward with its interior and
        // the band values of its inward exchange
        VetGdRingFlush();
        if (vgd_rg_state[i] < 0) {
          VgdRingOp o;
          o.kind = 3;
          o.i = i;
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
        VetGdRingFlush();
      }
      for (int g = 0; g < G; ++g) {
      const int d0 = (g*n + G - 1)/G, d1 = ((g + 1)*n + G - 1)/G;   // d*G/n == g
#if MPI_PARALLEL_ENABLED
      if (pipe && q > 0) {   // the band of the previous shell, this group
        Kokkos::Timer tf;
        xfinish();
        vgd_thalo += tf.seconds();
      }
#endif
      // flat (m, k, j, d) range, d fastest (as the MDRange before);
      // cstwin-geo: with the fused twin one thread computes the geometry once and the
      // values of both fields
      const int nkk = ke - ks + 1, njj = je - js + 1, ndd = d1 - d0;
      const int nper = usel ? nlist[pass] : (nmb1 + 1)*nkk*njj*ndd;
      const int lo = pass*nlt;
      // (cs early send): the body of the shell kernel for one (cell, direction)
      // as a device lambda, run by the shell kernel over the pass's list and (early send)
      // first over the shell's compact SEND list alone; cnt = false there: the clamped /
      // miss counters are counted once, by the shell kernel (which recomputes the same
      // values from the same inputs: bitwise)
#define M1_VGDCS_SHELL_NAME shell_
#include "rad_m1_vetgd_cs_shell.inc"
#undef M1_VGDCS_SHELL_NAME
      // cs early send (env VGD_CS_EARLY, default on; compact lists, G = 1, no NCCL): the
      // values this shell SENDS first, then their pack (pinned host buffer) and an event;
      // the full shell kernel runs while the host exchanges (finished below at the halo
      // point: event wait instead of a fence, MPI, unpack)
      if (early) {
        const VgdCsCl &E = vgd_cscl[vgd_cscur];
        const int np = static_cast<int>(vgd_csprk.size());
        const size_t sl = static_cast<size_t>(pass)*n1 + l;
        const int64_t *so = &E.so[sl*(np + 1)];
        const int64_t es0 = so[0], ens = so[np] - so[0];
        if (ens > 0) {
          auto si_ = E.si;
          auto S_ = vgd_cssl;
          vgd_esid = (vgd_esid >= std::numeric_limits<int>::max() - 1) ? 1
                                                                       : (vgd_esid + 1);
          if (vgd_esid == 1) {Kokkos::deep_copy(est_, 0);}   // (first use or wrapped)
          esid = vgd_esid;
          const int sid_ = esid;
          Kokkos::parallel_for("m1_vgdcs_estamp", VgdLight(DevExeSpace(), 0,
                                                            static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            est_(S_(e,0), S_(e,1) - ks, S_(e,2) - js, d) = sid_;
          });
          // the send cells and their pack on the second stream, concurrent
          // with the shell kernel (which skips them); the sweep stream waits for them
          // after the shell kernel (efork / eev)
          VgdEs2Make();
          efork.Record();
          efork.StreamWait(vgd_es2);
          Kokkos::parallel_for("m1_vgdcs_eshell",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            shell_(S_(e,0), S_(e,1), S_(e,2), d, true);
          });
          auto a_ = vgd_i;
          Real *sbp = vgd_cssb.data();
          Kokkos::parallel_for("m1_vgdcs_cpack",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
          });
        }
        eev.Record((esid != 0) ? vgd_es2 : DevExeSpace());
      }
      // the fused two-set shell kernel: this pass's list of the current set,
      // then the next set's, one launch; the next set's body is the include above with
      // that set's names (direction table, n . r_hat, intensities, read mask)
#define M1_VGDCS_LAUNCHN(F, N) \
      if (s_lb == 8) {                                                                 \
        Kokkos::parallel_for("m1_vgdcs_shell_xff", Kokkos::Experimental::require(      \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 8>>(DevExeSpace(), 0, N),    \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 6) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell_xff", Kokkos::Experimental::require(      \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 6>>(DevExeSpace(), 0, N),    \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 2) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell_xff", Kokkos::Experimental::require(      \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 2>>(DevExeSpace(), 0, N),    \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 1) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell_xff", Kokkos::Experimental::require(      \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 1>>(DevExeSpace(), 0, N),    \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else {                                                                         \
        Kokkos::parallel_for("m1_vgdcs_shell_xff", VgdLight(DevExeSpace(), 0, N), F);  \
      }
      const int esid_ = esid;
      esid = 0;
      if (xffs) {
        auto xsb_ = [&]() {
          auto dir_ = vgd_tabd.dir;
          auto mr_ = vgd_tabd.mr;
          auto vi_ = vgd_i_xb;
          auto msk_ = vgd_cscl[(xcurb >= 0) ? xcurb : vgd_cscur].msk;
#define M1_VGDCS_SHELL_NAME shellx_
#include "rad_m1_vetgd_cs_shell.inc"
#undef M1_VGDCS_SHELL_NAME
          return shellx_;
        }();
        const int nA = nlist[pass], nB = nlistb[pass];
        auto mainx_ = KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nA);
          int t = sb ? lstb_(lo + c - nA) : lst_(lo + c);
          const int d = d0 + t % ndd;
          t /= ndd;
          const int j = js + t % njj;
          t /= njj;
          const int k = ks + t % nkk;
          const int m = t/nkk;
          if (sb) {
            xsb_(m, k, j, d, true);
          } else {
            shell_(m, k, j, d, true);
          }
        };
        if (nA + nB > 0) {M1_VGDCS_LAUNCHN(mainx_, nA + nB)}
      } else {
#if defined(KOKKOS_ENABLE_CUDA)
      auto mainf_ = KOKKOS_LAMBDA(const int c) {
        int t = usel ? lst_(lo + c) : c;   // (d0 = 0, ndd = n with the list)
        const int d = d0 + t % ndd;
        t /= ndd;
        const int j = js + t % njj;
        t /= njj;
        const int k = ks + t % nkk;
        const int m = t/nkk;
        if (esid_ != 0 && est_(m, k - ks, j - js, d) == esid_) {return;}
        shell_(m, k, j, d, true);
      };
      // (the launch variants of VGD_CS_LB, below)
#define M1_VGDCS_LAUNCH(F) \
      if (s_lb == 8) {                                                                 \
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(          \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 8>>(DevExeSpace(), 0, nper), \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 6) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(          \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 6>>(DevExeSpace(), 0, nper), \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 2) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(          \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 2>>(DevExeSpace(), 0, nper), \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else if (s_lb == 1) {                                                          \
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(          \
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 1>>(DevExeSpace(), 0, nper), \
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);              \
      } else {                                                                         \
        Kokkos::parallel_for("m1_vgdcs_shell", VgdLight(DevExeSpace(), 0, nper), F);   \
      }
      M1_VGDCS_LAUNCH(mainf_)
#undef M1_VGDCS_LAUNCH
#else
      auto mainf_ = KOKKOS_LAMBDA(const int c) {
        int t = usel ? lst_(lo + c) : c;   // (d0 = 0, ndd = n with the list)
        const int d = d0 + t % ndd;
        t /= ndd;
        const int j = js + t % njj;
        t /= njj;
        const int k = ks + t % nkk;
        const int m = t/nkk;
        if (esid_ != 0 && est_(m, k - ks, j - js, d) == esid_) {return;}
        shell_(m, k, j, d, true);
      };
      // cslb: the shell kernel is latency bound at 25 % occupancy (114
      // registers, 4 blocks of 128 per SM; ncu: DRAM 4 %, issue 43 %); env VGD_CS_LB = 6
      // or 8 launches it with LaunchBounds<128, LB> (register cap, same code)
      if (s_lb == 8) {
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 8>>(DevExeSpace(), 0, nper),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      } else if (s_lb == 6) {
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 6>>(DevExeSpace(), 0, nper),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      } else if (s_lb == 2) {
        // on HIP the second bound is the minimum waves per EU, so 6 caps
        // the body at 80 VGPRs + 800 B scratch per thread (rocprofv3 13082414); 2 and 1
        // allow 256 / 512 (register allocation only, the same code)
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 2>>(DevExeSpace(), 0, nper),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      } else if (s_lb == 1) {
        Kokkos::parallel_for("m1_vgdcs_shell", Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::LaunchBounds<128, 1>>(DevExeSpace(), 0, nper),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      } else {
        Kokkos::parallel_for("m1_vgdcs_shell", VgdLight(DevExeSpace(), 0, nper), mainf_);
      }
#endif
      }   // !xffs
      if (esid_ != 0) {eev.StreamWait(DevExeSpace());}   // (the send cells)
#undef M1_VGDCS_LAUNCHN
#if MPI_PARALLEL_ENABLED
      if (pipe) {
        Kokkos::Timer tf;
        xstart(pass, l, g);
        vgd_thalo += tf.seconds();
      }
#endif
      }   // g
      if (rga) {
        // the previous shell of the pass has been read (stream order): its inward band
        // values saved (inward pass), its interior back to the compact slab; the next
        // shell lent a slab (run after this shell's exchange, below)
        if (q > 0) {
          const int ip = inw ? (i + 1) : (i - 1);
          VgdRingOp o;
          if (inw) {
            o.kind = 1;
            o.i = ip;
            vgd_rg_pend.push_back(o);
          }
          o.kind = 2;
          o.i = ip;
          vgd_rg_pend.push_back(o);
        }
        if (q + 1 < nsh) {
          VgdRingOp o;
          o.kind = 3;
          o.i = inw ? (i - 1) : (i + 1);
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
      }
      Kokkos::Timer th;
      if (pipe) {
        // this shell's exchanges are in flight (finished before the next
        // shell's kernels of the same group, or at the end of the pass)
      } else if (early) {
#if MPI_PARALLEL_ENABLED
        // cs early send: the pack has run (event), the shell kernel may still
        // be running; the messages of VetGdCsHaloC, then its unpack (stream order: after
        // the shell kernel)
        const VgdCsCl &E = vgd_cscl[vgd_cscur];
        const int np = static_cast<int>(vgd_csprk.size());
        const size_t sl = static_cast<size_t>(pass)*n1 + l;
        const int64_t *so = &E.so[sl*(np + 1)];
        const int64_t *ro = &E.ro[sl*(np + 1)];
        const int64_t s0 = so[0], ns = so[np] - so[0];
        const int64_t r0 = ro[0], nr = ro[np] - ro[0];
        if (static_cast<int64_t>(vgd_cssb.extent(0)) < ns ||
            static_cast<int64_t>(vgd_csrb.extent(0)) < nr) {
          VgdCsFatal("cs early send: message buffer too small");
        }
        eev.Wait();
        Kokkos::Timer tq;
        std::vector<MPI_Request> req;
        for (int p = 0; p < np; ++p) {
          const int rk = vgd_csprk[p];
          const int64_t cr = ro[p+1] - ro[p], cs = so[p+1] - so[p];
          if (cr > 0) {
            req.emplace_back();
            MPI_Irecv(vgd_csrb.data() + (ro[p] - r0), static_cast<int>(cr),
                      MPI_ATHENA_REAL, rk, 7103, MPI_COMM_WORLD, &req.back());
          }
          if (cs > 0) {
            req.emplace_back();
            MPI_Isend(vgd_cssb.data() + (so[p] - s0), static_cast<int>(cs),
                      MPI_ATHENA_REAL, rk, 7103, MPI_COMM_WORLD, &req.back());
            vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
          }
        }
        MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
        vgd_tmpi += tq.seconds();
        vgd_nexch += 1.0;
        vgd_csvsent += static_cast<double>(ns);
        vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;
        if (nr > 0) {
          auto a_ = vgd_i;
          auto R_ = vgd_csrl;
          auto ri_ = E.ri;
          const Real *rbp = vgd_csrb.data();
          Kokkos::parallel_for("m1_vgdcs_cunpack",
                               VgdLight(DevExeSpace(), 0, static_cast<int>(nr)),
          KOKKOS_LAMBDA(const int o) {
            const int x = ri_(r0 + o);
            const int e = x/n, d = x - e*n;
            a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rbp[o];
          });
        }
#endif
      } else if (ccl) {
        // in a fused sweep both sets' exchanges in one round (one message
        // per partner, VetGdCsHaloCX); else the two exchanges
        if (!(xffs && VetGdCsHaloCX(pass, i, xcurb))) {
          VetGdCsHaloC(pass, i);          // the read values only
          if (xffs) {
            // the next set's exchange (its own compact lists and array)
            const int ca = vgd_cscur;
            vgd_cscur = xcurb;
            std::swap(vgd_i, vgd_i_xb);
            VetGdCsHaloC(pass, i);
            std::swap(vgd_i, vgd_i_xb);
            vgd_cscur = ca;
          }
        }
      } else {
        VetGdCsHalo(vgd_i, n, i, i, false);  // remote band only (on-rank: read direct)
        if (xffs) {VetGdCsHalo(vgd_i_xb, n, i, i, false);}
      }
      if (rga && inw) {
        // the record of the values to save from this shell's inward band: the receive
        // list of its compact exchange (vgd_hxe[i], above)
        vgd_rg_sk[i] = 1;
        vgd_rg_se[i] = i;
      }
      if (rga) {VetGdRingFlush();}
      vgd_thalo += th.seconds();
    }
#if MPI_PARALLEL_ENABLED
    if (pipe) {   // the last shell's exchanges (its band: read by the outward pass)
      Kokkos::Timer tf;
      while (!xfl.empty()) {xfinish();}
      vgd_thalo += tf.seconds();
    }
#endif
    if (rga) {
      // the last shell of the pass: saved (inward) and back to its compact slab
      const int il = is + (inw ? lcut : (n1 - 1));
      if (vgd_rg_state[il] >= 0) {
        VgdRingOp o;
        if (inw) {
          o.kind = 1;
          o.i = il;
          vgd_rg_pend.push_back(o);
        }
        o.kind = 2;
        o.i = il;
        vgd_rg_pend.push_back(o);
      }
      VetGdRingFlush();
    }
  }
  VetGdRingSweepEnd();
  Kokkos::fence();
  if (ccl) {
    auto mh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_csmiss);
    if (mh(0) != 0) {
      VgdCsFatal("compact lists: " + std::to_string(mh(0)) + " remote band reads not "
                 "in the read mask (rank " + std::to_string(global_variable::my_rank) +
                 "); run with vet_gd_halo_compact = 0 (the dense band) and report");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CsTanFrameInit
//! \brief the per-cell orthonormal tangent frame cs_qf of Q: e1 = the normal of the x2
//! plane through the cell centre, e2 = the x3 plane's normal less its e1 part (both
//! planes hold the radial line, so both are tangential).  Cartesian components, exchanged
//! once through the tau_ten kind of exchange (a = at least 6 components, the same
//! no-resample map as the slots): a seam ghost holds its mirror cell's frame, a cube-
//! vertex corner ghost zeros, so Q contracts to zero there.

void RadiationM1::CsTanFrameInit(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                                 MeshBoundaryValuesCC *pb) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n3 = a.extent_int(2), n2 = a.extent_int(3);
  Kokkos::realloc(cs_qf, nmb1 + 1, 6, n3, n2);
  Kokkos::deep_copy(a, 0.0);
  auto a_ = a;
  auto cn2 = cs_cn2, cn3 = cs_cn3;
  par_for("m1_cstan_fr", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real e1[3], e2[3];
    for (int c = 0; c < 3; ++c) {e1[c] = cn2(m,j,c); e2[c] = cn3(m,k,c);}
    const Real d = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
    for (int c = 0; c < 3; ++c) {e2[c] -= d*e1[c];}
    const Real nn = 1.0/sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
    for (int c = 0; c < 3; ++c) {
      a_(m,c,k,j,is) = e1[c];
      a_(m,3+c,k,j,is) = e2[c]*nn;
    }
  });
  VetLatExchange(a, ac, pb);
  auto q_ = cs_qf;
  par_for("m1_cstan_fq", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int c = 0; c < 6; ++c) {q_(m,c,k,j) = a_(m,c,k,j,is);}
  });
  Kokkos::deep_copy(a, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdCsMoments
//! \brief part 0: D_rr = K_rr/J about the cell's r_hat -> LAT0 = D_rr - f_K (LAT1..3 the
//! Cartesian D_r,lat; LAT4, LAT5 = 0, or Q under vet_gd_tangential:);
//! part 1: the vet_gd_replace surface q of the top cell

namespace {
// see VgdMomLds (rad_m1_vetgd.cpp); env VGD_MOM_LDS=0: the par_for kernels
constexpr int VGD_LDS_TS = 64;
constexpr int VGD_LDS_DC = 16;
inline bool VgdMomLdsCs() {
#if defined(KOKKOS_ENABLE_HIP)
  static int on = -1;
  if (on < 0) {
    const char *e = std::getenv("VGD_MOM_LDS");
    on = (e != nullptr) ? std::atoi(e) : 1;
  }
  return on != 0;
#else
  return false;
#endif
}
}  // namespace

void RadiationM1::VetGdCsMoments(const int part) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int ilo = is + vlat_icut;
  auto vi_ = vgd_i;
  const int og = vgd_w - indcs.ng;
  auto dir_ = vgd_dir;
  auto tt_ = tau_ten;
  auto cr_ = vgd_csr;
  if (part == 0) {
    // vet_gd_tangential: Q in LAT4, LAT5 about the frame cs_qf (M1CsQab)
    const bool ctan = vgd_on && vgd_tan;
    auto qf_ = cs_qf;
    if (VgdMomLdsCs()) {
      // LDS-staged (rad_m1_vetgd.cpp VgdMomLds); same cells, sums, outputs
      const int nk = ke - ks + 1, nj = je - js + 1, ni = ie - is + 1;
      const int ncell = (nmb1 + 1)*ni*nj*nk;
      const int nlg = (ncell + VGD_LDS_TS - 1)/VGD_LDS_TS;
      const size_t scr = ScrArray1D<Real>::shmem_size(VGD_LDS_TS*VGD_LDS_DC);
      Kokkos::TeamPolicy<> pol(DevExeSpace(), nlg, VGD_LDS_TS, 1);
      Kokkos::parallel_for("m1_vgdcs_mom", pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
      KOKKOS_LAMBDA(const TeamMember_t &tm) {
        ScrArray1D<Real> sv(tm.team_scratch(0), VGD_LDS_TS*VGD_LDS_DC);
        const int c0 = tm.league_rank()*VGD_LDS_TS;
        const int tr = tm.team_rank();
        auto dec = [&](const int c, int &m, int &i, int &j, int &k) {
          int t = c;
          j = js + t % nj;
          t /= nj;
          k = ks + t % nk;
          t /= nk;
          i = is + t % ni;
          m = t/ni;
        };
        int m = 0, i = 0, j = 0, k = 0;
        const int c = c0 + tr;
        const bool val = (c < ncell);
        if (val) {dec(c, m, i, j, k);}
        const bool act = val && (i >= ilo);
        if (val) {
          for (int cc = 0; cc < M1_TT_NLAT; ++cc) {tt_(m,M1_TT_LAT0+cc,k,j,i) = 0.0;}
        }
        Real rx = 0.0, ry = 0.0, rz = 0.0;
        if (act) {
          rx = cr_(m,k+og,j+og,0); ry = cr_(m,k+og,j+og,1);
          rz = cr_(m,k+og,j+og,2);
        }
        Real jm = 0.0, rr = 0.0, kx = 0.0, ky = 0.0, kz = 0.0;
        Real q1m2 = 0.0, q12 = 0.0;
        for (int d0 = 0; d0 < n; d0 += VGD_LDS_DC) {
          const int dn = (n - d0 < VGD_LDS_DC) ? (n - d0) : VGD_LDS_DC;
          tm.team_barrier();
          for (int e = tr; e < VGD_LDS_TS*dn; e += VGD_LDS_TS) {
            const int cc = e/dn, dd = e - (e/dn)*dn;
            const int c2 = c0 + cc;
            if (c2 < ncell) {
              int m2, i2, j2, k2;
              dec(c2, m2, i2, j2, k2);
              if (i2 >= ilo) {sv(cc*VGD_LDS_DC + dd) = vi_(m2,d0+dd,k2+og,j2+og,i2);}
            }
          }
          tm.team_barrier();
          if (act) {
            for (int dd = 0; dd < dn; ++dd) {
              const int d = d0 + dd;
              const Real vv = sv(tr*VGD_LDS_DC + dd);
              const Real a = dir_(d,0)*rx + dir_(d,1)*ry + dir_(d,2)*rz;
              const Real wi = dir_(d,3)*vv;
              jm += wi;
              rr += wi*a*a;
              kx += wi*a*dir_(d,0);
              ky += wi*a*dir_(d,1);
              kz += wi*a*dir_(d,2);
              if (ctan) {
                const Real b1 = dir_(d,0)*qf_(m,0,k,j) + dir_(d,1)*qf_(m,1,k,j)
                                + dir_(d,2)*qf_(m,2,k,j);
                const Real b2 = dir_(d,0)*qf_(m,3,k,j) + dir_(d,1)*qf_(m,4,k,j)
                                + dir_(d,2)*qf_(m,5,k,j);
                q1m2 += wi*(b1*b1 - b2*b2);
                q12 += wi*b1*b2;
              }
            }
          }
        }
        if (act && jm > 0.0) {
          if (ctan) {
            tt_(m,M1_TT_LAT0+4,k,j,i) = 0.5*q1m2/jm;
            tt_(m,M1_TT_LAT0+5,k,j,i) = q12/jm;
          }
          tt_(m,M1_TT_LAT0,k,j,i) = rr/jm - tt_(m,0,k,j,i);
          const Real dr = rr/jm;
          tt_(m,M1_TT_LAT0+1,k,j,i) = kx/jm - dr*rx;
          tt_(m,M1_TT_LAT0+2,k,j,i) = ky/jm - dr*ry;
          tt_(m,M1_TT_LAT0+3,k,j,i) = kz/jm - dr*rz;
        }
      });
      return;
    }
    par_for("m1_vgdcs_mom", DevExeSpace(), 0, nmb1, is, ie, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
      for (int c = 0; c < M1_TT_NLAT; ++c) {tt_(m,M1_TT_LAT0+c,k,j,i) = 0.0;}
      if (i < ilo) {return;}
      const Real rx = cr_(m,k+og,j+og,0), ry = cr_(m,k+og,j+og,1);
      const Real rz = cr_(m,k+og,j+og,2);
      Real jm = 0.0, rr = 0.0, kx = 0.0, ky = 0.0, kz = 0.0;
      Real q1m2 = 0.0, q12 = 0.0;
      for (int d = 0; d < n; ++d) {
        const Real a = dir_(d,0)*rx + dir_(d,1)*ry + dir_(d,2)*rz;
        const Real wi = dir_(d,3)*vi_(m,d,k+og,j+og,i);
        jm += wi;
        rr += wi*a*a;
        kx += wi*a*dir_(d,0);
        ky += wi*a*dir_(d,1);
        kz += wi*a*dir_(d,2);
        if (ctan) {
          const Real b1 = dir_(d,0)*qf_(m,0,k,j) + dir_(d,1)*qf_(m,1,k,j)
                          + dir_(d,2)*qf_(m,2,k,j);
          const Real b2 = dir_(d,0)*qf_(m,3,k,j) + dir_(d,1)*qf_(m,4,k,j)
                          + dir_(d,2)*qf_(m,5,k,j);
          q1m2 += wi*(b1*b1 - b2*b2);
          q12 += wi*b1*b2;
        }
      }
      if (jm > 0.0) {
        if (ctan) {
          tt_(m,M1_TT_LAT0+4,k,j,i) = 0.5*q1m2/jm;
          tt_(m,M1_TT_LAT0+5,k,j,i) = q12/jm;
        }
        tt_(m,M1_TT_LAT0,k,j,i) = rr/jm - tt_(m,0,k,j,i);
        // v = D.rhat - D_rr rhat (CARTESIAN, tangential) -> LAT1..LAT3, the
        // D_r,lat of the cs operator (M1CsLat0 / M1CsLatF); read only with
        // vet_col_lat_offdiag = lagged | operator
        const Real dr = rr/jm;
        tt_(m,M1_TT_LAT0+1,k,j,i) = kx/jm - dr*rx;
        tt_(m,M1_TT_LAT0+2,k,j,i) = ky/jm - dr*ry;
        tt_(m,M1_TT_LAT0+3,k,j,i) = kz/jm - dr*rz;
      }
    });
    return;
  }
  const bool sqf = vcol_sqf;
  const Real jca = vcol_jca, jcb = vcol_jcb, jcap = vcol_jcap, q0 = marshak_q;
  const Real qlo = vcol_qmin, qhi = vcol_qmax;
  auto vq_ = vcol_q;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  par_for("m1_vgdcs_q", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const VgdRhatCs rh{cr_(m,k+og,j+og,0), cr_(m,k+og,j+og,1), cr_(m,k+og,j+og,2)};
    vq_(m,k,j) = VgdSurfaceQ(rh, dir_, vi_, n, m, k, j, og, ie, cx1v(m,ie)/cx1f(m,ie+1),
                             sqf, jca, jcb, jcap, qlo, qhi, q0);
  });
}

} // namespace radm1
