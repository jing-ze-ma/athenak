//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetgd_cs.cpp
//! \brief <rad_vet>/vet_gd on the cubed sphere: the sp sweep of rad_vet_vetgd.cpp
//! (shells, passes, segment and boundary rules) with the cs lateral stencil and band.
//! Band: the own panel's index space extended by vgd_w cells; beyond a seam the
//! neighbour's cells are stored unfolded (mirror image across the seam). No cell beyond
//! two edges or at depth >= N. Filled per shell by exact copies keyed on the global cell.
//! Read: the upwind point is located on its panel at fractional native indices:
//! bilinear inside the centre lattice; linear across a seam strip, each edge column
//! along the seam in its own chart; barycentric on the 3 cells at a cube vertex.
//! Weights >= 0, sum 1. Hand-over: D_rr - f_K into LAT0, the Cartesian tangential
//! vector D.rhat - D_rr rhat into LAT1..LAT3 (no tangential anisotropy on cs).

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
#include <thread>  // NOLINT(build/c++11)
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
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_vetseg.hpp"
#include "rad_vet/rad_vet_vgdray.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

// build generation per compact-list cache entry; not a member, so the dt-retry
// checkpoint keeps it and DtRetryAfterRestore drops entries rebuilt by the failed try
int64_t g_vgd_csgen[16] = {0};        // VET_VGD_CSCL_MAX
int64_t g_vgd_csgen_save[16] = {0};

namespace {

// stream event on the sweep's stream: record after a kernel, query from the host
// without a full fence (no device backend: fence)
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
  VgdFatalAt(__FILE__, "<rad_vet>/vet_gd (cubed sphere): ", msg);
}

KOKKOS_INLINE_FUNCTION
Real VgdCsDet3(const Real *a, const Real *b, const Real *c) {
  return a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0])
         + a[2]*(b[0]*c[1] - b[1]*c[0]);
}

// upwind stencil: up to 4 band positions (k, j) and their weights
struct VgdCsSt {
  int n;
  int k[4], j[4];
  Real w[4];
};

// stencil accessor of the shared upwind core (rad_vet_vgdray.hpp): weighted sum
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

// stencil of unit vector x for a block on panel p with first interior cell (j0, k0):
// band index = unfolded index - (j0, k0) + wb, extents c2 x c3. false if a cell is not
// in the band (dropped, rest renormalised; the caller counts a clamped read).
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
    // inside c's centre lattice: bilinear in native indices
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
    // seam strip: c's edge column and the neighbour's, linear across the seam
    const bool sa = !ina;                    // seam normal is c's a axis
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
    const bool na = (oa >= ot);             // seam normal is c2p's a axis
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
    // vertex corner: the three cells around the cube vertex, barycentric
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
  // global cells -> band positions of this block
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
//! \fn void RadiationVet::VetGdCsInit
//! \brief cs band setup (from VetGdInit): depth vgd_w, panel tables, unit vector of
//! every band cell, exact exchange lists

void RadiationVet::VetGdCsInit() {
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
  // blocks: panel and global first cell; a block's band comes from its own x1 level
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
  // unfolded position (J, K) of panel p -> the global cell (q, a, t) it holds
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
  // inverse maps (q's (a, t) -> p's unfolded (J, K)): signed permutation + shift
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
  // check: forward unfold and inverse map agree on every strip cell
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
  // band depth: sp formula with the narrowest angular spacing (+3: one more than sp for
  // the index kink across a seam); capped at N
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
      // shells next to another level: radii from the x1 ghosts
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
      std::cout << "<rad_vet> vet_gd (cubed sphere): N = " << nc << " per panel, "
                << "narrowest "
                << "cell spacing " << dmin << " rad, dense lateral band " << vgd_w
                << " cells (needed " << wneed << ")"
                << (vgd_capped ? ": capped at N, near-tangent reads beyond are clamped "
                                 "(counted)" : "") << std::endl;
    }
  }
  const int w = vgd_w;
  const int c2w = nx2 + 2*w, c3w = nx3 + 2*w;
  // device tables
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
  // unit vector of the cell at every band position (0: none)
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
  // exchange lists: receiver blocks in global order, band positions in (k, j) order;
  // the same iteration on every rank, so send and receive lists match in order
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
  // owner of every band position for direct on-rank reads; remote / none: itself
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
  // message buffers: one shell of every direction (source exchange chunked)
  const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
  const int nvb = std::max(vgd_n, 2*8);
  Kokkos::realloc(vgd_cssb, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_csrb, static_cast<size_t>(nent)*nvb);
  vgd_hmpi = (np > 0);
  // every side of every block carries the band: no direct neighbour reads
  Kokkos::realloc(vgd_hloc, 8*nmb);
  Kokkos::deep_copy(vgd_hloc, -1);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  // compact lists (vet_gd_halo_compact > 0, MPI band): receive entry of every band
  // position (-1: own interior, on-rank owner or no cell)
  vgd_cscomp = (vgd_hcomp > 0) && vgd_hmpi;
  Kokkos::realloc(vgd_csrix, nmb, c3w, c2w);
  {
    auto x_h = Kokkos::create_mirror_view(vgd_csrix);
    for (size_t q = 0; q < x_h.size(); ++q) {x_h.data()[q] = -1;}
    for (int e = 0; e < vgd_csnrcv; ++e) {x_h(rl[3*e], rl[3*e+1], rl[3*e+2]) = e;}
    Kokkos::deep_copy(vgd_csrix, x_h);
  }
  Kokkos::realloc(vgd_csmiss, 1);
  // one cache entry per set of vet_rotate_avg (the K sets cycle); entries beyond
  // vgd_csncl stay unallocated
  vgd_csncl = 2;
  if (vgd_cscomp) {
    vgd_csncl = std::min(std::max(2, vra_k), static_cast<int>(VET_VGD_CSCL_MAX));
  }
  vgd_csper = static_cast<int64_t>(vgd_n)*indcs.nx1;
  vgd_csper = 16*((vgd_csper + 15)/16);   // padded for the packed-word lists
  for (int c = 0; c < vgd_csncl; ++c) {
    vgd_cscl[c].alpha = -9.0;
    vgd_cscl[c].use = 0;
    vgd_cscl[c].msk = decltype(vgd_cscl[c].msk)();   // masks live in vgd_csmb
  }
  if (me == 0) {
    std::cout << "<rad_vet> vet_gd (cubed sphere): rank 0 band copies on-rank "
              << vgd_csnloc << ", sent " << vgd_csnsnd << ", received " << vgd_csnrcv
              << " (cells per shell), partners " << np << "; MPI band "
              << (vgd_cscomp ? "COMPACT (vet_gd_halo_compact > 0: per shell and pass "
                               "only the read (cell, direction) values)"
                             : (vgd_hmpi ? "dense (every remote band cell, every "
                                           "direction, per shell and pass)"
                                         : "none"))
              << "; list cache " << vgd_csncl << " sets" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsCompact
//! \brief compact lists of the current direction set, LRU cache of
//! max(2, vet_rotate_avg) entries. (1) read mask: replay the sweep's upwind geometry and
//! mark every remote band value read (1 inward, 2 outward; both cut variants, so the
//! lists depend on the set alone); (2) receivers send marks to owners; (3) both sides
//! build lists in the same (pass, shell, partner, entry, direction) order.

void RadiationVet::VetGdCsCompact() {
  int c = -1;
  for (int q = 0; q < vgd_csncl; ++q) {
    if (vgd_cscl[q].alpha == vgd_alpha) {c = q;}
  }
  if (c >= 0 && vgd_cscl[c].msk.data() != nullptr) {
    vgd_cscur = c;
    vgd_cscl[c].use = ++vgd_csuse;
    return;
  }
  // a cached set whose mask was reused gets its mask again (step 1 only)
  const bool remask = (c >= 0);
  Kokkos::Timer tb;
  if (!remask) {
    c = 0;   // least recently used entry (first of equal ones)
    for (int q = 1; q < vgd_csncl; ++q) {
      if (vgd_cscl[q].use < vgd_cscl[c].use) {c = q;}
    }
  }
  VgdCsCl &E = vgd_cscl[c];
  // two read masks only (vgd_csmb); the mask of the entry in use is kept
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
  const int64_t per_ = vgd_csper;   // n*n1 padded to 16
  Kokkos::Timer tsub;
  // (1) read mask (VetGdCsSweep's geometry; values never matter)
  par_for("vet_vgdcs_rmask", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, n - 1,
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
  // (3) the lists: counts per (pass, shell, partner), offsets, fill in (entry, d) order
  const size_t nsl = static_cast<size_t>(2)*n1;
  std::vector<int> sl_h, rl_h;
  int64_t ts = 0, tr = 0;
  if (per % 16 != 0) {VgdCsFatal("compact lists: n*n1 not padded to 16");}
  // lists from packed 2-bit words (partner blocks word aligned)
  const size_t wpe = per/16;   // words per entry
  // only the nonzero (word index, word) pairs, in increasing word order, are downloaded,
  // exchanged and scanned. Host copies use persistent buffers (grown, never freed): a
  // fresh buffer per build let UCX's registration cache send stale pages under
  // LD_PRELOAD profilers (MPI_ERR_TRUNCATE or hang).
  const int64_t nwr = static_cast<int64_t>(vgd_csnrcv)*static_cast<int64_t>(wpe);
  using RpL = Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>;
  int64_t nz = 0;
  Kokkos::parallel_scan("vet_vgdcs_nzn", RpL(DevExeSpace(), 0, nwr),
  KOKKOS_LAMBDA(const int64_t w, int64_t &u, const bool) {
    if (msk_(w) != 0u) {u += 1;}
  }, nz);
  DvceArray1D<int64_t> zi("vet_vgdcs_zi", std::max<int64_t>(nz, 1));
  DvceArray1D<uint32_t> zv("vet_vgdcs_zv", std::max<int64_t>(nz, 1));
  Kokkos::parallel_scan("vet_vgdcs_nzc", RpL(DevExeSpace(), 0, nwr),
  KOKKOS_LAMBDA(const int64_t w, int64_t &u, const bool fin) {
    const uint32_t x = msk_(w);
    if (x == 0u) {return;}
    if (fin) {zi(u) = w; zv(u) = x;}
    u += 1;
  });
  std::vector<int64_t> hzi(std::max<int64_t>(nz, 1));
  std::vector<uint32_t> hzv(std::max<int64_t>(nz, 1));
  Kokkos::deep_copy(Kokkos::View<int64_t *, Kokkos::HostSpace,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(hzi.data(), hzi.size()), zi);
  Kokkos::deep_copy(Kokkos::View<uint32_t *, Kokkos::HostSpace,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(hzv.data(), hzv.size()), zv);
  bool both = false;
  for (int64_t q = 0; q < nz; ++q) {
    const uint32_t x = hzv[q];
    if ((x & (x >> 1) & 0x55555555u) != 0u) {both = true;}
  }
  if (both) {VgdCsFatal("read mask: a value marked both inward and outward");}
  // own pairs per partner, into the persistent buffer
  std::vector<size_t> ho(np + 1, 0), hn(np, 0);
  std::vector<int64_t> hlo(np, 0);
  for (int p = 0; p < np; ++p) {
    const int64_t w0 = static_cast<int64_t>(vgd_cspr[p])*wpe;
    const int64_t w1 = w0 + static_cast<int64_t>(vgd_cspm[p])*wpe;
    if (w1 - w0 > static_cast<int64_t>(UINT32_MAX)) {
      VgdCsFatal("compact lists: a partner block of more than 2^32 mask words");
    }
    hlo[p] = std::lower_bound(hzi.begin(), hzi.begin() + nz, w0) - hzi.begin();
    hn[p] = (vgd_cspm[p] > 0)
            ? (std::lower_bound(hzi.begin(), hzi.begin() + nz, w1) - hzi.begin()) - hlo[p]
            : 0;
    ho[p+1] = ho[p] + 2*hn[p];
  }
  if (vgd_csmw_h.size() < std::max<size_t>(ho[np], 1)) {vgd_csmw_h.resize(ho[np]);}
  uint32_t *mp = vgd_csmw_h.data();
  for (int p = 0; p < np; ++p) {
    const int64_t w0 = static_cast<int64_t>(vgd_cspr[p])*wpe;
    for (size_t q = 0; q < hn[p]; ++q) {
      mp[ho[p] + 2*q] = static_cast<uint32_t>(hzi[hlo[p] + q] - w0);
      mp[ho[p] + 2*q + 1] = hzv[hlo[p] + q];
    }
  }
  vgd_cstb_d += tsub.seconds();
  tsub.reset();
  // owners' pairs (marks of the entries they send me), per partner
  std::vector<uint64_t> rn(np, 0), sn(np, 0);
  for (int p = 0; p < np; ++p) {sn[p] = hn[p];}
  std::vector<size_t> ro(np + 1, 0);
#if MPI_PARALLEL_ENABLED
  {
    // tags 7109 pair counts, 7107 pairs, 7108 checksums (7105: fused twin exchange);
    // a checksum mismatch is fatal
    auto hsum = [](const uint32_t *a, const size_t nw) {
      uint64_t h = 1469598103934665603ull;
      for (size_t w = 0; w < nw; ++w) {h = (h ^ a[w])*1099511628211ull;}
      return h;
    };
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_csprk[p];
      if (vgd_cspn[p] > 0) {
        req.emplace_back();
        MPI_Irecv(&rn[p], 1, MPI_UINT64_T, rk, 7109, MPI_COMM_WORLD, &req.back());
      }
      if (vgd_cspm[p] > 0) {
        req.emplace_back();
        MPI_Isend(&sn[p], 1, MPI_UINT64_T, rk, 7109, MPI_COMM_WORLD, &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    req.clear();
    for (int p = 0; p < np; ++p) {ro[p+1] = ro[p] + 2*rn[p];}
    if (vgd_csmw_r.size() < std::max<size_t>(ro[np], 1)) {vgd_csmw_r.resize(ro[np]);}
    std::vector<uint64_t> hs(np, 0), hr(np, 0);
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_csprk[p];
      if (vgd_cspn[p] > 0) {
        req.emplace_back();
        MPI_Irecv(vgd_csmw_r.data() + ro[p], static_cast<int>(2*rn[p]), MPI_UINT32_T,
                  rk, 7107, MPI_COMM_WORLD, &req.back());
        req.emplace_back();
        MPI_Irecv(&hr[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
      }
      if (vgd_cspm[p] > 0) {
        hs[p] = hsum(mp + ho[p], 2*hn[p]);
        req.emplace_back();
        MPI_Isend(mp + ho[p], static_cast<int>(2*hn[p]), MPI_UINT32_T, rk, 7107,
                  MPI_COMM_WORLD, &req.back());
        req.emplace_back();
        MPI_Isend(&hs[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    for (int p = 0; p < np; ++p) {
      if (vgd_cspn[p] <= 0) {continue;}
      if (hsum(vgd_csmw_r.data() + ro[p], 2*rn[p]) != hr[p]) {
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
  // each partner's counts and fills on its own host thread (disjoint ranges)
  auto parts = [np](const auto &f) {
    std::vector<std::thread> th;
    for (int p = 1; p < np; ++p) th.emplace_back(f, p);
    if (np > 0) f(0);
    for (auto &t : th) t.join();
  };
  // lists from partner p's pairs (pr + po[p], pn[p] of them), entries from st0[p]
  auto buildw = [&](const uint32_t *pr, const std::vector<size_t> &po,
                    const std::vector<int> &st0, std::vector<int64_t> &off,
                    std::vector<int> &lst) {
    std::vector<int64_t> nn(nsl*np, 0);
    parts([&](const int p) {
      const size_t q1 = (po[p+1] - po[p])/2;
      for (size_t q = 0; q < q1; ++q) {
        uint32_t x = pr[po[p] + 2*q + 1];
        const size_t r0 = static_cast<size_t>(pr[po[p] + 2*q])*16;   // in-block t
        while (x != 0u) {
          const int b = __builtin_ctz(x) >> 1;
          const uint32_t v = (x >> (2*b)) & 3u;
          x &= ~(3u << (2*b));
          const size_t r = (r0 + b) % per;
          const int sh = static_cast<int>(r % n1);
          nn[(static_cast<size_t>(v - 1)*n1 + sh)*np + p] += 1;
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
      const size_t q1 = (po[p+1] - po[p])/2;
      for (size_t q = 0; q < q1; ++q) {
        uint32_t x = pr[po[p] + 2*q + 1];
        const size_t r0 = static_cast<size_t>(pr[po[p] + 2*q])*16;
        while (x != 0u) {   // bits low to high = t increasing
          const int b = __builtin_ctz(x) >> 1;
          const uint32_t v = (x >> (2*b)) & 3u;
          x &= ~(3u << (2*b));
          const size_t rr = r0 + b;
          const int e = st0[p] + static_cast<int>(rr/per);
          const size_t r = rr % per;
          const int d = static_cast<int>(r/n1), sh = static_cast<int>(r % n1);
          lst[cur[(static_cast<size_t>(v - 1)*n1 + sh)*np + p]++] = e*n + d;
        }
      }
    });
    return t;
  };
  // send and receive lists concurrently
  std::thread tl([&]() {tr = buildw(mp, ho, vgd_cspr, E.ro, rl_h);});
  ts = buildw(vgd_csmw_r.data(), ro, vgd_csps, E.so, sl_h);
  tl.join();
  vgd_cstb_h += tsub.seconds();
  tsub.reset();
  // (entry, direction) order is the same on both sides; counts per (pass, shell,
  // partner) are checked by the message sizes below
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
  E.use = ++vgd_csuse;
  vgd_cscur = c;
  g_vgd_csgen[c] += 1;
  vgd_csnbuild += 1.0;
  vgd_cstbuild += tb.seconds();
  if (global_variable::my_rank == 0 && vgd_csnbuild <= 1.0) {
    const double dn = static_cast<double>(vgd_csnsnd)*n*n1*2.0;
    std::cout << "<rad_vet> vet_gd (cubed sphere) compact lists: rank 0 sends " << ts
              << " values per sweep (dense " << dn << ": "
              << ((dn > 0.0) ? ts/dn : 0.0) << "), receives " << tr << "; built in "
              << tb.seconds() << " s (once per direction set)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsHaloC
//! \brief compact exchange of shell i after pass `pass`: the listed (entry, direction)
//! values of every partner (VetGdCsCompact)

void RadiationVet::VetGdCsHaloC(const int pass, const int i) {
#if MPI_PARALLEL_ENABLED
  const VgdCsCl &E = vgd_cscl[vgd_cscur];
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1;
  const int n = vgd_n;
  const int np = static_cast<int>(vgd_csprk.size());
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
  Real *sbp = sb_.data();
  Real *rbp = rb_.data();
  if (ns > 0) {
    Kokkos::parallel_for("vet_vgdcs_cpack",
                         VgdLight(DevExeSpace(), 0, static_cast<int>(ns)),
    KOKKOS_LAMBDA(const int o) {
      const int x = si_(s0 + o);
      const int e = x/n, d = x - e*n;
      sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
    });
  }
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
  vgd_nexch += 1.0;
  vgd_csvsent += static_cast<double>(ns);
  vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;
  if (nr > 0) {
    Kokkos::parallel_for("vet_vgdcs_cunpack",
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

// per-partner offsets of the merged two-set exchange (VetGdCsHaloCX)
namespace {
constexpr int VGD_CSX_NS = 32;
struct VgdCsXIdx {int v[VGD_CSX_NS + 1];};
}   // namespace
//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdCsHaloCX
//! \brief both compact exchanges of shell i of a fused two-set sweep in one round (one
//! message per partner: current set's piece, then the next set's; tag 7107); bitwise two
//! VetGdCsHaloC calls. false (nothing done): the caller runs the two exchanges.

bool RadiationVet::VetGdCsHaloCX(const int pass, const int i, const int eb) {
#if MPI_PARALLEL_ENABLED
  if (eb < 0) {return false;}
  const VgdCsCl &EA = vgd_cscl[vgd_cscur];
  const VgdCsCl &EB = vgd_cscl[eb];
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
  // buffers hold both sets (grown once if needed)
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
    Kokkos::parallel_for("vet_vgdcs_cpackx", VgdLight(DevExeSpace(), 0, nsa + nsb),
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
    Kokkos::parallel_for("vet_vgdcs_cunpackx", VgdLight(DevExeSpace(), 0, nra + nrb),
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
//! \fn void RadiationVet::VetGdCsHalo
//! \brief fill every band position of a(m, v, k, j, i), v < nv, shells [i0, i1] with its
//! owner's interior value (exact copies, chunked). ws >= 0 (vet_gd_band_global,
//! depth-sorted lists): only positions within depth ws; ws < 0 (cubed sphere): all

template <class V>
void RadiationVet::VetGdCsHalo(V &a, const int nv, const int i0, const int i1,
                               const bool loc, const int ws) {
  auto a_ = a;
  const bool dp = (ws >= 0) && (vgd_gband || vgd_raw) && !vgd_gdnl.empty();
  const int wq = dp ? std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1) : 0;
  const int nl = loc ? (dp ? vgd_gdnl[wq] : vgd_csnloc) : 0;
  {
    const int ni = i1 - i0 + 1;
    if (nl > 0 && ni > 0) {
      auto L_ = vgd_csloc;
      par_for("vet_vgdcs_loc", DevExeSpace(), 0, nl - 1, 0, nv - 1, 0, ni - 1,
      KOKKOS_LAMBDA(const int e, const int v, const int ii) {
        a_(L_(e,0),v,L_(e,1),L_(e,2),i0+ii) = a_(L_(e,3),v,L_(e,4),L_(e,5),i0+ii);
      });
    }
  }
#if MPI_PARALLEL_ENABLED
  if (vgd_hmpi) {
    const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
    // full-sphere global band: device buffers (CUDA-aware MPI); cubed sphere: pinned host
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
        // depth-ws prefix of every partner segment (same offsets)
        par_for("vet_vgdgb_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DS_(e) > wq) {return;}
          gsb_(static_cast<size_t>(e)*nvi + v*ni + ii) =
              a_(S_(e,0),v,S_(e,1),S_(e,2),c0+ii);
        });
      } else if (vgd_csnsnd > 0) {
        par_for("vet_vgdcs_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
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
        par_for("vet_vgdgb_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0,
                ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DR_(e) > wq) {return;}
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              grb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      } else if (vgd_csnrcv > 0) {
        par_for("vet_vgdcs_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0,
                ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              rb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      }
    }
  }
#endif
}

template void RadiationVet::VetGdCsHalo<DvceArray5D<Real>>(DvceArray5D<Real> &a,
                                                         const int nv, const int i0,
                                                         const int i1, const bool loc,
                                                         const int ws);
template void RadiationVet::VetGdCsHalo<RadiationVet::VgdRag>(RadiationVet::VgdRag &a,
                                                              const int nv, const int i0,
                                                              const int i1,
                                                              const bool loc,
                                                              const int ws);

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsSource
//! \brief band of ln chi, ln S (vgd_cs) for every shell of the sweep

void RadiationVet::VetGdCsSource() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ilo = indcs.is, ie = indcs.ie;
  VetGdCsHalo(vgd_cs, 2, ilo, ie, true);
  Kokkos::fence();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsTab
//! \brief host tables of direction set d on cs: no walls, identity map, branch n . r_hat
//! of every band cell

void RadiationVet::VetGdCsTab(const std::vector<double> &d) {
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
//! \fn void RadiationVet::VetGdCsSweep
//! \brief one sweep (inward then outward) of every direction on the cubed sphere:
//! the VetGdSweep body with the cs upwind stencil

void RadiationVet::VetGdCsSweep() {
  if (vgd_xff_prep) {return;}   // nothing to prepare (lists made below)
  // vet_gd_rotate_xfade_fuse: the next set (own intensities vgd_i_xb, tables, lists)
  // is swept in the same shell kernels
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
  auto cnt_ = vray_cnt;
  auto mr_ = vgd_mr;
  auto cr_ = vgd_csr;
  auto bl_ = vgd_csblk;
  auto cm_ = vgd_cscm;
  auto nb_ = vgd_csnb;
  auto ow_ = vgd_csown;
  const int nc = vgd_csn;
  const Real h = 0.5*M_PI/nc;
  const bool bfc = vet_bc_face;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  const bool emi = VetBcEmit() && (lcut == 0) &&
                   (ibc_x1min != VET_IBC_REFLECT);
  const Real iemit = VetBathLo();
  const bool tbh = vgd_tbath && (ibc_x1max == VET_IBC_MARSHAK) && (iebath_x1max != 0.0);
  const Real itop = tbh ? iebath_x1max : 0.0;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const Real lcr = log(kVetCornerRatio);
  vgd_hsweep += xffs ? 2 : 1;
  // compact lists of this set (built on first sweep); a remote read missing from the
  // read mask is fatal below
  const bool ccl = vgd_cscomp;
  int xcurb = -1;   // next set's compact-list entry
  if (ccl) {
    VetGdCsCompact();
    if (xffs) {
      // current set first, then the next (the order of two sweeps: same LRU entry)
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
  const int64_t per_ = vgd_csper;   // mask stride per entry
  // vet_gd_band_ring (needs compact lists): the sp ring; the saved inward band values of
  // shell i are the receive list of its inward compact exchange (record vgd_hxe[i];
  // all addresses in the unused entry vgd_hxe[n1c]). Pure copies (bitwise).
  const bool rga = vgd_ring && ccl;
  if (vgd_ring && !ccl) {VgdCsFatal("vet_gd_band_ring without the compact lists");}
  // next set's records in vgd_hxe_alt; the ring's second array is its intensities
  auto ringrec = [&](std::vector<VgdHlEntry> &hxe, const VgdCsCl &E, const Real alr) {
    const int n1c = indcs.nx1 + 2*indcs.ng;
    if (static_cast<int>(hxe.size()) != 2*n1c) {hxe.resize(2*n1c);}
    const int np = static_cast<int>(vgd_csprk.size());
    const int64_t a0 = E.ro[0];
    const int64_t a1 = E.ro[static_cast<size_t>(n1 - 1)*(np + 1) + np];
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
      par_for("vet_vgdcs_ringrec", DevExeSpace(), 0, static_cast<int>(na) - 1,
      KOKKOS_LAMBDA(const int o) {
        const int x = ri_(a0 + o);
        const int e = x/n, d = x - e*n;
        bd_(o) = ((R_(e,0)*n + d)*c3b + R_(e,1))*c2b + R_(e,2);
      });
    }
    for (int l = 0; l < n1; ++l) {
      VgdHlEntry &X = hxe[is + l];
      const int64_t r0 = E.ro[static_cast<size_t>(l)*(np + 1)];
      const int64_t r1 = E.ro[static_cast<size_t>(l)*(np + 1) + np];
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
  if (rga) {vgd_cur = DevExeSpace();}   // ring copy kernel on this sweep's stream
  // per pass: codes ((m k) j) d of the (cell, direction) pairs of that pass's branch
  // (inward: n.r_hat < 0), in flat order, made each sweep from vgd_mr; the shell kernel
  // runs over this list only
  const int nkl = ke - ks + 1, njl = je - js + 1;
  const int nlt = (nmb1 + 1)*nkl*njl*n;
  int nlist[2] = {0, 0};
  if (vgd_cslst.extent_int(0) < 2*nlt) {
    vgd_cslst = DvceArray1D<int>();
    Kokkos::realloc(vgd_cslst, 2*nlt);
  }
  auto L_ = vgd_cslst;
  for (int p = 0; p < 2; ++p) {
    const bool inp = (p == 0);
    const int off = p*nlt;
    int cnt = 0;
    Kokkos::parallel_scan("vet_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                                  nlt),
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
  auto lst_ = vgd_cslst;
  // next set's pass lists (same scan over its own n . r_hat)
  int nlistb[2] = {0, 0};
  if (xffs) {
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
      Kokkos::parallel_scan("vet_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                                    nlt),
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
  bool early = ccl;
#if !MPI_PARALLEL_ENABLED
  early = false;
#endif
  if (xffs) {early = false;}   // early send: one set per shell kernel
#if defined(KOKKOS_ENABLE_HIP) && MPI_PARALLEL_ENABLED
  // HIP: > 2 ranks per GPU hang in the early send's event wait unless
  // GPU_MAX_HW_QUEUES <= 2 (same bits)
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
                 "use at most 2 ranks per GPU");
    }
  }
#endif
  // early-send cells are computed once: vet_vgdcs_estamp stamps the (m, k, j, d) that
  // vet_vgdcs_eshell computes (stamp vgd_esid); the shell kernel skips them (bitwise)
  const int nkb = ke - ks + 1, njb = je - js + 1;
  if (early && (vgd_est[0].extent_int(0) < nmb1 + 1 || vgd_est[0].extent_int(1) < nkb ||
                vgd_est[0].extent_int(2) < njb || vgd_est[0].extent_int(3) < n)) {
    vgd_est[0] = DvceArray4D<int>("vgd_est", nmb1 + 1, nkb, njb, n);   // zeroed
    vgd_esid = 0;
  }
  auto est_ = vgd_est[0];
  int esid = 0;   // 0: no early cells (never a stamp)
  VgdCsEv eev, efork;
  // shell kernel LaunchBounds<128, kLb> (register cap): HIP 2 (min waves per EU; more
  // spills), CUDA 8
#if defined(KOKKOS_ENABLE_HIP)
  constexpr unsigned kLb = 2;
#else
  constexpr unsigned kLb = 8;
#endif
  VetGdRingSweepBegin(xffs);   // the other set's array second
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    const int nsh = n1 - lcut;
    for (int q = 0; q < nsh; ++q) {
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      if (rga) {
        // this shell's wide slab: inward from scratch, outward with its interior and
        // its inward-exchange band values
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
      // flat (m, k, j, d) range, d fastest; with the fused twin one thread computes the
      // geometry once and the values of both fields
      const int nkk = ke - ks + 1, njj = je - js + 1;
      const int nper = nlist[pass];
      const int lo = pass*nlt;
      // shell body for one (cell, direction): run over the pass's list and (early send)
      // first over the compact send list with cnt = false (counters counted once)
#define VET_VGDCS_SHELL_NAME shell_
#include "rad_vet_vetgd_cs_shell.inc"
#undef VET_VGDCS_SHELL_NAME
      // early send: compute and pack this shell's send values first, record an event;
      // the full shell kernel overlaps the host exchange (finished at the halo point)
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
          if (vgd_esid == 1) {Kokkos::deep_copy(est_, 0);}   // first use or wrapped
          esid = vgd_esid;
          const int sid_ = esid;
          Kokkos::parallel_for("vet_vgdcs_estamp", VgdLight(DevExeSpace(), 0,
                                                            static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            est_(S_(e,0), S_(e,1) - ks, S_(e,2) - js, d) = sid_;
          });
          // send cells + pack on the second stream, concurrent with the shell kernel;
          // the sweep stream waits for them after it (efork / eev)
          VgdEs2Make();
          efork.Record();
          efork.StreamWait(vgd_es2);
          Kokkos::parallel_for("vet_vgdcs_eshell",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            shell_(S_(e,0), S_(e,1), S_(e,2), d, true);
          });
          auto a_ = vgd_i;
          Real *sbp = vgd_cssb.data();
          Kokkos::parallel_for("vet_vgdcs_cpack",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
          });
        }
        eev.Record((esid != 0) ? vgd_es2 : DevExeSpace());
      }
      // fused two-set shell kernel: this pass's list of the current set, then the next
      // set's, one launch (next set's body: the include with its names)
#define VET_VGDCS_LAUNCHN(F, N) \
      Kokkos::parallel_for("vet_vgdcs_shell_xff", Kokkos::Experimental::require(        \
          Kokkos::RangePolicy<Kokkos::LaunchBounds<128, kLb>>(DevExeSpace(), 0, N),    \
          Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);
      const int esid_ = esid;
      esid = 0;
      if (xffs) {
        auto xsb_ = [&]() {
          auto dir_ = vgd_tabd.dir;
          auto mr_ = vgd_tabd.mr;
          auto vi_ = vgd_i_xb;
          auto msk_ = vgd_cscl[(xcurb >= 0) ? xcurb : vgd_cscur].msk;
#define VET_VGDCS_SHELL_NAME shellx_
#include "rad_vet_vetgd_cs_shell.inc"
#undef VET_VGDCS_SHELL_NAME
          return shellx_;
        }();
        const int nA = nlist[pass], nB = nlistb[pass];
        auto mainx_ = KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nA);
          int t = sb ? lstb_(lo + c - nA) : lst_(lo + c);
          const int d = t % n;
          t /= n;
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
        if (nA + nB > 0) {VET_VGDCS_LAUNCHN(mainx_, nA + nB)}
      } else {
      auto mainf_ = KOKKOS_LAMBDA(const int c) {
        int t = lst_(lo + c);
        const int d = t % n;
        t /= n;
        const int j = js + t % njj;
        t /= njj;
        const int k = ks + t % nkk;
        const int m = t/nkk;
        if (esid_ != 0 && est_(m, k - ks, j - js, d) == esid_) {return;}
        shell_(m, k, j, d, true);
      };
      Kokkos::parallel_for("vet_vgdcs_shell", Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::LaunchBounds<128, kLb>>(DevExeSpace(), 0, nper),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      }   // !xffs
      if (esid_ != 0) {eev.StreamWait(DevExeSpace());}   // the send cells
#undef VET_VGDCS_LAUNCHN
      if (rga) {
        // previous shell read: save its inward band values (inward pass), interior back
        // to the compact slab; the next shell gets a slab after this exchange
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
      if (early) {
#if MPI_PARALLEL_ENABLED
        // early send: pack done (event); messages of VetGdCsHaloC, then the unpack
        // (stream-ordered after the shell kernel)
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
          Kokkos::parallel_for("vet_vgdcs_cunpack",
                               VgdLight(DevExeSpace(), 0, static_cast<int>(nr)),
          KOKKOS_LAMBDA(const int o) {
            const int x = ri_(r0 + o);
            const int e = x/n, d = x - e*n;
            a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rbp[o];
          });
        }
#endif
      } else if (ccl) {
        // fused sweep: both sets in one round (VetGdCsHaloCX), else two exchanges
        if (!(xffs && VetGdCsHaloCX(pass, i, xcurb))) {
          VetGdCsHaloC(pass, i);          // read values only
          if (xffs) {
            // next set's exchange (own compact lists and array)
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
        // values to save from this inward band: its compact receive list (vgd_hxe[i])
        vgd_rg_sk[i] = 1;
        vgd_rg_se[i] = i;
      }
      if (rga) {VetGdRingFlush();}
      vgd_thalo += th.seconds();
    }
    if (rga) {
      // last shell of the pass: saved (inward), back to its compact slab
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

} // namespace radvet
