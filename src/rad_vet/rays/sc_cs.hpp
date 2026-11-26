#ifndef RAD_VET_RAYS_SC_CS_HPP_
#define RAD_VET_RAYS_SC_CS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cs.hpp
//! \brief cubed-sphere rays: the lateral stencil.

#include <math.h>

#include <string>

#include "athena.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"

namespace radvet {

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

static inline void VgdCsFatal(const std::string &msg) {
  VgdFatalAt(__FILE__, "<rad_vet>/vet_gd (cubed sphere): ", msg);
}

KOKKOS_INLINE_FUNCTION
static Real VgdCsDet3(const Real *a, const Real *b, const Real *c) {
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
static KOKKOS_INLINE_FUNCTION
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

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_CS_HPP_
