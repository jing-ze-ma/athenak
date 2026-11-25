//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetgd.cpp
//! \brief <rad_vet>/vet_gd: the Eddington tensor on the spherical-polar grid from a
//! short-characteristics formal solution along globally fixed (Cartesian) directions.
//!
//! The grey equation n . grad I = chi (S - I) is swept shell by shell for a fixed set of
//! straight-ray directions (HEALPix or product set; weights integrate 1, n, nn exactly).
//! Two passes, one launch per shell, one thread per (block, cell, direction):
//!   inward (n . r_hat < 0), top-down: upwind point on shell i+1 (vacuum at the top);
//!   outward (n . r_hat >= 0), bottom-up: shell i-1, or beyond the tangent point on the
//!     same shell (inward value), or the inner face (diffusion intensity).
//! Upwind I, ln chi, ln S are bilinear in (theta, phi) index space; S is linear in tau.
//! Periodic theta/phi walls: ghost intensities are re-indexed by the local-frame map
//! n -> M n (nearest set direction). Block boundaries: lagged inflow (block-Jacobi).
//! Output: the intensities vgd_i (CleanClosureFromGd, CleanGdFaceSums).

#include <algorithm>
#include <array>
#include <chrono>  // NOLINT(build/c++11)
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif

#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "utils/repro_sum.hpp"
#include "rad_vet/rad_vet_vgdlg.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_vetseg.hpp"
#include "rad_vet/rad_vet_vgdray.hpp"
#include "utils/recipe.hpp"
#include "coordinates/grid_stretch.hpp"
#if defined(KOKKOS_ENABLE_CUDA)
#include <cub/cub.cuh>
#endif

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {
// <rad_vet>/vet_gd_twin_srw: weight the twin's pattern subtraction per cell by the
// symref weight w = min(J_twin/J, 1) (VgdSrW); w = 1 cells are bitwise. Parsed in
// rad_vet.cpp.
#if !defined(VET_VGD_SWEEP_TU)
// <rad_vet>/vet_rotate_advance: with vet_rotate_avg = K the angle set of block b does
// not repeat with period K, so the averaged window advances through SO(3) / the azimuth
// (parsed in VetRotInit)
bool g_vra_adv = false;
#endif

namespace {
void VgdFatal(const std::string &msg) {VgdFatalAt(__FILE__, "<rad_vet>/vet_gd: ", msg);}

KOKKOS_INLINE_FUNCTION
Real VgdLerp(const Real a, const Real b, const Real u) {return a + u*(b - a);}

// sp stencil accessor of the shared upwind core (rad_vet_vgdray.hpp): the four
// (theta, phi) corners q = 0..3: (k0,j0), (k0,j1), (k1,j0), (k1,j1); nested bilinear
struct VgdSt4 {
  int k0, k1, j0, j1;
  Real uj, uk;
  KOKKOS_INLINE_FUNCTION int N() const {return 4;}
  KOKKOS_INLINE_FUNCTION int K(const int e) const {return (e < 2) ? k0 : k1;}
  KOKKOS_INLINE_FUNCTION int J(const int e) const {return (e % 2 == 0) ? j0 : j1;}
  KOKKOS_INLINE_FUNCTION Real W(const int e) const {
    return ((e < 2) ? (1.0 - uk) : uk)*((e % 2 == 0) ? (1.0 - uj) : uj);
  }
  template <typename F>
  KOKKOS_INLINE_FUNCTION Real Interp(const F &f) const {
    return VgdLerp(VgdLerp(f(k0,j0), f(k0,j1), uj), VgdLerp(f(k1,j0), f(k1,j1), uj), uk);
  }
};

// vet_gd_halo_ipc: device flags written over NVLink by the peers: [0, NS) arrivals by
// partner slot, [NS, 2 NS) acks; a spin gives up after ~20 s and raises the error flag
constexpr int VGD_IPC_NS = 32;
constexpr long long VGD_IPC_SPIN = 40000000000LL;
// HIP: wall_clock64() runs at 100 MHz on gfx942: 2e9 ticks = 20 s
constexpr long long VGD_IPC_SPIN_HIP = 2000000000LL;
struct VgdIpcIdx {int v[VGD_IPC_NS + 1];};
// HIP IPC halo: flags in fine-grained device memory, system-scope atomics; the receive
// buffer is written by the partner's pack kernel and read by my scatter kernel after
// k2 has seen the partner's flag (system fences + kernel boundaries order the data)
#if defined(KOKKOS_ENABLE_HIP)
#define VETG(x) hip##x
#define VETG_STREAM(ex) (ex).hip_stream()
#elif defined(KOKKOS_ENABLE_CUDA)
#define VETG(x) cuda##x
#define VETG_STREAM(ex) (ex).cuda_stream()
#endif

KOKKOS_INLINE_FUNCTION
void VgdIpcWaitGE(const long long *f, const long long target, int *err) {
#if defined(__CUDA_ARCH__)
  const volatile long long *vf = f;
  const long long t0 = clock64();
  while (*vf < target) {
    if (clock64() - t0 > VGD_IPC_SPIN) {*err = 1; return;}
  }
  __threadfence_system();
#elif defined(__HIP_DEVICE_COMPILE__)
  const long long t0 = wall_clock64();
  while (__hip_atomic_load(f, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_SYSTEM) < target) {
    if (wall_clock64() - t0 > VGD_IPC_SPIN_HIP) {*err = 1; return;}
    __builtin_amdgcn_s_sleep(1);
  }
  __threadfence_system();
#else
  (void)f; (void)target; (void)err;
#endif
}

KOKKOS_INLINE_FUNCTION
void VgdIpcSet(long long *f, const long long v) {
#if defined(__CUDA_ARCH__)
  __threadfence_system();
  volatile long long *vf = f;
  *vf = v;
  __threadfence_system();
#elif defined(__HIP_DEVICE_COMPILE__)
  __threadfence_system();
  __hip_atomic_store(f, v, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_SYSTEM);
  __threadfence_system();
#else
  (void)f; (void)v;
#endif
}

// lean IPC protocol: one system fence (orders the pack's peer stores before the
// flags), then the arrival flag v into every partner's slot; the wait: every own
// arrival flag >= v, then one fence
KOKKOS_INLINE_FUNCTION
void VgdIpcSetN(const unsigned long long *psig, const int *pslot, const int np,
                const long long v) {
#if defined(__CUDA_ARCH__)
  __threadfence_system();
  for (int p = 0; p < np; ++p) {
    volatile long long *vf = reinterpret_cast<long long *>(psig[p]) + pslot[p];
    *vf = v;
  }
#elif defined(__HIP_DEVICE_COMPILE__)
  __threadfence_system();
  for (int p = 0; p < np; ++p) {
    long long *vf = reinterpret_cast<long long *>(psig[p]) + pslot[p];
    __hip_atomic_store(vf, v, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_SYSTEM);
  }
#else
  (void)psig; (void)pslot; (void)np; (void)v;
#endif
}

KOKKOS_INLINE_FUNCTION
void VgdIpcWaitN(const long long *f, const int np, const long long target, int *err) {
#if defined(__CUDA_ARCH__)
  const long long t0 = clock64();
  for (int p = 0; p < np; ++p) {
    const volatile long long *vf = f + p;
    while (*vf < target) {
      if (clock64() - t0 > VGD_IPC_SPIN) {*err = 1; return;}
    }
  }
  __threadfence_system();
#elif defined(__HIP_DEVICE_COMPILE__)
  const long long t0 = wall_clock64();
  for (int p = 0; p < np; ++p) {
    while (__hip_atomic_load(f + p, __ATOMIC_RELAXED, __HIP_MEMORY_SCOPE_SYSTEM) <
           target) {
      if (wall_clock64() - t0 > VGD_IPC_SPIN_HIP) {*err = 1; return;}
      __builtin_amdgcn_s_sleep(1);
    }
  }
  __threadfence_system();
#else
  (void)f; (void)np; (void)target; (void)err;
#endif
}

// vet_interp = auto: material edge in the upwind stencil, i.e. the ln chi contrast over
// the corners of nonzero weight > log(kVetCornerRatio)
template <class CS>
KOKKOS_INLINE_FUNCTION
bool VgdCornerEdge(const CS &csx, const int m, const int k0, const int k1, const int j0,
                   const int j1, const int iu, const Real uk, const Real uj) {
  Real lmn = 1.0e300, lmx = -1.0e300;
  for (int q = 0; q < 4; ++q) {
    const Real wq = ((q < 2) ? (1.0 - uk) : uk)*((q % 2 == 0) ? (1.0 - uj) : uj);
    if (!(wq > 0.0)) {continue;}
    const Real l = csx(m,0,(q < 2) ? k0 : k1,(q % 2 == 0) ? j0 : j1,iu);
    lmn = fmin(lmn, l);
    lmx = fmax(lmx, l);
  }
  return (lmx - lmn) > log(kVetCornerRatio);
}

// HEALPix ring-scheme centres (12 ns^2) with the moment-fixed weights
void VgdDirections(const int ns, std::vector<double> &d) {
  d.clear();
  for (int i = 1; i < 4*ns; ++i) {
    double z, ph0, dph;
    int cnt;
    if (i < ns) {
      z = 1.0 - static_cast<double>(i*i)/(3.0*ns*ns);
      cnt = 4*i; dph = M_PI/(2.0*i); ph0 = -0.5*dph;
    } else if (i <= 3*ns) {
      z = 4.0/3.0 - 2.0*i/(3.0*ns);
      const int s = (i - ns + 1) % 2;
      cnt = 4*ns; dph = M_PI/(2.0*ns); ph0 = -0.5*s*dph;
    } else {
      const int ii = 4*ns - i;
      z = -(1.0 - static_cast<double>(ii*ii)/(3.0*ns*ns));
      cnt = 4*ii; dph = M_PI/(2.0*ii); ph0 = -0.5*dph;
    }
    const double st = std::sqrt(std::max(1.0 - z*z, 0.0));
    for (int j = 1; j <= cnt; ++j) {
      const double p = ph0 + j*dph;
      d.push_back(st*std::cos(p));
      d.push_back(st*std::sin(p));
      d.push_back(z);
      d.push_back(0.0);
    }
  }
  const int n = static_cast<int>(d.size()/4);
  double m2 = 0.0, m4 = 0.0;
  for (int q = 0; q < n; ++q) {
    const double z2 = d[4*q+2]*d[4*q+2];
    m2 += z2/n;
    m4 += z2*z2/n;
  }
  const double c = (1.0/3.0 - m2)/(m4 - m2*m2);
  for (int q = 0; q < n; ++q) {
    d[4*q+3] = (1.0 + c*(d[4*q+2]*d[4*q+2] - m2))/n;
  }
}

// product set: Gauss-Legendre in n_z (nmu nodes) x nphi uniform azimuths per ring,
// alternate rings shifted by half a step (stagger); positive weights, sum w = 1,
// sum w n = 0, sum w nn = I/3 exactly
void VgdProduct(const int nmu, const int nphi, const bool stagger,
                std::vector<double> &d) {
  d.clear();
  for (int q = 0; q < nmu; ++q) {       // Newton on P_nmu
    double z = std::cos(M_PI*(q + 0.75)/(nmu + 0.5));
    double pp = 1.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int l = 1; l <= nmu; ++l) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0*l - 1.0)*z*p2 - (l - 1.0)*p3)/l;
      }
      pp = nmu*(z*p1 - p2)/(z*z - 1.0);
      const double dz = p1/pp;
      z -= dz;
      if (std::fabs(dz) < 1.0e-15) {break;}
    }
    const double wq = 2.0/((1.0 - z*z)*pp*pp);      // [-1, 1] weight, sum 2
    const double st = std::sqrt(std::max(1.0 - z*z, 0.0));
    const double off = (stagger && (q % 2 == 1)) ? 0.5 : 0.0;
    for (int k = 0; k < nphi; ++k) {
      const double ph = (k + off)*2.0*M_PI/nphi;
      d.push_back(st*std::cos(ph));
      d.push_back(st*std::sin(ph));
      d.push_back(z);
      d.push_back(0.5*wq/nphi);
    }
  }
}

// every stencil box f(iu, k_lo, k_hi, j_lo, j_hi) (inclusive band cells of shell iu)
// that the shell kernel's solve of (m, k, j, d) at shell i may read. Robust to rounding
// of another compilation: both branches near the threshold (rel 1e-9), every floor
// within 1e-8. A superset of the reads: used by the overlap split and the halo lists.
// sth / pol / tn / c2w: stretched-theta node table and full-sphere chart (sph2.inc)
template <class TD, class TB, class TX, class TF, class TT, class FN>
KOKKOS_INLINE_FUNCTION
void VgdCand(const int m, const int k, const int j, const int d, const int i, const int l,
             const bool inw, const int n1, const int lcut, const int js, const int ks,
             const int jlo, const int jhi, const int klo, const int khi, const TD &dir_,
             const TB &mbsize, const VgdLG &lg_, const TX &cx1v, const TF &cx1f,
             const bool sth,
             const bool pol, const TT &tn, const int c2w, const bool imir, const FN &f) {
  if (inw && l == n1 - 1) {return;}          // the top shell reads no intensity
  const Real twopi = 2.0*M_PI;
  const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
  const Real th = sth ? tn(m,j)
                      : lg_.Th(mbsize, m, j, js);
  const Real ph = lg_.Ph(mbsize, m, k, ks);
  const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
  const Real mr = nx*st*cp + ny*st*sp + nz*ct;
  const Real r = cx1v(m,i);
  const Real zr = r*fabs(mr);
  const Real p2 = fmax(r*r - zr*zr, 0.0);
  auto box = [&](const int iu, const Real ds, const int dq) {
    const Real xu = r*st*cp - ds*nx;
    const Real yu = r*st*sp - ds*ny;
    const Real zu = r*ct - ds*nz;
    const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
    Real thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
    Real dph = atan2(yu, xu) - ph;
    dph -= twopi*floor((dph + M_PI)/twopi);
    VgdPoleChart(pol, th, thu, dph);
    const Real fj = sth ? (VgdThFj(tn, m, j, thu, c2w, lg_.Sj(m, js)) - lg_.Sj(m, js))
                        : VgdLG::Fs(j, lg_.Sj(m, js), (thu - th)/lg_.Dx2(mbsize, m));
    const Real fk = VgdLG::Fs(k, lg_.Sk(m, ks), dph/lg_.Dx3(mbsize, m));
    const Real e = 1.0e-8;
    int ja = static_cast<int>(floor(fj - e)), jb = static_cast<int>(floor(fj + e));
    int ka = static_cast<int>(floor(fk - e)), kb = static_cast<int>(floor(fk + e));
    ja = (ja < jlo) ? jlo : ((ja > jhi) ? jhi : ja);
    jb = (jb < jlo) ? jlo : ((jb > jhi) ? jhi : jb);
    ka = (ka < klo) ? klo : ((ka > khi) ? khi : ka);
    kb = (kb < klo) ? klo : ((kb > khi) ? khi : kb);
    f(iu, ka, kb + 1, ja, jb + 1, dq);
  };
  if (inw) {
    const Real ru = cx1v(m,i+1);
    box(i + 1, sqrt(fmax(ru*ru - p2, 0.0)) - zr, d);
    return;
  }
  const Real rd = (l > lcut) ? cx1v(m,i-1) : cx1f(m,i);
  const Real t2 = rd*rd, tol = 1.0e-9*t2;
  if (l > lcut) {
    if (p2 <= t2*(1.0 + 1.0e-13) + tol) {box(i - 1, zr - sqrt(fmax(t2 - p2, 0.0)), d);}
    if (p2 > t2*(1.0 + 1.0e-13) - tol) {box(i, 2.0*zr, d);}
  } else if (p2 > t2 - tol) {
    box(i, 2.0*zr, d);                      // a turned ray
  }
  // reflecting inner face (clean scheme, imir): the first shell's inward intensity at
  // the face point in the specular image direction (every direction within 1e-9 of it)
  if (imir && l == lcut && p2 <= t2*(1.0 + 1.0e-13) + tol) {
    const Real ds = zr - sqrt(fmax(t2 - p2, 0.0));
    const Real xu = r*st*cp - ds*nx, yu = r*st*sp - ds*ny, zu = r*ct - ds*nz;
    const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
    const Real ux = xu/ru, uy = yu/ru, uz = zu/ru;
    const Real mu = nx*ux + ny*uy + nz*uz;
    const int nd = static_cast<int>(dir_.extent(0));
    Real best = -2.0;
    for (int pass = 0; pass < 2; ++pass) {
      for (int q = 0; q < nd; ++q) {
        const Real c = dir_(q,0)*(nx - 2.0*mu*ux) + dir_(q,1)*(ny - 2.0*mu*uy) +
                       dir_(q,2)*(nz - 2.0*mu*uz);
        if (dir_(q,0)*ux + dir_(q,1)*uy + dir_(q,2)*uz >= 0.0) {continue;}
        if (pass == 0) {
          best = fmax(best, c);
        } else if (c >= best - 1.0e-9) {
          box(i, ds, q);
        }
      }
    }
  }
}

// two in-place exclusive int scans (flag -> position, total at index n) in one pass:
// the two counts ride in the low / high 32 bits of one int64 sum (counts < 2^31: exact).
// Chunked and deterministic: per VGD_SCAN_CH chunk a team sum, a scan of the chunk
// sums, a team scan per chunk from its offset.
constexpr int VGD_SCAN_CH = 4096;
void VgdScan2(const DevExeSpace &ex, const DvceArray1D<int> &a, const size_t na,
              const DvceArray1D<int> &b, const size_t nb, const char *name,
              Kokkos::View<int64_t*, DevMemSpace> &sums) {
  auto a_ = a;
  auto b_ = b;
  const size_t n = (na > nb) ? na : nb;
  auto val = KOKKOS_LAMBDA(const size_t q) -> int64_t {
    const int64_t va = (q < na) ? static_cast<int64_t>(a_(q)) : 0;
    const int64_t vb = (q < nb) ? static_cast<int64_t>(b_(q)) : 0;
    return va + (vb << 32);
  };
  auto put = KOKKOS_LAMBDA(const size_t q, const int64_t acc) {
    if (q <= na) {a_(q) = static_cast<int>(acc & 0xffffffffLL);}
    if (q <= nb) {b_(q) = static_cast<int>(acc >> 32);}
  };
  // entries [0, n] (the last one gets the total); chunk c = [c CH, (c+1) CH)
  const size_t ntot = n + 1;
  const int nch = static_cast<int>((ntot + VGD_SCAN_CH - 1)/VGD_SCAN_CH);
  if (static_cast<int>(sums.extent(0)) < nch + 1) {
    Kokkos::realloc(Kokkos::WithoutInitializing, sums, nch + 1);
  }
  auto s_ = sums;
  Kokkos::parallel_for(name, Kokkos::TeamPolicy<>(ex, nch, Kokkos::AUTO),
  KOKKOS_LAMBDA(const TeamMember_t &tm) {
    const size_t q0 = static_cast<size_t>(tm.league_rank())*VGD_SCAN_CH;
    const size_t q1 = (q0 + VGD_SCAN_CH < ntot) ? (q0 + VGD_SCAN_CH) : ntot;
    int64_t sm = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, static_cast<int>(q1 - q0)),
    [&](const int t, int64_t &acc) {acc += val(q0 + t);}, sm);
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {s_(tm.league_rank()) = sm;});
  });
  Kokkos::parallel_scan(name, Kokkos::RangePolicy<>(ex, 0, nch),
  KOKKOS_LAMBDA(const int c, int64_t &acc, const bool fin) {
    const int64_t v = s_(c);
    if (fin) {s_(c) = acc;}
    acc += v;
  });
  Kokkos::parallel_for(name, Kokkos::TeamPolicy<>(ex, nch, Kokkos::AUTO),
  KOKKOS_LAMBDA(const TeamMember_t &tm) {
    const size_t q0 = static_cast<size_t>(tm.league_rank())*VGD_SCAN_CH;
    const size_t q1 = (q0 + VGD_SCAN_CH < ntot) ? (q0 + VGD_SCAN_CH) : ntot;
    const int64_t off = s_(tm.league_rank());
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(tm, static_cast<int>(q1 - q0)),
    [&](const int t, int64_t &acc, const bool fin) {
      const int64_t v = val(q0 + t);
      if (fin) {put(q0 + t, off + acc);}
      acc += v;
    });
  });
}

// in-place exclusive int scan of a[0, n] (a[n] = 0 first, then the total) with CUB on
// the instance's stream (integer sums: the positions of VgdScan2). Non-CUDA: VgdScan2.
void VgdScanOne(const DevExeSpace &ex, const DvceArray1D<int> &a, const size_t n,
                Kokkos::View<int64_t*, DevMemSpace> &sums, DvceArray1D<char> &tmp) {
#if defined(KOKKOS_ENABLE_CUDA)
  int *p = a.data();
  cudaStream_t st = ex.cuda_stream();
  (void)cudaMemsetAsync(p + n, 0, sizeof(int), st);
  size_t need = 0;
  cub::DeviceScan::ExclusiveSum(nullptr, need, p, p, static_cast<int>(n + 1), st);
  if (tmp.extent(0) < need) {
    ex.fence();
    Kokkos::realloc(Kokkos::WithoutInitializing, tmp, need + need/2 + 1024);
  }
  size_t have = tmp.extent(0);
  cub::DeviceScan::ExclusiveSum(tmp.data(), have, p, p, static_cast<int>(n + 1), st);
#else
  VgdScan2(ex, a, n, a, 0, "vet_vgd_hxb_scan", sums);
  (void)tmp;
#endif
}

// vet_gd_geo = n >= 1 (stored as vgd_ls = -n): the geodesic set, the vertices of the
// class-I subdivided icosahedron (10 n^2 + 2 directions); weight = solid angle of the
// vertex's spherical Voronoi cell / 4 pi. Full icosahedral symmetry: sum w = 1,
// sum w n = 0, sum w nn = I/3, 4th moments isotropic to round-off (VgdMomentCheck)
namespace {
double VgdSphTri(const double *a, const double *b, const double *c) {
  // the solid angle of the spherical triangle of unit vectors a, b, c (Van Oosterom)
  const double tp = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) +
                    a[2]*(b[0]*c[1] - b[1]*c[0]);
  const double ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  const double bc = b[0]*c[0] + b[1]*c[1] + b[2]*c[2];
  const double ca = c[0]*a[0] + c[1]*a[1] + c[2]*a[2];
  return 2.0*std::fabs(std::atan2(tp, 1.0 + ab + bc + ca));
}
void VgdUnit(double *v) {
  const double s = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= s; v[1] /= s; v[2] /= s;
}
}  // namespace

bool VgdGeodesic(const int nl, std::vector<double> &d) {
  if (nl < 1 || nl > 12) {return false;}
  const double g = 0.5*(1.0 + std::sqrt(5.0));
  const double iv[12][3] = {{-1, g, 0}, {1, g, 0}, {-1, -g, 0}, {1, -g, 0},
                            {0, -1, g}, {0, 1, g}, {0, -1, -g}, {0, 1, -g},
                            {g, 0, -1}, {g, 0, 1}, {-g, 0, -1}, {-g, 0, 1}};
  const int fc[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
                         {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                         {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
                         {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  // the vertices: barycentric points (i, j, k)/n of every face, merged by position
  std::vector<std::array<double, 3>> vx;
  std::vector<std::array<int, 3>> tr;
  auto find = [&](double *p) {
    VgdUnit(p);
    for (size_t q = 0; q < vx.size(); ++q) {
      const double dx = vx[q][0] - p[0], dy = vx[q][1] - p[1], dz = vx[q][2] - p[2];
      if (dx*dx + dy*dy + dz*dz < 1.0e-20) {return static_cast<int>(q);}
    }
    vx.push_back({p[0], p[1], p[2]});
    return static_cast<int>(vx.size()) - 1;
  };
  for (int f = 0; f < 20; ++f) {
    const double *A = iv[fc[f][0]], *B = iv[fc[f][1]], *C = iv[fc[f][2]];
    std::vector<int> id((nl + 1)*(nl + 2)/2);
    auto at = [&](const int i, const int j) {return i*(nl + 1) - i*(i - 1)/2 + j;};
    for (int i = 0; i <= nl; ++i) {
      for (int j = 0; j <= nl - i; ++j) {
        const int k = nl - i - j;
        double p[3];
        for (int c = 0; c < 3; ++c) {p[c] = (i*A[c] + j*B[c] + k*C[c])/nl;}
        id[at(i, j)] = find(p);
      }
    }
    for (int i = 0; i < nl; ++i) {
      for (int j = 0; j < nl - i; ++j) {
        tr.push_back({id[at(i, j)], id[at(i + 1, j)], id[at(i, j + 1)]});
        if (j + 1 < nl - i) {
          tr.push_back({id[at(i + 1, j)], id[at(i + 1, j + 1)], id[at(i, j + 1)]});
        }
      }
    }
  }
  const int nv = static_cast<int>(vx.size());
  if (nv != 10*nl*nl + 2) {return false;}
  std::vector<double> w(nv, 0.0);
  for (const auto &t : tr) {
    double cc[3];
    const double *a = vx[t[0]].data(), *b = vx[t[1]].data(), *c = vx[t[2]].data();
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    cc[0] = u[1]*v[2] - u[2]*v[1];
    cc[1] = u[2]*v[0] - u[0]*v[2];
    cc[2] = u[0]*v[1] - u[1]*v[0];
    if (cc[0]*a[0] + cc[1]*a[1] + cc[2]*a[2] < 0.0) {
      cc[0] = -cc[0]; cc[1] = -cc[1]; cc[2] = -cc[2];
    }
    VgdUnit(cc);
    for (int e = 0; e < 3; ++e) {
      const double *p = vx[t[e]].data();
      const double *q1 = vx[t[(e + 1) % 3]].data(), *q2 = vx[t[(e + 2) % 3]].data();
      double m1[3] = {p[0] + q1[0], p[1] + q1[1], p[2] + q1[2]};
      double m2[3] = {p[0] + q2[0], p[1] + q2[1], p[2] + q2[2]};
      VgdUnit(m1);
      VgdUnit(m2);
      w[t[e]] += VgdSphTri(p, m1, cc) + VgdSphTri(p, cc, m2);
    }
  }
  double ws = 0.0;
  for (int q = 0; q < nv; ++q) {ws += w[q];}
  d.clear();
  for (int q = 0; q < nv; ++q) {
    d.push_back(vx[q][0]);
    d.push_back(vx[q][1]);
    d.push_back(vx[q][2]);
    d.push_back(w[q]/ws);
  }
  return true;
}

// moment errors of a direction set (rows x y z w): max |sum w - 1|, |sum w n|,
// |sum w n_i n_j - delta_ij/3|, |sum w n_i n_j n_k n_l - (dd + dd + dd)/15|
void VgdMomentCheck(const std::vector<double> &d, double e[4]) {
  const int n = static_cast<int>(d.size()/4);
  double m0 = 0.0, m1[3] = {0, 0, 0}, m2[3][3] = {{0}}, m4[3][3][3][3] = {{{{0}}}};
  for (int q = 0; q < n; ++q) {
    const double *x = &d[4*q], w = d[4*q+3];
    m0 += w;
    for (int i = 0; i < 3; ++i) {
      m1[i] += w*x[i];
      for (int j = 0; j < 3; ++j) {
        m2[i][j] += w*x[i]*x[j];
        for (int k = 0; k < 3; ++k) {
          for (int l = 0; l < 3; ++l) {m4[i][j][k][l] += w*x[i]*x[j]*x[k]*x[l];}
        }
      }
    }
  }
  e[0] = std::fabs(m0 - 1.0);
  e[1] = e[2] = e[3] = 0.0;
  for (int i = 0; i < 3; ++i) {
    e[1] = std::max(e[1], std::fabs(m1[i]));
    for (int j = 0; j < 3; ++j) {
      e[2] = std::max(e[2], std::fabs(m2[i][j] - ((i == j) ? 1.0/3.0 : 0.0)));
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          const double iso = ((i == j && k == l) + (i == k && j == l) +
                              (i == l && j == k))/15.0;
          e[3] = std::max(e[3], std::fabs(m4[i][j][k][l] - iso));
        }
      }
    }
  }
}

// vet_gd_ls: level-symmetric LQ_N (Carlson; mu_1 of Lewis & Miller), positive weights
// per permutation class from the even moments (sum w nn = I/3 by octahedral symmetry)
bool VgdLevelSym(const int nl, std::vector<double> &d) {
  double m1;
  switch (nl) {
    case 6: m1 = 0.2666355; break;
    case 8: m1 = 0.2182179; break;
    case 10: m1 = 0.1893213; break;
    case 12: m1 = 0.1672126; break;
    default: return false;
  }
  const int h = nl/2;
  const double dl = 2.0*(1.0 - 3.0*m1*m1)/(nl - 2);
  std::vector<double> mu(h);
  for (int q = 0; q < h; ++q) {mu[q] = std::sqrt(m1*m1 + dl*q);}
  std::vector<std::vector<int>> cls;     // sorted index triples
  std::vector<int> pc;                   // class of each (unsigned) point
  std::vector<double> pt;
  for (int a = 0; a < h; ++a) {
    for (int b = 0; b < h; ++b) {
      const int c = h - 1 - a - b;
      if (c < 0) {continue;}
      std::vector<int> t = {a, b, c};
      std::sort(t.begin(), t.end());
      int id = -1;
      for (int q = 0; q < static_cast<int>(cls.size()); ++q) {
        if (cls[q] == t) {id = q;}
      }
      if (id < 0) {cls.push_back(t); id = static_cast<int>(cls.size()) - 1;}
      pc.push_back(id);
      pt.push_back(mu[a]); pt.push_back(mu[b]); pt.push_back(mu[c]);
    }
  }
  const int nc = static_cast<int>(cls.size());
  const int np = static_cast<int>(pc.size());
  // A w = rhs over one octant (weights per point sum to 1/8 there)
  std::vector<double> A(nc*nc, 0.0), r(nc);
  for (int row = 0; row < nc; ++row) {
    const int p2 = (row == 0) ? 0 : (2*(row + 1));
    for (int q = 0; q < np; ++q) {A[row*nc + pc[q]] += std::pow(pt[3*q], p2)*8.0;}
    r[row] = 1.0/(p2 + 1.0);
  }
  for (int c = 0; c < nc; ++c) {       // Gauss elimination with partial pivoting
    int pv = c;
    for (int q = c + 1; q < nc; ++q) {
      if (std::fabs(A[q*nc + c]) > std::fabs(A[pv*nc + c])) {pv = q;}
    }
    for (int q = 0; q < nc; ++q) {std::swap(A[c*nc + q], A[pv*nc + q]);}
    std::swap(r[c], r[pv]);
    for (int q = c + 1; q < nc; ++q) {
      const double f = A[q*nc + c]/A[c*nc + c];
      for (int t = c; t < nc; ++t) {A[q*nc + t] -= f*A[c*nc + t];}
      r[q] -= f*r[c];
    }
  }
  std::vector<double> w(nc);
  for (int c = nc - 1; c >= 0; --c) {
    double v = r[c];
    for (int t = c + 1; t < nc; ++t) {v -= A[c*nc + t]*w[t];}
    w[c] = v/A[c*nc + c];
  }
  d.clear();
  for (int q = 0; q < np; ++q) {
    for (int sg = 0; sg < 8; ++sg) {
      d.push_back(((sg & 1) ? -1.0 : 1.0)*pt[3*q]);
      d.push_back(((sg & 2) ? -1.0 : 1.0)*pt[3*q+1]);
      d.push_back(((sg & 4) ? -1.0 : 1.0)*pt[3*q+2]);
      d.push_back(w[pc[q]]);
    }
  }
  return true;
}
} // namespace
#if !defined(VET_VGD_SWEEP_TU)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdBuildSet
//! \brief the vet_gd direction set of the keys: product set, vet_gd_ls, vet_gd_geo or
//! HEALPix (vet_gd_nside); shared by vet_gd and vet_sc_set = gd

void RadiationVet::VgdBuildSet(std::vector<double> &d) const {
  if (vgd_gl_nmu > 0) {
    VgdProduct(vgd_gl_nmu, vgd_gl_nphi, vgd_gl_stag, d);
  } else if (vgd_ls > 0) {
    if (!VgdLevelSym(vgd_ls, d)) {VgdFatal("vet_gd_ls must be 6, 8, 10 or 12");}
  } else if (vgd_ls < 0) {   // vet_gd_geo
    if (!VgdGeodesic(-vgd_ls, d)) {VgdFatal("vet_gd_geo must be 1 .. 12");}
  } else {
    VgdDirections(vgd_nside, d);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdInit
//! \brief direction table, intensity array, exchange object, wall maps (from VetRayInit)

void RadiationVet::VetGdInit() {
  Mesh *pm = pmy_pack->pmesh;
  // nx3 = 1: axisymmetric; the phi band columns are the cell's own column with the
  // directions rotated by their phi offset (must map the set onto itself)
  if (!pm->multi_d) {VgdFatal("needs a 3-D wedge or a 2-D axisymmetric mesh");}
  if (!pm->three_d && (impl_beam_hr || impl_beam_fs)) {
    VgdFatal("2-D axisymmetric: no implicit_flux_beam");
  }
  if (ibc_x1max != VET_IBC_MARSHAK &&
      !(ibc_x1max == VET_IBC_REFLECT || ibc_x1max == VET_IBC_OUTFLOW)) {
    // the sweep has a vacuum top, or a mirror (VetGdTopMirror) for the clean scheme
    VgdFatal("outer x1 needs implicit_bc_x1max = marshak or "
             "reflect with scheme = clean");
  }
  if (ibc_x2min == VET_IBC_MARSHAK || ibc_x2max == VET_IBC_MARSHAK ||
      ibc_x3min == VET_IBC_MARSHAK || ibc_x3max == VET_IBC_MARSHAK) {
    // the gd sweep has no vacuum transverse face (vet_sc has)
    VgdFatal("implicit_bc_x2*/x3* = marshak is not wired into the vet_gd sweep");
  }
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  // cubed sphere (rad_vet_vetgd_cs.cpp): wedge-only options are refused; the cs band
  // is filled by its own exact copies (compact lists, see VetGdCsSweep); no walls
  vgd_csg = pm->use_cubed_sphere;
  // <mesh>/layout_invariant: the sweep's (theta, phi) from the global index
  // (rad_vet_vgdlg.hpp); a uniform single-level sp mesh
  vgd_lg.on = false;
  if (!vgd_csg && pm->layout_invariant) {
    vgd_lg.on = true;
    vgd_lg.x2a = pm->mesh_size.x2min;
    vgd_lg.dx2 = (pm->mesh_size.x2max - pm->mesh_size.x2min)/pm->mesh_indcs.nx2;
    vgd_lg.x3a = pm->mesh_size.x3min;
    vgd_lg.dx3 = (pm->mesh_size.x3max - pm->mesh_size.x3min)/pm->mesh_indcs.nx3;
    Kokkos::realloc(vgd_lg.o, pmy_pack->nmb_thispack, 2);
    vgd_lg.oh = Kokkos::create_mirror_view(vgd_lg.o);
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      const auto &ll = pm->lloc_eachmb[pmy_pack->gids + m];
      if (ll.level != pm->root_level) {
        VgdFatal("layout_invariant: the global vet_gd geometry needs a single-level "
                 "mesh");
      }
      vgd_lg.oh(m,0) = static_cast<int>(ll.lx2)*pm->mb_indcs.nx2;
      vgd_lg.oh(m,1) = static_cast<int>(ll.lx3)*pm->mb_indcs.nx3;
    }
    Kokkos::deep_copy(vgd_lg.o, vgd_lg.oh);
  }
  if (vgd_csg) {
    std::string bad;
    if (!bad.empty()) {VgdFatal("not available on the cubed sphere:" + bad);}
    vgd_wint = false;
    vgd_ipc_on = false;
    vgd_hx_on = false;
    // the unfused twin always uses this build's cut (as the fused twin); the previous
    // build's cut would make it decomposition-dependent
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd on the CUBED SPHERE: seam "
                << "strip + vertex triangle upwind stencil, exact per-shell band copies "
                << "(the sp halo lists / IPC / overlap keys do not apply)"
                << std::endl;
    }
  }
  // stretched theta / full sphere (rad_vet_vetgd_sph2.inc): options that assume the
  // uniform index formula are refused; the compact halo falls back to the dense one
  // (same values)
  if (vgd_sth) {
    std::string bad;
    if (vgd_thmir || vgd_thvac > 0) {bad += " vet_gd_theta_wall = mirror|vacuum|extrap;";}
    if (!bad.empty()) {
      VgdFatal("not available on a stretched theta / the full sphere:" +
               bad);
    }
    // with vet_gd_band_global the compact halo is the keep-tested global-band exchange
    // (VetGdGbHaloC); otherwise the dense halo (the slot keep test assumes uniform theta)
    if (vgd_hcomp > 0 && !vgd_gband) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd: vet_gd_halo_compact = " << vgd_hcomp
                  << " -> 0 (the dense per-shell slot halo; the compact keep test is "
                  << "ported only to vet_gd_band_global; speed only)" << std::endl;
      }
      vgd_hcomp = 0;
    }
    if (vgd_hcomp > 0) {
      // the slot-path speed machinery does not apply to the global band
      vgd_hx_on = false;
    }
    vgd_ipc_on = false;
  }
  // several MeshBlocks along x1: levels swept in turn, the boundary shell handed to the
  // next (VetGdX1Hand); flags 1 innermost, 2 outermost
  vgd_nb1 = 1;
  for (int g = 0; g < pm->nmb_total; ++g) {
    vgd_nb1 = std::max(vgd_nb1, pm->lloc_eachmb[g].lx1 + 1);
  }
  vgd_x1nb.assign(2*nmb, -1);
  vgd_x1lv.assign(nmb, 0);
  Kokkos::realloc(vgd_x1f, nmb);
  auto xf_h = Kokkos::create_mirror_view(vgd_x1f);
  for (int m = 0; m < nmb; ++m) {
    const LogicalLocation &l = pm->lloc_eachmb[pmy_pack->gids + m];
    vgd_x1lv[m] = l.lx1;
    xf_h(m) = (l.lx1 == 0 ? 1 : 0) + (l.lx1 == vgd_nb1 - 1 ? 2 : 0);
    if (l.lx1 == vgd_nb1 - 1 && ibc_x1max == VET_IBC_REFLECT) {xf_h(m) += 4;}
    xf_h(m) += 8*l.lx1;
    for (int g = 0; g < pm->nmb_total; ++g) {
      const LogicalLocation &o = pm->lloc_eachmb[g];
      if (o.level == l.level && o.panel == l.panel && o.lx2 == l.lx2 && o.lx3 == l.lx3 &&
          std::abs(o.lx1 - l.lx1) == 1) {
        vgd_x1nb[2*m + (o.lx1 > l.lx1 ? 1 : 0)] = g;
      }
    }
  }
  Kokkos::deep_copy(vgd_x1f, xf_h);
  Kokkos::realloc(vgd_x1mo, nmb);
  auto mo_h = Kokkos::create_mirror_view(vgd_x1mo);
  for (int m = 0; m < nmb; ++m) {mo_h(m) = m;}
  std::stable_sort(mo_h.data(), mo_h.data() + nmb,
                   [&](const int a, const int b) {return vgd_x1lv[a] < vgd_x1lv[b];});
  Kokkos::deep_copy(vgd_x1mo, mo_h);
  if (vgd_nb1 > 1) {
    if (vgd_csg || vgd_pol || indcs.nx1 < 2) {
      VgdFatal("several MeshBlocks along x1 need a wedge (no full sphere) and "
               "MeshBlock nx1 >= 2");
    }
    vgd_hcomp = 0;   // the dense per-shell halo (the compact lists assume one level)
  }
  // vgd_raw: the raw band where the exact-list band would clamp (same values where
  // nothing clamps)
  const bool raw_ok = !vgd_csg && !vgd_pol && !vgd_sth;
  vgd_raw = false;
  // nx3 = 1: the dense halo (the compact / exact lists mark 3-D band reads)
  if (!pm->three_d) {vgd_hcomp = 0;}
  auto raw_halo = [&]() {
    if (vgd_raw || !pm->three_d) {   // the compact halo: VetGdGbHaloC
      vgd_hx_on = false;
      vgd_ipc_on = vgd_ipc_on && vgd_raw && pm->three_d;   // VetGdGbHaloC's own IPC
    }
  };
  raw_halo();
  // vet_gd_wall_sym / vet_gd_theta_seam_sym do nothing without a phi wall or a
  // periodic theta seam
  const bool sp3 = !vgd_csg && !vgd_pol && pm->three_d;
  vgd_wsym = vgd_wsym && sp3 &&
             pm->mesh_size.x3max - pm->mesh_size.x3min < 2.0*M_PI - 1.0e-12;
  vgd_wst = vgd_wst && sp3 &&
            pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic;
  // the direction set: any table of unit vectors + weights with sum w = 1,
  // sum w nn = I/3
  if (vgd_gl_nmu > 0) {
    // product set: the azimuth step must divide the phi period (exact walls);
    // nphi >= 4 on the cubed sphere
    const double per = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const double m = per*vgd_gl_nphi/(2.0*M_PI);
    if (vgd_gl_nphi < 4 ||
        (!vgd_csg && !vgd_wsym && std::fabs(m - std::round(m)) > 1.0e-9)) {
      VgdFatal("vet_gd_gl_nphi must be >= 4 with 2 pi/nphi dividing the phi period of "
               "the wedge (exact wall maps), or vet_gd_wall_sym = true");
    }
  }
  VgdBuildSet(vgd_base);
  // the set's image under (x, y, z) -> (c x - t s y, s x + t c y, tz z): a turn about z
  // (t = 1), a meridian mirror (t = -1) or the equatorial mirror (tz = -1); -1 where
  // the image is not in the set
  const size_t nb = vgd_base.size()/4;
  auto img = [&](const double c, const double s, const double t, const double tz) {
    std::vector<int> p(nb, -1);
    for (size_t q = 0; q < nb; ++q) {
      const double *u = &vgd_base[4*q];
      const double mx = c*u[0] - t*s*u[1], my = s*u[0] + t*c*u[1];
      for (size_t r = 0; r < nb && p[q] < 0; ++r) {
        const double *v = &vgd_base[4*r];
        if (mx*v[0] + my*v[1] + tz*u[2]*v[2] > 1.0 - 1.0e-12 &&
            std::fabs(v[3] - u[3]) <= 1.0e-12*u[3]) {p[q] = static_cast<int>(r);}
      }
    }
    return p;
  };
  auto dev = [&](const std::vector<int> &p, const char *nm) {
    if (*std::min_element(p.begin(), p.end()) < 0) {
      VgdFatal(std::string(nm) + ": the set is not mirror-symmetric about the wedge "
               "mid-plane");
    }
    DvceArray1D<int> v(std::string(nm), nb);
    auto h = Kokkos::create_mirror_view(v);
    for (size_t q = 0; q < nb; ++q) {h(q) = p[q];}
    Kokkos::deep_copy(v, h);
    return v;
  };
  if (vgd_wsym && VgdRotOn()) {
    // a set turning about z on a wedge of whole azimuthal periods: the turn by the period
    // W maps every turned set onto itself (the two turns commute), so its phi wall maps
    // are exact and there is nothing to symmetrise (the key does nothing)
    const double w = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const std::vector<int> pw = img(std::cos(w), std::sin(w), 1.0, 1.0);
    if (!vgd_wint || vgd_rot3 || *std::min_element(pw.begin(), pw.end()) < 0) {
      VgdFatal("vet_gd_wall_sym with a rotating set needs vet_gd_wall_interp = true, a "
               "turn about z and a wedge of whole azimuthal periods of the set");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_wall_sym: the turning set's phi wall maps are exact "
                << "(no effect)" << std::endl;
    }
    vgd_wsym = false;
  }
  if (vgd_wsym) {
    // vet_gd_wall_sym: turn the set about z by the least angle that puts one of its
    // mirror planes (GL k pi/nphi, the octahedral sets k pi/4) on the wedge mid-plane
    if (!vgd_wint) {
      VgdFatal("vet_gd_wall_sym needs vet_gd_wall_interp = true");
    }
    const double pc = 0.5*(pm->mesh_size.x3min + pm->mesh_size.x3max);
    const double pp = (vgd_gl_nmu > 0) ? M_PI/vgd_gl_nphi : 0.25*M_PI;
    const double da = pc - pp*std::round(pc/pp);
    const double ca = std::cos(da), sa = std::sin(da);
    if (std::fabs(da) > 1.0e-12) {
      for (size_t q = 0; q < nb; ++q) {
        const double x = vgd_base[4*q], y = vgd_base[4*q+1];
        vgd_base[4*q] = ca*x - sa*y;
        vgd_base[4*q+1] = sa*x + ca*y;
      }
    }
    vgd_mir = dev(img(std::cos(2.0*pc), std::sin(2.0*pc), -1.0, 1.0), "vet_gd_wall_sym");
    // inexact phi wall map (the turn by the period is not a set map): the clean
    // scheme's seam faces in both frames (CleanGdFaceSums)
    vgd_wsw = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const std::vector<int> pw = img(std::cos(vgd_wsw), std::sin(vgd_wsw), 1.0, 1.0);
    vgd_wsx = *std::min_element(pw.begin(), pw.end()) < 0;
  }
  if (vgd_wst) {
    // vet_gd_theta_seam_sym: the clean scheme's theta seam faces in both frames, with the
    // equatorial mirror (a theta wedge centred on the equator).  A turn about z keeps the
    // set's equatorial mirror (the permutation of every turned set); a 3-D one does not
    if (!vgd_wint || (VgdRotOn() && vgd_rot3) ||
        std::fabs(pm->mesh_size.x2min + pm->mesh_size.x2max - M_PI) > 1.0e-12) {
      VgdFatal("vet_gd_theta_seam_sym needs a theta wedge centred on the equator, "
               "vet_gd_wall_interp = true and no 3-D set rotation");
    }
    vgd_mirz = dev(img(1.0, 0.0, 1.0, -1.0), "vet_gd_theta_seam_sym");
  }
  // the moments of the set, as built and after a fixed random 3-D rotation
  if (global_variable::my_rank == 0) {
    double e[4], er[4];
    VgdMomentCheck(vgd_base, e);
    const double ax[3] = {0.3713906763541037, -0.5570860145311556, 0.7427813527082074};
    const double an = 0.987654321, ca = std::cos(an), sa = std::sin(an);
    std::vector<double> rb(vgd_base);
    for (size_t q = 0; q < rb.size()/4; ++q) {   // Rodrigues
      const double *x = &vgd_base[4*q];
      const double kx = ax[1]*x[2] - ax[2]*x[1], ky = ax[2]*x[0] - ax[0]*x[2],
                   kz = ax[0]*x[1] - ax[1]*x[0];
      const double kd = ax[0]*x[0] + ax[1]*x[1] + ax[2]*x[2];
      rb[4*q] = x[0]*ca + kx*sa + ax[0]*kd*(1.0 - ca);
      rb[4*q+1] = x[1]*ca + ky*sa + ax[1]*kd*(1.0 - ca);
      rb[4*q+2] = x[2]*ca + kz*sa + ax[2]*kd*(1.0 - ca);
    }
    VgdMomentCheck(rb, er);
    std::cout << "<rad_vet> vet_gd direction set: " << vgd_base.size()/4
              << " directions; moment errors |sum w - 1| " << e[0] << ", |sum w n| "
              << e[1]
              << ", |sum w nn - I/3| " << e[2] << " (rotated " << er[2]
              << "), 4th-moment anisotropy " << e[3] << " (rotated " << er[3] << ")"
              << std::endl;
  }
  vgd_n = static_cast<int>(vgd_base.size()/4);
  // the set's z-rotation period: product 2 pi/nphi, geodesic pi, HEALPix / LQ_N pi/2
  vgd_zper = (vgd_gl_nmu > 0) ? 2.0*M_PI/vgd_gl_nphi : ((vgd_ls < 0) ? M_PI : 0.5*M_PI);
  if (!pm->three_d) {
    const double per = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const double c = std::cos(per), s = std::sin(per);
    for (int q = 0; q < vgd_n; ++q) {
      const double *u = &vgd_base[4*q];
      double bd = -2.0;
      for (int r = 0; r < vgd_n; ++r) {
        const double *v = &vgd_base[4*r];
        bd = std::max(bd, (c*u[0] - s*u[1])*v[0] + (s*u[0] + c*u[1])*v[1] + u[2]*v[2]);
      }
      if (bd < 1.0 - 1.0e-12) {
        VgdFatal("2-D axisymmetric: the phi extent of the mesh must rotate the direction "
                 "set onto itself (e.g. 2 pi/vet_gd_gl_nphi, pi/2 for HEALPix)");
      }
    }
  }
  if (vgd_n > VET_VGD_NMAX) {VgdFatal("too many directions");}
  const int n = vgd_n;
  Kokkos::realloc(vgd_dir, n, 4);
  Kokkos::realloc(vgd_x1b, nmb, std::max(n, 4)*indcs.nx3*indcs.nx2);
  Kokkos::realloc(vgd_x1r, nmb, std::max(n, 4)*indcs.nx3*indcs.nx2);
  // lateral band: the deepest reach of an upwind point in cells over every shell (cell
  // sizes r dtheta, r sin(theta) dphi); on the cubed sphere VetGdCsInit
  if (vgd_csg) {
    VetGdCsInit();
  } else {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto &mbs = pmy_pack->pmb->mb_size;
    const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
    const double smin = std::min(std::sin(t0), std::sin(t1));
    // stretched theta / full sphere: the narrowest row and the band latitude limit
    // (VetGdThetaDmin)
    const double dmin = vgd_sth ? VetGdThetaDmin()
                                : std::min(vgd_lg.Dx2H(mbs, 0),
                                           smin*vgd_lg.Dx3H(mbs, 0));
    // shell s is read (a) by s-1 inward, chord <= sqrt(r_s^2 - r_{s-1}^2), (b) by s+1
    // outward, chord <= sqrt(r_{s+1}^2 - r_s^2), (c) by s itself (turn / inner face),
    // chord <= 2 sqrt(r_s^2 - r_lo^2); angular offset <= chord / r_lo; + 2 cells
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    // pass 1: the raw band, where pass 0's exact lists would clamp
    for (int pass = 0; pass < 2; ++pass) {
      int nmax = vgd_raw ? std::numeric_limits<int>::max()
                         : (pm->three_d ? std::min(indcs.nx2, indcs.nx3) : indcs.nx2);
      // vet_gd_band_global: the band is filled from any block, so its depth is capped
      // by the sphere (half a ring in phi, pole-to-pole in theta): layout independent
      if (vgd_gband) {
        if (!vgd_pol) {
          VgdFatal("vet_gd_band_global needs the full sphere (mesh/use_polar_boundary)");
        }
        nmax = std::min(pm->mesh_indcs.nx3/2, pm->mesh_indcs.nx2);
      }
      vgd_wsh.assign(indcs.nx1 + 2*indcs.ng, 0);
      vgd_wsi.assign(indcs.nx1 + 2*indcs.ng, 0);
      vgd_wso.assign(indcs.nx1 + 2*indcs.ng, 0);
      int wneed = 0;
      for (int m = 0; m < ((vgd_nb1 > 1) ? nmb : 1); ++m) {
        // shells next to another level: radii from the x1 ghosts
        const bool lo = (xf_h(m) & 1), hi = (xf_h(m) & 2);
        for (int l = indcs.is; l <= indcs.ie; ++l) {
          const double r = x1v_h(m,l);
          const double rlo = (l > indcs.is || !lo) ? x1v_h(m,l-1) : x1f_h(m,l);
          const bool up = (l < indcs.ie || !hi);
          double a = 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
          if (up) {
            a = std::max(a, std::sqrt(x1v_h(m,l+1)*x1v_h(m,l+1) - r*r));
          }
          const int wl = static_cast<int>(std::ceil(a/(rlo*dmin))) + 2;
          wneed = std::max(wneed, wl);
          vgd_wsh[l] = std::min(wl, nmax);
          // inward-pass halo of shell s serves (a) and (c), the outward-pass halo (b)
          const double ain = 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
          const double aout = up ? std::sqrt(x1v_h(m,l+1)*x1v_h(m,l+1) - r*r) : 0.0;
          vgd_wsi[l] = std::min(static_cast<int>(std::ceil(ain/(rlo*dmin))) + 2, nmax);
          vgd_wso[l] = std::min(static_cast<int>(std::ceil(aout/(rlo*dmin))) + 2, nmax);
        }
      }
  #if MPI_PARALLEL_ENABLED
      if (vgd_nb1 > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &wneed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
  #endif
      vgd_w = std::min(wneed, nmax);
      if (pass == 0 && raw_ok && wneed > nmax) {
        vgd_raw = true;
        raw_halo();
        if (global_variable::my_rank == 0) {
          std::cout << "<rad_vet> vet_gd: the exact-list band needs "
                    << wneed << " cells, MeshBlocks of " << indcs.nx2 << " x "
                    << indcs.nx3 << " (theta x phi) hold " << nmax << ": the raw "
                    << "global band (no clamped reads)" << std::endl;
        }
        continue;
      }
      if (vgd_nb1 > 1) {   // every shell at the full band, the x1 ghosts too
        vgd_wsh.assign(c1, vgd_w);
        vgd_wsi = vgd_wsh;
        vgd_wso = vgd_wsh;
      }
      vgd_capped = (vgd_w < wneed);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd: exact per-shell lateral halo, band " << vgd_w
                  << " cells (needed " << wneed << "), per-shell depth mean "
                  << std::accumulate(vgd_wsh.begin(), vgd_wsh.end(), 0.0)/indcs.nx1
                  << (vgd_capped ? ": MeshBlocks too small, the deepest near-tangent "
                                   "reads are clamped (counted)" : "") << std::endl;
      }
      break;
    }
  }
  const int c2w = indcs.nx2 + 2*vgd_w, c3w = indcs.nx3 + 2*vgd_w;
  if (vgd_sth) {VetGdThetaInit();}   // rad_vet_vetgd_sph2.inc
  Kokkos::realloc(vgd_cs, nmb, 2, c3w, c2w, c1);
  Kokkos::realloc(vgd_wall, nmb, c3w, c2w);
  Kokkos::realloc(vgd_map, nmb, c3w, c2w, n);
  Kokkos::realloc(vgd_mr, nmb, c3w, c2w, n);
  if (vgd_wint) {
    Kokkos::realloc(vgd_m3, nmb, c3w, c2w, n, 3);
    Kokkos::realloc(vgd_w3, nmb, c3w, c2w, n, 3);
  }
  if (!vgd_csg) {VetGdHaloInit();}   // cs lists made by VetGdCsInit
  if (vgd_gband || vgd_raw) {VetGdGlobalBandInit();}   // sph2.inc: after the slots
  if (vgd_thmir) {
    // vet_gd_theta_wall = mirror: the theta neighbour across a wall is the block itself
    // (one MeshBlock in theta): corner slots are phi-edge slots
    if (vgd_csg) {
      VgdFatal("vet_gd_theta_wall = mirror: no theta walls on the cubed sphere");
    }
    if (vgd_nb2 != 1) {
      VgdFatal("vet_gd_theta_wall = mirror needs ONE MeshBlock along x2 (theta)");
    }
    if (vgd_ipc_on) {VgdFatal("vet_gd_theta_wall = mirror: vet_gd_halo_ipc not ported");}
  }
  if (vgd_thvac > 0) {
    // vet_gd_theta_wall = vacuum / extrap: only the theta ghosts beyond a physical wall
    // change (VgdTabCompute flags 3 / 4); vacuum allows any MeshBlock count in theta
    if (vgd_csg) {VgdFatal("vet_gd_theta_wall = vacuum/extrap: not on the cubed sphere");}
    if (vgd_ipc_on) {
      VgdFatal("vet_gd_theta_wall = vacuum/extrap: vet_gd_halo_ipc not ported");
    }
    if (vgd_thvac == 2 && vgd_nb2 != 1) {
      VgdFatal("vet_gd_theta_wall = extrap needs ONE MeshBlock along x2 (theta)");
    }
  }
  {
    // vet_gd_band_ring: only where every band value the sweep reads comes from this
    // sweep's exact-list (cs: compact; vgd_raw: own) exchange and some block has a
    // remote side; without the cs compact lists the clean sc rays keep the ragged band
    if (vgd_csg && vgd_hcomp == 0) {
      VgdFatal("vet_gd_halo_compact = 0 on the cubed sphere with the clean sc rays (no "
               "compact lists, no band ring: the full band); name vet_gd_halo_compact "
               "= 2");
    }
    const bool hxok = vgd_raw ? (vgd_hmpi && vgd_nb1 == 1)
                    : vgd_csg ? (vgd_cscomp && vgd_hmpi)
                              : (vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0));
    bool anyside = false;
    auto hl_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_hloc);
    for (size_t q = 0; q < hl_h.extent(0); ++q) {anyside = anyside || (hl_h(q) < 0);}
    vgd_ring = vgd_ring_par && hxok && anyside;
    // the fused two-set sweep keeps the ragged band on sp; on the cubed sphere the ring
    // records both sets' compact lists (VetGdCsSweep's ringrec)
  }
  // combinations the fused sweep never covers (VetGdXffOk) keep two sweeps and no
  // second array (the clean ring, CleanRingMoments, sweeps one set at a time)
  if (vgd_xff && ((cln_ring && !(cln_rfuse && cln_ring == 1 &&
                                          ibc_x1max != VET_IBC_REFLECT)) ||
                  (vgd_gband && !vgd_raw) || vgd_nb1 > 1)) {
    vgd_xff = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: not covered here ("
                << "xfade_m, ...): two sweeps per build" << std::endl;
    }
  }
  // the shells below the ray cut are filled just before a set's closure reads them
  // (CleanClosureFromGd): the fused sweep's two sets share their slabs.  From the cut of
  // the restart file, with a margin; a lower cut takes two sweeps (VetGdXffOk)
  vgd_rag_cut = (vgd_xff && vgd_nb1 == 1 && marena_on && !marena_sealed) ?
                std::max(cln_rcut - indcs.nx1/32, 0) : 0;
  VgdRagAlloc(vgd_i);   // after the halo init: band sides follow vgd_hloc
  if (vgd_xff) {
    // the other set's intensities, same layout and band sides
    VgdRagAlloc(vgd_i_xb);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: second-set intensities "
                << static_cast<double>(vgd_i_xb.d.extent(0))*sizeof(Real)/1.0e9
                << " GB on rank 0" << std::endl;
    }
  }
  VgdRagThick();
  // vet_gd_halo_exact marks the reads in the shell-list sweep
  if (vgd_hx_on && (!vgd_hmpi || vgd_hcomp <= 0)) {
    vgd_hx_on = false;         // does not apply (or nothing to exchange)
  }
  vgd_alpha = -1.0;
  // table prefetch (own device instance for the upload)
  if (VgdRotOn()) {
    vgd_tex = Kokkos::Experimental::partition_space(DevExeSpace(), 1)[0];
  }
  VetGdTables(VetGdAngle(pm->ncycle));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotInit
//! \brief rotation keys of the angle set, shared by vet_gd and vet_sc (one clock, one
//! key set); called from the closure parse, curved and Cartesian

void RadiationVet::VetRotInit(ParameterInput *pin) {
  const bool fresh = !global_variable::restart_run;
  // the rotated sets do not repeat (the averaged last-K window advances) on the curved
  // meshes and on the ring
  g_vra_adv = vgd_on || cln_ring > 0;
  {
    const bool rdt_named = pin->DoesParameterExist("rad_vet","vet_gd_rotate_dt");
    const bool rev_named = pin->DoesParameterExist("rad_vet","vet_gd_rotate_every");
    // vet_sc rotates on the same clock with the same keys, unless its set is pinned
    // by boundary beams (implicit_beam_n)
    const bool sc_pin = vet_sc &&
        (pin->DoesParameterExist("rad_vet","implicit_beam_n") &&
         pin->GetInteger("rad_vet","implicit_beam_n") > 0);
    // a pinned vet_sc set gets no default rotation and a named one is fatal (the beams
    // are injected along fixed rays of the set)
    const bool rot_cl = vgd_on || vet_sc;
    bool rot_phys = false;
    if (rot_cl && rdt_named) {
      const std::string rv = pin->GetString("rad_vet","vet_gd_rotate_dt");
      if (rv.compare("auto") == 0) {
        vgd_rot_auto = true;
      } else {
        vgd_rot_dt = pin->GetReal("rad_vet","vet_gd_rotate_dt");
      }
      rot_phys = vgd_rot_auto || (vgd_rot_dt > 0.0);
    } else if (vgd_on && !rev_named && fresh) {
      vgd_rot_auto = vgd_rot_unc = true;
      rot_phys = true;
    } else if (vet_sc && !sc_pin && !rev_named && fresh) {
      vgd_rot_auto = vgd_rot_unc = true;
      rot_phys = true;
    }
    if (rot_phys && vet_sc && sc_pin) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<rad_vet>/vet_gd_rotate_dt: the vet_sc set is pinned by "
                << "boundary beams (implicit_beam_n); no rotation"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rot_phys) {
      vsc_rot = vet_sc;
      vgd_rot = 0;
      pin_rot_ = pin;
      if (vet_sc && !cln_ring) {
        recipe::Rule(pin, "rad_vet", "vet_gd_rotate_slope", "correct");
        recipe::Rule(pin, "rad_vet", "vet_gd_rotate_slope_frac",
                     ParameterInput::RealToString(1.0 - 1.0/std::sqrt(2.0)));
      }
      if (vgd_rot_auto) {
        vgd_rot_auto_fac = 20.0;
        // recorded once resolved; until then the key reads "auto"
        pin->SetString("rad_vet","vet_gd_rotate_dt","auto");
      }
    }
  }
  if (VgdRotOn()) {
    recipe::Rule(pin, "rad_vet", "vet_rotate_advance", g_vra_adv ? "true" : "false");
  }
  // the axis: sp turns the set about the polar axis (exact wall maps for HEALPix), the
  // cubed sphere about a quasi-random axis per rotation block, so no direction stays
  // aligned with panel seams
  if (vgd_on) {
    vgd_rot3 = pmy_pack->pmesh->use_cubed_sphere;
    recipe::Rule(pin, "rad_vet", "vet_gd_rotate_axis", vgd_rot3 ? "3d" : "z");
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetGdAngle
//! \brief vet_gd_rotate_every = N > 0: the z-rotation of the set for cycle c, constant
//! per block of N cycles: alpha = P frac((c/N) g), g = (sqrt 5 - 1)/2, P the set's
//! period (deterministic: reruns and restarts bitwise)

Real RadiationVet::VetGdAngle(const int cyc) const {
  if (!VgdRotOn()) {return 0.0;}
  return VetGdAngleBlk(VetGdBlock(cyc));
}

Real RadiationVet::VetGdAngleBlk(const int64_t blk) const {
  // vet_gd_rotate_axis = 3d: the key of block b (decoded by VgdRot3) is zigzag(b), or
  // zigzag(b mod K) with vet_rotate_avg = K. vet_rotate_advance: always zigzag(b), so
  // any K consecutive blocks are consecutive R_3 points and the window moves by one set
  // per block. vet_rotate_avg_ema: a never-repeating sequence (3d: key of b; z: golden
  // phase of b)
  if (vra_ema) {
    if (vgd_rot3) {
      const int64_t z = (blk >= 0) ? 2*blk : (-2*blk - 1);
      return static_cast<Real>(z);
    }
    const double per = vgd_zper;
    const double b = static_cast<double>(blk);
    const double g = 0.5*(std::sqrt(5.0) - 1.0);
    const double a = per*(b*g - std::floor(b*g));
    return static_cast<Real>(a);
  }
  if (vgd_rot3) {
    int64_t r = blk;
    if (vra_on && vra_k >= 2 && !g_vra_adv) {r = ((blk % vra_k) + vra_k) % vra_k;}
    const int64_t z = (r >= 0) ? 2*r : (-2*r - 1);
    return static_cast<Real>(z);
  }
  // vet_rotate_avg = K: angle = period x RotPhase, (b mod K)/K, so any K consecutive
  // blocks are K equispaced rotations (the K-fold refined azimuth set)
  if (vra_on) {
    const double per = vgd_zper;
    return static_cast<Real>(per*RotPhase(blk));
  }
  const double b = static_cast<double>(blk);
  const double g = 0.5*(std::sqrt(5.0) - 1.0);
  // a product set keeps pi/2 here
  return ((vgd_gl_nmu > 0) ? 0.5*M_PI : vgd_zper)*(b*g - std::floor(b*g));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdRot3
//! \brief vet_gd_rotate_axis = 3d: the rotation matrix (row major) of block key `key`.
//! Shoemake (1992) uniform rotation from the R_3 point u_b = frac(b (1/g, 1/g^2,
//! 1/g^3)) (Roberts 2018), times the inverse of b = 0's (block 0 unrotated). A pure
//! function of the key: same bits on every rank and restart.

void RadiationVet::VgdRot3(const Real key, double r[9]) const {
  const int64_t z = static_cast<int64_t>(std::llround(static_cast<double>(key)));
  const int64_t b = (z % 2 == 0) ? (z/2) : (-(z + 1)/2);
  auto quat = [](const int64_t bb, double q[4]) {
    const double g = 1.2207440846057594754;   // the real root of x^4 = x + 1
    const double a[3] = {1.0/g, 1.0/(g*g), 1.0/(g*g*g)};
    double u[3];
    for (int c = 0; c < 3; ++c) {
      const double x = static_cast<double>(bb)*a[c];
      u[c] = x - std::floor(x);
    }
    const double s1 = std::sqrt(1.0 - u[0]), s2 = std::sqrt(u[0]);
    q[0] = s1*std::sin(2.0*M_PI*u[1]);
    q[1] = s1*std::cos(2.0*M_PI*u[1]);
    q[2] = s2*std::sin(2.0*M_PI*u[2]);
    q[3] = s2*std::cos(2.0*M_PI*u[2]);
  };
  auto mat = [](const double q[4], double m[9]) {   // q = (x, y, z, w), unit
    const double x = q[0], y = q[1], zz = q[2], w = q[3];
    m[0] = 1.0 - 2.0*(y*y + zz*zz); m[1] = 2.0*(x*y - zz*w);  m[2] = 2.0*(x*zz + y*w);
    m[3] = 2.0*(x*y + zz*w);  m[4] = 1.0 - 2.0*(x*x + zz*zz); m[5] = 2.0*(y*zz - x*w);
    m[6] = 2.0*(x*zz - y*w);  m[7] = 2.0*(y*zz + x*w);  m[8] = 1.0 - 2.0*(x*x + y*y);
  };
  if (b == 0) {
    for (int c = 0; c < 9; ++c) {r[c] = (c % 4 == 0) ? 1.0 : 0.0;}
    return;
  }
  double qb[4], q0[4], mb[9], m0[9];
  quat(b, qb);
  quat(0, q0);
  mat(qb, mb);
  mat(q0, m0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double v = 0.0;
      for (int k = 0; k < 3; ++k) {v += mb[3*i + k]*m0[3*j + k];}   // mb m0^T
      r[3*i + j] = v;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::RotPhase
//! \brief rotation of block blk as a fraction of the set's period in [0, 1): golden
//! frac(b g), or with vet_rotate_avg = K the stratified (b mod K)/K

double RadiationVet::RotPhase(const int64_t blk) const {
  // vet_rotate_advance: frac(b (1 + gam)/K), gam = g/(4 (K-1)): any K consecutive phases
  // are equispaced to within 0.16/K of the period, and the window drifts by gam/K per
  // block without repeating
  if (vra_on && vra_k >= 2 && g_vra_adv) {
    const double g = 0.5*(std::sqrt(5.0) - 1.0);
    const double kk = static_cast<double>(vra_k);
    const double x = static_cast<double>(blk)*(1.0 + g/(4.0*(kk - 1.0)))/kk;
    return x - std::floor(x);
  }
  if (vra_on && vra_k >= 2) {
    const int64_t r = ((blk % vra_k) + vra_k) % vra_k;
    return static_cast<double>(r)/static_cast<double>(vra_k);
  }
  const double b = static_cast<double>(blk);
  const double g = 0.5*(std::sqrt(5.0) - 1.0);
  return b*g - std::floor(b*g);
}

//----------------------------------------------------------------------------------------
//! \fn int64_t RadiationVet::VetGdBlock
//! \brief rotation block of the step of cycle cyc. Cycle clock (vet_gd_rotate_every =
//! N): cyc/N. Physical clock (vet_gd_rotate_dt > 0): floor((t + dt/2 - t0)/dt_rot),
//! the block of the step midpoint, shared by all stages and builds of the step. With
//! vet_gd_rotate_align every jump sits on a step boundary. Restart-safe.

int64_t RadiationVet::VetGdBlock(const int cyc) const {
  if (vgd_rot_dt > 0.0) {
    const Mesh *pm = pmy_pack->pmesh;
    const Real dt = (pm->dt > 0.0 && pm->dt < vgd_rot_dt) ? pm->dt : 0.0;
    return static_cast<int64_t>(std::floor((pm->time + 0.5*dt)/vgd_rot_dt));
  }
  return (vgd_rot > 0) ? static_cast<int64_t>(cyc/vgd_rot) : 0;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::RotLastStep
//! \brief is the current step the last of its rotation block? Physical clock: t + dt
//! reaches t0 + (b+1) dt_rot; cycle clock: (cyc+1) % N == 0. Shared with vet_sc.

bool RadiationVet::RotLastStep() const {
  const Mesh *pm = pmy_pack->pmesh;
  if (vgd_rot_dt > 0.0) {
    const int64_t b = VetGdBlock(pm->ncycle);
    const double tk = (static_cast<double>(b) + 1.0)*vgd_rot_dt;
    return pm->time + pm->dt >= tk - 1.0e-9*vgd_rot_dt;
  }
  if (vgd_rot > 0) {return ((pm->ncycle + 1) % vgd_rot) == 0;}
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdRotStep
//! \brief (Driver, before Time2RotCorrect) on the first step of a new rotation block
//! the stored hesdirk2 slope K1 belongs to the old set; be resets it (t2_ok = false).
//! The first call takes the block before t as the last one (restart on a boundary).

void RadiationVet::VetGdRotStep() {
  if (!(vgd_rot_dt > 0.0) || !(vgd_on || vsc_rot)) {return;}
  const int64_t b = VetGdBlock(pmy_pack->pmesh->ncycle);
  if (vgd_rot_lastblk == INT64_MIN) {
    const Mesh *pm = pmy_pack->pmesh;
    vgd_rot_lastblk = static_cast<int64_t>(
        std::floor(pm->time/vgd_rot_dt - 1.0e-9));
  }
  // a Cartesian set without the ring keeps the slope and asks Time2RotCorrect for its
  // jump; the cross-fade has no jump; clean_rot_ring = jump: a backward-Euler step
  if (vgd_xf && vgd_on) {
    vgd_rot_lastblk = b;
    return;
  }
  if (b != vgd_rot_lastblk && time_scheme == VET_TIME_HESDIRK2) {
    if (vet_sc && !cln_ring && Time2Active()) t2_rotcorr = true;
    if (cln_ring == 2) t2_ok = false;
  }
  vgd_rot_lastblk = b;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetGdAlignDt
//! \brief vet_gd_rotate_align (Mesh::NewTimeStep, after the global min): reach the next
//! boundary k dt_rot in n = ceil(rem/dt) equal steps, dt' = rem/n >= dt/2.
//! vet_gd_rotate_stretch = s > 0: n - 1 steps when dt' <= (1 + s) dt.

Real RadiationVet::VetGdAlignDt(const Real time, const Real dt) {
  // vet_gd_rotate_dt = auto: vgd_rot_auto_fac x the first global dt, recorded as a
  // number (restarts keep it)
  if (vgd_rot_auto && !(vgd_rot_dt > 0.0) && dt > 0.0) {
    vgd_rot_dt = vgd_rot_auto_fac*(vgd_rot_unc ? pmy_pack->pmesh->dt_uncapped : dt);
    vgd_rot_auto = false;
    if (pin_rot_ != nullptr) {pin_rot_->SetReal("rad_vet","vet_gd_rotate_dt",vgd_rot_dt);}
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_dt = auto: " << vgd_rot_auto_fac
                << " x the first "
                << (vgd_rot_unc ? "CFL dt (no dt_fixed / dt_retry cap) " : "dt ")
                << (vgd_rot_unc ? pmy_pack->pmesh->dt_uncapped : dt) << " = "
                << ParameterInput::RealToString(vgd_rot_dt)
                << " (recorded)" << std::endl;
    }
  }
  if (!(vgd_rot_dt > 0.0) || !(dt > 0.0)) {return dt;}
  const double x = time/vgd_rot_dt;
  double k = std::floor(x) + 1.0;
  // already on a boundary (to round-off): the next one
  if ((k - x)*vgd_rot_dt < 1.0e-9*dt) {k += 1.0;}
  const double rem = k*vgd_rot_dt - time;
  double n = std::ceil(rem/dt*(1.0 - 1.0e-9));
  if (n > 1.0e6) {return dt;}
  if (n > 1.0 && rem/(n - 1.0) <= 1.02*dt) n -= 1.0;   // (dt may grow by 2%)
  const Real dtn = static_cast<Real>(rem/std::max(n, 1.0));
  if (dtn < dt) {vgd_rot_nalign += 1.0;}
  return dtn;
}

//----------------------------------------------------------------------------------------
//! \fn static void VgdMirPick
//! \brief vet_gd_theta_wall = mirror: set directions standing in for direction m at a
//! cell with radial unit vector (rx, ry, rz): an exact match, else the nearest of m's
//! branch (sign of m . r_hat); with wint the 3 nearest, barycentric weights when all
//! >= 0, else inverse-angle (vet_gd_wall_interp's rule)

static void VgdMirPick(const std::vector<double> &d, const int n, const double mx,
                       const double my, const double mz, const double rx,
                       const double ry, const double rz, const bool wint, int r3[3],
                       double w3[3], const bool ronly = false) {
  const bool up = (mx*rx + my*ry + mz*rz) >= 0.0;
  double b3[3] = {-2.0, -2.0, -2.0};
  int best = 0;
  double bd = -2.0;
  r3[0] = r3[1] = r3[2] = 0;
  for (int r = 0; r < n; ++r) {
    const double dd = mx*d[4*r] + my*d[4*r+1] + mz*d[4*r+2];
    if (dd > bd) {bd = dd; best = r;}
    const double br = d[4*r]*rx + d[4*r+1]*ry + d[4*r+2]*rz;
    if ((br >= 0.0) != up) {continue;}
    if (dd > b3[0]) {
      b3[2] = b3[1]; r3[2] = r3[1]; b3[1] = b3[0]; r3[1] = r3[0];
      b3[0] = dd; r3[0] = r;
    } else if (dd > b3[1]) {
      b3[2] = b3[1]; r3[2] = r3[1]; b3[1] = dd; r3[1] = r;
    } else if (dd > b3[2]) {
      b3[2] = dd; r3[2] = r;
    }
  }
  w3[0] = 1.0;
  w3[1] = w3[2] = 0.0;
  if (bd > 1.0 - 1.0e-12) {r3[0] = r3[1] = r3[2] = best; return;}
  if (b3[0] <= -2.0) {r3[0] = r3[1] = r3[2] = best; return;}
  if (!wint || b3[2] <= -2.0) {r3[1] = r3[2] = r3[0]; return;}
  // m on a ring of the set (same n_z): linear in azimuth between the two bracketing
  // ring directions (n_z kept exactly)
  {
    int ra = -1, rb = -1;
    double da = -2.0, db = -2.0;
    for (int r = 0; r < n; ++r) {
      if (std::fabs(d[4*r+2] - mz) > 1.0e-10) {continue;}
      const double br = d[4*r]*rx + d[4*r+1]*ry + d[4*r+2]*rz;
      if ((br >= 0.0) != up) {continue;}
      const double dd = mx*d[4*r] + my*d[4*r+1] + mz*d[4*r+2];
      if (dd > da) {
        rb = ra; db = da; ra = r; da = dd;
      } else if (dd > db) {
        rb = r; db = dd;
      }
    }
    if (ra >= 0 && rb >= 0) {
      const double pt = std::atan2(my, mx);
      auto dang = [pt](const double x, const double y) {
        double a = std::atan2(y, x) - pt;
        a -= 2.0*M_PI*std::floor((a + M_PI)/(2.0*M_PI));
        return a;
      };
      const double aa = dang(d[4*ra], d[4*ra+1]), ab = dang(d[4*rb], d[4*rb+1]);
      if (aa*ab <= 0.0 && std::fabs(aa - ab) > 1.0e-12) {   // bracketing pair
        r3[0] = ra; r3[1] = rb; r3[2] = ra;
        w3[0] = std::fabs(ab)/std::fabs(aa - ab);
        w3[1] = 1.0 - w3[0];
        w3[2] = 0.0;
        return;
      }
    }
    // ronly (vet_gd_wall_sym): no bracketing pair -> the nearest same-branch ring
    // direction alone (paired under both mirrors)
    if (ronly && ra >= 0) {r3[0] = r3[1] = r3[2] = ra; return;}
  }
  const double *u = &d[4*r3[0]], *v = &d[4*r3[1]], *x = &d[4*r3[2]];
  auto det3 = [](const double *a, const double *b, const double *c) {
    return a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0])
           + a[2]*(b[0]*c[1] - b[1]*c[0]);
  };
  const double mv[3] = {mx, my, mz};
  const double dt0 = det3(u, v, x);
  bool ok = std::fabs(dt0) > 1.0e-12;
  double ww[3] = {0.0, 0.0, 0.0};
  if (ok) {
    ww[0] = det3(mv, v, x)/dt0;
    ww[1] = det3(u, mv, x)/dt0;
    ww[2] = det3(u, v, mv)/dt0;
    ok = (ww[0] >= -1.0e-12) && (ww[1] >= -1.0e-12) && (ww[2] >= -1.0e-12);
  }
  if (!ok) {
    for (int t = 0; t < 3; ++t) {
      ww[t] = 1.0/std::max(std::acos(std::min(b3[t], 1.0)), 1.0e-9);
    }
  }
  const double ws = std::max(ww[0], 0.0) + std::max(ww[1], 0.0) + std::max(ww[2], 0.0);
  for (int t = 0; t < 3; ++t) {w3[t] = std::max(ww[t], 0.0)/ws;}
}

//----------------------------------------------------------------------------------------
//! \fn std::vector<double> RadiationVet::VgdRotate
//! \brief the set at angle (or block key) alpha, packed as vgd_base (x, y, z, w):
//! rotated about z, or with vet_gd_rotate_axis = 3d by R_b (VgdRot3)

std::vector<double> RadiationVet::VgdRotate(const Real alpha) const {
  const int n = vgd_n;
  std::vector<double> d(vgd_base);
  if (vgd_rot3) {
    double rm[9];
    VgdRot3(alpha, rm);
    for (int q = 0; q < n; ++q) {
      const double x = vgd_base[4*q], y = vgd_base[4*q+1], z = vgd_base[4*q+2];
      double v[3];
      for (int c = 0; c < 3; ++c) {v[c] = rm[3*c]*x + rm[3*c+1]*y + rm[3*c+2]*z;}
      const double nv = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      for (int c = 0; c < 3; ++c) {d[4*q+c] = v[c]/nv;}
    }
  } else {
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    for (int q = 0; q < n; ++q) {
      const double x = vgd_base[4*q], y = vgd_base[4*q+1];
      d[4*q] = ca*x - sa*y;
      d[4*q+1] = sa*x + ca*y;
    }
  }
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdTabCompute
//! \brief host part of VetGdTables for angle alpha (rotated set, wall maps, branch
//! tables) into the pinned host buffers vgd_tabh; touches nothing else (safe on the
//! prefetch thread). nthr: host threads over (m, k) columns (result independent).

void RadiationVet::VgdTabCompute(const Real alpha, const int nthr) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c2 = indcs.nx2 + 2*vgd_w;
  const int c3 = indcs.nx3 + 2*vgd_w;
  const int n = vgd_n;
  VgdTabHost &H = vgd_tabh;
  std::vector<double> d = VgdRotate(alpha);
  auto &d_h = H.dir;
  for (int q = 0; q < n; ++q) {
    for (int c = 0; c < 4; ++c) {d_h(q,c) = d[4*q+c];}
  }
  if (vgd_csg) {   // no walls; the branch of every band cell
    VetGdCsTab(d);
    H.alpha = alpha;
    return;
  }
  // wall maps: lateral ghost cells outside the mesh in theta or phi hold the other
  // side's intensities; new(d) = old(map(d)), map = nearest set direction to M n_d
  auto &w_h = H.w;
  auto &mp_h = H.mp;
  const bool wint = vgd_wint;
  auto &m3_h = H.m3;
  auto &w3_h = H.w3;
  std::vector<double> &mt_v = H.mt;
  mt_v.assign(wint ? static_cast<size_t>(nmb)*c3*c2*n : 0, 0.0);
  auto mt_h = [&](int m, int k, int j, int q) -> double & {
    return mt_v[((static_cast<size_t>(m)*c3 + k)*c2 + j)*n + q];
  };
  // positions of the wall maps and branch tables: the face MIDPOINTS (x2v is the
  // volume centroid on sp; the periodic image of a midpoint is exact)
  auto &mbs = pmy_pack->pmb->mb_size;
  const int jsg = vgd_w, ksg = vgd_w;   // the band-index origin of the active cells
  auto x2_h = [&](const int m, const int j) {
    // stretched theta / full sphere: the node table (sph2.inc)
    if (vgd_sth) {return static_cast<Real>(vgd_thn_h[static_cast<size_t>(m)*c2 + j]);}
    return vgd_lg.ThH(mbs, m, j, jsg);
  };
  // full sphere: rows beyond a pole are the same points (extended chart), not a
  // periodic image: no theta wall
  const bool thwrap = !vgd_pol;
  auto x3_h = [&](const int m, const int k) {
    // vet_gd_band_global: the band column's phi from the global index
    if (vgd_gband) {return static_cast<Real>(vgd_phc_h[static_cast<size_t>(m)*c3 + k]);}
    return vgd_lg.PhH(mbs, m, k, ksg);
  };
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  const double p0 = pm->mesh_size.x3min, p1 = pm->mesh_size.x3max;
  int nwall = 0, nexact = 0, nmap = 0, nphx = 0, nthx = 0;
  // the (m, k) columns are independent: host threads over them (bitwise identical to
  // the serial loop; each thread keeps its own counters)
  auto column = [&](const int m, const int k, int &nwall, int &nexact, int &nmap,
                    int &nphx, int &nthx) {
    {
      for (int j = 0; j < c2; ++j) {
        const double th = x2_h(m,j), ph = x3_h(m,k);
        double thw = th, phw = ph;
        while (thwrap && thw < t0) {thw += t1 - t0;}
        while (thwrap && thw > t1) {thw -= t1 - t0;}
        while (phw < p0) {phw += p1 - p0;}
        while (phw > p1) {phw -= p1 - p0;}
        const bool wall = (thw != th) || (phw != ph);
        w_h(m,k,j) = wall ? 1 : 0;
        for (int q = 0; q < n; ++q) {mp_h(m,k,j,q) = q;}
        if (!wall) {continue;}
        ++nwall;
        const double st = std::sin(th), ct = std::cos(th);
        const double sp = std::sin(ph), cp = std::cos(ph);
        const double sw = std::sin(thw), cw = std::cos(thw);
        const double spw = std::sin(phw), cpw = std::cos(phw);
        for (int q = 0; q < n; ++q) {
          const double nx = d[4*q], ny = d[4*q+1], nz = d[4*q+2];
          const double a = nx*st*cp + ny*st*sp + nz*ct;
          const double b = nx*ct*cp + ny*ct*sp - nz*st;
          const double c = -nx*sp + ny*cp;
          const double mx = a*sw*cpw + b*cw*cpw - c*spw;
          const double my = a*sw*spw + b*cw*spw + c*cpw;
          const double mz = a*cw - b*sw;
          int best = 0;
          double bd = -2.0;
          for (int r = 0; r < n; ++r) {
            const double dd = mx*d[4*r] + my*d[4*r+1] + mz*d[4*r+2];
            if (dd > bd) {bd = dd; best = r;}
          }
          mp_h(m,k,j,q) = best;
          ++nmap;
          // vet_gd_band_global copies without re-indexing: exact = identity
          if (bd > 1.0 - 1.0e-12 && (!vgd_gband || best == q)) {
            ++nexact;
          } else if (thw == th) {
            ++nphx;   // a phi-only image the set does not map onto itself
          } else if (phw == ph) {
            ++nthx;   // a theta-only image (a turn about phi_hat: never exact)
          }
          if (wint) {
            // vet_gd_wall_interp: M n_q from the 3 nearest same-branch set directions at
            // the source cell (barycentric when all >= 0, else inverse-angle, sum 1); an
            // exact match keeps the single direction
            const double mt = nx*st*cp + ny*st*sp + nz*ct;   // target branch n.r_hat
            int r3[3] = {best, best, best};
            double w3[3] = {1.0, 0.0, 0.0};
            if (vgd_wsym && thw == th) {
              // vet_gd_wall_sym, phi wall: VgdMirPick, linear in azimuth on n_q's ring
              // (positive, mirror-paired)
              VgdMirPick(d, n, mx, my, mz, sw*cpw, sw*spw, cw, true, r3, w3, true);
            } else if (bd <= 1.0 - 1.0e-12) {
              double b3[3] = {-2.0, -2.0, -2.0};
              for (int r = 0; r < n; ++r) {
                const double br = d[4*r]*sw*cpw + d[4*r+1]*sw*spw + d[4*r+2]*cw;
                if ((br >= 0.0) != (mt >= 0.0)) {continue;}
                const double dd = mx*d[4*r] + my*d[4*r+1] + mz*d[4*r+2];
                if (dd > b3[0]) {
                  b3[2] = b3[1]; r3[2] = r3[1]; b3[1] = b3[0]; r3[1] = r3[0];
                  b3[0] = dd; r3[0] = r;
                } else if (dd > b3[1]) {
                  b3[2] = b3[1]; r3[2] = r3[1]; b3[1] = dd; r3[1] = r;
                } else if (dd > b3[2]) {
                  b3[2] = dd; r3[2] = r;
                }
              }
              if (b3[2] > -2.0) {
                // solve [n_r0 n_r1 n_r2] w = m (Cramer)
                const double *u = &d[4*r3[0]], *v = &d[4*r3[1]], *x = &d[4*r3[2]];
                auto det3 = [](const double *a, const double *b, const double *c) {
                  return a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0])
                         + a[2]*(b[0]*c[1] - b[1]*c[0]);
                };
                const double mv[3] = {mx, my, mz};
                const double dt0 = det3(u, v, x);
                bool ok = std::fabs(dt0) > 1.0e-12;
                double ww[3] = {0.0, 0.0, 0.0};
                if (ok) {
                  ww[0] = det3(mv, v, x)/dt0;
                  ww[1] = det3(u, mv, x)/dt0;
                  ww[2] = det3(u, v, mv)/dt0;
                  ok = (ww[0] >= -1.0e-12) && (ww[1] >= -1.0e-12) && (ww[2] >= -1.0e-12);
                }
                if (!ok) {
                  for (int t = 0; t < 3; ++t) {
                    ww[t] = 1.0/std::max(std::acos(std::min(b3[t], 1.0)), 1.0e-9);
                  }
                }
                const double ws = std::max(ww[0], 0.0) + std::max(ww[1], 0.0)
                                  + std::max(ww[2], 0.0);
                for (int t = 0; t < 3; ++t) {w3[t] = std::max(ww[t], 0.0)/ws;}
              } else {
                r3[0] = r3[1] = r3[2] = best;
              }
            }
            for (int t = 0; t < 3; ++t) {
              m3_h(m,k,j,q,t) = r3[t];
              w3_h(m,k,j,q,t) = w3[t];
            }
            mt_h(m,k,j,q) = mt;
          }
        }
      }
    }
  };
  {
    const int ncol = nmb*c3;
    const int hwc = static_cast<int>(std::thread::hardware_concurrency());
    const int nth = std::max(1, std::min(hwc, std::min(nthr, ncol)));
    std::vector<std::array<int, 5>> cnt(nth, {0, 0, 0, 0, 0});
    std::vector<std::thread> pool;
    for (int t = 0; t < nth; ++t) {
      pool.emplace_back([&, t]() {
        for (int c = t; c < ncol; c += nth) {
          column(c/c3, c % c3, cnt[t][0], cnt[t][1], cnt[t][2], cnt[t][3],
                 cnt[t][4]);
        }
      });
    }
    for (auto &th : pool) {th.join();}
    for (int t = 0; t < nth; ++t) {
      nwall += cnt[t][0];
      nexact += cnt[t][1];
      nmap += cnt[t][2];
      nphx += cnt[t][3];
      nthx += cnt[t][4];
    }
  }
  // branch of every stored value: n . r_hat of the direction stored at a lateral cell
  // (for a wall ghost the mapped one, at the source position)
  auto &mr_h = H.mr;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        double th = x2_h(m,j), ph = x3_h(m,k);
        if (thwrap && th < t0) {th += (t1 - t0);}
        if (thwrap && th > t1) {th -= (t1 - t0);}
        while (ph < p0) {ph += p1 - p0;}
        while (ph > p1) {ph -= p1 - p0;}
        const double st = std::sin(th), ct = std::cos(th);
        const double sp = std::sin(ph), cp = std::cos(ph);
        for (int q = 0; q < n; ++q) {
          const int r = mp_h(m,k,j,q);
          mr_h(m,k,j,q) = d[4*r]*st*cp + d[4*r+1]*st*sp + d[4*r+2]*ct;
          // interpolated wall ghost: the branch of the target direction (sources were
          // chosen on it)
          if (wint && w_h(m,k,j) != 0) {mr_h(m,k,j,q) = mt_h(m,k,j,q);}
        }
      }
    }
  }
  if (vgd_thmir) {
    // vet_gd_theta_wall = mirror: theta ghost rows (every k, phi unwrapped) read the
    // mirror row; their direction tables and branch (the target's: the local
    // reflection keeps n . r_hat)
    const int wb = vgd_w, nx2 = indcs.nx2;
    for (int m = 0; m < nmb; ++m) {
      for (int k = 0; k < c3; ++k) {
        const double ph = x3_h(m,k);
        const double sp = std::sin(ph), cp = std::cos(ph);
        for (int sd = 0; sd < 2; ++sd) {
          for (int g = 0; g < wb; ++g) {
            const int j = (sd == 1) ? (wb + nx2 + g) : (wb - 1 - g);
            const int jm = VgdMirRow(g, sd == 1, wb, nx2);
            const double sb = ((g/nx2) % 2 == 0) ? -1.0 : 1.0;   // even: no flip
            const double th = x2_h(m,j), thm = x2_h(m,jm);
            const double st = std::sin(th), ct = std::cos(th);
            const double sw = std::sin(thm), cw = std::cos(thm);
            const double rx = sw*cp, ry = sw*sp, rz = cw;
            w_h(m,k,j) = 2;
            for (int q = 0; q < n; ++q) {
              const double nx = d[4*q], ny = d[4*q+1], nz = d[4*q+2];
              const double a = nx*st*cp + ny*st*sp + nz*ct;
              const double b = nx*ct*cp + ny*ct*sp - nz*st;
              const double c = -nx*sp + ny*cp;
              // (a, -b, c) in the mirror cell's frame
              const double mx = a*sw*cp + sb*b*cw*cp - c*sp;
              const double my = a*sw*sp + sb*b*cw*sp + c*cp;
              const double mz = a*cw - sb*b*sw;
              int r3[3];
              double w3[3];
              VgdMirPick(d, n, mx, my, mz, rx, ry, rz, wint, r3, w3);
              // the theta ghost's own wall tables (flag 2 = mirror): one direction, or 3
              // with weights (wint)
              mp_h(m,k,j,q) = r3[0];
              if (wint) {
                for (int t = 0; t < 3; ++t) {
                  m3_h(m,k,j,q,t) = r3[t];
                  w3_h(m,k,j,q,t) = w3[t];
                }
              }
              // without wint the branch of the value read (r3[0] in the mirror row),
              // else a band corner would read that pass's other, unswept branch
              mr_h(m,k,j,q) = wint ? a : mr_h(m,k,jm,r3[0]);
            }
          }
        }
      }
    }
  }
  if (vgd_thvac > 0) {
    // vet_gd_theta_wall = vacuum (flag 3) / extrap (flag 4): theta ghost rows beyond a
    // physical wall; only the flag changes. The periodic maps stay for the twin (which
    // continues the uniform shell state beyond the wall) and its periodic halo reads
    const int wb = vgd_w, nx2 = indcs.nx2;
    const int fl = (vgd_thvac == 1) ? 3 : 4;
    for (int m = 0; m < nmb; ++m) {
      for (int k = 0; k < c3; ++k) {
        for (int sd = 0; sd < 2; ++sd) {
          for (int g = 0; g < wb; ++g) {
            const int j = (sd == 1) ? (wb + nx2 + g) : (wb - 1 - g);
            const double th = x2_h(m,j);
            if (th < t0 || th > t1) {w_h(m,k,j) = fl;}   // not an interior boundary
          }
        }
      }
    }
  }
  H.alpha = alpha;
  H.nwall = nwall;
  H.nexact = nexact;
  H.nmap = nmap;
  H.nphx = nphx;
  H.nthx = nthx;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdTabAlloc
//! \brief pinned host buffers of the tables (once, zero-filled) and, with dev, the
//! second set of device tables the prefetch uploads into

void RadiationVet::VgdTabAlloc(const bool dev) {
  VgdTabHost &H = vgd_tabh;
  if (H.dir.extent(0) != vgd_dir.extent(0)) {
    using PS = Kokkos::SharedHostPinnedSpace;
    H.dir = Kokkos::View<Real **, LayoutWrapper, PS>("vet_vgd_tab_dir", vgd_dir.extent(0),
                                                     vgd_dir.extent(1));
    H.w = Kokkos::View<int ***, LayoutWrapper, PS>("vet_vgd_tab_w", vgd_wall.extent(0),
                                                   vgd_wall.extent(1),
                                                   vgd_wall.extent(2));
    H.mp = Kokkos::View<int ****, LayoutWrapper, PS>("vet_vgd_tab_mp", vgd_map.extent(0),
                                                     vgd_map.extent(1), vgd_map.extent(2),
                                                     vgd_map.extent(3));
    H.mr = Kokkos::View<Real ****, LayoutWrapper, PS>("vet_vgd_tab_mr", vgd_mr.extent(0),
                                                      vgd_mr.extent(1), vgd_mr.extent(2),
                                                      vgd_mr.extent(3));
    if (vgd_wint) {
      H.m3 = Kokkos::View<int *****, LayoutWrapper, PS>("vet_vgd_tab_m3",
               vgd_m3.extent(0), vgd_m3.extent(1), vgd_m3.extent(2), vgd_m3.extent(3),
               vgd_m3.extent(4));
      H.w3 = Kokkos::View<Real *****, LayoutWrapper, PS>("vet_vgd_tab_w3",
               vgd_w3.extent(0), vgd_w3.extent(1), vgd_w3.extent(2), vgd_w3.extent(3),
               vgd_w3.extent(4));
    }
  }
  if (dev && vgd_tabd.dir.extent(0) != vgd_dir.extent(0)) {
    Kokkos::realloc(vgd_tabd.dir, vgd_dir.extent(0), vgd_dir.extent(1));
    Kokkos::realloc(vgd_tabd.wall, vgd_wall.extent(0), vgd_wall.extent(1),
                    vgd_wall.extent(2));
    Kokkos::realloc(vgd_tabd.map, vgd_map.extent(0), vgd_map.extent(1), vgd_map.extent(2),
                    vgd_map.extent(3));
    Kokkos::realloc(vgd_tabd.mr, vgd_mr.extent(0), vgd_mr.extent(1), vgd_mr.extent(2),
                    vgd_mr.extent(3));
    if (vgd_wint) {
      Kokkos::realloc(vgd_tabd.m3, vgd_m3.extent(0), vgd_m3.extent(1), vgd_m3.extent(2),
                      vgd_m3.extent(3), vgd_m3.extent(4));
      Kokkos::realloc(vgd_tabd.w3, vgd_w3.extent(0), vgd_w3.extent(1), vgd_w3.extent(2),
                      vgd_w3.extent(3), vgd_w3.extent(4));
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTabPrefetch
//! \brief compute the next rotation block's tables on a host thread (VgdTabCompute,
//! same numbers); VetGdTabPoll uploads them into the second device set, VetGdTables
//! swaps them in at the next rotation

void RadiationVet::VetGdTabPrefetch(const int cyc) {
  if (!VgdRotOn() || vgd_tab_state != 0 || vgd_xf) {return;}
  const Real an = VetGdAngleBlk(VetGdBlock(cyc) + 1);
  if (an == vgd_alpha) {return;}
  VgdTabAlloc(true);
  vgd_tab_state = 1;
  vgd_tabf = std::async(std::launch::async, [this, an]() {VgdTabCompute(an, 8);});
}

void RadiationVet::VetGdTabPoll() {
  if (vgd_tab_state != 1) {return;}
  if (vgd_tabf.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {return;}
  vgd_tabf.get();
  // the second device set is not read by anything in flight (swapped out at the last
  // rotation, fenced since): upload on the table instance
  VgdTabHost &H = vgd_tabh;
  Kokkos::deep_copy(vgd_tex, vgd_tabd.dir, H.dir);
  Kokkos::deep_copy(vgd_tex, vgd_tabd.wall, H.w);
  Kokkos::deep_copy(vgd_tex, vgd_tabd.map, H.mp);
  Kokkos::deep_copy(vgd_tex, vgd_tabd.mr, H.mr);
  if (vgd_wint) {
    Kokkos::deep_copy(vgd_tex, vgd_tabd.m3, H.m3);
    Kokkos::deep_copy(vgd_tex, vgd_tabd.w3, H.w3);
  }
  vgd_tab_state = 3;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTables
//! \brief the set rotated by alpha -> vgd_dir, wall maps and branch tables (from the
//! prefetched device set when it holds alpha, else computed now and copied)

// angle of the tables now in vgd_dir/wall/map/mr/m3/w3 (not a class member: the
// dt-retry checkpoint restores vgd_alpha but not these Views)
static Real g_vgd_tab_alpha = -9.0;

void RadiationVet::VetGdTables(const Real alpha) {
  VetGdTabPoll();
  VgdTabHost &H = vgd_tabh;
  if (vgd_tab_state == 3 && H.alpha == alpha) {
    vgd_tex.fence();
    std::swap(vgd_dir, vgd_tabd.dir);
    std::swap(vgd_wall, vgd_tabd.wall);
    std::swap(vgd_map, vgd_tabd.map);
    std::swap(vgd_mr, vgd_tabd.mr);
    if (vgd_wint) {
      std::swap(vgd_m3, vgd_tabd.m3);
      std::swap(vgd_w3, vgd_tabd.w3);
    }
    vgd_tab_nswap += 1.0;
  } else {
    if (vgd_tab_state == 1) {vgd_tabf.get();}     // another angle: not used
    if (vgd_tab_state == 3) {vgd_tex.fence();}    // the host buffers are being read
    VgdTabAlloc(false);
    if (H.alpha != alpha) {VgdTabCompute(alpha, 16);}
    Kokkos::deep_copy(vgd_dir, H.dir);
    Kokkos::deep_copy(vgd_mr, H.mr);
    if (vgd_wint) {
      Kokkos::deep_copy(vgd_m3, H.m3);
      Kokkos::deep_copy(vgd_w3, H.w3);
    }
    Kokkos::deep_copy(vgd_wall, H.w);
    Kokkos::deep_copy(vgd_map, H.mp);
  }
  vgd_tab_state = 0;
  const int n = vgd_n;
  const int nwall = H.nwall, nexact = H.nexact, nmap = H.nmap;
  if (vgd_gband && nexact != nmap) {
    VgdFatal("vet_gd_band_global: a wall direction map is not the identity (the global "
             "band copies do not re-index directions)");
  }
  if (H.nphx > 0 && !vgd_wsym) {
    // any set: an interpolated / nearest phi-wall map breaks the set's phi mirror
    // symmetry (a net F_phi), unless vet_gd_wall_sym pairs the maps
    VgdFatal("phi-periodic wedge width W must be a multiple of the ray set's azimuthal "
             "period for exact wall maps (" + std::to_string(H.nphx) + " of " +
             std::to_string(nmap) + " phi-wall direction maps inexact), or name "
             "rad_vet/vet_gd_wall_sym = true (with vet_gd_wall_interp = true) for "
             "mirror-paired phi walls");
  }
  if (H.nthx > 0 && !vgd_wst &&
      pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic) {
    // a theta-periodic mesh: its clean theta seam face is one-frame (a net F_theta); a
    // reflecting theta wall has no seam face
    VgdFatal("theta-periodic wedge: the theta wall maps are inexact (" +
             std::to_string(H.nthx) + " of " + std::to_string(nmap) + " direction maps); "
             "name rad_vet/vet_gd_theta_seam_sym = true (with vet_gd_wall_interp = true) "
             "for mirror-paired theta seams");
  }
  const bool first = (vgd_alpha < 0.0);
  // the compact global band's kept lists belong to a set: a new set remakes them
  if (alpha != vgd_alpha) {
    vgd_gbcf.clear();
    vgd_gbcb.clear();
    vgd_gbrec.clear();
  }
  vgd_alpha = alpha;
  g_vgd_tab_alpha = alpha;
  if (first && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> vet_gd: global-direction SC closure on the "
              << (vgd_csg ? "cubed sphere" : "sp wedge") << ", "
              << ((vgd_gl_nmu > 0) ? ("product GL" + std::to_string(vgd_gl_nmu) + "x"
                                      + std::to_string(vgd_gl_nphi)
                                      + (vgd_gl_stag ? " staggered" : ""))
                  : ((vgd_ls > 0) ? ("level-symmetric LQ" + std::to_string(vgd_ls))
                  : (vgd_ls < 0) ? ("geodesic level " + std::to_string(-vgd_ls))
                                  : ("HEALPix nside " + std::to_string(vgd_nside))))
              << " (" << n << " directions); "
              << (vgd_rot3 ? "3-D rotation (quasi-random axis; angle = block key)"
                           : "z-rotation") << " every "
              << ((vgd_rot_dt > 0.0) ? (std::to_string(vgd_rot_dt) + " in time")
                                     : (std::to_string(vgd_rot) + " cycle(s)"))
              << " (0 = fixed), first angle " << alpha << ";"
              << " rank 0 wall ghost columns " << nwall << ", direction maps exact "
              << nexact << " of " << nmap << " (the rest: nearest direction)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdXfTables
//! \brief vet_gd_rotate_xfade: the tables of alpha into vgd_dir/...; both cross-fade
//! sets stay on the device (current and vgd_tabd) and are swapped; another angle swaps
//! the current set out and builds the new one (VetGdTables)

void RadiationVet::VgdXfTables(const Real alpha) {
  if (alpha == vgd_alpha && g_vgd_tab_alpha == alpha) {return;}
  auto swp = [this]() {
    std::swap(vgd_dir, vgd_tabd.dir);
    std::swap(vgd_wall, vgd_tabd.wall);
    std::swap(vgd_map, vgd_tabd.map);
    std::swap(vgd_mr, vgd_tabd.mr);
    if (vgd_wint) {
      std::swap(vgd_m3, vgd_tabd.m3);
      std::swap(vgd_w3, vgd_tabd.w3);
    }
    // the exact halo lists (and the global band's kept lists) go with their set, so
    // the two sets do not rebuild each other's lists (bitwise). Not on the cubed sphere
    // (own per-set slots)
    if (!vgd_csg) {std::swap(vgd_hxe, vgd_hxe_alt);}
    VgdGbSwap();
    // the set's intensities and shell lists go with it (keyed caches)
    if (vgd_xff) {
      std::swap(vgd_i, vgd_i_xb);
      for (int p = 0; p < 2; ++p) {
        std::swap(vgd_shl[p], vgd_shl_xb[p]);
        std::swap(vgd_shn[p], vgd_shn_xb[p]);
      }
      std::swap(vgd_shl_alpha, vgd_shl_alpha_xb);
    }
  };
  Kokkos::fence();
  if (alpha == vgd_xf_alt && vgd_tabd.dir.extent(0) == vgd_dir.extent(0)) {
    swp();
    vgd_xf_alt = g_vgd_tab_alpha;
    vgd_alpha = alpha;
    g_vgd_tab_alpha = alpha;
    return;
  }
  VgdTabAlloc(true);
  swp();
  vgd_xf_alt = g_vgd_tab_alpha;
  VetGdTables(alpha);
}

namespace {
void VgdPackEvtRecord(const DevExeSpace &ex);   // defined with VetGdHlBegin
void VgdPackEvtWait(const DevExeSpace &ex);
}  // namespace

// second, high-priority stream of the cubed sphere's early send (VgdEs2Make)
namespace {
#if defined(KOKKOS_ENABLE_HIP)
DevExeSpace VgdPrioSpace() {
  int lo = 0, hi = 0;
  (void)hipDeviceGetStreamPriorityRange(&lo, &hi);
  hipStream_t s;
  (void)hipStreamCreateWithPriority(&s, hipStreamNonBlocking, hi);
  return DevExeSpace(s, Kokkos::Impl::ManageStream::yes);
}
#elif defined(KOKKOS_ENABLE_CUDA)
DevExeSpace VgdPrioSpace() {
  int lo = 0, hi = 0;
  (void)cudaDeviceGetStreamPriorityRange(&lo, &hi);
  cudaStream_t s;
  (void)cudaStreamCreateWithPriority(&s, cudaStreamNonBlocking, hi);
  return DevExeSpace(s, Kokkos::Impl::ManageStream::yes);
}
#else
DevExeSpace VgdPrioSpace() {return DevExeSpace();}
#endif
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdXffSwap
//! \brief swap the cross-fade's whole set context (tables, angles, exact lists,
//! intensities, shell lists) without a fence (handles only); between the two halo
//! exchanges of a fused shell when the merged message does not apply

void RadiationVet::VgdXffSwap() {
  std::swap(vgd_dir, vgd_tabd.dir);
  std::swap(vgd_wall, vgd_tabd.wall);
  std::swap(vgd_map, vgd_tabd.map);
  std::swap(vgd_mr, vgd_tabd.mr);
  if (vgd_wint) {
    std::swap(vgd_m3, vgd_tabd.m3);
    std::swap(vgd_w3, vgd_tabd.w3);
  }
  if (!vgd_csg) {std::swap(vgd_hxe, vgd_hxe_alt);}
  VgdGbSwap();
  std::swap(vgd_i, vgd_i_xb);
  for (int p = 0; p < 2; ++p) {
    std::swap(vgd_shl[p], vgd_shl_xb[p]);
    std::swap(vgd_shn[p], vgd_shn_xb[p]);
  }
  std::swap(vgd_shl_alpha, vgd_shl_alpha_xb);
  std::swap(vgd_alpha, vgd_xf_alt);
  g_vgd_tab_alpha = vgd_alpha;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdEs2Make
//! \brief the early send's second, high-priority stream vgd_es2 (made once)

void RadiationVet::VgdEs2Make() {
  if (!vgd_es2_made) {
    vgd_es2 = VgdPrioSpace();
    vgd_es2_made = true;
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdXffOk
//! \brief does the fused two-set sweep cover this configuration (same answer on every
//! rank)? Otherwise the two sweeps run

bool RadiationVet::VetGdXffOk() {
  if (!vgd_xff || !vgd_xf) {return false;}
  if (vray_icut < vgd_rag_cut) {return false;}   // (a cut below the shared slabs)
  if (vgd_csg) {
    // the cubed-sphere sweep (VetGdCsSweep): pass lists, no early send
    return !vgd_gband && (vgd_i_xb.dummy == vgd_i.dummy);
  }
  if (vgd_gband && !vgd_raw) {
    return false;
  }
  // IPC halo: one set per exchange (raw band: VetGdGbHaloX; CUDA lean: VetGdHlXff's
  // IPC form, no early send)
#if !defined(KOKKOS_ENABLE_CUDA)
  if (vgd_ipc_state == 1 && !vgd_raw) {return false;}
#endif
  if (vgd_i_xb.dummy != vgd_i.dummy) {return false;}
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHlXff
//! \brief exchange of shell i (pass vgd_hinw, depth ws) for both fused sets in one
//! message per partner (piece p = [A_p, B_p]); bitwise the two exact-list exchanges.
//! false (nothing done) unless both sets hold exact records of this exchange.

bool RadiationVet::VetGdHlXff(const int i, const int ws) {
#if MPI_PARALLEL_ENABLED
  if (!vgd_hmpi || vgd_hlp.on) {return false;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2_ = nx2 + 2*w, c3_ = nx3 + 2*w, nv_ = vgd_n;
  if (i < 0 || i >= n1c ||
      static_cast<int64_t>(nmb1 + 1)*nv_*c3_*c2_ >=
      static_cast<int64_t>(std::numeric_limits<int>::max())) {
    return false;
  }
  const int px = vgd_hinw ? 0 : 1;
  const int ex = px*n1c + i;
  if (ex >= static_cast<int>(vgd_hxe.size()) ||
      ex >= static_cast<int>(vgd_hxe_alt.size())) {
    return false;
  }
  const VgdHlEntry &XA = vgd_hxe[ex];
  const VgdHlEntry &XB = vgd_hxe_alt[ex];
  auto okx = [&](const VgdHlEntry &X, const Real al) {
    return X.gen != -1 && X.alpha == al && VetGdHxCut(ex, X) && X.ws == ws && X.exact;
  };
  if (!okx(XA, vgd_alpha) || !okx(XB, vgd_xf_alt)) {return false;}
  const int np = static_cast<int>(vgd_prk.size());
  if (np > VGD_IPC_NS || static_cast<int>(XA.sp.size()) != np + 1 ||
      static_cast<int>(XB.sp.size()) != np + 1 ||
      static_cast<int>(XA.rp.size()) != np + 1 ||
      static_cast<int>(XB.rp.size()) != np + 1) {
    return false;
  }
  // merged layouts: piece p at so[p] (send) / ro[p] (receive), A first
  std::vector<int> so(np + 1), ro(np + 1);
  VgdIpcIdx sa{}, sb{}, ra{}, rb{}, sm{}, rm{};
  for (int p = 0; p <= np; ++p) {
    so[p] = XA.sp[p] + XB.sp[p];
    ro[p] = XA.rp[p] + XB.rp[p];
    sa.v[p] = XA.sp[p];
    sb.v[p] = XB.sp[p];
    ra.v[p] = XA.rp[p];
    rb.v[p] = XB.rp[p];
    sm.v[p] = so[p];
    rm.v[p] = ro[p];
  }
  if (static_cast<int>(vgd_xff_sb.extent(0)) < std::max(so[np], 1)) {
    Kokkos::realloc(vgd_xff_sb, std::max(so[np], 1) + so[np]/8);
  }
  if (static_cast<int>(vgd_xff_rb.extent(0)) < std::max(ro[np], 1)) {
    Kokkos::realloc(vgd_xff_rb, std::max(ro[np], 1) + ro[np]/8);
  }
  auto sbuf_ = vgd_xff_sb;
  auto rbuf_ = vgd_xff_rb;
  auto aA_ = vgd_i;
  auto aB_ = vgd_i_xb;
  const int nsA = XA.ns, nsB = XB.ns;
  const int npk = np;
#if defined(KOKKOS_ENABLE_CUDA)
  // IPC form (CUDA lean, no early send): one kernel packs both sets' kept values of
  // partner p into its receive buffer, the last team sets the partners' flags, every
  // team waits for mine and scatters both sets (bitwise vet_vgd_xff_scatter)
  if (vgd_ipc_state == 1 &&
      vgd_rg_pend.empty() && static_cast<int>(XA.dsto.size()) == np &&
      static_cast<int>(XB.dsto.size()) == np) {
    if (!vgd_hlp.cnt.is_allocated()) {
      vgd_hlp.cnt = Kokkos::View<int, DevMemSpace>("vet_vgd_ipc_cnt");
    }
    const long long seq = ++vgd_ipc_seq, par = seq & 1;
    auto cnt_ = vgd_hlp.cnt;
    long long *sig = vgd_ipc_sig;
    int *err = vgd_ipc_err;
    auto psig = vgd_ipc_psig.d_view;
    auto pslot = vgd_ipc_pslot.d_view;
    auto prb = vgd_ipc_prb.d_view;
    auto phalf = vgd_ipc_phalf.d_view;
    auto pcs = vgd_ipc_pcs.d_view;
    const Real *rb_ = vgd_ipc_rbuf + par*vgd_ipc_half;
    const size_t cs_ = vgd_ipc_cs;
    VgdIpcIdx da{}, db{};
    for (int p = 0; p < np; ++p) {da.v[p] = XA.dsto[p]; db.v[p] = XB.dsto[p];}
    auto aA_ = vgd_i;
    auto aB_ = vgd_i_xb;
    auto slA_ = XA.sl;
    auto slB_ = XB.sl;
    const int nsA = XA.ns, nsB = XB.ns;
    const int nrA = XA.nr, nA = XA.nr + XA.nw, nrB = XB.nr, nB = XB.nr + XB.nw;
    auto dlA_ = XA.dl;
    auto srA_ = XA.sr;
    auto dwA_ = XA.dw;
    auto swA_ = XA.sw;
    auto dlB_ = XB.dl;
    auto srB_ = XB.sr;
    auto dwB_ = XB.dw;
    auto swB_ = XB.sw;
    auto w3A_ = vgd_w3;
    auto w3B_ = vgd_tabd.w3;
    constexpr int kT = 128;
    const int nt = std::max(nsA + nsB, nA + nB);
    const int nl = std::min((nt + kT - 1)/kT, 216);
    Kokkos::parallel_for("vet_vgd_ipc_xff", Kokkos::Experimental::require(
        Kokkos::TeamPolicy<>(vgd_cur, std::max(nl, 1), kT),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight),
    KOKKOS_LAMBDA(const TeamMember_t &tm) {
      const int e0 = tm.league_rank()*kT + tm.team_rank(), de = tm.league_size()*kT;
      for (int c = e0; c < nsA + nsB; c += de) {
        const bool b = (c >= nsA);
        const int cc = b ? (c - nsA) : c;
        const VgdIpcIdx &sx = b ? sb : sa;
        int p = 0;
        while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
        int ad = b ? slB_(cc) : slA_(cc);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + (b ? pcs(p) : 0) +
                    (b ? db.v[p] : da.v[p]) + (cc - sx.v[p]);
        *dst = b ? aB_(m,v,kb,jb,i) : aA_(m,v,kb,jb,i);
      }
#if defined(__CUDA_ARCH__)
      __threadfence_system();
#endif
      tm.team_barrier();
      if (tm.team_rank() == 0) {
        if (Kokkos::atomic_fetch_add(&cnt_(), 1) == tm.league_size() - 1) {
          cnt_() = 0;
          VgdIpcSetN(psig.data(), pslot.data(), npk, seq);
        }
        VgdIpcWaitN(sig, npk, seq, err);
      }
      tm.team_barrier();
      for (int c0 = e0; c0 < nA + nB; c0 += de) {
        const bool b = (c0 >= nA);
        const int c = b ? (c0 - nA) : c0;
        const int nr = b ? nrB : nrA;
        const Real *rx = rb_ + (b ? cs_ : 0);
        int ad = (c < nr) ? (b ? dlB_(c) : dlA_(c)) : (b ? dwB_(c - nr) : dwA_(c - nr));
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        Real val = 0.0;
        if (c < nr) {
          const int s = b ? srB_(c) : srA_(c);
          val = (s >= 0) ? rx[s] : 0.0;
        } else {
          const int cw = c - nr;
          for (int u = 0; u < 3; ++u) {
            const int s = b ? swB_(3*cw + u) : swA_(3*cw + u);
            const Real wq = b ? w3B_(m,kb,jb,v,u) : w3A_(m,kb,jb,v,u);
            val += wq*((s >= 0) ? rx[s] : 0.0);
          }
        }
        if (b) {
          aB_(m,v,kb,jb,i) = val;
        } else {
          aA_(m,v,kb,jb,i) = val;
        }
      }
    });
    vgd_ipc_n += 1.0;
    vgd_nexch += 2.0;
    vgd_hl_nuse += 2.0;
    vgd_xff_nmsg += 1.0;
    vgd_rg_lastx = true;
    vgd_rg_lastx_e = ex;
    return true;
  }
#endif
  if (nsA + nsB > 0) {
    auto slA_ = XA.sl;
    auto slB_ = XB.sl;
    Kokkos::parallel_for("vet_vgd_xff_pack", VgdLight(vgd_cur, 0, nsA + nsB),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nsA);
      const int cc = b ? (c - nsA) : c;
      const VgdIpcIdx &sx = b ? sb : sa;
      int p = 0;
      while (p + 1 < npk && sx.v[p+1] <= cc) {++p;}
      int ad = b ? slB_(cc) : slA_(cc);
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      const int pos = sm.v[p] + (b ? (sa.v[p+1] - sa.v[p]) : 0) + (cc - sx.v[p]);
      sbuf_(pos) = b ? aB_(m,v,kb,jb,i) : aA_(m,v,kb,jb,i);
    });
  }
  // band-ring copies of other shells go between pack and scatter (as VetGdHlBegin)
  VgdPackEvtRecord(vgd_cur);
  VetGdRingFlush();
  VgdPackEvtWait(vgd_cur);
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_prk[p];
    req.emplace_back();
    MPI_Irecv(rbuf_.data() + ro[p], ro[p+1] - ro[p], MPI_ATHENA_REAL, rk, 7005, vgd_comm,
              &req.back());
    req.emplace_back();
    MPI_Isend(sbuf_.data() + so[p], so[p+1] - so[p], MPI_ATHENA_REAL, rk, 7005, vgd_comm,
              &req.back());
  }
  vgd_tpost += tq.seconds();
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 2.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != ro[p+1] - ro[p]) {
      VgdFatal("vet_gd_rotate_xfade_fuse: sender and receiver lists disagree");
    }
    vgd_nbyte += static_cast<Real>(so[p+1] - so[p])*sizeof(Real);
  }
  // scatter both sets (A, then B), as VetGdHlEnd per set
  const int nrA = XA.nr, nwA = XA.nw, nrB = XB.nr, nwB = XB.nw;
  const int nA = nrA + nwA, nB = nrB + nwB;
  if (nA + nB > 0) {
    auto dlA_ = XA.dl;
    auto srA_ = XA.sr;
    auto dwA_ = XA.dw;
    auto swA_ = XA.sw;
    auto dlB_ = XB.dl;
    auto srB_ = XB.sr;
    auto dwB_ = XB.dw;
    auto swB_ = XB.sw;
    auto w3A_ = vgd_w3;
    auto w3B_ = vgd_tabd.w3;
    Kokkos::parallel_for("vet_vgd_xff_scatter", VgdLight(vgd_cur, 0, nA + nB),
    KOKKOS_LAMBDA(const int c0) {
      const bool b = (c0 >= nA);
      const int c = b ? (c0 - nA) : c0;
      const int nr = b ? nrB : nrA;
      const VgdIpcIdx &rx = b ? rb : ra;
      // merged receive position of compact source s of this set
      auto at = [&](const int s) -> int {
        int p = 0;
        while (p + 1 < npk && rx.v[p+1] <= s) {++p;}
        return rm.v[p] + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (s - rx.v[p]);
      };
      int ad = (c < nr) ? (b ? dlB_(c) : dlA_(c)) : (b ? dwB_(c - nr) : dwA_(c - nr));
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      if (c < nr) {
        const int s = b ? srB_(c) : srA_(c);
        const Real x = (s >= 0) ? rbuf_(at(s)) : 0.0;
        if (b) {
          aB_(m,v,kb,jb,i) = x;
        } else {
          aA_(m,v,kb,jb,i) = x;
        }
        return;
      }
      const int cw = c - nr;
      Real val = 0.0;
      for (int u = 0; u < 3; ++u) {
        const int s = b ? swB_(3*cw + u) : swA_(3*cw + u);
        const Real wq = b ? w3B_(m,kb,jb,v,u) : w3A_(m,kb,jb,v,u);
        val += wq*((s >= 0) ? rbuf_(at(s)) : 0.0);
      }
      if (b) {
        aB_(m,v,kb,jb,i) = val;
      } else {
        aA_(m,v,kb,jb,i) = val;
      }
    });
  }
  vgd_hl_nuse += 2.0;
  vgd_xff_nmsg += 1.0;
  vgd_rg_lastx = true;   // the band of this exchange came from both sets' lists
  vgd_rg_lastx_e = ex;
  return true;
#else
  (void) i;
  (void) ws;
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHaloInit
//! \brief the 8 lateral neighbour slots of every local MeshBlock (same x1 level; theta
//! and phi periodic by wrapping the logical location) and the message buffers

void RadiationVet::VetGdHaloInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  int nb2 = 1, nb3 = 1;
  for (int g = 0; g < pm->nmb_total; ++g) {
    nb2 = std::max(nb2, pm->lloc_eachmb[g].lx2 + 1);
    nb3 = std::max(nb3, pm->lloc_eachmb[g].lx3 + 1);
  }
  std::vector<int> gof(vgd_nb1*nb2*nb3, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx1*nb3 + l.lx3)*nb2 + l.lx2] = g;
  }
  for (int t = 0; t < vgd_nb1*nb2*nb3; ++t) {
    if (gof[t] < 0) {VgdFatal("incomplete MeshBlock grid (SMR/AMR not supported)");}
  }
  Kokkos::realloc(vgd_hloc, 8*nmb);
  auto hl_h = Kokkos::create_mirror_view(vgd_hloc);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  vgd_hmpi = false;
  vgd_nb2 = nb2;
  vgd_nb3 = nb3;
  Kokkos::realloc(vgd_lxy, nmb, 2);
  {
    auto lx_h = Kokkos::create_mirror_view(vgd_lxy);
    for (int m = 0; m < nmb; ++m) {
      lx_h(m,0) = pm->lloc_eachmb[g0 + m].lx2;
      lx_h(m,1) = pm->lloc_eachmb[g0 + m].lx3;
    }
    Kokkos::deep_copy(vgd_lxy, lx_h);
  }
  // full sphere: a slot beyond a pole takes the pole partner (lx2, lx3 + dk + nb3/2),
  // rows reversed (vgd_hpol 1 north, 2 south, else 0)
  if (vgd_pol && (nb3 % 2 != 0)) {
    VgdFatal("the full sphere needs an EVEN number of MeshBlocks in phi (pole partners)");
  }
  Kokkos::realloc(vgd_hpol, 8*nmb);
  auto hp_h = Kokkos::create_mirror_view(vgd_hpol);
  for (int m = 0; m < nmb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    for (int o = 0; o < 8; ++o) {
      const int oo = (o < 4) ? o : (o + 1);       // skip the centre of the 3x3
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int pol = VgdPoleSlot(o, l.lx2);
      hp_h(8*m + o) = pol;
      const int l2 = (pol != 0) ? l.lx2 : (l.lx2 + dj + nb2) % nb2;
      const int l3 = (pol != 0) ? (l.lx3 + dk + nb3/2 + nb3) % nb3
                                : (l.lx3 + dk + nb3) % nb3;
      const int g = gof[(l.lx1*nb3 + l3)*nb2 + l2];
      const int rk = pm->rank_eachmb[g];
      if (rk == me) {
        hl_h(8*m + o) = g - g0;
      } else {
        hl_h(8*m + o) = -1;
        vgd_hrank[8*m + o] = rk;
        vgd_hlid[8*m + o] = g - pm->gids_eachrank[rk];
        vgd_hmpi = true;
      }
    }
  }
  Kokkos::deep_copy(vgd_hloc, hl_h);
  Kokkos::deep_copy(vgd_hpol, hp_h);
  // the largest region (edge slot: nx x W) times the most variables per message
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nreg = std::max(indcs.nx2, indcs.nx3)*vgd_w;
  // one shell of all directions (the multi-shell ln chi, ln S exchange is chunked to
  // fit, VetGdHaloCs)
  (void) c1;
  vgd_maxcnt = nreg*std::max(vgd_n, 2);
  if (vgd_hmpi) {
    Kokkos::realloc(vgd_sbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    Kokkos::realloc(vgd_rbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    // message plan: pieces sorted by (rank, receiver block, receiver slot) on both
    // sides; piece size kn*jn with kn, jn in {nx, ws}
    std::vector<std::array<int, 4>> ps, pr;
    for (int m = 0; m < nmb; ++m) {
      for (int o = 0; o < 8; ++o) {
        if (hl_h(8*m + o) >= 0) {continue;}
        // the receiver's slot pointing back here: the opposite one, or across a pole
        // the one with the same theta side and the opposite phi side
        int orc = 7 - o;
        if (hp_h(8*m + o) != 0) {
          const int oo = (o < 4) ? o : (o + 1);
          const int dk = oo/3 - 1, dj = oo%3 - 1;
          const int ob = 3*(-dk + 1) + (dj + 1);
          orc = (ob < 4) ? ob : (ob - 1);
        }
        ps.push_back({vgd_hrank[8*m + o], vgd_hlid[8*m + o], orc, 8*m + o});
        pr.push_back({vgd_hrank[8*m + o], m, o, 8*m + o});
      }
    }
    std::sort(ps.begin(), ps.end());
    std::sort(pr.begin(), pr.end());
    vgd_prk.clear();
    for (const auto &x : ps) {
      if (vgd_prk.empty() || vgd_prk.back() != x[0]) {vgd_prk.push_back(x[0]);}
    }
    const int np = static_cast<int>(vgd_prk.size());
    Kokkos::realloc(vgd_soff, vgd_w + 1, 8*nmb);
    Kokkos::realloc(vgd_roff, vgd_w + 1, 8*nmb);
    auto so_h = Kokkos::create_mirror_view(vgd_soff);
    auto ro_h = Kokkos::create_mirror_view(vgd_roff);
    vgd_pdsp.assign(vgd_w + 1, std::vector<int>(np, 0));
    vgd_pscnt = vgd_pdsp;
    vgd_prdsp = vgd_pdsp;
    vgd_prcnt = vgd_pdsp;
    auto sz = [&](const int o, const int ws) {
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      return ((dk == 0) ? indcs.nx3 : ws)*((dj == 0) ? indcs.nx2 : ws);
    };
    for (int ws = 0; ws <= vgd_w; ++ws) {
      for (int t = 0; t < 8*nmb; ++t) {so_h(ws,t) = -1; ro_h(ws,t) = -1;}
      for (int side = 0; side < 2; ++side) {
        const auto &L = (side == 0) ? ps : pr;
        int off = 0, p = -1, prev = -1;
        for (const auto &x : L) {
          if (x[0] != prev) {
            ++p;
            prev = x[0];
            ((side == 0) ? vgd_pdsp : vgd_prdsp)[ws][p] = off;
          }
          const int o = x[3] % 8;
          ((side == 0) ? so_h : ro_h)(ws, x[3]) = off;
          // piece size: the receiver's slot 7 - o has the shape of my slot o
          off += sz(o, ws);
          ((side == 0) ? vgd_pscnt : vgd_prcnt)[ws][p] += sz(o, ws);
        }
      }
    }
    Kokkos::deep_copy(vgd_soff, so_h);
    Kokkos::deep_copy(vgd_roff, ro_h);
    if (vgd_hcomp > 0) {
      // partner boundaries of both dense layouts per depth ws (in units of nv*ni)
      const int npp = static_cast<int>(vgd_prk.size());
      Kokkos::realloc(vgd_pbd, vgd_w + 1, 2*(npp + 1));
      auto pb_h = Kokkos::create_mirror_view(vgd_pbd);
      for (int ws = 0; ws <= vgd_w; ++ws) {
        for (int p = 0; p < npp; ++p) {
          pb_h(ws,p) = vgd_pdsp[ws][p];
          pb_h(ws,npp + 1 + p) = vgd_prdsp[ws][p];
        }
        pb_h(ws,npp) = (npp > 0) ? (vgd_pdsp[ws][npp-1] + vgd_pscnt[ws][npp-1]) : 0;
        pb_h(ws,2*npp + 1) = (npp > 0) ? (vgd_prdsp[ws][npp-1] + vgd_prcnt[ws][npp-1])
                                       : 0;
      }
      Kokkos::deep_copy(vgd_pbd, pb_h);
      for (int sl = 0; sl < 2; ++sl) {
        Kokkos::realloc(vgd_hpb[sl], 2*(npp + 1));
      }
      auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       pmy_pack->pcoord->x1v);
      auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       pmy_pack->pcoord->xx1f);
      vgd_r1v.assign(x1v_h.extent(1), 0.0);
      vgd_r1f.assign(x1f_h.extent(1), 0.0);
      for (size_t q = 0; q < vgd_r1v.size(); ++q) {vgd_r1v[q] = x1v_h(0,q);}
      for (size_t q = 0; q < vgd_r1f.size(); ++q) {vgd_r1f[q] = x1f_h(0,q);}
      // exact per-shell sizes (ni = 1, nv = vgd_n) at the deepest band: sums of the
      // remote pieces (only the dense sbuf / rbuf carry the multi-shell exchange)
      const size_t ns = static_cast<size_t>(pb_h(vgd_w,npp))*vgd_n;
      const size_t nr = static_cast<size_t>(pb_h(vgd_w,2*npp + 1))*vgd_n;
      Kokkos::realloc(vgd_csb, std::max<size_t>(ns, 1));
      Kokkos::realloc(vgd_crb, std::max<size_t>(nr, 1));
      // mask slots: this shell's and the next one's
      for (int sl = 0; sl < 2; ++sl) {
        Kokkos::realloc(vgd_hfs[sl], ns + 1);
        Kokkos::realloc(vgd_hfr[sl], nr + 1);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHalo
//! \brief fill the lateral band (vgd_w cells) of a(m, v, k, j, i), v in [0, nv), i in
//! [i0, i1], from the 8 lateral neighbours' interiors: local ones copied directly (read
//! interiors, write bands: no race), remote ones by one message per (block, slot).
//! vgd_raw: global-cell lists, unmapped; intensities (mapd) only from remote owners

template <class V>
void RadiationVet::VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                             const bool mapd) {
  // vet_gd_band_global: every band position within depth ws from its owner's interior
  // (global-cell lists, sph2.inc); no direction re-indexing at the full sphere's walls
  if (vgd_gband || vgd_raw) {
    // compact per-shell exchange of the global band
    if (mapd && vgd_gbc && i1 == i0 && nv == vgd_n) {
      VetGdGbHaloC(a, i0, ws);
      return;
    }
    VetGdCsHalo(a, nv, i0, i1, !(vgd_raw && mapd), ws);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int ni = i1 - i0 + 1;
  const int mx = std::max(nx2, nx3)*ws*nv*ni;
  auto hl_ = vgd_hloc;
  auto hp_ = vgd_hpol;
  auto a_ = a;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const bool wint = vgd_wint;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  // region of slot o along one direction d = -1, 0, +1 with depth ws: (dest start, len,
  // source start) = (w - ws, ws, w + nx - ws), (w, nx, w), (w + nx, ws, w); wall
  // ghosts (mapd) take direction mp(v) of the source (new(d) = old(map(d)))
  const bool mpi = vgd_hmpi;
  if (mapd && !mpi) {return;}
  auto sb_ = vgd_sbuf;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  const int nvi = nv*ni;
  const bool hcomp = mapd && (vgd_hcomp > 0) && (ni == 1);
  // true when the compact halo's list path wrote the band itself
  bool done = false;
  if (hcomp) {done = VetGdHaloCompact(a, nv, i0, ws);}
  if (!hcomp) {
  par_for("vet_vgd_halo", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    const int cnt = nv*kn*jn*ni;
    if (t >= cnt) {return;}
    const int v = t/(kn*jn*ni);
    const int r1 = t - v*kn*jn*ni;
    const int kk = r1/(jn*ni);
    const int r2 = r1 - kk*jn*ni;
    const int jj = r2/ni;
    const int i = i0 + (r2 - jj*ni);
    const int nl = hl_(8*m + o);
    // a slot beyond a pole (full sphere): the partner's rows reversed
    const int pl = hp_(8*m + o);
    if (nl >= 0 && mapd) {return;}    // the sweep reads on-rank neighbours directly
    if (nl >= 0) {
      const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
      const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
      const int js_ = (dj < 0) ? (w + nx2 - ws) : w;
      // nx3 = 1: every phi band column is a copy of the one column
      const int ks_ = (nx3 == 1) ? w : (((dk < 0) ? (w + nx3 - ws) : w) + kk);
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const int jsrc = (pl != 0) ? VgdPoleRow(jd + jj, w, nx2) : (js_ + jj);
      a_(m,v,kd+kk,jd+jj,i) = a_(nl,vs,ks_,jsrc,i);
    } else if (mpi) {
      // my interior that the remote neighbour needs in its opposite slot (across a
      // pole: my rows reversed)
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (nx3 == 1) ? w : (((dk > 0) ? (w + nx3 - ws) : w) + kk);
      const int jsrc = (pl == 0) ? (js2 + jj)
                       : ((pl == 1) ? (w + ws - 1 - jj) : (w + nx2 - 1 - jj));
      sb_(static_cast<size_t>(so_(ws,8*m + o))*nvi + t) = a_(m,v,ks2,jsrc,i);
    }
  });
  }  // !hcomp
#if MPI_PARALLEL_ENABLED
  if (mpi && !hcomp) {
    vgd_cur.fence();
    Kokkos::Timer tq;
    std::vector<MPI_Request> req;
    auto rb_ = vgd_rbuf;
    // one message per partner rank and direction
    for (size_t p = 0; p < vgd_prk.size(); ++p) {
      const int rk = vgd_prk[p];
      req.emplace_back();
      MPI_Irecv(rb_.data() + static_cast<size_t>(vgd_prdsp[ws][p])*nvi,
                vgd_prcnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, vgd_comm,
                &req.back());
      req.emplace_back();
      MPI_Isend(sb_.data() + static_cast<size_t>(vgd_pdsp[ws][p])*nvi,
                vgd_pscnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, vgd_comm,
                &req.back());
    }
    vgd_tpost += tq.seconds();
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    vgd_tmpi += tq.seconds();
    vgd_nexch += 1.0;
    for (size_t p = 0; p < vgd_prk.size(); ++p) {
      vgd_nbyte += static_cast<Real>(vgd_pscnt[ws][p])*nvi*sizeof(Real);
    }
  }
  if (mpi && !done) {
    auto rb_ = vgd_rbuf;
    par_for("vet_vgd_unpack", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
    KOKKOS_LAMBDA(const int m, const int o, const int t) {
      if (hl_(8*m + o) >= 0) {return;}
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
      const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
      const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
      const int cnt = nv*kn*jn*ni;
      if (t >= cnt) {return;}
      const int v = t/(kn*jn*ni);
      const int r1 = t - v*kn*jn*ni;
      const int kk = r1/(jn*ni);
      const int r2 = r1 - kk*jn*ni;
      const int jj = r2/ni;
      const int ii = r2 - jj*ni;
      const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nvi;
      if (mapd && wint && wl_(m,kd+kk,jd+jj) != 0) {
        Real val = 0.0;
        for (int q = 0; q < 3; ++q) {
          const int vq = m3_(m,kd+kk,jd+jj,v,q);
          const size_t at = rb0 + ((vq*kn + kk)*jn + jj)*ni + ii;
          val += w3_(m,kd+kk,jd+jj,v,q)*rb_(at);
        }
        a_(m,v,kd+kk,jd+jj,i0+ii) = val;
        return;
      }
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + ((vs*kn + kk)*jn + jj)*ni + ii;
      a_(m,v,kd+kk,jd+jj,i0+ii) = rb_(at);
    });
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHaloCs
//! \brief VetGdHalo of ln chi, ln S over shells [i0, i1] at full band depth, in chunks
//! that fit vgd_maxcnt per slot (bitwise the one-call exchange)

void RadiationVet::VetGdHaloCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int per = std::max(indcs.nx2, indcs.nx3)*vgd_w*2;   // one shell, one slot
  const int nch = std::max(1, (per > 0) ? vgd_maxcnt/per : 1);
  for (int c0 = i0; c0 <= i1; c0 += nch) {
    VetGdHalo(a, 2, c0, std::min(c0 + nch - 1, i1), vgd_w, false);
  }
  if (vgd_thmir) {VetGdThMirCs(a, i0, i1);}
  if (vgd_thvac > 0) {VetGdThVacCs(a, i0, i1);}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdThVacCs
//! \brief vet_gd_theta_wall = vacuum / extrap: ln chi, ln S of the theta ghost rows
//! beyond a physical wall (flags 3, 4) = the edge row's (zero gradient)

void RadiationVet::VetGdThVacCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int wb = vgd_w, nx2 = indcs.nx2;
  const int c3 = indcs.nx3 + 2*vgd_w;
  const int c2 = nx2 + 2*wb;
  const Real t0 = pmy_pack->pmesh->mesh_size.x2min;
  const Real t1 = pmy_pack->pmesh->mesh_size.x2max;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const VgdLG lg_ = vgd_lg;
  auto a_ = a;
  auto wl_ = vgd_wall;
  par_for("vet_vgd_thvac_cs", vgd_cur, 0, nmb1, 0, 1, 0, c3 - 1, 0, 2*wb - 1, i0, i1,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int s, const int i) {
    const int g = (s < wb) ? s : (s - wb);
    const int j = (s < wb) ? (wb - 1 - g) : (wb + nx2 + g);
    if (wl_(m,k,j) < 3) {return;}
    // the band row next to the wall (the block's edge row with one theta MeshBlock)
    int je = (s < wb) ? static_cast<int>(lg_.Jw(mbsize, m, wb, t0))
                      : (static_cast<int>(lg_.Jw(mbsize, m, wb, t1)) - 1);
    je = (je < 0) ? 0 : ((je > c2 - 1) ? (c2 - 1) : je);
    a_(m,c,k,j,i) = a_(m,c,k,je,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdThMirCs
//! \brief vet_gd_theta_wall = mirror: ln chi, ln S of the theta ghost rows = their
//! mirror rows (after the exchange filled the phi band)

void RadiationVet::VetGdThMirCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int wb = vgd_w, nx2 = indcs.nx2;
  const int c3 = indcs.nx3 + 2*vgd_w;
  auto a_ = a;
  par_for("vet_vgd_thmir_cs", vgd_cur, 0, nmb1, 0, 1, 0, c3 - 1, 0, 2*wb - 1, i0, i1,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int s, const int i) {
    const int g = (s < wb) ? s : (s - wb);
    const int j = (s < wb) ? (wb - 1 - g) : (wb + nx2 + g);
    const int jm = VgdMirRow(g, s >= wb, wb, nx2);
    a_(m,c,k,j,i) = a_(m,c,k,jm,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn KOKKOS_INLINE_FUNCTION bool VgdHaloKeep
//! \brief vet_gd_halo_compact: is direction u at source cell (global gj, gk) at band
//! depth dep needed by the reads of shell s in this pass? Global indices only, so
//! sender and receiver agree bitwise. mode 1: the pass's branch (|n.r| < 1e-12 kept);
//! mode 2: also the ray's lateral reach >= dep - 3 cells. wall: branch only.

namespace {
// the keep test at a source-cell position (theta, phi): shared by the uniform slot halo
// (VgdHaloKeep) and the global band (VetGdGbHaloC, node table)
KOKKOS_INLINE_FUNCTION
bool VgdHaloKeepTP(const Real nx, const Real ny, const Real nz, const Real th,
                   const Real ph, const bool inw, const int mode, const bool wall,
                   const int dep, const Real rs, const Real rm, const Real rp,
                   const bool hasm, const bool hasp, const Real gam, const Real dmin) {
  const Real st = sin(th);
  const Real mu = nx*st*cos(ph) + ny*st*sin(ph) + nz*cos(th);
  if (fabs(mu) < 1.0e-12) {return true;}
  if (inw != (mu < 0.0)) {return false;}
  if (mode < 2 || wall) {return true;}
  // the reader is up to dep cells from the source cell (psi differs by up to that;
  // gam = 2 cell widths)
  const Real psi = acos(fmin(fabs(mu), 1.0)), gw = 0.5*gam*(dep + 2);
  const Real phi_ = fmin(psi + gw, 0.5*M_PI), plo_ = fmax(psi - gw, 0.0);
  const Real ph_ = rs*sin(phi_), pl_ = rs*sin(plo_);
  Real dl = 0.0;
  if (inw) {
    // NORM read by shell s-1 (rm); TURN read by shell s itself (r_lo = rm)
    if (hasm) {
      const Real pp = fmin(ph_, rm);
      dl = fmax(dl, asin(fmin(pp/rm, 1.0)) - asin(fmin(pp/rs, 1.0)));
      if (ph_ > rm) {dl = fmax(dl, 2.0*acos(fmin(fmax(pl_, rm)/rs, 1.0)));}
    }
  } else if (hasp) {
    // NORM read by shell s+1 (rp)
    const Real pp = fmin(ph_, rs);
    dl = fmax(dl, asin(fmin(pp/rs, 1.0)) - asin(fmin(pp/rp, 1.0)));
  }
  // the read point is up to dl from the cell centre: up to ceil(dl/dmin) cells away
  return dep <= static_cast<int>(ceil(dl/dmin)) + 3;
}

KOKKOS_INLINE_FUNCTION
bool VgdHaloKeep(const Real nx, const Real ny, const Real nz, const int gj, const int gk,
                 const Real t0, const Real dth, const Real p0, const Real dph,
                 const bool inw, const int mode, const bool wall, const int dep,
                 const Real rs, const Real rm, const Real rp, const bool hasm,
                 const bool hasp, const Real gam, const Real dmin) {
  const Real th = t0 + (gj + 0.5)*dth, ph = p0 + (gk + 0.5)*dph;
  return VgdHaloKeepTP(nx, ny, nz, th, ph, inw, mode, wall, dep, rs, rm, rp, hasm, hasp,
                       gam, dmin);
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHcPrep
//! \brief vet_gd_halo_compact = 1, 2: the keep masks of (shell, pass, depth) on the
//! dense layouts of both sides (VgdHaloKeep), their scans and partner bounds
void RadiationVet::VetGdHcPrep(const int slot, const int i0, const bool inw0,
                               const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int i = i0;
  const bool inw = inw0;
  const int mode = vgd_hcomp;
  const int nb2 = vgd_nb2, nb3 = vgd_nb3;
  // mirror / extrap theta walls: rows within ws of a wall keep their branch (gj global)
  const bool thm = vgd_thmir || (vgd_thvac == 2);
  const int nx2m = pm->mesh_indcs.nx2;
  const Real t0 = pm->mesh_size.x2min, p0 = pm->mesh_size.x3min;
  const Real dth = (pm->mesh_size.x2max - t0)/pm->mesh_indcs.nx2;
  const Real dph = (pm->mesh_size.x3max - p0)/pm->mesh_indcs.nx3;
  const Real smin = fmin(sin(t0), sin(pm->mesh_size.x2max));
  const Real dmin = fmin(dth, smin*dph);
  const Real gam = 2.0*fmax(dth, dph);
  const int lcut = indcs.is + vgd_scut;
  const Real rs = vgd_r1v[i];
  const bool hasm = (i >= lcut);
  const Real rm = (i > lcut) ? vgd_r1v[i-1] : vgd_r1f[i];
  const bool hasp = (i < indcs.ie);
  const Real rp = hasp ? vgd_r1v[i+1] : rs;
  auto hl_ = vgd_hloc;
  auto lx_ = vgd_lxy;
  auto dir_ = vgd_dir;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  auto fs_ = vgd_hfs[slot];
  auto fr_ = vgd_hfr[slot];
  const int nvi = nv;
  const size_t stot = static_cast<size_t>(vgd_pdsp[ws].empty() ? 0 :
                      (vgd_pdsp[ws].back() + vgd_pscnt[ws].back()))*nvi;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nvi;
  // (1) flags of both dense layouts (sender: my interior; receiver: my band)
  par_for("vet_vgd_hc_flag", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int lx2 = lx_(m,0), lx3 = lx_(m,1);
    const Real nx = dir_(v,0), ny = dir_(v,1), nz = dir_(v,2);
    {
      // send (my slot o, receiver slot direction -dj, -dk)
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
      const int gj = lx2*nx2 + (js2 + jj - w), gk = lx3*nx3 + (ks2 + kk - w);
      const int rdj = -dj, rdk = -dk;
      const int dpj = (rdj < 0) ? (ws - jj) : ((rdj > 0) ? (jj + 1) : 0);
      const int dpk = (rdk < 0) ? (ws - kk) : ((rdk > 0) ? (kk + 1) : 0);
      const int lr2 = (lx2 + dj + nb2) % nb2, lr3 = (lx3 + dk + nb3) % nb3;
      // mirror theta walls: rows within ws of a wall are read by mirrored rays too;
      // keep their whole branch (decided from gj alone)
      const bool wall = (lr2 + rdj < 0) || (lr2 + rdj >= nb2) || (lr3 + rdk < 0)
                        || (lr3 + rdk >= nb3) || (thm && (gj < ws || gj >= nx2m - ws));
      fs_(static_cast<size_t>(so_(ws,8*m + o))*nvi + t) =
        VgdHaloKeep(nx, ny, nz, gj, gk, t0, dth, p0, dph, inw, mode, wall,
                    (dpj > dpk) ? dpj : dpk, rs, rm, rp, hasm, hasp, gam, dmin) ? 1 : 0;
    }
    {
      // receive (my slot o): the source cell in the neighbour, global indices
      const int js_ = (dj < 0) ? (w + nx2 - ws) : w;
      const int ks_ = (dk < 0) ? (w + nx3 - ws) : w;
      const int ln2 = (lx2 + dj + nb2) % nb2, ln3 = (lx3 + dk + nb3) % nb3;
      const int gj = ln2*nx2 + (js_ + jj - w), gk = ln3*nx3 + (ks_ + kk - w);
      const int dpj = (dj < 0) ? (ws - jj) : ((dj > 0) ? (jj + 1) : 0);
      const int dpk = (dk < 0) ? (ws - kk) : ((dk > 0) ? (kk + 1) : 0);
      const bool wall = (lx2 + dj < 0) || (lx2 + dj >= nb2) || (lx3 + dk < 0)
                        || (lx3 + dk >= nb3) || (thm && (gj < ws || gj >= nx2m - ws));
      fr_(static_cast<size_t>(ro_(ws,8*m + o))*nvi + t) =
        VgdHaloKeep(nx, ny, nz, gj, gk, t0, dth, p0, dph, inw, mode, wall,
                    (dpj > dpk) ? dpj : dpk, rs, rm, rp, hasm, hasp, gam, dmin) ? 1 : 0;
    }
  });
  // (2) exclusive scans in place (flag -> position, total at n); flag = pos(q+1) -
  // pos(q); both in one launch (two 32-bit counts in one 64-bit sum)
  VgdScan2(vgd_cur, vgd_hfs[slot], stot, vgd_hfr[slot], rtot, "vet_vgd_hc_scan2",
           vgd_scan_sums);
  // (3) scan values at the partner boundaries, straight into pinned host memory (read
  // after the exchange's fence)
  const int np = static_cast<int>(vgd_prk.size());
  {
    auto pbd_ = vgd_pbd;
    auto pbv_ = vgd_hpb[slot];
    Kokkos::parallel_for("vet_vgd_hc_bnd",
                         Kokkos::RangePolicy<>(vgd_cur, 0, 2*(np + 1)),
    KOKKOS_LAMBDA(const int q) {
      const size_t at = static_cast<size_t>(pbd_(ws,q))*nvi;
      pbv_(q) = (q <= np) ? fs_(at) : fr_(at);
    });
  }
  vgd_htag[slot][0] = i0;
  vgd_htag[slot][1] = inw0 ? 1 : 0;
  vgd_htag[slot][2] = ws;
  vgd_htag[slot][3] = vgd_scut;
  vgd_htag[slot][4] = vgd_hsweep;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHlBuild
//! \brief vet_gd_halo_list_mb: entry e (pass inw, shell i, depth ws) from the compact
//! masks fs / fr and partner bounds pbh:
//!   sl     band address (m, v, k, j) of the c-th sent value,
//!   dl, sr band address of every receiver ghost with a sent source, and its position,
//!   dw, sw the same for vet_gd_wall_interp wall ghosts (3 sources, -1 = unsent).
//! Address = ((m nv + v) c3 + k) c2 + j at depth vgd_w. False when over the budget.

bool RadiationVet::VetGdHlBuild(const int e, const int i, const bool inw, const int ws,
                                const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                                const std::vector<int> &pbh) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  VgdHlEntry &E = vgd_hle[e];
  E.gen = -1;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  if (vgd_hl_fn.extent(0) < rtot + 1) {
    Kokkos::realloc(vgd_hl_fn, rtot + 1);
    Kokkos::realloc(vgd_hl_fw, rtot + 1);
  }
  auto hl_ = vgd_hloc;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  const bool wint = vgd_wint;
  auto fs_ = fs;
  auto fr_ = fr;
  auto fn_ = vgd_hl_fn;
  auto fw_ = vgd_hl_fw;
  // (1) receiver ghosts that the dense unpack fills from at least one sent value
  par_for("vet_vgd_hl_dflag", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    const size_t q = rb0 + t;
    int fn = 0, fw = 0;
    if (wint && wl_(m,kd+kk,jd+jj) != 0) {
      for (int u = 0; u < 3; ++u) {
        const size_t at = rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj;
        if (fr_(at + 1) > fr_(at)) {fw = 1;}
      }
    } else {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + (vs*kn + kk)*jn + jj;
      fn = (fr_(at + 1) > fr_(at)) ? 1 : 0;
    }
    fn_(q) = fn;
    fw_(q) = fw;
  });
  const DevExeSpace ex_ = vgd_cur;
  VgdScan2(ex_, vgd_hl_fn, rtot, vgd_hl_fw, rtot, "vet_vgd_hl_scan2", vgd_scan_sums);
  int nrw[2] = {0, 0};
  {
    Kokkos::View<int*, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(nrw, 1);
    Kokkos::deep_copy(ex_, h, Kokkos::subview(vgd_hl_fn, std::make_pair(rtot, rtot + 1)));
    Kokkos::View<int*, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h2(nrw + 1,
                                                                                 1);
    Kokkos::deep_copy(ex_, h2,
                      Kokkos::subview(vgd_hl_fw, std::make_pair(rtot, rtot + 1)));
    ex_.fence();
  }
  const int ns = pbh[np] - pbh[0];
  const int nr = nrw[0], nw = nrw[1];
  // every array at least its new size (capacity kept); budget on the sum
  const size_t nsz[5] = {static_cast<size_t>(std::max(ns, 1)),
                         static_cast<size_t>(std::max(nr, 1)),
                         static_cast<size_t>(std::max(nr, 1)),
                         static_cast<size_t>(std::max(nw, 1)),
                         static_cast<size_t>(std::max(3*nw, 1))};
  DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
  size_t have = 0, after = 0;
  for (int q = 0; q < 5; ++q) {
    have += arr[q]->extent(0);
    after += std::max(arr[q]->extent(0), nsz[q]);
  }
  if (after > have) {
    // budget: vet_gd_halo_list_mb, capped on a GPU at half the device memory free at
    // the first build
    if (vgd_hl_cap == 0) {
      vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;
      size_t fr = 0, tt = 0;
#if defined(KOKKOS_ENABLE_CUDA)
      if (cudaMemGetInfo(&fr, &tt) == cudaSuccess) {
        vgd_hl_cap = std::min(vgd_hl_cap, fr/2);
      }
#elif defined(KOKKOS_ENABLE_HIP)
      if (hipMemGetInfo(&fr, &tt) == hipSuccess) {
        vgd_hl_cap = std::min(vgd_hl_cap, fr/2);
      }
#endif
      (void) fr;
      (void) tt;
    }
    const size_t cap = vgd_hl_cap;
    if (vgd_hl_bytes + sizeof(int)*(after - have) > cap) {
      // over the budget: no new attempt for this mask set (dense path)
      E.alpha = vgd_alpha;
      E.scut = vgd_scut;
      E.ws = ws;
      E.nofit = true;
      return false;
    }
    for (int q = 0; q < 5; ++q) {
      if (arr[q]->extent(0) < nsz[q]) {Kokkos::realloc(*arr[q], nsz[q]);}
    }
    vgd_hl_bytes += sizeof(int)*(after - have);
  }
  auto sl_ = E.sl;
  auto dl_ = E.dl;
  auto sr_ = E.sr;
  auto dw_ = E.dw;
  auto sw_ = E.sw;
  // (2) the send list (compact pack addresses)
  par_for("vet_vgd_hl_slist", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const size_t q = static_cast<size_t>(so_(ws,8*m + o))*nv + t;
    if (fs_(q + 1) == fs_(q)) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
    const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
    sl_(fs_(q)) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
  });
  // (3) the receive lists
  par_for("vet_vgd_hl_rlist", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    const size_t q = rb0 + t;
    const int ad = ((m*nv + v)*c3 + kd + kk)*c2 + jd + jj;
    if (fw_(q + 1) > fw_(q)) {
      const int c = fw_(q);
      dw_(c) = ad;
      for (int u = 0; u < 3; ++u) {
        const size_t at = rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj;
        sw_(3*c + u) = (fr_(at + 1) > fr_(at)) ? fr_(at) : -1;
      }
    } else if (fn_(q + 1) > fn_(q)) {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      const size_t at = rb0 + (vs*kn + kk)*jn + jj;
      dl_(fn_(q)) = ad;
      sr_(fn_(q)) = fr_(at);
    }
  });
  E.sp.assign(pbh.begin(), pbh.begin() + np + 1);
  E.rp.assign(pbh.begin() + np + 1, pbh.begin() + 2*np + 2);
  E.ns = ns;
  E.nr = nr;
  E.nw = nw;
  E.alpha = vgd_alpha;
  E.scut = vgd_scut;
  E.ws = ws;
  E.gen = ++vgd_hl_gen;
  E.nofit = false;
  vgd_hl_nmade += 1.0;
  return true;
#else
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHaloCompact
//! \brief exchange of one shell: mask from VetGdHcPrep (prepared during the previous
//! shell's wait), compact pack, MPI, expansion; queue the next mask. false: the caller
//! unpacks the dense buffer. With a valid list entry X (vet_gd_halo_list_mb): gather
//! pack of X, MPI, zeros on the other pass's list Y in X's region, scatter of X (same
//! band values as the dense path; true). Defined below VetGdHlBegin / VetGdHlEnd.

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHlBegin
//! \brief first half of the list path of VetGdHaloCompact (gather pack, fence, posted
//! messages, zeros on Y's entries); false when the list path does not apply.
//! VetGdHlEnd completes it; the sweep overlap runs local items in between.

//----------------------------------------------------------------------------------------
// an event behind the per-shell halo pack on the sweep stream, so the host
// waits for the pack alone while the ring copies queued after it run (VetGdHlBegin)
namespace {
#if defined(KOKKOS_ENABLE_HIP)
hipEvent_t g_vgd_pkev = nullptr;
void VgdPackEvtRecord(const DevExeSpace &ex) {
  if (g_vgd_pkev == nullptr) {
    (void)hipEventCreateWithFlags(&g_vgd_pkev, hipEventDisableTiming);
  }
  (void)hipEventRecord(g_vgd_pkev, ex.hip_stream());
}
void VgdPackEvtWait(const DevExeSpace &) {(void)hipEventSynchronize(g_vgd_pkev);}
#elif defined(KOKKOS_ENABLE_CUDA)
cudaEvent_t g_vgd_pkev = nullptr;
void VgdPackEvtRecord(const DevExeSpace &ex) {
  if (g_vgd_pkev == nullptr) {
    (void)cudaEventCreateWithFlags(&g_vgd_pkev, cudaEventDisableTiming);
  }
  (void)cudaEventRecord(g_vgd_pkev, ex.cuda_stream());
}
void VgdPackEvtWait(const DevExeSpace &) {(void)cudaEventSynchronize(g_vgd_pkev);}
#else
void VgdPackEvtRecord(const DevExeSpace &) {}
void VgdPackEvtWait(const DevExeSpace &ex) {ex.fence();}
#endif
}  // namespace

template <class V>
bool RadiationVet::VetGdHlBegin(V &a, const int nv, const int i0, const int ws) {
#if MPI_PARALLEL_ENABLED
  if (vgd_hlp.on) {VgdFatal("VetGdHlBegin: an exchange is still pending");}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int i = i0;
  const bool inw = vgd_hinw;
  auto a_ = a;
  auto csb_ = vgd_csb;
  auto crb_ = vgd_crb;
  const int np = static_cast<int>(vgd_prk.size());
  // vet_gd_halo_list_mb: list entries X (this pass) and Y (the other pass) of shell i,
  // and the last writers of the band per pass
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  bool lst = (nv == vgd_n) && (i >= 0) && (i < n1c) &&
             (static_cast<int64_t>(nmb1 + 1)*nv*c3*c2 <
              static_cast<int64_t>(std::numeric_limits<int>::max()));
  const int px = inw ? 0 : 1;
  const int ex = px*n1c + i, ey = (1 - px)*n1c + i;
  std::array<std::vector<int64_t>, 2> *gwa = nullptr;
  if (lst) {
    if (static_cast<int>(vgd_hle.size()) != 2*n1c) {vgd_hle.resize(2*n1c);}
    auto gwget = [&](const void *p) {
      auto &g = vgd_hl_gw[p];
      if (static_cast<int>(g[0].size()) != n1c) {
        g[0].assign(n1c, 0);
        g[1].assign(n1c, 0);
      }
      return &g;
    };
    if constexpr (std::is_same<V, VgdIView>::value) {
      gwa = gwget(static_cast<const void *>(a.d.data()));
    } else {
      lst = false;
    }
  }
  auto xvalid = [&]() {
    const VgdHlEntry &X = vgd_hle[ex];
    return X.gen > 0 && X.alpha == vgd_alpha && X.scut == vgd_scut && X.ws == ws;
  };
  auto gwok = [&](std::array<std::vector<int64_t>, 2> *g) {
    const int64_t gx = (*g)[px][i], gy = (*g)[1 - px][i];
    return (gx == 0 || gx == vgd_hle[ex].gen) && (gy == 0 || (vgd_hle[ey].gen > 0 &&
                                                               gy == vgd_hle[ey].gen));
  };
  // vet_gd_halo_exact: only an exact record of this key (no last-writer condition, no
  // zeros: every ghost the sweep reads is on the list)
  const bool xe = vgd_hx_on;
  if (xe) {
    if (!(lst && VetGdHxKey(ex, ws) && vgd_hxe[ex].exact)) {return false;}
  } else if (!(lst && xvalid() && gwok(gwa))) {
    return false;
  }
  if constexpr (std::is_same<V, VgdIView>::value) {
    const VgdHlEntry &X = xe ? vgd_hxe[ex] : vgd_hle[ex];
    auto sl_ = X.sl;
    const int c2_ = c2, c3_ = c3, nv_ = nv;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    // vet_gd_halo_ipc: node-local and GPU-driven: wait for the partners' acks of this
    // parity, gather my values into their receive buffers, set their arrival flags,
    // wait for mine; all on the sweep stream (bitwise the MPI path)
    if (xe && vgd_ipc_state == 1 && static_cast<int>(X.dsto.size()) == np) {
      const long long seq = ++vgd_ipc_seq;
      const long long par = seq & 1;
      long long *sig = vgd_ipc_sig;
      int *err = vgd_ipc_err;
      auto psig = vgd_ipc_psig.d_view;
      auto pslot = vgd_ipc_pslot.d_view;
      auto prb = vgd_ipc_prb.d_view;
      auto phalf = vgd_ipc_phalf.d_view;
      const int np_ = np;
      // lean protocol, no ack kernel: partner p sets my arrival flag of exchange s in
      // its k2(s), after its scatter of s-1; so once my k2(s) has seen it, p has read
      // its parity (s-1)&1 buffer and my pack of s+1 may write there
      VgdIpcIdx spi, doi;
      for (int p = 0; p <= np; ++p) {spi.v[p] = X.sp[p];}
      for (int p = 0; p < np; ++p) {doi.v[p] = X.dsto[p];}
      // fused exchange (CUDA lean, no band-ring copy in between): teams pack as
      // vet_vgd_ipc_pack, the last one sets the partners' flags, every team waits for
      // mine and scatters as vet_vgd_hl_scatter (bitwise). At most 216 teams of 128: all
      // resident, so early spinning teams cannot starve the last one's pack
#if defined(KOKKOS_ENABLE_CUDA)
      if (vgd_hlp.now && vgd_rg_pend.empty() && X.ns + X.nr + X.nw > 0) {
        if (!vgd_hlp.cnt.is_allocated()) {
          vgd_hlp.cnt = Kokkos::View<int, DevMemSpace>("vet_vgd_ipc_cnt");
        }
        auto cnt_ = vgd_hlp.cnt;
        long long *sig = vgd_ipc_sig;
        int *err = vgd_ipc_err;
        const Real *rb_ = vgd_ipc_rbuf + par*vgd_ipc_half;
        const int ns_ = X.ns, nr = X.nr, nw = X.nw;
        auto dl_ = X.dl;
        auto sr_ = X.sr;
        auto dw_ = X.dw;
        auto sw_ = X.sw;
        auto w3_ = vgd_w3;
        constexpr int kT = 128;
        const int nt = std::max(ns_, nr + nw);
        const int nl = std::min((nt + kT - 1)/kT, 216);
        // light-weight launch: body as kernel arguments, no host wait
        Kokkos::parallel_for("vet_vgd_ipc_x", Kokkos::Experimental::require(
            Kokkos::TeamPolicy<>(vgd_cur, nl, kT),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        KOKKOS_LAMBDA(const TeamMember_t &tm) {
          const int c0 = tm.league_rank()*kT + tm.team_rank(), dc = tm.league_size()*kT;
          for (int c = c0; c < ns_; c += dc) {
            int p = 0;
            while (p + 1 < np_ && spi.v[p+1] <= c) {++p;}
            int ad = sl_(c);
            const int jb = ad % c2_;
            ad /= c2_;
            const int kb = ad % c3_;
            ad /= c3_;
            const int v = ad % nv_;
            const int m = ad/nv_;
            Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + doi.v[p] +
                        (c - spi.v[p]);
            *dst = a_(m,v,kb,jb,i);
          }
#if defined(__CUDA_ARCH__)
          __threadfence_system();
#endif
          tm.team_barrier();
          if (tm.team_rank() == 0) {
            if (Kokkos::atomic_fetch_add(&cnt_(), 1) == tm.league_size() - 1) {
              cnt_() = 0;
              VgdIpcSetN(psig.data(), pslot.data(), np_, seq);
            }
            VgdIpcWaitN(sig, np_, seq, err);
          }
          tm.team_barrier();
          for (int c = c0; c < nr + nw; c += dc) {
            int ad = (c < nr) ? dl_(c) : dw_(c - nr);
            const int jb = ad % c2_;
            ad /= c2_;
            const int kb = ad % c3_;
            ad /= c3_;
            const int v = ad % nv_;
            const int m = ad/nv_;
            if (c < nr) {
              const int s = sr_(c);
              a_(m,v,kb,jb,i) = (s >= 0) ? rb_[s] : 0.0;
              continue;
            }
            const int cw = c - nr;
            Real val = 0.0;
            for (int u = 0; u < 3; ++u) {
              const int s = sw_(3*cw + u);
              val += w3_(m,kb,jb,v,u)*((s >= 0) ? rb_[s] : 0.0);
            }
            a_(m,v,kb,jb,i) = val;
          }
        });
        vgd_hlp.done = true;
      } else
#endif
      if (X.ns > 0) {
        Kokkos::parallel_for("vet_vgd_ipc_pack", VgdLight(vgd_cur, 0, X.ns),
        KOKKOS_LAMBDA(const int c) {
          int p = 0;
          while (p + 1 < np_ && spi.v[p+1] <= c) {++p;}
          int ad = sl_(c);
          const int jb = ad % c2_;
          ad /= c2_;
          const int kb = ad % c3_;
          ad /= c3_;
          const int v = ad % nv_;
          const int m = ad/nv_;
          Real *dst = reinterpret_cast<Real *>(prb(p)) + par*phalf(p) + doi.v[p] +
                      (c - spi.v[p]);
          *dst = a_(m,v,kb,jb,i);
        });
      }
      // band-ring copies for this hook: after the pack that reads the slabs, before
      // the scatter that writes the lent slab
      VetGdRingFlush();
      if (!vgd_hlp.done) {
      Kokkos::parallel_for("vet_vgd_ipc_k2", Kokkos::RangePolicy<>(vgd_cur, 0, 1),
      KOKKOS_LAMBDA(const int) {
        VgdIpcSetN(psig.data(), pslot.data(), np_, seq);
        VgdIpcWaitN(sig, np_, seq, err);
      });
      }
      vgd_hlp.on = true;
      vgd_hlp.ipc = true;
      vgd_hlp.seq = seq;
      vgd_hlp.xe = xe;
      vgd_hlp.ex = ex;
      vgd_hlp.px = px;
      vgd_hlp.i = i;
      vgd_hlp.gwa = gwa;
      vgd_hlp.a = a;
      vgd_ipc_n += 1.0;
      return true;
    }
#endif
    if (X.ns > 0) {
      Kokkos::parallel_for("vet_vgd_hl_pack", VgdLight(vgd_cur, 0, X.ns),
      KOKKOS_LAMBDA(const int c) {
        int ad = sl_(c);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        csb_(c) = a_(m,v,kb,jb,i);
      });
    }
    const std::vector<int> &sp = X.sp, &rp_ = X.rp;
    // ring copies (other shells' slabs; stream order pack, ring, zero, scatter) are
    // queued on an event behind the pack, so they run while the host posts the
    // messages (the send buffer is complete at the event: bitwise)
    VgdPackEvtRecord(vgd_cur);
    VetGdRingFlush();
    VgdPackEvtWait(vgd_cur);
    vgd_hlp.tq.reset();
    auto &tq = vgd_hlp.tq;
    auto &req = vgd_hlp.req;
    req.clear();
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_prk[p];
      req.emplace_back();
      MPI_Irecv(crb_.data() + rp_[p], rp_[p+1] - rp_[p], MPI_ATHENA_REAL, rk, 7002,
                vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(csb_.data() + sp[p], sp[p+1] - sp[p], MPI_ATHENA_REAL, rk, 7002,
                vgd_comm, &req.back());
    }
    vgd_tpost += tq.seconds();
    // vet_gd_band_ring: the sweep's slab copies (other shells only) while the messages
    // travel
    VetGdRingFlush();
    // zeros on Y's receive entries inside X's region (the dense unpack writes 0 there)
    const VgdHlEntry &Y = vgd_hle[ey];
    const bool zy = !xe && ((*gwa)[1 - px][i] != 0);
    if (zy && (Y.nr > 0 || Y.nw > 0)) {
      const int ny = Y.nr, nyw = Y.nw;
      auto ydl_ = Y.dl;
      auto ydw_ = Y.dw;
      Kokkos::parallel_for("vet_vgd_hl_zero", VgdLight(vgd_cur, 0, ny + nyw),
      KOKKOS_LAMBDA(const int c) {
        int ad = (c < ny) ? ydl_(c) : ydw_(c - ny);
        const int jb = ad % c2_;
        ad /= c2_;
        const int kb = ad % c3_;
        ad /= c3_;
        const int v = ad % nv_;
        const int m = ad/nv_;
        if (jb < w - ws || jb >= w + nx2 + ws || kb < w - ws || kb >= w + nx3 + ws) {
          return;
        }
        a_(m,v,kb,jb,i) = 0.0;
      });
    }

    vgd_hlp.on = true;
    vgd_hlp.ipc = false;
    vgd_hlp.xe = xe;
    vgd_hlp.ex = ex;
    vgd_hlp.px = px;
    vgd_hlp.i = i;
    vgd_hlp.gwa = gwa;
    vgd_hlp.a = a;
    return true;
  }
#endif
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHlEnd
//! \brief second half of a list-path exchange started by VetGdHlBegin: MPI wait, size
//! check, scatter, last-writer ids

void RadiationVet::VetGdHlEnd() {
#if MPI_PARALLEL_ENABLED
  if (!vgd_hlp.on) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2_ = nx2 + 2*w, c3_ = nx3 + 2*w, nv_ = vgd_n;
  const int i = vgd_hlp.i, px = vgd_hlp.px;
  auto a_ = vgd_hlp.a;
  DvceArray1D<Real> crb_ = vgd_crb;
  auto &req = vgd_hlp.req;
  auto &tq = vgd_hlp.tq;
  auto *gwa = vgd_hlp.gwa;
  const int np = static_cast<int>(vgd_prk.size());
  const VgdHlEntry &X = vgd_hlp.xe ? vgd_hxe[vgd_hlp.ex] : vgd_hle[vgd_hlp.ex];
  const std::vector<int> &sp = X.sp, &rp_ = X.rp;
  if (vgd_hlp.ipc) {
    // vet_gd_halo_ipc: the values are in my IPC buffer of this parity (k2 waited for
    // them); the scatter reads it in place of crb
    const long long par = vgd_hlp.seq & 1;
    crb_ = DvceArray1D<Real>(vgd_ipc_rbuf + par*vgd_ipc_half, vgd_ipc_cs);
    vgd_nexch += 1.0;
    for (int p = 0; p < np; ++p) {
      vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
    }
  } else {
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != rp_[p+1] - rp_[p]) {
      VgdFatal("vet_gd_halo_list: sender and receiver lists disagree");
    }
    vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
  }
  }   // !ipc
  // the scatter (dense unpack of the expanded buffer, on the written entries)
  const int nr = X.nr, nw = X.nw;
  auto dl_ = X.dl;
  auto sr_ = X.sr;
  auto dw_ = X.dw;
  auto sw_ = X.sw;
  auto w3_ = vgd_w3;
  if (nr + nw > 0 && !vgd_hlp.done) {
    Kokkos::parallel_for("vet_vgd_hl_scatter",
                         VgdLight(vgd_cur, 0, nr + nw),
    KOKKOS_LAMBDA(const int c) {
      int ad = (c < nr) ? dl_(c) : dw_(c - nr);
      const int jb = ad % c2_;
      ad /= c2_;
      const int kb = ad % c3_;
      ad /= c3_;
      const int v = ad % nv_;
      const int m = ad/nv_;
      if (c < nr) {
        // exact lists: sr = -1 is a read ghost the masks do not send: 0, as the dense
        // path
        const int s = sr_(c);
        a_(m,v,kb,jb,i) = (s >= 0) ? crb_(s) : 0.0;
        return;
      }
      const int cw = c - nr;
      Real val = 0.0;
      for (int u = 0; u < 3; ++u) {
        const int s = sw_(3*cw + u);
        val += w3_(m,kb,jb,v,u)*((s >= 0) ? crb_(s) : 0.0);
      }
      a_(m,v,kb,jb,i) = val;
    });
  }
  vgd_rg_lastx = vgd_hlp.xe;          // the band of this exchange came by list
  vgd_rg_lastx_e = vgd_hlp.ex;
  if (!vgd_hlp.xe) {
    (*gwa)[px][i] = X.gen;
  }
  vgd_hl_nuse += 1.0;
  vgd_hlp.on = false;
  vgd_hlp.ipc = false;
  vgd_hlp.done = false;
  vgd_hlp.a = VgdRag();
#endif
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdHxKey
//! \brief vet_gd_halo_exact: record e (exact lists or fallback) was made for the
//! current direction set, cut and depth (same answer on every rank)

bool RadiationVet::VetGdHxKey(const int e, const int ws) const {
  if (e < 0 || e >= static_cast<int>(vgd_hxe.size())) {return false;}
  const VgdHlEntry &X = vgd_hxe[e];
  return X.gen != -1 && X.alpha == vgd_alpha && VetGdHxCut(e, X) && X.ws == ws;
}

// record X of (pass, shell l) = entry e holds for the current ray cut: made for it, or
// the move crosses neither l nor l + 1 (its reads are shells l - 1 .. l + 1), so a cut
// move remakes only the lists of the shells it crosses
bool RadiationVet::VetGdHxCut(const int e, const VgdHlEntry &X) const {
  if (X.scut == vgd_scut) {return true;}
  if (X.scut < 0 || vgd_scut < 0) {return false;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int l = e % (indcs.nx1 + 2*indcs.ng) - indcs.is;
  auto sd = [](const int a, const int c) {return (a > c) - (a < c);};
  return sd(l, X.scut) == sd(l, vgd_scut) && sd(l + 1, X.scut) == sd(l + 1, vgd_scut);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHxBuild2
//! \brief vet_gd_halo_exact = 2: the exact lists of record e from the read marks
//! vgd_hx_mk alone (no masks, no dense exchange); collective. (1) sources of marked
//! ghosts, (2) their compact positions, (3) per partner the count and ascending wanted
//! offsets, (4) a global vote (no mark outside depth ws, byte budget), (5) send
//! addresses and read ghosts. Every source of a read ghost is sent.

void RadiationVet::VetGdHxBuild2(const int e, const int i, const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nmb = nmb1 + 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  VgdHlEntry &E = vgd_hxe[e];
  auto hl_ = vgd_hloc;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  const bool wint = vgd_wint;
  auto mk_ = vgd_hx_mk;
  const int stamp = vgd_hx_stamp;
  const DevExeSpace ex_ = vgd_cur;
  auto grow = [](DvceArray1D<int> &a, const size_t n) {
    if (a.extent(0) < n) {Kokkos::realloc(Kokkos::WithoutInitializing, a, n);}
  };
  grow(vgd_hx_need, rtot + 1);
  grow(vgd_hx_t0, rtot + 1);
  grow(vgd_hx_t1, rtot + 1);
  auto nd_ = vgd_hx_need;
  auto f0_ = vgd_hx_t0;
  auto f1_ = vgd_hx_t1;
  Kokkos::deep_copy(ex_, Kokkos::subview(nd_, std::make_pair(size_t(0), rtot + 1)), 0);
  // (1) sources of the marked ghosts
  par_for("vet_vgd_hx_need", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    if (mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) != stamp) {return;}
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    if (wint && wl_(m,kd+kk,jd+jj) != 0) {
      for (int u = 0; u < 3; ++u) {
        nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj) = 1;
      }
    } else {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      nd_(rb0 + (vs*kn + kk)*jn + jj) = 1;
    }
  });
  int nout = 0;           // marked ghosts outside the region
  {
    Kokkos::View<int, HostMemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(&nout);
    Kokkos::deep_copy(ex_, h, vgd_hx_oc);
  }
  // read ghosts (non-wall f0, wall f1) over the receive layout
  par_for("vet_vgd_hx_dflag", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t q = static_cast<size_t>(ro_(ws,8*m + o))*nv + t;
    const bool rd = mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) == stamp;
    const bool wa = wint && wl_(m,kd+kk,jd+jj) != 0;
    f0_(q) = (rd && !wa) ? 1 : 0;
    f1_(q) = (rd && wa) ? 1 : 0;
  });
  // (2) compact receive positions of the needed sources and of the read ghosts
  VgdScan2(ex_, nd_, rtot, nd_, 0, "vet_vgd_hx_scan2", vgd_scan_sums);
  VgdScan2(ex_, f0_, rtot, f1_, rtot, "vet_vgd_hx_scan2", vgd_scan_sums);
  // per partner: needed counts (scanned need at the piece bounds), host
  std::vector<int> rpx(np + 1), spx(np + 1);
  int nrx = 0, nwx = 0;
  {
    if (static_cast<int>(vgd_hx_hb.extent(0)) < 4*np + 6) {
      vgd_hx_hb = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_hb",
                                                                    4*np + 6);
    }
    auto hb_ = vgd_hx_hb;
    const int o2 = 2*np + 4;
    for (int p = 0; p < np; ++p) {hb_(o2 + p) = vgd_prdsp[ws][p]*nv;}
    hb_(o2 + np) = static_cast<int>(rtot);
    const int np_ = np;
    const size_t rt = rtot;
    Kokkos::parallel_for("vet_vgd_hx_bnd", Kokkos::RangePolicy<>(ex_, 0, np + 3),
    KOKKOS_LAMBDA(const int q) {
      if (q <= np_) {
        hb_(q) = nd_(hb_(o2 + q));
      } else if (q == np_ + 1) {
        hb_(q) = f0_(rt);
      } else {
        hb_(q) = f1_(rt);
      }
    });
    ex_.fence();
    for (int p = 0; p <= np; ++p) {rpx[p] = hb_(p);}
    nrx = hb_(np + 1);
    nwx = hb_(np + 2);
  }
  const int nneed = rpx[np];
  // wanted offsets, ascending, relative to the partner's piece start
  grow(vgd_hx_nc, nneed + 1);
  auto rq_ = vgd_hx_nc;
  {
    std::vector<int> pst(np + 1);
    for (int p = 0; p < np; ++p) {pst[p] = vgd_prdsp[ws][p]*nv;}
    pst[np] = static_cast<int>(rtot);
    if (static_cast<int>(vgd_hx_pin.extent(0)) < np + 1) {
      vgd_hx_pin = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_pin",
                                                                     np + 1);
    }
    auto pst_d = vgd_hx_pin;
    for (int p = 0; p <= np; ++p) {pst_d(p) = pst[p];}
    const int np_ = np;
    Kokkos::parallel_for("vet_vgd_hx_req", Kokkos::RangePolicy<>(ex_, 0, rtot),
    KOKKOS_LAMBDA(const size_t q) {
      if (nd_(q + 1) == nd_(q)) {return;}
      int lo = 0, hi = np_ - 1;
      while (lo < hi) {             // partner piece holding q
        const int md = (lo + hi + 1)/2;
        if (static_cast<size_t>(pst_d(md)) <= q) {lo = md;} else {hi = md - 1;}
      }
      rq_(nd_(q)) = static_cast<int>(q) - pst_d(lo);
    });
    ex_.fence();
  }
  // (3) counts, then the offsets
  {
    std::vector<int> rc(np), sc(np);
    for (int p = 0; p < np; ++p) {rc[p] = rpx[p+1] - rpx[p];}
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      req.emplace_back();
      MPI_Irecv(&sc[p], 1, MPI_INT, vgd_prk[p], 7005, vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(&rc[p], 1, MPI_INT, vgd_prk[p], 7005, vgd_comm, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    spx[0] = 0;
    for (int p = 0; p < np; ++p) {spx[p+1] = spx[p] + sc[p];}
  }
  const int nsx = spx[np];
  grow(vgd_hx_ns, nsx + 1);
  auto so2_ = vgd_hx_ns;
  {
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      req.emplace_back();
      MPI_Irecv(so2_.data() + spx[p], spx[p+1] - spx[p], MPI_INT, vgd_prk[p], 7006,
                vgd_comm, &req.back());
      req.emplace_back();
      MPI_Isend(rq_.data() + rpx[p], rpx[p+1] - rpx[p], MPI_INT, vgd_prk[p], 7006,
                vgd_comm, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  }
  // (4) the vote
  const size_t nsz[5] = {static_cast<size_t>(std::max(nsx, 1)),
                         static_cast<size_t>(std::max(nrx, 1)),
                         static_cast<size_t>(std::max(nrx, 1)),
                         static_cast<size_t>(std::max(nwx, 1)),
                         static_cast<size_t>(std::max(3*nwx, 1))};
  DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
  size_t have = 0, after = 0;
  for (int q = 0; q < 5; ++q) {
    have += arr[q]->extent(0);
    after += std::max(arr[q]->extent(0), nsz[q]);
  }
  int ok = (nout == 0) ? 1 : 0;
  if (after > have) {
    if (vgd_hl_cap == 0) {vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;}
    if (vgd_hl_bytes + sizeof(int)*(after - have) > vgd_hl_cap) {ok = 0;}
  }
  int gok = 0;
  MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_MIN, vgd_comm);
  vgd_hx_nout += nout;
  E.alpha = vgd_alpha;
  E.scut = vgd_scut;
  E.ws = ws;
  E.exact = false;
  E.gen = -2;
  if (gok == 0) {vgd_hx_nfall += 1.0; return;}
  if (after > have) {
    for (int q = 0; q < 5; ++q) {
      if (arr[q]->extent(0) < nsz[q]) {Kokkos::realloc(*arr[q], nsz[q]);}
    }
    vgd_hl_bytes += sizeof(int)*(after - have);
  }
  auto sl_ = E.sl;
  auto dl_ = E.dl;
  auto sr_ = E.sr;
  auto dw_ = E.dw;
  auto sw_ = E.sw;
  // (5) send addresses: offset within my piece for partner p -> slot -> (v, k, j)
  {
    // my send pieces of depth ws sorted by start: (start, 8 m + o)
    auto so_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_soff);
    std::vector<std::array<int, 2>> pc;
    for (int t = 0; t < 8*nmb; ++t) {
      if (so_h(ws,t) >= 0) {pc.push_back({so_h(ws,t)*nv, t});}
    }
    std::sort(pc.begin(), pc.end());
    const int npc = static_cast<int>(pc.size());
    if (static_cast<int>(vgd_hx_pin2.extent(0)) < 2*npc + 2*np + 2) {
      vgd_hx_pin2 = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hx_pin2",
                                                                      2*npc + 2*np + 2);
    }
    auto ps_d = vgd_hx_pin2;
    auto ps_h = ps_d;
    for (int q = 0; q < npc; ++q) {ps_h(q) = pc[q][0]; ps_h(npc + q) = pc[q][1];}
    for (int p = 0; p <= np; ++p) {
      ps_h(2*npc + p) = spx[p];
      ps_h(2*npc + np + 1 + p) = (p < np) ? vgd_pdsp[ws][p]*nv : 0;
    }
    const int np_ = np;
    Kokkos::parallel_for("vet_vgd_hx_slist", Kokkos::RangePolicy<>(ex_, 0, nsx),
    KOKKOS_LAMBDA(const int c) {
      int lo = 0, hi = np_ - 1;     // partner of send entry c
      while (lo < hi) {
        const int md = (lo + hi + 1)/2;
        if (ps_d(2*npc + md) <= c) {lo = md;} else {hi = md - 1;}
      }
      const int q = ps_d(2*npc + np_ + 1 + lo) + so2_(c);
      int a = 0, b = npc - 1;       // piece holding q
      while (a < b) {
        const int md = (a + b + 1)/2;
        if (ps_d(md) <= q) {a = md;} else {b = md - 1;}
      }
      const int sid = ps_d(npc + a);
      const int m = sid/8, o = sid % 8;
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int jn = (dj == 0) ? nx2 : ws;
      const int kn = (dk == 0) ? nx3 : ws;
      const int t = q - ps_d(a);
      const int v = t/(kn*jn);
      const int r1 = t - v*kn*jn;
      const int kk = r1/jn;
      const int jj = r1 - kk*jn;
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
      sl_(c) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
    });
  }
  // read ghosts with their compact source positions (every source is sent)
  par_for("vet_vgd_hx_rlist", ex_, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
    const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
    const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
    const size_t q = rb0 + t;
    const int ad = ((m*nv + v)*c3 + kd + kk)*c2 + jd + jj;
    if (f1_(q + 1) > f1_(q)) {
      const int c = f1_(q);
      dw_(c) = ad;
      for (int u = 0; u < 3; ++u) {
        sw_(3*c + u) = nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj);
      }
    } else if (f0_(q + 1) > f0_(q)) {
      const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      dl_(f0_(q)) = ad;
      sr_(f0_(q)) = nd_(rb0 + (vs*kn + kk)*jn + jj);
    }
  });
  ex_.fence();
  E.sp = spx;
  E.rp = rpx;
  E.dsto.clear();
  E.ns = nsx;
  E.nr = nrx;
  E.nw = nwx;
  E.exact = true;
  E.gen = ++vgd_hl_gen;
  vgd_hx_nmade += 1.0;
  vgd_hx_vexact += rpx[np];
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdHxBatch
//! \brief the exact lists (vet_gd_halo_exact = 2) of every (pass, shell) job, built up
//! front in VetGdSweep in groups of G shells: per group (1) device-only marks, flags,
//! scans, counts and offsets, one fence; (2) one count message per partner, the
//! offsets, one vote (MPI_Allreduce of G flags); (3) send / receive lists, one fence.
//! Same expressions as VetGdHxBuild2 (bitwise). G <= 64 from a 256 MB pool (same on
//! every rank); the budget is checked on the group's sum (a shell over it: dense path).

void RadiationVet::VetGdHxBatch(const std::vector<std::array<int, 4>> &jobs) {
#if MPI_PARALLEL_ENABLED
  if (jobs.empty()) {return;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nmb = nmb1 + 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int nv = vgd_n, n = vgd_n;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  const int np = static_cast<int>(vgd_prk.size());
  const int n1 = indcs.nx1, is = indcs.is;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int lcut = vgd_scut;
  const bool imir_ = (ibc_x1min == VET_IBC_REFLECT) && lcut == 0;
  const int wb = w;
  const int js = wb, je = wb + nx2 - 1, ks = wb, ke = wb + nx3 - 1;
  const int jlo = 0, jhi = je + wb - 1, klo = 0, khi = ke + wb - 1;
  const bool sth_ = vgd_sth, pol_ = vgd_pol;
  auto thn_ = vgd_thn;
  const int c2w_ = indcs.nx2 + 2*wb;
  const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
  const int nx2b = nx2, nx3b = nx3;
  auto hl_ = vgd_hloc;
  auto ro_ = vgd_roff;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  auto m3_ = vgd_m3;
  auto dir_ = vgd_dir;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const VgdLG lg_ = vgd_lg;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const bool wint = vgd_wint;
  auto w3_ = vgd_w3;           // mirror theta walls: read marks
  const DevExeSpace ex_ = vgd_cur;
  if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
  Kokkos::Timer tb;
  auto rtot_of = [&](const int ws) -> size_t {
    return static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
           (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nv;
  };
  auto stot_of = [&](const int ws) -> size_t {
    return static_cast<size_t>(vgd_pdsp[ws].empty() ? 0 :
           (vgd_pdsp[ws].back() + vgd_pscnt[ws].back()))*nv;
  };
  // ---- setup (once): slot capacities, G, per-depth partner / piece tables ----
  if (vgd_hxb_g == 0) {
    size_t rcap = 1, scap = 1;
    for (int ws = 1; ws <= w; ++ws) {
      rcap = std::max(rcap, rtot_of(ws) + 1);
      scap = std::max(scap, stot_of(ws) + 1);
    }
    // pool 256 MB (exact lists built in groups, only at a new direction set)
    const double mb = 256.0;
    const int gmax = 64;
    const double per = (4.0*rcap + scap)*sizeof(int);
    int g = static_cast<int>(std::floor(mb*1048576.0/per));
    g = std::max(1, std::min(g, gmax));
    MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    vgd_hxb_g = g;
    vgd_hxb_rcap = rcap;
    vgd_hxb_scap = scap;
    const size_t G = static_cast<size_t>(g);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_nd, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_f0, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_f1, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_rq, G*rcap);
    Kokkos::realloc(Kokkos::WithoutInitializing, vgd_hxb_so, G*scap);
    vgd_hxb_oc = DvceArray1D<int>("vet_vgd_hxb_oc", G);
    vgd_hxb_hb = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hxb_hb",
                                                                   G*(np + 4));
    vgd_hxb_sp = Kokkos::View<int*, Kokkos::SharedHostPinnedSpace>("vet_vgd_hxb_sp",
                                                                   G*(np + 1));
    // per depth ws: the partners' receive piece starts (np + 1, last = rtot), my send
    // pieces sorted by start (start, 8 m + o), the partners' send starts
    auto so_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_soff);
    std::vector<std::array<int, 2>> pc0;
    for (int t = 0; t < 8*nmb; ++t) {
      if (so_h(1,t) >= 0) {pc0.push_back({0, t});}
    }
    const int npc = static_cast<int>(pc0.size());
    vgd_hxb_npc = npc;
    vgd_hxb_pst = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pst",
                                                                    w + 1, np + 1);
    vgd_hxb_pcs = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pcs",
                                                                    w + 1,
                                                                    std::max(npc, 1));
    vgd_hxb_pci = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pci",
                                                                    w + 1,
                                                                    std::max(npc, 1));
    vgd_hxb_pds = Kokkos::View<int**, LayoutWrapper, DevMemSpace>("vet_vgd_hxb_pds",
                                                                    w + 1, np + 1);
    auto pst_h = Kokkos::create_mirror_view(vgd_hxb_pst);
    auto pcs_h = Kokkos::create_mirror_view(vgd_hxb_pcs);
    auto pci_h = Kokkos::create_mirror_view(vgd_hxb_pci);
    auto pds_h = Kokkos::create_mirror_view(vgd_hxb_pds);
    for (int ws = 1; ws <= w; ++ws) {
      for (int p = 0; p < np; ++p) {
        pst_h(ws,p) = vgd_prdsp[ws][p]*nv;
        pds_h(ws,p) = vgd_pdsp[ws][p]*nv;
      }
      pst_h(ws,np) = static_cast<int>(rtot_of(ws));
      pds_h(ws,np) = 0;
      std::vector<std::array<int, 2>> pc;
      for (int t = 0; t < 8*nmb; ++t) {
        if (so_h(ws,t) >= 0) {pc.push_back({so_h(ws,t)*nv, t});}
      }
      std::sort(pc.begin(), pc.end());
      if (static_cast<int>(pc.size()) != npc) {
        VgdFatal("vet_gd_halo_exact batch: the send pieces differ between depths");
      }
      for (int q = 0; q < npc; ++q) {pcs_h(ws,q) = pc[q][0]; pci_h(ws,q) = pc[q][1];}
    }
    Kokkos::deep_copy(vgd_hxb_pst, pst_h);
    Kokkos::deep_copy(vgd_hxb_pcs, pcs_h);
    Kokkos::deep_copy(vgd_hxb_pci, pci_h);
    Kokkos::deep_copy(vgd_hxb_pds, pds_h);
    const size_t nall = static_cast<size_t>(nmb)*n*c3*c2;
    if (vgd_hx_mk.extent(0) < nall) {
      Kokkos::realloc(vgd_hx_mk, nall);
      vgd_hx_stamp = 0;
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_halo_exact lists made in groups of " << g
                << " shells (pool " << (G*per/1048576.0) << " MB per rank)" << std::endl;
    }
  }
  const int G = vgd_hxb_g;
  const size_t rcap = vgd_hxb_rcap, scap = vgd_hxb_scap;
  const int npc = vgd_hxb_npc;
  auto mk_ = vgd_hx_mk;
  auto ocv_ = vgd_hxb_oc;
  auto pst_ = vgd_hxb_pst;
  auto pcs_ = vgd_hxb_pcs;
  auto pci_ = vgd_hxb_pci;
  auto pds_ = vgd_hxb_pds;
  auto hb_ = vgd_hxb_hb;
  auto spp_ = vgd_hxb_sp;
  const int nb = np + 4;
  const int njob = static_cast<int>(jobs.size());
  for (int g0 = 0; g0 < njob; g0 += G) {
    const int gn = std::min(G, njob - g0);
    Kokkos::deep_copy(ex_, ocv_, 0);
    // ---- (1) device: marks, sources, flags, scans, bounds, offsets ----
    for (int s = 0; s < gn; ++s) {
      const int pass = jobs[g0+s][0], i = jobs[g0+s][1], l = jobs[g0+s][2];
      const int ws = jobs[g0+s][3];
      const bool inw = (pass == 0);
      const size_t rtot = rtot_of(ws);
      const int mx = std::max(nx2, nx3)*ws*nv;
      DvceArray1D<int> nd_(vgd_hxb_nd.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f0_(vgd_hxb_f0.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f1_(vgd_hxb_f1.data() + s*rcap, rtot + 1);
      DvceArray1D<int> rq_(vgd_hxb_rq.data() + s*rcap, rtot + 1);
      vgd_hx_stamp += 1;
      const int stamp = vgd_hx_stamp;
      const int wsx = ws;
      // marks (VetGdSweep hx_mark): ghosts of shell i's band read by the stencils of
      // the shells reading this exchange
      auto mark = [&](const int pr, const int ir, const int lr) {
        const bool inwr = (pr == 0);
        auto L_ = vgd_shl[pr];
        const int so = s;
        Kokkos::parallel_for("vet_vgd_hxb_mark", Kokkos::RangePolicy<>(ex_, 0,
                                                                       vgd_shn[pr]),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          VgdCand(m, k, j, d, ir, lr, inwr, n1, lcut, js, ks, jlo, jhi, klo, khi, dir_,
                  mbsize, lg_, cx1v, cx1f, sth_, pol_, thn_, c2w_, imir_,
          [&](const int iu, const int ka, const int kb, const int ja, const int jb,
              const int dr) {
            if (iu != i) {return;}
            for (int kk = ka; kk <= kb; ++kk) {
              for (int jj = ja; jj <= jb; ++jj) {
                const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
                const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
                if (ek == 0 && ej == 0) {continue;}
                if (ej != 0 && wl_(m,kk,jj) >= 2) {
                  // the main sweep reads the mirror row (flag 2) or the edge row
                  // (extrap, 4) through the phi-edge slot; vacuum (3) reads nothing.
                  // Flags 3, 4: the twin's periodic reads below too
                  const int wf = wl_(m,kk,jj);
                  const int o2 = 3*(ek + 1) + 1;
                  if (wf != 3 && ek != 0 && hl_(8*m + ((o2 < 4) ? o2 : (o2 - 1))) < 0) {
                    const bool hi = (ej > 0);
                    const int g = hi ? (jj - wb - nx2b) : (wb - 1 - jj);
                    const int jm = (wf == 4) ? (hi ? (wb + nx2b - 1) : wb)
                                             : VgdMirRow(g, hi, wb, nx2b);
                    for (int t = 0; t < ((wint && wf == 2) ? 3 : 1); ++t) {
                      if (t > 0 && w3_(m,kk,jj,dr,t) == 0.0) {continue;}
                      const int dq = (wf == 4) ? dr : (wint ? m3_(m,kk,jj,dr,t)
                                                           : mp_(m,kk,jj,dr));
                      mk_(((m*n + dq)*c3 + kk)*c2 + jm) = stamp;
                    }
                    if (kk < wb - wsx || kk >= wb + nx3b + wsx) {
                      Kokkos::atomic_add(&ocv_(so), 1);
                    }
                  }
                  if (wf == 2) {continue;}
                }
                const int oo = 3*(ek + 1) + (ej + 1);
                if (hl_(8*m + ((oo < 4) ? oo : (oo - 1))) >= 0) {continue;}
                mk_(((m*n + dr)*c3 + kk)*c2 + jj) = stamp;
                if (jj < wb - wsx || jj >= wb + nx2b + wsx || kk < wb - wsx ||
                    kk >= wb + nx3b + wsx) {
                  Kokkos::atomic_add(&ocv_(so), 1);
                }
              }
            }
          });
        });
      };
      if (inw) {
        if (l - 1 >= lcut) {mark(0, i - 1, l - 1);}
        mark(1, i, l);
      } else if (l + 1 <= n1 - 1) {
        mark(1, i + 1, l + 1);
      }
      (void) is;
      Kokkos::deep_copy(ex_, nd_, 0);
      // sources of the marked ghosts and read-ghost flags (VetGdHxBuild2 (1))
      par_for("vet_vgd_hxb_need", ex_, 0, nmb1, 0, 7, 0, mx - 1,
      KOKKOS_LAMBDA(const int m, const int o, const int t) {
        if (hl_(8*m + o) >= 0) {return;}
        const int oo = (o < 4) ? o : (o + 1);
        const int dk = oo/3 - 1, dj = oo%3 - 1;
        const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
        const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
        const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
        if (t >= nv*kn*jn) {return;}
        const int v = t/(kn*jn);
        const int r1 = t - v*kn*jn;
        const int kk = r1/jn;
        const int jj = r1 - kk*jn;
        const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
        const size_t q = rb0 + t;
        const bool rd = mk_(((m*nv + v)*c3 + kd + kk)*c2 + jd + jj) == stamp;
        const bool wa = wint && wl_(m,kd+kk,jd+jj) != 0;
        f0_(q) = (rd && !wa) ? 1 : 0;
        f1_(q) = (rd && wa) ? 1 : 0;
        if (!rd) {return;}
        if (wa) {
          for (int u = 0; u < 3; ++u) {
            nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj) = 1;
          }
        } else {
          const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
          nd_(rb0 + (vs*kn + kk)*jn + jj) = 1;
        }
      });
      VgdScanOne(ex_, nd_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      VgdScanOne(ex_, f0_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      VgdScanOne(ex_, f1_, rtot, vgd_scan_sums, vgd_hxb_tmp);
      // needed counts at the partner bounds, read-ghost counts, out count
      {
        const int np_ = np, so = s, nb_ = nb;
        const size_t rt = rtot;
        Kokkos::parallel_for("vet_vgd_hxb_bnd", Kokkos::RangePolicy<>(ex_, 0, np + 4),
        KOKKOS_LAMBDA(const int q) {
          int v = 0;
          if (q <= np_) {
            v = nd_(pst_(ws,q));
          } else if (q == np_ + 1) {
            v = f0_(rt);
          } else if (q == np_ + 2) {
            v = f1_(rt);
          } else {
            v = ocv_(so);
          }
          hb_(so*nb_ + q) = v;
        });
      }
      // wanted offsets, ascending, relative to the partner's piece start
      {
        const int np_ = np;
        Kokkos::parallel_for("vet_vgd_hxb_req", Kokkos::RangePolicy<>(ex_, 0, rtot),
        KOKKOS_LAMBDA(const size_t q) {
          if (nd_(q + 1) == nd_(q)) {return;}
          int lo = 0, hi = np_ - 1;
          while (lo < hi) {
            const int md = (lo + hi + 1)/2;
            if (static_cast<size_t>(pst_(ws,md)) <= q) {lo = md;} else {hi = md - 1;}
          }
          rq_(nd_(q)) = static_cast<int>(q) - pst_(ws,lo);
        });
      }
    }
    ex_.fence();
    // ---- (2) host + MPI: counts, offsets, vote ----
    std::vector<std::vector<int>> rpx(gn, std::vector<int>(np + 1));
    std::vector<std::vector<int>> spx(gn, std::vector<int>(np + 1));
    std::vector<int> nrx(gn), nwx(gn), nout(gn);
    // per partner one message of 2 gn ints: the counts I need from it and where its
    // values start in my receive layout (rpx: the IPC sender writes into my buffer)
    std::vector<int> rc(static_cast<size_t>(2*np)*gn), sc(static_cast<size_t>(2*np)*gn);
    for (int s = 0; s < gn; ++s) {
      for (int p = 0; p <= np; ++p) {rpx[s][p] = hb_(s*nb + p);}
      nrx[s] = hb_(s*nb + np + 1);
      nwx[s] = hb_(s*nb + np + 2);
      nout[s] = hb_(s*nb + np + 3);
      for (int p = 0; p < np; ++p) {
        rc[2*p*gn + s] = rpx[s][p+1] - rpx[s][p];
        rc[2*p*gn + gn + s] = rpx[s][p];
      }
    }
    {
      std::vector<MPI_Request> req;
      req.reserve(2*np);
      for (int p = 0; p < np; ++p) {
        req.emplace_back();
        MPI_Irecv(sc.data() + 2*p*gn, 2*gn, MPI_INT, vgd_prk[p], 7008, vgd_comm,
                  &req.back());
        req.emplace_back();
        MPI_Isend(rc.data() + 2*p*gn, 2*gn, MPI_INT, vgd_prk[p], 7008, vgd_comm,
                  &req.back());
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    }
    std::vector<std::vector<int>> dsto(gn, std::vector<int>(np));
    for (int s = 0; s < gn; ++s) {
      spx[s][0] = 0;
      for (int p = 0; p < np; ++p) {
        spx[s][p+1] = spx[s][p] + sc[2*p*gn + s];
        dsto[s][p] = sc[2*p*gn + gn + s];
      }
      if (static_cast<size_t>(spx[s][np]) + 1 > scap) {
        VgdFatal("vet_gd_halo_exact batch: more requested values than the send layout");
      }
    }
    {
      std::vector<MPI_Request> req;
      for (int s = 0; s < gn; ++s) {
        int *so2 = vgd_hxb_so.data() + s*scap;
        int *rq = vgd_hxb_rq.data() + s*rcap;
        for (int p = 0; p < np; ++p) {
          const int nrc = spx[s][p+1] - spx[s][p], nsd = rpx[s][p+1] - rpx[s][p];
          if (nrc > 0) {
            req.emplace_back();
            MPI_Irecv(so2 + spx[s][p], nrc, MPI_INT, vgd_prk[p], 7009, vgd_comm,
                      &req.back());
          }
          if (nsd > 0) {
            req.emplace_back();
            MPI_Isend(rq + rpx[s][p], nsd, MPI_INT, vgd_prk[p], 7009, vgd_comm,
                      &req.back());
          }
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    }
    // vote: no marked ghost outside the region, and the lists fit the budget (counted
    // over the group as if every earlier shell of it is made)
    std::vector<int> ok(gn), gok(gn);
    std::vector<std::array<size_t, 5>> tgt(gn);
    size_t grow = 0;
    if (vgd_hl_cap == 0) {vgd_hl_cap = static_cast<size_t>(vgd_hl_mb) << 20;}
    for (int s = 0; s < gn; ++s) {
      const int e = jobs[g0+s][0]*n1c + jobs[g0+s][1];
      VgdHlEntry &E = vgd_hxe[e];
      const size_t nsz[5] = {static_cast<size_t>(std::max(spx[s][np], 1)),
                             static_cast<size_t>(std::max(nrx[s], 1)),
                             static_cast<size_t>(std::max(nrx[s], 1)),
                             static_cast<size_t>(std::max(nwx[s], 1)),
                             static_cast<size_t>(std::max(3*nwx[s], 1))};
      DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
      size_t have = 0, after = 0;
      for (int q = 0; q < 5; ++q) {
        // growth with 1/4 headroom (sizes change with the angle)
        tgt[s][q] = (arr[q]->extent(0) >= nsz[q]) ? arr[q]->extent(0) :
                    (nsz[q] + nsz[q]/4 + 64);
        have += arr[q]->extent(0);
        after += tgt[s][q];
      }
      ok[s] = (nout[s] == 0) ? 1 : 0;
      if (after > have) {
        if (vgd_hl_bytes + grow + sizeof(int)*(after - have) > vgd_hl_cap) {
          ok[s] = 0;
        } else {
          grow += sizeof(int)*(after - have);
        }
      }
    }
    MPI_Allreduce(ok.data(), gok.data(), gn, MPI_INT, MPI_MIN, vgd_comm);
    // ---- (3) device: the lists ----
    for (int s = 0; s < gn; ++s) {
      const int pass = jobs[g0+s][0], i = jobs[g0+s][1], ws = jobs[g0+s][3];
      const int e = pass*n1c + i;
      VgdHlEntry &E = vgd_hxe[e];
      vgd_hx_nout += nout[s];
      E.alpha = vgd_alpha;
      E.scut = vgd_scut;
      E.ws = ws;
      E.exact = false;
      E.gen = -2;
      if (gok[s] == 0) {vgd_hx_nfall += 1.0; continue;}
      DvceArray1D<int> *arr[5] = {&E.sl, &E.dl, &E.sr, &E.dw, &E.sw};
      for (int q = 0; q < 5; ++q) {
        if (arr[q]->extent(0) < tgt[s][q]) {
          vgd_hl_bytes += sizeof(int)*(tgt[s][q] - arr[q]->extent(0));
          Kokkos::realloc(Kokkos::WithoutInitializing, *arr[q], tgt[s][q]);
          vgd_hxb_nrealloc += 1.0;
        }
      }
      const size_t rtot = rtot_of(ws);
      const int mx = std::max(nx2, nx3)*ws*nv;
      DvceArray1D<int> nd_(vgd_hxb_nd.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f0_(vgd_hxb_f0.data() + s*rcap, rtot + 1);
      DvceArray1D<int> f1_(vgd_hxb_f1.data() + s*rcap, rtot + 1);
      DvceArray1D<int> so2_(vgd_hxb_so.data() + s*scap, scap);
      auto sl_ = E.sl;
      auto dl_ = E.dl;
      auto sr_ = E.sr;
      auto dw_ = E.dw;
      auto sw_ = E.sw;
      const int nsx = spx[s][np];
      for (int p = 0; p <= np; ++p) {spp_(s*(np + 1) + p) = spx[s][p];}
      // send addresses (VetGdHxBuild2 (5))
      {
        const int np_ = np, so = s;
        Kokkos::parallel_for("vet_vgd_hxb_slist", Kokkos::RangePolicy<>(ex_, 0, nsx),
        KOKKOS_LAMBDA(const int c) {
          int lo = 0, hi = np_ - 1;     // partner of send entry c
          while (lo < hi) {
            const int md = (lo + hi + 1)/2;
            if (spp_(so*(np_ + 1) + md) <= c) {lo = md;} else {hi = md - 1;}
          }
          const int q = pds_(ws,lo) + so2_(c);
          int a = 0, b = npc - 1;       // piece holding q
          while (a < b) {
            const int md = (a + b + 1)/2;
            if (pcs_(ws,md) <= q) {a = md;} else {b = md - 1;}
          }
          const int sid = pci_(ws,a);
          const int m = sid/8, o = sid % 8;
          const int oo = (o < 4) ? o : (o + 1);
          const int dk = oo/3 - 1, dj = oo%3 - 1;
          const int jn = (dj == 0) ? nx2 : ws;
          const int kn = (dk == 0) ? nx3 : ws;
          const int t = q - pcs_(ws,a);
          const int v = t/(kn*jn);
          const int r1 = t - v*kn*jn;
          const int kk = r1/jn;
          const int jj = r1 - kk*jn;
          const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
          const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
          sl_(c) = ((m*nv + v)*c3 + ks2 + kk)*c2 + js2 + jj;
        });
      }
      // read ghosts with their compact source positions
      par_for("vet_vgd_hxb_rlist", ex_, 0, nmb1, 0, 7, 0, mx - 1,
      KOKKOS_LAMBDA(const int m, const int o, const int t) {
        if (hl_(8*m + o) >= 0) {return;}
        const int oo = (o < 4) ? o : (o + 1);
        const int dk = oo/3 - 1, dj = oo%3 - 1;
        const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
        const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
        const int jn = (dj == 0) ? nx2 : ws, kn = (dk == 0) ? nx3 : ws;
        if (t >= nv*kn*jn) {return;}
        const int v = t/(kn*jn);
        const int r1 = t - v*kn*jn;
        const int kk = r1/jn;
        const int jj = r1 - kk*jn;
        const size_t rb0 = static_cast<size_t>(ro_(ws,8*m + o))*nv;
        const size_t q = rb0 + t;
        const int ad = ((m*nv + v)*c3 + kd + kk)*c2 + jd + jj;
        if (f1_(q + 1) > f1_(q)) {
          const int c = f1_(q);
          dw_(c) = ad;
          for (int u = 0; u < 3; ++u) {
            sw_(3*c + u) = nd_(rb0 + (m3_(m,kd+kk,jd+jj,v,u)*kn + kk)*jn + jj);
          }
        } else if (f0_(q + 1) > f0_(q)) {
          const int vs = (wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
          dl_(f0_(q)) = ad;
          sr_(f0_(q)) = nd_(rb0 + (vs*kn + kk)*jn + jj);
        }
      });
      E.sp = spx[s];
      E.rp = rpx[s];
      E.dsto = dsto[s];
      E.ns = nsx;
      E.nr = nrx[s];
      E.nw = nwx[s];
      E.exact = true;
      E.gen = ++vgd_hl_gen;
      vgd_hx_nmade += 1.0;
      vgd_hx_vexact += rpx[s][np];
    }
    ex_.fence();
    vgd_hxb_ngroup += 1.0;
  }
  vgd_hxb_t += tb.seconds();
#endif
}

template <class V>
bool RadiationVet::VetGdHaloCompact(V &a, const int nv, const int i0, const int ws) {
#if MPI_PARALLEL_ENABLED
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int mx = std::max(nx2, nx3)*ws*nv;
  const int i = i0;
  const bool inw = vgd_hinw;
  auto hl_ = vgd_hloc;
  auto a_ = a;
  auto so_ = vgd_soff;
  auto csb_ = vgd_csb;
  auto crb_ = vgd_crb;
  auto rb_ = vgd_rbuf;
  const int nvi = nv;
  const size_t rtot = static_cast<size_t>(vgd_prdsp[ws].empty() ? 0 :
                      (vgd_prdsp[ws].back() + vgd_prcnt[ws].back()))*nvi;
  const int np = static_cast<int>(vgd_prk.size());
  // vet_gd_halo_list_mb: list entries X (this pass) and Y (the other pass) of shell i,
  // and the last writers of the band per pass
  const int n1c = indcs.nx1 + 2*indcs.ng;
  const int c2 = nx2 + 2*w, c3 = nx3 + 2*w;
  bool lst = (nv == vgd_n) && (i >= 0) && (i < n1c) &&
             (static_cast<int64_t>(nmb1 + 1)*nv*c3*c2 <
              static_cast<int64_t>(std::numeric_limits<int>::max()));
  const int px = inw ? 0 : 1;
  const int ex = px*n1c + i;
  std::array<std::vector<int64_t>, 2> *gwa = nullptr;
  if (lst) {
    if (static_cast<int>(vgd_hle.size()) != 2*n1c) {vgd_hle.resize(2*n1c);}
    auto gwget = [&](const void *p) {
      auto &g = vgd_hl_gw[p];
      if (static_cast<int>(g[0].size()) != n1c) {
        g[0].assign(n1c, 0);
        g[1].assign(n1c, 0);
      }
      return &g;
    };
    if constexpr (std::is_same<V, VgdIView>::value) {
      gwa = gwget(static_cast<const void *>(a.d.data()));
    } else {
      lst = false;
    }
  }
  auto xvalid = [&]() {
    const VgdHlEntry &X = vgd_hle[ex];
    return X.gen > 0 && X.alpha == vgd_alpha && X.scut == vgd_scut && X.ws == ws;
  };
  // ---- list path (begin + end; Begin decides: the mask lists X, or with
  // vet_gd_halo_exact X's exact record) ----
  vgd_hlp.now = true;
  if (lst && VetGdHlBegin(a, nv, i0, ws)) {
    vgd_hlp.now = false;
    VetGdHlEnd();
    return true;
  }
  vgd_hlp.now = false;
  // ---- the dense path ----
  // the slot holding this shell's mask (prepared by the previous exchange), else now
  auto tagged = [&](const int sl) {
    return vgd_htag[sl][0] == i && vgd_htag[sl][1] == (inw ? 1 : 0) &&
           vgd_htag[sl][2] == ws && vgd_htag[sl][3] == vgd_scut &&
           vgd_htag[sl][4] == vgd_hsweep;
  };
  int slot = -1;
  if (tagged(0)) {slot = 0;}
  if (tagged(1)) {slot = 1;}
  if (slot < 0) {
    slot = 1 - vgd_hlast;
    VetGdHcPrep(slot, i, inw, ws);
  }
  vgd_hlast = slot;
  auto fs_ = vgd_hfs[slot];
  auto fr_ = vgd_hfr[slot];
  // (3) compact pack
  par_for("vet_vgd_hc_pack", vgd_cur, 0, nmb1, 0, 7, 0, mx - 1,
  KOKKOS_LAMBDA(const int m, const int o, const int t) {
    if (hl_(8*m + o) >= 0) {return;}
    const int oo = (o < 4) ? o : (o + 1);
    const int dk = oo/3 - 1, dj = oo%3 - 1;
    const int jn = (dj == 0) ? nx2 : ws;
    const int kn = (dk == 0) ? nx3 : ws;
    if (t >= nv*kn*jn) {return;}
    const size_t q = static_cast<size_t>(so_(ws,8*m + o))*nvi + t;
    if (fs_(q + 1) == fs_(q)) {return;}
    const int v = t/(kn*jn);
    const int r1 = t - v*kn*jn;
    const int kk = r1/jn;
    const int jj = r1 - kk*jn;
    const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
    const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
    csb_(fs_(q)) = a_(m,v,ks2+kk,js2+jj,i);
  });
  vgd_cur.fence();
  std::vector<int> sp(np + 1), rp_(np + 1);
  {
    auto pbh = vgd_hpb[slot];
    for (int p = 0; p <= np; ++p) {
      sp[p] = pbh(p);
      rp_[p] = pbh(np + 1 + p);
    }
  }
  // vet_gd_halo_list_mb: X's lists from this exchange's masks (when X is not valid for
  // them), before the pipelined prep below reuses the masks' slot.
  // vet_gd_halo_exact: the exact lists when the record does not match (collective:
  // no rank has a usable record for this key)
  if (vgd_hx_on && lst) {
    if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
    vgd_hl_ndense += 1.0;
  }
  if (lst && !vgd_hx_on) {
    const VgdHlEntry &X0 = vgd_hle[ex];
    const bool nofit = X0.nofit && X0.alpha == vgd_alpha && X0.scut == vgd_scut &&
                       X0.ws == ws;
    if (!xvalid() && !nofit) {
      std::vector<int> pbv(sp);
      pbv.insert(pbv.end(), rp_.begin(), rp_.end());
      VetGdHlBuild(ex, i, inw, ws, fs_, fr_, pbv);
    }
    // the dense unpack writes X's whole region: the band follows X's lists when X is
    // valid, else it is unknown to the list path
    const int64_t g = xvalid() ? vgd_hle[ex].gen : -1;
    (*gwa)[px][i] = g;
    vgd_hl_ndense += 1.0;
  }
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_prk[p];
    req.emplace_back();
    MPI_Irecv(crb_.data() + rp_[p], rp_[p+1] - rp_[p], MPI_ATHENA_REAL, rk, 7002,
              vgd_comm, &req.back());
    req.emplace_back();
    MPI_Isend(csb_.data() + sp[p], sp[p+1] - sp[p], MPI_ATHENA_REAL, rk, 7002,
              vgd_comm, &req.back());
  }
  vgd_tpost += tq.seconds();
  // the next shell's mask on the device while the messages travel
  if (vgd_hnext[0] >= 0) {
    VetGdHcPrep(1 - slot, vgd_hnext[0], vgd_hnext[1] != 0, vgd_hnext[2]);
  }
  std::vector<MPI_Status> stat(req.size());
  MPI_Waitall(static_cast<int>(req.size()), req.data(), stat.data());
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  for (int p = 0; p < np; ++p) {
    int got = 0;
    MPI_Get_count(&stat[2*p], MPI_ATHENA_REAL, &got);
    if (got != rp_[p+1] - rp_[p]) {
      VgdFatal("vet_gd_halo_compact: sender and receiver masks disagree");
    }
    vgd_nbyte += static_cast<Real>(sp[p+1] - sp[p])*sizeof(Real);
  }
  // (5) expansion to the dense receive layout
  Kokkos::parallel_for("vet_vgd_hc_expand", Kokkos::RangePolicy<>(vgd_cur, 0, rtot),
  KOKKOS_LAMBDA(const size_t q) {
    rb_(q) = (fr_(q + 1) > fr_(q)) ? crb_(fr_(q)) : 0.0;
  });
  return false;
#else
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdIpcInit
//! \brief vet_gd_halo_ipc: node-local CUDA IPC halo of the exact lists: raw receive
//! buffer (2 parities x half, cs the partner stride) and flags, handles swapped with
//! every partner and opened. Off unless every partner is on this node. Collective.

void RadiationVet::VetGdIpcInit(const std::vector<int> &prk, size_t cs, size_t half) {
  vgd_ipc_state = -1;
#if MPI_PARALLEL_ENABLED && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  const int me = global_variable::my_rank;
  const int np = static_cast<int>(prk.size());
  int ok = (np <= VGD_IPC_NS) ? 1 : 0;
  MPI_Comm node;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, me, MPI_INFO_NULL, &node);
  int nn;
  MPI_Comm_size(node, &nn);
  std::vector<int> nodeg(nn);
  MPI_Allgather(&me, 1, MPI_INT, nodeg.data(), 1, MPI_INT, node);
  MPI_Comm_free(&node);
  // vgd_comm is MPI_COMM_WORLD or a duplicate (vet_gd_async): the same ranks
  for (int p = 0; p < np; ++p) {
    if (std::find(nodeg.begin(), nodeg.end(), prk[p]) == nodeg.end()) {ok = 0;}
  }
  vgd_ipc_cs = cs;
  vgd_ipc_half = half;
  if (ok) {
    if (VETG(Malloc)(reinterpret_cast<void **>(&vgd_ipc_rbuf),
                     2*vgd_ipc_half*sizeof(Real))
        != VETG(Success)) {ok = 0;}
#if defined(KOKKOS_ENABLE_HIP)
    if (ok && hipExtMallocWithFlags(reinterpret_cast<void **>(&vgd_ipc_sig),
                                    2*VGD_IPC_NS*sizeof(long long),
                                    hipDeviceMallocUncached) != hipSuccess) {ok = 0;}
#else
    if (ok && cudaMalloc(reinterpret_cast<void **>(&vgd_ipc_sig),
                         2*VGD_IPC_NS*sizeof(long long)) != cudaSuccess) {ok = 0;}
#endif
    if (ok && VETG(Malloc)(reinterpret_cast<void **>(&vgd_ipc_err), sizeof(int))
        != VETG(Success)) {ok = 0;}
    if (ok) {
      (void)VETG(Memset)(vgd_ipc_sig, 0, 2*VGD_IPC_NS*sizeof(long long));
      (void)VETG(Memset)(vgd_ipc_err, 0, sizeof(int));
      (void)VETG(Memset)(vgd_ipc_rbuf, 0, 2*vgd_ipc_half*sizeof(Real));
      (void)VETG(DeviceSynchronize)();
    }
  }
  {
    int g = ok;
    MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    ok = g;
  }
  if (ok) {
    VETG(IpcMemHandle_t) hb, hs;
    (void)VETG(IpcGetMemHandle)(&hb, vgd_ipc_rbuf);
    (void)VETG(IpcGetMemHandle)(&hs, vgd_ipc_sig);
    struct Info {VETG(IpcMemHandle_t) b, s; long long half, cs; int slot;};
    std::vector<Info> mine(np), theirs(np);
    std::vector<MPI_Request> rq(2*np);
    for (int p = 0; p < np; ++p) {
      mine[p] = Info{hb, hs, static_cast<long long>(vgd_ipc_half),
                     static_cast<long long>(vgd_ipc_cs), p};
      MPI_Irecv(&theirs[p], sizeof(Info), MPI_BYTE, prk[p], 7010, vgd_comm, &rq[p]);
    }
    for (int p = 0; p < np; ++p) {
      MPI_Isend(&mine[p], sizeof(Info), MPI_BYTE, prk[p], 7010, vgd_comm,
                &rq[np + p]);
    }
    MPI_Waitall(2*np, rq.data(), MPI_STATUSES_IGNORE);
    vgd_ipc_prb = DualArray1D<unsigned long long>("vet_vgd_ipc_prb", std::max(np, 1));
    vgd_ipc_psig = DualArray1D<unsigned long long>("vet_vgd_ipc_psig", std::max(np, 1));
    vgd_ipc_phalf = DualArray1D<long long>("vet_vgd_ipc_phalf", std::max(np, 1));
    vgd_ipc_pcs = DualArray1D<long long>("vet_vgd_ipc_pcs", std::max(np, 1));
    vgd_ipc_pslot = DualArray1D<int>("vet_vgd_ipc_pslot", std::max(np, 1));
    vgd_ipc_open.assign(2*np, nullptr);
    for (int p = 0; p < np; ++p) {
      void *pb = nullptr, *ps = nullptr;
      if (VETG(IpcOpenMemHandle)(&pb, theirs[p].b, VETG(IpcMemLazyEnablePeerAccess))
          != VETG(Success)) {ok = 0;}
      if (VETG(IpcOpenMemHandle)(&ps, theirs[p].s, VETG(IpcMemLazyEnablePeerAccess))
          != VETG(Success)) {ok = 0;}
      vgd_ipc_open[2*p] = pb;
      vgd_ipc_open[2*p + 1] = ps;
      vgd_ipc_prb.h_view(p) = reinterpret_cast<unsigned long long>(pb);
      vgd_ipc_psig.h_view(p) = reinterpret_cast<unsigned long long>(ps);
      vgd_ipc_phalf.h_view(p) = theirs[p].half;
      vgd_ipc_pcs.h_view(p) = theirs[p].cs;
      vgd_ipc_pslot.h_view(p) = theirs[p].slot;
    }
    (void)VETG(GetLastError)();
    vgd_ipc_prb.modify_host(); vgd_ipc_prb.sync_device();
    vgd_ipc_psig.modify_host(); vgd_ipc_psig.sync_device();
    vgd_ipc_phalf.modify_host(); vgd_ipc_phalf.sync_device();
    vgd_ipc_pcs.modify_host(); vgd_ipc_pcs.sync_device();
    vgd_ipc_pslot.modify_host(); vgd_ipc_pslot.sync_device();
    int g = ok;
    MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    ok = g;
    MPI_Barrier(vgd_comm);   // every peer's flags are zero before any exchange
  }
  vgd_ipc_state = ok ? 1 : -1;
  vgd_ipc_seq = 0;
  if (me == 0) {
    std::cout << "<rad_vet> vet_gd_halo_ipc: "
              << (ok ? "ON (node-local GPU IPC, GPU-driven per-shell halo of the exact "
                       "lists, no host wait)"
                     : "not possible (a partner off-node, or IPC failed): MPI path")
              << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdIpcCheck
//! \brief vet_gd_halo_ipc: after a (fenced) sweep, fatal if a spin timed out

void RadiationVet::VetGdIpcCheck() {
#if MPI_PARALLEL_ENABLED && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  if (vgd_ipc_state != 1) {return;}
  int h = 0;
  (void)VETG(Memcpy)(&h, vgd_ipc_err, sizeof(int), VETG(MemcpyDeviceToHost));
  if (h != 0) {
    std::cout << "### FATAL ERROR in rad_vet vet_gd_halo_ipc: a halo spin timed out "
                 "(rank "
              << global_variable::my_rank << ", exchange " << vgd_ipc_seq << ")"
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdX1Hand
//! \brief several MeshBlocks along x1: the boundary shell of every level to the x1
//! ghost shell of the level above (up) or below, its band exchanged as a sweep shell

void RadiationVet::VetGdX1Hand(const bool up) {
  if (vgd_scut != 0) {VgdFatal("several MeshBlocks along x1 need vray_icut = 0 (the "
                               "sweep from the inner face)");}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w;
  const int nx2 = indcs.nx2, nk = indcs.nx3*nx2;
  const int isrc = up ? indcs.ie : indcs.is, idst = up ? (indcs.is - 1) : (indcs.ie + 1);
  auto vi_ = vgd_i;
  auto b_ = vgd_x1b;
  auto r_ = vgd_x1r;
  par_for("vet_vgd_x1pack", vgd_cur, 0, nmb1, 0, n - 1, 0, nk - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int t) {
    b_(m,d*nk + t) = vi_(m,d,w + t/nx2,w + t % nx2,isrc);
  });
  VetGdX1Pass(up);
  par_for("vet_vgd_x1unpack", vgd_cur, 0, nmb1, 0, n - 1, 0, nk - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int t) {
    vi_(m,d,w + t/nx2,w + t % nx2,idst) = r_(m,d*nk + t);
  });
  VetGdHalo(vgd_i, n, idst, idst, w, true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTopMirror
//! \brief reflecting outer x1 (clean scheme): the inflow of every inward direction at
//! the top face = the top cell's last-sweep intensity in the specular image
//! n - 2 (n . r) r (nearest set direction), scaled so the inflow returns the outflow's
//! radial flux; the top cell's isotropic E when there is no outflow yet

void RadiationVet::VetGdTopMirror() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, ie = indcs.ie;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, ks = indcs.ks, js = indcs.js;
  auto dir_ = vgd_dir;
  auto tg_ = vgd_trig;
  auto mp_ = vgd_rmap;
  if (mp_.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(vgd_rmap, nmb1 + 1, nx3, nx2, n);
    mp_ = vgd_rmap;
  }
  // every sweep: a rotating set changes the directions
  par_for("vet_vgd_rmap", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1, 0, n - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
    Real best = -2.0;
    for (int q = 0; q < n; ++q) {
      Real c = 0.0;
      for (int a = 0; a < 3; ++a) {c += dir_(q,a)*(dir_(d,a) - 2.0*mr*r[a]);}
      if (c > best + 1.0e-12) {best = c; mp_(m,k,j,d) = q;}
    }
  });
  auto vi_ = vgd_i;
  auto iw_ = iw;
  auto xf_ = vgd_x1f;
  // clean sc rays: the last build's top shell from its copy cln_top (kept in restart
  // files; clean_rot_ring: the set's slot), else the top cell's E. cln_top holds either
  // every direction per cell (older files) or n1 per packed slot
  auto fi_ = cln_top;
  const int ns = cln_ring ? 2 : 1, n1 = fi_.extent_int(4);
  const bool all = (fi_.extent_int(1) == ns*n);
  const int off = (cln_ring ? cln_rqm : 0)*(all ? n : (n + n1 - 1)/std::max(n1, 1));
  const bool cp = cln_scmem && fi_.extent_int(0) == nmb1 + 1 && fi_.extent_int(1) > 0;
  const bool e1 = cln_scmem && !cp;
  par_for("vet_vgd_mirror", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (!(xf_(m) & 4)) {return;}
    auto top = [&](const int q) {
      if (!cp) {return vi_(m,q,k+w,j+w,ie);}
      return all ? fi_(m,off+q,k+ks,j+js,ie) : fi_(m,off+q/n1,k+ks,j+js,q % n1);
    };
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    Real hp = 0.0, hm = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
      if (mr >= 0.0) {
        hp += dir_(d,3)*mr*top(d);
      } else {
        hm -= dir_(d,3)*mr*top(mp_(m,k,j,d));
      }
    }
    const Real e0 = iw_(m,VET_IW_EN,k+ks,j+js,ie);
    for (int d = 0; d < n; ++d) {   // outward ones are not read
      vi_(m,d,k+w,j+w,ie+1) = (hp > 0.0 && hm > 0.0 && !e1) ? hp/hm*top(mp_(m,k,j,d))
                                                            : e0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTopCopy
//! \brief cln_scmem: the top shell of this build's sweep into its slot of cln_top
//! (blk mod 2 on clean_rot_ring), for the next build's VetGdTopMirror

void RadiationVet::VetGdTopCopy() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, ie = indcs.ie;
  const int ks = indcs.ks, js = indcs.js, n1 = cln_s.extent_int(4), nt = (n + n1 - 1)/n1;
  const int ns = cln_ring ? 2 : 1, off = (cln_ring ? cln_rqm : 0)*nt;
  if (cln_top.extent_int(0) != nmb1 + 1 || cln_top.extent_int(1) != ns*nt) {
    Kokkos::realloc(cln_top, nmb1 + 1, ns*nt, cln_s.extent_int(2), cln_s.extent_int(3),
                    n1);
  }
  auto t_ = cln_top;
  auto vi_ = vgd_i;
  par_for("vet_vgd_top", DevExeSpace(), 0, nmb1, 0, indcs.nx3 - 1, 0, indcs.nx2 - 1, 0,
          n - 1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    t_(m,off+d/n1,k+ks,j+js,d % n1) = vi_(m,d,k+w,j+w,ie);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdInnerMirror
//! \brief reflecting inner x1 (clean scheme), after the inward pass: per column of the
//! first shell, inward radial flux / outward flux of its specular images (nearest set
//! directions, robust to round-off): the scale of the mirrored inner-face inflow
//! (1 where either vanishes)

void RadiationVet::VetGdInnerMirror() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, is = indcs.is;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  if (vgd_isc.extent_int(0) != nmb1 + 1) {Kokkos::realloc(vgd_isc, nmb1 + 1, nx3, nx2);}
  auto dir_ = vgd_dir;
  auto tg_ = vgd_trig;
  auto vi_ = vgd_i;
  auto xf_ = vgd_x1f;
  auto sc_ = vgd_isc;
  par_for("vet_vgd_imirror", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    sc_(m,k,j) = 1.0;
    if (!(xf_(m) & 1)) {return;}
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    Real hi = 0.0, ho = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
      if (mr < 0.0) {
        hi -= dir_(d,3)*mr*vi_(m,d,k+w,j+w,is);
        continue;
      }
      int dm = 0;
      Real best = -2.0;
      for (int q = 0; q < n; ++q) {
        // the nearest inward set direction: outward ones are swept only after this
        Real c = 0.0;
        for (int a = 0; a < 3; ++a) {c += dir_(q,a)*(dir_(d,a) - 2.0*mr*r[a]);}
        const Real mq = dir_(q,0)*r[0] + dir_(q,1)*r[1] + dir_(q,2)*r[2];
        if (mq < 0.0 && c > best + 1.0e-12) {best = c; dm = q;}
      }
      ho += dir_(d,3)*mr*vi_(m,dm,k+w,j+w,is);
    }
    if (hi > 0.0 && ho > 0.0) {sc_(m,k,j) = hi/ho;}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdX1Pass
//! \brief vgd_x1b(m, :) of every MeshBlock to vgd_x1r of its x1 neighbour above (up)
//! or below (on-rank copy or one message per block)

void RadiationVet::VetGdX1Pass(const bool up) {
  Mesh *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack, g0 = pmy_pack->gids;
  const int cnt = vgd_x1b.extent_int(1);
  Kokkos::fence();
  auto bh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_x1b);
  auto rh = Kokkos::create_mirror_view(vgd_x1r);
#if MPI_PARALLEL_ENABLED
  std::vector<MPI_Request> req;
#endif
  for (int m = 0; m < nmb; ++m) {
    const int gr = vgd_x1nb[2*m + (up ? 1 : 0)], gs = vgd_x1nb[2*m + (up ? 0 : 1)];
    if (gr >= 0 && pm->rank_eachmb[gr] == global_variable::my_rank) {
      for (int t = 0; t < cnt; ++t) {rh(gr - g0,t) = bh(m,t);}
    }
#if MPI_PARALLEL_ENABLED
    if (gr >= 0 && pm->rank_eachmb[gr] != global_variable::my_rank) {
      req.emplace_back();
      MPI_Isend(&bh(m,0), cnt, MPI_ATHENA_REAL, pm->rank_eachmb[gr], gr, vgd_comm,
                &req.back());
    }
    if (gs >= 0 && pm->rank_eachmb[gs] != global_variable::my_rank) {
      req.emplace_back();
      MPI_Irecv(&rh(m,0), cnt, MPI_ATHENA_REAL, pm->rank_eachmb[gs], g0 + m, vgd_comm,
                &req.back());
    }
#endif
  }
#if MPI_PARALLEL_ENABLED
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
#endif
  Kokkos::deep_copy(vgd_x1r, rh);
}

#endif  // !VET_VGD_SWEEP_TU

// VetGdSweep and its shell kernels build in their own TU, rad_vet_vetgd_sweep.cpp (this
// file with VET_VGD_SWEEP_TU: includes, helpers above and this function only)
#if defined(VET_VGD_SWEEP_TU)
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdSweep
//! \brief one sweep (inward then outward) of every direction

void RadiationVet::VetGdSweep() {
  if (vgd_csg) {   // the cubed-sphere sweep (rad_vet_vetgd_cs.cpp)
    VetGdCsSweep();
    return;
  }
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  // band index space of vgd_i / vgd_cs / the tables: active cells at [wb, wb + nx)
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vgd_scut;
  auto cs_ = vgd_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vray_cnt;
  auto mr_ = vgd_mr;
  auto x1f_ = vgd_x1f;
  // outside the block: on-rank neighbours are read from their interior (re-indexed at
  // a wall), remote ones from the band filled per shell by VetGdHalo (vgd_raw: owner or
  // band, re-indexed at the read)
  auto hl_ = vgd_hloc;
  const bool raw_ = vgd_raw;
  auto own_ = vgd_own;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const int nx2b = indcs.nx2, nx3b = indcs.nx3;
  const bool wint = vgd_wint;
  auto m3_ = vgd_m3;
  auto w3_ = vgd_w3;
  // the wall kernel (bodym) for mirror, vacuum and extrap theta walls
  const bool thm = vgd_thmir || (vgd_thvac > 0);
  const bool bfc = vet_bc_face;   // vet_bc_diff = face
  const bool imir = (ibc_x1min == VET_IBC_REFLECT) && lcut == 0;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  // inner emission at the physical inner face (icut = 0): black sphere of flux F_in
  // (vet_bc_inner, VetBcEmit: emit, or auto at a Marshak face)
  const bool emi = VetBcEmit() && (lcut == 0) && ibc_x1min != VET_IBC_REFLECT;
  const Real iemit = VetBathLo();
  // top bath (vet_bc_bath, Marshak outer x1, implicit_ebath_x1max != 0): carried by the
  // top shell's inward rays, attenuated from the top face (as vet_sc)
  const bool tbh = vgd_tbath && (ibc_x1max == VET_IBC_MARSHAK) && (iebath_x1max != 0.0);
  const Real itop = tbh ? iebath_x1max : 0.0;
  const Real mt0 = pm->mesh_size.x2min, mtl = pm->mesh_size.x2max - mt0;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const VgdLG lg_ = vgd_lg;
  const int jlo = 0, jhi = je + wb - 1;
  const int klo = 0, khi = ke + wb - 1;
  const Real twopi = 2.0*M_PI;
  // stretched theta / full sphere (sph2.inc): pol_ 2 = vet_gd_band_global (cell-window
  // clamp), 1 = block band
  const bool sth_ = vgd_sth;
  const int pol_ = vgd_pol ? (vgd_gband ? 2 : 1) : 0;
  const bool gphb_ = vgd_gband; auto phcb_ = vgd_phc; const int gpwb_ = vgd_w;
  auto thn_ = vgd_thn;
  auto hp_ = vgd_hpol;
  auto wsd_ = vgd_wsd;
  const int c2w_ = indcs.nx2 + 2*wb;
  // the fused sweep counts as two sweeps, a list-only preparation as none
  if (!vgd_xff_prep) {vgd_hsweep += vgd_xff_now ? 2 : 1;}
  // vet_gd_shell_list: per pass the (m, k, j, d) whose n . r_hat (face midpoints) lies
  // on the pass's branch, widened by 1e-12 (the kernel's test still decides); remade
  // when the set changes
  if (vgd_shl_alpha != vgd_alpha) {
    const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
    const int ntot = (nmb1 + 1)*nkj*n;
    for (int p = 0; p < 2; ++p) {
      const bool inw = (p == 0);
      if (static_cast<int>(vgd_shl[p].extent(0)) < ntot) {
        Kokkos::realloc(vgd_shl[p], ntot);
      }
      auto L_ = vgd_shl[p];
      auto mo_ = vgd_x1mo;
      int cnt = 0;
      const bool gph_ = vgd_gband; auto phc_ = vgd_phc; const int gpw_ = vgd_w;
      Kokkos::parallel_scan("vet_vgd_shl", Kokkos::RangePolicy<>(vgd_cur, 0, ntot),
      KOKKOS_LAMBDA(const int c, int &acc, const bool fin) {
        int t = c;
        const int d = t % n;
        t /= n;
        const int m = mo_(t/nkj);
        t -= (t/nkj)*nkj;
        const int k = ks + t/nj;
        const int j = js + (t % nj);
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        const Real th = sth_ ? thn_(m,j)
                             : lg_.Th(mbsize, m, j, js);
        const Real ph = (gph_ ? phc_(m, gpw_ + k - ks) : lg_.Ph(mbsize, m, k, ks));
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        const Real mr = nx*st*cp + ny*st*sp + nz*ct;
        const bool keep = inw ? (mr < 1.0e-12) : (mr >= -1.0e-12);
        if (keep) {
          if (fin) {L_(acc) = (m*nkj + t)*n + d;}
          acc += 1;
        }
      }, cnt);
      vgd_shn[p] = cnt;
      // several MeshBlocks along x1: the code range of each level
      vgd_shlv[p].assign(vgd_nb1 + 1, cnt);
      vgd_shlv[p][0] = 0;
      if (vgd_nb1 > 1) {
        auto lh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), L_);
        for (int q = cnt - 1; q >= 0; --q) {vgd_shlv[p][vgd_x1lv[lh(q)/(n*nkj)]] = q;}
        for (int lv = vgd_nb1 - 1; lv > 0; --lv) {
          vgd_shlv[p][lv] = std::min(vgd_shlv[p][lv], vgd_shlv[p][lv + 1]);
        }
      }
    }
    vgd_shl_alpha = vgd_alpha;
  }
  // per-cell sin / cos of the face-midpoint (theta, phi) for the shell kernel (made
  // once; bitwise the kernel's own evaluation)
  if (static_cast<int>(vgd_trig.extent(0)) != nmb1 + 1) {
    Kokkos::realloc(vgd_trig, nmb1 + 1, ke - ks + 1, je - js + 1, 4);
    auto tw_ = vgd_trig;
    const bool gph_ = vgd_gband; auto phc_ = vgd_phc; const int gpw_ = vgd_w;
    par_for("vet_vgd_trig", vgd_cur, 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real th = sth_ ? thn_(m,j)
                           : lg_.Th(mbsize, m, j, js);
      const Real ph = (gph_ ? phc_(m, gpw_ + k - ks) : lg_.Ph(mbsize, m, k, ks));
      tw_(m,k-ks,j-js,0) = sin(th);
      tw_(m,k-ks,j-js,1) = cos(th);
      tw_(m,k-ks,j-js,2) = sin(ph);
      tw_(m,k-ks,j-js,3) = cos(ph);
    });
  }
  auto tg_ = vgd_trig;
  // the exact lists of every (pass, shell) whose record does not match (new set, cut or
  // depth), in a few batched collectives (VetGdHxBatch)
  if (vgd_ipc_on && vgd_ipc_state == 0 && vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0)) {
    const size_t cs = std::max<size_t>(vgd_crb.extent(0), 1);
    VetGdIpcInit(vgd_prk, cs, 2*cs);   // vet_gd_halo_ipc (collective, once)
  }
  if (vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0)) {
    const int n1c = indcs.nx1 + 2*indcs.ng;
    if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
    std::vector<std::array<int, 4>> jobs;
    const int nshb = n1 - lcut;
    for (int pass = 0; pass < 2; ++pass) {
      const bool inw = (pass == 0);
      for (int q = 0; q < nshb; ++q) {
        const int l = inw ? (n1 - 1 - q) : (lcut + q);
        const int i = is + l;
        const int wsx = inw ? vgd_wsi[i] : vgd_wso[i];
        if (wsx < 1 || wsx > vgd_w) {continue;}
        if (!VetGdHxKey(pass*n1c + i, wsx)) {jobs.push_back({pass, i, l, wsx});}
      }
    }
    VetGdHxBatch(jobs);
  }
  // a list-only preparation of this set (shell lists, exact halo lists)
  if (vgd_xff_prep) {return;}
  const bool xffs = vgd_xff_now;
  if (ibc_x1max == VET_IBC_REFLECT) {VetGdTopMirror();}
  // vet_gd_band_ring: shells are lent a wide slab while their band is needed (the
  // other set's array as the second one, VetGdRingFlush)
  VetGdRingSweepBegin(xffs);
  const bool rga = vgd_ring_act;
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    if (imir && !inw) {VetGdInnerMirror();}
    auto isc_ = vgd_isc;
#define VET_VGD_BODY_NAME body
#define VET_VGD_BODY_MIR 0
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#define VET_VGD_BODY_NAME bodym
#define VET_VGD_BODY_MIR 1
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#if defined(KOKKOS_ENABLE_CUDA)
    // CUDA gx-false kernels (rad_vet_vetgd_body.inc): bodyg / bodymg compile in
    // vpa = gx = false (nvcc has no generic device lambda for std::false_type)
#define VET_VGD_BODY_GX0
#define VET_VGD_BODY_NAME bodyg
#define VET_VGD_BODY_MIR 0
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#define VET_VGD_BODY_NAME bodymg
#define VET_VGD_BODY_MIR 1
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
#undef VET_VGD_BODY_GX0
#endif
    // shell bodies of the fused two-set kernel (run_xff), one per set, with its tables
    // and intensities bound to the body's names. Both come from one lambda expression
    // (one closure type), so the kernel holds one copy of the body and stack frame
#if defined(KOKKOS_ENABLE_CUDA)
#define VET_VGD_BODY_GX0
#endif
    using VgdDirV = decltype(vgd_dir);
    using VgdWallV = decltype(vgd_wall);
    using VgdMapV = decltype(vgd_map);
    using VgdMrV = decltype(vgd_mr);
    using VgdM3V = decltype(vgd_m3);
    using VgdW3V = decltype(vgd_w3);
    using VgdIV = decltype(vgd_i);
    auto mkx0 = [&](const VgdDirV dir_, const VgdWallV wl_, const VgdMapV mp_,
                    const VgdMrV mr_, const VgdM3V m3_, const VgdW3V w3_,
                    const VgdIV vi_) {
#define VET_VGD_BODY_NAME bodyx
#define VET_VGD_BODY_MIR 0
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
      return bodyx;
    };
    auto mkx1 = [&](const VgdDirV dir_, const VgdWallV wl_, const VgdMapV mp_,
                    const VgdMrV mr_, const VgdM3V m3_, const VgdW3V w3_,
                    const VgdIV vi_) {
#define VET_VGD_BODY_NAME bodyx
#define VET_VGD_BODY_MIR 1
#include "rad_vet_vetgd_body.inc"
#undef VET_VGD_BODY_NAME
#undef VET_VGD_BODY_MIR
      return bodyx;
    };
#if defined(KOKKOS_ENABLE_CUDA)
#undef VET_VGD_BODY_GX0
#endif
    // vet_gd_body_global: this pass's closure copied once to device memory; launches
    // carry only its address (bitwise)
    using VgdBodyT = decltype(body);
    using VgdBodymT = decltype(bodym);
    const VgdBodyT *bdp = nullptr;
    const VgdBodymT *bmp = nullptr;
    if (thm) {
      bmp = VgdClosureToDevice(bodym, vgd_bbuf1);
    } else {
      bdp = VgdClosureToDevice(body, vgd_bbuf0);
    }
#if defined(KOKKOS_ENABLE_CUDA)
    using VgdBodygT = decltype(bodyg);
    using VgdBodymgT = decltype(bodymg);
    const VgdBodygT *bgp = nullptr;
    const VgdBodymgT *bgq = nullptr;
    if (bmp != nullptr) {bgq = VgdClosureToDevice(bodymg, vgd_bbuf3);}
    if (bdp != nullptr) {bgp = VgdClosureToDevice(bodyg, vgd_bbuf2);}
#endif
    // both sets' gx-false closures in device memory (VetGdXffOk)
    using VgdXb0T = decltype(mkx0(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                  vgd_i));
    using VgdXb1T = decltype(mkx1(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                  vgd_i));
    const VgdXb0T *xa0 = nullptr, *xq0 = nullptr;
    const VgdXb1T *xa1 = nullptr, *xq1 = nullptr;
    if (xffs) {
#if defined(KOKKOS_ENABLE_CUDA)
      // the current set's closure: the buffer of its gx-false body
      auto &bba = vgd_bbuf2;
#else
      auto &bba = vgd_bbuf0;
#endif
      if (thm) {
        xa1 = VgdClosureToDevice(mkx1(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                      vgd_i), bba);
        xq1 = VgdClosureToDevice(mkx1(vgd_tabd.dir, vgd_tabd.wall, vgd_tabd.map,
                                      vgd_tabd.mr, vgd_tabd.m3, vgd_tabd.w3, vgd_i_xb),
                                 vgd_bbufx);
      } else {
        xa0 = VgdClosureToDevice(mkx0(vgd_dir, vgd_wall, vgd_map, vgd_mr, vgd_m3, vgd_w3,
                                      vgd_i), bba);
        xq0 = VgdClosureToDevice(mkx0(vgd_tabd.dir, vgd_tabd.wall, vgd_tabd.map,
                                      vgd_tabd.mr, vgd_tabd.m3, vgd_tabd.w3, vgd_i_xb),
                                 vgd_bbufx);
      }
    }
    const int nkj = (ke - ks + 1)*(je - js + 1), nj = je - js + 1;
    // fused shell kernel of (i, l): this pass's codes of the current set, then the
    // other set's, one launch; each thread runs its code's set (bitwise two launches)
    auto run_xff = [&](const int i, const int l) {
      const int nA = vgd_shn[pass], nB = vgd_shn_xb[pass];
      if (nA + nB <= 0) {return;}
      auto LA_ = vgd_shl[pass];
      auto LB_ = vgd_shl_xb[pass];
      const std::false_type f_;
      // current set's codes padded to a multiple of 1024: every block (power of two
      // <= 1024) runs one set only
      const int nP = ((nA + 1023)/1024)*1024;
      if (thm) {
        const auto *qa = xa1;
        const auto *qb = xq1;
        Kokkos::parallel_for("vet_vgd_shell_xff", VgdShellPol(vgd_cur, 0, nP + nB),
        KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nP);
          if (!sb && c >= nA) {return;}
          int t = sb ? LB_(c - nP) : LA_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*(sb ? qb : qa))(m, k, j, d, i, l, f_, f_);   // one body: one call site
        });
      } else {
        const auto *qa = xa0;
        const auto *qb = xq0;
        Kokkos::parallel_for("vet_vgd_shell_xff", VgdShellPol(vgd_cur, 0, nP + nB),
        KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nP);
          if (!sb && c >= nA) {return;}
          int t = sb ? LB_(c - nP) : LA_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*(sb ? qb : qa))(m, k, j, d, i, l, f_, f_);   // one body: one call site
        });
      }
    };
    // the shell kernel over codes cl(c), c in [c0, c1) (vet_gd_shell_list)
    auto run_list = [&](const DvceArray1D<int> &cl, const int c0, const int c1,
                        const int i, const int l, const char *name) {
      if (c1 <= c0) {return;}
      auto L_ = cl;
      // vet_gd_body_global: same codes, closure through a pointer; one kernel per
      // variant
      const std::false_type f_;
#if defined(KOKKOS_ENABLE_CUDA)
      // on CUDA the gx-false closures bodyg / bodymg
      const VgdBodygT *bp_ = bgp;
      const VgdBodymgT *bq_ = bgq;
#else
      const VgdBodyT *bp_ = bdp;
      const VgdBodymT *bq_ = bmp;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
      // launch-bounds variants of the CUDA gx-false kernel (VET_VGD_CUDA_LB; 0 = the
      // VgdShellPol launch): register allocation and block size only (same bits)
#define VET_VGD_GF_CU(...)                                                    \
      if (bq_ != nullptr) {                                               \
        Kokkos::parallel_for(name, __VA_ARGS__, KOKKOS_LAMBDA(const int c) { \
          int t = L_(c);                                                  \
          const int d = t % n;                                            \
          t /= n;                                                         \
          const int m = t/nkj;                                            \
          t -= m*nkj;                                                     \
          const int k = ks + t/nj;                                        \
          const int j = js + (t % nj);                                    \
          (*bq_)(m, k, j, d, i, l, f_, f_);                               \
        });                                                               \
      } else {                                                            \
        Kokkos::parallel_for(name, __VA_ARGS__, KOKKOS_LAMBDA(const int c) { \
          int t = L_(c);                                                  \
          const int d = t % n;                                            \
          t /= n;                                                         \
          const int m = t/nkj;                                            \
          t -= m*nkj;                                                     \
          const int k = ks + t/nj;                                        \
          const int j = js + (t % nj);                                    \
          (*bp_)(m, k, j, d, i, l, f_, f_);                               \
        });                                                               \
      }
      switch (VET_VGD_CUDA_LB) {
        case 1: VET_VGD_GF_CU(VgdLbPol<128, 4>(vgd_cur, c0, c1)); break;
        case 2: VET_VGD_GF_CU(VgdLbPol<128, 5>(vgd_cur, c0, c1)); break;
        case 3: VET_VGD_GF_CU(VgdLbPol<128, 6>(vgd_cur, c0, c1)); break;
        case 4: VET_VGD_GF_CU(VgdLbPol<256, 2>(vgd_cur, c0, c1)); break;
        case 5: VET_VGD_GF_CU(VgdLbPol<256, 3>(vgd_cur, c0, c1)); break;
        default: VET_VGD_GF_CU(VgdShellPol(vgd_cur, c0, c1)); break;
      }
#undef VET_VGD_GF_CU
#else
      if (bq_ != nullptr) {
        Kokkos::parallel_for(name, VgdShellPol(vgd_cur, c0, c1),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*bq_)(m, k, j, d, i, l, f_, f_);
        });
      } else {
        Kokkos::parallel_for(name, VgdShellPol(vgd_cur, c0, c1),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          (*bp_)(m, k, j, d, i, l, f_, f_);
        });
      }
#endif
    };
    const int nsh = n1 - lcut;
    // vet_gd_halo_exact: read marks of shell i's band in this pass (when its record
    // does not match; same decision on every rank): the ghosts read by inward shell i-1
    // and outward shell i's turned rays (inward exchange) or outward shell i+1 (outward)
    const int c2b = nx2b + 2*wb, c3b = nx3b + 2*wb;
    auto hx_mark = [&](const int i, const int l, const int wsx) -> bool {
      const int n1c = indcs.nx1 + 2*indcs.ng;
      if (static_cast<int>(vgd_hxe.size()) != 2*n1c) {vgd_hxe.resize(2*n1c);}
      if (VetGdHxKey(pass*n1c + i, wsx)) {return false;}
      const size_t nall = static_cast<size_t>(nmb1 + 1)*n*c3b*c2b;
      // stamped marks (no clearing): a ghost is marked when it holds the build's stamp
      if (vgd_hx_mk.extent(0) < nall) {
        Kokkos::realloc(vgd_hx_mk, nall);
        vgd_hx_stamp = 0;
      }
      if (!vgd_hx_oc.is_allocated()) {
        vgd_hx_oc = Kokkos::View<int, DevMemSpace>("vet_vgd_hx_oc");
      }
      vgd_hx_stamp += 1;
      const int stamp = vgd_hx_stamp;
      Kokkos::deep_copy(vgd_cur, vgd_hx_oc, 0);
      auto oc_ = vgd_hx_oc;
      auto mk_ = vgd_hx_mk;
      auto mark = [&](const int pr, const int ir, const int lr) {
        const bool inwr = (pr == 0);
        auto L_ = vgd_shl[pr];
        Kokkos::parallel_for("vet_vgd_hx_mark", Kokkos::RangePolicy<>(vgd_cur, 0,
                                                                      vgd_shn[pr]),
        KOKKOS_LAMBDA(const int c) {
          int t = L_(c);
          const int d = t % n;
          t /= n;
          const int m = t/nkj;
          t -= m*nkj;
          const int k = ks + t/nj;
          const int j = js + (t % nj);
          VgdCand(m, k, j, d, ir, lr, inwr, n1, lcut, js, ks, jlo, jhi, klo, khi, dir_,
                  mbsize, lg_, cx1v, cx1f, sth_, pol_, thn_, c2w_, imir,
          [&](const int iu, const int ka, const int kb, const int ja, const int jb,
              const int dr) {
            if (iu != i) {return;}
            for (int kk = ka; kk <= kb; ++kk) {
              for (int jj = ja; jj <= jb; ++jj) {
                const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
                const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
                if (ek == 0 && ej == 0) {continue;}
                if (ej != 0 && wl_(m,kk,jj) >= 2) {
                  // the main sweep reads the mirror row (flag 2) or the edge row
                  // (extrap, 4) through the phi-edge slot; vacuum (3) reads nothing.
                  // Flags 3, 4: the twin's periodic reads below too
                  const int wf = wl_(m,kk,jj);
                  const int o2 = 3*(ek + 1) + 1;
                  if (wf != 3 && ek != 0 && hl_(8*m + ((o2 < 4) ? o2 : (o2 - 1))) < 0) {
                    const bool hi = (ej > 0);
                    const int g = hi ? (jj - wb - nx2b) : (wb - 1 - jj);
                    const int jm = (wf == 4) ? (hi ? (wb + nx2b - 1) : wb)
                                             : VgdMirRow(g, hi, wb, nx2b);
                    for (int t = 0; t < ((wint && wf == 2) ? 3 : 1); ++t) {
                      if (t > 0 && w3_(m,kk,jj,dr,t) == 0.0) {continue;}
                      const int dq = (wf == 4) ? dr : (wint ? m3_(m,kk,jj,dr,t)
                                                           : mp_(m,kk,jj,dr));
                      mk_(((m*n + dq)*c3b + kk)*c2b + jm) = stamp;
                    }
                    if (kk < wb - wsx || kk >= wb + nx3b + wsx) {
                      Kokkos::atomic_add(&oc_(), 1);
                    }
                  }
                  if (wf == 2) {continue;}
                }
                const int oo = 3*(ek + 1) + (ej + 1);
                if (hl_(8*m + ((oo < 4) ? oo : (oo - 1))) >= 0) {continue;}
                mk_(((m*n + dr)*c3b + kk)*c2b + jj) = stamp;
                if (jj < wb - wsx || jj >= wb + nx2b + wsx || kk < wb - wsx ||
                    kk >= wb + nx3b + wsx) {
                  Kokkos::atomic_add(&oc_(), 1);
                }
              }
            }
          });
        });
      };
      if (inw) {
        if (l - 1 >= lcut) {mark(0, i - 1, l - 1);}
        mark(1, i, l);
      } else if (l + 1 <= n1 - 1) {
        mark(1, i + 1, l + 1);
      }
      return true;
    };
    const bool hxm = vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0);
    // several MeshBlocks along x1: levels in the pass's order, each from the boundary
    // shell handed on (wall maps couple all directions: no wavefront)
    const int nst = vgd_nb1;
    for (int qq = 0; qq < nst*nsh; ++qq) {
      const int q = qq % nsh, lv = inw ? (nst - 1 - qq/nsh) : qq/nsh;
      if (nst > 1 && q == 0 && qq > 0) {VetGdX1Hand(!inw);}
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      // the shell exchanged after this one (its mask is prepared while this shell's
      // messages travel)
      if (q + 1 < n1 - lcut) {
        const int in = is + (inw ? (n1 - 2 - q) : (lcut + q + 1));
        vgd_hnext[0] = in;
        vgd_hnext[1] = inw ? 1 : 0;
        vgd_hnext[2] = inw ? vgd_wsi[in] : vgd_wso[in];
      } else if (inw) {
        const int in = is + lcut;
        vgd_hnext[0] = in;
        vgd_hnext[1] = 0;
        vgd_hnext[2] = vgd_wso[in];
      } else {
        vgd_hnext[0] = -1;
      }
      if (rga) {
        // this shell's wide slab (queued during the previous exchange, else now):
        // inward from scratch, outward with its interior and inward band
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
      if (xffs) {
        // both cross-fade sets, one launch
        run_xff(i, l);
      } else {
        // shell list: only this pass's (and level's) (m, k, j, d)
        run_list(vgd_shl[pass], vgd_shlv[pass][lv], vgd_shlv[pass][lv + 1], i, l,
                 "vet_vgd_shell");
      }
      if (rga) {
        // the previous shell has been read (stream order): save its inward band
        // (inward pass), its interior back to the compact slab, lend the next shell a
        // slab; in this shell's exchange hook (VetGdHlBegin) or before the next kernel
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
      const int wsx = inw ? vgd_wsi[i] : vgd_wso[i];
      // the shell is complete on every block: its lateral band (exact, not lagged)
      Kokkos::Timer th;
      vgd_hinw = inw;
      if (hxm && hx_mark(i, l, wsx)) {
        // vet_gd_halo_exact = 2: lists from the marks alone, before the exchange
        VetGdHxBuild2(pass*(indcs.nx1 + 2*indcs.ng) + i, i, wsx);
      }
      if (xffs) {
        // one message per partner with both sets' exact lists (vgd_raw: VetGdGbHaloX);
        // else two exchanges, the other set's first (each writes only its own band)
        vgd_rg_lastx = false;
        const bool hx2 = vgd_raw ? VetGdGbHaloX(i, wsx)
                                 : VetGdHlXff(i, wsx);
        if (!hx2 && rga && !vgd_raw) {
          VgdFatal("vet_gd_rotate_xfade_fuse: the band ring needs the merged exchange");
        }
        if (rga && inw) {   // record of the saved inward band, as below; both sets
          vgd_rg_sk[i] = vgd_raw ? 3 : (vgd_rg_lastx ? 1 : 2);
          vgd_rg_se[i] = vgd_raw ? VetGdRawKey(i, true, wsx) : (vgd_rg_lastx ? i : -1);
        }
        if (!hx2) {
          VgdXffSwap();
          VetGdHalo(vgd_i, n, i, i, wsx, true);
          VgdXffSwap();
          VetGdHalo(vgd_i, n, i, i, wsx, true);
        }
      } else {
        vgd_rg_lastx = false;
        VetGdHalo(vgd_i, n, i, i, wsx, true);
        if (rga && inw) {
          // record of the values to save from this shell's inward band: the exchange's
          // exact list (vgd_raw: what it wrote), else the whole band
          const int ex = i;
          const bool xl = vgd_rg_lastx && (vgd_rg_lastx_e == ex);
          vgd_rg_sk[i] = vgd_raw ? 3 : (xl ? 1 : 2);
          vgd_rg_se[i] = vgd_raw ? VetGdRawKey(i, true, wsx) : (xl ? ex : -1);
        }
      }
      if (rga) {VetGdRingFlush();}
      // the next kernel is stream-ordered behind the unpack: no fence needed
      vgd_thalo += th.seconds();
    }
    if (rga) {
      // last shell of the pass: saved (inward) and back to its compact slab
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
}
#else

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdSource
//! \brief source of a gd build: ln chi, ln S of every cell (VetRaySource) into the band
//! index space of vgd_cs with its lateral band filled exactly (one multi-shell
//! exchange); with several MeshBlocks along x1 also the x1 ghost shells

void RadiationVet::VetGdSource() {
  VetRaySource();
  // from the first shell up (gd sweep), or every shell for the clean scheme, whose rays
  // and closure read the whole column
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int og = vgd_w - indcs.ng, ok = vgd_w - indcs.ks, x1g = (vgd_nb1 > 1) ? 1 : 0;
  const int cut = 0;
  const int ilo = indcs.is + cut - x1g, ie = indcs.ie + x1g;
  auto cs_ = vray_cs;
  auto cw_ = vgd_cs;
  par_for("vet_vgd_csw", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks,
          indcs.ke, indcs.js, indcs.je, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    cw_(m,0,k+ok,j+og,i) = cs_(m,0,k,j,i);
    cw_(m,1,k+ok,j+og,i) = cs_(m,1,k,j,i);
  });
  vgd_cur = DevExeSpace();
  if (vgd_csg) {   // the cs band copies
    VetGdCsSource();
    return;
  }
  VetGdHaloCs(vgd_cs, ilo, ie);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdBuild
//! \brief the gd build: source and first shell (VetRaySource),
//! sweep(s) with the lagged lateral inflow, moments; timed per part (fenced)

void RadiationVet::VetGdBuild() {
  Kokkos::Timer tm;
  VetGdSource();
  Kokkos::fence();
  vgd_tsrc += tm.seconds();
  VetGdTabPoll();      // the prefetched tables, uploaded when ready
  // vet_gd_rotate_every: a new z-angle at the first build of a block of N cycles; the
  // stored intensities belong to the old directions, so the lagged inflow re-converges
  const Real al = VetGdAngle(pmy_pack->pmesh->ncycle);
  // no vet_rotate_avg or cross-fade here: the clean scheme sweeps rotating sets on its
  // ring (CleanRingMoments); without clean_rot_ring they are refused (CleanInit)
  if (al != vgd_alpha) {
    Kokkos::Timer tt;
    VetGdTables(al);
    vgd_ttab += tt.seconds();
  }
  VetGdTabPrefetch(pmy_pack->pmesh->ncycle);
  VetGdSweepPost();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdSweepPost
//! \brief sweep part of a gd build with the current tables, and VetGdPost

void RadiationVet::VetGdSweepPost() {
  Kokkos::Timer tm;
  // the per-shell halo makes one sweep exact
  vgd_cur = DevExeSpace();
  vgd_scut = vray_icut;
  tm.reset();
  VetGdSweep();
  Kokkos::fence();
  vgd_tswp += tm.seconds();
  VetGdPost();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdPost
//! \brief after the sweep(s) of a build (also clean_rot_ring's): the GPU halo's error
//! flag, the clamp check over all ranks

void RadiationVet::VetGdPost() {
  Kokkos::fence();
  VetGdIpcCheck();
  vray_ncall += 1.0;
  // persistent pinned buffer
  if (vgd_pin_cnt.extent(0) != vray_cnt.extent(0)) {
    vgd_pin_cnt = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("vet_vgd_pin_cnt",
                                                                     vray_cnt.extent(0));
  }
  Kokkos::deep_copy(vgd_pin_cnt, vray_cnt);
  vray_nclamp = vgd_pin_cnt(0);
  // all ranks (sweep exactness: a read beyond the band is clamped and counted)
  vgd_nclamp_all = vray_nclamp;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &vgd_nclamp_all, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  // the band is capped at the block size: a read beyond it is clamped (inexact).
  // Fatal unless vet_gd_allow_clamp (warning); on the cubed sphere every clamped read
  // counts (the depth there is a bound with margin, not an exact reach)
  if ((vgd_capped || vgd_csg) && vgd_nclamp_all > vgd_nclamp_seen) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd: " << (vgd_nclamp_all - vgd_nclamp_seen)
                << " lateral reads were CLAMPED to the band (all ranks) at cycle "
                << pmy_pack->pmesh->ncycle << ": the MeshBlocks are smaller than the "
                << "near-tangent reach, the sweep is not exact" << std::endl;
    }
    if (!vgd_allow_clamp) {
      VgdFatal("clamped lateral reads (MeshBlocks too small for the reach); use larger "
               "MeshBlocks or <rad_vet>/vet_gd_allow_clamp = true");
    }
    vgd_nclamp_seen = vgd_nclamp_all;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdRagAlloc
//! \brief the ragged intensity array (VgdRag, rad_vet.hpp): shell i gets band depth
//! vgd_wsh[i] (its data's deepest read), ghost shells none; zero-filled

void RadiationVet::VgdRagAlloc(VgdRag &a) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  a.nmb = pmy_pack->nmb_thispack;
  a.n = vgd_n;
  a.nx2 = indcs.nx2;
  a.nx3 = indcs.nx3;
  a.w = vgd_w;
  a.c1 = c1;
  Kokkos::realloc(a.off, a.nmb, c1);
  Kokkos::realloc(a.wi, c1);
  Kokkos::realloc(a.sb, a.nmb, 4);
  Kokkos::realloc(a.oos, 1);
  auto off_h = Kokkos::create_mirror_view(a.off);
  auto wi_h = Kokkos::create_mirror_view(a.wi);
  auto sb_h = Kokkos::create_mirror_view(a.sb);
  // a side carries the band when any of its three neighbours (edge, two corners) is
  // not on this rank (vgd_hloc < 0)
  auto hl_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_hloc);
  for (int m = 0; m < a.nmb; ++m) {
    int side[4] = {0, 0, 0, 0};
    for (int o = 0; o < 8; ++o) {
      if (hl_h(8*m + o) >= 0) {continue;}
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      if (dk < 0) {side[0] = 1;}
      if (dk > 0) {side[1] = 1;}
      if (dj < 0) {side[2] = 1;}
      if (dj > 0) {side[3] = 1;}
    }
    for (int q = 0; q < 4; ++q) {sb_h(m,q) = side[q];}
  }
  if (vgd_raw) {
    // raw band (every slot remote): a side carries it only where a band position has
    // no on-rank owner (vgd_own < 0); on-rank owners are read directly
    auto ow_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_own);
    const int w = vgd_w;
    for (int m = 0; m < a.nmb; ++m) {
      for (int q = 0; q < 4; ++q) {sb_h(m,q) = 0;}
      for (int kk = 0; kk < a.nx3 + 2*w; ++kk) {
        for (int jj = 0; jj < a.nx2 + 2*w; ++jj) {
          if (ow_h(m,kk,jj,0) >= 0) {continue;}
          sb_h(m,0) |= (kk < w);
          sb_h(m,1) |= (kk >= w + a.nx3);
          sb_h(m,2) |= (jj < w);
          sb_h(m,3) |= (jj >= w + a.nx2);
        }
      }
    }
  }
  for (int i = 0; i < c1; ++i) {
    int wl = 0;
    if ((i >= indcs.is && i <= indcs.ie) || vgd_nb1 > 1) {
      wl = std::max(vgd_wsh[i], std::max(vgd_wsi[i], vgd_wso[i]));
      wl = std::min(wl, vgd_w);
    }
    wi_h(i) = wl;
  }
  int64_t tot = 0, tth = 0;
  const int cg = (vgd_rag_cut > 0) ? indcs.is + vgd_rag_cut : 0;
  if (vgd_ring) {
    // vet_gd_band_ring: every shell compact (depth 0) at [0, ...), then vgd_rg_k wide
    // slots of one shell each at the deepest band; shell i's full depth in
    // vgd_rg_wfull (same for both arrays)
    vgd_rg_wfull.assign(c1, 0);
    for (int i = 0; i < c1; ++i) {vgd_rg_wfull[i] = wi_h(i);}
    const int64_t cap = static_cast<int64_t>(a.nx3)*a.nx2*a.n;
    for (int m = 0; m < a.nmb; ++m) {
      for (int i = 0; i < c1; ++i) {
        off_h(m,i) = (i < cg) ? -1 - tth : tot;   // (below cg: VgdRagThick)
        ((i < cg) ? tth : tot) += cap;
      }
    }
    const int64_t ncomp = tot;
    vgd_rg_k = VGD_RING_K;
    vgd_rg_soff.assign(static_cast<size_t>(vgd_rg_k)*a.nmb, 0);
    for (int s = 0; s < vgd_rg_k; ++s) {
      for (int m = 0; m < a.nmb; ++m) {
        int64_t mx = 0;
        for (int i = 0; i < c1; ++i) {
          const int wl = wi_h(i);
          mx = std::max(mx, static_cast<int64_t>(a.n)*
                            (a.nx3 + (sb_h(m,0) + sb_h(m,1))*wl)*
                            (a.nx2 + (sb_h(m,2) + sb_h(m,3))*wl));
        }
        vgd_rg_soff[static_cast<size_t>(s)*a.nmb + m] = tot;
        tot += mx;
      }
    }
    for (int i = 0; i < c1; ++i) {wi_h(i) = 0;}
    vgd_rg_state.assign(c1, -1);
    vgd_rg_sown.assign(vgd_rg_k, -1);
    vgd_rg_so.assign(c1, 0);
    vgd_rg_sn.assign(c1, 0);
    vgd_rg_sk.assign(c1, 0);
    vgd_rg_se.assign(c1, -1);
    vgd_rg_pend.clear();
    Kokkos::realloc(a.offc, a.nmb, c1);
    Kokkos::realloc(vgd_rg_dsoff, vgd_rg_k, a.nmb);
    auto so_h = Kokkos::create_mirror_view(vgd_rg_dsoff);
    for (int m = 0; m < a.nmb; ++m) {
      for (int s = 0; s < vgd_rg_k; ++s) {
        so_h(s,m) = vgd_rg_soff[static_cast<size_t>(s)*a.nmb + m];
      }
    }
    Kokkos::deep_copy(a.offc, off_h);
    Kokkos::deep_copy(vgd_rg_dsoff, so_h);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_band_ring: compact intensities "
                << ncomp*sizeof(Real)/1.0e9 << " GB + "
                << vgd_rg_k << " band slots "
                << (tot - ncomp)*sizeof(Real)/1.0e9 << " GB on rank 0" << std::endl;
    }
  } else {
  for (int m = 0; m < a.nmb; ++m) {
    for (int i = 0; i < c1; ++i) {
      const int wl = wi_h(i);
      off_h(m,i) = (i < cg) ? -1 - tth : tot;
      ((i < cg) ? tth : tot) += static_cast<int64_t>(a.n)*
             (a.nx3 + (sb_h(m,0) + sb_h(m,1))*wl)*(a.nx2 + (sb_h(m,2) + sb_h(m,3))*wl);
    }
  }
  }
  a.dummy = tot;
  a.nth = tth;
  Kokkos::deep_copy(a.off, off_h);
  Kokkos::deep_copy(a.wi, wi_h);
  Kokkos::deep_copy(a.sb, sb_h);
  const int64_t nd = tot + 1;
  if (!MemArenaVgd(a, nd)) {Kokkos::realloc(a.d, nd);}
  Kokkos::deep_copy(a.d, 0.0);
  Kokkos::deep_copy(a.oos, 0.0);
  if (global_variable::my_rank == 0) {
    const double dense = static_cast<double>(a.nmb)*a.n*(a.nx3 + 2*vgd_w)*
                         (a.nx2 + 2*vgd_w)*c1;
    std::cout << "<rad_vet> vet_gd: ragged per-shell band intensity array "
              << nd*sizeof(Real)/1.0e9 << " GB on rank 0 (dense band "
              << dense*sizeof(Real)/1.0e9 << " GB)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VgdRagThick
//! \brief the slabs of the shells below vgd_rag_cut, one set of them for both fused
//! arrays, behind both in the arena: the arrays' offsets there (their own bases)

void RadiationVet::VgdRagThick() {
  if (vgd_i.nth == 0) return;
  MemArenaGrow(marena_sw + vgd_i.nth, false);   // (both arrays hold zeros)
  Real *t = marena.data() + marena_sw;
  marena_sw += vgd_i.nth + ((64 - (vgd_i.nth % 64)) % 64);
  for (VgdRag *a : {&vgd_i, &vgd_i_xb}) {
    const int64_t b = t - a->d.data();
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), a->off);
    for (size_t m = 0; m < h.extent(0); ++m) {
      for (size_t i = 0; i < h.extent(1); ++i) {
        if (h(m,i) < 0) h(m,i) = b + (-1 - h(m,i));
      }
    }
    Kokkos::deep_copy(a->off, h);
    if (vgd_ring) Kokkos::deep_copy(a->offc, h);
    a->d = DvceArray1D<Real>(a->d.data(), b + a->nth);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: the shells below " << vgd_rag_cut
              << " (ray cut " << cln_rcut << ") share one slab, "
              << vgd_i.nth*sizeof(Real)/1.0e9 << " GB on rank 0" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdRingSweepBegin / VetGdRingSweepEnd / VetGdRingFlush
//! \brief vet_gd_band_ring (rad_vet.hpp): per-sweep bookkeeping and the kernel of a
//! group of slab operations (save a shell's inward band; release = interior back to
//! the compact slab; acquire = interior in and saved band back, then both shells'
//! tables). Pure copies: the sweep reads the values of the dense-band layout.

// band_ring copy kernel (one thread per value; the last thread sets the tables)
struct VgdRingK {
  Real *da, *db, *st;
  int64_t *aoff, *boff;
  int *awi, *bwi;
  const int64_t *oc, *ocb, *so;
  const int *sb, *svdl, *svdw, *rsdl, *rsdw;
  int64_t offs0, offs1, ocs0, ocs1, sos0, sos1;
  int64_t nsv, nrl, ncp, nrs, svso, rsso, svn, rsn;
  int sbs0, sbs1, n, nx2, nx3, w, c2b, c3b, nm, svk, rsk, svnr, rsnr;
  int svs, svw, rli, rls, rlw, aqi, aqs, aqw;
  bool two, hrl, haq;
  // xb: db is the other cross-fade set's array; its saved / restored values follow its
  // own records (lists 2, counts svn2 / rsn2, after the first set's)
  const int *svdl2, *svdw2, *rsdl2, *rsdw2;
  int64_t svn2, rsn2;
  int svnr2, rsnr2;
  bool xb;
  // value of the wide slab of slot s (depth wl) of block m at band index (kb, jb);
  // -1 outside the shell's rectangle
  KOKKOS_INLINE_FUNCTION
  int64_t Wadr(const int s, const int wl, const int m, const int v, const int kb,
               const int jb) const {
    const int bk0 = sb[m*sbs0]*wl, bj0 = sb[m*sbs0 + 2*sbs1]*wl;
    const int c3 = nx3 + bk0 + sb[m*sbs0 + sbs1]*wl;
    const int c2 = nx2 + bj0 + sb[m*sbs0 + 3*sbs1]*wl;
    const int kk = kb - (w - bk0), jj = jb - (w - bj0);
    if (kk < 0 || kk >= c3 || jj < 0 || jj >= c2) {return -1;}
    return so[s*sos0 + m*sos1] + v +
           static_cast<int64_t>(n)*(kk + static_cast<int64_t>(c3)*jj);
  }
  // the t-th saved value of a record: its band index (m, v, kb, jb); false: none
  KOKKOS_INLINE_FUNCTION
  bool Rec(const int k, const int wl, const int *dl, const int *dw, const int nr,
           const int64_t t, int &m, int &v, int &kb, int &jb) const {
    if (k == 1) {
      int a = (t < nr) ? dl[t] : dw[t - nr];
      jb = a % c2b;
      a /= c2b;
      kb = a % c3b;
      a /= c3b;
      v = a % n;
      m = a/n;
      return true;
    }
    // dense: (m, v, kb, jb) over the band index rectangle of depth wl, interior out
    const int64_t c3w = nx3 + 2*wl, c2w = nx2 + 2*wl;
    int64_t r = t;
    v = static_cast<int>(r % n);
    r /= n;
    const int kw = static_cast<int>(r % c3w);
    r /= c3w;
    const int jw = static_cast<int>(r % c2w);
    m = static_cast<int>(r/c2w);
    kb = kw + (w - wl);
    jb = jw + (w - wl);
    return !(kb >= w && kb < w + nx3 && jb >= w && jb < w + nx2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t tt) const {
    int64_t t = tt;
    if (xb && t < nsv) {
      const bool b2 = (t >= svn);
      const int64_t tl = b2 ? (t - svn) : t;
      int m, v, kb, jb;
      if (!Rec(svk, svw, b2 ? svdl2 : svdl, b2 ? svdw2 : svdw, b2 ? svnr2 : svnr, tl,
               m, v, kb, jb)) {return;}
      const int64_t q = Wadr(svs, svw, m, v, kb, jb);
      if (q < 0) {return;}
      st[svso + t] = b2 ? db[q] : da[q];
      return;
    }
    if (t < nsv) {
      int m, v, kb, jb;
      if (!Rec(svk, svw, svdl, svdw, svnr, t, m, v, kb, jb)) {return;}
      const int64_t q = Wadr(svs, svw, m, v, kb, jb);
      if (q < 0) {return;}
      st[svso + t] = da[q];
      if (two) {st[svso + svn + t] = db[q];}
      return;
    }
    t -= nsv;
    if (t < nrl + ncp) {
      const bool out = (t < nrl);
      if (!out) {t -= nrl;}
      const int i = out ? rli : aqi;
      const int s = out ? rls : aqs;
      const int wl = out ? rlw : aqw;
      // 32-bit index arithmetic (nrl, ncp < 2^31, checked on the host)
      int r = static_cast<int>(t);
      const int v = r % n;
      r /= n;
      const int kk = r % nx3;
      r /= nx3;
      const int jj = r % nx2;
      const int m = r/nx2;
      const int64_t qi = v + static_cast<int64_t>(n)*(kk + static_cast<int64_t>(nx3)*jj);
      const int64_t qc = oc[m*ocs0 + i*ocs1] + qi, qb = ocb[m*ocs0 + i*ocs1] + qi;
      const int64_t qw = Wadr(s, wl, m, v, kk + w, jj + w);
      if (out) {
        da[qc] = da[qw];
        if (two) {db[qb] = db[qw];}
      } else {
        da[qw] = da[qc];
        if (two) {db[qw] = db[qb];}
      }
      return;
    }
    t -= nrl + ncp;
    if (xb && t < nrs) {
      const bool b2 = (t >= rsn);
      const int64_t tl = b2 ? (t - rsn) : t;
      int m, v, kb, jb;
      if (!Rec(rsk, aqw, b2 ? rsdl2 : rsdl, b2 ? rsdw2 : rsdw, b2 ? rsnr2 : rsnr, tl,
               m, v, kb, jb)) {return;}
      const int64_t q = Wadr(aqs, aqw, m, v, kb, jb);
      if (q < 0) {return;}
      if (b2) {
        db[q] = st[rsso + t];
      } else {
        da[q] = st[rsso + t];
      }
      return;
    }
    if (t < nrs) {
      int m, v, kb, jb;
      if (!Rec(rsk, aqw, rsdl, rsdw, rsnr, t, m, v, kb, jb)) {return;}
      const int64_t q = Wadr(aqs, aqw, m, v, kb, jb);
      if (q < 0) {return;}
      da[q] = st[rsso + t];
      if (two) {db[q] = st[rsso + rsn + t];}
      return;
    }
    // the tables (no other thread of this kernel reads them)
    for (int m = 0; m < nm; ++m) {
      if (hrl) {
        aoff[m*offs0 + rli*offs1] = oc[m*ocs0 + rli*ocs1];
        if (two) {boff[m*offs0 + rli*offs1] = ocb[m*ocs0 + rli*ocs1];}
      }
      if (haq) {
        aoff[m*offs0 + aqi*offs1] = so[aqs*sos0 + m*sos1];
        if (two) {boff[m*offs0 + aqi*offs1] = so[aqs*sos0 + m*sos1];}
      }
    }
    if (hrl) {
      awi[rli] = 0;
      if (two) {bwi[rli] = 0;}
    }
    if (haq) {
      awi[aqi] = aqw;
      if (two) {bwi[aqi] = aqw;}
    }
  }
};

static_assert(sizeof(VgdRingK) <= 512, "VgdRingK: keep below the by-value launch limit");

void RadiationVet::VetGdRingSweepBegin(const bool two) {
  vgd_ring_act = vgd_ring;
  if (!vgd_ring_act) {return;}
  for (size_t i = 0; i < vgd_rg_state.size(); ++i) {
    if (vgd_rg_state[i] >= 0) {VgdFatal("vet_gd_band_ring: a shell still holds a slot");}
  }
  for (int s = 0; s < vgd_rg_k; ++s) {vgd_rg_sown[s] = -1;}
  vgd_rg_pend.clear();
  vgd_rg_two = two;
  vgd_rg_xb = vgd_xff_now;   // the second array is the other set's
  vgd_rg_stop = 0;
  vgd_rg_nk += 1.0;
  // store of the saved inward band values: the exact lists of this sweep's inward
  // exchanges (made by the batch; grown below if needed)
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1c = indcs.nx1 + 2*indcs.ng;
  int64_t need = 0;
  if (static_cast<int>(vgd_hxe.size()) == 2*n1c) {
    for (int i = indcs.is; i <= indcs.ie; ++i) {
      const VgdHlEntry &X = vgd_hxe[i];
      if (X.gen != -1 && X.exact && X.alpha == vgd_alpha && VetGdHxCut(i, X)) {
        need += static_cast<int64_t>(X.nr) + X.nw;
      }
    }
  }
  need *= (two ? 2 : 1);
  if (static_cast<int64_t>(vgd_rg_store.extent(0)) < need) {
    Kokkos::fence();
    vgd_rg_store = DvceArray1D<Real>();
    Kokkos::realloc(vgd_rg_store, need + need/16 + 1024);
  }
}

void RadiationVet::VetGdRingSweepEnd() {
  if (!vgd_ring_act) {return;}
  VetGdRingFlush();
  for (size_t i = 0; i < vgd_rg_state.size(); ++i) {
    if (vgd_rg_state[i] >= 0) {VgdFatal("vet_gd_band_ring: a shell kept its slot");}
  }
  vgd_ring_act = false;
}

int RadiationVet::VetGdRingFreeSlot(const int avoid) const {
  for (int s = 0; s < vgd_rg_k; ++s) {
    if (vgd_rg_sown[s] < 0 && s != avoid) {return s;}
  }
  return -1;
}

void RadiationVet::VetGdRingFlush() {
  if (!vgd_ring_act || vgd_rg_pend.empty()) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int n = vgd_n, nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2b = nx2 + 2*w, c3b = nx3 + 2*w;
  const bool two = vgd_rg_two;
  const bool xb = two && vgd_rg_xb;
  // one group: at most one save, one release, one acquire
  size_t p = 0;
  while (p < vgd_rg_pend.size()) {
    VgdRingOp sv, rl, aq;
    while (p < vgd_rg_pend.size()) {
      const VgdRingOp &o = vgd_rg_pend[p];
      if ((o.kind == 1 && sv.kind) || (o.kind == 2 && rl.kind) ||
          (o.kind == 3 && aq.kind)) {
        break;
      }
      if (o.kind == 1) {
        sv = o;
        sv.slot = vgd_rg_state[o.i];
        if (sv.slot < 0) {VgdFatal("vet_gd_band_ring: save of a compact shell");}
      } else if (o.kind == 2) {
        rl = o;
        rl.slot = vgd_rg_state[o.i];
        if (rl.slot < 0) {VgdFatal("vet_gd_band_ring: release of a compact shell");}
        if (sv.kind && sv.i != o.i) {VgdFatal("vet_gd_band_ring: save / release differ");}
      } else {
        aq = o;
        if (vgd_rg_state[o.i] >= 0) {VgdFatal("vet_gd_band_ring: shell already lent");}
        aq.slot = VetGdRingFreeSlot(rl.kind ? rl.slot : -1);
        if (aq.slot < 0) {VgdFatal("vet_gd_band_ring: no free slot");}
      }
      ++p;
    }
    // save record of the shell saved here, and the store
    int64_t sv_so = 0, sv_n = 0, sv_n2 = 0;
    int sv_k = 0;
    DvceArray1D<int> sv_dl, sv_dw, sv_dl2, sv_dw2;
    int sv_nr = 0, sv_nr2 = 0;
    if (sv.kind) {
      const int i = sv.i;
      const int wl = vgd_rg_wfull[i];
      if (vgd_rg_sk[i] != 1 && vgd_rg_sk[i] != 3) {vgd_rg_sk[i] = 2;}   // whole band
      sv_k = vgd_rg_sk[i];
      if (xb && sv_k == 3) {   // the other set's own kept list
        VgdGbSwap();
        sv_dl2 = VetGdRawRec(vgd_rg_se[i]);
        VgdGbSwap();
        sv_nr2 = sv_dl2.extent_int(0);
      }
      if (sv_k == 3) {   // vgd_raw: the positions the raw exchange wrote
        sv_dl = VetGdRawRec(vgd_rg_se[i]);
        sv_nr = sv_dl.extent_int(0);
        sv_n = sv_nr;
        sv_k = 1;
      } else if (sv_k == 1) {
        const VgdHlEntry &X = vgd_hxe[vgd_rg_se[i]];
        sv_dl = X.dl;
        sv_dw = X.dw;
        sv_nr = X.nr;
        sv_n = static_cast<int64_t>(X.nr) + X.nw;
      } else {
        sv_n = static_cast<int64_t>(nmb)*n*(nx3 + 2*wl)*(nx2 + 2*wl);
        vgd_rg_ndense += 1.0;
      }
      sv_n2 = (xb && vgd_rg_sk[i] == 3) ? sv_nr2 : sv_n;
      if (xb && sv_k == 1 && vgd_rg_sk[i] != 3) {
        const VgdHlEntry &X2 = vgd_hxe_alt[vgd_rg_se[i]];
        sv_dl2 = X2.dl;
        sv_dw2 = X2.dw;
        sv_nr2 = X2.nr;
        sv_n2 = static_cast<int64_t>(X2.nr) + X2.nw;
      }
      if (static_cast<int>(vgd_rg_sn2.size()) != static_cast<int>(vgd_rg_sn.size())) {
        vgd_rg_sn2.assign(vgd_rg_sn.size(), 0);
      }
      vgd_rg_sn2[i] = sv_n2;
      const int64_t need = vgd_rg_stop + (xb ? (sv_n + sv_n2) : sv_n*(two ? 2 : 1));
      if (static_cast<int64_t>(vgd_rg_store.extent(0)) < need) {
        Kokkos::fence();
        Kokkos::resize(vgd_rg_store, need + need/8 + 1024);
      }
      sv_so = vgd_rg_stop;
      vgd_rg_so[i] = sv_so;
      vgd_rg_sn[i] = sv_n;
      vgd_rg_stop = need;
      // keep the saved values' lists with the record (the restore reads them)
    }
    // restore record of the acquired shell (outward pass)
    int64_t rs_so = 0, rs_n = 0, rs_n2 = 0;
    int rs_k = 0, rs_nr = 0, rs_nr2 = 0;
    DvceArray1D<int> rs_dl, rs_dw, rs_dl2, rs_dw2;
    if (aq.kind && aq.restore) {
      const int i = aq.i;
      rs_k = vgd_rg_sk[i];
      rs_so = vgd_rg_so[i];
      rs_n = vgd_rg_sn[i];
      if (xb && rs_k == 3) {   // the other set's raw record
        VgdGbSwap();
        rs_dl2 = VetGdRawRec(vgd_rg_se[i]);
        VgdGbSwap();
        rs_nr2 = rs_dl2.extent_int(0);
      }
      if (rs_k == 3) {
        rs_dl = VetGdRawRec(vgd_rg_se[i]);
        rs_nr = rs_dl.extent_int(0);
        rs_k = 1;
        if (rs_nr != rs_n) {VgdFatal("vet_gd_band_ring: the raw record changed");}
      } else if (rs_k == 1) {
        const VgdHlEntry &X = vgd_hxe[vgd_rg_se[i]];
        rs_dl = X.dl;
        rs_dw = X.dw;
        rs_nr = X.nr;
        if (static_cast<int64_t>(X.nr) + X.nw != rs_n) {
          VgdFatal("vet_gd_band_ring: the exact list changed between save and restore");
        }
      } else if (rs_k != 2) {
        VgdFatal("vet_gd_band_ring: restore without a save");
      }
      if (xb) {
        rs_n2 = vgd_rg_sn2[i];
        if (vgd_rg_sk[i] == 3 && rs_nr2 != rs_n2) {
          VgdFatal("vet_gd_band_ring: the other set's raw record changed");
        }
        if (rs_k == 1 && vgd_rg_sk[i] != 3) {
          const VgdHlEntry &X2 = vgd_hxe_alt[vgd_rg_se[i]];
          rs_dl2 = X2.dl;
          rs_dw2 = X2.dw;
          rs_nr2 = X2.nr;
          if (static_cast<int64_t>(X2.nr) + X2.nw != rs_n2) {
            VgdFatal("vet_gd_band_ring: the other set's exact list changed between save "
                     "and restore");
          }
        }
      }
    }
    auto a_ = vgd_i;
    auto b_ = xb ? vgd_i_xb : vgd_i;
    auto sb_ = vgd_i.sb;
    auto oc_ = a_.offc, ocb_ = b_.offc;
    auto so_ = vgd_rg_dsoff;
    auto st_ = vgd_rg_store;
    const int64_t nsv = sv.kind ? (xb ? (sv_n + sv_n2) : sv_n) : 0;
    const int64_t nin = static_cast<int64_t>(nmb)*n*nx3*nx2;
    const int64_t nrl = rl.kind ? nin : 0;
    const int64_t ncp = (aq.kind && aq.copyin) ? nin : 0;
    const int64_t nrs = (aq.kind && aq.restore) ? (xb ? (rs_n + rs_n2) : rs_n) : 0;
    const int64_t ntot = nsv + nrl + ncp + nrs + 1;
    if (nin >= static_cast<int64_t>(INT32_MAX)) {
      VgdFatal("vet_gd_band_ring: shell interior >= 2^31 values (32-bit copy index)");
    }
    // the copy kernel: a small POD functor of raw pointers (a functor of Views exceeds
    // the 512 B by-value launch limit and would go through constant memory)
    VgdRingK k;
    k.da = a_.d.data();
    k.db = b_.d.data();
    k.st = st_.data();
    k.aoff = a_.off.data();
    k.boff = b_.off.data();
    k.awi = a_.wi.data();
    k.bwi = b_.wi.data();
    k.oc = oc_.data();
    k.ocb = ocb_.data();
    k.so = so_.data();
    k.sb = sb_.data();
    k.offs0 = a_.off.stride(0);
    k.offs1 = a_.off.stride(1);
    k.ocs0 = oc_.stride(0);
    k.ocs1 = oc_.stride(1);
    k.sos0 = so_.stride(0);
    k.sos1 = so_.stride(1);
    k.sbs0 = static_cast<int>(sb_.stride(0));
    k.sbs1 = static_cast<int>(sb_.stride(1));
    if (b_.off.stride(0) != a_.off.stride(0) || b_.off.stride(1) != a_.off.stride(1)) {
      VgdFatal("vet_gd_band_ring: main / twin offset tables differ in layout");
    }
    k.svdl = sv_dl.data();
    k.svdw = sv_dw.data();
    k.rsdl = rs_dl.data();
    k.rsdw = rs_dw.data();
    k.nsv = nsv;
    k.nrl = nrl;
    k.ncp = ncp;
    k.nrs = nrs;
    k.svso = sv_so;
    k.rsso = rs_so;
    k.svn = sv_n;
    k.rsn = rs_n;
    k.n = n;
    k.nx2 = nx2;
    k.nx3 = nx3;
    k.w = w;
    k.c2b = c2b;
    k.c3b = c3b;
    k.nm = nmb;
    k.svk = sv_k;
    k.rsk = rs_k;
    k.svnr = sv_nr;
    k.rsnr = rs_nr;
    k.svs = sv.slot;
    k.svw = sv.kind ? vgd_rg_wfull[sv.i] : 0;
    k.rli = rl.i;
    k.rls = rl.slot;
    k.rlw = rl.kind ? vgd_rg_wfull[rl.i] : 0;
    k.aqi = aq.i;
    k.aqs = aq.slot;
    k.aqw = aq.kind ? vgd_rg_wfull[aq.i] : 0;
    k.two = two;
    k.hrl = rl.kind;
    k.haq = aq.kind;
    k.xb = xb;
    k.svdl2 = sv_dl2.data();
    k.svdw2 = sv_dw2.data();
    k.rsdl2 = rs_dl2.data();
    k.rsdw2 = rs_dw2.data();
    k.svn2 = sv_n2;
    k.rsn2 = rs_n2;
    k.svnr2 = sv_nr2;
    k.rsnr2 = rs_nr2;
    Kokkos::parallel_for("vet_vgd_ring", Kokkos::RangePolicy<>(vgd_cur, 0, ntot), k);
    if (rl.kind) {
      vgd_rg_sown[rl.slot] = -1;
      vgd_rg_state[rl.i] = -1;
    }
    if (aq.kind) {
      vgd_rg_sown[aq.slot] = aq.i;
      vgd_rg_state[aq.i] = aq.slot;
    }
  }
  vgd_rg_pend.clear();
}


//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::DtRetryQuiesce
//! \brief before the dt-retry checkpoint is taken (restoring = false) or put back
//! (true): join a pending symmetric-reference solve and, before a restore, the table
//! prefetch. Changes no number (vet_gd_async is refused with dt_retry).

extern int64_t g_vgd_csgen[16];       // rad_vet_vetgd_cs.cpp; VET_VGD_CSCL_MAX
extern int64_t g_vgd_csgen_save[16];

void RadiationVet::DtRetryQuiesce(bool restoring) {
  if (!restoring) {
    // compact-list build generations at the checkpoint
    for (int c = 0; c < vgd_csncl; ++c) {g_vgd_csgen_save[c] = g_vgd_csgen[c];}
    CleanRingPreShift();
    CleanDtrD(true);
    return;
  }
  if (vgd_tab_state == 1 && vgd_tabf.valid()) {vgd_tabf.wait();}
  if (vgd_on) {vgd_tex.fence();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::DtRetryAfterRestore
//! \brief after the checkpoint is put back: reset the table prefetch state (not
//! checkpointed) to "nothing prefetched"; the next rotation computes its tables
//! synchronously (VgdTabCompute: same numbers).

void RadiationVet::DtRetryAfterRestore() {
  CleanDtrD(false);
  if (vgd_tabf.valid()) {vgd_tabf.get();}
  if (vgd_on) {
    vgd_tex.fence();
    vgd_tab_state = 0;
    vgd_tabh.alpha = -9.0;
    // the sweep's halo caches (vgd_hle, vgd_hxe, mask slots vgd_htag, overlap groups)
    // have checkpointed keys but arrays remade in place: a failed attempt may leave a
    // restored key naming another cut's lists. Drop every record; the next sweep
    // remakes the ones it needs (same masks: same numbers).
    for (auto &x : vgd_hle) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    for (auto &x : vgd_hxe) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    for (auto &x : vgd_hxe_alt) {x.gen = -1; x.alpha = -9.0; x.nofit = false;}
    vgd_shl_alpha = vgd_shl_alpha_xb = -3.0;   // the shell lists (vgd_shl[_xb]) too
    for (int s = 0; s < 2; ++s) {vgd_htag[s][0] = -1;}
    // the cubed-sphere compact lists (vgd_cscl) are the same kind of record: drop an
    // entry the failed attempt rebuilt (LRU eviction), else the restored key names the
    // evicted set's arrays; the next sweep remakes it (same numbers)
    for (int c = 0; c < vgd_csncl; ++c) {
      if (g_vgd_csgen[c] != g_vgd_csgen_save[c]) {
        vgd_cscl[c].alpha = -9.0;
        vgd_cscl[c].use = 0;
      }
    }
    // the two read masks may hold other sets than the restored handles say: every
    // cached set remakes its mask at its next use (VetGdCsCompact)
    for (auto &e : vgd_cscl) {e.msk = decltype(e.msk)();}
    // a failed attempt on the last step of a rotation block left the tables on the
    // next set while the restored vgd_alpha names the current one: remake the tables
    // of the restored angle (VgdTabCompute: same numbers)
    if (vgd_alpha >= 0.0 && g_vgd_tab_alpha != vgd_alpha) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> dt_retry: vet_gd tables back from angle "
                  << g_vgd_tab_alpha << " to the restored " << vgd_alpha << std::endl;
      }
      VetGdTables(vgd_alpha);
    }
  }
}

// band exchanges called by VetGdSweep (its own TU)
template void RadiationVet::VetGdHalo<RadiationVet::VgdRag>(VgdRag &a, const int nv,
    const int i0, const int i1, const int ws, const bool mapd);
template bool RadiationVet::VetGdHlBegin<RadiationVet::VgdRag>(VgdRag &a, const int nv,
    const int i0, const int ws);
#endif  // !VET_VGD_SWEEP_TU
} // namespace radvet

#if !defined(VET_VGD_SWEEP_TU)
#include "rad_vet_vetgd_sph2.inc"
#endif
