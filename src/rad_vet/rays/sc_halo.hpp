#ifndef RAD_VET_RAYS_SC_HALO_HPP_
#define RAD_VET_RAYS_SC_HALO_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_halo.hpp
//! \brief curved-mesh rays: helpers of the sweep's halos (GPU flags, scans, candidate
//! cells).

#include <math.h>

#include <cstdint>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/grid_stretch.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#if defined(KOKKOS_ENABLE_CUDA)
#include <cub/cub.cuh>
#endif

namespace radvet {

// events behind the per-shell halo pack (rays/sc_sph_halo.cpp)
void ScPackEvtRecord(const DevExeSpace &ex);
void ScPackEvtWait(const DevExeSpace &ex);

static inline
void ScFatal(const std::string &msg) {ScFatalAt(__FILE__, "<rad_vet>/vet_gd: ", msg);}

KOKKOS_INLINE_FUNCTION
static Real ScLerp(const Real a, const Real b, const Real u) {return a + u*(b - a);}

// sp stencil accessor of the shared upwind core (rays/sc_ray.hpp): the four
// (theta, phi) corners q = 0..3: (k0,j0), (k0,j1), (k1,j0), (k1,j1); nested bilinear
struct ScSt4 {
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
    return ScLerp(ScLerp(f(k0,j0), f(k0,j1), uj), ScLerp(f(k1,j0), f(k1,j1), uj), uk);
  }
};

// vet_gd_halo_ipc: device flags written over NVLink by the peers: [0, NS) arrivals by
// partner slot, [NS, 2 NS) acks; a spin gives up after ~20 s and raises the error flag
constexpr int VGD_IPC_NS = 32;

constexpr long long VGD_IPC_SPIN = 40000000000LL;

// HIP: wall_clock64() runs at 100 MHz on gfx942: 2e9 ticks = 20 s
constexpr long long VGD_IPC_SPIN_HIP = 2000000000LL;

struct ScIpcIdx {int v[VGD_IPC_NS + 1];};

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
static void ScIpcWaitGE(const long long *f, const long long target, int *err) {
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
static void ScIpcSet(long long *f, const long long v) {
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
static KOKKOS_INLINE_FUNCTION
void ScIpcSetN(const unsigned long long *psig, const int *pslot, const int np,
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

static KOKKOS_INLINE_FUNCTION
void ScIpcWaitN(const long long *f, const int np, const long long target, int *err) {
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
static KOKKOS_INLINE_FUNCTION
bool ScCornerEdge(const CS &csx, const int m, const int k0, const int k1, const int j0,
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

// every stencil box f(iu, k_lo, k_hi, j_lo, j_hi) (inclusive band cells of shell iu)
// that the shell kernel's solve of (m, k, j, d) at shell i may read. Robust to rounding
// of another compilation: both branches near the threshold (rel 1e-9), every floor
// within 1e-8. A superset of the reads: used by the overlap split and the halo lists.
// sth / pol / tn / c2w: stretched-theta node table and full-sphere chart (sph2.inc)
template <class TD, class TB, class TX, class TF, class TT, class FN>
static KOKKOS_INLINE_FUNCTION
void ScCand(const int m, const int k, const int j, const int d, const int i, const int l,
            const bool inw, const int n1, const int lcut, const int js, const int ks,
            const int jlo, const int jhi, const int klo, const int khi, const TD &dir_,
            const TB &mbsize, const ScLG &lg_, const TX &cx1v, const TF &cx1f,
            const bool sth,
            const bool pol, const TT &tn, const int c2w, const bool imir, const int xr,
            const int slb, const FN &f) {
  if (inw && l == n1 - 1) {return;}          // the top shell reads no intensity
  // Cartesian slab (slb 1: 2-D, 2: 3-D): x1 plays r, the foot is a constant offset per
  // direction (as the shell kernel); the first shell's outward rays read the inner face
  if (slb > 0) {
    if (!inw && l == lcut) {return;}
    const int iu = inw ? (i + 1) : (i - 1);
    const Real ds = fabs(cx1v(m,iu) - cx1v(m,i))/fabs(dir_(d,0));
    const Real fj = ScLG::Fs(j, lg_.Sj(m, js), -ds*dir_(d,1)/lg_.Dx2(mbsize, m));
    const Real dk = (slb > 1) ? -ds*dir_(d,2)/lg_.Dx3(mbsize, m) : 0.0;
    const Real fk = ScLG::Fs(k, lg_.Sk(m, ks), dk);
    const Real e = 1.0e-8;
    int ja = static_cast<int>(floor(fj - e)), jb = static_cast<int>(floor(fj + e)) + 1;
    int ka = (slb > 1) ? static_cast<int>(floor(fk - e)) : k;   // (2-D: k1 = k0)
    int kb = (slb > 1) ? (static_cast<int>(floor(fk + e)) + 1) : k;
    ja = (ja - xr < jlo) ? jlo : (ja - xr);
    jb = (jb + xr > jhi + 1) ? (jhi + 1) : (jb + xr);
    ka = (ka - xr < klo) ? klo : (ka - xr);
    kb = (kb + xr > khi + 1) ? (khi + 1) : (kb + xr);
    f(iu, ka, kb, ja, jb, d);
    return;
  }
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
    ScPoleChart(pol, th, thu, dph);
    const Real fj = sth ? (ScThFj(tn, m, j, thu, c2w, lg_.Sj(m, js)) - lg_.Sj(m, js))
                        : ScLG::Fs(j, lg_.Sj(m, js), (thu - th)/lg_.Dx2(mbsize, m));
    const Real fk = ScLG::Fs(k, lg_.Sk(m, ks), dph/lg_.Dx3(mbsize, m));
    const Real e = 1.0e-8;
    int ja = static_cast<int>(floor(fj - e)), jb = static_cast<int>(floor(fj + e));
    int ka = static_cast<int>(floor(fk - e)), kb = static_cast<int>(floor(fk + e));
    ja = (ja < jlo) ? jlo : ((ja > jhi) ? jhi : ja);
    jb = (jb < jlo) ? jlo : ((jb > jhi) ? jhi : jb);
    ka = (ka < klo) ? klo : ((ka > khi) ? khi : ka);
    kb = (kb < klo) ? klo : ((kb > khi) ? khi : kb);
    // xr: reads beyond the 2x2 stencil (vgd_xr), inside the band
    const int kx0 = (ka - xr < klo) ? klo : (ka - xr);
    const int kx1 = (kb + 1 + xr > khi + 1) ? (khi + 1) : (kb + 1 + xr);
    const int jx0 = (ja - xr < jlo) ? jlo : (ja - xr);
    const int jx1 = (jb + 1 + xr > jhi + 1) ? (jhi + 1) : (jb + 1 + xr);
    f(iu, kx0, kx1, jx0, jx1, dq);
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

static inline
void ScScan2(const DevExeSpace &ex, const DvceArray1D<int> &a, const size_t na,
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
// the instance's stream (integer sums: the positions of ScScan2). Non-CUDA: ScScan2.
static inline
void ScScanOne(const DevExeSpace &ex, const DvceArray1D<int> &a, const size_t n,
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
  ScScan2(ex, a, n, a, 0, "vet_vgd_hxb_scan", sums);
  (void)tmp;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn KOKKOS_INLINE_FUNCTION bool ScHaloKeep
//! \brief vet_gd_halo_compact: is direction u at source cell (global gj, gk) at band
//! depth dep needed by the reads of shell s in this pass? Global indices only, so
//! sender and receiver agree bitwise. mode 1: the pass's branch (|n.r| < 1e-12 kept);
//! mode 2: also the ray's lateral reach >= dep - 3 cells. wall: branch only.

// the keep test at a source-cell position (theta, phi): shared by the uniform slot halo
// (ScHaloKeep) and the global band (ScGbHaloC, node table)
static KOKKOS_INLINE_FUNCTION
bool ScHaloKeepTP(const Real nx, const Real ny, const Real nz, const Real th,
                  const Real ph, const bool inw, const int mode, const bool wall,
                  const int dep, const Real rs, const Real rm, const Real rp,
                  const bool hasm, const bool hasp, const Real gam, const Real dmin,
                  const bool slb = false) {
  if (slb) {   // Cartesian slab: n . r_hat = n_1 everywhere; the branch alone
    return (fabs(nx) < 1.0e-12) || (inw == (nx < 0.0));
  }
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

// mesh theta face n (0 .. nx2), as Coordinates maps it (LeftEdgeX + ApplyThetaStretch)
static inline double ScThetaFace(const Mesh *pm, const int n) {
  const int nx = pm->mesh_indcs.nx2;
  const double x0 = pm->mesh_size.x2min, x1 = pm->mesh_size.x2max;
  const double rx = static_cast<double>(n)/static_cast<double>(nx);
  Real t = static_cast<Real>(rx*x1 - rx*x0 + x0);
  if (pm->use_grid_stretch_theta) {
    ApplyThetaStretch(pm->fStretchTheta, pm->fStretchThetaPlat, t);
  }
  return static_cast<double>(t);
}

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_HALO_HPP_
