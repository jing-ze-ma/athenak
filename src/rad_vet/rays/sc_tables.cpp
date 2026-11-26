//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_tables.cpp
//! \brief curved-mesh rays: wall maps and branch tables of a (rotated) direction set.
//!
//!  called from: rad_vet_closure_ring.cpp, rad_vet_restart.cpp, rays/sc_build.cpp,
//!    rays/sc_setup.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_cs_setup.cpp, rays/sc_directions.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <chrono>  // NOLINT(build/c++11)
#include <cmath>
#include <future>
#include <iostream>
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <utility>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_halo.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn static void ScMirPick
//! \brief vet_gd_theta_wall = mirror: set directions standing in for direction m at a
//! cell with radial unit vector (rx, ry, rz): an exact match, else the nearest of m's
//! branch (sign of m . r_hat); with wint the 3 nearest, barycentric weights when all
//! >= 0, else inverse-angle (vet_gd_wall_interp's rule)

static void ScMirPick(const std::vector<double> &d, const int n, const double mx,
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
//! \fn void RadiationVet::ScTabCompute
//! \brief host part of ScTables for angle alpha (rotated set, wall maps, branch
//! tables) into the pinned host buffers vgd_tabh; touches nothing else (safe on the
//! prefetch thread). nthr: host threads over (m, k) columns (result independent).

void RadiationVet::ScTabCompute(const Real alpha, const int nthr) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c2 = indcs.nx2 + 2*vgd_w;
  const int c3 = indcs.nx3 + 2*vgd_w;
  const int n = vgd_n;
  ScTabHost &H = vgd_tabh;
  std::vector<double> d = ScRotate(alpha);
  auto &d_h = H.dir;
  for (int q = 0; q < n; ++q) {
    for (int c = 0; c < 4; ++c) {d_h(q,c) = d[4*q+c];}
  }
  if (vgd_csg) {   // no walls; the branch of every band cell
    ScCsTab(d);
    H.alpha = alpha;
    return;
  }
  if (vgd_slab == 1) {
    // Cartesian slab: a periodic image is a translation (no wall, identity map) and the
    // branch n . r_hat is n_1 everywhere
    for (int m = 0; m < nmb; ++m) {
      for (int k = 0; k < c3; ++k) {
        for (int j = 0; j < c2; ++j) {
          H.w(m,k,j) = 0;
          for (int q = 0; q < n; ++q) {
            H.mp(m,k,j,q) = q;
            H.mr(m,k,j,q) = d[4*q];
          }
        }
      }
    }
    H.alpha = alpha;
    H.nwall = H.nexact = H.nmap = H.nphx = H.nthx = 0;
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
              // vet_gd_wall_sym, phi wall: ScMirPick, linear in azimuth on n_q's ring
              // (positive, mirror-paired)
              ScMirPick(d, n, mx, my, mz, sw*cpw, sw*spw, cw, true, r3, w3, true);
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
            const int jm = ScMirRow(g, sd == 1, wb, nx2);
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
              ScMirPick(d, n, mx, my, mz, rx, ry, rz, wint, r3, w3);
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
//! \fn void RadiationVet::ScTabAlloc
//! \brief pinned host buffers of the tables (once, zero-filled) and, with dev, the
//! second set of device tables the prefetch uploads into

void RadiationVet::ScTabAlloc(const bool dev) {
  ScTabHost &H = vgd_tabh;
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
//! \fn void RadiationVet::ScTabPrefetch
//! \brief compute the next rotation block's tables on a host thread (ScTabCompute,
//! same numbers); ScTabPoll uploads them into the second device set, ScTables
//! swaps them in at the next rotation

void RadiationVet::ScTabPrefetch(const int cyc) {
  if (!ScRotOn() || vgd_tab_state != 0 || vgd_xf) {return;}
  const Real an = RotAngleBlk(RotBlock(cyc) + 1);
  if (an == vgd_alpha) {return;}
  ScTabAlloc(true);
  vgd_tab_state = 1;
  vgd_tabf = std::async(std::launch::async, [this, an]() {ScTabCompute(an, 8);});
}

void RadiationVet::ScTabPoll() {
  if (vgd_tab_state != 1) {return;}
  if (vgd_tabf.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {return;}
  vgd_tabf.get();
  // the second device set is not read by anything in flight (swapped out at the last
  // rotation, fenced since): upload on the table instance
  ScTabHost &H = vgd_tabh;
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

// angle of the tables now in vgd_dir/wall/map/mr/m3/w3 (not a class member: the
// dt-retry checkpoint restores vgd_alpha but not these Views)
Real g_vgd_tab_alpha = -9.0;

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScTables
//! \brief the set rotated by alpha -> vgd_dir, wall maps and branch tables (from the
//! prefetched device set when it holds alpha, else computed now and copied)

void RadiationVet::ScTables(const Real alpha) {
  ScTabPoll();
  ScTabHost &H = vgd_tabh;
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
    ScTabAlloc(false);
    if (H.alpha != alpha) {ScTabCompute(alpha, 16);}
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
    ScFatal("vet_gd_band_global: a wall direction map is not the identity (the global "
            "band copies do not re-index directions)");
  }
  if (H.nphx > 0 && !vgd_wsym) {
    // any set: an interpolated / nearest phi-wall map breaks the set's phi mirror
    // symmetry (a net F_phi), unless vet_gd_wall_sym pairs the maps
    ScFatal("phi-periodic wedge width W must be a multiple of the ray set's azimuthal "
            "period for exact wall maps (" + std::to_string(H.nphx) + " of " +
            std::to_string(nmap) + " phi-wall direction maps inexact), or name "
            "rad_vet/vet_gd_wall_sym = true (with vet_gd_wall_interp = true) for "
            "mirror-paired phi walls");
  }
  if (H.nthx > 0 && !vgd_wst &&
      pmy_pack->pmesh->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic) {
    // a theta-periodic mesh: its clean theta seam face is one-frame (a net F_theta); a
    // reflecting theta wall has no seam face
    ScFatal("theta-periodic wedge: the theta wall maps are inexact (" +
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
              << (vgd_csg ? "cubed sphere"
                          : ((vgd_slab == 1) ? "Cartesian slab" : "sp wedge")) << ", "
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
//! \fn void RadiationVet::ScXfTables
//! \brief vet_gd_rotate_xfade: the tables of alpha into vgd_dir/...; both cross-fade
//! sets stay on the device (current and vgd_tabd) and are swapped; another angle swaps
//! the current set out and builds the new one (ScTables)

void RadiationVet::ScXfTables(const Real alpha) {
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
    ScGbSwap();
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
  ScTabAlloc(true);
  swp();
  vgd_xf_alt = g_vgd_tab_alpha;
  ScTables(alpha);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScXffSwap
//! \brief swap the cross-fade's whole set context (tables, angles, exact lists,
//! intensities, shell lists) without a fence (handles only); between the two halo
//! exchanges of a fused shell when the merged message does not apply

void RadiationVet::ScXffSwap() {
  std::swap(vgd_dir, vgd_tabd.dir);
  std::swap(vgd_wall, vgd_tabd.wall);
  std::swap(vgd_map, vgd_tabd.map);
  std::swap(vgd_mr, vgd_tabd.mr);
  if (vgd_wint) {
    std::swap(vgd_m3, vgd_tabd.m3);
    std::swap(vgd_w3, vgd_tabd.w3);
  }
  if (!vgd_csg) {std::swap(vgd_hxe, vgd_hxe_alt);}
  ScGbSwap();
  std::swap(vgd_i, vgd_i_xb);
  for (int p = 0; p < 2; ++p) {
    std::swap(vgd_shl[p], vgd_shl_xb[p]);
    std::swap(vgd_shn[p], vgd_shn_xb[p]);
  }
  std::swap(vgd_shl_alpha, vgd_shl_alpha_xb);
  std::swap(vgd_alpha, vgd_xf_alt);
  g_vgd_tab_alpha = vgd_alpha;
}

}  // namespace radvet
