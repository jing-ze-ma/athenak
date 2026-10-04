//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd.cpp
//! \brief <rad_m1>/vet_gd = true (m1-vet-gd; design rt_design_1003/SC_PROPER_SP.md,
//! candidate b): the Eddington tensor on the spherical-polar wedge from a
//! short-characteristics (SC) formal solution along GLOBALLY FIXED directions.
//!
//! WHAT IS SOLVED.  At every vet_col build (the same state, source S and extinction chi
//! as vet_col and vet_col_lat: VetLatSweep(0) fills ln chi, ln S and the first shell
//! vlat_icut), the grey transfer equation  n . grad I = chi (S - I)  is swept shell by
//! shell for the N = 12 vet_gd_nside^2 HEALPix ring centres n_d, fixed in the CARTESIAN
//! frame (straight rays: the curvature of the sp grid enters only through the positions,
//! there is no angular redistribution and no angular numerical diffusion).  The weights
//! are equal up to one global factor in n_z^2 that makes the set integrate 1, n and nn
//! exactly, so f = 1/3 holds for an isotropic or linear intensity in every cell frame.
//! Two passes, one launch per shell, one thread per (block, cell, direction):
//!   inward, top-down, for the (cell, direction) pairs with n . r_hat < 0: the upwind
//!     point is where the backward ray meets the sphere of shell i+1 (the top shell from
//!     the vacuum top face);
//!   outward, bottom-up, for n . r_hat >= 0: the sphere of shell i-1, or, for a ray that
//!     misses it, the same shell beyond the tangent point (an inward value of the first
//!     pass), or the inner face of the first shell (the diffusion intensity, as
//!     vet_col_lat).
//! For straight rays this order has no cycles (r is monotone on either side of the one
//! tangent point).  At the upwind point I (the SAME direction index), ln chi and ln S are
//! bilinear in (theta, phi) index space; the source is linear in tau along the segment
//! (vet_col's positive first-order weights).  This is rt_design_1003/sc_proper/gd_sc.py,
//! against which the C++ is gated.
//!
//! WALLS.  At a periodic theta or phi wall of the wedge the ghost cells hold the copied
//! intensities of the other side; they are re-indexed after every exchange by the LOCAL
//! FRAME convention of the run's own periodic boundaries (the (r, theta, phi)
//! components of a direction are kept, i.e. n -> M n with M the rotation between the
//! two frames, and the nearest set direction is read).  For a phi period of pi/2 the
//! HEALPix set is invariant and the map is an exact permutation.
//!
//! SEVERAL MESHBLOCKS.  An upwind point inside the block reads the current sweep; one in
//! the lateral ghost band reads the neighbour's intensity of the LAST sweep (exchanged
//! after every sweep): a lagged inflow across block boundaries (block-Jacobi; the first
//! build iterates vet_col_lat_init_iter times).  Reads beyond the ghost band are clamped
//! and counted (very oblique rays through laterally thin cells).
//!
//! HAND-OVER: the vet_col_lat interface (tau_ten slots M1_TT_LAT0.., the mesh basis):
//!   LAT0 = D_rr - f_K(vet_col), folded into slot 0 by VetLatBuild (so slot 0 = D_rr of
//!          the sweep above the first shell, vet_col's f_K below it),
//!   LAT1, LAT2 = D_r,theta and D_r,phi (the lagged Picard term M1SphLat when
//!          vet_col_lat_offdiag),
//!   LAT3, LAT4, LAT5 = D_tt - (1 - D_rr)/2, D_tp, D_pp - (1 - D_rr)/2 (the trace-free
//!          tangential anisotropy about the operator's isotropic tangential pressure:
//!          the lagged term M1SphTan of the transverse face equations,
//!          vet_gd_tangential).
//! All vet_col_lat keys apply (taucut, every, init_iter, offdiag, dump).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {
void VgdFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/vet_gd: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real VgdLerp(const Real a, const Real b, const Real u) {return a + u*(b - a);}

// the positive first-order SC weights (VcolW of rad_m1_vetcol.cpp, term for term)
KOKKOS_INLINE_FUNCTION
void VgdW(const Real dtau, Real &ex, Real &w0, Real &wu) {
  ex = exp(-dtau);
  if (dtau < 1.0e-3) {
    w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
    wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
  } else {
    const Real g = (1.0 - ex)/dtau;
    w0 = 1.0 - g;
    wu = g - ex;
  }
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

// level-symmetric LQ_N (Carlson; mu_1 of Lewis & Miller), positive weights per
// permutation class from the even moments 1, x^4, x^6, ... (sum w nn = I/3 exactly by
// the octahedral symmetry); the Athena++ / Jiang 2021 (Bruls 1999) family
// product set: Gauss-Legendre in n_z (nmu nodes on [-1, 1]) x nphi uniform azimuths per
// ring, alternate rings shifted by half an azimuth step (stagger); C_nphi-invariant about
// z, positive weights, sum w = 1, sum w n = 0, sum w nn = I/3 exactly (no fix)
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
      for (int q = 0; q < static_cast<int>(cls.size()); ++q) {if (cls[q] == t) {id = q;}}
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

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdInit
//! \brief direction table, intensity array, exchange object, wall maps (from VetLatInit)

void RadiationM1::VetGdInit() {
  Mesh *pm = pmy_pack->pmesh;
  if (vgd_nside < 1 || vgd_nside > 6) {VgdFatal("vet_gd_nside must be 1..6");}
  if (!pm->three_d) {VgdFatal("needs a 3-D wedge");}
  if (ibc_x1max != M1_IBC_MARSHAK) {
    // the sweep has a VACUUM top: a reflecting (or any non-Marshak) outer x1 of the M1
    // solve would be inconsistent with it near the top (GD_RESULTS.md, open item)
    VgdFatal("reflecting outer x1 not supported; use implicit_bc_x1max = marshak "
             "(+ vet_col_surface_q)");
  }
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  // the direction set (pluggable: any table of unit vectors + weights summing to 1 with
  // sum w nn = I/3; the sweep never looks inside it)
  if (vgd_gl_nmu > 0) {
    // the product set: its azimuth step must divide the wedge's phi period (exact walls)
    const double per = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const double m = per*vgd_gl_nphi/(2.0*M_PI);
    if (vgd_gl_nphi < 4 || std::fabs(m - std::round(m)) > 1.0e-9) {
      VgdFatal("vet_gd_gl_nphi must be >= 4 with 2 pi/nphi dividing the phi period of "
               "the wedge (exact wall maps)");
    }
    VgdProduct(vgd_gl_nmu, vgd_gl_nphi, vgd_gl_stag, vgd_base);
  } else if (vgd_ls > 0) {
    if (!VgdLevelSym(vgd_ls, vgd_base)) {VgdFatal("vet_gd_ls must be 6, 8, 10 or 12");}
  } else {
    VgdDirections(vgd_nside, vgd_base);
  }
  vgd_n = static_cast<int>(vgd_base.size()/4);
  if (vgd_n > M1_VGD_NMAX) {VgdFatal("too many directions");}
  const int n = vgd_n;
  Kokkos::realloc(vgd_dir, n, 4);
  // the lateral band: the deepest reach of an upwind point in cells, over every shell:
  // a NORM segment reaches sqrt(r_{l+1}^2 - r_l^2), a TURN segment 2 sqrt(r_l^2 -
  // r_{l-1}^2) sideways at most; cell sizes r dtheta and r sin(theta) dphi
  {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto &mbs = pmy_pack->pmb->mb_size;
    const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
    const double smin = std::min(std::sin(t0), std::sin(t1));
    const double dmin = std::min(mbs.h_view(0).dx2, smin*mbs.h_view(0).dx3);
    // the data of shell s is read (a) by shell s-1 in the inward pass, chord <=
    // sqrt(r_s^2 - r_{s-1}^2), (b) by shell s+1 in the outward pass, chord <=
    // sqrt(r_{s+1}^2 - r_s^2), (c) by shell s itself (TURN / inner face), chord <= 2
    // sqrt(r_s^2 - r_lo^2), r_lo = r_{s-1} (the inner face of the first shell); the
    // angular offset <= chord / r_lo; + 2 cells for the bilinear stencil and rounding
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    const int nmax = std::min(indcs.nx2, indcs.nx3);
    vgd_wsh.assign(indcs.nx1 + 2*indcs.ng, 0);
    int wneed = 0;
    for (int l = indcs.is; l <= indcs.ie; ++l) {
      const double r = x1v_h(0,l);
      const double rlo = (l > indcs.is) ? x1v_h(0,l-1) : x1f_h(0,l);
      double a = 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
      if (l < indcs.ie) {
        a = std::max(a, std::sqrt(x1v_h(0,l+1)*x1v_h(0,l+1) - r*r));
      }
      const int wl = static_cast<int>(std::ceil(a/(rlo*dmin))) + 2;
      wneed = std::max(wneed, wl);
      vgd_wsh[l] = std::min(wl, nmax);
    }
    vgd_w = std::min(wneed, nmax);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> vet_gd: exact per-shell lateral halo, band " << vgd_w
                << " cells (needed " << wneed << "), per-shell depth mean "
                << std::accumulate(vgd_wsh.begin(), vgd_wsh.end(), 0.0)/indcs.nx1
                << ((vgd_w < wneed) ? ": MeshBlocks too small, the deepest near-tangent "
                                      "reads are clamped (counted)" : "") << std::endl;
    }
  }
  const int c2w = indcs.nx2 + 2*vgd_w, c3w = indcs.nx3 + 2*vgd_w;
  Kokkos::realloc(vgd_i, nmb, n, c3w, c2w, c1);
  Kokkos::deep_copy(vgd_i, 0.0);
  Kokkos::realloc(vgd_cs, nmb, 2, c3w, c2w, c1);
  Kokkos::realloc(vgd_wall, nmb, c3w, c2w);
  Kokkos::realloc(vgd_map, nmb, c3w, c2w, n);
  Kokkos::realloc(vgd_mr, nmb, c3w, c2w, n);
  VetGdHaloInit();
  if (vgd_rbe > 0) {
    Kokkos::realloc(vgd_fk0, nmb, indcs.nx3 + 2*indcs.ng, indcs.nx2 + 2*indcs.ng, c1);
  }
  vgd_time_halo = (std::getenv("VGD_TIME_HALO") != nullptr);
  vgd_alpha = -1.0;
  VetGdTables(VetGdAngle(pm->ncycle));
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::VetGdAngle
//! \brief vet_gd_rotate_every = N > 0: the z-rotation of the set for cycle c, constant
//! over each block of N cycles: alpha = (pi/2) frac((c/N) g), g = (sqrt 5 - 1)/2 (a
//! deterministic low-discrepancy sequence of the cycle number: reruns and restarts
//! bitwise; pi/2 is the C4 period of the set)

Real RadiationM1::VetGdAngle(const int cyc) const {
  if (vgd_rot <= 0) {return 0.0;}
  const double b = static_cast<double>(cyc/vgd_rot);
  const double g = 0.5*(std::sqrt(5.0) - 1.0);
  return 0.5*M_PI*(b*g - std::floor(b*g));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdTables
//! \brief the set rotated by alpha about z -> vgd_dir, the wall maps and branch tables

void RadiationM1::VetGdTables(const Real alpha) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c2 = indcs.nx2 + 2*vgd_w;
  const int c3 = indcs.nx3 + 2*vgd_w;
  const int n = vgd_n;
  std::vector<double> d(vgd_base);
  const double ca = std::cos(alpha), sa = std::sin(alpha);
  for (int q = 0; q < n; ++q) {
    const double x = vgd_base[4*q], y = vgd_base[4*q+1];
    d[4*q] = ca*x - sa*y;
    d[4*q+1] = sa*x + ca*y;
  }
  auto d_h = Kokkos::create_mirror_view(vgd_dir);
  for (int q = 0; q < n; ++q) {
    for (int c = 0; c < 4; ++c) {d_h(q,c) = d[4*q+c];}
  }
  Kokkos::deep_copy(vgd_dir, d_h);
  // wall maps: lateral ghost cells outside the mesh in theta or phi hold the other
  // side's intensities; new(d) = old(map(d)), map = nearest set direction to M n_d
  auto w_h = Kokkos::create_mirror_view(vgd_wall);
  auto mp_h = Kokkos::create_mirror_view(vgd_map);
  // positions of the wall maps and branch tables: the face MIDPOINTS (x2v is the
  // volume centroid on sp; the periodic image of a midpoint is exact)
  auto &mbs = pmy_pack->pmb->mb_size;
  const int jsg = vgd_w, ksg = vgd_w;   // the band-index origin of the active cells
  auto x2_h = [&](const int m, const int j) {
    return mbs.h_view(m).x2min + (j - jsg + 0.5)*mbs.h_view(m).dx2;
  };
  auto x3_h = [&](const int m, const int k) {
    return mbs.h_view(m).x3min + (k - ksg + 0.5)*mbs.h_view(m).dx3;
  };
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  const double p0 = pm->mesh_size.x3min, p1 = pm->mesh_size.x3max;
  int nwall = 0, nexact = 0, nmap = 0;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        const double th = x2_h(m,j), ph = x3_h(m,k);
        double thw = th, phw = ph;
        if (th < t0) {thw = th + (t1 - t0);}
        if (th > t1) {thw = th - (t1 - t0);}
        if (ph < p0) {phw = ph + (p1 - p0);}
        if (ph > p1) {phw = ph - (p1 - p0);}
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
          if (bd > 1.0 - 1.0e-12) {++nexact;}
        }
      }
    }
  }
  // the branch of every stored value: n . r_hat of the direction actually stored at a
  // lateral cell (the mapped one, at the source position, for a wall ghost)
  auto mr_h = Kokkos::create_mirror_view(vgd_mr);
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        double th = x2_h(m,j), ph = x3_h(m,k);
        if (th < t0) {th += (t1 - t0);}
        if (th > t1) {th -= (t1 - t0);}
        if (ph < p0) {ph += (p1 - p0);}
        if (ph > p1) {ph -= (p1 - p0);}
        const double st = std::sin(th), ct = std::cos(th);
        const double sp = std::sin(ph), cp = std::cos(ph);
        for (int q = 0; q < n; ++q) {
          const int r = mp_h(m,k,j,q);
          mr_h(m,k,j,q) = d[4*r]*st*cp + d[4*r+1]*st*sp + d[4*r+2]*ct;
        }
      }
    }
  }
  Kokkos::deep_copy(vgd_mr, mr_h);
  Kokkos::deep_copy(vgd_wall, w_h);
  Kokkos::deep_copy(vgd_map, mp_h);
  const bool first = (vgd_alpha < 0.0);
  vgd_alpha = alpha;
  if (first && global_variable::my_rank == 0) {
    std::cout << "<rad_m1> vet_gd: global-direction SC closure on the sp wedge, "
              << ((vgd_gl_nmu > 0) ? ("product GL" + std::to_string(vgd_gl_nmu) + "x"
                                      + std::to_string(vgd_gl_nphi)
                                      + (vgd_gl_stag ? " staggered" : ""))
                  : ((vgd_ls > 0) ? ("level-symmetric LQ" + std::to_string(vgd_ls))
                                  : ("HEALPix nside " + std::to_string(vgd_nside))))
              << " (" << n << " directions); z-rotation every " << vgd_rot
              << " cycle(s) (0 = fixed), first angle " << alpha << ";"
              << " rank 0 wall ghost columns " << nwall << ", direction maps exact "
              << nexact << " of " << nmap << " (the rest: nearest direction)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdHaloInit
//! \brief the 8 lateral neighbour slots of every local MeshBlock (one x1 block, uniform
//! level; theta and phi periodic: the wrap of the logical location), message buffers

void RadiationM1::VetGdHaloInit() {
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
  std::vector<int> gof(nb2*nb3, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    if (l.lx1 != 0) {VgdFatal("needs ONE MeshBlock along x1");}
    gof[l.lx3*nb2 + l.lx2] = g;
  }
  for (int t = 0; t < nb2*nb3; ++t) {
    if (gof[t] < 0) {VgdFatal("incomplete MeshBlock grid (SMR/AMR not supported)");}
  }
  Kokkos::realloc(vgd_hloc, 8*nmb);
  auto hl_h = Kokkos::create_mirror_view(vgd_hloc);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  vgd_hmpi = false;
  for (int m = 0; m < nmb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    for (int o = 0; o < 8; ++o) {
      const int oo = (o < 4) ? o : (o + 1);       // skip the centre of the 3x3
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int l2 = (l.lx2 + dj + nb2) % nb2;
      const int l3 = (l.lx3 + dk + nb3) % nb3;
      const int g = gof[l3*nb2 + l2];
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
  // the largest region (an edge slot: nx x W) times the most variables per message
  // (all directions of one shell, or ln chi and ln S of every shell)
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nreg = std::max(indcs.nx2, indcs.nx3)*vgd_w;
  vgd_maxcnt = nreg*std::max(vgd_n, 2*c1);
  if (vgd_hmpi) {
    Kokkos::realloc(vgd_sbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    Kokkos::realloc(vgd_rbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    // the message plan: pieces sorted by (rank, receiver block, receiver slot) on both
    // sides; sizes per piece kn*jn with kn, jn in {nx, ws}
    std::vector<std::array<int, 4>> ps, pr;
    for (int m = 0; m < nmb; ++m) {
      for (int o = 0; o < 8; ++o) {
        if (hl_h(8*m + o) >= 0) {continue;}
        ps.push_back({vgd_hrank[8*m + o], vgd_hlid[8*m + o], 7 - o, 8*m + o});
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
          // the piece size: the receiver's slot (7 - o on the send side is the same
          // shape as my slot o)
          off += sz(o, ws);
          ((side == 0) ? vgd_pscnt : vgd_prcnt)[ws][p] += sz(o, ws);
        }
      }
    }
    Kokkos::deep_copy(vgd_soff, so_h);
    Kokkos::deep_copy(vgd_roff, ro_h);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdHalo
//! \brief fill the lateral band (vgd_w cells) of a band-index array a(m, v, k, j, i),
//! variables [0, nv), shells i in [i0, i1], from the 8 lateral neighbours' interiors:
//! a local neighbour is copied directly (it reads only interiors, it writes only bands:
//! no race), a remote one by one message per (block, slot)

void RadiationM1::VetGdHalo(DvceArray5D<Real> &a, const int nv, const int i0,
                            const int i1, const int ws, const bool mapd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int ni = i1 - i0 + 1;
  const int mx = std::max(nx2, nx3)*ws*nv*ni;
  auto hl_ = vgd_hloc;
  auto a_ = a;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  // region of slot o along one direction d = -1, 0, +1 with depth ws: (dest start, len,
  // neighbour source start) = (w - ws, ws, w + nx - ws), (w, nx, w), (w + nx, ws, w);
  // wall ghosts (mapd) take direction mp(v) of the source: the wall re-indexing folded
  // into the copy (new(d) = old(map(d)))
  const bool mpi = vgd_hmpi;
  if (mapd && !mpi) {return;}
  auto sb_ = vgd_sbuf;
  auto so_ = vgd_soff;
  auto ro_ = vgd_roff;
  const int nvi = nv*ni;
  par_for("m1_vgd_halo", DevExeSpace(), 0, nmb1, 0, 7, 0, mx - 1,
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
    if (nl >= 0 && mapd) {return;}    // the sweep reads on-rank neighbours directly
    if (nl >= 0) {
      const int jd = (dj < 0) ? (w - ws) : ((dj == 0) ? w : (w + nx2));
      const int kd = (dk < 0) ? (w - ws) : ((dk == 0) ? w : (w + nx3));
      const int js_ = (dj < 0) ? (w + nx2 - ws) : w;
      const int ks_ = (dk < 0) ? (w + nx3 - ws) : w;
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      a_(m,v,kd+kk,jd+jj,i) = a_(nl,vs,ks_+kk,js_+jj,i);
    } else if (mpi) {
      // my interior that the remote neighbour needs in ITS opposite slot
      const int js2 = (dj > 0) ? (w + nx2 - ws) : w;
      const int ks2 = (dk > 0) ? (w + nx3 - ws) : w;
      sb_(static_cast<size_t>(so_(ws,8*m + o))*nvi + t) = a_(m,v,ks2+kk,js2+jj,i);
    }
  });
#if MPI_PARALLEL_ENABLED
  if (mpi) {
    Kokkos::fence();
    Kokkos::Timer tq;
    std::vector<MPI_Request> req;
    auto rb_ = vgd_rbuf;
    // one message per partner rank and direction
    for (size_t p = 0; p < vgd_prk.size(); ++p) {
      const int rk = vgd_prk[p];
      req.emplace_back();
      MPI_Irecv(rb_.data() + static_cast<size_t>(vgd_prdsp[ws][p])*nvi,
                vgd_prcnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, MPI_COMM_WORLD,
                &req.back());
      req.emplace_back();
      MPI_Isend(sb_.data() + static_cast<size_t>(vgd_pdsp[ws][p])*nvi,
                vgd_pscnt[ws][p]*nvi, MPI_ATHENA_REAL, rk, 7001, MPI_COMM_WORLD,
                &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    vgd_tmpi += tq.seconds();
    par_for("m1_vgd_unpack", DevExeSpace(), 0, nmb1, 0, 7, 0, mx - 1,
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
      const int vs = (mapd && wl_(m,kd+kk,jd+jj) != 0) ? mp_(m,kd+kk,jd+jj,v) : v;
      a_(m,v,kd+kk,jd+jj,i0+ii) = rb_(static_cast<size_t>(ro_(ws,8*m + o))*nvi
                                      + ((vs*kn + kk)*jn + jj)*ni + ii);
    });
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdWall
//! \brief re-index the wall ghost intensities after an exchange (local-frame periodicity)

void RadiationM1::VetGdWall(const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int c2 = indcs.nx2 + 2*vgd_w, c3 = indcs.nx3 + 2*vgd_w;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  auto vi_ = vgd_i;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  par_for("m1_vgd_wall", DevExeSpace(), 0, nmb1, 0, c3 - 1, 0, c2 - 1, i0, i1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (wl_(m,k,j) == 0) {return;}
    Real buf[M1_VGD_NMAX];
    for (int q = 0; q < n; ++q) {buf[q] = vi_(m,q,k,j,i);}
    for (int q = 0; q < n; ++q) {vi_(m,q,k,j,i) = buf[mp_(m,k,j,q)];}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdSweep
//! \brief one sweep (inward then outward) of every direction

void RadiationM1::VetGdSweep() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  // the band index space of vgd_i / vgd_cs / the tables: active cells at [wb, wb + nx)
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vlat_icut;
  auto cs_ = vgd_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vlat_cnt;
  auto mr_ = vgd_mr;
  // intensities outside the block (band index) of an ON-RANK neighbour are read from its
  // interior directly (direction re-indexed at a wall); remote ones from the band, filled
  // per shell by VetGdHalo
  auto hl_ = vgd_hloc;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const int nx2b = indcs.nx2, nx3b = indcs.nx3;
  const bool bandx = vgd_bandx;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const int jlo = 0, jhi = je + wb - 1;
  const int klo = 0, khi = ke + wb - 1;
  const Real twopi = 2.0*M_PI;
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    for (int q = 0; q < n1 - lcut; ++q) {
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      par_for("m1_vgd_shell", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, n - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
        auto rd = [&](const int kk, const int jj, const int ii) -> Real {
          const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
          const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
          if (ek == 0 && ej == 0) {return vi_(m,d,kk,jj,ii);}
          const int oo = 3*(ek + 1) + (ej + 1);
          const int nl = hl_(8*m + ((oo < 4) ? oo : (oo - 1)));
          if (nl < 0) {return vi_(m,d,kk,jj,ii);}
          const int dd = (wl_(m,kk,jj) != 0) ? mp_(m,kk,jj,d) : d;
          return vi_(nl,dd,kk - ek*nx3b,jj - ej*nx2b,ii);
        };
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        // the face midpoints (uniform in index space, as the bilinear reads assume)
        const Real th = mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
        const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        const Real mr = nx*st*cp + ny*st*sp + nz*ct;
        if (inw == (mr >= 0.0)) {return;}
        const Real r = cx1v(m,i);
        const Real zr = r*fabs(mr);
        const Real p2 = fmax(r*r - zr*zr, 0.0);
        const Real ch0 = exp(cs_(m,0,k,j,i));
        const Real s0 = exp(cs_(m,1,k,j,i));
        Real iv;
        if (inw && l == n1 - 1) {
          // from the vacuum top face (the end cell's S extrapolated, as vet_col_lat)
          const Real rtop = cx1f(m,ie+1);
          const Real ds = sqrt(fmax(rtop*rtop - p2, 0.0)) - zr;
          const Real stp = fmax(1.5*s0 - 0.5*exp(cs_(m,1,k,j,i-1)), 1.0e-300);
          Real ex, w0, wu;
          VgdW(ch0*ds, ex, w0, wu);
          iv = wu*stp + w0*s0;
        } else if (!inw && zr < 1.0e-12*r) {
          iv = s0;                     // tangent at the cell centre (measure zero)
        } else {
          int typ = 1;                 // 1 upwind shell, 2 inner face, 3 own shell
          int iu = i + 1;
          Real ds;
          if (inw) {
            const Real ru = cx1v(m,i+1);
            ds = sqrt(ru*ru - p2) - zr;
          } else {
            const Real rd = (l > lcut) ? cx1v(m,i-1) : cx1f(m,i);
            if (l > lcut && p2 <= rd*rd*(1.0 + 1.0e-13)) {
              iu = i - 1;
              ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
            } else if (l == lcut && p2 <= rd*rd) {
              typ = 2;
              iu = i;
              ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
            } else {
              typ = 3;
              iu = i;
              ds = 2.0*zr;
            }
          }
          // the upwind point, its lateral cell and bilinear fractions (index space)
          const Real xu = r*st*cp - ds*nx;
          const Real yu = r*st*sp - ds*ny;
          const Real zu = r*ct - ds*nz;
          const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
          const Real thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
          Real dph = atan2(yu, xu) - ph;
          dph -= twopi*floor((dph + M_PI)/twopi);
          const Real fj = j + (thu - th)/mbsize.d_view(m).dx2;
          const Real fk = k + dph/mbsize.d_view(m).dx3;
          int j0 = static_cast<int>(floor(fj));
          int k0 = static_cast<int>(floor(fk));
          Real uj = fj - j0, uk = fk - k0;
          bool clp = false;
          if (j0 < jlo) {j0 = jlo; uj = 0.0; clp = true;}
          if (j0 > jhi) {j0 = jhi; uj = 1.0; clp = true;}
          if (k0 < klo) {k0 = klo; uk = 0.0; clp = true;}
          if (k0 > khi) {k0 = khi; uk = 1.0; clp = true;}
          if (clp) {Kokkos::atomic_add(&cnt_(0), 1.0);}
          const int j1 = j0 + 1, k1 = k0 + 1;
          bool bdone = false;
          if (clp && bandx && typ == 1) {
            // vet_gd_band_exit: the upwind point lies beyond the ghost band.  Instead of
            // reading the clamped (wrong) column, stop the segment where the ray leaves
            // the band (index linear in the path fraction s) and take I, chi, S there,
            // bilinear on the band-edge ghost columns of the two shells i and iu (the
            // lagged inflow of the last sweep) and linear in r between them
            const Real gj = fj - j, gk = fk - k;
            Real s = 1.0;
            if (fj < jlo) {s = fmin(s, (jlo - j)/gj);}
            if (fj > jhi + 1) {s = fmin(s, (jhi + 1 - j)/gj);}
            if (fk < klo) {s = fmin(s, (klo - k)/gk);}
            if (fk > khi + 1) {s = fmin(s, (khi + 1 - k)/gk);}
            s = fmax(s*(1.0 - 1.0e-9), 0.0);
            const Real dsx = s*ds;
            const Real xx = r*st*cp - dsx*nx, yx = r*st*sp - dsx*ny, zx = r*ct - dsx*nz;
            const Real rx = sqrt(xx*xx + yx*yx + zx*zx);
            const Real thx = acos(fmin(fmax(zx/rx, -1.0), 1.0));
            Real dpx = atan2(yx, xx) - ph;
            dpx -= twopi*floor((dpx + M_PI)/twopi);
            const Real fjx = fmin(fmax(j + (thx - th)/mbsize.d_view(m).dx2, Real(jlo)),
                                  Real(jhi + 1));
            const Real fkx = fmin(fmax(k + dpx/mbsize.d_view(m).dx3, Real(klo)),
                                  Real(khi + 1));
            int jx0 = static_cast<int>(floor(fjx)), kx0 = static_cast<int>(floor(fkx));
            jx0 = (jx0 > jhi) ? jhi : jx0;
            kx0 = (kx0 > khi) ? khi : kx0;
            const Real ujx = fjx - jx0, ukx = fkx - kx0;
            const Real wr = fmin(fmax((rx - r)/(cx1v(m,iu) - r), 0.0), 1.0);
            const bool want_in = inw;
            Real wsum = 0.0, isum = 0.0;
            for (int q = 0; q < 4; ++q) {
              const int kk = (q < 2) ? kx0 : (kx0 + 1);
              const int jj = (q % 2 == 0) ? jx0 : (jx0 + 1);
              const Real wq = ((q < 2) ? (1.0 - ukx) : ukx)
                              *((q % 2 == 0) ? (1.0 - ujx) : ujx);
              const Real mq = mr_(m,kk,jj,d);
              if (((want_in && mq < 0.0) || (!want_in && mq >= 0.0)) && wq > 0.0) {
                wsum += wq;
                isum += wq*((1.0 - wr)*rd(kk,jj,i) + wr*rd(kk,jj,iu));
              }
            }
            if (wsum > 0.0) {
              auto linx = [&](const int c, const int ii) {
                return VgdLerp(VgdLerp(cs_(m,c,kx0,jx0,ii), cs_(m,c,kx0,jx0+1,ii), ujx),
                               VgdLerp(cs_(m,c,kx0+1,jx0,ii), cs_(m,c,kx0+1,jx0+1,ii),
                                       ujx), ukx);
              };
              const Real chx = exp((1.0 - wr)*linx(0,i) + wr*linx(0,iu));
              const Real sx = exp((1.0 - wr)*linx(1,i) + wr*linx(1,iu));
              Real ex, w0, wu;
              VgdW(0.5*(chx + ch0)*dsx, ex, w0, wu);
              iv = fmax((isum/wsum)*ex + wu*sx + w0*s0, 0.0);
              bdone = true;
            }
          }
          if (bdone) {
            // done: the shortened segment from the band edge
          } else if (typ == 2) {
            // the diffusion intensity at the inner face of the first shell
            const Real rf = cx1f(m,i);
            const Real muf = sqrt(fmax(1.0 - p2/(rf*rf), 0.0));
            const int i1 = i + 1;
            auto ib = [&](const int kk, const int jj) {
              const Real sa = exp(cs_(m,1,kk,jj,i)), sb = exp(cs_(m,1,kk,jj,i1));
              const Real g0 = -(sb - sa)/(cx1v(m,i1) - cx1v(m,i))/exp(cs_(m,0,kk,jj,i));
              return sa + g0*muf;
            };
            auto sbt = [&](const int kk, const int jj) {
              return fmax(1.5*exp(cs_(m,1,kk,jj,i)) - 0.5*exp(cs_(m,1,kk,jj,i1)),
                          1.0e-300);
            };
            Real ex, w0, wu;
            VgdW(ch0*ds, ex, w0, wu);
            const Real bv = fmax(VgdLerp(VgdLerp(ib(k0,j0), ib(k0,j1), uj),
                                         VgdLerp(ib(k1,j0), ib(k1,j1), uj), uk), 0.0);
            const Real qv = VgdLerp(VgdLerp(sbt(k0,j0), sbt(k0,j1), uj),
                                    VgdLerp(sbt(k1,j0), sbt(k1,j1), uj), uk);
            iv = bv*ex + wu*qv + w0*s0;
          } else {
            // each (cell, direction) holds ONE intensity, of the pass its n_d . r_hat
            // selects; the ray read here is inward at the upwind point for the inward
            // pass and for a turned ray, outward otherwise.  Stencil cells of the other
            // pass hold the wrong branch (near the tangent band): they are dropped and
            // the bilinear weights renormalised (all valid: the plain bilinear)
            const bool want_in = inw || (typ == 3);
            bool ok[4];
            bool all = true;
            for (int q = 0; q < 4; ++q) {
              const int kk = (q < 2) ? k0 : k1;
              const int jj = (q % 2 == 0) ? j0 : j1;
              const Real mq = mr_(m,kk,jj,d);
              ok[q] = want_in ? (mq < 0.0) : (mq >= 0.0);
              all = all && ok[q];
            }
            Real ivu;
            if (all) {
              ivu = VgdLerp(VgdLerp(rd(k0,j0,iu), rd(k0,j1,iu), uj),
                            VgdLerp(rd(k1,j0,iu), rd(k1,j1,iu), uj), uk);
            } else {
              Real ws = 0.0, vs = 0.0;
              for (int q = 0; q < 4; ++q) {
                const int kk = (q < 2) ? k0 : k1;
                const int jj = (q % 2 == 0) ? j0 : j1;
                const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                *((q % 2 == 0) ? (1.0 - uj) : uj);
                if (ok[q] && wq > 0.0) {ws += wq; vs += wq*rd(kk,jj,iu);}
              }
              ivu = (ws > 0.0) ? vs/ws : -1.0;   // none: the upwind source (below)
            }
            auto lin = [&](const int c) {
              return VgdLerp(VgdLerp(cs_(m,c,k0,j0,iu), cs_(m,c,k0,j1,iu), uj),
                             VgdLerp(cs_(m,c,k1,j0,iu), cs_(m,c,k1,j1,iu), uj), uk);
            };
            const Real chu = exp(lin(0));
            const Real su = exp(lin(1));
            Real ex, w0, wu;
            VgdW(0.5*(chu + ch0)*ds, ex, w0, wu);
            iv = fmax(((ivu < 0.0) ? su : ivu)*ex + wu*su + w0*s0, 0.0);
          }
        }
        vi_(m,d,k,j,i) = iv;
      });
      // the shell is complete on every block: its lateral band, exact (not lagged)
      if (vgd_time_halo) {Kokkos::fence();}
      Kokkos::Timer th;
      VetGdHalo(vgd_i, n, i, i, vgd_wsh[i], true);
      Kokkos::fence();
      vgd_thalo += th.seconds();
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdMoments
//! \brief D = K/J in the local (r, theta, phi) frame -> tau_ten LAT slots (relative to
//! vet_col's uniaxial tensor of this build; VetLatBuild folds LAT0 into slot 0)

void RadiationM1::VetGdMoments() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int ilo = is + vlat_icut;
  auto vi_ = vgd_i;
  const int og = vgd_w - indcs.ng;   // mesh index -> band index of vgd_i
  auto dir_ = vgd_dir;
  auto tt_ = tau_ten;
  auto &mbsize = pmy_pack->pmb->mb_size;
  // vet_gd_replace with vet_col_surface_q: the outer Marshak q = H_r/J of the top cell
  const bool repq = vgd_replace && vcol_sq;
  auto vq_ = vcol_q;
  const Real qlo = vcol_qmin, qhi = vcol_qmax;
  par_for("m1_vgd_mom", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int c = 0; c < M1_TT_NLAT; ++c) {tt_(m,M1_TT_LAT0+c,k,j,i) = 0.0;}
    if (i < ilo) {return;}
    const Real th = mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
    const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
    const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
    Real jm = 0.0, rr = 0.0, rt = 0.0, rp = 0.0, tq = 0.0, tp = 0.0, pq = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      const Real a = nx*st*cp + ny*st*sp + nz*ct;
      const Real b = nx*ct*cp + ny*ct*sp - nz*st;
      const Real c = -nx*sp + ny*cp;
      const Real wi = dir_(d,3)*vi_(m,d,k+og,j+og,i);
      jm += wi;
      rr += wi*a*a;
      rt += wi*a*b;
      rp += wi*a*c;
      tq += wi*b*b;
      tp += wi*b*c;
      pq += wi*c*c;
    }
    if (jm > 0.0) {
      const Real fk = tt_(m,0,k,j,i);
      tt_(m,M1_TT_LAT0,k,j,i) = rr/jm - fk;
      tt_(m,M1_TT_LAT0+1,k,j,i) = rt/jm;
      tt_(m,M1_TT_LAT0+2,k,j,i) = rp/jm;
      tt_(m,M1_TT_LAT0+3,k,j,i) = 0.5*(tq - pq)/jm;    // D_tt - (1 - D_rr)/2
      tt_(m,M1_TT_LAT0+4,k,j,i) = tp/jm;
      tt_(m,M1_TT_LAT0+5,k,j,i) = -0.5*(tq - pq)/jm;   // D_pp - (1 - D_rr)/2
    }
  });
  // vet_gd_thin_taumin > 0: taper the LATERAL parts (LAT1..LAT5: D_r,lat and the
  // tangential anisotropy, both lagged terms of the face equations) to 0 in the far thin
  // top, log-linearly from 1 at tau_top = 10 taumin to 0 at tau_top = taumin (tau_top of
  // the column above the cell centre).  D_rr (LAT0, implicit) is kept.  The accuracy
  // region of the closure is below the photosphere; the lagged terms lose their fixed
  // point in thin cells at c dt >> dx (C1_RESULTS.md, GD_RESULTS.md sect. 5)
  if (vgd_taumin > 0.0) {
    auto cs_ = vlat_cs;
    auto cx1f_ = pmy_pack->pcoord->xx1f;
    const Real tlo = vgd_taumin;
    // vet_gd_thin_decades (default 1): the ramp from tau = taumin to taumin 10^decades;
    // vet_gd_thin_smooth (default false): smoothstep 3x^2 - 2x^3 instead of linear in
    // ln tau; vet_gd_thin_parts (bit 1 = D_r,lat LAT1-2, bit 2 = tangential LAT3-5,
    // default 3 = both)
    const Real dec = vgd_tdec;
    const Real lgw = dec*log(10.0);
    const Real thi = tlo*pow(10.0, dec);
    const bool smo = vgd_tsmooth;
    const int prt = vgd_tparts;
    par_for("m1_vgd_taper", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      Real tau = 0.0;
      for (int i = ie; i >= ilo; --i) {
        const Real dtc = exp(cs_(m,0,k,j,i))*(cx1f_(m,i+1) - cx1f_(m,i));
        const Real tc = tau + 0.5*dtc;
        tau += dtc;
        Real w = (tc <= tlo) ? 0.0 : ((tc >= thi) ? 1.0 : log(tc/tlo)/lgw);
        if (smo) {w = w*w*(3.0 - 2.0*w);}
        if (w < 1.0) {
          for (int c = 1; c < M1_TT_NLAT; ++c) {
            const int bit = (c <= 2) ? 1 : 2;
            if (prt & bit) {tt_(m,M1_TT_LAT0+c,k,j,i) *= w;}
          }
        }
      }
    });
  }
  if (!repq) {return;}
  // vet_gd_replace with vet_col_surface_q: the outer Marshak q in vet_col's FACE form,
  // q = H(face)/J_f, clamped to [vet_col_surface_qmin, _qmax]: H(face) = H_r of the top
  // cell carried to the face as r^2 H (free streaming over the half cell), J_f = J of the
  // top cell, or with vet_col_surface_face the same extrapolation as vet_col
  // (jca J_top - jcb J_below, limited to [0, jcap J_top])
  const bool sqf = vcol_sqf;
  const Real jca = vcol_jca, jcb = vcol_jcb, jcap = vcol_jcap, q0 = marshak_q;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  par_for("m1_vgd_q", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real th = mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
    const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
    const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
    Real jt = 0.0, jb = 0.0, ht = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real a = dir_(d,0)*st*cp + dir_(d,1)*st*sp + dir_(d,2)*ct;
      const Real w = dir_(d,3);
      jt += w*vi_(m,d,k+og,j+og,ie);
      jb += w*vi_(m,d,k+og,j+og,ie-1);
      ht += w*a*vi_(m,d,k+og,j+og,ie);
    }
    const Real rr = cx1v(m,ie)/cx1f(m,ie+1);
    const Real hf = ht*rr*rr;
    const Real jq = sqf ? fmin(fmax(jca*jt - jcb*jb, 0.0), jcap*jt) : jt;
    vq_(m,k,j) = (jq > 0.0) ? fmin(fmax(hf/jq, qlo), qhi) : q0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdSmooth
//! \brief vet_gd_smooth = N > 0 (diagnostic, default 0): N passes of a lateral 1-2-1
//! filter (theta, then phi) of the D_rr correction LAT0, ghosts exchanged before each
//! direction (periodic walls: scalar copy)

void RadiationM1::VetGdSmooth() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ilo = is + vlat_icut;
  if (vgd_sm.extent(0) == 0) {
    Kokkos::realloc(vgd_sm, pmy_pack->nmb_thispack, indcs.nx3 + 2*indcs.ng,
                    indcs.nx2 + 2*indcs.ng, indcs.nx1 + 2*indcs.ng);
  }
  auto tt_ = tau_ten;
  auto sm_ = vgd_sm;
  const int c0 = M1_TT_LAT0;
  for (int pass = 0; pass < vgd_smooth; ++pass) {
    for (int dir = 0; dir < 2; ++dir) {
      VetLatTTGhosts();
      par_for("m1_vgd_sm", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (dir == 0) {
          sm_(m,k,j,i) = 0.25*(tt_(m,c0,k,j-1,i) + 2.0*tt_(m,c0,k,j,i)
                               + tt_(m,c0,k,j+1,i));
        } else {
          sm_(m,k,j,i) = 0.25*(tt_(m,c0,k-1,j,i) + 2.0*tt_(m,c0,k,j,i)
                               + tt_(m,c0,k+1,j,i));
        }
      });
      par_for("m1_vgd_smc", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        tt_(m,c0,k,j,i) = sm_(m,k,j,i);
      });
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdBuild
//! \brief the gd build inside VetLatBuild: source and first shell (VetLatSweep(0)), the
//! sweep(s) with the lagged lateral inflow, the moments; timed per part (fenced)

void RadiationM1::VetGdBuild() {
  Kokkos::Timer tm;
  VetLatSweep(0);
  {
    // ln chi, ln S in the band index space, the lateral band filled exactly
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int og = vgd_w - indcs.ng;
    const int ilo = indcs.is + vlat_icut, ie = indcs.ie;
    auto cs_ = vlat_cs;
    auto cw_ = vgd_cs;
    par_for("m1_vgd_csw", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks,
            indcs.ke, indcs.js, indcs.je, ilo, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      cw_(m,0,k+og,j+og,i) = cs_(m,0,k,j,i);
      cw_(m,1,k+og,j+og,i) = cs_(m,1,k,j,i);
    });
    VetGdHalo(vgd_cs, 2, ilo, ie, vgd_w, false);
  }
  Kokkos::fence();
  vgd_tsrc += tm.seconds();
  // vet_gd_rotate_every: a new z-angle at the first build of a block of N cycles (the
  // build sits outside the Picard solve, once per step); the stored intensities belong
  // to the old directions, so the lagged inflow is re-converged as at the first build
  const Real al = VetGdAngle(pmy_pack->pmesh->ncycle);
  bool rotated = false;
  if (al != vgd_alpha) {
    Kokkos::Timer tt;
    VetGdTables(al);
    vgd_ttab += tt.seconds();
    rotated = true;
  }
  // the per-shell halo makes one sweep exact; vet_gd_iter > 1 only repeats it
  (void) rotated;
  const int nit = vgd_iter;
  for (int it = 0; it < nit; ++it) {
    tm.reset();
    VetGdSweep();
    Kokkos::fence();
    vgd_tswp += tm.seconds();
    vlat_ncall += 1.0;
  }
  tm.reset();
  VetGdMoments();
  if (vgd_smooth > 0) {VetGdSmooth();}
  VetLatOdMax();   // D_r,lat in the implicit operator (vet_col_lat_offdiag = operator)
  Kokkos::fence();
  vgd_tmom += tm.seconds();
  auto cnt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vlat_cnt);
  vlat_nclamp = cnt_h(0);
  // debug (gate): VGD_DUMP_I=<file>: rank 0's intensities incl. ghosts, first build
  const char *fi = std::getenv("VGD_DUMP_I");
  if (fi != nullptr && vlat_nbuild == 0 && global_variable::my_rank == 0) {
    auto vi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i);
    FILE *fp = std::fopen(fi, "wb");
    if (fp != nullptr) {
      int64_t sh[5];
      for (int c = 0; c < 5; ++c) {sh[c] = static_cast<int64_t>(vi_h.extent(c));}
      std::fwrite(sh, sizeof(int64_t), 5, fp);
      std::fwrite(vi_h.data(), sizeof(Real), vi_h.size(), fp);
      std::fclose(fp);
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdMms
//! \brief <rad_m1>/vet_gd_mms = true (test only; prints and exits): manufactured-solution
//! check of the discrete div(E D) the implicit transport assembles on the sp wedge, with
//! a smooth E and an anisotropic D whose six components are all nonzero.  The discrete
//! terms are evaluated with the SAME helpers and face formulas as the kernels
//! (M1SphDrr on x1 faces, M1EddDiag + the vet_gd deviation a on x2/x3 faces, M1SphLat,
//! M1SphTan, cell values averaged to the face), on host copies of the coordinates, and
//! compared per term with the analytic physical components
//!   r:  d_r P_rr + (2 P_rr - P_tt - P_pp)/r                         [R1 x1 face diag]
//!       (1/r) d_th P_rt + cot P_rt/r + (1/(r sin)) d_ph P_rp        [R2 M1SphLat d 0]
//!   th: (1/r) d_th P_tt                                             [T1 x2 face diag]
//!       (1/r^3) d_r (r^3 P_rt)                                      [T2 M1SphLat d 1]
//!       cot (P_tt - P_pp)/r + (1/(r sin)) d_ph P_tp                 [T3 M1SphTan d 1]
//!   ph: (1/(r sin)) d_ph P_pp                                       [P1 x3 face diag]
//!       (1/r^3) d_r (r^3 P_rp)                                      [P2 M1SphLat d 2]
//!       (1/r) d_th P_tp + 2 cot P_tp/r                              [P3 M1SphTan d 2]
//! (P_tt + P_pp = (1 - D_rr) E: the deviation a has no radial part.)  Errors: max over
//! the faces at least 2 cells from every block edge, relative to the max |component|.

void RadiationM1::VetGdMms() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int c1 = indcs.nx1 + 2*ng, c2 = indcs.nx2 + 2*ng, c3 = indcs.nx3 + 2*ng;
  const int nmb = pmy_pack->nmb_thispack;
  auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x1v);
  auto x2v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x2v);
  auto x3v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x3v);
  auto x1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx1f);
  auto &dxf = pmy_pack->pcoord->dxface;
  auto d1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dxf.x1f);
  auto d2f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dxf.x2f);
  auto d3f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dxf.x3f);
  auto &mbs = pmy_pack->pmb->mb_size;
  const double r0 = pm->mesh_size.x1min, r1 = pm->mesh_size.x1max;
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  const double p0 = pm->mesh_size.x3min, p1 = pm->mesh_size.x3max;
  // the manufactured fields (normalised coordinates u, v, w in [0, 1])
  auto U = [&](double r) {return (r - r0)/(r1 - r0);};
  auto V = [&](double t) {return (t - t0)/(t1 - t0);};
  auto W = [&](double p) {return (p - p0)/(p1 - p0);};
  auto fE = [&](double r, double t, double p) {
    return (1.0 + 0.5*U(r)*U(r))*(1.0 + 0.2*std::sin(3.0*V(t))*std::cos(2.0*W(p)));
  };
  auto fDrr = [&](double r, double t, double p) {
    return 0.5 + 0.1*std::sin(2.0*V(t) + 1.0)*std::cos(W(p))*(1.0 + U(r));
  };
  auto fA = [&](double r, double t, double p) {
    return 0.05*std::cos(V(t))*std::sin(2.0*W(p) + 0.3)*(1.0 + 0.5*U(r));
  };
  auto fDrt = [&](double r, double t, double p) {
    return 0.08*std::sin(V(t) + 0.2)*std::cos(W(p))*(1.0 + U(r)*U(r));
  };
  auto fDrp = [&](double r, double t, double p) {
    return 0.06*std::cos(2.0*V(t))*std::sin(W(p) + 0.4)*(1.0 + U(r));
  };
  auto fDtp = [&](double r, double t, double p) {
    return 0.04*std::sin(3.0*V(t) + 0.1)*std::cos(2.0*W(p))*(1.0 + 0.7*U(r));
  };
  // physical P components
  auto Prr = [&](double r, double t, double p) {return fDrr(r,t,p)*fE(r,t,p);};
  auto Ptt = [&](double r, double t, double p) {
    return (0.5*(1.0 - fDrr(r,t,p)) + fA(r,t,p))*fE(r,t,p);
  };
  auto Ppp = [&](double r, double t, double p) {
    return (0.5*(1.0 - fDrr(r,t,p)) - fA(r,t,p))*fE(r,t,p);
  };
  auto Prt = [&](double r, double t, double p) {return fDrt(r,t,p)*fE(r,t,p);};
  auto Prp = [&](double r, double t, double p) {return fDrp(r,t,p)*fE(r,t,p);};
  auto Ptp = [&](double r, double t, double p) {return fDtp(r,t,p)*fE(r,t,p);};
  // 4th-order central derivative of an analytic function along one coordinate
  auto der = [&](auto f, int c, double r, double t, double p) {
    const double h = 1.0e-3*((c == 0) ? (r1 - r0) : ((c == 1) ? (t1 - t0) : (p1 - p0)));
    auto g = [&](double s) {
      return (c == 0) ? f(r + s, t, p) : ((c == 1) ? f(r, t + s, p) : f(r, t, p + s));
    };
    return (8.0*(g(h) - g(-h)) - (g(2.0*h) - g(-2.0*h)))/(12.0*h);
  };
  // host work arrays as the kernels read them
  HostArray5D<Real> iwh("mms_iw", nmb, M1_NIW, c3, c2, c1);
  HostArray5D<Real> tth("mms_tt", nmb, M1_TT_LAT0 + M1_TT_NLAT, c3, c2, c1);
  const int c0 = M1_TT_LAT0;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        for (int i = 0; i < c1; ++i) {
          const double r = x1v(m,i), t = x2v(m,j), p = x3v(m,k);
          iwh(m,M1_IW_EP,k,j,i) = fE(r,t,p);
          iwh(m,M1_IW_WCHI,k,j,i) = fDrr(r,t,p);
          iwh(m,M1_IW_N1,k,j,i) = 1.0;
          iwh(m,M1_IW_N1+1,k,j,i) = 0.0;
          iwh(m,M1_IW_N1+2,k,j,i) = 0.0;
          tth(m,c0+1,k,j,i) = fDrt(r,t,p);
          tth(m,c0+2,k,j,i) = fDrp(r,t,p);
          tth(m,c0+3,k,j,i) = fA(r,t,p);
          tth(m,c0+4,k,j,i) = fDtp(r,t,p);
          tth(m,c0+5,k,j,i) = -fA(r,t,p);
        }
      }
    }
  }
  const int il = 0, iu = c1 - 1, jl = 0, ju = c2 - 1, kl = 0, ku = c3 - 1;
  const int ep = M1_IW_EP;
  double emax[8] = {0.0}, amax[3] = {0.0};
  const char *nm[8] = {"R1 x1-face diag", "R2 SphLat d0", "T1 x2-face diag",
                       "T2 SphLat d1", "T3 SphTan d1", "P1 x3-face diag", "P2 SphLat d2",
                       "P3 SphTan d2"};
  auto lat = [&](int m, int d, int k, int j, int i) {
    return M1SphLat(iwh, tth, c0, x1v, x2v, x3v, m, d, k, j, i, true, il, iu, jl, ju, kl,
                    ku, ep);
  };
  auto tan = [&](int m, int d, int k, int j, int i) {
    return M1SphTan(iwh, tth, c0, x1v, x2v, x3v, m, d, k, j, i, true, jl, ju, kl, ku, ep);
  };
  for (int m = 0; m < nmb; ++m) {
    for (int k = ks + 2; k <= ke - 2; ++k) {
      for (int j = js + 2; j <= je - 2; ++j) {
        for (int i = is + 2; i <= ie - 2; ++i) {
          const double r = x1v(m,i), t = x2v(m,j), p = x3v(m,k);
          // x1 face i+1/2 (between i and ip), as the S2 row
          {
            const int ip = i + 1;
            const double rf = x1f(m,ip);
            const double si = SQR(x1v(m,i)/rf), sp = SQR(x1v(m,ip)/rf);
            const double wi = M1SphDrr(iwh(m,M1_IW_WCHI,k,j,i), 1.0, si);
            const double wp = M1SphDrr(iwh(m,M1_IW_WCHI,k,j,ip), 1.0, sp);
            const double g1 = (wp*iwh(m,ep,k,j,ip) - wi*iwh(m,ep,k,j,i))/d1f(m,k,j,ip);
            const double g2 = 0.5*(lat(m,0,k,j,i) + lat(m,0,k,j,ip));
            const double a1 = der(Prr,0,rf,t,p) + (2.0*Prr(rf,t,p) - Ptt(rf,t,p)
                                                    - Ppp(rf,t,p))/rf;
            const double a2 = der(Prt,1,rf,t,p)/rf
                              + std::cos(t)/std::sin(t)*Prt(rf,t,p)/rf
                              + der(Prp,2,rf,t,p)/(rf*std::sin(t));
            emax[0] = std::max(emax[0], std::fabs(g1 - a1));
            emax[1] = std::max(emax[1], std::fabs(g2 - a2));
            amax[0] = std::max(amax[0], std::fabs(a1 + a2));
          }
          // x2 face j-1/2 (between jm and j), as m1_impl_f2face
          {
            const int jm = j - 1;
            const double tf = mbs.h_view(m).x2min + (j - js)*mbs.h_view(m).dx2;
            const double dl = M1EddDiag(iwh(m,M1_IW_WCHI,k,jm,i), 0.0)
                              + tth(m,c0+3,k,jm,i);
            const double dr = M1EddDiag(iwh(m,M1_IW_WCHI,k,j,i), 0.0)
                              + tth(m,c0+3,k,j,i);
            const double g1 = (dr*iwh(m,ep,k,j,i) - dl*iwh(m,ep,k,jm,i))/d2f(m,k,j,i);
            const double g2 = 0.5*(lat(m,1,k,jm,i) + lat(m,1,k,j,i));
            const double g3 = 0.5*(tan(m,1,k,jm,i) + tan(m,1,k,j,i));
            auto r3prt = [&](double rr, double tt, double pp) {
              return rr*rr*rr*Prt(rr,tt,pp);
            };
            const double ct = std::cos(tf)/std::sin(tf);
            const double a1 = der(Ptt,1,r,tf,p)/r;
            const double a2 = der(r3prt,0,r,tf,p)/(r*r*r);
            const double a3 = ct*(Ptt(r,tf,p) - Ppp(r,tf,p))/r
                              + der(Ptp,2,r,tf,p)/(r*std::sin(tf));
            emax[2] = std::max(emax[2], std::fabs(g1 - a1));
            emax[3] = std::max(emax[3], std::fabs(g2 - a2));
            emax[4] = std::max(emax[4], std::fabs(g3 - a3));
            amax[1] = std::max(amax[1], std::fabs(a1 + a2 + a3));
          }
          // x3 face k-1/2 (between km and k), as m1_impl_f3face
          {
            const int km = k - 1;
            const double pf = mbs.h_view(m).x3min + (k - ks)*mbs.h_view(m).dx3;
            const double dl = M1EddDiag(iwh(m,M1_IW_WCHI,km,j,i), 0.0)
                              + tth(m,c0+5,km,j,i);
            const double dr = M1EddDiag(iwh(m,M1_IW_WCHI,k,j,i), 0.0)
                              + tth(m,c0+5,k,j,i);
            const double g1 = (dr*iwh(m,ep,k,j,i) - dl*iwh(m,ep,km,j,i))/d3f(m,k,j,i);
            const double g2 = 0.5*(lat(m,2,km,j,i) + lat(m,2,k,j,i));
            const double g3 = 0.5*(tan(m,2,km,j,i) + tan(m,2,k,j,i));
            auto r3prp = [&](double rr, double tt, double pp) {
              return rr*rr*rr*Prp(rr,tt,pp);
            };
            const double ct = std::cos(t)/std::sin(t);
            const double a1 = der(Ppp,2,r,t,pf)/(r*std::sin(t));
            const double a2 = der(r3prp,0,r,t,pf)/(r*r*r);
            const double a3 = der(Ptp,1,r,t,pf)/r + 2.0*ct*Ptp(r,t,pf)/r;
            emax[5] = std::max(emax[5], std::fabs(g1 - a1));
            emax[6] = std::max(emax[6], std::fabs(g2 - a2));
            emax[7] = std::max(emax[7], std::fabs(g3 - a3));
            amax[2] = std::max(amax[2], std::fabs(a1 + a2 + a3));
          }
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, emax, 8, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, amax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    std::cout << "VGD_MMS nx = " << pm->mesh_indcs.nx1 << " " << pm->mesh_indcs.nx2
              << " " << pm->mesh_indcs.nx3 << std::endl;
    for (int t = 0; t < 8; ++t) {
      const int c = (t < 2) ? 0 : ((t < 5) ? 1 : 2);
      std::cout << "VGD_MMS " << nm[t] << " max|err|/max|div P_" << c << "| = "
                << emax[t]/amax[c] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdRealDiag
//! \brief DIAGNOSTIC (env VGD_RDIAG = N, every N-th implicit solve; nothing written to
//! the state): per band of the column optical depth tau_top (from the gd build's ln chi):
//!  (1) realizability of the closure against the new state: lambda_min of D - f f^T,
//!      D the tensor of tau_ten (D_rr slot 0; D_tt, D_pp = (1 - D_rr)/2 +- a; D_rt, D_rp,
//!      D_tp the LAT slots), f = F/(c E) (u0, physical components); |f| > 1 cells;
//!  (2) the sign structure of the last pass's 19-point stencil: rows with a POSITIVE
//!      off-diagonal (not an M-matrix row) and max sum(positive off-diag)/diag;
//!  plus the write-back |F| > c E scale-back counters (M1_POS_FCLIP, FCLIPM) since the
//!  last report.

void RadiationM1::VetGdRealDiag() {
  if (vgd_rdiag < 0) {
    const char *ev = std::getenv("VGD_RDIAG");
    vgd_rdiag = (ev != nullptr) ? std::max(0, std::atoi(ev)) : 0;
  }
  if (vgd_rdiag <= 0) {return;}
  if ((++vgd_rdcnt % vgd_rdiag) != 0) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto tt_ = tau_ten;
  auto u0_ = u0;
  auto cs_ = vlat_cs;
  auto cx1f_ = pmy_pack->pcoord->xx1f;
  auto st_ = ost;
  const bool hst = impl_stencil && (ost.extent(0) > 0);
  const int nst = hst ? static_cast<int>(ost.extent(1)) : 0;
  const Real cl = c_light;
  const int c0 = M1_TT_LAT0;
  constexpr int NB = 5, NQ = 8;
  // per band: ncell, n(lmin<0), min lmin, n(|f|>1), n(pos rows), max ratio, min Drr
  DvceArray2D<Real> acc("vgd_rd", NB, NQ);
  {
    auto a0 = Kokkos::create_mirror_view(acc);
    for (int b = 0; b < NB; ++b) {
      for (int q = 0; q < NQ; ++q) {a0(b,q) = (q == 2 || q == 6) ? 1.0e300 : 0.0;}
    }
    Kokkos::deep_copy(acc, a0);
  }
  auto ac_ = acc;
  par_for("m1_vgd_rdiag", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real tau = 0.0;
    for (int i = ie; i >= is; --i) {
      const Real dtc = exp(cs_(m,0,k,j,i))*(cx1f_(m,i+1) - cx1f_(m,i));
      const Real tc = tau + 0.5*dtc;
      tau += dtc;
      const int b = (tc >= 1.0) ? 0 : ((tc >= 0.1) ? 1 : ((tc >= 0.01) ? 2
                    : ((tc >= 1.0e-3) ? 3 : 4)));
      const Real drr = tt_(m,0,k,j,i);
      const Real a = tt_(m,c0+3,k,j,i);
      Real d[3][3];
      d[0][0] = drr;
      d[1][1] = 0.5*(1.0 - drr) + a;
      d[2][2] = 0.5*(1.0 - drr) - a;
      d[0][1] = d[1][0] = tt_(m,c0+1,k,j,i);
      d[0][2] = d[2][0] = tt_(m,c0+2,k,j,i);
      d[1][2] = d[2][1] = tt_(m,c0+4,k,j,i);
      const Real e = u0_(m,M1_E,k,j,i);
      Real f[3] = {0.0, 0.0, 0.0};
      if (e > 0.0) {
        for (int q = 0; q < 3; ++q) {f[q] = u0_(m,M1_F1+q,k,j,i)/(cl*e);}
      }
      Real mm[3][3];
      for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) {mm[p][q] = d[p][q] - f[p]*f[q];}
      }
      // smallest eigenvalue of the symmetric 3x3 (trigonometric closed form)
      const Real p1 = mm[0][1]*mm[0][1] + mm[0][2]*mm[0][2] + mm[1][2]*mm[1][2];
      const Real qm = (mm[0][0] + mm[1][1] + mm[2][2])/3.0;
      Real lmin;
      if (p1 < 1.0e-300) {
        lmin = fmin(mm[0][0], fmin(mm[1][1], mm[2][2]));
      } else {
        const Real p2 = (mm[0][0] - qm)*(mm[0][0] - qm) + (mm[1][1] - qm)*(mm[1][1] - qm)
                        + (mm[2][2] - qm)*(mm[2][2] - qm) + 2.0*p1;
        const Real pp = sqrt(p2/6.0);
        Real bm[3][3];
        for (int p = 0; p < 3; ++p) {
          for (int q = 0; q < 3; ++q) {bm[p][q] = (mm[p][q] - ((p == q) ? qm : 0.0))/pp;}
        }
        const Real det = bm[0][0]*(bm[1][1]*bm[2][2] - bm[1][2]*bm[2][1])
                         - bm[0][1]*(bm[1][0]*bm[2][2] - bm[1][2]*bm[2][0])
                         + bm[0][2]*(bm[1][0]*bm[2][1] - bm[1][1]*bm[2][0]);
        const Real r = fmin(fmax(0.5*det, -1.0), 1.0);
        const Real phi = acos(r)/3.0;
        lmin = qm + 2.0*pp*cos(phi + 2.0*M_PI/3.0);
      }
      const Real fn = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
      Kokkos::atomic_add(&ac_(b,0), 1.0);
      if (lmin < -1.0e-12) {Kokkos::atomic_add(&ac_(b,1), 1.0);}
      Kokkos::atomic_min(&ac_(b,2), lmin);
      if (fn > 1.0 + 1.0e-12) {Kokkos::atomic_add(&ac_(b,3), 1.0);}
      // the LATERAL part of f (|f_lat| > 0.5: the flux mostly sideways)
      if (f[1]*f[1] + f[2]*f[2] > 0.25) {Kokkos::atomic_add(&ac_(b,7), 1.0);}
      Kokkos::atomic_min(&ac_(b,6), drr - 1.0/3.0);
      if (hst) {
        const Real dg = st_(m,0,k,j,i);
        Real pos = 0.0;
        for (int q = 1; q < nst; ++q) {pos += fmax(st_(m,q,k,j,i), 0.0);}
        if (pos > 0.0) {
          Kokkos::atomic_add(&ac_(b,4), 1.0);
          Kokkos::atomic_max(&ac_(b,5), pos/fmax(fabs(dg), 1.0e-300));
        }
      }
    }
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acc);
  Real loc[NB*NQ];
  for (int b = 0; b < NB; ++b) {
    for (int q = 0; q < NQ; ++q) {loc[b*NQ + q] = h(b,q);}
  }
  auto pch = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
  Real fc[2] = {pch(M1_POS_FCLIP), pch(M1_POS_FCLIPM)};
#if MPI_PARALLEL_ENABLED
  {
    Real g[NB*NQ];
    Real gs[NB*NQ];
    for (int t = 0; t < NB*NQ; ++t) {gs[t] = loc[t];}
    MPI_Allreduce(gs, g, NB*NQ, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    Real gmn[NB*NQ], gmx[NB*NQ];
    MPI_Allreduce(gs, gmn, NB*NQ, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(gs, gmx, NB*NQ, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    for (int b = 0; b < NB; ++b) {
      for (int q = 0; q < NQ; ++q) {
        const int t = b*NQ + q;
        loc[t] = (q == 2 || q == 6) ? gmn[t] : ((q == 5) ? gmx[t] : g[t]);
      }
    }
    Real fg[2];
    MPI_Allreduce(fc, fg, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    fc[0] = fg[0];
    fc[1] = fg[1];
  }
#endif
  if (global_variable::my_rank == 0) {
    const char *bn[NB] = {">1", "0.1-1", "0.01-0.1", "1e-3-0.01", "<1e-3"};
    std::cout << "VGD_RDIAG cycle=" << pmy_pack->pmesh->ncycle << " time="
              << pmy_pack->pmesh->time << " Fclip(cum)=" << fc[0] << " sum(|f|-1)="
              << fc[1];
    for (int b = 0; b < NB; ++b) {
      const Real *v = &loc[b*NQ];
      std::cout << " | tau " << bn[b] << ": n=" << v[0] << " lmin<0 " << v[1]
                << " lmin " << v[2] << " |f|>1 " << v[3] << " posrow " << v[4]
                << " maxpos/diag " << v[5] << " min(Drr-1/3) " << v[6] << " |flat|>0.5 "
                << v[7];
    }
    std::cout << std::endl;
  }
}

} // namespace radm1
