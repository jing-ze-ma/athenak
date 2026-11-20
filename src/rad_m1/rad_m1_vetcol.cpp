//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetcol.cpp
//! \brief <rad_m1>/closure = vet_col: a VARIABLE EDDINGTON TENSOR for the multi-D
//! implicit M1 solve from a 1-D FORMAL SOLUTION PER RADIAL COLUMN (stage S5 of
//! docs/dev/rad_m1_curvilinear_design.md, option D; tests_m1/runs_5d_vetcol/README.md).
//!
//! WHAT IS SOLVED.  Once per hydro step (every vet_col_every steps), for every column
//! i = is..ie at fixed (m, k, j), the column is taken as LATERALLY HOMOGENEOUS and
//!     mu dI/dr + (1 - mu^2)/r dI/dmu = chi (S - I)      (sp: spherical)
//!     mu dI/dx = chi (S - I)                            (Cartesian: plane-parallel)
//! is solved with the column's own extinction chi = rho (kappa_F + kappa_s) (M1_IW_KT)
//! and the grey source of closure = vet_sc (rad_m1_vet.cpp header):
//!     S = eps_th a T^4 + (1 - eps_th) E^n,  eps_th = min(rho kappa_P / chi, 1),
//!     (vet_scatter = true, opt-in, read only when named: S = VetScatterMix, rad_m1.hpp,
//!     the Jiang 2021 eq. 6 source of Ma+2026 or emission by kappa_R - kappa_e only)
//! in E units (eps = 4 pi I / c, E = J).  The Eddington factor f_K = K/J of the solution
//! is handed to the solve as the FIXED uniaxial tensor
//!     D = diag(f_K, (1-f_K)/2, (1-f_K)/2)   about n = r_hat (x1 on Cartesian),
//! or about the M1 cell flux direction (vet_col_axis = flux), through the tau closure's
//! tau_ten (tau_closure is set too; only this build differs).  f_K is clamped to
//! [vet_col_fk_min, 1]: 1/3 (the old clamp, from vet_sc's eigen/flux projection, where
//! chi >= 1/3 by construction) or 0 (sp default since m1-sp-order2b).  About r_hat f_K
//! < 1/3 is physical (limb-brightened intensity: behind a radiating front, inside an
//! emitting shell) and D stays realizable for 0 <= f_K <= 1; the 1/3 clamp put a kink
//! into f_K where K/J crosses 1/3 (the vet_col pulse was order 0.8, runs_5q_sporder2b).
//! vet_col_axis = flux keeps 1/3.
//!
//! RAYS (sp): IMPACT-PARAMETER (p-ray) geometry.  A straight ray keeps its impact
//! parameter p; at a shell r >= p it has mu = sqrt(r^2 - p^2)/r.  Shells are the cell
//! centres x1v; the inner boundary sphere is the inner face r_in, the top the outer face.
//! The ray set, ordered by p:
//!   * vet_col_ncore + 1 CORE rays, p = r_in sqrt(1 - mu_f^2), mu_f = 1 - q/ncore at r_in
//!     (uniform in mu at the inner boundary);
//!   * per shell l, vet_col_nsub - 1 rays whose TANGENT POINT lies inside the interval
//!     (r_{l-1}, r_l) (r_in for l = 0), then the ray tangent AT shell l (mu = 0 there).
//! A shell sees every ray with p <= r_l (a prefix of the list), from mu = 1 (p = 0) to
//! mu = 0 (its own tangent ray): its angular grid follows the beam, so a free-streaming
//! cone of half-width r_in/r is resolved by the core rays at every radius.
//! SWEEPS: (1) incoming rays top down, from vacuum at the outer face; (2) outgoing rays
//! bottom up: a core ray starts at r_in with the DIFFUSION intensity E + 3 F_r mu_f/c of
//! the bottom cell's M1 state (as vet_sc), a ray tangent at shell l starts there with
//! its own incoming intensity (mu = 0: in = out), and a ray whose tangent point lies
//! inside the interval below shell l goes from its incoming intensity at shell l down to
//! the tangent point and back (the MIRRORED incoming intensity), with chi and S at the
//! tangent point interpolated linearly in r.  Along every segment chi is the mean of
//! the two ends and S is linear in tau: the positive first-order short-characteristics
//! weights of vet_sc.  The end segments (inner face to the first shell, last shell to
//! the outer face) use the end cell's chi and S at the face linearly extrapolated from
//! the two end cells (clipped at 0), as vet_sc.  There is no mu interpolation at all:
//! every intensity is carried along its own straight ray (the ray-curvature term is
//! exact).
//! GEOMETRY IN RELATIVE RADII: the path lengths are formed from x = r/r_top (the column
//! top, as the two-stream spherical code measures its areas relative to the column top,
//! the lesson of the He mode-3 fix) with z = sqrt((x-p)(x+p)) and the segment length
//! (x_b^2 - x_a^2)/(z_a + z_b), free of the cancellation of z_b - z_a near a tangent
//! point, and are multiplied by r_top once.
//! QUADRATURE at a shell: the trapezoid rule in mu over its rays, corrected per shell by
//! a factor (a + b mu + c mu^2) so that 1, mu, mu^2 are integrated EXACTLY on each
//! hemisphere: f_K = 1/3 exactly for an isotropic or a diffusion (a + b mu) intensity,
//! i.e. the thick limit is exact whatever the angular resolution.
//!
//! RAYS (Cartesian): vet_col_nmu Gauss-Legendre nodes per hemisphere, path dx1/mu per
//! cell, the same boundary segments: the 1-D plane-parallel form of vet_sc's sweep (on a
//! laterally uniform state both give the same K/J up to the order of the sums).
//!
//! WHEN.  ImplicitSolve calls VetColBuild where it calls vet_sc's formal solution: at the
//! start-of-step state (E^n, T^n, the M1 F of the bottom cell), before the predictor.
//! Under time_scheme = hesdirk2 both stage solves sit at t^{n+1}: time2_vet_col
//! (rad_m1_time2.cpp) = predict (sp default) builds once per step at U^n + dt K1, second
//! order in time; lag keeps D(U^n) (first order).
//!
//! REFLECTING OUTER x1 (vet_col_reflect_top, sp default with implicit_bc_x1max =
//! reflect, m1-sp-order2b).  Specular reflection at the top sphere keeps p, so the
//! incoming intensity of each ray at the top face is I_in = b + a I_in: b the outgoing
//! intensity at the top face of a sweep from a vacuum top, a the ray's round-trip
//! transmission (the product of its segment factors exp(-dtau); 0 for a core ray that
//! starts from the diffusion intensity).  A first sweep gives b and a, a second one runs
//! from I_in = b/(1 - min(a, vet_col_reflect_amax)).  Exact for the steady transfer
//! problem; in a nearly transparent closed cavity (a -> 1) the intensity is not set by
//! the instantaneous source at all (it remembers its history), which no formal
//! solution can give: the cap (default 0.99) bounds I_in there.  Twice the build cost.
//!
//! OUTER BOUNDARY (vet_col_surface_q = true, runs_5e).  The outgoing intensities of the
//! top shell are carried over the last half segment to the top face (the incoming
//! sweep's first segment backwards), H(face) = 1/2 sum wf mu_f I_f with the face
//! quadrature (sp: every ray at mu_f = sqrt(1 - p^2) plus a node mu = 0 of zero
//! intensity, the same moment-corrected trapezoid; Cartesian: the Gauss nodes), and the
//! outer Marshak face flux F = c q E(top cell) takes, per column,
//!     q = H(face)/J(top cell)      clamped to [vet_col_surface_qmin, _qmax]
//! (J at the top CELL, not at the face: with E(top cell) = J the face flux is then the
//! formal solution's own H).  Lagged like the tensor (same build); the BC branches of
//! rad_m1_implicit.cpp shadow marshak_q with it only at the outer x1 face.
//!
//! SECOND ORDER (vet_col_order2 = true, m1-sp-order2, tests_m1/runs_5o_sporder2; default
//! false, bitwise).  The sp build above is order 1.1-1.4 in f_K: (a) the trapezoid in mu
//! over p-ray nodes that crowd towards mu = 0 like sqrt(dr), (b) S and chi linear in the
//! PATH along a segment although r is quadratic in the path near a tangent point, (c) the
//! core rays started from the bottom CELL's E and F.  With the switch: (a) the base
//! weights come from the piecewise-quadratic interpolant in mu (VcolQuad2, never across
//! the core edge), then the same 1, mu, mu^2 correction; (b) S and chi are linear in r:
//! the mean gb of g(u) = (r(u) - r_a)/(r_b - r_a) per ray-segment (host tables vcol_gb,
//! vcol_gb0) sets the chi integral and the quadratic-in-u source weights of VcolW2; (c)
//! the core rays start from E extrapolated to the inner face and the stored face flux
//! f0x1, or from the mirrored incoming intensity at a reflecting inner x1.  f_K is then
//! second order (L1 2.1 on the T-S6 (ii) state).  Only the weights and the core start
//! change: the kernels' layout and cost are the same (+1 ms per build on the He wedge).
//!
//! COST AND LAYOUT.  vet_col_team = true (default): VetColBuildTeam, one team per column
//! with the rays over the threads and the moments summed per shell by one thread in the
//! same order, bitwise the kernel below on CPU.  vet_col_team = false: the original
//! kernel, one thread per column, the rays looped inside: per column and sweep
//! (ncore + 1 + nsub nx1) nx1 ray-segments on sp, nmu nx1 on Cartesian; the running
//! intensity of every ray is kept in vcol_buf (ray, column), coalesced over the columns.
//! The moments are summed in a fixed order: the tensor is deterministic and IDENTICAL on
//! every column of a spherically symmetric state.  Column-local: no communication.
//!
//! LIMITS (fatal at start-up): transport = implicit on a multi-D mesh, ONE MeshBlock
//! along x1 (every block then spans the same radial grid: the ray tables are shared), no
//! SMR/AMR, no periodic x1, no cubed sphere.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_vetseg.hpp"
#include "rad_m1/rad_m1_vetsrc.hpp"

namespace radm1 {

namespace {
void VcolFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/closure = vet_col: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

// Gauss-Legendre nodes and weights on [0, 1]
void VcolGauss01(const int n, std::vector<double> &x, std::vector<double> &w) {
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  for (int a = 0; a < n; ++a) {
    double z = std::cos(M_PI*(a + 0.75)/(n + 0.5));
    double pp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int q = 1; q <= n; ++q) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*q - 1.0)*z*p2 - (q - 1.0)*p3)/q;
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;
      if (std::fabs(z - z1) < 1.0e-15) break;
    }
    x[a] = 0.5*(1.0 - z);
    w[a] = 1.0/((1.0 - z*z)*pp*pp);
  }
}

// the trapezoid rule in mu over nodes mu[0] > mu[1] > ... (from 1 down to 0), times
// (a + b mu + c mu^2) so that 1, mu, mu^2 are exact on the hemisphere; returns false
// (and the plain trapezoid weights) when the correction fails or a weight turns negative
bool VcolQuad(const std::vector<double> &mu, std::vector<double> &wout) {
  const int kr = static_cast<int>(mu.size()) - 1;
  std::vector<double> w(kr + 1);
  for (int r = 0; r <= kr; ++r) {
    double mprev = (r == 0) ? mu[0] : mu[r-1];
    double mnext = (r == kr) ? mu[kr] : mu[r+1];
    w[r] = 0.5*(mprev - mnext);
  }
  double mm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int r = 0; r <= kr; ++r) {
    double mq = 1.0;
    for (int q = 0; q < 5; ++q) {mm[q] += w[r]*mq; mq *= mu[r];}
  }
  // solve [mm0 mm1 mm2; mm1 mm2 mm3; mm2 mm3 mm4] (a b c) = (1, 1/2, 1/3)
  double A[3][4] = {{mm[0], mm[1], mm[2], 1.0}, {mm[1], mm[2], mm[3], 0.5},
                    {mm[2], mm[3], mm[4], 1.0/3.0}
                   };
  for (int c = 0; c < 3; ++c) {
    int piv = c;
    for (int rr = c+1; rr < 3; ++rr) {
      if (std::fabs(A[rr][c]) > std::fabs(A[piv][c])) piv = rr;
    }
    for (int q = 0; q < 4; ++q) {std::swap(A[c][q], A[piv][q]);}
    for (int rr = 0; rr < 3; ++rr) {
      if (rr == c) continue;
      double f = A[rr][c]/A[c][c];
      for (int q = c; q < 4; ++q) {A[rr][q] -= f*A[c][q];}
    }
  }
  double ca = A[0][3]/A[0][0], cb = A[1][3]/A[1][1], cc = A[2][3]/A[2][2];
  bool ok = std::isfinite(ca) && std::isfinite(cb) && std::isfinite(cc);
  for (int r = 0; r <= kr && ok; ++r) {
    double m = mu[r];
    if (w[r]*(ca + cb*m + cc*m*m) < 0.0) {ok = false;}
  }
  wout.assign(kr + 1, 0.0);
  for (int r = 0; r <= kr; ++r) {
    double m = mu[r];
    wout[r] = ok ? w[r]*(ca + cb*m + cc*m*m) : w[r];
  }
  return ok;
}

// vet_col_quad = quadratic (m1-sp-order2, tests_m1/runs_5o_sporder2): the base weights
// of VcolQuad from the piecewise-QUADRATIC interpolant in mu instead of the trapezoid.
// The p-ray nodes crowd towards mu = 0 like sqrt(dr) (the rays tangent at the shells
// below), so the trapezoid leaves an O(dr^1.5) error in J, H, K there; per interval the
// quadratic through three neighbouring nodes (the mean of the two overlapping triples
// where both exist) leaves O(h^4) per interval, O(dr^2) summed.  The rule never
// straddles node `split` (the last core ray, mu_c: the core edge, where the outgoing
// intensity may jump), a panel of one interval is the trapezoid.  Then the moment
// correction of VcolQuad (1, mu, mu^2 exact).  Returns false (and VcolQuad's weights)
// if a weight turns negative.
double VcolLagInt(const double xp, const double xq, const double xr, const double lo,
                  const double hi) {
  // int_lo^hi (x - xq)(x - xr) dx / ((xp - xq)(xp - xr))
  auto prim = [&](double x) {
    return x*x*x/3.0 - 0.5*(xq + xr)*x*x + xq*xr*x;
  };
  return (prim(hi) - prim(lo))/((xp - xq)*(xp - xr));
}

bool VcolQuad2(const std::vector<double> &mu, const int split,
               std::vector<double> &wout) {
  const int kr = static_cast<int>(mu.size()) - 1;
  std::vector<double> w(kr + 1, 0.0);
  const int sp = std::min(std::max(split, 0), kr);
  const int pa[2] = {0, sp}, pb[2] = {sp, kr};
  for (int q = 0; q < 2; ++q) {
    const int a = pa[q], b = pb[q];
    if (b - a == 1) {
      w[a] += 0.5*(mu[a] - mu[a+1]);
      w[a+1] += 0.5*(mu[a] - mu[a+1]);
      continue;
    }
    for (int k = a; k < b; ++k) {
      const double lo = mu[k+1], hi = mu[k];
      int nt = 0;
      double wk[4] = {0.0, 0.0, 0.0, 0.0};   // weights of nodes k-1, k, k+1, k+2
      if (k - 1 >= a) {
        wk[0] += VcolLagInt(mu[k-1], mu[k], mu[k+1], lo, hi);
        wk[1] += VcolLagInt(mu[k], mu[k-1], mu[k+1], lo, hi);
        wk[2] += VcolLagInt(mu[k+1], mu[k-1], mu[k], lo, hi);
        ++nt;
      }
      if (k + 2 <= b) {
        wk[1] += VcolLagInt(mu[k], mu[k+1], mu[k+2], lo, hi);
        wk[2] += VcolLagInt(mu[k+1], mu[k], mu[k+2], lo, hi);
        wk[3] += VcolLagInt(mu[k+2], mu[k], mu[k+1], lo, hi);
        ++nt;
      }
      for (int t = 0; t < 4; ++t) {
        const int r = k - 1 + t;
        if (r >= 0 && r <= kr) {w[r] += wk[t]/nt;}
      }
    }
  }
  double mm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int r = 0; r <= kr; ++r) {
    double mq = 1.0;
    for (int q = 0; q < 5; ++q) {mm[q] += w[r]*mq; mq *= mu[r];}
  }
  double A[3][4] = {{mm[0], mm[1], mm[2], 1.0}, {mm[1], mm[2], mm[3], 0.5},
                    {mm[2], mm[3], mm[4], 1.0/3.0}
                   };
  for (int c = 0; c < 3; ++c) {
    int piv = c;
    for (int rr = c+1; rr < 3; ++rr) {
      if (std::fabs(A[rr][c]) > std::fabs(A[piv][c])) piv = rr;
    }
    for (int q = 0; q < 4; ++q) {std::swap(A[c][q], A[piv][q]);}
    for (int rr = 0; rr < 3; ++rr) {
      if (rr == c) continue;
      double f = A[rr][c]/A[c][c];
      for (int q = c; q < 4; ++q) {A[rr][q] -= f*A[c][q];}
    }
  }
  double ca = A[0][3]/A[0][0], cb = A[1][3]/A[1][1], cc = A[2][3]/A[2][2];
  bool ok = std::isfinite(ca) && std::isfinite(cb) && std::isfinite(cc);
  for (int r = 0; r <= kr && ok; ++r) {
    double m = mu[r];
    if (w[r]*(ca + cb*m + cc*m*m) < 0.0) {ok = false;}
  }
  if (!ok) {
    VcolQuad(mu, wout);
    return false;
  }
  wout.assign(kr + 1, 0.0);
  for (int r = 0; r <= kr; ++r) {
    double m = mu[r];
    wout[r] = w[r]*(ca + cb*m + cc*m*m);
  }
  return true;
}

// ray types in vcol_ray(r, 1)
constexpr int VC_CORE = 0;   // starts at the inner boundary sphere (face r_in)
constexpr int VC_TAN = 1;    // tangent exactly at its first shell L
constexpr int VC_SUB = 2;    // tangent point inside the interval below shell L

// vet_col_order2: the short-characteristics weights with S(u) = S_up + (S_dn - S_up)
// (u + a2 u (u - 1)) along the segment (u = the fraction of the optical path from the
// upstream end), a2 in [-1, 1]: a2 = 0 is VcolW (S linear in tau).  With
// W_k = int_0^1 u^k x e^{-x(1-u)} du: w0 = (1 - a2) W1 + a2 W2 (downstream end),
// wu = (1 - e^-x) - w0 (upstream end); both >= 0 (the profile is monotone in u).
KOKKOS_INLINE_FUNCTION
void VcolW2(const Real dtau, const Real a2, Real &ex, Real &w0, Real &wu) {
  ex = exp(-dtau);
  Real w1, w2;
  if (dtau < 1.0e-3) {
    w1 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
  } else {
    w1 = 1.0 - (1.0 - ex)/dtau;
  }
  if (dtau < 1.0e-2) {
    w2 = dtau*(1.0/3.0 - dtau*(1.0/12.0 - dtau*(1.0/60.0 - dtau*(1.0/360.0
                                                                 - dtau/2520.0))));
  } else {
    w2 = 1.0 - 2.0*w1/dtau;
  }
  w0 = (1.0 - a2)*w1 + a2*w2;
  wu = (1.0 - ex) - w0;
}

// the team build sums each shell's moments in VC_NG groups of consecutive rays
constexpr int VC_NG = 8;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetColInit
//! \brief checks, the ray tables (host, shared by every block), buffers

void RadiationM1::VetColInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  if (!trans_on) {VcolFatal("needs <rad_m1>/transport = implicit on a MULTI-D mesh");}
  if (pm->multilevel) {VcolFatal("SMR/AMR is not supported");}
  // (the clean closure from vet_gd sets up but never builds the columns)
  if (mindcs.nx1 != indcs.nx1 && !(clean && vgd_on)) {
    VcolFatal("needs ONE MeshBlock along x1 (meshblock/nx1 = mesh/nx1)");
  }
  if (ibc_x1min == M1_IBC_PERIODIC) {VcolFatal("a periodic x1 boundary has no top");}
  // STAGE CS2 (m1-cs-implicit): the cubed sphere with the RADIAL axis.  x1 is r on every
  // panel and one block spans the whole column, so the column is the radial line of a
  // panel cell and the p-ray tables of block 0 serve every block, as on sp.  The tensor
  // is uniaxial about r_hat; its tangential part (1-f_K)/2 is an isotropic pressure,
  // which the cs transverse rows read as the scalar D_tt (M1DDiag, CS1).  The flux axis
  // needs the cell F direction on the skewed panel basis: refused.
  if (pm->use_cubed_sphere && vcol_axis_flux) {
    VcolFatal("vet_col_axis = flux is not supported on the cubed sphere");
  }
  if (indcs.nx1 < 2) {VcolFatal("needs at least 2 cells along x1");}
  if (vcol_nc < 1 || vcol_np < 1 || vcol_nmu < 1 || vcol_every < 1) {
    VcolFatal("vet_col_ncore, vet_col_nsub, vet_col_nmu, vet_col_every must be >= 1");
  }
  vcol_sph = sph_geom;
  // vet_col_reflect_top: only on sp with a reflecting outer x1 of the M1 solve
  vcol_rtop = vcol_rtop && vcol_sph && (ibc_x1max == M1_IBC_REFLECT);
  vcol_sqf = vcol_sqf && vcol_sph && vcol_sq && impl_mface_lin;
  if (vcol_rtop && !(vcol_amax >= 0.0 && vcol_amax < 1.0)) {
    VcolFatal("vet_col_reflect_amax must be in [0, 1)");
  }
  const int n1 = indcs.nx1, is = indcs.is, ie = indcs.ie;

  // the radial grid of block 0: with one block along x1 and no SMR every block has it
  std::vector<double> rc(n1);
  double rin, rtop, dx1 = 0.0;
  if (vcol_sph) {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    for (int l = 0; l < n1; ++l) {rc[l] = x1v_h(0, is + l);}
    rin = x1f_h(0, is);
    rtop = x1f_h(0, ie + 1);
    // vet_col_surface_face: the face extrapolation of J (M1SphMarshakFaceE's weights)
    M1SphMarshakCoef(rc[n1-1], rc[n1-2], rtop, vcol_jca, vcol_jcb);
    vcol_jcap = 2.0*SQR(rc[n1-1]/rtop);
  } else {
    auto &sz = pmy_pack->pmb->mb_size.h_view(0);
    dx1 = sz.dx1;
    for (int l = 0; l < n1; ++l) {rc[l] = sz.x1min + (l + 0.5)*dx1;}
    rin = sz.x1min;
    rtop = sz.x1max;
  }
  vcol_rc = rc;

  // ---- the ray set
  std::vector<double> p, at, mf;
  std::vector<int> lr, ty;
  std::vector<int> klast(n1);
  const double xin = rin/rtop;
  std::vector<double> xr(n1 + 1);
  for (int l = 0; l < n1; ++l) {xr[l] = rc[l]/rtop;}
  xr[n1] = 1.0;
  if (vcol_sph) {
    for (int q = 0; q <= vcol_nc; ++q) {
      double muf = 1.0 - static_cast<double>(q)/vcol_nc;
      p.push_back(xin*std::sqrt((1.0 - muf)*(1.0 + muf)));
      lr.push_back(0);
      ty.push_back(VC_CORE);
      at.push_back(0.0);
      mf.push_back(muf);
    }
    for (int l = 0; l < n1; ++l) {
      double xlo = (l == 0) ? xin : xr[l-1];
      for (int s = 1; s < vcol_np; ++s) {
        double a = static_cast<double>(s)/vcol_np;
        p.push_back(xlo + (xr[l] - xlo)*a);
        lr.push_back(l);
        ty.push_back(VC_SUB);
        at.push_back(a);
        mf.push_back(0.0);
      }
      p.push_back(xr[l]);
      lr.push_back(l);
      ty.push_back(VC_TAN);
      at.push_back(1.0);
      mf.push_back(0.0);
      klast[l] = static_cast<int>(p.size()) - 1;
    }
  } else {
    std::vector<double> xg, wg;
    VcolGauss01(vcol_nmu, xg, wg);
    for (int q = 0; q < vcol_nmu; ++q) {
      p.push_back(0.0);
      lr.push_back(0);
      ty.push_back(VC_CORE);
      at.push_back(0.0);
      mf.push_back(xg[q]);
    }
    for (int l = 0; l < n1; ++l) {klast[l] = vcol_nmu - 1;}
  }
  const int nray = static_cast<int>(p.size());
  vcol_nray = nray;

  Kokkos::realloc(vcol_seg, nray, n1);
  Kokkos::realloc(vcol_mu, nray, n1);
  Kokkos::realloc(vcol_w, nray, n1);
  Kokkos::realloc(vcol_ray, nray, 4);
  const bool go2 = vcol_sph && vcol_o2;
  Kokkos::realloc(vcol_gb, go2 ? nray : 1, go2 ? n1 : 1);
  Kokkos::realloc(vcol_gb0, go2 ? nray : 1);
  auto gb_h = Kokkos::create_mirror_view(vcol_gb);
  auto gb0_h = Kokkos::create_mirror_view(vcol_gb0);
  for (int r = 0; r < static_cast<int>(gb_h.extent(0)); ++r) {
    gb0_h(r) = 0.5;
    for (int l = 0; l < static_cast<int>(gb_h.extent(1)); ++l) {gb_h(r, l) = 0.5;}
  }
  Kokkos::realloc(vcol_klast, n1);
  auto seg_h = Kokkos::create_mirror_view(vcol_seg);
  auto mu_h = Kokkos::create_mirror_view(vcol_mu);
  auto w_h = Kokkos::create_mirror_view(vcol_w);
  auto ray_h = Kokkos::create_mirror_view(vcol_ray);
  auto kl_h = Kokkos::create_mirror_view(vcol_klast);
  Kokkos::realloc(vcol_muf, nray);
  Kokkos::realloc(vcol_wf, nray);
  auto muf_h = Kokkos::create_mirror_view(vcol_muf);
  auto wf_h = Kokkos::create_mirror_view(vcol_wf);
  for (int r = 0; r < nray; ++r) {
    for (int l = 0; l < n1; ++l) {seg_h(r, l) = mu_h(r, l) = w_h(r, l) = 0.0;}
  }
  int nneg = 0;
  if (vcol_sph) {
    for (int r = 0; r < nray; ++r) {
      const double pr = p[r];
      auto zz = [&](double x) {return std::sqrt(std::fmax((x - pr)*(x + pr), 0.0));};
      for (int l = lr[r]; l < n1; ++l) {
        double za = (ty[r] == VC_TAN && l == lr[r]) ? 0.0 : zz(xr[l]);
        double zb = zz(xr[l+1]);
        seg_h(r, l) = rtop*((xr[l+1] - xr[l])*(xr[l+1] + xr[l])/(za + zb));
        mu_h(r, l) = za/xr[l];
      }
      double s0 = 0.0;
      if (ty[r] == VC_CORE) {
        double zf = zz(xin), z0 = zz(xr[0]);
        s0 = rtop*((xr[0] - xin)*(xr[0] + xin)/(zf + z0));
      } else if (ty[r] == VC_SUB) {
        s0 = rtop*zz(xr[lr[r]]);
      }
      ray_h(r, 0) = lr[r];
      ray_h(r, 1) = ty[r];
      ray_h(r, 2) = s0;
      ray_h(r, 3) = (ty[r] == VC_CORE) ? mf[r] : at[r];
    }
    if (vcol_o2) {
      // vet_col_order2: gb per segment, by 8-point Gauss in u.  r(u) = sqrt(p^2 +
      // (z_a + u (z_b - z_a))^2) in relative radii, z = sqrt(x^2 - p^2)
      std::vector<double> xg, wg;
      VcolGauss01(8, xg, wg);
      auto gbar = [&](double pr, double xa, double xb) {
        const double za = std::sqrt(std::fmax((xa - pr)*(xa + pr), 0.0));
        const double zb = std::sqrt(std::fmax((xb - pr)*(xb + pr), 0.0));
        if (!(xb > xa) || !(zb > za)) {return 0.5;}
        double g = 0.0;
        for (int q = 0; q < 8; ++q) {
          const double z = za + xg[q]*(zb - za);
          g += wg[q]*(std::sqrt(pr*pr + z*z) - xa)/(xb - xa);
        }
        return g;
      };
      for (int r = 0; r < nray; ++r) {
        for (int l = lr[r]; l < n1; ++l) {gb_h(r, l) = gbar(p[r], xr[l], xr[l+1]);}
        if (ty[r] == VC_CORE) {
          gb0_h(r) = gbar(p[r], xin, xr[0]);
        } else if (ty[r] == VC_SUB) {
          gb0_h(r) = gbar(p[r], p[r], xr[lr[r]]);
        }
      }
    }
    // quadrature per shell: trapezoid in mu over rays 0..klast (mu from 1 down to 0),
    // times (a + b mu + c mu^2) so that 1, mu, mu^2 are exact on the hemisphere
    for (int l = 0; l < n1; ++l) {
      const int kr = klast[l];
      std::vector<double> mus(kr + 1), wq;
      for (int r = 0; r <= kr; ++r) {mus[r] = mu_h(r, l);}
      const bool okq = vcol_o2 ? VcolQuad2(mus, vcol_nc, wq)
                                        : VcolQuad(mus, wq);
      if (!okq) {++nneg;}
      for (int r = 0; r <= kr; ++r) {w_h(r, l) = wq[r];}
    }
    // the top FACE (vet_col_surface_q): every ray, mu_f = sqrt(1 - p^2) (x = 1), and a
    // node mu = 0 whose outgoing intensity is 0 (no path inside the domain)
    std::vector<double> muf(nray + 1), wf;
    for (int r = 0; r < nray; ++r) {
      muf[r] = std::sqrt(std::fmax((1.0 - p[r])*(1.0 + p[r]), 0.0));
    }
    muf[nray] = 0.0;
    const bool okf = vcol_o2 ? VcolQuad2(muf, vcol_nc, wf)
                                      : VcolQuad(muf, wf);
    if (!okf) {++nneg;}
    for (int r = 0; r < nray; ++r) {muf_h(r) = muf[r]; wf_h(r) = wf[r];}
  } else {
    std::vector<double> xg, wg;
    VcolGauss01(vcol_nmu, xg, wg);
    for (int r = 0; r < nray; ++r) {
      for (int l = 0; l < n1; ++l) {
        seg_h(r, l) = ((l == n1 - 1) ? 0.5*dx1 : dx1)/xg[r];
        mu_h(r, l) = xg[r];
        w_h(r, l) = wg[r];
      }
      muf_h(r) = xg[r];
      wf_h(r) = wg[r];
      ray_h(r, 0) = 0;
      ray_h(r, 1) = VC_CORE;
      ray_h(r, 2) = 0.5*dx1/xg[r];
      ray_h(r, 3) = xg[r];
    }
  }
  for (int l = 0; l < n1; ++l) {kl_h(l) = klast[l];}
  Kokkos::deep_copy(vcol_seg, seg_h);
  Kokkos::deep_copy(vcol_mu, mu_h);
  Kokkos::deep_copy(vcol_w, w_h);
  Kokkos::deep_copy(vcol_ray, ray_h);
  Kokkos::deep_copy(vcol_gb, gb_h);
  Kokkos::deep_copy(vcol_gb0, gb0_h);
  Kokkos::deep_copy(vcol_klast, kl_h);
  // the team sweep's copy of seg and gb, (shell, ray, 2): the same values
  Kokkos::realloc(vcol_sgt, n1, nray, 2);
  {
    auto sgt_h = Kokkos::create_mirror_view(vcol_sgt);
    const bool hgb = (gb_h.extent_int(0) == nray);
    for (int l = 0; l < n1; ++l) {
      for (int r = 0; r < nray; ++r) {
        sgt_h(l, r, 0) = seg_h(r, l);
        sgt_h(l, r, 1) = hgb ? gb_h(r, l) : 0.5;
      }
    }
    Kokkos::deep_copy(vcol_sgt, sgt_h);
  }
  Kokkos::deep_copy(vcol_muf, muf_h);
  Kokkos::deep_copy(vcol_wf, wf_h);

  const int nmb = pmy_pack->nmb_thispack;
  const int ncol = nmb*indcs.nx2*indcs.nx3;
  int c1 = indcs.nx1 + 2*(indcs.ng);
  int c2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int c3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  if (vcol_team) {
    // shells per chunk: at most 16 and the team scratch (5 n1 + 2 nray + lc nray reals)
    // within 24 kB (runs_5e GPU sweep: 32 kB halves the occupancy, 3.3 ms at 16 vs 7.7
    // ms at 32 on the He wedge grid)
    const int avail = 3072 - 5*n1 - (vcol_rtop ? 3 : 2)*nray;
    // (m1-fast5-sp: 8 shells per chunk with the grouped moment sums, runs_5s job
    // 11980678: 4.95 ms per build at 8, 5.2-5.9 at 6/7/10/12, 7.5 at 16 on the He wedge)
    vcol_lc = std::max(1, std::min(std::min(n1, 8), avail/nray));
    if (vcol_lcin > 0) {vcol_lc = std::min(n1, vcol_lcin);}
  } else {
    Kokkos::realloc(vcol_buf, nray, ncol);
    Kokkos::deep_copy(vcol_buf, 0.0);
    if (vcol_rtop) {
      Kokkos::realloc(vcol_abuf, nray, ncol);
      Kokkos::deep_copy(vcol_abuf, 0.0);
    }
  }
  if (vcol_sq) {
    Kokkos::realloc(vcol_q, nmb, c3, c2);
    Kokkos::deep_copy(vcol_q, marshak_q);
    if (ibc_x1max != M1_IBC_MARSHAK && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> vet_col_surface_q: the outer x1 boundary is not marshak; "
                << "the q is built but unused" << std::endl;
    }
  }
  // vet_col_lat (rad_m1_vetlat.cpp): M1_TT_NLAT more slots, the lateral correction.
  // STAGE CS2: on the cubed sphere tau_ten carries the same 6 slots, M1_TT_LAT0..+5, the
  // LATERAL CORRECTION interface of the closure tensor (the planned short-
  // characteristics correction of vet_col, rt_design_1003/SP_SC_CLOSURE.md): D_ab in
  // the cs basis (r, xi, eta), covariant on the panel pair like F, order rr, r xi,
  // r eta, xi xi, xi eta, eta eta.  Zero, and read by nothing yet (CS2_RESULTS.md,
  // "interface"); the uniaxial (chi, n) of slots 0..3 is what runs.  Elsewhere 4 slots.
  Kokkos::realloc(tau_ten, nmb, (vlat_on || pm->use_cubed_sphere) ?
                  (M1_TT_LAT0 + M1_TT_NLAT) : 4, c3, c2, c1);
  Kokkos::deep_copy(tau_ten, 0.0);
  if (!vcol_dump.empty()) {
    Kokkos::realloc(vcol_mom, nmb, 5, c3, c2, c1);
    Kokkos::deep_copy(vcol_mom, 0.0);
  }
  tau_ready = true;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> closure = vet_col: per-column 1-D formal solution, "
              << (vcol_sph ? "spherical p-rays (ncore "
                           : "plane-parallel Gauss rays (nmu ")
              << (vcol_sph ? vcol_nc : vcol_nmu);
    if (vcol_sph) {std::cout << ", nsub " << vcol_np;}
    std::cout << ", " << nray << " rays, " << n1 << " shells), axis "
              << (vcol_axis_flux ? "M1 flux" : "radial") << ", rebuilt every "
              << vcol_every << " step(s)";
    if (vcol_o2) {std::cout << "; vet_col_order2 (quadratic mu, face core)";}
    std::cout << "; vet_bc_diff = " << (vet_bc_face ? "face" : "cell");
    if (vcol_fkmin != 1.0/3.0) {
      std::cout << "; f_K >= " << vcol_fkmin;
    }
    if (vcol_rtop) {
      std::cout << "; reflecting top (mirrored, a <= " << vcol_amax << ")";
    }
    if (nneg > 0) {std::cout << "; " << nneg << " shell(s) kept trapezoid weights";}
    if (vcol_team) {
      std::cout << "; team build, " << vcol_lc << " shells per chunk, team size "
                << ((vcol_ts > 0) ? std::to_string(vcol_ts) : std::string("AUTO"));
    }
    if (vcol_sq) {
      std::cout << "; outer Marshak q = H(face)/J(top cell) in [" << vcol_qmin << ", "
                << vcol_qmax << "]";
    }
    std::cout << std::endl;
  }
  if (vlat_on) {VetLatInit();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetColReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationM1::VetColReport() {
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_m1> closure = vet_col: " << vcol_ncall << " tensor builds ("
            << vcol_nskip << " skipped), " << vcol_time << " s (rank 0, fenced), "
            << ((vcol_ncall > 0.0) ? (1.0e3*vcol_time/vcol_ncall) : 0.0)
            << " ms per build" << std::endl;
  if (vlat_on) {
    std::cout << "<rad_m1> vet_col_lat: " << vlat_ncall << " sweeps (incl. the "
              << "first build's " << vlat_iinit << " inflow iterations), " << vlat_time
              << " s (rank 0, fenced), "
              << ((vlat_ncall > 0.0) ? (1.0e3*vlat_time/vlat_ncall) : 0.0)
              << " ms per sweep; lateral reads clamped to the ghost band (rank 0, "
              << "all sweeps): " << vlat_nclamp << std::endl;
    if (vgd_on) {
      std::cout << "<rad_m1> vet_gd: rank 0 fenced seconds: source " << vgd_tsrc
                << ", sweeps " << vgd_tswp << " (of which per-shell halo+walls "
                << vgd_thalo << ", of which MPI wait " << vgd_tmpi
                << "), tables (rotation) " << vgd_ttab << ", uniform twin " << vgd_ttwin
                << ", moments " << vgd_tmom << "; in-solve rebuilds " << vgd_nrb << " ("
                << vgd_trb << " s)" << std::endl;
      if (vra_on) {
        std::cout << "<rad_m1> vet_rotate_avg = " << vra_k << ": snapshot sweeps "
                  << vra_nsweep << " (rotations " << vra_nrot << ", inits " << vra_ninit
                  << "), " << vra_time << " s incl. offsets (rank 0, fenced)"
                  << std::endl;
      }
      if (vgd_symref) {
        std::cout << "<rad_m1> vet_gd_symref: host start " << vgd_tsr << " s, join "
                  << vgd_tsrj << " s (rank 0)" << std::endl;
      }
      std::cout << "<rad_m1> vet_gd halo (rank 0): exchanges " << vgd_nexch
                << ", bytes sent "
                << vgd_nbyte << " (" << ((vgd_nexch > 0.0) ? vgd_nbyte/vgd_nexch : 0.0)
                << " per exchange), post " << vgd_tpost << " s, partners "
                << vgd_prk.size() << "; band-clamped lateral reads (ALL ranks, all "
                << "sweeps) "
                << vgd_nclamp_all << std::endl;
      if (vgd_csg || vgd_gbc) {   // cs-s3b-1017 (share-1104: + the sp global band)
        std::cout << "<rad_m1> vet_gd " << (vgd_csg ? "cs" : "global") << " band (rank 0): "
                  << ((vgd_cscomp || vgd_gbc) ? "compact" : "dense") << ", intensity values sent "
                  << vgd_csvsent << " of dense " << vgd_csvdense << " ("
                  << ((vgd_csvdense > 0.0) ? vgd_csvsent/vgd_csvdense : 0.0)
                  << "); list builds " << vgd_csnbuild << " (" << vgd_cstbuild << " s)"
                  << std::endl;
      }
      if (vgd_csg && vgd_cscomp) {   // raycheap-1109: where the list builds go
        std::cout << "<rad_m1> raycheap-1109 cs list builds (rank 0 s): mask kernel "
                  << vgd_cstb_k << ", download+check " << vgd_cstb_d << ", mark exchange "
                  << vgd_cstb_m << ", host lists " << vgd_cstb_h << ", upload "
                  << vgd_cstb_u << "; mask stride " << vgd_csper
                  << (vgd_cspad ? " (padded)" : "") << std::endl;
      }
      std::cout << "<rad_m1> vet_gd vet_path = " << ((vet_path != 0) ? "cells" : "none")
                << ": subdivided steps (ALL ranks, all sweeps) " << vgd_npath << std::endl;
      {
        // accel-1009 ragged band: accesses outside a shell's slab (must be 0)
        Real no = 0.0, nt = 0.0;
        if (vgd_i.oos.extent(0) > 0) {
          auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i.oos);
          no = h(0);
        }
        if (vgd_itw.oos.extent(0) > 0) {
          auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_itw.oos);
          nt = h(0);
        }
        std::cout << "<rad_m1> vet_gd ragged band (rank 0): out-of-slab accesses main "
                  << no << " twin " << nt << std::endl;
        if (vgd_ring) {
          std::cout << "<rad_m1> vet_gd_band_ring (rank 0): sweeps " << vgd_rg_nk
                    << ", dense band saves " << vgd_rg_ndense << ", save store "
                    << vgd_rg_store.extent(0)*sizeof(Real)/1.0e9 << " GB"
                    << (vgd_ring_poison ? " (POISON test)" : "") << std::endl;
        }
        if (vgd_xff) {
          double oos2 = 0.0;
          if (vgd_i_xb.oos.extent(0) > 0) {
            auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i_xb.oos);
            oos2 = h(0);
          }
          std::cout << "<rad_m1> vet_gd_rotate_xfade_fuse (rank 0): fused builds "
                    << vgd_xff_n << ", two-sweep fallbacks " << vgd_xff_nfb
                    << ", merged exchanges " << vgd_xff_nmsg << ", out-of-slab accesses "
                    << "(other set) " << oos2 << ", early-sent shells " << vgd_xfe_n
                    << std::endl;
        }
      }
      if (vgd_hc_mb > 0 || vgd_twfuse || vgd_twdet || vgd_twseq) {
        std::cout << "<rad_m1> accel-1009 (rank 0): vet_gd_twin_fuse=" << vgd_twfuse
                  << " vet_gd_twin_lowmem=" << vgd_twseq
                  << " vet_gd_twin_det=" << vgd_twdet << " halo mask cache: "
                  << vgd_hc_nmade << " entries made, " << vgd_hc_bytes/1048576.0
                  << " MB of " << vgd_hc_mb << std::endl;
      }
      if (vgd_hcomp > 0 && vgd_hl_mb > 0) {
        std::cout << "<rad_m1> vgdspeed-1009 (rank 0): vet_gd halo lists "
                  << vgd_hl_nmade << " made, " << vgd_hl_bytes/1048576.0 << " MB of "
                  << vgd_hl_cap/1048576.0 << " (key " << vgd_hl_mb
                  << ", GPU: <= free/2); exchanges by list " << vgd_hl_nuse << ", dense "
                  << vgd_hl_ndense << std::endl;
      }
      if (vgd_hx_on) {
        std::cout << "<rad_m1> vgdfuse-1010 (rank 0): vet_gd_halo_exact lists made "
                  << vgd_hx_nmade << ", fallbacks " << vgd_hx_nfall
                  << " (marked reads outside the masks' region " << vgd_hx_nout
                  << "), read ghosts with an unsent source " << vgd_hx_nun
                  << "; received values per list: exact/mask "
                  << ((vgd_hx_vmask > 0.0) ? vgd_hx_vexact/vgd_hx_vmask : 0.0)
                  << std::endl;
      }
      std::cout << "<rad_m1> vgdrot-1011 (rank 0): exact-list batches " << vgd_hxb_ngroup
                << " (group " << vgd_hxb_g << ", " << vgd_hxb_t << " s, list reallocs "
                << vgd_hxb_nrealloc << "); table prefetch " << (vgd_tab_pre ? 1 : 0)
                << " swapped in " << vgd_tab_nswap << " of the rotations" << std::endl;
      if (vgd_ovl_g > 0 || vgd_fuse_h > 1) {
        std::cout << "<rad_m1> vgdfuse-1010 (rank 0): vet_gd_overlap=" << vgd_ovl_g
                  << " split shells " << vgd_ovl_nsplit << " (remote items "
                  << ((vgd_ovl_nall > 0.0) ? vgd_ovl_nrem/vgd_ovl_nall : 0.0)
                  << " of all), unsplit " << vgd_ovl_nfull
                  << (vgd_ovl_nofit ? " (over vet_gd_overlap_mb)" : "")
                  << "; vet_gd_fuse_shells=" << vgd_fuse_h << " skipped exchanges "
                  << vgd_fuse_nskip << std::endl;
      }
      if (vgd_async) {
        std::cout << "<rad_m1> vet_gd_async (rank 0): " << vgd_nasync
                  << " overlapped sweeps (" << (vgd_ainl ? "inline" : "thread")
                  << "), wait at the join " << vgd_tjoin << " s" << std::endl;
      }
    }
  }
}

} // namespace radm1
