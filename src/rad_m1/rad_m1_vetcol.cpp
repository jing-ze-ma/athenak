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
//! in E units (eps = 4 pi I / c, E = J).  The Eddington factor f_K = K/J of the solution
//! is handed to the solve as the FIXED uniaxial tensor
//!     D = diag(f_K, (1-f_K)/2, (1-f_K)/2)   about n = r_hat (x1 on Cartesian),
//! or about the M1 cell flux direction (vet_col_axis = flux), through the tau closure's
//! tau_ten (tau_closure is set too; only this build differs).  f_K is clamped to
//! [1/3, 1].
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
//! SMR/AMR, no periodic x1, no cubed sphere, time_scheme = be (the tau closure's rule).

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
#include "rad_m1/rad_m1_implicit.hpp"

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

// the positive first-order short-characteristics weights (vet_sc)
KOKKOS_INLINE_FUNCTION
void VcolW(const Real dtau, Real &ex, Real &w0, Real &wu) {
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

// vet_col_order2: a2 of VcolW2 from the mean gb of g(u) over an OUTGOING segment (u from
// the inner end): the quadratic g = c1 u + c2 u^2 with c1 + c2 = 1 and mean gb has
// c2 = 3 - 6 gb, clamped to [0, 1] (monotone).  Outgoing: a2 = c2; incoming (u from the
// outer end, g -> 1 - g(1 - u)): a2 = -c2.
KOKKOS_INLINE_FUNCTION
Real VcolC2(const Real gb) {
  return fmin(fmax(3.0 - 6.0*gb, 0.0), 1.0);
}
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
  if (mindcs.nx1 != indcs.nx1) {
    VcolFatal("needs ONE MeshBlock along x1 (meshblock/nx1 = mesh/nx1)");
  }
  if (ibc_x1min == M1_IBC_PERIODIC) {VcolFatal("a periodic x1 boundary has no top");}
  if (pm->use_cubed_sphere) {VcolFatal("the cubed sphere is not supported (stage S5)");}
  if (indcs.nx1 < 2) {VcolFatal("needs at least 2 cells along x1");}
  if (vcol_nc < 1 || vcol_np < 1 || vcol_nmu < 1 || vcol_every < 1) {
    VcolFatal("vet_col_ncore, vet_col_nsub, vet_col_nmu, vet_col_every must be >= 1");
  }
  vcol_sph = sph_geom;
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
    const int avail = 3072 - 5*n1 - 2*nray;
    vcol_lc = std::max(1, std::min(std::min(n1, 16), avail/nray));
    if (vcol_lcin > 0) {vcol_lc = std::min(n1, vcol_lcin);}
  } else {
    Kokkos::realloc(vcol_buf, nray, ncol);
    Kokkos::deep_copy(vcol_buf, 0.0);
  }
  if (vcol_sq) {
    Kokkos::realloc(vcol_q, nmb, c3, c2);
    Kokkos::deep_copy(vcol_q, marshak_q);
    if (ibc_x1max != M1_IBC_MARSHAK && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> vet_col_surface_q: the outer x1 boundary is not marshak; "
                << "the q is built but unused" << std::endl;
    }
  }
  Kokkos::realloc(tau_ten, nmb, 4, c3, c2, c1);
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
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetColBuild
//! \brief the formal solution of every column -> tau_ten (chi, n1, n2, n3)

void RadiationM1::VetColBuild() {
  if (!tau_ready) {VetColInit();}
  const int cyc = pmy_pack->pmesh->ncycle;
  if (vcol_built && vcol_every > 1 && (cyc % vcol_every) != 0) {
    vcol_nskip += 1.0;
    return;
  }
  Kokkos::fence();
  Kokkos::Timer timer;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int n1 = indcs.nx1, nj = indcs.nx2, nk = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const bool thermal = (pmy_pack->phydro != nullptr) && coupling && !opac_zero;
  const bool axf = vcol_axis_flux;
  const bool dmp = !vcol_dump.empty() && (vcol_dump_every > 0) &&
                   ((static_cast<int>(vcol_ncall) % vcol_dump_every) == 0);
  const Real cl = c_light, ar = arad, efl = e_floor;
  auto iw_ = iw;
  auto opac_ = opac;
  auto tt_ = tau_ten;
  auto mo_ = vcol_mom;
  auto seg_ = vcol_seg;
  auto mu_ = vcol_mu;
  auto w_ = vcol_w;
  auto ray_ = vcol_ray;
  auto kl_ = vcol_klast;
  auto buf_ = vcol_buf;
  auto muf_ = vcol_muf;
  auto wf_ = vcol_wf;
  auto vq_ = vcol_q;
  const bool sq = vcol_sq;
  const Real qlo = vcol_qmin, qhi = vcol_qmax, q0 = marshak_q;
  // vet_col_order2 (m1-sp-order2): the core rays start from E and F at the inner FACE
  // (E linearly extrapolated from the two bottom cells, F the stored face flux f0x1)
  const bool o2 = vcol_o2;
  auto f0f_ = f0x1;
  const bool g2 = vcol_o2 && vcol_sph;
  // vet_col_order2 with a reflecting inner x1 boundary: the core rays start from the
  // MIRRORED incoming intensity at the inner face (specular reflection, H = 0 there)
  const bool mir = vcol_o2 && (ibc_x1min == M1_IBC_REFLECT);
  auto gb_ = vcol_gb;
  auto gb0_ = vcol_gb0;
  const int nray = vcol_nray;
  if (vcol_team) {
    VetColBuildTeam(dmp);
  } else {
  par_for("m1_vcol", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int c = (m*nk + (k - ks))*nj + (j - js);
    // extinction and source of shell l
    auto chx = [&](const int l) {return fmax(iw_(m,M1_IW_KT,k,j,is+l), 1.0e-300);};
    auto src = [&](const int l) {
      const int i = is + l;
      Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
      Real s = e;
      if (thermal) {
        Real chi = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
        Real tg = iw_(m,M1_IW_TP,k,j,i);
        Real t2 = tg*tg;
        Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chi, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
      return s;
    };
    // the source at the two boundary faces, linearly extrapolated (vet_sc)
    const Real sbot = fmax(1.5*src(0) - 0.5*src(1), 0.0);
    const Real stop = fmax(1.5*src(n1-1) - 0.5*src(n1-2), 0.0);

    // (1) incoming rays, top down
    Real chu = 0.0, su = 0.0;
    for (int l = n1 - 1; l >= 0; --l) {
      const Real ch0 = chx(l), s0 = src(l);
      const bool top = (l == n1 - 1);
      const Real cseg = top ? ch0 : 0.5*(chu + ch0);
      const Real sup = top ? stop : su;
      const int kr = kl_(l);
      Real jj = 0.0, hh = 0.0, kk = 0.0;
      for (int r = 0; r <= kr; ++r) {
        const Real iu = top ? 0.0 : buf_(r,c);
        Real ex, w0, wu;
        VcolW(cseg*seg_(r,l), ex, w0, wu);
        if (g2) {
          // vet_col_order2: chi and S linear in r along the (incoming) segment
          const Real gb = gb_(r,l), chup = top ? ch0 : chu;
          VcolW2((chup + (ch0 - chup)*(1.0 - gb))*seg_(r,l), -VcolC2(gb), ex, w0, wu);
        }
        const Real iv = fmax(iu*ex + wu*sup + w0*s0, 0.0);
        buf_(r,c) = iv;
        const Real wq = w_(r,l)*iv, mq = mu_(r,l);
        jj += wq;
        hh -= wq*mq;
        kk += wq*mq*mq;
      }
      tt_(m,1,k,j,is+l) = jj;
      tt_(m,2,k,j,is+l) = hh;
      tt_(m,3,k,j,is+l) = kk;
      chu = ch0;
      su = s0;
    }

    // (2) outgoing rays, bottom up
    Real e0 = fmax(iw_(m,M1_IW_EN,k,j,is), efl);
    Real f0 = iw_(m,M1_IW_F1,k,j,is);
    if (o2) {
      e0 = fmax(1.5*iw_(m,M1_IW_EN,k,j,is) - 0.5*iw_(m,M1_IW_EN,k,j,is+1), efl);
      f0 = f0f_(m,k,j,is);
    }
    Real chd = 0.0, sd = 0.0, jtop = 0.0;
    for (int l = 0; l < n1; ++l) {
      const int i = is + l;
      const Real ch0 = chx(l), s0 = src(l);
      const Real chlo = (l == 0) ? ch0 : chd;     // chi below the shell (face: cell 0)
      const Real slo = (l == 0) ? sbot : sd;
      const int kr = kl_(l);
      Real jj = tt_(m,1,k,j,i), hh = tt_(m,2,k,j,i), kk = tt_(m,3,k,j,i);
      for (int r = 0; r <= kr; ++r) {
        const int lr = static_cast<int>(ray_(r,0));
        Real iv;
        Real ex, w0, wu;
        if (lr < l) {
          VcolW(0.5*(chd + ch0)*seg_(r,l-1), ex, w0, wu);
          if (g2) {
            const Real gb = gb_(r,l-1);
            VcolW2((chd + (ch0 - chd)*gb)*seg_(r,l-1), VcolC2(gb), ex, w0, wu);
          }
          iv = buf_(r,c)*ex + wu*sd + w0*s0;
        } else {
          const int ty = static_cast<int>(ray_(r,1));
          if (ty == VC_CORE) {
            // the diffusion intensity at the inner face, carried to the first shell
            Real ib = fmax(e0 + 3.0*f0*ray_(r,3)/cl, 0.0);
            VcolW(ch0*ray_(r,2), ex, w0, wu);
            if (g2) {VcolW2(ch0*ray_(r,2), -VcolC2(gb0_(r)), ex, w0, wu);}
            if (mir) {ib = fmax(buf_(r,c)*ex + wu*s0 + w0*sbot, 0.0);}
            if (g2) {VcolW2(ch0*ray_(r,2), VcolC2(gb0_(r)), ex, w0, wu);}
            iv = ib*ex + wu*sbot + w0*s0;
          } else if (ty == VC_TAN) {
            iv = buf_(r,c);            // mu = 0 at its tangent shell: out = in
          } else {
            // down to the tangent point inside the interval below and back up
            const Real a = ray_(r,3);
            const Real cht = (1.0 - a)*chlo + a*ch0;
            const Real st = (1.0 - a)*slo + a*s0;
            VcolW(0.5*(ch0 + cht)*ray_(r,2), ex, w0, wu);
            if (g2) {
              const Real gb = gb0_(r);
              const Real dt2 = (cht + (ch0 - cht)*gb)*ray_(r,2);
              VcolW2(dt2, -VcolC2(gb), ex, w0, wu);
              const Real it = fmax(buf_(r,c)*ex + wu*s0 + w0*st, 0.0);
              VcolW2(dt2, VcolC2(gb), ex, w0, wu);
              iv = it*ex + wu*st + w0*s0;
            } else {
            const Real it = fmax(buf_(r,c)*ex + wu*s0 + w0*st, 0.0);
            iv = it*ex + wu*st + w0*s0;
            }
          }
        }
        iv = fmax(iv, 0.0);
        buf_(r,c) = iv;
        const Real wq = w_(r,l)*iv, mq = mu_(r,l);
        jj += wq;
        hh += wq*mq;
        kk += wq*mq*mq;
      }
      jj *= 0.5;
      hh *= 0.5;
      kk *= 0.5;
      Real chi = 1.0/3.0;
      if (jj > 0.0) {chi = fmin(fmax(kk/jj, 1.0/3.0), 1.0);}
      Real n1v = 1.0, n2v = 0.0, n3v = 0.0;
      if (axf) {
        Real fa = iw_(m,M1_IW_F1,k,j,i), fb = iw_(m,M1_IW_F2,k,j,i);
        Real fc = thrd ? iw_(m,M1_IW_F3,k,j,i) : 0.0;
        Real fm = sqrt(fa*fa + fb*fb + fc*fc);
        if (fm > 1.0e-8*cl*fmax(iw_(m,M1_IW_EN,k,j,i), efl)) {
          n1v = fa/fm;
          n2v = fb/fm;
          n3v = fc/fm;
        }
      }
      tt_(m,0,k,j,i) = chi;
      tt_(m,1,k,j,i) = n1v;
      tt_(m,2,k,j,i) = n2v;
      tt_(m,3,k,j,i) = n3v;
      if (dmp) {
        mo_(m,0,k,j,i) = jj;
        mo_(m,1,k,j,i) = hh;
        mo_(m,2,k,j,i) = kk;
        mo_(m,3,k,j,i) = s0;
        mo_(m,4,k,j,i) = ch0;
      }
      chd = ch0;
      sd = s0;
      jtop = jj;
    }
    // vet_col_surface_q: the outgoing intensities carried from the top shell to the top
    // face (the incoming sweep's top segment), q = H(face)/J(top cell)
    if (sq) {
      const Real ch0 = chx(n1-1), s0 = src(n1-1);
      Real hf = 0.0;
      for (int r = 0; r < nray; ++r) {
        Real ex, w0, wu;
        VcolW(ch0*seg_(r,n1-1), ex, w0, wu);
        if (g2) {VcolW2(ch0*seg_(r,n1-1), VcolC2(gb_(r,n1-1)), ex, w0, wu);}
        const Real iv = fmax(buf_(r,c)*ex + wu*s0 + w0*stop, 0.0);
        hf += wf_(r)*muf_(r)*iv;
      }
      hf *= 0.5;
      vq_(m,k,j) = (jtop > 0.0) ? fmin(fmax(hf/jtop, qlo), qhi) : q0;
    }
  });
  }
  Kokkos::fence();
  vcol_time += timer.seconds();
  if (dmp) {VetColDumpColumn(static_cast<int>(vcol_ncall));}
  vcol_ncall += 1.0;
  vcol_built = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetColBuildTeam
//! \brief vet_col_team = true: the build of VetColBuild with ONE TEAM PER COLUMN.  The
//! column profile (chi, S) sits in team scratch; the rays of a shell are spread over the
//! team's threads (each ray's recurrence stays sequential in r, one thread per ray and
//! shell), the running intensity per ray in scratch.  The intensities of a CHUNK of
//! vcol_lc shells are kept in scratch and the moments of each shell are then summed by
//! one thread in the ray order r = 0, 1, ... of the one-thread-per-column kernel: every
//! operation and its order per ray and per moment is the same, so the tensor (and q) is
//! bitwise the one of the column kernel (up to the compiler's FMA contraction).

void RadiationM1::VetColBuildTeam(bool dmp) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const bool thermal = (pmy_pack->phydro != nullptr) && coupling && !opac_zero;
  const bool axf = vcol_axis_flux;
  const Real cl = c_light, ar = arad, efl = e_floor;
  auto iw_ = iw;
  auto opac_ = opac;
  auto tt_ = tau_ten;
  auto mo_ = vcol_mom;
  auto seg_ = vcol_seg;
  auto mu_ = vcol_mu;
  auto w_ = vcol_w;
  auto ray_ = vcol_ray;
  auto kl_ = vcol_klast;
  auto muf_ = vcol_muf;
  auto wf_ = vcol_wf;
  auto vq_ = vcol_q;
  const bool sq = vcol_sq;
  const Real qlo = vcol_qmin, qhi = vcol_qmax, q0 = marshak_q;
  // vet_col_order2 (m1-sp-order2): the core rays start from E and F at the inner FACE
  // (E linearly extrapolated from the two bottom cells, F the stored face flux f0x1)
  const bool o2 = vcol_o2;
  auto f0f_ = f0x1;
  const bool g2 = vcol_o2 && vcol_sph;
  // vet_col_order2 with a reflecting inner x1 boundary: the core rays start from the
  // MIRRORED incoming intensity at the inner face (specular reflection, H = 0 there)
  const bool mir = vcol_o2 && (ibc_x1min == M1_IBC_REFLECT);
  auto gb_ = vcol_gb;
  auto gb0_ = vcol_gb0;
  const int nray = vcol_nray, lc = vcol_lc;
  const size_t scr = ScrArray1D<Real>::shmem_size(5*n1) +
                     ScrArray1D<Real>::shmem_size(2*nray) +
                     ScrArray2D<Real>::shmem_size(lc, nray);

  const int nk = ke - ks + 1, nj = je - js + 1;
  const int nlg = (nmb1 + 1)*nk*nj;
  Kokkos::TeamPolicy<> pol =
      (vcol_ts > 0) ? Kokkos::TeamPolicy<>(DevExeSpace(), nlg, vcol_ts)
                    : Kokkos::TeamPolicy<>(DevExeSpace(), nlg, Kokkos::AUTO);
  Kokkos::parallel_for("m1_vcol_team", pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/(nk*nj);
    const int k = (tm.league_rank() - m*nk*nj)/nj + ks;
    const int j = (tm.league_rank() - m*nk*nj) % nj + js;
    ScrArray1D<Real> pr_(tm.team_scratch(0), 5*n1);    // chi, S, J_in, H_in, K_in
    ScrArray1D<Real> ir_(tm.team_scratch(0), 2*nray);  // running I; face I
    ScrArray2D<Real> ic_(tm.team_scratch(0), lc, nray);
    par_for_inner(tm, 0, n1-1, [&](const int l) {
      const int i = is + l;
      pr_(l) = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
      Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
      Real s = e;
      if (thermal) {
        Real chi = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
        Real tg = iw_(m,M1_IW_TP,k,j,i);
        Real t2 = tg*tg;
        Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chi, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
      pr_(n1 + l) = s;
    });
    tm.team_barrier();
    const Real sbot = fmax(1.5*pr_(n1) - 0.5*pr_(n1+1), 0.0);
    const Real stop = fmax(1.5*pr_(2*n1-1) - 0.5*pr_(2*n1-2), 0.0);

    // (1) incoming rays, top down, in chunks of lc shells
    for (int la = n1 - 1; la >= 0; la -= lc) {
      const int lb = (la - lc + 1 > 0) ? (la - lc + 1) : 0;
      for (int l = la; l >= lb; --l) {
        const Real ch0 = pr_(l), s0 = pr_(n1 + l);
        const bool top = (l == n1 - 1);
        const Real cseg = top ? ch0 : 0.5*(pr_(l+1) + ch0);
        const Real sup = top ? stop : pr_(n1 + l + 1);
        const int ll = la - l;
        par_for_inner(tm, 0, kl_(l), [&](const int r) {
          const Real iu = top ? 0.0 : ir_(r);
          Real ex, w0, wu;
          VcolW(cseg*seg_(r,l), ex, w0, wu);
          if (g2) {
            const Real gb = gb_(r,l), chup = top ? ch0 : pr_(l+1);
            VcolW2((chup + (ch0 - chup)*(1.0 - gb))*seg_(r,l), -VcolC2(gb), ex, w0, wu);
          }
          const Real iv = fmax(iu*ex + wu*sup + w0*s0, 0.0);
          ir_(r) = iv;
          ic_(ll, r) = iv;
        });
        tm.team_barrier();
      }
      par_for_inner(tm, lb, la, [&](const int l) {
        const int ll = la - l;
        Real jj = 0.0, hh = 0.0, kk = 0.0;
        const int kr = kl_(l);
        for (int r = 0; r <= kr; ++r) {
          const Real wq = w_(r,l)*ic_(ll, r), mq = mu_(r,l);
          jj += wq;
          hh -= wq*mq;
          kk += wq*mq*mq;
        }
        pr_(2*n1 + l) = jj;
        pr_(3*n1 + l) = hh;
        pr_(4*n1 + l) = kk;
      });
      tm.team_barrier();
    }

    // (2) outgoing rays, bottom up, in chunks of lc shells
    Real e0 = fmax(iw_(m,M1_IW_EN,k,j,is), efl);
    Real f0 = iw_(m,M1_IW_F1,k,j,is);
    if (o2) {
      e0 = fmax(1.5*iw_(m,M1_IW_EN,k,j,is) - 0.5*iw_(m,M1_IW_EN,k,j,is+1), efl);
      f0 = f0f_(m,k,j,is);
    }
    for (int la = 0; la < n1; la += lc) {
      const int lb = (la + lc - 1 < n1 - 1) ? (la + lc - 1) : (n1 - 1);
      for (int l = la; l <= lb; ++l) {
        const Real ch0 = pr_(l), s0 = pr_(n1 + l);
        const Real chd = (l == 0) ? 0.0 : pr_(l-1);
        const Real sd = (l == 0) ? 0.0 : pr_(n1 + l - 1);
        const Real chlo = (l == 0) ? ch0 : chd;
        const Real slo = (l == 0) ? sbot : sd;
        const int ll = l - la;
        par_for_inner(tm, 0, kl_(l), [&](const int r) {
          const int lr = static_cast<int>(ray_(r,0));
          Real iv;
          Real ex, w0, wu;
          if (lr < l) {
            VcolW(0.5*(chd + ch0)*seg_(r,l-1), ex, w0, wu);
            if (g2) {
              const Real gb = gb_(r,l-1);
              VcolW2((chd + (ch0 - chd)*gb)*seg_(r,l-1), VcolC2(gb), ex, w0, wu);
            }
            iv = ir_(r)*ex + wu*sd + w0*s0;
          } else {
            const int ty = static_cast<int>(ray_(r,1));
            if (ty == VC_CORE) {
              Real ib = fmax(e0 + 3.0*f0*ray_(r,3)/cl, 0.0);
              VcolW(ch0*ray_(r,2), ex, w0, wu);
              if (g2) {VcolW2(ch0*ray_(r,2), -VcolC2(gb0_(r)), ex, w0, wu);}
              if (mir) {ib = fmax(ir_(r)*ex + wu*s0 + w0*sbot, 0.0);}
              if (g2) {VcolW2(ch0*ray_(r,2), VcolC2(gb0_(r)), ex, w0, wu);}
              iv = ib*ex + wu*sbot + w0*s0;
            } else if (ty == VC_TAN) {
              iv = ir_(r);
            } else {
              const Real a = ray_(r,3);
              const Real cht = (1.0 - a)*chlo + a*ch0;
              const Real st = (1.0 - a)*slo + a*s0;
              VcolW(0.5*(ch0 + cht)*ray_(r,2), ex, w0, wu);
              if (g2) {
                const Real gb = gb0_(r);
                const Real dt2 = (cht + (ch0 - cht)*gb)*ray_(r,2);
                VcolW2(dt2, -VcolC2(gb), ex, w0, wu);
                const Real it = fmax(ir_(r)*ex + wu*s0 + w0*st, 0.0);
                VcolW2(dt2, VcolC2(gb), ex, w0, wu);
                iv = it*ex + wu*st + w0*s0;
              } else {
              const Real it = fmax(ir_(r)*ex + wu*s0 + w0*st, 0.0);
              iv = it*ex + wu*st + w0*s0;
              }
            }
          }
          iv = fmax(iv, 0.0);
          ir_(r) = iv;
          ic_(ll, r) = iv;
        });
        tm.team_barrier();
      }
      if (sq && lb == n1 - 1) {
        // vet_col_surface_q: the outgoing intensities carried to the top face
        const Real ch0 = pr_(n1-1), s0 = pr_(2*n1-1);
        par_for_inner(tm, 0, nray-1, [&](const int r) {
          Real ex, w0, wu;
          VcolW(ch0*seg_(r,n1-1), ex, w0, wu);
          if (g2) {VcolW2(ch0*seg_(r,n1-1), VcolC2(gb_(r,n1-1)), ex, w0, wu);}
          ir_(nray + r) = fmax(ir_(r)*ex + wu*s0 + w0*stop, 0.0);
        });
        tm.team_barrier();
      }
      par_for_inner(tm, la, lb, [&](const int l) {
        const int i = is + l;
        const int ll = l - la;
        Real jj = pr_(2*n1 + l), hh = pr_(3*n1 + l), kk = pr_(4*n1 + l);
        const int kr = kl_(l);
        for (int r = 0; r <= kr; ++r) {
          const Real wq = w_(r,l)*ic_(ll, r), mq = mu_(r,l);
          jj += wq;
          hh += wq*mq;
          kk += wq*mq*mq;
        }
        jj *= 0.5;
        hh *= 0.5;
        kk *= 0.5;
        Real chi = 1.0/3.0;
        if (jj > 0.0) {chi = fmin(fmax(kk/jj, 1.0/3.0), 1.0);}
        Real n1v = 1.0, n2v = 0.0, n3v = 0.0;
        if (axf) {
          Real fa = iw_(m,M1_IW_F1,k,j,i), fb = iw_(m,M1_IW_F2,k,j,i);
          Real fc = thrd ? iw_(m,M1_IW_F3,k,j,i) : 0.0;
          Real fm = sqrt(fa*fa + fb*fb + fc*fc);
          if (fm > 1.0e-8*cl*fmax(iw_(m,M1_IW_EN,k,j,i), efl)) {
            n1v = fa/fm;
            n2v = fb/fm;
            n3v = fc/fm;
          }
        }
        tt_(m,0,k,j,i) = chi;
        tt_(m,1,k,j,i) = n1v;
        tt_(m,2,k,j,i) = n2v;
        tt_(m,3,k,j,i) = n3v;
        if (dmp) {
          mo_(m,0,k,j,i) = jj;
          mo_(m,1,k,j,i) = hh;
          mo_(m,2,k,j,i) = kk;
          mo_(m,3,k,j,i) = pr_(n1 + l);
          mo_(m,4,k,j,i) = pr_(l);
        }
        if (sq && l == n1 - 1) {
          Real hf = 0.0;
          for (int r = 0; r < nray; ++r) {hf += wf_(r)*muf_(r)*ir_(nray + r);}
          hf *= 0.5;
          vq_(m,k,j) = (jj > 0.0) ? fmin(fmax(hf/jj, qlo), qhi) : q0;
        }
      });
      tm.team_barrier();
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetColDumpColumn
//! \brief rank 0, MeshBlock 0, column (ks, js): r, chi, S, J, H, K, f_K, E^n, F1 (the
//! diagnostic of T-S6).  Not timed.

void RadiationM1::VetColDumpColumn(int ncall) {
  if (global_variable::my_rank != 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, n1 = indcs.nx1;
  auto mo_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vcol_mom);
  auto iw_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), iw);
  auto tt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tau_ten);
  char nm[32];
  std::snprintf(nm, sizeof(nm), ".%05d.txt", ncall);
  std::ofstream f(vcol_dump + nm);
  f << "# vet_col column dump: cycle " << pmy_pack->pmesh->ncycle << " time "
    << std::setprecision(17) << pmy_pack->pmesh->time << "\n"
    << "# r chi S J H K fK chi_used E_n F1\n";
  f << std::setprecision(17);
  for (int l = 0; l < n1; ++l) {
    const int i = is + l;
    Real jj = mo_h(0,0,ks,js,i);
    f << vcol_rc[l] << " " << mo_h(0,4,ks,js,i) << " " << mo_h(0,3,ks,js,i) << " "
      << jj << " " << mo_h(0,1,ks,js,i) << " " << mo_h(0,2,ks,js,i) << " "
      << ((jj > 0.0) ? mo_h(0,2,ks,js,i)/jj : 0.0) << " " << tt_h(0,0,ks,js,i) << " "
      << iw_h(0,M1_IW_EN,ks,js,i) << " " << iw_h(0,M1_IW_F1,ks,js,i) << "\n";
  }
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
}

} // namespace radm1
