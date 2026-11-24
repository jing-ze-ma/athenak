//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
// ceiling_cell.cpp: cell-level unit test of the velocity ceiling / KE floor
// spherical polar (orthonormal), Cartesian and the cubed sphere (covariant m_a on the
// non-orthogonal unit tangents).  A known physical velocity V = A phi_hat + B r_hat +
// C theta_hat goes through the real code paths (SingleC2P_GeneralHyd,
// GnomonicRaiseVelFloors); the result is rebuilt as a Cartesian vector and checked:
//   E fixed; V'.phi_hat = A (L_z); |V'| = max(vceil, |A|); V' - A phi_hat parallel to
//   V - A phi_hat; eint' - eint = 1/2 rho (|V|^2 - |V'|^2) (true metric KE change);
//   the code's KE 1/2 m'_a v'^a = 1/2 rho |V'|^2; sp == cs == Cartesian V'.
// Header-only: build against a configured CPU build B (fastbuild slot, gcc 13 + openmpi)
//   mpicxx -std=c++17 -O2 -ffp-contract=off -fopenmp -I src -I B -I B/kokkos \
//     -I B/kokkos/core/src -I kokkos/core/src -I kokkos/tpls/desul/include \
//     -I kokkos/tpls/mdspan/include -I B/kokkos/containers/src \
//     -I kokkos/containers/src tests/vet_clean/unit/ceiling_cell.cpp \
//     B/kokkos/core/src/libkokkoscore.a -ldl -o ceiling_cell && ./ceiling_cell
#include <cmath>
#include <cstdio>
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_hyd.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/gnomonic_raisevel.hpp"

static double worst[16];
static const char *wname[16] = {"dE", "dLz(u_phi)", "|V'|-target", "perp parallel",
  "heat-dKE", "codeKE-|V'|^2", "sp-vs-cs V'", "cart-vs-cs V'", "metric |phi|-1",
  "cos_cell", "keyoff |V'|-vceil", "tfloor eint-etf", "tfloor dLz", "tfloor dE",
  "GnomonicPhiHat", "spare"};
static void rec(int k, double v) {
  v = std::fabs(v);
  if (v > worst[k]) worst[k] = v;
}

static double dot(const double *a, const double *b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

int main() {
  EOS_Data eos{};
  eos.gamma = 5.0/3.0;
  eos.is_ideal = true;
  eos.eos_type = EOSType::general;
  eos.dfloor = 1e-30; eos.pfloor = 1e-30; eos.tfloor = 1e-30; eos.sfloor = 0.0;
  eos.ekin_heat = true;
  eos.defer_cons_floors = false;
  const double d = 1.3;
  int ncase = 0;
  const double pts[][2] = {{0.0, 0.0}, {0.3, -0.2}, {0.7, 0.7}, {-0.75, 0.6},
                           {0.88, 0.1}, {-0.9, -0.85}, {0.2, 0.92}};   // incl. ghosts
  // (A, B, C) and vceil: |A| < vceil < |V|, |A| > vceil, tiny A
  const double vel[][4] = {{0.5, 1.0, -0.8, 1.0}, {2.0, 1.0, 0.8, 1.0},
                           {-1.7, -0.6, 1.4, 1.0}, {1e-3, 3.0, 2.0, 0.5},
                           {0.9, 0.0, 0.0, 0.5}};
  for (int p = 0; p < 6; ++p) for (auto &pt : pts) for (auto &vv : vel) {
    const double xi = pt[0], eta = pt[1];
    double q[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(p, xi, eta, q);
    cubed_sphere::PanelTangents(p, xi, eta, e1, e2);
    const double rr = std::sqrt(q[0]*q[0] + q[1]*q[1]);
    if (rr < 1e-3) continue;   // (the axis: no phi_hat)
    const double ph[3] = {-q[1]/rr, q[0]/rr, 0.0};
    const double th[3] = {ph[1]*q[2] - ph[2]*q[1], ph[2]*q[0] - ph[0]*q[2],
                          ph[0]*q[1] - ph[1]*q[0]};   // phi_hat x r_hat (south)
    const double x = std::tan(xi), y = std::tan(eta);
    const double c = -x*y/(std::sqrt(1 + x*x)*std::sqrt(1 + y*y));   // cos_cell
    rec(9, c - dot(e1, e2));
    Real f2, f3;   // (the code's covariant phi_hat, against phi_hat.e_a)
    GnomonicPhiHat(p, xi, eta, f2, f3);
    rec(14, std::fabs(f2 - dot(ph, e1)) + std::fabs(f3 - dot(ph, e2)));
    const double det = 1 - c*c;
    rec(8, (f2*f2 + f3*f3 - 2*c*f2*f3)/det - 1.0);
    const double A = vv[0], B = vv[1], C = vv[2], vc = vv[3];
    double V[3], Vp[3];
    for (int n = 0; n < 3; ++n) V[n] = A*ph[n] + B*q[n] + C*th[n];
    const double v2 = dot(V, V), ke = 0.5*d*v2;
    const double eint0 = 0.2*ke + 0.1, et0 = eint0 + ke;
    const double target = (std::fabs(A) > vc) ? std::fabs(A) : vc;
    const bool fires = v2 > vc*vc;
    // --- cubed sphere ---
    Real m1 = d*dot(V, q), m2 = d*dot(V, e1), m3 = d*dot(V, e2), et = et0;
    Real v1, w2, w3, eint, pn, g1, tp, de, vde;
    bool cu, fl, ts;
    eos.vceil = vc;
    GnomonicRaiseVelFloors(c, eos, true, false, vc, true, 1.0, 1.0, d, m1, m2, m3, et,
                           v1, w2, w3, eint, pn, g1, tp, cu, fl, ts, de, vde, f2, f3);
    for (int n = 0; n < 3; ++n) Vp[n] = v1*q[n] + w2*e1[n] + w3*e2[n];
    const double vp2 = dot(Vp, Vp);
    rec(0, (et - et0)/et0);
    rec(1, (dot(Vp, ph) - A)/std::sqrt(v2));
    if (fires) rec(2, (std::sqrt(vp2) - target)/target);
    {
      double a[3], b[3];
      for (int n = 0; n < 3; ++n) { a[n] = V[n] - A*ph[n]; b[n] = Vp[n] - A*ph[n]; }
      const double cr[3] = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                            a[0]*b[1] - a[1]*b[0]};
      rec(3, std::sqrt(dot(cr, cr))/(dot(a, a) + 1e-300));
    }
    rec(4, ((eint - eint0) - 0.5*d*(v2 - vp2))/ke);
    rec(5, (0.5*(m1*v1 + m2*w2 + m3*w3) - 0.5*d*vp2)/ke);
    // --- spherical polar: orthonormal (r, theta, phi) ---
    HydCons1D u{d, d*dot(V, q), d*dot(V, th), d*dot(V, ph), et0};
    HydPrim1D w;
    Real temp, pg, gg, efd = 0, dfv = 1, vcd = 0;
    bool b1 = false, b2 = false, b3 = false, ms = false, vu = false, vt = false,
         tc = false, tsu = false;
    eos.ekin_keep_x3 = true;
    SingleC2P_GeneralHyd(u, eos, w, 1.0, temp, pg, gg, b1, b2, b3, efd, ms, dfv, vu, vt,
                         tc, tsu, vcd);
    double Vs[3];
    for (int n = 0; n < 3; ++n) Vs[n] = w.vx*q[n] + w.vy*th[n] + w.vz*ph[n];
    double dd = 0;
    for (int n = 0; n < 3; ++n) dd = std::fmax(dd, std::fabs(Vs[n] - Vp[n]));
    rec(6, dd/std::sqrt(v2));
    rec(0, (u.e - et0)/et0);
    // --- Cartesian, the cell at R q ---
    const double R = 2.7;
    HydCons1D uc{d, d*V[0], d*V[1], d*V[2], et0};
    eos.ekin_keep_x3 = false;
    SingleC2P_GeneralHyd(uc, eos, w, 1.0, temp, pg, gg, b1, b2, b3, efd, ms, dfv, vu, vt,
                         tc, tsu, vcd, -R*q[1]/(R*rr), R*q[0]/(R*rr));
    dd = std::fmax(std::fabs(w.vx - Vp[0]), std::fmax(std::fabs(w.vy - Vp[1]),
                                                      std::fabs(w.vz - Vp[2])));
    rec(7, dd/std::sqrt(v2));
    // --- key off: the plain ceiling, |V'| = vceil ---
    eos.ekin_heat = false;
    m1 = d*dot(V, q); m2 = d*dot(V, e1); m3 = d*dot(V, e2); et = et0;
    GnomonicRaiseVelFloors(c, eos, true, false, vc, true, 1.0, 1.0, d, m1, m2, m3, et,
                           v1, w2, w3, eint, pn, g1, tp, cu, fl, ts, de, vde, f2, f3);
    for (int n = 0; n < 3; ++n) Vp[n] = v1*q[n] + w2*e1[n] + w3*e2[n];
    if (fires) rec(10, std::sqrt(dot(Vp, Vp)) - vc);
    eos.ekin_heat = true;
    // --- temperature floor paid from the KE without L_z (no ceiling) ---
    {
      const double etf = 0.4*ke;   // floor at T: e(rho, tfloor) = rho T/(gamma - 1)
      eos.tfloor = etf*(eos.gamma - 1.0)/d;
      const double ek_phi = 0.5*d*A*A;
      const double et1 = 0.95*etf + ke;   // eint below the floor
      m1 = d*dot(V, q); m2 = d*dot(V, e1); m3 = d*dot(V, e2); et = et1;
      GnomonicRaiseVelFloors(c, eos, true, false, 0.0, true, 1.0, 1.0, d, m1, m2, m3,
                             et, v1, w2, w3, eint, pn, g1, tp, cu, fl, ts, de, vde, f2,
                             f3);
      for (int n = 0; n < 3; ++n) Vp[n] = v1*q[n] + w2*e1[n] + w3*e2[n];
      if (et1 - ek_phi >= etf) {   // payable: E fixed, eint = etf
        rec(11, (eint - etf)/etf);
        rec(13, (et - et1)/et1);
      }
      rec(12, (dot(Vp, ph) - A)/std::sqrt(v2));
      eos.tfloor = 1e-30;
    }
    ++ncase;
  }
  std::printf("VCOV unit test: %d cells x velocities (6 panels, incl. ghost-extended "
              "and corner points)\n", ncase);
  for (int k = 0; k < 15; ++k) std::printf("  %-20s %.3e\n", wname[k], worst[k]);
  bool ok = true;
  for (int k = 0; k < 15; ++k) ok = ok && worst[k] < 1e-12;
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
