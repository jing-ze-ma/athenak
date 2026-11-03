#ifndef UTILS_DIVDISS_HPP_
#define UTILS_DIVDISS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file divdiss.hpp
//! \brief lowmach2-1102: multi-dimensional divergence dissipation for the low-Mach
//! Riemann solvers LHLLC / LHLLD (keys <hydro>/lhllc_divdiss, <mhd>/lhlld_divdiss).
//!
//! LHLLC/LHLLD scale the velocity-jump (acoustic) term of the contact pressure,
//!   p* = pbar - [m_L m_R/(m_L + m_R)] (u_R - u_L),       m_K = rho_K |S_K - u_K|,
//! by phi = O(M).  That term is the only part of the flux that damps a 2-cell
//! compressive mode of the normal velocity (the Rhie-Chow role), so at phi = O(M)
//! the mode survives (the reason for the checkerboard sensor, utils/cb_sensor.hpp).
//! Barsukow (2021, arXiv:2103.02621, eqs. 17-18 and 27-30; 2023, arXiv:2306.02847)
//! shows that a scheme is low-Mach compliant (stationarity preserving) when the
//! acoustic dissipation of the face pressure acts on a MULTI-DIMENSIONAL discrete
//! divergence, not on the 1-D normal jump: then it vanishes on discretely
//! divergence-free states and keeps the full acoustic rate on compressive ones.
//! Here the missing (1 - phi) of the acoustic term is restored acting on such a
//! divergence, of the MASS flux rho u (in a stratified medium div(rho u) ~ 0 is the
//! anelastic target, div u is not small):
//!   dp* = -(1 - phi) [m_L m_R/(m_L + m_R)] J / rhobar,   rhobar = (rho_L + rho_R)/2,
//! J = an undivided (x dx_n) face divergence of rho u built from the cell-centred w0.
//! J/rhobar replaces (u_R - u_L), so on a pure normal checkerboard the face gets
//! exactly the plain HLLC/HLLD dissipation.  Linear in the state, no switch, no
//! parameter; zero at phi = 1.  Four forms (key <..>_divdiss_form; default hpd):
//!
//!  raw (2)  Barsukow's corner divergence: J = dx_n div_face(rho u) with the normal
//!           difference averaged over the transverse neighbours with weights
//!           (1/4, 1/2, 1/4) per active transverse direction, and the transverse
//!           centred differences averaged over the two cells of the face (and smoothed
//!           (1/4,1/2,1/4) in the third direction).  = the average of the two
//!           node divergences at the face.  Kernel: node-divergence-free fields, which
//!           include u = (-1)^(i+j), v = 0.  FIRST ORDER on compressive smooth flow
//!           (damps sound at rate c k^2 dx).
//!  hp (1)   the Rhie-Chow high-pass of raw: J = dx_n [div_face - (div_c(L) +
//!           div_c(R))/2], div_c the cell divergence with the same transverse stencil.
//!           The transverse parts cancel identically and J reduces to the
//!           transverse-smoothed (1/4,1/2,1/4) odd part of the normal difference,
//!           J = (q_{-2} - 3 q_{-1} + 3 q_0 - q_1)/4 per row (the numerator of the
//!           checkerboard sensor), = O(dx^3 q''') on smooth flow (third order, the
//!           sound wave keeps its order).  Same kernel mode u = (-1)^(i+j) as raw.
//!  row (3)  hp without the transverse smoothing (face row only): damps every 2dx
//!           compressive mode including u = (-1)^(i+j); purely 1-D.
//!
//!  hpd (4) the normal-direction high-pass of the raw DIVERGENCE (not of the flux):
//!           J = (2 J_raw(f) - J_raw(f-1) - J_raw(f+1))/4 with f-1, f+1 the neighbouring
//!           faces along the normal.  Vanishes wherever the corner divergence does
//!           (discretely incompressible flow, at any smoothness: kinks in the flux of a
//!           solenoidal field cancel between the normal and transverse parts, which the
//!           flux high-pass hp cannot see), O(dx^3) on smooth compressive flow (sound keeps
//!           its order), and = J_raw for a normal checkerboard.  Kernel: that of raw
//!           plus the linear-in-x divergence fields.
//!
//!  hpdk (5) hpd + (row - hp): adds the transverse 2dx part of the face-row odd part, the
//!           only component the corner-averaged forms cannot see (kernel mode
//!           u = (-1)^(i+j), which grows under hpd/raw at M 1e-4 over a turnover).
//!
//! Cartesian prototype: on curvilinear grids w0(ivx) is the coordinate component and
//! raw would need the metric (refused there by the caller's mesh check: TODO).

#include "athena.hpp"

KOKKOS_INLINE_FUNCTION
Real DivDissQ(const DvceArray5D<Real> &w0, const int m, const int var, int p0, int p1,
              int p2) {
  const int n0 = static_cast<int>(w0.extent(4)) - 1;
  const int n1 = static_cast<int>(w0.extent(3)) - 1;
  const int n2 = static_cast<int>(w0.extent(2)) - 1;
  p0 = (p0 < 0) ? 0 : ((p0 > n0) ? n0 : p0);
  p1 = (p1 < 0) ? 0 : ((p1 > n1) ? n1 : p1);
  p2 = (p2 < 0) ? 0 : ((p2 > n2) ? n2 : p2);
  return w0(m,IDN,p2,p1,p0)*w0(m,var,p2,p1,p0);
}

//! J (units of rho*u) at face (k,j,i) of normal ivx; the face sits between cells
//! (normal index) -1 and 0.  act2/act3: x2/x3 active; dx[3] the cell widths.
KOKKOS_INLINE_FUNCTION
Real DivDissJ0(const DvceArray5D<Real> &w0, const int m, const int k, const int j,
               const int i, const int ivx, const int form, const bool act2,
               const bool act3, const Real dx1, const Real dx2, const Real dx3) {
  const int dn = ivx - IVX;
  const int t1 = (dn + 1) % 3, t2 = (dn + 2) % 3;
  const bool act[3] = {true, act2, act3};
  const Real dxa[3] = {dx1, dx2, dx3};
  const Real wt[3] = {0.25, 0.5, 0.25};
  const bool a1 = act[t1] && (form != 3);
  const bool a2 = act[t2] && (form != 3);
  Real jsum = 0.0;
  for (int b = -1; b <= 1; ++b) {
    if (!a1 && b != 0) continue;
    const Real wb = a1 ? wt[b+1] : 1.0;
    for (int c = -1; c <= 1; ++c) {
      if (!a2 && c != 0) continue;
      const Real wc = a2 ? wt[c+1] : 1.0;
      int p[3] = {i, j, k};
      p[t1] += b; p[t2] += c;
      Real q[4];
      for (int o = 0; o < 4; ++o) {
        int pp[3] = {p[0], p[1], p[2]};
        pp[dn] += o - 2;
        q[o] = DivDissQ(w0, m, ivx, pp[0], pp[1], pp[2]);
      }
      const Real jr = (form == 2) ? (q[2] - q[1])
                                  : 0.25*((q[0] - q[3]) + 3.0*(q[2] - q[1]));
      jsum += wb*wc*jr;
    }
  }
  if (form != 2) return jsum;
  // raw: transverse centred differences, averaged over the two cells of the face and
  // smoothed (1/4,1/2,1/4) along the remaining transverse direction
  for (int s = 0; s < 2; ++s) {
    const int td = (s == 0) ? t1 : t2;     // differenced direction
    const int to = (s == 0) ? t2 : t1;     // smoothed direction
    if (!act[td]) continue;
    const int ivt = IVX + td;
    Real tsum = 0.0;
    for (int a = -1; a <= 0; ++a) {
      for (int c = -1; c <= 1; ++c) {
        if (!act[to] && c != 0) continue;
        const Real wc = act[to] ? wt[c+1] : 1.0;
        int pp[3] = {i, j, k};
        pp[dn] += a; pp[to] += c;
        int pu[3] = {pp[0], pp[1], pp[2]}; pu[td] += 1;
        int pd[3] = {pp[0], pp[1], pp[2]}; pd[td] -= 1;
        tsum += 0.5*wc*0.5*(DivDissQ(w0, m, ivt, pu[0], pu[1], pu[2])
                            - DivDissQ(w0, m, ivt, pd[0], pd[1], pd[2]));
      }
    }
    jsum += (dxa[dn]/dxa[td])*tsum;
  }
  return jsum;
}


//! dispatcher: forms 1-3 directly, 4 (hpd) = normal second difference of raw (form 2)
KOKKOS_INLINE_FUNCTION
Real DivDissJ(const DvceArray5D<Real> &w0, const int m, const int k, const int j,
              const int i, const int ivx, const int form, const bool act2,
              const bool act3, const Real dx1, const Real dx2, const Real dx3) {
  if (form < 4) return DivDissJ0(w0, m, k, j, i, ivx, form, act2, act3, dx1, dx2, dx3);
  const int di = (ivx == IVX) ? 1 : 0, dj = (ivx == IVY) ? 1 : 0, dk = (ivx == IVZ) ? 1 : 0;
  const Real jc = DivDissJ0(w0, m, k, j, i, ivx, 2, act2, act3, dx1, dx2, dx3);
  const Real jm = DivDissJ0(w0, m, k-dk, j-dj, i-di, ivx, 2, act2, act3, dx1, dx2, dx3);
  const Real jp = DivDissJ0(w0, m, k+dk, j+dj, i+di, ivx, 2, act2, act3, dx1, dx2, dx3);
  Real jd = 0.25*(2.0*jc - jm - jp);
  // hpdk (5): + the TRANSVERSE high-pass of the face-row odd part, row - corner-smoothed
  // (= J_row for the kernel mode u = (-1)^(i+j), O(dx^2) along a face for smooth flow)
  if (form == 5) {
    jd += DivDissJ0(w0, m, k, j, i, ivx, 3, act2, act3, dx1, dx2, dx3)
        - DivDissJ0(w0, m, k, j, i, ivx, 1, act2, act3, dx1, dx2, dx3);
  }
  return jd;
}

#endif // UTILS_DIVDISS_HPP_
