//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ry_per_accretor.cpp
//! \brief The mass-gaining star of an Algol binary (RY Per first): the L1 stream hitting
//! a rigid, absorbing, rotating stellar surface, in the r-phi plane of the orbit.
//! Design: /viper/ptmp2/jinma/accretor_1006/DESIGN.md.
//! Build: -D PROBLEM=ry_per_accretor.
//!
//! GEOMETRY.  Spherical polar (mesh/use_spherical_polar), origin at the ACCRETOR centre,
//! x1 = r from R_acc to r_out (log-spaced: mesh/use_grid_stretch_r with f_stretch_r =
//! -ln(r_out/R_acc) gives r_i = R_acc (r_out/R_acc)^xi EXACTLY), x2 = a thin band of 4
//! cells about theta = pi/2 with REFLECTING faces, x3 = phi in [0, 2 pi) periodic.
//! phi = 0 points to the donor; phi > 0 is +y, the TRAILING side (Omega along +z).
//! Every force is evaluated IN THE ORBITAL PLANE (theta = pi/2): no theta
//! gravity, no theta Coriolis, sin(theta) = 1, so the band is a pure 2-D r-phi model.
//!
//! UNITS.  length = Rsun, velocity = km/s (time = Rsun/(km/s) = 6.957e5 s), density =
//! arbitrary: the run is SCALE-FREE in density (isothermal, no cooling), the code density
//! unit is the stream's peak density problem/rho_stream (default 1).  Every result
//! scales with Mdot; report ratios (j_acc = Jdot/Mdot, torque/Mdot).
//!
//! FRAME AND FORCES.  The frame co-rotates with the binary at Omega about the binary
//! centre of mass (CM), on the line of centres at x_cm = a M_d/(M_a + M_d) from
//! the accretor.  The CM is an inertial point, so a frame rotating about it at constant
//! Omega needs only the centrifugal and Coriolis terms.  Shifting the origin to
//! the accretor (a point at rest in that frame) is a pure translation, so the effective
//! potential is
//!     Phi(x,y) = -G M_a/r - G M_d/|x - a e_x| - Omega^2 [(x - x_cm)^2 + y^2]/2.
//! The INDIRECT TERM is already inside it: expanding the centrifugal part about the
//! accretor gives Omega^2 (x, y) - Omega^2 x_cm e_x, and Omega^2 x_cm = G M_d/a^2 because
//! Omega^2 = G (M_a + M_d)/a^3 (Kepler; this pgen sets Omega from the masses and a, so
//! the restricted three-body problem is exact and L1 is a true equilibrium).  The
//! constant -G M_d/a^2 e_x is minus the accretor's acceleration toward the donor, i.e.
//! the indirect term of an accretor-centred frame.  Coriolis -2 Omega x v is explicit:
//! a_r += 2 Omega v_phi, a_phi += -2 Omega v_r (no work).
//! The gravity source is the discrete gradient of Phi from its values on the cell's FACES
//! (r faces for g_r, phi faces for g_phi), precomputed once into (m,k,i) arrays.
//!
//! ENERGY.  Isothermal only (<hydro>/eos = isothermal; problem/thermo = isothermal).  The
//! key problem/thermo = adiabatic is RESERVED for the gamma = 5/3 option and is FATAL
//! today.  Building it: put Phi into phicc0/phi0 with <hydro>/etotgrav (he_star_m1.cpp
//! hs_phi fill), keep Coriolis explicit, do NOT add the centrifugal work as an explicit
//! energy source when the potential carries it (the 09-28 rot_potential double-count
//! guard, deep_hot_jupiter_rt.cpp `if (!(rotpot && use_etotgrav))`), and give the user
//! BCs a pressure.  With the isothermal EOS there is no energy equation, so neither the
//! etotgrav machinery nor that guard applies here.
//!
//! BOUNDARIES (user BCs on both x1 faces; mesh/ix1_bc = ox1_bc = user):
//!  inner (r = R_acc): rigid ABSORBING stellar surface.  v_r(ghost) = min(v_r(edge), 0)
//!    (a diode: gas only leaves into the star); density copied from the edge cell
//!    (problem/inner_rho = copy, absorbing: the surface gives no pressure support) or
//!    hydrostatically extrapolated in Phi (= hse; for the at-rest test); v_theta copied;
//!    problem/inner_vr = wall mirrors v_r instead of the diode (closed wall, test only);
//!    v_phi(ghost) = (spin - 1) Omega r (problem/inner_slip = noslip, star spin
//!    spin x Omega_orb, rotating frame) or copied (= free).
//!  outer (r = r_out): the STREAM enters through a phi window |phi - phi_s| < nsig sigma:
//!    rho = rho_stream exp(-(phi - phi_s)^2/(2 sigma^2)), v_r = v_r,s, v_phi = v_phi,s
//!    (rotating frame), from a ballistic L1 orbit integrated at start-up (same restricted
//!    three-body problem and launch as scripts/roche_stream.py: from L1 toward the
//!    accretor at the donor sound speed); sigma = w/r_out, w = problem/stream_width
//!    (default c_s,don/Omega).  Elsewhere: outflow without inflow, v_r = max(v_r, 0),
//!    v_theta, v_phi copied, density hydrostatically extrapolated (outer_rho = hse,
//!    default) or copied (= copy).
//!
//! HISTORY (problem/user_hist = true), all in code units, the plain sums over ranks:
//!   Mdom, Jdom           mass and INERTIAL z angular momentum about the accretor axis,
//!                        rho r (v_phi + Omega r), in the domain
//!   Min, Mout, Macc      CUMULATIVE mass in through the stream window, out through the
//!                        rest of r_out, into the star through r_in (time integrals of
//!                        Riemann fluxes with the RK weights of the integrator, rk1/2/3)
//!   Jacc                 cumulative inertial Lz into the star, -(F_M3 + Omega r F_D) r A
//!   Jstr                 the part of Jacc NOT carried by the edge cell's v_phi:
//!                        -(F_M3 - F_D v_phi(is)) r A.  With a no-slip wall this is the
//!                        surface stress (numerical viscosity of the Riemann solver vs
//!                        the ghost's wall velocity); approximate (cell-centred v_phi).
//!   Jout                 cumulative net inertial Lz out through r_out
//!   dMin ... dJout       the same six as RATES averaged over the last history interval
//! The cumulative integrals restart from zero on a restart.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "outputs/outputs.hpp"
#include "pgen.hpp"

namespace {

// physical constants (cgs) and the code units
constexpr Real kG = 6.674e-8, kMsun = 1.989e33, kRsun = 6.957e10;
constexpr Real kKB = 1.381e-16, kMH = 1.6726e-24;
constexpr Real kVel = 1.0e5;                 // km/s
constexpr Real kTime = kRsun/kVel;           // 6.957e5 s

struct RocheParams {
  Real gma, gmd;    // G M_a, G M_d   [(km/s)^2 Rsun]
  Real asep;        // separation     [Rsun]
  Real omega;       // orbital angular velocity [code 1/time]
  Real xcm;         // CM distance from the accretor [Rsun]
};

RocheParams rp_;
DvceArray3D<Real> gr_, gp_;    // discrete -dPhi/dr, -dPhi/(r dphi) at (m,k,i)
Real spin_ = 1.0;
bool noslip_ = true, inner_hse_ = false, outer_hse_ = true, stream_on_ = true;
bool inner_wall_ = false;
Real rho_s_ = 1.0, phi_s_ = 0.0, sig_ = 0.05, nsig_ = 3.0, vr_s_ = 0.0, vp_s_ = 0.0;
Real hse_cap_ = 40.0;
// flux accumulators (per rank), RK registers, stage counter
constexpr int kNacc = 6;
int nstages_ = 0;
Real rk_g0_[3], rk_g1_[3];
int stage_ctr_ = 0;
Real acc0_[kNacc] = {0.0}, acc1_[kNacc] = {0.0};
Real hist_prev_[kNacc] = {0.0};
Real hist_tprev_ = -1.0;

KOKKOS_INLINE_FUNCTION
Real RochePot(const RocheParams &p, const Real r, const Real phi) {
  const Real cp = cos(phi);
  const Real rd = sqrt(fmax(r*r + p.asep*p.asep - 2.0*p.asep*r*cp, 1.0e-30));
  return -p.gma/r - p.gmd/rd
         - 0.5*SQR(p.omega)*(r*r + p.xcm*p.xcm - 2.0*p.xcm*r*cp);
}

KOKKOS_INLINE_FUNCTION
Real WrapPhi(Real d) {
  const Real tp = 2.0*M_PI;
  d = fmod(d, tp);
  if (d > M_PI) d -= tp;
  if (d < -M_PI) d += tp;
  return d;
}

//----------------------------------------------------------------------------------------
// Ballistic L1 stream, restricted three-body problem in units a = 1, Omega = 1,
// G (M_a + M_d) = 1, CM origin, accretor at (-mu, 0), donor at (1 - mu, 0); the launch of
// scripts/roche_stream.py (from x_L1 - 1e-4 toward the accretor at speed eps).  RK4, a
// fixed step.  Returns the accretor-centred state where r first drops to r_cross, and the
// first periastron r_min.
struct BallisticOut {
  bool ok_cross = false;
  Real phi = 0.0, vr = 0.0, vt = 0.0;   // at r_cross (rad, units of a Omega)
  Real rmin = 0.0, phimin = 0.0;
  Real dl1 = 0.0;
};

void BallisticRhs(const Real mu, const Real *s, Real *ds) {
  const Real xa = -mu, xd = 1.0 - mu;
  const Real x = s[0], y = s[1], vx = s[2], vy = s[3];
  const Real r1 = std::hypot(x - xa, y), r2 = std::hypot(x - xd, y);
  const Real r13 = r1*r1*r1, r23 = r2*r2*r2;
  ds[0] = vx;
  ds[1] = vy;
  ds[2] = -(1.0 - mu)*(x - xa)/r13 - mu*(x - xd)/r23 + x + 2.0*vy;
  ds[3] = -(1.0 - mu)*y/r13 - mu*y/r23 + y - 2.0*vx;
}

BallisticOut IntegrateStream(const Real mu, const Real eps, const Real r_cross) {
  BallisticOut out;
  const Real xa = -mu, xd = 1.0 - mu;
  auto gx = [&](Real x) {
    return -(1.0 - mu)*(x - xa)/std::pow(std::fabs(x - xa), 3)
           - mu*(x - xd)/std::pow(std::fabs(x - xd), 3) + x;
  };
  // L1 by bisection (gx > 0 just right of the accretor? sign check both ways)
  Real lo = xa + 1.0e-3, hi = xd - 1.0e-3;
  Real flo = gx(lo);
  for (int it=0; it<200; ++it) {
    const Real mid = 0.5*(lo + hi);
    const Real fm = gx(mid);
    if ((fm > 0.0) == (flo > 0.0)) {
      lo = mid; flo = fm;
    } else {
      hi = mid;
    }
  }
  const Real xl1 = 0.5*(lo + hi);
  out.dl1 = xl1 - xa;
  Real s[4] = {xl1 - 1.0e-4, 0.0, -eps, 0.0};
  const Real h = 2.0e-6;
  Real rprev = std::hypot(s[0] - xa, s[1]);
  Real sp[4];
  bool falling = true;
  for (int n=0; n<20000000; ++n) {
    for (int q=0; q<4; ++q) sp[q] = s[q];
    Real k1[4], k2[4], k3[4], k4[4], t[4];
    BallisticRhs(mu, s, k1);
    for (int q=0; q<4; ++q) t[q] = s[q] + 0.5*h*k1[q];
    BallisticRhs(mu, t, k2);
    for (int q=0; q<4; ++q) t[q] = s[q] + 0.5*h*k2[q];
    BallisticRhs(mu, t, k3);
    for (int q=0; q<4; ++q) t[q] = s[q] + h*k3[q];
    BallisticRhs(mu, t, k4);
    for (int q=0; q<4; ++q) s[q] += h*(k1[q] + 2.0*k2[q] + 2.0*k3[q] + k4[q])/6.0;
    const Real r = std::hypot(s[0] - xa, s[1]);
    if (!out.ok_cross && rprev > r_cross && r <= r_cross) {
      const Real w = (rprev - r_cross)/(rprev - r);
      Real si[4];
      for (int q=0; q<4; ++q) si[q] = sp[q] + w*(s[q] - sp[q]);
      const Real dx = si[0] - xa, dy = si[1];
      const Real rr = std::hypot(dx, dy);
      const Real rx = dx/rr, ry = dy/rr;
      out.phi = std::atan2(dy, dx);
      out.vr = si[2]*rx + si[3]*ry;
      out.vt = -si[2]*ry + si[3]*rx;
      out.ok_cross = true;
    }
    if (falling && r > rprev) {   // first periastron
      out.rmin = rprev;
      out.phimin = std::atan2(sp[1], sp[0] - xa);
      falling = false;
      break;
    }
    rprev = r;
    if (r < 1.0e-3) break;
  }
  return out;
}

//----------------------------------------------------------------------------------------
//! the stream inflow / outflow / absorbing-surface BCs on the x1 faces
void RyPerBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x3v = pmbp->pcoord->x3v;
  auto u0 = pmbp->phydro->u0;
  const Real cs2 = SQR(pmbp->phydro->peos->eos_data.iso_cs);
  const Real dfl = pmbp->phydro->peos->eos_data.dfloor;
  const RocheParams p = rp_;
  const Real vwall = (spin_ - 1.0)*rp_.omega;
  const bool noslip = noslip_, ihse = inner_hse_, ohse = outer_hse_, strm = stream_on_;
  const bool iwall = inner_wall_;
  const Real rhos = rho_s_, phis = phi_s_, sig = sig_, wwin = nsig_*sig_;
  const Real vrs = vr_s_, vps = vp_s_, cap = hse_cap_;
  par_for("ryper_bc", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real ph = x3v(m,k);
    // ---- inner: absorbing rigid surface
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      const Real da = u0(m,IDN,k,j,is);
      const Real v1a = u0(m,IM1,k,j,is)/da;
      const Real v2a = u0(m,IM2,k,j,is)/da;
      const Real v3a = u0(m,IM3,k,j,is)/da;
      const Real pa = RochePot(p, x1v(m,is), ph);
      for (int g=0; g<ng; ++g) {
        const int ig = is - 1 - g;
        Real dg = da;
        if (ihse) {
          const Real ex = -(RochePot(p, x1v(m,ig), ph) - pa)/cs2;
          dg = da*exp(fmin(ex, cap));
        }
        dg = fmax(dg, dfl);
        // diode (absorbing), or a closed wall mirroring v_r (test only)
        const Real v1 = iwall ? -u0(m,IM1,k,j,is+g)/u0(m,IDN,k,j,is+g) : fmin(v1a, 0.0);
        const Real v3 = noslip ? vwall*x1v(m,ig) : v3a;
        u0(m,IDN,k,j,ig) = dg;
        u0(m,IM1,k,j,ig) = dg*v1;
        u0(m,IM2,k,j,ig) = dg*v2a;
        u0(m,IM3,k,j,ig) = dg*v3;
      }
    }
    // ---- outer: stream window, else outflow without inflow
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      const Real dph = WrapPhi(ph - phis);
      const bool win = strm && (fabs(dph) < wwin);
      const Real da = u0(m,IDN,k,j,ie);
      const Real v1a = u0(m,IM1,k,j,ie)/da;
      const Real v2a = u0(m,IM2,k,j,ie)/da;
      const Real v3a = u0(m,IM3,k,j,ie)/da;
      const Real pa = RochePot(p, x1v(m,ie), ph);
      for (int g=0; g<ng; ++g) {
        const int ig = ie + 1 + g;
        Real dg, v1, v2, v3;
        if (win) {
          dg = rhos*exp(-0.5*SQR(dph/sig));
          v1 = vrs; v2 = 0.0; v3 = vps;
        } else {
          dg = da;
          if (ohse) {
            const Real ex = -(RochePot(p, x1v(m,ig), ph) - pa)/cs2;
            dg = da*exp(fmin(ex, cap));
          }
          v1 = fmax(v1a, 0.0); v2 = v2a; v3 = v3a;
        }
        dg = fmax(dg, dfl);
        u0(m,IDN,k,j,ig) = dg;
        u0(m,IM1,k,j,ig) = dg*v1;
        u0(m,IM2,k,j,ig) = dg*v2;
        u0(m,IM3,k,j,ig) = dg*v3;
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! Roche gravity + Coriolis (in-plane), and the RK-weighted boundary flux integrals
void RyPerSrc(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto u0 = pmbp->phydro->u0;
  auto w0 = pmbp->phydro->w0;
  auto gr = gr_, gp = gp_;
  const Real om2 = 2.0*rp_.omega;
  par_for("ryper_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    const Real vr = w0(m,IVX,k,j,i), vp = w0(m,IVZ,k,j,i);
    u0(m,IM1,k,j,i) += bdt*d*(gr(m,k,i) + om2*vp);
    u0(m,IM3,k,j,i) += bdt*d*(gp(m,k,i) - om2*vr);
  });

  if (nstages_ <= 0) return;
  // boundary flux rates (into the domain counted positive for Min, into the star for
  // Macc/Jacc/Jstr, outward for Mout/Jout)
  auto &flx = pmbp->phydro->uflx.x1f;
  auto &area = pmbp->pcoord->area.x1f;
  auto &xf = pmbp->pcoord->xx1f;
  auto &x3v = pmbp->pcoord->x3v;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  const Real om = rp_.omega, phis = phi_s_, wwin = nsig_*sig_;
  const bool strm = stream_on_;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj;
  array_sum::GlobalSum fs;
  Kokkos::parallel_reduce("ryper_flx", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/(nk*nj);
    const int k = (idx - m*nk*nj)/nj + ks;
    const int j = idx%nj + js;
    array_sum::GlobalSum v;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      const Real a = area(m,k,j,is), r = xf(m,is);
      const Real fd = flx(m,IDN,k,j,is), f3 = flx(m,IM3,k,j,is);
      const Real v3a = w0(m,IVZ,k,j,is);
      v.the_array[2] = -fd*a;
      v.the_array[3] = -(f3 + om*r*fd)*r*a;
      v.the_array[4] = -(f3 - fd*v3a)*r*a;
    }
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      const Real a = area(m,k,j,ie+1), r = xf(m,ie+1);
      const Real fd = flx(m,IDN,k,j,ie+1), f3 = flx(m,IM3,k,j,ie+1);
      const bool win = strm && (fabs(WrapPhi(x3v(m,k) - phis)) < wwin);
      if (win) {
        v.the_array[0] = -fd*a;
      } else {
        v.the_array[1] = fd*a;
      }
      v.the_array[5] = (f3 + om*r*fd)*r*a;
    }
    sum += v;
  }, Kokkos::Sum<array_sum::GlobalSum>(fs));

  const int s = stage_ctr_ % nstages_;
  if (s == 0) {
    for (int n=0; n<kNacc; ++n) acc1_[n] = acc0_[n];
  }
  for (int n=0; n<kNacc; ++n) {
    acc0_[n] = rk_g0_[s]*acc0_[n] + rk_g1_[s]*acc1_[n] + bdt*fs.the_array[n];
  }
  stage_ctr_++;
}

//----------------------------------------------------------------------------------------
void RyPerHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = pmbp->phydro->u0;
  auto &vol = pmbp->pcoord->volume;
  auto &x1v = pmbp->pcoord->x1v;
  const Real om = rp_.omega;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  array_sum::GlobalSum s;
  Kokkos::parallel_reduce("ryper_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx%nx1 + is;
    array_sum::GlobalSum v;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
    const Real dv = vol(m,k,j,i), r = x1v(m,i), d = u0(m,IDN,k,j,i);
    v.the_array[0] = dv*d;
    v.the_array[1] = dv*r*(u0(m,IM3,k,j,i) + om*r*d);
    sum += v;
  }, Kokkos::Sum<array_sum::GlobalSum>(s));

  const char *lab[14] = {"Mdom", "Jdom", "Min", "Mout", "Macc", "Jacc", "Jstr", "Jout",
                         "dMin", "dMout", "dMacc", "dJacc", "dJstr", "dJout"};
  pdata->nhist = 14;
  for (int n=0; n<14; ++n) pdata->label[n] = lab[n];
  pdata->hdata[0] = s.the_array[0];
  pdata->hdata[1] = s.the_array[1];
  const Real t = pm->time;
  const Real dt = (hist_tprev_ >= 0.0) ? (t - hist_tprev_) : 0.0;
  for (int n=0; n<kNacc; ++n) {
    pdata->hdata[2+n] = acc0_[n];
    pdata->hdata[8+n] = (dt > 0.0) ? (acc0_[n] - hist_prev_[n])/dt : 0.0;
    hist_prev_[n] = acc0_[n];
  }
  hist_tprev_ = t;
}

//! release the namespace-scope Views before Kokkos::finalize
void RyPerFinal(ParameterInput *pin, Mesh *pm) {
  gr_ = DvceArray3D<Real>();
  gp_ = DvceArray3D<Real>();
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  const bool root = (global_variable::my_rank == 0);
  auto fatal = [&](const std::string &msg) {
    std::cout << "### FATAL ERROR in ry_per_accretor: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
  };
  if (pmbp->phydro == nullptr || pmbp->pmhd != nullptr) fatal("needs <hydro>, no <mhd>");
  if (!pmy_mesh_->use_spherical_polar) fatal("needs mesh/use_spherical_polar = true");
  std::string thermo = pin->GetOrAddString("problem", "thermo", "isothermal");
  if (thermo.compare("isothermal") != 0) {
    fatal("problem/thermo = '" + thermo + "': only isothermal is implemented (adiabatic "
          "is reserved; see the file header for what it needs)");
  }
  auto &eos = pmbp->phydro->peos->eos_data;
  if (eos.is_ideal) fatal("needs <hydro>/eos = isothermal");
  if (!user_srcs) fatal("set problem/user_srcs = true");
  if (!user_bcs) fatal("set mesh/ix1_bc = ox1_bc = user");

  // ---- binary (Barai+04 spectroscopic defaults; physical units)
  const Real ma = pin->GetOrAddReal("problem", "m_acc", 6.24);        // Msun
  const Real md = pin->GetOrAddReal("problem", "m_don", 1.69);        // Msun
  const Real asep = pin->GetOrAddReal("problem", "a_sep", 30.3);      // Rsun
  const Real pday = pin->GetOrAddReal("problem", "period", 6.863569); // d (reported only)
  const Real racc = pin->GetOrAddReal("problem", "r_acc", 4.06);      // Rsun
  const Real tdon = pin->GetOrAddReal("problem", "t_don", 6250.0);    // K
  const Real mudon = pin->GetOrAddReal("problem", "mu_don", 1.27);
  spin_ = pin->GetOrAddReal("problem", "spin", 1.0);   // Omega_* / Omega_orb
  std::string slip = pin->GetOrAddString("problem", "inner_slip", "noslip");
  noslip_ = (slip.compare("noslip") == 0);
  if (!noslip_ && slip.compare("free") != 0) fatal("inner_slip must be noslip | free");
  std::string irho = pin->GetOrAddString("problem", "inner_rho", "copy");
  std::string orho = pin->GetOrAddString("problem", "outer_rho", "hse");
  inner_hse_ = (irho.compare("hse") == 0);
  std::string ivr = pin->GetOrAddString("problem", "inner_vr", "diode");
  inner_wall_ = (ivr.compare("wall") == 0);
  if (!inner_wall_ && ivr.compare("diode") != 0) fatal("inner_vr must be diode | wall");
  outer_hse_ = (orho.compare("hse") == 0);
  hse_cap_ = pin->GetOrAddReal("problem", "hse_cap", 40.0);
  stream_on_ = pin->GetOrAddBoolean("problem", "stream", true);
  rho_s_ = pin->GetOrAddReal("problem", "rho_stream", 1.0);
  const Real rho_amb = pin->GetOrAddReal("problem", "rho_amb", 1.0e-5);
  nsig_ = pin->GetOrAddReal("problem", "stream_nsig", 3.0);

  const Real gu = kG*kMsun/(kRsun*kVel*kVel);   // G Msun in (km/s)^2 Rsun
  rp_.gma = gu*ma;
  rp_.gmd = gu*md;
  rp_.asep = asep;
  rp_.omega = std::sqrt(gu*(ma + md)/(asep*asep*asep));
  rp_.xcm = asep*md/(ma + md);
  const Real porb = 2.0*M_PI/rp_.omega;         // code time
  const Real r_in = pmy_mesh_->mesh_size.x1min, r_out = pmy_mesh_->mesh_size.x1max;

  // ---- ballistic stream
  const Real mu = md/(ma + md);
  const Real aom = asep*rp_.omega;               // km/s
  const Real csd = std::sqrt(kKB*tdon/(mudon*kMH))/kVel;   // km/s
  BallisticOut bo = IntegrateStream(mu, csd/aom, r_out/asep);
  BallisticOut bi = IntegrateStream(mu, csd/aom, racc/asep);
  const Real width = pin->GetOrAddReal("problem", "stream_width", csd/rp_.omega);
  if (!bo.ok_cross) fatal("the ballistic stream never reaches r_out = x1max");
  phi_s_ = bo.phi;
  vr_s_ = bo.vr*aom;
  vp_s_ = bo.vt*aom;
  if (pin->DoesParameterExist("problem", "stream_phi_deg")) {
    phi_s_ = pin->GetReal("problem", "stream_phi_deg")*M_PI/180.0;
  }
  if (pin->DoesParameterExist("problem", "stream_vr")) {
    vr_s_ = pin->GetReal("problem", "stream_vr");
  }
  if (pin->DoesParameterExist("problem", "stream_vphi")) {
    vp_s_ = pin->GetReal("problem", "stream_vphi");
  }
  sig_ = width/r_out;

  // ---- integrator weights for the boundary flux integrals
  std::string integ = pin->GetString("time", "integrator");
  nstages_ = 0;
  if (integ.compare("rk1") == 0) {
    nstages_ = 1; rk_g0_[0] = 0.0; rk_g1_[0] = 1.0;
  } else if (integ.compare("rk2") == 0) {
    nstages_ = 2; rk_g0_[0] = 0.0; rk_g1_[0] = 1.0; rk_g0_[1] = 0.5; rk_g1_[1] = 0.5;
  } else if (integ.compare("rk3") == 0) {
    nstages_ = 3; rk_g0_[0] = 0.0; rk_g1_[0] = 1.0; rk_g0_[1] = 0.25; rk_g1_[1] = 0.75;
    rk_g0_[2] = 2.0/3.0; rk_g1_[2] = 1.0/3.0;
  } else if (root) {
    std::cout << "ry_per_accretor: WARNING integrator " << integ << " not rk1/2/3: the "
              << "boundary flux integrals (Min..Jout) are switched off" << std::endl;
  }
  stage_ctr_ = 0;
  for (int n=0; n<kNacc; ++n) {
    acc0_[n] = acc1_[n] = hist_prev_[n] = 0.0;
  }
  hist_tprev_ = -1.0;

  if (root) {
    std::printf("ry_per_accretor: M_a %.4g M_d %.4g Msun, a %.4g Rsun, q %.4f, Omega %.6e"
                "/s (Kepler; 2pi/P_input %.6e), P_orb %.6f code = %.6f d\n", ma, md, asep,
                md/ma, rp_.omega/kTime, 2.0*M_PI/(pday*86400.0), porb,
                porb*kTime/86400.0);
    std::printf("ry_per_accretor: x_cm %.4f Rsun, d_L1 %.4f Rsun, r_in %.4f r_out %.4f "
                "Rsun = %.4f d_L1, R_acc key %.4f\n", rp_.xcm, bo.dl1*asep, r_in, r_out,
                r_out/(bo.dl1*asep), racc);
    std::printf("ry_per_accretor: c_s,don %.3f km/s, a Omega %.2f km/s, width %.4f Rsun, "
                "sigma_phi %.5f rad\n", csd, aom, width, sig_);
    std::printf("ry_per_accretor: ballistic at r_out: phi %.3f deg, v_r %.2f, v_phi %.2f "
                "km/s; stream state used: phi %.3f deg, v_r %.2f, v_phi %.2f km/s\n",
                bo.phi*180.0/M_PI, bo.vr*aom, bo.vt*aom, phi_s_*180.0/M_PI, vr_s_, vp_s_);
    std::printf("ry_per_accretor: ballistic r_min %.4f a = %.4f Rsun at phi %.2f deg; at "
                "R_acc: phi %.3f deg, v_r %.2f, v_phi %.2f km/s (angle from inward normal"
                " %.2f deg)\n", bo.rmin, bo.rmin*asep, bo.phimin*180.0/M_PI,
                bi.phi*180.0/M_PI, bi.vr*aom, bi.vt*aom,
                std::atan2(std::fabs(bi.vt), -bi.vr)*180.0/M_PI);
    std::printf("ry_per_accretor: spin %.3f (wall v_phi %.3f km/s at r_in, rotating "
                "frame), inner %s/%s/%s, outer %s, stream %s, iso_cs %.3f km/s\n", spin_,
                (spin_ - 1.0)*rp_.omega*r_in, slip.c_str(), irho.c_str(), ivr.c_str(),
                orho.c_str(), stream_on_ ? "on" : "off", eos.iso_cs);
    if (std::fabs(r_in - racc) > 1.0e-6*racc) {
      std::printf("ry_per_accretor: NOTE mesh x1min %.6f != problem/r_acc %.6f (the "
                  "absorbing surface sits at x1min)\n", r_in, racc);
    }
  }

  // ---- precomputed gravity (discrete gradients of Phi across the cell faces)
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmbp->nmb_thispack;
  Kokkos::realloc(gr_, nmb, n3, n1);
  Kokkos::realloc(gp_, nmb, n3, n1);
  {
    auto gr = gr_, gp = gp_;
    auto &x1v = pmbp->pcoord->x1v;
    auto &x1f = pmbp->pcoord->xx1f;
    auto &x3v = pmbp->pcoord->x3v;
    auto &x3f = pmbp->pcoord->xx3f;
    const RocheParams p = rp_;
    par_for("ryper_g", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      const Real r = x1v(m,i), rl = x1f(m,i), rr = x1f(m,i+1);
      const Real ph = x3v(m,k), pl = x3f(m,k), pr = x3f(m,k+1);
      gr(m,k,i) = -(RochePot(p, rr, ph) - RochePot(p, rl, ph))/(rr - rl);
      gp(m,k,i) = -(RochePot(p, r, pr) - RochePot(p, r, pl))/(r*(pr - pl));
    });
  }

  user_bcs_func = RyPerBC;
  user_srcs_func = RyPerSrc;
  if (user_hist) user_hist_func = RyPerHist;
  pgen_final_func = RyPerFinal;
  if (restart) return;

  // ---- initial state: at rest in the rotating frame
  std::string init = pin->GetOrAddString("problem", "init", "ambient");
  const bool ihse = (init.compare("hse") == 0);
  if (!ihse && init.compare("ambient") != 0) fatal("problem/init must be ambient | hse");
  const Real cs2 = SQR(eos.iso_cs), cap = hse_cap_, dfl = eos.dfloor;
  const RocheParams p = rp_;
  const Real phiref = RochePot(p, r_out, M_PI);
  auto u0 = pmbp->phydro->u0;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x3v = pmbp->pcoord->x3v;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  par_for("ryper_ic", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real d = rho_amb;
    if (ihse) {
      const Real ex = -(RochePot(p, x1v(m,i), x3v(m,k)) - phiref)/cs2;
      d = rho_amb*exp(fmin(ex, cap));
    }
    u0(m,IDN,k,j,i) = fmax(d, dfl);
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
  });
}
