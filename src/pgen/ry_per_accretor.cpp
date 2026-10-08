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
//! THERMO = GENERAL (problem/thermo = general, envelope stage only; accretor-rhd-1008,
//! design docs/dev/accretor_rhd_design.md stages S1/S2).  <hydro>/eos = general (table,
//! eos_radiation = false), <units> Rsun / km/s / a fixed density unit (mu 1), wb_option
//! = polytropic.  Every ideal-gas site of the envelope stage goes through the EOS: ghost
//! T and the isothermal ghost walk (WBAdvance), the WB pressure, the relaxation targets
//! (kelvin keys env_t_ph / env_t_amb / env_t_stream; envelope and atmosphere relax to
//! the column T(psi)), the ambient's P/rho, and the IC: T(psi), rho(psi) from the S0
//! column file problem/env_ic_file, balanced discretely (EnvICGeneral).  With thermo =
//! adiabatic (ideal gas) every one of these sites runs its original code verbatim.
//!
//! RADIATION (<rad_m1> present; stage S3 of docs/dev/accretor_rhd_design.md; needs
//! inner = envelope, thermo = general, env_ic = column).  Implicit grey M1 (closure
//! vet_col: vet_gd cannot run on the 4-cell theta band) with the run's Rosseland + Planck
//! tables (problem/env_opac_table, env_planck_table; the S0 column must be built from the
//! same Rosseland table and the same EOS).  The column's radiation force enters as a
//! REFERENCE (force_reference = wb_arad): Phi_eff = Phi_wb + G(psi), G = -int a_ref dr_eq
//! along the column, a_ref = kappa_t F/c (RadSetup), so the column balances in the
//! discrete WB sense with gas pressure only; the module applies the residual
//! kappa_t F/c - a_ref and this pgen adds the reference work (split).  Bottom: the
//! column's flux F(r_in) through the implicit face BC (implicit_bc_x1min = flux,
//! implicit_flux_x1min checked against the column to 1e-3); top: Marshak/vacuum
//! (implicit_bc_x1max = marshak, dark ghosts).  IC: E = a T_col^4, F_r = F_col(psi).
//! The hot ambient keeps its T relaxation and gets no absorption: <rad_m1>/
//! opac_abs_rho_max (module key) removes kappa_P, kappa_E below that density and keeps
//! electron scattering.  The relaxation/sponge energy is booked (history Erel).  The
//! relaxation inside the envelope and of the dense atmosphere should be off
//! (env_t_relax_env = env_t_relax_stream = 0).  History: the nine rate columns become
//! Lmeas Ltop Erad Etot KE Ebnd Erel Pic Picmax (RyPerHistEnv).
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "utils/wb_background.hpp"
#include "outputs/outputs.hpp"
#include "units/units.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
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
// stage 2 (problem/inner = envelope): the resolved stellar envelope, see the header
bool env_ = false, env_noslip_ = false;
Real env_np_ = 3.0, env_rhoph_ = 20.0, env_cph2_ = 240.25, env_rtop_ = 4.11;
Real env_tro_ = 1.0e-4, env_tri_ = 1.0e-4, env_trs_ = 1.0e-4, env_phis_ = 0.0, env_racc_ = 4.06;
Real env_dlt_ = 0.0;           // (spin^2 - 1) Omega^2 / 2
// problem/env_cool (default false): the envelope and stream/atmosphere relaxation times become
// the local radiative cooling time (diffusion + optically thin bridge, ES + Kramers opacity);
// the hot ambient keeps env_t_relax.  Constants in cgs, set in EnvSetup.
bool env_cool_ = false;
Real cool_rhou_ = 5.11e-8;      // g/cc per code density (problem/cool_rho_unit)
Real cool_kes_ = 0.34;          // cm^2/g electron scattering (problem/cool_kappa_es)
Real cool_kk0_ = 1.5e24;        // Kramers kappa = k0 rho T^-3.5 (problem/cool_kappa_k0)
Real cool_tfac_ = 0.0;          // K per (km/s)^2 of P/rho: mu m_H kVel^2 / k_B (problem/cool_mu)
Real cool_tmin_ = 0.0;          // floor of the cooling time [code time] (problem/cool_t_min)
Real env_rhoamb_ = 1.0e-6;     // floor (ambient) density = problem/rho_amb
Real env_camb2_ = 9.0e4;       // hot hydrostatic ambient P/rho (problem/env_cs_amb^2)
Real env_ramb_ = 1.0e-6;       // ambient density at r = R_acc (problem/env_amb_rho)
Real env_cstr2_ = 240.25;      // stream (far dense gas) P/rho, problem/env_cs_stream^2
Real env_rspin_ = 1.0e30;      // spin term / rotation ends at this radius (env_r_spin)
Real env_rhot_ = 1.0e30;       // dense gas outside the star relaxes to c_ph^2 below this
                               // radius, to env_cs_stream^2 above (problem/env_r_hot)
Real env_kamb_ = 3.0;          // gas denser than env_amb_k x the hydrostatic ambient is
                               // stream/star material and relaxes to c_ph^2
DvceArray3D<Real> rhoamb_;     // the hydrostatic ambient density at (m,k,i)
Real env_tsp_ = 1.0e-4;        // floor velocity damping time (problem/env_t_sponge)
Real env_fsp_ = 10.0;          // sponge acts where rho < env_sponge_rho x rho_amb
Real env_vcap_ = 100.0;        // |v| cap of sponge cells (problem/env_v_floor_max)
DvceArray3D<Real> cref2_;      // relaxation target P/rho at (m,k,i)
DvceArray3D<Real> gpl_;        // plain -dPhi_wb/dr (face difference), above hydro/wb_rmax
DvceArray1D<int> iracc_;       // per MeshBlock: index of the measuring x1 face (r_meas), or -1
// problem/env_top_mode = equipotential: Phi_wb is capped at env_phtop_ (the equipotential
// psi = -env_top_hp c_ph^2) instead of being flat above a spherical r_top; +huge = off
Real env_phtop_ = 1.0e300;
Real env_rmeas_ = 4.06;        // measuring face radius (problem/r_meas, default R_acc)
// problem/thermo = general (envelope only; default off = the ideal-gas code paths
// verbatim): <hydro>/eos = general (tabulated, gas only), temperatures in kelvin keys,
// the envelope T(psi) from a 1-D column file (problem/env_ic_file, S0 script
// docs/dev/accretor_rhd/make_ic_accretor_column.py), the discrete balance of the x1 WB
// pair (wb_option = polytropic) for the IC.  All T below are CODE temperatures of the
// general EOS (kelvin / eos.temp_cgs, mu_ref = 1).
bool gen_ = false;
Real gen_tph_ = 0.0;           // photospheric T (problem/env_t_ph, K)
Real gen_tamb_ = 0.0;          // hot ambient T (problem/env_t_amb, K)
Real gen_tstr_ = 0.0;          // stream T at the window (problem/env_t_stream, K)
Real gen_tfl_ = 0.0;           // ghost T floor (eos tfloor, or 1e-3 T_ph if unset)
DvceArray3D<Real> tcol_;       // column T(psi) at (m,k,i) [code T]
// flux accumulators (per rank), RK registers, stage counter.  The absorbing surface
// (stage 1) uses the first kNaccS, the envelope all kNacc.
constexpr int kNaccS = 6;
constexpr int kNacc = 9;
int nstages_ = 0;
Real rk_g0_[3], rk_g1_[3];
int stage_ctr_ = 0;
// <rad_m1> on the envelope (S3) adds four registers: radiation energy in at r_in (the
// imposed flux), radiation out at r_out (the comoving face flux f0x1), the hydro energy
// flux out at r_out (etotgrav: incl. rho Phi v), and the energy put in by the T
// relaxation and taken out by the floor sponge
constexpr int kNaccR = 13;
Real acc0_[kNaccR] = {0.0}, acc1_[kNaccR] = {0.0};
Real hist_prev_[kNaccR] = {0.0};
Real hist_tprev_ = -1.0;
// <rad_m1> on the envelope (stage S3, accretor-rhd-1008; see the header RADIATION)
bool rad_ = false;
Real rad_fin_ = 0.0;           // imposed x1min radiative flux (code), = F_col(r_in)
Real rad_lfac_ = 1.0;          // code band luminosity -> Lsun of the full sphere
Real rad_pn0_ = 0.0, rad_ps0_ = 0.0;   // Picard counters at the previous history output
std::vector<Real> rcol_psi_, rcol_g_, rcol_f_;   // column psi, G(psi), F(psi) [code]

KOKKOS_INLINE_FUNCTION
Real RochePot(const RocheParams &p, const Real r, const Real phi) {
  const Real cp = cos(phi);
  const Real rd = sqrt(fmax(r*r + p.asep*p.asep - 2.0*p.asep*r*cp, 1.0e-30));
  return -p.gma/r - p.gmd/rd
         - 0.5*SQR(p.omega)*(r*r + p.xcm*p.xcm - 2.0*p.xcm*r*cp);
}

//! stage 2: the potential of the x1 well-balanced pair.  Below r_top: Phi - Delta,
//! Delta = (s^2 - 1) Omega^2 r^2/2 (the envelope's effective potential incl. its own
//! spin).  Above r_top, where the gas is
//! unsupported ambient/stream, it is FLAT in r (the WB background then equals the cell
//! state: plain PLM and no WB force), and gravity there is the explicit -dD/dr,
//! D = Phi - Phi_wb (gr_).  The spin term stops growing at r_spin <= r_top
//! (problem/env_r_spin, default r_top): above it the atmosphere is at rest (synchronous).
KOKKOS_INLINE_FUNCTION
Real PhiWB(const RocheParams &p, const Real dlt, const Real rtop, const Real rspin,
           const Real phtop, const Real r, const Real phi) {
  const Real x = fmin(r, rtop), xs = fmin(x, rspin);
  return fmin(RochePot(p, x, phi) - dlt*xs*xs, phtop);
}

//! stage 2: envelope P/rho at depth psi = Phi_s - Phi_wb (n-polytrope shifted to the
//! photospheric P/rho c_ph^2 at psi = 0; isothermal c_ph^2 above)
KOKKOS_INLINE_FUNCTION
Real EnvC2(const Real psi, const Real cph2, const Real np) {
  return (psi > 0.0) ? cph2 + psi/(np + 1.0) : cph2;
}

KOKKOS_INLINE_FUNCTION
Real EnvRho(const Real psi, const Real cph2, const Real np, const Real rhoph) {
  return (psi > 0.0) ? rhoph*pow(1.0 + psi/((np + 1.0)*cph2), np)
                     : rhoph*exp(fmax(psi/cph2, -700.0));
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
    for (int n=0; n<kNaccS; ++n) acc1_[n] = acc0_[n];
  }
  for (int n=0; n<kNaccS; ++n) {
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
  for (int n=0; n<kNaccS; ++n) {
    pdata->hdata[2+n] = acc0_[n];
    pdata->hdata[8+n] = (dt > 0.0) ? (acc0_[n] - hist_prev_[n])/dt : 0.0;
    hist_prev_[n] = acc0_[n];
  }
  hist_tprev_ = t;
}

//----------------------------------------------------------------------------------------
//! stage 2 (problem/inner = envelope): closed wall at r_in under the envelope, stream
//! window / outflow at r_out.  Ideal gas with <hydro>/etotgrav: the conserved energy
//! carries rho Phi (phicc0), so the ghosts get it too.  Ghost density is the isothermal
//! hydrostatic extrapolation in Phi_wb at the edge cell's P/rho (the x1 WB background of
//! that cell, so the wall face sees equal states); velocities are mirrored (v_r) and
//! free-slip or held at the initial spin (v_phi, problem/env_wall_slip).
void RyPerBCEnv(Mesh *pm) {
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
  auto phicc = pmbp->phydro->phicc0;
  const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  const Real dfl = pmbp->phydro->peos->eos_data.dfloor;
  const RocheParams p = rp_;
  const Real vwall = (spin_ - 1.0)*rp_.omega;
  const bool noslip = env_noslip_, strm = stream_on_;
  const Real rhos = rho_s_, phis = phi_s_, sig = sig_, wwin = nsig_*sig_;
  const Real vrs = vr_s_, vps = vp_s_, cap = hse_cap_;
  const Real dlt = env_dlt_, rtop = env_rtop_, cph2 = env_cph2_, phtop = env_phtop_;
  const Real cstr2 = env_cstr2_, rspin = env_rspin_;
  const bool gen = gen_;
  const EOS_Data eosd = pmbp->phydro->peos->eos_data;
  const Real tfl = gen_tfl_, tstr = gen_tstr_;
  // <rad_m1> (S3): the inner ghosts walk along the EFFECTIVE potential of the x1 WB pair
  // (Phi_wb + G, the radiation-force reference, phicc_wb), and the radiation ghosts are
  // filled as he_star_m1's: r_in copy (the implicit solve imposes the flux there), r_out
  // vacuum (dark; Marshak through the implicit face BC)
  const bool rad = rad_;
  auto pwcc = pmbp->phydro->phicc_wb;
  DvceArray5D<Real> ur;
  Real rcl = 1.0, refl = 0.0;
  if (rad) {
    ur = pmbp->pradm1->u0;
    rcl = pmbp->pradm1->c_light;
    refl = pmbp->pradm1->e_floor;
  }
  par_for("ryper_bce", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real ph = x3v(m,k);
    if (gen) {
      // thermo = general: the same BCs through the EOS.  r_in: the local polytrope of
      // the two edge cells in Phi_wb (WBAdvance polytropic branch); r_out: ghost T = the
      // edge cell's (floored), density by the isothermal hydrostatic walk of the general
      // EOS in the TRUE Roche Phi; IEN = e + KE + rho Phi.
      if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        const Real da = u0(m,IDN,k,j,is);
        const Real kea = 0.5*(SQR(u0(m,IM1,k,j,is)) + SQR(u0(m,IM2,k,j,is))
                              + SQR(u0(m,IM3,k,j,is)))/da;
        Real ea = u0(m,IEN,k,j,is) - kea - da*phicc(m,k,j,is);
        Real ta = (ea > 0.0) ? eosd.Temperature(da, ea) : tfl;
        if (!(ta > tfl)) {
          ta = tfl;
          ea = eosd.EnergyFromTemperature(da, ta);
        }
        const Real pa = rad ? pwcc(m,k,j,is) : PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,is), ph);
        // the local polytrope of the two edge cells (wb_option = polytropic): T linear in
        // Phi_wb with their dT/dPhi, density by the WB polytropic walk, so the ghosts lie
        // on the edge cell's own WB background
        Real a = 0.0;
        {
          const int i1 = is + 1;
          const Real d1 = u0(m,IDN,k,j,i1);
          const Real e1 = u0(m,IEN,k,j,i1) - 0.5*(SQR(u0(m,IM1,k,j,i1))
                          + SQR(u0(m,IM2,k,j,i1)) + SQR(u0(m,IM3,k,j,i1)))/d1
                          - d1*phicc(m,k,j,i1);
          const Real dp1 = (rad ? pwcc(m,k,j,i1) :
                            PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,i1), ph)) - pa;
          if (e1 > 0.0 && dp1 != 0.0) a = (eosd.Temperature(d1, e1, ta) - ta)/dp1;
        }
        for (int g=0; g<ng; ++g) {
          const int ig = is - 1 - g, im = is + g;
          const Real pg = rad ? pwcc(m,k,j,ig) :
                          PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,ig), ph);
          Real dg = da, eg = ea, tg = ta;
          WBAdvance(eosd, 3, da, ea, pg - pa, dg, eg, tg, ta, a, ta, ta);
          if (!(Kokkos::isfinite(dg) && dg > 0.0 && Kokkos::isfinite(eg) && eg > 0.0 &&
                tg >= tfl)) {
            dg = da; eg = ea; tg = ta;
            WBAdvance(eosd, 1, da, ea, pg - pa, dg, eg, tg, ta, 0.0, ta, ta);
            if (!(Kokkos::isfinite(dg) && dg > 0.0)) dg = da;
            tg = ta;
          }
          dg = fmax(fmin(dg, da*exp(cap)), dfl);
          eg = eosd.EnergyFromTemperature(dg, tg);
          const Real dm = u0(m,IDN,k,j,im);
          const Real v1 = -u0(m,IM1,k,j,im)/dm;
          const Real v2 = u0(m,IM2,k,j,im)/dm;
          const Real v3 = noslip ? vwall*x1v(m,ig) : u0(m,IM3,k,j,im)/dm;
          u0(m,IDN,k,j,ig) = dg;
          u0(m,IM1,k,j,ig) = dg*v1;
          u0(m,IM2,k,j,ig) = dg*v2;
          u0(m,IM3,k,j,ig) = dg*v3;
          u0(m,IEN,k,j,ig) = eg + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
          if (rad) radm1::M1FillGhost(ur, m, k, j, ig, k, j, is, 1, 0, -1.0, rcl, refl);
        }
      }
      if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
        const Real dph = WrapPhi(ph - phis);
        const bool win = strm && (fabs(dph) < wwin);
        const Real da = u0(m,IDN,k,j,ie);
        const Real v1a = u0(m,IM1,k,j,ie)/da;
        const Real v2a = u0(m,IM2,k,j,ie)/da;
        const Real v3a = u0(m,IM3,k,j,ie)/da;
        Real ea = u0(m,IEN,k,j,ie) - 0.5*da*(v1a*v1a + v2a*v2a + v3a*v3a)
                  - da*phicc(m,k,j,ie);
        Real ta = (ea > 0.0) ? eosd.Temperature(da, ea) : tfl;
        if (!(ta > tfl)) {
          ta = tfl;
          ea = eosd.EnergyFromTemperature(da, ta);
        }
        const Real pa = RochePot(p, x1v(m,ie), ph);
        for (int g=0; g<ng; ++g) {
          const int ig = ie + 1 + g;
          Real dg, v1, v2, v3, tg;
          if (win) {
            dg = rhos*exp(-0.5*SQR(dph/sig));
            v1 = vrs; v2 = 0.0; v3 = vps; tg = tstr;
          } else {
            const Real pg = RochePot(p, x1v(m,ig), ph);
            Real eg = ea, tw = ta;
            dg = da;
            WBAdvance(eosd, 1, da, ea, pg - pa, dg, eg, tw, ta, 0.0, ta, ta);
            if (!(Kokkos::isfinite(dg) && dg > 0.0)) dg = da;
            dg = fmin(dg, da*exp(cap));
            v1 = fmax(v1a, 0.0); v2 = v2a; v3 = v3a; tg = ta;
          }
          dg = fmax(dg, dfl);
          u0(m,IDN,k,j,ig) = dg;
          u0(m,IM1,k,j,ig) = dg*v1;
          u0(m,IM2,k,j,ig) = dg*v2;
          u0(m,IM3,k,j,ig) = dg*v3;
          u0(m,IEN,k,j,ig) = eosd.EnergyFromTemperature(dg, tg)
                             + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
          if (rad) radm1::M1FillGhost(ur, m, k, j, ig, k, j, ie, 1, 2, 1.0, rcl, refl);
        }
      }
      return;
    }
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      const Real da = u0(m,IDN,k,j,is);
      const Real kea = 0.5*(SQR(u0(m,IM1,k,j,is)) + SQR(u0(m,IM2,k,j,is))
                            + SQR(u0(m,IM3,k,j,is)))/da;
      const Real ea = u0(m,IEN,k,j,is) - kea - da*phicc(m,k,j,is);
      const Real ta = fmax(gm1*ea/da, 1.0e-3*cph2);
      const Real pa = PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,is), ph);
      for (int g=0; g<ng; ++g) {
        const int ig = is - 1 - g, im = is + g;
        const Real pg = PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,ig), ph);
        const Real dg = fmax(da*exp(fmin(-(pg - pa)/ta, cap)), dfl);
        const Real dm = u0(m,IDN,k,j,im);
        const Real v1 = -u0(m,IM1,k,j,im)/dm;
        const Real v2 = u0(m,IM2,k,j,im)/dm;
        const Real v3 = noslip ? vwall*x1v(m,ig) : u0(m,IM3,k,j,im)/dm;
        u0(m,IDN,k,j,ig) = dg;
        u0(m,IM1,k,j,ig) = dg*v1;
        u0(m,IM2,k,j,ig) = dg*v2;
        u0(m,IM3,k,j,ig) = dg*v3;
        u0(m,IEN,k,j,ig) = dg*ta/gm1 + 0.5*dg*(v1*v1 + v2*v2 + v3*v3)
                           + dg*phicc(m,k,j,ig);
      }
    }
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      const Real dph = WrapPhi(ph - phis);
      const bool win = strm && (fabs(dph) < wwin);
      const Real da = u0(m,IDN,k,j,ie);
      const Real v1a = u0(m,IM1,k,j,ie)/da;
      const Real v2a = u0(m,IM2,k,j,ie)/da;
      const Real v3a = u0(m,IM3,k,j,ie)/da;
      const Real ea = u0(m,IEN,k,j,ie) - 0.5*da*(v1a*v1a + v2a*v2a + v3a*v3a)
                      - da*phicc(m,k,j,ie);
      const Real ta = fmax(gm1*ea/da, 1.0e-3*cph2);
      const Real pa = RochePot(p, x1v(m,ie), ph);
      for (int g=0; g<ng; ++g) {
        const int ig = ie + 1 + g;
        Real dg, v1, v2, v3, tg;
        if (win) {
          dg = rhos*exp(-0.5*SQR(dph/sig));
          v1 = vrs; v2 = 0.0; v3 = vps; tg = cstr2;
        } else {
          const Real pg = RochePot(p, x1v(m,ig), ph);
          dg = da*exp(fmin(-(pg - pa)/ta, cap));
          v1 = fmax(v1a, 0.0); v2 = v2a; v3 = v3a; tg = ta;
        }
        dg = fmax(dg, dfl);
        u0(m,IDN,k,j,ig) = dg;
        u0(m,IM1,k,j,ig) = dg*v1;
        u0(m,IM2,k,j,ig) = dg*v2;
        u0(m,IM3,k,j,ig) = dg*v3;
        u0(m,IEN,k,j,ig) = dg*tg/gm1 + 0.5*dg*(v1*v1 + v2*v2 + v3*v3)
                           + dg*phicc(m,k,j,ig);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! <rad_m1>: the volume sum of u(IEN) over the active cells of this rank
Real RyPerEsum(MeshBlockPack *pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.nx1, nji = indcs.nx2*indcs.nx1, nkji = indcs.nx3*nji;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = pmbp->phydro->u0;
  auto &volume = pmbp->pcoord->volume;
  Real es = 0.0;
  Kokkos::parallel_reduce("ryper_esum", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
    const int i = idx%ni + is;
    sum += u0(m,IEN,k,j,i)*volume(m,k,j,i);
  }, Kokkos::Sum<Real>(es));
  return es;
}

//----------------------------------------------------------------------------------------
//! stage 2 sources: x1 gravity in the well-balanced pressure form (the background of the
//! x1 WB pair, built with Phi_wb) plus -dDelta/dr, phi gravity and Coriolis explicit,
//! thermal relaxation toward P/rho = cref2 (exact exponential over the stage), and the
//! RK-weighted boundary flux integrals incl. the face r = R_acc and the stream AM.
//! The work of gravity is in the etotgrav energy flux (TRUE Phi incl. the orbital
//! centrifugal term); no explicit centrifugal work is added (09-28 double-count guard).
void RyPerSrcEnv(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto *ph = pmbp->phydro;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto wbq0 = ph->wbq0;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto gr = gr_, gp = gp_, c2r = cref2_, gpl = gpl_;
  auto &x1v = pmbp->pcoord->x1v;
  const Real rmax = ph->wb_rmax, phimax = ph->wb_phimax;
  auto pwcs = ph->phicc_wb;
  const Real gm1 = ph->peos->eos_data.gamma - 1.0;
  const Real om2 = 2.0*rp_.omega;
  const Real cph2 = env_cph2_, tro = env_tro_, tri = env_tri_, trs = env_trs_;
  const bool cool = env_cool_;
  const Real rhou = cool_rhou_, kes = cool_kes_, kk0 = cool_kk0_, tfac = cool_tfac_;
  const Real tcmin = cool_tmin_, gma = rp_.gma, lu = kRsun, vu = kVel, tu = kTime;
  const Real rfl = env_rhoamb_, tsp = env_tsp_, fsp = env_fsp_;
  const Real vcap = env_vcap_;
  const Real kamb = env_kamb_, camb2 = env_camb2_, cstr2 = env_cstr2_, rhot = env_rhot_;
  auto ramb = rhoamb_;
  const bool gen = gen_;
  const EOS_Data eosd = ph->peos->eos_data;
  const Real tstr = gen_tstr_, tamb = gen_tamb_;
  auto tcol = tcol_;
  // <rad_m1>: the energy the relaxation and the floor sponge put in/take out is booked
  // (the volume sum of u(IEN) before and after this kernel: no other term of it touches
  // the energy)
  const bool rad = rad_;
  const Real erel0 = rad ? RyPerEsum(pmbp) : 0.0;
  par_for("ryper_srce", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i), e = w0(m,IEN,k,j,i);
    const Real vr = w0(m,IVX,k,j,i), vp = w0(m,IVZ,k,j,i);
    // the WB pressure form where the x1 WB reconstruction is on (r <= wb_rmax), the
    // plain -rho dPhi_wb/dr above it (unsupported ambient: hydro/wb_rmax)
    Real fg;
    if ((rmax > 0.0 && x1v(m,i) > rmax) || pwcs(m,k,j,i) > phimax) {
      fg = d*gpl(m,k,i);
    } else {
      Real d1, pl, d2, pr, d3;
      WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      Real p;
      if (gen) {
        p = (e > 0.0) ? eosd.Pressure(d, e) : 0.5*(pl + pr);
      } else {
        p = (e > 0.0) ? gm1*e : 0.5*(pl + pr);
      }
      fg = (area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))/volume(m,k,j,i);
    }
    u0(m,IM1,k,j,i) += bdt*(fg + d*(gr(m,k,i) + om2*vp));
    u0(m,IM3,k,j,i) += bdt*d*(gp(m,k,i) - om2*vr);
    // relaxation target: the envelope profile inside; outside, material by density
    // relative to the hydrostatic ambient: denser than env_amb_k x rho_amb(r,phi) is
    // stream/atmosphere gas (c_ph^2), the rest is the hot ambient (c_amb^2).  Choosing
    // by temperature or by an absolute density instead heated the stream's thin
    // leading edge and the atmosphere's top to c_amb (precursor jets at 1700 km/s and
    // an evaporation wind; tests g_half, a2, 10-06/07).
    Real ct = c2r(m,k,i);
    // Time constants: env_t_relax_env inside, env_t_relax_stream for the stream/atmosphere
    // gas, env_t_relax for the hot ambient; t <= 0 switches the relaxation off (adiabatic).
    const bool inside = (ct > cph2*(1.0 + 1.0e-9));
    Real tr = tri;
    if (!inside) {
      const bool dense = (d > kamb*ramb(m,k,i));
      ct = dense ? ((x1v(m,i) <= rhot) ? cph2 : cstr2) : camb2;
      tr = dense ? trs : tro;
    }
    // local cooling time: t = e_gas/(c a T^4) (1/(kappa_abs rho) + 3 kappa rho H^2): thin
    // emission by the absorption (Kramers) opacity, diffusion by ES + Kramers; H = (P/rho)/g,
    // g = G M_a / r^2 (cgs inside, converted to code time); not for the hot ambient
    if (cool && (inside || d > kamb*ramb(m,k,i)) && e > 0.0 && d > 0.0) {
      const Real c2 = gm1*e/d;                                   // (km/s)^2
      const Real tk = tfac*c2;                                   // K
      const Real rc = d*rhou;                                    // g/cc
      const Real kab = kk0*rc*pow(tk, -3.5);                     // cm^2/g, Kramers
      const Real hl = fmin(c2*SQR(x1v(m,i))/gma, x1v(m,i))*lu;    // cm
      const Real eg = e*rhou*vu*vu;                              // erg/cc
      const Real at4 = 7.5657e-15*SQR(SQR(tk));                  // erg/cc
      const Real kr = (kes + kab)*rc;
      tr = eg/(2.99792458e10*at4)*(1.0/(kab*rc) + 3.0*kr*hl*hl)/tu;
      tr = fmax(tr, tcmin);
    }
    if (gen) {
      // thermo = general: the targets are temperatures.  The envelope and the dense gas
      // below env_r_hot (the stellar atmosphere) relax to the column T(psi), so the
      // static IC is not touched; the stream above env_r_hot to env_t_stream; the hot
      // ambient to env_t_amb.  (env_cool is refused with thermo = general.)
      if (tr > 0.0 && d > 0.0) {
        Real tt;
        if (inside) {
          tt = tcol(m,k,i);
        } else if (d > kamb*ramb(m,k,i)) {
          tt = (x1v(m,i) <= rhot) ? tcol(m,k,i) : tstr;
        } else {
          tt = tamb;
        }
        u0(m,IEN,k,j,i) += (eosd.EnergyFromTemperature(d, tt) - e)*(1.0 - exp(-bdt/tr));
      }
    } else if (tr > 0.0) {
      u0(m,IEN,k,j,i) += (d*ct/gm1 - e)*(1.0 - exp(-bdt/tr));
    }
    // floor sponge (outside the envelope, rho < env_sponge_rho rho_amb): the rotating-
    // frame velocity decays at 1/t_sponge and is capped at env_v_floor_max; the kinetic
    // energy removed leaves the total energy (no heating of the floor gas)
    const Real du = u0(m,IDN,k,j,i);
    if (!inside && d < fsp*rfl && du > 0.0) {
      const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
      const Real v = sqrt(m1*m1 + m2*m2 + m3*m3)/du;
      Real f = (tsp > 0.0) ? exp(-bdt/tsp) : 1.0;
      if (v*f > vcap) f = vcap/v;
      u0(m,IM1,k,j,i) = f*m1;
      u0(m,IM2,k,j,i) = f*m2;
      u0(m,IM3,k,j,i) = f*m3;
      u0(m,IEN,k,j,i) -= (1.0 - f*f)*0.5*(m1*m1 + m2*m2 + m3*m3)/du;
    }
  });
  Real derel = 0.0;
  if (rad) {
    derel = RyPerEsum(pmbp) - erel0;     // per rank, as the flux registers
    // the WORK of the reference radiation force a_ref (force_reference_work = split; the
    // module pays it from the radiation energy), as he_star_m1's HeStarGravity
    auto *pm1 = pmbp->pradm1;
    if (pm1->fref_wsplit) {
      auto aref = pm1->arad_ref;
      if (pm1->FrefWaccOn()) {
        auto wacc = pm1->FrefWacc();
        const Real g0 = pm1->FrefWaccGam0(bdt, pm->dt);
        par_for("ryper_fws_x", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const Real w = bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
          u0(m,IEN,k,j,i) += w;
          wacc(m,k,j,i) = g0*wacc(m,k,j,i) + w;
        });
      } else {
        par_for("ryper_fws", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          u0(m,IEN,k,j,i) += bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
        });
      }
    }
  }

  if (nstages_ <= 0) return;
  auto &flx = ph->uflx.x1f;
  auto &xf = pmbp->pcoord->xx1f;
  auto &x3v = pmbp->pcoord->x3v;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto ira = iracc_;
  const Real om = rp_.omega, phis = phi_s_, wwin = nsig_*sig_;
  const bool strm = stream_on_;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj;
  DvceArray4D<Real> ff0;
  if (rad) ff0 = pmbp->pradm1->f0x1;
  const bool haveff = rad && (ff0.extent_int(0) > 0);
  const Real fin = rad_fin_;
  array_sum::GlobalSum fs;
  Kokkos::parallel_reduce("ryper_flxe", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/(nk*nj);
    const int k = (idx - m*nk*nj)/nj + ks;
    const int j = idx%nj + js;
    array_sum::GlobalSum v;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      const Real a = area1(m,k,j,is), r = xf(m,is);
      const Real fd = flx(m,IDN,k,j,is), f3 = flx(m,IM3,k,j,is);
      v.the_array[2] = -fd*a;
      v.the_array[3] = -(f3 + om*r*fd)*r*a;
      v.the_array[4] = -(f3 - fd*w0(m,IVZ,k,j,is))*r*a;
      if (rad) v.the_array[9] = fin*a;
    }
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      const Real a = area1(m,k,j,ie+1), r = xf(m,ie+1);
      const Real fd = flx(m,IDN,k,j,ie+1), f3 = flx(m,IM3,k,j,ie+1);
      const bool win = strm && (fabs(WrapPhi(x3v(m,k) - phis)) < wwin);
      if (win) {
        v.the_array[0] = -fd*a;
        v.the_array[8] = -(f3 + om*r*fd)*r*a;
      } else {
        v.the_array[1] = fd*a;
      }
      v.the_array[5] = (f3 + om*r*fd)*r*a;
      if (rad) {
        if (haveff) v.the_array[10] = ff0(m,k,j,ie+1)*a;
        v.the_array[11] = flx(m,IEN,k,j,ie+1)*a;
      }
    }
    const int ir = ira(m);
    if (ir >= 0) {
      const Real a = area1(m,k,j,ir), r = xf(m,ir);
      const Real fd = flx(m,IDN,k,j,ir), f3 = flx(m,IM3,k,j,ir);
      v.the_array[6] = -fd*a;
      v.the_array[7] = -(f3 + om*r*fd)*r*a;
    }
    sum += v;
  }, Kokkos::Sum<array_sum::GlobalSum>(fs));

  const int s = stage_ctr_ % nstages_;
  const int nacc = rad ? kNaccR : kNacc;
  if (s == 0) {
    for (int n=0; n<nacc; ++n) acc1_[n] = acc0_[n];
  }
  for (int n=0; n<kNacc; ++n) {
    acc0_[n] = rk_g0_[s]*acc0_[n] + rk_g1_[s]*acc1_[n] + bdt*fs.the_array[n];
  }
  if (rad) {
    for (int n=kNacc; n<kNaccR - 1; ++n) {
      acc0_[n] = rk_g0_[s]*acc0_[n] + rk_g1_[s]*acc1_[n] + bdt*fs.the_array[n];
    }
    acc0_[kNaccR-1] = rk_g0_[s]*acc0_[kNaccR-1] + rk_g1_[s]*acc1_[kNaccR-1] + derel;
  }
  stage_ctr_++;
}

//----------------------------------------------------------------------------------------
//! stage 2 history: Mdom Jdom, Menv Jenv (cells with r < R_acc), then the nine cumulative
//! boundary integrals and their rates over the last history interval
void RyPerHistEnv(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = pmbp->phydro->u0;
  auto &vol = pmbp->pcoord->volume;
  auto &x1v = pmbp->pcoord->x1v;
  const Real om = rp_.omega, racc = env_rmeas_;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  // <rad_m1>: E_rad, gas + radiation energy, kinetic energy below r_meas, and the
  // comoving x1 face fluxes at r_meas and r_out
  const bool rad = rad_;
  const int ie = indcs.ie;
  DvceArray5D<Real> ur;
  DvceArray4D<Real> ff0;
  if (rad) {
    ur = pmbp->pradm1->u0;
    ff0 = pmbp->pradm1->f0x1;
  }
  const bool haveff = rad && (ff0.extent_int(0) > 0);
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto ira = iracc_;
  array_sum::GlobalSum s;
  Kokkos::parallel_reduce("ryper_histe", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx%nx1 + is;
    array_sum::GlobalSum v;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
    const Real dv = vol(m,k,j,i), r = x1v(m,i), d = u0(m,IDN,k,j,i);
    const Real jz = dv*r*(u0(m,IM3,k,j,i) + om*r*d);
    v.the_array[0] = dv*d;
    v.the_array[1] = jz;
    if (r < racc) {
      v.the_array[2] = dv*d;
      v.the_array[3] = jz;
    }
    if (rad) {
      const Real er = ur(m,radm1::M1_E,k,j,i);
      v.the_array[4] = dv*er;
      v.the_array[5] = dv*(u0(m,IEN,k,j,i) + er);
      if (r < racc) {
        v.the_array[6] = dv*0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                                 + SQR(u0(m,IM3,k,j,i)))/d;
      }
      if (haveff) {
        if (i == ira(m)) v.the_array[7] = ff0(m,k,j,i)*area1(m,k,j,i);
        if (i == ie && mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
          v.the_array[8] = ff0(m,k,j,i+1)*area1(m,k,j,i+1);
        }
      }
    }
    sum += v;
  }, Kokkos::Sum<array_sum::GlobalSum>(s));

  const char *lab[22] = {"Mdom", "Jdom", "Menv", "Jenv",
                         "Min", "Mout", "Mwal", "Jwal", "Jstr", "Jout", "MR", "JR", "Jin",
                         "dMin", "dMout", "dMwal", "dJwal", "dJstr", "dJout", "dMR",
                         "dJR", "dJin"};
  pdata->nhist = 22;
  for (int n=0; n<22; ++n) pdata->label[n] = lab[n];
  for (int n=0; n<4; ++n) pdata->hdata[n] = s.the_array[n];
  if (rad) {
    // <rad_m1>: the nine RATE columns are replaced (they are the differences of the
    // cumulative columns 4-12 over the history interval; NHISTORY_VARIABLES = 22 is
    // left alone so that no other build changes):
    //   Lmeas, Ltop  comoving x1 radiative flux through r_meas and r_out, in Lsun of the
    //                FULL sphere (band sum x 4 pi / band solid angle)
    //   Erad, Etot   radiation energy and gas (incl. KE and rho Phi) + radiation energy
    //   KE           kinetic energy (rotating frame) below r_meas
    //   Ebnd         CUMULATIVE energy through the boundaries: radiation in at r_in
    //                (imposed) - radiation out at r_out - hydro energy flux out at r_out
    //   Erel         CUMULATIVE energy put in by the T relaxation (ambient) and the floor
    //                sponge; energy conservation: Etot - Etot(0) = Ebnd + Erel (+ floors)
    //   Pic, Picmax  mean Picard passes per implicit solve since the previous output, and
    //                the largest pass count of any solve so far
    // (all energies code units of the band; cumulative columns restart from 0)
    const char *lr[9] = {"Lmeas", "Ltop", "Erad", "Etot", "KE", "Ebnd", "Erel", "Pic",
                         "Picmax"};
    for (int n=0; n<9; ++n) pdata->label[13+n] = lr[n];
    for (int n=0; n<kNacc; ++n) pdata->hdata[4+n] = acc0_[n];
    pdata->hdata[13] = s.the_array[7]*rad_lfac_;
    pdata->hdata[14] = s.the_array[8]*rad_lfac_;
    pdata->hdata[15] = s.the_array[4];
    pdata->hdata[16] = s.the_array[5];
    pdata->hdata[17] = s.the_array[6];
    pdata->hdata[18] = acc0_[9] - acc0_[10] - acc0_[11];
    pdata->hdata[19] = acc0_[12];
    pdata->hdata[20] = 0.0;
    pdata->hdata[21] = 0.0;
    auto *pm1 = pmbp->pradm1;
    if (global_variable::my_rank == 0) {
      const Real dn = pm1->impl_nstep - rad_pn0_, ds = pm1->impl_itsum - rad_ps0_;
      pdata->hdata[20] = (dn > 0.0) ? ds/dn : 0.0;
      pdata->hdata[21] = pm1->impl_itmax;
      rad_pn0_ = pm1->impl_nstep;
      rad_ps0_ = pm1->impl_itsum;
    }
    hist_tprev_ = pm->time;
    return;
  }
  const Real t = pm->time;
  const Real dt = (hist_tprev_ >= 0.0) ? (t - hist_tprev_) : 0.0;
  for (int n=0; n<kNacc; ++n) {
    pdata->hdata[4+n] = acc0_[n];
    pdata->hdata[13+n] = (dt > 0.0) ? (acc0_[n] - hist_prev_[n])/dt : 0.0;
    hist_prev_[n] = acc0_[n];
  }
  hist_tprev_ = t;
}

//! release the namespace-scope Views before Kokkos::finalize
void RyPerFinal(ParameterInput *pin, Mesh *pm) {
  gr_ = DvceArray3D<Real>();
  gp_ = DvceArray3D<Real>();
  cref2_ = DvceArray3D<Real>();
  tcol_ = DvceArray3D<Real>();
  rhoamb_ = DvceArray3D<Real>();
  gpl_ = DvceArray3D<Real>();
  iracc_ = DvceArray1D<int>();
}

//! thermo = general: log T, log rho of the column at psi (linear in psi, clamped)
void ColumnAt(const std::vector<Real> &ps, const std::vector<Real> &lt,
              const std::vector<Real> &lr, const Real psi, Real &t, Real &r) {
  const int n = static_cast<int>(ps.size());
  if (psi <= ps[0]) {
    t = lt[0]; r = lr[0];
    return;
  }
  if (psi >= ps[n-1]) {
    t = lt[n-1]; r = lr[n-1];
    return;
  }
  const int hi = static_cast<int>(std::upper_bound(ps.begin(), ps.end(), psi) - ps.begin());
  const int lo = hi - 1;
  const Real w = (psi - ps[lo])/(ps[hi] - ps[lo]);
  t = (1.0 - w)*lt[lo] + w*lt[hi];
  r = (1.0 - w)*lr[lo] + w*lr[hi];
}

//! <rad_m1>: y(psi) of the column (rcol_psi_ ascending), linear, clamped at the ends
Real RadColAt(const std::vector<Real> &y, const Real psi) {
  const auto &ps = rcol_psi_;
  const int n = static_cast<int>(ps.size());
  if (psi <= ps[0]) return y[0];
  if (psi >= ps[n-1]) return y[n-1];
  const int hi = static_cast<int>(std::upper_bound(ps.begin(), ps.end(), psi) - ps.begin());
  const int lo = hi - 1;
  const Real w = (psi - ps[lo])/(ps[hi] - ps[lo]);
  return (1.0 - w)*y[lo] + w*y[hi];
}

//! <rad_m1>: the merged opacity table format (he_star_m1's HsReadOpacityTable without the
//! low-density extension: below the grid the lookup holds the edge value)
void RyReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                        DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD,
                        Real &lt_lo, Real &lt_hi, Real &ld_lo, Real &ld_hi) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in ry_per_accretor: cannot open opacity table '" << fname
              << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  Real lt0 = 0.0, dlt = 0.0, ld0 = 0.0, dld = 0.0;
  bool have_grid = false;
  std::vector<Real> vals;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      if (!have_grid) {
        std::istringstream ss(line.substr(1));
        int a, b;
        Real c, d, e, g;
        if (ss >> a >> b >> c >> d >> e >> g) {
          nT = a; nD = b; lt0 = c; dlt = d; ld0 = e; dld = g;
          have_grid = true;
        }
      }
      continue;
    }
    vals.push_back(std::stod(line));
  }
  if (!have_grid || static_cast<int>(vals.size()) != nT*nD) {
    std::cout << "### FATAL ERROR in ry_per_accretor: opacity table '" << fname
              << "' has no grid line or the wrong value count" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Kokkos::realloc(tab, nT, nD);
  Kokkos::realloc(lT, nT);
  Kokkos::realloc(lD, nD);
  auto htab = Kokkos::create_mirror_view(tab);
  auto hlT = Kokkos::create_mirror_view(lT);
  auto hlD = Kokkos::create_mirror_view(lD);
  for (int i=0; i<nT; ++i) {
    hlT(i) = lt0 + i*dlt;
    for (int j=0; j<nD; ++j) htab(i,j) = vals[i*nD + j];
  }
  for (int j=0; j<nD; ++j) hlD(j) = ld0 + j*dld;
  Kokkos::deep_copy(tab, htab);
  Kokkos::deep_copy(lT, hlT);
  Kokkos::deep_copy(lD, hlD);
  lt_lo = lt0; lt_hi = lt0 + (nT - 1)*dlt;
  ld_lo = ld0; ld_hi = ld0 + (nD - 1)*dld;
}

//----------------------------------------------------------------------------------------
//! <rad_m1> on the envelope (stage S3): units and keys checks, the Rosseland + Planck
//! tables, and the RADIATION-FORCE REFERENCE of the column.
//! The column (S0 script) is in hydrostatic balance with P_gas + P_rad along psi:
//! dP_gas/dpsi = rho (1 - Gamma), Gamma = kappa_t F/(c g) its Eddington factor.  In the
//! diffusion limit on Roche equipotentials F is proportional to g = |grad Phi_wb| at fixed
//! psi (von Zeipel), so Gamma is a function of psi alone and the radiation force is
//! -grad G(psi), G(psi) = int Gamma dpsi = -int a_ref dr_eq along the column (a_ref =
//! kappa_t F/c, the column's own radiative acceleration, kappa_t from THIS run's table
//! lookup).  The x1 well-balanced pair then carries Phi_eff = Phi_wb + G(psi) (constant
//! above the envelope top, where psi is capped), the module the x1 reference a_ref =
//! -dG/dr (face difference, so the WB pair and the module subtract the same force) and
//! applies only the residual, and the plain gravity above hydro/wb_rmax / env_wb_depth
//! uses Phi_eff too.  gr_ (the -d(Phi - Phi_wb)/dr remainder) is unchanged.  The column
//! must come from the SAME opacity table: the run recomputes kappa_t, and the imposed
//! bottom flux must equal the column's F(r_in) (checked to 1e-3).
void RadSetup(ParameterInput *pin, MeshBlockPack *pmbp, const std::vector<Real> &cpsi,
              const std::vector<Real> &clt, const std::vector<Real> &clr,
              const std::vector<Real> &crq, const std::vector<Real> &cfl) {
  const bool root = (global_variable::my_rank == 0);
  auto fatal = [&](const std::string &msg) {
    std::cout << "### FATAL ERROR in ry_per_accretor (<rad_m1>): " << msg << std::endl;
    std::exit(EXIT_FAILURE);
  };
  auto *pm1 = pmbp->pradm1;
  auto *phd = pmbp->phydro;
  auto &eos = phd->peos->eos_data;
  // ---- units: Rsun, km/s, the fixed density unit; the module's constants must match
  const Real du = pmbp->punit->density_cgs(), lu = pmbp->punit->length_cgs();
  const Real vu = pmbp->punit->velocity_cgs(), tk = eos.temp_cgs;
  const Real eu = du*vu*vu, fu = eu*vu;          // erg/cc, erg/cm^2/s per code unit
  const Real cl = pm1->c_light, ar = pm1->arad;
  {
    std::ostringstream os;
    bool bad = false;
    auto chk = [&](const char *nm, const Real have, const Real want, const Real tol) {
      if (!(std::fabs(have/want - 1.0) <= tol)) {
        os << nm << " = " << have << " but this run's units need " << want << "; ";
        bad = true;
      }
    };
    chk("<rad_m1>/temp_unit_kelvin", pm1->otab.tunit, tk, 1.0e-5);
    chk("<rad_m1>/rho_unit_cgs", pm1->otab.dunit, du, 1.0e-5);
    chk("<rad_m1>/kappa_unit", pm1->otab.kunit, du*lu, 1.0e-5);
    chk("<rad_m1>/c_light", cl, 2.99792458e10/vu, 1.0e-6);
    chk("<rad_m1>/arad", ar, 7.5657e-15*SQR(SQR(tk))/eu, 1.0e-4);
    if (bad) fatal(os.str());
  }
  if (pm1->opacity_type != radm1::M1_OPAC_TABLE) fatal("needs <rad_m1>/opacity = table");
  if (pm1->force_ref != radm1::M1_FREF_WB_ARAD) {
    fatal("needs <rad_m1>/force_reference = wb_arad (Phi_eff carries the column's "
          "radiation force)");
  }
  if (pm1->transport == radm1::M1_TRANSPORT_EXPLICIT) fatal("needs implicit transport");
  if (pin->GetString("rad_m1", "implicit_bc_x1min").compare("flux") != 0 ||
      pin->GetString("rad_m1", "implicit_bc_x1max").compare("marshak") != 0) {
    fatal("needs implicit_bc_x1min = flux (the star's L at r_in) and implicit_bc_x1max = "
          "marshak (vacuum top)");
  }
  // ---- the tables (one grid), handed to the module
  Real lt_lo, lt_hi, ld_lo, ld_hi;
  {
    const std::string rt = pin->GetString("problem", "env_opac_table");
    const std::string pt = pin->GetString("problem", "env_planck_table");
    DvceArray2D<Real> krt, kpt;
    DvceArray1D<Real> mlT, mlD, plT, plD;
    int mnT = 0, mnD = 0, pnT = 0, pnD = 0;
    Real a1, a2, a3, a4;
    RyReadOpacityTable(rt, krt, mlT, mlD, mnT, mnD, lt_lo, lt_hi, ld_lo, ld_hi);
    RyReadOpacityTable(pt, kpt, plT, plD, pnT, pnD, a1, a2, a3, a4);
    if (mnT != pnT || mnD != pnD || std::fabs(a1 - lt_lo) + std::fabs(a2 - lt_hi) +
        std::fabs(a3 - ld_lo) + std::fabs(a4 - ld_hi) > 1.0e-9) {
      fatal("the Rosseland and Planck tables are not on the same grid");
    }
    pm1->SetOpacityTables(krt, kpt, mlT, mlD, mnT, mnD);
    if (root) {
      std::printf("ry_per_accretor: <rad_m1> tables %s + %s: log T %.3f .. %.3f, log rho "
                  "%.3f .. %.3f (edge-held outside); absorption mask below %.4g code "
                  "(%.4g g/cc), kappa_s there %.4g cm^2/g\n", rt.c_str(), pt.c_str(),
                  lt_lo, lt_hi, ld_lo, ld_hi, pm1->otab.amask_rho, pm1->otab.amask_rho*du,
                  pm1->otab.amask_kes/pm1->otab.kunit);
    }
  }
  // ---- kappa_t along the column with this run's lookup (device), a_ref, G(psi)
  const int nc = static_cast<int>(cpsi.size());
  std::vector<Real> kt(nc);
  {
    DvceArray1D<Real> dd("ryper_cd", nc), dt("ryper_ct", nc), dk("ryper_ck", nc);
    auto hd = Kokkos::create_mirror_view(dd);
    auto ht = Kokkos::create_mirror_view(dt);
    for (int n=0; n<nc; ++n) {
      hd(n) = std::exp(clr[n])/du;
      ht(n) = std::exp(clt[n])/tk;
    }
    Kokkos::deep_copy(dd, hd);
    Kokkos::deep_copy(dt, ht);
    radm1::M1OpacTab ot = pm1->otab;
    par_for("ryper_ckt", DevExeSpace(), 0, nc - 1, KOKKOS_LAMBDA(const int n) {
      Real op, oe, of, os;
      radm1::M1TableOpacities(ot, dd(n), dt(n), op, oe, of, os);
      dk(n) = of + os;
    });
    auto hk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dk);
    for (int n=0; n<nc; ++n) kt[n] = hk(n);
  }
  rcol_psi_ = cpsi;
  rcol_g_.assign(nc, 0.0);
  rcol_f_.assign(nc, 0.0);
  std::vector<Real> aref(nc);
  for (int n=0; n<nc; ++n) {
    rcol_f_[n] = cfl[n]/fu;
    aref[n] = kt[n]*rcol_f_[n]/cl;
  }
  for (int n=1; n<nc; ++n) {      // psi ascending = r_eq descending
    rcol_g_[n] = rcol_g_[n-1] + 0.5*(aref[n-1] + aref[n])*(crq[n-1] - crq[n]);
  }
  const RocheParams p = rp_;
  const Real dlt = env_dlt_, rtop = env_rtop_, rspin = env_rspin_, phtop = env_phtop_;
  const Real phis = env_phis_;
  // ---- the imposed bottom flux must be the column's F at psi(r_in, phi 90)
  const Real r_in = pmbp->pmesh->mesh_size.x1min;
  const Real psi_in = phis - PhiWB(p, dlt, rtop, rspin, phtop, r_in, 0.5*M_PI);
  rad_fin_ = pin->GetReal("rad_m1", "implicit_flux_x1min");
  const Real fcol = RadColAt(rcol_f_, psi_in);
  if (!(std::fabs(rad_fin_/fcol - 1.0) <= 1.0e-3)) {
    std::ostringstream os;
    os << "<rad_m1>/implicit_flux_x1min = " << rad_fin_ << " but the column has F(r_in) = "
       << fcol << " (code units, " << fcol*fu << " erg/cm^2/s)";
    fatal(os.str());
  }
  // band solid angle -> full sphere, code luminosity -> Lsun
  {
    const auto &ms = pmbp->pmesh->mesh_size;
    const Real omg = (std::cos(ms.x2min) - std::cos(ms.x2max))*(ms.x3max - ms.x3min);
    rad_lfac_ = (4.0*M_PI/omg)*fu*lu*lu/3.828e33;
  }
  // ---- Phi_eff into the x1 WB pair (centres, x1 faces), the plain-gravity force above
  // the WB region (gpl_), the module's x1 reference a_ref = -dG/dr
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmbp->nmb_thispack;
  auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
  auto hpc = Kokkos::create_mirror_view(phd->phicc_wb);
  auto hpf = Kokkos::create_mirror_view(phd->phi_wb_x1f);
  auto hgl = Kokkos::create_mirror_view(gpl_);
  DvceArray4D<Real> aref_d("ryper_aref", nmb, n3, n2, n1);
  auto ha = Kokkos::create_mirror_view(aref_d);
  Real amax = 0.0;
  for (int m=0; m<nmb; ++m) {
    for (int k=0; k<n3; ++k) {
      const Real ph = hx3(m,k);
      std::vector<Real> pf(n1 + 1), gf(n1 + 1);
      for (int i=0; i<=n1; ++i) {
        const Real w = PhiWB(p, dlt, rtop, rspin, phtop, hx1f(m,i), ph);
        gf[i] = RadColAt(rcol_g_, phis - w);
        pf[i] = w + gf[i];
      }
      for (int i=0; i<n1; ++i) {
        const Real w = PhiWB(p, dlt, rtop, rspin, phtop, hx1(m,i), ph);
        const Real pc = w + RadColAt(rcol_g_, phis - w);
        const Real dr = hx1f(m,i+1) - hx1f(m,i);
        const Real a = -(gf[i+1] - gf[i])/dr;
        hgl(m,k,i) = -(pf[i+1] - pf[i])/dr;
        amax = std::max(amax, a);
        for (int j=0; j<n2; ++j) {
          hpc(m,k,j,i) = pc;
          hpf(m,k,j,i) = pf[i];
          if (i == n1 - 1) hpf(m,k,j,i+1) = pf[i+1];
          ha(m,k,j,i) = a;
        }
      }
    }
  }
  Kokkos::deep_copy(phd->phicc_wb, hpc);
  Kokkos::deep_copy(phd->phi_wb_x1f, hpf);
  Kokkos::deep_copy(gpl_, hgl);
  Kokkos::deep_copy(aref_d, ha);
  pm1->fref_wsplit_ok = true;       // RyPerSrcEnv gives the reference work (split)
  pm1->SetForceReference(aref_d);
  if (root) {
    const Real g_in = RadColAt(rcol_g_, psi_in);
    std::printf("ry_per_accretor: <rad_m1> force reference: G(psi) over the column (psi "
                "%.5g .. %.5g), G(r_in) %.6g (km/s)^2 = %.4f psi(r_in); max a_ref %.5g; "
                "bottom flux %.6g code = %.6g erg/cm^2/s (L %.5g Lsun at r_in, phi 90)\n",
                cpsi.front(), cpsi.back(), g_in, g_in/psi_in, amax, rad_fin_,
                rad_fin_*fu, 4.0*M_PI*SQR(r_in*lu)*rad_fin_*fu/3.828e33);
  }
}

//----------------------------------------------------------------------------------------
//! <rad_m1>: the radiation IC (fresh start): E = a T_col^4 (T_col(psi): the column's,
//! clamped at its top; in the ambient the radiation of the top of the atmosphere), F_r =
//! the column's F(psi) on cells and on the x1 faces (f0x1), F_theta = F_phi = 0; then the
//! startup range check of every radiatively active cell (rho >= the absorption mask)
//! against the opacity grid and the EOS table (fatal).
void RadIC(ParameterInput *pin, MeshBlockPack *pmbp) {
  const bool root = (global_variable::my_rank == 0);
  auto *pm1 = pmbp->pradm1;
  auto *phd = pmbp->phydro;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmbp->nmb_thispack;
  const RocheParams p = rp_;
  const Real dlt = env_dlt_, rtop = env_rtop_, rspin = env_rspin_, phtop = env_phtop_;
  const Real phis = env_phis_;
  auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
  DvceArray3D<Real> fc("ryper_rfc", nmb, n3, n1), ffc("ryper_rff", nmb, n3, n1 + 1);
  auto hfc = Kokkos::create_mirror_view(fc);
  auto hff = Kokkos::create_mirror_view(ffc);
  for (int m=0; m<nmb; ++m) {
    for (int k=0; k<n3; ++k) {
      for (int i=0; i<n1; ++i) {
        hfc(m,k,i) = RadColAt(rcol_f_, phis - PhiWB(p, dlt, rtop, rspin, phtop, hx1(m,i),
                                                    hx3(m,k)));
      }
      for (int i=0; i<=n1; ++i) {
        hff(m,k,i) = RadColAt(rcol_f_, phis - PhiWB(p, dlt, rtop, rspin, phtop, hx1f(m,i),
                                                    hx3(m,k)));
      }
    }
  }
  Kokkos::deep_copy(fc, hfc);
  Kokkos::deep_copy(ffc, hff);
  auto ur = pm1->u0;
  auto tc = tcol_;
  const Real ar = pm1->arad, efl = pm1->e_floor, cl = pm1->c_light;
  par_for("ryper_ric", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = tc(m,k,i);
    const Real e = fmax(ar*SQR(SQR(t)), efl);
    ur(m,radm1::M1_E,k,j,i) = e;
    ur(m,radm1::M1_F1,k,j,i) = fmin(fc(m,k,i), 0.9*cl*e);
    ur(m,radm1::M1_F2,k,j,i) = 0.0;
    ur(m,radm1::M1_F3,k,j,i) = 0.0;
  });
  auto ff0 = pm1->f0x1;
  if (ff0.extent_int(0) >= nmb) {
    const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
    const int e1 = std::min(ff0.extent_int(3) - 1, n1);
    par_for("ryper_ricf", DevExeSpace(), 0, nmb - 1, 0, e3, 0, e2, 0, e1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      ff0(m,k,j,i) = ffc(m,k,i);
    });
  }
  // range check of the radiatively active cells (active zones)
  {
    const auto &eosd = phd->peos->eos_data;
    const EOS_Data eos = eosd;
    auto u0 = phd->u0;
    auto phicc = phd->phicc0;
    const Real thr = pm1->otab.amask_rho, du = pm1->otab.dunit, tk = pm1->otab.tunit;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int ni = indcs.nx1, nji = indcs.nx2*ni, nkji = indcs.nx3*nji;
    Real mn[2], mx[2];
    for (int q=0; q<2; ++q) {
      Real lo = 1.0e300, hi = -1.0e300;
      Kokkos::parallel_reduce("ryper_rchk_lo", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              nmb*nkji),
      KOKKOS_LAMBDA(const int idx, Real &l) {
        const int m = idx/nkji;
        const int k = (idx - m*nkji)/nji + ks;
        const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
        const int i = idx%ni + is;
        const Real d = u0(m,IDN,k,j,i);
        if (!(d >= thr)) return;
        const Real e = u0(m,IEN,k,j,i) - 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                       + SQR(u0(m,IM3,k,j,i)))/d - d*phicc(m,k,j,i);
        const Real v = (q == 0) ? log10(d*du) : log10(eos.Temperature(d, e)*tk);
        l = fmin(l, v);
      }, Kokkos::Min<Real>(lo));
      Kokkos::parallel_reduce("ryper_rchk_hi", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              nmb*nkji),
      KOKKOS_LAMBDA(const int idx, Real &l) {
        const int m = idx/nkji;
        const int k = (idx - m*nkji)/nji + ks;
        const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
        const int i = idx%ni + is;
        const Real d = u0(m,IDN,k,j,i);
        if (!(d >= thr)) return;
        const Real e = u0(m,IEN,k,j,i) - 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                       + SQR(u0(m,IM3,k,j,i)))/d - d*phicc(m,k,j,i);
        const Real v = (q == 0) ? log10(d*du) : log10(eos.Temperature(d, e)*tk);
        l = fmax(l, v);
      }, Kokkos::Max<Real>(hi));
      mn[q] = lo;
      mx[q] = hi;
    }
    auto &ot = pm1->otab;
    auto hlT = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ot.lT);
    auto hlD = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ot.lD);
    const Real tlo = hlT(0), thi = hlT(ot.nT - 1), dlo = hlD(0), dhi = hlD(ot.nD - 1);
    const bool ok = (mn[0] >= dlo && mx[0] <= dhi && mn[1] >= tlo && mx[1] <= thi);
    const Real elo_d = pin->GetReal("hydro", "eos_logd_min");
    const Real ehi_d = pin->GetReal("hydro", "eos_logd_max");
    const Real elo_t = pin->GetReal("hydro", "eos_logt_min");
    const Real ehi_t = pin->GetReal("hydro", "eos_logt_max");
    const bool oke = (mn[0] >= elo_d && mx[0] <= ehi_d && mn[1] >= elo_t && mx[1] <= ehi_t);
    if (root) {
      std::printf("ry_per_accretor: <rad_m1> radiatively active cells (rho >= %.4g code): "
                  "log rho %.4f .. %.4f, log T %.4f .. %.4f; opacity grid log rho %.3f .. "
                  "%.3f, log T %.3f .. %.3f: %s; EOS table: %s\n", thr, mn[0], mx[0], mn[1],
                  mx[1], dlo, dhi, tlo, thi, ok ? "inside" : "OUTSIDE",
                  oke ? "inside" : "OUTSIDE");
    }
    if (!(ok && oke)) {
      std::cout << "### FATAL ERROR in ry_per_accretor (<rad_m1>): radiatively active "
                << "cells outside the opacity grid or the EOS table" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------------
//! thermo = general: the envelope IC in DISCRETE hydrostatic balance of the x1 WB pair
//! (the he_star_m1 he_ic_balance march, per (m,k) column, host).  T is FIXED at the
//! column's T(psi) (tcol_); the x1 WB background of cell i (WBBackgroundStencil with the
//! run's wb_option and the cells' fixed T) walked to its faces gives PL_i(d_i) and
//! PR_i(d_i); the face condition PR_i(d_i) = PL_{i+1}(d_{i+1}) is solved cell by cell
//! (secant in ln d) from the ANCHOR, the active cell nearest the photosphere (psi = 0),
//! which takes the column's density there; upward and downward.  Gas pressure only (no
//! radiation in the EOS): below the photosphere the gas carries the whole weight, so the
//! deep density exceeds the column's (which was built with P_gas + P_rad) by the
//! integrated Eddington factor.  Upward, the hot hydrostatic ambient takes over as in
//! the ideal-gas IC (pressure below the ambient's, r > r_top or Phi_wb >= the top cap).
//! The inner ghost (i = 0) and the outer ghost (n1 - 1) copy the isothermal walk; the
//! BCs rewrite them anyway.
void EnvICGeneral(MeshBlockPack *pmbp, const Real racc, const Real rho_amb,
                  const bool eqtop, const bool icpoly, const std::vector<Real> &cpsi,
                  const std::vector<Real> &clt, const std::vector<Real> &clr) {
  const bool root = (global_variable::my_rank == 0);
  auto *phd = pmbp->phydro;
  const EOS_Data eos = phd->peos->eos_data;
  const WBOption wbo = phd->wb_option;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmbp->nmb_thispack;
  const int is = indcs.is, ie = indcs.ie;
  const RocheParams p = rp_;
  const Real dfl = eos.dfloor, tamb = gen_tamb_, camb2 = env_camb2_, ramb0 = env_ramb_;
  const Real rtop = env_rtop_, phtop = env_phtop_, rspin = env_rspin_, phis = env_phis_;
  const Real vsp = (spin_ - 1.0)*rp_.omega;
  const Real aeq = RochePot(p, racc, 0.5*M_PI);
  // the envelope top cap in the WB pair's own potential (Phi_eff under <rad_m1>)
  const Real phtop_c = rad_ ? (phtop + RadColAt(rcol_g_, phis - phtop)) : phtop;
  auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);
  auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), phd->phicc_wb);
  auto hpf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), phd->phi_wb_x1f);
  auto hpt = Kokkos::create_mirror_view_and_copy(HostMemSpace(), phd->phicc0);
  auto ht = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tcol_);
  DvceArray3D<Real> dd("ryper_icd", nmb, n3, n1), de("ryper_ice", nmb, n3, n1);
  DvceArray3D<Real> dv("ryper_icv", nmb, n3, n1);
  auto hd = Kokkos::create_mirror_view(dd);
  auto he = Kokkos::create_mirror_view(de);
  auto hv = Kokkos::create_mirror_view(dv);
  std::vector<Real> d(n1), e(n1), t(n1), pw(n1), pf(n1 + 1);
  Real rmax = 0.0, rhoanc_min = 1e300, rhoanc_max = -1e300;
  int nfail = 0;
  for (int m=0; m<nmb; ++m) {
    for (int k=0; k<n3; ++k) {
      for (int i=0; i<n1; ++i) {
        t[i] = ht(m,k,i);
        pw[i] = hpc(m,k,0,i);
        pf[i] = hpf(m,k,0,i);
      }
      pf[n1] = hpf(m,k,0,n1);
      // walk of cell i at density di: pressures at its inner and outer face
      auto walk = [&](const int i, const Real di, Real &pl, Real &pr) {
        Real ei, pp, cr, ct, cv;
        eos.ThermoAt(di, t[i], ei, pp, cr, ct, cv);
        WBState s0, s1, s2, s3, s4;
        WBBackgroundStencil(eos, wbo, di, di, di, ei, ei, ei, pw[i-1], pf[i], pw[i],
                            pf[i+1], pw[i+1], s0, s1, s2, s3, s4, t[i-1], t[i], t[i+1]);
        pl = s1.p;
        pr = s3.p;
      };
      // solve for d_i: side 0 = its inner-face pressure, 1 = its outer-face pressure
      auto solve = [&](const int i, const int side, const Real target, const Real dg) {
        Real pl, pr;
        Real x0 = std::log(dg), x1 = x0 - 0.05;
        walk(i, std::exp(x0), pl, pr);
        Real f0 = ((side == 0) ? pl : pr)/target - 1.0;
        walk(i, std::exp(x1), pl, pr);
        Real f1 = ((side == 0) ? pl : pr)/target - 1.0;
        for (int it=0; it<80 && std::fabs(f1) > 1.0e-14 && f1 != f0; ++it) {
          Real x2 = x1 - f1*(x1 - x0)/(f1 - f0);
          x2 = std::min(std::max(x2, x1 - 2.0), x1 + 2.0);
          x0 = x1; f0 = f1; x1 = x2;
          walk(i, std::exp(x1), pl, pr);
          f1 = ((side == 0) ? pl : pr)/target - 1.0;
        }
        if (!(std::fabs(f1) < 1.0e-8)) nfail++;
        rmax = std::max(rmax, std::fabs(f1));
        return std::exp(x1);
      };
      // anchor: polytrope: the deepest active cell at EnvRho (as the ideal-gas IC; the
      // profile is resolved there, so P stays a function of Phi_wb and the tidal phi force
      // is balanced); column: the active cell with |psi| smallest at the column's rho
      int ia = is;
      if (icpoly) {
        d[ia] = EnvRho(phis - pw[ia], env_cph2_, env_np_, env_rhoph_);
      } else {
        Real best = 1.0e300;
        for (int i=is; i<=ie; ++i) {
          // <rad_m1>: pw is Phi_eff = Phi_wb + G; psi is defined by Phi_wb alone
          const Real psi = rad_ ? (phis - PhiWB(p, env_dlt_, rtop, rspin, phtop, hx1(m,i),
                                                hx3(m,k)))
                                : (phis - pw[i]);
          if (std::fabs(psi) < best) {
            best = std::fabs(psi);
            ia = i;
          }
        }
        Real lt, lr;
        ColumnAt(cpsi, clt, clr, phis - pw[ia], lt, lr);
        d[ia] = std::exp(lr)/eos.dens_cgs;
      }
      rhoanc_min = std::min(rhoanc_min, d[ia]);
      rhoanc_max = std::max(rhoanc_max, d[ia]);
      // downward to cell 1
      for (int i=ia-1; i>=1; --i) {
        Real pl, pr;
        walk(i+1, d[i+1], pl, pr);
        d[i] = solve(i, 1, pl, d[i+1]);
      }
      // upward to n1 - 2; the hot ambient takes over
      bool amb = false;
      for (int i=ia; i<=n1-1; ++i) {
        const Real pha = hpt(m,k,0,i) - (eqtop ? aeq : RochePot(p, racc, hx3(m,k)));
        const Real da = ramb0*std::exp(std::max(-pha/camb2, -700.0));
        if (!amb && i > ia) {
          if (i <= n1 - 2) {
            Real pl, pr;
            walk(i-1, d[i-1], pl, pr);
            d[i] = solve(i, 0, pr, d[i-1]);
          } else {
            d[i] = d[i-1];
          }
        }
        Real pc = 0.0;
        if (!amb) {
          Real ei, cr, ct, cv;
          eos.ThermoAt(d[i], t[i], ei, pc, cr, ct, cv);
        }
        if (i > ia && (amb || hx1(m,i) > rtop || pw[i] >= phtop_c || pc < da*camb2 ||
                       i == n1 - 1)) {
          amb = true;
          d[i] = std::max(da, rho_amb);
          t[i] = tamb;
        }
      }
      // inner ghost i = 0: the isothermal walk from cell 1
      {
        Real e1 = eos.EnergyFromTemperature(d[1], t[1]);
        Real dg = d[1], eg = e1, tg = t[1];
        WBAdvance(eos, 1, d[1], e1, pw[0] - pw[1], dg, eg, tg, t[1], 0.0, t[1], t[1]);
        d[0] = dg;
        t[0] = t[1];
      }
      for (int i=0; i<n1; ++i) {
        d[i] = std::max(d[i], dfl);
        hd(m,k,i) = d[i];
        he(m,k,i) = eos.EnergyFromTemperature(d[i], t[i]);
        hv(m,k,i) = (hx1(m,i) < rspin && t[i] != tamb) ? vsp*hx1(m,i) : 0.0;
      }
    }
  }
  if (root) {
    std::printf("ry_per_accretor: general IC (rank 0): anchor rho %.5g .. %.5g code, "
                "balance residual max %.3e, %d unconverged cells\n", rhoanc_min, rhoanc_max,
                rmax, nfail);
  }
  Kokkos::deep_copy(dd, hd);
  Kokkos::deep_copy(de, he);
  Kokkos::deep_copy(dv, hv);
  auto u0 = phd->u0;
  auto phicc = phd->phicc0;
  par_for("ryper_icg", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dc = dd(m,k,i), v3 = dv(m,k,i);
    u0(m,IDN,k,j,i) = dc;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = dc*v3;
    u0(m,IEN,k,j,i) = de(m,k,i) + 0.5*dc*v3*v3 + dc*phicc(m,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! stage 2 set-up, on every start (restart included): keys, potentials (true Phi into
//! phicc0/phi0 for etotgrav, Phi_wb into the x1 WB pair), the -dDelta/dr force (into
//! gr_), the relaxation target, the R_acc face index, and on a fresh start the envelope
void EnvSetup(ParameterInput *pin, MeshBlockPack *pmbp, const Real racc,
              const Real rho_amb, const bool restart) {
  const bool root = (global_variable::my_rank == 0);
  auto *phd = pmbp->phydro;
  auto &eos = phd->peos->eos_data;
  const Real cph = pin->GetOrAddReal("problem", "env_cs_ph", 15.5);   // km/s
  env_cph2_ = cph*cph;
  env_np_ = pin->GetOrAddReal("problem", "env_npoly", 3.0);
  env_rhoph_ = pin->GetOrAddReal("problem", "env_rho_ph", 20.0);
  env_tro_ = pin->GetOrAddReal("problem", "env_t_relax", 1.0e-4);
  env_tri_ = pin->GetOrAddReal("problem", "env_t_relax_env", env_tro_);
  // stream/atmosphere gas outside the envelope (rho > env_amb_k rho_amb); default
  // env_t_relax (bitwise as before).  0 = adiabatic stream (shock heating kept).
  env_trs_ = pin->GetOrAddReal("problem", "env_t_relax_stream", env_tro_);
  env_cool_ = pin->GetOrAddBoolean("problem", "env_cool", false);
  if (env_cool_) {
    cool_rhou_ = pin->GetReal("problem", "cool_rho_unit");
    cool_kes_ = pin->GetOrAddReal("problem", "cool_kappa_es", 0.34);
    cool_kk0_ = pin->GetOrAddReal("problem", "cool_kappa_k0", 1.5e24);
    const Real mu = pin->GetOrAddReal("problem", "cool_mu", 0.62);
    cool_tfac_ = mu*1.67262192e-24*kVel*kVel/1.380649e-16;
    cool_tmin_ = pin->GetOrAddReal("problem", "cool_t_min", 0.0);
  }
  if (gen_ && env_cool_) {
    std::cout << "### FATAL ERROR in ry_per_accretor: env_cool is ideal-gas only (not with "
              << "thermo = general)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string ws = pin->GetOrAddString("problem", "env_wall_slip", "free");
  env_noslip_ = (ws.compare("noslip") == 0);
  if (!env_noslip_ && ws.compare("free") != 0) {
    std::cout << "### FATAL ERROR in ry_per_accretor: env_wall_slip must be free | noslip"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  env_racc_ = racc;
  env_rhoamb_ = rho_amb;
  const Real csa = pin->GetOrAddReal("problem", "env_cs_amb", 300.0);   // km/s
  env_camb2_ = csa*csa;
  env_ramb_ = pin->GetOrAddReal("problem", "env_amb_rho", 1.0e-6);
  if (gen_) {
    // thermo = general: temperatures in kelvin (defaults: the c_s keys at mu 0.62), code
    // T = kelvin/eos.temp_cgs.  The ambient's P/rho (its isothermal scale) comes from the
    // EOS at (env_amb_rho, T_amb); it is fully ionized, so this is ideal to < 1e-3.
    const Real kmu = 0.62*kMH*kVel*kVel/kKB;     // K per (km/s)^2 at mu 0.62
    const Real tk = eos.temp_cgs;
    gen_tph_ = pin->GetOrAddReal("problem", "env_t_ph", kmu*cph*cph)/tk;
    gen_tamb_ = pin->GetOrAddReal("problem", "env_t_amb", kmu*csa*csa)/tk;
    const Real cst0 = pin->GetOrAddReal("problem", "env_cs_stream", cph);
    gen_tstr_ = pin->GetOrAddReal("problem", "env_t_stream", kmu*cst0*cst0)/tk;
    gen_tfl_ = (eos.tfloor > 0.0) ? eos.tfloor : 1.0e-3*gen_tph_;
    const Real ea = eos.EnergyFromTemperature(env_ramb_, gen_tamb_);
    env_camb2_ = eos.Pressure(env_ramb_, ea, gen_tamb_)/env_ramb_;
    if (root) {
      std::printf("ry_per_accretor: thermo = general: T_ph %.1f K, T_amb %.4g K (P/rho %.5g"
                  " = (%.3f km/s)^2), T_stream %.1f K, ghost T floor %.1f K; %.5g K per "
                  "code T, %.5g g/cc per code density\n", gen_tph_*tk, gen_tamb_*tk,
                  env_camb2_, std::sqrt(env_camb2_), gen_tstr_*tk, gen_tfl_*tk, tk,
                  eos.dens_cgs);
    }
  }
  env_kamb_ = pin->GetOrAddReal("problem", "env_amb_k", 3.0);
  env_tsp_ = pin->GetOrAddReal("problem", "env_t_sponge", 1.0e-4);
  env_fsp_ = pin->GetOrAddReal("problem", "env_sponge_rho", 10.0);
  env_vcap_ = pin->GetOrAddReal("problem", "env_v_floor_max", 100.0);
  const Real ntop = pin->GetOrAddReal("problem", "env_top_hp", 15.0);
  // free slip: the x1 flux at r_in is the exact mirror (zero mass, energy and phi-
  // momentum flux) whatever the ghosts hold; noslip needs the ghosts' v_phi
  phd->wall_closed_ix1 = !env_noslip_;
  env_dlt_ = 0.5*(spin_*spin_ - 1.0)*SQR(rp_.omega);
  const RocheParams p = rp_;
  const Real dlt = env_dlt_;
  // r_spin: the star's rotation (and the spin term of Phi_wb) ends here; default =
  // r_top (read again below), huge while r_top is searched
  const Real rspk = pin->GetOrAddReal("problem", "env_r_spin", 1.0e30);
  env_phis_ = RochePot(p, racc, 0.5*M_PI) - dlt*SQR(fmin(racc, rspk));
  // r_top: the first log-grid x1 face above R_acc where the isothermal (c_ph) hydrostatic
  // atmosphere over the photosphere has dropped by e^env_top_hp (default 15 scale
  // heights; independent of the floor; problem/env_r_top overrides).
  // Above it Phi_wb is flat, the gas unsupported ambient at rest.
  {
    const Real x1a = pmbp->pmesh->mesh_size.x1min, x1b = pmbp->pmesh->mesh_size.x1max;
    const Real dlg = std::log(x1b/x1a)/pmbp->pmesh->mesh_indcs.nx1;
    env_rtop_ = x1b;
    for (int n=1; n<100000; ++n) {
      const Real r = racc*std::exp(n*dlg);
      if (r >= x1b) break;
      const Real ex = -(RochePot(rp_, r, 0.5*M_PI) - dlt*SQR(fmin(r, rspk))
                        - env_phis_)/env_cph2_;
      if (ex < -ntop) {
        env_rtop_ = r;
        break;
      }
    }
    env_rtop_ = pin->GetOrAddReal("problem", "env_r_top", env_rtop_);
  }
  // problem/env_top_mode = sphere (default: the spherical r_top above, found along
  // phi = 90 deg) | equipotential: the supported envelope ends on the equipotential
  // psi = -env_top_hp c_ph^2 at every phi (Phi_wb capped at that value, r_top -> r_out).
  // With spin 1 the Roche equipotentials bulge along the binary axis (photosphere 9.0
  // at phi 90 vs ~9.5 at phi 0 for Plaskett): a spherical cut leaves dense gas against
  // the hot ambient there (blow-out at start).  The spin term keeps its spherical end
  // (env_r_spin, default the phi = 90 r_top).
  const std::string topm = pin->GetOrAddString("problem", "env_top_mode", "sphere");
  const bool eqtop = (topm.compare("equipotential") == 0);
  if (!eqtop && topm.compare("sphere") != 0) {
    std::cout << "### FATAL ERROR in ry_per_accretor: env_top_mode must be sphere | "
              << "equipotential" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real rtop90 = env_rtop_;
  Real rtopmax = rtop90;           // largest radius of the envelope top over phi
  if (eqtop) {
    env_phtop_ = env_phis_ + ntop*env_cph2_;
    env_rtop_ = pmbp->pmesh->mesh_size.x1max;
    const Real rsp = fmin(rspk, rtop90);
    for (int q=0; q<=360; ++q) {           // phi = 0 .. 2 pi, 1 deg steps
      const Real ph = q*M_PI/180.0;
      for (Real r=racc; r<env_rtop_; r+=1.0e-4*racc) {
        if (RochePot(p, r, ph) - dlt*SQR(fmin(r, rsp)) >= env_phtop_) {
          rtopmax = fmax(rtopmax, r);
          break;
        }
      }
    }
  }
  const Real rtop = env_rtop_, phtop = env_phtop_;
  env_rspin_ = fmin(rspk, rtop90);
  const Real rspin = env_rspin_;
  // the stream's own temperature (default c_ph: unchanged runs) and the radius below
  // which dense gas outside the star belongs to the (hot) stellar atmosphere
  {
    const Real cst = pin->GetOrAddReal("problem", "env_cs_stream", std::sqrt(env_cph2_));
    env_cstr2_ = cst*cst;
    env_rhot_ = pin->GetOrAddReal("problem", "env_r_hot", rtopmax);
  }
  const Real phis = env_phis_, cph2 = env_cph2_, np = env_np_, rhoph = env_rhoph_;

  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmbp->nmb_thispack;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto &x3v = pmbp->pcoord->x3v;
  auto &x3f = pmbp->pcoord->xx3f;
  // the TRUE potential (Roche incl. the orbital centrifugal term) for etotgrav
  {
    auto phicc = phd->phicc0;
    auto ph1 = phd->phi0.x1f, ph2 = phd->phi0.x2f, ph3 = phd->phi0.x3f;
    const int n1m1 = n1 - 1, n2m1 = n2 - 1, n3m1 = n3 - 1;
    par_for("ryper_phi", DevExeSpace(), 0, nmb - 1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real pc = RochePot(p, x1v(m,i), x3v(m,k));
      phicc(m,k,j,i) = pc;
      ph1(m,k,j,i) = RochePot(p, x1f(m,i), x3v(m,k));
      if (i == n1m1) ph1(m,k,j,i+1) = RochePot(p, x1f(m,i+1), x3v(m,k));
      ph2(m,k,j,i) = pc;
      if (j == n2m1) ph2(m,k,j+1,i) = pc;
      ph3(m,k,j,i) = RochePot(p, x1v(m,i), x3f(m,k));
      if (k == n3m1) ph3(m,k+1,j,i) = RochePot(p, x1v(m,i), x3f(m,k+1));
    });
  }
  // ... and Phi_wb = Phi - Delta for the x1 well-balanced pair
  {
    phd->EnableWBEffectivePotential();
    auto pwc = phd->phicc_wb, pwf = phd->phi_wb_x1f;
    const int n1m1 = n1 - 1;
    par_for("ryper_phiwb", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real ph = x3v(m,k);
      pwc(m,k,j,i) = PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,i), ph);
      pwf(m,k,j,i) = PhiWB(p, dlt, rtop, rspin, phtop, x1f(m,i), ph);
      if (i == n1m1) {
        pwf(m,k,j,i+1) = PhiWB(p, dlt, rtop, rspin, phtop, x1f(m,i+1), ph);
      }
    });
  }
  // gr_ becomes the x1 force the WB pair does not carry: -dDelta/dr; target P/rho
  Kokkos::realloc(cref2_, nmb, n3, n1);
  Kokkos::realloc(rhoamb_, nmb, n3, n1);
  {
    auto ra = rhoamb_;
    const Real ramb0 = env_ramb_, camb2 = env_camb2_, dmin = rho_amb;
    // anchor: R_acc on the same phi (sphere), or the photospheric equipotential
    const Real aeq = RochePot(p, racc, 0.5*M_PI);
    par_for("ryper_ramb", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      const Real dph = RochePot(p, x1v(m,i), x3v(m,k))
                       - (eqtop ? aeq : RochePot(p, racc, x3v(m,k)));
      ra(m,k,i) = fmax(ramb0*exp(fmax(-dph/camb2, -700.0)), dmin);
    });
  }
  Kokkos::realloc(gpl_, nmb, n3, n1);
  {
    auto gr = gr_, c2r = cref2_, gpl = gpl_;
    par_for("ryper_gre", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      const Real rl = x1f(m,i), rr = x1f(m,i+1), ph = x3v(m,k);
      const Real wr = PhiWB(p, dlt, rtop, rspin, phtop, rr, ph);
      const Real wl = PhiWB(p, dlt, rtop, rspin, phtop, rl, ph);
      gr(m,k,i) = -((RochePot(p, rr, ph) - wr) - (RochePot(p, rl, ph) - wl))/(rr - rl);
      gpl(m,k,i) = -(wr - wl)/(rr - rl);
      const Real psi = phis - PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,i), x3v(m,k));
      c2r(m,k,i) = EnvC2(psi, cph2, np);
    });
  }
  // thermo = general: the envelope column T(psi), rho(psi) from problem/env_ic_file
  // (S0 script; ASCII, '#' comments, columns psi [(km/s)^2, = Phi_s - Phi_wb at phi 90]
  // T [K] rho [g/cc] [...], psi ascending), interpolated linearly in psi on log T and
  // log rho; clamped at the file's ends.  tcol_ = T at the cell centres.
  // problem/env_ic (thermo = general): polytrope (default) = the env13 structure, T(psi)
  // = T_ph EnvC2(psi)/c_ph^2 (the analytic n-polytrope's P/rho shape in kelvin), anchored
  // at the deepest cell at EnvRho (as the ideal-gas IC), gas-pressure hydrostatic: the S2
  // test of the general EOS on the env13 star.  column = the S0 radiative column
  // (env_ic_file).  NOTE: that column is in hydrostatic balance with P_gas + P_rad; with
  // gas pressure alone at its T(psi) the discrete balance makes the deep envelope ~5e4x
  // denser (the integrated Eddington factor), so column needs the radiation force
  // reference of stage S3 (a_ref in Phi_wb) before it is a sensible star.
  std::vector<Real> cpsi, clt, clr, crq, cfl;
  bool icpoly = true;
  if (gen_) {
    const std::string icm = pin->GetOrAddString("problem", "env_ic", "polytrope");
    icpoly = (icm.compare("polytrope") == 0);
    if (!icpoly && icm.compare("column") != 0) {
      std::cout << "### FATAL ERROR in ry_per_accretor: env_ic must be polytrope | column"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (gen_ && icpoly) {
    Kokkos::realloc(tcol_, nmb, n3, n1);
    auto tc = tcol_;
    const Real tph = gen_tph_;
    par_for("ryper_tcolp", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      const Real psi = phis - PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,i), x3v(m,k));
      tc(m,k,i) = tph*EnvC2(psi, cph2, np)/cph2;
    });
    if (root) {
      const Real psi_in = phis - PhiWB(p, dlt, rtop, rspin, phtop, pmbp->pmesh->mesh_size.x1min, 0.5*M_PI);
      std::printf("ry_per_accretor: env_ic = polytrope: T(r_in) %.4g K, T_ph %.1f K\n",
                  tph*EnvC2(psi_in, cph2, np)/cph2*eos.temp_cgs, tph*eos.temp_cgs);
    }
  }
  if (gen_ && !icpoly) {
    const std::string fn = pin->GetString("problem", "env_ic_file");
    std::ifstream f(fn);
    if (!f.good()) {
      std::cout << "### FATAL ERROR in ry_per_accretor: cannot open env_ic_file '" << fn
                << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::string ln;
    while (std::getline(f, ln)) {
      if (ln.empty() || ln[0] == '#') continue;
      std::istringstream is(ln);
      Real a, b, c;
      if (!(is >> a >> b >> c)) continue;
      // <rad_m1>: also r_eq [Rsun] (column 7) and F [erg/cm^2/s] (column 8)
      Real c4, c5, c6, c7, c8;
      if (rad_ && !(is >> c4 >> c5 >> c6 >> c7 >> c8)) {
        std::cout << "### FATAL ERROR in ry_per_accretor: <rad_m1> needs the 8+ column S0 "
                  << "file (psi T rho P_gas P_rad kappa_R r_eq F ...)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!cpsi.empty() && !(a > cpsi.back())) {
        std::cout << "### FATAL ERROR in ry_per_accretor: env_ic_file psi not strictly "
                  << "ascending at psi = " << a << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cpsi.push_back(a);
      clt.push_back(std::log(b));
      clr.push_back(std::log(c));
      if (rad_) {
        crq.push_back(c7);
        cfl.push_back(c8);
      }
    }
    if (cpsi.size() < 2) {
      std::cout << "### FATAL ERROR in ry_per_accretor: env_ic_file has < 2 rows"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // range check: the column inside the EOS table (T and rho, cgs)
    {
      const auto &tb = eos.tbl;
      Real tlo = 1e300, thi = -1e300, rlo = 1e300, rhi = -1e300;
      for (size_t n=0; n<cpsi.size(); ++n) {
        tlo = std::min(tlo, clt[n]); thi = std::max(thi, clt[n]);
        rlo = std::min(rlo, clr[n]); rhi = std::max(rhi, clr[n]);
      }
      const Real l10 = std::log(10.0);
      const bool okt = (tlo/l10 >= tb.ymin) && (thi/l10 <= tb.ymax);
      const bool okr = (rlo/l10 >= tb.xmin) && (rhi/l10 <= tb.xmax);
      const Real tak = gen_tamb_*eos.temp_cgs;
      const bool oka = (std::log10(tak) <= tb.ymax);
      if (root) {
        std::printf("ry_per_accretor: env_ic_file %s: %d rows, psi %.5g .. %.5g, log T "
                    "%.4f .. %.4f, log rho %.4f .. %.4f; EOS table log T %.3f .. %.3f, "
                    "log rho %.3f .. %.3f: %s\n", fn.c_str(), static_cast<int>(cpsi.size()),
                    cpsi.front(), cpsi.back(), tlo/l10, thi/l10, rlo/l10, rhi/l10,
                    tb.ymin, tb.ymax, tb.xmin, tb.xmax,
                    (okt && okr && oka) ? "inside" : "OUTSIDE");
      }
      if (!(okt && okr && oka)) {
        std::cout << "### FATAL ERROR in ry_per_accretor: the env_ic_file column or T_amb "
                  << "is outside the EOS table" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    Kokkos::realloc(tcol_, nmb, n3, n1);
    auto ht = Kokkos::create_mirror_view(tcol_);
    auto hx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1v);
    auto hx3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x3v);
    const Real tk = eos.temp_cgs;
    for (int m=0; m<nmb; ++m) {
      for (int k=0; k<n3; ++k) {
        for (int i=0; i<n1; ++i) {
          const Real psi = phis - PhiWB(p, dlt, rtop, rspin, phtop, hx1(m,i), hx3(m,k));
          Real lt, lr;
          ColumnAt(cpsi, clt, clr, psi, lt, lr);
          ht(m,k,i) = std::exp(lt)/tk;
        }
      }
    }
    Kokkos::deep_copy(tcol_, ht);
    if (rad_) RadSetup(pin, pmbp, cpsi, clt, clr, crq, cfl);
  }
  // the measuring x1 face in each MeshBlock (or -1): r = R_acc exactly (default), or,
  // with problem/r_meas, the face nearest r_meas (MR/JR, Menv/Jenv then refer to it; a
  // sphere outside the tidally bulged photosphere, so the envelope does not breathe
  // through it)
  const bool mexact = !pin->DoesParameterExist("problem", "r_meas");
  env_rmeas_ = pin->GetOrAddReal("problem", "r_meas", racc);
  Kokkos::realloc(iracc_, nmb);
  int nfound = 0;
  {
    auto hx = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1f);
    auto hi = Kokkos::create_mirror_view(iracc_);
    for (int m=0; m<nmb; ++m) {
      hi(m) = -1;
      Real best = 1.0e300;
      for (int i=indcs.is; i<=indcs.ie+1; ++i) {
        if (mexact) {
          if (std::fabs(hx(m,i) - racc) < 1.0e-6*racc) hi(m) = i;
        } else if (std::fabs(hx(m,i) - env_rmeas_) < best) {
          best = std::fabs(hx(m,i) - env_rmeas_);
          hi(m) = i;
        }
      }
      if (!mexact && hi(m) >= 0) env_rmeas_ = hx(m,hi(m));
      if (hi(m) >= 0) nfound++;
    }
    Kokkos::deep_copy(iracc_, hi);
  }
  // problem/env_wb_depth > 0: the x1 well-balanced scheme is switched off (as above
  // hydro/wb_rmax) outside the equipotential through r = R_acc - env_wb_depth at
  // phi = 90 deg, i.e. at the same depth below the photosphere at every phi
  {
    const Real wdep = pin->GetOrAddReal("problem", "env_wb_depth", 0.0);
    if (wdep > 0.0) {
      const Real rw = racc - wdep;
      phd->wb_phimax = RochePot(p, rw, 0.5*M_PI) - dlt*SQR(fmin(rw, rspin));
      // <rad_m1>: phicc_wb holds Phi_eff = Phi_wb + G(psi)
      if (rad_) phd->wb_phimax += RadColAt(rcol_g_, env_phis_ - phd->wb_phimax);
      if (root) {
        std::printf("ry_per_accretor: x1 WB off outside the equipotential through r = "
                    "%.4f (phi 90), Phi_wb > %.6e\n", rw, phd->wb_phimax);
      }
    }
  }
  if (root && !mexact) {
    std::printf("ry_per_accretor: measuring face (MR, JR, Menv, Jenv) at r = %.5f\n",
                env_rmeas_);
  }
  if (root && eqtop) {
    std::printf("ry_per_accretor: envelope top on the equipotential psi = -%.1f c_ph^2 "
                "(Phi_wb cap %.6e), r_top %.4f at phi 90, max over phi %.4f\n", ntop,
                phtop, rtop90, rtopmax);
  }
  const Real r_in = pmbp->pmesh->mesh_size.x1min;
  if (root) {
    const Real psi_in = phis - PhiWB(p, dlt, rtop, rspin, phtop, r_in, 0.5*M_PI);
    std::printf("ry_per_accretor: ENVELOPE n %.2f, c_ph %.3f km/s, rho_ph %.4g, r_in %.5f"
                " = %.4f R_acc; at r_in P/rho %.4g (c_iso %.2f km/s), rho %.4e; r_top "
                "%.4f, t_relax %.3g / env %.3g / stream %.3g, wall %s; R_acc face found in "
                "%d of %d local blocks\n", np, cph, rhoph, r_in, r_in/racc,
                EnvC2(psi_in, cph2, np),
                std::sqrt(EnvC2(psi_in, cph2, np)), EnvRho(psi_in, cph2, np, rhoph), rtop,
                env_tro_, env_tri_, env_trs_, ws.c_str(), nfound, nmb);
  }
  if (nfound == 0 && root) {
    std::printf("ry_per_accretor: WARNING no x1 face at r = R_acc = %.6f on rank 0: "
                "MR/JR are zero there (choose x1min so that R_acc is a log-grid face)\n",
                racc);
  }
  if (restart) return;

  // ---- initial state: the envelope in DISCRETE hydrostatic equilibrium of the x1 WB
  // pair, column by column.  T = P/rho is the target profile at the cell centres; the
  // anchor is the deepest active cell at the analytic polytrope density (there the
  // profile is resolved, so P stays a function of Phi_wb and the tidal phi force is
  // balanced); from it the pressure is carried face by face with each cell's own
  // isothermal background,
  //   p_f = p_i exp(-(Phi_wb,f - Phi_wb,i)/T_i),  p_{i+1} = p_f exp((Phi_wb,f -
  //   Phi_wb,i+1)/T_{i+1}),
  // which is the exact rest state of the isothermal-option WB reconstruction and source
  // even where H_p is not resolved (the photosphere).  Upward the c_ph atmosphere ends at
  // the hot ambient's pressure or r_top; above, the hot hydrostatic ambient at rest.
  // Rotation (spin - 1) Omega r
  // (rotating frame) in the supported part (r < r_top).
  if (gen_) {
    EnvICGeneral(pmbp, racc, rho_amb, eqtop, icpoly, cpsi, clt, clr);
    if (rad_) RadIC(pin, pmbp);
    return;
  }
  const Real gm1 = eos.gamma - 1.0, dfl = eos.dfloor;
  const Real vsp = (spin_ - 1.0)*rp_.omega;
  const Real camb2 = env_camb2_, ramb = env_ramb_;
  const Real aeq = RochePot(p, racc, 0.5*M_PI);
  auto u0 = phd->u0;
  auto phicc = phd->phicc0;
  auto pwc = phd->phicc_wb, pwf = phd->phi_wb_x1f;
  auto c2r = cref2_;
  const int isa = indcs.is;
  {
    auto hi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), iracc_);
    for (int m=0; m<nmb; ++m) {
      if (hi(m) < 0) {
        std::cout << "### FATAL ERROR in ry_per_accretor: inner = envelope needs every "
                  << "MeshBlock to span r = R_acc (meshblock/nx1 = mesh/nx1, R_acc a "
                  << "log-grid face)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  par_for("ryper_ice", DevExeSpace(), 0, nmb - 1, 0, n3 - 1,
  KOKKOS_LAMBDA(const int m, const int k) {
    const int ia = isa;                   // the deepest active cell
    const Real psia = phis - pwc(m,k,0,ia);
    // pass 1 (down from the anchor, incl. the ghosts), pass 2 (up)
    for (int pass=0; pass<2; ++pass) {
      Real pp = EnvRho(psia, cph2, np, rhoph)*c2r(m,k,ia);
      Real tp = c2r(m,k,ia);
      bool amb = false;
      const int i0 = (pass == 0) ? ia : ia + 1;
      const int i1 = (pass == 0) ? -1 : n1;
      const int di = (pass == 0) ? -1 : 1;
      for (int i=i0; i!=i1; i+=di) {
        const Real t = c2r(m,k,i);
        Real pc = pp;
        if (i != ia) {
          const int ip = i - di;                       // previous cell
          const int f = (di > 0) ? i : i + 1;          // face between ip and i
          const Real pf = pp*exp(-(pwf(m,k,0,f) - pwc(m,k,0,ip))/tp);
          pc = pf*exp((pwf(m,k,0,f) - pwc(m,k,0,i))/t);
        }
        Real d = pc/t, tt = t;
        Real v3 = (x1v(m,i) < rspin) ? vsp*x1v(m,i) : 0.0;
        // the hot ambient: isothermal (c_amb) at rest, hydrostatic in the TRUE Roche
        // potential, env_amb_rho at R_acc; it takes over where the cold column's
        // pressure falls below its own (or above r_top)
        const Real pha = phicc(m,k,0,i) - (eqtop ? aeq : RochePot(p, racc, x3v(m,k)));
        const Real da = ramb*exp(fmax(-pha/camb2, -700.0));
        if (pass == 1 && (amb || x1v(m,i) > rtop || pwc(m,k,0,i) >= phtop ||
                          pc < da*camb2)) {
          amb = true;
          d = fmax(da, rho_amb); tt = camb2; v3 = 0.0;
        }
        d = fmax(d, dfl);
        for (int j=0; j<n2; ++j) {
          u0(m,IDN,k,j,i) = d;
          u0(m,IM1,k,j,i) = 0.0;
          u0(m,IM2,k,j,i) = 0.0;
          u0(m,IM3,k,j,i) = d*v3;
          u0(m,IEN,k,j,i) = d*tt/gm1 + 0.5*d*v3*v3 + d*phicc(m,k,j,i);
        }
        pp = pc; tp = t;
      }
    }
  });
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
  // stage 2: problem/inner = envelope (read only when named, so a stage-1 input's
  // parameter dump is unchanged)
  env_ = false;
  if (pin->DoesParameterExist("problem", "inner")) {
    std::string inner = pin->GetString("problem", "inner");
    env_ = (inner.compare("envelope") == 0);
    if (!env_ && inner.compare("surface") != 0) {
      fatal("problem/inner must be surface | envelope");
    }
  }
  auto &eos = pmbp->phydro->peos->eos_data;
  gen_ = false;
  // <rad_m1> (stage S3): only with inner = envelope, thermo = general, env_ic = column
  rad_ = (pmbp->pradm1 != nullptr);
  if (rad_) {
    if (!env_ || thermo.compare("general") != 0 ||
        pin->GetOrAddString("problem", "env_ic", "polytrope").compare("column") != 0) {
      fatal("<rad_m1> needs problem/inner = envelope, thermo = general and env_ic = column");
    }
  }
  if (env_) {
    auto *phd = pmbp->phydro;
    gen_ = (thermo.compare("general") == 0);
    if (gen_) {
      // thermo = general: tabulated general EOS, gas only, in Rsun / km/s units
      if (pin->GetString("hydro", "eos").compare("general") != 0 || !eos.tbl.active) {
        fatal("thermo = general needs <hydro>/eos = general with general_eos = table");
      }
      if (eos.tbl.radiation) fatal("thermo = general needs <hydro>/eos_radiation = false");
      if (phd->wb_option != WBOption::polytropic) {
        fatal("thermo = general needs <hydro>/wb_option = polytropic");
      }
      if (pmbp->punit == nullptr) fatal("thermo = general needs a <units> block");
      const Real lu = pmbp->punit->length_cgs(), vu = pmbp->punit->velocity_cgs();
      if (std::fabs(lu/kRsun - 1.0) > 1.0e-6 || std::fabs(vu/kVel - 1.0) > 1.0e-6 ||
          std::fabs(pmbp->punit->mu() - 1.0) > 1.0e-12) {
        fatal("thermo = general needs <units> length_cgs = 6.957e10 (Rsun), time_cgs = "
              "6.957e5 (velocity km/s) and mu = 1");
      }
      if (root) {
        std::printf("ry_per_accretor: thermo = general, density unit %.6g g/cc\n",
                    pmbp->punit->density_cgs());
      }
    } else if (thermo.compare("adiabatic") != 0) {
      fatal("inner = envelope needs thermo = adiabatic | general");
    }
    if (!gen_ && !eos.is_ideal) fatal("inner = envelope needs <hydro>/eos = ideal");
    if (!(phd->use_etotgrav && phd->use_wellbalance_dynamic && phd->use_wb_x1)) {
      fatal("inner = envelope needs <hydro>/etotgrav, wellbalance_dynamic and wb_x1 = "
            "true");
    }
  } else {
    if (thermo.compare("isothermal") != 0) {
      fatal("problem/thermo = '" + thermo + "': only isothermal is implemented for the "
            "absorbing surface (adiabatic needs problem/inner = envelope)");
    }
    if (eos.is_ideal) fatal("needs <hydro>/eos = isothermal");
  }
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
    if (env_) {
      std::printf("ry_per_accretor: spin %.3f (envelope v_phi %.3f km/s at R_acc, "
                  "rotating frame), stream %s\n", spin_, (spin_ - 1.0)*rp_.omega*racc,
                  stream_on_ ? "on" : "off");
    } else {
      std::printf("ry_per_accretor: spin %.3f (wall v_phi %.3f km/s at r_in, rotating "
                  "frame), inner %s/%s/%s, outer %s, stream %s, iso_cs %.3f km/s\n",
                  spin_, (spin_ - 1.0)*rp_.omega*r_in, slip.c_str(), irho.c_str(),
                  ivr.c_str(),
                  orho.c_str(), stream_on_ ? "on" : "off", eos.iso_cs);
    }
    if (!env_ && std::fabs(r_in - racc) > 1.0e-6*racc) {
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
  if (env_) {
    EnvSetup(pin, pmbp, racc, rho_amb, restart);
    user_bcs_func = RyPerBCEnv;
    user_srcs_func = RyPerSrcEnv;
    if (user_hist) user_hist_func = RyPerHistEnv;
    return;
  }
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
