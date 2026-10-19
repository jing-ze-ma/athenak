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
//! implicit_flux_x1min checked against env_lin_band_factor x the column's F(r_in) to
//! 1e-3; the factor, default 1, is the von Zeipel band equivalent of the 2-D slab and
//! must be 1 in 3-D); top: Marshak/vacuum
//! (implicit_bc_x1max = marshak, dark ghosts).  IC: E = a T_col^4, F_r = F_col(psi).
//! The hot ambient keeps its T relaxation and gets no absorption: <rad_m1>/
//! opac_abs_rho_max (module key) removes kappa_P, kappa_E below that density and keeps
//! electron scattering.  The relaxation/sponge energy is booked (history Erel).  The
//! relaxation inside the envelope and of the dense atmosphere should be off
//! (env_t_relax_env = env_t_relax_stream = 0).  History: the nine rate columns become
//! Lmeas Ltop Erad Etot KE Ebnd Erel Pic Picmax (RyPerHistEnv).
//!
//! PHYSICAL AMBIENT (problem/env_amb_mode = radeq; default hot = the hot hydrostatic ambient
//! above; S5 of accretor-rhd-1008, user 10-09).  The gas between the stars is tenuous gas at
//! radiative equilibrium: absorption ON (no opac_abs_rho_max), no T relaxation (env_t_relax
//! = 0; the M1 coupling itself sets T_gas = T_rad), held by the density floor rho_amb (=
//! hydro/dfloor; it cannot be hydrostatic at ~2e4 K over the domain).  IC: the column's
//! atmosphere ends where its density falls to rho_amb; above, rho_amb at T_col (the
//! radiation IC's own T, so gas and radiation start in equilibrium).  The floor is applied
//! HERE (RyPerSrcEnv, after the sources; the added gas has the cell's velocity and T), so
//! its mass is booked (history column 8 "Mfl", cumulative, in place of the wall stress Jstr)
//! and its energy goes into Erel with the sponge; the EOS floor then never acts.
//!
//! STREAM RADIATION (problem/stream_rad = true; default false: dark window).  The optically
//! thick stream enters with its radiation in equilibrium with its gas: in the window the M1
//! ghosts hold E = a T_s^4 and the lab-frame flux of optically thick advection F = (4/3) E v
//! (env_t_stream, the stream velocity), and since the implicit x1 solve uses a FACE-FLUX
//! boundary (the ghosts are not read), the outer face of the window columns gets the same
//! state as an INFLOW face (RadiationM1::SetX1maxBathColumns): comoving flux 0 (= the
//! ghosts' F - (4/3) v E), lab-frame inflow (1 + chi) v E_s (= (4/3) v E_s, thick limit).
//! Elsewhere the outer face is vacuum/Marshak, unchanged.  The advective radiation energy
//! through r_out (module's ibadv_x1max) is booked in Ebnd with the comoving face flux.
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
//! STREAM AS A VOLUME SOURCE + DONOR SINK (accretor-bigdom-1015, design sect. 13; keys read
//! only when given, defaults = the boundary window above):
//!  problem/stream_inject = volume: the domain extends beyond the old r_out (e.g. 1.5 d_L1)
//!    and the stream is injected in a NOZZLE, a shell of problem/stream_ncell (2) cells
//!    whose inner face is the x1 face nearest problem/stream_r (default x1max), inside
//!    the same phi window: after the sources of every stage those cells are RESET to the
//!    window state of the old ghosts (rho_s Gaussian, ballistic v at stream_r, T(phi),
//!    E = a T^4, F = (4/3) E v) -- an internal Dirichlet region, so the face below it sees
//!    exactly what the old r_out face saw.  The mass, inertial Lz, gas energy and
//!    radiation energy the reset adds are booked as Min, Jin and (with negative sign) in
//!    the hydro / radiation energy-out registers of Ebnd.  The outer x1 boundary has no
//!    window then: outflow everywhere (v_r = max(v_r, 0); problem/outer_rho = copy gives
//!    a zero-gradient density and T, honoured with thermo = general), radiation Marshak
//!    with a dark (zero) bath, as the old top outside the window.  problem/env_r_top_max
//!    caps the WB top radius (set it to the old r_out: beyond L1 the Roche potential falls
//!    below the envelope's top cap again and would look like envelope).
//!  problem/donor_mask = roche: the donor's Roche lobe in the band (Phi < Phi(L1) and x >
//!    d_L1) is a sink: after the sources of every stage its cells are reset to the floor
//!    ambient (rho_amb, at rest in the rotating frame, T = the IC's ambient T_col); the
//!    removed mass and gas energy are booked in history columns Mmsk and Emsk (in place
//!    of Mwal, Jwal: identically 0 with the closed free-slip wall).  Radiation is not
//!    touched, unless problem/donor_irr = true: then the lobe also holds E = a T_d^4,
//!    F = 0 (T_d = problem/donor_teff, default 26200 K), a black-body surface whose light
//!    enters the grid by the M1 transport (booked in Ebnd as radiation in).
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
//! The cumulative integrals continue across a restart (acctrc-1019c: stored in the rst's
//! pgen block, RyPerRstWrite; a pre-1019c rst restarts them from zero).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
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
#include "driver/dt_retry.hpp"
#include "units/units.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
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
// accwall-1019 booking fix: the radiation the stage resets put into the M1 state (nozzle
// stream light zs[3], donor_irr / sector light zs[6]) as a PLAIN sum.  The M1 u0 is not
// RK-combined, so the RK-weighted booking in register 10 kept only half of the stage-1
// increment (EBUD-1019: the whole 'unbooked energy gain').  History only.
Real erinj_ = 0.0;
// <rad_m1> on the envelope (stage S3, accretor-rhd-1008; see the header RADIATION)
bool rad_ = false;
Real rad_fin_ = 0.0;           // imposed x1min radiative flux (code), = f_band F_col(r_in)
Real rad_lfac_ = 1.0;          // code band luminosity -> Lsun of the full sphere
Real rad_pn0_ = 0.0, rad_ps0_ = 0.0;   // Picard counters at the previous history output
Real picmax_rst_ = 0.0;        // Picmax carried by the restart file (acctrc-1019c)
// S5 (user 10-09): problem/env_amb_mode = radeq (floor ambient at radiative equilibrium,
// floor booked) and problem/stream_rad (the window's radiation injected); see the header
bool amb_radeq_ = false;
bool strad_ = false;
Real strad_e_ = 0.0;           // a T_s^4 of the window (code)
bool stprof_ = false;          // problem/stream_t_profile = diffusion (StreamT)
Real gen_tsurf_ = 0.0;         // problem/stream_t_surf (code T): the profile's floor
Real strad_ar_ = 0.0;          // <rad_m1> arad, for E = a T(phi)^4 in the window ghosts
DvceArray3D<Real> rbath_;      // per-column outer Marshak bath (m,k,j), handed to <rad_m1>
std::vector<Real> rcol_psi_, rcol_g_, rcol_f_;   // column psi, G(psi), F(psi) [code]
// accretor-bigdom-1015 (design sect. 13): the stream as a VOLUME source inside a larger
// domain (problem/stream_inject = volume), the donor's Roche lobe as a sink
// (problem/donor_mask = roche); see the header STREAM AS A VOLUME SOURCE
bool stvol_ = false;           // stream_inject = volume
Real str_r_ = 0.0;             // problem/stream_r: radius of the injection (default x1max)
Real nz_rlo_ = 0.0, nz_rhi_ = 0.0;   // the nozzle shell: cells with x1v in (rlo, rhi)
bool msk_ = false;             // donor_mask = roche
Real msk_phl1_ = 0.0, msk_xl1_ = 0.0;   // Phi(L1) and the L1 distance from the accretor
bool mirr_ = false;            // problem/donor_irr: the mask radiates at T_d
Real mirr_e_ = 0.0;            // a T_d^4 (code), x donor_irr_scale
// accretor-bigdom-v4 (10-16): donor_irr_scale (2-D band: auto = physical 3-D dilution),
// donor_gdark (von Zeipel beta: E = s a T_loc^4, T_loc^4 = C g_eff^(4 beta), floored),
// sink_zero_v (zero v in the sink cells after the M1 coupling)
Real mirr_sc_ = 1.0;
Real mirr_beta_ = 0.0;
Real mirr_ec_ = 0.0, mirr_ef_ = 0.0;   // s a C (code E per g^(4 beta)) and s a T_floor^4
bool szv_ = false;
bool stl1_ = false;            // stream_inject = l1 (the nozzle at L1, Lubow-Shu state)
int nz_if_ = -1;               // local x1 face index of the nozzle's inner face
Real l1_c_ = 0.0, l1_t_ = 0.0; // L1 nozzle: speed along -x (km/s) and uniform T (code)
// two extra history registers (columns 23/24, only with a nozzle or the mask, so older
// decks keep 22 columns): Mnzf = cumulative mass through the nozzle's inner face (RK flux
// integral, the stream that actually leaves the nozzle toward the star), Eirr = cumulative
// radiation energy put into the lobe by donor_irr
bool accx_on_ = false;
bool nzback_ = false;          // problem/stream_back_sink: the sector behind the nozzle is sink
// accretor-sector-light-1017: problem/stream_sector_light, the radiation held in the back-sink
// sector (stream_back_sink with stream_rad): 0 full (v4 default: the stream's own light,
// unscaled), 1 scaled (x donor_irr_scale), 2 surface (the donor lobe's surface light, as
// the lobe; v3 rule), 3 dark (left to the transport, as the ambient)
int ssl_ = 0;
bool flmom_ = false;           // problem/floor_keep_momentum (radeq floor: v = m/rho_amb, capped)
bool dtrep_ = false;           // problem/dt_report: print the hydro dt-setting cell at hst outputs
Real accx0_[2] = {0.0, 0.0}, accx1_[2] = {0.0, 0.0};
// accretor-3d-1016: a WIDE theta range (half-width about the orbital plane > 2 deg, e.g.
// +-20 deg) switches the nozzle, the donor mask and the nozzle-face register to 3-D: the
// nozzle is Gaussian in z = r cos(theta) with sigma_z (problem/stream_sigma_z, Rsun), its
// peak density normalised to the FULL stream rate stream_mdot over the 3-D window, and the
// mask is the exact 3-D Roche lobe.  The thin 2-D band (+-0.35 deg) never enters these
// branches (bitwise unchanged).  Forces, IC and WB stay those of the band (theta-invariant).
bool g3d_ = false;
Real sz_ = 0.0;                // sigma_z of the nozzle (Rsun)
// accretor-3d-1017 (problem/phys_3d = true, 3-D mode only, envelope + thermo = general; read
// only when given, so every 2-D deck and the 049c1f39 3-D benchmark deck run the old code):
// the PHYSICS in 3-D.  (1) every potential is the exact 3-D Roche potential RochePot3(r,
// theta, phi) (centrifugal about the z axis through the CM, the same function as the 3-D
// lobe mask): the true Phi (phicc0, phi0 on all three face sets, etotgrav), Phi_wb of the
// x1 WB pair (PhiWB3, the spin term with the cylindrical radius r sin theta), the explicit
// remainders g_r = -d(Phi - Phi_wb)/dr, g_theta = -dPhi/(r dtheta), g_phi = -dPhi/(r sin
// theta dphi) from face differences, all on (m,k,j,i); (2) the full rotating-frame Coriolis
// -2 Omega x v (a_r += 2 Omega v_phi sin th, a_th += 2 Omega v_phi cos th, a_phi -= 2 Omega
// (v_r sin th + v_th cos th); no work), the centrifugal term inside Phi (its work is in the
// etotgrav flux; no explicit centrifugal work: the deep_hot_jupiter_rt rot_potential +
// etotgrav double-count guard); (3) the envelope IC, the column T(psi), the radiation
// reference G(psi), the ambient rho_amb and the relaxation targets on the 3-D equipotentials
// psi = Phi_s - Phi_wb(r, theta, phi) (column by (m,k,j)); (4) problem/theta_bc = reflect
// (default; mesh ix2/ox2_bc = reflect) | outflow | hse (mesh ix2/ox2_bc = user: ghosts by
// the isothermal hydrostatic walk in theta in the true Phi at the edge cell's T; outflow:
// v_theta a diode (outflow only), v_r, v_phi copied; hse: a closed wall, v from the mirror
// cell with v_theta flipped; M1 ghosts reflecting, dark only with <rad_m1>/implicit_bc_x2min
// /max = marshak; the mass, Lz and energy through the theta faces join Mout, Jout and the
// hydro energy out), and problem/theta_sponge_deg / theta_sponge_tau (a velocity sponge in a
// layer next to each theta wall, linear ramp, outside the envelope, KE removed from E);
// (5) the inertial Lz bookkeeping with the cylindrical radius r sin theta.
bool p3_ = false;
bool thout_ = false;           // theta_bc = outflow | hse (user theta faces of this pgen)
bool thhse_ = false;           // theta_bc = hse: closed hydrostatic wall (v_theta mirrored)
// accwall-1019: problem/theta_wall_closed = true (theta_bc = hse only; read only when named,
// default false = bitwise unchanged): the hse theta walls are EXACTLY closed.  After the
// whole hydro flux task (RyPerWallFlux, the Hydro::Fluxes hook) the x2 face fluxes at the
// domain's theta walls carry no mass, energy (incl. the etotgrav rho Phi v term), x1 or x3
// momentum or scalars; only the x2 momentum flux is kept = the Riemann solver's interface
// pressure between the edge cell and its hydrostatic ghost (the hse ghosts exist to make
// exactly this the column's hydrostatic wall pressure; it also holds the solver's
// rho v_n^2 part of the old leak, O(v_n^2) with v_n ~ 1e-3 km/s).  Without it the hse
// ghost (hydrostatic density, mirrored velocity) and the edge cell form a non-mirror
// Riemann problem whose mass flux leaked ~2x the stream rate INTO the star (EBUD-1019).
bool thwc_ = false;
Real thsp_w_ = 0.0, thsp_tau_ = 0.0;   // theta sponge width (rad), time (code)
DvceArray4D<Real> gr4_, gt4_, gp4_, gpl4_;   // 3-D: g_r remainder, g_theta, g_phi, plain g_r
DvceArray4D<Real> cref24_, rhoamb4_, tcol4_; // 3-D: relaxation target, ambient, column T
// <rad_m1>/implicit_bc_x2min|max = flux (phys_3d): the imposed equilibrium lateral flux at the
// theta walls (RadSetup) and its total outward rate on this rank (booked as radiation out)
DvceArray3D<Real> fx2lo_, fx2hi_;
Real fx2leak_ = 0.0;
// acctrc-1019 (USER 10-19): two passive advective tracers, <hydro>/nscalars = 2 (0 = off: no
// branch below runs, every kernel and register of the default is unchanged).  Conserved
// u(nhydro + n) = rho s_n, advected by the core with the upwind mass flux (hydro_fluxes.cpp).
//   s0 STREAM:          1 in the gas the nozzle writes (stream_inject = volume; the r_out
//                       window ghosts with boundary injection), 0 elsewhere at the start.
//   s1 BOUNDARY INFLOW: 1 in all gas that enters through any physical boundary of the mesh
//                       other than the nozzle (r_in, r_out, the theta walls).
// IC envelope / ambient 0.  Mass added by any floor (EOS dfloor, radeq floor) carries s = 0:
// the conserved rho s is left as it is, so the concentration is diluted by the added mass.
// The labelling of boundary inflow is EXACT: after the hydro update of every stage the
// scalar flux through each user boundary face is replaced by F_rho (0, 1) where F_rho points
// INTO the domain (TrcFluxFix; the ghosts themselves carry the zero-gradient / mirror
// concentration, so the interior slopes and the outflow are not biased).  Diagnostics: the
// separate file <basename>.trc.hst (ryper.user.hst is unchanged); registers per rank, RK
// weighted like acc0_, cumulative, carried across restarts (1019c), in the dt_retry checkpoint.
bool trc_ = false;
int trc_n0_ = 5;               // index of s0 in u0 (= Hydro::nhydro); s1 = trc_n0_ + 1
constexpr int kNtr = 9;
// 0 Nz0  s0 mass written by the nozzle reset (dn - (rho s0)_old) dV
// 1 Nz1  s1 mass removed by the nozzle reset ((rho s1)_old dV)
// 2 Sk0  s0 mass removed by the sinks (donor mask, back sink: reset to s = 0)
// 3 Sk1  s1 mass removed by the sinks
// 4 Fi0, 5 Fi1  s0, s1 mass OUT through r_in (outward positive)
// 6 Fo0, 7 Fo1  s0, s1 mass OUT through r_out and the user theta walls (outward positive)
// 8 Bin  mass that ENTERED through the non-nozzle user boundaries (= the s1 source)
Real trc0_[kNtr] = {0.0}, trc1_[kNtr] = {0.0};
std::string trc_file_;

KOKKOS_INLINE_FUNCTION
Real RochePot(const RocheParams &p, const Real r, const Real phi) {
  const Real cp = cos(phi);
  const Real rd = sqrt(fmax(r*r + p.asep*p.asep - 2.0*p.asep*r*cp, 1.0e-30));
  return -p.gma/r - p.gmd/rd
         - 0.5*SQR(p.omega)*(r*r + p.xcm*p.xcm - 2.0*p.xcm*r*cp);
}

//! accretor-3d-1016: the exact 3-D Roche potential at (r, theta, phi), sth = sin(theta)
//! (x = r sth cos phi, y = r sth sin phi, z = r cos theta; centrifugal about the z axis
//! through the CM).  Used only by the 3-D donor mask.
KOKKOS_INLINE_FUNCTION
Real RochePot3(const RocheParams &p, const Real r, const Real sth, const Real phi) {
  const Real xs = r*sth*cos(phi);
  const Real rd = sqrt(fmax(r*r + p.asep*p.asep - 2.0*p.asep*xs, 1.0e-30));
  return -p.gma/r - p.gmd/rd
         - 0.5*SQR(p.omega)*(r*r*sth*sth + p.xcm*p.xcm - 2.0*p.xcm*xs);
}

//! accretor-bigdom-v4 (10-16): the Roche potential and |grad Phi| (effective gravity) at a
//! Cartesian point (accretor at the origin, donor at (a, 0, 0), centrifugal about the z
//! axis through the CM at (xcm, 0, 0)), and the donor's von Zeipel gravity: |grad Phi| at
//! the point where the ray from the donor centre through (x, y, z) leaves the L1 lobe
//! (Phi = Phi_L1).  A ray that never leaves it (exactly through the L1 neck) gives 0.
KOKKOS_INLINE_FUNCTION
Real RochePotXYZ(const RocheParams &p, const Real x, const Real y, const Real z) {
  const Real r = sqrt(fmax(x*x + y*y + z*z, 1.0e-30));
  const Real rd = sqrt(fmax(SQR(x - p.asep) + y*y + z*z, 1.0e-30));
  return -p.gma/r - p.gmd/rd - 0.5*SQR(p.omega)*(SQR(x - p.xcm) + y*y);
}

KOKKOS_INLINE_FUNCTION
Real RocheGeffXYZ(const RocheParams &p, const Real x, const Real y, const Real z,
                  Real *gx = nullptr, Real *gy = nullptr, Real *gz = nullptr) {
  const Real r = sqrt(fmax(x*x + y*y + z*z, 1.0e-30));
  const Real rd = sqrt(fmax(SQR(x - p.asep) + y*y + z*z, 1.0e-30));
  const Real ar = p.gma/(r*r*r), ad = p.gmd/(rd*rd*rd), o2 = SQR(p.omega);
  const Real dx = ar*x + ad*(x - p.asep) - o2*(x - p.xcm);
  const Real dy = ar*y + ad*y - o2*y;
  const Real dz = ar*z + ad*z;
  if (gx != nullptr) {*gx = dx; *gy = dy; *gz = dz;}
  return sqrt(dx*dx + dy*dy + dz*dz);
}

//! surface point along the unit ray (nx, ny, nz) from the donor centre: march out in steps
//! of 0.005 a from t0, bisect the first crossing of Phi_L1.  Returns t (0 = none).  The
//! march stops at the L1 distance a - x_L1 (+0.1 %), the farthest point of the lobe: a ray
//! that slips through the thin neck would otherwise find the gainer side
KOKKOS_INLINE_FUNCTION
Real DonorSurfT(const RocheParams &p, const Real phl1, const Real xl1, const Real nx,
                const Real ny, const Real nz, const Real t0) {
  const Real tmx = 1.001*(p.asep - xl1);
  const Real st = 0.005*p.asep;
  Real ta = t0, tb = -1.0;
  for (int n=0; n<400; ++n) {
    const Real t = ta + st;
    if (t > tmx) break;
    if (RochePotXYZ(p, p.asep + t*nx, t*ny, t*nz) > phl1) {tb = t; break;}
    ta = t;
  }
  if (tb < 0.0) return 0.0;
  for (int n=0; n<40; ++n) {
    const Real tm = 0.5*(ta + tb);
    if (RochePotXYZ(p, p.asep + tm*nx, tm*ny, tm*nz) > phl1) tb = tm; else ta = tm;
  }
  return 0.5*(ta + tb);
}

KOKKOS_INLINE_FUNCTION
Real DonorSurfGeff(const RocheParams &p, const Real phl1, const Real xl1, const Real x,
                   const Real y, const Real z) {
  Real nx = x - p.asep, ny = y, nz = z;
  Real t0 = sqrt(nx*nx + ny*ny + nz*nz);
  if (t0 < 1.0e-6*p.asep) {nx = 0.0; ny = 1.0; nz = 0.0; t0 = 0.0;} else {
    nx /= t0; ny /= t0; nz /= t0;
  }
  // a cell centre a little outside the lobe surface (the mask uses its own Phi test) is
  // started from the centre
  if (RochePotXYZ(p, x, y, z) > phl1) t0 = 0.0;
  const Real t = DonorSurfT(p, phl1, xl1, nx, ny, nz, t0);
  if (t <= 0.0) return 0.0;
  return RocheGeffXYZ(p, p.asep + t*nx, t*ny, t*nz);
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

//! accretor-3d-1017 (phys_3d): PhiWB in 3-D, sth = sin(theta): the exact 3-D Roche
//! potential capped at r_top and at the top equipotential, the spin term with the
//! cylindrical radius min(r, r_spin) sin(theta)
KOKKOS_INLINE_FUNCTION
Real PhiWB3(const RocheParams &p, const Real dlt, const Real rtop, const Real rspin,
            const Real phtop, const Real r, const Real sth, const Real phi) {
  const Real x = fmin(r, rtop), xs = fmin(x, rspin)*sth;
  return fmin(RochePot3(p, x, sth, phi) - dlt*xs*xs, phtop);
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

//! problem/stream_t_profile = diffusion (S5, user 10-09): the stream's T across the window.
//! The core T_c (env_t_stream, the one-zone value) falls toward the stream surface as the
//! lowest diffusion (cooling) mode of a slab in optical-depth coordinate, T^4 = T_c^4
//! cos(pi/2 tau(x)/tau_half), with tau(x)/tau_half = erf(|x|/sqrt 2) for the Gaussian
//! density at constant opacity (x = dphi/sigma), clamped below at the thin
//! radiative-equilibrium T_eq (problem/stream_t_surf).  On = false: T_c everywhere.
KOKKOS_INLINE_FUNCTION
Real StreamT(const bool on, const Real tc, const Real teq, const Real x) {
  if (!on) return tc;
  const Real c = cos(0.5*M_PI*erf(fabs(x)/sqrt(2.0)));
  return fmax(tc*sqrt(sqrt(fmax(c, 0.0))), teq);
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
  const bool stvol = stvol_, ohse = outer_hse_ || !gen_;
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
  const bool strad = strad_;
  const Real sre = strad_e_;
  const bool stp = stprof_;
  const Real tsf = gen_tsurf_, sar = strad_ar_;
  auto pwcc = pmbp->phydro->phicc_wb;
  DvceArray5D<Real> ur;
  Real rcl = 1.0, refl = 0.0;
  if (rad) {
    ur = pmbp->pradm1->u0;
    rcl = pmbp->pradm1->c_light;
    refl = pmbp->pradm1->e_floor;
  }
  // accretor-3d-1017 (phys_3d): the 3-D potentials and the cylindrical spin radius
  const bool p3 = p3_;
  auto &x2v = pmbp->pcoord->x2v;
  par_for("ryper_bce", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const Real ph = x3v(m,k);
    const Real sth = p3 ? sin(x2v(m,j)) : 1.0;
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
        const Real pa = rad ? pwcc(m,k,j,is) :
                        (p3 ? PhiWB3(p, dlt, rtop, rspin, phtop, x1v(m,is), sth, ph)
                            : PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,is), ph));
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
                            (p3 ? PhiWB3(p, dlt, rtop, rspin, phtop, x1v(m,i1), sth, ph)
                                : PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,i1), ph))) - pa;
          if (e1 > 0.0 && dp1 != 0.0) a = (eosd.Temperature(d1, e1, ta) - ta)/dp1;
        }
        for (int g=0; g<ng; ++g) {
          const int ig = is - 1 - g, im = is + g;
          const Real pg = rad ? pwcc(m,k,j,ig) :
                          (p3 ? PhiWB3(p, dlt, rtop, rspin, phtop, x1v(m,ig), sth, ph)
                              : PhiWB(p, dlt, rtop, rspin, phtop, x1v(m,ig), ph));
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
          const Real v3 = noslip ? vwall*x1v(m,ig)*sth : u0(m,IM3,k,j,im)/dm;
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
        const bool win = strm && !stvol && (fabs(dph) < wwin);
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
        const Real pa = p3 ? RochePot3(p, x1v(m,ie), sth, ph) : RochePot(p, x1v(m,ie), ph);
        for (int g=0; g<ng; ++g) {
          const int ig = ie + 1 + g;
          Real dg, v1, v2, v3, tg;
          if (win) {
            dg = rhos*exp(-0.5*SQR(dph/sig));
            v1 = vrs; v2 = 0.0; v3 = vps; tg = StreamT(stp, tstr, tsf, dph/sig);
          } else {
            dg = da;
            if (ohse) {
              const Real pg = p3 ? RochePot3(p, x1v(m,ig), sth, ph) : RochePot(p, x1v(m,ig), ph);
              Real eg = ea, tw = ta;
              WBAdvance(eosd, 1, da, ea, pg - pa, dg, eg, tw, ta, 0.0, ta, ta);
              if (!(Kokkos::isfinite(dg) && dg > 0.0)) dg = da;
              dg = fmin(dg, da*exp(cap));
            }
            v1 = fmax(v1a, 0.0); v2 = v2a; v3 = v3a; tg = ta;
          }
          dg = fmax(dg, dfl);
          u0(m,IDN,k,j,ig) = dg;
          u0(m,IM1,k,j,ig) = dg*v1;
          u0(m,IM2,k,j,ig) = dg*v2;
          u0(m,IM3,k,j,ig) = dg*v3;
          u0(m,IEN,k,j,ig) = eosd.EnergyFromTemperature(dg, tg)
                             + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
          if (rad) {
            if (win && strad) {
              // problem/stream_rad: the stream's radiation in equilibrium with its gas,
              // E = a T_s^4, lab-frame flux of optically thick advection (4/3) E v
              const Real eg = stp ? sar*SQR(SQR(tg)) : sre;
              ur(m,radm1::M1_E,k,j,ig) = eg;
              ur(m,radm1::M1_F1,k,j,ig) = (4.0/3.0)*eg*v1;
              ur(m,radm1::M1_F2,k,j,ig) = 0.0;
              ur(m,radm1::M1_F3,k,j,ig) = (4.0/3.0)*eg*v3;
            } else {
              radm1::M1FillGhost(ur, m, k, j, ig, k, j, ie, 1, 2, 1.0, rcl, refl);
            }
          }
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
  // acctrc-1019: the tracer ghosts of the x1 faces, from the ghost densities just written.
  // r_in (closed wall, the face states are mirrored by hydro wall_closed_ix1): the mirror
  // cell's concentrations; r_out: the edge cell's (zero gradient: an unbiased outflow);
  // the r_out stream window (boundary injection, stream_inject != volume): s0 = 1, s1 = 0.
  // Inflow through these faces is labelled s1 = 1 exactly by TrcFluxFix, not here.
  if (trc_) {
    const int n0 = trc_n0_;
    par_for("ryper_bctr1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int g=0; g<ng; ++g) {
          const int ig = is - 1 - g, im = is + g;
          const Real dm = u0(m,IDN,k,j,im), dg = u0(m,IDN,k,j,ig);
          for (int n=n0; n<n0+2; ++n) {
            u0(m,n,k,j,ig) = (dm > 0.0) ? dg*u0(m,n,k,j,im)/dm : 0.0;
          }
        }
      }
      if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
        const bool win = strm && !stvol && (fabs(WrapPhi(x3v(m,k) - phis)) < wwin);
        const Real da = u0(m,IDN,k,j,ie);
        for (int g=0; g<ng; ++g) {
          const int ig = ie + 1 + g;
          const Real dg = u0(m,IDN,k,j,ig);
          if (win) {
            u0(m,n0,k,j,ig) = dg;
            u0(m,n0+1,k,j,ig) = 0.0;
          } else {
            for (int n=n0; n<n0+2; ++n) {
              u0(m,n,k,j,ig) = (da > 0.0) ? dg*u0(m,n,k,j,ie)/da : 0.0;
            }
          }
        }
      }
    });
  }
  // accretor-3d-1017: problem/theta_bc = outflow | hse (mesh ix2/ox2_bc = user): the theta
  // ghosts continue the edge cell isothermally (its T) and hydrostatically in the TRUE 3-D
  // Roche potential along theta (the outer-x1 hse walk of the general EOS, capped at
  // e^hse_cap), so the wall face carries the hydrostatic pressure (a plain reflecting wall
  // gives the face the edge cell's own pressure and leaves the theta component of the tidal
  // and centrifugal force unbalanced there: v_theta ~ 7 km/s in the deep envelope at the
  // walls, Raven s3dA 10-17).  outflow: v_theta a diode (outflow only), v_r, v_phi copied
  // from the edge cell; hse: a closed wall, v (all three, v_theta sign-flipped) from the
  // mirror cell.  Radiation ghosts: reflecting (as the implicit default for these faces;
  // the star's interior reaches the walls, a dark face there drains it), dark only with
  // <rad_m1>/implicit_bc_x2min/max = marshak.
  if (!(p3 && thout_ && gen)) return;
  const int js = indcs.js, je = indcs.je;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const bool thh = thhse_;
  int rm2lo = 1, rm2hi = 1;
  if (rad) {
    // acc3d-1017: implicit_bc_x2 = outflow (zero-gradient outgoing flux) -> copied ghosts
    auto rmode = [](const int b) {
      return (b == radm1::M1_IBC_MARSHAK) ? 2 : ((b == radm1::M1_IBC_OUTFLOW) ? 0 : 1);
    };
    rm2lo = rmode(pmbp->pradm1->ibc_x2min);
    rm2hi = rmode(pmbp->pradm1->ibc_x2max);
  }
  par_for("ryper_bcth", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    const Real ph = x3v(m,k), r = x1v(m,i);
    for (int sd=0; sd<2; ++sd) {
      if (sd == 0 && mbbcs.d_view(m, BoundaryFace::inner_x2) != BoundaryFlag::user) continue;
      if (sd == 1 && mbbcs.d_view(m, BoundaryFace::outer_x2) != BoundaryFlag::user) continue;
      const int ja = (sd == 0) ? js : je;
      const Real da = u0(m,IDN,k,ja,i);
      const Real v1a = u0(m,IM1,k,ja,i)/da;
      const Real v2a = u0(m,IM2,k,ja,i)/da;
      const Real v3a = u0(m,IM3,k,ja,i)/da;
      Real ea = u0(m,IEN,k,ja,i) - 0.5*da*(v1a*v1a + v2a*v2a + v3a*v3a)
                - da*phicc(m,k,ja,i);
      Real ta = (ea > 0.0) ? eosd.Temperature(da, ea) : tfl;
      if (!(ta > tfl)) {
        ta = tfl;
        ea = eosd.EnergyFromTemperature(da, ta);
      }
      // v3: the walk follows the NET potential of the column's balance, Phi + G(psi) (the
      // true Phi plus the radiation-force reference G = Phi_eff - Phi_wb; G = 0 without
      // <rad_m1>, constant above the envelope top), with the local polytrope of the edge
      // cell and its inner neighbour along theta (T linear in that potential, the inner-x1
      // WB ghost rule; isothermal fallback), so the ghosts lie on the column's own P(psi),
      // T(psi) and the wall face carries its hydrostatic gas pressure
      auto pnet = [&](const int jj) {
        const Real sj = sin(x2v(m,jj));
        const Real gref = rad ? (pwcc(m,k,jj,i) - PhiWB3(p, dlt, rtop, rspin, phtop, r, sj, ph))
                              : 0.0;
        return RochePot3(p, r, sj, ph) + gref;
      };
      const Real pa = pnet(ja);
      Real aT = 0.0;
      {
        const int j1 = (sd == 0) ? js + 1 : je - 1;
        const Real d1 = u0(m,IDN,k,j1,i);
        const Real e1 = u0(m,IEN,k,j1,i) - 0.5*(SQR(u0(m,IM1,k,j1,i)) + SQR(u0(m,IM2,k,j1,i))
                        + SQR(u0(m,IM3,k,j1,i)))/d1 - d1*phicc(m,k,j1,i);
        const Real dp1 = pnet(j1) - pa;
        if (e1 > 0.0 && d1 > 0.0 && dp1 != 0.0) aT = (eosd.Temperature(d1, e1, ta) - ta)/dp1;
      }
      for (int g=0; g<ng; ++g) {
        const int jg = (sd == 0) ? js - 1 - g : je + 1 + g;
        const int jm = (sd == 0) ? js + g : je - g;     // the mirror cell
        const Real pg = pnet(jg);
        Real dg = da, eg = ea, tw = ta;
        WBAdvance(eosd, 3, da, ea, pg - pa, dg, eg, tw, ta, aT, ta, ta);
        if (!(Kokkos::isfinite(dg) && dg > 0.0 && Kokkos::isfinite(eg) && eg > 0.0 &&
              tw >= tfl)) {
          dg = da; eg = ea; tw = ta;
          WBAdvance(eosd, 1, da, ea, pg - pa, dg, eg, tw, ta, 0.0, ta, ta);
          if (!(Kokkos::isfinite(dg) && dg > 0.0)) dg = da;
          tw = ta;
        }
        dg = fmax(fmin(dg, da*exp(cap)), dfl);
        Real v1 = v1a, v2 = (sd == 0) ? fmin(v2a, 0.0) : fmax(v2a, 0.0), v3 = v3a;
        if (thh) {
          const Real dm = u0(m,IDN,k,jm,i);
          v1 = u0(m,IM1,k,jm,i)/dm;
          v2 = -u0(m,IM2,k,jm,i)/dm;
          v3 = u0(m,IM3,k,jm,i)/dm;
        }
        u0(m,IDN,k,jg,i) = dg;
        u0(m,IM1,k,jg,i) = dg*v1;
        u0(m,IM2,k,jg,i) = dg*v2;
        u0(m,IM3,k,jg,i) = dg*v3;
        u0(m,IEN,k,jg,i) = eosd.EnergyFromTemperature(dg, tw)
                           + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,jg,i);
        if (rad) {
          const int md = (sd == 0) ? rm2lo : rm2hi;
          radm1::M1FillGhost(ur, m, k, jg, i, k, (md == 1) ? jm : ja, i, 2, md,
                             (sd == 0) ? -1.0 : 1.0, rcl, refl);
        }
      }
    }
  });
  // acctrc-1019: the tracer ghosts of the user theta walls: hse (closed wall) the mirror
  // cell's concentrations, outflow the edge cell's; inflow is labelled by TrcFluxFix
  if (trc_) {
    const int n0 = trc_n0_;
    par_for("ryper_bctr2", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int i) {
      for (int sd=0; sd<2; ++sd) {
        if (sd == 0 && mbbcs.d_view(m, BoundaryFace::inner_x2) != BoundaryFlag::user) continue;
        if (sd == 1 && mbbcs.d_view(m, BoundaryFace::outer_x2) != BoundaryFlag::user) continue;
        const int ja = (sd == 0) ? js : je;
        for (int g=0; g<ng; ++g) {
          const int jg = (sd == 0) ? js - 1 - g : je + 1 + g;
          const int jr = thh ? ((sd == 0) ? js + g : je - g) : ja;
          const Real dr = u0(m,IDN,k,jr,i), dg = u0(m,IDN,k,jg,i);
          for (int n=n0; n<n0+2; ++n) {
            u0(m,n,k,jg,i) = (dr > 0.0) ? dg*u0(m,n,k,jr,i)/dr : 0.0;
          }
        }
      }
    });
  }
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
//----------------------------------------------------------------------------------------
//! accretor-bigdom-v4 (10-16): the sink cells (the donor lobe and the back-sink sector) are
//! reset to floor gas AT REST after the other sources of every stage, but the M1 implicit
//! coupling that follows gives that 1e-8 floor gas the radiation momentum of the light it
//! holds (Raven 10-16: v_r 2e5-1e6 km/s in the first back-sink row behind the nozzle, in
//! every output and in Emsk).  This hook (RadiationM1::HydroConToPrim, after the coupling,
//! before the fluid refresh) zeroes their momentum again and books the kinetic energy it
//! removes in Emsk.  The hydro dt is unchanged (computed before the M1 stage, and these
//! cells never set it).  Active only while a sink can hold light or follow the stream
//! (donor_irr on, or the back sink with the stream on); problem/sink_zero_v = false: off.
//! accwall-1019 (problem/theta_wall_closed): the hse theta walls exactly closed.  Called
//! by Hydro::Fluxes after every flux contribution (Riemann + FOFC, etotgrav, WB removal),
//! so the RK update, the spherical-polar L_z source (it reads x2 IM3 fluxes), the
//! <hydro>/am_transport module and the boundary registers of RyPerSrcEnv all see zeros:
//! every x2 flux component except IM2 is set to 0 on the domain-wall faces (j = js of a
//! block with a user inner_x2 face, j = je + 1 with a user outer_x2 face), over the whole
//! k, i extent of the flux array.
void RyPerWallFlux(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int js = indcs.js, je = indcs.je;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto flx2 = pmbp->phydro->uflx.x2f;
  const int nv = flx2.extent_int(1), n3 = flx2.extent_int(2), n1 = flx2.extent_int(4);
  auto &mbbcs = pmbp->pmb->mb_bcs;
  par_for("ryper_wallflx", DevExeSpace(), 0, nmb1, 0, n3 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int i) {
    for (int sd=0; sd<2; ++sd) {
      if (sd == 0 && mbbcs.d_view(m, BoundaryFace::inner_x2) != BoundaryFlag::user) continue;
      if (sd == 1 && mbbcs.d_view(m, BoundaryFace::outer_x2) != BoundaryFlag::user) continue;
      const int jf = (sd == 0) ? js : je + 1;
      for (int n=0; n<nv; ++n) {
        if (n != IM2) flx2(m,n,k,jf,i) = 0.0;
      }
    }
  });
}

void RyPerM1Post(Mesh *pm) {
  if (!szv_ || !msk_) return;
  if (!(mirr_ || (nzback_ && stvol_ && stream_on_))) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto u0 = pmbp->phydro->u0;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x2c = pmbp->pcoord->x2v;
  auto &x3c = pmbp->pcoord->x3v;
  auto &volume = pmbp->pcoord->volume;
  const bool nzon = stvol_ && stream_on_, nzbk = nzback_ && nzon, g3 = g3d_;
  const Real nlo = nz_rlo_, nhi = nz_rhi_, phs = phi_s_, sg = sig_, ww = nsig_*sig_;
  const Real szz = sz_, wz = nsig_*sz_;
  const Real phl1 = msk_phl1_, xl1 = msk_xl1_;
  const RocheParams pr = rp_;
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Real eke = 0.0;
  Kokkos::parallel_reduce("ryper_m1post", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &sum) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
    const int i = idx%ni + is;
    const Real r = x1v(m,i), phv = x3c(m,k);
    const Real dph = WrapPhi(phv - phs);
    const Real zc = g3 ? nlo*cos(x2c(m,j)) : 0.0;
    const Real sth = g3 ? sin(x2c(m,j)) : 1.0;
    const bool innz = nzon && r > nlo && r < nhi && fabs(dph) < ww && (!g3 || fabs(zc) < wz);
    if (innz) return;
    bool inmk;
    if (g3) {
      inmk = (r*sth*cos(phv) > xl1) && (RochePot3(pr, r, sth, phv) < phl1);
    } else {
      inmk = (r*cos(phv) > xl1) && (RochePot(pr, r, phv) < phl1);
    }
    if (nzbk && r > nhi && fabs(dph) < ww + sg && (!g3 || fabs(zc) < wz + szz)) inmk = true;
    if (!inmk) return;
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
    if (!(d > 0.0)) return;
    const Real ke = 0.5*(m1*m1 + m2*m2 + m3*m3)/d;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) -= ke;
    sum += ke*volume(m,k,j,i);
  }, Kokkos::Sum<Real>(eke));
  acc0_[3] += eke;
}

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
  // accretor-3d-1017 (phys_3d): the (m,k,j,i) forces/targets, full Coriolis, theta sponge
  const bool p3 = p3_;
  auto gr4 = gr4_, gt4 = gt4_, gp4 = gp4_, gpl4 = gpl4_;
  auto c2r4 = cref24_, ramb4 = rhoamb4_, tcol4 = tcol4_;
  auto &x2v = pmbp->pcoord->x2v;
  const Real thsw = thsp_w_, thtau = thsp_tau_;
  const Real th0 = pm->mesh_size.x2min, th1 = pm->mesh_size.x2max;
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
      fg = d*(p3 ? gpl4(m,k,j,i) : gpl(m,k,i));
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
    if (p3) {
      // 3-D: gravity remainders on (m,k,j,i), the full Coriolis -2 Omega x v
      const Real vt = w0(m,IVY,k,j,i);
      const Real st = sin(x2v(m,j)), cth = cos(x2v(m,j));
      u0(m,IM1,k,j,i) += bdt*(fg + d*(gr4(m,k,j,i) + om2*vp*st));
      u0(m,IM2,k,j,i) += bdt*d*(gt4(m,k,j,i) + om2*vp*cth);
      u0(m,IM3,k,j,i) += bdt*d*(gp4(m,k,j,i) - om2*(vr*st + vt*cth));
    } else {
      u0(m,IM1,k,j,i) += bdt*(fg + d*(gr(m,k,i) + om2*vp));
      u0(m,IM3,k,j,i) += bdt*d*(gp(m,k,i) - om2*vr);
    }
    const Real rambc = p3 ? ramb4(m,k,j,i) : ramb(m,k,i);
    const Real tcolc = p3 ? tcol4(m,k,j,i) : tcol(m,k,i);
    // relaxation target: the envelope profile inside; outside, material by density
    // relative to the hydrostatic ambient: denser than env_amb_k x rho_amb(r,phi) is
    // stream/atmosphere gas (c_ph^2), the rest is the hot ambient (c_amb^2).  Choosing
    // by temperature or by an absolute density instead heated the stream's thin
    // leading edge and the atmosphere's top to c_amb (precursor jets at 1700 km/s and
    // an evaporation wind; tests g_half, a2, 10-06/07).
    Real ct = p3 ? c2r4(m,k,j,i) : c2r(m,k,i);
    // Time constants: env_t_relax_env inside, env_t_relax_stream for the stream/atmosphere
    // gas, env_t_relax for the hot ambient; t <= 0 switches the relaxation off (adiabatic).
    const bool inside = (ct > cph2*(1.0 + 1.0e-9));
    Real tr = tri;
    if (!inside) {
      const bool dense = (d > kamb*rambc);
      ct = dense ? ((x1v(m,i) <= rhot) ? cph2 : cstr2) : camb2;
      tr = dense ? trs : tro;
    }
    // local cooling time: t = e_gas/(c a T^4) (1/(kappa_abs rho) + 3 kappa rho H^2): thin
    // emission by the absorption (Kramers) opacity, diffusion by ES + Kramers; H = (P/rho)/g,
    // g = G M_a / r^2 (cgs inside, converted to code time); not for the hot ambient
    if (cool && (inside || d > kamb*rambc) && e > 0.0 && d > 0.0) {
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
          tt = tcolc;
        } else if (d > kamb*rambc) {
          tt = (x1v(m,i) <= rhot) ? tcolc : tstr;
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
    // accretor-3d-1017: theta sponge (problem/theta_sponge_deg > 0): outside the envelope,
    // within the layer next to each theta wall, the rotating-frame velocity decays at a rate
    // rising linearly from 0 at the layer's inner edge to 1/theta_sponge_tau at the wall;
    // the kinetic energy removed leaves E (booked in Erel with the floor sponge)
    if (p3 && thsw > 0.0 && !inside) {
      const Real dw = fmin(x2v(m,j) - th0, th1 - x2v(m,j));
      const Real dd = u0(m,IDN,k,j,i);
      if (dw < thsw && dd > 0.0) {
        const Real f = exp(-bdt*(1.0 - dw/thsw)/thtau);
        const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
        u0(m,IM1,k,j,i) = f*m1;
        u0(m,IM2,k,j,i) = f*m2;
        u0(m,IM3,k,j,i) = f*m3;
        u0(m,IEN,k,j,i) -= (1.0 - f*f)*0.5*(m1*m1 + m2*m2 + m3*m3)/dd;
      }
    }
  });
  // problem/env_amb_mode = radeq: the density floor applied here (the EOS floor then never
  // acts), added gas at the cell's velocity and temperature; its mass is booked (Mfl), its
  // energy is inside the Erel difference below
  Real dmfl = 0.0;
  if (amb_radeq_) {
    auto phc = ph->phicc0;
    const Real tfl = gen_tfl_;
    const bool fmom = flmom_;
    const Real vcf = env_vcap_;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("ryper_flr", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &sum) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = idx%ni + is;
      const Real du = u0(m,IDN,k,j,i);
      if (!(du < rfl)) return;
      const Real dw = w0(m,IDN,k,j,i), ew = w0(m,IEN,k,j,i);
      Real v1 = w0(m,IVX,k,j,i), v2 = w0(m,IVY,k,j,i), v3 = w0(m,IVZ,k,j,i);
      Real tc = (dw > 0.0 && ew > 0.0) ? eosd.Temperature(dw, ew) : tfl;
      if (du > 0.0) {
        v1 = u0(m,IM1,k,j,i)/du; v2 = u0(m,IM2,k,j,i)/du; v3 = u0(m,IM3,k,j,i)/du;
        const Real ei = u0(m,IEN,k,j,i) - 0.5*du*(v1*v1 + v2*v2 + v3*v3)
                        - du*phc(m,k,j,i);
        if (ei > 0.0) tc = eosd.Temperature(du, ei, tc);
      }
      if (fmom) {
        // bigdom-1015 (nozzle/mask runs only): the added gas is at REST, i.e. momentum is
        // kept (v = m/rho_amb), and |v| is capped at env_v_floor_max.  Next to the dense
        // L1 nozzle a near-vacuum cell can lose mass to the HLLC flux while the pressure
        // term gives it momentum: m/du then reached 7e5 km/s (dt collapse, g7 10-15)
        v1 = u0(m,IM1,k,j,i)/rfl; v2 = u0(m,IM2,k,j,i)/rfl; v3 = u0(m,IM3,k,j,i)/rfl;
        const Real vv = sqrt(v1*v1 + v2*v2 + v3*v3);
        if (vv > vcf) {
          v1 *= vcf/vv; v2 *= vcf/vv; v3 *= vcf/vv;
        }
      }
      tc = fmax(tc, tfl);
      u0(m,IDN,k,j,i) = rfl;
      u0(m,IM1,k,j,i) = rfl*v1;
      u0(m,IM2,k,j,i) = rfl*v2;
      u0(m,IM3,k,j,i) = rfl*v3;
      u0(m,IEN,k,j,i) = eosd.EnergyFromTemperature(rfl, tc)
                        + 0.5*rfl*(v1*v1 + v2*v2 + v3*v3) + rfl*phc(m,k,j,i);
      sum += (rfl - du)*volume(m,k,j,i);
    }, Kokkos::Sum<Real>(dmfl));
  }
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

  // acctrc-1019: TrcFluxFix -- the exact boundary-inflow label.  The core advected the
  // scalars with the upwind flux F_s = F_rho s_upwind; at every user boundary face where
  // F_rho points INTO the domain the scalar flux is replaced by F_rho (s0, s1) = F_rho (0, 1)
  // (the r_out stream window of boundary injection: (1, 0)), and the edge cell gets the
  // difference with this stage's weight, -bdt (F_new - F_old) A / V on the divergence of
  // the face (hydro_update.cpp form), i.e. exactly the update the new flux would have made.
  // Outflow faces keep the core's upwind flux.  The corrected fluxes are summed (per rank)
  // into the registers Fi*, Fo*, Bin below.
  array_sum::GlobalSum tfs;
  for (int n=0; n<NREDUCTION_VARIABLES; ++n) tfs.the_array[n] = 0.0;
  if (trc_) {
    const int n0 = trc_n0_;
    auto &tf1 = ph->uflx.x1f;
    auto &tf2 = ph->uflx.x2f;
    auto &ta2 = pmbp->pcoord->area.x2f;
    auto &tx3 = pmbp->pcoord->x3v;
    auto &tbcs = pmbp->pmb->mb_bcs;
    const bool twin = stream_on_ && !stvol_;
    const Real tphs = phi_s_, tww = nsig_*sig_;
    const int tnj = je - js + 1, tnk = ke - ks + 1, tni = ie - is + 1;
    array_sum::GlobalSum t1, t2;
    Kokkos::parallel_reduce("ryper_trfx1", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*tnk*tnj),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
      const int m = idx/(tnk*tnj);
      const int k = (idx - m*tnk*tnj)/tnj + ks;
      const int j = idx%tnj + js;
      array_sum::GlobalSum v;
      for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
      if (tbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        const Real a = area1(m,k,j,is), fd = tf1(m,IDN,k,j,is);
        const Real fn[2] = {0.0, (fd > 0.0) ? fd : 0.0};
        for (int q=0; q<2; ++q) {
          if (fd > 0.0) {
            u0(m,n0+q,k,j,is) += bdt*(fn[q] - tf1(m,n0+q,k,j,is))*a/volume(m,k,j,is);
            tf1(m,n0+q,k,j,is) = fn[q];
          }
          v.the_array[q] = -tf1(m,n0+q,k,j,is)*a;
        }
        if (fd > 0.0) v.the_array[4] = fd*a;
      }
      if (tbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
        const Real a = area1(m,k,j,ie+1), fd = tf1(m,IDN,k,j,ie+1);
        const bool win = twin && (fabs(WrapPhi(tx3(m,k) - tphs)) < tww);
        const Real fn[2] = {win ? fd : 0.0, win ? 0.0 : fd};
        for (int q=0; q<2; ++q) {
          if (fd < 0.0) {
            u0(m,n0+q,k,j,ie) -= bdt*(fn[q] - tf1(m,n0+q,k,j,ie+1))*a/volume(m,k,j,ie);
            tf1(m,n0+q,k,j,ie+1) = fn[q];
          }
          v.the_array[2+q] = tf1(m,n0+q,k,j,ie+1)*a;
        }
        if (fd < 0.0 && !win) v.the_array[4] -= fd*a;
      }
      sum += v;
    }, Kokkos::Sum<array_sum::GlobalSum>(t1));
    for (int n=0; n<5; ++n) tfs.the_array[n] = t1.the_array[n];
    if (p3 && thout_) {
      Kokkos::parallel_reduce("ryper_trfx2", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1 + 1)*tnk*tni),
      KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
        const int m = idx/(tnk*tni);
        const int k = (idx - m*tnk*tni)/tni + ks;
        const int i = idx%tni + is;
        array_sum::GlobalSum v;
        for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
        for (int sd=0; sd<2; ++sd) {
          if (sd == 0 && tbcs.d_view(m, BoundaryFace::inner_x2) != BoundaryFlag::user) continue;
          if (sd == 1 && tbcs.d_view(m, BoundaryFace::outer_x2) != BoundaryFlag::user) continue;
          const int jf = (sd == 0) ? js : je + 1, jc = (sd == 0) ? js : je;
          const Real sg = (sd == 0) ? -1.0 : 1.0;     // outward sign of the face flux
          const Real a = ta2(m,k,jf,i), fd = tf2(m,IDN,k,jf,i);
          const bool inflow = (sg*fd < 0.0);
          const Real fn[2] = {0.0, fd};
          for (int q=0; q<2; ++q) {
            if (inflow) {
              u0(m,n0+q,k,jc,i) -= sg*bdt*(fn[q] - tf2(m,n0+q,k,jf,i))*a/volume(m,k,jc,i);
              tf2(m,n0+q,k,jf,i) = fn[q];
            }
            v.the_array[2+q] += sg*tf2(m,n0+q,k,jf,i)*a;
          }
          if (inflow) v.the_array[4] -= sg*fd*a;
        }
        sum += v;
      }, Kokkos::Sum<array_sum::GlobalSum>(t2));
      for (int n=2; n<5; ++n) tfs.the_array[n] += t2.the_array[n];
    }
  }

  // accretor-bigdom-1015: the stream nozzle (stream_inject = volume) and the donor sink
  // (donor_mask = roche), reset after every other source of the stage; sums: 0 nozzle
  // mass added, 1 its inertial Lz, 2 its gas energy (incl. KE, rho Phi), 3 its radiation
  // energy, 4 mask mass removed, 5 mask gas energy removed (per rank)
  array_sum::GlobalSum zs;
  for (int n=0; n<NREDUCTION_VARIABLES; ++n) zs.the_array[n] = 0.0;
  if ((stvol_ && stream_on_) || msk_) {
    const bool nzon = stvol_ && stream_on_, mskon = msk_;
    const Real nlo = nz_rlo_, nhi = nz_rhi_;
    const Real rhos = rho_s_, phs = phi_s_, sg = sig_, ww = nsig_*sig_;
    const Real vrs = vr_s_, vps = vp_s_, om = rp_.omega;
    const Real tcs = gen_tstr_, tsf = gen_tsurf_, sar = strad_ar_, sre = strad_e_;
    const bool stp = stprof_, srad = rad && strad_;
    const Real phl1 = msk_phl1_, xl1 = msk_xl1_, dfl = eosd.dfloor;
    const bool mirr = rad && mirr_;
    const Real med = mirr_e_;
    const Real gdb = 4.0*mirr_beta_, gdc = mirr_ec_, gdf = mirr_ef_;
    const bool gdk = mirr_beta_ > 0.0;
    const bool l1 = stl1_, nzbk = nzback_ && nzon;
    const int ssl = ssl_;
    const Real ssc = mirr_sc_;
    const Real l1c = l1_c_, l1t = l1_t_;
    const RocheParams pr = rp_;
    const bool g3 = g3d_;
    const Real szz = sz_, wz = nsig_*sz_;
    const bool trc = trc_;
    const int n0 = trc_n0_;
    auto phc = ph->phicc0;
    auto &x3c = pmbp->pcoord->x3v;
    auto &x2c = pmbp->pcoord->x2v;
    DvceArray5D<Real> urd;
    if (rad) urd = pmbp->pradm1->u0;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("ryper_nzmsk", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
      const int m = idx/nkji;
      const int k = (idx - m*nkji)/nji + ks;
      const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
      const int i = idx%ni + is;
      const Real r = x1v(m,i), phv = x3c(m,k), dv = volume(m,k,j,i);
      const Real dph = WrapPhi(phv - phs);
      // 3-D: the nozzle's height coordinate z = r_f cos(theta) at its INNER FACE radius r_f
      // (every nozzle cell of a theta row has the same state, as every cell of a phi column,
      // so the upwind face flux equals the normalisation sum exactly); sin(theta) for the mask
      const Real zc = g3 ? nlo*cos(x2c(m,j)) : 0.0;
      const Real sth = g3 ? sin(x2c(m,j)) : 1.0;
      const bool innz = nzon && r > nlo && r < nhi && fabs(dph) < ww && (!g3 || fabs(zc) < wz);
      bool inmk = false, inbk = false;
      if (mskon && !innz) {
        if (g3) {
          inmk = (r*sth*cos(phv) > xl1) && (RochePot3(pr, r, sth, phv) < phl1);
        } else {
          inmk = (r*cos(phv) > xl1) && (RochePot(pr, r, phv) < phl1);
        }
        // stream_inject = l1: the sector BEHIND the nozzle (r above its outer face, |dphi|
        // < 3 sigma + 1 sigma) joins the sink.  The sonic L1 state expands backward into
        // the near-vacuum beyond L1 (v + 3 c_s > 0), where the floor gas got garbage
        // energies and 1e5-1e9 km/s (dt collapse, g8/g9 10-15); the donor's atmosphere
        // that really lies there is not modelled, so it is a sink like the lobe.
        if (!inmk && nzbk && r > nhi && fabs(dph) < ww + sg && (!g3 || fabs(zc) < wz + szz)) {
          inmk = true; inbk = true;
        }
      }
      if (!(innz || inmk)) return;
      array_sum::GlobalSum v;
      for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
      const Real d0 = u0(m,IDN,k,j,i), m30 = u0(m,IM3,k,j,i), e0 = u0(m,IEN,k,j,i);
      Real dn, v1, v3, tn;
      if (innz && l1) {
        // L1 nozzle (Lubow & Shu 1975): isothermal sonic outflow along -x at c_T, Gaussian
        // across the stream; in the rotating frame v = -c_T e_x
        dn = fmax(rhos*exp(-0.5*SQR(dph/sg)), dfl);
        v1 = -l1c*cos(phv); v3 = l1c*sin(phv); tn = l1t;
      } else if (innz && g3) {
        // 3-D: Gaussian in phi (arc sigma) and in z (sigma_z); the diffusion T profile in
        // the elliptical radius of the cross-section
        const Real q2 = SQR(dph/sg) + SQR(zc/szz);
        dn = fmax(rhos*exp(-0.5*q2), dfl);
        v1 = vrs; v3 = vps; tn = StreamT(stp, tcs, tsf, sqrt(q2));
      } else if (innz) {
        dn = fmax(rhos*exp(-0.5*SQR(dph/sg)), dfl);
        v1 = vrs; v3 = vps; tn = StreamT(stp, tcs, tsf, dph/sg);
      } else {
        dn = rfl; v1 = 0.0; v3 = 0.0; tn = p3 ? tcol4(m,k,j,i) : tcol(m,k,i);
      }
      const Real en = eosd.EnergyFromTemperature(dn, tn) + 0.5*dn*(v1*v1 + v3*v3)
                      + dn*phc(m,k,j,i);
      u0(m,IDN,k,j,i) = dn;
      u0(m,IM1,k,j,i) = dn*v1;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = dn*v3;
      u0(m,IEN,k,j,i) = en;
      // acctrc-1019: the nozzle writes stream gas (s0 = 1, s1 = 0); a sink cell is reset to
      // s = 0 (the floor gas at rest it holds is not tracked); the tracer mass each writes
      // or removes goes to registers 14-17 (-> Nz0, Nz1, Sk0, Sk1)
      if (trc) {
        const Real s0o = u0(m,n0,k,j,i), s1o = u0(m,n0+1,k,j,i);
        u0(m,n0,k,j,i) = innz ? dn : 0.0;
        u0(m,n0+1,k,j,i) = 0.0;
        if (innz) {
          v.the_array[14] = (dn - s0o)*dv;
          v.the_array[15] = s1o*dv;
        } else {
          v.the_array[16] = s0o*dv;
          v.the_array[17] = s1o*dv;
        }
      }
      if (innz) {
        v.the_array[0] = (dn - d0)*dv;
        // phys_3d: the inertial Lz with the cylindrical radius r sin(theta)
        const Real rj = p3 ? r*sth : r;
        v.the_array[1] = ((dn*v3 - m30) + om*rj*(dn - d0))*rj*dv;
        v.the_array[2] = (en - e0)*dv;
        if (srad) {
          const Real er = (stp || l1) ? sar*SQR(SQR(tn)) : sre;
          v.the_array[3] = (er - urd(m,radm1::M1_E,k,j,i))*dv;
          urd(m,radm1::M1_E,k,j,i) = er;
          urd(m,radm1::M1_F1,k,j,i) = (4.0/3.0)*er*v1;
          urd(m,radm1::M1_F2,k,j,i) = 0.0;
          urd(m,radm1::M1_F3,k,j,i) = (4.0/3.0)*er*v3;
        }
      } else {
        v.the_array[4] = (d0 - dn)*dv;
        v.the_array[5] = (e0 - en)*dv;
        // acc3d-held-1018: stream_sector_light = neck -- the lobe rule (as surface) only in the
        // sector's physical part, between the nozzle and L1 (r < d_L1) inside the nozzle
        // window; the rest of the 3-D sector (ambient beside and beyond L1) is a dark sink
        const bool neck = (ssl == 4) && (r < xl1) && fabs(dph) < ww && (!g3 || fabs(zc) < wz);
        if (inbk && ssl != 2 && !neck) {
          // v4 (user 10-16): the back-sink sector between the nozzle and L1 holds the
          // STREAM's own light, the nozzle's radiation state continued backward (core T
          // env_t_stream, the same diffusion profile in phi (and z), floor stream_t_surf),
          // at rest (F = 0), within the nozzle's width (3 sigma); outside it, its radiation
          // is left to the transport (dark, as the ambient).  Not scaled (not donor light);
          // booked with the nozzle's radiation (register 3).
          // stream_sector_light (10-17): full = this; scaled = x donor_irr_scale (the band's
          // 2-D -> 3-D factor of the held donor light); dark = no light held (all of the
          // sector left to the transport); surface = the lobe rule below (donor light)
          const bool inw = fabs(dph) < ww && (!g3 || fabs(zc) < wz);
          if (srad && inw && ssl != 3 && ssl != 4) {
            const Real tb = g3 ? StreamT(stp, tcs, tsf, sqrt(SQR(dph/sg) + SQR(zc/szz)))
                               : StreamT(stp, tcs, tsf, dph/sg);
            Real er = stp ? sar*SQR(SQR(tb)) : sre;
            if (ssl == 1) er *= ssc;
            v.the_array[3] = (er - urd(m,radm1::M1_E,k,j,i))*dv;
            urd(m,radm1::M1_E,k,j,i) = er;
            urd(m,radm1::M1_F1,k,j,i) = 0.0;
            urd(m,radm1::M1_F2,k,j,i) = 0.0;
            urd(m,radm1::M1_F3,k,j,i) = 0.0;
          }
        } else if (mirr) {
          // donor_irr: the lobe is a black body (isotropic, comoving F = 0); its light
          // leaves through the mask surface by the M1 transport.  donor_gdark > 0: the
          // local T of the von Zeipel surface along the ray from the donor centre
          Real eh = med;
          if (gdk) {
            const Real cph = cos(phv), sph = sin(phv);
            const Real xx = r*sth*cph, yy = r*sth*sph, zz = g3 ? r*cos(x2c(m,j)) : 0.0;
            const Real ge = DonorSurfGeff(pr, phl1, xl1, xx, yy, zz);
            eh = fmax(gdc*pow(ge, gdb), gdf);
          }
          v.the_array[6] = (eh - urd(m,radm1::M1_E,k,j,i))*dv;
          urd(m,radm1::M1_E,k,j,i) = eh;
          urd(m,radm1::M1_F1,k,j,i) = 0.0;
          urd(m,radm1::M1_F2,k,j,i) = 0.0;
          urd(m,radm1::M1_F3,k,j,i) = 0.0;
        }
      }
      sum += v;
    }, Kokkos::Sum<array_sum::GlobalSum>(zs));
  }

  if (nstages_ <= 0) return;
  auto &flx = ph->uflx.x1f;
  auto &xf = pmbp->pcoord->xx1f;
  auto &x3v = pmbp->pcoord->x3v;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto ira = iracc_;
  const Real om = rp_.omega, phis = phi_s_, wwin = nsig_*sig_;
  const bool strm = stream_on_ && !stvol_;   // the r_out window (boundary injection only)
  const int nzif = (stvol_ && stream_on_) ? nz_if_ : -1;
  const bool g3f = g3d_;
  const Real wzf = nsig_*sz_;
  auto &x2vf = pmbp->pcoord->x2v;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj;
  DvceArray4D<Real> ff0;
  if (rad) ff0 = pmbp->pradm1->f0x1;
  const bool haveff = rad && (ff0.extent_int(0) > 0);
  const Real fin = rad_fin_;
  // per-column outer bath (env_amb_mode = radeq or stream_rad): the advective radiation
  // energy through r_out (implicit_bc_advect) joins the comoving face flux in register 10
  const bool hbdv = rad && pmbp->pradm1->iebath_hi_col_on;
  DvceArray3D<Real> bdv;
  if (hbdv) bdv = pmbp->pradm1->ibadv_x1max;
  // phys_3d: Lz with the cylindrical radius; theta_bc = outflow: the theta faces
  const bool thf = p3 && thout_;
  auto &flx2 = ph->uflx.x2f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &x2ff = pmbp->pcoord->xx2f;
  array_sum::GlobalSum fs;
  Kokkos::parallel_reduce("ryper_flxe", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
  KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
    const int m = idx/(nk*nj);
    const int k = (idx - m*nk*nj)/nj + ks;
    const int j = idx%nj + js;
    array_sum::GlobalSum v;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) v.the_array[n] = 0.0;
    const Real sj = p3 ? sin(x2vf(m,j)) : 1.0;
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      const Real a = area1(m,k,j,is), r = p3 ? xf(m,is)*sj : xf(m,is);
      const Real fd = flx(m,IDN,k,j,is), f3 = flx(m,IM3,k,j,is);
      v.the_array[2] = -fd*a;
      v.the_array[3] = -(f3 + om*r*fd)*r*a;
      v.the_array[4] = -(f3 - fd*w0(m,IVZ,k,j,is))*r*a;
      if (rad) v.the_array[9] = fin*a;
    }
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      const Real a = area1(m,k,j,ie+1), r = p3 ? xf(m,ie+1)*sj : xf(m,ie+1);
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
        if (hbdv) v.the_array[10] += bdv(m,k,j)*a;
        v.the_array[11] = flx(m,IEN,k,j,ie+1)*a;
      }
    }
    // theta_bc = outflow: mass, inertial Lz and hydro energy (incl. rho Phi) out through the
    // theta walls join Mout (1), Jout (5) and the hydro energy out (11)
    if (thf) {
      for (int sd=0; sd<2; ++sd) {
        if (sd == 0 && !(j == js && mbbcs.d_view(m, BoundaryFace::inner_x2) ==
                         BoundaryFlag::user)) continue;
        if (sd == 1 && !(j == je && mbbcs.d_view(m, BoundaryFace::outer_x2) ==
                         BoundaryFlag::user)) continue;
        const int jf = (sd == 0) ? js : je + 1;
        const Real sg = (sd == 0) ? -1.0 : 1.0;
        const Real sf = sin(x2ff(m,jf));
        for (int i=is; i<=ie; ++i) {
          const Real a = area2(m,k,jf,i), rr = x1v(m,i)*sf;
          const Real fd = flx2(m,IDN,k,jf,i), f3 = flx2(m,IM3,k,jf,i);
          v.the_array[1] += sg*fd*a;
          v.the_array[5] += sg*(f3 + om*rr*fd)*rr*a;
          if (rad) v.the_array[11] += sg*flx2(m,IEN,k,jf,i)*a;
        }
      }
    }
    const int ir = ira(m);
    if (ir >= 0) {
      const Real a = area1(m,k,j,ir), r = p3 ? xf(m,ir)*sj : xf(m,ir);
      const Real fd = flx(m,IDN,k,j,ir), f3 = flx(m,IM3,k,j,ir);
      v.the_array[6] = -fd*a;
      v.the_array[7] = -(f3 + om*r*fd)*r*a;
    }
    // the nozzle's inner face (Mnzf), inside its phi window
    if (nzif >= 0 && fabs(WrapPhi(x3v(m,k) - phis)) < wwin &&
        (!g3f || fabs(xf(m,nzif)*cos(x2vf(m,j))) < wzf)) {
      v.the_array[13] = -flx(m,IDN,k,j,nzif)*area1(m,k,j,nzif);
    }
    sum += v;
  }, Kokkos::Sum<array_sum::GlobalSum>(fs));

  const int s = stage_ctr_ % nstages_;
  const int nacc = rad ? kNaccR : kNacc;
  // env_amb_mode = radeq: register 4 (the wall-stress part Jstr of Jacc, meaningless deep in
  // the envelope) holds the cumulative floor mass Mfl instead
  if (amb_radeq_) fs.the_array[4] = 0.0;
  // acc3d-1017: the imposed x2-wall radiative flux (implicit_bc_x2 = flux) leaves as radiation
  if (rad) fs.the_array[10] += fx2leak_;
  // donor_mask: registers 2, 3 (Mwal, Jwal: zero with the closed free-slip wall) hold the
  // cumulative mass and gas energy removed by the mask
  if (msk_) fs.the_array[2] = fs.the_array[3] = 0.0;
  if (s == 0) {
    for (int n=0; n<nacc; ++n) acc1_[n] = acc0_[n];
  }
  for (int n=0; n<kNacc; ++n) {
    acc0_[n] = rk_g0_[s]*acc0_[n] + rk_g1_[s]*acc1_[n] + bdt*fs.the_array[n];
  }
  if (amb_radeq_) acc0_[4] += dmfl;
  // the nozzle's additions (Min, Jin) and the mask's removals (Mmsk, Emsk), as increments
  // (they act on the stage's u like the floor, after its RK combination)
  acc0_[0] += zs.the_array[0];
  acc0_[8] += zs.the_array[1];
  if (msk_) {
    acc0_[2] += zs.the_array[4];
    acc0_[3] += zs.the_array[5];
  }
  if (accx_on_) {
    if (s == 0) {
      accx1_[0] = accx0_[0];
      accx1_[1] = accx0_[1];
    }
    accx0_[0] = rk_g0_[s]*accx0_[0] + rk_g1_[s]*accx1_[0] + bdt*fs.the_array[13];
    accx0_[1] += zs.the_array[6];   // accwall-1019: plain sum (see erinj_)
  }
  if (rad) {
    for (int n=kNacc; n<kNaccR - 1; ++n) {
      acc0_[n] = rk_g0_[s]*acc0_[n] + rk_g1_[s]*acc1_[n] + bdt*fs.the_array[n];
    }
    acc0_[kNaccR-1] = rk_g0_[s]*acc0_[kNaccR-1] + rk_g1_[s]*acc1_[kNaccR-1] + derel;
    // the nozzle's energy (and the donor_irr radiation) enters Ebnd = in - rad out (10)
    // - hydro out (11)
    erinj_ += zs.the_array[3] + zs.the_array[6];
    acc0_[11] -= zs.the_array[2];
  }
  // acctrc-1019: the tracer registers (fluxes RK weighted as acc0_, resets as increments)
  if (trc_) {
    if (s == 0) {
      for (int n=0; n<kNtr; ++n) trc1_[n] = trc0_[n];
    }
    for (int n=0; n<kNtr; ++n) {
      trc0_[n] = rk_g0_[s]*trc0_[n] + rk_g1_[s]*trc1_[n] + ((n >= 4) ? bdt*tfs.the_array[n-4] : 0.0);
    }
    for (int n=0; n<4; ++n) trc0_[n] += zs.the_array[14+n];
  }
  stage_ctr_++;
}

//----------------------------------------------------------------------------------------
//! stage 2 history: Mdom Jdom, Menv Jenv (cells with r < R_acc), then the nine cumulative
//! boundary integrals and their rates over the last history interval
//----------------------------------------------------------------------------------------
//! acctrc-1019: the tracer diagnostics, one row per history output in <basename>.trc.hst
//! (appended; ryper.user.hst is not touched).  Columns: time, Ms0, Ms1 (domain integrals of
//! rho s0, rho s1), the nine cumulative registers Nz0 Nz1 Sk0 Sk1 Fi0 Fi1 Fo0 Fo1 Bin (see
//! kNtr), min/max of s0 and s1 and max(s0 + s1) over the active cells.  Budgets (to the
//! round-off of the c2p clip of negative rho s): Ms0 - Ms0(start) = Nz0 - Sk0 - Fi0 - Fo0,
//! Ms1 - Ms1(start) = - Nz1 - Sk1 - Fi1 - Fo1; the registers continue across a restart
//! (acctrc-1019c rst block; from 0 only after a restart from a pre-1019c file).
void RyPerTrcHist(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  auto u0 = pmbp->phydro->u0;
  auto &vol = pmbp->pcoord->volume;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const int n0 = trc_n0_;
  Real ms[2] = {0.0, 0.0};
  Real mn0 = 0.0, mx0 = 0.0, mn1 = 0.0, mx1 = 0.0, mxs = 0.0;
  Kokkos::parallel_reduce("ryper_trh", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a0, Real &a1, Real &l0, Real &h0, Real &l1, Real &h1,
                Real &hs) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    const int i = idx%nx1 + is;
    const Real dv = vol(m,k,j,i), d = u0(m,IDN,k,j,i);
    const Real q0 = u0(m,n0,k,j,i), q1 = u0(m,n0+1,k,j,i);
    a0 += q0*dv;
    a1 += q1*dv;
    const Real s0 = q0/d, s1 = q1/d;
    l0 = fmin(l0, s0); h0 = fmax(h0, s0);
    l1 = fmin(l1, s1); h1 = fmax(h1, s1);
    hs = fmax(hs, s0 + s1);
  }, Kokkos::Sum<Real>(ms[0]), Kokkos::Sum<Real>(ms[1]), Kokkos::Min<Real>(mn0),
     Kokkos::Max<Real>(mx0), Kokkos::Min<Real>(mn1), Kokkos::Max<Real>(mx1),
     Kokkos::Max<Real>(mxs));
  Real sb[2 + kNtr];
  sb[0] = ms[0]; sb[1] = ms[1];
  for (int n=0; n<kNtr; ++n) sb[2+n] = trc0_[n];
  Real lo[2] = {mn0, mn1}, hi[3] = {mx0, mx1, mxs};
#if MPI_PARALLEL_ENABLED
  if (global_variable::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, sb, 2 + kNtr, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, lo, 2, MPI_ATHENA_REAL, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, hi, 3, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(sb, sb, 2 + kNtr, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(lo, lo, 2, MPI_ATHENA_REAL, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(hi, hi, 3, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
  }
#endif
  if (global_variable::my_rank != 0) return;
  bool hdr = false;
  {
    std::ifstream f(trc_file_);
    hdr = !f.good() || f.peek() == std::ifstream::traits_type::eof();
  }
  FILE *fp = std::fopen(trc_file_.c_str(), "a");
  if (fp == nullptr) return;
  if (hdr) {
    std::fprintf(fp, "# ry_per_accretor tracers (acctrc-1019): s0 stream, s1 boundary inflow; "
                 "cumulative columns continue across restarts (rst block, acctrc-1019c)\n#  [1]=time  [2]=Ms0  "
                 "[3]=Ms1  [4]=Nz0  [5]=Nz1  [6]=Sk0  [7]=Sk1  [8]=Fi0  [9]=Fi1  [10]=Fo0  "
                 "[11]=Fo1  [12]=Bin  [13]=s0min  [14]=s0max  [15]=s1min  [16]=s1max  "
                 "[17]=smax\n");
  }
  std::fprintf(fp, " %.15e", pm->time);
  for (int n=0; n<2 + kNtr; ++n) std::fprintf(fp, " %.15e", sb[n]);
  std::fprintf(fp, " %.15e %.15e %.15e %.15e %.15e\n", lo[0], hi[0], lo[1], hi[1], hi[2]);
  std::fclose(fp);
}

void RyPerHistEnv(HistoryData *pdata, Mesh *pm) {
  if (trc_) RyPerTrcHist(pm);
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
  const bool p3 = p3_;
  auto &x2v = pmbp->pcoord->x2v;
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
    // phys_3d: the cylindrical radius about the accretor's z axis
    const Real rz = p3 ? r*sin(x2v(m,j)) : r;
    const Real jz = dv*rz*(u0(m,IM3,k,j,i) + om*rz*d);
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

  // problem/dt_report (bigdom-1015): which cell sets the hydro dt, at every hst output.
  // accretor-3d-1017: the GLOBAL minimum over the ranks (MPI MINLOC of each rank's own
  // minimum dtnew; was rank 0's local minimum); the owning rank broadcasts the cell, rank 0
  // prints it (the old line, then the rank and theta of the cell)
  if (dtrep_) {
    auto *phd = pmbp->phydro;
    const int dm = phd->dtnew_m, dk = phd->dtnew_k, dj = phd->dtnew_j, di = phd->dtnew_i;
    const bool have = (dm >= 0);
    double hb[12];
    for (int n=0; n<12; ++n) hb[n] = 0.0;
    if (have) {
      DvceArray1D<Real> dv("ryper_dtrep", 8);
      auto w0d = phd->w0;
      auto &x3d = pmbp->pcoord->x3v;
      auto &x2d = pmbp->pcoord->x2v;
      par_for("ryper_dtrep_f", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int n) {
        dv(0) = x1v(dm,di);
        dv(1) = x3d(dm,dk);
        dv(2) = w0d(dm,IDN,dk,dj,di);
        dv(3) = w0d(dm,IVX,dk,dj,di);
        dv(4) = w0d(dm,IVY,dk,dj,di);
        dv(5) = w0d(dm,IVZ,dk,dj,di);
        dv(6) = w0d(dm,IEN,dk,dj,di);
        dv(7) = x2d(dm,dj);
      });
      auto hv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dv);
      for (int n=0; n<8; ++n) hb[n] = hv(n);
      hb[8] = dm; hb[9] = dk; hb[10] = dj; hb[11] = di;
    }
    struct {double v; int r;} loc{have ? static_cast<double>(phd->dtnew) : 1.0e300,
                                  global_variable::my_rank};
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &loc, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    MPI_Bcast(hb, 12, MPI_DOUBLE, loc.r, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0 && loc.v < 1.0e299) {
      std::printf("ry_per_accretor dt_report: t %.6e dt %.4e cell (m,k,j,i) (%d,%d,%d,%d) r %.5f "
                  "phi %.4f rad: rho %.4g v_r %.4g v_th %.4g v_phi %.4g km/s e/rho %.4g "
                  "(km/s)^2 [rank %d, theta %.4f deg]\n", pm->time, pm->dt,
                  static_cast<int>(hb[8]), static_cast<int>(hb[9]), static_cast<int>(hb[10]),
                  static_cast<int>(hb[11]), hb[0], hb[1], hb[2], hb[3], hb[4], hb[5],
                  hb[6]/hb[2], loc.r, hb[7]*180.0/M_PI);
    }
  }
  const char *lab[22] = {"Mdom", "Jdom", "Menv", "Jenv",
                         "Min", "Mout", "Mwal", "Jwal", "Jstr", "Jout", "MR", "JR", "Jin",
                         "dMin", "dMout", "dMwal", "dJwal", "dJstr", "dJout", "dMR",
                         "dJR", "dJin"};
  pdata->nhist = 22;
  for (int n=0; n<22; ++n) pdata->label[n] = lab[n];
  if (amb_radeq_) pdata->label[8] = "Mfl";   // cumulative floor mass (env_amb_mode = radeq)
  if (msk_) {                                  // donor_mask = roche (bigdom-1015)
    pdata->label[6] = "Mmsk";
    pdata->label[7] = "Emsk";
  }
  if (accx_on_) {                              // columns 23/24 (bigdom-1015)
    pdata->nhist = 24;
    pdata->label[22] = "Mnzf";
    pdata->label[23] = "Eirr";
    pdata->hdata[22] = accx0_[0];
    pdata->hdata[23] = accx0_[1];
  }
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
    // (all energies code units of the band; cumulative columns and Picmax continue across
    // restarts, acctrc-1019c)
    const char *lr[9] = {"Lmeas", "Ltop", "Erad", "Etot", "KE", "Ebnd", "Erel", "Pic",
                         "Picmax"};
    for (int n=0; n<9; ++n) pdata->label[13+n] = lr[n];
    for (int n=0; n<kNacc; ++n) pdata->hdata[4+n] = acc0_[n];
    pdata->hdata[13] = s.the_array[7]*rad_lfac_;
    pdata->hdata[14] = s.the_array[8]*rad_lfac_;
    pdata->hdata[15] = s.the_array[4];
    pdata->hdata[16] = s.the_array[5];
    pdata->hdata[17] = s.the_array[6];
    pdata->hdata[18] = acc0_[9] - acc0_[10] - acc0_[11] + erinj_;
    pdata->hdata[19] = acc0_[12];
    pdata->hdata[20] = 0.0;
    pdata->hdata[21] = 0.0;
    auto *pm1 = pmbp->pradm1;
    if (global_variable::my_rank == 0) {
      const Real dn = pm1->impl_nstep - rad_pn0_, ds = pm1->impl_itsum - rad_ps0_;
      pdata->hdata[20] = (dn > 0.0) ? ds/dn : 0.0;
      pdata->hdata[21] = fmax(picmax_rst_, pm1->impl_itmax);   // global (1019c)
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
  rbath_ = DvceArray3D<Real>();
  fx2lo_ = DvceArray3D<Real>();
  fx2hi_ = DvceArray3D<Real>();
  gr4_ = DvceArray4D<Real>();
  gt4_ = DvceArray4D<Real>();
  gp4_ = DvceArray4D<Real>();
  gpl4_ = DvceArray4D<Real>();
  cref24_ = DvceArray4D<Real>();
  rhoamb4_ = DvceArray4D<Real>();
  tcol4_ = DvceArray4D<Real>();
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
  // ---- the stream window under radiation (S4): the injected gas carries its T through the
  // EOS (EnergyFromTemperature, RyPerBCEnv), its own radiation field is NOT injected (dark
  // r_out ghosts, Marshak window; user 10-08), and its temperature must not be relaxed
  // (env_t_relax_stream = 0).  Report the peak state so the no-injection premise
  // (P_rad << P_gas) and the stream's optical depth can be checked from the log.
  if (stream_on_) {
    if (env_trs_ != 0.0) {
      fatal("stream on under <rad_m1> needs problem/env_t_relax_stream = 0 (the stream's T "
            "is set by the EOS state at the window and by radiation, not relaxed)");
    }
    const Real d = rho_s_, t = stl1_ ? l1_t_ : gen_tstr_;
    const Real e = eos.EnergyFromTemperature(d, t);
    const Real pg = eos.Pressure(d, e, t);
    const Real prad = ar*SQR(SQR(t))/3.0;
    DvceArray1D<Real> sk("ryper_sk", 2);
    radm1::M1OpacTab ot = pm1->otab;
    par_for("ryper_skt", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int n) {
      Real op, oe, of, os;
      radm1::M1TableOpacities(ot, d, t, op, oe, of, os);
      sk(0) = of + os;
      sk(1) = op;
    });
    auto hsk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sk);
    const Real wa = sig_*str_r_;      // arc sigma at the injection radius (code)
    if (root) {
      std::printf("ry_per_accretor: <rad_m1> stream window peak: rho %.4g g/cc, T %.1f K, "
                  "P_rad/P_gas %.4g (%s), kappa_R %.4g "
                  "kappa_P %.4g cm^2/g, tau across one arc sigma (%.4g Rsun) %.4g\n",
                  d*du, t*tk, prad/pg, strad_ ? "radiation INJECTED: E = a T^4 = "
                  "window ghosts and outer Marshak bath, F = (4/3) E v" :
                  "radiation NOT injected: dark ghosts", hsk(0)/(du*lu),
                  hsk(1)/(du*lu), wa, hsk(0)*d*wa);
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
  // problem/env_lin_band_factor (default 1): L_in = band equivalent of the theta-thin
  // equatorial slab.  With T = T(psi) and F ~ g (von Zeipel) the band of the Roche-
  // distorted star radiates <(g/g90)(r/r90)^2/cos alpha>_phi of the phi = 90 column's L at
  // the photosphere (docs/dev/accretor_rhd/band_luminosity.py: 0.960504 for the Plaskett
  // gainer, geometric only), and the slab cannot shed flux to the poles; the imposed
  // bottom flux must then be the factor times F_col(r_in) so that the band radiates what
  // it receives at the column's T_eff.  A 3-D run (full theta) must use 1.
  rad_fin_ = pin->GetReal("rad_m1", "implicit_flux_x1min");
  const Real lbf = pin->GetOrAddReal("problem", "env_lin_band_factor", 1.0);
  if (!(lbf > 0.0 && lbf <= 1.0)) fatal("problem/env_lin_band_factor must be in (0, 1]");
  const Real fcol = lbf*RadColAt(rcol_f_, psi_in);
  if (!(std::fabs(rad_fin_/fcol - 1.0) <= 1.0e-3)) {
    std::ostringstream os;
    os << "<rad_m1>/implicit_flux_x1min = " << rad_fin_ << " but env_lin_band_factor ("
       << lbf << ") x the column's F(r_in) = " << fcol << " (code units, " << fcol*fu
       << " erg/cm^2/s)";
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
  // accretor-3d-1017 (phys_3d): Phi_eff = Phi_wb + G(psi) on the 3-D equipotentials, column
  // by (m,k,j); the plain-gravity force above the WB region into gpl4_
  if (p3_) {
    auto hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
    auto hgl4 = Kokkos::create_mirror_view(gpl4_);
    std::vector<Real> pf(n1 + 1), gf(n1 + 1);
    for (int m=0; m<nmb; ++m) {
      for (int k=0; k<n3; ++k) {
        const Real ph = hx3(m,k);
        for (int j=0; j<n2; ++j) {
          const Real st = std::sin(hx2(m,j));
          for (int i=0; i<=n1; ++i) {
            const Real w = PhiWB3(p, dlt, rtop, rspin, phtop, hx1f(m,i), st, ph);
            gf[i] = RadColAt(rcol_g_, phis - w);
            pf[i] = w + gf[i];
          }
          for (int i=0; i<n1; ++i) {
            const Real w = PhiWB3(p, dlt, rtop, rspin, phtop, hx1(m,i), st, ph);
            const Real dr = hx1f(m,i+1) - hx1f(m,i);
            const Real a = -(gf[i+1] - gf[i])/dr;
            hgl4(m,k,j,i) = -(pf[i+1] - pf[i])/dr;
            amax = std::max(amax, a);
            hpc(m,k,j,i) = w + RadColAt(rcol_g_, phis - w);
            hpf(m,k,j,i) = pf[i];
            if (i == n1 - 1) hpf(m,k,j,i+1) = pf[i+1];
            ha(m,k,j,i) = a;
          }
        }
      }
    }
    Kokkos::deep_copy(gpl4_, hgl4);
    // acc3d-1017: <rad_m1>/implicit_bc_x2min|max = flux: the wall faces carry the column's
    // EQUILIBRIUM lateral flux, F_theta = F_col(psi) (dPhi/(r dtheta))/|grad Phi| at the face
    // (F along grad Phi of the true 3-D Roche potential, i.e. normal to the tilted
    // equipotentials; F_col at the face's psi, no von Zeipel g factor), so the wall neither
    // blocks the poleward flux of the tilted photosphere (reflect) nor drains the star
    // (outflow/Marshak).  Explicit (imposed); its outward total is booked in Ebnd.
    fx2leak_ = 0.0;
    if (pm1->ibc_x2min == radm1::M1_IBC_X2FLUX || pm1->ibc_x2max == radm1::M1_IBC_X2FLUX) {
      Kokkos::realloc(fx2lo_, nmb, n3, n1);
      Kokkos::realloc(fx2hi_, nmb, n3, n1);
      auto hlo = Kokkos::create_mirror_view(fx2lo_);
      auto hhi = Kokkos::create_mirror_view(fx2hi_);
      auto hx2f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx2f);
      auto ha2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x2f);
      auto hbc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pmb->mb_bcs.d_view);
      const int js = indcs.js, je = indcs.je, is = indcs.is, ie = indcs.ie, ks = indcs.ks;
      const int ke = indcs.ke;
      // |grad Phi| of the true 3-D Roche potential at the column's r_eq (phi = 90, plane)
      std::vector<Real> gcol(cpsi.size());
      for (size_t n=0; n<cpsi.size(); ++n) {
        const Real e = 1.0e-6, rq = crq[n], h = 0.5*M_PI;
        const Real a = (RochePot3(p, rq*(1.0 + e), 1.0, h) - RochePot3(p, rq*(1.0 - e), 1.0, h))
                       /(2.0*e*rq);
        const Real b = (RochePot3(p, rq, 1.0, h + e) - RochePot3(p, rq, 1.0, h - e))/(2.0*e*rq);
        gcol[n] = std::sqrt(a*a + b*b);
      }
      auto fth = [&](const Real r, const Real th, const Real ph) {
        const Real e = 1.0e-6;
        const Real gt = (RochePot3(p, r, std::sin(th + e), ph) -
                         RochePot3(p, r, std::sin(th - e), ph))/(2.0*e*r);
        const Real gr = (RochePot3(p, r*(1.0 + e), std::sin(th), ph) -
                         RochePot3(p, r*(1.0 - e), std::sin(th), ph))/(2.0*e*r);
        const Real gp = (RochePot3(p, r, std::sin(th), ph + e) -
                         RochePot3(p, r, std::sin(th), ph - e))/(2.0*e*r*std::sin(th));
        const Real g = std::sqrt(gt*gt + gr*gr + gp*gp);
        const Real w = PhiWB3(p, dlt, rtop, rspin, phtop, r, std::sin(th), ph);
        // only inside the photosphere (psi > 0, the diffusion regime the von Zeipel flux
        // describes); above it the wall face carries no imposed flux (reflecting): the
        // column's F capped at the top over the whole outer wall drained the star (the first
        // try imposed 1.4x L_r out sideways; eint <= 0 retries at the walls)
        if (!(phis - w > 0.0)) return 0.0;
        // von Zeipel (10-18): on an equipotential F ~ g_eff, so the column's F (taken at
        // its own position r_eq(psi), theta = pi/2, phi = pi/2) scales by g/g_col(psi); the
        // lateral component is then F_col (g/g_col) gt/g = F_col gt/g_col
        return RadColAt(rcol_f_, phis - w)*gt/std::fmax(RadColAt(gcol, phis - w), 1.0e-30);
      };
      Real lk = 0.0;
      for (int m=0; m<nmb; ++m) {
        const Real tlo = hx2f(m,js), thi = hx2f(m,je+1);
        const bool plo = (hbc(m, BoundaryFace::inner_x2) != BoundaryFlag::block &&
                          hbc(m, BoundaryFace::inner_x2) != BoundaryFlag::periodic);
        const bool phi2 = (hbc(m, BoundaryFace::outer_x2) != BoundaryFlag::block &&
                           hbc(m, BoundaryFace::outer_x2) != BoundaryFlag::periodic);
        for (int k=0; k<n3; ++k) {
          for (int i=0; i<n1; ++i) {
            hlo(m,k,i) = fth(hx1(m,i), tlo, hx3(m,k));
            hhi(m,k,i) = fth(hx1(m,i), thi, hx3(m,k));
            if (k >= ks && k <= ke && i >= is && i <= ie) {
              if (plo && pm1->ibc_x2min == radm1::M1_IBC_X2FLUX) lk -= hlo(m,k,i)*ha2(m,k,js,i);
              if (phi2 && pm1->ibc_x2max == radm1::M1_IBC_X2FLUX) {
                lk += hhi(m,k,i)*ha2(m,k,je+1,i);
              }
            }
          }
        }
      }
      fx2leak_ = lk;
      Kokkos::deep_copy(fx2lo_, hlo);
      Kokkos::deep_copy(fx2hi_, hhi);
      pm1->SetX2FaceFlux(fx2lo_, fx2hi_);
      if (root) {
        std::printf("ry_per_accretor: <rad_m1> x2 walls: IMPOSED equilibrium lateral flux "
                    "F_col(psi) (g/g_col) dPhi/(r dth)/|grad Phi| (von Zeipel); outward total (rank 0) %.6g code "
                    "(= %.4g Lsun of the full sphere scale)\n", lk, lk*rad_lfac_);
      }
    }
  }
  for (int m=0; m<nmb && !p3_; ++m) {
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
                "bottom flux %.6g code = %.6g erg/cm^2/s (L %.5g Lsun at r_in, phi 90; "
                "env_lin_band_factor %.6g)\n",
                cpsi.front(), cpsi.back(), g_in, g_in/psi_in, amax, rad_fin_,
                rad_fin_*fu, 4.0*M_PI*SQR(r_in*lu)*rad_fin_*fu/3.828e33,
                pin->GetReal("problem", "env_lin_band_factor"));
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
  // accretor-3d-1017 (phys_3d): F_col(psi) and T_col(psi) on the 3-D equipotentials
  const bool p3 = p3_;
  DvceArray4D<Real> fc4, ffc4;
  if (p3) {
    fc4 = DvceArray4D<Real>("ryper_rfc4", nmb, n3, n2, n1);
    ffc4 = DvceArray4D<Real>("ryper_rff4", nmb, n3, n2, n1 + 1);
    auto hfc4 = Kokkos::create_mirror_view(fc4);
    auto hff4 = Kokkos::create_mirror_view(ffc4);
    auto hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
    for (int m=0; m<nmb; ++m) {
      for (int k=0; k<n3; ++k) {
        for (int j=0; j<n2; ++j) {
          const Real st = std::sin(hx2(m,j));
          for (int i=0; i<n1; ++i) {
            hfc4(m,k,j,i) = RadColAt(rcol_f_, phis - PhiWB3(p, dlt, rtop, rspin, phtop,
                                                             hx1(m,i), st, hx3(m,k)));
          }
          for (int i=0; i<=n1; ++i) {
            hff4(m,k,j,i) = RadColAt(rcol_f_, phis - PhiWB3(p, dlt, rtop, rspin, phtop,
                                                             hx1f(m,i), st, hx3(m,k)));
          }
        }
      }
    }
    Kokkos::deep_copy(fc4, hfc4);
    Kokkos::deep_copy(ffc4, hff4);
  }
  for (int m=0; m<nmb && !p3; ++m) {
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
  auto tc4 = tcol4_;
  const Real ar = pm1->arad, efl = pm1->e_floor, cl = pm1->c_light;
  par_for("ryper_ric", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = p3 ? tc4(m,k,j,i) : tc(m,k,i);
    const Real e = fmax(ar*SQR(SQR(t)), efl);
    ur(m,radm1::M1_E,k,j,i) = e;
    ur(m,radm1::M1_F1,k,j,i) = fmin(p3 ? fc4(m,k,j,i) : fc(m,k,i), 0.9*cl*e);
    ur(m,radm1::M1_F2,k,j,i) = 0.0;
    ur(m,radm1::M1_F3,k,j,i) = 0.0;
  });
  auto ff0 = pm1->f0x1;
  if (ff0.extent_int(0) >= nmb) {
    const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
    const int e1 = std::min(ff0.extent_int(3) - 1, n1);
    par_for("ryper_ricf", DevExeSpace(), 0, nmb - 1, 0, e3, 0, e2, 0, e1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      ff0(m,k,j,i) = p3 ? ffc4(m,k,j,i) : ffc(m,k,i);
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
  // accretor-3d-1017 (phys_3d): the same march for every (m,k,j) column on the 3-D
  // equipotentials (Phi_wb, Phi_eff, T_col and the ambient all depend on theta)
  const bool p3 = p3_;
  const int nj3 = p3 ? n2 : 1;
  DvceArray4D<Real> dd4, de4, dv4;
  decltype(Kokkos::create_mirror_view(dd4)) hd4, he4, hv4;
  decltype(Kokkos::create_mirror_view_and_copy(HostMemSpace(), tcol4_)) ht4;
  decltype(Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v)) hx2;
  if (p3) {
    dd4 = DvceArray4D<Real>("ryper_icd4", nmb, n3, n2, n1);
    de4 = DvceArray4D<Real>("ryper_ice4", nmb, n3, n2, n1);
    dv4 = DvceArray4D<Real>("ryper_icv4", nmb, n3, n2, n1);
    hd4 = Kokkos::create_mirror_view(dd4);
    he4 = Kokkos::create_mirror_view(de4);
    hv4 = Kokkos::create_mirror_view(dv4);
    ht4 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), tcol4_);
    hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
  }
  std::vector<Real> d(n1), e(n1), t(n1), pw(n1), pf(n1 + 1);
  std::vector<char> isa(n1, 0);
  const bool radeq = amb_radeq_;
  Real rmax = 0.0, rhoanc_min = 1e300, rhoanc_max = -1e300;
  int nfail = 0;
  for (int m=0; m<nmb; ++m) {
    for (int k=0; k<n3; ++k) {
    for (int jc=0; jc<nj3; ++jc) {
      const Real st = p3 ? std::sin(hx2(m,jc)) : 1.0;
      for (int i=0; i<n1; ++i) {
        t[i] = p3 ? ht4(m,k,jc,i) : ht(m,k,i);
        pw[i] = hpc(m,k,jc,i);
        pf[i] = hpf(m,k,jc,i);
      }
      pf[n1] = hpf(m,k,jc,n1);
      // the Phi_wb of cell i (psi = Phi_s - Phi_wb; under <rad_m1> pw holds Phi_eff)
      auto pwbc = [&](const int i) {
        return p3 ? PhiWB3(p, env_dlt_, rtop, rspin, phtop, hx1(m,i), st, hx3(m,k))
                  : PhiWB(p, env_dlt_, rtop, rspin, phtop, hx1(m,i), hx3(m,k));
      };
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
          const Real psi = rad_ ? (phis - pwbc(i)) : (phis - pw[i]);
          if (std::fabs(psi) < best) {
            best = std::fabs(psi);
            ia = i;
          }
        }
        Real lt, lr;
        const Real psia = rad_ ? (phis - pwbc(ia)) : (phis - pw[ia]);
        ColumnAt(cpsi, clt, clr, psia, lt, lr);
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
      // upward to n1 - 2; the hot ambient takes over (env_amb_mode = radeq: the floor at the
      // column's T, where the column's density has fallen to rho_amb)
      bool amb = false;
      for (int i=0; i<n1; ++i) isa[i] = 0;
      for (int i=ia; i<=n1-1; ++i) {
        const Real pha = hpt(m,k,jc,i) - (eqtop ? aeq : (p3 ? RochePot3(p, racc, st, hx3(m,k))
                                                            : RochePot(p, racc, hx3(m,k))));
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
        if (i > ia && (amb || hx1(m,i) > rtop || pw[i] >= phtop_c ||
                       (radeq ? (d[i] < rho_amb) : (pc < da*camb2)) || i == n1 - 1)) {
          amb = true;
          isa[i] = 1;
          if (radeq) {
            d[i] = rho_amb;          // t[i] stays T_col: E = a T_col^4 in RadIC
          } else {
            d[i] = std::max(da, rho_amb);
            t[i] = tamb;
          }
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
        const Real vi = (hx1(m,i) < rspin && (radeq ? (isa[i] == 0) : (t[i] != tamb))) ?
                        vsp*hx1(m,i) : 0.0;
        if (p3) {
          hd4(m,k,jc,i) = d[i];
          he4(m,k,jc,i) = eos.EnergyFromTemperature(d[i], t[i]);
          hv4(m,k,jc,i) = vi*st;
        } else {
          hd(m,k,i) = d[i];
          he(m,k,i) = eos.EnergyFromTemperature(d[i], t[i]);
          hv(m,k,i) = vi;
        }
      }
    }
    }
  }
  if (root) {
    std::printf("ry_per_accretor: general IC (rank 0): anchor rho %.5g .. %.5g code, "
                "balance residual max %.3e, %d unconverged cells\n", rhoanc_min, rhoanc_max,
                rmax, nfail);
  }
  auto u0 = phd->u0;
  auto phicc = phd->phicc0;
  if (p3) {
    Kokkos::deep_copy(dd4, hd4);
    Kokkos::deep_copy(de4, he4);
    Kokkos::deep_copy(dv4, hv4);
    par_for("ryper_icg3", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dc = dd4(m,k,j,i), v3 = dv4(m,k,j,i);
      u0(m,IDN,k,j,i) = dc;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = dc*v3;
      u0(m,IEN,k,j,i) = de4(m,k,j,i) + 0.5*dc*v3*v3 + dc*phicc(m,k,j,i);
    });
    return;
  }
  Kokkos::deep_copy(dd, hd);
  Kokkos::deep_copy(de, he);
  Kokkos::deep_copy(dv, hv);
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
  // ---- S5 (user 10-09): the physical ambient and the stream's radiation (read only when
  // named, so older inputs' parameter dumps are unchanged)
  {
    auto fatal5 = [&](const std::string &msg) {
      std::cout << "### FATAL ERROR in ry_per_accretor: " << msg << std::endl;
      std::exit(EXIT_FAILURE);
    };
    amb_radeq_ = false;
    if (pin->DoesParameterExist("problem", "env_amb_mode")) {
      const std::string am = pin->GetString("problem", "env_amb_mode");
      amb_radeq_ = (am.compare("radeq") == 0);
      if (!amb_radeq_ && am.compare("hot") != 0) fatal5("env_amb_mode must be hot | radeq");
    }
    strad_ = pin->DoesParameterExist("problem", "stream_rad") &&
             pin->GetBoolean("problem", "stream_rad");
    stprof_ = false;
    if (pin->DoesParameterExist("problem", "stream_t_profile")) {
      const std::string sp = pin->GetString("problem", "stream_t_profile");
      stprof_ = (sp.compare("diffusion") == 0);
      if (!stprof_ && sp.compare("uniform") != 0) {
        fatal5("stream_t_profile must be uniform | diffusion");
      }
      if (stprof_) {
        gen_tsurf_ = pin->GetReal("problem", "stream_t_surf")/eos.temp_cgs;
        if (!gen_ || !(gen_tsurf_ > 0.0 && gen_tsurf_ <= gen_tstr_)) {
          fatal5("stream_t_profile = diffusion needs thermo = general and 0 < "
                 "stream_t_surf <= env_t_stream");
        }
      }
    }
    if ((amb_radeq_ || strad_) && !rad_) {
      fatal5("env_amb_mode = radeq and stream_rad need <rad_m1>");
    }
    if (amb_radeq_) {
      if (env_tro_ != 0.0) {
        fatal5("env_amb_mode = radeq needs problem/env_t_relax = 0 (the floor gas reaches "
               "radiative equilibrium through the M1 coupling, no T relaxation)");
      }
      if (pmbp->pradm1->otab.amask_rho > 0.0) {
        fatal5("env_amb_mode = radeq needs absorption ON: remove <rad_m1>/opac_abs_rho_max");
      }
      if (eos.dfloor != rho_amb) {
        fatal5("env_amb_mode = radeq needs hydro/dfloor = problem/rho_amb (the floor is "
               "applied and booked by this pgen)");
      }
    }
    if (rad_ && (amb_radeq_ || strad_)) {
      // the per-column outer bath: a T_s^4 in the stream window (stream_rad, stream on),
      // 0 elsewhere (vacuum, = the scalar implicit_ebath_x1max default)
      auto *pm1 = pmbp->pradm1;
      if (pm1->iebath_x1max != 0.0) fatal5("stream_rad / radeq need implicit_ebath_x1max = 0");
      strad_e_ = pm1->arad*SQR(SQR(gen_tstr_));
      auto &ind = pmbp->pmesh->mb_indcs;
      const int nb = pmbp->nmb_thispack;
      const int m3 = (ind.nx3 > 1) ? (ind.nx3 + 2*ind.ng) : 1;
      const int m2 = (ind.nx2 > 1) ? (ind.nx2 + 2*ind.ng) : 1;
      Kokkos::realloc(rbath_, nb, m3, m2);
      auto rb = rbath_;
      auto &x3c = pmbp->pcoord->x3v;
      const bool on = strad_ && stream_on_ && !stvol_;   // volume: dark top everywhere
      strad_ar_ = pm1->arad;
      const Real es = strad_e_, phs = phi_s_, ww = nsig_*sig_, sg = sig_;
      const bool stp = stprof_;
      const Real tc = gen_tstr_, tsf = gen_tsurf_, ar = pm1->arad;
      par_for("ryper_bath", DevExeSpace(), 0, nb - 1, 0, m3 - 1, 0, m2 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        const Real dph = WrapPhi(x3c(m,k) - phs);
        const Real t = StreamT(stp, tc, tsf, dph/sg);
        rb(m,k,j) = (on && fabs(dph) < ww) ? (stp ? ar*SQR(SQR(t)) : es) : 0.0;
      });
      if (root && stp) {
        std::printf("ry_per_accretor: stream_t_profile = diffusion: T at x = dphi/sigma = 0, 1, "
                    "2, 2.5, 3: %.0f %.0f %.0f %.0f %.0f K (floor stream_t_surf %.0f K)\n",
                    StreamT(true, tc, tsf, 0.0)*eos.temp_cgs,
                    StreamT(true, tc, tsf, 1.0)*eos.temp_cgs,
                    StreamT(true, tc, tsf, 2.0)*eos.temp_cgs,
                    StreamT(true, tc, tsf, 2.5)*eos.temp_cgs,
                    StreamT(true, tc, tsf, 3.0)*eos.temp_cgs, tsf*eos.temp_cgs);
      }
      pm1->SetX1maxBathColumns(rbath_);
      if (root) {
        std::printf("ry_per_accretor: S5: ambient %s; stream radiation %s (E = a T_s^4 = "
                    "%.6g code at T_s %.1f K, F = (4/3) E v; outer Marshak bath per column)"
                    "\n", amb_radeq_ ? "radeq (floor at radiative equilibrium, absorption "
                    "on, floor booked as Mfl)" : "hot", on ? "INJECTED" : "off (dark)",
                    strad_e_, gen_tstr_*eos.temp_cgs);
      }
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
  // accretor-bigdom-1015: problem/env_r_top_max caps r_top (read only when given): in a
  // domain beyond L1 the Roche potential drops below the top cap again inside the donor's
  // lobe, which the uncapped Phi_wb would take for envelope (psi > 0)
  if (pin->DoesParameterExist("problem", "env_r_top_max")) {
    env_rtop_ = fmin(env_rtop_, pin->GetReal("problem", "env_r_top_max"));
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
    const bool p3 = p3_;
    auto &x2c = pmbp->pcoord->x2v;
    auto &x2fc = pmbp->pcoord->xx2f;
    par_for("ryper_phi", DevExeSpace(), 0, nmb - 1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (p3) {
        // accretor-3d-1017: the exact 3-D Roche potential on the centres and all faces
        const Real st = sin(x2c(m,j)), phk = x3v(m,k);
        const Real pc3 = RochePot3(p, x1v(m,i), st, phk);
        phicc(m,k,j,i) = pc3;
        ph1(m,k,j,i) = RochePot3(p, x1f(m,i), st, phk);
        if (i == n1m1) ph1(m,k,j,i+1) = RochePot3(p, x1f(m,i+1), st, phk);
        ph2(m,k,j,i) = RochePot3(p, x1v(m,i), sin(x2fc(m,j)), phk);
        if (j == n2m1) ph2(m,k,j+1,i) = RochePot3(p, x1v(m,i), sin(x2fc(m,j+1)), phk);
        ph3(m,k,j,i) = RochePot3(p, x1v(m,i), st, x3f(m,k));
        if (k == n3m1) ph3(m,k+1,j,i) = RochePot3(p, x1v(m,i), st, x3f(m,k+1));
        return;
      }
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
    const bool p3 = p3_;
    auto &x2c = pmbp->pcoord->x2v;
    par_for("ryper_phiwb", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real ph = x3v(m,k);
      if (p3) {
        const Real st = sin(x2c(m,j));
        pwc(m,k,j,i) = PhiWB3(p, dlt, rtop, rspin, phtop, x1v(m,i), st, ph);
        pwf(m,k,j,i) = PhiWB3(p, dlt, rtop, rspin, phtop, x1f(m,i), st, ph);
        if (i == n1m1) {
          pwf(m,k,j,i+1) = PhiWB3(p, dlt, rtop, rspin, phtop, x1f(m,i+1), st, ph);
        }
        return;
      }
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
  // accretor-3d-1017 (phys_3d): the same on (m,k,j,i) with the 3-D potentials, plus the
  // explicit g_theta and g_phi (face differences of the true Phi); the ambient anchored at
  // R_acc on the same (theta, phi) ray (sphere) or the photospheric equipotential
  if (p3_) {
    Kokkos::realloc(gr4_, nmb, n3, n2, n1);
    Kokkos::realloc(gt4_, nmb, n3, n2, n1);
    Kokkos::realloc(gp4_, nmb, n3, n2, n1);
    Kokkos::realloc(gpl4_, nmb, n3, n2, n1);
    Kokkos::realloc(cref24_, nmb, n3, n2, n1);
    Kokkos::realloc(rhoamb4_, nmb, n3, n2, n1);
    auto gr = gr4_, gt = gt4_, gp = gp4_, gpl = gpl4_, c2r = cref24_, ra = rhoamb4_;
    auto &x2c = pmbp->pcoord->x2v;
    auto &x2fc = pmbp->pcoord->xx2f;
    const Real ramb0 = env_ramb_, camb2 = env_camb2_, dmin = rho_amb;
    const Real aeq = RochePot(p, racc, 0.5*M_PI);
    par_for("ryper_gre3", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real rl = x1f(m,i), rr = x1f(m,i+1), r = x1v(m,i), ph = x3v(m,k);
      const Real th = x2c(m,j), st = sin(th);
      const Real tl = x2fc(m,j), tr = x2fc(m,j+1);
      const Real wr = PhiWB3(p, dlt, rtop, rspin, phtop, rr, st, ph);
      const Real wl = PhiWB3(p, dlt, rtop, rspin, phtop, rl, st, ph);
      gr(m,k,j,i) = -((RochePot3(p, rr, st, ph) - wr) - (RochePot3(p, rl, st, ph) - wl))
                    /(rr - rl);
      gpl(m,k,j,i) = -(wr - wl)/(rr - rl);
      gt(m,k,j,i) = -(RochePot3(p, r, sin(tr), ph) - RochePot3(p, r, sin(tl), ph))
                    /(r*(tr - tl));
      gp(m,k,j,i) = -(RochePot3(p, r, st, x3f(m,k+1)) - RochePot3(p, r, st, x3f(m,k)))
                    /(r*st*(x3f(m,k+1) - x3f(m,k)));
      const Real psi = phis - PhiWB3(p, dlt, rtop, rspin, phtop, r, st, ph);
      c2r(m,k,j,i) = EnvC2(psi, cph2, np);
      const Real dph = RochePot3(p, r, st, ph) - (eqtop ? aeq : RochePot3(p, racc, st, ph));
      ra(m,k,j,i) = fmax(ramb0*exp(fmax(-dph/camb2, -700.0)), dmin);
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
    if (p3_) {
      Kokkos::realloc(tcol4_, nmb, n3, n2, n1);
      auto tc4 = tcol4_;
      auto &x2c = pmbp->pcoord->x2v;
      par_for("ryper_tcolp3", DevExeSpace(), 0, nmb - 1, 0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real psi = phis - PhiWB3(p, dlt, rtop, rspin, phtop, x1v(m,i), sin(x2c(m,j)),
                                       x3v(m,k));
        tc4(m,k,j,i) = tph*EnvC2(psi, cph2, np)/cph2;
      });
    }
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
    // acc3d-stream-1018 (USER 10-18 "go with a"): problem/env_ic_vz (read only when named,
    // default off = the barotropic column T(psi) on every equipotential, bitwise): the von
    // Zeipel temperature in the atmosphere, T = T_col(psi) (g_ph/g_col(0))^(w/4), g_ph =
    // |grad Phi| of the true Roche potential at the ray's photosphere (3-D: theta, phi;
    // band: the plane), g_col(0) = the same at the column's photosphere (phi 90 cut), w = 1 above the photosphere (psi <= 0), tapering linearly to 0
    // at psi = env_ic_vz_depth c_ph^2 below it (default 3: deeper the barotropic T(psi) of a
    // radiative interior is kept).  Density follows from the IC's own hydrostatic march
    // along each ray (anchored at the photospheric cell at the column's rho), the radeq
    // floor/ambient and RadIC E = a T^4 use the same T, so the IC stays in HSE per ray.
    const bool icvz = pin->DoesParameterExist("problem", "env_ic_vz") &&
                      pin->GetBoolean("problem", "env_ic_vz");
    std::vector<Real> gcolv;
    Real vzdep = 0.0;
    if (icvz) {
      if (!rad_) {
        std::cout << "### FATAL ERROR in ry_per_accretor: problem/env_ic_vz needs <rad_m1>" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      vzdep = pin->GetOrAddReal("problem", "env_ic_vz_depth", 3.0)*env_cph2_;
      gcolv.resize(cpsi.size());
      for (size_t n=0; n<cpsi.size(); ++n) {
        const Real ee = 1.0e-6, rq = crq[n], h = 0.5*M_PI;
        const Real a = (RochePot3(p, rq*(1.0 + ee), 1.0, h) -
                        RochePot3(p, rq*(1.0 - ee), 1.0, h))/(2.0*ee*rq);
        gcolv[n] = std::fabs(a);   // the plane at phi 90: dPhi/dtheta = dPhi/dphi = 0
      }
    }
    Real vzmin = 1.0e300, vzmax = -1.0e300;
    // g_eff is taken at the ray's own PHOTOSPHERE (psi = 0 on this theta, phi ray; v10: the
    // local g of v9 fell to 0.24 g_col in the far ambient, cooling the whole radeq floor by
    // 20-50 %, which is not von Zeipel), so every cell of the ray carries the factor of its
    // surface patch, (g_ph/g_col(0))^(w/4); cached per ray
    Real vzth = -1.0, vzph = -1.0, vzrat = 1.0;
    auto vzray = [&](const Real th, const Real ph) {
      if (th == vzth && ph == vzph) return vzrat;
      vzth = th; vzph = ph; vzrat = 1.0;
      const Real sth = std::sin(th);
      auto psir = [&](const Real r) {
        return phis - (p3_ ? PhiWB3(p, dlt, rtop, rspin, phtop, r, sth, ph)
                           : PhiWB(p, dlt, rtop, rspin, phtop, r, ph));
      };
      Real a = 0.7*racc, b = std::min(1.5*racc, rtop);
      if (!(psir(a) > 0.0 && psir(b) < 0.0)) return vzrat;
      for (int it=0; it<60; ++it) {
        const Real c = 0.5*(a + b);
        if (psir(c) > 0.0) a = c; else b = c;
      }
      const Real r = 0.5*(a + b), ee = 1.0e-6;
      const Real gt = (RochePot3(p, r, std::sin(th + ee), ph) -
                       RochePot3(p, r, std::sin(th - ee), ph))/(2.0*ee*r);
      const Real gr = (RochePot3(p, r*(1.0 + ee), sth, ph) -
                       RochePot3(p, r*(1.0 - ee), sth, ph))/(2.0*ee*r);
      const Real gp = (RochePot3(p, r, sth, ph + ee) - RochePot3(p, r, sth, ph - ee))
                      /(2.0*ee*r*sth);
      // g_col at psi = 0 (linear in the column's psi)
      const size_t n = std::upper_bound(cpsi.begin(), cpsi.end(), 0.0) - cpsi.begin();
      Real gc = gcolv.back();
      if (n == 0) {
        gc = gcolv.front();
      } else if (n < cpsi.size()) {
        const Real x = (0.0 - cpsi[n-1])/(cpsi[n] - cpsi[n-1]);
        gc = (1.0 - x)*gcolv[n-1] + x*gcolv[n];
      }
      vzrat = std::sqrt(gt*gt + gr*gr + gp*gp)/gc;
      return vzrat;
    };
    auto vzfac = [&](const Real r, const Real th, const Real ph, const Real psi) {
      (void) r;
      if (!icvz) return static_cast<Real>(1.0);
      const Real w = (psi <= 0.0) ? 1.0 : std::max(0.0, 1.0 - psi/vzdep);
      if (!(w > 0.0)) return static_cast<Real>(1.0);
      const Real f = std::pow(vzray(th, ph), 0.25*w);
      vzmin = std::min(vzmin, f);
      vzmax = std::max(vzmax, f);
      return f;
    };
    for (int m=0; m<nmb; ++m) {
      for (int k=0; k<n3; ++k) {
        for (int i=0; i<n1; ++i) {
          const Real psi = phis - PhiWB(p, dlt, rtop, rspin, phtop, hx1(m,i), hx3(m,k));
          Real lt, lr;
          ColumnAt(cpsi, clt, clr, psi, lt, lr);
          ht(m,k,i) = std::exp(lt)/tk;
          if (icvz) ht(m,k,i) *= vzfac(hx1(m,i), 0.5*M_PI, hx3(m,k), psi);
        }
      }
    }
    Kokkos::deep_copy(tcol_, ht);
    if (p3_) {
      // accretor-3d-1017: T(psi) of the column on the 3-D equipotentials
      Kokkos::realloc(tcol4_, nmb, n3, n2, n1);
      auto ht4 = Kokkos::create_mirror_view(tcol4_);
      auto hx2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
      for (int m=0; m<nmb; ++m) {
        for (int k=0; k<n3; ++k) {
          for (int j=0; j<n2; ++j) {
            const Real st = std::sin(hx2(m,j));
            for (int i=0; i<n1; ++i) {
              const Real psi = phis - PhiWB3(p, dlt, rtop, rspin, phtop, hx1(m,i), st,
                                             hx3(m,k));
              Real lt, lr;
              ColumnAt(cpsi, clt, clr, psi, lt, lr);
              ht4(m,k,j,i) = std::exp(lt)/tk;
              if (icvz) ht4(m,k,j,i) *= vzfac(hx1(m,i), hx2(m,j), hx3(m,k), psi);
            }
          }
        }
      }
      Kokkos::deep_copy(tcol4_, ht4);
    }
    if (icvz && root) {
      std::printf("ry_per_accretor: env_ic_vz: von Zeipel IC T x (g/g_col)^(w/4), taper depth "
                  "%.4g (psi), factor %.5f .. %.5f (rank 0)\n", vzdep, vzmin, vzmax);
    }
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

//----------------------------------------------------------------------------------------
//! acctrc-1019c (USER 10-19, Caltech NOTE eb0b2c7a): the CUMULATIVE HISTORY REGISTERS in
//! the restart file, so that Min..Jin, Mfl/Mmsk/Emsk, Ebnd/Erel, Mnzf/Eirr (user.hst) and
//! the tracer sums Nz0..Bin (trc.hst) continue across a restart and a restarted run is
//! bitwise the continuous one.  Carried in the generic marked pgen state block
//! (kPgenRstMagic, pgen.hpp; ProblemGenerator::pgen_rst_write_func), so the hydro records
//! and every other pgen's file are untouched.  Layout (all ranks write the same bytes):
//!   char[8] "RYPACC01", int32 nranks, int32 kRstPer, int32 kRstGlob, int32 pad,
//!   nranks x kRstPer Reals (rank-ordered: acc0_[13], hist_prev_[13], accx0_[2], erinj_,
//!   trc0_[9]; the registers are PER RANK partial sums, MPI-summed at each output, so they
//!   are stored per rank to keep the sums bitwise), then kRstGlob Reals (rank 0's
//!   hist_tprev_, the Picard passes/solves since the last output, Picmax).
//! The stage registers acc1_/accx1_/trc1_ are copied from acc0_ at stage 1 and need no
//! storing; stage_ctr_ is a multiple of nstages_ at a cycle end.  With the per-column outer
//! bath (<rad_m1> iebath_hi_col_on: env_amb_mode = radeq or stream_rad) the block also holds
//! the module's ibadv_x1max (the previous solve's advective radiation flux through r_out,
//! read by the first hydro stages of the next cycle into Ebnd, register 10; not M1 state and
//! not in the rst otherwise): int32 hdr[3] = nk*nj per MeshBlock, then nmb_total x nk*nj
//! Reals in global MeshBlock order (restored by gid, so also onto other rank counts).
//! Without it the restarted Ebnd is offset by the first cycle's advective term.  Restart:
//!   * same nranks: every register restored exactly (bitwise continuation);
//!   * other nranks: the per-rank values are summed onto rank 0 (continuous, not bitwise);
//!   * a file without the block (any pre-1019c rst): registers from 0, as before (logged).
//! Picmax is GLOBAL (the largest pass count since t = 0, max(stored, <rad_m1> impl_itmax)).
//! problem/rst_registers = false (read only when named) writes no block (the 4960ac49 file).
constexpr char kRyRstMagic[8] = {'R', 'Y', 'P', 'A', 'C', 'C', '0', '1'};
constexpr int kRstPer = 2*kNaccR + 2 + 1 + kNtr;   // 38
constexpr int kRstGlob = 4;
MeshBlockPack *rst_pmbp_ = nullptr;
std::vector<Real> rst_badv_;   // ibadv_x1max read from the restart file (gid order), applied
int rst_badv_n_ = -1;          // after EnvSetup (RyPerRstBadv); nk*nj per MB (-1: no block)

std::vector<char> RyPerRstWrite() {
  Real loc[kRstPer];
  int o = 0;
  for (int n=0; n<kNaccR; ++n) loc[o++] = acc0_[n];
  for (int n=0; n<kNaccR; ++n) loc[o++] = hist_prev_[n];
  for (int n=0; n<2; ++n) loc[o++] = accx0_[n];
  loc[o++] = erinj_;
  for (int n=0; n<kNtr; ++n) loc[o++] = trc0_[n];
  const int nr = global_variable::nranks;
  radm1::RadiationM1 *pm1 = (rst_pmbp_ != nullptr) ? rst_pmbp_->pradm1 : nullptr;
  // the outer bath's advective flux (see above), gathered in gid order
  int nbv = 0;
  std::vector<Real> bv;
  if (pm1 != nullptr && pm1->iebath_hi_col_on && pm1->ibadv_x1max.extent_int(0) > 0) {
    auto &bd = pm1->ibadv_x1max;
    nbv = bd.extent_int(1)*bd.extent_int(2);
    Mesh *pmsh = rst_pmbp_->pmesh;
    const int nmbt = pmsh->nmb_total, nml = bd.extent_int(0);
    auto hb = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), bd);
    std::vector<Real> lb(static_cast<size_t>(nml)*nbv);
    const int n1 = bd.extent_int(1), n2 = bd.extent_int(2);
    for (int m=0; m<nml; ++m) for (int k=0; k<n1; ++k) for (int j=0; j<n2; ++j)
      lb[(static_cast<size_t>(m)*n1 + k)*n2 + j] = hb(m,k,j);
    bv.assign(static_cast<size_t>(nmbt)*nbv, 0.0);
#if MPI_PARALLEL_ENABLED
    std::vector<int> cnt(nr), dsp(nr);
    for (int r=0; r<nr; ++r) {
      cnt[r] = pmsh->nmb_eachrank[r]*nbv;
      dsp[r] = pmsh->gids_eachrank[r]*nbv;
    }
    MPI_Allgatherv(lb.data(), nml*nbv, MPI_ATHENA_REAL, bv.data(), cnt.data(), dsp.data(),
                   MPI_ATHENA_REAL, MPI_COMM_WORLD);
#else
    for (size_t n=0; n<lb.size(); ++n) bv[n] = lb[n];
#endif
  }
  std::vector<Real> all(static_cast<size_t>(nr)*kRstPer + kRstGlob + bv.size(), 0.0);
#if MPI_PARALLEL_ENABLED
  MPI_Allgather(loc, kRstPer, MPI_ATHENA_REAL, all.data(), kRstPer, MPI_ATHENA_REAL,
                MPI_COMM_WORLD);
#else
  for (int n=0; n<kRstPer; ++n) all[n] = loc[n];
#endif
  Real glb[kRstGlob] = {hist_tprev_, 0.0, 0.0, picmax_rst_};
  if (pm1 != nullptr) {
    glb[1] = pm1->impl_nstep - rad_pn0_;
    glb[2] = pm1->impl_itsum - rad_ps0_;
    glb[3] = fmax(picmax_rst_, pm1->impl_itmax);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(glb, kRstGlob, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif
  for (int n=0; n<kRstGlob; ++n) all[static_cast<size_t>(nr)*kRstPer + n] = glb[n];
  for (size_t n=0; n<bv.size(); ++n) all[static_cast<size_t>(nr)*kRstPer + kRstGlob + n] = bv[n];
  const std::int32_t hdr[4] = {static_cast<std::int32_t>(nr), kRstPer, kRstGlob,
                               static_cast<std::int32_t>(nbv)};
  std::vector<char> out(sizeof(kRyRstMagic) + sizeof(hdr) + all.size()*sizeof(Real));
  std::memcpy(out.data(), &(kRyRstMagic[0]), sizeof(kRyRstMagic));
  std::memcpy(out.data() + sizeof(kRyRstMagic), &(hdr[0]), sizeof(hdr));
  std::memcpy(out.data() + sizeof(kRyRstMagic) + sizeof(hdr), all.data(),
              all.size()*sizeof(Real));
  return out;
}

void RyPerRstRestore(const std::vector<char> &blk) {
  const bool root = (global_variable::my_rank == 0);
  const size_t h = sizeof(kRyRstMagic) + 4*sizeof(std::int32_t);
  if (blk.size() < h || std::memcmp(blk.data(), &(kRyRstMagic[0]), sizeof(kRyRstMagic))) {
    if (root) {
      std::printf("ry_per_accretor: this restart file carries no history registers (pre-"
                  "acctrc-1019c): the cumulative user.hst/trc.hst columns restart from 0\n");
    }
    return;
  }
  std::int32_t hdr[4];
  std::memcpy(&(hdr[0]), blk.data() + sizeof(kRyRstMagic), sizeof(hdr));
  const int nrf = hdr[0], nbv = hdr[3];
  const size_t nreg = static_cast<size_t>(nrf)*kRstPer + kRstGlob;
  const size_t nrest = blk.size() - h;
  const size_t nbvt = (nbv > 0 && nrest > nreg*sizeof(Real)) ?
                      (nrest/sizeof(Real) - nreg) : 0;
  if (nrf < 1 || hdr[1] != kRstPer || hdr[2] != kRstGlob || nbv < 0 ||
      (nbv > 0 && nbvt % nbv != 0) || blk.size() != h + (nreg + nbvt)*sizeof(Real)) {
    if (root) {
      std::printf("ry_per_accretor: WARNING the restart file's history-register block has "
                  "an unknown layout (%d ranks x %d + %d, %zu bytes): ignored, the "
                  "cumulative columns restart from 0\n", nrf, hdr[1], hdr[2], blk.size());
    }
    return;
  }
  std::vector<Real> all(nreg + nbvt);
  std::memcpy(all.data(), blk.data() + h, all.size()*sizeof(Real));
  rst_badv_n_ = nbv;
  rst_badv_.assign(all.begin() + nreg, all.end());
  const int nr = global_variable::nranks, me = global_variable::my_rank;
  Real loc[kRstPer];
  if (nrf == nr) {
    for (int n=0; n<kRstPer; ++n) loc[n] = all[static_cast<size_t>(me)*kRstPer + n];
  } else {
    for (int n=0; n<kRstPer; ++n) {
      Real sum = 0.0;
      for (int r=0; r<nrf; ++r) sum += all[static_cast<size_t>(r)*kRstPer + n];
      loc[n] = root ? sum : 0.0;
    }
  }
  int o = 0;
  for (int n=0; n<kNaccR; ++n) acc0_[n] = acc1_[n] = loc[o++];
  for (int n=0; n<kNaccR; ++n) hist_prev_[n] = loc[o++];
  for (int n=0; n<2; ++n) accx0_[n] = accx1_[n] = loc[o++];
  erinj_ = loc[o++];
  for (int n=0; n<kNtr; ++n) trc0_[n] = trc1_[n] = loc[o++];
  const Real *glb = all.data() + static_cast<size_t>(nrf)*kRstPer;
  hist_tprev_ = glb[0];
  picmax_rst_ = glb[3];
  if (root) {
    radm1::RadiationM1 *pm1 = (rst_pmbp_ != nullptr) ? rst_pmbp_->pradm1 : nullptr;
    const Real pn = (pm1 != nullptr) ? pm1->impl_nstep : 0.0;
    const Real ps = (pm1 != nullptr) ? pm1->impl_itsum : 0.0;
    rad_pn0_ = pn - glb[1];
    rad_ps0_ = ps - glb[2];
    std::printf("ry_per_accretor: history registers restored from the restart file (%d "
                "ranks%s)\n", nrf, (nrf == nr) ? ", bitwise continuation" :
                ": summed onto rank 0, continuous but not bitwise");
  }
}

//! the stored ibadv_x1max into <rad_m1> (after EnvSetup: SetX1maxBathColumns zeroes it)
void RyPerRstBadv() {
  radm1::RadiationM1 *pm1 = (rst_pmbp_ != nullptr) ? rst_pmbp_->pradm1 : nullptr;
  const bool on = (pm1 != nullptr) && pm1->iebath_hi_col_on &&
                  (pm1->ibadv_x1max.extent_int(0) > 0);
  if (rst_badv_n_ < 0) return;                  // no block (pre-1019c file): logged already
  if (!on && rst_badv_.empty()) return;
  const bool root = (global_variable::my_rank == 0);
  Mesh *pmsh = rst_pmbp_->pmesh;
  if (!on || rst_badv_.empty() ||
      rst_badv_n_ != pm1->ibadv_x1max.extent_int(1)*pm1->ibadv_x1max.extent_int(2) ||
      rst_badv_.size() != static_cast<size_t>(pmsh->nmb_total)*rst_badv_n_) {
    if (root) {
      std::printf("ry_per_accretor: WARNING the restart file's outer-bath advective flux "
                  "(%zu values) does not match this run: not restored (Ebnd offset by the "
                  "first cycle's advective term)\n", rst_badv_.size());
    }
    rst_badv_.clear();
    return;
  }
  auto &bd = pm1->ibadv_x1max;
  auto hb = Kokkos::create_mirror_view(bd);
  const int g0 = pmsh->gids_eachrank[global_variable::my_rank];
  const int n1 = bd.extent_int(1), n2 = bd.extent_int(2);
  for (int m=0; m<bd.extent_int(0); ++m) for (int k=0; k<n1; ++k) for (int j=0; j<n2; ++j)
    hb(m,k,j) = rst_badv_[(static_cast<size_t>(g0 + m)*n1 + k)*n2 + j];
  Kokkos::deep_copy(bd, hb);
  rst_badv_.clear();
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
  // accretor-bigdom-1015: stream_inject = boundary (default: the r_out window) | volume
  // (the nozzle at stream_r inside the domain); keys read only when given
  stvol_ = false;
  stl1_ = false;
  nz_if_ = -1;
  if (pin->DoesParameterExist("problem", "stream_inject")) {
    const std::string si = pin->GetString("problem", "stream_inject");
    stl1_ = (si.compare("l1") == 0);
    stvol_ = (si.compare("volume") == 0) || stl1_;
    if (!stvol_ && si.compare("boundary") != 0) {
      fatal("stream_inject must be boundary | volume | l1");
    }
  }
  str_r_ = r_out;
  if (pin->DoesParameterExist("problem", "stream_r")) {
    str_r_ = pin->GetReal("problem", "stream_r");
    if (!stvol_ && str_r_ != r_out) fatal("stream_r != x1max needs stream_inject = volume");
  }
  // (stream_inject = l1: the trajectory only gives d_L1 and the reported path; the
  // reported crossing radius is then R_acc)
  BallisticOut bo = IntegrateStream(mu, csd/aom, (stl1_ ? racc : str_r_)/asep);
  BallisticOut bi = IntegrateStream(mu, csd/aom, racc/asep);
  const Real width = pin->GetOrAddReal("problem", "stream_width", csd/rp_.omega);
  if (!bo.ok_cross) fatal("the ballistic stream never reaches stream_r (default x1max)");
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
  sig_ = width/str_r_;
  // stream_inject = l1 (bigdom-1015, user 10-15): the nozzle sits just inside L1 with the
  // Lubow & Shu (1975) state: isothermal, sonic flow along -x at c_T (problem/stream_l1_v,
  // default the launch speed c_s,don of t_don), uniform T (problem/stream_l1_t [K], default
  // T_L1 155300 K of stream_physics_1009), Gaussian across the stream with the in-plane
  // hydrostatic width sigma_y = c_T/(Omega sqrt(B)), B = A - 1, A = (1 - mu)/d_a^3 + mu/d_d^3
  // the potential curvature at L1 (a = 1 units; problem/stream_l1_width overrides), and
  // the peak density from the band's mass rate problem/stream_mdot (code mass/time, the
  // discrete face flux of the window, set after the nozzle faces below)
  Real l1a = 0.0, l1b = 0.0;
  if (stl1_) {
    const Real da = bo.dl1, dd = 1.0 - bo.dl1;
    l1a = (1.0 - mu)/(da*da*da) + mu/(dd*dd*dd);
    l1b = l1a - 1.0;
    l1_c_ = pin->GetOrAddReal("problem", "stream_l1_v", csd);
    const Real wl1 = pin->GetOrAddReal("problem", "stream_l1_width",
                                       l1_c_/(rp_.omega*std::sqrt(l1b)));
    str_r_ = bo.dl1*asep;
    sig_ = wl1/str_r_;
    phi_s_ = 0.0;
    vr_s_ = -l1_c_;
    vp_s_ = 0.0;
    if (gen_) {
      l1_t_ = pin->GetOrAddReal("problem", "stream_l1_t", 155300.0)/eos.temp_cgs;
    }
  }
  // accretor-3d-1016: 3-D mode from the theta range (the 2-D band is +-0.35 deg)
  {
    const auto &ms = pmy_mesh_->mesh_size;
    g3d_ = (std::fmax(std::fabs(ms.x2min - 0.5*M_PI), std::fabs(ms.x2max - 0.5*M_PI))
            > 2.0*M_PI/180.0);
    sz_ = 0.0;
    if (g3d_ && (stvol_ || pin->DoesParameterExist("problem", "donor_mask"))) {
      if (stvol_) {
        if (!pin->DoesParameterExist("problem", "stream_sigma_z") ||
            !pin->DoesParameterExist("problem", "stream_mdot")) {
          fatal("3-D (wide theta) nozzle needs problem/stream_sigma_z [Rsun] and "
                "problem/stream_mdot (the FULL stream rate, code mass/time)");
        }
        sz_ = pin->GetReal("problem", "stream_sigma_z");
        if (!(sz_ > 0.0)) fatal("problem/stream_sigma_z must be > 0");
      }
      if (root) {
        std::printf("ry_per_accretor: 3-D mode (theta %.3f..%.3f deg): nozzle Gaussian in z, "
                    "exact 3-D Roche-lobe mask; forces/IC/WB stay theta-invariant (band)\n",
                    ms.x2min*180.0/M_PI, ms.x2max*180.0/M_PI);
      }
    }
    // accretor-3d-1017: problem/phys_3d (the 3-D physics; see the namespace header at p3_),
    // problem/theta_bc, problem/theta_sponge_deg / theta_sponge_tau; all read only when
    // given, so older decks (2-D and the 049c1f39 3-D benchmark) keep their parameter dump
    p3_ = false;
    thout_ = false;
    thhse_ = false;
    thsp_w_ = 0.0;
    thsp_tau_ = 0.0;
    if (pin->DoesParameterExist("problem", "phys_3d")) {
      p3_ = pin->GetBoolean("problem", "phys_3d");
      if (p3_ && !g3d_) {
        fatal("problem/phys_3d needs the 3-D mode (theta range beyond 2 deg from the plane)");
      }
      if (p3_ && !(env_ && gen_)) {
        fatal("problem/phys_3d needs problem/inner = envelope and thermo = general");
      }
    }
    if (pin->DoesParameterExist("problem", "theta_bc")) {
      const std::string tb = pin->GetString("problem", "theta_bc");
      thhse_ = (tb.compare("hse") == 0);
      thout_ = (tb.compare("outflow") == 0) || thhse_;
      if (!thout_ && tb.compare("reflect") != 0) {
        fatal("theta_bc must be reflect | outflow | hse");
      }
      if (thout_ && !p3_) fatal("theta_bc = outflow | hse needs problem/phys_3d = true");
    }
    thwc_ = false;
    if (pin->DoesParameterExist("problem", "theta_wall_closed")) {
      thwc_ = pin->GetBoolean("problem", "theta_wall_closed");
      if (thwc_ && !thhse_) fatal("problem/theta_wall_closed needs theta_bc = hse");
      if (thwc_ && root) {
        std::printf("ry_per_accretor: theta_wall_closed: the hse theta-wall faces carry "
                    "only the x2 momentum flux (interface pressure); mass, energy, x1/x3 "
                    "momentum and scalar fluxes are zero there\n");
      }
    }
    {
      const bool u2 = (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::user) ||
                      (pmy_mesh_->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::user);
      const bool u2b = (pmy_mesh_->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::user) &&
                       (pmy_mesh_->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::user);
      if (thout_ && !u2b) fatal("theta_bc = outflow | hse needs mesh/ix2_bc = ox2_bc = user");
      if (!thout_ && u2) fatal("mesh ix2_bc/ox2_bc = user needs problem/theta_bc = outflow | hse");
    }
    if (pin->DoesParameterExist("problem", "theta_sponge_deg")) {
      thsp_w_ = pin->GetReal("problem", "theta_sponge_deg")*M_PI/180.0;
      if (thsp_w_ > 0.0) {
        if (!p3_) fatal("theta_sponge_deg > 0 needs problem/phys_3d = true");
        thsp_tau_ = pin->GetReal("problem", "theta_sponge_tau");
        if (!(thsp_tau_ > 0.0)) fatal("problem/theta_sponge_tau must be > 0");
      }
    }
    if (root && p3_) {
      std::printf("ry_per_accretor: phys_3d: 3-D Roche potential (gravity r/theta/phi, Phi_wb, "
                  "etotgrav), full Coriolis, IC/WB/ambient on 3-D equipotentials; theta walls "
                  "%s; theta sponge %s (width %.3f deg, tau %.4g)\n",
                  thhse_ ? "HSE (closed, hydrostatic theta ghosts, v_theta mirrored)" :
                  (thout_ ? "OUTFLOW (hydrostatic theta ghosts, v_theta diode)" : "reflecting"),
                  (thsp_w_ > 0.0) ? "on" : "off", thsp_w_*180.0/M_PI, thsp_tau_);
    }
  }
  // the nozzle shell: stream_ncell cells above the x1 face nearest stream_r
  msk_ = false;
  if (stvol_ || pin->DoesParameterExist("problem", "donor_mask")) {
    if (!env_ || !gen_) fatal("stream_inject = volume / donor_mask need inner = envelope "
                              "and thermo = general");
  }
  if (stvol_) {
    if (pmy_mesh_->mesh_indcs.nx1 != pmy_mesh_->mb_indcs.nx1) {
      fatal("stream_inject = volume needs one MeshBlock in x1");
    }
    const int nzc = pin->GetOrAddInteger("problem", "stream_ncell", 2);
    auto hxf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
    const int is = pmy_mesh_->mb_indcs.is, ie = pmy_mesh_->mb_indcs.ie;
    int fb = is;
    for (int i=is; i<=ie+1; ++i) {
      if (std::fabs(hxf(0,i) - str_r_) < std::fabs(hxf(0,fb) - str_r_)) fb = i;
    }
    if (stl1_) fb -= nzc;     // l1: the nozzle's OUTER face is the face nearest d_L1
    if (!(nzc >= 1 && fb > is && fb + nzc <= ie + 1)) {
      fatal("the nozzle (stream_r, stream_ncell) must lie inside the domain");
    }
    nz_rlo_ = hxf(0,fb);
    nz_rhi_ = hxf(0,fb + nzc);
    nz_if_ = fb;
    if (stl1_ || pin->DoesParameterExist("problem", "stream_mdot")) {
      // peak density from the band's mass rate (problem/stream_mdot; required for l1, optional
      // for volume, where it replaces rho_stream): the discrete sum over the window cells of
      // the inner face's rho |v_r| A (A = r^2 (cos th1 - cos th2) dphi), so a supersonic
      // nozzle delivers exactly stream_mdot (upwind face flux)
      if (!pin->DoesParameterExist("problem", "stream_mdot")) {
        fatal("stream_inject = l1 needs problem/stream_mdot (code mass/time in the band)");
      }
      const Real mdot = pin->GetReal("problem", "stream_mdot");
      const auto &ms = pmy_mesh_->mesh_size;
      const int nk = pmy_mesh_->mesh_indcs.nx3;
      const Real dphk = (ms.x3max - ms.x3min)/nk, rf = nz_rlo_;
      const Real ath = rf*rf*(std::cos(ms.x2min) - std::cos(ms.x2max))*dphk;
      Real sum = 0.0;
      if (g3d_) {
        // 3-D: stream_mdot is the FULL stream; the discrete sum over the window cells
        // (|dphi| < nsig sigma_phi, |z| < nsig sigma_z, z = r_f cos theta) of the inner
        // face's e^{-q^2/2} |v_r| A, A = r_f^2 (cos th_j - cos th_j+1) dphi
        if (stl1_) fatal("stream_inject = l1 is 2-D only (use volume in 3-D)");
        const int nj = pmy_mesh_->mesh_indcs.nx2;
        const Real dth = (ms.x2max - ms.x2min)/nj;
        Real sum0 = 0.0;          // per unit height at z = 0 (central slice), for the report
        for (int k=0; k<nk; ++k) {
          const Real phk = ms.x3min + (k + 0.5)*dphk;
          const Real dph = WrapPhi(phk - phi_s_);
          if (std::fabs(dph) >= nsig_*sig_) continue;
          sum0 += std::exp(-0.5*SQR(dph/sig_))*std::fabs(vr_s_)*rf*dphk;
          for (int j=0; j<nj; ++j) {
            const Real t0 = ms.x2min + j*dth, t1 = t0 + dth;
            const Real zj = rf*std::cos(0.5*(t0 + t1));
            if (std::fabs(zj) >= nsig_*sz_) continue;
            sum += std::exp(-0.5*(SQR(dph/sig_) + SQR(zj/sz_)))*std::fabs(vr_s_)*
                   rf*rf*(std::cos(t0) - std::cos(t1))*dphk;
          }
        }
        rho_s_ = mdot/sum;
        if (global_variable::my_rank == 0) {
          std::printf("ry_per_accretor: 3-D nozzle: sigma_z %.4f Rsun (|z| < %.4f, theta "
                      "%.3f..%.3f deg), full rate stream_mdot %.6g code -> peak rho %.6g "
                      "code; central slice %.6g code/time per Rsun of height = %.6g in the "
                      "2-D band of height r_f 0.0122718 (r_f %.5f)\n", sz_, nsig_*sz_,
                      90.0 - std::asin(std::fmin(nsig_*sz_/rf, 1.0))*180.0/M_PI,
                      90.0 + std::asin(std::fmin(nsig_*sz_/rf, 1.0))*180.0/M_PI, mdot,
                      rho_s_, rho_s_*sum0, rho_s_*sum0*rf*0.0122718, rf);
        }
        sum = 0.0;
      }
      for (int k=0; k<nk; ++k) {
        const Real phk = ms.x3min + (k + 0.5)*dphk;
        const Real dph = WrapPhi(phk - phi_s_);
        if (std::fabs(dph) < nsig_*sig_) {
          const Real vr = stl1_ ? l1_c_*std::cos(phk) : std::fabs(vr_s_);
          sum += std::exp(-0.5*SQR(dph/sig_))*vr*ath;
        }
      }
      if (!g3d_) rho_s_ = mdot/sum;
    }
  }
  // the two stabilisers added for the L1 nozzle (10-15), keys default: on for l1, off for
  // volume (stream_back_sink); on with any nozzle/mask (floor_keep_momentum)
  nzback_ = stvol_ && pin->GetOrAddBoolean("problem", "stream_back_sink", true);
  // accretor-sector-light-1017 (user 10-17): problem/stream_sector_light (read, not added:
  // the key stays out of the PAR_DUMP of runs that do not set it; default full = v4)
  ssl_ = 0;
  if (pin->DoesParameterExist("problem", "stream_sector_light")) {
    const std::string sl = pin->GetString("problem", "stream_sector_light");
    if (sl.compare("full") == 0) {
      ssl_ = 0;
    } else if (sl.compare("scaled") == 0) {
      ssl_ = 1;
    } else if (sl.compare("surface") == 0) {
      ssl_ = 2;
    } else if (sl.compare("dark") == 0) {
      ssl_ = 3;
    } else if (sl.compare("neck") == 0) {
      ssl_ = 4;
    } else {
      fatal("stream_sector_light must be full | scaled | surface | dark | neck");
    }
    if (root) {
      std::printf("ry_per_accretor: stream_sector_light %s (back-sink sector radiation)\n",
                  sl.c_str());
    }
  }
  dtrep_ = pin->DoesParameterExist("problem", "dt_report") &&
           pin->GetBoolean("problem", "dt_report");
  flmom_ = (stvol_ || pin->DoesParameterExist("problem", "donor_mask")) &&
           pin->GetOrAddBoolean("problem", "floor_keep_momentum", true);
  if (pin->DoesParameterExist("problem", "donor_mask")) {
    const std::string dm = pin->GetString("problem", "donor_mask");
    msk_ = (dm.compare("roche") == 0);
    if (!msk_ && dm.compare("none") != 0) fatal("donor_mask must be none | roche");
    msk_xl1_ = bo.dl1*asep;
    msk_phl1_ = RochePot(rp_, msk_xl1_, 0.0);
  }
  // problem/donor_irr (default false; needs donor_mask = roche and <rad_m1>): the masked
  // lobe holds E = a T_d^4, F = 0 every stage, T_d = problem/donor_teff [K] (default 26200,
  // the donor T_eff of the Wade+2026 MESA track used by the S5 deck)
  mirr_ = pin->DoesParameterExist("problem", "donor_irr") &&
          pin->GetBoolean("problem", "donor_irr");
  if ((ssl_ == 1 || ssl_ == 2 || ssl_ == 4) && !mirr_) {
    fatal("stream_sector_light = scaled | surface needs donor_irr = true");
  }
  if (mirr_) {
    if (!msk_ || !rad_) fatal("donor_irr needs donor_mask = roche and <rad_m1>");
    // T_d from the donor luminosity problem/donor_lum [Lsun] (default 10^4.84 = 69183, the
    // Wade+2026 track value of the S5 deck) over the Roche lobe's volume-equivalent sphere,
    // 4 pi R_L^2 (Eggleton 1983 at q = M_d/M_a): the per-area flux of the full 3-D donor,
    // which the 2-D band samples (the band's luminosities scale by 4 pi / band solid angle,
    // as Lmeas/Ltop).  problem/donor_teff [K], if given, overrides.
    const Real qd = md/ma, q3 = std::cbrt(qd);
    const Real rlv = asep*0.49*q3*q3/(0.6*q3*q3 + std::log(1.0 + q3));
    const Real ldon = pin->GetOrAddReal("problem", "donor_lum", std::pow(10.0, 4.84));
    Real tdk = std::pow(ldon*3.828e33/(4.0*M_PI*SQR(rlv*kRsun)*5.6704e-5), 0.25);
    if (pin->DoesParameterExist("problem", "donor_teff")) {
      tdk = pin->GetReal("problem", "donor_teff");
    }
    if (root) {
      std::printf("ry_per_accretor: donor_irr: L_d %.6g Lsun, R_L (Eggleton, volume) %.4f "
                  "Rsun -> T_d %.1f K%s\n", ldon, rlv, tdk,
                  pin->DoesParameterExist("problem", "donor_teff") ? " (donor_teff key)" : "");
    }
    const Real tunit = pmbp->phydro->peos->eos_data.temp_cgs;
    const Real tdc = tdk/tunit;
    mirr_e_ = pmbp->pradm1->arad*SQR(SQR(tdc));
    if (root) {
      std::printf("ry_per_accretor: donor_irr: the donor lobe radiates at T_d %.1f K (E = "
                  "a T^4 = %.6g code, F = 0); booked in Ebnd as radiation in\n", tdk,
                  mirr_e_);
    }
    // ---- v4 (user 10-16) gravity darkening, problem/donor_gdark = beta (default 0.25,
    // von Zeipel, radiative envelope; 0 = the uniform T_d above): T_loc^4 = C g^(4 beta),
    // g = |grad Phi| (rotating-frame Roche) on the L1 lobe surface, floored at
    // donor_gdark_tfloor (0.3) x T_mean, C from  int sigma T_loc^4 dA = L_d  over the FULL
    // 3-D lobe surface (also in the 2-D band: the band holds the equatorial cut of the 3-D
    // donor; T_mean^4 = L_d/(sigma A_lobe)).  With donor_teff given, L_d = 4 pi R_L^2
    // sigma T_d^4.  The surface integral and the in-plane line integral are by rays from
    // the donor centre (DonorSurfT), polar axis toward L1.
    mirr_beta_ = pin->GetOrAddReal("problem", "donor_gdark", 0.25);
    const Real tfr = pin->GetOrAddReal("problem", "donor_gdark_tfloor", 0.3);
    const RocheParams rpp = rp_;
    const Real q4 = SQR(SQR(tdk))*4.0*M_PI*SQR(rlv);   // L_d/sigma [K^4 Rsun^2]
    const int na = 720, nq = 48;   // alpha (0..pi), psi quarter (0..pi/2), symmetric in y, z
    std::vector<Real> ga(na*nq), da(na*nq);
    std::vector<Real> gl(2*na), dl(2*na);
    Real atot = 0.0, ltot = 0.0;
    for (int ia=0; ia<na; ++ia) {
      const Real al = (ia + 0.5)*M_PI/na, sa = std::sin(al), ca = std::cos(al);
      for (int iq=0; iq<nq; ++iq) {
        const Real ps = (iq + 0.5)*0.5*M_PI/nq;
        const Real nx = -ca, ny = sa*std::cos(ps), nz = sa*std::sin(ps);
        const Real t = DonorSurfT(rpp, msk_phl1_, msk_xl1_, nx, ny, nz, 0.0);
        Real gx, gy, gz;
        const Real g = (t > 0.0) ? RocheGeffXYZ(rpp, rpp.asep + t*nx, t*ny, t*nz,
                                                &gx, &gy, &gz) : 0.0;
        const Real cf = (t > 0.0) ? std::fmax(std::fabs(nx*gx + ny*gy + nz*gz)/
                                              std::fmax(g, 1.0e-30), 0.05) : 1.0;
        ga[ia*nq + iq] = g;
        da[ia*nq + iq] = 4.0*SQR(t)*sa*(M_PI/na)*(0.5*M_PI/nq)/cf;
        atot += da[ia*nq + iq];
      }
      // the in-plane (z = 0) cut, both y signs (symmetric): the band's 2-D source
      const Real nx = -ca, ny = sa;
      const Real t = DonorSurfT(rpp, msk_phl1_, msk_xl1_, nx, ny, 0.0, 0.0);
      Real gx, gy, gz;
      const Real g = (t > 0.0) ? RocheGeffXYZ(rpp, rpp.asep + t*nx, t*ny, 0.0,
                                              &gx, &gy, &gz) : 0.0;
      const Real cf = (t > 0.0) ? std::fmax(std::fabs(nx*gx + ny*gy)/std::fmax(g, 1.0e-30),
                                            0.05) : 1.0;
      gl[ia] = g;
      dl[ia] = 2.0*t*(M_PI/na)/cf;
      ltot += dl[ia];
    }
    const Real tm4 = q4/atot, fl4 = SQR(SQR(tfr))*tm4;
    const Real b4 = 4.0*mirr_beta_;
    Real cc = tm4;   // C [K^4 per g^(4 beta)], solved for int max(C g^4b, floor) dA = q4
    if (mirr_beta_ > 0.0) {
      Real gm = 0.0;
      for (int n=0; n<na*nq; ++n) gm += ga[n]*da[n];
      gm /= atot;
      Real clo = 1.0e-6*tm4/std::pow(gm, b4), chi = 1.0e6*tm4/std::pow(gm, b4);
      for (int it=0; it<200; ++it) {
        const Real cm = std::sqrt(clo*chi);
        Real lsum = 0.0;
        for (int n=0; n<na*nq; ++n) lsum += std::fmax(cm*std::pow(ga[n], b4), fl4)*da[n];
        if (lsum > q4) chi = cm; else clo = cm;
      }
      cc = std::sqrt(clo*chi);
    }
    auto t4loc = [&](const Real g) {
      return (mirr_beta_ > 0.0) ? std::fmax(cc*std::pow(g, b4), fl4) : SQR(SQR(tdk));
    };
    // T^4-weighted lobe area and in-plane perimeter: the 3-D luminosity and the band's
    // per-length luminosity
    Real lw3 = 0.0, lw2 = 0.0;
    for (int n=0; n<na*nq; ++n) lw3 += t4loc(ga[n])*da[n];
    for (int ia=0; ia<na; ++ia) lw2 += t4loc(gl[ia])*dl[ia];
    // ---- v4 (user 10-16) problem/donor_irr_scale: the factor on the held donor light.
    // Default: 1 in 3-D; AUTO in the 2-D band, s = s_geo x s_m1:
    //  s_geo = (L_3D/(4 pi D^2)) / (L_2D/(2 pi D)) = int T^4 dA / (2 D oint T^4 dl): the
    //    band (reflecting theta walls = a donor infinitely extended in z) dilutes its light
    //    as 1/D (cylinder), the 3-D donor as 1/D^2 (sphere); for a sphere s_geo = R/D.
    //    D = a_sep - problem/donor_irr_rcmp (default 9.6, the gainer's facing photosphere)
    //  s_m1 = problem/donor_irr_m1fac: the band M1/vet_col transport's excess over the
    //    ideal cylindrical dilution at the facing point (CALIBRATED on Caltech, v4 static
    //    check; see accretor_rhd_design.md 15)
    const Real rcmp = pin->GetOrAddReal("problem", "donor_irr_rcmp", 9.6);
    const Real dcmp = asep - rcmp;
    const Real sgeo = lw3/(2.0*dcmp*lw2);
    const Real sm1 = pin->GetOrAddReal("problem", "donor_irr_m1fac", 1.0);
    if (pin->DoesParameterExist("problem", "donor_irr_scale")) {
      mirr_sc_ = pin->GetReal("problem", "donor_irr_scale");
    } else {
      mirr_sc_ = g3d_ ? 1.0 : sgeo*sm1;
    }
    const Real arc = pmbp->pradm1->arad/SQR(SQR(tunit));
    mirr_ec_ = mirr_sc_*arc*cc;
    mirr_ef_ = mirr_sc_*arc*fl4;
    mirr_e_ *= mirr_sc_;
    if (root) {
      std::printf("ry_per_accretor: donor_gdark beta %.4g: lobe area %.5g Rsun^2 (4 pi R_L^2 "
                  "%.5g), in-plane perimeter %.5g Rsun, T_mean %.1f K, floor %.3g T_mean; "
                  "C %.6g K^4/g^%.3g\n", mirr_beta_, atot, 4.0*M_PI*SQR(rlv), ltot,
                  std::pow(tm4, 0.25), tfr, cc, b4);
      const int ias[8] = {0, na/36, na/12, na/6, na/3, na/2, 2*na/3, na - 1};
      for (int n=0; n<8; ++n) {
        const int ia = ias[n];
        std::printf("ry_per_accretor: donor_gdark: in-plane alpha %6.1f deg from L1: g %.5g, "
                    "T_loc %.1f K\n", (ia + 0.5)*180.0/na, gl[ia], std::pow(t4loc(gl[ia]),
                    0.25));
      }
      {
        const Real t = DonorSurfT(rpp, msk_phl1_, msk_xl1_, 0.0, 0.0, 1.0, 0.0);
        const Real g = RocheGeffXYZ(rpp, rpp.asep, 0.0, t);
        std::printf("ry_per_accretor: donor_gdark: pole (z, R_pole %.4f): g %.5g, T_loc "
                    "%.1f K\n", t, g, std::pow(t4loc(g), 0.25));
      }
      std::printf("ry_per_accretor: donor_irr_scale %.6g (%s; s_geo %.6g = int T^4 dA / "
                  "(2 D oint T^4 dl), D = %.4f, s_m1 %.4g); 3-D target flux at D: "
                  "L_d/(4 pi D^2) = %.6g code\n", mirr_sc_,
                  pin->DoesParameterExist("problem", "donor_irr_scale") ? "key" :
                  (g3d_ ? "3-D default" : "2-D auto"), sgeo, dcmp, sm1,
                  0.25*pmbp->pradm1->c_light*arc*lw3/(4.0*M_PI*SQR(dcmp)));
    }
  }
  // (read, not added: the key stays out of the PAR_DUMP of runs that do not set it)
  szv_ = !pin->DoesParameterExist("problem", "sink_zero_v") ||
         pin->GetBoolean("problem", "sink_zero_v");

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
  // accretor-rhd-1015: the per-cycle state of the source/flux accumulators joins the
  // <time>/dt_retry checkpoint.  The boundary-flux registers acc0_/acc1_ and the stage
  // counter advance in every stage of a cycle, so a failed attempt would otherwise be
  // booked twice in Min..Ebnd/Erel (history only; the state is not affected).  The outer
  // face's advective radiation flux ibadv_x1max (S5 per-column bath), read by the
  // hydro-stage accumulator from the PREVIOUS cycle's solve like f0x1, is in the m1
  // minimal list (driver/dt_retry.cpp kDtrMinimal), not registered here (an alias entry
  // would be poisoned after its restore by dt_retry_poison).
  // History-output state (hist_prev_, hist_tprev_, rad_pn0_/ps0_) changes only at outputs,
  // which a failed attempt never reaches; rbath_, tcol_, rhoamb_, cref2_, gpl_ are set
  // at setup only.
  dtretry::RegisterValue("ry_stage_ctr", &stage_ctr_);
  erinj_ = 0.0;
  dtretry::RegisterValue("ry_erinj", &erinj_);
  for (int n=0; n<kNaccR; ++n) {
    dtretry::RegisterValue(("ry_acc0_" + std::to_string(n)).c_str(), &acc0_[n]);
    dtretry::RegisterValue(("ry_acc1_" + std::to_string(n)).c_str(), &acc1_[n]);
  }
  // bigdom-1015: the two extra registers (only with a nozzle or the mask)
  accx_on_ = stvol_ || msk_;
  for (int n=0; n<2; ++n) accx0_[n] = accx1_[n] = 0.0;
  if (accx_on_) {
    for (int n=0; n<2; ++n) {
      dtretry::RegisterValue(("ry_accx0_" + std::to_string(n)).c_str(), &accx0_[n]);
      dtretry::RegisterValue(("ry_accx1_" + std::to_string(n)).c_str(), &accx1_[n]);
    }
  }
  // acctrc-1019: the passive tracers (<hydro>/nscalars = 2; 0 = off, nothing registered)
  trc_ = false;
  for (int n=0; n<kNtr; ++n) trc0_[n] = trc1_[n] = 0.0;
  if (pmbp->phydro->nscalars != 0) {
    if (pmbp->phydro->nscalars != 2) {
      fatal("<hydro>/nscalars must be 0 or 2 (acctrc: s0 stream, s1 boundary inflow)");
    }
    if (!env_ || !gen_) fatal("the tracers (<hydro>/nscalars = 2) need inner = envelope and "
                              "thermo = general");
    if (pin->DoesParameterExist("hydro", "scalars_rst_missing") &&
        pin->GetString("hydro", "scalars_rst_missing").compare("zero") != 0) {
      fatal("<hydro>/scalars_rst_missing must be zero (or absent)");
    }
    trc_ = true;
    trc_n0_ = pmbp->phydro->nhydro;
    trc_file_ = pin->GetString("job", "basename") + ".trc.hst";
    for (int n=0; n<kNtr; ++n) {
      dtretry::RegisterValue(("ry_trc0_" + std::to_string(n)).c_str(), &trc0_[n]);
      dtretry::RegisterValue(("ry_trc1_" + std::to_string(n)).c_str(), &trc1_[n]);
    }
    if (root) {
      std::printf("ry_per_accretor: tracers on (<hydro>/nscalars = 2): s0 stream (nozzle), "
                  "s1 boundary inflow (r_in, r_out, theta walls); diagnostics in %s\n",
                  trc_file_.c_str());
    }
  }

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
    if (stl1_) {
      std::printf("ry_per_accretor: stream_inject = l1 (Lubow-Shu nozzle): L1 curvature A "
                  "%.4f B %.4f; c_T %.3f km/s along -x, T_L1 %.1f K, sigma_y %.4f Rsun "
                  "(sigma_phi %.5f rad at d_L1), peak rho %.6g code (%.4g g/cc) for "
                  "stream_mdot %.6g; nozzle cells r %.5f .. %.5f (outer face at d_L1 %.5f), "
                  "|phi| < %.4f rad; outer boundary: outflow (outer_rho %s), dark Marshak "
                  "top\n", l1a, l1b, l1_c_, l1_t_*eos.temp_cgs, sig_*str_r_, sig_, rho_s_,
                  rho_s_*(pmbp->punit ? pmbp->punit->density_cgs() : 1.0),
                  pin->GetReal("problem", "stream_mdot"), nz_rlo_, nz_rhi_, str_r_,
                  nsig_*sig_, orho.c_str());
    } else if (stvol_) {
      std::printf("ry_per_accretor: stream_inject = volume: the ballistic state above is at "
                  "stream_r %.4f Rsun = %.4f d_L1; nozzle cells r %.5f .. %.5f, |dphi| < "
                  "%.4f rad; peak rho %.6g code%s; back sink %s, floor_keep_momentum %s; "
                  "outer boundary: outflow (outer_rho %s), dark Marshak top\n",
                  str_r_, str_r_/(bo.dl1*asep), nz_rlo_, nz_rhi_, nsig_*sig_, rho_s_,
                  pin->DoesParameterExist("problem", "stream_mdot") ? " (from stream_mdot)" : "",
                  nzback_ ? "on" : "off", flmom_ ? "on" : "off", orho.c_str());
    }
    if (msk_) {
      std::printf("ry_per_accretor: donor_mask = roche: sink where Phi < Phi(L1) = %.6g "
                  "(km/s)^2 and x > d_L1 = %.4f Rsun; Mmsk/Emsk in hst columns 7/8\n",
                  msk_phl1_, msk_xl1_);
    }
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

  // acctrc-1019c: the cumulative history registers in the restart file (RyPerRstWrite)
  rst_pmbp_ = pmbp;
  picmax_rst_ = 0.0;
  if (!pin->DoesParameterExist("problem", "rst_registers") ||
      pin->GetBoolean("problem", "rst_registers")) {
    pgen_rst_write_func = RyPerRstWrite;
    if (restart) RyPerRstRestore(pgen_rststate);
  }
  user_bcs_func = RyPerBC;
  user_srcs_func = RyPerSrc;
  if (user_hist) user_hist_func = RyPerHist;
  pgen_final_func = RyPerFinal;
  if (env_) {
    EnvSetup(pin, pmbp, racc, rho_amb, restart);
    if (restart && pgen_rst_write_func != nullptr) RyPerRstBadv();
    user_bcs_func = RyPerBCEnv;
    user_srcs_func = RyPerSrcEnv;
    if (user_hist) user_hist_func = RyPerHistEnv;
    if (msk_) user_m1post_func = RyPerM1Post;
    if (thwc_) user_hydro_flux_func = RyPerWallFlux;
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
