//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file red_giant.cpp
//! \brief A self-luminous stellar envelope: the outer shell of a red giant (or any cool
//! giant) between an inner wall that injects the stellar luminosity and an outer wall
//! above the photosphere.  Point-mass gravity, stellar Rosseland opacities, and the
//! radiative-conduction operator carrying the flux; the initial column is the
//! radiative-convective hydrostatic equilibrium of that same operator.
//!
//! Runs on the cubed sphere and on spherical polar (x1 = r), and as a Cartesian x1
//! column for one-dimensional tests (x1 is then r - rin + mesh x1min).
//!
//! problem/rt_use_cons (default FALSE, here and in deep_hot_jupiter_rt.cpp; the FeCZ box
//! pgen box_convection defaults it TRUE).  Where the two-stream sweep reads the cell's
//! thermodynamic state from.  With it FALSE the sweep reads w0, the primitives the
//! PREVIOUS ConToPrim wrote, i.e. the state at the START of the stage: the RK update,
//! the explicit source terms and the implicit radial conduction have all moved the gas
//! since, so the radiative source is computed against a state that is one whole stage
//! stale.  That lag is not harmless -- it is a delayed thermostat, and a delayed
//! thermostat is a numerical FORCING whose amplitude is linear in dt.  It was measured
//! driving the organ-pipe acoustic modes of the closed FeCZ boxes (saturated v_rms
//! strictly monotone in dt and in nothing else; see the memory note
//! rt-source-dt-forcing.md), which is why box_convection hands the sweep the post-RK
//! state and why TRUE is the recommended setting for new runs here as well.
//!
//! It is left FALSE by default deliberately: flipping the default would silently change
//! every red-giant and deep-hot-Jupiter run ever restarted from these inputs, and the
//! two states differ by O(dt) in the radiative source.  Opt in from the input file.
//! (The conserved state is not guaranteed positive after an RK stage; the sweep's
//! EintFromCons path carries its own non-positive guard -- see utils/two_stream_rt.hpp.)
//!
//!   problem/mstar        mass inside the inner wall [g] (point mass; g = G M / r^2)
//!   problem/lstar        luminosity injected through the inner wall [erg/s]
//!   problem/teff         effective temperature [K]: the top of the column has
//!                        T^4 = Teff^4/2 (Eddington grey at tau = 0)
//!   problem/ptop         pressure at the outer wall [dyn/cm^2]
//!   problem/opac_table   log10 kappa_R on (log10 T, log10 rho): the file written by
//!                        tools/stellar_opac/merge_rosseland.py
//!   problem/mu           mean molecular weight for an IDEAL gas (ignored by a general
//!                        EOS, which carries its own composition)
//!   problem/vpert        velocity seed, in units of the local sound speed (0 = rest)
//!   problem/kappa_fac    scales kappa_R everywhere (the column here AND the conduction
//!                        operator, which is told).  A TEST knob: a giant envelope is
//!                        convective from tau ~ 10 down, so radiation carries only
//!                        grad_ad/grad_rad of L there and a 1-D column, which cannot
//!                        convect, cannot be steady.  kappa_fac ~ 1e-3 makes the
//!                        whole envelope radiative, and then the column must sit still.
//!   problem/kappa_const  > 0: a CONSTANT kappa_R [cm^2/g] replaces the table (column and
//!                        operator).  The grey Eddington atmosphere with constant
//!                        opacity has kappa p = g tau exactly, so grad_rad = tau/(4 tau +
//!                        8/3) < 0.25 < grad_ad: the whole ideal-gas column is RADIATIVE
//!                        and its equilibrium is known in closed form,
//!                        T^4 = (3/4) Teff^4 (tau + 2/3), p = g tau/kappa -- the 1-D
//!                        validation of the operator, the walls and the relaxation.
//!                        (kappa_fac cannot do this: scaling a T-dependent kappa leaves
//!                        kappa p / (g tau), and with it the convection zone, unchanged.)
//!   problem/mlt_alpha    > 0: a MIXING-LENGTH convective flux on the radial faces, with
//!                        l = alpha H_p (Boehm-Vitense; Kippenhahn et al. eq. 7.6),
//!                        wherever grad > grad_ad.  A giant envelope is convective from
//!                        just below the photosphere down, so a column that cannot
//!                        convect cannot be steady; this carries L - F_rad in 1-D, and
//!                        is a sub-grid option in 3-D (default 0 = off: resolve it).
//!   problem/rt_ck        true: the CORRELATED-K two-stream carries the optically thin
//!                        layers instead of the grey relaxation below.  The tables are
//!                        the Exo-FMS ones (problem/ck_table, ck_data_dir), the stellar
//!                        sweep is off (Teq = 0) and the internal temperature is teff,
//!                        so the band solver sees a self-luminous atmosphere.  The tau
//!                        blend hands over to radiative diffusion underneath exactly as
//!                        it does for the hot Jupiter.
//!   problem/inner_bc     wall (default) or open.  THE INNER WALL SITS INSIDE THE
//!                        CONVECTION ZONE -- for this star the radiative-convective
//!                        boundary is at 0.97 of the outer edge, so 94 % of the domain
//!                        convects and the wall is 12 density scale heights inside it.
//!                        A reflecting wall there bounces plumes, forces the convective
//!                        flux to zero and pins the envelope entropy at its initial
//!                        value.  `open` replaces it with the Stein-Nordlund / CO5BOLD
//!                        (Freytag+ 2012) treatment that solar_convection.cpp already
//!                        uses: the ghosts continue hydrostatically at the interior
//!                        temperature with the velocity copied, so plumes pass through,
//!                        and the lowest active layer is relaxed each step -- upflows
//!                        toward a prescribed deep adiabat (this is the energy input and
//!                        it sets the emergent flux), pressure toward the shell mean,
//!                        the mean density restored, and the net mass flux driven to
//!                        zero.  With `open` the luminosity is an OUTPUT, so
//!                        rad_flux_inner must be 0: the inflow entropy carries the
//!                        energy, not a diffusive flux through a wall.
//!   problem/outer_bc     wall (default) or open.  A wall fills its ghosts from the
//!                        INITIAL column, which stops being the interior's own state as
//!                        soon as the thin top adjusts; the standing jump that leaves at
//!                        the top face drains the top layer.  `open` continues the top
//!                        active cell hydrostatically at its own temperature, outflow
//!                        only unless problem/open_outer_noinflow is false.  Mass and
//!                        energy can then leave; problem/face_budget reports how much.
//!   problem/s_relax_cs, problem/s_relax_cp   the two relaxation rates (0.1, 0.3)
//!   problem/column_dump  if set, rank 0 writes the initial column to this file
//!   problem/ic_profile   if set, the initial stratification is READ from this file
//!                        ("r[cm] rho eint" ascending in r) instead of marched, and
//!                        problem/ptop and the outer-wall gas-solution check are unused
//!   problem/nfine        nodes of the fine radial grid the column lives on
//!   problem/user_srcs    must be true (the gravity source lives here)
//!
//! <hydro|mhd>/isotropic_conduction = radiative with rad_kappa_src = table_rho makes the
//! conduction operator read the same opacity table; rad_flux_inner < 0 lets this file
//! set the inner flux to L/(4 pi rin^2).  Use ix1_bc = ox1_bc = user: each radial
//! boundary is then either a reflecting wall whose ghosts come from the initial column,
//! or open (problem/inner_bc, problem/outer_bc).
//!
//! THE INITIAL COLUMN.  From the outer wall inward (and outward through the ghosts):
//!     d ln p / dr = -rho g / p,     d ln T / d ln p = min(grad_rad, grad_ad),
//!     grad_rad = 3 kappa_R p L / (16 pi a c G M T^4),
//! with rho(p, T) and grad_ad from the run's own EOS and kappa_R from the table.  The
//! diffusion form of grad_rad is exactly the Eddington grey slope, so the only thing the
//! optically thin top adds is the boundary value T^4(tau = 0) = Teff^4/2.  Where the
//! radiative gradient exceeds the adiabatic one the column follows the adiabat -- the
//! Schwarzschild criterion, mixing length omitted, which is what a resolved convection
//! zone is supposed to supply.  The result is the steady state of the conduction operator
//! in the radiative layers and of adiabatic convection below, so the run starts from
//! its own equilibrium rather than relaxing to it.

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/conduction.hpp"
#include "utils/rosseland.hpp"
#include "utils/wb_background.hpp"
#include "units/units.hpp"
#include "pgen.hpp"
#include "pgen_eos_utils.hpp"
#include "utils/correlated_k.hpp"
#include "utils/two_stream_rt.hpp"
#include "utils/runaway_scan.hpp"
#include "utils/eint_from_cons.hpp"

using pgen_eos::DensFromPT;
using pgen_eos::EintFromDensT;
using pgen_eos::PresFromEint;
using pgen_eos::GradAd;
using pgen_eos::TempKelvin;

void RedGiantGravity(Mesh *pm, Real bdt);
// problem/rt_strang, problem/rt_once_per_cycle: the radiation operator applied OUTSIDE
// the RK stages (ProblemGenerator::user_split_func).  See the declarations of rt_strang_
void RedGiantRTSplit(Mesh *pm, Real bdt);
void RedGiantGravityLedger(Mesh *pm, Real bdt);
void RedGiantBC(Mesh *pm);
void RedGiantFinal(ParameterInput *pin, Mesh *pm);

namespace {
// physical constants, cgs
constexpr Real kGrav = 6.674e-8;
constexpr Real kBoltz = 1.380649e-16;
constexpr Real kMH = 1.6726e-24;
constexpr Real kArad = 7.5657e-15;
constexpr Real kClight = 2.99792458e10;
constexpr Real kSigmaSB = 5.670374419e-5;

// the star and the grid, code units unless said otherwise
Real gm_ = 0.0;          // G M, code units
Real rin_ = 1.0;         // radius of the inner wall
Real x1min_ = 0.0;       // mesh x1min (the Cartesian column maps x1 -> r)
bool curv_ = false;      // spherical polar or cubed sphere: x1 IS r
Real rgas_ = 1.0;        // k/(mu m_H) for an ideal gas, code units
Real gm1_ = 0.4;
bool etotgrav_ = false;
// <problem>/wb_grav_source = background (default) | plain.  Under wellbalance_dynamic the
// radial gravity source is the BACKGROUND's own pressure difference across the cell,
// which cancels the well-balanced reconstruction exactly.  With `plain` the source
// reverts to -rho G M / r^2 (energy handled exactly as in the non-WB path) while the
// reconstruction and the fluxes keep the deviation form: the cancellation is deliberately
// broken, which separates an overstability living in the SOURCE (a delayed or weakened
// restoring force) from one living in the RECONSTRUCTION.  Not a production option.
bool wb_grav_plain_ = false;
// <problem>/wb_ramp [cm], 0 = a hard edge.  The STATIC well-balanced background is set to
// zero outside [<hydro>/wb_rmin, <hydro>/wb_rmax], which is exactly the plain scheme
// there -- a zero background makes the deviation the full state, the face background adds
// nothing, the flux removal removes nothing and the gravity source is -rho g.  But the
// deviation is reconstructed with a THREE-CELL stencil, so a hard step in the background
// puts a step of the background's own size into the slope of the cell beside it.  wb_ramp
// takes the background to zero over a raised cosine of this width instead; a few radial
// cells is enough.  The scheme stays exactly consistent for any background field -- the
// tapered one is simply less completely cancelled, ending at the plain scheme.
Real wb_ramp_ = 0.0;

// the initial column on a fine uniform grid in r: ln p [code] and T [K]
DvceArray1D<Real> lnp_d_, tk_d_;
// The MIXING-LENGTH convective velocity on the same fine column, for the velocity seed.
// problem/vpert_mlt (default true when the column is built with MLT): the seed amplitude
// is v_c(r) rather than a Mach number.  A red giant's surface convection is transonic --
// v_c ~ (F/rho)^(1/3) is 1.4e6 cm/s against c_s = 5.4e5 at the photosphere -- so a Mach
// 1e-3 seed starts three orders of magnitude below the amplitude the flow has to reach,
// and the layers cool faster than that many e-foldings take.  problem/vpert then
// multiplies v_c instead of c_s.
DvceArray1D<Real> vc_d_;
bool vpert_mlt_ = true;
// problem/ic_tau_rad: the optical depth above which the INITIAL COLUMN follows grad_rad
// rather than the convective gradient.  It used to be <hydro>/rad_tau_hi, which made the
// initial star depend on where the RT/diffusion handover was put: moving the handover to
// tau = 300 for stability also built the column radiative down to tau = 300, where
// grad_rad/grad_ad ~ 1e3, i.e. a different and far hotter star.  Defaults to rad_tau_hi,
// so nothing changes unless it is set.
Real ic_tau_rad_ = -1.0;
// problem/vpert_mach_max (default 0.3): cap the seed at this fraction of the local sound
// speed.  v_c peaks just under the photosphere, where convection is least efficient and
// the density is lowest, and there it is TRANSONIC -- seeding a single smooth mode at
// that amplitude drove the top 50 radial cells non-finite within 3400 cycles.  Deep down
// v_c/c_s ~ 1e-4 and the cap never binds, so this only tames the surface.
Real vpert_mach_max_ = 0.3;
// problem/vpert_rmin [cm], 0 = off: seed nothing below this radius, with a five-cell
// cosine ramp above it.  The seed is a STANDING radial flow, and where it crosses the
// radiative-convective boundary it advects entropy across an entropy STEP at a constant
// rate -- which shows up as a superadiabaticity at the two cells straddling the boundary
// growing LINEARLY from t = 0, at a rate no closure changes.  Measured: 2.5e-4 per 1e6 s
// at i = 5/6, identical to four digits between a run with the subgrid convective flux
// and one without.  Surface convection is what the seed is for; the radiative interior
// needs no seed at all.
Real vpert_rmin_ = 0.0;
// problem/vpert_cart (default false): seed the velocity from a CHART-INDEPENDENT function
// of the Cartesian direction instead of the meshblock-local (x2,x3) sinusoid.  The
// default seed is a sinusoid in MESHBLOCK-local coordinates, so it is the SAME field in
// all 96 blocks and is discontinuous at every panel seam: it imprints the chart on the
// flow, and the three chart copies of a cube vertex are then inequivalent by
// construction.  With
// vpert_cart the angular pattern is a cubic-harmonic polynomial of the unit direction,
// identical for every chart, so any chart signature that appears afterwards is the
// SCHEME's and not the initial condition's.
bool vpert_cart_ = false;
Real rlo_ = 0.0, drf_ = 1.0;
int nfine_ = 0;
// the opacity table (for the thin-layer relaxation) and the star's Teff
DvceArray2D<Real> ktab_;
DvceArray1D<Real> klT_, klD_;
int knT_ = 0, knD_ = 0;
Real teff_ = 0.0;
Real kfac_ = 1.0;        // problem/kappa_fac, applied to every table lookup
// The window in logR = log10 rho - 3 log10 T + 18 over which the opacity table is DATA.
// Both of its sources are tabulated in that variable over a finite range; outside it the
// merge fills with the edge value, which is constant in density and is not physics.  The
// table's header records the window and the initial column is checked against it below.
Real opac_lR_lo_ = -1.0e30, opac_lR_hi_ = 1.0e30;
// the open inner boundary (problem/inner_bc = open) and the deep adiabat it relaxes to
bool open_inner_ = false;
// problem/outer_bc = wall (default) or open.  A `wall` outer boundary fills its ghosts
// from the INITIAL column, so once the top layers have adjusted -- and they do, the
// optically thin top is exactly where the atmosphere is free to move -- the ghost holds
// a t = 0 density above a cell that no longer has it.  Dense over rarefied, at a face
// the Riemann solver then reads as a standing pressure jump: it pumps the top layer,
// which is the drain and the supersonic top cell seen in every run so far.  `open`
// continues the interior hydrostatically at ITS OWN temperature instead, so a settled
// atmosphere gives no flux through the face at all and the top is free to find its own
// stratification.
bool open_outer_ = false;
// problem/open_outer_noinflow (default true): with an open top, clamp the ghost's radial
// velocity to outflow.  The ghost is an extrapolation, not a reservoir; letting it push
// mass back down invents material the domain never had.  False copies the velocity as the
// inner boundary does, which is the right choice only if something above is being
// modelled.
bool open_outer_noinflow_ = true;
// problem/open_outer_lid (default false): keep the open ghost's hydrostatic (rho, e)
// continuation but MIRROR the radial velocity about the top face, so no mass crosses it.
// For a star whose radial pulsation lifts BOUND gas (a few per cent of v_esc) through a
// top that sits just above the photosphere: outflow-only drains it irreversibly.
bool open_outer_lid_ = false;
Real p_base_ = 0.0, t_base_ = 0.0;      // the initial column AT the inner wall, cgs
// --- AMBIENT MEDIUM outside the star (problem/bg_rho > 0).  A constant (rho, T) filling
// the domain above the point where the hydrostatic column thins to bg_rho, so that x1max
// can be pushed well beyond the photosphere without the star inflating to meet it.
//
// SETTING bg_rho AND bg_temp.  The whole region sits INSIDE the Bondi radius of the star
// (r_B = GM/c_s^2 > x1max for any T below ~3.9e5 K), so bg_rho is the density at infinity
// of a Bondi flow and the accretion rate is
//     Mdot = 4 pi lambda (GM)^2 rho_bg / c_bg^3,   lambda = 1/4 for gamma = 5/3,
// which RISES as the gas gets colder.  A local free-fall estimate, 4 pi r^2 rho v_ff at
// x1max, is ~5000x more permissive and is the wrong bound to use.  At 400 K, holding
// Mdot below 1e-10 Msun/yr wants bg_rho <~ 4e-22 g/cm^3; the mass constraint is
// irrelevant beside it (that background weighs ~2e-18 of the envelope).
Real bg_rho_ = 0.0;                     // problem/bg_rho [g/cm^3], 0 = no background
Real bg_temp_ = 400.0;                  // problem/bg_temp [K]
Real bg_rjoin_ = 0.0;                   // radius of the star/background join [code]
// problem/bg_hydrostatic: what sits above the star.  false (default) is the constant
// (bg_rho, bg_temp) medium below; true replaces it with a HYDROSTATIC isothermal corona
// at bg_temp, pressure-matched to the last stellar level.  A cold medium is not in
// balance with anything -- its scale height is ~1e-3 R -- so it free-falls onto the
// star's thin atmosphere, piles up, and its weight crushes the atmosphere onto the
// photosphere.  A hot corona at 6e5 K has H ~ 2r: it holds itself up, and since the
// join is at equal PRESSURE the density jumps DOWN by T_star/T_corona across it, so the
// light gas sits on top of the heavy gas and the interface is Rayleigh-Taylor stable.
bool bg_hydrostatic_ = false;           // problem/bg_hydrostatic
// problem/opac_floor [cm^2/g]: what the opacity becomes below the table's valid logR
// window.  NOT zero.  Zero makes the cell unable to emit, and then the semi-implicit
// radiative step has no equilibrium to relax to -- e_eq = e (A/E)^(1/4) diverges as
// E -> 0 and the update degenerates to an unbounded explicit one, which is exactly the
// state the rt_de_max limiter was papering over.  A small FINITE floor keeps emission and
// absorption consistent (both carry the same kappa), so the cell relaxes to its local
// radiative equilibrium temperature at a slow but well-defined rate.  It stays physically
// transparent: at 1e-5 cm^2/g and 1e-22 g/cm^3 the optical depth across the whole ambient
// shell is ~1e-15.
Real opac_floor_ = 1.0e-5;
// problem/opac_tmin [K]: if > 0, the table is CLAMPED in temperature -- every row below
// it is a copy of the first row at or above it.  The merged table is AESOPUS with grains,
// and kappa_R at 1e-11 g/cm^3 rises from 2e-5 cm^2/g at 2600 K to 1e-2 at 2000 K; this
// switches that condensation step off for a test (opac_tmin = 2500 leaves the molecular
// opacity above it untouched).  Default off.
Real opac_tmin_ = -1.0;
Real cs_change_ = 0.1, cp_change_ = 0.3;
// problem/open_budget (cycles, 0 = off): the open boundary is meant to be the star's
// ENERGY SOURCE, but each of its three passes rewrites the base shell's total energy and
// none of them is conservative by construction.  Accumulate what each one adds, in erg,
// and report it against L so a drain of the size of the luminosity cannot hide.
int open_budget_ = 0;
// problem/open_conserve (default true): make the density-restoring pass energy neutral.
// That pass adds the SAME drho to every cell of the base shell at the cell's own
// specific energy, so it moves mass between hot and cold columns and the total energy
// changes by drho * sum_i e_i V_i -- a sink of thousands of L once the shell has any
// entropy contrast (measured: -4e3 L, against the +1 L the boundary is meant to supply).
// With this on, whatever it adds is removed again uniformly over the shell.
bool open_conserve_ = true;
// problem/face_budget (cycles, 0 = off): the energy and mass the RADIAL FACES of the
// domain carry, integrated over each shell.  The open boundary's own passes are only
// half of its budget -- the other half is what the Riemann solver advects across the
// same face, and nothing else in the code reports it.  Sampled once per cycle, at the
// first stage, times dt: an O(dt) estimate of the cycle's transfer, which is all that
// is needed to compare a suspected sink against L.
// problem/open_nomassflux (default true with inner_bc = open): hold the NET mass flux
// through the inner FACE at zero.  The boundary's fourth pass drives the shell-mean
// radial velocity at the lowest cell CENTRE to zero, but the mass leaves across the
// face, and measured (problem/face_budget) that face drained 3.0e30 g and 2.0e44 erg --
// the whole of the run's energy loss -- in 2.3e4 s.  Each stage this removes the net
// mass the face just delivered, uniformly over the base shell, together with the energy
// and momentum that mass carried at the shell-mean specific enthalpy and velocity.  A
// zero-net-mass convective enthalpy flux, which is the luminosity the boundary is meant
// to supply, passes through untouched.
bool open_nomassflux_ = true;
// problem/wall_noflux (default true, inner_bc = wall): make the reflecting wall exactly
// impermeable.  The ghost carries the INITIAL COLUMN's density and pressure with only
// v_r mirrored -- which is what keeps it in balance with the well-balanced background --
// so once the interior drifts from the column the Riemann problem at the wall is no
// longer symmetric, HLLC's contact speed is not zero, and the face advects mass.
// Measured on the cubed sphere it reached -100 L and was still growing.  Each stage this
// removes the mass the face just delivered, cell by cell (a wall is impermeable
// locally, not on average), together with the energy that mass carried at the cell's own
// specific total enthalpy.  It does NOT touch the momentum flux, whose pressure term is
// what holds the star up, and it does NOT touch the conduction flux: the wall's imposed
// luminosity is added to the same flx1(IEN) and must survive.
bool wall_noflux_ = true;
// Kinetic energy density from the CONSERVED momenta.  On the cubed sphere u0(IM2,IM3)
// are COVARIANT components on a non-orthogonal basis whose angle has cos = cos_cell, so
// the kinetic energy is 0.5 m_i v^i = 0.5 [m1^2 + (m2^2 + m3^2 - 2 c m2 m3)/(1 - c^2)]/d,
// the same form history.cpp uses; sp and Cartesian have c = 0 and recover the plain sum.
// The sponge below used the orthogonal sum: subtracting the wrong kinetic energy and
// adding back fac^2 of it changes the TRUE internal energy by (KE_true - KE_orth)
// (1 - fac^2) every step, a systematic drain that is largest where the basis is most
// skewed -- the panel corners -- and grows with v^2.  Once the open-top atmosphere went
// transonic (~8 km/s at i >= 300, t ~ 1e6) the 24 cube-vertex cells inside the sponge
// were drained to the temperature floor within 2e4 s, became pressureless holes with
// 40 km/s in them, and the run NaN'd (T7/T8/T9, 2026-09-09).  prod11's vertex cells
// show the same drain at the 5-8 % level.
KOKKOS_INLINE_FUNCTION
Real CsKinetic(const Real d, const Real m1, const Real m2, const Real m3, const Real c) {
  return 0.5*(m1*m1 + (m2*m2 + m3*m3 - 2.0*c*m2*m3)/(1.0 - c*c))/d;
}
// problem/sponge (default true), sponge_zbot (0.96, as a fraction of the radial domain),
// sponge_c (0.1, in units of cs/dr): an absorbing layer under the outer wall.  Both
// radial walls are reflecting, so the acoustic flux convection drives up through the
// photosphere has nowhere to go; a wave amplitude grows as rho^-1/2 and this envelope is
// far more stratified than the solar box that needed the same layer.  The default zbot
// sits ABOVE tau = 2/3 (0.952 of the domain for the standard input) so the layer never
// damps the convection zone itself.  The kinetic energy removed is DISCARDED, not
// thermalised -- the layer stands in for an atmosphere that carries the flux away, and
// depositing it locally would keep inflating the top, which is the problem being fixed.
bool sponge_on_ = true;
Real sponge_zbot_ = 0.96, sponge_c_ = 0.1;
int open_debug_ = 0;          // problem/open_debug: boundary calls to trace, 0 = off
// THE OPEN-GHOST GUARD.  R9_lid_open32 died at t = 1.99e5 with the hydro dt going from
// 37 s to 7e-15 in ONE step and 1.64e6 non-finite cells spanning i ~ 200..323 on every
// rank within 26 cycles.  Nothing local spreads that far that fast; the two-stream's
// down-sweep from the top ghost does, and the ghost it starts from is built by
// fill_open_out, whose DensFromPT answers a failed inversion with the 1e6 g/cm^3
// sentinel (and a NaN temperature).  One bad ghost then poisons the whole column.
// These count the cells where the hydrostatic continuation was rejected and the initial
// column was used instead, and record the first one so the offending state is visible.
// POINTERS, not Views: a file-scope View is destroyed after Kokkos::finalize and aborts
// the run on the way out -- which is worse than useless, because it kills the final
// restart flush.  The two-stream's own scratch uses the same `new` and never-delete
// pattern for exactly this reason.
DvceArray1D<int> *bcguard_cnt_ = nullptr;   // [0] running total of fallbacks, atomic
DvceArray1D<Real> *bcguard_rec_ = nullptr;  // the first offender: m,k,j,i,d_i,e_i,t_i,p_g
int bcguard_seen_ = 0;              // total at the last readback
bool bcguard_reported_ = false;     // has the first offender been printed
// THE SPONGE'S POSITIVITY GUARD.  Same pointer-not-View reason as above.  [0] counts the
// cells the sponge refused to damp because their conserved energy was not yet physical;
// the record holds the first one: m,k,j,i,d,ei,ke,phi.
DvceArray1D<int> *spgguard_cnt_ = nullptr;
DvceArray1D<Real> *spgguard_rec_ = nullptr;
bool spgguard_reported_ = false;    // has the first offender been printed
// --- problem/nan_report (default false, so a default run is bit-identical).  After each
// operator in this file that writes u0(IEN), scan the active cells for a conserved state
// that is already unphysical -- non-finite d, momentum or energy, or E <= 0 -- and print
// the first one, then keep running.  The point is to say WHICH operator writes the NaN
// energy at a standing cold spike (the prod11_catch death), which the sponge guard can
// only report after the fact.  One reduction per operator, only when the switch is on.
// Pointers, not Views, for the same reason as the guards above.
bool nan_report_ = false;
// (problem/bc_use_cons is gone: the ghost fills ALWAYS continue the interior cell from
// u0.  Reading the previous stage's w0 made a restart non-reproducible -- w0 does not
// live in the restart file and is still zero at the first physical-BC call of
// Driver::Initialize -- so the old default is not a legal option.  See state_i in
// RedGiantBC.)
DvceArray1D<int> *nanrep_cnt_ = nullptr;
DvceArray1D<Real> *nanrep_rec_ = nullptr;
int nanrep_lines_ = 0;              // lines printed on this rank
const int nanrep_maxlines_ = 400;   // then stay quiet; the run is unaffected
int open_dbg_calls_ = 0;
int face_budget_ = 0;
int face_cycle_ = -1;
Real face_E_in_ = 0.0, face_E_out_ = 0.0, face_M_in_ = 0.0, face_M_out_ = 0.0;
Real face_E_in_l_ = 0.0, face_E_out_l_ = 0.0, face_t_last_ = 0.0;
Real open_dE_ent_ = 0.0, open_dE_prs_ = 0.0, open_dE_den_ = 0.0;
Real open_lstar_ = 0.0, open_t_last_ = 0.0;
Real open_dE_ent_l_ = 0.0, open_dE_prs_l_ = 0.0, open_dE_den_l_ = 0.0;
// --- OPERATOR-SPLIT SCHEDULING, ported from box_convection.cpp (:1922-2110).
// problem/rt_strang: the radiation operator (the two-stream column solve, its radiative
// force and, with rt_split_transverse, the implicit transverse diffusion) is applied as
// a Strang-split operator around the time integrator -- dt/2 before, dt/2 after -- and
// NOT inside the RK stages.  problem/rt_once_per_cycle is the Lie variant: once with the
// full dt after the last stage.  Both are enrolled through user_split_func.
bool rt_strang_ = false;
bool rt_once_ = false;
// problem/rt_split_transverse: with the column out of the stage, take the implicit
// transverse (ADI) operator out with it, over the same bdt.  Default ON under a split.
bool rt_split_tr_ = false;

// --- problem/vdamp_top_tau, problem/vdamp_top_time (0 = OFF, the default, and bitwise
// inert): the TAU-BASED top sponge ported from box_convection.cpp (:1691-1692 and the
// kernel at :2774).  Each stage the RADIAL momentum is relaxed as
//     m1 -> m1/(1 + f bdt/vdamp_top_time),
// with the ramp f = 1 where the COLUMN optical depth tau <= vdamp_top_tau/3, f = 0 where
// tau >= vdamp_top_tau and a raised cosine in log tau between (RadBlendWeight, the same
// shape the vertical blend uses).  It reads the per-column tau the conduction module
// builds each sweep (Conduction::rad_tauf), so it FOLLOWS the photosphere, which is what
// a radiation-pressure-dominated envelope needs; red_giant's own sponge (problem/sponge)
// is positioned by radius fraction instead and both may run together.
Real vdamp_tau_ = 0.0, vdamp_time_ = 20.0;
bool vdamp_printed_ = false;
// --- problem/vdamp_bot_cells, vdamp_bot_time, vdamp_bot_mean_only (0 = OFF, bitwise
// inert): the bottom sponge (box_convection.cpp:1708-1710, kernel :2822), over the
// lowest N ACTIVE RADIAL cells i = is .. is+N-1, with the ramp
//     f = (1 + cos(pi (i-is)/N))/2,   g = 1 - exp(-f bdt/vdamp_bot_time),
//     m1 -> m1 - g m1          (vdamp_bot_mean_only = false, the default), or
//     m1 -> m1 - g rho <v1>_shell   (vdamp_bot_mean_only = true; the shell mean is taken
//                                    over every angular cell at that radius, all ranks).
// The internal energy is conserved: the kinetic-energy change is written back to IEN.
// Insurance against the deep g-mode cavity under an open inner boundary.
// problem/vdamp_all_until (code time, 0 = off) and problem/vdamp_all_time: damp v1 in
// EVERY cell at the rate 1/vdamp_all_time while time < vdamp_all_until, then release.
// A relaxed-start test: does radial motion regrow from rest (instability) or not (IC)?
Real vda_until_ = 0.0, vda_time_ = 20.0;
int vdb_cells_ = 0;
Real vdb_time_ = 20.0;
bool vdb_mean_ = false;
DvceArray1D<Real> vdb_d_;        // (N): the shell sums of v1, then the shell means
HostArray1D<Real> vdb_h_;

// --- problem/rt_surface_dt, rt_surface_file: the EMERGENT TOP-FACE FLUX per angular
// cell, the spherical reading of box_convection.cpp's rt_surface dump (:889-955).  Every
// MeshBlock spans the whole radius, so the i = ie+1 face of every block IS the top of the
// star and every rank owns whole columns.  RECORD LAYOUT (little-endian, appended):
//     int64  n              number of angular cells in this record (all ranks, all panels
//                           concatenated in rank then block then (k,j) order)
//     double t              simulation time
//     double row[n][4]      panel, x2v, x3v, F_top   -- panel is the cubed-sphere panel
//                           index 0..5 as a double (-1 off the cubed sphere), x2v/x3v the
//                           cell-centre gnomonic coordinates on that panel in [-1,1], and
//                           F_top the net radial two-stream flux on the cell's top face.
// The box writes 3 doubles per row (x2, x3, F) with no panel column, so a reader such
// as analysis_0916/fmode_w7_t10 needs one change: read 4 columns and carry the panel.
Real surf_dt_ = 0.0;             // <= 0 disables it
Real surf_next_ = -1.0;          // next dump time; armed at the first source call
char surf_file_[256] = "rt_surface.bin";
DvceArray2D<Real> surf_d_;       // (ncol_local, 4): panel, x2, x3, F_top
HostArray2D<Real> surf_h_;
bool surf_alloc_ = false;
// --- problem/rt_profile_dt, rt_profile_file: SHELL-AVERAGED RADIAL PROFILES, the
// spherical reading of box_convection.cpp's rt_profile dump (:968-1050).  The average is
// AREA WEIGHTED over every angular cell at a radial index -- on the cubed sphere the
// gnomonic cells differ in area by ~4 between panel centre and corner, so a plain cell
// count would be a different average at every radius.  Every MeshBlock spans the whole
// radius, so the local index i - is IS the global radial index.  RECORD LAYOUT:
//     double t
//     int32  nx1, int32 nvar (= 8)
//     double x1v[nx1]                      the radial cell centres
//     double q[nvar][nx1]                  the shell means, slot order
//        0 rho   1 v1   2 rho v1   3 v1^2   4 v2^2+v3^2   5 T[K]   6 eint   7 v1(eint+p)
// with, PRECISELY (fixed 2026-09-17, tests_r9 task 3; every slot area weighted):
//        0  rho                      w0(IDN)                          [g/cm^3]
//        1  v1                       w0(IVX)                          [cm/s]
//        2  rho v1                                                    [g/cm^2/s]
//        3  v1^2                                                      [cm^2/s^2]
//        4  |v_perp|^2               v2^2+v3^2+2 cos(a) v2 v3, the GNOMONIC metric norm
//        5  T                        KELVIN: Hydro::wtemp (the CODE temperature the EOS
//                                    root find returned) times EOS_Data::temp_cgs
//        6  eint                     the INTERNAL ENERGY DENSITY [erg/cm^3], as
//                                    EintFromCons: u0(IEN) - KE(metric) - emag - rho Phi
//        7  v1 (eint + p)            with p the GAS+TAPERED-RADIATION PRESSURE, i.e.
//                                    wder(IDPR) under a general EOS
// WHAT WAS WRONG BEFORE.  Slot 5 was wtemp RAW, i.e. the code temperature, which on this
// star is 8.3145e7 times the kelvin the comment promises (T = 1.19e5 K at 0.94 R read out
// as 9.88e12).  Slot 6 was u0(IEN) - KE with NO rho*Phi subtraction, so under
// <hydro>/etotgrav = true (which every he4_presn input sets) it was eint + rho*GM(1/rin
// - 1/r): right at the inner wall, and 80x the true eint at the top of the domain, where
// eint is 1e-8 of the base and rho*Phi is not.  Slot 7 used w0(IPR) as `p`, but IPR ==
// IEN and under a GENERAL EOS the primitive slot holds the INTERNAL ENERGY, not the
// pressure (general_c2p_hyd.hpp:56 `w.e = u.e - e_k`; the pressure is in wder(IDPR)), so
// the enthalpy flux was v1*2*eint on the ideal branch's arithmetic.  Both are diagnostic
// only -- nothing reads them back into the evolution -- but tests_r6/an6.py's energy
// columns do, and were wrong in the envelope by those factors.
// Identical to the box's record except that x1v is the STRETCHED grid's own cell centres
// rather than a uniform reconstruction, so a reader must take x1v from the file (the box
// reader already does), and except for the two corrections above, which the box's own
// writer (box_convection.cpp:BoxConvProfileDump) does NOT have: it is plane-parallel with
// etotgrav off in every input that uses it, where slot 6 is already right, and its slot 7
// carries the same general-EOS defect.  Left untouched there on purpose: the box's
// productions are compared against their own archived records.
constexpr int kNProf = 8;        // see the slot list above
Real prof_dt_ = 0.0;             // <= 0 disables it
Real prof_next_ = -1.0;
char prof_file_[256] = "rt_profile.bin";
DvceArray2D<Real> prof_d_;       // (kNProf+1, nx1): the weighted shell SUMS, then means
HostArray2D<Real> prof_h_;
bool prof_alloc_ = false;

// --- problem/vpert_var: "v1" (default) seeds the RADIAL VELOCITY, exactly as before;
// "eint" instead scales the internal energy at fixed density, i.e. an ENTROPY seed.  The
// box's two options (box_convection.cpp:1383-1393), on red_giant's chart-free angular
// pattern.
int vpert_var_ = 0;
// --- problem/vpert_tau_lo, problem/vpert_tau_hi: the optical-depth WINDOW of the seed,
// the spherical reading of the box's height window vpert_zlo/vpert_zhi (which are quoted
// in the He inputs as the tau values they mean).  tau_hi is the DEEP edge, tau_lo the
// shallow one; the envelope is sin(pi ln(tau_hi/tau)/ln(tau_hi/tau_lo)) inside the
// window and 0 outside, so the seed vanishes at both edges.  tau comes from the INITIAL
// column, on its own fine radial grid.  Both 0 (the default) = no window, bitwise the
// old seed.
Real vpert_tau_lo_ = 0.0, vpert_tau_hi_ = 0.0;
DvceArray1D<Real> tau_d_;        // the initial column's tau on the fine radial grid

bool relax_ = false;     // the optically thin relaxation is on (radiative + tau blend)
bool rt_ck_ = false;     // problem/rt_ck: the band solver replaces that relaxation
// problem/rt_grey: the GREY two-stream on the same machinery, sharing the conduction
// module's opacity table.  Mutually exclusive with rt_ck; it replaces the same
// relaxation.  See two_stream_rt::rt_grey for why a grey control exists at all.
bool rt_grey_ = false;
// problem/ck_dump_t2, ck_dump_file2: re-arm the solver's one-shot column dump once the
// run reaches t2, so the emergent flux can be compared BEFORE and AFTER the atmosphere
// has adjusted to the band opacities.  The global energy budget cannot do this: the
// envelope holds ~1e48 erg while the flux imbalance integrates to ~1e40 over a short run.
Real ck_dump_t2_ = -1.0;
std::string ck_dump_file2_;
bool ck_dumped2_ = false;
Real mlt_alpha_ = 0.0;   // problem/mlt_alpha; > 0 turns the convective flux on
// problem/mlt_rmax [cm], 0 = off: no MLT flux above this radius.  The column is built
// RADIATIVE above tau = ic_tau_rad (10), where the two-stream / tau-blend carries L, so
// the subgrid flux must hand over there too.  Left unbounded it runs up into the thin
// layers, ends at whatever face the 10 %-of-e cap stops it, and dumps everything into
// that one cell: measured 43-140 % hotter photosphere and atmosphere within 6e5 s,
// Mach ~1.5 outflow at tau ~ 1 and a dt collapse (R1_mlt, 2026-09-09).
Real mlt_rmax_ = 0.0;
// problem/mlt_tau_lo, mlt_tau_hi: the SMOOTH hand-over to radiation, in the run's own
// evolving Rosseland tau.  A hard radial cut (mlt_rmax) does not work: at t = 0 F_mlt
// falls from 0.79 L to 0.03 L across the two cells at tau ~ 30, where the radiative
// solution still carries only 6e-4 L, so ~0.9 L is deposited in those two cells; the
// photosphere ran 45-150 % hot and the run NaN'd at the cut radius by t = 6.2e5
// (R1b_mlt, 2026-09-09).  Radiation only takes over as tau falls, so the subgrid flux
// is weighted by w = clamp(log10(tau/tau_lo)/log10(tau_hi/tau_lo), 0, 1): full MLT
// below tau_hi, nothing above tau_lo, and three decades of grid cells in between over
// which the two fluxes exchange the load.
Real mlt_tau_lo_ = 30.0;
Real mlt_tau_hi_ = 1.0e4;
Real lstar_cgs_ = 0.0;   // L [erg/s], for the physical cap F_mlt <= L/(4 pi r^2)
// problem/mlt_ramp_time [s], 0 = off: the applied MLT flux is multiplied by
// min(1, t/mlt_ramp_time).  Switching a flux of order L on in a single step drives an
// acoustic transient -- the R1d switch-on launched a Mach ~1.5 outflow through the
// atmosphere -- because the cut region is handed its whole heating rate before it can
// expand.  Ramping over a few sound-crossing times of the convection zone lets the
// stratification follow the source instead of being shocked by it.
Real mlt_ramp_ = 0.0;
// problem/mlt_ramp_down_time [s], 0 = off, and problem/mlt_hold_time [s]: the HAND-OVER
// to resolved convection.  mlt_ramp_time only ramps the subgrid flux UP, which is what a
// 1-D column that has to carry L for ever needs; a 3-D run wants the opposite -- the
// subgrid flux carries the convective luminosity from t = 0, while the seed grows, and
// is then withdrawn so that the resolved flow has to take the load.  The applied flux is
// multiplied by 1 for t <= mlt_hold_time and falls linearly to 0 at
// mlt_hold_time + mlt_ramp_down_time.  (The shell-mean closure already hands over
// PHYSICALLY -- it asks only for the deficit F_req - F_cond - F_2s - F_resolved, so it
// shrinks on its own as the resolved flux grows -- but that hand-over is a fixed point of
// the run, not a controlled one: if the resolved flow stalls, the subgrid flux simply
// keeps carrying L for ever and the run never tests whether convection has started.  The
// ramp-down forces the question at a time of the user's choosing.)
Real mlt_ramp_dn_ = 0.0;
Real mlt_hold_ = 0.0;
// problem/mlt_x_thr: the superadiabaticity above which the MLT AMPLITUDE is trusted.
// On the shell mean the deep convection zone sits on the adiabat to ~1e-6 in grad, and
// the SIGN of that difference is round-off: face by face it alternates, and the old
// closure (which returned zero wherever x <= 0) therefore left 14 % of L uncarried at
// the radiative-convective boundary, i = 7.  Above the threshold x is a real number
// (1e-2..0.5 near the surface) and F_MLT(x) binds; below it only the DEFICIT is known,
// and that is what gets carried.
Real mlt_x_thr_ = 1.0e-4;
// problem/mlt_flux_fix (default false), shell-mean path only: on faces the closure leaves
// alone (radiatively stable, or a negative deficit) carry the SIGNED shell residual
// F_req - F_rad - F_res - F_2s, capped at mlt_flux_fix_cap*F_req, below mlt_flux_fix_rmax
// [cm].  The initial column is in radiative equilibrium with ITS OWN transfer, not with
// the run's two-stream on this grid: the mismatch is 0.1-0.3 % of L (tests_r12), which
// on a Gamma ~ 0.95 envelope heats the layers under 0.93 R, cools those above and
// drives the secular expansion.  This makes L(r) = L on every interior face in the mean.
bool mlt_flux_fix_ = false;
// problem/mlt_split_deposit (needs problem/rt_strang): deposit div F_conv in the SAME
// Strang half steps as the two-stream, not inside the RK stages.  In the MLT zone
// div F_conv and div F_2s are each up to 6e-3 of eint per step and cancel; with one in
// the stages and the other split, the stages see a pressure error ~ dt div F_conv, a
// force dipole LINEAR in dt at the top of the MLT zone (tests_r12 fc3/fb0/fc07: +8.5 /
// +4.1 / +1.0 % of rho g at cfl 0.3 / 0.15 / 0.0375; sign flips with mlt_alpha = 0).
// The stage call still builds fconv_.  mlt_owed_ carries the half step made before the
// first build.
Real cs_max_ = 0.0, cs_mass_added_ = 0.0;   // problem/cs_max, see RgCsCeiling
bool rt_no_heat_ = false; // problem/rt_no_heat, see RedGiantRTSplit
int eledger_ = 0;         // problem/e_ledger, see RedGiantRTSplit
bool mlt_split_dep_ = false, mlt_fconv_ready_ = false;
Real mlt_owed_ = 0.0;
Real mlt_flux_fix_cap_ = 0.02, mlt_flux_fix_rmax_ = 0.0;
// problem/mlt_relax_time [s], 0 = off: relax the applied 1-D profile toward the freshly
// computed flux by dt/relax_time each call.  The deficit now contains the RESOLVED
// convective flux, which is a turbulent correlation and flickers from dump to dump; the
// applied subgrid flux should follow its mean, not its noise.  Seeded on the first call
// so that t = 0 already carries the full closure.
Real mlt_relax_ = 1.0e4;
bool mlt_relax_seeded_ = false;
std::string mlt_dump_ = "";  // problem/mlt_dump: write the faces of one column once
bool mlt_dumped_ = false;
// The dump has to wait for the two-stream.  On the FIRST source call of a run the RT
// solver has not run yet, so its face-flux array does not exist and F_2s reads zero --
// which made the t = 0 budget check show a full L unaccounted for above the photosphere
// where the two-stream is the only carrier.  Dump at the first call that has it (still
// the first cycle), or after a few calls if this run has no band solver at all.
int mlt_dump_calls_ = 0;
// problem/mlt_mean: evaluate the MLT closure ONCE PER RADIAL SHELL, on the horizontal
// mean state, instead of once per column.  Two measured defects of the per-column form
// force this.  (1) DEEP FLICKER: at i ~ 8-30 (tau ~ 1e12) the superadiabaticity that
// carries L is x ~ 1e-5, while the cell-to-cell noise of the discrete
// d ln T / d ln p across two cells is +-5e-6.  The 3/2 power and the L cap then turn
// that noise into a flux that flips between 0 and the full L on ADJACENT faces
// (drop/L = +-1 face to face in R1c_mlt's mltfaces_t0.txt): a checkerboard heating
// source of order L, which is not a subgrid model of anything.  (2) TOP DUMP: any
// shutoff of a per-column flux -- a radius cut or a tau taper -- deposits whatever it
// still carries into the last few cells: 0.9 L (R1b_mlt, photosphere +45-150 %, NaN)
// and 0.17 L per cell at tau ~ 50-150 (R1c_mlt).  The 1D profile is smooth (the noise
// averages down by sqrt(N_columns)), and the hand-over is done PHYSICALLY instead of by
// a taper: the shell only asks MLT for the part of L/(4 pi r^2) that radiative
// diffusion is not already carrying, so as the diffusive flux rises to L near the
// surface the convective one falls to zero on its own.
bool mlt_mean_ = true;
// (i, 13): shell sums of rho, e, T, p, the count, the CONDUCTION module's own discrete
// face flux and its blend weight, (rank 0's block 0, first column) that flux for one
// column so the shell mean can be checked against a real column, and the four sums the
// RESOLVED convective flux needs: rho v_r, rho v_r h, h and the kinetic term, and the
// TWO-STREAM's own net face flux
DvceArray2D<Real> shell_;
DvceArray1D<Real> fmlt1d_;  // (i): the shell's convective flux on face i, code units
DvceArray2D<Real> fmean_;   // (i, 17): the 1D closure's diagnostics, cgs
DvceArray2D<Real> fdiag_;   // (i, 8): grad, grad_ad, H_p, c_p, v, F, F_cap, F_used
DvceArray4D<Real> fconv_; // its radial face flux, code units, (m,k,j,i) on x1 faces
DvceArray4D<Real> taumlt_;  // cell-centred Rosseland tau from the top, (m,k,j,i)

//----------------------------------------------------------------------------------------
//! \brief the MLT closure's RELAXED STATE for the restart file (kPgenRstMagic in
//! pgen.hpp).  fmlt1d and mlt_relax_seeded_ are the only state the closure carries
//! between calls -- shell_, fmean_, fdiag_, taumlt_ and fconv_ are all rebuilt from the
//! current solution every call -- and without them a restart re-seeds fmlt1d to the
//! instantaneous target and the applied subgrid flux jumps by a factor 1.5-7 in one step,
//! which kicks the star (tests_r6/README.md section 3).  Layout: int32 nface, int32
//! seeded, then nface Reals, in the pgen's own code flux units.  Enrolled only when
//! problem/mlt_alpha > 0, so no other configuration writes a block at all.

std::vector<char> RedGiantMltRestartState() {
  std::vector<char> out;
  if (!(mlt_alpha_ > 0.0) || fmlt1d_.extent_int(0) <= 0) return out;
  const int nf = fmlt1d_.extent_int(0);
  auto hf = Kokkos::create_mirror_view(fmlt1d_);
  Kokkos::deep_copy(hf, fmlt1d_);
  const std::int32_t hdr[2] = {static_cast<std::int32_t>(nf),
                               static_cast<std::int32_t>(mlt_relax_seeded_ ? 1 : 0)};
  out.resize(sizeof(hdr) + nf*sizeof(Real));
  std::memcpy(out.data(), &(hdr[0]), sizeof(hdr));
  std::vector<Real> tmp(nf);
  for (int i = 0; i < nf; ++i) tmp[i] = hf(i);
  std::memcpy(out.data() + sizeof(hdr), tmp.data(), nf*sizeof(Real));
  return out;
}

//----------------------------------------------------------------------------------------
//! \brief the other half: put a restart's stored profile back into fmlt1d_.  Tolerates
//! its absence (a file written with problem/mlt_alpha = 0, or by any older binary) and
//! refuses a profile of the wrong length rather than reading past it -- a restart onto a
//! different radial grid has no business reusing a face-by-face profile.

void RedGiantMltRestartRestore(const std::vector<char> &blk) {
  if (!(mlt_alpha_ > 0.0) || fmlt1d_.extent_int(0) <= 0) return;
  const int nf = fmlt1d_.extent_int(0);
  if (blk.size() < 2*sizeof(std::int32_t)) {
    if (global_variable::my_rank == 0) {
      std::cout << "### red_giant: this restart file carries no MLT closure state; the "
                << "relaxed profile is re-seeded (expect a kick if "
                << "problem/mlt_relax_time > 0)" << std::endl;
    }
    return;
  }
  std::int32_t hdr[2];
  std::memcpy(&(hdr[0]), blk.data(), sizeof(hdr));
  if (hdr[0] != nf || blk.size() != sizeof(hdr) + nf*sizeof(Real)) {
    std::cout << "### FATAL ERROR in red_giant: the restart file's MLT closure profile "
              << "has " << hdr[0] << " faces and " << blk.size() << " bytes, this run "
              << "has " << nf << " faces.  Restarting a problem/mlt_alpha > 0 run onto a "
              << "different radial grid is not supported: set problem/mlt_relax_time = 0 "
              << "or start from scratch." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::vector<Real> tmp(nf);
  std::memcpy(tmp.data(), blk.data() + sizeof(hdr), nf*sizeof(Real));
  auto hf = Kokkos::create_mirror_view(fmlt1d_);
  for (int i = 0; i < nf; ++i) hf(i) = tmp[i];
  Kokkos::deep_copy(fmlt1d_, hf);
  mlt_relax_seeded_ = (hdr[1] != 0);
  if (global_variable::my_rank == 0) {
    std::cout << "### red_giant: MLT closure state read from the restart file ("
              << nf << " faces, seeded = " << mlt_relax_seeded_ << ")" << std::endl;
  }
}

KOKKOS_INLINE_FUNCTION Real GravAt(const Real gm, const Real r) {
  return gm/(r*r);
}
// increasing outward, zero at the inner wall (the code's sign convention for phi)
KOKKOS_INLINE_FUNCTION Real PotAt(const Real gm, const Real rin, const Real r) {
  return gm*(1.0/rin - 1.0/r);
}
KOKKOS_INLINE_FUNCTION Real RadiusOf(const bool curv, const Real x1, const Real rin,
                                     const Real x1min) {
  return curv ? x1 : (rin + x1 - x1min);
}
// The STATIC well-balanced background's radial window: 1 inside [rmin, rmax], 0 outside,
// with a raised-cosine ramp of width `ramp` on each side (ramp = 0: a hard edge).  A pure
// function of r, so cells and faces see the same field and the scheme stays consistent.
KOKKOS_INLINE_FUNCTION Real WbWindow(const Real r, const Real rmin, const Real rmax,
                                     const Real ramp) {
  Real w = 1.0;
  if (rmax > 0.0) {
    if (r >= rmax) return 0.0;
    if (ramp > 0.0 && r > rmax - ramp) {
      w *= 0.5*(1.0 - cos(M_PI*(rmax - r)/ramp));
    }
  }
  if (rmin > 0.0) {
    if (r <= rmin) return 0.0;
    if (ramp > 0.0 && r < rmin + ramp) {
      w *= 0.5*(1.0 - cos(M_PI*(r - rmin)/ramp));
    }
  }
  return w;
}

// the STATIC background state at radius r: the initial column's (rho, e), tapered by the
// window.  A free function, not a lambda captured inside the fill kernel: a device lambda
// may not capture another host-defined lambda.
KOKKOS_INLINE_FUNCTION
void WbBgAt(const EOS_Data &eos, const Real rgas, const Real igm1,
            const DvceArray1D<Real> &lnp, const DvceArray1D<Real> &tk, const int nf,
            const Real rlo, const Real drf, const Real rmin, const Real rmax,
            const Real ramp, const Real r, Real &dbg, Real &ebg);

// the column at radius r: ln p (code) and T (K), linear in r between fine nodes
KOKKOS_INLINE_FUNCTION void ColumnAt(const DvceArray1D<Real> &lnp,
                                     const DvceArray1D<Real> &tk, const int nfine,
                                     const Real rlo, const Real drf, const Real r,
                                     Real &lp, Real &t) {
  Real s = (r - rlo)/drf;
  int ii = static_cast<int>(s);
  ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
  const Real f = s - ii;
  lp = lnp(ii)*(1.0 - f) + lnp(ii+1)*f;
  t = tk(ii)*(1.0 - f) + tk(ii+1)*f;
}
KOKKOS_INLINE_FUNCTION
void WbBgAt(const EOS_Data &eos, const Real rgas, const Real igm1,
            const DvceArray1D<Real> &lnp, const DvceArray1D<Real> &tk, const int nf,
            const Real rlo, const Real drf, const Real rmin, const Real rmax,
            const Real ramp, const Real r, Real &dbg, Real &ebg) {
  Real lp, t;
  ColumnAt(lnp, tk, nf, rlo, drf, r, lp, t);
  const Real d = DensFromPT(eos, rgas, exp(lp), t);
  const Real w = WbWindow(r, rmin, rmax, ramp);
  dbg = w*d;
  ebg = w*EintFromDensT(eos, rgas, igm1, d, t);
}

// the same as ColumnAt, on host mirrors, for the start-up checks
template <typename V1>
void ColumnAtHost(const V1 &lnp, const V1 &tk, const int nfine, const Real rlo,
                  const Real drf, const Real r, Real &lp, Real &t) {
  Real s = (r - rlo)/drf;
  int ii = static_cast<int>(s);
  ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
  const Real f = s - ii;
  lp = lnp(ii)*(1.0 - f) + lnp(ii+1)*f;
  t = tk(ii)*(1.0 - f) + tk(ii+1)*f;
}

// log-bilinear kappa_R(T, rho) on the (log10 T, log10 rho) table, clamped
KOKKOS_INLINE_FUNCTION Real KappaTab(const DvceArray2D<Real> &tab,
                                     const DvceArray1D<Real> &lT,
                                     const DvceArray1D<Real> &lD, const int nT,
                                     const int nD, const Real tk, const Real rho) {
  const Real x = log10(tk), y = log10(rho);
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (x > lT(0)) {
    if (x >= lT(nT-1)) {
      i = nT-2; fx = 1.0;
    } else {
      i = static_cast<int>((x - lT(0))/(lT(1) - lT(0)));
      i = (i < 0) ? 0 : ((i > nT-2) ? nT-2 : i);
      fx = (x - lT(i))/(lT(i+1) - lT(i));
    }
  }
  if (y > lD(0)) {
    if (y >= lD(nD-1)) {
      j = nD-2; fy = 1.0;
    } else {
      j = static_cast<int>((y - lD(0))/(lD(1) - lD(0)));
      j = (j < 0) ? 0 : ((j > nD-2) ? nD-2 : j);
      fy = (y - lD(j))/(lD(j+1) - lD(j));
    }
  }
  const Real lk = (1.0-fx)*((1.0-fy)*tab(i,j) + fy*tab(i,j+1))
                +      fx *((1.0-fy)*tab(i+1,j) + fy*tab(i+1,j+1));
  return pow(10.0, lk);
}

//----------------------------------------------------------------------------------------
//! \fn ReadOpacityTable
//! \brief the merged Rosseland table: comment lines, one of them "# nT nD lTmin dlT
//! lDmin dlD", then nT*nD values of log10 kappa_R with T slowest

void ReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                      DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in red_giant: cannot open problem/opac_table '"
              << fname << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  Real lt0 = 0.0, dlt = 0.0, ld0 = 0.0, dld = 0.0;
  bool have_grid = false;
  std::vector<Real> vals;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      // "# valid_logR lo hi": outside that window the table is edge-filled, not data
      const std::size_t vp = line.find("valid_logR");
      if (vp != std::string::npos) {
        std::istringstream vs(line.substr(vp + 10));
        Real lo, hi;
        if (vs >> lo >> hi) { opac_lR_lo_ = lo; opac_lR_hi_ = hi; }
      }
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
    std::cout << "### FATAL ERROR in red_giant: opacity table '" << fname
              << "' has no grid line or " << vals.size() << " values for "
              << nT << " x " << nD << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // EXTEND THE DENSITY AXIS DOWN TO log10 rho = -24, AT THE FLOOR.  The table starts at
  // log10 rho = -14, and every consumer of it (the column here, the conduction operator,
  // the grey two-stream) CLAMPS to the table edge.  So the ambient medium at
  // rho = 2e-22 and the star's thin atmosphere at 1e-18 were both evaluated at
  // rho = 1e-14, where at 400-3000 K the node logR is -3.8..-6.6 -- inside the valid
  // window -- and the clamp handed back REAL molecular opacity, 1e-3..1e-2, for a cell
  // eight decades thinner than the node it was read from.  dtau per cell came out 1e2-1e3
  // too large and the background heated 400 -> 1800 K in 1000 s, when at opac_floor it
  // should take 1e5-1e6 s.  Prepending rows AT THE FLOOR puts the clamp edge below any
  // density the run can reach, so the lookup sees the true density and returns the floor.
  const Real lD_ext = -24.0;
  const int npre = (dld > 0.0 && ld0 > lD_ext)
                   ? static_cast<int>(std::lround((ld0 - lD_ext)/dld)) : 0;
  const int nD_raw = nD;
  const Real ld0_raw = ld0;
  nD += npre;
  ld0 -= npre*dld;
  Kokkos::realloc(tab, nT, nD);
  Kokkos::realloc(lT, nT);
  Kokkos::realloc(lD, nD);
  auto htab = Kokkos::create_mirror_view(tab);
  auto hlT = Kokkos::create_mirror_view(lT);
  auto hlD = Kokkos::create_mirror_view(lD);
  const Real lkfl = log10(opac_floor_);
  for (int i=0; i<nT; ++i) {
    hlT(i) = lt0 + i*dlt;
    for (int j=0; j<npre; ++j) htab(i,j) = lkfl;
    for (int j=0; j<nD_raw; ++j) htab(i,npre+j) = vals[i*nD_raw + j];
  }
  for (int j=0; j<nD; ++j) hlD(j) = ld0 + j*dld;
  // FLOOR THE OPACITY BELOW THE TABLE'S VALID logR WINDOW.  Both sources are tabulated in
  // logR = log rho - 3 log T + 18, and the merge EDGE-FILLED everything outside it: a
  // constant-in-density value that is not data.  For the THIN side (logR < lo) there is
  // essentially no material there to absorb, so replace the fabricated number with
  // problem/opac_floor -- small enough to be transparent, finite so that the cell can
  // still emit and the semi-implicit source has an equilibrium to relax to (see the note
  // on opac_floor_).  The ramp this puts on the last interpolation stencil inside the
  // window is the price, and it falls on layers whose tau is already negligible.  The
  // DENSE side (logR > hi) is left alone: there the edge-fill really would matter, and
  // the start-up check still refuses to run there.
  int nzero = 0;
  if (opac_lR_lo_ > -1.0e29) {
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nD; ++j) {
        // Floor ONLY the prepended rows (rho below the file's range: genuinely thin gas).
        // The file's own nodes below the logR window keep their EDGE-FILLED values: a
        // photospheric cell that has been shock-heated to 6e6 K at 2e-8 g/cm^3 also has
        // logR < -8, but it sits at tau ~ 100 in the diffusion operator, and handing it
        // the thin-gas floor made kappa_rad = 16 sigma T^3/(3 kappa rho) 1e5-1e6 x too
        // large and took the conduction timestep to 1e-7 s (R8_lid_open, 2026-09-09).
        // The edge value (electron scattering / free-free at the window's low-logR end)
        // is the physically sensible extrapolation for such a state.
        if (j < npre && (hlD(j) - 3.0*hlT(i) + 18.0) < opac_lR_lo_) {
          htab(i,j) = lkfl;
          ++nzero;
        }
      }
    }
  }
  // problem/opac_tmin: clamp the table in T (see opac_tmin_).  Done after the floor so
  // the prepended rows keep their floor value in every copied row too.
  int nclampT = 0;
  if (opac_tmin_ > 0.0) {
    const Real ltmin = log10(opac_tmin_);
    int iref = 0;
    while (iref < nT-1 && hlT(iref) < ltmin) ++iref;
    for (int i=0; i<iref; ++i) {
      for (int j=0; j<nD; ++j) htab(i,j) = htab(iref,j);
      ++nclampT;
    }
    std::cout << "red_giant: problem/opac_tmin = " << opac_tmin_ << " K: " << nclampT
              << " table rows below log10 T = " << hlT(iref) << " replaced by that row"
              << " (grain opacity OFF)" << std::endl;
  }
  Kokkos::deep_copy(tab, htab);
  Kokkos::deep_copy(lT, hlT);
  Kokkos::deep_copy(lD, hlD);
  std::cout << "red_giant: opacity table '" << fname << "', " << nT << " x " << nD
            << " nodes, log10 T " << lt0 << ".." << lt0 + (nT-1)*dlt
            << ", log10 rho " << ld0 << ".." << ld0 + (nD-1)*dld << " (" << npre
            << " rows prepended below the file's " << ld0_raw << " at the floor, so the"
            << " clamp cannot hand back real opacity for a thin cell); " << nzero
            << " prepended nodes (rho below the file) set to problem/opac_floor"
            << " = " << opac_floor_ << " cm^2/g (optically thin; the merge edge-filled"
            << " them)" << std::endl;
}
//----------------------------------------------------------------------------------------
//! \fn void RGNanScan
//! \brief problem/nan_report: scan u0 over the active cells for an unphysical conserved
//! state and print the first one, tagged with the operator that just ran.
//!
//! Called after each of this file's energy-updating operators; a no-op unless
//! problem/nan_report is set, so the default path is untouched.

void RGNanScan(Mesh *pm, const char *opname) {
  if (!nan_report_ || nanrep_lines_ >= nanrep_maxlines_) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &ccell_n = pmbp->pcoord->cos_cell;    // see CsKinetic
  const bool cs_n = pm->use_cubed_sphere;
  if (nanrep_cnt_ == nullptr) {
    nanrep_cnt_ = new DvceArray1D<int>("rg_nanrep_cnt", 1);
    nanrep_rec_ = new DvceArray1D<Real>("rg_nanrep_rec", 8);
  }
  auto ncnt = *nanrep_cnt_;
  auto nrec = *nanrep_rec_;
  Kokkos::deep_copy(ncnt, 0);
  Kokkos::deep_copy(nrec, 0.0);
  par_for("rg_nanscan", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = u0(m,IDN,k,j,i);
    const Real m1 = u0(m,IM1,k,j,i);
    const Real m2 = u0(m,IM2,k,j,i);
    const Real m3 = u0(m,IM3,k,j,i);
    const Real e = u0(m,IEN,k,j,i);
    const bool bad = !isfinite(d) || !isfinite(m1) || !isfinite(m2) || !isfinite(m3)
                     || !isfinite(e) || !(e > 0.0) || !(d > 0.0);
    if (!bad) return;
    if (Kokkos::atomic_fetch_add(&ncnt(0), 1) == 0) {
      const Real cc = cs_n ? ccell_n(m,k,j) : 0.0;
      nrec(0) = static_cast<Real>(m);
      nrec(1) = static_cast<Real>(k);
      nrec(2) = static_cast<Real>(j);
      nrec(3) = static_cast<Real>(i);
      nrec(4) = d;
      nrec(5) = e;
      nrec(6) = CsKinetic(d, m1, m2, m3, cc);
      nrec(7) = m1;
    }
  });
  auto hc = Kokkos::create_mirror_view(ncnt);
  Kokkos::deep_copy(hc, ncnt);
  if (hc(0) <= 0) return;
  auto hr = Kokkos::create_mirror_view(nrec);
  Kokkos::deep_copy(hr, nrec);
  ++nanrep_lines_;
  const int mb = static_cast<int>(hr(0));
  std::cout << "### rg nan_report [" << opname << "] rank " << global_variable::my_rank
            << " cycle " << pm->ncycle << " t = " << pm->time << ": " << hc(0)
            << " bad cell(s); first (m,k,j,i) = (" << mb << ","
            << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
            << static_cast<int>(hr(3)) << ") gid = " << (pmbp->gids + mb)
            << " d = " << hr(4) << " u(IEN) = " << hr(5) << " ke = " << hr(6)
            << " u(IM1) = " << hr(7) << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn RedGiantSurfaceDump
//! \brief append one record of the EMERGENT TOP-FACE FLUX per angular cell.
//!
//! The spherical reading of box_convection.cpp's rt_surface dump (:899-955): there the
//! columns are vertical and the record is F_top(x2, x3); here they are radial and the
//! record is F_top per angular cell of every cubed-sphere panel, all panels concatenated.
//! Every MeshBlock spans the whole radius (the mode-3 and mlt_mean checks enforce it), so
//! the i = ie+1 face of every block IS the top of the star and every rank owns whole
//! columns; each row therefore carries its OWN panel and (x2, x3) cell centre and the
//! reader never has to reconstruct the decomposition.  Record layout: see surf_file_.
//!
//! Called immediately after the two-stream, so what is written is the flux the solver has
//! just produced, on the state it was handed.  Diagnostic only.

void RedGiantSurfaceDump(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int ncol = nmb*nx3*nx2;
  if (!surf_alloc_ || surf_d_.extent_int(0) != ncol) {
    Kokkos::realloc(surf_d_, ncol, 4);
    Kokkos::realloc(surf_h_, ncol, 4);
    surf_alloc_ = true;
  }
  auto fb = two_stream_rt::rt_face_flux();
  const int nblk = two_stream_rt::rt_face_nblk();
  auto &size = pmbp->pmb->mb_size;
  auto &mbpanel = pmbp->pmb->mb_panel;
  const bool cs_s = pm->use_cubed_sphere;
  auto sd = surf_d_;
  par_for("rg_surf", DevExeSpace(), 0, ncol-1, KOKKOS_LAMBDA(const int idx) {
    const int m = idx/(nx3*nx2);
    const int kj = idx - m*(nx3*nx2);
    const int kk = kj/nx2;
    const int jj = kj - kk*nx2;
    Real ft = 0.0;
    for (int b=0; b<nblk; ++b) ft += fb(m,b,ie+1,ks+kk,js+jj);
    sd(idx,0) = cs_s ? static_cast<Real>(mbpanel.d_view(m)) : -1.0;
    sd(idx,1) = CellCenterX(jj, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    sd(idx,2) = CellCenterX(kk, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    sd(idx,3) = ft;
  });
  Kokkos::deep_copy(surf_h_, surf_d_);

  std::vector<double> mine(4*ncol);
  for (int n=0; n<ncol; ++n) {
    for (int c=0; c<4; ++c) mine[4*n+c] = static_cast<double>(surf_h_(n,c));
  }
  std::vector<double> all;
  int ntot = ncol;
#if MPI_PARALLEL_ENABLED
  const int nr = global_variable::nranks;
  std::vector<int> cnt(nr, 0), disp(nr, 0);
  int mycnt = 4*ncol;
  MPI_Gather(&mycnt, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  int tot = 0;
  for (int r=0; r<nr; ++r) {
    disp[r] = tot;
    tot += cnt[r];
  }
  all.resize((global_variable::my_rank == 0) ? tot : 1);
  MPI_Gatherv(mine.data(), mycnt, MPI_DOUBLE, all.data(), cnt.data(), disp.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ntot = tot/4;
#else
  all.swap(mine);
#endif
  if (global_variable::my_rank != 0) return;
  FILE *pf = std::fopen(surf_file_, "ab");
  if (pf == nullptr) {
    std::cout << "### FATAL ERROR in red_giant: cannot append to "
              << "problem/rt_surface_file '" << surf_file_ << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int64_t n64 = static_cast<int64_t>(ntot);
  const double tnow = static_cast<double>(pm->time);
  std::fwrite(&n64, sizeof(int64_t), 1, pf);
  std::fwrite(&tnow, sizeof(double), 1, pf);
  std::fwrite(all.data(), sizeof(double), 4*ntot, pf);
  std::fclose(pf);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn RedGiantSurfaceTick
//! \brief the cadence guard for the surface dump (problem/rt_surface_dt).  pm->time does
//! not move between RK stages, so advancing surf_next_ PAST it here is what makes the
//! dump fire once per cycle rather than once per stage.

void RedGiantSurfaceTick(Mesh *pm) {
  if (!(surf_dt_ > 0.0) || !two_stream_rt::rt_face_flux_ready()) return;
  if (surf_next_ < 0.0) surf_next_ = pm->time;
  if (pm->time >= surf_next_) {
    RedGiantSurfaceDump(pm);
    while (surf_next_ <= pm->time) surf_next_ += surf_dt_;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn RedGiantProfileDump
//! \brief append one record of the SHELL-AVERAGED radial profile.
//!
//! The spherical reading of box_convection.cpp's rt_profile dump (:981-1050).  The
//! average is AREA WEIGHTED, with the weight the cell's own radial face area: on the
//! cubed sphere the gnomonic cells differ in area by a factor ~4 between panel centre and
//! corner, so a plain cell count would be a different average at every radius.  One team
//! per radial index reduces over that shell's (m, k, j); the weighted SUMS are gathered
//! with MPI_Reduce and rank 0 divides by the summed weight (slot kNProf).
//!
//! Diagnostic only: it reads the primitives w0, the conserved energy u0(IEN) and the
//! temperature Hydro::wtemp that ConsToPrim left, and writes nothing back.  Record
//! layout: see prof_file_.

void RedGiantProfileDump(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  const int nkj = nmb*nx3*nx2;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (!prof_alloc_ || prof_d_.extent_int(1) != nx1) {
    Kokkos::realloc(prof_d_, kNProf+1, nx1);
    Kokkos::realloc(prof_h_, kNProf+1, nx1);
    prof_alloc_ = true;
  }
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto wt = is_mhd ? pmbp->pmhd->wtemp : pmbp->phydro->wtemp;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &ccellp = pmbp->pcoord->cos_cell;
  const bool cs_p = pm->use_cubed_sphere;
  // slots 6 and 7: the internal energy has to come out of u0 the way ConToPrim takes it
  // out -- with the gnomonic kinetic energy, the magnetic energy and, under
  // <hydro>/etotgrav, rho*Phi -- and the pressure out of wder, because under a general
  // EOS the primitive slot IPR holds the internal energy instead.  See the layout above.
  auto phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  auto bcc = is_mhd ? pmbp->pmhd->bcc0 : pmbp->phydro->w0;   // unread unless is_mhd
  auto wder_p = is_mhd ? pmbp->pmhd->wder : pmbp->phydro->wder;
  const bool etg_p = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool cgen_p = pmbp->phydro != nullptr ? pmbp->phydro->peos->eos_data.IsGeneral()
                                              : pmbp->pmhd->peos->eos_data.IsGeneral();
  const Real gm1_p = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                             : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  // slot 5 in KELVIN.  Hydro::wtemp is the CODE temperature the root find returned --
  // k_B T/(mu_ref m_u) in code energy per unit mass -- and EOS_Data::temp_cgs is the
  // kelvin per unit of it (eos.cpp:281, = temperature_cgs()/mu, so the <units>/mu cancels
  // and on a cgs <units> block it is m_u/k_B = 1.2027e-8).  Dumping wtemp raw made slot 5
  // 8.31e7 times the kelvin its comment promised.
  const Real tcgs_p = is_mhd ? pmbp->pmhd->peos->eos_data.temp_cgs
                             : pmbp->phydro->peos->eos_data.temp_cgs;
  auto pd = prof_d_;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nx1, Kokkos::AUTO);
  Kokkos::parallel_for("rg_prof", policy,
  KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tmember) {
    const int i = is + tmember.league_rank();
    array_sum::GlobalSum sum;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkj),
    [&](const int idx, array_sum::GlobalSum &ls) {
      const int m = idx/(nx3*nx2);
      const int kj = idx - m*(nx3*nx2);
      const int k = ks + kj/nx2;
      const int j = js + (kj - (kj/nx2)*nx2);
      const Real aw = area1(m,k,j,i);          // the area weight of this angular cell
      const Real d = w0(m,IDN,k,j,i);
      const Real v1 = w0(m,IVX,k,j,i);
      const Real v2 = w0(m,IVY,k,j,i);
      const Real v3 = w0(m,IVZ,k,j,i);
      const Real pg = cgen_p ? wder_p(m,IDPR,k,j,i) : w0(m,IPR,k,j,i)*gm1_p;
      // the kinetic energy on THIS grid's basis (the gnomonic tangents are not
      // orthogonal to each other), the same form CsKinetic uses
      const Real cc = cs_p ? ccellp(m,k,j) : 0.0;
      // the internal energy exactly as ConToPrim extracts it: EintFromCons subtracts the
      // metric kinetic energy, the magnetic energy and rho*Phi (the last only under
      // etotgrav, which is what slot 6 used to be wrong by)
      const Real emag = is_mhd ? MagEnergyCC(bcc, m, k, j, i) : 0.0;
      const Real ei = EintFromCons(u0, m, k, j, i, cc, cs_p, etg_p,
                                   etg_p ? phicc(m,k,j,i) : 0.0, emag);
      ls.the_array[0] += aw*d;
      ls.the_array[1] += aw*v1;
      ls.the_array[2] += aw*d*v1;
      ls.the_array[3] += aw*v1*v1;
      ls.the_array[4] += aw*(v2*v2 + v3*v3 + 2.0*cc*v2*v3);
      ls.the_array[5] += aw*(tcgs_p*wt(m,k,j,i));
      ls.the_array[6] += aw*ei;
      ls.the_array[7] += aw*v1*(ei + pg);
      ls.the_array[8] += aw;
    }, Kokkos::Sum<array_sum::GlobalSum>(sum));
    Kokkos::single(Kokkos::PerTeam(tmember), [&]() {
      for (int n=0; n<kNProf+1; ++n) pd(n, i-is) = sum.the_array[n];
    });
  });
  Kokkos::fence();
  Kokkos::deep_copy(prof_h_, prof_d_);

  const int nq = (kNProf+1)*nx1;
  std::vector<double> buf(nq);
  for (int n=0; n<kNProf+1; ++n) {
    for (int i=0; i<nx1; ++i) buf[n*nx1+i] = static_cast<double>(prof_h_(n,i));
  }
#if MPI_PARALLEL_ENABLED
  std::vector<double> rbuf((global_variable::my_rank == 0) ? nq : 1);
  MPI_Reduce(buf.data(), rbuf.data(), nq, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) buf.swap(rbuf);
#endif
  // the radial cell centres: the grid may be stretched, so they are read off the
  // coordinates rather than reconstructed from the mesh extent as the box does
  std::vector<double> x1v(nx1);
  {
    auto hx = Kokkos::create_mirror_view(pmbp->pcoord->x1v);
    Kokkos::deep_copy(hx, pmbp->pcoord->x1v);
    for (int i=0; i<nx1; ++i) x1v[i] = static_cast<double>(hx(0, is+i));
  }
  if (global_variable::my_rank != 0) return;
  for (int i=0; i<nx1; ++i) {
    const double wsum = buf[kNProf*nx1+i];
    const double fn = (wsum > 0.0) ? 1.0/wsum : 0.0;
    for (int n=0; n<kNProf; ++n) buf[n*nx1+i] *= fn;
  }
  FILE *pfp = std::fopen(prof_file_, "ab");
  if (pfp == nullptr) {
    std::cout << "### FATAL ERROR in red_giant: cannot append to "
              << "problem/rt_profile_file '" << prof_file_ << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const double tnow = static_cast<double>(pm->time);
  const int32_t n1 = static_cast<int32_t>(nx1);
  const int32_t nv = static_cast<int32_t>(kNProf);
  std::fwrite(&tnow, sizeof(double), 1, pfp);
  std::fwrite(&n1, sizeof(int32_t), 1, pfp);
  std::fwrite(&nv, sizeof(int32_t), 1, pfp);
  std::fwrite(x1v.data(), sizeof(double), nx1, pfp);
  std::fwrite(buf.data(), sizeof(double), kNProf*nx1, pfp);
  std::fclose(pfp);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn RedGiantProfileTick
//! \brief the cadence guard for the profile dump (problem/rt_profile_dt).  Called from
//! the top of the source term, so the record is the start-of-cycle state; pm->time does
//! not move between stages, so advancing prof_next_ past it fires this once per cycle.

void RedGiantProfileTick(Mesh *pm) {
  if (!(prof_dt_ > 0.0)) return;
  if (prof_next_ < 0.0) prof_next_ = pm->time;
  if (pm->time >= prof_next_) {
    RedGiantProfileDump(pm);
    while (prof_next_ <= pm->time) prof_next_ += prof_dt_;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantRTSweep
//! \brief the band/grey two-stream call itself, shared by the in-stage source term and
//! the operator-split entry point (problem/rt_strang, rt_once_per_cycle) so that the two
//! are the SAME code with a different dt.

void RedGiantRTSweep(Mesh *pm, Real bdt) {
  if (!ck_dumped2_ && ck_dump_t2_ >= 0.0 && !ck_dump_file2_.empty() &&
      pm->time >= ck_dump_t2_) {
    ck_dumped2_ = true;
    two_stream_rt::rt_dump_file = ck_dump_file2_;
    two_stream_rt::rt_dump_done = false;      // re-arm the one-shot dump
  }
  two_stream_rt::picket_fence_two_stream_RT(pm, bdt);
  RGNanScan(pm, "RT_two_stream");
  runaway_scan::Scan(pm, "RT_two_stream");
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  eledger_ = pin->GetOrAddInteger("problem", "e_ledger", 0);
  rt_no_heat_ = pin->GetOrAddBoolean("problem", "rt_no_heat", false);
  cs_max_ = pin->GetOrAddReal("problem", "cs_max", 0.0);
  user_srcs_func = (eledger_ > 0 || cs_max_ > 0.0) ? RedGiantGravityLedger
                                                   : RedGiantGravity;
  user_bcs_func = RedGiantBC;
  pgen_final_func = RedGiantFinal;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  const Real gamma = eos.gamma;
  const Real gm1 = gamma - 1.0, igm1 = 1.0/gm1;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;
  if (pmbp->punit == nullptr) {
    std::cout << "### FATAL ERROR in red_giant: a <units> block is required" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real lunit = pmbp->punit->length_cgs();
  const Real dunit = pmbp->punit->density_cgs();
  const Real punit_ = pmbp->punit->pressure_cgs();
  const Real vunit = pmbp->punit->velocity_cgs();
  const bool curv = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
  if (!curv && !(pmy_mesh_->one_d)) {
    std::cout << "### FATAL ERROR in red_giant: on a Cartesian mesh only a 1D x1 column "
              << "is supported (or use spherical polar / the cubed sphere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // --- the star
  const Real mstar = pin->GetReal("problem", "mstar");
  const Real lstar = pin->GetReal("problem", "lstar");
  lstar_cgs_ = lstar;
  const Real teff = pin->GetReal("problem", "teff");
  const Real ptop_cgs = pin->GetReal("problem", "ptop");
  const Real mu = pin->GetOrAddReal("problem", "mu", 0.62);
  const Real vpert = pin->GetOrAddReal("problem", "vpert", 0.0);
  vpert_mlt_ = pin->GetOrAddBoolean("problem", "vpert_mlt", true);
  ic_tau_rad_ = pin->GetOrAddReal("problem", "ic_tau_rad", -1.0);
  bg_rho_ = pin->GetOrAddReal("problem", "bg_rho", 0.0);
  bg_temp_ = pin->GetOrAddReal("problem", "bg_temp", 400.0);
  bg_hydrostatic_ = pin->GetOrAddBoolean("problem", "bg_hydrostatic", false);
  // problem/bg_rtop: the radius at which the column is anchored at problem/ptop.  The
  // column is integrated DOWNWARD from here, so this -- not x1max -- is what fixes where
  // the star's photosphere lands.  Default x1max reproduces the old behaviour exactly;
  // with a background, set it to the top of the STAR and let x1max sit far above it.
  const Real bg_rtop = pin->GetOrAddReal("problem", "bg_rtop", 0.0);
  vpert_mach_max_ = pin->GetOrAddReal("problem", "vpert_mach_max", 0.3);
  vpert_rmin_ = pin->GetOrAddReal("problem", "vpert_rmin", 0.0);
  vpert_cart_ = pin->GetOrAddBoolean("problem", "vpert_cart", false);
  // problem/vpert_var, problem/vpert_tau_lo/hi: WHICH variable the seed perturbs and the
  // optical-depth window it lives in (see their declarations).  Both are the box's
  // options on red_giant's chart-free angular pattern; the defaults are the old seed.
  {
    const std::string pvs = pin->GetOrAddString("problem", "vpert_var", "v1");
    if (pvs.compare("v1") == 0) {
      vpert_var_ = 0;
    } else if (pvs.compare("eint") == 0) {
      vpert_var_ = 1;
    } else {
      std::cout << "### FATAL ERROR in red_giant: problem/vpert_var must be \"v1\" or "
                << "\"eint\"" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  vpert_tau_lo_ = pin->GetOrAddReal("problem", "vpert_tau_lo", 0.0);
  vpert_tau_hi_ = pin->GetOrAddReal("problem", "vpert_tau_hi", 0.0);
  if ((vpert_tau_lo_ > 0.0) != (vpert_tau_hi_ > 0.0) ||
      (vpert_tau_hi_ > 0.0 && !(vpert_tau_hi_ > vpert_tau_lo_))) {
    std::cout << "### FATAL ERROR in red_giant: the seed window needs both "
              << "problem/vpert_tau_lo and problem/vpert_tau_hi > 0 with tau_hi (the "
              << "DEEP edge) > tau_lo" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real kfac = pin->GetOrAddReal("problem", "kappa_fac", 1.0);
  const Real kconst = pin->GetOrAddReal("problem", "kappa_const", 0.0);
  mlt_alpha_ = pin->GetOrAddReal("problem", "mlt_alpha", 0.0);
  mlt_rmax_ = pin->GetOrAddReal("problem", "mlt_rmax", 0.0);
  mlt_tau_lo_ = pin->GetOrAddReal("problem", "mlt_tau_lo", 30.0);
  mlt_tau_hi_ = pin->GetOrAddReal("problem", "mlt_tau_hi", 1.0e4);
  // problem/mlt_alpha_ic: the mixing length used to build the INITIAL COLUMN, which is
  // a separate question from whether the MLT flux runs as a source term.  Defaults to
  // mlt_alpha, so setting one knob does both; set it alone (with mlt_alpha = 0) to start
  // a RESOLVED-convection run from a stratification that already carries L.
  const Real mlt_alpha_ic = pin->GetOrAddReal("problem", "mlt_alpha_ic", mlt_alpha_);
  mlt_dump_ = pin->GetOrAddString("problem", "mlt_dump", "");
  mlt_mean_ = pin->GetOrAddBoolean("problem", "mlt_mean", true);
  mlt_ramp_ = pin->GetOrAddReal("problem", "mlt_ramp_time", 0.0);
  mlt_ramp_dn_ = pin->GetOrAddReal("problem", "mlt_ramp_down_time", 0.0);
  mlt_hold_ = pin->GetOrAddReal("problem", "mlt_hold_time", 0.0);
  mlt_x_thr_ = pin->GetOrAddReal("problem", "mlt_x_thr", 1.0e-4);
  mlt_flux_fix_ = pin->GetOrAddBoolean("problem", "mlt_flux_fix", false);
  mlt_split_dep_ = pin->GetOrAddBoolean("problem", "mlt_split_deposit", false);
  mlt_flux_fix_cap_ = pin->GetOrAddReal("problem", "mlt_flux_fix_cap", 0.02);
  mlt_flux_fix_rmax_ = pin->GetOrAddReal("problem", "mlt_flux_fix_rmax", 0.0);
  mlt_relax_ = pin->GetOrAddReal("problem", "mlt_relax_time", 1.0e4);
  // the shell mean indexes cells by their GLOBAL radial position, so a MeshBlock that
  // covers only part of the radius would average cells at different radii together
  if (mlt_alpha_ > 0.0 && mlt_mean_ &&
      pmy_mesh_->mesh_indcs.nx1 != pmy_mesh_->mb_indcs.nx1) {
    std::cout << "### FATAL ERROR in red_giant: problem/mlt_mean needs every MeshBlock "
              << "to span the whole radial range, but mesh/nx1 = "
              << pmy_mesh_->mesh_indcs.nx1 << " and meshblock/nx1 = "
              << pmy_mesh_->mb_indcs.nx1 << ". Set meshblock/nx1 = mesh/nx1, or set "
              << "problem/mlt_mean = false." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  {
    std::string wbgs = pin->GetOrAddString("problem", "wb_grav_source", "background");
    if (wbgs.compare("background") == 0) {
      wb_grav_plain_ = false;
    } else if (wbgs.compare("plain") == 0) {
      wb_grav_plain_ = true;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<problem> wb_grav_source = '" << wbgs
                << "' must be 'background' or 'plain'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // the ramp that takes the STATIC background to zero at wb_rmin / wb_rmax [cm]
  wb_ramp_ = pin->GetOrAddReal("problem", "wb_ramp", 0.0)/lunit;
  sponge_on_ = pin->GetOrAddBoolean("problem", "sponge", true);
  nan_report_ = pin->GetOrAddBoolean("problem", "nan_report", false);
  runaway_scan::on = pin->GetOrAddBoolean("problem","runaway_scan",false);
  runaway_scan::rmin = pin->GetOrAddReal("problem","runaway_rmin",3.3e12);
  runaway_scan::ratio_print = pin->GetOrAddReal("problem","runaway_ratio",3.0);
  runaway_scan::ratio_state = pin->GetOrAddReal("problem","runaway_ratio_state",10.0);
  two_stream_rt::rt_nan_report = nan_report_;
  sponge_zbot_ = pin->GetOrAddReal("problem", "sponge_zbot", 0.96);
  sponge_c_ = pin->GetOrAddReal("problem", "sponge_c", 0.1);
  opac_floor_ = pin->GetOrAddReal("problem", "opac_floor", 1.0e-5);
  opac_tmin_ = pin->GetOrAddReal("problem", "opac_tmin", -1.0);
  const std::string opac = pin->GetString("problem", "opac_table");
  const std::string dump = pin->GetOrAddString("problem", "column_dump", "");
  // problem/ic_profile: an EXTERNALLY SUPPLIED 1-D stratification, three columns
  // "r[cm] rho[g/cm^3] eint[erg/cm^3]" in ASCENDING r (the same file format
  // box_convection reads, with r in place of z).  When it is set the pgen's own
  // hydrostatic+MLT march is SKIPPED ENTIRELY and the column below is built from the
  // file instead.  It exists because that march cannot build a RADIATION-DOMINATED
  // envelope: it anchors the top at the grey tau = 0 temperature, where a T^4/3 is
  // within a few per cent of the total pressure, refuses that state outright (the "no
  // gas solution at the outer wall" abort), knows nothing about the EOS radiation taper
  // (<hydro>/eos_rad_rho_hi/lo) or about the two-stream radiative force
  // (problem/rt_rad_force) that carries the pressure the taper removed, and the error
  // those make at the top runs away inward.  A file built OFFLINE against this run's own
  // table -- solving p_gas(rho,T) + a T^4/3 = p_column(r) node by node and then applying
  // the taper at fixed (rho,T) -- has none of those problems.  rho may be NON-MONOTONE
  // (a density inversion is physical in a super-Eddington envelope), which is why the
  // interpolation below searches on r and never on rho.
  const std::string icprof = pin->GetOrAddString("problem", "ic_profile", "");
  x1min_ = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  rin_ = curv ? x1min_ : pin->GetReal("problem", "rin")/lunit;
  const Real rout = curv ? x1max : rin_ + (x1max - x1min_);
  gm_ = kGrav*mstar/(lunit*lunit*lunit)*(pmbp->punit->time_cgs()
                                         *pmbp->punit->time_cgs());
  // ideal gas only: p = rho Rgas T with T in KELVIN (pgen_eos_utils convention), so
  // Rgas carries velocity^2 per kelvin in code units
  rgas_ = kBoltz/(mu*kMH)/(vunit*vunit);
  curv_ = curv; gm1_ = gm1; etotgrav_ = etotgrav;
  const Real gm = gm_, rin = rin_, x1min = x1min_, rgas = rgas_;

  // --- the opacity table, for the column here and for the conduction operator
  DvceArray2D<Real> ktab;
  DvceArray1D<Real> klT, klD;
  int knT = 0, knD = 0;
  ReadOpacityTable(opac, ktab, klT, klD, knT, knD);
  if (kconst > 0.0) {
    auto h = Kokkos::create_mirror_view(ktab);
    for (int i=0; i<knT; ++i) {
      for (int j=0; j<knD; ++j) h(i,j) = log10(kconst);
    }
    Kokkos::deep_copy(ktab, h);
    std::cout << "red_giant: kappa_const = " << kconst << " cm^2/g REPLACES the table"
              << std::endl;
  }
  Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
  if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
    if (pc->rad_kappa_tab && !pc->rad_kappa_rho) {
      std::cout << "### FATAL ERROR in red_giant: use rad_kappa_src = table_rho, the "
                << "stellar table is on (T, rho)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    pc->rad_kappa_fac = kfac;
    if (pc->rad_kappa_rho) {
      pc->rad_kr_tab = ktab;
      pc->rad_kr_lT = klT;
      pc->rad_kr_lP = klD;
      pc->rad_kr_nT = knT;
      pc->rad_kr_nP = knD;
    }
    if (pc->rad_flux_inner < 0.0) {
      const Real fin = lstar/(4.0*M_PI*SQR(rin*lunit));       // erg/cm^2/s
      pc->rad_flux_inner = fin/(punit_*vunit);
      std::cout << "red_giant: inner flux L/(4 pi rin^2) = " << fin
                << " erg/cm^2/s = sigma (" << pow(fin/kSigmaSB, 0.25) << " K)^4"
                << std::endl;
    }
    // THE OPTICALLY THIN LAYERS.  Diffusion is the wrong physics above tau ~ 1 and,
    // worse, explicitly stiff there: kappa_rad = 16 sigma T^3/(3 kappa_R rho) diverges
    // as kappa_R rho -> 0 while the layer's heat capacity vanishes, so the explicit
    // diffusion step collapses to microseconds even with the flux limiter (which does
    // nothing where the gradient is flat).  So the tau blend is REQUIRED: diffusion
    // carries the flux with weight w(tau) and, with weight 1 - w, each thin cell relaxes
    // toward the grey Eddington temperature of its own optical depth,
    //     T_eq^4 = (3/4) Teff^4 (tau + 2/3),
    // on its radiative time  t_rad = (rho c_v) / (4 kappa_P rho sigma T^3), applied as
    // the exact exponential decay so it is stable at any timestep.  In radiative
    // equilibrium that profile IS what the diffusion delivers below, so the two hand
    // over consistently for the luminosity Teff encodes.  See RedGiantGravity.
    if (!pc->rad_tau_mode) {
      std::cout << "### FATAL ERROR in red_giant: radiative conduction needs the tau "
                << "blend here (set rad_tau_lo/rad_tau_hi, e.g. 1 and 10): the optically "
                << "thin layers are relaxed toward the Eddington profile, not diffused"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    relax_ = true;
  }

  // --- the TAU-BASED TOP SPONGE (problem/vdamp_top_tau, vdamp_top_time), ported from
  // box_convection.cpp:1691-1706.  It reads the per-column optical depth the conduction
  // module builds, so the tau blend has to be on.  See the declaration of vdamp_tau_.
  vdamp_tau_ = pin->GetOrAddReal("problem", "vdamp_top_tau", 0.0);
  vdamp_time_ = pin->GetOrAddReal("problem", "vdamp_top_time", 20.0);
  if (vdamp_tau_ > 0.0) {
    if (pc == nullptr || !pc->rad_tau_mode) {
      std::cout << "### FATAL ERROR in red_giant: problem/vdamp_top_tau needs the "
                << "per-column optical depth, i.e. the tau blend (<hydro>/rad_tau_hi > 0)"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(vdamp_time_ > 0.0)) {
      std::cout << "### FATAL ERROR in red_giant: problem/vdamp_top_time must be > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // --- the BOTTOM SPONGE (problem/vdamp_bot_cells, vdamp_bot_time,
  // vdamp_bot_mean_only), ported from box_convection.cpp:1707-1733.  The box additionally
  // demands a CLOSED bottom wall; here the sponge exists precisely to tame the deep
  // g-mode cavity under an OPEN inner boundary, so that check is dropped -- but the
  // shell-mean variant needs every MeshBlock to span the whole radius, as the MLT shell
  // mean does, because it indexes cells by their global radial position.
  vda_until_ = pin->GetOrAddReal("problem", "vdamp_all_until", 0.0);
  vda_time_ = pin->GetOrAddReal("problem", "vdamp_all_time", 20.0);
  vdb_cells_ = pin->GetOrAddInteger("problem", "vdamp_bot_cells", 0);
  vdb_time_ = pin->GetOrAddReal("problem", "vdamp_bot_time", 20.0);
  vdb_mean_ = pin->GetOrAddBoolean("problem", "vdamp_bot_mean_only", false);
  if (vdb_cells_ > 0) {
    if (2*vdb_cells_ >= pmy_mesh_->mesh_indcs.nx1) {
      std::cout << "### FATAL ERROR in red_giant: problem/vdamp_bot_cells = "
                << vdb_cells_ << " must be < nx1/2 = "
                << pmy_mesh_->mesh_indcs.nx1/2 << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(vdb_time_ > 0.0)) {
      std::cout << "### FATAL ERROR in red_giant: problem/vdamp_bot_time must be > 0"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (vdb_mean_ && pmy_mesh_->mesh_indcs.nx1 != pmy_mesh_->mb_indcs.nx1) {
      std::cout << "### FATAL ERROR in red_giant: problem/vdamp_bot_mean_only averages "
                << "over a shell at a GLOBAL radial index, so every MeshBlock must span "
                << "the whole radius: set meshblock/nx1 = mesh/nx1." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // --- THE TWO DIAGNOSTIC DUMPS (problem/rt_surface_dt, problem/rt_profile_dt), ported
  // from box_convection.cpp:2153-2178 and reinterpreted for a sphere.  See the
  // declarations of surf_file_ and prof_file_ for the record layouts.  Both are OFF at
  // dt <= 0, which is the default, and neither writes back to the state.
  surf_dt_ = pin->GetOrAddReal("problem", "rt_surface_dt", 0.0);
  {
    const std::string sf = pin->GetOrAddString("problem", "rt_surface_file",
                                               "rt_surface.bin");
    if (sf.size() >= sizeof(surf_file_)) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_surface_file is longer than "
                << "255 characters" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::snprintf(surf_file_, sizeof(surf_file_), "%s", sf.c_str());
  }
  surf_next_ = -1.0;
  surf_alloc_ = false;
  prof_dt_ = pin->GetOrAddReal("problem", "rt_profile_dt", 0.0);
  {
    const std::string pf = pin->GetOrAddString("problem", "rt_profile_file",
                                               "rt_profile.bin");
    if (pf.size() >= sizeof(prof_file_)) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_profile_file is longer than "
                << "255 characters" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::snprintf(prof_file_, sizeof(prof_file_), "%s", pf.c_str());
  }
  prof_next_ = -1.0;
  prof_alloc_ = false;
  // the surface dump takes the i = ie+1 face of every block as the top of the star, and
  // the profile indexes cells by their global radial position: both need whole columns
  if ((surf_dt_ > 0.0 || prof_dt_ > 0.0) &&
      pmy_mesh_->mesh_indcs.nx1 != pmy_mesh_->mb_indcs.nx1) {
    std::cout << "### FATAL ERROR in red_giant: problem/rt_surface_dt and "
              << "rt_profile_dt need every MeshBlock to span the whole radius: set "
              << "meshblock/nx1 = mesh/nx1." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // --- the CORRELATED-K two-stream for the optically thin layers
  // (utils/two_stream_rt.hpp)
  rt_ck_ = pin->GetOrAddBoolean("problem", "rt_ck", false);
  rt_grey_ = pin->GetOrAddBoolean("problem", "rt_grey", false);
  if (rt_ck_ && rt_grey_) {
    std::cout << "### FATAL ERROR in red_giant: problem/rt_ck and problem/rt_grey are "
              << "two solvers for the same layers. Pick one." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (rt_grey_) {
    namespace ts = two_stream_rt;
    ts::rt_grey = true;
    ts::rt_split = true;             // implied: the grey kernel lives on the split path
    // THE MERGED IMPLICIT COLUMN SOLVE (problem/rt_implicit_column), ported from
    // box_convection.cpp:1844-1921.  3 = the EXACT block-tridiagonal solve of the whole
    // radial column; 0 (the default) is the old per-cell semi-implicit apply and leaves
    // this problem generator bitwise unchanged.  It is read FIRST because the mode-3
    // gate (two_stream_rt.hpp:1386-1414) demands rt_use_cons, rt_src_direct,
    // rt_semi_implicit, !rt_explicit, !rt_layer_legacy, !rt_top_re and rt_outer_iter = 1,
    // i.e. the OPPOSITE of four of red_giant's legacy defaults: so every one of those
    // defaults is keyed on col3 below.  Each is still a GetOrAdd, so the input file
    // always wins and a mode-0 run reads exactly the values it read before.
    ts::rt_implicit_column = pin->GetOrAddInteger("problem", "rt_implicit_column", 0);
    const bool col3 = (ts::rt_implicit_column == 3);
    ts::rt_de_max = pin->GetOrAddReal("problem", "rt_de_max", 0.5);
    // these five default to the OLD (pre-fix) solver under rt_implicit_column = 0, as
    // they must for every problem generator sharing two_stream_rt.hpp; the red-giant
    // inputs turn them on explicitly.  Under mode 3 they default to the FIXED solver,
    // which is what that gate requires anyway (see col3 above).
    ts::rt_semi_lin = pin->GetOrAddBoolean("problem", "rt_semi_lin", !col3);
    ts::rt_explicit = pin->GetOrAddBoolean("problem", "rt_explicit", false);
    ts::rt_newton = pin->GetOrAddBoolean("problem", "rt_newton", col3);
    ts::rt_rescue_eq = pin->GetOrAddBoolean("problem", "rt_rescue_eq", col3);
    ts::rt_ali_diag = pin->GetOrAddBoolean("problem", "rt_ali_diag", true);
    ts::rt_relax_sub = pin->GetOrAddInteger("problem", "rt_relax_sub", 1);
    ts::rt_relax_xcrit = pin->GetOrAddReal("problem", "rt_relax_xcrit", 1.0);
    ts::rt_relax_submax = pin->GetOrAddInteger("problem", "rt_relax_submax", 32);
    ts::rt_src_direct = pin->GetOrAddBoolean("problem", "rt_src_direct", col3);
    ts::rt_top_clamp = pin->GetOrAddBoolean("problem", "rt_top_clamp", col3);
    ts::rt_use_cons = pin->GetOrAddBoolean("problem", "rt_use_cons", col3);
    ts::rt_bface = pin->GetOrAddBoolean("problem", "rt_bface", false);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ts::rt_int_at_cut = false;       // the same handover argument as rt_ck below
    // The column above the domain radiates to space instead of mirroring the ghost back
    // down.  Default ON here because a star has nothing above it: with outer_bc = open
    // the mirroring column sealed the atmosphere -- 6650 K isothermal, emergent flux
    // 0.3 % of L.  Set problem/rt_top_re = false to get the old boundary back.
    ts::rt_top_re = pin->GetOrAddBoolean("problem", "rt_top_re", !col3);
    // problem/rt_top_vacuum (default false = the unresolved-column model above): nothing
    // above the domain, I_down = 0 at the top face.  The right boundary when the domain
    // top sits in optically thin gas (the He4 presupernova star: tau = 1e-2 at x1max),
    // and the box's production setting; here it was never read, so the mode-3 column,
    // with rt_top_re forced off, always used the column model (2026-09-18).
    ts::rt_top_vacuum = pin->GetOrAddBoolean("problem", "rt_top_vacuum", false);
    // see two_stream_rt.hpp, rt_cut_bc_legacy: the grey sweep used to start the upward
    // intensity at the cut isotropic at B, which is exactly half the flux the column
    // supports.  The fix is the default here too; true restores the old boundary.
    ts::rt_cut_bc_legacy = pin->GetOrAddBoolean("problem", "rt_cut_bc_legacy", false);
    // see two_stream_rt.hpp, rt_layer_legacy: every layer used to span a whole
    // cell in optical depth but carry a source running between two cell CENTRES.
    // The fix is the default here too; true restores the old layers.
    ts::rt_layer_legacy = pin->GetOrAddBoolean("problem", "rt_layer_legacy", false);
    ts::rt_cell_report = pin->GetOrAddBoolean("problem", "rt_cell_report", false);
    ts::rt_report_r = pin->GetOrAddReal("problem", "rt_report_r", 3.887e12);
    ts::rt_report_every = pin->GetOrAddInteger("problem", "rt_report_every", 100);
    ts::rt_tint_override = teff;
    ts::rt_star_teff = 0.0;
    ts::rt_dump_file = pin->GetOrAddString("problem", "ck_dump_file", "");
    ts::rt_dump_m = pin->GetOrAddInteger("problem", "ck_dump_m", 0);
    ts::rt_dump_j = pin->GetOrAddInteger("problem", "ck_dump_j", -1);
    ts::rt_dump_k = pin->GetOrAddInteger("problem", "ck_dump_k", -1);
    ts::rt_use_cons = pin->GetOrAddBoolean("problem", "rt_use_cons", col3);
    ts::rt_bface = pin->GetOrAddBoolean("problem", "rt_bface", false);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ck_dump_t2_ = pin->GetOrAddReal("problem", "ck_dump_t2", -1.0);
    ck_dump_file2_ = pin->GetOrAddString("problem", "ck_dump_file2", "");
    // ck_nquad picks the angular quadrature for the grey sweep too: 1 = hemispheric
    // mean, 2 = two-point Gauss
    correlated_k::ck_nq = pin->GetOrAddInteger("problem", "ck_nquad", 1);
    // ---- the column solve's own controls, ported verbatim from box_convection.cpp
    // (:1844-1921).  Every name, and every default, is the box's, so the same input
    // line means the same thing in both problem generators.
    ts::rt_semi_implicit = pin->GetOrAddBoolean("problem", "rt_semi_implicit", true);
    ts::rt_outer_iter = pin->GetOrAddInteger("problem", "rt_outer_iter", 1);
    ts::rt_outer_verbose = pin->GetOrAddBoolean("problem", "rt_outer_verbose", false);
    ts::rt_impl_tol = pin->GetOrAddReal("problem", "rt_impl_tol", col3 ? 1.0e-8 : 1.0e-6);
    ts::rt_col3_ex_iter = pin->GetOrAddBoolean("problem", "rt_col3_ex_iter", false);
    // mode 3 converges in 2-4 Newton steps and is cheap per step, so it gets one more
    // than the linearised modes by default
    ts::rt_impl_maxit = pin->GetOrAddInteger("problem", "rt_impl_maxit", col3 ? 8 : 5);
    ts::rt_impl_exjac = pin->GetOrAddBoolean("problem", "rt_impl_exjac", true);
    ts::rt_impl_norm = pin->GetOrAddInteger("problem", "rt_impl_norm", 1);
    ts::rt_impl_norm_eps = pin->GetOrAddReal("problem", "rt_impl_norm_eps", 1.0e-3);
    ts::rt_impl_dstop = pin->GetOrAddBoolean("problem", "rt_impl_dstop", true);
    ts::rt_impl_rescheck = pin->GetOrAddBoolean("problem", "rt_impl_rescheck", true);
    ts::rt_impl_ablate = pin->GetOrAddInteger("problem", "rt_impl_ablate", 0);
    ts::rt_impl_fixit = pin->GetOrAddBoolean("problem", "rt_impl_fixit", false);
    if (ts::rt_impl_ablate != 0 || ts::rt_impl_fixit) {
      std::cout << "### WARNING in red_giant: problem/rt_impl_ablate or rt_impl_fixit "
                << "is set.  These are TIMING INSTRUMENTATION and the mode-3 solve they "
                << "produce is NOT a correct solve." << std::endl;
    }
    ts::rt_impl_cvfreeze = pin->GetOrAddInteger("problem", "rt_impl_cvfreeze", 0);
    // problem/rt_impl_mixed: 0 = double (default), 1 = single-precision Newton
    // correction, 2 = and single-precision stored factors.  See two_stream_rt.hpp.
    ts::rt_impl_mixed = pin->GetOrAddInteger("problem", "rt_impl_mixed", 0);
    if (ts::rt_impl_mixed < 0 || ts::rt_impl_mixed > 2) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_impl_mixed must be 0, 1 or 2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/rt_impl_reuse: reuse the block factorisation across Newton passes, with a
    // contraction check (1) or unconditionally (2).  0 (the default) is the old code.
    ts::rt_impl_reuse = pin->GetOrAddInteger("problem", "rt_impl_reuse", 0);
    ts::rt_impl_reuse_rho = pin->GetOrAddReal("problem", "rt_impl_reuse_rho", 0.3);
    if (ts::rt_impl_reuse < 0 || ts::rt_impl_reuse > 2) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_impl_reuse must be 0, 1 or 2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // how the mode-3 block system is solved: the serial block Thomas (one thread per
    // column) or the team-partitioned solve (two_stream_column_partition.hpp)
    {
      const std::string sv = pin->GetOrAddString("problem", "rt_impl_solver", "thomas");
      if (sv.compare("thomas") == 0) {
        ts::rt_impl_solver = 0;
      } else if (sv.compare("pcr") == 0) {
        ts::rt_impl_solver = 1;
      } else {
        std::cout << "### FATAL ERROR in red_giant: problem/rt_impl_solver = " << sv
                  << " is not one of thomas, pcr" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    ts::rt_impl_redpar = pin->GetOrAddBoolean("problem", "rt_impl_redpar", false);
    ts::rt_col3_hybrid_tau = pin->GetOrAddReal("problem", "rt_col3_hybrid_tau", 0.0);
    ts::rt_col3_split_deep = pin->GetOrAddBoolean("problem", "rt_col3_split_deep", false);
    ts::rt_col3_split_w = pin->GetOrAddInteger("problem", "rt_col3_split_w", 8);
    // segments per column = the Kokkos team size of the PCR partitioned solve
    ts::rt_impl_nseg = pin->GetOrAddInteger("problem", "rt_impl_nseg", 64);
    if (ts::rt_impl_nseg < 1) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_impl_nseg must be >= 1"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/rt_impl_warm: warm-start the mode-3 Newton from the previous call's
    // converged Planck function (1) or from a linear extrapolation of the last two (2).
    ts::rt_impl_warm = pin->GetOrAddInteger("problem", "rt_impl_warm", 0);
    if (ts::rt_impl_warm < 0 || ts::rt_impl_warm > 2) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_impl_warm must be 0, 1 or 2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_impl_tau_min = pin->GetOrAddReal("problem", "rt_impl_tau_min", 1.0);
    ts::rt_impl_dtmax = pin->GetOrAddReal("problem", "rt_impl_dtmax", 0.25);
    ts::rt_impl_tau_blend = pin->GetOrAddReal("problem", "rt_impl_tau_blend", 1.0);
    ts::rt_col3_skip_sweep = pin->GetOrAddBoolean("problem", "rt_col3_skip_sweep", false);
    if (ts::rt_col3_skip_sweep && !(col3 && ts::rt_col3_ex_iter && ts::rt_src_direct)) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_col3_skip_sweep needs "
                << "rt_implicit_column = 3, rt_col3_ex_iter and rt_src_direct (the "
                << "column solve then supplies the face flux from its own converged "
                << "intensities)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (col3) {
      // the column solve is radial and private to a team, so one MeshBlock must hold
      // the whole x1 extent, and the private intensity column is sized at compile time
      // (box_convection.cpp:1778-1793, the same two checks)
      if (pmy_mesh_->mb_indcs.nx1 != pmy_mesh_->mesh_indcs.nx1) {
        std::cout << "### FATAL ERROR in red_giant: problem/rt_implicit_column = 3 "
                  << "solves a whole RADIAL column inside one MeshBlock, but mesh/nx1 = "
                  << pmy_mesh_->mesh_indcs.nx1 << " and meshblock/nx1 = "
                  << pmy_mesh_->mb_indcs.nx1 << ". Set meshblock/nx1 = mesh/nx1 and "
                  << "decompose in x2/x3 only." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (indcs.nx1 + 2*ng > 520) {
        std::cout << "### FATAL ERROR in red_giant: the grey column sweep dispatches on "
                  << "compile-time column tiers, the largest of which is 520, but nx1 = "
                  << indcs.nx1 << " and nghost = " << ng << " give "
                  << (indcs.nx1 + 2*ng) << "." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // ---- the thin-region RADIATIVE FORCE (box_convection.cpp:2185-2320) -------------
    // The other half of the EOS's radiation taper: the taper removes (1-w) of the LTE
    // radiation pressure from the gas, and this puts the force that pressure carried
    // back as an explicit momentum source, (1-w) rho kappa F/c + Prad grad w.  Neither
    // half is meaningful alone, so both are required together.  Mandatory for the He
    // star, where kappa F/(c g) = 1.23 -- the envelope is super-Eddington.
    //
    // The g(r) the force works against is NOT a parameter of the source: the momentum
    // term is built from the two-stream's own face flux and the taper weight, and the
    // point-mass gravity it opposes is applied by this file's own rg_grav kernel.  The
    // only place the solver needs a gravity is EffGravAt, for the unresolved column
    // above the top face, and that already reads hot_jupiter_param.grav with
    // grav_point_mass = true (set below), i.e. the true g(r) = G M / r^2.  ts::
    // rt_force_grav is used ONLY to normalise the rt_force_verbose print, so it is set
    // to g at the inner wall, the same reference GravAt uses.
    ts::rt_rad_force = pin->GetOrAddBoolean("problem", "rt_rad_force", false);
    ts::rt_force_verbose = pin->GetOrAddInteger("problem", "rt_force_verbose", 0);
    ts::rt_force_grav = kGrav*mstar/SQR(rin_*lunit);
    if (ts::rt_rad_force) {
      if (!pmbp->phydro->peos->eos_data.tbl.rad_taper) {
        std::cout << "### FATAL ERROR in red_giant: problem/rt_rad_force is the momentum "
                  << "source that goes with the EOS radiation taper, and is only defined "
                  << "when the taper is on: set <hydro>/eos_rad_rho_hi and "
                  << "eos_rad_rho_lo (with eos_radiation = true), or drop rt_rad_force."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (global_variable::my_rank == 0) {
        std::cout << "red_giant: RT radiative momentum source ON -- (1-w) rho kappa F/c "
                  << "+ Prad grad w, with w from the EOS taper; g(r) = GM/r^2"
                  << std::endl;
      }
    }
    // the TIME CENTRING of the two coupling terms mode 3 leaves first order.  Both are
    // default-off and bitwise off; both need the exact column solve.
    ts::rt_force_center = pin->GetOrAddInteger("problem", "rt_force_center", 0);
    if (ts::rt_force_center < 0 || ts::rt_force_center > 2) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_force_center must be 0 "
                << "(entry flux), 1 (converged flux) or 2 (the average)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ts::rt_force_center > 0 && !col3) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_force_center needs "
                << "problem/rt_implicit_column = 3 (only the exact column solve has a "
                << "converged flux to centre against)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ts::rt_src_theta = pin->GetOrAddReal("problem", "rt_src_theta", 1.0);
    if (ts::rt_src_theta < 0.0 || ts::rt_src_theta > 1.0) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_src_theta must be in [0, 1]"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (ts::rt_src_theta != 1.0 && !col3) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_src_theta != 1 needs "
                << "problem/rt_implicit_column = 3" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // ---- problem/rt_bottom_flux (box_convection.cpp:2214) --------------------------
    // Hand the internal flux to the TWO-STREAM's own lower boundary instead of to the
    // conduction wall face.  Only consistent when the tau blend leaves the sweep the
    // whole column, i.e. <hydro>/rad_tau_lo deeper than the domain so the blend weight
    // is 0 on every face; that cannot be checked here (the weights are built per cycle),
    // so it is stated, not enforced.  Meaningless with inner_bc = open, which injects
    // the luminosity as the entropy of the inflow instead.
    {
      const bool rt_botflux = pin->GetOrAddBoolean("problem", "rt_bottom_flux", false);
      if (rt_botflux) {
        if (pc == nullptr || !pc->rad_tau_mode) {
          std::cout << "### FATAL ERROR in red_giant: problem/rt_bottom_flux puts the "
                    << "internal flux on the two-stream's lower boundary and takes it "
                    << "off the conduction wall face, which is only consistent when the "
                    << "sweep reaches that wall: <hydro>/rad_tau_hi must be set."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        ts::rt_bot_flux = pc->rad_flux_inner;
        pc->rad_flux_inner = 0.0;
        if (global_variable::my_rank == 0) {
          std::printf("red_giant: F_bot = %.5e is carried by the TWO-STREAM's lower "
                      "boundary (rt_bottom_flux); the conduction wall face injects "
                      "nothing.  This assumes <hydro>/rad_tau_lo = %.3e is DEEPER than "
                      "the domain, so the blend weight is 0 on every face.\n",
                      ts::rt_bot_flux, pc->rad_tau_lo);
        }
      }
    }
    // ---- OPERATOR-SPLIT SCHEDULING (box_convection.cpp:1922-2110) -------------------
    // problem/rt_strang: the whole radiation operator is Strang-split around the time
    // integrator (dt/2 before, dt/2 after) instead of being applied once per RK stage;
    // problem/rt_once_per_cycle applies it once with the full dt after the last stage.
    // Both take the source OUT of the stage, so the in-stage call below is skipped.
    // problem/rt_imex is deliberately NOT ported.
    rt_strang_ = pin->GetOrAddBoolean("problem", "rt_strang", false);
    rt_once_ = pin->GetOrAddBoolean("problem", "rt_once_per_cycle", false);
    if (rt_strang_ && rt_once_) {
      std::cout << "### FATAL ERROR in red_giant: problem/rt_strang and "
                << "problem/rt_once_per_cycle are mutually exclusive" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_strang_ || rt_once_) {
      user_split_func = RedGiantRTSplit;
      user_split_once = rt_once_;
      // problem/rt_split_transverse: with the column out of the stage, the implicit
      // transverse (ADI) operator moves WITH it, over the same bdt and on the state the
      // column solve has just relaxed.  Default ON whenever a split is active.
      rt_split_tr_ = pin->GetOrAddBoolean("problem", "rt_split_transverse", true);
      if (rt_split_tr_ && pc != nullptr && pc->rad_implicit_ang) {
        pc->rad_tr_split_out = true;
      }
      if (global_variable::my_rank == 0) {
        std::cout << "red_giant: the two-stream is applied "
                  << (rt_strang_ ? "STRANG-SPLIT around the time integrator (dt/2 "
                                   "before, dt/2 after)"
                                 : "ONCE per cycle with the full dt after the last RK "
                                   "stage")
                  << " and NOT inside the RK stages"
                  << (rt_split_tr_ ? "; the transverse operator moves with it" : "")
                  << std::endl;
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: GREY two-stream (problem/rt_grey), opacity from the "
                << "conduction module's table, " << pin->GetOrAddInteger("problem",
                   "ck_nquad", 1) << "-point angular quadrature" << std::endl;
      std::cout << "           column above the domain: "
                << (two_stream_rt::rt_top_re
                    ? "radiative equilibrium, it radiates to space"
                    : "the ghost's own Planck function (it mirrors the top cell back)")
                << std::endl;
    }
  }
  if (rt_ck_) {
    namespace ck = correlated_k;
    namespace ts = two_stream_rt;
    ts::rt_ck = true;
    ts::rt_split = true;             // implied by rt_ck
    ts::rt_ck_pcut = pin->GetOrAddReal("problem", "ck_pcut_bar", 1.0e30);
    ts::rt_de_max = pin->GetOrAddReal("problem", "rt_de_max", 0.5);
    // these five default to the OLD (pre-fix) solver, as they must for every problem
    // generator sharing two_stream_rt.hpp; the red-giant inputs turn them on explicitly
    ts::rt_semi_lin = pin->GetOrAddBoolean("problem", "rt_semi_lin", true);
    ts::rt_explicit = pin->GetOrAddBoolean("problem", "rt_explicit", false);
    ts::rt_newton = pin->GetOrAddBoolean("problem", "rt_newton", false);
    ts::rt_rescue_eq = pin->GetOrAddBoolean("problem", "rt_rescue_eq", false);
    ts::rt_ali_diag = pin->GetOrAddBoolean("problem", "rt_ali_diag", true);
    ts::rt_relax_sub = pin->GetOrAddInteger("problem", "rt_relax_sub", 1);
    ts::rt_relax_xcrit = pin->GetOrAddReal("problem", "rt_relax_xcrit", 1.0);
    ts::rt_relax_submax = pin->GetOrAddInteger("problem", "rt_relax_submax", 32);
    ts::rt_src_direct = pin->GetOrAddBoolean("problem", "rt_src_direct", false);
    ts::rt_top_clamp = pin->GetOrAddBoolean("problem", "rt_top_clamp", false);
    ts::rt_use_cons = pin->GetOrAddBoolean("problem", "rt_use_cons", false);
    ts::rt_bface = pin->GetOrAddBoolean("problem", "rt_bface", false);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    // THE SOLVER INJECTS NO ENERGY, and cannot: its upward intensity at the cut is the
    // local band Planck function, which is the right handover to an optically thick,
    // diffusive interior.  The internal-flux term is off twice over -- set false here,
    // and forced false by the tau blend, which this problem generator requires above.
    // So the only energy input is the inner boundary: rad_flux_inner through a wall, or
    // the entropy of the inflow through an open one.
    ts::rt_int_at_cut = false;
    // T_int reaches only the grey picket-fence path and the albedo, neither of which runs
    // with rt_ck and the blend on, and the albedo multiplies a stellar flux that is zero.
    // It is set to a sane value rather than left at 0 so those helpers, which are
    // evaluated unconditionally, are not handed a degenerate temperature.
    ts::rt_tint_override = teff;
    ts::rt_dump_file = pin->GetOrAddString("problem", "ck_dump_file", "");
    ts::rt_dump_m = pin->GetOrAddInteger("problem", "ck_dump_m", 0);
    ts::rt_dump_j = pin->GetOrAddInteger("problem", "ck_dump_j", -1);
    ts::rt_dump_k = pin->GetOrAddInteger("problem", "ck_dump_k", -1);
    ts::rt_use_cons = pin->GetOrAddBoolean("problem", "rt_use_cons", false);
    ts::rt_bface = pin->GetOrAddBoolean("problem", "rt_bface", false);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ck_dump_t2_ = pin->GetOrAddReal("problem", "ck_dump_t2", -1.0);
    ck_dump_file2_ = pin->GetOrAddString("problem", "ck_dump_file2", "");
    ts::rt_star_teff = 0.0;          // no host: the stellar band fractions are unused
    ts::rt_nchain = ck::CK_NB*ck::CK_NG*pin->GetOrAddInteger("problem", "ck_nquad", 1);
    ck::ck_nq = pin->GetOrAddInteger("problem", "ck_nquad", 1);
    ck::read_ck_table(pin->GetString("problem", "ck_table"), ts::rt_ck_pcut);
    ck::build_planck_fractions(ts::rt_ck_pcut);
    ck::read_ck_continuum(pin->GetString("problem", "ck_data_dir"),
                          pin->GetOrAddString("problem", "ck_swflux",
                                              "sw_band_flux_W121_11.txt"), 0.0);
    ck::ck_selftest();
    ck::ck_rt_selftest();
    // problem/opac_compare: write the Rosseland mean DERIVED FROM THE CK TABLE on its own
    // (T, p) grid, so it can be compared against the stellar table the diffusion uses.
    // The blend hands over between the two, and they are different data for the same gas:
    // if they disagree where they overlap, the handover is not conservative.  Built into
    // the conduction object and then undone, since the run must keep the stellar table.
    {
      const std::string ocmp = pin->GetOrAddString("problem", "opac_compare", "");
      if (!ocmp.empty() && pc != nullptr && pc->rad_kappa_tab) {
        auto sv_tab = pc->rad_kr_tab;
        auto sv_lT = pc->rad_kr_lT;
        auto sv_lP = pc->rad_kr_lP;
        const int sv_nT = pc->rad_kr_nT, sv_nP = pc->rad_kr_nP;
        const bool sv_rho = pc->rad_kappa_rho;
        pc->rad_kappa_rho = false;
        ck::ck_build_rosseland_table(pc);
        auto hk = Kokkos::create_mirror_view(pc->rad_kr_tab);
        auto hT = Kokkos::create_mirror_view(pc->rad_kr_lT);
        auto hP = Kokkos::create_mirror_view(pc->rad_kr_lP);
        Kokkos::deep_copy(hk, pc->rad_kr_tab);
        Kokkos::deep_copy(hT, pc->rad_kr_lT);
        Kokkos::deep_copy(hP, pc->rad_kr_lP);
        if (global_variable::my_rank == 0) {
          std::ofstream f(ocmp);
          f.precision(10);
          f << "# Rosseland mean from the CORRELATED-K table + continuum\n"
            << "# nT nP, then nT values of log10 T[K], nP of log10 p[dyn/cm^2], then\n"
            << "# nT*nP values of log10 kappa_R [cm^2/g], T slowest\n"
            << pc->rad_kr_nT << " " << pc->rad_kr_nP << "\n";
          for (int i = 0; i < pc->rad_kr_nT; ++i) f << hT(i) << "\n";
          for (int j = 0; j < pc->rad_kr_nP; ++j) f << hP(j) << "\n";
          for (int i = 0; i < pc->rad_kr_nT; ++i) {
            for (int j = 0; j < pc->rad_kr_nP; ++j) f << hk(i, j) << "\n";
          }
          std::cout << "red_giant: ck-derived Rosseland table written to '" << ocmp
                    << "'" << std::endl;
        }
        pc->rad_kr_tab = sv_tab;
        pc->rad_kr_lT = sv_lT;
        pc->rad_kr_lP = sv_lP;
        pc->rad_kr_nT = sv_nT;
        pc->rad_kr_nP = sv_nP;
        pc->rad_kappa_rho = sv_rho;
      }
    }
    std::cout << "red_giant: correlated-k two-stream ON, T_int = " << teff
              << " K, no irradiation" << std::endl;
  }
  // The star and grid numbers BOTH band solvers read; Teq = 0 switches the stellar sweep
  // off.  This used to be set only on the correlated-k path, which left the grey path
  // with grav = 0: the top-ghost column's dtau = kappa p / g_eff was then +inf, which the
  // exponential absorbed silently, and became 0/0 = NaN the moment a transparent top cell
  // made kappa zero.  The solver needs the gravity whichever band scheme is running.
  if (rt_ck_ || rt_grey_) {
    hot_jupiter_param.Teq = 0.0;
    hot_jupiter_param.omega = 0.0;
    hot_jupiter_param.grav = kGrav*mstar/SQR(rin*lunit);
    hot_jupiter_param.ap = rin;
    hot_jupiter_param.Rgas = rgas_;
    hot_jupiter_param.met = pin->GetOrAddReal("problem", "met", 0.0);
    hot_jupiter_param.grav_point_mass = true;
    hot_jupiter_param.stellar_tide = false;
    hot_jupiter_param.rot_potential = false;
  }
  ktab_ = ktab; klT_ = klT; klD_ = klD; knT_ = knT; knD_ = knD; teff_ = teff;
  kfac_ = kfac;
  if (mlt_alpha_ > 0.0) {
    Kokkos::realloc(fconv_, pmbp->nmb_thispack, n3m1+1, n2m1+1, n1m1+2);
    Kokkos::realloc(taumlt_, pmbp->nmb_thispack, n3m1+1, n2m1+1, n1m1+2);
    Kokkos::realloc(fdiag_, n1m1+2, 8);
    Kokkos::realloc(shell_, n1m1+2, 13);
    Kokkos::realloc(fmlt1d_, n1m1+2);
    Kokkos::realloc(fmean_, n1m1+2, 17);
    // the relaxed profile goes into, and comes out of, the restart file (see
    // RedGiantMltRestartState).  Enrolled only here, so a run with problem/mlt_alpha = 0
    // writes no block and every other pgen's restart format is untouched.
    pgen_rst_write_func = RedGiantMltRestartState;
    if (restart) RedGiantMltRestartRestore(pgen_rststate);
  }

  // --- the initial column: fine grid in r from below the inner ghosts to above the
  // outer ones (10 % margins cover any stretch), integrated from the outer wall
  const int nfine = pin->GetOrAddInteger("problem", "nfine",
                                        40*pmy_mesh_->mesh_indcs.nx1 + 2*ng*40);
  const Real rlo = rin - 0.1*(rout - rin), rhi = rout + 0.1*(rout - rin);
  const Real drf = (rhi - rlo)/(nfine - 1);
  // The column is integrated DOWNWARD from itop, so itop -- not the grid -- is what
  // fixes the star's size.  Without a background it is the outer wall, as before.
  const Real ranch = (bg_rtop > 0.0) ? bg_rtop/lunit : rout;
  const int itop = static_cast<int>((ranch - rlo)/drf + 0.5);
  DvceArray1D<Real> lnp("rg_lnp", nfine), tk("rg_tk", nfine);
  DvceArray1D<Real> kap("rg_kap", nfine), grad("rg_grad", nfine), tau("rg_tau", nfine);
  DvceArray1D<Real> vcc("rg_vc", nfine), rhoa("rg_rho", nfine);
  // --- problem/ic_profile: (rho, eint) on the fine grid, read on the host.  Both are
  // interpolated LOG-LINEARLY in r, as box_convection does: the profile spans decades in
  // both and a linear interpolant on the coarse file nodes would not be hydrostatic to
  // the accuracy the file itself has.  Nodes of the fine grid outside the file's range
  // -- the 10 % margins, which no cell and no ghost ever reaches -- continue the end
  // segment; ColumnAt clamps there in any case.
  const bool useprof = !icprof.empty();
  DualArray1D<Real> pfd("rg_pfd", nfine), pfe("rg_pfe", nfine);
  if (useprof) {
    std::ifstream pf(icprof);
    if (!pf.good()) {
      std::cout << "### FATAL ERROR in red_giant: cannot open problem/ic_profile '"
                << icprof << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::vector<Real> rf, df, ef;
    std::string line;
    while (std::getline(pf, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      Real a, b, c;
      if (!(is >> a >> b >> c)) continue;
      rf.push_back(a/lunit);
      df.push_back(b/dunit);
      ef.push_back(c/punit_);
    }
    if (rf.size() < 2) {
      std::cout << "### FATAL ERROR in red_giant: problem/ic_profile '" << icprof
                << "' has " << rf.size() << " usable rows" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (std::size_t q = 1; q < rf.size(); ++q) {
      if (!(rf[q] > rf[q-1])) {
        std::cout << "### FATAL ERROR in red_giant: problem/ic_profile must be in "
                  << "ASCENDING r; row " << q << " is r = " << rf[q]*lunit
                  << " after " << rf[q-1]*lunit << " cm" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // COVERAGE is checked on the RADII THE GRID ACTUALLY HAS, ghosts included -- not on
    // the fine grid's 10 % margins, which are scratch.  The ghosts matter because the
    // wall and open radial boundaries fill them from this column.
    {
      auto &x1v_ic = pmbp->pcoord->x1v;
      auto hx1v_ic = Kokkos::create_mirror_view(x1v_ic);
      Kokkos::deep_copy(hx1v_ic, x1v_ic);
      Real rgmin = std::numeric_limits<Real>::max();
      Real rgmax = -std::numeric_limits<Real>::max();
      for (int m = 0; m <= nmb1; ++m) {
        for (int i = 0; i <= n1m1; ++i) {
          const Real rr = RadiusOf(curv, hx1v_ic(m,i), rin, x1min_);
          rgmin = fmin(rgmin, rr);
          rgmax = fmax(rgmax, rr);
        }
      }
      if (rf.front() > rgmin || rf.back() < rgmax) {
        std::cout << "### FATAL ERROR in red_giant: problem/ic_profile spans ["
                  << rf.front()*lunit << ", " << rf.back()*lunit << "] cm, which does "
                  << "not cover the radial grid plus its ghosts ["
                  << rgmin*lunit << ", " << rgmax*lunit << "] cm" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    std::size_t kk = 0;
    for (int i = 0; i < nfine; ++i) {
      const Real r = rlo + i*drf;
      while (kk + 2 < rf.size() && rf[kk+1] < r) ++kk;
      const Real w = (r - rf[kk])/(rf[kk+1] - rf[kk]);
      pfd.h_view(i) = std::exp(std::log(df[kk])*(1.0 - w) + std::log(df[kk+1])*w);
      pfe.h_view(i) = std::exp(std::log(ef[kk])*(1.0 - w) + std::log(ef[kk+1])*w);
    }
    pfd.modify_host();  pfd.sync_device();
    pfe.modify_host();  pfe.sync_device();
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: initial stratification READ FROM " << icprof << " ("
                << rf.size() << " nodes, " << rf.front()*lunit << " .. "
                << rf.back()*lunit << " cm); the hydrostatic+MLT march, problem/ptop "
                << "and the outer-wall gas-solution check are all OVERRIDDEN"
                << std::endl;
    }
  }
  {
    const Real ttop = teff*pow(0.5, 0.25);
    const Real ptop = ptop_cgs/punit_;
    const Real lum = lstar, mass = mstar, lun = lunit, dun = dunit, pun = punit_;
    const Real alpha_ic = mlt_alpha_ic;
    const Real vunit_c = vunit;
    const Real tun = pmbp->punit->temperature_cgs();
    // one serial sweep on the device (the EOS and the table live there); RK2 in r
    // where the thin-layer relaxation acts (tau < rad_tau_hi) the column is RADIATIVE
    // whatever the Schwarzschild criterion says: the relaxation targets the Eddington
    // profile there, and an adiabatic start would be cooled toward it from the first
    // step (in a real star convection is inefficient that far out anyway)
    const Real tauhi_blend = (pc != nullptr && pc->rad_tau_mode) ? pc->rad_tau_hi : 0.0;
    const Real tauhi = (ic_tau_rad_ >= 0.0) ? ic_tau_rad_ : tauhi_blend;
    const Real dbg = bg_rho_/dunit, tbg = bg_temp_, gm1c = gm1;
    const bool bghse = bg_hydrostatic_;
    auto pfd_d = pfd.d_view;
    auto pfe_d = pfe.d_view;
    const bool useprof_k = useprof;
    par_for("rg_column", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int dummy) {
      // gradient at (p [code], T [K], tau): d ln T / d ln p and the state
      auto nabla = [&](const Real r, const Real p, const Real t, const Real ta,
                       Real &rho, Real &kr, Real &gr, Real &vc) {
        vc = 0.0;
        rho = DensFromPT(eos, rgas, p, t);
        kr = kfac*KappaTab(ktab, klT, klD, knT, knD, t, rho*dun);
        const Real grad_rad = 3.0*kr*(p*pun)*lum
                              /(16.0*M_PI*kArad*kClight*kGrav*mass*t*t*t*t);
        const Real grad_ad = GradAd(eos, gamma, rgas, p, t);
        if (ta < tauhi || grad_rad < grad_ad) {
          gr = grad_rad;                       // radiative, or the two-stream's domain
          return gr;
        }
        if (!(alpha_ic > 0.0)) {
          gr = grad_ad;             // perfectly efficient convection: no driving at all
          return gr;
        }
        // MIXING LENGTH.  grad_ad is the limit of infinitely efficient convection: it
        // carries the flux at zero superadiabaticity, so a column built on it is
        // marginally stable, has no convective flux, and cannot start convecting.  The
        // gradient that actually carries F = L/(4 pi r^2) splits it between radiation
        // and the same MLT closure the source term uses (see the rg_mlt_flux kernel),
        //   F (1 - grad/grad_rad) = rho cp T sqrt(g delta) l^2 x^(3/2)
        //                           / (4 sqrt2 Hp^(3/2)),   x = grad - grad_ad,
        // whose left side falls and right side rises in grad, so the root in
        // [grad_ad, grad_rad] is unique; bisect for it.  Deep down the convection is
        // efficient and it sits a hair above grad_ad, which is why this changes nothing
        // there and everything in the superadiabatic layer below the photosphere.
        const Real rcm = r*lun;
        const Real g = kGrav*mass/(rcm*rcm);                    // cgs
        const Real dcgs = rho*dun, pcgs = p*pun;
        const Real ftot = lum/(4.0*M_PI*rcm*rcm);               // erg/cm^2/s
        const Real hp = pcgs/(dcgs*g);
        const Real ell = alpha_ic*hp;
        const Real e = EintFromDensT(eos, rgas, igm1, rho, t);
        Real cv = (EintFromDensT(eos, rgas, igm1, rho, 1.01*t) - e)/(0.01*t)*pun/dcgs;
        if (!(cv > 0.0)) cv = 1.5*pcgs/(dcgs*t);
        Real delta = 1.0, cp = gamma*cv;
        if (eos.IsGeneral()) {
          const Real chit = eos.ChiT(rho, e), chir = eos.ChiRho(rho, e);
          delta = chit/chir;
          cp = cv*eos.Gamma1(rho, e)/chir;
        }
        const Real cmlt = dcgs*cp*t*sqrt(g*delta)*ell*ell/(4.0*sqrt(2.0)*hp*sqrt(hp));
        if (!(cmlt > 0.0) || !(ftot > 0.0)) {
          gr = grad_ad;
          return gr;
        }
        Real lo = grad_ad, hi = grad_rad;
        for (int it = 0; it < 60; ++it) {
          const Real mid = 0.5*(lo + hi);
          const Real x = mid - grad_ad;
          if (ftot*(1.0 - mid/grad_rad) > cmlt*x*sqrt(x)) {
            lo = mid;
          } else {
            hi = mid;
          }
        }
        gr = 0.5*(lo + hi);
        // the convective velocity that carries the rest of the flux, the same expression
        // the rg_mlt_flux source term uses.  This is what the seed is scaled to.
        vc = sqrt(g*delta*ell*ell*(gr - grad_ad)/(8.0*hp))/vunit_c;
        return gr;
      };
      // a hydrostatic RK2 step from (r, lnp, T, tau) over dr (signed; tau grows inward)
      auto step = [&](const Real r, const Real lp0, const Real t0, const Real dr,
                      Real &lp1, Real &t1, Real &ta) {
        Real rho, kr, gr, vcd;
        const Real p0 = exp(lp0);
        nabla(r, p0, t0, ta, rho, kr, gr, vcd);
        const Real dlnp_a = -rho*GravAt(gm, r)/p0;
        const Real lpm = lp0 + 0.5*dr*dlnp_a;
        const Real tm = t0*exp(0.5*dr*dlnp_a*gr);
        const Real rm = r + 0.5*dr;
        nabla(rm, exp(lpm), tm, ta, rho, kr, gr, vcd);
        const Real dlnp_m = -rho*GravAt(gm, rm)/exp(lpm);
        lp1 = lp0 + dr*dlnp_m;
        t1 = t0*exp(dr*dlnp_m*gr);
        ta -= kr*rho*dun*dr*lun;        // dr < 0 inward: tau increases
        if (ta < 0.0) ta = 0.0;
      };
      if (useprof_k) {
        // problem/ic_profile: the column IS the file.  Everything downstream -- the
        // initial state, the radial ghosts of both boundary modes, the static
        // well-balanced background, the sponges, the seed's tau window and column_dump
        // -- reads (lnp, tk) through ColumnAt, so converting the supplied (rho, eint)
        // here with THIS run's own EOS is all that is needed, and it is the run's own
        // (p, T) by construction.  No march, no anchor, no ptop.
        for (int i = 0; i < nfine; ++i) {
          Real pp_, tt_;
          pgen_eos::PresTempFromEint(eos, gm1c, rgas, pfd_d(i), pfe_d(i), -1.0, pp_, tt_);
          lnp(i) = log(pp_);
          tk(i) = tt_;
        }
      } else {
        lnp(itop) = log(ptop);
        tk(itop) = ttop;
        Real ta = 0.0;
        for (int i = itop+1; i < nfine; ++i) {
          step(rlo + (i-1)*drf, lnp(i-1), tk(i-1), drf, lnp(i), tk(i), ta);
        }
        ta = 0.0;
        for (int i = itop-1; i >= 0; --i) {
          step(rlo + (i+1)*drf, lnp(i+1), tk(i+1), -drf, lnp(i), tk(i), ta);
        }
      }
      // --- THE AMBIENT MEDIUM.  Everything above the point where the star's own column
      // has thinned to bg_rho becomes either a constant (bg_rho, bg_temp) medium or, with
      // bg_hydrostatic, an isothermal corona at bg_temp; bg_rho is only the threshold
      // that locates the join in the second case.  For the constant medium, joining on
      // DENSITY
      // rather than on pressure is deliberate: an equal-pressure join would put gas that
      // is (T_star/T_bg) times heavier on top of the atmosphere and be Rayleigh-Taylor
      // unstable from the first step.  Matching density instead leaves the background
      // under-pressured by that same ratio, which is the physical situation -- the
      // atmosphere of a giant is not confined by its surroundings -- and it expands
      // slightly into the background rather than the background falling into it.
      if (dbg > 0.0) {
        const Real ebg = EintFromDensT(eos, rgas, igm1, dbg, tbg);
        const Real pbg = PresFromEint(eos, gm1c, dbg, ebg);
        // Scan UPWARD, not down from the top.  Above the point where the hydrostatic
        // column runs out of atmosphere the tabulated EOS stops inverting and
        // DensFromPT returns its 1e6 g/cm^3 failure sentinel with T = NaN; a downward
        // scan breaks on that sentinel immediately and writes no background at all.
        // Going up, the first level that is not both finite and denser than bg_rho IS
        // the join, whichever of the two ends the star.
        // from the ANCHOR upward: below itop the column is the star by construction, and
        // the bottom of the fine grid sits at a negative radius where the inversion
        // returns the same sentinel the top does.
        int ijoin = nfine;
        for (int i = itop; i < nfine; ++i) {
          const Real rho_i = DensFromPT(eos, rgas, exp(lnp(i)), tk(i));
          if (!(rho_i > dbg) || !(rho_i < 1.0e5) || !(tk(i) > 0.0)) { ijoin = i; break; }
        }
        if (bghse) {
          // THE HOT CORONA.  Continue the SAME hydrostatic integration upward, but
          // isothermally at T = bg_temp and starting from the last stellar level's own
          // pressure, so the join is pressure-matched and the only jump across it is in
          // temperature (density falls by T_star/T_corona).  RK2 in r on
          //   d ln p / dr = -rho(p, T_c) g(r) / p,
          // the same form `step` integrates, with the gradient fixed at zero instead of
          // taken from `nabla`: there is no flux to carry up here and no opacity to carry
          // it with (see rad_kappa_rmax), so the corona's only job is to stand still.
          Real lp = lnp(ijoin-1);
          for (int i = ijoin; i < nfine; ++i) {
            const Real r0 = rlo + (i-1)*drf;
            const Real p0 = exp(lp);
            const Real k1 = -DensFromPT(eos, rgas, p0, tbg)*GravAt(gm, r0)/p0;
            const Real pm = exp(lp + 0.5*drf*k1);
            const Real rm = r0 + 0.5*drf;
            const Real k2 = -DensFromPT(eos, rgas, pm, tbg)*GravAt(gm, rm)/pm;
            lp += drf*k2;
            lnp(i) = lp;
            tk(i) = tbg;
          }
        } else {
          for (int i = ijoin; i < nfine; ++i) {
            lnp(i) = log(pbg);
            tk(i) = tbg;
          }
        }
      }
      // diagnostics: tau from the top down, kappa, and the gradient actually used
      Real tsum = 0.0;
      for (int i = nfine-1; i >= 0; --i) {
        Real rho, kr, gr, vcd;
        nabla(rlo + i*drf, exp(lnp(i)), tk(i), tsum, rho, kr, gr, vcd);
        kap(i) = kr; grad(i) = gr; vcc(i) = vcd; rhoa(i) = rho;
        tau(i) = tsum;
        tsum += kr*rho*dun*drf*lun;
      }
      (void) tun;
    });
  }
  lnp_d_ = lnp; tk_d_ = tk; vc_d_ = vcc; rlo_ = rlo; drf_ = drf; nfine_ = nfine;
  tau_d_ = tau;            // problem/vpert_tau_lo/hi: the seed's optical-depth window

  // --- report the column, and dump it if asked
  {
    auto hlnp = Kokkos::create_mirror_view(lnp);
    auto htk = Kokkos::create_mirror_view(tk);
    auto hkap = Kokkos::create_mirror_view(kap);
    auto hgrad = Kokkos::create_mirror_view(grad);
    auto htau = Kokkos::create_mirror_view(tau);
    Kokkos::deep_copy(hlnp, lnp);
    Kokkos::deep_copy(htk, tk);
    Kokkos::deep_copy(hkap, kap);
    Kokkos::deep_copy(hgrad, grad);
    Kokkos::deep_copy(htau, tau);
    // where the star stops and the ambient medium begins.  The background block wrote
    // tk verbatim, so the first point from the top that is NOT at bg_temp is the join.
    bg_rjoin_ = 0.0;
    if (bg_rho_ > 0.0) {
      int ij = 0;
      for (int i = nfine-1; i >= 0; --i) {
        if (htk(i) != bg_temp_) { ij = i+1; break; }
      }
      bg_rjoin_ = rlo + ij*drf;
      if (bg_hydrostatic_) {
        // The corona is only useful if it can hold itself up on the grid it is given, so
        // print the join in full and the scale height it is supposed to have there.
        auto hrho = Kokkos::create_mirror_view(rhoa);
        Kokkos::deep_copy(hrho, rhoa);
        // mu from the corona's OWN state, not from <units>/mu: the gas is ionized at
        // 6e5 K and the table knows it, so read it back as rho k T/(p m_H).
        const Real rj = bg_rjoin_*lunit;
        const Real mu_c = hrho(ij)*dunit*kBoltz*htk(ij)
                          /(exp(hlnp(ij))*punit_*kMH);
        const Real hscale = kBoltz*bg_temp_*rj*rj/(mu_c*kMH*kGrav*mstar);
        std::cout << "red_giant: HYDROSTATIC corona at T = " << bg_temp_ << " K, join at "
                  << "r = " << rj << " cm, H = " << hscale << " cm = " << hscale/rj
                  << " r (mu = " << mu_c << ")" << std::endl;
        const int itopf = nfine-1;
        const int idx[4] = {ij-1, ij, ij+5, itopf};
        for (int q = 0; q < 4; ++q) {
          const int i = (idx[q] < nfine) ? idx[q] : nfine-1;
          std::cout << "red_giant:   i = " << i << " r = " << (rlo + i*drf)*lunit
                    << " p = " << exp(hlnp(i))*punit_ << " rho = " << hrho(i)*dunit
                    << " T = " << htk(i) << std::endl;
        }
      }
      std::cout << "red_giant: ambient medium rho = " << bg_rho_ << " g/cm^3, T = "
                << bg_temp_ << " K, from r = " << bg_rjoin_*lunit << " cm ("
                << bg_rjoin_*lunit/(rout*lunit) << " of x1max) outward; the star's "
                << "column is anchored at p = problem/ptop at r = " << ranch*lunit
                << " cm" << std::endl;
    }
    int iphot = nfine-1, iconv = -1;
    for (int i = nfine-1; i >= 0; --i) {
      if (htau(i) >= 2.0/3.0 && iphot == nfine-1) iphot = i;
    }
    // the outermost point (below the wall) where the column is on the adiabat, i.e.
    // where the radiative gradient EXCEEDS the one used (which is then grad_ad)
    for (int i = itop; i >= 0; --i) {
      Real p = exp(hlnp(i)), t = htk(i);
      Real grad_rad = 3.0*hkap(i)*(p*punit_)*lstar
                      /(16.0*M_PI*kArad*kClight*kGrav*mstar*t*t*t*t);
      if (grad_rad > hgrad(i)*(1.0 + 1.0e-6) && iconv < 0) {
        iconv = i;
      }
    }
    const int ibot = static_cast<int>((rin - rlo)/drf + 0.5);
    // THE OPACITY TABLE MUST BE DATA WHERE THIS STAR LIVES.  Outside the recorded logR
    // window the merge filled with the edge value, constant in density; a run there is
    // not using an opacity at all.  Checked on the ACTUAL radial grid including the
    // ghosts, since the boundary cells feed the wall flux -- not merely on the domain.
    // A first line of defence, not a guarantee: cells can still wander out as the
    // atmosphere relaxes.
    // ...and it is about the TABLE, so a run whose opacity is problem/kappa_const
    // never touches it: the window is then not a statement about anything.
    if (opac_lR_lo_ > -1.0e29 && !(kconst > 0.0)) {
      auto &x1v_chk = pmbp->pcoord->x1v;
      auto hx1v = Kokkos::create_mirror_view(x1v_chk);
      Kokkos::deep_copy(hx1v, x1v_chk);
      int nbad = 0;
      Real wl = 0.0, wh = 0.0, wr = 0.0;
      for (int m = 0; m <= nmb1; ++m) {
        for (int i = 0; i <= n1m1; ++i) {
          const Real rr_ = curv ? hx1v(m,i) : rin;
          if (!(rr_ > 0.0)) continue;
          // The AMBIENT MEDIUM is outside the tabulated stellar regime by construction
          // (logR ~ -11 for a cold, thin background), and is meant to be: at 1e-22 g/cm^3
          // its optical depth across the whole shell is ~1e-10 whatever kappa comes back,
          // so an edge-filled opacity there transports nothing.  The check below is about
          // the STAR, so skip the background rather than let it abort the run.
          if (bg_rho_ > 0.0 && rr_ >= bg_rjoin_) continue;
          Real lp, tt;
          ColumnAtHost(hlnp, htk, nfine, rlo, drf, rr_, lp, tt);
          const Real dd = DensFromPT(eos, rgas, exp(lp), tt)*dunit;
          if (!(dd > 0.0)) continue;
          const Real lR = log10(dd) - 3.0*log10(tt) + 18.0;
          // Only the DENSE side is a refusal now.  Below the window the table has been
          // zeroed above, which is the physically right answer for gas too thin to
          // absorb, so a cell there is legal -- it simply does not participate in the
          // radiative transfer.
          if (lR > opac_lR_hi_) {
            if (nbad == 0) { wl = lR; wh = lR; wr = rr_*lunit; }
            wl = fmin(wl, lR);
            wh = fmax(wh, lR);
            ++nbad;
          }
        }
      }
      if (nbad > 0) {
        std::cout << "### FATAL ERROR in red_giant: " << nbad << " cells (ghosts "
                  << "included) lie ABOVE the opacity table's valid window logR = ["
                  << opac_lR_lo_ << ", " << opac_lR_hi_ << "]: they reach logR " << wl
                  << " .. " << wh << ", first at r = " << wr << " cm. There the table is "
                  << "edge-filled, not data. Move the domain, or extend the tables."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::cout << "red_giant: no cell is denser than the opacity table's valid window "
                << "logR = [" << opac_lR_lo_ << ", " << opac_lR_hi_ << "]; cells below "
                << "it carry problem/opac_floor" << std::endl;
    }

    // --- the inner boundary, and the deep adiabat an OPEN one relaxes toward
    {
      const std::string ibc = pin->GetOrAddString("problem", "inner_bc", "wall");
      if (ibc.compare("open") == 0) {
        open_inner_ = true;
      } else if (ibc.compare("wall") != 0) {
        std::cout << "### FATAL ERROR in red_giant: problem/inner_bc must be wall or open"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      const std::string obc = pin->GetOrAddString("problem", "outer_bc", "wall");
      if (obc.compare("open") == 0) {
        open_outer_ = true;
      } else if (obc.compare("wall") != 0) {
        std::cout << "### FATAL ERROR in red_giant: problem/outer_bc must be wall or open"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      open_outer_noinflow_ = pin->GetOrAddBoolean("problem", "open_outer_noinflow", true);
      open_outer_lid_ = pin->GetOrAddBoolean("problem", "open_outer_lid", false);
      if (open_outer_ && global_variable::my_rank == 0) {
        std::cout << "red_giant: OPEN outer boundary -- the ghosts continue the top "
                  << "active cell hydrostatically at its own temperature"
                  << (open_outer_noinflow_ ? ", outflow only." : ", velocity copied.")
                  << " Mass and energy CAN leave the domain; problem/face_budget "
                  << "reports how much." << std::endl;
      }
      cs_change_ = pin->GetOrAddReal("problem", "s_relax_cs", 0.1);
      cp_change_ = pin->GetOrAddReal("problem", "s_relax_cp", 0.3);
      open_budget_ = pin->GetOrAddInteger("problem", "open_budget", 0);
      face_budget_ = pin->GetOrAddInteger("problem", "face_budget", 0);
      open_debug_ = pin->GetOrAddInteger("problem", "open_debug", 0);
      wall_noflux_ = pin->GetOrAddBoolean("problem", "wall_noflux", true);
      open_nomassflux_ = pin->GetOrAddBoolean("problem", "open_nomassflux", true);
      open_conserve_ = pin->GetOrAddBoolean("problem", "open_conserve", true);
      open_lstar_ = lstar;
      p_base_ = exp(hlnp(ibot))*punit_;
      t_base_ = htk(ibot);
      if (open_inner_) {
        if (pc != nullptr && pc->rad_flux_inner != 0.0) {
          std::cout << "### FATAL ERROR in red_giant: inner_bc = open carries the energy "
                    << "in as the ENTROPY of the inflow, so <block>/rad_flux_inner must "
                    << "be 0. It is " << pc->rad_flux_inner << ". With both, the "
                    << "luminosity is double counted." << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::cout << "red_giant: OPEN inner boundary (Stein-Nordlund / CO5BOLD). The "
                  << "deep adiabat is anchored at p = " << p_base_ << " dyn/cm^2, T = "
                  << t_base_ << " K; the emergent luminosity is now an OUTPUT."
                  << std::endl;
      }
    }
    // the top state must have a gas solution: under a general EOS with radiation the
    // total pressure cannot fall below a T^4/3, and SolveDensity then hands back the
    // table floor -- a column of floor density that looks like a run and is not one
    if (!useprof) {
      const Real p_t = exp(hlnp(itop)), t_t = htk(itop);
      const Real rho_t = DensFromPT(eos, rgas, p_t, t_t)*dunit;
      const Real prad = kArad*t_t*t_t*t_t*t_t/3.0;
      const bool eos_rad = eos.IsGeneral() && eos.tbl.radiation;
      if (!(rho_t > 1.0e-30) || !std::isfinite(rho_t)
          || (eos_rad && prad > 0.9*ptop_cgs)) {
        std::cout << "### FATAL ERROR in red_giant: no gas solution at the outer wall: "
                  << "problem/ptop = " << ptop_cgs << " dyn/cm^2 against a radiation "
                  << "pressure a T^4/3 = " << prad << " at T = " << t_t << " K (rho = "
                  << rho_t << " g/cm^3, eos_radiation "
                  << (eos_rad ? "on" : "off")
                  << "). Raise ptop above a T^4/3, or switch eos_radiation off."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: M = " << mstar << " g, L = " << lstar << " erg/s, Teff = "
                << teff << " K, rin = " << rin*lunit << " cm, rout = " << rout*lunit
                << " cm" << std::endl
                << "           wall p/T: outer " << ptop_cgs << " / " << htk(itop)
                << " K, inner " << exp(hlnp(ibot))*punit_ << " / " << htk(ibot) << " K"
                << std::endl
                << "           tau = 2/3 at r = " << (rlo + iphot*drf)*lunit
                << " cm, T = " << htk(iphot)
                << " K; radiative-convective boundary at r = "
                << ((iconv >= 0) ? (rlo + iconv*drf)*lunit : -1.0) << " cm"
                << std::endl;
      if (!dump.empty()) {
        std::ofstream df(dump);
        df.precision(10);
        df << std::scientific;
        auto hvc = Kokkos::create_mirror_view(vcc);
        Kokkos::deep_copy(hvc, vcc);
        df << "# red_giant initial column\n# r[cm] p[dyn/cm2] T[K] rho[g/cm3] kappa_R "
           << "grad tau v_mlt[cm/s]\n";
        for (int i = 0; i < nfine; ++i) {
          Real p = exp(hlnp(i)), t = htk(i);
          df << (rlo + i*drf)*lunit << " " << p*punit_ << " " << t << " "
             << DensFromPT(eos, rgas, p, t)*dunit << " " << hkap(i) << " " << hgrad(i)
             << " " << htau(i) << " " << hvc(i)*vunit << "\n";
        }
      }
    }
  }

  // --- the potential, needed on restarts too (rebuilt here, not stored)
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray4D<Real> ph2 = is_mhd ? pmbp->pmhd->phi0.x2f : pmbp->phydro->phi0.x2f;
  DvceArray4D<Real> ph3 = is_mhd ? pmbp->pmhd->phi0.x3f : pmbp->phydro->phi0.x3f;
  const bool have_phi = (etotgrav || wbdyn);
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  if (have_phi) {
    par_for("rg_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
      const Real xl = curv ? x1f_(m,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
      const Real xr = curv ? x1f_(m,i+1) : LeftEdgeX(i+1-is, indcs.nx1, x1lo, x1hi);
      const Real phi_c = PotAt(gm, rin, RadiusOf(curv, xc, rin, x1min));
      phicc(m,k,j,i) = phi_c;
      ph1(m,k,j,i) = PotAt(gm, rin, RadiusOf(curv, xl, rin, x1min));
      if (i == n1m1) ph1(m,k,j,i+1) = PotAt(gm, rin, RadiusOf(curv, xr, rin, x1min));
      ph2(m,k,j,i) = phi_c;
      ph3(m,k,j,i) = phi_c;
      if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
      if (k == n3m1) ph3(m,k+1,j,i) = phi_c;
    });
  }
  // --- the STATIC well-balanced background (<hydro>/wellbalance_static), needed on
  // restarts too: it is the problem generator's own initial column, evaluated over the
  // FULL arrays including every ghost zone, and it never evolves.
  //
  // WHY IT IS THE COLUMN.  The DYNAMIC scheme rebuilds a hydrostatic background from the
  // current state in every stencil, so the background responds to the perturbation it is
  // meant to be measured against; on this star that feedback drives an exponentially
  // growing buoyancy mode in the deep radiative zone (RG_fofc_long2/mode).  A background
  // fixed at the initial column cannot feed back, and it is a genuine hydrostatic
  // solution of this problem's own EOS and gravity, so the deviation stays small
  // wherever the star has not moved.
  //
  // WHAT MAKES IT BALANCED.  The deviation reconstruction hands the Riemann solver the
  // face background wherever the deviation vanishes; RemoveWbFlux then subtracts exactly
  // that face pressure from the momentum flux and the geometric source subtracts pwb, so
  // a cell whose state IS the background gets zero momentum change from the flux
  // divergence and the geometry together.  The gravity source below must therefore also
  // vanish there, which is what -g*(rho - rho_bg) does -- to round-off, for any
  // background field, however the background was built.
  const bool wbstat = is_mhd ? pmbp->pmhd->use_wellbalance_static
                             : pmbp->phydro->use_wellbalance_static;
  if (wbstat && wbdyn) {
    std::cout << "### FATAL ERROR in red_giant: wellbalance_static and "
              << "wellbalance_dynamic are two backgrounds for one scheme; pick one"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (wbstat) {
    const Real wbrmax = is_mhd ? pmbp->pmhd->wb_rmax : pmbp->phydro->wb_rmax;
    const Real wbrmin = is_mhd ? pmbp->pmhd->wb_rmin : pmbp->phydro->wb_rmin;
    const Real wbramp = wb_ramp_;
    auto u0wb = is_mhd ? pmbp->pmhd->u0wb : pmbp->phydro->u0wb;
    auto w0wb = is_mhd ? pmbp->pmhd->w0wb : pmbp->phydro->w0wb;
    auto wf1 = is_mhd ? pmbp->pmhd->w0facewb.x1f : pmbp->phydro->w0facewb.x1f;
    auto wf2 = is_mhd ? pmbp->pmhd->w0facewb.x2f : pmbp->phydro->w0facewb.x2f;
    auto wf3 = is_mhd ? pmbp->pmhd->w0facewb.x3f : pmbp->phydro->w0facewb.x3f;
    auto lnp_w = lnp_d_;
    auto tk_w = tk_d_;
    const int nfw = nfine_;
    const Real rlo_w = rlo_, drf_w = drf_;
    par_for("rg_wb_static", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
      const Real rc = RadiusOf(curv, xc, rin, x1min);
      Real dbg, ebg;
      WbBgAt(eos, rgas, igm1, lnp_w, tk_w, nfw, rlo_w, drf_w,
             wbrmin, wbrmax, wbramp, rc, dbg, ebg);
      w0wb(m,IDN,k,j,i) = dbg;
      w0wb(m,IVX,k,j,i) = 0.0;
      w0wb(m,IVY,k,j,i) = 0.0;
      w0wb(m,IVZ,k,j,i) = 0.0;
      w0wb(m,IEN,k,j,i) = ebg;
      u0wb(m,IDN,k,j,i) = dbg;
      u0wb(m,IM1,k,j,i) = 0.0;
      u0wb(m,IM2,k,j,i) = 0.0;
      u0wb(m,IM3,k,j,i) = 0.0;
      // with etotgrav the conserved energy carries the potential, exactly as the initial
      // state below does
      u0wb(m,IEN,k,j,i) = ebg + (etotgrav ? dbg*PotAt(gm, rin, rc) : 0.0);
      // the ANGULAR faces sit at the cell's own radius: the background is a function of
      // r alone, so their state is the cell's
      wf2(m,IDN,k,j,i) = dbg; wf2(m,IEN,k,j,i) = ebg;
      wf3(m,IDN,k,j,i) = dbg; wf3(m,IEN,k,j,i) = ebg;
      wf2(m,IVX,k,j,i) = 0.0; wf2(m,IVY,k,j,i) = 0.0; wf2(m,IVZ,k,j,i) = 0.0;
      wf3(m,IVX,k,j,i) = 0.0; wf3(m,IVY,k,j,i) = 0.0; wf3(m,IVZ,k,j,i) = 0.0;
      if (j == n2m1) {
        wf2(m,IDN,k,j+1,i) = dbg; wf2(m,IEN,k,j+1,i) = ebg;
        wf2(m,IVX,k,j+1,i) = 0.0; wf2(m,IVY,k,j+1,i) = 0.0; wf2(m,IVZ,k,j+1,i) = 0.0;
      }
      if (k == n3m1) {
        wf3(m,IDN,k+1,j,i) = dbg; wf3(m,IEN,k+1,j,i) = ebg;
        wf3(m,IVX,k+1,j,i) = 0.0; wf3(m,IVY,k+1,j,i) = 0.0; wf3(m,IVZ,k+1,j,i) = 0.0;
      }
      // the RADIAL faces take their own radius
      const Real xl = curv ? x1f_(m,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
      Real dfl, efl;
      WbBgAt(eos, rgas, igm1, lnp_w, tk_w, nfw, rlo_w, drf_w, wbrmin, wbrmax, wbramp,
             RadiusOf(curv, xl, rin, x1min), dfl, efl);
      wf1(m,IDN,k,j,i) = dfl; wf1(m,IEN,k,j,i) = efl;
      wf1(m,IVX,k,j,i) = 0.0; wf1(m,IVY,k,j,i) = 0.0; wf1(m,IVZ,k,j,i) = 0.0;
      if (i == n1m1) {
        const Real xr = curv ? x1f_(m,i+1) : LeftEdgeX(i+1-is, indcs.nx1, x1lo, x1hi);
        Real dfr, efr;
        WbBgAt(eos, rgas, igm1, lnp_w, tk_w, nfw, rlo_w, drf_w, wbrmin, wbrmax, wbramp,
               RadiusOf(curv, xr, rin, x1min), dfr, efr);
        wf1(m,IDN,k,j,i+1) = dfr; wf1(m,IEN,k,j,i+1) = efr;
        wf1(m,IVX,k,j,i+1) = 0.0; wf1(m,IVY,k,j,i+1) = 0.0; wf1(m,IVZ,k,j,i+1) = 0.0;
      }
    });
  }

  if (restart) return;

  // --- the initial state
  const int nf = nfine;
  auto vcd = vc_d_;
  const bool vpert_mlt = vpert_mlt_ && (mlt_alpha_ic > 0.0);
  const Real vmachmax = vpert_mach_max_;
  const Real vrmin = vpert_rmin_/lunit;             // code radius
  // problem/vpert_var and the tau window (see their declarations).  tau comes from the
  // INITIAL column on the same fine radial grid the state does, so the window is
  // evaluated with ColumnAt's interpolation and nothing extra is stored.
  const int pvar = vpert_var_;
  auto taud = tau_d_;
  const bool twin = (vpert_tau_hi_ > 0.0);
  const Real tlo = vpert_tau_lo_, thi = vpert_tau_hi_;
  const Real itlr = twin ? 1.0/log(thi/tlo) : 0.0;
  const bool vcart = vpert_cart_ && pmy_mesh_->use_cubed_sphere;
  auto &mbpanel_ic = pmbp->pmb->mb_panel;
  // The conserved momenta are COVARIANT on the cubed sphere's non-orthogonal angular
  // basis (GnomonicEquiangleLowerMom); c = 0 off that grid recovers the orthogonal form.
  auto &ccell_ic = pmbp->pcoord->cos_cell;
  const bool cs_ic = pmy_mesh_->use_cubed_sphere;
  // --- RANDOM, SYMMETRY-BREAKING PART OF THE CHART-FREE SEED --------------------------
  // problem/vpert_cart on its own leaves the flow with the EXACT 48-fold symmetry of the
  // cube: measured on V4, the stripe amplitudes along j and k agree to every printed
  // digit, every cold-cell count is a multiple of 12, and the eight cube-vertex copies
  // agree to 1e-6 out to t = 3.6e5.  A production run must not start on that symmetry,
  // so a random field is superposed on the cubic-harmonic one.
  //
  // It is a sum of vpert_nrand PLANE WAVES in the CARTESIAN position of the cell,
  // x = r qhat, and so is -- exactly like the cubic-harmonic part -- a pure function of
  // the physical point: no MeshBlock index, no panel index, no rank, and no per-cell
  // draw.  Per-cell white noise would be grid locked, and grid-scale noise is precisely
  // what the chart-free seed exists to avoid; a band-limited superposition is smooth on
  // the cell scale and resolved by construction.  Each wave carries a polarisation
  // PERPENDICULAR to its wave vector, e_n . k_n = 0, so the Cartesian field is
  // solenoidal (the spherical shell and the radial taper break that at the two radial
  // ends, as they already do for the structured part).
  //
  // Every draw comes from ONE host-side mt19937_64 with a fixed seed, run once here, so
  // all ranks build the same table and the seed is bit-reproducible under any
  // decomposition.  The generator is never called inside a kernel.
  const int vrnum = pin->GetOrAddInteger("problem", "vpert_nrand", 32);
  const int vrseed = pin->GetOrAddInteger("problem", "vpert_seed", 1234);
  const Real vrandamp = pin->GetOrAddReal("problem", "vpert_rand_amp", 1.0);
  // The default band is the two length scales the structured chart-free seed already
  // contains, so the random part adds no scale the seed did not have:
  //   vpert_lmin = its RADIAL wavelength.  That part is sin(3 . 2pi (x1 - x1lo)/
  //     (x1hi - x1lo)), i.e. three waves across a MeshBlock, so one third of a
  //     MeshBlock's radial extent (the mean one, if the radial grid is stretched).
  //   vpert_lmax = the ANGULAR wavelength of its l = 4 cubic harmonic at the top of the
  //     domain, 2 pi r_out / 4.
  // Both are given in CGS, like problem/vpert_rmin, and converted to code length here.
  const Real x1lo_g = pmy_mesh_->mesh_size.x1min, x1hi_g = pmy_mesh_->mesh_size.x1max;
  const int nmbx1 = pmy_mesh_->mesh_indcs.nx1/indcs.nx1;
  const Real vlmin_def = (x1hi_g - x1lo_g)/(3.0*static_cast<Real>(nmbx1));
  const Real vlmax_def = 0.5*M_PI*x1hi_g;
  const Real vlmin = pin->GetOrAddReal("problem","vpert_lmin",vlmin_def*lunit)/lunit;
  const Real vlmax = pin->GetOrAddReal("problem","vpert_lmax",vlmax_def*lunit)/lunit;
  // vpert_rand_amp = 0 switches the random part off and reproduces the structured
  // chart-free seed bit for bit: nothing below runs and the kernel branch is skipped.
  const int nrw = (vcart && vrandamp != 0.0 && vrnum > 0) ? vrnum : 0;
  if (nrw > 0 && (!(vlmin > 0.0) || vlmax < vlmin)) {
    std::cout << "### FATAL ERROR in red_giant: need 0 < problem/vpert_lmin <= "
              << "problem/vpert_lmax" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  DualArray2D<Real> vrwave("rg_vpert_rand", (nrw > 0) ? nrw : 1, 7);
  if (nrw > 0) {
    std::mt19937_64 rng(vrseed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double lnrat = log(static_cast<double>(vlmax/vlmin));
    for (int n=0; n<nrw; ++n) {
      // wave vector: direction uniform on the sphere, wavelength log-uniform in the band
      const double cz = 2.0*uni(rng) - 1.0;
      const double sz = sqrt((1.0 - cz*cz > 0.0) ? (1.0 - cz*cz) : 0.0);
      const double az = 2.0*M_PI*uni(rng);
      const double kh[3] = {sz*cos(az), sz*sin(az), cz};
      const double lam = static_cast<double>(vlmin)*exp(lnrat*uni(rng));
      const double kmag = 2.0*M_PI/lam;
      // an orthonormal pair spanning the plane perpendicular to kh, then a random
      // direction inside it: e . kh = 0 is what makes the superposition divergence free
      double t1[3];
      if (fabs(kh[2]) < 0.9) {
        t1[0] = -kh[1]; t1[1] = kh[0]; t1[2] = 0.0;
      } else {
        t1[0] = 0.0; t1[1] = -kh[2]; t1[2] = kh[1];
      }
      const double nt1 = sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
      for (int c=0; c<3; ++c) { t1[c] /= nt1; }
      const double t2[3] = {kh[1]*t1[2] - kh[2]*t1[1], kh[2]*t1[0] - kh[0]*t1[2],
                            kh[0]*t1[1] - kh[1]*t1[0]};
      const double psi = 2.0*M_PI*uni(rng);
      for (int c=0; c<3; ++c) {
        vrwave.h_view(n,c) = static_cast<Real>(kmag*kh[c]);
        vrwave.h_view(n,3+c) = static_cast<Real>(cos(psi)*t1[c] + sin(psi)*t2[c]);
      }
      vrwave.h_view(n,6) = static_cast<Real>(2.0*M_PI*uni(rng));
    }
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: chart-free seed + " << nrw << " random plane waves, "
                << "lambda in [" << vlmin*lunit << ", " << vlmax*lunit << "] cm, "
                << "seed " << vrseed << ", amp " << vrandamp << std::endl;
    }
  }
  vrwave.modify_host();
  vrwave.sync_device();
  auto vrwd = vrwave.d_view;
  // problem/vpert_sp_rand (off the cubed sphere, e.g. the spherical-polar WEDGE): a
  // random superposition of vpert_nrand modes PERIODIC over the mesh in x2 and x3,
  //     sum_n cos(2 pi (a_n u + b_n w) + phi_n) sin(pi c_n s + chi_n),
  // u, w, s the mesh-wide fractional coordinates, a_n, b_n integers in [-kmax, kmax] (not
  // both zero), c_n in 1..3.  Without it the seed off the cubed sphere is ONE mode whose
  // period is the MeshBlock, and everything else has to grow from round-off.  Radial
  // component only, so under vpert_var = eint it is a scalar entropy pattern.
  const int nsp = (!vcart && vrnum > 0) ? vrnum : 0;
  const bool sprand = pin->GetOrAddBoolean("problem", "vpert_sp_rand", false)
                      && (nsp > 0);
  const int spkmax = pin->GetOrAddInteger("problem", "vpert_sp_kmax", 8);
  DualArray2D<Real> spwave("rg_vpert_sp", (nsp > 0) ? nsp : 1, 5);
  if (sprand) {
    std::mt19937_64 rng2(vrseed + 7919);
    std::uniform_real_distribution<double> uni2(0.0, 1.0);
    for (int n=0; n<nsp; ++n) {
      int an = 0, bn = 0;
      while (an == 0 && bn == 0) {
        an = static_cast<int>(floor((2*spkmax + 1)*uni2(rng2))) - spkmax;
        bn = static_cast<int>(floor((2*spkmax + 1)*uni2(rng2))) - spkmax;
      }
      spwave.h_view(n,0) = static_cast<Real>(an);
      spwave.h_view(n,1) = static_cast<Real>(bn);
      spwave.h_view(n,2) = static_cast<Real>(1 + static_cast<int>(floor(3.0*uni2(rng2))));
      spwave.h_view(n,3) = static_cast<Real>(2.0*M_PI*uni2(rng2));
      spwave.h_view(n,4) = static_cast<Real>(2.0*M_PI*uni2(rng2));
    }
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: periodic random seed, " << nsp << " modes, |a|,|b| <= "
                << spkmax << ", seed " << vrseed << std::endl;
    }
  }
  spwave.modify_host();
  spwave.sync_device();
  auto spwd = spwave.d_view;
  const Real mx1lo = pmy_mesh_->mesh_size.x1min, mx1hi = pmy_mesh_->mesh_size.x1max;
  const Real mx2lo = pmy_mesh_->mesh_size.x2min, mx2hi = pmy_mesh_->mesh_size.x2max;
  const Real mx3lo = pmy_mesh_->mesh_size.x3min, mx3hi = pmy_mesh_->mesh_size.x3max;
  const Real spnorm = (nsp > 0) ? 2.0/sqrt(static_cast<Real>(nsp)) : 0.0;  // unit rms
  // rms matching, so that vpert_rand_amp = 1 gives the two parts the same rms speed.
  // The structured part is the single RADIAL component vpert cs sin(3 . 2pi u) G(qhat)
  // with G = 4(K4 - 1/5) + 6 K6.  Over the sphere <K4> = 1/5 (which is exactly what the
  // -1/5 removes, so G has zero mean), <K4^2> = 1/21, <K6^2> = 1/105 and <K4 K6> = 0 by
  // parity, giving <G^2> = 16(1/21 - 1/25) + 36/105 = 244/525; with <sin^2> = 1/2 its
  // rms speed is sqrt(244/525) vpert cs / sqrt(2).  For nrw waves of amplitude A with
  // random unit polarisation and phase the rms speed is A sqrt(nrw/2), so the two are
  // equal at A = sqrt(244/525) vpert cs vpert_rand_amp / sqrt(nrw) -- the sqrt(2)s
  // cancel.  (Both means verified numerically against 8e6 uniform directions.)
  const Real vrwamp = (nrw > 0)
                      ? sqrt(244.0/525.0)*vrandamp/sqrt(static_cast<Real>(nrw)) : 0.0;
  par_for("rg_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real x2lo = size.d_view(m).x2min, x2hi = size.d_view(m).x2max;
    const Real x3lo = size.d_view(m).x3min, x3hi = size.d_view(m).x3max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    Real lp, t;
    ColumnAt(lnp, tk, nf, rlo, drf, r, lp, t);
    const Real p = exp(lp);
    const Real d = DensFromPT(eos, rgas, p, t);
    const Real e = EintFromDensT(eos, rgas, igm1, d, t);
    Real v1 = 0.0, v2 = 0.0, v3 = 0.0;
    if (vpert > 0.0) {
      const Real x2v = CellCenterX(j-js, indcs.nx2, x2lo, x2hi);
      const Real x3v = CellCenterX(k-ks, indcs.nx3, x3lo, x3hi);
      // the amplitude: the MLT convective velocity where the column has one, else the
      // sound speed (the historical Mach-number seed)
      const Real csl = sqrt(gamma*p/d);
      // With the MLT seed, a cell the column gives no convective velocity -- the
      // radiative interior below the RCB, and the optically thin top -- gets NO seed.
      // Falling back to the sound speed there means seeding the densest material in the
      // star at Mach vpert, which with vpert = 1 is a transonic slug in the radiative
      // core: 8 orders of magnitude more kinetic energy than the rest of the seed put
      // together, and the upper domain went non-finite within 3400 cycles.
      Real amp = vpert_mlt ? 0.0 : csl;
      if (vpert_mlt) {
        const Real xf = (r - rlo)/drf;
        int i0 = static_cast<int>(xf);
        i0 = (i0 < 0) ? 0 : ((i0 > nf-2) ? nf-2 : i0);
        const Real f = xf - static_cast<Real>(i0);
        const Real vcl = vcd(i0)*(1.0 - f) + vcd(i0+1)*f;
        if (vcl > 0.0) amp = fmin(vcl, vmachmax*csl);
      }
      // the radial taper: 0 below vpert_rmin, a raised cosine over the five cells above
      if (vrmin > 0.0) {
        const Real wramp = 5.0*(x1hi - x1lo)/indcs.nx1;
        const Real srmp = (r - vrmin)/wramp;
        amp *= (srmp <= 0.0) ? 0.0
               : ((srmp >= 1.0) ? 1.0 : 0.5*(1.0 - cos(M_PI*srmp)));
      }
      // problem/vpert_tau_lo/hi: the OPTICAL-DEPTH window, the spherical reading of the
      // box's vpert_zlo/vpert_zhi.  sin(pi ln(tau_hi/tau)/ln(tau_hi/tau_lo)) inside,
      // zero outside, so the seed vanishes at both edges of the driving layer.
      if (twin) {
        const Real xf = (r - rlo)/drf;
        int i0 = static_cast<int>(xf);
        i0 = (i0 < 0) ? 0 : ((i0 > nf-2) ? nf-2 : i0);
        const Real ff = xf - static_cast<Real>(i0);
        const Real tac = taud(i0)*(1.0 - ff) + taud(i0+1)*ff;
        amp *= (tac > tlo && tac < thi) ? sin(M_PI*log(thi/tac)*itlr) : 0.0;
      }
      // problem/vpert_var = eint: the pattern below is then DIMENSIONLESS (the Mach
      // number of the velocity seed it replaces), and multiplies the internal energy at
      // fixed density instead of setting a velocity -- an ENTROPY seed.
      const Real cs = (pvar == 0) ? amp : amp/csl;
      const Real tp = 2.0*M_PI;
      if (vcart) {
        // CHART-INDEPENDENT seed: purely RADIAL, so it needs no tangential basis and is
        // the same physical field seen from any panel.  The angular pattern is built from
        // the two lowest cubic harmonics of the unit direction q, which are smooth on the
        // whole sphere and therefore continuous across every seam.
        Real q[3];
        cubed_sphere::PanelToCart(mbpanel_ic.d_view(m), 0.25*M_PI*x2v,
                                  0.25*M_PI*x3v, q);
        const Real q2x = q[0]*q[0], q2y = q[1]*q[1], q2z = q[2]*q[2];
        const Real k4 = q2x*q2y + q2y*q2z + q2z*q2x;         // l=4 cubic harmonic
        const Real k6 = q[0]*q[1]*q[2];                      // l=3 cubic harmonic
        v1 = vpert*cs*sin(3.0*tp*(xc - x1lo)/(x1hi - x1lo))
             *(4.0*(k4 - 0.2) + 6.0*k6);
        v2 = 0.0;
        v3 = 0.0;
        if (nrw > 0) {
          // the random plane waves, evaluated at the CARTESIAN position r qhat
          const Real px = r*q[0], py = r*q[1], pz = r*q[2];
          Real vr1 = 0.0, vr2 = 0.0, vr3 = 0.0;
          for (int n=0; n<nrw; ++n) {
            const Real cw = cos(vrwd(n,0)*px + vrwd(n,1)*py + vrwd(n,2)*pz + vrwd(n,6));
            vr1 += vrwd(n,3)*cw;
            vr2 += vrwd(n,4)*cw;
            vr3 += vrwd(n,5)*cw;
          }
          const Real aw = vrwamp*vpert*cs;   // the same envelope and taper as above
          vr1 *= aw; vr2 *= aw; vr3 *= aw;
          // Split the Cartesian vector on the basis the code stores velocities in:
          // V = v1 rhat + v2 e1 + v3 e2 with e1, e2 the UNIT gnomonic tangents, which
          // are perpendicular to rhat but NOT to each other (e1 . e2 = cos_cell).  So
          // the radial part is a plain dot product and the tangential pair has to be
          // raised with the 2x2 Gram matrix, exactly as cubed_sphere::TransformMomentum
          // does.  Taking the components straight off e1, e2 would be wrong by O(1) at
          // the panel corners.
          Real e1[3], e2[3];
          cubed_sphere::PanelTangents(mbpanel_ic.d_view(m), 0.25*M_PI*x2v,
                                      0.25*M_PI*x3v, e1, e2);
          const Real cg = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
          const Real dtg = 1.0 - cg*cg;
          const Real md2 = vr1*e1[0] + vr2*e1[1] + vr3*e1[2];
          const Real md3 = vr1*e2[0] + vr2*e2[1] + vr3*e2[2];
          v1 += vr1*q[0] + vr2*q[1] + vr3*q[2];
          v2 += (md2 - cg*md3)/dtg;
          v3 += (md3 - cg*md2)/dtg;
        }
      } else if (sprand) {
        const Real uu = (x2v - mx2lo)/(mx2hi - mx2lo), ww = (x3v - mx3lo)/(mx3hi - mx3lo);
        const Real ss = (xc - mx1lo)/(mx1hi - mx1lo);
        Real sum = 0.0;
        for (int n=0; n<nsp; ++n) {
          sum += cos(tp*(spwd(n,0)*uu + spwd(n,1)*ww) + spwd(n,3))
                *sin(M_PI*spwd(n,2)*ss + spwd(n,4));
        }
        v1 = vpert*cs*spnorm*sum;
        v2 = 0.0;
        v3 = 0.0;
      } else {
        v1 = vpert*cs*sin(3.0*tp*(xc - x1lo)/(x1hi - x1lo))
             *cos(2.0*tp*(x2v - x2lo)/(x2hi - x2lo))
             *cos(tp*(x3v - x3lo)/(x3hi - x3lo));
        v2 = vpert*cs*sin(2.0*tp*(x2v - x2lo)/(x2hi - x2lo))
             *cos(tp*(x3v - x3lo)/(x3hi - x3lo));
        v3 = vpert*cs*sin(tp*(x3v - x3lo)/(x3hi - x3lo))
             *cos(2.0*tp*(x2v - x2lo)/(x2hi - x2lo));
      }
    }
    // problem/vpert_var = eint: the RADIAL component of the pattern (which carries the
    // chart-free angular harmonics and, with vpert_cart, the random plane waves) becomes
    // a relative internal-energy perturbation at fixed density; the gas starts at rest.
    Real efac = 1.0;
    if (pvar == 1) {
      efac = 1.0 + v1;      // v1 already carries the vpert amplitude
      v1 = 0.0;
      v2 = 0.0;
      v3 = 0.0;
    }
    const Real cc = cs_ic ? ccell_ic(m,k,j) : 0.0;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*(v2 + cc*v3);
    u0(m,IM3,k,j,i) = d*(v3 + cc*v2);
    u0(m,IEN,k,j,i) = e*efac + 0.5*d*(v1*v1 + v2*v2 + v3*v3 + 2.0*cc*v2*v3);
    if (etotgrav) u0(m,IEN,k,j,i) += d*PotAt(gm, rin, r);
  });
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    par_for("rg_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      b.x1f(m,k,j,i) = 0.0;
      b.x2f(m,k,j,i) = 0.0;
      b.x3f(m,k,j,i) = 0.0;
      if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
      if (j == n2m1) b.x2f(m,k,j+1,i) = 0.0;
      if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantFaceBudget
//! \brief what the two radial faces of the domain carry, in erg/s and g/s against L.
//! The update adds +area(is) F(is) at the bottom and -area(ie+1) F(ie+1) at the top, so
//! both numbers below are signed as a GAIN by the domain.  Called from the source term,
//! which runs after the fluxes and before the next stage overwrites them.

void RedGiantFaceBudget(Mesh *pm) {
  if (face_budget_ <= 0 || pm->ncycle == face_cycle_) return;
  face_cycle_ = pm->ncycle;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &flx1 = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const Real dt = pm->dt;
  Real sEi = 0.0, sEo = 0.0, sMi = 0.0, sMo = 0.0;
  Kokkos::parallel_reduce("rg_facebud", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
  KOKKOS_LAMBDA(const int m, Real &aEi, Real &aEo, Real &aMi, Real &aMo) {
    const bool inb = (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user);
    const bool oub = (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        if (inb) {
          aEi += area1(m,k,j,is)*flx1(m,IEN,k,j,is);
          aMi += area1(m,k,j,is)*flx1(m,IDN,k,j,is);
        }
        if (oub) {
          aEo -= area1(m,k,j,ie+1)*flx1(m,IEN,k,j,ie+1);
          aMo -= area1(m,k,j,ie+1)*flx1(m,IDN,k,j,ie+1);
        }
      }
    }
  }, sEi, sEo, sMi, sMo);
#if MPI_PARALLEL_ENABLED
  {
    Real g[4] = {sEi, sEo, sMi, sMo};
    MPI_Allreduce(MPI_IN_PLACE, g, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    sEi = g[0]; sEo = g[1]; sMi = g[2]; sMo = g[3];
  }
#endif
  const Real eunit = pmbp->punit->pressure_cgs()*SQR(pmbp->punit->length_cgs())
                     *pmbp->punit->length_cgs();
  const Real munit = pmbp->punit->density_cgs()*SQR(pmbp->punit->length_cgs())
                     *pmbp->punit->length_cgs();
  face_E_in_  += sEi*dt*eunit;
  face_E_out_ += sEo*dt*eunit;
  face_M_in_  += sMi*dt*munit;
  face_M_out_ += sMo*dt*munit;
  // THE RADIATIVE LUMINOSITY.  The two-stream is a volumetric source, so none of it
  // shows up in the flux array above; a single dumped column cannot give L either, since
  // the emergent flux varies from column to column by more than the mean.  Sum the
  // solver's own net face flux instead: over the TOP face that IS L_out, and over each
  // column's cut face it is the flux the diffusion operator hands the two-stream from
  // below.  Fb is in code flux units on faces is..ie+1, so it needs the same area
  // weighting as the hydro fluxes.  The conduction module keeps no face flux of its own
  // (it accumulates straight into uflx), so its share of the top face is not separable
  // here -- it is already inside the "outer" number above.
  Real sRad = 0.0, sCut = 0.0;
  const bool rt_on = two_stream_rt::rt_face_flux_ready();
  if (rt_on) {
    auto Fb_g = two_stream_rt::rt_face_flux();
    auto icut_g = two_stream_rt::rt_cut_index();
    const int nblk = two_stream_rt::rt_face_nblk();
    Kokkos::parallel_reduce("rg_radbud", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &aR, Real &aC) {
      const bool oub = (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user);
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          if (oub) {
            Real ft = 0.0;
            for (int b=0; b<nblk; ++b) ft += Fb_g(m,b,ie+1,k,j);
            aR += area1(m,k,j,ie+1)*ft;
          }
          int ic = icut_g(m,k,j);
          if (ic < is) ic = is;
          if (ic <= ie) {
            Real fc = 0.0;
            for (int b=0; b<nblk; ++b) fc += Fb_g(m,b,ic,k,j);
            aC += area1(m,k,j,ic)*fc;
          }
        }
      }
    }, sRad, sCut);
#if MPI_PARALLEL_ENABLED
    {
      Real g2[2] = {sRad, sCut};
      MPI_Allreduce(MPI_IN_PLACE, g2, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      sRad = g2[0]; sCut = g2[1];
    }
#endif
  }
  const Real funit_rad = pmbp->punit->pressure_cgs()*pmbp->punit->velocity_cgs()
                         *SQR(pmbp->punit->length_cgs());
  if (pm->ncycle % face_budget_ == 0 && global_variable::my_rank == 0) {
    const Real tnow = pm->time*pmbp->punit->time_cgs();
    const Real dtw = tnow - face_t_last_;
    if (dtw > 0.0) {
      const Real iL = 1.0/open_lstar_;
      std::cout << "face budget: gain/L  inner=" << (face_E_in_ - face_E_in_l_)/dtw*iL
                << " outer=" << (face_E_out_ - face_E_out_l_)/dtw*iL
                << " | erg in=" << face_E_in_ << " out=" << face_E_out_
                << " | g in=" << face_M_in_ << " out=" << face_M_out_;
      if (rt_on) {
        std::cout << " | L_rad,out/L = " << sRad*funit_rad*iL
                  << " L_rad,cut/L = " << sCut*funit_rad*iL;
      }
      std::cout << " (t = " << tnow << " s)" << std::endl;
    }
    face_t_last_ = tnow;
    face_E_in_l_ = face_E_in_;
    face_E_out_l_ = face_E_out_;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantGravity
//! \brief the radial gravity source, -rho G M / r^2, or under wellbalance_dynamic the
//! background's own pressure difference across the cell (face-sum form on the
//! curvilinear grids, which carries the geometric term; plain difference on the column)

void RedGiantGravity(Mesh *pm, Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  RedGiantFaceBudget(pm);
  // problem/rt_profile_dt: the shell-averaged radial profile, once per cycle, on the
  // state this routine inherits (the start-of-cycle state in stage 1)
  RedGiantProfileTick(pm);
  // problem/nan_report call sites: one scan after each operator below that writes
  // u0(IEN).  The first is the state this routine INHERITS, i.e. everything RKUpdate
  // applied -- the hydro flux divergence and the conduction heat flux with it.
  RGNanScan(pm, "entry_after_RKUpdate+conduction");
  runaway_scan::Scan(pm, "entry_after_RKUpdate+conduction");
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;
  const bool wbx1 = is_mhd ? pmbp->pmhd->use_wb_x1 : pmbp->phydro->use_wb_x1;
  // the STATIC well-balanced background (filled from the initial column in UserProblem).
  // Its gravity source is -g*(rho - rho_bg): the background's own pressure difference has
  // already been taken out of the momentum RHS by RemoveWbFlux and the geometric source,
  // so a cell sitting exactly on the background must receive nothing here.  Where the
  // background was windowed to zero (outside [wb_rmin, wb_rmax]) this IS -rho g, so the
  // radius cutoffs need no branch of their own.
  const bool wbstat = is_mhd ? pmbp->pmhd->use_wellbalance_static
                             : pmbp->phydro->use_wellbalance_static;
  auto w0wb = is_mhd ? pmbp->pmhd->w0wb : pmbp->phydro->w0wb;
  // outside [wb_rmin, wb_rmax] the reconstruction dropped the well-balanced background,
  // so the source term must drop it too: the two are only well balanced TOGETHER.  Both
  // modules carry the pair, so take it from whichever one is evolving the fluid --
  // <mhd>/wb_r* under MHD, <hydro>/wb_r* otherwise.
  const Real wbrmax = is_mhd ? pmbp->pmhd->wb_rmax : pmbp->phydro->wb_rmax;
  const Real wbrmin = is_mhd ? pmbp->pmhd->wb_rmin : pmbp->phydro->wb_rmin;
  // <problem>/wb_grav_source = plain: keep the WB reconstruction, drop the WB source
  const bool wbplain = wb_grav_plain_;
  const WBOption wbo = is_mhd ? pmbp->pmhd->wb_option : pmbp->phydro->wb_option;
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray5D<Real> wbq0 = is_mhd ? pmbp->pmhd->wbq0 : pmbp->phydro->wbq0;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto &ccell_g = pmbp->pcoord->cos_cell;    // see CsKinetic
  const bool cs_g = pm->use_cubed_sphere;
  const bool curv = curv_;
  const Real gm = gm_, rin = rin_, x1min = x1min_;
  par_for("rg_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    const Real d = w0(m,IDN,k,j,i);
    const Real g = GravAt(gm, r);
    Real src = -bdt*g*d;
    if (!etotgrav) u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
    if (wbstat && !wbplain) {
      src = -bdt*g*(d - w0wb(m,IDN,k,j,i));
    }
    if (wbdyn && !wbplain && !(wbrmax > 0.0 && r > wbrmax) &&
        !(wbrmin > 0.0 && r < wbrmin)) {
      Real pl, pr, d1, d2, d3;
      if (wbx1) {
        WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      } else {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i), w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i), w0(m,IEN,k,j,i+1),
            phicc(m,k,j,i-1), ph1(m,k,j,i), phicc(m,k,j,i), ph1(m,k,j,i+1),
            phicc(m,k,j,i+1), d1, pl, d2, pr, d3);
      }
      if (curv) {
        // Same trap the sponge guards against (see "THE SPONGE RUNS BEFORE ConsToPrim"
        // below): the tabulated EOS takes log10(e), so a non-positive or non-finite
        // primitive energy returns NaN, src is NaN, and IM1 goes NaN while the density
        // stays finite -- the prod11 death at t=3.2106e7.  w0(IEN) can be non-positive
        // transiently here because the unclipped RT energy update precedes the floor.
        // Fall back to the background itself, p = (pl+pr)/2: with the cell pressure
        // taken as the background the curvature source reduces to the pure background
        // gradient, which is the well-balanced answer for a cell whose own state is not
        // yet physical.  (The non-curv branch below needs no p at all.)
        const Real e_ = w0(m,IEN,k,j,i);
        const Real p = (e_ > 0.0) ? eos.Pressure(d, e_) : 0.5*(pl + pr);
        src = bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))/volume(m,k,j,i);
      } else {
        src = bdt*(pr - pl)/((x1hi - x1lo)/indcs.nx1);
      }
    }
    u0(m,IM1,k,j,i) += src;
  });
  RGNanScan(pm, "gravity_WB_source");
  runaway_scan::Scan(pm, "gravity_WB_source");

  // --- the mixing-length convective flux (see the header): on each interior radial
  // face, where the face's d ln T / d ln p exceeds grad_ad,
  //     F_conv = rho c_p T sqrt(g delta) l^2 (grad - grad_ad)^(3/2)
  //              / (4 sqrt2 H_p^(3/2)),
  // l = alpha H_p, H_p = p/(rho g), delta = chi_T/chi_rho, c_p = c_v Gamma_1/chi_rho, all
  // in cgs then converted; zero on both walls.  Capped so one step moves at most a tenth
  // of the smaller neighbour's internal energy.
  if (mlt_alpha_ > 0.0) {
    auto fconv = fconv_;
    auto fdiag = fdiag_;
    auto taumlt = taumlt_;
    const Real alpha = mlt_alpha_, rgas = rgas_, igm1 = 1.0/gm1_, gamma = gm1_ + 1.0;
    const Real rmax_mlt = mlt_rmax_, lstar_c = lstar_cgs_;
    const bool mltmean = mlt_mean_;
    auto fmlt1d = fmlt1d_;
    const Real lgtlo = log10(mlt_tau_lo_), lgthi = log10(mlt_tau_hi_);
    const Real dunit = pmbp->punit->density_cgs();
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real vunit = pmbp->punit->velocity_cgs();
    const Real lunit = pmbp->punit->length_cgs();
    const Real tunit = pmbp->punit->time_cgs();
    auto &x1f_ = pmbp->pcoord->xx1f;
    // the run's OWN Rosseland optical depth, integrated DOWN each column from the top
    // active cell: tau(ie) = 0.5 kappa rho dx, tau(i) = tau(i+1) + half the cell above
    // plus half of this one.  The taper has to follow the gas the run actually has, not
    // the tau of the initial column: the photosphere moves.  320 cells serially per
    // (m,k,j) is nothing next to the flux kernel it feeds.
    if (!mltmean) {
      auto ktab = ktab_;
      auto klT = klT_;
      auto klD = klD_;
      const int knT = knT_, knD = knD_;
      const Real kfac = kfac_;
      par_for("rg_mlt_tau", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
        const Real dxu = (x1hi - x1lo)/indcs.nx1;
        Real prev = 0.0;
        for (int i = ie; i >= is; --i) {
          const Real d = w0(m,IDN,k,j,i);
          const Real e = w0(m,IEN,k,j,i);
          const Real p = eos.Pressure(d, e);
          const Real t = TempKelvin(eos, rgas, d, e, p);
          const Real kp = kfac*KappaTab(ktab, klT, klD, knT, knD, t, d*dunit);
          const Real dxc = (curv ? (x1f_(m,i+1) - x1f_(m,i)) : dxu)*lunit;
          const Real dtau = 0.5*kp*(d*dunit)*dxc;
          taumlt(m,k,j,i) = prev + dtau;
          prev = taumlt(m,k,j,i) + dtau;
        }
      });
    }
    // --- THE SHELL-MEAN CLOSURE (problem/mlt_mean; see the declaration for the two
    // measured defects it replaces).  Pass 1: horizontal sums of rho, e, T and p per
    // radial cell over every column of every rank.  T and p are evaluated PER CELL and
    // then averaged, not taken from the averaged state, so the mean pressure is the one
    // the gas actually has.  Every MeshBlock spans the whole radius (checked at
    // start-up), so the local index i IS the global radial index.
    if (mltmean) {
      auto shell = shell_;
      // THE HAND-OVER HAS TO USE THE OPERATOR'S OWN FLUX, not an estimate of it.  The
      // closure asks convection for whatever radiation leaves undone, so its idea of
      // "what radiation carries" must be the number the conduction module will actually
      // put on that face this stage.  16 sigma T^3/(3 kappa rho) |dT/dr| on the mean
      // profile is not that number: it has neither the free-streaming limiter nor the
      // tau-blend weight, and at the RCB it read 0.99 L where the operator delivered
      // 0.969 L, so 1-3 % of L was deposited there every step with nothing to carry it
      // and the superadiabatic dipole at i = 5/6 grew linearly.  Conduction stores no
      // face flux (it accumulates straight into uflx), so replicate it exactly, per
      // column, from the same w0 the flux task used -- the pgen source runs in the same
      // stage, after `flux` and before `w0` is touched again, so these are the identical
      // inputs -- and shell-average the RESULT rather than the profile.
      Conduction *pcnd = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
      const bool condrad = (pcnd != nullptr &&
                            pcnd->iso_cond_type.compare("radiative") == 0);
      const bool cgen = eos.IsGeneral();
      auto wtemp_c = is_mhd ? pmbp->pmhd->wtemp : pmbp->phydro->wtemp;
      auto wder_c = is_mhd ? pmbp->pmhd->wder : pmbp->phydro->wder;
      const Real tunitK = pmbp->punit->temperature_cgs();
      const Real cgm1 = gm1_;
      const bool cktab = condrad && pcnd->rad_kappa_tab && pcnd->rad_kr_nT > 0;
      const bool ckrho = cktab && pcnd->rad_kappa_rho;
      auto ckt = cktab ? pcnd->rad_kr_tab : DvceArray2D<Real>("d",1,1);
      auto cklT = cktab ? pcnd->rad_kr_lT : DvceArray1D<Real>("d",1);
      auto cklP = cktab ? pcnd->rad_kr_lP : DvceArray1D<Real>("d",1);
      const int cknT = cktab ? pcnd->rad_kr_nT : 0;
      const int cknP = cktab ? pcnd->rad_kr_nP : 0;
      const Real ckfac = condrad ? pcnd->rad_kappa_fac : 1.0;
      const Real cmet = condrad ? pcnd->rad_met : 0.0;
      const bool climit = condrad && pcnd->rad_flux_limit;
      const bool ctaum = condrad && pcnd->rad_tau_mode;
      const bool cbldr = condrad && pcnd->rad_blend_radial;
      const Real cpcut = (condrad && pcnd->rad_tau_mode) ? -1.0
                         : (condrad ? pcnd->rad_pcut : 0.0);
      auto cwf = (condrad && pcnd->rad_tau_mode) ? pcnd->rad_w
                 : DvceArray4D<Real>("d",1,1,1,1);
      const Real csig = 5.670374419e-5;
      const bool rank0 = (global_variable::my_rank == 0);
      // the TWO-STREAM's own net face flux, the third carrier.  Above the tau blend the
      // conduction weight is zero and the two-stream carries the whole luminosity; a
      // deficit that counted only conduction therefore asked the subgrid model for a
      // second full L in optically thin gas, and the top cell went to Mach ~1.
      const bool rt2s = two_stream_rt::rt_face_flux_ready();
      auto Fb2s = rt2s ? two_stream_rt::rt_face_flux()
                       : DvceArray5D<Real>("d",1,1,1,1,1);
      const int n2sb = rt2s ? two_stream_rt::rt_face_nblk() : 0;
      par_for("rg_mlt_shell", DevExeSpace(), is, ie, KOKKOS_LAMBDA(const int i) {
        Real sd = 0.0, se = 0.0, st = 0.0, sp = 0.0, sn = 0.0;
        Real sfc = 0.0, sw = 0.0, fc0 = 0.0;
        // the RESOLVED convective flux needs a CORRELATION, so the products have to be
        // summed per cell and decorrelated afterwards, not formed from the means
        Real srv = 0.0, srvh = 0.0, sh = 0.0, ske = 0.0, s2s = 0.0;
        for (int m = 0; m <= nmb1; ++m) {
          for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
              const Real d = w0(m,IDN,k,j,i);
              const Real e = w0(m,IEN,k,j,i);
              const Real p = eos.Pressure(d, e);
              sd += d;
              se += e;
              sp += p;
              st += TempKelvin(eos, rgas, d, e, p);
              sn += 1.0;
              const Real vr = w0(m,IVX,k,j,i);
              const Real v2 = SQR(vr) + SQR(w0(m,IVY,k,j,i)) + SQR(w0(m,IVZ,k,j,i));
              const Real hh = (e + p)/d;                 // specific enthalpy
              if (rt2s) {
                for (int b = 0; b < n2sb; ++b) s2s += Fb2s(m,b,i,k,j);
              }
              srv += d*vr;
              srvh += d*vr*hh;
              sh += hh;
              ske += 0.5*d*vr*v2;
              if (!condrad || i == is) continue;
              // --- conduction's radiative x1-face flux, term for term (see
              // Conduction::AddIsotropicHeatFluxRadiative) ---
              const Real tl = cgen ? wtemp_c(m,k,j,i-1)
                                   : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*cgm1;
              const Real tr = cgen ? wtemp_c(m,k,j,i)
                                   : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*cgm1;
              const Real pl = cgen ? wder_c(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*cgm1;
              const Real pr = cgen ? wder_c(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*cgm1;
              const Real dlc = curv ? (x1v_(m,i) - x1v_(m,i-1))
                                    : (size.d_view(m).x1max
                                       - size.d_view(m).x1min)/indcs.nx1;
              const Real wt = (ctaum && cbldr) ? cwf(m,k,j,i) : 1.0;
              const Real pfc = 0.5*(pl + pr);
              Real fcond = 0.0;
              if (pfc >= cpcut) {
                const Real tkf = 0.5*(tl + tr)*tunitK;
                const Real rhof = 0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dunit;
                const Real kr = cktab
                    ? RosselandTable(ckt, cklT, cklP, cknT, cknP, tkf,
                                     ckrho ? rhof : pfc*punit_)
                    : RosselandFreedman2014(tkf, pfc*punit_, cmet);
                const Real kapr = 16.0*csig*tkf*tkf*tkf/(3.0*ckfac*kr*rhof);
                fcond = -kapr*((tr - tl)/dlc)*tunitK/lunit;
                if (climit) {
                  const Real ffree = csig*tkf*tkf*tkf*tkf;
                  fcond /= sqrt(1.0 + SQR(fcond/ffree));
                }
                fcond *= wt;
              }
              sfc += fcond;                       // cgs already: eflx_unit cancels below
              sw += wt;
              if (rank0 && m == 0 && k == ks && j == js) fc0 = fcond;
            }
          }
        }
        shell(i,0) = sd; shell(i,1) = se; shell(i,2) = st; shell(i,3) = sp;
        shell(i,4) = sn; shell(i,5) = sfc; shell(i,6) = sw; shell(i,7) = fc0;
        shell(i,8) = srv; shell(i,9) = srvh; shell(i,10) = sh; shell(i,11) = ske;
        shell(i,12) = s2s;
      });
#if MPI_PARALLEL_ENABLED
      {
        const int ni = ie - is + 1;
        auto hs = Kokkos::create_mirror_view(shell);
        Kokkos::deep_copy(hs, shell);
        MPI_Allreduce(MPI_IN_PLACE, &hs(is,0), 13*ni, MPI_ATHENA_REAL, MPI_SUM,
                      MPI_COMM_WORLD);
        Kokkos::deep_copy(shell, hs);
      }
#endif
      // Pass 2, CLOSURE v3.  The star has to carry F_req = L/(4 pi r^2) through every
      // shell, and three things can carry it: the conduction operator, the resolved
      // flow, and this subgrid model.  v2 asked for f = min(F_MLT(x), F_req - F_rad)
      // and got both halves wrong -- it double counted the resolved flux (the run
      // radiated 1.0-1.5 L against the control's 0.3-0.6 L) and it gated the whole
      // thing on the sign of a discrete x that is round-off on the adiabat (14 % of L
      // uncarried at the RCB face, i = 7, where x came out at -3e-7).  v3 decides
      // instability with grad_rad > grad_ad, which differentiates nothing, and asks for
      // the DEFICIT the other two carriers leave.  No tau taper and no mlt_rmax here:
      // the deficit IS the hand-over and it vanishes wherever radiation takes the load.
      auto fmean = fmean_;
      const Real xthr = mlt_x_thr_, relaxt = mlt_relax_;
      const bool ffix = mlt_flux_fix_;
      const Real ffcap = mlt_flux_fix_cap_, ffrmax = mlt_flux_fix_rmax_;
      const bool seed = !mlt_relax_seeded_;
      auto ktab = ktab_;
      auto klT = klT_;
      auto klD = klD_;
      const int knT = knT_, knD = knD_;
      const Real kfac = kfac_;
      const Real sigsb = 0.25*kArad*kClight;
      par_for("rg_mlt_mean_face", DevExeSpace(), is+1, ie, KOKKOS_LAMBDA(const int i) {
        // fmlt1d is NOT cleared here: it is the relaxed state and carries over.
        for (int q = 0; q < 17; ++q) fmean(i,q) = 0.0;
        const Real nl = shell(i-1,4), nr = shell(i,4);
        if (!(nl > 0.0) || !(nr > 0.0)) return;
        const Real dl = shell(i-1,0)/nl, dr_ = shell(i,0)/nr;
        const Real el = shell(i-1,1)/nl, er = shell(i,1)/nr;
        const Real tl = shell(i-1,2)/nl, tr = shell(i,2)/nr;
        const Real pl = shell(i-1,3)/nl, pr = shell(i,3)/nr;
        const Real x1lo = size.d_view(0).x1min, x1hi = size.d_view(0).x1max;
        const Real dxu = (x1hi - x1lo)/indcs.nx1;
        const Real xf = curv ? x1f_(0,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
        const Real rf = RadiusOf(curv, xf, rin, x1min);
        const Real rcm = rf*lunit;
        const Real freq = lstar_c/(4.0*M_PI*rcm*rcm);
        const Real pf = 0.5*(pl + pr), tf = 0.5*(tl + tr), df = 0.5*(dl + dr_);
        const Real ef = 0.5*(el + er);
        const Real dxcm = (curv ? (x1v_(0,i) - x1v_(0,i-1)) : dxu)*lunit;
        // what radiation carries.  fdiff is the plain diffusion estimate, kept only as
        // a diagnostic; frad is the CONDUCTION OPERATOR'S OWN discrete face flux,
        // shell-averaged, which is the flux that will actually be deposited.
        const Real kp = kfac*KappaTab(ktab, klT, klD, knT, knD, tf, df*dunit);
        Real fdiff = 0.0;
        if (kp > 0.0 && dxcm > 0.0) {
          fdiff = 16.0*sigsb*tf*tf*tf/(3.0*kp*df*dunit)*fabs(tl - tr)/dxcm;
        }
        const Real frad = condrad ? shell(i,5)/nr : fdiff;
        // THE RESOLVED CONVECTIVE FLUX.  The star is a 3-D run: most of the convective
        // luminosity below the photosphere is CARRIED BY THE FLOW ALREADY, and a closure
        // that asks the subgrid model for F_req - F_rad hands it a second copy of the
        // same energy.  Measured: with that closure the run radiated 1.0-1.5 L from the
        // top face against the 0.3-0.6 L of the no-MLT control.  So subtract what the
        // resolved flow carries -- the enthalpy correlation plus the kinetic term,
        //   F_res = <rho v_r h> - <h><rho v_r> + <(1/2) rho v_r v^2>,
        // formed per cell (the products were summed, not the means) and averaged to the
        // face.  Removing <h><rho v_r> is what makes this a FLUX and not a mass flux
        // times an enthalpy: the shell carries no net mass, but round-off in <rho v_r>
        // times a large <h> would otherwise swamp the correlation.
        Real fres = 0.0;
        {
          const Real rvl = shell(i-1,8)/nl, rvr_ = shell(i,8)/nr;
          const Real rhl = shell(i-1,9)/nl, rhr = shell(i,9)/nr;
          const Real hl = shell(i-1,10)/nl, hr = shell(i,10)/nr;
          const Real kl = shell(i-1,11)/nl, kr_ = shell(i,11)/nr;
          const Real fl = rhl - hl*rvl + kl;
          const Real fr = rhr - hr*rvr_ + kr_;
          fres = 0.5*(fl + fr)*punit_*vunit;                           // cgs
        }
        fmean(i,0) = rcm; fmean(i,1) = tf; fmean(i,2) = pf*punit_;
        fmean(i,7) = frad; fmean(i,9) = freq;
        fmean(i,10) = fdiff; fmean(i,11) = condrad ? shell(i,6)/nr : 0.0;
        fmean(i,12) = condrad ? shell(i,7) : 0.0;
        // the two-stream's share of this face, shell-averaged.  It is written only at
        // faces at or above each column's cut and is exactly zero below, so it can be
        // summed as it stands.
        const Real f2s = shell(i,12)/nr*punit_*vunit;                  // cgs
        // the blend weight of this face: the share of the radiative flux the CONDUCTION
        // operator already carries, and therefore the share of F_2s that frad repeats
        const Real wbl = condrad ? shell(i,6)/nr : 0.0;
        fmean(i,13) = fres; fmean(i,16) = f2s;
        // thermodynamics at the face, cgs.  Needed by both branches below, so it is no
        // longer behind the x > 0 test the old closure returned on.
        const Real g = GravAt(gm, rf)*lunit/(tunit*tunit);
        const Real hp = pf*punit_/(df*dunit*g);
        const Real ell = alpha*hp;
        Real cv = (EintFromDensT(eos, rgas, igm1, df, 1.01*tf)
                   - EintFromDensT(eos, rgas, igm1, df, tf))/(0.01*tf)*punit_/(df*dunit);
        if (!(cv > 0.0)) cv = 1.5*pf*punit_/(df*dunit*tf);
        Real delta = 1.0, cp = gamma*cv;
        if (eos.IsGeneral()) {
          const Real chit = eos.ChiT(df, ef), chir = eos.ChiRho(df, ef);
          delta = chit/chir;
          cp = cv*eos.Gamma1(df, ef)/chir;
        }
        const Real grad_ad = GradAd(eos, gamma, rgas, pf, tf);
        Real grad = 0.0, x = -1.0;
        if (pl > pr && tl > tr) {
          grad = log(tl/tr)/log(pl/pr);
          x = grad - grad_ad;
        }
        fmean(i,3) = grad; fmean(i,4) = grad_ad; fmean(i,5) = x;
        // (a) IS THE SHELL CONVECTIVE?  Not "is the discrete x positive" -- that
        // question is unanswerable on an adiabat resolved to 1e-6 -- but Schwarzschild
        // in the form that does not differentiate the solution at all:
        //   grad_rad = 3 kappa p F_req / (16 sigma g T^4)  >  grad_ad.
        // kappa here is the same Rosseland table lookup the conduction operator makes,
        // at the same face state, so the two criteria cannot disagree about the opacity.
        const Real grad_rad = 3.0*kp*(pf*punit_)*freq
                              /(16.0*sigsb*g*tf*tf*tf*tf);
        fmean(i,14) = grad_rad;
        Real ftar = 0.0;                                               // cgs
        if (grad_rad > grad_ad) {
          // (b) what is left for the SUBGRID model, after the diffusion operator, the
          // two-stream and the resolved flow have each been counted once
          // THE DOUBLE COUNT.  frad is the conduction module's BLENDED face flux,
          // w F_raddiff, while f2s is the two-stream's UNBLENDED one; the apply kernel
          // hands the shell (1-w) F_2s, not F_2s.  Subtracting the whole F_2s therefore
          // counts the blended part twice, the deficit comes out too small wherever
          // w > 0, and the subgrid model is asked for a flux the star is already
          // carrying (tests_3d/handover/README.md section 3).  Exactly f2s when w = 0,
          // which is the shipped whole-column blend, so that configuration is bitwise.
          Real f2se = f2s;
          if (wbl > 0.0) f2se = (1.0 - wbl)*f2s;
          const Real dfc = fmax(0.0, freq - frad - fres - f2se);
          fmean(i,15) = dfc;
          if (dfc > 0.0) {
            // (c) efficient vs inefficient.  Where x is a real number the mixing-length
            // amplitude is the physics and it binds; where it is jitter, only the
            // deficit is known, and the deficit is what gets carried.  D goes to zero
            // continuously at the RCB because F_rad_cond goes to F_req there, so this
            // does not create a jump at the boundary the way the sign test did.
            Real f = dfc;
            if (x > xthr) {
              const Real fmlt = df*dunit*cp*tf*sqrt(g*delta)*ell*ell*x*sqrt(x)
                                /(4.0*sqrt(2.0)*hp*sqrt(hp));
              fmean(i,6) = fmlt;
              f = fmin(fmlt, dfc);
            }
            // NO chi LIMITER HERE.  The per-column path applies one because there the
            // flux really is a nonlinear diffusion of the cell's own gradient, with
            // chi = 1.5 (F/x) H_p/(T rho c_p) unbounded as x -> 0.  In the shell-mean
            // path the applied flux is an IMPOSED 1-D profile: it is set from the shell
            // budget, capped at F_req, ramped in over mlt_ramp_time and relaxed over
            // mlt_relax_time, and it does not respond to the local gradient at all --
            // so there is no explicit-diffusion eigenvalue to bound.  Applying the cap
            // anyway, at x_eff = x_thr, throttled the flux to 0.16 F_req over
            // i = 160-300 while nothing else carried it, and left 65-82 % of the
            // luminosity uncarried through the upper envelope.
            ftar = fmax(0.0, f);
          }
        }
        // problem/mlt_flux_fix: the signed residual where the closure carried nothing
        if (ffix && ftar == 0.0 && (ffrmax <= 0.0 || rcm < ffrmax)) {
          const Real f2sf = (wbl > 0.0) ? (1.0 - wbl)*f2s : f2s;
          const Real res = freq - frad - fres - f2sf;
          ftar = fmax(-ffcap*freq, fmin(ffcap*freq, res));
        }
        // (d) relax the applied profile toward the target
        const Real fnew = ftar/(punit_*vunit);                         // code flux
        if (relaxt > 0.0 && !seed) {
          fmlt1d(i) += fmin(1.0, bdt*tunit/relaxt)*(fnew - fmlt1d(i));
        } else {
          fmlt1d(i) = fnew;
        }
        fmean(i,8) = fmlt1d(i)*punit_*vunit;
      });
      mlt_relax_seeded_ = true;
    }
    // problem/mlt_ramp_time: fade the applied flux in, so the cut region heats over
    // that time instead of in one step (see the declaration).
    Real ramp = (mlt_ramp_ > 0.0)
                ? fmin(1.0, pm->time*tunit/mlt_ramp_) : 1.0;
    // problem/mlt_ramp_down_time / mlt_hold_time: the hand-over to resolved convection
    // (see the declaration).  Off by default, so nothing above changes.
    if (mlt_ramp_dn_ > 0.0) {
      const Real tdn = pm->time*tunit - mlt_hold_;
      ramp *= (tdn <= 0.0) ? 1.0
              : ((tdn >= mlt_ramp_dn_) ? 0.0 : (1.0 - tdn/mlt_ramp_dn_));
    }
    if (global_variable::my_rank == 0 && mlt_ramp_dn_ > 0.0
        && pm->ncycle % 100 == 0) {
      std::cout << "### red_giant: MLT subgrid flux fraction = " << ramp
                << " at t = " << pm->time*tunit << " s" << std::endl;
    }
    par_for("rg_mlt_flux", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fconv(m,k,j,i) = 0.0;
      if (m == 0 && k == ks && j == js) {
        for (int q = 0; q < 8; ++q) fdiag(i,q) = 0.0;
      }
      if (i == is || i == ie+1) return;
      // the shell-mean path already holds the whole closure in fmlt1d(i); all that is
      // left per column is the local safety cap (a tenth of the smaller neighbour's
      // internal energy per step), so a cold column cannot be over-drained
      if (mltmean) {
        const Real x1lo_ = size.d_view(m).x1min, x1hi_ = size.d_view(m).x1max;
        const Real dxl = curv ? (x1f_(m,i) - x1f_(m,i-1))
                              : (x1hi_ - x1lo_)/indcs.nx1;
        const Real fmax_ = 0.1*fmin(w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i))*dxl/bdt;
        fconv(m,k,j,i) = ramp*fmin(fmlt1d(i), fmax_);
        return;
      }
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xf = curv ? x1f_(m,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
      const Real rf = RadiusOf(curv, xf, rin, x1min);
      if (rmax_mlt > 0.0 && rf*lunit > rmax_mlt) return;    // radiation's layers
      const Real dl = w0(m,IDN,k,j,i-1), dr_ = w0(m,IDN,k,j,i);
      const Real el = w0(m,IEN,k,j,i-1), er = w0(m,IEN,k,j,i);
      const Real pl = eos.Pressure(dl, el), pr = eos.Pressure(dr_, er);
      const Real tl = TempKelvin(eos, rgas, dl, el, pl);
      const Real tr = TempKelvin(eos, rgas, dr_, er, pr);
      if (!(pl > pr) || !(tl > tr)) return;             // not stratified the right way
      const Real grad = log(tl/tr)/log(pl/pr);
      const Real pf = 0.5*(pl + pr), tf = 0.5*(tl + tr), df = 0.5*(dl + dr_);
      const Real grad_ad = GradAd(eos, gamma, rgas, pf, tf);
      if (!(grad > grad_ad)) return;
      const Real ef = 0.5*(el + er);
      // thermodynamics at the face, cgs
      const Real g = GravAt(gm, rf)*lunit/(tunit*tunit);
      const Real hp = pf*punit_/(df*dunit*g);
      const Real ell = alpha*hp;
      Real cv = (EintFromDensT(eos, rgas, igm1, df, 1.01*tf)
                 - EintFromDensT(eos, rgas, igm1, df, tf))/(0.01*tf)*punit_/(df*dunit);
      if (!(cv > 0.0)) cv = 1.5*pf*punit_/(df*dunit*tf);
      Real delta = 1.0, cp = gamma*cv;
      if (eos.IsGeneral()) {
        const Real chit = eos.ChiT(df, ef), chir = eos.ChiRho(df, ef);
        delta = chit/chir;
        cp = cv*eos.Gamma1(df, ef)/chir;
      }
      const Real x = grad - grad_ad;
      const Real vc = sqrt(g*delta*ell*ell*x/(8.0*hp));
      Real f = df*dunit*cp*tf*sqrt(g*delta)*ell*ell*x*sqrt(x)
               /(4.0*sqrt(2.0)*hp*sqrt(hp));
      // the subgrid flux can never exceed what the star has to carry.  Without this a
      // layer that steepened while nothing carried its flux (sup ~ 0.5 below the
      // photosphere) hands MLT a flux of many L the moment it is switched on.
      const Real rcm = rf*lunit;
      f = fmin(f, lstar_c/(4.0*M_PI*rcm*rcm));
      // the hand-over to radiation, on the face's own optical depth (geometric mean of
      // the two cells, which is the linear mean in log tau the weight is built on)
      const Real taul = taumlt(m,k,j,i-1), taur = taumlt(m,k,j,i);
      if (!(taul > 0.0) || !(taur > 0.0)) return;
      const Real lgtf = 0.5*(log10(taul) + log10(taur));
      Real wtau = (lgthi > lgtlo) ? (lgtf - lgtlo)/(lgthi - lgtlo)
                                  : ((lgtf >= lgthi) ? 1.0 : 0.0);
      wtau = (wtau < 0.0) ? 0.0 : ((wtau > 1.0) ? 1.0 : wtau);
      if (!(wtau > 0.0)) return;
      f *= wtau;
      f /= (punit_*vunit);                                // code flux
      // STABILITY.  This is an explicit nonlinear diffusion of T with diffusivity
      //   chi = (dF/d(dT/dr)) / (rho c_p) = (3/2) (F/x) (H_p/T) / (rho c_p),
      // which is unbounded as the superadiabaticity x grows and rho falls (the
      // tabulated column blew up 1e6-fold this way).  An explicit step is stable for
      // chi <= dx^2/(2 dt); the flux is scaled down to hold chi at half that, so a
      // strongly superadiabatic layer under-transports rather than explodes.  A tenth
      // of the smaller internal energy per step caps it as well.
      const Real dxl = curv ? (x1f_(m,i) - x1f_(m,i-1)) : (x1hi - x1lo)/indcs.nx1;
      const Real dxc = dxl*lunit;
      const Real chi = 1.5*(f*punit_*vunit/x)*(hp/tf)/(df*dunit*cp);   // cm^2/s
      const Real chimax = 0.25*dxc*dxc/(bdt*tunit);
      if (chi > chimax) f *= chimax/chi;
      const Real fmax_ = 0.1*fmin(el, er)*dxl/bdt;
      fconv(m,k,j,i) = ramp*fmin(f, fmax_);
      if (m == 0 && k == ks && j == js) {
        fdiag(i,0) = grad; fdiag(i,1) = grad_ad; fdiag(i,2) = hp; fdiag(i,3) = cp;
        fdiag(i,4) = vc; fdiag(i,5) = f*punit_*vunit; fdiag(i,6) = fmax_*punit_*vunit;
        fdiag(i,7) = fconv(m,k,j,i)*punit_*vunit;
      }
    });
    ++mlt_dump_calls_;
    if (!mlt_dumped_ && !mlt_dump_.empty() && global_variable::my_rank == 0 &&
        (two_stream_rt::rt_face_flux_ready() || mlt_dump_calls_ > 4)) {
      mlt_dumped_ = true;
      std::ofstream df(mlt_dump_);
      df.precision(6);
      df << std::scientific;
      if (mltmean) {
        auto hm = Kokkos::create_mirror_view(fmean_);
        Kokkos::deep_copy(hm, fmean_);
        df << "# red_giant MLT faces, SHELL MEAN (problem/mlt_mean), first call (cgs)\n"
           << "# i r T_f p_f grad grad_ad x F_mlt F_rad F_used F_req"
           << " F_raddiff w_blend F_rad_col0 F_conv_res grad_rad D F_2s\n";
        for (int i = is+1; i <= ie; ++i) {
          df << i;
          for (int q = 0; q < 17; ++q) df << " " << hm(i,q);
          df << "\n";
        }
      } else {
        auto hd = Kokkos::create_mirror_view(fdiag);
        Kokkos::deep_copy(hd, fdiag);
        df << "# red_giant MLT faces, block 0, first column, first source call (cgs)\n"
           << "# i grad grad_ad H_p c_p v_conv F_mlt F_cap F_used\n";
        for (int i = is; i <= ie+1; ++i) {
          df << i;
          for (int q = 0; q < 8; ++q) df << " " << hd(i,q);
          df << "\n";
        }
      }
    }
    mlt_fconv_ready_ = true;
    if (!(mlt_split_dep_ && rt_strang_))
    par_for("rg_mlt_div", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      Real div;
      if (curv) {
        div = (area1(m,k,j,i)*fconv(m,k,j,i) - area1(m,k,j,i+1)*fconv(m,k,j,i+1))
              /volume(m,k,j,i);
      } else {
        div = (fconv(m,k,j,i) - fconv(m,k,j,i+1))/((x1hi - x1lo)/indcs.nx1);
      }
      u0(m,IEN,k,j,i) += bdt*div;
    });
    RGNanScan(pm, "MLT_flux_divergence");
  }

  // --- OPEN INNER BOUNDARY: the CO5BOLD (Freytag+ 2012) relaxation of the lowest active
  // layer, the same treatment solar_convection.cpp uses, on a spherical shell.  Four
  // steps: relax UPFLOWING gas toward the deep adiabat (this is the energy input and it
  // sets the emergent luminosity), damp pressure toward the shell mean, restore the mean
  // density, and remove the net radial mass flux.  The means are global over the whole
  // shell at i = is, so they need a reduction across ranks.
  if (open_inner_) {
    auto &mb_bcs = pmbp->pmb->mb_bcs;
    const Real rgas = rgas_, igm1 = 1.0/gm1_, gamma = gm1_ + 1.0;
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real p_base = p_base_/punit_, t_base = t_base_;
    const Real ccs = cs_change_, ccp = cp_change_;
    auto &x1f_r = pmbp->pcoord->xx1f;
    // Pass 0: cancel the NET mass the inner face delivered this stage.  The fluxes the
    // RKUpdate just applied are still in uflx, so this is exact per stage rather than a
    // correction chasing a drift.
    if (open_nomassflux_) {
      auto &flx1 = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
      auto &area1_ = pmbp->pcoord->area.x1f;
      auto &vol0_ = pmbp->pcoord->volume;
      Real sM = 0.0, sV = 0.0, sH = 0.0, sD = 0.0, sPv = 0.0;
      Kokkos::parallel_reduce("rg_co50", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
      KOKKOS_LAMBDA(const int m, Real &aM, Real &aV, Real &aH, Real &aD, Real &aPv) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real vol = vol0_(m,k,j,is);
            const Real d = u0(m,IDN,k,j,is);
            const Real cc = cs_g ? ccell_g(m,k,j) : 0.0;
            Real ei = u0(m,IEN,k,j,is)
                      - CsKinetic(d, u0(m,IM1,k,j,is), u0(m,IM2,k,j,is),
                                  u0(m,IM3,k,j,is), cc);
            if (etotgrav) ei -= d*phicc(m,k,j,is);
            aM += area1_(m,k,j,is)*flx1(m,IDN,k,j,is);
            aV += vol;
            aH += (u0(m,IEN,k,j,is) + eos.Pressure(d, ei))*vol;   // total enthalpy
            aD += d*vol;
            aPv += u0(m,IM1,k,j,is)*vol;
          }
        }
      }, sM, sV, sH, sD, sPv);
#if MPI_PARALLEL_ENABLED
      {
        Real g[5] = {sM, sV, sH, sD, sPv};
        MPI_Allreduce(MPI_IN_PLACE, g, 5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
        sM = g[0]; sV = g[1]; sH = g[2]; sD = g[3]; sPv = g[4];
      }
#endif
      if (sV > 0.0 && sD > 0.0) {
        const Real dd = -bdt*sM/sV;      // code density, the same in every base cell
        const Real hbar = sH/sD;         // shell-mean specific total enthalpy
        const Real vbar = sPv/sD;        // shell-mean radial velocity
        par_for("rg_co50b", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
          if (u0(m,IDN,k,j,is) + dd <= 0.0) return;
          u0(m,IDN,k,j,is) += dd;
          u0(m,IM1,k,j,is) += dd*vbar;
          u0(m,IEN,k,j,is) += dd*hbar;
        });
      }
    }
    // Pass A: shell means of density and pressure
    Real srho0 = 0.0, sumP = 0.0, sN = 0.0;
    Kokkos::parallel_reduce("rg_co5A", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &a, Real &b, Real &c) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real r = w0(m,IDN,k,j,is);
            a += r;
            b += eos.Pressure(r, w0(m,IEN,k,j,is));
            c += 1.0;
          }
        }
      }
    }, srho0, sumP, sN);
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {srho0, sumP, sN};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      srho0 = g[0]; sumP = g[1]; sN = g[2];
    }
#endif
    const Real meanP = (sN > 0.0) ? sumP/sN : 0.0;
    // Pass B: the inflow entropy relaxation, then the pressure damping.  Both rewrite
    // u0(IEN); the budget (problem/open_budget) accumulates each one separately, which
    // is the only way to see which of them is the star's energy source and which is not.
    auto &vol_ = pmbp->pcoord->volume;
    const Real eunit = pmbp->punit->pressure_cgs()
                       *SQR(pmbp->punit->length_cgs())*pmbp->punit->length_cgs();
    Real dEent = 0.0, dEprs = 0.0;
    Kokkos::parallel_reduce("rg_co5B", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &sEent, Real &sEprs) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      const Real dr_ = x1f_r(m,is+1) - x1f_r(m,is);
      Real r = u0(m,IDN,k,j,is);
      const Real v1 = u0(m,IM1,k,j,is)/r, v2 = u0(m,IM2,k,j,is)/r;
      const Real v3 = u0(m,IM3,k,j,is)/r;
      const Real Ein = u0(m,IEN,k,j,is);
      Real ei = u0(m,IEN,k,j,is) - 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) ei -= r*phicc(m,k,j,is);
      Real es = ei/r;
      Real P = eos.Pressure(r, ei);
      Real g1 = eos.Gamma1(r, ei);
      const Real cs = sqrt(g1*P/r);
      if (v1 > 0.0) {
        // the point on the adiabat through the base state AT THIS CELL'S PRESSURE; a
        // fractional step toward it at constant pressure IS the entropy relaxation
        const Real gad = GradAd(eos, gamma, rgas, p_base*punit_, t_base);
        const Real T_ad = t_base*pow(P/p_base, gad);
        const Real r_ad = DensFromPT(eos, rgas, P, T_ad);
        if (r_ad > 0.0) {
          const Real es_ad = EintFromDensT(eos, rgas, igm1, r_ad, T_ad)/r_ad;
          const Real rlx = ccs*bdt*cs/dr_;
          r  += rlx*(r_ad - r);
          es += rlx*(es_ad - es);
        }
      }
      {
        Real Em = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
        if (etotgrav) Em += r*phicc(m,k,j,is);
        sEent += (Em - Ein)*vol_(m,k,j,is)*eunit;
        sEprs -= Em*vol_(m,k,j,is)*eunit;    // completed after the pressure damping
      }
      // damp pressure toward the shell mean, adiabatically
      Real P1 = eos.Pressure(r, r*es);
      Real g1p = eos.Gamma1(r, r*es);
      const Real cs2 = g1p*P1/r;
      const Real rlxp = ccp*bdt*sqrt(cs2)/dr_;
      r  += rlxp*(1.0/cs2)*(meanP - P1);
      es += rlxp*(1.0/(g1p*r))*(meanP - P1);
      u0(m,IDN,k,j,is) = r;
      u0(m,IM1,k,j,is) = r*v1;
      u0(m,IM2,k,j,is) = r*v2;
      u0(m,IM3,k,j,is) = r*v3;
      Real E = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) E += r*phicc(m,k,j,is);
      u0(m,IEN,k,j,is) = E;
      sEprs += E*vol_(m,k,j,is)*eunit;
      }}
    }, dEent, dEprs);
    // Pass C: the shell means again, after those two steps
    Real srho2 = 0.0, srho2v = 0.0, sv = 0.0;
    Kokkos::parallel_reduce("rg_co5C", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &a, Real &b, Real &c) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real r = u0(m,IDN,k,j,is);
            const Real v1 = u0(m,IM1,k,j,is)/r;
            a += r; b += r*v1; c += v1;
          }
        }
      }
    }, srho2, srho2v, sv);
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {srho2, srho2v, sv};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      srho2 = g[0]; srho2v = g[1]; sv = g[2];
    }
#endif
    const Real drho4 = (sN > 0.0) ? (srho0 - srho2)/sN : 0.0;
    const Real cvel = (srho0 > 0.0) ? (srho2v + drho4*sv)/srho0 : 0.0;
    // Pass D: restore the mean density, and drive the net radial mass flux to zero
    Real dEden = 0.0;
    Kokkos::parallel_reduce("rg_co5D", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &sEden) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      const Real Ein = u0(m,IEN,k,j,is);
      const Real r0 = u0(m,IDN,k,j,is);
      const Real v1 = u0(m,IM1,k,j,is)/r0 - cvel;
      const Real v2 = u0(m,IM2,k,j,is)/r0, v3 = u0(m,IM3,k,j,is)/r0;
      Real ei = u0(m,IEN,k,j,is) - 0.5*r0*(SQR(u0(m,IM1,k,j,is)/r0) + v2*v2 + v3*v3);
      if (etotgrav) ei -= r0*phicc(m,k,j,is);
      const Real es = ei/r0;
      const Real r = r0 + drho4;
      if (!(r > 0.0)) continue;
      u0(m,IDN,k,j,is) = r;
      u0(m,IM1,k,j,is) = r*v1;
      u0(m,IM2,k,j,is) = r*v2;
      u0(m,IM3,k,j,is) = r*v3;
      Real E = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) E += r*phicc(m,k,j,is);
      u0(m,IEN,k,j,is) = E;
      sEden += (E - Ein)*vol_(m,k,j,is)*eunit;
      }}
    }, dEden);
    // The report.  Rates over the interval since the last one, and the running totals,
    // both in units of L: the boundary is supposed to feed the star, so anything that
    // is not O(L) is a bug in one of the three passes above.
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {dEent, dEprs, dEden};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      dEent = g[0]; dEprs = g[1]; dEden = g[2];
    }
#endif
    // Pass E: give back what Pass D took.  Uniformly in volume, so the shell's own
    // structure is untouched and only the spurious global term goes away.
    if (open_conserve_ && dEden != 0.0) {
      Real svol = 0.0;
      Kokkos::parallel_reduce("rg_co5Evol", Kokkos::RangePolicy<>(DevExeSpace(),0,nmb1+1),
      KOKKOS_LAMBDA(const int m, Real &a) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) { a += vol_(m,k,j,is); }
          }
        }
      }, svol);
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &svol, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
      if (svol > 0.0) {
        const Real de = -dEden/(svol*eunit);      // code energy density
        par_for("rg_co5E", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
          u0(m,IEN,k,j,is) += de;
        });
        dEden = 0.0;
      }
    }
    if (open_budget_ > 0) {
      open_dE_ent_ += dEent;
      open_dE_prs_ += dEprs;
      open_dE_den_ += dEden;
      if (pm->ncycle % open_budget_ == 0 && global_variable::my_rank == 0) {
        const Real tnow = pm->time*pmbp->punit->time_cgs();
        const Real dtw = tnow - open_t_last_;
        const Real iL = 1.0/open_lstar_;
        if (dtw > 0.0) {
          std::cout << "open BC budget: rate/L  entropy="
                    << (open_dE_ent_ - open_dE_ent_l_)/dtw*iL << " pressure="
                    << (open_dE_prs_ - open_dE_prs_l_)/dtw*iL << " density="
                    << (open_dE_den_ - open_dE_den_l_)/dtw*iL
                    << " | total erg = " << open_dE_ent_ + open_dE_prs_ + open_dE_den_
                    << " (t = " << tnow << " s)" << std::endl;
        }
        open_t_last_ = tnow;
        open_dE_ent_l_ = open_dE_ent_;
        open_dE_prs_l_ = open_dE_prs_;
        open_dE_den_l_ = open_dE_den_;
      }
    }
    RGNanScan(pm, "open_inner_BC");
  }

  // --- THE WALL IS IMPERMEABLE.  See wall_noflux_ above: cancel the mass the inner face
  // advected this stage and the energy it carried, cell by cell, leaving the momentum
  // flux (its pressure term holds the star up) and the conduction flux (it carries the
  // imposed luminosity through the same face) alone.
  if (!open_inner_ && wall_noflux_) {
    auto &flx1w = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
    auto &area1w = pmbp->pcoord->area.x1f;
    auto &volw = pmbp->pcoord->volume;
    auto &mb_bcs = pmbp->pmb->mb_bcs;
    par_for("rg_wallflux", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      const Real dm = bdt*area1w(m,k,j,is)*flx1w(m,IDN,k,j,is)/volw(m,k,j,is);
      if (dm == 0.0) return;
      const Real d = u0(m,IDN,k,j,is);
      if (!(d - dm > 0.0)) return;
      // Cancel EXACTLY the energy the face flux applied this stage (the hydro flux
      // array holds the advected total-energy + pressure-work flux; the conduction
      // flux is a separate operator and is untouched). Rebuilding it as dm*h from the
      // cell state was not exact: with etotgrav the cell holds rho*Phi at its centre
      // while the face flux carries the reconstructed face energy, and the mismatch
      // (either sign, depending on the reconstruction) drifted tot-E by 1e-4 in 2e6 s
      // (RG_fofc_long3, +8.8e-5) or -7e-5 in 2e5 s (RG_v4 with the h = gas-only variant).
      const Real dE = bdt*area1w(m,k,j,is)*flx1w(m,IEN,k,j,is)/volw(m,k,j,is);
      u0(m,IDN,k,j,is) -= dm;
      u0(m,IEN,k,j,is) -= dE;
    });
    RGNanScan(pm, "wall_noflux");
  }

  // --- THE SPONGE under the outer wall.  See sponge_on_ above.
  if (sponge_on_) {
    const Real r0s = pm->mesh_size.x1min, r1s = pm->mesh_size.x1max;
    const Real zs = r0s + sponge_zbot_*(r1s - r0s);
    const Real izw = (r1s > zs) ? 1.0/(r1s - zs) : 0.0;
    const Real sc = sponge_c_;
    auto &x1vs = pmbp->pcoord->x1v;
    auto &dx1s = pmbp->pcoord->dx1;
    auto &ccell_s = pmbp->pcoord->cos_cell;    // see CsKinetic
    const bool cs_s = pm->use_cubed_sphere;
    const Real gm1l = gm1_, gamma_ = gm1_ + 1.0;   // locals: no host global in the lambda
    if (spgguard_cnt_ == nullptr) {
      spgguard_cnt_ = new DvceArray1D<int>("rg_spgguard_cnt", 1);
      spgguard_rec_ = new DvceArray1D<Real>("rg_spgguard_rec", 8);
      Kokkos::deep_copy(*spgguard_cnt_, 0);
      Kokkos::deep_copy(*spgguard_rec_, 0.0);
    }
    auto scnt = *spgguard_cnt_;
    auto srec = *spgguard_rec_;
    par_for("rg_sponge", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xc = curv ? x1vs(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
      if (!(xc > zs)) return;
      const Real ramp = SQR((xc - zs)*izw);
      const Real d = u0(m,IDN,k,j,i);
      const Real m1 = u0(m,IM1,k,j,i), m2 = u0(m,IM2,k,j,i), m3 = u0(m,IM3,k,j,i);
      // the kinetic energy on THIS grid's basis (CsKinetic), not the orthogonal sum
      const Real cc = cs_s ? ccell_s(m,k,j) : 0.0;
      const Real ke = CsKinetic(d, m1, m2, m3, cc);
      Real ei = u0(m,IEN,k,j,i) - ke;
      if (etotgrav) ei -= d*phicc(m,k,j,i);
      // THE SPONGE RUNS BEFORE ConsToPrim, so what it reads out of u0 has had no floor
      // applied to it.  The atmosphere it damps is hypersonic (Mach 470 measured at
      // i = 309 at the death of T7/T8/T9), so E - E_kin - d phi can land non-positive
      // after a single flux update -- which is ordinary, and is exactly what the
      // pressure floor a few tasks later exists to repair.  The tabulated EOS, however,
      // takes log10(e): handed a non-positive one it returns NaN, cs2 is NaN, and the
      // cell is NaN in one step with no precursor of any kind.  This is the ONLY place
      // in the run that evaluates the EOS on an unfloored energy.  A cell whose state
      // is not yet physical simply does not get damped this step.
      if (!(d > 0.0) || !(ei > 0.0)) {
        if (Kokkos::atomic_fetch_add(&scnt(0), 1) == 0) {
          srec(0) = static_cast<Real>(m);   srec(1) = static_cast<Real>(k);
          srec(2) = static_cast<Real>(j);   srec(3) = static_cast<Real>(i);
          srec(4) = d;                      srec(5) = ei;
          srec(6) = ke;
          srec(7) = etotgrav ? d*phicc(m,k,j,i) : 0.0;
        }
        return;
      }
      const Real cs2 = eos.IsGeneral() ? eos.Gamma1(d, ei)*eos.Pressure(d, ei)/d
                                       : gamma_*gm1l*ei/d;
      const Real dz = curv ? dx1s(m,k,j,i) : (x1hi - x1lo)/indcs.nx1;
      const Real fac = 1.0/(1.0 + sc*ramp*bdt*sqrt(cs2)/dz);
      // scaling the whole momentum vector by fac scales the kinetic energy by fac^2 on
      // any basis, so the internal energy is left exactly where it was
      u0(m,IM1,k,j,i) = fac*m1;
      u0(m,IM2,k,j,i) = fac*m2;
      u0(m,IM3,k,j,i) = fac*m3;
      Real E = ei + fac*fac*ke;
      if (etotgrav) E += d*phicc(m,k,j,i);
      u0(m,IEN,k,j,i) = E;
    });
    // Say so the first time the guard fires, with the state that tripped it, and stay
    // quiet afterwards.  One device sync on that first call only.
    if (!spgguard_reported_) {
      auto hs = Kokkos::create_mirror_view(scnt);
      Kokkos::deep_copy(hs, scnt);
      if (hs(0) > 0) {
        auto hr = Kokkos::create_mirror_view(srec);
        Kokkos::deep_copy(hr, srec);
        std::cout << "### red_giant sponge guard FIRED (rank " << global_variable::my_rank
                  << ", cycle " << pm->ncycle << ", t = " << pm->time << "): "
                  << hs(0) << " cell(s) not damped; first (m,k,j,i) = ("
                  << static_cast<int>(hr(0)) << "," << static_cast<int>(hr(1)) << ","
                  << static_cast<int>(hr(2)) << "," << static_cast<int>(hr(3))
                  << ") d = " << hr(4) << " ei = " << hr(5) << " ke = " << hr(6)
                  << " dphi = " << hr(7) << std::endl;
        spgguard_reported_ = true;
      }
    }
    RGNanScan(pm, "sponge");
  runaway_scan::Scan(pm, "sponge");
  }

  // --- THE TAU-BASED TOP SPONGE (problem/vdamp_top_tau; see the declaration of
  // vdamp_tau_), ported from box_convection.cpp:2703-2816.  It damps the RADIAL momentum
  // only and writes the kinetic-energy change back to IEN, so the internal energy is left
  // exactly where it was.  On the cubed sphere IM1 is orthogonal to the angular pair (the
  // cos_cell cross term couples IM2 and IM3 alone), so 0.5 (m1n^2 - m1o^2)/rho is the
  // exact kinetic-energy change here as well as on a Cartesian grid.
  {
    Conduction *pcsp = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
    // rad_w_built guards the one stage in which the column tau does not exist yet: an
    // all-zero tau would read as "thin everywhere" and damp the whole envelope
    const bool vdamp_on = (vdamp_tau_ > 0.0) && (pcsp != nullptr) && pcsp->rad_w_built;
    if (vdamp_on) {
      const Real vd_hi = vdamp_tau_, vd_lo = vdamp_tau_/3.0;
      const Real vd_rate = bdt/vdamp_time_;
      auto vtauf = pcsp->rad_tauf;
      par_for("rg_vdtop", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        // the cell's LOWER face carries the larger tau, the conservative (weaker) choice
        const Real f = 1.0 - RadBlendWeight(vtauf(m,k,j,i), vd_lo, vd_hi);
        if (f > 0.0) {
          const Real dc = u0(m,IDN,k,j,i);
          const Real m1o = u0(m,IM1,k,j,i);
          const Real m1n = m1o/(1.0 + f*vd_rate);
          u0(m,IM1,k,j,i) = m1n;
          u0(m,IEN,k,j,i) += 0.5*(SQR(m1n) - SQR(m1o))/dc;
        }
      });
      // the sponge's radial range, once, at the first stage it runs: the column tau only
      // exists after BuildRadWeights, so it cannot be known at start-up
      if (!vdamp_printed_) {
        vdamp_printed_ = true;
        Real rlow = std::numeric_limits<Real>::max();
        Kokkos::parallel_reduce("rg_vdtop_range",
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                               {nmb1+1,ke+1,je+1,ie+1}),
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmin) {
          if (1.0 - RadBlendWeight(vtauf(m,k,j,i), vd_lo, vd_hi) > 0.0) {
            const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
            const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
            lmin = (xc < lmin) ? xc : lmin;
          }
        }, Kokkos::Min<Real>(rlow));
#if MPI_PARALLEL_ENABLED
        Real rg_ = rlow;
        MPI_Allreduce(&rg_, &rlow, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
        if (global_variable::my_rank == 0) {
          if (rlow < pm->mesh_size.x1max) {
            std::printf("  red_giant top v1 sponge active: r >= %.5e, f = 1 at tau <= "
                        "%.4e, timescale %.4e s\n", rlow, vd_lo, vdamp_time_);
          } else {
            std::printf("  red_giant top v1 sponge: NO cell has tau < %.4e -- inert\n",
                        vd_hi);
          }
        }
      }
      RGNanScan(pm, "vdamp_top");
    }
  }

  // --- THE BOTTOM SPONGE (problem/vdamp_bot_cells; see the declaration of vdb_cells_),
  // ported from box_convection.cpp:2818-2879.  It runs on the u0 the sponges above left
  // and BEFORE the radiation source, so nothing injected this stage is cancelled by it.
  if (vdb_cells_ > 0) {
    const int nb = vdb_cells_;
    if (vdb_d_.extent_int(0) != nb) {
      Kokkos::realloc(vdb_d_, nb);
      Kokkos::realloc(vdb_h_, nb);
    }
    auto vdb = vdb_d_;
    if (vdb_mean_) {
      // one team per radial index of the layer reduces v1 over that SHELL's (m,k,j); the
      // sums are Allreduced (every rank applies the mean) and divided by the global count
      const int lnx2 = indcs.nx2, lnx3 = indcs.nx3;
      const int nkj = (nmb1+1)*lnx3*lnx2;
      const int gnx2 = pm->mesh_indcs.nx2, gnx3 = pm->mesh_indcs.nx3;
      const int npanel = pm->use_cubed_sphere ? 6 : 1;
      Kokkos::TeamPolicy<> vpol(DevExeSpace(), nb, Kokkos::AUTO);
      Kokkos::parallel_for("rg_vdbot_mean", vpol,
      KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tmember) {
        const int i = is + tmember.league_rank();
        Real vs = 0.0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkj),
        [&](const int idx, Real &ls) {
          const int m = idx/(lnx3*lnx2);
          const int kj = idx - m*(lnx3*lnx2);
          const int k = ks + kj/lnx2;
          const int j = js + (kj - (kj/lnx2)*lnx2);
          ls += u0(m,IM1,k,j,i)/u0(m,IDN,k,j,i);
        }, Kokkos::Sum<Real>(vs));
        Kokkos::single(Kokkos::PerTeam(tmember), [&]() { vdb(i-is) = vs; });
      });
      Kokkos::fence();
      Kokkos::deep_copy(vdb_h_, vdb_d_);
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, vdb_h_.data(), nb, MPI_ATHENA_REAL, MPI_SUM,
                    MPI_COMM_WORLD);
#endif
      const Real fpl = 1.0/static_cast<Real>(gnx2*gnx3*npanel);
      for (int q=0; q<nb; ++q) vdb_h_(q) *= fpl;
      Kokkos::deep_copy(vdb_d_, vdb_h_);
    }
    const Real vb_rate = bdt/vdb_time_;
    const bool vb_mean = vdb_mean_;
    par_for("rg_vdbot", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, is+nb-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real f = 0.5*(1.0 + cos(M_PI*static_cast<Real>(i-is)/nb));
      const Real gg = 1.0 - exp(-f*vb_rate);
      const Real dc = u0(m,IDN,k,j,i);
      const Real m1o = u0(m,IM1,k,j,i);
      const Real m1n = m1o - gg*(vb_mean ? dc*vdb(i-is) : m1o);
      u0(m,IM1,k,j,i) = m1n;
      u0(m,IEN,k,j,i) += 0.5*(SQR(m1n) - SQR(m1o))/dc;
    });
    RGNanScan(pm, "vdamp_bot");
  }
  if (vda_until_ > 0.0 && pm->time < vda_until_) {
    const Real ga = 1.0 - exp(-bdt/vda_time_);
    par_for("rg_vdall", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dc = u0(m,IDN,k,j,i);
      const Real m1o = u0(m,IM1,k,j,i);
      const Real m1n = m1o*(1.0 - ga);
      u0(m,IM1,k,j,i) = m1n;
      u0(m,IEN,k,j,i) += 0.5*(SQR(m1n) - SQR(m1o))/dc;
    });
  }

  // --- the optically thin layers: the correlated-k or the grey two-stream if one of
  // them is on, else the old grey Eddington relaxation.  problem/rt_strang and
  // problem/rt_once_per_cycle take this source OUT of the stage: RedGiantRTSplit then
  // makes exactly the same call, with a different dt and at a different point in the
  // cycle, and nothing else moves.
  if ((rt_ck_ || rt_grey_) && !rt_strang_ && !rt_once_) {
    RedGiantRTSweep(pm, bdt);
    // problem/rt_surface_dt: the emergent top-face flux the sweep has just produced
    RedGiantSurfaceTick(pm);
  }
  // --- the grey relaxation (see UserProblem): with weight 1 - w each cell
  // decays toward the Eddington temperature of its optical depth on its radiative time
  if (relax_ && !rt_ck_ && !rt_grey_) {
    Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
    auto &wf = pc->rad_w;
    auto &tf = pc->rad_tauf;
    auto ktab = ktab_;
    auto klT = klT_;
    auto klD = klD_;
    const int knT = knT_, knD = knD_;
    const Real teff4 = SQR(SQR(teff_));
    const Real kfac = kfac_;
    const Real rgas = rgas_, igm1 = 1.0/gm1_, gm1 = gm1_;
    const Real dunit = pmbp->punit->density_cgs();
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real dt_cgs = bdt*pmbp->punit->time_cgs();
    par_for("rg_relax", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real wc = 0.5*(wf(m,k,j,i) + wf(m,k,j,i+1));
      if (wc >= 1.0) return;
      const Real tauc = 0.5*(tf(m,k,j,i) + tf(m,k,j,i+1));
      const Real d = w0(m,IDN,k,j,i);
      const Real e = w0(m,IEN,k,j,i);
      const Real p = eos.Pressure(d, e);
      const Real t = TempKelvin(eos, rgas, d, e, p);
      const Real teq = pow(0.75*teff4*(tauc + 2.0/3.0), 0.25);
      const Real eeq = EintFromDensT(eos, rgas, igm1, d, teq);
      // rho c_v [erg/cm^3/K] by a finite difference of the EOS's own e(rho, T)
      const Real e1 = EintFromDensT(eos, rgas, igm1, d, 1.01*t);
      const Real e0 = EintFromDensT(eos, rgas, igm1, d, t);
      // fall back to the ideal-gas value if the table's finite difference misbehaves
      Real rcv = (e1 - e0)/(0.01*t)*punit_;
      if (!(rcv > 0.0)) rcv = 1.5*p*punit_/t;
      const Real kp = kfac*KappaTab(ktab, klT, klD, knT, knD, t, d*dunit);
      const Real trad = rcv/(4.0*kp*(d*dunit)*kSigmaSB*t*t*t);
      if (!(trad > 0.0)) return;
      const Real de = -(1.0 - wc)*(e - eeq)*(1.0 - exp(-dt_cgs/trad));
      u0(m,IEN,k,j,i) += de;
      (void) gm1;
    });
    RGNanScan(pm, "grey_relax");
  runaway_scan::Scan(pm, "grey_relax");
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantRTSplit
//! \brief problem/rt_strang (and problem/rt_once_per_cycle): the two-stream as an
//! OPERATOR-SPLIT source around the time integrator instead of an in-stage one.
//! Hydro::RTStrangSplit calls this twice a cycle with bdt = dt/2 (or once with the full
//! dt under rt_once_per_cycle), on the state at the start and at the end of the cycle,
//! and runs ConToPrim after each.  It is exactly the call RedGiantGravity makes when the
//! split is off, with a different dt and at a different point in the cycle.
//!
//! With problem/rt_split_transverse the implicit transverse (ADI) operator travels with
//! it, over the same bdt and on the state the column solve has just relaxed: its x2/x3
//! ghosts are the last exchange's, exactly as they are for the in-stage task it replaces.

// problem/e_ledger (cycles, 0 = off; single rank only): which OPERATOR changes the
// domain's sum of u0(IEN) dV.  el_(0) = the two-stream sweep (with its radiative force),
// el_(1) = the split MLT deposit, el_(2) = everything in RedGiantGravity; the rest of
// the change since the first call is the hydro update, i.e. the two radial faces.
Real el_[3] = {0.0, 0.0, 0.0}, el_e0_ = 0.0, el_t0_ = -1.0;

Real RgTotE(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nji = nx2*nx1, nkji = nx3*nx2*nx1;
  auto &volume = pmbp->pcoord->volume;
  Real sum = 0.0;
  Kokkos::parallel_reduce("rg_tote", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &acc) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji;
    const int j = (idx - m*nkji - k*nji)/nx1;
    const int i = (idx - m*nkji - k*nji - j*nx1) + is;
    acc += u0(m,IEN,k+ks,j+js,i)*volume(m,k+ks,j+js,i);
  }, Kokkos::Sum<Real>(sum));
  return sum;
}

// problem/cs_max [cm/s], 0 = off: a SOUND-SPEED CEILING applied as a density floor,
//     rho >= (10/9) e_int / cs_max^2
// (c_s^2 = 4 e/(9 rho) radiation dominated, 10 e/(9 rho) for a gamma = 5/3 gas).  The
// 1-D He4 column evacuates a cavity under the Fe bump (0.85-0.87 R, rho 4e-9 -> 1.6e-13
// at T = 1.5e5 K); the EOS keeps LTE radiation pressure in a hot cell whatever its
// density, so c_s = 4e9 cm/s there and dt falls from 17 s to 0.3 s (tests_r12/sd3).
// The added mass is at rest and carries no energy (see the kernel).

void RgCsCeiling(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nji = nx2*nx1, nkji = nx3*nx2*nx1;
  auto &volume = pmbp->pcoord->volume;
  auto &size = pmbp->pmb->mb_size;
  auto &x1v_ = pmbp->pcoord->x1v;
  const bool curv = curv_, etotgrav = etotgrav_;
  const Real gm = gm_, rin = rin_, x1min = x1min_;
  const Real fac = (10.0/9.0)/SQR(cs_max_/pmbp->punit->velocity_cgs());
  Real dm = 0.0;
  Kokkos::parallel_reduce("rg_csmax", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmbp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &acc) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k-ks)*nji)/nx1 + js;
    const int i = (idx - m*nkji - (k-ks)*nji - (j-js)*nx1) + is;
    const Real d = u0(m,IDN,k,j,i);
    const Real dmin = fac*w0(m,IEN,k,j,i);
    if (d > 0.0 && d < dmin) {
      // the mass is added AT REST: momentum and total energy are kept, so the cell's
      // kinetic energy falls by d/dmin and the difference becomes heat.  Keeping the
      // VELOCITY instead injected kinetic energy (tests_r12/cm_3.0e8: tot-E +40 %).
      const Real dd = dmin - d;
      u0(m,IDN,k,j,i) = dmin;
      if (etotgrav) {
        const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
        const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, nx1, x1lo, x1hi);
        u0(m,IEN,k,j,i) += dd*PotAt(gm, rin, RadiusOf(curv, xc, rin, x1min));
      }
      acc += dd*volume(m,k,j,i);
    }
  }, Kokkos::Sum<Real>(dm));
  cs_mass_added_ += dm;
}

void RedGiantGravityLedger(Mesh *pm, Real bdt) {
  const Real e0 = (eledger_ > 0) ? RgTotE(pm) : 0.0;
  if (eledger_ > 0 && el_t0_ < 0.0) { el_t0_ = pm->time; el_e0_ = e0; }
  RedGiantGravity(pm, bdt);
  if (eledger_ > 0) el_[2] += RgTotE(pm) - e0;
  if (cs_max_ > 0.0) RgCsCeiling(pm);
}

// the same ledger PER RADIAL SHELL, for the sweep alone: sign*sum(u0(IEN) dV) into
// slot (0,i), and with sign = 0 the telescoped face fluxes the sweep reports,
// bdt*(A_i F2s_i - A_i+1 F2s_i+1), into slot (1,i).  Written to e_ledger_shell.txt.
DvceArray2D<Real> el_sh_;
Real el_tsum_ = 0.0;      // sum of the sweep dt the shell ledger has integrated
void RgShellAcc(Mesh *pm, const Real sign, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  if (el_sh_.extent(0) == 0) Kokkos::realloc(el_sh_, 2, indcs.nx1+1);
  auto sh = el_sh_;
  auto &volume = pmbp->pcoord->volume;
  auto &area1 = pmbp->pcoord->area.x1f;
  if (sign != 0.0) {
    par_for("rg_elsh", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Kokkos::atomic_add(&sh(0,i-is), sign*u0(m,IEN,k,j,i)*volume(m,k,j,i));
    });
  } else if (two_stream_rt::rt_face_flux_ready()) {
    // slot (1,i): the time integral of the FACE luminosity, sum_kj A_i F_i bdt
    auto fb = two_stream_rt::rt_face_flux();
    const int nblk = two_stream_rt::rt_face_nblk();
    par_for("rg_elshf", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real ft = 0.0;
      for (int b=0; b<nblk; ++b) ft += fb(m,b,i,k,j);
      Kokkos::atomic_add(&sh(1,i-is), bdt*area1(m,k,j,i)*ft);
      // slot (0,nx1) is otherwise unused: the summed face AREA of the bottom face, and
      // of the top face times 1e-30 packed apart is not worth it -- bottom only
      if (i == is) Kokkos::atomic_add(&sh(0,ie+1-is), bdt*area1(m,k,j,i));
    });
  }
}

// problem/rt_no_heat (diagnostic): the split operator keeps its MOMENTUM change (the
// radiative force) but deposits no heat: u0(IEN) is put back to its entry value plus the
// kinetic-energy change.  An adiabatic star with the force law intact.
DvceArray5D<Real> rt_nh_save_;
void RedGiantRTSplitBody(Mesh *pm, Real bdt);
void RedGiantRTSplitLedger(Mesh *pm, Real bdt);
void RedGiantRTSplit(Mesh *pm, Real bdt) {
  if (!rt_no_heat_) { RedGiantRTSplitLedger(pm, bdt); return; }
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &u0 = (pmbp->pmhd != nullptr) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  if (rt_nh_save_.extent(0) == 0) {
    Kokkos::realloc(rt_nh_save_, u0.extent(0), 2, u0.extent(2), u0.extent(3),
                    u0.extent(4));
  }
  auto sv = rt_nh_save_;
  const int n1 = u0.extent(4) - 1, n2 = u0.extent(3) - 1, n3 = u0.extent(2) - 1;
  const int nm = u0.extent(0) - 1;
  par_for("rg_nh0", DevExeSpace(), 0, nm, 0, n3, 0, n2, 0, n1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    sv(m,0,k,j,i) = u0(m,IEN,k,j,i);
    sv(m,1,k,j,i) = u0(m,IM1,k,j,i);
  });
  RedGiantRTSplitLedger(pm, bdt);
  par_for("rg_nh1", DevExeSpace(), 0, nm, 0, n3, 0, n2, 0, n1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = u0(m,IDN,k,j,i);
    const Real dke = (d > 0.0)
                   ? 0.5*(SQR(u0(m,IM1,k,j,i)) - SQR(sv(m,1,k,j,i)))/d : 0.0;
    u0(m,IEN,k,j,i) = sv(m,0,k,j,i) + dke;
  });
}

void RedGiantRTSplitLedger(Mesh *pm, Real bdt) {
  if (eledger_ <= 0) { RedGiantRTSplitBody(pm, bdt); return; }
  static int ncall = 0;
  const Real e0 = RgTotE(pm);
  if (el_t0_ < 0.0) { el_t0_ = pm->time; el_e0_ = e0; }
  const Real mlt0 = el_[1];
  RgShellAcc(pm, -1.0, bdt);
  RedGiantRTSplitBody(pm, bdt);
  RgShellAcc(pm, 1.0, bdt);
  RgShellAcc(pm, 0.0, bdt);
  el_tsum_ += bdt;
  el_[0] += RgTotE(pm) - e0 - (el_[1] - mlt0);
  if ((++ncall % 2 == 0) && (pm->ncycle % eledger_ == 0)) {
    const Real de = RgTotE(pm) - el_e0_;
    std::cout << std::scientific << std::setprecision(6) << "### e_ledger t=" << pm->time
              << " dE=" << de << " rt=" << el_[0] << " mlt=" << el_[1] << " srcs="
              << el_[2] << " hydro=" << de - el_[0] - el_[1] - el_[2] << std::endl;
    auto hs = Kokkos::create_mirror_view(el_sh_);
    Kokkos::deep_copy(hs, el_sh_);
    std::ofstream sf("e_ledger_shell.txt");
    sf << "# t=" << el_tsum_ << " : i  dE_sweep(+mlt)  bdt*(A F2s)_i-(A F2s)_i+1 (code)\n"
       << std::scientific << std::setprecision(8);
    for (int i = 0; i <= pm->mb_indcs.nx1; ++i) {
      sf << i << " " << hs(0,i) << " " << hs(1,i) << "\n";
    }
  }
}

void RedGiantRTSplitBody(Mesh *pm, Real bdt) {
  if (!(rt_ck_ || rt_grey_)) return;
  // problem/mlt_split_deposit: div F_conv with the sweep's own dt (see the declaration),
  // BEFORE the sweep: the implicit column solve then starts from e + dt S, i.e. it is
  // backward Euler WITH the source, and a state with div F_2s + div F_conv = 0 is its
  // fixed point.  Deposited after the sweep it is not (tests_r12/fs0: dipole 1.25x).
  if (mlt_split_dep_ && rt_strang_ && mlt_alpha_ > 0.0) {
    if (!mlt_fconv_ready_) {
      mlt_owed_ += bdt;
    } else {
      MeshBlockPack *pmbp = pm->pmb_pack;
      const bool is_mhd = (pmbp->pmhd != nullptr);
      auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
      auto &indcs = pm->mb_indcs;
      const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
      const int ks = indcs.ks, ke = indcs.ke;
      const int nmb1 = pmbp->nmb_thispack - 1;
      auto &size = pmbp->pmb->mb_size;
      auto &area1 = pmbp->pcoord->area.x1f;
      auto &volume = pmbp->pcoord->volume;
      const bool curv = curv_;
      auto fconv = fconv_;
      const Real ddt = bdt + mlt_owed_;
      mlt_owed_ = 0.0;
      const Real em0 = (eledger_ > 0) ? RgTotE(pm) : 0.0;
      par_for("rg_mlt_div_split", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
        Real div;
        if (curv) {
          div = (area1(m,k,j,i)*fconv(m,k,j,i) - area1(m,k,j,i+1)*fconv(m,k,j,i+1))
                /volume(m,k,j,i);
        } else {
          div = (fconv(m,k,j,i) - fconv(m,k,j,i+1))/((x1hi - x1lo)/indcs.nx1);
        }
        u0(m,IEN,k,j,i) += ddt*div;
      });
      if (eledger_ > 0) el_[1] += RgTotE(pm) - em0;
      RGNanScan(pm, "MLT_flux_divergence_split");
    }
  }
  RedGiantRTSweep(pm, bdt);
  if (rt_split_tr_) {
    MeshBlockPack *pmbp = pm->pmb_pack;
    const bool is_mhd = (pmbp->pmhd != nullptr);
    Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
    auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
    const EOS_Data &eos = is_mhd ? pmbp->pmhd->peos->eos_data
                                 : pmbp->phydro->peos->eos_data;
    if (pc != nullptr && pc->rad_implicit_ang) {
      if (pc->rad_sts_all) {
        pc->StsConductionUpdate(u0, eos, bdt);
      } else {
        pc->ImplicitTransverseUpdate(u0, eos, bdt);
      }
      RGNanScan(pm, "RT_transverse_split");
    }
  }
  // the per-angular-cell surface dump, once per cycle (pm->time does not move between
  // the two Strang half steps, so the arming below fires on the first of them only)
  RedGiantSurfaceTick(pm);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantBC
//! \brief both radial walls: the ghost cells hold the initial column at their own
//! radius, with the radial velocity mirrored (a reflecting wall whose ghost state is in
//! hydrostatic balance with the interior, so the well-balanced scheme sees no jump)

void RedGiantBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  auto &x1v_ = pmbp->pcoord->x1v;
  // see the note at rg_ic: u0(IM2,IM3) are COVARIANT on the cubed sphere.  cos_cell is
  // allocated and filled over the FULL padded angular range, ghosts included.
  auto &ccell_b = pmbp->pcoord->cos_cell;
  const bool cs_b = pm->use_cubed_sphere;
  const bool curv = curv_, etotgrav = etotgrav_;
  const Real gm = gm_, rin = rin_, x1min = x1min_, rgas = rgas_, igm1 = 1.0/gm1_;
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  const Real rlo = rlo_, drf = drf_;
  const int nfine = nfine_;
  auto lnp = lnp_d_;
  auto tk = tk_d_;
  // OPEN inner boundary: continue the interior hydrostatically at ITS OWN temperature and
  // copy the velocity, so a plume crosses the boundary instead of bouncing off it.  The
  // pressure follows d ln p = -(rho/p) dPhi from the lowest active cell, which keeps the
  // ghost in balance with the interior rather than with the initial column -- tying it to
  // the column would fight the relaxation and pressurise the envelope.
  const bool open_in = open_inner_;
  // OPEN outer boundary: the same continuation at the other end.  See open_outer_.
  const bool open_out = open_outer_;
  const bool open_out_noin = open_outer_noinflow_;
  const bool open_out_lid = open_outer_lid_;
  // The AMBIENT MEDIUM is not hydrostatic -- a 400 K gas has a scale height of ~1e-3 R
  // at the outer wall -- so it falls in at the free-fall speed and the reflecting ghost
  // below (v1 mirrored) lets the top cells evacuate with nothing coming back.  Measured:
  // the outermost 18 layers went to vacuum in 5e3 s and the hole ate inward ~1 cell per
  // 1e3 s until the run died.  Where the ghost sits IN the background, copy the interior
  // velocity instead of mirroring it, so the boundary passes the flow and the fixed
  // (bg_rho, bg_temp) ghost REPLENISHES what falls in.  That makes the outer edge a
  // steady accretion boundary at the free-fall rate -- 4e-11 Msun/yr for the default
  // background, which is the rate the medium was chosen to deliver in the first place.
  const Real bgrj = (bg_rho_ > 0.0) ? bg_rjoin_ : -1.0;
  // problem/open_debug: print the open ghost's inputs and outputs for one column on the
  // first few boundary calls.  A ghost that comes out non-finite says nothing about
  // which of the four EOS calls did it.
  const bool open_dbg = (open_debug_ > 0) && (open_dbg_calls_++ < open_debug_);
  // the guard's counters, allocated once and never zeroed: the total is cumulative and is
  // read back only every 1000 cycles, which keeps a device sync out of every BC call
  if (bcguard_cnt_ == nullptr) {
    bcguard_cnt_ = new DvceArray1D<int>("rg_bcguard_cnt", 1);
    bcguard_rec_ = new DvceArray1D<Real>("rg_bcguard_rec", 8);
    Kokkos::deep_copy(*bcguard_cnt_, 0);
    Kokkos::deep_copy(*bcguard_rec_, 0.0);
  }
  auto gcnt = *bcguard_cnt_;
  auto grec = *bcguard_rec_;
  // The temperature window the tabulated EOS can actually invert in.  A continuation that
  // hands DensFromPT a T outside it is extrapolated, and an extrapolated inversion is
  // exactly what produces the sentinel; reject it before it reaches the ghost.
  const bool eostab = eos.tbl.active;
  const Real t_lo = eostab ? pow(10.0, eos.tbl.ymin) : 0.0;
  const Real t_hi = eostab ? pow(10.0, eos.tbl.ymax) : 1.0e30;
  // The INITIAL COLUMN at a ghost's own radius -- what `fill` uses, factored out so the
  // open ghosts can fall back to it cell by cell when their continuation is rejected.
  auto column_state = KOKKOS_LAMBDA(const int m, const int i, Real &d, Real &e) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    Real lp, t;
    ColumnAt(lnp, tk, nfine, rlo, drf, r, lp, t);
    d = DensFromPT(eos, rgas, exp(lp), t);
    e = EintFromDensT(eos, rgas, igm1, d, t);
  };
  // Is a hydrostatic continuation usable?  Everything finite and positive, the density
  // within three decades of the cell it continues (the sentinel is 14 decades away, and a
  // real ghost never moves by more than the ~0.1 % that one cell of scale height buys),
  // and the temperature inside the table's own inversion window.
  auto ghost_ok = KOKKOS_LAMBDA(const Real d_g, const Real e_g, const Real d_i,
                                const Real t_i, const Real p_g) {
    return isfinite(d_g) && isfinite(e_g) && isfinite(p_g) && isfinite(t_i)
           && (d_g > 0.0) && (e_g > 0.0) && (p_g > 0.0)
           && (d_g < 1.0e3*d_i) && (d_g > 1.0e-3*d_i)
           && (t_i > t_lo) && (t_i < t_hi);
  };
  // record the first rejection so the offending state is not lost, and count the rest
  auto guard_hit = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                                 const Real d_i, const Real e_i, const Real t_i,
                                 const Real p_g) {
    if (Kokkos::atomic_fetch_add(&gcnt(0), 1) == 0) {
      grec(0) = static_cast<Real>(m); grec(1) = static_cast<Real>(k);
      grec(2) = static_cast<Real>(j); grec(3) = static_cast<Real>(i);
      grec(4) = d_i; grec(5) = e_i; grec(6) = t_i; grec(7) = p_g;
    }
  };
  // WHERE THE GHOST FILLS READ THE INTERIOR CELL THEY CONTINUE: always from u0, never
  // from w0.  Two reasons, and the second is fatal.
  // (1) ApplyPhysicalBCs runs AFTER RKUpdate, the user source terms (gravity, WB, the
  // two-stream) and the implicit radial conduction, and BEFORE ConToPrim -- so w0 here
  // is the previous stage's inversion and predates every one of those edits to u0.  In a
  // smooth cell the difference is O(dt); at the photosphere, which is exactly where these
  // runs die, it is a systematic jump at the boundary face.
  // (2) w0 IS NOT RESTART STATE.  restart.cpp writes u0 only (ghosts included), and
  // Driver::InitBoundaryValuesAndPrimitives (driver.cpp) calls ApplyPhysicalBCs BEFORE
  // ConToPrim, so at the first boundary call of a restarted run w0 is still zero: the
  // `w0(...,is/ie) > 0` gates below then took the initial-column branch and overwrote the
  // ghosts the file had restored with a state the running boundary never produces.  That
  // made every radial ghost differ on a restart and ~17k active cells differ one cycle
  // later.  Reading u0 -- which IS restored bitwise, and which the running boundary reads
  // at exactly the same point of the update -- makes the fill idempotent, so a restart is
  // a bitwise continuation, and it is the physically correct state as well.
  // The velocity is RAISED with the gnomonic metric (the stored momentum is covariant,
  // see GnomonicEquiangleRaiseVel), which is the inverse of the lowering the ghost write
  // below performs, so the round trip is exact.
  auto state_i = [=] (const int m, const int k, const int j, const int im,
                      Real &d_i, Real &e_i, Real &v1, Real &v2, Real &v3) {
    const Real cc = cs_b ? ccell_b(m,k,j) : 0.0;
    d_i = u0(m,IDN,k,j,im);
    e_i = EintFromCons(u0, m, k, j, im, cc, cs_b, etotgrav,
                       etotgrav ? phicc(m,k,j,im) : 0.0);
    const Real q1 = u0(m,IM1,k,j,im), q2 = u0(m,IM2,k,j,im), q3 = u0(m,IM3,k,j,im);
    const Real det = 1.0 - cc*cc;
    v1 = q1/d_i;
    v2 = cs_b ? (q2 - cc*q3)/(d_i*det) : q2/d_i;
    v3 = cs_b ? (q3 - cc*q2)/(d_i*det) : q3/d_i;
  };
  auto fill_open = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                                 const int im) {
    Real d_i, e_i, v1, v2, v3;
    state_i(m, k, j, im, d_i, e_i, v1, v2, v3);
    const Real p_i = eos.Pressure(d_i, e_i);
    const Real t_i = TempKelvin(eos, rgas, d_i, e_i, p_i);
    const Real dphi = phicc(m,k,j,i) - phicc(m,k,j,im);
    const Real p_g = p_i*exp(-(d_i/p_i)*dphi);
    Real d_g = DensFromPT(eos, rgas, p_g, t_i);
    Real e_g = EintFromDensT(eos, rgas, igm1, d_g, t_i);
    if (!ghost_ok(d_g, e_g, d_i, t_i, p_g)) {
      guard_hit(m, k, j, i, d_i, e_i, t_i, p_g);
      column_state(m, i, d_g, e_g);
    }
    if (open_dbg && m == 0 && k == 2 && j == 2) {
      Kokkos::printf("rg_open_bc i=%d im=%d: d_i=%.6e e_i=%.6e p_i=%.6e t_i=%.6e "
                     "dphi=%.6e p_g=%.6e d_g=%.6e e_g=%.6e v1=%.6e\n",
                     i, im, d_i, e_i, p_i, t_i, dphi, p_g, d_g, e_g, v1);
    }
    w0(m,IDN,k,j,i) = d_g;
    w0(m,IEN,k,j,i) = e_g;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    const Real cc = cs_b ? ccell_b(m,k,j) : 0.0;
    u0(m,IDN,k,j,i) = d_g;
    u0(m,IM1,k,j,i) = d_g*v1;
    u0(m,IM2,k,j,i) = d_g*(v2 + cc*v3);
    u0(m,IM3,k,j,i) = d_g*(v3 + cc*v2);
    Real et = e_g + 0.5*d_g*(v1*v1 + v2*v2 + v3*v3 + 2.0*cc*v2*v3);
    if (etotgrav) et += d_g*phicc(m,k,j,i);
    u0(m,IEN,k,j,i) = et;
  };
  // The OUTER open ghost.  The same hydrostatic continuation as fill_open, kept separate
  // because the two ends want different things of the velocity: the inner boundary copies
  // it so plumes cross, this one clamps it to outflow so an extrapolated ghost cannot
  // push mass back in.  dphi is positive going up, so the exponential thins the ghost.
  auto fill_open_out = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                                     const int im) {
    Real d_i, e_i, v1, v2, v3;
    state_i(m, k, j, im, d_i, e_i, v1, v2, v3);
    const Real p_i = eos.Pressure(d_i, e_i);
    const Real t_i = TempKelvin(eos, rgas, d_i, e_i, p_i);
    const Real dphi = phicc(m,k,j,i) - phicc(m,k,j,im);
    const Real p_g = p_i*exp(-(d_i/p_i)*dphi);
    Real d_g = DensFromPT(eos, rgas, p_g, t_i);
    Real e_g = EintFromDensT(eos, rgas, igm1, d_g, t_i);
    if (!ghost_ok(d_g, e_g, d_i, t_i, p_g)) {
      guard_hit(m, k, j, i, d_i, e_i, t_i, p_g);
      column_state(m, i, d_g, e_g);
    }
    if (open_out_lid) {
      Real d_m, e_m, v2m, v3m;       // the mirror cell about the top face
      state_i(m, k, j, 2*im + 1 - i, d_m, e_m, v1, v2m, v3m);
      v1 = -v1;
    } else if (open_out_noin && v1 < 0.0) {
      v1 = 0.0;
    }
    if (open_dbg && m == 0 && k == 2 && j == 2) {
      Kokkos::printf("rg_open_bc_out i=%d im=%d: d_i=%.6e e_i=%.6e p_i=%.6e t_i=%.6e "
                     "dphi=%.6e p_g=%.6e d_g=%.6e e_g=%.6e v1=%.6e\n",
                     i, im, d_i, e_i, p_i, t_i, dphi, p_g, d_g, e_g, v1);
    }
    w0(m,IDN,k,j,i) = d_g;
    w0(m,IEN,k,j,i) = e_g;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    const Real cc = cs_b ? ccell_b(m,k,j) : 0.0;
    u0(m,IDN,k,j,i) = d_g;
    u0(m,IM1,k,j,i) = d_g*v1;
    u0(m,IM2,k,j,i) = d_g*(v2 + cc*v3);
    u0(m,IM3,k,j,i) = d_g*(v3 + cc*v2);
    Real et = e_g + 0.5*d_g*(v1*v1 + v2*v2 + v3*v3 + 2.0*cc*v2*v3);
    if (etotgrav) et += d_g*phicc(m,k,j,i);
    u0(m,IEN,k,j,i) = et;
  };
  auto fill = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            const int im) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    Real lp, t;
    ColumnAt(lnp, tk, nfine, rlo, drf, r, lp, t);
    const Real p = exp(lp);
    const Real d = DensFromPT(eos, rgas, p, t);
    const Real e = EintFromDensT(eos, rgas, igm1, d, t);
    // in the background: pass the flow (accretion boundary); elsewhere: reflecting wall
    Real d_wi, e_wi, vi1, v2, v3;
    state_i(m, k, j, im, d_wi, e_wi, vi1, v2, v3);
    const Real v1 = ((bgrj > 0.0) && (r >= bgrj)) ? vi1 : -vi1;
    w0(m,IDN,k,j,i) = d;
    w0(m,IEN,k,j,i) = e;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    const Real cc = cs_b ? ccell_b(m,k,j) : 0.0;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*(v2 + cc*v3);
    u0(m,IM3,k,j,i) = d*(v3 + cc*v2);
    Real et = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3 + 2.0*cc*v2*v3);
    if (etotgrav) et += d*PotAt(gm, rin, r);
    u0(m,IEN,k,j,i) = et;
  };
  par_for("rg_bc_x1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int n) {
    // AN OPEN GHOST NEEDS AN INTERIOR TO CONTINUE.  The boundary function runs before the
    // problem generator has filled the state -- twice, on the calls that set up the mesh
    // -- and the conserved variables are then exactly zero.  Continuing THAT
    // hydrostatically gives a pressure of NaN, and SolveDensity answers a NaN pressure
    // with its 1e6 g/cm^3
    // sentinel: a ghost 17 orders of magnitude denser than the cell it sits against,
    // which destroys the outermost active cell on the first step and collapses the
    // timestep by 2e4 before the run has produced a single output.  Fall back to the
    // initial column whenever the neighbouring active cell has no state yet; once it
    // has, the continuation is used and at t = 0 it reproduces the column anyway.
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      if (open_in && u0(m,IDN,k,j,is) > 0.0) {
        fill_open(m, k, j, is-1-n, is);      // always from the lowest ACTIVE cell
      } else {
        fill(m, k, j, is-1-n, is+n);
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      if (open_out && u0(m,IDN,k,j,ie) > 0.0) {
        fill_open_out(m, k, j, ie+1+n, ie);   // always from the highest ACTIVE cell
      } else {
        fill(m, k, j, ie+1+n, ie-n);
      }
    }
  });
  // Report the guard.  The counter is read back only on the first call of every 1000th
  // cycle, so a healthy run pays one device sync per 1000 cycles and prints nothing at
  // all; a run that starts rejecting ghosts says so, with the state that did it.
  if (pm->ncycle % 1000 == 0) {
    auto hc = Kokkos::create_mirror_view(gcnt);
    Kokkos::deep_copy(hc, gcnt);
    if (hc(0) > bcguard_seen_) {
      if (!bcguard_reported_) {
        auto hr = Kokkos::create_mirror_view(grec);
        Kokkos::deep_copy(hr, grec);
        std::cout << "### red_giant open-ghost guard FIRED (rank "
                  << global_variable::my_rank << "): first at (m,k,j,i) = ("
                  << static_cast<int>(hr(0)) << "," << static_cast<int>(hr(1)) << ","
                  << static_cast<int>(hr(2)) << "," << static_cast<int>(hr(3))
                  << ") d_i=" << hr(4) << " e_i=" << hr(5) << " t_i=" << hr(6)
                  << " p_g=" << hr(7) << "; the initial column is used there instead"
                  << std::endl;
        bcguard_reported_ = true;
      } else {
        std::cout << "red_giant open-ghost guard: " << (hc(0) - bcguard_seen_)
                  << " fallbacks since the last report, " << hc(0) << " total (rank "
                  << global_variable::my_rank << ", cycle " << pm->ncycle << ")"
                  << std::endl;
      }
      bcguard_seen_ = hc(0);
    }
  }
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    const int n1m1 = indcs.nx1 + 2*ng - 1;
    par_for("rg_bc_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (i < is || i > ie) {
        b.x1f(m,k,j,i) = 0.0;
        b.x2f(m,k,j,i) = 0.0;
        b.x3f(m,k,j,i) = 0.0;
        if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
        if (j == n2m1) b.x2f(m,k,j+1,i) = 0.0;
        if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
      }
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantFinal

void RedGiantFinal(ParameterInput *pin, Mesh *pm) {
  // every file-scope View must be released before Kokkos::finalize, or it aborts
  lnp_d_ = DvceArray1D<Real>();
  tk_d_ = DvceArray1D<Real>();
  vc_d_ = DvceArray1D<Real>();   // "rg_vc": missing here aborted the run at exit
  ktab_ = DvceArray2D<Real>();
  fconv_ = DvceArray4D<Real>();
  taumlt_ = DvceArray4D<Real>();
  fdiag_ = DvceArray2D<Real>();
  shell_ = DvceArray2D<Real>();
  fmlt1d_ = DvceArray1D<Real>();
  fmean_ = DvceArray2D<Real>();
  klT_ = DvceArray1D<Real>();
  klD_ = DvceArray1D<Real>();
  tau_d_ = DvceArray1D<Real>();  // the seed's tau window (problem/vpert_tau_lo/hi)
  vdb_d_ = DvceArray1D<Real>();  // the bottom sponge's shell means...
  vdb_h_ = HostArray1D<Real>();  // ...and its HOST mirror: Kokkos tracks host
  surf_d_ = DvceArray2D<Real>(); // allocations too, and a default-constructed View
  surf_h_ = HostArray2D<Real>(); // reallocated later carries the empty label "", which
  prof_d_ = DvceArray2D<Real>(); // is what the abort at exit named
  prof_h_ = HostArray2D<Real>();
  return;
}
