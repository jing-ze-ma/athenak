#ifndef UTILS_TWO_STREAM_RT_HPP_
#define UTILS_TWO_STREAM_RT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt.hpp
//! \brief The two-stream radiative-transfer SOLVER: the picket-fence and correlated-k
//! column sweeps, their configuration, their scratch arrays, and the explicit-source
//! limiter.
//!
//! EXTRACTED from src/pgen/deep_hot_jupiter_rt.cpp, where it grew, so that a second
//! problem generator can use it.  The arithmetic is untouched; the extraction is gated on
//! a bitwise-identical problem/ck_dump_file column.  What changed: the enclosing
//! namespace and `inline` on the namespace-scope definitions.
//!
//! HOW IT IS CONFIGURED.  Through the namespace-scope variables below (rt_ck, rt_split,
//! rt_nchain, rt_ck_pcut, rt_de_max, rt_int_at_cut, the dump controls) and through
//! `pm->pgen->hot_jupiter_param`, which despite its name is a generic carrier for the
//! star and planet numbers this solver needs: Teq, grav, ap, omega, Rgas, met and the
//! gravity flags.  A self-luminous object sets Teq = 0.  Tidying that into an explicit
//! config struct is the obvious next step; it was deliberately NOT done in the same
//! change as the move, so that the move could be proved byte-for-byte.
//!
//! WHAT IT NEEDS.  correlated_k for the band opacities and atm_column for the geometry
//! and gravity; both are included here.  The scratch arrays are the rt_*_ptr pointers,
//! allocated by the caller once the mesh exists.

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "utils/eint_from_cons.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"
#include "utils/correlated_k.hpp"
#include "utils/atm_column.hpp"
#include "pgen/pgen_eos_utils.hpp"

namespace two_stream_rt {

//! \fn void par_reduce_clip4 / par_reduce_clip3
//! \brief par_for with an int sum reduction bolted on, flattened exactly the way
//! athena.hpp's par_for flattens its ranges. These exist only so the source limiter can
//! report how many cells it clipped without a second pass over the grid; that is why they
//! live here rather than in athena.hpp.

template <typename Function>
inline void par_reduce_clip4(const std::string &name, const int ml, const int mu,
                             const int kl, const int ku, const int jl, const int ju,
                             const int il, const int iu, int &nclip,
                             const Function &function) {
  const int nk = ku-kl+1, nj = ju-jl+1, ni = iu-il+1;
  const int nkji = nk*nj*ni, nji = nj*ni;
  int cnt = 0;
  Kokkos::parallel_reduce(name,
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (mu-ml+1)*nkji),
  KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    function(m+ml, k+kl, j+jl, i, sum);
  }, cnt);
  nclip += cnt;
}

template <typename Function>
inline void par_reduce_clip3(const std::string &name, const int ml, const int mu,
                             const int kl, const int ku, const int jl, const int ju,
                             int &nclip, const Function &function) {
  const int nk = ku-kl+1, nj = ju-jl+1;
  const int nkj = nk*nj;
  int cnt = 0;
  Kokkos::parallel_reduce(name,
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (mu-ml+1)*nkj),
  KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = idx/nkj;
    int k = (idx - m*nkj)/nj;
    int j = (idx - m*nkj - k*nj);
    function(m+ml, k+kl, j+jl, sum);
  }, cnt);
  nclip += cnt;
}

using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;
using pgen_eos::EintFromDensT;

using correlated_k::CK_NB;
using correlated_k::CK_NG;
using correlated_k::CK_DIFFUSIVITY;
using correlated_k::ck_nq;
using correlated_k::ck_nT;
using correlated_k::ck_nP;
using correlated_k::ck_lT_ptr;
using correlated_k::ck_lP_ptr;
using correlated_k::ck_lk_ptr;
using correlated_k::ck_gw_ptr;
using correlated_k::ck_wl_ptr;
using correlated_k::ck_swf_ptr;
using correlated_k::ck_pf_ptr;
using correlated_k::ck_pf_lTmin;
using correlated_k::ck_pf_idlT;
using correlated_k::ce_nT;
using correlated_k::ce_nP;
using correlated_k::ce_lT_ptr;
using correlated_k::ce_lP_ptr;
using correlated_k::ce_ptr;
using correlated_k::cia_nT_ptr;
using correlated_k::cia_T_ptr;
using correlated_k::cia_k_ptr;
using correlated_k::ray_x_ptr;
using correlated_k::ck_tp_index;
using correlated_k::ck_kappa;
using correlated_k::ck_continuum;
using correlated_k::ck_planck_frac;
using correlated_k::ck_planck_bands;
using atm_column::TGuess;
using atm_column::CSCellAngles;
using atm_column::GravAccAt;
using atm_column::GravPotAt;
using atm_column::TideAccR;
using atm_column::TideAccT;
using atm_column::TideAccP;
using atm_column::EffGravAt;

// Number of (band, quadrature) chains stepped together in the IR sweep. This is the
// instruction-level parallelism the sweep gets; it also fixes the private I_ir_down_c
// footprint at NC radial columns, independent of how many chains are requested.
#ifndef RT_NB
#define RT_NB 4
#endif

// Radial extent of the private intensity column in the GREY split kernel. The
// correlated-k kernel no longer uses this -- it instantiates itself at several sizes and
// dispatches the smallest that fits n1 at run time, because an oversized column is not
// free: at n1 = 68 the chain kernel costs 454 ms with a 72-deep column against 540 with a
// 272-deep one. The grey split path is a test path and keeps the fixed size, guarded.
#ifndef RT_NNC
#define RT_NNC 72
#endif

// Cache the per-layer two-stream coefficients between the two sweeps of the split
// chain kernel. The down-sweep at layer i-1 and the up-sweep at layer i are the SAME
// layer -- both use dtau = tau[i-1]-tau[i] and the same (T,p) proxy -- so they compute
// the identical e0, and alp(down) == gm(up), bet(down) == bet(up). Caching removes one
// expm1, one divide and (with the table on) one k-table lookup and two logs per cell
// per chain, at the price of three more private [NC][NNC] arrays.
#ifndef RT_CACHE
#define RT_CACHE 0
#endif

// Run the correlated-k chain kernel's recurrence in single precision. Default OFF, but
// the
// reason has changed and is worth stating precisely.
//
// Measured while the kernel was still starved on scattered loads, this bought 2.8 %.
// Re-measured once the array layout was fixed it buys 1.42x on rt_chain, 457 -> 325 ms
// per
// 100 cycles: with the memory fed, the FP64 transcendentals really are the next thing in
// the way. That is a 25 % saving on RT.
//
// It is off anyway because 25 % of RT is only 6 % of the run. RT is 41 % of kernel time
// and kernel time is less than wall, so the whole trade is 8.62 -> 8.10 s of integration
// for a 1.7e-4 relative change in the net longwave flux. Judge it on the total, not on
// the
// kernel.
//
// Accuracy if it is ever wanted: the recurrence is a contraction and single precision
// handles it fine. The exposure is the NET flux, which deep down is the difference of
// I_up and I_down when both are close to B. On a real column the outgoing flux at the top
// agrees to 1e-7 and the flux at the cut to 1e-5; the worst level is 1.7e-4.
#ifndef RT_FP32
#define RT_FP32 0
#endif
#if RT_FP32
using RtF = float;
#define RT_EXPM1(x) expm1f(x)
#define RT_EXP(x)   expf(x)
#else
using RtF = Real;
#define RT_EXPM1(x) expm1(x)
#define RT_EXP(x)   exp(x)
#endif

// Column solves ("chains") the RT kernel steps per cell: 4 for the grey picket fence
// (two IR channels x the two Gauss angles), CK_NB*CK_NG*ck_nquad for correlated-k.
// Derived from the scheme, not an input.
inline int rt_nchain = 4;

// problem/rt_split: run the RT as three kernels with the chain-block index promoted to
// a parallel dimension, instead of one kernel that loops over chains serially. Same
// arithmetic, same block-summation order -- see the use site.
inline bool rt_split = false;
// PER-CYCLE RT DIAGNOSTIC (deep_hot_jupiter_rt's problem/diag_gid).  Off by default and
// costing nothing when off: the array is not even allocated.  (m,slot,k,j,i) with
// slot 0 = the source RATE src [code units, erg/cm^3/s] BEFORE the limiter and before
// multiplying by bdt, slot 1 = the de actually added to u0(IEN) after the semi-implicit
// relaxation and the clip, slot 2 = 1.0 if the limiter changed de on this call,
// slot 3 = Ft, the block-summed net flux on the TOP face i+1 after the taublend factor,
// slot 4 = Fb, the same on the BOTTOM face i -- exactly the two values that make up
// src = -(Ft-Fb)/dx1 -- slot 5 = Qs, the block-summed stellar heating (0 below icut),
// slot 6 = Em, the block-summed emission rate with the semi-implicit branch's taublend
// factor, computed for the diagnostic even when the source is applied explicitly.
inline bool rt_diag = false;
inline DvceArray5D<Real> *rt_diag_ptr = nullptr;
// number of cells the RT source limiter clipped on the LAST call (host side, free)
inline int rt_nclip_last = 0;
// (m,k,j,i) face optical depth from the top
inline DvceArray4D<Real> *rt_tau_ptr = nullptr;
inline DvceArray4D<Real> *rt_B_ptr = nullptr;     // (m,k,j,i) Planck function
inline DvceArray4D<Real> *rt_Qv_ptr = nullptr;    // (m,k,j,i) stellar heating rate
inline DvceArray4D<Real> *rt_cf_ptr = nullptr;    // (m,k,j,{gamir1,gamir2,beta})
// (m,blk,k,j,i) IR flux, one slot per block
inline DvceArray5D<Real> *rt_Fb_ptr = nullptr;

// Correlated-k per-cell fields, filled by rt_pre and read by the chain kernel. Laid out
// (m,band,k,j,i) so that the radial index is contiguous for a thread sweeping one band --
// (m,k,j,i,band) would stride every read of the sweep by CK_NB.
inline DvceArray5D<Real> *rt_kc_ptr = nullptr;    // (m,b,k,j,i) continuum kappa [cm^2/g]
inline DvceArray5D<Real> *rt_Bb_ptr = nullptr;    // (m,b,k,j,i) band Planck intensity
inline DvceArray4D<Real> *rt_T_ptr = nullptr;     // (m,k,j,i) temperature [K]
inline DvceArray4D<Real> *rt_pb_ptr = nullptr;    // (m,k,j,i) pressure [bar]
// (m,k,j,i) continuous k-table index in T
inline DvceArray4D<Real> *rt_xT_ptr = nullptr;
// (m,k,j,i) continuous k-table index in p
inline DvceArray4D<Real> *rt_xP_ptr = nullptr;
// (m,k,j) deepest cell doing correlated-k
inline DvceArray3D<int>  *rt_icut_ptr = nullptr;
// (m,blk,k,j,i) stellar heating, one per block
inline DvceArray5D<Real> *rt_Qb_ptr = nullptr;
// (m,blk,i,k,j) the LOCAL thermal emission rate per unit volume [erg/cm^3/s] that each
// chain block puts into cell i.  Summed over blocks this is 4 sigma kappa_P rho T^4, the
// only part of the radiative source that depends on the cell's OWN temperature.  It is
// what turns the explicit update into a stable exponential one; see the apply kernel.
inline DvceArray5D<Real> *rt_Em_ptr = nullptr;
inline DvceArray5D<Real> *rt_Src_ptr = nullptr;   // per-cell net source, formed DIRECTLY

// --- correlated-k (Lee/Exo-FMS premixed tables, Kataria+2013 11-band grid) --------
// problem/rt_ck turns it on; problem/ck_table is the path to the premixed table and
// problem/ck_pcut_bar the pressure below which correlated-k is used at all (deeper than
// that the atmosphere is optically thick -- grey tau ~ 9e3 at 10 bar for this setup --
// and the interior is far outside any molecular table, up to 12000 K).
//
// Table layout and the traps it carries are documented in data/exo_fms_ck/PROVENANCE.md.
// The two that matter here: the data block runs the band index BACKWARDS, and kappa is
// cgs cm^2/g, which is what this code already works in.
inline bool rt_ck = false;
// problem/ck_star_teff: host effective temperature. If > 0 the stellar spectrum is a
// blackbody at this temperature, which needs no external data since only the SHAPE is
// used. If <= 0, problem/ck_swflux is read instead.
//
// Default 6000 K, the middle of the range the ultra-hot Jupiter GCM literature actually
// uses: Tan et al. (2024, MNRAS 528, 1016) grid over host T_eff of 5500, 6000 and 6500 K,
// and Parmentier et al. (2018, A&A 617, A110) model WASP-121b at 6460 K. That range also
// keeps the flux falling bluer than the grid's 0.26 um edge down to 1-3 %, which is what
// makes the 11-band Kataria structure defensible in the first place -- an A-type host
// would put 18-32 % outside it.
inline Real rt_star_teff = 6000.0;
// problem/ck_dump_file: write one column's RT solution -- level pressures, temperatures,
// net longwave flux and stellar heating -- straight out of the production kernel, once,
// at
// the first RT call. This exists so the kernel ITSELF can be compared against an external
// code on the same profile, rather than a transcription of it. See
// docs/correlated_k_rt.md.
inline std::string rt_dump_file = "";
inline int rt_dump_m = 0;
inline int rt_dump_j = -1;
inline int rt_dump_k = -1;
inline bool rt_dump_done = false;
inline Real rt_ck_pcut = 10.0;                    // bar
// problem/rt_grey: the GREY two-stream, built on the correlated-k machinery rather than
// on the old picket-fence path.  One band, one opacity, taken from the conduction
// module's own kappa table (<hydro|mhd>/rad_kappa_src) so that the two-stream and the
// radiative diffusion it blends with cannot disagree about the opacity, and the full
// Planck function sigma T^4/pi as the source.
//
// WHY IT EXISTS, and why it is not the old grey path.  The old one is the Parmentier
// picket fence with the Freedman fit: two IR channels whose split (gamma_1, gamma_2,
// beta) is a fit to IRRADIATED giant planets, and an analytic opacity that is 1.7x off
// the stellar table on a cool giant.  This one shares every structural fix the
// correlated-k path got -- the layer-integrated two-stream coefficients, the flux
// DIFFERENCE deposition, the tau blend, the icut handover, the semi-implicit source --
// and differs from it in exactly one way: emission and absorption use the SAME opacity,
// so kappa cancels out of the local radiative balance and the equilibrium temperature of
// a thin layer is a fixed point no matter how the opacity varies.  A correlated-k layer
// weights its absorption by the incident spectrum and its emission by its own, so a
// layer that drifts cold can absorb less than it emits and keep drifting; that is the red
// giant top-layer runaway.  Grey cannot do that.  It also cannot represent line
// blanketing, which is the reason the band solver exists -- so this is the control, and
// the diagnosis, not automatically the production choice.
inline bool rt_grey = false;

// --- read-only access to the band solver's radial face flux, for diagnostics ---------
// rt_Fb holds the NET longwave flux on the radial faces, (m, block, i, k, j), in code
// flux units (pressure x velocity), for faces i = is..ie+1.  It is a namespace-scope
// pointer and lives for the run once the split band path has allocated it, so a problem
// generator can sum it without re-running the solver.  rt_face_flux_ready() is the guard
// every caller must use: the monolithic and non-band paths never fill it.
inline bool rt_face_flux_ready() {
  return (rt_Fb_ptr != nullptr) && (rt_ck || rt_grey);
}
inline int rt_face_nblk() {
  return (rt_Fb_ptr != nullptr) ? rt_Fb_ptr->extent_int(1) : 0;
}
inline DvceArray5D<Real> rt_face_flux() { return *rt_Fb_ptr; }
// (m,k,j) deepest cell the band solver integrates; the face below it is where the
// interior flux is handed in from the diffusion operator.
inline DvceArray3D<int> rt_cut_index() { return *rt_icut_ptr; }
// problem/rt_de_max: the cap in LimitRTSource, as a fraction of the cell's internal
// energy per RT application. Applies to every EXPLICIT radiative update -- grey and
// correlated-k, split and monolithic. Set <= 0 to disable the limiter entirely.
inline Real rt_de_max = 0.5;
// problem/rt_semi_lin: recover the OLD semi-implicit step, which linearized the emission
// about the current state and relaxed at lambda = 4E/e.  See the long note in rt_apply:
// that form bounds cooling but leaves heating explicit and unbounded, so a cold optically
// thin cell overshoots its equilibrium badly.  Kept only to reproduce runs made before
// the fix.  DEFAULT TRUE, i.e. the old linearization, so that every problem generator
// sharing this header reproduces its pre-fix runs bit for bit; the red-giant inputs set
// it false to select the new equilibrium-relaxation step.
inline bool rt_semi_lin = true;
// problem/rt_use_cons: take the cell's internal energy and density from the CONSERVED
// state u0 instead of from w0.  w0 is the previous stage's ConToPrim output, and by the
// time the user source function runs RKUpdate, the explicit source terms and the
// implicit radial conduction have all moved u0.  In a smooth cell the two agree to
// O(dt); in a runaway cell the hydro step changes e by ~100 % per stage, so the
// equilibrium, the Newton iterate, the positivity guard and LimitRTSource are all
// evaluated on a number that no longer exists while the step is applied to u0 anyway --
// which is how a cell is driven to a negative internal energy, repaired by the floor,
// and handed to the next Riemann solve as a 0.1 K cell beside a 1e4 K one.  Default
// FALSE so every existing input reproduces bit-for-bit; the red-giant runs set it true.
inline bool rt_use_cons = false;
// problem/rt_bface: the emissivity-weighted far-endpoint Planck source (see BFace).
// It is a red-giant fix -- it exists because the corona/star join put a cell's emission
// on a neighbour's Planck function 1e9 times its own -- and this header is shared with
// solar_convection and the hot-Jupiter problems, which have no such join.  DEFAULT FALSE
// therefore returns the old endpoint b_far bit for bit at every one of the seven call
// sites; the red-giant inputs set it true.
inline bool rt_bface = false;
// problem/rt_explicit: take the source EXPLICITLY, de = src*bdt, and skip the whole
// equilibrium block -- no closed form, no linearization, no Newton.  With the direct
// (cancellation-free) source the explicit step's stability limit is the local radiative
// time e/|src|, which is ~1e3 s at tau 1-10 against a 30 s timestep, so this is a usable
// control and not just a diagnostic: it removes the semi-implicit update from the picture
// entirely for the R9 death bisection.  rt_de_max still applies if it is set (> 0).
inline bool rt_explicit = false;
// how many cells the Newton loop had to be rescued from a non-positive internal energy
// (see rt_apply).  A device counter, read back where the clip count is reported.
inline DvceArray1D<int> *rt_efix_ptr = nullptr;
// problem/rt_newton: refine the semi-implicit step with a Newton solve of the exact
// backward-Euler balance, seeded by the closed form.  The closed form already lands on
// the right equilibrium under e ~ T; this drops that assumption and uses the EOS's own
// T(e) and c_v(e), which matters wherever H2 or H is partly dissociated.  Costs one
// Temperature() and one SpecificHeatCv() per iteration, and from that seed it is
// normally one or two.  General EOS only; with an ideal gas e ~ T is exact and the
// closed form is already the answer.  DEFAULT FALSE so shared pgens are unchanged; the
// red-giant inputs set it true.
inline bool rt_newton = false;
// problem/rt_rescue_eq: when the Newton step would leave e <= 0, land the cell on the
// radiative equilibrium it actually sees -- deq = ei((A/Em)^(1/4) - 1), i.e. Em(T_eq) = A
// with A = src + Em the absorption, which in the thin limit is (pi/mu) kappa rho
// (I_up + I_dn) so that sigma T_eq^4 = pi (I_up + I_dn)/2 up to the diffusivity factor --
// instead of the unconditional 99.9 % drop.  Floored at that same 99.9 %, so it can only
// make the rescue less violent.  DEFAULT FALSE = the old 99.9 % policy, so shared pgens
// are unchanged; the red-giant inputs set it true.
inline bool rt_rescue_eq = false;
// problem/rt_src_direct: form the per-cell radiative source DIRECTLY as absorption minus
// emission during the sweeps, instead of as the difference of the two face fluxes.  The
// two are the same number algebraically -- the stream update across a layer is
//     I_out = (1 - e0) I_in + (alp B_far + bet B_near),   alp + bet = e0,
// so the flux change over the layer is e0 (I_in - B) exactly, and summing those over a
// column telescopes to F_top - F_bot -- but not numerically.  In a transparent cell the
// two face fluxes agree to every digit and their difference is round-off: measured
// 1.7e-11 relative on 5.67e9, i.e. 5.8e-12 erg/cm^3/s against a cell whose whole
// internal energy is 6.3e-12.  That noise, of either sign, is what heated the ambient
// medium to 3000 K in two steps and then drove it to zero; no limiter or opacity floor
// touches it, because it is not a physical term at all.  Formed directly the source is
// O(dtau) with no large numbers cancelling, and it vanishes as the cell goes transparent
// the way the physics says it must.  Off reproduces the flux-difference form bit for bit,
// and is the DEFAULT for that reason; the red-giant inputs set it true.
inline bool rt_src_direct = false;
// problem/rt_top_clamp: the band solver's top slot i = ie+1 takes its fluid state from
// the top ACTIVE cell ie instead of from the hydro ghost (see rt_pre_tp).  It decouples
// the whole column from whatever the outer boundary condition put in the ghost, which is
// what a poisoned ghost needs, but it changes numbers in any run whose top ghost is a
// legitimate state.  DEFAULT FALSE = read the ghost exactly as before; the red-giant
// inputs set it true.
inline bool rt_top_clamp = false;
// problem/rt_semi_implicit: apply the split two-stream source semi-implicitly, relaxing
// the cell toward radiative equilibrium instead of stepping explicitly (see the apply
// kernel).  Default true, which is the behaviour since 048dff30.  The flag exists
// because the semi-implicit step is NOT answer-preserving in general: on the cubed
// sphere deep_hot_jupiter 09-07 configuration it moves the kinetic energy by 3 % in 20
// cycles (A/B 2026-09-09).  Runs that must reproduce pre-048dff30 results set it false,
// which restores the plain explicit de = src*bdt.
// NOTE: this and problem/rt_explicit above are two switches onto the same branch, kept
// separate because their defaults differ.  The semi-implicit block runs only when
// rt_semi_implicit is true AND rt_explicit is false, so each default (true / false)
// reproduces its own side's behaviour and either flag alone selects the explicit step.
inline bool rt_semi_implicit = true;
// problem/ck_int_at_cut: deliver the planet's internal flux sigma T_int^4 as an extra
// upward source at the correlated-k cut (the historical behaviour, true). Set false when
// the layers below the cut carry it themselves -- <mhd|hydro>/isotropic_conduction =
// radiative with rad_flux_inner at the bottom wall -- so the cut's upward intensity is
// just the thermalised Planck function and nothing is counted twice.
inline bool rt_int_at_cut = true;
// problem/rt_top_re: what the unresolved column ABOVE the domain sends back down.
// false (historical) makes it radiate at the ghost cell's own temperature. That is safe
// only while the ghost is pinned to something outside the solution: with an open outer
// boundary the ghost IS the top active cell, the column hands back exactly what that cell
// emitted, and the net flux through the top face goes to zero as the column thickens --
// the atmosphere isothermalises and the object stops radiating (measured on the red
// giant: emergent flux 0.3 % of L, top density 22x). true treats the column as a slab in
// RADIATIVE EQUILIBRIUM instead: it absorbs (1-e^-dtau) I_up from below and re-emits half
// up to space and half back down, so its source is I_up/2 whatever its optical depth and
// the face always keeps at least half of its outgoing flux. Grey path only so far.
// NOTE, and the reason for the startup warning below: the back-radiation is delivered
// as (1 - exp(-dtau/mu)) bsrc with dtau = kappa(ghost) p/g.  Above rad_kappa_rmax the
// ghost's opacity is rad_kappa_above, so with rad_kappa_above = 0 -- the exact setting
// every corona run uses -- dtau is 0 and this boundary hands back EXACTLY ZERO however
// carefully bsrc was probed.  Measured in V8_corona: I_dn at the top face of the active
// column is 0.000000e+00 at every cycle.  Behaviour unchanged; it is only documented.
inline bool rt_top_re = false;
// SELF-LUMINOUS objects: > 0 uses this internal temperature directly instead of the
// Thorngren+2019 T_int(T_eq) relation, which is a fit for IRRADIATED giant planets and
// floors at 100 K.  A star or a brown dwarf sets it to its own T_eff (and Teq = 0, which
// zeroes the stellar sweep).  0 keeps the historical behaviour exactly.
inline Real rt_tint_override = 0.0;
// problem/ad_dump_file: write the INITIAL (p, T) profile with its actual and adiabatic
// logarithmic gradients, once, before and after adjust_ad_pT_arr. A diagnostic for
// whether the starting atmosphere is convectively unstable where the adjustment did not
// reach.
// adjust_ad_pT_arr scans from the bottom, breaks at the FIRST crossing and enforces a
// single adiabat below it, so with a general EOS -- where grad_ad dips at H2 dissociation
// as well as at H ionization -- a second unstable zone could in principle survive it.
// Measured on the deep hot Jupiter setup it does not: after the adjustment there are zero
// super-adiabatic levels, H2 on or off, because there is only one crossing to begin with.
inline std::string ad_dump_file = "";
// problem/rt_apply_debug: print the energy balance the applied source actually sees, for
// the top rt_apply_debug_n cells of the (rt_dump_m, rt_dump_k, rt_dump_j) column, on the
// first N calls.  The column dump says what the SWEEP produced; this says what the CELL
// then did with it, which is where a thin top layer that will not sit at its equilibrium
// temperature has to be caught.
inline int rt_apply_debug = 0;
inline int rt_apply_debug_n = 8;
// problem/rt_cell_report: the radiative-balance report for the cells that decide the
// top of the active column.  Prints, once per RT call, the FIRST cell that trips the
// Newton positivity rescue -- its location, state, opacity, tau to the top, the up and
// down streams at its two faces, the absorption and emission terms, e_eq, dt and every
// Newton iterate -- and, every rt_report_every cycles, the same for the fixed cell that
// contains rt_report_r on the (ks, js) column.  Off by default: a default run neither
// allocates the two stream arrays nor prints.
inline bool rt_cell_report = false;
inline Real rt_report_r = 3.887e12;               // problem/rt_report_r [cm]
inline int rt_report_every = 100;                 // problem/rt_report_every [cycles]
inline DvceArray4D<Real> *rt_idn_ptr = nullptr;   // downward stream at face i
inline DvceArray4D<Real> *rt_iup_ptr = nullptr;   // upward stream at face i
inline bool rt_srclim_warned = false;             // the one-time warning has been issued
// problem/nan_report: catch the cell whose conserved energy the GREY apply makes
// non-finite or non-positive, IN the kernel, with the inputs that produced it.  Off by
// default, so a default run is bit-identical.  Pointers, not Views, for the same reason
// the pgen's guards use pointers: a file-scope View outlives Kokkos::finalize.
inline bool rt_nan_report = false;
inline DvceArray1D<int> *rt_nanrep_cnt = nullptr;
inline DvceArray1D<Real> *rt_nanrep_rec = nullptr;
inline int rt_nanrep_lines = 0;                   // printed so far on this rank
inline int rt_nanrep_maxlines = 400;              // then stay quiet (the run continues)

// problem/rt_use_cons NaN GUARD.  EintFromCons returns whatever E - KE gives, INCLUDING
// a negative number: in a runaway cell the total energy is kinetic dominated by many
// orders of magnitude and the difference is catastrophic cancellation.  A non-positive
// e makes TempKelvin/PresTempFromEint return NaN, and one NaN Planck function poisons
// the whole column (measured: T(i-1) = nan, B(i-1) = nan, I_dn = nan ahead of the
// red-giant dt collapse).  eiN clamps such a cell to e(rho, tfloor) -- the same floor
// ConsToPrim would apply -- so the sweep never sees a NaN, and counts the clamps here.
// Active only under rt_use_cons; with it off eiN returns w0 and nothing changes.
inline DvceArray1D<int> *rt_eiclamp_cnt = nullptr;
inline bool rt_eiclamp_warned = false;

// The SAME hazard one level up: a (p, T) pair the EOS could not form -- a non-positive
// or non-finite conserved DENSITY in the cell a slot reads, or a solve that did not
// converge -- gives a NaN Planck function and a NaN top-of-column optical depth, and a
// single NaN poisons the WHOLE column below it.  The down-sweep carries the stream as
// I(i) = (1-e0)*I(i+1) + ..., and at large dtau (1-e0) underflows to EXACTLY zero, so
// 0*NaN = NaN survives to arbitrary depth: measured in I8_vceil, 348 rescued cells at
// tau_to_top = 75-90 with I_dn = -nan while their own B, T and kappa were finite.  The
// two guards below are the density and the state analogues of the eiN clamp, counted in
// one counter, and they are no-ops on any state the EOS returned cleanly.
inline DvceArray1D<int> *rt_stclamp_cnt = nullptr;
inline bool rt_stclamp_warned = false;

//----------------------------------------------------------------------------------------
//! \fn bool RTBadState
//! \brief true when a (p, T) pair cannot be used to form a Planck function or a
//! top-of-column optical depth: NaN, infinite, or non-positive.
KOKKOS_INLINE_FUNCTION
bool RTBadState(const Real p, const Real t) {
  return !(t > 0.0) || !(p > 0.0) || !Kokkos::isfinite(t) || !Kokkos::isfinite(p);
}

//----------------------------------------------------------------------------------------
//! \fn Real RTTopDtau
//! \brief the optical depth of the unresolved hydrostatic column above the domain,
//! kappa p / g_eff, with the g_eff = 0 case resolved.
//!
//! g_eff > 0 returns the old expression BIT FOR BIT.  Where the effective gravity is
//! zero (a plane-parallel run that never set one) the old form was kappa*p/0, which is
//! +inf for an opaque ghost -- an opaque lid, harmless -- but 0/0 = NaN as soon as the
//! ghost's opacity is zero, which is exactly what the density gate makes it above the
//! star.  The limit is taken here instead: no opacity above the domain means no column.
KOKKOS_INLINE_FUNCTION
Real RTTopDtau(const Real kap, const Real pcgs, const Real geff) {
  const Real kp = kap*pcgs;
  if (geff > 0.0) return kp/geff;           // the old expression, unchanged
  return (kp > 0.0) ? 1.0e30 : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RTSourceLimiterWarn
//! \brief say ONCE, from rank 0, that LimitRTSource has clipped cells.
//!
//! Once, not per call: a run that trips this trips it every cycle, and the message is
//! about the configuration, not about the individual step.
inline void RTSourceLimiterWarn(const int nclip) {
  if (nclip <= 0 || rt_srclim_warned) return;
  rt_srclim_warned = true;
  if (global_variable::my_rank == 0) {
    std::cout << "### WARNING in deep_hot_jupiter_rt: the explicit radiative source was "
              << "clipped in " << nclip << " cell(s) by problem/rt_de_max = " << rt_de_max
              << ".\n    The radiative time e/|src| is shorter than the timestep there, "
              << "so the radiation is outside\n    the regime this operator-split scheme "
              << "is valid in. The usual cause is a pressure or\n    density floor "
              << "pinning the top of the atmosphere; see docs/correlated_k_rt.md.\n"
              << "    Reported once per run." << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real LimitRTSource
//! \brief cap one explicit radiative energy update at a fraction of the cell's internal
//! energy.
//!
//! WHY THIS EXISTS. The radiation is operator split and applied EXPLICITLY at the
//! hydrodynamic timestep: u0(IEN) += src*bdt, with src the net flux divergence plus the
//! stellar heating. That is stable only while the local radiative time e/|src| exceeds
//! bdt, and nothing in the code enforces it. Measured on the column that first blew up in
//! the ideal-gas correlated-k run (see docs/correlated_k_rt.md): a healthy column has
//! e/|src| between 86 and 570 timesteps, but once the pressure floor pins the top of the
//! atmosphere at p = pfloor -- which fixes e = pfloor/(gamma-1), 2.1 erg/cm^3 there -- an
//! unremarkable flux divergence of 29 erg/cm^3/s gives e/|src| = 0.014 timesteps. The
//! explicit update then drives e far negative, the floor rescues it, and the next step
//! overshoots harder: adjacent cells at 650, 2081 and 9806 K, and a NaN a few cycles
//! later.
//!
//! WHY A HARD CLAMP, NOT A SMOOTH ONE. A smooth limiter perturbs every cell a little; a
//! clamp is the IDENTITY wherever |de| < de_max*e, so it cannot move an answer in any
//! regime where the explicit step was legitimate in the first place. The correlated-k
//! fluxes are validated against Exo-FMS and that validation has to survive this. At the
//! default de_max = 0.5 the clamp is 40x looser than the worst healthy cell measured
//! above and 140x tighter than the runaway, so it separates the two cleanly.
//!
//! This bounds the damage; it does not make the step accurate. A run that trips it is
//! reporting that its floors, or its timestep, put the radiation outside the regime the
//! scheme is valid in -- which is why tripping it is warned about exactly once.
// The emissivity ratio below which BFace stops trusting the neighbour's Planck function
// as this layer's far source endpoint.  1/10: at V8_corona's join the ratio is exactly 0
// (rad_kappa_above = 0) and in a smooth stellar column it never falls below ~0.5, so this
// separates the two by 5x on the smooth side and infinitely on the pathological one.
#define RT_BFACE_R 0.1

//----------------------------------------------------------------------------------------
//! \fn Real BFace
//! \brief the Planck source at a layer's FAR endpoint, weighted by emitting matter.
//!
//! WHY THIS EXISTS.  Both sweeps take the source function of a layer as B interpolated
//! linearly between the two cell CENTRES -- alp*B_far + bet*B_near -- and Em averages the
//! same two.  Nothing in that asks whether the far cell radiates at all.  Across the join
//! at rad_kappa_rmax it does not: measured in V8_corona at t = 0, the last active cell
//! (i = 439, T = 3363.6 K, B = 2.310e9) sits directly under the first corona cell
//! (T = 6.0e5 K, kappa = 0, B = 2.339e18), so its emission was formed with a Planck
//! function 1.01e9 times its own: 3.60e-8 erg/cm^3/s of emission against an absorption of
//! 5.50e-17 and an internal energy of 7.46e-10, i.e. 1480 times its own energy removed in
//! one 30.6 s step.  All 192 cells of that shell hit the temperature floor on the FIRST
//! RT call and the corona then accreted onto the cold sink.
//!
//! WHAT IT DOES.  It is SURGICAL: while the far cell's emissivity kappa*rho is within a
//! factor RT_BFACE_R of this layer's, the old endpoint B_far is returned unchanged, so
//! every column whose opacity varies smoothly is BIT-IDENTICAL to the pre-fix code
//! (measured on the V7 grey column: all 39 active faces identical).  Below that ratio the
//! endpoint is blended into the emissivity-weighted face source
//!     S = (j_own B_own + j_far B_far)/(j_own + j_far),   j = kappa rho,
//! with a weight that reaches S exactly as the far cell stops radiating -- so an inert
//! neighbour contributes B_own, i.e. the layer emits with its OWN Planck function.  The
//! blend is continuous at the threshold (S is entered with weight 0 there), so a cell
//! cannot jump between the two forms.
KOKKOS_INLINE_FUNCTION
Real BFace(const Real k_own, const Real k_far, const Real b_own, const Real b_far,
           const bool on) {
  if (!on) return b_far;                  // problem/rt_bface off: the pre-fix endpoint
  const Real kt = RT_BFACE_R*k_own;
  if (k_far >= kt) return b_far;          // the old expression, bit for bit
  if (!(k_own > 0.0)) return b_far;       // this layer does not emit at all
  if (!(k_far > 0.0)) return b_own;       // the far cell is radiatively inert
  const Real w = k_far/kt;                // 1 at the threshold, 0 at k_far = 0
  const Real sem = (k_own*b_own + k_far*b_far)/(k_own + k_far);
  return w*b_far + (1.0 - w)*sem;
}

KOKKOS_INLINE_FUNCTION
Real LimitRTSource(const Real de, const Real eint, const Real de_max) {
  const Real cap = de_max*eint;
  return (de > cap) ? cap : ((de < -cap) ? -cap : de);
}

KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A) {
  // Parmentier+2015
  Real X = log10(Teff);
  Real g = gg*0.01;
  Real a,b;
  if (Teff < 250.0) {
    a = -0.335*pow(g,0.07);
    b = 0.0;
  } else if (Teff < 750.0) {
    a = -0.335*pow(g,0.07) + 2.149*pow(g,0.135);
    b = -0.896*pow(g,0.135);
  } else if (Teff < 1250.0) {
    a = -0.335*pow(g,0.07) - 0.428*pow(g,0.135);
    b = 0.0;
  } else {
    a = 16.947 - 3.174*pow(g,0.07) - 4.051*pow(g,0.135);
    b = -5.472 + 0.917*pow(g,0.07) + 1.170*pow(g,0.135);
  }
  Real log10A = a + b*X;
  A = pow(10.0,log10A);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2,
                            Real &gamv3, Real &beta, Real &gamir1, Real &gamir2) {
  // Parmentier & Giollot 2014; Parmentier+2015; Roth+2024
  Real X = log10(Teff);
  Real a3, a2, a1, b3, b2, b1, ab, bb;

  Real ap = -2.36;
  Real bp = 13.92;
  Real cp = -19.38;
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ap = -12.45;
    bp = 82.25;
    cp = -134.42;
  }

  if (Teff < 2000.0) {
    ab = 0.84;
    bb = 0.0;
  } else {
    ab = 6.21;
    bb = -1.63;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ab = 3.0;
    bb = -0.69;
  }

  if (Teff < 200.0) {
    a3 = -3.03;
    b3 = -0.2;
    a2 = -7.37;
    b2 = 2.53;
    a1 = -5.51;
    b1 = 2.48;
  } else if (Teff < 300.0) {
    a3 = -13.87;
    b3 = 4.51;
    a2 = 13.99;
    b2 = -6.75;
    a1 = 1.23;
    b1 = -0.45;
  } else if (Teff < 600.0) {
    a3 = -11.95;
    b3 = 3.74;
    a2 = -15.18;
    b2 = 5.02;
    a1 = 8.65;
    b1 = -3.45;
  } else if (Teff < 1400.0) {
    a3 = -6.97;
    b3 = 1.94;
    a2 = -10.41;
    b2 = 3.31;
    a1 = -12.96;
    b1 = 4.33;
  } else if (Teff < 2000.0) {
    a3 = -3.65;
    b3 = 0.89;
    a2 = -19.95;
    b2 = 6.34;
    a1 = -23.75;
    b1 = 7.76;
  } else {
    a3 = -6.02;
    b3 = 1.61;
    a2 = 13.56;
    b2 = -3.81;
    a1 = 12.65;
    b1 = -3.27;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    if (Teff < 2000.0) {
      a3 = 0.02;
      b3 = -0.28;
      a2 = 6.96;
      b2 = -2.21;
      a1 = -1.68;
      b1 = 0.75;
    } else {
      a3 = -16.54;
      b3 = 4.74;
      a2 = -2.4;
      b2 = 0.62;
      a1 = 10.37;
      b1 = -2.91;
    }
  }
  Real log10gamv1 = a1 + b1*X;
  Real log10gamv2 = a2 + b2*X;
  Real log10gamv3 = a3 + b3*X;
  Real log10gamp = ap*SQR(X) + bp*X + cp;
  beta = ab + bb*X;

  gamv1 = pow(10.0,log10gamv1);
  gamv2 = pow(10.0,log10gamv2);
  gamv3 = pow(10.0,log10gamv3);
  Real gamp = pow(10.0,log10gamp);
  Real dum = (gamp-1.0)/(2.0*beta*(1.0-beta));
  Real R = 1.0 + dum + sqrt(SQR(dum)+dum);
  gamir1 = beta + R - beta*R;
  gamir2 = gamir1/R;

  return;
}

KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr) {
  // Freedman+2014, shared with the conduction module (utils/rosseland.hpp)
  kapr = RosselandFreedman2014(T, p, met);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint) {
  // Thorngren+2019 +Erratum
  Real boltz_sigma = 5.6704e-5;
  Real F = 4.0*boltz_sigma*SQR(SQR(Teq));
  Tint = 0.39*Teq*exp(-SQR(log10(F)-9.0-0.14)/1.095);
  Tint = (Tint < 100.0) ? 100.0 : Tint;
  return;
}

inline void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    // Noti+2023; Lee+2021

    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;

    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    // problem/rt_use_cons: where the solver reads the cell's thermodynamic state from.
    // eiN/rhoN are the ONLY way this routine touches e and rho below, so the switch
    // reaches the precompute, the sweeps, the equilibrium, the Newton, the positivity
    // guard, the rescue floor and LimitRTSource together -- there is no half-converted
    // path.  With the switch off they return exactly w0(...), so the default is
    // bit-identical by construction.
    const bool usecons_ = rt_use_cons;
    const bool bface_on = rt_bface;
    const bool cs_uc_ = pm->use_cubed_sphere;
    auto cosc_uc_ = pmbp->pcoord->cos_cell;
    const bool etg_uc_ = (pmbp->phydro != nullptr) ? pmbp->phydro->use_etotgrav
                       : ((pmbp->pmhd != nullptr) ? pmbp->pmhd->use_etotgrav : false);
    DvceArray4D<Real> phicc_uc_("rt_phi_dummy", 1, 1, 1, 1);
    if (pmbp->phydro != nullptr) {
      phicc_uc_ = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      phicc_uc_ = pmbp->pmhd->phicc0;
    }
    auto u0_uc_ = u0;
    auto w0_uc_ = w0;
    // see the note on rt_eiclamp_cnt: e <= 0 (or NaN) out of the conserved state is
    // clamped to e(rho, tfloor) before any temperature or Planck function is formed
    if (rt_eiclamp_cnt == nullptr) {
      rt_eiclamp_cnt = new DvceArray1D<int>("rt_eiclamp", 1);
      Kokkos::deep_copy(*rt_eiclamp_cnt, 0);
    }
    auto eicl_g = *rt_eiclamp_cnt;
    auto eos_uc_ = eos;
    auto eiN = [=] (const int m, const int k, const int j, const int i) {
      if (!usecons_) return w0_uc_(m,IEN,k,j,i);
      const Real ei_uc = EintFromCons(u0_uc_, m, k, j, i,
                                      cs_uc_ ? cosc_uc_(m,k,j) : 0.0, cs_uc_,
                                      etg_uc_, etg_uc_ ? phicc_uc_(m,k,j,i) : 0.0);
      if (ei_uc > 0.0) return ei_uc;           // false for NaN too, which is the point
      Kokkos::atomic_fetch_add(&eicl_g(0), 1);
      Real ei_fl = eos_uc_.EnergyFromTemperature(u0_uc_(m,IDN,k,j,i), eos_uc_.tfloor);
      if (!(ei_fl > 0.0)) ei_fl = 1.0e-300;
      return ei_fl;
    };
    if (usecons_ && !rt_eiclamp_warned) {
      auto eicl_h = Kokkos::create_mirror_view(eicl_g);
      Kokkos::deep_copy(eicl_h, eicl_g);
      if (eicl_h(0) > 0) {
        rt_eiclamp_warned = true;
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING in two_stream_rt: rt_use_cons gave a non-positive "
                    << "internal energy in " << eicl_h(0) << " cell read(s); clamped to "
                    << "e(rho,tfloor) before the temperature/Planck evaluation. "
                    << "Reported once." << std::endl;
        }
      }
    }
    // see rt_stclamp_cnt: a non-positive or non-finite density, and a (p,T) the EOS
    // could not form from it, are clamped before any Planck function or optical depth
    if (rt_stclamp_cnt == nullptr) {
      rt_stclamp_cnt = new DvceArray1D<int>("rt_stclamp", 1);
      Kokkos::deep_copy(*rt_stclamp_cnt, 0);
    }
    auto stcl_g = *rt_stclamp_cnt;
    if (!rt_stclamp_warned) {
      auto stcl_h = Kokkos::create_mirror_view(stcl_g);
      Kokkos::deep_copy(stcl_h, stcl_g);
      if (stcl_h(0) > 0) {
        rt_stclamp_warned = true;
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING in two_stream_rt: a non-finite or non-positive "
                    << "density/(p,T) state was read in " << stcl_h(0) << " cell(s); "
                    << "clamped to the floor state before the Planck function and the "
                    << "optical depth, so the column below it is not poisoned. "
                    << "Reported once." << std::endl;
        }
      }
    }
    const Real dfl_uc_ = eos.dfloor;
    auto rhoN = [=] (const int m, const int k, const int j, const int i) {
      const Real d = usecons_ ? u0_uc_(m,IDN,k,j,i) : w0_uc_(m,IDN,k,j,i);
      if (d > 0.0) return d;                 // false for NaN too, as in eiN
      Kokkos::atomic_fetch_add(&stcl_g(0), 1);
      return dfl_uc_;
    };
    {
      static bool uc_announced = false;
      if (!uc_announced && global_variable::my_rank == 0) {
        uc_announced = true;
        std::cout << "### two_stream_rt: thermodynamic state read from "
                  << (usecons_ ? "u0 (CONSERVED, problem/rt_use_cons = true)"
                               : "w0 (previous ConToPrim; rt_use_cons = false)")
                  << std::endl;
      }
    }

    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;

    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;

//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real met = 0.0;

//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    const bool tide = pm->pgen->hot_jupiter_param.stellar_tide;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    Real met = pm->pgen->hot_jupiter_param.met;

    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;

    Real Tirr = Teq*sqrt(2);
    Real Tirr4 = SQR(SQR(Tirr));
    Real Fstar = boltz_sigma*Tirr4;
    Real Tint;
    if (rt_tint_override > 0.0) {
      Tint = rt_tint_override;
    } else {
      get_Tint(Teq, Tint);
    }
    Real Tint4 = SQR(SQR(Tint));
    bool int_at_cut = rt_int_at_cut;
    const bool top_re = rt_top_re;
    Real Iint = boltz_sigma/M_PI*Tint4;

    const int nchain_rt = rt_nchain;

    Real mug[2];
    Real wg[2];
    mug[0] = 0.21132487;
    mug[1] = 0.78867513;
    wg[0] = 0.5;
    wg[1] = 0.5;

    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);


    // ================================================================================
    // CHAIN-PARALLEL SPLIT PATH  (problem/rt_split)
    // --------------------------------------------------------------------------------
    // The monolithic kernel below is parallel over (m,k,j) only -- 8192 columns here,
    // 128 wavefronts on 1216 SIMDs -- so it runs at MeanOccupancyPerCU ~0.5 and
    // VALUBusy 4-6 %, and the frequency chains are a SERIAL loop inside each thread.
    // That is why the cost is exactly linear in the chain count: it is exp() latency on
    // a dependency chain with no other wave to hide it, not a throughput limit.
    //
    // This path splits the same arithmetic into three kernels so the chain-block index
    // becomes a parallel dimension:
    //   A "rt_pre"    (m,k,j)       chain-independent: tau, B, Q_v and the picket-fence
    //                               coefficients, written to global arrays.
    //   B "rt_chain"  (m,blk,k,j)   one thread per column per block of RT_NB chains.
    //                               nblk x more threads than the monolithic kernel.
    //   C "rt_apply"  (m,k,j,i)     sum the per-block fluxes in block order and apply.
    //
    // Kernel C sums the blocks in the SAME order the serial loop accumulated them,
    // starting from the same 0.0, so with one block (production, 4 chains) the result is
    // bitwise identical to the monolithic path. With the harness active the extra blocks
    // contribute exactly 0.0, so it stays bitwise identical at any chain count -- which
    // is what makes the speed-up measurable against an unchanged answer.
    if (rt_split) {
      constexpr int NC = RT_NB;
      // grey runs one chain: one band, one column sweep, whatever the angular quadrature
      const int nblk = rt_grey ? 1 : (nchain_rt + NC - 1)/NC;
      if (rt_tau_ptr == nullptr) {
        const int nmb = pmbp->nmb_thispack;
        // Deliberately leaked, like the k-table: a namespace-scope View would outlive
        // Kokkos::finalize(). rt_Fb is the only large one, nmb*nblk*n3*n2*n1 Reals.
        rt_tau_ptr = new DvceArray4D<Real>("rt_tau", nmb, n3, n2, n1);
        rt_B_ptr   = new DvceArray4D<Real>("rt_B",   nmb, n3, n2, n1);
        rt_Qv_ptr  = new DvceArray4D<Real>("rt_Qv",  nmb, n3, n2, n1);
        rt_cf_ptr  = new DvceArray4D<Real>("rt_cf",  nmb, n3, n2, 4);
        // (m, slot, i, k, j). Threads in a wave vary in j -- see the flattening in
        // par_for -- so j has to be the FASTEST array index. With j second-slowest, as
        // (m,slot,k,j,i) had it, adjacent lanes were 544 bytes apart and each one pulled
        // its own cache line: 8 useful bytes out of every 64 fetched.
        rt_Fb_ptr  = new DvceArray5D<Real>("rt_Fb",  nmb, nblk, n1, n3, n2);
        rt_Em_ptr  = new DvceArray5D<Real>("rt_Em",  nmb, nblk, n1, n3, n2);
        rt_Src_ptr = new DvceArray5D<Real>("rt_Src", nmb, nblk, n1, n3, n2);
        if (rt_ck || rt_grey) {
          const int nb_a = rt_ck ? CK_NB : 1;
          rt_kc_ptr = new DvceArray5D<Real>("rt_kc", nmb, nb_a, n1, n3, n2);
          rt_Bb_ptr = new DvceArray5D<Real>("rt_Bb", nmb, nb_a, n1, n3, n2);
          rt_T_ptr  = new DvceArray4D<Real>("rt_T",  nmb, n3, n2, n1);
          rt_pb_ptr = new DvceArray4D<Real>("rt_pb", nmb, n3, n2, n1);
          rt_xT_ptr = new DvceArray4D<Real>("rt_xT", nmb, n3, n2, n1);
          rt_xP_ptr = new DvceArray4D<Real>("rt_xP", nmb, n3, n2, n1);
          rt_icut_ptr = new DvceArray3D<int>("rt_icut", nmb, n3, n2);
          rt_Qb_ptr = new DvceArray5D<Real>("rt_Qb", nmb, nblk, n1, n3, n2);
        }
        if (rt_diag) {
          rt_diag_ptr = new DvceArray5D<Real>("rt_diag", nmb, 7, n3, n2, n1);
        }
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: RT split path ON, " << nblk
                    << " chain block(s) of " << NC << " -> "
                    << nblk << "x the thread count of the monolithic kernel"
                    << std::endl;
        }
      }
      auto tau_g = *rt_tau_ptr;
      auto B_g   = *rt_B_ptr;
      auto Qv_g  = *rt_Qv_ptr;
      auto cf_g  = *rt_cf_ptr;
      auto Fb_g  = *rt_Fb_ptr;
      auto Em_g  = *rt_Em_ptr;
      auto Src_g = *rt_Src_ptr;
      const bool ck_on = rt_ck;
      // the grey path shares the correlated-k scaffolding: the per-cell (T, p) and
      // opacity precompute, the cut, the tau blend and the semi-implicit application.
      // band_on says "one of the two band solvers is running", ck_on says which.
      const bool grey_on = rt_grey;
      const bool band_on = ck_on || grey_on;
      // problem/rt_top_clamp: top slot reads cell ie instead of the hydro ghost
      const bool topclamp = rt_top_clamp;
      // the optical-depth blend with the conduction module's radiative diffusion: its
      // x1-face weight w (rad_w) says how much of each face's longwave flux the
      // two-stream still owns (1 - w); the column's RT bottom is the first face with
      // w = 1, and no internal flux is injected there (the diffusion carries it)
      Conduction *pcond_rt = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->pcond
                                                            : pm->pmb_pack->phydro->pcond;
      const bool taublend = (band_on && pcond_rt != nullptr &&
                             pcond_rt->rad_tau_mode);
      auto w_g = taublend ? pcond_rt->rad_w : DvceArray4D<Real>("rt_w_dummy",1,1,1,1);
      auto tauf_g = taublend ? pcond_rt->rad_tauf
                             : DvceArray4D<Real>("rt_tau_dummy",1,1,1,1);
      if (taublend) int_at_cut = false;
      // rt_top_re is silently a no-op when the ghost is radiatively inert: see the note
      // on rt_top_re.  Say so once, from rank 0; the behaviour is unchanged.
      {
        static bool topre_warned = false;
        // NOTE the DENSITY gate is deliberately not part of this test: with
        // rad_gate_rho the top ghost's opacity is G(rho_top)*kappa, which is zero only
        // where the gate is actually closed, so rt_top_re is a no-op there and live
        // everywhere else.  Warning on it would be wrong more often than right.
        if (top_re && !topre_warned && pcond_rt != nullptr &&
            pcond_rt->rad_kappa_rmax > 0.0 && pcond_rt->rad_kappa_above == 0.0 &&
            global_variable::my_rank == 0) {
          topre_warned = true;
          std::cout << "### WARNING in two_stream_rt: problem/rt_top_re is ON but "
                    << "rad_kappa_above = 0, so the top ghost has zero optical depth and "
                    << "the back-radiation it hands down is exactly zero. The switch has "
                    << "no effect in this configuration." << std::endl;
        }
      }
      auto kc_g   = (band_on) ? *rt_kc_ptr : Fb_g;
      auto Bb_g   = (band_on) ? *rt_Bb_ptr : Fb_g;
      auto T_g    = (band_on) ? *rt_T_ptr  : tau_g;
      auto pb_g   = (band_on) ? *rt_pb_ptr : tau_g;
      auto xT_g   = (band_on) ? *rt_xT_ptr : tau_g;
      auto xP_g   = (band_on) ? *rt_xP_ptr : tau_g;
      auto icut_g = (band_on) ? *rt_icut_ptr : DvceArray3D<int>("dummy",1,1,1);
      auto Qb_g   = (band_on) ? *rt_Qb_ptr : Fb_g;
      // see rt_cell_report: the per-face streams the report needs, q = 0 only
      const bool report_on = rt_cell_report && band_on;
      if (report_on && rt_idn_ptr == nullptr) {
        const int nmb_r = pmbp->nmb_thispack;
        rt_idn_ptr = new DvceArray4D<Real>("rt_idn", nmb_r, n3, n2, n1);
        rt_iup_ptr = new DvceArray4D<Real>("rt_iup", nmb_r, n3, n2, n1);
      }
      auto idn_g = report_on ? *rt_idn_ptr : DvceArray4D<Real>("rt_idn_d", 1, 1, 1, 1);
      auto iup_g = report_on ? *rt_iup_ptr : DvceArray4D<Real>("rt_iup_d", 1, 1, 1, 1);
      // the grey opacity: the conduction module's own table if it has one, else the
      // Freedman fit, which is what the old grey path used unconditionally
      const bool grey_ktab = grey_on && pcond_rt != nullptr &&
                             pcond_rt->rad_kappa_tab && pcond_rt->rad_kr_nT > 0;
      const bool grey_krho = grey_ktab && pcond_rt->rad_kappa_rho;
      auto grey_kt  = grey_ktab ? pcond_rt->rad_kr_tab : DvceArray2D<Real>("d",1,1);
      auto grey_klT = grey_ktab ? pcond_rt->rad_kr_lT : DvceArray1D<Real>("d",1);
      auto grey_klP = grey_ktab ? pcond_rt->rad_kr_lP : DvceArray1D<Real>("d",1);
      const int grey_nT = grey_ktab ? pcond_rt->rad_kr_nT : 0;
      const int grey_nP = grey_ktab ? pcond_rt->rad_kr_nP : 0;
      const Real grey_kfac = (pcond_rt != nullptr) ? pcond_rt->rad_kappa_fac : 1.0;
      // the radiatively inert region the conduction module defines (rad_kappa_rmax): the
      // two-stream has to go quiet over exactly the same cells, or the corona the
      // diffusion refuses to touch would still be cooled by the band solver.  Zero
      // opacity is exact here, not a limit: dtau = 0 gives e0 = -expm1(0) = 0, so alp,
      // bet and gm all vanish, Src and Em pick up nothing, and rt_apply's Newton branch
      // is skipped because Em <= 0.  With tau flat through the corona the blend weight w
      // is 0 there too, so the w*F handover term adds nothing either.
      const Real grey_krmax = (pcond_rt != nullptr) ? pcond_rt->rad_kappa_rmax : 0.0;
      const Real grey_kabove = (pcond_rt != nullptr) ? pcond_rt->rad_kappa_above : 0.0;
      // ...and its DENSITY form (rad_gate_rho), which is the one to use when the
      // artificial medium and the star exchange gas: kappa_eff = G kappa + (1-G) kabove
      // with G = RadGate(rho).  Applied to kc_g, which is the ONE array every sweep
      // below reads -- the top-ghost dtau, the BFace emissivity weighting and the
      // up/down sweeps all go through it, so gating it here gates the whole solver.
      const Real gate_rho = (pcond_rt != nullptr) ? pcond_rt->rad_gate_rho : 0.0;
      const Real gate_dex = (pcond_rt != nullptr) ? pcond_rt->rad_gate_dex : 0.5;
      auto ckswf  = (ck_on) ? *ck_swf_ptr : DvceArray1D<Real>("d",1);
      auto cklk = (ck_on) ? *ck_lk_ptr : DvceArray4D<Real>("d",1,1,1,1);
      auto cklT = (ck_on) ? *ck_lT_ptr : DvceArray1D<Real>("d",1);
      auto cklP = (ck_on) ? *ck_lP_ptr : DvceArray1D<Real>("d",1);
      auto ckgw = (ck_on) ? *ck_gw_ptr : DvceArray1D<Real>("d",1);
      auto ckwl = (ck_on) ? *ck_wl_ptr : DvceArray1D<Real>("d",1);
      auto ckpf = (ck_on) ? *ck_pf_ptr : DvceArray2D<Real>("d",1,1);
      auto cece = (ck_on) ? *ce_ptr : DvceArray3D<Real>("d",1,1,1);
      auto celT = (ck_on) ? *ce_lT_ptr : DvceArray1D<Real>("d",1);
      auto celP = (ck_on) ? *ce_lP_ptr : DvceArray1D<Real>("d",1);
      auto cian = (ck_on) ? *cia_nT_ptr : DvceArray1D<int>("d",1);
      auto ciaT = (ck_on) ? *cia_T_ptr : DvceArray2D<Real>("d",1,1);
      auto ciak = (ck_on) ? *cia_k_ptr : DvceArray3D<Real>("d",1,1,1);
      auto rayx = (ck_on) ? *ray_x_ptr : DvceArray2D<Real>("d",1,1);
      const int ckNT = ck_nT;
      const int ckNP = ck_nP;
      const int ceNT = ce_nT;
      const int ceNP = ce_nP;
      const Real pfl0 = ck_pf_lTmin;
      const Real pfid = ck_pf_idlT;
      const Real pcut = rt_ck_pcut;
      const int ck_nq_ = ck_nq;

      if (!band_on && n1 > RT_NNC) {
        // the correlated-k path dispatches its column size at run time; the grey split
        // path below still uses the fixed RT_NNC, so it has to be checked
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_split with grey "
                  << "RT needs RT_NNC >= n1, but RT_NNC = " << RT_NNC
                  << " and n1 = " << n1
                  << ". Rebuild with -DRT_NNC=" << n1 << ", or use the monolithic path "
                  << "(problem/rt_split=false), which sizes its own arrays." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // ---- A (correlated-k): the per-cell work has NO radial dependency, so it does
      // not belong in a per-column kernel. Splitting it out takes it from 8192 threads at
      // 0.52 waves/CU and 2.8 % VALUBusy -- the same starvation the monolithic kernel
      // had -- to one thread per cell. Only mu0 and the cut index are per column, and
      // with correlated-k on, the grey optical depth sweep, the grey Planck function and
      // the three-band Q_v that the old rt_pre computed are all dead: nothing reads them.
      if (band_on) {
        par_for("rt_pre_geom", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          const Real x2v = x2v_(m,j);
          const Real x3v = x3v_(m,k);
          Real lam, phi, theta;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -theta+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            lam = x3v*iap;
            theta = -lam+M_PI/2.0;
            phi = x2v*iap;
          }
          Real mu0 = sin(theta)*cos(phi);
          if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
          cf_g(m,k,j,3) = mu0;
        });
        par_for("rt_pre_tp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          // THE TOP SLOT READS THE TOP ACTIVE CELL, NOT THE HYDRO GHOST.  i = ie+1 is the
          // unresolved column above the domain, and physically it IS the continuation of
          // cell ie -- taking its state from ie rather than from w0's ghost decouples the
          // whole band solver from whatever the boundary condition put there.  R9's death
          // (t = 1.99e5, hydro dt 37 s -> 7e-15 in one step, 1.64e6 non-finite cells over
          // i ~ 200..323 on every rank in 26 cycles) is what a single poisoned ghost does
          // once the down-sweep starts from it; nothing else spreads that fast.
          const int ii = (topclamp && i > ie) ? ie : i;
          Real pp, TT;
          const Real dd = rhoN(m,k,j,ii);
          PresTempFromEint(eos,gm1,Rgas,dd,eiN(m,k,j,ii),
                           TGuess(wtemp_, m, k, j, ii),pp,TT);
          if (RTBadState(pp, TT)) {          // see RTBadState: one NaN kills the column
            Kokkos::atomic_fetch_add(&stcl_g(0), 1);
            TT = 0.0;                        // the marker the opacity kernels read
            pp = 0.0;
          }
          T_g(m,k,j,i) = TT;
          pb_g(m,k,j,i) = pp*1.0e-6;
        });
        par_for("rt_pre_cut", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          // i = is is the bottom, so pressure falls as i rises: the cut is the deepest
          // cell still shallower than pcut, i.e. the first one scanning up
          int icut = ie+1;
          if (taublend) {
            // the RT reaches down to the first face it still owns a share of
            for (int i=is; i<ie+2; ++i) {
              if (w_g(m,k,j,i) < 1.0) { icut = i; break; }
            }
          } else {
            for (int i=is; i<ie+2; ++i) {
              if (pb_g(m,k,j,i) < pcut) { icut = i; break; }
            }
          }
          icut_g(m,k,j) = icut;
        });
        if (grey_on) {
          par_for("rt_pre_opac_grey", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            if (i < icut_g(m,k,j)) return;        // deeper than the cut: never read
            const int ii = (topclamp && i > ie) ? ie : i;  // top slot: see rt_pre_tp
            // T = 0 is rt_pre_tp's marker for a state the EOS could not form (see
            // RTBadState).  Such a cell takes NO part in this RT call: zero opacity and
            // zero emission make the layer transparent (e0 = 0, so alp = bet = 0, the
            // streams pass through untouched, Src and Em pick up nothing and rt_apply's
            // Newton branch is skipped) -- the same "exact" inert state the corona above
            // rad_kappa_rmax already has.  Giving it a FLOOR state instead would hand a
            // 0.1 K blackbody to a cell with 8000 K neighbours, which collapses dt.
            if (!(T_g(m,k,j,i) > 0.0)) {
              kc_g(m,0,i,k,j) = 0.0;
              Bb_g(m,0,i,k,j) = 0.0;
              return;
            }
            const Real TT = T_g(m,k,j,i);
            const Real pcgs = pb_g(m,k,j,i)*1.0e6;
            const Real rho = rhoN(m,k,j,ii);
            Real kr = (grey_krmax > 0.0 && x1v_(m,i) > grey_krmax)
                ? grey_kabove
                : (grey_ktab
                   ? RosselandTable(grey_kt, grey_klT, grey_klP, grey_nT, grey_nP, TT,
                                    grey_krho ? rho : pcgs)
                   : RosselandFreedman2014(TT, pcgs, met));
            if (gate_rho > 0.0) {
              const Real g = RadGate(rho, gate_rho, gate_dex);
              kr = g*kr + (1.0 - g)*grey_kabove;
            }
            kc_g(m,0,i,k,j) = grey_kfac*kr;
            Bb_g(m,0,i,k,j) = boltz_sigma/M_PI*SQR(SQR(TT));
          });
        } else {
        par_for("rt_pre_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          if (i < icut_g(m,k,j)) return;          // deeper than the cut: never read
          const int ii = (topclamp && i > ie) ? ie : i;  // top slot: see rt_pre_tp
          if (!(T_g(m,k,j,i) > 0.0)) {            // unusable state: inert, see the grey
            xT_g(m,k,j,i) = 0.0;                  // kernel's note
            xP_g(m,k,j,i) = 0.0;
            for (int b=0; b<CK_NB; ++b) {
              kc_g(m,b,i,k,j) = 0.0;
              Bb_g(m,b,i,k,j) = 0.0;
            }
            return;
          }
          const Real TT = T_g(m,k,j,i);
          const Real pbar = pb_g(m,k,j,i);
          int iT, iP;
          Real fT, fP;
          ck_tp_index(cklT, ckNT, log10(TT), iT, fT);
          ck_tp_index(cklP, ckNP, log10(pbar), iP, fP);
          xT_g(m,k,j,i) = static_cast<Real>(iT) + fT;
          xP_g(m,k,j,i) = static_cast<Real>(iP) + fP;
          Real kcb[CK_NB];
          ck_continuum(cece, celT, celP, ceNT, ceNP, cian, ciaT, ciak, rayx, ckwl,
                       TT, pbar, rhoN(m,k,j,ii), kcb);
          // the density gate.  The correlated-k path never carried the rad_kappa_rmax
          // radius test -- the inert corona is a grey-path feature -- but the gate is a
          // property of the GAS, so it must reach every band here as well or a ck run
          // would keep heating the medium the grey run refuses to touch.
          const Real gk = (gate_rho > 0.0)
              ? RadGate(rhoN(m,k,j,ii), gate_rho, gate_dex) : 1.0;
          const Real sigT4_pi = boltz_sigma/M_PI*SQR(SQR(TT));
          for (int b=0; b<CK_NB; ++b) {
            kc_g(m,b,i,k,j) = (gate_rho > 0.0)
                ? (gk*kcb[b] + (1.0 - gk)*grey_kabove) : kcb[b];
            Bb_g(m,b,i,k,j) = sigT4_pi*ck_planck_frac(ckpf, pfl0, pfid, TT, b);
          }
        });
        }
      } else {
      // ---- A: chain-independent per-column precompute -----------------------------
      par_for("rt_pre", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        // Subviews so the body below is the monolithic kernel's code unchanged: these
        // stand in for the private tau_down_r_f[NN], B[NN], Q_v[NN] arrays.
        auto tau_down_r_f = Kokkos::subview(tau_g, m, k, j, Kokkos::ALL);
        auto B            = Kokkos::subview(B_g,   m, k, j, Kokkos::ALL);
        auto Q_v          = Kokkos::subview(Qv_g,  m, k, j, Kokkos::ALL);
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);

        Real rtop = x1v_(m,ie+1);
        Real rbot = x1v_(m,is);

        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);

        Real mus = (mu0 > 0.0) ? mu0 : 0.0;
        Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
        Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
        get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);

        // 3 V Bands
        // top
        // the top slot is the continuation of the top ACTIVE cell, not the hydro ghost
        // (see rt_pre_tp): the band solver must not be able to read a poisoned ghost
        const int itop = topclamp ? ie : (ie+1);
        Real rho = rhoN(m,k,j,itop);
        Real p = PresFromEint(eos,gm1,rho,eiN(m,k,j,itop));
        Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,itop),p);
        bool badtop = RTBadState(p, T);      // see RTBadState: one NaN kills the column
        if (badtop) {
          Kokkos::atomic_fetch_add(&stcl_g(0), 1);
          T = 0.0;
          p = 0.0;
        }
        B[ie+1] = badtop ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr = 0.0;
        if (!badtop) get_kapr(T, p, met, kapr);
        // tau = kappa p / g for the UNRESOLVED column above the domain. g must be the
        // value at the top, not the surface value: at r/ap = 1.5 they differ by 2.3x.
        // The tidal term belongs here too: it is the EFFECTIVE gravity that sets how
        // much mass the unresolved column above the domain holds. See EffGravAt.
        Real tau_r_f = RTTopDtau(kapr, p,
                                 EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0,
                                           tide));
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = (kapr*rho != 0.0) ? tau_r_f/(kapr*rho) : 0.0;
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
        // beam transmission at the face ABOVE the cell being filled, carried down the
        // sweep so that differencing the flux across a cell costs no extra exp
        Real trp1 = trans1;
        Real trp2 = trans2;
        Real trp3 = trans3;
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = rhoN(m,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,eiN(m,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          bool badcell = RTBadState(p, T);   // see RTBadState: one NaN kills the column
          if (badcell) {
            Kokkos::atomic_fetch_add(&stcl_g(0), 1);
            T = 0.0;
            p = 0.0;
          }
          B[i] = badcell ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr = 0.0;
          if (!badcell) get_kapr(T, p, met, kapr);
          Real dr = dx1(m,k,j,i);
          tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
          Real r = x1f_(m,i);
////          Real delta = (drtop+(rtop-r))/r;
//          Real delta = dr/r;
//          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
//          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
          Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
          Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
//          Real trans1 = exp(-gamv1*tausl);
//          Real trans2 = exp(-gamv2*tausl);
//          Real trans3 = exp(-gamv3*tausl);
//          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
          Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
          // Deposit the flux DIFFERENCE across the cell. The old form,
          // kappa rho F exp(-tau) with tau at the lower face, is right only for a thin
          // layer: it returns u e^-u / (1 - e^-u) of what the cell actually absorbs, with
          // u = dtau/mu, which is 0.95 at u = 0.1 but 0.58 at u = 1. On this grid, about
          // 0.46 scale heights per cell, that put dtau/mu near one wherever it mattered
          // and lost about 24 % of the incident stellar flux. Found by comparing the
          // correlated-k version of the same expression against Exo-FMS on an identical
          // column. Written this way the column integral telescopes to
          // mu F (1 - e^-tau_total) exactly, and it still reduces to the old expression
          // as dtau -> 0.
          Real Qv = (1.0-albedo)*Fstar*(1.0/3.0)
                  * ((trp1-trans1)+(trp2-trans2)+(trp3-trans3))/(fac*dr);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
          trp1 = trans1;
          trp2 = trans2;
          trp3 = trans3;
        }
        cf_g(m,k,j,0) = gamir1;
        cf_g(m,k,j,1) = gamir2;
        cf_g(m,k,j,2) = beta;
        cf_g(m,k,j,3) = mu0;
      });
      }

      // ---- B (correlated-k): one thread per (column, chain block) ------------------
      // Chain c decodes as c = ((band*CK_NG) + g)*2 + quad, so a block of RT_NB = 4
      // consecutive chains is two g-points x two angular quadrature points of ONE band,
      // which is what lets the block share a band's continuum and Planck function.
      //
      // The longwave needs no cumulative optical depth: dtau is kappa*rho*dr, purely
      // local, so going from a grey tau scaled by gamma to a per-chain kappa is a lookup
      // and nothing structural. Only the sweep's lower limit changes, from is to icut.
      if (grey_on) {
        // ---- B (grey): one thread per column.  One band, one opacity, one sweep.
        //
        // Structurally this is the correlated-k kernel with the chain loop collapsed to
        // the angular quadrature: the same layer-integrated coefficients (alpha, beta on
        // the linear-in-tau source), the same column-above-the-domain top boundary, the
        // same thermalised intensity at the cut, and the same flux accumulation.  What it
        // does NOT carry is the band structure, and that is the point: with one opacity
        // the cell's emission and its absorption of the field scale together, so the
        // local balance has a stable fixed point.
        //
        // Em is the cell's own emission per unit volume and it is written EXACTLY here,
        // 4 sigma kappa rho T^4, rather than through the hemispheric-mean quadrature the
        // band kernel uses (17 % low).  It only sets the implicit relaxation rate, but
        // the rate is what decides whether a thin cell can be integrated at the
        // hydrodynamic timestep, so it is worth having right.
        auto launch_grey_chain = [&](auto nn_tag) {
          constexpr int NN = decltype(nn_tag)::value;
          par_for("rt_chain_grey", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            for (int i=is; i<ie+2; ++i) {
              Fb_g(m,0,i,k,j) = 0.0;
              Qb_g(m,0,i,k,j) = 0.0;
              Em_g(m,0,i,k,j) = 0.0;
              Src_g(m,0,i,k,j) = 0.0;
            }
            const int icut = icut_g(m,k,j);
            if (icut > ie) return;                  // whole column deeper than the cut
            // ck_nquad = 1 is the hemispheric mean (mu = 1/1.66), 2 the two-point
            // Gauss-Legendre quadrature the band solver offers on the same switch
            const int nq = (ck_nq_ > 1) ? 2 : 1;
            Real muq[2], wfq[2];
            if (nq == 1) {
              muq[0] = 1.0/CK_DIFFUSIVITY;
              wfq[0] = M_PI;                        // F = pi I
            } else {
              for (int q=0; q<2; ++q) {
                muq[q] = mug[q];
                wfq[q] = 2.0*M_PI*wg[q]*mug[q];
              }
            }
            Real I_down[2][NN];
            // Top: the unresolved hydrostatic column above the domain, p/g of it, at the
            // top cell's opacity -- the same construction the band solver uses.
            {
              const Real kap = kc_g(m,0,ie+1,k,j);
              const Real mu0 = cf_g(m,k,j,3);
              const Real dtau = RTTopDtau(kap, pb_g(m,k,j,ie+1)*1.0e6,
                                          EffGravAt(grav, ap, x1v_(m,ie+1), grav_pmass,
                                                    omega, mu0, tide));
              // The source of that column: the ghost's own Planck function, or -- with
              // rt_top_re -- the radiative-equilibrium value I_up/2, which is what a slab
              // with nothing but space above it must emit downward.  See rt_top_re.
              //
              // I_up at the top face does not depend on I_down anywhere: the up-sweep
              // reads only the intensity it starts with at the cut and the Planck
              // function.  So the equilibrium source can be had EXACTLY, with a scalar
              // probe sweep here and no lag and no stored state -- one extra pass over
              // the column, run only when this boundary is selected.
              Real bsrc[2];
              for (int q=0; q<nq; ++q) {
                bsrc[q] = Bb_g(m,0,ie+1,k,j);
              }
              if (top_re) {
                for (int q=0; q<nq; ++q) {
                  Real ip = Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0);
                  for (int i=icut+1; i<ie+2; ++i) {
                    const Real krb = kc_g(m,0,i-1,k,j)*rhoN(m,k,j,i-1);
                    const Real x = krb*dx1(m,k,j,i-1)/muq[q];
                    const Real e0 = -expm1(-x);
                    const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                    const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                                  : (x/2.0 - SQR(x)/3.0);
                    // emissivity-weighted far endpoint: see BFace
                    const int iir = (i > ie) ? ie : i;
                    const Real bfar = BFace(krb, kc_g(m,0,i,k,j)*rhoN(m,k,j,iir),
                                            Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
                    ip = (1.0-e0)*ip + bet*bfar + gm*Bb_g(m,0,i-1,k,j);
                  }
                  bsrc[q] = 0.5*ip;
                }
              }
              for (int q=0; q<nq; ++q) {
                I_down[q][ie+1] = (1.0 - exp(-dtau/muq[q]))*bsrc[q];
              }
              if (report_on) idn_g(m,k,j,ie+1) = I_down[0][ie+1];
            }
            // down-sweep
            for (int i=ie; i>icut-1; --i) {
              const Real krb_d = kc_g(m,0,i,k,j)*rhoN(m,k,j,i);
              const Real dtau_i = krb_d*dx1(m,k,j,i);
              // the far endpoint of the layer's source function, weighted by emitting
              // matter: above rad_kappa_rmax the neighbour has kappa = 0 and a Planck
              // function 1e9x this cell's, which used to drain it to the floor in one
              // step (see BFace).  Identical to Bb_g(i+1) at equal opacity.
              const int iip = (i+1 > ie) ? ie : i+1;
              const Real bfar_d = BFace(krb_d, kc_g(m,0,i+1,k,j)*rhoN(m,k,j,iip),
                                        Bb_g(m,0,i,k,j), Bb_g(m,0,i+1,k,j), bface_on);
              for (int q=0; q<nq; ++q) {
                const Real x = dtau_i/muq[q];
                const Real e0 = -expm1(-x);
                const Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
                const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                // direct source: what this stream leaves in cell i, absorbed minus
                // emitted
                Src_g(m,0,i,k,j) += wfq[q]/dx1(m,k,j,i)
                                  *(e0*I_down[q][i+1]
                                    - (alp*bfar_d + bet*Bb_g(m,0,i,k,j)));
                I_down[q][i] = (1.0-e0)*I_down[q][i+1]
                             + alp*bfar_d + bet*Bb_g(m,0,i,k,j);
              }
              if (report_on) idn_g(m,k,j,i) = I_down[0][i];
            }
            // Bottom of the RT domain: thermalised, plus the internal flux if the layers
            // below are not carrying it themselves (see rt_int_at_cut).
            Real I_up[2];
            for (int q=0; q<nq; ++q) {
              I_up[q] = Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0);
              Fb_g(m,0,icut,k,j) += wfq[q]*(I_up[q] - I_down[q][icut]);
            }
            if (report_on) iup_g(m,k,j,icut) = I_up[0];
            // up-sweep
            for (int i=icut+1; i<ie+2; ++i) {
              const Real kap = kc_g(m,0,i-1,k,j);
              const Real rho = rhoN(m,k,j,i-1);
              const Real dtau_i = kap*rho*dx1(m,k,j,i-1);
              // same emissivity weighting for the upward stream and for Em: see BFace
              const int iiu = (i > ie) ? ie : i;
              const Real bfar_u = BFace(kap*rho, kc_g(m,0,i,k,j)*rhoN(m,k,j,iiu),
                                        Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
              for (int q=0; q<nq; ++q) {
                const Real x = dtau_i/muq[q];
                const Real e0 = -expm1(-x);
                const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
                const Real Iup_in = I_up[q];
                // the same for the upward stream through layer i-1
                Src_g(m,0,i-1,k,j) += wfq[q]/dx1(m,k,j,i-1)
                                    *(e0*Iup_in
                                      - (bet*bfar_u
                                         + gm*Bb_g(m,0,i-1,k,j)));
                I_up[q] = (1.0-e0)*Iup_in
                        + bet*bfar_u + gm*Bb_g(m,0,i-1,k,j);
                Fb_g(m,0,i,k,j) += wfq[q]*(I_up[q] - I_down[q][i]);
              }
              if (report_on) iup_g(m,k,j,i) = I_up[0];
              Em_g(m,0,i-1,k,j) = 4.0*M_PI*kap*rho
                                * 0.5*(bfar_u + Bb_g(m,0,i-1,k,j));
            }
          });
        };
        if (n1 <= 72) {
          launch_grey_chain(std::integral_constant<int, 72>{});
        } else if (n1 <= 136) {
          launch_grey_chain(std::integral_constant<int, 136>{});
        } else if (n1 <= 264) {
          launch_grey_chain(std::integral_constant<int, 264>{});
        } else if (n1 <= 520) {
          launch_grey_chain(std::integral_constant<int, 520>{});
        } else {
          std::cout << "### FATAL ERROR in two_stream_rt: n1 = " << n1
                    << " exceeds the largest grey radial tier (520). Add a tier to the "
                    << "dispatch in picket_fence_two_stream_RT." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else if (ck_on) {
        // The private intensity column has to be sized at COMPILE time, but the radial
        // extent is only known at run time, and an oversized one is not free: at n1 = 68
        // the chain kernel costs 454 ms with a 72-deep column, 503 at 136 and 540 at 272,
        // because the scratch footprint grows with it. So instead of one generous ceiling
        // that taxes every run, the kernel is instantiated at a few sizes and the
        // smallest
        // one that fits is dispatched at run time.
        auto launch_ck_chain = [&](auto nn_tag) {
          constexpr int NN = decltype(nn_tag)::value;
          par_for("rt_chain_ck", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
            constexpr int NC = RT_NB;
            for (int i=is; i<ie+2; ++i) {
              Fb_g(m,blk,i,k,j) = 0.0;
              Qb_g(m,blk,i,k,j) = 0.0;
              Em_g(m,blk,i,k,j) = 0.0;
              Src_g(m,blk,i,k,j) = 0.0;
            }
            const int icut = icut_g(m,k,j);
            if (icut > ie) return;                  // whole column deeper than the cut
            // Shortwave. This is the one part of the scheme that genuinely restructures:
            // the longwave only ever needs a LAYER optical depth, which is local, but the
            // direct stellar beam needs the CUMULATIVE depth from the top, so each chain
            // carries its own downward recurrence and it cannot live in the per-column
            // rt_pre. It rides along in the longwave down-sweep because the two share the
            // same kappa lookup at every cell -- doing it in a second kernel would pay
            // for
            // that lookup twice.
            const Real mu0 = cf_g(m,k,j,3);
            const Real facsw = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
            const bool lit = (mu0 > 0.0);
            Real tausw[NC];
            Real transw[NC];      // beam transmission at the face above the current cell

            int bandc[NC], gc[NC];
            Real muc[NC], wfc[NC], wgc[NC];
            for (int cc=0; cc<NC; ++cc) {
              const int c = blk*NC + cc;
              if (ck_nq_ == 1) {
                gc[cc] = c % CK_NG;
                bandc[cc] = c/CK_NG;
                muc[cc] = 1.0/CK_DIFFUSIVITY;
                wfc[cc] = M_PI*ckgw(gc[cc]);          // F = pi I
              } else {
                const int nq = c % 2;
                gc[cc] = (c/2) % CK_NG;
                bandc[cc] = c/(2*CK_NG);
                muc[cc] = mug[nq];
                wfc[cc] = 2.0*M_PI*wg[nq]*mug[nq]*ckgw(gc[cc]);
              }
              // the shortwave weights by the g-point alone: it is a direct beam, not an
              // angular quadrature, and with nquad = 2 each angular point would otherwise
              // double-count the incident flux
              wgc[cc] = ckgw(gc[cc])/static_cast<Real>(ck_nq_);
            }
            // A WARNING ABOUT MEASURING THIS KERNEL. Seven optimisations were measured
            // against it while its per-cell arrays were laid out (m,slot,k,j,i), which
            // put
            // adjacent lanes 544 bytes apart. All seven failed, and several of those
            // verdicts
            // were artefacts of that: with the wave starved on scattered loads, nothing
            // done
            // to the arithmetic could show up. Re-measured on the (m,slot,i,k,j) layout:
            //
            //                                     starved layout    coalesced layout
            //   FP32 recurrence                        +2.8 %            +1.42x
            //   RT_NB = 2 instead of 4                  -16 %      faster on rt_chain,
            //                                                      slower on the total
            //   dropping this private column            -22 %            neutral
            //   RT_NB = 8                              slower            slower
            //
            // So: do not trust a null result on this kernel without checking that memory
            // is
            // not the thing in the way. Still genuinely useless, both layouts: the
            // k-table
            // layout (0.5 %), blocking the lookup by band (0.4 %), and precomputing kappa
            // per (cell, chain) (-1 %, because those four table loads are cache hits).
            //
            // Storing the intensity column costs 2304 bytes of scratch per thread, and
            // accumulating each sweep's flux contribution separately instead would remove
            // it entirely. On the starved layout that cost 22 % (1585 -> 1927 ms),
            // because
            // the extra flux traffic it adds was uncoalesced. On the fixed layout it is
            // neutral (454 vs 457 ms), which confirms the diagnosis. Neutral is not a
            // reason to change it, so the column stays.
            RtF I_down[NC][NN];

            // Top: the column above the domain, using the top cell's opacity over the
            // hydrostatic column p/g -- the same construction the grey scheme uses.
            {
              const Real Ttop = T_g(m,k,j,ie+1);
              const Real ptop = pb_g(m,k,j,ie+1);
              int iT, iP;
              Real fT, fP;
              const Real xTv = xT_g(m,k,j,ie+1);
              const Real xPv = xP_g(m,k,j,ie+1);
              iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
              iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
              for (int cc=0; cc<NC; ++cc) {
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, bandc[cc], gc[cc])
                               + kc_g(m,bandc[cc],ie+1,k,j);
                const Real dtau = RTTopDtau(kap, ptop*1.0e6,
                                            EffGravAt(grav, ap, x1v_(m,ie+1),
                                                      grav_pmass, omega, mu0, tide));
                const RtF trans = RT_EXP(-static_cast<RtF>(dtau/muc[cc]));
                I_down[cc][ie+1] = (static_cast<RtF>(1.0)-trans)
                                 * static_cast<RtF>(Bb_g(m,bandc[cc],ie+1,k,j));
                tausw[cc] = dtau;               // beam already crossed the column above
                transw[cc] = RT_EXP(-static_cast<RtF>(dtau*facsw));
              }
            }
            // down-sweep
            for (int i=ie; i>icut-1; --i) {
              const Real rho = rhoN(m,k,j,i);
              const Real drho = rho*dx1(m,k,j,i);
              int iT, iP;
              Real fT, fP;
              const Real xTv = xT_g(m,k,j,i);
              const Real xPv = xP_g(m,k,j,i);
              iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
              iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
              // the far endpoint of this layer's source function, weighted by emitting
              // matter (see BFace).  The neighbour's kappa costs a second table lookup,
              // so it is paid only where the two Planck functions differ by more than 4x
              // -- a smooth column takes the old expression and is bit-identical.
              const int ir1 = (i+1 > ie) ? ie : i+1;
              const Real rhf = rhoN(m,k,j,ir1);
              const Real xTf = xT_g(m,k,j,i+1);
              const Real xPf = xP_g(m,k,j,i+1);
              const int iTf = static_cast<int>(xTf);
              const int iPf = static_cast<int>(xPf);
              const Real fTf = xTf - static_cast<Real>(iTf);
              const Real fPf = xPf - static_cast<Real>(iPf);
              for (int cc=0; cc<NC; ++cc) {
                const int b = bandc[cc];
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                                 + kc_g(m,b,i,k,j);
                const Real bown = Bb_g(m,b,i,k,j);
                Real bfar = Bb_g(m,b,i+1,k,j);
                if (bfar > 4.0*bown || bown > 4.0*bfar) {
                  const Real kapf = ck_kappa(cklk, iTf, fTf, iPf, fPf, b, gc[cc])
                                  + kc_g(m,b,i+1,k,j);
                  bfar = BFace(kap*rho, kapf*rhf, bown, bfar, bface_on);
                }
                const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
                const RtF e0 = -RT_EXPM1(-x);
                const RtF one = static_cast<RtF>(1.0);
                const RtF alp = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                              : (x/2 - x*x/3);
                const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                              : (x/2 - x*x/6);
                Src_g(m,blk,i,k,j) += wfc[cc]/dx1(m,k,j,i)
                    *static_cast<Real>(e0*I_down[cc][i+1]
                      - (alp*static_cast<RtF>(bfar)
                         + bet*static_cast<RtF>(bown)));
                I_down[cc][i] = (one-e0)*I_down[cc][i+1]
                              + alp*static_cast<RtF>(bfar)
                              + bet*static_cast<RtF>(bown);
                // Direct beam. Deposit the flux DIFFERENCE across the cell, not
                // kappa rho F exp(-tau) evaluated at one face. The latter is what the
                // grey
                // picket fence does, and it under-deposits badly once a layer is not
                // thin:
                // the ratio of deposited to absorbed is u e^-u/(1 - e^-u) with u =
                // dtau/mu,
                // which is 0.95 at u = 0.1 but 0.58 at u = 1 and 0.31 at u = 2. Summed
                // down
                // a column with u ~ 0.5 it loses a quarter of the incident flux --
                // measured
                // against Exo-FMS on an identical column, which is how this was found.
                // Written this way the column integral is F mu (1 - e^-tau_total) by
                // construction, and it still reduces to the old form as dtau -> 0.
                tausw[cc] += kap*drho;
                if (lit) {
                  const Real tnew = RT_EXP(-static_cast<RtF>(tausw[cc]*facsw));
                  Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(b)*wgc[cc]/facsw
                            * (transw[cc] - tnew)/dx1(m,k,j,i);
                  transw[cc] = tnew;
                }
              }
            }
            // Bottom of the CORRELATED-K DOMAIN, not of the column. At the cut the grey
            // optical depth is of order 1e4, so the layer is thermalised to e^-tau and
            // the
            // upward intensity is its own Planck function; the planet's internal flux is
            // delivered here as an extra band-weighted source. Below the cut nothing
            // radiative is applied -- that region is optically thick and convective, and
            // the flux simply passes through it.
            RtF I_up[NC];
            for (int cc=0; cc<NC; ++cc) {
              const int b = bandc[cc];
              const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                                * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
              I_up[cc] = static_cast<RtF>(Bb_g(m,b,icut,k,j) + Iint_b);
              Fb_g(m,blk,icut,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][icut]);
            }
            // up-sweep
            for (int i=icut+1; i<ie+2; ++i) {
              const Real rho = rhoN(m,k,j,i-1);
              const Real drho = rho*dx1(m,k,j,i-1);
              int iT, iP;
              Real fT, fP;
              const Real xTv = xT_g(m,k,j,i-1);
              const Real xPv = xP_g(m,k,j,i-1);
              iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
              iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
              // same emissivity-weighted far endpoint, and the same 4x guard: see BFace
              const int ir1 = (i > ie) ? ie : i;
              const Real rhf = rhoN(m,k,j,ir1);
              const Real xTf = xT_g(m,k,j,i);
              const Real xPf = xP_g(m,k,j,i);
              const int iTf = static_cast<int>(xTf);
              const int iPf = static_cast<int>(xPf);
              const Real fTf = xTf - static_cast<Real>(iTf);
              const Real fPf = xPf - static_cast<Real>(iPf);
              for (int cc=0; cc<NC; ++cc) {
                const int b = bandc[cc];
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                               + kc_g(m,b,i-1,k,j);
                const Real bown = Bb_g(m,b,i-1,k,j);
                Real bfar = Bb_g(m,b,i,k,j);
                if (bfar > 4.0*bown || bown > 4.0*bfar) {
                  const Real kapf = ck_kappa(cklk, iTf, fTf, iPf, fPf, b, gc[cc])
                                  + kc_g(m,b,i,k,j);
                  bfar = BFace(kap*rho, kapf*rhf, bown, bfar, bface_on);
                }
                const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
                const RtF e0 = -RT_EXPM1(-x);
                const RtF one = static_cast<RtF>(1.0);
                const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                              : (x/2 - x*x/6);
                const RtF gm = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                             : (x/2 - x*x/3);
                const RtF Iup_in = I_up[cc];
                Src_g(m,blk,i-1,k,j) += wfc[cc]/dx1(m,k,j,i-1)
                    *static_cast<Real>(e0*Iup_in
                      - (bet*static_cast<RtF>(bfar)
                         + gm*static_cast<RtF>(bown)));
                I_up[cc] = (one-e0)*Iup_in
                         + bet*static_cast<RtF>(bfar)
                         + gm*static_cast<RtF>(bown);
                Fb_g(m,blk,i,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][i]);
                // the cell's OWN emission, both hemispheres: in the thin limit each
                // stream adds wfc*(kap*rho*dr/mu)*B over the layer, so per unit volume
                // the two together give 2*(wfc/mu)*kap*rho*B.  Summed over bands and g
                // this is 2*CK_DIFFUSIVITY*sigma*kappa_P*rho*T^4, i.e. the exact
                // 4 sigma kappa_P rho T^4 with 1.66 in place of 2 for the hemispheric
                // mean.  17 % low, which is well inside what a rate estimate needs.
                Em_g(m,blk,i-1,k,j) += 2.0*(wfc[cc]/muc[cc])*kap*rho
                                     * 0.5*(bfar + bown);
              }
            }
          });
        };
        if (n1 <= 72) {
          launch_ck_chain(std::integral_constant<int, 72>{});
        } else if (n1 <= 136) {
          launch_ck_chain(std::integral_constant<int, 136>{});
        } else if (n1 <= 264) {
          launch_ck_chain(std::integral_constant<int, 264>{});
        } else if (n1 <= 520) {
          launch_ck_chain(std::integral_constant<int, 520>{});
        } else {
          std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: n1 = " << n1
                    << " exceeds the largest correlated-k radial tier (520). Add a tier "
                        << "to "
                    << "the dispatch in picket_fence_two_stream_RT." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
      par_for("rt_chain", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
        constexpr int NN = RT_NNC;
        auto tau_down_r_f = Kokkos::subview(tau_g, m, k, j, Kokkos::ALL);
        auto B            = Kokkos::subview(B_g,   m, k, j, Kokkos::ALL);
        auto F_ir_f       = Kokkos::subview(Fb_g,  m, blk, k, j, Kokkos::ALL);
        const Real gamir1 = cf_g(m,k,j,0);
        const Real gamir2 = cf_g(m,k,j,1);
        const Real beta   = cf_g(m,k,j,2);
        // this block owns its own flux slot, so it starts from zero exactly as the
        // serial accumulator did
        for (int i=is; i<ie+2; ++i) {
          Fb_g(m,blk,i,k,j) = 0.0;
          Em_g(m,blk,i,k,j) = 0.0;
          Src_g(m,blk,i,k,j) = 0.0;
        }
          Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC];
          for (int cc=0; cc<NC; ++cc) {
            const int c = blk*NC + cc;
            const int n = (c/2) % 2;
            const int vir = c % 2;
            muggc[cc] = mug[n];
            wggc[cc] = wg[n];
            gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
            fbc[cc] = (vir == 0) ? beta : (1.0-beta);
          }
          Real I_ir_down_c[NC][NN];
#if RT_CACHE
          Real e0c[NC][NN], alpc[NC][NN], betc[NC][NN];
#endif

          // top
          for (int cc=0; cc<NC; ++cc) {
            Real dtauir = gamirc[cc]*tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/muggc[cc]);
            I_ir_down_c[cc][ie+1] = (1.0-trans)*(fbc[cc]*B[ie+1]);
          }
          // down-sweep
          for (int i=ie; i>is-1; --i) {
            Real dtau_i = tau_down_r_f[i]-tau_down_r_f[i+1];
            // emissivity-weighted far endpoint of the layer's source function: dtau/dr
            // is kappa rho, so the two layers' optical thicknesses per unit length are
            // the weights (see BFace).  The ghost above ie has no dtau in the array, so
            // the topmost layer keeps the old endpoint exactly.
            const Real kro_g = dtau_i/dx1(m,k,j,i);
            const Real krf_g = (i < ie) ? (tau_down_r_f[i+1]-tau_down_r_f[i+2])
                                        / dx1(m,k,j,i+1) : kro_g;
            const Real bfr_g = BFace(kro_g, krf_g, B[i], B[i+1], bface_on);
            for (int cc=0; cc<NC; ++cc) {
              Real dtauir = gamirc[cc]*dtau_i;
              Real x = dtauir/muggc[cc];
              Real e0 = -expm1(-x);
              Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              // direct source: what this stream leaves in cell i, absorbed minus emitted
              Src_g(m,blk,i,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/dx1(m,k,j,i)
                                  *(e0*I_ir_down_c[cc][i+1]
                                    - fbc[cc]*(alp*bfr_g + bet*B[i]));
              I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                                 + alp*fbc[cc]*bfr_g + bet*fbc[cc]*B[i];
#if RT_CACHE
              e0c[cc][i] = e0;
              alpc[cc][i] = alp;
              betc[cc][i] = bet;
#endif
            }
          }

          // bottom
          Real I_ir_up_c[NC];
          for (int cc=0; cc<NC; ++cc) {
            I_ir_up_c[cc] = Iint + I_ir_down_c[cc][is];
            Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is];
            Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
            Fb_g(m,blk,is,k,j) += (F_ir_up_f - F_ir_down_f);
          }
          // up-sweep, accumulating the band flux as it goes
          for (int i=is+1; i<ie+2; ++i) {
#if !RT_CACHE
            Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
#endif
            // the same emissivity-weighted far endpoint for the upward stream and for
            // Em; the topmost layer keeps the old endpoint (see the down-sweep)
            const Real krou = (tau_down_r_f[i-1]-tau_down_r_f[i])/dx1(m,k,j,i-1);
            const Real krfu = (i < ie+1) ? (tau_down_r_f[i]-tau_down_r_f[i+1])
                                         / dx1(m,k,j,i) : krou;
            const Real bfru = BFace(krou, krfu, B[i-1], B[i], bface_on);
            for (int cc=0; cc<NC; ++cc) {
#if RT_CACHE
              // layer i-1, already solved on the way down
              const Real e0 = e0c[cc][i-1];
              const Real bet = betc[cc][i-1];
              const Real gm = alpc[cc][i-1];
#else
              Real dtauir = gamirc[cc]*dtau_i;
              Real x = dtauir/muggc[cc];
              Real e0 = -expm1(-x);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
#endif
              const Real Iup_in = I_ir_up_c[cc];
              Src_g(m,blk,i-1,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/dx1(m,k,j,i-1)
                                    *(e0*Iup_in - fbc[cc]*(bet*bfru + gm*B[i-1]));
              I_ir_up_c[cc] = (1.0-e0)*Iup_in
                            + bet*fbc[cc]*bfru + gm*fbc[cc]*B[i-1];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
              Fb_g(m,blk,i,k,j) += (F_ir_up_f - F_ir_down_f);
              // the cell's own emission per unit volume, both hemispheres.  dtau_i is
              // kappa rho dr for this layer, so dtau_i/dr = kappa rho and the layer
              // thickness cancels out of the volumetric rate.
              const Real dtau_ly = tau_down_r_f[i-1]-tau_down_r_f[i];
              Em_g(m,blk,i-1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_ly
                                   / dx1(m,k,j,i-1)*fbc[cc]*0.5*(bfru+B[i-1]);
            }
          }
      });
      }

      // ---- C: reduce over blocks in order, then apply ------------------------------
      int nclip = 0;
      const Real demax = rt_de_max;
      const bool semilin = rt_semi_lin;
      const bool explicit_on = rt_explicit;
      const bool newton_on = rt_newton;
      if (rt_efix_ptr == nullptr) {
        rt_efix_ptr = new DvceArray1D<int>("rt_efix", 3);
        Kokkos::deep_copy(*rt_efix_ptr, 0);
      }
      auto efix_g = *rt_efix_ptr;
      const bool resc_eq = rt_rescue_eq;
      // rt_cell_report: claimed once per RT call, so only the FIRST rescued cell prints
      DvceArray1D<int> repc_g(std::string("rt_repc"), 1);
      Kokkos::deep_copy(repc_g, 0);
      const bool fixed_on = report_on && (rt_report_every > 0) &&
                            (pm->ncycle % rt_report_every == 0);
      const Real rep_r = rt_report_r;
      const int rep_k = ks, rep_j = js;
      const int rep_cyc = pm->ncycle;
      const Real rep_time = pm->time;
      const bool direct_on = rt_src_direct && (ck_on || grey_on);
      const bool semi_imp = rt_semi_implicit;
      // see rt_apply_debug: which column, and how many calls are left to print
      const bool dbg_on = (rt_apply_debug > 0);
      const int dbg_m = rt_dump_m;
      const int dbg_j = (rt_dump_j >= 0) ? rt_dump_j : (js + je)/2;
      const int dbg_k = (rt_dump_k >= 0) ? rt_dump_k : (ks + ke)/2;
      const int dbg_n = rt_apply_debug_n;
      if (dbg_on) --rt_apply_debug;
      // per-cycle diagnostic: an empty View captures fine, so the lambda needs no
      // branch on the pointer itself
      const bool diag = rt_diag;
      DvceArray5D<Real> dg;
      if (diag) dg = *rt_diag_ptr;
      par_reduce_clip4("rt_apply", 0, nmb1, ks, ke, js, je, is, ie, nclip,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nc) {
        Real Ft = 0.0, Fb = 0.0;
        for (int b=0; b<nblk; ++b) {
          Ft += Fb_g(m,b,i+1,k,j);
          Fb += Fb_g(m,b,i,k,j);
        }
        // the two-stream's share of each face in the tau blend
        Real src;
        if (direct_on) {
          // The direct source is the divergence of the FULL two-stream flux.  What the
          // two-stream actually deposits is the divergence of its blended share,
          //     -d[(1-w)F]/dr = -dF/dr + d[wF]/dr,
          // and the second term is the flux it hands to the conduction operator across
          // the blend layer -- an O(F/dr) exchange, not a correction.  Keep it, formed
          // from the raw face fluxes.  Where w = 0 (the thin layers this rewrite is for)
          // it vanishes identically, so no round-off of F - F re-enters there; where it
          // is non-zero the cell is at tau > 10 and its energy dwarfs that round-off.
          src = 0.0;
          for (int b=0; b<nblk; ++b) src += Src_g(m,b,i,k,j);
          if (taublend) {
            src += (w_g(m,k,j,i+1)*Ft - w_g(m,k,j,i)*Fb)/dx1(m,k,j,i);
          }
        } else {
          if (taublend) {
            Ft *= (1.0 - w_g(m,k,j,i+1));
            Fb *= (1.0 - w_g(m,k,j,i));
          }
          src = -(Ft-Fb)/dx1(m,k,j,i);
        }
        Real Qs_d = 0.0;   // the stellar heating that entered src, for the diagnostic
        if (band_on) {
          // deeper than the cut nothing radiative is applied: that region is optically
          // thick and convective, and the stellar beam died decades of optical depth
          // above
          if (i < icut_g(m,k,j)) {
            src = 0.0;
          } else {
            Real Qs = 0.0;
            for (int b=0; b<nblk; ++b) Qs += Qb_g(m,b,i,k,j);
            src += Qs;
            Qs_d = Qs;
          }
        } else {
          Qs_d = Qv_g(m,k,j,i);
          src += Qs_d;
        }
        // SEMI-IMPLICIT APPLICATION.  The source splits as src = A - E(T), A being the
        // absorption of the field from elsewhere, fixed on this step, and E the cell's
        // own emission, which is 4 sigma kappa_P rho T^4 and so scales as T^4.  The step
        // has to be right in two limits at once: it must reduce to the explicit rate
        // src*bdt as bdt -> 0, and it must land on the state where the cell stops
        // exchanging energy, E(T) = A, as bdt -> infinity.
        //
        // WHY NOT LINEARIZE.  Treating E implicitly by linearizing about the CURRENT
        // state gives a relaxation rate lambda = dE/de = (4E/T)/(rho c_v) ~ 4E/e and the
        // step (src/lambda)(1 - exp(-lambda bdt)).  That is only valid while |src| <~ E.
        // Its asymptote is e(A-E)/(4E), which for a cold optically thin cell in a hotter
        // field (A >> E) DIVERGES as E -> 0: the emission is damped but the absorption
        // stays explicit and unbounded.  A 400 K ambient cell sitting in a field whose
        // equilibrium temperature is 1500 K has A/E = (1500/400)^4 = 197 and is handed
        // de/e ~ 49 in one step, where the true answer is T_eq/T - 1 = 2.75 -- an
        // overshoot of ~18x, and the reason rt_de_max used to clip essentially every cell
        // above the photosphere.  Cooling was safe (A = 0 gives -e/4); heating was not.
        //
        // WHAT THIS DOES INSTEAD.  Relax toward the TRUE fixed point.  E ~ T^4 and, under
        // the same e ~ rho c_v T the old rate assumed, e ~ T, so E(T_eq) = A puts the
        // equilibrium at e_eq = e (A/E)^(1/4) and deq = e_eq - e.  Then take the rate
        // from the initial slope rather than from lambda -- x = src*bdt/deq -- so that
        //     de = deq (1 - exp(-x))
        // is exactly src*bdt for small bdt and exactly deq for large bdt.  src and deq
        // always share a sign (A > E <=> e_eq > e), so x >= 0 and the step is monotone
        // and can never cross the equilibrium from either side.  Where the departure is
        // small it agrees with the old linearized form, since deq -> src/lambda there.
        //
        // The e ~ T behind e_eq is the same approximation the old rate made, and it
        // UNDER-estimates e_eq wherever H2 or H is partly dissociated, i.e. it errs
        // toward a smaller step.  problem/rt_semi_lin recovers the old linearization.
        Real de = src*bdt;
        // rt_cell_report bookkeeping: the pieces of the step, kept for the report below
        Real dg_A = 0.0, dg_Em = 0.0, dg_deq = 0.0;
        Real dg_it[8];
        int dg_nit = 0;
        bool dg_resc = false;
        for (int q=0; q<8; ++q) dg_it[q] = 0.0;
        // problem/rt_explicit (or problem/rt_semi_implicit = false): nothing else
        // touches de.  See rt_explicit and rt_semi_implicit.
        if (!explicit_on && semi_imp) {
          Real Em = 0.0;
          for (int b=0; b<nblk; ++b) Em += Em_g(m,b,i,k,j);
          if (taublend) {
            Em *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (band_on && i < icut_g(m,k,j)) Em = 0.0;
          const Real ei = eiN(m,k,j,i);
          if (Em > 0.0 && ei > 0.0) {
            const Real sdt = src*bdt;
            if (semilin) {
              const Real lam = 4.0*Em/ei;
              const Real x = lam*bdt;
              de = (x > 1.0e-4) ? (src/lam)*(-expm1(-x)) : sdt;
            } else {
              const Real absn = src + Em;              // A, held fixed over the step
              // e_eq - e.  With nothing arriving the equilibrium is T = 0, i.e. -e.
              const Real deq = (absn > 0.0) ? ei*(sqrt(sqrt(absn/Em)) - 1.0) : -ei;
              dg_A = absn; dg_Em = Em; dg_deq = deq;
              if (deq != 0.0) {
                const Real x = sdt/deq;                // >= 0: src and deq share a sign
                de = (x > 1.0e-4) ? deq*(-expm1(-x)) : sdt;
              } else {
                de = sdt;
              }
              // NEWTON REFINEMENT, seeded by the closed form above.  That estimate is
              // already the right asymptote, so this only has to correct the e ~ T it
              // assumed -- with H2 dissociating or H ionizing, e(T) is far steeper than
              // linear and the equilibrium moves.  Solve the exact backward-Euler
              // balance with the emission at the NEW temperature,
              //     F(de) = de - A dt + E(T_old) (T_new/T_old)^4 dt = 0,
              // which is the same E ~ T^4 the band solver emits with, but with T(e) and
              // c_v(e) taken from the EOS.  F is monotone in de (both terms increase),
              // so Newton from a bracketing-quality guess converges in a step or two;
              // F(0) = -src dt recovers the explicit answer if it stops immediately.
              if (newton_on && eos.IsGeneral()) {
                const Real t0 = T_g(m,k,j,i);
                const Real d0 = rhoN(m,k,j,i);
                const Real abdt = absn*bdt, embdt = Em*bdt;
                if (t0 > 0.0 && d0 > 0.0) {
                  for (int it=0; it<8; ++it) {
                    const Real e1 = ei + de;
                    if (!(e1 > 0.0)) break;
                    const Real tc = eos.Temperature(d0, e1);
                    const Real t1 = tc*eos.temp_cgs;
                    if (!(t1 > 0.0)) break;
                    const Real cv = d0*eos.SpecificHeatCv(d0, e1, tc);
                    if (!(cv > 0.0)) break;
                    const Real r4 = SQR(SQR(t1/t0));
                    const Real fx = de - abdt + embdt*r4;
                    // d(r4)/de = 4 r4/T dT/de, with dT/de = temp_cgs/(d c_v)
                    const Real dfx = 1.0 + embdt*4.0*r4/t1*(eos.temp_cgs/cv);
                    if (!(dfx > 0.0)) break;
                    const Real step = fx/dfx;
                    de -= step;
                    if (dg_nit < 8) dg_it[dg_nit++] = de;
                    if (fabs(step) <= 1.0e-8*(fabs(de) + fabs(ei))) break;
                  }
                }
                // THE LOOP CHECKS e1 > 0 AT THE TOP, NOT AT THE BOTTOM.  The last
                // `de -= step` is never validated, so Newton can exit having pushed the
                // cell to ei + de <= 0 -- a one-step NaN with no counterpart in the
                // closed form, which guarantees e1 = ei*exp(-x) > 0.  R9's all-column
                // NaN at t = 1.99e5 is under bisection; this closes the only path in
                // this branch to a non-positive energy.  Falling back to a 99.9 %
                // drop keeps the cell cooling hard without ever crossing zero.
                if (!(ei + de > 0.0)) {
                  // problem/rt_rescue_eq: land on the equilibrium the cell actually
                  // sees rather than on a fixed 99.9 % drop.  deq is e_eq - e from the
                  // closed form above, i.e. Em(T_eq) = A = src + Em; the old floor is
                  // kept underneath it, so this can only make the rescue gentler.
                  const Real defl = -(1.0 - 1.0e-3)*ei;
                  if (resc_eq && deq < 0.0 && deq > defl) {
                    de = deq;
                    Kokkos::atomic_fetch_add(&efix_g(1), 1);
                  } else {
                    de = defl;
                    Kokkos::atomic_fetch_add(&efix_g(2), 1);
                  }
                  Kokkos::atomic_fetch_add(&efix_g(0), 1);
                  dg_resc = true;
                }
              }
            }
          }
        }
        const Real de_pre = de;
        if (demax > 0.0) {
          const Real dl = LimitRTSource(de, eiN(m,k,j,i), demax);
          if (dl != de) { ++nc; de = dl; }
        }
        if (diag) {
          Real Em_d = 0.0;
          for (int b=0; b<nblk; ++b) Em_d += Em_g(m,b,i,k,j);
          if (taublend) {
            Em_d *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (band_on && i < icut_g(m,k,j)) Em_d = 0.0;
          dg(m,0,k,j,i) = src;
          dg(m,1,k,j,i) = de;
          dg(m,2,k,j,i) = (de != de_pre) ? 1.0 : 0.0;
          dg(m,3,k,j,i) = Ft;
          dg(m,4,k,j,i) = Fb;
          dg(m,5,k,j,i) = Qs_d;
          dg(m,6,k,j,i) = Em_d;
        }
        if (dbg_on && m == dbg_m && k == dbg_k && j == dbg_j && i > ie - dbg_n) {
          Real Em = 0.0, Qs = 0.0;
          for (int b=0; b<nblk; ++b) { Em += Em_g(m,b,i,k,j); Qs += Qb_g(m,b,i,k,j); }
          // the raw direct source, before the blend handover and the beam are added
          Real srcraw = 0.0;
          for (int b=0; b<nblk; ++b) srcraw += Src_g(m,b,i,k,j);
          const Real ei = eiN(m,k,j,i);
          Kokkos::printf("rt_apply i=%d T=%.4e d=%.4e e=%.4e Fb=%.6e Ft=%.6e "
                         "divF=%.14e srcraw=%.14e src=%.4e Qs=%.4e Em=%.4e "
                         "lamdt=%.4e de=%.4e "
                         "de/e=%.4e tau=%.4e w=%.4e dx=%.4e\n",
                         i, T_g(m,k,j,i), rhoN(m,k,j,i), ei, Fb, Ft,
                         -(Ft-Fb)/dx1(m,k,j,i), srcraw, src, Qs, Em,
                         (Em > 0.0 && ei > 0.0) ? 4.0*Em/ei*bdt : 0.0,
                         de, de/ei, taublend ? tauf_g(m,k,j,i) : 0.0,
                         taublend ? w_g(m,k,j,i) : 0.0, dx1(m,k,j,i));
        }
        // ---- rt_cell_report ---------------------------------------------------
        if (report_on) {
          const bool fixedcell = fixed_on && (k == rep_k) && (j == rep_j) &&
              (fabs(x1v_(m,i) - rep_r) < 0.5*dx1(m,k,j,i));
          bool doprint = fixedcell;
          if (dg_resc) {
            if (Kokkos::atomic_fetch_add(&repc_g(0), 1) == 0) doprint = true;
          }
          if (doprint && i >= icut_g(m,k,j)) {
            const Real kap = kc_g(m,0,i,k,j);
            const Real rho = rhoN(m,k,j,i);
            const Real dtc = kap*rho*dx1(m,k,j,i);
            Real tautop = 0.0;
            for (int i2=i; i2<=ie; ++i2) {
              tautop += kc_g(m,0,i2,k,j)*rhoN(m,k,j,i2)*dx1(m,k,j,i2);
            }
            // the hemispheric-mean weights the sweeps used (ck_nquad = 1); with
            // ck_nquad = 2 the split into absorption and emission below is indicative
            const Real mu1 = 1.0/CK_DIFFUSIVITY, wf = M_PI;
            const Real xq = dtc/mu1;
            const Real e0 = -expm1(-xq);
            const Real alp = (xq > 1.0e-3) ? (e0 - 1.0 + e0/xq) : (xq/2.0 - SQR(xq)/3.0);
            const Real bet = (xq > 1.0e-3) ? (1.0 - e0/xq) : (xq/2.0 - SQR(xq)/6.0);
            const Real gmq = alp;
            const Real bi = Bb_g(m,0,i,k,j);
            const Real bip = Bb_g(m,0,i+1,k,j);
            const Real bim = (i > is) ? Bb_g(m,0,i-1,k,j) : bi;
            const Real idn = idn_g(m,k,j,i+1), iup = iup_g(m,k,j,i);
            const Real absdn = wf/dx1(m,k,j,i)*e0*idn;
            const Real absup = wf/dx1(m,k,j,i)*e0*iup;
            const Real emidn = wf/dx1(m,k,j,i)*(alp*bip + bet*bi);
            const Real emiup = wf/dx1(m,k,j,i)*(bet*bip + gmq*bi);
            const Real ei = eiN(m,k,j,i);
            // the staleness itself, so a report says whether w0 and u0 had parted
            const Real ei_w0_r = w0_uc_(m,IEN,k,j,i);
            const Real ei_u0_r = EintFromCons(u0_uc_, m, k, j, i,
                                              cs_uc_ ? cosc_uc_(m,k,j) : 0.0, cs_uc_,
                                              etg_uc_,
                                              etg_uc_ ? phicc_uc_(m,k,j,i) : 0.0);
            const Real rho_w0_r = w0_uc_(m,IDN,k,j,i);
            const Real rho_u0_r = u0_uc_(m,IDN,k,j,i);
            Kokkos::printf(
              "### rt_cell_report %s ncycle=%d t=%.6e m=%d k=%d j=%d i=%d r=%.6e\n"
              "    rho=%.4e T=%.6e kap=%.4e kaprho=%.4e dtau=%.4e tau_to_top=%.4e "
              "dt=%.4e\n"
              "    I_dn(i+1)=%.6e I_dn(i)=%.6e I_up(i)=%.6e I_up(i+1)=%.6e "
              "F(i)=%.6e F(i+1)=%.6e\n"
              "    B(i-1)=%.6e B(i)=%.6e B(i+1)=%.6e T(i-1)=%.4e T(i+1)=%.4e "
              "kap(i-1)=%.3e kap(i+1)=%.3e\n"
              "    abs_dn=%.6e abs_up=%.6e emi_dn=%.6e emi_up=%.6e src=%.6e Em=%.6e "
              "A=%.6e\n"
              "    e=%.6e deq=%.6e e_eq=%.6e de=%.6e de/e=%.6e rescued=%d nit=%d\n"
              "    ei_w0=%.6e ei_u0=%.6e stale=%.4e rho_w0=%.6e rho_u0=%.6e\n"
              "    newton de: %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
              dg_resc ? "RESCUE" : "FIXED", rep_cyc, rep_time, m, k, j, i, x1v_(m,i),
              rho, T_g(m,k,j,i), kap, kap*rho, dtc, tautop, bdt,
              idn, idn_g(m,k,j,i), iup, iup_g(m,k,j,i+1),
              Fb_g(m,0,i,k,j), Fb_g(m,0,i+1,k,j),
              bim, bi, bip, (i > is) ? T_g(m,k,j,i-1) : 0.0, T_g(m,k,j,i+1),
              (i > is) ? kc_g(m,0,i-1,k,j) : 0.0, kc_g(m,0,i+1,k,j),
              absdn, absup, emidn, emiup, src, dg_Em, dg_A,
              ei, dg_deq, ei + dg_deq, de, de/ei, dg_resc ? 1 : 0, dg_nit,
              ei_w0_r, ei_u0_r,
              (ei_w0_r != 0.0) ? (ei_u0_r - ei_w0_r)/ei_w0_r : 0.0,
              rho_w0_r, rho_u0_r,
              dg_it[0], dg_it[1], dg_it[2], dg_it[3],
              dg_it[4], dg_it[5], dg_it[6], dg_it[7]);
          }
        }
        u0(m,IEN,k,j,i) += de;
      });
      rt_nclip_last = nclip;
      RTSourceLimiterWarn(nclip);
      // The Newton positivity rescue should never fire.  Say so the first time it does,
      // with the running total, and stay quiet afterwards.
      {
        static bool efix_warned = false;
        static int efix_seen = 0;
        auto he = Kokkos::create_mirror_view(efix_g);
        Kokkos::deep_copy(he, efix_g);
        if (he(0) > efix_seen && !efix_warned) {
          std::cout << "### two_stream_rt: the Newton step left e <= 0 in " << he(0)
                    << " cell(s): " << he(1) << " rescued to the radiative equilibrium "
                    << "(problem/rt_rescue_eq), " << he(2) << " to the 99.9 % floor; "
                    << "this is reported once" << std::endl;
          efix_warned = true;
        }
        efix_seen = he(0);
      }

      // ---- one-shot column dump, for cross-code comparison -------------------------
      if (band_on && !rt_dump_file.empty() && !rt_dump_done) {
        rt_dump_done = true;
        const int jd = (rt_dump_j >= 0) ? rt_dump_j : (js + je)/2;
        const int kd = (rt_dump_k >= 0) ? rt_dump_k : (ks + ke)/2;
        const int md = (rt_dump_m <= nmb1) ? rt_dump_m : 0;
        DvceArray2D<Real> col("rt_col", n1, 9);
        par_for("rt_dumpcol", DevExeSpace(), is, ie+1, KOKKOS_LAMBDA(const int i) {
          Real Fs = 0.0;
          Real Qs = 0.0;
          for (int b=0; b<nblk; ++b) {
            Fs += Fb_g(md,b,i,kd,jd);
            Qs += Qb_g(md,b,i,kd,jd);
          }
          col(i,0) = pb_g(md,kd,jd,i);
          col(i,1) = T_g(md,kd,jd,i);
          col(i,2) = Fs;
          col(i,3) = Qs;
          col(i,4) = x1f_(md,i);
          // Gamma_1 and grad_ad of the CURRENT state. Both collapse to the ideal values
          // when the EOS is ideal; under the tabulated EOS they dip hard through the H2
          // dissociation and H ionization bands, and how steeply they vary across a cell
          // is what the reconstruction has to cope with.
          const Real dd = rhoN(md,kd,jd,i);
          const Real ee = eiN(md,kd,jd,i);
          col(i,5) = eos.IsGeneral() ? eos.Gamma1(dd, ee) : eos.gamma;
          col(i,6) = GradAd(eos, eos.gamma, Rgas, pb_g(md,kd,jd,i)*1.0e6,
                            T_g(md,kd,jd,i));
          col(i,7) = taublend ? tauf_g(md,kd,jd,i) : 0.0;
          col(i,8) = taublend ? w_g(md,kd,jd,i) : 0.0;
        });
        auto hc = Kokkos::create_mirror_view(col);
        Kokkos::deep_copy(hc, col);
        auto hcut = Kokkos::create_mirror_view(icut_g);
        Kokkos::deep_copy(hcut, icut_g);
        auto hcf = Kokkos::create_mirror_view(cf_g);
        Kokkos::deep_copy(hcf, cf_g);
        if (global_variable::my_rank == 0) {
          std::ofstream f(rt_dump_file);
          f.precision(10);
          f << std::scientific;
          f << (grey_on ? "# two_stream_rt GREY column dump\n"
                        : "# deep_hot_jupiter_rt correlated-k column dump\n")
            << "# meshblock " << md << ", k = " << kd << ", j = " << jd
            << ", mu0 = " << hcf(md,kd,jd,3) << ", icut = " << hcut(md,kd,jd)
            << " (is = " << is << ", ie = " << ie << ")\n"
            << "# T_int = " << Tint << " K, T_irr = " << Tirr
            << " K, grav = " << grav << " cm/s^2\n"
            << "# fluxes are cgs: 1 erg/s/cm^2 = 1e-3 W/m^2. F_lw is NET (up minus "
                << "down)\n"
            << "# i  r_face[cm]  p[bar]  T[K]  F_lw_net[erg/s/cm2]  Q_sw[erg/s/cm3]"
            << "  Gamma_1  grad_ad  tau_R(face)  w_diff(face)\n";
          for (int i=is; i<ie+2; ++i) {
            f << i << " " << hc(i,4) << " " << hc(i,0) << " " << hc(i,1)
              << " " << hc(i,2) << " " << hc(i,3)
              << " " << hc(i,5) << " " << hc(i,6)
              << " " << hc(i,7) << " " << hc(i,8) << "\n";
          }
          f.close();
          std::cout << "two_stream_rt: wrote "
                    << (grey_on ? "grey" : "correlated-k") << " column dump to '"
                    << rt_dump_file << "' (k = " << kd << ", j = " << jd
                    << ", mu0 = " << hcf(md,kd,jd,3) << ")" << std::endl;
        }
      }
      return;
    }

//    size_t scr_size = 8 * ScrArray1D<Real>::shmem_size(n1);
//    int scr_level = 0;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, scr_level,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    // The private radial arrays here were a flat 270, with no check. Production runs at
    // nx1 = 256, i.e. n1 = 260, so they were ten cells from silently overrunning per-
    // thread
    // memory -- and nghost = 4 with the same nx1 would have gone over. Dispatched at run
    // time over compile-time tiers instead, as the correlated-k chain kernel is.
    int nclip_grey = 0;
    const Real demax_grey = rt_de_max;
    // --- problem/nan_report (see rt_nan_report): the in-kernel catcher.  One int
    // counter and one 16-slot record, zeroed per call; the FIRST offending cell on the
    // rank wins the record.  Everything it needs is recomputed inside the `bad` branch,
    // so the fast path pays nothing but a bool test.
    const bool nanrep_g = rt_nan_report;
    if (nanrep_g && rt_nanrep_cnt == nullptr) {
      rt_nanrep_cnt = new DvceArray1D<int>("rt_nanrep_cnt", 1);
      rt_nanrep_rec = new DvceArray1D<Real>("rt_nanrep_rec", 16);
    }
    auto nrcnt = nanrep_g ? *rt_nanrep_cnt : DvceArray1D<int>("d", 1);
    auto nrrec = nanrep_g ? *rt_nanrep_rec : DvceArray1D<Real>("d", 1);
    if (nanrep_g) {
      Kokkos::deep_copy(nrcnt, 0);
      Kokkos::deep_copy(nrrec, 0.0);
    }
    auto launch_grey_rt = [&](auto nn_tag) {
      constexpr int NN = decltype(nn_tag)::value;
      par_reduce_clip3("2stream_rt", 0, nmb1, ks, ke, js, je, nclip_grey,
      KOKKOS_LAMBDA(const int m, const int k, const int j, int &nc) {
  //        ScrArray1D<Real> tau_down_r_f(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> F_v_down_f(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> B(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> I_ir_down_f(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> I_ir_up_f(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> F_ir_f(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> Q_v(member.team_scratch(scr_level), n1);
  //        ScrArray1D<Real> kapJ_ir(member.team_scratch(scr_level), n1);
          Real tau_down_r_f[NN];
  //        Real F_v_down_f[NN];
          Real B[NN];
          Real F_ir_f[NN];
          Real Q_v[NN];
  //        Real kapJ_ir[NN];

          Real x2v = x2v_(m,j);
          Real x3v = x3v_(m,k);

          Real rtop = x1v_(m,ie+1);
          Real rbot = x1v_(m,is);

          Real lam, phi, theta;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -theta+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            lam = x3v*iap;
            theta = -lam+M_PI/2.0;
            phi = x2v*iap;
          }
          Real ex = sin(theta)*cos(phi);
          Real ex0 = 1.0;
          Real mu0 = ex*ex0;
          if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);

          Real mus = (mu0 > 0.0) ? mu0 : 0.0;
          Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
          Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
          get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);

          // 3 V Bands
          // top
          Real rho = rhoN(m,k,j,ie+1);
          Real p = PresFromEint(eos,gm1,rho,eiN(m,k,j,ie+1));
          Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,ie+1),p);
          bool badtop = RTBadState(p, T);    // see RTBadState: one NaN kills the column
          if (badtop) {
            Kokkos::atomic_fetch_add(&stcl_g(0), 1);
            T = 0.0;
            p = 0.0;
          }
          B[ie+1] = badtop ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr = 0.0;
          if (!badtop) get_kapr(T, p, met, kapr);
          Real tau_r_f = RTTopDtau(kapr, p,
                                   EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0,
                                             tide));
          tau_down_r_f[ie+1] = tau_r_f;
          Real drtop = (kapr*rho != 0.0) ? tau_r_f/(kapr*rho) : 0.0;
          Real delta = drtop/rtop;
          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
          fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          Real tausl = tau_r_f*fac;
          Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
          Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
          Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
          // beam transmission at the face ABOVE the cell being filled, carried down the
          // sweep so that differencing the flux across a cell costs no extra exp
          Real trp1 = trans1;
          Real trp2 = trans2;
          Real trp3 = trans3;
  //        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
  //        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
          // down-sweep
          for (int i=ie; i>is-1; --i) {
            Real rho = rhoN(m,k,j,i);
            Real p, T;
            PresTempFromEint(eos,gm1,Rgas,rho,eiN(m,k,j,i),
                             TGuess(wtemp_, m, k, j, i),p,T);
            bool badcell = RTBadState(p, T); // see RTBadState: one NaN kills the column
            if (badcell) {
              Kokkos::atomic_fetch_add(&stcl_g(0), 1);
              T = 0.0;
              p = 0.0;
            }
            B[i] = badcell ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
            Real kapr = 0.0;
            if (!badcell) get_kapr(T, p, met, kapr);
            Real dr = dx1(m,k,j,i);
            tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
            Real r = x1f_(m,i);
  ////          Real delta = (drtop+(rtop-r))/r;
  //          Real delta = dr/r;
  //          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
  //          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
            Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
            Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
            Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
            Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
  //          Real trans1 = exp(-gamv1*tausl);
  //          Real trans2 = exp(-gamv2*tausl);
  //          Real trans3 = exp(-gamv3*tausl);
  //          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
            Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
            // Deposit the flux DIFFERENCE across the cell. The old form,
            // kappa rho F exp(-tau) with tau at the lower face, is right only for a thin
            // layer: it returns u e^-u / (1 - e^-u) of what the cell actually absorbs,
            // with
            // u = dtau/mu, which is 0.95 at u = 0.1 but 0.58 at u = 1. On this grid,
            // about
            // 0.46 scale heights per cell, that put dtau/mu near one wherever it mattered
            // and lost about 24 % of the incident stellar flux. Found by comparing the
            // correlated-k version of the same expression against Exo-FMS on an identical
            // column. Written this way the column integral telescopes to
            // mu F (1 - e^-tau_total) exactly, and it still reduces to the old expression
            // as dtau -> 0.
            Real Qv = (1.0-albedo)*Fstar*(1.0/3.0)
                    * ((trp1-trans1)+(trp2-trans2)+(trp3-trans3))/(fac*dr);
            Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
            trp1 = trans1;
            trp2 = trans2;
            trp3 = trans3;
          }

          // 2 IR Bands x two quadrature points, interleaved, in blocks of NC.
          //
          // Each (band, quadrature) combination is an independent pair of linear
          // recurrences in radius, and running them one after another leaves the
          // wavefront
          // stalled on a single dependency chain: this kernel has only nmb*nx3*nx2/64
          // wavefronts for 912 SIMDs, so there is no other wave to hide that latency and
          // VALUBusy sits near 3 %. Stepping NC combinations inside one radial loop gives
          // the chain NC independent strands to overlap, and blocking keeps the private
          // I_ir_down_c footprint at NC columns however many chains are requested.
          //
          // This is the grey path, so nchain_rt is 4 and there is exactly one block; the
          // blocking survives because the correlated-k kernel shares the structure.
          constexpr int NC = RT_NB;
          const int nblk_rt = (nchain_rt + NC - 1)/NC;
          for (int i=is; i<ie+2; ++i) {
            F_ir_f[i] = 0.0;
          }

          for (int blk=0; blk<nblk_rt; ++blk) {
            Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC];
            for (int cc=0; cc<NC; ++cc) {
              const int c = blk*NC + cc;
              const int n = (c/2) % 2;
              const int vir = c % 2;
              muggc[cc] = mug[n];
              wggc[cc] = wg[n];
              gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
              fbc[cc] = (vir == 0) ? beta : (1.0-beta);
            }
            Real I_ir_down_c[NC][NN];

            // top
            for (int cc=0; cc<NC; ++cc) {
              Real dtauir = gamirc[cc]*tau_down_r_f[ie+1];
              Real trans = exp(-dtauir/muggc[cc]);
              I_ir_down_c[cc][ie+1] = (1.0-trans)*(fbc[cc]*B[ie+1]);
            }
            // down-sweep
            for (int i=ie; i>is-1; --i) {
              Real dtau_i = tau_down_r_f[i]-tau_down_r_f[i+1];
              for (int cc=0; cc<NC; ++cc) {
                Real dtauir = gamirc[cc]*dtau_i;
                Real x = dtauir/muggc[cc];
                Real e0 = -expm1(-x);
                Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
                Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
                I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                                   + alp*fbc[cc]*B[i+1] + bet*fbc[cc]*B[i];
              }
            }

            // bottom
            Real I_ir_up_c[NC];
            for (int cc=0; cc<NC; ++cc) {
              I_ir_up_c[cc] = Iint + I_ir_down_c[cc][is];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is];
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
              F_ir_f[is] += (F_ir_up_f - F_ir_down_f);
            }
            // up-sweep, accumulating the band flux as it goes
            for (int i=is+1; i<ie+2; ++i) {
              Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
              for (int cc=0; cc<NC; ++cc) {
                Real dtauir = gamirc[cc]*dtau_i;
                Real x = dtauir/muggc[cc];
                Real e0 = -expm1(-x);
                Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
                Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
                I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]
                              + bet*fbc[cc]*B[i] + gm*fbc[cc]*B[i-1];
                Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
                Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
                F_ir_f[i] += (F_ir_up_f - F_ir_down_f);
              }
            }
          }

  //        // Sync all threads in the team so that scratch memory is consistent
  //        member.team_barrier();

  //        par_for_inner(member, is, ie, [&](const int i) {
          for (int i=is; i<ie+1; ++i) {
            // source term as flux divergence
            Real area_t = area1(m,k,j,i+1);
            Real area_b = area1(m,k,j,i);
            Real vol = volume(m,k,j,i);
              Real Ft = F_ir_f[i+1];//-F_v_down_f(i+1);
              Real Fb = F_ir_f[i];//-F_v_down_f(i);
            Real src = -(Ft-Fb)/dx1(m,k,j,i);
            if (correct_spherical) {
              src = -(Ft*area_t-Fb*area_b)/vol;
            }
              src += Q_v[i];
            Real du_flux = src*bdt;

  //          // source term semi-implicit
  //          Real p = PresFromEint(eos,gm1,rhoN(m,k,j,i),eiN(m,k,j,i));
  //          Real rho = rhoN(m,k,j,i);
  //          Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,i),p);
  //          Real kapr;
  //          get_kapr(T, p, met, kapr);
  //          Real cv = Rgas*rho*igm1;
  //          Real e0 = eos.IsGeneral() ? eiN(m,k,j,i) : cv*T;
  //          Real kk = 0.0;
  //          Real bb = du_flux + e0;
  ////          Real bb = Q_v(i)*bdt + e0;
  //          for (int vir=0; vir<2; ++vir) {
  //            Real gamir, fb;
  //            if (vir == 0) {
  //              gamir = gamir1;
  //              fb = beta;
  //            } else {
  //              gamir = gamir2;
  //              fb = 1.0-beta;
  //            }
  //            kk += -4.0*M_PI*gamir*kapr*rho*fb*boltz_sigma/M_PI*bdt;
  //            bb += 4.0*M_PI*gamir*kapr*rho*fb*B[i]*bdt;
  //          }
  ////          bb += 4.0*M_PI*rho*kapJ_ir(i)*bdt;
  //          int ierr=0;
  //          Real e = e0;
  //          // Newton-Raphson. A general EOS has no e = c_v T with constant c_v, so the
  //          // iteration runs on the internal energy directly rather than on T:
  //          // F(e) = e - kk T(e)^4 - bb, with dT/de = temp_cgs/(d c_v) since T is in K.
  //          for (int n=0; n<100; ++n) {
  //            Real de;
  //            if (eos.IsGeneral()) {
  //              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e));
  //              de = e - kk*SQR(SQR(T)) - bb;
  //              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
  //              T = eos.Temperature(rho,e)*eos.temp_cgs;
  //            } else {
  //              e = cv*T;
  //              de = e - kk*SQR(SQR(T)) - bb;
  //              T -= de / (cv - 4.0*kk*T*T*T);
  //            }
  //            if (T < 0.0) {
  //              e = e0;
  //              ierr = 1;
  //              break;
  //            }
  //            if (fabs(de) <= 1.0e-10*e)
  //              break;
  //          }
  //          Real du_src = e-e0;
  //
  //          Real du = (fabs(du_flux) < e0 && ierr == 1) ? du_flux : du_src;
            Real du = du_flux;
            if (demax_grey > 0.0) {
              const Real dl = LimitRTSource(du, eiN(m,k,j,i), demax_grey);
              if (dl != du) { ++nc; du = dl; }
            }
            u0(m,IEN,k,j,i) += du;
            // --- the nan_report catcher.  Record the cell the grey apply just made
            // non-finite or non-positive, with the state that produced it.
            if (nanrep_g) {
              const Real enew = u0(m,IEN,k,j,i);
              if (!(enew > 0.0) || !isfinite(enew) || !isfinite(du)) {
                if (Kokkos::atomic_fetch_add(&nrcnt(0), 1) == 0) {
                  Real pc, tc, pm1, tm1, pp1, tp1, kc;
                  PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i), eiN(m,k,j,i),
                                   TGuess(wtemp_, m, k, j, i), pc, tc);
                  PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i-1), eiN(m,k,j,i-1),
                                   TGuess(wtemp_, m, k, j, i-1), pm1, tm1);
                  PresTempFromEint(eos, gm1, Rgas, rhoN(m,k,j,i+1), eiN(m,k,j,i+1),
                                   TGuess(wtemp_, m, k, j, i+1), pp1, tp1);
                  get_kapr(tc, pc, met, kc);
                  nrrec(0) = static_cast<Real>(m);
                  nrrec(1) = static_cast<Real>(k);
                  nrrec(2) = static_cast<Real>(j);
                  nrrec(3) = static_cast<Real>(i);
                  nrrec(4) = rhoN(m,k,j,i);
                  nrrec(5) = eiN(m,k,j,i);
                  nrrec(6) = tc;
                  nrrec(7) = kc*rhoN(m,k,j,i);
                  nrrec(8) = du;
                  nrrec(9) = du_flux;
                  nrrec(10) = enew;
                  nrrec(11) = eiN(m,k,j,i-1);
                  nrrec(12) = tm1;
                  nrrec(13) = eiN(m,k,j,i+1);
                  nrrec(14) = tp1;
                  nrrec(15) = tau_down_r_f[i];
                }
              }
            }
          }
  //        });

  //        // Sync all threads in the team so that scratch memory is consistent
  //        member.team_barrier();

      });
    };
    if (n1 <= 72) {
      launch_grey_rt(std::integral_constant<int, 72>{});
    } else if (n1 <= 136) {
      launch_grey_rt(std::integral_constant<int, 136>{});
    } else if (n1 <= 264) {
      launch_grey_rt(std::integral_constant<int, 264>{});
    } else if (n1 <= 520) {
      launch_grey_rt(std::integral_constant<int, 520>{});
    } else if (n1 <= 1032) {
      launch_grey_rt(std::integral_constant<int, 1032>{});
    } else {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: n1 = " << n1
                << " exceeds the largest radial tier (1032). Add a tier to the dispatch "
                << "in picket_fence_two_stream_RT." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    RTSourceLimiterWarn(nclip_grey);
    if (nanrep_g && rt_nanrep_lines < rt_nanrep_maxlines) {
      auto hgc = Kokkos::create_mirror_view(nrcnt);
      Kokkos::deep_copy(hgc, nrcnt);
      if (hgc(0) > 0) {
        auto hgr = Kokkos::create_mirror_view(nrrec);
        Kokkos::deep_copy(hgr, nrrec);
        ++rt_nanrep_lines;
        std::cout << "### rg nan_report [RT_grey_apply] rank "
                  << global_variable::my_rank << " cycle " << pm->ncycle
                  << " t = " << pm->time << ": " << hgc(0) << " cell(s); first (m,k,j,i)"
                  << " = (" << static_cast<int>(hgr(0)) << ","
                  << static_cast<int>(hgr(1)) << "," << static_cast<int>(hgr(2)) << ","
                  << static_cast<int>(hgr(3)) << ")"
                  << " d = " << hgr(4) << " e_pre = " << hgr(5) << " T = " << hgr(6)
                  << " rho_kap = " << hgr(7) << " du = " << hgr(8)
                  << " du_flux = " << hgr(9) << " u_new = " << hgr(10)
                  << " e[i-1] = " << hgr(11) << " T[i-1] = " << hgr(12)
                  << " e[i+1] = " << hgr(13) << " T[i+1] = " << hgr(14)
                  << " tau = " << hgr(15) << std::endl;
      }
    }

    return;
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_RT_HPP_
