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
#include <cstdio>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "coordinates/cell_locations.hpp"
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
#include "utils/rad_taper.hpp"
#include "utils/two_stream_column_implicit.hpp"
#include "utils/two_stream_column_partition.hpp"

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

// problem/rt_plane_parallel: run the GREY sweep on a PLANE-PARALLEL Cartesian mesh.
//
// Everything else in this file assumes a RADIAL mesh: the cell centres and faces come
// from Coordinates::x1v / xx1f and the radial width from Coordinates::dx1, and all three
// are 1x1 placeholder Views on a Cartesian mesh (coordinates.cpp reallocates them only
// for spherical-polar and cubed-sphere grids), so reading them there is an out-of-bounds
// access.  With this flag the three are taken from the MeshBlock's own RegionSize
// instead -- uniform dx1, centres from CellCenterX, faces from LeftEdgeX -- which is
// exactly what a plane-parallel column needs and costs the radial path nothing: the
// substitutions are runtime branches on a `const bool` that is false in every existing
// run, and the spherical expressions are untouched.
//
// The sweep needs no other geometric change, because with this flag on the area and
// volume helpers of the radial path (AFC/ACC/VLS/VLA below, and the SPHERICAL
// DILUTION note that introduces them) collapse to A = 1 and V = dx1, which is exactly
// the plane-parallel operator -(F_top - F_bot)/dx1 this path used to carry
// everywhere.  What DOES have to be right is the
// gravity model -- a plane-parallel box has constant g, so the caller sets
// grav_point_mass = false and stellar_tide = false, and EffGravAt then returns g
// whatever radius it is handed.
//
// RESTRICTED to the grey split path (rt_grey && rt_split): the monolithic picket-fence
// kernel and the correlated-k kernel carry stellar-beam geometry (the substellar angle,
// the slant path) that has no plane-parallel meaning.  The guard is a fatal error.
inline bool rt_plane_parallel = false;

//----------------------------------------------------------------------------------------
// problem/rt_rad_force -- the RADIATIVE MOMENTUM SOURCE that goes with the EOS's
// thin-region radiation taper (<eos>/eos_rad_rho_hi, see utils/rad_taper.hpp).
//
// WHY.  With the taper on, the hydro's pressure carries only w(rho) of the LTE radiation
// pressure, so the momentum equation loses (1-w) of -grad Prad -- and in the thin layer
// that is precisely where the real force lives: MEASURED on the relaxed B-star column,
// grad Prad/(rho g) is 0.10 at tau = 1e-2, 0.15 at tau = 0.1, 0.23 at tau = 2/3 and 0.30
// at tau = 10, and at the very top it agrees with kappa_R F/(c g) to 1.8x and below tau 1
// to a few per cent.  Dropping it would make the top of the box lighter than it is.
//
// WHAT IS ADDED, per cell:
//     f = (1 - w) rho kappa_R F_net/c   +   Prad grad w,     Prad = a T^4/3,
// with F_net the face-averaged net two-stream flux (positive upward, so the force is
// along +x1, against gravity).  The FIRST term is the direct momentum deposition the
// pressure gradient no longer supplies.  The SECOND is a correction, not a force: with
// p = pgas + w(rho) Prad the hydro differences
//     -grad(w Prad) = -w grad Prad - Prad grad w,
// and the last piece is an artefact of the taper: inside the ramp dw/dln rho is of
// order one, so it is of the same order as grad Prad itself, up to ~20 % of g.  Adding
// Prad grad w back cancels it exactly, and the total radiative force is then -grad Prad
// in the deep limit (w = 1), kappa rho F/c in the thin limit (w = 0), and continuous
// through the ramp.  The x1 component of the first term is the only one the column-wise
// two-stream knows; the second is a genuine pressure gradient and is applied in all three
// directions.  The work v.f is added to the total energy alongside.
//
// The step is EXPLICIT and lives in the same place as the tau-blend handover: after the
// semi-implicit relaxation, which must not damp it.  Default OFF; box_convection requires
// the EOS taper to be on with it.
inline bool rt_rad_force = false;
// problem/rt_force_center (default 0, bitwise off; needs rt_implicit_column = 3): WHICH
// face flux the radiative momentum source above is driven by.  In mode 3 the energy is
// solved exactly and implicitly over the stage, but Fb is deliberately left at the
// ENTRY sweep's (lambda-iteration) values, so the force alone still responds to the
// stage-START radiation field: for a radiatively relaxed surface mode that is a phase
// lag linear in dt, and the work f.v it does over a cycle is then an O(dt) forcing.
//   0 = the entry sweep's flux (as now),
//   1 = the flux of THIS stage's converged column solve (the column writes its own face
//       flux, exactly as rt_col3_skip_sweep already makes it do),
//   2 = the average of the two, i.e. a trapezoidal centring in the stage.
// With the switch on the column writes Fb, so rad_f2s, the emergent-flux history and the
// surface dumps read the CONVERGED flux too; the force itself is what mode 1/2 differ in.
inline int rt_force_center = 0;
inline DvceArray4D<Real> *rt_fbsave_ptr = nullptr;  // rt_force_center: the entry Fb sum
// problem/rt_src_theta (default 1.0, bitwise off; needs rt_implicit_column = 3): the
// time centring of the column's own energy deposit.  The solve is backward Euler within
// the stage (theta = 1, bitwise the old code), which lags T' by O(bdt/t_rad); theta < 1
// blends that deposit with the EXACT exponential relaxation of the same implied rate, in
// factor space -- see RTCol3ThetaDe in two_stream_column_implicit.hpp for why a literal
// theta*implicit + (1-theta)*explicit blend is unusable (it blew the He column up by a
// factor 1e5 at theta = 0.5).  theta = 0.5 is the Crank-Nicolson-like centring, theta = 0
// the exact linear relaxation.  Exact only in the linear limit, and NOT conservative
// against the column's own flux divergence: a diagnostic switch, not a production one.
// Needs the entry sweep, so it turns problem/rt_col3_skip_sweep off.
inline Real rt_src_theta = 1.0;
// problem/rt_budget_verbose (box_convection.cpp): when non-null, the radiative force's
// WORK term v.f is accumulated, box-integrated over the step, into slot 10 of this
// array.  Diagnostic only: nothing here changes a source term.
inline DvceArray1D<Real> *rt_bud_ptr = nullptr;
// problem/rt_src_dump: dump the grey sweep's per-cell source assembly for ONE column
// (m = 0, k = ks, j = js) on the next N apply-kernel calls, every cell from the band cut
// to the top: the raw face fluxes, both forms of the divergence, the assembled source
// and the source actually applied.  Diagnostic only.
inline int rt_src_dump = 0;
// problem/rt_force_verbose -- print the hydrostatic balance of every cell inside the
// taper ramp for this many RT calls, then stop.  a_p + a_g + a_f normalised by g: inside
// the ramp this is what says whether the Prad grad w correction is doing its job.
inline int rt_force_verbose = 0;
inline Real rt_force_grav = 0.0;    // |g| used only to normalise that print
// problem/rt_top_vacuum: nothing above the top of the domain.  The default top boundary
// is the UNRESOLVED HYDROSTATIC COLUMN, of optical depth kappa p / g_eff, which is the
// right model for a star whose atmosphere continues above x1max.  A local box is cut out
// of a stratification and the cells above it are not part of the problem; more to the
// point, a box with weak gravity would be handed p/g = a huge column and sealed.  With
// this flag the downward intensity entering the top face is exactly zero, which is the
// standard grey-atmosphere boundary and the one solar_convection.cpp's local solver uses.
inline bool rt_top_vacuum = false;

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
// ...and the STATE of the first cell it happened to, which is the only way to tell a
// genuine stiff-cooling cell from a poisoned neighbour without re-running.  Filled once
// per run by whichever cell claims efix(0) == 0, printed with the warning.
inline DvceArray1D<Real> *rt_efix_rec = nullptr;
// the column energy budget of the source application, accumulated when rt_cell_report is
// on: slot 0 is the EXPLICIT deposit sum(src*bdt*dx1), which telescopes to F_bot - F_top
// over a column, and slot 1 the deposit actually applied, sum(de*dx1).  The two differ by
// whatever the relaxation, the sub-cycle and LimitRTSource did, so their ratio is the
// conservation statement for problem/rt_relax_sub.
inline DvceArray1D<Real> *rt_desum_ptr = nullptr;
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
// problem/rt_relax_sub: SUB-CYCLE the local relaxation of the two-stream source.  The
// single-step form damps each cell toward the equilibrium of the un-relaxed column by a
// per-cell factor (1 - exp(-x))/x; neighbours at different x are damped by different
// amounts, so the non-local exchange between them no longer cancels and what is left is
// a pressure perturbation ~ dt x the compression rate, which pumps the standing acoustic
// modes of a closed box (the He-star box's saturated v_rms was LINEAR in dt).  With this
// set > 1 the local balance is integrated in nsub sub-steps of bdt/nsub instead, the
// absorbed field A frozen for the stage and the cell's own emission re-formed from the
// running energy after each sub-step, so every sub-step is taken at small x.  The value
// is the FLOOR on nsub; the stiffness rule below can raise it.  DEFAULT 1 = the
// single-step form, BIT FOR BIT.
// problem/rt_ali_diag: the ACCELERATED-LAMBDA (ALI) DIAGONAL in the per-cell
// semi-implicit apply.  Without it the step relaxes the cell toward A, the absorption
// the sweep computed, at the rate 4E/e set by the cell's OWN emission -- and in an
// optically thick cell that rate is enormous, so x = src*bdt/deq reaches 1e4-1e5 and the
// factor (1 - exp(-x))/x delivers 1/x of the source to the gas.  MEASURED on the B-star
// relaxed column (bench/bstar_fecz/leakdiag/dump0, arm d, w = 0 everywhere, dtau 5-134):
// the assembled source sums to +0.2931 F_bot = F_bot - F_top exactly, and the APPLIED
// sum is +0.000089 F_bot.  Nothing is clipped, no guard fires, the Newton converges in
// 1-2 steps -- the source is simply damped away, cell by cell, and the energy the sweep
// took out of the radiation field is never given to the gas.  That is the 0.29 F_bot
// deficit: Ftop/F_bot froze at 0.707 because the column can never re-relax.
//
// The mistake is treating A as an EXTERNAL field.  In a thick cell almost all of A is
// the cell's own emission coming back from its immediate neighbours, so the true
// response of the net exchange to a change in this cell's source function is not 4E/e
// but (1 - Lambda*_ii) 4E/e, with Lambda*_ii the diagonal of the Lambda operator.  That
// is textbook accelerated Lambda iteration; the sub-cycle and the Newton refinement are
// both fixes to the wrong rate rather than to the rate itself.
//
// Lambda*_ii here is the local escape-probability diagonal the sweep's own layer
// coefficients already carry, 1 - (1 - e^-x)/x at x = dtau/mu, averaged over the
// quadrature.  It goes to 1 as dtau -> infinity, which sends the relaxation rate to zero
// and the step to the EXACT explicit one, de = src*bdt -- conservative, and stable
// because a thick cell's net imbalance relaxes on the diffusion time, not the thermal
// one -- and to 0 as dtau -> 0, where the old stiff relaxation is what is wanted and is
// recovered bitwise.  rt_ali_diag = false restores the old behaviour exactly.
inline bool rt_ali_diag = true;
inline int rt_relax_sub = 1;
// problem/rt_relax_xcrit: the per-sub-step stiffness the sub-cycle aims for.  nsub is
// ceil(x/rt_relax_xcrit) with x = src_relax*bdt/deq the first step's stiffness, floored
// at rt_relax_sub.  Inert unless rt_relax_sub > 1.
inline Real rt_relax_xcrit = 1.0;
// problem/rt_relax_submax: the cap on nsub, so a single pathological cell cannot cost an
// unbounded number of EOS calls.  Inert unless rt_relax_sub > 1.
inline int rt_relax_submax = 32;
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
// problem/rt_outer_iter: how many times the WHOLE two-stream source step is repeated
// per stage.  Default 1 = the behaviour of every run so far, and bitwise so.
//
// WHAT IT TESTS.  With one pass the absorbed field A that each cell relaxes toward is
// formed from ONE sweep of the UN-relaxed column: the intensities are those of the state
// at the start of the stage, while the cell's own energy is then moved by the full stage
// increment.  The non-local half of the exchange therefore lags the local half by one
// stage, and that lag is proportional to bdt -- exactly the signature of the He-star
// box's dt-LINEAR saturated v_rms.  Sub-cycling (rt_relax_sub) cannot see this: it
// re-forms the cell's own emission but keeps A frozen, which is the whole point of it.
//
// WHAT IT DOES.  Pass k re-runs the column sweep on the CURRENT running state -- the
// partially relaxed column, e = e^n + de_{k-1} -- and then recomputes the TOTAL stage
// increment de_k from the ORIGINAL e^n against that updated field.  It is a fixed-point
// iteration for the implicit balance e^{n+1} = e^n + bdt S(e^{n+1}), not an accumulation
// of extra increments: de_k replaces de_{k-1} rather than adding to it, and u0 carries
// only the difference.  The handover term src_ex is re-formed from the updated face
// fluxes each pass too, since it is explicit and exact in flux form and should therefore
// use the CONVERGED field.  k passes cost ~k times the RT time.
//
// REQUIREMENTS (checked at the top of the wrapper).  The sweep has to be able to SEE the
// running state, so it needs problem/rt_use_cons (u0 is what the passes update; w0 is
// stale until the next ConToPrim), the grey split path, and the semi-implicit apply.
inline int rt_outer_iter = 1;
// the running total stage increment de_k, (m,k,j,i); allocated only when rt_outer_iter>1
inline DvceArray4D<Real> *rt_deacc_ptr = nullptr;
// the fixed-point convergence of the pass: (max |de_k - de_{k-1}|/|de_k|, max |de_k|)
inline DvceArray1D<Real> *rt_oconv_ptr = nullptr;
// problem/rt_outer_verbose: print that convergence every rt_report_every cycles even
// without problem/rt_cell_report.  Inert unless rt_outer_iter > 1.
inline bool rt_outer_verbose = false;
// ---- problem/rt_implicit_column: THE MERGED IMPLICIT COLUMN SOLVE -------------------
//
// WHAT IS WRONG WITH THE PER-CELL RELAXATION.  The grey split sweep is exactly LINEAR in
// the cell Planck functions B_j at frozen opacity: the net source is
//     Src_i = kappa_i rho_i [(M B)_i + g_i] - 4 pi kappa_i rho_i B_i,
// M the (dense, but strongly banded) exchange matrix the two sweeps build and g the
// boundary terms.  The apply block then relaxes each cell SEPARATELY toward the
// equilibrium of that one sweep, damping it by its own factor (1 - e^-x)/x.  Two
// neighbours exchanging O(1e3) F of radiation with a net of O(1) F are damped by
// DIFFERENT factors, so what survives the near-cancellation is not the net but a
// residual proportional to dt -- a pump, not a relaxation.  That is the measured
// dt-LINEAR saturated v_rms of the He-star box (5.29e4 / 1.61e5 / 2.82e5 cm/s over cfl
// 0.15 / 0.30 / 0.45), and lagging the neighbours in an outer fixed point only contracts
// it at ~0.9 per pass, which is useless.
//
// WHAT THIS DOES.  Per column the implicit balance is
//     F(T) = C_v (T - T*)/(beta dt) - kappa rho (M - I) B(T) = 0,
// whose Jacobian J = C_v/(beta dt) + kappa rho (I - M) dB/dT is an M-matrix.  Its
// NEAREST-NEIGHBOUR part is accumulated during the sweep itself, out of the e0/alp/bet
// layer quantities already in registers, and folded straight into the radial implicit
// conduction tridiagonal (Conduction::ImplicitRadialUpdate), which is solved for the
// same column and is already an M-matrix in dT.  Everything two cells away or further,
// the stellar beam and the tau-blend handover stay in the EXPLICIT residual R_i, so the
// fixed point of the outer iteration is the exact backward-Euler balance whatever the
// truncated Jacobian gets wrong -- and the radiative exchange and the radiative
// diffusion are then applied by ONE conservative column update instead of two split ones.
//
// With the switch on the two-stream applies NO de of its own: it stores R_i, the
// three-point Jacobian and dB_i/dT_i on the Conduction object and calls the tridiagonal
// itself, once per outer pass (rt_outer_iter), so Hydro/MHD::ImplicitConduction is a
// no-op.  0 = the old per-cell relaxation, bitwise.
inline int rt_implicit_column = 0;
// convergence tolerance and pass cap of that Newton, used when rt_outer_iter is not set
inline Real rt_impl_tol = 1.0e-6;
inline int rt_impl_maxit = 5;
// problem/rt_col3_ex_iter: under rt_implicit_column = 3, re-form the tau-blend handover
// src_ex from the COLUMN'S OWN converged face flux instead of the entry sweep's, so that
// the applied source is the divergence of one field, div[(1-w)F_3], and telescopes.  See
// two_stream_column_implicit.hpp, step 4a'.  Off = bitwise the frozen handover.
inline bool rt_col3_ex_iter = false;
// problem/rt_col3_skip_sweep: under rt_implicit_column = 3 with the tau blend weight w
// = 0 on EVERY face (rad_tau_lo deeper than the whole box, i.e. the rt_bottom_flux
// production configuration), SKIP the explicit entry sweep, which the column no longer
// consumes.  With w = 0 the handover src_ex is 0 on every cell and, under
// rt_col3_ex_iter, the column re-forms it from its own flux anyway: mode 3 then reads
// only kc, Bb, icut, Qb (identically 0 on the grey path) and the frozen top-face
// intensity, all of which the PRE-kernels and rt_c3_top build without the sweep.  The
// sweep's own products -- Fb (the face flux every diagnostic, rad_f2s, rt_rad_force and
// the surface/history dumps read), Src and Em -- are then supplied differently: Fb is
// written by the COLUMN SOLVE from its own converged D/U intensities (which is the field
// mode 3 actually applies, not the entry-state one), while Src and Em stay at zero, so
// every consumer of THOSE two must be off.  box_convection.cpp enforces the whole list.
// Default false = bitwise off.
inline bool rt_col3_skip_sweep = false;
// ---- mode 3 CONVERGENCE, the switches of the acceleration pass ----------------------
// problem/rt_impl_exjac: carry the rt_col3_ex_iter handover in the JACOBIAN instead of
// lagging it.  With ex_iter on the applied source is A_i = Src_i + div[w F_3]_i + Q_i,
// and div[w F_3] is LINEAR in the very unknowns the block system already carries, so
// there is no excuse for lagging it: off, the missing term makes the whole Newton a
// PICARD iteration wherever w is neither 0 nor 1, and the measured convergence is linear
// at ~0.4 per pass (13-16 passes to 1e-12, bench/bstar_fecz/m3tol).  On, the block rows
// are the exact derivative.  Off = the old lagged assembly, bitwise.
//
// MEASURED, bench/bstar_fecz/m3acc/sweep, the relaxed 15 Msun B-star column at cycle 0,
// passes taken to reach a given residual (norm 0, no step-size stop, maxit 30):
//     tol     1e-2   1e-4   1e-6   1e-8   1e-10  1e-12
//   exjac 0     2      4      7     10     13     16      LINEAR, rate 0.215 per pass
//   exjac 1     2      3      3      3      4      4      QUADRATIC: 1.6e-4, 4.3e-9,
//                                                         5.4e-14, i.e. round-off at 4
// Both reach the SAME fixed point (max|du/u| 6.350328e-02 either way at tol 1e-14), so
// this is purely the path.  Cost of the whole run, 2000 cycles on 1 MI300A:
//   old default (tol 1e-6, maxit 6, never converged, exits on maxit)   1.53e5 zc/s
//   old converged (tol 1e-12, maxit 20, 13-16 passes)                  8.36e4 zc/s
//   new default (exjac, tol 1e-8, maxit 8, 3 passes, resid 5e-10)      2.15e5 zc/s
// against 3.96e5 zc/s for the standard scheme (rt_implicit_column = 0), so the CONVERGED
// mode 3 now costs 1.84x the standard scheme where the unconverged one cost 2.58x.
inline bool rt_impl_exjac = true;
// problem/rt_impl_norm: the residual norm the stopping rule measures.  0 = the old
// max_i |R_i|/e_i, which the optically thin TOP of a stellar column dominates simply
// because its e is 5e-6 of the base's -- a fixed relative tolerance there demands an
// absolute precision nothing else in the scheme has.  1 = max_i |R_i|/(e_i + eps e_max),
// eps = problem/rt_impl_norm_eps, so a cell holding a negligible share of the column's
// energy is measured against that share and not against itself.
inline int rt_impl_norm = 1;
inline Real rt_impl_norm_eps = 1.0e-3;
// problem/rt_impl_dstop: also stop when the Newton STEP is small, max_i |db_i/b_i| < tol.
// A converged Newton takes a negligible step; insisting on the residual as well costs a
// whole extra pass whose only effect is to confirm it.
inline bool rt_impl_dstop = true;
// problem/rt_impl_rescheck: form the energy-row residual BEFORE the block factorisation
// so that the pass which only confirms convergence never factorises.  Bitwise identical
// to the standard path (same test, same value, same iterate); see the note on RTCol3.
inline bool rt_impl_rescheck = false;
// problem/rt_impl_ablate, problem/rt_impl_fixit: TIMING INSTRUMENTATION ONLY, serial
// (thomas) path only.  They do not produce a correct solve -- see RTCol3.
inline int rt_impl_ablate = 0;
inline bool rt_impl_fixit = false;
// problem/rt_impl_cvfreeze: freeze de/dT (the only nonlinear entry of the block) after
// this many passes, 0 = never.  The radiative part of the diagonal is exact and frozen
// already, so this only quasi-Newtons the thin cells.
inline int rt_impl_cvfreeze = 0;
// problem/rt_impl_reuse: REUSE THE BLOCK FACTORISATION across Newton passes.  The
// Jacobian of the column system is constant in the iterate except for the single entry
// de/db (the gas heat capacity), and the right-hand side is one number per cell -- the
// energy-row residual -- so a pass that keeps pass 1's factors costs a residual
// evaluation and two small matrix-vector products instead of the whole assembly,
// factorisation and reduced elimination.  1 = reuse with a contraction check (the pass
// refactorises when the residual failed to fall by a factor 0.3), 2 = reuse always
// (diagnostic).  0 (the default) is the old code, bitwise.  The converged answer is the
// same root of the same residual to the same tolerance; only the iterates between the
// first and the last differ, at O(the frozen de/db), like rt_impl_cvfreeze.  The
// PARTITIONED path (rt_impl_solver = pcr) only; thomas is refused.
inline int rt_impl_reuse = 0;
// problem/rt_impl_reuse_rho: the contraction a reuse pass must show for the NEXT pass to
// keep the factorisation.  The residual of a full Newton pass falls much faster than of
// a frozen-Jacobian one, so a loose threshold buys cheap passes at the price of more of
// them; tighten it to refactorise sooner.
inline Real rt_impl_reuse_rho = 0.3;
// problem/rt_impl_tau_min: THE TWO-LEVEL SPLIT.  Only a cell whose OWN Rosseland optical
// depth kappa rho dr reaches this goes into the tridiagonal.  Linearising the emission
// about the current state gives dT ~ (T/4)(A/E), which diverges as the cell's own
// emission E -> 0: measured on the He-star box the first merged step emptied 33 cells and
// collapsed dt by three decades, which is the divergence the apply block's own note
// ("WHY NOT LINEARIZE ... DIVERGES as E -> 0") warns about.  In a cell with dtau >~ 1 the
// field is within a factor of its own B, A/E = O(1), and the linearisation is excellent
// -- and those are exactly the cells the (1 - e^-x)/x damping pumps.  Thin cells keep
// the nonlinear per-cell relaxation, which relaxes toward the TRUE fixed point and is
// bounded.
inline Real rt_impl_tau_min = 1.0;
// problem/rt_impl_dtmax: cap on |dT|/T per tridiagonal pass, the belt for the thick rows
inline Real rt_impl_dtmax = 0.25;
// problem/rt_impl_tau_blend: 1 = the hard switch above.  b > 1 makes the thick weight
// rise LINEARLY IN log dtau from tau_min/b to tau_min*b, and the cell then does both:
// the tridiagonal carries w R and w J, the per-cell relaxation applies (1-w) de, and a
// neighbour's coupling to it is w-weighted with the remaining (1-w) handed over as the
// known dB of rt_col_dtex.  At w = 0 and w = 1 this is bitwise the two pure branches.
inline Real rt_impl_tau_blend = 1.0;
// ---- rt_implicit_column = 3: THE EXACT IMPLICIT COLUMN SOLVE ------------------------
// The intensities become unknowns alongside the gas energy, so the column system is
// exactly block-tridiagonal (5x5) and a block Thomas solves it with no Jacobi/
// Gauss-Seidel iteration at all -- which is what modes 1 and 2 could not do.  See
// utils/two_stream_column_implicit.hpp for the rows, the M-matrix argument and what is
// frozen.  Mode 3 REPLACES the two-stream's own semi-implicit apply and nothing else:
// the tau-blend handover, the radial conduction tridiagonal (rad_implicit_x1) and the
// transverse operator are untouched, and rt_col_active stays false so the conduction
// wrapper task still runs.  Grey only; fatal under correlated-k or picket fence.
// problem/rt_impl_solver: HOW the block-tridiagonal column system is solved.
//   thomas (0, the default)  one thread per column, the serial block Thomas of
//                            two_stream_column_implicit.hpp -- the reference bit pattern
//   pcr    (1)               one TEAM per column: the column is cut into
//                            problem/rt_impl_nseg segments, each thread eliminates its
//                            own segment with the incoming unknown carried symbolically,
//                            and the resulting 5x5 block system over the segment
//                            BOUNDARIES is solved by one thread and back-substituted in
//                            parallel.  Same system, same Newton, same clamps; the two
//                            agree to round-off.  See two_stream_column_partition.hpp.
inline int rt_impl_solver = 0;
// problem/rt_impl_mixed: MIXED PRECISION in the partitioned column solve's block algebra.
//   0 (default)  everything double -- bit for bit the code before the switch existed
//   1            the Newton CORRECTION is formed in SINGLE precision: the 5x5 inverses
//                (RTCol3Inv5X), the segment forward eliminations, the p/Q/R recurrence
//                and the back-substitutions.  The residual, the convergence test, the
//                layer/Planck/kappa coefficients, the reduced system and the u0(IEN)
//                update stay double, so the outer Newton is an iterative refinement of
//                the float factor and converges to the same rt_impl_tol.
//   2            ... and the STORED factors (G, H, d, and the reuse pair v, M) live in
//                their own float workspace, which halves their memory traffic.
// A block whose float Gauss-Jordan loses its pivots is redone in double and counted in
// stat(24) (nmixfb under rt_outer_verbose).  Partitioned path only.
inline int rt_impl_mixed = 0;
inline int rt_impl_nseg = 64;
// problem/rt_impl_redpar: on the pcr path, the REDUCED block-tridiagonal system over the
// segment boundaries (nsg 5x5 rows, non-periodic) is by default eliminated serially by
// one lane while the other nsg-1 idle -- ~29 % of a Newton pass at nseg = 32.  With this
// switch it is solved by PARALLEL CYCLIC REDUCTION instead: lane s owns block row s and
// log2(nsg) rounds of neighbour eliminations leave y_s = B_s^-1 r_s on every lane at
// once, with no reduced back-substitution.  Same system, different elimination order, so
// the answer is equal to ROUND-OFF, not bitwise.  Falls back to the serial solve unless
// 2 <= nsg <= 64 and nsg is a power of two.  Costs the larger per-segment workspace
// (RTCOL3_NRDP), so the array is only enlarged when the switch is on.
inline bool rt_impl_redpar = false;
// problem/rt_col3_hybrid_tau: THE HYBRID COLUMN.  Below this COLUMN optical depth the
// cells are solved as DIFFUSION -- one unknown per cell, F = (4 pi/3) dB/dtau across the
// faces, which is the exact deep limit of the same two-stream -- instead of as the full
// five-unknown two-stream, all inside the SAME Newton system and the same block
// elimination (1x1 blocks deep, 5x5 thin, a 1x5/5x1 pair at the interface).  The two
// segments are coupled by flux continuity: the thin segment's lower boundary is the deep
// limit U_q = b_f + mu_q dB/dtau with b_f and dB/dtau carried implicitly by the deep
// unknowns, and the deep segment's top face loses exactly the two-stream's own net flux
// at that face.  0 = off, bitwise the whole-column solve.  See
// utils/two_stream_column_implicit.hpp, RTCol3Hyb.  Supported on BOTH solvers: with
// rt_impl_solver = pcr the deep segments run a scalar partitioned Thomas and the reduced
// system has mixed block sizes.  See two_stream_column_partition.hpp.
inline Real rt_col3_hybrid_tau = 0.0;
// problem/rt_col3_split_deep: BALANCE the partitioned path's segments by WORK.  The
// hybrid makes a deep cell a scalar row, but the team's segments are equal in CELLS and
// every phase ends on a team barrier, so the lanes that hold the thin cells still carry
// nc/nseg full 5x5 rows and the pass costs what it did with the hybrid off -- the 1.09x.
// With this switch a thin cell counts rt_col3_split_w deep cells when the boundaries are
// laid out, the thin cells spread over all nseg lanes, and the critical path falls to
// about n_thin/nseg 5x5 rows.  Partition only: the answer moves by round-off.
// MEASURED SLOWER, and kept off: see the note on RTCol3::split_deep.  The team's lanes
// are one wavefront, so a cell loop costs deep-plus-thin and scales with the most cells
// any lane holds; the equal partition already minimises that.  tau_hyb 30 on the B star:
// 223.7 ms/call plain, 258.4 balanced, 293.1 with the hybrid off.
inline bool rt_col3_split_deep = false;
inline int rt_col3_split_w = 8;   // problem/rt_col3_split_w, thin cost / deep cost
// problem/rt_impl_warm: WARM-START the mode-3 Newton from the previous call's converged
// Planck function instead of from the entry state's.  0 = off (bitwise the old code),
// 1 = the previous converged b per cell, 2 = linear extrapolation in time from the last
// two.  Only the ITERATE moves; the residual, the Jacobian, the tolerance and the entry
// state are untouched, so the converged root is the same one and the answer is equal to
// round-off rather than bitwise.  Costs one (warm = 1) or two (warm = 2) extra 4D arrays.
inline int rt_impl_warm = 0;
inline DvceArray5D<Real> *rt_c3wk_ptr = nullptr;
inline DvceArray5D<float> *rt_c3wkf_ptr = nullptr;
inline DvceArray4D<Real> *rt_c3rd_ptr = nullptr;
inline DvceArray4D<Real> *rt_c3top_ptr = nullptr;
inline DvceArray1D<Real> *rt_c3stat_ptr = nullptr;
inline DvceArray4D<Real> *rt_c3bp_ptr = nullptr;   // the warm-start history, level n
inline DvceArray4D<Real> *rt_c3bp2_ptr = nullptr;  // the warm-start history, level n-1
inline Real rt_c3_bdt_prev = 0.0;                  // the previous call's bdt
inline int rt_c3_lines = 0;
// problem/ck_int_at_cut: deliver the planet's internal flux sigma T_int^4 as an extra
// upward source at the correlated-k cut (the historical behaviour, true). Set false when
// the layers below the cut carry it themselves -- <mhd|hydro>/isotropic_conduction =
// radiative with rad_flux_inner at the bottom wall -- so the cut's upward intensity is
// just the thermalised Planck function and nothing is counted twice.
inline bool rt_int_at_cut = true;
// problem/rt_cut_bc_legacy: the UPWARD intensity the grey sweep starts with at the cut.
//
// The sweep used to start with the isotropic thermalised value I_up(mu) = B(icut) for
// every mu.  That is the zeroth term of the deep expansion and it drops the first one,
// which is the ONLY term that carries a flux: with I_up = B and the down-sweep already
// carrying its own gradient, sum_q w_q (I_up - I_down) comes out at EXACTLY half the
// diffusion flux the same column supports.  MEASURED on the He-star 1-D column
// (bench/hestar_fecz/instab1d): F_2s(cut)/F_diff(cut) = 2.0002, healing over ~3 cells
// above the cut, and the residual of that mismatch is what the tau blend then deposits
// in the ramp.
//
// The first term is the standard diffusion limit,
//     I_up(mu) = B + mu dB/dtau,
// with tau increasing DOWNWARD (into the star, i.e. towards smaller i), so dB/dtau is
// formed here from the two cell-centre Planck functions straddling the cut face and the
// Rosseland optical depth between those two CENTRES,
//     dB/dtau = (B_icut - B_{icut+1}) / (0.5 (kappa rho dz)_icut
//                                        + 0.5 (kappa rho dz)_{icut+1}),
// which is positive on a star (deeper is hotter) and so raises I_up above B.  With the
// two-point Gauss-Legendre quadrature (rt_nquad = 2) this reproduces F = (4 pi/3) dB/dtau
// exactly; with the hemispheric mean (rt_nquad = 1) it gives 2 pi mu_H dB/dtau = 0.904 of
// it, which is the quadrature's own error and not this boundary's.
//
// Default false, i.e. THE FIX IS ON: the old behaviour was a bug.  Set true only to
// reproduce a pre-fix run bitwise.  Grey split path only -- the correlated-k sweep keeps
// its own cut boundary, which is never used under a tau blend.
inline bool rt_cut_bc_legacy = false;
// problem/rt_bottom_flux: THE IMPOSED FLUX ON THE TWO-STREAM'S OWN LOWER BOUNDARY.
//
// With the tau blend switched off (rad_tau_lo/hi beyond the bottom of the box, so the
// blend weight is 0 on every face) the cut falls on the inner x1 wall and the two-stream
// owns the WHOLE column: the radiative-diffusion operator has nothing left to carry, and
// the handover that splits the two is gone.  The energy then has to enter through the
// two-stream's own lower boundary instead of through the conduction wall face.
//
// This is the code-unit flux that boundary carries.  The cut's upward intensity is the
// deep limit I_up(mu) = B + mu dB/dtau with the gradient set BY THE FLUX rather than read
// off the two cells above it,
//     dB/dtau = 3 F_bot/(4 pi),
// which with the two-point Gauss-Legendre quadrature returns sum_q w_q (I_up - I_down) =
// F_bot exactly (rt_nquad = 1 returns 0.904 of it, the hemispheric mean's own error).
// 0 = off, the local gradient of rt_cut_bc_legacy = false.  box_convection.cpp sets it
// from <hydro>/rad_flux_inner under <problem>/rt_bottom_flux, and zeroes rad_flux_inner
// at the same time so the flux is not also injected at the conduction wall face.
inline Real rt_bot_flux = 0.0;
// problem/rt_layer_legacy: the STAGGERED layer the grey and correlated-k sweeps used.
//
// Every layer was integrated over a WHOLE cell in optical depth, dtau = kappa rho dz of
// cell i, but with its source running linearly in tau from that cell's CENTRE Planck
// function to the NEIGHBOUR CELL'S CENTRE one.  The interval the optical depth measures
// and the interval the source runs over are then offset by half a cell, and the two
// endpoints of the source are two cells apart while the layer is one cell thick.  The
// error is O(dtau^2) and smooth, but it is not small where it matters: MEASURED against
// an accurate formal solution on the He-star column, the emergent flux came out ~4 % low
// at dtau/cell of 1-3 and wandered between 0.73 and 1.06 of the true flux at dtau/cell
// 3-6, which is exactly the handover ramp of the He and B-star boxes.
//
// THE FIX.  The layers now run between cell CENTRES: layer (i, i+1) has
//     dtau = 0.5 (kappa rho dz)_i + 0.5 (kappa rho dz)_{i+1},
// each half with its OWN opacity, and a source linear in tau between B_i and B_{i+1} --
// which is exact for that layer, since the source really is (piecewise) linear there.
// The layer is integrated as its two halves, split at the face, with the face source
// taken at its own tau within the layer; composing the two halves reproduces the whole
// layer exactly, so this costs nothing in accuracy and hands back the FACE intensity as a
// by-product.  That face intensity is what Fb_g reports, so the fluxes the tau blend and
// the heating read are the formal solution AT the face, not an interpolation of it.  The
// column is closed by two HALF layers: the upper half of cell ie, from the top face to
// its centre, with the source held at B(ie); and the lower half of cell icut, between its
// centre and the cut face, with B continued at the deep-limit gradient dB/dtau (zero
// under rt_cut_bc_legacy, so that flag still recovers its own old behaviour).
//
// Each half layer lies entirely inside ONE cell, so the absorbed and emitted amounts that
// build Src_g are attributed unambiguously: a cell's emission is simply split between its
// two halves.  Em_g, the relaxation rate, becomes the cell's own 4 sigma kappa rho T^4
// instead of the centre-to-centre average of two Planck functions.
//
// Default false, i.e. THE FIX IS ON.  Set true to reproduce a pre-fix run bitwise.
inline bool rt_layer_legacy = false;
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

//----------------------------------------------------------------------------------------
//! \fn Real BFaceW
//! \brief d(BFace)/d(b_far): BFace is a CONVEX COMBINATION of b_far and b_own at frozen
//! opacity, so one weight describes it completely and the own-side weight is 1 - this.
//! Used only by the merged column solve (rt_implicit_column) to linearise the sweep.
KOKKOS_INLINE_FUNCTION
Real BFaceW(const Real k_own, const Real k_far, const bool on) {
  if (!on) return 1.0;
  const Real kt = RT_BFACE_R*k_own;
  if (k_far >= kt) return 1.0;
  if (!(k_own > 0.0)) return 1.0;
  if (!(k_far > 0.0)) return 0.0;
  const Real w = k_far/kt;
  return w + (1.0 - w)*k_far/(k_own + k_far);
}

//----------------------------------------------------------------------------------------
//! \fn void RTLayer
//! \brief one exact short-characteristic step through a layer (see rt_layer_legacy).
//!
//! The layer has optical thickness dtau along the ray (mu the direction cosine) and a
//! source linear in tau from s_in, where the ray ENTERS, to s_out, where it LEAVES.  The
//! exponential coefficients are the ones the sweeps have always used, with their series
//! forms below x = 1e-3; what changes with the new layers is only which interval they are
//! handed.  absorbed and emitted come back separately because the per-cell source needs
//! them apart, and because their difference is exactly the change in the intensity, which
//! is what keeps Src_g and the divergence of the face fluxes the same number.
KOKKOS_INLINE_FUNCTION
void RTLayer(const Real dtau, const Real mu, const Real s_in, const Real s_out,
             Real &intens, Real &absorbed, Real &emitted) {
  const Real x = dtau/mu;
  const Real e0 = -expm1(-x);
  const Real c_in  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
  const Real c_out = (x > 1.0e-3) ? (1.0 - e0/x)      : (x/2.0 - SQR(x)/6.0);
  absorbed = e0*intens;
  emitted  = c_in*s_in + c_out*s_out;
  intens   = (1.0 - e0)*intens + emitted;
}

//! \fn void RTLayerCoef
//! \brief the three exponential coefficients RTLayer forms internally, without the
//! transport.  Used only by rt_implicit_column, to linearise a layer in its endpoints.
KOKKOS_INLINE_FUNCTION
void RTLayerCoef(const Real dtau, const Real mu, Real &e0, Real &c_in, Real &c_out) {
  const Real x = dtau/mu;
  e0 = -expm1(-x);
  c_in  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
  c_out = (x > 1.0e-3) ? (1.0 - e0/x)      : (x/2.0 - SQR(x)/6.0);
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

inline void picket_fence_two_stream_RT_pass(Mesh *pm, Real bdt, const int oit,
                                           const int nit);

//----------------------------------------------------------------------------------------
//! \fn picket_fence_two_stream_RT
//! \brief the two-stream source step.  One pass unless problem/rt_outer_iter > 1, in
//! which case the whole step (sweep + apply) is repeated as a fixed-point iteration for
//! the implicit balance; see rt_outer_iter.

// ---- problem/work_hist: THE PER-OPERATOR PROBE (box_convection) ---------------------
// A diagnostic hook, null by default.  The only call site in this header is at the end
// of the mode-3 column solve, with tag 2: at that point the column has applied all of
// its energy to u0 and the ONLY thing left in the call is the radiative momentum force
// (rt_apply's energy branch is skipped under mode 3), so the interval that closes there
// is exactly "column heating" and the interval that closes when the call returns is
// exactly "radiative force".  Left null nothing is called and nothing changes.
inline void (*rt_probe)(const int tag) = nullptr;

// ---- problem/rt_kappa_frozen: THE MODE SEES NO OPACITY PERTURBATION -----------------
// With this on, the grey opacity kc_g is replaced, after it has been built from the
// local (rho,T) as usual, by its HORIZONTAL (x2,x3) MEAN in each x1 row.  So the mean
// opacity profile still follows the mean state stage by stage, while delta kappa (the
// part that carries the horizontal mode) is identically zero -- which kills the
// kappa-mechanism channel in the column solve AND in the radiative force, both of which
// read kc_g and nothing else for the opacity.  (It is the mean of kappa, not kappa of
// the mean state; the two differ only at second order in the perturbation, which is
// below the linear-mode test's resolution.)  The horizontal mean is GLOBAL: every rank
// reduces its own plane sums and the result is broadcast back, because a per-rank mean
// of a cos(k x) pattern is not zero.  The transverse ADI operator builds its own
// conductances in Conduction and is NOT covered.  Grey path only.
inline bool rt_kappa_frozen = false;
inline DvceArray1D<Real> *rt_kfrz_d_ptr = nullptr;
inline HostArray1D<Real> *rt_kfrz_h_ptr = nullptr;

inline void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
  const int nit = (rt_outer_iter > 1) ? rt_outer_iter : 1;
  if (nit > 1) {
    static bool checked = false;
    if (!checked) {
      checked = true;
      if (!(rt_grey && rt_split) || !rt_use_cons || rt_explicit || !rt_semi_implicit) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_outer_iter > 1 needs "
                  << "the GREY SPLIT sweep (rt_grey + rt_split), problem/rt_use_cons "
                  << "(the passes update u0, and w0 is stale until the next ConToPrim) "
                  << "and the semi-implicit apply (rt_semi_implicit, !rt_explicit). "
                  << "Got rt_grey=" << rt_grey << " rt_split=" << rt_split
                  << " rt_use_cons=" << rt_use_cons << " rt_explicit=" << rt_explicit
                  << " rt_semi_implicit=" << rt_semi_implicit << "." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  // ---- problem/rt_implicit_column: the merged tridiagonal Newton --------------------
  // Requirements, checked once.  The column solve linearises the GREY SPLIT sweep, reads
  // the running state out of u0 (rt_use_cons) between passes, and replaces -- not
  // supplements -- the semi-implicit per-cell relaxation, which therefore has to be the
  // step it is replacing.
  Conduction *pc = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->pcond
                                                   : pm->pmb_pack->phydro->pcond;
  if (rt_implicit_column == 3) {
    // ---- mode 3, the exact block-tridiagonal column solve ---------------------------
    // It needs the grey split sweep (it IS that sweep, differentiated), the conserved
    // state (it applies de to u0 and must read the post-RK state), the direct per-cell
    // source (its rows ARE the per-half-layer absorbed-minus-emitted balance), and the
    // centre-to-centre layers.  It does NOT need rad_implicit_x1: the radial conduction
    // stays a separate implicit tridiagonal, applied by its own task.
    static bool c3checked = false;
    if (!c3checked) {
      c3checked = true;
      if (!(rt_grey && rt_split) || rt_ck || !rt_use_cons || rt_explicit ||
          !rt_semi_implicit || !rt_src_direct || rt_layer_legacy || rt_top_re ||
          rt_outer_iter > 1) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_implicit_column = 3 "
                  << "needs the GREY SPLIT sweep (rt_grey + rt_split, and NOT "
                  << "correlated-k or picket fence), problem/rt_use_cons, "
                  << "rt_src_direct, the centre-to-centre layers (rt_layer_legacy = "
                  << "false), rt_top_re = false, rt_outer_iter = 1 and the semi-implicit "
                  << "apply it replaces (rt_semi_implicit, !rt_explicit).  Got rt_grey="
                  << rt_grey << " rt_split=" << rt_split << " rt_ck=" << rt_ck
                  << " rt_use_cons=" << rt_use_cons << " rt_explicit=" << rt_explicit
                  << " rt_semi_implicit=" << rt_semi_implicit
                  << " rt_src_direct=" << rt_src_direct
                  << " rt_layer_legacy=" << rt_layer_legacy
                  << " rt_top_re=" << rt_top_re
                  << " rt_outer_iter=" << rt_outer_iter << "." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  } else if (rt_implicit_column > 0) {
    static bool cchecked = false;
    if (!cchecked) {
      cchecked = true;
      if (!(rt_grey && rt_split) || !rt_use_cons || rt_explicit || !rt_semi_implicit ||
          pc == nullptr || !pc->rad_implicit_x1) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_implicit_column "
                  << "needs the GREY SPLIT sweep (rt_grey + rt_split), "
                  << "problem/rt_use_cons, the semi-implicit apply (rt_semi_implicit, "
                  << "!rt_explicit) and <hydro>/ or <mhd>/rad_implicit_x1 = true.  Got "
                  << "rt_grey=" << rt_grey << " rt_split=" << rt_split
                  << " rt_use_cons=" << rt_use_cons << " rt_explicit=" << rt_explicit
                  << " rt_semi_implicit=" << rt_semi_implicit
                  << " rad_implicit_x1="
                  << ((pc != nullptr) ? pc->rad_implicit_x1 : false) << "." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  for (int oit=0; oit<nit; ++oit) {
    picket_fence_two_stream_RT_pass(pm, bdt, oit, nit);
    // mode 3 solved and applied the column inside the pass; the radial conduction stays
    // a separate operator, applied by the ImplicitConduction task as usual.
    if (rt_implicit_column > 0 && rt_implicit_column != 3 && pc != nullptr) {
      // the sweep wrote R, the Jacobian and dB/dT for the CURRENT state; solve the
      // column for dT and let the tridiagonal apply the energy, radiative exchange and
      // radiative diffusion together
      MeshBlockPack *pp = pm->pmb_pack;
      if (pp->pmhd != nullptr) {
        pc->ImplicitRadialUpdate(pp->pmhd->u0, pp->pmhd->peos->eos_data, bdt, true, oit);
      } else {
        pc->ImplicitRadialUpdate(pp->phydro->u0, pp->phydro->peos->eos_data, bdt, true,
                                 oit);
      }
    }
  }
}

inline void picket_fence_two_stream_RT_pass(Mesh *pm, Real bdt, const int oit,
                                            const int nit) {
  const bool outer_on = (nit > 1);
  const bool outer_last = (oit == nit - 1);
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

    // --- PLANE-PARALLEL geometry (problem/rt_plane_parallel).  See the flag's note: the
    // three radial geometry Views above are placeholders on a Cartesian mesh, so in this
    // mode every read of them goes through X1V / X1F / DX1 instead, which rebuild the
    // same quantities from the MeshBlock's RegionSize.  With the flag off each helper
    // returns the View element it replaced, so the radial path is unchanged bit for bit.
    const bool pp_ = rt_plane_parallel;
    if (pp_ && !(rt_grey && rt_split)) {
      std::cout << "### FATAL ERROR in two_stream_rt: rt_plane_parallel is implemented "
                << "for the GREY split sweep only (rt_grey = true, which implies "
                << "rt_split); the picket-fence and correlated-k kernels carry radial "
                << "stellar-beam geometry." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pp_ && (use_cubed_sphere_ || use_spherical_polar)) {
      std::cout << "### FATAL ERROR in two_stream_rt: rt_plane_parallel is for a "
                << "Cartesian mesh, but this mesh is spherical-polar or cubed-sphere."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const int is_pp = is;
    const int nx1_pp = indcs.nx1;
    auto X1V = [=] (const int m, const int i) {
      return pp_ ? CellCenterX(i-is_pp, nx1_pp, size.d_view(m).x1min,
                               size.d_view(m).x1max)
                 : x1v_(m,i);
    };
    auto X1F = [=] (const int m, const int i) {
      return pp_ ? LeftEdgeX(i-is_pp, nx1_pp, size.d_view(m).x1min,
                             size.d_view(m).x1max)
                 : x1f_(m,i);
    };
    auto dx1_ = pmbp->pcoord->dx1;
    auto DX1 = [=] (const int m, const int k, const int j, const int i) {
      return pp_ ? size.d_view(m).dx1 : dx1_(m,k,j,i);
    };

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
    // MHD: u0(IEN) also holds the magnetic energy, so the extraction has to subtract it
    // or every temperature this solver forms is the temperature of e + ME.  bcc0 is the
    // cell-centred form of the CURRENT b0 (see MagEnergyCC): the RT apply and its
    // precompute run from MHD::MHDSrcTerms, inside the stage but BEFORE MHD::CT, so b0
    // has not moved since the ConToPrim that filled bcc0 and re-averaging the faces here
    // would return the same numbers.  On the cubed sphere bcc0 is already in the
    // orthonormal frame, so no metric enters.
    const bool mhd_uc_ = (pmbp->pmhd != nullptr);
    DvceArray5D<Real> bcc_uc_("rt_bcc_dummy", 1, 1, 1, 1, 1);
    if (pmbp->pmhd != nullptr) {
      bcc_uc_ = pmbp->pmhd->bcc0;
    }
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
                                      etg_uc_, etg_uc_ ? phicc_uc_(m,k,j,i) : 0.0,
                                      mhd_uc_ ? MagEnergyCC(bcc_uc_,m,k,j,i) : 0.0);
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

    // ================================================================================
    // SPHERICAL DILUTION: WHY THE DEPOSIT *AND* THE PROPAGATION CARRY THE AREA
    //
    // The sweep integrates, per ray direction mu, the plane-parallel transfer equation
    //     mu dI/dr = -kappa rho (I - B),
    // and used to deposit -(F_top - F_bot)/dx1 with F = 2 pi w mu I.  On a RADIAL mesh
    // that is wrong twice over, and the two errors have to be fixed together:
    //   (a) the deposit is a plane-parallel divergence.  What the shell actually loses
    //       is (A_top F_top - A_bot F_bot)/V, and over r_out/r_in = 2.5 the area ratio
    //       is 6.25 -- the emergent luminosity of a global star was wrong by that much;
    //   (b) the PROPAGATION does not dilute.  A transparent shell returns I_out = I_in,
    //       i.e. F_out = F_in, so L = A F grows like r^2 -- and had only (a) been fixed
    //       that same transparent shell would be handed a spurious cooling
    //       -(A_top - A_bot) F/V = -2 F/r.  NEITHER FORM ALONE IS CONSERVATIVE.
    //
    // The exact statement is the zeroth moment, div F = kappa rho (4 pi B - c E), i.e.
    // (1/r^2) d(r^2 F)/dr on the left.  Multiply the transfer equation by r^2 and write
    // J = r^2 I (equivalently J = A I; the column's solid angle cancels):
    //     mu dJ/dr = -kappa rho (J - A B),
    // which is the SAME equation the layer solve already integrates exactly, with the
    // source function A(r) B(r) in place of B(r).  So the sweep keeps its form and
    // carries the AREA-WEIGHTED intensity J instead of I:
    //   * every layer source endpoint is multiplied by the area AT THAT ENDPOINT -- a
    //     face endpoint by area.x1f (AFC), a cell-centre endpoint by the cell's mean
    //     area V/dx1 (ACC);
    //   * the boundary data are scaled the same way (the top incoming intensity by
    //     A(ie+1), the internal flux entering the bottom face by A(is));
    //   * the face flux handed to the rest of the code stays a flux PER UNIT AREA,
    //     recovered as F = 2 pi w mu J / A;
    //   * the deposit is (J_in - J_out), summed over rays, divided by the cell VOLUME.
    //     It telescopes to exactly -(A_top F_top - A_bot F_bot)/V, so the column
    //     conserves energy cell by cell and L = A F is constant wherever the local
    //     balance vanishes -- which is the 1-D radiative-equilibrium gate.
    // The ray-curvature term of the true spherical transfer equation, (1-mu^2)/r dI/dmu,
    // is NOT represented: this is the standard radial-ray approximation, exact in the
    // two limits that matter here (diffusion, and a radially streaming flux) and
    // conservative everywhere.
    //
    // PLANE-PARALLEL INERTNESS.  Under problem/rt_plane_parallel the helpers return
    // A = 1 exactly and V = the same dx1 the expression used before, so every new
    // factor is a multiplication or a division by 1.0 -- exact in IEEE arithmetic --
    // and the Cartesian box is bitwise what it was.  The helpers are also the only
    // reads of the Coordinates area/volume Views, which are 1x1x1x1 placeholders on a
    // Cartesian mesh (coordinates.cpp:63-94) exactly like x1v/dx1, so the ternary is
    // what keeps them unread there.
    //
    // WHAT IS NOT CONVERTED: the monolithic (rt_split = false) sweep, the correlated-k
    // kernel and the mode-1/2 nearest-neighbour Jacobian still carry the plane-parallel
    // divergence.  They serve the box and the r_out/r_in < 1.3 hot-Jupiter/red-giant
    // domains, where the error is at the per-cent level.  The GREY SPLIT sweep
    // (rt_grey + rt_split), the apply kernel and the mode-3 column solve -- the path a
    // global star runs on -- are area/volume correct.
    // ================================================================================
    // the face area, the cell's MEAN area V/dx1 (the weight of a cell-centre source
    // endpoint), and the cell volume.  THE TWO VOLUME HELPERS DIFFER ONLY IN THEIR
    // PLANE-PARALLEL BRANCH, and each kernel must use the one matching the radial width
    // it already spells: VLS for the picket-fence chain kernel, which indexes the raw
    // dx1 View, and VLA wherever the code goes through DX1 (the grey sweep and the
    // apply kernel).  That is what makes the substitution bitwise inert on the box --
    // getting it the wrong way round reads the 1x1x1x1 dx1 placeholder there.
    auto AFC = [=] (const int m, const int k, const int j, const int i) {
      return pp_ ? 1.0 : area1(m,k,j,i);
    };
    auto ACC = [=] (const int m, const int k, const int j, const int i) {
      return pp_ ? 1.0 : volume(m,k,j,i)/dx1(m,k,j,i);
    };
    auto VLS = [=] (const int m, const int k, const int j, const int i) {
      return pp_ ? dx1(m,k,j,i) : volume(m,k,j,i);
    };
    auto VLA = [=] (const int m, const int k, const int j, const int i) {
      return pp_ ? DX1(m,k,j,i) : volume(m,k,j,i);
    };

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
    const bool cut_legacy = rt_cut_bc_legacy;   // see rt_cut_bc_legacy
    const Real bot_flux = rt_bot_flux;          // see rt_bottom_flux
    const bool layer_legacy = rt_layer_legacy;  // see rt_layer_legacy
    const bool top_re = rt_top_re;
    const bool top_vac = rt_top_vacuum;
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
      // ---- problem/rt_implicit_column: the merged tridiagonal Newton -----------------
      // The requirements are checked once, in the wrapper.  The arrays are allocated on
      // the Conduction object the first time through; with the switch off implcol_ is
      // false everywhere and jac_g is a 1-element dummy that is captured and never read.
      const bool implcol_ = (rt_implicit_column > 0) && (rt_implicit_column != 3) &&
                            grey_on && rt_split;
      // rt_implicit_column = 3: the exact block-tridiagonal column solve.  It runs as a
      // separate kernel right after the sweep and owns the energy update, so the apply
      // block below skips de entirely (skip_de) and keeps only its diagnostics and the
      // radiative momentum source.
      const bool mode3_ = (rt_implicit_column == 3) && grey_on && rt_split;
      // problem/rt_col3_skip_sweep: the entry sweep is dead weight under mode 3 when the
      // blend weight is 0 everywhere.  The REQUIREMENTS are checked in the problem
      // generator; the one thing it cannot check there is that w really is 0 on every
      // face (the weights are built per cycle), so that is checked here, once.
      const bool skipsweep_ = mode3_ && rt_col3_skip_sweep;
      if (skipsweep_) {
        static bool wchecked = false;
        if (!wchecked) {
          wchecked = true;
          Real wmax = 0.0;
          if (taublend) {
            Kokkos::parallel_reduce("rt_c3_wchk",
              Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*(ke-ks+1)*(je-js+1)
                                                      *(ie+2-is)),
              KOKKOS_LAMBDA(const int idx, Real &mx) {
                const int nn1 = ie + 2 - is, nn2 = je - js + 1, nn3 = ke - ks + 1;
                const int i = is + (idx % nn1);
                const int j = js + ((idx/nn1) % nn2);
                const int k = ks + ((idx/(nn1*nn2)) % nn3);
                const int m = idx/(nn1*nn2*nn3);
                const Real w = fabs(w_g(m,k,j,i));
                if (w > mx) mx = w;
              }, Kokkos::Max<Real>(wmax));
          }
          if (wmax > 0.0) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_col3_skip_sweep "
                      << "needs the tau-blend weight to be 0 on EVERY x1 face (the "
                      << "handover the skipped sweep would feed is then identically "
                      << "zero), but max|w| = " << wmax << ". Set <hydro>/rad_tau_lo "
                      << "deeper than the bottom of the box, or drop the switch."
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
      }
      if (implcol_ && pcond_rt != nullptr && !pcond_rt->rt_col_alloc) {
        pcond_rt->EnableRTColumn();
      }
      auto jac_g = (implcol_ && pcond_rt != nullptr)
                 ? pcond_rt->rt_col_jac : DvceArray5D<Real>("rt_jac_dummy",1,1,1,1,1);
      auto res_g = (implcol_ && pcond_rt != nullptr)
                 ? pcond_rt->rt_col_res : DvceArray4D<Real>("rt_res_dummy",1,1,1,1);
      auto dbt_g = (implcol_ && pcond_rt != nullptr)
                 ? pcond_rt->rt_col_dbdt : DvceArray4D<Real>("rt_dbt_dummy",1,1,1,1);
      auto dtx_g = (implcol_ && pcond_rt != nullptr)
                 ? pcond_rt->rt_col_dtex : DvceArray4D<Real>("rt_dtx_dummy",1,1,1,1);
      auto tn_g = (implcol_ && pcond_rt != nullptr)
                 ? pcond_rt->rt_col_tn : DvceArray4D<Real>("rt_tn_dummy",1,1,1,1);
      const Real taumin_ = rt_impl_tau_min;
      const Real taublnd_ = (rt_impl_tau_blend > 1.0) ? rt_impl_tau_blend : 0.0;
      // the one-shot assembly dump, see rtcol_asm below
      const bool rtdbg_ = implcol_ && rt_outer_verbose && (pm->ncycle == 0);
      if (implcol_ && pcond_rt != nullptr) pcond_rt->rt_col_verbose = rtdbg_;
      if (implcol_ && pcond_rt != nullptr) pcond_rt->rt_col_dtmax = rt_impl_dtmax;
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
          // PLANE-PARALLEL: there is no substellar direction, and x2v/x3v are 1x1
          // placeholder Views on a Cartesian mesh, so they must not be read at all.
          // mu0 reaches only EffGravAt's tidal term, which is off here.
          if (pp_) {
            cf_g(m,k,j,3) = 0.0;
            return;
          }
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
            Real kr = (grey_krmax > 0.0 && X1V(m,i) > grey_krmax)
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
          // ---- problem/rt_kappa_frozen: kc_g -> its horizontal mean in each x1 row ---
          // See the note at the top of picket_fence_two_stream_RT.  Only the cells this
          // call actually reads take part (i >= icut, T > 0), so a row that straddles a
          // column cut is still the mean of the cells that are used; the cells below the
          // cut are left alone because nothing reads them.  Off by default = no kernel.
          if (rt_kappa_frozen) {
            const int nr = ie + 2 - is;          // the rows i = is .. ie+1
            if (rt_kfrz_d_ptr == nullptr) {
              rt_kfrz_d_ptr = new DvceArray1D<Real>("rt_kfrz_d", 2*nr);
              rt_kfrz_h_ptr = new HostArray1D<Real>("rt_kfrz_h", 2*nr);
            }
            auto kfd = *rt_kfrz_d_ptr;
            auto kfh = *rt_kfrz_h_ptr;
            const int nkj_f = (nmb1 + 1)*(ke - ks + 1)*(je - js + 1);
            const int nx2_f = je - js + 1, nx3_f = ke - ks + 1;
            Kokkos::TeamPolicy<> pol_f(DevExeSpace(), nr, Kokkos::AUTO);
            Kokkos::parallel_for("rt_kfrz_sum", pol_f,
            KOKKOS_LAMBDA(Kokkos::TeamPolicy<>::member_type tm) {
              const int i = is + tm.league_rank();
              array_sum::GlobalSum sm;
              Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nkj_f),
              [&](const int idx, array_sum::GlobalSum &ls) {
                const int m = idx/(nx3_f*nx2_f);
                const int kj = idx - m*(nx3_f*nx2_f);
                const int k = ks + kj/nx2_f;
                const int j = js + (kj - (kj/nx2_f)*nx2_f);
                const bool ok = (i >= icut_g(m,k,j)) && (T_g(m,k,j,i) > 0.0);
                array_sum::GlobalSum lv;
                for (int n=0; n<NREDUCTION_VARIABLES; ++n) lv.the_array[n] = 0.0;
                lv.the_array[0] = ok ? kc_g(m,0,i,k,j) : 0.0;
                lv.the_array[1] = ok ? 1.0 : 0.0;
                ls += lv;
              }, Kokkos::Sum<array_sum::GlobalSum>(sm));
              Kokkos::single(Kokkos::PerTeam(tm), [&]() {
                kfd(tm.league_rank()) = sm.the_array[0];
                kfd(nr + tm.league_rank()) = sm.the_array[1];
              });
            });
            Kokkos::fence();
            Kokkos::deep_copy(kfh, kfd);
#if MPI_PARALLEL_ENABLED
            MPI_Allreduce(MPI_IN_PLACE, kfh.data(), 2*nr, MPI_ATHENA_REAL, MPI_SUM,
                          MPI_COMM_WORLD);
#endif
            for (int r=0; r<nr; ++r) {
              kfh(r) = (kfh(nr + r) > 0.0) ? (kfh(r)/kfh(nr + r)) : 0.0;
            }
            Kokkos::deep_copy(kfd, kfh);
            par_for("rt_kfrz_put", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
            KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
              if (i < icut_g(m,k,j)) return;
              if (!(T_g(m,k,j,i) > 0.0)) return;
              kc_g(m,0,i,k,j) = kfd(i - is);
            });
          }
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
            // THE DEEP-LIMIT GRADIENT AT THE CUT (see rt_cut_bc_legacy).  dB/dtau with
            // tau increasing downward, from the two cell centres straddling the cut face
            // and the Rosseland tau between those centres.  Zero in the legacy mode, and
            // zero if the cut sits on the top cell (no i+1 to difference against).
            Real dbdtau_cut = 0.0;
            if (!cut_legacy && icut + 1 <= ie) {
              const Real krb0 = kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut);
              const Real krb1 = kc_g(m,0,icut+1,k,j)*rhoN(m,k,j,icut+1);
              const Real dtc = 0.5*(krb0*DX1(m,k,j,icut) + krb1*DX1(m,k,j,icut+1));
              if (dtc > 0.0) {
                dbdtau_cut = (Bb_g(m,0,icut,k,j) - Bb_g(m,0,icut+1,k,j))/dtc;
              }
            }
            // problem/rt_bottom_flux: the cut IS the bottom wall and carries the
            // imposed internal flux, so its gradient is set by that flux, not by the
            // column's own two deepest cells.  See rt_bot_flux.
            if (bot_flux > 0.0) dbdtau_cut = 3.0*bot_flux/(4.0*M_PI);
            // ck_nquad = 1 is the hemispheric mean (mu = 1/1.66), 2 the two-point
            // Gauss-Legendre quadrature the band solver offers on the same switch
            const int nq = (ck_nq_ > 1) ? 2 : 1;
            // THE LAYER between the centres of cells iL and iL+1 (see rt_layer_legacy):
            // its two half thicknesses, the Planck source at each centre, and the source
            // at the face that separates them, taken at its own tau inside the layer.
            // BFace is applied to BOTH endpoints and symmetrically, so a radiatively
            // inert neighbour leaves a layer emitting with its own Planck function, as it
            // did before -- and with equal opacities on both sides this is just the
            // straight line from B(iL) to B(iL+1).
            auto rt_layer = [&](const int iL, Real &dt_l, Real &dt_u,
                                Real &s_l, Real &s_u, Real &s_f) {
              const Real kl = kc_g(m,0,iL,k,j)*rhoN(m,k,j,iL);
              const Real ku = kc_g(m,0,iL+1,k,j)*rhoN(m,k,j,iL+1);
              const Real bl = Bb_g(m,0,iL,k,j);
              const Real bu = Bb_g(m,0,iL+1,k,j);
              dt_l = 0.5*kl*DX1(m,k,j,iL);
              dt_u = 0.5*ku*DX1(m,k,j,iL+1);
              s_l = BFace(ku, kl, bu, bl, bface_on);
              s_u = BFace(kl, ku, bl, bu, bface_on);
              const Real dtc = dt_l + dt_u;
              s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc)) : (0.5*(s_l + s_u));
            };
            // ---- rt_implicit_column: the SAME layer, differentiated -----------------
            // Every source above is a convex combination of the two centre Planck
            // functions, so four weights describe the layer completely:
            //   ds_l/dB_l = p_l,  ds_l/dB_u = 1 - p_l,
            //   ds_u/dB_u = p_u,  ds_u/dB_l = 1 - p_u
            // and s_f inherits them through the same dt_l/(dt_l+dt_u) interpolation.
            const bool jac_on = implcol_;
            auto rt_layer_w = [&](const int iL, Real &dl_l, Real &dl_u, Real &du_l,
                                  Real &du_u, Real &df_l, Real &df_u) {
              const Real kl = kc_g(m,0,iL,k,j)*rhoN(m,k,j,iL);
              const Real ku = kc_g(m,0,iL+1,k,j)*rhoN(m,k,j,iL+1);
              const Real pl = BFaceW(ku, kl, bface_on);     // ds_l/dB_l
              const Real pu = BFaceW(kl, ku, bface_on);     // ds_u/dB_u
              const Real dtl = 0.5*kl*DX1(m,k,j,iL);
              const Real dtu = 0.5*ku*DX1(m,k,j,iL+1);
              const Real dtc = dtl + dtu;
              const Real f = (dtc > 0.0) ? (dtl/dtc) : 0.5;
              dl_l = pl;
              dl_u = 1.0 - pl;
              du_l = 1.0 - pu;
              du_u = pu;
              df_l = dl_l*(1.0 - f) + du_l*f;
              df_u = dl_u*(1.0 - f) + du_u*f;
            };
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
            // AREA-WEIGHTED INTENSITIES from here on: J = A I, the SPHERICAL
            // DILUTION of the caller's note.  The sweep integrates the same layer
            // solve, mu dJ/dr = -kappa rho (J - A B), so every source endpoint is
            // multiplied by the area AT ITS OWN POSITION (a face endpoint by AFC, a
            // cell-centre endpoint by the cell's mean area ACC), every deposit is per
            // VOLUME (VLS), and a face flux handed back is J/A again.  Under
            // rt_plane_parallel all of that is A = 1, V = dx1 and bitwise the old code.
            const Real aft_ = AFC(m,k,j,ie+1);
            const Real afc_ = AFC(m,k,j,icut);
            Real I_down[2][NN];
            // Top: the unresolved hydrostatic column above the domain, p/g of it, at the
            // top cell's opacity -- the same construction the band solver uses.
            {
              const Real kap = kc_g(m,0,ie+1,k,j);
              const Real mu0 = cf_g(m,k,j,3);
              // rt_top_vacuum: nothing above the domain, so the column has no optical
              // depth and hands back no intensity whatever its source function is.
              const Real dtau = top_vac ? 0.0
                  : RTTopDtau(kap, pb_g(m,k,j,ie+1)*1.0e6,
                              EffGravAt(grav, ap, X1V(m,ie+1), grav_pmass,
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
              if (top_re && !layer_legacy) {
                // the same probe on the centre-to-centre layers (see rt_layer_legacy):
                // the cut half, then every layer, then the top half
                const Real dcut = 0.5*kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut)
                                * DX1(m,k,j,icut);
                const Real bcut = Bb_g(m,0,icut,k,j);
                const Real dtop = 0.5*kc_g(m,0,ie,k,j)*rhoN(m,k,j,ie)*DX1(m,k,j,ie);
                const Real btop = Bb_g(m,0,ie,k,j);
                // the probe is the up-sweep, so it dilutes exactly as the up-sweep
                // does: it carries J and is divided back by the top face's area
                const Real acc_ = ACC(m,k,j,icut), act_ = ACC(m,k,j,ie);
                for (int q=0; q<nq; ++q) {
                  const Real bfc = bcut + dbdtau_cut*dcut;
                  Real ip = (bfc + (int_at_cut ? Iint : 0.0)
                             + muq[q]*dbdtau_cut)*afc_;
                  Real ab, em;
                  RTLayer(dcut, muq[q], bfc*afc_, bcut*acc_, ip, ab, em);
                  for (int i=icut; i<ie; ++i) {
                    Real dt_l, dt_u, s_l, s_u, s_f;
                    rt_layer(i, dt_l, dt_u, s_l, s_u, s_f);
                    const Real acl = ACC(m,k,j,i), afm = AFC(m,k,j,i+1);
                    const Real acu = ACC(m,k,j,i+1);
                    RTLayer(dt_l, muq[q], s_l*acl, s_f*afm, ip, ab, em);
                    RTLayer(dt_u, muq[q], s_f*afm, s_u*acu, ip, ab, em);
                  }
                  RTLayer(dtop, muq[q], btop*act_, btop*aft_, ip, ab, em);
                  bsrc[q] = 0.5*ip/aft_;
                }
              } else if (top_re) {
                for (int q=0; q<nq; ++q) {
                  Real ip = (Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0)
                          + muq[q]*dbdtau_cut)*afc_;
                  for (int i=icut+1; i<ie+2; ++i) {
                    const Real krb = kc_g(m,0,i-1,k,j)*rhoN(m,k,j,i-1);
                    const Real x = krb*DX1(m,k,j,i-1)/muq[q];
                    const Real e0 = -expm1(-x);
                    const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                    const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                                  : (x/2.0 - SQR(x)/3.0);
                    // emissivity-weighted far endpoint: see BFace
                    const int iir = (i > ie) ? ie : i;
                    const Real bfar = BFace(krb, kc_g(m,0,i,k,j)*rhoN(m,k,j,iir),
                                            Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
                    ip = (1.0-e0)*ip + bet*bfar*ACC(m,k,j,i)
                       + gm*Bb_g(m,0,i-1,k,j)*ACC(m,k,j,i-1);
                  }
                  bsrc[q] = 0.5*ip/aft_;
                }
              }
              for (int q=0; q<nq; ++q) {
                I_down[q][ie+1] = (1.0 - exp(-dtau/muq[q]))*(bsrc[q]*aft_);
              }
              if (report_on) idn_g(m,k,j,ie+1) = I_down[0][ie+1];
            }
            // THE LAYERS.  See rt_layer_legacy: the default branch runs layers
            // between cell CENTRES; the legacy one the old staggered whole-cell
            // layers, bit for bit.
            Real I_up[2];
            if (layer_legacy) {
              // down-sweep on the OLD staggered whole-cell layers
              for (int i=ie; i>icut-1; --i) {
                const Real krb_d = kc_g(m,0,i,k,j)*rhoN(m,k,j,i);
                const Real dtau_i = krb_d*DX1(m,k,j,i);
                // the far endpoint of the layer's source function, weighted by emitting
                // matter: above rad_kappa_rmax the neighbour has kappa = 0 and a Planck
                // function 1e9x this cell's, which used to drain it to the floor in one
                // step (see BFace).  Identical to Bb_g(i+1) at equal opacity.
                const int iip = (i+1 > ie) ? ie : i+1;
                const Real bfar_d = BFace(krb_d, kc_g(m,0,i+1,k,j)*rhoN(m,k,j,iip),
                                          Bb_g(m,0,i,k,j), Bb_g(m,0,i+1,k,j), bface_on);
                const Real bfa_d = bfar_d*ACC(m,k,j,i+1);
                const Real bwn_d = Bb_g(m,0,i,k,j)*ACC(m,k,j,i);
                for (int q=0; q<nq; ++q) {
                  const Real x = dtau_i/muq[q];
                  const Real e0 = -expm1(-x);
                  const Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                               : (x/2.0 - SQR(x)/3.0);
                  const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                  // direct source: what this stream leaves in cell i, absorbed minus
                  // emitted
                  Src_g(m,0,i,k,j) += wfq[q]/VLA(m,k,j,i)
                                    *(e0*I_down[q][i+1]
                                      - (alp*bfa_d + bet*bwn_d));
                  I_down[q][i] = (1.0-e0)*I_down[q][i+1]
                               + alp*bfa_d + bet*bwn_d;
                }
                if (report_on) idn_g(m,k,j,i) = I_down[0][i];
              }
              // Bottom of the RT domain: thermalised, plus the internal flux if the
              // layers below are not carrying it themselves (see rt_int_at_cut).
              for (int q=0; q<nq; ++q) {
                I_up[q] = (Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0)
                        + muq[q]*dbdtau_cut)*afc_;
                Fb_g(m,0,icut,k,j) += wfq[q]*(I_up[q] - I_down[q][icut])/afc_;
              }
              if (report_on) iup_g(m,k,j,icut) = I_up[0];
              // up-sweep
              for (int i=icut+1; i<ie+2; ++i) {
                const Real kap = kc_g(m,0,i-1,k,j);
                const Real rho = rhoN(m,k,j,i-1);
                const Real dtau_i = kap*rho*DX1(m,k,j,i-1);
                // same emissivity weighting for the upward stream and for Em: see BFace
                const int iiu = (i > ie) ? ie : i;
                const Real bfar_u = BFace(kap*rho, kc_g(m,0,i,k,j)*rhoN(m,k,j,iiu),
                                          Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
                const Real bfa_u = bfar_u*ACC(m,k,j,i);
                const Real bwn_u = Bb_g(m,0,i-1,k,j)*ACC(m,k,j,i-1);
                const Real afi_ = AFC(m,k,j,i);
                for (int q=0; q<nq; ++q) {
                  const Real x = dtau_i/muq[q];
                  const Real e0 = -expm1(-x);
                  const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                  const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                               : (x/2.0 - SQR(x)/3.0);
                  const Real Iup_in = I_up[q];
                  // the same for the upward stream through layer i-1
                  Src_g(m,0,i-1,k,j) += wfq[q]/VLA(m,k,j,i-1)
                                      *(e0*Iup_in
                                        - (bet*bfa_u
                                           + gm*bwn_u));
                  I_up[q] = (1.0-e0)*Iup_in
                          + bet*bfa_u + gm*bwn_u;
                  Fb_g(m,0,i,k,j) += wfq[q]*(I_up[q] - I_down[q][i])/afi_;
                }
                if (report_on) iup_g(m,k,j,i) = I_up[0];
                Em_g(m,0,i-1,k,j) = 4.0*M_PI*kap*rho
                                  * 0.5*(bfar_u + Bb_g(m,0,i-1,k,j));
              }
            } else {
              // ---- the centre-to-centre layers (see rt_layer_legacy) ----
              // The top half layer: the upper half of cell ie, entered at the top face,
              // with the source held at the cell's own centre value -- there is nothing
              // above it to interpolate towards.
              const Real dt_top = 0.5*kc_g(m,0,ie,k,j)*rhoN(m,k,j,ie)*DX1(m,k,j,ie);
              const Real b_top = Bb_g(m,0,ie,k,j);
              // The cut half layer: the lower half of cell icut.  B is continued from
              // the centre to the cut face at the deep-limit gradient, the same
              // linear-in-tau behaviour the cut boundary itself assumes (and zero under
              // rt_cut_bc_legacy).  A thick half layer forgets its entry source anyway,
              // since that coefficient falls off as 1/x.
              const Real dt_cut = 0.5*kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut)
                                * DX1(m,k,j,icut);
              const Real b_cut = Bb_g(m,0,icut,k,j);
              const Real b_cutf = b_cut + dbdtau_cut*dt_cut;
              Real ab, em;
              // down-sweep: centre to centre, recording the face intensity in between
              Real Idn[2];
              // dI_down(current face)/dB of the two cells nearest below it, per stream
              Real ddn_c[2], ddn_p[2];
              if (jac_on) {
                for (int i=is; i<ie+2; ++i) {
                  jac_g(m,0,k,j,i) = 0.0;
                  jac_g(m,1,k,j,i) = 0.0;
                  jac_g(m,2,k,j,i) = 0.0;
                }
                for (int q=0; q<2; ++q) {
                  ddn_c[q] = 0.0;
                  ddn_p[q] = 0.0;
                }
              }
              // the top half layer runs from the TOP FACE to the CENTRE of cell ie
              const Real act_ = ACC(m,k,j,ie), acc_ = ACC(m,k,j,icut);
              for (int q=0; q<nq; ++q) {
                Idn[q] = I_down[q][ie+1];
                RTLayer(dt_top, muq[q], b_top*aft_, b_top*act_, Idn[q], ab, em);
                Src_g(m,0,ie,k,j) += wfq[q]/VLA(m,k,j,ie)*(ab - em);
                if (jac_on) {
                  // the top half layer: both endpoints are cell ie's own B
                  Real e0t, cit, cot;
                  RTLayerCoef(dt_top, muq[q], e0t, cit, cot);
                  const Real W = wfq[q]/DX1(m,k,j,ie);
                  jac_g(m,1,k,j,ie) -= W*(cit + cot);
                  ddn_c[q] = cit + cot;      // dIdn/dB(ie), the cell it just crossed
                  ddn_p[q] = 0.0;
                }
              }
              for (int i=ie; i>icut; --i) {
                Real dt_l, dt_u, s_l, s_u, s_f;
                rt_layer(i-1, dt_l, dt_u, s_l, s_u, s_f);
                Real wl_l = 0.0, wl_u = 0.0, wu_l = 0.0, wu_u = 0.0;
                Real wf_l = 0.0, wf_u = 0.0;
                if (jac_on) rt_layer_w(i-1, wl_l, wl_u, wu_l, wu_u, wf_l, wf_u);
                const Real acd_u = ACC(m,k,j,i), afd_m = AFC(m,k,j,i);
                const Real acd_l = ACC(m,k,j,i-1);
                for (int q=0; q<nq; ++q) {
                  RTLayer(dt_u, muq[q], s_u*acd_u, s_f*afd_m, Idn[q], ab, em);
                  Src_g(m,0,i,k,j) += wfq[q]/VLA(m,k,j,i)*(ab - em);
                  I_down[q][i] = Idn[q];
                  RTLayer(dt_l, muq[q], s_f*afd_m, s_l*acd_l, Idn[q], ab, em);
                  Src_g(m,0,i-1,k,j) += wfq[q]/VLA(m,k,j,i-1)*(ab - em);
                  if (jac_on) {
                    // upper half: enters at s_u (centre i), leaves at s_f; deposits in i
                    Real e0u, ciu, cou, e0l, cil, col;
                    RTLayerCoef(dt_u, muq[q], e0u, ciu, cou);
                    RTLayerCoef(dt_l, muq[q], e0l, cil, col);
                    const Real Wi = wfq[q]/DX1(m,k,j,i);
                    const Real Wm = wfq[q]/DX1(m,k,j,i-1);
                    const Real eu_u = ciu*wu_u + cou*wf_u;   // d(emitted upper)/dB_i
                    const Real eu_l = ciu*wu_l + cou*wf_l;   // ... /dB_{i-1}
                    jac_g(m,1,k,j,i)   += Wi*(e0u*ddn_c[q] - eu_u);
                    jac_g(m,2,k,j,i)   += Wi*e0u*ddn_p[q];
                    jac_g(m,0,k,j,i)   -= Wi*eu_l;
                    // the intensity at the face between the two half layers
                    const Real dm_u = (1.0-e0u)*ddn_c[q] + eu_u;
                    const Real dm_l = eu_l;
                    const Real el_l = cil*wf_l + col*wl_l;   // d(emitted lower)/dB_{i-1}
                    const Real el_u = cil*wf_u + col*wl_u;   // ... /dB_i
                    jac_g(m,1,k,j,i-1) += Wm*(e0l*dm_l - el_l);
                    jac_g(m,2,k,j,i-1) += Wm*(e0l*dm_u - el_u);
                    // what leaves the layer, for the next pair one cell down
                    const Real do_l = (1.0-e0l)*dm_l + el_l;
                    const Real do_u = (1.0-e0l)*dm_u + el_u;
                    ddn_c[q] = do_l;
                    ddn_p[q] = do_u;
                  }
                }
                if (report_on) idn_g(m,k,j,i) = I_down[0][i];
              }
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_cut, muq[q], b_cut*acc_, b_cutf*afc_, Idn[q], ab, em);
                Src_g(m,0,icut,k,j) += wfq[q]/VLA(m,k,j,icut)*(ab - em);
                I_down[q][icut] = Idn[q];
              }
              if (report_on) idn_g(m,k,j,icut) = I_down[0][icut];
              // Bottom of the RT domain, now AT the cut face: thermalised, plus the
              // internal flux if the layers below are not carrying it themselves (see
              // rt_int_at_cut), plus the deep-limit gradient (see rt_cut_bc_legacy).
              for (int q=0; q<nq; ++q) {
                I_up[q] = (b_cutf + (int_at_cut ? Iint : 0.0)
                           + muq[q]*dbdtau_cut)*afc_;
                Fb_g(m,0,icut,k,j) += wfq[q]*(I_up[q] - I_down[q][icut])/afc_;
              }
              if (report_on) iup_g(m,k,j,icut) = I_up[0];
              // up-sweep
              // dI_up(current face)/dB of the two cells nearest below it, per stream
              Real dup_c[2], dup_m[2];
              if (jac_on) {
                for (int q=0; q<2; ++q) {
                  dup_c[q] = 0.0;
                  dup_m[q] = 0.0;
                }
              }
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_cut, muq[q], b_cutf*afc_, b_cut*acc_, I_up[q], ab, em);
                Src_g(m,0,icut,k,j) += wfq[q]/VLA(m,k,j,icut)*(ab - em);
                if (jac_on) {
                  // the cut half layer: both endpoints are cell icut's own B (b_cutf
                  // adds only the frozen deep gradient), and I_up entered it from the
                  // thermalised boundary, which is B(icut) as well
                  Real e0c, cic, coc;
                  RTLayerCoef(dt_cut, muq[q], e0c, cic, coc);
                  const Real W = wfq[q]/DX1(m,k,j,icut);
                  jac_g(m,1,k,j,icut) += W*(e0c - (cic + coc));
                  dup_c[q] = (1.0-e0c) + cic + coc;
                  dup_m[q] = 0.0;
                }
              }
              for (int i=icut; i<ie; ++i) {
                Real dt_l, dt_u, s_l, s_u, s_f;
                rt_layer(i, dt_l, dt_u, s_l, s_u, s_f);
                Real wl_l = 0.0, wl_u = 0.0, wu_l = 0.0, wu_u = 0.0;
                Real wf_l = 0.0, wf_u = 0.0;
                if (jac_on) rt_layer_w(i, wl_l, wl_u, wu_l, wu_u, wf_l, wf_u);
                const Real acu_l = ACC(m,k,j,i), afu_m = AFC(m,k,j,i+1);
                const Real acu_u = ACC(m,k,j,i+1);
                for (int q=0; q<nq; ++q) {
                  RTLayer(dt_l, muq[q], s_l*acu_l, s_f*afu_m, I_up[q], ab, em);
                  Src_g(m,0,i,k,j) += wfq[q]/VLA(m,k,j,i)*(ab - em);
                  Fb_g(m,0,i+1,k,j) += wfq[q]*(I_up[q] - I_down[q][i+1])/afu_m;
                  if (report_on && q == 0) iup_g(m,k,j,i+1) = I_up[0];
                  RTLayer(dt_u, muq[q], s_f*afu_m, s_u*acu_u, I_up[q], ab, em);
                  Src_g(m,0,i+1,k,j) += wfq[q]/VLA(m,k,j,i+1)*(ab - em);
                  if (jac_on) {
                    Real e0l, cil, col, e0u, ciu, cou;
                    RTLayerCoef(dt_l, muq[q], e0l, cil, col);
                    RTLayerCoef(dt_u, muq[q], e0u, ciu, cou);
                    const Real Wi = wfq[q]/DX1(m,k,j,i);
                    const Real Wp = wfq[q]/DX1(m,k,j,i+1);
                    const Real el_l = cil*wl_l + col*wf_l;   // d(emitted lower)/dB_i
                    const Real el_u = cil*wl_u + col*wf_u;   // ... /dB_{i+1}
                    jac_g(m,1,k,j,i) += Wi*(e0l*dup_c[q] - el_l);
                    jac_g(m,0,k,j,i) += Wi*e0l*dup_m[q];
                    jac_g(m,2,k,j,i) -= Wi*el_u;
                    const Real dm_c = (1.0-e0l)*dup_c[q] + el_l;
                    const Real dm_u = el_u;
                    const Real eu_u = ciu*wf_u + cou*wu_u;   // d(emitted upper)/dB_{i+1}
                    const Real eu_l = ciu*wf_l + cou*wu_l;   // ... /dB_i
                    jac_g(m,1,k,j,i+1) += Wp*(e0u*dm_u - eu_u);
                    jac_g(m,0,k,j,i+1) += Wp*(e0u*dm_c - eu_l);
                    const Real do_u = (1.0-e0u)*dm_u + eu_u;
                    const Real do_c = (1.0-e0u)*dm_c + eu_l;
                    dup_c[q] = do_u;
                    dup_m[q] = do_c;
                  }
                }
              }
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_top, muq[q], b_top*act_, b_top*aft_, I_up[q], ab, em);
                Src_g(m,0,ie,k,j) += wfq[q]/VLA(m,k,j,ie)*(ab - em);
                Fb_g(m,0,ie+1,k,j) += wfq[q]*(I_up[q] - I_down[q][ie+1])/aft_;
                if (jac_on) {
                  Real e0t, cit, cot;
                  RTLayerCoef(dt_top, muq[q], e0t, cit, cot);
                  const Real W = wfq[q]/DX1(m,k,j,ie);
                  jac_g(m,1,k,j,ie) += W*(e0t*dup_c[q] - (cit + cot));
                  jac_g(m,0,k,j,ie) += W*e0t*dup_m[q];
                }
              }
              if (report_on) iup_g(m,k,j,ie+1) = I_up[0];
              // the cell's own emission, exactly 4 sigma kappa rho T^4: each cell now
              // owns both of its half layers, so nothing here is a two-centre average
              for (int i=icut; i<ie+1; ++i) {
                Em_g(m,0,i,k,j) = 4.0*M_PI*kc_g(m,0,i,k,j)*rhoN(m,k,j,i)
                                * Bb_g(m,0,i,k,j);
              }
            }
          });
        };
        if (skipsweep_) {
          // NOTHING to launch: kc_g/Bb_g/icut_g come from the pre-kernels, Qb_g is
          // identically 0 on the grey path (only the correlated-k chain accumulates
          // into it) and stays at its allocation zero, Fb_g is written by the column
          // solve below, and Src_g/Em_g stay zero with no consumer left.
        } else if (n1 <= 72) {
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
        // ---- rt_implicit_column = 3: THE EXACT IMPLICIT COLUMN SOLVE ----------------
        // The sweep above has just filled kc_g/Bb_g (the FROZEN opacity and the entry
        // Planck function), Fb_g/Src_g (which the handover term is formed from) and
        // Qb_g.  Everything the block solve needs is therefore in place; it re-solves
        // the same column with the intensities as unknowns and applies the energy
        // itself.  Fb_g/Src_g/Em_g are deliberately LEFT at the entry-state values, so
        // rad_f2s, rt_rad_force and every diagnostic below see exactly the numbers mode
        // 0 hands them.
        if (mode3_) {
          const int nmb_c3 = pmbp->nmb_thispack;
          // THE PARTITIONED SOLVER (problem/rt_impl_solver = pcr).  A team per column
          // needs the spike H and the two sweep factors on top of the serial layout, and
          // a per-segment array for the reduced system.  On a host backend a Kokkos team
          // is one thread, so the partition is forced to a single segment there -- which
          // is the serial block Thomas, to round-off.
          const bool c3par = (rt_impl_solver == 1);
          const bool c3rp = c3par && rt_impl_redpar;
          // problem/rt_impl_reuse enlarges BOTH workspaces (the per-cell factors v, M
          // and the reduced solve's own stored factors); with the switch off not a word
          // is added.  The serial (thomas) path does not implement it.
          if (rt_impl_mixed > 0 && !c3par) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_impl_mixed "
                      << "requires problem/rt_impl_solver = pcr (the partitioned column "
                      << "solver); the thomas path is double only." << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (rt_impl_reuse > 0 && !c3par) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_impl_reuse "
                      << "requires problem/rt_impl_solver = pcr (the partitioned column "
                      << "solver); the thomas path does not store the factorisation."
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
          const int c3nrd = (rt_impl_reuse > 0) ? (c3rp ? RTCOL3_NRDPU : RTCOL3_NRDU)
                                                : (c3rp ? RTCOL3_NRDP : RTCOL3_NRD);
          int c3nseg = c3par ? rt_impl_nseg : 1;
          if (std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
            c3nseg = 1;
          }
          const int c3nw = c3par ? ((rt_impl_reuse > 0) ? RTCOL3_NWPU : RTCOL3_NWP)
                                 : RTCOL3_NW;
          if (rt_c3wk_ptr == nullptr) {
            // the partitioned path TRANSPOSES the workspace so that the fast index is
            // the one the threads of a team differ in; see RTCol3::Wk
            rt_c3wk_ptr = c3par
                ? new DvceArray5D<Real>("rt_c3wk", nmb_c3, c3nw, n3, n2, n1)
                : new DvceArray5D<Real>("rt_c3wk", nmb_c3, c3nw, n1, n3, n2);
            rt_c3top_ptr = new DvceArray4D<Real>("rt_c3top", nmb_c3, 2, n3, n2);
            rt_c3stat_ptr = new DvceArray1D<Real>("rt_c3stat", 26);
            // the warm-start history.  Zero-initialised, so "no history" is the state of
            // every cell on the first call and the fallback fires there by construction.
            const int wn1 = (rt_impl_warm > 0) ? n1 : 1;
            const int wn2 = (rt_impl_warm > 0) ? n2 : 1;
            const int wn3 = (rt_impl_warm > 0) ? n3 : 1;
            const int wnm = (rt_impl_warm > 0) ? nmb_c3 : 1;
            rt_c3bp_ptr = new DvceArray4D<Real>("rt_c3bp", wnm, wn3, wn2, wn1);
            rt_c3bp2_ptr = new DvceArray4D<Real>("rt_c3bp2",
                                                 (rt_impl_warm > 1) ? wnm : 1,
                                                 (rt_impl_warm > 1) ? wn3 : 1,
                                                 (rt_impl_warm > 1) ? wn2 : 1,
                                                 (rt_impl_warm > 1) ? wn1 : 1);
            rt_c3rd_ptr = new DvceArray4D<Real>("rt_c3rd", nmb_c3, n3, n2,
                                                c3par ? c3nseg*c3nrd : 1);
            // problem/rt_impl_mixed = 2: the compact SINGLE-precision factor workspace.
            // Sized 1 when the switch is not 2, so nothing is added otherwise.
            const int c3nwf = (rt_impl_mixed == 2)
                ? ((rt_impl_reuse > 0) ? RTCOL3_NWFU : RTCOL3_NWF) : 1;
            rt_c3wkf_ptr = (rt_impl_mixed == 2)
                ? new DvceArray5D<float>("rt_c3wkf", nmb_c3, c3nwf, n3, n2, n1)
                : new DvceArray5D<float>("rt_c3wkf", 1, 1, 1, 1, 1);
            if (global_variable::my_rank == 0) {
              const double wmb = static_cast<double>(nmb_c3)*c3nw*n1*n2*n3
                                 *sizeof(Real)/1.0e6;
              const double rmb = c3par ? (static_cast<double>(nmb_c3)*c3nseg*c3nrd
                                          *n2*n3*sizeof(Real)/1.0e6) : 0.0;
              std::cout << "### two_stream_rt: rt_implicit_column = 3, the EXACT "
                        << "block-tridiagonal column solve (intensities as unknowns); "
                        << "solver " << (c3par ? "pcr" : "thomas")
                        << " nseg " << c3nseg
                        << (c3rp ? " redpar" : "")
                        << " hybrid_tau " << rt_col3_hybrid_tau
                        << (rt_col3_split_deep ? " split_deep" : "")
                        << (rt_impl_reuse > 0
                            ? ((rt_impl_reuse > 1) ? " reuse(always)" : " reuse") : "")
                        << (rt_impl_mixed > 0
                            ? ((rt_impl_mixed > 1) ? " mixed(2)" : " mixed(1)") : "")
                        << "; workspace " << wmb << " + " << rmb << " MB" << std::endl;
            }
          }
          auto c3wk = *rt_c3wk_ptr;
          auto c3wkf = *rt_c3wkf_ptr;
          auto c3rd = *rt_c3rd_ptr;
          auto c3top = *rt_c3top_ptr;
          auto c3stat = *rt_c3stat_ptr;
          auto c3bp = *rt_c3bp_ptr;
          auto c3bp2 = *rt_c3bp2_ptr;
          Kokkos::deep_copy(c3stat, 0.0);
          // the angular quadrature, the same one the sweep just used
          const int nq3 = (ck_nq_ > 1) ? 2 : 1;
          Real mu3[2], wf3[2];
          if (nq3 == 1) {
            mu3[0] = 1.0/CK_DIFFUSIVITY;
            wf3[0] = M_PI;
            mu3[1] = 1.0;
            wf3[1] = 0.0;
          } else {
            for (int q=0; q<2; ++q) {
              mu3[q] = mug[q];
              wf3[q] = 2.0*M_PI*wg[q]*mug[q];
            }
          }
          // the frozen top-face downward intensity: the sweep's unresolved-column model
          const Real mu3a = mu3[0], mu3b = mu3[1];
          const bool tv3 = top_vac;
          par_for("rt_c3_top", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            const Real kap = kc_g(m,0,ie+1,k,j);
            const Real mu0 = cf_g(m,k,j,3);
            const Real dtau = tv3 ? 0.0
                : RTTopDtau(kap, pb_g(m,k,j,ie+1)*1.0e6,
                            EffGravAt(grav, ap, X1V(m,ie+1), grav_pmass, omega, mu0,
                                      tide));
            const Real bsrc = Bb_g(m,0,ie+1,k,j);
            c3top(m,0,k,j) = (1.0 - exp(-dtau/mu3a))*bsrc;
            c3top(m,1,k,j) = (1.0 - exp(-dtau/mu3b))*bsrc;
          });
          // ---- problem/rt_force_center: centre the flux the radiative force uses --
          // The column below writes its OWN converged face flux into Fb (the
          // rt_col3_skip_sweep path, RTCol3::wrflux); mode 2 needs the ENTRY flux as
          // well, so it is summed over the blocks and saved here, before the solve.
          const int fcen_ = rt_force_center;
          if (fcen_ == 2 && rt_fbsave_ptr == nullptr) {
            rt_fbsave_ptr = new DvceArray4D<Real>("rt_fbsave", nmb_c3, n3, n2, n1);
          }
          auto fbsv_ = (rt_fbsave_ptr != nullptr) ? *rt_fbsave_ptr
                     : DvceArray4D<Real>("rt_fbsave_d", 1, 1, 1, 1);
          if (fcen_ == 2) {
            const int nblk_s = nblk;
            par_for("rt_fbsave", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
            KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
              Real ff = 0.0;
              for (int b=0; b<nblk_s; ++b) ff += Fb_g(m,b,i,k,j);
              fbsv_(m,k,j,i) = ff;
            });
          }
          RTCol3 c3;
          c3.u0 = u0;
          c3.bcc = bcc_uc_;
          c3.phicc = phicc_uc_;
          c3.cosc = cosc_uc_;
          c3.kc = kc_g;
          c3.Bb = Bb_g;
          c3.Fb = Fb_g;
          c3.Src = Src_g;
          c3.Qb = Qb_g;
          c3.Tg = T_g;
          c3.wblend = taublend ? w_g : DvceArray4D<Real>("rt_c3w_d",1,1,1,1);
          c3.icut = icut_g;
          c3.wk = c3wk;
          c3.wkf = c3wkf;
          c3.mixed = rt_impl_mixed;
          c3.rd = c3rd;
          c3.nseg = c3nseg;
          c3.redpar = c3rp;
          c3.dtop = c3top;
          c3.stat = c3stat;
          c3.size = size;
          c3.dx1 = dx1_;
          // the radial geometry of the SPHERICAL DILUTION (see the note above): the
          // column solve carries the same area-weighted intensities the sweep does
          c3.vol = volume;
          c3.area1 = area1;
          c3.eos = eos;
          c3.bdt = bdt;
          c3.sigma = boltz_sigma;
          c3.Iint = Iint;
          c3.bot_flux = bot_flux;
          c3.mu[0] = mu3[0];
          c3.mu[1] = mu3[1];
          c3.wf[0] = wf3[0];
          c3.wf[1] = wf3[1];
          c3.tol = rt_impl_tol;
          c3.norm = rt_impl_norm;
          c3.norm_eps = rt_impl_norm_eps;
          c3.exjac = rt_impl_exjac;
          c3.dstop = rt_impl_dstop;
          c3.rescheck = rt_impl_rescheck;
          c3.ablate = rt_impl_ablate;
          c3.fixit = rt_impl_fixit;
          c3.cvfreeze = rt_impl_cvfreeze;
          c3.reuse = rt_impl_reuse;
          c3.reuse_rho = rt_impl_reuse_rho;
          c3.warm = rt_impl_warm;
          c3.bprev = c3bp;
          c3.bprev2 = c3bp2;
          // the extrapolation ratio: this call's bdt over the last one's.  The first
          // call of a run has no history at all, so the value is irrelevant there.
          c3.dtr = (rt_c3_bdt_prev > 0.0) ? (bdt/rt_c3_bdt_prev) : 1.0;
          rt_c3_bdt_prev = bdt;
          c3.dfloor = eos.dfloor;
          c3.rgas = Rgas;
          c3.gm1 = gm1;
          c3.nq = nq3;
          c3.maxit = (rt_impl_maxit > 0) ? rt_impl_maxit : 8;
          c3.is = is;
          c3.ie = ie;
          c3.is_pp = is_pp;
          c3.nx1_pp = nx1_pp;
          c3.pp = pp_;
          c3.mhd = mhd_uc_;
          c3.etg = etg_uc_;
          c3.cs = cs_uc_;
          c3.bface = bface_on;
          c3.taublend = taublend;
          c3.int_at_cut = int_at_cut;
          c3.cut_legacy = cut_legacy;
          c3.direct = rt_src_direct;
          c3.ex_iter = rt_col3_ex_iter;
          c3.hyb_tau = rt_col3_hybrid_tau;
          c3.split_deep = rt_col3_split_deep;
          c3.split_w = rt_col3_split_w;
          c3.wrflux = skipsweep_ || (fcen_ > 0);
          c3.theta = rt_src_theta;
          c3.dump = rt_outer_verbose && (pm->ncycle == 0);
          if (c3par) {
            RTCol3TeamLaunch(c3, nmb1, ks, ke, js, je);
          } else {
            RTCol3Launch(c3, nmb1, ks, ke, js, je);
          }
          if (rt_outer_verbose && global_variable::my_rank == 0 &&
              rt_c3_lines < 4000 &&
              (rt_report_every <= 0 || pm->ncycle % rt_report_every == 0)) {
            ++rt_c3_lines;
            auto hs = Kokkos::create_mirror_view(c3stat);
            Kokkos::deep_copy(hs, c3stat);
            const Real ncol = (hs(1) > 0.0) ? hs(1) : 1.0;
            std::cout << "### rt_col3 ncycle=" << pm->ncycle << " t=" << pm->time
                      << " dt=" << bdt << " ncol=" << static_cast<int>(hs(1))
                      << " it_mean=" << hs(0)/ncol
                      << " it_max=" << static_cast<int>(hs(2))
                      << " nclamp=" << static_cast<int>(hs(3))
                      << " max|db/b|_it=" << hs(4)
                      << " max|du/u|_tot=" << hs(18)
                      << " eos_roundtrip=" << hs(5)
                      << " nsing=" << static_cast<int>(hs(6))
                      << " nbadE=" << static_cast<int>(hs(7))
                      << " nfloor=" << static_cast<int>(hs(8))
                      << " sum_de_dx=" << hs(9)
                      << " sum|de|dx=" << hs(10)
                      << " ndiagviol=" << static_cast<int>(hs(11))
                      << " nwarmfb=" << static_cast<int>(hs(21))
                      << " nrefac=" << hs(22)/ncol
                      << " nreusepass=" << hs(23)/ncol
                      << " nmixfb=" << static_cast<int>(hs(24))
                      << " resid_max=" << hs(19)
                      << " resid_mean=" << hs(20)/ncol
                      << " budget_rel="
                      << ((hs(10) > 0.0) ? (hs(9) - hs(12))/hs(10) : 0.0)
                      << " telescope_rel="
                      << ((hs(15) > 0.0) ? (hs(13) - hs(14))/hs(15) : 0.0)
                      << " Ftop_impl/Ftop_sweep="
                      << ((hs(17) != 0.0) ? hs(16)/hs(17) : 0.0) << std::endl;
          }
          // problem/work_hist: CLOSE THE COLUMN-HEATING INTERVAL.  The column has just
          // applied all of its energy to u0 and the only operator left in this call is
          // the radiative momentum force.  Null hook = no-op.
          if (rt_probe != nullptr) rt_probe(2);
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
            if (layer_legacy) {
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
                // so it is paid only where the two Planck functions differ by more than
                // 4x
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
              // upward intensity is its own Planck function; the planet's internal flux
              // is
              // delivered here as an extra band-weighted source. Below the cut nothing
              // radiative is applied -- that region is optically thick and convective,
              // and
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
            } else {
              // ---- the centre-to-centre layers (see rt_layer_legacy) ----
              // Every layer now runs from one cell CENTRE to the next, split at the face
              // into its two halves, so the interval the optical depth measures and the
              // interval the source runs over are the same one.  The column is closed by
              // two half layers, the upper half of cell ie and the lower half of cell
              // icut, each with its own cell's Planck function at both ends.
              //
              // NO EXTRA TABLE LOOKUP.  The far cell's kappa rho is CARRIED in kfar from
              // the previous iteration of the sweep, so this pays exactly one ck_kappa
              // per cell and chain -- one FEWER than the staggered layers, which looked
              // the neighbour up again whenever the two Planck functions differed by 4x.
              // kfar and the running intensity Idn are the only new per-chain state.
              Real kfar[NC];
              RtF Idn[NC];
              // one half-layer step, in the chain's own precision: the same exponential
              // coefficients the staggered layers used, handed the half interval.  dsrc
              // comes back as absorbed minus emitted, which is what Src_g wants and is
              // exactly the change in the intensity.
              auto step = [&](const Real dtau, const Real mu, const Real s_in,
                              const Real s_out, RtF &I, Real &dsrc) {
                const RtF x = static_cast<RtF>(dtau/mu);
                const RtF e0 = -RT_EXPM1(-x);
                const RtF one = static_cast<RtF>(1.0);
                const RtF cin = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                               : (x/2 - x*x/3);
                const RtF cout = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                                : (x/2 - x*x/6);
                const RtF em = cin*static_cast<RtF>(s_in)
                             + cout*static_cast<RtF>(s_out);
                dsrc = static_cast<Real>(e0*I - em);
                I = (one - e0)*I + em;
              };
              for (int cc=0; cc<NC; ++cc) {
                Idn[cc] = I_down[cc][ie+1];
              }
              // down-sweep, recording the FACE intensity in between the two halves
              for (int i=ie; i>icut-1; --i) {
                const Real rho = rhoN(m,k,j,i);
                const Real dz = dx1(m,k,j,i);
                const Real drho = rho*dz;
                const Real dzf = (i < ie) ? dx1(m,k,j,i+1) : dz;
                const Real xTv = xT_g(m,k,j,i);
                const Real xPv = xP_g(m,k,j,i);
                const int iT = static_cast<int>(xTv);
                const int iP = static_cast<int>(xPv);
                const Real fT = xTv - static_cast<Real>(iT);
                const Real fP = xPv - static_cast<Real>(iP);
                for (int cc=0; cc<NC; ++cc) {
                  const int b = bandc[cc];
                  const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                                 + kc_g(m,b,i,k,j);
                  const Real kro = kap*rho;
                  const Real bown = Bb_g(m,b,i,k,j);
                  Real dsrc;
                  if (i == ie) {
                    // the top half layer: the upper half of cell ie, entered at the top
                    // face, its source held at the cell's own value -- there is nothing
                    // above it to interpolate towards
                    step(0.5*kro*dz, muc[cc], bown, bown, Idn[cc], dsrc);
                    Src_g(m,blk,i,k,j) += wfc[cc]/dz*dsrc;
                  } else {
                    // the layer between the centres of cells i+1 and i.  BFace is
                    // applied to both endpoints and symmetrically, so a radiatively
                    // inert neighbour leaves a layer emitting with its own Planck
                    // function; with equal opacities this is the straight line.
                    const Real kru = kfar[cc];
                    const Real bfar = Bb_g(m,b,i+1,k,j);
                    const Real dt_u = 0.5*kru*dzf;
                    const Real dt_l = 0.5*kro*dz;
                    const Real s_l = BFace(kru, kro, bfar, bown, bface_on);
                    const Real s_u = BFace(kro, kru, bown, bfar, bface_on);
                    const Real dtc = dt_l + dt_u;
                    const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                                 : (0.5*(s_l + s_u));
                    step(dt_u, muc[cc], s_u, s_f, Idn[cc], dsrc);
                    Src_g(m,blk,i+1,k,j) += wfc[cc]/dzf*dsrc;
                    I_down[cc][i+1] = Idn[cc];
                    step(dt_l, muc[cc], s_f, s_l, Idn[cc], dsrc);
                    Src_g(m,blk,i,k,j) += wfc[cc]/dz*dsrc;
                  }
                  kfar[cc] = kro;
                  // Direct beam, UNCHANGED: it crosses whole cells and carries no
                  // source, so the layer construction does not touch it.  Deposit the
                  // flux DIFFERENCE across the cell, not kappa rho F exp(-tau) at one
                  // face: the latter is what the grey picket fence does, and it
                  // under-deposits badly once a layer is not thin (the ratio of
                  // deposited to absorbed is u e^-u/(1 - e^-u), 0.95 at u = 0.1 but 0.58
                  // at u = 1), losing a quarter of the incident flux down a column with
                  // u ~ 0.5 -- measured against Exo-FMS on an identical column.
                  tausw[cc] += kap*drho;
                  if (lit) {
                    const Real tnew = RT_EXP(-static_cast<RtF>(tausw[cc]*facsw));
                    Qb_g(m,blk,i,k,j) += (1.0-albedo)*Fstar*ckswf(b)*wgc[cc]/facsw
                              * (transw[cc] - tnew)/dz;
                    transw[cc] = tnew;
                  }
                }
              }
              // the lower half of cell icut, from its centre down to the cut face
              {
                const Real dz = dx1(m,k,j,icut);
                for (int cc=0; cc<NC; ++cc) {
                  const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
                  Real dsrc;
                  step(0.5*kfar[cc]*dz, muc[cc], bcut, bcut, Idn[cc], dsrc);
                  Src_g(m,blk,icut,k,j) += wfc[cc]/dz*dsrc;
                  I_down[cc][icut] = Idn[cc];
                }
              }
              // Bottom of the CORRELATED-K DOMAIN, not of the column, and now AT the cut
              // face.  There the grey optical depth is of order 1e4, so the layer is
              // thermalised to e^-tau and the upward intensity is its own Planck
              // function; the planet's internal flux is delivered here as an extra
              // band-weighted source.  Below the cut nothing radiative is applied.
              RtF I_up[NC];
              for (int cc=0; cc<NC; ++cc) {
                const int b = bandc[cc];
                const Real Iint_b = (int_at_cut ? boltz_sigma/M_PI*Tint4
                                  * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
                I_up[cc] = static_cast<RtF>(Bb_g(m,b,icut,k,j) + Iint_b);
                Fb_g(m,blk,icut,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][icut]);
              }
              // up-sweep: the cut half layer first, on the opacity the down-sweep left
              // in kfar, then centre to centre.  Em_g is now the cell's OWN emission,
              // 2 (wfc/mu) kappa rho B summed over the chains, instead of the
              // centre-to-centre average of two Planck functions.
              {
                const Real dz = dx1(m,k,j,icut);
                for (int cc=0; cc<NC; ++cc) {
                  const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
                  Real dsrc;
                  step(0.5*kfar[cc]*dz, muc[cc], bcut, bcut, I_up[cc], dsrc);
                  Src_g(m,blk,icut,k,j) += wfc[cc]/dz*dsrc;
                  Em_g(m,blk,icut,k,j) += 2.0*(wfc[cc]/muc[cc])*kfar[cc]*bcut;
                }
              }
              for (int i=icut; i<ie; ++i) {
                const Real rhou = rhoN(m,k,j,i+1);
                const Real dzu = dx1(m,k,j,i+1);
                const Real dzl = dx1(m,k,j,i);
                const Real xTv = xT_g(m,k,j,i+1);
                const Real xPv = xP_g(m,k,j,i+1);
                const int iT = static_cast<int>(xTv);
                const int iP = static_cast<int>(xPv);
                const Real fT = xTv - static_cast<Real>(iT);
                const Real fP = xPv - static_cast<Real>(iP);
                for (int cc=0; cc<NC; ++cc) {
                  const int b = bandc[cc];
                  const Real kapu = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                                  + kc_g(m,b,i+1,k,j);
                  const Real kru = kapu*rhou;
                  const Real krl = kfar[cc];
                  const Real bl = Bb_g(m,b,i,k,j);
                  const Real bu = Bb_g(m,b,i+1,k,j);
                  const Real dt_l = 0.5*krl*dzl;
                  const Real dt_u = 0.5*kru*dzu;
                  const Real s_l = BFace(kru, krl, bu, bl, bface_on);
                  const Real s_u = BFace(krl, kru, bl, bu, bface_on);
                  const Real dtc = dt_l + dt_u;
                  const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc))
                                               : (0.5*(s_l + s_u));
                  Real dsrc;
                  step(dt_l, muc[cc], s_l, s_f, I_up[cc], dsrc);
                  Src_g(m,blk,i,k,j) += wfc[cc]/dzl*dsrc;
                  Fb_g(m,blk,i+1,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][i+1]);
                  step(dt_u, muc[cc], s_f, s_u, I_up[cc], dsrc);
                  Src_g(m,blk,i+1,k,j) += wfc[cc]/dzu*dsrc;
                  Em_g(m,blk,i+1,k,j) += 2.0*(wfc[cc]/muc[cc])*kru*bu;
                  kfar[cc] = kru;
                }
              }
              // the top half of cell ie, back out through the top face
              {
                const Real dz = dx1(m,k,j,ie);
                for (int cc=0; cc<NC; ++cc) {
                  const Real btop = Bb_g(m,bandc[cc],ie,k,j);
                  Real dsrc;
                  step(0.5*kfar[cc]*dz, muc[cc], btop, btop, I_up[cc], dsrc);
                  Src_g(m,blk,ie,k,j) += wfc[cc]/dz*dsrc;
                  Fb_g(m,blk,ie+1,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][ie+1]);
                }
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
        // ---- rt_implicit_column: the nearest-neighbour linearisation of this sweep ---
        // dSrc_i/dB_{i-1,i,i+1} at frozen opacity, summed over the NC quadrature chains
        // with their weights.  nblk == 1 on the grey path, so one thread owns the whole
        // column and these are ordinary stores.  Nothing here reads back into the sweep.
        const bool jac_on = implcol_;
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
        if (jac_on) {
          for (int i=is; i<ie+2; ++i) {
            jac_g(m,0,k,j,i) = 0.0;
            jac_g(m,1,k,j,i) = 0.0;
            jac_g(m,2,k,j,i) = 0.0;
          }
        }
        // d I_down(entering face of the next layer)/dB(that layer's upper neighbour),
        // and the two upward counterparts; carried along the sweeps, per chain
        Real djd[NC], djup[NC], djuo[NC];
        if (jac_on) {
          for (int cc=0; cc<NC; ++cc) {
            djd[cc] = 0.0;
            djup[cc] = 0.0;
            djuo[cc] = 0.0;
          }
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

          // top.  Everything carried from here on is the AREA-WEIGHTED intensity
          // J = A I of the dilution note above, so the incoming beam enters scaled by
          // the top face's area.
          const Real aft_ = AFC(m,k,j,ie+1);
          for (int cc=0; cc<NC; ++cc) {
            Real dtauir = gamirc[cc]*tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/muggc[cc]);
            I_ir_down_c[cc][ie+1] = (1.0-trans)*(fbc[cc]*B[ie+1]*aft_);
          }
          if (layer_legacy) {
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
              // the source function of the J equation is A B at the endpoint's own
              // position; this layer runs between the two cell centres' mean areas
              const Real bfa_g = bfr_g*ACC(m,k,j,i+1);
              const Real bwn_g = B[i]*ACC(m,k,j,i);
              for (int cc=0; cc<NC; ++cc) {
                Real dtauir = gamirc[cc]*dtau_i;
                Real x = dtauir/muggc[cc];
                Real e0 = -expm1(-x);
                Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
                Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
                // direct source: what this stream leaves in cell i, absorbed minus
                // emitted
                Src_g(m,blk,i,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i)
                                    *(e0*I_ir_down_c[cc][i+1]
                                      - fbc[cc]*(alp*bfa_g + bet*bwn_g));
                I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                                   + alp*fbc[cc]*bfa_g + bet*fbc[cc]*bwn_g;
#if RT_CACHE
                e0c[cc][i] = e0;
                alpc[cc][i] = alp;
                betc[cc][i] = bet;
#endif
              }
            }

            // bottom
            Real I_ir_up_c[NC];
            const Real afb_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              I_ir_up_c[cc] = Iint*afb_ + I_ir_down_c[cc][is];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is]/afb_;
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afb_;
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
              const Real bfa_u = bfru*ACC(m,k,j,i);
              const Real bwn_u = B[i-1]*ACC(m,k,j,i-1);
              const Real afi_ = AFC(m,k,j,i);
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
                Src_g(m,blk,i-1,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i-1)
                                      *(e0*Iup_in - fbc[cc]*(bet*bfa_u + gm*bwn_u));
                I_ir_up_c[cc] = (1.0-e0)*Iup_in
                              + bet*fbc[cc]*bfa_u + gm*fbc[cc]*bwn_u;
                Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i]/afi_;
                Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afi_;
                Fb_g(m,blk,i,k,j) += (F_ir_up_f - F_ir_down_f);
                // the cell's own emission per unit volume, both hemispheres.  dtau_i is
                // kappa rho dr for this layer, so dtau_i/dr = kappa rho and the layer
                // thickness cancels out of the volumetric rate.
                const Real dtau_ly = tau_down_r_f[i-1]-tau_down_r_f[i];
                Em_g(m,blk,i-1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_ly
                                     / dx1(m,k,j,i-1)*fbc[cc]*0.5*(bfru+B[i-1]);
              }
            }
          } else {
            // ---- the centre-to-centre layers (see rt_layer_legacy) ----
            // Same construction as the grey and correlated-k sweeps: a layer runs from
            // one cell CENTRE to the next, split at the face into two halves, each half
            // lying entirely inside one cell.  Here the optical depths come ready-made
            // from the face-cumulative tau array, so a half layer is simply half a
            // cell's dtau and nothing extra is looked up.  The coefficient cache is not
            // used on this path: both halves of a cell have the SAME thickness, but the
            // cache is indexed by the whole-cell layer the legacy branch walks.
            Real Idn[NC];
            // one half-layer step; dsrc is absorbed minus emitted, which is exactly the
            // change in the intensity, so Src_g and the divergence of Fb_g agree
            auto step = [&](const Real dtau, const Real mu, const Real s_in,
                            const Real s_out, Real &I, Real &dsrc) {
              const Real x = dtau/mu;
              const Real e0 = -expm1(-x);
              const Real cin = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              const Real cout = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              const Real em = cin*s_in + cout*s_out;
              dsrc = e0*I - em;
              I = (1.0 - e0)*I + em;
            };
            for (int cc=0; cc<NC; ++cc) {
              Idn[cc] = I_ir_down_c[cc][ie+1];
            }
            // down-sweep, recording the FACE intensity in between the two halves
            for (int i=ie; i>is-1; --i) {
              const Real dtau_l = tau_down_r_f[i]-tau_down_r_f[i+1];
              const Real dzl = dx1(m,k,j,i);
              if (i == ie) {
                // the top half layer: the upper half of cell ie, its source held at the
                // cell's own Planck function -- there is nothing above it
                // the half layer runs from the TOP FACE to the cell CENTRE, so its
                // two source endpoints take those two areas
                const Real act_ = ACC(m,k,j,ie);
                for (int cc=0; cc<NC; ++cc) {
                  const Real s = fbc[cc]*B[ie];
                  Real dsrc;
                  step(0.5*gamirc[cc]*dtau_l, muggc[cc], s*aft_, s*act_, Idn[cc], dsrc);
                  Src_g(m,blk,ie,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,ie)*dsrc;
                }
              } else {
                // the layer between the centres of cells i+1 and i.  dtau/dr is kappa
                // rho, so the two cells' optical thicknesses per unit length are the
                // emissivity weights BFace wants, applied symmetrically to both
                // endpoints.
                const Real dtau_u = tau_down_r_f[i+1]-tau_down_r_f[i+2];
                const Real dzu = dx1(m,k,j,i+1);
                const Real krl = dtau_l/dzl;
                const Real kru = dtau_u/dzu;
                const Real s_l = BFace(kru, krl, B[i+1], B[i], bface_on);
                const Real s_u = BFace(krl, kru, B[i], B[i+1], bface_on);
                const Real dtc = dtau_l + dtau_u;
                const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dtau_l/dtc))
                                             : (0.5*(s_l + s_u));
                const Real acu_ = ACC(m,k,j,i+1), afm_ = AFC(m,k,j,i+1);
                const Real acl_ = ACC(m,k,j,i);
                for (int cc=0; cc<NC; ++cc) {
                  Real dsrc;
                  step(0.5*gamirc[cc]*dtau_u, muggc[cc], fbc[cc]*s_u*acu_,
                       fbc[cc]*s_f*afm_, Idn[cc], dsrc);
                  Src_g(m,blk,i+1,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i+1)
                                        *dsrc;
                  I_ir_down_c[cc][i+1] = Idn[cc];
                  step(0.5*gamirc[cc]*dtau_l, muggc[cc], fbc[cc]*s_f*afm_,
                       fbc[cc]*s_l*acl_, Idn[cc], dsrc);
                  Src_g(m,blk,i,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i)*dsrc;
                }
              }
            }
            // the lower half of cell is, down to the bottom face
            const Real dtau_is = tau_down_r_f[is]-tau_down_r_f[is+1];
            const Real acs_ = ACC(m,k,j,is), afs_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              const Real s = fbc[cc]*B[is];
              Real dsrc;
              step(0.5*gamirc[cc]*dtau_is, muggc[cc], s*acs_, s*afs_, Idn[cc], dsrc);
              Src_g(m,blk,is,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,is)*dsrc;
              I_ir_down_c[cc][is] = Idn[cc];
            }

            // bottom, now AT the bottom face
            Real I_ir_up_c[NC];
            const Real afb_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              I_ir_up_c[cc] = Iint*afb_ + I_ir_down_c[cc][is];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is]/afb_;
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afb_;
              Fb_g(m,blk,is,k,j) += (F_ir_up_f - F_ir_down_f);
            }
            // up-sweep: the bottom half layer first, then centre to centre, accumulating
            // the band flux as it goes.  Em_g is the cell's OWN emission now, since each
            // cell owns both of its halves.
            for (int cc=0; cc<NC; ++cc) {
              const Real s = fbc[cc]*B[is];
              Real dsrc;
              step(0.5*gamirc[cc]*dtau_is, muggc[cc], s*afs_, s*acs_, I_ir_up_c[cc],
                   dsrc);
              Src_g(m,blk,is,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,is)*dsrc;
              Em_g(m,blk,is,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_is
                                  / dx1(m,k,j,is)*fbc[cc]*B[is];
            }
            for (int i=is; i<ie; ++i) {
              const Real dtau_l = tau_down_r_f[i]-tau_down_r_f[i+1];
              const Real dtau_u = tau_down_r_f[i+1]-tau_down_r_f[i+2];
              const Real dzl = dx1(m,k,j,i);
              const Real dzu = dx1(m,k,j,i+1);
              const Real krl = dtau_l/dzl;
              const Real kru = dtau_u/dzu;
              const Real s_l = BFace(kru, krl, B[i+1], B[i], bface_on);
              const Real s_u = BFace(krl, kru, B[i], B[i+1], bface_on);
              const Real dtc = dtau_l + dtau_u;
              const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dtau_l/dtc))
                                           : (0.5*(s_l + s_u));
              const Real acl2_ = ACC(m,k,j,i), afm2_ = AFC(m,k,j,i+1);
              const Real acu2_ = ACC(m,k,j,i+1);
              for (int cc=0; cc<NC; ++cc) {
                const Real pref = 2.0*M_PI*wggc[cc]*muggc[cc];
                Real dsrc;
                step(0.5*gamirc[cc]*dtau_l, muggc[cc], fbc[cc]*s_l*acl2_,
                     fbc[cc]*s_f*afm2_, I_ir_up_c[cc], dsrc);
                Src_g(m,blk,i,k,j) += pref/VLS(m,k,j,i)*dsrc;
                Fb_g(m,blk,i+1,k,j) += pref*(I_ir_up_c[cc]
                                             - I_ir_down_c[cc][i+1])/afm2_;
                step(0.5*gamirc[cc]*dtau_u, muggc[cc], fbc[cc]*s_f*afm2_,
                     fbc[cc]*s_u*acu2_, I_ir_up_c[cc], dsrc);
                Src_g(m,blk,i+1,k,j) += pref/VLS(m,k,j,i+1)*dsrc;
                Em_g(m,blk,i+1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_u
                                     / dzu*fbc[cc]*B[i+1];
              }
            }
            // the top half of cell ie, back out through the top face
            {
              const Real dtau_ie = tau_down_r_f[ie]-tau_down_r_f[ie+1];
              const Real act2_ = ACC(m,k,j,ie);
              for (int cc=0; cc<NC; ++cc) {
                const Real pref = 2.0*M_PI*wggc[cc]*muggc[cc];
                const Real s = fbc[cc]*B[ie];
                Real dsrc;
                step(0.5*gamirc[cc]*dtau_ie, muggc[cc], s*act2_, s*aft_,
                     I_ir_up_c[cc], dsrc);
                Src_g(m,blk,ie,k,j) += pref/VLS(m,k,j,ie)*dsrc;
                Fb_g(m,blk,ie+1,k,j) += pref*(I_ir_up_c[cc]
                                              - I_ir_down_c[cc][ie+1])/aft_;
              }
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
      // see rt_ali_diag.  Only the band paths carry a per-cell continuum opacity in
      // kc_g, so the diagonal is available there and nowhere else.
      const bool ali_on = rt_ali_diag && band_on;
      const int ali_nq = (ck_nq_ > 1) ? 2 : 1;
      if (rt_efix_ptr == nullptr) {
        rt_efix_ptr = new DvceArray1D<int>("rt_efix", 3);
        Kokkos::deep_copy(*rt_efix_ptr, 0);
        rt_efix_rec = new DvceArray1D<Real>("rt_efix_rec", 16);
        Kokkos::deep_copy(*rt_efix_rec, 0.0);
      }
      auto efix_g = *rt_efix_ptr;
      auto efrec_g = *rt_efix_rec;
      if (rt_desum_ptr == nullptr) {
        rt_desum_ptr = new DvceArray1D<Real>("rt_desum", 2);
      }
      auto dsum_g = *rt_desum_ptr;
      Kokkos::deep_copy(dsum_g, 0.0);
      // ---- the outer fixed-point iteration, see rt_outer_iter --------------------
      // deacc holds de_{k-1}, the total stage increment this cell has already been
      // given.  u0 carries it, so the next pass's sweep reads the partially relaxed
      // column (rt_use_cons is required for exactly that reason) and the apply below
      // recovers e^n as eiN - deacc.  With rt_outer_iter = 1 nothing here is allocated
      // and every use of deacc is compiled behind outer_on, so the path is bitwise.
      DvceArray4D<Real> deacc_g;
      DvceArray1D<Real> oconv_g;
      if (outer_on) {
        if (rt_deacc_ptr == nullptr ||
            rt_deacc_ptr->extent(0) != static_cast<size_t>(nmb1+1) ||
            rt_deacc_ptr->extent(3) != static_cast<size_t>(n1)) {
          if (rt_deacc_ptr != nullptr) delete rt_deacc_ptr;
          rt_deacc_ptr = new DvceArray4D<Real>("rt_deacc", nmb1+1, n3, n2, n1);
        }
        if (rt_oconv_ptr == nullptr) {
          rt_oconv_ptr = new DvceArray1D<Real>("rt_oconv", 6);
        }
        deacc_g = *rt_deacc_ptr;
        oconv_g = *rt_oconv_ptr;
        if (oit == 0) Kokkos::deep_copy(deacc_g, 0.0);
        Kokkos::deep_copy(oconv_g, 0.0);
      }
      const int efix_cyc = pm->ncycle;
      const bool resc_eq = rt_rescue_eq;
      // the sub-cycled local relaxation; see rt_relax_sub
      const int nsub_in = rt_relax_sub;
      const int nsub_max = (rt_relax_submax > 1) ? rt_relax_submax : 1;
      const Real xcrit = rt_relax_xcrit;
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
      if (dbg_on && outer_last) --rt_apply_debug;
      // per-cycle diagnostic: an empty View captures fine, so the lambda needs no
      // branch on the pointer itself
      const bool diag = rt_diag;
      DvceArray5D<Real> dg;
      if (diag) dg = *rt_diag_ptr;
      // ---- the radiative momentum source (problem/rt_rad_force) --------------------
      const bool radforce = rt_rad_force && grey_on;
      // problem/rt_force_center: 0 = the entry sweep's flux (Fb as it stands), 1 = the
      // converged column flux (mode 3 has already overwritten Fb with it above), 2 = the
      // average of the two, with the entry flux taken from the saved copy.
      const int fcen_a = (rt_implicit_column == 3) ? rt_force_center : 0;
      auto fbs_a = (rt_fbsave_ptr != nullptr) ? *rt_fbsave_ptr
                 : DvceArray4D<Real>("rt_fbs_dummy", 1, 1, 1, 1);
      const Real inv_c = 1.0/2.99792458e10;      // cgs: this path runs in cgs code units
      const Real arad_f = eos.tbl.arad;
      const Real xlo_f = eos.tbl.rad_lrho_lo, xhi_f = eos.tbl.rad_lrho_hi;
      const bool md_f = pm->multi_d, td_f = pm->three_d;
      const bool fverb = radforce && (rt_force_verbose > 0);
      if (fverb) --rt_force_verbose;
      // problem/rt_force_verbose normalises the cell balance by g.  rt_force_grav is a
      // CONSTANT (the pgen's GM/r_in^2), which on a global radial domain is not the
      // gravity any cell feels -- at the He star's photosphere it is 4x the local value,
      // so both the printed accelerations and the residual were mis-scaled.  Use the
      // LOCAL g(r) that the force itself works against, and keep rt_force_grav only as
      // the fallback for a mesh with no radius (the plane-parallel box).
      const Real gver = rt_force_grav;
      auto cfv_ = (rt_cf_ptr != nullptr) ? *rt_cf_ptr
                : DvceArray4D<Real>("rt_cf_dummy", 1, 1, 1, 1);
      const bool cfv_on = (rt_cf_ptr != nullptr);
      const Real apo_ = ap;             // the orbit's a; the print shadows `ap` below
      const int vcyc = pm->ncycle;
      auto eos_f = eos;
      // problem/rt_budget_verbose: the v.f work accumulator (see rt_bud_ptr)
      // problem/rt_src_dump: one column, this call only (see rt_src_dump)
      const bool sdump_ = (rt_src_dump > 0);
      if (sdump_) --rt_src_dump;
      const int sdcyc_ = pm->ncycle;
      const bool budg_ = (rt_bud_ptr != nullptr) && radforce;
      auto bud_ = budg_ ? *rt_bud_ptr : DvceArray1D<Real>("rtbuddummy", 1);
      // ---- rad_blend_use_2s: hand the conduction module this sweep's face flux -----
      // (see conduction.hpp).  The NET radial face flux of the band solver, summed over
      // its blocks, in the same code flux units the conduction operator works in.  It is
      // COPIED OUT rather than read back through rt_face_flux(), so the conduction module
      // needs no compile-time knowledge of this header, and it is written before the
      // source is applied so that the number the diffusion operator scales by w is
      // exactly the number this call scales by 1 - w.
      if (taublend && pcond_rt->rad_blend_use_2s > 0 && outer_last) {
        auto f2s_out = pcond_rt->rad_f2s;
        const int nblk_f = nblk;
        par_for("rt_f2s", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          Real ff = 0.0;
          for (int b=0; b<nblk_f; ++b) ff += Fb_g(m,b,i,k,j);
          f2s_out(m,k,j,i) = ff;
        });
        pcond_rt->rad_f2s_ready = true;
      }
      par_reduce_clip4("rt_apply", 0, nmb1, ks, ke, js, je, is, ie, nclip,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nc) {
        Real Ft = 0.0, Fb = 0.0;
        for (int b=0; b<nblk; ++b) {
          Ft += Fb_g(m,b,i+1,k,j);
          Fb += Fb_g(m,b,i,k,j);
        }
        // the face fluxes are PER UNIT AREA; every divergence below is therefore
        // (A_top F_top - A_bot F_bot)/V (the dilution note above), which under
        // rt_plane_parallel is A = 1, V = DX1 and bitwise the old expression
        const Real aft_a = AFC(m,k,j,i+1), afb_a = AFC(m,k,j,i);
        // the two-stream's share of each face in the tau blend
        Real dg_srcd = 0.0;     // problem/rt_src_dump: the raw absorbed-minus-emitted sum
        Real src;
        // THE BLEND HANDOVER IS NOT AN EMISSION, AND MUST NOT BE RELAXED.  src below is
        // the divergence of the two-stream's BLENDED share, and inside the tau ramp it
        // is dominated by the handover term d[w F]/dr: the flux the two-stream gives up
        // to the conduction operator (and the conduction operator, face by face, takes
        // up in full).  The semi-implicit step further down damps whatever it is handed
        // by (1 - e^-x)/x per cell, which is the right thing to do to the LOCAL
        // absorption-emission balance -- it is stiff -- and the wrong thing to do to a
        // handover: one half of the exchange is then damped and the other is not, the
        // two no longer cancel, and the residual is deposited in the blend layer for as
        // long as the run lasts.  MEASURED on the 1-D grey atmosphere with tau_bottom
        // = 1000 and rad_tau_lo/hi = 10/100, started FROM the analytic steady state:
        // +0.098 F of spurious heating, 68 % of it inside tau 10-100, against +0.006 F
        // with problem/rt_semi_implicit = false.
        //
        // So the source is split.  src_relax, the (1-w)-weighted local balance, goes
        // through the relaxation and is paired with the (1-w)-weighted Em it must
        // balance; src_ex = src - src_relax, which is the handover and nothing else, is
        // applied EXPLICITLY and therefore exactly.  The explicit limit of the whole
        // step is src*bdt either way.  Outside the blend (w = 0 on both faces, i.e.
        // every cell the two-stream owns alone, and every run with no tau blend at all)
        // src_relax == src and src_ex == 0 identically, so this is inert there.
        Real src_relax;
        const Real wbar = taublend ? 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1)) : 0.0;
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
          dg_srcd = src;
          src_relax = src;
          if (taublend) {
            src += (w_g(m,k,j,i+1)*(Ft*aft_a)
                    - w_g(m,k,j,i)*(Fb*afb_a))/VLA(m,k,j,i);
            src_relax *= (1.0 - wbar);
          }
        } else {
          const Real dvg = -(Ft*aft_a-Fb*afb_a)/VLA(m,k,j,i);
          if (taublend) {
            Ft *= (1.0 - w_g(m,k,j,i+1));
            Fb *= (1.0 - w_g(m,k,j,i));
          }
          src = -(Ft*aft_a-Fb*afb_a)/VLA(m,k,j,i);
          src_relax = taublend ? (1.0 - wbar)*dvg : src;
        }
        Real Qs_d = 0.0;   // the stellar heating that entered src, for the diagnostic
        if (band_on) {
          // deeper than the cut nothing radiative is applied: that region is optically
          // thick and convective, and the stellar beam died decades of optical depth
          // above
          if (i < icut_g(m,k,j)) {
            src = 0.0;
            src_relax = 0.0;
          } else {
            Real Qs = 0.0;
            for (int b=0; b<nblk; ++b) Qs += Qb_g(m,b,i,k,j);
            src += Qs;
            src_relax += Qs;
            Qs_d = Qs;
          }
        } else {
          Qs_d = Qv_g(m,k,j,i);
          src += Qs_d;
          src_relax += Qs_d;
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
        // ---- rt_implicit_column: hand the column solve the residual and the Jacobian
        // R_i is the FULL source at this state -- handover, stellar beam and all -- so
        // the outer iteration's fixed point is the exact backward-Euler balance.  The
        // Jacobian is scaled by the same (1 - wbar) the relaxed part of the source
        // carries, and zeroed below the band cut where no source is applied at all; the
        // handover's own B-dependence is left in R.
        // ---- THE TWO-LEVEL SPLIT, see rt_impl_tau_min ------------------------------
        // dtau of the cell itself, straight out of the sweep's own tau array.
        Real wthk_ = 0.0;              // this cell's share of the tridiagonal
        bool thick_ = false;
        if (implcol_) {
          // kappa_R rho dr of the cell itself, from the SAME opacity cache the sweep
          // built its layers with.  (tau_g is not filled on a box with no stellar beam,
          // so it cannot be used for this.)
          const Real dtau_c = kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*DX1(m,k,j,i);
          if (taublnd_ > 0.0) {
            if (!(dtau_c > 0.0)) {
              wthk_ = 0.0;
            } else {
              const Real u = log(dtau_c*taublnd_/taumin_)/log(taublnd_*taublnd_);
              wthk_ = (u <= 0.0) ? 0.0 : ((u >= 1.0) ? 1.0 : u);
            }
          } else {
            wthk_ = (dtau_c >= taumin_) ? 1.0 : 0.0;
          }
          thick_ = (wthk_ > 0.0);
          Real jsc = taublend ? (1.0 - wbar) : 1.0;
          if (band_on && i < icut_g(m,k,j)) jsc = 0.0;
          jsc *= wthk_;             // a thin row carries NO two-stream term at all
          if (jsc != 1.0) {
            jac_g(m,0,k,j,i) *= jsc;
            jac_g(m,1,k,j,i) *= jsc;
            jac_g(m,2,k,j,i) *= jsc;
          }
          const Real tkc = T_g(m,k,j,i);
          // dbdt doubles as the "how much of this cell is an unknown" weight: zero here
          // removes the row's own diagonal AND a thick neighbour's coupling to it
          dbt_g(m,k,j,i) = (tkc > 0.0) ? wthk_*4.0*boltz_sigma/M_PI*tkc*tkc*tkc : 0.0;
          res_g(m,k,j,i) = wthk_*src;
          dtx_g(m,k,j,i) = 0.0;     // filled below by the share the cell relaxes itself
        }
        // the share the tridiagonal owns is not applied here; the rest takes the old
        // nonlinear relaxation, unchanged
        const bool skip_de = (implcol_ && (wthk_ >= 1.0)) || mode3_;
        // ---- ONE-SHOT ASSEMBLY DUMP (problem/rt_outer_verbose, cycle 0, one column) --
        // A and E are taken with EXACTLY the pre-existing apply block's definitions:
        // E = Em, the per-volume emission the sweep subtracted, and A = src_relax + Em,
        // so R = src = A - E + src_ex with the handover counted once and only once.
        if (implcol_ && rtdbg_ && m == 0 && k == ks && j == js) {
          Real Emd = 0.0;
          for (int b=0; b<nblk; ++b) Emd += Em_g(m,b,i,k,j);
          if (taublend) Emd *= 1.0 - wbar;
          if (band_on && i < icut_g(m,k,j)) Emd = 0.0;
          const Real dtau_c = kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*DX1(m,k,j,i);
          const Real tkd = T_g(m,k,j,i);
          const Real dbd = (tkd > 0.0) ? 4.0*boltz_sigma/M_PI*tkd*tkd*tkd : 0.0;
          Kokkos::printf("### rtcol_asm i=%d dtau=%.4e thick=%d T=%.5e rho=%.4e "
                         "A=%.6e E=%.6e R=%.6e src_ex=%.4e jm=%.6e j0=%.6e jp=%.6e "
                         "dBdT=%.6e sink=%.6e A/E=%.4e\n",
                         i, dtau_c, thick_ ? 1 : 0, tkd, rhoN(m,k,j,i),
                         src_relax + Emd, Emd, src, src - src_relax,
                         jac_g(m,0,k,j,i), jac_g(m,1,k,j,i), jac_g(m,2,k,j,i), dbd,
                         4.0*kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*boltz_sigma*tkd*tkd*tkd,
                         (Emd != 0.0) ? (src_relax + Emd)/Emd : 0.0);
        }
        // ---- rt_ali_diag: 1 - Lambda*_ii, the local escape factor -----------------
        Real ome_ = 1.0;
        if (ali_on) {
          const Real dtau_a = kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*DX1(m,k,j,i);
          if (dtau_a > 0.0) {
            Real lst = 0.0;
            if (ali_nq > 1) {
              const Real mq[2] = {0.21132487, 0.78867513};
              for (int q=0; q<2; ++q) {
                const Real xq = dtau_a/mq[q];
                lst += 0.5*((xq > 1.0e-3) ? (1.0 - (-expm1(-xq))/xq)
                                          : (0.5*xq - xq*xq/6.0));
              }
            } else {
              const Real xq = dtau_a*1.66;
              lst = (xq > 1.0e-3) ? (1.0 - (-expm1(-xq))/xq) : (0.5*xq - xq*xq/6.0);
            }
            ome_ = 1.0 - lst;
            if (!(ome_ > 1.0e-8)) ome_ = 1.0e-8;   // never divide by zero below
            if (ome_ > 1.0) ome_ = 1.0;
          }
        }
        Real de = src*bdt;
        // rt_cell_report bookkeeping: the pieces of the step, kept for the report below
        Real dg_A = 0.0, dg_Em = 0.0, dg_deq = 0.0;
        Real dg_it[8];
        int dg_nit = 0;
        bool dg_resc = false;
        int dg_nsub = 1;   // rt_src_dump: the sub-cycle count actually used
        for (int q=0; q<8; ++q) dg_it[q] = 0.0;
        // problem/rt_explicit (or problem/rt_semi_implicit = false): nothing else
        // touches de.  See rt_explicit and rt_semi_implicit.
        if (!explicit_on && semi_imp && !skip_de) {
          Real Em = 0.0;
          for (int b=0; b<nblk; ++b) Em += Em_g(m,b,i,k,j);
          if (taublend) {
            Em *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (band_on && i < icut_g(m,k,j)) Em = 0.0;
          // e^n, the energy this stage STARTED from.  Without the outer iteration that
          // is just the state the solver reads; with it, u0 already carries de_{k-1}
          // and has to be walked back, because de below is the TOTAL stage increment
          // and not an addition to what the previous pass left.  T_g and Em, by
          // contrast, are deliberately the UPDATED state's: Em(T^{k-1}) together with
          // t0 = T^{k-1} is what makes the Newton's Em*(T(e^n+de)/t0)^4 the emission at
          // the NEW energy, i.e. F(de) = de - bdt(A - Em(e^n+de)) exactly.
          const Real ei = outer_on ? (eiN(m,k,j,i) - deacc_g(m,k,j,i)) : eiN(m,k,j,i);
          if (Em > 0.0 && ei > 0.0) {
            const Real src_ex = src - src_relax;     // the handover, applied exactly
            const Real sdt = src_relax*bdt;
            if (semilin) {
              const Real lam = ome_*4.0*Em/ei;     // see rt_ali_diag
              const Real x = lam*bdt;
              de = (x > 1.0e-4) ? (src_relax/lam)*(-expm1(-x)) : sdt;
              de += src_ex*bdt;
            } else {
              // SUB-CYCLED RELAXATION (problem/rt_relax_sub).  The step above relaxes
              // toward the equilibrium of the UN-relaxed column: A is one Jacobi pass,
              // formed from the intensities of the state the sweep saw, and the cell is
              // then damped toward it by a per-cell factor (1 - exp(-x))/x.  That factor
              // is what breaks the problem: the non-local exchange between two cells
              // cancels exactly only while both are damped by the SAME amount, and two
              // neighbours at different x are not.  What survives the near-cancellation
              // is a pressure perturbation proportional to dt times the local compression
              // rate, i.e. a term that pumps rather than damps, and in a closed box it
              // feeds the organ-pipe modes: the He-star box's saturated v_rms came out
              // LINEAR in dt, which no converged solution can be.
              //
              // The fix is to stop taking one large damped step.  A, the absorbed field,
              // stays frozen for the stage -- it is what the sweep computed and re-doing
              // the sweep is what this is trying to avoid -- but the LOCAL balance is
              // integrated in nsub sub-steps of bdt/nsub, re-forming the cell's own
              // emission from the running energy after each one (Em(T) ~ T^4, with T(e)
              // from the EOS where there is one).  Each sub-step is then taken at small
              // x, where (1 - exp(-x))/x -> 1 and the damping that broke the cancellation
              // is gone, while the sum still lands on the same fixed point E(T) = A.
              //
              // nsub is chosen per cell from the stiffness of the first sub-step,
              // ceil(x/rt_relax_xcrit), floored at rt_relax_sub and capped at
              // rt_relax_submax.  rt_relax_sub <= 1 disables the whole thing and is
              // BITWISE the single-step form above.
              // A_eff: the equilibrium the cell is relaxed toward is the one where the
              // NET exchange vanishes, (1 - Lambda*)(E(T) - E) = src, i.e. E_eq = E +
              // src/(1 - Lambda*).  With rt_ali_diag off ome_ is 1 and this is the old
              // A = src + Em, bit for bit.  See rt_ali_diag.
              const Real absn = src_relax/ome_ + Em;  // A, held fixed over the step
              // e_eq - e.  With nothing arriving the equilibrium is T = 0, i.e. -e.
              const Real deq0 = (absn > 0.0) ? ei*(sqrt(sqrt(absn/Em)) - 1.0) : -ei;
              int nsub = 1;
              if (nsub_in > 1) {
                nsub = nsub_in;
                if (deq0 != 0.0 && xcrit > 0.0) {
                  const Real xn = ceil((sdt/deq0)/xcrit);
                  if (xn > static_cast<Real>(nsub)) {
                    nsub = (xn >= static_cast<Real>(nsub_max)) ? nsub_max
                                                              : static_cast<int>(xn);
                  }
                }
                if (nsub > nsub_max) nsub = nsub_max;
                if (nsub < 1) nsub = 1;
              }
              dg_nsub = nsub;
              const Real sub_bdt = bdt/static_cast<Real>(nsub);
              const Real sub_sdt = src_relax*sub_bdt;
              const Real t0 = T_g(m,k,j,i);
              const Real d0 = rhoN(m,k,j,i);
              Real ec = ei;                            // the running internal energy
              Real Emc = Em;                           // its emission, re-formed below
              de = 0.0;
              for (int s=0; s<nsub; ++s) {
                const Real deq = (absn > 0.0) ? ec*(sqrt(sqrt(absn/Emc)) - 1.0) : -ec;
                dg_A = absn; dg_Em = Emc; dg_deq = deq;
                Real des;
                if (deq != 0.0) {
                  const Real x = sub_sdt/deq;          // >= 0: src and deq share a sign
                  des = (x > 1.0e-4) ? deq*(-expm1(-x)) : sub_sdt;
                } else {
                  des = sub_sdt;
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
                // Under sub-cycling T_old and E(T_old) stay the values of the cell at the
                // START of the stage, so r4 is measured against the same fixed A.
                if (newton_on && eos.IsGeneral()) {
                  // both terms carry the same (1 - Lambda*), so F(de) = de -
                  // bdt*src_relax + ome*bdt*Em*(r4 - 1): explicit at r4 = 1, zero at the
                  // ALI equilibrium E(T) = A_eff.  ome_ = 1 is the old expression.
                  const Real abdt = ome_*absn*sub_bdt, embdt = ome_*Em*sub_bdt;
                  // SAFEGUARDED, because a bare Newton here does not converge.  F is
                  // MONOTONE INCREASING in de -- both de and E(T(e+de)) rise with de --
                  // so its root is unique and bracketing is available for free: every
                  // iterate with F < 0 is a lower bound on the root and every one with F
                  // > 0 an upper bound.  Without that, in an optically thin
                  // radiation-dominated cell the step is very nearly the EXPLICIT one
                  // (the damping denominator dfx -> 1 as the cell's own emission stops
                  // controlling its temperature) and it overshoots the root by orders of
                  // magnitude.  MEASURED in RG_v4/out.txt, where rt_cell_report printed
                  // the iterates: on a cell with e = 6.02e-3 the sequence was de =
                  // -1.63e+1, +1.06e+1, -3.33e+1 -- oscillating, each iterate 3-4 decades
                  // past e -- and on another, e = 3.06e-1 was handed de = -3.37e-1 on the
                  // FIRST step.  Both left e + de <= 0 and fell through to the rescue,
                  // which can only leave the cell at ~1e-3 e, from which the table EOS
                  // returns a floor temperature and the sound speed collapses the
                  // timestep.  That is the mechanism behind the He-star top-cell "dt
                  // COLLAPSE" events and the red giant's 152 eos_tclamp family.
                  //
                  // The lower bracket nlo is the 99.9 % floor the rescue used, and F(nlo)
                  // < 0 whenever the absorption A is non-negative, so it is a true lower
                  // bound on the root.  ONLY a step that would cross zero energy is
                  // replaced, by a bisection of [nlo, des]; every other step is taken
                  // exactly as before, so this is inert on every cell the bare iteration
                  // handled.  A wider safeguard was tried first -- bracketing from both
                  // sides and clamping the converged answer into the bracket -- and is
                  // NOT in the tree: it perturbed healthy cells (the 1-D He-star column's
                  // timestep fell 0.107 -> 0.034 s with no rescue anywhere in the run),
                  // because the Newton legitimately walks past the e ~ T equilibrium
                  // wherever the EOS is far from that, which is the reason the refinement
                  // exists at all.
                  const Real nlo = -(1.0 - 1.0e-3)*ec;
                  if (t0 > 0.0 && d0 > 0.0) {
                    for (int it=0; it<8; ++it) {
                      const Real e1 = ec + des;
                      if (!(e1 > 0.0)) break;
                      const Real tc = eos.Temperature(d0, e1);
                      const Real t1 = tc*eos.temp_cgs;
                      if (!(t1 > 0.0)) break;
                      const Real cv = d0*eos.SpecificHeatCv(d0, e1, tc);
                      if (!(cv > 0.0)) break;
                      const Real r4 = SQR(SQR(t1/t0));
                      const Real fx = des - abdt + embdt*r4;
                      // d(r4)/de = 4 r4/T dT/de, with dT/de = temp_cgs/(d c_v)
                      const Real dfx = 1.0 + embdt*4.0*r4/t1*(eos.temp_cgs/cv);
                      if (!(dfx > 0.0)) break;
                      const Real step = fx/dfx;
                      const Real dn = des - step;
                      // THE ONE INTERVENTION.  A step that would leave a non-positive
                      // internal energy is replaced by a bisection of [nlo, des], which
                      // brackets the root: F is monotone and F(des) > 0 is what makes the
                      // step negative in the first place, while F(nlo) < 0 for any
                      // non-negative absorption.  Every other step is taken EXACTLY as
                      // before -- dn is `des - step`, the same expression `des -= step`
                      // evaluated -- so a cell the bare iteration handled is untouched,
                      // bit for bit, and `ec + des > 1e-3 ec > 0` now holds
                      // unconditionally after the loop.
                      des = (ec + dn > 0.0) ? dn : 0.5*(nlo + des);
                      if (dg_nit < 8) dg_it[dg_nit++] = des;
                      if (fabs(step) <= 1.0e-8*(fabs(des) + fabs(ec))) break;
                    }
                  }
                  // THE LOOP CHECKS e1 > 0 AT THE TOP, NOT AT THE BOTTOM.  The last `des
                  // -= step` is never validated, so Newton can exit having pushed the
                  // cell to ec + des <= 0 -- a one-step NaN with no counterpart in the
                  // closed form, which guarantees e1 = ec*exp(-x) > 0.  R9's all-column
                  // NaN at t = 1.99e5 is under bisection; this closes the only path in
                  // this branch to a non-positive energy.  Falling back to a 99.9 % drop
                  // keeps the cell cooling hard without ever crossing zero.
                  if (!(ec + des > 0.0)) {
                    const Real de_nt = des;        // what Newton left, for the report
                    // problem/rt_rescue_eq: land on the equilibrium the cell actually
                    // sees rather than on a fixed 99.9 % drop.  deq is e_eq - e from the
                    // closed form above, i.e. Em(T_eq) = A = src + Em; the old floor is
                    // kept underneath it, so this can only make the rescue gentler.
                    const Real defl = -(1.0 - 1.0e-3)*ec;
                    if (resc_eq && deq < 0.0 && deq > defl) {
                      des = deq;
                      Kokkos::atomic_fetch_add(&efix_g(1), 1);
                    } else {
                      des = defl;
                      Kokkos::atomic_fetch_add(&efix_g(2), 1);
                    }
                    if (Kokkos::atomic_fetch_add(&efix_g(0), 1) == 0) {
                      efrec_g(0) = static_cast<Real>(m);
                      efrec_g(1) = static_cast<Real>(k);
                      efrec_g(2) = static_cast<Real>(j);
                      efrec_g(3) = static_cast<Real>(i);
                      efrec_g(4) = rhoN(m,k,j,i);
                      efrec_g(5) = ei;
                      efrec_g(6) = t0;
                      efrec_g(7) = src;
                      efrec_g(8) = src_relax;
                      efrec_g(9) = Em;
                      efrec_g(10) = absn;
                      efrec_g(11) = deq;
                      efrec_g(12) = de_nt;
                      efrec_g(13) = bdt;
                      efrec_g(14) = des;
                      efrec_g(15) = static_cast<Real>(efix_cyc);
                    }
                    dg_resc = true;
                  }
                }
                de += des;
                ec = ei + de;
                // Re-form the emission of the UPDATED state for the next sub-step.  Em
                // scales as T^4 and T comes from the EOS wherever the Newton block has
                // one; with e ~ T (ideal gas, and the closed form's own assumption) the
                // ratio is just ec/ei.  This is the whole point of the sub-cycle: without
                // it every sub-step would relax toward the same stale equilibrium.
                if (s + 1 < nsub) {
                  if (!(ec > 0.0)) break;
                  Real trat;
                  if (newton_on && eos.IsGeneral() && t0 > 0.0 && d0 > 0.0) {
                    const Real t1 = eos.Temperature(d0, ec)*eos.temp_cgs;
                    trat = (t1 > 0.0) ? t1/t0 : 0.0;
                  } else {
                    trat = ec/ei;
                  }
                  Emc = Em*SQR(SQR(trat));
                  if (!(Emc > 0.0)) break;
                }
              }
              // ...and the handover the relaxation was deliberately not shown.  Added
              // after the Newton refinement and the rescue, both of which are statements
              // about the LOCAL balance alone; src_ex is zero unless this cell sits
              // inside the tau ramp.  It is OUTSIDE the sub-cycle: it is an exact
              // explicit term, not part of the local balance being relaxed.
              de += src_ex*bdt;
            }
          }
        }
        const Real de_pre = de;
        if (demax > 0.0 && !skip_de) {
          const Real dl = LimitRTSource(de,
              outer_on ? (eiN(m,k,j,i) - deacc_g(m,k,j,i)) : eiN(m,k,j,i), demax);
          if (dl != de) { ++nc; de = dl; }
        }
        if (report_on && outer_last) {
          const Real dxb = DX1(m,k,j,i);
          Kokkos::atomic_add(&dsum_g(0), src*bdt*dxb);
          Kokkos::atomic_add(&dsum_g(1), de*dxb);
        }
        if (diag && outer_last) {
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
        if (dbg_on && outer_last && m == dbg_m && k == dbg_k && j == dbg_j &&
            i > ie - dbg_n) {
          Real Em = 0.0, Qs = 0.0;
          for (int b=0; b<nblk; ++b) { Em += Em_g(m,b,i,k,j); Qs += Qb_g(m,b,i,k,j); }
          // the raw direct source, before the blend handover and the beam are added
          Real srcraw = 0.0;
          for (int b=0; b<nblk; ++b) srcraw += Src_g(m,b,i,k,j);
          const Real ei = eiN(m,k,j,i);
          Kokkos::printf("rt_apply i=%d T=%.4e d=%.4e e=%.4e Fb=%.6e Ft=%.6e "
                         "divF=%.14e srcraw=%.14e src=%.4e Qs=%.4e Em=%.4e "
                         "lamdt=%.4e de=%.4e "
                         "de/e=%.4e tau=%.4e w=%.4e dx=%.4e "
                         "Ab=%.6e At=%.6e Ac=%.6e V=%.6e\n",
                         i, T_g(m,k,j,i), rhoN(m,k,j,i), ei, Fb, Ft,
                         -(Ft*aft_a-Fb*afb_a)/VLA(m,k,j,i), srcraw, src, Qs, Em,
                         (Em > 0.0 && ei > 0.0) ? 4.0*Em/ei*bdt : 0.0,
                         de, de/ei, taublend ? tauf_g(m,k,j,i) : 0.0,
                         taublend ? w_g(m,k,j,i) : 0.0, DX1(m,k,j,i),
                         afb_a, aft_a, ACC(m,k,j,i), VLA(m,k,j,i));
        }
        // ---- rt_cell_report ---------------------------------------------------
        if (report_on && outer_last) {
          const bool fixedcell = fixed_on && (k == rep_k) && (j == rep_j) &&
              (fabs(X1V(m,i) - rep_r) < 0.5*DX1(m,k,j,i));
          bool doprint = fixedcell;
          if (dg_resc) {
            if (Kokkos::atomic_fetch_add(&repc_g(0), 1) == 0) doprint = true;
          }
          if (doprint && i >= icut_g(m,k,j)) {
            const Real kap = kc_g(m,0,i,k,j);
            const Real rho = rhoN(m,k,j,i);
            const Real dtc = kap*rho*DX1(m,k,j,i);
            Real tautop = 0.0;
            for (int i2=i; i2<=ie; ++i2) {
              tautop += kc_g(m,0,i2,k,j)*rhoN(m,k,j,i2)*DX1(m,k,j,i2);
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
            const Real absdn = wf/DX1(m,k,j,i)*e0*idn;
            const Real absup = wf/DX1(m,k,j,i)*e0*iup;
            const Real emidn = wf/DX1(m,k,j,i)*(alp*bip + bet*bi);
            const Real emiup = wf/DX1(m,k,j,i)*(bet*bip + gmq*bi);
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
              dg_resc ? "RESCUE" : "FIXED", rep_cyc, rep_time, m, k, j, i, X1V(m,i),
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
        // the part of de this cell applies ITSELF: all of it when the tridiagonal owns
        // none of the cell, none when it owns all
        const Real de_app = implcol_ ? (1.0 - wthk_)*de : de;
        if (sdump_ && m == 0 && k == ks && j == js &&
            (!band_on || i >= icut_g(m,k,j))) {
          const Real dxc = DX1(m,k,j,i);
          const Real wb_ = taublend ? w_g(m,k,j,i) : 0.0;
          const Real wt_ = taublend ? w_g(m,k,j,i+1) : 0.0;
          const Real dtauc = kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*dxc;
          Kokkos::printf("### rt_srcdump cyc=%d i=%d w_b=%.6f w_t=%.6f dtau=%.4e "
                         "T=%.6e dx=%.6e Fb=%.8e Ft=%.8e divf=%.8e divw=%.8e "
                         "srcdir=%.8e src=%.8e src_relax=%.8e de_o_dt=%.8e "
                         "A=%.8e Em=%.8e AoE_m1=%.8e deq=%.8e x=%.8e nsub=%d "
                         "nit=%d resc=%d clip=%d ei=%.8e de=%.8e\n",
                         sdcyc_, i, wb_, wt_, dtauc, T_g(m,k,j,i), dxc, Fb, Ft,
                         -(Ft - Fb)/dxc,
                         -((1.0 - wt_)*Ft - (1.0 - wb_)*Fb)/dxc,
                         dg_srcd, src, src_relax, de/bdt,
                         dg_A, dg_Em, (dg_Em != 0.0) ? (dg_A/dg_Em - 1.0) : 0.0,
                         dg_deq, (dg_deq != 0.0) ? (src_relax*bdt/dg_deq) : 0.0,
                         dg_nsub, dg_nit, dg_resc ? 1 : 0,
                         (de != de_pre) ? 1 : 0, eiN(m,k,j,i), de);
        }
        if (implcol_ && wthk_ < 1.0) {
          // Record the change in this cell's Planck function that its own relaxation is
          // applying, so a neighbouring row can take the exchange with the part of it
          // that is NOT an unknown as a known right-hand-side term (see rt_col_dtex).
          //
          // IT IS MEASURED AGAINST THE STATE THE SWEEP SAW, not against e^n: R was formed
          // from the running column, so the only thing the rows have not been shown is
          // the step about to be taken.  de is the TOTAL stage increment (rt_outer_iter
          // re-derives it from e^n every pass), hence the walk back through deacc.
          const Real t0x = T_g(m,k,j,i);
          const Real ei0 = outer_on ? (eiN(m,k,j,i) - deacc_g(m,k,j,i)) : eiN(m,k,j,i);
          if (t0x > 0.0 && ei0 > 0.0 && (ei0 + de_app) > 0.0) {
            Real t1x;
            if (newton_on && eos.IsGeneral()) {
              t1x = eos.Temperature(rhoN(m,k,j,i), ei0 + de_app)*eos.temp_cgs;
            } else {
              t1x = t0x*(ei0 + de_app)/eiN(m,k,j,i);
            }
            if (t1x > 0.0) {
              dtx_g(m,k,j,i) = boltz_sigma/M_PI*(SQR(SQR(t1x)) - SQR(SQR(t0x)));
              // AND RE-ANCHOR T^n.  The tridiagonal's heat-capacity term measures
              // T^{k-1} - T^n to keep the conduction operator from being applied in full
              // on every outer pass.  This cell's OWN relaxation also moved T^{k-1}, and
              // that part is not an unconverged conduction increment: without moving the
              // anchor with it the next pass's conduction row would try to undo it.  THIS
              // IS THE OUTER-LOOP DEFECT that made the fixed point diverge.
              if (outer_on && oit > 0 && eos.temp_cgs > 0.0) {
                tn_g(m,k,j,i) += (t1x - t0x)/eos.temp_cgs;
              }
            }
          }
        }
        if (skip_de) {
          // NOTHING is applied here: the tridiagonal owns the energy update, and it is
          // called by the wrapper as soon as this kernel is done.  de is kept only for
          // the diagnostics above, where it is the EXPLICIT rate this state would give.
          if (outer_on) {
            const Real dprev = deacc_g(m,k,j,i);
            deacc_g(m,k,j,i) = de;
            const Real ade = fabs(de);
            Kokkos::atomic_max(&oconv_g(0), (ade > 0.0) ? fabs(de - dprev)/ade : 0.0);
            Kokkos::atomic_max(&oconv_g(1), ade);
          }
        } else if (outer_on) {
          // de REPLACES de_{k-1}: u0 carries the running total, never a sum of passes
          const Real dprev = deacc_g(m,k,j,i);
          u0(m,IEN,k,j,i) += de_app - dprev;
          deacc_g(m,k,j,i) = de_app;
          const Real ade = fabs(de_app);
          Kokkos::atomic_max(&oconv_g(0), (ade > 0.0) ? fabs(de_app - dprev)/ade : 0.0);
          Kokkos::atomic_max(&oconv_g(1), ade);
          // ABSOLUTE per class: the ratio above is a max over cells and is therefore
          // owned by whichever cell has the SMALLEST |de_k|, which says nothing about
          // convergence.  Slots 2/3 are the cells this kernel applies (thin, or the
          // relaxed share of a blended one), 4/5 the ones that touch a tridiagonal row.
          const int cls = (implcol_ && thick_) ? 4 : 2;
          Kokkos::atomic_max(&oconv_g(cls), fabs(de_app - dprev));
          Kokkos::atomic_max(&oconv_g(cls+1), ade);
        } else {
          u0(m,IEN,k,j,i) += de_app;
        }
        // ---- radiative momentum source, see rt_rad_force --------------------------
        if (radforce && outer_last) {
          const Real rho = rhoN(m,k,j,i);
          Real wr, dwdx;
          rad_taper::Weight(log10(rho), xlo_f, xhi_f, wr, dwdx);
          // the cell's net two-stream flux, positive upward (rt_force_center: Ft/Fb
          // are this call's Fb array, so 0/1 differ only in what the column wrote; 2
          // averages it with the saved entry flux)
          Real ftc_ = Ft, fbc_ = Fb;
          if (fcen_a == 2) {
            ftc_ = 0.5*(Ft + fbs_a(m,k,j,i+1));
            fbc_ = 0.5*(Fb + fbs_a(m,k,j,i));
          }
          const Real fnet = 0.5*(ftc_ + fbc_);
          Real f1 = (1.0 - wr)*rho*kc_g(m,0,i,k,j)*fnet*inv_c;
          Real f2 = 0.0, f3 = 0.0, prgw = 0.0;
          if (dwdx != 0.0) {
            // Prad grad w, with grad w = w'(x) grad rho/(rho ln10) -- the SAME derivative
            // the EOS put into chi_rho, so the two cancel in the continuum limit
            const Real tk = T_g(m,k,j,i);
            const Real cg = (arad_f*tk*tk*tk*tk/3.0)*dwdx*M_LOG10E/rho;
            prgw = cg*(rhoN(m,k,j,i+1) - rhoN(m,k,j,i-1))/(X1V(m,i+1) - X1V(m,i-1));
            f1 += prgw;
            if (md_f) {
              f2 = cg*(rhoN(m,k,j+1,i) - rhoN(m,k,j-1,i))/(2.0*size.d_view(m).dx2);
            }
            if (td_f) {
              f3 = cg*(rhoN(m,k+1,j,i) - rhoN(m,k-1,j,i))/(2.0*size.d_view(m).dx3);
            }
          }
          const Real dinv = 1.0/u0(m,IDN,k,j,i);
          const Real v1 = u0(m,IM1,k,j,i)*dinv;
          const Real v2 = u0(m,IM2,k,j,i)*dinv;
          const Real v3 = u0(m,IM3,k,j,i)*dinv;
          u0(m,IM1,k,j,i) += f1*bdt;
          if (md_f) u0(m,IM2,k,j,i) += f2*bdt;
          if (td_f) u0(m,IM3,k,j,i) += f3*bdt;
          u0(m,IEN,k,j,i) += (v1*f1 + v2*f2 + v3*f3)*bdt;
          if (budg_) {
            const Real dvb = DX1(m,k,j,i)*size.d_view(m).dx2*size.d_view(m).dx3;
            Kokkos::atomic_add(&bud_(10), (v1*f1 + v2*f2 + v3*f3)*bdt*dvb);
          }
          if (fverb && wr > 0.0 && wr < 1.0) {
            // the hydrostatic balance of this cell: pressure gradient, gravity, source
            EOSThermoState sm, sp;
            eos_f.tbl.EvalNoMu(rhoN(m,k,j,i-1), T_g(m,k,j,i-1), sm);
            eos_f.tbl.EvalNoMu(rhoN(m,k,j,i+1), T_g(m,k,j,i+1), sp);
            const Real ap = -(sp.p - sm.p)/((X1V(m,i+1) - X1V(m,i-1))*rho);
            const Real gloc = !pp_
                ? EffGravAt(grav, apo_, X1V(m,i), grav_pmass, omega,
                            cfv_on ? cfv_(m,k,j,3) : 0.0, tide)
                : gver;
            const Real gg = (gloc > 0.0) ? gloc : ((gver > 0.0) ? gver : 1.0);
            Kokkos::printf("### rt_force ncycle=%d i=%d rho=%.4e T=%.4e w=%.4f "
                           "a_p/g=%.6e a_f/g=%.6e a_prgw/g=%.6e resid/g=%.6e "
                           "F=%.4e kap=%.4e\n",
                           vcyc, i, rho, T_g(m,k,j,i), wr, ap/gg, f1/(rho*gg),
                           prgw/(rho*gg), (ap - gg + f1/rho)/gg,
                           0.5*(Ft + Fb), kc_g(m,0,i,k,j));
          }
        }
      });
      rt_nclip_last = nclip;
      if (outer_last) RTSourceLimiterWarn(nclip);
      // ---- the per-cycle clip/rescue CENSUS (problem/rt_outer_verbose) -------------
      // How many cells the per-STEP caps actually touched this stage: nclip is the
      // LimitRTSource count reduced over the apply kernel, efix(0..2) the running Newton
      // positivity rescue totals.  Diagnostic only; nothing here changes an answer.
      if (rt_outer_verbose && global_variable::my_rank == 0 && pm->ncycle % 200 == 0) {
        auto hcen = Kokkos::create_mirror_view(efix_g);
        Kokkos::deep_copy(hcen, efix_g);
        std::cout << "### rt_clip ncycle=" << pm->ncycle << " t=" << pm->time
                  << " nclip=" << nclip << " efix_tot=" << hcen(0)
                  << " resc_eq=" << hcen(1) << " resc_floor=" << hcen(2) << std::endl;
      }
      // ---- the fixed-point convergence of this pass, see rt_outer_iter ------------
      if (outer_on && (rt_outer_verbose || report_on) && rt_report_every > 0 &&
          (pm->ncycle % rt_report_every == 0) && global_variable::my_rank == 0) {
        auto hc = Kokkos::create_mirror_view(oconv_g);
        Kokkos::deep_copy(hc, oconv_g);
        std::cout << "### rt_outer ncycle=" << pm->ncycle << " pass " << (oit+1)
                  << "/" << nit << "  max|de_k-de_k-1|/|de_k| = " << hc(0)
                  << "  max|de_k| = " << hc(1)
                  << "  | relaxed: max|dde| = " << hc(2) << " max|de| = " << hc(3)
                  << "  | interface: max|dde| = " << hc(4) << " max|de| = " << hc(5)
                  << std::endl;
      }
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
          auto hr = Kokkos::create_mirror_view(efrec_g);
          Kokkos::deep_copy(hr, efrec_g);
          std::cout << "    first such cell: cycle " << static_cast<int>(hr(15))
                    << " (m,k,j,i) = (" << static_cast<int>(hr(0)) << ","
                    << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
                    << static_cast<int>(hr(3)) << ")  rho = " << hr(4)
                    << "  e = " << hr(5) << "  T = " << hr(6) << " K" << std::endl
                    << "      src = " << hr(7) << "  src_relax = " << hr(8)
                    << "  Em = " << hr(9) << "  A = src_relax+Em = " << hr(10)
                    << "  bdt = " << hr(13) << std::endl
                    << "      src*bdt/e = " << (hr(5) != 0.0 ? hr(7)*hr(13)/hr(5) : 0.0)
                    << "  Em*bdt/e = " << (hr(5) != 0.0 ? hr(9)*hr(13)/hr(5) : 0.0)
                    << "  deq/e = " << (hr(5) != 0.0 ? hr(11)/hr(5) : 0.0)
                    << "  de_newton/e = " << (hr(5) != 0.0 ? hr(12)/hr(5) : 0.0)
                    << "  de_used/e = " << (hr(5) != 0.0 ? hr(14)/hr(5) : 0.0)
                    << std::endl;
          efix_warned = true;
        }
        efix_seen = he(0);
      }
      // the source's energy budget: what the flux divergence asked for against what the
      // (sub-cycled, relaxed, limited) application actually deposited.  See rt_desum_ptr.
      if (report_on && fixed_on && outer_last) {
        auto hd = Kokkos::create_mirror_view(dsum_g);
        Kokkos::deep_copy(hd, dsum_g);
        std::cout << "### rt_desum ncycle=" << pm->ncycle << " sum(src*dt*dx) = " << hd(0)
                  << "  sum(de*dx) = " << hd(1) << "  rel = "
                  << ((hd(0) != 0.0) ? (hd(1) - hd(0))/hd(0) : 0.0) << std::endl;
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
          col(i,4) = X1F(md,i);
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
