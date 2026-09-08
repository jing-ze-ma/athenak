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
// problem/rt_de_max: the cap in LimitRTSource, as a fraction of the cell's internal
// energy per RT application. Applies to every EXPLICIT radiative update -- grey and
// correlated-k, split and monolithic. Set <= 0 to disable the limiter entirely.
inline Real rt_de_max = 0.5;
// problem/ck_int_at_cut: deliver the planet's internal flux sigma T_int^4 as an extra
// upward source at the correlated-k cut (the historical behaviour, true). Set false when
// the layers below the cut carry it themselves -- <mhd|hydro>/isotropic_conduction =
// radiative with rad_flux_inner at the bottom wall -- so the cut's upward intensity is
// just the thermalised Planck function and nothing is counted twice.
inline bool rt_int_at_cut = true;
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
inline bool rt_srclim_warned = false;             // the one-time warning has been issued

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
      const int nblk = (nchain_rt + NC - 1)/NC;
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
        if (rt_ck) {
          rt_kc_ptr = new DvceArray5D<Real>("rt_kc", nmb, CK_NB, n1, n3, n2);
          rt_Bb_ptr = new DvceArray5D<Real>("rt_Bb", nmb, CK_NB, n1, n3, n2);
          rt_T_ptr  = new DvceArray4D<Real>("rt_T",  nmb, n3, n2, n1);
          rt_pb_ptr = new DvceArray4D<Real>("rt_pb", nmb, n3, n2, n1);
          rt_xT_ptr = new DvceArray4D<Real>("rt_xT", nmb, n3, n2, n1);
          rt_xP_ptr = new DvceArray4D<Real>("rt_xP", nmb, n3, n2, n1);
          rt_icut_ptr = new DvceArray3D<int>("rt_icut", nmb, n3, n2);
          rt_Qb_ptr = new DvceArray5D<Real>("rt_Qb", nmb, nblk, n1, n3, n2);
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
      const bool ck_on = rt_ck;
      // the optical-depth blend with the conduction module's radiative diffusion: its
      // x1-face weight w (rad_w) says how much of each face's longwave flux the
      // two-stream still owns (1 - w); the column's RT bottom is the first face with
      // w = 1, and no internal flux is injected there (the diffusion carries it)
      Conduction *pcond_rt = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->pcond
                                                            : pm->pmb_pack->phydro->pcond;
      const bool taublend = (ck_on && pcond_rt != nullptr && pcond_rt->rad_tau_mode);
      auto w_g = taublend ? pcond_rt->rad_w : DvceArray4D<Real>("rt_w_dummy",1,1,1,1);
      auto tauf_g = taublend ? pcond_rt->rad_tauf
                             : DvceArray4D<Real>("rt_tau_dummy",1,1,1,1);
      if (taublend) int_at_cut = false;
      auto kc_g   = (ck_on) ? *rt_kc_ptr : Fb_g;
      auto Bb_g   = (ck_on) ? *rt_Bb_ptr : Fb_g;
      auto T_g    = (ck_on) ? *rt_T_ptr  : tau_g;
      auto pb_g   = (ck_on) ? *rt_pb_ptr : tau_g;
      auto xT_g   = (ck_on) ? *rt_xT_ptr : tau_g;
      auto xP_g   = (ck_on) ? *rt_xP_ptr : tau_g;
      auto icut_g = (ck_on) ? *rt_icut_ptr : DvceArray3D<int>("dummy",1,1,1);
      auto Qb_g   = (ck_on) ? *rt_Qb_ptr : Fb_g;
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

      if (!rt_ck && n1 > RT_NNC) {
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
      if (ck_on) {
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
          Real pp, TT;
          PresTempFromEint(eos,gm1,Rgas,w0(m,IDN,k,j,i),w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),pp,TT);
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
        par_for("rt_pre_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          if (i < icut_g(m,k,j)) return;          // deeper than the cut: never read
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
                       TT, pbar, w0(m,IDN,k,j,i), kcb);
          const Real sigT4_pi = boltz_sigma/M_PI*SQR(SQR(TT));
          for (int b=0; b<CK_NB; ++b) {
            kc_g(m,b,i,k,j) = kcb[b];
            Bb_g(m,b,i,k,j) = sigT4_pi*ck_planck_frac(ckpf, pfl0, pfid, TT, b);
          }
        });
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
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr;
        get_kapr(T, p, met, kapr);
        // tau = kappa p / g for the UNRESOLVED column above the domain. g must be the
        // value at the top, not the surface value: at r/ap = 1.5 they differ by 2.3x.
        // The tidal term belongs here too: it is the EFFECTIVE gravity that sets how
        // much mass the unresolved column above the domain holds. See EffGravAt.
        Real tau_r_f = kapr*p/EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0, tide);
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = tau_r_f/(kapr*rho);
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
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr;
          get_kapr(T, p, met, kapr);
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
      if (ck_on) {
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
                const Real dtau = kap*ptop*1.0e6/EffGravAt(grav, ap, x1v_(m,ie+1),
                                                           grav_pmass, omega, mu0, tide);
                const RtF trans = RT_EXP(-static_cast<RtF>(dtau/muc[cc]));
                I_down[cc][ie+1] = (static_cast<RtF>(1.0)-trans)
                                 * static_cast<RtF>(Bb_g(m,bandc[cc],ie+1,k,j));
                tausw[cc] = dtau;               // beam already crossed the column above
                transw[cc] = RT_EXP(-static_cast<RtF>(dtau*facsw));
              }
            }
            // down-sweep
            for (int i=ie; i>icut-1; --i) {
              const Real rho = w0(m,IDN,k,j,i);
              const Real drho = rho*dx1(m,k,j,i);
              int iT, iP;
              Real fT, fP;
              const Real xTv = xT_g(m,k,j,i);
              const Real xPv = xP_g(m,k,j,i);
              iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
              iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
              for (int cc=0; cc<NC; ++cc) {
                const int b = bandc[cc];
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                                 + kc_g(m,b,i,k,j);
                const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
                const RtF e0 = -RT_EXPM1(-x);
                const RtF one = static_cast<RtF>(1.0);
                const RtF alp = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                              : (x/2 - x*x/3);
                const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                              : (x/2 - x*x/6);
                I_down[cc][i] = (one-e0)*I_down[cc][i+1]
                              + alp*static_cast<RtF>(Bb_g(m,b,i+1,k,j))
                              + bet*static_cast<RtF>(Bb_g(m,b,i,k,j));
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
              const Real rho = w0(m,IDN,k,j,i-1);
              const Real drho = rho*dx1(m,k,j,i-1);
              int iT, iP;
              Real fT, fP;
              const Real xTv = xT_g(m,k,j,i-1);
              const Real xPv = xP_g(m,k,j,i-1);
              iT = static_cast<int>(xTv); fT = xTv - static_cast<Real>(iT);
              iP = static_cast<int>(xPv); fP = xPv - static_cast<Real>(iP);
              for (int cc=0; cc<NC; ++cc) {
                const int b = bandc[cc];
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc])
                               + kc_g(m,b,i-1,k,j);
                const RtF x = static_cast<RtF>(kap*drho/muc[cc]);
                const RtF e0 = -RT_EXPM1(-x);
                const RtF one = static_cast<RtF>(1.0);
                const RtF bet = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x)
                                                              : (x/2 - x*x/6);
                const RtF gm = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x)
                                                             : (x/2 - x*x/3);
                I_up[cc] = (one-e0)*I_up[cc]
                         + bet*static_cast<RtF>(Bb_g(m,b,i,k,j))
                         + gm*static_cast<RtF>(Bb_g(m,b,i-1,k,j));
                Fb_g(m,blk,i,k,j) += wfc[cc]*(I_up[cc] - I_down[cc][i]);
                // the cell's OWN emission, both hemispheres: in the thin limit each
                // stream adds wfc*(kap*rho*dr/mu)*B over the layer, so per unit volume
                // the two together give 2*(wfc/mu)*kap*rho*B.  Summed over bands and g
                // this is 2*CK_DIFFUSIVITY*sigma*kappa_P*rho*T^4, i.e. the exact
                // 4 sigma kappa_P rho T^4 with 1.66 in place of 2 for the hemispheric
                // mean.  17 % low, which is well inside what a rate estimate needs.
                Em_g(m,blk,i-1,k,j) += 2.0*(wfc[cc]/muc[cc])*kap*rho
                                     * 0.5*(Bb_g(m,b,i,k,j) + Bb_g(m,b,i-1,k,j));
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
            for (int cc=0; cc<NC; ++cc) {
              Real dtauir = gamirc[cc]*dtau_i;
              Real x = dtauir/muggc[cc];
              Real e0 = -expm1(-x);
              Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                                 + alp*fbc[cc]*B[i+1] + bet*fbc[cc]*B[i];
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
              I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]
                            + bet*fbc[cc]*B[i] + gm*fbc[cc]*B[i-1];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
              Fb_g(m,blk,i,k,j) += (F_ir_up_f - F_ir_down_f);
              // the cell's own emission per unit volume, both hemispheres.  dtau_i is
              // kappa rho dr for this layer, so dtau_i/dr = kappa rho and the layer
              // thickness cancels out of the volumetric rate.
              const Real dtau_ly = tau_down_r_f[i-1]-tau_down_r_f[i];
              Em_g(m,blk,i-1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_ly
                                   / dx1(m,k,j,i-1)*fbc[cc]*0.5*(B[i]+B[i-1]);
            }
          }
      });
      }

      // ---- C: reduce over blocks in order, then apply ------------------------------
      int nclip = 0;
      const Real demax = rt_de_max;
      // see rt_apply_debug: which column, and how many calls are left to print
      const bool dbg_on = (rt_apply_debug > 0);
      const int dbg_m = rt_dump_m;
      const int dbg_j = (rt_dump_j >= 0) ? rt_dump_j : (js + je)/2;
      const int dbg_k = (rt_dump_k >= 0) ? rt_dump_k : (ks + ke)/2;
      const int dbg_n = rt_apply_debug_n;
      if (dbg_on) --rt_apply_debug;
      par_reduce_clip4("rt_apply", 0, nmb1, ks, ke, js, je, is, ie, nclip,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nc) {
        Real Ft = 0.0, Fb = 0.0;
        for (int b=0; b<nblk; ++b) {
          Ft += Fb_g(m,b,i+1,k,j);
          Fb += Fb_g(m,b,i,k,j);
        }
        // the two-stream's share of each face in the tau blend
        if (taublend) {
          Ft *= (1.0 - w_g(m,k,j,i+1));
          Fb *= (1.0 - w_g(m,k,j,i));
        }
        Real src = -(Ft-Fb)/dx1(m,k,j,i);
        if (ck_on) {
          // deeper than the cut nothing radiative is applied: that region is optically
          // thick and convective, and the stellar beam died decades of optical depth
          // above
          if (i < icut_g(m,k,j)) {
            src = 0.0;
          } else {
            Real Qs = 0.0;
            for (int b=0; b<nblk; ++b) Qs += Qb_g(m,b,i,k,j);
            src += Qs;
          }
        } else {
          src += Qv_g(m,k,j,i);
        }
        // SEMI-IMPLICIT APPLICATION.  The source splits as src = A - E(T), A being the
        // absorption of the field from elsewhere, fixed on this step, and E the cell's
        // own emission, which is 4 sigma kappa_P rho T^4 and so scales as T^4.  Treating
        // E implicitly and A explicitly, and using de/dT = rho c_v ~ e/T (exact for an
        // ideal gas, and an UNDER-estimate of rho c_v wherever H2 or H is partly
        // dissociated, which only makes the step more damped), gives a linear relaxation
        // with rate lambda = dE/de = (4E/T)/(rho c_v) = 4E/e.  The exact solution over
        // bdt is the exponential below; it reduces to src*bdt when lambda*bdt << 1 and to
        // the equilibrium offset src/lambda when lambda*bdt >> 1, and it can never
        // overshoot the equilibrium.  Pure cooling is then bounded by e/4 per step
        // whatever the timestep, which is what rt_de_max used to impose by hand.
        Real de = src*bdt;
        {
          Real Em = 0.0;
          for (int b=0; b<nblk; ++b) Em += Em_g(m,b,i,k,j);
          if (taublend) {
            Em *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (ck_on && i < icut_g(m,k,j)) Em = 0.0;
          const Real ei = w0(m,IEN,k,j,i);
          if (Em > 0.0 && ei > 0.0) {
            const Real lam = 4.0*Em/ei;
            const Real x = lam*bdt;
            de = (x > 1.0e-4) ? (src/lam)*(-expm1(-x)) : src*bdt;
          }
        }
        if (demax > 0.0) {
          const Real dl = LimitRTSource(de, w0(m,IEN,k,j,i), demax);
          if (dl != de) { ++nc; de = dl; }
        }
        if (dbg_on && m == dbg_m && k == dbg_k && j == dbg_j && i > ie - dbg_n) {
          Real Em = 0.0, Qs = 0.0;
          for (int b=0; b<nblk; ++b) { Em += Em_g(m,b,i,k,j); Qs += Qb_g(m,b,i,k,j); }
          const Real ei = w0(m,IEN,k,j,i);
          Kokkos::printf("rt_apply i=%d T=%.4e d=%.4e e=%.4e Fb=%.6e Ft=%.6e "
                         "divF=%.4e Qs=%.4e Em=%.4e lamdt=%.4e de=%.4e de/e=%.4e "
                         "tau=%.4e w=%.4e dx=%.4e\n",
                         i, T_g(m,k,j,i), w0(m,IDN,k,j,i), ei, Fb, Ft,
                         -(Ft-Fb)/dx1(m,k,j,i), Qs, Em,
                         (Em > 0.0 && ei > 0.0) ? 4.0*Em/ei*bdt : 0.0,
                         de, de/ei, taublend ? tauf_g(m,k,j,i) : 0.0,
                         taublend ? w_g(m,k,j,i) : 0.0, dx1(m,k,j,i));
        }
        u0(m,IEN,k,j,i) += de;
      });
      RTSourceLimiterWarn(nclip);

      // ---- one-shot column dump, for cross-code comparison -------------------------
      if (ck_on && !rt_dump_file.empty() && !rt_dump_done) {
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
          const Real dd = w0(md,IDN,kd,jd,i);
          const Real ee = w0(md,IEN,kd,jd,i);
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
          f << "# deep_hot_jupiter_rt correlated-k column dump\n"
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
          std::cout << "deep_hot_jupiter_rt: wrote correlated-k column dump to '"
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
          Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
          Real rho = w0(m,IDN,k,j,ie+1);
          Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
          B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr;
          get_kapr(T, p, met, kapr);
          Real tau_r_f = kapr*p/EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0, tide);
          tau_down_r_f[ie+1] = tau_r_f;
          Real drtop = tau_r_f/(kapr*rho);
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
            Real rho = w0(m,IDN,k,j,i);
            Real p, T;
            PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                             TGuess(wtemp_, m, k, j, i),p,T);
            B[i] = boltz_sigma/M_PI*SQR(SQR(T));
            Real kapr;
            get_kapr(T, p, met, kapr);
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
  //          Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,i),w0(m,IEN,k,j,i));
  //          Real rho = w0(m,IDN,k,j,i);
  //          Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,i),p);
  //          Real kapr;
  //          get_kapr(T, p, met, kapr);
  //          Real cv = Rgas*rho*igm1;
  //          Real e0 = eos.IsGeneral() ? w0(m,IEN,k,j,i) : cv*T;
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
              const Real dl = LimitRTSource(du, w0(m,IEN,k,j,i), demax_grey);
              if (dl != du) { ++nc; du = dl; }
            }
            u0(m,IEN,k,j,i) += du;
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

    return;
}

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_RT_HPP_
