//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter.cpp
//! \brief Problem generator for the deep hot Jupiter.
//!
//! REFERENCE: Heng, Menou, Phillipps, MNRAS, 413, 2380 (2011); Deitrick, Mendonça, Schroffenegger, Grimm, Tsai, Heng, ApJS, 248, 30 (2020)

// C++ headers
#include <cmath>
#include <iostream> // cout
#include <fstream>  // ifstream, for the correlated-k table
#include <string>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"
#include "pgen_eos_utils.hpp"
#include "diffusion/resistivity.hpp"

#include <Kokkos_Random.hpp>

// EOS-aware conversions shared with the other stratified problem generators
using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt);
void double_gray_two_stream_RT(Mesh *pm, Real bdt);

void picket_fence_two_stream_RT(Mesh *pm, Real bdt);
KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2);
KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr);
KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real &taulim, const Real &A, const Real &B, const Real (&C)[3], const Real (&D)[3], const Real (&E)[3], const Real (&gamv)[3], const Real &tau, Real &T);
template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr);

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);


// Number of (band, quadrature) chains stepped together in the IR sweep. This is the
// instruction-level parallelism the sweep gets; it also fixes the private I_ir_down_c
// footprint at NC radial columns, independent of how many chains are requested.
#ifndef RT_NB
#define RT_NB 4
#endif

// Radial extent of the private I_ir_down_c column in the chain-parallel split kernel.
// The monolithic kernel uses NN = 270 for every private array, which is 18 kB of
// scratch per thread at 88 chains. That does not bind today because the kernel is
// occupancy-starved anyway, but the split kernel is meant to put thousands of waves in
// flight, and then it would. Sized to the real n1 (checked at runtime) instead.
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

namespace {
// --- blocked-band RT scaling harness (measurement only) ---------------------------
// problem/rt_nchain: total (band, quadrature) chains to step, 4 = production picket
// fence. problem/rt_ktab: give the synthetic chains a per-cell table lookup, so the
// memory traffic of a correlated-k scheme is measured and not just the recurrences.
// Chains past the first four carry weight zero, so a harness run must remain bitwise
// identical to production -- see the sweep for why.
int rt_nchain = 4;
bool rt_ktab = false;
constexpr int RT_KT_NT = 30;
constexpr int RT_KT_NP = 20;
// Deliberately leaked: a Kokkos View at namespace scope would be destroyed after
// Kokkos::finalize(). One 88x30x20 table is 422 kB, the size of a real k-table.
DvceArray3D<Real> *rt_ktab_ptr = nullptr;
DvceArray1D<Real> *rt_wgt_ptr = nullptr;

// problem/rt_split: run the RT as three kernels with the chain-block index promoted to
// a parallel dimension, instead of one kernel that loops over chains serially. Same
// arithmetic, same block-summation order -- see the use site.
bool rt_split = false;
DvceArray4D<Real> *rt_tau_ptr = nullptr;   // (m,k,j,i) face optical depth from the top
DvceArray4D<Real> *rt_B_ptr = nullptr;     // (m,k,j,i) Planck function
DvceArray4D<Real> *rt_Qv_ptr = nullptr;    // (m,k,j,i) stellar heating rate
DvceArray4D<Real> *rt_cf_ptr = nullptr;    // (m,k,j,{gamir1,gamir2,beta})
DvceArray5D<Real> *rt_Fb_ptr = nullptr;    // (m,blk,k,j,i) IR flux, one slot per block

// Correlated-k per-cell fields, filled by rt_pre and read by the chain kernel. Laid out
// (m,band,k,j,i) so that the radial index is contiguous for a thread sweeping one band --
// (m,k,j,i,band) would stride every read of the sweep by CK_NB.
DvceArray5D<Real> *rt_kc_ptr = nullptr;    // (m,b,k,j,i) continuum kappa [cm^2/g]
DvceArray5D<Real> *rt_Bb_ptr = nullptr;    // (m,b,k,j,i) band Planck intensity
DvceArray4D<Real> *rt_T_ptr = nullptr;     // (m,k,j,i) temperature [K]
DvceArray4D<Real> *rt_pb_ptr = nullptr;    // (m,k,j,i) pressure [bar]
DvceArray4D<Real> *rt_xT_ptr = nullptr;    // (m,k,j,i) continuous k-table index in T
DvceArray4D<Real> *rt_xP_ptr = nullptr;    // (m,k,j,i) continuous k-table index in p
DvceArray3D<int>  *rt_icut_ptr = nullptr;  // (m,k,j) deepest cell doing correlated-k
DvceArray5D<Real> *rt_Qb_ptr = nullptr;    // (m,blk,k,j,i) stellar heating, one per block

// --- correlated-k (Lee/Exo-FMS premixed tables, Kataria+2013 11-band grid) --------
// problem/rt_ck turns it on; problem/ck_table is the path to the premixed table and
// problem/ck_pcut_bar the pressure below which correlated-k is used at all (deeper than
// that the atmosphere is optically thick -- grey tau ~ 9e3 at 10 bar for this setup --
// and the interior is far outside any molecular table, up to 12000 K).
//
// Table layout and the traps it carries are documented in data/exo_fms_ck/PROVENANCE.md.
// The two that matter here: the data block runs the band index BACKWARDS, and kappa is
// cgs cm^2/g, which is what this code already works in.
bool rt_ck = false;
Real rt_ck_pcut = 10.0;                    // bar
// problem/ck_nquad: angular treatment of the longwave. 1 uses a single diffusivity
// factor, mu = 1/1.66, which is what GCMs normally do and what the chain count in the
// correlated-k literature assumes. 2 keeps the 2-point Gauss quadrature the picket fence
// uses, which doubles the number of column solves. Both normalise to F = pi <I>: for the
// Gauss pair, 2 pi sum w mu = 2 pi (0.5) = pi.
int ck_nq = 1;
constexpr Real CK_DIFFUSIVITY = 1.66;
constexpr int CK_NB = 11;                  // bands, fixed by the Kataria grid
constexpr int CK_NG = 8;                   // g-points per band, fixed by the table
int ck_nT = 0;
int ck_nP = 0;
// Deliberately leaked, as with the other tables here: a namespace-scope Kokkos View
// would be destroyed after Kokkos::finalize().
DvceArray1D<Real> *ck_lT_ptr = nullptr;    // log10 T grid [K]
DvceArray1D<Real> *ck_lP_ptr = nullptr;    // log10 p grid [bar]
// Layout is (band, g, iT, iP): for one chain the entire (T,p) plane is contiguous and only
// 38*34*8 = 10 kB, so the four bilinear corners sit 8 and 272 bytes apart rather than the
// 704 bytes and 24 kB that (iT,iP,band,g) would give.
//
// This was expected to be worth about a factor of two in the chain kernel. It was MEASURED
// AND IT IS NOT: 2725 -> 2711 ms per 100 cycles, half a per cent. The k-table access
// pattern is simply not the bottleneck, which is the same answer an earlier experiment
// with a synthetic table gave when it blocked the lookup by band and gained nothing.
// The layout is kept because it is the more natural one, not because it is faster.
DvceArray4D<Real> *ck_lk_ptr = nullptr;    // (b,g,iT,iP) log10 kappa [cm^2/g]
DvceArray1D<Real> *ck_gw_ptr = nullptr;    // g-point weights, sum to 1
DvceArray1D<Real> *ck_wl_ptr = nullptr;    // band edges [um], descending, CK_NB+1 of them
DvceArray1D<Real> *ck_swf_ptr = nullptr;   // stellar flux FRACTION per band, sums to 1

// Band-integrated Planck function. The grey scheme uses B = sigma T^4 / pi; correlated-k
// needs B_b(T) = (sigma T^4 / pi) * f_b(T), the fraction of the Planck function falling in
// band b. f_b is smooth in log T, so it is tabulated once on a uniform log10 T grid and
// interpolated -- evaluating the series per cell per band would be absurd.
constexpr int CK_NPF = 512;
constexpr Real CK_PF_TMIN = 50.0;
constexpr Real CK_PF_TMAX = 20000.0;
DvceArray2D<Real> *ck_pf_ptr = nullptr;    // (iT, band) fractional Planck function
Real ck_pf_lTmin = 0.0;
Real ck_pf_idlT = 0.0;                     // 1 / grid spacing in log10 T

//----------------------------------------------------------------------------------------
//! \fn Real planck_fraction_below()
//  \brief fraction of blackbody emission at wavelengths SHORTER than lam, as a function of
//  lam*T in um K. The standard series (Chang & Rhee 1984): with xi = c2/(lam T),
//    F = 15/pi^4 sum_n e^{-n xi}/n (xi^3 + 3 xi^2/n + 6 xi/n^2 + 6/n^3).
//  Host side only, evaluated once per table entry at startup.

Real planck_fraction_below(const Real lamT) {
  const Real c2 = 1.4387769e4;             // hc/k in um K
  if (lamT <= 0.0) return 0.0;
  const Real xi = c2/lamT;
  if (xi > 7.0e2) return 0.0;              // exp underflow: nothing below this lam
  Real sum = 0.0;
  for (int n=1; n<=500; ++n) {
    const Real nx = n*xi;
    if (nx > 7.0e2) break;
    const Real e = std::exp(-nx);
    const Real rn = 1.0/static_cast<Real>(n);
    sum += e*rn*(xi*xi*xi + 3.0*xi*xi*rn + 6.0*xi*rn*rn + 6.0*rn*rn*rn);
  }
  const Real pi4 = M_PI*M_PI*M_PI*M_PI;
  return 15.0/pi4*sum;
}

//----------------------------------------------------------------------------------------
//! \fn void build_planck_fractions()
//  \brief tabulate f_b(T) on a uniform log10 T grid, given the band edges already read.
//
//  The Kataria grid spans 0.26 to 324.68 um and does NOT capture the whole Planck function:
//  at high T a real fraction escapes past the blue edge. That flux has to go somewhere or
//  the scheme silently loses energy, so the two outermost bands are extended to 0 and
//  infinity -- the sub-0.26 um tail joins the bluest band and the super-324.68 um tail the
//  reddest. It is an approximation, since those tails get their host band's kappa, but it
//  is a small one wherever the tails are small, and it makes sum_b f_b = 1 exactly.

void build_planck_fractions() {
  ck_pf_ptr = new DvceArray2D<Real>("ck_pf", CK_NPF, CK_NB);
  auto hpf = Kokkos::create_mirror_view(*ck_pf_ptr);
  auto hwl = Kokkos::create_mirror_view(*ck_wl_ptr);
  Kokkos::deep_copy(hwl, *ck_wl_ptr);

  const Real lTmin = std::log10(CK_PF_TMIN);
  const Real lTmax = std::log10(CK_PF_TMAX);
  const Real dlT = (lTmax - lTmin)/static_cast<Real>(CK_NPF-1);
  ck_pf_lTmin = lTmin;
  ck_pf_idlT = 1.0/dlT;

  Real blue1pc_T = -1.0;                   // T above which >1% of the flux is bluer than
                                           // the grid: the band structure stops capturing
                                           // the spectrum and the fold-in stops being small
  Real worst_sum_err = 0.0;
  for (int i=0; i<CK_NPF; ++i) {
    const Real T = std::pow(10.0, lTmin + i*dlT);
    // hwl is descending, so band b runs from hwl(b) (long) to hwl(b+1) (short) and
    // F(long) > F(short).
    for (int b=0; b<CK_NB; ++b) {
      hpf(i,b) = planck_fraction_below(hwl(b)*T) - planck_fraction_below(hwl(b+1)*T);
    }
    const Real red_tail = 1.0 - planck_fraction_below(hwl(0)*T);
    const Real blue_tail = planck_fraction_below(hwl(CK_NB)*T);
    hpf(i,0) += red_tail;
    hpf(i,CK_NB-1) += blue_tail;
    if (blue1pc_T < 0.0 && blue_tail > 0.01) blue1pc_T = T;
    Real sum = 0.0;
    for (int b=0; b<CK_NB; ++b) sum += hpf(i,b);
    worst_sum_err = std::max(worst_sum_err, std::fabs(sum - 1.0));
  }
  Kokkos::deep_copy(*ck_pf_ptr, hpf);

  if (global_variable::my_rank == 0) {
    std::cout << "  Planck fractions tabulated on " << CK_NPF << " points, T = "
              << CK_PF_TMIN << " .. " << CK_PF_TMAX << " K; worst |sum_b f_b - 1| = "
              << worst_sum_err << std::endl
              << "  >1% of the Planck function falls bluer than " << hwl(CK_NB)
              << " um above T = " << blue1pc_T << " K; beyond that the 11-band grid stops"
              << std::endl
              << "  capturing the spectrum and folding the tail into the bluest band is no"
              << " longer a small correction" << std::endl;
    if (blue1pc_T > 0.0 && blue1pc_T < CK_PF_TMAX) {
      std::cout << "  (with the p < " << rt_ck_pcut << " bar cut this setup tops out near "
                << "4800 K, where the tail is 3e-3)" << std::endl;
    }
    if (worst_sum_err > 1.0e-12) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: Planck fractions do not sum "
                << "to 1 (worst error " << worst_sum_err << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ck_planck_bands()
//  \brief B_b(T) for all bands at one cell, in the same units as the grey B, i.e. the
//  Planck INTENSITY sigma T^4 / pi split across the bands. Linear in log10 T on a uniform
//  grid, so the index is analytic and there is no search.

//----------------------------------------------------------------------------------------
//! \fn Real ck_planck_frac()
//  \brief f_b(T) for a single band. Uniform log10 T grid, so the index is arithmetic.

KOKKOS_INLINE_FUNCTION
Real ck_planck_frac(const DvceArray2D<Real> &pf, const Real lTmin, const Real idlT,
                    const Real &T, const int b) {
  Real x = (log10(T) - lTmin)*idlT;
  x = (x < 0.0) ? 0.0 : ((x > CK_NPF-1.0) ? CK_NPF-1.0 : x);
  int i = static_cast<int>(x);
  i = (i > CK_NPF-2) ? CK_NPF-2 : i;
  const Real f = x - static_cast<Real>(i);
  return (1.0-f)*pf(i,b) + f*pf(i+1,b);
}

KOKKOS_INLINE_FUNCTION
void ck_planck_bands(const DvceArray2D<Real> &pf, const Real lTmin, const Real idlT,
                     const Real &sigT4_pi, const Real &T, Real (&Bb)[CK_NB]) {
  Real x = (log10(T) - lTmin)*idlT;
  x = (x < 0.0) ? 0.0 : ((x > CK_NPF-1.0) ? CK_NPF-1.0 : x);
  int i = static_cast<int>(x);
  i = (i > CK_NPF-2) ? CK_NPF-2 : i;
  const Real f = x - static_cast<Real>(i);
  for (int b=0; b<CK_NB; ++b) {
    Bb[b] = sigT4_pi*((1.0-f)*pf(i,b) + f*pf(i+1,b));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void read_ck_table()
//  \brief host-side reader for the Exo-FMS premixed correlated-k table (the HELIOS-k
//  format, ck_form == 2 in Exo-FMS's src/ck_opacity_mod.f90). Everything past the first
//  line is whitespace-separated numbers, so it is read as one token stream.
//
//  kappa is stored as log10 with a 1e-99 floor, and the T and p grids as log10, because
//  that is the space the interpolation has to happen in: kappa spans ~40 decades across
//  the table and linear interpolation of it is meaningless.

void read_ck_table(const std::string &fname) {
  std::ifstream f(fname);
  if (!f.is_open()) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open correlated-k "
              << "table '" << fname << "'. Set problem/ck_table." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string species;
  std::getline(f, species);                       // line 1: species list, text
  int nT, nP, nb, ng;
  f >> nT >> nP >> nb >> ng;
  if (nb != CK_NB || ng != CK_NG) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: correlated-k table is "
              << nb << " bands x " << ng << " g-points, but this build is compiled for "
              << CK_NB << " x " << CK_NG << ". Use the 11-band g8 table." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  ck_nT = nT;
  ck_nP = nP;
  ck_lT_ptr = new DvceArray1D<Real>("ck_lT", nT);
  ck_lP_ptr = new DvceArray1D<Real>("ck_lP", nP);
  ck_lk_ptr = new DvceArray4D<Real>("ck_lk", CK_NB, CK_NG, nT, nP);
  ck_gw_ptr = new DvceArray1D<Real>("ck_gw", CK_NG);
  ck_wl_ptr = new DvceArray1D<Real>("ck_wl", CK_NB+1);
  auto hlT = Kokkos::create_mirror_view(*ck_lT_ptr);
  auto hlP = Kokkos::create_mirror_view(*ck_lP_ptr);
  auto hlk = Kokkos::create_mirror_view(*ck_lk_ptr);
  auto hgw = Kokkos::create_mirror_view(*ck_gw_ptr);
  auto hwl = Kokkos::create_mirror_view(*ck_wl_ptr);

  Real v;
  for (int i=0; i<nT; ++i) { f >> v; hlT(i) = std::log10(v); }
  for (int j=0; j<nP; ++j) { f >> v; hlP(j) = std::log10(v); }
  for (int b=0; b<CK_NB+1; ++b) { f >> v; hwl(b) = v; }        // um, descending
  for (int b=0; b<CK_NB+1; ++b) { f >> v; }                    // wavenumbers, unused
  for (int g=0; g<CK_NG; ++g) { f >> v; }                      // g nodes, unused: each
                                                               // g-point is its own
                                                               // column solve
  for (int g=0; g<CK_NG; ++g) { f >> v; hgw(g) = v; }
  // ORDERING. Records run in the same order as the wl edges, i.e. DESCENDING wavelength:
  // the first record of each (T,p) block is 324.68-20 um and the last is 0.26-0.42 um.
  //
  // Note this is the opposite of what Exo-FMS's own reader appears to do -- its loop is
  // `do b = nwl, 1, -1` -- so it was checked against the data instead of trusted. The
  // discriminator is condensation: with band 10 = 0.26-0.42 um the optical opacity at
  // 0.1 bar goes 1.3e-8, 1.3e-7, 0.71, 31, 60 cm^2/g at T = 300, 800, 1500, 2500,
  // 3500 K, which is TiO/VO/Fe/Na/K appearing as they vaporise, while band 0 =
  // 20-324.68 um falls from 14 cm^2/g with T, which is the H2O rotational band. Reversed,
  // both are physically impossible. Getting this wrong silently swaps the optical and the
  // far infrared, which looks plausible in a plot and is completely wrong.
  for (int i=0; i<nT; ++i) {
    for (int j=0; j<nP; ++j) {
      for (int b=0; b<CK_NB; ++b) {
        for (int g=0; g<CK_NG; ++g) {
          f >> v;
          hlk(b,g,i,j) = std::log10((v > 1.0e-99) ? v : 1.0e-99);
        }
      }
    }
  }
  if (!f) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: correlated-k table '" << fname
              << "' ended early; expected " << nT*nP*CK_NB*CK_NG << " kappa values."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Kokkos::deep_copy(*ck_lT_ptr, hlT);
  Kokkos::deep_copy(*ck_lP_ptr, hlP);
  Kokkos::deep_copy(*ck_lk_ptr, hlk);
  Kokkos::deep_copy(*ck_gw_ptr, hgw);
  Kokkos::deep_copy(*ck_wl_ptr, hwl);

  if (global_variable::my_rank == 0) {
    Real wsum = 0.0;
    for (int g=0; g<CK_NG; ++g) wsum += hgw(g);
    std::cout << "deep_hot_jupiter_rt: correlated-k table '" << fname << "'" << std::endl
              << "  " << nT << " T x " << nP << " p x " << CK_NB << " bands x " << CK_NG
              << " g, T = " << std::pow(10.0, hlT(0)) << " .. "
              << std::pow(10.0, hlT(nT-1)) << " K, p = " << std::pow(10.0, hlP(0))
              << " .. " << std::pow(10.0, hlP(nP-1)) << " bar" << std::endl
              << "  bands " << hwl(CK_NB) << " .. " << hwl(0)
              << " um, g weights sum to " << wsum << std::endl
              << "  correlated-k applied where p < " << rt_ck_pcut << " bar" << std::endl;
    // the grids must be increasing: the lookup does a bracketing search on them
    for (int i=1; i<nT; ++i) {
      if (hlT(i) <= hlT(i-1)) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: ck table T grid is not "
                  << "increasing at index " << i << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    for (int j=1; j<nP; ++j) {
      if (hlP(j) <= hlP(j-1)) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: ck table p grid is not "
                  << "increasing at index " << j << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (std::fabs(wsum - 1.0) > 1.0e-10) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: ck g weights sum to " << wsum
                << ", not 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return;
}

// --- continuum: CIA, Rayleigh, and the equilibrium composition they need ----------
// The k-table is line opacity only. CIA and Rayleigh are grey WITHIN a band, so they add
// to every g-point rather than making new chains:
//     kappa_tot(b,g) = kappa_ck(b,g) + kappa_CIA(b) + kappa_Ray(b).
// Both need number densities, so they need the equilibrium composition. That comes from
// the FastChem table shipped alongside the k-table -- same chemistry the premixed
// opacities were built with, which is the point: mixing a different chemistry into the
// continuum than into the lines would be inconsistent.
constexpr int CK_NCIA = 4;                 // H2-H2, H2-He, H2-H, He-H
constexpr int CK_CIA_NTMAX = 512;
constexpr int CK_NRAY = 4;                 // H2, He, H, e-
constexpr int CK_NCE = 6;                  // mu, then VMR of H2, He, H, e-, H-
int ce_nT = 0;
int ce_nP = 0;
DvceArray1D<Real> *ce_lT_ptr = nullptr;
DvceArray1D<Real> *ce_lP_ptr = nullptr;
DvceArray3D<Real> *ce_ptr = nullptr;       // (iT,iP,CK_NCE)
DvceArray1D<int>  *cia_nT_ptr = nullptr;   // per-pair grid length
DvceArray2D<Real> *cia_T_ptr = nullptr;    // (pair,iT) -- the four grids differ wildly,
DvceArray3D<Real> *cia_k_ptr = nullptr;    // (pair,iT,band)   200-3000 K to 200-9900 K
DvceArray2D<Real> *ray_x_ptr = nullptr;    // (species,band) cm^2/molecule, T independent

//----------------------------------------------------------------------------------------
//! \fn void read_ck_continuum()
//  \brief read the FastChem composition table, the four CIA pair tables and the Rayleigh
//  cross sections. All are whitespace-separated numbers after a short header.

void read_ck_continuum(const std::string &dir, const std::string &swfile) {
  // ---- FastChem composition: "nT nP nrec nspecies", species names, T grid, p grid,
  // then nrec records of {mu, VMR(H2), VMR(He), VMR(H), VMR(e-), VMR(H-)}. Note SIX
  // columns for five species: mu is prepended.
  {
    const std::string fn = dir + "/CE_tables/FastChem_ck_1x_int.txt";
    std::ifstream f(fn);
    if (!f.is_open()) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open '" << fn
                << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int nrec, nsp;
    f >> ce_nT >> ce_nP >> nrec >> nsp;
    std::string nm;
    for (int n=0; n<nsp; ++n) f >> nm;
    ce_lT_ptr = new DvceArray1D<Real>("ce_lT", ce_nT);
    ce_lP_ptr = new DvceArray1D<Real>("ce_lP", ce_nP);
    ce_ptr = new DvceArray3D<Real>("ce", ce_nT, ce_nP, CK_NCE);
    auto hT = Kokkos::create_mirror_view(*ce_lT_ptr);
    auto hP = Kokkos::create_mirror_view(*ce_lP_ptr);
    auto hC = Kokkos::create_mirror_view(*ce_ptr);
    Real v;
    for (int i=0; i<ce_nT; ++i) { f >> v; hT(i) = std::log10(v); }
    for (int j=0; j<ce_nP; ++j) { f >> v; hP(j) = std::log10(v); }
    for (int i=0; i<ce_nT; ++i) {
      for (int j=0; j<ce_nP; ++j) {
        for (int n=0; n<CK_NCE; ++n) { f >> v; hC(i,j,n) = v; }
      }
    }
    if (!f) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: '" << fn
                << "' ended early" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Kokkos::deep_copy(*ce_lT_ptr, hT);
    Kokkos::deep_copy(*ce_lP_ptr, hP);
    Kokkos::deep_copy(*ce_ptr, hC);
  }
  // ---- CIA pairs: "nT nband", T grid, band wavenumbers, then nT rows of nband values.
  {
    const char *pf[CK_NCIA] = {"H2-H2", "H2-He", "H2-H", "He-H"};
    cia_nT_ptr = new DvceArray1D<int>("cia_nT", CK_NCIA);
    cia_T_ptr = new DvceArray2D<Real>("cia_T", CK_NCIA, CK_CIA_NTMAX);
    cia_k_ptr = new DvceArray3D<Real>("cia_k", CK_NCIA, CK_CIA_NTMAX, CK_NB);
    auto hn = Kokkos::create_mirror_view(*cia_nT_ptr);
    auto hT = Kokkos::create_mirror_view(*cia_T_ptr);
    auto hk = Kokkos::create_mirror_view(*cia_k_ptr);
    Kokkos::deep_copy(hT, 0.0);
    Kokkos::deep_copy(hk, 0.0);
    for (int s=0; s<CK_NCIA; ++s) {
      const std::string fn = dir + "/cia/" + pf[s] + "_reform_11.txt";
      std::ifstream f(fn);
      if (!f.is_open()) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open '" << fn
                  << "'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      int nT, nb;
      f >> nT >> nb;
      if (nb != CK_NB || nT > CK_CIA_NTMAX) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: CIA table '" << fn
                  << "' is " << nT << " x " << nb << ", need <= " << CK_CIA_NTMAX
                  << " x " << CK_NB << std::endl;
        std::exit(EXIT_FAILURE);
      }
      hn(s) = nT;
      Real v;
      for (int i=0; i<nT; ++i) { f >> v; hT(s,i) = v; }
      for (int b=0; b<CK_NB; ++b) { f >> v; }          // band wavenumbers, unused
      for (int i=0; i<nT; ++i) {
        for (int b=0; b<CK_NB; ++b) { f >> v; hk(s,i,b) = v; }
      }
      if (!f) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: '" << fn << "' ended early"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    Kokkos::deep_copy(*cia_nT_ptr, hn);
    Kokkos::deep_copy(*cia_T_ptr, hT);
    Kokkos::deep_copy(*cia_k_ptr, hk);
  }
  // ---- stellar flux per band. Ordered ASCENDING in wavelength, i.e. OPPOSITE to the
  // k-table and to wl -- verified by reversing it and recovering a 6460 K blackbody to
  // about 1 % per band, where as listed it is exactly backwards. Only the SHAPE is taken:
  // the values are renormalised to sum to one and multiplied by the code's own
  // sigma T_irr^4, so the total insolation matches the grey scheme it replaces and the
  // file's absolute normalisation never has to be pinned down.
  {
    const std::string fn = dir + "/sw_flux/" + swfile;
    std::ifstream f(fn);
    if (!f.is_open()) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open '" << fn
                << "'. Set problem/ck_swflux." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ck_swf_ptr = new DvceArray1D<Real>("ck_swf", CK_NB);
    auto hs = Kokkos::create_mirror_view(*ck_swf_ptr);
    Real v[CK_NB];
    Real tot = 0.0;
    for (int b=0; b<CK_NB; ++b) { f >> v[b]; tot += v[b]; }
    if (!f || tot <= 0.0) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: bad stellar flux file '"
                << fn << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int b=0; b<CK_NB; ++b) hs(b) = v[CK_NB-1-b]/tot;    // reverse into wl order
    Kokkos::deep_copy(*ck_swf_ptr, hs);
    if (global_variable::my_rank == 0) {
      std::cout << "  stellar spectrum '" << swfile << "', band fractions "
                << hs(CK_NB-1) << " (bluest) .. " << hs(0) << " (reddest)" << std::endl;
    }
  }
  // ---- Rayleigh: one species name, then CK_NB cross sections in cm^2/molecule.
  {
    const char *rf[CK_NRAY] = {"H2", "He", "H", "e-"};
    ray_x_ptr = new DvceArray2D<Real>("ray_x", CK_NRAY, CK_NB);
    auto hr = Kokkos::create_mirror_view(*ray_x_ptr);
    for (int s=0; s<CK_NRAY; ++s) {
      const std::string fn = dir + "/ray/Ray_" + rf[s] + "_11.txt";
      std::ifstream f(fn);
      if (!f.is_open()) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open '" << fn
                  << "'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string nm;
      std::getline(f, nm);
      for (int b=0; b<CK_NB; ++b) { f >> hr(s,b); }
      if (!f) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: '" << fn << "' ended early"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    Kokkos::deep_copy(*ray_x_ptr, hr);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "  continuum: FastChem composition " << ce_nT << " T x " << ce_nP
              << " p, " << CK_NCIA << " CIA pairs, " << CK_NRAY << " Rayleigh species"
              << std::endl;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ck_tp_index()
//  \brief bracket x in a monotonically increasing grid, returning the lower index and the
//  interpolation fraction. CLAMPED, not extrapolated, at both ends.
//
//  The clamping is the whole point. Clamping only the index and letting the fraction run
//  free is what made the earlier synthetic-table harness produce a negative opacity, a
//  sign-flipped optical depth increment, and an intensity recurrence that diverged to inf.
//  A table lookup that silently extrapolates 40 decades of kappa is not a lookup.

KOKKOS_INLINE_FUNCTION
void ck_tp_index(const DvceArray1D<Real> &lg, const int n, const Real &x,
                 int &i, Real &f) {
  if (x <= lg(0)) { i = 0; f = 0.0; return; }
  if (x >= lg(n-1)) { i = n-2; f = 1.0; return; }
  int lo = 0;
  int hi = n-1;
  while (hi - lo > 1) {
    const int mid = (lo + hi)/2;
    if (x < lg(mid)) { hi = mid; } else { lo = mid; }
  }
  i = lo;
  f = (x - lg(lo))/(lg(lo+1) - lg(lo));
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ck_kappa()
//  \brief kappa [cm^2/g] for one (band, g-point), bilinear in log10 kappa over
//  (log10 T, log10 p). The table is stored as log10 because kappa spans about forty
//  decades; interpolating it linearly would be meaningless.
//
//  (iT, fT, iP, fP) come from ck_tp_index and depend only on the cell, so a caller
//  handling several chains in one cell computes them once. Measurement says the index
//  arithmetic is not the cost either way.

KOKKOS_INLINE_FUNCTION
Real ck_kappa(const DvceArray4D<Real> &lk, const int iT, const Real &fT,
              const int iP, const Real &fP, const int b, const int g) {
  const Real k00 = lk(b, g, iT  , iP  );
  const Real k01 = lk(b, g, iT  , iP+1);
  const Real k10 = lk(b, g, iT+1, iP  );
  const Real k11 = lk(b, g, iT+1, iP+1);
  const Real lkap = (1.0-fT)*((1.0-fP)*k00 + fP*k01)
                  +      fT *((1.0-fP)*k10 + fP*k11);
  return exp(2.302585092994046*lkap);       // 10^lkap
}

//----------------------------------------------------------------------------------------
//! \fn void ck_continuum()
//  \brief grey-within-band continuum opacity [cm^2/g] for every band at one cell: the four
//  CIA pairs plus Rayleigh scattering off H2, He, H and e-.
//
//  CIA scales as the PRODUCT of the two collider number densities, Rayleigh as one, so
//  both need the equilibrium composition, taken from the FastChem table. The simulation's
//  own rho is used for the mass conversion rather than the table's, since that is the
//  density the rest of the scheme works with.
//
//  Every table index is clamped, and the four CIA pairs have grids that stop in very
//  different places -- 200-3000 K for H2-H2 but 200-9900 K for H2-He -- so clamping is not
//  an edge case here, it is the normal state of affairs above 3000 K. Held-flat
//  extrapolation of CIA is what Exo-FMS does too.
//
//  NOT included: H- bound-free and free-free. Those dominate the continuum above about
//  3000 K and are the reason the table's Rosseland mean still falls short of the grey
//  Freedman opacity there -- see the note in read_ck_table's validation.

KOKKOS_INLINE_FUNCTION
void ck_continuum(const DvceArray3D<Real> &ce, const DvceArray1D<Real> &celT,
                  const DvceArray1D<Real> &celP, const int ceNT, const int ceNP,
                  const DvceArray1D<int> &ciaN, const DvceArray2D<Real> &ciaT,
                  const DvceArray3D<Real> &ciak, const DvceArray2D<Real> &rayx,
                  const DvceArray1D<Real> &wl,
                  const Real &T, const Real &pbar, const Real &rho, Real (&kc)[CK_NB]) {
  // composition at this (T,p)
  int iT, iP;
  Real fT, fP;
  ck_tp_index(celT, ceNT, log10(T), iT, fT);
  ck_tp_index(celP, ceNP, log10(pbar), iP, fP);
  Real vmr[CK_NCE];
  for (int n=0; n<CK_NCE; ++n) {
    vmr[n] = (1.0-fT)*((1.0-fP)*ce(iT  ,iP,n) + fP*ce(iT  ,iP+1,n))
           +      fT *((1.0-fP)*ce(iT+1,iP,n) + fP*ce(iT+1,iP+1,n));
  }
  // total number density, ideal gas. vmr[0] is mu and is not a species.
  const Real kboltz = 1.380649e-16;
  const Real ntot = pbar*1.0e6/(kboltz*T);
  const Real irho = 1.0/rho;

  for (int b=0; b<CK_NB; ++b) kc[b] = 0.0;

  // CIA: pair (i1,i2) indexes into vmr, where 1=H2 2=He 3=H 4=e- 5=H-
  const int i1[CK_NCIA] = {1, 1, 1, 2};
  const int i2[CK_NCIA] = {1, 2, 3, 3};
  for (int s=0; s<CK_NCIA; ++s) {
    const int n = ciaN(s);
    int it = 0;
    Real ft = 0.0;
    if (T <= ciaT(s,0)) {
      it = 0; ft = 0.0;
    } else if (T >= ciaT(s,n-1)) {
      it = n-2; ft = 1.0;
    } else {
      int lo = 0;
      int hi = n-1;
      while (hi - lo > 1) {
        const int mid = (lo + hi)/2;
        if (T < ciaT(s,mid)) { hi = mid; } else { lo = mid; }
      }
      it = lo;
      ft = (T - ciaT(s,lo))/(ciaT(s,lo+1) - ciaT(s,lo));
    }
    const Real nn = vmr[i1[s]]*ntot*vmr[i2[s]]*ntot*irho;
    for (int b=0; b<CK_NB; ++b) {
      kc[b] += ((1.0-ft)*ciak(s,it,b) + ft*ciak(s,it+1,b))*nn;
    }
  }
  // Rayleigh: cross section per molecule, no temperature dependence
  for (int s=0; s<CK_NRAY; ++s) {
    const Real nn = vmr[s+1]*ntot*irho;
    for (int b=0; b<CK_NB; ++b) {
      kc[b] += rayx(s,b)*nn;
    }
  }

  // H- bound-free and free-free, John (1988). This is the dominant continuum above about
  // 3000 K, which under the p < 10 bar cut is the BOTTOM of the correlated-k region, so it
  // is not a refinement. n(H-) is taken straight from the FastChem table rather than
  // reconstructed from Saha, which keeps it consistent with the line opacities.
  {
    const Real lam0 = 1.6419;              // um, H- photodetachment threshold
    const Real Cbf[6] = {152.519, 49.534, -118.858, 92.536, -34.194, 4.982};
    // free-free coefficients: set 1 for 0.1823 < lam < 0.3645 um, set 2 for lam >= 0.3645
    const Real Aff1[6] = {518.1021, 472.2636, -482.2089, 115.5291, 0.0, 0.0};
    const Real Bff1[6] = {-734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0};
    const Real Cff1[6] = {1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0};
    const Real Dff1[6] = {-479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0};
    const Real Eff1[6] = {93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0};
    const Real Fff1[6] = {-6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0};
    const Real Aff2[6] = {0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830};
    const Real Bff2[6] = {0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170};
    const Real Cff2[6] = {0.0, -2054.2910, 8746.5230, -13651.1050, 8642.9700, -1863.8640};
    const Real Dff2[6] = {0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880};
    const Real Eff2[6] = {0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880};
    const Real Fff2[6] = {0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850};
    const Real T5040 = 5040.0/T;
    const Real nHm = vmr[5]*ntot;                    // H- number density
    const Real Pe_nH = vmr[4]*ntot*vmr[3]*ntot*kboltz*T;   // P(e-) * n(H)
    for (int b=0; b<CK_NB; ++b) {
      // the band's representative wavelength is the one at the mean WAVENUMBER of its
      // edges, which is what the binned CIA tables were built on
      const Real lam = 2.0/(1.0/wl(b) + 1.0/wl(b+1));
      // bound-free: zero longward of the detachment threshold
      Real xbf = 0.0;
      if (lam < lam0) {
        const Real dk = 1.0/lam - 1.0/lam0;
        // exponents are n/2, so walk them with a running sqrt instead of six pow()s
        const Real sdk = sqrt(dk);
        Real dp = 1.0;
        Real fbf = 0.0;
        for (int n=0; n<6; ++n) { fbf += Cbf[n]*dp; dp *= sdk; }
        xbf = 1.0e-18*lam*lam*lam*(dk*sdk)*fbf;
      }
      // free-free
      Real sff = 0.0;
      if (lam >= 0.3645 || (lam > 0.1823 && lam < 0.3645)) {
        const bool set2 = (lam >= 0.3645);
        const Real st = sqrt(T5040);
        Real tp = T5040;                     // exponent (n+2)/2, walked by sqrt(T5040)
        for (int n=0; n<6; ++n) {
          const Real An = set2 ? Aff2[n] : Aff1[n];
          const Real Bn = set2 ? Bff2[n] : Bff1[n];
          const Real Cn = set2 ? Cff2[n] : Cff1[n];
          const Real Dn = set2 ? Dff2[n] : Dff1[n];
          const Real En = set2 ? Eff2[n] : Eff1[n];
          const Real Fn = set2 ? Fff2[n] : Fff1[n];
          sff += tp*(lam*lam*An + Bn + Cn/lam
                 + Dn/(lam*lam) + En/(lam*lam*lam) + Fn/(lam*lam*lam*lam));
          tp *= st;
        }
      }
      // xbf is cm^2 per H-, sff*1e-29 is cm^4/dyne and multiplies P(e-) n(H)
      kc[b] += (xbf*nHm + 1.0e-29*sff*Pe_nH)*irho;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ck_selftest()
//  \brief exercise the k-table on the DEVICE at startup and check the one thing that can
//  go wrong silently: the band ordering.
//
//  The discriminator is condensation. TiO, VO, Fe, Na and K are condensed out below about
//  1000 K and in the gas phase by 2500 K, so the bluest band's opacity has to climb by
//  orders of magnitude between the two, while the reddest band -- the H2O rotational band
//  -- falls. If the band index is reversed these swap, which looks entirely plausible in a
//  plot and is completely wrong. Exo-FMS's own reader loops the other way, so this is not
//  a hypothetical.

void ck_selftest() {
  auto lk = *ck_lk_ptr;
  auto lT = *ck_lT_ptr;
  auto lP = *ck_lP_ptr;
  auto gw = *ck_gw_ptr;
  const int nT = ck_nT;
  const int nP = ck_nP;
  auto ce = *ce_ptr;
  auto celT = *ce_lT_ptr;
  auto celP = *ce_lP_ptr;
  auto ciaN = *cia_nT_ptr;
  auto ciaT = *cia_T_ptr;
  auto ciak = *cia_k_ptr;
  auto rayx = *ray_x_ptr;
  auto wl = *ck_wl_ptr;
  const int ceNT = ce_nT;
  const int ceNP = ce_nP;
  DvceArray1D<Real> out("ck_selftest", 4+CK_NB);
  par_for("ck_selftest", DevExeSpace(), 0, 3, KOKKOS_LAMBDA(const int n) {
    const Real Tv = (n % 2 == 0) ? 800.0 : 2500.0;
    const int b = (n < 2) ? 0 : (CK_NB-1);      // 0 = reddest, CK_NB-1 = bluest
    int iT, iP;
    Real fT, fP;
    ck_tp_index(lT, nT, log10(Tv), iT, fT);
    ck_tp_index(lP, nP, log10(0.1), iP, fP);    // 0.1 bar
    Real km = 0.0;
    for (int g=0; g<CK_NG; ++g) {
      km += gw(g)*ck_kappa(lk, iT, fT, iP, fP, b, g);
    }
    out(n) = km;
    if (n == 0) {
      // continuum at the same reference point, with rho = 1 so the printed number is the
      // volumetric coefficient and can be checked against an independent parse
      Real kc[CK_NB];
      ck_continuum(ce, celT, celP, ceNT, ceNP, ciaN, ciaT, ciak, rayx, wl,
                   3500.0, 0.1, 1.0, kc);
      for (int b=0; b<CK_NB; ++b) out(4+b) = kc[b];
    }
  });
  auto h = Kokkos::create_mirror_view(out);
  Kokkos::deep_copy(h, out);
  if (global_variable::my_rank == 0) {
    std::cout << "  band-mean kappa at 0.1 bar [cm^2/g]:  reddest band  "
              << h(0) << " (800 K) -> " << h(1) << " (2500 K)" << std::endl
              << "                                       bluest band   "
              << h(2) << " (800 K) -> " << h(3) << " (2500 K)" << std::endl;
    std::cout << "  continuum at 3500 K, 0.1 bar, rho=1 [cm^-1], reddest to bluest:"
              << std::endl << "   ";
    for (int b=0; b<CK_NB; ++b) std::cout << " " << h(4+b);
    std::cout << std::endl;
    if (!(h(3) > 100.0*h(2) && h(0) > h(1))) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: the correlated-k band ordering "
                << "looks reversed." << std::endl
                << "  Expected the bluest band to climb steeply from 800 to 2500 K "
                << "(TiO/VO/Fe/Na/K leaving condensation)" << std::endl
                << "  and the reddest band (H2O rotational) to fall. See "
                << "data/exo_fms_ck/PROVENANCE.md." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ktab_lookup()
//  \brief bilinear interpolation in a synthetic correlated-k table. The cell coordinates
//  are computed once per cell by the caller, as a real scheme would compute log T and
//  log p once per cell rather than once per chain.

KOKKOS_INLINE_FUNCTION
Real ktab_lookup(const DvceArray3D<Real> &kt, const int c,
                 const Real &ct, const Real &cp) {
  const int nt = RT_KT_NT;
  const int np = RT_KT_NP;
  // Clamp the COORDINATE, not just the index. Clamping only the index leaves the
  // fraction as a huge extrapolation weight, which drives the interpolated opacity
  // negative, flips the sign of the optical depth increment, and makes the intensity
  // recurrence diverge to inf -- which then poisons the zero-weight sum as 0*inf = NaN.
  Real xt = (ct < 0.0) ? 0.0 : ((ct > nt-1.0) ? nt-1.0 : ct);
  Real xp = (cp < 0.0) ? 0.0 : ((cp > np-1.0) ? np-1.0 : cp);
  int it = static_cast<int>(xt);
  int ip = static_cast<int>(xp);
  it = (it > nt-2) ? nt-2 : it;
  ip = (ip > np-2) ? np-2 : ip;
  const Real ft = xt - static_cast<Real>(it);
  const Real fp = xp - static_cast<Real>(ip);
  const Real k00 = kt(c,it  ,ip  );
  const Real k01 = kt(c,it  ,ip+1);
  const Real k10 = kt(c,it+1,ip  );
  const Real k11 = kt(c,it+1,ip+1);
  return (1.0-ft)*((1.0-fp)*k00 + fp*k01) + ft*((1.0-fp)*k10 + fp*k11);
}

// problem/bc_outer_maxwell: whether the outer-x1 ghost extrapolation carries the
// divergence of the Maxwell stress.
//
// The extrapolation is hydrostatic, and this term was added so the ghost sees the
// magnetic force as well as gravity. It enters as e0 -= e_i*dM1mag/(rho_i*grav_acc), so
// its size relative to the hydrostatic term is (B^2/dr)/(rho g). At the outermost
// shell of the standard setup (rho = 4.0e-10, dr = 5.6e7, g = 942) that is 0.03 at 3 G,
// 0.33 at 10 G, 3.0 at 30 G and 33 at 100 G. Past O(1) it is no longer a correction: it
// swings the ghost density by its own size and the last two active radial shells go.
//
// It is now applied as an effective gravity inside the hydrostatic solve instead, and
// clamped, so it is bounded at any field strength -- see the use site. DEFAULT ON, since
// in that form it is the better physics; set false to drop the magnetic force entirely.
bool bc_outer_maxwell = true;
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  bool use_etotgrav = false;
  bool use_wellbalance_static = false;
  bool use_wellbalance_dynamic = false;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  // read before anything restart-sensitive: the outer BC needs it on restarts too
  bc_outer_maxwell = pin->GetOrAddBoolean("problem","bc_outer_maxwell",true);
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: outer-x1 Maxwell-stress term in the ghost "
              << "extrapolation is " << (bc_outer_maxwell ? "ON" : "off") << std::endl;
  }
  // blocked-band RT scaling harness: sized once, before any RT call
  rt_nchain = pin->GetOrAddInteger("problem","rt_nchain",4);
  rt_ktab = pin->GetOrAddBoolean("problem","rt_ktab",false);
  rt_split = pin->GetOrAddBoolean("problem","rt_split",false);
  rt_ck = pin->GetOrAddBoolean("problem","rt_ck",false);
  rt_ck_pcut = pin->GetOrAddReal("problem","ck_pcut_bar",10.0);
  ck_nq = pin->GetOrAddInteger("problem","ck_nquad",1);
  if (ck_nq != 1 && ck_nq != 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_nquad must be 1 "
              << "(diffusivity factor) or 2 (Gauss), got " << ck_nq << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (rt_ck && ck_lk_ptr == nullptr) {
    read_ck_table(pin->GetOrAddString("problem","ck_table",
                                      "data/exo_fms_ck/ck/Premixed_1x_g8_11.txt"));
    build_planck_fractions();
    read_ck_continuum(pin->GetOrAddString("problem","ck_data_dir",
                                          "data/exo_fms_ck"),
                      pin->GetOrAddString("problem","ck_swflux",
                                          "sw_band_flux_W121_11.txt"));
    ck_selftest();
    // The chain set is fixed by the table once correlated-k is on: every (band, g-point)
    // for each of the two Gauss points of the two-stream angular quadrature. Note this is
    // TWICE the 88 usually quoted for an 11 x 8 scheme -- 88 counts band x g, and this
    // kernel also carries the 2-point angular quadrature the picket fence uses.
    rt_nchain = CK_NB*CK_NG*ck_nq;
    if (global_variable::my_rank == 0) {
      std::cout << "  " << CK_NB << " bands x " << CK_NG << " g-points x " << ck_nq
                << " angular point(s) = " << rt_nchain << " column solves per cell";
      if (ck_nq == 1) {
        std::cout << " (diffusivity factor " << CK_DIFFUSIVITY << ")" << std::endl;
      } else {
        std::cout << " (2-point Gauss)" << std::endl;
      }
    }
  }
  if (rt_nchain < 4) rt_nchain = 4;
  if (rt_ktab_ptr == nullptr) {
    rt_ktab_ptr = new DvceArray3D<Real>("rt_ktab", rt_nchain, RT_KT_NT, RT_KT_NP);
    rt_wgt_ptr = new DvceArray1D<Real>("rt_wgt", rt_nchain);
    auto ktab = *rt_ktab_ptr;
    auto rtwgt = *rt_wgt_ptr;
    par_for("rt_ktab_init", DevExeSpace(), 0, rt_nchain-1, 0, RT_KT_NT-1, 0, RT_KT_NP-1,
    KOKKOS_LAMBDA(const int c, const int it, const int ip) {
      // order-unity, smooth in (T,p), different per chain
      ktab(c,it,ip) = 0.1 + 0.45*(1.0 + sin(0.11*c + 0.3*it)*cos(0.2*ip + 0.07*c));
    });
    const bool ckon = rt_ck;
    const int cknq = ck_nq;
    auto ckgw = (rt_ck) ? *ck_gw_ptr : DvceArray1D<Real>("dummy", 1);
    par_for("rt_wgt_init", DevExeSpace(), 0, rt_nchain-1,
    KOKKOS_LAMBDA(const int c) {
      rtwgt(c) = (c < 4) ? 1.0 : 0.0;
      if (ckon) rtwgt(c) = ckgw((c/cknq) % CK_NG);  // correlated-k: the g-point weight
    });
    if (global_variable::my_rank == 0 && rt_nchain > 4) {
      std::cout << "deep_hot_jupiter_rt: RT scaling harness active, " << rt_nchain
                << " chains in blocks of " << RT_NB << ", table lookup "
                << (rt_ktab ? "ON" : "off")
                << " (chains past 4 carry zero weight)" << std::endl;
    }
  }
  user_bcs_func = HydrostaticEquilibrium;
  // NOTE: this function must NOT return early on a restart. Only the initial condition
  // (the "probini" kernel) and the magnetic field at the end are genuinely one-off;
  // everything in between builds BACKGROUND state that lives in memory only and is
  // therefore gone after a restart -- the gravitational potential phicc0/phi0 consumed by
  // etotgrav and by the well-balanced reconstruction, and the well-balanced reference
  // atmosphere u0wb/w0wb/w0facewb. Returning here left phi identically zero, which
  // silently changed the physics at the restart point for any run with etotgrav = true.
  // The two one-off blocks are guarded with `if (!restart)` individually instead.
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pin->GetInteger("mesh", "nx1") != pin->GetInteger("meshblock", "nx1")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only allows one meshblock in r direction for the RT to work properly" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    
  Real r0, r1;
  r0 = pmy_mesh_->mesh_size.x1min;
  r1 = pmy_mesh_->mesh_size.x1max;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_;
  DvceArray5D<Real> w0_;
    
  DvceArray4D<Real> phi0_x1f;
  DvceArray4D<Real> phi0_x2f;
  DvceArray4D<Real> phi0_x3f;
  DvceArray4D<Real> phicc0;
  DvceArray5D<Real> u0wb;
  DvceArray5D<Real> w0wb;
  DvceArray5D<Real> w0facewb_x1f;
  DvceArray5D<Real> w0facewb_x2f;
  DvceArray5D<Real> w0facewb_x3f;
    
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &area3 = pmbp->pcoord->area.x3f;
  auto &dxe1 = pmbp->pcoord->dxedge.x1e;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;

  Real gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    w0_ = pmbp->phydro->w0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
    use_etotgrav = pmbp->phydro->use_etotgrav;
    use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
    phi0_x1f = pmbp->phydro->phi0.x1f;
    phi0_x2f = pmbp->phydro->phi0.x2f;
    phi0_x3f = pmbp->phydro->phi0.x3f;
    phicc0 = pmbp->phydro->phicc0;
    u0wb = pmbp->phydro->u0wb;
    w0wb = pmbp->phydro->w0wb;
    w0facewb_x1f = pmbp->phydro->w0facewb.x1f;
    w0facewb_x2f = pmbp->phydro->w0facewb.x2f;
    w0facewb_x3f = pmbp->phydro->w0facewb.x3f;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    w0_ = pmbp->pmhd->w0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
    use_etotgrav = pmbp->pmhd->use_etotgrav;
    use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
    phi0_x1f = pmbp->pmhd->phi0.x1f;
    phi0_x2f = pmbp->pmhd->phi0.x2f;
    phi0_x3f = pmbp->pmhd->phi0.x3f;
    phicc0 = pmbp->pmhd->phicc0;
    u0wb = pmbp->pmhd->u0wb;
    w0wb = pmbp->pmhd->w0wb;
    w0facewb_x1f = pmbp->pmhd->w0facewb.x1f;
    w0facewb_x2f = pmbp->pmhd->w0facewb.x2f;
    w0facewb_x3f = pmbp->pmhd->w0facewb.x3f;
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
  // by-value copy, capturable in the device lambdas below. Every pressure -> internal
  // energy conversion goes through it: for a general EOS e is not p/(gamma-1), and the
  // well-balanced background arrays built here feed the static scheme directly.
  auto eos = (pmbp->phydro != nullptr) ? pmbp->phydro->peos->eos_data
                                       : pmbp->pmhd->peos->eos_data;
    
//  Real Teq = 1469.0;
//  Real grav_acc = -942.0;
//  Real ap = 9.44e9;
//  Real omega = 2.06e-5;
//  Real Rgas = 4.593e7;
//  Real met = 0.0;
    
    Real Teq = pin->GetReal("problem","Teq");
    Real grav_acc = -pin->GetReal("problem","grav");
    Real ap = pin->GetReal("problem","ap");
    Real omega = pin->GetReal("problem","omega");
    Real Rgas = pin->GetReal("problem","Rgas");
    Real met = pin->GetReal("problem","met");

    // <problem>/Rgas is the ideal-gas R/mu that fixes this atmosphere's mean molecular
    // weight. Under a general EOS composition lives in the EOS instead, so Rgas is unused
    // and the Kelvin conversion is eos.temp_cgs = (pres_cgs/dens_cgs) m_u/k_B, which
    // carries NO mu. The two therefore agree only at mu = 1, i.e. Rgas*temp_cgs = 1. That
    // is not required for a physical run -- but it is required for the general path to
    // reproduce the ideal one, which is how this pgen is verified, so say so loudly.
    if (eos.IsGeneral() && global_variable::my_rank == 0) {
      Real mu_implied = 1.0/(Rgas*eos.temp_cgs);
      if (fabs(mu_implied - 1.0) > 1.0e-4) {
        std::cout << std::endl << "### WARNING! in " << __FILE__ << " at line "
                  << __LINE__ << std::endl
                  << "<problem>/Rgas = " << Rgas << " implies mean molecular weight "
                  << mu_implied << ", but the general EOS supplies its own composition "
                  << "and is being asked for temperature directly." << std::endl
                  << "Rgas is now unused; this run will NOT reproduce the eos=ideal run. "
                  << "Set Rgas = " << 1.0/eos.temp_cgs << " for that comparison."
                  << std::endl << std::endl;
      }
    }

    // <problem>/met is [M/H] in dex and feeds the opacity fit in get_kapr(); the EOS
    // scales its metal electron donors by eos_metal_mh, which DEFAULTS to met and so
    // normally agrees automatically. It can only differ if it was set explicitly, and the
    // result would be an atmosphere opaque at one metallicity and conducting at another,
    // with the electron fraction -- hence the Ohmic resistivity -- wrong by
    // 10^(met - eos_metal_mh). Refuse rather than let that pass.
    if (eos.IsGeneral() && eos.MetalIonization()) {
      if (fabs(met - eos.MetalMetallicity()) > 1.0e-6) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<problem>/met = " << met << " but the EOS metal donors were built "
                  << "with eos_metal_mh = " << eos.MetalMetallicity()
                  << "; both are [M/H] in dex and must agree. Remove the explicit "
                  << "eos_metal_mh and it will follow met." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    Real iap = 1.0/ap;
    
    Real grav = -grav_acc;
    Real Tirr = Teq*sqrt(2);
    Real Tint;
    get_Tint(Teq, Tint);
    Real mus = 1.0;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_wb_eos_arr(eos, Rgas, grav_acc, N, (r1-r0)*1.1, zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    
    DualArray1D<Real> Tarr_init("Tarrinit", N);
    DualArray1D<Real> lgparr_init("lgparrinit", N);
    get_picket_fence_pT_arr(eos, Rgas, gamma, Tint, Tirr, met, grav, mus, N, Tarr_init.h_view, lgparr_init.h_view);
    
    Tarr_init.modify_host();
    Tarr_init.sync_device();
    lgparr_init.modify_host();
    lgparr_init.sync_device();
    
//    DvceArray1D<Real> zinitarr("zinitarr", N);
//    DvceArray1D<Real> logpinitarr("logpinitarr", N);
    DualArray1D<Real> zarr_init("zarrinit", N);
    DualArray1D<Real> logparr_init("logparrinit", N);
    get_init_eos_arr(eos, Rgas, grav_acc, Tarr_init.h_view, lgparr_init.h_view, N, (r1-r0)*1.1, zarr_init.h_view, logparr_init.h_view);
    
    zarr_init.modify_host();
    zarr_init.sync_device();
    logparr_init.modify_host();
    logparr_init.sync_device();
  
    // one-off: the initial condition. Skipped on a restart, where u0/w0 come from file.
    if (!restart) {
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
        
      Real x1v, x2v, x3v;
      if (use_spherical_polar) {
        x1v = x1v_(m,i);
        x2v = x2v_(m,j);
        x3v = x3v_(m,k);
      } else {
        x1v = CellCenterX(i-is, nx1, x1min, x1max);
        x2v = CellCenterX(j-js, nx2, x2min, x2max);
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }
      Real r = x1v;
      if (use_spherical_polar) x1v -= ap;
        
      Real lam, theta, phi;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -x2v+M_PI/2.0;
        phi = x3v-M_PI;
      } else {
        theta = -(x2v*iap-M_PI/2.0);
        lam = x2v*iap;
        phi = x1v*iap;
      }
        
      Real pwb, denwb;
      get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
      Real p, den;
      get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
//      get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//      get_init_eos(zinitarr,logpinitarr,x3v,lam,phi,den,p);
//      p = pwb;
//      den = denwb;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
      w0_(m,IDN,k,j,i) = den;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
        if (use_etotgrav) {
            u0_(m,IEN,k,j,i) += den*phicc;
        }
    });
    }  // end of !restart guard on the initial condition

    par_for("probwb", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        }
        Real r = x1v;
        if (use_spherical_polar) x1v -= ap;
          
        Real lam, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -x2v+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          theta = -(x3v*iap-M_PI/2.0);
          lam = x3v*iap;
          phi = x2v*iap;
        }
          
        Real pwb, denwb;
        get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
        Real p, den;
        get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
//        get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//        get_init_eos(zinitarr,logpinitarr,x1v,lam,phi,den,p);
//        p = pwb;
//        den = denwb;
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
        
      if (use_etotgrav || use_wellbalance_dynamic) {
          if (use_spherical_polar) {
            x1v = x1f_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              if (use_spherical_polar) {
                x1v = x1f_(m,i+1);
              } else {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              }
              r = x1v;
              if (use_spherical_polar) x1v -= ap;
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2f_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x2f(m,k,j,i) = phicc;
          if (j == je) {
              if (use_spherical_polar) {
                x2v = x2f_(m,j+1);
              } else {
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3f_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              if (use_spherical_polar) {
                x3v = x3f_(m,k+1);
              } else {
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x3f(m,k+1,j,i) = phicc;
          }
      }
        if (use_wellbalance_static) {
            if (use_spherical_polar) {
              x1v = x1f_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            Real denwb;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (i == ie) {
                if (use_spherical_polar) {
                  x1v = x1f_(m,i+1);
                } else {
                  x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                }
                if (use_spherical_polar) x1v -= ap;
                get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2f_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (j == je) {
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3f_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (k == ke) {
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = EintFromP(eos, igm1, denwb, pwb);
            }
        }
    });
    if (use_etotgrav || use_wellbalance_dynamic) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            
            Real x1v, x2v, x3v;
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            Real r = x1v;
            if (use_spherical_polar) x1v -= ap;
              
            Real lam, theta, phi;
            if (use_spherical_polar) {
              theta = x2v;
              lam = -x2v+M_PI/2.0;
              phi = x3v-M_PI;
            } else {
              theta = -(x3v*iap-M_PI/2.0);
              lam = x3v*iap;
              phi = x2v*iap;
            }
            
            Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
            phicc0(m,k,j,i) = phicc;
        });
//        par_for("wbgravbc", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1,
//        KOKKOS_LAMBDA(int m, int k, int j) {
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect) {
//            for (int i=0; i<ng; ++i) {
//              phicc0(m,k,j,ie+i+1) = phicc0(m,k,j,ie-i);
//              u0_(m,IDN,k,j,ie+i+1) = u0_(m,IDN,k,j,ie-i);
//              u0_(m,IEN,k,j,ie+i+1) = u0_(m,IEN,k,j,ie-i);
//              phi0_x1f(m,k,j,ie+i+1) = phicc0(m,k,j,ie);
//            }
//          }
//        });
    }
    if (use_wellbalance_static) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          int nx1 = indcs.nx1;
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
            
          Real x1v, x2v, x3v;
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          Real r = x1v;
          if (use_spherical_polar) x1v -= ap;
            
          Real pwb, denwb;
          get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
          if (use_etotgrav) {
              Real phicc = - grav_acc * x1v;
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
        });
    }

    // initialize magnetic fields if MHD. One-off like the initial condition above: on a
    // restart b0/bcc0 are read from file.
    if (!restart && pmbp->pmhd != nullptr) {
      // Read magnetic field strength
      Real bbot = pin->GetReal("problem","bbot");
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),0, n3m1, 0, n2m1, 0, n1m1, //ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
          if (use_spherical_polar) {
            Real x1v = x1v_(m,i);
            Real x2v = x2v_(m,j);
            Real x3v = x3v_(m,k);
            Real x1fl = x1f_(m,i);
            Real x2fl = x2f_(m,j);
            Real x3fl = x3f_(m,k);
            Real x1fr = x1f_(m,i+1);
            Real x2fr = x2f_(m,j+1);
            Real x3fr = x3f_(m,k+1);
              
            Real A1 = 0.0;
            Real A2 = 0.0;
            Real A2ip = 0.0;
            Real A2kp = 0.0;
            Real A2ipkp = 0.0;
            Real A3 = 0.5*bbot*r0*sin(x2fl)/SQR(x1fl/r0);
            Real A3ip = 0.5*bbot*r0*sin(x2fl)/SQR(x1fr/r0);
            Real A3jp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fl/r0);
            Real A3ipjp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fr/r0);
//            Real A3 = 0.5*bbot*x1fl*sin(x2fl);
//            Real A3ip = 0.5*bbot*x1fr*sin(x2fl);
//            Real A3jp = 0.5*bbot*x1fl*sin(x2fr);
//            Real A3ipjp = 0.5*bbot*x1fr*sin(x2fr);
              
            b0.x1f(m,k,j,i) = (dxe3(m,k,j+1,i)*A3jp - dxe3(m,k,j,i)*A3)/area1(m,k,j,i) - (dxe2(m,k+1,j,i)*A2kp - dxe2(m,k,j,i)*A2)/area1(m,k,j,i);
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1)) {
              b0.x2f(m,k,j,i) = - ((x1fr*(x3fr-x3fl))*A3ip - (x1fl*(x3fr-x3fl))*A3) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0.x2f(m,k,j,i) = - (dxe3(m,k,j,i+1)*A3ip - dxe3(m,k,j,i)*A3)/area2(m,k,j,i);
            }
            b0.x3f(m,k,j,i) = (dxe2(m,k,j,i+1)*A2ip - dxe2(m,k,j,i)*A2)/area3(m,k,j,i);
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j,i+1)*A3ip)/area1(m,k,j,i+1) - (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k,j,i+1)*A2ip)/area1(m,k,j,i+1);
            Real b0x2fjp;
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js-1) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je)) {
              b0x2fjp = - ((x1fr*(x3fr-x3fl))*A3ipjp - (x1fl*(x3fr-x3fl))*A3jp) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0x2fjp = - (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j+1,i)*A3jp)/area2(m,k,j+1,i);
            }
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k+1,j,i)*A2kp)/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;
              
            Real lw, rw;
            lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0x2fjp;
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0x3fkp;
              
            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))+SQR(bcc0(m,IBY,k,j,i))+SQR(bcc0(m,IBZ,k,j,i)));
          }
      });
    }

    // On a restart u0 and b0 are restored exactly, ghost zones included, but bcc0 is a
    // DERIVED array that the restart file does not carry, so it starts at zero. The
    // inner-x1 user boundary swaps the magnetic energy of its ghost cells by subtracting
    // 0.5*bcc0^2 and adding it back from the current face fields; with bcc0 still zero it
    // subtracts nothing and so double counts the magnetic energy on the very first step.
    // Fill bcc0 from the restored face fields here, which is exactly the state the
    // initialisation above leaves on a fresh start. The interpolation must match the one
    // ConsToPrim uses, hence the split on use_spherical_polar.
    if (restart && pmbp->pmhd != nullptr) {
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_bcc_restart", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
              0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        if (use_spherical_polar) {
          Real lw, rw;
          lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0.x2f(m,k,j+1,i);
          lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0.x3f(m,k+1,j,i);
        } else {
          bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
          bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        }
      });
    }

  return;
}


void HydrostaticEquilibrium(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  int nmb = pm->pmb_pack->nmb_thispack;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

    // wtemp holds the temperature ConsToPrim solved for; the outer-x1 ghost extrapolation
    // warm starts its hydrostatic solve from the last active cell's value. General EOS
    // only -- a zero-size View otherwise, captured and never read.
    DvceArray4D<Real> wtemp_;
    DvceArray5D<Real> u0_;
    DvceArray5D<Real> w0_;
    Real gamma;
    bool use_etotgrav = false;
    bool use_wellbalance_dynamic = false;
    
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    DvceArray5D<Real> bcc0;
    DvceArray4D<Real> b0_x1f;
    DvceArray4D<Real> b0_x2f;
    DvceArray4D<Real> b0_x3f;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &area2 = pmbp->pcoord->area.x2f;
    auto &area3 = pmbp->pcoord->area.x3f;
    auto &volume = pmbp->pcoord->volume;
    auto &z_ov_rE = pmbp->pcoord->z_ov_rE;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;

    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      bcc0 = pmbp->pmhd->bcc0;
      b0_x1f = pmbp->pmhd->b0.x1f;
      b0_x2f = pmbp->pmhd->b0.x2f;
      b0_x3f = pmbp->pmhd->b0.x3f;
    }
    Real bbot = pm->pgen->hot_jupiter_param.bbot;
    
    int nvar = u0_.extent_int(1);
    
    Real igm1 = 1.0/(gamma-1.0);
    Real gigm1 = gamma*igm1;
    Real gm1ig = (gamma-1.0)/gamma;
    Real ig = 1.0/gamma;
    
    if (pmbp->pmhd != nullptr) {
      par_for("usrboundaryx1_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = is+i;
              int iex = is-i-1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin+1)/area1(m,k,j,iex);
              b0_x1f(m,k,j,iex) = b0_x1f(m,k,j,iin+1)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = is-i-1;
//              Real div_rest = b0_x1f(m,k,j,iex+1)*area1(m,k,j,iex+1) + (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) + (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex) = div_rest/area1(m,k,j,iex);
//            }
          }
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = ie-i;
              int iex = ie+i+1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin)/area1(m,k,j,iex+1);
              b0_x1f(m,k,j,iex+1) = b0_x1f(m,k,j,iin)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = ie+i+1;
//              Real div_rest = b0_x1f(m,k,j,iex)*area1(m,k,j,iex) - (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) - (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex+1) = div_rest/area1(m,k,j,iex+1);
//            }
//              for (int i=0; i<ng; ++i) {
//                b0_x1f(m,k,j,ie+i+2) = -b0_x1f(m,k,j,ie-i);
//                b0_x2f(m,k,j,ie+i+1) =  b0_x2f(m,k,j,ie-i);
//                if (j == n2-1) {b0_x2f(m,k,j+1,ie+i+1) = b0_x2f(m,k,j+1,ie-i);}
//                b0_x3f(m,k,j,ie+i+1) =  b0_x3f(m,k,j,ie-i);
//                if (k == n3-1) {b0_x3f(m,k+1,j,ie+i+1) = b0_x3f(m,k+1,j,ie-i);}
//              }
        }
      });
    }
    
//    if (pmbp->phydro != nullptr) {
//      if (use_etotgrav) {
//        pmbp->phydro->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->phydro->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
//    else if (pmbp->pmhd != nullptr) {
//      auto b0 = pmbp->pmhd->b0;
//      if (use_etotgrav) {
//        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->pmhd->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
    if (pmbp->phydro != nullptr) {
      if (use_etotgrav) {
        pmbp->phydro->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->phydro->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }
    else if (pmbp->pmhd != nullptr) {
      auto b0 = pmbp->pmhd->b0;
      if (use_etotgrav) {
        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->pmhd->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }

//    par_for("usrboundaryx1", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n2-1),
//    KOKKOS_LAMBDA(int m, int n, int k, int j) {
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
//          for (int i=0; i<ng; ++i) {
//            u0_(m,n,k,j,ie+i+1) = u0_(m,n,k,j,ie);
//          }
//        }
//    });
    
    // Local copy of the file-scope flag. Reading the global directly from the kernel is
    // a reference to a __host__ variable in device code, which hipcc rejects outright --
    // the switch has to be captured by value like any other host state.
    const bool bc_outer_maxwell_ = bc_outer_maxwell;

    // The cell-centred field in the OUTER-x1 ghost zones is built here, in its own kernel,
    // and not in the extrapolation kernel below.  It used to be computed inline there, in
    // the same launch that reads bcc0 at (k,j+1) and (k+1,j) for the Maxwell stress -- cells
    // owned by OTHER threads of that launch.  Whether a thread saw its neighbour's new value
    // or the previous cycle's stale one depended on how the wavefronts happened to be
    // scheduled, which made the whole run non-reproducible on a GPU: two runs of the same
    // binary on one rank diverged within a single cycle, worst in the outer half of the
    // domain, where the two-stream RT then spread it along each radial column.  On a CPU the
    // ascending j loop always lost the race the same way, so this was deterministic (and
    // deterministically wrong) until the problem first ran on an accelerator.
    //
    // Splitting the write from the read is the whole fix: every bcc0 the stress term reads
    // is now a completed write from a previous kernel.
    //
    // The INNER-x1 ghosts are deliberately NOT hoisted: that branch subtracts the magnetic
    // energy using the OLD bcc0 and adds it back with the new one, so moving the update
    // ahead of it would change what it subtracts.  It has no cross-thread read and was
    // measured deterministic on its own.
    par_for("usrboundaryx1_bcc_outer", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user &&
            pmbp->pmhd != nullptr) {
          for (int i=0; i<ng; ++i) {
            Real lw, rw;
            lw = (x1f_(m,(ie+i+1)+1)-x1v_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            rw = (x1v_(m,(ie+i+1))-x1f_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            bcc0(m,IBX,k,j,(ie+i+1)) = lw*b0_x1f(m,k,j,(ie+i+1))
                                     + rw*b0_x1f(m,k,j,(ie+i+1)+1);
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,(ie+i+1)) = lw*b0_x2f(m,k,j,(ie+i+1))
                                     + rw*b0_x2f(m,k,j+1,(ie+i+1));
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,(ie+i+1)) = lw*b0_x3f(m,k,j,(ie+i+1))
                                     + rw*b0_x3f(m,k+1,j,(ie+i+1));
          }
        }
    });

    par_for("usrboundaryx1_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real rho_i = u0_(m,IDN,k,j,is);
//          Real e_i = w0_(m,IEN,k,j,is);
//          Real phi_i = phicc0(m,k,j,is);
//          Real q0_i = log(e_i);
//          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            if (pmbp->pmhd != nullptr) {
              u0_(m,IEN,k,j,is-i-1) -= 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
              Real lw, rw;
              lw = (x1f_(m,(is-i-1)+1)-x1v_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              rw = (x1v_(m,(is-i-1))-x1f_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              bcc0(m,IBX,k,j,(is-i-1)) = lw*b0_x1f(m,k,j,(is-i-1)) + rw*b0_x1f(m,k,j,(is-i-1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(is-i-1)) = lw*b0_x2f(m,k,j,(is-i-1)) + rw*b0_x2f(m,k,j+1,(is-i-1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(is-i-1)) = lw*b0_x3f(m,k,j,(is-i-1)) + rw*b0_x3f(m,k+1,j,(is-i-1));
              u0_(m,IEN,k,j,is-i-1) += 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
            }
//              Real rho0_ip = u0_(m,IDN,k,j,(is-i-1));
//              u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM1,k,j,(is-i-1)) = 0.0;//-u0_(m,IM1,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IEN,k,j,(is-i-1)) = w0_(m,IEN,k,j,(is-i-1)) + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
////            Real dphi_i = phicc0(m,k,j,(is-i-1))-phi_i;
////            Real q0_ip = q0_i - factor_i * dphi_i;
////            Real e0_ip = exp(q0_ip);
////            if (e0_ip < 0.0) e0_ip = e_i;
////            Real rho0_ip = e0_ip/e_i*rho_i;
////            u0_(m,IDN,k,j,(is-i-1)) = rho0_ip;
////            u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
////            u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
////            Real mom = u0_(m,IM1,k,j,is)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,is)/rho_i*rho0_ip); //
////            u0_(m,IM1,k,j,(is-i-1)) = mom;
////            u0_(m,IEN,k,j,(is-i-1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
//            if (use_etotgrav) u0_(m,IEN,k,j,(is-i-1)) += rho0_ip*phicc0(m,k,j,(is-i-1));
//            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(is-i-1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(is-i-1)))+SQR(bcc0(m,IBY,k,j,(is-i-1)))+SQR(bcc0(m,IBZ,k,j,(is-i-1))));
          }
        }
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
          Real rho_i = w0_(m,IDN,k,j,ie);
//          Real e_i = u0_(m,IEN,k,j,ie) - 0.5*(SQR(u0_(m,IM1,k,j,ie))+SQR(u0_(m,IM2,k,j,ie))+SQR(u0_(m,IM3,k,j,ie)))/rho_i;
//          if (use_etotgrav) e_i -= rho_i*phicc0(m,k,j,ie);
//          if (pmbp->pmhd != nullptr) e_i -= 0.5*(SQR(bcc0(m,IBX,k,j,ie))+SQR(bcc0(m,IBY,k,j,ie))+SQR(bcc0(m,IBZ,k,j,ie)));
          Real e_i = w0_(m,IEN,k,j,ie);
          Real phi_i = phicc0(m,k,j,ie);
          Real q0_i = log(e_i);
          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            Real dM1mag = 0.0;
            if (pmbp->pmhd != nullptr) {
              // bcc0 in these ghost cells was built by usrboundaryx1_bcc_outer above.
              if (bc_outer_maxwell_) {
              // (k,j+1) and (k+1,j) are read here. They belong to other threads, which is
              // why the write had to move to its own kernel -- see the note on that kernel.
              // At the outermost ghost row j+1 and k+1 leave the cell-centred array (the
              // FACE arrays have the extra slot, the cell-centred one does not), so the
              // index is held at the edge: the correction degenerates to a one-sided
              // difference in that row rather than reading past the end, which is what the
              // previous version did.
              int jp1 = (j+1 < n2) ? (j+1) : j;
              int kp1 = (k+1 < n3) ? (k+1) : k;
              // Same one past the end in x1: on the LAST ghost cell (i = ng-1) the cell
              // centred index ie+i+2 is n1, one beyond the array. The face array b0_x1f has
              // the extra slot and is indexed as before; only bcc0 is held back.
              int ip2 = (ie+i+2 < n1) ? (ie+i+2) : (n1-1);
              Real pb = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
              Real pbp1 = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+2)))+SQR(bcc0(m,IBY,k,j,ip2))+SQR(bcc0(m,IBZ,k,j,ip2)));
              Real M11 = pb - SQR(b0_x1f(m,k,j,(ie+i+1)));
              Real M11p1 = pbp1 - SQR(b0_x1f(m,k,j,(ie+i+2)));
              Real M12 = - b0_x2f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M12p1 = - b0_x2f(m,k,j+1,(ie+i+1)) * bcc0(m,IBX,k,jp1,(ie+i+1));
              Real M13 = - b0_x3f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M13p1 = - b0_x3f(m,k+1,j,(ie+i+1)) * bcc0(m,IBX,kp1,j,(ie+i+1));
              dM1mag = -( (M11p1*area1(m,k,j,(ie+i+2))-M11*area1(m,k,j,(ie+i+1))) + (M12p1*area2(m,k,j+1,(ie+i+1))-M12*area2(m,k,j,(ie+i+1))) + (M13p1*area3(m,k+1,j,(ie+i+1))-M13*area3(m,k,j,(ie+i+1))) )/volume(m,k,j,(ie+i+1));
              dM1mag += z_ov_rE(m,k,j,(ie+i+1)) * 0.5*SQR(bcc0(m,IBX,k,j,(ie+i+1)));
              }
            }
            // Hydrostatic extrapolation into the ghost zone, at fixed temperature. The
            // closed form below is the ideal-gas isothermal background; for a general EOS
            // the same statement is WBAdvance's isothermal branch, which integrates
            // dln d/dPhi = -d/(p chi_rho) instead of assuming p = (gamma-1)e.
            Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
            // The magnetic force enters as an EFFECTIVE GRAVITY, not as an additive
            // shift of the extrapolated energy. dM1mag/rho is an acceleration, so
            // dphi -> dphi*(1 - dM1mag/(rho |g|)) is the same to first order but stays
            // inside the exponential (and inside WBAdvance's EOS-consistent
            // integration), so a magnetic force comparable to gravity changes the SCALE
            // HEIGHT instead of swinging the answer linearly through zero.
            //
            // The clamp is what makes it usable at high field. gmag -> 1 is the
            // force-free limit where the atmosphere stops falling off, and gmag > 1 is
            // net outward, which an outward-decaying ghost cannot represent; gmag < -1
            // would compress the ghost without bound. Both ends are held back, so the
            // ghost degrades to "very extended" rather than to nonsense.
            if (bc_outer_maxwell_) {
              Real gmag = dM1mag/(rho_i*fabs(grav_acc));
              gmag = fmin(fmax(gmag, -1.0), 0.9);
              dphi_i *= (1.0 - gmag);
            }
            Real e0_hyd, rho0_hyd;
            if (eos.IsGeneral()) {
              rho0_hyd = rho_i;
              e0_hyd = e_i;
              Real t_hyd = -1.0;   // WBAdvance's temperature hand-off; unused here
              // Warm start from cell ie's own temperature, which the ConsToPrim call on
              // this column above has already solved for and left in wtemp. This runs per
              // ghost cell per stage, and the isothermal branch inverts twice.
              WBAdvance(eos, 1, rho_i, e_i, dphi_i, rho0_hyd, e0_hyd, t_hyd,
                        wtemp_(m,k,j,ie));
            } else {
              e0_hyd = exp(q0_i - factor_i * dphi_i);
              rho0_hyd = e0_hyd/e_i*rho_i;
            }
            Real e0_ip = e0_hyd;
            if (e0_ip < 0.0) e0_ip = e_i;
            // density and energy now come from the SAME hydrostatic solve, so they are
            // thermodynamically consistent by construction; the old code had to rescale
            // the density by the energy's relative shift because the magnetic term was
            // bolted on afterwards
            Real rho0_ip = rho0_hyd;
//            rho0_ip = rho_i;
//            e0_ip = e_i;
            u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
            Real mom = u0_(m,IM1,k,j,ie)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,ie)/rho_i*rho0_ip); // 
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(ie+i+1)))+SQR(u0_(m,IM2,k,j,(ie+i+1)))+SQR(u0_(m,IM3,k,j,(ie+i+1))))/rho0_ip;
            if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(ie+i+1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
          }
        }
    });

////    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n1-1),
////    KOKKOS_LAMBDA(int m, int n, int k, int i) {
////        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,js-j-1,i) = -u0_(m,n,k,js+j,i);
////            } else {
////              u0_(m,n,k,js-j-1,i) = u0_(m,n,k,js+j,i);
////            }
////          }
////        }
////        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,je+j+1,i) = -u0_(m,n,k,je-j,i);
////            } else {
////              u0_(m,n,k,je+j+1,i) = u0_(m,n,k,je-j,i);
////            }
////          }
////        }
////    });
//    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
//    KOKKOS_LAMBDA(int m, int k, int i) {
//        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,js-j-1,i) = u0_(m,IDN,k,js,i);
//            u0_(m,IM1,k,js-j-1,i) = u0_(m,IM1,k,js,i);
//            u0_(m,IM3,k,js-j-1,i) = u0_(m,IM3,k,js,i);
//            u0_(m,IM2,k,js-j-1,i) = fmin(0.0,u0_(m,IM2,k,js,i));
//            u0_(m,IEN,k,js-j-1,i) = u0_(m,IEN,k,js,i)-0.5*SQR(u0_(m,IM2,k,js,i))+0.5*SQR(u0_(m,IM2,k,js-j-1,i));
//          }
//        }
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,je+j+1,i) = u0_(m,IDN,k,je,i);
//            u0_(m,IM1,k,je+j+1,i) = u0_(m,IM1,k,je,i);
//            u0_(m,IM3,k,je+j+1,i) = u0_(m,IM3,k,je,i);
//            u0_(m,IM2,k,je+j+1,i) = fmax(0.0,u0_(m,IM2,k,je,i));
//            u0_(m,IEN,k,je+j+1,i) = u0_(m,IEN,k,je,i)-0.5*SQR(u0_(m,IM2,k,je,i))+0.5*SQR(u0_(m,IM2,k,je+j+1,i));
//          }
//        }
//    });
//    if (pmbp->pmhd != nullptr) {
//      par_for("usrboundaryx2_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//      KOKKOS_LAMBDA(int m, int k, int i) {
//          if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = js+j;
//                int jex = js-j-1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin+1,i)/area2(m,k,jex,i);
//                b0_x2f(m,k,jex,i) = b0_x2f(m,k,jin+1,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = js-j-1;
////                Real div_rest = b0_x2f(m,k,jex+1,i)*area2(m,k,jex+1,i) + (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) + (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex,i) = div_rest/area2(m,k,jex,i);
////              }
//          }
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = je-j;
//                int jex = je+j+1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin,i)/area2(m,k,jex+1,i);
//                b0_x2f(m,k,jex+1,i) = b0_x2f(m,k,jin,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = je+j+1;
////                Real div_rest = b0_x2f(m,k,jex,i)*area2(m,k,jex,i) - (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) - (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex+1,i) = div_rest/area2(m,k,jex+1,i);
////              }
//          }
//      });
//        par_for("usrboundaryx2_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//        KOKKOS_LAMBDA(int m, int k, int i) {
//            if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js,i) = lw*b0_x1f(m,k,js,i) + rw*b0_x1f(m,k,js,i+1);
//                lw = (x2f_(m,js+1)-x2v_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                rw = (x2v_(m,js)-x2f_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                bcc0(m,IBY,k,js,i) = lw*b0_x2f(m,k,js,i) + rw*b0_x2f(m,k,js+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js,i) = lw*b0_x3f(m,k,js,i) + rw*b0_x3f(m,k+1,js,i);
//
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,je,i) = lw*b0_x1f(m,k,je,i) + rw*b0_x1f(m,k,je,i+1);
//                lw = (x2f_(m,je+1)-x2v_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                rw = (x2v_(m,je)-x2f_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                bcc0(m,IBY,k,je,i) = lw*b0_x2f(m,k,je,i) + rw*b0_x2f(m,k,je+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,je,i) = lw*b0_x3f(m,k,je,i) + rw*b0_x3f(m,k+1,je,i);
//
//              for (int j=0; j<ng; ++j) {
////                Real lw, rw;
////                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                bcc0(m,IBX,k,js+j,i) = lw*b0_x1f(m,k,js+j,i) + rw*b0_x1f(m,k,js+j,i+1);
////                lw = (x2f_(m,js+j+1)-x2v_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                rw = (x2v_(m,js+j)-x2f_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                bcc0(m,IBY,k,js+j,i) = lw*b0_x2f(m,k,js+j,i) + rw*b0_x2f(m,k,js+j+1,i);
////                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                bcc0(m,IBZ,k,js+j,i) = lw*b0_x3f(m,k,js+j,i) + rw*b0_x3f(m,k+1,js+j,i);
//                u0_(m,IEN,k,js-j-1,i) -= 0.5*(SQR(bcc0(m,IBX,k,js,i))+SQR(bcc0(m,IBY,k,js,i))+SQR(bcc0(m,IBZ,k,js,i)));
//
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js-j-1,i) = lw*b0_x1f(m,k,js-j-1,i) + rw*b0_x1f(m,k,js-j-1,i+1);
//                lw = (x2f_(m,js-j-1+1)-x2v_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                rw = (x2v_(m,js-j-1)-x2f_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                bcc0(m,IBY,k,js-j-1,i) = lw*b0_x2f(m,k,js-j-1,i) + rw*b0_x2f(m,k,js-j-1+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js-j-1,i) = lw*b0_x3f(m,k,js-j-1,i) + rw*b0_x3f(m,k+1,js-j-1,i);
//                u0_(m,IEN,k,js-j-1,i) += 0.5*(SQR(bcc0(m,IBX,k,js-j-1,i))+SQR(bcc0(m,IBY,k,js-j-1,i))+SQR(bcc0(m,IBZ,k,js-j-1,i)));
//              }
//            }
//            if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
////                  Real lw, rw;
////                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  bcc0(m,IBX,k,je-j,i) = lw*b0_x1f(m,k,je-j,i) + rw*b0_x1f(m,k,je-j,i+1);
////                  lw = (x2f_(m,je-j+1)-x2v_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  rw = (x2v_(m,je-j)-x2f_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  bcc0(m,IBY,k,je-j,i) = lw*b0_x2f(m,k,je-j,i) + rw*b0_x2f(m,k,je-j+1,i);
////                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  bcc0(m,IBZ,k,je-j,i) = lw*b0_x3f(m,k,je-j,i) + rw*b0_x3f(m,k+1,je-j,i);
//                  u0_(m,IEN,k,je+j+1,i) -= 0.5*(SQR(bcc0(m,IBX,k,je,i))+SQR(bcc0(m,IBY,k,je,i))+SQR(bcc0(m,IBZ,k,je,i)));
//
//                  Real lw, rw;
//                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  bcc0(m,IBX,k,je+j+1,i) = lw*b0_x1f(m,k,je+j+1,i) + rw*b0_x1f(m,k,je+j+1,i+1);
//                  lw = (x2f_(m,je+j+1+1)-x2v_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  rw = (x2v_(m,je+j+1)-x2f_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  bcc0(m,IBY,k,je+j+1,i) = lw*b0_x2f(m,k,je+j+1,i) + rw*b0_x2f(m,k,je+j+1+1,i);
//                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  bcc0(m,IBZ,k,je+j+1,i) = lw*b0_x3f(m,k,je+j+1,i) + rw*b0_x3f(m,k+1,je+j+1,i);
//                  u0_(m,IEN,k,je+j+1,i) += 0.5*(SQR(bcc0(m,IBX,k,je+j+1,i))+SQR(bcc0(m,IBY,k,je+j+1,i))+SQR(bcc0(m,IBZ,k,je+j+1,i)));
//              }
//            }
//        });
//    }
  return;
}


void SourceFunc(Mesh *pm, Real bdt) {
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

    DvceArray5D<Real> u0, w0, w0wb;
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    bool use_etotgrav = false;
    bool use_wellbalance_static = false;
    bool use_wellbalance_dynamic = false;
    Real gamma;
    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
      w0wb = pmbp->phydro->w0wb;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      w0wb = pmbp->pmhd->w0wb;
    }
    
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    
//    Real ap = 9.44e9;
//    Real omega = 2.06e-5;
//    Real grav_acc = -942.0;
//    Real Rgas = 4.593e7;
    
//    ParameterInput* pin;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    
    Real time = pm->time;
    
    picket_fence_two_stream_RT(pm, bdt);

    par_for("usrsource", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
        
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        }
        
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
          z = x1v-ap;
          r = x1v;
        } else {
          lam = x3v*iap;
          phi = x2v*iap;
          z = x1v;
        }
        Real rho = w0(m,IDN,k,j,i);
        Real p = eos.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,i),p);
        
        Real area_r = area1(m,k,j,i+1);
        Real area_l = area1(m,k,j,i);
        Real vol = volume(m,k,j,i);
        
        // gravity
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance_static) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_dynamic) {
          // The source is the background's own pressure difference across the cell, so it
          // has to come from the same entry point the reconstruction uses: getWBq0, which
          // dispatches to the ideal-gas closed forms or to the general-EOS background and
          // returns pressure directly, not an energy to be multiplied by (gamma-1).
          Real pl,pr,dum1,dum2,dum3;
          if (pmbp->phydro != nullptr) {
              pmbp->phydro->getWBq0(eos, WBVar::wb_pres,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,pl,dum2,pr,dum3);
          } else if (pmbp->pmhd != nullptr) {
              pmbp->pmhd->getWBq0(eos, WBVar::wb_pres,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,pl,dum2,pr,dum3);
          }
          src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol;
        }
        u0(m,IM1,k,j,i) += src;
        
        // Forces in the corotating frame
        if (use_spherical_polar) {
          Real vtheta = w0(m,IVY,k,j,i);
          Real vphi = w0(m,IVZ,k,j,i);
          Real vr = w0(m,IVX,k,j,i);
          Real sine = sin(theta);
          Real cosine = cos(theta);
          Real oor = SQR(omega)*r*sine;
          Real cor = 2.0*omega*vphi;
          u0(m,IM2,k,j,i) += rho*(cor+oor)*cosine*bdt;
          u0(m,IM3,k,j,i) += -rho*2.0*omega*(vr*sine+vtheta*cosine)*bdt;
          u0(m,IM1,k,j,i) += rho*(cor+oor)*sine*bdt;
//          if (!use_etotgrav)
          u0(m,IEN,k,j,i) += rho*oor*(vr*sine+vtheta*cosine)*bdt;
        } else {
          // corotating beta-plane approximation e.g. Fromang+2016
          Real omega1 = omega*lam;
          u0(m,IM2,k,j,i) += -2.0*rho*omega1*(-w0(m,IVZ,k,j,i))*bdt;
          u0(m,IM3,k,j,i) += -2.0*rho*omega1*w0(m,IVY,k,j,i)*bdt;
        }

//        // Newtonian cooling
//        Real Teq, itrad;
//        get_eq_Tp(lam, phi, p, Teq);
//        get_itrad(p,itrad);
////        Real t0 = 2.16e6;
////        Real ff = (t0-time)/t0;
////        ff = (ff < 0.0) ? 0.0 : ff;
////        Real gg = 2.0*ff;
////        itrad /= pow(10.0,gg);
//        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
//        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;

        // Rayleigh drag (initial relaxation)
        Real dyntime = 2.0*M_PI/omega;
        Real tau1 = dyntime / 10.0;
        Real tau2 = dyntime;
        Real t1   = 2.0 * dyntime;
        Real t2   = 5.0 * dyntime;
        Real itdrag, fredux;
        if(time < t1) {
          itdrag = 1.0 / tau1;
        } else if(time < t2)
          {
            Real alp = (time - t1) / (t2 - t1);
            itdrag   = 1.0 / tau1 * pow(tau1 / tau2, alp);
          }
        else {
          itdrag = 0.0;
        }
        if (time < t2) {
          fredux = itdrag*bdt;///(1.0+itdrag*bdt);
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }

        // Top sponge layer
        Real bar = 1.0e6;
        Real logpl = log(1.0e-6*bar);
        Real logpt = log(1.0e-7*bar);
        Real logp = log(p);
        Real fdrag = 1.0 - (logp-logpt)/(logpl-logpt); // high p = 0, low p = 1
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        // Bottom sponge layer
        logpl = log(1.0e2*bar);
        logpt = log(5.0e1*bar);
        logp = log(p);
        fdrag = (logp-logpt)/(logpl-logpt); // high p = 1, low p = 0
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
//        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
    });

    return;
}

void double_gray_two_stream_RT(Mesh *pm, Real bdt) {
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
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 100.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 5;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
        constexpr int NN = 270;
        Real tau_ir_down_f[NN];
        Real F_v_down_f[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real B[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real rtop = x1f_(m,ie+1);
        Real rbot = x1f_(m,is);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = -1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real pm1 = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie),w0(m,IEN,k,j,ie));
        Real pf = exp((log(p)+log(pm1))/2.0);
        Real tau_v_f = 0.0;//pf/(grav/kap_v);
        Real tau_ir_f = 0.0;//pf/(grav/kap_ir);
        
        tau_ir_down_f[ie+1] = tau_ir_f;
        F_v_down_f[ie+1] = Fstar*mu0;//*exp(-tau_v_f/muf);
        F_v_down_f[ie+1] = (mu0 > 0.0) ? F_v_down_f[ie+1] : 0.0;
        
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
          Real kap_ir = 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
            
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);

//          Real muf = sqrt(1.0-SQR(ap/rf)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r1));
          Real muf = (mu0 < mucr) ? mucr : mu0;
          if (test_oned) muf = mu0;
          Real dtau_v = kap_v*rho*dr;
          Real dtau_ir = kap_ir*rho*dr;
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf); // Zhang+2023
          Real trans = exp(-dtau_v/muf);
          F_v_down_f[i] = F_v_down_f[i+1]*trans*fac;
          tau_ir_down_f[i] = tau_ir_down_f[i+1] + dtau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        I_ir_down_f[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf);
          Real dtau = tau_ir_down_f[i]-tau_ir_down_f[i+1];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i]+B[i+1])/2.0;
          I_ir_down_f[i] = (I_ir_down_f[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        I_ir_up_f[is] = Iint;
        // up-sweep
        for (int i=is+1; i<ie+2; ++i) {
          Real rf = x1f_(m,i);
          Real rfm1 = x1f_(m,i-1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rfm1/rf);
          Real dtau = tau_ir_down_f[i-1]-tau_ir_down_f[i];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i-1]+B[i])/2.0;
          I_ir_up_f[i] = (I_ir_up_f[i-1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        // flux divergence
        for (int i=is; i<ie+1; ++i) {
          Real Ft = 2.0*M_PI*mu1*(I_ir_up_f[i+1]-I_ir_down_f[i+1])-F_v_down_f[i+1];
          Real Fb = 2.0*M_PI*mu1*(I_ir_up_f[i]-I_ir_down_f[i])-F_v_down_f[i];
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) src = -(Ft*area_t-Fb*area_b)/vol;
          u0(m,IEN,k,j,i) += src*bdt;
        }
        
    });
    
    return;
}

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt) {
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
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real grav = 942.0;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 500.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 4;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(0), n1);
        constexpr int NN = 270;
        Real tau_ir_down[NN];
        Real I_down[NN];
        Real B[NN];
        Real Q_v[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(50.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real tau_v = 0.0;//p/(grav/kap_v);
        Real tau_ir = 0.0;//p/(grav/kap_ir);
//        if (test_oned) {
//          tau_v = p/(grav/kap_v)/mu0;
//          tau_ir = p/(grav/kap_ir)/mu1;
//        }
        tau_ir_down[ie+1] = tau_ir;
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
            Real kap_ir = 1.0e-2; // 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
            
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
////          Real mu = sqrt(1.0-SQR(r1/r)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r));
//          Real mu = (mu0 < mucr) ? mucr : mu0;
//          if (test_oned) mu = mu0;
//          Real dtau_v = kap_v*rho*dr;
          Real delta = dr/r;
          Real drcor = r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          if (test_oned) drcor = dr/mu0;
          Real dtau_v = kap_v*rho*drcor;
          Real dtau_ir = kap_ir*rho*dr;
            
          tau_v += dtau_v;
          Real fac = 1.0;
//          if (correct_spherical) fac = SQR(r1/r); // Zhang+2023
//          Real Q_v = kap_v*rho*Fstar*fac*exp(-tau_v/mu); // Zhang+2023
//          Q_v = (mu0 > 0.0) ? Q_v : 0.0;
          Real Qv = kap_v*rho*Fstar*fac*exp(-tau_v);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
//          u0(m,IEN,k,j,i) += Q_v*bdt;
        
          tau_ir_down[i] = tau_ir;
          tau_ir += dtau_ir/mu1;
          if (i==is) tau_ir_down[i-1] = tau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        Real rtop = x1v_(m,ie+1);
        I_down[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real r = x1v_(m,i);
          Real rp1 = x1v_(m,i+1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rp1/r);
          Real dtau = tau_ir_down[i]-tau_ir_down[i+1];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i)+B(i+1))/2.0;
          I_down[i] = (I_down[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        Real rbot = x1v_(m,is-1);
        Real I_up = Iint;
        // up-sweep
        for (int i=is; i<ie+1; ++i) {
          Real r = x1v_(m,i);
          Real rm1 = x1v_(m,i-1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rm1/r);
          Real dtau = tau_ir_down[i-1]-tau_ir_down[i];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i-1)+B(i))/2.0;
          I_up = (I_up*trans + Bavg*(1.0-trans))*fac;
          Real J = (I_up+I_down[i])/2.0;
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
            Real kap_ir = 1.0e-2; //2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
//          Real Q_ir = 4.0*M_PI*kap_ir*rho*(J-B(i));
//          u0(m,IEN,k,j,i) += Q_ir*bdt;
            
          Real kk = -4.0*M_PI*kap_ir*rho*boltz_sigma/M_PI*bdt;
          Real e0, e;
          if (eos.IsGeneral()) {
            // Same implicit balance e - kk T^4 = bb, but a general EOS has no e = c_v T
            // with constant c_v, so Newton-Raphson runs on the internal energy directly:
            // F(e) = e - kk T(e)^4 - bb, with dT/de = 1/(d c_v). c_v is per unit mass in
            // CODE units while T here is in Kelvin, hence the temp_cgs factor.
            e0 = w0(m,IEN,k,j,i);
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            e = e0;
            // `tc` is T in CODE temperature, carried alongside the kelvin `T` so that
            // neither EOS call has to solve for it. The two-argument SpecificHeatCv(d,e)
            // and Temperature(d,e) are the COLD-START forms -- the first is documented
            // "setup-time use only" because it inverts e(d,T) itself, and the second
            // brackets from scratch -- so using them here cost two full root finds per
            // Newton step, per cell, per stage. T and e are consistent at the top of
            // every iteration, so c_v can be evaluated at the temperature already in
            // hand, and the refresh below only needs the previous T as a warm start.
            Real tc = T/eos.temp_cgs;
            for (int n=0; n<100; ++n) {
              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e,tc));
              Real de = e - kk*SQR(SQR(T)) - bb;
              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
              if (fabs(de) <= 1.0e-10*e)
                break;
              tc = eos.Temperature(rho,e,tc);
              T = tc*eos.temp_cgs;
            }
          } else {
            Real cv = Rgas*rho*igm1;
            e0 = cv*T;
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            // Newton-Raphson
            for (int n=0; n<100; ++n) {
              e = cv*T;
              Real de = e - kk*SQR(SQR(T)) - bb;
              T -= de / (cv - 4.0*kk*T*T*T);
              if (fabs(de) <= 1.0e-10*e)
                break;
            }
          }
          u0(m,IEN,k,j,i) += e-e0;
        }
        
    });
    
    return;
}


KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td) {
    
    Real bar = 1.0e6;
    Real pl = 1.0e-3*bar;
    Real pt = log10(p/bar);
    Real ptl = log10(pl/bar);
    
    Real fn[13];
    Real fd[14];
    
    fn[0] = 1388.77348;
    fn[1] = 279.575848;
    fn[2] = -213.835822;
    fn[3] = 21.0010475;
    fn[4] = 100.938036;
    fn[5] = 12.7972336;
    fn[6] = -13.9266925;
    fn[7] = -3.70783272;
    fn[8] = 0.522370269;
    fn[9] = 0.320837882;
    fn[10]= 0.0451831612;
    fn[11]= 2.18195583e-3;
    fn[12]= 3.98938097e-6;

    fd[0] = 2152.06036;
    fd[1] = 29.3485512;
    fd[2] = -183.318696;
    fd[3] = 46.3893130;
    fd[4] = 19.8116485;
    fd[5] = -28.5473177;
    fd[6] = -2.52726545;
    fd[7] = 8.43627538;
    fd[8] = 2.62945375;
    fd[9] = -0.297098168;
    fd[10]= -0.286871487;
    fd[11]= -0.0590629443;
    fd[12]= -5.38679474e-3;
    fd[13]= -1.89972415e-4;
    
    Real Tnstar = 0.0;
    Real Tdstar = 0.0;
    Real Tnstar_pl = 0.0;
    Real Tdstar_pl = 0.0;
    Real ptn = 1.0;
    Real ptln = 1.0;
    for(int ilogp=0; ilogp<13; ilogp++) {
        Tnstar += fn[ilogp]*ptn;
        Tnstar_pl += fn[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }
    ptn = 1.0;
    ptln = 1.0;
    for(int ilogp=0; ilogp<14; ilogp++) {
        Tdstar += fd[ilogp]*ptn;
        Tdstar_pl += fd[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }

    Tn = Tnstar;
    Td = Tdstar;
//    if (p < pl) {
//        Tn = Tnstar_pl * exp(0.1*log10(p/pl));
//        if (Tn < 250) {
//            Tn = 250;
//        }
//        Td = Tdstar_pl * exp(0.015*log10(p/pl));
//        if (Td < 1000) {
//            Td = 1000;
//        }
//    }
//    if (Td < Tn) {
//        Td = Tn;
//    }
    Real x = log10(p / pl);

    Real Tn_new = Tnstar_pl;// * exp(0.1   * x);
    Real Td_new = Tdstar_pl;// * exp(0.015 * x);

    Tn_new = fmax(Tn_new, 250.0);
    Td_new = fmax(Td_new, 1000.0);

    // blend
    Tn = (p < pl) ? Tn_new : Tn;
    Td = (p < pl) ? Td_new : Td;

    // enforce Td >= Tn
    Td = fmax(Td, Tn);

    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
//    Real Tn, Td;
//    get_daynight_Tp(p, Tn, Td);
//
////    T = Tn;
//
//    Real Tn2 = Tn*Tn;
//    Real Tn4 = Tn2*Tn2;
//    Real Td2 = Td*Td;
//    Real Td4 = Td2*Td2;
//    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;
//
//    Real Tmid = sqrt(sqrt(Tmid4));
//    T = Tmid;
//
////    Real bar = 1.0e6;
////    Real Teq = 1469.0;
//////    Real Tirr = Teq*sqrt(2);
//////    Real Tint = 100.0;
//////    Real g = 942.0;
//////    Real mus = cos(50.0/90.0*M_PI/2.0);
//////    Real fH = 0.5;
//////    Real fK = 1.0/3.0;
//////    Real kap_v = 4.0e-3;
//////    Real kap_ir = 1.0e-2;
//////    Real gam = kap_v/kap_ir;
//////    Real tau = p/(g/kap_ir);
//////    tau = (tau < 0.0) ? 0.0 : tau;
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
////////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
//////    T = sqrt(sqrt(T4)); //2581.5574;
////    T = Teq;
    
    
    return;
}

template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
//    T = Tn;
    
    Real Tn2 = Tn*Tn;
    Real Tn4 = Tn2*Tn2;
    Real Td2 = Td*Td;
    Real Td4 = Td2*Td2;
    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;

    Real Tmid = sqrt(sqrt(Tmid4));
    T = Tmid;
    
//    Real bar = 1.0e6;
//    Real Teq = 1469.0;
////    Real Tirr = Teq*sqrt(2);
////    Real Tint = 100.0;
////    Real g = 942.0;
////    Real mus = cos(50.0/90.0*M_PI/2.0);
////    Real fH = 0.5;
////    Real fK = 1.0/3.0;
////    Real kap_v = 4.0e-3;
////    Real kap_ir = 1.0e-2;
////    Real gam = kap_v/kap_ir;
////    Real tau = p/(g/kap_ir);
////    tau = (tau < 0.0) ? 0.0 : tau;
////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
////    T = sqrt(sqrt(T4)); //2581.5574;
//    T = Teq;
    
    return;
}


template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr) {
    
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_wb_Tp(p,T);

        if (eos.IsGeneral()) {
          // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
          // branch is the same thing with rho = p/(Rgas T), kept in its original form.
          Real rho = DensFromPT(eos, Rgas, p, T);
          logparr(n+1) = logparr(n) + grav_acc*dz*rho/p;
        } else {
          logparr(n+1) = logparr(n) + fac/T;
        }
    }

    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_wb_Tp(p,T);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_wb_Tp(p0,T0);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr) {

//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_init_Tp_host(N, Tarr, lgparr, p, T);

        if (eos.IsGeneral()) {
          // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
          // branch is the same thing with rho = p/(Rgas T), kept in its original form.
          Real rho = DensFromPT(eos, Rgas, p, T);
          logparr(n+1) = logparr(n) + grav_acc*dz*rho/p;
        } else {
          logparr(n+1) = logparr(n) + fac/T;
        }
    }


    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {

//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tp(N, Tarr, lgparr, p,T);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_init_Tp(N, Tarr, lgparr, p0,T0);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
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
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2) {
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
  // Freedman+2014
  Real T1 = T;
  Real p1 = p;
  if (T > 4000.0) T1 = 4000.0;
  if (T < 75.0) T1 = 75.0;
  if (p > 3.0e8) p1 = 3.0e8;
  if (p < 1.0) p1 = 1.0;
  Real lgT = log10(T1);
  Real lgp = log10(p1);
  Real c1 = 10.602;
  Real c2 = 2.882;
  Real c3 = 6.09e-15;
  Real c4 = 2.954;
  Real c5 = -2.526;
  Real c6 = 0.843;
  Real c7 = -5.490;
  Real c8, c9, c10, c11, c12, c13;
  if (T1 < 800.0) {
    c8 = -14.051;
    c9 = 3.055;
    c10 = 0.024;
    c11 = 1.877;
    c12 = -0.445;
    c13 = 0.8321;
  } else {
    c8 = 82.241;
    c9 = -55.456;
    c10 = 8.754;
    c11 = 0.7048;
    c12 = -0.0414;
    c13 = 0.8321;
  }
    
  Real lgkl = c1*atan(lgT-c2) - c3/(lgp+c4)*exp(SQR(lgT-c5)) + c6*met + c7;
  Real lgkh = c8 + c9*lgT + c10*SQR(lgT) + lgp*(c11+c12*lgT) + c13*met*(0.5+1.0/M_PI*atan((lgT-2.5)/0.2));
  Real kl = pow(10.0,lgkl);
  Real kh = pow(10.0,lgkh);
  kapr = kl + kh;
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

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real Teq = Tirr/sqrt(2);
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);
    
    Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
    Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
    get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
    
    Real R = gamir1/gamir2;
    Real gamp = gamir1 + gamir2 - SQR(gamir2)*R;
    taulim = sqrt(R/3.0)*sqrt(beta*SQR(R-1.0)-SQR(beta*(R-1.0))+R)/SQR(gamir1);
    Real At1 = SQR(gamir1)*log(1.0+1.0/(taulim*gamir1));
    Real At2 = SQR(gamir2)*log(1.0+1.0/(taulim*gamir2));
    
    Real a0 = 1.0/gamir1 + 1.0/gamir2;
    Real a1 = -1.0/(3.0*SQR(taulim))*(gamp/(1.0-gamp)*(gamir1+gamir2-2.0)/(gamir1+gamir2) + (gamir1+gamir2)*taulim - (At1+At2)*SQR(taulim));
    Real b0 = 1.0/(gamir1*gamir2/(gamir1-gamir2)*(At1-At2)/3.0 - SQR(gamir1*gamir2)/sqrt(3.0*gamp) - SQR(gamir1*gamir2)*gamir1*gamir2/(1.0-gamir1)/(1.0-gamir2)/(gamir1+gamir2));
    A = 1.0/3.0*(a0+a1*b0);
    B = -1.0/3.0*SQR(gamir1*gamir2)/gamp*b0;
    
//    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    
    for (int iv=0; iv<3; ++iv) {
        Real gamv;
        if (iv==0) gamv = gamv1;
        if (iv==1) gamv = gamv2;
        if (iv==2) gamv = gamv3;
        gamvv[iv] = gamv;
        
        Real longf = (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(gamir1+gamir2) - 3.0*gamv*(6.0*SQR(gamir1*gamir2)-SQR(gamv)*(SQR(gamir1)+SQR(gamir2)));
        Real a2 = SQR(taulim)/(gamp*SQR(gamv)) * longf/(1.0-SQR(gamv*taulim));
        Real Av1 = SQR(gamir1)*log(1.0+gamv/gamir1);
        Real Av2 = SQR(gamir2)*log(1.0+gamv/gamir2);
        Real a3 = -SQR(taulim) * (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(Av1+Av2) / (gamp*SQR(gamv)*gamv*(1.0-SQR(gamv*taulim)));
        Real b1 = gamir1*gamir2*(3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*SQR(taulim) / (gamp*SQR(gamv)*(SQR(gamv*taulim)-1.0));
        Real b2 = 3.0*(gamir1+gamir2)*SQR(gamv)*gamv / ((3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv)));
        Real b3 = (Av2-Av1) / (gamv*(gamir1-gamir2));
        
        C[iv] = -1.0/3.0*(b0*b1*(1.0+b2+b3)*a1 + a2 + a3);
        D[iv] = -B*b1*(1.0+b2+b3);
        E[iv] = (3.0-SQR(gamv/gamir1))*(3.0-SQR(gamv/gamir2)) / (9.0*gamv*(SQR(gamv*taulim)-1.0));
//        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C + D*exp(-tau/taulim) + E*exp(-gamv*tau));
    }
    
//    T = sqrt(sqrt(T4));
    return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real taulim, const Real &A, const Real &B, const Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamv)[3], const Real &tau, Real &T) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    for (int iv=0; iv<3; ++iv) {
        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C[iv] + D[iv]*exp(-tau/taulim) + E[iv]*exp(-gamv[iv]*tau));
    }
    T = sqrt(sqrt(T4));
    return;
}

template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr) {
    Real bar = 1.0e6;
    Real tautop = 1.0e-6;
    Real Ttop;
    Real taulim, A, B, C[3], D[3], E[3], gamv[3];
    get_picket_fence_Ttau_coeff(Tint, Tirr, met, grav, mus, taulim, A, B, C, D, E, gamv);
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    Real kapr, ptop;
    ptop = 2.0; // 2e-6 bar
    get_kapr(Ttop, ptop, met, kapr);
    tautop = kapr/grav*ptop;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    
    Real lgptop = log10(ptop);
    Real pbot = 300.0*bar;
    Real lgpbot = log10(pbot);
    Real dlgp = (lgpbot-lgptop)/(N-1);
    Real lgp = lgptop;
    Real p = pow(10.0,lgp);
    Real tau = tautop;
    Real lgtau = log10(tau);
    Real T;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
    get_kapr(Ttop, p, met, kapr);
    for (int ip=0; ip<N; ++ ip) {
        Tarr(ip) = T;
        lgparr(ip) = lgp;
        Real K = p*kapr/grav*log(10.0);
        tau += K*dlgp;
        lgp += dlgp;
        p = pow(10.0,lgp);
//        lgp += dlgp;
//        Real pp = pow(10.0,lgp);
//        Real dp = pp-p;
//        Real K = kapr/grav;
//        tau += K*dp;
//        p = pp;
        get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
        get_kapr(T, p, met, kapr);
    }
    
    adjust_ad_pT_arr(eos, Rgas, gamma, N, Tarr, lgparr);
    
    return;
}

template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr) {

  // --- find convective boundary (search from bottom) ---
  int ic = -1;

  for (int ip = 0; ip < N-1; ++ip) {
    int i_inv = N - 1 - ip;
    int i_inv_p1 = N - 1 - (ip + 1);

    Real lgT  = log10(Tarr(i_inv));
    Real lgT1 = log10(Tarr(i_inv_p1));

    Real lgp  = lgparr(i_inv);
    Real lgp1 = lgparr(i_inv_p1);

    Real nabla = (lgT - lgT1) / (lgp - lgp1);

    Real T = Tarr(i_inv);
    Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgp), T);

    if (nabla < nabla_ad) {
      ic = i_inv_p1;  // map back to original indexing
      break;
    }
  }

  // fallback if no crossing found
  if (ic > 0) {
    // --- enforce adiabat downward ---
    for (int ip = ic; ip < N-1; ++ip) {
      Real T = Tarr(ip);
      Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgparr(ip)), T);

      Tarr(ip+1) = pow(10.0,
        nabla_ad * (lgparr(ip+1) - lgparr(ip)) + log10(T)
      );
    }
  }
    
  return;
}

void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
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
    get_Tint(Teq, Tint);
    Real Tint4 = SQR(SQR(Tint));
    Real Iint = boltz_sigma/M_PI*Tint4;
    
    const int nchain_rt = rt_nchain;
    const bool use_ktab_ = rt_ktab;
    auto ktab_ = *rt_ktab_ptr;
    auto rt_wgt_ = *rt_wgt_ptr;

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
        if (n1 > RT_NNC) {
          std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_split needs "
                    << "RT_NNC >= n1, but RT_NNC = " << RT_NNC << " and n1 = " << n1
                    << ". Rebuild with -DRT_NNC=" << n1 << " or larger." << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // Deliberately leaked, like the k-table: a namespace-scope View would outlive
        // Kokkos::finalize(). rt_Fb is the only large one, nmb*nblk*n3*n2*n1 Reals.
        rt_tau_ptr = new DvceArray4D<Real>("rt_tau", nmb, n3, n2, n1);
        rt_B_ptr   = new DvceArray4D<Real>("rt_B",   nmb, n3, n2, n1);
        rt_Qv_ptr  = new DvceArray4D<Real>("rt_Qv",  nmb, n3, n2, n1);
        rt_cf_ptr  = new DvceArray4D<Real>("rt_cf",  nmb, n3, n2, 4);
        rt_Fb_ptr  = new DvceArray5D<Real>("rt_Fb",  nmb, nblk, n3, n2, n1);
        if (rt_ck) {
          rt_kc_ptr = new DvceArray5D<Real>("rt_kc", nmb, CK_NB, n3, n2, n1);
          rt_Bb_ptr = new DvceArray5D<Real>("rt_Bb", nmb, CK_NB, n3, n2, n1);
          rt_T_ptr  = new DvceArray4D<Real>("rt_T",  nmb, n3, n2, n1);
          rt_pb_ptr = new DvceArray4D<Real>("rt_pb", nmb, n3, n2, n1);
          rt_xT_ptr = new DvceArray4D<Real>("rt_xT", nmb, n3, n2, n1);
          rt_xP_ptr = new DvceArray4D<Real>("rt_xP", nmb, n3, n2, n1);
          rt_icut_ptr = new DvceArray3D<int>("rt_icut", nmb, n3, n2);
          rt_Qb_ptr = new DvceArray5D<Real>("rt_Qb", nmb, nblk, n3, n2, n1);
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
      const bool ck_on = rt_ck;
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
        Real tau_r_f = kapr*p/grav;
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = tau_r_f/(kapr*rho);
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
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
          Real Qv = kapr*rho*(1.0-albedo)*Fstar*1.0/3.0*(gamv1*trans1+gamv2*trans2+gamv3*trans3);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
        }
        cf_g(m,k,j,0) = gamir1;
        cf_g(m,k,j,1) = gamir2;
        cf_g(m,k,j,2) = beta;
        cf_g(m,k,j,3) = mu0;

        if (ck_on) {
          // Everything the chain kernel needs that depends only on the cell: T, p, the
          // band Planck functions and the continuum. The continuum in particular returns
          // all CK_NB bands at once, so computing it here rather than per chain block
          // avoids repeating it nblk times per cell.
          // i = is is the BOTTOM of the column and i = ie+1 the top, so pressure falls
          // as i rises. The cut is the deepest cell still shallower than pcut, i.e. the
          // first one found scanning up from the bottom.
          int icut = ie+1;
          bool found = false;
          for (int i=is; i<ie+2; ++i) {
            Real rho = w0(m,IDN,k,j,i);
            Real pp, TT;
            PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),pp,TT);
            const Real pbar = pp*1.0e-6;
            T_g(m,k,j,i) = TT;
            pb_g(m,k,j,i) = pbar;
            if (!found && pbar < pcut) { icut = i; found = true; }
          }
          icut_g(m,k,j) = icut;
          // Second pass, only over the correlated-k region. The continuum is by far the
          // most expensive thing here -- four CIA searches, H- and eleven bands -- so
          // building it for the deep column too would be a third of it wasted. The table
          // index is hoisted here as well: it depends only on the cell, and recomputing
          // it inside the chain kernel would repeat both searches once per block.
          for (int i=icut; i<ie+2; ++i) {
            const Real rho = w0(m,IDN,k,j,i);
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
                         TT, pbar, rho, kcb);
            const Real sigT4_pi = boltz_sigma/M_PI*SQR(SQR(TT));
            for (int b=0; b<CK_NB; ++b) {
              kc_g(m,b,k,j,i) = kcb[b];
              Bb_g(m,b,k,j,i) = sigT4_pi*ck_planck_frac(ckpf, pfl0, pfid, TT, b);
            }
          }
        }
      });

      // ---- B (correlated-k): one thread per (column, chain block) ------------------
      // Chain c decodes as c = ((band*CK_NG) + g)*2 + quad, so a block of RT_NB = 4
      // consecutive chains is two g-points x two angular quadrature points of ONE band,
      // which is what lets the block share a band's continuum and Planck function.
      //
      // The longwave needs no cumulative optical depth: dtau is kappa*rho*dr, purely
      // local, so going from a grey tau scaled by gamma to a per-chain kappa is a lookup
      // and nothing structural. Only the sweep's lower limit changes, from is to icut.
      if (ck_on) {
        par_for("rt_chain_ck", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
          constexpr int NC = RT_NB;
          constexpr int NN = RT_NNC;
          auto F_ir_f = Kokkos::subview(Fb_g, m, blk, k, j, Kokkos::ALL);
          auto Q_v_f  = Kokkos::subview(Qb_g, m, blk, k, j, Kokkos::ALL);
          for (int i=is; i<ie+2; ++i) { F_ir_f[i] = 0.0; Q_v_f[i] = 0.0; }
          const int icut = icut_g(m,k,j);
          if (icut > ie) return;                  // whole column deeper than the cut
          // Shortwave. This is the one part of the scheme that genuinely restructures:
          // the longwave only ever needs a LAYER optical depth, which is local, but the
          // direct stellar beam needs the CUMULATIVE depth from the top, so each chain
          // carries its own downward recurrence and it cannot live in the per-column
          // rt_pre. It rides along in the longwave down-sweep because the two share the
          // same kappa lookup at every cell -- doing it in a second kernel would pay for
          // that lookup twice.
          const Real mu0 = cf_g(m,k,j,3);
          const Real facsw = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          const bool lit = (mu0 > 0.0);
          Real tausw[NC];

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
          Real I_down[NC][NN];

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
                             + kc_g(m,bandc[cc],k,j,ie+1);
              const Real dtau = kap*ptop*1.0e6/grav;
              const Real trans = exp(-dtau/muc[cc]);
              I_down[cc][ie+1] = (1.0-trans)*Bb_g(m,bandc[cc],k,j,ie+1);
              tausw[cc] = dtau;               // beam already crossed the column above
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
              const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, gc[cc]) + kc_g(m,b,k,j,i);
              const Real x = kap*drho/muc[cc];
              const Real e0 = -expm1(-x);
              const Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              I_down[cc][i] = (1.0-e0)*I_down[cc][i+1]
                            + alp*Bb_g(m,b,k,j,i+1) + bet*Bb_g(m,b,k,j,i);
              // direct beam: accumulate to this cell's lower face, then deposit
              tausw[cc] += kap*drho;
              if (lit) {
                Q_v_f[i] += kap*rho*(1.0-albedo)*Fstar*ckswf(b)*wgc[cc]
                          * exp(-tausw[cc]*facsw);
              }
            }
          }
          // Bottom of the CORRELATED-K DOMAIN, not of the column. At the cut the grey
          // optical depth is of order 1e4, so the layer is thermalised to e^-tau and the
          // upward intensity is its own Planck function; the planet's internal flux is
          // delivered here as an extra band-weighted source. Below the cut nothing
          // radiative is applied -- that region is optically thick and convective, and
          // the flux simply passes through it.
          Real I_up[NC];
          for (int cc=0; cc<NC; ++cc) {
            const int b = bandc[cc];
            const Real Iint_b = boltz_sigma/M_PI*Tint4
                              * ck_planck_frac(ckpf, pfl0, pfid, Tint, b);
            I_up[cc] = Bb_g(m,b,k,j,icut) + Iint_b;
            F_ir_f[icut] += wfc[cc]*(I_up[cc] - I_down[cc][icut]);
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
                             + kc_g(m,b,k,j,i-1);
              const Real x = kap*drho/muc[cc];
              const Real e0 = -expm1(-x);
              const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              const Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              I_up[cc] = (1.0-e0)*I_up[cc]
                       + bet*Bb_g(m,b,k,j,i) + gm*Bb_g(m,b,k,j,i-1);
              F_ir_f[i] += wfc[cc]*(I_up[cc] - I_down[cc][i]);
            }
          }
        });
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
          F_ir_f[i] = 0.0;
        }
          Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC], wtc[NC];
          bool synth[NC];
          for (int cc=0; cc<NC; ++cc) {
            const int c = blk*NC + cc;
            const int n = (c/2) % 2;
            const int vir = c % 2;
            muggc[cc] = mug[n];
            wggc[cc] = wg[n];
            gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
            fbc[cc] = (vir == 0) ? beta : (1.0-beta);
            wtc[cc] = rt_wgt_(c);
            synth[cc] = (c >= 4);
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
            Real ct = 0.0, cp = 0.0;
            if (use_ktab_) {
              ct = 1.2*log(1.0 + B[i]);
              cp = 0.7*log(1.0 + tau_down_r_f[i]);
            }
            for (int cc=0; cc<NC; ++cc) {
              Real gam = gamirc[cc];
              if (use_ktab_ && synth[cc]) {
                // Per-cell correlated-k style lookup: bilinear in a (T, p) proxy pair.
                // This measures the table traffic a real k-table would generate -- the
                // indices are data dependent and differ per cell, per chain -- it is not
                // a physical opacity, and it only ever feeds a zero-weight chain.
                gam = ktab_lookup(ktab_, blk*NC + cc, ct, cp);
              }
              Real dtauir = gam*dtau_i;
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
            F_ir_f[is] += wtc[cc]*(F_ir_up_f - F_ir_down_f);
          }
          // up-sweep, accumulating the band flux as it goes
          for (int i=is+1; i<ie+2; ++i) {
#if !RT_CACHE
            Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
            Real ct = 0.0, cp = 0.0;
            if (use_ktab_) {
              ct = 1.2*log(1.0 + B[i-1]);
              cp = 0.7*log(1.0 + tau_down_r_f[i-1]);
            }
#endif
            for (int cc=0; cc<NC; ++cc) {
#if RT_CACHE
              // layer i-1, already solved on the way down
              const Real e0 = e0c[cc][i-1];
              const Real bet = betc[cc][i-1];
              const Real gm = alpc[cc][i-1];
#else
              Real gam = gamirc[cc];
              if (use_ktab_ && synth[cc]) {
                gam = ktab_lookup(ktab_, blk*NC + cc, ct, cp);
              }
              Real dtauir = gam*dtau_i;
              Real x = dtauir/muggc[cc];
              Real e0 = -expm1(-x);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
#endif
              I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]
                            + bet*fbc[cc]*B[i] + gm*fbc[cc]*B[i-1];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
              F_ir_f[i] += wtc[cc]*(F_ir_up_f - F_ir_down_f);
            }
          }
      });
      }

      // ---- C: reduce over blocks in order, then apply ------------------------------
      par_for("rt_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real Ft = 0.0, Fb = 0.0;
        for (int b=0; b<nblk; ++b) {
          Ft += Fb_g(m,b,k,j,i+1);
          Fb += Fb_g(m,b,k,j,i);
        }
        Real src = -(Ft-Fb)/dx1(m,k,j,i);
        if (ck_on) {
          // deeper than the cut nothing radiative is applied: that region is optically
          // thick and convective, and the stellar beam died decades of optical depth above
          if (i < icut_g(m,k,j)) {
            src = 0.0;
          } else {
            Real Qs = 0.0;
            for (int b=0; b<nblk; ++b) Qs += Qb_g(m,b,k,j,i);
            src += Qs;
          }
        } else {
          src += Qv_g(m,k,j,i);
        }
        u0(m,IEN,k,j,i) += src*bdt;
      });
      return;
    }

//    size_t scr_size = 8 * ScrArray1D<Real>::shmem_size(n1);
//    int scr_level = 0;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, scr_level,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_down_r_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> B(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_ir_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> kapJ_ir(member.team_scratch(scr_level), n1);
        constexpr int NN = 270;
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
        Real tau_r_f = kapr*p/grav;
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = tau_r_f/(kapr*rho);
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
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
          Real Qv = kapr*rho*(1.0-albedo)*Fstar*1.0/3.0*(gamv1*trans1+gamv2*trans2+gamv3*trans3);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
        }
        
        // 2 IR Bands x two quadrature points, interleaved, in blocks of NC.
        //
        // Each (band, quadrature) combination is an independent pair of linear
        // recurrences in radius, and running them one after another leaves the wavefront
        // stalled on a single dependency chain: this kernel has only nmb*nx3*nx2/64
        // wavefronts for 912 SIMDs, so there is no other wave to hide that latency and
        // VALUBusy sits near 3 %. Stepping NC combinations inside one radial loop gives
        // the chain NC independent strands to overlap, and blocking keeps the private
        // I_ir_down_c footprint at NC columns however many chains are requested.
        //
        // Chains beyond the first four are the SCALING HARNESS (problem/rt_nchain, see
        // the note at the head of this function): they do the same lookup, the same two
        // recurrences and the same flux accumulation, but their weight rt_wgt is zero at
        // runtime, so the atmosphere is untouched and the run stays on the trajectory of
        // the production scheme. Multiplying by an exact 1.0 or 0.0 cannot perturb the
        // sum, so a harness run must stay bitwise identical -- that is the check that
        // the measurement is honest.
        constexpr int NC = RT_NB;
        const int nblk_rt = (nchain_rt + NC - 1)/NC;
        for (int i=is; i<ie+2; ++i) {
          F_ir_f[i] = 0.0;
        }

        for (int blk=0; blk<nblk_rt; ++blk) {
          Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC], wtc[NC];
          bool synth[NC];
          for (int cc=0; cc<NC; ++cc) {
            const int c = blk*NC + cc;
            const int n = (c/2) % 2;
            const int vir = c % 2;
            muggc[cc] = mug[n];
            wggc[cc] = wg[n];
            gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
            fbc[cc] = (vir == 0) ? beta : (1.0-beta);
            wtc[cc] = rt_wgt_(c);
            synth[cc] = (c >= 4);
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
            Real ct = 0.0, cp = 0.0;
            if (use_ktab_) {
              ct = 1.2*log(1.0 + B[i]);
              cp = 0.7*log(1.0 + tau_down_r_f[i]);
            }
            for (int cc=0; cc<NC; ++cc) {
              Real gam = gamirc[cc];
              if (use_ktab_ && synth[cc]) {
                // Per-cell correlated-k style lookup: bilinear in a (T, p) proxy pair.
                // This measures the table traffic a real k-table would generate -- the
                // indices are data dependent and differ per cell, per chain -- it is not
                // a physical opacity, and it only ever feeds a zero-weight chain.
                gam = ktab_lookup(ktab_, blk*NC + cc, ct, cp);
              }
              Real dtauir = gam*dtau_i;
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
            F_ir_f[is] += wtc[cc]*(F_ir_up_f - F_ir_down_f);
          }
          // up-sweep, accumulating the band flux as it goes
          for (int i=is+1; i<ie+2; ++i) {
            Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
            Real ct = 0.0, cp = 0.0;
            if (use_ktab_) {
              ct = 1.2*log(1.0 + B[i-1]);
              cp = 0.7*log(1.0 + tau_down_r_f[i-1]);
            }
            for (int cc=0; cc<NC; ++cc) {
              Real gam = gamirc[cc];
              if (use_ktab_ && synth[cc]) {
                gam = ktab_lookup(ktab_, blk*NC + cc, ct, cp);
              }
              Real dtauir = gam*dtau_i;
              Real x = dtauir/muggc[cc];
              Real e0 = -expm1(-x);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]
                            + bet*fbc[cc]*B[i] + gm*fbc[cc]*B[i-1];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i];
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc];
              F_ir_f[i] += wtc[cc]*(F_ir_up_f - F_ir_down_f);
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
          u0(m,IEN,k,j,i) += du;
        }
//        });
        
//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();
        
    });
    
    return;
}

