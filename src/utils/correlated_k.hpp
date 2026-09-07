#ifndef UTILS_CORRELATED_K_HPP_
#define UTILS_CORRELATED_K_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file correlated_k.hpp
//! \brief The Exo-FMS correlated-k OPACITY layer: the k-table, the continuum (CIA,
//! Rayleigh and H-), the equilibrium-chemistry table, the band-integrated Planck
//! fractions, the device-side lookups, the start-up self-tests, and the Rosseland mean
//! that couples all of it to the radiative-conduction module.
//!
//! EXTRACTED from src/pgen/deep_hot_jupiter_rt.cpp, where it grew, so that a second
//! problem generator can use the same bands and the same tables.  The arithmetic is
//! untouched; the extraction is gated on a bitwise-identical problem/ck_dump_file
//! column.  What changed: the enclosing namespace, `inline` on the namespace-scope
//! definitions, and three signatures that now take as ARGUMENTS the two solver
//! parameters the reports used to read from problem-generator globals
//! (the correlated-k pressure cut and the host effective temperature), so that this
//! module owns no solver state.
//!
//! WHAT IS HERE AND WHAT IS NOT.  The opacity layer only: given (T, p, rho) it returns
//! kappa per band and g-point, and the band Planck fractions.  The two-stream SOLVER
//! that consumes them stays in the hot Jupiter generator, because it is written around
//! that problem's irradiation and arrays.
//!
//! MERGING WITH RADIATIVE CONDUCTION.  ck_build_rosseland_table() tabulates the
//! Rosseland mean of this same k-table plus continuum on the table's own (T, p) grid and
//! installs it in a Conduction object (<block>/rad_kappa_src = table).  That is what
//! makes the optical-depth blend consistent: the two-stream above and the diffusion
//! below then use ONE opacity, so the handover conserves flux instead of jumping at the
//! seam.  Call it once at start-up, after the tables are read.  A Planck mean, which an
//! optically thin relaxation wants and the Rosseland mean is a poor stand-in for, is the
//! obvious next addition: the weighting is the only thing that differs.
//!
//! USE.  On the host at start-up:
//!     correlated_k::read_ck_table(file, pcut_bar);
//!     correlated_k::build_planck_fractions(pcut_bar);
//!     correlated_k::read_ck_continuum(dir, swflux, star_teff);
//!     correlated_k::ck_selftest();   correlated_k::ck_rt_selftest();
//!     correlated_k::ck_build_rosseland_table(pcond);      // if blending with diffusion
//! then, in a kernel, copy the Views out of the pointers FIRST -- a device lambda must
//! not dereference a host pointer -- and call ck_tp_index, ck_kappa, ck_planck_bands and
//! ck_continuum.
//!
//! The tables are third-party and are not in git: see data/exo_fms_ck/PROVENANCE.md.

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "diffusion/conduction.hpp"
#include "utils/rosseland.hpp"

namespace correlated_k {


// problem/ck_nquad: angular treatment of the longwave. 1 uses a single diffusivity
// factor, mu = 1/1.66, which is what GCMs normally do and what the chain count in the
// correlated-k literature assumes. 2 keeps the 2-point Gauss quadrature the picket fence
// uses, which doubles the number of column solves. Both normalise to F = pi <I>: for the
// Gauss pair, 2 pi sum w mu = 2 pi (0.5) = pi.
inline int ck_nq = 1;
constexpr Real CK_DIFFUSIVITY = 1.66;
constexpr int CK_NB = 11;                  // bands, fixed by the Kataria grid
constexpr int CK_NG = 8;                   // g-points per band, fixed by the table
inline int ck_nT = 0;
inline int ck_nP = 0;
// Deliberately leaked, as with the other tables here: a namespace-scope Kokkos View
// would be destroyed after Kokkos::finalize().
inline DvceArray1D<Real> *ck_lT_ptr = nullptr;    // log10 T grid [K]
inline DvceArray1D<Real> *ck_lP_ptr = nullptr;    // log10 p grid [bar]
// Layout is (band, g, iT, iP): for one chain the entire (T,p) plane is contiguous and
// only
// 38*34*8 = 10 kB, so the four bilinear corners sit 8 and 272 bytes apart rather than the
// 704 bytes and 24 kB that (iT,iP,band,g) would give.
//
// This was expected to be worth about a factor of two in the chain kernel. It was
// MEASURED
// AND IT IS NOT: 2725 -> 2711 ms per 100 cycles, half a per cent. The k-table access
// pattern is simply not the bottleneck, which is the same answer an earlier experiment
// with a synthetic table gave when it blocked the lookup by band and gained nothing.
// The layout is kept because it is the more natural one, not because it is faster.
inline DvceArray4D<Real> *ck_lk_ptr = nullptr;    // (b,g,iT,iP) log10 kappa [cm^2/g]
inline DvceArray1D<Real> *ck_gw_ptr = nullptr;    // g-point weights, sum to 1
// band edges [um], descending, CK_NB+1 of them
inline DvceArray1D<Real> *ck_wl_ptr = nullptr;
// stellar flux FRACTION per band, sums to 1
inline DvceArray1D<Real> *ck_swf_ptr = nullptr;

// Band-integrated Planck function. The grey scheme uses B = sigma T^4 / pi; correlated-k
// needs B_b(T) = (sigma T^4 / pi) * f_b(T), the fraction of the Planck function falling
// in
// band b. f_b is smooth in log T, so it is tabulated once on a uniform log10 T grid and
// interpolated -- evaluating the series per cell per band would be absurd.
constexpr int CK_NPF = 512;
constexpr Real CK_PF_TMIN = 50.0;
constexpr Real CK_PF_TMAX = 20000.0;
inline DvceArray2D<Real> *ck_pf_ptr = nullptr;    // (iT, band) fractional Planck function
inline Real ck_pf_lTmin = 0.0;
inline Real ck_pf_idlT = 0.0;                     // 1 / grid spacing in log10 T

//----------------------------------------------------------------------------------------
//! \fn Real planck_fraction_below()
//  \brief fraction of blackbody emission at wavelengths SHORTER than lam, as a function
// of
//  lam*T in um K. The standard series (Chang & Rhee 1984): with xi = c2/(lam T),
//    F = 15/pi^4 sum_n e^{-n xi}/n (xi^3 + 3 xi^2/n + 6 xi/n^2 + 6/n^3).
//  Host side only, evaluated once per table entry at startup.

inline Real planck_fraction_below(const Real lamT) {
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
//  The Kataria grid spans 0.26 to 324.68 um and does NOT capture the whole Planck
// function:
//  at high T a real fraction escapes past the blue edge. That flux has to go somewhere or
//  the scheme silently loses energy, so the two outermost bands are extended to 0 and
//  infinity -- the sub-0.26 um tail joins the bluest band and the super-324.68 um tail
// the
//  reddest. It is an approximation, since those tails get their host band's kappa, but it
//  is a small one wherever the tails are small, and it makes sum_b f_b = 1 exactly.

inline void build_planck_fractions(const Real pcut_bar) {
  const Real rt_ck_pcut = pcut_bar;   // only for the advisory message below
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
                                           // the spectrum and the fold-in stops being
                                           // small
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
              << "  capturing the spectrum and folding the tail into the bluest "
              << "band is no"
              << " longer a small correction" << std::endl;
    if (blue1pc_T > 0.0 && blue1pc_T < CK_PF_TMAX) {
      std::cout << "  (with the p < " << rt_ck_pcut
                << " bar cut this setup tops out near "
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

// Templated on the view type so the SAME lookup serves the device kernels and the
// host-side photosphere diagnostic, whose mirrors live in a different memory space.
template <typename PFView>
KOKKOS_INLINE_FUNCTION
Real ck_planck_frac(const PFView &pf, const Real lTmin, const Real idlT,
                    const Real &T, const int b) {
  Real x = (log10(T) - lTmin)*idlT;
  // NaN DEFEATS AN ORDINARY CLAMP.  `NaN < 0` and `NaN > hi` are BOTH false, so a
  // two-sided ternary passes NaN straight through, and static_cast<int>(NaN) is INT_MIN:
  // the lookup then indexes the table at -2147483648 and SEGFAULTS.  Measured -- a
  // cubed-sphere run whose upper atmosphere underflowed below the EOS table's 71 K floor
  // died exactly here, while spherical polar merely produced NaN and kept going.  Writing
  // the low test as !(x > 0) makes NaN take the low branch, which turns a crash into a
  // clamped value the surrounding NaN checks can still notice.
  x = !(x > 0.0) ? 0.0 : ((x > CK_NPF-1.0) ? CK_NPF-1.0 : x);
  int i = static_cast<int>(x);
  i = (i > CK_NPF-2) ? CK_NPF-2 : i;
  const Real f = x - static_cast<Real>(i);
  return (1.0-f)*pf(i,b) + f*pf(i+1,b);
}

template <typename PFView>
KOKKOS_INLINE_FUNCTION
void ck_planck_bands(const PFView &pf, const Real lTmin, const Real idlT,
                     const Real &sigT4_pi, const Real &T, Real (&Bb)[CK_NB]) {
  Real x = (log10(T) - lTmin)*idlT;
  // NaN-safe, for the reason spelled out in ck_planck_frac above.
  x = !(x > 0.0) ? 0.0 : ((x > CK_NPF-1.0) ? CK_NPF-1.0 : x);
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

inline void read_ck_table(const std::string &fname, const Real pcut_bar) {
  const Real rt_ck_pcut = pcut_bar;   // only for the advisory message below
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
inline int ce_nT = 0;
inline int ce_nP = 0;
inline DvceArray1D<Real> *ce_lT_ptr = nullptr;
inline DvceArray1D<Real> *ce_lP_ptr = nullptr;
inline DvceArray3D<Real> *ce_ptr = nullptr;       // (iT,iP,CK_NCE)
inline DvceArray1D<int>  *cia_nT_ptr = nullptr;   // per-pair grid length
// (pair,iT) -- the four grids differ wildly,
inline DvceArray2D<Real> *cia_T_ptr = nullptr;
// (pair,iT,band)   200-3000 K to 200-9900 K
inline DvceArray3D<Real> *cia_k_ptr = nullptr;
// (species,band) cm^2/molecule, T independent
inline DvceArray2D<Real> *ray_x_ptr = nullptr;

//----------------------------------------------------------------------------------------
//! \fn void read_ck_continuum()
//  \brief read the FastChem composition table, the four CIA pair tables and the Rayleigh
//  cross sections. All are whitespace-separated numbers after a short header.

inline void read_ck_continuum(const std::string &dir, const std::string &swfile,
                              const Real star_teff) {
  const Real rt_star_teff = star_teff;   // the stellar band fractions are built here
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
    if (rt_star_teff <= 0.0 && !f.is_open()) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: could not open '" << fn
                << "'. Set problem/ck_swflux." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ck_swf_ptr = new DvceArray1D<Real>("ck_swf", CK_NB);
    auto hs = Kokkos::create_mirror_view(*ck_swf_ptr);
    auto hwl = Kokkos::create_mirror_view(*ck_wl_ptr);
    Kokkos::deep_copy(hwl, *ck_wl_ptr);
    Real blue_tail = 0.0;
    if (rt_star_teff > 0.0) {
      // Blackbody host at problem/ck_star_teff. Only the SHAPE of the stellar spectrum is
      // used -- it is renormalised to the code's own sigma T_irr^4 -- so this needs the
      // host's effective temperature and nothing else, and it needs no external data.
      // Same tail handling as the Planck fractions: flux outside the grid is folded into
      // the outermost bands so that all of the insolation is deposited somewhere.
      Real tot = 0.0;
      for (int b=0; b<CK_NB; ++b) {
        hs(b) = planck_fraction_below(hwl(b)*rt_star_teff)
              - planck_fraction_below(hwl(b+1)*rt_star_teff);
        tot += hs(b);
      }
      blue_tail = planck_fraction_below(hwl(CK_NB)*rt_star_teff);
      hs(0) += 1.0 - planck_fraction_below(hwl(0)*rt_star_teff);
      hs(CK_NB-1) += blue_tail;
    } else {
      Real v[CK_NB];
      Real tot = 0.0;
      for (int b=0; b<CK_NB; ++b) { f >> v[b]; tot += v[b]; }
      if (!f || tot <= 0.0) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: bad stellar flux file '"
                  << fn << "'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      for (int b=0; b<CK_NB; ++b) hs(b) = v[CK_NB-1-b]/tot;  // reverse into wl order
    }
    Kokkos::deep_copy(*ck_swf_ptr, hs);
    if (global_variable::my_rank == 0) {
      if (rt_star_teff > 0.0) {
        std::cout << "  stellar spectrum: blackbody at T_eff = " << rt_star_teff << " K"
                  << std::endl;
        std::cout << "    " << 100.0*blue_tail << " % of it is bluer than the grid's "
                  << hwl(CK_NB) << " um edge and is folded into the bluest band";
        if (blue_tail > 0.05) {
          std::cout << " -- THAT IS A LOT. The band grid does not cover this host; its "
                    << "UV flux will be deposited with near-UV opacities, hence "
                    << "too deep.";
        }
        std::cout << std::endl;
      } else {
        std::cout << "  stellar spectrum '" << swfile << "' (file), band fractions "
                  << hs(CK_NB-1) << " (bluest) .. " << hs(0) << " (reddest)" << std::endl;
      }
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
//  sign-flipped optical depth increment, and an intensity recurrence that diverged to
// inf.
//  A table lookup that silently extrapolates 40 decades of kappa is not a lookup.

KOKKOS_INLINE_FUNCTION
void ck_tp_index(const DvceArray1D<Real> &lg, const int n, const Real &x,
                 int &i, Real &f) {
  // !(x > lg(0)) rather than x <= lg(0), so a NaN takes this branch instead of falling
  // through to the bisection, which would return a valid index with a NaN FRACTION and
  // quietly poison kappa. See the note in ck_planck_frac.
  if (!(x > lg(0))) { i = 0; f = 0.0; return; }
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

template <typename LKView>
KOKKOS_INLINE_FUNCTION
Real ck_kappa(const LKView &lk, const int iT, const Real &fT,
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
//! \fn void ck_rt_selftest()
//  \brief two exact limits of the two-stream solver, run on the device at startup against
//  a synthetic column, using the same recurrence the chain kernel uses.
//
//  1. ISOTHERMAL INVARIANCE. An isothermal atmosphere bathed in its own Planck function
//     must stay in equilibrium: the net flux has to vanish identically at every level.
//     This works because alp + bet = (e0 - 1 + e0/x) + (1 - e0/x) = e0 exactly, so
//     I = (1-e0) I + e0 B has B as a fixed point. It is the sharpest check there is on
// the
//     recurrence coefficients and on the small-x branch, and it catches sign and
//     normalisation errors in them that a plausible-looking profile would hide.
//
//  2. TRANSPARENT SLAB. With the layers made optically thin and a blackbody floor, the
//     emergent flux must be exactly sigma T^4. That is a different statement: it tests
// the
//     WEIGHTS rather than the recurrence -- the g-point weights summing to one, the band
//     Planck fractions summing to one, and the flux prefactor being pi and not 2 pi.
//     Test 1 cannot see any of those, since it holds whatever the weights are.

inline void ck_rt_selftest() {
  auto lk = *ck_lk_ptr;
  auto pf = *ck_pf_ptr;
  auto gwv = *ck_gw_ptr;
  const Real l0 = ck_pf_lTmin;
  const Real id = ck_pf_idlT;
  const int nq = ck_nq;
  Real mugl[2];
  Real wgl[2];
  mugl[0] = 0.21132487;  mugl[1] = 0.78867513;
  wgl[0] = 0.5;          wgl[1] = 0.5;
  const Real boltz = 5.6704e-5;
  DvceArray1D<Real> out("ck_rt_selftest", 2);
  par_for("ck_rt_selftest", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int dummy) {
    const int NL = 40;                       // synthetic layers
    const Real Ttest = 2000.0;
    const Real sigT4_pi = boltz/M_PI*SQR(SQR(Ttest));
    Real worst_iso = 0.0;
    Real Fthin = 0.0;
    for (int tst=0; tst<2; ++tst) {
      const Real dtau = (tst == 0) ? 0.3 : 1.0e-12;   // thick-ish, then transparent
      Real Fnet[NL+2];
      for (int i=0; i<NL+2; ++i) Fnet[i] = 0.0;
      const int nch = CK_NB*CK_NG*nq;
      for (int c=0; c<nch; ++c) {
        int b, g;
        Real mu, wf;
        if (nq == 1) {
          g = c % CK_NG;  b = c/CK_NG;
          mu = 1.0/CK_DIFFUSIVITY;  wf = M_PI*gwv(g);
        } else {
          const int q = c % 2;
          g = (c/2) % CK_NG;  b = c/(2*CK_NG);
          mu = mugl[q];  wf = 2.0*M_PI*wgl[q]*mugl[q]*gwv(g);
        }
        const Real Bb = sigT4_pi*ck_planck_frac(pf, l0, id, Ttest, b);
        const Real x = dtau/mu;
        const Real e0 = -expm1(-x);
        const Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
        const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
        // test 1 bathes the top in B, test 2 leaves it dark
        Real Id[NL+2];
        Id[NL+1] = (tst == 0) ? Bb : 0.0;
        for (int i=NL; i>=0; --i) {
          Id[i] = (1.0-e0)*Id[i+1] + alp*Bb + bet*Bb;
        }
        Real Iu = Bb;                        // blackbody floor, no extra internal flux
        Fnet[0] += wf*(Iu - Id[0]);
        for (int i=1; i<NL+2; ++i) {
          Iu = (1.0-e0)*Iu + bet*Bb + alp*Bb;
          Fnet[i] += wf*(Iu - Id[i]);
        }
      }
      if (tst == 0) {
        for (int i=0; i<NL+2; ++i) {
          const Real a = fabs(Fnet[i])/(boltz*SQR(SQR(Ttest)));
          worst_iso = (a > worst_iso) ? a : worst_iso;
        }
      } else {
        Fthin = Fnet[NL+1]/(boltz*SQR(SQR(Ttest)));
      }
    }
    out(0) = worst_iso;
    out(1) = Fthin;
  });
  auto h = Kokkos::create_mirror_view(out);
  Kokkos::deep_copy(h, out);
  if (global_variable::my_rank == 0) {
    std::cout << "  two-stream self-test: isothermal net flux |F|/sigmaT^4 <= " << h(0)
              << ", transparent slab F/sigmaT^4 = " << h(1) << std::endl;
    if (!(h(0) < 1.0e-12)) {
      std::cout << "### FATAL ERROR in correlated_k: the two-stream recurrence does "
                << "not hold an isothermal atmosphere in equilibrium (worst "
                << "|F|/sigmaT^4 = "
                << h(0) << "). Check alp, bet and the small-x branch." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(fabs(h(1) - 1.0) < 1.0e-10)) {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: a transparent slab over a "
                << "blackbody floor emits " << h(1) << " sigma T^4, not 1. The g-point "
                << "weights, the band Planck fractions or the flux prefactor are wrong."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ck_continuum()
//  \brief grey-within-band continuum opacity [cm^2/g] for every band at one cell: the
// four
//  CIA pairs plus Rayleigh scattering off H2, He, H and e-.
//
//  CIA scales as the PRODUCT of the two collider number densities, Rayleigh as one, so
//  both need the equilibrium composition, taken from the FastChem table. The simulation's
//  own rho is used for the mass conversion rather than the table's, since that is the
//  density the rest of the scheme works with.
//
//  Every table index is clamped, and the four CIA pairs have grids that stop in very
//  different places -- 200-3000 K for H2-H2 but 200-9900 K for H2-He -- so clamping is
// not
//  an edge case here, it is the normal state of affairs above 3000 K. Held-flat
//  extrapolation of CIA is what Exo-FMS does too.
//
//  H- bound-free and free-free (John 1988) ARE included below, from the FastChem n(H-),
//  n(e-) and n(H); they dominate the continuum above about 3000 K. With them and CIA in,
//  the table's Rosseland mean is 1.0-1.8x Freedman+2014 at 3500 K
// (docs/correlated_k_rt.md).
//  NOT included: the H2- and He- free-free files that ship in cia/ (unread).

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
  // 3000 K, which under the p < 10 bar cut is the BOTTOM of the correlated-k region, so
  // it
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
    const Real Dff2[6] = {0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300,
                          2095.2880};
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
//  -- falls. If the band index is reversed these swap, which looks entirely plausible in
// a
//  plot and is completely wrong. Exo-FMS's own reader loops the other way, so this is not
//  a hypothetical.

inline void ck_selftest() {
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
      std::cout << "### FATAL ERROR in correlated_k: the correlated-k band ordering "
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
//! \fn void ck_build_rosseland_table()
//  \brief the Rosseland mean of the correlated-k table + continuum, tabulated ONCE on the
//  table's own (T, p) grid and handed to the conduction module, so the radiative
//  diffusion below the two-stream region uses the SAME opacity as the two-stream above
//  it (<mhd>/rad_kappa_src = table). 1/kappa_R = sum_b w_b sum_g gw_g/(k_bg + kc_b) /
//  sum_b w_b with w_b = d(sigma T^4 f_b)/dT from the band Planck-fraction table (centred
//  difference at +-1 % in T); the continuum (CIA, Rayleigh, H-) counts as extinction and
//  is evaluated at the equilibrium density p mu m_H/(k T) with mu from the FastChem
// table.
//  Outside the grid the conduction module holds the edge value.

inline void ck_build_rosseland_table(Conduction *pc) {
  if (pc == nullptr || !pc->rad_kappa_tab) return;
  if (ck_lk_ptr == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "rad_kappa_src = table needs the correlated-k tables: call "
              << "read_ck_table/read_ck_continuum first" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nT = ck_nT, nP = ck_nP;
  DvceArray2D<Real> tab("ck_kR", nT, nP);
  DvceArray1D<Real> lPcgs("ck_kR_lP", nP);
  auto cklk = *ck_lk_ptr;
  auto cklT = *ck_lT_ptr;
  auto cklP = *ck_lP_ptr;
  auto ckgw = *ck_gw_ptr;
  auto ckwl = *ck_wl_ptr;
  auto ckpf = *ck_pf_ptr;
  auto cece = *ce_ptr;
  auto celT = *ce_lT_ptr;
  auto celP = *ce_lP_ptr;
  auto cian = *cia_nT_ptr;
  auto ciaT = *cia_T_ptr;
  auto ciak = *cia_k_ptr;
  auto rayx = *ray_x_ptr;
  const int ceNT = ce_nT, ceNP = ce_nP;
  const Real pfl0 = ck_pf_lTmin, pfid = ck_pf_idlT;
  const Real boltz_sigma = 5.670374419e-5;
  par_for("ck_rosseland_tab", DevExeSpace(), 0, nT-1, 0, nP-1,
  KOKKOS_LAMBDA(const int it, const int ip) {
    const Real TT = pow(10.0, cklT(it));
    const Real pbar = pow(10.0, cklP(ip));
    if (it == 0) lPcgs(ip) = cklP(ip) + 6.0;
    // equilibrium mean molecular weight at this grid point -> density for the continuum
    int jT, jP;
    Real gT, gP;
    ck_tp_index(celT, ceNT, log10(TT), jT, gT);
    ck_tp_index(celP, ceNP, log10(pbar), jP, gP);
    const Real mu = (1.0-gT)*((1.0-gP)*cece(jT,jP,0) + gP*cece(jT,jP+1,0))
                  +      gT *((1.0-gP)*cece(jT+1,jP,0) + gP*cece(jT+1,jP+1,0));
    const Real rho = pbar*1.0e6*mu*1.6726e-24/(1.380649e-16*TT);
    Real kcb[CK_NB];
    ck_continuum(cece, celT, celP, ceNT, ceNP, cian, ciaT, ciak, rayx, ckwl,
                 TT, pbar, rho, kcb);
    const Real Tp = 1.01*TT, Tm = 0.99*TT;
    const Real sp = boltz_sigma*SQR(SQR(Tp)), sm = boltz_sigma*SQR(SQR(Tm));
    Real num = 0.0, den = 0.0;
    for (int b=0; b<CK_NB; ++b) {
      const Real wb = sp*ck_planck_frac(ckpf, pfl0, pfid, Tp, b)
                    - sm*ck_planck_frac(ckpf, pfl0, pfid, Tm, b);
      if (!(wb > 0.0)) continue;
      Real inv = 0.0;
      for (int g=0; g<CK_NG; ++g) {
        inv += ckgw(g)/(ck_kappa(cklk, it, 0.0, ip, 0.0, b, g) + kcb[b]);
      }
      num += wb*inv;
      den += wb;
    }
    tab(it,ip) = log10(den/num);
  });
  pc->rad_kr_tab = tab;
  pc->rad_kr_lT = cklT;
  pc->rad_kr_lP = lPcgs;
  pc->rad_kr_nT = nT;
  pc->rad_kr_nP = nP;
  // a few values against the Freedman fit, so a wrong band order or unit shows here
  auto htab = Kokkos::create_mirror_view(tab);
  Kokkos::deep_copy(htab, tab);
  auto hlT = Kokkos::create_mirror_view(cklT);
  Kokkos::deep_copy(hlT, cklT);
  auto hlP = Kokkos::create_mirror_view(cklP);
  Kokkos::deep_copy(hlP, cklP);
  if (global_variable::my_rank == 0) {
    std::cout << "  Rosseland table for the radiative diffusion (ck + continuum), "
              << nT << " T x " << nP << " p; kappa_R / Freedman+2014:" << std::endl;
    const int its[4] = {9, 19, 25, 31}, ips[2] = {24, 28};   // 1000/2500/3700/4900 K
    for (int a=0; a<4; ++a) {
      for (int c=0; c<2; ++c) {
        const Real T = std::pow(10.0, hlT(its[a])), pb = std::pow(10.0, hlP(ips[c]));
        const Real kr = std::pow(10.0, htab(its[a],ips[c]));
        std::cout << "    T = " << T << " K, p = " << pb << " bar: " << kr << " / "
                  << RosselandFreedman2014(T, pb*1.0e6, 0.0) << " = "
                  << kr/RosselandFreedman2014(T, pb*1.0e6, 0.0) << std::endl;
      }
    }
  }
  return;
}

}  // namespace correlated_k

#endif  // UTILS_CORRELATED_K_HPP_
