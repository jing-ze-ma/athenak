#ifndef EOS_EOS_COMPOSITION_HPP_
#define EOS_EOS_COMPOSITION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_composition.hpp
//! \brief analytic equation of state for a partially dissociated, partially ionized
//! hydrogen/helium mixture, evaluated in (rho,T) and in cgs units.
//!
//! This is the PHYSICS model behind the general EOS. It is HOST-ONLY and deliberately
//! makes no attempt to be fast: it is evaluated once per node when the interpolation
//! table is built (see eos_table.{hpp,cpp}) and never on the hot path. Everything the
//! solver actually calls comes from the table.
//!
//! Species: H2, H, H+, He, He+, He++, free electrons, and a "metal" component that
//! contributes particles and, optionally, singly ionizing electron donors. The
//! composition is
//! obtained from
//!   - the Saha equation for H -> H+ + e, He -> He+ + e and He+ -> He++ + e, and
//!   - the dissociation equilibrium H2 -> 2H,
//! closed by charge neutrality, which is solved as a single 1D root find on the electron
//! number density.
//!
//! The internal energy per unit mass is the sum of
//!   - translational energy (3/2)nkT of every free particle, electrons included,
//!   - the rotational and vibrational energy of H2, from its partition function, and
//!   - the chemical energy stored in dissociated and ionized material, measured from a
//!     zero point of "all hydrogen in ground-state H2, all helium neutral".
//!
//! NOT included, and the main limitations of this model: Coulomb (non-ideal) corrections,
//! electron degeneracy and pressure ionization, excited bound states beyond the ground
//! state statistical weights, and any ionization of the metal component. It is therefore
//! valid for the weakly coupled, non-degenerate regime -- gas-giant and stellar envelope
//! conditions -- and not for the deep interiors of stars or degenerate objects.
//!
//! Radiation is deliberately NOT included here. It is added analytically on top of the
//! tabulated gas contribution (see eos_table.hpp), which makes the runtime radiation
//! switch exact and keeps the tabulated surfaces smooth.

#include <math.h>

#include <algorithm>

//----------------------------------------------------------------------------------------
//! \namespace eos_cgs
//! \brief physical constants in cgs, and the atomic data used by the composition model

namespace eos_cgs {

constexpr double kboltz = 1.380649e-16;         // Boltzmann constant, erg/K
constexpr double hplanck = 6.62607015e-27;      // Planck constant, erg s
constexpr double m_u = 1.66053906660e-24;       // atomic mass unit, g
constexpr double m_el = 9.1093837015e-28;       // electron mass, g
constexpr double a_rad = 7.5657332503e-15;      // radiation constant, erg/cm^3/K^4
constexpr double ev = 1.602176634e-12;          // electron volt, erg

constexpr double a_hyd = 1.008;                 // atomic weight of H
constexpr double a_hel = 4.002602;              // atomic weight of He

constexpr double chi_h = 13.598434*ev;          // H ionization potential
constexpr double chi_he1 = 24.587389*ev;        // He   -> He+  ionization potential
constexpr double chi_he2 = 54.417765*ev;        // He+  -> He++ ionization potential
constexpr double chi_d = 4.478007*ev;           // H2 dissociation energy from (v=0,J=0)

constexpr double theta_rot = 85.4;              // H2 rotational temperature, K
constexpr double theta_vib = 6332.0;            // H2 vibrational temperature, K

//----------------------------------------------------------------------------------------
//! \struct MetalDonor
//! \brief one singly-ionizing metal, as an electron source
//!
//! In a cool, weakly ionized atmosphere the free electrons do NOT come from hydrogen:
//! 13.6 eV is far too high to matter below ~4000 K, while these species sit at 4-8 eV and
//! are substantially ionized by 2500 K. They set the electrical conductivity, and hence
//! the Ohmic resistivity, over the whole range where a hot-Jupiter atmosphere lives.
//!
//! `abun` is the solar abundance by NUMBER relative to hydrogen nuclei; `gfac` is
//! 2 g_ion/g_neutral from ground-state statistical weights. Only the first ionization is
//! tracked -- the second potentials are all above 11 eV and irrelevant here.
//!
//! Ground-state weights are a good approximation for the alkalis, whose first excited
//! states lie ~2 eV up, and a poorer one for Fe, whose many low-lying levels make the
//! true neutral partition function ~1.5-2x the ground-state 25 near 3000-5000 K. That
//! biases n_e by tens of percent where iron dominates, small against the orders of
//! magnitude at stake, but it is not exact.
//! `tc_a` and `tc_b` are the condensation curve, in the Clausius-Clapeyron form
//!
//!     T_cond(p, [M/H]) = tc_b / (tc_a - log10(p/1 bar) - 0.5 [M/H])
//!
//! which is the shape these curves actually have: T_cond rises with pressure and, more
//! weakly, with metallicity, because both raise the partial pressure of the condensible.
//!
//! THE COEFFICIENTS ARE APPROXIMATE. They are anchored to the widely quoted 1 bar
//! condensation temperatures -- Fe ~1800 K, Al (corundum) ~1700 K, Ca ~1650 K, Mg
//! (forsterite/enstatite) ~1600 K, Na (Na2S) ~1200 K, K (KCl) ~1000 K -- with slopes
//! chosen to give the usual ~80-150 K per decade of pressure, and a common metallicity
//! coefficient of 0.5. They are NOT a specific published fit, and the real chemistry is a
//! network (Na2S and KCl form by reaction with H2S and HCl, not by simple vaporization),
//! so treat the resulting T_cond as good to perhaps a hundred kelvin. Replace the
//! constants here with a Lodders/Visscher fit if that accuracy matters; the framework
//! does not change.
struct MetalDonor {
  double abun;
  double chi;
  double gfac;
  double tc_a;
  double tc_b;
};
constexpr int n_metal_donor = 6;
constexpr MetalDonor metal_donor[n_metal_donor] = {
  {2.0e-6, 5.139*ev, 1.0,     13.0, 15600.0},   // Na, as Na2S,  ~1200 K at 1 bar
  {1.2e-7, 4.341*ev, 1.0,     12.5, 12500.0},   // K,  as KCl,   ~1000 K at 1 bar
  {2.2e-6, 6.113*ev, 4.0,     12.1, 19965.0},   // Ca,           ~1650 K at 1 bar
  {3.0e-6, 5.986*ev, 1.0/3.0, 12.0, 20400.0},   // Al, corundum, ~1700 K at 1 bar
  {3.8e-5, 7.646*ev, 4.0,     12.2, 19520.0},   // Mg, silicate, ~1600 K at 1 bar
  {3.2e-5, 7.902*ev, 2.4,     12.0, 21600.0},   // Fe,           ~1800 K at 1 bar
};
constexpr double metal_tc_mh = 0.5;   // rough d/d[M/H] coefficient, common to all

}  // namespace eos_cgs

//----------------------------------------------------------------------------------------
//! \struct EOSCompositionState
//! \brief everything the composition model returns for one (rho,T) point, in cgs

struct EOSCompositionState {
  double e_spec;    // specific internal energy of the gas, erg/g
  double p_spec;    // p/rho for the gas, erg/g (so that p = rho*p_spec)
  double mu;        // mean molecular weight, mass per free particle in units of m_u
  double xh2;       // fraction of H nuclei locked in H2 (0 = fully atomic)
  double xhii;      // fraction of H nuclei ionized
  double xe;        // free electrons per particle, n_e/n_tot -- what resistivity needs
};

//----------------------------------------------------------------------------------------
//! \struct EOSCompositionModel
//! \brief the analytic (rho,T) EOS, parameterised by composition and by which pieces of
//! physics are switched on. Host-only.

struct EOSCompositionModel {
  double xhyd = 0.7381;     // hydrogen mass fraction
  double yhel = 0.2485;     // helium mass fraction
  double a_metal = 16.0;    // mean atomic weight of the inert metal component
  bool include_h2 = true;   // include H2 formation, dissociation and rot/vib energy
  bool include_ion = true;  // include H and He ionization

  //! Metal ionization. OFF by default, so every existing run is unchanged. The tracked
  //! species are already inside the inert metal lump and contribute one heavy particle
  //! each whether ionized or not, so switching this on adds ONLY their electrons: mu
  //! moves by ~0.01% in a molecular atmosphere and ~0.4% in the fully ionized limit.
  //! What changes by orders of magnitude is n_e.
  bool include_metal_ion = false;
  //! Metallicity of the tracked donors, [M/H] in DEX, so that it is the same number and
  //! the same convention as the metallicity a problem generator feeds its opacity -- the
  //! two must agree or the atmosphere is opaque at one metallicity and conducting at
  //! another. Scales the electron donors only; the inert lump is still set by
  //! Z = 1 - X - Y, so adjust xhyd and yhel too for a large departure from solar.
  double metal_mh = 0.0;
  //! Per-species rainout on the condensation curves in eos_cgs::metal_donor. Each donor
  //! is removed from the gas phase below its own T_cond(p, [M/H]), which is where it
  //! stops supplying electrons. Off by default because it costs a second pass through the
  //! composition model and most runs never get cold enough to need it.
  //!
  //! Only the ELECTRON DONATION is suppressed. The condensed material is not removed from
  //! the mass or the particle count: whether it rains out of the column or stays as cloud
  //! is a transport question this local model cannot answer, and the mass involved is
  //! ~1e-4 of the gas, so leaving it in place is both simpler and the smaller error.
  bool include_metal_cond = false;
  //! Crude alternative to the above: a single temperature, below which ALL tracked metals
  //! are removed over one decade in T. Overrides include_metal_cond when positive.
  //!
  //! Real condensation is species-specific and pressure-dependent -- Fe near 1800 K,
  //! Na2S near 1200 K, KCl near 1000 K at 1 bar, all shifting with pressure -- and this
  //! single knob captures none of that structure. It exists so that a run which reaches
  //! into the condensation regime is not silently given full gas-phase abundances, which
  //! would overstate n_e by orders of magnitude there. Leave it off if the domain stays
  //! hot.
  double metal_tcond = 0.0;

  //! \fn double MetalFraction
  //! \brief metal mass fraction Z = 1 - X - Y
  double MetalFraction() const { return (1.0 - xhyd - yhel); }

  //--------------------------------------------------------------------------------------
  //! \fn void RotPartition
  //! \brief rotational partition function of H2 and the associated internal energy per
  //! molecule, including the 1:3 para:ortho nuclear spin weights.
  //!
  //! The nuclear spin degeneracy (4 in the high temperature limit) is carried HERE rather
  //! than being dropped from both sides, so that the statistical weights used in the
  //! dissociation constant below are the complete ones for H2 and for two H atoms and no
  //! convention-dependent factor is left over.
  //!
  //! Below 10*theta_rot the sum over J is short and is done directly; above it the sum is
  //! replaced by the Mulholland asymptotic expansion, which is accurate to better than
  //! 1e-8 there and avoids summing thousands of terms at high temperature (where H2 has
  //! long since dissociated and the value is irrelevant anyway).
  static void RotPartition(const double t, double &zrot, double &erot) {
    const double th = eos_cgs::theta_rot;
    if (t < 10.0*th) {
      double zsum = 0.0, esum = 0.0;
      for (int jj=0; jj<200; ++jj) {
        double gns = (jj % 2 == 0) ? 1.0 : 3.0;         // para : ortho = 1 : 3
        double ej = static_cast<double>(jj)*(jj + 1.0)*th;
        double term = gns*(2.0*jj + 1.0)*exp(-ej/t);
        zsum += term;
        esum += term*ej;
        if (jj > 4 && term < 1.0e-18*zsum) break;
      }
      zrot = zsum;
      erot = eos_cgs::kboltz*esum/zsum;
    } else {
      // Mulholland expansion for a homonuclear diatomic (symmetry number 2), times the
      // total nuclear spin degeneracy 4
      double u = th/t;
      double series = 1.0 + u/3.0 + u*u/15.0 + 4.0*u*u*u/315.0;
      double dseries = (1.0/3.0 + 2.0*u/15.0 + 12.0*u*u/315.0)*(-u/t);
      zrot = 2.0*(t/th)*series;
      // E = k T^2 dln(zrot)/dT, with zrot = C*T*series(T)
      erot = eos_cgs::kboltz*t*t*(1.0/t + dseries/series);
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn double VibPartition
  //! \brief harmonic vibrational partition function of H2, measured from the v=0 level,
  //! and the corresponding energy per molecule
  static void VibPartition(const double t, double &zvib, double &evib) {
    double u = eos_cgs::theta_vib/t;
    if (u > 300.0) {                 // exp(-u) underflows; the mode is frozen out
      zvib = 1.0;
      evib = 0.0;
    } else {
      double ex = exp(-u);
      zvib = 1.0/(1.0 - ex);
      evib = eos_cgs::kboltz*eos_cgs::theta_vib*ex/(1.0 - ex);
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn double SahaFactor
  //! \brief the (2 pi m_e k T/h^2)^(3/2) factor common to all three ionization
  //! equilibria, in cm^-3
  static double SahaFactor(const double t) {
    double hh = eos_cgs::hplanck*eos_cgs::hplanck;
    double x = 2.0*M_PI*eos_cgs::m_el*eos_cgs::kboltz*t/hh;
    return pow(x, 1.5);
  }

  //--------------------------------------------------------------------------------------
  //! \fn double DissociationConstant
  //! \brief K_D(T) = n_H^2/n_H2 for H2 <-> 2H, in cm^-3.
  //!
  //! K_D = (pi m_H k T/h^2)^(3/2) * g_H^2/(g_H2 z_rot z_vib) * exp(-chi_D/kT), with
  //! g_H = 4 (electron spin 2 times proton spin 2) and the H2 electronic weight 1; the
  //! nuclear spin degeneracy of H2 is already inside z_rot.
  static double DissociationConstant(const double t) {
    double zrot, erot, zvib, evib;
    RotPartition(t, zrot, erot);
    VibPartition(t, zvib, evib);
    double mh = eos_cgs::a_hyd*eos_cgs::m_u;
    double x = M_PI*mh*eos_cgs::kboltz*t/(eos_cgs::hplanck*eos_cgs::hplanck);
    double expo = -eos_cgs::chi_d/(eos_cgs::kboltz*t);
    if (expo < -700.0) return 0.0;                 // fully molecular; underflow guard
    return pow(x, 1.5)*(16.0/(zrot*zvib))*exp(expo);
  }

  //--------------------------------------------------------------------------------------
  //! \fn double ChargeResidual
  //! \brief n_H+ + n_He+ + 2 n_He++ - n_e at a trial electron density, i.e. the function
  //! whose root is the equilibrium electron density. Monotonically decreasing in n_e.
  //! Also returns the species densities implied by the trial value.
  void Species(const double nel, const double nh_tot, const double nhe_tot,
               const bool has_h2, const double kdis,
               const double kh, const double khe1, const double khe2,
               double &nh2, double &nh, double &nhii,
               double &nhe, double &nheii, double &nheiii,
               const double *kmet = nullptr, const double nmet_base = 0.0,
               double *nmion = nullptr, const double *fcond = nullptr) const {
    // hydrogen: n_Htot = 2 n_H2 + n_H + n_H+, with n_H+ = (K_H/n_e) n_H and
    // n_H2 = n_H^2/K_D. Solved as a quadratic in n_H, in the form that stays accurate
    // when the molecular term is negligible. H2 is then recovered from the nuclei budget
    // rather than from n_H^2/K_D, so that hydrogen nuclei are conserved exactly even
    // where K_D underflows and the atomic fraction is unresolvably small.
    double yion = (kh > 0.0) ? kh/nel : 0.0;
    double bcoef = 1.0 + yion;
    if (!has_h2) {
      nh = nh_tot/bcoef;
      nh2 = 0.0;
    } else if (kdis > 0.0) {
      double disc = bcoef*bcoef + 8.0*nh_tot/kdis;
      nh = 2.0*nh_tot/(bcoef + sqrt(disc));
      nh2 = 0.5*(nh_tot - nh*bcoef);
    } else {
      nh = 0.0;                    // K_D -> 0: everything not ionized is molecular
      nh2 = 0.5*nh_tot;
    }
    nhii = yion*nh;

    // helium: two successive ionizations off the same electron pool
    double r1 = (khe1 > 0.0) ? khe1/nel : 0.0;
    double r2 = (khe2 > 0.0) ? khe2/nel : 0.0;
    double denom = 1.0 + r1 + r1*r2;
    nhe = nhe_tot/denom;
    nheii = r1*nhe;
    nheiii = r1*r2*nhe;

    // metals: each singly ionizing off the same electron pool, and each contributing a
    // term that decreases with n_e, so the charge residual stays monotonic and the
    // bisection bracket below is still valid
    if (nmion != nullptr) {
      double sum = 0.0;
      if (kmet != nullptr) {
        for (int s=0; s<eos_cgs::n_metal_donor; ++s) {
          double rm = kmet[s]/nel;
          double fc = (fcond != nullptr) ? fcond[s] : 1.0;
          sum += eos_cgs::metal_donor[s].abun*nmet_base*fc*rm/(1.0 + rm);
        }
      }
      *nmion = sum;
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn void CondensationFactors
  //! \brief gas-phase fraction of each metal donor at a given temperature and pressure.
  //!
  //! Each species is removed below its own T_cond(p, [M/H]) from the Clausius-Clapeyron
  //! curve in eos_cgs::metal_donor. The transition is a tanh over 5% in temperature,
  //! not a step: the result is tabulated and then differenced to build the interpolation
  //! table, and a discontinuity would give node derivatives that make the bicubic Hermite
  //! patch ring.
  void CondensationFactors(const double t, const double p_cgs, double *fcond) const {
    if (metal_tcond > 0.0) {              // crude single-temperature override
      double f = log10(t/(0.1*metal_tcond));
      f = (f < 0.0) ? 0.0 : ((f > 1.0) ? 1.0 : f);
      for (int s=0; s<eos_cgs::n_metal_donor; ++s) fcond[s] = f;
      return;
    }
    const double lgp_bar = log10(p_cgs/1.0e6);          // 1 bar = 1e6 barye
    for (int s=0; s<eos_cgs::n_metal_donor; ++s) {
      double denom = eos_cgs::metal_donor[s].tc_a - lgp_bar
                   - eos_cgs::metal_tc_mh*metal_mh;
      // a non-positive denominator means the curve has run off its range of validity;
      // treat it as "no condensation" rather than producing a negative T_cond
      double tc = (denom > 0.0) ? eos_cgs::metal_donor[s].tc_b/denom : 0.0;
      fcond[s] = 0.5*(1.0 + tanh((t - tc)/(0.05*(tc > 0.0 ? tc : 1.0))));
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn EOSCompositionState Evaluate
  //! \brief composition, pressure and internal energy at (rho,T) in cgs.
  //!
  //! Condensation needs the pressure, and the pressure is an output, so it takes two
  //! passes: one at full gas-phase abundance to get p, then the real one. That is exact
  //! enough because the donors are ~1e-4 of the particles, so p is insensitive to whether
  //! they have condensed -- and this is host-side table-build code, where a factor of two
  //! costs nothing.
  EOSCompositionState Evaluate(const double rho, const double t) const {
    if (include_metal_ion && (include_metal_cond || metal_tcond > 0.0)) {
      EOSCompositionState s0 = EvaluateAt(rho, t, nullptr);
      double fcond[eos_cgs::n_metal_donor];
      CondensationFactors(t, rho*s0.p_spec, fcond);
      return EvaluateAt(rho, t, fcond);
    }
    return EvaluateAt(rho, t, nullptr);
  }

  //--------------------------------------------------------------------------------------
  //! \fn EOSCompositionState EvaluateAt
  //! \brief the model itself, at a prescribed gas-phase fraction for each metal donor
  EOSCompositionState EvaluateAt(const double rho, const double t,
                                 const double *fcond) const {
    using namespace eos_cgs;  // NOLINT(build/namespaces)

    // nuclei per unit volume
    double nh_tot = xhyd*rho/(a_hyd*m_u);
    double nhe_tot = yhel*rho/(a_hel*m_u);
    double nz = (a_metal > 0.0) ? MetalFraction()*rho/(a_metal*m_u) : 0.0;

    // equilibrium constants
    double kdis = include_h2 ? DissociationConstant(t) : 0.0;
    double sf = SahaFactor(t);
    double kh = 0.0, khe1 = 0.0, khe2 = 0.0;
    if (include_ion) {
      double bh = -chi_h/(kboltz*t);
      double b1 = -chi_he1/(kboltz*t);
      double b2 = -chi_he2/(kboltz*t);
      kh   = (bh > -700.0) ? sf*exp(bh) : 0.0;         // g factor 2*(1/2) = 1
      khe1 = (b1 > -700.0) ? 4.0*sf*exp(b1) : 0.0;     // g factor 2*(2/1) = 4
      khe2 = (b2 > -700.0) ? sf*exp(b2) : 0.0;         // g factor 2*(1/2) = 1
    }

    // metal donors. The T-dependent factors are hoisted out of the bisection below, so
    // the extra cost inside it is one divide per species per iteration.
    double kmet[eos_cgs::n_metal_donor];
    double nmet_base = 0.0;
    const bool has_metal = include_metal_ion && (nh_tot > 0.0);
    if (has_metal) {
      nmet_base = nh_tot*pow(10.0, metal_mh);
      for (int s=0; s<eos_cgs::n_metal_donor; ++s) {
        double bm = -eos_cgs::metal_donor[s].chi/(kboltz*t);
        kmet[s] = (bm > -700.0) ? eos_cgs::metal_donor[s].gfac*sf*exp(bm) : 0.0;
      }
    } else {
      for (int s=0; s<eos_cgs::n_metal_donor; ++s) kmet[s] = 0.0;
    }

    // Solve charge neutrality for the electron density by bisection in log n_e. The
    // residual n_H+ + n_He+ + 2n_He++ - n_e is positive at small n_e (where the Saha
    // ratios drive everything to the ionized side) and negative at n_e = n_max, so the
    // bracket below always contains the root. Forty decades of headroom is far more than
    // the ionized fraction can ever need to be relevant: below that the electrons make no
    // measurable contribution to mu, p or e, and the bracket floor is returned instead.
    double nmax = nh_tot + 2.0*nhe_tot;
    double nh2, nh, nhii, nhe, nheii, nheiii, nmion = 0.0;
    double nel = 0.0;
    bool any_ion = (kh > 0.0 || khe1 > 0.0);
    if (has_metal) {
      for (int s=0; s<eos_cgs::n_metal_donor; ++s) any_ion = any_ion || (kmet[s] > 0.0);
    }
    if (any_ion) {
      double lo = log(nmax) - 40.0*M_LN10;
      double hi = log(nmax);
      for (int it=0; it<80; ++it) {
        double mid = 0.5*(lo + hi);
        double ntry = exp(mid);
        Species(ntry, nh_tot, nhe_tot, include_h2, kdis, kh, khe1, khe2,
                nh2, nh, nhii, nhe, nheii, nheiii, kmet, nmet_base, &nmion, fcond);
        double res = nhii + nheii + 2.0*nheiii + nmion - ntry;
        if (res > 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
        if ((hi - lo) < 1.0e-13) break;
      }
      nel = exp(0.5*(lo + hi));
    } else {
      nel = 1.0e-40*nmax;
    }
    Species(nel, nh_tot, nhe_tot, include_h2, kdis, kh, khe1, khe2,
            nh2, nh, nhii, nhe, nheii, nheiii, kmet, nmet_base, &nmion, fcond);
    // the free electrons are whatever the ions supply, not the bisection iterate
    nel = nhii + nheii + 2.0*nheiii + nmion;

    // total free particles, and the pressure and mean molecular weight that follow
    double ntot = nh2 + nh + nhii + nhe + nheii + nheiii + nz + nel;
    double pres = ntot*kboltz*t;

    // internal energy: translational, then the H2 internal modes, then chemical
    double eint = 1.5*pres;
    if (nh2 > 0.0) {
      double zrot, erot, zvib, evib;
      RotPartition(t, zrot, erot);
      VibPartition(t, zvib, evib);
      eint += nh2*(erot + evib);
    }
    // Chemical energy, measured from a zero point of "all hydrogen in ground state H2,
    // all helium neutral". With H2 switched off there is no molecular reservoir to
    // measure from, so the zero point moves to atomic H and the dissociation energy drops
    // out; that leaves a pure monatomic ideal gas when ionization is off too, which is
    // what the table machinery is validated against.
    double e_href = include_h2 ? 0.5*chi_d : 0.0;
    eint += nh*e_href + nhii*(e_href + chi_h)
          + nheii*chi_he1 + nheiii*(chi_he1 + chi_he2);
    // metal ionization energy, measured from neutral. Worth ~1e-5 of the total, far too
    // small to matter for p or e, but it is what makes c_v consistent across the metal
    // ionization zone -- and c_v is differenced out of this model to build the table.
    if (has_metal) {
      for (int s=0; s<eos_cgs::n_metal_donor; ++s) {
        double rm = kmet[s]/nel;
        double fc = (fcond != nullptr) ? fcond[s] : 1.0;
        eint += eos_cgs::metal_donor[s].abun*nmet_base*fc*(rm/(1.0 + rm))
                *eos_cgs::metal_donor[s].chi;
      }
    }

    EOSCompositionState s;
    s.e_spec = eint/rho;
    s.p_spec = pres/rho;
    s.mu = rho/(ntot*m_u);
    s.xh2 = (nh_tot > 0.0) ? 2.0*nh2/nh_tot : 0.0;
    s.xhii = (nh_tot > 0.0) ? nhii/nh_tot : 0.0;
    s.xe = (ntot > 0.0) ? nel/ntot : 0.0;
    return s;
  }
};

#endif // EOS_EOS_COMPOSITION_HPP_
