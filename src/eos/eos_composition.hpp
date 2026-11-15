#ifndef EOS_EOS_COMPOSITION_HPP_
#define EOS_EOS_COMPOSITION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_composition.hpp
//! \brief analytic equation of state of a partially dissociated, partially ionized H/He
//! mixture in (rho, T), cgs.  Host-only: evaluated once per node to build the table
//! (eos_table.{hpp,cpp}).
//!
//! Species: H2, H, H+, He, He+, He++, free electrons and a metal component (particles and
//! optional singly ionizing electron donors).  Composition from Saha for H, He, He+ and
//! the dissociation equilibrium H2 <-> 2H, closed by charge neutrality (one root find on
//! n_e). Internal energy: translational (3/2) nkT of all free particles, H2 rotation and
//! vibration from its partition functions, and the chemical energy measured from all H in
//! ground-state H2 and all He neutral.
//!
//! Not included: Coulomb corrections, degeneracy, pressure ionization and excited bound
//! states, so the model is valid for weakly coupled, non-degenerate gas (planetary and
//! stellar envelopes).  Radiation is added on top of the table (eos_table.hpp).

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
//! \brief one singly ionizing metal as an electron source.  Below ~4000 K the free
//! electrons come from these 4-8 eV species, not from hydrogen.  abun is the solar
//! abundance by number relative to H nuclei, gfac = 2 g_ion/g_neutral from ground-state
//! weights (a poorer approximation for Fe, whose neutral partition function is ~1.5-2x
//! larger near 3000-5000 K).  tc_a..tc_d give the condensation curve of the phase that
//! removes the element from the gas,
//!
//!     10^4 / T_cond = tc_a - tc_b log10(p_T/bar) + tc_d (log10 p_T)^2 - tc_c [M/H]
//!
//! with p_T the total pressure.  Literature coefficients, per species:
//!
//!   Na  Na2S       Morley et al. 2012, ApJ 756, 172 (arXiv:1206.4313), Table 1
//!                    10^4/T_cond(Na2S)  = 10.045 - 0.72  log p_t - 1.08  [Fe/H]
//!   K   KCl        Morley et al. 2012, ApJ 756, 172 (arXiv:1206.4313), Table 1
//!                    10^4/T_cond(KCl)   = 12.479 - 0.879 log p_t - 0.879 [Fe/H]
//!   Ca  CaAl4O7    Wakeford et al. 2017, MNRAS 464, 4247 (arXiv:1610.03325)
//!                    10^4/T_cond(CaAl4O7) = 4.990 - 0.2394 log P_T
//!                                            + 1.398e-3 (log P_T)^2 - 0.595 [Fe/H]
//!   Al  Al2O3      Wakeford et al. 2017, MNRAS 464, 4247 (arXiv:1610.03325)
//!                    10^4/T_cond(Al2O3)   = 5.014 - 0.2179 log P_T
//!                                            + 2.264e-3 (log P_T)^2 - 0.580 [Fe/H]
//!   Mg  Mg2SiO4    Visscher, Lodders & Fegley 2010, ApJ 716, 1060 (arXiv:1001.3639),
//!                  eq 18
//!                    10^4/T_cond(Mg2SiO4) = 5.89 - 0.37 log P_T - 0.73 [Fe/H]
//!   Fe  Fe metal   Visscher, Lodders & Fegley 2010, ApJ 716, 1060 (arXiv:1001.3639) eq 2
//!                    10^4/T_cond(Fe)      = 5.44 - 0.48 log P_T - 0.48 [Fe/H]
//!
//! Visscher+2010 is valid for 800-2500 K and [Fe/H] <= +0.5.  Ca (grossite) and Mg
//! (forsterite) are removed slightly early; below ~1700 K n_e is set by Na and K anyway.
struct MetalDonor {
  double abun;
  double chi;
  double gfac;
  double tc_a;
  double tc_b;
  double tc_c;
  double tc_d;
};
constexpr int kNumMetalDonor = 6;
constexpr MetalDonor metal_donor[kNumMetalDonor] = {
  // abundance  chi        2g+/g0     tc_a   tc_b    tc_c   tc_d   phase / source / 1 bar
  {2.0e-6, 5.139*ev, 1.0,     10.045, 0.72,   1.08,  0.0},      // Na Na2S    Mo12   995 K
  {1.2e-7, 4.341*ev, 1.0,     12.479, 0.879,  0.879, 0.0},      // K  KCl     Mo12   801 K
  {2.2e-6, 6.113*ev, 4.0,      4.990, 0.2394, 0.595, 1.398e-3}, // Ca CaAl4O7 Wa17  2004 K
  {3.0e-6, 5.986*ev, 1.0/3.0,  5.014, 0.2179, 0.580, 2.264e-3}, // Al Al2O3   Wa17  1994 K
  {3.8e-5, 7.646*ev, 4.0,      5.89,  0.37,   0.73,  0.0},      // Mg Mg2SiO4 Vi10  1698 K
  {3.2e-5, 7.902*ev, 2.4,      5.44,  0.48,   0.48,  0.0},      // Fe metal   Vi10  1838 K
  // Mo12 = Morley et al. 2012; Wa17 = Wakeford et al. 2017; Vi10 = Visscher, Lodders &
  // Fegley 2010 (equations above)
};

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
//! \brief the analytic (rho, T) EOS, parameterised by composition and physics switches

struct EOSCompositionModel {
  double xhyd = 0.7381;     // hydrogen mass fraction
  double yhel = 0.2485;     // helium mass fraction
  double a_metal = 16.0;    // mean atomic weight of the inert metal component
  bool include_h2 = true;   // include H2 formation, dissociation and rot/vib energy
  bool include_ion = true;  // include H and He ionization

  //! metal ionization (default off): adds the donors' electrons only (they are already
  //! particles of the inert metal lump), so mu barely changes while n_e changes by orders
  //! of magnitude
  bool include_metal_ion = false;
  //! [M/H] of the donors in dex (same convention as the opacity metallicity); the inert
  //! lump stays Z = 1 - X - Y
  double metal_mh = 0.0;
  //! per-species condensation on the curves of eos_cgs::metal_donor (default off): only
  //! the electron donation is removed below T_cond, the mass and particle count are kept
  bool include_metal_cond = false;

  //! \fn double MetalFraction
  //! \brief metal mass fraction Z = 1 - X - Y
  double MetalFraction() const { return (1.0 - xhyd - yhel); }

  //--------------------------------------------------------------------------------------
  //! \fn void RotPartition
  //! \brief rotational partition function of H2 and its internal energy per molecule,
  //! with the para:ortho nuclear spin weights included (complete weights in K_D).  Direct
  //! sum below 10 theta_rot, Mulholland expansion above.
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
        for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) {
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
  //! \brief gas-phase fraction of each metal donor at (T, p): a tanh over 5% in T around
  //! T_cond, smooth so the differenced node derivatives do not ring
  void CondensationFactors(const double t, const double p_cgs, double *fcond) const {
    // clamp p to 1e-8..1e3 bar, where the quadratic fits are valid
    double lgp = log10(p_cgs/1.0e6);                    // 1 bar = 1e6 barye
    lgp = (lgp < -8.0) ? -8.0 : ((lgp > 3.0) ? 3.0 : lgp);
    for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) {
      const eos_cgs::MetalDonor &m = eos_cgs::metal_donor[s];
      double denom = m.tc_a - m.tc_b*lgp + m.tc_d*lgp*lgp - m.tc_c*metal_mh;
      // a non-positive denominator means the curve has run off its range of validity;
      // treat that as "no condensation" rather than produce a negative T_cond
      double tc = (denom > 0.0) ? 1.0e4/denom : 0.0;
      fcond[s] = 0.5*(1.0 + tanh((t - tc)/(0.05*(tc > 0.0 ? tc : 1.0))));
    }
  }

  //--------------------------------------------------------------------------------------
  //! \fn EOSCompositionState Evaluate
  //! \brief composition, pressure and internal energy at (rho, T) in cgs; with
  //! condensation, a first pass at full gas-phase abundance supplies p for T_cond
  EOSCompositionState Evaluate(const double rho, const double t) const {
    if (include_metal_ion && include_metal_cond) {
      EOSCompositionState s0 = EvaluateAt(rho, t, nullptr);
      double fcond[eos_cgs::kNumMetalDonor];
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
    double kmet[eos_cgs::kNumMetalDonor];
    double nmet_base = 0.0;
    const bool has_metal = include_metal_ion && (nh_tot > 0.0);
    if (has_metal) {
      nmet_base = nh_tot*pow(10.0, metal_mh);
      for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) {
        double bm = -eos_cgs::metal_donor[s].chi/(kboltz*t);
        kmet[s] = (bm > -700.0) ? eos_cgs::metal_donor[s].gfac*sf*exp(bm) : 0.0;
      }
    } else {
      for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) kmet[s] = 0.0;
    }

    // charge neutrality by bisection in log n_e: the residual is positive at small n_e
    // and negative at n_max; the 40-decade floor is returned where electrons do not
    // matter
    double nmax = nh_tot + 2.0*nhe_tot;
    double nh2, nh, nhii, nhe, nheii, nheiii, nmion = 0.0;
    double nel = 0.0;
    bool any_ion = (kh > 0.0 || khe1 > 0.0);
    if (has_metal) {
      for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) any_ion = any_ion || (kmet[s] > 0.0);
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
    // chemical energy from the zero point all H in ground-state H2, all He neutral
    // (atomic H without H2, so the gas is monatomic ideal with ionization off as well)
    double e_href = include_h2 ? 0.5*chi_d : 0.0;
    eint += nh*e_href + nhii*(e_href + chi_h)
          + nheii*chi_he1 + nheiii*(chi_he1 + chi_he2);
    // metal ionization energy from neutral (keeps c_v consistent across the metal
    // ionization zone)
    if (has_metal) {
      for (int s=0; s<eos_cgs::kNumMetalDonor; ++s) {
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
