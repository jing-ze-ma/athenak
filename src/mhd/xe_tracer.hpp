#ifndef MHD_XE_TRACER_HPP_
#define MHD_XE_TRACER_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file xe_tracer.hpp
//! \brief advected electron fraction (xetracer-1101), <mhd>/xe_tracer = true.
//!
//! One extra MHD passive scalar (the LAST one, index nmhd+nscalars-1) carries
//! z = x_e = n_e/n_tot, the same quantity ResistivityEOS() and SetAmbipolar() otherwise
//! take from the local Saha table (EOS_Data::ElectronFractionKelvin).  With the key on,
//! both use z instead.  Between the advection steps z relaxes towards the local Saha
//! value n_eq = x_eq n_tot by single-carrier radiative recombination
//!
//!     dn_e/dt = alpha(T) (n_eq^2 - n_e^2),
//!
//! operator split once per cycle, pointwise, backward Euler over the cycle dt with
//! n_tot frozen.  Backward Euler is a quadratic in n_e^{n+1}; its positive root is
//!
//!     n_e^{n+1} = 2c/(1 + sqrt(1 + 4 a c)),   a = alpha dt,  c = n_e^n + a n_eq^2,
//!
//! (the cancellation-free form), which is exact for the discrete scheme, unconditionally
//! stable and monotone, and returns n_eq for a -> infinity (dayside, Damkoehler >> 1).
//! Afterwards z is capped at max(x_eq, x_don), x_don = the gas-phase metal donors fully
//! singly ionized (the same abundances, [M/H] and condensation curves the composition
//! table was built with, eos_composition.hpp), so condensation removes the carriers.
//!
//! alpha = Verner & Ferland (1996) radiative recombination of Na+ (rrfit, Na-like Na I;
//! the same fit the offline nightxe_1101/chem.py uses), times <mhd>/xe_alpha_fac.  K+
//! (2.762e-13 (T/1e4)^-0.80, rrfit power law) is within 30 % of it at 500-3000 K, and
//! Na has 17x the abundance of K.  Night carriers are Na+ and K+ (closed-shell cores:
//! dielectronic recombination negligible).

#include <cmath>
#include "athena.hpp"
#include "eos/eos_composition.hpp"

namespace xetracer {

constexpr int kNDonor = eos_cgs::n_metal_donor;

//! everything the source kernel needs, captured by value (no host pointers)
struct XeTracerParams {
  Real alpha_fac = 1.0;     // <mhd>/xe_alpha_fac
  Real time_cgs = 1.0;      // s per code time unit
  Real xh = 0.7381;         // <mhd>/eos_xh (hydrogen mass fraction)
  Real metal_mh = 0.0;      // [M/H] of the donors
  Real metal_tcond = 0.0;   // single-temperature rainout override
  bool metal_cond = false;  // per-species condensation curves
  bool metal_ion = false;   // donors present at all
  Real abun[kNDonor];
  Real tca[kNDonor], tcb[kNDonor], tcc[kNDonor], tcd[kNDonor];
};

//! fill the donor tables from eos_cgs::metal_donor (host)
inline void FillDonors(XeTracerParams &p) {
  for (int s=0; s<kNDonor; ++s) {
    p.abun[s] = eos_cgs::metal_donor[s].abun;
    p.tca[s] = eos_cgs::metal_donor[s].tc_a;
    p.tcb[s] = eos_cgs::metal_donor[s].tc_b;
    p.tcc[s] = eos_cgs::metal_donor[s].tc_c;
    p.tcd[s] = eos_cgs::metal_donor[s].tc_d;
  }
}

//! radiative recombination coefficient of Na+ [cm^3/s], Verner & Ferland 1996
KOKKOS_INLINE_FUNCTION
Real AlphaRR(const Real t) {
  const Real a = 5.641e-12, b = 0.1749, t0 = 3.077e2, t1 = 2.617e6;
  const Real s0 = sqrt(t/t0), s1 = sqrt(t/t1);
  return a/(s0*pow(1.0 + s0, 1.0 - b)*pow(1.0 + s1, 1.0 + b));
}

//! electron density [cm^-3] if every gas-phase metal donor were singly ionized.  Same
//! condensation factors as EOSCompositionModel::CondensationFactors.
KOKKOS_INLINE_FUNCTION
Real DonorCap(const XeTracerParams &p, const Real t, const Real p_cgs, const Real nh) {
  if (!p.metal_ion) return 0.0;
  const Real nbase = nh*pow(10.0, p.metal_mh);
  Real sum = 0.0;
  if (p.metal_tcond > 0.0) {
    Real f = log10(t/(0.1*p.metal_tcond));
    f = (f < 0.0) ? 0.0 : ((f > 1.0) ? 1.0 : f);
    for (int s=0; s<kNDonor; ++s) sum += p.abun[s]*f;
    return nbase*sum;
  }
  Real lgp = log10(p_cgs/1.0e6);
  lgp = (lgp < -8.0) ? -8.0 : ((lgp > 3.0) ? 3.0 : lgp);
  for (int s=0; s<kNDonor; ++s) {
    Real fc = 1.0;
    if (p.metal_cond) {
      const Real den = p.tca[s] - p.tcb[s]*lgp + p.tcd[s]*lgp*lgp - p.tcc[s]*p.metal_mh;
      const Real tc = (den > 0.0) ? 1.0e4/den : 0.0;
      fc = 0.5*(1.0 + tanh((t - tc)/(0.05*(tc > 0.0 ? tc : 1.0))));
    }
    sum += p.abun[s]*fc;
  }
  return nbase*sum;
}

//! backward-Euler step of dn/dt = alpha (neq^2 - n^2) over a = alpha dt (exact root)
KOKKOS_INLINE_FUNCTION
Real Relax(const Real n, const Real neq, const Real a) {
  const Real np = (n > 0.0) ? n : 0.0;
  const Real c = np + a*neq*neq;
  return 2.0*c/(1.0 + sqrt(1.0 + 4.0*a*c));
}

}  // namespace xetracer

#endif  // MHD_XE_TRACER_HPP_
