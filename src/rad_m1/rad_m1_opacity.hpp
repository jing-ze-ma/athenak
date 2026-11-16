#ifndef RAD_M1_RAD_M1_OPACITY_HPP_
#define RAD_M1_RAD_M1_OPACITY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_opacity.hpp
//! \brief the opacity laws of <rad_m1> as device-inline functions.  All return opacities
//! per unit mass in code units; the caller stores rho*kappa in RadiationM1::opac.

#include <math.h>
#include "athena.hpp"
#include "rad_m1/rad_m1.hpp"
#include "diffusion/conduction.hpp"   // RosselandTable: ONE tabulated-opacity lookup

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1TableOpacities
//! \brief <rad_m1>/opacity = table: stellar Rosseland and Planck means looked up in
//! Kelvin and g/cm^3 on the shared (log10 T, log10 rho) grid of M1OpacTab.  The Rosseland
//! total already contains electron scattering, so kappa_F = kappa_R with kappa_s = 0; the
//! Planck table (absorption only) gives kappa_P = kappa_E.  With vet_scatter (tab.scat !=
//! 0) the total is split into of = kappa_R - kappa_e and os = kappa_e (Ma+2026 eq. 7).
//! \param[in]  tab   the two tables and their unit conversions
//! \param[in]  d, t  density and temperature of the cell, code units
//! \param[out] op,oe,of,os  the four opacities per unit mass, code units

KOKKOS_INLINE_FUNCTION
void M1TableOpacities(const M1OpacTab &tab, const Real d, const Real t,
                      Real &op, Real &oe, Real &of, Real &os) {
  const Real tk = t*tab.tunit;
  const Real dc = d*tab.dunit;
  // both tables share the grid: logarithms, bisections and weights are formed once
  const Real x = log10(tk), y = log10(dc);
  const int nT = tab.nT, nD = tab.nD;
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (!(x > tab.lT(0))) {
    i = 0; fx = 0.0;
  } else if (x >= tab.lT(nT-1)) {
    i = nT-2; fx = 1.0;
  } else {
    int hi = nT-2;
    while (i < hi) {
      const int mid = (i + hi + 1)/2;
      if (tab.lT(mid) <= x) {
        i = mid;
      } else {
        hi = mid - 1;
      }
    }
    fx = (x - tab.lT(i))/(tab.lT(i+1) - tab.lT(i));
  }
  if (!(y > tab.lD(0))) {
    j = 0; fy = 0.0;
  } else if (y >= tab.lD(nD-1)) {
    j = nD-2; fy = 1.0;
  } else {
    int hj = nD-2;
    while (j < hj) {
      const int mid = (j + hj + 1)/2;
      if (tab.lD(mid) <= y) {
        j = mid;
      } else {
        hj = mid - 1;
      }
    }
    fy = (y - tab.lD(j))/(tab.lD(j+1) - tab.lD(j));
  }
  const Real lr = (1.0-fx)*((1.0-fy)*tab.kr(i,j) + fy*tab.kr(i,j+1))
                +      fx *((1.0-fy)*tab.kr(i+1,j) + fy*tab.kr(i+1,j+1));
  const Real lp = (1.0-fx)*((1.0-fy)*tab.kp(i,j) + fy*tab.kp(i,j+1))
                +      fx *((1.0-fy)*tab.kp(i+1,j) + fy*tab.kp(i+1,j+1));
  const Real kr = pow(10.0, lr);
  const Real kp = pow(10.0, lp);
  op = kp*tab.kunit;
  oe = op;
  of = kr*tab.kunit;
  os = 0.0;
  // vet_scatter: split the total into kappa_a + kappa_e, guarded as in Ma+2026; scat = 2
  // (EOS kappa_e) starts from the constant split and M1ScatterEos refines it
  if (tab.scat != 0) {
    const Real kt = of;
    os = M1ScatterGuard(kt, tab.kes, tab.kfl);
    of = kt - os;
  }
  // <rad_m1>/opac_abs_rho_max: below this density the gas only scatters (0 = off)
  if (d < tab.amask_rho) {
    op = 0.0;
    oe = 0.0;
    of = tab.amask_kes;
    os = 0.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1ScatterEos
//! \brief vet_scatter_kappa_e = eos: split kappa_R = of + os with the
//! ionisation-dependent electron-scattering opacity
//!     kappa_e = sigma_T n_e/rho = (sigma_T/m_u) x_e/mu,  x_e = n_e/n_tot,
//! from the general EOS composition, guarded as in Ma+2026.  No-op unless tab.scat == 2;
//! d, t in code units.
template <class EosT>
KOKKOS_INLINE_FUNCTION
void M1ScatterEos(const M1OpacTab &tab, const EosT &eos, const Real d, const Real t,
                  Real &of, Real &os) {
  if (tab.scat != 2) return;
  const Real kt = of + os;
  const Real xe = eos.ElectronFraction(d, 0.0, t);
  const Real mu = eos.MeanMolecularWeight(d, 0.0, t);
  os = M1ScatterGuard(kt, tab.kxe*xe/fmax(mu, 1.0e-30), tab.kfl);
  of = kt - os;
}

//----------------------------------------------------------------------------------------
//! \fn M1UserOpacity
//! \brief opacity law of <rad_m1>/opacity = user, inline because device lambdas cannot
//! call host function pointers; a problem with its own law edits this function.
//! \param[in]  d, t        density and temperature of the cell, code units
//! \param[in]  kp,ke,kf,ks the <rad_m1> constants, as defaults
//! \param[out] op,oe,of,os the four opacities per unit mass

KOKKOS_INLINE_FUNCTION
void M1UserOpacity(const Real d, const Real t, const Real kp, const Real ke,
                   const Real kf, const Real ks, const Real rho_ref, const Real t_ref,
                   const Real aa, const Real bb,
                   Real &op, Real &oe, Real &of, Real &os) {
  // default: a peaked kappa(T), the four constants times
  //   (d/rho_ref)^a [1 + b exp(-ln^2(t/t_ref)/(2 s^2))],  s = 0.2;
  // b = 0 (opac_b default) and a = 0 give the constant law
  const Real fd = (aa == 0.0) ? 1.0 : pow(fmax(d, 0.0)/rho_ref, aa);
  Real fb = 1.0;
  if (bb != 0.0) {
    const Real x = log(fmax(t, 1.0e-300)/t_ref);
    fb = 1.0 + bb*exp(-x*x/(2.0*0.04));
  }
  const Real sc = fd*fb;
  op = kp*sc;
  oe = ke*sc;
  of = kf*sc;
  os = ks*sc;
}

//----------------------------------------------------------------------------------------
//! \fn M1Opacities
//! \brief the four opacities per unit mass at (d, t) under the law otype:
//!   const     kappa = kappa0
//!   powerlaw  kappa = kappa0 (d/rho_ref)^a (t/t_ref)^b, one pair of exponents for all
//! four   user      M1UserOpacity above

KOKKOS_INLINE_FUNCTION
void M1Opacities(const int otype, const Real d, const Real t,
                 const Real kp, const Real ke, const Real kf, const Real ks,
                 const Real rho_ref, const Real t_ref, const Real aa, const Real bb,
                 Real &op, Real &oe, Real &of, Real &os) {
  if (otype == M1_OPAC_POWERLAW) {
    Real fd = (aa == 0.0) ? 1.0 : pow(fmax(d, 0.0)/rho_ref, aa);
    Real ft = (bb == 0.0) ? 1.0 : pow(fmax(t, 0.0)/t_ref, bb);
    Real sc = fd*ft;
    op = kp*sc;
    oe = ke*sc;
    of = kf*sc;
    os = ks*sc;
  } else if (otype == M1_OPAC_USER) {
    M1UserOpacity(d, t, kp, ke, kf, ks, rho_ref, t_ref, aa, bb, op, oe, of, os);
  } else {
    op = kp;
    oe = ke;
    of = kf;
    os = ks;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1CellOpacities
//! \brief the four opacities per unit mass at (d, t) of the <rad_m1> law: the tables with
//! the EOS electron scattering (M1ScatterEos), or M1Opacities

template <class EosT>
KOKKOS_INLINE_FUNCTION
void M1CellOpacities(const int otype, const M1OpacTab &ot, const EosT &eos, const Real d,
                     const Real t, const Real kp, const Real ke, const Real kf,
                     const Real ks, const Real rho_ref, const Real t_ref, const Real aa,
                     const Real bb, Real &op, Real &oe, Real &of, Real &os) {
  if (otype == M1_OPAC_TABLE) {
    M1TableOpacities(ot, d, t, op, oe, of, os);
    M1ScatterEos(ot, eos, d, t, of, os);
  } else {
    M1Opacities(otype, d, t, kp, ke, kf, ks, rho_ref, t_ref, aa, bb, op, oe, of, os);
  }
}

} // namespace radm1
#endif // RAD_M1_RAD_M1_OPACITY_HPP_
