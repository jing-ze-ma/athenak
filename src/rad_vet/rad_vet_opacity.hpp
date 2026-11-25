#ifndef RAD_VET_RAD_VET_OPACITY_HPP_
#define RAD_VET_RAD_VET_OPACITY_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_opacity.hpp
//! \brief <rad_vet> opacity laws (device inline): opacities per unit mass, code units;
//! the caller stores rho*kappa in RadiationVet::opac.

#include <math.h>
#include "athena.hpp"
#include "rad_vet/rad_vet.hpp"
#include "diffusion/conduction.hpp"   // RosselandTable: ONE tabulated-opacity lookup

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn M1TableOpacities
//! \brief opacity = table: Rosseland and Planck means, bilinear in (log10 T, log10 rho)
//! (K, g/cm^3). kappa_F = kappa_R (includes e-scattering), kappa_P = kappa_E; with
//! vet_scatter (tab.scat != 0): of = kappa_R - kappa_e, os = kappa_e (Ma+2026 eq. 7).
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
  // vet_scatter: split into kappa_a + kappa_e; scat = 2 is refined by M1ScatterEos
  if (tab.scat != 0) {
    const Real kt = of;
    os = M1ScatterGuard(kt, tab.kes, tab.kfl);
    of = kt - os;
  }
  // <rad_vet>/opac_abs_rho_max: below this density the gas only scatters (0 = off)
  if (d < tab.amask_rho) {
    op = 0.0;
    oe = 0.0;
    of = tab.amask_kes;
    os = 0.0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn M1ScatterEos
//! \brief vet_scatter_kappa_e = eos: split kappa_R = of + os with the EOS electron
//! scattering kappa_e = (sigma_T/m_u) x_e/mu, x_e = n_e/n_tot. No-op unless scat == 2.
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
//! \brief opacity = user law (inline: device lambdas cannot call host function
//! pointers); a problem with its own law edits this function.
//! \param[in]  d, t        density and temperature of the cell, code units
//! \param[in]  kp,ke,kf,ks the <rad_vet> constants, as defaults
//! \param[out] op,oe,of,os the four opacities per unit mass

KOKKOS_INLINE_FUNCTION
void M1UserOpacity(const Real d, const Real t, const Real kp, const Real ke,
                   const Real kf, const Real ks, const Real rho_ref, const Real t_ref,
                   const Real aa, const Real bb,
                   Real &op, Real &oe, Real &of, Real &os) {
  // default: constants times (d/rho_ref)^a [1 + b exp(-ln^2(t/t_ref)/(2 s^2))], s = 0.2
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
//!   const: kappa0; powerlaw: kappa0 (d/rho_ref)^a (t/t_ref)^b; user: M1UserOpacity

KOKKOS_INLINE_FUNCTION
void M1Opacities(const int otype, const Real d, const Real t,
                 const Real kp, const Real ke, const Real kf, const Real ks,
                 const Real rho_ref, const Real t_ref, const Real aa, const Real bb,
                 Real &op, Real &oe, Real &of, Real &os) {
  if (otype == VET_OPAC_POWERLAW) {
    Real fd = (aa == 0.0) ? 1.0 : pow(fmax(d, 0.0)/rho_ref, aa);
    Real ft = (bb == 0.0) ? 1.0 : pow(fmax(t, 0.0)/t_ref, bb);
    Real sc = fd*ft;
    op = kp*sc;
    oe = ke*sc;
    of = kf*sc;
    os = ks*sc;
  } else if (otype == VET_OPAC_USER) {
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
//! \brief the four opacities per unit mass at (d, t): tables (+ M1ScatterEos) or
//! M1Opacities

template <class EosT>
KOKKOS_INLINE_FUNCTION
void M1CellOpacities(const int otype, const M1OpacTab &ot, const EosT &eos, const Real d,
                     const Real t, const Real kp, const Real ke, const Real kf,
                     const Real ks, const Real rho_ref, const Real t_ref, const Real aa,
                     const Real bb, Real &op, Real &oe, Real &of, Real &os) {
  if (otype == VET_OPAC_TABLE) {
    M1TableOpacities(ot, d, t, op, oe, of, os);
    M1ScatterEos(ot, eos, d, t, of, os);
  } else {
    M1Opacities(otype, d, t, kp, ke, kf, ks, rho_ref, t_ref, aa, bb, op, oe, of, os);
  }
}

} // namespace radvet
#endif // RAD_VET_RAD_VET_OPACITY_HPP_
