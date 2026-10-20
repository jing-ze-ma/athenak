#ifndef RAD_M1_RAD_M1_VETSRC_HPP_
#define RAD_M1_RAD_M1_VETSRC_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetsrc.hpp
//! \brief the relaxed / blended VET source of one cell, shared by every formal solution
//! (vet_col, vet_col_lat, vet_gd, vet_sc): vet_source = relaxed and
//! time2_vet_source = blend (unify-sc-gd-1012 #3).

#include "athena.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn VetRelaxedSource
//! unify-sc-gd-1012 (#3): ONE copy for vet_col, vet_col_lat / vet_gd and vet_sc (was
//! VcolRelaxedSource / VlatRelaxedSource, term for term the same).
//! \brief vet_col_source = relaxed: the thermal source eth a T*^4 + (1 - eth) E* of one
//! cell, T* and E* the state after a LOCAL backward-Euler gas-radiation exchange over the
//! step dt from the build state (rho, T_b, E_b):
//!   E* - E_b = chat dt (rho kappa_P a T*^4 - rho kappa_E E*)
//!   rho e(T*) - rho e(T_b) = -(c/chat) (E* - E_b)
//!   rho e(T*) + c dt rho kappa_P a T*^4/(1+x) = rho e(T_b) + c dt rho kappa_E E_b/(1+x)
//! with x = chat dt rho kappa_E, by the safeguarded root find of the implicit gas solve.
//! WHY: the tensor is lagged, built from the start-of-step state, where T_b is the gas
//! temperature after the operator-split hydro step (gas-only EOS).  Where the exchange is
//! stiff (dt/t_eq >> 1) that T_b is a transient the solve removes at once: in a
//! radiation-dominated cell it carries the gas-only adiabatic response to the step's
//! compression, (gamma_gas - 1) drho/rho, of which the coupled state keeps only the gas
//! heat-capacity share ~ rho c_v/(4 a T^3) (~1/800 at P_gas/P = 0.01).  The per-column
//! (laterally homogeneous) formal solution turned that transient into an Eddington factor
//! perturbation on cell scales wherever the cells are optically thin to thick (tau_cell
//! ~ 0.1-1), and a transverse checkerboard of it drives grid-scale radial flows through
//! -(c/chi) E d(f)/dr (the he_presn wedge instability, 0929).  Where the exchange is weak
//! (x, c dt rho kappa_P small) T* -> T_b and E* -> E_b: the old source.
//!
//! vetsrc-1020, vet_source = relaxed_esc (qe > 0; qe = 0 is relaxed, bitwise): the same
//! exchange with the cell's photon ESCAPE to its neighbours, which hold the radiation at
//! the current J ~ E_b over the step:
//!   E* - E_b = chat dt (rho kappa_P a T*^4 - rho kappa_E E*) - y (E* - E_b)
//!   rho e(T*) + c dt rho kappa_P a T*^4 (1+y)/(1+x+y) = rho e(T_b) + c dt rho kappa_E
//!               E_b (1+y)/(1+x+y),      E* = [(1+y) E_b + chat dt rho kappa_P a T*^4]/(1+x+y)
//! y = chat dt qe, qe = 1/(h (2 + chi h/3)) (VetEscapeRate) the escape rate of a cell of
//! width h along the column axis, per unit chat: c/(2h) free streaming out of a thin
//! slab, c/(3 chi h^2) ~ the diffusion rate pi^2 c/(3 chi h^2) of the cell-scale mode in a
//! thick one.  x/y = tau_a (2 + tau/3), tau = chi h, tau_a = rho kappa_E h: the cell's
//! thermalization depth, ~ (h/L_th)^2 with L_th = 1/sqrt(3 chi rho kappa_a) when thick.
//! Thin or effectively thin (scattering-dominated, h < L_th): E* -> E_b, the formal
//! solution keeps the NON-LOCAL J (the NLTE atmosphere, nlte2/nltek: relaxed drove E*
//! -> B(T) there, +10-13 % in tau < 1); thick (h >> L_th): y << x, relaxed.  T* keeps
//! relaxing wherever the gas exchange is stiff (its rate c dt rho kappa_P 4 a T^3
//! (1+y)/((1+x+y) rho c_v) is >= that of relaxed), so the T_b transient of the 0929
//! wedge instability is still removed; held gas (huge rho c_v) gives T* = T_b and the
//! gas source exactly in the thin limit.

template <class EosT>
KOKKOS_INLINE_FUNCTION
Real VetRelaxedSourceQ(const EosT &eos, const Real d, const Real tb, const Real eb,
                       const Real rkp, const Real rke, const Real eth, const Real ar,
                       const Real cl, const Real ch, const Real dt,
                       const int sf, const Real chi, const Real rkes,
                       const Real qe = 0.0) {
  // vetsrc-1020 relaxed_esc: y = chat dt qe (the escape over the step, see above)
  const bool esc = (qe > 0.0);
  const Real y = ch*dt*qe;
  const Real ix = esc ? (1.0 + y)/(1.0 + ch*dt*rke + y) : 1.0/(1.0 + ch*dt*rke);
  M1EosDirect<EosT> th{eos};
  Real ee, cv;
  th(d, tb, ee, cv);
  Real ts = tb;
  bool ok = true;
  // vetsrc-1020: cl = 0 from the callers when the gas is a PRESCRIBED background
  // (<rad_m1>/gas_feedback = false: the solve does not heat it), T* = T_b
  if (cl > 0.0) {
    (void) M1ImplTemperatureT(th, d, tb, ee, cl*dt*rkp*ar*ix, cl*dt*rke*eb*ix, ts, ok);
    if (!(ok && ts > 0.0)) {ts = tb;}
  }
  const Real t2 = ts*ts;
  const Real t4 = t2*t2;
  const Real es = esc ? ((1.0 + y)*eb + ch*dt*rkp*ar*t4)/(1.0 + ch*dt*rke + y)
                      : (eb + ch*dt*rkp*ar*t4)*ix;
  // vet_scatter: B* of the relaxed gas with (sf = 1, 2; vet_scatter_j = current) the
  // CURRENT J = eb (E^n or the Picard iterate) or (sf = 3, 4; vet_scatter_j = relaxed)
  // the relaxed E*, the consistent end-of-step pair (S - E* = (E* - eb)/(c dt chi))
  if (sf != 0) {
    return VetScatterMix((sf > 2) ? (sf - 2) : sf, rkp, chi, rkes, ar*t4,
                         (sf > 2) ? es : eb);
  }
  return eth*ar*t4 + (1.0 - eth)*es;
}

// sporder-1012, time2_vet_col_source = blend: S = S*(0) + w (S*(dt) - S*(0)),
// w = x^2/(x0^2 + x^2), x = x_E + x_g the cell's exchange stiffness over the step,
// x_E = chat dt rho kappa_E, x_g = c dt rho kappa_P 4 a T^3/((1 + x_E) rho c_v).  The
// relaxed source of the PREDICTED state (~ t^{n+1}) relaxes it by a further dt: O(dt)
// error where the exchange is resolved (first order in time, vcpred-1010).  blend is the
// state's own source there (O(x^3) from it) and the relaxed source where the exchange is
// stiff (x >> x0, w -> 1: no current-J vet_scatter source, which state amplifies).
template <class EosT>
KOKKOS_INLINE_FUNCTION
Real VetRelaxedSourceBQ(const EosT &eos, const Real d, const Real tb, const Real eb,
                       const Real rkp, const Real rke, const Real eth, const Real ar,
                       const Real cl, const Real ch, const Real dt,
                       const int sf, const Real chi, const Real rkes, const Real x0,
                       const Real qe = 0.0) {
  const Real s1 = VetRelaxedSourceQ(eos, d, tb, eb, rkp, rke, eth, ar, cl, ch, dt, sf,
                                    chi, rkes, qe);
  const Real s0 = VetRelaxedSourceQ(eos, d, tb, eb, rkp, rke, eth, ar, cl, ch, 0.0, sf,
                                    chi, rkes);
  M1EosDirect<EosT> th{eos};
  Real ee, cv;
  th(d, tb, ee, cv);
  const Real xe = ch*dt*rke;
  const Real xg = cl*dt*rkp*ar*4.0*tb*tb*tb/((1.0 + xe)*fmax(d*cv, 1.0e-300));
  const Real x = xe + xg;
  const Real w = (x*x)/(x0*x0 + x*x);
  return s0 + w*(s1 - s0);
}

//! amdef-1019: the 5864f189 relaxed source VERBATIM (VetRelaxedSource / VetRelaxedSourceB
//! without qe), called whenever relaxed_esc is not named, so the relaxed default keeps its
//! exact device code generation (A100 FMA contraction: the qe-generalised kernel changed
//! vet_source = relaxed runs at round-off, hg145 gbit); relaxed_esc calls the ...Q forms
template <class EosT>
KOKKOS_INLINE_FUNCTION
Real VetRelaxedSource(const EosT &eos, const Real d, const Real tb, const Real eb,
                       const Real rkp, const Real rke, const Real eth, const Real ar,
                       const Real cl, const Real ch, const Real dt,
                       const int sf, const Real chi, const Real rkes) {
  const Real ix = 1.0/(1.0 + ch*dt*rke);
  M1EosDirect<EosT> th{eos};
  Real ee, cv;
  th(d, tb, ee, cv);
  Real ts = tb;
  bool ok = true;
  (void) M1ImplTemperatureT(th, d, tb, ee, cl*dt*rkp*ar*ix, cl*dt*rke*eb*ix, ts, ok);
  if (!(ok && ts > 0.0)) {ts = tb;}
  const Real t2 = ts*ts;
  const Real t4 = t2*t2;
  const Real es = (eb + ch*dt*rkp*ar*t4)*ix;
  // vet_scatter: B* of the relaxed gas with (sf = 1, 2; vet_scatter_j = current) the
  // CURRENT J = eb (E^n or the Picard iterate) or (sf = 3, 4; vet_scatter_j = relaxed)
  // the relaxed E*, the consistent end-of-step pair (S - E* = (E* - eb)/(c dt chi))
  if (sf != 0) {
    return VetScatterMix((sf > 2) ? (sf - 2) : sf, rkp, chi, rkes, ar*t4,
                         (sf > 2) ? es : eb);
  }
  return eth*ar*t4 + (1.0 - eth)*es;
}

// sporder-1012, time2_vet_col_source = blend: S = S*(0) + w (S*(dt) - S*(0)),
// w = x^2/(x0^2 + x^2), x = x_E + x_g the cell's exchange stiffness over the step,
// x_E = chat dt rho kappa_E, x_g = c dt rho kappa_P 4 a T^3/((1 + x_E) rho c_v).  The
// relaxed source of the PREDICTED state (~ t^{n+1}) relaxes it by a further dt: O(dt)
// error where the exchange is resolved (first order in time, vcpred-1010).  blend is the
// state's own source there (O(x^3) from it) and the relaxed source where the exchange is
// stiff (x >> x0, w -> 1: no current-J vet_scatter source, which state amplifies).
template <class EosT>
KOKKOS_INLINE_FUNCTION
Real VetRelaxedSourceB(const EosT &eos, const Real d, const Real tb, const Real eb,
                       const Real rkp, const Real rke, const Real eth, const Real ar,
                       const Real cl, const Real ch, const Real dt,
                       const int sf, const Real chi, const Real rkes, const Real x0) {
  const Real s1 = VetRelaxedSource(eos, d, tb, eb, rkp, rke, eth, ar, cl, ch, dt, sf,
                                    chi, rkes);
  const Real s0 = VetRelaxedSource(eos, d, tb, eb, rkp, rke, eth, ar, cl, ch, 0.0, sf,
                                    chi, rkes);
  M1EosDirect<EosT> th{eos};
  Real ee, cv;
  th(d, tb, ee, cv);
  const Real xe = ch*dt*rke;
  const Real xg = cl*dt*rkp*ar*4.0*tb*tb*tb/((1.0 + xe)*fmax(d*cv, 1.0e-300));
  const Real x = xe + xg;
  const Real w = (x*x)/(x0*x0 + x*x);
  return s0 + w*(s1 - s0);
}

//----------------------------------------------------------------------------------------
//! \fn VetEscapeRate
//! \brief vetsrc-1020 (vet_source = relaxed_esc): the photon escape rate per unit chat of
//! a cell of width h (the column axis: dr on sp / cs, dx1 on Cartesian) and extinction
//! chi, qe = 1/(h (2 + chi h/3)); 0 when off (the relaxed source, bitwise)
KOKKOS_INLINE_FUNCTION
Real VetEscapeRate(const bool on, const Real chi, const Real h) {
  return on ? 1.0/(h*(2.0 + chi*h*(1.0/3.0))) : 0.0;
}

} // namespace radm1

#endif // RAD_M1_RAD_M1_VETSRC_HPP_
