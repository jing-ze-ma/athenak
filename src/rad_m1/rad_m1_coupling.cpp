//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_coupling.cpp
//! \brief the implicit local matter coupling of design sect. 4, called inside BOTH
//! PD-ARS stages, and the gas half of the PD-ARS combination.
//!
//! Per cell, with (*) the state left by the explicit transport update of this stage:
//!
//!  (a) ENERGY EXCHANGE, implicit.  With E' = [E* + chat dt (rho kappa_P a T'^4
//!      - rho kappa_E dE0)]/(1 + chat dt rho kappa_E), solve
//!         y(T') = e_gas(rho,T') + (c/chat) E'(T') - [e_gas* + (c/chat) E*] = 0
//!      by a bracketed Newton with bisection fallback.  T'^4 is NOT linearised: that is
//!      what makes a tabulated EOS converge.  The gas energy is then SET algebraically
//!      from the invariant rather than integrated, so e_gas + (c/chat) E is conserved to
//!      round-off and not to the solver tolerance.  dE0 = E0* - E* is the O(beta)
//!      velocity correction of the comoving energy density, taken EXPLICITLY at (*).
//!  (b) FLUX, linear backward Euler:
//!         F' = [F* + chat dt rho k_t (v E' + v.P') - chat c dt beta g0']/(1+chat dt k_t)
//!         delta(rho v) = -(F' - F*)/(chat c),     g0' = (E* - E')/(chat dt)
//!  (c) WORK, explicit: W = vbar.delta(rho v) with vbar the mean of the old and the new
//!      velocity -- exactly the kinetic-energy change -- so u(IEN) += W leaves the gas
//!      INTERNAL energy of (a) untouched, and E' -= (chat/c) W keeps the invariant.
//!
//! PD-ARS (design sect. 5, corrected for the gas).  With L the transport operator and S
//! this solve, stage 1 is U1 = U^n + dt L(U^n) + dt S(U1) and stage 2
//!   U^{n+1} = U^n + dt/2[L(U^n) + L(U1)] + dt/2 S(U1) + dt/2 S(U^{n+1}).
//! The module's Heun weights already put U^n + dt/2[L(U^n)+L(U1)] + dt/2 S(U1) into the
//! radiation array at the start of stage 2, because the stage-1 source increment is
//! carried inside U1 and halved by the average.  The GAS has no transport update here,
//! so the same average has to be formed explicitly: ugas1 holds the gas state at the
//! start of the substep and stage 2 replaces g by (g + g^n)/2 before solving.  Without
//! that the two halves of the invariant get different weights and conservation of
//! e_gas + (c/chat) E breaks at O(dt).

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::Coupling
//! \brief design sect. 4, one cell-local implicit solve per active cell.

TaskStatus RadiationM1::Coupling(Driver *pdrive, int stage) {
  if (!coupling) return TaskStatus::complete;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto u0_ = u0;
  auto opac_ = opac;
  auto ug1 = ugas1;
  auto uh = pmy_pack->phydro->u0;
  auto eos = pmy_pack->phydro->peos->eos_data;
  auto cnt_ = cnt;

  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  Real ar = arad;
  bool edd = eddington;
  bool feedback = gas_feedback;
  // stage 1 integrates the source over dt, stage 2 over dt/2 (the PD-ARS tableau)
  Real dti = (stage == 1) ? dt_sub : (0.5*dt_sub);
  bool stage1 = (stage == 1);

  par_for("m1_coupling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // ---- the gas half of the PD-ARS combination
    Real gi[4];
    if (stage1) {
      ug1(m,0,k,j,i) = uh(m,IEN,k,j,i);
      ug1(m,1,k,j,i) = uh(m,IM1,k,j,i);
      ug1(m,2,k,j,i) = uh(m,IM2,k,j,i);
      ug1(m,3,k,j,i) = uh(m,IM3,k,j,i);
      gi[0] = uh(m,IEN,k,j,i);
      gi[1] = uh(m,IM1,k,j,i);
      gi[2] = uh(m,IM2,k,j,i);
      gi[3] = uh(m,IM3,k,j,i);
    } else {
      gi[0] = 0.5*(uh(m,IEN,k,j,i) + ug1(m,0,k,j,i));
      gi[1] = 0.5*(uh(m,IM1,k,j,i) + ug1(m,1,k,j,i));
      gi[2] = 0.5*(uh(m,IM2,k,j,i) + ug1(m,2,k,j,i));
      gi[3] = 0.5*(uh(m,IM3,k,j,i) + ug1(m,3,k,j,i));
      if (feedback) {
        uh(m,IEN,k,j,i) = gi[0];
        uh(m,IM1,k,j,i) = gi[1];
        uh(m,IM2,k,j,i) = gi[2];
        uh(m,IM3,k,j,i) = gi[3];
      }
    }

    Real rkp = opac_(m,M1_OP_P,k,j,i);
    Real rke = opac_(m,M1_OP_E,k,j,i);
    Real rkt = opac_(m,M1_OP_T,k,j,i);
    if (rkp == 0.0 && rke == 0.0 && rkt == 0.0) return;

    Real dd = uh(m,IDN,k,j,i);
    Real idd = 1.0/fmax(dd, 1.0e-300);
    Real v1 = gi[1]*idd, v2 = gi[2]*idd, v3 = gi[3]*idd;
    Real ekin = 0.5*dd*(v1*v1 + v2*v2 + v3*v3);
    Real eg = gi[0] - ekin;                       // gas INTERNAL energy density

    Real es = fmax(u0_(m,M1_E,k,j,i), efl);
    Real fs1 = u0_(m,M1_F1,k,j,i);
    Real fs2 = u0_(m,M1_F2,k,j,i);
    Real fs3 = u0_(m,M1_F3,k,j,i);

    // closure of the (*) state, used for P* and kept as the flux DIRECTION in (b)
    Real r1, r2, r3, rn;
    M1ReducedFlux(cl, es, fs1, fs2, fs3, r1, r2, r3, rn);
    Real chi = edd ? (1.0/3.0) : M1Chi(rn);
    Real p11, p21, p31, p12, p22, p32, p13, p23, p33;
    M1PressureCol(1, es, r1, r2, r3, rn, chi, p11, p21, p31);
    M1PressureCol(2, es, r1, r2, r3, rn, chi, p12, p22, p32);
    M1PressureCol(3, es, r1, r2, r3, rn, chi, p13, p23, p33);

    Real b1 = v1/cl, b2 = v2/cl, b3 = v3/cl;
    Real b2sq = b1*b1 + b2*b2 + b3*b3;
    Real bdotf = (b1*fs1 + b2*fs2 + b3*fs3)/cl;
    Real bpb = b1*(p11*b1 + p12*b2 + p13*b3) + b2*(p21*b1 + p22*b2 + p23*b3)
             + b3*(p31*b1 + p32*b2 + p33*b3);
    // E0 = (1 + beta^2) E - 2 beta.F/c + beta.P.beta; only the correction is explicit
    Real de0 = b2sq*es - 2.0*bdotf + bpb;

    Real ctc = cl/ch;
    Real inv = 1.0/(1.0 + ch*dti*rke);
    Real tot = eg + ctc*es;                      // the invariant, per unit volume
    Real ep = es;

    // ---------------------------------------------------------------- (a) energy
    if (rkp > 0.0 || rke > 0.0) {
      Real ca = ch*dti*rkp*ar;                   // E'(T) = (es + ca T^4 + cb)*inv
      Real cb = -ch*dti*rke*de0;
      Real tg = eos.Temperature(dd, fmax(eg, 1.0e-300));
      // y(T) = e_gas(d,T) + (c/chat) E'(T) - tot, strictly increasing in T
      Real tlo = tg, thi = tg;
      Real ylo, yhi;
      {
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(dd, tg, ee, pp, cr, ct, cv);
        Real y0 = ee + ctc*fmax((es + ca*tg*tg*tg*tg + cb)*inv, 0.0) - tot;
        ylo = y0;
        yhi = y0;
      }
      int nit = 0;
      bool ok = true;
      if (ylo > 0.0) {
        // T' < T: walk the low end down
        for (int it=0; it<80 && ylo > 0.0; ++it) {
          tlo *= 0.5;
          Real ee, pp, cr, ct, cv;
          eos.ThermoAt(dd, tlo, ee, pp, cr, ct, cv);
          ylo = ee + ctc*fmax((es + ca*tlo*tlo*tlo*tlo + cb)*inv, 0.0) - tot;
          ++nit;
        }
        ok = (ylo <= 0.0);
      } else if (yhi < 0.0) {
        for (int it=0; it<80 && yhi < 0.0; ++it) {
          thi *= 2.0;
          Real ee, pp, cr, ct, cv;
          eos.ThermoAt(dd, thi, ee, pp, cr, ct, cv);
          yhi = ee + ctc*fmax((es + ca*thi*thi*thi*thi + cb)*inv, 0.0) - tot;
          ++nit;
        }
        ok = (yhi >= 0.0);
      }
      Real tp = tg;
      if (!ok) {
        // no bracket: leave the state alone and count it
        Kokkos::atomic_add(&cnt_.d_view(M1_CNT_NBRAK), 1.0);
        Kokkos::atomic_add(&cnt_.d_view(M1_CNT_NFAIL), 1.0);
      } else {
        // Safeguarded Newton (Press et al. `rtsafe`): a Newton step is taken only when
        // it stays inside the bracket AND is shrinking the interval at least as fast as
        // bisection would; otherwise the step IS a bisection.  The plain "Newton, fall
        // back to bisection only when it leaves the bracket" version stalls on a
        // TABULATED EOS, whose c_v is a separate interpolated surface and therefore not
        // exactly the slope of its own e(T): Newton then converges from one side while
        // the far end of the bracket never moves.
        tp = 0.5*(tlo + thi);
        if (tlo == thi) {tp = tlo;}
        Real dxold = fabs(thi - tlo);
        Real dx = dxold;
        int it = 0;
        bool conv = (tlo == thi);
        for (; it<M1_MAXIT && !conv; ++it) {
          Real ee, pp, cr, ct, cv;
          eos.ThermoAt(dd, tp, ee, pp, cr, ct, cv);
          Real t3 = tp*tp*tp;
          Real y = ee + ctc*fmax((es + ca*tp*t3 + cb)*inv, 0.0) - tot;
          if (y > 0.0) {thi = tp;} else {tlo = tp;}
          Real dy = dd*cv + cl*dti*rkp*4.0*ar*t3*inv;
          dxold = dx;
          Real tn = (dy > 0.0) ? (tp - y/dy) : tp;
          if (!(tn > tlo && tn < thi) || (fabs(2.0*y) > fabs(dxold*dy))) {
            tn = 0.5*(tlo + thi);
          }
          dx = fabs(tn - tp);
          conv = (dx <= M1_RTOL*fabs(tn)) ||
                 (fabs(thi - tlo) <= M1_RTOL*fabs(tp));
          tp = tn;
        }
        nit += it;
        if (!conv) {Kokkos::atomic_add(&cnt_.d_view(M1_CNT_NFAIL), 1.0);}
      }
      Kokkos::atomic_add(&cnt_.d_view(M1_CNT_NSOLVE), 1.0);
      Kokkos::atomic_add(&cnt_.d_view(M1_CNT_ITSUM), static_cast<Real>(nit));
      Kokkos::atomic_max(&cnt_.d_view(M1_CNT_ITMAX), static_cast<Real>(nit));

      ep = fmax((es + ca*tp*tp*tp*tp + cb)*inv, efl);
      // the gas energy is SET from the invariant: conservation is algebraic
      eg = tot - ctc*ep;
    }

    // ---------------------------------------------------------------- (b) flux
    Real fp1 = fs1, fp2 = fs2, fp3 = fs3;
    if (rkt > 0.0) {
      Real g0p = (dti > 0.0) ? ((es - ep)/(ch*dti)) : 0.0;
      // P' = P(E', f*): the closure is homogeneous of degree one in E at fixed f
      Real sc = (es > 0.0) ? (ep/es) : 0.0;
      Real a1 = v1*ep + sc*(v1*p11 + v2*p21 + v3*p31);
      Real a2 = v2*ep + sc*(v1*p12 + v2*p22 + v3*p32);
      Real a3 = v3*ep + sc*(v1*p13 + v2*p23 + v3*p33);
      Real den = 1.0/(1.0 + ch*dti*rkt);
      fp1 = (fs1 + ch*dti*rkt*a1 - ch*dti*v1*g0p)*den;
      fp2 = (fs2 + ch*dti*rkt*a2 - ch*dti*v2*g0p)*den;
      fp3 = (fs3 + ch*dti*rkt*a3 - ch*dti*v3*g0p)*den;
    }
    Real dm1 = -(fp1 - fs1)/(ch*cl);
    Real dm2 = -(fp2 - fs2)/(ch*cl);
    Real dm3 = -(fp3 - fs3)/(ch*cl);

    // ---------------------------------------------------------------- (c) work
    Real work = 0.0;
    if (feedback) {
      Real w1 = (gi[1] + dm1)*idd;
      Real w2 = (gi[2] + dm2)*idd;
      Real w3 = (gi[3] + dm3)*idd;
      work = 0.5*((v1 + w1)*dm1 + (v2 + w2)*dm2 + (v3 + w3)*dm3);
      ep -= (ch/cl)*work;
    }

    // limits, then write back
    M1ApplyLimits(cl, efl, ep, fp1, fp2, fp3);
    u0_(m,M1_E,k,j,i) = ep;
    u0_(m,M1_F1,k,j,i) = fp1;
    u0_(m,M1_F2,k,j,i) = fp2;
    u0_(m,M1_F3,k,j,i) = fp3;
    if (feedback) {
      uh(m,IM1,k,j,i) = gi[1] + dm1;
      uh(m,IM2,k,j,i) = gi[2] + dm2;
      uh(m,IM3,k,j,i) = gi[3] + dm3;
      // e_gas(new) + KE(old) + W, and W is exactly KE(new) - KE(old)
      uh(m,IEN,k,j,i) = eg + ekin + work;
    }
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::HydroConToPrim
//! \brief refresh hydro's GHOST ZONES and primitives after the coupling has written its
//! conserved u0, once per substep stage.
//!
//! The ConToPrim alone (which is all design sect. 5 asks for) is not enough, and the
//! reason cost a debugging round in T5: the hydro stage chain is
//! CopyCons -> Fluxes -> ... -> SendU -> RecvU -> BCs -> ConToPrim, so its FLUXES read
//! ghost-zone primitives that were exchanged BEFORE this module last wrote u0.  A
//! uniform single zone whose gas energy the coupling had just dropped by two decades
//! came back with a density spread of 7e-8 and a spurious momentum, because the first
//! hydro flux of the next cycle saw stale ghosts: a jump of the full energy change at
//! the MeshBlock boundary.  So the conserved state is restricted, exchanged and
//! bounded here, exactly as Driver::InitBoundaryValuesAndPrimitives does it, before
//! the inversion.

TaskStatus RadiationM1::HydroConToPrim(Driver *pdrive, int stage) {
  if (!coupling || !gas_feedback) return TaskStatus::complete;
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph == nullptr) return TaskStatus::complete;
  (void) ph->RestrictU(pdrive, 0);
  (void) ph->InitRecv(pdrive, -1);    // stage < 0 suppresses InitFluxRecv
  (void) ph->SendU(pdrive, 0);
  (void) ph->ClearSend(pdrive, -1);   // stage = -1: only clear SendU
  (void) ph->ClearRecv(pdrive, -1);   // stage = -1: only clear RecvU
  (void) ph->RecvU(pdrive, 0);
  (void) ph->ApplyPhysicalBCs(pdrive, 0);
  (void) ph->Prolongate(pdrive, 0);
  return ph->ConToPrim(pdrive, 0);
}

} // namespace radm1
