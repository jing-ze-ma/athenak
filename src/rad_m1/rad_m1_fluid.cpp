//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_fluid.cpp
//! \brief the fluid side of the M1 module: the ghost exchange and inversion after a write
//! of the conserved fluid (HydroConToPrim, FluidRefresh), the near-floor gas cap
//! (FloorGasCap) and |B|^2/2 (EmagBuild).
//!
//! THE REFERENCE FORCE, <rad_m1>/force_reference = wb_arad (milestone 2a).  When an
//! external well-balanced scheme already carries part of the radiative support -- the
//! effective potential Phi_eff of box_convection, whose gravity source delivers
//! -rho*g + rho*arad_ref to the gas momentum and NO energy term, because a_rad is not in
//! the etotgrav potential -- the coupling must not hand the gas that part a second time.
//! The bookkeeping that makes the split exact is: MOMENTUM residual, WORK full.
//!
//!   radiation field:  F' - F*                 (the full force; unchanged)
//!   gas momentum:     dm_res = dm - dt rho arad_ref,   dm = -(F' - F*)/(chat c)
//!   gas energy:       u(IEN) += W with W = vbar.dm, the work of the FULL force
//!   radiation energy: E' -= (chat/c) W
//!
//! Momentum is conserved because the WB source supplies exactly the dt rho arad_ref that
//! was taken out here.  Energy: writing u(IEN) = e_gas + KE_old + W while the momentum
//! actually applied is dm_res raises the gas INTERNAL energy by vbar.(dm - dm_res) =
//! vbar.dt rho arad_ref -- precisely the amount the WB momentum source will then remove
//! from it when it adds dt rho arad_ref to IM1 without an energy term.  Over the cycle
//! the internal energy is left alone and the radiation field has paid for the whole of
//! the kinetic work, exactly once.  The two vbar are a stage apart, so the cancellation
//! is exact to O(dt^2), the order of the operator split itself.
//!

#include <math.h>

#include <iostream>

#include "utils/eint_ledger.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_mergelat.hpp"
#include "mhd/mhd.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "pgen/pgen.hpp"

namespace radm1 {

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
  // <hydro>/merge_lat_rad = true: the radiation E, F of the merged lateral
  // groups replaced by their group values after the solve (hydro_mergelat.hpp)
  {
    hydro::Hydro *phm = pmy_pack->phydro;
    if (phm != nullptr && phm->pmlat != nullptr && phm->pmlat->on && phm->pmlat->rad) {
      DvceArray4D<Real> nophi;
      phm->pmlat->Average(u0, M1_NVAR, M1_NVAR, 1, nophi, false);
    }
  }
  // <rad_m1>/polar_ring_avg: E, F and the gas the coupling wrote to their
  // phi-group averages (rad_m1_ringavg.cpp), before the fluid refresh below
  if (rav_on) {
    RavState(coupling && gas_feedback && fl_on && pmy_pack->phydro != nullptr);
  }
  eint_ledger::Mark(pmy_pack->pmesh, "m1_ravg");
  if (!coupling || !gas_feedback) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  if (pmy_pack->pmesh->pgen->user_m1post_func != nullptr) {
    (pmy_pack->pmesh->pgen->user_m1post_func)(pmy_pack->pmesh);
  }
  eint_ledger::Mark(pmy_pack->pmesh, "m1_post");
  // not on the state of a failed solve, which dt_retry discards (its
  // inf / NaN cells made the cap's energy counters inf for the rest of the run)
  if (fcap_rho > 0.0 && !dtr_fail) {FloorGasCap();}
  FluidRefresh(pdrive);
  eint_ledger::Mark(pmy_pack->pmesh, "m1_c2p");
  TmrMark(9);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::FloorGasCap
//! \brief cap |v| and T_gas of near-floor gas (rho < fcap_rho) after
//! the implicit solve, active cells, on the conserved state (FluidRefresh follows).
//! The cap acts where the gas carries a negligible share of the energy (rho 4+ decades
//! below the photosphere) and only there sets the CFL dt (hot, fast floor pockets).
//! T <= T_rad mimics Compton equilibration, which the model has no term for.  Energy:
//! the clipped kinetic energy (fcap_ke_rad) and the excess internal energy go to E_rad,
//! so gas + radiation energy is conserved per cell (chat = c).  Momentum is removed.

void RadiationM1::FloorGasCap() {
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph == nullptr) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u0_ = u0;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  const bool etg = fl.use_etotgrav;
  auto phicc = fl.phicc0;
  // cell volume for the energy counters (spherical polar; 1 elsewhere)
  const bool hv = pmy_pack->pmesh->use_spherical_polar;
  auto &x1f = pmy_pack->pcoord->xx1f;
  auto &x2f = pmy_pack->pcoord->xx2f;
  auto &x3f = pmy_pack->pcoord->xx3f;
  auto cnt_ = fcap_cnt;
  const Real rc = fcap_rho, vc = fcap_v;
  const bool tc = fcap_t, kr = fcap_ke_rad;
  const Real tmin = fcap_tmin/otab.tunit;   // code temperature
  const Real ar = arad;
  par_for("m1_floorcap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dd = uh(m,IDN,k,j,i);
    if (!(dd < rc) || !(dd > 0.0)) return;
    const Real wv = hv ? ((x1f(m,i+1)*x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i)*x1f(m,i))
                          /3.0*(cos(x2f(m,j)) - cos(x2f(m,j+1)))*(x3f(m,k+1) - x3f(m,k)))
                       : 1.0;
    Real m1 = uh(m,IM1,k,j,i), m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
    const Real msq = m1*m1 + m2*m2 + m3*m3;
    Real ekin = 0.5*msq/dd;
    const Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
    Real etot = uh(m,IEN,k,j,i);
    Real erad = u0_(m,M1_E,k,j,i);
    if (vc > 0.0 && msq > SQR(vc*dd)) {
      const Real fs = vc*dd/sqrt(msq);
      const Real dke = (1.0 - fs*fs)*ekin;
      // (the momentum removed, booked: |dp| V and, spherical polar, its L_z)
      Kokkos::atomic_add(&cnt_.d_view(4), (1.0 - fs)*sqrt(msq)*wv);
      if (hv) {
        Kokkos::atomic_add(&cnt_.d_view(5), (1.0 - fs)*m3*0.5*(x1f(m,i) + x1f(m,i+1))*
                           sin(0.5*(x2f(m,j) + x2f(m,j+1)))*wv);
      }
      m1 *= fs; m2 *= fs; m3 *= fs;
      uh(m,IM1,k,j,i) = m1; uh(m,IM2,k,j,i) = m2; uh(m,IM3,k,j,i) = m3;
      ekin -= dke;
      if (kr) {
        etot -= dke;
        erad += dke;
      }
      Kokkos::atomic_add(&cnt_.d_view(0), 1.0);
      Kokkos::atomic_add(&cnt_.d_view(2), dke*wv);
    }
    if (tc) {
      const Real eg = etot - ekin - egrv;
      if (eg > 0.0) {
        const Real tg = eos.Temperature(dd, eg);
        const Real trad = sqrt(sqrt(fmax(erad, 0.0)/ar));
        // never below the EOS temperature floor: the C2P that follows would raise e_int
        // back to e(rho, tfloor) and create the energy this cap moved to E_rad
        const Real tcap = fmax(fmax(trad, tmin), eos.tfloor);
        if (tg > tcap) {
          const Real enew = eos.EnergyFromTemperature(dd, tcap);
          const Real de = eg - enew;
          if (de > 0.0) {
            etot -= de;
            erad += de;
            Kokkos::atomic_add(&cnt_.d_view(1), 1.0);
            Kokkos::atomic_add(&cnt_.d_view(3), de*wv);
          }
        }
      }
    }
    uh(m,IEN,k,j,i) = etot;
    u0_(m,M1_E,k,j,i) = erad;
  });
  const int rep = fcap_report;
  const int cyc = pmy_pack->pmesh->ncycle;
  if (rep > 0 && (cyc % rep) == 0) {
    fcap_cnt.modify_device();
    fcap_cnt.sync_host();
    Real c[6];
    for (int n = 0; n < 6; ++n) c[n] = fcap_cnt.h_view(n);
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, c, 6, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> floorcap cycle " << cyc << " (cumulative, all calls): n_v "
                << c[0] << " n_t " << c[1] << " dKE " << c[2] << " dEint " << c[3]
                << " (erg, to E_rad" << (kr ? "" : "; KE thermalised") << ")"
                << " momentum removed |dp| " << c[4] << " dL_z " << c[5] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::FluidRefresh
//! \brief restrict, exchange, bound and invert the fluid after a host-side write of its
//! conserved u0 (see HydroConToPrim above for why).  Hydro: the hydro task calls
//! verbatim (bitwise).  MHD (m1-mhd): M1 never writes b0, so only u0 is exchanged, and
//! through pbval_u directly: MHD::InitRecv(-1) also posts the B receives and
//! MHD::ClearRecv(-1) MPI_Waits on them, and no B is sent here -- a deadlock on more
//! than one rank.  MHD::ApplyPhysicalBCs re-applies the B BCs to an unchanged b0
//! (harmless) and the pgen's user BC, Prolongate and ConToPrim rebuild bcc0 and w0.

void RadiationM1::FluidRefresh(Driver *pdrive) {
  if (pmy_pack->pmhd != nullptr) {
    mhd::MHD *pm = pmy_pack->pmhd;
    (void) pm->RestrictU(pdrive, 0);
    (void) pm->pbval_u->InitRecv(pm->nmhd + pm->nscalars);
    (void) pm->pbval_u->PackAndSendCC(pm->u0, pm->coarse_u0);
    (void) pm->pbval_u->ClearSend();
    (void) pm->pbval_u->ClearRecv();
    (void) pm->pbval_u->RecvAndUnpackCC(pm->u0, pm->coarse_u0);
    (void) pm->ApplyPhysicalBCs(pdrive, 0);
    (void) pm->Prolongate(pdrive, 0);
    (void) pm->ConToPrim(pdrive, 0);
    return;
  }
  hydro::Hydro *ph = pmy_pack->phydro;
  if (ph == nullptr) return;
  // the gas the solve wrote, back to group values (hydro_mergelat.hpp)
  if (ph->pmlat != nullptr && ph->pmlat->on) {
    ph->pmlat->Average(ph->u0, ph->nhydro, ph->nhydro + ph->nscalars, 0, ph->phicc0,
                       ph->use_etotgrav);
  }
  (void) ph->RestrictU(pdrive, 0);
  (void) ph->InitRecv(pdrive, -1);    // stage < 0 suppresses InitFluxRecv
  (void) ph->SendU(pdrive, 0);
  (void) ph->ClearSend(pdrive, -1);   // stage = -1: only clear SendU
  (void) ph->ClearRecv(pdrive, -1);   // stage = -1: only clear RecvU
  (void) ph->RecvU(pdrive, 0);
  (void) ph->ApplyPhysicalBCs(pdrive, 0);
  (void) ph->Prolongate(pdrive, 0);
  (void) ph->ConToPrim(pdrive, 0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::EmagBuild
//! \brief m1-mhd: |B|^2/2 per cell (incl. ghosts) of b0 (un = false -> emag0) or of the
//! U^n copy b1 (un = true -> emag1), with the MHD C2P's face average (M1EmagCell).
//! No-op for hydro.

void RadiationM1::EmagBuild(bool un) {
  mhd::MHD *pm = pmy_pack->pmhd;
  if (pm == nullptr) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1 + 2*(indcs.ng);
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto b1f = un ? pm->b1.x1f : pm->b0.x1f;
  auto b2f = un ? pm->b1.x2f : pm->b0.x2f;
  auto b3f = un ? pm->b1.x3f : pm->b0.x3f;
  auto em = un ? emag1 : emag0;
  const bool sph = pmy_pack->pmesh->use_spherical_polar;
  const bool csr = pmy_pack->pmesh->use_cubed_sphere;
  auto &x1v = pmy_pack->pcoord->x1v;
  auto &x1f = pmy_pack->pcoord->xx1f;
  auto &x2v = pmy_pack->pcoord->x2v;
  auto &x2f = pmy_pack->pcoord->xx2f;
  auto &x3v = pmy_pack->pcoord->x3v;
  auto &x3f = pmy_pack->pcoord->xx3f;
  if (csr) {
    // the cubed sphere's |B|^2 is that of the ORTHONORMAL triple the MHD
    // C2P builds (Coordinates::GnomonicEquiangleRaiseVelMHD): the faces hold the
    // projections B.nhat (nhat_xi.nhat_eta = -c), bx = B.rhat (CellCenteredRadialFld),
    // by = B.e_xi = (b_xi + c b_eta)/s, bz = B.(e_eta - c e_xi)/s = b_eta; the plain
    // sum of the face-normal means (M1EmagCell's generic branch) is wrong by O(c)
    auto ccl = pmy_pack->pcoord->cos_cell;
    auto csn = pmy_pack->pcoord->sin_cell;
    par_for("m1_emag_cs", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real c = ccl(m,k,j), sn = csn(m,k,j);
      const Real bx = CellCenteredRadialFld(b1f(m,k,j,i), b1f(m,k,j,i+1), x1f(m,i),
                                            x1f(m,i+1), x1v(m,i));
      const Real byn = 0.5*(b2f(m,k,j,i) + b2f(m,k,j+1,i));
      const Real bz = 0.5*(b3f(m,k,j,i) + b3f(m,k+1,j,i));
      const Real by = (byn + c*bz)/sn;
      em(m,k,j,i) = 0.5*(bx*bx + by*by + bz*bz);
    });
    return;
  }
  par_for("m1_emag", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    em(m,k,j,i) = M1EmagCell(b1f, b2f, b3f, sph, csr, x1v, x1f, x2v, x2f, x3v, x3f,
                             m, k, j, i);
  });
}

} // namespace radm1
