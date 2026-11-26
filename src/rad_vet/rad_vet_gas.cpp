//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_gas.cpp
//! \brief fluid side of the VET module: ghost exchange and inversion after a write of
//! the conserved fluid, near-floor gas cap, |B|^2/2, the opacity task and RSLA check.
//!
//!  called from: rad_vet_imex.cpp, rad_vet_steady.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_imex.cpp, rad_vet_init.cpp, rad_vet_ledger.cpp

#include <math.h>

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/gnomonic_raisevel.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_hyd.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_mergelat.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "utils/eint_ledger.hpp"
#include "utils/scratch_view.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::GasUpdate
//! \brief refresh hydro ghost zones and primitives after the coupling wrote u0, once
//! per stage: the next hydro fluxes read ghost primitives, so ConToPrim alone would
//! leave stale ghosts (restrict, exchange, BCs as InitBoundaryValuesAndPrimitives).

TaskStatus RadiationVet::GasUpdate(Driver *pdrive, int stage) {
  Time2GasTime(pdrive, 2);   // sdirk2: the gas of the stage-1 solve at t^{n+1}
  // <hydro>/merge_lat_rad: radiation E, F of merged lateral groups -> group values
  {
    hydro::Hydro *phm = pmy_pack->phydro;
    if (phm != nullptr && phm->pmlat != nullptr && phm->pmlat->on && phm->pmlat->rad) {
      DvceArray4D<Real> nophi;
      phm->pmlat->Average(u0, VET_NVAR, VET_NVAR, 1, nophi, false);
    }
  }
  eint_ledger::Mark(pmy_pack->pmesh, "vet_ravg");
  if (!coupling || !gas_feedback) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  if (pmy_pack->pmesh->pgen->user_vetpost_func != nullptr) {
    (pmy_pack->pmesh->pgen->user_vetpost_func)(pmy_pack->pmesh);
  }
  eint_ledger::Mark(pmy_pack->pmesh, "vet_post");
  // skip on a failed solve (dt_retry discards it; inf/NaN would poison the counters)
  if (fcap_rho > 0.0 && !dtr_fail) {FloorGasCap();}
  eint_ledger::Mark(pmy_pack->pmesh, "floorcap");
  FluidRefresh(pdrive);
  eint_ledger::Mark(pmy_pack->pmesh, "vet_c2p");
  TmrMark(9);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void FloorCapMomentum
//! \brief the floor-gas cap's |m|^2 msq of a cell momentum p (the cubed sphere: m2, m3
//! covariant on the unit tangents, cos_cell c, so msq = m1^2 + (m2^2 - 2 c m2 m3 +
//! m3^2)/(1 - c^2), the metric norm of the ceiling's |v|) and, with lz, its phi_hat f
//! (spherical polar e3, Cartesian (-y, x)/R of the cell, the cubed sphere its covariant
//! (0, f2, f3), GnomonicPhiHat), the part pf = m.phi_hat that carries L_z (the cubed
//! sphere m_a phi^a) and the Cartesian arm R of the cell (else 0)

KOKKOS_INLINE_FUNCTION
void FloorCapMomentum(const bool sph, const bool cub, const bool lz, const int pnl,
                      const RegionSize &sz, const Real c, const int ii, const int jj,
                      const int kk, const int nx1, const int nx2, const int nx3,
                      const Real p[3], Real f[3], Real &msq, Real &pf, Real &r) {
  const Real det = 1.0 - c*c;
  msq = cub ? p[0]*p[0] + (p[1]*p[1] - 2.0*c*p[1]*p[2] + p[2]*p[2])/det
            : p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  f[0] = f[1] = 0.0;
  f[2] = sph ? 1.0 : 0.0;
  pf = r = 0.0;
  if (!lz) return;
  if (cub) {
    GnomonicPhiHat(pnl, 0.25*M_PI*CellCenterX(jj, nx2, sz.x2min, sz.x2max),
                   0.25*M_PI*CellCenterX(kk, nx3, sz.x3min, sz.x3max), f[1], f[2]);
  } else if (!sph) {
    const Real cx = CellCenterX(ii, nx1, sz.x1min, sz.x1max);
    const Real cy = CellCenterX(jj, nx2, sz.x2min, sz.x2max);
    r = sqrt(cx*cx + cy*cy);
    if (r > 0.0) {
      f[0] = -cy/r;
      f[1] = cx/r;
    }
  }
  pf = cub ? (p[1]*(f[1] - c*f[2]) + p[2]*(f[2] - c*f[1]))/det
           : p[0]*f[0] + p[1]*f[1] + p[2]*f[2];
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FloorGasCap
//! \brief cap |v| and T_gas of near-floor gas (rho < fcap_rho) after
//! the implicit solve, active cells, on the conserved state (FluidRefresh follows).
//! The cap acts where the gas carries a negligible share of the energy (rho 4+ decades
//! below the photosphere) and only there sets the CFL dt (hot, fast floor pockets).
//! T <= T_rad mimics Compton equilibration, which the model has no term for.  Energy:
//! the clipped kinetic energy and the excess internal energy go to E_rad,
//! so gas + radiation energy is conserved per cell (chat = c).  Momentum is removed.
//! floorcap_lz_spin: the whole velocity is capped (so |v| <= floorcap_v), and the L_z
//! X it removes goes to the gas at or above floorcap_rho as one rigid rotation
//! dOmega = X/I, I = sum rho V R^2 (both summed per cell, then in a fixed order:
//! SumCells), its kinetic energy paid from the cell's E_rad: L_z and E exact.

void RadiationVet::FloorGasCap() {
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
  auto cnt_ = fcap_cnt;
  const Real rc = fcap_rho, vc = fcap_v;
  const Real ar = arad;
  const bool flz = fcap_lz, sph = pmy_pack->pmesh->use_spherical_polar;
  auto mbs = pmy_pack->pmb->mb_size;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  // the cubed sphere: the metric norm of the covariant m2, m3 (FloorCapMomentum)
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  auto cosc = cub ? pmy_pack->pcoord->cos_cell : DvceArray3D<Real>("fcap_c",1,1,1);
  auto pnl = pmy_pack->pmb->mb_panel.d_view;
  // the counters by the mesh cell volume, |dp| and L_z in the frame of y (CellVolume,
  // LzArm)
  SolveArgs x = MakeArgs();
  // (floorcap_lz_spin: dOmega = X/I, from the state before the cap)
  const bool spn = fcap_spin;
  Real dom = 0.0;
  if (spn) {
    auto lxx = ScratchView<DvceArray4D<Real>>("fcap_lxx", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    auto lxi = ScratchView<DvceArray4D<Real>>("fcap_lxi", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    par_for("vet_floorcap_lz", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dd = uh(m,IDN,k,j,i);
      const Real p[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
      Real f[3], msq, pf, r;
      FloorCapMomentum(sph, cub, true, pnl(m), mbs.d_view(m), cub ? cosc(m,k,j) : 0.0,
                       i - is, j - js, k - ks, nx1, nx2, nx3, p, f, msq, pf, r);
      if (sph || cub) r = LzArm(x, m, k, j, i);
      const Real wv = CellVolume(x, m, k, j, i);
      const bool cap = (dd < rc) && (dd > 0.0) && (vc > 0.0) && (msq > SQR(vc*dd));
      lxx(m,k,j,i) = cap ? (1.0 - vc*dd/sqrt(msq))*pf*r*wv : 0.0;
      lxi(m,k,j,i) = (dd >= rc) ? dd*wv*r*r : 0.0;
    });
    const Real lx = SumCells(lxx), li = SumCells(lxi);
    if (lx != 0.0 && !(li > 0.0)) {
      std::cout << "### FATAL ERROR: <rad_vet>/floorcap_lz_spin: the cap removed L_z "
                << lx << " but no gas is at or above floorcap_rho to take it (name "
                << "floorcap_lz_spin = false)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (lx != 0.0) dom = lx/li;
  }
  auto sl = ScratchView<DvceArray1D<Real>>("cln_csl", 6*CLN_NSL);
  Kokkos::deep_copy(sl, 0.0);
  par_for("vet_floorcap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dd = uh(m,IDN,k,j,i);
    if (!(dd > 0.0) || (!(dd < rc) && dom == 0.0)) return;
    const Real wv = CellVolume(x, m, k, j, i);
    Real m1 = uh(m,IM1,k,j,i), m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
    const Real c = cub ? cosc(m,k,j) : 0.0;
    const Real pm[3] = {m1, m2, m3};
    Real f[3], msq, pf, rf;
    FloorCapMomentum(sph, cub, flz, pnl(m), mbs.d_view(m), c, i - is, j - js, k - ks, nx1,
                     nx2, nx3, pm, f, msq, pf, rf);
    if (!(dd < rc)) {   // (floorcap_lz_spin: dOmega rho R phi_hat, E_rad pays its KE)
      const Real r = (sph || cub) ? LzArm(x, m, k, j, i) : rf, w = dom*dd*r;
      const Real q[3] = {m1 + w*f[0], m2 + w*f[1], m3 + w*f[2]};
      Real g[3], qsq, qf, qr;
      FloorCapMomentum(sph, cub, false, 0, mbs.d_view(m), c, 0, 0, 0, 1, 1, 1, q, g, qsq,
                       qf, qr);
      const Real dk = 0.5*(qsq - msq)/dd;
      uh(m,IM1,k,j,i) = q[0]; uh(m,IM2,k,j,i) = q[1]; uh(m,IM3,k,j,i) = q[2];
      uh(m,IEN,k,j,i) += dk;
      u0_(m,VET_E,k,j,i) -= dk;
      SlotAdd(sl, 5, m, k, j, i, -w*r*wv);
      return;
    }
    Real ekin = 0.5*msq/dd;
    const Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
    Real etot = uh(m,IEN,k,j,i);
    Real erad = u0_(m,VET_E,k,j,i);
    if (vc > 0.0 && msq > SQR(vc*dd) && flz && !spn) {
      // (floorcap_lz: the part pf = m.phi_hat that carries L_z never scaled, EkinScale)
      Real fs;
      const Real ek = EkinScale(ekin, 0.5*pf*pf/dd, 0.5*vc*vc*dd, fs);
      const Real dke = ekin - ek;
      SlotAdd(sl, 4, m, k, j, i, (1.0 - fs)*sqrt(fmax(msq - pf*pf, 0.0))*wv);
      m1 = fs*(m1 - pf*f[0]) + pf*f[0];
      m2 = fs*(m2 - pf*f[1]) + pf*f[1];
      m3 = fs*(m3 - pf*f[2]) + pf*f[2];
      uh(m,IM1,k,j,i) = m1; uh(m,IM2,k,j,i) = m2; uh(m,IM3,k,j,i) = m3;
      ekin = ek;
      etot -= dke;
      erad += dke;
      SlotAdd(sl, 0, m, k, j, i, 1.0);
      SlotAdd(sl, 2, m, k, j, i, dke*wv);
    } else if (vc > 0.0 && msq > SQR(vc*dd)) {
      const Real fs = vc*dd/sqrt(msq);
      const Real dke = (1.0 - fs*fs)*ekin;
      // book the removed momentum |dp| V and its L_z
      Real p[3] = {m1, m2, m3};
      GasMomentum(x, uh, m, k, j, i, p);
      SlotAdd(sl, 4, m, k, j, i, (1.0 - fs)*sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])*wv);
      SlotAdd(sl, 5, m, k, j, i, (1.0 - fs)*LzArm(x, m, k, j, i)*p[2]*wv);
      m1 *= fs; m2 *= fs; m3 *= fs;
      uh(m,IM1,k,j,i) = m1; uh(m,IM2,k,j,i) = m2; uh(m,IM3,k,j,i) = m3;
      ekin -= dke;
      etot -= dke;
      erad += dke;
      SlotAdd(sl, 0, m, k, j, i, 1.0);
      SlotAdd(sl, 2, m, k, j, i, dke*wv);
    }
    const Real eg = etot - ekin - egrv;
    if (eg > 0.0) {
      const Real tg = eos.Temperature(dd, eg);
      const Real trad = sqrt(sqrt(fmax(erad, 0.0)/ar));
      // not below the EOS floor: the C2P would raise e_int back and create energy
      const Real tcap = fmax(trad, eos.tfloor);
      if (tg > tcap) {
        const Real enew = eos.EnergyFromTemperature(dd, tcap);
        const Real de = eg - enew;
        if (de > 0.0) {
          etot -= de;
          erad += de;
          SlotAdd(sl, 1, m, k, j, i, 1.0);
          SlotAdd(sl, 3, m, k, j, i, de*wv);
        }
      }
    }
    uh(m,IEN,k,j,i) = etot;
    u0_(m,VET_E,k,j,i) = erad;
  });
  SlotSum(sl, cnt_.d_view, 6);
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
      std::cout << "<rad_vet> floorcap cycle " << cyc << " (cumulative, all calls): n_v "
                << c[0] << " n_t " << c[1] << " dKE " << c[2] << " dEint " << c[3]
                << " (erg, to E_rad)"
                << " momentum removed |dp| " << c[4] << " dL_z " << c[5] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluidRefresh
//! \brief restrict, exchange, bound and invert the fluid after a write of u0.
//! Hydro: the hydro task calls. MHD: only u0 is exchanged, through pbval_u directly
//! (MHD::InitRecv/ClearRecv(-1) would wait on B receives never sent: deadlock).
//! lat = false: skip the merged lateral groups' average.

void RadiationVet::FluidRefresh(Driver *pdrive, bool lat) {
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
  // the gas the solve wrote, back to group values
  if (lat && ph->pmlat != nullptr && ph->pmlat->on) {
    ph->pmlat->Average(ph->u0, ph->nhydro, ph->nhydro + ph->nscalars, 0, ph->phicc0,
                       ph->use_etotgrav);
    eint_ledger::Mark(pmy_pack->pmesh, "mrglat");
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

}  // namespace radvet
