//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_fluid.cpp
//! \brief the fluid side of the M1 module: the ghost exchange and inversion after a write
//! of the conserved fluid (HydroConToPrim, FluidRefresh), the near-floor gas cap
//! (FloorGasCap) and |B|^2/2 (EmagBuild).
//!
//! THE REFERENCE FORCE, <rad_vet>/force_reference = wb_arad (milestone 2a).  When an
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

#include <cstdlib>
#include <iostream>
#include <math.h>

#include "utils/eint_ledger.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_mergelat.hpp"
#include "mhd/mhd.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "pgen/pgen.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "utils/eint_from_cons.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::HydroConToPrim
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

TaskStatus RadiationVet::HydroConToPrim(Driver *pdrive, int stage) {
  Time2GasTime(pdrive, 2);   // (sdirk2: the gas of the stage-1 solve at t^{n+1})
  // <hydro>/merge_lat_rad = true: the radiation E, F of the merged lateral
  // groups replaced by their group values after the solve (hydro_mergelat.hpp)
  {
    hydro::Hydro *phm = pmy_pack->phydro;
    if (phm != nullptr && phm->pmlat != nullptr && phm->pmlat->on && phm->pmlat->rad) {
      DvceArray4D<Real> nophi;
      phm->pmlat->Average(u0, M1_NVAR, M1_NVAR, 1, nophi, false);
    }
  }
  eint_ledger::Mark(pmy_pack->pmesh, "vet_ravg");
  if (!coupling || !gas_feedback) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  if (pmy_pack->pmesh->pgen->user_m1post_func != nullptr) {
    (pmy_pack->pmesh->pgen->user_m1post_func)(pmy_pack->pmesh);
  }
  eint_ledger::Mark(pmy_pack->pmesh, "vet_post");
  // not on the state of a failed solve, which dt_retry discards (its
  // inf / NaN cells made the cap's energy counters inf for the rest of the run)
  if (fcap_rho > 0.0 && !dtr_fail) {FloorGasCap();}
  FluidRefresh(pdrive);
  eint_ledger::Mark(pmy_pack->pmesh, "vet_c2p");
  TmrMark(9);
  return TaskStatus::complete;
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
  // cell volume for the energy counters (spherical polar; 1 elsewhere)
  const bool hv = pmy_pack->pmesh->use_spherical_polar;
  auto &x1f = pmy_pack->pcoord->xx1f;
  auto &x2f = pmy_pack->pcoord->xx2f;
  auto &x3f = pmy_pack->pcoord->xx3f;
  auto cnt_ = fcap_cnt;
  const Real rc = fcap_rho, vc = fcap_v;
  const Real ar = arad;
  const bool mv = fcap_mvol;
  CleanArgs x;
  if (mv) x = CleanMakeArgs();
  par_for("vet_floorcap", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dd = uh(m,IDN,k,j,i);
    if (!(dd < rc) || !(dd > 0.0)) return;
    const Real wv = mv ? CleanVol(x, m, k, j, i) :
        (hv ? ((x1f(m,i+1)*x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i)*x1f(m,i))
               /3.0*(cos(x2f(m,j)) - cos(x2f(m,j+1)))*(x3f(m,k+1) - x3f(m,k))) : 1.0);
    Real m1 = uh(m,IM1,k,j,i), m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
    const Real msq = m1*m1 + m2*m2 + m3*m3;
    Real ekin = 0.5*msq/dd;
    const Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
    Real etot = uh(m,IEN,k,j,i);
    Real erad = u0_(m,M1_E,k,j,i);
    if (vc > 0.0 && msq > SQR(vc*dd)) {
      const Real fs = vc*dd/sqrt(msq);
      const Real dke = (1.0 - fs*fs)*ekin;
      // (the momentum removed, booked: |dp| V and, spherical polar or with
      // floorcap_mesh_vol on every mesh, its L_z)
      Real p[3] = {m1, m2, m3};
      if (mv) CleanGasMomentum(x, uh, m, k, j, i, p);
      Kokkos::atomic_add(&cnt_.d_view(4), (1.0 - fs)*(mv ? sqrt(p[0]*p[0] + p[1]*p[1] +
                         p[2]*p[2]) : sqrt(msq))*wv);
      if (mv) {
        Kokkos::atomic_add(&cnt_.d_view(5), (1.0 - fs)*CleanLzArm(x, m, k, j, i)*p[2]*wv);
      } else if (hv) {
        Kokkos::atomic_add(&cnt_.d_view(5), (1.0 - fs)*m3*0.5*(x1f(m,i) + x1f(m,i+1))*
                           sin(0.5*(x2f(m,j) + x2f(m,j+1)))*wv);
      }
      m1 *= fs; m2 *= fs; m3 *= fs;
      uh(m,IM1,k,j,i) = m1; uh(m,IM2,k,j,i) = m2; uh(m,IM3,k,j,i) = m3;
      ekin -= dke;
      etot -= dke;
      erad += dke;
      Kokkos::atomic_add(&cnt_.d_view(0), 1.0);
      Kokkos::atomic_add(&cnt_.d_view(2), dke*wv);
    }
    const Real eg = etot - ekin - egrv;
    if (eg > 0.0) {
      const Real tg = eos.Temperature(dd, eg);
      const Real trad = sqrt(sqrt(fmax(erad, 0.0)/ar));
      // never below the EOS temperature floor: the C2P that follows would raise e_int
      // back to e(rho, tfloor) and create the energy this cap moved to E_rad
      const Real tcap = fmax(trad, eos.tfloor);
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
      std::cout << "<rad_vet> floorcap cycle " << cyc << " (cumulative, all calls): n_v "
                << c[0] << " n_t " << c[1] << " dKE " << c[2] << " dEint " << c[3]
                << " (erg, to E_rad)"
                << " momentum removed |dp| " << c[4] << " dL_z " << c[5] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FluidRefresh
//! \brief restrict, exchange, bound and invert the fluid after a host-side write of its
//! conserved u0 (see HydroConToPrim above for why).  Hydro: the hydro task calls
//! verbatim (bitwise).  MHD (m1-mhd): M1 never writes b0, so only u0 is exchanged, and
//! through pbval_u directly: MHD::InitRecv(-1) also posts the B receives and
//! MHD::ClearRecv(-1) MPI_Waits on them, and no B is sent here -- a deadlock on more
//! than one rank.  MHD::ApplyPhysicalBCs re-applies the B BCs to an unchanged b0
//! (harmless) and the pgen's user BC, Prolongate and ConToPrim rebuild bcc0 and w0.
//! lat = false: without the merged lateral groups' average (Time2GasTime's flux state).

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
  // the gas the solve wrote, back to group values (hydro_mergelat.hpp)
  if (lat && ph->pmlat != nullptr && ph->pmlat->on) {
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
//! \fn void RadiationVet::EmagBuild
//! \brief m1-mhd: |B|^2/2 per cell (incl. ghosts) of b0 (un = false -> emag0) or of the
//! U^n copy b1 (un = true -> emag1), with the MHD C2P's face average (M1EmagCell).
//! No-op for hydro.

void RadiationVet::EmagBuild(bool un) {
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
    par_for("vet_emag_cs", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
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
  par_for("vet_emag", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    em(m,k,j,i) = M1EmagCell(b1f, b2f, b3f, sph, csr, x1v, x1f, x2v, x2f, x3v, x3f,
                             m, k, j, i);
  });
}


//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::Opacity

TaskStatus RadiationVet::Opacity(Driver *pdrive, int stage) {
  // first gas reader of both M1 chains: refresh |B|^2/2 of the stage's field
  if (fl_mhd && !emag_hold) EmagBuild(false);
  if (opac_zero) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  if (opacity_type == M1_OPAC_TABLE && otab.nT <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/opacity = table, but no tables were handed to the "
      << "module: the problem generator must call RadiationVet::SetOpacityTables"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto opac_ = opac;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  // with <hydro>/etotgrav the conserved energy includes rho*Phi: subtract it as well
  const bool etg = fl.use_etotgrav;
  auto phicc = fl.phicc0;
  const bool mhd = fl_mhd;      // m1-mhd: and |B|^2/2
  auto emag_ = emag0;
  // cubed sphere: the momentum is covariant, KE = 0.5 m_a v^a (gnomonic_raisevel.hpp)
  const bool csk = pmy_pack->pmesh->use_cubed_sphere;
  auto ccl = pmy_pack->pcoord->cos_cell;
  int otype = opacity_type;
  Real kp = kappa_p, ke = kappa_e, kf = kappa_f, ks = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  bool need_t = (otype != M1_OPAC_CONST);
  // copy the tables (held by value) into the closure instead of capturing this
  M1OpacTab ot = otab;
  const bool osc = vscat;   // vet_scatter: rho kappa_e -> opac(M1_OP_S)
  if (ot.scat == 2 && !eos.tbl.active) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/vet_scatter_kappa_e = eos needs a general (tabulated) "
      << "EOS: the electron fraction comes from its Saha composition" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  par_for_lb("vet_opacity", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real d = uh(m,IDN,k,j,i);
    Real t = 0.0;
    if (need_t) {
      Real ke_dens = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                          SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        ke_dens = CsKeMom(uh(m,IM1,k,j,i), m2, m3, c)
                  /fmax(d, 1.0e-300);
      }
      Real eint = uh(m,IEN,k,j,i) - ke_dens;
      if (etg) eint -= d*phicc(m,k,j,i);
      if (mhd) eint -= emag_(m,k,j,i);
      t = eos.Temperature(d, fmax(eint, 0.0));
    }
    Real op, oe, of, os;
    M1CellOpacities(otype, ot, eos, d, t, kp, ke, kf, ks, rref, tref, aa, bb, op, oe, of,
                    os);
    opac_(m,M1_OP_P,k,j,i) = d*op;
    opac_(m,M1_OP_E,k,j,i) = d*oe;
    opac_(m,M1_OP_T,k,j,i) = d*(of + os);
    if (osc) opac_(m,M1_OP_S,k,j,i) = d*os;
  });

  if (!rsla_done) RSLACheck();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RSLACheck
//! \brief reduced-speed-of-light admissibility check, run once on the first filled
//! opacity array: prints K = chat/(v_max tau_max), with tau_max the largest x1 column
//! optical depth over the active cells and v_max = <rad_vet>/rsla_vmax if positive, else
//! the measured max |v|

void RadiationVet::RSLACheck() {
  rsla_done = true;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto opac_ = opac;
  auto uh = FluidRef::Get(pmy_pack).u0;
  auto &size = pmy_pack->pmb->mb_size;

  Real taumax = 0.0;
  Kokkos::parallel_reduce("vet_rsla_tau",
  Kokkos::MDRangePolicy<Kokkos::Rank<3>>(DevExeSpace(), {0,ks,js}, {nmb1+1,ke+1,je+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, Real &lmax) {
    Real tau = 0.0;
    for (int i=is; i<=ie; ++i) {
      tau += opac_(m,M1_OP_T,k,j,i)*size.d_view(m).dx1;
    }
    lmax = (tau > lmax) ? tau : lmax;
  }, Kokkos::Max<Real>(taumax));

  Real vmax = 0.0;
  Kokkos::parallel_reduce("vet_rsla_v",
  Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                         {nmb1+1,ke+1,je+1,ie+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
    Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
    Real v = sqrt(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                  SQR(uh(m,IM3,k,j,i)))/d;
    lmax = (v > lmax) ? v : lmax;
  }, Kokkos::Max<Real>(vmax));

#if MPI_PARALLEL_ENABLED
  {
    Real loc[2] = {taumax, vmax}, glob[2];
    MPI_Allreduce(loc, glob, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    taumax = glob[0];
    vmax = glob[1];
  }
#endif
  Real vmeas = vmax;
  if (rsla_vmax > 0.0) vmax = rsla_vmax;
  // chat = c: nothing to check
  if (chat_over_c >= 1.0) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> RSLA check: chat = c, no reduced speed of light (tau_max="
                << taumax << ", measured v_max=" << vmeas << ")" << std::endl;
    }
    return;
  }
  Real denom = vmax*((taumax > 1.0) ? taumax : 1.0);
  Real ratio = (chat > 0.0) ? (denom/chat) : 0.0;
  Real kk = (denom > 0.0) ? (chat/denom) : 0.0;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> RSLA check: tau_max=" << taumax
              << " v_max=" << vmax << " (measured " << vmeas << ")"
              << " chat=" << chat << " -> K = chat/(v_max tau_max) = " << kk
              << ", v_max tau_max/chat = " << ratio << std::endl;
  }
  if (ratio > 1.0 && !rsla_force) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> reduced speed of light is inadmissible: "
      << "v_max*tau_max/chat = " << ratio << " > 1 (design sect. 2).  Raise "
      << "chat_over_c, or set <rad_vet>/rsla_force = true to run anyway" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ratio > rsla_warn && global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_vet> v_max*tau_max/chat = " << ratio
              << " exceeds <rad_vet>/rsla_warn = " << rsla_warn
              << "; the reduced speed of light is marginal (design sect. 2)"
              << std::endl;
  }
  return;
}

} // namespace radvet
