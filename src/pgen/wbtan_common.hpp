#ifndef PGEN_WBTAN_COMMON_HPP_
#define PGEN_WBTAN_COMMON_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wbtan_common.hpp
//! \brief roche-wbtan-1101: the WELL-BALANCED TANGENTIAL potential force, shared by
//! problem generators (deep_hot_jupiter_rt problem/roche_wb_tan; the accretor next), for
//! the cubed sphere AND spherical polar, Hydro AND MHD.
//!
//! The x2/x3 part of -rho grad Phi_tot (centrifugal, tide, ...) is replaced by
//!     S_t = -(rho/rho_eq) R_t(W_eq) = rho g_t,
//! R_t the code's OWN discrete tangential momentum RHS of the problem's equipotential
//! equilibrium W_eq (v = 0; MHD: B = 0, see ProbeT).  Then that equilibrium is balanced in
//! x2/x3 to round-off (seams and sp poles included), whatever the reconstruction,
//! Riemann solver, metric terms or angular-momentum corrections, and for any departure
//! from it the force is rho g_t, g_t = -grad_t Phi_tot + O(h^2).
//!
//! API (one WbTan per problem generator, file scope):
//!   WbTan w;                        on/ready flags, W_eq, g_t, the private exchange
//!   Alloc(w, pmbp)                  W_eq (nmb, 2, n3, n2, n1: rho_eq, eint_eq, ALL cells
//!                                   incl. ghosts) and g_t (nmb, 2, ...: covariant g_2, g_3)
//!   [the pgen fills *w.weq from its own IC formula, fresh start AND restart]
//!   Probe(pm, w, extra, who)        once per run segment, at the first user-source call
//!                                   (w.ready false): runs the fluid's own operators on
//!                                   W_eq and fills g_t; the live state is saved/restored.
//!                                   extra (ExtraFn or nullptr): adds the pgen's own stage
//!                                   terms of W_eq with bdt = 1 to us (IM2/IM3), e.g. a pgen
//!                                   AM fix; scratch is a free u0-shaped array.
//!   ApplyTan(u0, gt, m,k,j,i, rho, a2, a3, bdt)
//!                                   in the user source, ONE call for sp and cs:
//!                                   u0(IM2|IM3) += rho*(a_c + g_c)*bdt, a_c = the pgen's
//!                                   remaining analytic tangential acceleration (Coriolis)
//!                                   projected on the cell's basis by its geometry front-end.
//! Energy: none here; the etotgrav flux books the potential work (the caller requires
//! etotgrav), so E_tot stays exactly conservative.

#include <iostream>
#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "bvals/bvals.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "diffusion/resistivity.hpp"

namespace wbtan {

struct WbTan {
  bool on = false;
  bool ready = false;
  DvceArray5D<Real> *weq = nullptr;   // (nmb, 2, n3, n2, n1): rho_eq, eint_eq (all cells)
  DvceArray5D<Real> *gt = nullptr;    // (nmb, 2, n3, n2, n1): g_2, g_3 (covariant accel.)
  // optional (accwbtan): the FULL conserved W_eq (nmb, 5: IDN..IEN, n3, n2, n1; E as the
  // fluid stores it, etotgrav potential and KE included): a ROTATING equilibrium.  When
  // set it replaces weq in the probe (AllocU).
  DvceArray5D<Real> *ueq = nullptr;
  MeshBoundaryValuesCC *pbv = nullptr;   // the probe's private exchange (own comms)
};

//! the problem generator's own stage terms of W_eq: add them (bdt = 1) to us(IM2/IM3);
//! scratch is a free u0-shaped array (the fluid's u0 itself, restored after the probe)
using ExtraFn = void (*)(Mesh *pm, DvceArray5D<Real> &us, DvceArray5D<Real> &scratch);

inline void Alloc(WbTan &w, MeshBlockPack *pmbp) {
  if (w.weq != nullptr) return;
  auto &ic = pmbp->pmesh->mb_indcs;
  const int n1 = ic.nx1 + 2*ic.ng;
  const int n2 = (ic.nx2 > 1) ? (ic.nx2 + 2*ic.ng) : 1;
  const int n3 = (ic.nx3 > 1) ? (ic.nx3 + 2*ic.ng) : 1;
  w.weq = new DvceArray5D<Real>("wbtan_weq", pmbp->nmb_thispack, 2, n3, n2, n1);
  w.gt = new DvceArray5D<Real>("wbtan_gt", pmbp->nmb_thispack, 2, n3, n2, n1);
}

inline void AllocU(WbTan &w, MeshBlockPack *pmbp) {
  Alloc(w, pmbp);
  if (w.ueq != nullptr) return;
  auto &ic = pmbp->pmesh->mb_indcs;
  const int n1 = ic.nx1 + 2*ic.ng;
  const int n2 = (ic.nx2 > 1) ? (ic.nx2 + 2*ic.ng) : 1;
  const int n3 = (ic.nx3 > 1) ? (ic.nx3 + 2*ic.ng) : 1;
  w.ueq = new DvceArray5D<Real>("wbtan_ueq", pmbp->nmb_thispack, 5, n3, n2, n1);
}

//! wbtan-rstg-1101: on a restart, g_t from the file's per-cell pgen slabs instead of a new
//! probe (a restart must not depend on rebuilding W_eq bit for bit, nor on the probe being
//! side-effect free on a developed state).  True (and w.ready) when the file had exactly
//! the two g_t slabs of this mesh; false (re-probe at the first source call) otherwise.
inline bool LoadFromRestart(ProblemGenerator *pg, WbTan &w, MeshBlockPack *pmbp,
                            const char *who) {
  Alloc(w, pmbp);
  auto &in = pg->pgen_rst_slabs_in;
  const bool ok = (in.extent_int(1) == 2) && (in.extent(0) == w.gt->extent(0)) &&
                  (in.extent(2) == w.gt->extent(2)) && (in.extent(3) == w.gt->extent(3)) &&
                  (in.extent(4) == w.gt->extent(4));
  if (ok) {
    Kokkos::deep_copy(*w.gt, in);
    w.ready = true;
    pg->pgen_rst_slabs = w.gt;
  }
  if (global_variable::my_rank == 0) {
    std::cout << who << ": roche_wb_tan: " << (ok ? "g_t read from the restart file (no "
              "re-probe)" : "the restart file has no g_t: re-probe at the first source "
              "call") << std::endl;
  }
  return ok;
}

//! the force, one code path for every geometry: a2/a3 are the caller's other tangential
//! accelerations on the cell's (covariant) basis
KOKKOS_INLINE_FUNCTION
void ApplyTan(const DvceArray5D<Real> &u0, const DvceArray5D<Real> &gt, const int m,
              const int k, const int j, const int i, const Real rho, const Real a2,
              const Real a3, const Real bdt) {
  u0(m,IM2,k,j,i) += rho*(a2 + gt(m,0,k,j,i))*bdt;
  u0(m,IM3,k,j,i) += rho*(a3 + gt(m,1,k,j,i))*bdt;
}

//----------------------------------------------------------------------------------------
//! \fn ProbeT
//! \brief g_t from the fluid F (hydro::Hydro or mhd::MHD).  The live state (u0, w0, wder,
//! wtemp, uflx, the x1 WB cache; MHD also b0, bcc0) is saved; W_eq fills the active cells
//! and the ghosts are made exactly as the stage makes them (a private
//! MeshBoundaryValuesCC with pbval_u's seam settings, ApplyPhysicalBCs, ConToPrim with no
//! temperature guess so g is restart-invariant); F::Fluxes runs on it (stage 0: no
//! driver; FOFC held off -- it acts only on bad cells, none in an equilibrium); the
//! panel-seam fluxes are reconciled as SendFlux/RecvFlux do; the tangential flux
//! divergence is formed exactly as RKUpdate forms it; the coordinate source is added by
//! the very routine the stage calls (bdt = 1, into a zeroed array), then the fluid's AM
//! module (AmModule::AmTransportProbe) and the caller's extra terms.  Then
//!     g_c = -R_c/rho_eq,  R_c = -div F_c + S_c  (c = IM2, IM3),
//! and everything is restored.
//! MHD: the SAME equipotential state with B = 0 (b0, bcc0 zeroed for the probe, also
//! after the BCs and C2P).  With the field, "the equilibrium" would include the IC's
//! Lorentz force and g_t would freeze its negative into a permanent force although B
//! evolves; the cancelled part is the pressure/potential balance only, and the Lorentz
//! force stays entirely the field's own.

template <class F>
void ProbeT(Mesh *pm, F *ph, WbTan &w, ExtraFn extra, const char *who) {
  constexpr bool kMHD = std::is_same<F, mhd::MHD>::value;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const bool multi_d = pm->multi_d, three_d = pm->three_d;

  auto dup5 = [](const DvceArray5D<Real> &a) {
    DvceArray5D<Real> b(Kokkos::view_alloc("wbtan_sv", Kokkos::WithoutInitializing),
                        a.layout());
    Kokkos::deep_copy(b, a);
    return b;
  };
  auto dup4 = [](const DvceArray4D<Real> &a) {
    DvceArray4D<Real> b(Kokkos::view_alloc("wbtan_sv", Kokkos::WithoutInitializing),
                        a.layout());
    Kokkos::deep_copy(b, a);
    return b;
  };
  const bool hasder = (ph->wder.size() > 0);
  const bool hastmp = (ph->wtemp.size() > 0);
  const bool haswbq = (ph->wbq0.size() > 0);
  DvceArray5D<Real> su0 = dup5(ph->u0), sw0 = dup5(ph->w0);
  DvceArray5D<Real> sder, swbq, sf1 = dup5(ph->uflx.x1f), sf2 = dup5(ph->uflx.x2f);
  DvceArray5D<Real> sf3 = dup5(ph->uflx.x3f);
  DvceArray4D<Real> stmp;
  if (hasder) sder = dup5(ph->wder);
  if (hastmp) stmp = dup4(ph->wtemp);
  if (haswbq) swbq = dup5(ph->wbq0);
  // MHD: the probe is the B = 0 equilibrium (see the note in WbTanProbe)
  DvceArray4D<Real> sb1, sb2, sb3, seta, sadc;
  DvceArray4D<Real> se31, se21, se12, se32, se23, se13;
  DvceArray4D<Real> sr1, sr2, sr3;
  DvceArray5D<Real> sbcc;
  bool sxe = false;
  if constexpr (kMHD) {
    // tracer-1101: a pending x_e Saha start (fresh start / restart without the tracer) is
    // for the live state, not the probe's: kept pending across the probe's ConToPrim
    sxe = ph->xe_init_pending;
    sb1 = dup4(ph->b0.x1f); sb2 = dup4(ph->b0.x2f); sb3 = dup4(ph->b0.x3f);
    // wbtan-etafix-1101: MHD::ConToPrim at stage 0 (as the probe calls it) recomputes the
    // resistivity arrays (Ohmic eta_b, ambipolar ad_coef) from the probe's state, whose
    // passive scalars -- the <mhd>/xe_tracer x_e among them -- are zero: an eta of the
    // probe state.  Saved here and restored below (E_B x450 in the first cycle with
    // xe_tracer + ambipolar = eos, integ-test-1101 gate (b)).
    // wbtan-emffix-1101: MHD::Fluxes also writes the face electric fields e3x1 ... e1x3,
    // which the stage's CornerE/CT read AFTER the user source (where the probe runs): the
    // probe's B = 0 fields replaced this stage's EMF (the stage-1 CT of the first cycle
    // after a fresh start or a restart: restart continuation broke).  Saved and restored.
    se31 = dup4(ph->e3x1); se21 = dup4(ph->e2x1); se12 = dup4(ph->e1x2);
    se32 = dup4(ph->e3x2); se23 = dup4(ph->e2x3); se13 = dup4(ph->e1x3);
    if (ph->presist != nullptr) {
      if (ph->presist->eta_b.size() > 0) seta = dup4(ph->presist->eta_b);
      if (ph->presist->ad_coef.size() > 0) sadc = dup4(ph->presist->ad_coef);
      // wbtan-efldfix-1102: without STS, MHD::Fluxes -> AddResistiveFluxes clears and
      // rebuilds the resistive edge EMF efld_resist from the probe's B = 0; the stage's
      // EField adds it AFTER the user source, so the first stage after a fresh start or
      // a restart had no resistive EMF (no-STS restart continuation broke).
      if (!ph->presist->use_rkg_sts) {
        sr1 = dup4(ph->presist->efld_resist.x1e);
        sr2 = dup4(ph->presist->efld_resist.x2e);
        sr3 = dup4(ph->presist->efld_resist.x3e);
      }
    }
    sbcc = dup5(ph->bcc0);
  }
  auto zero_b = [&]() {
    if constexpr (kMHD) {
      Kokkos::deep_copy(ph->b0.x1f, 0.0); Kokkos::deep_copy(ph->b0.x2f, 0.0);
      Kokkos::deep_copy(ph->b0.x3f, 0.0); Kokkos::deep_copy(ph->bcc0, 0.0);
    }
  };
  const bool sv_built = ph->wb_cache_built;
  const bool sv_fofc = ph->use_fofc;
  // followup-1102 (port of Caltech 0f9dc36ec): two more pieces the probe would leave
  // behind.  (1) Mesh::ecounter: the probe's ConToPrim of W_eq adds its own floor hits,
  //     which would enter the event-log interval of the first cycle of every run segment
  //     (diagnostic only; restored for a clean log).  (2) the conduction module: Fluxes
  //     calls AddHeatFluxes, which (rad_cap_ang / rad_implicit_ang / rad_sts_*) rebuilds
  //     the frozen angular coefficients and caps from the probe state; the heat flux is
  //     energy only and g_t takes the momentum components, so the probe runs without it.
  const EventCounters sv_ec = pm->ecounter;
  Conduction *sv_cond = ph->pcond;
  ph->pcond = nullptr;

  // W_eq -> u0 in the active cells (v = 0; E = eint (+ rho Phi_tot with etotgrav), as
  // the IC), then the ghosts EXACTLY as the live stage gets them: the MeshBlock/panel
  // exchange (own MeshBoundaryValuesCC, own communicators, same seam resample and rho
  // guard as pbval_u), the physical + user x1 BCs, Hydro::ConToPrim
  auto u0 = ph->u0;
  auto weq = *w.weq;
  auto phi0 = ph->phicc0;
  const bool etg = ph->use_etotgrav;
  int nvar;
  if constexpr (kMHD) {nvar = ph->nmhd + ph->nscalars;} else {nvar = ph->nhydro + ph->nscalars;}
  zero_b();
  const bool full = (w.ueq != nullptr);
  DvceArray5D<Real> ueq;
  if (full) ueq = *w.ueq;
  par_for("wbtan_u", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (full) {
      for (int n=IDN; n<=IEN; ++n) u0(m,n,k,j,i) = ueq(m,n,k,j,i);
    } else {
      const Real d = weq(m,0,k,j,i);
      u0(m,IDN,k,j,i) = d;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = weq(m,1,k,j,i) + (etg ? d*phi0(m,k,j,i) : 0.0);
    }
    for (int n=IEN+1; n<nvar; ++n) u0(m,n,k,j,i) = 0.0;
  });
  if (w.pbv == nullptr) {
    w.pbv = new MeshBoundaryValuesCC(pmbp, nullptr, false);
    w.pbv->cs_rho_guard = ph->pbval_u->cs_rho_guard;
    w.pbv->cs_lin_resample = ph->pbval_u->cs_lin_resample;
    w.pbv->cs_blockend_ext_ = ph->pbval_u->cs_blockend_ext_;   // layoutbug-1103
    w.pbv->InitializeBuffers(nvar);
  }
  MeshBoundaryValuesCC *pbv = w.pbv;
  while (pbv->InitRecv(nvar) == TaskStatus::incomplete) {}
  if (pm->use_cubed_sphere) {
    while (pbv->InitFluxSeamRecv(nvar) == TaskStatus::incomplete) {}
  }
  while (pbv->PackAndSendCC(ph->u0, ph->coarse_u0) == TaskStatus::incomplete) {}
  while (pbv->RecvAndUnpackCC(ph->u0, ph->coarse_u0) == TaskStatus::incomplete) {}
  while (pbv->ClearSend() == TaskStatus::incomplete) {}
  while (pbv->ClearRecv() == TaskStatus::incomplete) {}
  (void) ph->ApplyPhysicalBCs(nullptr, 0);
  zero_b();   // a user BC may impose a ghost field; the probe state has none
  // no temperature guess (tguess <= 0: the table solve starts from its own bracket), so
  // W_eq's (T, p, Gamma_1) do not depend on the live wtemp: g is the same on a fresh
  // start and on every restart
  if (hastmp) Kokkos::deep_copy(ph->wtemp, 0.0);
  const bool sv_trk = ph->c2p_track;   // C2PTrack needs the driver
  ph->c2p_track = false;
  (void) ph->ConToPrim(nullptr, 0);
  ph->c2p_track = sv_trk;
  zero_b();
  ph->wb_cache_built = false;   // the x1 WB background of W_eq
  ph->use_fofc = false;
  ph->Fluxes(nullptr, 0);
  // the panel-seam flux reconciliation of SendFlux/RecvFlux
  if (pm->use_cubed_sphere) {
    while (pbv->PackAndSendFluxSeamCC(ph->uflx) == TaskStatus::incomplete) {}
    while (pbv->RecvAndUnpackFluxSeamCC(ph->uflx) == TaskStatus::incomplete) {}
    while (pbv->ClearFluxSend() == TaskStatus::incomplete) {}
    while (pbv->ClearFluxRecv() == TaskStatus::incomplete) {}
  }
  // the coordinate source of W_eq alone, with bdt = 1, into a zeroed array
  DvceArray5D<Real> us("wbtan_us", ph->u0.extent(0), ph->u0.extent(1),
                       ph->u0.extent(2), ph->u0.extent(3), ph->u0.extent(4));
  if constexpr (kMHD) {
    if (pm->use_cubed_sphere) {
      pmbp->pcoord->SrcTermsGnomonicEquiangleMHD(ph->w0, ph->bcc0, ph->wder, ph->pwb,
                                                 ph->uflx, ph->peos->eos_data, 1.0, us);
    } else {
      pmbp->pcoord->SrcTermsSphericalPolarMHD(ph->w0, ph->bcc0, ph->pwb, ph->uflx,
                                              ph->peos->eos_data, 1.0, us);
    }
  } else {
    if (pm->use_cubed_sphere) {
      pmbp->pcoord->SrcTermsGnomonicEquiangle(ph->w0, ph->wder, ph->pwb, ph->uflx,
                                              ph->peos->eos_data, 1.0, us);
    } else {
      pmbp->pcoord->SrcTermsSphericalPolarHydro(ph->w0, ph->pwb, ph->uflx,
                                                ph->peos->eos_data, 1.0, us);
    }
  }
  // the angular-momentum corrections of the same stage, from W_eq's fluxes: the hydro
  // module's (<hydro>/am_transport | am_vector, HydroSrcTerms) and the pgen's own
  // (problem/coriolis_am | am_transport, AmFix in SourceFunc), both with bdt = 1 into
  // `us` (AmFix adds to the fluid's u0, so u0 is pointed at a zeroed scratch for it)
  if (ph->pam != nullptr) {
    DvceArray5D<Real> bcc_;
    if constexpr (kMHD) bcc_ = ph->bcc0;
    hydro::AmFluid fl{us, ph->u1, ph->w0, ph->wder, bcc_, ph->uflx,
                      ph->peos->eos_data, kMHD, ph->use_wellbalance_static, ph->pwb};
    ph->pam->AmTransportProbe(fl);
    // a rotating W_eq: the frame force (Coriolis) of ITS velocity field belongs to the
    // static rotating equilibrium (the am_vector mass-flux torque correction does not)
    ph->pam->FrameForceProbe(fl);
  }
  // the problem generator's own stage terms of W_eq (bdt = 1) into `us`; u0 is free
  // scratch here (restored below)
  if (extra != nullptr) extra(pm, us, ph->u0);
  auto flx1 = ph->uflx.x1f, flx2 = ph->uflx.x2f, flx3 = ph->uflx.x3f;
  auto area1 = pmbp->pcoord->area.x1f, area2 = pmbp->pcoord->area.x2f;
  auto area3 = pmbp->pcoord->area.x3f;
  auto volume = pmbp->pcoord->volume;
  auto w0 = ph->w0;
  auto gt = *w.gt;
  Kokkos::deep_copy(gt, 0.0);
  par_for("wbtan_gt", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    for (int c=0; c<2; ++c) {
      const int n = IM2 + c;
      Real divf = (flx1(m,n,k,j,i+1)*area1(m,k,j,i+1)
                   - flx1(m,n,k,j,i)*area1(m,k,j,i))/volume(m,k,j,i);
      if (multi_d) {
        divf += (flx2(m,n,k,j+1,i)*area2(m,k,j+1,i)
                 - flx2(m,n,k,j,i)*area2(m,k,j,i))/volume(m,k,j,i);
      }
      if (three_d) {
        divf += (flx3(m,n,k+1,j,i)*area3(m,k+1,j,i)
                 - flx3(m,n,k,j,i)*area3(m,k,j,i))/volume(m,k,j,i);
      }
      const Real rhs = -divf + us(m,n,k,j,i);
      gt(m,c,k,j,i) = -rhs/w0(m,IDN,k,j,i);
    }
  });

  // restore the live state
  Kokkos::deep_copy(ph->u0, su0);
  Kokkos::deep_copy(ph->w0, sw0);
  Kokkos::deep_copy(ph->uflx.x1f, sf1);
  Kokkos::deep_copy(ph->uflx.x2f, sf2);
  Kokkos::deep_copy(ph->uflx.x3f, sf3);
  if (hasder) Kokkos::deep_copy(ph->wder, sder);
  if (hastmp) Kokkos::deep_copy(ph->wtemp, stmp);
  if (haswbq) Kokkos::deep_copy(ph->wbq0, swbq);
  if constexpr (kMHD) {
    ph->xe_init_pending = sxe;
    Kokkos::deep_copy(ph->e3x1, se31); Kokkos::deep_copy(ph->e2x1, se21);
    Kokkos::deep_copy(ph->e1x2, se12); Kokkos::deep_copy(ph->e3x2, se32);
    Kokkos::deep_copy(ph->e2x3, se23); Kokkos::deep_copy(ph->e1x3, se13);
    if (seta.size() > 0) Kokkos::deep_copy(ph->presist->eta_b, seta);
    if (sadc.size() > 0) Kokkos::deep_copy(ph->presist->ad_coef, sadc);
    if (sr1.size() > 0) {
      Kokkos::deep_copy(ph->presist->efld_resist.x1e, sr1);
      Kokkos::deep_copy(ph->presist->efld_resist.x2e, sr2);
      Kokkos::deep_copy(ph->presist->efld_resist.x3e, sr3);
    }
    Kokkos::deep_copy(ph->b0.x1f, sb1); Kokkos::deep_copy(ph->b0.x2f, sb2);
    Kokkos::deep_copy(ph->b0.x3f, sb3); Kokkos::deep_copy(ph->bcc0, sbcc);
  }
  ph->wb_cache_built = sv_built;
  ph->use_fofc = sv_fofc;
  ph->pcond = sv_cond;
  pm->ecounter = sv_ec;
  w.ready = true;
  // wbtan-rstg-1101: every restart file from now on carries g_t (kPgenSlabRstMagic)
  if (pm->pgen != nullptr) pm->pgen->pgen_rst_slabs = w.gt;

  // one line: the size of g_t (max |g_c| over the rank's cells, rank 0)
  Real gmax = 0.0;
  Kokkos::parallel_reduce("wbtan_gmax",
      Kokkos::RangePolicy<DevExeSpace>(0, static_cast<int>(gt.size())),
      KOKKOS_LAMBDA(const int q, Real &mx) {
        mx = fmax(mx, fabs(gt.data()[q]));
      }, Kokkos::Max<Real>(gmax));
  if (global_variable::my_rank == 0) {
    std::cout << who << ": roche_wb_tan: equilibrium residual probe "
              << "done at cycle " << pm->ncycle << ", max |g_t| (rank 0) = " << gmax
              << " cm/s^2" << std::endl;
  }
}


//! dispatch on the pack's fluid (MHD first, as the problem generators do)
inline void Probe(Mesh *pm, WbTan &w, ExtraFn extra, const char *who) {
  if (pm->pmb_pack->pmhd != nullptr) {
    ProbeT(pm, pm->pmb_pack->pmhd, w, extra, who);
  } else {
    ProbeT(pm, pm->pmb_pack->phydro, w, extra, who);
  }
}

}  // namespace wbtan

#endif  // PGEN_WBTAN_COMMON_HPP_
