//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file he_strip_bg.cpp
//! \brief stripbg-1102 (USER 11-02, He giant global): a ONE-SHOT restart-time transform
//! of the he_star_m1 problem that removes the gas above the photosphere (lifted gas
//! falling back, filling the volume to the outer boundary) and replaces it with a
//! static, optically thin ambient in radiative equilibrium.
//!
//! Keys (<problem>, all read only when strip_bg_on_restart = true; absent = bitwise):
//!   strip_bg_on_restart (false)  apply the transform once, at this restart
//!   strip_bg_time       (none)   if named: apply only when the restart time equals it
//!                                (relative 1e-9), so an inherited key never re-fires
//!   strip_bg_tau        (2/3)    the cut: per column, the radial Rosseland optical depth
//!                                tau(r) = int_r^{r_top} (kappa_F + kappa_s) rho dr from
//!                                the outer mesh edge inward, at the cell centre; cell
//!                                i_c = the outermost active cell with tau >= strip_bg_tau.
//!                                Every cell i <= i_c is kept untouched; every cell
//!                                i > i_c (and the outer x1 ghosts) is REPLACED.
//!   strip_bg_rho_fac    (1e-3)   rho_amb(r) = max(dfloor, rho_0 (r_c/r)^2),
//!   strip_bg_rho_max    (2e-14)  rho_0 = min(fac rho(i_c), rho_max) (code density)
//!   strip_bg_nblend     (3)      cells i_c+1..i_c+nb: log-linear blend in rho and T from
//!                                the kept cell's (rho, T) to the ambient (rho_amb(r), T_amb)
//!                                -- a NEW profile (no old gas above i_c is kept or
//!                                averaged), only so that P has no jump at r_c
//! Replaced cells: v = 0, T = T_amb = max(tfloor, (E/a_rad)^(1/4)) of the restart's own
//! M1 E (radiative equilibrium with the star's light; E and F are kept as they are),
//! e = e(rho, T) from the EOS, IEN = e + rho Phi (etotgrav), passive scalars 0, am spin 0;
//! the general-EOS restart caches (wtemp, wder p/Gamma_1, eint_rst) of those cells are
//! set to the new state so the first (frozen) ConToPrim is consistent and the kept
//! cells stay bitwise the unmodified restart's.
//! Columns: every (k,j) of each MeshBlock incl. the x2/x3 ghost columns (the restart keeps
//! ghosts), the corner columns (both ghost) excepted.  A column whose tau never reaches
//! strip_bg_tau is left untouched (counted).  After the transform the key is set false in
//! the ParameterInput, so the restart files this run writes do not carry it.
//========================================================================================

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "utils/eint_from_cons.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

void HsStripBackground(ParameterInput *pin, Mesh *pm) {
  if (!pin->DoesParameterExist("problem", "strip_bg_on_restart") ||
      !pin->GetBoolean("problem", "strip_bg_on_restart")) return;   // absent: nothing added
  if (pin->DoesParameterExist("problem", "strip_bg_time")) {
    const Real ts = pin->GetReal("problem", "strip_bg_time");
    if (std::fabs(pm->time - ts) > 1.0e-9*std::max(std::fabs(ts), 1.0)) {
      if (global_variable::my_rank == 0) {
        std::cout << "strip_bg: restart time " << pm->time << " != strip_bg_time " << ts
                  << ": NOT applied" << std::endl;
      }
      pin->SetBoolean("problem", "strip_bg_on_restart", false);
      return;
    }
  }
  const Real taucut = pin->GetOrAddReal("problem", "strip_bg_tau", 2.0/3.0);
  const Real rfac = pin->GetOrAddReal("problem", "strip_bg_rho_fac", 1.0e-3);
  const Real rmax = pin->GetOrAddReal("problem", "strip_bg_rho_max", 2.0e-14);
  const int nb = pin->GetOrAddInteger("problem", "strip_bg_nblend", 3);

  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  auto *pm1 = pmbp->pradm1;
  if (ph == nullptr || pm1 == nullptr) {
    std::cout << "### FATAL strip_bg: needs <hydro> and <rad_m1>" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 != pm->mesh_indcs.nx1) {
    std::cout << "### FATAL strip_bg: needs ONE MeshBlock across x1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = indcs.nx2 + 2*ng - 1, n3m1 = indcs.nx3 + 2*ng - 1;
  const int nmb = pmbp->nmb_thispack;
  const bool hcs = pm->use_cubed_sphere;
  const bool etg = ph->use_etotgrav;
  const int nhyd = ph->nhydro, nscal = ph->nscalars;

  auto u0 = ph->u0;
  auto wtemp = ph->wtemp;
  auto wder = ph->wder;
  const bool hwt = (wtemp.extent_int(0) > 0), hwd = (wder.extent_int(0) > 0);
  const bool hei = ph->c2p_eint_rst && (ph->eint_rst.extent_int(0) > 0);
  auto eir = ph->eint_rst;
  auto phi = ph->phicc0;
  auto ur = pm1->u0;
  auto ot = pm1->otab;
  const Real ar = pm1->arad;
  auto eos = ph->peos->eos_data;
  const Real dfl = eos.dfloor, tfl = eos.tfloor;
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto cosc = pmbp->pcoord->cos_cell;
  auto &volume = pmbp->pcoord->volume;
  hydro::AmModule *pam = hydro::AmModuleOf(pmbp);
  const bool hsp = (pam != nullptr) && (pam->am_spin.extent_int(0) > 0);
  DvceArray4D<Real> spin;
  if (hsp) spin = pam->am_spin;

  // per column: index of the kept cell (-1 = untouched), r_c, mass removed / added
  DvceArray3D<Real> rc("sbg_rc", nmb, n3m1 + 1, n2m1 + 1);
  DvceArray3D<int> icut("sbg_ic", nmb, n3m1 + 1, n2m1 + 1);
  DvceArray3D<Real> dmo("sbg_dmo", nmb, n3m1 + 1, n2m1 + 1);
  DvceArray3D<Real> dmn("sbg_dmn", nmb, n3m1 + 1, n2m1 + 1);

  par_for("strip_bg", DevExeSpace(), 0, nmb - 1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    icut(m,k,j) = -2;
    rc(m,k,j) = 0.0;
    dmo(m,k,j) = 0.0;
    dmn(m,k,j) = 0.0;
    const bool gj = (j < js || j > je), gk = (k < ks || k > ke);
    if (gj && gk) return;                       // corner column: not used
    const Real c = hcs ? cosc(m,k,j) : 0.0;
    auto tgas = [&](const int i) {
      const Real d = u0(m,IDN,k,j,i);
      Real ke2 = hcs ? CsKeMom(u0(m,IM1,k,j,i), u0(m,IM2,k,j,i), u0(m,IM3,k,j,i), c)
                     : 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                            + SQR(u0(m,IM3,k,j,i)));
      Real e = u0(m,IEN,k,j,i) - ke2/d - (etg ? d*phi(m,k,j,i) : 0.0);
      const Real tg = hwt ? wtemp(m,k,j,i) : -1.0;
      return eos.Temperature(d, e, tg);
    };
    // tau from the outer edge inward, at cell centres
    Real tau = 0.0;
    int ic = -1;
    Real tc = 0.0;
    for (int i = ie; i >= is; --i) {
      const Real d = u0(m,IDN,k,j,i);
      const Real t = tgas(i);
      Real op, oe, of, os;
      radm1::M1TableOpacities(ot, d, t, op, oe, of, os);
      const Real dt = (of + os)*d*(x1f(m,i+1) - x1f(m,i));
      if (tau + 0.5*dt >= taucut) {
        ic = i;
        tc = t;
        break;
      }
      tau += dt;
    }
    icut(m,k,j) = ic;
    if (ic < 0) return;                         // never optically thick: untouched
    const Real r0 = x1v(m,ic);
    rc(m,k,j) = r0;
    const Real dc = u0(m,IDN,k,j,ic);
    const Real d0 = fmin(rfac*dc, rmax);
    for (int i = ic + 1; i <= n1m1; ++i) {
      const Real r = x1v(m,i);
      const Real damb = fmax(dfl, d0*SQR(r0/r));
      const Real er = ur(m,0,k,j,i);
      const Real tamb = fmax(tfl, (er > 0.0) ? sqrt(sqrt(er/ar)) : 0.0);
      Real d = damb, t = tamb;
      const int n = i - ic;
      if (n <= nb) {
        const Real w = static_cast<Real>(n)/static_cast<Real>(nb + 1);
        d = exp((1.0 - w)*log(dc) + w*log(damb));
        t = exp((1.0 - w)*log(tc) + w*log(tamb));
      }
      if (i <= ie && !gj && !gk) {
        dmo(m,k,j) += u0(m,IDN,k,j,i)*volume(m,k,j,i);
        dmn(m,k,j) += d*volume(m,k,j,i);
      }
      const Real e = eos.EnergyFromTemperature(d, t);
      u0(m,IDN,k,j,i) = d;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      u0(m,IEN,k,j,i) = e + (etg ? d*phi(m,k,j,i) : 0.0);
      for (int s = 0; s < nscal; ++s) u0(m,nhyd+s,k,j,i) = 0.0;
      if (hsp) spin(m,k,j,i) = 0.0;
      if (hwt) wtemp(m,k,j,i) = t;
      if (hwd) {
        Real p, g1;
        eos.PressureAndGamma1(d, e, t, p, g1);
        wder(m,IDPR,k,j,i) = p;
        wder(m,IDG1,k,j,i) = g1;
      }
      if (hei) eir(m,k,j,i) = e;
    }
  });

  // report (active columns only)
  auto hic = Kokkos::create_mirror_view_and_copy(HostMemSpace(), icut);
  auto hrc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rc);
  auto hmo = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dmo);
  auto hmn = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dmn);
  std::vector<double> rcs;
  double sums[4] = {0.0, 0.0, 0.0, 0.0};   // cols, untouched cols, M removed, M added
  for (int m = 0; m < nmb; ++m) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        sums[0] += 1.0;
        if (hic(m,k,j) < 0) {
          sums[1] += 1.0;
          continue;
        }
        rcs.push_back(hrc(m,k,j));
        sums[2] += hmo(m,k,j);
        sums[3] += hmn(m,k,j);
      }
    }
  }
  double rmn = 1.0e300, rmx = -1.0e300;
  for (double r : rcs) {
    rmn = std::min(rmn, r);
    rmx = std::max(rmx, r);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &rmn, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &rmx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    std::cout << "strip_bg (stripbg-1102): APPLIED at t = " << pm->time << ", cycle "
              << pm->ncycle << ": tau_cut = " << taucut << ", rho_fac = " << rfac
              << ", rho_max = " << rmax << ", nblend = " << nb << std::endl
              << "  columns " << sums[0] << " (untouched, never tau >= tau_cut: "
              << sums[1] << "), r_cut min " << rmn << " max " << rmx
              << ", mass removed " << sums[2] << " added " << sums[3] << std::endl;
  }
  // one-shot: never inherited by the restart files this run writes
  pin->SetBoolean("problem", "strip_bg_on_restart", false);
}

//----------------------------------------------------------------------------------------
//! rstpert-1105 (USER 11-05, HG prodG stall at t 4.4232e7): a ONE-SHOT restart-time tiny
//! relative perturbation of the M1 radiation energy E of every active cell,
//!   E -> E (1 + amp (2 h - 1)),  h = hash(gid, k, j, i) in [0, 1) (layout-independent),
//! to move the trajectory off a deterministic stall.  Keys (<problem>):
//!   rst_pert_erad (absent = nothing, bitwise)  the relative amplitude amp (>= 0)
//!   rst_pert_time (required with it)           apply only when the restart time equals it
//!                                              (relative 1e-9): an inherited key never re-fires
void HsRestartPerturb(ParameterInput *pin, Mesh *pm) {
  if (!pin->DoesParameterExist("problem", "rst_pert_erad")) return;   // absent: bitwise
  const Real amp = pin->GetReal("problem", "rst_pert_erad");
  if (!(amp > 0.0)) return;
  if (!pin->DoesParameterExist("problem", "rst_pert_time")) {
    std::cout << "### FATAL rst_pert_erad needs problem/rst_pert_time" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real ts = pin->GetReal("problem", "rst_pert_time");
  if (std::fabs(pm->time - ts) > 1.0e-9*std::max(std::fabs(ts), 1.0)) {
    if (global_variable::my_rank == 0) {
      std::cout << "rst_pert: restart time " << pm->time << " != rst_pert_time " << ts
                << ": NOT applied" << std::endl;
    }
    return;
  }
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto *pm1 = pmbp->pradm1;
  if (pm1 == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack, gids = pmbp->gids;
  auto ur = pm1->u0;
  par_for("rst_pert_erad", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    uint64_t x = (static_cast<uint64_t>(m + gids) << 48) ^
                 (static_cast<uint64_t>(k) << 32) ^ (static_cast<uint64_t>(j) << 16) ^
                 static_cast<uint64_t>(i);
    x += 0x9E3779B97F4A7C15ULL;   // splitmix64
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    x = x ^ (x >> 31);
    const Real h = static_cast<Real>(x >> 11) * (1.0/9007199254740992.0);
    ur(m,0,k,j,i) *= (1.0 + amp*(2.0*h - 1.0));
  });
  Kokkos::fence();
  if (global_variable::my_rank == 0) {
    std::cout << "rst_pert: M1 E of every active cell perturbed by relative " << amp
              << " (hash noise) at t = " << pm->time << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! eosreinit-1112 (USER 11-12, RSG eos_h2 switch at a restart): problem/eos_reinit_from_t
//! (absent or false = nothing, bitwise).  ONE-SHOT, on a restart: re-set the internal
//! energy from the restart's (rho, T) with THIS run's EOS table, keeping T.  The general
//! EOS measures e from a composition-dependent zero point (with <hydro>/eos_h2 = true it
//! is "all H in ground-state H2", so every atomic H carries +chi_d/2), so restarting with
//! a different table at the file's e jumps T (eos_h2 false -> true: T collapses into the
//! molecular branch).  T is the restart file's temperature cache (Hydro::wtemp, restored
//! by the restart reader before UserProblem; the run needs it, the call fatals without
//! it).  Per cell (ghosts included): u0(IEN) += e_new(rho, T) - e_old, with e_old the
//! file's internal energy (the eint_rst cache when the file carried it, else u0(IEN)
//! - KE - rho Phi with the cubed-sphere covariant KE); KE, rho Phi, rho, momenta and the
//! M1 radiation (E, F) are untouched.  The derived cache (p, Gamma_1) is re-evaluated
//! with the new table at (rho, e_new, T) and eint_rst set to e_new, so the frozen first
//! ConToPrim of the restart (Hydro::c2p_freeze_derived) carries the NEW table's state.
//! The total internal-energy change over the active cells (= the total-energy change,
//! etotgrav and radiation unchanged) is volume-integrated and printed, with the pressure
//! change at fixed (rho, T) (the molecular-region mu change).  Hydro only (the pgen
//! refuses <mhd>); refused under <hydro>/dual_energy (its pressure slots belong to the
//! old table).  The key is cleared after use, so the restart files this run writes never
//! re-fire it.
void HsEosReinitFromT(ParameterInput *pin, Mesh *pm) {
  if (!pin->DoesParameterExist("problem", "eos_reinit_from_t") ||
      !pin->GetBoolean("problem", "eos_reinit_from_t")) return;   // absent: bitwise
  MeshBlockPack *pmbp = pm->pmb_pack;
  hydro::Hydro *ph = pmbp->phydro;
  if (ph == nullptr || pmbp->pmhd != nullptr) {
    std::cout << "### FATAL eos_reinit_from_t: hydro only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto eos = ph->peos->eos_data;
  auto wtemp = ph->wtemp;
  auto wder = ph->wder;
  if (!eos.IsGeneral() || wtemp.extent_int(0) == 0 || wder.extent_int(0) == 0 ||
      !ph->c2p_freeze_derived) {
    std::cout << "### FATAL eos_reinit_from_t needs the general EOS and a restart file "
              << "carrying its T and derived caches" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ph->dual.on) {
    std::cout << "### FATAL eos_reinit_from_t is not supported with <hydro>/dual_energy"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb = pmbp->nmb_thispack;
  const bool hcs = pm->use_cubed_sphere;
  const bool etg = ph->use_etotgrav;
  const bool hei = ph->c2p_eint_rst && (ph->eint_rst.extent_int(0) > 0);
  auto eir = ph->eint_rst;
  auto u0 = ph->u0;
  auto phi = ph->phicc0;
  auto cosc = pmbp->pcoord->cos_cell;
  auto &volume = pmbp->pcoord->volume;
  // sums over the ACTIVE cells: 0 dE_int, 1 E_int old, 2 cells changed, 3 int p_old dV,
  // 4 int p_new dV, 5 cells p_new/p_old < 0.9, 6 int rho dV of those, 7 int rho dV;
  // 8 min p_new/p_old (as -max of -ratio), 9 skipped cells (no T or e)
  constexpr int NS = 10;
  DvceArray1D<Real> sums("eosre_sums", NS);
  Kokkos::deep_copy(sums, 0.0);
  Kokkos::deep_copy(Kokkos::subview(sums, 8), -1.0e300);
  par_for("eos_reinit_from_t", DevExeSpace(), 0, nmb - 1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
    const Real d = u0(m,IDN,k,j,i);
    const Real t = wtemp(m,k,j,i);
    Real eold;
    if (hei) {
      eold = eir(m,k,j,i);
    } else {
      const Real c = hcs ? cosc(m,k,j) : 0.0;
      const Real ke2 = hcs ? CsKeMom(u0(m,IM1,k,j,i), u0(m,IM2,k,j,i),
                                     u0(m,IM3,k,j,i), c)
                           : 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                                  + SQR(u0(m,IM3,k,j,i)));
      eold = u0(m,IEN,k,j,i) - ke2/d - (etg ? d*phi(m,k,j,i) : 0.0);
    }
    if (!(d > 0.0) || !(t > 0.0) || !(eold > 0.0)) {
      if (act) Kokkos::atomic_add(&sums(9), 1.0);
      return;
    }
    const Real enew = eos.EnergyFromTemperature(d, t);
    const Real pold = wder(m,IDPR,k,j,i);
    u0(m,IEN,k,j,i) += enew - eold;
    Real p, g1;
    eos.PressureAndGamma1(d, enew, t, p, g1);
    wder(m,IDPR,k,j,i) = p;
    wder(m,IDG1,k,j,i) = g1;
    if (hei) eir(m,k,j,i) = enew;
    if (act) {
      const Real v = volume(m,k,j,i);
      Kokkos::atomic_add(&sums(0), (enew - eold)*v);
      Kokkos::atomic_add(&sums(1), eold*v);
      Kokkos::atomic_add(&sums(2), 1.0);
      Kokkos::atomic_add(&sums(3), pold*v);
      Kokkos::atomic_add(&sums(4), p*v);
      Kokkos::atomic_add(&sums(7), d*v);
      const Real rt = (pold > 0.0) ? p/pold : 1.0;
      if (rt < 0.9) {
        Kokkos::atomic_add(&sums(5), 1.0);
        Kokkos::atomic_add(&sums(6), d*v);
      }
      Kokkos::atomic_max(&sums(8), -rt);
    }
  });
  Kokkos::fence();
  auto hs = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sums);
  double s[NS];
  for (int n = 0; n < NS; ++n) s[n] = static_cast<double>(hs(n));
  double smx = s[8];
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, s, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &s[9], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &smx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0) {
    std::printf("eos_reinit_from_t (eosreinit-1112): APPLIED at t = %.10e, cycle %d "
                "(e from the restart (rho, T) with this run's EOS table; T, rho, "
                "momenta, M1 E/F kept; e_old from %s)\n", pm->time, pm->ncycle,
                hei ? "the restart's eint cache" : "the conserved state");
    std::printf("eos_reinit_from_t: active cells %.0f (skipped %.0f); int E_int dV old "
                "%.10e new %.10e erg; dE_total = dE_int = %.10e erg (%.4e of E_int "
                "old)\n",
                s[2], s[9], s[1], s[1] + s[0], s[0], (s[1] != 0.0) ? s[0]/s[1] : 0.0);
    std::printf("eos_reinit_from_t: int p dV old %.10e new %.10e; cells with p_new/p_old "
                "< 0.9: %.0f (mass %.6e of %.6e g); min p_new/p_old %.6f\n",
                s[3], s[4], s[5], s[6], s[7], -smx);
  }
  // one-shot: never inherited by the restart files this run writes
  pin->SetBoolean("problem", "eos_reinit_from_t", false);
}
