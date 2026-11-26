//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_imex.cpp
//! \brief time_scheme = hesdirk2: the IMEX stages (Heun hydro, ESDIRK implicit
//! radiation), their state and its restart channels.
//! Implicit part: A = [[0,0,0],[1-g,g,0],[1/2,b2,g]], g = 1 - 1/sqrt 2, b2 = g/(2(1-g)),
//! first slope FSAL.  A stage is the backward-Euler ImplicitSolve with dt -> g dt and the
//! old vector rhs = (stage start) + t2inc, t2inc = (1-g) dt K1 (stage 1) or
//! dt (g/2 K1 + (b2 - g/2) K2) (stage 2); its slope is K = (Y - rhs)/(g dt).  Without a
//! valid K1 (first step, restart without the slope, failed stage) the step is backward
//! Euler and stores K1 = (Y - rhs)/dt.
//!
//!  called from: rad_vet_gas.cpp, rays/sc_cart_sweep.cpp, rays/sc_cart_sweep_mb.cpp,
//!    rays/sc_directions.cpp, rays/sc_source.cpp, driver/driver.cpp
//!  calls into: rad_vet_fct.cpp, rad_vet_gas.cpp, rad_vet_ledger.cpp, rad_vet_opacity.cpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "utils/deep_copy_across.hpp"
#include "utils/eint_from_cons.hpp"
#include "utils/imex_sdirk2.hpp"

namespace radvet {

using deep_copy_across::DeepCopyAcross;

namespace {

const Real kT2G = imex_sdirk2::Gamma();   // shared tableau

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2Active
//! \brief true when the next step runs the two stage solves

bool RadiationVet::Time2Active() {
  if (time_scheme != VET_TIME_HESDIRK2) return false;
  if (static_cast<int>(t2k1.extent(0)) < pmy_pack->nmb_thispack) {t2_ok = false;}
  return t2_ok;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2GasTime
//! \brief time2_tableau = sdirk2: one time per stage for gas and radiation, which keeps
//! the IMEX coupling conditions b_E.c_I = b_I.c_E = 1/2 with the Heun hydro.
//!   which = 1 (before the stage-1 solve): the Heun-stage-1 gas U* is kept in t2hs; the
//!     solve takes S = U^n + g (U* - U^n), the gas at t^n + g dt.
//!   which = 2 (first in GasUpdate, gas Y = S + d): the gas becomes U* + d.
//!   which = 3 (after the stage tasks): w0 is refreshed from the hydro's stage-2 flux
//!     state U* + f d (f = 1/g non-stiff, -> 1 stiff), then u0 is U* + d again.
//! time2_gtime_mom: momentum's f = 1/g in every cell; IEN takes its kinetic-energy
//! change.  The flux state's eint increment (hs_ IEN, raw) is (f - 1) d's, or with
//! time2_gtime_slope (1-g) dt minmod(K1, K2), K1 the last step's stage-2 slope of eint
//! and K2 = d/(g dt) (f not used for the energy).  At which = 3 it is added to the gas
//! the caps and floors left, in the Patankar form where it is negative (positive).  Its
//! floors act on the primitives only: u0 is restored to U* + d as it was, so they
//! create no energy and nothing is paid.
//! time2_gtime_pos = false: the old flux state U* + f d,
//! formed at which = 2 on the gas before the caps, its floors kept in u0.
//! <mhd>: the solve takes B from S as it takes the momentum, B^n + g (B* - B^n), through
//! its |B|^2/2 (emag0, held from which = 1 to 2); d of B is 0, so the gas U* + d and the
//! flux state have B*, and the CT update of the hydro's stage 2 takes its EMFs from the
//! flux state as the gas fluxes (Heun weights).  b0 is not touched.

//----------------------------------------------------------------------------------------

void RadiationVet::Time2GasTime(Driver *pdrive, int which) {
  if (t2_tab != 2 || !coupling || !fl_on || (which > 1 && t2_gt != which - 1)) return;
  t2_gt = (which == 3) ? 0 : which;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto u0_ = fl.u0;
  auto u1_ = fl.u1;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nmb = u0_.extent_int(0), nv = u0_.extent_int(1);
  const int n3 = u0_.extent_int(2), n2 = u0_.extent_int(3), n1 = u0_.extent_int(4);
  const bool all = (which == 1);
  const int il = all ? 0 : id.is, iu = all ? n1-1 : id.ie, jl = all ? 0 : id.js;
  const int ju = all ? n2-1 : id.je, kl = all ? 0 : id.ks, ku = all ? n3-1 : id.ke;
  if (t2hs.extent_int(0) != nmb || t2hs.extent_int(1) != nv) {
    Kokkos::realloc(t2hs, nmb, nv, n3, n2, n1);
  }
  auto hs_ = t2hs;
  // flux-state factor on d: 1/g + O(k dt), -> 1 where the exchange is stiff, k = ch rho
  // kappa_E + 4 cl rho kappa_P a T^3/(rho c_v)
  auto op_ = opac;
  auto eos = fl.eos;
  auto phc = fl.phicc0;
  auto em_ = emag0;
  auto ccl = pmy_pack->pcoord->cos_cell;
  const bool etg = fl.use_etotgrav, csk = cs_geom, mh = fl.mhd;
  const Real g = kT2G, cl = c_light, ch = chat, ar = arad, dt = pmy_pack->pmesh->dt;
  if (which == 1 && mh) {
    EmagBuild(g);   // the solve's field: S's
    emag_hold = true;
  }
  const bool gs = t2_gts, ps = t2_gtp;
  auto k1_ = t2k1;
  par_for("vet_t2_gtime", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = 0.0, xs = 0.0;
    if (which == 2) {
      const Real dd = u0_(m,IDN,k,j,i);
      const Real eg = u0_(m,IEN,k,j,i) - CsKeMom(u0_(m,IM1,k,j,i), u0_(m,IM2,k,j,i),
                      u0_(m,IM3,k,j,i), csk ? ccl(m,k,j) : 0.0)/dd -
                      (etg ? dd*phc(m,k,j,i) : 0.0) - (mh ? em_(m,k,j,i) : 0.0);
      const Real tk = (eg > 0.0) ? eos.Temperature(dd, eg) : 0.0;
      const Real z = dt*(ch*op_(m,VET_OP_E,k,j,i) + ((eg > 0.0) ?
                     4.0*cl*op_(m,VET_OP_P,k,j,i)*ar*SQR(SQR(tk))/eg : 0.0));
      e = (1.0 + g*z)/(g*(1.0 + z)) - 1.0;
      if (gs) {
        const Real cc = csk ? ccl(m,k,j) : 0.0, r = u1_(m,IDN,k,j,i);
        Real p[3], q[3], s[3];   // momenta of U^n, U^n + g dt K1 and S
        for (int n = 0; n < 3; ++n) {
          p[n] = u1_(m,IM1+n,k,j,i);
          q[n] = p[n] + g*dt*k1_(m,VET_T2_M1+n,k,j,i);
          s[n] = p[n] + g*(hs_(m,IM1+n,k,j,i) - p[n]);
        }
        // the eint changes g dt K1 and g dt K2 = d
        const Real a = g*dt*k1_(m,VET_T2_EN,k,j,i) -
                       (CsKeMom(q[0], q[1], q[2], cc) - CsKeMom(p[0], p[1], p[2], cc))/r;
        const Real es = u1_(m,IEN,k,j,i) + g*(hs_(m,IEN,k,j,i) - u1_(m,IEN,k,j,i)) -
                        CsKeMom(s[0], s[1], s[2], cc)/dd -
                        (etg ? dd*phc(m,k,j,i) : 0.0) - (mh ? em_(m,k,j,i) : 0.0);
        const Real b = eg - es;
        xs = (a*b > 0.0) ? (1.0/g - 1.0)*((fabs(a) < fabs(b)) ? a : b) : 0.0;
        e = 0.0;
      }
    }
    if (which == 3 && ps) {
      const Real c = csk ? ccl(m,k,j) : 0.0, r = u0_(m,IDN,k,j,i);
      const Real k0 = CsKeMom(u0_(m,IM1,k,j,i), u0_(m,IM2,k,j,i), u0_(m,IM3,k,j,i), c)/r;
      const Real k1 = CsKeMom(u0_(m,IM1,k,j,i) + hs_(m,IM1,k,j,i), u0_(m,IM2,k,j,i) +
                              hs_(m,IM2,k,j,i), u0_(m,IM3,k,j,i) + hs_(m,IM3,k,j,i), c)/r;
      const Real eb = u0_(m,IEN,k,j,i) - k0 - (etg ? r*phc(m,k,j,i) : 0.0) -
                      (mh ? em_(m,k,j,i) : 0.0);
      xs = hs_(m,IEN,k,j,i);
      hs_(m,IEN,k,j,i) = k1 - k0 + ((eb > 0.0) ? ((xs < 0.0) ? xs*eb/(eb - xs) : xs) :
                                    0.0);
    }
    for (int n = 0; n < nv; ++n) {
      const Real un = u1_(m,n,k,j,i), us = hs_(m,n,k,j,i), u = u0_(m,n,k,j,i);
      if (which == 1) {
        hs_(m,n,k,j,i) = u;
        u0_(m,n,k,j,i) = un + g*(u - un);
      } else if (which == 2) {
        const Real d = u - (un + g*(us - un));   // (0 where the solve left it)
        u0_(m,n,k,j,i) = us + d;
        hs_(m,n,k,j,i) = ((n >= IM1 && n <= IM3) ? 1.0/g - 1.0 : e)*d;
      } else {
        u0_(m,n,k,j,i) = u + us;
        if (ps) {hs_(m,n,k,j,i) = u;}   // (the register to restore)
      }
    }
    if (which == 2 && gs) {
      hs_(m,IEN,k,j,i) = xs;
    } else if (which == 2) {
      // ps: (f - 1) d of the total energy less that of the kinetic; else IEN's with the
      // momentum's f
      const Real c = csk ? ccl(m,k,j) : 0.0, r = u0_(m,IDN,k,j,i) + hs_(m,IDN,k,j,i);
      const Real s = e*g/(1.0 - g);
      auto ke = [&](const Real q) {
        return CsKeMom(u0_(m,IM1,k,j,i) + q*hs_(m,IM1,k,j,i), u0_(m,IM2,k,j,i) +
                       q*hs_(m,IM2,k,j,i), u0_(m,IM3,k,j,i) + q*hs_(m,IM3,k,j,i), c)/r;
      };
      hs_(m,IEN,k,j,i) += ps ? ke(0.0) - ke(s) : ke(1.0) - ke(s);
    }
  });
  if (which == 2 && mh) {
    emag_hold = false;
    EmagBuild(1.0);   // B* (U* + d)
  }
  if (which != 3) return;
  fluid::DualEnergy &du = mh ? fl.pm->dual : fl.ph->dual;
  MeshBoundaryValuesCC *pbu = mh ? fl.pm->pbval_u : fl.ph->pbval_u;
  auto cu0 = mh ? fl.pm->coarse_u0 : fl.ph->coarse_u0;
  const bool dpay = du.pay;
  if (ps) {du.pay = du.book = false;}
  FluidRefresh(pdrive, false);   // ghosts + w0 of the flux state
  du.pay = dpay;
  du.book = true;
  // u0 back to U* + d, in the ghosts the neighbour's (exchanged; the flux state at the
  // physical boundaries): a ghost is then its neighbour's cell bitwise, so the stage-2
  // dfloor_flux_limit takes the same donor state on both sides of a MeshBlock face (the
  // face flux does not depend on the MeshBlock layout)
  par_for("vet_t2_gtime3g", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (k >= kl && k <= ku && j >= jl && j <= ju && i >= il && i <= iu) return;
    for (int n = 0; n < nv; ++n) {hs_(m,n,k,j,i) = ps ? u0_(m,n,k,j,i) : 0.0;}
  });
  (void) pbu->InitRecv(nv);
  (void) pbu->PackAndSendCC(hs_, cu0);
  (void) pbu->ClearSend();
  (void) pbu->ClearRecv();
  (void) pbu->RecvAndUnpackCC(hs_, cu0);
  par_for("vet_t2_gtime3", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int n = 0; n < nv; ++n) {
      u0_(m,n,k,j,i) = ps ? hs_(m,n,k,j,i) : u0_(m,n,k,j,i) - hs_(m,n,k,j,i);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2FormStage
//! \brief after the hydro's stage `stage`: the radiation stage start (stage 1 U^n, stage
//! 2 the Heun average of U^n and Y2) and t2inc, which ImplicitSolve adds to form the old
//! vector of the stage solve

void RadiationVet::Time2FormStage(int stage, Real dt) {
  const Real g = kT2G;
  const Real w0 = t2_w0, w1 = t2_w1;
  // hesdirk2 and sdirk2 take the Heun stage start; trbdf2 the w0/w1 average
  const bool heun = (t2_tab == 0 || t2_tab == 2);
  Real c1, c2;
  auto u0_ = u0;
  auto u1_ = u1;
  const int nf1 = static_cast<int>(f0x1.extent(3));
  if (stage == 1) {
    Kokkos::deep_copy(DevExeSpace(), u1, u0);
    Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
    if (trans_on) {
      Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
      if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
    }
    c1 = t2_a21*dt;
    c2 = 0.0;
  } else {
    int nmb1 = static_cast<int>(u0.extent(0)) - 1;
    int n3 = static_cast<int>(u0.extent(2)), n2 = static_cast<int>(u0.extent(3));
    int n1 = static_cast<int>(u0.extent(4));
    par_for("vet_t2_avg", DevExeSpace(), 0, nmb1, 0, VET_NVAR-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      u0_(m,n,k,j,i) = heun ? 0.5*(u1_(m,n,k,j,i) + u0_(m,n,k,j,i))
                            : (w1*u1_(m,n,k,j,i) + w0*u0_(m,n,k,j,i));
    });
    auto avgf = [&](DvceArray4D<Real> &a, DvceArray4D<Real> &b) {
      auto a_ = a;
      auto b_ = b;
      par_for("vet_t2_avgf", DevExeSpace(), 0, static_cast<int>(a.extent(0))-1,
              0, static_cast<int>(a.extent(1))-1, 0, static_cast<int>(a.extent(2))-1,
              0, static_cast<int>(a.extent(3))-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        a_(m,k,j,i) = heun ? 0.5*(b_(m,k,j,i) + a_(m,k,j,i))
                           : (w1*b_(m,k,j,i) + w0*a_(m,k,j,i));
      });
    };
    avgf(f0x1, t2f1);
    if (trans_on) {
      avgf(f0x2, t2f2);
      if (trans_x3) {avgf(f0x3, t2f3);}
    }
    c1 = t2_c31*dt;
    c2 = t2_c32*dt;
  }
  (void) nf1;
  auto k1_ = t2k1;
  auto k2_ = t2k2;
  auto in_ = t2inc;
  int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
  int n3 = static_cast<int>(t2inc.extent(2)), n2 = static_cast<int>(t2inc.extent(3));
  int n1 = static_cast<int>(t2inc.extent(4));
  const bool two = (stage != 1);
  par_for("vet_t2_inc", DevExeSpace(), 0, nmb1, 0, VET_T2_NK-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    Real v = c1*k1_(m,n,k,j,i);
    if (two) {v += c2*k2_(m,n,k,j,i);}
    in_(m,n,k,j,i) = v;
  });
  cln_lzc = c1*cln_lzs1 + (two ? c2*cln_lzs2 : 0.0);   // (clean_cs_lz_abs = 2)
  FluxMemoryInc(two, c1, c2, heun ? 0.5 : w1, heun ? 0.5 : w0);
  Time2FaceFct(stage, c1, c2);
  BoundaryEnergyAcc(stage, c1, c2, (stage == 1) ? 0.0 : (heun ? 0.5 : w0), dt);
  t2_solve = (stage == 1) ? VET_T2S_STAGE1 : VET_T2S_STAGE2;
  t2_fail = false;
  nsub = 1;
  dt_sub = g*dt;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Stash
//! \brief stage-loop start (before Heun stage 1): checkpoints the clean energy-flux
//! memory and ray history

void RadiationVet::Time2Stash() {
  FluxMemoryKeep(true);
  RayHistory(true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RayHistory
//! \brief save: keep the rays' history (RayKeep); else restore it (a stage failed)

void RadiationVet::RayHistory(const bool save) {
  DvceArray5D<Real> *a[6] = {&cln_fvi, &cln_fh, &cln_top, &cln_rr, &cln_rpd, &cln_rnew};
  auto &k = cln_rk;
  auto cp = [](auto &d, const auto &s) {
    if (s.size() == 0) {
      d = std::remove_reference_t<decltype(d)>();
      return;
    }
    bool same = d.size() == s.size();
    for (unsigned r = 0; r < s.rank(); ++r) same = same && d.extent(r) == s.extent(r);
    if (!same) {
      d = std::remove_reference_t<decltype(d)>(
          Kokkos::view_alloc(std::string("cln_rk"), Kokkos::WithoutInitializing),
          s.layout());
    }
    Kokkos::deep_copy(d, s);
  };
  // (cln_rfix: no build of the cycle shifts the ring, and each writes the pending Delta
  // before it reads it, so neither is needed back: as the dt_retry checkpoint)
  const bool fix = save ? cln_rfix : k.rfix;
  for (int q = 0; q < 6; ++q) {
    if (fix && (a[q] == &cln_rr || a[q] == &cln_rpd)) {
      k.v[q] = DvceArray5D<Real>();
    } else {
      save ? cp(k.v[q], *a[q]) : cp(*a[q], k.v[q]);
    }
  }
  save ? cp(k.fvd, cln_fvd) : cp(cln_fvd, k.fvd);
  save ? k.fvdm.resize(cln_fvdm.size()) : cln_fvdm.resize(k.fvdm.size());
  for (std::size_t q = 0; q < k.fvdm.size(); ++q) {
    save ? cp(k.fvdm[q], cln_fvdm[q]) : cp(cln_fvdm[q], k.fvdm[q]);
  }
  auto sw = [&](auto &x, auto &y) {if (save) {y = x;} else {x = y;}};
  sw(cln_fvblk, k.fvblk); sw(cln_rb, k.rb); sw(cln_rpb, k.rpb); sw(cln_rbld, k.rbld);
  sw(cln_rpt, k.rpt); sw(cln_rs, k.rs); sw(cln_rguard_sum, k.rg); sw(cln_roff_sum, k.ro);
  sw(cln_rns, k.rns); sw(cln_rok, k.rok); sw(cln_rfix, k.rfix); sw(cln_rlast, k.rlast);
  sw(cln_rjump, k.rjump); sw(cln_fvn, k.fvn); sw(cln_fvkey, k.fvkey);
  sw(cln_rcut, k.rcut); sw(cln_rcut_n, k.rcut_n); sw(cln_rfdone, k.rfdone);
  sw(cln_rhist, k.rhist);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Restore
//! \brief a stage failed: put U^n back (radiation, faces, gas, B, spin) and refresh the
//! hydro; the Driver reruns the Heun stages and takes the step with backward Euler

void RadiationVet::Time2Restore(Driver *pdrive) {
  t2_gt = 0;   // (Time2GasTime)
  emag_hold = false;
  Kokkos::deep_copy(DevExeSpace(), u0, u1);
  Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (fl.on) {
    Kokkos::deep_copy(DevExeSpace(), fl.u0, fl.u1);
    // mhd: the field back to U^n (b1 is the stage-1 copy of b0 incl. ghosts)
    if (fl.mhd) {
      mhd::MHD *pm = pmy_pack->pmhd;
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x1f, pm->b1.x1f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x2f, pm->b1.x2f);
      Kokkos::deep_copy(DevExeSpace(), pm->b0.x3f, pm->b1.x3f);
    }
    FluidRefresh(pdrive);
    if (fl.mhd) EmagBuild(1.0);
  }
  // am_vector: the spin reservoir back to U^n (am_spin1)
  if (am_vec) {
    hydro::AmModule *pam = hydro::AmModuleOf(pmy_pack);
    Kokkos::deep_copy(DevExeSpace(), pam->am_spin, pam->am_spin1);
  }
  FluxMemoryKeep(false);
  RayHistory(false);
  cln_lzp = 0.0;   // (the L_z ledger: the failed stages' pending part, LzBook)
  t2_fail = false;
  t2_ok = false;
  t2_nfall += 1.0;
  // the backward-Euler redo is reported, never silent
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage solve failed, the step is redone "
                "with backward Euler (stage fallback %.0f)\n", pmy_pack->pmesh->ncycle,
                t2_solve, t2_nfall);
    if (cln_flrun > 0.0 && cln_flrid == 4*pmy_pack->pmesh->ncycle + t2_solve) {
      std::printf(
          "CLEAN WARNING cycle %d t2s %d: the implicit_pos_floor creation %.6e of "
                  "the failed attempt is discarded with it\n", pmy_pack->pmesh->ncycle,
                  t2_solve, cln_flrun);
    }
  }
  cln_flrun = 0.0;
  // the predictor increments were stored by the failed solve: the redo starts cold from
  // U^n (extrapolating them again can grow without bound)
  pred_ok = false;
  pred2_ok = false;
  t2_solve = VET_T2S_NONE;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2RotCorrect
//! \brief vet_gd_rotate_slope = correct, on the first step of a new vet_gd rotation block
//! with a valid K1 of the old set: one stage-type solve of dt_c =
//! vet_gd_rotate_slope_frac dt from U^n with the new set and old vector U^n - dt_c K1.
//! Its slope K1 + (I - dt_c A)^{-1} Delta f takes the closure jump where dt_c A is small
//! and keeps K1 where stiff; a failed solve keeps K1.  Returns true when a solve ran.

bool RadiationVet::Time2RotCorrect(Driver *pdrive, Real dt) {
  if (!t2_rotcorr) return false;
  t2_rotcorr = false;
  if (time_scheme != VET_TIME_HESDIRK2 || !Time2Active()) return false;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const Real dtc = (1.0 - 1.0/std::sqrt(2.0))*dt;
  // U^n into the stage-start copies (Time2FormStage(1) overwrites them again)
  Kokkos::deep_copy(DevExeSpace(), u1, u0);
  Kokkos::deep_copy(DevExeSpace(), t2f1, f0x1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), t2f2, f0x2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), t2f3, f0x3);}
  }
  if (fl.on) {Kokkos::deep_copy(DevExeSpace(), fl.u1, fl.u0);}
  // mhd: U^n = (u1, b1) for the solve; b1 still holds the last step's start
  if (fl.on && fl.mhd) {
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x1f, fl.pm->b0.x1f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x2f, fl.pm->b0.x2f);
    Kokkos::deep_copy(DevExeSpace(), fl.pm->b1.x3f, fl.pm->b0.x3f);
  }
  // the old vector U^n - dt_c K1 (every channel)
  {
    auto k1_ = t2k1;
    auto in_ = t2inc;
    const Real c = -dtc;
    int nmb1 = static_cast<int>(t2inc.extent(0)) - 1;
    int n3 = static_cast<int>(t2inc.extent(2)), n2 = static_cast<int>(t2inc.extent(3));
    int n1 = static_cast<int>(t2inc.extent(4));
    par_for("vet_t2_rcinc", DevExeSpace(), 0, nmb1, 0, VET_T2_NK-1, 0, n3-1, 0, n2-1,
            0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      in_(m,n,k,j,i) = c*k1_(m,n,k,j,i);
    });
    cln_lzc = c*cln_lzs1;   // (clean_cs_lz_abs = 2)
  }
  Mesh *pmesh = pmy_pack->pmesh;
  t2_inprobe = true;
  t2_solve = VET_T2S_STAGE1;
  t2_fail = false;
  nsub = 1;
  dt_sub = dtc;
  const Real its0 = impl_itsum;
  pdrive->ExecuteTaskList(pmesh, "vet_before_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "vet_stagen", 1);
  pdrive->ExecuteTaskList(pmesh, "vet_after_stagen", 1);
  vgd_rot_ncpass += static_cast<double>(impl_itsum - its0);
  t2_inprobe = false;
  const bool ok = !t2_fail;
  if (ok) {
    Kokkos::deep_copy(DevExeSpace(), t2k1, t2k2);
    cln_lzs1 = cln_lzs2;   // (clean_cs_lz_abs = 2)
    if (t2_pfst && t2pf1.extent_int(0) > 0) {
      Kokkos::deep_copy(DevExeSpace(), t2pf1, t2pf2);
      t2pf1_ok = t2pf2_ok;
    }
    vgd_rot_ncorr += 1.0;
  } else {
    vgd_rot_ncfail += 1.0;
  }
  // U^n back
  Kokkos::deep_copy(DevExeSpace(), u0, u1);
  Kokkos::deep_copy(DevExeSpace(), f0x1, t2f1);
  if (trans_on) {
    Kokkos::deep_copy(DevExeSpace(), f0x2, t2f2);
    if (trans_x3) {Kokkos::deep_copy(DevExeSpace(), f0x3, t2f3);}
  }
  if (fl.on) {
    Kokkos::deep_copy(DevExeSpace(), fl.u0, fl.u1);
    FluidRefresh(pdrive);
    if (fl.mhd) EmagBuild(1.0);   // b0 untouched; same faces
  }
  t2_solve = VET_T2S_NONE;
  t2_fail = false;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2AutoBeBegin
//! \brief time2_auto_be: called by the Driver at the cycle start (after RotNewBlock).
//! A skip step is put in the state the stage fallback leaves (Time2Restore): no valid
//! slope (Heun hydro, then the BESTORE solve), cold predictor.  A retried cycle calls it
//! again with the same state.

bool RadiationVet::Time2AutoBeBegin() {
  t2_abe_now = false;
  // t2_beonly: every step (the first too) is a skip step
  const bool bo = t2_beonly && (time_scheme == VET_TIME_HESDIRK2);
  if (!bo && (t2_abe <= 0 || time_scheme != VET_TIME_HESDIRK2 || t2_abe_left <= 0)) {
    return false;
  }
  if (!bo && !Time2Active()) return false;   // no slope anyway (the BE path follows)
  t2_abe_now = true;
  t2_ok = false;
  t2_fail = false;
  t2_solve = VET_T2S_NONE;
  pred_ok = false;
  pred2_ok = false;
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2AutoBeEnd
//! \brief an accepted cycle of a hesdirk2 run.  kind 0 stage step, 1 stage fallback +
//! BE redo, 2 time2_auto_be step, 3 other BE step, 4 discarded dt_retry attempt (timing
//! only); tw_a / tw_b: fenced wall seconds of the attempt or step, and of the redo

void RadiationVet::Time2AutoBeEnd(int kind, Real tw_a, Real tw_b) {
  if (kind == 0) {t2w_clean += tw_a;}
  if (kind == 1) {t2w_fail += tw_a; t2w_redo += tw_b;}
  if (kind == 2) {t2w_abe += tw_a;}
  if (kind == 3) {t2w_be += tw_a;}
  if (kind == 4) {t2w_retry += tw_a; return;}
  const bool was_skip = t2_abe_now;
  t2_abe_now = false;
  if (t2_abe <= 0) return;
  if (kind == 0) {
    if (t2_abe_len > 0) {t2_abe_ntry += 1.0; t2_abe_nok += 1.0;}
    t2_abe_len = 0;
    t2_abe_left = 0;
  } else if (kind == 1) {
    if (t2_abe_len > 0) {t2_abe_ntry += 1.0;}
    t2_abe_len = (t2_abe_len == 0) ? 1 : std::min(2*t2_abe_len, t2_abe);
    t2_abe_left = t2_abe_len;
  } else if (kind == 2 && was_skip) {
    t2_abe_n += 1.0;
    t2_abe_left -= 1;
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VsrcNoes
//! \brief vet_source_noesrc in force for the current solve: every backward-Euler solve
//! (under hesdirk2 too), never a stage solve

bool RadiationVet::VsrcNoes() const {
  if (!vsrc_noes) return false;
  if (time_scheme == VET_TIME_BE) return true;
  return (t2_solve != VET_T2S_STAGE1) && (t2_solve != VET_T2S_STAGE2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RunStage
//! \brief one implicit stage: the vet_* task lists (rad_vet_tasks.cpp)

void RadiationVet::RunStage(Driver *pdrive, const int stage) {
  Mesh *pm = pmy_pack->pmesh;
  pdrive->ExecuteTaskList(pm, "vet_before_stagen", stage);
  pdrive->ExecuteTaskList(pm, "vet_stagen", stage);
  pdrive->ExecuteTaskList(pm, "vet_after_stagen", stage);
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::Time2Stage
//! \brief after the hydro's stage `stage`: the implicit radiation solve of the stage;
//! false: it failed (the step goes back to U^n and through BackwardEuler)

bool RadiationVet::Time2Stage(Driver *pdrive, const int stage) {
  if (stage == 1) {Time2GasTime(pdrive, 1);}   // (sdirk2 only)
  Time2FormStage(stage, pmy_pack->pmesh->dt);
  RunStage(pdrive, 1);
  if (t2_fail) return false;
  if (stage == 1) {
    Time2GasTime(pdrive, 3);
    EbndStage1();
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::Time2Accept
//! \brief both stages accepted: the end of the step's radiation

void RadiationVet::Time2Accept(Driver *pdrive) {
  // hesdirk2: U^{n+1} = Y3 is in place and K1 is stored by the stage-2 solve
  t2_solve = radvet::VET_T2S_NONE;
  t2_nstep += 1.0;
  (void) StageTimers(pdrive, 1);
  (void) NewTimeStep(pdrive, 1);
  t2_dtprev = pmy_pack->pmesh->dt;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::BackwardEuler
//! \brief the step's radiation by backward Euler after the hydro stages: no valid slope
//! (first step), a failed stage, or no hesdirk2; under hesdirk2 it stores K1

void RadiationVet::BackwardEuler(Driver *pdrive) {
  // hesdirk2 without a valid slope: backward Euler, which stores K1
  const bool vetbe2 = (time_scheme == VET_TIME_HESDIRK2);
  t2_solve = vetbe2 ? VET_T2S_BESTORE : VET_T2S_NONE;
  int nsub = SetSubsteps(pmy_pack->pmesh->dt);
  // (history: boundary energy)
  BoundaryEnergyAcc(0, 0.0, 0.0, 0.0, pmy_pack->pmesh->dt);
  for (int nst=0; nst<nsub; ++nst) {
    for (int stage=1; stage<=(nstage); ++stage) {
      RunStage(pdrive, stage);
    }
    (void) StageTimers(pdrive, 1);
  }
  (void) NewTimeStep(pdrive, 1);
  if (vetbe2) {
    t2_solve = VET_T2S_NONE;
    t2_ok = true;
    t2_nbe += 1.0;
    t2_dtprev = pmy_pack->pmesh->dt;
  }
}

}  // namespace radvet
