//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file he_drive.cpp
//! \brief hg262-drive-1113: the driven inner boundary of he_star_m1 (see he_drive.hpp
//! for the keys): the inflow table, the live monopole gravity and the RK-exact boundary
//! ledgers.  State lives in this file; the pgen calls the functions of he_drive.hpp.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "utils/repro_sum.hpp"
#include "utils/time_table.hpp"
#include "pgen/he_drive.hpp"

namespace hsdrive {
namespace {

constexpr char kTag[8] = {'H', 'S', 'D', 'R', 'V', '0', '0', '1'};
// the ledgers: 0 mass in (inner face, inward), 1 mass out (outer face, outward),
// 2 energy in, 3 energy out (hydro total energy incl. rho Phi), 4 sponge (sum of the
// stage changes of IEN V, negative)
constexpr int kNL = 5;

bool on_ = false, inf_ = false, live_ = false;
timetable::Table itab_, mtab_;
bool mtab_on_ = false;
Real t0_ = 0.0;
SetupArgs a_;
int every_ = 1, logev_ = 1;
char logfn_[1024] = {0};   // <basename>.hsdrive.txt
// RK stage coefficients of the integrator (driver.cpp), for the stage ledger
int nst_ = 0;
Real gam0_[4] = {0.0, 0.0, 0.0, 0.0}, gam1_[4] = {0.0, 0.0, 0.0, 0.0};
Real beta_[4] = {0.0, 0.0, 0.0, 0.0};
Real led_[kNL] = {0.0, 0.0, 0.0, 0.0, 0.0}, ledn_[kNL] = {0.0, 0.0, 0.0, 0.0, 0.0};
int led_cyc_ = -1;
Real erad_ = 0.0, wgrav_ = 0.0, tprev_ = 0.0, ain_ = 0.0;
// open items 11-13: the radiation energy through the physical x1 faces AS THE SCHEME
// APPLIES IT (RadiationM1 Ebnd with ebnd_split: the hesdirk2 tableau weights of the
// slopes' face fluxes, incl. implicit_bc_advect's end-face A E): in at x1min, out at
// x1max (the Marshak top).  ebd_ = false: the transport cannot book it (columns 0).
bool ebd_ = false;
// problem/he_drive_closure_lin = const | table: the L_in the adaptive MLT closure targets
bool cltab_ = false;
Mesh *pmesh_ = nullptr;
// the BC time: t_n + dt inside a cycle (set by StageLedger), else the mesh time
bool bct_set_ = false;
Real bct_ = 0.0;
// live gravity: reference shell-mean densities (active x1 cells), the table mass at the
// reference, and the APPLIED increments (cell shell-mass changes dM_i, interior dM_int)
bool ref_ = false;
Real mref_ = 0.0, dmint_ = 0.0;
std::vector<Real> rhoref_, ddm_;
// a restart's stored numbers (RestartStrip -> RestartApply)
bool rst_have_ = false;
std::vector<Real> rst_;
int nx1_ = 0;

void Fatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in he_drive.cpp (hg262-drive-1113)" << std::endl << msg
            << std::endl;
  std::exit(EXIT_FAILURE);
}

//! rank-ordered sum of v[0..n) over all ranks (deterministic for a given layout)
void DetSum(Real *v, const int n) {
#if MPI_PARALLEL_ENABLED
  const int nr = global_variable::nranks;
  if (nr <= 1 || n <= 0) return;
  std::vector<Real> mine(v, v + n);
  std::vector<Real> all(static_cast<size_t>(n)*nr);
  MPI_Allgather(mine.data(), n, MPI_ATHENA_REAL, all.data(), n, MPI_ATHENA_REAL,
                MPI_COMM_WORLD);
  for (int i = 0; i < n; ++i) {
    Real s = 0.0;
    for (int r = 0; r < nr; ++r) {s += all[static_cast<size_t>(r)*n + i];}
    v[i] = s;
  }
#else
  (void) v; (void) n;
#endif
}

//! shell sums over the active cells, per x1 index (absolute, 0 .. ncells1-1):
//! q = 0 rho V, 1 V, 2 IEN V, 3 E_rad V.  Order-independent under layout_invariant
//! (repro::ShellSums), else a fixed serial order per rank + DetSum.
void ShellSums(Mesh *pm, std::vector<Real> *out) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nc1 = indcs.nx1 + 2*indcs.ng;
  const int nmb = pmbp->nmb_thispack;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto u0 = fl.u0;
  auto ur = pmbp->pradm1->u0;
  auto &vol = pmbp->pcoord->volume;
  constexpr int nq = 4;
  out->assign(static_cast<size_t>(nq)*nc1, 0.0);
  if (pm->layout_invariant) {
    std::vector<Real> rs;
    repro::ShellSums("hsdrv_shell_r", nq, nc1, nmb, ks, ke, js, je, is, ie,
                     KOKKOS_LAMBDA(const int q, const int m, const int k, const int j,
                                   const int i) {
      const Real v = vol(m,k,j,i);
      if (q == 0) {return u0(m,IDN,k,j,i)*v;}
      if (q == 1) {return v;}
      if (q == 2) {return u0(m,IEN,k,j,i)*v;}
      return ur(m,radm1::M1_E,k,j,i)*v;
    }, rs);
    for (int q = 0; q < nq; ++q) {
      for (int i = is; i <= ie; ++i) {
        (*out)[static_cast<size_t>(q)*nc1 + i] = rs[static_cast<size_t>(q)*nc1 + i];
      }
    }
    return;
  }
  DvceArray2D<Real> sh("hsdrv_shell", nq, nc1);
  par_for("hsdrv_shell", DevExeSpace(), is, ie, KOKKOS_LAMBDA(const int i) {
    Real s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (int m = 0; m < nmb; ++m) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          const Real v = vol(m,k,j,i);
          s0 += u0(m,IDN,k,j,i)*v;
          s1 += v;
          s2 += u0(m,IEN,k,j,i)*v;
          s3 += ur(m,radm1::M1_E,k,j,i)*v;
        }
      }
    }
    sh(0,i) = s0; sh(1,i) = s1; sh(2,i) = s2; sh(3,i) = s3;
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sh);
  for (int q = 0; q < nq; ++q) {
    for (int i = 0; i < nc1; ++i) {(*out)[static_cast<size_t>(q)*nc1 + i] = h(q,i);}
  }
  DetSum(out->data(), nq*nc1);
}

//! the radial grid of the (single) x1 MeshBlock row: faces 0 .. nc1, centres 0 .. nc1-1
void RadialGrid(Mesh *pm, std::vector<Real> *rf, std::vector<Real> *rc) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const int nc1 = pm->mb_indcs.nx1 + 2*pm->mb_indcs.ng;
  auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  rf->resize(nc1 + 1);
  rc->resize(nc1);
  for (int i = 0; i <= nc1; ++i) {(*rf)[i] = hf(0,i);}
  for (int i = 0; i < nc1; ++i) {(*rc)[i] = hc(0,i);}
}

//! dPhi on the cells and x1 faces (absolute indices, ghosts included) of the mass
//! increments dmi (inside r_in) and ddm[i - is] (shell i, uniform density in the cell);
//! zero at the face r_in = x1f(is), like every potential of this pgen
void DPhi(Mesh *pm, const Real dmi, const std::vector<Real> &ddm, std::vector<Real> *pc,
          std::vector<Real> *pf) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int nc1 = indcs.nx1 + 2*indcs.ng;
  std::vector<Real> rf, rc;
  RadialGrid(pm, &rf, &rc);
  const Real g = a_.gcode;
  pc->assign(nc1, 0.0);
  pf->assign(nc1 + 1, 0.0);
  const Real rin = rf[is];
  // below r_in: the point mass dM_int
  for (int i = 0; i < is; ++i) {
    (*pf)[i] = g*dmi*(1.0/rin - 1.0/rf[i]);
    (*pc)[i] = g*dmi*(1.0/rin - 1.0/rc[i]);
  }
  Real dm = dmi, ph = 0.0;
  for (int i = is; i <= ie; ++i) {
    const Real r0 = rf[i], r1 = rf[i+1], d3 = r1*r1*r1 - r0*r0*r0;
    const Real b = ddm[i - is]/d3, a = dm - b*r0*r0*r0;
    (*pf)[i] = ph;
    const Real x = rc[i];
    (*pc)[i] = ph + g*(a*(1.0/r0 - 1.0/x) + 0.5*b*(x*x - r0*r0));
    ph += g*(a*(1.0/r0 - 1.0/r1) + 0.5*b*(r1*r1 - r0*r0));
    dm += ddm[i - is];
  }
  (*pf)[ie+1] = ph;
  const Real rt = rf[ie+1];
  for (int i = ie + 1; i < nc1; ++i) {
    (*pc)[i] = ph + g*dm*(1.0/rt - 1.0/rc[i]);
    (*pf)[i+1] = ph + g*dm*(1.0/rt - 1.0/rf[i+1]);
  }
}

//! the reference potentials (device clones of every potential array when the live
//! gravity reference is taken; on a restart the rebuilt arrays ARE the reference): the
//! arrays are always SET to reference + dPhi (never accumulated), so a restart that
//! rebuilds the reference and applies the stored dPhi reproduces them bitwise
DvceArray4D<Real> rcc_, rx1_, rx2_, rx3_, rwc_, rwf_;

void TakeRef(Mesh *pm) {
  radm1::FluidRef fl = radm1::FluidRef::Get(pm->pmb_pack);
  auto clone = [](const DvceArray4D<Real> &a, DvceArray4D<Real> *c) {
    *c = DvceArray4D<Real>(Kokkos::view_alloc(std::string("hsdrv_ref"),
                                              Kokkos::WithoutInitializing), a.layout());
    Kokkos::deep_copy(*c, a);
  };
  clone(fl.phicc0, &rcc_);
  clone(fl.phi0->x1f, &rx1_);
  clone(fl.phi0->x2f, &rx2_);
  clone(fl.phi0->x3f, &rx3_);
  if (fl.phicc_wb.data() != fl.phicc0.data()) {clone(fl.phicc_wb, &rwc_);}
  if (fl.phi_wb_x1f.data() != fl.phi0->x1f.data()) {clone(fl.phi_wb_x1f, &rwf_);}
}

//! set every potential array to reference + (dpc, dpf); energy: E += rho (new - old)
void SetPhi(Mesh *pm, const std::vector<Real> &dpc, const std::vector<Real> &dpf,
            const bool energy) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int nc1 = indcs.nx1 + 2*indcs.ng;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng - 1) : 0;
  const int n1m1 = nc1 - 1;
  const int nmb1 = pmbp->nmb_thispack - 1;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  DvceArray1D<Real> dc("hsdrv_dpc", nc1), df("hsdrv_dpf", nc1 + 1);
  auto hdc = Kokkos::create_mirror_view(dc);
  auto hdf = Kokkos::create_mirror_view(df);
  for (int i = 0; i < nc1; ++i) {hdc(i) = dpc[i];}
  for (int i = 0; i <= nc1; ++i) {hdf(i) = dpf[i];}
  Kokkos::deep_copy(dc, hdc);
  Kokkos::deep_copy(df, hdf);
  auto phicc = fl.phicc0;
  auto ph1 = fl.phi0->x1f, ph2 = fl.phi0->x2f, ph3 = fl.phi0->x3f;
  auto pwc = fl.phicc_wb, pwf = fl.phi_wb_x1f;
  // the effective pair is its own allocation only after EnableWBEffectivePotential
  const bool wbs = (pwc.data() != phicc.data());
  const bool wfs = (pwf.data() != ph1.data());
  auto rcc = rcc_, rx1 = rx1_, rx2 = rx2_, rx3 = rx3_, rwc = rwc_, rwf = rwf_;
  auto u0 = fl.u0;
  par_for("hsdrv_setphi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = dc(i);
    const Real pn = rcc(m,k,j,i) + d;
    if (energy) {u0(m,IEN,k,j,i) += u0(m,IDN,k,j,i)*(pn - phicc(m,k,j,i));}
    phicc(m,k,j,i) = pn;
    ph1(m,k,j,i) = rx1(m,k,j,i) + df(i);
    if (i == n1m1) {ph1(m,k,j,i+1) = rx1(m,k,j,i+1) + df(i+1);}
    ph2(m,k,j,i) = rx2(m,k,j,i) + d;
    ph3(m,k,j,i) = rx3(m,k,j,i) + d;
    if (j == n2m1) {ph2(m,k,j+1,i) = rx2(m,k,j+1,i) + d;}
    if (k == n3m1) {ph3(m,k+1,j,i) = rx3(m,k+1,j,i) + d;}
    if (wbs) {pwc(m,k,j,i) = rwc(m,k,j,i) + d;}
    if (wfs) {
      pwf(m,k,j,i) = rwf(m,k,j,i) + df(i);
      if (i == n1m1) {pwf(m,k,j,i+1) = rwf(m,k,j,i+1) + df(i+1);}
    }
  });
  if (fl.ph != nullptr) {fl.ph->wb_cache_built = false;}
  if (fl.pm != nullptr) {fl.pm->wb_cache_built = false;}
}

Real MtabAt(const Real t) {return mtab_on_ ? mtab_.At(0, t - t0_) : 0.0;}

//! the live-gravity refresh at the current state (energy: E += rho dPhi_increment).
//! sh = the shell sums of this state (ShellSums).  Returns the work sum rho d(dPhi) V.
Real Refresh(Mesh *pm, const std::vector<Real> &sh, const bool energy) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int nc1 = indcs.nx1 + 2*indcs.ng;
  std::vector<Real> rf, rc;
  RadialGrid(pm, &rf, &rc);
  if (!ref_) {
    rhoref_.assign(nx1_, 0.0);
    for (int i = is; i <= ie; ++i) {rhoref_[i - is] = sh[i]/sh[nc1 + i];}
    mref_ = MtabAt(pm->time);
    dmint_ = 0.0;
    ddm_.assign(nx1_, 0.0);
    ref_ = true;
    TakeRef(pm);
    if (global_variable::my_rank == 0) {
      std::cout << "he_drive: live gravity reference taken at cycle " << pm->ncycle
                << ", t = " << pm->time << " (M_int table " << mref_ << ")" << std::endl;
    }
    return 0.0;
  }
  std::vector<Real> nddm(nx1_);
  for (int i = is; i <= ie; ++i) {
    const Real r0 = rf[i], r1 = rf[i+1];
    const Real drho = sh[i]/sh[nc1 + i] - rhoref_[i - is];
    nddm[i - is] = (4.0*M_PI/3.0)*drho*(r1*r1*r1 - r0*r0*r0);
  }
  const Real ndmi = MtabAt(pm->time) - mref_;
  std::vector<Real> pc0, pf0, pc1, pf1;
  DPhi(pm, dmint_, ddm_, &pc0, &pf0);
  DPhi(pm, ndmi, nddm, &pc1, &pf1);
  std::vector<Real> dc(nc1);
  Real w = 0.0;
  for (int i = 0; i < nc1; ++i) {dc[i] = pc1[i] - pc0[i];}
  for (int i = is; i <= ie; ++i) {w += dc[i]*sh[i];}
  SetPhi(pm, pc1, pf1, energy);
  ddm_ = nddm;
  dmint_ = ndmi;
  return w;
}

//! the four boundary-face sums of this stage (mass in at is, mass out at ie+1, energy
//! in, energy out), with the effective fluxes of the pgen's wall corrections (hs_zflux)
void FaceSums(Mesh *pm, const bool binf, const bool zflux, Real *f) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb = pmbp->nmb_thispack;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto flx1 = fl.uflx->x1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  array_sum::GlobalSum s;
  Kokkos::parallel_reduce("hsdrv_face", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nx3*nx2),
  KOKKOS_LAMBDA(const int n, array_sum::GlobalSum &acc) {
    const int m = n/(nx3*nx2);
    const int k = (n/nx2)%nx3 + ks, j = n%nx2 + js;
    array_sum::GlobalSum h;
    for (int q = 0; q < NREDUCTION_VARIABLES; ++q) {h.the_array[q] = 0.0;}
    if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      // the wall correction removes the inner mass and energy fluxes (wall mode)
      if (binf || !zflux) {
        h.the_array[0] = flx1(m,IDN,k,j,is)*area1(m,k,j,is);
        h.the_array[2] = flx1(m,IEN,k,j,is)*area1(m,k,j,is);
      }
    } else if (mbbcs.d_view(m, BoundaryFace::inner_x1) != BoundaryFlag::block) {
      h.the_array[0] = flx1(m,IDN,k,j,is)*area1(m,k,j,is);
      h.the_array[2] = flx1(m,IEN,k,j,is)*area1(m,k,j,is);
    }
    if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      // the outer correction removes inflow (flux <= 0)
      if (!zflux || flx1(m,IDN,k,j,ie+1) > 0.0) {
        h.the_array[1] = flx1(m,IDN,k,j,ie+1)*area1(m,k,j,ie+1);
        h.the_array[3] = flx1(m,IEN,k,j,ie+1)*area1(m,k,j,ie+1);
      }
    } else if (mbbcs.d_view(m, BoundaryFace::outer_x1) != BoundaryFlag::block) {
      h.the_array[1] = flx1(m,IDN,k,j,ie+1)*area1(m,k,j,ie+1);
      h.the_array[3] = flx1(m,IEN,k,j,ie+1)*area1(m,k,j,ie+1);
    }
    acc += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(s));
  for (int q = 0; q < 4; ++q) {f[q] = s.the_array[q];}
#if MPI_PARALLEL_ENABLED
  Real g[4];
  MPI_Allreduce(f, g, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < 4; ++q) {f[q] = g[q];}
#endif
}

void WriteLog(Mesh *pm, const std::vector<Real> &sh) {
  if (global_variable::my_rank != 0) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int nc1 = indcs.nx1 + 2*indcs.ng;
  Real mt = 0.0, et = 0.0;
  for (int i = is; i <= ie; ++i) {
    mt += sh[i];
    et += sh[2*nc1 + i] + sh[3*nc1 + i];
  }
  FILE *fp = std::fopen(logfn_, "a");
  if (fp == nullptr) return;
  const Real gmi = a_.gm + a_.gcode*dmint_;
  // the applied potential increment and mass increment at the top face (live gravity)
  Real ptop = 0.0, gdm = 0.0;
  if (ref_) {
    std::vector<Real> pc, pf;
    DPhi(pm, dmint_, ddm_, &pc, &pf);
    ptop = pf[ie+1];
    Real dm = dmint_;
    for (int i = 0; i < nx1_; ++i) {dm += ddm_[i];}
    gdm = a_.gcode*dm;
  }
  radm1::RadiationM1 *pm1 = pm->pmb_pack->pradm1;
  const Real ebi = ebd_ ? pm1->ebnd_acc_in : 0.0;
  const Real ebo = ebd_ ? pm1->ebnd_acc_out : 0.0;
  const Real esa = ebd_ ? pm1->ebnd_acc_src : 0.0;
  std::fprintf(fp, "%d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e "
               "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", pm->ncycle, pm->time,
               mt, et, led_[0], led_[1], led_[2], led_[3], led_[4], erad_, wgrav_, gmi,
               InnerFlux(pm, a_.fconst), ain_, ptop, gdm, ebi, ebo, esa);
  std::fclose(fp);
}

}  // namespace

void CycleBody(Mesh *pm, const bool log);

bool On() {return on_;}
bool InflowOn() {return inf_;}
bool LiveOn() {return live_;}

void Setup(ParameterInput *pin, Mesh *pm, const bool restart, const SetupArgs &a) {
  a_ = a;
  const bool hinf = pin->DoesParameterExist("problem","he_bc_inner_table");
  const bool hlive = pin->DoesParameterExist("problem","he_gm_live") &&
                     pin->GetBoolean("problem","he_gm_live");
  radm1::RadiationM1 *pm1 = pm->pmb_pack->pradm1;
  const bool hftab = (pm1 != nullptr) && pm1->ftab_on;
  bool led = hinf || hlive || hftab;
  if (pin->DoesParameterExist("problem","he_drive_ledger")) {
    led = pin->GetBoolean("problem","he_drive_ledger");
    if (!led && (hinf || hlive)) {
      Fatal("problem/he_drive_ledger = false is not allowed with he_bc_inner_table or "
            "he_gm_live (their state travels with the ledger)");
    }
  }
  // open items 11-13 (item 2): the adaptive MLT closure targets L_in.  const (default):
  // the deck's constant implicit_flux_x1min (the relaxed state's target); table: the
  // current table flux (the first version; it re-targets the closure flux at EVERY
  // radius to L_in(t) at once, which takes the extra injection out of the base cells
  // convectively and drops the radiative L_bot, see DRIVE_CODE_1113.md)
  {
    const bool cln = pin->DoesParameterExist("problem","he_drive_closure_lin");
    const std::string cl = (hftab && cln) ?
        pin->GetString("problem","he_drive_closure_lin") : std::string("const");
    if (cl != "const" && cl != "table") {
      Fatal("problem/he_drive_closure_lin must be const or table");
    }
    cltab_ = (cl == "table");
  }
  on_ = led;
  pmesh_ = pm;
  if (!on_) {
    if (rst_have_ && global_variable::my_rank == 0) {
      std::cout << "he_drive: the restart file carries a drive trailer but no drive key "
                << "is named: ignored" << std::endl;
    }
    return;
  }
  auto &indcs = pm->mb_indcs;
  nx1_ = indcs.nx1;
  if (pm->mesh_indcs.nx1 != indcs.nx1) {
    Fatal("the drive needs ONE MeshBlock along x1 (mesh nx1 = meshblock nx1)");
  }
  t0_ = pin->GetOrAddReal("problem","he_drive_t0",0.0);
  every_ = pin->GetOrAddInteger("problem","he_gm_live_every",1);
  logev_ = pin->GetOrAddInteger("problem","he_drive_log_every",1);
  if (every_ < 1 || logev_ < 1) {
    Fatal("problem/he_gm_live_every and he_drive_log_every must be >= 1");
  }
  std::snprintf(logfn_, sizeof(logfn_), "%s.hsdrive.txt",
                pin->GetString("job","basename").c_str());
  // the integrator's stage coefficients (driver.cpp), for the stage ledger
  const std::string ig = pin->GetString("time","integrator");
  if (ig == "rk1") {
    nst_ = 1; gam0_[0] = 0.0; gam1_[0] = 1.0; beta_[0] = 1.0;
  } else if (ig == "rk2") {
    nst_ = 2; gam0_[0] = 0.0; gam1_[0] = 1.0; beta_[0] = 1.0;
    gam0_[1] = 0.5; gam1_[1] = 0.5; beta_[1] = 0.5;
    if (pin->DoesBlockExist("rad_m1") &&
        pin->DoesParameterExist("rad_m1", "time2_tableau") &&
        pin->GetString("rad_m1", "time2_tableau").compare("trbdf2") == 0) {
      const Real c2 = 2.0 - std::sqrt(2.0);
      beta_[0] = c2;
      gam0_[1] = 0.25;
      gam1_[1] = 0.75;
      beta_[1] = 0.5/c2;
    }
  } else if (ig == "rk3") {
    nst_ = 3; gam0_[0] = 0.0; gam1_[0] = 1.0; beta_[0] = 1.0;
    gam0_[1] = 0.25; gam1_[1] = 0.75; beta_[1] = 0.25;
    gam0_[2] = 2.0/3.0; gam1_[2] = 1.0/3.0; beta_[2] = 2.0/3.0;
  } else {
    Fatal("the drive ledger supports <time>/integrator = rk1 | rk2 | rk3");
  }
  if (hinf) {
    if (!a.binf) Fatal("problem/he_bc_inner_table needs problem/he_bc_inner = inflow");
    const std::string fn = pin->GetString("problem","he_bc_inner_table");
    const std::string err = itab_.Read(fn, 3);
    if (!err.empty()) Fatal("problem/he_bc_inner_table: " + err);
    if (!(a.rho_rin > 0.0 && a.eint_rin > 0.0)) {
      Fatal("problem/he_bc_inner_table: the initial column at r_in is not positive");
    }
    inf_ = true;
  }
  if (hlive) {
    if (a.pointmass) Fatal("problem/he_gm_live needs he_gravity = wb");
    if (a.bc_hse) Fatal("problem/he_gm_live does not support he_bc_outer = hse");
    if (!(a.gcode > 0.0)) Fatal("problem/he_gm_live needs <units> (G in code units)");
    if (pin->DoesParameterExist("problem","he_gm_table")) {
      const std::string fn = pin->GetString("problem","he_gm_table");
      const std::string err = mtab_.Read(fn, 1);
      if (!err.empty()) Fatal("problem/he_gm_table: " + err);
      mtab_on_ = true;
    }
    live_ = true;
  }
  // open items 11-13: the scheme's own x1-face radiation energy (Ebnd, split), so the
  // Marshak top loss is in the ledger; only where Ebnd is defined (implicit transport,
  // no multi-rate, no coupling_split)
  ebd_ = (pm1 != nullptr) && (pm1->transport == radm1::M1_TRANSPORT_IMPLICIT) &&
         !pm1->MRActive() && (pm1->csplit == 0);
  if (ebd_) {
    pm1->ebnd_on = true;
    pm1->ebnd_adv = true;   // EbndFluxX1 adds A E only where implicit_bc_advect is on
    pm1->ebnd_split = true;
  }
  // the inner face area of the whole mesh (Erad ledger)
  {
    MeshBlockPack *pmbp = pm->pmb_pack;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nx2 = indcs.nx2, nx3 = indcs.nx3, nmb = pmbp->nmb_thispack;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &mbbcs = pmbp->pmb->mb_bcs;
    Real s = 0.0;
    Kokkos::parallel_reduce("hsdrv_ain", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            nmb*nx3*nx2),
    KOKKOS_LAMBDA(const int n, Real &acc) {
      const int m = n/(nx3*nx2);
      const int k = (n/nx2)%nx3 + ks, j = n%nx2 + js;
      if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        acc += area1(m,k,j,is);
      }
    }, Kokkos::Sum<Real>(s));
    DetSum(&s, 1);
    ain_ = s;
  }
  tprev_ = pm->time;
  if (global_variable::my_rank == 0) {
    std::cout << "he_drive (hg262-drive-1113): ledger on (log " << logfn_ << " every "
              << logev_ << " cycles), inflow table " << (inf_ ? "on" : "off")
              << ", live gravity " << (live_ ? "on" : "off");
    if (live_) {
      std::cout << " (every " << every_ << " cycles, M_int table "
                << (mtab_on_ ? "on" : "off") << ")";
    }
    std::cout << ", rad_m1 flux table " << (hftab ? "on" : "off") << ", A_in " << ain_
              << ", table t0 " << t0_ << ", Ebnd faces " << (ebd_ ? "on" : "off")
              << ", closure L_in " << (cltab_ ? "table" : "const")
              << (restart ? " (restart)" : "") << std::endl;
  }
  if (!restart && global_variable::my_rank == 0) {
    FILE *fp = std::fopen(logfn_, "w");
    if (fp != nullptr) {
      std::fprintf(fp, "# hg262-drive-1113 ledger: [1]cycle [2]time [3]M_tot [4]E_tot "
                   "(IEN + E_rad, active) [5]Mb_in [6]Mb_out [7]Eb_in [8]Eb_out "
                   "[9]E_sponge [10]Erad_in (A_in int F dt) [11]W_grav [12]GM_int "
                   "[13]F_in [14]A_in [15]dPhi(x1f(ie+1)) [16]G dm(x1f(ie+1)) "
                   "[17]ErB_in [18]ErB_out (radiation energy in at x1min / OUT at x1max "
                   "as the scheme applies it, RadiationM1 Ebnd) [19]E_src (the explicit "
                   "radiation source esrc = the MLT scaffold, as the solves apply it). "
                   "Budget: [4]-[4]_0 = "
                   "[7]-[8]+[9]+[11]+[17]-[18]+[19]\n");
      std::fclose(fp);
    }
  }
}

void RestartStrip(std::vector<char> *blk) {
  const std::size_t nt = sizeof(kTag) + sizeof(Real);
  if (blk->size() < nt ||
      std::memcmp(blk->data() + blk->size() - sizeof(kTag), &(kTag[0]), sizeof(kTag))
      != 0) {
    return;
  }
  Real rn;
  std::memcpy(&rn, blk->data() + blk->size() - nt, sizeof(Real));
  const std::size_t n = static_cast<std::size_t>(rn);
  const std::size_t nb = nt + n*sizeof(Real);
  if (blk->size() < nb) Fatal("the restart file's HSDRV001 trailer is truncated");
  rst_.assign(n, 0.0);
  std::memcpy(rst_.data(), blk->data() + blk->size() - nb, n*sizeof(Real));
  blk->resize(blk->size() - nb);
  rst_have_ = true;
}

void RstAppend(std::vector<char> *out) {
  if (!on_) return;
  // payload: nx1, ledgers (kNL), ledgers at the stage-loop start (kNL), their cycle,
  // erad, wgrav, tprev, ref flag, mref, dmint, rhoref (nx1), ddm (nx1), [open items
  // 11-13] ebnd_acc_in, ebnd_acc_out, ebnd_acc_src, a spare
  std::vector<Real> p;
  p.push_back(static_cast<Real>(nx1_));
  for (int q = 0; q < kNL; ++q) {p.push_back(led_[q]);}
  for (int q = 0; q < kNL; ++q) {p.push_back(ledn_[q]);}
  p.push_back(static_cast<Real>(led_cyc_));
  p.push_back(erad_);
  p.push_back(wgrav_);
  p.push_back(tprev_);
  p.push_back(ref_ ? 1.0 : 0.0);
  p.push_back(mref_);
  p.push_back(dmint_);
  for (int i = 0; i < nx1_; ++i) {p.push_back(ref_ ? rhoref_[i] : 0.0);}
  for (int i = 0; i < nx1_; ++i) {p.push_back(ref_ ? ddm_[i] : 0.0);}
  radm1::RadiationM1 *pm1 = (pmesh_ != nullptr) ? pmesh_->pmb_pack->pradm1 : nullptr;
  p.push_back(ebd_ ? pm1->ebnd_acc_in : 0.0);
  p.push_back(ebd_ ? pm1->ebnd_acc_out : 0.0);
  p.push_back(ebd_ ? pm1->ebnd_acc_src : 0.0);
  p.push_back(0.0);   // spare
  const Real rn = static_cast<Real>(p.size());
  const std::size_t o = out->size();
  out->resize(o + (p.size() + 1)*sizeof(Real) + sizeof(kTag));
  std::memcpy(out->data() + o, p.data(), p.size()*sizeof(Real));
  std::memcpy(out->data() + o + p.size()*sizeof(Real), &rn, sizeof(Real));
  std::memcpy(out->data() + o + (p.size() + 1)*sizeof(Real), &(kTag[0]), sizeof(kTag));
}

void RestartApply(Mesh *pm) {
  if (!on_) return;
  bool have = false;
  if (rst_have_) {
    const std::size_t nexp = 1 + 2*kNL + 7 + 2*static_cast<std::size_t>(nx1_);
    // the first version's trailer (caf4fe6bb) has no Ebnd pair: those start at 0
    if ((rst_.size() != nexp && rst_.size() != nexp + 4) ||
        static_cast<int>(rst_[0]) != nx1_) {
      Fatal("the restart file's HSDRV001 trailer does not match this radial grid");
    }
    std::size_t o = 1;
    for (int q = 0; q < kNL; ++q) {led_[q] = rst_[o++];}
    for (int q = 0; q < kNL; ++q) {ledn_[q] = rst_[o++];}
    led_cyc_ = static_cast<int>(rst_[o++]);
    erad_ = rst_[o++];
    wgrav_ = rst_[o++];
    tprev_ = rst_[o++];
    ref_ = (rst_[o++] != 0.0);
    mref_ = rst_[o++];
    dmint_ = rst_[o++];
    rhoref_.assign(rst_.begin() + o, rst_.begin() + o + nx1_);
    o += nx1_;
    ddm_.assign(rst_.begin() + o, rst_.begin() + o + nx1_);
    o += nx1_;
    if (rst_.size() == nexp + 4) {
      if (ebd_) {
        radm1::RadiationM1 *pm1 = pm->pmb_pack->pradm1;
        pm1->ebnd_acc_in = rst_[o];
        pm1->ebnd_acc_out = rst_[o + 1];
        pm1->ebnd_acc = pm1->ebnd_acc_in - pm1->ebnd_acc_out;
        pm1->ebnd_acc_src = rst_[o + 2];
      }
    }
    have = true;
  }
  if (!live_ && have && ref_ &&
      (dmint_ != 0.0 || ddm_ != std::vector<Real>(nx1_, 0.0))) {
    Fatal("the restart state carries a live-gravity increment: keep problem/he_gm_live "
          "= true on this restart");
  }
  // the straight run's end-of-cycle work after it wrote this restart file (Cycle: the
  // Erad ledger of the last cycle and the live-gravity refresh, no log line)
  if (have) {
    if (live_ && ref_) {
      std::vector<Real> pc, pf;
      DPhi(pm, dmint_, ddm_, &pc, &pf);
      TakeRef(pm);   // the rebuilt potentials are the reference
      SetPhi(pm, pc, pf, false);
    }
    CycleBody(pm, false);
    return;
  }
  if (!live_) return;
  std::vector<Real> sh;
  ShellSums(pm, &sh);
  Refresh(pm, sh, false);   // no trailer: the reference is this restart state
}

void InflowAt(Mesh *pm, Real *d, Real *t, Real *v) {
  const Real tt = (bct_set_ ? bct_ : pm->time) - t0_;
  *d = itab_.At(0, tt);
  *t = itab_.At(1, tt)/a_.tk;
  *v = itab_.At(2, tt);
}

void StageLedger(Mesh *pm, const Real bdt, const bool binf, const bool zflux,
                 const Real dsp) {
  int s = -1;
  for (int q = 0; q < nst_; ++q) {
    if (bdt == beta_[q]*pm->dt) {s = q; break;}
  }
  if (s < 0) Fatal("the stage ledger cannot identify the RK stage from beta dt");
  if (s == 0 && pm->ncycle != led_cyc_) {
    for (int q = 0; q < kNL; ++q) {ledn_[q] = led_[q];}
    led_cyc_ = pm->ncycle;
  }
  Real f[4];
  FaceSums(pm, binf, zflux, f);
  const Real g0 = gam0_[s], g1 = gam1_[s];
  for (int q = 0; q < 4; ++q) {led_[q] = g0*led_[q] + g1*ledn_[q] + bdt*f[q];}
  led_[4] = g0*led_[4] + g1*ledn_[4] + dsp;
  bct_set_ = true;
  bct_ = pm->time + pm->dt;
}

//! the end-of-cycle work; log = false on a restart (the straight run wrote that line)
void CycleBody(Mesh *pm, const bool log) {
  radm1::RadiationM1 *pm1 = pm->pmb_pack->pradm1;
  const Real t = pm->time;
  if (pm1 != nullptr && pm1->ftab_on) {
    erad_ += ain_*pm1->InnerFluxIntegral(tprev_, t);
  } else {
    erad_ += ain_*a_.fconst*(t - tprev_);
  }
  tprev_ = t;
  const bool logd = log && (pm->ncycle % logev_ == 0);
  const bool refr = live_ && (pm->ncycle % every_ == 0 || !ref_);
  if (!logd && !refr) return;
  std::vector<Real> sh;
  ShellSums(pm, &sh);
  if (refr) {
    wgrav_ += Refresh(pm, sh, true);
    if (logd) {ShellSums(pm, &sh);}
  }
  if (logd) {WriteLog(pm, sh);}
}

void Cycle(Mesh *pm) {
  if (on_) {CycleBody(pm, true);}
}

Real ClosureFlux(Mesh *pm, const Real fconst) {
  return cltab_ ? InnerFlux(pm, fconst) : fconst;
}

Real InnerFlux(Mesh *pm, const Real fconst) {
  radm1::RadiationM1 *pm1 = pm->pmb_pack->pradm1;
  if (pm1 != nullptr && pm1->ftab_on) {return pm1->InnerFluxAt(pm->time);}
  return fconst;
}

void Final() {
  rcc_ = DvceArray4D<Real>(); rx1_ = DvceArray4D<Real>(); rx2_ = DvceArray4D<Real>();
  rx3_ = DvceArray4D<Real>(); rwc_ = DvceArray4D<Real>(); rwf_ = DvceArray4D<Real>();
  itab_ = timetable::Table();
  mtab_ = timetable::Table();
}

}  // namespace hsdrive
