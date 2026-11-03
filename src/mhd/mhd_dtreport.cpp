//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_dtreport.cpp
//! \brief dtdiag-1103: opt-in report of WHAT sets the MHD CFL time step.
//!
//! <mhd>/dt_report = N (default 0 = off; read only when named): every N-th call of
//! MHD::NewTimeStep (last stage, i.e. once per cycle) one extra pass over all cells
//! prints, per rank,
//!   - the fast-magnetosonic dt (the one NewTimeStep uses; the copy here must equal
//!     dtnew, which is printed alongside as a check), over ALL cells and over cells with
//!     p >= dt_report_p1 and p >= dt_report_p2 (code pressure; defaults 1 and 1e3 =
//!     1e-6 and 1e-3 bar in cgs code units);
//!   - the same with the SOUND speed in place of the fast speed (B ignored) and with the
//!     FLOW speed alone, so the signal can be named;
//!   - the number of cells within 1.2x of the minimum in the three pressure bands, and of
//!     cells at the <mhd>/valfven_max_cgs cap (B^2/rho >= 0.99 v_A,max^2);
//!   - the state of the winning cell in each band (r, p, rho, T, B^2, v_A, c_s, the
//!     signal and flow speed and dx in the winning direction, sin_cell, panel, the unit
//!     position vector).
//! Nothing here writes state: with the key absent nothing runs (bitwise, zero cost).
//! Ideal or general EOS, Newtonian MHD, no <mhd>/rad_signal_speed (refused when on).

#include <math.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "diffusion/resistivity.hpp"
#include "mhd.hpp"

namespace mhd {

namespace {
int rep_every = 0;
Real rep_p1 = 1.0, rep_p2 = 1.0e3;
int rep_calls = 0;
}  // namespace

void DtReportInit(ParameterInput *pin, MeshBlockPack *pp, const bool rss) {
  if (pin->DoesParameterExist("mhd", "dt_report")) {
    rep_every = pin->GetInteger("mhd", "dt_report");
    // prod1103: the report reads the curvilinear cell widths (pcoord dx1/dx2/dx3, sin of
    // the cs cell angle), so it is 3-D cubed sphere / spherical polar only; refuse the key
    // elsewhere at start-up instead of switching it off silently.  Hydro has no dt_report.
    if (rep_every > 0 && (!pp->pmesh->three_d || rss ||
        !(pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar))) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<mhd>/dt_report needs a 3-D cubed-sphere or spherical-"
                << "polar mesh and rad_signal_speed off" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    rep_p1 = pin->GetOrAddReal("mhd", "dt_report_p1", 1.0);
    rep_p2 = pin->GetOrAddReal("mhd", "dt_report_p2", 1.0e3);
  }
}

void DtReport(MHD *pm, MeshBlockPack *pp) {
  if (rep_every <= 0) return;
  if ((rep_calls++) % rep_every != 0) return;
  Mesh *pmesh = pp->pmesh;
  if (pm->rad_signal_speed || !pmesh->three_d) {
    if (rep_calls == 1) std::cout << "### dtrep: needs 3-D and rad_signal_speed off; "
                                     "report disabled" << std::endl;
    rep_every = 0;
    return;
  }
  auto &indcs = pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1, js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w0_ = pm->w0;
  auto bcc_ = pm->bcc0;
  auto wder_ = pm->wder;
  auto eos = pm->peos->eos_data;
  const bool gen = eos.IsGeneral();
  const int ipr = eos.is_ideal ? static_cast<int>(IEN) : static_cast<int>(IDN);
  const bool curv = (pmesh->use_cubed_sphere || pmesh->use_spherical_polar);
  const bool cs_ = pmesh->use_cubed_sphere;
  if (!curv) { rep_every = 0; return; }
  auto dx1_ = pp->pcoord->dx1, dx2_ = pp->pcoord->dx2, dx3_ = pp->pcoord->dx3;
  auto sn_ = pp->pcoord->sin_cell;
  const Real p1 = rep_p1, p2 = rep_p2;
  const Real va2i = eos.va2max_inv;

  // per-cell: fast / sound / flow dt (min over the three directions), pressure
  auto cell = KOKKOS_LAMBDA (int m, int k, int j, int i, Real &dtf, Real &dts,
                                   Real &dtv, Real &pr, int &dir, Real &vab) {
    const Real d = w0_(m,IDN,k,j,i);
    const Real bx = bcc_(m,IBX,k,j,i), by = bcc_(m,IBY,k,j,i), bz = bcc_(m,IBZ,k,j,i);
    Real p, g1;
    if (gen) { p = wder_(m,IDPR,k,j,i); g1 = wder_(m,IDG1,k,j,i); }
    else { p = eos.IdealGasPressure(w0_(m,ipr,k,j,i)); g1 = eos.gamma; }
    pr = p;
    const Real cs = sqrt(g1*p/d);
    const Real sn = cs_ ? sn_(m,k,j) : 1.0;
    Real cf[3];
    cf[0] = eos.FastSpeedFromP(d, p, g1, bx, by, bz);
    cf[1] = eos.FastSpeedFromP(d, p, g1, by, bz, bx)/sn;
    cf[2] = eos.FastSpeedFromP(d, p, g1, bz, bx, by)/sn;
    const Real v[3] = {fabs(w0_(m,IVX,k,j,i)), fabs(w0_(m,IVY,k,j,i)),
                       fabs(w0_(m,IVZ,k,j,i))};
    const Real h[3] = {dx1_(m,k,j,i), dx2_(m,k,j,i), dx3_(m,k,j,i)};
    const Real csd[3] = {cs, cs/sn, cs/sn};
    dtf = dts = dtv = std::numeric_limits<Real>::max(); dir = 0;
    for (int a = 0; a < 3; ++a) {
      const Real f = h[a]/(v[a] + cf[a]);
      if (f < dtf) { dtf = f; dir = a; }
      dts = fmin(dts, h[a]/(v[a] + csd[a]));
      if (v[a] > 0.0) dtv = fmin(dtv, h[a]/v[a]);
    }
    vab = (bx*bx + by*by + bz*bz)/d;   // v_A^2
  };

  using VL = Kokkos::ValLocScalar<Real, int>;
  VL ra, rb, rc;
  Real sa, sb, fa;
  Kokkos::parallel_reduce("mhd_dtrep1", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, VL &a, VL &b, VL &c, Real &s_a, Real &s_b, Real &f_a) {
    int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1 + is; k += ks; j += js;
    Real dtf, dts, dtv, pr, vab; int dir;
    cell(m, k, j, i, dtf, dts, dtv, pr, dir, vab);
    if (dtf < a.val) { a.val = dtf; a.loc = idx; }
    if (pr >= p1 && dtf < b.val) { b.val = dtf; b.loc = idx; }
    if (pr >= p2 && dtf < c.val) { c.val = dtf; c.loc = idx; }
    s_a = fmin(s_a, dts);
    if (pr >= p1) s_b = fmin(s_b, dts);
    f_a = fmin(f_a, dtv);
  }, Kokkos::MinLoc<Real, int>(ra), Kokkos::MinLoc<Real, int>(rb),
     Kokkos::MinLoc<Real, int>(rc), Kokkos::Min<Real>(sa), Kokkos::Min<Real>(sb),
     Kokkos::Min<Real>(fa));

  const Real thr = 1.2*ra.val, thr1 = 1.2*rb.val;
  int n0 = 0, n1 = 0, n2 = 0, ncap = 0, ncapt = 0, nb1 = 0;
  Kokkos::parallel_reduce("mhd_dtrep2", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, int &c0, int &c1, int &c2, int &cc, int &cct, int &cb1) {
    int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1 + is; k += ks; j += js;
    Real dtf, dts, dtv, pr, vab; int dir;
    cell(m, k, j, i, dtf, dts, dtv, pr, dir, vab);
    if (dtf < thr) {
      if (pr < p1) ++c0; else if (pr < p2) ++c1; else ++c2;
    }
    if (pr >= p1 && dtf < thr1) ++cb1;
    if (va2i > 0.0 && vab*va2i >= 0.99) { ++cc; if (pr < p1) ++cct; }
  }, Kokkos::Sum<int>(n0), Kokkos::Sum<int>(n1), Kokkos::Sum<int>(n2),
     Kokkos::Sum<int>(ncap), Kokkos::Sum<int>(ncapt), Kokkos::Sum<int>(nb1));

  // winners' state
  constexpr int nv = 18;
  DualArray2D<Real> wv("dtrep_w", 3, nv);
  const int locs[3] = {ra.loc, rb.loc, rc.loc};
  const int l0 = locs[0], l1 = locs[1], l2 = locs[2];
  auto x1v_ = pp->pcoord->x1v, x2v_ = pp->pcoord->x2v, x3v_ = pp->pcoord->x3v;
  auto pan_ = pp->pmb->mb_panel.d_view;
  auto tem_ = pm->wtemp;
  auto wvd = wv.d_view;
  Kokkos::parallel_for("mhd_dtrep3", Kokkos::RangePolicy<>(DevExeSpace(), 0, 3),
  KOKKOS_LAMBDA(const int n) {
    const int idx = (n == 0) ? l0 : ((n == 1) ? l1 : l2);
    for (int q = 0; q < nv; ++q) wvd(n,q) = -1.0;
    if (idx < 0 || idx >= nmkji) return;
    int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1 + is; k += ks; j += js;
    Real dtf, dts, dtv, pr, vab; int dir;
    cell(m, k, j, i, dtf, dts, dtv, pr, dir, vab);
    const Real d = w0_(m,IDN,k,j,i);
    const Real bx = bcc_(m,IBX,k,j,i), by = bcc_(m,IBY,k,j,i), bz = bcc_(m,IBZ,k,j,i);
    Real g1 = gen ? wder_(m,IDG1,k,j,i) : eos.gamma;
    const Real sn = cs_ ? sn_(m,k,j) : 1.0;
    Real b[3] = {bx, by, bz};
    const Real cfd = eos.FastSpeedFromP(d, pr, g1, b[dir], b[(dir+1)%3], b[(dir+2)%3])
                     /(dir > 0 ? sn : 1.0);
    const Real vd = fabs(w0_(m,IVX+dir,k,j,i));
    const Real hd = (dir == 0) ? dx1_(m,k,j,i) : ((dir == 1) ? dx2_(m,k,j,i)
                                                            : dx3_(m,k,j,i));
    Real q3[3] = {0.0, 0.0, 0.0};
    const int pnl = pan_(m);
    if (cs_) cubed_sphere::PanelToCart(pnl, 0.25*M_PI*x2v_(m,j), 0.25*M_PI*x3v_(m,k), q3);
    wvd(n,0) = x1v_(m,i); wvd(n,1) = pr; wvd(n,2) = d;
    wvd(n,3) = gen ? tem_(m,k,j,i) : pr/d;
    wvd(n,4) = bx*bx + by*by + bz*bz; wvd(n,5) = sqrt(vab);
    wvd(n,6) = sqrt(g1*pr/d); wvd(n,7) = cfd; wvd(n,8) = vd; wvd(n,9) = hd;
    wvd(n,10) = dir + 1; wvd(n,11) = sn; wvd(n,12) = pnl;
    wvd(n,13) = q3[0]; wvd(n,14) = q3[1]; wvd(n,15) = q3[2];
    wvd(n,16) = static_cast<Real>(m);
    wvd(n,17) = static_cast<Real>(((k-ks)*1000 + (j-js))*1000 + (i-is));
  });
  wv.modify_device();
  wv.sync_host();

  Real dtd = -1.0; bool fe = false;
  if (pm->presist != nullptr && pm->presist->use_rkg_sts) {
    dtd = pm->presist->dtnew; fe = pm->presist->rkg_dt_fe;
  }
  const Real cfl = pmesh->cfl_no;
  char buf[1024];
  std::snprintf(buf, sizeof(buf), "### dtrep rank=%d cycle=%d time=%.10e cfl=%.3f "
      "dtnew=%.6e | fast all=%.6e p>=p1=%.6e p>=p2=%.6e | sound all=%.6e p>=p1=%.6e | "
      "flow all=%.6e | near1.2x top=%d fe=%d deep=%d near1.2x(p>=p1)=%d | vacap all=%d "
      "top=%d | rkg_dtnew=%.6e fe=%d",
      global_variable::my_rank, pmesh->ncycle, pmesh->time, cfl, pm->dtnew, ra.val,
      rb.val, rc.val, sa, sb, fa, n0, n1, n2, nb1, ncap, ncapt, dtd, fe ? 1 : 0);
  std::cout << buf << std::endl;
  const char *tag[3] = {"all", "p>=p1", "p>=p2"};
  for (int n = 0; n < 3; ++n) {
    auto &h = wv.h_view;
    const int kji = static_cast<int>(h(n,17));
    std::snprintf(buf, sizeof(buf), "### dtrep   %-6s m=%d k=%d j=%d i=%d panel=%d "
        "r=%.6e p=%.4e rho=%.4e T=%.1f B2=%.4e vA=%.4e cs=%.4e dir=%d sig=%.4e "
        "|v|=%.4e dx=%.4e sin=%.4f q=(%.4f,%.4f,%.4f)",
        tag[n], static_cast<int>(h(n,16)), kji/1000000, (kji/1000)%1000, kji%1000,
        static_cast<int>(h(n,12)), h(n,0), h(n,1), h(n,2), h(n,3), h(n,4), h(n,5),
        h(n,6), static_cast<int>(h(n,10)), h(n,7), h(n,8), h(n,9), h(n,11), h(n,13),
        h(n,14), h(n,15));
    std::cout << buf << std::endl;
  }
}

}  // namespace mhd
