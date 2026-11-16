//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean.cpp
//! \brief Stage solve of the clean implicit radiation scheme (<rad_m1>/scheme = clean).
//! One call solves one implicit stage (or a backward-Euler step) for E, coupled to the
//! gas by Newton iteration on the exchange and a fixed-point iteration on the velocity.
//! See rad_m1_clean.hpp for the discrete equations.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn void CleanFatal
//! \brief stops the run with a message about an unsupported clean-scheme setup

static void CleanFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPreflight
//! \brief the production keys of <rad_m1> (and <time>) an input names with scheme =
//! clean (not the defaults a restart file lists): a key the clean solve replaces names
//! its clean counterpart, the others it does not read are listed once, and a key whose
//! loss would change the physics stops the run unless it is off.  The keys the clean
//! solve reads (boundaries, beams, the time step, multi-rate, the hesdirk2 tableau and
//! its backward-Euler steps) pass.

void RadiationM1::CleanPreflight(ParameterInput *pin) {
  // block/key (* ends a prefix): nullptr read, "" no effect, "!" fatal, else clean key
  static const char *tab[][2] = {
    {"rad_m1/implicit_bc_*", nullptr}, {"rad_m1/implicit_flux_x1m*", nullptr},
    {"rad_m1/implicit_ebath_*", nullptr}, {"rad_m1/implicit_beam*", nullptr},
    {"rad_m1/implicit_cfl", nullptr}, {"rad_m1/implicit_allow_multid", nullptr},
    {"rad_m1/implicit_mr_*", nullptr}, {"rad_m1/implicit_det_trace", nullptr},
    {"rad_m1/implicit_tol", "clean_tol_nl"}, {"rad_m1/implicit_lin_tol", "clean_tol"},
    {"rad_m1/implicit_maxit", "clean_maxnl"},
    {"rad_m1/implicit_lin_maxit", "clean_maxit"},
    {"rad_m1/implicit_recon", "clean_recon"},
    {"rad_m1/implicit_closure_lag", "clean_closure"},
    {"rad_m1/time2_vet", "clean_closure"}, {"rad_m1/time2_vet_col", "clean_closure"},
    {"rad_m1/time2_vet_sc", "clean_closure"},
    {"rad_m1/time2_vet_extrap", "clean_closure"}, {"rad_m1/time2_lin_tol*", "clean_tol"},
    {"rad_m1/implicit_ghost_x1min", "!"}, {"rad_m1/implicit_*", ""},
    {"rad_m1/time2_accept_fallback", ""}, {"rad_m1/time2_posfix*", ""},
    {"rad_m1/time2_k1_*", ""}, {"rad_m1/time2_face_fct", ""},
    {"time/dt_retry_t2_be", ""}};
  const bool say = (global_variable::my_rank == 0);
  std::string none;
  for (const InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_m1" && ib.block_name != "time") continue;
    for (const InputLine &ln : ib.line) {
      if (ln.param_comment == ParameterInput::default_comment) continue;
      const std::string k = ib.block_name + "/" + ln.param_name, &v = ln.param_value;
      for (const auto &t : tab) {
        const std::string p(t[0]);
        const size_t n = p.size() - (p.back() == '*');
        if (k.compare(0, n, p, 0, n) != 0 || (n == p.size() && k.size() != n)) continue;
        if (t[1] == nullptr) break;
        const std::string c(t[1]);
        if (c == "!" && v != "false" && v != "0") {
          CleanFatal("<" + k + "> = " + v + " has no counterpart in scheme = clean");
        }
        if (c.empty()) none += " " + k;
        if (say && !c.empty() && c != "!") {
          std::cout << "### WARNING (scheme = clean): <" << k << "> has no effect, the "
                    << "clean solve reads <rad_m1>/" << c << std::endl;
        }
        break;
      }
    }
  }
  if (say && !none.empty()) {
    std::cout << "### WARNING (scheme = clean): no effect:" << none << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanInit
//! \brief allocates the arrays of the clean solve and checks that the setup is supported

void RadiationM1::CleanInit() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb = pmy_pack->nmb_thispack;
  Kokkos::realloc(cln_s, nmb, CLN_N, n3, n2, n1);
  Kokkos::realloc(cln_v, nmb, 2, n3, n2, n1);
  Kokkos::realloc(cln_w, nmb, CLW_N, n3, n2, n1);
  const int nfc = cln_milne ? CLF_KK + 6 : (cln_consist ? CLF_HK + 2 : CLF_N);
  Kokkos::realloc(cln_fc1, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, nfc, n3, n2, n1);
  Kokkos::deep_copy(cln_w, 0.0);
  Kokkos::deep_copy(cln_fc1, 0.0);
  Kokkos::deep_copy(cln_fc2, 0.0);
  Kokkos::deep_copy(cln_fc3, 0.0);
  if (cln_cell) {
    Kokkos::realloc(cln_x, nmb, CLC_N, n3, n2, n1);
    Kokkos::deep_copy(cln_x, 0.0);
  }
  if (force_ref == M1_FREF_CLEAN && cln_fref.extent_int(0) != nmb) {
    Kokkos::realloc(cln_fref, nmb, 3, n3, n2, n1);   // (a restart has read it)
  }
  // multigrid levels of the cell layout: halve every direction with an even number
  // (>= 4) of cells, down to the level where no direction can be halved
  if (cln_cell) {
    std::array<int, 3> nc = {indcs.nx1, indcs.nx2, indcs.nx3};
    while (true) {
      bool coarser = false;
      for (int d = 0; d < 3; ++d) {
        if (nc[d] >= 4 && nc[d] % 2 == 0) {
          nc[d] /= 2;
          coarser = true;
        }
      }
      if (!coarser) break;
      cln_mgn.push_back(nc);
      cln_mg.emplace_back("cln_mg", nmb, CLM_N, nc[2], nc[1], nc[0]);
    }
  }
  if (cln_extrap > 0) {
    // (a restart may have read t^n and the end of the last step)
    if (cln_mn.extent_int(0) != nmb) Kokkos::realloc(cln_mn, nmb, CLN_N, n3, n2, n1);
    Kokkos::realloc(cln_mp, nmb, CLN_N, n3, n2, n1);
    if (cln_mc.extent_int(0) != nmb) Kokkos::realloc(cln_mc, nmb, CLN_N, n3, n2, n1);
    if (cln_extrap == 2) Kokkos::realloc(cln_m1, nmb, CLN_N, n3, n2, n1);
  }
  Mesh *pm = pmy_pack->pmesh;
  if (pm->multilevel || (pm->nmb_total > 1 && !cln_cell) || chat != c_light) {
    CleanFatal("<rad_m1>/scheme = clean currently needs a uniform mesh (several "
               "MeshBlocks: clean_layout = cell) and chat = c");
  }
  // the boundaries of the clean solve: open (bath, beam; outflow without a bath),
  // reflect, periodic and, in the cell layout, an imposed flux (x1, and x2 with the
  // problem generator's face fluxes)
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  for (int s = 0; s < 6; ++s) {
    const bool flux = bc[s] == M1_IBC_FLUX || (bc[s] == M1_IBC_X2FLUX && ifx2_on);
    if (!(bc[s] == M1_IBC_MARSHAK || bc[s] == M1_IBC_OUTFLOW || bc[s] == M1_IBC_REFLECT ||
          bc[s] == M1_IBC_PERIODIC || (flux && cln_cell))) {
      CleanFatal("<rad_m1>/implicit_bc_* with scheme = clean: auto | marshak | outflow | "
                 "reflect | periodic, flux (x1, x2 with the problem generator's face "
                 "fluxes) with clean_layout = cell; not efix");
    }
  }
  if (esrc_on && !cln_cell) {
    CleanFatal("<rad_m1>/esrc (problem/mlt_flux_frozen) with scheme = clean needs "
               "clean_layout = cell");
  }
  if (impl_lbath || iebath_hi_col_on) {
    CleanFatal("scheme = clean has no ghost-built inner face (problem/he_bc_linhao) and "
               "no per-column outer bath");
  }
  if ((force_ref != M1_FREF_NONE && !cln_cell) ||
      (force_ref == M1_FREF_WB_ARAD && !fref_wsplit)) {
    CleanFatal("<rad_m1>/force_reference with scheme = clean needs clean_layout = cell "
               "and, for wb_arad, the reference work given by the problem generator "
               "(force_reference_work = split)");
  }
  if (pm->use_cubed_sphere) {
    // the rays are vet_gd's directions transported conservatively (CleanFvSweep), or
    // the analytic sphere (eddington: of infinite radius); an even panel resolution
    // keeps the cell centres off the z axis (the frame of y, CleanCsFrame)
    if (!cln_cell || !(cln_sphere || cln_edd || vgd_on) || (cln_consist && !vgd_on) ||
        vgd_twin || VgdRotOn() || vra_on || vgd_async || pm->mesh_indcs.nx2 % 2 != 0) {
      CleanFatal("<rad_m1>/scheme = clean on the cubed sphere needs clean_layout = cell, "
                 "clean_closure = sphere | eddington (no clean_consistency) or closure = "
                 "vet_col with vet_gd (clean_rays = fv1 | fv, no twin, rotation or "
                 "async) and an even <mesh>/nx2");
    }
    if (vgd_on) {
      cln_rfv = std::max(cln_rfv, 1);
      vlat_taucut = std::numeric_limits<Real>::max();
      marena_on = false;
    } else {
      if (cln_edd) cln_sph_r = std::numeric_limits<Real>::max();
      cln_sphere = true;
    }
    CleanCsGeometry();
    cln_init = true;
    return;
  }
  if (pm->use_spherical_polar) {
    // the rays are vet_gd's (closure = vet_col with vet_gd) on its fixed direction set,
    // swept from the first shell, or the analytic sphere (clean_closure = sphere)
    if (!cln_cell || !(cln_sphere || vgd_on) || (cln_consist && !vgd_on) ||
        vgd_twin || VgdRotOn() || vra_on || vgd_async) {
      CleanFatal("<rad_m1>/scheme = clean on spherical polar needs clean_layout = "
                 "cell, clean_closure = sphere (no clean_consistency) or closure = "
                 "vet_col with vet_gd (vet_gd_twin = false, vet_gd_rotate_dt = 0, no "
                 "vet_rotate_avg or vet_gd_async)");
    }
    // the rays' theta walls follow the mesh: a reflecting wall mirrors (the conservative
    // transport reflects specularly, CleanFvWalls), a periodic one wraps
    if (vgd_on && !pm->use_polar_boundary) {
      vgd_thmir = (cln_rfv == 0) &&
                  (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::reflect);
      vgd_thvac = 0;
    }
    // (the sweep starts at the first shell; its intensities keep their own memory)
    vlat_taucut = std::numeric_limits<Real>::max();
    marena_on = false;
    cln_init = true;
    return;
  }
  if (cln_sphere || !vet_sc || !vet_full) {
    CleanFatal("<rad_m1>/scheme = clean on a Cartesian mesh needs closure = vet_sc with "
               "vet_tensor = full (clean_closure = sphere is for spherical polar)");
  }
  // the formal solution must provide the half-range sums and C1
  impl_beam_hr = true;
  thsw_c1 = true;
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCsGeometry
//! \brief cubed sphere: the cell geometry of the finite-volume divergence with Cartesian
//! y.  A lateral face lies in the great-circle plane through the origin and its two
//! edge points P, P' on the cube surface: normal along P x P', area (r_r^2 - r_l^2) th/2
//! with th the angle between P and P'.  A radial face takes r_f^2 S with S the vector
//! solid angle, minus half the sum of th n over the lateral faces (outward), so that the
//! six vector areas of a cell sum to zero to round-off; the centroid is
//! (r_r^4 - r_l^4) S/(4 V).  The cube points come from the global index on the panel,
//! exactly +-1 on its edges, so that both panels at a seam build the face from the same
//! Cartesian numbers (single-valued to the bit).  A lateral face's first moment int x dA
//! is (r_r^3 - r_l^3)/3 2 sin(th/2) along the bisector of P, P'.  The frame of y is that
//! of the hydro (hydro_am.hpp): r_hat of the cell centre and the unit panel tangents
//! e_xi, e_eta (cubed_sphere::PanelTangents).  The centroids and r_hat take the halo
//! exchange (a seam ghost is the neighbour's cell), then the centroid distances.

void RadiationM1::CleanCsGeometry() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, ks = indcs.ks, nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*indcs.ng, n2 = indcs.nx2 + 2*indcs.ng;
  const int n3 = indcs.nx3 + 2*indcs.ng, np = pmy_pack->pmesh->mesh_indcs.nx2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  for (int d = 0; d < 3; ++d) {
    Kokkos::realloc(cln_ga[d], nmb1 + 1, n3, n2, n1);
    Kokkos::realloc(cln_gd[d], nmb1 + 1, n3, n2, n1);
  }
  Kokkos::realloc(cln_gn, nmb1 + 1, CLG_NV, n3, n2, n1);
  auto ga0 = cln_ga[0], ga1 = cln_ga[1], ga2 = cln_ga[2];
  auto gd0 = cln_gd[0], gd1 = cln_gd[1], gd2 = cln_gd[2];
  auto gn = cln_gn;
  auto size = pmy_pack->pmb->mb_size.d_view;
  auto pnl = pmy_pack->pmb->mb_panel.d_view;
  auto rf = pmy_pack->pcoord->xx1f;
  auto vol = pmy_pack->pcoord->volume;
  par_for("cln_csgeo", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real fa[3], fb[3], fn[3];
    cubed_sphere::PanelFrame(pnl(m), fa, fb, fn);
    // tangent of the panel angle at face q along x2 (d = 0) or x3 (d = 1)
    auto tq = [&](const int d, const int q) {
      const Real x0 = (d == 0) ? size(m).x2min : size(m).x3min;
      const int s = 2*(q - ((d == 0) ? js : ks) +
                       static_cast<int>(floor(0.5*(x0 + 1.0)*np + 0.5))) - np;
      return (s == np || s == -np) ? static_cast<Real>(s/np)
                                   : tan(0.25*M_PI*static_cast<Real>(s)/np);
    };
    // angle and unit normal (+index) of face q along x2 (d = 0) or x3 (d = 1) across
    // cell o of the other direction: P(eta_r) x P(eta_l), P(xi_l) x P(xi_r)
    auto lat = [&](const int d, const int q, const int o, Real &th, Real nv[3],
                   Real bv[3]) {
      const Real t = tq(d, q), u0 = tq(1 - d, o), u1 = tq(1 - d, o + 1);
      Real p[3], pp[3];
      for (int c = 0; c < 3; ++c) {
        p[c] = (d == 0) ? fa[c]*t + fb[c]*u1 + fn[c] : fa[c]*u0 + fb[c]*t + fn[c];
        pp[c] = (d == 0) ? fa[c]*t + fb[c]*u0 + fn[c] : fa[c]*u1 + fb[c]*t + fn[c];
      }
      nv[0] = p[1]*pp[2] - p[2]*pp[1];
      nv[1] = p[2]*pp[0] - p[0]*pp[2];
      nv[2] = p[0]*pp[1] - p[1]*pp[0];
      const Real s = std::sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
      th = atan2(s, p[0]*pp[0] + p[1]*pp[1] + p[2]*pp[2]);
      const Real ip = 1.0/std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
      const Real ipp = 1.0/std::sqrt(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2]);
      for (int c = 0; c < 3; ++c) {
        nv[c] /= s;
        bv[c] = p[c]*ip + pp[c]*ipp;
      }
      const Real ib = 2.0*sin(0.5*th)/std::sqrt(bv[0]*bv[0] + bv[1]*bv[1] + bv[2]*bv[2]);
      for (int c = 0; c < 3; ++c) bv[c] *= ib;
    };
    Real th[4], nv[4][3], bv[4][3], sv[3];
    lat(0, j, k, th[0], nv[0], bv[0]);
    lat(0, j + 1, k, th[1], nv[1], bv[1]);
    lat(1, k, j, th[2], nv[2], bv[2]);
    lat(1, k + 1, j, th[3], nv[3], bv[3]);
    for (int c = 0; c < 3; ++c) {
      sv[c] = -0.5*(th[1]*nv[1][c] - th[0]*nv[0][c] + th[3]*nv[3][c] - th[2]*nv[2][c]);
    }
    const Real ss = std::sqrt(sv[0]*sv[0] + sv[1]*sv[1] + sv[2]*sv[2]);
    const Real rl = rf(m,i), rr = rf(m,i+1), al = 0.5*(rr*rr - rl*rl);
    ga0(m,k,j,i) = rl*rl*ss;
    ga1(m,k,j,i) = al*th[0];
    ga2(m,k,j,i) = al*th[2];
    for (int c = 0; c < 3; ++c) {
      gn(m,CLG_N+c,k,j,i) = sv[c]/ss;
      gn(m,CLG_N+3+c,k,j,i) = nv[0][c];
      gn(m,CLG_N+6+c,k,j,i) = nv[2][c];
      gn(m,CLG_X+c,k,j,i) = 0.25*(rr*rr*rr*rr - rl*rl*rl*rl)*sv[c]/vol(m,k,j,i);
      gn(m,CLG_M+c,k,j,i) = (rr*rr*rr - rl*rl*rl)/3.0*bv[0][c];
      gn(m,CLG_M+3+c,k,j,i) = (rr*rr*rr - rl*rl*rl)/3.0*bv[2][c];
    }
    const Real xc = 0.25*M_PI*CellCenterX(j - js, nx2, size(m).x2min, size(m).x2max);
    const Real ec = 0.25*M_PI*CellCenterX(k - ks, nx3, size(m).x3min, size(m).x3max);
    Real qc[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(pnl(m), xc, ec, qc);
    cubed_sphere::PanelTangents(pnl(m), xc, ec, e1, e2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    for (int c = 0; c < 3; ++c) {
      gn(m,CLG_R+c,k,j,i) = qc[c]*qn;
      gn(m,CLG_E+c,k,j,i) = e1[c];
      gn(m,CLG_E+3+c,k,j,i) = e2[c];
    }
  });
  CleanGhostFill(cln_gn, CLG_X, 6, false);
  par_for("cln_csdist", DevExeSpace(), 0, nmb1, 1, n3-1, 1, n2-1, 1, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int d = 0; d < 3; ++d) {
      const int kl = k - (d == 2), jl = j - (d == 1), il = i - (d == 0);
      Real s = 0.0;
      for (int c = 0; c < 3; ++c) {
        s += SQR(gn(m,CLG_X+c,k,j,i) - gn(m,CLG_X+c,kl,jl,il));
      }
      auto gd = (d == 0) ? gd0 : ((d == 1) ? gd1 : gd2);
      gd(m,k,j,i) = std::sqrt(s);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanSteadyStart
//! \brief <rad_m1>/clean_init = steady (fresh starts; Driver::Initialize at the first
//! dt): the initial radiation relaxed to the initial matter and sources.  The gas is
//! frozen (no hydro stage, no feedback: its emission and opacity stay those of the
//! initial state) and backward-Euler radiation steps of clean_init_dtfac x dt are taken
//! until the state stops changing, max(|dE|, |dF|/c) < clean_init_tol max E (the largest
//! E of the passes).  A step of dt -> infinity drops the time term, so the fixed point is
//! the steady solution of the scheme itself (rays, closure, consistency and kinetic face
//! sums of that state); the closure is a function of the state, so it has converged
//! with it.  A pass whose solve fails or leaves E < 0 is undone and the step shrinks by
//! 1000, down to dt.  The large step weights the RC term (theta h) and the moment rows
//! over the flux rows, so the passes then go on at the run's dt, to the tolerance or
//! until the change contracts by less than 1/2 per pass: what is left is the slow
//! evolution of the scheme at that dt, which the run carries on.  The radiation energy
//! and momentum are thereby redefined (reported as volume sums, F in the local basis on
//! spherical polar); the gas is not touched, except its internal energy with
//! force_reference = clean.  The closure history is then dropped, so the first step
//! starts from this state alone.

void RadiationM1::CleanSteadyStart(Driver *pdrive) {
  if (!cln_steady) return;
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  const Real dt0 = pm->dt, cl = c_light;
  const bool sfb = gas_feedback;
  const int st2 = t2_solve;
  const Real nc0 = vet_ncall;
  // gas frozen, plain backward-Euler steps without substeps; big: the large step.  With
  // force_reference = clean the gas internal energy relaxes with E (rho and momentum
  // frozen), so that the reference state is in radiative equilibrium
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool relax = (force_ref == M1_FREF_CLEAN) && sfb && fl.on;
  auto uh = fl.u0;
  gas_feedback = relax;
  cln_vzero = relax;
  t2_solve = M1_T2S_NONE;
  nsub = 1;
  bool big = true, ok = false;
  int n2 = 0;
  Real dp = 0.0, fac = cln_steady_dtfac, emax = 0.0;
  auto u0_ = u0;
  DvceArray5D<Real> uo("cln_steady_u", nmb1 + 1, 4, u0.extent(2), u0.extent(3),
                       u0.extent(4));
  DvceArray5D<Real> ug("cln_steady_g", relax ? nmb1 + 1 : 1, 4, u0.extent(2),
                       u0.extent(3), u0.extent(4));
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj*ni;
  auto &sz = pmy_pack->pmb->mb_size.h_view(0);
  const Real dv = sz.dx1*sz.dx2*sz.dx3;
  const bool sph = pm->use_spherical_polar;
  auto vol = pmy_pack->pcoord->volume;
  // totals of E and F (code units x volume)
  auto sums = [&](Real t[4]) {
    for (int c = 0; c < 4; ++c) {
      Real a = 0.0;
      Kokkos::parallel_reduce("cln_steady_sum", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              ntot),
      KOKKOS_LAMBDA(const int q, Real &sl) {
        const int m = q/(ni*nj*nk), k = ks + (q/(ni*nj))%nk;
        const int j = js + (q/ni)%nj, i = is + q%ni;
        sl += (sph ? vol(m,k,j,i) : dv)*u0_(m,M1_E+c,k,j,i);
      }, Kokkos::Sum<Real>(a));
      t[c] = CleanReduce(a, false);
    }
  };
  // u0 <-> uo, ghosts included
  auto copy = [&](const bool save) {
    Kokkos::parallel_for("cln_steady_copy", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         uo.size()),
    KOKKOS_LAMBDA(const int q) {
      const int n4 = uo.extent(4), n3 = uo.extent(3), n2 = uo.extent(2);
      const int i = q%n4, j = (q/n4)%n3, k = (q/(n4*n3))%n2, c = (q/(n4*n3*n2))%4;
      const int m = q/(n4*n3*n2*4);
      const int g = (c == 0) ? IEN : IM1 + c - 1;
      if (save) {
        uo(m,c,k,j,i) = u0_(m,M1_E+c,k,j,i);
        if (relax) ug(m,c,k,j,i) = uh(m,g,k,j,i);
      } else {
        u0_(m,M1_E+c,k,j,i) = uo(m,c,k,j,i);
        if (relax) uh(m,g,k,j,i) = ug(m,c,k,j,i);
      }
    });
    if (relax && !save) FluidRefresh(pdrive);
  };
  Real t0[4], t1[4];
  sums(t0);
  int n = 0;
  Real d = 0.0;
  for (; n < cln_steady_n; ++n) {
    pm->dt = big ? fac*dt0 : dt0;
    dt_sub = pm->dt;
    copy(true);
    const Real nf = cln_nfail + cln_nneg;
    // one radiation step: the ordinary m1 task lists (closure, solve, ghost exchange)
    for (int stage = 1; stage <= nstage; ++stage) {
      pdrive->ExecuteTaskList(pm, "m1_before_stagen", stage);
      pdrive->ExecuteTaskList(pm, "m1_stagen", stage);
      pdrive->ExecuteTaskList(pm, "m1_after_stagen", stage);
    }
    (void) ApplyClosureLimits(pdrive, 1);
    if (cln_nfail + cln_nneg > nf) {
      copy(false);
      dtr_fail = false;
      cln_have_n = cln_have_p = cln_have_c = cln_have_1 = false;
      if (global_variable::my_rank == 0) {
        std::printf("<rad_m1> clean_init steady pass %d: solve failed or E < 0 at %.0e "
                    "dt, undone\n", n, big ? fac : 1.0);
      }
      if (!big || fac <= 1.0) break;
      fac = fmax(1.0e-3*fac, 1.0);
      continue;
    }
    // its largest change of E and F/c over max E
    Real dm = 0.0, em = 0.0;
    Kokkos::parallel_reduce("cln_steady_chg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            ntot),
    KOKKOS_LAMBDA(const int q, Real &dl, Real &el) {
      const int m = q/(ni*nj*nk), k = ks + (q/(ni*nj))%nk;
      const int j = js + (q/ni)%nj, i = is + q%ni;
      Real x = fabs(u0_(m,M1_E,k,j,i) - uo(m,0,k,j,i));
      for (int c = 1; c < 4; ++c) {
        x = fmax(x, fabs(u0_(m,M1_E+c,k,j,i) - uo(m,c,k,j,i))/cl);
      }
      dl = fmax(dl, x);
      el = fmax(el, fmax(u0_(m,M1_E,k,j,i), uo(m,0,k,j,i)));
    }, Kokkos::Max<Real>(dm), Kokkos::Max<Real>(em));
    emax = fmax(emax, CleanReduce(em, true));
    d = CleanReduce(dm, true)/fmax(emax, 1.0e-300);
    if (global_variable::my_rank == 0) {
      std::printf("<rad_m1> clean_init steady pass %d: max(|dE|, |dF|/c)/max E = %.3e\n",
                  n, d);
    }
    if (!big && (d < cln_steady_tol || (++n2 > 2 && d > 0.5*dp))) break;
    if (d < cln_steady_tol) big = false;
    ok = ok || !big;
    dp = d;
  }
  // restore the run's setup; the first step starts from this state without history
  sums(t1);
  gas_feedback = sfb;
  cln_vzero = false;
  t2_solve = st2;
  t2_ok = false;
  dtr_fail = false;
  cln_have_n = cln_have_p = cln_have_c = cln_have_1 = false;
  pm->dt = dt0;
  (void) SetSubsteps(dt0);
  if (global_variable::my_rank == 0) {
    std::printf("<rad_m1> clean_init steady: %s after %d passes (%.0f formal solutions); "
                "radiation E %.10e -> %.10e, F/c %.4e %.4e %.4e -> %.4e %.4e %.4e "
                "(volume sums)\n", ok ? "converged" : "NOT converged",
                std::min(n + 1, cln_steady_n), vet_ncall - nc0, t0[0], t1[0], t0[1]/cl,
                t0[2]/cl, t0[3]/cl, t1[1]/cl, t1[2]/cl, t1[3]/cl);
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::CleanSolve
//! \brief the implicit stage solve of the clean scheme (replaces ImplicitSolve), in the
//! staggered layout (E in cells, F on faces) or the cell layout (E and F in cells)

TaskStatus RadiationM1::CleanSolve(Driver *pdriver, int stage) {
  if (!cln_init) CleanInit();
  // the vet_gd direction set before any use (the consistency term reads it at the first
  // solve of a restart, before any ray build)
  if (vgd_on && !tau_ready) VetColInit();
  FluidRef fl = FluidRef::Get(pmy_pack);
  cln_feedback = gas_feedback && fl.on;
  const bool cell = cln_cell;

  // 1. starting face fluxes on the first cycle (restarts carry them)
  if (!cln_faces_set) {
    cln_faces_set = true;
    if (pmy_pack->pmesh->ncycle == 0 && !cell) CleanInitFaceFluxes();
  }

  // 2. closure (from the gas before the stage increments are added), old vector, face
  //    coefficients
  if (cln_extrap > 0) {
    CleanClosureExtrap();
  } else if (!cln_closit) {
    CleanClosure();
  }
  CleanOldVector();
  if (cell) {
    CleanCellOldVector();
    CleanCellFaceCoefficients();
  } else {
    CleanFaceCoefficients();
  }

  // 3. Newton iteration on the gas coupling and Picard passes on the lagged quantities:
  //    the slopes (reconstruction) and, with clean_closure = iterate, the closure from
  //    the rays of the iterate (at most clean_maxpic, then they are frozen and only the
  //    gas coupling iterates on); one linear solve without any of these.  With
  //    clean_closure = extrap(2) the stage-1 solve, whose state Y1 differs from
  //    U(t^{n+1}) at O(dt^2), takes the closure of its own first iterate once (pass 1):
  //    the extrapolated closure alone makes the step first order in time.  In the cell
  //    layout a solution with cells E < 0 switches their faces to the low-order flux
  //    and is solved again (at most clean_maxpos times).
  const bool clos1 = (cln_extrap > 0) && (t2_solve == M1_T2S_STAGE1);
  const bool iterate = cln_feedback || (cell && cln_recon) || cln_closit || clos1;
  const bool lagged = (cell && cln_recon) || cln_closit || clos1;
  const bool t2stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  int nit = 0, nit_sum = 0, nsolve = 0, npic = 0, npos = 0, nmark = 0;
  bool all_low = false;
  cln_all_low = false;
  Real dpic = 0.0;
  bool converged = !iterate;
  const int maxnl = iterate ? cln_maxnl : 1;
  // the lagged slopes settle only to the accuracy of the linear solve
  const Real tol_nl = cell ? fmax(cln_tol_nl, 10.0*cln_tol) : cln_tol_nl;
  // a change that stopped decreasing within 100 times both its round-off level (cln_fe,
  // cln_fv, from the magnitudes of the terms of the solve) and tol_nl is converged too
  auto settled = [&](const Real d, const Real d0, const Real f) {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  };
  for (int ipos = 0; ; ++ipos) {
    CleanFirstIterate();
    if (cell) CleanCellFirstIterate();
    if (cln_feedback) CleanAdvectiveCoefficients();
    converged = !iterate;
    npic = 0;
    Real de0 = 1.0e300, dv0 = 1.0e300;
    for (int inl = 0; inl < maxnl; ++inl) {
      if (lagged && inl < cln_maxpic) ++npic;
      const bool clos = (cln_closit && inl < cln_maxpic) || (clos1 && inl == 1);
      if (clos) {
        CleanClosureIterate();
        if (cell) {
          CleanCellFaceSpeeds();
        } else {
          CleanFaceCoefficients();
        }
        if (cln_feedback) CleanAdvectiveCoefficients();
      }
      if (cell) {
        if (cln_consist && (inl == 0 || clos)) CleanCellConsistency();
        if (cln_recon && inl < cln_maxpic) CleanCellSlopeWeights(CLC_K);
        CleanCellRightHandSide();
        CleanCellPreconditionerSetup();
      } else {
        CleanRightHandSide();
        CleanPreconditionerSetup();
      }
      nit = CleanBiCGStab();
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (!iterate || nit < 0) break;
      const Real de = cell ? CleanCellNewtonUpdate() : CleanNewtonUpdate();
      if (lagged && inl < cln_maxpic) dpic = de;
      Real dv = 0.0;
      if (cln_feedback) dv = cell ? CleanCellVelocityUpdate() : CleanVelocityUpdate();
      if (settled(de, de0, cln_fe) && settled(dv, dv0, cln_fv) && !(clos1 && inl == 0)) {
        converged = true;
        break;
      }
      de0 = de;
      dv0 = dv;
    }
    if (iterate) {
      if (cell) {
        CleanCellCopyIterate();
      } else {
        CleanCopyIterate();
      }
    }
    if (!cell || !cln_poslim) break;
    const int nbad = CleanCellMarkFaces();
    if (nbad == 0) break;
    nmark += nbad;
    if (ipos >= cln_maxpos) {
      // a stage gives up (the step is redone by backward Euler); a backward-Euler step
      // is solved once more with every face low-order
      if (t2stage || all_low) break;
      CleanCellAllLowOrder();
      all_low = true;
    }
    ++npos;
  }
  if (!converged) nit = -std::abs(nit) - 1;

  // 4. positivity of the solution
  const int nneg = cell ? CleanCellCountNegative() : CleanCountNegative();
  if (cell && cln_diagcells && all_low) CleanCellDiagnose();
  if (nit < 0) cln_nfail += 1.0;
  if (nneg > 0) cln_nneg += 1.0;
  // a backward-Euler solve that did not converge fails the cycle: with <time>/dt_retry
  // the cycle is taken again with a smaller dt (a failed stage is redone by backward
  // Euler through t2_fail below)
  if (nit < 0 && !t2stage) {
    if (dtr_on) {
      if (!dtr_fail) {
        dtr_why = 1;
        dtr_resid = -1.0;
      }
      dtr_fail = true;
    }
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d: backward-Euler solve not converged (%s, %d "
                  "so far)\n", pmy_pack->pmesh->ncycle,
                  dtr_on ? "dt_retry" : "no dt_retry: accepted",
                  static_cast<int>(cln_nfail));
    }
  }
  if (global_variable::my_rank == 0) {
    if (cell) {
      std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pic %d %.2e pos %d mark %d "
                  "all %d neg %d nreal %d guard %d solves %d\n", pmy_pack->pmesh->ncycle,
                  t2_solve, nit, nit_sum, npic, dpic, npos, nmark, all_low ? 1 : 0, nneg,
                  cln_nnreal, cln_nguard, nsolve);
    } else {
      std::printf("CLEAN cycle %d t2s %d it %d neg %d guard %d\n",
                  pmy_pack->pmesh->ncycle, t2_solve, nit, nneg, cln_nguard);
    }
  }

  // 5. new state: face fluxes, E, cell fluxes, gas, stage slopes
  if (cell) {
    CleanCellWriteBack();
  } else {
    CleanWriteBack();
  }
  // a write-back that leaves gas eint < 0 (energy for the fluid floor to create) is
  // rejected like a failed solve (the stage is redone by backward Euler, a
  // backward-Euler step by dt_retry)
  if (cell && cln_feedback && cln_neint > 0 && !t2stage && dtr_on) dtr_fail = true;
  if (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2) {
    t2_fail = (nit < 0) || (nneg > 0) || (cell && cln_feedback && cln_neint > 0);
  }

  // 6. clean_closure = extrap(2): ray moments of the new state.  The last solve of a
  //    step also re-evaluates its slope K1 = f(U^{n+1}) with the closure of these rays
  //    (the difference quotient of the solve carries the extrapolated closure, whose
  //    error the stiff coupling turns into an O(dt) error of the slope; extrap1 keeps the
  //    slope, whose stage-2 closure is that of stage 1)
  if (cln_extrap > 0 && !t2_fail) {
    const bool reslope = cell && cln_extrap != 3 && t2_solve != M1_T2S_NONE &&
                         t2_solve != M1_T2S_STAGE1;
    if (reslope) CleanCellClosureSlope(-1.0);
    CleanClosureAfterSolve();
    if (reslope) {
      CleanClosureFill();
      CleanGhostFill(cln_s, 0, CLN_N);
      CleanCellFaceSpeeds();
      if (cln_consist) CleanCellConsistency();
      if (cln_milne) CleanCellMilneSetup();
      CleanCellClosureSlope(1.0);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanInitFaceFluxes
//! \brief face fluxes as the mean of the two adjacent cell fluxes

void RadiationM1::CleanInitFaceFluxes() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u0_ = u0;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("cln_f0", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      f0(m,k,j,i) = 0.5*(u0_(m,M1_F1+dir,k,j,i) + u0_(m,M1_F1+dir,k-dk,j-dj,i-di));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanOldVector
//! \brief E and gas state the stage starts from (stage increments added), with the gas
//! internal energy and temperature

void RadiationM1::CleanOldVector() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  const bool feedback = cln_feedback;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fluid = fl.on;
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto u0_ = u0;
  auto inc = t2inc;
  auto w = cln_w;
  CleanArgs x = CleanMakeArgs();

  par_for("cln_old", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real t = 0.0;
    if (fluid) {
      if (feedback && stage) {
        for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += inc(m,M1_T2_M1+n,k,j,i);
        uh(m,IEN,k,j,i) += inc(m,M1_T2_EN,k,j,i);
      }
      const Real d = uh(m,IDN,k,j,i);
      Real p[3];
      CleanGasMomentum(x, uh, m, k, j, i, p);
      const Real ekin = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/d;
      const Real eint = fmax(uh(m,IEN,k,j,i) - ekin, 0.0);
      t = eos.Temperature(d, eint);
      w(m,CLW_EI,k,j,i) = eint;
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,M1_E,k,j,i) + (stage ? inc(m,M1_T2_E,k,j,i) : 0.0);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanAdvectiveCoefficients
//! \brief a_d = v_d + D_db v_b, the advective flux per unit E, with v the gas velocity of
//! the old vector plus the momentum change of the current iterate

void RadiationM1::CleanAdvectiveCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto w = cln_w;
  CleanArgs x = CleanMakeArgs();

  par_for("cln_advcoef", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i);
    Real v[3], a[3];
    CleanGasMomentum(x, uh, m, k, j, i, v);
    for (int n = 0; n < 3; ++n) v[n] = (v[n] + w(m,CLW_DM+n,k,j,i))/d;
    CleanAdvect(x, m, k, j, i, v, a);
    for (int n = 0; n < 3; ++n) w(m,CLW_A1+n,k,j,i) = a[n];
  });
  CleanGhostFill(cln_w, CLW_A1, 3);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFirstIterate
//! \brief first Newton iterate: E of the old vector, no momentum change, and the gas
//! temperature in local exchange equilibrium with that E

void RadiationM1::CleanFirstIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool feedback = cln_feedback;
  const Real hc = dt_sub*c_light, arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;

  par_for("cln_iter0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = fmax(w(m,CLW_EO,k,j,i), 0.0);
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) w(m,CLW_DM+n,k,j,i) = 0.0;
    if (feedback) {
      w(m,CLW_T,k,j,i) = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                             hc*op(m,M1_OP_P,k,j,i)*arad_,
                             w(m,CLW_EI,k,j,i) + hc*op(m,M1_OP_E,k,j,i)*e,
                             w(m,CLW_T,k,j,i));
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRightHandSide
//! \brief B = E_old + h c kP a T^4 - r h c kE E_k - W - h div F_source, with the exchange
//! linearised about the iterate (r = s/(d cv + s), s = 4 h c kP a T^3), W the work of
//! the radiation force of the iterate on the gas, and the initial guess X = E_k

void RadiationM1::CleanRightHandSide() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const Real arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;

  par_for("cln_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hckp = x.h*x.c*op(m,M1_OP_P,k,j,i);
    const Real hcke = x.h*x.c*op(m,M1_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real e, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
      const Real s = 4.0*hckp*arad_*t*t*t;
      rr = s/(d*cv + s);
      // work of the radiation force on the gas, (v_old + dm/(2 rho)).dm: exactly the
      // kinetic energy the gas gains, taken from the radiation
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (uh(m,IM1+n,k,j,i) + 0.5*dm)*dm/d;
      }
    }
    x.w(m,CLW_RR,k,j,i) = rr;
    x.w(m,CLW_B,k,j,i) = x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
                         - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work
                         - CleanSourceDivergence(x, m, k, j, i);
    x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanNewtonUpdate
//! \brief new iterate E_k = X and the gas temperature in exchange equilibrium with it
//! (the absorption takes max(X, 0); the staggered layout also keeps that as the iterate,
//! the cell layout keeps X, whose E < 0 its positivity limiter treats).
//! Returns max |dE|/E, |dT|/T.

Real RadiationM1::CleanNewtonUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const Real hc = dt_sub*c_light, arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  const bool cell = cln_cell;

  Real change = 0.0;
  Kokkos::parallel_reduce("cln_newton", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dmax) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real e_old = w(m,CLW_EK,k,j,i);
    const Real t_old = w(m,CLW_T,k,j,i);
    const Real x = w(m,CLW_X,k,j,i);
    const Real ex = fmax(x, 0.0);
    const Real e = cell ? x : ex;
    const Real t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                       hc*op(m,M1_OP_P,k,j,i)*arad_,
                       w(m,CLW_EI,k,j,i) + hc*op(m,M1_OP_E,k,j,i)*ex, t_old);
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_T,k,j,i) = t;
    const Real de = fabs(e - e_old)/fmax(fmax(e, e_old), 1.0e-300);
    const Real dt = fabs(t - t_old)/fmax(fmax(t, t_old), 1.0e-300);
    dmax = fmax(dmax, fmax(de, dt));
  }, Kokkos::Max<Real>(change));
  return CleanReduce(change, true);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanVelocityUpdate
//! \brief momentum change of the new iterate and the advective coefficients of the new
//! velocity (implicit velocity).  Returns max |dv| over max |v|.

Real RadiationM1::CleanVelocityUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanGhostFill(cln_w, CLW_EK, 1);
  CleanArgs x = CleanMakeArgs();
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;

  Real dvmax = 0.0, vmax = 0.0;
  Kokkos::parallel_reduce("cln_velocity", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real d = uh(m,IDN,k,j,i);
    Real v2 = 0.0;
    for (int dir = 0; dir < 3; ++dir) {
      Real dm = 0.0;
      if (dir == 0 || (dir == 1 && x.multi_d) || (dir == 2 && x.three_d)) {
        dm = CleanMomentumChange(x, CLW_EK, dir, m, k, j, i);
      }
      dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+dir,k,j,i))/d);
      x.w(m,CLW_DM+dir,k,j,i) = dm;
      v2 += SQR((uh(m,IM1+dir,k,j,i) + dm)/d);
    }
    v_max = fmax(v_max, std::sqrt(v2));
  }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax));
  CleanAdvectiveCoefficients();
  return CleanReduce(dvmax, true)/fmax(CleanReduce(vmax, true), 1.0e-300);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCopyIterate
//! \brief the solution is the last Newton iterate

void RadiationM1::CleanCopyIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  par_for("cln_copy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_X,k,j,i) = w(m,CLW_EK,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanCountNegative
//! \brief sets values within the solver tolerance of zero (|E| <= 10 tol max|B|) to zero
//! and returns the number of cells that are still negative

int RadiationM1::CleanCountNegative() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w = cln_w;

  Real bmax = 0.0;
  Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(w(m,CLW_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));

  const Real zero = 10.0*cln_tol*CleanReduce(bmax, true);
  Real nneg = 0.0;
  Kokkos::parallel_reduce("cln_neg", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    Real &e = w(m,CLW_X,k+ks,j+js,i+is);
    if (e < 0.0 && e >= -zero) e = 0.0;
    if (!(e >= 0.0)) n += 1.0;
  }, Kokkos::Sum<Real>(nneg));
  return static_cast<int>(CleanReduce(nneg, false));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanWriteBack
//! \brief stores the solution: face fluxes, E, cell fluxes (face means), the gas momentum
//! and energy, and the stage slopes K = (U - U_old)/h.  The gas takes exactly the momentum
//! and energy the radiation gives up, so both totals change only by face fluxes.

void RadiationM1::CleanWriteBack() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool slopes = (t2_solve != M1_T2S_NONE);
  const bool feedback = cln_feedback;
  const Real h = dt_sub;
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanArgs x = CleanMakeArgs();
  auto kk = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;

  // face fluxes of the solution
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = M1_T2_F1 + dir;
    par_for("cln_wbface", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real f = CleanFaceFlux(x, CLW_X, dir, m, k, j, i, true);
      if (slopes) kk(m,tf,k,j,i) = (f - fc(m,CLF_G,k,j,i))/h;
      f0(m,k,j,i) = f;
    });
  }

  // E, cell fluxes, gas
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto u0_ = u0;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  par_for("cln_wbcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = x.w(m,CLW_X,k,j,i);
    Real dm[3] = {0.0, 0.0, 0.0};
    Real den = 0.0;
    if (feedback) {
      dm[0] = CleanMomentumChange(x, CLW_X, 0, m, k, j, i);
      if (multi_d) dm[1] = CleanMomentumChange(x, CLW_X, 1, m, k, j, i);
      if (three_d) dm[2] = CleanMomentumChange(x, CLW_X, 2, m, k, j, i);
      for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += dm[n];
      // gas energy: exactly what the radiation gives up besides its flux divergence, so
      // that total energy changes only by the face fluxes
      Real div = x.idx[0]*(f1(m,k,j,i+1) - f1(m,k,j,i));
      if (multi_d) div += x.idx[1]*(f2(m,k,j+1,i) - f2(m,k,j,i));
      if (three_d) div += x.idx[2]*(f3(m,k+1,j,i) - f3(m,k,j,i));
      den = -(e - x.w(m,CLW_EO,k,j,i)) - h*div;
      uh(m,IEN,k,j,i) += den;
    }
    if (slopes) {
      kk(m,M1_T2_E,k,j,i) = (e - x.w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,M1_T2_M1+n,k,j,i) = dm[n]/h;
      kk(m,M1_T2_EN,k,j,i) = den/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    u0_(m,M1_F1,k,j,i) = 0.5*(f1(m,k,j,i) + f1(m,k,j,i+1));
    u0_(m,M1_F1+1,k,j,i) = multi_d ? 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i)) : 0.0;
    u0_(m,M1_F1+2,k,j,i) = three_d ? 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i)) : 0.0;
  });
}

}  // namespace radm1
