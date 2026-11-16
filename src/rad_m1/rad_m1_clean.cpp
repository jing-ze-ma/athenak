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
extern bool g_vra_adv;   // vet_rotate_advance (rad_m1_vetgd.cpp)
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
    {"rad_m1/implicit_ghost_x1min", "!"}, {"rad_m1/ic_rad*", "clean_init"},
    {"rad_m1/implicit_*", ""},
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
  Kokkos::realloc(cln_v, nmb, 3, n3, n2, n1);
  Kokkos::realloc(cln_w, nmb, CLW_N, n3, n2, n1);
  // the cubed sphere takes the half ranges about the face normals from the FV rays
  cln_hk = cln_consist || (pmy_pack->pmesh->use_cubed_sphere && vgd_on);
  const int nfc = std::max(cln_milne ? CLF_KK + 6 : (cln_hk ? CLF_HK + 2 : CLF_N),
                           cln_fsolve ? CLF_N + 4 : 0);
  // the opacities depend on T: tables, or the T exponent (bump) of powerlaw (user)
  cln_opt = cln_cell && !opac_freeze && !opac_zero && (opacity_type == M1_OPAC_TABLE ||
            (opacity_type != M1_OPAC_CONST && opac_b != 0.0));
  Kokkos::realloc(cln_fc1, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, nfc, n3, n2, n1);
  Kokkos::deep_copy(cln_w, 0.0);
  Kokkos::deep_copy(cln_fc1, 0.0);
  Kokkos::deep_copy(cln_fc2, 0.0);
  Kokkos::deep_copy(cln_fc3, 0.0);
  for (int d = 0; d < 3 && cln_rstd.extent_int(0) == nmb && cln_fsolve; ++d) {
    auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);   // (from a restart)
    auto AL = Kokkos::ALL;
    Kokkos::deep_copy(Kokkos::subview(fc, AL, CLF_N + 3, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, d, AL, AL, AL));
  }
  if (cln_cell) {
    Kokkos::realloc(cln_x, nmb, CLC_N, n3, n2, n1);
    Kokkos::deep_copy(cln_x, 0.0);
  }
  if (force_ref == M1_FREF_CLEAN && cln_fref.extent_int(0) != nmb) {
    Kokkos::realloc(cln_fref, nmb, 3, n3, n2, n1);   // (a restart has read it)
  }
  // multigrid levels of the cell layout: halve every direction with an even number
  // (>= 4) of cells, down to the level where no direction can be halved (clean_pc =
  // prodmg: x1 is never halved, two coarse levels)
  if (cln_cell && cln_pc != CLN_PC_X1LINE) {
    const bool pm = (cln_pc == CLN_PC_PRODMG);
    std::array<int, 3> nc = {indcs.nx1, indcs.nx2, indcs.nx3};
    while (!(pm && cln_mgn.size() == 2)) {
      bool coarser = false;
      for (int d = (pm ? 1 : 0); d < 3; ++d) {
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
  if (pm->multilevel || chat != c_light) {
    CleanFatal("<rad_m1>/scheme = clean currently needs a uniform mesh and chat = c");
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
  if ((force_ref == M1_FREF_WB_ARAD && !cln_cell) ||
      (force_ref == M1_FREF_WB_ARAD && !fref_wsplit)) {
    CleanFatal("<rad_m1>/force_reference with scheme = clean needs clean_layout = cell "
               "and, for wb_arad, the reference work given by the problem generator "
               "(force_reference_work = split)");
  }
  if (vgd_on && VgdRotOn()) {
    // the rotating FV rays (CleanFvSets): the closure is the J-weighted mean of the
    // sets, continuous across the blocks with the cross-fade; no slope correction
    const bool fixed = vra_on && !g_vra_adv;
    if (vra_ema || vra_eq || (vgd_xf && !(vgd_rot_dt > 0.0 || vgd_rot_auto)) ||
        (cln_rfv == 2 && !fixed) || (vgd_rot3 && !pm->three_d)) {
      CleanFatal("<rad_m1>/scheme = clean with a rotating vet_gd set: vet_rotate_avg = K "
                 "(J-weighted, no _ema), vet_gd_rotate_xfade with vet_gd_rotate_dt, "
                 "clean_rays = fv only for a fixed union (vet_rotate_avg without "
                 "vet_rotate_advance), and on a 2-D axisymmetric mesh the z axis");
    }
    vgd_rot_slope = 0;
  }
  const bool ring_ok = vgd_on && cln_rfv > 0 && vra_on && g_vra_adv &&
                       (vgd_rot_dt > 0.0 || vgd_rot_auto) && vgd_rot_align &&
                       vra_k <= M1_RA_NF;
  if (!cln_ring_named && ring_ok && vgd_xf) cln_ring = 1;
  if (cln_ring && !(ring_ok && vgd_xf == (cln_ring == 1))) {
    CleanFatal("<rad_m1>/clean_rot_ring needs vet_gd FV rays rotating in time: "
               "vet_rotate_avg = K with vet_rotate_advance, vet_gd_rotate_dt and "
               "_align, with vet_gd_rotate_xfade for live, without for jump");
  }
  // jump: the first step of a block is a backward-Euler step (no slope of the old sets)
  if (cln_ring == 2) vgd_rot_slope = 2;
  if (pm->use_cubed_sphere) {
    // the rays are vet_gd's directions transported conservatively (CleanFvSweep), or
    // the analytic sphere (eddington: of infinite radius); an even panel resolution
    // keeps the cell centres off the z axis (the frame of y, CleanCsFrame)
    // (staggered: without gas coupling, clean_consistency and clean_stag_fsolve so far)
    if ((!cln_cell && (gas_feedback || cln_consist || cln_fsolve)) ||
        !(cln_sphere || cln_edd || vgd_on) || (cln_consist && !vgd_on) ||
        vgd_twin || vgd_async || pm->mesh_indcs.nx2 % 2 != 0) {
      CleanFatal("<rad_m1>/scheme = clean on the cubed sphere needs clean_layout = cell "
                 "(staggered: no gas_feedback, clean_consistency, clean_stag_fsolve), "
                 "clean_closure = sphere | eddington (no clean_consistency) or closure = "
                 "vet_col with vet_gd (clean_rays = fv1 | fv, no twin or async) and an "
                 "even <mesh>/nx2");
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
    // (staggered: also clean_closure = eddington)
    const bool stg = !cln_cell;
    if ((!cln_cell && !stg) || !(cln_sphere || vgd_on || (stg && cln_edd)) ||
        (cln_consist && !vgd_on) || vgd_twin || vgd_async ||
        (VgdRotOn() && cln_rfv == 0)) {
      CleanFatal("<rad_m1>/scheme = clean on spherical polar needs clean_layout = "
                 "cell, clean_closure = sphere (no clean_consistency) or closure = "
                 "vet_col with vet_gd (vet_gd_twin = false, no vet_gd_async; a rotating "
                 "set needs clean_rays = fv1 | fv)");
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
  if (cln_rfv > 0) CleanCartGeometry();
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCartGeometry
//! \brief Cartesian mesh with the FV rays (clean_rays = fv1 | fv | lld): the cell
//! geometry the rays read through the curvilinear hooks (volumes, face areas, centre and
//! face coordinates)

void RadiationM1::CleanCartGeometry() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto *pc = pmy_pack->pcoord;
  const int nmb = pmy_pack->nmb_thispack, ng = indcs.ng;
  const int c1 = indcs.nx1 + 2*ng, c2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int c3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  Kokkos::realloc(pc->volume, nmb, c3, c2, c1);
  Kokkos::realloc(pc->area.x1f, nmb, c3, c2, c1 + 1);
  Kokkos::realloc(pc->area.x2f, nmb, c3, c2 + 1, c1);
  Kokkos::realloc(pc->area.x3f, nmb, c3 + 1, c2, c1);
  Kokkos::realloc(pc->x1v, nmb, c1);
  Kokkos::realloc(pc->x2v, nmb, c2);
  Kokkos::realloc(pc->x3v, nmb, c3);
  Kokkos::realloc(pc->xx1f, nmb, c1 + 1);
  Kokkos::realloc(pc->xx2f, nmb, c2 + 1);
  Kokkos::realloc(pc->xx3f, nmb, c3 + 1);
  auto &size = pmy_pack->pmb->mb_size;
  auto vol = pc->volume;
  auto a1 = pc->area.x1f, a2 = pc->area.x2f, a3 = pc->area.x3f;
  auto x1v = pc->x1v, x2v = pc->x2v, x3v = pc->x3v;
  auto x1f = pc->xx1f, x2f = pc->xx2f, x3f = pc->xx3f;
  par_for("cln_cart_geo", DevExeSpace(), 0, nmb - 1, 0, c3, 0, c2, 0, c1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d1 = size.d_view(m).dx1, d2 = size.d_view(m).dx2, d3 = size.d_view(m).dx3;
    const bool i_ = i < c1, j_ = j < c2, k_ = k < c3;
    if (j_ && k_) a1(m,k,j,i) = d2*d3;
    if (i_ && k_) a2(m,k,j,i) = d1*d3;
    if (i_ && j_) a3(m,k,j,i) = d1*d2;
    if (!(i_ && j_ && k_)) return;
    vol(m,k,j,i) = d1*d2*d3;
    if (j == 0 && k == 0) {
      x1f(m,i) = size.d_view(m).x1min + (i - is)*d1;
      x1v(m,i) = x1f(m,i) + 0.5*d1;
      if (i == c1 - 1) x1f(m,c1) = x1f(m,i) + d1;
    }
    if (i == 0 && k == 0) {
      x2f(m,j) = size.d_view(m).x2min + (j - js)*d2;
      x2v(m,j) = x2f(m,j) + 0.5*d2;
      if (j == c2 - 1) x2f(m,c2) = x2f(m,j) + d2;
    }
    if (i == 0 && j == 0) {
      x3f(m,k) = size.d_view(m).x3min + (k - ks)*d3;
      x3v(m,k) = x3f(m,k) + 0.5*d3;
      if (k == c3 - 1) x3f(m,c3) = x3f(m,k) + d3;
    }
  });
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
    cln_tc = -1.0e300;   // (predict: the closure of this pass's state, not of the last dt's)
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
  cln_tc = -1.0e300;
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
  if (vgd_on && cln_rfv > 0 && cln_fvn == 0) {
    CleanFvSets(pmy_pack->pmesh->time, pmy_pack->pmesh->ncycle - 1);
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  cln_feedback = gas_feedback && fl.on;
  cln_poslim = cln_poslim0 && (!cln_posgas || cln_feedback);
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
  if (cln_ring == 2) CleanRingSolveSet();
  CleanOldVector();
  if (cell) {
    CleanCellOldVector();
    CleanCellFaceCoefficients();
  } else {
    CleanFaceCoefficients();
    if (cln_consist) CleanStagBracket();
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
  std::array<int, 2> nfirst = {0, 0};
  cln_nclip = 0;
  bool all_low = false;
  cln_all_low = false;
  Real dpic = 0.0;
  bool converged = !iterate;
  const int maxnl = iterate ? cln_maxnl : 1;
  // the lagged slopes (and the staggered velocity) settle only to the accuracy of the
  // linear solve
  const Real tol_nl = fmax(cln_tol_nl, 10.0*cln_tol);
  // a change that stopped decreasing within 100 times both its round-off level (cln_fe,
  // cln_fv, from the magnitudes of the terms of the solve) and tol_nl is converged too
  auto settled = [&](const Real d, const Real d0, const Real f) {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  };
  if (cell) Kokkos::deep_copy(Kokkos::subview(cln_v, Kokkos::ALL, 2, Kokkos::ALL,
                                              Kokkos::ALL, Kokkos::ALL), 0.0);
  for (int ipos = 0; ; ++ipos) {
    CleanFirstIterate();
    if (cell) CleanCellFirstIterate();
    if (cln_opt && cln_feedback) {
      CleanCellOpacity();
      CleanCellFaceSpeeds();
    }
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
          if (cln_consist) CleanStagBracket();
        }
        if (cln_feedback) CleanAdvectiveCoefficients();
      }
      if (cln_opt && cln_feedback && inl > 0) {   // the opacities of the iterate's T
        CleanCellOpacity();
        CleanCellFaceSpeeds();
      }
      if (cell) {
        if (cln_hk && (inl == 0 || clos || (cln_opt && cln_feedback))) {
          CleanCellConsistency();
        }
        if (cln_recon && inl < cln_maxpic) CleanCellSlopeWeights(CLC_K);
        CleanCellRightHandSide();
        if (inl % cln_pcreuse == 0) CleanCellPreconditionerSetup();
      } else {
        CleanRightHandSide();
        CleanPreconditionerSetup();
      }
      nit = CleanBiCGStab();
      if (cell && cln_yres && cln_feedback) CleanCellResidual();
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (!iterate || nit < 0) break;
      const Real de = cell ? CleanCellNewtonUpdate() : CleanNewtonUpdate();
      if (lagged && inl < cln_maxpic) dpic = de;
      Real dv = 0.0;
      if (cln_feedback) dv = cell ? CleanCellVelocityUpdate() : CleanVelocityUpdate();
      if (cln_nlprint && global_variable::my_rank == 0) {
        std::printf("CLEAN pass cycle %d t2s %d pos %d pass %d it %d de %.3e dv %.3e\n",
                    pmy_pack->pmesh->ncycle, t2_solve, ipos, inl, nit, de, dv);
      }
      // the remaining change of a contraction q < 1/2 is q/(1-q) times that of the pass
      auto est = [&](const Real d, const Real d0) {
        const Real q = d/d0;
        return (q < 0.5) ? q/(1.0 - q)*d : d;
      };
      if (cln_conv_est && inl > 0 && est(de, de0) < tol_nl && est(dv, dv0) < tol_nl &&
          !(clos1 && inl == 0)) {
        converged = true;
        break;
      }
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
    // (staggered: R1, R2 of the solution; with clean_pos_limit the faces of their cells
    // and the R2 faces take the rays' half-range flux, w = 1)
    const int nbad = cell ? (cln_poslim ? CleanCellMarkFaces() : 0) :
                     CleanStagRealizable(cln_poslim && ipos < cln_maxpos);
    if (!cell && ipos == 0) nfirst = {cln_nr1, cln_nr2};
    if (!cln_poslim || nbad == 0) break;
    nmark += nbad;
    if (ipos >= cln_maxpos) {
      // a stage gives up (the step is redone by backward Euler); a backward-Euler step
      // is solved once more with every face low-order
      if (t2stage || all_low) break;
      if (cell) CleanCellAllLowOrder();
      if (!cell) CleanStagRealizable(2);
      all_low = true;
    }
    ++npos;
  }
  if (!converged) nit = -std::abs(nit) - 1;

  // 4. positivity of the solution
  const int nneg = cell ? CleanCellCountNegative() : CleanCountNegative();
  if (cell && cln_diagcells && all_low) CleanCellDiagnose();
  if (!cell && cln_diagforce) CleanStagForceDiag();
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
      std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pos %d mark %d all %d neg %d "
                  "r1 %d r2 %d first %d %d guard %d solves %d\n", pmy_pack->pmesh->ncycle,
                  t2_solve, nit, nit_sum, npos, nmark, all_low ? 1 : 0, nneg, cln_nr1,
                  cln_nr2, nfirst[0], nfirst[1], cln_nguard, nsolve);
      if (cln_split) std::printf("CLEAN SPLIT cycle %d t2s %d r3 %d max %.3e clip %d "
                                 "Emin/Emax %.3e Eneg/E %.3e\n", pmy_pack->pmesh->ncycle,
                                 t2_solve, cln_nr3, cln_fh, cln_nclip, cln_eneg[0],
                                 cln_eneg[1]);
    }
  }

  // 5. new state: face fluxes, E, cell fluxes, gas, stage slopes
  if (cell) {
    CleanCellWriteBack();
  } else {
    cln_dok = (nit >= 0);   // (a failed solve does not change the energy flux's memory)
    CleanWriteBack();
  }
  // a write-back that leaves gas eint < 0 (energy for the fluid floor to create) is
  // rejected like a failed solve (the stage is redone by backward Euler, a
  // backward-Euler step by dt_retry)
  if (cln_feedback && cln_neint > 0 && !t2stage && dtr_on) dtr_fail = true;
  if (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2) {
    t2_fail = (nit < 0) || (nneg > 0) || (cln_feedback && cln_neint > 0) ||
              (!cell && cln_poslim && cln_nr2 > 0);
  }
  // (force_reference = clean: the reference is the force of the converged, accepted
  // first coupled solve, every iterate of which recaptures it; a failed one is redone)
  if (force_ref == M1_FREF_CLEAN && cln_feedback && !cln_vzero && nit >= 0 &&
      nneg == 0 && !(t2stage && t2_fail) && !dtr_fail) cln_fref_set = true;

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
      if (cln_hk) CleanCellConsistency();
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
      const Real eint = fmax(CleanGasEint(x, uh, ekin, m, k, j, i), 0.0);
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
//! the radiation force of the iterate on the gas, and the initial guess X = E_k.  With
//! gas coupling the beta terms of Ma et al. (2025, eq. 56; see rad_m1_clean_cell.hpp):
//! the emission takes G = 1 - beta^2/2 (CLW_GM) and B adds the lagged A - hc (kP a T^4 G
//! - kE E) = hc kE [(beta^2/2 - beta.a/c) E_k + 2 beta.y_k] (y_k = F/c of the iterate,
//! mean of the two faces; the gas eint takes its opposite, CLW_EP), with beta = v/c of
//! the iterate (CLW_BT) and hc kP a T^4 (CLW_RE) for the face fluxes (eq. 57,
//! CleanFaceWeights)

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

  if (x.feedback) {
    par_for("cln_beta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = 0; n < 3; ++n) {
        x.w(m,CLW_BT+n,k,j,i) = (uh(m,IM1+n,k,j,i) + x.w(m,CLW_DM+n,k,j,i))/
                                (uh(m,IDN,k,j,i)*x.c);
      }
      x.w(m,CLW_RE,k,j,i) = x.h*x.c*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,i)));
    });
    CleanGhostFill(cln_w, CLW_BT, 4);
    CleanGhostFill(cln_w, CLW_EK, 1);
  }
  CleanArgs xh = x;   // (the face momentum state, clean_stag_split)
  xh.hi = cln_split;
  par_for("cln_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hcke = x.h*x.c*op(m,M1_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0, gm = 1.0, sb = 0.0;
    if (x.feedback) {
      Real b2 = 0.0, ba = 0.0, by = 0.0;
      for (int n = 0; n < 3; ++n) {
        const Real bt = x.w(m,CLW_BT+n,k,j,i);
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        if (n < x.nf) {
          by += bt*0.5*(CleanFaceFlux(xh, CLW_EK, n, m, k, j, i, true) +
                        CleanFaceFlux(xh, CLW_EK, n, m, k+dk, j+dj, i+di, true))/x.c;
        }
        b2 += SQR(bt);
        ba += bt*x.w(m,CLW_A1+n,k,j,i)/x.c;
      }
      gm = 1.0 - 0.5*b2;
      sb = hcke*((0.5*b2 - ba)*x.w(m,CLW_EK,k,j,i) + 2.0*by);
    }
    const Real hckp = x.h*x.c*op(m,M1_OP_P,k,j,i)*gm;
    x.w(m,CLW_GM,k,j,i) = gm;
    x.w(m,CLW_EP,k,j,i) = -sb;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real e, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
      const Real s = 4.0*hckp*arad_*t*t*t;
      rr = (s > 0.0) ? s/(d*cv + s) : 0.0;   // (no exchange: none, also where T = 0)
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
                         - CleanSourceDivergence(x, m, k, j, i) + sb;
    x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanNewtonUpdate
//! \brief new iterate E_k = X and the gas temperature in exchange equilibrium with it
//! (the absorption takes max(X, 0); the staggered layout also keeps that as the iterate,
//! the cell layout keeps X, whose E < 0 its positivity limiter treats).  The cell layout
//! takes the gas eint change -A of eq. 56 of Ma et al. (2025) (CLW_GM, CLW_EP of
//! CleanCellRightHandSide, CleanCellNewtonUpdate).
//! Returns max |dE|/E, |dT|/T (clean_nl_norm = cell: |dE| over max(|E|, the E scale
//! cln_nl_escale)).

Real RadiationM1::CleanNewtonUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const Real hc = dt_sub*c_light, arad_ = arad;
  const bool ncell = cln_nl_cell;
  const Real esc = fmax(cln_nl_escale, 1.0e-300);
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  const bool cell = cln_cell, poslim = cln_poslim;

  Real change = 0.0, femax = 0.0;
  Kokkos::parallel_reduce("cln_newton", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dmax, Real &fe_max) {
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
    const Real e = (cell || poslim) ? x : ex, g = w(m,CLW_GM,k,j,i);
    const Real t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                       hc*op(m,M1_OP_P,k,j,i)*g*arad_,
                       w(m,CLW_EI,k,j,i) + hc*op(m,M1_OP_E,k,j,i)*ex
                       + w(m,CLW_EP,k,j,i), t_old);
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_T,k,j,i) = t;
    const Real de = ncell ? fabs(e - e_old)/fmax(fmax(fabs(e), fabs(e_old)), esc)
                          : fabs(e - e_old)/fmax(fmax(fabs(e), fabs(e_old)), 1.0e-300);
    const Real dt = fabs(t - t_old)/fmax(fmax(t, t_old), 1.0e-300);
    dmax = fmax(dmax, fmax(de, dt));
    // staggered: round-off level of E (its terms over E; cell: CleanCellVelocityUpdate)
    if (!cell) fe_max = fmax(fe_max, (fabs(e)*(1.0 + hc*op(m,M1_OP_E,k,j,i)) +
        fabs(w(m,CLW_EO,k,j,i)) + hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(t)))/
        fmax(fabs(e), 1.0e-300));
  }, Kokkos::Max<Real>(change), Kokkos::Max<Real>(femax));
  if (!cell) cln_fe = std::numeric_limits<Real>::epsilon()*CleanReduce(femax, true);
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
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;
  if (cln_fwork && cln_fsolve && cln_split) {   // (clean_stag_fwork: the momentum of the
    CleanStagFSolve(CLW_EK);                     //  F solve, which the write-back gives
    xh.fs = true;                                //  the gas, for the work W of the E row)
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto op = opac;
  const Real arad_ = arad;
  const bool fref = (force_ref == M1_FREF_CLEAN);
  const bool cap = fref && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;

  Real dvmax = 0.0, vmax = 0.0, fvmax = 0.0;
  Kokkos::parallel_reduce("cln_velocity", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
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
        Real mag = 0.0;
        dm = CleanMomentumChange(xh, CLW_EK, dir, m, k, j, i, &mag);
        // (force_reference = clean: the first coupled solve's force is the reference,
        // withheld from the gas; nothing while the gas momentum is frozen)
        if (cap) fr(m,dir,k,j,i) = dm/x.h;
        if (fref) dm -= cap ? dm : x.h*fr(m,dir,k,j,i);
        if (x.vzero) dm = 0.0;
        // (times that of E in the stress, as CleanCellVelocityUpdate)
        const Real e = fmax(fabs(x.w(m,CLW_EK,k,j,i)), 1.0e-300), hc = x.h*x.c;
        fv_max = fmax(fv_max, mag/d*(2.0 + (e*hc*op(m,M1_OP_E,k,j,i) + fabs(x.w(m,
                 CLW_EO,k,j,i)) + hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,
                 i))))/e));
      }
      dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+dir,k,j,i))/d);
      x.w(m,CLW_DM+dir,k,j,i) = dm;
      v2 += SQR((uh(m,IM1+dir,k,j,i) + dm)/d);
    }
    v_max = fmax(v_max, std::sqrt(v2));
  }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  CleanAdvectiveCoefficients();
  if (cap) cln_fref_set = true;
  vmax = fmax(CleanReduce(vmax, true), 1.0e-300);
  cln_fv = std::numeric_limits<Real>::epsilon()*CleanReduce(fvmax, true)/vmax;
  return CleanReduce(dvmax, true)/vmax;
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
  // (staggered clean_stag_tau_pos: E < 0 within tau_pos max E is kept, not counted)
  const Real zt = cln_cell ? 0.0 : cln_taupos*cln_emax;
  Real nneg = 0.0;
  Kokkos::parallel_reduce("cln_neg", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    Real &e = w(m,CLW_X,k+ks,j+js,i+is);
    if (e < 0.0 && e >= -zero) e = 0.0;
    if (!(e >= -zt)) n += 1.0;
  }, Kokkos::Sum<Real>(nneg));
  return static_cast<int>(CleanReduce(nneg, false));
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanStagRealizable
//! \brief staggered layout: the realizability of the solution X beyond the solver
//! tolerance (10 clean_tol max|B|, CleanCountNegative), R1: cells with E < 0 (cln_nr1),
//! R2: interior faces whose flux leaves more than c E of the donor cell (cln_nr2).
//! mark = 1 makes the faces of the R1 cells and the R2 faces low-order, mark = 2 every
//! interior face: w = 1, the rays' half-range flux c (h+ E_L + h- E_R) (until the next
//! CleanFaceCoefficients).  With every face low-order a backward-Euler solve has an
//! M-matrix (E >= 0) and |F_f| <= c E_donor (|h+-| <= 1).  Returns cln_nr1 + cln_nr2.

int RadiationM1::CleanStagRealizable(const int mark) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanArgs x = CleanMakeArgs();
  Real bmax = 0.0, emax = 0.0, emin = 0.0, eneg = 0.0, esum = 0.0;
  Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max, Real &e_max, Real &e_min, Real &e_neg,
                Real &e_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(x.w(m,CLW_B,k+ks,j+js,i+is)));
    const Real e = x.w(m,CLW_X,k+ks,j+js,i+is);
    e_max = fmax(e_max, e);
    e_min = fmin(e_min, e);
    e_neg += fmin(e, 0.0);
    e_sum += fabs(e);
  }, Kokkos::Max<Real>(bmax), Kokkos::Max<Real>(emax), Kokkos::Min<Real>(emin),
  Kokkos::Sum<Real>(eneg), Kokkos::Sum<Real>(esum));
  emax = CleanReduce(emax, true);
  cln_emax = emax;
  cln_eneg[0] = -CleanReduce(-emin, true)/fmax(emax, 1.0e-300);   // min E/max E
  cln_eneg[1] = CleanReduce(eneg, false)/fmax(CleanReduce(esum, false), 1.0e-300);
  // (clean_stag_tau_pos: R1, R2 within tau_pos max E)
  const Real zero = fmax(10.0*cln_tol*CleanReduce(bmax, true), cln_taupos*emax);
  Real n1 = 0.0, n2 = 0.0, n3 = 0.0, n5 = 0.0;
  cln_fh = 0.0;
  // clean_stag_bclip: the passes limit the bracket instead of making faces low-order
  // (1: clip to the donor bound at X; 2: an R1 cell's outgoing bracket scaled to what
  // it has, E_old + incoming bracket + emission, as CleanCellConsistency's FV rays)
  const bool clip = (mark == 1) && x.consist && (cln_bclip & 1);
  const bool fct = (mark == 1) && x.consist && (cln_bclip & 2);
  const bool keep = (cln_bclip & 4);   // (a low-order face keeps the limited bracket)
  auto lam = cln_v;
  if (fct) {
    const Real arad_ = arad;
    auto op = opac;
    par_for("cln_bfct", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is, is+nx1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real pout = 0.0, pin = 0.0;
      for (int d = 0; d < x.nf; ++d) {
        const Real cl = x.fc[d](m,CLF_RC,k,j,i);
        const Real ch = x.fc[d](m,CLF_RC,k+(d == 2),j+(d == 1),i+(d == 0));
        pout += x.h*x.idx[d]*(fmax(-cl, 0.0) + fmax(ch, 0.0));
        pin += x.h*x.idx[d]*(fmax(cl, 0.0) + fmax(-ch, 0.0));
      }
      const Real a = fmax(x.w(m,CLW_EO,k,j,i), 0.0) + pin + x.h*x.c*op(m,M1_OP_P,k,j,i)*
                     arad_*SQR(SQR(x.w(m,CLW_T,k,j,i)));
      lam(m,0,k,j,i) = (!(x.w(m,CLW_X,k,j,i) >= -zero) && pout > a) ? a/pout : 1.0;
    });
    CleanGhostFill(cln_v, 0, 1);
  }
  CleanArgs xh = x;
  xh.hi = true;
  const bool split = cln_split;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0, r5 = 0.0;
    Kokkos::parallel_reduce("cln_real", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &q1, Real &q2, Real &q3, Real &q4, Real &q5) {
      int m = idx/fkji;
      int k = (idx - m*fkji)/fji;
      int j = (idx - m*fkji - k*fji)/f1;
      int i = idx - m*fkji - k*fji - j*f1;
      const bool cell = (i < nx1 && j < nx2 && k < nx3);   // (its low face)
      k += ks;
      j += js;
      i += is;
      const Real el = x.w(m,CLW_X,k-dk,j-dj,i-di), er = x.w(m,CLW_X,k,j,i);
      if (dir == 0 && cell && !(er >= -zero)) q1 += 1.0;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const Real f = CleanFaceFlux(x, CLW_X, dir, m, k, j, i, true);
      const Real ed = (f > 0.0) ? el : er;
      const bool r2 = !(fabs(f) <= (1.0 + 1.0e-12)*x.c*ed + x.c*zero);
      if (r2) q2 += 1.0;
      if (split) {   // (diagnostic: the momentum state)
        const Real fh = CleanFaceFlux(xh, CLW_X, dir, m, k, j, i, true);
        const Real eh = (fh > 0.0) ? el : er;
        if (!(fabs(fh) <= (1.0 + 1.0e-12)*x.c*eh + x.c*zero)) q3 += 1.0;
        q4 = fmax(q4, fabs(fh)/(x.c*fmax(eh, 0.0) + x.c*zero + 1.0e-300));
      }
      if (clip || fct) {
        Real &cb = x.fc[dir](m,CLF_RC,k,j,i);
        const Real c0 = cb, p0 = f - cb;
        if (fct) cb *= (cb > 0.0) ? lam(m,0,k-dk,j-dj,i-di) : lam(m,0,k,j,i);
        if (clip) cb = fmin(fmax(cb, -x.c*fmax(er, 0.0) - p0), x.c*fmax(el, 0.0) - p0);
        if (cb != c0) q5 += 1.0;
      } else if (mark == 2 || (mark == 1 && (r2 || !(el >= -zero) || !(er >= -zero)))) {
        x.fc[dir](m,CLF_W,k,j,i) = 1.0;
        if (x.consist && !keep) x.fc[dir](m,CLF_RC,k,j,i) = 0.0;   // (no bracket)
      }
    }, Kokkos::Sum<Real>(r1), Kokkos::Sum<Real>(r2), Kokkos::Sum<Real>(r3),
    Kokkos::Max<Real>(r4), Kokkos::Sum<Real>(r5));
    n5 += r5;
    n1 += r1;
    n2 += r2;
    n3 += r3;
    cln_fh = fmax(cln_fh, r4);
  }
  cln_nr1 = static_cast<int>(CleanReduce(n1, false));
  cln_nr2 = static_cast<int>(CleanReduce(n2, false));
  cln_nr3 = static_cast<int>(CleanReduce(n3, false));
  cln_fh = CleanReduce(cln_fh, true);
  cln_nclip += static_cast<int>(CleanReduce(n5, false));
  return cln_nr1 + cln_nr2;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagBracket
//! \brief staggered clean_consistency (with clean_stag_split, Cartesian): the bracket
//! C_f = w (F_ray,f - Phi_f(J)) of every interior face (CLF_RC), F_ray the rays' energy
//! flux through it (CLN_FX, CLN_F2, CLN_F3 of its upper cell) and Phi_f(J) the energy
//! flux of the rays' own J, so that Phi(J) + C = (1 - w) Phi(J) + w F_ray.  Lagged (fixed
//! in the solve), energy flux only (CleanFaceSource); a low-order face drops it.

void RadiationM1::CleanStagBracket() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.consist = false;
  x.feedback = false;   // (the gas terms read the last iterate's work vectors)
  par_for("cln_rayj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x.w(m,CLW_Z,k,j,i) = x.cs(m,CLN_J,k,j,i);
  });
  CleanGhostFill(cln_w, CLW_Z, 1);
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int n = (dir == 0) ? CLN_FX : ((dir == 1) ? CLN_F2 : CLN_F3);
    par_for("cln_sbracket", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // (Phi(J) at its steady value, -(c/chi_f) div(D J) in the share 1 - w: the bracket,
      // and with it the steady state of the scheme, does not depend on h or on G)
      const Real w = x.fc[dir](m,CLF_W,k,j,i);
      const Real chi = 0.5*(x.op(m,M1_OP_T,k,j,i) + x.op(m,M1_OP_T,k-dk,j-dj,i-di));
      const Real cd = (chi > 0.0) ? x.c/chi : x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c;
      const Real pj = -(1.0 - w)*cd*CleanFaceStress(x, CLW_Z, dir, m, k, j, i, true) +
          w*x.c*(x.fc[dir](m,CLF_HP,k,j,i)*x.w(m,CLW_Z,k-dk,j-dj,i-di) +
                 x.fc[dir](m,CLF_HM,k,j,i)*x.w(m,CLW_Z,k,j,i));
      x.fc[dir](m,CLF_RC,k,j,i) = CleanBoundaryFace(x, dir, m, k, j, i) ? 0.0 :
          w*(x.c*x.cs(m,n,k,j,i) - pj);
    });
  }
  if (!(cln_bclip & 4)) return;
  // clean_stag_bclip & 4: every cell's outgoing bracket at most what it has, E_old +
  // h c kP a T^4 (FCT, factor of the donor): the right-hand side stays >= 0, so the
  // all-low-order backward-Euler solve keeps E >= 0 with the bracket
  const Real arad_ = arad;
  auto op = opac;
  auto lam = cln_v;
  const bool cln8 = (cln_bclip & 8);   // (budget with the incoming bracket)
  par_for("cln_bfct0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real pout = 0.0, pin = 0.0;
    for (int d = 0; d < x.nf; ++d) {
      const Real cl = x.fc[d](m,CLF_RC,k,j,i);
      const Real ch = x.fc[d](m,CLF_RC,k+(d == 2),j+(d == 1),i+(d == 0));
      pout += x.h*x.idx[d]*(fmax(-cl, 0.0) + fmax(ch, 0.0));
      pin += x.h*x.idx[d]*(fmax(cl, 0.0) + fmax(-ch, 0.0));
    }
    const Real a = fmax(x.w(m,CLW_EO,k,j,i), 0.0) + x.h*x.c*op(m,M1_OP_P,k,j,i)*arad_*
                   SQR(SQR(x.w(m,CLW_T,k,j,i))) + ((cln8) ? pin : 0.0);
    lam(m,0,k,j,i) = (pout > a) ? a/pout : 1.0;
  });
  CleanGhostFill(cln_v, 0, 1);
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("cln_bfct1", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real &cb = x.fc[dir](m,CLF_RC,k,j,i);
      cb *= (cb > 0.0) ? lam(m,0,k-dk,j-dj,i-di) : lam(m,0,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagFSolve
//! \brief clean_stag_fsolve (with clean_stag_split): with E of the solution, the face
//! momentum F on its own, one unknown per face, from
//!   (F - G)/h + c^2 [P_R - P_L]/dx_f = -c chi_f (F - A) (+ the beta terms),
//!   P_i = (D E)_i - alpha_i [F]_i/(2c),  [F]_i the jump of the two faces of cell i,
//! the cell layout's dissipation (fY = c (D E) - al (c/2) [y]) on the stress: per face
//!   F (1/theta + h c (al_L + al_R)/(2 dx_f)) - h c (al_L F_{f-1} + al_R F_{f+1})/(2 dx_f)
//!     = F^H/theta,
//! al_i = max(A(tau_i)^2, w of its faces), 1 where a face of the cell has
//! |F^H| > c E_donor.  al = 0 gives F^H (a face with al_L = al_R = 0 keeps it exactly).
//! The lines along each direction are tridiagonal (Thomas per MeshBlock line; the end
//! faces next to other MeshBlocks or a periodic wrap take the neighbours' F lagged,
//! iterated with the ghost exchange).  CLF_N holds F, CLF_N + 1 the dissipation term
//! (al_R [F]_R - al_L [F]_L)/(2 c dx_f) the gas momentum takes (CleanMomentumChange).
//! clean_stag_fkin: the streaming share w_i (C1 of the cell, the faces' ramp) of every
//! cell's momentum flux is the rays' kinetic upwind flux of the face momentum,
//!   P_i = (1 - w_i) [(D E)_i - al_i [F]_i/(2c)] + w_i (h+_i F_lo + h-_i F_hi)/c,
//! the staggered form of the cell layout's c (h+ y_L + h- y_R) (y = F/c; a beam, F = c E,
//! h+ = 1, gives D E = E), and on spherical polar the w share of the curvature terms of
//! D E is dropped (the cell layout's streaming share has the rays' own).  It telescopes,
//! the lines stay tridiagonal (M-matrix: h+ >= 0 >= h-), and w = al = 0 keeps F^H.

void RadiationM1::CleanStagFSolve(const int e) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  auto op = opac;
  auto al = cln_v;
  const bool fk = cln_fkin, dok = cln_dok;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    // F^H of every face (CLF_N) and al of every cell (cln_v slot 0)
    par_for("cln_fs0", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fc(m,CLF_N,k,j,i) = CleanFaceFlux(xh, e, dir, m, k, j, i, true);
      fc(m,CLF_N+1,k,j,i) = 0.0;
      if (!fk || CleanBoundaryFace(x, dir, m, k, j, i)) return;
      // (fkin: the w share of the transverse stresses D_{dir,b} E of F^H; on spherical
      // polar also that of its curvature terms)
      for (int b = 0; b < x.nf; ++b) {
        Real p[4], q[4], cp, cm;
        if (b == dir) continue;
        CleanStagEdge(x, e, dir, b, m, k, j, i, 1, c1lo, c1hi, p);
        CleanStagEdge(x, e, dir, b, m, k, j, i, -1, c1lo, c1hi, q);
        CleanStagEdgeWeights(x, dir, b, m, k, j, i, cp, cm);
        fc(m,CLF_N+1,k,j,i) += cp*p[0]*p[3] - cm*q[0]*q[3];
      }
      if (!x.sph) return;
      auto wq = [&](const int q, const int kk, const int jj, const int ii) {
        return CleanSmoothStep((x.cs(m,CLN_C1,kk,jj,ii) - c1lo)/(c1hi - c1lo))*
               CleanStagSph(x, q, dir, m, k, j, i)*x.w(m,e,kk,jj,ii);
      };
      fc(m,CLF_N+1,k,j,i) += 0.5*CleanStagSph(x, 1, dir, m, k, j, i)*
                             (wq(2, k-dk, j-dj, i-di) + wq(3, k, j, i));
    });
    par_for("cln_fs1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real rc;
      const Real tau = op(m,M1_OP_T,k,j,i)*CleanCellFaceDistance(x, dir, m, k, j, i);
      // (fkin: the streaming share is kinetic, so only the R2 faces need the dissipation)
      Real a = fk ? 0.0 : CleanDissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
      for (int s = 0; s < 2; ++s) {
        const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
        if (!fk) a = fmax(a, fc(m,CLF_W,kf,jf,iff));
        const Real f = fc(m,CLF_N,kf,jf,iff);
        const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
        if (!CleanBoundaryFace(x, dir, m, kf, jf, iff) &&
            !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
      }
      // (fkin: w of the cell, slot 2, and the part of P - D E without F, -w D E, slot 1)
      const Real wc = fk ? CleanSmoothStep((x.cs(m,CLN_C1,k,j,i) - c1lo)/(c1hi - c1lo))
                         : 0.0;
      al(m,0,k,j,i) = (1.0 - wc)*a;
      // (and the stress the kinetic form misses for the rays' own state, F = c s J:
      // (D - s^2) J, 0 for a beam, J/3 for an isotropic field, with J of the rays where
      // the cell holds at least that much, lagged)
      const Real dd = x.cs(m,CLN_D11+dir,k,j,i), ec = x.w(m,e,k,j,i);
      const Real sr = x.cs(m,CLN_HP+dir,k,j,i) + x.cs(m,CLN_HM+dir,k,j,i);
      al(m,1,k,j,i) = wc*((dd - sr*sr)*fmin(x.cs(m,CLN_J,k,j,i), fmax(ec, 0.0)) - dd*ec);
      al(m,2,k,j,i) = wc;
      // (F of the low face, for the exchange; with transverse terms the last F solve's)
      x.w(m,CLW_Y,k,j,i) = fc(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,j,i);
    });
    CleanGhostFill(cln_v, 0, 3);
    CleanGhostFill(cln_w, CLW_Y, 1);
    const int n = (dir == 0) ? ie - is + 1 : ((dir == 1) ? je - js + 1 : ke - ks + 1);
    const int na = (dir == 0) ? je - js + 1 : ie - is + 1;
    const int nb = (dir == 2) ? je - js + 1 : ke - ks + 1;
    for (int it = 0; it < 64; ++it) {
      Real dmax = 0.0, fmx = 0.0;
      if (fk && x.multi_d) {   // (the other rows' F of the last pass, transverse terms)
        Kokkos::deep_copy(Kokkos::subview(cln_w, Kokkos::ALL, CLW_Z, Kokkos::ALL,
                                          Kokkos::ALL, Kokkos::ALL),
                          Kokkos::subview(cln_w, Kokkos::ALL, CLW_Y, Kokkos::ALL,
                                          Kokkos::ALL, Kokkos::ALL));
      }
      Kokkos::parallel_reduce("cln_fs2", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1 + 1)*na*nb),
      KOKKOS_LAMBDA(const int idx, Real &d_max, Real &f_max) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
        // the line: faces t = 0..n along dir
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        Real cp = 0.0, dp = 0.0;
        for (int t = 0; t <= n; ++t) {
          int k, j, i;
          face(t, k, j, i);
          Real lo = 0.0, hi = 0.0, dg = 1.0, r = fc(m,CLF_N,k,j,i);
          const bool bd = CleanBoundaryFace(x, dir, m, k, j, i);
          const Real all = al(m,0,k-dk,j-dj,i-di), alr = al(m,0,k,j,i);
          const Real wl = al(m,2,k-dk,j-dj,i-di), wr = al(m,2,k,j,i);
          // (fkin: the rays' kinetic transport of F along b != dir through the dual faces
          // of CleanStagEdge, the other rows' F lagged; started from the last F solve's)
          Real tg = 0.0, tr = 0.0;
          for (int b = 0; b < x.nf && fk && !bd; ++b) {
            Real p[4], q[4], cp, cm;
            if (b == dir) continue;
            CleanStagEdge(x, e, dir, b, m, k, j, i, 1, c1lo, c1hi, p);
            CleanStagEdge(x, e, dir, b, m, k, j, i, -1, c1lo, c1hi, q);
            CleanStagEdgeWeights(x, dir, b, m, k, j, i, cp, cm);
            const int bk = (b == 2), bj = (b == 1), bi = (b == 0);
            tg += x.h*x.c*(cp*p[0]*p[1] - cm*q[0]*q[2]);
            tr -= x.h*x.c*(cp*p[0]*p[2]*x.w(m,CLW_Z,k+bk,j+bj,i+bi) -
                           cm*q[0]*q[1]*x.w(m,CLW_Z,k-bk,j-bj,i-bi));
          }
          if (!bd && (all > 0.0 || alr > 0.0 || wl > 0.0 || wr > 0.0 || tg > 0.0)) {
            const Real g = 0.5*x.h*x.c/CleanCellFaceDistance(x, dir, m, k, j, i);
            const Real it_ = 1.0/fc(m,CLF_TH,k,j,i);
            // (fkin: h c^2 [w P_kin]_L^R/dx, P_kin = (h+ F_lo + h- F_hi)/c of a cell)
            const Real q = x.h*x.c*(x.sph ? CleanStagSph(x, 0, dir, m, k, j, i)
                                          : x.idx[dir]);
            const Real hpl = x.cs(m,CLN_HP+dir,k-dk,j-dj,i-di);
            const Real hml = x.cs(m,CLN_HM+dir,k-dk,j-dj,i-di);
            const Real hpr = x.cs(m,CLN_HP+dir,k,j,i), hmr = x.cs(m,CLN_HM+dir,k,j,i);
            lo = -g*all - q*wl*hpl;
            hi = -g*alr + q*wr*hmr;
            dg = it_ - lo - hi + q*(wr*(hpr + hmr) - wl*(hpl + hml)) + tg;
            r = fc(m,CLF_N,k,j,i)*it_ + q*x.c*(al(m,1,k-dk,j-dj,i-di) - al(m,1,k,j,i)) +
                x.h*x.c*x.c*fc(m,CLF_N+1,k,j,i) + tr;
            // (the faces beyond a MeshBlock line end: the neighbours', lagged)
            if (t == 0) {r -= lo*x.w(m,CLW_Y,k-dk,j-dj,i-di); lo = 0.0;}
            if (t == n) {r -= hi*x.w(m,CLW_Y,k+dk,j+dj,i+di); hi = 0.0;}
          }
          const Real den = dg - lo*cp;
          cp = hi/den;
          dp = (r - lo*dp)/den;
          x.w(m,CLW_CP,k,j,i) = cp;
          x.w(m,CLW_R,k,j,i) = dp;
        }
        Real fn = 0.0;
        for (int t = n; t >= 0; --t) {
          int k, j, i;
          face(t, k, j, i);
          const Real fv = x.w(m,CLW_R,k,j,i) - ((t < n) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
          fn = fv;
          // (the last face's slot is a ghost: the neighbour's or a physical face)
          if (t < n) d_max = fmax(d_max, fabs(fv - x.w(m,CLW_Y,k,j,i)));
          f_max = fmax(f_max, fabs(fv));
          x.w(m,CLW_Y,k,j,i) = fv;
        }
      }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(fmx));
      CleanGhostFill(cln_w, CLW_Y, 1);
      if (CleanReduce(dmax, true) <= 1.0e-14*CleanReduce(fmx, true)) break;
    }
    // F and the dissipation term of every face
    par_for("cln_fs3", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);   // (F^H: the boundary's, kept)
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const Real f = x.w(m,CLW_Y,k,j,i);
      // (fkin, the solution's F solve: d = Phi^H - F, Phi^H = F^H + theta d the energy
      // flux's high-order part, whose memory is G + d (CleanFaceSource): its own history,
      // so that its steady state, F^H of G = F^H, does not depend on h)
      if (fk && e == CLW_X && dok) {
        fc(m,CLF_N+3,k,j,i) += fc(m,CLF_N,k,j,i) - f - (1.0 - fc(m,CLF_TH,k,j,i))*
                               fc(m,CLF_N+3,k,j,i);
      }
      fc(m,CLF_N,k,j,i) = f;
      fc(m,CLF_N+2,k,j,i) = f;
      const Real all = al(m,0,k-dk,j-dj,i-di), alr = al(m,0,k,j,i);
      const Real wl = al(m,2,k-dk,j-dj,i-di), wr = al(m,2,k,j,i);
      Real tk = 0.0;   // (fkin: the transverse kinetic stress)
      for (int b = 0; b < x.nf && fk; ++b) {
        Real p[4], q[4], cp, cm;
        if (b == dir) continue;
        CleanStagEdge(x, e, dir, b, m, k, j, i, 1, c1lo, c1hi, p);
        CleanStagEdge(x, e, dir, b, m, k, j, i, -1, c1lo, c1hi, q);
        CleanStagEdgeWeights(x, dir, b, m, k, j, i, cp, cm);
        const int bk = (b == 2), bj = (b == 1), bi = (b == 0);
        tk += (cp*p[0]*(p[1]*f + p[2]*x.w(m,CLW_Y,k+bk,j+bj,i+bi)) -
               cm*q[0]*(q[1]*x.w(m,CLW_Y,k-bk,j-bj,i-bi) + q[2]*f))/x.c;
      }
      fc(m,CLF_N+1,k,j,i) -= tk;
      if (!(all > 0.0 || alr > 0.0 || wl > 0.0 || wr > 0.0)) return;
      const Real fm = x.w(m,CLW_Y,k-dk,j-dj,i-di), fp = x.w(m,CLW_Y,k+dk,j+dj,i+di);
      auto pk = [&](const int kk, const int jj, const int ii, const Real l, const Real h) {
        return al(m,2,kk,jj,ii)*(x.cs(m,CLN_HP+dir,kk,jj,ii)*l +
               x.cs(m,CLN_HM+dir,kk,jj,ii)*h)/x.c + al(m,1,kk,jj,ii);
      };
      fc(m,CLF_N+1,k,j,i) += (alr*(fp - f) - all*(f - fm))/
                             (2.0*x.c*CleanCellFaceDistance(x, dir, m, k, j, i)) -
          (x.sph ? CleanStagSph(x, 0, dir, m, k, j, i) : x.idx[dir])*
          (pk(k, j, i, f, fp) - pk(k-dk, j-dj, i-di, fm, f));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagForceDiag
//! \brief diagnostic (clean_diag_force, staggered): on the low-order faces (w = 1) the
//! momentum a face hands to the gas, -[(F - G)/c^2 + h d_n(D E)] (CleanMomentumChange),
//! against the physical one, h chi_f (F - A)/c = (1/theta - 1)(F - A)/c^2 (A the face
//! advective flux).  Prints sum |mismatch| over those faces, the part with h c chi_f > 1,
//! both over sum |physical| of all interior faces, the largest |mismatch| over the
//! largest |physical|, and sum |mismatch| over sum |F|/c^2 of all interior faces.

void RadiationM1::CleanStagForceDiag() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;   // (the flux the gas sees)
  // phys, mis, thick, n, |F|/c^2, max mis, max phys
  Real r[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Kokkos::parallel_reduce("cln_fdiag", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &p, Real &mi, Real &th, Real &n, Real &fm,
                  Real &mm, Real &mp) {
      int m = idx/fkji;
      int k = (idx - m*fkji)/fji;
      int j = (idx - m*fkji - k*fji)/f1;
      int i = idx - m*fkji - k*fji - j*f1;
      k += ks;
      j += js;
      i += is;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const Real f = CleanFaceFlux(xh, CLW_X, dir, m, k, j, i, true), c2 = x.c*x.c;
      const Real t = x.fc[dir](m,CLF_TH,k,j,i);
      const Real a = x.feedback ? 0.5*(x.w(m,CLW_A1+dir,k-dk,j-dj,i-di)*
                     x.w(m,CLW_X,k-dk,j-dj,i-di) + x.w(m,CLW_A1+dir,k,j,i)*
                     x.w(m,CLW_X,k,j,i)) : 0.0;
      const Real ph = (1.0/t - 1.0)*(f - a)/c2;
      p += fabs(ph);
      fm += fabs(f)/c2;
      mp = fmax(mp, fabs(ph));
      if (x.fc[dir](m,CLF_W,k,j,i) < 1.0) return;
      const Real d = fabs(-(f - x.fc[dir](m,CLF_G,k,j,i))/c2 -
                          x.h*CleanFaceStress(x, CLW_X, dir, m, k, j, i, true) - ph);
      mi += d;
      if (t < 0.5) th += d;
      n += 1.0;
      mm = fmax(mm, d);
    }, Kokkos::Sum<Real>(q[0]), Kokkos::Sum<Real>(q[1]), Kokkos::Sum<Real>(q[2]),
    Kokkos::Sum<Real>(q[3]), Kokkos::Sum<Real>(q[4]), Kokkos::Max<Real>(q[5]),
    Kokkos::Max<Real>(q[6]));
    for (int n = 0; n < 5; ++n) r[n] += q[n];
    for (int n = 5; n < 7; ++n) r[n] = fmax(r[n], q[n]);
  }
  for (int n = 0; n < 7; ++n) r[n] = CleanReduce(r[n], n >= 5);
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN FORCE cycle %d t2s %d low %d mis %.3e thick %.3e max %.3e "
                "mom %.3e\n", pmy_pack->pmesh->ncycle, t2_solve, static_cast<int>(r[3]),
                r[1]/fmax(r[0], 1.0e-300), r[2]/fmax(r[0], 1.0e-300),
                r[5]/fmax(r[6], 1.0e-300), r[1]/fmax(r[4], 1.0e-300));
  }
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
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;   // (clean_stag_split: the faces store the momentum state)
  auto kk = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;
  if (cln_fsolve && cln_split) {   // the face F of the F solve and its dissipation
    CleanStagFSolve(CLW_X);
    xh.fs = true;
  }

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
      const Real f = xh.fs ? fc(m,CLF_N,k,j,i) : CleanFaceFlux(xh, CLW_X, dir, m, k, j, i, true);
      if (slopes) kk(m,tf,k,j,i) = (f - fc(m,CLF_G,k,j,i))/h;
      f0(m,k,j,i) = f;
    });
  }

  // E, cell fluxes, gas
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto u0_ = u0;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  const bool fref = (force_ref == M1_FREF_CLEAN);
  auto fr = cln_fref;
  const Real arad_ = arad;
  par_for("cln_wbcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = x.w(m,CLW_X,k,j,i);
    Real dm[3] = {0.0, 0.0, 0.0};
    Real den = 0.0;
    if (feedback) {
      dm[0] = CleanMomentumChange(xh, CLW_X, 0, m, k, j, i);
      if (multi_d) dm[1] = CleanMomentumChange(xh, CLW_X, 1, m, k, j, i);
      if (three_d) dm[2] = CleanMomentumChange(xh, CLW_X, 2, m, k, j, i);
      for (int n = 0; n < 3; ++n) {   // (force_reference = clean, frozen gas momentum)
        if (fref) dm[n] -= h*fr(m,n,k,j,i);
        if (x.vzero) dm[n] = 0.0;
      }
      // (the kinetic energy the gas gains, the work W of the E row at convergence)
      Real wk = 0.0;
      for (int n = 0; n < 3; ++n) {
        wk += (uh(m,IM1+n,k,j,i) + 0.5*dm[n])*dm[n]/uh(m,IDN,k,j,i);
      }
      for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += dm[n];
      // gas energy: exactly what the radiation gives up besides its flux divergence, so
      // that total energy changes only by the face fluxes
      Real div = x.idx[0]*(f1(m,k,j,i+1) - f1(m,k,j,i));
      if (multi_d) div += x.idx[1]*(f2(m,k,j+1,i) - f2(m,k,j,i));
      if (three_d) div += x.idx[2]*(f3(m,k+1,j,i) - f3(m,k,j,i));
      // (clean_stag_split: the faces store F^H, E moved by the energy flux)
      den = -(e - x.w(m,CLW_EO,k,j,i)) - ((xh.hi || x.sph || x.cub) ?
            CleanFluxDivergence(x, CLW_X, m, k, j, i, true) : h*div);
      if (x.vzero) {   // (steady start: the eint of the iterate's T, not the residual of
        Real ei, p, chi_rho, chi_t, cv;   //  a step h -> infinity, as the cell layout)
        eos.ThermoAt(uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), ei, p, chi_rho, chi_t, cv);
        den = ei - x.w(m,CLW_EI,k,j,i);
      } else {
        // the gas takes the exchange of the solution with the iterate's T, the beta terms
        // and W, E the rest (its fluxes): the linear solve's residual goes to E, not eint
        const Real t = x.w(m,CLW_T,k,j,i);
        const Real ex = h*x.c*(x.op(m,M1_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_*SQR(SQR(t)) -
                               x.op(m,M1_OP_E,k,j,i)*fmax(e, 0.0));
        // (the residual over the row's local diagonal: all of it in a scattering cell,
        // little where the exchange is stiff, whose X the solve fixes to tol/diagonal)
        const Real rd = (den - x.w(m,CLW_EP,k,j,i) + ex - wk)/
                        (1.0 + x.h*x.c*x.op(m,M1_OP_E,k,j,i)*(1.0 - x.w(m,CLW_RR,k,j,i)));
        e += rd;
        den -= rd;
      }
      uh(m,IEN,k,j,i) += den;
      x.w(m,CLW_JV,k,j,i) = den - wk;   // (its change of the gas eint, cln_neint)
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
    // (clean_stag_split: the cell flux is that of the energy flux Phi, the faces keep F^H)
    for (int d = 0; d < 3 && xh.hi; ++d) {
      if ((d == 1 && !multi_d) || (d == 2 && !three_d)) continue;
      u0_(m,M1_F1+d,k,j,i) = 0.5*(CleanFaceFlux(x, CLW_X, d, m, k, j, i, true) +
          CleanFaceFlux(x, CLW_X, d, m, k+(d == 2), j+(d == 1), i+(d == 0), true));
    }
  });
  // cells the write-back left with gas eint < 0 by lowering it (the solve is then not
  // accepted, as in the cell layout; an eint the hydro left <= 0 is not the solve's)
  cln_neint = 0;
  if (!feedback) return;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nx1 = id.nx1, nx2 = id.nx2, nx3 = id.nx3, nji = nx2*nx1, nkji = nx3*nji;
  Real nb = 0.0;
  Kokkos::parallel_reduce("cln_neint", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji + ks;
    const int j = (idx%nji)/nx1 + js, i = idx%nx1 + is;
    const Real ke = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                         SQR(uh(m,IM3,k,j,i)))/uh(m,IDN,k,j,i);
    if (!(CleanGasEint(x, uh, ke, m, k, j, i) > 0.0) && !(x.w(m,CLW_JV,k,j,i) >= 0.0) &&
        x.w(m,CLW_EI,k,j,i) > 0.0) n += 1.0;
  }, Kokkos::Sum<Real>(nb));
  cln_neint = static_cast<int>(CleanReduce(nb, false));
}

}  // namespace radm1
