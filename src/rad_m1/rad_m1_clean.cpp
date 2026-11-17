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
#include <utility>

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
// the policy of the heavy per-cell reductions: Kokkos' HIP backend compiles a plain
// RangePolicy reduction for up to 1024 threads per block, which caps it at 128 VGPRs on
// gfx942 and spilled them (cln_velocity 1005 VGPRs); 256 threads per block as athena's
// par_for and old production's m1-fast4.  Max reductions and counts: bitwise.
using CleanRp256 = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 1>>;
//----------------------------------------------------------------------------------------
//! \fn void CleanFatal
//! \brief stops the run with a message about an unsupported clean-scheme setup

static void CleanFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagGasMomentum
//! \brief staggered: the gas momentum of cell (k,j,i) in the components of its face
//! fluxes, those of the mesh; on the cubed sphere Cartesian (from the frame of y of
//! CleanGasMomentum)

KOKKOS_INLINE_FUNCTION
static void CleanStagGasMomentum(const CleanArgs &x, const DvceArray5D<Real> &uh,
                                 const int m, const int k, const int j, const int i,
                                 Real p[3]) {
  CleanGasMomentum(x, uh, m, k, j, i, p);
  if (!x.cub) return;
  Real q[3][3];
  const Real y[3] = {p[0], p[1], p[2]};
  CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) p[b] = y[0]*q[0][b] + y[1]*q[1][b] + y[2]*q[2][b];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagRefForce
//! \brief staggered: the reference force density force_reference withholds from the gas
//! of cell (k,j,i) (CleanRefForce of the cell layout) in the components of
//! CleanStagGasMomentum: clean the stored force, wb_arad d arad_ref along x1 (the radius
//! on the cubed sphere: no torque), wb_arad_md d times the stored acceleration (every
//! direction, cln_fref 0..2), else 0.  With wb_arad the gas takes the residual of
//! the radiation force, as the production kick (rad_m1.hpp:360-366, the residual dv of
//! ImplicitWorkRow, rad_m1_implicit.cpp:12353)

KOKKOS_INLINE_FUNCTION
static void CleanStagRefForce(const CleanArgs &x, const int fref,
                              const DvceArray5D<Real> &fr, const DvceArray4D<Real> &aref,
                              const Real d, const int m, const int k, const int j,
                              const int i, Real rf[3]) {
  Real q[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  if (fref == M1_FREF_WB_ARAD && x.cub) CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) {
    rf[b] = (fref == M1_FREF_CLEAN) ? fr(m,b,k,j,i) :
            ((fref == M1_FREF_WB_ARAD) ? d*aref(m,k,j,i)*q[0][b] :
            ((fref == M1_FREF_WB_MD) ? d*fr(m,b,k,j,i) : 0.0));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagMomentum
//! \brief staggered: the momentum the radiation gives to cell (k,j,i) over the stage
//! (CleanMomentumChange of each face family, mag its terms' size) in the components of
//! CleanStagGasMomentum.  On the cubed sphere the faces' changes (along their normals,
//! shared by wt) give the vector of CleanStagCsVector, of which the gas takes the r and
//! theta components; its phi component is the torque row of the cell layout's y_phi: the
//! change of L_z = V R F_phi/c^2 of the radiation's cell vector F (CleanStagCsVector of
//! the face momenta) and the torque h sum_f (M_f x (D E)_f n_f)_z of the stress through
//! the faces (M_f = int x dA, (D E)_f the mean of the two cells; none on physical
//! boundaries), so that gas plus radiation L_z changes only by round-off.

KOKKOS_INLINE_FUNCTION
static void CleanStagMomentum(const CleanArgs &x, const int e, const int m, const int k,
                              const int j, const int i, Real dm[3], Real mag[3]) {
  Real cf[6], wf[6];
  for (int d = 0; d < 3; ++d) {
    mag[d] = 0.0;
    dm[d] = (d < x.nf) ? CleanMomentumChange(x, e, d, m, k, j, i, &mag[d], &cf[2*d],
                                             &wf[2*d]) : 0.0;
  }
  if (!x.cub) return;
  Real dv[3], df[6], one[6], fd[3], q[3][3];
  CleanStagCsVector(x, m, k, j, i, cf, wf, dv);
  Real tz = 0.0;
  for (int d = 0; d < 3; ++d) {
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
      const Real f = x.fs ? x.fc[d](m,CLF_N,kf,jf,iff)
                          : CleanFaceFlux(x, e, d, m, kf, jf, iff, true);
      df[2*d+s] = f - x.fc[d](m,CLF_G,kf,jf,iff);
      one[2*d+s] = 1.0;
      if (CleanBoundaryFace(x, d, m, kf, jf, iff)) continue;
      Real pn[3], mf[3];
      for (int b = 0; b < 3; ++b) {
        pn[b] = 0.0;
        for (int c = 0; c < 3; ++c) {
          const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
          pn[b] += 0.5*(x.cs(m,n,kf-dk,jf-dj,iff-di)*x.w(m,e,kf-dk,jf-dj,iff-di) +
                        x.cs(m,n,kf,jf,iff)*x.w(m,e,kf,jf,iff))*
                   x.g(m,CLG_N+3*d+c,kf,jf,iff);
        }
        mf[b] = (d == 0) ? x.x1f(m,iff)*x.area[0](m,kf,jf,iff)*x.g(m,CLG_N+b,kf,jf,iff)
                         : x.g(m,CLG_M+3*(d-1)+b,kf,jf,iff);
      }
      tz += (s ? 1.0 : -1.0)*(mf[0]*pn[1] - mf[1]*pn[0]);
    }
  }
  CleanStagCsVector(x, m, k, j, i, df, one, fd);
  const Real rc = CleanCsFrame(x, m, k, j, i, q);
  const Real dr = q[0][0]*dv[0] + q[0][1]*dv[1] + q[0][2]*dv[2];
  const Real dt = q[1][0]*dv[0] + q[1][1]*dv[1] + q[1][2]*dv[2];
  const Real dp = -(q[2][0]*fd[0] + q[2][1]*fd[1] + q[2][2]*fd[2])/(x.c*x.c) -
                  x.h*tz/(x.vol(m,k,j,i)*rc);
  for (int b = 0; b < 3; ++b) dm[b] = dr*q[0][b] + dt*q[1][b] + dp*q[2][b];
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
    {"rad_m1/implicit_predictor", nullptr}, {"rad_m1/implicit_predictor_tcap", nullptr},
    {"rad_m1/implicit_resid_fatal", nullptr},
    {"rad_m1/implicit_pos_gas*", nullptr}, {"rad_m1/implicit_*", ""},
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
  // (staggered on the cubed sphere: no clean_stag_vimp yet, CleanStagDrag's half-cell
  // masses and D_nn would be the resampled gas ghosts' and the Cartesian D_dd at a seam)
  if (pmy_pack->pmesh->use_cubed_sphere && !cln_cell && cln_vimp) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> staggered on the cubed sphere: clean_stag_vimp = false"
                << std::endl;
    }
    cln_vimp = false;
  }
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb = pmy_pack->nmb_thispack;
  Kokkos::realloc(cln_s, nmb, CleanNs(cln_cell, pmy_pack->pmesh->use_cubed_sphere), n3,
                  n2, n1);
  // (staggered on the cubed sphere: the cell terms and faces of CleanStagFSolveCs)
  const bool stcs = pmy_pack->pmesh->use_cubed_sphere && !cln_cell;
  Kokkos::realloc(cln_v, nmb, stcs ? 18 : 3, n3, n2, n1);
  Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
  Kokkos::realloc(cln_w.o, nmb, (cln_cell ? CLW_NC : CLW_N) - CLW_NS, n3, n2, n1);
  cln_w.ns = CLW_NS;
  // the cubed sphere takes the half ranges about the face normals from the FV rays
  cln_hk = cln_consist || (pmy_pack->pmesh->use_cubed_sphere && vgd_on);
  const int nfc = std::max(cln_milne ? CLF_KK + 6 : (cln_hk ? CLF_HK + 2 : CLF_RC + 1),
                           cln_cell ? (cln_fsolve ? CLF_N + 4 : 0) :
                           ((cln_fsolve && cln_fkin) ? CLF_KD + 2 :
                           (cln_vimp ? CLF_DQ + 1 : CLF_QL + 3)));
  // the opacities depend on T: tables, or the T exponent (bump) of powerlaw (user)
  cln_opt = cln_cell && !opac_freeze && !opac_zero &&
            (opacity_type == M1_OPAC_TABLE ||
            (opacity_type != M1_OPAC_CONST && opac_b != 0.0));
  Kokkos::realloc(cln_fc1, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, nfc, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, nfc, n3, n2, n1);
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  Kokkos::deep_copy(cln_fc1, 0.0);
  Kokkos::deep_copy(cln_fc2, 0.0);
  Kokkos::deep_copy(cln_fc3, 0.0);
  for (int d = 0; d < 3 && cln_rstd.extent_int(0) == nmb && cln_fsolve; ++d) {
    auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);   // (from a restart)
    auto AL = Kokkos::ALL;
    Kokkos::deep_copy(Kokkos::subview(fc, AL, CLF_N + 3, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, d, AL, AL, AL));
    Kokkos::deep_copy(Kokkos::subview(fc, AL, CLF_KD, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, 3 + d, AL, AL, AL));
  }
  if (cln_cell) {
    const int nf = (cln_hk || cln_milne) ? CLC_NF : CLC_U;
    const int nc = (cln_pc != CLN_PC_X1LINE) ? CLC_N : (cln_pcfac ? nf : CLC_NL);
    Kokkos::realloc(cln_x, nmb, nc, n3, n2, n1);
    Kokkos::deep_copy(cln_x, 0.0);
    if (cln_yres == M1_YRES_LOCAL) {
      Kokkos::realloc(cln_dy, nmb, 3, n3, n2, n1);
      Kokkos::deep_copy(cln_dy, 0.0);
    }
  }
  if ((force_ref == M1_FREF_CLEAN || force_ref == M1_FREF_WB_MD) &&
      cln_fref.extent_int(0) != nmb) {   // (a restart has read it)
    Kokkos::realloc(cln_fref, nmb, (force_ref == M1_FREF_WB_MD) ? 9 : 3, n3, n2, n1);
  }
  if (cln_stst) Kokkos::realloc(cln_st, nmb, CLN_ST, indcs.nx3, indcs.nx2, indcs.nx1);
  // multigrid levels of the cell layout: halve every direction with an even number
  // (>= 4) of cells, down to the level where no direction can be halved (clean_pc =
  // prodmg, and staggered clean_stag_pc = mg: x1 is never halved, two coarse levels)
  if ((cln_cell && cln_pc != CLN_PC_X1LINE) || cln_stagmg) {
    const bool pm = (cln_cell && cln_pc == CLN_PC_PRODMG) || cln_stagmg;
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
      cln_mg.emplace_back("cln_mg", nmb, cln_stagmg ? CLS_N : CLM_N, nc[2], nc[1],
                          nc[0]);
    }
  }
  if (cln_extrap > 0) {
    // (a restart may have read t^n and the end of the last step)
    if (cln_mn.extent_int(0) != nmb) Kokkos::realloc(cln_mn, nmb, CLN_N, n3, n2, n1);
    if (cln_mc.extent_int(0) != nmb) Kokkos::realloc(cln_mc, nmb, CLN_N, n3, n2, n1);
    if (cln_mshare) {
      cln_mp = cln_mc;
    } else {
      Kokkos::realloc(cln_mp, nmb, CLN_N, n3, n2, n1);
    }
    if (cln_extrap == 2) Kokkos::realloc(cln_m1, nmb, CLN_N, n3, n2, n1);
  }
  Mesh *pm = pmy_pack->pmesh;
  if (pm->multilevel || chat != c_light) {
    CleanFatal("<rad_m1>/scheme = clean currently needs a uniform mesh and chat = c");
  }
  // the boundaries of the clean solve: open (bath, beam; outflow without a bath),
  // reflect, periodic and an imposed flux (x1; in the cell layout also x2 with the
  // problem generator's face fluxes).  Staggered: as the production row, where the
  // imposed flux is the face's whole energy flux and its face flux (rad_m1_implicit.cpp
  // ImplicitSolve, lines 17518, 18062 and 18638), the face takes it with no E weights
  // (CleanFaceWeights, CleanFaceSource), its F solve keeps it as the face F (no F
  // equation) and the inner half cell absorbs its momentum (CleanMomentumChange), as at
  // the open faces
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  for (int s = 0; s < 6; ++s) {
    const bool flux = bc[s] == M1_IBC_FLUX || (bc[s] == M1_IBC_X2FLUX && ifx2_on);
    if (!(bc[s] == M1_IBC_MARSHAK || bc[s] == M1_IBC_OUTFLOW || bc[s] == M1_IBC_REFLECT ||
          bc[s] == M1_IBC_PERIODIC || (flux && (cln_cell || s < 2)))) {
      CleanFatal("<rad_m1>/implicit_bc_* with scheme = clean: auto | marshak | outflow | "
                 "reflect | periodic, flux (x1; x2 with the problem generator's face "
                 "fluxes with clean_layout = cell); not efix");
    }
  }
  if (impl_lbath || iebath_hi_col_on) {
    CleanFatal("scheme = clean has no ghost-built inner face (problem/he_bc_linhao) and "
               "no per-column outer bath");
  }
  if (t2_theta > 0.0) {   // (the drop's re-entry and carry are the old write-back's)
    CleanFatal("<rad_m1>/time2_stage_local_be or time2_stage_theta with scheme = clean: "
               "its dropped stage increments are neither re-entered nor carried (energy "
               "not conserved)");
  }
  if ((force_ref == M1_FREF_WB_ARAD || force_ref == M1_FREF_WB_MD) && !fref_wsplit) {
    CleanFatal("<rad_m1>/force_reference = wb_arad(_md) with scheme = clean needs the "
               "reference work given by the problem generator (force_reference_work = "
               "split)");
  }
  if (force_ref == M1_FREF_WB_MD && (cln_cell || (!cln_steady && !cln_fref_set))) {
    CleanFatal("<rad_m1>/force_reference = wb_arad_md needs clean_layout = staggered and "
               "clean_init = steady (its reference is the steady start's force)");
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
  }
  const bool ring_ok = vgd_on && vra_on && g_vra_adv &&
                       (vgd_rot_dt > 0.0 || vgd_rot_auto) && vgd_rot_align &&
                       vra_k <= M1_RA_NF;
  if (cln_ring && !(ring_ok && vgd_xf == (cln_ring == 1))) {
    CleanFatal("<rad_m1>/clean_rot_ring needs vet_gd FV rays rotating in time: "
               "vet_rotate_avg = K with vet_rotate_advance, vet_gd_rotate_dt and "
               "_align, with vet_gd_rotate_xfade for live, without for jump");
  }
  // the sc rays are validated with the staggered layout only (the cell layout's ray face
  // sums of sc are the balance march: not conservative on the cubed sphere, L1(E) 0.4 on
  // the RAYSET2 sp decks)
  if (vgd_on && cln_rfv == 0 && cln_cell) {
    CleanFatal("<rad_m1>/clean_rays = sc needs clean_layout = staggered (the cell layout "
               "takes clean_rays = fv1)");
  }
  // the sc rays sweep one set at a time: a mean or a cross-fade of sets is the ring's
  if (vgd_on && VgdRotOn() && cln_rfv == 0 && (vra_on || vgd_xf) && !cln_ring) {
    CleanFatal("<rad_m1>/clean_rays = sc with vet_rotate_avg or vet_gd_rotate_xfade "
               "needs clean_rot_ring");
  }
  // jump: the first step of a block is a backward-Euler step (no slope of the old sets)
  if (cln_ring == 2) vgd_rot_slope = 2;
  if (pm->use_cubed_sphere) {
    // the rays are vet_gd's directions transported conservatively (CleanFvSweep) or
    // swept by vet_gd itself (clean_rays = sc: its band, VetGdCsSweep), or
    // the analytic sphere (eddington: of infinite radius); an even panel resolution
    // keeps the cell centres off the z axis (the frame of y, CleanCsFrame)
    if (!(cln_sphere || cln_edd || vgd_on) || (cln_consist && !vgd_on) ||
        vgd_twin || vgd_async || pm->mesh_indcs.nx2 % 2 != 0) {
      CleanFatal("<rad_m1>/scheme = clean on the cubed sphere needs "
                 "clean_closure = sphere | eddington (no clean_consistency) or closure = "
                 "vet_col with vet_gd (clean_rays = fv1 | fv | sc, no twin or async) and "
                 "an even <mesh>/nx2");
    }
    if (!vgd_on) {
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
        (cln_consist && !vgd_on) || vgd_twin || vgd_async) {
      CleanFatal("<rad_m1>/scheme = clean on spherical polar needs clean_layout = "
                 "cell, clean_closure = sphere (no clean_consistency) or closure = "
                 "vet_col with vet_gd (vet_gd_twin = false, no vet_gd_async)");
    }
    // the rays' theta walls follow the mesh: a reflecting wall mirrors (the conservative
    // transport reflects specularly, CleanFvWalls), a periodic one wraps
    if (vgd_on && !pm->use_polar_boundary) {
      vgd_thmir = (cln_rfv == 0) &&
                  (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::reflect);
      vgd_thvac = 0;
    }
    // (the sweep starts at the first shell, the constructor's vlat_taucut)
    // (implicit_flux_beam = halfrange has no effect here: no clean solve reads the
    // half-range ratios vgd_hr, so they are neither made nor checkpointed)
    impl_beam_hr = false;
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
//! exchange (a seam ghost is the neighbour's cell), then the centroid distances; so do
//! the owner's x2 and x3 normals, volume and MeshBlock gid (CLG_F), with which the
//! staggered layout finds a seam ghost's face families (CleanCsSlot).

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
  auto gid = pmy_pack->pmb->mb_gid.d_view;
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
      gn(m,CLG_F+c,k,j,i) = nv[0][c];
      gn(m,CLG_F+4+c,k,j,i) = nv[2][c];
    }
    gn(m,CLG_F+3,k,j,i) = vol(m,k,j,i);
    gn(m,CLG_F+7,k,j,i) = gid(m);
    gn(m,CLG_A,k,j,i) = 0.5*(rl*rl + rr*rr)*ss;
    gn(m,CLG_A+1,k,j,i) = 0.5*al*(th[0] + th[1]);
    gn(m,CLG_A+2,k,j,i) = 0.5*al*(th[2] + th[3]);
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
  CleanGhostFill(cln_gn, CLG_F, 11, false);
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
//! force_reference = clean or icr_fmeas.  The closure history is then dropped, so the
//! first step starts from this state alone.

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
  // force_reference = clean, or the problem generator's equilibrium IC (icr_fmeas,
  // he_star_m1 he_ic_equil), the gas internal energy relaxes with E (rho and momentum
  // frozen), so that the reference state is in radiative equilibrium; icr_fmeas also
  // keeps the x1 radiation force of each pass less the reference, dm - h rho a_ref
  // (CleanWriteBack), in icr_dm (icr_dm_dt = h), as the old module's force pass
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool relax = (force_ref == M1_FREF_CLEAN || icr_fmeas) && sfb && fl.on;
  if (relax && icr_fmeas) {
    Kokkos::realloc(icr_dm, nmb1 + 1, u0.extent(2), u0.extent(3), u0.extent(4));
  }
  auto uh = fl.u0;
  gas_feedback = relax;
  cln_vzero = relax;
  cln_mdcap = (force_ref == M1_FREF_WB_MD);   // (every pass's write-back: the last wins)
  t2_solve = M1_T2S_NONE;
  nsub = 1;
  bool big = true, ok = false;
  int n2 = 0;
  // (icr_fmeas: the large step is reached by a ladder x10 per pass from dt, its smaller
  // steps the thermal relaxation of the gas, whose solves fail when started at the top)
  const bool lad = relax && icr_fmeas;
  Real dp = 0.0, fac = lad ? 1.0 : cln_steady_dtfac, emax = 0.0;
  auto u0_ = u0;
  DvceArray5D<Real> uo("cln_steady_u", nmb1 + 1, 4, u0.extent(2), u0.extent(3),
                       u0.extent(4));
  DvceArray5D<Real> ug("cln_steady_g", relax ? nmb1 + 1 : 1, 4, u0.extent(2),
                       u0.extent(3), u0.extent(4));
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj*ni;
  auto &sz = pmy_pack->pmb->mb_size.h_view(0);
  const Real dv = sz.dx1*sz.dx2*sz.dx3;
  const bool sph = pm->use_spherical_polar || pm->use_cubed_sphere;   // (volume sums)
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
    if (icr_fmeas) icr_dm_dt = pm->dt;
    copy(true);
    cln_tc = -1.0e300;   // (predict: the closure of this pass, not of the last dt)
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
    if (d < cln_steady_tol && (!lad || fac >= cln_steady_dtfac)) big = false;
    if (lad) fac = fmin(10.0*fac, cln_steady_dtfac);
    ok = ok || !big;
    dp = d;
  }
  // restore the run's setup; the first step starts from this state without history
  sums(t1);
  if (cln_mdcap) {   // (wb_arad_md: the reference in the ghosts, the hydro walks read it)
    CleanGhostFill(cln_fref, 0, 9, true, 0);
    if (pm->use_spherical_polar && pm->multi_d) {   // (across a pole: no theta halves)
      auto fr = cln_fref;
      auto bc = pmy_pack->pmb->mb_bcs.d_view;
      const int n1 = fr.extent_int(4), n3 = fr.extent_int(2), ng = indcs.ng;
      par_for("cln_mdpole", DevExeSpace(), 0, nmb1, 0, n3-1, 0, ng-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int g, const int i) {
        for (int s = 0; s < 2; ++s) {
          if (bc(m,s ? BoundaryFace::outer_x2 : BoundaryFace::inner_x2) !=
              BoundaryFlag::polar) continue;
          const int j = s ? je + 1 + g : js - 1 - g;
          fr(m,5,k,j,i) = fr(m,6,k,j,i) = 0.0;
        }
      });
    }
    if (pm->use_cubed_sphere) {
      // (the cubed sphere: the lateral walks run along (e_xi, e_eta), not the face normals
      // of the shares, and a ghost across a panel seam holds the neighbour's axes: each
      // lateral half is a.(X_face - X_centre), the acceleration's line integral)
      auto fr = cln_fref;
      auto pnl = pmy_pack->pmb->mb_panel.d_view;
      auto x1v = pmy_pack->pcoord->x1v, x2v = pmy_pack->pcoord->x2v;
      auto x3v = pmy_pack->pcoord->x3v;
      auto x2f = pmy_pack->pcoord->xx2f, x3f = pmy_pack->pcoord->xx3f;
      const int n1 = fr.extent_int(4), n2 = fr.extent_int(3), n3 = fr.extent_int(2);
      par_for("cln_mdcs", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real c[3], f[3];
        const Real q = 0.25*M_PI;
        cubed_sphere::PanelToCart(pnl(m), q*x2v(m,j), q*x3v(m,k), c);
        for (int sd = 0; sd < 4; ++sd) {   // (xi -, +, eta -, +)
          const Real a2 = (sd < 2) ? x2f(m,j+sd) : x2v(m,j);
          const Real a3 = (sd < 2) ? x3v(m,k) : x3f(m,k+sd-2);
          cubed_sphere::PanelToCart(pnl(m), q*a2, q*a3, f);
          Real p = 0.0;
          for (int b = 0; b < 3; ++b) p += fr(m,b,k,j,i)*x1v(m,i)*(f[b] - c[b]);
          fr(m,5+sd,k,j,i) = (sd % 2) ? p : -p;
        }
      });
    }
    cln_fref_set = true;
    cln_mdcap = false;
  }
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
  // (the sc rays' ring on the cubed sphere: the union's directions and weights too)
  if (vgd_on && (cln_rfv > 0 || (cln_cssc && cln_ring)) && cln_fvn == 0) {
    CleanFvSets(pmy_pack->pmesh->time, pmy_pack->pmesh->ncycle - 1);
    // (the sc ring: its two sets, the intensities a restart file holds, or none from a
    // file of another ray method: zero, until the first build)
    if (cln_cssc) {
      CleanScUnion2();
      if (cln_fvi.extent_int(1) != cln_fvn) {
        Kokkos::realloc(cln_fvi, pmy_pack->nmb_thispack, cln_fvn, cln_s.extent_int(2),
                        cln_s.extent_int(3), cln_s.extent_int(4));
      }
    }
  }
  if (cln_cssc && !cln_ring) {   // (clean_rays = sc: the sweep's set, its cln_fvi)
    cln_fvd = vgd_dir;
    cln_fvn = vgd_n;
  }
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
  if (cln_ring == 2) CleanRingSolveSet();
  CleanOldVector();
  if (cell) {
    CleanCellOldVector();
    CleanCellFaceCoefficients();
  } else {
    CleanFaceCoefficients();
    if (cln_consist) CleanStagBracket();
    if (cln_gflux) CleanGhostFill(cln_w, CLW_EO, 1);   // (the donors' E_old, memory)
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
  bool all_low = false;
  cln_all_low = false;
  Real dpic = 0.0;
  bool converged = !iterate;
  const int maxnl = iterate ? cln_maxnl : 1;
  // the lagged slopes (and the staggered velocity) settle only to the accuracy of the
  // linear solve
  // (clean_conv = budget: of the pass, CleanNewtonUpdate; the first solve's from the
  // floor, and no state from the last solve, so that restarts stay bitwise)
  Real tol_nl = fmax(cln_tol_nl, 10.0*cln_tol);
  cln_tolb = tol_nl;
  // (the budget's norm and time error only in Newton passes of the gas coupling whose
  // K1 = f(U^n) is known; a single solve, a cold one without the time-error estimate,
  // or one without gas feedback, whose E row's residual is the energy error (no gas
  // row takes it), the strict test, whose answer is that of the converged step)
  cln_bon = cln_budget && iterate && cln_feedback && (t2_ok || cln_k1ok) &&
            t2_solve != M1_T2S_NONE && !cln_vzero &&
            t2k1.extent_int(0) >= pmy_pack->nmb_thispack;
  cln_ewr = -1.0;
  // a change that stopped decreasing within 100 times both its round-off level (cln_fe,
  // cln_fv, from the magnitudes of the terms of the solve) and tol_nl is converged too
  auto settled = [&](const Real d, const Real d0, const Real f) {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  };
  // <rad_m1>/dbg_rst_fp_cycle: fingerprints of the first solve of that cycle
  const bool fpc = (pmy_pack->pmesh->ncycle == dbg_rst_fp_cycle) && !dbg_rst_fp_done;
  if (fpc) RstFingerprint("cln_setup");
  if (cell) Kokkos::deep_copy(Kokkos::subview(cln_v, Kokkos::ALL, 2, Kokkos::ALL,
                                              Kokkos::ALL, Kokkos::ALL), 0.0);
  Real dnl = 0.0;   // (the last pass's change, implicit_resid_fatal)
  for (int ipos = 0; ; ++ipos) {
    CleanFirstIterate();
    cln_fsk = false;
    if (cell) CleanCellFirstIterate();
    if (cln_opt && cln_feedback) {
      CleanCellOpacity();
      CleanCellFaceSpeeds();
    }
    if (cln_feedback) CleanAdvectiveCoefficients();
    // (the first iterate's beta and emission for the drag's face flux, not the last
    // solve's: a restart, failed stage or discarded attempt does not carry those)
    if (cln_vimp && cln_feedback) CleanRightHandSide(true);
    if (cln_vimp && cln_feedback) CleanStagDrag(true);
    converged = !iterate;
    npic = 0;
    Real de0 = 1.0e300, dv0 = 1.0e300, trel = 0.0;
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
        // (clean_stag_sbimp: the first iterate's T in exchange with E_k and the beta
        // terms sb(E_k) the RHS just built, then the RHS of that T)
        if (inl == 0 && cln_sbimp && cln_feedback) {
          CleanNewtonUpdate();
          CleanRightHandSide();
        }
        CleanPreconditionerSetup();
      }
      if (fpc) RstFingerprint(("cln_pass" + std::to_string(inl) + "_pre").c_str());
      // (budget: a pass whose iterate the next pass takes the closure of is solved in
      // the strict norm, so that closure does not carry the inexact solve: with
      // clean_closure = extrap1/2 or iterate the step's answer depends on it)
      const bool cnext = (cln_closit && inl + 1 < cln_maxpic) || (clos1 && inl == 0);
      if (cnext) cln_ewr = -2.0;
      nit = CleanBiCGStab();
      if (cnext) cln_ewr = -1.0;
      if (fpc) RstFingerprint(("cln_pass" + std::to_string(inl) + "_post").c_str());
      // (clean_yres = local: the local exchange of the solution X, kept while a pass
      // leaves X as it is, so that the gas does not take the round-off of a new residual)
      if (cell && cln_yres != M1_YRES_OFF && cln_feedback &&
          (cln_yres == M1_YRES_POST || inl == 0 || nit != 0)) {
        CleanCellResidual();
      }
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (nit < 0) dnl = std::numeric_limits<Real>::infinity();
      if (!iterate || nit < 0) break;
      const Real de = cell ? CleanCellNewtonUpdate() : CleanNewtonUpdate();
      if (cln_bon) tol_nl = cln_tolb;
      if (lagged && inl < cln_maxpic) dpic = de;
      Real dv = 0.0;
      if (cln_feedback) dv = cell ? CleanCellVelocityUpdate() : CleanVelocityUpdate();
      dnl = fmax(de, dv);
      // (clean_tol_rel, as the old production's implicit_tol_rel: the strict test stops at
      // kappa times the first pass's change, at most clean_tol_rel_max)
      if (cln_tolrel > 0.0 && !cln_bon && inl == 0) trel = fmin(cln_tolrel*dnl, cln_tolrmx);
      tol_nl = fmax(tol_nl, trel);
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
    // (staggered: R1, R2 of the solution, counted)
    const int nbad = cell ? (cln_poslim ? CleanCellMarkFaces() : 0) :
                     CleanStagRealizable();
    if (!cell && ipos == 0) nfirst = {cln_nr1, cln_nr2};
    if (!cln_poslim || nbad == 0) break;
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
  // Euler through t2_fail below); without, as the old production's Picard solve, the run
  // stops where the last change (any after a failed linear solve) exceeds
  // implicit_resid_fatal, and goes on with the iterate below it
  if (nit < 0 && !t2stage) {
    const bool fatal = !dtr_on && impl_res_fatal > 0.0 && !(dnl <= impl_res_fatal);
    if (dtr_on) {
      if (!dtr_fail) {
        dtr_why = 1;
        dtr_resid = -1.0;
      }
      dtr_fail = true;
    }
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d: backward-Euler solve not converged (%s, %d "
                  "so far; last change %.3e)\n", pmy_pack->pmesh->ncycle,
                  dtr_on ? "dt_retry" : (fatal ? "no dt_retry: > implicit_resid_fatal" :
                  "no dt_retry: accepted"), static_cast<int>(cln_nfail), dnl);
    }
    if (fatal) CleanFatal("<rad_m1> clean: a backward-Euler solve did not converge and "
                          "its last change exceeds <rad_m1>/implicit_resid_fatal");
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
    }
  }

  // 5. new state: face fluxes, E, cell fluxes, gas, stage slopes
  if (cell) {
    CleanCellWriteBack();
  } else {
    // (a failed solve or a probe does not change the energy flux's memory)
    cln_dok = (nit >= 0) && !t2_inprobe;
    CleanWriteBack();
  }
  // a write-back that leaves gas eint < 0 (energy for the fluid floor to create) is
  // rejected like a failed solve (the stage is redone by backward Euler, a
  // backward-Euler step by dt_retry)
  if (cln_feedback && cln_neint > 0 && !t2stage && dtr_on) {
    if (!dtr_fail) {
      dtr_why = 4;   // (reason eint<=0, as the old production's)
      dtr_resid = cln_neint;
    }
    dtr_fail = true;
  }
  if (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2) {
    t2_fail = (nit < 0) || (nneg > 0) || (cln_feedback && cln_neint > 0);
    // time2_reject = retry: the cycle again at a smaller dt instead of the backward-Euler
    // redo; after dt_retry_t2_be such cuts the redo (as ImplicitSolve's t2retry-1020)
    const bool rt = t2_fail && t2_reject == 1 && dtr_on && !t2_inprobe && !t2_pfix_want;
    if (rt && !dtr_t2be) {
      if (!dtr_fail) {dtr_why = 5; dtr_resid = -1.0;}
      dtr_fail = true;
    }
    if (rt && dtr_t2be) dtr_t2hit = true;
  }
  // (K1 = f(U^n) for the next solve's predictor and time error: what an accepted
  // backward-Euler or stage-2 solve wrote; a failed one leaves none, the backward-Euler
  // redo and auto-BE steps of hesdirk2 keep the last)
  if (t2_solve == M1_T2S_BESTORE || t2_solve == M1_T2S_STAGE2) {
    cln_k1ok = (nit >= 0) && !dtr_fail && !(t2_solve == M1_T2S_STAGE2 && t2_fail);
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
      CleanGhostFill(cln_s, 0, cln_s.extent_int(1));
      CleanCellFaceSpeeds();
      if (cln_hk) CleanCellConsistency();
      if (cln_milne) CleanCellMilneSetup();
      CleanCellClosureSlope(1.0);
    }
  }
  if (fpc) dbg_rst_fp_done = true;
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
//! \fn void RadiationM1::CleanRstFaces
//! \brief a staggered restart from a file without the clean scheme's state (the old
//! production scheme): its faces hold the comoving normal flux F0 and its cells F =
//! mean(F0) + a E, a = v + v.D (ImplicitSolve's write-back, rad_m1_implicit.cpp:20287,
//! on the cubed sphere in covariant tangential components, :20282-20289), while the
//! staggered faces hold the lab-frame flux.  Every face not on a physical boundary takes
//! the old scheme's own face enthalpy flux a_f E_f (implicit_enthalpy, M1EnthCorr: plm,
//! central, or the donor cell for upwind, picked by the sign of a_f), with each cell's
//! a = (F - mean F0)/E from the file (its closure, no other enters); on the cubed sphere
//! a is Cartesian (exchanged across the seams, CleanStagN takes the face's component).
//! The faces of the physical boundaries keep theirs

void RadiationM1::CleanRstFaces() {
  if (!cln_rstold) return;
  if (!cln_init) CleanInit();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto u0_ = u0;
  auto w = cln_w;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  auto cw = pmy_pack->pcoord->cos_cell, sw = pmy_pack->pcoord->sin_cell;
  const Real efl = e_floor;
  const int enm = (dop_face == 1) ? M1_IENTH_CENTRAL : impl_enth;
  const bool fre = cln_freal;
  DvceArray1D<Real> rc("cln_rstr", 2);
  // a of every cell (CLW_DM, 3) and its E (CLW_B), as the old scheme's write-back
  par_for("cln_rsta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = fmax(u0_(m,M1_E,k,j,i), efl);
    Real g[3] = {0.5*(f1(m,k,j,i) + f1(m,k,j,i+1)), 0.0, 0.0};
    if (x.nf > 1) g[1] = 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i));
    if (x.nf > 2) g[2] = 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i));
    w(m,CLW_B,k,j,i) = e;
    for (int d = 0; d < 3; ++d) w(m,CLW_DM+d,k,j,i) = (u0_(m,M1_F1+d,k,j,i) - g[d])/e;
    if (!x.cub) return;
    // (covariant e_xi, e_eta components, those of the face normals g = sin (a^xi, a^eta))
    const Real c = cw(m,k,j), sn = sw(m,k,j), s2 = 1.0 - c*c;
    const Real a2 = (u0_(m,M1_F2,k,j,i) - c*u0_(m,M1_F3,k,j,i))/s2 - g[1]/sn;
    const Real a3 = (u0_(m,M1_F3,k,j,i) - c*u0_(m,M1_F2,k,j,i))/s2 - g[2]/sn;
    const Real a1 = u0_(m,M1_F1,k,j,i) - g[0];
    for (int b = 0; b < 3; ++b) {
      w(m,CLW_DM+b,k,j,i) = (a1*x.g(m,CLG_R+b,k,j,i) + a2*x.g(m,CLG_E+b,k,j,i) +
                             a3*x.g(m,CLG_E+3+b,k,j,i))/e;
    }
  });
  CleanGhostFill(cln_w, CLW_DM, 4);
  for (int dir = 0; dir < x.nf; ++dir) {
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("cln_rstf", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      // (cell o = -2, -1 | 0, 1 of the face: its a along the face normal, and E)
      auto a = [&](const int o) {
        return CleanStagN(x, CLW_DM, dir, m, k, j, i, k+o*dk, j+o*dj, i+o*di);
      };
      auto e = [&](const int o) {return w(m,CLW_B,k+o*dk,j+o*dj,i+o*di);};
      Real af = 0.5*(a(-1) + a(0)), ef = 0.5*(e(-1) + e(0)), dum;
      if (enm == M1_IENTH_UPWIND) {
        ef = (af > 0.0) ? e(-1) : e(0);
        af = (af > 0.0) ? a(-1) : a(0);
      } else if (enm >= M1_IENTH_PLM &&
                 !CleanBoundaryFace(x, dir, m, k-dk, j-dj, i-di) &&
                 !CleanBoundaryFace(x, dir, m, k+dk, j+dj, i+di)) {
        Real afl, afr;
        PLM(a(-2), a(-1), a(0), afl, dum);
        PLM(a(-1), a(0), a(1), dum, afr);
        af = 0.5*(afl + afr);
        if (af > 0.0) PLM(e(-2), e(-1), e(0), ef, dum);
        if (af < 0.0) PLM(e(-1), e(0), e(1), dum, ef);
      }
      f0(m,k,j,i) += af*ef;
      // (clean_stag_freal: realisable, at most c E of its donor.  In an E-floor hole of
      // the old scheme (its floor raised E after the write-back) the file's F and F0 are
      // not those of its E: the old scheme's next write-back would have clipped them
      // (M1ApplyLimits), so the excess is no momentum of its state; counted and its
      // momentum sum |dF| V/c^2 reported, not given to the floor gas)
      const Real fc = fmin(fmax(f0(m,k,j,i), -x.c*e(0)), x.c*e(-1));
      if (fre && fc != f0(m,k,j,i)) {
        const Real v = (x.sph || x.cub) ? x.vol(m,k,j,i) :
                       1.0/(x.idx[0]*x.idx[1]*x.idx[2]);
        Kokkos::atomic_add(&rc(0), 1.0);
        Kokkos::atomic_add(&rc(1), fabs(f0(m,k,j,i) - fc)*v/(x.c*x.c));
        f0(m,k,j,i) = fc;
      }
    });
  }
  if (fre) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rc);
    const Real n = CleanReduce(h(0), false), p = CleanReduce(h(1), false);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> clean_stag_freal RESTART CONVERSION (one time, not a "
                << "deposit): faces of the old scheme's file above c E_donor taken at it="
                << n << ", radiation momentum removed sum |dF| V/c^2=" << p
                << " (its energy equivalent c sum |dF| V/c^2=" << p*x.c
                << "; E unchanged)" << std::endl;
    }
  }
  // the cell flux of these faces, as CleanWriteBack stores it (on the cubed sphere the
  // face vector in the frame of y, not the old scheme's covariant panel components)
  par_for("cln_rstc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), (x.nf > 1) ? f2(m,k,j,i) : 0.0,
                       (x.nf > 1) ? f2(m,k,j+1,i) : 0.0, (x.nf > 2) ? f3(m,k,j,i) : 0.0,
                       (x.nf > 2) ? f3(m,k+1,j,i) : 0.0};
    Real fv[3] = {0.5*(c[0] + c[1]), 0.5*(c[2] + c[3]), 0.5*(c[4] + c[5])}, q[3][3];
    if (x.cub) {
      const Real wc[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      Real fy[3];
      CleanStagCsVector(x, m, k, j, i, c, wc, fy);
      CleanCsFrame(x, m, k, j, i, q);
      for (int d = 0; d < 3; ++d) fv[d] = q[d][0]*fy[0] + q[d][1]*fy[1] + q[d][2]*fy[2];
    }
    for (int d = 0; d < 3; ++d) u0_(m,M1_F1+d,k,j,i) = fv[d];
  });
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> clean restart from a file without the clean state: the faces "
              << "take its enthalpy flux a E (lab-frame F)" << std::endl;
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
  CleanArgs x = CleanMakeArgs(), xa = x;
  // (staggered on the cubed sphere: Cartesian, CleanStagGasMomentum)
  const bool sc = x.cub && !cln_cell;
  xa.cub = x.cub && !sc;

  par_for("cln_advcoef", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i);
    Real v[3], a[3];
    if (sc) {
      CleanStagGasMomentum(x, uh, m, k, j, i, v);
    } else {
      CleanGasMomentum(x, uh, m, k, j, i, v);
    }
    for (int n = 0; n < 3; ++n) v[n] = (v[n] + w(m,CLW_DM+n,k,j,i))/d;
    CleanAdvect(xa, m, k, j, i, v, a);
    for (int n = 0; n < 3; ++n) w(m,CLW_A1+n,k,j,i) = a[n];
  });
  // (staggered: their ghosts are first read after the next exchange, the RHS's)
  if (cln_cell) {
    CleanGhostFill(cln_w, CLW_A1, 3);
  } else {
    CleanGhostDefer(cln_w, CLW_A1, 3);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFirstIterate
//! \brief first Newton iterate: E of the old vector, no momentum change, and the gas
//! temperature in local exchange equilibrium with that E.  implicit_predictor (step, as
//! the old production): E + h K1_E instead, K1 the slope of U^n (hesdirk2, the last
//! solve's), where that moves the radiation temperature by less than the factor
//! implicit_predictor_tcap (staggered).  The starting point only: the solution does not
//! move

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
  auto k1 = t2k1;
  const bool pred = cln_kpred && !cln_cell && (t2_ok || cln_k1ok) &&
                    t2_solve != M1_T2S_NONE && !cln_vzero && t2k1.extent_int(0) > nmb1;
  const Real h = dt_sub, ecap = SQR(SQR(impl_pred_tcap)), rtol = cln_trtol;

  const Real tgf = cln_gtf ? eos.tfloor : 0.0;   // (clean_gas_tfloor)
  par_for("cln_iter0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real eo = fmax(w(m,CLW_EO,k,j,i), 0.0);
    const Real ep = pred ? eo + h*k1(m,M1_T2_E,k,j,i) : eo;
    const Real e = (ep > eo/ecap && ep < eo*ecap) ? ep : eo;
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) w(m,CLW_DM+n,k,j,i) = 0.0;
    if (feedback) {
      w(m,CLW_EP,k,j,i) = 0.0;   // (the beta terms of this T, CleanRightHandSide)
      w(m,CLW_T,k,j,i) = fmax(CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                             hc*op(m,M1_OP_P,k,j,i)*arad_,
                             w(m,CLW_EI,k,j,i) + hc*op(m,M1_OP_E,k,j,i)*e,
                             w(m,CLW_T,k,j,i), rtol), tgf);
    }
  });
  // (staggered: the ghosts of E_k for the gas terms, with the next exchange)
  if (feedback && !cln_cell) CleanGhostDefer(cln_w, CLW_EK, 1);
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagAdvFace
//! \brief the advective flux a E at face (k,j,i) normal to n, as the face fluxes carry
//! it: clean_stag_upwind the share phi as a_f E_f, donor cell plus the lagged rest,
//! CleanFaceWeights; the central mean there left phi (a_f E_f - mean) in F - a E, a
//! comoving flux the gas does not have, which 2 hc kE beta/c fed back lagged at every
//! front (clean_stag_gflux, CleanRightHandSide; xh the face momentum state)

KOKKOS_INLINE_FUNCTION
static Real CleanStagAdvFace(const CleanArgs &x, const CleanArgs &xh, const int n,
                             const int m, const int kk, const int jj, const int ii) {
  const int kl = kk - (n == 2), jl = jj - (n == 1), il = ii - (n == 0);
  const Real ul = x.w(m,CLW_A1+n,kl,jl,il), ur = x.w(m,CLW_A1+n,kk,jj,ii);
  const Real el = x.w(m,CLW_EK,kl,jl,il), er = x.w(m,CLW_EK,kk,jj,ii);
  const Real ac = 0.5*(ul*el + ur*er);
  if (!x.upw || x.cub || CleanBoundaryFace(x, n, m, kk, jj, ii)) return ac;
  Real al, ar, phi = 0.0;
  CleanFaceWeights(xh, n, m, kk, jj, ii, al, ar, &phi);
  if (!(phi > 0.0)) return ac;
  const Real ah = CleanStagAdvectCorr(x, n, m, kk, jj, ii) +
                  ((ul + ur > 0.0) ? ul*el : ur*er);
  return ac + phi*(ah - ac);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanRightHandSide
//! \brief B = E_old + h c kP a T^4 - r h c kE E_k - W - h div F_source + h esrc, with the
//! exchange linearised about the iterate (r = s/(d cv + s), s = 4 h c kP a T^3), W the
//! work of the radiation force of the iterate on the gas (with wb_arad and the reference
//! work h d arad_ref.v at the iterate's velocity, the row term of ImplicitWorkRow,
//! rad_m1_implicit.cpp:12438-12465), esrc the explicit energy source of <rad_m1> (the
//! production old-vector deposit, rad_m1_implicit.cpp:15335-15351), and the initial
//! guess X = E_k.  With gas coupling the beta terms of Ma et al. (2025, eq. 56; see
//! rad_m1_clean_cell.hpp):
//! the emission takes G = 1 - beta^2/2 (CLW_GM) and B adds the lagged A - hc (kP a T^4 G
//! - kE E) = hc kE [(beta^2/2 - beta.a/c) E_k + 2 beta.y_k] (y_k = F/c of the iterate,
//! mean of the two faces; the gas eint takes its opposite, CLW_EP), with beta = v/c of
//! the iterate (CLW_BT) and hc kP a T^4 (CLW_RE) for the face fluxes (eq. 57,
//! CleanFaceWeights).  clean_stag_sbimp takes (1 - r) sb(E) into the operator
//! (CleanApply, CleanStagBetaE) and B the same at E_k: the solution is that of the lagged
//! terms

void RadiationM1::CleanRightHandSide(const bool beta) {
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
  // (wb_arad: the reference work rides the row, as ImplicitWorkRow's force_reference_work
  // = split)
  const int fref = force_ref;
  const bool wsp = (fref == M1_FREF_WB_ARAD || fref == M1_FREF_WB_MD) && !x.vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  // (a backward-Euler step, FrefWaccOn: E pays exactly the reference work the gas got in
  // the hydro stages, its share h/dt of fref_wacc, as production's be write-back,
  // rad_m1_implicit.cpp m1_impl_fws_x)
  const bool wpa = wsp && FrefWaccOn();
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : dt_sub/pmy_pack->pmesh->dt;
  // (the explicit energy source esrc, problem/mlt_flux_frozen: h esrc joins the old E,
  // as ImplicitSolve's dt (chat/c) esrc; the slope of E carries it, CleanWriteBack)
  const bool eso = esrc_on;
  auto es = esrc;
  const bool bud = cln_bon;
  // (CleanSteadyStart's equilibrium relax, icr_fmeas: h -> infinity)
  const bool eqr = icr_fmeas && x.vzero && x.feedback && !x.sbi;

  if (x.feedback) {
    par_for("cln_beta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real p[3];
      CleanStagGasMomentum(x, uh, m, k, j, i, p);
      for (int n = 0; n < 3; ++n) {
        x.w(m,CLW_BT+n,k,j,i) = (p[n] + x.w(m,CLW_DM+n,k,j,i))/(uh(m,IDN,k,j,i)*x.c);
      }
      x.w(m,CLW_RE,k,j,i) = x.h*x.c*op(m,M1_OP_P,k,j,i)*arad_*
                            SQR(SQR(x.w(m,CLW_T,k,j,i)));
      // (r of the exchange, for the face beta terms: CleanFaceWeights, CleanFaceSource)
      const Real d = uh(m,IDN,k,j,i), t = x.w(m,CLW_T,k,j,i);
      Real b2 = 0.0, e, pg, chi_rho, chi_t, cv;
      for (int n = 0; n < 3; ++n) b2 += SQR(x.w(m,CLW_BT+n,k,j,i));
      eos.ThermoAt(d, t, e, pg, chi_rho, chi_t, cv);
      const Real s = 4.0*x.h*x.c*op(m,M1_OP_P,k,j,i)*(1.0 - 0.5*b2)*arad_*t*t*t;
      x.w(m,CLW_GX,k,j,i) = (s > 0.0) ? s/(d*cv + s) : 0.0;
    });
    CleanGhostFill(cln_w, CLW_BT, 4);   // (E_k's ghosts: those of its last change)
    CleanStagCsLateral(CLW_EK);
    CleanGhostFill(cln_w, CLW_GX, 1);
  }
  if (beta) return;
  CleanArgs xh = x;   // (the face momentum state, clean_stag_split)
  xh.hi = cln_split;
  xh.fs = cln_fsk;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const auto *peos = CleanDev(eos, 2);
  const Real tgf = cln_gtf ? eos.tfloor : 0.0;   // (clean_gas_tfloor)
  // each face's source once, for its two cells (CLF_SR), and with gflux its a E (CLF_PT)
  const int nf = x.nf;
  const bool af = x.feedback && x.gf && !x.cub;
  par_for("cln_rhsaf", DevExeSpace(), 0, nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    for (int n = 0; n < nf; ++n) {
      if (i > ie+(n == 0) || j > je+(n == 1) || k > ke+(n == 2)) continue;
      x.fc[n](m,CLF_SR,k,j,i) = CleanFaceSource(x, n, m, k, j, i);
      if (af) x.fc[n](m,CLF_PT,k,j,i) = CleanStagAdvFace(x, xh, n, m, k, j, i);
    }
  });
  par_for("cln_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    const auto &eos = *peos;
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hcke = x.h*x.c*op(m,M1_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0, gm = 1.0, sb = 0.0, dop = 0.0, ek = 0.0, gres = 0.0;
    Real eg = 0.0;
    bool gr = false;
    if (x.feedback) {
      Real b2 = 0.0, ba = 0.0, by = 0.0, bp = 0.0, y2 = 0.0, fs[3] = {0.0, 0.0, 0.0};
      // (clean_stag_fwork: the F solve's F, which the gas takes)
      auto ff = [&](const int d, const int kk, const int jj, const int ii) {
        return xh.fs ? xh.fc[d](m,CLF_N,kk,jj,ii)
                     : CleanFaceFluxCall(xh, CLW_EK, d, m, kk, jj, ii, true);
      };
      // (gflux: an open physical boundary face's flux, absorbed in an opaque half cell,
      // enters the cell mean as its mean over the half cell, (1 - e^-q)/q, q = chi dx/2,
      // as the momentum it gives, CleanMomentumChange; an imposed flux runs through it)
      Real fw[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      for (int n = 0; n < 3; ++n) {   // (trip 3: fw stays in registers)
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        for (int sd = 0; sd < 2 && x.gf && n < x.nf; ++sd) {
          if (!CleanBoundaryFace(x, n, m, k+sd*dk, j+sd*dj, i+sd*di) ||
              !CleanOpen(x, m, 2*n + sd)) continue;
          const Real q = 0.5*op(m,M1_OP_T,k,j,i)/CleanDivW(x, n, m, k, j, i, sd);
          fw[2*n+sd] = (q > 1.0e-8) ? -expm1(-q)/q : 1.0;
        }
      }
      for (int n = 0; n < 3; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        if (n < x.nf) fs[n] = fw[2*n]*ff(n, k, j, i) + fw[2*n+1]*ff(n, k+dk, j+dj, i+di);
      }
      if (x.cub) {   // (the cell's flux vector, CleanStagCsVector)
        Real c[6], w[6];
        for (int d = 0; d < 6; ++d) w[d] = 1.0;
        for (int d = 0; d < 3; ++d) {
          c[2*d] = fw[2*d]*ff(d, k, j, i);
          c[2*d+1] = fw[2*d+1]*ff(d, k+(d == 2), j+(d == 1), i+(d == 0));
        }
        CleanStagCsVector(x, m, k, j, i, c, w, fs);
        for (int d = 0; d < 3; ++d) fs[d] *= 2.0;
      }
      for (int n = 0; n < 3; ++n) {
        const Real bt = x.w(m,CLW_BT+n,k,j,i);
        if (n < x.nf) by += bt*0.5*fs[n]/x.c;
        y2 += SQR(0.5*fs[n]/x.c);
        // (gflux: a E face by face; the cubed sphere's face families are not the vector's
        // components, there the cell's a E, the form without gflux)
        if (n < x.nf && x.gf && !x.cub) {
          const int dk = (n == 2), dj = (n == 1), di = (n == 0);
          bp += bt*0.5*(x.fc[n](m,CLF_PT,k,j,i) + x.fc[n](m,CLF_PT,k+dk,j+dj,i+di))/x.c;
        }
        b2 += SQR(bt);
        ba += bt*x.w(m,CLW_A1+n,k,j,i)/x.c;
      }
      if (x.gf && x.cub) bp = ba*x.w(m,CLW_EK,k,j,i);
      // (y realisable, |y| <= E, as implicit_realisable_coupling: the gas eint takes the
      // same sb, CLW_EP, so the energy stays exact)
      ek = fmax(x.w(m,CLW_EK,k,j,i), 0.0);
      if (y2 > SQR(ek)) by *= ek/std::sqrt(y2);
      gm = 1.0 - 0.5*b2;
      sb = hcke*((0.5*b2 - ba)*x.w(m,CLW_EK,k,j,i) + 2.0*by);
      // gflux: sb = hcke [dop E + 2 beta.F'/c] with the comoving flux F' = F - a E taken
      // face by face (F - A of each face, ~ theta in an opaque cell), not as the faces'
      // mean F less the cell's a E: at a front the two means differ by the neighbour's
      // a E, which hcke >> 1 turned into a lagged source the iteration could not settle
      dop = 0.5*b2 + ba;
      if (x.gf) sb = hcke*(dop*x.w(m,CLW_EK,k,j,i) + 2.0*(by - bp));
      // (clean_stag_freal: the comoving part realisable, |E0 - E| <= 2 |beta| E as the
      // old production's implicit_realisable_coupling: a face F' of a hot neighbour
      // through a cold opaque cell is no flux of its E; the gas takes the same sb)
      if (x.fre) {
        const Real s0 = hcke*dop*x.w(m,CLW_EK,k,j,i), sl = 2.0*hcke*std::sqrt(b2)*ek;
        sb = s0 + fmin(fmax(sb - s0, -sl), sl);
      }
      if (x.gf) x.w(m,CLW_JV+1,k,j,i) = hcke*dop;   // (d sb/dE, CleanNewtonUpdate)
      gr = !x.sbi && sb != -x.w(m,CLW_EP,k,j,i);   // (T took the last pass's sb)
    }
    const Real hckp = x.h*x.c*op(m,M1_OP_P,k,j,i)*gm;
    x.w(m,CLW_GM,k,j,i) = gm;
    x.w(m,CLW_EP,k,j,i) = -sb;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, eg, p, chi_rho, chi_t, cv);
      const Real s = 4.0*hckp*arad_*t*t*t;
      rr = (s > 0.0) ? s/(d*cv + s) : 0.0;   // (no exchange: none, also where T = 0)
      // (1 - r = d cv/(d cv + s) without cancellation, as implicit_src_stable)
      x.w(m,CLW_RY,k,j,i) = (s > 0.0) ? d*cv/(d*cv + s) : 1.0;
      // the gas row's residual at the iterate with these sb (its T is the root with those
      // of the last pass), share r, as the old production's SRCR: a change of sb between
      // passes enters E with 1 - r, the exchange's share (radiation-dominated opaque
      // gas, r ~ 1, hcke beta F'/c >> eint: taken whole it contracted only ~0.5 per pass)
      if (gr) gres = rr*(x.w(m,CLW_EI,k,j,i) + hcke*ek - sb - eg - hckp*arad_*SQR(SQR(t)));
      // (clean_gas_tfloor: a gas at its floor T holds e(rho, tfloor), its row whole in E,
      // r = 1: the radiation takes exactly the energy the gas gives, CleanWriteBack)
      if (tgf > 0.0 && !(t > tgf)) {
        rr = 1.0;
        x.w(m,CLW_RY,k,j,i) = 0.0;
        gres = x.w(m,CLW_EI,k,j,i) + hcke*ek - sb - eg - hckp*arad_*SQR(SQR(t));
      }
      // gflux: the Doppler part hcke dop E of sb is a change of the absorption, hcke
      // (1 - dop): taken in the Newton step as the exchange (r -> r + dop (1 - r), sb
      // keeps its value at E_k), not lagged; where hcke dop >> 1 (an opaque cell with a
      // moving gas) the lagged term drove the iteration away
      if (x.gf) {
        const Real dc = fmin(fmax(dop, 0.0), 1.0);
        rr += dc*(1.0 - rr);
        x.w(m,CLW_RY,k,j,i) *= 1.0 - dc;
      }
      // work of the radiation force on the gas, (v_old + dm/(2 rho)).dm: exactly the
      // kinetic energy the gas gains, taken from the radiation; with wb_arad also the
      // work h rho arad_ref.v of the reference force at the iterate's velocity, which the
      // problem generator gives the gas
      Real pg[3], rf[3] = {0.0, 0.0, 0.0};
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      if (wsp) CleanStagRefForce(x, fref, fr, aref, d, m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (pg[n] + 0.5*dm)*dm/d;
        if (wsp && !wpa) work += x.h*rf[n]*(pg[n] + dm)/d;
      }
      if (wpa) work += fra*wc(m,k,j,i);
    }
    x.w(m,CLW_RR,k,j,i) = rr;
    // (eqr: the exchange through the gas row of the root T, hckp a T^4 - r hcke E =
    // EI - e(T) - sb + (1 - r) hcke E, without the cancellation of its two terms, each
    // ~ hcke E, whose round-off at h -> infinity exceeds E in opaque cells)
    x.w(m,CLW_B,k,j,i) = eqr ?
        x.w(m,CLW_EO,k,j,i) + x.w(m,CLW_EI,k,j,i) - eg - work
        + x.w(m,CLW_RY,k,j,i)*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0)
        - CleanSourceDivergence(x, m, k, j, i) + (eso ? x.h*es(m,k,j,i) : 0.0) :
        x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
        - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work
        - CleanSourceDivergence(x, m, k, j, i) + sb + gres + (eso ? x.h*es(m,k,j,i) : 0.0);
    // (clean_stag_sbimp: (1 - r) sb(E) in the operator as the cell layout's, Newton form)
    if (x.sbi) x.w(m,CLW_B,k,j,i) -= (1.0 - rr)*CleanStagBetaE(xh, CLW_EK, m, k, j, i);
    x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
    if (!bud) return;   // (clean_conv = budget: the gross terms of the E and gas rows)
    const Real xc = hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) + hckp*arad_*SQR(SQR(t)) +
                    fabs(work) + fabs(sb);
    // (a face's flux: the old vector's, at least what the cell's E can stream through it
    // in h, c theta E, theta = 1/(1 + h c chi): free streaming starts from F = 0; a
    // physical boundary face, whose old F and ghost opacity a fresh start has not set
    // yet, the cell's streaming alone)
    Real fl = eso ? fabs(es(m,k,j,i)) : 0.0;
    const Real eo = fabs(x.w(m,CLW_EO,k,j,i));
    const Real cs = x.c*eo/(1.0 + x.h*x.c*op(m,M1_OP_T,k,j,i));
    for (int n = 0; n < x.nf; ++n) {
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*(n == 2), jf = j + sd*(n == 1), iff = i + sd*(n == 0);
        const auto &f = x.fc[n];
        fl += CleanDivW(x, n, m, k, j, i, sd)*(CleanBoundaryFace(x, n, m, kf, jf, iff) ?
              cs : fmax(fabs(f(m,CLF_G,kf,jf,iff)), x.c*f(m,CLF_TH,kf,jf,iff)*eo));
      }
    }
    x.w(m,CLW_SE,k,j,i) = eo + xc + x.h*fl;
    x.w(m,CLW_SG,k,j,i) = fabs(x.w(m,CLW_EI,k,j,i)) + xc;
  });
  if (!bud) return;
  // (an empty cell, no energy, no flux: its rows at 1e-12 of the largest scale)
  auto w = cln_w;
  Real smax[2] = {0.0, 0.0};
  Kokkos::parallel_reduce("cln_smax", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &se, Real &sg) {
    se = fmax(se, w(m,CLW_SE,k,j,i));
    sg = fmax(sg, w(m,CLW_SG,k,j,i));
  }, Kokkos::Max<Real>(smax[0]), Kokkos::Max<Real>(smax[1]));
  const Real fse = 1.0e-12*CleanReduce(smax[0], true);
  const Real fsg = fmax(1.0e-12*CleanReduce(smax[1], true), 1.0e-300);
  par_for("cln_sfloor", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_SE,k,j,i) = fmax(w(m,CLW_SE,k,j,i), fmax(fse, 1.0e-300));
    w(m,CLW_SG,k,j,i) = fmax(w(m,CLW_SG,k,j,i), fsg);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanNewtonUpdate
//! \brief new iterate E_k = X and the gas temperature in exchange equilibrium with it
//! (the absorption takes max(X, 0); the staggered layout also keeps that as the iterate
//! without its positivity limiter and clean_stag_gflux, the cell layout keeps X, whose
//! E < 0 its positivity limiter treats; clean_stag_sbimp: the gas eint takes -sb(X),
//! CleanStagBetaE).  The cell layout
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
  const bool cell = cln_cell, poslim = cln_poslim || cln_gflux;
  const bool gfd = !cell && cln_gflux && !cln_sbimp && cln_feedback;
  // (clean_conv = budget: the changes over the rows' gross terms, and the time error h/2
  // |f(Y) - f(U^n)| of the iterate, Heun against Euler, K1 = f(U^n))
  const bool bud = cln_bon;
  const bool kok = bud && (t2_ok || cln_k1ok) && t2_solve != M1_T2S_NONE &&
                   !cln_vzero && t2k1.extent_int(0) >= nmb;
  const Real h = dt_sub, rtol = cln_trtol;
  auto k1 = t2k1;
  if (cln_sbimp && cln_feedback) {   // (clean_stag_sbimp: the gas takes -sb of X)
    CleanGhostFill(cln_w, CLW_X, 1);
    CleanArgs xh = CleanMakeArgs();
    xh.hi = cln_split;
    const CleanArgs *pxh = CleanDevArgs(xh, 1);
    par_for("cln_sbx", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is, is+nx1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &xh = *pxh;   // (in device memory: no per-thread copy)
      w(m,CLW_EP,k,j,i) -= CleanStagBetaE(xh, CLW_X, m, k, j, i) -
                           CleanStagBetaE(xh, CLW_EK, m, k, j, i);
    });
  }

  Real change = 0.0, femax = 0.0, ermax = 0.0;
  // (256 threads per block: the exchange temperature needs more than 128 registers)
  using Rp256 = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 1>>;
  const Real tgf = cln_gtf ? eos.tfloor : 0.0;   // (clean_gas_tfloor)
  Kokkos::parallel_reduce("cln_newton", Rp256(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dmax, Real &fe_max, Real &er_max) {
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
    // (gflux: the gas takes the Doppler part hc kE dop E of sb at X, as the E row takes
    // it in r; at E_k it lagged a source of hc kE dop, in a moving opaque cell more than
    // the net exchange (1 - r) hc kE, and T and E drifted apart pass by pass; sbimp:
    // the gas takes all of sb at X)
    if (gfd) w(m,CLW_EP,k,j,i) -= w(m,CLW_JV+1,k,j,i)*(x - e_old);
    const Real t = fmax(CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                       hc*op(m,M1_OP_P,k,j,i)*g*arad_,
                       w(m,CLW_EI,k,j,i) + hc*op(m,M1_OP_E,k,j,i)*ex
                       + w(m,CLW_EP,k,j,i), t_old, rtol), tgf);
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_T,k,j,i) = t;
    const Real de = ncell ? fabs(e - e_old)/fmax(fmax(fabs(e), fabs(e_old)), esc)
                          : fabs(e - e_old)/fmax(fmax(fabs(e), fabs(e_old)), 1.0e-300);
    const Real dt = fabs(t - t_old)/fmax(fmax(t, t_old), 1.0e-300);
    if (bud) {
      const Real d = uh(m,IDN,k,j,i), se = w(m,CLW_SE,k,j,i);
      Real e1, e0, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e1, p, chi_rho, chi_t, cv);
      eos.ThermoAt(d, t_old, e0, p, chi_rho, chi_t, cv);
      dmax = fmax(dmax, fmax(fabs(e - e_old)/se, fabs(e1 - e0)/w(m,CLW_SG,k,j,i)));
      if (kok) er_max = fmax(er_max, 0.5*fabs(e - w(m,CLW_EO,k,j,i) -
                                              h*k1(m,M1_T2_E,k,j,i))/se);
    } else {
      dmax = fmax(dmax, fmax(de, dt));
    }
    // staggered: round-off level of E (its terms over E; cell: CleanCellVelocityUpdate)
    if (!cell) fe_max = fmax(fe_max, (fabs(e)*(1.0 + hc*op(m,M1_OP_E,k,j,i)) +
        fabs(w(m,CLW_EO,k,j,i)) + hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(t)))/
        fmax(fabs(e), 1.0e-300));
  }, Kokkos::Max<Real>(change), Kokkos::Max<Real>(femax), Kokkos::Max<Real>(ermax));
  if (!cell) cln_fe = std::numeric_limits<Real>::epsilon()*CleanReduce(femax, true);
  if (bud) {   // (the tolerance of the next test and solve; for settled, the level of a
    // change the linear solve fixes at the floor, its 0.1 tfl in the budget's norm)
    const Real tfl = fmax(cln_tol_nl, 10.0*cln_tol);
    cln_tolb = kok ? fmin(fmax(cln_bkap*CleanReduce(ermax, true), tfl), cln_bmax) : tfl;
    cln_fe = 0.1*tfl;
  }
  // (clean_stag_sbimp, fwork: the face beta terms' emission of the new T, hc kP a T^4,
  // for the F solve of the velocity update, whose F the next pass's y and W take)
  if ((cln_sbimp || cln_fwork) && cln_feedback) {
    par_for("cln_sbre", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is,
            is+nx1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_RE,k,j,i) = hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(w(m,CLW_T,k,j,i)));
    });
    CleanGhostDefer(cln_w, CLW_RE, 1);
  }
  // (staggered: the ghosts of the new E_k, for the velocity update and the next RHS)
  if (!cell && cln_feedback) CleanGhostFill(cln_w, CLW_EK, 1);
  return CleanReduce(change, true);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanVelocityUpdate
//! \brief momentum change of the new iterate, less the reference force of
//! force_reference (CleanStagRefForce), and the advective coefficients of the new
//! velocity (implicit velocity).  Returns max |dv| over max |v|.

Real RadiationM1::CleanVelocityUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanStagCsLateral(CLW_EK);   // (E_k's ghosts: CleanNewtonUpdate)
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;
  xh.mcl = true;
  if (cln_fwork) {   // (clean_stag_fwork: the momentum of the
    CleanStagFSolve(CLW_EK);                     //  F solve, which the write-back gives
    xh.fs = true;                                //  the gas, for the work W of the E row)
    cln_fsk = true;
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto op = opac;
  const Real arad_ = arad;
  const int fref = force_ref;
  const bool cap = (fref == M1_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;

  // (with the F solve's F: the faces' terms and the drag first, once per face, then the
  // cells' momentum and their advective coefficients a_d, CleanAdvectiveCoefficients)
  const bool fa = xh.fs;
  if (fa) CleanStagDrag(false, true);
  xh.fa = fa;
  auto xa = x;
  xa.cub = false;
  Real dvmax = 0.0, vmax = 0.0, fvmax = 0.0;
  Kokkos::parallel_reduce("cln_velocity", CleanRp256(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real d = uh(m,IDN,k,j,i);
    Real v2 = 0.0, dmv[3], mg[3], pg[3], rf[3];
    CleanStagMomentum(xh, CLW_EK, m, k, j, i, dmv, mg);
    CleanStagGasMomentum(x, uh, m, k, j, i, pg);
    CleanStagRefForce(x, fref, fr, aref, d, m, k, j, i, rf);
    for (int dir = 0; dir < 3; ++dir) {
      Real dm = 0.0;
      if (dir == 0 || (dir == 1 && x.multi_d) || (dir == 2 && x.three_d)) {
        const Real mag = mg[dir];
        dm = dmv[dir];
        v_max = fmax(v_max, mag/d);   // (|v| floor: the size of the stage's force terms)
        // (force_reference = clean: the first coupled solve's force is the reference,
        // withheld from the gas, as is wb_arad's; nothing while the gas momentum is
        // frozen)
        if (cap) fr(m,dir,k,j,i) = dm/x.h;
        if (fref != M1_FREF_NONE) dm -= cap ? dm : x.h*rf[dir];
        if (x.vzero) dm = 0.0;
        // (times that of E in the stress, as CleanCellVelocityUpdate)
        const Real e = fmax(fabs(x.w(m,CLW_EK,k,j,i)), 1.0e-300), hc = x.h*x.c;
        fv_max = fmax(fv_max, mag/d*(2.0 + (e*hc*op(m,M1_OP_E,k,j,i) + fabs(x.w(m,
                 CLW_EO,k,j,i)) + hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,
                 i))))/e));
      }
      dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+dir,k,j,i))/d);
      x.w(m,CLW_DM+dir,k,j,i) = dm;
      v2 += SQR((pg[dir] + dm)/d);
    }
    v_max = fmax(v_max, std::sqrt(v2));
    if (!fa) return;
    Real v[3], a[3];
    for (int n = 0; n < 3; ++n) v[n] = (pg[n] + x.w(m,CLW_DM+n,k,j,i))/d;
    CleanAdvect(xa, m, k, j, i, v, a);
    for (int n = 0; n < 3; ++n) x.w(m,CLW_A1+n,k,j,i) = a[n];
  }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  if (fa) {
    CleanGhostDefer(cln_w, CLW_A1, 3);   // (as CleanAdvectiveCoefficients)
  } else {
    if (cln_vimp) CleanStagDrag(false);   // (with the coefficients of this dm)
    CleanAdvectiveCoefficients();
  }
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
//! \brief the number of cells whose E is negative beyond the solver tolerance (10 tol
//! max|B|); those within it CleanWriteBack raises to zero, counted

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
    if (!(w(m,CLW_X,k+ks,j+js,i+is) >= -zero)) n += 1.0;
  }, Kokkos::Sum<Real>(nneg));
  return static_cast<int>(CleanReduce(nneg, false));
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanStagRealizable
//! \brief staggered layout: the realizability of the solution X beyond the solver
//! tolerance (10 clean_tol max|B|, CleanCountNegative), R1: cells with E < 0 (cln_nr1),
//! R2: interior faces whose flux leaves more than c E of the donor cell (cln_nr2).
//! Returns cln_nr1 + cln_nr2.

int RadiationM1::CleanStagRealizable() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanStagCsLateral(CLW_X);
  CleanArgs x = CleanMakeArgs();
  Real bmax = 0.0;
  Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(x.w(m,CLW_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));
  const Real zero = 10.0*cln_tol*CleanReduce(bmax, true);
  Real n1 = 0.0, n2 = 0.0;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real r1 = 0.0, r2 = 0.0;
    Kokkos::parallel_reduce("cln_real", CleanRp256(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &q1, Real &q2) {
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
    }, Kokkos::Sum<Real>(r1), Kokkos::Sum<Real>(r2));
    n1 += r1;
    n2 += r2;
  }
  cln_nr1 = static_cast<int>(CleanReduce(n1, false));
  cln_nr2 = static_cast<int>(CleanReduce(n2, false));
  return cln_nr1 + cln_nr2;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagBracket
//! \brief staggered clean_consistency: the bracket C_f = w (F_ray,f - Phi_f(J)) of every
//! interior face, F_ray the rays' energy flux through it (CLN_FX, CLN_F2, CLN_F3 of its
//! upper cell) and Phi_f(J) the energy flux of the rays' own J, so that Phi(J) + C =
//! (1 - w) Phi(J) + w F_ray.  Lagged (fixed in the solve), energy flux only.  It enters
//! the streaming share as flux factors (clean_stag_gflux),
//! Phi_w = c (g+ E_L + g- E_R) with c (g+ J_L + g- J_R) w = c (h+ J_L + h- J_R) w + C_f
//! (CLF_HP, CLF_HM), so that the face flux of E = J is the bracket's, (1 - w) Phi(J) +
//! w F_ray: the half ranges H+-_f = h+-_f J of CleanFaceCoefficients, tilted by the
//! defect dd = C_f/(w c) in proportion to their sizes,
//!   H+-' = H+-_f (1 +- dd/s),  s = H+_f - H-_f,
//! the shrinking one at most to 0.  Linear in E with the upwind signs, so the streaming
//! share keeps the M-matrix for any E (an empty cell sends nothing, whatever the rays
//! hold; the partial-current form of HOLO methods).  g+- are not capped at 1: the sweep's
//! face flux can exceed c J of the donor at a beam edge (its foot interpolation reaches
//! past the two cells), and the cap there breaks the consistency of the beam.
//! On the cubed sphere the high x2, x3 faces of a MeshBlock take the rays' flux from its
//! last cell (CLN_FH; a seam ghost's low face is another face), as CleanCellConsistency.

void RadiationM1::CleanStagBracket() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.consist = false;
  x.feedback = false;   // (the gas terms read the last iterate's work vectors)
  const Real efl = e_floor;
  par_for("cln_rayj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x.w(m,CLW_Z,k,j,i) = x.cs(m,CLN_J,k,j,i);
  });
  CleanGhostFill(cln_w, CLW_Z, 1);
  CleanStagCsLateral(CLW_Z);
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
      Real gl = 0.0, gr = 0.0;
      CleanStagArea(x, dir, m, k, j, i, w, gl, gr);
      const Real pj = -(1.0 - w)*cd*CleanFaceStress(x, CLW_Z, dir, m, k, j, i, true) +
          w*x.c*(x.fc[dir](m,CLF_HP,k,j,i)*x.w(m,CLW_Z,k-dk,j-dj,i-di) +
                 x.fc[dir](m,CLF_HM,k,j,i)*x.w(m,CLW_Z,k,j,i)) +
          gl*x.w(m,CLW_Z,k-dk,j-dj,i-di) + gr*x.w(m,CLW_Z,k,j,i);
      // (cubed sphere: the high x2, x3 face of a MeshBlock from its last cell, CLN_FH)
      const bool up = x.cub && dir > 0 && (j > x.je || k > x.ke);
      const Real fy = up ? x.cs(m,CLN_FH+4*(dir-1),k-dk,j-dj,i-di) : x.cs(m,n,k,j,i);
      const bool bf = CleanBoundaryFace(x, dir, m, k, j, i);
      if (bf) return;
      // (the rays' J also where it is below clean_ray_min E: their face flux passes
      // through cells they hold almost nothing in, the side cells of a diagonal beam)
      const Real jl = x.cs(m,CLN_JR,k-dk,j-dj,i-di), jr = x.cs(m,CLN_JR,k,j,i);
      const Real hl = x.fc[dir](m,CLF_HP,k,j,i)*jl, hr = x.fc[dir](m,CLF_HM,k,j,i)*jr;
      const Real dd = fy - pj/x.c, fr = hl + hr + dd, s = hl - hr;
      // (H+', and H-' = fr - H+'; without half ranges the cell's h+- stay)
      const Real hp = !(s > 0.0) ? hl : ((dd < 0.0) ? fmax(hl + dd*hl/s, 0.0) :
                                        fr - fmin(hr - dd*hr/s, 0.0));
      // (rays' J at or below the radiation floor e_floor are no radiation field, J ~
      // 1e-100 in the E-floor holes of an evolved envelope: there the factors take the
      // realisable half range |H+-| <= J; unbounded they made rows of the operator ~1e100)
      const Real hpj = hp/jl, hmj = ((s > 0.0) ? fr - hp : hr)/jr;
      if (jl > 0.0) x.fc[dir](m,CLF_HP,k,j,i) = (jl > efl) ? hpj : fmin(hpj, 1.0);
      if (jr > 0.0) x.fc[dir](m,CLF_HM,k,j,i) = (jr > efl) ? hmj : fmax(hmj, -1.0);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanDMem
//! \brief clean_stag_fkin: the energy flux's memory d (CLF_N + 3) of the stage-loop start
//! into CLF_DS (save) or back (a rejected stage: its write-back, and that of a stage 1
//! before it, belong to a solution the backward-Euler redo discards)

void RadiationM1::CleanDMem(const bool save) {
  if (!clean || !cln_fsolve || !cln_fkin) return;
  if (!cln_init) CleanInit();
  auto AL = Kokkos::ALL;
  for (auto fc : {cln_fc1, cln_fc2, cln_fc3}) {
    auto d = Kokkos::subview(fc, AL, CLF_N + 3, AL, AL, AL);
    auto s = Kokkos::subview(fc, AL, CLF_DS, AL, AL, AL);
    if (save) {Kokkos::deep_copy(s, d);} else {Kokkos::deep_copy(d, s);}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanDInc
//! \brief clean_stag_fkin: d of a hesdirk2 stage's old vector, from the stage start state
//! (stage 2: w1 d^n + w0 Y1, as Time2FormStage) plus the stage increment c1 K1 + c2 K2 of
//! d's own slopes (CLF_KD, CLF_KD + 1; the write-back's), as every other component: the
//! energy flux's memory is part of the state the scheme integrates in time

void RadiationM1::CleanDInc(const bool two, const Real c1, const Real c2, const Real w1,
                            const Real w0) {
  if (!clean || !cln_fsolve || !cln_fkin) return;
  for (auto fc : {cln_fc1, cln_fc2, cln_fc3}) {
    par_for("cln_dinc", DevExeSpace(), 0, fc.extent_int(0) - 1, 0, fc.extent_int(2) - 1,
            0, fc.extent_int(3) - 1, 0, fc.extent_int(4) - 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real d = fc(m,CLF_N+3,k,j,i) + c1*fc(m,CLF_KD,k,j,i);
      if (two) d += w1*fc(m,CLF_DS,k,j,i) + (w0 - 1.0)*fc(m,CLF_N+3,k,j,i) +
                    c2*fc(m,CLF_KD+1,k,j,i);
      fc(m,CLF_N+3,k,j,i) = d;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanDtrD
//! \brief dt_retry: the memory d and its slope K1 (CLF_N + 3, CLF_KD) are the only face
//! values a cycle takes from the last one; the checkpoint keeps them in cln_rstd (save,
//! the restart's layout) instead of the face arrays

void RadiationM1::CleanDtrD(const bool save) {
  const int nmb = cln_fc1.extent_int(0);
  if (!clean || !cln_fsolve || !cln_fkin || nmb == 0) return;
  auto AL = Kokkos::ALL;
  if (save && (cln_rstd.extent_int(0) != nmb || cln_rstd.extent_int(1) != 6)) {
    Kokkos::realloc(cln_rstd, nmb, 6, cln_fc1.extent_int(2), cln_fc1.extent_int(3),
                    cln_fc1.extent_int(4));
  }
  for (int e = 0; e < cln_rstd.extent_int(1) && cln_rstd.extent_int(0) == nmb; ++e) {
    auto d = Kokkos::subview((e%3 == 0) ? cln_fc1 : ((e%3 == 1) ? cln_fc2 : cln_fc3), AL,
                             (e < 3) ? CLF_N + 3 : CLF_KD, AL, AL, AL);
    auto r = Kokkos::subview(cln_rstd, AL, e, AL, AL, AL);
    if (save) {Kokkos::deep_copy(r, d);} else {Kokkos::deep_copy(d, r);}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagFSolve
//! \brief clean_stag_fsolve (with clean_stag_split): with E of the solution, the face
//! momentum F on its own, one unknown per face, from
//!   (F - G)/h + c^2 [P_R - P_L]/dx_f = -c chi_f (F - A) (+ the beta terms),
//!   P_i = (D E)_i - alpha_i [F]_i/(2c),  [F]_i the jump of the two faces of cell i,
//! the cell layout's dissipation (fY = c (D E) - al (c/2) [y]) on the stress: per face
//!   F (1/theta + h c (al_L + al_R)/(2 dx_f))
//!     - h c (al_L F_{f-1} + al_R F_{f+1})/(2 dx_f)
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
//! h+ = 1, gives D E = E), on spherical polar in area form, (1/dx +- s) P_R,L (s of
//! CleanStagSph: the w share keeps 2 P_rr/r, F r^2 of a beam is steady).  It telescopes,
//! the lines stay tridiagonal (M-matrix: h+ >= 0 >= h-), and w = al = 0 keeps F^H.

void RadiationM1::CleanStagFSolve(const int e) {
  if (pmy_pack->pmesh->use_cubed_sphere) {
    CleanStagFSolveCs(e);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  auto op = opac;
  auto AL = Kokkos::ALL;
  const bool fk = cln_fkin, dok = cln_dok;
  const int kd = (t2_solve == M1_T2S_NONE) ? 0 :
                 ((t2_solve == M1_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD);
  const Real c1lo = cln_c0, c1hi = cln_c1;
  // (the F of every direction, CLW_Y, CLW_Z, CLW_S: one exchange)
  auto yfill = [&]() {
    if (x.nf == 3) CleanGhostDefer(cln_w, CLW_S, 1);
    CleanGhostFill(cln_w, CLW_Y, std::min(x.nf, 2));
  };
  // the directions in three phases (F^H, al; the line solves; F), each direction with
  // its own al (CLW_L + 3 dir) and F (cy), so that each exchange serves all of them
  for (int q = 0; q < 3*x.nf; ++q) {
    const int dir = q%x.nf, ph = q/x.nf;
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    auto al = Kokkos::subview(cln_w.o, AL, std::make_pair(CLW_L - CLW_NS + 3*dir,
                              CLW_L - CLW_NS + 3*dir + 3), AL, AL, AL);
    const int cy = (dir == 2) ? CLW_S : CLW_Y + dir;
    if (ph == 0) {
      // F^H of every face (CLF_N) and al of every cell (al slot 0)
      const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
      par_for("cln_fs0", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
        fc(m,CLF_N,k,j,i) = CleanFaceFlux(xh, e, dir, m, k, j, i, true);
        fc(m,CLF_N+1,k,j,i) = 0.0;
        if (!fk || !x.sph || dir == 2 || CleanBoundaryFace(x, dir, m, k, j, i)) return;
        // (spherical polar: the w share's transverse pressures, the curvature of D E but
        // that of the normal stress (q = 4), are the rays' own (CLN_G, the cell layout's
        // CleanCellSphSource) where the cell holds their J, faded as CleanRayShare)
        const Real q1 = CleanStagSph(x, 1, dir, m, k, j, i);
        const Real qn = 2.0*CleanStagSph(x, 4, dir, m, k, j, i)/q1;
        auto wq = [&](const int q, const int kk, const int jj, const int ii) {
          const Real ec = x.w(m,e,kk,jj,ii);
          const Real rs = CleanSmoothStep(2.5*fmax(ec, 0.0)/
                                          fmax(x.cs(m,CLN_J,kk,jj,ii), 1.0e-300) - 0.25);
          return CleanSmoothStep((x.cs(m,CLN_C1,kk,jj,ii) - c1lo)/(c1hi - c1lo))*rs*
                 (q1*(CleanStagSph(x, q, dir, m, k, j, i) -
                      qn*x.cs(m,CLN_D11+dir,kk,jj,ii))*ec -
                  x.cs(m,CLN_G+dir,kk,jj,ii)/x.vol(m,kk,jj,ii));
        };
        fc(m,CLF_N+1,k,j,i) += 0.5*(wq(2, k-dk, j-dj, i-di) + wq(3, k, j, i));
      });
      par_for("cln_fs1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rc;
        const Real tau = op(m,M1_OP_T,k,j,i)*CleanCellFaceDistance(x, dir, m, k, j, i);
        // (fkin: the streaming share is kinetic, without dissipation)
        Real a = fk ? 0.0 : CleanDissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
        for (int s = 0; s < 2; ++s) {
          const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
          if (!fk) a = fmax(a, fc(m,CLF_W,kf,jf,iff));
          const Real f = fc(m,CLF_N,kf,jf,iff);
          const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
          if (!fk && !CleanBoundaryFace(x, dir, m, kf, jf, iff) &&
              !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
        }
        // (fkin: w of the cell, slot 2, and the part of P - D E without F, -w D E,
        // slot 1)
        const Real wc = fk ? CleanSmoothStep((x.cs(m,CLN_C1,k,j,i) - c1lo)/(c1hi - c1lo))
                           : 0.0;
        al(m,0,k,j,i) = (1.0 - wc)*a;
        // (and the stress the kinetic form misses for the rays' own state, F = c s J:
        // (D - s^2) J, 0 for a beam, J/3 for an isotropic field, with J of the rays where
        // the cell holds at least that much, lagged)
        const Real dd = x.cs(m,CLN_D11+dir,k,j,i), ec = x.w(m,e,k,j,i);
        const Real sr = x.cs(m,CLN_HP+dir,k,j,i) + x.cs(m,CLN_HM+dir,k,j,i);
        al(m,1,k,j,i) = wc*((dd - sr*sr)*fmin(x.cs(m,CLN_J,k,j,i), fmax(ec, 0.0)) -
                            dd*ec);
        al(m,2,k,j,i) = wc;
        // (F of the low face, for the exchange; multi-D: the last F solve's)
        x.w(m,cy,k,j,i) = fc(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,j,i);
      });
      if (dir == x.nf - 1) {
        CleanGhostDefer(cln_w, CLW_L, 3*x.nf);
        yfill();
      }
    }
    const int n = (dir == 0) ? ie - is + 1 : ((dir == 1) ? je - js + 1 : ke - ks + 1);
    const int na = (dir == 0) ? je - js + 1 : ie - is + 1;
    const int nb = (dir == 2) ? je - js + 1 : ke - ks + 1;
    if (ph == 1) {
      // the lines' coefficients are the same in every pass below, only the neighbours'
      // F at the line ends changes: once, lo (CLW_AM), hi at the last face (CLW_AP), the
      // Thomas factors den (CLW_A0) and c' (CLW_CP), r (CLW_R0) and the faces with an
      // equation (CLW_P; scratch slots of the E solve)
      Kokkos::parallel_for("cln_fs2f", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                           (nmb1 + 1)*na*nb),
      KOKKOS_LAMBDA(const int idx) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        Real cp = 0.0;
        for (int t = 0; t <= n; ++t) {
          int k, j, i;
          face(t, k, j, i);
          Real lo = 0.0, hi = 0.0, dg = 1.0, r = fc(m,CLF_N,k,j,i);
          const bool bd = CleanBoundaryFace(x, dir, m, k, j, i);
          const Real all = al(m,0,k-dk,j-dj,i-di), alr = al(m,0,k,j,i);
          const Real wl = al(m,2,k-dk,j-dj,i-di), wr = al(m,2,k,j,i);
          const bool eq = !bd && (all > 0.0 || alr > 0.0 || wl > 0.0 || wr > 0.0);
          if (eq) {
            const Real g = 0.5*x.h*x.c/CleanCellFaceDistance(x, dir, m, k, j, i);
            const Real it_ = 1.0/fc(m,CLF_TH,k,j,i);
            // (fkin: h c^2 [w P_kin]_L^R/dx, P_kin = (h+ F_lo + h- F_hi)/c of a cell)
            // (spherical polar: area form, (1/dx -+ s) of cell L, R)
            const Real q = x.h*x.c*(x.sph ? CleanStagSph(x, 0, dir, m, k, j, i)
                                          : x.idx[dir]);
            const Real qs = x.sph ? x.h*x.c*CleanStagSph(x, 4, dir, m, k, j, i) : 0.0;
            const Real hpl = x.cs(m,CLN_HP+dir,k-dk,j-dj,i-di);
            const Real hml = x.cs(m,CLN_HM+dir,k-dk,j-dj,i-di);
            const Real hpr = x.cs(m,CLN_HP+dir,k,j,i), hmr = x.cs(m,CLN_HM+dir,k,j,i);
            lo = -g*all - q*wl*hpl + qs*wl*hpl;
            hi = -g*alr + q*wr*hmr + qs*wr*hmr;
            dg = it_ - lo - hi + q*(wr*(hpr + hmr) - wl*(hpl + hml)) +
                 qs*(wr*(hpr + hmr) + wl*(hpl + hml));
            r = fc(m,CLF_N,k,j,i)*it_ + q*x.c*(al(m,1,k-dk,j-dj,i-di) - al(m,1,k,j,i)) -
                qs*x.c*(al(m,1,k-dk,j-dj,i-di) + al(m,1,k,j,i)) +
                x.h*x.c*x.c*fc(m,CLF_N+1,k,j,i);
          }
          // (the faces beyond a MeshBlock line end: the neighbours', lagged)
          x.w(m,CLW_AM,k,j,i) = lo;
          x.w(m,CLW_AP,k,j,i) = hi;
          if (eq && t == 0) lo = 0.0;
          if (eq && t == n) hi = 0.0;
          const Real den = dg - lo*cp;
          cp = hi/den;
          x.w(m,CLW_A0,k,j,i) = den;
          x.w(m,CLW_CP,k,j,i) = cp;
          x.w(m,CLW_R0,k,j,i) = r;
          x.w(m,CLW_P,k,j,i) = eq ? 1.0 : 0.0;
        }
      });
      // one exact Thomas solve per line; the faces beyond a line end of a MeshBlock or
      // a periodic wrap are the last F solve's (lagged)
      Kokkos::parallel_for("cln_fs2", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                           (nmb1 + 1)*na*nb),
      KOKKOS_LAMBDA(const int idx) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        Real dp = 0.0;
        for (int t = 0; t <= n; ++t) {
          int k, j, i;
          face(t, k, j, i);
          Real r = x.w(m,CLW_R0,k,j,i);
          if (x.w(m,CLW_P,k,j,i) > 0.0) {
            if (t == 0) r -= x.w(m,CLW_AM,k,j,i)*x.w(m,cy,k-dk,j-dj,i-di);
            if (t == n) r -= x.w(m,CLW_AP,k,j,i)*x.w(m,cy,k+dk,j+dj,i+di);
          }
          dp = (r - ((t == 0) ? 0.0 : x.w(m,CLW_AM,k,j,i))*dp)/x.w(m,CLW_A0,k,j,i);
          x.w(m,CLW_R,k,j,i) = dp;
        }
        Real fn = 0.0;
        for (int t = n; t >= 0; --t) {   // (the last face's slot is a ghost: the
          int k, j, i;                   //  neighbour's or a physical face)
          face(t, k, j, i);
          fn = x.w(m,CLW_R,k,j,i) - ((t < n) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
          x.w(m,cy,k,j,i) = fn;
        }
      });
      if (dir == x.nf - 1) yfill();
    }
    if (ph == 2) {
      // F and the dissipation term of every face
      par_for("cln_fs3", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);   // (F^H: the boundary's, kept)
        if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
        const Real f = x.w(m,cy,k,j,i);
        // (fkin, the solution's F solve: d = Phi^H - F, Phi^H = F^H + theta d the energy
        // flux's high-order part, whose memory is G + d (CleanFaceSource): its own
        // history, so that its steady state, F^H of G = F^H, does not depend on h)
        if (fk && e == CLW_X && dok) {   // (and its slope, CleanDInc)
          const Real dd = fc(m,CLF_N,k,j,i) - f - (1.0 - fc(m,CLF_TH,k,j,i))*
                          fc(m,CLF_N+3,k,j,i);
          fc(m,CLF_N+3,k,j,i) += dd;
          if (kd > 0) fc(m,kd,k,j,i) = dd/x.h;
        }
        fc(m,CLF_N,k,j,i) = f;
        fc(m,CLF_N+2,k,j,i) = f;
        const Real all = al(m,0,k-dk,j-dj,i-di), alr = al(m,0,k,j,i);
        const Real wl = al(m,2,k-dk,j-dj,i-di), wr = al(m,2,k,j,i);
        if (!(all > 0.0 || alr > 0.0 || wl > 0.0 || wr > 0.0)) return;
        // (a physical high face keeps F^H; the ghost fill replaced its slot)
        const Real fm = x.w(m,cy,k-dk,j-dj,i-di);
        const Real fp = CleanBoundaryFace(x, dir, m, k+dk, j+dj, i+di) ?
                        fc(m,CLF_N,k+dk,j+dj,i+di) : x.w(m,cy,k+dk,j+dj,i+di);
        auto pk = [&](const int kk, const int jj, const int ii, const Real l,
                      const Real h) {
          return al(m,2,kk,jj,ii)*(x.cs(m,CLN_HP+dir,kk,jj,ii)*l +
                 x.cs(m,CLN_HM+dir,kk,jj,ii)*h)/x.c + al(m,1,kk,jj,ii);
        };
        const Real q = x.sph ? CleanStagSph(x, 0, dir, m, k, j, i) : x.idx[dir];
        const Real qs = x.sph ? CleanStagSph(x, 4, dir, m, k, j, i) : 0.0;
        fc(m,CLF_N+1,k,j,i) += (alr*(fp - f) - all*(f - fm))/
                               (2.0*x.c*CleanCellFaceDistance(x, dir, m, k, j, i)) -
            q*(pk(k, j, i, f, fp) - pk(k-dk, j-dj, i-di, fm, f)) -
            qs*(pk(k, j, i, f, fp) + pk(k-dk, j-dj, i-di, fm, f));
      });
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagFSolveCs
//! \brief clean_stag_fsolve on the cubed sphere: the F solve of CleanStagFSolve for the
//! face-normal momentum F_n.  The cell terms of a line (the dissipation weight al, the
//! kinetic share w (h+ F_lo + h- F_hi)/c and its lagged pressure (D - s^2) J) are those
//! of the cell's face family along the line, D taken about the family's mean normal
//! (CleanCsFamily).  The kinetic stress w (F_kin/c - s^2 J) of the family's normal, 0 for
//! an isotropic field, enters in area form, (A_R P_R - A_L P_L)/V_f, A the mean area of
//! the cell's two faces of the family, so that free streaming F keeps r^2 F; w D (J - E)
//! as a gradient (A_f/V_f).  A line that ends at another MeshBlock takes the far
//! face of the next cell from that cell's owner, lagged (cln_v, in the owner's axes,
//! CleanCsSlot: across a seam an x2 line may continue as an x3 line), so the lines of all
//! directions iterate together.  A face shared by two MeshBlocks then takes the solution
//! of the one with the lower gid: both see the same F to the bit.  The lateral stress
//! keeps D E (no transverse kinetic stress of clean_stag_fkin on the cubed sphere: its
//! lateral Gauss sums are also the non-orthogonal correction of the normal stress, which
//! a kinetic transport of F does not reproduce, and a static envelope would lose its
//! balance on the panel-tangent faces).

void RadiationM1::CleanStagFSolveCs(const int e) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanStagCsLateral(e);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  auto op = opac;
  auto v = cln_v;   // (cell terms of family q: 3 q + 0..2, 15 + q; F: 9 + 2 q + hi)
  const bool fk = cln_fkin, dok = cln_dok;
  const int kd = (t2_solve == M1_T2S_NONE) ? 0 :
                 ((t2_solve == M1_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD);
  const Real c1lo = cln_c0, c1hi = cln_c1;
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    const CleanArgs *pxh = CleanDevArgs(xh, 1);
    par_for("cln_fcs0", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &xh = *pxh;   // (in device memory: no per-thread copy)
      fc(m,CLF_N,k,j,i) = CleanFaceFlux(xh, e, dir, m, k, j, i, true);
      fc(m,CLF_N+1,k,j,i) = 0.0;
    });
  }
  // the cell terms of every family, started from the last F solve's faces
  const CleanArgs *px = CleanDevArgs(x, 0);
  par_for("cln_fcs1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
    Real nd[3][3];
    CleanCsFamily(x, m, k, j, i, nd);
    const Real ec = x.w(m,e,k,j,i);
    const Real wc = fk ? CleanSmoothStep((x.cs(m,CLN_C1,k,j,i) - c1lo)/(c1hi - c1lo))
                       : 0.0;
    for (int q = 0; q < 3; ++q) {
      const int dk = (q == 2), dj = (q == 1), di = (q == 0);
      Real rc;
      const Real tau = op(m,M1_OP_T,k,j,i)*CleanCellFaceDistance(x, q, m, k, j, i);
      Real a = fk ? 0.0 : CleanDissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
      for (int s = 0; s < 2; ++s) {
        const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
        if (!fk) a = fmax(a, x.fc[q](m,CLF_W,kf,jf,iff));
        const Real f = x.fc[q](m,CLF_N,kf,jf,iff);
        const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
        if (!fk && !CleanBoundaryFace(x, q, m, kf, jf, iff) &&
            !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
        v(m,9+2*q+s,k,j,i) = x.fc[q](m,CLF_N+2,kf,jf,iff);
      }
      Real dd = 0.0;
      for (int b = 0; b < 3; ++b) {
        for (int c = 0; c < 3; ++c) {
          dd += nd[q][b]*nd[q][c]*x.cs(m,CLN_D11+((b == c) ? b : 2 + b + c),k,j,i);
        }
      }
      const Real sr = x.cs(m,CLN_HP+q,k,j,i) + x.cs(m,CLN_HM+q,k,j,i);
      // (the kinetic stress w (F_kin/c - s^2 J), 0 for an isotropic field, times the mean
      // area of the family's two faces: area form; w D (J - E) stays a gradient)
      const Real ac = 0.5*(x.area[q](m,k,j,i) + x.area[q](m,k+dk,j+dj,i+di));
      const Real jc = fmin(x.cs(m,CLN_J,k,j,i), fmax(ec, 0.0));
      v(m,3*q,k,j,i) = (1.0 - wc)*a;
      v(m,3*q+1,k,j,i) = wc*dd*(jc - ec);
      v(m,3*q+2,k,j,i) = ac*wc;
      v(m,15+q,k,j,i) = -ac*wc*sr*sr*jc;
    }
  });
  CleanGhostFill(cln_v, 0, 18);
  // the cell terms of cell (k,j,i) at its face (kf,jf,iff) normal to dir, along dir: al,
  // the lagged pressure, w and the half ranges h+, h- (of the owner's family, oriented)
  auto cell = KOKKOS_LAMBDA(const int dir, const int m, const int kf, const int jf,
                            const int iff, const int k, const int j, const int i,
                            Real t[6]) {
    Real s;
    const int q = CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2;
    for (int n = 0; n < 3; ++n) t[n] = v(m,3*q+n,k,j,i);
    const Real hp = x.cs(m,CLN_HP+q,k,j,i), hm = x.cs(m,CLN_HM+q,k,j,i);
    t[3] = (s > 0.0) ? hp : -hm;
    t[4] = (s > 0.0) ? hm : -hp;
    t[5] = v(m,15+q,k,j,i);
  };
  // F of the other face of the family of cell (k,j,i) along face (kf,jf,iff), ours
  auto far = KOKKOS_LAMBDA(const int dir, const int m, const int kf, const int jf,
                           const int iff, const int k, const int j, const int i) {
    Real s;
    const int n = CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, true, s);
    return s*v(m,9+n,k,j,i);
  };
  const int nn[3] = {ie - is + 1, je - js + 1, ke - ks + 1};
  for (int it = 0; it < 64; ++it) {
    Real dmax = 0.0, fmx = 0.0;
    for (int dir = 0; dir < 3; ++dir) {
      const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
      const int nb = (dir == 2) ? nn[1] : nn[2];
      auto fc = x.fc[dir];
      Real d0 = 0.0, f0 = 0.0;
      Kokkos::parallel_reduce("cln_fcs2", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              (nmb1 + 1)*na*nb),
      KOKKOS_LAMBDA(const int idx, Real &d_max, Real &f_max) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        Real cp = 0.0, dp = 0.0;
        for (int t = 0; t <= n; ++t) {
          int k, j, i;
          face(t, k, j, i);
          const int kl = k - dk, jl = j - dj, il = i - di;
          Real lo = 0.0, hi = 0.0, dg = 1.0, r = fc(m,CLF_N,k,j,i);
          Real tl[6], tr[6];
          cell(dir, m, k, j, i, kl, jl, il, tl);
          cell(dir, m, k, j, i, k, j, i, tr);
          if (!CleanBoundaryFace(x, dir, m, k, j, i) &&
              (tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0)) {
            const Real g = 0.5*x.h*x.c/CleanCellFaceDistance(x, dir, m, k, j, i);
            const Real it_ = 1.0/fc(m,CLF_TH,k,j,i);
            const Real q = x.h*x.c/(0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,il)));
            lo = -g*tl[0] - q*tl[2]*tl[3];
            hi = -g*tr[0] + q*tr[2]*tr[4];
            dg = it_ - lo - hi + q*(tr[2]*(tr[3] + tr[4]) - tl[2]*(tl[3] + tl[4]));
            r = fc(m,CLF_N,k,j,i)*it_ + q*x.c*(x.area[dir](m,k,j,i)*(tl[1] - tr[1]) +
                                               tl[5] - tr[5]);
            // (the faces beyond a MeshBlock line end: their owners', lagged)
            if (t == 0) {r -= lo*far(dir, m, k, j, i, kl, jl, il); lo = 0.0;}
            if (t == n) {r -= hi*far(dir, m, k, j, i, k, j, i); hi = 0.0;}
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
          const int kl = k - dk, jl = j - dj, il = i - di;
          fn = x.w(m,CLW_R,k,j,i) - ((t < n) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
          const Real fo = (t > 0) ? v(m,9+2*dir+1,kl,jl,il) : v(m,9+2*dir,k,j,i);
          d_max = fmax(d_max, fabs(fn - fo));
          f_max = fmax(f_max, fmax(fabs(fn), 1.0e-8*x.c*fabs(x.w(m,e,k,j,i))));
          if (t > 0) v(m,9+2*dir+1,kl,jl,il) = fn;
          if (t < n) v(m,9+2*dir,k,j,i) = fn;
        }
      }, Kokkos::Max<Real>(d0), Kokkos::Max<Real>(f0));
      dmax = fmax(dmax, d0);
      fmx = fmax(fmx, f0);
    }
    CleanGhostFill(cln_v, 9, 6);
    if (CleanReduce(dmax, true) <= 1.0e-14*CleanReduce(fmx, true)) break;
  }
  // F of every face (CLF_N + 2), a face shared with another MeshBlock that of the lower
  // gid, and the cells' slots of it
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    par_for("cln_fcs3", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) {
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);
        return;
      }
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      Real f = lin ? v(m,9+2*dir+1,kl,jl,il) : v(m,9+2*dir,k,j,i);
      if (!(lin && rin)) {   // (the other MeshBlock's solution, in our axes)
        const int kg = lin ? k : kl, jg = lin ? j : jl, ig = lin ? i : il;
        const int ko = lin ? kl : k, jo = lin ? jl : j, io = lin ? il : i;
        Real s;
        const int n = CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s);
        if (x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io)) f = s*v(m,9+n,kg,jg,ig);
      }
      fc(m,CLF_N+2,k,j,i) = f;
    });
  }
  par_for("cln_fcs4", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int q = 0; q < 3; ++q) {
      v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
      v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
    }
  });
  CleanGhostFill(cln_v, 9, 6);
  // the energy flux's memory d, F and the dissipation term of every face
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    par_for("cln_fcs5", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real f = fc(m,CLF_N+2,k,j,i);
      if (fk && e == CLW_X && dok) {   // (and its slope, CleanDInc)
        const Real dd = fc(m,CLF_N,k,j,i) - f - (1.0 - fc(m,CLF_TH,k,j,i))*
                        fc(m,CLF_N+3,k,j,i);
        fc(m,CLF_N+3,k,j,i) += dd;
        if (kd > 0) fc(m,kd,k,j,i) = dd/x.h;
      }
      fc(m,CLF_N,k,j,i) = f;
      Real tl[6], tr[6];
      cell(dir, m, k, j, i, kl, jl, il, tl);
      cell(dir, m, k, j, i, k, j, i, tr);
      if (!(tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0)) return;
      const Real fm = far(dir, m, k, j, i, kl, jl, il);
      const Real fp = far(dir, m, k, j, i, k, j, i);
      const Real pl = tl[2]*(tl[3]*fm + tl[4]*f)/x.c + tl[5];
      const Real pr = tr[2]*(tr[3]*f + tr[4]*fp)/x.c + tr[5];
      fc(m,CLF_N+1,k,j,i) = (tr[0]*(fp - f) - tl[0]*(f - fm))/
                            (2.0*x.c*CleanCellFaceDistance(x, dir, m, k, j, i)) -
          (pr - pl + x.area[dir](m,k,j,i)*(tr[1] - tl[1]))/
          (0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,kl,jl,il)));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagCsLateral
//! \brief staggered on the cubed sphere: for E in work vector e, the lateral stress sums
//! of every cell and face family q, Y_q = 1/4 sum over its two faces of A n.(P' - P),
//! P = D E (Cartesian, n outward, P' of the cell across).  The momentum F_n of a face
//! changes by n.div P of its two half cells: A n.(P_R - P_L) and the Y of the families
//! other than the face's (CleanFaceStress), exact for a uniform P.  The owner computes Y
//! and the ghosts take it, so that a seam face sees the neighbour's true geometry.

void RadiationM1::CleanStagCsLateral(const int e) {
  if (!pmy_pack->pmesh->use_cubed_sphere || cln_cell) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x = CleanMakeArgs();
  par_for("cln_cslat", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int q = 0; q < 3; ++q) {
      const int dk = (q == 2), dj = (q == 1), di = (q == 0);
      Real y[3] = {0.0, 0.0, 0.0};
      for (int s = -1; s <= 1; s += 2) {
        const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
        if (CleanBoundaryFace(x, q, m, kf, jf, iff)) continue;
        const Real a = 0.25*s*x.area[q](m,kf,jf,iff);
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
            y[b] += a*x.g(m,CLG_N+3*q+c,kf,jf,iff)*(x.cs(m,n,k+s*dk,j+s*dj,i+s*di)*
                    x.w(m,e,k+s*dk,j+s*dj,i+s*di) - x.cs(m,n,k,j,i)*x.w(m,e,k,j,i));
          }
        }
      }
      for (int b = 0; b < 3; ++b) x.w(m,CLW_L+3*q+b,k,j,i) = y[b];
    }
  });
  CleanGhostFill(cln_w, CLW_L, 9);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagDrag
//! \brief clean_stag_vimp (staggered, gas coupling): the gas velocity implicit in the
//! face momentum, eliminated face by face so that the E solve stays scalar.  A change
//! dm_f of the momentum face f gives up moves its two half cells (mass rho_c/2, share
//! wt_c = chi_c/(chi_L + chi_R)) by dv_c = 2 wt_c dm_f/rho_c, and the advective flux by
//! (1 - theta) dA_f = beta c^2 dm_f with
//!   beta = (1 - theta) sum_c wt_c (1 + D_nn,c) E_c/(rho_c c^2),
//! E of the iterate.  With dm_f = -[(F - F^k)/c^2 + h d_n(D (E - E^k))] relative to the
//! flux of the iterate's velocity,
//!   F = (1 - b) F^H(E) + b [Q - h c^2 d_n(D E)],  b = beta/(1 + beta),
//!   Q = F^k + h c^2 d_n(D E^k) = G - c^2 dm_f^k (CLF_DB, CLF_DQ; CleanFaceWeights,
//!   CleanFaceCross, CleanFaceSource)
//! and the energy flux takes the same increment.  When the loop has converged F = F^H:
//! the solution is that of the lagged velocity.  first: the first iterate, whose
//! velocity is that of the old vector: dm_f^k = 0, Q = G; with force_reference = clean
//! (or wb_arad) the gas holds the reference force, which balances the radiation force of
//! the start, so Q is that of E^k (no drag kick); b = 0 while the gas momentum is frozen
//! and before the reference force is set.  vel (CleanVelocityUpdate with the F solve's
//! F): first the face's term of CleanMomentumChange less its share wt, with the stress
//! once per face for both (CLF_PT); without clean_stag_vimp only that.

void RadiationM1::CleanStagDrag(const bool first, const bool vel) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanStagCsLateral(CLW_EK);   // (cubed sphere: the lateral stress of the iterate)
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  xh.drag = x.drag && !first;   // (the flux the iterate's velocity came from)
  auto uh = FluidRef::Get(pmy_pack).u0;
  const bool off = cln_vzero || (first && force_ref == M1_FREF_CLEAN && !cln_fref_set);
  const bool lag = !first || (force_ref != M1_FREF_NONE);
  const bool dr = cln_vimp;
  const int nf = x.nf;   // (one launch for the faces of all directions)
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  par_for("cln_drag", DevExeSpace(), 0, nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    for (int dir = 0; dir < nf; ++dir) {
      const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      if (i > ie+di || j > je+dj || k > ke+dk) continue;
      auto &fc = x.fc[dir];
      if (CleanBoundaryFace(x, dir, m, k, j, i)) continue;   // (no F equation)
      const Real s = (vel || lag) ? CleanFaceStress(x, CLW_EK, dir, m, k, j, i, true)
                                  : 0.0;
      if (vel) {
        fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) + x.h*s -
                             x.h*fc(m,CLF_N+1,k,j,i);
      }
      if (!dr) continue;
      const Real q = !lag ? fc(m,CLF_G,k,j,i) :
                     CleanFaceFlux(xh, CLW_EK, dir, m, k, j, i, true) + x.h*x.c*x.c*s;
      const Real xl = x.op(m,M1_OP_T,k-dk,j-dj,i-di), xs = xl + x.op(m,M1_OP_T,k,j,i);
      const Real wl = (xs > 0.0) ? xl/xs : 0.5;
      auto ce = [&](const int kk, const int jj, const int ii) {
        return (1.0 + x.cs(m,CLN_D11+dir,kk,jj,ii))*fmax(x.w(m,CLW_EK,kk,jj,ii), 0.0)/
               uh(m,IDN,kk,jj,ii);
      };
      const Real bt = off ? 0.0 : (1.0 - fc(m,CLF_TH,k,j,i))*(wl*ce(k-dk, j-dj, i-di) +
                      (1.0 - wl)*ce(k, j, i))/(x.c*x.c);
      fc(m,CLF_DB,k,j,i) = bt/(1.0 + bt);
      fc(m,CLF_DQ,k,j,i) = q;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanWriteBack
//! \brief stores the solution: face fluxes, E, cell fluxes (face means), the gas momentum
//! and energy, and the stage slopes K = (U - U_old)/h.  The gas takes exactly the momentum
//! and energy the radiation gives up, so both totals change only by face fluxes; with
//! force_reference less the reference force, and with wb_arad (force_reference_work =
//! split) E pays the reference work h d arad_ref.v at the velocity after the kick, which
//! the problem generator gives the gas (rad_m1_implicit.cpp:12438-12465); the explicit
//! source esrc is the radiation's (the gas takes the exchange only).  On the
//! cubed sphere the gas takes the vector of CleanStagMomentum in the hydro's components
//! (CleanGasChange), and the cell flux of the outputs is the face momenta's vector
//! (CleanStagCsVector) in the frame of y: total energy and L_z (of that vector) are
//! conserved to round-off.

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
  CleanStagCsLateral(CLW_X);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;   // (clean_stag_split: the faces store the momentum state)
  auto kk = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;
  // h div Phi of the solution (CLW_Q, a slot the F solve leaves alone) before the F solve
  // below moves the energy flux's memory d (clean_stag_fkin): E takes the residual of the
  // operator it was solved with; with clean_stag_split also the cell flux, that of Phi
  // (the faces keep F^H)
  auto u0_ = u0;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  par_for("cln_wbdiv", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    x.w(m,CLW_Q,k,j,i) = CleanFluxDivergence(x, CLW_X, m, k, j, i, true);
    for (int d = 0; d < 3 && xh.hi; ++d) {
      u0_(m,M1_F1+d,k,j,i) = ((d == 1 && !multi_d) || (d == 2 && !three_d)) ? 0.0 :
          0.5*(CleanFaceFlux(x, CLW_X, d, m, k, j, i, true) +
               CleanFaceFlux(x, CLW_X, d, m, k+(d == 2), j+(d == 1), i+(d == 0), true));
    }
  });
  // (time2_face_fct: the energy flux Phi of every face, which moved E, with its slope)
  if (t2_pfst && slopes && t2pf1.extent_int(0) > nmb1) {
    auto pf = (t2_solve == M1_T2S_STAGE1) ? t2pf2 : t2pf1;
    for (int d = 0; d < x.nf; ++d) {
      const int dk = (d == 2), dj = (d == 1), di = (d == 0);
      par_for("cln_wbpf", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        pf(m,d,k,j,i) = CleanFaceFlux(x, CLW_X, d, m, k, j, i, true);
      });
    }
    if (t2_solve == M1_T2S_STAGE1) {t2pf2_ok = true;} else {t2pf1_ok = true;}
  }
  // (the history's boundary energy: + h A Phi of the physical x1 faces, CleanBndAcc)
  if (cln_lacc.extent_int(0) > nmb1) {
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;
    const bool stg = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
    par_for("cln_wbbnd", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      la(m,s,k,j) = !CleanBoundaryFace(x, 0, m, k, j, i) ? 0.0 :
                    la(m,stg ? 2+s : s,k,j) + h*ar(m,k,j,i)*CleanFaceFlux(x, CLW_X, 0, m,
                                                                     k, j, i, true);
    });
  }
  CleanStagFSolve(CLW_X);   // (the face F of the F solve and its dissipation)
  xh.fs = true;

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
      const Real f = fc(m,CLF_N,k,j,i);
      if (slopes) kk(m,tf,k,j,i) = (f - fc(m,CLF_G,k,j,i))/h;
      f0(m,k,j,i) = f;
    });
  }

  // E, cell fluxes, gas
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  const int fref = force_ref;
  const bool wsp = (fref == M1_FREF_WB_ARAD || fref == M1_FREF_WB_MD);
  auto fr = cln_fref;
  auto aref = arad_ref;
  const bool wpa = wsp && FrefWaccOn();   // (fref_wacc: CleanRightHandSide)
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : h/pmy_pack->pmesh->dt;
  const bool wex = wsp && fref_wacc_t2 && !fref_pay_post && !cln_cell &&
                   time_scheme == M1_TIME_HESDIRK2;   // (no slope carries a payment)
  const bool fm = icr_fmeas && cln_vzero;   // (CleanSteadyStart: icr_dm)
  auto fmd = icr_dm;
  const bool eso = esrc_on;
  auto es = esrc;
  const Real arad_ = arad;
  // (wb_arad_md, clean_init = steady: the reference of the steady state, cln_fref 0..2
  // its acceleration in the components of dm, 3 + 2 dir + side the potential P that the
  // hydro walks of <hydro>/wb_dirs = all take along dir: the line integral over each half
  // cell (side 0 face to centre, 1 centre to face) of its acceleration, 2 cf/(h rho), the
  // half cell's share of the face family's momentum (the cubed sphere's lateral halves:
  // CleanSteadyStart).  None outside the Roche lobe, problem/wb_lobe_cut)
  const bool mdc = cln_mdcap;
  const Real phx = fl.ph ? fl.ph->wb_phimax : (fl.pm ? fl.pm->wb_phimax : 1.0e300);
  auto phc = fl.phicc_wb;
  auto cx1 = pmy_pack->pcoord->x1v, cx2 = pmy_pack->pcoord->x2v;
  auto cx3 = pmy_pack->pcoord->x3v, cf1 = pmy_pack->pcoord->xx1f;
  auto cf2 = pmy_pack->pcoord->xx2f, cf3 = pmy_pack->pcoord->xx3f;
  const Real tgf = cln_gtf ? eos.tfloor : 0.0;   // (clean_gas_tfloor)
  const bool gtf = cln_gtf, pgas = impl_pos_gas;
  const Real pgf = impl_pos_gas_frac, efl = e_floor;
  auto pc_ = pos_cnt_d;
  par_for("cln_wbcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = x.w(m,CLW_X,k,j,i);
    if (mdc) {
      Real dr[3], mg[3], cf[2], wf[2];
      CleanStagMomentum(xh, CLW_X, m, k, j, i, dr, mg);
      const Real d = uh(m,IDN,k,j,i), s = (phc(m,k,j,i) > phx) ? 0.0 : 1.0/(h*d);
      for (int n = 0; n < 3; ++n) {
        fr(m,n,k,j,i) = s*dr[n];
        fr(m,3+2*n,k,j,i) = fr(m,4+2*n,k,j,i) = 0.0;
        if (n >= x.nf) continue;
        CleanMomentumChange(xh, CLW_X, n, m, k, j, i, nullptr, cf, wf);
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        for (int sd = 0; sd < 2; ++sd) {
          const int kf = k + sd*dk, jf = j + sd*dj, iff = i + sd*di;
          const Real l = !(x.sph || x.cub) ? 0.5/x.idx[n] : ((n == 0) ?
              fabs(cf1(m,iff) - cx1(m,i)) : ((n == 1) ? cx1(m,i)*fabs(cf2(m,jf) -
              cx2(m,j)) : cx1(m,i)*sin(cx2(m,j))*fabs(cf3(m,kf) - cx3(m,k))));
          fr(m,3+2*n+sd,k,j,i) = 2.0*s*cf[sd]*l;
        }
      }
    }
    Real dm[3] = {0.0, 0.0, 0.0}, dh[3] = {0.0, 0.0, 0.0};
    Real den = 0.0, wsl = 0.0;
    if (feedback) {
      Real mg[3], pg[3], rf[3];
      CleanStagMomentum(xh, CLW_X, m, k, j, i, dm, mg);
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      CleanStagRefForce(x, fref, fr, aref, uh(m,IDN,k,j,i), m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {   // (force_reference, frozen gas momentum)
        if (fref != M1_FREF_NONE) dm[n] -= h*rf[n];
        if (fm && n == 0) fmd(m,k,j,i) = dm[0];
        if (x.vzero) dm[n] = 0.0;
      }
      // (the kinetic energy the gas gains, the work W of the E row at convergence; with
      // wb_arad wr, the reference work at the velocity after the kick, which E pays and
      // the problem generator gives the gas)
      Real wk = 0.0, wr = 0.0;
      for (int n = 0; n < 3; ++n) {
        wk += (pg[n] + 0.5*dm[n])*dm[n]/uh(m,IDN,k,j,i);
        if (wsp && !x.vzero) wr += h*rf[n]*(pg[n] + dm[n])/uh(m,IDN,k,j,i);
        dh[n] = dm[n];
      }
      if (x.cub) {   // (the hydro's components of dm, from the frame of y)
        Real q[3][3], dy[3];
        CleanCsFrame(x, m, k, j, i, q);
        for (int a = 0; a < 3; ++a) dy[a] = q[a][0]*dm[0] + q[a][1]*dm[1] + q[a][2]*dm[2];
        CleanGasChange(x, m, k, j, i, dy, dh);
      }
      for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += dh[n];
      // gas energy: exactly what the radiation gives up besides its flux divergence, so
      // that total energy changes only by the face fluxes
      Real div = x.idx[0]*(f1(m,k,j,i+1) - f1(m,k,j,i));
      if (multi_d) div += x.idx[1]*(f2(m,k,j+1,i) - f2(m,k,j,i));
      if (three_d) div += x.idx[2]*(f3(m,k+1,j,i) - f3(m,k,j,i));
      // (clean_stag_split: the faces store F^H, E moved by the energy flux)
      // (and esrc, the radiation's external source: the gas takes the exchange only)
      den = -(e - x.w(m,CLW_EO,k,j,i)) - ((xh.hi || x.sph || x.cub) ?
            x.w(m,CLW_Q,k,j,i) : h*div) + (eso ? h*es(m,k,j,i) : 0.0);
      if (x.vzero) {   // (steady start: the eint of the iterate's T, not the residual of
        Real ei, p, chi_rho, chi_t, cv;   //  a step h -> infinity, as the cell layout)
        eos.ThermoAt(uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), ei, p, chi_rho, chi_t, cv);
        den = ei - x.w(m,CLW_EI,k,j,i);
      } else {
        // the gas takes the exchange of the solution with the iterate's T, the beta terms
        // and W, E the rest (its fluxes): the linear solve's residual goes to E, not eint
        const Real t = x.w(m,CLW_T,k,j,i);
        const Real ex = h*x.c*(x.op(m,M1_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_*
                               SQR(SQR(t)) -
                               x.op(m,M1_OP_E,k,j,i)*fmax(e, 0.0));
        // (the residual over the row's local diagonal: all of it in a scattering cell,
        // little where the exchange is stiff, whose X the solve fixes to tol/diagonal)
        // (wpa: E pays fref_wacc; the slope keeps no payment (wex: FrefWaccStage) or the
        // payment wr, which the next step's stage increments re-apply)
        const Real wp = wpa ? fra*wc(m,k,j,i) : wr;
        // (clean_gas_tfloor: the gas eint changes by that of its T root (>= tfloor),
        // e(rho, T) - eint, not by the exchange, a difference of two terms ~E >> eint
        // where r ~ 1: positive by construction, equal to it at a converged root; E
        // takes all the rest, the solve's residual and any work the gas's kinetic
        // energy change did not receive)
        Real eg = x.w(m,CLW_EP,k,j,i) - ex;
        if (gtf) {
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(uh(m,IDN,k,j,i), fmax(t, tgf), ei, p, chi_rho, chi_t, cv);
          eg = ei - x.w(m,CLW_EI,k,j,i);
        }
        const Real rd = (den - eg - wk - wp)/(gtf ? 1.0 :
                        1.0 + x.h*x.c*x.op(m,M1_OP_E,k,j,i)*x.w(m,CLW_RY,k,j,i));
        e += rd;
        den -= rd + wp;
        wsl = wex ? wp : wp - wr;
        // (implicit_pos_gas, as the old production's write-back: a gas eint below e_min =
        // max(e(rho, tfloor), implicit_pos_gas_frac x its old eint) is raised to it with
        // energy from the cell's E above e_floor, gas + E of the cell kept; counted)
        if (pgas) {
          Real ef = 0.0, pr, chi_rho, chi_t, cv;
          if (eos.tfloor > 0.0) {
            eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
          }
          const Real ei = x.w(m,CLW_EI,k,j,i), emin = fmax(ef, pgf*fmax(ei, 0.0));
          const Real take = fmin(emin - ei - (den - wk), e - efl);
          if (take > 1.0e-12*emin) {   // (beyond eint's round-off at e_min)
            const Real vol = (x.sph || x.cub) ? x.vol(m,k,j,i) :
                             1.0/(x.idx[0]*x.idx[1]*x.idx[2]);
            den += take;
            e -= take;
            Kokkos::atomic_add(&pc_(M1_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(M1_POS_GAS_DE), take*vol);
          }
        }
      }
      uh(m,IEN,k,j,i) += den;
      x.w(m,CLW_JV,k,j,i) = den - wk;   // (its change of the gas eint, cln_neint)
    }
    // (an E below zero within the solver tolerance, as the old production's
    // implicit_pos_floor: raised to zero with the energy of the cell's gas above e_min,
    // counted; what the gas cannot cover stays, and the solve is not accepted
    // (cln_neint); without the gas the raise is counted as created)
    if (e < 0.0) {
      const Real vol = (x.sph || x.cub) ? x.vol(m,k,j,i) :
                       1.0/(x.idx[0]*x.idx[1]*x.idx[2]);
      const bool gas = feedback && !x.vzero;
      Real take = gas ? 0.0 : -e;
      if (gas) {
        Real ef = 0.0, pr, chi_rho, chi_t, cv;
        if (eos.tfloor > 0.0) {
          eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
        }
        const Real ei = x.w(m,CLW_EI,k,j,i), emin = fmax(ef, pgf*fmax(ei, 0.0));
        take = fmin(-e, fmax(ei + x.w(m,CLW_JV,k,j,i) - emin, 0.0));
        den -= take;
        uh(m,IEN,k,j,i) -= take;
        x.w(m,CLW_JV,k,j,i) -= take;
        Kokkos::atomic_add(&pc_(M1_POS_FLR_DE), take*vol);
      } else {
        Kokkos::atomic_add(&pc_(M1_POS_FLR_UN), take*vol);
      }
      e += take;
      Kokkos::atomic_add(&pc_(M1_POS_FLR), 1.0);
    }
    if (slopes) {
      kk(m,M1_T2_E,k,j,i) = (e + wsl - x.w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,M1_T2_M1+n,k,j,i) = dh[n]/h;
      kk(m,M1_T2_EN,k,j,i) = den/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    if (!xh.hi) {   // (clean_stag_split: the cell flux of Phi, above)
      u0_(m,M1_F1,k,j,i) = 0.5*(f1(m,k,j,i) + f1(m,k,j,i+1));
      u0_(m,M1_F1+1,k,j,i) = multi_d ? 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i)) : 0.0;
      u0_(m,M1_F1+2,k,j,i) = three_d ? 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i)) : 0.0;
    }
    if (x.cub) {   // (cubed sphere: the face momenta's vector, CleanStagCsVector, in the
      Real c[6], w[6], fv[3], q[3][3];   //  frame of y: its L_z is CleanStagMomentum's)
      for (int d = 0; d < 6; ++d) w[d] = 1.0;
      c[0] = f1(m,k,j,i);
      c[1] = f1(m,k,j,i+1);
      c[2] = f2(m,k,j,i);
      c[3] = f2(m,k,j+1,i);
      c[4] = f3(m,k,j,i);
      c[5] = f3(m,k+1,j,i);
      CleanStagCsVector(x, m, k, j, i, c, w, fv);
      CleanCsFrame(x, m, k, j, i, q);
      for (int d = 0; d < 3; ++d) {
        u0_(m,M1_F1+d,k,j,i) = q[d][0]*fv[0] + q[d][1]*fv[1] + q[d][2]*fv[2];
      }
    }
  });
  // cells the write-back left with gas eint < 0, or below the EOS floor, by lowering it
  // (the solve is then not accepted, as in the cell layout; an eint the hydro left there
  // is not the solve's)
  cln_neint = 0;
  if (!feedback) return;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nx1 = id.nx1, nx2 = id.nx2, nx3 = id.nx3, nji = nx2*nx1, nkji = nx3*nji;
  Real nb = 0.0;
  Kokkos::parallel_reduce("cln_neint", CleanRp256(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji + ks;
    const int j = (idx%nji)/nx1 + js, i = idx%nx1 + is;
    Real p[3];
    CleanGasMomentum(x, uh, m, k, j, i, p);
    const Real ke = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/uh(m,IDN,k,j,i);
    // (and below the eint of the EOS temperature floor: the hydro's floor would add
    // energy the radiation's budget does not see)
    Real ef = 0.0;
    if (eos.tfloor > 0.0) {
      Real pr, chi_rho, chi_t, cv;
      eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
    }
    const Real tz = 1.0e-12*(fabs(uh(m,IEN,k,j,i)) + ke + ef);   // (eint's round-off)
    if ((!(CleanGasEint(x, uh, ke, m, k, j, i) > ef - tz) &&
         !(x.w(m,CLW_JV,k,j,i) >= 0.0) && x.w(m,CLW_EI,k,j,i) > ef + tz) ||
        u0_(m,M1_E,k,j,i) < 0.0) n += 1.0;
  }, Kokkos::Sum<Real>(nb));
  cln_neint = static_cast<int>(CleanReduce(nb, false));
}

}  // namespace radm1
