//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_clean.cpp
//! \brief Stage solve of the clean implicit radiation scheme: one implicit stage (or a
//! backward-Euler step) for E, coupled to the gas by Newton iteration on the exchange and
//! a fixed-point iteration on the velocity.  Discrete equations: rad_vet_clean.hpp.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_clean.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "utils/repro_sum.hpp"
#include "utils/scratch_view.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
// VET_TU_PART: compiled as 3 TUs (this file = part 0, rad_vet_clean_p<k>.cpp = part k).
#ifndef VET_TU_PART_N
#define VET_TU_PART_N 0
#endif
#define VET_TU_PART(k) (VET_TU_PART_N == (k))

namespace radvet {
extern bool g_vra_adv;   // vet_rotate_advance (rad_vet_vetgd.cpp)
// heavy per-cell reductions at 256 threads per block (HIP's 1024 default caps VGPRs at
// 128 on gfx942 and spills); max reductions and counts are bitwise.
using CleanRp256 = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 1>>;
// heavy per-cell passes: 2 blocks of 256 per SM (<= 128 registers) on CUDA, 1 on HIP
#if defined(KOKKOS_ENABLE_CUDA)
static constexpr int CleanLbMb(int mb) {return mb;}
#else
static constexpr int CleanLbMb(int) {return 1;}
#endif
using CleanRpLb2 =
    Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, CleanLbMb(2)>>;
// 4-D par_for (same index flattening) at MB blocks of 256 per SM
template <int MB, typename F>
static void CleanParForLb(const char *name, const int nmb1, const int ks, const int ke,
                          const int js, const int je, const int is, const int ie,
                          const F &f) {
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int mb = 2;   // HIP: cln_rhs at 2 waves per SIMD (faster despite spills)
#else
  constexpr int mb = CleanLbMb(MB);
#endif
  Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,
                       mb>>(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/ni;
    f(m, ks + k, js + j, is + idx - m*nkji - k*nji - j*ni);
  });
}
//----------------------------------------------------------------------------------------
//! \fn void CleanFatal
//! \brief stops the run with a message about an unsupported clean-scheme setup

static void CleanFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  // _Exit: no static destructors behind live Kokkos views
  std::_Exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanReduceMax
//! \brief the maxima of a, b, c over all ranks, one reduction (CleanReduce)

static void CleanReduceMax(Real &a, Real &b, Real &c) {
  Real v[3] = {a, b, c};
  CleanMaxRanks(v, 3);
  a = v[0];
  b = v[1];
  c = v[2];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagGasMomentum
//! \brief staggered: the gas momentum of cell (k,j,i) in the mesh's face-flux components
//! (Cartesian on the cubed sphere)

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

//! \brief staggered: v from CleanStagGasMomentum's components to the hydro's (changed on
//! the cubed sphere only)

KOKKOS_INLINE_FUNCTION
static void CleanStagHydro(const CleanArgs &x, const int m, const int k, const int j,
                           const int i, Real v[3]) {
  if (!x.cub) return;
  Real q[3][3];
  const Real y[3] = {v[0], v[1], v[2]};
  CleanCsFrame(x, m, k, j, i, q);
  for (int d = 0; d < 3; ++d) v[d] = q[d][0]*y[0] + q[d][1]*y[1] + q[d][2]*y[2];
}

//! \brief staggered: the radial component of a vector v in the components of
//! CleanStagGasMomentum (spherical polar: v[0]; the cubed sphere: r_hat . v)

KOKKOS_INLINE_FUNCTION
static Real CleanStagRadial(const CleanArgs &x, const int m, const int k, const int j,
                            const int i, const Real v[3]) {
  if (!x.cub) return v[0];
  return x.g(m,CLG_R,k,j,i)*v[0] + x.g(m,CLG_R+1,k,j,i)*v[1] + x.g(m,CLG_R+2,k,j,i)*v[2];
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagRefForce
//! \brief staggered: the reference force density withheld from the gas of cell (k,j,i),
//! in CleanStagGasMomentum's components.  clean: the stored force; wb_arad: d arad_ref
//! along x1 (radial on the cubed sphere: no torque); wb_arad_md: d times the stored
//! acceleration; else 0.  The gas then takes the residual of the radiation force.

KOKKOS_INLINE_FUNCTION
static void CleanStagRefForce(const CleanArgs &x, const int fref,
                              const DvceArray5D<Real> &fr, const DvceArray4D<Real> &aref,
                              const Real d, const int m, const int k, const int j,
                              const int i, Real rf[3]) {
  Real q[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  if (fref == VET_FREF_WB_ARAD && x.cub) CleanCsFrame(x, m, k, j, i, q);
  for (int b = 0; b < 3; ++b) {
    rf[b] = (fref == VET_FREF_CLEAN) ? fr(m,b,k,j,i) :
            ((fref == VET_FREF_WB_ARAD) ? d*aref(m,k,j,i)*q[0][b] :
            ((fref == VET_FREF_WB_MD) ? d*fr(m,b,k,j,i) : 0.0));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagMomentum
//! \brief staggered: the radiation's momentum change of cell (k,j,i) over the stage, the
//! cell's shares of its faces' absorption (CleanMomentumChange; mag: the terms' size), in
//! CleanStagGasMomentum's components.
//! Cubed sphere: the face changes combine by CleanStagCsVector (weights 1/2).  Lateral
//! faces are no phi faces, so L_z does not telescope: dlz = -V R (dF_phi/c^2 + dp_phi) is
//! booked as the radiation's (CleanLzBook).  3-D spherical polar: dlz of the phi faces,
//! which telescopes to the L_z crossing the boundaries.

template <bool kFs = false>
KOKKOS_INLINE_FUNCTION
static void CleanStagMomentum(const CleanArgs &x, const int e, const int m, const int k,
                              const int j, const int i, Real dm[3], Real mag[3],
                              Real *dlz = nullptr) {
  Real cf[6], wf[6];
  for (int d = 0; d < 3; ++d) {
    mag[d] = 0.0;
    dm[d] = (d < x.nf) ? CleanMomentumChange<kFs>(x, e, d, m, k, j, i, &mag[d],
                                                  &cf[2*d], &wf[2*d]) : 0.0;
  }
  if (!x.cub) {
    if (dlz && x.sph && x.nf > 2) {
      Real df[2];
      for (int s = 0; s < 2; ++s) {
        const Real f = (kFs || x.fs) ? x.fc[2](m,CLF_N,k+s,j,i)
                                     : CleanFaceFlux(x, e, 2, m, k+s, j, i, true);
        df[s] = f - x.fc[2](m,CLF_G,k+s,j,i);
      }
      *dlz = (-0.5*(df[0] + df[1])/(x.c*x.c) - dm[2])*x.vol(m,k,j,i)*x.x1v(m,i)*
             x.tg(m,2,j);
    }
    return;
  }
  Real dv[3], q[3][3];
  const Real hw[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  CleanStagCsVector(x, m, k, j, i, cf, hw, dv);
  const Real rc = CleanCsFrame(x, m, k, j, i, q);
  const Real dr = q[0][0]*dv[0] + q[0][1]*dv[1] + q[0][2]*dv[2];
  const Real dt = q[1][0]*dv[0] + q[1][1]*dv[1] + q[1][2]*dv[2];
  const Real dp = q[2][0]*dv[0] + q[2][1]*dv[1] + q[2][2]*dv[2];
  if (dlz) {
    const Real one[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Real df[6], fd[3];
    for (int d = 0; d < 6; ++d) {
      const int n = d/2, kf = k + (d%2)*(n == 2), jf = j + (d%2)*(n == 1);
      const int iff = i + (d%2)*(n == 0);
      const Real f = (kFs || x.fs) ? x.fc[n](m,CLF_N,kf,jf,iff)
                                   : CleanFaceFlux(x, e, n, m, kf, jf, iff, true);
      df[d] = f - x.fc[n](m,CLF_G,kf,jf,iff);
    }
    CleanStagCsVector(x, m, k, j, i, df, one, fd);
    *dlz = (-(q[2][0]*fd[0] + q[2][1]*fd[1] + q[2][2]*fd[2])/(x.c*x.c) - dp)*
           x.vol(m,k,j,i)*rc;
  }
  for (int b = 0; b < 3; ++b) dm[b] = dr*q[0][b] + dt*q[1][b] + dp*q[2][b];
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanCsTorque
//! \brief cubed sphere, clean_cs_lz_tel: the gas L_z change of cell (k,j,i) over the
//! stage in telescoping form; rl: the radiation's (halves of V F_f lambda_f/c^2 of the
//! lateral faces, lambda_f = (M_f x n_f)_z/A_f); xs: the stress part (exact less faces').
//! The gas takes the faces' absorption torque V cf lambda plus the exact stress torque of
//! P = D E: fluxes r_f A (n x P n)_z (radial) and (M_f x P n)_z (lateral), P the two-cell
//! mean (one-sided at a radial boundary).  Interior faces cancel; CleanLzBook books the
//! boundary torque.  bf (clean_cs_lz_bnd_free): radial boundary faces take no stress
//! torque (lever arm along n), so gas and radiation L_z are exact with no ledger.

KOKKOS_INLINE_FUNCTION
static Real CleanCsTorque(const CleanArgs &x, const int e, const int m, const int k,
                          const int j, const int i, Real &rl, Real &xs,
                          const bool bf) {
  const Real v = x.vol(m,k,j,i), ic2 = 1.0/(x.c*x.c);
  Real gl = 0.0, gs = 0.0, tc = 0.0, p0[6];
  rl = 0.0;
  for (int q = 1; q < 3; ++q) {   // lateral faces: absorption and momenta
    Real cf[2], wf[2], tq[2];
    CleanMomentumChange(x, e, q, m, k, j, i, nullptr, cf, wf, tq);
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1);
      const Real f = x.fs ? x.fc[q](m,CLF_N,kf,jf,i)
                          : CleanFaceFlux(x, e, q, m, kf, jf, i, true);
      const Real df = (f - x.fc[q](m,CLF_G,kf,jf,i))*ic2;
      const Real lam = CleanCsArm(x, q, m, kf, jf, i);
      gl += v*lam*cf[s];
      gs += 0.5*v*lam*(tq[s] - df);
      rl += 0.5*v*lam*df;
    }
  }
  for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
  for (int q = 0; q < 3; ++q) {   // stress torque through the six faces
    const int dk = (q == 2), dj = (q == 1), di = (q == 0);
    for (int s = 0; s < 2; ++s) {
      const int o = 2*s - 1, kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
      const int ko = k + o*dk, jo = j + o*dj, io = i + o*di;
      const bool bnd = CleanBoundaryFace(x, q, m, kf, jf, iff);
      if (bf && bnd && q == 0) continue;
      Real pf[6], a[3], t[2];
      for (int n = 0; n < 6; ++n) {
        pf[n] = bnd ? p0[n] : 0.5*(p0[n] + x.cs(m,CLN_D11+n,ko,jo,io)*x.w(m,e,ko,jo,io));
      }
      for (int c = 0; c < 3; ++c) {
        a[c] = (q == 0) ? x.x1f(m,iff)*x.area[0](m,kf,jf,iff)*x.g.At<CLG_N>(m,c,kf,jf,iff)
                        : x.g.At<CLG_M>(m,3*(q-1)+c,kf,jf,iff);
      }
      for (int b = 0; b < 2; ++b) {   // P n, outward: x and y components
        t[b] = 0.0;
        for (int c = 0; c < 3; ++c) {
          t[b] += o*pf[(b == c) ? b : 2 + b + c]*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff);
        }
      }
      tc += a[0]*t[1] - a[1]*t[0];
    }
  }
  xs = gs - x.h*tc;
  return gl + xs;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanVelocityCell
//! \brief CleanVelocityUpdate's cell (k,j,i): the momentum change (CLW_DM), the advective
//! coefficients of the new velocity if fa (CLW_A1), and the maxima of |dv|, |v| and the
//! force scale; kFs: F known from the F solve (xh.fs, xh.fa)

template <bool kFs>
KOKKOS_INLINE_FUNCTION
static void CleanVelocityCell(const CleanArgs &x, const CleanArgs &xh,
                              const CleanArgs &xa,
                              const DvceArray5D<Real> &uh,
                              const int fref, const DvceArray5D<Real> &fr,
                              const DvceArray4D<Real> &aref, const bool cap,
                              const DvceArray5D<Real> &op, const Real arad_,
                              const bool fa,
                              const int m, const int k, const int j, const int i,
                              Real &dv_max, Real &v_max, Real &fv_max) {
  const Real d = uh(m,IDN,k,j,i);
  Real v2 = 0.0, dmv[3], mg[3], pg[3], rf[3];
  CleanStagMomentum<kFs>(xh, CLW_EK, m, k, j, i, dmv, mg);
  CleanStagGasMomentum(x, uh, m, k, j, i, pg);
  CleanStagRefForce(x, fref, fr, aref, d, m, k, j, i, rf);
  for (int dir = 0; dir < 3; ++dir) {
    Real dm = 0.0;
    if (dir == 0 || (dir == 1 && x.multi_d) || (dir == 2 && x.three_d)) {
      const Real mag = mg[dir];
      dm = dmv[dir];
      v_max = fmax(v_max, mag/d);   // |v| floor: size of the stage's force terms
      // force_reference = clean: the first coupled solve's force is the reference,
      // withheld from the gas (as wb_arad's)
      if (cap) fr(m,dir,k,j,i) = dm/x.h;
      if (fref != VET_FREF_NONE) dm -= cap ? dm : x.h*rf[dir];
      if (x.vzero) dm = 0.0;
      // force scale relative to E in the stress
      const Real e = fmax(fabs(x.w(m,CLW_EK,k,j,i)), 1.0e-300), hc = x.h*x.c;
      fv_max = fmax(fv_max, mag/d*(2.0 + (e*hc*op(m,VET_OP_E,k,j,i) + fabs(x.w(m,
               CLW_EO,k,j,i)) + hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,
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
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanPreflight
//! \brief checks the <rad_vet> / <time> keys an input names under scheme = clean: a
//! replaced key names its clean counterpart, unread keys are listed once, and a key that
//! would change the physics stops the run unless off.

#if VET_TU_PART(0)
void RadiationVet::CleanPreflight(ParameterInput *pin) {
  // block/key (* ends a prefix): nullptr read, "" no effect, "!" fatal unless off, "x"
  // removed (fatal when named), "=v" fatal unless v (old restart value), else clean key
  static const char *tab[][2] = {
    {"rad_vet/implicit_bc_advect", ""}, {"rad_vet/implicit_bc_*", nullptr},
    // inner-flux time hooks: removed
    {"rad_vet/implicit_flux_x1min_mod_*", "x"},
    {"rad_vet/implicit_flux_x1min_ramp_*", "x"},
    {"rad_vet/implicit_flux_x1min_table*", "x"},
    {"rad_vet/implicit_flux_x1m*", nullptr},
    {"rad_vet/implicit_ebath_*", nullptr}, {"rad_vet/implicit_beam*", nullptr},
    {"rad_vet/implicit_cfl", nullptr},
    {"rad_vet/implicit_mr_*", nullptr}, {"rad_vet/implicit_det_trace", nullptr},
    {"rad_vet/implicit_tol", "clean_tol_nl"}, {"rad_vet/implicit_lin_tol", "clean_tol"},
    // read by ImplicitInit; a value other than the clean key's prints a NOTE
    {"rad_vet/implicit_maxit", nullptr},
    {"rad_vet/implicit_tol_rel_max", nullptr},
    {"rad_vet/implicit_lin_maxit", nullptr},
    {"rad_vet/implicit_recon", ""},
    {"rad_vet/implicit_closure_lag", "clean_closure"},
    {"rad_vet/time2_vet", "clean_closure"}, {"rad_vet/time2_vet_col", "clean_closure"},
    {"rad_vet/time2_vet_sc", "clean_closure"},
    {"rad_vet/time2_vet_extrap", "clean_closure"},
    {"rad_vet/time2_lin_tol*", "clean_tol"},
    {"rad_vet/implicit_ghost_x1min", "!"}, {"rad_vet/ic_rad*", "clean_init"},
    {"rad_vet/implicit_predictor", nullptr}, {"rad_vet/implicit_predictor_tcap", nullptr},
    {"rad_vet/implicit_resid_fatal", nullptr},
    // ImplicitInit members the clean steps use: enthalpy flux, conservative gas source,
    // Marshak face, predictor order, hold
    {"rad_vet/implicit_enthalpy", nullptr}, {"rad_vet/implicit_gas_src", nullptr},
    {"rad_vet/implicit_marshak_face", nullptr},
    {"rad_vet/implicit_predictor_order", nullptr},
    {"rad_vet/implicit_hold", nullptr},
    {"rad_vet/implicit_pos_gas*", nullptr}, {"rad_vet/implicit_pos_floor", nullptr},
    {"rad_vet/implicit_pos_floor_cons", nullptr},
    {"rad_vet/implicit_pos_floor_gas_keep", nullptr},
    {"rad_vet/implicit_pos_gas_cons", nullptr},
    {"rad_vet/implicit_*", ""},
    // hesdirk2 stage formation (Time2FormStage / Time2Stash) and dt_retry: read;
    // ImplicitSolve's acceptance, posfix and KE cap: none
    {"rad_vet/time2_k1_reset", "!"},   // K1 reset of held cells: removed
    {"rad_vet/time2_k1_*", nullptr}, {"rad_vet/time2_face_fct*", nullptr},
    {"time/dt_retry_t2_be", nullptr}, {"time/dt_retry_max", ""},
    {"rad_vet/time2_accept_fallback", ""}, {"rad_vet/time2_posfix*", ""},
    {"rad_vet/time2_inc_kecap", ""}, {"rad_vet/time2_pred_cap", ""},
    // implicit solve's position limiters, Zalesak round cap, vet_col lateral posfix and
    // x2 vet fields: none
    {"rad_vet/time2_adm_poslim", ""}, {"rad_vet/time2_dtr_poslim", ""},
    {"rad_vet/zalesak_cap_rounds", ""}, {"rad_vet/vet_col_lat_posfix", ""},
    {"rad_vet/time2_vet_x2_fields", ""},
    // cell-layout keys: removed
    {"rad_vet/clean_alpha", "x"}, {"rad_vet/clean_recon", "x"},
    {"rad_vet/clean_face_model", "x"}, {"rad_vet/clean_consistency_weight", "x"},
    {"rad_vet/clean_psi_eps", "x"}, {"rad_vet/clean_pos_*", "x"},
    {"rad_vet/clean_diag_cells", "x"}, {"rad_vet/clean_maxpos", "x"},
    {"rad_vet/clean_vfreeze", "x"}, {"rad_vet/clean_gas_exact", "x"},
    {"rad_vet/clean_yres", "x"}, {"rad_vet/clean_nl_efloor", "x"},
    {"rad_vet/clean_pc", "x"}, {"rad_vet/clean_pc_sweeps", "x"},
    {"rad_vet/clean_pc_factor*", "x"}, {"rad_vet/clean_pc_reuse", "x"},
    {"rad_vet/clean_dissip_power", "x"},
    {"rad_vet/clean_line_teams", "x"}, {"rad_vet/clean_line_solver", "x"},
    {"rad_vet/clean_line_chunk", "x"},
    {"rad_vet/clean_fv_*", "x"},   // FV rays: removed
    // explicit transport's face-flux and source options: only their old defaults
    {"rad_vet/thick_flux", "=none"}, {"rad_vet/ap_form", "=alpha2"},
    {"rad_vet/reconstruct", "=plm"}, {"rad_vet/split_vel", "=recon"},
    {"rad_vet/source_form", "=full"}, {"rad_vet/f_source", "=cell"},
    {"rad_vet/dbg_tensor", "=none"}, {"rad_vet/dbg_tensor_tilt", ""},
    // old solver's diagnostics, multi-rate, coupling split, K1 FCT: removed
    {"rad_vet/dbg_energy_audit", "!"}, {"rad_vet/dbg_cs_metric_selftest", "!"},
    {"rad_vet/implicit_det_trace", "!"}, {"rad_vet/dbg_hydro_off", "!"},
    {"rad_vet/time2_k1_fct", "!"}, {"rad_vet/implicit_mr_every", "!"},
    {"rad_vet/coupling_split", "!"},
    // Zalesak fixed-point and local-round variants: only off
    {"rad_vet/zalesak_converge", "!"}, {"rad_vet/zalesak_local_rounds", "!"},
    {"rad_vet/zalesak_colsweep", "!"},
    // old step history: never kept, its restart channels are skipped
    {"rad_vet/clean_skip_time2_hist", "=true"},
    // hesdirk2 start probe: only be
    {"rad_vet/time2_start", "=be"},
    // extrapolated vet_sc builds: only 1
    {"rad_vet/vet_sc_every", "=1"},
    // restart fingerprints (debug hook): removed
    {"rad_vet/dbg_rst_fp*", "x"},
    // explicit transport's sub-cycling: removed
    {"rad_vet/cfl_rad", ""}, {"rad_vet/subcycle", ""},
    // vet_gd twin and symmetric reference: inert here (with vet_gd, a twin named true
    // stops in rad_vet.cpp)
    {"rad_vet/vet_gd_twin*", ""},
    {"rad_vet/vet_gd_symref*", ""},
    {"rad_vet/vet_gd_async", "!"}, {"rad_vet/time2_stage_theta", "!"},
    {"rad_vet/time2_stage_local_be", "!"},
    // study keys: removed
    {"rad_vet/vet_gd_smooth", "!"}, {"rad_vet/vet_dump", "=none"},
    {"rad_vet/vet_dump_every", ""}, {"rad_vet/vet_gd_halo_cache_mb", "!"},
    {"rad_vet/vet_gd_beam", "!"}, {"rad_vet/vet_gd_beam_*", "x"},
    {"rad_vet/vet_gd_beam2_*", "x"}, {"rad_vet/vet_gd_mms", "!"},
    {"rad_vet/vet_gd_hr_split", "=none"}, {"rad_vet/vet_gd_hr_split_*", ""},
    {"rad_vet/vet_gd_hr_latfit", "=-1"}, {"rad_vet/vet_gd_hr_latfit_*", ""},
    {"rad_vet/sp_face_sync", "!"}, {"rad_vet/sp_face_sync_relax", ""},
    // old solver's diagnostics and vet_col build keys: no effect
    {"rad_vet/report_newton_fb", ""}, {"rad_vet/dbg_t2_admiss", ""},
    {"rad_vet/dbg_opac_part", ""}, {"rad_vet/doppler_vel_dbg", ""},
    {"rad_vet/dbg_trans_memory", ""}, {"rad_vet/dbg_cell_lo", ""},
    {"rad_vet/dbg_cell_hi", ""}, {"rad_vet/vet_col_lat_posfix_force", ""},
    {"rad_vet/vet_col_lat_dump", ""}, {"rad_vet/vet_col_lat_dump_every", ""},
    {"rad_vet/polar_ring_avg_lump", "!"}, {"rad_vet/polar_ring_avg_gz", "!"},
    {"rad_vet/time2_vstage", ""}, {"rad_vet/vet_col_order2", ""},
    {"rad_vet/vet_col_team", ""}, {"rad_vet/vet_col_team_size", ""},
    {"rad_vet/vet_col_chunk", ""}, {"rad_vet/vet_col_spread", ""},
    // vet_sc on several MeshBlocks: ray kernel, angle sweep and lag removed
    {"rad_vet/vet_mb_kernel", "x"}, {"rad_vet/vet_mb_tblock*", "x"},
    {"rad_vet/vet_mb_mom_fuse", "x"}, {"rad_vet/vet_mb_mom_batch", "x"},
    {"rad_vet/vet_mb_halo", "x"}, {"rad_vet/vet_mb_overlap", "x"},
    {"rad_vet/vet_mb_angles", "x"}, {"rad_vet/vet_mb_lag", "x"},
    {"rad_vet/vet_mb_agroup", "=1"}, {"rad_vet/vet_mb_agg", "=true"},
    // vet_sc Milne diagnostic: removed
    {"rad_vet/vet_milne", "!"},
    // vet_path = cells: removed (decks and restart files store none)
    {"rad_vet/vet_path", "=none"},
    // vet_seg_cross: removed
    {"rad_vet/vet_seg_cross", "!"},
    // vet_col_lat depth cut: the clean rays always start at the first shell
    {"rad_vet/vet_col_lat_taucut", ""},
    // vet_gd upwind profile shift and dense shell launch: removed
    {"rad_vet/vet_gd_upwind", "=bilinear"}, {"rad_vet/vet_gd_upwind_*", ""},
    {"rad_vet/vet_gd_shell_list", "=true"},
    // cs sweep's pipelined direction groups: removed
    {"rad_vet/vet_gd_cs_xgroups", "=1"},
    // inward-only band ring slabs: removed
    {"rad_vet/vet_gd_band_ring_inward", "!"},
    // staggered increment and cross-flux study forms: removed
    {"rad_vet/clean_stag_kpos", "!"}, {"rad_vet/clean_stag_kpat", "=false"},
    {"rad_vet/clean_stag_kreal", "!"}, {"rad_vet/clean_stag_xsweep", "!"},
    {"rad_vet/clean_stag_xpat", "!"}, {"rad_vet/clean_thin_donor", "!"},
    {"rad_vet/clean_stag_stencil", "=true"}, {"rad_vet/clean_stag_sbimp", "!"},
    {"rad_vet/clean_lz_bnd_free", "!"}, {"rad_vet/clean_stag_thinf_held", "!"},
    {"rad_vet/clean_face_hclat", "!"},
    {"rad_vet/clean_ray_source_relaxed", "!"},
    // removed, named by decks or restart files: seam clamp, pass print, fixed set's band
    // reach (covered by the bound's band); halo depth cap and ray cut only at kept values
    {"rad_vet/cs_cross_mono", ""}, {"rad_vet/clean_nl_print", ""},
    {"rad_vet/vet_gd_reach", ""}, {"rad_vet/clean_halo_depth", "=0"},
    {"rad_vet/clean_ray_cut", "=true"}};
  const bool say = (global_variable::my_rank == 0);
  std::string none;
  std::vector<std::string> drop;
  for (InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_vet" && ib.block_name != "time") continue;
    for (InputLine &ln : ib.line) {
      const std::string k = ib.block_name + "/" + ln.param_name, &v = ln.param_value;
      for (const auto &t : tab) {
        const std::string p(t[0]);
        const size_t n = p.size() - (p.back() == '*');
        if (k.compare(0, n, p, 0, n) != 0 || (n == p.size() && k.size() != n)) continue;
        if (t[1] == nullptr) break;
        // handled here: never reported as unread (a default-echo line: no message)
        ln.accessed = true;
        if (ln.param_comment == ParameterInput::default_comment) break;
        const std::string c(t[1]);
        if (c == "x") {
          CleanFatal("<" + k + "> belonged to a removed path (the cell layout, the FV "
                     "rays, a debug hook removed or a retired study key); remove it "
                     "from the input");
        }
        if (c[0] == '=' && v != c.substr(1) && !(c == "=false" && v == "0")) {
          CleanFatal("<" + k + "> = " + v + ": removed from this build (only " +
                     c.substr(1) + " is accepted)");
        }
        if (c == "!" && v != "false" && v != "0" && v != "off") {
          CleanFatal("<" + k + "> = " + v + " has no counterpart in scheme = clean");
        }
        if (c.empty() || c == "!" || c[0] == '=') none += " " + k;   // (accepted value)
        if (c != "!" && c[0] != '=' && k.compare(0, 17, "rad_vet/implicit_") == 0) {
          drop.push_back(ln.param_name);
        }
        if (say && !c.empty() && c != "!" && c[0] != '=') {
          std::cout << "### WARNING (scheme = clean): <" << k << "> has no effect, the "
                    << "clean solve reads <rad_vet>/" << c << std::endl;
        }
        break;
      }
    }
  }
  if (say && !none.empty()) {
    std::cout << "### WARNING (scheme = clean): no effect:" << none << std::endl;
  }
  // implicit_* keys the clean solve does not read leave a fresh run's input, so
  // ImplicitInit parses its defaults; a restart keeps its file's keys (its layout)
  for (InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_vet" || global_variable::restart_run) continue;
    ib.line.remove_if([&](const InputLine &ln) {
      return std::find(drop.begin(), drop.end(), ln.param_name) != drop.end();
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanInit
//! \brief allocates the arrays of the clean solve and checks that the setup is supported

void RadiationVet::CleanInit() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb = pmy_pack->nmb_thispack;
  // clean_closure = iterate: not the closure a restart file holds (Time2RstSet)
  if (!cln_cip) {
    Kokkos::realloc(cln_s, nmb, CleanNs(pmy_pack->pmesh->use_cubed_sphere, cln_rfh),
                    n3, n2, n1);
  }
  if (cln_phic.extent_int(0) != nmb) {
    Kokkos::realloc(cln_phic, nmb, 3, n3, n2, n1);
  }
  // staggered on the cubed sphere: cell terms and faces of CleanStagFSolve
  const bool stcs = pmy_pack->pmesh->use_cubed_sphere;
  Kokkos::realloc(cln_v, nmb, stcs ? 18 : 3, n3, n2, n1);
  // mem_arena: the shared slots follow cln_w.o (CleanWoBind), unless the closure is built
  // inside the solve (clean_closure = iterate | extrap*)
  cln_w.s = DvceArray5D<Real>();
  if (!marena_on || cln_closit || cln_extrap > 0) {
    Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
  }
  const int nwo = (cln_xte > 0) ? CLW_NX : CLW_N;
  // mem_arena: cln_w.o is made at the first solve, over the intensities (CleanWoBind)
  cln_wo_nc = nwo - CLW_NS;
  cln_w.o = DvceArray5D<Real>();
  if (!marena_on) Kokkos::realloc(cln_w.o, nmb, cln_wo_nc, n3, n2, n1);
  cln_w.ns = CLW_NS;
  const int nfc = (cln_fsolve && cln_fkin) ? CLF_KD + 2 :
                  (cln_vimp ? CLF_DQ + 1 : CLF_QL + 3);
  // clean_stag_opac_update: only for T-dependent opacities (tables, powerlaw T exponent)
  cln_opt = cln_opt && !opac_zero && (opacity_type == VET_OPAC_TABLE ||
            (opacity_type != VET_OPAC_CONST && opac_b != 0.0));
  // staggered bracket's lagged flux DG: one more slot (CleanStagBracket)
  cln_fdg = cln_consist ? nfc : -1;
  cln_xqs = cln_xd ? nfc + (cln_fdg >= 0) : -1;   // clean_stag_xdonor's q
  int nfx = nfc + (cln_fdg >= 0) + (cln_xd ? 5 : 0);
  cln_fwbs = cln_fwbst ? nfx++ : -1;
  Kokkos::realloc(cln_fc1, nmb, nfx, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, nfx, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, nfx, n3, n2, n1);
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  Kokkos::deep_copy(cln_fc1, 0.0);
  Kokkos::deep_copy(cln_fc2, 0.0);
  Kokkos::deep_copy(cln_fc3, 0.0);
  for (int d = 0; d < 3 && cln_rstd.extent_int(0) == nmb && cln_fsolve; ++d) {
    auto fc = (d == 0) ? cln_fc1 : ((d == 1) ? cln_fc2 : cln_fc3);   // from a restart
    auto AL = Kokkos::ALL;
    Kokkos::deep_copy(Kokkos::subview(fc, AL, CLF_N + 3, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, d, AL, AL, AL));
    Kokkos::deep_copy(Kokkos::subview(fc, AL, CLF_KD, AL, AL, AL),
                      Kokkos::subview(cln_rstd, AL, 3 + d, AL, AL, AL));
  }
  if ((force_ref == VET_FREF_CLEAN || force_ref == VET_FREF_WB_MD) &&
      cln_fref.extent_int(0) != nmb) {   // unless a restart has read it
    Kokkos::realloc(cln_fref, nmb, (force_ref == VET_FREF_WB_MD) ? 9 : 3 + cln_lzt, n3,
                    n2,
                    n1);
  }
  Kokkos::realloc(cln_st, nmb, CLN_ST, indcs.nx3, indcs.nx2, indcs.nx1);
  if (pmy_pack->pmesh->use_cubed_sphere) {
    Kokkos::realloc(cln_csy, nmb, CLN_CSY, indcs.nx3, indcs.nx2, indcs.nx1);
  }
  if (cln_ysum) Kokkos::realloc(cln_yl, nmb, 9, n3, n2, n1);
  // multigrid levels of the preconditioner: halve x2 and x3 where the tile's are even
  // (>= 4), clean_stag_mg_levels coarse levels at most (x1 is never halved)
  std::array<int, 3> nc = {indcs.nx1, indcs.nx2, indcs.nx3};
  if (cln_mgt[0] == 0) cln_mgt = nc;   // layout_tile: the tiles' cells
  std::array<int, 3> tc = cln_mgt;
  while (static_cast<int>(cln_mgn.size()) < cln_mgnl) {
    bool coarser = false;
    for (int d = 1; d < 3; ++d) {
      if (tc[d] >= 4 && tc[d] % 2 == 0) {
        nc[d] /= 2;
        tc[d] /= 2;
        coarser = true;
      }
    }
    if (!coarser) break;
    cln_mgn.push_back(nc);
    cln_mg.emplace_back("cln_mg", nmb, CLS_N, nc[2], nc[1], nc[0]);
  }
  if (cln_extrap > 0) {
    // a restart may have read t^n and the end of the last step
    const int ns = cln_s.extent_int(1), c3 = indcs.nx3, c2 = indcs.nx2, c1 = indcs.nx1;
    if (cln_mn.extent_int(0) != nmb) Kokkos::realloc(cln_mn, nmb, ns, c3, c2, c1);
    if (cln_mc.extent_int(0) != nmb) Kokkos::realloc(cln_mc, nmb, ns, c3, c2, c1);
    if (!CleanExtrapSx()) Kokkos::realloc(cln_mp, nmb, ns, c3, c2, c1);
    if (cln_extrap == 2) Kokkos::realloc(cln_m1, nmb, ns, c3, c2, c1);
  }
  Mesh *pm = pmy_pack->pmesh;
  if (pm->multilevel || chat != c_light) {
    CleanFatal("<rad_vet>/scheme = clean currently needs a uniform mesh and chat = c");
  }
  // boundaries: open (bath, beam; outflow without a bath), reflect, periodic and an
  // imposed x1 flux.  The imposed flux is the face's whole energy flux, taken with no E
  // weights (CleanFaceWeights, CleanFaceSource), kept as the face F by the F solve; the
  // inner half cell absorbs its momentum (CleanMomentumChange), as at open faces.
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  for (int s = 0; s < 6; ++s) {
    const bool flux = bc[s] == VET_IBC_FLUX || (bc[s] == VET_IBC_X2FLUX && ifx2_on);
    if (!(bc[s] == VET_IBC_MARSHAK || bc[s] == VET_IBC_OUTFLOW || bc[s] == VET_IBC_REFLECT
          ||
          bc[s] == VET_IBC_PERIODIC || (flux && s < 2))) {
      CleanFatal("<rad_vet>/implicit_bc_* with scheme = clean: auto | marshak | outflow "
                 "| "
                 "reflect | periodic, flux (x1); not efix");
    }
  }
  if (iebath_hi_col_on) {
    CleanFatal("scheme = clean has no per-column outer bath");
  }
  if ((force_ref == VET_FREF_WB_ARAD || force_ref == VET_FREF_WB_MD) && !fref_wsplit) {
    CleanFatal("<rad_vet>/force_reference = wb_arad(_md) with scheme = clean needs the "
               "reference work given by the problem generator (force_reference_work = "
               "split)");
  }
  if (force_ref == VET_FREF_WB_MD && !cln_steady && !cln_fref_set) {
    CleanFatal("<rad_vet>/force_reference = wb_arad_md needs clean_init = steady (its "
               "reference is the steady start's force)");
  }
  if ((vgd_on || (vsc_rot && cln_ring)) && VgdRotOn()) {
    // rotating rays (CleanFvSets): the closure is the J-weighted mean of the sets, made
    // continuous across blocks by the cross-fade, no slope correction.
    // vet_rotate_avg_ema: only with clean_rot_ring = live and without _dema, _tau, _lc,
    // _diag (CleanRingEma)
    const bool ema_ok = cln_ring == 1 && !vra_ema_dema && vra_ema_tau <= 0.0 &&
                        vra_ema_lc <= 0.0 && vra_ema_diag == 0;
    if ((vra_ema && !ema_ok) || vra_eq ||
        (vgd_xf && !(vgd_rot_dt > 0.0 || vgd_rot_auto)) ||
        (vgd_rot3 && !pm->three_d)) {
      CleanFatal("<rad_vet>/scheme = clean with a rotating vet_gd set: vet_rotate_avg = "
                 "K "
                 "(J-weighted; _ema only with clean_rot_ring = live, no _dema, _tau, "
                 "_lc, _diag), vet_gd_rotate_xfade with vet_gd_rotate_dt, and on a "
                 "2-D axisymmetric mesh the z axis");
    }
  }
  const bool ring_ok = (vgd_on || vsc_rot) && vra_on && g_vra_adv &&
                       (vgd_rot_dt > 0.0 || vgd_rot_auto) && vgd_rot_align &&
                       vra_k <= VET_RA_NF;
  if (cln_ring && !(ring_ok && vgd_xf == (cln_ring == 1))) {
    CleanFatal("<rad_vet>/clean_rot_ring needs vet_gd or vet_sc rays rotating in time: "
               "vet_rotate_avg = K with vet_rotate_advance, vet_gd_rotate_dt and "
               "_align, with vet_gd_rotate_xfade for live, without for jump");
  }
  // sc rays sweep one set at a time: a mean or cross-fade of sets needs the ring
  if (vgd_on && VgdRotOn() && (vra_on || vgd_xf) && !cln_ring) {
    CleanFatal("<rad_vet>/clean_rays = sc with vet_rotate_avg or vet_gd_rotate_xfade "
               "needs clean_rot_ring");
  }
  // jump: a block's first step is backward Euler (no slope of the old sets)
  if (cln_ring == 2) vgd_rot_slope = 2;
  if (pm->use_cubed_sphere) {
    // rays: vet_gd's (clean_rays = sc: its band, VetGdCsSweep) or the analytic sphere
    // (eddington: infinite radius); an even panel resolution keeps cell centres off the z
    // axis (CleanCsFrame)
    if (!(cln_sphere || cln_edd || vgd_on) || (cln_consist && !vgd_on) ||
        pm->mesh_indcs.nx2 % 2 != 0) {
      CleanFatal("<rad_vet>/scheme = clean on the cubed sphere needs "
                 "clean_closure = sphere | eddington (no clean_consistency) or closure = "
                 "vet_col with vet_gd (clean_rays = sc) and "
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
    // rays: vet_gd's fixed direction set swept from the first shell, or the analytic
    // sphere (clean_closure = sphere | eddington)
    if (!(cln_sphere || vgd_on || cln_edd) ||
        (cln_consist && !vgd_on)) {
      CleanFatal("<rad_vet>/scheme = clean on spherical polar needs clean_closure = "
                 "sphere | eddington (no clean_consistency) or closure = "
                 "vet_col with vet_gd");
    }
    // theta walls: vet_gd_theta_wall or the mesh's (rad_vet.cpp); the sweep starts at the
    // first shell or the ray cut (vlat_icut).  No clean solve reads the half-range ratios
    // vgd_hr, so implicit_flux_beam = halfrange has no effect.
    impl_beam_hr = false;
    cln_init = true;
    return;
  }
  if (cln_sphere || !vet_sc || !vet_full) {
    CleanFatal("<rad_vet>/scheme = clean on a Cartesian mesh needs closure = vet_sc with "
               "vet_tensor = full (clean_closure = sphere is for spherical polar)");
  }
  // the formal solution provides the half-range sums and C1
  impl_beam_hr = true;
  thsw_c1 = true;
  cln_init = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanCsGeometry
//! \brief cubed sphere: the cell geometry of the finite-volume divergence in Cartesian y.
//! Lateral face: in the great-circle plane of its edge points P, P' on the cube surface,
//! normal along P x P', area (r_r^2 - r_l^2) th/2 (th the angle P, P'), first moment
//! (r_r^3 - r_l^3)/3 2 sin(th/2) along the bisector.  Radial face: r_f^2 S (S the vector
//! solid angle) minus half the sum of outward th n of the lateral faces, so a cell's six
//! vector areas sum to zero to round-off; centroid (r_r^4 - r_l^4) S/(4 V).
//! Cube points come from the global panel index (exactly +-1 on edges), so both panels at
//! a seam build a face from the same numbers, bitwise.  The frame of y is the hydro's:
//! r_hat and the panel tangents e_xi, e_eta.  Centroids, r_hat, the owner's x2/x3
//! normals, volume and gid (CLG_F, for CleanCsSlot) take the halo exchange.

void RadiationVet::CleanCsGeometry() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, ks = indcs.ks, nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*indcs.ng, n2 = indcs.nx2 + 2*indcs.ng;
  const int n3 = indcs.nx3 + 2*indcs.ng, np = pmy_pack->pmesh->mesh_indcs.nx2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  for (int d = 0; d < 3; ++d) {
    Kokkos::realloc(cln_ga[d], nmb1 + 1, n3, n2, n1);
    Kokkos::realloc(cln_gd[d], nmb1 + 1, n3, n2, n1);
  }
  // all CLG_NV in gn first (exchanges, distances), then into cln_gn's two arrays
  DvceArray5D<Real> gn("cln_csgn", nmb1 + 1, CLG_NV, n3, n2, n1);
  auto ga0 = cln_ga[0], ga1 = cln_ga[1], ga2 = cln_ga[2];
  auto gd0 = cln_gd[0], gd1 = cln_gd[1], gd2 = cln_gd[2];
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
  CleanGhostFill(gn, CLG_X, 6, false);
  CleanGhostFill(gn, CLG_F, 11, false);
  par_for("cln_csq", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int b = 0;
    for (int d = 1; d < 3; ++d) {
      for (int lo = 0; lo < 2 && (lo == 0 || (d == 1 ? j : k) > 0); ++lo) {
        Real s;
        const int q = CleanCsSlotQ(gn, d, m, k, j, i, k - lo*(d == 2), j - lo*(d == 1),
                                   i, s);
        b += ((q - 1) + 2*(s > 0.0)) << (4*(d - 1) + 2*lo);
      }
    }
    gn(m,CLG_Q,k,j,i) = b;
  });
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
  // the directions of CleanGeo::a are bitwise the same along each x1 range of a line
  const int is = indcs.is, ie = indcs.ie;
  Kokkos::realloc(cln_gn.r, nmb1 + 1, CLG_NV - CLG_NA, n3, n2, n1);
  Kokkos::realloc(cln_gn.a, nmb1 + 1, CLG_NA, n3, n2, 3);
  cln_gn.is = is;
  cln_gn.ie = ie;
  auto gr = cln_gn.r;
  auto gt = cln_gn.a;
  int nbad = 0;
  Kokkos::parallel_reduce("cln_csgeo2",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, DevExeSpace>(
      {0, 0, 0, 0}, {nmb1 + 1, n3, n2, n1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nb) {
    for (int n = 0; n < CLG_NV; ++n) {
      const int s = CleanGeo::Slot(n);
      const int c = (i < is) ? 0 : ((i > ie) ? 2 : 1), i0 = (c == 0) ? 0 : ((c == 1) ? is
                                                                        : ie + 1);
      const Real v = gn(m,n,k,j,i), v0 = gn(m,n,k,j,i0);
      if (s < 0) {
        gr(m,-1-s,k,j,i) = v;
      } else if (v != v0 || Kokkos::signbit(v) != Kokkos::signbit(v0)) {
        nb += 1;
      } else if (i == i0) {
        gt(m,s,k,j,c) = v;
      }
    }
  }, nbad);
  if (nbad > 0) {
    CleanFatal("<rad_vet> clean cubed sphere: " + std::to_string(nbad) +
               " geometry directions change along x1 (CleanGeo)");
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanSteadyStart
//! \brief <rad_vet>/clean_init = steady (fresh start, at the first dt): relaxes the
//! initial radiation to the initial matter and sources.
//! The gas is frozen (no hydro, no feedback) and backward-Euler radiation steps of
//! clean_init_dtfac x dt are taken until max(|dE|, |dF|/c) < clean_init_tol max E; as
//! dt -> infinity the fixed point is the scheme's own steady solution (closure included).
//! A failed pass (or E < 0) is undone and the step shrinks by 1000, down to dt.  The
//! passes then go on at the run's dt, to the tolerance or until the change contracts by
//! less than 1/2 per pass.  Only E and F change, except the gas internal energy with
//! force_reference = clean or icr_fmeas; the closure history is dropped afterwards.

void RadiationVet::CleanSteadyStart(Driver *pdrive) {
  if (!cln_steady) return;
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  const Real dt0 = pm->dt, cl = c_light;
  const bool sfb = gas_feedback;
  const int st2 = t2_solve;
  const Real nc0 = vet_ncall;
  // gas frozen, backward-Euler steps without substeps; big: the large step.  With
  // force_reference = clean or icr_fmeas the gas internal energy relaxes with E (rho,
  // momentum frozen) to radiative equilibrium; icr_fmeas also keeps each pass's x1 force
  // less the reference, dm - h rho a_ref, in icr_dm (icr_dm_dt = h, CleanWriteBack)
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool relax = (force_ref == VET_FREF_CLEAN || icr_fmeas) && sfb && fl.on;
  if (relax && icr_fmeas) {
    Kokkos::realloc(icr_dm, nmb1 + 1, u0.extent(2), u0.extent(3), u0.extent(4));
  }
  auto uh = fl.u0;
  gas_feedback = relax;
  cln_vzero = relax;
  cln_mdcap = (force_ref == VET_FREF_WB_MD);   // every pass's write-back: the last wins
  t2_solve = VET_T2S_NONE;
  nsub = 1;
  bool big = true, ok = false;
  int n2 = 0;
  // icr_fmeas: reach the large step by x10 per pass from dt (the gas solves fail when
  // started at the top)
  const bool lad = relax && icr_fmeas;
  Real dp = 0.0, fac = lad ? 1.0 : cln_steady_dtfac, emax = 0.0;
  cln_insteady = true;
  auto u0_ = u0;
  DvceArray5D<Real> uo("cln_steady_u", nmb1 + 1, 4, u0.extent(2), u0.extent(3),
                       u0.extent(4));
  DvceArray5D<Real> ug("cln_steady_g", relax ? nmb1 + 1 : 1, 4, u0.extent(2),
                       u0.extent(3), u0.extent(4));
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj*ni;
  const Real dv = MeshBlock0Dx(0)*MeshBlock0Dx(1)*MeshBlock0Dx(2);
  const bool sph = pm->use_spherical_polar || pm->use_cubed_sphere;   // volume sums
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
        sl += (sph ? vol(m,k,j,i) : dv)*u0_(m,VET_E+c,k,j,i);
      }, Kokkos::Sum<Real>(a));
      t[c] = CleanReduce(a, false);
    }
  };
  // u0 <-> uo with ghosts, plus face F and the write-back captures (cln_fref, icr_dm)
  DvceArray4D<Real> fo[4];
  DvceArray5D<Real> ro;
  auto copy = [&](const bool save) {
    DvceArray4D<Real> fv[4] = {f0x1, f0x2, f0x3, icr_dm};
    for (int d = 0; d < 4; ++d) {
      if (fo[d].size() != fv[d].size()) fo[d] = DvceArray4D<Real>("cln_steady_f",
                                                                  fv[d].layout());
      if (save) Kokkos::deep_copy(fo[d], fv[d]);
      if (!save) Kokkos::deep_copy(fv[d], fo[d]);
    }
    if (cln_mdcap && ro.size() != cln_fref.size()) {
      ro = DvceArray5D<Real>("cln_steady_r", cln_fref.layout());
    }
    if (cln_mdcap) Kokkos::deep_copy(save ? ro : cln_fref, save ? cln_fref : ro);
    Kokkos::parallel_for("cln_steady_copy", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         uo.size()),
    KOKKOS_LAMBDA(const int q) {
      const int n4 = uo.extent(4), n3 = uo.extent(3), n2 = uo.extent(2);
      const int i = q%n4, j = (q/n4)%n3, k = (q/(n4*n3))%n2, c = (q/(n4*n3*n2))%4;
      const int m = q/(n4*n3*n2*4);
      const int g = (c == 0) ? IEN : IM1 + c - 1;
      if (save) {
        uo(m,c,k,j,i) = u0_(m,VET_E+c,k,j,i);
        if (relax) ug(m,c,k,j,i) = uh(m,g,k,j,i);
      } else {
        u0_(m,VET_E+c,k,j,i) = uo(m,c,k,j,i);
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
    cln_tc = -1.0e300;   // predict this pass's closure, not the last dt's
    const Real nf = cln_nfail + cln_nneg;
    // one radiation step: the ordinary m1 task lists (closure, solve, ghost exchange)
    for (int stage = 1; stage <= nstage; ++stage) {
      pdrive->ExecuteTaskList(pm, "vet_before_stagen", stage);
      pdrive->ExecuteTaskList(pm, "vet_stagen", stage);
      pdrive->ExecuteTaskList(pm, "vet_after_stagen", stage);
    }
    (void) ApplyClosureLimits(pdrive, 1);
    if (cln_nfail + cln_nneg > nf) {
      copy(false);
      dtr_fail = false;
      cln_have_n = cln_have_p = cln_have_c = cln_have_1 = false;
      if (global_variable::my_rank == 0) {
        std::printf("<rad_vet> clean_init steady pass %d: solve failed or E < 0 at %.0e "
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
      Real x = fabs(u0_(m,VET_E,k,j,i) - uo(m,0,k,j,i));
      for (int c = 1; c < 4; ++c) {
        x = fmax(x, fabs(u0_(m,VET_E+c,k,j,i) - uo(m,c,k,j,i))/cl);
      }
      dl = fmax(dl, x);
      el = fmax(el, fmax(u0_(m,VET_E,k,j,i), uo(m,0,k,j,i)));
    }, Kokkos::Max<Real>(dm), Kokkos::Max<Real>(em));
    emax = fmax(emax, CleanReduce(em, true));
    d = CleanReduce(dm, true)/fmax(emax, 1.0e-300);
    if (global_variable::my_rank == 0) {
      std::printf("<rad_vet> clean_init steady pass %d: max(|dE|, |dF|/c)/max E = %.3e\n",
                  n, d);
    }
    if (!big && (d < cln_steady_tol || (++n2 > 2 && d > 0.5*dp))) break;
    if (d < cln_steady_tol && (!lad || fac >= cln_steady_dtfac)) big = false;
    if (lad) fac = fmin(10.0*fac, cln_steady_dtfac);
    ok = ok || !big;
    dp = d;
  }
  cln_insteady = false;
  // restore the run's setup; the first step starts from this state without history
  sums(t1);
  if (cln_mdcap) {   // wb_arad_md: the reference in the ghosts, read by the hydro walks
    CleanGhostFill(cln_fref, 0, 9, true, 0);
    if (pm->use_spherical_polar && pm->multi_d) {   // across a pole: no theta halves
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
      // cubed sphere: the lateral walks run along (e_xi, e_eta) and a seam ghost holds
      // the neighbour's axes, so each lateral half is a.(X_face - X_centre)
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
        for (int sd = 0; sd < 4; ++sd) {   // xi -, +, eta -, +
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
    std::printf("<rad_vet> clean_init steady: %s after %d passes, residual %.3e "
                "(tol %.1e; %.0f formal solutions); radiation E %.10e -> %.10e, F/c "
                "%.4e %.4e %.4e -> %.4e %.4e %.4e (volume sums)\n",
                ok ? "converged" : "NOT converged",
                std::min(n + 1, cln_steady_n), d, cln_steady_tol, vet_ncall - nc0, t0[0],
                t1[0], t0[1]/cl, t0[2]/cl, t0[3]/cl, t1[1]/cl, t1[2]/cl, t1[3]/cl);
    if (!ok) std::printf("CLEAN WARNING clean_init steady: not converged (residual %.3e, "
                         "tol %.1e)\n", d, cln_steady_tol);
  }
  // a last change that is not small is a diverged start, not an unconverged one
  if (!(d < 1.0e-2)) {
    CleanFatal("clean_init steady: the last pass changed E or F/c by " +
               std::to_string(d) + " of max E (diverged: clean_init_dtfac, _n)");
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::CleanSolve
//! \brief the implicit stage solve of the clean scheme, in the staggered layout (E in
//! cells, F on faces)

TaskStatus RadiationVet::CleanSolve(Driver *pdriver, int stage) {
  if (!cln_init) CleanInit();
  tau_ten = DvceArray5D<Real>();   // (an old file's, RssMode: cln_s from here on)
  // the vet_gd direction set first (a restart's first solve reads it before any build)
  if (vgd_on && !tau_ready) VetColInit();
  // sc rays' ring on the cubed sphere: the union's directions and weights too
  if (vgd_on && cln_cssc && cln_ring && cln_fvn == 0) {
    CleanFvSets(pmy_pack->pmesh->time, pmy_pack->pmesh->ncycle - 1);
    // sc ring: its two sets' intensities from the restart file (zero from a file of
    // another ray method, until the first build)
    if (cln_cssc) CleanScUnion2();
  }
  if (cln_cssc && !cln_ring) {   // clean_rays = sc: the sweep's set, its cln_fvi
    cln_fvd = vgd_dir;
    cln_fvn = vgd_n;
  }
  // restart from a file with the intensities (version 5): their face half ranges, set by
  // set
  if (cln_cssc && cln_fh.extent_int(0) == 0 && cln_fvi.extent_int(0) > 0) {
    for (int q = 0; q < cln_fvi.extent_int(1)/vgd_n; ++q) {
      CleanFaceHalf(cln_fvd, q*vgd_n, q == 0, q*vgd_n);
    }
    cln_fvi = DvceArray5D<Real>();
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  cln_feedback = gas_feedback && fl.on;

  // 1. starting face fluxes on the first cycle (restarts carry them)
  if (!cln_faces_set) {
    cln_faces_set = true;
    if (pmy_pack->pmesh->ncycle == 0) CleanInitFaceFluxes();
  }

  // 2. closure (from the gas before the stage increments), old vector, face coefficients.
  //    clean_closure = iterate with no closure in place takes the predicted one
  if (cln_extrap > 0) {
    CleanClosureExtrap();
  } else if (!cln_closit || !cln_cip) {
    CleanClosure();
  }
  if (cln_ring == 2) CleanRingSolveSet();
  CleanWoBind();
  cln_xdh = false;   // clean_stag_xdonor_freeze: no hold before this solve's passes
  CleanOldVector();
  CleanFaceCoefficients();
  // donors' E_old; the bracket's gate q
  if (cln_gflux || cln_consist) CleanGhostFill(cln_w, CLW_EO, 1);
  if (cln_consist) CleanStagBracket();
  // 3. Newton iteration on the gas coupling, Picard passes on the lagged slopes and (with
  //    clean_closure = iterate, at most clean_maxpic) the closure; one linear solve
  //    without these.  clean_closure = extrap(2): the stage-1 solve takes the closure of
  //    its own first iterate once (pass 1), else the step is first order in time.
  const bool clos1 = (cln_extrap > 0) && (t2_solve == VET_T2S_STAGE1);
  const bool iterate = cln_feedback || cln_closit || clos1;
  const bool t2stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  int nit = 0, nit_sum = 0, nsolve = 0;
  std::array<int, 2> nfirst = {0, 0};
  cln_ncwin = 0;
  bool converged = !iterate;
  const int maxnl = iterate ? cln_maxnl : 1;
  // lagged slopes settle only to the linear tolerance.  clean_conv = budget: the pass's
  // tolerance (CleanNewtonUpdate), the first solve's from the floor (restarts bitwise)
  Real tol_nl = fmax(cln_tol_nl, 10.0*cln_tol);
  cln_tolb = tol_nl;
  // budget norm and time error only in gas-coupled Newton passes with K1 = f(U^n) known;
  // otherwise the strict test
  cln_bon = cln_budget && iterate && cln_feedback && (t2_ok || cln_k1ok) &&
            t2_solve != VET_T2S_NONE && !cln_vzero &&
            t2k1.extent_int(0) >= pmy_pack->nmb_thispack;
  cln_ewr = -1.0;
  // a change that stopped decreasing within 100 x both its round-off level (cln_fe,
  // cln_fv) and tol_nl is converged too
  auto settled = [&](const Real d, const Real d0, const Real f) {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  };
  Real dnl = 0.0, dnl0 = -1.0;   // the last pass's change and the one before
  {
    CleanFirstIterate();
    cln_fsk = false;
    cln_fwbok = false;
    cln_csyok = false;   // lateral sums' coefficients once per closure
    if (cln_opt && cln_feedback) CleanStageOpacity();
    if (cln_feedback) CleanAdvectiveCoefficients();
    // the drag's face flux from the first iterate's beta and emission (a restart, failed
    // stage or discarded attempt does not carry the last solve's)
    if (cln_vimp && cln_feedback) CleanRightHandSide(true);
    if (cln_vimp && cln_feedback) CleanStagDrag(true);
    // clean_stag_xdonor: q of the first iterate with its drag share b (CLW_L holds E_k's
    // lateral sums after the first drag)
    CleanStagXDonor(true, !(cln_vimp && cln_feedback));
    converged = !iterate;
    cln_gpic = false;   // clean_gas_newton_switch: the solve's history
    std::vector<Real> rhist;
    Real de0 = 1.0e300, dv0 = 1.0e300, trel = 0.0;
    for (int inl = 0; inl < maxnl; ++inl) {
      const bool clos = (cln_closit && inl < cln_maxpic) || (clos1 && inl == 1);
      if (clos) {
        cln_csyok = false;
        // sdirk2 / trbdf2: fwb's s of the stage's own closure (pass 0's is the last
        // step's, CleanClosureExtrap)
        if (clos1 && t2_tab != 0) cln_fwbok = false;
        CleanClosureIterate();
        CleanFaceCoefficients();
        if (cln_consist) CleanStagBracket();
        if (cln_feedback) CleanAdvectiveCoefficients();
      }
      // opacities of the iterate's T, in passes < clean_stag_opac_passes
      const bool opr = cln_opt && cln_feedback && inl > 0 && inl < cln_optn;
      if (opr) {
        CleanStageOpacity();
        if (cln_optdr && cln_vimp) CleanStagDrag(false);
        CleanAdvectiveCoefficients();
      }
      // clean_stag_xdonor: q of the iterate; CLW_L holds E_k's lateral sums unless the
      // closure, bracket or opacity refresh rewrote them (and on the cubed sphere)
      if (inl > 0) CleanStagXDonor(false, clos || opr || !cln_feedback,
                                   cln_xdfz > 0 && inl >= cln_xdfz);
      // after the xdonor: the lateral sums of E_k are its
      CleanRightHandSide(false, !(inl > 0 && cln_xd_on));
      // clean_xt_cons: the first iterate's T from the pass's gas row, then its RHS
      if (inl == 0 && cln_xtcons && cln_feedback) {
        CleanNewtonUpdate();
        CleanRightHandSide();
      }
      CleanPreconditionerSetup();
      // budget: a pass whose iterate gives the next closure is solved in the strict norm
      const bool cnext = (cln_closit && inl + 1 < cln_maxpic) || (clos1 && inl == 0);
      if (cnext) cln_ewr = (cln_c1b && !dtr_redo) ? -3.0 : -2.0;
      nit = CleanBiCGStab();
      if (cnext) cln_ewr = -1.0;
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (nit < 0) dnl = std::numeric_limits<Real>::infinity();
      // implicit_lin_maxit: a linear solve that ends at its cap leaves its iterate to the
      // pass, which is then not the last one
      const bool lp = cln_linpass && nit < 0 && inl + 1 < maxnl;
      if (!iterate || (nit < 0 && !lp)) break;
      const bool fbl = lp;
      if (!t2stage && inl + 1 == maxnl) {   // for CleanWorstCell
        auto sv = Kokkos::subview(cln_w.s, Kokkos::ALL, std::make_pair(CLW_EK,
                                  CLW_EK + 1), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        Kokkos::realloc(cln_wci, sv.extent(0), 1, sv.extent(2), sv.extent(3),
                        sv.extent(4));
        Kokkos::deep_copy(cln_wci, sv);
      }
      const Real de = CleanNewtonUpdate();
      if (cln_bon) tol_nl = cln_tolb;
      Real dv = 0.0;
      // (clean_vel_settled: a pass that cannot end the solve keeps v and F; CLW_L the
      // lateral sums of E_k, as the velocity update leaves them)
      const bool vs = cln_vset && !lp &&
                      (!settled(de, de0, cln_fe) || (clos1 && inl == 0));
      if (cln_feedback && vs) CleanStagCsLateral(CLW_EK);
      if (cln_feedback && !vs) dv = CleanVelocityUpdate();
      dnl0 = (inl + 1 == maxnl) ? dnl : -1.0;
      dnl = fmax(de, dv);
      // clean_gas_newton_switch = W: a stalled solve (pass + 1 >= _min and > W, the
      // minimum over the last W passes >= clean_stall_fac x the one before, up at least
      // twice in the window) drops the gas Newton slope for its remaining passes: the gas
      // T takes the bracketed root of its row, a contraction (Newton cycles where c_v
      // jumps, e.g. H ionisation).  Deterministic in the passes: restarts bitwise.
      if (cln_gnsw > 0 && !cln_gpic) {
        rhist.push_back(dnl);
        const int nw = cln_gnsw;
        if (inl + 1 >= cln_gnswmin && inl + 1 > nw) {
          Real mwin = 1.0e300, mpre = 1.0e300;
          int nup = 0;
          for (int q = 0; q <= inl; ++q) {
            if (q > inl - nw) {
              mwin = fmin(mwin, rhist[q]);
              if (rhist[q] > rhist[q-1]) ++nup;
            } else {
              mpre = fmin(mpre, rhist[q]);
            }
          }
          if (mwin >= cln_stfac*mpre && nup >= 2) {
            cln_gpic = true;
            cln_ngnsw += 1.0;
            if (global_variable::my_rank == 0) {
              std::printf("CLEAN GNSW cycle %d t2s %d pass %d: stalled, the gas Newton "
                          "slope off (%.0f so far)\n", pmy_pack->pmesh->ncycle, t2_solve,
                          inl, cln_ngnsw);
            }
          }
        }
      }
      // strict test: stop at 0.1 x the first pass's change, at most 1e-4
      if (!cln_bon && inl == 0) trel = fmin(0.1*dnl, 1.0e-4);
      tol_nl = fmax(tol_nl, trel);
      // clean_verbose >= 2: each pass's changes against the tolerance it is tested with
      if (cln_verbose > 1 && global_variable::my_rank == 0) {
        std::printf("CLEAN PASS cycle %d t2s %d pass %d it %d de %.3e dv %.3e tol %.3e "
                    "fe %.3e fv %.3e%s%s\n", pmy_pack->pmesh->ncycle, t2_solve, inl,
                    nit, de, dv, tol_nl, cln_fe, cln_fv, clos ? " clos" : "",
                    opr ? " opac" : "");
      }
      if (fbl) {
        de0 = de;
        dv0 = dv;
        continue;
      }
      if (settled(de, de0, cln_fe) && settled(dv, dv0, cln_fv) && !(clos1 && inl == 0)) {
        converged = true;
        break;
      }
      de0 = de;
      if (!vs) dv0 = dv;
    }
    if (iterate) CleanCopyIterate();
    // R1, R2 of the solution, counted
    CleanStagRealizable();
    nfirst = {cln_nr1, cln_nr2};
  }
  if (!converged) {
    nit = -std::abs(nit) - 1;
  }

  // 4. positivity of the solution
  const int nneg = CleanCountNegative();
  if (nit < 0) cln_nfail += 1.0;
  if (nneg > 0) cln_nneg += 1.0;
  // an unconverged backward-Euler solve fails the cycle: with <time>/dt_retry it is
  // retaken at a smaller dt (a failed stage is redone by backward Euler via t2_fail);
  // without, the run stops if the last change exceeds implicit_resid_fatal
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
    CleanWorstCell(dnl, dnl0);
    if (fatal) CleanFatal("<rad_vet> clean: a backward-Euler solve did not converge and "
                          "its last change exceeds <rad_vet>/implicit_resid_fatal");
  }
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pos 0 mark 0 all 0 neg %d "
                "r1 %d r2 %d first %d %d guard %d solves %d cwin %d\n",
                pmy_pack->pmesh->ncycle, t2_solve, nit, nit_sum, nneg, cln_nr1, cln_nr2,
                nfirst[0], nfirst[1], cln_nguard, nsolve, cln_ncwin);
  }
  if (cln_mclc.extent_int(0) == 1) {   // last pass's face momentum terms at the cap
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_mclc);
    const Real n = CleanReduce(h(0), false);
    if (global_variable::my_rank == 0 && cln_verbose > 0 && n > 0.0) {
      std::printf(
          "CLEAN MCAP cycle %d t2s %d: %.0f face momentum terms at the cap h chi_f "
          "E in the last pass\n", pmy_pack->pmesh->ncycle, t2_solve, n);
    }
  }

  // 5. new state: face fluxes, E, cell fluxes, gas, stage slopes
  // a failed solve or a probe does not change the energy flux's memory
  cln_dok = (nit >= 0) && !t2_inprobe;
  // clean_wb_floor_fail: E floor energy not from the cell's gas in this write-back
  const bool wff = cln_wbff && !t2_inprobe;
  auto ftk = [&]() {
    auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pos_cnt_d);
    return hc(VET_POS_FLR_MV) + hc(VET_POS_FLR_UN);
  };
  const Real ft0 = wff ? ftk() : 0.0;
  CleanWriteBack();
  const Real ftw = wff ? CleanReduce(ftk() - ft0, false) : 0.0;
  const bool ffl = ftw > 0.0 && (t2stage || (dtr_on && dtr_ncut < dtr_maxcut));
  if (ftw > 0.0 && !t2stage && dtr_on && ffl) {
    if (!dtr_fail) {
      dtr_why = 8;   // efloor_column
      dtr_resid = ftw;
    }
    dtr_fail = true;
  }
  if (ftw > 0.0 && !ffl) cln_wbffn += 1.0;
  if (ftw > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN %s cycle %d t2s %d: the E floor took %.6e not from the cells' "
                "gas (%s)\n", ffl ? "WBFLOOR" : "WARNING", pmy_pack->pmesh->ncycle,
                t2_solve, ftw, ffl ? (t2stage ? "stage fails" : "dt_retry") :
                "accepted: last cut or no dt_retry");
  }
  // a write-back that leaves gas eint < 0 is rejected like a failed solve (stage redone
  // by backward Euler, a backward-Euler step by dt_retry)
  if (cln_feedback && cln_neint > 0 && !t2stage && dtr_on) {
    if (!dtr_fail) {
      dtr_why = 4;   // reason eint<=0
      dtr_resid = cln_neint;
    }
    dtr_fail = true;
  }
  if (cln_feedback && cln_neint > 0 && !t2stage && !dtr_on && !t2_inprobe &&
      global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: backward-Euler step accepted with %d "
                "cells of gas eint below 0 / its floor, E < 0 or not finite (no "
                "dt_retry)\n", pmy_pack->pmesh->ncycle, t2_solve, cln_neint);
  }
  if (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2) {
    // implicit_pos_floor(_cons): an E < 0 raised by the conservative write-back floor
    // does not fail the stage (counted)
    const bool flr = impl_pos_floor && impl_pos_flr_cons;
    t2_fail = (nit < 0) || (nneg > 0 && !flr) || (cln_feedback && cln_neint > 0) ||
              (t2_solve == VET_T2S_STAGE1 && pmy_pack->pmesh->ncycle == t2_dbg_fail) ||
              ffl;
    // time2_reject = retry: redo the cycle at a smaller dt instead of the backward-Euler
    // redo; after dt_retry_t2_be such cuts, the redo
    const bool rt = t2_fail && t2_reject == 1 && dtr_on && !t2_inprobe;
    if (rt && !dtr_t2be) {
      if (!dtr_fail) {dtr_why = 5; dtr_resid = -1.0;}
      dtr_fail = true;
    }
    if (rt && dtr_t2be) dtr_t2hit = true;
    if (cln_lzbook && cln_feedback && t2_solve == VET_T2S_STAGE2 && !t2_fail) {
      CleanLzBook(0.0, false);   // the accepted stage 2's L_z, pending until here
    }
  }
  // K1 = f(U^n) for the next predictor and time error: from an accepted backward-Euler
  // or stage-2 solve; a failed one leaves none (BE redo and auto-BE steps keep the last)
  if (t2_solve == VET_T2S_BESTORE || t2_solve == VET_T2S_STAGE2) {
    cln_k1ok = (nit >= 0) && !dtr_fail && !(t2_solve == VET_T2S_STAGE2 && t2_fail);
  }
  // force_reference = clean: the reference is the force of the converged, accepted first
  // coupled solve (each iterate recaptures it; a failed one is redone)
  if (force_ref == VET_FREF_CLEAN && cln_feedback && !cln_vzero && nit >= 0 &&
      nneg == 0 && !(t2stage && t2_fail) && !dtr_fail && !t2_inprobe) cln_fref_set = true;

  cln_xd_on = false;   // clean_stag_xdonor: q is of this solve's passes only

  // 6. clean_closure = extrap(2): ray moments of the new state
  if (cln_extrap > 0 && !t2_fail) CleanClosureAfterSolve();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanInitFaceFluxes
//! \brief face fluxes as the mean of the two adjacent cell fluxes

void RadiationVet::CleanInitFaceFluxes() {
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
    CleanFor("cln_f0", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      f0(m,k,j,i) = 0.5*(u0_(m,VET_F1+dir,k,j,i) + u0_(m,VET_F1+dir,k-dk,j-dj,i-di));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRstFaces
//! \brief staggered restart from a file without the clean state (the old scheme): its
//! faces hold the comoving normal flux F0 and its cells F = mean(F0) + a E, a = v + v.D,
//! while the staggered faces hold the lab-frame flux.  Each interior face adds the old
//! scheme's enthalpy flux a_f E_f (implicit_enthalpy: plm, central, or upwind donor) with
//! a = (F - mean F0)/E per cell (Cartesian on the cubed sphere, exchanged across seams);
//! physical-boundary faces keep theirs.

void RadiationVet::CleanRstFaces() {
  if (!cln_rstold) return;
  if (!cln_init) CleanInit();
  CleanWoBind();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto u0_ = u0;
  auto w = cln_w;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  auto cw = pmy_pack->pcoord->cos_cell, sw = pmy_pack->pcoord->sin_cell;
  const Real efl = e_floor;
  const int enm = (dop_face == 1) ? VET_IENTH_CENTRAL : impl_enth;
  const bool fre = cln_freal;
  DvceArray1D<Real> rc("cln_rstr", 2);
  // a of every cell (CLW_DM, 3) and its E (CLW_B), as the old scheme's write-back
  par_for("cln_rsta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = fmax(u0_(m,VET_E,k,j,i), efl);
    Real g[3] = {0.5*(f1(m,k,j,i) + f1(m,k,j,i+1)), 0.0, 0.0};
    if (x.nf > 1) g[1] = 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i));
    if (x.nf > 2) g[2] = 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i));
    w(m,CLW_B,k,j,i) = e;
    for (int d = 0; d < 3; ++d) w(m,CLW_DM+d,k,j,i) = (u0_(m,VET_F1+d,k,j,i) - g[d])/e;
    if (!x.cub) return;
    // covariant e_xi, e_eta components (face normals g = sin (a^xi, a^eta))
    const Real c = cw(m,k,j), sn = sw(m,k,j), s2 = 1.0 - c*c;
    const Real a2 = (u0_(m,VET_F2,k,j,i) - c*u0_(m,VET_F3,k,j,i))/s2 - g[1]/sn;
    const Real a3 = (u0_(m,VET_F3,k,j,i) - c*u0_(m,VET_F2,k,j,i))/s2 - g[2]/sn;
    const Real a1 = u0_(m,VET_F1,k,j,i) - g[0];
    for (int b = 0; b < 3; ++b) {
      w(m,CLW_DM+b,k,j,i) = (a1*x.g.At<CLG_R>(m,b,k,j,i) + a2*x.g.At<CLG_E>(m,b,k,j,i) +
                             a3*x.g.At<CLG_E>(m,3+b,k,j,i))/e;
    }
  });
  CleanGhostFill(cln_w, CLW_DM, 4);
  for (int dir = 0; dir < x.nf; ++dir) {
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("cln_rstf", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      // cell o = -2, -1 | 0, 1 of the face: its a along the face normal, and E
      auto a = [&](const int o) {
        return CleanStagN(x, CLW_DM, dir, m, k, j, i, k+o*dk, j+o*dj, i+o*di);
      };
      auto e = [&](const int o) {return w(m,CLW_B,k+o*dk,j+o*dj,i+o*di);};
      Real af = 0.5*(a(-1) + a(0)), ef = 0.5*(e(-1) + e(0)), dum;
      if (enm == VET_IENTH_UPWIND) {
        ef = (af > 0.0) ? e(-1) : e(0);
        af = (af > 0.0) ? a(-1) : a(0);
      } else if (enm >= VET_IENTH_PLM &&
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
      // clean_stag_freal: realisable, |F| <= c E of its donor.  In an old E-floor hole
      // the file's F, F0 are not those of its E (the old scheme would have clipped them),
      // so the excess is dropped: counted and its momentum sum |dF| V/c^2 reported
      const Real fc = fmin(fmax(f0(m,k,j,i), -x.c*e(0)), x.c*e(-1));
      if (fre && fc != f0(m,k,j,i)) {
        const Real v = CleanVol(x, m, k, j, i);
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
      std::cout << "<rad_vet> clean_stag_freal RESTART CONVERSION (one time, not a "
                << "deposit): faces of the old scheme's file above c E_donor taken at it="
                << n << ", radiation momentum removed sum |dF| V/c^2=" << p
                << " (its energy equivalent c sum |dF| V/c^2=" << p*x.c
                << "; E unchanged)" << std::endl;
    }
  }
  // the cell flux of these faces, as CleanWriteBack stores it (cubed sphere: the frame of
  // y, not the old covariant panel components)
  par_for("cln_rstc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), (x.nf > 1) ? f2(m,k,j,i) : 0.0,
                       (x.nf > 1) ? f2(m,k,j+1,i) : 0.0, (x.nf > 2) ? f3(m,k,j,i) : 0.0,
                       (x.nf > 2) ? f3(m,k+1,j,i) : 0.0};
    Real fv[3];
    CleanStagCellVector(x, m, k, j, i, c, fv);
    CleanStagHydro(x, m, k, j, i, fv);
    for (int d = 0; d < 3; ++d) u0_(m,VET_F1+d,k,j,i) = fv[d];
  });
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> clean restart from a file without the clean state: the faces "
              << "take its enthalpy flux a E (lab-frame F)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanOldVector
//! \brief E and gas state the stage starts from (stage increments added), with the gas
//! internal energy and temperature

void RadiationVet::CleanOldVector() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  const bool feedback = cln_feedback, cry = cln_eicry;
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fluid = fl.on;
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto u0_ = u0;
  auto inc = t2inc;
  auto w = cln_w;
  CleanArgs x = CleanMakeArgs();

  // counted: gas eint < 0 at the stage start (taken as 0 by the exchange), E < 0, and a
  // non-finite start (a NaN max(e0, 0) would hide: the solve fails, the step is redone).
  // A flat range; the counts are integers, exact in any order.
  Real nn[3] = {0.0, 0.0, 0.0};
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("cln_old", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &ng, Real &ne, Real &nf) {
    const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
    const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
    Real t = 0.0;
    if (fluid) {
      if (feedback && stage) {
        for (int n = 0; n < 3; ++n) uh(m,IM1+n,k,j,i) += inc(m,VET_T2_M1+n,k,j,i);
        uh(m,IEN,k,j,i) += inc(m,VET_T2_EN,k,j,i);
      }
      const Real d = uh(m,IDN,k,j,i);
      Real p[3];
      CleanGasMomentum(x, uh, m, k, j, i, p);
      const Real ekin = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/d;
      // clean_stag_eicarry: a deficit e0 < 0 stays in the gas row, paid by its implicit
      // exchange (e(T) > 0); else taken as 0, the deficit in CLW_JV for the write-back
      const Real e0 = CleanGasEint(x, uh, ekin, m, k, j, i);
      const Real eint = cry ? e0 : fmax(e0, 0.0);
      if (e0 < 0.0) ng += 1.0;
      if (!(Kokkos::isfinite(e0) && Kokkos::isfinite(d))) nf += 1.0;
      t = eos.Temperature(d, fmax(eint, 0.0));
      w(m,CLW_EI,k,j,i) = eint;
      w(m,CLW_JV,k,j,i) = e0 - eint;   // the gas's deficit below 0: CleanWriteBack
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,VET_E,k,j,i) + (stage ? inc(m,VET_T2_E,k,j,i) : 0.0);
    if (w(m,CLW_EO,k,j,i) < 0.0) ne += 1.0;
    if (!Kokkos::isfinite(w(m,CLW_EO,k,j,i))) nf += 1.0;
  }, Kokkos::Sum<Real>(nn[0]), Kokkos::Sum<Real>(nn[1]), Kokkos::Sum<Real>(nn[2]));
  nn[0] = CleanReduce(nn[0], false);
  nn[1] = CleanReduce(nn[1], false);
  nn[2] = CleanReduce(nn[2], false);
  // clean_cs_lz_abs = 2: the increment's L_z, of the slopes' sums (Time2FormStage)
  if (cln_lzbook && fluid && feedback && stage) CleanLzBook(cln_lzc, true);
  if (global_variable::my_rank == 0 && nn[2] > 0.0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage start not finite in %.0f cells "
                "(gas or E; the solve fails and the step is redone)\n",
                pmy_pack->pmesh->ncycle, t2_solve, nn[2]);
  }
  if (global_variable::my_rank == 0 && cln_verbose > 0 && (nn[0] > 0.0 || nn[1] > 0.0)) {
    std::printf("CLEAN OLDVEC cycle %d t2s %d eint<0 %.0f (%s) E<0 %.0f\n",
                pmy_pack->pmesh->ncycle, t2_solve, nn[0], cry ? "carried" : "taken as 0",
                nn[1]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagXDonor
//! \brief clean_stag_xdonor (staggered): the off-diagonal stress flux of each face,
//! C = -(1-w) theta h c^2 d_t(D_nt E) (CleanFaceCross), of the iterate E_k enters the
//! pass as a donor-cell coefficient, C(E) -> q E_donor, q = C(E_k)/E_k,donor, |q| <= c.
//! q adds to al (>= 0) or ar (<= 0) only, so the operator stays an M-matrix for any D, h:
//! a non-negative right-hand side gives a non-negative solution, without a mask.  At the
//! fixed point q E_donor = C(E) where |C| <= c E_donor; beyond, the cross flux is limited
//! at free streaming.  The gas momentum keeps the central stress.
//! Lagged C alternates pass to pass where nothing else holds E (optically thin), so a
//! face whose change turns back, d = l x the last (l < 0), takes d/(1 - l), the fixed
//! point of that mode.  first: a solve's first iterate, no history.

void RadiationVet::CleanStagXDonor(const bool first, const bool lat, const bool hold) {
  if (!cln_xd || !pmy_pack->pmesh->multi_d) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  cln_xd_on = false;
  CleanGhostDefer(cln_w, CLW_EK, 1);   // one exchange
  CleanGhostFill(cln_w, CLW_EO, 1);   // donors' E_old of the lag split
  if (lat) CleanStagCsLateral(CLW_EK);
  // hold (clean_stag_xdonor_freeze): keep C, q and the explicit part of the last pass
  cln_xdh = hold;
  if (hold) {cln_xd_on = true; return;}
  CleanArgs x = CleanMakeArgs();   // x.xq = -1: the central cross flux
  const int xq = cln_xqs;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    // zero every face (one flat range); the pass's split (CleanStagXSplit) sets the
    // active ones
    const int n4 = fc.extent_int(4), n34 = fc.extent_int(3)*n4;
    const int n234 = fc.extent_int(2)*n34;
    AthFor("cln_xdonor", Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0,
                         (nmb1 + 1)*n234),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/n234, k = (idx - m*n234)/n34;
      const int j = (idx - m*n234 - k*n34)/n4, i = idx%n4;
      for (int q = 0; q < 3; ++q) fc(m,xq+q,k,j,i) = 0.0;
      if (k < ks || k > ke+dk || j < js || j > je+dj || i < is || i > ie+di) return;
      const Real cn = CleanFaceCross(x, CLW_EK, dir, m, k, j, i);
      const Real cp = fc(m,xq+3,k,j,i), d = cn - cp, dp = fc(m,xq+4,k,j,i);
      // l < 0 where C turns back: step d/(1 - l), the fixed point of a mode C -> l C;
      // continuous in l
      const Real l = (d*dp < 0.0) ? d/dp : 0.0;
      const Real cu = first ? cn : cp + d/(1.0 - l);
      fc(m,xq+3,k,j,i) = cu;
      fc(m,xq+4,k,j,i) = first ? 0.0 : cu - cp;
    });
  }
  cln_xd_on = true;
}

#endif  // VET_TU_PART(0)

#if VET_TU_PART(1)
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanWorstCell
//! \brief CLEAN WORST diagnostic of an unconverged backward-Euler solve: its kind from
//! the last two changes (STALL: fell less than 2x) and the cell of the largest |dE|/S

void RadiationVet::CleanWorstCell(const Real dnl, const Real dnl0) {
  const bool ok = std::isfinite(dnl) && dnl0 >= 0.0 && cln_bon;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int ni = id.nx1, nji = id.nx2*ni, nkji = id.nx3*nji, is = id.is, js = id.js;
  const int ks = id.ks;
  auto w = cln_w;
  auto ci = cln_wci;
  auto op = opac;
  Kokkos::MaxLoc<Real,int>::value_type ml;
  Kokkos::parallel_reduce("cln_worst", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          ok ? pmy_pack->nmb_thispack*nkji : 0),
  KOKKOS_LAMBDA(const int q, Kokkos::MaxLoc<Real,int>::value_type &v) {
    const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/ni, i = is + q%ni;
    const Real d = fabs(w(m,CLW_EK,k,j,i) - ci(m,0,k,j,i))/w(m,CLW_SE,k,j,i);
    if (d > v.val) {v.val = d; v.loc = q;}
  }, Kokkos::MaxLoc<Real,int>(ml));
  const Real gmx = CleanReduce(ml.val, true);
  if (ok ? !(ml.loc >= 0 && ml.val >= gmx) : global_variable::my_rank != 0) return;
  const int q = std::max(ml.loc, 0), m = q/nkji, k = ks + (q%nkji)/nji;
  const int j = js + (q%nji)/ni, i = is + q%ni;
  DvceArray1D<Real> cv("cln_worst", 5);
  if (ok) {
    par_for("cln_worst_cell", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      cv(0) = w(m,CLW_EK,k,j,i);
      cv(1) = w(m,CLW_EO,k,j,i);
      cv(2) = w(m,CLW_T,k,j,i);
      cv(3) = w(m,CLW_SE,k,j,i);
      cv(4) = op(m,VET_OP_T,k,j,i);
    });
  }
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cv);
  std::printf("CLEAN WORST cycle %d %s last change %.3e (before %.3e) gid %d k %d j %d "
              "i %d E %.6e E_old %.6e T %.6e S %.3e chi_T %.3e\n",
              pmy_pack->pmesh->ncycle,
              !std::isfinite(dnl) ? "BREAKDOWN" : (dnl0 < 0.0 ? "?" : (dnl > 2.0*dnl0 ?
              "DIVERGE" : (dnl > 0.5*dnl0 ? "STALL" : "SLOW"))), dnl, dnl0,
              ok ? pmy_pack->pmb->mb_gid.h_view(m) : -1, k, j, i, h(0), h(1), h(2), h(3),
              h(4));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanAdvectiveCoefficients
//! \brief a_d = v_d + D_db v_b, the advective flux per unit E, v the gas velocity of the
//! old vector plus the iterate's momentum change

void RadiationVet::CleanAdvectiveCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto w = cln_w;
  CleanArgs x = CleanMakeArgs(), xa = x;
  // staggered on the cubed sphere: Cartesian (CleanStagGasMomentum)
  const bool sc = x.cub;
  xa.cub = false;

  CleanFor("cln_advcoef", nmb1, ks, ke, js, je, is, ie,
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
  // their ghosts are first read after the next exchange (the RHS's)
  CleanGhostDefer(cln_w, CLW_A1, 3);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStageOpacity
//! \brief clean_stag_opac_update: the active cells' opacities at the Newton iterate's gas
//! temperature (stage density), then the face coefficients and bracket that use them

void RadiationVet::CleanStageOpacity() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  const int otype = opacity_type;
  const VetOpacTab ot = otab;
  const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, ks_ = kappa_s;
  const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  const bool osc = vscat;
  par_for("cln_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i);
    Real o[4];
    VetCellOpacities(otype, ot, eos, d, w(m,CLW_T,k,j,i), kp, kev, kf, ks_, rref, tref,
                     aa,
                    bb, o[0], o[1], o[2], o[3]);
    op(m,VET_OP_P,k,j,i) = d*o[0];
    op(m,VET_OP_E,k,j,i) = d*o[1];
    op(m,VET_OP_T,k,j,i) = d*(o[2] + o[3]);
    if (osc) op(m,VET_OP_S,k,j,i) = d*o[3];
  });
  CleanGhostFill(opac, 0, opac.extent_int(1), false);
  CleanFaceCoefficients();
  if (cln_consist) CleanStagBracket();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanFirstIterate
//! \brief first Newton iterate: E of the old vector, no momentum change, gas T in local
//! exchange equilibrium with that E.  implicit_predictor: E + h K1_E instead (K1 the
//! slope of U^n), where that changes the radiation temperature by less than
//! implicit_predictor_tcap.  Only the starting point: the solution does not move.

void RadiationVet::CleanFirstIterate() {
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
  // clean_pred_stage: sdirk2's stage 2 from this step's stage-1 slope, E^n + dt K1
  auto k1 = (cln_pst && t2_tab == 2 && t2_solve == VET_T2S_STAGE2) ? t2k2 : t2k1;
  // clean_pred_xch: sdirk2's stage 2 from (K1 + z K2^{n-1})/(1 + z), z = h (ch kE +
  // 4 c kP a T^4/e) the exchange stiffness: stiff cells follow the last step's stage-2
  // slope (the slow manifold's drift), slow ones the latest slope
  const bool pxc = cln_pxc && t2_tab == 2 && t2_solve == VET_T2S_STAGE2;
  auto ks1 = t2k2, ks2 = t2k1;   // sdirk2 stage 2: K1 of this step, K2 of the last
  const Real ch = chat;
  const bool pred = cln_kpred && (t2_ok || cln_k1ok) &&
                    t2_solve != VET_T2S_NONE && !cln_vzero && t2k1.extent_int(0) > nmb1;
  const Real h = dt_sub, ecap = SQR(SQR(impl_pred_tcap)), rtol = cln_trtol;
  const bool xtn = cln_xtnorm;   // clean_xt_norm: E responds, 1/(1 + hc kE) folded
  const bool xtg = cln_xtcons;   // clean_xt_cons: G = 1 - beta^2/2 of the gas
  const Real cl = c_light;

  CleanFor("cln_iter0", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real eo = fmax(w(m,CLW_EO,k,j,i), 0.0);
    Real ep = pred ? eo + h*k1(m,VET_T2_E,k,j,i) : eo;
    if (pred && pxc) {
      const Real ei = feedback ? w(m,CLW_EI,k,j,i) : 0.0;
      const Real tk = (ei > 0.0) ? eos.Temperature(uh(m,IDN,k,j,i), ei) : 0.0;
      const Real z = h*ch*op(m,VET_OP_E,k,j,i) +
                     ((ei > 0.0) ? 4.0*hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(tk))/ei
                                 : 0.0);
      ep = eo + h*(ks1(m,VET_T2_E,k,j,i) + z*ks2(m,VET_T2_E,k,j,i))/(1.0 + z);
    }
    const Real e = (ep > eo/ecap && ep < eo*ecap) ? ep : eo;
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) w(m,CLW_DM+n,k,j,i) = 0.0;
    if (feedback) {
      w(m,CLW_EP,k,j,i) = 0.0;   // beta terms of this T (CleanRightHandSide)
      const Real fe = xtn ? 1.0/(1.0 + hc*op(m,VET_OP_E,k,j,i)) : 1.0;
      const Real gm = xtg ? 1.0 - 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                                       SQR(uh(m,IM3,k,j,i)))/SQR(uh(m,IDN,k,j,i)*cl)
                          : 1.0;
      w(m,CLW_T,k,j,i) = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                             hc*op(m,VET_OP_P,k,j,i)*arad_*fe*gm,
                             w(m,CLW_EI,k,j,i) + hc*op(m,VET_OP_E,k,j,i)*e*fe,
                             w(m,CLW_T,k,j,i), rtol);
    }
  });
  // ghosts of E_k for the gas terms, with the next exchange
  if (feedback) CleanGhostDefer(cln_w, CLW_EK, 1);
}
#endif  // VET_TU_PART(1)

//----------------------------------------------------------------------------------------
//! \fn Real CleanStagAdvFace
//! \brief the advective flux a E at face (k,j,i) normal to n as the face fluxes carry it:
//! clean_stag_upwind takes the share phi as donor-cell a_f E_f plus the lagged rest
//! (CleanFaceWeights), else the central mean; xh the face momentum state

KOKKOS_INLINE_FUNCTION
static Real CleanStagAdvFace(const CleanArgs &x, const CleanArgs &xh, const int n,
                             const int m, const int kk, const int jj, const int ii) {
  const int kl = kk - (n == 2), jl = jj - (n == 1), il = ii - (n == 0);
  const Real ul = x.w(m,CLW_A1+n,kl,jl,il), ur = x.w(m,CLW_A1+n,kk,jj,ii);
  const Real el = x.w(m,CLW_EK,kl,jl,il), er = x.w(m,CLW_EK,kk,jj,ii);
  const Real ac = 0.5*(ul*el + ur*er);
  if (!x.upw || CleanBoundaryFace(x, n, m, kk, jj, ii)) return ac;
  Real al, ar, phi = 0.0;
  CleanFaceWeights(xh, n, m, kk, jj, ii, al, ar, &phi);
  if (!(phi > 0.0)) return ac;
  const Real ah = CleanStagAdvectCorr(x, n, m, kk, jj, ii) +
                  ((ul + ur > 0.0) ? ul*el : ur*er);
  return ac + phi*(ah - ac);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRightHandSide
//! \brief B = E_old + h c kP a T^4 - r h c kE E_k - W - h div F_source + h esrc, the
//! exchange linearised about the iterate (r = s/(d cv + s), s = 4 h c kP a T^3), and the
//! initial guess X = E_k.  W: the work of the iterate's radiation force on the gas (with
//! wb_arad also the reference work h d arad_ref.v); esrc: the explicit energy source.
//! With gas coupling, the beta terms of Ma et al. (2025, eq. 56): emission takes
//! G = 1 - beta^2/2 (CLW_GM) and B adds the lagged hc kE [(beta^2/2 - beta.a/c) E_k +
//! 2 beta.y_k] (y_k = F/c of the iterate; the gas eint takes its opposite, CLW_EP), with
//! beta = v/c (CLW_BT) and hc kP a T^4 (CLW_RE) for the face fluxes (eq. 57).

#if VET_TU_PART(1)
void RadiationVet::CleanRightHandSide(const bool beta, const bool lat) {
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
  // wb_arad: the reference work rides the row (force_reference_work = split)
  const int fref = force_ref;
  const bool wsp = (fref == VET_FREF_WB_ARAD || fref == VET_FREF_WB_MD) && !x.vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  // backward-Euler step, FrefWaccOn: E pays exactly the reference work the gas got in the
  // hydro stages, its share h/dt of fref_wacc
  const bool wpa = wsp && FrefWaccOn();
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : dt_sub/pmy_pack->pmesh->dt;
  // explicit energy source esrc: h esrc joins the old E (the slope of E carries it)
  const bool eso = esrc_on;
  auto es = esrc;
  const bool bud = cln_bon;
  // CleanSteadyStart's equilibrium relax (icr_fmeas): h -> infinity
  const bool eqr = icr_fmeas && x.vzero && x.feedback;

  if (x.feedback) {
    const bool rx = x.cub && cln_vimp;
    CLEAN_HIPX(const CleanArgs *pxb = CleanDevArgs(x, 0);)
    CleanFor<0, 3>("cln_beta", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *pxb;)   // in device memory: no per-thread copy
      Real p[3];
      CleanStagGasMomentum(x, uh, m, k, j, i, p);
      for (int n = 0; n < 3; ++n) {
        x.w(m,CLW_BT+n,k,j,i) = (p[n] + x.w(m,CLW_DM+n,k,j,i))/(uh(m,IDN,k,j,i)*x.c);
      }
      x.w(m,CLW_RE,k,j,i) = x.h*x.c*op(m,VET_OP_P,k,j,i)*arad_*
                            SQR(SQR(x.w(m,CLW_T,k,j,i)));
      // r of the exchange, for the face beta terms (CleanFaceWeights, CleanFaceSource)
      const Real d = uh(m,IDN,k,j,i), t = x.w(m,CLW_T,k,j,i);
      Real b2 = 0.0, e, pg, chi_rho, chi_t, cv;
      for (int n = 0; n < 3; ++n) b2 += SQR(x.w(m,CLW_BT+n,k,j,i));
      eos.ThermoAt(d, t, e, pg, chi_rho, chi_t, cv);
      const Real s = 4.0*x.h*x.c*op(m,VET_OP_P,k,j,i)*(1.0 - 0.5*b2)*arad_*t*t*t;
      x.w(m,CLW_GX,k,j,i) = (s > 0.0) ? s/(d*cv + s) : 0.0;
      if (rx) x.w(m,CLW_JV+2,k,j,i) = d;   // CleanStagDrag's rho: a seam ghost's own
    });
    // E_k's ghosts are those of its last change; BT, RE, GX in one exchange before the
    // lateral sums
    CleanGhostDefer(cln_w, CLW_BT, 4);
    CleanGhostFill(cln_w, rx ? CLW_JV+2 : CLW_GX, rx ? 2 : 1);
    if (lat) CleanStagCsLateral(CLW_EK);
  }
  if (beta) return;
  CleanArgs xh = x;   // face momentum state (clean_stag_split)
  xh.hi = cln_split;
  // a pass without the F solve's F: its F^H, each face once in cln_rhsaf into CLF_N
  const bool fh = x.feedback && !cln_fsk;
  xh.fs = true;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const auto *peos = CleanDev(eos, 2);
  CleanArgs xr = x;   // clean_stag_xdonor: the face sources as they are
  xr.xq = -1;
  const CleanArgs *pxr = CleanDevArgs(xr, 3);
  // each face's source once, for its two cells (CLF_SR), and with gflux its a E (CLF_PT)
  const int nf = x.nf;
  const bool af = x.feedback && x.gf && !x.cub, xsp = (x.xq >= 0) && !cln_xdh;
  // one launch per direction, its face code for a compile-time direction (CleanDirCase)
  for (int dir = 0; dir < nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    CleanFor<4>("cln_rhsaf", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
      CleanDirCase(dir, [&](auto d) {
        constexpr int n = decltype(d)::value;
        if (xsp) CleanStagXSplit(x, *pxr, n, m, k, j, i);
        x.fc[n](m,CLF_SR,k,j,i) = CleanFaceSource(x, n, m, k, j, i);
        if (fh) x.fc[n](m,CLF_N,k,j,i) = CleanFaceFlux<n>(xh, CLW_EK, n, m, k, j, i,
                                                          true);
        if (af) x.fc[n](m,CLF_PT,k,j,i) = CleanStagAdvFace(x, xh, n, m, k, j, i);
      });
    });
  }
  const bool pdo = x.pd;   // clean_stag_central: the rows' diagonal CLW_PD
  const bool xte = (cln_xte > 0) && x.feedback;   // clean_xt_exact: CLW_XA, CLW_XW
  const bool gpic = cln_gpic;   // clean_gas_newton_switch
  // clean_gas_secant: CLW_RR holds the last pass's T
  const bool gsc = cln_gsec == 1 || cln_gsec == 3 ||
                   ((cln_gsec == 2 || cln_gsec == 4) && t2_solve == VET_T2S_STAGE2);
  const bool gmx = cln_gsec >= 3;   // max: the chord only where it is the steeper
  const Real eps_ = std::numeric_limits<Real>::epsilon();
  CleanParForLb<3>("cln_rhs", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
    const auto &eos = *peos;
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hcke = x.h*x.c*op(m,VET_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0, gm = 1.0, sb = 0.0, dop = 0.0, ek = 0.0, gres = 0.0;
    Real eg = 0.0, wq = 0.0;
    bool gr = false;
    if (x.feedback) {
      Real b2 = 0.0, ba = 0.0, by = 0.0, bp = 0.0, y2 = 0.0, fs[3] = {0.0, 0.0, 0.0};
      // clean_stag_fwork: the F solve's F, which the gas takes
      auto ff = [&](const int d, const int kk, const int jj, const int ii) {
        return xh.fc[d](m,CLF_N,kk,jj,ii);
      };
      // gflux: an open boundary face's flux, absorbed in an opaque half cell, enters the
      // cell mean as its half-cell mean, (1 - e^-q)/q, q = chi dx/2 (CleanMomentumChange)
      Real fw[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      for (int n = 0; n < 3; ++n) {   // trip 3: fw stays in registers
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        for (int sd = 0; sd < 2 && x.gf && n < x.nf; ++sd) {
          if (!CleanBoundaryFace(x, n, m, k+sd*dk, j+sd*dj, i+sd*di) ||
              !CleanOpen(x, m, 2*n + sd)) continue;
          const Real q = 0.5*op(m,VET_OP_T,k,j,i)/CleanDivW(x, n, m, k, j, i, sd);
          fw[2*n+sd] = (q > 1.0e-8) ? -expm1(-q)/q : 1.0;
        }
      }
      for (int n = 0; n < 3; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        if (n < x.nf) fs[n] = fw[2*n]*ff(n, k, j, i) + fw[2*n+1]*ff(n, k+dk, j+dj, i+di);
      }
      if (x.cub) {   // the cell's flux vector (CleanStagCellVector)
        Real c[6];
        for (int d = 0; d < 3; ++d) {
          c[2*d] = fw[2*d]*ff(d, k, j, i);
          c[2*d+1] = fw[2*d+1]*ff(d, k+(d == 2), j+(d == 1), i+(d == 0));
        }
        CleanStagCellVector(x, m, k, j, i, c, fs);
        for (int d = 0; d < 3; ++d) fs[d] *= 2.0;
      }
      for (int n = 0; n < 3; ++n) {
        const Real bt = x.w(m,CLW_BT+n,k,j,i);
        if (n < x.nf) by += bt*0.5*fs[n]/x.c;
        y2 += SQR(0.5*fs[n]/x.c);
        // gflux: a E face by face (cubed sphere: the cell's a E, its face families are
        // not the vector's components)
        if (n < x.nf && x.gf && !x.cub) {
          const int dk = (n == 2), dj = (n == 1), di = (n == 0);
          bp += bt*0.5*(x.fc[n](m,CLF_PT,k,j,i) + x.fc[n](m,CLF_PT,k+dk,j+dj,i+di))/x.c;
        }
        b2 += SQR(bt);
        ba += bt*x.w(m,CLW_A1+n,k,j,i)/x.c;
      }
      if (x.gf && x.cub) bp = ba*x.w(m,CLW_EK,k,j,i);
      // y realisable, |y| <= E; the gas eint takes the same sb (CLW_EP): energy exact
      ek = fmax(x.w(m,CLW_EK,k,j,i), 0.0);
      if (y2 > SQR(ek)) by *= ek/std::sqrt(y2);
      gm = 1.0 - 0.5*b2;
      sb = hcke*((0.5*b2 - ba)*x.w(m,CLW_EK,k,j,i) + 2.0*by);
      // gflux: sb = hcke [dop E + 2 beta.F'/c] with the comoving flux F' = F - a E taken
      // face by face, not as the faces' mean F less the cell's a E (at a front these
      // differ by the neighbour's a E, a lagged source the iteration could not settle)
      dop = 0.5*b2 + ba;
      if (x.gf) sb = hcke*(dop*x.w(m,CLW_EK,k,j,i) + 2.0*(by - bp));
      // clean_stag_freal: the comoving part realisable, |E0 - E| <= 2 |beta| E (a face F'
      // of a hot neighbour through a cold opaque cell is no flux of its E); gas: same sb
      if (x.fre) {
        const Real s0 = hcke*dop*x.w(m,CLW_EK,k,j,i), sl = 2.0*hcke*std::sqrt(b2)*ek;
        sb = s0 + fmin(fmax(sb - s0, -sl), sl);
      }
      if (x.gf) x.w(m,CLW_JV+1,k,j,i) = hcke*dop;   // d sb/dE (CleanNewtonUpdate)
      gr = sb != -x.w(m,CLW_EP,k,j,i);   // T took the last pass's sb
    }
    const Real hckp = x.h*x.c*op(m,VET_OP_P,k,j,i)*gm;
    x.w(m,CLW_GM,k,j,i) = gm;
    x.w(m,CLW_EP,k,j,i) = -sb;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, eg, p, chi_rho, chi_t, cv);
      // exchange slopes: the tangent at the iterate's T (Newton, quadratic) where the
      // exchange in B, EO + (1 - r) hc kP a T^4 - r (e(T) - eint), is >= 0; else the
      // secant through T_o (the old vector's gas T) and T, exact at T, so the solution is
      // the same and the exchange in B is EO + (1 - r) hc kP a T_o^4 >= 0 by construction
      // (the tangent overshoots the root of the concave E -> e map when radiation heats
      // cold gas, giving E < 0).
      Real s = 4.0*hckp*arad_*t*t*t, cg = d*cv;
      // clean_stag_eicarry: a carried deficit eint < 0 takes the chord's T_o at 0
      const Real ei = fmax(x.w(m,CLW_EI,k,j,i), 0.0), eo = x.w(m,CLW_EO,k,j,i);
      const bool sec = eg > ei && eo*(cg + s) + cg*hckp*arad_*SQR(SQR(t)) < s*(eg - ei);
      const Real to = sec ? eos.Temperature(d, ei, t) : t;
      if (sec && t - to > 1.0e-6*t) {
        cg = (eg - ei)/(t - to);
        s = hckp*arad_*(t + to)*(t*t + to*to);
      }
      // clean_gas_secant: elsewhere the chord through the last pass's T (the tangent at
      // a c_v jump, e.g. H ionisation, cycles across the jump); the root is unchanged
      const Real tp = x.w(m,CLW_RR,k,j,i);
      // max: the chord only where it exceeds the tangent d cv by more than the round-off
      // of e(T) - e(T_p) (Newton where e(T) is smooth; bitwise tangent for an ideal gas)
      if (gsc && !sec && tp > 0.0 && tp != t) {
        Real ep, pp, crp, ctp, cvp;
        eos.ThermoAt(d, tp, ep, pp, crp, ctp, cvp);
        const Real dx = (eg - ep - cg*(t - tp))*(t - tp);
        if (!gmx ||
            dx > 16.0*eps_*(fabs(eg) + fabs(ep) + fabs(cg*(t - tp)))*fabs(t - tp)) {
          cg = (eg - ep)/(t - tp);
          s = hckp*arad_*(t + tp)*(t*t + tp*tp);
        }
      }
      rr = (s > 0.0) ? s/(cg + s) : 0.0;   // no exchange where T = 0
      // 1 - r = cg/(cg + s) without cancellation
      x.w(m,CLW_RY,k,j,i) = (s > 0.0) ? cg/(cg + s) : 1.0;
      // clean_gas_newton_switch: a stalled solve's remaining passes drop the gas Newton
      // slope (Picard on the bracketed root T of the gas row)
      if (gpic) {
        rr = 0.0;
        x.w(m,CLW_RY,k,j,i) = 1.0;
      }
      // the gas row's residual at the iterate with these sb, share r: a change of sb
      // between passes enters E with 1 - r (taken whole it contracts only ~0.5 per pass
      // in radiation-dominated opaque gas)
      if (gr) gres = rr*(x.w(m,CLW_EI,k,j,i) + hcke*ek - sb - eg -
                         hckp*arad_*SQR(SQR(t)));
      // gflux: the Doppler part hcke dop E of sb is taken in the Newton step as the
      // exchange (r -> r + dop (1 - r), sb kept at E_k), not lagged: where hcke dop >> 1
      // the lagged term diverges
      if (x.gf) {
        const Real dc = fmin(fmax(dop, 0.0), 1.0);
        rr += dc*(1.0 - rr);
        x.w(m,CLW_RY,k,j,i) *= 1.0 - dc;
      }
      // work of the radiation force on the gas, (v_old + dm/(2 rho)).dm: exactly the
      // kinetic energy the gas gains; with wb_arad also the reference force's work
      // h rho arad_ref.v at the iterate's velocity
      Real pg[3], rf[3] = {0.0, 0.0, 0.0};
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      if (wsp) CleanStagRefForce(x, fref, fr, aref, d, m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (pg[n] + 0.5*dm)*dm/d;
        if (wsp && !wpa) work += x.h*rf[n]*(pg[n] + dm)/d;
      }
      if (wpa) work += fra*wc(m,k,j,i);
      // clean_stag_central: a drain W > 0 goes on the diagonal as q X, q = W/E_k (same
      // fixed point, b >= 0, X > 0; the write-back books q X); a source stays in B
      if (x.ucp && work > 0.0 && ek > 0.0) {
        wq = work/ek;
        work = 0.0;
      }
    }
    x.w(m,CLW_RR,k,j,i) = gsc ? t : rr;
    // eqr: the exchange via the gas row of the root T, hckp a T^4 - r hcke E =
    // EI - e(T) - sb + (1 - r) hcke E, avoiding the cancellation of two ~ hcke E terms
    // (round-off at h -> infinity exceeds E in opaque cells)
    const Real sdv = CleanSourceDivergence(x, m, k, j, i);
    const Real hes = eso ? x.h*es(m,k,j,i) : 0.0;
    x.w(m,CLW_B,k,j,i) = eqr ?
        x.w(m,CLW_EO,k,j,i) + x.w(m,CLW_EI,k,j,i) - eg - work
        + x.w(m,CLW_RY,k,j,i)*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - sdv + hes :
        x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
        - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work - sdv + sb + gres + hes;
    if (pdo) x.w(m,CLW_PD,k,j,i) = wq;
    // clean_xt_exact: the row's exchange as linearised (B less its rest) and its work
    if (xte) {
      x.w(m,CLW_XA,k,j,i) = eqr ? 0.0 : x.w(m,CLW_B,k,j,i) - (x.w(m,CLW_EO,k,j,i) - work
          - sdv + hes);
      x.w(m,CLW_XW,k,j,i) = work;
    }
    x.w(m,CLW_X,k,j,i) = x.w(m,CLW_EK,k,j,i);
    if (!bud) return;   // clean_conv = budget: the gross terms of the E and gas rows
    const Real xc = hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) + hckp*arad_*SQR(SQR(t)) +
                    fabs(work) + fabs(sb);
    // face flux scale: the old vector's, at least what the cell's E can stream through
    // it in h, c theta E, theta = 1/(1 + h c chi); a physical boundary face takes the
    // cell's streaming alone
    Real fl = eso ? fabs(es(m,k,j,i)) : 0.0;
    const Real eo = fabs(x.w(m,CLW_EO,k,j,i));
    const Real cs = x.c*eo/(1.0 + x.h*x.c*op(m,VET_OP_T,k,j,i));
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
  // an empty cell (no energy, no flux): its rows at 1e-12 of the largest scale
  auto w = cln_w;
  Real smax[2] = {0.0, 0.0};
  // flat range reduction
  const int n1 = ie - is + 1, n2 = je - js + 1, nji = n2*n1, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("cln_smax", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &se, Real &sg) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    se = fmax(se, w(m,CLW_SE,k,j,i));
    sg = fmax(sg, w(m,CLW_SG,k,j,i));
  }, Kokkos::Max<Real>(smax[0]), Kokkos::Max<Real>(smax[1]));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, smax, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  const Real fse = 1.0e-12*smax[0];
  const Real fsg = fmax(1.0e-12*smax[1], 1.0e-300);
  CleanFor("cln_sfloor", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_SE,k,j,i) = fmax(w(m,CLW_SE,k,j,i), fmax(fse, 1.0e-300));
    w(m,CLW_SG,k,j,i) = fmax(w(m,CLW_SG,k,j,i), fsg);
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanNewtonUpdate
//! \brief new iterate E_k = X and the gas temperature in exchange equilibrium with it
//! (absorption of max(X, 0); also the iterate without clean_stag_gflux).  Returns
//! max |dE|/E, |dT|/T (clean_nl_norm = cell: |dE| over max(|E|, cln_nl_escale)).

Real RadiationVet::CleanNewtonUpdate() {
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
  const bool poslim = cln_gflux;
  const bool gfd = cln_gflux && cln_feedback;
  // clean_conv = budget: changes over the rows' gross terms, and the time error
  // h/2 |f(Y) - f(U^n)| of the iterate (Heun against Euler, K1 = f(U^n))
  const bool bud = cln_bon;
  const bool kok = bud && (t2_ok || cln_k1ok) && t2_solve != VET_T2S_NONE &&
                   !cln_vzero && t2k1.extent_int(0) >= nmb;
  const Real h = dt_sub, rtol = cln_trtol;
  auto k1 = t2k1;

  // clean_xt_exact x2: T from the local 2x2 of the E and gas rows (S = X diag - XA fixed)
  const bool x2 = (cln_xte == 2) && cln_feedback;
  Real change = 0.0, femax = 0.0, ermax = 0.0;
  // CUDA: 4 blocks of 256 per SM; the reduction keeps its blocks and grid (bitwise)
  using NewtonRp =
      Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, CleanLbMb(4)>>;
  AthReduce("cln_newton", NewtonRp(DevExeSpace(), 0, nmb*nkji),
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
    Real x = w(m,CLW_X,k,j,i);
    const Real xr = poslim ? x : fmax(x, 0.0);   // the change: undamped
    const Real ex = fmax(x, 0.0);
    const Real e = poslim ? x : ex, g = w(m,CLW_GM,k,j,i);
    // gflux: the gas takes the Doppler part hc kE dop E of sb at X, as the E row in r
    // (lagged at E_k, T and E drifted apart pass by pass in moving opaque cells)
    if (gfd) w(m,CLW_EP,k,j,i) -= w(m,CLW_JV+1,k,j,i)*(x - e_old);
    const Real hcke = hc*op(m,VET_OP_E,k,j,i);
    Real t;
    if (x2) {
      const Real dg = 1.0 + hcke*w(m,CLW_RY,k,j,i);
      const Real sr = x*dg - w(m,CLW_XA,k,j,i);
      t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
              hc*op(m,VET_OP_P,k,j,i)*g*arad_/(1.0 + hcke),
              w(m,CLW_EI,k,j,i) + (hcke*sr + w(m,CLW_EP,k,j,i))/(1.0 + hcke), t_old,
              rtol);
    } else {
      t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i), hc*op(m,VET_OP_P,k,j,i)*g*arad_,
              w(m,CLW_EI,k,j,i) + hcke*ex + w(m,CLW_EP,k,j,i), t_old, rtol);
    }
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_T,k,j,i) = t;
    const Real de = ncell ? fabs(xr - e_old)/fmax(fmax(fabs(xr), fabs(e_old)), esc)
                          : fabs(xr - e_old)/fmax(fmax(fabs(xr), fabs(e_old)), 1.0e-300);
    const Real dt = fabs(t - t_old)/fmax(fmax(t, t_old), 1.0e-300);
    if (bud) {
      const Real d = uh(m,IDN,k,j,i), se = w(m,CLW_SE,k,j,i);
      Real e1, e0, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e1, p, chi_rho, chi_t, cv);
      eos.ThermoAt(d, t_old, e0, p, chi_rho, chi_t, cv);
      dmax = fmax(dmax, fmax(fabs(xr - e_old)/se, fabs(e1 - e0)/w(m,CLW_SG,k,j,i)));
      if (kok) er_max = fmax(er_max, 0.5*fabs(e - w(m,CLW_EO,k,j,i) -
                                              h*k1(m,VET_T2_E,k,j,i))/se);
    } else {
      dmax = fmax(dmax, fmax(de, dt));
    }
    // round-off level of E (its terms over E)
    fe_max = fmax(fe_max, (fabs(e)*(1.0 + hc*op(m,VET_OP_E,k,j,i)) +
        fabs(w(m,CLW_EO,k,j,i)) + hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(t)))/
        fmax(fabs(e), 1.0e-300));
  }, Kokkos::Max<Real>(change), Kokkos::Max<Real>(femax), Kokkos::Max<Real>(ermax));
  // one reduction over the ranks
  Real mx[3] = {change, femax, ermax};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, mx, 3, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  change = mx[0]; femax = mx[1]; ermax = mx[2];
  cln_fe = std::numeric_limits<Real>::epsilon()*femax;
  if (bud) {   // (the tolerance of the next test and solve; for settled, the level of a
    // change the linear solve fixes at the floor, its 0.1 tfl in the budget's norm)
    // (a step dt_retry redoes, whose state just failed, is solved strictly too)
    // (clean_clos1_budget, clean_bud_maxnorm, clean_vel_settled: the time error down to
    // clean_tol_nl, as each relaxes the passes to it; at 10 clean_tol a step whose time
    // error is below, a 1e-6 wave's 1e-11 at a change of 1e-9, lost the wave's phase)
    const Real tfl = fmax(cln_tol_nl, 10.0*cln_tol);
    const Real tlo = (cln_c1b || cln_bmn || cln_vset) ? cln_tol_nl : tfl;
    cln_tolb = (kok && !dtr_redo) ?
               fmin(fmax(cln_bkap*ermax, tlo), cln_bmax) : tfl;
    cln_fe = 0.1*tfl;
  }
  // fwork: emission hc kP a T^4 of the new T for the face beta terms of the velocity
  // update's F solve
  if (cln_fwork && cln_feedback) {
    par_for("cln_sbre", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is,
            is+nx1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_RE,k,j,i) = hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(w(m,CLW_T,k,j,i)));
    });
    CleanGhostDefer(cln_w, CLW_RE, 1);
  }
  // staggered: the ghosts of the new E_k, for the velocity update and the next RHS
  if (cln_feedback) CleanGhostFill(cln_w, CLW_EK, 1);
  return change;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanVelocityUpdate
//! \brief momentum change of the new iterate less force_reference's reference force, and
//! the new velocity's advective coefficients.  Returns max |dv| over max |v|.

Real RadiationVet::CleanVelocityUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanStagCsLateral(CLW_EK);   // E_k's ghosts: CleanNewtonUpdate
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;
  xh.mcl = true;
  if (cln_mclc.extent_int(0) != 1) Kokkos::realloc(cln_mclc, 1);
  Kokkos::deep_copy(cln_mclc, 0.0);
  xh.mclc = cln_mclc;
  // the F solve's face tie after the drag (same face T: CLF_PT)
  const bool td = cln_fwork;
  // clean_stag_fwork: the F solve's momentum, which the write-back gives the gas (work W)
  if (cln_fwork) {
    CleanStagFSolve(CLW_EK, !td, false);
    xh.fs = true;
    xh.tf = !x.cub && !cln_tgid;
    cln_fsk = true;
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto op = opac;
  const Real arad_ = arad;
  const int fref = force_ref;
  const bool cap = (fref == VET_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;

  // with the F solve's F: the faces' terms and drag first, once per face, then the cells'
  // momentum and advective coefficients
  const bool fa = xh.fs;
  if (fa) CleanStagDrag(false, true, false);   // lateral sums of E_k: above
  if (td) CleanStagFaceTie(CLW_EK, (t2_solve == VET_T2S_NONE) ? 0 : ((t2_solve ==
                           VET_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD), true);
  xh.fa = fa;
  auto xa = x;
  xa.cub = false;
  Real dvmax = 0.0, vmax = 0.0, fvmax = 0.0;
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
             const CleanArgs *pxa = CleanDevArgs(xa, 4);)
  // fa: CleanStagMomentum<true> reads the faces' terms, in a kernel of its own (HIP: 256
  // VGPRs, 2 waves per SIMD)
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int vmb = 2;
#else
  constexpr int vmb = CleanLbMb(2);
#endif
  auto cell = KOKKOS_LAMBDA(const int idx, int &m, int &k, int &j, int &i) {
    m = idx/nkji;
    k = (idx - m*nkji)/nji + ks;
    j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
  };
  if (fa) {
    AthReduce("cln_velocity", Kokkos::RangePolicy<DevExeSpace,
              Kokkos::LaunchBounds<256, vmb>>(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // device memory
      int m, k, j, i;
      cell(idx, m, k, j, i);
      CleanVelocityCell<true>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, true, m,
                              k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  } else {
    AthReduce("cln_velocity", CleanRpLb2(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // device memory
      int m, k, j, i;
      cell(idx, m, k, j, i);
      CleanVelocityCell<false>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, false, m,
                               k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  }
  CleanReduceMax(dvmax, vmax, fvmax);
  if (fa) {
    CleanGhostDefer(cln_w, CLW_A1, 3);   // as CleanAdvectiveCoefficients
  } else {
    if (cln_vimp) CleanStagDrag(false, false, false);   // this dm's coefficients
    CleanAdvectiveCoefficients();
  }
  vmax = fmax(vmax, 1.0e-300);
  cln_fv = std::numeric_limits<Real>::epsilon()*fvmax/vmax;
  return dvmax/vmax;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanCopyIterate
//! \brief the solution is the last Newton iterate

void RadiationVet::CleanCopyIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  CleanFor("cln_copy", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_X,k,j,i) = w(m,CLW_EK,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CleanCountNegative
//! \brief the number of cells with E < 0 beyond the solver tolerance (10 tol max|B|);
//! CleanWriteBack raises those within it to zero, counted

int RadiationVet::CleanCountNegative() {
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
//! \fn int RadiationVet::CleanStagRealizable
//! \brief staggered: realizability of the solution X beyond the solver tolerance (10
//! clean_tol max|B|).  R1: cells with E < 0 (cln_nr1); R2: interior faces whose flux
//! exceeds c E of the donor (cln_nr2).  Returns cln_nr1 + cln_nr2.

int RadiationVet::CleanStagRealizable() {
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
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real r1 = 0.0, r2 = 0.0;
    Kokkos::parallel_reduce("cln_real", CleanRpLb2(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &q1, Real &q2) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // in device memory: no per-thread copy
      int m = idx/fkji;
      int k = (idx - m*fkji)/fji;
      int j = (idx - m*fkji - k*fji)/f1;
      int i = idx - m*fkji - k*fji - j*f1;
      const bool cell = (i < nx1 && j < nx2 && k < nx3);   // its low face
      k += ks;
      j += js;
      i += is;
      const Real el = x.w(m,CLW_X,k-dk,j-dj,i-di), er = x.w(m,CLW_X,k,j,i);
      if (dir == 0 && cell && !(er >= -zero)) q1 += 1.0;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const Real f = CleanFaceFlux(x, CLW_X, dir, m, k, j, i, true);
      x.fc[dir](m,CLF_PT,k,j,i) = f;   // for CleanWriteBack, next (QL is free)
      const Real ed = (f > 0.0) ? el : er;
      const bool r2 = !(fabs(f) <= (1.0 + 1.0e-12)*x.c*ed + x.c*zero);
      if (r2) q2 += 1.0;
    }, Kokkos::Sum<Real>(r1), Kokkos::Sum<Real>(r2));
    n1 += r1;
    n2 += r2;
  }
  Real nr[2] = {n1, n2};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, nr, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  cln_nr1 = static_cast<int>(nr[0]);
  cln_nr2 = static_cast<int>(nr[1]);
  return cln_nr1 + cln_nr2;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagBracket
//! \brief staggered clean_consistency: the bracket C_f = w (F_ray,f - Phi_f(J)) of every
//! interior face, F_ray the rays' energy flux (CLN_FX, CLN_F2, CLN_F3 of the upper cell),
//! Phi_f(J) the flux of the rays' own J, so Phi(J) + C = (1 - w) Phi(J) + w F_ray.
//! Lagged (fixed in the solve), energy flux only.  It enters the streaming share as flux
//! factors, Phi_w = c (g+ E_L + g- E_R) with
//!   c (g+ J_L + g- J_R) w = c (h+ J_L + h- J_R) w + C_f  (CLF_HP, CLF_HM):
//! the half ranges H+-_f = h+-_f J tilted by dd = C_f/(w c) in proportion to their sizes,
//!   H+-' = H+-_f (1 +- dd/s),  s = H+_f - H-_f,
//! the shrinking one at most to 0.  Linear in E with upwind signs, so the M-matrix holds
//! for any E.  g+- are not capped at 1 (a beam edge's face flux can exceed c J_donor).
//! Cubed sphere: a panel seam's high x2, x3 faces take the rays' flux from their last cell
//! (CLN_FH, CleanRayFace).  clean_stag_thinf takes the flux of the cells' half ranges there (the face
//! sum with fixed sets): the two MeshBlocks of a face mix or extrapolate its sum with
//! their own cells' weights, and the write-back's face tie would move the difference
//! into one side's E.  clean_stag_rayratio: J -> E in Phi and the half ranges, F_ray ->
//! (F_ray/J_d) E_d of the upwind donor d (no bracket where J_d <= e_floor).

void RadiationVet::CleanStagBracket() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.feedback = false;   // the gas terms read the last iterate's work vectors
  const Real efl = e_floor;
  const bool rat = cln_rayratio, thf = cln_thinf;
  const bool hct = cln_hctau;
  const bool lat = thf && pmy_pack->pmesh->use_cubed_sphere;
  const Real rmn = cln_raymin;
  auto iw_ = iw;
  // clean_stag_thinf: R, what the rays hold of the cell: its J and its faces' fluxes
  CleanFor("cln_rayj", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x.w(m,CLW_Z,k,j,i) = rat ? iw_(m,VET_IW_EN,k,j,i) : x.cs(m,CLN_J,k,j,i);
    if (!thf) return;
    Real r = fmax(x.cs(m,CLN_JR,k,j,i), 0.0);
    for (int n = 0; n < x.nf; ++n) {
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*(n == 2), jf = j + sd*(n == 1), iff = i + sd*(n == 0);
        if (sd > 0 && CleanBoundaryFace(x, n, m, kf, jf, iff)) continue;
        r = fmax(r, fabs(CleanRayFace(x, n, m, kf, jf, iff)));
      }
    }
    x.w(m,CLW_Y,k,j,i) = r;
  });
  CleanGhostFill(cln_w, thf ? CLW_Y : CLW_Z, thf ? 2 : 1);
  CleanStagCsLateral(CLW_Z);
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
  const bool hrm = cln_hface;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    CleanFor("cln_sbracket", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // in device memory: no per-thread copy
      // Phi(J) at its steady value, -(c/chi_f) div(D J), in the share 1 - w, so the
      // bracket (and the scheme's steady state) does not depend on h or G.  Where the
      // diffusion flux exceeds free streaming, lambda = |div(D J)|/(chi_f J_f) > 1, it
      // takes c J_f along -div(D J) (Larsen's max limiter)
      const Real w = x.fc[dir](m,CLF_W,k,j,i);
      const Real chi = CleanFaceChi(x.op, x.x1f, x.mbc, hrm, dir, m, k, j, i, is, ie);
      const Real st = CleanFaceStress(x, CLW_Z, dir, m, k, j, i, true);
      const Real jf = 0.5*(x.w(m,CLW_Z,k,j,i) + x.w(m,CLW_Z,k-dk,j-dj,i-di));
      const Real cd = (jf > 0.0 && fabs(st) > chi*jf) ? x.c*jf/fabs(st) :
                      ((chi > 0.0) ? x.c/chi : x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c);
      // cells' half ranges realisable, h+ in [0, 1], h- in [-1, 0]; the rest is DG below
      x.fc[dir](m,CLF_HP,k,j,i) = fmin(fmax(x.fc[dir](m,CLF_HP,k,j,i), 0.0), 1.0);
      x.fc[dir](m,CLF_HM,k,j,i) = fmax(fmin(x.fc[dir](m,CLF_HM,k,j,i), 0.0), -1.0);
      Real gl = 0.0, gr = 0.0;
      CleanStagArea(x, dir, m, k, j, i, w, gl, gr);
      const Real pj = -(1.0 - w)*cd*st +
          w*x.c*(x.fc[dir](m,CLF_HP,k,j,i)*x.w(m,CLW_Z,k-dk,j-dj,i-di) +
                 x.fc[dir](m,CLF_HM,k,j,i)*x.w(m,CLW_Z,k,j,i)) +
          gl*x.w(m,CLW_Z,k-dk,j-dj,i-di) + gr*x.w(m,CLW_Z,k,j,i);
      const Real fr0 = CleanRayFace(x, dir, m, k, j, i);
      const bool bf = CleanBoundaryFace(x, dir, m, k, j, i);
      if (bf) return;
      // the rays' J also below clean_ray_min E (a diagonal beam's side cells);
      // rayratio: E, and the rays' flux per J of its donor times the donor's E
      const Real jl = rat ? x.w(m,CLW_Z,k-dk,j-dj,i-di) : x.cs(m,CLN_JR,k-dk,j-dj,i-di);
      const Real jr = rat ? x.w(m,CLW_Z,k,j,i) : x.cs(m,CLN_JR,k,j,i);
      const Real jd = (fr0 > 0.0) ? x.cs(m,CLN_JR,k-dk,j-dj,i-di) : x.cs(m,CLN_JR,k,j,i);
      if (rat && !(jd > efl)) {
        x.fc[dir](m,x.fdg,k,j,i) = 0.0;
        return;
      }
      // clean_face_hctau: the rays' face flux interpolated linearly in optical depth,
      // H_f = (tau_R H_L + tau_L H_R)/(tau_L + tau_R), tau the half-cell depths (H r^2 on
      // sp / cs radial faces; x1 faces there, every face on Cartesian).  F_ray (x
      // (r_d/r_f)^2 on sp / cs) keeps the share e = (1 - min/max of the cells' r^2 H)
      // e^-tau_d, so a jump in H (beam edge, shadow) keeps the transport solution and
      // smooth or opaque regions the interpolation; a convex blend, single-valued per
      // face.
      Real fy = rat ? fr0/jd*((fr0 > 0.0) ? jl : jr) : fr0;
      if (lat && dir > 0) {
        fy = x.fc[dir](m,CLF_HP,k,j,i)*jl + x.fc[dir](m,CLF_HM,k,j,i)*jr;
      }
      if (hct && (dir == 0 || !(x.sph || x.cub))) {
        const int il = i - di;
        auto hc = [&](const int kk, const int jj, const int ii) {
          return (x.cs(m,CLN_HP+dir,kk,jj,ii) + x.cs(m,CLN_HM+dir,kk,jj,ii))*
                 x.cs(m,CLN_JR,kk,jj,ii);
        };
        Real wl, wr, tg[4];
        CleanHcTauW(x, dir, m, k, j, i, wl, wr, tg);
        const Real hl = hc(k-dk, j-dj, il), hr = hc(k, j, i);
        const Real al = fabs(tg[2]*hl), ar = fabs(tg[3]*hr);
        const Real e = ((hl*hr > 0.0) ? 1.0 - fmin(al, ar)/fmax(al, ar) : 1.0)*
                       exp(-tg[(fr0 > 0.0) ? 0 : 1]);
        const Real ga = (tg[(fr0 > 0.0) ? 2 : 3] - 1.0)*fy;
        fy = (1.0 - e)*(wl*hl + wr*hr) + e*(fy + ga);
      }
      // clean_stag_thinf: the rest of the rays' face flux, dd = F_y - h+ J_L - h- J_R,
      // tilted into the face half ranges in proportion to their sizes (the shrinking one
      // at most to 0), per J of the donor (at least min(clean_ray_min E_old, R)), g+ in
      // [0, 1], g- in [-1, 0]: the thin share w c (g+ E_L + g- E_R) is within [-c E_R,
      // c E_L], realisable by construction, no DG.  A flux the cap cuts is carried by a
      // higher E (a diagonal beam's side cells: J ~ 0, R ~ E).  In a ray hole (R <<
      // clean_ray_min E_old) E streams along the rays' own direction: per clean_ray_min
      // E_old, g -> 0 as J -> 0, so E would stay, its flux independent of it.
      if (thf) {
        const Real hl = x.fc[dir](m,CLF_HP,k,j,i)*jl, hr = x.fc[dir](m,CLF_HM,k,j,i)*jr;
        const Real dd = fy - (hl + hr), fr = hl + hr + dd, s = hl - hr;
        const Real tl = hl + dd*hl/s, tr = hr - dd*hr/s;
        const bool cut = (tl < 0.0) || (tr > 0.0);
        const Real hp = !(s > 0.0) ? hl : (cut ? ((tl < 0.0) ? 0.0 : fmax(fr, 0.0)) : tl);
        const Real hm = !(s > 0.0) ? hr : (cut ? ((tl < 0.0) ? fmin(fr, 0.0) : 0.0) : tr);
        const Real el = fmax(x.w(m,CLW_EO,k-dk,j-dj,i-di), 0.0);
        const Real er = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
        const Real dl = fmax(jl, fmin(rmn*el, x.w(m,CLW_Y,k-dk,j-dj,i-di)));
        const Real dr = fmax(jr, fmin(rmn*er, x.w(m,CLW_Y,k,j,i)));
        const Real gp = fmin(hp/dl, 1.0), gm = fmax(hm/dr, -1.0);
        if (jl > 0.0) x.fc[dir](m,CLF_HP,k,j,i) = gp;
        if (jr > 0.0) x.fc[dir](m,CLF_HM,k,j,i) = gm;
        x.fc[dir](m,x.fdg,k,j,i) = 0.0;
        return;
      }
      // default form: half ranges h+- realisable, the rest the lagged flux DG/w = q c
      // (F_y - h+ J_L - h- J_R) (CleanThinFaceDG), F_y clamped to [-J_R, J_L], q = r^2
      // the quasi-steadiness gate, r the worse of the two cells' min(E_old/J, J/E_old)
      auto rq = [&](const Real a, const Real b) {
        return (fmax(a, b) > 0.0) ? fmin(a, b)/fmax(a, b) : 1.0;
      };
      const Real r = fmin(rq(fmax(x.w(m,CLW_EO,k-dk,j-dj,i-di), 0.0), fmax(jl, 0.0)),
                          rq(fmax(x.w(m,CLW_EO,k,j,i), 0.0), fmax(jr, 0.0)));
      const Real hup = x.fc[dir](m,CLF_HP,k,j,i)*jl + x.fc[dir](m,CLF_HM,k,j,i)*jr;
      x.fc[dir](m,x.fdg,k,j,i) = r*r*x.c*(fmin(fmax(fy, -jr), jl) - hup);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanDMem
//! \brief clean_stag_fkin: the energy flux memory d (CLF_N + 3) at the stage-loop start
//! into CLF_DS (save) or back (a rejected stage's write-backs are discarded)

void RadiationVet::CleanDMem(const bool save) {
  if (!cln_fsolve || !cln_fkin) return;
  if (!cln_init) CleanInit();
  auto AL = Kokkos::ALL;
  for (auto fc : {cln_fc1, cln_fc2, cln_fc3}) {
    auto d = Kokkos::subview(fc, AL, CLF_N + 3, AL, AL, AL);
    auto s = Kokkos::subview(fc, AL, CLF_DS, AL, AL, AL);
    if (save) {Kokkos::deep_copy(s, d);} else {Kokkos::deep_copy(d, s);}
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanDInc
//! \brief clean_stag_fkin: d of a hesdirk2 stage's old vector: the stage start state
//! (stage 2: w1 d^n + w0 Y1) plus the increment c1 K1 + c2 K2 of d's own slopes (CLF_KD,
//! CLF_KD + 1), as every other component of the integrated state

void RadiationVet::CleanDInc(const bool two, const Real c1, const Real c2, const Real w1,
                             const Real w0) {
  if (!cln_fsolve || !cln_fkin) return;
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
#endif  // VET_TU_PART(1)

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanDtrD
//! \brief dt_retry: the memory d and its slope K1 (CLF_N + 3, CLF_KD), the only face
//! values a cycle takes from the last one, to / from cln_rstd (the restart's layout)

#if VET_TU_PART(2)
void RadiationVet::CleanDtrD(const bool save) {
  const int nmb = cln_fc1.extent_int(0);
  if (!cln_fsolve || !cln_fkin || nmb == 0) return;
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
#endif  // VET_TU_PART(2)

//----------------------------------------------------------------------------------------
//! \struct CleanFLine
//! \brief the F solve's lines on every mesh (CleanStagFSolve): the cell terms t of cell
//! (k,j,i) at its face (kf,jf,iff) normal to dir (al, lagged pressure p, w, h+, h-, area
//! pressure pa) and the F of the line's faces.  Spherical polar, Cartesian: t in
//! CLW_L + 3 dir, F of each cell's low face in CLW_Y + dir (CLW_S).  Cubed sphere: cln_v
//! in the owner's axes (an x2 line may continue as an x3 line across a seam, so a cell
//! keeps both faces of each family; family q: terms 3 q + 0..2, 15 + q, F 9 + 2 q + hi)

struct CleanFLine {
  DvceArray5D<Real> v;
  bool cub;
  KOKKOS_INLINE_FUNCTION
  static int Cy(const int dir) {return (dir == 2) ? CLW_S : CLW_Y + dir;}
  KOKKOS_INLINE_FUNCTION
  void Cell(const CleanArgs &x, const int dir, const int m, const int kf, const int jf,
            const int iff, const int k, const int j, const int i, Real t[6]) const {
    if (!cub) {
      for (int n = 0; n < 3; ++n) t[n] = x.w(m,CLW_L+3*dir+n,k,j,i);
      t[3] = x.cs(m,CLN_HP+dir,k,j,i);
      t[4] = x.cs(m,CLN_HM+dir,k,j,i);
      t[5] = 0.0;
      return;
    }
    Real s;
    const int q = CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, false, s)/2;
    for (int n = 0; n < 3; ++n) t[n] = v(m,3*q+n,k,j,i);
    const Real hp = x.cs(m,CLN_HP+q,k,j,i), hm = x.cs(m,CLN_HM+q,k,j,i);
    t[3] = (s > 0.0) ? hp : -hm;
    t[4] = (s > 0.0) ? hm : -hp;
    t[5] = v(m,15+q,k,j,i);
  }
  // slot of the F of face (kf,jf,iff) (far: of the family's other face of cell (k,j,i))
  // in that cell's storage, s its orientation along ours
  KOKKOS_INLINE_FUNCTION
  Real &Slot(const CleanArgs &x, const int dir, const int m, const int kf, const int jf,
             const int iff, const int k, const int j, const int i, const bool far,
             Real &s) const {
    s = 1.0;
    if (cub) return v(m,9+CleanCsSlot(x, dir, m, kf, jf, iff, k, j, i, far, s),k,j,i);
    if (!far) return x.w(m,Cy(dir),kf,jf,iff);
    const int up = (k == kf && j == jf && i == iff) ? 1 : 0;
    return x.w(m,Cy(dir),k+up*(dir == 2),j+up*(dir == 1),i+up*(dir == 0));
  }
  // F of the family's other face of cell (k,j,i) along face (kf,jf,iff), in our axes
  // (sp, Cartesian: a physical face's F^H); Near: of that face
  KOKKOS_INLINE_FUNCTION
  Real Far(const CleanArgs &x, const int dir, const int m, const int kf, const int jf,
           const int iff, const int k, const int j, const int i) const {
    if (!cub) {
      const int up = (k == kf && j == jf && i == iff) ? 1 : 0;
      const int ko = k + up*(dir == 2), jo = j + up*(dir == 1), io = i + up*(dir == 0);
      if (CleanBoundaryFace(x, dir, m, ko, jo, io)) return x.fc[dir](m,CLF_N,ko,jo,io);
    }
    Real s;
    const Real f = Slot(x, dir, m, kf, jf, iff, k, j, i, true, s);
    return s*f;
  }
  KOKKOS_INLINE_FUNCTION
  Real Near(const CleanArgs &x, const int dir, const int m, const int kf, const int jf,
            const int iff, const int k, const int j, const int i) const {
    Real s;
    const Real f = Slot(x, dir, m, kf, jf, iff, k, j, i, false, s);
    return s*f;
  }
  // whether the MeshBlock solves face (k,j,i), number t of the n + 1 of its line (tie:
  // the block whose low face it is; cubed sphere with clean_stag_fschur: the lower
  // gid's), else its identity row
  KOKKOS_INLINE_FUNCTION
  bool Own(const CleanArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t, const int n, const bool tie) const {
    if (!cub) return !(tie && t == n && !CleanBoundaryFace(x, dir, m, k, j, i));
    if (!tie || (t > 0 && t < n) || CleanBoundaryFace(x, dir, m, k, j, i)) return true;
    const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0), o = (t == 0) ? 1 : 0;
    const Real gc = x.g(m,CLG_F+7,k-(1-o)*dk,j-(1-o)*dj,i-(1-o)*di);
    const Real go = x.g(m,CLG_F+7,k-o*dk,j-o*dj,i-o*di);
    return gc < go || (gc == go && t == 0);
  }
  // F of face (k,j,i), number t of the n + 1 faces of its line along dir: get, put
  KOKKOS_INLINE_FUNCTION
  Real Get(const CleanArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t) const {
    if (!cub) return x.w(m,Cy(dir),k,j,i);
    return (t > 0) ? v(m,9+2*dir+1,k-(dir == 2),j-(dir == 1),i-(dir == 0))
                   : v(m,9+2*dir,k,j,i);
  }
  KOKKOS_INLINE_FUNCTION
  void Put(const CleanArgs &x, const int dir, const int m, const int k, const int j,
           const int i, const int t, const int n, const Real f) const {
    if (!cub) {
      x.w(m,Cy(dir),k,j,i) = f;
      return;
    }
    if (t > 0) v(m,9+2*dir+1,k-(dir == 2),j-(dir == 1),i-(dir == 0)) = f;
    if (t < n) v(m,9+2*dir,k,j,i) = f;
  }
};

//----------------------------------------------------------------------------------------
//! \fn bool CleanStagFRow
//! \brief the F solve's row lo F_{f-1} + dg F_f + hi F_{f+1} = r of face (k,j,i) normal
//! to dir; false where the face has no equation (F = F^H) or another block owns it (own
//! false): the identity row.  sm (clean_stag_fschur): the row sum, without cancellation.
//! The kinetic stress over h c is q (P_R - P_L) + qs (P_R + P_L): spherical polar's area
//! form (1/dx +- s), Cartesian 1/dx, cubed sphere (A_R P_R - A_L P_L)/V_f.

template <class FC>
KOKKOS_INLINE_FUNCTION
bool CleanStagFRow(const CleanArgs &x, const CleanFLine &ln, const FC &fc, const int dir,
                   const int m, const int k, const int j, const int i, Real &lo, Real &hi,
                   Real &dg, Real &r, const bool own = true, Real *sm = nullptr) {
  const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
  lo = 0.0, hi = 0.0, dg = 1.0, r = fc(m,CLF_N,k,j,i);
  if (sm) *sm = 1.0;
  if (!own) return false;
  const bool bd = CleanBoundaryFace(x, dir, m, k, j, i);
  Real tl[6], tr[6];
  ln.Cell(x, dir, m, k, j, i, k-dk, j-dj, i-di, tl);
  ln.Cell(x, dir, m, k, j, i, k, j, i, tr);
  const bool eq = !bd && (tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0);
  if (eq) {
    const Real g = 0.5*x.h*x.c/CleanCellFaceDistance(x, dir, m, k, j, i);
    const Real it_ = 1.0/fc(m,CLF_TH,k,j,i);
    // fkin: h c^2 [w P_kin]_L^R/dx, P_kin = (h+ F_lo + h- F_hi)/c of a cell; the cubed
    // sphere's cells carry their area pressures A P (CleanFaceMet)
    const CleanFaceMet fm = CleanFaceMetric(x, dir, m, k, j, i);
    const Real q = x.h*x.c*(x.cub ? 1.0 : fm.a)/fm.v, qs = x.h*x.c*fm.s;
    const Real ag = x.cub ? fm.a : 1.0;
    lo = -g*tl[0] - q*tl[2]*tl[3] + qs*tl[2]*tl[3];
    hi = -g*tr[0] + q*tr[2]*tr[4] + qs*tr[2]*tr[4];
    const Real sq = q*(tr[2]*(tr[3] + tr[4]) - tl[2]*(tl[3] + tl[4]));
    const Real ss = qs*(tr[2]*(tr[3] + tr[4]) + tl[2]*(tl[3] + tl[4]));
    dg = it_ - lo - hi + sq + ss;
    if (sm) *sm = it_ + sq + ss;
    r = fc(m,CLF_N,k,j,i)*it_ + q*x.c*(ag*(tl[1] - tr[1]) + tl[5] - tr[5]) -
        qs*x.c*(tl[1] + tr[1]) + x.h*x.c*x.c*fc(m,CLF_N+1,k,j,i);
  }
  return eq;
}
//----------------------------------------------------------------------------------------
//! \fn Real CleanStagFRay
//! \brief clean_stag_fwb: the rays' own lab-frame momentum of face (kf,jf,iff): their
//! flux factor (clamped to |g/J| <= 1) times the donor's E (work vector e), plus the
//! upwind advective a E.  Where the donor holds no rays (CLN_J = 0) the thin share's own
//! energy flux c (h+ E_L + h- E_R): a function of E on every face, never the last F, so
//! that the well-balanced rows of a line hold no lagged neighbour

KOKKOS_INLINE_FUNCTION
Real CleanStagFRay(const CleanArgs &x, const int e, const int dir, const int m,
                   const int kf, const int jf, const int iff) {
  const int ek = (dir == 2), ej = (dir == 1), ei = (dir == 0);
  const Real g = CleanRayFace(x, dir, m, kf, jf, iff);
  auto ec = [&](const int l) {return fmax(x.w(m,e,kf-l*ek,jf-l*ej,iff-l*ei), 0.0);};
  const int u = (g > 0.0) ? 1 : 0;
  const Real jd = x.cs(m,CLN_J,kf-u*ek,jf-u*ej,iff-u*ei);
  const Real a = 0.5*(CleanStagN(x, CLW_A1, dir, m, kf, jf, iff, kf-ek, jf-ej, iff-ei) +
                      CleanStagN(x, CLW_A1, dir, m, kf, jf, iff, kf, jf, iff));
  const Real f = (jd > 0.0) ? x.c*fmin(fmax(g/jd, -1.0), 1.0)*ec(u) :
                 x.c*(x.fc[dir](m,CLF_HP,kf,jf,iff)*ec(1) +
                      x.fc[dir](m,CLF_HM,kf,jf,iff)*ec(0));
  return f + a*ec(a > 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagFOwn
//! \brief clean_stag_fschur: a MeshBlock's line (m, a, b) along dir (n cells) as segment
//! s of the F line: first own face x = s0 + s1 yL + s2 xR, last own face y = s4 + s5 yL +
//! s6 xR (yL the face before, xR the face after), s3 = 1 - s1 - s2 and s7 = 1 - s5 - s6
//! from the row sums (no cancellation on nearly singular lines); Thomas sweeps up for y,
//! down for x.  layout_invariant_mode = tile: the tile of n cells from cell o, il / ih
//! its ends inside the MeshBlock, owned as between MeshBlocks.

template <class FC>
KOKKOS_INLINE_FUNCTION
void CleanStagFOwn(const CleanArgs &x, const CleanFLine &ln, const FC &fc, const int dir,
                   const int m, const int a, const int b, const int n, Real s[8],
                   const int o = 0, const bool il = false, const bool ih = false) {
  auto face = [&](const int t, int &k, int &j, int &i) {
    k = (dir == 2) ? x.ks + o + t : x.ks + b;
    j = (dir == 1) ? x.js + o + t : ((dir == 2) ? x.js + b : x.js + a);
    i = (dir == 0) ? x.is + o + t : x.is + a;
  };
  int k, j, i;
  Real lo, hi, dg, r, sm;
  face(0, k, j, i);
  const int t0 = (il ? !x.cub : ln.Own(x, dir, m, k, j, i, 0, n, true)) ? 0 : 1;
  face(n, k, j, i);
  const bool bh = !ih && CleanBoundaryFace(x, dir, m, k, j, i);   // F^H
  const int t1 = (x.cub && !bh && (ih || ln.Own(x, dir, m, k, j, i, n, n, true))) ? n
                                                                                : n - 1;
  const Real fh = bh ? fc(m,CLF_N,k,j,i) : 0.0, wh = bh ? 1.0 : 0.0;
  Real cp = 0.0, dc = 0.0, dl = 1.0, dw = 0.0;   // F_{-1} = yL
  for (int t = t0; t <= t1; ++t) {
    face(t, k, j, i);
    CleanStagFRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, true, &sm);
    const Real den = dg - lo*cp;
    dc = (r - lo*dc)/den;
    dl = -lo*dl/den;
    dw = (sm - lo*dw)/den;
    cp = hi/den;
  }
  s[4] = dc - cp*fh, s[5] = dl, s[6] = bh ? 0.0 : -cp, s[7] = dw - cp*wh;
  Real cq = 0.0, ec = fh, er = 1.0 - wh, ew = wh;   // F_{t1+1} = xR or F^H
  for (int t = t1; t >= t0; --t) {
    face(t, k, j, i);
    CleanStagFRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, true, &sm);
    const Real den = dg - hi*cq;
    ec = (r - hi*ec)/den;
    er = -hi*er/den;
    ew = (sm - hi*ew)/den;
    cq = lo/den;
  }
  s[0] = ec, s[1] = -cq, s[2] = er, s[3] = ew;
}

//----------------------------------------------------------------------------------------
//! \fn void CleanStagFJoin
//! \brief clean_stag_fschur: the segment c of two adjacent ones p, q (CleanStagFOwn), and
//! p's last face ya and q's first xb in terms of c's outer faces; s3, s7 as sums (no
//! cancellation)

KOKKOS_INLINE_FUNCTION
void CleanStagFJoin(const Real p[8], const Real q[8], Real c[8], Real ya[3],
                    Real xb[3]) {
  const Real d = p[5] + p[7] + p[6]*(q[2] + q[3]);
  ya[0] = (p[4] + p[6]*q[0])/d, ya[1] = p[5]/d, ya[2] = p[6]*q[2]/d;
  xb[0] = (q[0] + q[1]*p[4])/d, xb[1] = q[1]*p[5]/d, xb[2] = q[2]/d;
  c[0] = p[0] + p[2]*xb[0], c[1] = p[1] + p[2]*xb[1], c[2] = p[2]*xb[2];
  c[3] = p[3] + p[2]*(q[3] + q[1]*p[7])/d;
  c[4] = q[4] + q[5]*ya[0], c[5] = q[5]*ya[1], c[6] = q[6] + q[5]*ya[2];
  c[7] = q[7] + q[5]*(p[7] + p[6]*q[3])/d;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFSolve
//! \brief clean_stag_fsolve: with E of the solution, the face momentum F alone, one
//! unknown per face:
//!   (F - G)/h + c^2 [P_R - P_L]/dx_f = -c chi_f (F - A) (+ beta terms),
//!   P_i = (D E)_i - al_i [F]_i/(2c),  [F]_i the jump of cell i's two faces,
//! al_i = max(A(tau_i)^2, w of its faces) (1 where |F^H| > c E_donor, not with fkin); al
//! = 0 gives F^H.  Lines are tridiagonal (Thomas per MeshBlock line, ends lagged and
//! iterated with the ghost exchange).  CLF_N holds F, CLF_N + 1 the dissipation term
//! (al_R [F]_R - al_L [F]_L)/(2 c dx_f) the gas takes.
//! clean_stag_fkin: the streaming share w_i takes the rays' kinetic upwind flux,
//!   P_i = (1 - w_i) [(D E)_i - al_i [F]_i/(2c)] + w_i (h+_i F_lo + h-_i F_hi)/c,
//! in area form on spherical polar (F r^2 of a beam steady); telescoping, M-matrix
//! (h+ >= 0 >= h-).  Cubed sphere: a line's cell terms are the face family's (CleanFLine,
//! CleanCsFamily), area form (A_R P_R - A_L P_L)/V_f; seams iterate (up to 64 sweeps,
//! 1e-14), a shared face takes the lower gid's F (bitwise); lateral stress stays D E.

#if VET_TU_PART(2)
// layout_invariant_mode = tile: per MeshBlock line of nk tiles, the tiles' segments, the
// folds from the left and right (8 nk each) and the tiles' ends (2 nk)
static void CleanFsjt(DvceArray3D<Real> &v, const int nmb, const int nl, const int nk) {
  if (v.extent_int(0) == nmb && v.extent_int(1) == nl && v.extent_int(2) == 26*nk) return;
  v = DvceArray3D<Real>("cln_fsjt", nmb, nl, 26*nk);
}

// the F slots 9 + 2 q + hi of cln_v the ghosts of a face of family q take
// (CleanGSeg::fm): its own family's; f0 = 2: without the x1 family's
static std::array<uint32_t, 4> CleanFMask(const int f0) {
  return (f0 > 0) ? std::array<uint32_t, 4>{0, 3, 014, 0}
                  : std::array<uint32_t, 4>{3, 014, 060, 0};
}

void RadiationVet::CleanStagFSolve(const int e, const bool tie_now, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // clean_stag_fwb_stage: fwb's s of E_old (sm, the call with CLW_EO) before a round's
  // first solve (clean_init steady: every solve).  clean_stag_fwb_t2: s of the stage-time
  // state E* = E_old + h K (K the last slope, positive, within the first iterate's cap
  // implicit_predictor_tcap^4, else E_old), so F_ray and s are second order in time;
  // E_old is kept in CLW_V and copied back
  const int sw = (cln_fwb && !cln_insteady) ? cln_fwbs : -1;
  const bool sm = sw >= 0 && e == CLW_EO;
  const int ef = cln_fwbst ? CLW_EO : e;   // F_ray's E: fwb_stage's E_old
  if (sw >= 0 && !cln_fwbok) {
    cln_fwbok = true;
    const bool st = cln_fwbt2 && (t2_solve == VET_T2S_STAGE2 ||
                                  (t2_solve == VET_T2S_STAGE1 && t2_ok));
    auto w = cln_w;
    auto kv = (t2_solve == VET_T2S_STAGE1) ? t2k1 : t2k2;
    const Real h = dt_sub, ecap = SQR(SQR(impl_pred_tcap));
    const int ng = indcs.ng, n1 = indcs.nx1 + 2*ng;
    const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
    const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
    if (st) {
      CleanFor("cln_fwbt2", nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real eo = w(m,CLW_EO,k,j,i);
        w(m,CLW_V,k,j,i) = eo;
        if (k < ks || k > ke || j < js || j > je || i < is || i > ie) return;
        const Real ep = CleanPredictE(eo, h*kv(m,VET_T2_E,k,j,i));
        w(m,CLW_EO,k,j,i) = (ep > eo/ecap && ep < eo*ecap) ? ep : eo;
      });
      CleanGhostFill(cln_w, CLW_EO, 1);
    }
    CleanStagFSolve(CLW_EO, false, true);
    if (st) {
      CleanFor("cln_fwbt2r", nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        w(m,CLW_EO,k,j,i) = w(m,CLW_V,k,j,i);
      });
    }
    if (!lat) CleanStagCsLateral(e);
  }
  // lat false: the caller's lateral sums are those of e
  if (lat) CleanStagCsLateral(e);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  auto op = opac;
  const CleanFLine ln{cln_v, x.cub};
  const bool fk = cln_fkin, dok = cln_dok, cub = x.cub;
  const bool rms = cln_rmsh && cln_rstold && !cln_rmdone && fk && e == CLW_X && dok;
  DvceArray1D<Real> rmc("cln_rmc", 1);
  const bool fwb = cln_fwb;
  if (fwb && cln_fwbc.extent_int(0) != 3) Kokkos::realloc(cln_fwbc, 3);
  if (fwb && (sw < 0 || sm)) Kokkos::deep_copy(cln_fwbc, 0.0);
  auto fwc = cln_fwbc;
  const bool wbc = e == CLW_X;
  const int kd = (t2_solve == VET_T2S_NONE) ? 0 :
                 ((t2_solve == VET_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD);
  const Real c1lo = cln_c0, c1hi = cln_c1;
  // clean_stag_fs_line = pcr; on the cubed sphere with fschur its x1 lines, in panel
  const bool fpcr = cln_fspcr && (!cub || cln_fschur) &&
      !std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>;
  // the F of every direction (CLW_Y, CLW_Z, CLW_S or cln_v) in one exchange, started
  // after the last direction.  Cubed sphere, one MeshBlock along x1: the x1 family's
  // slots (0..2, 9, 10, 15) have no line across a MeshBlock face: not exchanged
  const int f0 = (cub && pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel)
                 ? 2 : 0;
  auto yfill = [&](const int part) {
    if (cub) {   // CleanFLine reads the faces' own cells, of the face's family (fm)
      CleanGhostFill(cln_v, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));
      return;
    }
    if (x.nf == 3) CleanGhostDefer(cln_w, CLW_S, 1);
    CleanGhostFill(cln_w, CLW_Y, std::min(x.nf, 2), part);
  };
  // three phases (F^H and al; line solves; F), each direction with its own al and F, so
  // each exchange serves all.  clean_stag_fsweeps N: up to N line solves, each with the
  // line ends' F of the last (converging across MeshBlock faces and the periodic wrap);
  // the cubed sphere without clean_stag_fschur iterates its seams (64 by default)
  const int nsw = std::max(cln_fsweeps, 1);
  const bool tie = !cub || cln_fschur;
  const bool fov = !cub && nsw == 1;
  Real fsd = 0.0, fsf = 0.0;
  for (int q = 0; q < (2 + nsw)*x.nf; ++q) {
    const int dir = q%x.nf, ph = std::min(q/x.nf, 1) + (q >= (1 + nsw)*x.nf);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    if (ph == 0) {
      // F^H of every face (CLF_N) and the cell terms of every cell (CleanFLine)
      const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
      CleanFor<4>("cln_fs0", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
        CleanDirCase(dir, [&](auto d) {
          constexpr int n = decltype(d)::value;
          fc(m,CLF_N,k,j,i) = CleanFaceFlux<n>(xh, e, n, m, k, j, i, true);
        });
        fc(m,CLF_PT,k,j,i) = fc(m,CLF_N,k,j,i);   // F^H for CleanStagDrag (QL is free)
        fc(m,CLF_N+1,k,j,i) = (sw >= 0 && !sm) ? fc(m,sw,k,j,i) : 0.0;
        if (sm) fc(m,sw,k,j,i) = 0.0;
        if (!fk || !x.sph || dir == 2 || CleanBoundaryFace(x, dir, m, k, j, i)) return;
        // spherical polar: the w share's transverse pressures (but not the normal
        // stress's, q = 4) are the rays' own (CLN_G) where the cell holds their J
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
      auto v = cln_v;
      auto f1 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rc;
        const Real tau = op(m,VET_OP_T,k,j,i)*CleanCellFaceDistance(x, dir, m, k, j, i);
        // fkin: the streaming share is kinetic, without dissipation
        Real a = fk ? 0.0 : CleanDissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
        for (int s = 0; s < 2; ++s) {
          const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
          if (!fk) a = fmax(a, fc(m,CLF_W,kf,jf,iff));
          const Real f = fc(m,CLF_N,kf,jf,iff);
          const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
          if (!fk && !CleanBoundaryFace(x, dir, m, kf, jf, iff) &&
              !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
          // cubed sphere faces: the last F solve's
          if (cub) v(m,9+2*dir+s,k,j,i) = fc(m,CLF_N+2,kf,jf,iff);
        }
        // fkin: w of the cell; the part of P - D E without F, -w D E, and the stress the
        // kinetic form misses for the rays' own state, (D - s^2) J (0 for a beam, J/3
        // isotropic), lagged; cubed sphere: D about the family's mean normal, area form
        const Real wc = fk ? CleanSmoothStep((x.cs(m,CLN_C1,k,j,i) - c1lo)/(c1hi - c1lo))
                           : 0.0;
        const Real ec = x.w(m,e,k,j,i);
        const Real sr = x.cs(m,CLN_HP+dir,k,j,i) + x.cs(m,CLN_HM+dir,k,j,i);
        const Real jc = fmin(x.cs(m,CLN_J,k,j,i), fmax(ec, 0.0));
        if (!cub) {
          const Real dd = x.cs(m,CLN_D11+dir,k,j,i);
          x.w(m,CLW_L+3*dir,k,j,i) = (1.0 - wc)*a;
          x.w(m,CLW_L+3*dir+1,k,j,i) = wc*((dd - sr*sr)*jc - dd*ec);
          x.w(m,CLW_L+3*dir+2,k,j,i) = wc;
          // F of the low face, for the exchange (multi-D: the last F solve's)
          x.w(m,CleanFLine::Cy(dir),k,j,i) = fc(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,
                                                j,i);
          return;
        }
        Real nd[3][3], dd = 0.0;
        CleanCsFamily(x, m, k, j, i, nd);
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            dd += nd[dir][b]*nd[dir][c]*x.cs(m,CLN_D11+((b == c) ? b : 2 + b + c),k,j,i);
          }
        }
        const Real ac = 0.5*(x.area[dir](m,k,j,i) + x.area[dir](m,k+dk,j+dj,i+di));
        v(m,3*dir,k,j,i) = (1.0 - wc)*a;
        v(m,3*dir+1,k,j,i) = wc*dd*(jc - ec);
        v(m,3*dir+2,k,j,i) = ac*wc;
        v(m,15+dir,k,j,i) = -ac*wc*sr*sr*jc;
      };
      // cubed sphere, one MeshBlock along x1: one depth-1 exchange of the family terms
      // (fm), started once the last direction's sent cells are done
      const bool pv = cub && f0 > 0 && dir == x.nf - 1;
      const std::function<void()> post = [&]() {
        CleanGhostDefer(cln_v, 3, 6, 1);
        cln_hpend.back().fm = {0, 7, 070, 0};
        CleanGhostDefer(cln_v, 11, 4, 1);
        cln_hpend.back().fm = {0, 3, 014, 0};
        CleanGhostFill(cln_v, 16, 2, true, 1, 1, {0, 1, 2, 0});
      };
      const std::function<void()> fin = [&]() {CleanGhostFill(cln_v, 16, 2, true, 2, 1);};
      if (pv) {
        CleanOverlapFor("cln_fs1", pmy_pack, true, f1, fin, -1, post);
      } else {
        CleanFor("cln_fs1", nmb1, ks, ke, js, je, is, ie, f1);
      }
      if (dir == x.nf - 1 && !pv) {
        if (cub) {
          CleanGhostFill(cln_v, 0, 18, true, 0, 1);
        } else {
          CleanGhostDefer(cln_w, CLW_L, 3*x.nf);
          yfill(0);
        }
      }
    }
    // clean_stag_fschur: the faces beyond the line ends exact, not lagged.
    // clean_stag_fwb: each row less its streaming share w of its residual at the rays'
    // own state F_ray, so F_ray is the thin rows' exact steady state.  F_ray (lab frame):
    // the rays' flux factor times the donor's old E plus the upwind advective a E; a
    // physical face its F^H; also the memory G.  The thick share keeps its row, s takes
    // 1/(1 + tau_f), tau_f = chi_f dx_f.  The residual s enters CLF_N + 1: the gas takes
    // the absorption of F, not s, whose s V/c^2 is booked as the radiation's (L_z: dlz).
    // E and the energy flux memory G + d are unchanged.  clean_stag_fwb_stage: s once
    // per round (sm) into slot sw, added in phase 0
    if (fwb && q == x.nf && (sw < 0 || sm)) {
      const CleanArgs *pw = CleanDevArgs(x, 0);
      // clean_stag_fwb_line, _ring: s and its weight of the projected faces
      // (CleanStagFwbProject)
      DvceArray4D<Real> rs, rw;
      if (cln_fwbpj[0] || cln_fwbpj[2]) {
        rs = ScratchView<DvceArray4D<Real>>("cln_fwrs", nmb1 + 1, u0.extent_int(2),
                                            u0.extent_int(3), u0.extent_int(4));
        rw = ScratchView<DvceArray4D<Real>>("cln_fwrw", nmb1 + 1, u0.extent_int(2),
                                            u0.extent_int(3), u0.extent_int(4));
      }
      auto sl = ScratchView<DvceArray1D<Real>>("cln_csl", 6*CLN_NSL);
      if (wbc || sm) Kokkos::deep_copy(sl, 0.0);
      // the F slots of the lines (CleanFLine) hold F_ray, then the solve's start again
      for (int pass = 0; pass < 2; ++pass) {
        for (int d = 0; d < x.nf; ++d) {
          const int ek = (d == 2), ej = (d == 1), ei = (d == 0);
          auto fd = x.fc[d];
          CleanFor("cln_fwbs", nmb1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            const CleanArgs &x = *pw;
            for (int s = 0; s < (cub ? 2 : 1); ++s) {
              const int kf = k + s*ek, jf = j + s*ej, iff = i + s*ei;
              Real f = cub ? fd(m,CLF_N+2,kf,jf,iff)
                           : fd(m,(fk && x.multi_d) ? CLF_N+2 : CLF_N,k,j,i);
              if (pass == 0 && !CleanBoundaryFace(x, d, m, kf, jf, iff)) {
                f = CleanStagFRay(x, ef, d, m, kf, jf, iff);
              }
              if (cub) {
                ln.v(m,9+2*d+s,k,j,i) = f;
              } else {
                x.w(m,CleanFLine::Cy(d),k,j,i) = f;
              }
            }
          });
        }
        yfill(0);
        if (pass == 1) break;
        for (int d = 0; d < x.nf; ++d) {
          const int ek = (d == 2), ej = (d == 1), ei = (d == 0);
          auto fd = x.fc[d];
          const bool pj = cln_fwbpj[d];
          if (pj) {
            Kokkos::deep_copy(rs, 0.0);
            Kokkos::deep_copy(rw, 0.0);
          }
          CleanFor("cln_fwb", nmb1, ks, ke+ek, js, je+ej, is, ie+ei,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            const CleanArgs &x = *pw;
            Real lo, hi, dg, r;
            if (CleanBoundaryFace(x, d, m, k, j, i)) return;
            if (!CleanStagFRow(x, ln, fd, d, m, k, j, i, lo, hi, dg, r)) return;
            const Real f0 = CleanStagFRay(x, ef, d, m, k, j, i);
            const Real fm = ln.Far(x, d, m, k, j, i, k-ek, j-ej, i-ei);
            const Real fp = ln.Far(x, d, m, k, j, i, k, j, i);
            const Real qd = x.drag ? fd(m,CLF_DB,k,j,i) : 0.0;
            const Real tau = 0.5*(x.op(m,VET_OP_T,k,j,i) +
                                  x.op(m,VET_OP_T,k-ek,j-ej,i-ei))*
                             CleanCellFaceDistance(x, d, m, k, j, i);
            const Real s = fd(m,CLF_W,k,j,i)/(1.0 + tau)*
                (lo*fm + dg*f0 + hi*fp - r - (1.0 - qd)*(f0 - fd(m,CLF_G,k,j,i)));
            if (pj) {
              rs(m,k,j,i) = s;
              rw(m,k,j,i) = fd(m,CLF_W,k,j,i)/(1.0 + tau);
              return;
            }
            fd(m,sm ? sw : CLF_N+1,k,j,i) += s/(x.h*x.c*x.c);
            const Real p = fd(m,CLF_TH,k,j,i)*s*0.5*(CleanVol(x, m, k, j, i) +
                           CleanVol(x, m, k-ek, j-ej, i-ei))/(x.c*x.c);
            if (!wbc && !sm) return;   // counters: the printed s, the write-back's
            CleanSlotAdd(sl, 0, m, k, j, i, 1.0);
            CleanSlotAdd(sl, 1, m, k, j, i, fabs(p));
            if (d == 0) CleanSlotAdd(sl, 2, m, k, j, i, p);
          });
          if (pj) CleanStagFwbProject(d, sm ? sw : CLF_N + 1, rs, rw, wbc || sm);
        }
      }
      if (wbc || sm) CleanSlotSum(sl, fwc, 3);
    }
    if (sm && q == x.nf) return;
    if (cln_fschur && q == x.nf) CleanStagFSchur();
    const int n = (dir == 0) ? ie - is + 1 : ((dir == 1) ? je - js + 1 : ke - ks + 1);
    const int na = (dir == 0) ? je - js + 1 : ie - is + 1;
    const int nb = (dir == 2) ? je - js + 1 : ke - ks + 1;
    Real dmax = 0.0, fmx = 0.0;
    if (ph == 1) {
      // per face (coalesced): lo (CLW_AM), hi (CLW_AP), dg (CLW_A0, then the Thomas den),
      // r (CLW_R0), faces with an equation (CLW_P), and in the sweep c' (CLW_CP); E-solve
      // scratch slots shared by the directions (the x1 PCR reads them too)
      CleanFor("cln_fs2f", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const int t = (dir == 0) ? i - is : ((dir == 1) ? j - js : k - ks);
        // face tie: a line's high face shared with another MeshBlock or the periodic
        // wrap is its owner's (whose low face it is): an identity row (CleanStagFaceTie)
        const bool own = ln.Own(x, dir, m, k, j, i, t, n, tie);
        Real lo, hi, dg, r;
        const bool eq = CleanStagFRow(x, ln, fc, dir, m, k, j, i, lo, hi, dg, r, own);
        x.w(m,CLW_AM,k,j,i) = lo;
        x.w(m,CLW_AP,k,j,i) = hi;
        x.w(m,CLW_A0,k,j,i) = dg;
        x.w(m,CLW_R0,k,j,i) = r;
        x.w(m,CLW_P,k,j,i) = !own ? -1.0 : (eq ? 1.0 : 0.0);
      });
    }
    if (ph == 1 && dir == 0 && fpcr) {
      // clean_stag_fs_pcr: on a GPU the x1 lines by parallel cyclic reduction, one team
      // per line (as CleanStagLineSolve); F to round-off of the Thomas solve
      const int nx = n + 1;
      int nround = 0;
      while ((1 << nround) < nx) ++nround;
      const size_t scr = ScrArray1D<Real>::shmem_size(8*nx);
      Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), (nmb1 + 1)*na*nb, 256);
      const int lv0 = TeamScratchLevel(scr, 0, pol.team_size());
      Kokkos::parallel_for("cln_fs2p", pol.set_scratch_size(lv0, Kokkos::PerTeam(scr)),
      KOKKOS_LAMBDA(TeamMember_t tm) {
        const int m = tm.league_rank()/(na*nb), a = tm.league_rank()%na;
        const int k = ks + (tm.league_rank()/na)%nb, j = js + a;
        ScrArray1D<Real> sw(tm.team_scratch(lv0), 8*nx);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int t) {
          const int i = is + t;
          // rows of cln_fs2f; faces beyond a MeshBlock line end: the neighbours', lagged;
          // ftie: a shared high face its owner's, the identity row
          const Real lo = x.w(m,CLW_AM,k,j,i), hi = x.w(m,CLW_AP,k,j,i);
          const Real dg = x.w(m,CLW_A0,k,j,i), p = x.w(m,CLW_P,k,j,i);
          Real r = x.w(m,CLW_R0,k,j,i);
          if (p > 0.0) {
            if (t == 0) r -= lo*ln.Far(x, dir, m, k, j, i, k, j, i-1);
            if (t == n) r -= hi*ln.Far(x, dir, m, k, j, i, k, j, i);
          }
          if (p < 0.0) r = ln.Near(x, dir, m, k, j, i, k, j, i-(t == 0));
          sw(t) = (t > 0) ? lo : 0.0;
          sw(nx + t) = dg;
          sw(2*nx + t) = (t < n) ? hi : 0.0;
          sw(3*nx + t) = r;
        });
        tm.team_barrier();
        int o = 0;
        for (int rd = 0, s = 1; rd < nround; ++rd, s *= 2) {
          const int d = 4*nx - o;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
            Real a = 0.0, b = sw(o + nx + q), c = 0.0, r = sw(o + 3*nx + q);
            if (q >= s) {
              const Real f = -sw(o + q)/sw(o + nx + q - s);
              a = f*sw(o + q - s);
              b += f*sw(o + 2*nx + q - s);
              r += f*sw(o + 3*nx + q - s);
            }
            if (q + s < nx) {
              const Real g = -sw(o + 2*nx + q)/sw(o + nx + q + s);
              c = g*sw(o + 2*nx + q + s);
              b += g*sw(o + q + s);
              r += g*sw(o + 3*nx + q + s);
            }
            sw(d + q) = a;
            sw(d + nx + q) = b;
            sw(d + 2*nx + q) = c;
            sw(d + 3*nx + q) = r;
          });
          tm.team_barrier();
          o = d;
        }
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int t) {
          ln.Put(x, dir, m, k, j, is+t, t, n, sw(o + 3*nx + t)/sw(o + nx + t));
        });
      });
    } else if (ph == 1) {
      // one exact Thomas solve per line, the faces beyond a MeshBlock line end or
      // periodic wrap the last F solve's (lagged); cubed sphere: also the largest face
      // change and largest F.  layout_invariant_mode = tile: each tile of a lateral line
      // solved as a MeshBlock line, its inner ends from cln_fsjt, joined
      // (CleanStagFSchurTile) or the last sweep's; a face two tiles share is its
      // non-owner's identity row
      const int nk = (dir > 0 && repro::tile_n > 0) ? n/repro::tile_n : 1, nt = n/nk;
      if (nk > 1) CleanFsjt(cln_fsjt[dir], nmb1 + 1, na*nb, nk);
      auto jt = cln_fsjt[dir];
      const bool fsc = cln_fschur, olo = !(cub && tie), ohi = cub || !tie;
      auto line = KOKKOS_LAMBDA(const int idx, Real &d_max, Real &f_max) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        for (int u = 1; u < nk && !fsc; ++u) {
          int k, j, i;
          face(u*nt, k, j, i);
          jt(m,l,25*nk+u-1) = ohi ? ln.Far(x, dir, m, k, j, i, k, j, i)
                                  : ln.Near(x, dir, m, k, j, i, k, j, i);
          jt(m,l,24*nk+u) = olo ? ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di)
                                : ln.Near(x, dir, m, k, j, i, k-dk, j-dj, i-di);
        }
        for (int u = 0; u < nk; ++u) {
          const int o = u*nt;
          const bool il = u > 0, ih = u < nk - 1;
          Real dp = 0.0, cp = 0.0;
          for (int t = 0; t <= nt; ++t) {
            int k, j, i;
            face(o + t, k, j, i);
            const bool ie = (il && t == 0) || (ih && t == nt);
            const bool io = ie && !((t == 0) ? olo : ohi);
            const Real ee = ie ? jt(m,l,((t == 0) ? 24*nk : 25*nk)+u) : 0.0;
            const Real am = io ? 0.0 : x.w(m,CLW_AM,k,j,i);
            const Real ap = io ? 0.0 : x.w(m,CLW_AP,k,j,i);
            const Real pp = io ? -1.0 : x.w(m,CLW_P,k,j,i);
            Real den = io ? 1.0 : x.w(m,CLW_A0,k,j,i);   // dg
            // faces beyond a MeshBlock line end: the neighbours', lagged
            const bool eq = pp > 0.0;
            den -= ((eq && t == 0) ? 0.0 : am)*cp;
            cp = ((eq && t == nt) ? 0.0 : ap)/den;
            x.w(m,CLW_CP,k,j,i) = cp;
            const int t0 = (t == 0) ? 1 : 0;
            Real r = (pp < 0.0) ? (ie ? ee : ln.Near(x, dir, m, k, j, i, k-t0*dk, j-t0*dj,
                                                     i-t0*di)) :
                x.w(m,CLW_R0,k,j,i);
            if (pp > 0.0) {
              if (t == 0) {
                r -= am*(il ? ee : ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di));
              }
              if (t == nt) r -= ap*(ih ? ee : ln.Far(x, dir, m, k, j, i, k, j, i));
            }
            dp = (r - ((t == 0) ? 0.0 : am)*dp)/den;
            x.w(m,CLW_R,k,j,i) = dp;
          }
          Real fn = 0.0;
          // the last face's slot is a ghost: the neighbour's or a physical face
          for (int t = nt; t >= 0; --t) {
            int k, j, i;
            face(o + t, k, j, i);
            fn = x.w(m,CLW_R,k,j,i) - ((t < nt) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
            if (cub) {
              d_max = fmax(d_max, fabs(fn - ln.Get(x, dir, m, k, j, i, t)));
              f_max = fmax(f_max, fmax(fabs(fn), 1.0e-8*x.c*fabs(x.w(m,e,k,j,i))));
            }
            ln.Put(x, dir, m, k, j, i, t, nt, fn);
          }
        }
      };
      // one warp or wavefront per block: few x1 lines, one thread each, latency bound
#if defined(KOKKOS_ENABLE_HIP)
      using LineBounds = Kokkos::LaunchBounds<64,1>;
#else
      using LineBounds = Kokkos::LaunchBounds<32,1>;
#endif
      const auto pol = Kokkos::RangePolicy<DevExeSpace, LineBounds>(DevExeSpace(), 0,
                                                                    (nmb1 + 1)*na*nb);
      if (cub && nsw > 1) {
        AthReduce("cln_fs2", pol, line, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(fmx));
      } else {
        AthFor("cln_fs2", pol, KOKKOS_LAMBDA(const int idx) {
          Real d_max = 0.0, f_max = 0.0;
          line(idx, d_max, f_max);
        });
      }
    }
    if (ph == 1) {
      fsd = (dir == 0) ? dmax : fmax(fsd, dmax);
      fsf = (dir == 0) ? fmx : fmax(fsf, fmx);
    }
    if (ph == 1 && dir == x.nf - 1) {
      yfill(fov ? 1 : 0);
      // sweeps: the cubed sphere until no face changes by 1e-14 of the largest F, else
      // until the owners' high-face F changes by < 1e-12 of itself
      const bool done = nsw > 1 &&
          (cub ? CleanReduce(fsd, true) <= 1.0e-14*CleanReduce(fsf, true)
               : (CleanStagFSweepChange() < 1.0e-12 && q >= 2*x.nf));
      if (done) q = (1 + nsw)*x.nf - 1;
    }
    // cubed sphere: a face shared with another MeshBlock takes the lower gid's F, and the
    // cells' slots of it
    if (cub && q == (1 + nsw)*x.nf) CleanStagFCsShared();
    if (ph == 2) {
      // F and the dissipation term of every face
      const CleanArgs *px = CleanDevArgs(x, 0);
      auto f3 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const CleanArgs &x = *px;   // in device memory: no per-thread copy
        const Real f = cub ? fc(m,CLF_N+2,k,j,i) : x.w(m,CleanFLine::Cy(dir),k,j,i);
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);   // F^H: the boundary's, kept
        if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
        // fkin, the solution's F solve: d = Phi^H - F, Phi^H = F^H + theta d the energy
        // flux's high-order part, memory G + d (CleanFaceSource): its own history, so its
        // steady state does not depend on h
        if (fk && e == CLW_X && dok) {   // and its slope (CleanDInc)
          const Real th = fc(m,CLF_TH,k,j,i);
          const Real dd = fc(m,CLF_N,k,j,i) - f - (1.0 - th)*fc(m,CLF_N+3,k,j,i);
          fc(m,CLF_N+3,k,j,i) += dd;
          if (kd > 0) fc(m,kd,k,j,i) = dd/x.h;
          // clean_rst_mshare, the first solution after a restart from an old scheme's
          // file: the memory G + d at the fixed point of Phi^H = theta (G + d) + B of
          // this E, the VET share's steady flux B/(1 - theta) (above c E of its donor:
          // counted), so the restart has no relaxation transient
          if (rms && th < 1.0) {
            const Real ms = (fc(m,CLF_N,k,j,i) - th*fc(m,CLF_G,k,j,i))/(1.0 - th);
            if (fabs(ms) > x.c*fmax(x.w(m,e,k,j,i), x.w(m,e,k-dk,j-dj,i-di))) {
              Kokkos::atomic_add(&rmc(0), 1.0);
            }
            fc(m,CLF_N+3,k,j,i) = ms - f;
          }
        }
        fc(m,CLF_N,k,j,i) = f;
        fc(m,CLF_N+2,k,j,i) = f;
        Real tl[6], tr[6];
        ln.Cell(x, dir, m, k, j, i, k-dk, j-dj, i-di, tl);
        ln.Cell(x, dir, m, k, j, i, k, j, i, tr);
        if (!(tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0)) return;
        // a physical high face keeps F^H; the ghost fill replaced its slot
        const Real fm = ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di);
        const Real fp = ln.Far(x, dir, m, k, j, i, k, j, i);
        auto pk = [&](const Real t[6], const Real l, const Real h) {
          return t[2]*(t[3]*l + t[4]*h)/x.c + (cub ? t[5] : t[1]);
        };
        const Real pl = pk(tl, fm, f), pr = pk(tr, f, fp);
        // the kinetic stress: the cubed sphere's area form over V_f, else the F row's
        Real dp, ds = 0.0;
        if (cub) {
          dp = (pr - pl + CleanCsNormA(x, dir, m, k, j, i)*(tr[1] - tl[1]))/
               (0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,k-dk,j-dj,i-di)));
        } else {
          dp = (x.sph ? CleanStagSph(x, 0, dir, m, k, j, i) : x.idx[dir])*(pr - pl);
          ds = x.sph ? CleanStagSph(x, 4, dir, m, k, j, i)*(pr + pl) : 0.0;
        }
        const Real dc = 2.0*x.c*CleanCellFaceDistance(x, dir, m, k, j, i);
        fc(m,CLF_N+1,k,j,i) += (tr[0]*(fp - f) - tl[0]*(f - fm))/dc - dp - ds;
      };
      // x1: the faces that read no ghost while the exchange of F finishes
      const bool ov = fov && dir == 0;
      const std::function<void()> fin = [&]() {
        if (ov) CleanGhostFill(cln_w, CLW_Y, 1, 2);
      };
      CleanOverlapFor("cln_fs3", pmy_pack, ov, f3, fin, dir);
    }
  }
  if (cln_fcone) CleanStagFCone(e, kd);
  if (tie_now) CleanStagFaceTie(e, kd);
  if (fwb && e == CLW_X) {   // clean_stag_fwb: radiation momentum of the residuals
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_fwbc);
    const Real n = CleanReduce(h(0), false), p = CleanReduce(h(1), false);
    const Real p1 = CleanReduce(h(2), false);
    if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(p + p1))) {
      std::printf("CLEAN %sFWB cycle %d t2s %d faces %.0f radiation momentum sum |dp| "
                  "%.6e x1 faces sum dp %.6e\n", std::isfinite(p + p1) ? "" : "WARNING ",
                  pmy_pack->pmesh->ncycle, t2_solve, n, p, p1);
    }
  }
  if (rms) {   // clean_rst_mshare: once, at the old file's first solution
    cln_rmdone = true;
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rmc);
    const Real n = CleanReduce(h(0), false);
    if (global_variable::my_rank == 0) {
      std::printf(
          "<rad_vet> clean_rst_mshare RESTART CONVERSION (one time): the VET share's "
          "memory at its steady flux (faces above c E of the donor: %.0f)\n", n);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFCone
//! \brief clean_stag_fcone: the cell flux vector (CleanStagCellVector, F_i = sum_f F_f
//! u_f) within the realisable cone |F_i| <= c E_i (E of work vector e), by construction.
//! Each face is scaled towards 0 by l_f = min(l_L, l_R), l_i = min(1, c E_i/S_i), S_i =
//! sum_f |F_f| |u_f| >= |F_i|, so |F_i| <= l_i S_i <= c E_i.  One value per face (cubed
//! sphere: the lower gid's) keeps flux form.  The gas takes the absorption of the scaled
//! F (CLF_N + 1), G + d is kept, the momentum cut is booked and printed (CLEAN FCONE).

void RadiationVet::CleanStagFCone(const int e, const int kd) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const bool cub = x.cub, dd = cln_fkin && e == CLW_X && cln_dok, wbc = e == CLW_X;
  const bool fg = cln_fcgas;   // (clean_stag_fcone_gas)
  // (clean_stag_fcone_lz: l_i of the faces without L_z arm (fa), CLW_CP >= 0 with the
  // others' 1, or -l_i of the others with fa's 0, where those alone exceed c E;
  // clean_lz_spin: the others' cut withheld, its L_z X the gas's rigid rotation; the
  // cubed sphere's gas phi momentum is lz_tel's torque, which holds the cut already)
  const bool lz = cln_fclz, spn = cln_lzsp, sph = x.sph;
  auto fa = KOKKOS_LAMBDA(const int q) {return !CleanLzFamily(sph, cub, q);};
  auto lf = KOKKOS_LAMBDA(const Real cp, const int q) {
    return !lz ? cp : (fa(q) ? fmax(cp, 0.0) : ((cp < 0.0) ? -cp : 1.0));
  };
  const bool emx = cln_fcemx;   // (clean_stag_fcone_emax)
  if (cln_fconc.extent_int(0) != 5) Kokkos::realloc(cln_fconc, 5);
  if (spn && cln_fcx.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(cln_fcx, nmb1 + 1, 3, cln_w.s.extent_int(2) + 1,
                    cln_w.s.extent_int(3) + 1, cln_w.s.extent_int(4) + 1);
  }
  auto fcx = cln_fcx;
  Kokkos::deep_copy(cln_fconc, 0.0);
  auto fcc = cln_fconc;
  CleanFor("cln_fcone_l", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = 0.0, st = 0.0;   // (st: the faces with L_z arm, clean_stag_fcone_lz)
    for (int d = 0; d < 2*x.nf; ++d) {
      const int q = d/2, h = d%2;
      Real c[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, u[3];
      c[d] = 1.0;
      CleanStagCellVector(x, m, k, j, i, c, u);
      const Real a = fabs(x.fc[q](m,CLF_N,k+h*(q == 2),j+h*(q == 1),i+h*(q == 0)))*
                     sqrt(SQR(u[0]) + SQR(u[1]) + SQR(u[2]));
      if (lz && !fa(q)) {
        st += a;
      } else {
        s += a;
      }
    }
    const Real ce = x.c*fmax(x.w(m,e,k,j,i), 0.0);
    x.w(m,CLW_CP,k,j,i) = (st > ce) ? -ce/st : ((s + st > ce) ? (ce - st)/s : 1.0);
  });
  CleanGhostFill(cln_w, CLW_CP, 1);
  for (int d = 0; d < x.nf; ++d) {
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    auto fc = x.fc[d];
    CleanFor("cln_fcone_f", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // a physical boundary face: its inner cell's l, after the cubed sphere's tie; an
      // opaque face (tau_f = chi_f dx_f) keeps its share tau/(1 + tau) of F (a diffusion
      // front's face mean is no free-streaming excess)
      const bool bd = CleanBoundaryFace(x, d, m, k, j, i);
      const Real tau = 0.5*(x.op(m,VET_OP_T,k,j,i) + x.op(m,VET_OP_T,k-dk,j-dj,i-di))*
                       CleanCellFaceDistance(x, d, m, k, j, i);
      const Real l = fmin(lf(x.w(m,CLW_CP,k-dk,j-dj,i-di), d),
                          lf(x.w(m,CLW_CP,k,j,i), d));
      Real fs = (bd ? 1.0 : 1.0 - (1.0 - l)/(1.0 + tau))*fc(m,CLF_N,k,j,i);
      if (emx && !bd) {   // (clean_stag_fcone_emax: |F_f| <= c max(E_L, E_R))
        const Real ce = x.c*fmax(fmax(x.w(m,e,k,j,i), x.w(m,e,k-dk,j-dj,i-di)), 0.0);
        fs = fmin(fmax(fs, -ce), ce);
      }
      fc(m,CLF_N+2,k,j,i) = fs;
    });
  }
  if (cub) {   // a face two MeshBlocks share: its lower gid's value
    auto v = cln_v;
    CleanFor("cln_fcone_v", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int q = 0; q < 3; ++q) {
        v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
        v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
      }
    });
    const int f0 = (pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel)
                   ? 2 : 0;
    CleanGhostFill(cln_v, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));
    CleanStagFCsShared();
  }
  auto sl = ScratchView<DvceArray1D<Real>>("cln_csl", 6*CLN_NSL);
  if (wbc) Kokkos::deep_copy(sl, 0.0);
  for (int d = 0; d < x.nf; ++d) {
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    auto fc = x.fc[d];
    CleanFor("cln_fcone_b", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // a physical boundary face has no F equation: the gas takes its half cell's
      // absorption of fs (CleanMomentumChange)
      const bool bd = CleanBoundaryFace(x, d, m, k, j, i);
      const bool lin = (d == 0) ? (i > is) : ((d == 1) ? (j > js) : (k > ks));
      const Real f = fc(m,CLF_N,k,j,i);
      const Real fs = !bd ? fc(m,CLF_N+2,k,j,i) :
                      lf(x.w(m,CLW_CP,k-lin*dk,j-lin*dj,i-lin*di), d)*f;
      fc(m,CLF_N+2,k,j,i) = fs;
      if (spn) fcx(m,d,k,j,i) = 0.0;
      if (fs == f) return;
      const Real it = 1.0/fc(m,CLF_TH,k,j,i);
      // (clean_lz_spin: the arm faces' cut withheld from the gas per face, CleanWriteBack)
      if (spn && !bd && !fa(d) && !fg && !cub) fcx(m,d,k,j,i) = (f - fs)*it/(x.c*x.c);
      if (!bd && !fg && !(lz && !fa(d) && !(spn && !cub))) {
        fc(m,CLF_N+1,k,j,i) -= (f - fs)*it/(x.h*x.c*x.c);
      }
      if (dd && !bd) {   // (the energy flux's memory G + d kept)
        fc(m,CLF_N+3,k,j,i) += f - fs;
        if (kd > 0) fc(m,kd,k,j,i) += (f - fs)/x.h;
      }
      fc(m,CLF_N,k,j,i) = fs;
      if (!wbc) return;   // counters: the write-back's, as printed
      CleanSlotAdd(sl, 0, m, k, j, i, 1.0);
      const Real dp = (f - fs)*it*0.5*(CleanVol(x, m, k, j, i) +
                      CleanVol(x, m, k-dk, j-dj, i-di))/(x.c*x.c);
      CleanSlotAdd(sl, 1, m, k, j, i, fabs(dp));
      // (fcone_lz: the momentum removed along each face family, signed, the gas's not)
      if (lz && !bd && fa(d)) CleanSlotAdd(sl, 2 + d, m, k, j, i, dp);
    });
  }
  if (wbc) CleanSlotSum(sl, fcc, 5);
  if (e == CLW_X) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_fconc);
    const Real n = CleanReduce(h(0), false), p = CleanReduce(h(1), false);
    if (global_variable::my_rank == 0 && cln_verbose > 0 && n > 0.0) {
      std::printf("CLEAN FCONE cycle %d t2s %d faces %.0f momentum cut sum |dF| V/c^2 "
                  "%.6e\n", pmy_pack->pmesh->ncycle, t2_solve, n, p);
    }
    if (n > 0.0 && cln_fclz) {
      Real q[3];   // (every rank: the sums are collective)
      for (int d = 0; d < 3; ++d) q[d] = CleanReduce(h(2+d), false);
      if (global_variable::my_rank == 0) {
        std::printf("CLEAN FCONE_LZ cycle %d t2s %d momentum removed (signed, V/c^2) x1 "
                    "%.6e x2 %.6e x3 %.6e\n", pmy_pack->pmesh->ncycle, t2_solve, q[0],
                    q[1], q[2]);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRingMatesInit
//! \brief the MeshBlocks of each local block's line along d, position by position
//! (CleanStagFwbProject): src(m, p) >= 0 the local block, else -1 - slot of a received
//! one; snd / rcv the blocks exchanged (block or slot, rank, gid), one message each

void RadiationVet::CleanRingMatesInit(const int d, CleanRingMates &rm) {
  Mesh *pm = pmy_pack->pmesh;
  const int nr[3] = {pm->nmb_rootx1, pm->nmb_rootx2, pm->nmb_rootx3};
  const int nmb = pmy_pack->nmb_thispack, me = global_variable::my_rank;
  rm.nmb = nmb;
  rm.nt = pm->nmb_total;
  rm.snd.clear();
  rm.rcv.clear();
  std::map<int, int> slot;   // (received gid -> slot)
  std::set<std::array<int, 3>> sent;
  Kokkos::realloc(rm.src, nmb, nr[d]);
  auto h = Kokkos::create_mirror_view(rm.src);
  auto loc = [&](const int g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    return std::array<int, 4>{static_cast<int>(l.lx1), static_cast<int>(l.lx2),
                              static_cast<int>(l.lx3), l.level};
  };
  for (int m = 0; m < nmb; ++m) {
    const int gm = pmy_pack->pmb->mb_gid.h_view(m);
    std::array<int, 4> w = loc(gm);
    const int p0 = w[d];
    for (int p = 0; p < nr[d]; ++p) {
      w[d] = p;
      int g = 0;
      while (g < pm->nmb_total && loc(g) != w) ++g;
      if (g == pm->nmb_total) CleanFatal("fwb: no MeshBlock at a line's position");
      const int rk = pm->rank_eachmb[g];
      if (rk == me) {
        h(m,p) = g - pm->gids_eachrank[me];
        continue;
      }
      if (slot.find(g) == slot.end()) {
        const int k = static_cast<int>(slot.size());
        slot[g] = k;
        rm.rcv.push_back({k, rk, g});
      }
      h(m,p) = -1 - slot[g];
      if (p != p0) sent.insert({m, rk, gm});
    }
  }
  rm.snd.assign(sent.begin(), sent.end());
  rm.nslot = static_cast<int>(slot.size());
  Kokkos::deep_copy(rm.src, h);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFwbProject
//! \brief clean_stag_fwb_line (d = 0), clean_stag_fwb_ring (d = 2): fwb's residual s of
//! the faces along d (rs, its weight w/(1 + tau) rw, CleanStagFSolve) less mu rw,
//! mu = sum V_f s/sum V_f rw over the whole line through the mesh (a radial line, a phi
//! ring), added to face slot slot.  A radial line's faces share their direction, a
//! ring's their lever arm R sin(theta): the net s V of each line is 0, so fwb gives no
//! momentum along it and no L_z.  A line across MeshBlocks: each block's sums at its
//! position along d in one mesh-wide array, summed over the ranks (every slot from one
//! block, the others 0: exact), then over the positions in order, so mu does not depend
//! on the rank count (one MeshBlock along d: the block's own sums).  Tile mode: the
//! sums of each tile along a lateral d, then the tiles in order, so mu does not depend on
//! the MeshBlocks either

void RadiationVet::CleanStagFwbProject(const int d, const int slot,
                                       const DvceArray4D<Real> &rs,
                                       const DvceArray4D<Real> &rw, const bool cnt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, nmb = pmy_pack->nmb_thispack;
  Mesh *pm = pmy_pack->pmesh;
  const int nx[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nr[3] = {pm->nmb_rootx1, pm->nmb_rootx2, pm->nmb_rootx3};
  const int da = (d == 0) ? 1 : 0, db = (d == 2) ? 1 : 2;   // transverse dirs
  const int n = nx[d], na = nx[da], nb = nx[db];
  const int dk = (d == 2), dj = (d == 1), di = (d == 0);
  const bool wrap = nr[d] == 1 && pm->mesh_bcs[2*d] == BoundaryFlag::periodic;
  CleanArgs x = CleanMakeArgs();
  auto fc = x.fc[d];
  auto fwc = cln_fwbc;
  const int tn = (d > 0 && repro::tile_n > 0) ? repro::tile_n : n, nc = n/tn;
  auto mu = ScratchView<DvceArray4D<Real>>("cln_fwmu" + std::to_string(d), nmb, nb, na,
                                           2*nc);
  auto face = KOKKOS_LAMBDA(const int t, const int a, const int b, int &k, int &j,
                            int &i) {
    i = is + ((d == 0) ? t : ((da == 0) ? a : b));
    j = js + ((d == 1) ? t : ((da == 1) ? a : b));
    k = ks + ((d == 2) ? t : b);
  };
  par_for("cln_fwbp", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1,
  KOKKOS_LAMBDA(const int m, const int b, const int a) {
    for (int c = 0; c < nc; ++c) {
      Real ss = 0.0, sw = 0.0;
      for (int t = c*tn; t < (c + 1)*tn; ++t) {   // (V_f: the face's weight, own faces)
        int k, j, i;
        face(t, a, b, k, j, i);
        const Real vf = CleanVol(x, m, k, j, i) + CleanVol(x, m, k-dk, j-dj, i-di);
        ss += vf*rs(m,k,j,i);
        sw += vf*rw(m,k,j,i);
      }
      mu(m,b,a,2*c) = ss;
      mu(m,b,a,2*c+1) = sw;
    }
  });
  int nq = nc;
  auto mq = mu;
  if (nr[d] > 1) {   // (the other MeshBlocks' sums, then all positions in order)
    CleanRingMates &rm = cln_fwr[d];
    if (rm.nmb != nmb || rm.nt != pm->nmb_total) CleanRingMatesInit(d, rm);
    auto rb = ScratchView<DvceArray4D<Real>>("cln_fwrb" + std::to_string(d),
                                             std::max(rm.nslot, 1), nb, na, 2*nc);
    mq = ScratchView<DvceArray4D<Real>>("cln_fwmq" + std::to_string(d), nmb, nb, na, 2);
#if MPI_PARALLEL_ENABLED
    const int nw = nb*na*2*nc;
    std::vector<MPI_Request> rq(rm.rcv.size() + rm.snd.size());
    DevExeSpace().fence();
    size_t r = 0;
    for (auto &v : rm.rcv) {
      MPI_Irecv(rb.data() + static_cast<size_t>(v[0])*nw, nw, MPI_ATHENA_REAL, v[1],
                4250 + v[2], MPI_COMM_WORLD, &rq[r++]);
    }
    for (auto &v : rm.snd) {
      MPI_Isend(mu.data() + static_cast<size_t>(v[0])*nw, nw, MPI_ATHENA_REAL, v[1],
                4250 + v[2], MPI_COMM_WORLD, &rq[r++]);
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
#endif
    const int np = nr[d], ncc = nc;
    auto src = rm.src;
    auto mo = mq;
    nq = 1;
    par_for("cln_fwbr", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1, 0, 1,
    KOKKOS_LAMBDA(const int m, const int b, const int a, const int c) {
      Real s = 0.0;
      for (int p = 0; p < np; ++p) {
        const int sl = src(m,p);
        for (int t = 0; t < ncc; ++t) {
          s += (sl >= 0) ? mu(sl,b,a,2*t+c) : rb(-1-sl,b,a,2*t+c);
        }
      }
      mo(m,b,a,c) = s;
    });
  }
  par_for("cln_fwbq", DevExeSpace(), 0, nmb - 1, 0, nb - 1, 0, na - 1,
  KOKKOS_LAMBDA(const int m, const int b, const int a) {
    Real ss = 0.0, sw = 0.0;
    for (int c = 0; c < nq; ++c) {
      ss += mq(m,b,a,2*c);
      sw += mq(m,b,a,2*c+1);
    }
    const Real q = (sw > 0.0) ? ss/sw : 0.0;
    Real cn = 0.0, cp = 0.0, c1 = 0.0;   // (the line's counters: one atomic each)
    for (int t = 0; t <= n; ++t) {   // (one MeshBlock round a periodic wrap: low face)
      int k, j, i, kr, jr, ir;
      face(t, a, b, k, j, i);
      face((wrap && t == n) ? 0 : t, a, b, kr, jr, ir);
      const Real s = rs(m,kr,jr,ir) - q*rw(m,kr,jr,ir);
      fc(m,slot,k,j,i) += s/(x.h*x.c*x.c);
      if (!cnt || t == n) continue;   // counters: the printed s, the write-back's
      const Real p = fc(m,CLF_TH,k,j,i)*s*0.5*(CleanVol(x, m, k, j, i) +
                     CleanVol(x, m, k-dk, j-dj, i-di))/(x.c*x.c);
      cn += 1.0;
      cp += fabs(p);
      c1 += p;
    }
    if (!cnt) return;
    Kokkos::atomic_add(&fwc(0), cn);
    Kokkos::atomic_add(&fwc(1), cp);
    if (d == 0) Kokkos::atomic_add(&fwc(2), c1);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFaceTie
//! \brief the face tie (spherical polar, Cartesian): a face shared by two MeshBlocks (or
//! one across a periodic wrap) takes its owner's values (the block whose low face it is),
//! as its F does: the momentum exchange T = (F - G)/c^2 + h div P - h N1
//! (CleanMomentumChange), the memory d and its slope.  Both sides then give the gas the
//! halves of one T.  Each block computes its low faces' into CLW_L + 3 dir (+ 1 d, + 2
//! slope); the exchange hands them to the other side.  Cubed sphere (and
//! clean_stag_tie_gid): T in CLF_PT, the lower gid's (CleanStagCsTie).

void RadiationVet::CleanStagFaceTie(const int e, const int kq, const bool pt) {
  // d's slope: clean_stag_fkin only (otherwise the face arrays end before CLF_KD)
  const int kd = cln_fkin ? kq : 0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.fs = true;
  const int nf = x.nf;
  if (x.cub || cln_tgid) {
    for (int d = 0; d < nf && !pt; ++d) {   // as CleanStagDrag; none at physical faces
      auto fc = x.fc[d];
      CleanFor("cln_tiet", nmb1, ks, ke+(d == 2), js, je+(d == 1), is, ie+(d == 0),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (CleanBoundaryFace(x, d, m, k, j, i)) return;
        fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) +
            x.h*CleanFaceStress(x, e, d, m, k, j, i, true) - x.h*fc(m,CLF_N+1,k,j,i);
      });
    }
    CleanStagCsTie((kd > 0) ? 3 : 2, CLF_PT, CLF_N + 3, kd);
    return;
  }
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
  CleanFor("cln_tie", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    CLEAN_HIPX(const CleanArgs &x = *px;)   // in device memory: no per-thread copy
    for (int d = 0; d < nf; ++d) {
      const Real f = x.fc[d](m,CLF_N,k,j,i), g = x.fc[d](m,CLF_G,k,j,i);
      // pt: CleanStagDrag just took it from the same F (CLF_PT), except on physical faces
      x.w(m,CLW_L+3*d,k,j,i) = (pt && !CleanBoundaryFace(x, d, m, k, j, i)) ?
          x.fc[d](m,CLF_PT,k,j,i) : (f - g)/(x.c*x.c) + x.h*CleanFaceStress(x, e, d, m, k,
                               j, i, true) - x.h*x.fc[d](m,CLF_N+1,k,j,i);
      x.w(m,CLW_L+3*d+1,k,j,i) = x.fc[d](m,CLF_N+3,k,j,i);
      x.w(m,CLW_L+3*d+2,k,j,i) = (kd > 0) ? x.fc[d](m,kd,k,j,i) : 0.0;
    }
  });
  CleanGhostFill(cln_w, CLW_L, 3*nf);
  for (int d = 0; d < nf; ++d) {   // the block's high faces: their owner's d
    const int dk = (d == 2), dj = (d == 1), di = (d == 0);
    auto fc = x.fc[d];
    par_for("cln_tie_hi", DevExeSpace(), 0, nmb1, dk ? ke+1 : ks, ke+dk, dj ? je+1 : js,
            je+dj, di ? ie+1 : is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, d, m, k, j, i)) return;
      fc(m,CLF_N+3,k,j,i) = x.w(m,CLW_L+3*d+1,k,j,i);
      if (kd > 0) fc(m,kd,k,j,i) = x.w(m,CLW_L+3*d+2,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanStagFSweepChange
//! \brief clean_stag_fsweeps: the largest relative change of the F the lines' high faces
//! took from their neighbours (CLW_Y, CLW_Z, CLW_S ghosts) since the last sweep (kept in
//! CLW_TM + dir; the first sweep only stores it)

Real RadiationVet::CleanStagFSweepChange() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack, nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nx3 = indcs.nx3, n12 = nx1*nx2, n13 = nx1*nx3, n23 = nx2*nx3;
  const int nfc = n23 + n13 + n12;
  CleanArgs x = CleanMakeArgs();
  const int nf = x.nf;
  Real dmax = 0.0;
  Kokkos::parallel_reduce("cln_fsweep", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nfc),
  KOKKOS_LAMBDA(const int idx, Real &mx) {
    const int m = idx/nfc, q = idx%nfc;
    const int d = (q < n23) ? 0 : ((q < n23 + n13) ? 1 : 2);
    if (d >= nf) return;
    int k, j, i;
    if (d == 0) {i = ie + 1; j = js + q%nx2; k = ks + q/nx2;}
    if (d == 1) {j = je + 1; i = is + (q - n23)%nx1; k = ks + (q - n23)/nx1;}
    if (d == 2) {k = ke + 1; i = is + (q - n23 - n13)%nx1; j = js + (q - n23 - n13)/nx1;}
    if (CleanBoundaryFace(x, d, m, k, j, i)) return;
    const int cy = (d == 2) ? CLW_S : CLW_Y + d;
    const Real f = x.w(m,cy,k,j,i), f0 = x.w(m,CLW_TM+d,k,j,i);
    mx = fmax(mx, fabs(f - f0)/fmax(fabs(f), 1.0e-300));
    x.w(m,CLW_TM+d,k,j,i) = f;
  }, Kokkos::Max<Real>(dmax));
  return CleanReduce(dmax, true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFSchur
//! \brief clean_stag_fschur (uniform mesh; with the face tie or the cubed sphere's owner
//! rule): the F solve's lines exact across MeshBlocks and periodic wraps (cubed sphere:
//! round the four-panel rings).  Each MeshBlock line is a segment, its end faces affine
//! in the faces beyond (CleanStagFOwn).  Exchange rounds join segments (CleanStagFJoin):
//! a MeshBlock sends right its own joined to the one from the left (and mirrored) until
//! they span the line, or once round a periodic one; the faces beyond its ends then
//! replace the lagged ghosts in the line solves.  A segment goes out in three slots of
//! three cells at each end; the six free slots serve two directions at a time.  Across
//! a cubed-sphere seam a reversed neighbour line (CleanCsSlot) is read reversed.

void RadiationVet::CleanStagFSchur() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const CleanFLine ln{cln_v, x.cub};
  const int nn[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nm[3] = {pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3};
  const int fs[6] = {CLW_TM, CLW_TM+1, CLW_TM+2, CLW_TM+3, CLW_BW, CLW_PW};
  // rounds to the left and right: a line of nb MeshBlocks nb - 1 each, a periodic one
  // nb/2 and (nb - 1)/2; the group of directions exchanged together, and their slots
  // (cubed sphere: x2 and x3 lines in one group, as a seam may turn one into the other)
  const bool cub = x.cub;
  int nr[3][2], gr[3], gs[3] = {0, 0, 0}, ngd = 0;
  bool per[3];
  for (int dd = 0; dd < x.nf; ++dd) {
    const int d = cub ? (dd + 1)%3 : dd;
    const int nb = (cub && d > 0) ? 4*nm[d]/nn[d] : nm[d]/nn[d];
    per[d] = (cub && d > 0) || (pm->mesh_bcs[2*d] == BoundaryFlag::periodic);
    nr[d][0] = per[d] ? nb/2 : nb - 1;
    nr[d][1] = per[d] ? (nb - 1)/2 : nb - 1;
    gr[d] = (nr[d][0] > 0) ? ngd/2 : ((per[d]) ? 0 : -1);
    gs[d] = (nr[d][0] > 0) ? 3*(ngd++%2) : 0;
  }
  // layout_invariant_mode = tile, a lateral line of more than one tile
  const int tn = repro::tile_n;
  if (tn > 0 && (nn[1] > tn || (x.nf == 3 && nn[2] > tn))) {
    CleanStagFSchurTile(gr, gs, per, ngd);
    return;
  }
  for (int g = 0; g == 0 || 2*g < ngd; ++g) {
    int nrt = 0, nd = 0;
    for (int d = 0; d < x.nf; ++d) {
      if (gr[d] == g && nr[d][0] > 0) nrt = std::max(nrt, nr[d][0]), ++nd;
    }
    for (int r = 1; r <= nrt + 1; ++r) {
      const bool fin = (r == nrt + 1);
      for (int dir = 0; dir < x.nf; ++dir) {
        const int rl = nr[dir][0], rr = nr[dir][1];
        if (gr[dir] != g || (!fin && (rl == 0 || (r > rl && r != rr + 1)))) continue;
        const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
        const int nb = (dir == 2) ? nn[1] : nn[2];
        auto fc = x.fc[dir];
        const int c0 = fs[gs[dir]], c1 = fs[gs[dir]+1], c2 = fs[gs[dir]+2];
        const int g1 = gs[1], g2 = gs[2];
        const bool pd = per[dir];
        const int rs = cub ? rl : rr;   // rounds the low end sends
        AthFor("cln_fsch", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             (nmb1 + 1)*na*nb),
        KOKKOS_LAMBDA(const int idx) {
          const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
          auto cell = [&](const int t, int &k, int &j, int &i) {
            k = (dir == 2) ? x.ks + t : x.ks + b;
            j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
            i = (dir == 0) ? x.is + t : x.is + a;
          };
          // number q of a segment: slot q/3, cell t0 + sg (q%3)
          auto slot = [&](const int t0, const int sg, const int q) -> Real & {
            int k, j, i;
            cell(t0 + sg*(q%3), k, j, i);
            return x.w(m,(q < 3) ? c0 : ((q < 6) ? c1 : c2),k,j,i);
          };
          auto keepo = [&](const int q) -> Real & {
            int k, j, i;
            cell(q/4, k, j, i);
            return x.w(m,((gs[dir] == 0) ? CLW_AM : CLW_R) + q%4,k,j,i);
          };
          // a segment in the ghosts beyond the low or high end; cubed sphere: in the
          // neighbour family's slots, reversed (F -> -F) where it runs the other way
          auto ghost = [&](const bool he, Real g[8]) {
            int k, j, i, kc, jc, ic;
            cell(he ? n : -1, k, j, i);
            cell(he ? n : 0, kc, jc, ic);
            Real sg = 1.0;
            const int fq = cub ? CleanCsSlot(x, dir, m, kc, jc, ic, k, j, i, false, sg)/2
                               : dir;
            const int o = (fq == dir) ? gs[dir] : ((fq == 1) ? g1 : g2);
            for (int q = 0; q < 8; ++q) {
              cell(he ? n + q%3 : -1 - q%3, k, j, i);
              g[q] = x.w(m,fs[o + q/3],k,j,i);
            }
            if (sg > 0.0) return;
            const Real h[8] = {g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]};
            g[0] = -h[4], g[1] = h[6], g[2] = h[5], g[3] = h[7];
            g[4] = -h[0], g[5] = h[2], g[6] = h[1], g[7] = h[3];
          };
          int k, j, i;
          cell(0, k, j, i);
          const bool bl = CleanBoundaryFace(x, dir, m, k, j, i);
          cell(n, k, j, i);
          const bool bh = CleanBoundaryFace(x, dir, m, k, j, i);
          const bool hl = r > 1 && rl > 0 && !bl, hr = r > 1 && rr > 0 && !bh;
          Real o[8], sl[8], sr[8], c[8], cr[8], ya[3], xb[3], yr[3], xr[3];
          // the own segment once, kept for the rounds in cells 0, 1 of four E-solve
          // scratch slots, a set per direction of the group (a one-MeshBlock periodic
          // line sends nothing: only in the last round)
          const bool r1 = r == 1 || rl == 0;
          for (int q = 0; q < 8 && !r1; ++q) o[q] = keepo(q);
          if (r1) CleanStagFOwn(x, ln, fc, dir, m, a, b, n, o);
          for (int q = 0; q < 8 && r1 && !fin; ++q) keepo(q) = o[q];
          // the segment to the left in the ghosts; to the right kept after rr (cubed
          // sphere: in the faces' free CLF_N + 2, both ends send for rl rounds)
          auto keep = [&](const int q) -> Real & {
            int k, j, i;
            cell(q, k, j, i);
            return fc(m,CLF_N+2,k,j,i);
          };
          if (hl) ghost(false, sl);
          if (hr && (r - 1 <= rr || cub)) ghost(true, sr);
          for (int q = 0; q < 8; ++q) {
            if (!hl) sl[q] = 0.0;
            if (hr && r - 1 > rr && !cub) sr[q] = slot(0, 1, q);
            if (!hr) sr[q] = 0.0;
            c[q] = o[q], cr[q] = o[q];
          }
          if (hl) CleanStagFJoin(sl, o, c, ya, xb);
          if (!fin) {
            if (hr && r <= rs) CleanStagFJoin(o, sr, cr, yr, xr);
            for (int q = 0; q < 8; ++q) {
              if (r <= rl) slot(n - 1, -1, q) = c[q];
              if (r <= rs + ((hr && !cub) ? 1 : 0)) {
                slot(0, 1, q) = (r <= rs) ? cr[q] : sr[q];
              }
              if (cub && hr && r == rr + 1) keep(q) = sr[q];
            }
            return;
          }
          for (int q = 0; q < 8 && cub && hr && r - 1 > rr; ++q) sr[q] = keep(q);
          // the whole line (periodic: determinant (1 - c2)(1 - c5) - c1 c6 as a sum)
          for (int q = 0; q < 8; ++q) cr[q] = c[q];
          if (hr) CleanStagFJoin(c, sr, cr, yr, xr);
          Real xo = 0.0, yo = 0.0;
          if (pd) {
            const Real det = cr[3]*(cr[6] + cr[7]) + cr[1]*cr[7];
            xo = (cr[0]*(cr[6] + cr[7]) + cr[1]*cr[4])/det;
            yo = (cr[4]*(cr[1] + cr[3]) + cr[6]*cr[0])/det;
          }
          const Real gr = hr ? xr[0] + xr[1]*yo + xr[2]*xo : xo;
          const Real gl = hl ? ya[0] + ya[1]*yo + ya[2]*gr : yo;
          // the faces beyond the own ones (next block's first xR, last block's last yL),
          // in the slots the line solve reads
          Real sg;
          int kf, jf, iff;
          if (!bh) {   // cell n
            const bool on = ln.Own(x, dir, m, k, j, i, n, n, true);
            Real &f = ln.Slot(x, dir, m, k, j, i, k, j, i, on, sg);
            f = sg*gr;
          }
          cell(0, kf, jf, iff);
          cell(-1, k, j, i);
          if (!bl) {
            const bool on = ln.Own(x, dir, m, kf, jf, iff, 0, n, true);
            Real &f = ln.Slot(x, dir, m, kf, jf, iff, k, j, i, on, sg);
            f = sg*gl;
          }
        });
      }
      // the ghosts of a face take its line family's slots (CleanGSeg::fm)
      if (r <= nrt) {
        std::array<uint32_t, 4> fa = {0, 0, 0, 0}, fb = fa;
        for (int d = 0; d < x.nf; ++d) {
          if (gr[d] == g) fa[d] = (gs[d] > 0) ? 3 : 0, fb[d] = (gs[d] > 0) ? 8 : 7;
        }
        if (nd == 2) CleanGhostDefer(cln_w, CLW_BW, 2), cln_hpend.back().fm = fa;
        CleanGhostDefer(cln_w, CLW_TM, (nd == 2) ? 4 : 3);
        cln_hpend.back().fm = fb;
        CleanGhostFlush();
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFSchurTile
//! \brief layout_invariant_mode = tile: CleanStagFSchur with the nk tiles of a lateral
//! MeshBlock line as MeshBlocks, so the F solve is the same for every whole-tile layout.
//! Tile g's whole line is join(join(L, o_g), R), L the left fold of the nt/2 segments
//! before it (periodic; else all), R the right fold of the (nt - 1)/2 after.  Joins do
//! not associate in round-off, so the folds are built in that order: in pass p each
//! MeshBlock sends right the left fold from its tile p and left the right fold to its
//! tile nk - 1 - p, extended by each block they pass.  Inner tile ends go to cln_fsjt,
//! the MeshBlock's to the ghosts.

void RadiationVet::CleanStagFSchurTile(const int gr[3], const int gs[3],
                                       const bool per[3], const int ngd) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const CleanFLine ln{cln_v, x.cub};
  const bool cub = x.cub;
  const int nn[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nm[3] = {pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3};
  const int fs[6] = {CLW_TM, CLW_TM+1, CLW_TM+2, CLW_TM+3, CLW_BW, CLW_PW};
  const int g1 = gs[1], g2 = gs[2];
  int nk[3] = {1, 1, 1}, nbk[3] = {1, 1, 1};
  for (int d = 0; d < x.nf; ++d) {
    if (d > 0) nk[d] = nn[d]/repro::tile_n;
    nbk[d] = ((cub && d > 0) ? 4 : 1)*nm[d]/nn[d];
  }
  // the tiles' segments
  for (int dir = 0; dir < x.nf; ++dir) {
    if (nk[dir] == 1 && gr[dir] < 0) continue;
    const int nkd = nk[dir], nt = nn[dir]/nkd, na = (dir == 0) ? nn[1] : nn[0];
    const int nb = (dir == 2) ? nn[1] : nn[2];
    CleanFsjt(cln_fsjt[dir], nmb1 + 1, na*nb, nkd);
    auto jt = cln_fsjt[dir];
    auto fc = x.fc[dir];
    AthFor("cln_fsjo", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
      for (int u = 0; u < nkd; ++u) {
        Real s[8];
        CleanStagFOwn(x, ln, fc, dir, m, a, b, nt, s, u*nt, u > 0, u < nkd - 1);
        for (int q = 0; q < 8; ++q) jt(m,l,8*u+q) = s[q];
      }
    });
  }
  // the passes of rounds of each group of directions (CleanStagFSchur)
  for (int g = 0; g == 0 || 2*g < ngd; ++g) {
    int np = 0, nrt = 0, nd = 0;
    for (int d = 0; d < x.nf; ++d) {
      if (gr[d] != g || nbk[d] < 2) continue;
      np = std::max(np, per[d] ? nk[d] : 1);
      nrt = std::max(nrt, per[d] ? nbk[d]/2 : nbk[d] - 1);
      ++nd;
    }
    for (int p = 0; p < np; ++p) {
      for (int r = 1; r <= nrt + 1; ++r) {
        for (int dir = 0; dir < x.nf; ++dir) {
          const bool pd = per[dir];
          const int nkd = nk[dir], ns = pd ? nbk[dir]/2 : nbk[dir] - 1;
          if (gr[dir] != g || nbk[dir] < 2 || p >= (pd ? nkd : 1) || r > ns + 1) continue;
          const int ul = pd ? p : 0, ur = pd ? (nkd - p)%nkd : 0;
          // the round the right fold is kept (tile 0 of a periodic line has (nt - 1)/2
          // tiles after it)
          const int rk = (pd && p == 0) ? ((nbk[dir]*nkd - 1)/2)/nkd + 1 : ns + 1;
          const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
          const int nb = (dir == 2) ? nn[1] : nn[2];
          const int c0 = fs[gs[dir]], c1 = fs[gs[dir]+1], c2 = fs[gs[dir]+2];
          const int gd = gs[dir];
          auto jt = cln_fsjt[dir];
          AthFor("cln_fsjr", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
          KOKKOS_LAMBDA(const int idx) {
            const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
            auto cell = [&](const int t, int &k, int &j, int &i) {
              k = (dir == 2) ? x.ks + t : x.ks + b;
              j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
              i = (dir == 0) ? x.is + t : x.is + a;
            };
            auto slot = [&](const int t0, const int sg, const int q) -> Real & {
              int k, j, i;
              cell(t0 + sg*(q%3), k, j, i);
              return x.w(m,(q < 3) ? c0 : ((q < 6) ? c1 : c2),k,j,i);
            };
            auto ghost = [&](const bool he, Real g[8]) {   // as CleanStagFSchur
              int k, j, i, kc, jc, ic;
              cell(he ? n : -1, k, j, i);
              cell(he ? n : 0, kc, jc, ic);
              Real sg = 1.0;
              const int fq = cub ? CleanCsSlot(x, dir, m, kc, jc, ic, k, j, i, false,
                                               sg)/2 : dir;
              const int o = (fq == dir) ? gd : ((fq == 1) ? g1 : g2);
              for (int q = 0; q < 8; ++q) {
                cell(he ? n + q%3 : -1 - q%3, k, j, i);
                g[q] = x.w(m,fs[o + q/3],k,j,i);
              }
              if (sg > 0.0) return;
              const Real h[8] = {g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]};
              g[0] = -h[4], g[1] = h[6], g[2] = h[5], g[3] = h[7];
              g[4] = -h[0], g[5] = h[2], g[6] = h[1], g[7] = h[3];
            };
            int k, j, i;
            cell(0, k, j, i);
            const bool xl = r > 1 && !CleanBoundaryFace(x, dir, m, k, j, i);
            cell(n, k, j, i);
            const bool xr = r > 1 && !CleanBoundaryFace(x, dir, m, k, j, i);
            Real sl[8], sr[8], f[8], h[8], s[8], ya[3], xb[3];
            if (xl) ghost(false, sl);
            if (xr) ghost(true, sr);
            for (int q = 0; q < 8; ++q) {
              if (xl && r == ns + 1) jt(m,l,8*nkd+8*ul+q) = sl[q];
              if (xr && r == rk) jt(m,l,16*nkd+8*ur+q) = sr[q];
            }
            if (r > ns) return;
            // to the right: the left fold from tile ul (or the one from the left
            // extended), to the left: the right fold to tile nk - 1 - ul (or extended)
            const int u0 = xl ? 0 : ul, u1 = xr ? nkd - 1 : nkd - 1 - ul;
            for (int q = 0; q < 8; ++q) f[q] = xl ? sl[q] : jt(m,l,8*u0+q);
            for (int u = xl ? 0 : u0 + 1; u < nkd; ++u) {
              for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q);
              CleanStagFJoin(f, s, h, ya, xb);
              for (int q = 0; q < 8; ++q) f[q] = h[q];
            }
            for (int q = 0; q < 8; ++q) slot(n - 1, -1, q) = f[q];
            for (int q = 0; q < 8; ++q) f[q] = xr ? sr[q] : jt(m,l,8*u1+q);
            for (int u = xr ? nkd - 1 : u1 - 1; u >= 0; --u) {
              for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q);
              CleanStagFJoin(s, f, h, ya, xb);
              for (int q = 0; q < 8; ++q) f[q] = h[q];
            }
            for (int q = 0; q < 8; ++q) slot(0, 1, q) = f[q];
          });
        }
        if (r <= nrt) {   // as CleanStagFSchur
          std::array<uint32_t, 4> fa = {0, 0, 0, 0}, fb = fa;
          for (int d = 0; d < x.nf; ++d) {
            if (gr[d] == g) fa[d] = (gs[d] > 0) ? 3 : 0, fb[d] = (gs[d] > 0) ? 8 : 7;
          }
          if (nd == 2) CleanGhostDefer(cln_w, CLW_BW, 2), cln_hpend.back().fm = fa;
          CleanGhostDefer(cln_w, CLW_TM, (nd == 2) ? 4 : 3);
          cln_hpend.back().fm = fb;
          CleanGhostFlush();
        }
      }
    }
  }
  // each tile's whole line, its ends
  for (int dir = 0; dir < x.nf; ++dir) {
    if (nk[dir] == 1 && gr[dir] < 0) continue;
    const bool pd = per[dir];
    const int nkd = nk[dir], n = nn[dir], nbd = nbk[dir], na = (dir == 0) ? nn[1] : nn[0];
    const int nb = (dir == 2) ? nn[1] : nn[2];
    auto jt = cln_fsjt[dir];
    AthFor("cln_fsjf", Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*na*nb),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb, l = idx%(na*nb);
      auto cell = [&](const int t, int &k, int &j, int &i) {
        k = (dir == 2) ? x.ks + t : x.ks + b;
        j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
        i = (dir == 0) ? x.is + t : x.is + a;
      };
      int k, j, i, kf, jf, iff;
      cell(0, k, j, i);
      const bool bl = CleanBoundaryFace(x, dir, m, k, j, i);
      cell(n, k, j, i);
      const bool bh = CleanBoundaryFace(x, dir, m, k, j, i);
      for (int u = 0; u < nkd; ++u) {
        Real lf[8], rf[8], s[8], h[8], c[8], cr[8], ya[3], xb[3], yr[3], xr[3];
        Real yd[3], xd[3];
        bool hl = false, hr = false;
        auto addl = [&](const int o) {   // L = join(L, the segment at o)
          for (int q = 0; q < 8; ++q) s[q] = jt(m,l,o+q);
          if (hl) CleanStagFJoin(lf, s, h, yd, xd);
          for (int q = 0; q < 8; ++q) lf[q] = hl ? h[q] : s[q];
          hl = true;
        };
        auto addr = [&](const int o) {   // R = join(the segment at o, R)
          for (int q = 0; q < 8; ++q) s[q] = jt(m,l,o+q);
          if (hr) CleanStagFJoin(s, rf, h, yd, xd);
          for (int q = 0; q < 8; ++q) rf[q] = hr ? h[q] : s[q];
          hr = true;
        };
        if (pd && nbd == 1) {   // round the wrap of one MeshBlock
          for (int t = u - nkd/2; t < u; ++t) addl(8*((t + nkd)%nkd));
          for (int t = u + (nkd - 1)/2; t > u; --t) addr(8*(t%nkd));
        } else {
          if (pd || !bl) addl(8*nkd + 8*(pd ? u : 0));
          for (int t = 0; t < u; ++t) addl(8*t);
          if (pd ? (u > 0 || (nbd*nkd - 1)/2 >= nkd) : !bh) addr(16*nkd + 8*(pd ? u : 0));
          for (int t = nkd - 1; t > u; --t) addr(8*t);
        }
        for (int q = 0; q < 8; ++q) s[q] = jt(m,l,8*u+q), c[q] = s[q];
        if (hl) CleanStagFJoin(lf, s, c, ya, xb);
        for (int q = 0; q < 8; ++q) cr[q] = c[q];
        if (hr) CleanStagFJoin(c, rf, cr, yr, xr);
        Real xo = 0.0, yo = 0.0;
        if (pd) {
          const Real det = cr[3]*(cr[6] + cr[7]) + cr[1]*cr[7];
          xo = (cr[0]*(cr[6] + cr[7]) + cr[1]*cr[4])/det;
          yo = (cr[4]*(cr[1] + cr[3]) + cr[6]*cr[0])/det;
        }
        const Real gR = hr ? xr[0] + xr[1]*yo + xr[2]*xo : xo;
        const Real gL = hl ? ya[0] + ya[1]*yo + ya[2]*gR : yo;
        if (u > 0) jt(m,l,24*nkd+u) = gL;
        if (u < nkd - 1) jt(m,l,25*nkd+u) = gR;
        Real sg;
        cell(n, k, j, i);
        if (u == nkd - 1 && !bh) {   // cell n, as CleanStagFSchur
          const bool on = ln.Own(x, dir, m, k, j, i, n, n, true);
          Real &f = ln.Slot(x, dir, m, k, j, i, k, j, i, on, sg);
          f = sg*gR;
        }
        cell(0, kf, jf, iff);
        cell(-1, k, j, i);
        if (u == 0 && !bl) {
          const bool on = ln.Own(x, dir, m, kf, jf, iff, 0, n, true);
          Real &f = ln.Slot(x, dir, m, kf, jf, iff, k, j, i, on, sg);
          f = sg*gL;
        }
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagFCsShared
//! \brief cubed sphere, the F solve's last step: a face shared with another MeshBlock
//! takes the lower gid's solution (in our axes, CleanCsSlot), bitwise the same on both
//! sides (CLF_N + 2), and the cells' slots of cln_v take it

void RadiationVet::CleanStagFCsShared() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto v = cln_v;
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor("cln_fcs3", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      if (CleanBoundaryFace(x, dir, m, k, j, i)) {
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);
        return;
      }
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      Real f = lin ? v(m,9+2*dir+1,kl,jl,il) : v(m,9+2*dir,k,j,i);
      if (!(lin && rin)) {   // the other MeshBlock's solution, in our axes
        const int kg = lin ? k : kl, jg = lin ? j : jl, ig = lin ? i : il;
        const int ko = lin ? kl : k, jo = lin ? jl : j, io = lin ? il : i;
        Real s;
        const int n = CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s);
        if (x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io)) f = s*v(m,9+n,kg,jg,ig);
      }
      fc(m,CLF_N+2,k,j,i) = f;
    });
  }
  CleanFor("cln_fcs4", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int q = 0; q < 3; ++q) {
      v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
      v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
    }
  });
  // one MeshBlock along x1: not the x1 family's F (CleanStagFSolve)
  const int f0 =
      (pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel) ? 2 : 0;
  CleanGhostFill(cln_v, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));   // (as the F solve)
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagCsTie
//! \brief cubed sphere, the face tie: the nc face slots c0, c1, c2 of a face shared with
//! another MeshBlock take the lower gid's values in our axes (CleanCsSlot), bitwise one
//! value on both sides.  Each cell puts its boundary faces' in cln_v slot 3 n + q, which
//! the ghosts take.  clean_stag_tie_gid on spherical polar / Cartesian: the same in
//! CLW_L + 3 n + q; a periodic self-neighbour takes its low face's (the owner's).

void RadiationVet::CleanStagCsTie(const int nc, const int c0, const int c1,
                                  const int c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const bool cub = x.cub;
  const int v0 = cub ? 0 : CLW_L, nf = x.nf;
  const CleanW v = cub ? CleanW(cln_v) : cln_w;
  auto gid = pmy_pack->pmb->mb_gid.d_view;
  auto nb = pmy_pack->pmb->nghbr.d_view;
  CleanFor("cln_cstie", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int c[3] = {c0, c1, c2};
    for (int q = 0; q < nf; ++q) {
      const int u = (q == 0) ? (i == ie) : ((q == 1) ? (j == je) : (k == ke));
      for (int n = 0; n < nc; ++n) {
        v(m,v0+3*n+q,k,j,i) = x.fc[q](m,c[n],k+u*(q == 2),j+u*(q == 1),i+u*(q == 0));
      }
    }
  });
  if (cub) {
    CleanGhostFill(cln_v, 0, 3*nc, true, 0, 1);
  } else {
    CleanGhostFill(cln_w, CLW_L, 3*nc, 0, 1);
  }
  for (int dir = 0; dir < nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor("cln_cstie2", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      if ((lin && rin) || CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const int kg = lin ? k : k - dk, jg = lin ? j : j - dj, ig = lin ? i : i - di;
      const int ko = lin ? k - dk : k, jo = lin ? j - dj : j, io = lin ? i - di : i;
      // sp, C: NeighborIndex of the face, 0/4 x1, 8/12 x2, 24/28 x3
      const int gn = cub ? 0 : nb(m,((dir == 2) ? 24 : 8*dir) + 4*lin).gid;
      if (cub ? !(x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io))
              : !(gn < gid(m) || (gn == gid(m) && lin))) return;
      const int c[3] = {c0, c1, c2};
      Real s = 1.0;
      const int q = cub ? CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s)/2 : dir;
      for (int n = 0; n < nc; ++n) fc(m,c[n],k,j,i) = s*v(m,v0+3*n+q,kg,jg,ig);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagCsLateral
//! \brief staggered, cubed sphere: for E in work vector e, the lateral stress sums of
//! every cell and face family q, Y_q = 1/4 sum over its two faces of A n.(P' - P), P = D
//! E (Cartesian, n outward, P' of the cell across; one-sided where a side's ghosts do not
//! carry the field, CleanTransSide).  A face's F_n changes by n.div P of its half cells:
//! A n.(P_R - P_L) plus the other families' Y (CleanFaceStress), exact for uniform P.
//! The owner computes Y and the ghosts take it (a seam face sees the true geometry).

void RadiationVet::CleanStagCsLateral(const int e, const bool st, const int part) {
  if (!pmy_pack->pmesh->use_cubed_sphere && !cln_ysum) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x = CleanMakeArgs();
  const int ks = indcs.ks, js = indcs.js, is = indcs.is;
  auto cy = cln_csy;
  // the sums' array and first slot (off the cubed sphere cln_yl, x.yl)
  const bool cub = pmy_pack->pmesh->use_cubed_sphere, wl = (CLW_L < cln_w.ns);
  auto ya = !cub ? cln_yl : (wl ? cln_w.s : cln_w.o);
  const int y0 = !cub ? 0 : (wl ? CLW_L : CLW_L - cln_w.ns);
  // st: from the pass's coefficients (cln_csyb, CleanPreconditionerSetup)
  if (st) {
    par_for("cln_cslst", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
            indcs.js, indcs.je, indcs.is, indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int q = 0; q < 3; ++q) {
        const int dk = (q == 2), dj = (q == 1), di = (q == 0);
        // a seam side's slot holds the cell 2 away in panel (cln_csyb)
        const int ol = CleanTransSide(x, q, m, k, j, i, 0) ? -1 : 2;
        const int oh = CleanTransSide(x, q, m, k, j, i, 1) ? 1 : -2;
        for (int b = 0; b < 3; ++b) {
          x.w(m,CLW_L+3*q+b,k,j,i) = cy(m,18+3*q+b,k-ks,j-js,i-is)*x.w(m,e,k,j,i) +
              ((cy(m,6*q+b,k-ks,j-js,i-is) != 0.0) ?   // none past a boundary
               cy(m,6*q+b,k-ks,j-js,i-is)*x.w(m,e,k+ol*dk,j+ol*dj,i+ol*di) : 0.0) +
              ((cy(m,6*q+3+b,k-ks,j-js,i-is) != 0.0) ?
               cy(m,6*q+3+b,k-ks,j-js,i-is)*x.w(m,e,k+oh*dk,j+oh*dj,i+oh*di) : 0.0);
        }
      }
    });
    // read by the faces' own cells: depth 1, the two families other than the face's
    // (CleanStagCsLat; part 1: CleanApply finishes it)
    CleanGhostFill(cln_w, CLW_L, 9, part, 1, {0770, 0707, 0077, 0});
    return;
  }
  auto f = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real p0[6];   // the cell's own P, once
    for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
    for (int q = 0; q < 3; ++q) {
      const int dk = (q == 2), dj = (q == 1), di = (q == 0);
      Real y[3] = {0.0, 0.0, 0.0};
      // a side whose ghosts do not carry the field drops out, the other side's weight
      // doubled (one-sided, CleanTransSide)
      const int vl = CleanTransSide(x, q, m, k, j, i, 0);
      const int vh = CleanTransSide(x, q, m, k, j, i, 1);
      // a seam side keeps its face, the ghost's P extrapolated in panel, clamped
      // (CleanCsSeamCoef); an open boundary drops its side, the other's weight doubled
      const bool o2 = (vl + vh == 1) && x.mbc(m,2*q+vl) == BoundaryFlag::panel;
      const Real wq = (vl + vh == 2 || o2) ? 0.25 : 0.5;
      // clean_stag_lat_sums off the cubed sphere: A/V of a Cartesian face 1/dx, V = 1;
      // no family past the mesh's dimension
      for (int s = -1; s <= 1 && q < x.nf; s += 2) {
        const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
        const bool ex = o2 && !CleanTransSide(x, q, m, k, j, i, s > 0);
        if (!ex && !CleanTransSide(x, q, m, k, j, i, s > 0)) continue;
        const Real a = wq*s*((x.cub || x.sph) ? x.area[q](m,kf,jf,iff) : x.idx[q]);
        if (!ex) {   // each component's P(s) - P(0) and the normal once
          Real dp[6], nf[3];
          for (int n = 0; n < 6; ++n) {
            dp[n] = x.cs(m,CLN_D11+n,k+s*dk,j+s*dj,i+s*di)*x.w(m,e,k+s*dk,j+s*dj,i+s*di) -
                    p0[n];
          }
          for (int c = 0; c < 3; ++c) nf[c] = CleanNf(x, c, q, m, kf, jf, iff);
          for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) y[b] += a*nf[c]*dp[(b == c) ? b : 2 + b + c];
          }
          continue;
        }
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
            auto p = [&](const int o) {
              return x.cs(m,n,k+o*dk,j+o*dj,i+o*di)*x.w(m,e,k+o*dk,j+o*dj,i+o*di);
            };
            Real cf[3] = {0.0, 0.0, 0.0};
            if (ex) CleanCsSeamCoef(x, n, q, m, k, j, i, -s, cf);
            y[b] += a*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff)*(ex ? cf[0]*p(0) + cf[1]*p(-s) +
                                                         cf[2]*p(-2*s) : p(s) - p(0));
          }
        }
      }
      for (int b = 0; b < 3; ++b) ya(m,y0+3*q+b,k,j,i) = y[b];
    }
  };
  // the cells it sends first, its exchange in flight behind the inner ones
  const bool ov = (cln_hmax >= 9);
  const std::function<void()> post = [&]() {
    CleanGhostFill(ya, y0, 9, true, ov ? 1 : 0, 1, {0770, 0707, 0077, 0});   // as above
  };
  const std::function<void()> fin = [&]() {
    if (ov) CleanGhostFill(ya, y0, 9, true, 2, 1);
  };
  CleanOverlapFor("cln_cslat", pmy_pack, ov, f, fin, -1, post);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanStagDrag
//! \brief clean_stag_vimp (staggered, gas coupling): the gas velocity implicit in the
//! face momentum, eliminated face by face so the E solve stays scalar.  A change dm_f of
//! face f's momentum moves its two half cells (mass rho_c/2, share wt_c = chi_c/(chi_L +
//! chi_R)) by dv_c = 2 wt_c dm_f/rho_c, and the advective flux by (1 - theta) dA_f =
//! beta c^2 dm_f, beta = (1 - theta) sum_c wt_c (1 + D_nn,c) E_c/(rho_c c^2).  So
//!   F = (1 - b) F^H(E) + b [Q - h c^2 d_n(D E)],  b = beta/(1 + beta),
//!   Q = F^k + h c^2 d_n(D E^k) = G - c^2 dm_f^k (CLF_DB, CLF_DQ),
//! and the energy flux takes the same increment; at convergence F = F^H (the lagged
//! velocity's solution).  first: dm_f^k = 0, Q = G (with a reference force: Q of E^k, no
//! drag kick); b = 0 while the gas momentum is frozen.  vel (CleanVelocityUpdate with the
//! F solve's F): first the face's CleanMomentumChange term less its share wt (CLF_PT).

void RadiationVet::CleanStagDrag(const bool first, const bool vel, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // cubed sphere: the lateral stress of the iterate (lat false: the caller's)
  if (lat) CleanStagCsLateral(CLW_EK);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  xh.drag = x.drag && !first;   // the flux the iterate's velocity came from
  auto uh = FluidRef::Get(pmy_pack).u0;
  const bool off = cln_vzero || (first && force_ref == VET_FREF_CLEAN && !cln_fref_set);
  const bool lag = !first || (force_ref != VET_FREF_NONE);
  const bool dr = cln_vimp;
  const int nf = x.nf;   // one launch for the faces of all directions
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  CleanFor("cln_drag", nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // in device memory: no per-thread copy
    for (int dir = 0; dir < nf; ++dir) {
      const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      if (i > ie+di || j > je+dj || k > ke+dk) continue;
      auto &fc = x.fc[dir];
      if (CleanBoundaryFace(x, dir, m, k, j, i)) continue;   // no F equation
      const Real s = (vel || lag) ? CleanFaceStress(x, CLW_EK, dir, m, k, j, i, true)
                                  : 0.0;
      // vel: the F solve's fs0 just took F^H of E_k with these weights (CLF_PT)
      const Real fq = (vel && lag) ? fc(m,CLF_PT,k,j,i) : 0.0;
      if (vel) {
        fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) + x.h*s -
                             x.h*fc(m,CLF_N+1,k,j,i);
      }
      if (!dr) continue;
      const Real q = !lag ? fc(m,CLF_G,k,j,i) :
                     (vel ? fq : CleanFaceFlux(xh, CLW_EK, dir, m, k, j, i, true)) +
                     x.h*x.c*x.c*s;
      // cubed sphere: n.D.n, the owners' half-cell V (as CleanMomentumChange) and the
      // exchange's rho (CLW_JV + 2: a seam ghost's own cell), so both sides build one b
      auto vo = [&](const int kk, const int jj, const int ii) {
        return x.cub ? x.g(m,CLG_F+3,kk,jj,ii) : 1.0;
      };
      const Real vl = vo(k-dk, j-dj, i-di), vr = vo(k, j, i);
      const Real xl = x.op(m,VET_OP_T,k-dk,j-dj,i-di)*vl;
      const Real xs = xl + x.op(m,VET_OP_T,k,j,i)*vr, wl = (xs > 0.0) ? xl/xs : 0.5;
      auto ce = [&](const int kk, const int jj, const int ii) {
        const Real d = CleanStagDnn(x, dir, m, k, j, i, kk, jj, ii);
        return (1.0 + d)*fmax(x.w(m,CLW_EK,kk,jj,ii), 0.0)*0.5*(vl + vr)/
               (vo(kk, jj, ii)*(x.cub ? x.w(m,CLW_JV+2,kk,jj,ii) : uh(m,IDN,kk,jj,ii)));
      };
      const Real bt = off ? 0.0 : (1.0 - fc(m,CLF_TH,k,j,i))*(wl*ce(k-dk, j-dj, i-di) +
                      (1.0 - wl)*ce(k, j, i))/(x.c*x.c);
      fc(m,CLF_DB,k,j,i) = bt/(1.0 + bt);
      fc(m,CLF_DQ,k,j,i) = q;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanLzBook
//! \brief clean_cs_lz_abs = 2: books the cells' sum s of the L_z the gas did not take (a
//! write-back, or with start the stage start's increment) as its gas momentum enters the
//! step: stage-1 part weighted as in the stage-2 start (Heun 1/2), stage-2 parts pending
//! until stage 2 is accepted (s = 0; a backward-Euler redo drops them, Time2Restore).
//! Total and K1's sum in <rad_vet> (restart).  3-D spherical polar: s is the L_z that
//! crossed the boundaries (CleanStagMomentum's dlz), <rad_vet>/clean_lz_ledger.

void RadiationVet::CleanLzBook(const Real s, const bool start) {
  if (t2_inprobe) return;   // a helper solve, put back (Time2RotCorrect)
  if (t2_solve == VET_T2S_STAGE1) {
    cln_lzp = (start ? 0.0 : cln_lzp) + ((t2_tab == 0 || t2_tab == 2) ? 0.5 : t2_w0)*s;
    return;
  }
  if (t2_solve == VET_T2S_STAGE2 && start) {
    cln_lzp += s;
    return;
  }
  cln_lzl += (t2_solve == VET_T2S_STAGE2) ? cln_lzp : cln_lzp + s;
  cln_lzp = 0.0;
  const bool cs = (cln_lza == 2);
  cln_lzpin->SetReal("rad_vet", VetLzKey(cln_lzpin, cs, false), cln_lzl);
  cln_lzpin->SetReal("rad_vet", VetLzKey(cln_lzpin, cs, true), cln_lzs1);
  if (global_variable::my_rank == 0 && cln_verbose > 0 && !cln_insteady) {
    std::printf("CLEAN LZ%s cycle %d t2s %d ledger %.16e\n", cs ? "" : "B",
                pmy_pack->pmesh->ncycle, t2_solve, cln_lzl);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanWbFct
//! \brief clean_wb_fct: Zalesak limit of the written faces' energy flux Phi (cln_wbdiv's,
//! fp 0..2) so that E > e_floor after the write-back by construction.  A cell's E is
//! Q - P: Q its E without face transport (cln_wbcell probe, fp 3, less e_floor and
//! round-off) plus inflow, P its outflow.  A cell with Q - R P not > 0 takes R = min(1,
//! (1 - 1e-12) Q/P) (0, counted 'left', where Q <= 0); each face is scaled by its donor's
//! R; rounds until no R falls (an R that fell carries its sign bit, so the next round
//! evaluates only its neighbours).  The limited faces' (1 - s) Phi leave qf (CLW_Q), pf
//! and the x1 boundary ledger.  Face F is not touched: momentum and L_z exact.

void RadiationVet::CleanWbFct(const CleanArgs &x, const CleanArgs *px,
                              const DvceArray5D<Real> &fp, const bool pfo,
                              const DvceArray5D<Real> &pf) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, ie = id.ie, js = id.js, je = id.je, ks = id.ks, ke = id.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1, nf = x.nf;
  if (iw.extent_int(1) <= Time2FctSlot(1)) {
    CleanFatal("<rad_vet>/clean_wb_fct needs time2_face_fct's iw slots");
  }
  auto iw_ = iw;
  const int ca = Time2FctSlot(0), cb = Time2FctSlot(1);
  for (const int c : {ca, cb}) {   // R = 1, ghosts too: physical ones stay unlimited
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, c, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), 1.0);
  }
  // a cell's Q and P at the donors' R in slot c; a face's scale, its donor's R in ca
  auto qp = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, const int c,
                          Real &q, Real &po) {
    const CleanArgs &x = *px;
    q = fp(m,3,k,j,i);
    po = 0.0;
    for (int n = 0; n < nf; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      const Real gl = x.h*CleanDivW(x, n, m, k, j, i, 0)*fp(m,n,k,j,i);   // + = in
      const Real gh = x.h*CleanDivW(x, n, m, k, j, i, 1)*fp(m,n,k+dk,j+dj,i+di);
      if (gl > 0.0) {q += fabs(iw_(m,c,k-dk,j-dj,i-di))*gl;} else {po -= gl;}
      if (gh > 0.0) {po += gh;} else {q -= fabs(iw_(m,c,k+dk,j+dj,i+di))*gh;}
    }
  };
  auto sf = KOKKOS_LAMBDA(const int n, const int m, const int k, const int j,
                          const int i) {
    const Real f = fp(m,n,k,j,i);
    return (f > 0.0) ? iw_(m,ca,k-(n == 2),j-(n == 1),i-(n == 0)) :
                       ((f < 0.0) ? iw_(m,ca,k,j,i) : 1.0);
  };
  const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
  int nr = 0, sa = ca, sb = cb;
  for (; nr < 4096; ++nr) {   // R only falls: until none newly limited
    Real nv = 0.0;
    const bool all = (nr == 0);
    Kokkos::parallel_reduce("cln_wbfr", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &a) {
      const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
      const int i = is + idx%n1;
      const Real rc = fabs(iw_(m,sa,k,j,i));
      bool act = all;
      for (int n = 0; n < nf && !act; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        act = Kokkos::signbit(iw_(m,sa,k-dk,j-dj,i-di)) ||
              Kokkos::signbit(iw_(m,sa,k+dk,j+dj,i+di));
      }
      Real r = rc;
      if (act) {
        Real q, po;
        qp(m, k, j, i, sa, q, po);
        const Real rq = (q > 0.0 && po > 0.0) ? fmin(1.0, (1.0 - 1.0e-12)*q/po) : 0.0;
        if (!(q - rc*po > 0.0) && rq < rc) {
          r = -rq;
          a += 1.0;
        }
      }
      iw_(m,sb,k,j,i) = r;
    }, Kokkos::Sum<Real>(nv));
    Time2FctHalo(sb);
    std::swap(sa, sb);
    if (!(CleanReduce(nv, false) > 0.0)) break;
  }
  if (nr == 4096) {
    ZalesakCapWarn("clean_wb_fct", pmy_pack->pmesh->ncycle, t2_solve);
  }
  if (sa != ca) {   // the last round's R, no sign bit set, ghosts too
    Kokkos::deep_copy(DevExeSpace(), Kokkos::subview(iw, Kokkos::ALL, ca, Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), Kokkos::subview(iw, Kokkos::ALL, sa,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (nr == 0) return;   // nothing limited: bitwise
  Real lv[2] = {0.0, 0.0};   // 'left': not > 0 at R; the outflow the limit kept
  Kokkos::parallel_reduce("cln_wbfl", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a, Real &b) {
    const CleanArgs &x = *px;
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
    const int i = is + idx%n1;
    Real q, po;
    qp(m, k, j, i, ca, q, po);
    const Real rc = iw_(m,ca,k,j,i);
    if (!(q - rc*po > 0.0)) a += 1.0;
    b += (1.0 - rc)*po*CleanVol(x, m, k, j, i);
    Real dq = 0.0;   // qf less the limited faces' share
    for (int n = 0; n < nf; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      dq += CleanDivF(x, n, m, k, j, i, (1.0 - sf(n, m, k, j, i))*fp(m,n,k,j,i),
                      (1.0 - sf(n, m, k+dk, j+dj, i+di))*fp(m,n,k+dk,j+dj,i+di));
    }
    x.w(m,CLW_Q,k,j,i) -= x.h*dq;
  }, Kokkos::Sum<Real>(lv[0]), Kokkos::Sum<Real>(lv[1]));
  lv[0] = CleanReduce(lv[0], false);
  lv[1] = CleanReduce(lv[1], false);
  if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(lv[1]))) {
    std::printf("CLEAN %sWBFCT cycle %d t2s %d rounds %d left %.0f kept %.6e\n",
                std::isfinite(lv[1]) ? "" : "WARNING ", pmy_pack->pmesh->ncycle,
                t2_solve, nr, lv[0], lv[1]);
  }
  for (int n = 0; n < nf && pfo; ++n) {   // the stored Phi: what moved E
    CleanFor("cln_wbfs", nmb1, ks, ke+(n == 2), js, je+(n == 1), is, ie+(n == 0),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pf(m,n,k,j,i) -= (1.0 - sf(n, m, k, j, i))*fp(m,n,k,j,i);
    });
  }
  if (cln_lacc.extent_int(0) > nmb1) {   // the history's boundary energy (cln_wbbnd)
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;
    const Real h = x.h;
    par_for("cln_wbfb", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      if (CleanBoundaryFace(*px, 0, m, k, j, i)) {
        la(m,s,k,j) -= h*ar(m,k,j,i)*(1.0 - sf(0, m, k, j, i))*fp(m,0,k,j,i);
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanWriteBack
//! \brief stores the solution: face fluxes, E, cell fluxes, gas momentum and energy, and
//! the stage slopes K = (U - U_old)/h.  The gas takes exactly the momentum and energy the
//! radiation gives up (less the reference force with force_reference), so both totals
//! change only by face fluxes; with wb_arad (force_reference_work = split) E pays the
//! reference work h d arad_ref.v; esrc is the radiation's alone.  Cubed sphere: the gas
//! takes CleanStagMomentum's vector in the hydro's components, the output cell flux is
//! the face momenta's vector (frame of y): energy and L_z conserved to round-off.

void RadiationVet::CleanWriteBack() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool slopes = (t2_solve != VET_T2S_NONE);
  const bool feedback = cln_feedback;
  const Real h = dt_sub;
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanStagCsLateral(CLW_X);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = cln_split;   // clean_stag_split: the faces store the momentum state
  auto kk = (t2_solve == VET_T2S_STAGE1) ? t2k2 : t2k1;
  // h div Phi of the solution (CLW_Q) before the F solve below moves the memory d
  // (clean_stag_fkin): E takes the residual of the operator it was solved with; with
  // clean_stag_split also the cell flux of Phi (the faces keep F^H)
  auto u0_ = u0;
  // face tie: a high face's Phi is its owner's (CLW_L + dir of the next cell), so energy
  // telescopes across MeshBlock faces and the theta-periodic wrap; cubed sphere: each
  // shared face's CLF_PT the lower gid's (CleanStagCsTie)
  const bool tie = !x.cub && !cln_tgid;
  if (x.cub || cln_tgid) CleanStagCsTie(1, CLF_PT);
  if (tie) {
    const int nf = x.nf;
    CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
    CleanFor("cln_wbphi", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // in device memory: no per-thread copy
      // CleanStagRealizable's Phi except on physical faces
      for (int d = 0; d < nf; ++d) {
        x.w(m,CLW_L+d,k,j,i) = CleanBoundaryFace(x, d, m, k, j, i) ?
            CleanFaceFlux(x, CLW_X, d, m, k, j, i, true) : x.fc[d](m,CLF_PT,k,j,i);
      }
    });
    CleanGhostFill(cln_w, CLW_L, nf);
  }
  // clean_stag_xdonor: the cell flux (radiation momentum) of the central cross flux and
  // face sources, as the gas momentum's stress: Phi's donor form equals it only at the
  // exact fixed point of q
  CleanArgs xc = x;
  xc.xq = -1;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const CleanArgs *pxc = CleanDevArgs(xc, 2);
  // time2_face_fct: the energy flux Phi of every face, which moved E, with its slope
  const bool pfo = t2_pfst && slopes && t2pf1.extent_int(0) > nmb1;
  auto pf = (t2_solve == VET_T2S_STAGE1) ? t2pf2 : t2pf1;
  // xdonor: each face's central cross flux once, for both its cells, in cln_xf's scratch
  // (its users write it before they read it)
  auto phq = cln_phic;   // the cell flux of the faces' Phi (output vet_phi)
  DvceArray5D<Real> xcv;
  if (xh.hi && x.xq >= 0) {
    xcv = ScratchView<DvceArray5D<Real>>("cln_xf", nmb1 + 1, 3, cln_w.s.extent_int(2),
                                         cln_w.s.extent_int(3), cln_w.s.extent_int(4));
    for (int n = 0; n < x.nf; ++n) {
      CleanFor("cln_wbxc", nmb1, ks, ke+(n == 2), js, je+(n == 1), is, ie+(n == 0),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        CleanDirCase(n, [&](auto d) {
          constexpr int nd = decltype(d)::value;
          xcv(m,nd,k,j,i) = CleanFaceFlux<nd>(*pxc, CLW_X, nd, m, k, j, i, true);
        });
      });
    }
  }
  // clean_wb_fct: each face's Phi for CleanWbFct, components 0..2 of cln_wbfp
  const bool wf = cln_wbfct && feedback && !x.vzero;
  const DvceArray5D<Real> wfp = !wf ? DvceArray5D<Real>() :
      ScratchView<DvceArray5D<Real>>("cln_wbfp", nmb1 + 1, 4, u0.extent_int(2),
                                     u0.extent_int(3), u0.extent_int(4));
  CleanFor("cln_wbdiv", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // device memory: no copy
    // h div Phi with the faces' Phi of CleanStagRealizable (CLF_PT), CleanFaceFlux on
    // physical faces; with the tie the low faces' CLW_L and the high faces' owner's
    Real dv = 0.0, pc[6];   // cubed sphere: each face's (CleanStagCellVector)
    for (int n = 0; n < 3; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      if (n >= x.nf) {
        if (xh.hi) phq(m,n,k,j,i) = 0.0;
        continue;
      }
      const Real fl = tie ? x.w(m,CLW_L+n,k,j,i) : (CleanBoundaryFace(x, n, m, k, j, i) ?
          CleanFaceFlux(x, CLW_X, n, m, k, j, i, true) : x.fc[n](m,CLF_PT,k,j,i));
      const bool bh = CleanBoundaryFace(x, n, m, k+dk, j+dj, i+di);
      const Real fu = bh ? CleanFaceFlux(x, CLW_X, n, m, k+dk, j+dj, i+di, true) :
                           x.fc[n](m,CLF_PT,k+dk,j+dj,i+di);
      const Real fh = (tie && !bh) ? x.w(m,CLW_L+n,k+dk,j+dj,i+di) : fu;
      dv += CleanDivF(x, n, m, k, j, i, fl, fh);
      if (wf) {
        wfp(m,n,k,j,i) = fl;
        if (n == 0 ? i == x.ie : (n == 1 ? j == x.je : k == x.ke)) {
          wfp(m,n,k+dk,j+dj,i+di) = fh;
        }
      }
      if (xh.hi && x.xq >= 0) {
        phq(m,n,k,j,i) = 0.5*(xcv(m,n,k,j,i) + xcv(m,n,k+dk,j+dj,i+di));
      } else if (xh.hi) {
        phq(m,n,k,j,i) = 0.5*(fl + fu);
      }
      if (xh.hi && x.cub) {
        pc[2*n] = (x.xq >= 0) ? xcv(m,n,k,j,i) : fl;
        pc[2*n+1] = (x.xq >= 0) ? xcv(m,n,k+dk,j+dj,i+di) : fu;
      }
      if (pfo) pf(m,n,k,j,i) = fl;
      // the MeshBlock's high face: the Phi that moved E (the owner's with the tie), so
      // time2_face_fct limits one flux on both sides
      if (pfo && (n == 0 ? i == x.ie : (n == 1 ? j == x.je : k == x.ke))) {
        pf(m,n,k+dk,j+dj,i+di) = fh;
      }
    }
    if (xh.hi && x.cub) {   // cubed sphere: the faces' Phi's vector, frame of y
      Real fv[3];
      CleanStagCellVector(x, m, k, j, i, pc, fv);
      CleanStagHydro(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) phq(m,d,k,j,i) = fv[d];
    }
    x.w(m,CLW_Q,k,j,i) = x.h*dv;
  });
  if (pfo) {
    if (t2_solve == VET_T2S_STAGE1) {t2pf2_ok = true;} else {t2pf1_ok = true;}
  }
  // the history's boundary energy: + h A Phi of the physical x1 faces (CleanBndAcc)
  if (cln_lacc.extent_int(0) > nmb1) {
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;   // x1 face area of the clean fluxes
    const bool stg = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
    par_for("cln_wbbnd", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      la(m,s,k,j) = !CleanBoundaryFace(x, 0, m, k, j, i) ? 0.0 :
                    la(m,stg ? 2+s : s,k,j) + h*ar(m,k,j,i)*CleanFaceFlux(x, CLW_X, 0, m,
                                                                     k, j, i, true);
    });
  }
  // face F of the F solve and its dissipation (lateral sums of X: above)
  CleanStagFSolve(CLW_X, true, false);
  xh.fs = true;
  xh.tf = !x.cub && !cln_tgid;
  xh.fa = x.cub || cln_tgid;   // the tied T in CLF_PT (CleanStagFaceTie)

  // clean_cs_lz_tel: CleanCsTorque's stress torque xs is zero-sum noise per cell; the gas
  // takes only its net X = sum xs, as the solid-body rotation dOmega = X/I of least
  // kinetic energy, I = sum rho V R^2.  tqv: the gas's L_z gl and the ledger's; lxs, lxi:
  // xs and rho V R^2, lzs their sums X and I
  const bool lzt = cln_lzt && feedback;
  DvceArray5D<Real> tqv;
  Real lzs[2] = {0.0, 0.0};
  // (clean_lz_spin: X, the L_z of fcone's cut of the faces with an arm (CleanStagFCone)
  // and of clean_wb_mcap's cap (a pass of the arm families' terms), the gas's as the
  // rigid rotation dOmega = X/I, I = sum rho V R^2; the cubed sphere: none, its gas phi
  // momentum is lz_tel's torque of the uncut, uncapped face terms)
  const bool spn = cln_lzsp && feedback && !x.vzero && !x.cub;
  Real lzx = 0.0;
  if (spn) {   // (per cell, then summed in a fixed order: each cell's own halves)
    CleanArgs xm = xh;
    xm.mcl = cln_wbmc;
    xm.mclc = DvceArray1D<Real>("cln_lzmc", 1);
    auto fcx = cln_fcx;
    auto lxx = ScratchView<DvceArray4D<Real>>("cln_lxx", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    CLEAN_HIPX(const CleanArgs *pxm = CleanDevArgs(xm, 0), *pxh = CleanDevArgs(xh, 1);)
    CleanFor("cln_lzmc", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &xm = *pxm, &xh = *pxh;)
      Real s = 0.0;
      for (int d = 0; d < xh.nf; ++d) {
        if (!CleanLzFamily(xh.sph, xh.cub, d)) continue;
        Real cf[2], wf[2], a = 0.0;
        const Real u = CleanMomentumChange(xh, CLW_X, d, m, k, j, i, nullptr, cf, wf);
        for (int sd = 0; sd < 2; ++sd) {   // (fcone's cut, as the gas's share of it)
          a += wf[sd]*fcx(m,d,k+sd*(d == 2),j+sd*(d == 1),i+sd*(d == 0));
        }
        if (xm.mcl) a += u - CleanMomentumChange(xm, CLW_X, d, m, k, j, i);
        s += a*CleanDirLzArm(xh, d, m, k, j, i);
      }
      lxx(m,k,j,i) = s*CleanVol(xh, m, k, j, i);
    });
    lzx = CleanSumCells(lxx);
    lzs[0] = lzx;
    if (!lzt) {
      auto ug = FluidRef::Get(pmy_pack).u0;
      auto lxi = ScratchView<DvceArray4D<Real>>("cln_lxi", nmb1 + 1, u0.extent_int(2),
                                                u0.extent_int(3), u0.extent_int(4));
      CleanFor("cln_lzsi", nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        lxi(m,k,j,i) = ug(m,IDN,k,j,i)*CleanVol(x, m, k, j, i)*
                       SQR(CleanLzArm(x, m, k, j, i));
      });
      lzs[1] = CleanSumCells(lxi);
    }
  }
  if (lzt) {
    tqv = ScratchView<DvceArray5D<Real>>("cln_tqv", nmb1 + 1, 3, u0.extent_int(2),
                                         u0.extent_int(3), u0.extent_int(4));
    auto lxs = ScratchView<DvceArray4D<Real>>("cln_lxs", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    auto lxi = ScratchView<DvceArray4D<Real>>("cln_lxi", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    px = CleanDevArgs(x, 0);   // xh changed: the F solve's F
    pxh = CleanDevArgs(xh, 1);
    auto ug = FluidRef::Get(pmy_pack).u0;
    const bool bf = cln_lzbf;
    CleanFor("cln_wbtq", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh;)
      Real rl, xs, q[3][3];
      const Real gl = CleanCsTorque(xh, CLW_X, m, k, j, i, rl, xs, bf);
      const Real rc = CleanCsFrame(x, m, k, j, i, q);
      lxs(m,k,j,i) = xs;
      lxi(m,k,j,i) = ug(m,IDN,k,j,i)*x.vol(m,k,j,i)*rc*rc;
      tqv(m,1,k,j,i) = gl - xs;   // the absorption torque
      tqv(m,2,k,j,i) = -(gl + rl);
    });
    lzs[0] += CleanSumCells(lxs);
    lzs[1] = CleanSumCells(lxi);
  }
  const Real lzw_ = (lzs[1] > 0.0) ? lzs[0]/lzs[1] : 0.0;   // (dOmega)
  if (spn && global_variable::my_rank == 0 && lzx != 0.0) {
    std::printf("CLEAN LZSPIN cycle %d t2s %d: L_z of the cut and capped momentum X "
                "%.6e, the gas's as dOmega %.6e (I %.6e)\n", pmy_pack->pmesh->ncycle,
                t2_solve, lzx, lzw_, lzs[1]);
  }

  // face fluxes of the solution
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = VET_T2_F1 + dir;
    CleanFor("cln_wbface", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real f = fc(m,CLF_N,k,j,i);
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
  const bool wsp = (fref == VET_FREF_WB_ARAD || fref == VET_FREF_WB_MD);
  const bool fcap = (fref == VET_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  const bool wpa = wsp && FrefWaccOn();   // fref_wacc: CleanRightHandSide
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : h/pmy_pack->pmesh->dt;
  const bool wex = wsp && fref_wacc_t2 &&
                   time_scheme == VET_TIME_HESDIRK2;   // no slope carries a payment
  const bool fm = icr_fmeas && cln_vzero;   // CleanSteadyStart: icr_dm
  auto fmd = icr_dm;
  const bool eso = esrc_on;
  auto es = esrc;
  const Real arad_ = arad;
  // wb_arad_md, clean_init = steady: the steady state's reference; cln_fref 0..2 its
  // acceleration (components of dm), 3 + 2 dir + side the potential the hydro walks of
  // <hydro>/wb_dirs = all take along dir: the line integral of 2 cf/(h rho) over each
  // half cell (side 0 face to centre, 1 centre to face).  None outside the Roche lobe
  // (problem/wb_lobe_cut)
  const bool mdc = cln_mdcap;
  const Real phx = fl.ph ? fl.ph->wb_phimax : (fl.pm ? fl.pm->wb_phimax : 1.0e300);
  auto phc = fl.phicc_wb;
  auto cx1 = pmy_pack->pcoord->x1v, cx2 = pmy_pack->pcoord->x2v;
  auto cx3 = pmy_pack->pcoord->x3v, cf1 = pmy_pack->pcoord->xx1f;
  auto cf2 = pmy_pack->pcoord->xx2f, cf3 = pmy_pack->pcoord->xx3f;
  const Real tgf = cln_gtf ? eos.tfloor : 0.0;   // clean_gas_tfloor
  const bool gtf = cln_gtf, pgas = impl_pos_gas, lex = cln_wbxl && feedback && gtf;
  if (lex && !(tgf > 0.0)) CleanFatal("<rad_vet>/clean_wb_exch_local needs a tfloor > 0");
  const Real pgf = impl_pos_gas_frac, efl = e_floor, pgk = impl_pos_flr_gk;
  // clean_wb_floor_diag N: the first N floored cells of each write-back, their terms
  const int fdn = cln_wbfdg;
  DvceArray1D<Real> fdc;
  if (fdn > 0) fdc = DvceArray1D<Real>("cln_wbfdc", 1);
  const int fcy = pmy_pack->pmesh->ncycle, fts = t2_solve;
  auto pc_ = pos_cnt_d;
  const bool pflr = impl_pos_floor, pfc = pflr && impl_pos_flr_cons;
  // the E floor's raise neither gas nor column covered (created, VET_POS_FLR_UN): its
  // change in this write-back, CLEAN WARNING below
  Real un0 = 0.0;
  if (pflr) Kokkos::deep_copy(un0, Kokkos::subview(pc_, VET_POS_FLR_UN));
  if (pfc && pfc_def.extent_int(0) != nmb1 + 1) {
    auto &icp = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(pfc_def, nmb1 + 1, icp.nx3 + 2*icp.ng, icp.nx2 + 2*icp.ng,
                    icp.nx1 + 2*icp.ng);
  }
  if (pfc) Kokkos::deep_copy(pfc_def, 0.0);
  auto pfd = pfc_def;
  const bool pgc = pfc && pgas && impl_pos_gas_cons;   // implicit_pos_gas_cons
  if (pgc && pgc_def.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(pgc_def, nmb1 + 1, pfc_def.extent_int(1), pfc_def.extent_int(2),
                    pfc_def.extent_int(3));
  }
  if (pgc) Kokkos::deep_copy(pgc_def, 0.0);
  auto pgd = pgc_def;
  const bool xwb = (cln_xte > 0) && feedback && !cln_vzero;
  const bool cls = xwb;   // the closure in CLW_PW, printed below
  // (clean_wb_mcap: the gas's momentum at the passes' cap, that of the E row's work; the
  // terms capped and their momentum withheld per face family, CLEAN MCAP below)
  const bool wmc = cln_wbmc && feedback && !x.vzero;
  if (wmc) {
    xh.mcl = true;
    xh.mclc = DvceArray1D<Real>("cln_wbmcc", 4);
  }
  px = CleanDevArgs(x, 0);   // the F solve took the slots; xh changed
  pxh = CleanDevArgs(xh, 1);
  const bool lzb = cln_lzbook && feedback;
  DvceArray4D<Real> lzw;
  if (lzb) lzw = ScratchView<DvceArray4D<Real>>("cln_lzw", nmb1 + 1, u0.extent_int(2),
                                                 u0.extent_int(3), u0.extent_int(4));
  // (prb, clean_wb_fct's probe: no side effects, each cell's E without its face
  // transport, e + qf, into cln_wbfp 3 for CleanWbFct, before the cells' write-back)
  // (clean_dbg_inf: once, in the first stage-1 write-back of that cycle)
  const bool dinf = (cln_dbginf == fcy) && t2_solve == VET_T2S_STAGE1 && !t2_inprobe &&
                    global_variable::my_rank == 0;
  if (dinf) cln_dbginf = -1;
  auto wbc = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                           const bool prb) CLEAN_HIP_INLINE {
    CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh;)   // in device memory
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
    Real dg[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // clean_wb_floor_diag: eg wk wp dw rd
    const Real d0 = x.w(m,CLW_JV,k,j,i);   // the stage start's eint deficit, <= 0
    if (feedback) {
      Real mg[3], pg[3], rf[3], dlz = 0.0;
      CleanStagMomentum(xh, CLW_X, m, k, j, i, dm, mg, &dlz);
      if (spn) {   // (clean_lz_spin: dOmega rho z x r; the ledger's withheld L_z back)
        const Real w = lzw_*uh(m,IDN,k,j,i), r = CleanLzArm(x, m, k, j, i);
        if (x.sph) {
          dm[2] += w*r;
        } else {
          dm[0] -= w*x.x2v(m,j);
          dm[1] += w*x.x1v(m,i);
        }
        dlz -= w*r*r*CleanVol(x, m, k, j, i);
      }
      if (lzt) {   // (clean_cs_lz_tel: the gas phi momentum of CleanCsTorque, cln_wbtq)
        Real q[3][3];
        const Real rc = CleanCsFrame(x, m, k, j, i, q);
        const Real dp = tqv(m,1,k,j,i)/(x.vol(m,k,j,i)*rc) + lzw_*uh(m,IDN,k,j,i)*rc -
                        (q[2][0]*dm[0] + q[2][1]*dm[1] + q[2][2]*dm[2]);
        // force_reference = clean: the reference withholds this phi row too (cln_fref 3),
        // its torque booked as that of cln_fref 0..2
        const Real dw = (fref != VET_FREF_CLEAN) ? 0.0 : (fcap ? dp : h*fr(m,3,k,j,i));
        if (fcap) fr(m,3,k,j,i) = dp/h;
        for (int n = 0; n < 3; ++n) dm[n] += (dp - dw)*q[2][n];
        dlz = tqv(m,2,k,j,i) + x.vol(m,k,j,i)*rc*dw;
      }
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      CleanStagRefForce(x, fref, fr, aref, uh(m,IDN,k,j,i), m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {   // force_reference
        if (fref != VET_FREF_NONE) dm[n] -= h*rf[n];
      }
      if (lzb && fref == VET_FREF_CLEAN) {   // the reference torque withheld: booked
        Real q[3][3];
        const Real r = x.cub ? CleanCsFrame(x, m, k, j, i, q) : x.x1v(m,i)*x.tg(m,2,j);
        dlz += h*x.vol(m,k,j,i)*r*(x.cub ? q[2][0]*rf[0] + q[2][1]*rf[1] + q[2][2]*rf[2]
                                         : rf[2]);
      }
      if (fm) fmd(m,k,j,i) = CleanStagRadial(x, m, k, j, i, dm);
      if (x.vzero) dm[0] = dm[1] = dm[2] = dlz = 0.0;   // frozen gas momentum
      if (lzb) lzw(m,k,j,i) = dlz;   // clean_cs_lz_abs = 2: CleanLzBook
      // the kinetic energy the gas gains (the E row's W at convergence); with wb_arad wr,
      // the reference work at the velocity after the kick, which E pays
      Real wk = 0.0, wr = 0.0;
      for (int n = 0; n < 3; ++n) {
        wk += (pg[n] + 0.5*dm[n])*dm[n]/uh(m,IDN,k,j,i);
        if (wsp && !x.vzero) wr += h*rf[n]*(pg[n] + dm[n])/uh(m,IDN,k,j,i);
        dh[n] = dm[n];
      }
      if (x.cub) {   // the hydro's components of dm, from the frame of y
        Real dy[3] = {dm[0], dm[1], dm[2]};
        CleanStagHydro(x, m, k, j, i, dy);
        CleanGasChange(x, m, k, j, i, dy, dh);
      }
      for (int n = 0; n < 3 && !prb; ++n) uh(m,IM1+n,k,j,i) += dh[n];
      // gas energy: exactly what the radiation gives up besides its flux divergence, so
      // total energy changes only by the face fluxes
      Real div = x.idx[0]*(f1(m,k,j,i+1) - f1(m,k,j,i));
      if (multi_d) div += x.idx[1]*(f2(m,k,j+1,i) - f2(m,k,j,i));
      if (three_d) div += x.idx[2]*(f3(m,k+1,j,i) - f3(m,k,j,i));
      // clean_stag_split: the faces store F^H, E moved by the energy flux; esrc is the
      // radiation's external source (the gas takes the exchange only)
      den = -(e - x.w(m,CLW_EO,k,j,i)) - ((xh.hi || x.sph || x.cub) ?
            x.w(m,CLW_Q,k,j,i) : h*div) + (eso ? h*es(m,k,j,i) : 0.0);
      // steady start: the eint of the iterate's T, not the residual of a step h -> inf
      if (x.vzero) {
        Real ei, p, chi_rho, chi_t, cv;
        eos.ThermoAt(uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), ei, p, chi_rho, chi_t, cv);
        den = ei - x.w(m,CLW_EI,k,j,i);
      } else {
        // the gas takes the exchange of the solution at the iterate's T, the beta terms
        // and W; E the rest (its fluxes): the linear solve's residual goes to E, not eint
        const Real t = x.w(m,CLW_T,k,j,i);
        const Real ex = h*x.c*(x.op(m,VET_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_*
                               SQR(SQR(t)) -
                               x.op(m,VET_OP_E,k,j,i)*fmax(e, 0.0));
        // the residual over the row's local diagonal (all of it in a scattering cell,
        // little where the exchange is stiff).  wpa: E pays fref_wacc; the slope keeps no
        // payment (wex: FrefWaccStage) or the payment wr, re-applied by the next stages
        const Real wp = wpa ? fra*wc(m,k,j,i) : wr;
        // clean_gas_tfloor: the gas eint changes by that of its T root (>= tfloor),
        // e(rho, T) - eint, not by the exchange (two terms ~E >> eint where r ~ 1):
        // positive by construction; E takes all the rest.  The floor only here (inside
        // the Newton passes it caused 2-cycles); a root below it raises the gas to
        // e(rho, tfloor) with E's energy, counted as implicit_pos_gas
        Real eg = x.w(m,CLW_EP,k,j,i) - ex;
        // clean_wb_exch_local: the cell's exchange solved exactly from the transported
        // E* = X + den - wk - wp at the solve's opacities, G and s = EP: the root
        // T >= tfloor of (1 + b)(e(rho, T) - eint) + c4 T^4 - b E* - s (b = h c kE,
        // c4 = h c kP G a), increasing in T, so unique; E = (E* - s + c4 T^4)/(1 + b) > 0
        // for E* - s > 0 and the gas takes E* - E; the face FCT's probe keeps the
        // iterate's T
        const Real es = x.w(m,CLW_X,k,j,i) + den - wk - wp;
        const Real sp = x.w(m,CLW_EP,k,j,i);
        const Real lb = x.h*x.c*x.op(m,VET_OP_E,k,j,i);
        if (gtf && lex && !prb) {
          const Real d = uh(m,IDN,k,j,i), e0 = x.w(m,CLW_EI,k,j,i);
          const Real c4 = x.h*x.c*x.op(m,VET_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_;
          auto fx = [&](const Real tt, Real &df) {
            Real ee, p, chi_rho, chi_t, cv;
            eos.ThermoAt(d, tt, ee, p, chi_rho, chi_t, cv);
            df = (1.0 + lb)*d*cv + 4.0*c4*tt*tt*tt;
            return (1.0 + lb)*(ee - e0) + c4*SQR(SQR(tt)) - lb*es - sp;
          };
          Real df, tl = tgf, th = fmax(t, tgf), tr = tgf;
          const bool flo = !(fx(tl, df) < 0.0);
          if (!flo) {
            for (int n = 0; n < 2100 && fx(th, df) < 0.0; ++n) {tl = th; th *= 2.0;}
            tr = th;
            Real fr = fx(tr, df);
            for (int n = 0; n < 200 && fr != 0.0 && th - tl > 4.0e-16*th; ++n) {
              if (fr > 0.0) {th = tr;} else {tl = tr;}
              const Real tn = tr - fr/df;
              tr = (tn > tl && tn < th) ? tn : 0.5*(tl + th);
              fr = fx(tr, df);
            }
          }
          // a root at tfloor: the exchange at T = tfloor; the gas's rest to e(rho,
          // tfloor) is implicit_pos_gas's raise, below
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(d, tr, ei, p, chi_rho, chi_t, cv);
          eg = flo ? (lb*es + sp - c4*SQR(SQR(tgf)))/(1.0 + lb) : ei - e0;
        } else if (gtf) {
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(uh(m,IDN,k,j,i), fmax(t, tgf), ei, p, chi_rho, chi_t, cv);
          eg = ei - x.w(m,CLW_EI,k,j,i);
          if (t < tgf && !prb) {
            Real e0;
            eos.ThermoAt(uh(m,IDN,k,j,i), t, e0, p, chi_rho, chi_t, cv);
            const Real vol = CleanVol(x, m, k, j, i);
            Kokkos::atomic_add(&pc_(VET_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), (ei - e0)*vol);
          }
        }
        // clean_xt_exact: the work the row took but the gas's kinetic energy did not
        // receive goes to E whole, eint the heating only; clean_stag_central: and q X
        const Real dw = xwb ? x.w(m,CLW_XW,k,j,i) - wk - wp + (x.ucp ?
                        x.w(m,CLW_PD,k,j,i)*x.w(m,CLW_X,k,j,i) : 0.0) : 0.0;
        const Real rd = (den - eg - wk - wp - dw)/(gtf ? 1.0 :
                        1.0 + x.h*x.c*x.op(m,VET_OP_E,k,j,i)*x.w(m,CLW_RY,k,j,i)) + dw;
        if (cls) x.w(m,CLW_PW,k,j,i) = fabs(den - eg - wk - wp)/fmax(fabs(e) +
            fabs(x.w(m,CLW_EO,k,j,i)) + fabs(x.w(m,CLW_EI,k,j,i)), 1.0e-300);
        e += rd;
        if (fdn > 0) {dg[0] = eg; dg[1] = wk; dg[2] = wp; dg[3] = dw; dg[4] = rd;}
        if (prb) {   // less E's round-off and e_floor: a limited cell stays above it
          const Real q = x.w(m,CLW_Q,k,j,i);
          wfp(m,3,k,j,i) = e + q - efl -
              1.0e-12*(fabs(e) + fabs(q) + fabs(eg) + fabs(wk) +
              fabs(wp) + fabs(x.w(m,CLW_X,k,j,i)) + fabs(x.w(m,CLW_EO,k,j,i)));
          return;
        }
        den -= rd + wp;
        wsl = wex ? wp : wp - wr;
        // implicit_pos_gas: a gas eint below e_min = max(e(rho, tfloor),
        // implicit_pos_gas_frac x its old eint) is raised to it with energy from the
        // cell's E above e_floor (gas + E kept; counted)
        if (pgas) {
          Real ef = 0.0, pr, chi_rho, chi_t, cv;
          if (eos.tfloor > 0.0) {
            eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
          }
          // from the gas's true eint, the stage start's deficit below 0 included
          const Real ei = x.w(m,CLW_EI,k,j,i), emin = fmax(ef, pgf*fmax(ei, 0.0));
          const Real need = emin - ei - d0 - (den - wk);
          const Real take = fmin(need, e - efl);
          const Real vol = CleanVol(x, m, k, j, i);
          const bool tk = take > 1.0e-12*emin;   // beyond eint's round-off at e_min
          if (tk) {
            den += take;
            e -= take;
            Kokkos::atomic_add(&pc_(VET_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), take*vol);
          }
          // implicit_pos_gas_cons: what E could not cover, from the column, below
          const Real un = need - (tk ? take : 0.0);
          if (pgc && un > 1.0e-12*emin) pgd(m,k,j,i) = un*vol;
        }
      }
      uh(m,IEN,k,j,i) += den;
      x.w(m,CLW_JV,k,j,i) = den - wk;   // its change of the gas eint (cln_neint)
    }
    // (implicit_pos_floor: an E at or below e_floor is
    // raised to it with the energy of the cell's gas above e_min; what the gas cannot
    // cover is taken from the radiation of the cell's column (implicit_pos_floor_cons,
    // below), else counted as created; counted).  A non-finite E (a diverged iterate) is
    // never raised: it stays non-finite, and the check below rejects the solve
    if (dinf && !prb && m == 0 && k == ks && j == js && i == is) e = -INFINITY;
    if (pflr && e <= efl && Kokkos::isfinite(e)) {
      const Real vol = CleanVol(x, m, k, j, i);
      const bool gas = feedback && !x.vzero;
      Real take = 0.0;
      if (gas) {
        Real ef = 0.0, pr, chi_rho, chi_t, cv;
        if (eos.tfloor > 0.0) {
          eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
        }
        const Real ei = x.w(m,CLW_EI,k,j,i), en = ei + d0 + x.w(m,CLW_JV,k,j,i);
        // implicit_pos_floor_gas_keep: the gas keeps pgk of its eint
        const Real emin = fmax(fmax(ef, pgf*fmax(ei, 0.0)), pgk*fmax(en, 0.0));
        take = fmin(efl - e, fmax(en - emin, 0.0));
        den -= take;
        uh(m,IEN,k,j,i) -= take;
        x.w(m,CLW_JV,k,j,i) -= take;
        Kokkos::atomic_add(&pc_(VET_POS_FLR_DE), take*vol);
      }
      if (pfc) {
        pfd(m,k,j,i) = (efl - e - take)*vol;
      } else {
        Kokkos::atomic_add(&pc_(VET_POS_FLR_UN), (efl - e - take)*vol);
      }
      if (fdn > 0 && Kokkos::atomic_fetch_add(&fdc(0), 1.0) < fdn) {
        printf("CLEAN WBFDIAG cycle %d t2s %d m %d kji %d %d %d rho %.4e T %.4e EO %.4e "
               "X %.4e E %.4e Q %.4e EI %.4e d0 %.4e eg %.4e wk %.4e wp %.4e dw %.4e "
               "rd %.4e den %.4e take %.4e dm %.4e %.4e %.4e vol %.4e\n", fcy, fts, m, k,
               j, i, uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), x.w(m,CLW_EO,k,j,i),
               x.w(m,CLW_X,k,j,i), e, x.w(m,CLW_Q,k,j,i), x.w(m,CLW_EI,k,j,i), d0, dg[0],
               dg[1], dg[2], dg[3], dg[4], den, take, dm[0], dm[1], dm[2], vol);
      }
      e = efl;
      Kokkos::atomic_add(&pc_(VET_POS_FLR), 1.0);
    }
    if (slopes) {
      kk(m,VET_T2_E,k,j,i) = (e + wsl - x.w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,VET_T2_M1+n,k,j,i) = dh[n]/h;
      kk(m,VET_T2_EN,k,j,i) = den/h;
    }
    u0_(m,VET_E,k,j,i) = e;
    // the cell flux: the face momenta's cell vector (frame of y on the cubed sphere);
    // Phi's in cln_phic, above
    {
      const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), multi_d ? f2(m,k,j,i) : 0.0,
                         multi_d ? f2(m,k,j+1,i) : 0.0, three_d ? f3(m,k,j,i) : 0.0,
                         three_d ? f3(m,k+1,j,i) : 0.0};
      Real fv[3];
      CleanStagCellVector(x, m, k, j, i, c, fv);
      CleanStagHydro(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) {
        u0_(m,VET_F1+d,k,j,i) = fv[d];
        if (!xh.hi) phq(m,d,k,j,i) = fv[d];   // no split: Phi is the face momentum
      }
    }
  };
  if (wf) {
    CleanFor("cln_wbprb", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      wbc(m, k, j, i, true);
    });
    CleanWbFct(x, px, wfp, pfo, pf);
  }
  if (wmc) Kokkos::deep_copy(xh.mclc, 0.0);   // (the write-back's own, not the probe's)
  CleanFor("cln_wbcell", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    wbc(m, k, j, i, false);
  });
  if (wmc) {
    auto hm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), xh.mclc);
    Real r[4];
    for (int n = 0; n < 4; ++n) r[n] = CleanReduce(hm(n), false);
    if (global_variable::my_rank == 0 && r[0] > 0.0) {
      std::printf("CLEAN MCAP cycle %d t2s %d: %.0f face momentum terms at the cap in the "
                  "write-back, gas momentum withheld sum dp V (x1 x2 x3) %.6e %.6e %.6e\n",
                  pmy_pack->pmesh->ncycle, t2_solve, r[0], r[1], r[2], r[3]);
    }
  }
  if (lzb) {   // (the L_z ledger: the stage's L_z as its slope, then booked)
    const Real s = CleanSumCells(lzw);
    if (slopes) (t2_solve == VET_T2S_STAGE1 ? cln_lzs2 : cln_lzs1) = s/h;
    CleanLzBook(s, t2_solve == VET_T2S_STAGE2);
  }
  // xt_exact: the stage's closure |dE + de + dK + h div Phi| over E + E_old + eint, max
  // over the cells
  if (cls) {
    auto w = cln_w;
    const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
    Real cm = 0.0;
    Kokkos::parallel_reduce("cln_wbclose", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &c) {
      const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
      const Real v = w(m,CLW_PW,k,j,is + idx%n1);
      c = fmax(c, Kokkos::isfinite(v) ? v : INFINITY);   // Max drops a NaN
    }, Kokkos::Max<Real>(cm));
    cm = CleanReduce(cm, true);
    if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(cm))) {
      std::printf("CLEAN %sXT cycle %d t2s %d write-back closure %.3e\n",
                  std::isfinite(cm) ? "" : "WARNING ",
                  pmy_pack->pmesh->ncycle, t2_solve, cm);
    }
  }
  // implicit_pos_floor_cons: the uncovered raise of a cell from the radiation above
  // e_floor of the nearest cells of its radial column, i -+ d, split by their energy
  // above the floor, with their stage slopes.  clean_stag_colflux: the stored Phi carries
  // each move as x1 face fluxes, so time2_face_fct limits the moves.
  const bool pfk = pfo && cln_colf;
  if (pfc) {
    par_for("cln_pos_flr_cons", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      CLEAN_HIPX(const CleanArgs &x = *px;)
      auto vol = [&](const int q) {
        return CleanVol(x, m, k, j, q);
      };
      Real mv = 0.0, un = 0.0;
      for (int i = is; i <= ie; ++i) {
        Real def = pfd(m,k,j,i);
        for (int d = 1; d <= ie - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(u0_(m,VET_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,VET_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,VET_E,k,j,q) = fmax(u0_(m,VET_E,k,j,q) - de, efl);
            if (slopes) kk(m,VET_T2_E,k,j,q) -= de/h;
            for (int f = (sd ? i + 1 : qa + 1); pfk && f <= (sd ? qb : i); ++f) {
              pf(m,0,k,j,f) += (sd ? -1.0 : 1.0)*de*vol(q)/
                               (h*CleanDivW(x, 0, m, k, j, f, 0)*vol(f));
            }
          }
          def -= tk;
          mv += tk;
        }
        un += fmax(def, 0.0);
      }
      if (mv > 0.0) Kokkos::atomic_add(&pc_(VET_POS_FLR_MV), mv);
      if (un > 0.0) Kokkos::atomic_add(&pc_(VET_POS_FLR_UN), un);
      if (!pgc) return;
      // implicit_pos_gas_cons: the gas raise E could not cover, the same way, into the
      // cell's gas; what the column cannot give stays below e_min (the EOS floor ledger)
      Real gm = 0.0;
      for (int i = is; i <= ie; ++i) {
        Real def = pgd(m,k,j,i);
        for (int d = 1; d <= ie - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(u0_(m,VET_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,VET_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,VET_E,k,j,q) = fmax(u0_(m,VET_E,k,j,q) - de, efl);
            if (slopes) kk(m,VET_T2_E,k,j,q) -= de/h;
            for (int f = (sd ? i + 1 : qa + 1); pfk && f <= (sd ? qb : i); ++f) {
              pf(m,0,k,j,f) += (sd ? -1.0 : 1.0)*de*vol(q)/
                               (h*CleanDivW(x, 0, m, k, j, f, 0)*vol(f));
            }
          }
          const Real dg = tk/vol(i);
          uh(m,IEN,k,j,i) += dg;
          x.w(m,CLW_JV,k,j,i) += dg;
          if (slopes) kk(m,VET_T2_EN,k,j,i) += dg/h;
          def -= tk;
          gm += tk;
        }
      }
      if (gm > 0.0) Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), gm);
    });
  }
  if (pflr) {
    Real un = 0.0;
    Kokkos::deep_copy(un, Kokkos::subview(pc_, VET_POS_FLR_UN));
    un = CleanReduce(un - un0, false);
    cln_flrun = un;
    cln_flrid = 4*pmy_pack->pmesh->ncycle + t2_solve;
    if (un > 0.0 && global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d t2s %d: implicit_pos_floor created %.6e of E "
                  "(neither the gas nor the column covered the raise)\n",
                  pmy_pack->pmesh->ncycle, t2_solve, un);
    }
  }
  // cells the write-back left with gas eint < 0 or below the EOS floor (the solve is then
  // not accepted; an eint the hydro left there is not the solve's)
  cln_neint = 0;
  if (!feedback) return;
  const bool cry = cln_eicry;   // a carried deficit left below: its E could not pay
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nx1 = id.nx1, nx2 = id.nx2, nx3 = id.nx3, nji = nx2*nx1, nkji = nx3*nji;
  Real nb = 0.0, nh = 0.0;
  Kokkos::parallel_reduce("cln_neint", CleanRp256(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n, Real &nf) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji + ks;
    const int j = (idx%nji)/nx1 + js, i = idx%nx1 + is;
    Real p[3];
    CleanGasMomentum(x, uh, m, k, j, i, p);
    const Real ke = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/uh(m,IDN,k,j,i);
    // and below the eint of the EOS temperature floor (the hydro floor would add energy
    // the radiation's budget does not see)
    Real ef = 0.0;
    if (eos.tfloor > 0.0) {
      Real pr, chi_rho, chi_t, cv;
      eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
    }
    const Real tz = 1.0e-12*(fabs(uh(m,IEN,k,j,i)) + ke + ef);   // eint's round-off
    // and any non-finite E, F or gas state, also at the floor: rejected and retried,
    // never accepted
    const bool fin = Kokkos::isfinite(u0_(m,VET_E,k,j,i) + u0_(m,VET_F1,k,j,i) +
                                      u0_(m,VET_F1+1,k,j,i) + u0_(m,VET_F1+2,k,j,i) +
                                      uh(m,IDN,k,j,i) + uh(m,IEN,k,j,i) + ke);
    const bool low = !(CleanGasEint(x, uh, ke, m, k, j, i) > ef - tz);
    if ((low && !(x.w(m,CLW_JV,k,j,i) >= 0.0) && x.w(m,CLW_EI,k,j,i) > ef + tz) ||
        (low && cry && x.w(m,CLW_EI,k,j,i) < 0.0) ||
        u0_(m,VET_E,k,j,i) < 0.0 || !fin) n += 1.0;
      // without implicit_pos_gas: an eint the hydro's floor will raise, counted
    if (low && !pgas) nf += 1.0;
  }, Kokkos::Sum<Real>(nb), Kokkos::Sum<Real>(nh));
  cln_neint = static_cast<int>(CleanReduce(nb, false));
  nh = CleanReduce(nh, false);
  if (nh > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: %.0f cells left with gas eint below the "
                "EOS floor (implicit_pos_gas off: the hydro's floor raises them, energy "
                "created)\n", pmy_pack->pmesh->ncycle, t2_solve, nh);
  }
}
#endif  // VET_TU_PART(2)

}  // namespace radvet
