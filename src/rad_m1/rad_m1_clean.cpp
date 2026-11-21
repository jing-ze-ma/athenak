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
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
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
#include "utils/scratch_view.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
// M1_TU_PART: this file is compiled as 3 translation units (rad_m1_clean.cpp = part 0, rad_m1_clean_p<k>.cpp
// = part k, see split_tu.py); every member function below is in exactly one part.
#ifndef M1_TU_PART_N
#define M1_TU_PART_N 0
#endif
#define M1_TU_PART(k) (M1_TU_PART_N == (k))

namespace radm1 {
extern bool g_vra_adv;   // vet_rotate_advance (rad_m1_vetgd.cpp)
// the policy of the heavy per-cell reductions: Kokkos' HIP backend compiles a plain
// RangePolicy reduction for up to 1024 threads per block, which caps it at 128 VGPRs on
// gfx942 and spilled them (cln_velocity 1005 VGPRs); 256 threads per block as athena's
// par_for and old production's m1-fast4.  Max reductions and counts: bitwise.
using CleanRp256 = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 1>>;
// (the heavy per-cell passes at 2 blocks of 256 per SM, <= 128 registers, as the old
// production's par_for_lb (m1-fast4): 164-213 registers left 1 block per SM; the
// arithmetic is unchanged)
// (CUDA only: on HIP (wave64) the bounds stay at 1 block, as measured on MI300A, but
// for cln_rhs, CleanParForLb)
#if defined(KOKKOS_ENABLE_CUDA)
static constexpr int CleanLbMb(int mb) {return mb;}
#else
static constexpr int CleanLbMb(int) {return 1;}
#endif
using CleanRpLb2 =
    Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, CleanLbMb(2)>>;
// (the 4-D par_for, same flattening and index arithmetic, at MB blocks of 256 per SM)
template <int MB, typename F>
static void CleanParForLb(const char *name, const int nmb1, const int ks, const int ke,
                          const int js, const int je, const int is, const int ie,
                          const F &f) {
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int mb = 2;   // (cln_rhs at 256 VGPRs, 2 waves per SIMD: MI300A RSG -45 ms
                          // per 8 cycles, though it spills more; CleanRpLb2's passes at 2
                          // measured no gain)
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
  // (_Exit: no static destructors behind live Kokkos views, as VgdFatalAt)
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

//! \brief staggered: v of CleanStagGasMomentum's components in the hydro's (the cubed
//! sphere: from the frame of y; else unchanged)

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
//! CleanStagGasMomentum: the sum of the cell's shares of its faces' absorption.  On the
//! cubed sphere the faces' changes (along their normals) give the vector of
//! CleanStagCsVector of weights 1/2, their sum as on spherical polar.  Its lateral faces
//! are no phi faces, so L_z does not telescope: dlz, minus the L_z V R (dF_phi/c^2 +
//! dp_phi) of the radiation's cell vector F (CleanStagCsVector of the face momenta) and
//! the gas, is booked as the radiation's (CleanLzBook; the stress torques through the
//! inner faces cancel in its sum).  3-D spherical polar: dlz the same of the phi faces'
//! F; inside they telescope, so its sum is the L_z that crossed the boundaries, at an open
//! or imposed radial boundary the one-sided torque form of the transverse stress
//! (CleanFaceStress) less the F solve's dissipation there.

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
//! \brief cubed sphere, clean_cs_lz_tel: the L_z the gas of cell (k,j,i) takes over the
//! stage in telescoping form, and in rl the change of the radiation's (its faces' momenta:
//! the cell's halves of V F_f lambda_f/c^2, lambda_f = (M_f x n_f)_z/A_f of the lateral
//! faces, 0 of the radial ones).  The gas takes the faces' absorption as their torque,
//! V cf lambda, and in place of the faces' own stress torque, 1/2 V lambda (T - dF/c^2),
//! the exact one of P = D E over the cell: the angular-momentum fluxes r_f A (n x P n)_z
//! of its radial faces and (M_f x P n)_z of its lateral ones, P the two cells' mean (one
//! sided at a radial boundary).  The interior faces cancel in the sum, so gas and
//! radiation L_z change only by the radial boundaries' torque (CleanLzBook).  In xs the
//! stress part, the exact torque less the faces' own: its net the gas's (CleanWriteBack)

KOKKOS_INLINE_FUNCTION
static Real CleanCsTorque(const CleanArgs &x, const int e, const int m, const int k,
                          const int j, const int i, Real &rl, Real &xs) {
  const Real v = x.vol(m,k,j,i), ic2 = 1.0/(x.c*x.c);
  Real gl = 0.0, gs = 0.0, tc = 0.0, p0[6];
  rl = 0.0;
  for (int q = 1; q < 3; ++q) {   // (the lateral faces' absorption and momenta)
    Real cf[2], wf[2], tq[2];
    CleanMomentumChange(x, e, q, m, k, j, i, nullptr, cf, wf, tq);
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1);
      const Real f = x.fs ? x.fc[q](m,CLF_N,kf,jf,i) : CleanFaceFlux(x, e, q, m, kf, jf, i, true);
      const Real df = (f - x.fc[q](m,CLF_G,kf,jf,i))*ic2;
      const Real lam = CleanCsArm(x, q, m, kf, jf, i);
      gl += v*lam*cf[s];
      gs += 0.5*v*lam*(tq[s] - df);
      rl += 0.5*v*lam*df;
    }
  }
  for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
  for (int q = 0; q < 3; ++q) {   // (the stress torque through the cell's six faces)
    const int dk = (q == 2), dj = (q == 1), di = (q == 0);
    for (int s = 0; s < 2; ++s) {
      const int o = 2*s - 1, kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
      const int ko = k + o*dk, jo = j + o*dj, io = i + o*di;
      const bool bnd = CleanBoundaryFace(x, q, m, kf, jf, iff);
      Real pf[6], a[3], t[2];
      for (int n = 0; n < 6; ++n) {
        pf[n] = bnd ? p0[n] : 0.5*(p0[n] + x.cs(m,CLN_D11+n,ko,jo,io)*x.w(m,e,ko,jo,io));
      }
      for (int c = 0; c < 3; ++c) {
        a[c] = (q == 0) ? x.x1f(m,iff)*x.area[0](m,kf,jf,iff)*x.g.At<CLG_N>(m,c,kf,jf,iff)
                        : x.g.At<CLG_M>(m,3*(q-1)+c,kf,jf,iff);
      }
      for (int b = 0; b < 2; ++b) {   // (P n, outward: its x and y components)
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
//! \brief CleanVelocityUpdate's cell (k,j,i): the momentum change (CLW_DM), with fa the
//! advective coefficients of the new velocity (CLW_A1), and the maxima of |dv|, |v| and
//! of the force scale; kFs: the F solve's F known set (xh.fs, xh.fa)

template <bool kFs>
KOKKOS_INLINE_FUNCTION
static void CleanVelocityCell(const CleanArgs &x, const CleanArgs &xh, const CleanArgs &xa,
                              const DvceArray5D<Real> &uh,
                              const int fref, const DvceArray5D<Real> &fr,
                              const DvceArray4D<Real> &aref, const bool cap,
                              const DvceArray5D<Real> &op, const Real arad_, const bool fa,
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
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPreflight
//! \brief the production keys of <rad_m1> (and <time>) an input names with scheme =
//! clean (not the defaults a restart file lists): a key the clean solve replaces names
//! its clean counterpart, the others it does not read are listed once, and a key whose
//! loss would change the physics stops the run unless it is off.  The keys the clean
//! solve reads (boundaries, beams, the time step, multi-rate, the hesdirk2 tableau and
//! its backward-Euler steps) pass.

#if M1_TU_PART(0)
void RadiationM1::CleanPreflight(ParameterInput *pin) {
  // block/key (* ends a prefix): nullptr read, "" no effect, "!" fatal, else clean key
  static const char *tab[][2] = {
    {"rad_m1/implicit_bc_*", nullptr}, {"rad_m1/implicit_flux_x1m*", nullptr},
    {"rad_m1/implicit_ebath_*", nullptr}, {"rad_m1/implicit_beam*", nullptr},
    {"rad_m1/implicit_cfl", nullptr}, {"rad_m1/implicit_allow_multid", nullptr},
    {"rad_m1/implicit_mr_*", nullptr}, {"rad_m1/implicit_det_trace", nullptr},
    {"rad_m1/implicit_tol", "clean_tol_nl"}, {"rad_m1/implicit_lin_tol", "clean_tol"},
    {"rad_m1/implicit_maxit", nullptr},   // (read by the old module's ImplicitInit; a
    {"rad_m1/implicit_tol_rel_max", nullptr},   // named other value than the clean key's
    {"rad_m1/implicit_lin_maxit", nullptr},   // prints a NOTE, rad_m1.cpp)
    {"rad_m1/implicit_recon", "clean_recon"},
    {"rad_m1/implicit_closure_lag", "clean_closure"},
    {"rad_m1/time2_vet", "clean_closure"}, {"rad_m1/time2_vet_col", "clean_closure"},
    {"rad_m1/time2_vet_sc", "clean_closure"},
    {"rad_m1/time2_vet_extrap", "clean_closure"}, {"rad_m1/time2_lin_tol*", "clean_tol"},
    {"rad_m1/implicit_ghost_x1min", "!"}, {"rad_m1/ic_rad*", "clean_init"},
    {"rad_m1/implicit_predictor", nullptr}, {"rad_m1/implicit_predictor_tcap", nullptr},
    {"rad_m1/implicit_resid_fatal", nullptr},
    // (members of ImplicitInit the clean steps use: the enthalpy flux, the conservative gas
    // source, the Marshak face, the predictor order, the hold, the fused kernels, the floor
    // solve's pairing with implicit_pos_floor)
    {"rad_m1/implicit_enthalpy", nullptr}, {"rad_m1/implicit_gas_src", nullptr},
    {"rad_m1/implicit_marshak_face", nullptr}, {"rad_m1/implicit_predictor_order", nullptr},
    {"rad_m1/implicit_hold", nullptr}, {"rad_m1/implicit_fast_kernels", nullptr},
    {"rad_m1/implicit_pos_floor_solve", nullptr},
    {"rad_m1/implicit_pos_gas*", nullptr}, {"rad_m1/implicit_pos_floor", nullptr},
    {"rad_m1/implicit_pos_floor_cons", nullptr},
    {"rad_m1/implicit_pos_floor_gas_keep", nullptr}, {"rad_m1/implicit_pos_gas_cons", nullptr},
    {"rad_m1/implicit_*", ""},
    // (the hesdirk2 stage formation the clean solve shares, Time2FormStage / Time2Stash,
    // and dt_retry: read; ImplicitSolve's acceptance and posfix and the KE cap: none)
    {"rad_m1/time2_k1_*", nullptr}, {"rad_m1/time2_face_fct*", nullptr},
    {"time/dt_retry_t2_be", nullptr},
    {"rad_m1/time2_accept_fallback", ""}, {"rad_m1/time2_posfix*", ""},
    {"rad_m1/time2_inc_kecap", ""},
    // (the implicit solve's admissibility and dt_retry position limiters, its Zalesak
    // round cap, the vet_col lateral posfix and the x2 vet fields: none under clean)
    {"rad_m1/time2_adm_poslim", ""}, {"rad_m1/time2_dtr_poslim", ""},
    {"rad_m1/zalesak_cap_rounds", ""}, {"rad_m1/vet_col_lat_posfix", ""},
    {"rad_m1/time2_vet_x2_fields", ""}};
  const bool say = (global_variable::my_rank == 0);
  std::string none;
  std::vector<std::string> drop;
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
        if (c != "!" && k.compare(0, 16, "rad_m1/implicit_") == 0) {
          drop.push_back(ln.param_name);
        }
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
  // (the old module's implicit_* keys the clean solve does not read leave the input of a
  // fresh run, so that its ImplicitInit parses its defaults and none of its cross-checks
  // sees them; a restart keeps its file's keys, which also set the layout of its data)
  for (InputBlock &ib : pin->block) {
    if (ib.block_name != "rad_m1" || global_variable::restart_run) continue;
    ib.line.remove_if([&](const InputLine &ln) {
      return std::find(drop.begin(), drop.end(), ln.param_name) != drop.end();
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanInit
//! \brief allocates the arrays of the clean solve and checks that the setup is supported

void RadiationM1::CleanInit() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  // (the L_z ledger follows the gas momentum's increments as the slopes' sums: not the
  // cell-wise scalings of these)
  if (cln_lzbook && (t2_theta > 0.0 || t2_inckc > 0.0 || t2_flr_rho > 0.0 || cln_kreal ||
                       cln_kpos || cln_cell)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<rad_m1> the L_z ledger (cubed sphere, 3-D spherical polar) needs no "
              << "time2_stage_theta, time2_inc_kecap, "
              << "time2_floor_eint, clean_stag_kreal, clean_stag_kpos (staggered)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  if (!cln_cell && cln_phic.extent_int(0) != nmb) {
    Kokkos::realloc(cln_phic, nmb, 3, n3, n2, n1);
  }
  // (staggered on the cubed sphere: the cell terms and faces of CleanStagFSolve)
  const bool stcs = pmy_pack->pmesh->use_cubed_sphere && !cln_cell;
  Kokkos::realloc(cln_v, nmb, stcs ? 18 : 3, n3, n2, n1);
  // (mem_arena, staggered, closure outside the solve: the shared slots too, after
  // cln_w.o, CleanWoBind; clean_closure = iterate | extrap* builds inside the solve)
  cln_w.s = DvceArray5D<Real>();
  if (cln_cell || !marena_on || cln_closit || cln_extrap > 0) {
    Kokkos::realloc(cln_w.s, nmb, CLW_NS, n3, n2, n1);
  }
  const int nwo = cln_cell ? CLW_NC : (cln_kreal ? CLW_NR : (cln_kpos ? CLW_NKP :
                  ((cln_xte > 0) ? CLW_NX : ((cln_kpos || cln_xpat) ? CLW_NK : CLW_N))));
  // (mem_arena, staggered: made at the first solve, over the intensities, CleanWoBind)
  cln_wo_nc = nwo - CLW_NS;
  cln_w.o = DvceArray5D<Real>();
  if (cln_cell || !marena_on) Kokkos::realloc(cln_w.o, nmb, cln_wo_nc, n3, n2, n1);
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
  // (the staggered bracket's lagged flux DG in one more slot, CleanStagBracket)
  cln_fdg = (!cln_cell && cln_consist) ? nfc : -1;
  cln_xqs = cln_xd ? nfc + (cln_fdg >= 0) : -1;   // (clean_stag_xdonor's q)
  cln_kqs = cln_kpos ? nfc + (cln_fdg >= 0) + (cln_xd ? 5 : 0) : -1;   // (kpos: q)
  const int nfx = nfc + (cln_fdg >= 0) + (cln_xd ? 5 : 0) + (cln_kpos ? 2 : 0);
  if (cln_kpos || cln_xpat) {
    Kokkos::realloc(cln_kcnt, 2);
    Kokkos::deep_copy(cln_kcnt, 0.0);
  }
  Kokkos::realloc(cln_fc1, nmb, nfx, n3, n2, n1);
  Kokkos::realloc(cln_fc2, nmb, nfx, n3, n2, n1);
  Kokkos::realloc(cln_fc3, nmb, nfx, n3, n2, n1);
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
    Kokkos::realloc(cln_fref, nmb, (force_ref == M1_FREF_WB_MD) ? 9 : 3 + cln_lzt, n3, n2,
                    n1);
  }
  if (cln_stst) Kokkos::realloc(cln_st, nmb, CLN_ST, indcs.nx3, indcs.nx2, indcs.nx1);
  if (cln_stst && pmy_pack->pmesh->use_cubed_sphere) {
    Kokkos::realloc(cln_csy, nmb, CLN_CSY, indcs.nx3, indcs.nx2, indcs.nx1);
  }
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
    // (vet_rotate_avg_ema: production's EMA, only on clean_rot_ring = live; without
    // _dema, _tau, _lc and _diag, which CleanRingEma does not carry)
    const bool ema_ok = cln_ring == 1 && !vra_ema_dema && vra_ema_tau <= 0.0 &&
                        vra_ema_lc <= 0.0 && vra_ema_diag == 0;
    if ((vra_ema && !ema_ok) || vra_eq ||
        (vgd_xf && !(vgd_rot_dt > 0.0 || vgd_rot_auto)) ||
        (cln_rfv == 2 && !fixed) || (vgd_rot3 && !pm->three_d)) {
      CleanFatal("<rad_m1>/scheme = clean with a rotating vet_gd set: vet_rotate_avg = K "
                 "(J-weighted; _ema only with clean_rot_ring = live, no _dema, _tau, "
                 "_lc, _diag), vet_gd_rotate_xfade with vet_gd_rotate_dt, clean_rays = "
                 "fv only for a fixed union (vet_rotate_avg without vet_rotate_advance), "
                 "and on a 2-D axisymmetric mesh the z axis");
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
    // (the rays' theta walls: rad_m1.cpp, vet_gd_theta_wall or, unnamed, the mesh's walls)
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
  // (all CLG_NV in gn first: the exchanges, the distances; then into cln_gn's two arrays)
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
  // the directions of CleanGeo::a are the same in every cell of an x1 range of a line,
  // to the bit
  const int is = indcs.is, ie = indcs.ie;
  Kokkos::realloc(cln_gn.r, nmb1 + 1, CLG_NV - CLG_NA, n3, n2, n1);
  Kokkos::realloc(cln_gn.a, nmb1 + 1, CLG_NA, n3, n2, 3);
  cln_gn.is = is;
  cln_gn.ie = ie;
  auto gr = cln_gn.r;
  auto gt = cln_gn.a;
  int nbad = 0;
  Kokkos::parallel_reduce("cln_csgeo2", Kokkos::MDRangePolicy<Kokkos::Rank<4>, DevExeSpace>(
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
    CleanFatal("<rad_m1> clean cubed sphere: " + std::to_string(nbad) +
               " geometry directions change along x1 (CleanGeo)");
  }
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
  cln_insteady = true;
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
  cln_insteady = false;
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
    std::printf("<rad_m1> clean_init steady: %s after %d passes, residual %.3e "
                "(tol %.1e; %.0f formal solutions); radiation E %.10e -> %.10e, F/c "
                "%.4e %.4e %.4e -> %.4e %.4e %.4e (volume sums)\n",
                ok ? "converged" : "NOT converged",
                std::min(n + 1, cln_steady_n), d, cln_steady_tol, vet_ncall - nc0, t0[0],
                t1[0], t0[1]/cl, t0[2]/cl, t0[3]/cl, t1[1]/cl, t1[2]/cl, t1[3]/cl);
  }
  // a last change that is not small is a diverged start, not an unconverged one
  if (!(d < 1.0e-2)) {
    CleanFatal("clean_init steady: the last pass changed E or F/c by " +
               std::to_string(d) + " of max E (diverged: clean_init_dtfac, _n)");
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
      if (cln_cell && cln_fvi.extent_int(1) != cln_fvn) {   // (staggered: cln_fh)
        Kokkos::realloc(cln_fvi, pmy_pack->nmb_thispack, cln_fvn, cln_s.extent_int(2),
                        cln_s.extent_int(3), cln_s.extent_int(4));
      }
    }
  }
  if (cln_cssc && !cln_ring) {   // (clean_rays = sc: the sweep's set, its cln_fvi)
    cln_fvd = vgd_dir;
    cln_fvn = vgd_n;
  }
  // (a restart from a file with the intensities (version 5), staggered: the face half
  // ranges they give, set by set as the builds fold them)
  if (cln_cssc && !cln_cell && cln_fh.extent_int(0) == 0 && cln_fvi.extent_int(0) > 0) {
    for (int q = 0; q < cln_fvi.extent_int(1)/vgd_n; ++q) {
      CleanFaceHalf(cln_fvd, q*vgd_n, q == 0, q*vgd_n);
    }
    cln_fvi = DvceArray5D<Real>();
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
  CleanWoBind();
  CleanStagKReal();   // (clean_stag_kreal: the kick E can pay, the rest deferred)
  cln_xdh = false;   // (clean_stag_xdonor_freeze: no hold before this solve's passes)
  CleanOldVector();
  CleanStagKPos();   // (clean_stag_kpos: the increment in Patankar form)
  if (cell) {
    CleanCellOldVector();
    CleanCellFaceCoefficients();
  } else {
    CleanFaceCoefficients();
    CleanStagKDelta();   // (clean_stag_kpos = delta: the faces' rest dG)
    // (the donors' E_old, memory; the bracket's gate q)
    if (cln_gflux || cln_consist) CleanGhostFill(cln_w, CLW_EO, 1);
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
  bool all_low = false;
  cln_all_low = false;
  cln_ncwin = 0;
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
  Real dnl = 0.0, dnl0 = -1.0;   // (the last pass's change and the one before)
  for (int ipos = 0; ; ++ipos) {
    cln_pdfst = true;   // (clean_pass_damp: no history before a solve's first pass)
    CleanFirstIterate();
    cln_fsk = false;
    cln_csyok = false;   // (spdb-1119: the lateral sums' coefficients once per closure)
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
    // (clean_stag_xdonor: q of the first iterate, with its drag share b)
    // (after the first drag: CLW_L holds E_k's lateral sums, CleanStagDrag)
    CleanStagXDonor(true, !(cln_vimp && cln_feedback));
    converged = !iterate;
    npic = 0;
    Real de0 = 1.0e300, dv0 = 1.0e300, trel = 0.0;
    for (int inl = 0; inl < maxnl; ++inl) {
      if (lagged && inl < cln_maxpic) ++npic;
      const bool clos = (cln_closit && inl < cln_maxpic) || (clos1 && inl == 1);
      if (clos) {
        cln_csyok = false;
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
        // (clean_stag_xdonor: q of the iterate; CLW_L holds E_k's lateral sums, those of
        // the velocity update, but where the closure moved D or the bracket took CLW_L)
        if (inl > 0) CleanStagXDonor(false, clos || !cln_feedback,
                                     cln_xdfz > 0 && inl >= cln_xdfz);
        // (after the xdonor: the lateral sums of E_k are its, CleanStagXDonor)
        CleanRightHandSide(false, !(inl > 0 && cln_xd_on));
        // (clean_stag_sbimp: the first iterate's T in exchange with E_k and the beta
        // terms sb(E_k) the RHS just built, then the RHS of that T)
        // (clean_xt_cons: the same, the first iterate's T from the pass's gas row)
        if (inl == 0 && (cln_sbimp || cln_xtcons) && cln_feedback) {
          CleanNewtonUpdate();
          cln_pdfst = true;   // (clean_pass_damp: that update's X is the last solve's)
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
      if (cln_xsw && !cell && nit >= 0) CleanStagXSweep();
      if (cln_xd && !cell && cln_nlprint) CleanStagXCount(inl, nit);   // (diagnostic)
      nit_sum += std::abs(nit);
      ++nsolve;
      cln_nsolve += 1.0;
      cln_nit += std::abs(nit);
      if (nit < 0) dnl = std::numeric_limits<Real>::infinity();
      // (implicit_lin_maxit named, as the old production: a linear solve that ends at
      // its cap leaves its iterate to the pass, which is then not the last one)
      const bool lp = cln_linpass && !cell && nit < 0 && inl + 1 < maxnl;
      if (!iterate || (nit < 0 && !lp)) break;
      const bool fbl = lp;
      if (!cell && !t2stage && inl + 1 == maxnl) {   // (CleanWorstCell)
        auto sv = Kokkos::subview(cln_w.s, Kokkos::ALL, std::make_pair(CLW_EK,
                                  CLW_EK + 1), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        Kokkos::realloc(cln_wci, sv.extent(0), 1, sv.extent(2), sv.extent(3),
                        sv.extent(4));
        Kokkos::deep_copy(cln_wci, sv);
      }
      const Real de = cell ? CleanCellNewtonUpdate() : CleanNewtonUpdate();
      if (cln_bon) tol_nl = cln_tolb;
      if (lagged && inl < cln_maxpic) dpic = de;
      Real dv = 0.0;
      if (cln_feedback) dv = cell ? CleanCellVelocityUpdate() : CleanVelocityUpdate();
      dnl0 = (inl + 1 == maxnl) ? dnl : -1.0;
      dnl = fmax(de, dv);
      // (clean_tol_rel, as the old production's implicit_tol_rel: the strict test stops at
      // kappa times the first pass's change, at most clean_tol_rel_max)
      if (cln_tolrel > 0.0 && !cln_bon && inl == 0) trel = fmin(cln_tolrel*dnl, cln_tolrmx);
      tol_nl = fmax(tol_nl, trel);
      if (cln_nlprint && global_variable::my_rank == 0) {
        std::printf("CLEAN pass cycle %d t2s %d pos %d pass %d it %d de %.3e dv %.3e\n",
                    pmy_pack->pmesh->ncycle, t2_solve, ipos, inl, nit, de, dv);
      }
      if (fbl) {
        de0 = de;
        dv0 = dv;
        continue;
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
    // (clean_stag_xpat: the counts of the passes)
    if (cln_xpat && !cell) CleanStagKEnd();
    // (clean_stag_kpos = redo: X <= 0 left, the stage again with those donors' faces
    // moved)
    if (!cell && CleanStagKRedo()) continue;
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
  if (!converged || (pmy_pack->pmesh->ncycle == dbg_nonconv && dtr_ncut == 0)) {
    nit = -std::abs(nit) - 1;
  }

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
    if (!cell) CleanWorstCell(dnl, dnl0);
    if (fatal) CleanFatal("<rad_m1> clean: a backward-Euler solve did not converge and "
                          "its last change exceeds <rad_m1>/implicit_resid_fatal");
  }
  if (cln_rtb && dtr_redo && cln_bon && nit >= 0 && global_variable::my_rank == 0) {
    std::printf("CLEAN RETRY cycle %d t2s %d: the redone solve converged to the budget "
                "tolerance %.3e, last change %.3e\n", pmy_pack->pmesh->ncycle, t2_solve,
                cln_tolb, dnl);
  }
  if (global_variable::my_rank == 0) {
    if (cell) {
      std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pic %d %.2e pos %d mark %d "
                  "all %d neg %d nreal %d guard %d solves %d\n", pmy_pack->pmesh->ncycle,
                  t2_solve, nit, nit_sum, npic, dpic, npos, nmark, all_low ? 1 : 0, nneg,
                  cln_nnreal, cln_nguard, nsolve);
    } else {
      std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pos %d mark %d all %d neg %d "
                  "r1 %d r2 %d first %d %d guard %d solves %d cwin %d\n",
                  pmy_pack->pmesh->ncycle, t2_solve, nit, nit_sum, npos, nmark,
                  all_low ? 1 : 0, nneg, cln_nr1, cln_nr2, nfirst[0], nfirst[1],
                  cln_nguard, nsolve, cln_ncwin);
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
    // (as the old production with implicit_pos_floor(_cons): an E < 0 that the
    // conservative write-back floor raises does not fail the stage; the floor is counted)
    const bool flr = impl_pos_floor && impl_pos_flr_cons;
    t2_fail = (nit < 0) || (nneg > 0 && !flr) || (cln_feedback && cln_neint > 0) ||
              (t2_solve == M1_T2S_STAGE1 && pmy_pack->pmesh->ncycle == t2_dbg_fail);
    // time2_reject = retry: the cycle again at a smaller dt instead of the backward-Euler
    // redo; after dt_retry_t2_be such cuts the redo (as ImplicitSolve's t2retry-1020)
    const bool rt = t2_fail && t2_reject == 1 && dtr_on && !t2_inprobe && !t2_pfix_want;
    if (rt && !dtr_t2be) {
      if (!dtr_fail) {dtr_why = 5; dtr_resid = -1.0;}
      dtr_fail = true;
    }
    if (rt && dtr_t2be) dtr_t2hit = true;
    if (cln_lzbook && cln_feedback && t2_solve == M1_T2S_STAGE2 && !t2_fail) {
      CleanLzBook(0.0, false);   // (the accepted stage 2's L_z, pending until here)
    }
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

  cln_xd_on = false;   // (clean_stag_xdonor: q is of this solve's passes only)
  cln_kq_on = false;   // (clean_stag_kpos: q, pd of this stage)

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
    CleanFor("cln_f0", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
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
  if (cln_fvgf) {   // (the ghosts of the restart's cln_fvi, CleanFvRstCopy)
    CleanGhostFill(cln_fvi, 0, cln_fvi.extent_int(1));
    cln_fvgf = false;
  }
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
    Real fv[3];
    CleanStagCellVector(x, m, k, j, i, c, fv);
    CleanStagHydro(x, m, k, j, i, fv);
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
//! \fn void RadiationM1::CleanStagKReal
//! \brief clean_stag_kreal (staggered, hesdirk2 stage solves): the stage increment's gas
//! momentum kick D (t2inc) gives the gas the kinetic energy W(1), W(s) = KE(m_s + s D) -
//! KE(m_s) = a s + b s^2 (a = p_s.p_D/rho, b = |p_D|^2/(2 rho), p the momentum in the
//! frame of y, CleanGasVector), which its slopes took from E.  A cell whose E cannot pay
//! it, W(1) > (1 - th) E_a (E_a = E_old + W(1), th = CLN_KPOS_TH), takes the kick at the
//! largest s in [0, 1] with W(s) <= (1 - th) max(E_a, 0) (closed form; s = 0 is always
//! feasible): the old vector's momentum m_s + s D, its E_old + W(1) - W(s) >= th E_a and
//! its gas energy less the same, so E + gas and the gas eint are as formed; the rest
//! (1 - s) D (CLW_KR, frame of y) is the stage's: the momentum change of the iterate and
//! the write-back take it, the E row its work (CleanRightHandSide, CleanWriteBack).
//! Momentum and energy exact, L_z too (the gas takes all of D); where s = 1 the
//! increment as is.  Counted (CLEAN KREAL: cells, min s).

void RadiationM1::CleanStagKReal() {
  cln_krel_on = false;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  FluidRef fl = FluidRef::Get(pmy_pack);
  // (with clean_stag_kpos the bound is on kpos's E_old, CleanStagKPos)
  if (!cln_kreal || cln_kpos || cln_cell || !stage || !cln_feedback || !fl.on) return;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int is = id.is, js = id.js, ks = id.ks, nx1 = id.nx1, nx2 = id.nx2;
  const int nji = nx2*nx1, nkji = id.nx3*nji, nmb = pmy_pack->nmb_thispack;
  auto uh = fl.u0;
  auto u0_ = u0;
  auto inc = t2inc;
  auto w = cln_w;
  CleanArgs x = CleanMakeArgs();
  const Real th = CLN_KPOS_TH;
  const bool kps = cln_kpos;
  Real nc = 0.0, smin = 1.0;
  Kokkos::parallel_reduce("cln_kreal", CleanRp256(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n, Real &sm) {
    const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/nx1;
    const int i = is + idx%nx1;
    // (clean_stag_kpos's applied share CLW_KS: 1, t2inc holds s D)
    if (kps) w(m,CLW_KS,k,j,i) = 1.0;
    const Real d = uh(m,IDN,k,j,i);
    const Real us[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
    const Real ud[3] = {inc(m,M1_T2_M1,k,j,i), inc(m,M1_T2_M1+1,k,j,i),
                        inc(m,M1_T2_M1+2,k,j,i)};
    Real ps[3], pd[3], a = 0.0, b = 0.0;
    CleanGasVector(x, m, k, j, i, us, ps);
    CleanGasVector(x, m, k, j, i, ud, pd);
    for (int q = 0; q < 3; ++q) {
      a += ps[q]*pd[q]/d;
      b += 0.5*SQR(pd[q])/d;
      w(m,CLW_KR+q,k,j,i) = 0.0;
    }
    const Real w1 = a + b, ea = u0_(m,M1_E,k,j,i) + inc(m,M1_T2_E,k,j,i) + w1;
    const Real bd = (1.0 - th)*fmax(ea, 0.0);
    if (!(w1 > bd) || !(b > 0.0)) return;
    // (the largest root of b s^2 + a s - bd = 0, without cancellation)
    const Real sq = std::sqrt(a*a + 4.0*b*bd);
    Real s = (a >= 0.0) ? ((a + sq > 0.0) ? 2.0*bd/(a + sq) : 0.0) : (sq - a)/(2.0*b);
    s = fmin(fmax(s, 0.0), 1.0);
    const Real dw = w1 - (a + b*s)*s;
    for (int q = 0; q < 3; ++q) {
      inc(m,M1_T2_M1+q,k,j,i) = s*ud[q];
      w(m,CLW_KR+q,k,j,i) = (1.0 - s)*pd[q];
    }
    inc(m,M1_T2_E,k,j,i) += dw;
    inc(m,M1_T2_EN,k,j,i) -= dw;
    n += 1.0;
    sm = fmin(sm, s);
  }, Kokkos::Sum<Real>(nc), Kokkos::Min<Real>(smin));
  nc = CleanReduce(nc, false);
  smin = -CleanReduce(-smin, true);
  cln_krel_on = (nc > 0.0);
  if (cln_krel_on && global_variable::my_rank == 0) {
    std::printf("CLEAN KREAL cycle %d t2s %d cells %.0f min s %.3e\n",
                pmy_pack->pmesh->ncycle, t2_solve, nc, smin);
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

  // (counted: gas eint < 0 of the stage start, which the exchange takes as 0 and the
  // write-back's implicit_pos_gas meets, and E < 0)
  // (a flat range, as cln_newton: A100 ran the 4-D MDRange reduction at 39 ms; the
  // counts are whole numbers, exact in any order)
  // (and a stage start that is not finite: a NaN the hydro stage left, which eint = max(e0,
  // 0) would hide; the solve then fails and the step is redone, CLEAN WARNING, rule 9)
  Real nn[3] = {0.0, 0.0, 0.0};
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Kokkos::parallel_reduce("cln_old", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &ng, Real &ne, Real &nf) {
    const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
    const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
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
      const Real e0 = CleanGasEint(x, uh, ekin, m, k, j, i), eint = fmax(e0, 0.0);
      if (e0 < 0.0) ng += 1.0;
      if (!(Kokkos::isfinite(e0) && Kokkos::isfinite(d))) nf += 1.0;
      t = eos.Temperature(d, eint);
      w(m,CLW_EI,k,j,i) = eint;
      w(m,CLW_JV,k,j,i) = e0 - eint;   // (the gas's deficit below 0: CleanWriteBack)
    }
    w(m,CLW_T,k,j,i) = t;
    w(m,CLW_EO,k,j,i) = u0_(m,M1_E,k,j,i) + (stage ? inc(m,M1_T2_E,k,j,i) : 0.0);
    if (w(m,CLW_EO,k,j,i) < 0.0) ne += 1.0;
    if (!Kokkos::isfinite(w(m,CLW_EO,k,j,i))) nf += 1.0;
  }, Kokkos::Sum<Real>(nn[0]), Kokkos::Sum<Real>(nn[1]), Kokkos::Sum<Real>(nn[2]));
  nn[0] = CleanReduce(nn[0], false);
  nn[1] = CleanReduce(nn[1], false);
  nn[2] = CleanReduce(nn[2], false);
  // (clean_cs_lz_abs = 2: the increment's L_z, of the slopes' sums, Time2FormStage)
  if (cln_lzbook && fluid && feedback && stage) CleanLzBook(cln_lzc, true);
  if (global_variable::my_rank == 0 && nn[2] > 0.0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: stage start not finite in %.0f cells (gas or "
                "E; the solve fails and the step is redone)\n", pmy_pack->pmesh->ncycle,
                t2_solve, nn[2]);
  }
  if (global_variable::my_rank == 0 && (nn[0] > 0.0 || nn[1] > 0.0)) {
    std::printf("CLEAN OLDVEC cycle %d t2s %d eint<0 %.0f (taken as 0) E<0 %.0f\n",
                pmy_pack->pmesh->ncycle, t2_solve, nn[0], nn[1]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagXDonor
//! \brief clean_stag_xdonor (staggered): the off-diagonal stress flux of each face,
//! C = -(1-w) theta h c^2 d_t(D_nt E) (CleanFaceCross; cubed sphere: the lateral sums),
//! of the iterate E_k enters the pass as a donor-cell coefficient of the face,
//!   C(E) -> q E_donor,  q = C(E_k)/E_k,donor,  |q| <= c,
//! the donor the cell the flux leaves.  q adds to al (>= 0) or ar (<= 0) only, so the
//! operator stays an M-matrix (column dominant, conservative) for any D, h: every pass's
//! solution of a non-negative right-hand side is non-negative, without a mask.  At the
//! fixed point q E_donor = C(E) where |C| <= c E_donor (the central discretisation);
//! beyond, the cross flux is limited at free streaming.  The gas momentum keeps the
//! central stress (CleanMomentumChange: momentum and L_z by telescoping).
//! C of the iterate is a lagged (Jacobi) coupling of the face to the cells beside it,
//! with a gain ~ the light CFL where nothing else holds E (optically thin): there its
//! checkerboard mode alternates from pass to pass (RSG floor gas, tau ~ 1e-9: the
//! passes stall at the round-off of the solve and the floor-density gas velocity grows,
//! S3P4).  A face whose C turns back, its change d = l times the last one (l < 0),
//! takes d/(1 - l), the fixed point of that mode (half the change for l = -1); C is
//! unchanged at the fixed point.  first: the first iterate of a solve, no history.

void RadiationM1::CleanStagXDonor(const bool first, const bool lat, const bool hold) {
  if (!cln_xd || cln_cell || !pmy_pack->pmesh->multi_d) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  cln_xd_on = false;
  CleanGhostDefer(cln_w, CLW_EK, 1);   // (one exchange)
  CleanGhostFill(cln_w, CLW_EO, 1);   // (the donors' E_old of the lag split)
  if (lat) CleanStagCsLateral(CLW_EK);
  // (hold, clean_stag_xdonor_freeze: C and the right-hand side's split of it, its q and
  // explicit part, kept as the last pass set them: donor and cap no longer switch)
  cln_xdh = hold;
  if (hold) {cln_xd_on = true; return;}
  CleanArgs x = CleanMakeArgs();   // (x.xq = -1: the central cross flux)
  const int xq = cln_xqs;
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    // (faces outside the active range: none, every solve; the pass's split,
    // CleanStagXSplit in the right-hand side, sets the active ones) (one flat range over
    // every face: A100 R136 ran the 5-D subview fill and the 4-D par_for at 0.55 ms each)
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
      // (the change over the last one, l < 0 where C turns back: the step d/(1 - l), the
      // fixed point of a mode C -> l C exactly; continuous in l, l >= 0 as it is)
      const Real l = (d*dp < 0.0) ? d/dp : 0.0;
      const Real cu = first ? cn : cp + d/(1.0 - l);
      fc(m,xq+3,k,j,i) = cu;
      fc(m,xq+4,k,j,i) = first ? 0.0 : cu - cp;
    });
  }
  cln_xd_on = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagKPos
//! \brief clean_stag_kpos (staggered, hesdirk2 stage solves): the explicit increment of
//! the old vector, E_old = E_s + c1 K1_E (+ c2 K2_E), E_s the stage start, in the
//! modified-Patankar form of MPRK22 (Burchard, Deleersnijder & Meister 2003), by
//! construction positive for any h, without selection or rounds.  The increment is
//! -div G + Q, G = c1 Phi1 + c2 Phi2 the energy its slopes' solves moved across each
//! face (t2pf1, t2pf2: Phi without these terms, CleanWriteBack), Q the cell part.
//! - every interior face's G enters the stage as q X_donor, q = |G|/(h E_ref,donor)
//!   (CleanFaceWeights, the E rows and the write-back's energy only, never the momentum);
//! - a physical x1 boundary's outflow and the removal R = -Q > 0 (radiation the slope
//!   gave the gas: the eint its increment gained, which the gas gives back at the stage
//!   start; the work on its kinetic energy stays explicit, counted) enter the row's
//!   diagonal as pd X, pd = (R_b + R_g)/E_ref (CLW_KM; CLW_KF the
//!   boundary's part, which the write-back books as boundary energy; the gas row takes
//!   pd_g X; without gas feedback the removal into the fixed gas on the diagonal whole);
//! - E_old = E_s + Q+ + the boundary inflows - the gas-limited rest (counted), >= 0.
//! E_ref = max(E_s + tau K1_E, th E_s), tau = c1 + c2 + h, the predictor of the stage
//! value (X/E_ref = 1 + O(h^2): 2nd order), from stored state only.  q and pd are
//! constants of the stage: the operator stays an M-matrix (al >= 0 >= ar, column
//! dominant), b >= 0, so X >= 0; one flux per face and the gas taking pd_g X: E + gas
//! exact for any iterate.  The slopes and t2pf leave these terms out.  q, pd at most
//! 1e12 c, 1e12 (E_ref -> 0), counted (CLEAN KPOS).

void RadiationM1::CleanStagKPos() {
  cln_kq_on = false;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  if (!cln_kpos || cln_cell || !stage) return;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool two = (t2_solve == M1_T2S_STAGE2);   // (stage 2: K1 and K2)
  if (!t2pf1_ok || (two && !t2pf2_ok) || t2pf1.extent_int(0) <= nmb1) {
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN KPOS cycle %d t2s %d: no stored face slope, increment as is\n",
                  pmy_pack->pmesh->ncycle, t2_solve);
    }
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto w = cln_w;
  auto u0_ = u0;
  auto k1v = t2k1, k2v = t2k2;
  CleanArgs x = CleanMakeArgs();
  const Real c1 = cln_kc1, c2 = cln_kc2, th = CLN_KPOS_TH;
  const Real qcap = 1.0e12*x.c;
  auto p1 = t2pf1, p2 = t2pf2;
  // the references, the E each slope was evaluated at (MPRK22): K1 = f(U^n) at E_s
  // (CLW_PW), K2 = f(Y1) at the stage-1 solution (stage 2: CLW_KF, cln_ky1), with the
  // donors' ghosts
  const bool y1 = two && cln_ky1.extent_int(0) > nmb1;
  auto ky = cln_ky1;
  CleanFor("cln_kpos_ref", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_PW,k,j,i) = fmax(u0_(m,M1_E,k,j,i), 0.0);
    w(m,CLW_KF,k,j,i) = y1 ? fmax(ky(m,k,j,i), 0.0) : w(m,CLW_PW,k,j,i);
  });
  // (the donors whose faces move: all with faces, none yet with cell / redo)
  const Real kt0 = (cln_kmode == 1 || cln_kmode == 4) ? 1.0 : 0.0;
  const bool kdl = (cln_kmode == 4);   // (delta: the faces' q later, CleanStagKDelta)
  CleanFor("cln_kpos_kt", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_KT,k,j,i) = kt0;
  });
  CleanGhostFill(cln_w, CLW_PW, 1);
  if (two) CleanGhostFill(cln_w, CLW_KF, 1);
  CleanGhostFill(cln_w, CLW_KT, 1);
  cln_kredo = 0;
  cln_kredo_neg = 0.0;
  // 1. the faces' q
  const int kq = cln_kqs;
  Real nq[2] = {0.0, 0.0};   // (faces moved, at the guard)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    Real nd[2] = {0.0, 0.0};
    Kokkos::parallel_reduce("cln_kpos_q", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
        DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + dk + 1, je + dj + 1, ie + di + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &n, Real &nc) {
      Real ql = 0.0, qr = 0.0;
      // (each slope's G over its own reference of the donor)
      auto put = [&](const Real g, const int sl) {
        if (g == 0.0) return;
        const bool up = (g > 0.0);
        const Real ed = up ? w(m,sl,k-dk,j-dj,i-di) : w(m,sl,k,j,i);
        if (!((up ? w(m,CLW_KT,k-dk,j-dj,i-di) : w(m,CLW_KT,k,j,i)) > 0.0)) return;
        Real q = qcap;
        if (fabs(g) < qcap*x.h*ed) {
          q = fabs(g)/(x.h*ed);
        } else {
          nc += 1.0;
        }
        if (up) ql += q;
        if (!up) qr -= q;
      };
      if (!kdl && !CleanBoundaryFace(x, dir, m, k, j, i)) {
        put(c1*p1(m,dir,k,j,i), CLW_PW);
        if (two) put(c2*p2(m,dir,k,j,i), CLW_KF);
        if (ql != 0.0 || qr != 0.0) n += 1.0;
      }
      fc(m,kq,k,j,i) = ql;
      fc(m,kq+1,k,j,i) = qr;
    }, Kokkos::Sum<Real>(nd[0]), Kokkos::Sum<Real>(nd[1]));
    nq[0] += nd[0];
    nq[1] += nd[1];
  }
  // 2. the old vector, pd, the gas's give-back
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto inc = t2inc;
  const bool gas = cln_feedback, krk = cln_kreal && gas;
  // (cell / redo: a boundary outflow explicit unless it takes E_old below 0)
  const bool kbt = (cln_kmode == 2 || cln_kmode == 3);
  // (cells with pd, gas-limited, pd guard, E_old < 0, of which the stage start's eint < 0)
  Real nn[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  Real nkr = 0.0, skr = 1.0;   // (clean_stag_kreal: cells limited, min s)
  Kokkos::parallel_reduce("cln_kpos", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &np, Real &ng,
                Real &nc, Real &ne, Real &nj, Real &lk, Real &ls) {
    const Real es = u0_(m,M1_E,k,j,i);
    if (krk) for (int n = 0; n < 3; ++n) w(m,CLW_KR+n,k,j,i) = 0.0;
    Real e = es, dg = 0.0, dgs[2] = {0.0, 0.0}, rbs[2] = {0.0, 0.0};
    for (int d = 0; d < x.nf; ++d) {
      const int dk = (d == 2), dj = (d == 1), di = (d == 0);
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*dk, jf = j + sd*dj, iff = i + sd*di;
        const bool bf = CleanBoundaryFace(x, d, m, kf, jf, iff);
        for (int s = 0; s < 1 + two; ++s) {
          const Real g = s ? c2*p2(m,d,kf,jf,iff) : c1*p1(m,d,kf,jf,iff);
          const Real a = CleanDivW(x, d, m, k, j, i, sd)*(sd ? -g : g);   // (into i)
          dg += a;
          dgs[s] += a;
          if (!bf) {   // (a face whose donor's faces do not move: explicit)
            const int sn = 2*sd - 1;
            const Real fd = (a < 0.0) ? w(m,CLW_KT,k,j,i) :
                            w(m,CLW_KT,k+sn*dk,j+sn*dj,i+sn*di);
            if (!(fd > 0.0)) e += a;
            continue;
          }
          if (a >= 0.0) {
            e += a;        // (a boundary's inflow explicit)
          } else if (d == 0) {
            rbs[s] -= a;   // (an x1 boundary's outflow on the diagonal)
          } else {
            e += a;        // (another physical boundary's outflow explicit, counted)
          }
        }
      }
    }
    const Real q = w(m,CLW_EO,k,j,i) - es - dg;
    Real rg = 0.0, wke = 0.0;   // (wke: the kinetic-energy work removed explicitly)
    if (q >= 0.0) {
      e += q;
    } else if (gas && w(m,CLW_JV,k,j,i) >= 0.0) {   // (a stage start eint >= 0)
      // (the gas gives back the eint its increment gained, de = inc_EN less the kinetic
      // energy of the momentum increment, at most -Q: eint stays >= that of the stage
      // start; the rest, the work on the gas's kinetic energy, a removal of E, counted)
      Real m2 = 0.0, m2s = 0.0;
      for (int n = 0; n < 3; ++n) {
        m2 += SQR(uh(m,IM1+n,k,j,i));
        m2s += SQR(uh(m,IM1+n,k,j,i) - inc(m,M1_T2_M1+n,k,j,i));
      }
      const Real de = inc(m,M1_T2_EN,k,j,i) - 0.5*(m2 - m2s)/uh(m,IDN,k,j,i);
      rg = fmin(fmax(de, 0.0), fmin(-q, fmax(w(m,CLW_EI,k,j,i), 0.0)));
      Real rem = -q - rg;   // (the explicit rest: the kinetic-energy work W(1))
      // (clean_stag_kreal: the increment's kick D at the largest share s whose work this
      // E_old, KW = e, can pay with the rest: W(s) <= W(1) - (rem - (1 - th) KW), W(s) =
      // a s + b s^2 (CleanStagKReal's form); the gas takes m_s + s D and its energy less
      // W(1) - W(s) (eint as formed), the stage the rest (1 - s) D (CLW_KR))
      if (krk && rem > 0.0) {
        const Real d = uh(m,IDN,k,j,i);
        const Real ud[3] = {inc(m,M1_T2_M1,k,j,i), inc(m,M1_T2_M1+1,k,j,i),
                            inc(m,M1_T2_M1+2,k,j,i)};
        const Real us[3] = {uh(m,IM1,k,j,i) - ud[0], uh(m,IM2,k,j,i) - ud[1],
                            uh(m,IM3,k,j,i) - ud[2]};
        Real ps[3], pv[3], a = 0.0, b = 0.0;
        CleanGasVector(x, m, k, j, i, us, ps);
        CleanGasVector(x, m, k, j, i, ud, pv);
        for (int n = 0; n < 3; ++n) {
          a += ps[n]*pv[n]/d;
          b += 0.5*SQR(pv[n])/d;
        }
        const Real w1 = a + b, exc = rem - (1.0 - th)*fmax(e, 0.0);
        if (exc > 0.0 && b > 0.0 && w1 > 0.0) {
          const Real bd = fmax(w1 - exc, 0.0);
          const Real sq = std::sqrt(a*a + 4.0*b*bd);
          Real s = (a >= 0.0) ? ((a + sq > 0.0) ? 2.0*bd/(a + sq) : 0.0) :
                   (sq - a)/(2.0*b);
          s = fmin(fmax(s, 0.0), 1.0);
          const Real dw = w1 - (a + b*s)*s;
          for (int n = 0; n < 3; ++n) {
            uh(m,IM1+n,k,j,i) -= (1.0 - s)*ud[n];
            w(m,CLW_KR+n,k,j,i) = (1.0 - s)*pv[n];
          }
          uh(m,IEN,k,j,i) -= dw;
          rem -= dw;
          lk += 1.0;
          ls = fmin(ls, s);
        }
      }
      if (rem > 0.0) {
        e -= rem;
        wke = rem;
        ng += 1.0;
      }
    } else if (!gas) {
      rg = -q;   // (no gas feedback: the removal into a fixed gas, on the diagonal whole)
    } else {
      e += q;
      ng += 1.0;
    }
    if (gas && rg > 0.0) {
      uh(m,IEN,k,j,i) -= rg;
      w(m,CLW_EI,k,j,i) -= rg;
      w(m,CLW_T,k,j,i) = eos.Temperature(uh(m,IDN,k,j,i), w(m,CLW_EI,k,j,i));
    }
    // (pd: each slope's part over its reference; the cell part split by the slopes'
    // own cell removals, -(c K_E - dg) of each)
    Real qs[2] = {-(c1*k1v(m,M1_T2_E,k,j,i) - dgs[0]),
                  two ? -(c2*k2v(m,M1_T2_E,k,j,i) - dgs[1]) : 0.0};
    qs[0] = fmax(qs[0], 0.0);
    qs[1] = fmax(qs[1], 0.0);
    const Real qt = qs[0] + qs[1];
    // (an x1 boundary's outflow stays explicit where E_old can pay it, rb <= (1 - th)
    // E_old: the well-balanced decks' photosphere keeps its explicit flux; only past that
    // on the diagonal)
    if (rbs[0] + rbs[1] <= (kbt ? 1.0 : 1.0 - th)*fmax(e, 0.0)) {
      e -= rbs[0] + rbs[1];
      rbs[0] = rbs[1] = 0.0;
    }
    Real pd = 0.0, pb = 0.0;
    for (int s = 0; s < 1 + two; ++s) {
      const Real er = w(m,s ? CLW_KF : CLW_PW,k,j,i);
      const Real rs = rbs[s] + ((qt > 0.0) ? rg*qs[s]/qt : (s ? 0.0 : rg));
      if (!(rs > 0.0)) continue;
      if (rs < 1.0e12*er) {
        pd += rs/er;
        pb += rbs[s]/er;
      } else {
        pd += 1.0e12;
        pb += 1.0e12*rbs[s]/rs;
        nc += 1.0;
      }
    }
    if (pd > 0.0) np += 1.0;
    if (e < 0.0) ne += 1.0;
    if (e < 0.0 && gas && w(m,CLW_JV,k,j,i) < 0.0) nj += 1.0;
    w(m,CLW_KF,k,j,i) = pb;
    w(m,CLW_KM,k,j,i) = pd;
    w(m,CLW_EO,k,j,i) = e;
  }, Kokkos::Sum<Real>(nn[0]), Kokkos::Sum<Real>(nn[1]), Kokkos::Sum<Real>(nn[2]),
     Kokkos::Sum<Real>(nn[3]), Kokkos::Sum<Real>(nn[4]), Kokkos::Sum<Real>(nkr),
     Kokkos::Min<Real>(skr));
  for (int n = 0; n < 5; ++n) nn[n] = CleanReduce(nn[n], false);
  // (the history's boundary energy: a moved x1 boundary outflow (CLW_KF > 0) leaves the
  // old vector; the write-back books pd_b X)
  if (cln_lacc.extent_int(0) > nmb1) {
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;   // (the x1 face area of the clean fluxes)
    par_for("cln_kpos_bnd", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      if (!CleanBoundaryFace(x, 0, m, k, j, i)) return;
      const Real g = ar(m,k,j,i)*(c1*p1(m,0,k,j,i) + (two ? c2*p2(m,0,k,j,i) : 0.0));
      if (s ? !(g > 0.0) : !(g < 0.0)) return;   // (an outflow: its donor is the cell)
      if (!(w(m,CLW_KF,k,j,s ? ie : is) > 0.0)) return;   // (kept explicit)
      la(m,2+s,k,j) -= g;
    });
  }
  nkr = CleanReduce(nkr, false);
  skr = -CleanReduce(-skr, true);
  cln_krel_on = (nkr > 0.0);
  if (cln_krel_on && global_variable::my_rank == 0) {
    std::printf("CLEAN KREAL cycle %d t2s %d cells %.0f min s %.3e (on the kpos E_old)\n",
                pmy_pack->pmesh->ncycle, t2_solve, nkr, skr);
  }
  nq[0] = CleanReduce(nq[0], false);
  nq[1] = CleanReduce(nq[1], false);
  cln_kq_on = true;
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN KPOS cycle %d t2s %d faces %.0f cells pd %.0f guard %.0f/%.0f "
                "gas-limited %.0f Eo<0 %.0f (eint<0 at the stage start %.0f)\n",
                pmy_pack->pmesh->ncycle, t2_solve, nq[0], nn[0], nq[1], nn[2], nn[1], nn[3],
                nn[4]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagKDelta
//! \brief clean_stag_kpos = delta (after CleanFaceCoefficients): each interior face's slope
//! flux G = c1 Phi1 + c2 Phi2 splits into tau (al E_L + ar E_R), tau = c1 + c2, with the
//! stage's own coefficients (al >= 0 >= ar, CleanFaceWeights), taken implicit in X
//! (x.kt = tau/h scales al, ar of the E rows), and the rest dG = G - tau (al E_s,L + ar
//! E_s,R), moved as q X_donor (q = |dG|/(h E_s,donor), the donor by dG's sign).  The
//! old vector has no face part (as faces), so b >= 0 and the operator an M-matrix: X >= 0
//! by construction; at X = E_s the face carries G exactly, and near equilibrium dG -> 0
//! (the moved part vanishes; the well-balanced and steady states as the stage's own
//! flux).  Counted (CLEAN KPOS delta).

void RadiationM1::CleanStagKDelta() {
  if (!cln_kq_on || cln_kmode != 4) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool two = (t2_solve == M1_T2S_STAGE2);
  const Real c1 = cln_kc1, c2 = cln_kc2, tau = c1 + (two ? c2 : 0.0);
  cln_ktau = tau/dt_sub;
  CleanArgs x = CleanMakeArgs(), y = x;
  y.kq = -1;
  y.xq = -1;
  y.kt = 0.0;
  auto w = cln_w;
  auto p1 = t2pf1, p2 = t2pf2;
  const int kq = cln_kqs;
  const Real qcap = 1.0e12*x.c;
  Real nq[2] = {0.0, 0.0};   // (faces moved, at the guard)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    Real nd[2] = {0.0, 0.0};
    Kokkos::parallel_reduce("cln_kdelta", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
        DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + dk + 1, je + dj + 1, ie + di + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &n, Real &nc) {
      Real ql = 0.0, qr = 0.0;
      if (!CleanBoundaryFace(y, dir, m, k, j, i)) {
        Real al, ar;
        CleanFaceWeights(y, dir, m, k, j, i, al, ar);
        const Real el = w(m,CLW_PW,k-dk,j-dj,i-di), er = w(m,CLW_PW,k,j,i);
        const Real g = c1*p1(m,dir,k,j,i) + (two ? c2*p2(m,dir,k,j,i) : 0.0);
        const Real dg = g - tau*(al*el + ar*er);
        if (dg != 0.0) {
          const bool up = (dg > 0.0);
          const Real ed = up ? el : er;
          Real q = qcap;
          if (fabs(dg) < qcap*y.h*ed) {
            q = fabs(dg)/(y.h*ed);
          } else {
            nc += 1.0;
          }
          if (up) ql = q;
          if (!up) qr = -q;
          n += 1.0;
        }
      }
      fc(m,kq,k,j,i) = ql;
      fc(m,kq+1,k,j,i) = qr;
    }, Kokkos::Sum<Real>(nd[0]), Kokkos::Sum<Real>(nd[1]));
    nq[0] += nd[0];
    nq[1] += nd[1];
  }
  nq[0] = CleanReduce(nq[0], false);
  nq[1] = CleanReduce(nq[1], false);
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN KPOS delta cycle %d t2s %d tau/h %.3e faces dG %.0f guard %.0f\n",
                pmy_pack->pmesh->ncycle, t2_solve, cln_ktau, nq[0], nq[1]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::CleanStagKRedo
//! \brief clean_stag_kpos = redo: after the passes, the cells the solve left at X <= 0
//! (globally counted) become donors whose slope face fluxes move (q X_donor, as faces);
//! their explicit part leaves the old vector on both sides of each face, and the stage
//! is solved again (true).  At most CLN_KPOS_REDO rounds; the last moves every donor's
//! faces (= faces: b >= 0, an M-matrix, X >= 0 by construction).  Deterministic (flags
//! from stored state, no in-place neighbour reads), counted (CLEAN KPOS redo).

bool RadiationM1::CleanStagKRedo() {
  if (cln_kmode != 3 || !cln_kq_on || cln_kredo >= CLN_KPOS_REDO) return false;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  auto u0_ = u0;
  // (the last round, all faces: the cap, or a round that did not lower the count of
  // X <= 0 cells (a receiver cascade; cheaper than more rounds))
  const bool cap = (cln_kredo + 1 >= CLN_KPOS_REDO);
  // 1. the cells X <= 0 (counted), then the new donors (KT = 2): X <= 0 there or at a face
  //    neighbour (one layer: the receivers that would lose a moved inflow move with it),
  //    not moved yet; the last round all
  CleanGhostFill(cln_w, CLW_X, 1);
  CleanArgs xa = CleanMakeArgs();
  const int nf = xa.nf;
  Real nn[2] = {0.0, 0.0};   // (cells X <= 0, new donors)
  Kokkos::parallel_reduce("cln_kredo_c", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a) {
    if (!(w(m,CLW_X,k,j,i) > 0.0)) a += 1.0;
  }, Kokkos::Sum<Real>(nn[0]));
  nn[0] = CleanReduce(nn[0], false);
  if (!(nn[0] > 0.0)) return false;   // (positive: no round)
  const bool last = cap || (cln_kredo > 0 && nn[0] >= cln_kredo_neg);
  cln_kredo_neg = nn[0];
  Kokkos::parallel_reduce("cln_kredo_n", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &b) {
    const bool neg = !(w(m,CLW_X,k,j,i) > 0.0);
    bool nb = neg;
    for (int d = 0; d < nf && !nb; ++d) {
      const int dk = (d == 2), dj = (d == 1), di = (d == 0);
      nb = (!CleanBoundaryFace(xa, d, m, k, j, i) && !(w(m,CLW_X,k-dk,j-dj,i-di) > 0.0)) ||
           (!CleanBoundaryFace(xa, d, m, k+dk, j+dj, i+di) &&
            !(w(m,CLW_X,k+dk,j+dj,i+di) > 0.0));
    }
    if ((nb || last) && !(w(m,CLW_KT,k,j,i) > 0.0)) {
      w(m,CLW_KT,k,j,i) = 2.0;
      b += 1.0;
    }
  }, Kokkos::Sum<Real>(nn[1]));
  nn[1] = CleanReduce(nn[1], false);
  ++cln_kredo;
  // 2. the references (K1: E_s in CLW_PW; K2: the stage-1 solution in CLW_KS)
  const bool two = (t2_solve == M1_T2S_STAGE2);
  const bool y1 = two && cln_ky1.extent_int(0) > nmb1;
  auto ky = cln_ky1;
  CleanFor("cln_kredo_ref", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_PW,k,j,i) = fmax(u0_(m,M1_E,k,j,i), 0.0);
    w(m,CLW_KS,k,j,i) = y1 ? fmax(ky(m,k,j,i), 0.0) : w(m,CLW_PW,k,j,i);
  });
  CleanGhostFill(cln_w, CLW_PW, 1);
  CleanGhostFill(cln_w, CLW_KS, 1);
  CleanGhostFill(cln_w, CLW_KT, 1);
  CleanArgs x = CleanMakeArgs();
  const Real c1 = cln_kc1, c2 = cln_kc2, qcap = 1.0e12*x.c;
  auto p1 = t2pf1, p2 = t2pf2;
  const int kq = cln_kqs;
  // 3. the new donors' faces: q added, their explicit part out of E_old
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    par_for("cln_kredo_q", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      auto put = [&](const Real g, const int sl) {
        if (g == 0.0) return;
        const bool up = (g > 0.0);
        if (!((up ? w(m,CLW_KT,k-dk,j-dj,i-di) : w(m,CLW_KT,k,j,i)) > 1.5)) return;
        const Real ed = up ? w(m,sl,k-dk,j-dj,i-di) : w(m,sl,k,j,i);
        const Real q = (fabs(g) < qcap*x.h*ed) ? fabs(g)/(x.h*ed) : qcap;
        if (up) fc(m,kq,k,j,i) += q;
        if (!up) fc(m,kq+1,k,j,i) -= q;
      };
      put(c1*p1(m,dir,k,j,i), CLW_PW);
      if (two) put(c2*p2(m,dir,k,j,i), CLW_KS);
    });
  }
  CleanFor("cln_kredo_eo", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real de = 0.0;
    for (int d = 0; d < x.nf; ++d) {
      const int dk = (d == 2), dj = (d == 1), di = (d == 0);
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*dk, jf = j + sd*dj, iff = i + sd*di, sn = 2*sd - 1;
        if (CleanBoundaryFace(x, d, m, kf, jf, iff)) continue;
        for (int s = 0; s < 1 + two; ++s) {
          const Real g = s ? c2*p2(m,d,kf,jf,iff) : c1*p1(m,d,kf,jf,iff);
          const Real a = CleanDivW(x, d, m, k, j, i, sd)*(sd ? -g : g);   // (into i)
          const Real fd = (a < 0.0) ? w(m,CLW_KT,k,j,i) :
                          w(m,CLW_KT,k+sn*dk,j+sn*dj,i+sn*di);
          if (fd > 1.5) de -= a;
        }
      }
    }
    w(m,CLW_EO,k,j,i) += de;
  });
  CleanFor("cln_kredo_1", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (w(m,CLW_KT,k,j,i) > 1.5) w(m,CLW_KT,k,j,i) = 1.0;
  });
  CleanGhostFill(cln_w, CLW_EO, 1);
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN KPOS redo cycle %d t2s %d round %d X<=0 %.0f new donors %.0f%s\n",
                pmy_pack->pmesh->ncycle, t2_solve, cln_kredo, nn[0], nn[1],
                last ? " (last: all)" : "");
  }
  return true;
}
#endif  // M1_TU_PART(0)

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagKEnd
//! \brief clean_stag_xpat after the passes: the counts of the solve, exchange rows in
//! Patankar form (last pass; sum over passes)

#if M1_TU_PART(1)
void RadiationM1::CleanStagKEnd() {
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_kcnt);
  Real v[2] = {h(0), h(1)};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank == 0 && v[1] > 0.0) {
    std::printf("CLEAN XPAT cycle %d t2s %d end: rows %.0f (sum %.0f)\n",
                pmy_pack->pmesh->ncycle, t2_solve, v[0], v[1]);
  }
  Kokkos::deep_copy(cln_kcnt, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagXSweep
//! \brief clean_stag_xsweep: one line-Jacobi sweep from the linear solve's X,
//!   X = X+ + L^{-1} (B - A X+),  X+ = max(X, 0),
//! L the x1-line tridiagonal part of A with A's diagonal (CleanPreconditionerSetup).  A
//! an M-matrix (clean_stag_xdonor), A = L - N is a regular splitting (L^{-1} >= 0,
//! N >= 0), so X = L^{-1} (B + N X+) >= 0 where B >= 0, and the sweep contracts the
//! error of the solution.

void RadiationM1::CleanStagXSweep() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  par_for("cln_xsw0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_V,k,j,i) = fmax(w(m,CLW_X,k,j,i), 0.0);
  });
  CleanApply(CLW_V, CLW_S);
  par_for("cln_xsw1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_R,k,j,i) = w(m,CLW_B,k,j,i) - w(m,CLW_S,k,j,i);
  });
  CleanStagLineSolve(CleanStagLevel(0), CLW_R, CLW_Z, -1);
  par_for("cln_xsw2", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_X,k,j,i) = w(m,CLW_V,k,j,i) + w(m,CLW_Z,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanWorstCell
//! \brief CLEAN WORST (diagnostic) of a backward-Euler solve that did not converge: its
//! kind from the last two changes (STALL: fell less than 2x, a cut dt does not help) and
//! the cell of the largest |dE|/S of the last pass, printed by its rank

void RadiationM1::CleanWorstCell(const Real dnl, const Real dnl0) {
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
      cv(4) = op(m,M1_OP_T,k,j,i);
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
//! \fn void RadiationM1::CleanStagXCount
//! \brief clean_stag_xdonor: the pass's cells with X <= 0 and faces at the cap (CLEAN
//! XDONOR, diagnostic)

void RadiationM1::CleanStagXCount(const int inl, const int nit) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  auto w = cln_w;
  Real nn = 0.0, xmin = 0.0, nb = 0.0, nbx = 0.0, neo = 0.0, nsd = 0.0;
  CleanArgs x = CleanMakeArgs();
  // (of the X <= 0 cells: the row term most negative relative to E_old)
  Real tw = 0.0, te = 0.0, ts = 0.0, tg = 0.0;
  // (flat range reductions, as cln_old: A100 R136 ran the 4-D MDRange ones at 6.7 ms
  // (cln_xcount) and 1.6 ms per launch (cln_xfaces, per direction); every count is a
  // whole number and the minimum is exact, in any order)
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  {
    auto xdg = cln_xdg;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("cln_xcount", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &a, Real &b, Real &c, Real &d, Real &e, Real &f,
                  Real &tw_, Real &te_, Real &ts_, Real &tg_) {
      const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
      const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
      const Real v = w(m,CLW_X,k,j,i);
      if (!(v > 0.0)) {
        int n = 0;
        for (int q = 1; q < 4; ++q) if (xdg(m,q,k,j,i) < xdg(m,n,k,j,i)) n = q;
        if (xdg(m,n,k,j,i) < 0.0) {
          if (n == 0) tw_ += 1.0;
          if (n == 1) te_ += 1.0;
          if (n == 2) ts_ += 1.0;
          if (n == 3) tg_ += 1.0;
        }
      }
      if (!(v > 0.0)) a += 1.0;
      if (w(m,CLW_B,k,j,i) < 0.0) c += 1.0;
      if (!(v > 0.0) && w(m,CLW_B,k,j,i) < 0.0) d += 1.0;
      // (of the X <= 0 cells: the old vector < 0, or its face sources' outflow above it)
      const Real eo = w(m,CLW_EO,k,j,i);
      if (!(v > 0.0) && eo < 0.0) e += 1.0;
      if (!(v > 0.0) && eo >= 0.0) {   // (CLF_SR holds ar now: the sources anew)
        Real sd = 0.0;
        for (int q = 0; q < x.nf; ++q) {
          const int dk = (q == 2), dj = (q == 1), di = (q == 0);
          sd += CleanDivW(x, q, m, k, j, i, 1)*
                    CleanFaceSource(x, q, m, k+dk, j+dj, i+di) -
                CleanDivW(x, q, m, k, j, i, 0)*CleanFaceSource(x, q, m, k, j, i);
        }
        if (eo < x.h*sd) f += 1.0;
      }
      b = fmin(b, v/fmax(fabs(eo), 1.0e-300));
    }, Kokkos::Sum<Real>(nn), Kokkos::Min<Real>(xmin), Kokkos::Sum<Real>(nb),
       Kokkos::Sum<Real>(nbx), Kokkos::Sum<Real>(neo), Kokkos::Sum<Real>(nsd),
       Kokkos::Sum<Real>(tw), Kokkos::Sum<Real>(te), Kokkos::Sum<Real>(ts),
       Kokkos::Sum<Real>(tg));
  }
  // (faces: with a donor part, at the cap c; the directions' faces in one range)
  Real nd = 0.0, ncp = 0.0;
  if (cln_xd_on) {
    const int nf = x.nf, xq = cln_xqs;
    const Real c = x.c;
    auto fc0 = x.fc[0], fc1 = x.fc[1], fc2 = x.fc[2];
    const int ni = ie - is + 2, nji = (je - js + 2)*ni, nkji = (ke - ks + 2)*nji;
    Kokkos::parallel_reduce("cln_xfaces", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &a, Real &b) {
      const int m = idx/nkji, k = ks + (idx - m*nkji)/nji;
      const int j = js + (idx - m*nkji - (k - ks)*nji)/ni, i = is + idx%ni;
      for (int dir = 0; dir < nf; ++dir) {
        if (i > ie + (dir == 0) || j > je + (dir == 1) || k > ke + (dir == 2)) continue;
        const auto &fc = (dir == 0) ? fc0 : ((dir == 1) ? fc1 : fc2);
        const Real ql = fc(m,xq,k,j,i), qr = fc(m,xq+1,k,j,i);
        if (ql != 0.0 || qr != 0.0) a += 1.0;
        if (ql >= c || -qr >= c) b += 1.0;
      }
    }, Kokkos::Sum<Real>(nd), Kokkos::Sum<Real>(ncp));
  }
  // (the counts in one reduction over the ranks: whole numbers, exact in any order)
  Real cn[11] = {tw, te, ts, tg, nn, nb, nbx, neo, nsd, ncp, nd};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, cn, 11, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  tw = cn[0]; te = cn[1]; ts = cn[2]; tg = cn[3]; nn = cn[4]; nb = cn[5]; nbx = cn[6];
  neo = cn[7]; nsd = cn[8]; cln_xd_ncap = cn[9]; nd = cn[10];
  xmin = -CleanReduce(-xmin, true);
  if (global_variable::my_rank == 0 && (nn > 0.0 || nd > 0.0)) {
    std::printf("CLEAN XDONOR cycle %d t2s %d pass %d xle0 %.0f minX/Eo %.3e donor %.0f "
                "cap %.0f b<0 %.0f (with X<=0 %.0f; Eo<0 %.0f, Eo<hdivS %.0f) most-neg "
                "work %.0f exch %.0f sb %.0f gres %.0f it %d\n",
                pmy_pack->pmesh->ncycle, t2_solve, inl, nn, xmin, nd, cln_xd_ncap, nb, nbx,
                neo, nsd, tw, te, ts, tg, nit);
  }
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
  const bool xtn = cln_xtnorm;   // (clean_xt_norm: E responds, 1/(1 + hc kE) folded)
  const bool xtg = cln_xtcons;   // (clean_xt_cons: G = 1 - beta^2/2 of the gas)
  const Real cl = c_light;
  const bool kpg = cln_kq_on;   // (clean_stag_kpos: the gas row's pd_g E)

  CleanFor("cln_iter0", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real eo = fmax(w(m,CLW_EO,k,j,i), 0.0);
    const Real ep = pred ? eo + h*k1(m,M1_T2_E,k,j,i) : eo;
    const Real e = (ep > eo/ecap && ep < eo*ecap) ? ep : eo;
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) w(m,CLW_DM+n,k,j,i) = 0.0;
    if (feedback) {
      w(m,CLW_EP,k,j,i) = 0.0;   // (the beta terms of this T, CleanRightHandSide)
      const Real fe = xtn ? 1.0/(1.0 + hc*op(m,M1_OP_E,k,j,i)) : 1.0;
      const Real gm = xtg ? 1.0 - 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                                       SQR(uh(m,IM3,k,j,i)))/SQR(uh(m,IDN,k,j,i)*cl) : 1.0;
      const Real pg = kpg ? w(m,CLW_KM,k,j,i) - w(m,CLW_KF,k,j,i) : 0.0;
      w(m,CLW_T,k,j,i) = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
                             hc*op(m,M1_OP_P,k,j,i)*arad_*fe*gm,
                             w(m,CLW_EI,k,j,i) + (hc*op(m,M1_OP_E,k,j,i) + pg)*e*fe,
                             w(m,CLW_T,k,j,i), rtol);
    }
  });
  // (staggered: the ghosts of E_k for the gas terms, with the next exchange)
  if (feedback && !cln_cell) CleanGhostDefer(cln_w, CLW_EK, 1);
}
#endif  // M1_TU_PART(1)

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
  if (!x.upw || CleanBoundaryFace(x, n, m, kk, jj, ii)) return ac;
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

#if M1_TU_PART(1)
void RadiationM1::CleanRightHandSide(const bool beta, const bool lat) {
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
    CLEAN_HIPX(const CleanArgs *pxb = CleanDevArgs(x, 0);)
    CleanFor<0, 3>("cln_beta", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *pxb;)   // (in device memory: no per-thread copy)
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
    // (E_k's ghosts: those of its last change; BT, RE and GX in one exchange, before the
    // lateral sums, which read what an earlier fill deferred)
    CleanGhostDefer(cln_w, CLW_BT, 4);
    CleanGhostFill(cln_w, CLW_GX, 1);
    if (lat) CleanStagCsLateral(CLW_EK);
  }
  if (beta) return;
  CleanArgs xh = x;   // (the face momentum state, clean_stag_split)
  xh.hi = cln_split;
  xh.kq = -1;
  // (a pass without the F solve's F of the iterate: its F^H, each face once in cln_rhsaf
  // into CLF_N, which the next F solve sets anew)
  const bool fh = x.feedback && !cln_fsk;
  xh.fs = true;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const auto *peos = CleanDev(eos, 2);
  CleanArgs xr = x;   // (clean_stag_xdonor: the face sources as they are)
  xr.xq = -1;
  const CleanArgs *pxr = CleanDevArgs(xr, 3);
  // each face's source once, for its two cells (CLF_SR), and with gflux its a E (CLF_PT)
  const int nf = x.nf;
  const bool af = x.feedback && x.gf && !x.cub, xsp = (x.xq >= 0) && !cln_xdh;
  CleanFor<4>("cln_rhsaf", nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    for (int n = 0; n < nf; ++n) {
      if (i > ie+(n == 0) || j > je+(n == 1) || k > ke+(n == 2)) continue;
      if (xsp) CleanStagXSplit(x, *pxr, n, m, k, j, i);
      x.fc[n](m,CLF_SR,k,j,i) = CleanFaceSource(x, n, m, k, j, i);
      if (fh) x.fc[n](m,CLF_N,k,j,i) = CleanFaceFlux(xh, CLW_EK, n, m, k, j, i, true);
      if (af) x.fc[n](m,CLF_PT,k,j,i) = CleanStagAdvFace(x, xh, n, m, k, j, i);
    }
  });
  // (clean_stag_xdonor diagnostic: the row's terms, CleanStagXCount; clean_nl_print only)
  const bool xdd = cln_xd && cln_nlprint;
  if (xdd && cln_xdg.extent_int(0) < nmb1 + 1) {
    Kokkos::realloc(cln_xdg, nmb1 + 1, 4, cln_w.s.extent_int(2), cln_w.s.extent_int(3),
                    cln_w.s.extent_int(4));
  }
  auto xdg = cln_xdg;
  // (clean_stag_kpos / xpat: the rows' Patankar diagonal CLW_PD; counts of this pass)
  const bool pdo = x.pd, xpd = cln_xpat;
  const bool xte = (cln_xte > 0) && x.feedback;   // (clean_xt_exact: CLW_XA, CLW_XW)
  const bool krl = cln_krel_on && x.feedback;   // (clean_stag_kreal: CLW_KR)
  auto kc = cln_kcnt;
  if (pdo) Kokkos::deep_copy(Kokkos::subview(kc, 0), 0.0);
  CleanParForLb<3>("cln_rhs", nmb1, ks, ke, js, je, is, ie,
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
        return xh.fc[d](m,CLF_N,kk,jj,ii);
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
      if (x.cub) {   // (the cell's flux vector, CleanStagCellVector)
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
      // the exchange's slopes: the tangent at the iterate's T, Newton's quadratic rate,
      // where it takes the exchange into B as EO + (1 - r) hc kP a T^4 - r (e(T) - eint)
      // >= 0; where that is < 0, the secant through T_o, the gas T of the old vector,
      // and T (emission and gas energy): exact at T as the tangent, so the solution is
      // the same, and the exchange in B is EO + (1 - r) hc kP a T_o^4 >= 0 by
      // construction.  The tangent there (e(T) - eint >> E_old, r ~ 1: radiation
      // heating cold gas under the photosphere) gave B << 0 and E < 0, a Newton step
      // past the root of the concave E -> e map (the old production's SRCR ~ q (e^n -
      // e(T_k)) << 0).  The secant everywhere above T_o converged only linearly where
      // the chord is far from the tangent (R136: +25 % passes from thin outer cells)
      Real s = 4.0*hckp*arad_*t*t*t, cg = d*cv;
      const Real ei = x.w(m,CLW_EI,k,j,i), eo = x.w(m,CLW_EO,k,j,i);
      const bool sec = eg > ei && eo*(cg + s) + cg*hckp*arad_*SQR(SQR(t)) < s*(eg - ei);
      const Real to = sec ? eos.Temperature(d, ei, t) : t;
      if (sec && t - to > 1.0e-6*t) {
        cg = (eg - ei)/(t - to);
        s = hckp*arad_*(t + to)*(t*t + to*to);
      }
      rr = (s > 0.0) ? s/(cg + s) : 0.0;   // (no exchange: none, also where T = 0)
      // (1 - r = cg/(cg + s) without cancellation, as implicit_src_stable)
      x.w(m,CLW_RY,k,j,i) = (s > 0.0) ? cg/(cg + s) : 1.0;
      // the gas row's residual at the iterate with these sb (its T is the root with those
      // of the last pass), share r, as the old production's SRCR: a change of sb between
      // passes enters E with 1 - r, the exchange's share (radiation-dominated opaque
      // gas, r ~ 1, hcke beta F'/c >> eint: taken whole it contracted only ~0.5 per pass)
      // (clean_stag_kpos: the gas row takes pd_g E, CleanStagKPos)
      const Real pdg = (x.kq >= 0) ? x.w(m,CLW_KM,k,j,i) - x.w(m,CLW_KF,k,j,i) : 0.0;
      if (gr) gres = rr*(x.w(m,CLW_EI,k,j,i) + (hcke + pdg)*ek - sb - eg -
                         hckp*arad_*SQR(SQR(t)));
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
        // (clean_stag_kreal: with the stage's deferred kick, CleanStagKReal)
        const Real dm = x.w(m,CLW_DM+n,k,j,i) + (krl ? x.w(m,CLW_KR+n,k,j,i) : 0.0);
        work += (pg[n] + 0.5*dm)*dm/d;
        if (wsp && !wpa) work += x.h*rf[n]*(pg[n] + dm)/d;
      }
      if (wpa) work += fra*wc(m,k,j,i);
    }
    x.w(m,CLW_RR,k,j,i) = rr;
    // (eqr: the exchange through the gas row of the root T, hckp a T^4 - r hcke E =
    // EI - e(T) - sb + (1 - r) hcke E, without the cancellation of its two terms, each
    // ~ hcke E, whose round-off at h -> infinity exceeds E in opaque cells)
    if (xdd) {
      xdg(m,0,k,j,i) = -work;
      xdg(m,1,k,j,i) = eqr ? 0.0 : -rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0);
      xdg(m,2,k,j,i) = eqr ? 0.0 : sb;
      xdg(m,3,k,j,i) = eqr ? 0.0 : gres;
    }
    // (the face sources' divergence and the source term once: each call is inlined)
    const Real sdv = CleanSourceDivergence(x, m, k, j, i);
    const Real hes = eso ? x.h*es(m,k,j,i) : 0.0;
    x.w(m,CLW_B,k,j,i) = eqr ?
        x.w(m,CLW_EO,k,j,i) + x.w(m,CLW_EI,k,j,i) - eg - work
        + x.w(m,CLW_RY,k,j,i)*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - sdv + hes :
        x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
        - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work - sdv + sb + gres + hes;
    // (clean_stag_sbimp: (1 - r) sb(E) in the operator as the cell layout's, Newton form)
    if (x.sbi) x.w(m,CLW_B,k,j,i) -= (1.0 - rr)*CleanStagBetaE(xh, CLW_EK, m, k, j, i);
    // (clean_stag_xpat: the exchange's removal R = r hcke E_k explicit while it leaves B >= th E_old (th = CLN_KPOS_TH),
    // past that, R - (B0 - th E_old) (B0 = B + R) on the diagonal as its ratio to E_k:
    // B >= min(B0, th E_old) >= 0 where B0 >= 0; continuous; the same fixed point)
    if (pdo) {
      Real pd = x.kq >= 0 ? x.w(m,CLW_KM,k,j,i) : 0.0;   // (clean_stag_kpos)
      const Real rx = rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0);
      if (xpd && !eqr && rx > 0.0) {
        const Real b0 = x.w(m,CLW_B,k,j,i) + rx;
        const Real ex = fmax(b0 - CLN_KPOS_TH*fmax(x.w(m,CLW_EO,k,j,i), 0.0), 0.0);
        if (rx > ex) {
          pd += rr*hcke*(1.0 - ex/rx);
          x.w(m,CLW_B,k,j,i) = b0 - ex;
          Kokkos::atomic_add(&kc(0), 1.0);
          Kokkos::atomic_add(&kc(1), 1.0);
        }
      }
      x.w(m,CLW_PD,k,j,i) = pd;
    }
    // (clean_xt_exact: the row's exchange as linearised, B less its rest, and its work)
    if (xte) {
      x.w(m,CLW_XA,k,j,i) = eqr ? 0.0 : x.w(m,CLW_B,k,j,i) - (x.w(m,CLW_EO,k,j,i) - work
          - sdv + hes);
      x.w(m,CLW_XW,k,j,i) = work;
    }
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
  // (a flat range: the 4-D MDRange reduction took 1.3 ms on the A100)
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
    xh.kq = -1;
    const CleanArgs *pxh = CleanDevArgs(xh, 1);
    par_for("cln_sbx", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is, is+nx1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &xh = *pxh;   // (in device memory: no per-thread copy)
      w(m,CLW_EP,k,j,i) -= CleanStagBetaE(xh, CLW_X, m, k, j, i) -
                           CleanStagBetaE(xh, CLW_EK, m, k, j, i);
    });
  }

  // (clean_xt_exact: the closure of the row's linearised exchange at the new T; x2: T
  // from the local 2x2 of the E and gas rows, the E row's rest S = X diag - XA fixed)
  const bool xte = (cln_xte > 0) && !cell && cln_feedback, x2 = xte && cln_xte == 2;
  const bool pdx = x2 && cln_xpat, pxe = xte && cln_xpat;
  const bool kpg = cln_kq_on;   // (clean_stag_kpos: the gas row's pd_g X)
  const bool dmp = cln_pdmp && !cell, dfst = cln_pdfst;
  if (dmp && cln_pdh.extent_int(0) < nmb) {
    Kokkos::realloc(cln_pdh, nmb, cln_w.s.extent_int(2), cln_w.s.extent_int(3),
                    cln_w.s.extent_int(4));
  }
  auto hd = cln_pdh;
  cln_pdfst = false;
  Real change = 0.0, femax = 0.0, ermax = 0.0, xcl = 0.0;
  // (CUDA: 4 blocks of 256 per SM, A100 R136 -1.5 ms/cycle; the reduction keeps its
  // blocks and grid: bitwise)
  using NewtonRp =
      Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, CleanLbMb(4)>>;
  AthReduce("cln_newton", NewtonRp(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dmax, Real &fe_max, Real &er_max, Real &xc) {
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
    const Real xr = (cell || poslim) ? x : fmax(x, 0.0);   // (the change: undamped)
    // (clean_pass_damp: a change that turns back, l = d/d_last < 0, takes d/(1 - l), the
    // fixed point of that mode, as CleanStagXDonor's C; the fixed point is unchanged;
    // an iterate E <= 0 takes the solution whole, so that it leaves E < 0 at once)
    if (dmp) {
      const Real d = x - e_old, dp = dfst ? 0.0 : hd(m,k,j,i);
      const Real l = (d*dp < 0.0 && e_old > 0.0) ? d/dp : 0.0;
      x = e_old + d/(1.0 - l);
      hd(m,k,j,i) = x - e_old;
    }
    const Real ex = fmax(x, 0.0);
    const Real e = (cell || poslim) ? x : ex, g = w(m,CLW_GM,k,j,i);
    // (gflux: the gas takes the Doppler part hc kE dop E of sb at X, as the E row takes
    // it in r; at E_k it lagged a source of hc kE dop, in a moving opaque cell more than
    // the net exchange (1 - r) hc kE, and T and E drifted apart pass by pass; sbimp:
    // the gas takes all of sb at X)
    if (gfd) w(m,CLW_EP,k,j,i) -= w(m,CLW_JV+1,k,j,i)*(x - e_old);
    const Real hcke = hc*op(m,M1_OP_E,k,j,i);
    Real t;
    if (x2) {
      const Real dg = 1.0 + hcke*w(m,CLW_RY,k,j,i) + (pdx ? w(m,CLW_PD,k,j,i) : 0.0);
      const Real sr = x*dg - w(m,CLW_XA,k,j,i);
      t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i),
              hc*op(m,M1_OP_P,k,j,i)*g*arad_/(1.0 + hcke),
              w(m,CLW_EI,k,j,i) + (hcke*sr + w(m,CLW_EP,k,j,i))/(1.0 + hcke), t_old, rtol);
    } else {
      const Real pg = kpg ? w(m,CLW_KM,k,j,i) - w(m,CLW_KF,k,j,i) : 0.0;   // (kpos)
      t = CleanExchangeTemperature(eos, uh(m,IDN,k,j,i), hc*op(m,M1_OP_P,k,j,i)*g*arad_,
              w(m,CLW_EI,k,j,i) + (hcke + pg)*ex + w(m,CLW_EP,k,j,i), t_old, rtol);
    }
    if (xte) {
      const Real dg = 1.0 + hcke*w(m,CLW_RY,k,j,i) +
                      (pxe ? w(m,CLW_PD,k,j,i) : 0.0);
      const Real ga = hc*op(m,M1_OP_P,k,j,i)*g*arad_*SQR(SQR(t)) - hcke*x -
                      w(m,CLW_EP,k,j,i);
      xc = fmax(xc, fabs(w(m,CLW_XA,k,j,i) - (dg - 1.0)*x - ga)/fmax(fabs(x) +
                    fabs(w(m,CLW_EO,k,j,i)) + fabs(w(m,CLW_EI,k,j,i)), 1.0e-300));
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
                                              h*k1(m,M1_T2_E,k,j,i))/se);
    } else {
      dmax = fmax(dmax, fmax(de, dt));
    }
    // staggered: round-off level of E (its terms over E; cell: CleanCellVelocityUpdate)
    if (!cell) fe_max = fmax(fe_max, (fabs(e)*(1.0 + hc*op(m,M1_OP_E,k,j,i)) +
        fabs(w(m,CLW_EO,k,j,i)) + hc*op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(t)))/
        fmax(fabs(e), 1.0e-300));
  }, Kokkos::Max<Real>(change), Kokkos::Max<Real>(femax), Kokkos::Max<Real>(ermax),
     Kokkos::Max<Real>(xcl));
  // (one reduction over the ranks)
  Real mx[4] = {change, femax, ermax, xcl};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, mx, (xte && cln_nlprint) ? 4 : 3, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
#endif
  change = mx[0]; femax = mx[1]; ermax = mx[2]; xcl = mx[3];
  if (xte && cln_nlprint) {   // (diagnostic: CLEAN XT with clean_nl_print only)
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN XT cycle %d t2s %d pass closure %.3e\n", pmy_pack->pmesh->ncycle,
                  t2_solve, xcl);
    }
  }
  if (!cell) cln_fe = std::numeric_limits<Real>::epsilon()*femax;
  if (bud) {   // (the tolerance of the next test and solve; for settled, the level of a
    // change the linear solve fixes at the floor, its 0.1 tfl in the budget's norm)
    // (a step dt_retry redoes, whose state just failed, is solved strictly too)
    const Real tfl = fmax(cln_tol_nl, 10.0*cln_tol);
    cln_tolb = (kok && (!dtr_redo || cln_rtb)) ?
               fmin(fmax(cln_bkap*ermax, tfl), cln_bmax) : tfl;
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
  return change;
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
  xh.kq = -1;
  xh.mcl = true;
  // (the F solve's face tie after the drag, which evaluates the same face T: CLF_PT)
  const bool td = cln_fwork && cln_ftie;
  if (cln_fwork) {   // (clean_stag_fwork: the momentum of the
    CleanStagFSolve(CLW_EK, !td, false);         //  F solve, which the write-back gives
    xh.fs = true;                                //  the gas, for the work W of the E row)
    xh.tf = cln_ftie && !x.cub;
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
  if (fa) CleanStagDrag(false, true, false);   // (the lateral sums of E_k: above)
  if (td) CleanStagFaceTie(CLW_EK, (t2_solve == M1_T2S_NONE) ? 0 : ((t2_solve ==
                           M1_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD), true);
  xh.fa = fa;
  auto xa = x;
  xa.cub = false;
  Real dvmax = 0.0, vmax = 0.0, fvmax = 0.0;
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
             const CleanArgs *pxa = CleanDevArgs(xa, 4);)
  // (with the F solve's F, fa: CleanStagMomentum<true> reads the faces' terms and
  // computes none, in a kernel of its own; on HIP at 256 VGPRs, 2 waves per SIMD)
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
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // (device memory)
      int m, k, j, i;
      cell(idx, m, k, j, i);
      CleanVelocityCell<true>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, true, m,
                              k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  } else {
    AthReduce("cln_velocity", CleanRpLb2(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // (device memory)
      int m, k, j, i;
      cell(idx, m, k, j, i);
      CleanVelocityCell<false>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, false, m,
                               k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  }
  CleanReduceMax(dvmax, vmax, fvmax);
  if (fa) {
    CleanGhostDefer(cln_w, CLW_A1, 3);   // (as CleanAdvectiveCoefficients)
  } else {
    if (cln_vimp) CleanStagDrag(false);   // (with the coefficients of this dm)
    CleanAdvectiveCoefficients();
  }
  vmax = fmax(vmax, 1.0e-300);
  cln_fv = std::numeric_limits<Real>::epsilon()*fvmax/vmax;
  return dvmax/vmax;
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
  CleanFor("cln_copy", nmb1, ks, ke, js, je, is, ie,
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
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real r1 = 0.0, r2 = 0.0;
    Kokkos::parallel_reduce("cln_real", CleanRpLb2(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &q1, Real &q2) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
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
      x.fc[dir](m,CLF_PT,k,j,i) = f;   // (for CleanWriteBack, next; QL is free)
      const Real fr = f - CleanKposFlux(x, CLW_X, dir, m, k, j, i);   // (kpos: Phi's own)
      const Real ed = (fr > 0.0) ? el : er;
      const bool r2 = !(fabs(fr) <= (1.0 + 1.0e-12)*x.c*ed + x.c*zero);
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
//! clean_stag_rayratio (the old production's form: the rays enter only as ratios times
//! the state's E): J -> E in Phi and in the half ranges, and F_ray -> (F_ray/J_d) E_d of
//! the face's donor d (the upwind cell of F_ray; none where J_d <= e_floor: no bracket).
//! The same for E = J; a ray set off the state's scale (a source hole, a lagged ring
//! ratio) cannot reach the operator.

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
  const bool rat = cln_rayratio, thf = cln_thinf;
  const Real rmn = cln_raymin;
  auto iw_ = iw;
  CleanFor("cln_rayj", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x.w(m,CLW_Z,k,j,i) = rat ? iw_(m,M1_IW_EN,k,j,i) : x.cs(m,CLN_J,k,j,i);
  });
  CleanGhostFill(cln_w, CLW_Z, 1);
  CleanStagCsLateral(CLW_Z);
  // (clean_stag_thinf_held: theta = min(1, max(J, P)/(clean_ray_min E_old)) of every
  // cell in CLW_V, a Krylov scratch slot, P the largest rays' flux out of one of its
  // faces: 1 where the rays carry the cell (as its closure, a beam's transit cell
  // included), down to 0 where they hold less than clean_ray_min of its E)
  const bool thh = cln_thinh, thg = cln_thhg;
  if (thh) {
    CleanFor("cln_thheld", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real p = 0.0;
      for (int d = 0; d < x.nf; ++d) {
        const int n = (d == 0) ? CLN_FX : ((d == 1) ? CLN_F2 : CLN_F3);
        const bool up = x.cub && d > 0 && ((d == 1) ? j == x.je : k == x.ke);
        const Real fh = up ? x.cs(m,CLN_FH+4*(d-1),k,j,i) :
                        x.cs(m,n,k+(d == 2),j+(d == 1),i+(d == 0));
        p = fmax(p, fmax(fh, -x.cs(m,n,k,j,i)));
      }
      const Real e = rmn*x.w(m,CLW_EO,k,j,i), h = fmax(x.cs(m,CLN_JR,k,j,i), p);
      // (clean_stag_thinf_held_gate: times r^2, r = E_old/H where the rays bring more than
      // the cell holds, H = max(J, P) > E_old (the old production's quasi-steadiness gate,
      // on the front side only): a dark floor cell a beam just hit moves with its own half
      // range; a cell the rays under-fill keeps C1's theta)
      // weighted by the cell's absorption depth over the stage, tau_a = h c kappa_P: where
      // the gas takes the radiation within the stage the step-start rays do not predict
      // the cell being filled; a free-streaming front (tau_a ~ 0, a beam) keeps C1's theta)
      const Real eo = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
      const Real r = (thg && h > eo) ? eo/h : 1.0;
      const Real ta = thg ? fmin(x.h*x.c*x.op(m,M1_OP_P,k,j,i), 1.0) : 0.0;
      x.w(m,CLW_V,k,j,i) = ((e > h) ? h/e : 1.0)*(1.0 - (1.0 - r*r)*ta);
    });
    CleanGhostFill(cln_w, CLW_V, 1);
    CleanStagCsLateral(CLW_V);
  }
  CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int n = (dir == 0) ? CLN_FX : ((dir == 1) ? CLN_F2 : CLN_F3);
    CleanFor("cln_sbracket", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
      // (Phi(J) at its steady value, -(c/chi_f) div(D J) in the share 1 - w: the bracket,
      // and with it the steady state of the scheme, does not depend on h or on G.  The
      // diffusion model holds where its flux is at most the free-streaming c J_f, lambda =
      // |div(D J)|/(chi_f J_f) <= 1; where the face is thinner (an outer envelope's floor
      // gas at w < 1, chi_f 7e-22: c/chi_f made Phi/c ~ 5e3 J, h- ~ -9e3, E < 0) the model
      // takes that bound, c J_f along -div(D J) (Larsen's max limiter; Wilson's
      // c/(chi_f + |div(D J)|/J_f) also changed it by 1/(1 + lambda) where it is valid)
      const Real w = x.fc[dir](m,CLF_W,k,j,i);
      const Real chi = 0.5*(x.op(m,M1_OP_T,k,j,i) + x.op(m,M1_OP_T,k-dk,j-dj,i-di));
      const Real st = CleanFaceStress(x, CLW_Z, dir, m, k, j, i, true);
      const Real jf = 0.5*(x.w(m,CLW_Z,k,j,i) + x.w(m,CLW_Z,k-dk,j-dj,i-di));
      const Real cd = (jf > 0.0 && fabs(st) > chi*jf) ? x.c*jf/fabs(st) :
                      ((chi > 0.0) ? x.c/chi : x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c);
      // (the cells' half ranges realisable, h+ in [0, 1], h- in [-1, 0], as the old
      // production's; the rest of the rays' flux is DG below)
      x.fc[dir](m,CLF_HP,k,j,i) = fmin(fmax(x.fc[dir](m,CLF_HP,k,j,i), 0.0), 1.0);
      x.fc[dir](m,CLF_HM,k,j,i) = fmax(fmin(x.fc[dir](m,CLF_HM,k,j,i), 0.0), -1.0);
      Real gl = 0.0, gr = 0.0;
      CleanStagArea(x, dir, m, k, j, i, w, gl, gr);
      const Real pj = -(1.0 - w)*cd*st +
          w*x.c*(x.fc[dir](m,CLF_HP,k,j,i)*x.w(m,CLW_Z,k-dk,j-dj,i-di) +
                 x.fc[dir](m,CLF_HM,k,j,i)*x.w(m,CLW_Z,k,j,i)) +
          gl*x.w(m,CLW_Z,k-dk,j-dj,i-di) + gr*x.w(m,CLW_Z,k,j,i);
      // (cubed sphere: the high x2, x3 face of a MeshBlock from its last cell, CLN_FH)
      const bool up = x.cub && dir > 0 && (j > x.je || k > x.ke);
      const Real fr0 = up ? x.cs(m,CLN_FH+4*(dir-1),k-dk,j-dj,i-di) : x.cs(m,n,k,j,i);
      const bool bf = CleanBoundaryFace(x, dir, m, k, j, i);
      if (bf) return;
      // (the rays' J also where it is below clean_ray_min E: their face flux passes
      // through cells they hold almost nothing in, the side cells of a diagonal beam;
      // rayratio: E, and the rays' flux per J of its donor times the donor's E)
      const Real jl = rat ? x.w(m,CLW_Z,k-dk,j-dj,i-di) : x.cs(m,CLN_JR,k-dk,j-dj,i-di);
      const Real jr = rat ? x.w(m,CLW_Z,k,j,i) : x.cs(m,CLN_JR,k,j,i);
      const Real jd = (fr0 > 0.0) ? x.cs(m,CLN_JR,k-dk,j-dj,i-di) : x.cs(m,CLN_JR,k,j,i);
      if (rat && !(jd > efl)) {
        x.fc[dir](m,x.fdg,k,j,i) = 0.0;
        return;
      }
      const Real fy = rat ? fr0/jd*((fr0 > 0.0) ? jl : jr) : fr0;
      // (clean_stag_thinf: the rest of the rays' face flux, dd = F_y - h+ J_L - h- J_R, in
      // the face half ranges, tilted in proportion to their sizes (the shrinking one at
      // most to 0, the rest then on the other), per J of the donor, at least clean_ray_min
      // E_old (a cell the rays hold almost nothing in), within g+ in [0, 1], g- in [-1, 0]:
      // the thin share w c (g+ E_L + g- E_R) is linear in E with the upwind signs and
      // within [-c E_R, c E_L], realisable by construction; no DG.  A face flux the cap
      // cuts (a beam narrower than a cell, F_y > J_L) is carried by a higher E)
      if (thf) {
        const Real hl = x.fc[dir](m,CLF_HP,k,j,i)*jl, hr = x.fc[dir](m,CLF_HM,k,j,i)*jr;
        const Real dd = fy - (hl + hr), fr = hl + hr + dd, s = hl - hr;
        const Real tl = hl + dd*hl/s, tr = hr - dd*hr/s;
        const bool cut = (tl < 0.0) || (tr > 0.0);
        const Real hp = !(s > 0.0) ? hl : (cut ? ((tl < 0.0) ? 0.0 : fmax(fr, 0.0)) : tl);
        const Real hm = !(s > 0.0) ? hr : (cut ? ((tl < 0.0) ? fmin(fr, 0.0) : 0.0) : tr);
        const Real el = fmax(x.w(m,CLW_EO,k-dk,j-dj,i-di), 0.0);
        const Real er = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
        // (clean_stag_thinf_held: the share theta of the donor's E the rays carry moves
        // with their flux factor, the rest with the donor's own half range: a convex
        // blend, realisable; a ray-shadowed donor (J << E) is no dam)
        const Real gp = fmin(hp/fmax(jl, rmn*el), 1.0), gm = fmax(hm/fmax(jr, rmn*er), -1.0);
        const Real ql = thh ? x.w(m,CLW_V,k-dk,j-dj,i-di) : 1.0;
        const Real qr = thh ? x.w(m,CLW_V,k,j,i) : 1.0;
        auto &fp = x.fc[dir](m,CLF_HP,k,j,i), &fm = x.fc[dir](m,CLF_HM,k,j,i);
        if (jl > 0.0) fp = thh ? ql*gp + (1.0 - ql)*fp : gp;
        if (jr > 0.0) fm = thh ? qr*gm + (1.0 - qr)*fm : gm;
        x.fc[dir](m,x.fdg,k,j,i) = 0.0;
        return;
      }
      // (the old production's form, M1ThinFace implicit_thin_corr = hc: the cells' half
      // ranges h+- realisable, the rest of the rays' face flux the lagged flux DG/w = q c
      // (F_y - h+ J_L - h- J_R) of the face (CleanThinFaceDG), F_y clamped to [-J_R, J_L]
      // (|F| <= c J of the donor), q = r^2 the quasi-steadiness gate, r the worse of the
      // two cells' min(E_old/J, J/E_old) (E of the old vector, ghosts filled); tilted into
      // h+- the defect reached F/J_R ~ 1e3 .. 4e20 at a ray-shadowed donor below a lit
      // cell of the evolved HG envelope hegiant.00184)
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
#endif  // M1_TU_PART(1)

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanDtrD
//! \brief dt_retry: the memory d and its slope K1 (CLF_N + 3, CLF_KD) are the only face
//! values a cycle takes from the last one; the checkpoint keeps them in cln_rstd (save,
//! the restart's layout) instead of the face arrays

#if M1_TU_PART(2)
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
#endif  // M1_TU_PART(2)

//----------------------------------------------------------------------------------------
//! \struct CleanFLine
//! \brief the F solve's lines on every mesh (CleanStagFSolve): the cell terms t of cell
//! (k,j,i) at its face (kf,jf,iff) normal to dir, of the cell's face family along the
//! line (al, the lagged pressure p, w, h+, h- and the area pressure pa), and the F of
//! the line's faces.  Spherical polar and Cartesian: t in CLW_L + 3 dir, F of each
//! cell's low face in CLW_Y + dir (CLW_S); cubed sphere: cln_v in the owner's axes
//! (CleanCsSlot: across a seam an x2 line may continue as an x3 line, so a cell keeps
//! both faces of each family; cell terms of family q: 3 q + 0..2, 15 + q; F: 9 + 2 q +
//! hi)

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
  // the slot of the F of face (kf,jf,iff) (far: of the other face of the family of
  // cell (k,j,i) along it) in that cell's storage, s its orientation along ours
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
  // F of the other face of the family of cell (k,j,i) along its face (kf,jf,iff), in
  // our axes (spherical polar, Cartesian: a physical face's F^H); Near: of that face
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
  // whether the MeshBlock solves face (k,j,i), number t of the n + 1 of its line along
  // dir (tie: clean_stag_ftie, the block whose low face it is; on the cubed sphere, with
  // clean_stag_fschur, the lower gid's, as CleanStagFCsShared), else its identity row
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
//! \brief the row lo F_{f-1} + dg F_f + hi F_{f+1} = r of face (k,j,i) normal to dir in
//! the F solve (CleanStagFSolve); false where the face has no equation (F = F^H) or (own
//! false, clean_stag_ftie) another block owns it: the identity row, r the caller's;
//! sm (clean_stag_fschur): the row sum lo + dg + hi, without the cancellation.  The
//! geometry: the kinetic stress over h c is q (P_R - P_L) + qs (P_R + P_L), spherical
//! polar's area form (1/dx +- s) (CleanStagSph), Cartesian 1/dx, the cubed sphere's
//! (A_R P_R - A_L P_L)/V_f (the cells' area pressures pa, the lagged one a gradient
//! A_f/V_f)

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
    // (fkin: h c^2 [w P_kin]_L^R/dx, P_kin = (h+ F_lo + h- F_hi)/c of a cell)
    Real q, qs = 0.0, ag = 1.0;
    if (x.cub) {
      q = x.h*x.c/(0.5*(x.g(m,CLG_F+3,k,j,i) + x.g(m,CLG_F+3,k-dk,j-dj,i-di)));
      ag = CleanCsNormA(x, dir, m, k, j, i);
    } else {
      q = x.h*x.c*(x.sph ? CleanStagSph(x, 0, dir, m, k, j, i) : x.idx[dir]);
      qs = x.sph ? x.h*x.c*CleanStagSph(x, 4, dir, m, k, j, i) : 0.0;
    }
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
//! \fn void CleanStagFOwn
//! \brief clean_stag_fschur: the line (m, a, b) along dir of a MeshBlock (n cells) as a
//! segment s of the F line: its first own face x = s0 + s1 yL + s2 xR and its last own
//! face y = s4 + s5 yL + s6 xR, yL the face before the first (another block's or the
//! wrap's last own), xR the face after the last (the next block's first own: the shared
//! high face, the identity row, or on the cubed sphere, where the lower gid owns a shared
//! face (CleanFLine::Own), the one beyond it), and s3 = 1 - s1 - s2,
//! s7 = 1 - s5 - s6 from the row sums (the response to a uniform F: near 0 on a nearly
//! singular line, h c/dx >> 1/theta, which 1 - s1 - s2 would lose); by a Thomas sweep
//! up for y and one down for x

template <class FC>
KOKKOS_INLINE_FUNCTION
void CleanStagFOwn(const CleanArgs &x, const CleanFLine &ln, const FC &fc, const int dir,
                   const int m, const int a, const int b, const int n, Real s[8]) {
  auto face = [&](const int t, int &k, int &j, int &i) {
    k = (dir == 2) ? x.ks + t : x.ks + b;
    j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
    i = (dir == 0) ? x.is + t : x.is + a;
  };
  int k, j, i;
  Real lo, hi, dg, r, sm;
  face(0, k, j, i);
  const int t0 = ln.Own(x, dir, m, k, j, i, 0, n, true) ? 0 : 1;
  face(n, k, j, i);
  const bool bh = CleanBoundaryFace(x, dir, m, k, j, i);   // (F^H)
  const int t1 = (x.cub && !bh && ln.Own(x, dir, m, k, j, i, n, n, true)) ? n : n - 1;
  const Real fh = bh ? fc(m,CLF_N,k,j,i) : 0.0, wh = bh ? 1.0 : 0.0;
  Real cp = 0.0, dc = 0.0, dl = 1.0, dw = 0.0;   // (F_{-1} = yL)
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
  Real cq = 0.0, ec = fh, er = 1.0 - wh, ew = wh;   // (F_{t1+1} = xR or F^H)
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
//! \brief clean_stag_fschur: the segment c of two adjacent ones, p then q (CleanStagFOwn),
//! and p's last face ya = ya0 + ya1 yL + ya2 xR and q's first xb (same) in the faces
//! beyond c's ends; 1 - p6 q1 and c's s3, s7 as sums of the parts' (no cancellation)

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
//! |F^H| > c E_donor (with fkin only under clean_stag_flim: al on the (1 - w) share).
//! al = 0 gives F^H (a face with al_L = al_R = 0 keeps it exactly).
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
//! One solve on every mesh, its lines through CleanFLine: on the cubed sphere the cell
//! terms of a line (al, the lagged pressure, w and the half ranges) are those of the
//! cell's face family along it, D about the family's mean normal
//! (CleanCsFamily), the kinetic stress w (F_kin/c - s^2 J) of the family's normal, 0 for
//! an isotropic field, in area form, (A_R P_R - A_L P_L)/V_f, A the mean area of the
//! cell's two faces of the family, so that free streaming F keeps r^2 F, and w D (J - E)
//! a gradient (A_f/V_f); its line ends iterate the seams (up to 64 sweeps, to 1e-14 of
//! the largest F), and a face two MeshBlocks share takes the solution of the one with
//! the lower gid, so that both see the same F to the bit.  Its lateral stress keeps D E
//! (no transverse kinetic stress: its lateral Gauss sums are also the non-orthogonal
//! correction of the normal stress, which a kinetic transport of F does not reproduce,
//! and a static envelope would lose its balance on the panel-tangent faces).

#if M1_TU_PART(2)
// (the F slots 9 + 2 q + hi of cln_v that the ghosts of a face of family q take
// (CleanGSeg::fm): its own family's, CleanFLine; f0 = 2: without the x1 family's)
static std::array<uint32_t, 4> CleanFMask(const int f0) {
  return (f0 > 0) ? std::array<uint32_t, 4>{0, 3, 014, 0}
                  : std::array<uint32_t, 4>{3, 014, 060, 0};
}

void RadiationM1::CleanStagFSolve(const int e, const bool tie_now, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // (lat false: the caller's are those of e, CLW_L as they were)
  if (lat) CleanStagCsLateral(e);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  xh.kq = -1;
  auto op = opac;
  const CleanFLine ln{cln_v, x.cub};
  const bool fk = cln_fkin, dok = cln_dok, cub = x.cub, flm = cln_flim;
  if (flm && cln_flimc.extent_int(0) != 2) Kokkos::realloc(cln_flimc, 2);
  if (flm) Kokkos::deep_copy(cln_flimc, 0.0);
  auto flc = cln_flimc;
  const int kd = (t2_solve == M1_T2S_NONE) ? 0 :
                 ((t2_solve == M1_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD);
  const Real c1lo = cln_c0, c1hi = cln_c1;
  // (the cubed sphere with clean_stag_fs_pcr_cs and fschur: its x1 lines, in panel, too)
  const bool fpcr = cln_fspcr && (!cub || (cln_fspcs && cln_fschur)) &&
      !std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>;
  // (the F of every direction, CLW_Y, CLW_Z, CLW_S or cln_v: one exchange; part 1:
  // started, the last one, which x1's F and dissipation term overlap)
  // (cubed sphere, one MeshBlock along x1: the x1 family's slots (0..2, 9, 10, 15) have no
  // line across a MeshBlock face and their ghosts are not read: not exchanged)
  const int f0 = (cub && pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel)
                 ? 2 : 0;
  auto yfill = [&](const int part) {
    if (cub) {   // (CleanFLine reads the faces' own cells, of the face's family, fm)
      CleanGhostFill(cln_v, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));
      return;
    }
    if (x.nf == 3) CleanGhostDefer(cln_w, CLW_S, 1);
    CleanGhostFill(cln_w, CLW_Y, std::min(x.nf, 2), part);
  };
  // the directions in three phases (F^H, al; the line solves; F), each direction with
  // its own al and F, so that each exchange serves all of them (clean_stag_fsweeps N:
  // the line solves up to N times, each with the line ends' F of the last, so that the
  // lines converge across MeshBlock faces and the periodic wrap)
  // (the cubed sphere without clean_stag_fschur iterates its seams: up to 64 sweeps)
  const int nsw = (cub && !cln_fschur) ? 64 : std::max(cln_fsweeps, 1);
  const bool tie = cub ? cln_fschur : cln_ftie;
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
        const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
        fc(m,CLF_N,k,j,i) = CleanFaceFlux(xh, e, dir, m, k, j, i, true);
        fc(m,CLF_PT,k,j,i) = fc(m,CLF_N,k,j,i);   // (F^H for CleanStagDrag; QL is free)
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
      auto v = cln_v;
      auto f1 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rc;
        const Real tau = op(m,M1_OP_T,k,j,i)*CleanCellFaceDistance(x, dir, m, k, j, i);
        // (fkin: the streaming share is kinetic, without dissipation)
        Real a = fk ? 0.0 : CleanDissipationWeight(CLA_TAU, tau, 0.0, 0.0, rc);
        for (int s = 0; s < 2; ++s) {
          const int kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
          if (!fk) a = fmax(a, fc(m,CLF_W,kf,jf,iff));
          const Real f = fc(m,CLF_N,kf,jf,iff);
          const Real ed = (f > 0.0) ? x.w(m,e,kf-dk,jf-dj,iff-di) : x.w(m,e,kf,jf,iff);
          // (clean_stag_flim with fkin: the dissipation of a face's excess over its
          // donor's c E, a stress, telescopes: momentum and L_z exact)
          if ((!fk || flm) && !CleanBoundaryFace(x, dir, m, kf, jf, iff) &&
              !(fabs(f) <= (1.0 + 1.0e-12)*x.c*fmax(ed, 0.0))) a = 1.0;
          // (the cubed sphere's faces: the last F solve's)
          if (cub) v(m,9+2*dir+s,k,j,i) = fc(m,CLF_N+2,kf,jf,iff);
        }
        // (fkin: w of the cell; the part of P - D E without F, -w D E, and the stress
        // the kinetic form misses for the rays' own state, F = c s J: (D - s^2) J, 0 for
        // a beam, J/3 for an isotropic field, with J of the rays where the cell holds at
        // least that much, lagged; the cubed sphere: D about the family's mean normal,
        // s^2 J in area form)
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
          // (F of the low face, for the exchange; multi-D: the last F solve's)
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
      // (the cubed sphere, one MeshBlock along x1: one exchange, depth 1, the ghosts of a
      // face its family's terms (fm), started once the last direction's cells it sends
      // are done, in flight behind its inner cells)
      const bool pv = cub && f0 > 0 && cln_hmax >= 12 && dir == x.nf - 1;
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
    // (clean_stag_fschur: the faces beyond the line ends exact, not lagged)
    if (cln_fschur && q == x.nf) CleanStagFSchur();
    const int n = (dir == 0) ? ie - is + 1 : ((dir == 1) ? je - js + 1 : ke - ks + 1);
    const int na = (dir == 0) ? je - js + 1 : ie - is + 1;
    const int nb = (dir == 2) ? je - js + 1 : ke - ks + 1;
    Real dmax = 0.0, fmx = 0.0;
    if (ph == 1) {
      // per face (coalesced) lo (CLW_AM), hi (CLW_AP), dg (CLW_A0, then the Thomas
      // factor den), r (CLW_R0) and the faces with an equation (CLW_P; scratch slots of
      // the E solve, which the directions share: every pass; the x1 PCR reads them too),
      // and in the line sweep the Thomas factors den and c' (CLW_CP)
      CleanFor("cln_fs2f", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const int t = (dir == 0) ? i - is : ((dir == 1) ? j - js : k - ks);
        // (clean_stag_ftie: the high face of a line that another MeshBlock or the
        // periodic wrap shares is its owner's, the block whose low face it is: an
        // identity row, the value its CLW_Y ghost brings, CleanStagFaceTie)
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
      // (clean_stag_fs_pcr: on a GPU the x1 lines by parallel cyclic reduction, one team
      // per line coalesced along x1, as the preconditioner's (CleanStagLineSolve); the
      // same system, F to round-off of the Thomas solve)
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
          // (the rows of cln_fs2f; the faces beyond a MeshBlock line end: the
          // neighbours', lagged; ftie: a shared high face its owner's, the identity row)
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
      // one exact Thomas solve per line; the faces beyond a line end of a MeshBlock or
      // a periodic wrap are the last F solve's (lagged); on the cubed sphere the largest
      // change of a face and the largest F (a reduction only there: no fence elsewhere)
      auto line = KOKKOS_LAMBDA(const int idx, Real &d_max, Real &f_max) {
        const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
        auto face = [&](const int t, int &k, int &j, int &i) {
          if (dir == 0) {i = is + t; j = js + a; k = ks + b;}
          if (dir == 1) {j = js + t; i = is + a; k = ks + b;}
          if (dir == 2) {k = ks + t; i = is + a; j = js + b;}
        };
        Real dp = 0.0, cp = 0.0;
        for (int t = 0; t <= n; ++t) {
          int k, j, i;
          face(t, k, j, i);
          Real den = x.w(m,CLW_A0,k,j,i);   // (dg)
          // (the faces beyond a MeshBlock line end: the neighbours', lagged)
          const bool eq = x.w(m,CLW_P,k,j,i) > 0.0;
          den -= ((eq && t == 0) ? 0.0 : x.w(m,CLW_AM,k,j,i))*cp;
          cp = ((eq && t == n) ? 0.0 : x.w(m,CLW_AP,k,j,i))/den;
          x.w(m,CLW_CP,k,j,i) = cp;
          Real r = (x.w(m,CLW_P,k,j,i) < 0.0) ?
              ln.Near(x, dir, m, k, j, i, k-(t == 0)*dk, j-(t == 0)*dj, i-(t == 0)*di) :
              x.w(m,CLW_R0,k,j,i);
          if (x.w(m,CLW_P,k,j,i) > 0.0) {
            if (t == 0) r -= x.w(m,CLW_AM,k,j,i)*ln.Far(x, dir, m, k, j, i, k-dk, j-dj,
                                                         i-di);
            if (t == n) r -= x.w(m,CLW_AP,k,j,i)*ln.Far(x, dir, m, k, j, i, k, j, i);
          }
          dp = (r - ((t == 0) ? 0.0 : x.w(m,CLW_AM,k,j,i))*dp)/den;
          x.w(m,CLW_R,k,j,i) = dp;
        }
        Real fn = 0.0;
        for (int t = n; t >= 0; --t) {   // (the last face's slot is a ghost: the
          int k, j, i;                   //  neighbour's or a physical face)
          face(t, k, j, i);
          fn = x.w(m,CLW_R,k,j,i) - ((t < n) ? x.w(m,CLW_CP,k,j,i)*fn : 0.0);
          if (cub) {
            d_max = fmax(d_max, fabs(fn - ln.Get(x, dir, m, k, j, i, t)));
            f_max = fmax(f_max, fmax(fabs(fn), 1.0e-8*x.c*fabs(x.w(m,e,k,j,i))));
          }
          ln.Put(x, dir, m, k, j, i, t, n, fn);
        }
      };
      // (one warp or wavefront per block: the x1 lines are few, one thread each,
      // latency bound)
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
      // (sweeps: the cubed sphere until no face changes by 1e-14 of the largest F, else
      // until the high faces' F, their owners', changes by < 1e-12 of itself)
      const bool done = nsw > 1 &&
          (cub ? CleanReduce(fsd, true) <= 1.0e-14*CleanReduce(fsf, true)
               : (CleanStagFSweepChange() < 1.0e-12 && q >= 2*x.nf));
      if (done) q = (1 + nsw)*x.nf - 1;
    }
    // (the cubed sphere: a face shared with another MeshBlock takes the F of the lower
    // gid, and the cells' slots of it)
    if (cub && q == (1 + nsw)*x.nf) CleanStagFCsShared();
    if (ph == 2) {
      // F and the dissipation term of every face
      const CleanArgs *px = CleanDevArgs(x, 0);
      auto f3 = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
        const Real f = cub ? fc(m,CLF_N+2,k,j,i) : x.w(m,CleanFLine::Cy(dir),k,j,i);
        fc(m,CLF_N+2,k,j,i) = fc(m,CLF_N,k,j,i);   // (F^H: the boundary's, kept)
        if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
        // (clean_stag_flim: the solution F within 5/4 of its donor's realisable range,
        // [-c E_R, c E_L]: in floor gas, chi ~ 1e-19, the lateral curl of a ray D at odds
        // with the moment F drives F secularly, only chi damps it (|F| ~ 2e4 c E,
        // T19HOT), and the dissipation above, on the (1 - w) share, does not reach it.
        // The 1/4 is the bound's own discretisation (a face against its donor's mean E;
        // free streaming in fs: up to 0.12): a resolved or streaming flow keeps its F.
        // The cut f - fs moved F's equation by -(f - fs)/theta; the gas term
        // (F - G)/c^2 + h div P - h CLF_N+1 then takes h chi (F - A)/c of fs, the
        // absorption of the realisable flux, not the cut: CLF_N+1 -= (f - fs)/(theta h
        // c^2).  Its radiation momentum leaves the total: counted and printed.  E and
        // the gas energy stay exact (flux form; the work is that of the gas term).  The
        // kinetic stress below keeps the solution's f, of the equation it satisfied)
        Real fs = f;
        if (flm) {
          fs = fmin(fmax(f, -1.25*x.c*fmax(x.w(m,e,k,j,i), 0.0)),
                    1.25*x.c*fmax(x.w(m,e,k-dk,j-dj,i-di), 0.0));
          if (fs != f) {
            const Real it = 1.0/fc(m,CLF_TH,k,j,i);
            fc(m,CLF_N+1,k,j,i) -= (f - fs)*it/(x.h*x.c*x.c);
            Kokkos::atomic_add(&flc(0), 1.0);
            Kokkos::atomic_add(&flc(1), fabs(f - fs)*it*0.5*(CleanVol(x, m, k, j, i) +
                               CleanVol(x, m, k-dk, j-dj, i-di))/(x.c*x.c));
          }
        }
        // (fkin, the solution's F solve: d = Phi^H - F, Phi^H = F^H + theta d the energy
        // flux's high-order part, whose memory is G + d (CleanFaceSource): its own
        // history, so that its steady state, F^H of G = F^H, does not depend on h)
        if (fk && e == CLW_X && dok) {   // (and its slope, CleanDInc)
          const Real dd = fc(m,CLF_N,k,j,i) - fs - (1.0 - fc(m,CLF_TH,k,j,i))*
                          fc(m,CLF_N+3,k,j,i);
          fc(m,CLF_N+3,k,j,i) += dd;
          if (kd > 0) fc(m,kd,k,j,i) = dd/x.h;
        }
        fc(m,CLF_N,k,j,i) = fs;
        fc(m,CLF_N+2,k,j,i) = fs;
        Real tl[6], tr[6];
        ln.Cell(x, dir, m, k, j, i, k-dk, j-dj, i-di, tl);
        ln.Cell(x, dir, m, k, j, i, k, j, i, tr);
        if (!(tl[0] > 0.0 || tr[0] > 0.0 || tl[2] > 0.0 || tr[2] > 0.0)) return;
        // (a physical high face keeps F^H; the ghost fill replaced its slot)
        const Real fm = ln.Far(x, dir, m, k, j, i, k-dk, j-dj, i-di);
        const Real fp = ln.Far(x, dir, m, k, j, i, k, j, i);
        auto pk = [&](const Real t[6], const Real l, const Real h) {
          return t[2]*(t[3]*l + t[4]*h)/x.c + (cub ? t[5] : t[1]);
        };
        const Real pl = pk(tl, fm, f), pr = pk(tr, f, fp);
        // (the kinetic stress: the cubed sphere's area form over V_f, else that of the
        // F row, CleanStagFRow)
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
      // (x1: the faces that read no ghost while the exchange of F finishes)
      const bool ov = fov && dir == 0;
      const std::function<void()> fin = [&]() {
        if (ov) CleanGhostFill(cln_w, CLW_Y, 1, 2);
      };
      CleanOverlapFor("cln_fs3", pmy_pack, ov, f3, fin, dir);
    }
  }
  if (cln_ftie && tie_now) CleanStagFaceTie(e, kd);
  if (flm) {   // (clean_stag_flim: the faces cut and their radiation momentum, booked)
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_flimc);
    const Real n = CleanReduce(h(0), false), p = CleanReduce(h(1), false);
    if (global_variable::my_rank == 0 && n > 0.0) {
      std::printf("CLEAN FLIM cycle %d t2s %d faces %.0f momentum cut sum |dF| V/c^2 "
                  "%.6e\n", pmy_pack->pmesh->ncycle, t2_solve, n, p);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagFaceTie
//! \brief clean_stag_ftie (spherical polar, Cartesian): a face shared by two MeshBlocks,
//! or by one with itself across a periodic wrap, takes the values of its owner, the block
//! whose low face it is, as its F does (the F solve's identity row, CleanStagFSolve):
//! the momentum exchange T = (F - G)/c^2 + h div P - h N1 of each face
//! (CleanMomentumChange), the energy flux's memory d and its slope.  The two sides of a
//! seam then give the gas the two halves of one T, its face's absorption where the F
//! solve converged; across the theta-periodic wrap of a wedge, whose two copies have
//! opposite cot(theta) terms, the owner's.  Each block computes them for its low faces
//! into CLW_L + 3 dir (+ 1 d, + 2 slope; free after the F solve), the exchange hands the
//! high faces' to their other side.  The cubed sphere: T in CLF_PT (pt: the drag's),
//! each face's three of the lower gid (CleanStagCsTie).

void RadiationM1::CleanStagFaceTie(const int e, const int kq, const bool pt) {
  // (d's slope: clean_stag_fkin only; without it the face arrays end before CLF_KD)
  const int kd = cln_fkin ? kq : 0;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.fs = true;
  const int nf = x.nf;
  if (x.cub) {
    for (int d = 0; d < nf && !pt; ++d) {   // (as CleanStagDrag; none at physical faces)
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
    CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
    for (int d = 0; d < nf; ++d) {
      const Real f = x.fc[d](m,CLF_N,k,j,i), g = x.fc[d](m,CLF_G,k,j,i);
      // (pt: CleanStagDrag's of the same F just took it, CLF_PT, but on physical faces)
      x.w(m,CLW_L+3*d,k,j,i) = (pt && !CleanBoundaryFace(x, d, m, k, j, i)) ?
          x.fc[d](m,CLF_PT,k,j,i) : (f - g)/(x.c*x.c) + x.h*CleanFaceStress(x, e, d, m, k,
                               j, i, true) - x.h*x.fc[d](m,CLF_N+1,k,j,i);
      x.w(m,CLW_L+3*d+1,k,j,i) = x.fc[d](m,CLF_N+3,k,j,i);
      x.w(m,CLW_L+3*d+2,k,j,i) = (kd > 0) ? x.fc[d](m,kd,k,j,i) : 0.0;
    }
  });
  CleanGhostFill(cln_w, CLW_L, 3*nf);
  for (int d = 0; d < nf; ++d) {   // (the high faces of the block: their owner's d)
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
//! \fn Real RadiationM1::CleanStagFSweepChange
//! \brief clean_stag_fsweeps: the largest change, relative to itself, of the F the high
//! faces of the F solve's lines took from their neighbours (CLW_Y, CLW_Z, CLW_S ghosts)
//! since the last sweep (kept in CLW_TM + dir, free in the write-back; the first sweep
//! only stores it: CLW_TM holds another solve's or the arena's bytes before)

Real RadiationM1::CleanStagFSweepChange() {
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
//! \fn void RadiationM1::CleanStagFSchur
//! \brief clean_stag_fschur (uniform mesh; with clean_stag_ftie, the cubed sphere with
//! its owner rule, CleanFLine::Own): the F solve's lines exact across MeshBlocks and
//! periodic wraps, on the cubed sphere round the rings of four panels' lateral lines.
//! The line of each MeshBlock is a segment of the whole line, its first and last face
//! affine in the faces beyond its ends (CleanStagFOwn).  Rounds of exchanges join the
//! segments (CleanStagFJoin): in each, a MeshBlock sends to its right neighbour its own
//! joined to the segment it got from the left (and the mirror to the left), until with
//! its own they span the line, or once round a periodic one.  The faces beyond its ends
//! then follow from that system of two unknowns and take the place of the lagged ones
//! (the CLW_Y, CLW_Z, CLW_S ghosts) in the line solves (Thomas or PCR).  A segment goes
//! out in three slots of three cells at each end of the line, which the exchange gives
//! the neighbour's ghosts: the six slots the F solve leaves free serve two directions at
//! a time.  Across a cubed-sphere seam the neighbour's line may be its other lateral
//! family and run the other way (CleanCsSlot): its segment is read from that family's
//! slots and reversed, so both ends send for the same rounds.

void RadiationM1::CleanStagFSchur() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const CleanFLine ln{cln_v, x.cub};
  const int nn[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int nm[3] = {pm->mesh_indcs.nx1, pm->mesh_indcs.nx2, pm->mesh_indcs.nx3};
  const int fs[6] = {CLW_TM, CLW_TM+1, CLW_TM+2, CLW_TM+3, CLW_BW, CLW_PW};
  // (the rounds of the segments to the left and to the right: a line of nb MeshBlocks
  // with ends nb - 1 each, a periodic one nb/2 and (nb - 1)/2; the group of directions
  // whose exchanges go together, and the slots of a direction in it)
  // (the cubed sphere: a lateral line is a ring round four panels, its x2 and x3 lines
  // in one group, as a seam may turn one into the other)
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
  for (int g = 0; g == 0 || 2*g < ngd; ++g) {
    int nrt = 0, nd = 0;
    for (int d = 0; d < x.nf; ++d) {
      if (gr[d] == g && nr[d][0] > 0) nrt = std::max(nrt, nr[d][0]), ++nd;
    }
    for (int r = 1; r <= nrt + 1; ++r) {
      const bool fin = (r == nrt + 1);
      for (int dir = 0; dir < x.nf; ++dir) {
        const int rl = nr[dir][0], rr = nr[dir][1];
        if (gr[dir] != g || (!fin && r > rl && r != rr + 1)) continue;
        const int n = nn[dir], na = (dir == 0) ? nn[1] : nn[0];
        const int nb = (dir == 2) ? nn[1] : nn[2];
        auto fc = x.fc[dir];
        const int c0 = fs[gs[dir]], c1 = fs[gs[dir]+1], c2 = fs[gs[dir]+2];
        const int g1 = gs[1], g2 = gs[2];
        const bool pd = per[dir];
        const int rs = cub ? rl : rr;   // (rounds the low end sends)
        AthFor("cln_fsch", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                             (nmb1 + 1)*na*nb),
        KOKKOS_LAMBDA(const int idx) {
          const int m = idx/(na*nb), a = idx%na, b = (idx/na)%nb;
          auto cell = [&](const int t, int &k, int &j, int &i) {
            k = (dir == 2) ? x.ks + t : x.ks + b;
            j = (dir == 1) ? x.js + t : ((dir == 2) ? x.js + b : x.js + a);
            i = (dir == 0) ? x.is + t : x.is + a;
          };
          // (number q of a segment: slot q/3, cell t0 + sg (q%3))
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
          // (a segment in the ghosts beyond the low or high end; on the cubed sphere in
          // the slots of the neighbour's family along the line, reversed, F -> -F, where
          // it runs the other way: CleanCsSlot)
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
          // (the own segment once, kept for the rounds in cells 0, 1 of four scratch
          // slots of the E solve, a set for each direction of the group)
          for (int q = 0; q < 8 && r > 1; ++q) o[q] = keepo(q);
          if (r == 1) CleanStagFOwn(x, ln, fc, dir, m, a, b, n, o);
          for (int q = 0; q < 8 && r == 1 && !fin; ++q) keepo(q) = o[q];
          // (the segment to the left, in the ghosts; to the right, kept after rr: at the
          // start, on the cubed sphere in the line's faces' CLF_N + 2, free here, as both
          // ends send for rl rounds there, for a neighbour that runs the other way)
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
          // the whole line (once round a periodic one, whose ends are each other's: its
          // determinant (1 - c2)(1 - c5) - c1 c6 as a sum)
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
          // (the faces beyond the own ones: the next block's first own face xR and the
          // last block's last yL, in the slots the line solve reads, CleanFLine)
          Real sg;
          int kf, jf, iff;
          if (!bh) {   // (cell n)
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
      if (r <= nrt) {   // (the ghosts of a face take its line family's slots, CleanGSeg::fm)
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
//! \fn void RadiationM1::CleanStagFCsShared
//! \brief cubed sphere, the F solve's last step before its faces' F (CleanStagFSolve):
//! a face shared with another MeshBlock takes the solution of the one with the lower
//! gid (in our axes, CleanCsSlot), so that both see the same F to the bit (CLF_N + 2),
//! and the cells' slots of cln_v take it

void RadiationM1::CleanStagFCsShared() {
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
  CleanFor("cln_fcs4", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int q = 0; q < 3; ++q) {
      v(m,9+2*q,k,j,i) = x.fc[q](m,CLF_N+2,k,j,i);
      v(m,9+2*q+1,k,j,i) = x.fc[q](m,CLF_N+2,k+(q == 2),j+(q == 1),i+(q == 0));
    }
  });
  // (one MeshBlock along x1: not the x1 family's F, CleanStagFSolve)
  const int f0 = (pmy_pack->pmesh->nmb_rootx1 == 1 && !pmy_pack->pmesh->multilevel) ? 2 : 0;
  CleanGhostFill(cln_v, 9 + f0, 6 - f0, true, 0, 1, CleanFMask(f0));   // (as the F solve)
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagCsTie
//! \brief cubed sphere, clean_stag_ftie: the nc face slots c0, c1, c2 of a face shared
//! with another MeshBlock, in a panel or across a seam, take the values of the lower gid
//! (as its F, CleanStagFCsShared), in our axes (CleanCsSlot), so that both sides use one
//! value to the bit.  Each cell puts those of its face of family q on its MeshBlock's
//! boundary (else its high face) in cln_v slot 3 n + q (the F solve's cell terms, free
//! after it), which the ghosts take.

void RadiationM1::CleanStagCsTie(const int nc, const int c0, const int c1, const int c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto v = cln_v;
  CleanFor("cln_cstie", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int c[3] = {c0, c1, c2};
    for (int q = 0; q < 3; ++q) {
      const int u = (q == 0) ? (i == ie) : ((q == 1) ? (j == je) : (k == ke));
      for (int n = 0; n < nc; ++n) {
        v(m,3*n+q,k,j,i) = x.fc[q](m,c[n],k+u*(q == 2),j+u*(q == 1),i+u*(q == 0));
      }
    }
  });
  CleanGhostFill(cln_v, 0, 3*nc, true, 0, 1);
  for (int dir = 0; dir < 3; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    CleanFor("cln_cstie2", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool lin = (dir == 0) ? (i > is) : ((dir == 1) ? (j > js) : (k > ks));
      const bool rin = (dir == 0) ? (i <= ie) : ((dir == 1) ? (j <= je) : (k <= ke));
      if ((lin && rin) || CleanBoundaryFace(x, dir, m, k, j, i)) return;
      const int kg = lin ? k : k - dk, jg = lin ? j : j - dj, ig = lin ? i : i - di;
      const int ko = lin ? k - dk : k, jo = lin ? j - dj : j, io = lin ? i - di : i;
      if (!(x.g(m,CLG_F+7,kg,jg,ig) < x.g(m,CLG_F+7,ko,jo,io))) return;
      const int c[3] = {c0, c1, c2};
      Real s;
      const int q = CleanCsSlot(x, dir, m, k, j, i, kg, jg, ig, false, s)/2;
      for (int n = 0; n < nc; ++n) fc(m,c[n],k,j,i) = s*v(m,3*n+q,kg,jg,ig);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagCsLateral
//! \brief staggered on the cubed sphere: for E in work vector e, the lateral stress sums
//! of every cell and face family q, Y_q = 1/4 sum over its two faces of A n.(P' - P),
//! P = D E (Cartesian, n outward, P' of the cell across; one sided, 1/2 of the other
//! face, where a side's ghosts do not carry the field, CleanTransSide).  The momentum
//! F_n of a face changes by n.div P of its two half cells: A n.(P_R - P_L) and the Y of
//! the families
//! other than the face's (CleanFaceStress), exact for a uniform P.  The owner computes Y
//! and the ghosts take it, so that a seam face sees the neighbour's true geometry.

void RadiationM1::CleanStagCsLateral(const int e, const bool st, const int part) {
  if (!pmy_pack->pmesh->use_cubed_sphere || cln_cell) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  CleanArgs x = CleanMakeArgs();
  const int ks = indcs.ks, js = indcs.js, is = indcs.is;
  auto cy = cln_csy;
  // (st: from the coefficients of the pass, cln_csyb, CleanPreconditionerSetup)
  if (st) {
    par_for("cln_cslst", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
            indcs.js, indcs.je, indcs.is, indcs.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int q = 0; q < 3; ++q) {
        const int dk = (q == 2), dj = (q == 1), di = (q == 0);
        // (a seam side's slot holds the cell 2 away in panel, cln_csyb)
        const int ol = CleanTransSide(x, q, m, k, j, i, 0) ? -1 : 2;
        const int oh = CleanTransSide(x, q, m, k, j, i, 1) ? 1 : -2;
        for (int b = 0; b < 3; ++b) {
          x.w(m,CLW_L+3*q+b,k,j,i) = cy(m,18+3*q+b,k-ks,j-js,i-is)*x.w(m,e,k,j,i) +
              ((cy(m,6*q+b,k-ks,j-js,i-is) != 0.0) ?   // (none past a boundary)
               cy(m,6*q+b,k-ks,j-js,i-is)*x.w(m,e,k+ol*dk,j+ol*dj,i+ol*di) : 0.0) +
              ((cy(m,6*q+3+b,k-ks,j-js,i-is) != 0.0) ?
               cy(m,6*q+3+b,k-ks,j-js,i-is)*x.w(m,e,k+oh*dk,j+oh*dj,i+oh*di) : 0.0);
        }
      }
    });
    // (read by the faces' own cells: depth 1, the two families other than the face's,
    // CleanStagCsLat; part 1: CleanApply finishes it)
    CleanGhostFill(cln_w, CLW_L, 9, part, 1, {0770, 0707, 0077, 0});
    return;
  }
  auto f = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real p0[6];   // (the cell's own P, once)
    for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
    for (int q = 0; q < 3; ++q) {
      const int dk = (q == 2), dj = (q == 1), di = (q == 0);
      Real y[3] = {0.0, 0.0, 0.0};
      // (a side whose ghosts do not carry the field drops out, the other side's weight
      // doubled: one sided, CleanTransSide)
      const int vl = CleanTransSide(x, q, m, k, j, i, 0);
      const int vh = CleanTransSide(x, q, m, k, j, i, 1);
      // (a seam side keeps its face, the ghost's P extrapolated in panel, clamped, as
      // production's M1CsTanDer, CleanCsSeamCoef; an open boundary drops its side, the
      // other's weight doubled)
      const bool o2 = (vl + vh == 1) && x.mbc(m,2*q+vl) == BoundaryFlag::panel;
      const Real wq = (vl + vh == 2 || o2) ? 0.25 : 0.5;
      for (int s = -1; s <= 1; s += 2) {
        const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
        const bool ex = o2 && !CleanTransSide(x, q, m, k, j, i, s > 0);
        if (!ex && !CleanTransSide(x, q, m, k, j, i, s > 0)) continue;
        const Real a = wq*s*x.area[q](m,kf,jf,iff);
        if (!ex) {   // (each component's P(s) - P(0) and the normal once)
          Real dp[6], nf[3];
          for (int n = 0; n < 6; ++n) {
            dp[n] = x.cs(m,CLN_D11+n,k+s*dk,j+s*dj,i+s*di)*x.w(m,e,k+s*dk,j+s*dj,i+s*di) -
                    p0[n];
          }
          for (int c = 0; c < 3; ++c) nf[c] = x.g.At<CLG_N>(m,3*q+c,kf,jf,iff);
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
      for (int b = 0; b < 3; ++b) x.w(m,CLW_L+3*q+b,k,j,i) = y[b];
    }
  };
  // (the cells it sends first, its exchange in flight behind the inner ones)
  const bool ov = (cln_hmax >= 9);
  const std::function<void()> post = [&]() {
    CleanGhostFill(cln_w, CLW_L, 9, ov ? 1 : 0, 1, {0770, 0707, 0077, 0});   // (as above)
  };
  const std::function<void()> fin = [&]() {
    if (ov) CleanGhostFill(cln_w, CLW_L, 9, 2, 1);
  };
  CleanOverlapFor("cln_cslat", pmy_pack, ov, f, fin, -1, post);
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

void RadiationM1::CleanStagDrag(const bool first, const bool vel, const bool lat) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // (cubed sphere: the lateral stress of the iterate; lat false: the caller's)
  if (lat) CleanStagCsLateral(CLW_EK);
  CleanArgs x = CleanMakeArgs(), xh = x;
  xh.hi = true;
  xh.kq = -1;
  xh.drag = x.drag && !first;   // (the flux the iterate's velocity came from)
  auto uh = FluidRef::Get(pmy_pack).u0;
  const bool off = cln_vzero || (first && force_ref == M1_FREF_CLEAN && !cln_fref_set);
  const bool lag = !first || (force_ref != M1_FREF_NONE);
  const bool dr = cln_vimp;
  const int nf = x.nf;   // (one launch for the faces of all directions)
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  CleanFor("cln_drag", nmb1, ks, ke+(nf > 2), js, je+(nf > 1), is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    for (int dir = 0; dir < nf; ++dir) {
      const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      if (i > ie+di || j > je+dj || k > ke+dk) continue;
      auto &fc = x.fc[dir];
      if (CleanBoundaryFace(x, dir, m, k, j, i)) continue;   // (no F equation)
      const Real s = (vel || lag) ? CleanFaceStress(x, CLW_EK, dir, m, k, j, i, true)
                                  : 0.0;
      // (vel: the F solve's fs0 just took F^H of E_k with these weights, CLF_PT)
      const Real fq = (vel && lag) ? fc(m,CLF_PT,k,j,i) : 0.0;
      if (vel) {
        fc(m,CLF_PT,k,j,i) = (fc(m,CLF_N,k,j,i) - fc(m,CLF_G,k,j,i))/(x.c*x.c) + x.h*s -
                             x.h*fc(m,CLF_N+1,k,j,i);
      }
      if (!dr) continue;
      const Real q = !lag ? fc(m,CLF_G,k,j,i) :
                     (vel ? fq : CleanFaceFlux(xh, CLW_EK, dir, m, k, j, i, true)) +
                     x.h*x.c*x.c*s;
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
//! \fn void RadiationM1::CleanLzBook
//! \brief clean_cs_lz_abs = 2: books the cells' sum s of the L_z the gas did not take (of
//! a write-back, or with start of the stage start's increment) as the gas momentum it
//! belongs to enters the step: a stage-1 part with the weight of the stage-1 state in the
//! stage-2 start (Heun 1/2), the stage-2 start and write-back pending until the solve
//! accepts stage 2 (s = 0; a backward-Euler redo from U^n drops them, Time2Restore); the
//! total and K1's sum in <rad_m1> (restart).
//! 3-D spherical polar: s is the L_z that crossed the boundaries (CleanStagMomentum's dlz),
//! <rad_m1>/clean_sp_lz_bnd

void RadiationM1::CleanLzBook(const Real s, const bool start) {
  if (t2_solve == M1_T2S_STAGE1) {
    cln_lzp = (start ? 0.0 : cln_lzp) + ((t2_tab == 0 || t2_tab == 2) ? 0.5 : t2_w0)*s;
    return;
  }
  if (t2_solve == M1_T2S_STAGE2 && start) {
    cln_lzp += s;
    return;
  }
  cln_lzl += (t2_solve == M1_T2S_STAGE2) ? cln_lzp : cln_lzp + s;
  cln_lzp = 0.0;
  const bool cs = (cln_lza == 2);
  cln_lzpin->SetReal("rad_m1", cs ? "clean_cs_lz_ledger" : "clean_sp_lz_bnd", cln_lzl);
  cln_lzpin->SetReal("rad_m1", cs ? "clean_cs_lz_k1" : "clean_sp_lz_bnd_k1", cln_lzs1);
  if (global_variable::my_rank == 0 && !cln_insteady) {
    std::printf("CLEAN LZ%s cycle %d t2s %d ledger %.16e\n", cs ? "" : "B",
                pmy_pack->pmesh->ncycle, t2_solve, cln_lzl);
  }
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
  xh.kq = -1;   // (clean_stag_kpos: the momentum without the moved flux)
  auto kk = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;
  // h div Phi of the solution (CLW_Q, a slot the F solve leaves alone) before the F solve
  // below moves the energy flux's memory d (clean_stag_fkin): E takes the residual of the
  // operator it was solved with; with clean_stag_split also the cell flux, that of Phi
  // (the faces keep F^H)
  auto u0_ = u0;
  // (clean_stag_ftie: a high face's Phi its owner's, CLW_L + dir of the next cell, so the
  // energy telescopes across MeshBlock faces and the theta-periodic wrap; the solve's
  // own Phi there differs by its copy's geometry, a residual E takes; the cubed sphere:
  // each shared face's CLF_PT the lower gid's, CleanStagCsTie)
  const bool tie = cln_ftie && !x.cub;
  if (cln_ftie && x.cub) CleanStagCsTie(1, CLF_PT);
  if (tie) {
    const int nf = x.nf;
    CLEAN_HIPX(const CleanArgs *px = CleanDevArgs(x, 0);)
    CleanFor("cln_wbphi", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
      for (int d = 0; d < nf; ++d) {   // (CleanStagRealizable's Phi but on physical faces)
        x.w(m,CLW_L+d,k,j,i) = CleanBoundaryFace(x, d, m, k, j, i) ?
            CleanFaceFlux(x, CLW_X, d, m, k, j, i, true) : x.fc[d](m,CLF_PT,k,j,i);
      }
    });
    CleanGhostFill(cln_w, CLW_L, nf);
  }
  // (clean_stag_xdonor: the cell flux, the radiation's momentum, of the central cross flux
  // and face sources, as the gas momentum's stress (CleanMomentumChange): Phi's donor
  // form equals it only at the exact fixed point of q, so total momentum would carry
  // the passes' tolerance)
  CleanArgs xc = x;
  xc.xq = -1;
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);
  const CleanArgs *pxc = CleanDevArgs(xc, 2);
  // (time2_face_fct: the energy flux Phi of every face, which moved E, with its slope)
  const bool pfo = t2_pfst && slopes && t2pf1.extent_int(0) > nmb1;
  auto pf = (t2_solve == M1_T2S_STAGE1) ? t2pf2 : t2pf1;
  // (xdonor: each face's central cross flux once, for the two cells that read it, in
  // cln_xf's scratch (dead here: its users write it before they read it); A100 R136
  // cln_wbdiv 1.6 -> 14.5 ms per launch with it evaluated twice per face)
  auto phq = cln_phic;   // (the cell flux of the faces' Phi, output m1_phi)
  DvceArray5D<Real> xcv;
  if (xh.hi && x.xq >= 0) {
    xcv = ScratchView<DvceArray5D<Real>>("cln_xf", nmb1 + 1, 3, cln_w.s.extent_int(2),
                                         cln_w.s.extent_int(3), cln_w.s.extent_int(4));
    for (int n = 0; n < x.nf; ++n) {
      CleanFor("cln_wbxc", nmb1, ks, ke+(n == 2), js, je+(n == 1), is, ie+(n == 0),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        xcv(m,n,k,j,i) = CleanFaceFlux(*pxc, CLW_X, n, m, k, j, i, true);
      });
    }
  }
  CleanFor("cln_wbdiv", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (device memory: no copy)
    // (CleanFluxDivergence(Tie) with the faces' Phi of CleanStagRealizable, CLF_PT, but
    // on physical faces; with the tie the low faces' CLW_L and the high faces' owner's)
    Real dv = 0.0, pc[6];   // (the cubed sphere: each face's, CleanStagCellVector)
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
      // (clean_stag_kpos: the cell flux without the moved increment flux)
      if (xh.hi && x.xq >= 0) {
        phq(m,n,k,j,i) = 0.5*(xcv(m,n,k,j,i) + xcv(m,n,k+dk,j+dj,i+di) -
                              CleanKposFlux(x, CLW_X, n, m, k, j, i) -
                              CleanKposFlux(x, CLW_X, n, m, k+dk, j+dj, i+di));
      } else if (xh.hi) {
        phq(m,n,k,j,i) = 0.5*(fl + fu - CleanKposFlux(x, CLW_X, n, m, k, j, i) -
                              CleanKposFlux(x, CLW_X, n, m, k+dk, j+dj, i+di));
      }
      if (xh.hi && x.cub) {
        pc[2*n] = ((x.xq >= 0) ? xcv(m,n,k,j,i) : fl) -
                  CleanKposFlux(x, CLW_X, n, m, k, j, i);
        pc[2*n+1] = ((x.xq >= 0) ? xcv(m,n,k+dk,j+dj,i+di) : fu) -
                    CleanKposFlux(x, CLW_X, n, m, k+dk, j+dj, i+di);
      }
      // (clean_stag_kpos: the stored Phi without the moved increment flux)
      if (pfo) pf(m,n,k,j,i) = fl - CleanKposFlux(x, CLW_X, n, m, k, j, i);
      // (the high face of the MeshBlock: the Phi that moved E, the owner's with the tie,
      // so time2_face_fct limits one flux on both sides of a seam or the theta wrap)
      if (pfo && (n == 0 ? i == x.ie : (n == 1 ? j == x.je : k == x.ke))) {
        pf(m,n,k+dk,j+dj,i+di) = fh - CleanKposFlux(x, CLW_X, n, m, k+dk, j+dj, i+di);
      }
    }
    if (xh.hi && x.cub) {   // (the cubed sphere: the faces' Phi's vector, frame of y)
      Real fv[3];
      CleanStagCellVector(x, m, k, j, i, pc, fv);
      CleanStagHydro(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) phq(m,d,k,j,i) = fv[d];
    }
    x.w(m,CLW_Q,k,j,i) = x.h*dv;
  });
  if (pfo) {
    if (t2_solve == M1_T2S_STAGE1) {t2pf2_ok = true;} else {t2pf1_ok = true;}
  }
  // (the history's boundary energy: + h A Phi of the physical x1 faces, CleanBndAcc)
  if (cln_lacc.extent_int(0) > nmb1) {
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;   // (the x1 face area of the clean fluxes)
    const bool stg = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
    par_for("cln_wbbnd", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      // (clean_stag_kpos: the boundary outflow of the increment, pd_b X of the cell)
      const int ic = s ? ie : is;
      const Real kb = (x.kq >= 0) ? (s ? 1.0 : -1.0)*x.w(m,CLW_KF,k,j,ic)*
                      fmax(x.w(m,CLW_X,k,j,ic), 0.0)*x.vol(m,k,j,ic) : 0.0;
      la(m,s,k,j) = !CleanBoundaryFace(x, 0, m, k, j, i) ? 0.0 :
                    la(m,stg ? 2+s : s,k,j) + h*ar(m,k,j,i)*CleanFaceFlux(x, CLW_X, 0, m,
                                                                     k, j, i, true) + kb;
    });
  }
  // (the face F of the F solve and its dissipation; the lateral sums of X: above)
  CleanStagFSolve(CLW_X, true, false);
  xh.fs = true;
  xh.tf = cln_ftie && !x.cub;
  xh.fa = cln_ftie && x.cub;   // (the cubed sphere's tied T: CLF_PT, CleanStagFaceTie)

  // (clean_cs_lz_tel: the stress torque CleanCsTorque's xs, a defect of the radiation's
  // momentum equation, is zero-sum noise per cell (sum |xs| ~1e5 its net); the gas takes
  // its net X = sum xs only, as the solid-body rotation dOmega = X/I of least kinetic
  // energy, I = sum rho V R^2, so no cell takes the per-cell defect (LZTEL3: a cell's
  // share spun thin gas up, and in the faces it made F unrealisable).  tqv: the gas's
  // L_z gl and the ledger's; lxs, lxi: xs and rho V R^2, lzs their sums X and I)
  const bool lzt = cln_lzt && feedback;
  DvceArray5D<Real> tqv;
  Real lzs[2] = {0.0, 0.0};
  if (lzt) {
    tqv = ScratchView<DvceArray5D<Real>>("cln_tqv", nmb1 + 1, 3, u0.extent_int(2),
                                         u0.extent_int(3), u0.extent_int(4));
    auto lxs = ScratchView<DvceArray4D<Real>>("cln_lxs", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    auto lxi = ScratchView<DvceArray4D<Real>>("cln_lxi", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    px = CleanDevArgs(x, 0);   // (xh changed: the F solve's F)
    pxh = CleanDevArgs(xh, 1);
    auto ug = FluidRef::Get(pmy_pack).u0;
    CleanFor("cln_wbtq", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh;)
      Real rl, xs, q[3][3];
      const Real gl = CleanCsTorque(xh, CLW_X, m, k, j, i, rl, xs);
      const Real rc = CleanCsFrame(x, m, k, j, i, q);
      lxs(m,k,j,i) = xs;
      lxi(m,k,j,i) = ug(m,IDN,k,j,i)*x.vol(m,k,j,i)*rc*rc;
      tqv(m,1,k,j,i) = gl - xs;   // (the absorption torque)
      tqv(m,2,k,j,i) = -(gl + rl);
    });
    lzs[0] = CleanSumCells(lxs);
    lzs[1] = CleanSumCells(lxi);
  }
  const Real lzw_ = (lzs[1] > 0.0) ? lzs[0]/lzs[1] : 0.0;   // (dOmega)

  // face fluxes of the solution
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = M1_T2_F1 + dir;
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
  const bool wsp = (fref == M1_FREF_WB_ARAD || fref == M1_FREF_WB_MD);
  const bool fcap = (fref == M1_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
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
  const Real pgf = impl_pos_gas_frac, efl = e_floor, pgk = impl_pos_flr_gk;
  auto pc_ = pos_cnt_d;
  const bool pflr = impl_pos_floor, pfc = pflr && impl_pos_flr_cons;
  if (pfc && pfc_def.extent_int(0) != nmb1 + 1) {
    auto &icp = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(pfc_def, nmb1 + 1, icp.nx3 + 2*icp.ng, icp.nx2 + 2*icp.ng,
                    icp.nx1 + 2*icp.ng);
  }
  if (pfc) Kokkos::deep_copy(pfc_def, 0.0);
  auto pfd = pfc_def;
  const bool pgc = pfc && pgas && impl_pos_gas_cons;   // (implicit_pos_gas_cons)
  if (pgc && pgc_def.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(pgc_def, nmb1 + 1, pfc_def.extent_int(1), pfc_def.extent_int(2),
                    pfc_def.extent_int(3));
  }
  if (pgc) Kokkos::deep_copy(pgc_def, 0.0);
  auto pgd = pgc_def;
  const bool ky1s = cln_kpos && t2_solve == M1_T2S_STAGE1;
  if (ky1s && cln_ky1.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(cln_ky1, nmb1 + 1, u0.extent_int(2), u0.extent_int(3), u0.extent_int(4));
  }
  auto ky1v = cln_ky1;
  const bool xwb = (cln_xte > 0) && feedback && !cln_vzero;
  const bool krl = cln_krel_on && feedback;   // (clean_stag_kreal: CLW_KR)
  const bool kpb = cln_kpb && feedback;
  const bool cls = xwb || kpb;   // (the closure in CLW_PW, printed below)
  px = CleanDevArgs(x, 0);   // (the F solve took the slots; xh changed)
  pxh = CleanDevArgs(xh, 1);
  const bool lzb = cln_lzbook && feedback;
  DvceArray4D<Real> lzw;
  if (lzb) lzw = ScratchView<DvceArray4D<Real>>("cln_lzw", nmb1 + 1, u0.extent_int(2),
                                                 u0.extent_int(3), u0.extent_int(4));
  CleanFor("cln_wbcell", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    CLEAN_HIPX(const CleanArgs &x = *px, &xh = *pxh;)   // (in device memory)
    Real e = x.w(m,CLW_X,k,j,i);
    // (clean_stag_kpos: h div of the moved increment flux, the diagonal's boundary part
    // pd_b X (boundary energy) and its gas part pd_g X)
    Real hdq = 0.0, pbx = 0.0, pgx = 0.0;
    if (x.kq >= 0) {
      for (int n = 0; n < x.nf; ++n) {
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        const Real ql = CleanKposFlux(x, CLW_X, n, m, k, j, i);
        const Real qh = CleanKposFlux(x, CLW_X, n, m, k+dk, j+dj, i+di);
        hdq += h*CleanDivF(x, n, m, k, j, i, ql, qh);
      }
      pbx = x.w(m,CLW_KF,k,j,i)*fmax(e, 0.0);
      pgx = (x.w(m,CLW_KM,k,j,i) - x.w(m,CLW_KF,k,j,i))*fmax(e, 0.0);
    }
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
    Real dm[3] = {0.0, 0.0, 0.0}, dh[3] = {0.0, 0.0, 0.0}, dkr[3] = {0.0, 0.0, 0.0};
    Real den = 0.0, wsl = 0.0, wkr = 0.0;
    const Real d0 = x.w(m,CLW_JV,k,j,i);   // (the stage start's eint deficit, <= 0)
    if (feedback) {
      Real mg[3], pg[3], rf[3], dlz = 0.0;
      CleanStagMomentum(xh, CLW_X, m, k, j, i, dm, mg, &dlz);
      if (lzt) {   // (clean_cs_lz_tel: the gas phi momentum of CleanCsTorque, cln_wbtq)
        Real q[3][3];
        const Real rc = CleanCsFrame(x, m, k, j, i, q);
        const Real dp = tqv(m,1,k,j,i)/(x.vol(m,k,j,i)*rc) + lzw_*uh(m,IDN,k,j,i)*rc -
                        (q[2][0]*dm[0] + q[2][1]*dm[1] + q[2][2]*dm[2]);
        // (force_reference = clean: the reference withholds this phi row too, cln_fref 3,
        // its torque booked as that of cln_fref 0..2)
        const Real dw = (fref != M1_FREF_CLEAN) ? 0.0 : (fcap ? dp : h*fr(m,3,k,j,i));
        if (fcap) fr(m,3,k,j,i) = dp/h;
        for (int n = 0; n < 3; ++n) dm[n] += (dp - dw)*q[2][n];
        dlz = tqv(m,2,k,j,i) + x.vol(m,k,j,i)*rc*dw;
      }
      CleanStagGasMomentum(x, uh, m, k, j, i, pg);
      CleanStagRefForce(x, fref, fr, aref, uh(m,IDN,k,j,i), m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {   // (force_reference)
        if (fref != M1_FREF_NONE) dm[n] -= h*rf[n];
      }
      if (lzb && fref == M1_FREF_CLEAN) {   // (the reference torque withheld: booked)
        Real q[3][3];
        const Real r = x.cub ? CleanCsFrame(x, m, k, j, i, q) : x.x1v(m,i)*x.tg(m,2,j);
        dlz += h*x.vol(m,k,j,i)*r*(x.cub ? q[2][0]*rf[0] + q[2][1]*rf[1] + q[2][2]*rf[2]
                                         : rf[2]);
      }
      if (fm) fmd(m,k,j,i) = CleanStagRadial(x, m, k, j, i, dm);
      if (x.vzero) dm[0] = dm[1] = dm[2] = dlz = 0.0;   // (frozen gas momentum)
      if (lzb) lzw(m,k,j,i) = dlz;   // (clean_cs_lz_abs = 2: CleanLzBook)
      // (clean_stag_kreal: the stage's deferred kick, CleanStagKReal; its own kinetic
      // energy wkr leaves the slopes, which are the stage's f)
      if (krl && !x.vzero) {
        for (int n = 0; n < 3; ++n) {
          const Real kr = x.w(m,CLW_KR+n,k,j,i);
          wkr += (pg[n] + 0.5*kr)*kr/uh(m,IDN,k,j,i);
          dm[n] += kr;
        }
        Real kv[3] = {x.w(m,CLW_KR,k,j,i), x.w(m,CLW_KR+1,k,j,i), x.w(m,CLW_KR+2,k,j,i)};
        CleanGasChange(x, m, k, j, i, kv, dkr);
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
        Real dy[3] = {dm[0], dm[1], dm[2]};
        CleanStagHydro(x, m, k, j, i, dy);
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
      // (kpos: the Cartesian face F carries no moved flux, the boundary's pd_b X leaves)
      den = -(e - x.w(m,CLW_EO,k,j,i)) - ((xh.hi || x.sph || x.cub) ?
            x.w(m,CLW_Q,k,j,i) : h*div + hdq) + (eso ? h*es(m,k,j,i) : 0.0) - pbx;
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
        // energy change did not receive.  The floor only here: inside the Newton passes
        // a root below it switched the gas row whole into E and back, a 2-cycle of the
        // passes; a root below it raises the gas to e(rho, tfloor) with E's energy,
        // counted as implicit_pos_gas)
        Real eg = x.w(m,CLW_EP,k,j,i) - ex + pgx;   // (kpos: the gas's pd_g X)
        if (gtf) {
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(uh(m,IDN,k,j,i), fmax(t, tgf), ei, p, chi_rho, chi_t, cv);
          eg = ei - x.w(m,CLW_EI,k,j,i);
          if (t < tgf) {
            Real e0;
            eos.ThermoAt(uh(m,IDN,k,j,i), t, e0, p, chi_rho, chi_t, cv);
            const Real vol = CleanVol(x, m, k, j, i);
            Kokkos::atomic_add(&pc_(M1_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(M1_POS_GAS_DE), (ei - e0)*vol);
          }
        }
        // (clean_xt_exact: the work the row took and the gas's kinetic energy did not
        // receive to E whole, eint the heating only)
        const Real dw = xwb ? x.w(m,CLW_XW,k,j,i) - wk - wp : 0.0;
        const Real rd = (den - eg - wk - wp - dw)/(gtf ? 1.0 :
                        1.0 + x.h*x.c*x.op(m,M1_OP_E,k,j,i)*x.w(m,CLW_RY,k,j,i)) + dw;
        if (cls) x.w(m,CLW_PW,k,j,i) = fabs(den - eg - wk - wp)/fmax(fabs(e) +
            fabs(x.w(m,CLW_EO,k,j,i)) + fabs(x.w(m,CLW_EI,k,j,i)), 1.0e-300);
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
          // (from the gas's true eint: the stage start's deficit below 0 too, which the
          // exchange took as 0 and the hydro's floor would have created)
          const Real ei = x.w(m,CLW_EI,k,j,i), emin = fmax(ef, pgf*fmax(ei, 0.0));
          const Real need = emin - ei - d0 - (den - wk);
          const Real take = fmin(need, e - efl);
          const Real vol = CleanVol(x, m, k, j, i);
          const bool tk = take > 1.0e-12*emin;   // (beyond eint's round-off at e_min)
          if (tk) {
            den += take;
            e -= take;
            Kokkos::atomic_add(&pc_(M1_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(M1_POS_GAS_DE), take*vol);
          }
          // (implicit_pos_gas_cons: what E could not cover, from the column, below)
          const Real un = need - (tk ? take : 0.0);
          if (pgc && un > 1.0e-12*emin) pgd(m,k,j,i) = un*vol;
        }
      }
      uh(m,IEN,k,j,i) += den;
      x.w(m,CLW_JV,k,j,i) = den - wk;   // (its change of the gas eint, cln_neint)
    }
    // (implicit_pos_floor, the old production's write-back: an E at or below e_floor is
    // raised to it with the energy of the cell's gas above e_min; what the gas cannot
    // cover is taken from the radiation of the cell's column (implicit_pos_floor_cons,
    // below), else counted as created; counted)
    if (pflr && e <= efl) {
      const Real vol = CleanVol(x, m, k, j, i);
      const bool gas = feedback && !x.vzero;
      Real take = 0.0;
      if (gas) {
        Real ef = 0.0, pr, chi_rho, chi_t, cv;
        if (eos.tfloor > 0.0) {
          eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
        }
        const Real ei = x.w(m,CLW_EI,k,j,i), en = ei + d0 + x.w(m,CLW_JV,k,j,i);
        // (implicit_pos_floor_gas_keep: the gas keeps pgk of its eint)
        const Real emin = fmax(fmax(ef, pgf*fmax(ei, 0.0)), pgk*fmax(en, 0.0));
        take = fmin(efl - e, fmax(en - emin, 0.0));
        den -= take;
        uh(m,IEN,k,j,i) -= take;
        x.w(m,CLW_JV,k,j,i) -= take;
        Kokkos::atomic_add(&pc_(M1_POS_FLR_DE), take*vol);
      }
      if (pfc) {
        pfd(m,k,j,i) = (efl - e - take)*vol;
      } else {
        Kokkos::atomic_add(&pc_(M1_POS_FLR_UN), (efl - e - take)*vol);
      }
      e = efl;
      Kokkos::atomic_add(&pc_(M1_POS_FLR), 1.0);
    }
    if (slopes) {
      // (clean_stag_kpos: the slopes f(U) without the moved increment, hdq + pd X, and
      // the gas's pd_g X: K_E + K_EN = -div of the stored Phi; clean_stag_kreal: without
      // the deferred kick and its kinetic energy wkr)
      kk(m,M1_T2_E,k,j,i) = (e + wsl - x.w(m,CLW_EO,k,j,i) + hdq + pbx + pgx + wkr)/h;
      for (int n = 0; n < 3; ++n) kk(m,M1_T2_M1+n,k,j,i) = (dh[n] - dkr[n])/h;
      kk(m,M1_T2_EN,k,j,i) = (den - (feedback ? pgx : 0.0) - wkr)/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    if (ky1s) ky1v(m,k,j,i) = e;   // (clean_stag_kpos: K2's reference, CleanStagKPos)
    {   // (the cell flux: the face momenta's cell vector, the radiation momentum on every
        //  mesh, in the frame of y on the cubed sphere; Phi's in cln_phic, above)
      const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), multi_d ? f2(m,k,j,i) : 0.0,
                         multi_d ? f2(m,k,j+1,i) : 0.0, three_d ? f3(m,k,j,i) : 0.0,
                         three_d ? f3(m,k+1,j,i) : 0.0};
      Real fv[3];
      CleanStagCellVector(x, m, k, j, i, c, fv);
      CleanStagHydro(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) {
        u0_(m,M1_F1+d,k,j,i) = fv[d];
        if (!xh.hi) phq(m,d,k,j,i) = fv[d];   // (no split: Phi is the face momentum)
      }
    }
  });
  if (lzb) {   // (the L_z ledger: the stage's L_z as its slope, then booked)
    const Real s = CleanSumCells(lzw);
    if (slopes) (t2_solve == M1_T2S_STAGE1 ? cln_lzs2 : cln_lzs1) = s/h;
    CleanLzBook(s, t2_solve == M1_T2S_STAGE2);
  }
  if (cls) {   // (xt_exact / kpos: the stage's closure |dE + de + dK + h div Phi|
               //  over E + E_old + eint, max over the cells)
    auto w = cln_w;
    const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
    Real cm = 0.0;
    Kokkos::parallel_reduce("cln_wbclose", CleanRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &c) {
      const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
      c = fmax(c, w(m,CLW_PW,k,j,is + idx%n1));
    }, Kokkos::Max<Real>(cm));
    cm = CleanReduce(cm, true);
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN %s cycle %d t2s %d write-back closure %.3e\n", xwb ? "XT" : "KPOS",
                  pmy_pack->pmesh->ncycle, t2_solve, cm);
    }
  }
  cln_kpb = false;
  // (implicit_pos_floor_cons, the old production's column pass: the uncovered raise of a
  // cell from the radiation above e_floor of the nearest cells of its radial column, i -+
  // d, the pair split by their energy above the floor; their stage slopes with it.
  // clean_stag_colflux: the stored Phi carries each move as x1 face fluxes from q to i,
  // so the slopes stay -div Phi + their cell pairs and time2_face_fct limits the moves)
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
          const Real aa = (qa >= is) ? fmax(u0_(m,M1_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,M1_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,M1_E,k,j,q) = fmax(u0_(m,M1_E,k,j,q) - de, efl);
            if (slopes) kk(m,M1_T2_E,k,j,q) -= de/h;
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
      if (mv > 0.0) Kokkos::atomic_add(&pc_(M1_POS_FLR_MV), mv);
      if (un > 0.0) Kokkos::atomic_add(&pc_(M1_POS_FLR_UN), un);
      if (!pgc) return;
      // (implicit_pos_gas_cons: the gas raise E could not cover, the same way, into the
      // cell's gas; what the column cannot give stays below e_min: the hydro EOS floor
      // ledger books it)
      Real gm = 0.0;
      for (int i = is; i <= ie; ++i) {
        Real def = pgd(m,k,j,i);
        for (int d = 1; d <= ie - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(u0_(m,M1_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,M1_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,M1_E,k,j,q) = fmax(u0_(m,M1_E,k,j,q) - de, efl);
            if (slopes) kk(m,M1_T2_E,k,j,q) -= de/h;
            for (int f = (sd ? i + 1 : qa + 1); pfk && f <= (sd ? qb : i); ++f) {
              pf(m,0,k,j,f) += (sd ? -1.0 : 1.0)*de*vol(q)/
                               (h*CleanDivW(x, 0, m, k, j, f, 0)*vol(f));
            }
          }
          const Real dg = tk/vol(i);
          uh(m,IEN,k,j,i) += dg;
          x.w(m,CLW_JV,k,j,i) += dg;
          if (slopes) kk(m,M1_T2_EN,k,j,i) += dg/h;
          def -= tk;
          gm += tk;
        }
      }
      if (gm > 0.0) Kokkos::atomic_add(&pc_(M1_POS_GAS_DE), gm);
    });
  }
  // cells the write-back left with gas eint < 0, or below the EOS floor, by lowering it
  // (the solve is then not accepted, as in the cell layout; an eint the hydro left there
  // is not the solve's)
  cln_neint = 0;
  if (!feedback) return;
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
    // (and below the eint of the EOS temperature floor: the hydro's floor would add
    // energy the radiation's budget does not see)
    Real ef = 0.0;
    if (eos.tfloor > 0.0) {
      Real pr, chi_rho, chi_t, cv;
      eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
    }
    const Real tz = 1.0e-12*(fabs(uh(m,IEN,k,j,i)) + ke + ef);   // (eint's round-off)
    // (and any non-finite E, F or gas state, also of a cell at the floor: rejected and
    // retried as the others, never accepted)
    const bool fin = Kokkos::isfinite(u0_(m,M1_E,k,j,i) + u0_(m,M1_F1,k,j,i) +
                                      u0_(m,M1_F1+1,k,j,i) + u0_(m,M1_F1+2,k,j,i) +
                                      uh(m,IDN,k,j,i) + uh(m,IEN,k,j,i) + ke);
    const bool low = !(CleanGasEint(x, uh, ke, m, k, j, i) > ef - tz);
    if ((low && !(x.w(m,CLW_JV,k,j,i) >= 0.0) && x.w(m,CLW_EI,k,j,i) > ef + tz) ||
        u0_(m,M1_E,k,j,i) < 0.0 || !fin) n += 1.0;
    // (without implicit_pos_gas: an eint the hydro's floor will raise, counted)
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
#endif  // M1_TU_PART(2)

}  // namespace radm1
