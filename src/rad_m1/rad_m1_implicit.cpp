//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_implicit.cpp
//! \brief IMPLICIT transport along x1 columns, <rad_m1>/transport = implicit_x1
//! (milestone 3a of docs/dev/rad_m1_implicit_design.md).
//!
//! ONE backward-Euler solve per hydro step, at the TRUE speed of light, no sub-cycling
//! and no PD-ARS.  The normal flux lives on the x1 FACES and is eliminated there, which
//! turns the coupled (E, F) system into one tridiagonal M-matrix system for E per column:
//!
//!   F0'_f = theta_f [ F0^n_f - chat c dt (w_i E'_i - w_{i-1} E'_{i-1})/dx
//!                     - chat dt v_f g0_f ],    theta_f = 1/(1 + chat dt (rho k_t)_f)
//!   E'_i + (dt/dx)(chat/c)[ (F0' + A')_{i+1/2} - (F0' + A')_{i-1/2} ]
//!        = E^n_i + dt chat (rho kappa_P a T'^4 - rho kappa_E E0'_i)
//!
//! with w = P_11/E from the LAGGED closure (= chi in 1-D), A = a E the enthalpy flux
//! (a = v1 (1 + w)) upwinded with the face velocity, and the emission term linearised in
//! T and eliminated into the diagonal (design sect. 2).  Both the diffusion and the
//! upwind-advection off-diagonals are non-positive and the source adds
//! dt chat rho kappa_E rho c_v/B >= 0 to the diagonal, so the matrix is an M-matrix and
//! E' > 0 at any dt.
//!
//! OUTER LOOP: Picard on (w, theta, a, de0, g0, the upwind directions, optionally the
//! opacities) plus the safeguarded scalar root find for T' at fixed E'.  At convergence
//! the linearised emission correction vanishes identically, so the converged state
//! satisfies the NONLINEAR backward-Euler equations (see the note above
//! M1ImplTemperature).
//!
//! RESTRICTIONS of 3a, all fatal:
//!   * exactly ONE MeshBlock along x1 (each column is solved by a plain Thomas sweep
//!     inside one block; the partitioned line solve of two_stream_column_partition.hpp
//!     is NOT used);
//!   * nx2 = nx3 = 1 unless <rad_m1>/implicit_allow_multid = true, in which case the
//!     columns are INDEPENDENT: there is no transport along x2/x3 in this mode and
//!     F_2 = F_3 = 0 always;
//!   * no SMR/AMR.

#include <float.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn M1EnthIdx
//! \brief implicit_enthalpy: the cell index a face stencil may read along one direction,
//! for the raw (unwrapped) index ii.  With a periodic wrap inside the block (cyc) the
//! index is wrapped; otherwise it must lie in [lo - hlo, hi + hhi], hlo/hhi being the
//! number of ghost layers filled on that side (0 at a physical boundary).

KOKKOS_INLINE_FUNCTION
int M1EnthIdx(const int ii, const int lo, const int hi, const bool cyc, const int hlo,
              const int hhi, bool &ok) {
  if (cyc) {
    const int n = hi - lo + 1;
    int r = ii;
    if (r < lo) {r += n;}
    if (r > hi) {r -= n;}
    ok = true;
    return r;
  }
  ok = (ii >= lo - hlo) && (ii <= hi + hhi);
  return ok ? ii : lo;
}

//----------------------------------------------------------------------------------------
//! \fn M1EnthCorr
//! \brief implicit_enthalpy: (high-order face enthalpy flux) - (donor-cell face flux the
//! matrix carries), both at the lagged iterate.  el2, el, er, er2 are E and al2, al, ar,
//! ar2 the advective coefficients a at the cells L-1, L, R, R+1 of the face, and ok says
//! whether L-1 AND R+1 exist; vf is the face velocity whose sign picks the donor cell of
//! the matrix.  Requiring BOTH outer cells for the plm form makes the choice the same for
//! the two MeshBlocks that share a face (each of them misses one of the two when the halo
//! is too thin), so the face flux stays single valued.
//!   central: a_f = (a_L + a_R)/2, E_f = (E_L + E_R)/2.
//!   plm:     a_f = the mean of the two van Leer plm face values of a (it reduces to the
//!            central mean at an extremum of a and is a 4-point interpolation where a is
//!            smooth and monotone), E_f = the van Leer plm value of E from the side
//!            upwind of a_f (monotone, E_f <= 2 E_donor).  Central where !ok.

KOKKOS_INLINE_FUNCTION
Real M1EnthCorr(const int mode, const Real el2, const Real el, const Real er,
                const Real er2, const bool ok, const Real al2, const Real al,
                const Real ar, const Real ar2, const Real vf) {
  const Real alow = (vf > 0.0) ? (al*el) : (ar*er);
  Real af = 0.5*(al + ar);
  Real ef = 0.5*(el + er);
  if (mode == M1_IENTH_PLM && ok) {
    Real dum, afl, afr;
    PLM(al2, al, ar, afl, dum);
    PLM(al, ar, ar2, dum, afr);
    af = 0.5*(afl + afr);
    if (af > 0.0) {
      PLM(el2, el, er, ef, dum);
    } else if (af < 0.0) {
      PLM(el, er, er2, dum, ef);
    }
  }
  return af*ef - alow;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn ImplBCFromString
//! \brief one x1 boundary type of the implicit solve, from the input string, or from the
//! mesh boundary flag when the input does not name one.

int ImplBCFromString(const std::string &s, const BoundaryFlag mbc) {
  if (s.compare("auto") == 0) {
    if (mbc == BoundaryFlag::periodic) return M1_IBC_PERIODIC;
    if (mbc == BoundaryFlag::reflect) return M1_IBC_REFLECT;
    return M1_IBC_MARSHAK;
  }
  if (s.compare("marshak") == 0) return M1_IBC_MARSHAK;
  if (s.compare("flux") == 0) return M1_IBC_FLUX;
  if (s.compare("reflect") == 0) return M1_IBC_REFLECT;
  if (s.compare("periodic") == 0) return M1_IBC_PERIODIC;
  if (s.compare("efix") == 0) return M1_IBC_EFIX;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    << std::endl << "<rad_m1>/implicit_bc_x1min|max = '" << s << "' is not a valid "
    << "choice (auto | marshak | flux | reflect | periodic | efix)" << std::endl;
  std::exit(EXIT_FAILURE);
  return M1_IBC_MARSHAK;
}

//----------------------------------------------------------------------------------------
//! \fn ImplFatal
//! \brief one fatal-error exit with a message, used by the 3a2 option parsers

void ImplFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitInit
//! \brief read the <rad_m1> parameters of the implicit solver, check the restrictions of
//! 3a and allocate the face array and the work array.  Called at the END of the
//! constructor, from the delimited hook there, and a no-op in explicit mode.

void RadiationM1::ImplicitInit(ParameterInput *pin) {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;

  // ---- MILESTONE 3b phase B.  `implicit` = the full 7-point solve; `implicit_x1` keeps
  // every branch below on the 3a/3a2/3c arithmetic, bit for bit.
  const bool full = (transport == M1_TRANSPORT_IMPLICIT);
  trans_on = false;
  trans_x3 = false;
  impl_cfl = pin->GetOrAddReal("rad_m1","implicit_cfl",-1.0);
  impl_tol = pin->GetOrAddReal("rad_m1","implicit_tol",1.0e-8);
  // a lagged transverse coupling needs more outer passes than a pure column solve, so
  // the ceiling is raised (not the count: the loop still stops at convergence).
  impl_maxit = pin->GetOrAddInteger("rad_m1","implicit_maxit", full ? 200 : 30);
  impl_lin_tol = pin->GetOrAddReal("rad_m1","implicit_lin_tol",1.0e-10);
  impl_lin_maxit = pin->GetOrAddInteger("rad_m1","implicit_lin_maxit",200);
  {std::string sv = pin->GetOrAddString("rad_m1","implicit_solver","line_jacobi");
  if (sv.compare("line_jacobi") == 0) {
    impl_solver = M1_ISOLV_LINE_JACOBI;
  } else if (sv.compare("bicgstab") == 0) {
    impl_solver = M1_ISOLV_BICGSTAB;
  } else {
    ImplFatal("<rad_m1>/implicit_solver = '" + sv
              + "' is not a choice (line_jacobi | bicgstab)");
  }
  }
  // ---- MILESTONE 3b phase D: what happens to the off-diagonal Eddington terms, and how
  // fast the closure is allowed to move between Picard passes.  The defaults reproduce
  // phase C on every 1-D / implicit_x1 configuration (there are no off-diagonal terms and
  // no transverse closure there), and `operator` is the default in multi-D, where lagging
  // them has no fixed point in an optically thin cell.
  std::string sod = pin->GetOrAddString("rad_m1","implicit_offdiag","auto");
  bool od_auto = (sod.compare("auto") == 0);
  if (od_auto || sod.compare("lagged") == 0) {
    impl_offdiag = M1_OD_LAGGED;
  } else if (sod.compare("operator") == 0) {
    impl_offdiag = M1_OD_OPERATOR;
  } else if (sod.compare("none") == 0) {
    impl_offdiag = M1_OD_NONE;
  } else {
    ImplFatal("<rad_m1>/implicit_offdiag = '" + sod
              + "' is not a choice (auto | lagged | operator | none)");
  }
  impl_crelax = pin->GetOrAddReal("rad_m1","implicit_closure_relax",1.0);
  if (!(impl_crelax > 0.0) || impl_crelax > 1.0) {
    ImplFatal("<rad_m1>/implicit_closure_relax must lie in (0,1]");
  }
  impl_crelax_thin = pin->GetOrAddBoolean("rad_m1","implicit_closure_relax_thin",false);
  {std::string sc = pin->GetOrAddString("rad_m1","implicit_closure_lag","pass");
  if (sc.compare("pass") == 0) {
    impl_clag_step = false;
  } else if (sc.compare("step") == 0) {
    impl_clag_step = true;
  } else {
    ImplFatal("<rad_m1>/implicit_closure_lag = '" + sc
              + "' is not a choice (pass | step)");
  }
  }
  // ---- the TRANSVERSE realizability limiter.  `none` is not merely the default: it
  // takes the OLD expression for theta everywhere (see ImplicitTransTheta), so an input
  // file that does not name it is bitwise unchanged.
  {std::string st = pin->GetOrAddString("rad_m1","implicit_trans_limit","none");
  if (st.compare("none") == 0) {
    impl_tlim = M1_TLIM_NONE;
  } else if (st.compare("lp") == 0) {
    impl_tlim = M1_TLIM_LP;
  } else {
    ImplFatal("<rad_m1>/implicit_trans_limit = '" + st
              + "' is not a choice (none | lp)");
  }
  }
  impl_tfmax = pin->GetOrAddReal("rad_m1","implicit_trans_fmax",1.0);
  if (!(impl_tfmax > 0.0)) {
    ImplFatal("<rad_m1>/implicit_trans_fmax must be positive");
  }
  // ---- MILESTONE 3e: ANDERSON acceleration of the outer (Picard) iteration.  `none` is
  // not merely the default: nothing below is allocated and ImplicitAccelSave/Apply are
  // never called, so an input file that does not name it is bitwise unchanged.
  {std::string sa = pin->GetOrAddString("rad_m1","implicit_accel","none");
  if (sa.compare("none") == 0) {
    impl_accel = M1_IACC_NONE;
  } else if (sa.compare("anderson") == 0) {
    impl_accel = M1_IACC_ANDERSON;
  } else {
    ImplFatal("<rad_m1>/implicit_accel = '" + sa
              + "' is not a choice (none | anderson)");
  }
  }
  // ---- MILESTONE 3g: the GAS-RADIATION energy coupling.  Both default false; neither
  // allocates anything nor is referenced then, so an input file that does not name them
  // is bitwise unchanged.
  impl_gas_newton = pin->GetOrAddBoolean("rad_m1","implicit_gas_newton",false);
  impl_eos_cache = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache",false);
  impl_ecnt = pin->GetOrAddInteger("rad_m1","implicit_eos_cache_nt",2);
  impl_eccheck = pin->GetOrAddBoolean("rad_m1","implicit_eos_cache_check",true);
  // the x1 LINE SOLVE (preconditioner and line-Jacobi pass): thomas = one thread per
  // column, serial recurrence (the original); pcr = one team per column, parallel
  // cyclic reduction in team scratch (GPU).  Same system; the answers agree to
  // round-off, not bitwise.  The gathered stack sweep (part_nblk > 1) is Thomas always.
  // DEFAULT pcr since bench/m1_defaults_0923 (15x faster line solve on the GPU, +1 %
  // wall on the host; thomas reproduces the earlier default bitwise).
  {std::string ls = pin->GetOrAddString("rad_m1","implicit_line_solver","pcr");
  if (ls.compare("thomas") == 0) {
    impl_line_solver = 0;
  } else if (ls.compare("pcr") == 0) {
    impl_line_solver = 1;
  } else {
    ImplFatal("<rad_m1>/implicit_line_solver = '" + ls
              + "' is not a choice (thomas | pcr)");
  }
  }
  impl_pcr_team = pin->GetOrAddInteger("rad_m1","implicit_pcr_team",0);
  impl_pcr_check = pin->GetOrAddBoolean("rad_m1","implicit_pcr_check",false);
  if (impl_pcr_team < 0 || impl_pcr_team > 1024) {
    ImplFatal("<rad_m1>/implicit_pcr_team must lie in [0,1024]");
  }
  // SPEED-UP 3: the host synchronisations of the BiCGStab loop.  0 = the original loop;
  // 1 = the same recurrence with fused reductions and vector updates (3 blocking
  // reductions per iteration instead of 5); 2 = 1 with alpha kept on the device, so
  // rhat.v does not block either (1 rank only; with more ranks 2 acts as 1).  Levels
  // 1 and 2 sum in a different order than 0: same answer to round-off, not bitwise.
  // DEFAULT 1 since bench/m1_defaults_0923 (0 reproduces the earlier default).
  impl_bcg_sync = pin->GetOrAddInteger("rad_m1","implicit_bcg_sync",1);
  if (impl_bcg_sync < 0 || impl_bcg_sync > 2) {
    ImplFatal("<rad_m1>/implicit_bcg_sync must be 0, 1 or 2");
  }
  if (impl_bcg_sync == 2) {
    bcg_rvd = Kokkos::View<Real, DevMemSpace>("m1_bcg_rvd");
  }
  // GPU COST OF THE KRYLOV ITERATION (bench/m1_fast_0923).  Off, nothing is allocated
  // or called and the path is bitwise the pre-0923 code.  Since
  // tests_m1/runs_3p_fastdefault they DEFAULT ON (the `FAST` set of
  // tests_m1/runs_3k_gpu3d/README_FAST.md: halo_direct, od_cache, krylov_fuse = 3,
  // op_stencil, precond = rbgs_fwd; 2.0-2.4x on the GPU, the converged answer moved at
  // the +-1-ulp level) for the closures whose tensor is fixed within a step (eddington,
  // vet_sc, tau) whenever the configuration admits them: transport = implicit with
  // bicgstab, implicit_bcg_sync = 1, implicit_line_solver = pcr, ONE MeshBlock along x1,
  // no implicit_lin_cnorm, and (op_stencil) no periodic x1 wrap.  m1 / minerbo /
  // kershaw keep them off (not gated there).  A key the input names keeps its value.
  //  implicit_halo_direct = true   the implicit exchanges as ONE on-rank copy kernel
  //    when every neighbour of every block is on its own rank at the same level
  //    (otherwise the ordinary exchange, as before).  Bitwise.
  //  implicit_od_cache = true      the off-diagonal Eddington operator from a per-cell
  //    cache of sum_e d_e P_de, fused with the 7-point row into one kernel.  Bitwise.
  //  implicit_krylov_fuse = 1      the pcr preconditioner reads and writes the Krylov
  //    vectors itself and carries the p / s updates (bitwise); = 2 also puts the
  //    (rhat,v) and (t,s),(t,t) reductions in the operator kernel (round-off: the sums
  //    are ordered differently).  Needs implicit_bcg_sync = 1, implicit_line_solver =
  //    pcr, one block along x1; 2 needs implicit_od_cache.  3 = 2 with TWO blocking
  //    reductions per iteration (ImplicitBiCGStabTwo; round-off).
  //  implicit_op_stencil = true    the frozen operator of each Picard pass written out
  //    once as a 19-point stencil (7-point row + off-diagonal Eddington terms), applied
  //    by one kernel per Krylov product (round-off: the terms are grouped differently).
  //    Needs implicit_bcg_sync = 1; not with a periodic x1 wrap.
  //  implicit_precond_float = true  the line solves of the fused path's preconditioner
  //    in float (the operator, the vectors and every reduction stay double: a
  //    preconditioner only sets the convergence rate).  Needs implicit_krylov_fuse >= 1.
  //  implicit_precond = line | rbgs | rbgs_fwd   the preconditioner of the fused path:
  //    x1 line Jacobi (the original), symmetric / forward red-black transverse line
  //    Gauss-Seidel, block-local (ImplicitPrecondX).  Needs implicit_krylov_fuse >= 1.
  const bool fixcl = eddington || vet_sc || tau_closure;
  bool fdef = fixcl && full && (impl_solver == M1_ISOLV_BICGSTAB) &&
              (impl_bcg_sync == 1) && (impl_line_solver == 1) &&
              (pmy_pack->pmesh->mesh_indcs.nx1 == pmy_pack->pmesh->mb_indcs.nx1);
  if (pin->DoesParameterExist("rad_m1","implicit_lin_cnorm") &&
      pin->GetReal("rad_m1","implicit_lin_cnorm") > 0.0) {
    fdef = false;
  }
  // the x1 wrap as the boundary parsing below will read it (same keys, same defaults)
  const bool x1per = (ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      pmy_pack->pmesh->mesh_bcs[static_cast<int>(BoundaryFace::inner_x1)])
      == M1_IBC_PERIODIC);
  impl_halo_direct = pin->GetOrAddBoolean("rad_m1","implicit_halo_direct",fdef);
  impl_odc = pin->GetOrAddBoolean("rad_m1","implicit_od_cache",fdef);
  impl_kfuse = pin->GetOrAddInteger("rad_m1","implicit_krylov_fuse",fdef ? 3 : 0);
  impl_stencil = pin->GetOrAddBoolean("rad_m1","implicit_op_stencil",fdef && !x1per);
  impl_prec_float = pin->GetOrAddBoolean("rad_m1","implicit_precond_float",false);
  if (impl_stencil && impl_bcg_sync != 1) {
    ImplFatal("<rad_m1>/implicit_op_stencil needs implicit_bcg_sync = 1");
  }
  {std::string pc = pin->GetOrAddString("rad_m1","implicit_precond",
                                        fdef ? "rbgs_fwd" : "line");
  if (pc.compare("line") == 0) {
    impl_prec = 0;
  } else if (pc.compare("rbgs") == 0) {
    impl_prec = 1;
  } else if (pc.compare("rbgs_fwd") == 0) {
    impl_prec = 2;
  } else {
    ImplFatal("<rad_m1>/implicit_precond = '" + pc
              + "' is not a choice (line | rbgs | rbgs_fwd)");
  }
  }
  if (impl_kfuse < 0 || impl_kfuse > 3) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse must be 0, 1, 2 or 3");
  }
  if (impl_kfuse == 3 && pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0) > 0.0) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse = 3 does not take implicit_lin_cnorm");
  }
  if (impl_kfuse > 0 && (impl_bcg_sync != 1 || impl_line_solver != 1)) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse needs implicit_bcg_sync = 1 and "
              "implicit_line_solver = pcr");
  }
  if (impl_kfuse >= 2 && !impl_odc) {
    ImplFatal("<rad_m1>/implicit_krylov_fuse >= 2 needs implicit_od_cache = true");
  }
  if (impl_prec_float && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond_float needs implicit_krylov_fuse >= 1");
  }
  if (impl_prec > 0 && impl_kfuse == 0) {
    ImplFatal("<rad_m1>/implicit_precond = rbgs needs implicit_krylov_fuse >= 1");
  }
  // multi-rank Krylov (rad_m1_krylov.cpp, tests_m1/runs_3w_krylov): both default OFF
  impl_kpipe = pin->GetOrAddBoolean("rad_m1","implicit_krylov_pipe",false);
  impl_halo_mpi = pin->GetOrAddBoolean("rad_m1","implicit_halo_mpi",false);
  hm_state = 0;
  hm_comm = nullptr;
  if (impl_kpipe && impl_kfuse != 3) {
    ImplFatal("<rad_m1>/implicit_krylov_pipe needs implicit_krylov_fuse = 3");
  }
  if (impl_halo_mpi && !impl_halo_direct) {
    ImplFatal("<rad_m1>/implicit_halo_mpi needs implicit_halo_direct = true");
  }
  // the Picard pass count (bench/m1_picard_0923): a per-pass log, off by default
  impl_plog = pin->GetOrAddInteger("rad_m1","implicit_picard_log",0);
  // ...and the options that cut it.  Since bench/m1_defaults_0923 they DEFAULT ON for
  // the closures that do not read the iterate (eddington, vet_sc, tau: 1.6x on the GPU,
  // statistics moved at the solver-tolerance level) and stay OFF for m1 / minerbo /
  // kershaw, whose closure moves with the iterate: there the lresid test is NOT
  // redundant and the inexact (Eisenstat-Walker) pass 0 raised NON-CONVERGED steps
  // (tests_m1/runs_3k_gpu3d/README_DEFAULTS.md).  The earlier default path is
  // implicit_lres_test = true, implicit_conv_est = false, implicit_lin_ew_max = 0,
  // implicit_predictor = none (with implicit_bcg_sync = 0 and implicit_line_solver =
  // thomas it is bitwise the pre-flip code):
  //  implicit_lres_test = false  drop the pass-to-pass transverse-change test under
  //    bicgstab, where the transverse coupling is IN the operator and is solved to
  //    implicit_lin_tol on every pass (the test lags the Picard test by one pass);
  //  implicit_conv_est = true    stop once q/(1-q) times the last change (q = the
  //    measured contraction of the last two passes, required < 1/2) is below
  //    implicit_tol, i.e. without the confirming pass;
  //  implicit_lin_ew_max > 0     Eisenstat-Walker (choice 2) inner tolerance
  //    max(lin_tol max|b|, eta_k max|r0|), eta_0 = ew_max,
  //    eta_k = min(ew_max, gamma (|r0_k|/|r0_{k-1}|)^2) with the gamma eta_{k-1}^2
  //    safeguard (bcg_sync >= 1 only);
  //  implicit_predictor = step   start the loop from the previous step's implicit
  //    increment scaled by dt/dt_prev (closures that do not read the iterate only;
  //    the default for them since 0923, restart-safe).
  //  The Eisenstat-Walker default is 1e-2 with bcg_sync >= 1 and 0 (off) with
  //  bcg_sync = 0, so an input that asks for the original BiCGStab loop still runs.
  impl_lres_test = pin->GetOrAddBoolean("rad_m1","implicit_lres_test",!fixcl);
  impl_conv_est = pin->GetOrAddBoolean("rad_m1","implicit_conv_est",fixcl);
  impl_ew_max = pin->GetOrAddReal("rad_m1","implicit_lin_ew_max",
                                  (fixcl && impl_bcg_sync >= 1) ? 1.0e-2 : 0.0);
  impl_ew_gam = pin->GetOrAddReal("rad_m1","implicit_lin_ew_gamma",0.9);
  //  implicit_lin_cnorm > 0     the inner test on max_i |r_i|/(s_i E^k_i) < cnorm,
  //    s_i = 1 + SRCB_i the row excess: a bound on the per-cell relative error of E
  //    the residual implies (bcg_sync >= 1 only).
  impl_lin_cnorm = pin->GetOrAddReal("rad_m1","implicit_lin_cnorm",0.0);
  if (impl_lin_cnorm < 0.0 || (impl_lin_cnorm > 0.0 && impl_bcg_sync == 0)) {
    ImplFatal("<rad_m1>/implicit_lin_cnorm must be >= 0 and needs "
              "implicit_bcg_sync >= 1");
  }
  if (impl_ew_max < 0.0 || impl_ew_max >= 1.0 || !(impl_ew_gam > 0.0)) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max must lie in [0,1), ew_gamma > 0");
  }
  if (impl_ew_max > 0.0 && impl_bcg_sync == 0) {
    ImplFatal("<rad_m1>/implicit_lin_ew_max needs implicit_bcg_sync >= 1");
  }
  {std::string pr = pin->GetOrAddString("rad_m1","implicit_predictor",
                                        fixcl ? "step" : "none");
  // predictor = step is the default for the closures that do not read the iterate: its
  // state (ipred, pred_ok, pred_dt) travels in the restart file (radm1::kM1PredRstMagic,
  // tests_m1/runs_3k_gpu3d/README_PREDRST.md), so a restart stays bitwise.
  if (pr.compare("none") == 0) {
    impl_pred = false;
  } else if (pr.compare("step") == 0) {
    impl_pred = true;
  } else {
    ImplFatal("<rad_m1>/implicit_predictor = '" + pr
              + "' is not a choice (none | step)");
  }
  }
  if (impl_ecnt < 0 || impl_ecnt > M1_EC_NTMAX) {
    ImplFatal("<rad_m1>/implicit_eos_cache_nt must lie in [0,8]");
  }
  impl_and_m = pin->GetOrAddInteger("rad_m1","implicit_anderson_m",5);
  impl_and_beta = pin->GetOrAddReal("rad_m1","implicit_anderson_beta",1.0);
  impl_and_start = pin->GetOrAddInteger("rad_m1","implicit_anderson_start",1);
  if (impl_accel == M1_IACC_ANDERSON) {
    if (impl_and_m < 1 || impl_and_m > M1_AND_MMAX) {
      ImplFatal("<rad_m1>/implicit_anderson_m must lie in [1,10]");
    }
    if (!(impl_and_beta > 0.0) || impl_and_beta > 1.0) {
      ImplFatal("<rad_m1>/implicit_anderson_beta must lie in (0,1]");
    }
    if (impl_and_start < 1) {
      ImplFatal("<rad_m1>/implicit_anderson_start must be >= 1 (pass 0 has no history)");
    }
  }
  impl_opac_update = pin->GetOrAddBoolean("rad_m1","implicit_opac_update",false);
  impl_allow_multid = pin->GetOrAddBoolean("rad_m1","implicit_allow_multid",false);
  marshak_q = pin->GetOrAddReal("rad_m1","marshak_q",0.5);
  // ---- milestone 3a2 options.  All three default to the 3a behaviour, so an input file
  // that does not name them reproduces RESULTS.txt of runs_3a exactly.
  std::string sfx = pin->GetOrAddString("rad_m1","implicit_flux","central");
  if (sfx.compare("central") == 0) {
    impl_flux = M1_IFLUX_CENTRAL;
  } else if (sfx.compare("ap_hll") == 0) {
    impl_flux = M1_IFLUX_APHLL;
  } else if (sfx.compare("berthon") == 0) {
    impl_flux = M1_IFLUX_BERTHON;
  } else if (sfx.compare("blend") == 0) {
    impl_flux = M1_IFLUX_BLEND;
  } else {
    ImplFatal("<rad_m1>/implicit_flux = '" + sfx
              + "' is not a choice (central | ap_hll | berthon | blend)");
  }
  // 3b phase B: the transverse operator is built for the face-eliminated (central) form
  // only.  The HLL/berthon/blend face fluxes carry per-face coefficients (ifw) that exist
  // for the x1 faces alone, so anything but `central` would silently be central in x2/x3.
  if (full && impl_flux != M1_IFLUX_CENTRAL) {
    ImplFatal("<rad_m1>/transport = implicit supports implicit_flux = central only "
              "(the asymptotic-preserving forms are built for the x1 faces); use "
              "transport = implicit_x1 for ap_hll | berthon | blend");
  }
  // ---- milestone 3c.  The weight of implicit_flux = blend.  Inert for every other
  // flux, and the two ends of the blend are BITWISE central and berthon.
  std::string sbl = pin->GetOrAddString("rad_m1","implicit_blend","tau_f");
  if (sbl.compare("tau") == 0) {
    impl_blend = M1_IBLEND_TAU;
  } else if (sbl.compare("f") == 0) {
    impl_blend = M1_IBLEND_F;
  } else if (sbl.compare("tau_f") == 0) {
    impl_blend = M1_IBLEND_TAUF;
  } else {
    ImplFatal("<rad_m1>/implicit_blend = '" + sbl
              + "' is not a choice (tau | f | tau_f)");
  }
  std::string sbm = pin->GetOrAddString("rad_m1","implicit_blend_fmode","max");
  if (sbm.compare("max") == 0) {
    impl_blend_fmode = M1_IBFM_MAX;
  } else if (sbm.compare("mean") == 0) {
    impl_blend_fmode = M1_IBFM_MEAN;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_fmode = '" + sbm
              + "' is not a choice (max | mean)");
  }
  std::string sbw = pin->GetOrAddString("rad_m1","implicit_blend_mode","flux");
  if (sbw.compare("flux") == 0) {
    impl_blend_mode = M1_IBMODE_FLUX;
  } else if (sbw.compare("dissipation") == 0) {
    impl_blend_mode = M1_IBMODE_DISSIP;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_mode = '" + sbw
              + "' is not a choice (flux | dissipation)");
  }
  impl_blend_tau0 = pin->GetOrAddReal("rad_m1","implicit_blend_tau0",1.0);
  impl_blend_flo = pin->GetOrAddReal("rad_m1","implicit_blend_flo",0.6);
  impl_blend_fhi = pin->GetOrAddReal("rad_m1","implicit_blend_fhi",0.9);
  if (!(impl_blend_tau0 > 0.0) || !(impl_blend_fhi > impl_blend_flo)) {
    ImplFatal("<rad_m1>: implicit_blend_tau0 must be positive and implicit_blend_fhi "
              "must exceed implicit_blend_flo");
  }
  std::string srn = pin->GetOrAddString("rad_m1","implicit_recon","dc");
  if (srn.compare("dc") == 0) {
    impl_recon = M1_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = M1_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_m1>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // implicit_enthalpy (see rad_m1_implicit.hpp).  Read only when it is named, so that
  // the parameter dump of an input that does not name it is unchanged.
  impl_enth = M1_IENTH_UPWIND;
  if (pin->DoesParameterExist("rad_m1","implicit_enthalpy")) {
    std::string sen = pin->GetString("rad_m1","implicit_enthalpy");
    if (sen.compare("upwind") == 0) {
      impl_enth = M1_IENTH_UPWIND;
    } else if (sen.compare("central") == 0) {
      impl_enth = M1_IENTH_CENTRAL;
    } else if (sen.compare("plm") == 0) {
      impl_enth = M1_IENTH_PLM;
    } else {
      ImplFatal("<rad_m1>/implicit_enthalpy = '" + sen
                + "' is not a choice (upwind | central | plm)");
    }
  }
  // LIMIT 4 of the 3a findings is NOT implemented in 3a2: a column still has to live
  // inside one MeshBlock along x1 (the fatal below).  The option is parsed so that the
  // input files and the gate scripts can already name it, and `gather` fatals rather
  // than silently doing something else.
  std::string spt = pin->GetOrAddString("rad_m1","implicit_partition","none");
  if (spt.compare("none") == 0) {
    impl_part = M1_IPART_NONE;
  } else if (spt.compare("gather") == 0) {
    impl_part = M1_IPART_GATHER;
  } else {
    ImplFatal("<rad_m1>/implicit_partition = '" + spt
              + "' is not a choice (none | gather)");
  }
  impl_recon_w = pin->GetOrAddReal("rad_m1","implicit_recon_w",-1.0);
  impl_res_floor = pin->GetOrAddReal("rad_m1","implicit_res_floor",0.0);
  std::string slg = pin->GetOrAddString("rad_m1","implicit_recon_lag","picard");
  impl_recon_freeze = (slg.compare("step") == 0);
  // MILESTONE 3c: freeze the deferred correction, and with it the plm limiter's choice,
  // after this many Picard passes.  The 3a2 finding is that the limiter keeps switching
  // on a handful of cells and the iteration is a small limit cycle that never meets
  // implicit_tol, so plm_dc costs implicit_maxit passes per step; freezing the
  // correction after a few passes leaves an ordinary linear system to converge.
  // <= 0 (the default) never freezes, i.e. reproduces 3a2.
  impl_recon_npass = pin->GetOrAddInteger("rad_m1","implicit_recon_npass",-1);
  if (!impl_recon_freeze && slg.compare("picard") != 0) {
    ImplFatal("<rad_m1>/implicit_recon_lag = '" + slg
              + "' is not a choice (step | picard)");
  }
  // LIMIT 3 of the 3a findings: the DEFAULT is now `true`.  It is a bug fix, not a
  // tuning knob (the boundary face used to hand its whole momentum to one interior cell,
  // which gave that cell 1.5 face-shares of radiative force: He column bottom-cell |v1|
  // 10.3 -> 0.47 v_MLT).  The key is kept so that `false` reproduces runs_3a/RESULTS.txt.
  impl_bmom_half = pin->GetOrAddBoolean("rad_m1","implicit_bmom_half",true);
  if (impl_lin_maxit < 1) {
    ImplFatal("<rad_m1>/implicit_lin_maxit must be >= 1");
  }
  if (!(impl_tol > 0.0) || impl_maxit < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/implicit_tol must be positive and implicit_maxit >= 1"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &mbcs = pmy_pack->pmesh->mesh_bcs;
  ibc_x1min = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      mbcs[static_cast<int>(BoundaryFace::inner_x1)]);
  ibc_x1max = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1max","auto"),
      mbcs[static_cast<int>(BoundaryFace::outer_x1)]);
  iflux_x1min = pin->GetOrAddReal("rad_m1","implicit_flux_x1min",0.0);
  iflux_x1max = pin->GetOrAddReal("rad_m1","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_m1","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_m1","implicit_ebath_x1max",0.0);
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- the restrictions of 3a
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mindcs = pmy_pack->pmesh->mesh_indcs;
  part_nblk = 1;
  if (mindcs.nx1 != indcs.nx1) {
    if (impl_part != M1_IPART_GATHER) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/transport = implicit_x1 with implicit_partition = none "
        << "needs exactly ONE MeshBlock along x1: <meshblock>/nx1 must equal <mesh>/nx1 ("
        << indcs.nx1 << " vs " << mindcs.nx1 << ").  Set <rad_m1>/implicit_partition = "
        << "gather for the line solve partitioned over MeshBlocks and ranks" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((mindcs.nx1 % indcs.nx1) != 0) {
      ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nx1 to be an exact "
                "multiple of <meshblock>/nx1 (uniform mesh only)");
    }
    part_nblk = mindcs.nx1/indcs.nx1;
  }
  if (part_nblk > 1 && (ibc_x1min == M1_IBC_PERIODIC)) {
    ImplFatal("<rad_m1>/implicit_partition = gather does not support PERIODIC x1 across "
              "more than one MeshBlock (the cyclic Thomas sweep of 3a wraps inside one "
              "block).  Use one MeshBlock along x1, or a non-periodic x1 boundary pair");
  }
  if (part_nblk > 1 && indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nghost >= 2");
  }
  if (full) {
    trans_on = pmy_pack->pmesh->multi_d;
    trans_x3 = pmy_pack->pmesh->three_d;
    if (!trans_on) {
      // a 1-D mesh has no transverse direction: the solve IS the x1 column solve, and
      // gate G3 is exactly this statement.
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1>: transport=implicit on a 1-D mesh is the x1 column solve"
                  << std::endl;
      }
    }
  }
  if ((indcs.nx2 > 1 || indcs.nx3 > 1) && !impl_allow_multid && !full) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does NOT transport along x2/x3. "
      << "Set <rad_m1>/implicit_allow_multid = true to run a set of INDEPENDENT x1 "
      << "columns (F_2 = F_3 = 0 everywhere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does not support SMR/AMR"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // one solve per hydro step: no sub-cycling, one stage
  subcycle = false;
  nstage = 1;
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (impl_cfl > 0.0) || (!other_sets_dt);

  // ---- arrays
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(f0x1, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1, 0.0);
  Kokkos::realloc(f0x1n, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1n, 0.0);
  // MILESTONE 3b phase C: the BiCGStab wrapper needs the four/six transverse off-diagonal
  // coefficients and ten Krylov vectors on top of the phase-B work array.  They are
  // allocated only when it is selected AND the mesh is multi-D, so line_jacobi keeps the
  // footprint (and, on a 1-D mesh, the solve IS the column solve and there is no system
  // to wrap).
  bicg_on = trans_on && (impl_solver == M1_ISOLV_BICGSTAB);
  // MILESTONE 3b phase D.  `auto` = the off-diagonal Eddington terms go INTO the operator
  // wherever there is a Krylov solver to carry them, and stay LAGGED otherwise, which is
  // what every line_jacobi and implicit_x1 configuration did before phase D.
  if (od_auto) {
    impl_offdiag = bicg_on ? M1_OD_OPERATOR : M1_OD_LAGGED;
  }
  if (impl_offdiag == M1_OD_OPERATOR && trans_on && !bicg_on) {
    ImplFatal("<rad_m1>/implicit_offdiag = operator needs implicit_solver = bicgstab: "
              "the off-diagonal Eddington terms are a 9-/19-point coupling and the "
              "x1 line solve of line_jacobi cannot carry them");
  }
  od_now = impl_offdiag;
  int niw = full ? (bicg_on ? M1_NIW_K : M1_NIW) : M1_NIW_X1;
  // MILESTONE 3g: the five per-cell components of the gas coupling are APPENDED, so
  // every index above keeps the value it had and the array grows only when asked for.
  iw_gas = -1;
  if (impl_gas_newton || impl_eos_cache) {
    iw_gas = niw;
    niw += M1_NIW_GAS;
  }
  if (impl_eos_cache) {
    impl_nec = M1EosCacheNComp(impl_ecnt);
    Kokkos::realloc(ecache, nmb, impl_nec, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ecache, -1.0);
  }
  Kokkos::realloc(iw, nmb, niw, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(iw, 0.0);
  if (impl_pred) {
    Kokkos::realloc(ipred, nmb, 3, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(ipred, 0.0);
    pred_ok = false;
  }
  Kokkos::realloc(ifw, nmb, M1_NIFW, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(ifw, 0.0);
  if (trans_on) {
    Kokkos::realloc(f0x2, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2, 0.0);
    Kokkos::realloc(f0x2n, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::deep_copy(f0x2n, 0.0);
    if (trans_x3) {
      Kokkos::realloc(f0x3, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3, 0.0);
      Kokkos::realloc(f0x3n, nmb, ncells3+1, ncells2, ncells1);
      Kokkos::deep_copy(f0x3n, 0.0);
    }
    if (impl_tlim != M1_TLIM_NONE) {
      Kokkos::realloc(thx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(thx2, 0.0);
      Kokkos::realloc(klx2, nmb, ncells3, ncells2+1, ncells1);
      Kokkos::deep_copy(klx2, 0.0);
      if (trans_x3) {
        Kokkos::realloc(thx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(thx3, 0.0);
        Kokkos::realloc(klx3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::deep_copy(klx3, 0.0);
      }
    }
    // the transverse halo rides the module's ORDINARY cell-centred boundary machinery on
    // a scratch array, which is what gives it periodic wrap, corner/edge neighbours and
    // MPI for free; the hand-rolled x1 halo of 3b is then not used at all under
    // transport = implicit (this exchange carries its six quantities and seven more).
    Kokkos::realloc(thw, nmb, M1_NHALO_T, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thw, 0.0);
    Kokkos::realloc(thw_c, nmb, M1_NHALO_T, 1, 1, 1);
    pbval_th = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_th->InitializeBuffers(M1_NHALO_T);
    // the NARROW exchange of the same list: the M1_NHALO_Q components a Picard pass can
    // move once the closure is frozen.  A separate array because the exchange takes the
    // variable count from the array's second extent, and a prefix subview of a
    // LayoutRight array is not one.
    Kokkos::realloc(thq, nmb, M1_NHALO_Q, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(thq, 0.0);
    Kokkos::realloc(thq_c, nmb, M1_NHALO_Q, 1, 1, 1);
    pbval_tq = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    pbval_tq->InitializeBuffers(M1_NHALO_Q);
    // the deep interior of the scratch arrays is neither read by a send nor written by a
    // receive when every neighbour is at the SAME level (a same-level buffer reaches ng
    // cells in from the active boundary, buffs_cc.cpp) and neither the cubed-sphere
    // resample nor the polar transform is in play; then the copies to and from iw can
    // skip it.  SMR/AMR is already a fatal above.
    halo_shell = !(pmy_pack->pmesh->multilevel || pmy_pack->pmesh->use_cubed_sphere ||
                   pmy_pack->pmesh->use_polar_boundary);
    {
      // ONE more exchange object, for the single Krylov vector the operator application
      // needs in its ghost zones -- and for the one-component exchange of E alone, which
      // is why it is allocated under every solver.  It is used strictly SEQUENTIALLY
      // with pbval_th and pbval_tq (each exchange runs its
      // InitRecv/Send/Recv/Clear chain to completion before the next
      // starts) and every rank issues the identical SEQUENCE of exchanges -- the Picard
      // count, the BiCGStab count and every breakdown decision are taken from GLOBAL
      // reductions -- so MPI's non-overtaking guarantee keeps the two streams apart even
      // though they share the tag space.
      Kokkos::realloc(krw, nmb, 1, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(krw, 0.0);
      Kokkos::realloc(krw_c, nmb, 1, 1, 1, 1);
      pbval_kr = new MeshBoundaryValuesCC(pmy_pack, pin, false);
      pbval_kr->InitializeBuffers(1);
    }
    if (impl_halo_direct) {ImplicitHaloDirectInit();}
    if (impl_odc) {
      Kokkos::realloc(odc, nmb, 3, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(odc, 0.0);
    }
    if (impl_stencil) {
      if (ibc_x1min == M1_IBC_PERIODIC) {
        ImplFatal("<rad_m1>/implicit_op_stencil does not take a periodic x1 wrap");
      }
      Kokkos::realloc(ost, nmb, 19, ncells3, ncells2, ncells1);
      Kokkos::deep_copy(ost, 0.0);
    }
  }
  // MILESTONE 3e: the Anderson histories.  Allocated ONLY when the acceleration is on.
  if (impl_accel != M1_IACC_NONE) {
    aa_nc = trans_on ? 4 : 2;
    Kokkos::realloc(aa_xc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xc, 0.0);
    Kokkos::realloc(aa_fc, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fc, 0.0);
    Kokkos::realloc(aa_xp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_xp, 0.0);
    Kokkos::realloc(aa_fp, nmb, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_fp, 0.0);
    Kokkos::realloc(aa_dx, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_dx, 0.0);
    Kokkos::realloc(aa_df, nmb, impl_and_m, aa_nc, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_df, 0.0);
    Kokkos::realloc(aa_sc, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(aa_sc, 1.0);
  }
  ImplicitPartitionInit();

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: transport=" << (full ? "implicit" : "implicit_x1")
              << (full ? (trans_x3 ? " (3-D)" : (trans_on ? " (2-D)" : " (1-D)")) : "")
              << " solver=" << (full ? (bicg_on ? "bicgstab" : "line_jacobi") : "thomas")
              << " lin_tol=" << impl_lin_tol
              << " lin_maxit=" << impl_lin_maxit << std::endl;
    std::cout << "         implicit_cfl="
              << impl_cfl << " tol=" << impl_tol << " maxit=" << impl_maxit
              << " opac_update=" << (impl_opac_update ? "true" : "false")
              << " marshak_q=" << marshak_q << std::endl;
    std::cout << "         implicit_flux="
              << ((impl_flux == M1_IFLUX_APHLL) ? "ap_hll" :
                  ((impl_flux == M1_IFLUX_BERTHON) ? "berthon" :
                   ((impl_flux == M1_IFLUX_BLEND) ? "blend" : "central")))
              << ((impl_flux == M1_IFLUX_BLEND) ? (":" + sbl + "/" + sbm + "/" + sbw)
                                                : std::string(""))
              << " implicit_recon="
              << ((impl_recon == M1_IRECON_PLMDC) ? "plm_dc" : "dc")
              << " implicit_partition="
              << ((impl_part == M1_IPART_GATHER) ? "gather" : "none") << std::endl;
    if (impl_enth != M1_IENTH_UPWIND) {
      std::cout << "         implicit_enthalpy="
                << ((impl_enth == M1_IENTH_PLM) ? "plm" : "central")
                << " (deferred correction)" << std::endl;
    }
    if (trans_on) {
      std::cout << "         implicit_offdiag="
                << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                    ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
                << (od_auto ? " (auto)" : "")
                << " closure_relax=" << impl_crelax
                << (impl_crelax_thin ? " (thin faces only)" : "")
                << " closure_lag=" << (impl_clag_step ? "step" : "pass")
                << " trans_limit="
                << ((impl_tlim == M1_TLIM_LP) ? "lp" : "none")
                << ((impl_tlim == M1_TLIM_LP)
                    ? (" fmax=" + std::to_string(impl_tfmax)) : std::string(""))
                << std::endl;
    }
    std::cout << "         implicit_accel="
              << ((impl_accel == M1_IACC_ANDERSON) ? "anderson" : "none")
              << ((impl_accel == M1_IACC_ANDERSON)
                  ? (" m=" + std::to_string(impl_and_m)
                     + " beta=" + std::to_string(impl_and_beta)
                     + " start=" + std::to_string(impl_and_start)
                     + " ncomp=" + std::to_string(aa_nc))
                  : std::string("")) << std::endl;
    if (impl_gas_newton || impl_eos_cache) {
      std::cout << "         implicit_gas_newton="
                << (impl_gas_newton ? "true" : "false")
                << " implicit_eos_cache=" << (impl_eos_cache ? "true" : "false")
                << (impl_eos_cache ? (" nt=" + std::to_string(impl_ecnt)
                                      + " ncomp=" + std::to_string(impl_nec)
                                      + " check="
                                      + (impl_eccheck ? "true" : "false"))
                                   : std::string("")) << std::endl;
    }
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
  if (vet_sc) {VetInit(pin);}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPartitionInit
//! \brief milestone 3b, LIMIT 4: build the topology of the x1 stacks and allocate the
//! gather/scatter and halo buffers.  A "stack" is the set of MeshBlocks that share the
//! same (x2,x3) footprint, ordered by their x1 logical location; its ROOT is the block
//! with the lowest one.  Uniform mesh only (the caller has already fatalled on SMR/AMR
//! and on a non-integer block count along x1), so the stack of a block is found by a
//! scan of the global LogicalLocation list.
//!
//! Nothing is allocated and nothing is communicated when part_nblk == 1: the 3a path is
//! then bitwise what it was.

void RadiationM1::ImplicitPartitionInit() {
  part_nroot = 0;
  part_nx1g = 0;
  part_nlay = 0;
  part_nqa = M1_NHALO_A;
  part_any_mpi = false;

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int g0 = pmy_pack->gids;
  // the three per-block tables are ALWAYS allocated: ImplicitSolve reads part_pos in
  // every kernel, and with one MeshBlock per column the defaults (position 0 of a stack
  // of one) select exactly the 3a branches.
  Kokkos::realloc(part_pos, nmb);
  Kokkos::realloc(part_slot, nmb);
  Kokkos::realloc(part_nbr, 2*nmb);
  for (int m=0; m<nmb; ++m) {
    part_pos.h_view(m) = 0;
    part_slot.h_view(m) = -1;
    part_nbr.h_view(2*m) = -1;
    part_nbr.h_view(2*m+1) = -1;
  }
  if (part_nblk <= 1) {
    part_pos.modify_host();
    part_pos.sync_device();
    part_slot.modify_host();
    part_slot.sync_device();
    part_nbr.modify_host();
    part_nbr.sync_device();
    return;
  }
  part_nx1g = part_nblk*indcs.nx1;
  part_nlay = std::min(indcs.ng, 2);

  // (1) for every LOCAL block: its position in the stack, the global ids of the stack
  // members and of its two x1 neighbours.
  part_rootgid.assign(nmb, -1);
  part_rootrank.assign(nmb, -1);
  part_nbrrank.assign(2*nmb, -1);
  part_nbrgid.assign(2*nmb, -1);
  std::vector<int> stack(part_nblk);
  for (int m=0; m<nmb; ++m) {
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {stack[p] = -1;}
    for (int g=0; g<pm->nmb_total; ++g) {
      LogicalLocation &lg = pm->lloc_eachmb[g];
      if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level) {
        if (lg.lx1 >= 0 && lg.lx1 < part_nblk) {stack[lg.lx1] = g;}
      }
    }
    for (int p=0; p<part_nblk; ++p) {
      if (stack[p] < 0) {
        ImplFatal("<rad_m1>/implicit_partition = gather: the x1 stack of a MeshBlock is "
                  "incomplete (a non-uniform mesh?)");
      }
    }
    int pos = static_cast<int>(lm.lx1);
    part_pos.h_view(m) = pos;
    part_rootgid[m] = stack[0];
    part_rootrank[m] = pm->rank_eachmb[stack[0]];
    if (pos > 0) {
      part_nbrgid[2*m] = stack[pos-1];
      part_nbrrank[2*m] = pm->rank_eachmb[stack[pos-1]];
    }
    if (pos < part_nblk-1) {
      part_nbrgid[2*m+1] = stack[pos+1];
      part_nbrrank[2*m+1] = pm->rank_eachmb[stack[pos+1]];
    }
    for (int s=0; s<2; ++s) {
      int gn = part_nbrgid[2*m+s];
      bool loc = (gn >= 0) && (part_nbrrank[2*m+s] == global_variable::my_rank);
      part_nbr.h_view(2*m+s) = loc ? (gn - g0) : -1;
      if (gn >= 0 && part_nbrrank[2*m+s] != global_variable::my_rank) {
        part_any_mpi = true;
      }
    }
    if (part_rootrank[m] != global_variable::my_rank) {part_any_mpi = true;}
  }

  // (2) the local ROOTS and their member lists
  part_mrank.clear();
  part_mgid.clear();
  std::vector<int> rootmb;
  for (int m=0; m<nmb; ++m) {
    if (part_pos.h_view(m) != 0) continue;
    rootmb.push_back(m);
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {
      int gfound = -1;
      for (int g=0; g<pm->nmb_total; ++g) {
        LogicalLocation &lg = pm->lloc_eachmb[g];
        if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level && lg.lx1 == p) {
          gfound = g;
        }
      }
      part_mgid.push_back(gfound);
      part_mrank.push_back(pm->rank_eachmb[gfound]);
      if (pm->rank_eachmb[gfound] != global_variable::my_rank) {part_any_mpi = true;}
    }
  }
  part_nroot = static_cast<int>(rootmb.size());
  for (int m=0; m<nmb; ++m) {
    int sl = -1;
    if (part_rootrank[m] == global_variable::my_rank) {
      for (int s=0; s<part_nroot; ++s) {
        if (rootmb[s] + g0 == part_rootgid[m]) {sl = s;}
      }
    }
    part_slot.h_view(m) = sl;
  }
  part_pos.modify_host();
  part_pos.sync_device();
  part_slot.modify_host();
  part_slot.sync_device();
  part_nbr.modify_host();
  part_nbr.sync_device();

  // (3) buffers
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(part_sys, std::max(part_nroot,1), 6, ncells3, ncells2, part_nx1g);
  Kokkos::deep_copy(part_sys, 0.0);
  int nrow = 4*ncells3*ncells2*indcs.nx1;
  Kokkos::realloc(part_sbuf, nmb, nrow);
  Kokkos::realloc(part_rbuf, std::max(part_nroot*part_nblk,1), nrow);
  Kokkos::realloc(part_sbuf_h, nmb, nrow);
  Kokkos::realloc(part_rbuf_h, std::max(part_nroot*part_nblk,1), nrow);
  int nhal = M1_NHALO_A*part_nlay*ncells3*ncells2;
  Kokkos::realloc(part_hbuf, nmb, 4, nhal);
  Kokkos::realloc(part_hbuf_h, nmb, 4, nhal);

  if (global_variable::my_rank == 0) {
    std::cout << "         implicit_partition=gather: " << part_nblk
              << " MeshBlocks per x1 column, " << part_nx1g << " rows per gathered line, "
              << part_nlay << " halo layers" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitX1Halo
//! \brief exchange the x1 ghost layers of the work array `iw` between the MeshBlocks of
//! one x1 stack.  `eponly` picks the one-quantity set (M1_IW_EP, the new iterate, needed
//! by the face update and by the next pass) instead of the six LAGGED quantities
//! (M1HaloCompA).
//!
//! What makes the partitioned solve BITWISE identical to a single-block solve is that
//! the ghost value a block reads is the VERY NUMBER its neighbour computed, copied, and
//! never a quantity recomputed from a hydro ghost: this routine moves nothing else.

void RadiationM1::ImplicitX1Halo(bool eponly) {
  if (part_nblk <= 1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nl = part_nlay;
  int nq = eponly ? 1 : M1_NHALO_A;
  const bool ep1 = eponly;
  auto iw_ = iw;
  auto nbr_ = part_nbr;
  int nj = je - js + 1, nk = ke - ks + 1;

  // (1) same-rank neighbours: a plain device copy, neighbour ACTIVE -> my GHOST
  par_for("m1_impl_halo_loc", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int mlo = nbr_.d_view(2*m);
    int mhi = nbr_.d_view(2*m+1);
    if (mlo >= 0) {iw_(m,nc,k,j,is-1-l) = iw_(mlo,nc,k,j,ie-l);}
    if (mhi >= 0) {iw_(m,nc,k,j,ie+1+l) = iw_(mhi,nc,k,j,is+l);}
  });
  if (!part_any_mpi) return;

#if MPI_PARALLEL_ENABLED
  // (2) remote neighbours.  Buffer slots: 0 = send to lo, 1 = send to hi,
  // 2 = recv from lo, 3 = recv from hi.
  auto hb_ = part_hbuf;
  par_for("m1_impl_halo_pack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    hb_(m,0,idx) = iw_(m,nc,k,j,is+l);
    hb_(m,1,idx) = iw_(m,nc,k,j,ie-l);
  });
  int nbuf = nq*nl*nk*nj;
  Kokkos::deep_copy(part_hbuf_h, part_hbuf);
  std::vector<MPI_Request> req;
  // the tag identifies (receiving local block, receiving side, which halo set); the
  // source rank is named in the receive, so it need only be unique per rank pair.
  int te = eponly ? 1 : 0;
  int *gr = pmy_pack->pmesh->gids_eachrank;
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_hbuf_h(m,2+s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*m + 2*s + te, MPI_COMM_WORLD, &req.back());
    }
  }
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      // the neighbour receives this message into ITS slot 2+(1-s)
      int lidn = part_nbrgid[2*m+s] - gr[rk];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_hbuf_h(m,s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*lidn + 2*(1-s) + te, MPI_COMM_WORLD, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  Kokkos::deep_copy(part_hbuf, part_hbuf_h);
  auto nbrk = part_nbr;
  par_for("m1_impl_halo_unpack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    if (nbrk.d_view(2*m) < 0) {iw_(m,nc,k,j,is-1-l) = hb_(m,2,idx);}
    if (nbrk.d_view(2*m+1) < 0) {iw_(m,nc,k,j,ie+1+l) = hb_(m,3,idx);}
  });
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransverseHalo
//! \brief milestone 3b phase B: put the M1_NHALO_T LAGGED quantities the x2/x3 faces (and
//! the lagged off-diagonal Eddington terms of every face) read at a neighbouring cell
//! into the scratch array `thw`, exchange it with ALL six neighbours through the module's
//! ORDINARY cell-centred boundary machinery -- which is what supplies periodic wrap,
//! edge/corner neighbours and MPI without another hand-rolled protocol -- and copy the
//! ghost zones back into `iw`.
//!
//! The values a block reads in its ghost zones are therefore the VERY NUMBERS its
//! neighbour computed, so the two blocks that share a face assemble that face's flux from
//! bit-identical inputs: the face flux is single-valued and the scheme is conservative
//! across a MeshBlock boundary (gate G4).
//!
//! PHYSICAL (non-periodic) boundaries are NOT filled here: every kernel below branches on
//! the boundary flag instead and imposes F = 0 on a physical x2/x3 face (reflecting).

void RadiationM1::ImplicitTransverseHalo(int nq) {
  if (!trans_on) return;
  ImplicitHaloExchange(nq, -1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloCopy
//! \brief copy `nq` components between iw and the scratch halo array `sc` (`topack`
//! picks the direction), over the HALO SHELL only.
//!
//! What the exchange reads out of `sc` is the outermost ng ACTIVE cells, and what it
//! writes back into it is the ghost zones (buffs_cc.cpp, isame); with same-level
//! neighbours only -- which is all the implicit solve allows -- nothing else in `sc` is
//! ever touched.  The box [is+ng,ie-ng] x [js+ng,je-ng] x [ks+ng,ke-ng] is therefore
//! neither sent nor received and need not be copied in either direction: on the way in
//! its value is never read, and on the way back it is a copy of what the pack put there.
//! At a 2-D 84 x 32 block with ng = 2 that is 480 cells per component instead of 3168.
//! `halo_shell` falls back to the full copy for the exchanges that reach deeper (the
//! cubed-sphere resample and the polar transform read the whole strip).
//!
//! The loop carries (m,n,k,j) and runs the contiguous i direction inside, so the index
//! arithmetic of the flattened range policy (an integer division per index) is paid once
//! per ROW instead of once per cell.

void RadiationM1::ImplicitHaloCopy(DvceArray5D<Real> &sc, int nq, int c0, bool topack) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // the deep interior, as an index box.  A direction with no neighbour (a degenerate
  // dimension) is interior everywhere; a direction too thin to have one (nx <= 2 ng), or
  // a run whose exchange reaches deeper, makes the box empty and the copy full.
  int ng = indcs.ng;
  int ilo = halo_shell ? (indcs.is + ng) : n1;
  int ihi = halo_shell ? (indcs.ie - ng) : (n1 - 1);
  int jlo = (indcs.nx2 > 1) ? (indcs.js + ng) : 0;
  int jhi = (indcs.nx2 > 1) ? (indcs.je - ng) : 0;
  int klo = (indcs.nx3 > 1) ? (indcs.ks + ng) : 0;
  int khi = (indcs.nx3 > 1) ? (indcs.ke - ng) : 0;
  if (ilo > ihi) {
    ilo = n1;
    ihi = n1 - 1;
  }
  const int il_ = ilo, iu_ = ihi, jl_ = jlo, ju_ = jhi, kl_ = klo, ku_ = khi;
  const int nc0 = c0, nn1 = n1;
  const bool pack_ = topack;
  auto iw_ = iw;
  auto sc_ = sc;
  if (impl_halo_direct) {
    // implicit_halo_direct on a mesh whose neighbours are not all on this rank: the
    // same shell copy with ONE thread per cell (coalesced along i) instead of one per
    // row -- the row loop runs ~13k threads with strided rows on the 3-D box.
    par_for("m1_impl_hcpyf", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_) && (i >= il_) &&
          (i <= iu_)) {
        return;
      }
      const int nc = (nc0 >= 0) ? nc0 : M1HaloCompT(n);
      if (pack_) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      } else {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    });
    return;
  }
  par_for("m1_impl_hcpy", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j) {
    const int nc = (nc0 >= 0) ? nc0 : M1HaloCompT(n);
    // the two i runs this row copies: the whole row unless the row is interior
    int ia = 0, ib = nn1 - 1, ic = nn1, id = nn1 - 1;
    if ((k >= kl_) && (k <= ku_) && (j >= jl_) && (j <= ju_)) {
      ib = il_ - 1;
      ic = iu_ + 1;
    }
    if (pack_) {
      for (int i=ia; i<=ib; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        sc_(m,n,k,j,i) = iw_(m,nc,k,j,i);
      }
    } else {
      for (int i=ia; i<=ib; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
      for (int i=ic; i<=id; ++i) {
        iw_(m,nc,k,j,i) = sc_(m,n,k,j,i);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloExchange
//! \brief the common core: pack, exchange through the ordinary cell-centred boundary
//! machinery on the PERSISTENT scratch array whose width matches `nq`, unpack.  The
//! three widths (M1_NHALO_T, M1_NHALO_Q, 1) have one array and one boundary object each;
//! they are used strictly sequentially, so they may share the MPI tag space.

void RadiationM1::ImplicitHaloExchange(int nq, int c0) {
  if (halo_direct_on) {   // implicit_halo_direct, every neighbour on this rank
    ImplicitHaloDirect(nq, c0);
    return;
  }
  if (impl_halo_mpi) {     // implicit_halo_mpi: on-rank copy + one message per rank
    if (hm_state == 0) {ImplicitHaloMPIInit();}
    if (hm_state == 1) {
      ImplicitHaloMPI(nq, c0);
      return;
    }
  }
  DvceArray5D<Real> *pa, *pc;
  MeshBoundaryValuesCC *pb;
  if (nq == M1_NHALO_T) {
    pa = &thw; pc = &thw_c; pb = pbval_th;
  } else if (nq == M1_NHALO_Q) {
    pa = &thq; pc = &thq_c; pb = pbval_tq;
  } else {
    pa = &krw; pc = &krw_c; pb = pbval_kr;
  }
  ImplicitHaloCopy(*pa, nq, c0, true);
  while (pb->InitRecv(nq) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(*pa, *pc) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
  ImplicitHaloCopy(*pa, nq, c0, false);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransTheta
//! \brief the TRANSVERSE realizability limiter: fill thx2/thx3, the ONE face theta every
//! use of the transverse face elimination reads (the face-flux kernels, the cell-terms
//! kernel, ImplicitOffDiagOp and the post-solve reconstruction).
//!
//! The face-eliminated transverse flux
//!   F_f' = theta_f [F_f^n - c^2 dt G_f - c dt v_f g0_f],  G_f = gr_f + off_f,
//!   theta_f = 1/(1 + c dt kt_f)
//! has NO free-streaming bound.  Its steady state is F_f = -c G_f/kt_f, the unlimited
//! diffusive flux, which at c dt/dx ~ 7e3 is super-luminal wherever rho kappa is tiny:
//! in the optically thin top of the 2-D He slab |F_2|/(c E) saturates at the post-solve
//! clip of 1 and E goes horizontally unphysical within 2 s.
//!
//! Under implicit_trans_limit = lp each face gets a LAGGED limiter opacity
//!   klim_f = |G_f| / (phi_f E_f),   E_f = (E_L + E_R)/2 of the lagged iterate,
//!   phi_f  = fmax sqrt(max(1 - f1_f^2, 0.01)),  f1_f the face mean of the lagged cell
//!            x1 reduced flux F1/(c E), clipped to [-1,1],
//! and theta_f = 1/(1 + c dt (kt_f + klim_f)).  Then in steady state
//!   |F_f| = c |G_f| / (kt_f + |G_f|/(phi_f E_f)) <= c phi_f E_f,
//! i.e. the transverse flux can never exceed the room the x1 flux leaves in the
//! realizability cone, while where R = |G_f|/(kt_f E_f) << 1 -- every optically thick
//! face -- the change is O(R) and the diffusion limit is untouched.  This is a flux
//! limiter of exactly the Levermore-Pomraning form, written as an opacity so that it
//! enters the ALREADY linear face elimination.
//!
//! klim only ever makes theta_f SMALLER, and theta_f multiplies the transverse
//! off-diagonals of the 7-point row (and its own diagonal contribution by the same
//! factor), so the row stays diagonally dominant with the same signs: the M-matrix
//! property that guarantees E' > 0 at any dt is preserved.
//!
//! klim follows implicit_closure_lag: with `step` it is computed ONCE per step, from the
//! entry state (newk is then true only on the first Picard pass), and frozen for every
//! later pass and every Krylov application of that step; with `pass` it is recomputed at
//! each pass, from that pass' lagged iterate.  Either way it is a lagged coefficient and
//! the linear system stays linear.

void RadiationM1::ImplicitTransTheta(bool newk) {
  if (impl_tlim == M1_TLIM_NONE) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto th2_ = thx2;
  auto kl2_ = klx2;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  auto th3_ = thx3;
  auto kl3_ = klx3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  Real cl = c_light, ch = chat, dt = dt_sub;
  Real efl = e_floor;
  Real fmx = impl_tfmax;
  const bool thrd = trans_x3;
  const bool nk = newk;
  const int odm = od_now;

  par_for("m1_impl_tlim2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      // a physical x2 face is reflecting: its flux is zero and its theta is never read
      kl2_(m,k,j,i) = 0.0;
      th2_(m,k,j,i) = 1.0;
      return;
    }
    int jm = j - 1;
    if (nk) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag a6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int kl = (!thrd || ((a5 != BoundaryFlag::block) &&
                          (a5 != BoundaryFlag::periodic))) ? ks : ks-1;
      int ku = (!thrd || ((a6 != BoundaryFlag::block) &&
                          (a6 != BoundaryFlag::periodic))) ? ke : ke+1;
      int jl = plo ? js : js-1;
      int ju = phi ? je : je+1;
      Real el = fmax(iw_(m,M1_IW_EP,k,jm,i), efl);
      Real er = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real dl = M1DDiag(iw_,vd_,dfull,m,1,k,jm,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
      Real gf = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/dx2;
      if (odm != M1_OD_NONE) {
        gf += 0.5*(M1OffDiv(iw_,m,1,k,jm,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      Real rl = iw_(m,M1_IW_F1,k,jm,i)/(cl*el);
      Real rr = iw_(m,M1_IW_F1,k,j,i)/(cl*er);
      rl = fmin(fmax(rl, -1.0), 1.0);
      rr = fmin(fmax(rr, -1.0), 1.0);
      Real f1f = fmin(fmax(0.5*(rl + rr), -1.0), 1.0);
      Real phif = fmx*sqrt(fmax(1.0 - f1f*f1f, 0.01));
      Real ef = fmax(0.5*(el + er), 1.0e-300);
      kl2_(m,k,j,i) = fabs(gf)/(phif*ef);
    }
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,jm,i) + iw_(m,M1_IW_KT,k,j,i));
    th2_(m,k,j,i) = 1.0/(1.0 + ch*dt*(ktf + kl2_(m,k,j,i)));
  });

  if (!thrd) return;
  par_for("m1_impl_tlim3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) {
      kl3_(m,k,j,i) = 0.0;
      th3_(m,k,j,i) = 1.0;
      return;
    }
    int km = k - 1;
    if (nk) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
      BoundaryFlag a4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      bool pj1 = (a3 != BoundaryFlag::block) && (a3 != BoundaryFlag::periodic);
      bool pj2 = (a4 != BoundaryFlag::block) && (a4 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int jl = pj1 ? js : js-1;
      int ju = pj2 ? je : je+1;
      int kl = plo ? ks : ks-1;
      int ku = phi ? ke : ke+1;
      Real el = fmax(iw_(m,M1_IW_EP,km,j,i), efl);
      Real er = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real dl = M1DDiag(iw_,vd_,dfull,m,2,km,j,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      Real gf = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/dx3;
      if (odm != M1_OD_NONE) {
        gf += 0.5*(M1OffDiv(iw_,m,2,km,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      Real rl = iw_(m,M1_IW_F1,km,j,i)/(cl*el);
      Real rr = iw_(m,M1_IW_F1,k,j,i)/(cl*er);
      rl = fmin(fmax(rl, -1.0), 1.0);
      rr = fmin(fmax(rr, -1.0), 1.0);
      Real f1f = fmin(fmax(0.5*(rl + rr), -1.0), 1.0);
      Real phif = fmx*sqrt(fmax(1.0 - f1f*f1f, 0.01));
      Real ef = fmax(0.5*(el + er), 1.0e-300);
      kl3_(m,k,j,i) = fabs(gf)/(phif*ef);
    }
    Real ktf = 0.5*(iw_(m,M1_IW_KT,km,j,i) + iw_(m,M1_IW_KT,k,j,i));
    th3_(m,k,j,i) = 1.0/(1.0 + ch*dt*(ktf + kl3_(m,k,j,i)));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTransverseTerms
//! \brief milestone 3b phase B: one line-Jacobi evaluation of the TRANSVERSE (x2, x3)
//! part of the operator at the current Picard iterate.
//!
//! Per direction d and face f, exactly the face-eliminated form the x1 faces use:
//!
//!   F0_f' = theta_f [ F0_f^n - c^2 dt (P_dd,R - P_dd,L)/dx_d - c dt v_f g0_f
//!                     - c^2 dt (sum_{e != d} d_e P_de)_f ]
//!   theta_f = 1/(1 + c dt (rho kappa_t)_f),  (rho kappa_t)_f the arithmetic face mean
//!
//! with P_dd = D_dd E of the LAGGED closure, D_dd = (1-chi)/2 + (3 chi - 1)/2 n_d^2, and
//! the off-diagonal divergence fully lagged (centred differences of the previous pass' E
//! and closure), so it is a pure right-hand-side term.  The advective enthalpy flux
//! A_d = v_d E + (v.P)_d is upwinded with the face velocity, exactly as in x1.
//!
//! The result is split into
//!   M1_IW_TDIA = dT/dE_c >= 0, which goes on the matrix DIAGONAL, and
//!   M1_IW_TRHS = -(T(E^k) - TDIA E^k_c), the neighbours' lagged contribution.
//! Keeping the diagonal part on the matrix is what preserves the full 7-point M-matrix:
//! the column sums stay >= 1 and E' > 0 at any dt.
//!
//! It also measures the TRUE residual of the full 7-point linear system.  After the line
//! solve E^{k+1} satisfies  D1(E^{k+1}) + TDIA E^{k+1} + U(E^k) = b  with
//! U(E) = T(E) - TDIA E_c, so the residual of the FULL system at E^{k+1} is exactly
//! U(E^k) - U(E^{k+1}): the change of the lagged off-diagonal term between two passes.
//! That is what M1_IW_LRES holds, and it is tested separately from the Picard residual --
//! a lagged coupling can stall and look converged on |dE|/E alone.

void RadiationM1::ImplicitTransverseTerms(bool first) {
  if (!trans_on) return;
  // the transverse realizability limiter: one evaluation of theta for the whole pass,
  // which everything below and the Krylov operator then READ.  klim itself follows the
  // closure lag (frozen for the step under implicit_closure_lag = step).
  ImplicitTransTheta(!impl_clag_step || first);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto f2_ = f0x2;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  auto f2n_ = f0x2n;
  auto f3_ = f0x3;
  auto f3n_ = f0x3n;
  // the transverse realizability limiter.  Under `none` the OLD expression for theta is
  // taken, term for term, so the arithmetic of every earlier configuration is bitwise
  // unchanged; under `lp` every use below reads the SAME thx2/thx3 the Krylov operator
  // reads.
  const bool lm = (impl_tlim != M1_TLIM_NONE);
  auto th2_ = thx2;
  auto th3_ = thx3;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  Real cl = c_light, ch = chat, dt = dt_sub;
  const bool thrd = trans_x3;
  const bool fst = first;
  const Real wmem = dbg_trans_memory;
  // MILESTONE 3b phase C: also store the transverse OFF-DIAGONAL coefficients of the
  // frozen row, which is what the BiCGStab operator applies and what turns TRHS back into
  // the right-hand side of the full system.
  const bool bcg = bicg_on;
  // MILESTONE 3b phase D: implicit_offdiag.  The face flux stored here is the PHYSICAL
  // one and keeps the off-diagonal term at the current iterate under `lagged` and
  // `operator` alike; what changes between the two is the row the solve is given (see
  // ImplicitSolve step (e), where the term is subtracted from the right-hand side and
  // handed to the operator instead).  Under `none` it is dropped everywhere.
  const int odm = od_now;

  // (1) the x2 face fluxes
  par_for("m1_impl_f2face", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      f2_(m,k,j,i) = 0.0;
      return;
    }
    Real dx1 = mbsize.d_view(m).dx1;
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
    BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
    BoundaryFlag a5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag a6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
    bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
    int il = pi1 ? is : is-1;
    int iu = pi2 ? ie : ie+1;
    int kl = (!thrd || ((a5 != BoundaryFlag::block) &&
                        (a5 != BoundaryFlag::periodic))) ? ks : ks-1;
    int ku = (!thrd || ((a6 != BoundaryFlag::block) &&
                        (a6 != BoundaryFlag::periodic))) ? ke : ke+1;
    int jl = plo ? js : js-1;
    int ju = phi ? je : je+1;
    int jm = j - 1;
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,jm,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
    Real dl = M1DDiag(iw_,vd_,dfull,m,1,k,jm,i);
    Real dr = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
    Real gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,k,jm,i))/dx2;
    Real vf = 0.5*(iw_(m,M1_IW_V2,k,jm,i) + iw_(m,M1_IW_V2,k,j,i));
    Real g0f = 0.5*(iw_(m,M1_IW_G0,k,jm,i) + iw_(m,M1_IW_G0,k,j,i));
    Real off = 0.0;
    if (odm != M1_OD_NONE) {
      off = 0.5*(M1OffDiv(iw_,m,1,k,jm,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                          dfull)
                 + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                            dfull));
    }
    f2_(m,k,j,i) = th*(wmem*f2n_(m,k,j,i) - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*off);
  });

  // (2) the x3 face fluxes
  if (thrd) {
    par_for("m1_impl_f3face", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
      bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
      if ((k == ks && plo) || (k == ke+1 && phi)) {
        f3_(m,k,j,i) = 0.0;
        return;
      }
      Real dx1 = mbsize.d_view(m).dx1;
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag a1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
      BoundaryFlag a2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
      BoundaryFlag a3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
      BoundaryFlag a4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
      bool pi1 = (a1 != BoundaryFlag::block) && (a1 != BoundaryFlag::periodic);
      bool pi2 = (a2 != BoundaryFlag::block) && (a2 != BoundaryFlag::periodic);
      bool pj1 = (a3 != BoundaryFlag::block) && (a3 != BoundaryFlag::periodic);
      bool pj2 = (a4 != BoundaryFlag::block) && (a4 != BoundaryFlag::periodic);
      int il = pi1 ? is : is-1;
      int iu = pi2 ? ie : ie+1;
      int jl = pj1 ? js : js-1;
      int ju = pj2 ? je : je+1;
      int kl = plo ? ks : ks-1;
      int ku = phi ? ke : ke+1;
      int km = k - 1;
      Real ktf = 0.5*(iw_(m,M1_IW_KT,km,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real dl = M1DDiag(iw_,vd_,dfull,m,2,km,j,i);
      Real dr = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      Real gr = (dr*iw_(m,M1_IW_EP,k,j,i) - dl*iw_(m,M1_IW_EP,km,j,i))/dx3;
      Real vf = 0.5*(iw_(m,M1_IW_V3,km,j,i) + iw_(m,M1_IW_V3,k,j,i));
      Real g0f = 0.5*(iw_(m,M1_IW_G0,km,j,i) + iw_(m,M1_IW_G0,k,j,i));
      Real off = 0.0;
      if (odm != M1_OD_NONE) {
        off = 0.5*(M1OffDiv(iw_,m,2,km,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                            vd_,dfull)
                   + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                              vd_,dfull));
      }
      f3_(m,k,j,i) = th*(wmem*f3n_(m,k,j,i)
                         - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*off);
    });
  }

  // (3) the cell terms: the diagonal part, the lagged right-hand side and the residual
  // of the full 7-point system.  implicit_enthalpy adds its deferred correction to the
  // face fluxes fp/fm/gp/gm, which reach only the right-hand side (TRHS), never the row.
  const int enm = impl_enth;
  const bool enth2 = (enm != M1_IENTH_UPWIND);
  const int ngh = indcs.ng;
  par_for("m1_impl_tcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dx2 = mbsize.d_view(m).dx2;
    Real dx3 = mbsize.d_view(m).dx3;
    Real cr = ch/cl;
    Real ec = iw_(m,M1_IW_EP,k,j,i);
    Real dia = 0.0, tt = 0.0;
    BoundaryFlag b3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag b4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool p2lo = (b3 != BoundaryFlag::block) && (b3 != BoundaryFlag::periodic);
    bool p2hi = (b4 != BoundaryFlag::block) && (b4 != BoundaryFlag::periodic);
    Real nu2 = dt/dx2;
    Real d2c = M1DDiag(iw_,vd_,dfull,m,1,k,j,i);
    Real a2c = iw_(m,M1_IW_A2,k,j,i);
    Real fp = 0.0, fm = 0.0;
    Real cjp = 0.0, cjm = 0.0;
    if (!(j == je && p2hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
      Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V2,k,j,i) + iw_(m,M1_IW_V2,k,j+1,i));
      fp = f2_(m,k,j+1,i);
      if (vf > 0.0) {
        fp += a2c*ec;
        dia += nu2*cr*a2c;
      } else {
        fp += iw_(m,M1_IW_A2,k,j+1,i)*iw_(m,M1_IW_EP,k,j+1,i);
        if (bcg) {cjp += nu2*cr*iw_(m,M1_IW_A2,k,j+1,i);}
      }
      if (enth2) {
        bool o0, o3;
        int j0 = M1EnthIdx(j-1, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o0);
        int j3 = M1EnthIdx(j+2, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o3);
        fp += M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j0,i), ec, iw_(m,M1_IW_EP,k,j+1,i),
                         iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                         a2c, iw_(m,M1_IW_A2,k,j+1,i), iw_(m,M1_IW_A2,k,j3,i), vf);
      }
      dia += nu2*th*ch*ch*dt*d2c/dx2;
      if (bcg) {
        Real d2p = M1DDiag(iw_,vd_,dfull,m,1,k,j+1,i);
        cjp -= nu2*th*ch*ch*dt*d2p/dx2;
      }
    }
    if (!(j == js && p2lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real vf = 0.5*(iw_(m,M1_IW_V2,k,j-1,i) + iw_(m,M1_IW_V2,k,j,i));
      fm = f2_(m,k,j,i);
      if (vf > 0.0) {
        fm += iw_(m,M1_IW_A2,k,j-1,i)*iw_(m,M1_IW_EP,k,j-1,i);
        if (bcg) {cjm -= nu2*cr*iw_(m,M1_IW_A2,k,j-1,i);}
      } else {
        fm += a2c*ec;
        dia -= nu2*cr*a2c;
      }
      if (enth2) {
        bool o0, o3;
        int j0 = M1EnthIdx(j-2, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o0);
        int j3 = M1EnthIdx(j+1, js, je, false, p2lo ? 0 : ngh, p2hi ? 0 : ngh, o3);
        fm += M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j0,i), iw_(m,M1_IW_EP,k,j-1,i), ec,
                         iw_(m,M1_IW_EP,k,j3,i), o0 && o3, iw_(m,M1_IW_A2,k,j0,i),
                         iw_(m,M1_IW_A2,k,j-1,i), a2c, iw_(m,M1_IW_A2,k,j3,i), vf);
      }
      dia += nu2*th*ch*ch*dt*d2c/dx2;
      if (bcg) {
        Real d2m = M1DDiag(iw_,vd_,dfull,m,1,k,j-1,i);
        cjm -= nu2*th*ch*ch*dt*d2m/dx2;
      }
    }
    tt += nu2*cr*(fp - fm);
    if (bcg) {
      iw_(m,M1_IW_CJM,k,j,i) = cjm;
      iw_(m,M1_IW_CJP,k,j,i) = cjp;
      iw_(m,M1_IW_CKM,k,j,i) = 0.0;
      iw_(m,M1_IW_CKP,k,j,i) = 0.0;
    }
    if (thrd) {
      BoundaryFlag b5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag b6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool p3lo = (b5 != BoundaryFlag::block) && (b5 != BoundaryFlag::periodic);
      bool p3hi = (b6 != BoundaryFlag::block) && (b6 != BoundaryFlag::periodic);
      Real nu3 = dt/dx3;
      Real d3c = M1DDiag(iw_,vd_,dfull,m,2,k,j,i);
      Real a3c = iw_(m,M1_IW_A3,k,j,i);
      Real gp = 0.0, gm = 0.0;
      Real ckp = 0.0, ckm = 0.0;
      if (!(k == ke && p3hi)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
        Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V3,k,j,i) + iw_(m,M1_IW_V3,k+1,j,i));
        gp = f3_(m,k+1,j,i);
        if (vf > 0.0) {
          gp += a3c*ec;
          dia += nu3*cr*a3c;
        } else {
          gp += iw_(m,M1_IW_A3,k+1,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
          if (bcg) {ckp += nu3*cr*iw_(m,M1_IW_A3,k+1,j,i);}
        }
        if (enth2) {
          bool o0, o3;
          int k0 = M1EnthIdx(k-1, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o0);
          int k3 = M1EnthIdx(k+2, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o3);
          gp += M1EnthCorr(enm, iw_(m,M1_IW_EP,k0,j,i), ec, iw_(m,M1_IW_EP,k+1,j,i),
                           iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                           iw_(m,M1_IW_A3,k0,j,i), a3c, iw_(m,M1_IW_A3,k+1,j,i),
                           iw_(m,M1_IW_A3,k3,j,i), vf);
        }
        dia += nu3*th*ch*ch*dt*d3c/dx3;
        if (bcg) {
          Real d3p = M1DDiag(iw_,vd_,dfull,m,2,k+1,j,i);
          ckp -= nu3*th*ch*ch*dt*d3p/dx3;
        }
      }
      if (!(k == ks && p3lo)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
        Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V3,k-1,j,i) + iw_(m,M1_IW_V3,k,j,i));
        gm = f3_(m,k,j,i);
        if (vf > 0.0) {
          gm += iw_(m,M1_IW_A3,k-1,j,i)*iw_(m,M1_IW_EP,k-1,j,i);
          if (bcg) {ckm -= nu3*cr*iw_(m,M1_IW_A3,k-1,j,i);}
        } else {
          gm += a3c*ec;
          dia -= nu3*cr*a3c;
        }
        if (enth2) {
          bool o0, o3;
          int k0 = M1EnthIdx(k-2, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o0);
          int k3 = M1EnthIdx(k+1, ks, ke, false, p3lo ? 0 : ngh, p3hi ? 0 : ngh, o3);
          gm += M1EnthCorr(enm, iw_(m,M1_IW_EP,k0,j,i), iw_(m,M1_IW_EP,k-1,j,i), ec,
                           iw_(m,M1_IW_EP,k3,j,i), o0 && o3,
                           iw_(m,M1_IW_A3,k0,j,i), iw_(m,M1_IW_A3,k-1,j,i), a3c,
                           iw_(m,M1_IW_A3,k3,j,i), vf);
        }
        dia += nu3*th*ch*ch*dt*d3c/dx3;
        if (bcg) {
          Real d3m = M1DDiag(iw_,vd_,dfull,m,2,k-1,j,i);
          ckm -= nu3*th*ch*ch*dt*d3m/dx3;
        }
      }
      tt += nu3*cr*(gp - gm);
      if (bcg) {
        iw_(m,M1_IW_CKM,k,j,i) = ckm;
        iw_(m,M1_IW_CKP,k,j,i) = ckp;
      }
    }
    Real unew = tt - dia*ec;
    Real uold = -iw_(m,M1_IW_TRHS,k,j,i);
    iw_(m,M1_IW_LRES,k,j,i) = fst ? 1.0e30 : fabs(unew - uold);
    iw_(m,M1_IW_TDIA,k,j,i) = dia;
    iw_(m,M1_IW_TRHS,k,j,i) = -unew;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGatherSolve
//! \brief gather the assembled rows (a,b,c,r) of every x1 stack onto its root block, run
//! there the IDENTICAL serial Thomas sweep the single-block solve runs, and scatter the
//! solution back into iw(M1_IW_S2).  Two messages per Picard iteration per non-root
//! block (a gather of 4*nx1 reals per column and a scatter of nx1).
//!
//! Cost: the root sweeps part_nblk*nx1 rows serially per column.  That is acceptable
//! here (the columns are independent and the kernel is parallel over (root,k,j), and
//! these columns are 84-512 cells long), but it is the serial bottleneck of the scheme;
//! the scalable successor is Schur condensation to the block-interface unknowns (one
//! reduced tridiagonal system of part_nblk rows per column, solved after two local
//! sweeps) or cyclic reduction, neither of which can be bitwise.

void RadiationM1::ImplicitGatherSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nx1 = indcs.nx1, nx1g = part_nx1g;
  int nj = je - js + 1, nk = ke - ks + 1;
  auto iw_ = iw;
  auto sys_ = part_sys;
  auto pos_ = part_pos;
  auto slot_ = part_slot;

  // (1) same-rank members: copy their rows straight into the root's gathered system
  par_for("m1_impl_gth_loc", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    sys_(sl,n,k,j,pos_.d_view(m)*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
  });

#if MPI_PARALLEL_ENABLED
  int nrow = 4*nk*nj*nx1;
  std::vector<MPI_Request> req;
  if (part_any_mpi) {
    auto sb_ = part_sbuf;
    par_for("m1_impl_gth_pack", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      sb_(m,((n*nk + (k-ks))*nj + (j-js))*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
    });
    Kokkos::deep_copy(part_sbuf_h, part_sbuf);
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(&part_rbuf_h(s*part_nblk+p,0), nrow, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]), MPI_COMM_WORLD, &req.back());
      }
    }
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_sbuf_h(m,0), nrow, MPI_ATHENA_REAL, part_rootrank[m],
                2*m, MPI_COMM_WORLD, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    req.clear();
    if (part_nroot > 0) {
      Kokkos::deep_copy(part_rbuf, part_rbuf_h);
      auto rb_ = part_rbuf;
      // which (slot,position) pairs are remote: encoded as a host loop over kernels
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_gth_unp", DevExeSpace(), 0, 3, ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int n, const int k, const int j, const int i) {
            sys_(ss,n,k,j,pp*nx1 + i) =
                rb_(row,((n*nk + (k-ks))*nj + (j-js))*nx1 + i);
          });
        }
      }
    }
  }
#endif

  // (2) the gathered Thomas sweep: the same arithmetic, over nx1g rows
  if (part_nroot > 0) {
    par_for("m1_impl_gth_thomas", DevExeSpace(), 0, part_nroot-1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int s, const int k, const int j) {
      Real bet = sys_(s,1,k,j,0);
      sys_(s,5,k,j,0) = sys_(s,3,k,j,0)/bet;
      for (int i=1; i<nx1g; ++i) {
        sys_(s,4,k,j,i) = sys_(s,2,k,j,i-1)/bet;
        bet = sys_(s,1,k,j,i) - sys_(s,0,k,j,i)*sys_(s,4,k,j,i);
        sys_(s,5,k,j,i) = (sys_(s,3,k,j,i) - sys_(s,0,k,j,i)*sys_(s,5,k,j,i-1))/bet;
      }
      for (int i=nx1g-2; i>=0; --i) {
        sys_(s,5,k,j,i) -= sys_(s,4,k,j,i+1)*sys_(s,5,k,j,i+1);
      }
    });
  }

  // (3) scatter
#if MPI_PARALLEL_ENABLED
  if (part_any_mpi) {
    int nsol = nk*nj*nx1;
    if (part_nroot > 0) {
      auto rb_ = part_rbuf;
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_sct_pack", DevExeSpace(), ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int k, const int j, const int i) {
            rb_(row,((k-ks)*nj + (j-js))*nx1 + i) = sys_(ss,5,k,j,pp*nx1 + i);
          });
        }
      }
      Kokkos::deep_copy(part_rbuf_h, part_rbuf);
    }
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_sbuf_h(m,0), nsol, MPI_ATHENA_REAL, part_rootrank[m],
                2*m + 1, MPI_COMM_WORLD, &req.back());
    }
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&part_rbuf_h(s*part_nblk+p,0), nsol, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]) + 1, MPI_COMM_WORLD,
                  &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    Kokkos::deep_copy(part_sbuf, part_sbuf_h);
    auto sb_ = part_sbuf;
    par_for("m1_impl_sct_unp", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      iw_(m,M1_IW_S2,k,j,i) = sb_(m,((k-ks)*nj + (j-js))*(ie-is+1) + (i-is));
    });
  }
#endif
  par_for("m1_impl_sct_loc", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    iw_(m,M1_IW_S2,k,j,i) = sys_(sl,5,k,j,pos_.d_view(m)*(ie-is+1) + (i-is));
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitTridiagSolve
//! \brief the x1 LINE SOLVE of the assembled rows: read (M1_IW_TA, TB, TC, TR) and leave
//! the solution in M1_IW_S2.  Plain Thomas, cyclic Thomas (Sherman-Morrison) when the x1
//! boundaries are periodic inside one MeshBlock, or the gathered stack sweep when the
//! column spans several blocks.
//!
//! It is BOTH the line-Jacobi pass (called once per Picard iteration with the lagged
//! right-hand side) and the PRECONDITIONER of the BiCGStab wrapper (called twice per
//! inner iteration with a Krylov vector in M1_IW_TR).  Extracting it changed no
//! arithmetic: the two sweeps below are verbatim what the Picard loop used to run inline.

void RadiationM1::ImplicitTridiagSolve() {
  if (part_nblk > 1) {
    ImplicitGatherSolve();
    return;
  }
  if (!impl_pcr_check) {
    if (impl_line_solver == 1) {
      ImplicitPCRSolve();
    } else {
      ImplicitThomasSolve();
    }
    return;
  }
  // implicit_pcr_check: run the OTHER solver first, keep its answer, then the selected
  // one (whose answer stays in M1_IW_S2), and record max|dx|/max|x| over the pack.
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  if (pcr_chk.extent(0) == 0) {
    Kokkos::realloc(pcr_chk, nmb1+1, indcs.nx3 + 2*indcs.ng*(indcs.nx3 > 1 ? 1 : 0),
                    indcs.nx2 + 2*indcs.ng*(indcs.nx2 > 1 ? 1 : 0), indcs.nx1+2*indcs.ng);
  }
  auto ck_ = pcr_chk;
  if (impl_line_solver == 1) {
    ImplicitThomasSolve();
  } else {
    ImplicitPCRSolve();
  }
  par_for("m1_impl_pcrck_cp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    ck_(m,k,j,i) = iw_(m,M1_IW_S2,k,j,i);
  });
  if (impl_line_solver == 1) {
    ImplicitPCRSolve();
  } else {
    ImplicitThomasSolve();
  }
  const int nkji = (ke-ks+1)*(je-js+1)*(ie-is+1);
  const int nji = (je-js+1)*(ie-is+1), ni = ie-is+1;
  Real dmax = 0.0, xmax = 0.0;
  Kokkos::parallel_reduce("m1_impl_pcrck_red",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dm, Real &xm) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + is;
    k += ks;
    j += js;
    dm = fmax(dm, fabs(iw_(m,M1_IW_S2,k,j,i) - ck_(m,k,j,i)));
    xm = fmax(xm, fabs(ck_(m,k,j,i)));
  }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(xmax));
  Real rel = (xmax > 0.0) ? dmax/xmax : dmax;
  pcr_chk_max = fmax(pcr_chk_max, rel);
  pcr_chk_n += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitThomasSolve
//! \brief the Thomas / cyclic-Thomas line solve, one thread per (m,k,j) column (the
//! original ImplicitTridiagSolve body, unchanged)

void RadiationM1::ImplicitThomasSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  par_for("m1_impl_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (!cyclic) {
      Real bet = iw_(m,M1_IW_TB,k,j,is);
      iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
      for (int i=is+1; i<=ie; ++i) {
        iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
        bet = iw_(m,M1_IW_TB,k,j,i) - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
        iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
      }
      for (int i=ie-1; i>=is; --i) {
        iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
      }
    } else {
      // cyclic tridiagonal, Sherman-Morrison (Press et al. `cyclic`).  alpha is the
      // BOTTOM-LEFT corner, c(ie) (row ie coupling to cell is), and beta the TOP-RIGHT
      // one, a(is) (row is coupling to cell ie) -- that is the convention u and v below
      // are built for, and swapping them is invisible on a symmetric matrix but wrong
      // as soon as upwind advection makes the two off-diagonals differ: it produced a
      // spurious dipole across the seam (-14 % in the first cell, +11 % in the last)
      // on the very first step of T4b, where the exact answer is "nothing moves".
      Real alpha = iw_(m,M1_IW_TC,k,j,ie);
      Real beta = iw_(m,M1_IW_TA,k,j,is);
      Real gam = -iw_(m,M1_IW_TB,k,j,is);
      Real bb0 = iw_(m,M1_IW_TB,k,j,is) - gam;
      Real bbn = iw_(m,M1_IW_TB,k,j,ie) - alpha*beta/gam;
      // solve A' y = r and A' z = u with u = (gam,0,...,0,alpha)
      Real bet = bb0;
      iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
      iw_(m,M1_IW_S3,k,j,is) = gam/bet;
      for (int i=is+1; i<=ie; ++i) {
        Real bd = (i == ie) ? bbn : iw_(m,M1_IW_TB,k,j,i);
        iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
        bet = bd - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
        iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
        Real uu = (i == ie) ? alpha : 0.0;
        iw_(m,M1_IW_S3,k,j,i) = (uu
                                 - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S3,k,j,i-1))/bet;
      }
      for (int i=ie-1; i>=is; --i) {
        iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
        iw_(m,M1_IW_S3,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S3,k,j,i+1);
      }
      // x = y - z (v.y)/(1 + v.z),  v = (1,0,...,0,beta/gam)
      Real vy = iw_(m,M1_IW_S2,k,j,is) + (beta/gam)*iw_(m,M1_IW_S2,k,j,ie);
      Real vz = iw_(m,M1_IW_S3,k,j,is) + (beta/gam)*iw_(m,M1_IW_S3,k,j,ie);
      Real fac = vy/(1.0 + vz);
      for (int i=is; i<=ie; ++i) {
        iw_(m,M1_IW_S2,k,j,i) -= fac*iw_(m,M1_IW_S3,k,j,i);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPCRSolve
//! \brief implicit_line_solver = pcr: the same x1 line solve as the Thomas sweep of
//! ImplicitTridiagSolve (non-cyclic, or cyclic by the same Sherman-Morrison split), by
//! PARALLEL CYCLIC REDUCTION with one Kokkos team per (m,k,j) column.  The rows are
//! loaded into team scratch with the team's threads running along i (coalesced), then
//! ceil(log2 nx1) PCR rounds (double-buffered, one team barrier each) decouple every
//! row, and x_i = r_i/b_i.  Out-of-range neighbours of a round are the identity row.
//! Under the cyclic split the second right-hand side u = (gam,0,...,0,alpha) rides the
//! same elimination.  Writes only M1_IW_S2 (the Thomas scratch S1/S3 is not touched and
//! nothing else reads it).  O(n log n) work, so on a CPU (team size 1) it is slower than
//! Thomas; it is meant for the GPU, where the Thomas sweep has one thread per column.

void RadiationM1::ImplicitPCRSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = ie - is + 1;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const int nkj = nk*nj;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const int nv = cyclic ? 5 : 4;   // a, b, c, r (+ u) per buffer
  size_t scr_size = ScrArray1D<Real>::shmem_size(2*nv*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  Kokkos::TeamPolicy<DevExeSpace> policy;
  if (!std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    // a GPU: an explicit team size (the threads run along i)
    int ts = impl_pcr_team;
    if (ts == 0) {
      ts = 1;
      while (ts < nx && ts < 256) ts *= 2;
    }
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, ts);
  } else {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, Kokkos::AUTO);
  }
  Kokkos::parallel_for("m1_impl_pcr",
                       policy.set_scratch_size(0, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/nkj;
    const int k = (tm.league_rank() - m*nkj)/nj + ks;
    const int j = (tm.league_rank() - m*nkj)%nj + js;
    ScrArray1D<Real> sw(tm.team_scratch(0), 2*nv*nx);
    // buffer q (0/1), variable v: sw(q*nv*nx + v*nx + i); v = 0 a, 1 b, 2 c, 3 r, 4 u
    Real alpha = 0.0, beta = 0.0, gam = 1.0;
    if (cyclic) {
      alpha = iw_(m,M1_IW_TC,k,j,ie);
      beta = iw_(m,M1_IW_TA,k,j,is);
      gam = -iw_(m,M1_IW_TB,k,j,is);
    }
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      const int ii = i + is;
      Real bd = iw_(m,M1_IW_TB,k,j,ii);
      if (cyclic) {
        if (i == 0) bd -= gam;
        if (i == nx-1) bd -= alpha*beta/gam;
        sw(4*nx + i) = (i == 0) ? gam : ((i == nx-1) ? alpha : 0.0);
      }
      sw(i) = (i == 0) ? 0.0 : iw_(m,M1_IW_TA,k,j,ii);
      sw(nx + i) = bd;
      sw(2*nx + i) = (i == nx-1) ? 0.0 : iw_(m,M1_IW_TC,k,j,ii);
      sw(3*nx + i) = iw_(m,M1_IW_TR,k,j,ii);
    });
    tm.team_barrier();
    int src = 0;
    for (int rd=0, s=1; rd<nround; ++rd, s*=2) {
      const int o = src*nv*nx, d = (1-src)*nv*nx;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        const int im = i - s, ip = i + s;
        Real ai = sw(o + i), bi = sw(o + nx + i), ci = sw(o + 2*nx + i);
        Real ri = sw(o + 3*nx + i);
        Real ui = cyclic ? sw(o + 4*nx + i) : 0.0;
        Real an = 0.0, cn = 0.0;
        if (im >= 0) {
          Real f = -ai/sw(o + nx + im);
          an = f*sw(o + im);
          bi += f*sw(o + 2*nx + im);
          ri += f*sw(o + 3*nx + im);
          if (cyclic) ui += f*sw(o + 4*nx + im);
        }
        if (ip < nx) {
          Real g = -ci/sw(o + nx + ip);
          cn = g*sw(o + 2*nx + ip);
          bi += g*sw(o + ip);
          ri += g*sw(o + 3*nx + ip);
          if (cyclic) ui += g*sw(o + 4*nx + ip);
        }
        sw(d + i) = an;
        sw(d + nx + i) = bi;
        sw(d + 2*nx + i) = cn;
        sw(d + 3*nx + i) = ri;
        if (cyclic) sw(d + 4*nx + i) = ui;
      });
      tm.team_barrier();
      src = 1 - src;
    }
    const int o = src*nv*nx;
    if (!cyclic) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        iw_(m,M1_IW_S2,k,j,i+is) = sw(o + 3*nx + i)/sw(o + nx + i);
      });
    } else {
      // x = y - z (v.y)/(1 + v.z),  v = (1,0,...,0,beta/gam)
      Real y0 = sw(o + 3*nx)/sw(o + nx);
      Real yn = sw(o + 4*nx - 1)/sw(o + 2*nx - 1);
      Real z0 = sw(o + 4*nx)/sw(o + nx);
      Real zn = sw(o + 5*nx - 1)/sw(o + 2*nx - 1);
      Real fac = (y0 + (beta/gam)*yn)/(1.0 + z0 + (beta/gam)*zn);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        Real bi = sw(o + nx + i);
        iw_(m,M1_IW_S2,k,j,i+is) = sw(o + 3*nx + i)/bi - fac*(sw(o + 4*nx + i)/bi);
      });
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPCRSolveX
//! \brief implicit_krylov_fuse >= 1: the pcr line solve of ImplicitPCRSolve, term for
//! term, with the right-hand side taken from component `rc` of iw (upd = 0) or MADE in
//! the load phase as the BiCGStab update p = r + c1 (p - c2 v) (upd = 1, written to
//! M1_IW_KP) or s = r - c1 v (upd = 2, written to M1_IW_KS), and the answer written
//! straight into component `zc`.  That removes the staging copy into M1_IW_TR, the copy
//! out of M1_IW_S2 and the separate p / s kernels: 1 launch where there were 3.
//!
//! implicit_precond = rbgs: `col` = 0 / 1 solves only the columns (k,j) whose parity
//! (k-ks)+(j-js) is `col` (one launch per colour, half the teams), and `sub` >= 0
//! subtracts the transverse 5-point coupling sum_nb C_nb z_nb of component `sub` from
//! the right-hand side, over the neighbours INSIDE the MeshBlock only (the other colour,
//! which the previous half-sweep has just written).  col < 0: every column, sub ignored.

namespace {
//! the kernel of ImplicitPCRSolveX, with the elimination carried in T (Real, or float
//! under implicit_precond_float: a preconditioner only has to be a FIXED linear map,
//! so its precision sets the convergence rate, not the converged answer)
template <typename T>
void M1PCRX(const DvceArray5D<Real> &iw_, Kokkos::TeamPolicy<DevExeSpace> policy,
            const int is, const int ie, const int js, const int je, const int ks,
            const int ke, const int nkj, const int njl, const bool colr, const int cl_,
            const int cs_, const bool thrd, const bool cyclic, const int cr,
            const int cz, const int up, const Real a1, const Real a2) {
  const int nx = ie - is + 1;
  const int nv = cyclic ? 5 : 4;   // a, b, c, r (+ u) per buffer
  size_t scr_size = ScrArray1D<T>::shmem_size(2*nv*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  Kokkos::parallel_for("m1_impl_pcrx",
                       policy.set_scratch_size(0, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/nkj;
    const int kk = (tm.league_rank() - m*nkj)/njl;
    const int jj = (tm.league_rank() - m*nkj)%njl;
    const int k = kk + ks;
    const int j = colr ? (js + 2*jj + ((cl_ + kk) & 1)) : (jj + js);
    if (j > je) return;   // team-uniform
    ScrArray1D<T> sw(tm.team_scratch(0), 2*nv*nx);
    Real alpha = 0.0, beta = 0.0, gam = 1.0;
    if (cyclic) {
      alpha = iw_(m,M1_IW_TC,k,j,ie);
      beta = iw_(m,M1_IW_TA,k,j,is);
      gam = -iw_(m,M1_IW_TB,k,j,is);
    }
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      const int ii = i + is;
      Real bd = iw_(m,M1_IW_TB,k,j,ii);
      if (cyclic) {
        if (i == 0) bd -= gam;
        if (i == nx-1) bd -= alpha*beta/gam;
        sw(4*nx + i) = static_cast<T>((i == 0) ? gam : ((i == nx-1) ? alpha : 0.0));
      }
      sw(i) = static_cast<T>((i == 0) ? 0.0 : iw_(m,M1_IW_TA,k,j,ii));
      sw(nx + i) = static_cast<T>(bd);
      sw(2*nx + i) = static_cast<T>((i == nx-1) ? 0.0 : iw_(m,M1_IW_TC,k,j,ii));
      Real rr;
      if (up == 1) {
        rr = iw_(m,M1_IW_KR,k,j,ii)
             + a1*(iw_(m,M1_IW_KP,k,j,ii) - a2*iw_(m,M1_IW_KV,k,j,ii));
        iw_(m,M1_IW_KP,k,j,ii) = rr;
      } else if (up == 2) {
        rr = iw_(m,M1_IW_KR,k,j,ii) - a1*iw_(m,M1_IW_KV,k,j,ii);
        iw_(m,M1_IW_KS,k,j,ii) = rr;
      } else {
        rr = iw_(m,cr,k,j,ii);
      }
      if (cs_ >= 0) {
        if (j > js) rr -= iw_(m,M1_IW_CJM,k,j,ii)*iw_(m,cs_,k,j-1,ii);
        if (j < je) rr -= iw_(m,M1_IW_CJP,k,j,ii)*iw_(m,cs_,k,j+1,ii);
        if (thrd) {
          if (k > ks) rr -= iw_(m,M1_IW_CKM,k,j,ii)*iw_(m,cs_,k-1,j,ii);
          if (k < ke) rr -= iw_(m,M1_IW_CKP,k,j,ii)*iw_(m,cs_,k+1,j,ii);
        }
      }
      sw(3*nx + i) = static_cast<T>(rr);
    });
    tm.team_barrier();
    int src = 0;
    for (int rd=0, s=1; rd<nround; ++rd, s*=2) {
      const int o = src*nv*nx, d = (1-src)*nv*nx;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        const int im = i - s, ip = i + s;
        T ai = sw(o + i), bi = sw(o + nx + i), ci = sw(o + 2*nx + i);
        T ri = sw(o + 3*nx + i);
        T ui = cyclic ? sw(o + 4*nx + i) : static_cast<T>(0.0);
        T an = 0.0, cn = 0.0;
        if (im >= 0) {
          T f = -ai/sw(o + nx + im);
          an = f*sw(o + im);
          bi += f*sw(o + 2*nx + im);
          ri += f*sw(o + 3*nx + im);
          if (cyclic) ui += f*sw(o + 4*nx + im);
        }
        if (ip < nx) {
          T g = -ci/sw(o + nx + ip);
          cn = g*sw(o + 2*nx + ip);
          bi += g*sw(o + ip);
          ri += g*sw(o + 3*nx + ip);
          if (cyclic) ui += g*sw(o + 4*nx + ip);
        }
        sw(d + i) = an;
        sw(d + nx + i) = bi;
        sw(d + 2*nx + i) = cn;
        sw(d + 3*nx + i) = ri;
        if (cyclic) sw(d + 4*nx + i) = ui;
      });
      tm.team_barrier();
      src = 1 - src;
    }
    const int o = src*nv*nx;
    if (!cyclic) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        iw_(m,cz,k,j,i+is) = static_cast<Real>(sw(o + 3*nx + i)/sw(o + nx + i));
      });
    } else {
      T y0 = sw(o + 3*nx)/sw(o + nx);
      T yn = sw(o + 4*nx - 1)/sw(o + 2*nx - 1);
      T z0 = sw(o + 4*nx)/sw(o + nx);
      T zn = sw(o + 5*nx - 1)/sw(o + 2*nx - 1);
      T bg = static_cast<T>(beta/gam);
      T fac = (y0 + bg*yn)/(static_cast<T>(1.0) + z0 + bg*zn);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        T bi = sw(o + nx + i);
        iw_(m,cz,k,j,i+is) = static_cast<Real>(sw(o + 3*nx + i)/bi
                                               - fac*(sw(o + 4*nx + i)/bi));
      });
    }
  });
}
} // namespace

void RadiationM1::ImplicitPCRSolveX(int rc, int zc, int upd, Real c1, Real c2, int col,
                                    int sub) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = ie - is + 1;
  const int nj = je - js + 1, nk = ke - ks + 1;
  const bool colr = (col >= 0);
  const int njl = colr ? (nj + 1)/2 : nj;   // league columns per (m,k)
  const int nkj = nk*njl;
  const int cl_ = colr ? col : 0;
  const int cs_ = (colr && sub >= 0) ? sub : -1;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  Kokkos::TeamPolicy<DevExeSpace> policy;
  if (!std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    int ts = impl_pcr_team;
    if (ts == 0) {
      ts = 1;
      while (ts < nx && ts < 256) ts *= 2;
    }
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, ts);
  } else {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, Kokkos::AUTO);
  }
  if (impl_prec_float) {
    M1PCRX<float>(iw, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                  trans_x3, cyclic, rc, zc, upd, c1, c2);
  } else {
    M1PCRX<Real>(iw, policy, is, ie, js, je, ks, ke, nkj, njl, colr, cl_, cs_,
                 trans_x3, cyclic, rc, zc, upd, c1, c2);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPrecondX
//! \brief the preconditioner of the implicit_krylov_fuse path, z = M^{-1} r with r from
//! `rc` or made by the p / s update `upd` (ImplicitPCRSolveX):
//!  implicit_precond = line: M = the x1 line part of the row (as ImplicitPrecond);
//!  implicit_precond = rbgs: ONE symmetric red-black transverse line Gauss-Seidel sweep
//!    from z = 0 (red, black, red; each colour one line solve of its columns against
//!    the just-updated other colour), block-local: the x2/x3 couplings across a
//!    MeshBlock face are left out of M (they stay in the operator), so it needs no
//!    communication.  M is symmetric positive definite whenever the 5-point part is, and
//!    stronger than line Jacobi on the transverse coupling (an approximate inverse of
//!    the whole 5-point-per-line system instead of its x1 part alone).
//!  implicit_precond = rbgs_fwd: the forward half only (red, black).

void RadiationM1::ImplicitPrecondX(int rc, int zc, int upd, Real c1, Real c2) {
  if (impl_prec == 0) {
    ImplicitPCRSolveX(rc, zc, upd, c1, c2, -1, -1);
    return;
  }
  // the right-hand side the third half-sweep re-reads: the update the first one wrote
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 0, -1);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 1, zc);
  if (impl_prec == 1) {
    ImplicitPCRSolveX(rr, zc, 0, 0.0, 0.0, 0, zc);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitODCache
//! \brief implicit_od_cache: odc(m,d,k,j,i) = M1OffDiv(x, d) at every cell a face of
//! ImplicitOffDiagOp reads it from -- the active box plus one layer in x1, x2 (and x3)
//! -- with exactly the index limits that routine passes, so every face value
//! 0.5*(odc_L + odc_R) is the very number ImplicitOffDiagOp forms (bitwise), from 3
//! evaluations per cell instead of 12.  Cells the faces never read are computed and
//! ignored.

void RadiationM1::ImplicitODCache(int xc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  const bool thrd = trans_x3;
  const int e3 = thrd ? 1 : 0;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto od_ = odc;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  const int cx = xc;
  par_for("m1_impl_odc", DevExeSpace(), 0, nmb1, ks-e3, ke+e3, js-1, je+1, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dx1 = mbsize(m).dx1;
    Real dx2 = mbsize(m).dx2;
    Real dx3 = mbsize(m).dx3;
    BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
    BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
    if (!p2lo) {jl = js-1;}
    if (!p2hi) {ju = je+1;}
    if (thrd && !p3lo) {kl = ks-1;}
    if (thrd && !p3hi) {ku = ke+1;}
    od_(m,0,k,j,i) = M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull);
    od_(m,1,k,j,i) = M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull);
    if (thrd) {
      od_(m,2,k,j,i) = M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirectInit
//! \brief implicit_halo_direct: tabulate, for every MeshBlock of the pack and each of the
//! 26 directions (ox1,ox2,ox3), the LOCAL index of the same-level neighbour that fills
//! that ghost region, or -1 where there is none (a physical boundary: its ghost zones
//! are not filled by the ordinary exchange either).  The direct copy is used only when
//! EVERY rank finds every neighbour on its own rank at the same level, outside the
//! cubed-sphere and polar transforms (so all ranks take the same branch and no MPI
//! exchange is left half-posted); otherwise the ordinary exchange runs, as before.

void RadiationM1::ImplicitHaloDirectInit() {
  auto *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack;
  hd_src = DualArray2D<int>("m1_hd_src", nmb, 27);
  int ok = (pm->multilevel || pm->use_cubed_sphere || pm->use_polar_boundary) ? 0 : 1;
  auto &nb = pmy_pack->pmb->nghbr;
  auto &lev = pmy_pack->pmb->mb_lev;
  const int e2 = pm->multi_d ? 1 : 0;
  const int e3 = pm->three_d ? 1 : 0;
  for (int m = 0; m < nmb; ++m) {
    for (int d = 0; d < 27; ++d) {hd_src.h_view(m,d) = -1;}
    for (int o3 = -e3; o3 <= e3; ++o3) {
      for (int o2 = -e2; o2 <= e2; ++o2) {
        for (int o1 = -1; o1 <= 1; ++o1) {
          if (o1 == 0 && o2 == 0 && o3 == 0) continue;
          int n = NeighborIndex(o1, o2, o3, 0, 0);
          if (n < 0 || n >= pmy_pack->pmb->nnghbr) {ok = 0; continue;}
          const NeighborBlock &q = nb.h_view(m,n);
          if (q.gid < 0) continue;
          if (q.rank != global_variable::my_rank || q.lev != lev.h_view(m)) {
            ok = 0;
            continue;
          }
          hd_src.h_view(m, (o1+1) + 3*(o2+1) + 9*(o3+1)) = q.gid - pmy_pack->gids;
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  {int g = ok;
  MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  ok = g;}
#endif
  halo_direct_on = (ok == 1);
  hd_src.modify_host();
  hd_src.sync_device();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitHaloDirect
//! \brief implicit_halo_direct: the ghost zones of `nq` components (c0 >= 0: that one
//! component; else the M1HaloCompT list) filled by ONE kernel that copies each ghost
//! cell from the active cell of the neighbour that owns it (all neighbours are on this
//! rank at the same level, ImplicitHaloDirectInit).  The same numbers the pack /
//! exchange / unpack chain delivers: 1 launch where there were 4, and no host work.

void RadiationM1::ImplicitHaloDirect(int nq, int c0) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int n1 = nx1 + 2*indcs.ng;
  const int n2 = (nx2 > 1) ? (nx2 + 2*indcs.ng) : 1;
  const int n3 = (nx3 > 1) ? (nx3 + 2*indcs.ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool md = (nx2 > 1), td = (nx3 > 1);
  auto iw_ = iw;
  auto tab = hd_src.d_view;
  const int nc0 = c0;
  par_for("m1_impl_hdir", DevExeSpace(), 0, nmb1, 0, nq-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    const int o1 = (i < is) ? -1 : ((i > ie) ? 1 : 0);
    const int o2 = md ? ((j < js) ? -1 : ((j > je) ? 1 : 0)) : 0;
    const int o3 = td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
    if (o1 == 0 && o2 == 0 && o3 == 0) return;
    const int src = tab(m, (o1+1) + 3*(o2+1) + 9*(o3+1));
    if (src < 0) return;
    const int nc = (nc0 >= 0) ? nc0 : M1HaloCompT(n);
    iw_(m,nc,k,j,i) = iw_(src,nc,k - o3*nx3,j - o2*nx2,i - o1*nx1);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitKrylovHalo
//! \brief milestone 3b phase C: put ONE component of the work array into the scratch
//! array `krw`, exchange it with all six neighbours through the module's ordinary
//! cell-centred boundary machinery (which is what supplies periodic wrap, edge/corner
//! neighbours and MPI), and copy the ghost zones back.  PHYSICAL boundaries are not
//! filled: the row's coefficient towards such a neighbour is identically zero, so the
//! ghost value is multiplied by zero and never read in anger.

void RadiationM1::ImplicitKrylovHalo(int comp) {
  ImplicitHaloExchange(1, comp);
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1OdFaces
//! \brief implicit_od_cache: the six face terms of ImplicitOffDiagOp at cell (m,k,j,i),
//! in the same order and with the same arithmetic, from the per-cell cache od(m,d,...)
//! = M1OffDiv(x, d) instead of re-evaluating M1OffDiv at both cells of every face.
//! The caller has already returned on an M1_IBC_EFIX row.

KOKKOS_INLINE_FUNCTION
Real M1OdFaces(const DvceArray5D<Real> &iw_, const DvceArray5D<Real> &od_,
               const DvceArray4D<Real> &th2_, const DvceArray4D<Real> &th3_,
               const bool lm, const int m, const int k, const int j, const int i,
               const int is, const int ie, const int js, const int je, const int ks,
               const int ke, const bool cyclic, const bool botb, const bool topb,
               const bool p2lo, const bool p2hi, const bool p3lo, const bool p3hi,
               const bool thrd, const Real dx1, const Real dx2, const Real dx3,
               const Real ch, const Real cl, const Real dt) {
  Real cr = ch/cl;
  Real kk = ch*cl*dt;
  Real y = 0.0;
  if (i < ie || cyclic || !topb) {
    int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,0,k,j,i) + od_(m,0,k,j,ip));
    y -= (dt/dx1)*cr*th*kk*od;
  }
  if (i > is || cyclic || !botb) {
    int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
    Real th = 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,0,k,j,im) + od_(m,0,k,j,i));
    y += (dt/dx1)*cr*th*kk*od;
  }
  Real nu2 = dt/dx2;
  if (!(j == je && p2hi)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
    Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,1,k,j,i) + od_(m,1,k,j+1,i));
    y -= nu2*cr*th*kk*od;
  }
  if (!(j == js && p2lo)) {
    Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
    Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
    Real od = 0.5*(od_(m,1,k,j-1,i) + od_(m,1,k,j,i));
    y += nu2*cr*th*kk*od;
  }
  if (thrd) {
    Real nu3 = dt/dx3;
    if (!(k == ke && p3hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
      Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(od_(m,2,k,j,i) + od_(m,2,k+1,j,i));
      y -= nu3*cr*th*kk*od;
    }
    if (!(k == ks && p3lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(od_(m,2,k-1,j,i) + od_(m,2,k,j,i));
      y += nu3*cr*th*kk*od;
    }
  }
  return y;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOffDiagOpC
//! \brief implicit_od_cache: ImplicitOffDiagOp (y += sgn L_off(x)) from the cache, and,
//! with `with7`, the whole operator y = A x in ONE kernel (the 7-point row of
//! ImplicitApplyOp, then + L_off(x) exactly as ImplicitApplyOp + ImplicitOffDiagOp
//! form it).  The reduction rides in the operator kernel (implicit_krylov_fuse >= 2):
//! red = 1: out[0] = (rhat,y);  red = 2: out[0] = (y,s), out[1] = (y,y);
//! red = 3: as 2 plus out[2] = (rhat,y);  red = 4: as 1 plus out[3] = max|r|.

void RadiationM1::ImplicitOffDiagOpC(int xc, int yc, Real sgn, bool with7, int red,
                                     Real *out) {
  ImplicitODCache(xc);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto od_ = odc;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  const bool lm = (impl_tlim != M1_TLIM_NONE);
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const int cx = xc, cy = yc;
  const Real sg = sgn;
  const bool w7 = with7;
  const int rm = red;
  // the standalone call (with7 = false) is made only under the operator form
  const bool odon = w7 ? (od_now == M1_OD_OPERATOR) : true;
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  // one cell's result; returns y (the value stored)
  auto row = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) -> Real {
    Real y7 = 0.0;
    if (w7) {
      int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
      int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
      y7 = iw_(m,M1_IW_TB,k,j,i)*iw_(m,cx,k,j,i)
           + iw_(m,M1_IW_TA,k,j,i)*iw_(m,cx,k,j,im)
           + iw_(m,M1_IW_TC,k,j,i)*iw_(m,cx,k,j,ip)
           + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,cx,k,j-1,i)
           + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,cx,k,j+1,i);
      if (thrd) {
        y7 += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,cx,k-1,j,i)
              + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,cx,k+1,j,i);
      }
    } else {
      y7 = iw_(m,cy,k,j,i);
    }
    if (!odon) {
      iw_(m,cy,k,j,i) = y7;
      return y7;
    }
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                    (i == ie && topb && bchi == M1_IBC_EFIX))) {
      iw_(m,cy,k,j,i) = y7;
      return y7;
    }
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    Real y = M1OdFaces(iw_, od_, th2_, th3_, lm, m, k, j, i, is, ie, js, je, ks, ke,
                       cyclic, botb, topb, p2lo, p2hi, p3lo, p3hi, thrd,
                       mbsize(m).dx1, mbsize(m).dx2, mbsize(m).dx3, ch, cl, dt);
    Real out = y7 + sg*y;
    iw_(m,cy,k,j,i) = out;
    return out;
  };
  if (rm == 0) {
    par_for("m1_impl_opc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      row(m, k, j, i);
    });
    return;
  }
  // 256-thread blocks: the default 1024-thread block of a reduction with a 4-Real
  // value takes 33 kB of LDS, i.e. ONE block per CU (measured 97 us vs 47 us for the
  // same stencil as a par_for)
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  Real a0 = 0.0, a1 = 0.0, a2 = 0.0, amx = 0.0;
  Kokkos::parallel_reduce("m1_impl_opcr", pol,
  KOKKOS_LAMBDA(const int idx, Real &l0, Real &l1, Real &l2, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real y = row(m, k, j, i);
    if (rm == 1 || rm == 4) {
      l0 += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        lmx = (a > lmx) ? a : lmx;
      }
    } else {
      l0 += y*iw_(m,M1_IW_KS,k,j,i);
      l1 += y*y;
      if (rm == 3) {l2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, a0, a1, a2, Kokkos::Max<Real>(amx));
  out[0] = a0;
  out[1] = a1;
  out[2] = a2;
  out[3] = amx;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1StIdx
//! \brief implicit_op_stencil: the slot of the 19-point stencil for the offset
//! (di,dj,dk) in {-1,0,1}^3 with at most two non-zero entries: 0 centre; 1..6 the faces
//! i-1,i+1,j-1,j+1,k-1,k+1; 7..10 the (i,j) edges, 11..14 the (i,k) edges, 15..18 the
//! (j,k) edges, each ordered (-,-),(+,-),(-,+),(+,+) in (first, second) axis.

KOKKOS_INLINE_FUNCTION
int M1StIdx(const int di, const int dj, const int dk) {
  if (dk == 0) {
    if (dj == 0) {return (di == 0) ? 0 : ((di < 0) ? 1 : 2);}
    if (di == 0) {return (dj < 0) ? 3 : 4;}
    return 7 + ((di > 0) ? 1 : 0) + ((dj > 0) ? 2 : 0);
  }
  if (dj == 0) {
    if (di == 0) {return (dk < 0) ? 5 : 6;}
    return 11 + ((di > 0) ? 1 : 0) + ((dk > 0) ? 2 : 0);
  }
  return 15 + ((dj > 0) ? 1 : 0) + ((dk > 0) ? 2 : 0);
}

//! D_de of the frozen closure at one cell (the coefficient M1POff multiplies x by)
KOKKOS_INLINE_FUNCTION
Real M1DOffC(const DvceArray5D<Real> &iw, const DvceArray5D<Real> &vd, const bool full,
             const int m, const int a, const int b, const int k, const int j,
             const int i) {
  if (full) {return vd(m,M1_VET_D11+2+a+b,k,j,i);}
  return M1EddOff(iw(m,M1_IW_WCHI,k,j,i), iw(m,M1_IW_N1+a,k,j,i),
                  iw(m,M1_IW_N1+b,k,j,i));
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitStencilBuild
//! \brief implicit_op_stencil: the frozen operator of this Picard pass -- the 7-point
//! row (TA,TB,TC,CJM..CKP) plus, under implicit_offdiag = operator, the off-diagonal
//! Eddington terms of ImplicitOffDiagOp -- written out ONCE as a 19-point stencil
//! st(m,0..18,k,j,i).  Every face term of ImplicitOffDiagOp,
//!   -/+ (dt/dx_d) (chat/c) theta_f chat c dt * 0.5 [OD_d(L) + OD_d(R)],
//!   OD_d(q) = sum_{e!=d} [D_de x](q_a) - [D_de x](q_b)) / ((a-b) dx_e)
//! with the same one-sided clamps at physical faces, is linear in x with coefficients
//! frozen over the pass, so it distributes onto the centre, the 6 face and the 12 edge
//! neighbours.  Applying the stencil (ImplicitStencilOp) is then one read of 19
//! coefficients per cell instead of re-deriving D_de and theta at every face in every
//! Krylov iteration.  Same operator; the sums are grouped differently (round-off).
//! Not for a periodic x1 wrap (cyclic), which keeps the od_cache path.

void RadiationM1::ImplicitStencilBuild() {
  if (ibc_x1min == M1_IBC_PERIODIC) {   // a problem generator may set it late
    ImplFatal("<rad_m1>/implicit_op_stencil does not take a periodic x1 wrap");
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto st_ = ost;
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool thrd = trans_x3;
  const bool lm = (impl_tlim != M1_TLIM_NONE);
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const bool odon = (od_now == M1_OD_OPERATOR);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Real emax = 0.0;
  Kokkos::parallel_reduce("m1_impl_stb",
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>(DevExeSpace(), 0,
                                                                (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real c[19];
    for (int o = 0; o < 19; ++o) {c[o] = 0.0;}
    c[0] = iw_(m,M1_IW_TB,k,j,i);
    c[1] = iw_(m,M1_IW_TA,k,j,i);
    c[2] = iw_(m,M1_IW_TC,k,j,i);
    c[3] = iw_(m,M1_IW_CJM,k,j,i);
    c[4] = iw_(m,M1_IW_CJP,k,j,i);
    if (thrd) {
      c[5] = iw_(m,M1_IW_CKM,k,j,i);
      c[6] = iw_(m,M1_IW_CKP,k,j,i);
    }
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    bool efix = (i == is && botb && bclo == M1_IBC_EFIX) ||
                (i == ie && topb && bchi == M1_IBC_EFIX);
    if (odon && !efix) {
      Real dxv[3] = {mbsize(m).dx1, mbsize(m).dx2, mbsize(m).dx3};
      BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
      BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
      BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
      BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
      BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
      BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
      bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
      bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
      bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
      bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
      // the index limits M1OffDiv is given, per axis
      int lo[3] = {is, js, ks}, hi[3] = {ie, je, ke};
      if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {lo[0] = is-1;}
      if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {hi[0] = ie+1;}
      if (!p2lo) {lo[1] = js-1;}
      if (!p2hi) {hi[1] = je+1;}
      if (thrd && !p3lo) {lo[2] = ks-1;}
      if (thrd && !p3hi) {hi[2] = ke+1;}
      const int cc[3] = {i, j, k};
      const Real cr = ch/cl, kk = ch*cl*dt;
      const int nd = thrd ? 3 : 2;
      for (int d = 0; d < nd; ++d) {
        for (int sd = -1; sd <= 1; sd += 2) {
          // does this face carry the term (the conditions of ImplicitOffDiagOp)?
          bool has;
          if (d == 0) {
            has = (sd > 0) ? (i < ie || !topb) : (i > is || !botb);
          } else if (d == 1) {
            has = (sd > 0) ? !(j == je && p2hi) : !(j == js && p2lo);
          } else {
            has = (sd > 0) ? !(k == ke && p3hi) : !(k == ks && p3lo);
          }
          if (!has) continue;
          // the face theta, from the two cells' transport opacities (or the limiter)
          int nb[3] = {i, j, k};
          nb[d] += sd;
          Real th;
          if (d > 0 && lm) {
            th = (d == 1) ? th2_(m,k,(sd > 0) ? j+1 : j,i)
                          : th3_(m,(sd > 0) ? k+1 : k,j,i);
          } else {
            Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,nb[2],nb[1],nb[0]));
            th = 1.0/(1.0 + ch*dt*ktf);
          }
          // upper face: y -= w od; lower face: y += w od; od = 0.5 (OD(c) + OD(nb))
          const Real w = -static_cast<Real>(sd)*(dt/dxv[d])*cr*th*kk*0.5;
          for (int qs = 0; qs <= 1; ++qs) {
            int q[3] = {i, j, k};
            if (qs == 1) {q[d] += sd;}
            for (int e = 0; e < nd; ++e) {
              if (e == d) continue;
              int qa[3] = {q[0], q[1], q[2]}, qb[3] = {q[0], q[1], q[2]};
              if (q[e] + 1 <= hi[e]) {qa[e] = q[e] + 1;}
              if (q[e] - 1 >= lo[e]) {qb[e] = q[e] - 1;}
              if (qa[e] == qb[e]) continue;
              const Real f = w/((qa[e] - qb[e])*dxv[e]);
              const Real da = M1DOffC(iw_, vd_, dfull, m, d, e, qa[2], qa[1], qa[0]);
              const Real db = M1DOffC(iw_, vd_, dfull, m, d, e, qb[2], qb[1], qb[0]);
              c[M1StIdx(qa[0]-cc[0], qa[1]-cc[1], qa[2]-cc[2])] += f*da;
              c[M1StIdx(qb[0]-cc[0], qb[1]-cc[1], qb[2]-cc[2])] -= f*db;
            }
          }
        }
      }
    }
    for (int o = 0; o < 19; ++o) {st_(m,o,k,j,i) = c[o];}
    for (int o = 7; o < 19; ++o) {lmx = fmax(lmx, fabs(c[o]));}
  }, Kokkos::Max<Real>(emax));
  // no edge coefficient anywhere (the Eddington closure: D_ab = 0 off the diagonal, or
  // implicit_offdiag not operator): ImplicitStencilOp reads the 7 face/centre slots only,
  // which adds the same non-zero terms in the same order
#if MPI_PARALLEL_ENABLED
  {Real g = emax;
  MPI_Allreduce(&emax, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  emax = g;}
#endif
  st_edges = (emax > 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitStencilOp
//! \brief implicit_op_stencil: y = A x from the 19-point stencil of ImplicitStencilBuild
//! (the caller has filled the ghost zones of x).  `red` and `out` as ImplicitOffDiagOpC.

void RadiationM1::ImplicitStencilOp(int xc, int yc, int red, Real *out) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto st_ = ost;
  const bool thrd = trans_x3;
  const int cx = xc, cy = yc;
  const int rm = red;
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  const bool edg = st_edges;
  auto row = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) -> Real {
    Real y = st_(m,0,k,j,i)*iw_(m,cx,k,j,i)
             + st_(m,1,k,j,i)*iw_(m,cx,k,j,i-1) + st_(m,2,k,j,i)*iw_(m,cx,k,j,i+1)
             + st_(m,3,k,j,i)*iw_(m,cx,k,j-1,i) + st_(m,4,k,j,i)*iw_(m,cx,k,j+1,i);
    if (edg) {
      y += st_(m,7,k,j,i)*iw_(m,cx,k,j-1,i-1) + st_(m,8,k,j,i)*iw_(m,cx,k,j-1,i+1)
           + st_(m,9,k,j,i)*iw_(m,cx,k,j+1,i-1) + st_(m,10,k,j,i)*iw_(m,cx,k,j+1,i+1);
    }
    if (thrd) {
      y += st_(m,5,k,j,i)*iw_(m,cx,k-1,j,i) + st_(m,6,k,j,i)*iw_(m,cx,k+1,j,i);
      if (edg) {
        y += st_(m,11,k,j,i)*iw_(m,cx,k-1,j,i-1) + st_(m,12,k,j,i)*iw_(m,cx,k-1,j,i+1)
             + st_(m,13,k,j,i)*iw_(m,cx,k+1,j,i-1) + st_(m,14,k,j,i)*iw_(m,cx,k+1,j,i+1)
             + st_(m,15,k,j,i)*iw_(m,cx,k-1,j-1,i) + st_(m,16,k,j,i)*iw_(m,cx,k-1,j+1,i)
             + st_(m,17,k,j,i)*iw_(m,cx,k+1,j-1,i) + st_(m,18,k,j,i)*iw_(m,cx,k+1,j+1,i);
      }
    }
    iw_(m,cy,k,j,i) = y;
    return y;
  };
  if (rm == 0) {
    par_for("m1_impl_sto", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      row(m, k, j, i);
    });
    return;
  }
  // 256-thread blocks: the default 1024-thread block of a reduction with a 4-Real
  // value takes 33 kB of LDS, i.e. ONE block per CU (measured 97 us vs 47 us for the
  // same stencil as a par_for)
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  Real a0 = 0.0, a1 = 0.0, a2 = 0.0, amx = 0.0;
  Kokkos::parallel_reduce("m1_impl_stor", pol,
  KOKKOS_LAMBDA(const int idx, Real &l0, Real &l1, Real &l2, Real &lmx) {
    int m = idx/nkji;
    int r = idx - m*nkji;
    int k = r/nji;
    r -= k*nji;
    int j = r/ni;
    int i = r - j*ni;
    k += ks; j += js; i += is;
    Real y = row(m, k, j, i);
    if (rm == 1 || rm == 4) {
      l0 += iw_(m,M1_IW_KRH,k,j,i)*y;
      if (rm == 4) {
        Real a = fabs(iw_(m,M1_IW_KR,k,j,i));
        lmx = (a > lmx) ? a : lmx;
      }
    } else {
      l0 += y*iw_(m,M1_IW_KS,k,j,i);
      l1 += y*y;
      if (rm == 3) {l2 += iw_(m,M1_IW_KRH,k,j,i)*y;}
    }
  }, a0, a1, a2, Kokkos::Max<Real>(amx));
  out[0] = a0;
  out[1] = a1;
  out[2] = a2;
  out[3] = amx;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOpX
//! \brief y = A x (+ the reductions `red` of ImplicitOffDiagOpC) on the fast path: from
//! the 19-point stencil under implicit_op_stencil, else from the od cache.  The caller
//! has filled the ghost zones of x.

void RadiationM1::ImplicitOpX(int xc, int yc, int red, Real *out) {
  if (impl_stencil) {
    ImplicitStencilOp(xc, yc, red, out);
  } else {
    ImplicitOffDiagOpC(xc, yc, 1.0, true, red, out);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOffDiagOp
//! \brief MILESTONE 3b phase D: accumulate sgn * L_off(x) into the component yc, where
//! L_off is the contribution the OFF-DIAGONAL Eddington terms of the face equations make
//! to the cell row:
//!
//!   L_off(x)_c = - sum_d (dt/dx_d) (chat/c) chat c dt
//!                  [ theta_{d,+} od_{d,+}(x) - theta_{d,-} od_{d,-}(x) ],
//!   od_{d,f}(x) = (1/2) [ (sum_{e!=d} d_e D_de x)_L + (sum_{e!=d} d_e D_de x)_R ],
//!
//! i.e. exactly the term the face fluxes of ImplicitTransverseTerms and of the x1
//! assembly carry, with the LAGGED energy replaced by the argument x.  The closure
//! (chi, n, hence D_ab) is frozen inside a Picard pass, so this IS a linear operator: a
//! 9-point stencil in 2-D and a 19-point one in 3-D, built from the same one-layer halo
//! the 7-point operator already exchanges (a cell-centred exchange fills the edge and
//! corner neighbours, which is what the cross derivatives read).
//!
//! It is used twice per Picard pass under implicit_offdiag = operator: once with x = E^k
//! to REMOVE the lagged term from the right-hand side the assembly produced, and inside
//! every operator application to put it back on the left.  Both use the same routine, so
//! the two cannot drift apart.
//!
//! The boundary logic mirrors the assembly face by face: a physical x1 face carries an
//! imposed flux and no off-diagonal term, a physical x2/x3 face is reflecting (F = 0),
//! and an M1_IBC_EFIX Dirichlet row is replaced whole and gets nothing at all.

void RadiationM1::ImplicitOffDiagOp(int xc, int yc, Real sgn) {
  if (impl_odc) {   // the same numbers from the per-cell cache
    ImplicitOffDiagOpC(xc, yc, sgn, false, 0, nullptr);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  // capture the DEVICE Views only: whole DualViews pushed this functor past 512 bytes,
  // i.e. onto Kokkos HIP's constant-memory launch, which waits on the previous such
  // launch (hip_event_synchronize) -- one hidden host stall per call, ~290 per cycle.
  auto mbsize = pmy_pack->pmb->mb_size.d_view;
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
  auto pos_ = part_pos.d_view;
  auto vd_ = vet_cell;   // vet_tensor = full (M1DDiag, M1OffDiv)
  const bool dfull = vet_full;
  const int nblkx1 = part_nblk;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  // the TRANSVERSE face theta of the operator must be the very number the face fluxes
  // and the cell terms used, limiter or not (ImplicitTransTheta).  The x1 faces are not
  // limited.
  const bool lm = (impl_tlim != M1_TLIM_NONE);
  auto th2_ = thx2;
  auto th3_ = thx3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  Real cl = c_light, ch = chat, dt = dt_sub;
  const int cx = xc, cy = yc;
  const Real sg = sgn;
  par_for("m1_impl_odop", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int ipos = pos_(m);
    bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                    (i == ie && topb && bchi == M1_IBC_EFIX))) {
      return;
    }
    Real dx1 = mbsize(m).dx1;
    Real dx2 = mbsize(m).dx2;
    Real dx3 = mbsize(m).dx3;
    BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
    BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
    int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
    if (!p2lo) {jl = js-1;}
    if (!p2hi) {ju = je+1;}
    if (thrd && !p3lo) {kl = ks-1;}
    if (thrd && !p3hi) {ku = ke+1;}
    Real cr = ch/cl;
    Real kk = ch*cl*dt;
    Real y = 0.0;
    // ---- the two x1 faces
    if (i < ie || cyclic || !topb) {
      int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,0,k,j,ip,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y -= (dt/dx1)*cr*th*kk*od;
    }
    if (i > is || cyclic || !botb) {
      int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
      Real th = 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,0,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y += (dt/dx1)*cr*th*kk*od;
    }
    // ---- the two x2 faces
    Real nu2 = dt/dx2;
    if (!(j == je && p2hi)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
      Real th = lm ? th2_(m,k,j+1,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,1,k,j+1,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y -= nu2*cr*th*kk*od;
    }
    if (!(j == js && p2lo)) {
      Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
      Real th = lm ? th2_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
      Real od = 0.5*(M1OffDiv(iw_,m,1,k,j-1,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                              dfull)
                     + M1OffDiv(iw_,m,1,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull));
      y += nu2*cr*th*kk*od;
    }
    // ---- the two x3 faces
    if (thrd) {
      Real nu3 = dt/dx3;
      if (!(k == ke && p3hi)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
        Real th = lm ? th3_(m,k+1,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real od = 0.5*(M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull)
                       + M1OffDiv(iw_,m,2,k+1,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,
                                  vd_,dfull));
        y -= nu3*cr*th*kk*od;
      }
      if (!(k == ks && p3lo)) {
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
        Real th = lm ? th3_(m,k,j,i) : 1.0/(1.0 + ch*dt*ktf);
        Real od = 0.5*(M1OffDiv(iw_,m,2,k-1,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                dfull)
                       + M1OffDiv(iw_,m,2,k,j,i,dx1,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,cx,vd_,
                                  dfull));
        y += nu3*cr*th*kk*od;
      }
    }
    iw_(m,cy,k,j,i) += sg*y;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitApplyOp
//! \brief milestone 3b phase C: y = A x for the FROZEN 7-point operator of the current
//! Picard pass -- the assembled tridiagonal row (M1_IW_TA, TB, TC, with TB already
//! carrying the transverse diagonal M1_IW_TDIA) plus the four/six transverse
//! off-diagonals (M1_IW_CJM..CKP).  One halo exchange of x, then one stencil kernel.
//!
//! The x1 wrap of a periodic single-block column is handled exactly as the assembly
//! handles it; with several blocks stacked along x1 the coupling to the neighbouring
//! block is through the ghost cell the halo just filled, which is the same cell the
//! gathered line solve treats as an interior row.

void RadiationM1::ImplicitApplyOp(int xc, int yc) {
  ImplicitKrylovHalo(xc);
  if (impl_stencil) {   // the 19-point stencil of this pass
    ImplicitStencilOp(xc, yc, 0, nullptr);
    return;
  }
  if (impl_odc) {   // 7-point row + off-diagonal Eddington terms in one kernel
    ImplicitOffDiagOpC(xc, yc, 1.0, true, 0, nullptr);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const bool cyclic = (ibc_x1min == M1_IBC_PERIODIC);
  const bool thrd = trans_x3;
  const int cx = xc, cy = yc;
  par_for("m1_impl_op", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
    int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
    Real y = iw_(m,M1_IW_TB,k,j,i)*iw_(m,cx,k,j,i)
             + iw_(m,M1_IW_TA,k,j,i)*iw_(m,cx,k,j,im)
             + iw_(m,M1_IW_TC,k,j,i)*iw_(m,cx,k,j,ip)
             + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,cx,k,j-1,i)
             + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,cx,k,j+1,i);
    if (thrd) {
      y += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,cx,k-1,j,i)
           + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,cx,k+1,j,i);
    }
    iw_(m,cy,k,j,i) = y;
  });
  // MILESTONE 3b phase D: the off-diagonal Eddington coupling, when it is part of the
  // operator.  The halo of x above is a full cell-centred exchange, so the edge and
  // corner neighbours the cross derivatives read are already filled: the 9-/19-point
  // stencil costs one extra kernel and no communication at all.
  if (od_now == M1_OD_OPERATOR) {
    ImplicitOffDiagOp(xc, yc, 1.0);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPrecond
//! \brief milestone 3b phase C: z = M^{-1} r, with M the x1 line part of the row (the
//! exact tridiagonal solve including the full diagonal).  It stages r through M1_IW_TR,
//! which the assembly has already been read out of by the time the Krylov loop runs.

void RadiationM1::ImplicitPrecond(int rc, int zc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const int cr = rc, cz = zc;
  if (cr >= 0) {
    par_for("m1_impl_prein", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_TR,k,j,i) = iw_(m,cr,k,j,i);
    });
  }
  ImplicitTridiagSolve();
  par_for("m1_impl_preout", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,cz,k,j,i) = iw_(m,M1_IW_S2,k,j,i);
  });
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn M1GlobalSum2
//! \brief two global SUMs in one MPI_Allreduce.  The per-rank part is a Kokkos reduction
//! over the same fixed index range every time, so its summation order is reproducible on
//! a given rank count, and the cross-rank part is one collective.

void M1GlobalSum2(Real &a, Real &b) {
#if MPI_PARALLEL_ENABLED
  Real loc[2] = {a, b}, glb[2];
  MPI_Allreduce(loc, glb, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  a = glb[0];
  b = glb[1];
#else
  (void)a;
  (void)b;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1GlobalSumArr
//! \brief n global SUMs in one MPI_Allreduce (milestone 3e: the row of the Anderson
//! normal equations).  Same reproducibility argument as M1GlobalSum2.

void M1GlobalSumArr(Real *a, const int n) {
#if MPI_PARALLEL_ENABLED
  Real glb[NREDUCTION_VARIABLES];
  MPI_Allreduce(a, glb, n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int q = 0; q < n; ++q) {a[q] = glb[q];}
#else
  (void)a;
  (void)n;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1SolveSPD
//! \brief solve the small SYMMETRIC system A g = b in place by Gaussian elimination with
//! partial pivoting (n <= M1_AND_MMAX).  Returns false if the matrix is numerically
//! singular, in which case the caller falls back to a plain Picard pass.

bool M1SolveSPD(Real *a, Real *b, const int n) {
  for (int p = 0; p < n; ++p) {
    int piv = p;
    Real amx = fabs(a[p*n + p]);
    for (int r = p+1; r < n; ++r) {
      Real v = fabs(a[r*n + p]);
      if (v > amx) {amx = v; piv = r;}
    }
    if (!(amx > 0.0)) {return false;}
    if (piv != p) {
      for (int q = 0; q < n; ++q) {
        Real t = a[p*n + q];
        a[p*n + q] = a[piv*n + q];
        a[piv*n + q] = t;
      }
      Real t = b[p];
      b[p] = b[piv];
      b[piv] = t;
    }
    Real ip = 1.0/a[p*n + p];
    for (int r = p+1; r < n; ++r) {
      Real fct = a[r*n + p]*ip;
      if (fct == 0.0) {continue;}
      for (int q = p; q < n; ++q) {a[r*n + q] -= fct*a[p*n + q];}
      b[r] -= fct*b[p];
    }
  }
  for (int r = n-1; r >= 0; --r) {
    Real s = b[r];
    for (int q = r+1; q < n; ++q) {s -= a[r*n + q]*b[q];}
    b[r] = s/a[r*n + r];
  }
  for (int r = 0; r < n; ++r) {
    if (!std::isfinite(b[r])) {return false;}
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn M1GlobalMax
void M1GlobalMax(Real &a) {
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&a, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  a = g;
#else
  (void)a;
#endif
}

//----------------------------------------------------------------------------------------
//! \struct M1BcgVal, M1BcgRed
//! \brief implicit_bcg_sync > 0: ONE reduction kernel returns up to three sums and one
//! max (e.g. rhat.r and max|r| straight from the kernel that updates r).  M1BcgRed is a
//! Kokkos reducer of the same shape as Kokkos::Sum; |r| >= 0, so 0 is the identity of
//! the max slot.

struct M1BcgVal {
  Real s0, s1, s2, mx;
};

template <class Space>
struct M1BcgRed {
 public:
  using reducer = M1BcgRed<Space>;
  using value_type = M1BcgVal;
  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool refs;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit M1BcgRed(value_type &v) : value(&v), refs(true) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    d.s0 += s.s0;
    d.s1 += s.s1;
    d.s2 += s.s2;
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    v.s0 = 0.0;
    v.s1 = 0.0;
    v.s2 = 0.0;
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return refs;}
};

//! the flattened (m,k,j,i) index of the fused reductions, as par_for flattens it
KOKKOS_INLINE_FUNCTION
void M1BcgIdx(const int idx, const int nkji, const int nji, const int ni,
              int &m, int &k, int &j, int &i) {
  m = idx/nkji;
  int r = idx - m*nkji;
  k = r/nji;
  r -= k*nji;
  j = r/ni;
  i = r - j*ni;
}

#if MPI_PARALLEL_ENABLED
//! the MPI operation of M1GlobalBcg: sum the first three slots, max the fourth
void M1BcgOpFn(void *in, void *inout, int *len, MPI_Datatype *) {
  Real *a = static_cast<Real *>(in);
  Real *b = static_cast<Real *>(inout);
  for (int q = 0; q < *len; ++q) {
    b[4*q] += a[4*q];
    b[4*q + 1] += a[4*q + 1];
    b[4*q + 2] += a[4*q + 2];
    b[4*q + 3] = (a[4*q + 3] > b[4*q + 3]) ? a[4*q + 3] : b[4*q + 3];
  }
}
#endif

//----------------------------------------------------------------------------------------
//! \fn M1GlobalBcg
//! \brief the three sums and the max of an M1BcgVal over all ranks in ONE MPI_Allreduce
//! (a 4-Real contiguous type with a user operation).  Nothing to do on one rank.

void M1GlobalBcg(M1BcgVal &v) {
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks == 1) {return;}
  static MPI_Datatype typ = MPI_DATATYPE_NULL;
  static MPI_Op op = MPI_OP_NULL;
  if (op == MPI_OP_NULL) {
    MPI_Type_contiguous(4, MPI_ATHENA_REAL, &typ);
    MPI_Type_commit(&typ);
    MPI_Op_create(&M1BcgOpFn, 1, &op);
  }
  Real loc[4] = {v.s0, v.s1, v.s2, v.mx}, glb[4];
  MPI_Allreduce(loc, glb, 1, typ, op, MPI_COMM_WORLD);
  v.s0 = glb[0];
  v.s1 = glb[1];
  v.s2 = glb[2];
  v.mx = glb[3];
#else
  (void)v;
#endif
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitAccelSave
//! \brief MILESTONE 3e: snapshot the SCALED state x_k entering a Picard pass.
//!
//! The fixed-point vector is the COMPLETE lagged state one pass maps to the next: the
//! solve unknown E and the three cell-centred fluxes F_1, F_2, F_3, which are what the
//! top of the next pass rebuilds the closure (chi, n, hence the whole Eddington tensor)
//! from.  Nothing else a pass reads is independent state: the face fluxes f0x1/f0x2/f0x3
//! are recomputed from E at the end of every pass, the velocities and opacities are
//! frozen over the step, and the temperature is slaved to E by the scalar root find.
//!
//! SCALING.  E and F are not commensurate (F ~ c E), and an unscaled 2-norm would be
//! dominated by the optically thick base, where E is ~ 10 orders above the thin top --
//! which is exactly where the iteration does NOT need help.  Each cell is therefore
//! divided by its OWN energy scale S = max(E^n, e_floor), fixed over the whole step (so
//! the least-squares problem stays linear across passes), and the flux components by c S:
//!   x = ( E/S , F_1/(c S) , F_2/(c S) , F_3/(c S) ) ,
//! i.e. the Anderson residual is the RELATIVE fixed-point residual, cell by cell, and the
//! thin top carries the same weight as the base.

void RadiationM1::ImplicitAccelSave() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto xc_ = aa_xc;
  auto sc_ = aa_sc;
  const Real cl = c_light;
  const int nc = aa_nc;
  par_for("m1_acc_save", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
    Real w = 1.0/(sc_(m,k,j,i)*((c == 0) ? 1.0 : cl));
    xc_(m,c,k,j,i) = iw_(m,M1AccComp(c),k,j,i)*w;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitAccelApply
//! \brief MILESTONE 3e: ANDERSON ACCELERATION (Walker & Ni 2011, SIAM J. Numer.
//! Anal. 49, 1715) of the Picard map, applied at the END of pass `it`, where the iw
//! state is G(x_k).
//!
//! With g_k = G(x_k) - x_k and the histories dX_t = x_{t+1} - x_t and dG_t = g_{t+1} -
//! g_t of the last m passes, the accelerated iterate is
//!
//!   gamma = argmin_gamma || g_k - dG gamma ||_2 ,
//!   x_{k+1} = x_k + beta g_k - (dX + beta dG) gamma ,
//!
//! i.e. the beta-mixed map evaluated at the point of the affine span of the last m
//! iterates whose linear residual model is smallest.  beta = 1 is the plain (undamped)
//! form.  The least-squares problem is solved on the HOST from the NORMAL equations
//! (dG^T dG + lambda I) gamma = dG^T g_k, with lambda = M1_AND_REG tr(dG^T dG)/m a
//! Tikhonov term that makes a nearly dependent history harmless; every inner product is a
//! GLOBAL sum (one MPI_Allreduce of m+1 doubles per history row), so the answer does not
//! depend on the decomposition.
//!
//! REALIZABILITY is re-applied to the accelerated iterate (E >= e_floor, |F| <= c E)
//! before it is written back: the extrapolation is a linear combination of realizable
//! states and the M1 admissible set is convex in (E,F), but the E floor and the per-cell
//! scaling break exact convexity, and an inadmissible lagged state would give the next
//! pass a closure with chi outside [1/3,1].
//!
//! SAFEGUARD: if ||g_k|| exceeds 10x the previous pass', the history is dropped and the
//! pass falls back to plain Picard (x_{k+1} = G(x_k)); the event is counted.  The
//! convergence test of the outer loop is untouched -- it still measures |dE|/E and the
//! true linear residual of the UNACCELERATED map, which is the fixed-point residual.

void RadiationM1::ImplicitAccelApply(int it) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto xc_ = aa_xc;
  auto fc_ = aa_fc;
  auto xp_ = aa_xp;
  auto fp_ = aa_fp;
  auto dx_ = aa_dx;
  auto df_ = aa_df;
  auto sc_ = aa_sc;
  const Real cl = c_light;
  const Real efl = e_floor;
  const int nc = aa_nc;
  const int mm = impl_and_m;

  // (1) the fixed-point residual of this pass, in the scaled variables
  par_for("m1_acc_res", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
    Real w = 1.0/(sc_(m,k,j,i)*((c == 0) ? 1.0 : cl));
    fc_(m,c,k,j,i) = iw_(m,M1AccComp(c),k,j,i)*w - xc_(m,c,k,j,i);
  });

  // (2) its GLOBAL 2-norm
  Real fn2 = 0.0;
  Kokkos::parallel_reduce("m1_acc_nrm",
  Kokkos::MDRangePolicy<Kokkos::Rank<5>>(DevExeSpace(), {0,0,ks,js,is},
                                         {nmb1+1,nc,ke+1,je+1,ie+1}),
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i,
                Real &lsum) {
    lsum += SQR(fc_(m,c,k,j,i));
  }, Kokkos::Sum<Real>(fn2));
  {Real dum = 0.0;
  M1GlobalSum2(fn2, dum);}
  Real fnorm = sqrt(fmax(fn2, 0.0));

  // (3) the divergence safeguard
  bool restart = false;
  if (aa_fnp > 0.0 && fnorm > 10.0*aa_fnp) {
    aa_nh = 0;
    aa_head = 0;
    aa_hasp = false;
    aa_nrst += 1.0;
    restart = true;
  }

  // (4) push (dX, dG) onto the ring, then remember this pass
  if (aa_hasp) {
    const int q = aa_head;
    par_for("m1_acc_push", DevExeSpace(), 0, nmb1, 0, nc-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i) {
      dx_(m,q,c,k,j,i) = xc_(m,c,k,j,i) - xp_(m,c,k,j,i);
      df_(m,q,c,k,j,i) = fc_(m,c,k,j,i) - fp_(m,c,k,j,i);
    });
    aa_head = (aa_head + 1) % mm;
    aa_nh = std::min(aa_nh + 1, mm);
  }
  Kokkos::deep_copy(DevExeSpace(), aa_xp, aa_xc);
  Kokkos::deep_copy(DevExeSpace(), aa_fp, aa_fc);
  aa_hasp = true;
  aa_fnp = fnorm;

  // nothing to extrapolate from (or the safeguard fired, or the pass is before
  // implicit_anderson_start): leave the iw state as G(x_k), which IS the Picard iterate
  if (restart || it < impl_and_start || aa_nh < 1) {return;}

  // (5) the m x m normal equations, one GLOBAL reduction of m+1 numbers per row.  The
  // history columns in use are the last aa_nh pushed onto the ring.
  const int nh = aa_nh;
  int cid[M1_AND_MMAX];
  for (int t = 0; t < nh; ++t) {cid[t] = (aa_head - nh + t + mm) % mm;}
  Real amat[M1_AND_MMAX*M1_AND_MMAX], bvec[M1_AND_MMAX];
  Real trc = 0.0;
  for (int t = 0; t < nh; ++t) {
    const int qr = cid[t];
    array_sum::GlobalSum row;
    Kokkos::parallel_reduce("m1_acc_dot",
    Kokkos::MDRangePolicy<Kokkos::Rank<5>>(DevExeSpace(), {0,0,ks,js,is},
                                           {nmb1+1,nc,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int c, const int k, const int j, const int i,
                  array_sum::GlobalSum &lsum) {
      Real dq = df_(m,qr,c,k,j,i);
      for (int r = 0; r < nh; ++r) {
        lsum.the_array[r] += dq*df_(m,cid[r],c,k,j,i);
      }
      lsum.the_array[nh] += dq*fc_(m,c,k,j,i);
    }, Kokkos::Sum<array_sum::GlobalSum>(row));
    M1GlobalSumArr(row.the_array, nh+1);
    for (int r = 0; r < nh; ++r) {amat[t*nh + r] = row.the_array[r];}
    bvec[t] = row.the_array[nh];
    trc += row.the_array[t];
  }
  const Real lam = M1_AND_REG*fmax(trc, 1.0e-300)/static_cast<Real>(nh);
  for (int t = 0; t < nh; ++t) {amat[t*nh + t] += lam;}
  if (!M1SolveSPD(amat, bvec, nh)) {return;}

  // (6) the accelerated iterate, with the realizability limits re-applied
  Real gam[M1_AND_MMAX];
  for (int t = 0; t < nh; ++t) {gam[t] = bvec[t];}
  const Real bta = impl_and_beta;
  par_for("m1_acc_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real xn[4] = {0.0, 0.0, 0.0, 0.0};
    for (int c = 0; c < nc; ++c) {
      Real v = xc_(m,c,k,j,i) + bta*fc_(m,c,k,j,i);
      for (int t = 0; t < nh; ++t) {
        const int q = cid[t];
        v -= gam[t]*(dx_(m,q,c,k,j,i) + bta*df_(m,q,c,k,j,i));
      }
      xn[c] = v;
    }
    Real s = sc_(m,k,j,i);
    Real e = xn[0]*s;
    if (!(e > efl)) {e = efl;}
    Real f1 = xn[1]*cl*s;
    Real f2 = (nc > 2) ? (xn[2]*cl*s) : 0.0;
    Real f3 = (nc > 3) ? (xn[3]*cl*s) : 0.0;
    Real fm = sqrt(f1*f1 + f2*f2 + f3*f3);
    Real fmx = cl*e;
    if (fm > fmx) {
      Real sf = fmx/fm;
      f1 *= sf;
      f2 *= sf;
      f3 *= sf;
    }
    iw_(m,M1_IW_EP,k,j,i) = e;
    iw_(m,M1_IW_F1,k,j,i) = f1;
    if (nc > 2) {
      iw_(m,M1_IW_F2,k,j,i) = f2;
      iw_(m,M1_IW_F3,k,j,i) = f3;
    }
  });
  aa_nacc += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStab
//! \brief milestone 3b phase C: solve the FROZEN 7-point linear system of one Picard pass
//! by matrix-free BiCGStab, RIGHT-preconditioned by the exact x1 line solve, and leave
//! the answer in M1_IW_S2 -- exactly where the line-Jacobi pass leaves it, so nothing
//! downstream of step (e) of ImplicitSolve knows which solver ran.
//!
//! The system is A x = b with
//!   A x = tridiag(TA,TB,TC) x + CJM x_{j-1} + CJP x_{j+1} + CKM x_{k-1} + CKP x_{k+1},
//!   b   = TR + sum_nb C_nb E^k_nb,
//! i.e. the right-hand side line Jacobi uses PLUS the lagged off-diagonal term it had
//! moved there.  One line-Jacobi pass is x <- M^{-1}(b - sum C x) with the same M and the
//! same A, so the two solvers have the SAME fixed point and converge to the same answer;
//! only the number of passes differs.  The initial guess is the Picard iterate itself.
//!
//! Stopping is on the TRUE residual max|b - A x| relative to max|b|, below
//! implicit_lin_tol: the recursive residual is monitored every iteration (it is free)
//! and, when it passes, ONE extra operator application checks the true one.  If the true
//! residual disagrees the recurrence is RESTARTED from it, which is the standard cure for
//! the drift between the two.
//!
//! BREAKDOWN (rho or rhat.v underflowing, omega vanishing) restarts the recurrence once
//! from the current iterate with a fresh shadow residual; a second breakdown in the same
//! pass falls back to ONE line-Jacobi update, which always exists, and is counted.
//!
//! Three global reductions per iteration (rhat.r with r.r, rhat.v, and t.s with t.t),
//! five scalars in all.

int RadiationM1::ImplicitBiCGStab(Real rhsmax) {
  if (impl_bcg_sync > 0) {return ImplicitBiCGStabFused(rhsmax);}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  Kokkos::MDRangePolicy<Kokkos::Rank<4>> rng(DevExeSpace(), {0,ks,js,is},
                                             {nmb1+1,ke+1,je+1,ie+1});

  // x0 = the Picard iterate; r0 = b - A x0
  par_for("m1_impl_bcg_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  par_for("m1_impl_bcg_r0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
  });
  Real rnorm = 0.0;
  Kokkos::parallel_reduce("m1_impl_bcg_rn", rng,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
    Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
    lmax = (r > lmax) ? r : lmax;
  }, Kokkos::Max<Real>(rnorm));
  M1GlobalMax(rnorm);
  bcg_nred += 1.0;
  bcg_r0rel = rnorm/bscale;

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = (rnorm/bscale < tol);
  bool fell_back = false;
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    // rho_new = (rhat, r)
    Real rhon = 0.0;
    Kokkos::parallel_reduce("m1_impl_bcg_rho", rng,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
      ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KR,k,j,i);
    }, rhon);
    Real dummy = 0.0;
    M1GlobalSum2(rhon, dummy);
    bcg_nred += 1.0;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      const Real bt = beta, om = omega;
      par_for("m1_impl_bcg_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        iw_(m,M1_IW_KP,k,j,i) = iw_(m,M1_IW_KR,k,j,i)
            + bt*(iw_(m,M1_IW_KP,k,j,i) - om*iw_(m,M1_IW_KV,k,j,i));
      });
      ImplicitPrecond(M1_IW_KP, M1_IW_KY);
      ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
      Real rv = 0.0;
      Kokkos::parallel_reduce("m1_impl_bcg_rv", rng,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
        ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KV,k,j,i);
      }, rv);
      Real d2 = 0.0;
      M1GlobalSum2(rv, d2);
      bcg_nred += 1.0;
      if (!(fabs(rv) > M1_BCG_EPS)) {
        breakdown = true;
      } else {
        alpha = rhon/rv;
        const Real al = alpha;
        par_for("m1_impl_bcg_s", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          iw_(m,M1_IW_KS,k,j,i) = iw_(m,M1_IW_KR,k,j,i) - al*iw_(m,M1_IW_KV,k,j,i);
        });
        ImplicitPrecond(M1_IW_KS, M1_IW_KZ);
        ImplicitApplyOp(M1_IW_KZ, M1_IW_KTT);
        Real ts = 0.0, tt2 = 0.0;
        Kokkos::parallel_reduce("m1_impl_bcg_ts", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
          ls += iw_(m,M1_IW_KTT,k,j,i)*iw_(m,M1_IW_KS,k,j,i);
        }, ts);
        Kokkos::parallel_reduce("m1_impl_bcg_tt", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
          ls += iw_(m,M1_IW_KTT,k,j,i)*iw_(m,M1_IW_KTT,k,j,i);
        }, tt2);
        M1GlobalSum2(ts, tt2);
        bcg_nred += 1.0;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real ow = omega;
        par_for("m1_impl_bcg_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
          iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
        });
        rho = rhon;
        rnorm = 0.0;
        Kokkos::parallel_reduce("m1_impl_bcg_rn2", rng,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
          Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
          lmax = (r > lmax) ? r : lmax;
        }, Kokkos::Max<Real>(rnorm));
        M1GlobalMax(rnorm);
        bcg_nred += 1.0;
        if (rnorm/bscale < tol) {
          // the TRUE residual, which is what the tolerance is about
          ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
          par_for("m1_impl_bcg_true", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
          });
          Real rt = 0.0;
          Kokkos::parallel_reduce("m1_impl_bcg_rnt", rng,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
            Real r = fabs(iw_(m,M1_IW_KR,k,j,i));
            lmax = (r > lmax) ? r : lmax;
          }, Kokkos::Max<Real>(rt));
          M1GlobalMax(rt);
          bcg_nred += 1.0;
          if (rt/bscale < tol) {
            done = true;
          } else {
            breakdown = true;   // restart the recurrence from the true residual
          }
        }
        if (!done && !(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        // FALL BACK: one line-Jacobi update, which is always available -- its right-hand
        // side is b minus the lagged off-diagonal term, i.e. the M1_IW_TR the assembly
        // produced, rebuilt here because the preconditioner has overwritten it.
        fell_back = true;
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      par_for("m1_impl_bcg_rs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
      });
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }

  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitBiCGStabEnd
//! \brief the common end of both BiCGStab loops: the line-Jacobi fallback or the copy of
//! the iterate into M1_IW_S2, and the inner-iteration statistics.

void RadiationM1::ImplicitBiCGStabEnd(int nit, bool fell_back) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  if (fell_back) {
    bcg_nfall += 1.0;
    const bool thrd = trans_x3;
    par_for("m1_impl_bcg_lj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real r = iw_(m,M1_IW_KB,k,j,i)
               - iw_(m,M1_IW_CJM,k,j,i)*iw_(m,M1_IW_EP,k,j-1,i)
               - iw_(m,M1_IW_CJP,k,j,i)*iw_(m,M1_IW_EP,k,j+1,i);
      if (thrd) {
        r -= iw_(m,M1_IW_CKM,k,j,i)*iw_(m,M1_IW_EP,k-1,j,i)
             + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
      }
      iw_(m,M1_IW_TR,k,j,i) = r;
    });
    // MILESTONE 3b phase D: the right-hand side of a LINE-JACOBI update is the assembled
    // TR, which under implicit_offdiag = operator still carries -L_off(E^k) -- and the
    // caller has just added L_off(E^k) to b.  Take it out again, so that the fallback is
    // the same update it was before phase D.
    if (od_now == M1_OD_OPERATOR) {
      ImplicitOffDiagOp(M1_IW_EP, M1_IW_TR, -1.0);
    }
    ImplicitTridiagSolve();
  } else {
    par_for("m1_impl_bcg_out", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_S2,k,j,i) = iw_(m,M1_IW_KX,k,j,i);
    });
  }
  bcg_nsolve += 1.0;
  bcg_itsum += static_cast<Real>(nit);
  bcg_itmax = std::max(bcg_itmax, static_cast<Real>(nit));
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabFused
//! \brief implicit_bcg_sync = 1 | 2: the SAME right-preconditioned BiCGStab recurrence as
//! ImplicitBiCGStab (same breakdown, restart, true-residual and fallback logic), with its
//! host synchronisations and kernel launches cut.  On the GPU every reduction into a
//! host scalar blocks the host until the queue is empty, and the next launch then starts
//! on an idle device (measured: ~40 us per blocking reduction).  Per iteration:
//!  * rho_{k+1} = (rhat, r_{k+1}) and max|r_{k+1}| come out of the kernel that updates x
//!    and r, i.e. ONE reduction where the original has three kernels and two syncs
//!    (upd, max|r|, and (rhat,r) at the top of the next iteration);
//!  * (t,s) and (t,t) come out of one kernel;
//!  * the vector updates p and s also stage the preconditioner input M1_IW_TR, which
//!    removes the copy kernel of ImplicitPrecond;
//!  * sync level 2 on ONE rank: rhat.v is reduced into a device scalar and s = r - alpha
//!    v reads alpha from there, so rhat.v does not block; the host learns it from the
//!    (t,s),(t,t) reduction and only then applies the |rhat.v| breakdown test.  If that
//!    test fails, s/z/t are discarded and the recurrence restarts from x, which the
//!    iteration had not touched -- exactly the original branch.
//! Blocking reductions per iteration: 5 kernels / 4 MPI_Allreduce (level 0), 3 / 3
//! (level 1), 2 / 0 (level 2, one rank).  The recurrence is identical in exact
//! arithmetic; the sums are ordered differently, so the iterates agree to round-off only.

int RadiationM1::ImplicitBiCGStabFused(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace> pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  using HRed = M1BcgRed<Kokkos::HostSpace>;
  const bool devrv = (impl_bcg_sync == 2) && (global_variable::nranks == 1);
  const int kf = impl_kfuse;   // implicit_krylov_fuse (needs bcg_sync = 1, pcr)
  if (impl_stencil) {ImplicitStencilBuild();}   // the operator of this pass, once
  if (kf == 3 && impl_kpipe) {return ImplicitBiCGStabPipe(rhsmax);}
  if (kf == 3) {return ImplicitBiCGStabTwo(rhsmax);}
  auto rvd_ = bcg_rvd;
  // implicit_lin_cnorm > 0: the max norm is taken of r_i/(s_i E^k_i), s_i = 1 + SRCB_i
  // the row EXCESS of the M-matrix (diagonal minus the off-diagonal moduli: the
  // transport rows sum to zero, the absorption does not), instead of r_i/max|b|.  For
  // an M-matrix A s >= s componentwise, so |dE_i| <= max_j |r_j|/s_j: the norm bounds
  // the error of E in every cell relative to the local E; the test is r < lin_cnorm.
  const bool cn = (impl_lin_cnorm > 0.0);
  const Real efl = e_floor;

  // x0 = the Picard iterate; r0 = b - A x0, with max|r0| and (r0,r0) in the same kernel
  par_for("m1_impl_bcgf_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  M1BcgVal red;
  Kokkos::parallel_reduce("m1_impl_bcgf_r0", pol,
  KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
    int m, k, j, i;
    M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
    k += ks; j += js; i += is;
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
    v.s0 += r*r;
    Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                            *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
    v.mx = (a > v.mx) ? a : v.mx;
  }, HRed(red));
  M1GlobalBcg(red);
  bcg_nred += 1.0;
  Real rnorm = red.mx;
  Real rhon = red.s0;   // (rhat, r) of the NEXT iteration, always known on entry
  bcg_r0rel = cn ? rnorm : (rnorm/bscale);
  // Eisenstat-Walker (implicit_lin_ew_max > 0): max|r0| IS the nonlinear residual of
  // the Picard iterate in the max norm (the system was re-linearised about it), so the
  // forcing term needs nothing that is not already here.  Off: the fixed test, as is.
  const bool ew = (impl_ew_max > 0.0);
  Real tabs = cn ? impl_lin_cnorm : (tol*bscale);
  if (ew) {
    Real eta = impl_ew_max;
    if (ew_fprev > 0.0) {
      eta = impl_ew_gam*SQR(rnorm/ew_fprev);
      const Real sg = impl_ew_gam*SQR(ew_etaprev);
      if (sg > 0.1) {eta = fmax(eta, sg);}
      eta = fmin(eta, impl_ew_max);
    }
    ew_fprev = rnorm;
    ew_etaprev = eta;
    tabs = fmax(tabs, eta*rnorm);
  }
  auto lin_done = [=](const Real r) {
    return (ew || cn) ? (r < tabs) : (r/bscale < tol);
  };

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      const Real bt = beta, om = omega;
      Real rv = 0.0;
      bool rvdone = false;
      if (kf > 0) {
        // implicit_krylov_fuse: the p update rides in the preconditioner's load phase
        ImplicitPrecondX(-1, M1_IW_KY, 1, bt, om);
        if (kf == 2) {
          ImplicitKrylovHalo(M1_IW_KY);
          Real o4[4];
          ImplicitOpX(M1_IW_KY, M1_IW_KV, 1, o4);
          rv = o4[0];
          Real d2 = 0.0;
          M1GlobalSum2(rv, d2);
          bcg_nred += 1.0;
          rvdone = true;
        } else {
          ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
        }
      } else {
      par_for("m1_impl_bcgf_p", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real p = iw_(m,M1_IW_KR,k,j,i)
                 + bt*(iw_(m,M1_IW_KP,k,j,i) - om*iw_(m,M1_IW_KV,k,j,i));
        iw_(m,M1_IW_KP,k,j,i) = p;
        iw_(m,M1_IW_TR,k,j,i) = p;
      });
      ImplicitPrecond(-1, M1_IW_KY);
      ImplicitApplyOp(M1_IW_KY, M1_IW_KV);
      }
      auto rvf = KOKKOS_LAMBDA(const int idx, Real &ls) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        ls += iw_(m,M1_IW_KRH,k,j,i)*iw_(m,M1_IW_KV,k,j,i);
      };
      if (devrv) {
        // no host sync: alpha stays on the device until the (t,s) reduction
        Kokkos::parallel_reduce("m1_impl_bcgf_rv", pol, rvf,
                                Kokkos::Sum<Real, DevMemSpace>(rvd_));
      } else {
        if (!rvdone) {
        Kokkos::parallel_reduce("m1_impl_bcgf_rv", pol, rvf, rv);
        Real d2 = 0.0;
        M1GlobalSum2(rv, d2);
        bcg_nred += 1.0;
        }
        if (!(fabs(rv) > M1_BCG_EPS)) {
          breakdown = true;
        } else {
          alpha = rhon/rv;
        }
      }
      if (!breakdown) {
        const Real al = alpha, rh = rhon;
        const bool dv = devrv;
        if (kf > 0) {
          ImplicitPrecondX(-1, M1_IW_KZ, 2, al, 0.0);
        } else {
        par_for("m1_impl_bcgf_s", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          Real a = dv ? (rh/rvd_()) : al;
          Real s = iw_(m,M1_IW_KR,k,j,i) - a*iw_(m,M1_IW_KV,k,j,i);
          iw_(m,M1_IW_KS,k,j,i) = s;
          iw_(m,M1_IW_TR,k,j,i) = s;
        });
        ImplicitPrecond(-1, M1_IW_KZ);
        }
        if (kf == 2) {
          ImplicitKrylovHalo(M1_IW_KZ);
          red.s2 = 0.0;
          red.mx = 0.0;
          Real o4[4];
          ImplicitOpX(M1_IW_KZ, M1_IW_KTT, 2, o4);
          red.s0 = o4[0];
          red.s1 = o4[1];
        } else {
        ImplicitApplyOp(M1_IW_KZ, M1_IW_KTT);
        Kokkos::parallel_reduce("m1_impl_bcgf_ts", pol,
        KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
          int m, k, j, i;
          M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
          k += ks; j += js; i += is;
          Real t = iw_(m,M1_IW_KTT,k,j,i);
          v.s0 += t*iw_(m,M1_IW_KS,k,j,i);
          v.s1 += t*t;
          if (dv && idx == 0) {v.s2 += rvd_();}   // carries rhat.v to the host, exactly
        }, HRed(red));
        }
        M1GlobalBcg(red);
        bcg_nred += 1.0;
        if (devrv) {
          rv = red.s2;
          if (!(fabs(rv) > M1_BCG_EPS)) {
            breakdown = true;
          } else {
            alpha = rhon/rv;
          }
        }
      }
      if (!breakdown) {
        const Real ts = red.s0, tt2 = red.s1;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real al = alpha, ow = omega;
        Kokkos::parallel_reduce("m1_impl_bcgf_upd", pol,
        KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
          int m, k, j, i;
          M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
          k += ks; j += js; i += is;
          iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
          Real r = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
          iw_(m,M1_IW_KR,k,j,i) = r;
          v.s0 += iw_(m,M1_IW_KRH,k,j,i)*r;
          Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                                  *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
          v.mx = (a > v.mx) ? a : v.mx;
        }, HRed(red));
        M1GlobalBcg(red);
        bcg_nred += 1.0;
        rho = rhon;
        rhon = red.s0;
        rnorm = red.mx;
        if (lin_done(rnorm)) {
          // the TRUE residual, which is what the tolerance is about
          ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
          Kokkos::parallel_reduce("m1_impl_bcgf_true", pol,
          KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
            int m, k, j, i;
            M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
            k += ks; j += js; i += is;
            Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
            iw_(m,M1_IW_KR,k,j,i) = r;
            Real a = cn ? (fabs(r)/((1.0 + fmax(iw_(m,M1_IW_SRCB,k,j,i), 0.0))
                                    *fmax(iw_(m,M1_IW_EP,k,j,i), efl))) : fabs(r);
            v.mx = (a > v.mx) ? a : v.mx;
          }, HRed(red));
          M1GlobalBcg(red);
          bcg_nred += 1.0;
          if (lin_done(red.mx)) {
            done = true;
          } else {
            breakdown = true;   // restart the recurrence from the true residual
          }
        }
        if (!done && !(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        fell_back = true;   // one line-Jacobi update (ImplicitBiCGStabEnd)
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      Kokkos::parallel_reduce("m1_impl_bcgf_rs", pol,
      KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
        v.s0 += r*r;
      }, HRed(red));
      M1GlobalBcg(red);
      bcg_nred += 1.0;
      rhon = red.s0;
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::ImplicitBiCGStabTwo
//! \brief implicit_krylov_fuse = 3: the right-preconditioned BiCGStab of
//! ImplicitBiCGStabFused with TWO blocking reductions per iteration instead of three.
//!  * (rhat, v) comes out of the operator kernel v = A y, together with max|r| of the
//!    CURRENT r (the one the previous iteration made): the convergence test is taken
//!    one half-iteration late, and on success that half-iteration (p, y, v; x and r are
//!    untouched by it) is discarded;
//!  * (t,s), (t,t) and (rhat,t) come out of the operator kernel t = A z;
//!  * the update x += alpha y + omega z, r = s - omega t is then a plain kernel, and
//!    rho_{k+1} = (rhat, r_{k+1}) = rho_k - alpha (rhat,v) - omega (rhat,t) follows by
//!    recurrence (exact in exact arithmetic, since (rhat,s) = rho_k - alpha (rhat,v)).
//! Breakdown, restart, true-residual and fallback rules are those of the fused loop.
//! Same iterates in exact arithmetic; round-off differs (the recurrence for rho).

int RadiationM1::ImplicitBiCGStabTwo(Real rhsmax) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  const Real tol = impl_lin_tol;
  const Real bscale = fmax(rhsmax, 1.0e-300);
  const int ni = ie - is + 1;
  const int nji = (je - js + 1)*ni;
  const int nkji = (ke - ks + 1)*nji;
  Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256,1>>
      pol(DevExeSpace(), 0, (nmb1 + 1)*nkji);
  using HRed = M1BcgRed<Kokkos::HostSpace>;

  // x0 = the Picard iterate; r0 = b - A x0, with max|r0| and (r0,r0) in the same kernel
  par_for("m1_impl_bcg2_x0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,M1_IW_KX,k,j,i) = iw_(m,M1_IW_EP,k,j,i);
  });
  ImplicitApplyOp(M1_IW_KX, M1_IW_KV);
  M1BcgVal red;
  Kokkos::parallel_reduce("m1_impl_bcg2_r0", pol,
  KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
    int m, k, j, i;
    M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
    k += ks; j += js; i += is;
    Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KV,k,j,i);
    iw_(m,M1_IW_KR,k,j,i) = r;
    iw_(m,M1_IW_KRH,k,j,i) = r;
    iw_(m,M1_IW_KP,k,j,i) = 0.0;
    iw_(m,M1_IW_KV,k,j,i) = 0.0;
    v.s0 += r*r;
    Real a = fabs(r);
    v.mx = (a > v.mx) ? a : v.mx;
  }, HRed(red));
  M1GlobalBcg(red);
  bcg_nred += 1.0;
  Real rnorm = red.mx;
  Real rhon = red.s0;
  bcg_r0rel = rnorm/bscale;
  const bool ew = (impl_ew_max > 0.0);
  Real tabs = tol*bscale;
  if (ew) {
    Real eta = impl_ew_max;
    if (ew_fprev > 0.0) {
      eta = impl_ew_gam*SQR(rnorm/ew_fprev);
      const Real sg = impl_ew_gam*SQR(ew_etaprev);
      if (sg > 0.1) {eta = fmax(eta, sg);}
      eta = fmin(eta, impl_ew_max);
    }
    ew_fprev = rnorm;
    ew_etaprev = eta;
    tabs = fmax(tabs, eta*rnorm);
  }
  auto lin_done = [=](const Real r) {
    return ew ? (r < tabs) : (r/bscale < tol);
  };
  // the true residual of x, which is what the tolerance is about
  auto true_ok = [&]() -> bool {
    ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
    M1BcgVal tr;
    Kokkos::parallel_reduce("m1_impl_bcg2_true", pol,
    KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
      int m, k, j, i;
      M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
      k += ks; j += js; i += is;
      Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
      iw_(m,M1_IW_KR,k,j,i) = r;
      Real a = fabs(r);
      v.mx = (a > v.mx) ? a : v.mx;
    }, HRed(tr));
    M1GlobalBcg(tr);
    bcg_nred += 1.0;
    return lin_done(tr.mx);
  };

  int nit = 0;
  int nrestart = 0;
  Real rho = 1.0, alpha = 1.0, omega = 1.0;
  bool done = lin_done(rnorm);
  bool fell_back = false;
  bool pend = false;   // r was updated by the last iteration; its max is not known yet
  while (!done && nit < impl_lin_maxit) {
    ++nit;
    bool breakdown = !(fabs(rhon) > M1_BCG_EPS) || !(fabs(omega) > M1_BCG_EPS);
    if (!breakdown) {
      Real beta = (rhon/rho)*(alpha/omega);
      ImplicitPrecondX(-1, M1_IW_KY, 1, beta, omega);
      ImplicitKrylovHalo(M1_IW_KY);
      Real o4[4];
      ImplicitOpX(M1_IW_KY, M1_IW_KV, 4, o4);
      M1BcgVal a;
      a.s0 = o4[0]; a.s1 = 0.0; a.s2 = 0.0; a.mx = o4[3];
      M1GlobalBcg(a);
      bcg_nred += 1.0;
      if (pend) {
        pend = false;
        if (lin_done(a.mx)) {
          --nit;   // this half-iteration is discarded: x and r are those it started from
          if (true_ok()) {
            done = true;
            break;
          }
          breakdown = true;   // restart the recurrence from the true residual
        }
      }
      const Real rv = a.s0;
      if (!breakdown && !(fabs(rv) > M1_BCG_EPS)) {breakdown = true;}
      if (!breakdown) {
        alpha = rhon/rv;
        ImplicitPrecondX(-1, M1_IW_KZ, 2, alpha, 0.0);
        ImplicitKrylovHalo(M1_IW_KZ);
        ImplicitOpX(M1_IW_KZ, M1_IW_KTT, 3, o4);
        red.s0 = o4[0]; red.s1 = o4[1]; red.s2 = o4[2]; red.mx = 0.0;
        M1GlobalBcg(red);
        bcg_nred += 1.0;
        const Real ts = red.s0, tt2 = red.s1, rt = red.s2;
        omega = (tt2 > 0.0) ? (ts/tt2) : 0.0;
        const Real al = alpha, ow = omega;
        par_for("m1_impl_bcg2_upd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          iw_(m,M1_IW_KX,k,j,i) += al*iw_(m,M1_IW_KY,k,j,i) + ow*iw_(m,M1_IW_KZ,k,j,i);
          iw_(m,M1_IW_KR,k,j,i) = iw_(m,M1_IW_KS,k,j,i) - ow*iw_(m,M1_IW_KTT,k,j,i);
        });
        rho = rhon;
        rhon = rhon - alpha*rv - omega*rt;
        pend = true;
        if (!(fabs(omega) > M1_BCG_EPS)) {breakdown = true;}
      }
    }
    if (breakdown && !done) {
      pend = false;
      ++nrestart;
      bcg_nbreak += 1.0;
      if (nrestart > 2) {
        fell_back = true;
        break;
      }
      ImplicitApplyOp(M1_IW_KX, M1_IW_KTT);
      Kokkos::parallel_reduce("m1_impl_bcg2_rs", pol,
      KOKKOS_LAMBDA(const int idx, M1BcgVal &v) {
        int m, k, j, i;
        M1BcgIdx(idx, nkji, nji, ni, m, k, j, i);
        k += ks; j += js; i += is;
        Real r = iw_(m,M1_IW_KB,k,j,i) - iw_(m,M1_IW_KTT,k,j,i);
        iw_(m,M1_IW_KR,k,j,i) = r;
        iw_(m,M1_IW_KRH,k,j,i) = r;
        iw_(m,M1_IW_KP,k,j,i) = 0.0;
        iw_(m,M1_IW_KV,k,j,i) = 0.0;
        v.s0 += r*r;
      }, HRed(red));
      M1GlobalBcg(red);
      bcg_nred += 1.0;
      rhon = red.s0;
      rho = 1.0;
      alpha = 1.0;
      omega = 1.0;
    }
  }
  // the cap reached with an update whose max is not known yet: test it once
  if (!done && !fell_back && pend) {
    if (true_ok()) {done = true;}
  }
  ImplicitBiCGStabEnd(nit, fell_back);
  return nit;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::SetImplicitX1BC
//! \brief let a problem generator name the x1 boundary types (and the imposed fluxes) of
//! the implicit solve, for a mesh whose x1 flags are `user` and whose own BC routine the
//! module cannot interpret.  Alternative to <rad_m1>/implicit_bc_x1min|max.

void RadiationM1::SetImplicitX1BC(int lo_type, Real lo_flux, int hi_type, Real hi_flux) {
  ibc_x1min = lo_type;
  iflux_x1min = lo_flux;
  ibc_x1max = hi_type;
  iflux_x1max = hi_flux;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitReport
//! \brief one line at the end of the run with the Picard statistics

void RadiationM1::ImplicitReport() {
  if (transport < M1_TRANSPORT_IMPLICIT_X1) return;
  // the Picard iteration count is MPI_MAX-reduced every step (see ImplicitSolve), so
  // every rank holds the same three numbers and no reduction is needed here
  if (global_variable::my_rank != 0) return;
  Real mean = (impl_nstep > 0.0) ? (impl_itsum/impl_nstep) : 0.0;
  std::cout << "<rad_m1> implicit transport: solves=" << impl_nstep
            << " Picard iterations mean=" << mean << " max=" << impl_itmax
            << " NON-CONVERGED=" << impl_nfail << std::endl;
  if (impl_accel == M1_IACC_ANDERSON) {
    Real apst = (impl_nstep > 0.0) ? (aa_nacc/impl_nstep) : 0.0;
    std::cout << "<rad_m1> anderson: m=" << impl_and_m
              << " beta=" << impl_and_beta
              << " start=" << impl_and_start
              << " accelerated passes=" << aa_nacc
              << " (" << apst << " per solve)"
              << " history restarts=" << aa_nrst << std::endl;
  }
  if (trans_on) {
    Real lmean = (impl_nstep > 0.0) ? (impl_linsum/impl_nstep) : 0.0;
    std::cout << "<rad_m1> implicit transverse ("
              << (bicg_on ? "bicgstab" : "line_jacobi")
              << "): final 7-point linear "
              << "residual mean=" << lmean << " max=" << impl_linmax
              << " tol=" << impl_lin_tol << std::endl;
  }
  if (bicg_on) {
    Real imean = (bcg_nsolve > 0.0) ? (bcg_itsum/bcg_nsolve) : 0.0;
    Real rper = (bcg_itsum > 0.0) ? (bcg_nred/bcg_itsum) : 0.0;
    std::cout << "<rad_m1> bicgstab: outer passes=" << impl_itsum
              << " linear solves=" << bcg_nsolve
              << " inner iterations mean=" << imean << " max=" << bcg_itmax
              << " total=" << bcg_itsum << std::endl;
    std::cout << "<rad_m1> bicgstab: breakdowns=" << bcg_nbreak
              << " line_jacobi fallbacks=" << bcg_nfall
              << " global reductions=" << bcg_nred
              << " (" << rper << " per inner iteration)" << std::endl;
  }
  if (impl_gas_newton || impl_eos_cache) {
    Real fpc = (gas_ncell > 0.0) ? (newt_nfb/gas_ncell) : 0.0;
    Real mpc = (gas_ncell > 0.0) ? (ec_nmiss/gas_ncell) : 0.0;
    std::cout << "<rad_m1> gas coupling: newton="
              << (impl_gas_newton ? "true" : "false")
              << " eos_cache=" << (impl_eos_cache ? "true" : "false")
              << " cell-passes=" << gas_ncell
              << " newton fallbacks=" << newt_nfb << " (" << fpc << " per cell-pass)"
              << std::endl;
    if (impl_eos_cache) {
      std::cout << "<rad_m1> eos_cache: nt=" << impl_ecnt
                << " misses=" << ec_nmiss << " (" << mpc << " per cell-pass)"
                << " max |de|/e vs the table=" << ec_emax
                << " max |dq|/q of the exchanged energy=" << ec_tmax << std::endl;
    }
  }
  if (impl_pcr_check) {
    std::cout << "<rad_m1> line solver=" << ((impl_line_solver == 1) ? "pcr" : "thomas")
              << " pcr_check: calls=" << pcr_chk_n
              << " max over calls of max|x_pcr - x_thomas|/max|x|=" << pcr_chk_max
              << " (rank 0)" << std::endl;
  }
  if (trans_on) {
    std::cout << "<rad_m1> offdiag="
              << ((impl_offdiag == M1_OD_OPERATOR) ? "operator" :
                  ((impl_offdiag == M1_OD_NONE) ? "none" : "lagged"))
              << " closure_relax=" << impl_crelax
              << " closure_lag=" << (impl_clag_step ? "step" : "pass")
              << " positivity fallbacks=" << od_nfall
              << " min E from the solve=" << od_emin << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPicardLog
//! \brief DIAGNOSTIC (<rad_m1>/implicit_picard_log): one line per Picard pass with the
//! E and T parts of the Picard residual (and the x1 index of the cell that owns each),
//! the transverse change lresid, the inner iterations and the inner starting residual.

void RadiationM1::ImplicitPicardLog(int it, int nin, Real resid, Real lresid, bool srct) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  using MaxLoc = Kokkos::MaxLoc<Real,int>;
  Real vals[3] = {0.0, 0.0, 0.0};
  int locs[3] = {-1, -1, -1};
  const int comp[3] = {M1_IW_S1, M1_IW_S3, M1_IW_LRES};
  for (int q = 0; q < 3; ++q) {
    if (q == 1 && !srct) {continue;}
    if (q == 2 && !trans_on) {continue;}
    if (q == 0 && !srct) {
      // without the gas coupling RES is the E part alone
      locs[0] = 0;
    }
    const int c = (q == 0 && !srct) ? M1_IW_RES : comp[q];
    MaxLoc::value_type mloc;
    Kokkos::parallel_reduce("m1_impl_plog",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                  MaxLoc::value_type &lmx) {
      Real r = iw_(m,c,k,j,i);
      if (r > lmx.val) {
        lmx.val = r;
        lmx.loc = i;
      }
    }, MaxLoc(mloc));
    vals[q] = mloc.val;
    locs[q] = mloc.loc;
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> plog step=" << static_cast<int>(impl_nstep) << " pass=" << it
              << " res=" << resid << " resE=" << vals[0] << " iE=" << locs[0]
              << " resT=" << vals[1] << " iT=" << locs[1]
              << " lres=" << lresid << " iL=" << locs[2]
              << " nin=" << nin << " r0=" << bcg_r0rel << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ImplicitSolve
//! \brief the whole backward-Euler step: the Picard loop, the tridiagonal column solves,
//! the write-back into u0 and into the gas.

TaskStatus RadiationM1::ImplicitSolve(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  // closure = vet_sc: the cost of the whole solve is timed against the formal solution
  const bool vetsc = vet_sc;
  auto vc_ = vet_cell;
  // vet_tensor = full: every D_ab of the solve is read from vet_cell (M1DDiag, M1POff)
  auto vd_ = vet_cell;
  const bool dfull = vet_full;
  Kokkos::Timer vtimer;
  if (vetsc) {Kokkos::fence(); vtimer.reset();}

  auto u0_ = u0;
  auto iw_ = iw;
  auto ifw_ = ifw;
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);
  const bool berth = (impl_flux == M1_IFLUX_BERTHON) || (impl_flux == M1_IFLUX_BLEND);
  // MILESTONE 3c: the smooth per-face convex blend.  `bkind` etc. are read only when
  // `blend` is true, and `blend` is false for every 3a/3a2 flux, so those paths keep the
  // arithmetic they had (the weight enters as an exact multiplication by 1.0).
  const bool blend = (impl_flux == M1_IFLUX_BLEND);
  const int bkind = impl_blend, bfm = impl_blend_fmode;
  const bool bdis = blend && (impl_blend_mode == M1_IBMODE_DISSIP);
  const Real btau0 = impl_blend_tau0;
  const Real bflo = impl_blend_flo, bfhi = impl_blend_fhi;
  const bool plmdc = (impl_recon == M1_IRECON_PLMDC);
  const Real rwin = impl_recon_w;
  const Real rfl_ = impl_res_floor;
  const bool rfreeze = impl_recon_freeze;
  const int rnpass = impl_recon_npass;
  auto f0_ = f0x1;
  // F0^n, the face flux at the START of the step: the Picard loop overwrites f0x1 with
  // each new iterate, so the backward-Euler right-hand side needs its own copy
  auto f0n_ = f0x1n;
  Kokkos::deep_copy(DevExeSpace(), f0x1n, f0x1);
  // MILESTONE 3b phase B: the transverse couplings.  `trans` is false for every 3a/3a2/3c
  // configuration and for a 1-D mesh, so those paths keep their arithmetic bit for bit.
  const bool trans = trans_on;
  const bool thrd = trans_x3;
  // MILESTONE 3b phase C: implicit_solver = bicgstab.  False for line_jacobi and for
  // every 3a/3a2/3c configuration, so those paths keep their arithmetic bit for bit.
  const bool bicg = bicg_on;
  // MILESTONE 3b phase D.  Every step STARTS in the configured off-diagonal mode; the
  // positivity fallback below may drop this step to `none` (the operator is not an
  // M-matrix, so E' > 0 is no longer guaranteed by construction).
  od_now = impl_offdiag;
  // the closure under-relaxation and the start-of-step closure freeze.  Both are inert
  // at their defaults (w = 1, lag = pass), so phase C arithmetic is untouched.
  const Real crw = impl_crelax;
  const bool crthin = impl_crelax_thin;
  const bool clagst = impl_clag_step;
  // DIAGNOSTIC dbg_tensor = frozen | tilt: the stored tensor axis survives the step reset
  const bool tpers = (dbg_tensor == 1 || dbg_tensor == 2) && dbg_tensor_init;
  // closure = tau (rad_m1_tau.cpp): (chi, n) from the column optical depth, built on the
  // first pass of the step and read by step (b) on every pass
  const bool tauc = tau_closure;
  if (tauc && !tau_ready) {TauClosureInit();}
  auto tt_ = tau_ten;
  auto f2_ = f0x2;
  auto f3_ = f0x3;
  if (trans) {
    Kokkos::deep_copy(DevExeSpace(), f0x2n, f0x2);
    if (thrd) {Kokkos::deep_copy(DevExeSpace(), f0x3n, f0x3);}
  }
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  auto opac_ = opac;
  auto &mbsize = pmy_pack->pmb->mb_size;
  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  Real ar = arad;
  Real dt = dt_sub;
  bool edd = eddington;
  const int chk = chi_kind;
  bool ovc = source_ovc;
  bool feedback = gas_feedback;
  // MILESTONE 3b phase E: the DEBUG switches of the explicit coupling are honoured here
  // too (they used to be hard-wired to `true` on this branch).  Both default to `true`,
  // so every earlier configuration is bitwise unchanged; what they buy is the ability to
  // ask a multi-D run which HALF of the gas coupling drives a flow -- the radiative FORCE
  // (dbg_gas_force) or the energy exchange (dbg_gas_heat) -- exactly as the 1-D
  // pulsation diagnosis did.  `opac_freeze` is still not on this branch.
  const bool dbgf = dbg_gas_force;
  const bool dbgh = dbg_gas_heat;
  const bool dbgft = dbg_gas_force_trans;
  // LIMIT 3 of the 3a findings.  A physical boundary face hands its whole
  // dt (rho k_t)_f F0_f/c to its ONE interior cell in 3a, so that cell receives 1.5
  // face-shares of radiative force where every interior cell receives 1.0; the residual
  // is a steady force the well-balanced reference a_rad_ref does not carry, and it drives
  // the 10.3 v_MLT bottom-cell flow of I8.  With implicit_bmom_half the boundary face
  // gives HALF, like any other face: the cell-averaged radiative force is then
  // (rho k_t F/c) with F the mean of the cell's two faces, everywhere.  The other half
  // leaves the domain with the radiation, which is where it physically goes.
  const bool bmhalf = impl_bmom_half;
  bool fref = (force_ref == M1_FREF_WB_ARAD);
  auto aref_ = arad_ref;
  Real mq = marshak_q;
  int bclo = ibc_x1min, bchi = ibc_x1max;
  Real fxlo = iflux_x1min, fxhi = iflux_x1max;
  Real eblo = iebath_x1min, ebhi = iebath_x1max;
  bool cyclic = (bclo == M1_IBC_PERIODIC);
  // MILESTONE 3b, LIMIT 4.  With more than one MeshBlock along x1 a block is at a
  // PHYSICAL x1 boundary only when it sits at the corresponding end of its stack; the
  // faces it shares with a stack neighbour are ordinary interior faces whose other cell
  // is a GHOST cell, filled by ImplicitX1Halo with the very numbers the neighbour
  // computed.  With part_nblk == 1 every block is both ends and nothing below changes.
  const int nblkx1 = part_nblk;
  auto pos_ = part_pos;
  const int nlay_ = (part_nblk > 1) ? part_nlay : 0;

  const bool have_hydro = (pmy_pack->phydro != nullptr);
  const bool src_on = have_hydro && coupling && dbgh && !opac_zero;
  // MILESTONE 3g: the gas-radiation energy coupling.  Both are false by default and
  // every branch they guard is then dead, so the pre-3g arithmetic is untouched.
  const bool gnewt = impl_gas_newton && src_on;
  const bool usec = impl_eos_cache && have_hydro;
  const bool gasx = (iw_gas >= 0);
  const int igb = gasx ? (iw_gas + M1_IWG_BK) : 0;
  const int igr = gasx ? (iw_gas + M1_IWG_RK) : 0;
  const int igy = gasx ? (iw_gas + M1_IWG_YR) : 0;
  const int igf = gasx ? (iw_gas + M1_IWG_FB) : 0;
  const int igm = gasx ? (iw_gas + M1_IWG_MS) : 0;
  const int ecnt = impl_ecnt;
  auto ec_ = ecache;
  auto uh = have_hydro ? pmy_pack->phydro->u0 : u0;
  const bool etg = have_hydro ? pmy_pack->phydro->use_etotgrav : false;
  auto phicc = have_hydro ? pmy_pack->phydro->phicc0 : arad_ref;

  //-------------------------------------------------------------------------- start state
  par_for("m1_impl_i0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = fmax(u0_(m,M1_E,k,j,i), efl);
    iw_(m,M1_IW_EN,k,j,i) = e;
    iw_(m,M1_IW_EP,k,j,i) = e;
    iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
    iw_(m,M1_IW_V1,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
    iw_(m,M1_IW_DE0,k,j,i) = 0.0;
    iw_(m,M1_IW_G0,k,j,i) = 0.0;
    iw_(m,M1_IW_TP,k,j,i) = 0.0;
    iw_(m,M1_IW_EGN,k,j,i) = 0.0;
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
    if (gasx) {
      iw_(m,igy,k,j,i) = 0.0;
      iw_(m,igf,k,j,i) = 0.0;
      iw_(m,igm,k,j,i) = 0.0;
    }
    if (trans) {
      iw_(m,M1_IW_V2,k,j,i) = 0.0;
      iw_(m,M1_IW_V3,k,j,i) = 0.0;
      iw_(m,M1_IW_F2,k,j,i) = u0_(m,M1_F2,k,j,i);
      iw_(m,M1_IW_F3,k,j,i) = u0_(m,M1_F3,k,j,i);
      if (!tpers) {
        iw_(m,M1_IW_N1,k,j,i) = 0.0;
        iw_(m,M1_IW_N2,k,j,i) = 0.0;
        iw_(m,M1_IW_N3,k,j,i) = 0.0;
      }
      iw_(m,M1_IW_A2,k,j,i) = 0.0;
      iw_(m,M1_IW_A3,k,j,i) = 0.0;
      iw_(m,M1_IW_TDIA,k,j,i) = 0.0;
      iw_(m,M1_IW_TRHS,k,j,i) = 0.0;
      iw_(m,M1_IW_LRES,k,j,i) = 0.0;
    }
  });

  if (have_hydro) {
    auto eos = pmy_pack->phydro->peos->eos_data;
    par_for("m1_impl_i1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      Real ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                       SQR(uh(m,IM3,k,j,i)))*idd;
      Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      Real eg = uh(m,IEN,k,j,i) - ekin - egrv;
      iw_(m,M1_IW_EGN,k,j,i) = eg;
      iw_(m,M1_IW_TP,k,j,i) = eos.Temperature(dd, fmax(eg, 1.0e-300));
      iw_(m,M1_IW_V1,k,j,i) = uh(m,IM1,k,j,i)*idd;
      if (trans) {
        iw_(m,M1_IW_V2,k,j,i) = uh(m,IM2,k,j,i)*idd;
        iw_(m,M1_IW_V3,k,j,i) = uh(m,IM3,k,j,i)*idd;
      }
    });
    // MILESTONE 3g: the FROZEN-DENSITY e(T) cache.  rho does not move over the step, so
    // the density direction of the tabulated energy surface is collapsed ONCE here and
    // every pass of the Picard loop reads a 1-D cubic in ln T instead of the 2-D table.
    if (usec) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      par_for("m1_impl_ecb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        M1EosCacheBuild(eos, ec_, m, k, j, i, ecnt, uh(m,IDN,k,j,i),
                        iw_(m,M1_IW_TP,k,j,i));
      });
    }
  }

  // The SCALE of the Picard convergence test.  With implicit_res_floor = 0 (the 3a
  // default) the residual is the pure relative change |dE|/E, which in a run with a large
  // dynamic range is dominated by cells many orders below the peak: on gate I6 with
  // implicit_recon = plm_dc the tail cells sit 9 orders under the maximum and keep the
  // reported residual above the tolerance for ever, although the SOLUTION is converged
  // (maxit 30 and maxit 100 give the same amplitude to five digits).  A positive
  // implicit_res_floor scales those cells by the column peak instead,
  // res = |dE|/max(E, implicit_res_floor*max(E)).
  Real emax0 = 0.0;
  if (rfl_ > 0.0) {
    Kokkos::parallel_reduce("m1_impl_emax",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = iw_(m,M1_IW_EN,k,j,i);
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(emax0));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&emax0, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    emax0 = g;}
#endif
  }
  const Real escale = rfl_*emax0;

  // the SCALE the TRUE linear residual is measured against: the max norm of the
  // right-hand side of the 7-point system, which is E^n plus the (small) source terms.
  Real rhsmax = 1.0;
  if (trans) {
    Real rm = 0.0;
    Kokkos::parallel_reduce("m1_impl_rhsmax",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = fabs(iw_(m,M1_IW_EN,k,j,i));
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(rm));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&rm, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    rm = g;}
#endif
    rhsmax = fmax(rm, 1.0e-300);
  }
  Real lresid = 0.0;

  // the x1 ghost layers the partitioned solve reads.  E of the start-of-step state is
  // sent once here; the six lagged quantities and the new E are sent inside the loop.
  // Under transport = implicit the SIX-neighbour exchange of ImplicitTransverseHalo
  // carries all of that and seven quantities more, so the hand-rolled x1 halo is not run.
  if (trans) {
    ImplicitTransverseHalo(M1_NHALO_T);
  } else {
    ImplicitX1Halo(true);
  }

  // closure = vet_sc: the formal solution of the start-of-step state.  (chi, n) are
  // then read by step (b) on every Picard pass: the tensor is lagged by one hydro step.
  if (vetsc) {VetShortChar();}

  // implicit_predictor = step: start the Picard loop from the previous step's implicit
  // increment, scaled by dt/dt_prev.  Only the STARTING POINT moves: E^n (M1_IW_EN),
  // e^n, the EOS cache window and every scale below are those of the step, and the
  // closures it is allowed with do not read the iterate (the Eddington D_ab does not
  // depend on n, vet_sc and tau read their own arrays), so the fixed point is unchanged.
  // It runs AFTER the vet_sc formal solution, which reads T^n (M1_IW_TP) as its source,
  // and the moved E is then sent to the ghost cells.
  const bool pred = impl_pred && (edd || vetsc || tauc) && !(aphll && rfreeze);
  if (impl_pred && (static_cast<int>(ipred.extent(0)) != nmb1 + 1)) {
    pred_ok = false;   // the pack changed size (AMR): start cold
  }
  auto pd_ = ipred;
  if (pred) {
    const bool pok = pred_ok && (pred_dt > 0.0);
    const Real rat = pok ? (dt/pred_dt) : 0.0;
    const bool hh = have_hydro;
    par_for("m1_impl_pred", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real tn = iw_(m,M1_IW_TP,k,j,i);
      pd_(m,2,k,j,i) = tn;
      if (!pok) {return;}
      Real ep = iw_(m,M1_IW_EN,k,j,i) + rat*pd_(m,0,k,j,i);
      if (ep > efl) {iw_(m,M1_IW_EP,k,j,i) = ep;}
      if (hh) {
        Real tp = tn + rat*pd_(m,1,k,j,i);
        if (tp > 0.5*tn && tp < 2.0*tn) {iw_(m,M1_IW_TP,k,j,i) = tp;}
      }
    });
    if (pok) {
      if (trans) {ImplicitTransverseHalo(1);} else {ImplicitX1Halo(true);}
    }
  }

  // MILESTONE 3e: the Anderson histories start empty at every step, and the per-cell
  // scale of the fixed-point vector is frozen at the start-of-step energy (see
  // ImplicitAccelSave).  Nothing here runs under implicit_accel = none.
  const bool accel = (impl_accel != M1_IACC_NONE);
  if (accel) {
    aa_nh = 0;
    aa_head = 0;
    aa_hasp = false;
    aa_fnp = -1.0;
    auto sc_ = aa_sc;
    par_for("m1_acc_scale", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      sc_(m,k,j,i) = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
    });
  }

  //--------------------------------------------------------------------- the Picard loop
  int it = 0;
  Real resid = 0.0;
  bool converged = false;
  // the per-pass log (implicit_picard_log): E and T parts of the Picard residual, the
  // transverse change, the inner iterations and the inner starting residual
  const bool plog = (impl_plog > 0) && (impl_nstep < static_cast<Real>(impl_plog));
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  Real rprev = -1.0;
  for (it = 0; it < impl_maxit && !converged; ++it) {
    int nin = -1;
    // MILESTONE 3e: x_k, the state this pass maps
    if (accel) {ImplicitAccelSave();}
    // the off-diagonal mode of THIS pass (the positivity fallback can change it)
    const int odm = od_now;
    // the closure moves only after the first pass, and not at all under
    // implicit_closure_lag = step
    const bool dorel = (crw < 1.0) && (it > 0);
    const bool dofreeze = clagst && (it > 0);
    // DIAGNOSTIC dbg_tensor (VET scaffolding): tkeep = read the stored tensor; ttau =
    // rebuild it from the optical depth (first pass of every step); ttilt = rotate the
    // axis of the tensor computed on the very first pass of the run
    const int tmode = dbg_tensor;
    const bool tkeep = (tmode != 0) && ((it > 0) || (tmode != 3 && dbg_tensor_init));
    const bool ttau = (tmode == 3) && (it == 0);
    const bool ttilt = (tmode == 2) && (it == 0) && !dbg_tensor_init;
    const Real talp = dbg_tensor_tilt;
    const Real tx2min = pmy_pack->pmesh->mesh_size.x2min;
    const Real tx2len = pmy_pack->pmesh->mesh_size.x2max - tx2min;
    if (tmode == 3 && pmy_pack->pmesh->mesh_indcs.nx1 != indcs.nx1) {
      ImplFatal("<rad_m1>/dbg_tensor = tau needs ONE MeshBlock along x1");
    }
    // (a) optional opacity re-evaluation at the current temperature iterate
    if (impl_opac_update && it > 0 && have_hydro && !opac_zero) {
      int otype = opacity_type;
      Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kscat = kappa_s;
      Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
      M1OpacTab ot = otab;
      par_for("m1_impl_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real d = uh(m,IDN,k,j,i);
        Real t = iw_(m,M1_IW_TP,k,j,i);
        Real op, oe, of, os;
        if (otype == M1_OPAC_TABLE) {
          M1TableOpacities(ot, d, t, op, oe, of, os);
        } else {
          M1Opacities(otype, d, t, kp, kev, kf, kscat, rref, tref, aa, bb, op, oe,
                      of, os);
        }
        opac_(m,M1_OP_P,k,j,i) = d*op;
        opac_(m,M1_OP_E,k,j,i) = d*oe;
        opac_(m,M1_OP_T,k,j,i) = d*(of + os);
        iw_(m,M1_IW_KT,k,j,i) = d*(of + os);
      });
    }

    if (tauc && it == 0) {TauClosureBuild();}
    // (b) the lagged closure, the enthalpy-flux coefficient, de0 and g0
    par_for("m1_impl_lag", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real e = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real f1 = iw_(m,M1_IW_F1,k,j,i);
      Real rf = f1/(cl*e);
      if (rf > 1.0) {rf = 1.0;}
      if (rf < -1.0) {rf = -1.0;}
      Real chi = edd ? (1.0/3.0) : M1Chi(fabs(rf), chk);
      Real v1 = iw_(m,M1_IW_V1,k,j,i);
      Real de0;
      if (!trans) {
        iw_(m,M1_IW_WCHI,k,j,i) = chi;
        iw_(m,M1_IW_ADV,k,j,i) = v1*(1.0 + chi);
        // E0 - E with F_2 = F_3 = 0: P_11 = chi E, P_22 = P_33 = (1-chi) E/2
        Real b1 = v1/cl;
        de0 = ovc ? (-2.0*b1*f1/cl) : (b1*b1*e - 2.0*b1*f1/cl + b1*b1*chi*e);
      } else {
        // MILESTONE 3b phase B: the closure of the MULTI-DIMENSIONAL solve.  The reduced
        // flux is the MAGNITUDE |F|/(c E) and the Eddington tensor is built around the
        // unit flux direction n, D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b.  In
        // 1-D this is the branch above with n = (+-1,0,0), which is why the 1-D mesh is
        // handled there (gate G3) and nothing about implicit_x1 changes.
        Real f2c = iw_(m,M1_IW_F2,k,j,i);
        Real f3c = iw_(m,M1_IW_F3,k,j,i);
        Real fm = sqrt(f1*f1 + f2c*f2c + f3c*f3c);
        Real rfm = fm/(cl*e);
        if (rfm > 1.0) {rfm = 1.0;}
        chi = edd ? (1.0/3.0) : M1Chi(rfm, chk);
        Real ifm = 1.0/fmax(fm, 1.0e-300);
        Real n1 = f1*ifm, n2 = f2c*ifm, n3 = f3c*ifm;
        Real v2 = iw_(m,M1_IW_V2,k,j,i), v3 = iw_(m,M1_IW_V3,k,j,i);
        // MILESTONE 3b phase D: how fast the LAGGED closure is allowed to move.
        //   implicit_closure_lag = step  freezes (chi, n) at the start-of-step state for
        //     the whole step (the Eddington tensor is then explicit in time, as in a VET
        //     code that reuses the previous step's tensor), leaving ONE linear solve plus
        //     the temperature nonlinearity per step;
        //   implicit_closure_relax = w   under-relaxes them between passes, optionally
        //     only where the cell is optically thin (theta > 1/2), which is where the
        //     closure feeds back on the solve through (c dt/dx)^2.
        if (ttau) {
          // column optical depth to the top of the block at j-1, j, j+1 (cell centres)
          Real dx1 = mbsize.d_view(m).dx1;
          Real dx2 = mbsize.d_view(m).dx2;
          Real tm = 0.5*iw_(m,M1_IW_KT,k,j-1,i)*dx1;
          Real tc = 0.5*iw_(m,M1_IW_KT,k,j,i)*dx1;
          Real tq = 0.5*iw_(m,M1_IW_KT,k,j+1,i)*dx1;
          for (int ii = i+1; ii <= ie; ++ii) {
            tm += iw_(m,M1_IW_KT,k,j-1,ii)*dx1;
            tc += iw_(m,M1_IW_KT,k,j,ii)*dx1;
            tq += iw_(m,M1_IW_KT,k,j+1,ii)*dx1;
          }
          // exact grey plane-parallel K/J = (tau + q_inf)/(3 (tau + q(tau))), q = Hopf
          Real qh = 0.710446 - 0.133054*exp(-3.4488*tc);
          chi = (tc + 0.710446)/(3.0*(tc + qh));
          Real g1 = -iw_(m,M1_IW_KT,k,j,i);
          Real g2 = (tq - tm)/(2.0*dx2);
          Real ign = 1.0/fmax(sqrt(g1*g1 + g2*g2), 1.0e-300);
          n1 = -g1*ign;
          n2 = -g2*ign;
          n3 = 0.0;
        } else if (ttilt) {
          Real dx2 = mbsize.d_view(m).dx2;
          Real x2v = mbsize.d_view(m).x2min + (static_cast<Real>(j - js) + 0.5)*dx2;
          Real al = talp*sin(2.0*M_PI*(x2v - tx2min)/tx2len);
          Real ca = cos(al), sa = sin(al);
          Real r1 = ca*n1 - sa*n2, r2 = sa*n1 + ca*n2;
          n1 = r1;
          n2 = r2;
        }
        if (vetsc) {
          chi = vc_(m,M1_VET_CHI,k,j,i);
          n1 = vc_(m,M1_VET_N1,k,j,i);
          n2 = vc_(m,M1_VET_N1+1,k,j,i);
          n3 = vc_(m,M1_VET_N1+2,k,j,i);
        } else if (tkeep) {
          chi = iw_(m,M1_IW_WCHI,k,j,i);
          n1 = iw_(m,M1_IW_N1,k,j,i);
          n2 = iw_(m,M1_IW_N2,k,j,i);
          n3 = iw_(m,M1_IW_N3,k,j,i);
        } else if (dofreeze) {
          chi = iw_(m,M1_IW_WCHI,k,j,i);
          n1 = iw_(m,M1_IW_N1,k,j,i);
          n2 = iw_(m,M1_IW_N2,k,j,i);
          n3 = iw_(m,M1_IW_N3,k,j,i);
        } else if (dorel && (!crthin || (ch*dt*iw_(m,M1_IW_KT,k,j,i) < 1.0))) {
          Real w1 = 1.0 - crw;
          chi = w1*iw_(m,M1_IW_WCHI,k,j,i) + crw*chi;
          n1 = w1*iw_(m,M1_IW_N1,k,j,i) + crw*n1;
          n2 = w1*iw_(m,M1_IW_N2,k,j,i) + crw*n2;
          n3 = w1*iw_(m,M1_IW_N3,k,j,i) + crw*n3;
          Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
          if (nn > 0.0) {
            Real inn = 1.0/nn;
            n1 *= inn;
            n2 *= inn;
            n3 *= inn;
          }
        }
        if (tauc) {
          chi = tt_(m,0,k,j,i);
          n1 = tt_(m,1,k,j,i);
          n2 = tt_(m,2,k,j,i);
          n3 = tt_(m,3,k,j,i);
        }
        iw_(m,M1_IW_WCHI,k,j,i) = chi;
        iw_(m,M1_IW_N1,k,j,i) = n1;
        iw_(m,M1_IW_N2,k,j,i) = n2;
        iw_(m,M1_IW_N3,k,j,i) = n3;
        // a_d = v_d + (v.D)_d, so that the enthalpy flux A_d = v_d E + (v.P)_d = a_d E
        Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2), d33 = M1EddDiag(chi,n3);
        Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
        Real d23 = M1EddOff(chi,n2,n3);
        if (dfull) {
          // vet_tensor = full: the guarded K/J of the formal solution, all six components
          d11 = vd_(m,M1_VET_D11,k,j,i);
          d22 = vd_(m,M1_VET_D11+1,k,j,i);
          d33 = vd_(m,M1_VET_D11+2,k,j,i);
          d12 = vd_(m,M1_VET_D11+3,k,j,i);
          d13 = vd_(m,M1_VET_D11+4,k,j,i);
          d23 = vd_(m,M1_VET_D11+5,k,j,i);
        }
        iw_(m,M1_IW_ADV,k,j,i) = v1 + (v1*d11 + v2*d12 + v3*d13);
        iw_(m,M1_IW_A2,k,j,i) = v2 + (v1*d12 + v2*d22 + v3*d23);
        iw_(m,M1_IW_A3,k,j,i) = v3 + (v1*d13 + v2*d23 + v3*d33);
        // E0 - E to O(beta^2), with the full pressure tensor
        Real b1 = v1/cl, b2 = v2/cl, b3 = v3/cl;
        Real bf = (b1*f1 + b2*f2c + b3*f3c)/cl;
        Real bpb = (b1*b1*d11 + b2*b2*d22 + b3*b3*d33
                    + 2.0*(b1*b2*d12 + b1*b3*d13 + b2*b3*d23))*e;
        Real b2sq = b1*b1 + b2*b2 + b3*b3;
        de0 = ovc ? (-2.0*bf) : (b2sq*e - 2.0*bf + bpb);
      }
      iw_(m,M1_IW_DE0,k,j,i) = de0;
      Real rkev = opac_(m,M1_OP_E,k,j,i);
      Real rkpv = opac_(m,M1_OP_P,k,j,i);
      Real tp = iw_(m,M1_IW_TP,k,j,i);
      Real t2 = tp*tp;
      iw_(m,M1_IW_G0,k,j,i) = rkev*(e + de0) - rkpv*ar*t2*t2;
      // The COMOVING reduced flux of the iterate, which is what the HLL part of the
      // ap_hll flux lags (the HLL acts on F0 only; the enthalpy flux A is added back
      // upwinded, exactly as in the explicit advective split).  The LAB f above still
      // drives the closure chi, as it does in the explicit scheme.
      //
      // It is NOT F0_cell/(c E) with F0_cell the arithmetic mean of the two faces.  That
      // is design risk R4 and it is fatal here: in free streaming the upwind face flux is
      // c E_{i-1}, so the cell mean is c (E_{i-1}+E_i)/2 and the derived f is
      // (1 + E_{i-1}/E_i)/2, i.e. 0.5 rather than 1 on the steep side of a pulse.  The
      // wave speeds then reopen to +-c/sqrt(3), the HLL flux turns CENTRED, and the
      // I6 pulse is flattened to its box mean in one crossing (amplitude ratio 0.0014).
      // Each FACE flux is therefore normalised by the E of the cell it comes FROM, which
      // is exactly 1 for an upwind free-streaming face, and the cell value is the mean of
      // the two face ratios.  On a cold start (f0x1 is zero-initialised and the problem
      // generator's state lives in u0) the cell flux is used instead.
      Real r0;
      Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
      if (fabs(fl) + fabs(fr) > 0.0) {
        int iml = (i > is) ? (i-1)
                  : (cyclic ? ie : ((pos_.d_view(m) > 0) ? (is-1) : is));
        int ipr = (i < ie) ? (i+1)
                  : (cyclic ? is : ((pos_.d_view(m) < nblkx1-1) ? (ie+1) : ie));
        Real eul = fmax((fl > 0.0) ? iw_(m,M1_IW_EP,k,j,iml) : e, efl);
        Real eur = fmax((fr > 0.0) ? e : iw_(m,M1_IW_EP,k,j,ipr), efl);
        r0 = 0.5*(fl/(cl*eul) + fr/(cl*eur));
      } else {
        r0 = (f1 - iw_(m,M1_IW_ADV,k,j,i)*e)/(cl*e);
      }
      if (r0 > 1.0) {r0 = 1.0;}
      if (r0 < -1.0) {r0 = -1.0;}
      iw_(m,M1_IW_RF0,k,j,i) = r0;
    });
    if ((tmode == 1 || tmode == 2) && it == 0) {dbg_tensor_init = true;}

    // the x1 halo of the LAGGED quantities (w, a, g0, v1, the comoving reduced flux and
    // the transport opacity).  Every face of the stack is then assembled by both of its
    // blocks from bit-identical numbers.
    if (trans) {
      // ...and only the M1_NHALO_Q of them a pass can still move when the
      // closure is frozen for the step (see M1HaloCompT).
      ImplicitTransverseHalo(dofreeze ? M1_NHALO_Q : M1_NHALO_T);
      // (b1) the lagged transverse operator: the x2/x3 face fluxes of this iterate, their
      // diagonal contribution to the matrix and their lagged right-hand side.
      ImplicitTransverseTerms(it == 0);
    } else {
      ImplicitX1Halo(false);
    }

    // (b2) the FACE coefficients of the asymptotic-preserving HLL blend.  Nothing here
    // runs under implicit_flux = central, where ifw stays identically zero and the row
    // assembled below is bitwise the 3a one.
    // The deferred correction is by default evaluated ONCE per step, at the start-of-step
    // state (implicit_recon_lag = step).  Recomputing it every Picard pass
    // (= picard) makes the loop a limit cycle: the plm limiter keeps switching on a few
    // cells and the strict tolerance is never reached, at 30 iterations per step against
    // 2, although the answer is the same to five digits.  The correction is a lagged,
    // explicit term in any case, so evaluating it at E^n costs nothing in order.
    const int iter = it;
    const bool doface = aphll && (it == 0 || !rfreeze);
    if (doface) {
      const bool dodg = plmdc && (it == 0
                                  || (!rfreeze && (rnpass <= 0 || it < rnpass)));
      par_for("m1_impl_aphll", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        int ipos = pos_.d_view(m);
        bool phys = (((i == is) && (ipos == 0)) ||
                     ((i == ie+1) && (ipos == nblkx1-1))) && !cyclic;
        if (phys) {
          // a physical boundary face: the flux is IMPOSED there (flux / Marshak /
          // reflect / efix), so there is no Riemann problem and no blend.
          ifw_(m,M1_IFW_AL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCR,k,j,i) = 0.0;
          ifw_(m,M1_IFW_DG,k,j,i) = 0.0;
          return;
        }
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real dx = mbsize.d_view(m).dx1;
        Real rfl = iw_(m,M1_IW_RF0,k,j,im);
        Real rfr = iw_(m,M1_IW_RF0,k,j,ip);
        // closed-form M1 wave speeds of the two LAGGED states (1-D: mu = sign f)
        Real bl, br;
        if (edd) {
          br = ch/sqrt(3.0);
          bl = -br;
        } else {
          Real lml, lpl, lmr, lpr;
          M1WaveSpeeds(fabs(rfl), (rfl >= 0.0) ? 1.0 : -1.0, lml, lpl);
          M1WaveSpeeds(fabs(rfr), (rfr >= 0.0) ? 1.0 : -1.0, lmr, lpr);
          bl = ch*fmin(fmin(lml, lmr), 0.0);
          br = ch*fmax(fmax(lpl, lpr), 0.0);
        }
        // alpha: Bloch et al. (2021) eq. 25 with the (1-f^2) guard and the arithmetic
        // face mean of the CELL optical depth.  lp*lm <= 0, so den >= 1 and alpha <= 1.
        Real tauf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,ip))*dx;
        Real al = 1.0;
        if (tauf > 0.0) {
          Real fbar = 0.5*(fabs(rfl) + fabs(rfr));
          Real guard = fmax(1.0 - fbar*fbar, 0.0);
          Real lp = br/ch, lm = bl/ch;
          Real den = 1.0 - 3.0*tauf*guard*lp*lm/(lp - lm + 1.0e-300);
          al = 1.0/fmax(den, 1.0);
        }
        // F_HLL = [b_R c_h f_L E'_L - b_L c_h f_R E'_R + b_R b_L (E'_R - E'_L)]/(b_R-b_L)
        // is linear in E'.  Split into the ADVECTIVE part (the two physical fluxes) and
        // the DISSIPATION (the jump term), because the two carry different weights:
        //
        //   F = alpha F_adv + alpha^2 F_dis + (1 - alpha) F_diff.
        //
        // The dissipation must carry alpha^2 and not alpha.  At piecewise-constant states
        // -- which is what the MATRIX is built from, whatever implicit_recon says -- the
        // HLL dissipation IS the physical diffusion (b_R b_L dE/(b_R-b_L) -> -c dE/3),
        // so weighting it alpha and adding (1-alpha) F_diff on top counts the diffusion
        // TWICE: measured on gate I1 at tau_cell = 1e3, d(sigma^2)/dt came out 2.0022 x
        // the analytic 2D at every CFL.  alpha^2 ~ 1/tau^2 kills it against F_diff
        // ~ 1/tau and leaves the thin limit (alpha -> 1) exactly the plain HLL flux.
        // This is the `alpha2` form of the explicit scheme (rad_m1_closure.hpp), reached
        // here for the same reason.
        //
        // The two fmax()/fmin() are the M-MATRIX GUARDS.  At alpha = 1 they are provably
        // inactive (the HLL consistency condition b_L <= c_h f <= b_R holds on the M1
        // admissible set), but alpha < 1 rescales the two parts differently and the
        // E'_R coefficient can turn positive; the clamp then drops it to zero, which
        // only makes the face flux more upwind and leaves conservation exact (it is one
        // number per face, used with opposite signs by the two cells).
        Real invb = 1.0/(br - bl + 1.0e-300);
        Real adl = br*ch*rfl*invb;        // E'_L coefficient of F_adv
        Real adr = -bl*ch*rfr*invb;       // E'_R coefficient of F_adv
        Real dk = -br*bl*invb;            // >= 0, the dissipation coefficient
        //
        // implicit_flux = berthon drops F_diff altogether and weights BOTH parts of the
        // HLL flux by alpha, which is what alpha was constructed for: alpha (F_adv +
        // F_dis) is the physical diffusion to first order in 1/tau, the advective part
        // supplying the 1/(0.866 tau) that the dissipation alone is short of.  The
        // F_diff weight is then zero, which is what storing AL = 1 below means.
        Real wdis = berth ? al : (al*al);
        // MILESTONE 3c.  w_f in [0,1] from the LAGGED face quantities: the whole face
        // flux is (1 - w_f) F_central + w_f F_berthon (implicit_blend_mode = flux), or
        // the FULL central flux plus w_f times the HLL DISSIPATION alone
        // (= dissipation).  w_f = 1 with mode = flux reproduces `berthon` bitwise and
        // w_f = 0 reproduces `central` bitwise: the multiplications below are by exactly
        // 1.0 or exactly 0.0.  Every input of the weight lives in iw, which the x1 halo
        // of the partitioned solve already carries.
        Real wf = 1.0;
        if (blend) {
          wf = M1BlendWeight(bkind, bfm, tauf, btau0, rfl, rfr, bflo, bfhi);
        }
        // In `dissipation` mode the advective part of the HLL flux is NOT added (the
        // central flux already carries the transport); only the jump term is.
        Real wadv = bdis ? 0.0 : al;
        Real wdsq = bdis ? al : wdis;
        Real ccl = wf*fmax(wadv*adl + wdsq*dk, 0.0);
        Real ccr = wf*fmin(wadv*adr - wdsq*dk, 0.0);
        // how much of the CENTRAL (face-eliminated) flux the row keeps is 1 - AL.
        ifw_(m,M1_IFW_AL,k,j,i) = bdis ? 0.0 : (berth ? wf : al);
        ifw_(m,M1_IFW_HCL,k,j,i) = ccl;
        ifw_(m,M1_IFW_HCR,k,j,i) = ccr;
        // the plm DEFERRED CORRECTION: the difference between the plm and the dc HLL
        // flux at the PREVIOUS iterate.  It goes to the right-hand side, so the matrix
        // stays the low-order M-matrix.  Both E and the comoving reduced flux are
        // reconstructed, with the same limiter the explicit scheme uses; a face whose
        // 4-cell stencil leaves the block falls back to dc (zero correction).
        Real dg = 0.0;
        if (dodg && al > 0.0) {
          int ilo = is - ((ipos > 0) ? nlay_ : 0);
          int ihi = ie + ((ipos < nblkx1-1) ? nlay_ : 0);
          int imm = (im > ilo) ? (im-1) : (cyclic ? ie : -1);
          int ipp = (ip < ihi) ? (ip+1) : (cyclic ? is : -1);
          if (imm >= 0 && ipp >= 0) {
            Real dum;
            Real elp, erp, flp, frp;
            PLM(iw_(m,M1_IW_EP,k,j,imm), iw_(m,M1_IW_EP,k,j,im),
                iw_(m,M1_IW_EP,k,j,ip), elp, dum);
            PLM(iw_(m,M1_IW_EP,k,j,im), iw_(m,M1_IW_EP,k,j,ip),
                iw_(m,M1_IW_EP,k,j,ipp), dum, erp);
            PLM(iw_(m,M1_IW_RF0,k,j,imm), iw_(m,M1_IW_RF0,k,j,im),
                iw_(m,M1_IW_RF0,k,j,ip), flp, dum);
            PLM(iw_(m,M1_IW_RF0,k,j,im), iw_(m,M1_IW_RF0,k,j,ip),
                iw_(m,M1_IW_RF0,k,j,ipp), dum, frp);
            Real ecl = iw_(m,M1_IW_EP,k,j,im), ecr = iw_(m,M1_IW_EP,k,j,ip);
            // the deferred correction applies to the UPWIND part only, so it carries the
            // same weight w_f the upwind part carries (3c).
            Real gp = wf*(wadv*(br*ch*flp*elp - bl*ch*frp*erp)*invb
                          + wdsq*(-dk)*(erp - elp));
            Real gc = ccl*ecl + ccr*ecr;
            // ADMISSIBILITY of the corrected face flux against the DONOR cell.  The
            // reconstructed face energy may exceed the donor cell's own (plm puts
            // E_i (r-1)/(r+1) on top of E_i for a geometric ratio r), and c times that is
            // then faster than the donor can physically emit: the cell drains below what
            // it receives and, on an exponentially falling background, the drain
            // cascades.  Measured on I6 before this clamp: the 1e-4 background of the
            // free-streaming pulse collapsed onto the floor and the peak grew 7x
            // (amplitude ratio 6.98, Picard never converging).  The low-order flux gc
            // already satisfies this bound, so the clamp never removes the whole
            // correction, only the inadmissible part of it.
            Real gmax = wf*ch*ecl, gmin = -wf*ch*ecr;
            gp = fmin(fmax(gp, gmin), gmax);
            // The DEFERRED-CORRECTION WEIGHT.  A deferred correction is a fixed-point
            // iteration x <- A_low^-1 (b + (A_low - A_high) x), and for advection its
            // contraction factor is ~ 2 nu/(1 + nu) with nu = chat dt/dx: it converges
            // only below CFL ~ 1 and DIVERGES above it (measured: Picard never converges
            // at implicit_cfl = 10 and the I6 pulse amplitude comes out 3.85).  The
            // correction is therefore weighted by w = 1/(1 + nu) unless
            // <rad_m1>/implicit_recon_w names a fixed value.  That makes the contraction
            // factor 2 nu/(1 + nu)^2 <= 1/2 at EVERY CFL, and the fixed point a convex
            // blend of the dc and plm fluxes -- still a monotone flux, second-order where
            // w -> 1 (nu << 1, which is where a propagating front is resolved in time at
            // all) and dc where the step is so long that the front is not resolved.
            Real wdc = (rwin > 0.0) ? rwin : (1.0/(1.0 + ch*dt/dx));
            dg = wdc*(gp - gc);
            // UNDER-RELAXATION across the Picard passes.  The plm limiter keeps switching
            // on a handful of cells and the un-relaxed iteration is a small-amplitude
            // limit cycle that never meets the tolerance (30 passes per step against 2,
            // with the answer already right to five digits).  Averaging with the previous
            // pass leaves the fixed point untouched and breaks the cycle.
            if (iter > 0) {dg = 0.5*(dg + ifw_(m,M1_IFW_DG,k,j,i));}
          }
        }
        ifw_(m,M1_IFW_DG,k,j,i) = dg;
      });
    }

    // (c) the emission/absorption source, linearised in T about the iterate
    if (src_on) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      par_for("m1_impl_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rkpv = opac_(m,M1_OP_P,k,j,i);
        Real rkev = opac_(m,M1_OP_E,k,j,i);
        if (rkpv == 0.0 && rkev == 0.0) {
          iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
          iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
          return;
        }
        Real dd = uh(m,IDN,k,j,i);
        Real tk = iw_(m,M1_IW_TP,k,j,i);
        Real nmiss = 0.0;
        M1EosCached<decltype(eos), decltype(ec_)> thc{eos, ec_, m, k, j, i, ecnt,
                                                      &nmiss};
        M1EosDirect<decltype(eos)> thd{eos};
        Real ee, cv;
        if (usec) {
          thc(dd, tk, ee, cv);
        } else {
          thd(dd, tk, ee, cv);
        }
        Real t3 = tk*tk*tk;
        Real t4 = t3*tk;
        Real de0 = iw_(m,M1_IW_DE0,k,j,i);
        Real bk = dd*cv + 4.0*cl*dt*rkpv*ar*t3;
        Real rk = iw_(m,M1_IW_EGN,k,j,i) - ee - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
        // MILESTONE 3g, the Newton SAFEGUARD.  The temperature this pass starts from was
        // produced by the Newton step of the previous pass, whose linearisation dropped
        // the curvature of e(T) and of T^4.  Here -- where e(T_k) has just been
        // evaluated anyway, so the test is FREE -- the exact nonlinear gas residual
        //   y(T) = rho e(T) + c dt rho kappa_P a T^4 - rho e^n - c dt rho kappa_E E0'
        //        = -(R_k + c dt rho kappa_E E')
        // is compared with the value it had BEFORE that step, at the SAME E' (the solve
        // has not moved E since).  If it did not decrease, the Newton step is discarded
        // and the bracketed root find of the pre-3g scheme is run for this cell.
        if (gnewt) {
          Real yprev = iw_(m,igy,k,j,i);
          if (yprev > 0.0) {
            Real ep = iw_(m,M1_IW_EP,k,j,i);
            Real ynow = fabs(rk + cl*dt*rkev*ep);
            // ...but only where the residual still MEANS something.  Once the cell has
            // converged, y is a difference of numbers that cancel to round-off and it
            // stops decreasing monotonically; without this floor every converged cell
            // buys a bracketed root find in every remaining pass (measured: 13 % of all
            // cell-passes fell back, against 0.6 % with it).
            Real ysc = fmax(fabs(iw_(m,M1_IW_EGN,k,j,i)), cl*dt*rkev*fmax(ep, 0.0));
            if (!(ynow < yprev) && ynow > M1_IMPL_TRTOL*ysc) {
              Real tn = tk;
              bool ok = true;
              if (usec) {
                (void) M1ImplTemperatureT(thc, dd, tk, iw_(m,M1_IW_EGN,k,j,i),
                                          cl*dt*rkpv*ar, cl*dt*rkev*(ep + de0), tn, ok);
              } else {
                (void) M1ImplTemperatureT(thd, dd, tk, iw_(m,M1_IW_EGN,k,j,i),
                                          cl*dt*rkpv*ar, cl*dt*rkev*(ep + de0), tn, ok);
              }
              if (ok && tn > 0.0) {
                tk = tn;
                iw_(m,M1_IW_TP,k,j,i) = tk;
                if (usec) {thc(dd, tk, ee, cv);} else {thd(dd, tk, ee, cv);}
                t3 = tk*tk*tk;
                t4 = t3*tk;
                bk = dd*cv + 4.0*cl*dt*rkpv*ar*t3;
                rk = iw_(m,M1_IW_EGN,k,j,i) - ee - cl*dt*rkpv*ar*t4
                     + cl*dt*rkev*de0;
              }
              iw_(m,igf,k,j,i) += 1.0;
            }
          }
          iw_(m,igb,k,j,i) = bk;
          iw_(m,igr,k,j,i) = rk;
        }
        if (gasx) {
          iw_(m,igm,k,j,i) += nmiss;
        }
        Real emis = dt*ch*rkpv*ar;
        Real kk = (bk > 0.0) ? (emis*4.0*t3*cl*dt*rkev/bk) : 0.0;
        iw_(m,M1_IW_SRCB,k,j,i) = dt*ch*rkev - kk;
        iw_(m,M1_IW_SRCR,k,j,i) = emis*t4 - dt*ch*rkev*de0
                                  + ((bk > 0.0) ? (emis*4.0*t3*rk/bk) : 0.0);
      });
    }

    // (d) assemble the tridiagonal system of every column
    // implicit_enthalpy: the deferred correction of the x1 enthalpy flux (header)
    const int enm = impl_enth;
    const bool enth2 = (enm != M1_IENTH_UPWIND);
    const int ngh = indcs.ng;
    par_for("m1_impl_asm", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      Real nu = dt/dx;
      Real cr = ch/cl;
      int ipos = pos_.d_view(m);
      bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
      Real wi = iw_(m,M1_IW_WCHI,k,j,i);
      // in 1-D the stored closure IS the diagonal Eddington component (P_11/E = chi);
      // in multi-D it is chi and D_11 has to be built from the lagged flux direction.
      if (trans) {wi = M1DDiag(iw_,vd_,dfull,m,0,k,j,i);}
      Real ai = iw_(m,M1_IW_ADV,k,j,i);
      Real vi = iw_(m,M1_IW_V1,k,j,i);
      Real aa = 0.0, bb = 1.0, cc = 0.0;
      Real rr = iw_(m,M1_IW_EN,k,j,i) + iw_(m,M1_IW_SRCR,k,j,i);
      bb += iw_(m,M1_IW_SRCB,k,j,i);
      // MILESTONE 3b phase B: the LINE-JACOBI transverse couplings.  Their diagonal part
      // stays on the diagonal (the 7-point M-matrix), the neighbours' lagged part goes to
      // the right-hand side.
      int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
      if (trans) {
        bb += iw_(m,M1_IW_TDIA,k,j,i);
        rr += iw_(m,M1_IW_TRHS,k,j,i);
        BoundaryFlag q1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
        BoundaryFlag q2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
        BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
        BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
        BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
        BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
        if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
        if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
        if ((q3 == BoundaryFlag::block) || (q3 == BoundaryFlag::periodic)) {jl = js-1;}
        if ((q4 == BoundaryFlag::block) || (q4 == BoundaryFlag::periodic)) {ju = je+1;}
        if (thrd && ((q5 == BoundaryFlag::block) ||
                     (q5 == BoundaryFlag::periodic))) {kl = ks-1;}
        if (thrd && ((q6 == BoundaryFlag::block) ||
                     (q6 == BoundaryFlag::periodic))) {ku = ke+1;}
      }
      Real dx2 = mbsize.d_view(m).dx2;
      Real dx3 = mbsize.d_view(m).dx3;
      // implicit_enthalpy: the x1 ghost layers the plm stencil may read on each side
      int hxl = 0, hxh = 0;
      if (enth2) {
        if (trans) {
          hxl = (il < is) ? ngh : 0;
          hxh = (iu > ie) ? ngh : 0;
        } else {
          hxl = botb ? 0 : nlay_;
          hxh = topb ? 0 : nlay_;
        }
      }

      // ---- face i+1/2
      if (i < ie || cyclic || !topb) {
        int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i+1);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        Real wp = iw_(m,M1_IW_WCHI,k,j,ip);
        if (trans) {wp = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);}
        cc -= nu*df*wp;
        Real vf = 0.5*(vi + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,i) + iw_(m,M1_IW_G0,k,j,ip));
        // the OFF-diagonal Eddington terms of the x1 flux equation, d_2 P_12 + d_3 P_13,
        // fully lagged and therefore a right-hand-side term
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,i,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,vd_,
                             dfull)
                    + M1OffDiv(iw_,m,0,k,j,ip,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,
                               M1_IW_EP,vd_,dfull));
        }
        rr -= nu*cr*om*th*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f - ch*cl*dt*od);
        // the HLL part: its E'_L coefficient is >= 0 (diagonal) and its E'_R coefficient
        // <= 0 (upper off-diagonal), so the blend keeps the M-matrix.
        bb += nu*ifw_(m,M1_IFW_HCL,k,j,i+1);
        cc += nu*ifw_(m,M1_IFW_HCR,k,j,i+1);
        rr -= nu*ifw_(m,M1_IFW_DG,k,j,i+1);
        if (vf > 0.0) {
          bb += nu*cr*ai;
        } else {
          cc += nu*cr*iw_(m,M1_IW_ADV,k,j,ip);
        }
        if (enth2) {
          bool o0, o3;
          int i0 = M1EnthIdx(i-1, is, ie, cyclic, hxl, hxh, o0);
          int i3 = M1EnthIdx(i+2, is, ie, cyclic, hxl, hxh, o3);
          rr -= nu*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,i),
                                 iw_(m,M1_IW_EP,k,j,ip), iw_(m,M1_IW_EP,k,j,i3),
                                 o0 && o3, iw_(m,M1_IW_ADV,k,j,i0), ai,
                                 iw_(m,M1_IW_ADV,k,j,ip), iw_(m,M1_IW_ADV,k,j,i3), vf);
        }
      } else if (bchi == M1_IBC_MARSHAK) {
        bb += nu*ch*mq;
        rr += nu*ch*mq*ebhi;
      } else if (bchi == M1_IBC_FLUX) {
        rr -= nu*cr*fxhi;
      }

      // ---- face i-1/2
      if (i > is || cyclic || !botb) {
        int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        Real wm = iw_(m,M1_IW_WCHI,k,j,im);
        if (trans) {wm = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);}
        aa -= nu*df*wm;
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + vi);
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,i));
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                             vd_,dfull)
                    + M1OffDiv(iw_,m,0,k,j,i,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                               vd_,dfull));
        }
        rr += nu*cr*om*th*(f0n_(m,k,j,i) - ch*dt*vf*g0f - ch*cl*dt*od);
        aa -= nu*ifw_(m,M1_IFW_HCL,k,j,i);
        bb -= nu*ifw_(m,M1_IFW_HCR,k,j,i);
        rr += nu*ifw_(m,M1_IFW_DG,k,j,i);
        if (vf > 0.0) {
          aa -= nu*cr*iw_(m,M1_IW_ADV,k,j,im);
        } else {
          bb -= nu*cr*ai;
        }
        if (enth2) {
          bool o0, o3;
          int i0 = M1EnthIdx(i-2, is, ie, cyclic, hxl, hxh, o0);
          int i3 = M1EnthIdx(i+1, is, ie, cyclic, hxl, hxh, o3);
          rr += nu*cr*M1EnthCorr(enm, iw_(m,M1_IW_EP,k,j,i0), iw_(m,M1_IW_EP,k,j,im),
                                 iw_(m,M1_IW_EP,k,j,i), iw_(m,M1_IW_EP,k,j,i3),
                                 o0 && o3, iw_(m,M1_IW_ADV,k,j,i0),
                                 iw_(m,M1_IW_ADV,k,j,im), ai, iw_(m,M1_IW_ADV,k,j,i3),
                                 vf);
        }
      } else if (bclo == M1_IBC_MARSHAK) {
        bb += nu*ch*mq;
        rr += nu*ch*mq*eblo;
      } else if (bclo == M1_IBC_FLUX) {
        rr += nu*cr*fxlo;
      }

      // a Dirichlet end cell: the whole row is replaced, which keeps the matrix an
      // M-matrix and anchors the level of E (see M1_IBC_EFIX)
      if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                      (i == ie && topb && bchi == M1_IBC_EFIX))) {
        aa = 0.0;
        bb = 1.0;
        cc = 0.0;
        rr = iw_(m,M1_IW_EN,k,j,i);
        // the WHOLE row is replaced, transverse couplings included, or the BiCGStab
        // operator would carry an off-diagonal the preconditioner's row does not have.
        if (bicg) {
          iw_(m,M1_IW_CJM,k,j,i) = 0.0;
          iw_(m,M1_IW_CJP,k,j,i) = 0.0;
          iw_(m,M1_IW_CKM,k,j,i) = 0.0;
          iw_(m,M1_IW_CKP,k,j,i) = 0.0;
        }
      }
      iw_(m,M1_IW_TA,k,j,i) = aa;
      iw_(m,M1_IW_TB,k,j,i) = bb;
      iw_(m,M1_IW_TC,k,j,i) = cc;
      iw_(m,M1_IW_TR,k,j,i) = rr;
    });

    // (e) SOLVE the linear system of this pass.  With implicit_solver = line_jacobi
    // that is ONE x1 line solve with the lagged transverse term already on the
    // right-hand side (the outer Picard loop then IS the Jacobi iteration); with
    // bicgstab the very same system -- the same matrix and the same right-hand side --
    // is solved to implicit_lin_tol by a Krylov iteration preconditioned by that line
    // solve, and the outer loop is left with the nonlinearity alone.
    if (bicg) {
      // b = TR + sum_nb C_nb E^k_nb: undo the move of the lagged off-diagonal term to
      // the right-hand side, so that the operator and the right-hand side describe the
      // same system.  EP still holds E^k here, ghosts included.
      par_for("m1_impl_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real b = iw_(m,M1_IW_TR,k,j,i)
                 + iw_(m,M1_IW_CJM,k,j,i)*iw_(m,M1_IW_EP,k,j-1,i)
                 + iw_(m,M1_IW_CJP,k,j,i)*iw_(m,M1_IW_EP,k,j+1,i);
        if (thrd) {
          b += iw_(m,M1_IW_CKM,k,j,i)*iw_(m,M1_IW_EP,k-1,j,i)
               + iw_(m,M1_IW_CKP,k,j,i)*iw_(m,M1_IW_EP,k+1,j,i);
        }
        iw_(m,M1_IW_KB,k,j,i) = b;
      });
      // MILESTONE 3b phase D.  Under implicit_offdiag = operator the assembly has put
      // -L_off(E^k) on the right-hand side (it is inside TR, through the face fluxes);
      // adding L_off(E^k) back takes it out again, and the operator application adds
      // L_off(x) on the LEFT.  The two changes cancel at x = E^k by construction, so the
      // residual the linear solver measures is the residual of the same system the
      // lagged form measures -- what changes is where the term is solved.
      if (odm == M1_OD_OPERATOR) {
        ImplicitOffDiagOp(M1_IW_EP, M1_IW_KB, 1.0);
      }
      nin = ImplicitBiCGStab(rhsmax);
      if (odm == M1_OD_OPERATOR) {
        // POSITIVITY.  The cross-derivative coefficients have mixed signs, so the
        // 9-/19-point operator is not an M-matrix and E' > 0 is no longer guaranteed.
        // Measure the smallest E the solve produced and, if any cell is non-positive,
        // drop the REST of this step to implicit_offdiag = NONE and count the event.
        // `none` is the M-matrix form that survives: `lagged` is not an alternative
        // here, because lagging these terms in an optically thin cell is exactly what
        // has no fixed point (measured: the seeded He slab blows up in 14 steps with
        // lagged + a frozen closure, and runs with none + a frozen closure).
        Real emin = 1.0e300;
        Kokkos::parallel_reduce("m1_impl_odmin",
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                               {nmb1+1,ke+1,je+1,ie+1}),
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmin) {
          Real r = iw_(m,M1_IW_S2,k,j,i);
          lmin = (r < lmin) ? r : lmin;
        }, Kokkos::Min<Real>(emin));
#if MPI_PARALLEL_ENABLED
        {Real g;
        MPI_Allreduce(&emin, &g, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
        emin = g;}
#endif
        od_emin = std::min(od_emin, emin);
        if (!(emin > 0.0)) {
          od_now = M1_OD_NONE;
          od_nfall += 1.0;
        }
      }
    } else {
      ImplicitTridiagSolve();
    }

    // (f) accept E', solve for T' and measure the Picard residual
    if (src_on) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      par_for("m1_impl_tsolve", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        Real rkpv = opac_(m,M1_OP_P,k,j,i);
        Real rkev = opac_(m,M1_OP_E,k,j,i);
        Real told = iw_(m,M1_IW_TP,k,j,i);
        Real tnew = told;
        if (rkpv > 0.0 || rkev > 0.0) {
          Real dd = uh(m,IDN,k,j,i);
          Real de0 = iw_(m,M1_IW_DE0,k,j,i);
          bool ok = true;
          // MILESTONE 3g.  The gas has ALREADY been eliminated locally to build the row
          // (step (c)): the linearised energy equation is B_k dT = R_k + c dt rho
          // kappa_E E', whose dT is exactly what put -4 a T_k^3 c dt rho kappa_E/B_k on
          // the diagonal and the rest on the right-hand side.  With
          // implicit_gas_newton the SAME relation supplies T', at no table evaluation at
          // all, instead of re-solving the nonlinear equation from scratch in every
          // pass.  It is one Newton step of that equation, so the Picard loop is now a
          // Newton iteration on the coupled (E,T) system, and its fixed point -- where
          // the loop stops, |dT|/T < implicit_tol -- satisfies
          // R_k + c dt rho kappa_E E' = 0, i.e. the EXACT nonlinear backward-Euler gas
          // equation with e(T) and T^4 evaluated (not linearised) at the final T.
          // The elimination only ADDS to the diagonal of the radiation row (the
          // coefficient 4 a T^3 c dt rho kappa_E emis/B_k is >= 0 whenever B_k > 0), so
          // the M-matrix property of sect. 7 is untouched by it.
          bool done = false;
          if (gnewt) {
            Real bk = iw_(m,igb,k,j,i);
            Real rk = iw_(m,igr,k,j,i);
            Real yk = rk + cl*dt*rkev*enew;
            Real dtk = (bk > 0.0) ? (yk/bk) : 0.0;
            if (bk > 0.0 && fabs(dtk) <= M1_NEWT_TRUST*told && (told + dtk) > 0.0) {
              tnew = told + dtk;
              iw_(m,igy,k,j,i) = fabs(yk);
              done = true;
            }
          }
          if (!done) {
            // the pre-3g bracketed root find: also the per-cell FALLBACK of the Newton
            // update (c_v <= 0, a step outside the trust region, a non-positive T).
            if (usec) {
              Real nmiss = 0.0;
              M1EosCached<decltype(eos), decltype(ec_)> thc{eos, ec_, m, k, j, i, ecnt,
                                                            &nmiss};
              (void) M1ImplTemperatureT(thc, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                        cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew,
                                        ok);
              if (gasx) {iw_(m,igm,k,j,i) += nmiss;}
            } else {
              (void) M1ImplTemperature(eos, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                       cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew, ok);
            }
            if (!ok) {tnew = told;}
            if (gnewt) {
              iw_(m,igy,k,j,i) = 0.0;
              iw_(m,igf,k,j,i) += 1.0;
            }
          }
        }
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_TP,k,j,i) = tnew;
        Real re = fabs(enew - eold)/fmax(fmax(fabs(enew), escale), 1.0e-300);
        Real rt = fabs(tnew - told)/fmax(fabs(tnew), 1.0e-300);
        iw_(m,M1_IW_RES,k,j,i) = fmax(re, rt);
        if (plog) {
          iw_(m,M1_IW_S1,k,j,i) = re;
          iw_(m,M1_IW_S3,k,j,i) = rt;
        }
      });
    } else {
      par_for("m1_impl_accept", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_RES,k,j,i) = fabs(enew - eold)
                                 /fmax(fmax(fabs(enew), escale), 1.0e-300);
      });
    }

    // the new iterate's E has to reach the ghost cells before the FACE update, or the
    // two blocks that share a face would build it from different states.
    if (trans) {
      // E is the ONLY halo quantity the pass changed since the exchange
      // above, so one component goes, not the whole list.
      ImplicitTransverseHalo(1);
    } else {
      ImplicitX1Halo(true);
    }

    // (g) the face fluxes of the new iterate, and the derived cell flux
    par_for("m1_impl_face", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      int ipos = pos_.d_view(m);
      bool lo = (i == is) && (ipos == 0);
      bool hi = (i == ie+1) && (ipos == nblkx1-1);
      if ((lo || hi) && !cyclic) {
        int bc = lo ? bclo : bchi;
        int ic = lo ? is : ie;
        Real sgn = lo ? -1.0 : 1.0;
        Real fb = 0.0;
        if (bc == M1_IBC_MARSHAK) {
          fb = sgn*cl*mq*(iw_(m,M1_IW_EP,k,j,ic) - (lo ? eblo : ebhi));
        } else if (bc == M1_IBC_FLUX) {
          fb = lo ? fxlo : fxhi;
        } else if (bc == M1_IBC_EFIX) {
          // the Dirichlet row does not define a face flux; the adjacent interior face is
          // copied so that the cell flux and the momentum deposit stay finite.  The
          // boundary cell is held fixed anyway, so nothing downstream depends on it.
          fb = f0_(m,k,j,lo ? (is+1) : ie);
        }
        f0_(m,k,j,i) = fb;
      } else {
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,ip));
        Real wp = iw_(m,M1_IW_WCHI,k,j,ip);
        Real wm = iw_(m,M1_IW_WCHI,k,j,im);
        if (trans) {
          wp = M1DDiag(iw_,vd_,dfull,m,0,k,j,ip);
          wm = M1DDiag(iw_,vd_,dfull,m,0,k,j,im);
        }
        Real gr = (wp*iw_(m,M1_IW_EP,k,j,ip) - wm*iw_(m,M1_IW_EP,k,j,im))/dx;
        Real od = 0.0;
        if (trans && odm != M1_OD_NONE) {
          Real dx2 = mbsize.d_view(m).dx2;
          Real dx3 = mbsize.d_view(m).dx3;
          int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
          BoundaryFlag q1 = mbbcs.d_view(m,BoundaryFace::inner_x1);
          BoundaryFlag q2 = mbbcs.d_view(m,BoundaryFace::outer_x1);
          BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
          BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
          BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
          BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
          if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
          if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
          if ((q3 == BoundaryFlag::block) || (q3 == BoundaryFlag::periodic)) {jl = js-1;}
          if ((q4 == BoundaryFlag::block) || (q4 == BoundaryFlag::periodic)) {ju = je+1;}
          if (thrd && ((q5 == BoundaryFlag::block) ||
                       (q5 == BoundaryFlag::periodic))) {kl = ks-1;}
          if (thrd && ((q6 == BoundaryFlag::block) ||
                       (q6 == BoundaryFlag::periodic))) {ku = ke+1;}
          od = 0.5*(M1OffDiv(iw_,m,0,k,j,im,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,M1_IW_EP,
                             vd_,dfull)
                    + M1OffDiv(iw_,m,0,k,j,ip,dx,dx2,dx3,thrd,il,iu,jl,ju,kl,ku,
                               M1_IW_EP,vd_,dfull));
        }
        Real fn = th*(f0n_(m,k,j,(i == ie+1 && cyclic) ? is : i)
                      - ch*cl*dt*gr - ch*dt*vf*g0f - ch*cl*dt*od);
        // the SAME blend the row was assembled with, so that the stored comoving face
        // flux (which the restart file carries, which the momentum deposit uses and which
        // the next iterate's reduced flux is built from) is the flux the solve applied.
        // The F0^n MEMORY term sits entirely in the F_diff branch: in the thin limit
        // alpha -> 1 and it must NOT survive, or the lagged flux would fight the upwind
        // HLL flux and the front would be damped exactly as it is under `central`.
        // (3c) the test is on the FLUX FORM, not on al: with implicit_blend_mode =
        // dissipation the central weight is 1 (al = 0) and the HLL part is still there.
        Real al = ifw_(m,M1_IFW_AL,k,j,i);
        if (aphll) {
          Real gh = ifw_(m,M1_IFW_HCL,k,j,i)*iw_(m,M1_IW_EP,k,j,im)
                    + ifw_(m,M1_IFW_HCR,k,j,i)*iw_(m,M1_IW_EP,k,j,ip)
                    + ifw_(m,M1_IFW_DG,k,j,i);
          fn = (1.0 - al)*fn + (cl/ch)*gh;
        }
        f0_(m,k,j,i) = fn;
      }
    });
    if (cyclic) {
      // the wrapped face is stored twice; keep the two copies identical
      par_for("m1_impl_facewrap", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        f0_(m,k,j,ie+1) = f0_(m,k,j,is);
      });
    }
    par_for("m1_impl_f1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_F1,k,j,i) = 0.5*(f0_(m,k,j,i) + f0_(m,k,j,i+1))
                              + iw_(m,M1_IW_ADV,k,j,i)*iw_(m,M1_IW_EP,k,j,i);
      if (trans) {
        Real ep = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_F2,k,j,i) = 0.5*(f2_(m,k,j,i) + f2_(m,k,j+1,i))
                                + iw_(m,M1_IW_A2,k,j,i)*ep;
        if (thrd) {
          iw_(m,M1_IW_F3,k,j,i) = 0.5*(f3_(m,k,j,i) + f3_(m,k+1,j,i))
                                  + iw_(m,M1_IW_A3,k,j,i)*ep;
        }
      }
    });

    // (h) convergence
    resid = 0.0;
    Kokkos::parallel_reduce("m1_impl_res",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = iw_(m,M1_IW_RES,k,j,i);
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(resid));
#if MPI_PARALLEL_ENABLED
    // the convergence test must be GLOBAL: with a partitioned column the ranks would
    // otherwise take different numbers of Picard passes and the gather would deadlock,
    // and even with rank-local columns a per-rank test makes the answer depend on the
    // decomposition.  One MPI_MAX of one double per pass.
    {Real rg;
    MPI_Allreduce(&resid, &rg, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    resid = rg;}
#endif
    // MILESTONE 3b phase B.  The Picard test alone is NOT enough once the transverse
    // couplings are lagged: |dE|/E can stall while the off-diagonal terms are still
    // moving.  The TRUE residual of the full 7-point system is measured separately (see
    // ImplicitTransverseTerms) and both have to be met.
    lresid = 0.0;
    if (trans) {
      Kokkos::parallel_reduce("m1_impl_lres",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                             {nmb1+1,ke+1,je+1,ie+1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
        Real r = iw_(m,M1_IW_LRES,k,j,i);
        lmax = (r > lmax) ? r : lmax;
      }, Kokkos::Max<Real>(lresid));
#if MPI_PARALLEL_ENABLED
      {Real lg;
      MPI_Allreduce(&lresid, &lg, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
      lresid = lg;}
#endif
      lresid /= rhsmax;
    }
    converged = (resid < impl_tol) && (!trans || (lresid < impl_lin_tol));
    if (impl_conv_est || !impl_lres_test) {
      // the Picard test, optionally without the confirming pass: q is the contraction
      // of the last two passes, and q/(1-q) times this change bounds what is left
      bool pc = (resid < impl_tol);
      if (impl_conv_est && !pc && rprev > 0.0) {
        Real q = resid/rprev;
        pc = (q < 0.5) && (resid*q/(1.0 - q) < impl_tol);
      }
      bool lc = !trans || (lresid < impl_lin_tol) || (!impl_lres_test && bicg);
      converged = pc && lc;
    }
    rprev = resid;
    if (plog) {ImplicitPicardLog(it, nin, resid, lresid, src_on);}
    // MILESTONE 3e: ACCELERATE.  Only on a pass that is followed by another one: the
    // state the step ENDS on must be the one the face fluxes of step (g) were built
    // from, so a converged pass -- and the last pass of a non-converged step -- keeps
    // the plain Picard iterate.  The accelerated E then has to reach the ghost cells
    // before the next pass assembles a row from it.
    if (accel && !converged && (it + 1 < impl_maxit)) {
      ImplicitAccelApply(it);
      if (trans) {
        // the acceleration rewrites the cell flux F1 as well as E
        ImplicitTransverseHalo(M1_NHALO_T);
      } else {
        ImplicitX1Halo(true);
      }
    }
  }
  if (pred) {
    const bool hh = have_hydro;
    par_for("m1_impl_pstore", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pd_(m,0,k,j,i) = iw_(m,M1_IW_EP,k,j,i) - iw_(m,M1_IW_EN,k,j,i);
      pd_(m,1,k,j,i) = hh ? (iw_(m,M1_IW_TP,k,j,i) - pd_(m,2,k,j,i)) : 0.0;
    });
    pred_ok = true;
    pred_dt = dt;
  }
  if (trans) {
    // refresh the stored transverse face fluxes with the CONVERGED E, so that the
    // persistent state the restart file carries, the momentum deposit and the derived
    // F_2/F_3 are the fluxes the solve ended on (exactly what step (g) does for x1).
    ImplicitTransverseTerms(false);
    impl_linmax = std::max(impl_linmax, lresid);
    impl_linsum += lresid;
  }
#if MPI_PARALLEL_ENABLED
  {
    // a column never crosses a rank here, but the ITERATION COUNT must be global or the
    // ranks would take different numbers of Picard passes and diverge
    int ilocal = it, iglob;
    MPI_Allreduce(&ilocal, &iglob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    it = iglob;
  }
#endif
  impl_nstep += 1.0;
  impl_itsum += static_cast<Real>(it);
  impl_itmax = std::max(impl_itmax, static_cast<Real>(it));
  if (!converged) {
    impl_nfail += 1.0;
    // MILESTONE 3e: WHERE the outer iteration stalled.  One line per non-converged step
    // with the cell that owns the largest Picard residual -- the depth index i is what
    // says whether the stall is in the optically thin top or at the base.  Nothing here
    // changes any number; a step that converges (every step of every earlier gate)
    // prints nothing.
    using MaxLoc = Kokkos::MaxLoc<Real,int>;
    MaxLoc::value_type mloc;
    const int nj = je - js + 1, ni = ie - is + 1;
    Kokkos::parallel_reduce("m1_impl_resloc",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                  MaxLoc::value_type &lmx) {
      Real r = iw_(m,M1_IW_RES,k,j,i);
      if (r > lmx.val) {
        lmx.val = r;
        lmx.loc = ((m*(ke-ks+1) + (k-ks))*nj + (j-js))*ni + (i-is);
      }
    }, MaxLoc(mloc));
    int lmb = mloc.loc/((ke-ks+1)*nj*ni);
    int lrem = mloc.loc - lmb*(ke-ks+1)*nj*ni;
    int lk = lrem/(nj*ni);
    int lj = (lrem - lk*nj*ni)/ni;
    int li = lrem - lk*nj*ni - lj*ni;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> Picard NON-CONVERGED after " << it << " passes:"
                << " resid=" << resid << " (tol " << impl_tol << ")"
                << " lin_resid=" << lresid
                << " worst cell of rank 0 (m,k,j,i)=(" << lmb << "," << (lk+ks) << ","
                << (lj+js) << "," << (li+is) << ")" << std::endl;
    }
  }

  // MILESTONE 3g: the per-cell counters of the gas coupling, reduced ONCE per step.
  if (gasx) {
    Real sfb = 0.0, sms = 0.0;
    Kokkos::parallel_reduce("m1_impl_gcnt",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lf) {
      lf += iw_(m,igf,k,j,i);
    }, Kokkos::Sum<Real>(sfb));
    Kokkos::parallel_reduce("m1_impl_gms",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ls) {
      ls += iw_(m,igm,k,j,i);
    }, Kokkos::Sum<Real>(sms));
    Real ncell = static_cast<Real>(nmb1+1)*static_cast<Real>(ke-ks+1)
                 *static_cast<Real>(je-js+1)*static_cast<Real>(ie-is+1);
#if MPI_PARALLEL_ENABLED
    {Real lo[3] = {sfb, sms, ncell}, gl[3];
    MPI_Allreduce(lo, gl, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    sfb = gl[0];
    sms = gl[1];
    ncell = gl[2];}
#endif
    newt_nfb += sfb;
    ec_nmiss += sms;
    gas_ncell += ncell*static_cast<Real>(it);
  }

  // MILESTONE 3g: ONE TRUE-TABLE evaluation per cell at the state the step ends on, to
  // MEASURE what the cache cost.  Two numbers are taken: the relative error of e(T')
  // itself, and the relative error of the energy q = SRCR - SRCB E' the gas actually
  // exchanges with the radiation -- the only channel through which the cache can reach
  // the evolved state.
  //
  // Nothing is CORRECTED here, deliberately.  The gas energy is set from the ASSEMBLED
  // row (see the write-back below), which is what makes e_gas + (c/chat) E change by the
  // fluxes and the work term alone to round-off; overwriting q after the solve with a
  // re-evaluated one would break exactly that algebraic balance.  Consistency with the
  // real table is instead structural: the evolved gas state is (rho, e_gas), T' is not
  // persistent, and the next step re-inverts e_gas through the real table.
  if (usec && impl_eccheck && src_on) {
    auto eos = pmy_pack->phydro->peos->eos_data;
    Real emx = 0.0, qmx = 0.0;
    Kokkos::parallel_reduce("m1_impl_eck",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      iw_(m,igm,k,j,i) = 0.0;
      Real rkpv = opac_(m,M1_OP_P,k,j,i);
      Real rkev = opac_(m,M1_OP_E,k,j,i);
      if (rkpv == 0.0 && rkev == 0.0) {return;}
      Real dd = uh(m,IDN,k,j,i);
      Real tk = iw_(m,M1_IW_TP,k,j,i);
      Real ce, ccv;
      if (!M1EosCacheEval(eos, ec_, m, k, j, i, ecnt, dd, tk, ce, ccv)) {return;}
      Real te, pp, cr, ct, tcv;
      eos.ThermoAt(dd, tk, te, pp, cr, ct, tcv);
      lmax = fmax(lmax, fabs(ce - te)/fmax(fabs(te), 1.0e-300));
      // the same row the pass assembled, re-made with the TRUE table
      Real t3 = tk*tk*tk, t4 = t3*tk;
      Real de0 = iw_(m,M1_IW_DE0,k,j,i);
      Real ep = iw_(m,M1_IW_EP,k,j,i);
      Real emis = dt*ch*rkpv*ar;
      Real bkc = dd*ccv + 4.0*cl*dt*rkpv*ar*t3;
      Real bkt = dd*tcv + 4.0*cl*dt*rkpv*ar*t3;
      Real rkc = iw_(m,M1_IW_EGN,k,j,i) - ce - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
      Real rkt = iw_(m,M1_IW_EGN,k,j,i) - te - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
      Real qc = emis*t4 - dt*ch*rkev*de0 - (dt*ch*rkev)*ep;
      Real qt = qc;
      if (bkc > 0.0) {qc += emis*4.0*t3*(rkc + cl*dt*rkev*ep)/bkc;}
      if (bkt > 0.0) {qt += emis*4.0*t3*(rkt + cl*dt*rkev*ep)/bkt;}
      iw_(m,igm,k,j,i) = fabs(qt - qc)/fmax(fabs(qt), 1.0e-300);
    }, Kokkos::Max<Real>(emx));
    Kokkos::parallel_reduce("m1_impl_eckq",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      lmax = fmax(lmax, iw_(m,igm,k,j,i));
    }, Kokkos::Max<Real>(qmx));
#if MPI_PARALLEL_ENABLED
    {Real lo[2] = {emx, qmx}, gl[2];
    MPI_Allreduce(lo, gl, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    emx = gl[0];
    qmx = gl[1];}
#endif
    ec_emax = std::max(ec_emax, emx);
    ec_tmax = std::max(ec_tmax, qmx);
  }

  //------------------------------------------------------------------------- write back
  par_for("m1_impl_wb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real ep = iw_(m,M1_IW_EP,k,j,i);
    Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
    Real fp1 = 0.5*(fl + fr) + iw_(m,M1_IW_ADV,k,j,i)*ep;
    // MILESTONE 3b phase B: the derived cell-centred transverse fluxes, the face means
    // plus the enthalpy flux, exactly as in x1.
    Real fp2 = 0.0, fp3 = 0.0;
    Real g2l = 0.0, g2r = 0.0, g3l = 0.0, g3r = 0.0;
    if (trans) {
      g2l = f2_(m,k,j,i);
      g2r = f2_(m,k,j+1,i);
      fp2 = 0.5*(g2l + g2r) + iw_(m,M1_IW_A2,k,j,i)*ep;
      if (thrd) {
        g3l = f3_(m,k,j,i);
        g3r = f3_(m,k+1,j,i);
        fp3 = 0.5*(g3l + g3r) + iw_(m,M1_IW_A3,k,j,i)*ep;
      }
    }

    Real work = 0.0, dm1 = 0.0, dmref = 0.0, eg = 0.0, ekin = 0.0, egrv = 0.0;
    Real dm2 = 0.0, dm3 = 0.0;
    Real dd = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
    if (have_hydro) {
      dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      v1 = uh(m,IM1,k,j,i)*idd;
      v2 = uh(m,IM2,k,j,i)*idd;
      v3 = uh(m,IM3,k,j,i)*idd;
      ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                  SQR(uh(m,IM3,k,j,i)))*idd;
      egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      eg = iw_(m,M1_IW_EGN,k,j,i);
      // (a) the energy the RADIATION gained from the gas over the step.  It is taken
      // from the ASSEMBLED row, q = SRCR - SRCB E', and not from rho kappa_P a T'^4:
      // the two differ by the Picard remainder of the linearisation, and only the first
      // is the amount the solved E actually received.  Setting the gas energy from it
      // makes e_gas + (c/chat) E change by the face fluxes and the work term ALONE, to
      // round-off -- measured: the T'^4 form drifted 2.8e-11 over 2000 steps of T5, this
      // one 0.  At convergence the two agree, so T' stays the consistent temperature.
      if (coupling && dbgh) {
        Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*ep;
        eg -= (cl/ch)*qq;
      }
      // (b) MOMENTUM.  Each x1 face hands dt (rho k_t)_f F0'_f/c to the gas, half to
      // each of its two cells (a physical boundary face gives all of it to its one
      // interior cell), which is exactly what the implicit face source removed from the
      // radiation: sum_cells dm = sum_faces dt (rho k_t)_f F0'_f/c.
      if (coupling && dbgf) {
        Real ktl, ktr;
        Real wl = 0.5, wr = 0.5;
        int ipos = pos_.d_view(m);
        if (i == is && ipos == 0 && !cyclic) {
          ktl = iw_(m,M1_IW_KT,k,j,i);
          wl = bmhalf ? 0.5 : 1.0;
        } else {
          int im = (cyclic && i == is) ? ie : (i-1);
          ktl = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
        }
        if (i == ie && ipos == nblkx1-1 && !cyclic) {
          ktr = iw_(m,M1_IW_KT,k,j,i);
          wr = bmhalf ? 0.5 : 1.0;
        } else {
          int ip = (cyclic && i == ie) ? is : (i+1);
          ktr = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
        }
        dm1 = (dt/cl)*(wl*ktl*fl + wr*ktr*fr);
        dmref = fref ? (dt*dd*aref_(m,k,j,i)) : 0.0;
        if (trans && dbgft) {
          // the same rule per transverse direction: each face hands
          // dt (rho k_t)_f F0_f/c to the gas, half to each of its two cells, and a
          // PHYSICAL boundary face (where F0 is zero anyway) half as well under
          // implicit_bmom_half.
          BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
          BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
          bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
          bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
          Real ktc = iw_(m,M1_IW_KT,k,j,i);
          Real kl2 = (j == js && p2lo) ? ktc
                     : 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + ktc);
          Real kr2 = (j == je && p2hi) ? ktc
                     : 0.5*(ktc + iw_(m,M1_IW_KT,k,j+1,i));
          Real u2 = ((j == js && p2lo) && !bmhalf) ? 1.0 : 0.5;
          Real w2 = ((j == je && p2hi) && !bmhalf) ? 1.0 : 0.5;
          dm2 = (dt/cl)*(u2*kl2*g2l + w2*kr2*g2r);
          if (thrd) {
            BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
            BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
            bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
            bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
            Real kl3 = (k == ks && p3lo) ? ktc
                       : 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + ktc);
            Real kr3 = (k == ke && p3hi) ? ktc
                       : 0.5*(ktc + iw_(m,M1_IW_KT,k+1,j,i));
            Real u3 = ((k == ks && p3lo) && !bmhalf) ? 1.0 : 0.5;
            Real w3 = ((k == ke && p3hi) && !bmhalf) ? 1.0 : 0.5;
            dm3 = (dt/cl)*(u3*kl3*g3l + w3*kr3*g3r);
          }
        }
        if (feedback) {
          Real idg = 1.0/fmax(dd, 1.0e-300);
          Real w1 = (uh(m,IM1,k,j,i) + dm1)*idg;
          work = 0.5*(v1 + w1)*dm1;
          if (trans && dbgft) {
            Real w2n = (uh(m,IM2,k,j,i) + dm2)*idg;
            work += 0.5*(v2 + w2n)*dm2;
            if (thrd) {
              Real w3n = (uh(m,IM3,k,j,i) + dm3)*idg;
              work += 0.5*(v3 + w3n)*dm3;
            }
          }
          ep -= (ch/cl)*work;
        }
      }
    }

    M1ApplyLimits(cl, efl, ep, fp1, fp2, fp3);
    u0_(m,M1_E,k,j,i) = ep;
    u0_(m,M1_F1,k,j,i) = fp1;
    u0_(m,M1_F2,k,j,i) = fp2;
    u0_(m,M1_F3,k,j,i) = fp3;
    if (have_hydro && feedback) {
      uh(m,IM1,k,j,i) = uh(m,IM1,k,j,i) + dm1 - dmref;
      if (trans && dbgft) {
        uh(m,IM2,k,j,i) = uh(m,IM2,k,j,i) + dm2;
        if (thrd) {uh(m,IM3,k,j,i) = uh(m,IM3,k,j,i) + dm3;}
      }
      uh(m,IEN,k,j,i) = eg + ekin + egrv + work;
    }
  });

  if (vetsc) {Kokkos::fence(); vet_itime += vtimer.seconds();}
  return TaskStatus::complete;
}

} // namespace radm1
