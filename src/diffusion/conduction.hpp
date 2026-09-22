#ifndef DIFFUSION_CONDUCTION_HPP_
#define DIFFUSION_CONDUCTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.hpp
//! \brief Contains data and functions that implement various formulations for conduction.
//  Currently only isotropic conduction implemented

#include <cstdint>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "parameter_input.hpp"

// the dedicated one-variable ghost exchange the implicit transverse solve drives
class MeshBoundaryValuesCC;

//----------------------------------------------------------------------------------------
//! \class Conduction
//! \brief data and functions that implement thermal conduction in Hydro and MHD

class Conduction {
 public:
  Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~Conduction();

  // data
  Real dtnew;
  // the cell that set dtnew, for the collapse report in Mesh::NewTimeStep
  int dtnew_m = -1, dtnew_k = -1, dtnew_j = -1, dtnew_i = -1;
  // ...and its state, so the report says WHY that cell is slow rather than only where
  // it is.  Filled by one extra single-cell kernel, and only when dtnew has just
  // collapsed (or on the first call), so it costs nothing in a healthy run.
  static constexpr int ndtdiag = 16;
  DualArray1D<Real> dt_diag;
  bool dt_diag_valid = false;
  Real dtnew_prev = -1.0;
  std::string iso_cond_type; // "constant", "spitzer", "spitzer_limited", "radiative"
  // radiative conductivity kappa_rad = 16 sigma T^3/(3 kappa_R rho), kappa_R the
  // Freedman+2014 Rosseland mean (utils/rosseland.hpp), applied on faces whose pressure
  // is at or above rad_pcut (the layers the two-stream RT does not reach), flux-limited
  // to
  // sigma T^4, with an imposed heat flux through the inner x1 wall (the planet's internal
  // flux). rad_kappa_fac scales kappa_R (tests).
  Real rad_met = 0.0;          // [M/H] for the Rosseland fit
  Real rad_pcut = 0.0;         // code units; faces with p < rad_pcut get no flux
  Real rad_flux_inner = 0.0;   // code units; heat flux through the inner x1 wall
  Real rad_kappa_fac = 1.0;
  // <problem>/nan_report: record the first x1 face whose radiative heat flux comes out
  // non-finite, with the inputs that made it.  Default false; nothing runs when off.
  bool nan_report = false;
  bool rad_flux_limit = true;
  // The FREE-STREAMING FLUX the limiter saturates at.  The limiter is
  //     F = -K grad T / sqrt(1 + (K grad T/F_free)^2),
  // and F_free must be the largest flux the radiation field can carry, which for an
  // isotropic LTE field of energy density E = a T^4 streaming freely is
  //     F_free = c E = c a T^4 = 4 sigma T^4.
  // Every limiter in this module (the explicit face flux, RadFaceKCode, the implicit
  // radial solve, and both timestep estimates) uses this factor of 4, so they cannot
  // drift.
  // OPTICAL-DEPTH BLEND (rad_tau_hi > 0): instead of the pressure cut, each x1 face
  // carries a weight w(tau_R) rising smoothly from 0 at rad_tau_lo to 1 at rad_tau_hi,
  // tau_R the column's Rosseland depth from the top; the diffusion flux is multiplied
  // by w and a two-stream RT that reads rad_w multiplies its own by 1 - w, so the two
  // operators overlap and hand over conservatively.  rad_w and rad_tauf sit on x1 faces.
  Real rad_tau_lo = 0.0, rad_tau_hi = 0.0;
  bool rad_tau_mode = false;
  // rad_tmax_kappa (0 = off): a CEILING on the temperature that enters kappa_rad, both
  // its T^3 and the kappa_R lookup.  kappa_rad ~ T^3/(kappa_R rho) explodes when a
  // single drained cell runs away to 1e5-1e7 K, and the explicit diffusive dt collapses
  // with it -- conduction reports the runaway rather than causing it, and capping T here
  // keeps the timestep survivable while the column refills.  The temperature GRADIENT
  // and the flux limiter's sigma T^4 are left at the true temperature.
  Real rad_tmax = 0.0;         // K; <= 0 disables the cap
  // rad_blend_radial (default true): the blend weight applies to the x1 faces too.
  // False keeps the radial diffusion at full weight everywhere and applies w only to
  // the angular faces -- for a SELF-LUMINOUS object with no two-stream above the
  // blend, where flux-limited radial diffusion is the surface treatment and the
  // horizontal exchange in the optically thin layers is what is switched off.
  bool rad_blend_radial = true;
  // rad_blend_transverse (default true): the blend weight also multiplies the x2/x3
  // (TRANSVERSE) face conductances -- of the explicit flux, of the rad_cap_ang cap and
  // of the super-time-stepped operator's coefficients alike, as a cell-centred 4-point
  // average 0.25*(wf(i)+wf(i+1)) of the two cells sharing the face.  True is the
  // historical behaviour and is bitwise unchanged.
  //
  // SET IT FALSE WHENEVER THE GREY TWO-STREAM OWNS THE WHOLE COLUMN.  In that
  // configuration (rad_tau_lo/rad_tau_hi pushed below the bottom of the box so that
  // w = 0 on every x1 face, with problem/rt_bottom_flux imposing rad_flux_inner on the
  // sweep's own lower boundary) the radial diffusion is deliberately inert -- but the
  // two-stream is a plane-parallel, column-by-column solver and carries NO horizontal
  // radiative transport at all.  With the weight applied transversely as well, w = 0
  // would switch the transverse operator off everywhere too, and the box would have no
  // horizontal radiative exchange of any kind, which is precisely the transport this
  // operator exists to supply.  False decouples the transverse conductance from the
  // vertical blend: the x2/x3 faces use weight 1 while the x1 faces keep w
  // (rad_blend_radial is untouched, and so is cap_c1 under rad_sts_all).  The flux
  // limiter and the rad_gate_rho density gate still apply on every transverse face.
  bool rad_blend_transverse = true;

  // rad_tr_tau_lo / rad_tr_tau_hi (default 0 = OFF, and bitwise inert then): the SMOOTH
  // TRANSVERSE TAPER.  Every x2/x3 face conductance -- the frozen cap_c2/cap_c3 of
  // BuildAngularCoeffs, the explicit x2/x3 face flux, and the transverse conduction dt --
  // is multiplied by a weight that rises from 0 at rad_tr_tau_lo to 1 at rad_tr_tau_hi in
  // the COLUMN optical depth tau measured down from the top face (the rad_tauf array
  // BuildRadWeights fills, face-averaged over the four x1 faces of the transverse face,
  // exactly as the vertical blend averages its weight there).  It is the same raised
  // cosine in log tau the vertical blend uses (RadBlendWeight), so the two tapers cannot
  // disagree about the shape of a handover.
  //
  // WHY.  The transverse operator is a DIFFUSION approximation of horizontal radiative
  // exchange, and that approximation is invalid exactly where tau < 1: a photon there
  // crosses the box horizontally without being absorbed, so there is no local flux
  // -K grad_h T to speak of.  Those same planes are also the stiffest -- K ~ 1/(kappa
  // rho) diverges as the gas thins -- so they set the RKL1 substage count of the whole
  // operator (145-197 substages against ~15 for the deep planes in the B-star box).
  // Tapering them off removes physics that was never right and buys most of the cost.
  //
  // The DENSITY gate rad_gate_rho is INDEPENDENT of this and stays: it is a safety
  // against the artificial low-density medium, selected by what the gas is rather than by
  // how deep it sits, and the two multiply.
  //
  // PRODUCTION (bstar_fecz/prod_whole): rad_tr_tau_lo = 3, rad_tr_tau_hi = 10.  The
  // photosphere sits at tau = 2/3, i.e. at weight 0, so no horizontal radiative diffusion
  // survives at or above the surface; the operator is at full strength by tau = 10.
  Real rad_tr_tau_lo = 0.0, rad_tr_tau_hi = 0.0;
  // rad_kappa_src = freedman (default) | table: with table, kappa_R(T,p) is a bilinear
  // lookup of log10 kappa_R over (log10 T, log10 p[cgs]) in a table the problem
  // generator hands over ONCE at start-up (deep_hot_jupiter_rt tabulates the Rosseland
  // mean of its own correlated-k table + continuum on that table's grid, so the
  // diffusion and the two-stream share one opacity). Same cost as the Freedman fit.
  // Until the table is set (rad_kr_nT == 0) the Freedman fit is used.
  bool rad_kappa_tab = false;
  // rad_kappa_src = table_rho: the same lookup with log10 rho [g/cm^3] as the second
  // axis instead of log10 p.  Stellar opacity tables (OPAL/OPLIB, AESOPUS, Ferguson) are
  // tabulated in (T, rho) -- via logR = log rho - 3 log T6 -- and a (T, p) axis would
  // need the equation of state to convert, on the table's grid, once per node.  With
  // this mode the pgen hands over (T, rho) directly and every lookup site already has
  // the face or cell density.  rad_kr_lP then holds log10 rho.
  bool rad_kappa_rho = false;
  // rad_kappa_rmax > 0 [code length]: a RADIATIVELY INERT region above this radius.  The
  // red-giant runs put a hot hydrostatic corona above the star, and that corona must
  // neither cool nor set the timestep: it is a numerical lid, not a stellar layer, and a
  // 6e5 K gas at 1e-24 g/cm^3 has a conduction time short enough to stop the run dead.
  // Above rmax the radiative conduction flux is zero, the cell is dropped from the
  // conduction timestep, and the optical depth integrated down from the top accumulates
  // at rad_kappa_above instead of the table value -- so tau does not grow through the
  // corona and the blend weight w stays 0 there.  The two-stream reads the same two
  // numbers (see two_stream_rt.hpp) so both radiative operators go quiet together.
  Real rad_kappa_rmax = 0.0;
  // rad_kappa_above [cm^2/g]: the opacity used above rad_kappa_rmax.  0 (the default)
  // makes the corona perfectly transparent; the grey chain handles kappa = 0 exactly
  // (dtau = 0 gives e0 = 0, so both the source and the emission vanish identically).
  Real rad_kappa_above = 0.0;
  // rad_gate_rho [g/cm^3] > 0: the DENSITY form of the same idea, and the one to use
  // whenever the artificial medium can exchange gas with the star.  The radius test
  // above is a fixed shell: it makes the corona inert, but it ALSO strips the opacity
  // from stellar gas that inflates or is ejected past that radius -- which is exactly
  // what the open-top red-giant runs do from t ~ 5.7e5 on (mass above 3.6e12 grows by
  // 4.5-14x, the photosphere crosses the cutoff, and the two-stream then radiates 9000 K
  // gas as a bare blackbody at 2-4 L).  Gating on the LOCAL DENSITY instead makes the
  // opacity follow the gas: the 1e-15 corona stays inert wherever it is, and stellar
  // material stays opaque wherever it goes.  The gate is a logistic in log10 rho (see
  // RadGate) rising from 0.1 to 0.9 over rad_gate_dex decades, and the effective opacity
  // is  kappa_eff = G kappa_table + (1 - G) rad_kappa_above.  Mutually exclusive with
  // rad_kappa_rmax (fatal if both are set): one criterion or the other.
  Real rad_gate_rho = 0.0;
  Real rad_gate_dex = 0.5;
  // rad_implicit_x1 (<hydro>/ or <mhd>/rad_implicit_x1, default false): solve the RADIAL
  // radiative diffusion IMPLICITLY (backward Euler) instead of adding it to the
  // face fluxes.  The explicit radial operator has dt ~ dx^2 rho c_v/kappa_rad with
  // kappa_rad ~ T^3/(kappa_R rho): in an evacuated cell just above the photosphere
  // that limit falls below 1e-4 s while the run takes 30 s steps, and because dt is
  // evaluated from the PREVIOUS state the overshoot is unbounded -- the cell runs away
  // to 1e10 K in a single step.  With this on, the interior x1 faces carry no explicit
  // flux; a per-column tridiagonal solve (Conduction::ImplicitRadialUpdate) applies the
  // same flux-limited operator unconditionally stably, and the x1 conduction timestep
  // constraint is dropped.  The whole radial extent must live in ONE MeshBlock (the
  // solve is column-local, no MPI), which the constructor checks.  x2/x3 stay explicit.
  // In MHD the frozen internal energy also has the magnetic energy taken out of it
  // (MagEnergyCC); the operator itself is identical.
  bool rad_implicit_x1 = false;
  // scratch for the tridiagonal solve, allocated once: slots (e*, T*, 1/(rho c_v), p,
  // A_f K_f/dl_f, c', d') per cell/face of every column.  See imp_* below.
  static constexpr int nimpw = 7;
  static constexpr int IMPE = 0, IMPT = 1, IMPA = 2, IMPP = 3;
  static constexpr int IMPC = 4, IMPCP = 5, IMPDP = 6;
  DvceArray5D<Real> imp_wrk;
  DvceArray1D<int> imp_flag;   // 1 element: has the first bad cell been recorded?
  DvceArray1D<Real> imp_rec;   // 8 elements: that cell's identity and state
  // ---- <hydro>/<mhd>/rad_x1_uform and rad_x1_kiter: the two cures for the frozen-K
  // error linearising the radial conduction in T with the face conductance frozen makes.
  // Both default to off/1, and with them off not one
  // expression of the solve changes (bitwise).
  //
  // rad_x1_uform.  At frozen opacity the radiative flux is EXACTLY linear in u = T^4,
  //   F = -(a c/(3 kappa_R rho)) du/dl = -(K/(4 T^3)) du/dl,
  // so solving for du instead of dT removes the T^3 part of the nonlinearity outright
  // (measured: it halves the frozen-K flux error).  The face recipe is untouched -- the
  // off-diagonal is the SAME C_f divided by 4 T_f^3 with T_f the frozen face
  // temperature, so the two forms coincide where T is uniform -- the diagonal becomes
  // (rho c_v V)_i/(4 T_i^3 beta_dt), and the energy increment comes back as
  // de_i = (rho c_v)_i du_i/(4 T_i^3), through the same conservative positivity clip.
  // Still symmetric, still an M-matrix.
  //
  // rad_x1_kiter.  The residue the u-form leaves is the frozen d ln kappa_R/d ln T
  // (~0 in a He-star convection zone, -2.4 at the B-star Fe bump).  kiter > 1 re-runs
  // the whole assembly -- T, c_v, kappa_R(T,p), the limiter and the face T^3 -- at the
  // state the previous pass wrote, and re-solves: a Picard iteration on the backward-
  // Euler balance.  Passes 2..k keep the heat-capacity term anchored on T^n (imp_tn),
  // exactly as the merged two-stream passes do, so each pass CORRECTS the step already
  // taken instead of taking another full one.  Both switches are refused together with
  // the merged two-stream column solve, which carries its own linearisation.
  bool rad_x1_uform = false;
  int rad_x1_kiter = 1;
  DvceArray4D<Real> imp_tn;    // T^n, the state pass 0 started from (kiter > 1 only)
  int imp_lines = 0;           // lines printed so far by the debug report
  void ImplicitRadialUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                            const Real beta_dt, const bool rt_on = false,
                            const int rt_pass = 0);
  // ---- THE MERGED TWO-STREAM COLUMN SOLVE (<problem>/rt_implicit_column) -------------
  // The grey two-stream source is exactly linear in the cell Planck functions B_j at
  // frozen opacity, so its nearest-neighbour linearisation folds straight into the
  // tridiagonal above: the column is then solved for the radiative exchange AND the
  // radiative diffusion at once, by one backward-Euler Newton step, instead of the
  // two_stream's per-cell relaxation which damps each cell by its own (1 - e^-x)/x and
  // so breaks the O(1e3)-to-O(1) cancellation between neighbours (the dt-linear velocity
  // pump; see two_stream_rt.hpp, rt_implicit_column).
  //
  // two_stream_rt fills these three arrays instead of applying its own de, and then
  // calls ImplicitRadialUpdate itself, once per outer pass.  The wrapper task
  // ImplicitConduction is a no-op while rt_col_active: the solve has already happened.
  //   rt_col_res  : R_i, the FULL explicit two-stream source at the current state, in
  //                 code energy density per code time (what the old apply block would
  //                 have multiplied by bdt).
  //   rt_col_jac  : dR_i/dB_{i-1}, dR_i/dB_i, dR_i/dB_{i+1} at frozen opacity, summed
  //                 over the quadrature chains.  Everything >= 2 cells away stays in R.
  //   rt_col_dbdt : dB_i/dT_i with T in KELVIN, i.e. 4 B_i/T_i.  The Kelvin -> code
  //                 conversion is applied here, where temperature_cgs() lives.
  //   rt_col_tn   : T* at the START of the outer iteration, needed by passes 2..k so the
  //                 heat-capacity term stays anchored on T^n and the conduction operator
  //                 is not re-applied in full on every pass.
  //   rt_col_dtex : for a THIN cell (see rt_impl_tau_min), the CHANGE in its own Planck
  //                 function that its per-cell relaxation has just applied.  A thick row
  //                 next to it keeps its Jacobian entry dR_i/dB_thin, but as a KNOWN
  //                 right-hand-side contribution rather than an unknown, so the exchange
  //                 across a thick/thin interface is still counted exactly once.
  bool rt_col_active = false;       // <problem>/rt_implicit_column, set by two_stream_rt
  Real rt_col_dtmax = 0.25;         // <problem>/rt_impl_dtmax: cap on |dT|/T per pass
  bool rt_col_verbose = false;      // one-shot assembly dump of one column at cycle 0
  DvceArray4D<Real> rt_col_dtex;
  bool rt_col_alloc = false;
  DvceArray4D<Real> rt_col_res;
  DvceArray5D<Real> rt_col_jac;
  DvceArray4D<Real> rt_col_dbdt;
  DvceArray4D<Real> rt_col_tn;
  // the M-matrix audit of the merged rows: how many rows had |off-diagonals| exceeding
  // the diagonal, and the worst excess ratio.  Reported once per rad_col_report cycles.
  // 6: nviol, worst ratio, sum V de, expected sum V de, dT caps, worst capped |dT|/T
  DvceArray1D<Real> rt_col_diag;
  int rt_col_lines = 0;
  void EnableRTColumn();
  // rad_angular (<hydro>/ or <mhd>/rad_angular, default true): the MASTER SWITCH of the
  // HORIZONTAL (x2/x3, "angular"/"transverse") part of the radiative conduction.  True is
  // everything below and is bitwise the historical arithmetic.  FALSE removes the
  // horizontal operator outright:
  //   * AddIsotropicHeatFluxRadiative returns after the x1 faces, so NO radiative heat
  //     flux is added to any x2/x3 face -- on every mesh type, cubed-sphere seam and
  //     spherical-polar paths included, since those live inside the same two kernels;
  //   * BuildAngularCoeffs returns immediately, so no frozen transverse conductance,
  //     no cap_x/cap_c2/cap_c3 and no rad_cap_ang report are built (this also covers the
  //     problem generators that call it themselves under problem/rt_split_transverse);
  //   * NewTimeStep drops the transverse conduction limits dt2 and dt3.
  // Everything RADIAL is untouched: the x1 face flux, rad_implicit_x1, the tau blend
  // (rad_tau_lo/hi, rad_blend_radial), rad_flux_inner and the two-stream hand-over
  // (rad_blend_use_2s, rad_f2s).  Ordinary (constant / spitzer) thermal conduction is not
  // affected at all -- the switch is read and acted on only for
  // isotropic_conduction = radiative.
  //
  // WHAT IT IS FOR.  The horizontal radiative diffusion is the stiff, expensive and
  // physically dubious half of the operator (see rad_tr_tau_lo below); a run that wants
  // no horizontal radiative transport at all -- because a column solver owns that
  // transport, or because the arm is a deliberate ablation -- had to fake it by setting
  // the smooth taper past the top of the domain (rad_tr_tau_lo = 1e30).  That route is
  // bitwise equivalent to this one but still pays for the angular kernels, and it needs
  // the tau blend to be on.  This switch is the direct statement.
  //
  // COMBINATIONS.  rad_angular = false is a FATAL with rad_implicit_ang, rad_sts_all,
  // rad_sts_split and rad_ang_solver = adi: those ARE the transverse treatment, and
  // silently running them with no transverse physics would be a trap.  rad_cap_ang > 0
  // is instead accepted SILENTLY and becomes a no-op (the production inputs carry
  // rad_cap_ang = 0.5 and must not have to be edited): nothing angular is built for it.
  // Note that the dt-collapse report in Mesh::NewTimeStep keys its "capped" label on
  // rad_cap_ang alone, so with rad_angular = false and rad_cap_ang = 0 it prints the
  // sentinel dt2 = dt3 = -1 instead of a word; -1 means "no transverse constraint".
  bool rad_angular = true;
  // rad_cap_ang (<hydro>/ or <mhd>/rad_cap_ang, default 0 = off): a CONSERVATIVE cap on
  // the explicit ANGULAR (x2/x3) radiative diffusion.  With the radial direction made
  // implicit the same evacuated-cell runaway simply migrates to the angular faces (the
  // I1 test died at t = 5.76e5 with T = 3e10 K in a cell 60x below its shell median),
  // because the angular operator is still explicit and its limit dt ~ dx^2 rho c_v/kappa
  // is evaluated from the PREVIOUS state.  Define per cell the explicit angular
  // stiffness
  //     x_i = beta_dt alpha_i (sum over its 4 angular faces of A_f K_f/dl_f)/V_i
  //         = dt/dt_cond,angular,
  // with alpha_i = 1/(rho_i c_v,i) and K_f the SAME frozen face conductivity the x2/x3
  // kernels use (RadFaceKappa x flux limiter x tau-blend weight, divided by sin(alpha)
  // on the cubed sphere, which is the coefficient of (T_j - T_i)/dl in the face flux).
  // Every angular face flux is then multiplied by
  //     f_f = min(1, rad_cap_ang/max(x_i, x_j)),
  // one number per face applied to both of its cells: exactly conservative, purely
  // local, and it makes the effective x of every cell <= rad_cap_ang, so at 0.5 the
  // operator is unconditionally stable and monotone (the checkerboard mode's
  // amplification factor stays >= 0).  Requires rad_implicit_x1; drops dt2/dt3.
  Real rad_cap_ang = 0.0;
  // rad_implicit_ang (<hydro>/ or <mhd>/rad_implicit_ang, default false): solve the
  // TRANSVERSE (x2/x3) radiative diffusion with an unconditionally stable operator-split
  // update instead of adding it to the face fluxes.  rad_cap_ang slows the transverse
  // exchange down to whatever an explicit step can carry, which is stable and
  // conservative but WRONG: in the iron-convection-zone boxes it throttles the very
  // horizontal radiative exchange that sets the structure.  With this on, the x2/x3
  // faces carry no explicit flux, the x2/x3 conduction timestep constraints are dropped,
  // and Conduction::ImplicitTransverseUpdate (conduction_transverse.cpp) advances the
  // per-layer 2D operator over the stage with an RKL1 super-time-stepping loop.
  // v1 is CARTESIAN and uniform-grid only, and mutually exclusive with rad_cap_ang.
  // rad_implicit_x1 is NOT required: the two operators are split from each other and
  // from the hydro, and the radial direction keeps whatever treatment it was given.
  bool rad_implicit_ang = false;
  // rad_tr_split_out (set by a problem generator, not by the input file): the TRANSVERSE
  // operator has been taken OUT of the RK stage and is run by the same operator-split
  // step that runs the radiation column (box_convection's problem/rt_strang and
  // problem/rt_imex, through problem/rt_split_transverse).  Hydro::/MHD::
  // ImplicitTransverseConduction is then a no-op: solving here as well would apply the
  // operator twice a stage.
  bool rad_tr_split_out = false;
  int rad_ang_maxit = 200;      // ceiling on the RKL1 substage count of one call
  bool rad_ang_verbose = false; // report the substage count and the conservation residual
  static constexpr int ntrs = 2;
  static constexpr int TRST = 0, TRSA = 1;   // frozen T*, and alpha = 1/(rho c_v)
  DvceArray5D<Real> tr_st;
  // the RKL1 registers Y_{j-2}, Y_{j-1}, Y_j on the energy increment, one variable each
  // so the dedicated MeshBoundaryValuesCC below can exchange whichever one is current
  DvceArray5D<Real> tr_ya, tr_yb, tr_yc, tr_ycoar;
  MeshBoundaryValuesCC *pbval_tr = nullptr;
  int ang_lines = 0;
  // ---- rad_ang_solver (<hydro>/ or <mhd>/rad_ang_solver, "sts" by default): WHICH
  // solver advances the transverse system.  "sts" is the RKL1 super-time-stepping loop
  // and is the default, so the switch is bitwise inert unless it is set.  "adi" is an
  // ALTERNATING-DIRECTION IMPLICIT step: the SAME frozen coefficients, the SAME flux-form
  // 5-point stencil and the SAME conservation check, but the linear ODE
  //     dy/dt = b + A y,   A = A2 + A3 (the x2 and the x3 faces),   y(0) = 0
  // is advanced over the whole stage tau = beta_dt in ONE step by a sequential
  // backward-Euler (LOD) splitting -- two tridiagonal LINE solves instead of s explicit
  // substages -- which is stiffly accurate, unlike the theta-weighted (Douglas) scheme
  // this module used to also offer.
  //
  // COST.  Two line solves and NO cell halo exchange at all: the only communication is
  // the T*/alpha refresh exchange the operator already does, the Gershgorin Allreduce,
  // and one tiny ring gather of interface coefficients per sweep.  RKL1 pays one stencil
  // + one halo per substage, and s grows as the square root of the stiffness.
  //
  // THE LINE SOLVE CROSSES MeshBlocks AND RANKS -- a line along x2 runs through all
  // nblk2 = <mesh>/nx2 / <meshblock>/nx2 blocks of its row, periodically.  It is solved
  // by PARTITION (the SPIKE / reduced-system pattern of the mode-3 column partition in
  // utils/two_stream_column_partition.hpp): each block runs Thomas on the INTERIOR of its
  // own piece of the line for three right-hand sides -- the data, and the two unit
  // vectors on its first and last row -- which gives its two interface unknowns as an
  // affine function of its two neighbours' facing interface unknowns.  Those 6 numbers
  // per line are gathered around the ring of blocks (nblk-1 rounds of a neighbour shift,
  // a device copy when the neighbour is on this rank and an MPI message when it is not),
  // the resulting 2*nblk reduced system is solved REDUNDANTLY and in one canonical
  // (absolute-block-index) order on every block of the ring -- so all of them get
  // bitwise-identical interface values -- and each block back-substitutes locally.
  //
  // WHAT IS REFUSED: rad_sts_split (an RKL1 stability construction with no meaning for
  // a direct solve) and rad_sts_all (ADI is the transverse operator only).  A direction
  // decomposed into more than one MeshBlock must be PERIODIC (an open chain of blocks
  // would need a second gather pass that is not implemented).  rad_tr_halo_every is
  // IGNORED with a warning: there is no substage recurrence whose ghost skin it could run
  // down.  rad_sts_once is orthogonal and works unchanged.
  bool rad_ang_adi = false;
  // rad_adi_scheme (lod | lod2 | lod2a, "lod" by default): WHICH splitting the two
  // sweeps implement.  See the long note at the schemes in conduction_transverse.cpp.
  //
  // lod is only FIRST order, and under-damps a moderately stiff mode (1/(1+z) against
  // the exact e^-z): on the Gaussian test its L1 is 5.6x RKL1's at 32^2 and on the B-star
  // gate box dT_tau10/dT_tau1 come out 7-8x RKL1's.
  //   lod2  -- RICHARDSON extrapolation in the step size: y = 2 LOD(tau/2)^2 - LOD(tau),
  //            the two half-steps taken with the sweep order alternated.  SECOND order,
  //            and still L-stable: R(z) = 2/(1+z/2)^2 - 1/(1+z) -> 0 as z -> infinity and
  //            |R| <= 1 on z >= 0.  THREE sweep pairs per stage and TWO factorisation
  //            sets (tau and tau/2 share nothing), i.e. ~3x the lod cost.
  static constexpr int ADISCM_LOD = 1;
  static constexpr int ADISCM_LOD2 = 2;
  static constexpr int ADISCM_LOD2A = 4;   // lod2 with the half-step order alternated
  int rad_adi_scm = ADISCM_LOD;
  // rad_adi_cross_iter (cubed sphere only, default 1): outer iterations of the ADI step
  // in which the EXPLICIT metric cross term is re-evaluated at the current answer
  // Th = T* + alpha y instead of at T*.  The fixed-point contraction factor is the
  // cross/diagonal symbol ratio, <= 1/4 in the stiff mid-band and -> 0 at the Nyquist
  // checkerboard, so 2 iterations remove the mid-band ringing described in the design
  // note.  1 = the plain lagged cross term.  Ignored off the cubed sphere.
  int rad_adi_cross_iter = 1;
  // rad_adi_seam_w (cubed sphere + rad_ang_solver = adi, default 0.5): where the
  // pair-implicit PANEL-SEAM sub-step evaluates the temperature, between the state the
  // directional sweeps started from (0) and the one they ended at (1).  0.5 is the
  // trapezoidal rule and is the default; see conduction_transverse.cpp.
  Real rad_adi_seam_w = 0.5;
  static constexpr int NADIB = 8;      // ceiling on blocks per line (reduced size 2*N)
  DvceArray5D<Real> tr_aw, tr_acp;     // the second spike, and the Thomas scratch
  DvceArray5D<Real> tr_yf;             // lod2 only: the full-step answer for Richardson
  DvceArray5D<Real> tr_ylg;            // rad_adi_cross_iter > 1: the lagged answer the
                                       // explicit metric cross term is evaluated at
  DvceArray4D<Real> tr_ared;           // (m, block, line, 6): the interface coefficients
  DvceArray1D<int> tr_aact, tr_ab0;    // per plane: active?  per block: its ring index
  HostArray1D<int> tr_aact_h, tr_ab0_h;
  int adi_nb2 = 1, adi_nb3 = 1;        // blocks along a line, per direction
  int adi_nlmax = 1;
#if MPI_PARALLEL_ENABLED
  MPI_Comm adi_comm;                   // its own communicator for the ring gather
  bool adi_comm_set = false;
#endif
  void ImplicitTransverseUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                const Real beta_dt);
  // rad_sts_all (<hydro>/ or <mhd>/rad_sts_all, default false): ONE super-time-stepped
  // radiative conduction operator over ALL THREE directions, replacing the split of an
  // implicit tridiagonal x1 solve (rad_implicit_x1) plus the RKL1 transverse operator
  // (rad_implicit_ang).  The stencil of ImplicitTransverseUpdate is extended to 7 points
  // by the x1 faces, every direction is linearised about the SAME refreshed state, and
  // the whole thing is advanced by one RKL1 loop -- so there is no splitting error
  // between the radial and the horizontal operator at all.  That splitting error is not
  // academic: linearising the transverse operator about the state the radial solve had
  // already moved is what grew max|de| by ~1.5x per cycle in the He-star FeCZ box (see
  // the file comment of conduction_transverse.cpp).  Implies rad_implicit_ang (it IS
  // the transverse treatment), is mutually exclusive with rad_implicit_x1 and
  // rad_cap_ang, and inherits every restriction of rad_implicit_ang: Cartesian,
  // uniform grid, no SMR/AMR, nghost >= 2, plus the whole x1 extent in one MeshBlock.
  //
  // WHAT IT IS FOR, AND WHAT IT IS NOT FOR.  An RKL1 super-step costs s substages for
  // s^2+s explicit steps, i.e. the SQUARE ROOT of the stiffness ratio -- but ONE s is
  // chosen for the whole mesh from the global maximum of the Gershgorin radius, so the
  // cost is set by the single stiffest cell in every direction at once.  In an
  // ISOTROPIC CARTESIAN BOX (the FeCZ boxes, the solar box) all three directions have
  // the same dx and the same K, the x1 faces raise the row radius by a factor of order
  // 3/2, and the substage count grows by only ~20 %: the unified operator is then
  // strictly better than the split, since it costs almost nothing extra and removes a
  // first-order-in-dt error term.
  //
  // On a STRETCHED SPHERICAL (or cubed-sphere) grid it is the WRONG CHOICE, and the
  // curvilinear guards below refuse it outright.  There the radial direction is both the
  // finest (dr/r ~ 1e-3 near the photosphere against a degree in angle) and the one
  // carrying the whole stellar flux, so its stiffness runs orders of magnitude above the
  // transverse one; folding it into the same RKL1 loop would set s from the radial rows
  // and multiply the cost of the transverse operator by that ratio's square root, where
  // the tridiagonal solve of ImplicitRadialUpdate handles exactly that direction in ONE
  // unconditionally stable sweep because it is column-local.  Radially stiff, angularly
  // mild -> keep rad_implicit_x1 + rad_implicit_ang.  Isotropic -> rad_sts_all.
  bool rad_sts_all = false;
  void StsConductionUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                           const Real beta_dt);
  // rad_sts_split (<hydro>/ or <mhd>/rad_sts_split, default false): the STIFFNESS SPLIT
  // of the super-time-stepped operator.  Every face conductance is written as
  //     C_f = C_exp,f + C_sts,f,   C_exp,f = min(C_f, C_max,f),   C_sts,f = C_f - C_exp,f
  // where C_max,f is the largest conductance the CURRENT hydro step can carry
  // EXPLICITLY on that face.  C_exp is added to the ordinary face fluxes in the flux
  // kernels of AddIsotropicHeatFluxRadiative, i.e. inside the RK stages, and C_sts alone
  // goes into the RKL1 loop.  Two things follow, and they are the whole point:
  //   * the non-stiff part of the operator -- which is most of it away from the stiff
  //     cells (37 % of cells in the B-star FeCZ box, nearly all in the He box) -- is
  //     advanced BY THE RK INTEGRATOR, so it is coupled to the hydro at the integrator's
  //     own order instead of being operator-split from it at first order, and
  //   * the RKL1 loop sees a smaller Gershgorin radius, so it takes fewer substages, and
  //     a MeshBlock in which C_sts vanishes on every face is skipped outright
  //     (sts_blk below).
  // Both parts are in FLUX form, so each is separately conservative.  The dt limiter is
  // untouched: C_exp is bounded BY CONSTRUCTION, never by dt.
  //
  // C_max,f.  The explicit stability limit this code already uses in NewTimeStep is, per
  // cell, dt <= cfl/(2 ndim) dx^2 rho c_v/K per direction, which is exactly the row-sum
  // statement
  //     x_i = beta_dt alpha_i (sum over the cell's faces of C_f)/V_i <= cfl,
  //     alpha_i = 1/(rho_i c_v,i)
  // (see the derivation in BuildAngularCoeffs).  The explicit part is therefore given the
  // budget x_i <= rad_sts_split_x (0.5 by default, the value at which the operator is
  // monotone), distributed EVENLY over the nf = 2 x ndim faces of the row, and a face is
  // held to the smaller of the budgets of the two cells that share it:
  //     C_max,f = rad_sts_split_x min(V_i/alpha_i, V_j/alpha_j)/(nf beta_dt).
  // Taking the min makes the bound hold for BOTH rows, and the even split makes the row
  // sum at most nf x (budget/nf) = budget.  Symmetric in i and j, so the two cells (and
  // the two MeshBlocks at a block face) form bitwise the same number and the explicit
  // part stays exactly conservative.
  //
  // A face the RKL1 loop treats as CLOSED (a physical, non-periodic x2/x3 boundary, and
  // the two physical x1 faces) gets no explicit part either: the split only redistributes
  // the operator that is already there, it does not open faces.
  //
  // NOTE for rad_implicit_ang WITHOUT rad_sts_all: the x1 faces are then not part of this
  // operator at all and keep their own explicit treatment, whose share of the row sum the
  // conduction dt limiter still bounds by cfl/ndim.  The total explicit row sum is then
  // at most rad_sts_split_x + cfl/ndim rather than rad_sts_split_x.
  bool rad_sts_split = false;
  Real rad_sts_split_x = 0.5;
  // the explicit FRACTION C_exp,f/C_f of each face, 0 when the face carries nothing;
  // allocated only when rad_sts_split is set.  cap_c1/2/3 then hold C_sts, not C_f.
  DvceArray4D<Real> cap_f1, cap_f2, cap_f3;
  // per MeshBlock: does ANY face of it still carry C_sts > 0?  A block where the split
  // left nothing behind is skipped by the RKL1 stencil (its increment is zero there).
  DvceArray1D<int> sts_blk;
  bool sts_blk_used = false;
  // rad_sts_once (default false): apply the RKL1 operator ONCE per cycle, after the LAST
  // RK stage and with the FULL dt, instead of once per stage with beta_dt.  The substage
  // count grows as sqrt(tau), so one call over dt costs sqrt(2) x one call over dt/2 and
  // the two-stage integrator saves ~30 % of the operator.  The splitting between this
  // operator and the hydro is FIRST ORDER either way -- running it every stage does not
  // make it second order, because the operator is not part of the RK right-hand side --
  // so this trades nothing but the size of the first-order term.
  bool rad_sts_once = false;
  // The round-off margin on the RKL1 stability limit is a fixed 0.10 (see
  // conduction_transverse.cpp): an s-substage RKL1 super-step covers |lambda| tau <=
  // s^2+s, and the substage count is chosen from R = 0.5 (1 + 0.10) max_i z_i.
  // rad_sts_perplane (default false): PER-PLANE SUBSTAGE COUNTS.  The transverse
  // operator couples cells only inside a horizontal plane of fixed x1 index i -- its
  // stencil has no x1 face at all -- so the planes are INDEPENDENT linear systems that
  // happen to be integrated in one kernel.  One global substage count therefore makes
  // every plane pay the stiffness of the stiffest one (145-197 substages in the
  // whole-column B-star box, where the deep planes need ~15).
  //
  // With this on, the Gershgorin radius is reduced PER PLANE (a global max over m, k, j
  // at fixed i, reduced across MPI), each plane gets its own substage count s_i and its
  // own RKL1 w1 = 2/(s_i^2+s_i), and the loop runs s_max = max_i s_i substages during
  // which a plane with j > s_i is skipped by the stencil -- its register is already its
  // final value and is simply carried forward.  mu_j and nu_j do not depend on s, so the
  // only per-plane coefficient is w1.  Each plane is then exactly the s_i-stage RKL1
  // scheme it would have got on its own, so the answer is a converged RKL1 solution of
  // the same operator plane by plane, and conservation is untouched: the flux form is
  // per plane, and the two cells sharing an x2/x3 face are always in the SAME plane and
  // multiply the same mu~_j.
  //
  // BITWISE identical to the global-s loop when every plane has the same s_i (a state
  // uniform in x1, e.g. inputs/tests/rad_transverse_gauss.athinput).
  //
  // NOT available with rad_sts_all: the 7-point stencil has x1 faces, which couple the
  // planes, and the substage count then has to be global.  It is silently ignored there.
  // Needs the whole x1 extent in one MeshBlock, so that a local index i IS a plane.
  bool rad_sts_perplane = false;
  // rad_tr_window (default true): restrict the three transverse halo exchanges to the
  // x1 plane bracket that the substage actually reads (see conduction_transverse.cpp).
  // false = exchange the FULL x1 range every time, which is the pre-window behaviour
  // and a control for any placement- or window-dependent difference.
  bool rad_tr_window = true;
  // ---- COMMUNICATION COST OF THE RKL1 TRANSVERSE OPERATOR.  This switch changes only
  // WHICH GHOST CELLS ARE EXCHANGED AND WHEN; the arithmetic every ACTIVE cell performs
  // is untouched, so it is bitwise-inert at its default and meant to stay bitwise when
  // on.  See the bookkeeping note in conduction_transverse.cpp.
  //
  // rad_tr_halo_every (int, default 1 = exchange every substage): exchange the RKL1
  // increment only every N substages and compute the intermediate substages on a GHOST
  // SKIN as well, letting the skin shrink by one layer per substage.  N is clamped to
  // nghost.  N >= 3 also exchanges the SECOND register at each refill (see the note in
  // conduction_transverse.cpp: the two-iterate recurrence otherwise runs the skin out).
  int rad_tr_halo_every = 1;
  DvceArray1D<Real> tr_zpl;     // the Gershgorin radius of each x1 plane
  DvceArray1D<int> tr_spl;      // ... its substage count s_i ...
  DvceArray1D<Real> tr_w1pl;    // ... and its RKL1 w1 = 2/(s_i^2 + s_i)
  HostArray1D<Real> tr_zpl_h;
  HostArray1D<int> tr_spl_h;
  HostArray1D<Real> tr_w1pl_h;
  // cumulative substage count of this operator, reported alongside the per-call one so
  // that the COST of a run can be read off a single line
  std::int64_t sts_nsub_tot = 0;
  std::int64_t sts_ncall = 0;
  void BuildAngularCoeffs(const DvceArray5D<Real> &w0, const EOS_Data &eos,
                          const Real beta_dt);
  // beta_dt of the CURRENT stage.  The angular fluxes are added in Hydro::Fluxes, which
  // runs before the RK update, so the cap has no other way of knowing the step it is
  // capping; Hydro::Fluxes sets this immediately before calling AddHeatFluxes.
  Real stage_beta_dt = 0.0;
  DvceArray4D<Real> cap_x;    // the stiffness x_i, active cells + one angular ghost
  DvceArray4D<Real> cap_c1;   // ... and on the x1 faces (rad_sts_all only)
  DvceArray4D<Real> cap_c2;   // A_f K_f/dl_f on the x2 faces
  DvceArray4D<Real> cap_c3;   // ... and on the x3 faces
  // CUBED SPHERE (and rad_cs_exact): the METRIC CROSS-TERM coefficient of the same
  // faces, G_f = C_f dl_f cos(alpha).  The gnomonic xi/eta lines meet at an angle, so
  // the total flux through an x2 face is NOT -C_f (T_j - T_i) but
  //     Phi_f = -C_f (T_j - T_i) + G_f ge_f,
  // with ge_f the face-tangential (eta) derivative, the same two-cell average the
  // explicit face flux of AddIsotropicHeatFluxRadiative uses.  The cross term couples
  // j to k and is therefore NOT part of the tridiagonal/5-point structure the solvers
  // invert: it is carried EXPLICITLY (see docs/dev/cs_implicit_transverse.md).  Empty
  // on a Cartesian mesh, and every read of it is guarded by the flag.
  DvceArray4D<Real> cap_g2;
  DvceArray4D<Real> cap_g3;
  DvceArray1D<int> cap_cnt;   // 2: cells with x_i > cap, faces actually capped
  DvceArray1D<Real> cap_rec;  // 6: max x_i and the cell that carries it
  int cap_lines = 0;
  Real cap_diag_x = 0.0;      // the largest x_i of the last call (rank-local)
  int cap_diag_over = 0;      // cells with x_i > rad_cap_ang in the last call
  // rad_cs_exact (default true): the exact face-normal derivative on the cubed sphere;
  // false drops the metric cross term and the 1/sin(alpha) -- DIAGNOSTIC only
  bool rad_cs_exact = true;
  int rad_kr_nT = 0, rad_kr_nP = 0;
  DvceArray2D<Real> rad_kr_tab;            // (iT, iP) log10 kappa_R [cm^2/g]
  DvceArray1D<Real> rad_kr_lT, rad_kr_lP;  // log10 T [K], log10 p [dyn/cm^2], ascending
  DvceArray4D<Real> rad_w, rad_tauf;
  // rad_w is filled by BuildRadWeights, which runs as a task inside the stage.  Until
  // it has, every weight is zero and the timestep below reads that as "no face carries
  // any diffusive flux" -- so the FIRST cycle would run at the hydro timestep with the
  // conduction operator fully on.  NewTimeStep builds the weights itself if this is
  // still false, which is the case at initialisation.
  bool rad_w_built = false;
  // rad_blend_use_2s: REMOVED as a runtime switch (measured a net loss on every arm it
  // was tried on -- see git history).  This module can no longer turn it on, but the
  // two members below survive because src/utils/two_stream_rt.hpp still names them
  // (guarded by rad_blend_use_2s > 0, which is now permanently false): rad_f2s stays
  // unallocated and rad_f2s_ready stays false forever.
  int rad_blend_use_2s = 0;
  DvceArray4D<Real> rad_f2s;
  bool rad_f2s_ready = false;
  void BuildRadWeights(const DvceArray5D<Real> &w, const EOS_Data &eos);
  // PER-CYCLE DIAGNOSTIC (deep_hot_jupiter_rt's problem/diag_gid).  Off by default, and
  // cond_diag is not allocated until EnableDiag is called, so it costs nothing when off.
  // (m,slot,k,j,i) with slot 0/1/2 = the energy flux this module ADDS to
  // flx1/flx2/flx3(m,IEN,k,j,i) on the x1/x2/x3 face of index i/j/k, slot 3 = the
  // radiative conductivity kappa_rad at the cell exactly as NewTimeStep forms it (code
  // units), slot 4 = the flux-limited effective diffusivity keff on the x1 direction,
  // slot 5 = the cell's own conduction dt candidate, the minimum over the three
  // directions BEFORE the dimensional factor fac and before the CFL number.
  bool diag = false;
  DvceArray5D<Real> cond_diag;
  void EnableDiag(int nmb, int n3, int n2, int n1);
  Real kappa_iso;            // isotropic thermal conductivity
  Real kappa_iso_limit;      // limit to isotropic thermal conductivity

  // functions
  void AddHeatFluxes(const DvceArray5D<Real> &w, const EOS_Data &eos,
                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxRadiative(const DvceArray5D<Real> &w, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx);
  void AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                       DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);

 private:
  // the ONE RKL1 loop both super-time-stepped entry points run: `with_x1` = false is the
  // 5-point transverse operator (rad_implicit_ang), true the 7-point unified one
  // (rad_sts_all).  See conduction_transverse.cpp.
  void RklConductionUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                           const Real beta_dt, const bool with_x1);
  MeshBlockPack* pmy_pack;
  // "hydro" or "mhd": identifies which physics module owns this Conduction object, so
  // that the general EOS derived-variable arrays (cached temperature and pressure) of the
  // right module can be found.
  std::string my_block;
};
//! \fn Real RosselandTable
//! \brief log-bilinear lookup of a tabulated kappa_R(T,p), indices clamped to the grid
KOKKOS_INLINE_FUNCTION
Real RosselandTable(const DvceArray2D<Real> &tab, const DvceArray1D<Real> &lT,
                    const DvceArray1D<Real> &lP, const int nT, const int nP,
                    const Real tk, const Real pcgs) {
  const Real x = log10(tk), y = log10(pcgs);
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (!(x > lT(0))) {
    i = 0; fx = 0.0;
  } else if (x >= lT(nT-1)) {
    i = nT-2; fx = 1.0;
  } else {
    while (i < nT-2 && lT(i+1) <= x) ++i;
    fx = (x - lT(i))/(lT(i+1) - lT(i));
  }
  if (!(y > lP(0))) {
    j = 0; fy = 0.0;
  } else if (y >= lP(nP-1)) {
    j = nP-2; fy = 1.0;
  } else {
    while (j < nP-2 && lP(j+1) <= y) ++j;
    fy = (y - lP(j))/(lP(j+1) - lP(j));
  }
  const Real lk = (1.0-fx)*((1.0-fy)*tab(i,j) + fy*tab(i,j+1))
                +      fx *((1.0-fy)*tab(i+1,j) + fy*tab(i+1,j+1));
  return pow(10.0, lk);
}

//! \fn Real KappaTemp
//! \brief the temperature used in the radiative conductivity: the cell temperature,
//! capped at tmax when tmax > 0 (see Conduction::rad_tmax).  Bitwise identity when off.
KOKKOS_INLINE_FUNCTION
Real KappaTemp(const Real tk, const Real tmax) {
  return (tmax > 0.0) ? fmin(tk, tmax) : tk;
}

//! \fn Real RadGate
//! \brief the density gate G(rho): 1 where the gas is dense enough to carry its
//! tabulated opacity, 0 in the artificial low-density medium, a logistic in log10 rho
//! between.  Exactly
//!     G = 1/(1 + exp(-(log10(rho) - log10(rho_gate))/s)),   s = dex/(2 ln 9),
//! so G = 1/2 at the threshold and G runs from 0.1 to 0.9 over exactly `dex` decades
//! (G = 0.9 at u/s = ln 9, and the width in u is 2 s ln 9 = dex).  rho <= 0 gives 0; the
//! exponent is clamped so a 20-decade excursion cannot overflow.
//! ONE definition, used at every opacity site in the conduction module and the
//! two-stream, so the two operators cannot disagree about which gas is radiative.
KOKKOS_INLINE_FUNCTION
Real RadGate(const Real rho_cgs, const Real rho_gate, const Real dex) {
  if (!(rho_gate > 0.0)) return 1.0;      // gate off
  if (!(rho_cgs > 0.0)) return 0.0;
  const Real s = dex/(2.0*log(9.0));
  Real u = (log10(rho_cgs) - log10(rho_gate))/s;
  if (u > 40.0) return 1.0;
  if (u < -40.0) return 0.0;
  return 1.0/(1.0 + exp(-u));
}

//! \fn Real RadBlendWeight
//! \brief the diffusion weight of a face: 0 for tau <= lo, 1 for tau >= hi, a raised
//! cosine in log tau between
KOKKOS_INLINE_FUNCTION
Real RadBlendWeight(const Real tau, const Real lo, const Real hi) {
  if (tau <= lo) return 0.0;
  if (tau >= hi) return 1.0;
  const Real s = log(tau/lo)/log(hi/lo);
  return 0.5*(1.0 - cos(M_PI*s));
}

//! \fn Real RadTaperWeight
//! \brief the TRANSVERSE taper weight of a face (Conduction::rad_tr_tau_lo/hi): the same
//! raised cosine in log tau as RadBlendWeight, and identically 1 when the taper is off,
//! so every call site is bitwise inert until rad_tr_tau_lo is set.
KOKKOS_INLINE_FUNCTION
Real RadTaperWeight(const Real tau, const Real lo, const Real hi) {
  return (lo > 0.0) ? RadBlendWeight(tau, lo, hi) : 1.0;
}

#endif // DIFFUSION_CONDUCTION_HPP_
