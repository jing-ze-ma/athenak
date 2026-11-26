#ifndef RAD_VET_RAD_VET_HPP_
#define RAD_VET_RAD_VET_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet.hpp
//! \brief the grey two-moment (M1 / VET) radiation module <rad_vet>: the RadiationVet
//! class, its constants (array slots, key enums) and small inline helpers.

#include <array>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <utility>
#include <vector>
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"
#include "utils/time_table.hpp"
#include "rad_vet/rad_vet_defs.hpp"
#include "rad_vet/rays/sc_band.hpp"

// forward declarations
class Driver;

namespace radvet {
//----------------------------------------------------------------------------------------
//! \class RadiationVet
//! \brief the grey M1 / VET radiation module

class RadiationVet {
 public:
  RadiationVet(MeshBlockPack *ppack, ParameterInput *pin);
  ~RadiationVet();
  // the constructor's parts
  void ReadSchemeKeys(ParameterInput *pin);    // rad_vet_params.cpp
  void CheckGeometry(ParameterInput *pin);     // rad_vet_geometry.cpp
  void ReadMatterKeys(ParameterInput *pin);    // rad_vet_params_matter.cpp
  void ReadClosureKeys(ParameterInput *pin);   // rad_vet_params_closure.cpp

  // <time>/dt_retry: dtr_on: a failed cycle is retaken with a smaller dt; a solve that
  // would stop the run raises dtr_fail with reason dtr_why instead
  bool dtr_on = false;
  bool dtr_fail = false;
  // 1 NON-CONVERGED, 2 DIVERGED, 3 DIVERGED (masked), 4 eint<=0
  int dtr_why = 0;
  Real dtr_resid = 0.0;
  template <class V> void DtRetryVisit(V &v);   // rad_vet_dtretry_visit.hpp (generated)
  void DtRetryQuiesce(bool restoring);          // no host thread writes while restoring
  void DtRetryAfterRestore();                   // make the async prefetch state safe

  // parameters read from <rad_vet>
  Real c_light{};      // speed of light in code units (required)
  Real chat{};         // = c_light (no reduced speed of light)
  Real e_floor{};      // floor on E

  // opacities per unit mass (code units): kappa0, or kappa0 (rho/rho_ref)^a (T/t_ref)^b
  int opacity_type{};           // VET_OPAC_*
  Real kappa_p{}, kappa_e{}, kappa_f{}, kappa_s{};
  Real opac_rho_ref{}, opac_t_ref{}, opac_a{}, opac_b{};
  bool opac_zero{};             // all four opacities are identically zero
  VetOpacTab otab;               // opacity = table: filled by SetOpacityTables

  // <rad_vet>/force_reference: with VET_FREF_WB_ARAD the gas gets the residual between
  // the radiative force and the reference acceleration arad_ref of an external
  // well-balanced scheme (filled by the pgen through SetForceReference)
  int force_ref{};
  DvceArray4D<Real> arad_ref;   // (m,k,j,i), a reference x1 acceleration, code units
  // optional explicit radiation energy source (code units), set by a pgen with
  // esrc_on = true (implicit only): adds dt*(chat/c)*esrc to E's old vector and slope
  DvceArray4D<Real> esrc;
  bool esrc_on = false;
  // vet_source_noesrc: the ray source reads E^n without the
  // explicit esrc deposit (be only)
  bool vsrc_noes = false;
  bool VsrcNoes() const;

  // the fluid is <hydro> or <mhd> (FluidRef). Under MHD the gas internal energy subtracts
  // emag0 = |B|^2/2 of the current field (b0) or emag1 of the U^n field (b1), EmagBuild;
  // hydro: unread dummies
  bool fl_on = false;           // a fluid exists (<hydro> or <mhd>)
  bool fl_mhd = false;          // it is <mhd>
  bool emag_hold = false;       // Opacity must not rebuild emag0 (the U^n swap)
  DvceArray4D<Real> emag0, emag1;
  void EmagBuild(bool un);      // un = false: emag0 from b0; true: emag1 from b1
  // after a host write of the fluid's u0: restrict, exchange, BCs, prolongate, ConToPrim.
  // MHD exchanges u only through pbval_u (MHD::InitRecv would wait on B forever).
  void FluidRefresh(Driver *pdrive, bool lat = true);
  // floorcap_rho: after the implicit solve, cells with rho < floorcap_rho get
  // |v| <= floorcap_v and T_gas <= T_rad; the removed energy goes to E_rad (FloorGasCap)
  Real fcap_rho = 0.0, fcap_v = 0.0;
  int fcap_report = 0;
  bool fcap_lz = false;     // floorcap_lz: the velocity cap keeps the L_z momentum
  bool fcap_spin = false;   // floorcap_lz_spin: its L_z X the dense gas's dOmega = X/I
  DualArray1D<Real> fcap_cnt;   // removed: 0 n_v, 1 n_t, 2 dKE, 3 dEint, 4 |dp|V, 5 dL_zV
  void FloorGasCap();


  // matter coupling
  bool coupling{};     // run the implicit local solve at all
  bool gas_feedback{}; // write the solve's reaction back to the gas (off: static medium)
  Real arad{};         // radiation constant in CODE units: the equilibrium energy
                       // density is arad*T^4 with T the EOS's code temperature

  // debug: dbg_gas_force = false gives the gas no momentum and no work
  bool dbg_gas_force{};

  // ---- Cartesian short characteristics (rays/sc_cart_*.cpp): Eddington tensor K/J of a
  // grey
  // formal solution, once per hydro step from the start-of-step state
  bool vet_sc{};             // closure = sc on a Cartesian mesh
  int vet_nmu = 0, vet_nphi = 0, vet_nray = 0;  // mu nodes per hemisphere, azimuths, rays
  bool vsc_rayk = false;     // vet_sc_kernel = ray: ray-parallel layer launch
  DvceArray4D<Real> vet_rayg;   // its per-ray face sums (block, 5 nray, k, j)
  bool vet_bc_bath = false;  // a Marshak x1 end gives the rays its bath
  // x1 boundary intensity of vet_sc and vet_gd (VetBcInit). vet_bc_inner: entering
  // intensity at a Marshak / flux inner face, 0 diff, 1 emit, 2 auto (emit at Marshak,
  // diff at flux). vet_bc_reflect = mirror: reflecting ends mirror the leaving rays.
  int vet_bc_inner = 0;
  bool vet_bc_mirror = false;   // mirror at a reflecting x1min (vet_bc_reflect = mirror)
  bool vet_bc_mirror_hi = false;   // ... at a reflecting x1max (mirror or auto, vet_sc)
  void VetBcInit(ParameterInput *pin, const bool gd);
  bool VetBcEmit() const;    // the inner rule at the current implicit_bc_x1min is emit
  Real VetBathLo() const;    // vet_sc: entering inner intensity of emit (-1: diff)
  bool vet_bc_face = false;  // vet_bc_diff = face: inner diffusion I at the face
  bool vet_x1per = false;    // vet_x1_periodic: periodic x1 sweep
  int vet_x1npass = 0;       // its number of passes
  // vet_x1_ptol: stop the periodic sweep when d_p/(1 - q) < vet_x1_ptol (0: npass)
  Real vet_x1ptol = 0.0;
  DvceArray4D<Real> vet_x1ent;   // (nmb, nray, k, j): the previous pass's exit plane
  // a build has run, so the next adaptive sweep starts warm from vet_x1ent; both are in
  // the restart file (kVetX1eRstMagic; a file without them restarts the sweep cold)
  bool vet_x1_warm = false;
  int X1eRstNslab(const int n1) const {
    return vet_x1_warm ? (vet_x1ent.extent_int(1) + n1 - 1)/n1 : 0;
  }
  void X1eRstPack(DvceArray5D<Real> &dst, const int nmb) const;
  bool X1eRstUnpack(const DvceArray5D<Real> &src, const int nmb, const int nsl);
  Real vet_x1_nbuild = 0.0, vet_x1_npsum = 0.0, vet_x1_ncap = 0.0;
  int vet_x1_npmax = 0;
  int vet_mirror_npass = 4;  // passes of a cavity mirrored at both x1 ends
  bool vet_x1src = false;    // vet_x1_periodic_inflow = source
  // a non-periodic x2 (x3) mesh boundary is a vacuum face (I = 0) for the formal solution
  bool vet_x2open = false, vet_x3open = false;
  // vet_sc angle set and its rotation on the vet_gd clock (keys vet_gd_rotate_*).
  // Double-Gauss: azimuths about x1 turn by (2 pi/nphi) q; ring: in-plane angles by
  // 0.8 (q - 1/2)(2 pi/N); q in [0, 1) the block's low-discrepancy fraction.
  int vet_ring_n = 0, vet_ring_g = 1;   // vet_ring, vet_ring_groups (0: double-Gauss)
  bool vsc_rot = false;                 // vet_sc rotates (physical clock)
  Real vsc_beta = 0.0;                  // current rotation of the vet_sc set
  bool vsc_gd = false;       // vet_sc_set = gd: the vet_gd direction set, pole on x1
  int vsc_gdmir = 0;         // ray r + nray/2 is the x1 mirror of ray r (1) or not (-1)
  void ScCartAngles(const Real beta);
  void ScGdAngles(const Real beta);
  void ScBuildSet(std::vector<double> &d) const;
  void VetRotInit(ParameterInput *pin);   // (rays/sc_directions.cpp)
  void VetRotUpdate();
  Real VetRotBeta(const int64_t blk) const;   // rotation of block blk
  void VetRotSet(const Real beta);
  void VetScRotAvgSetup();     // vet_sc's averaged fields (vet_cell)
  void ScCartSweep();    // build with the current set
  bool vet_full = false;     // vet_tensor = full: the solve reads all six D_ab = K_ab/J
  static constexpr Real vet_eig_min = 0.0;   // full: eigenvalue floor of D
  Real vet_nguard{}, vet_ncell{}; // full: cell-calls the guard changed D / all cell-calls
  Real vet_guard_max{};      // full: largest |D_guarded - D_raw| (max norm)
  Real vet_d11_min{};        // full: smallest guarded D_11
  Real vet_time{}, vet_itime{}, vet_ncall{}; // SC seconds, implicit-solve seconds, calls
  DvceArray2D<Real> vet_ang;   // (nray, 4): mu1, mu2, mu3, weight (sum of weights = 1)
  DvceArray5D<Real> vet_cell;  // (nmb, VETRAY_NC, k, j, i), VETRAY_*
  DvceArray5D<Real> vet_ipl;   // (nmb, 2, nray, k, j): intensity of the last two layers
  // multi-block sweep state (plane buffers, neighbour tables, exchange); null on one
  // MeshBlock
  ScCartMBState *vet_mbs;

  // fold dtnew into Mesh::NewTimeStep (no sub-cycling, or no other module sets dt)
  bool sets_mesh_dt{};

  // implicit transport (rad_vet_implicit.hpp)
  int nstage{};                 // 1: one stage solve per task-list pass
  Real impl_cfl{};              // <rad_vet>/implicit_cfl: dtnew = impl_cfl*dx/c (<=0 off)
  DvceArray1D<Real> opn_nskip_d;
  // an energy-conserving positivity limiter of the written-back gas eint;
  // implicit_pos_floor: the E floor takes its energy from the gas. pos_cnt_d: device
  // counters (VET_POS_*), pos_cnt the all-rank totals
  bool impl_pos_floor = false;
  DvceArray1D<Real> pos_cnt_d;
  Real pos_cnt[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // E-floor energy created by discarded attempts (rank-local), not in pos_cnt
  Real pos_un_drop = 0.0;
  int dtr_ncut = 0, dtr_maxcut = 8;   // DtRetry cuts of this cycle, limit
  bool dtr_redo = false;   // this cycle's failed step is redone (CLEAN RETRY)
  Real marshak_q{};             // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min{}, ibc_x1max{}; // VET_IBC_*
  // implicit_bc_x2min|x2max|x3min|x3max: auto (periodic on a periodic mesh, else
  // reflect) | reflect | marshak (F_f = c marshak_q E outward; formal solution I = 0)
  int ibc_x2min = 2, ibc_x2max = 2, ibc_x3min = 2, ibc_x3max = 2;
  Real iflux_x1min{}, iflux_x1max{}; // imposed face flux of VET_IBC_FLUX
  Real iebath_x1min{}, iebath_x1max{}; // VET_IBC_MARSHAK: incident bath
  // per-column inflow of the outer x1 Marshak face (SetX1maxBathColumns, (nmb, n3, n2)):
  // E_in > 0 is a thick inflow face, comoving flux 0 plus the enthalpy inflow
  // v1 (1 + chi) E_in; E_in = 0 keeps the Marshak face. ibadv_x1max: lab-frame advective
  // energy flux through the outer face.
  bool iebath_hi_col_on = false;
  DvceArray3D<Real> iebath_hi_col;
  // implicit_bc_x2min|max = flux: imposed lab-frame F_theta(m,k,i) (+ toward larger
  // theta) from SetX2FaceFlux, explicit in the solve
  bool ifx2_on = false;
  DvceArray3D<Real> ifx2_lo, ifx2_hi;
  DvceArray3D<Real> ibadv_x1max;
  // implicit_beam_n > 0 (vet_sc): up to 2 collimated beams entering x1min on an x2
  // patch, one SC ray each. ibm_i(m,q,k,j): entering intensity of beam ray q; ibm_f:
  // its face flux c sum w mu1 I, which replaces the Marshak row where > 0
  int ibm_n = 0, ibm_nr = 0;
  int ibm_ray[4] = {-1, -1, -1, -1};
  DvceArray4D<Real> ibm_i;
  DvceArray3D<Real> ibm_f;
  DvceArray4D<Real> ibm_s;   // beams' stress column c sum w Omega_b mu1 I (clean)
  void ImplBeamInit(ParameterInput *pin);
  DvceArray4D<Real> f0x1;       // face-normal comoving x1 flux (m,k,j,i), restart state
  DvceArray5D<Real> iw;         // per-cell work array of the solve, VET_NIW components
  Real impl_nstep{}, impl_itsum{}, impl_itmax{};  // Picard counters
  // implicit_thin_freeze: cells with c dt rho kappa_P below it keep their start-of-step
  // opacities for the whole solve
  DvceArray4D<Real> thin_frz;   // 1 = frozen cell of this solve
  // implicit_resid_fatal (0 off): stop when a solve ends with resid above it or
  // non-finite, before the write-back
  Real impl_res_fatal = 0.0;
  // implicit_timers = N: fenced timers of the stage tasks and ImplicitSolve parts over
  // cycles >= N, printed by ImplicitReport (diagnostic)
  int tmr_c0 = 0;
  bool tmr_on = false;
  double tmr_acc[12] = {}, tmr_cnt[8] = {}, tmr_last = 0.0;
  Kokkos::Timer tmr_t;
  void TmrMark(int c);
  bool dbg_hydro_off = false;   // diagnostic: the Driver skips hydro stages
  // force_reference_work (wb_arad): full = the solve gives the gas the work of the full
  // force; split = the WB source gives rho arad_ref's work at its stage, the solve only
  // the residual's; the radiation pays both
  bool fref_wsplit = false;
  bool fref_wsplit_ok = false;  // set by a pgen whose WB source gives the reference work
                                // (before SetForceReference, which resolves auto)
  bool fref_wsplit_auto = false;  // force_reference_work = auto, not yet resolved
  // split with be: the pgen adds its increment to fref_wacc (acc = gam0_s acc + incr,
  // as the gas energy combines) and the be write-back makes E pay exactly fref_wacc
  // (times dt_solve/dt_mesh)
  DvceArray4D<Real> fref_wacc;
  int fref_hnst = 0;
  Real fref_hgam0[4] = {0.0, 0.0, 0.0, 0.0}, fref_hbeta[4] = {0.0, 0.0, 0.0, 0.0};
  bool FrefWaccOn() const;
  // clean staggered hesdirk2 (time2_fref_wacc): each stage solve pays the reference work
  // of the hydro stage before it (gam0 = 0); no slope carries a payment
  bool FrefWaccStage() const;
  // gam0 of the hydro stage whose beta dt is bdt (FATAL if none matches)
  Real FrefWaccGam0(Real bdt, Real dt_mesh) const;
  // the accumulator, allocated to the shape of arad_ref
  DvceArray4D<Real> FrefWacc();

  // ---- implicit face fluxes (constants: rad_vet_implicit.hpp)
  // half-range / thin-switch face options (read only when named):
  //  implicit_hr_wscale_x1 / _lat: scale the half-range face coefficients (diagnostic)
  //  implicit_hr_fcimp: dissipation form F = F_c(E') + w (c/2)(s_L E_L - s_R E_R),
  //    s = h+ - h-, AL = 0, HCL >= 0, HCR <= 0 (M-matrix)
  //  implicit_hr_recon_sfgate: plm on h only off the beams (Jh > J: donor)
  //  implicit_hr_recon_lim = smooth0: slope 0 where the smoothness test failed
  //  implicit_hr_recon_hjump = X: donor where |dh+| + |dh-| > X between neighbours;
  //    _hjump_w = W: gate reach; implicit_hr_recon_kjump = K: donor where chi differs by
  //    > K within W
  //  implicit_hr_recon_hsmooth = C: donor where |dh+ - dh-|/(|dh+| + |dh-| + eps) > C
  //    within W (eps = _hsmooth_eps)
  //  implicit_hr_recon_hlim: sig = min(sig_E, sig_H+, sig_H-)
  //  implicit_hr_recon_qe_tau = T: on cells with chi dx_d > T the slope is h_c dE sig/4
  //  implicit_thin_gface = max | mean | jump (mean + (g_R - g_L)^2/_gface_d): face
  //    thinness g_f of the switch weight
  //  implicit_thin_gface_shell: x1-face g from shell means (repro::ShellSums)
  //  implicit_hr_fsc: in the blend band the central branch takes c (h_L E_L + h_R E_R)/2
  //  implicit_thin_gface_filter = L: x1-face g from a hat filter of half-width L
  //  implicit_thin_gw_dump: g_L, g_R and weight AL of the last build (output vet_gw)
  DvceArray5D<Real> gwd;
  // time2_face_fct: the restart file carries K1's face energy flux t2pf1 (one slab per
  // direction, kVetFaxRstMagic). Old files' implicit_hr_cmem_rst slabs: read, ignored.
  int FaxRstNslab() const {
    const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
    return (t2_pfst && t2pf1_ok && t2pf1.extent_int(0) > 0) ? nd : 0;
  }
  void FaxRstPack(DvceArray5D<Real> out, int nmb);
  bool FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl);
  // transverse (x2, x3; cs seams) faces, fmm_t(m,n,k,j,i) at the face: x2 [0] Llo,
  // [1] Lhi (E'-linear coefficients, chat units), [2] G = cr F(E^k) - Llo E_lo -
  // Lhi E_hi, [3] lam; x3 [4..7]; [8..11] F_L coefficients of a limited face (x2, x3),
  // F_L = chat ((1 + s) E_lo - (1 - s) E_hi)/2; [12], [13] lagged off-diagonal part of
  // the x2 / x3 flux (lab units). iw iw_fmm+3..+6: dt A/V of x2+, x2-, x3+, x3- faces.
  DvceArray5D<Real> fmm_t;
  // Zalesak rounds exchange R only when a cell within 2 ng of a face
  // changed. Counters: rounds, exchanges done.
  // A round evaluates only cells whose stencil changed in the
  // previous round or its R exchange (others cannot be newly limited: bitwise)
  Real zl_rounds = 0.0, zl_halos = 0.0;
  // implicit_thin_switch (VET_THSW_*): dt-free thin/thick weight frozen per solve
  // (formal-solution |H|/J, C1, E^n Knudsen number); implicit_thin_corr adds the lagged
  // formal-solution flux
  int impl_thsw = 0;
  bool thsw_c1 = false;           // the sweep accumulates VETRAY_NL
  Real thsw_h0 = 0.3, thsw_h1 = 0.6;
  // implicit_thin_h0_x1 / _h1_x1: thin-switch band of the x1 faces (_hx1_curv: only on
  // curvilinear meshes); implicit_thin_recon = js (VetThinJs: ratios from the face mean
  // of H+-); implicit_thin_bq = lin (outer Marshak q from the face-extrapolated H+);
  // implicit_thin_hcsf = jh: hc/hcs take the sf donor scale Jh (vet_sc), j: the cell J
  bool impl_beam_fs = false;    // implicit_flux_beam = fs: |f| = |H|/J of vet_sc
  // implicit_flux_beam = halfrange: F = (1 - w) F_central + w c (r+_L E_L + r-_R E_R),
  // r+- = H+-/J the half-range ratios along the normal, w = exp(-(tau_f/tau0)^2)
  // (implicit_blend = tau) or 1 (berthon)
  bool impl_beam_hr = false;
  // implicit_hr_face = sf (default with vet_sc): face ratios from the face mean of the
  // half-range fluxes over the donor scale Jh = max(J, outgoing SC face fluxes) where J
  // does not carry the SC flux (VetHrSf); cell: the cell ratios
  bool impl_hr_sf = false;
  // implicit_hr_excess = iso (vet_gd): the part of E^n above the formal-solution J takes
  // the isotropic half-range shape (+-1/4) (VetHrExcess)
  double hsf_ncell = 0.0, hsf_nsf = 0.0;   // cell-solves, of which Jh > J (rank 0)
  DvceArray5D<Real> vgd_hr;   // (m, 6, k, j, i): r+_r r-_r r+_t r-_t r+_p r-_p
  // with implicit_thin_switch = fs* or implicit_thin_corr vgd_hr has VET_VGD_HRNX
  // components: + C1, eps, J
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_pin_cnt;
  // ---- implicit transverse (x2/x3) coupling
  bool trans_on{};              // transport == implicit and the mesh is multi-D
  bool trans_x3{};              // ...and three-dimensional
  DvceArray4D<Real> f0x2, f0x3;   // face-normal comoving x2 / x3 fluxes (m,k,j,i),
                                  // restart state (trans_on only)

  // ---- the linear solve
  int impl_prec_fb{};           // implicit_precond_fallback: -1 none, else impl_prec code
  bool bcg_capped{};            // the last BiCGStab stopped at lin_maxit, test unmet
  // ---- Krylov GPU cost options (off: nothing allocated)
  DvceArray5D<Real> odc;        // (m,3,k,j,i) operator cache
  DvceArray5D<Real> ost;        // (m,19,k,j,i) stencil
  bool st_edges{};              // any edge coefficient is non-zero
  int impl_onep{};              // implicit_one_pass: check period N (0 = off)
  Real onep_qa[3]{}, onep_qb[3]{};  // last two measured contractions (be, stage 1, 2)
  Real onep_cnt[3]{};           // solves since the last measurement
  bool impl_onep_auto{};        // implicit_one_pass_auto: switch one_pass off per kind
                                // of solve when it does not pay
  Real onep_off[3]{};           // 1: one_pass switched off for this kind (be, s1, s2)
  Real onep_actr[3]{};          // on: eligible solves since the last one-pass acceptance;
                                // off: switched-off solves since it was switched off
  Real onep_prb[3]{};           // 1: on as a re-probe, not yet confirmed by an acceptance
  int impl_pord = 2;            // implicit_predictor_order of an old restart file
  int impl_prec{};              // implicit_precond: 0 x1 line Jacobi, 1 symmetric
                                // red-black transverse line Gauss-Seidel (block-local),
                                // 2 its forward half, 3 mg (rbgs_fwd + block-local
                                // semicoarsening), 4 mg_gc (mg + global coarse space)
  Real bcg_r0rel{};             // max|b - A x0|/max|b| of the last BiCGStab call
  DvceArray1D<Real> lay_part;   // layout_invariant_mode = block: per-MeshBlock
                                // partials of VetLayReduce (sized on first use)
  Kokkos::View<Real*, Kokkos::HostSpace> lay_part_h;   // their host copy
  Real ew_fprev{}, ew_etaprev{}; // max|r0| and eta of the previous pass of this step
  bool pred_ok{};               // a previous-step increment is stored
  Real pred_dt{};               // dt of that step
  DvceArray5D<Real> ipred;      // (m,3,k,j,i): dE_prev, dT_prev, T^n of this step

  // ---- time_scheme = hesdirk2: two implicit stage solves inside the Heun stages
  int time_scheme{};            // VET_TIME_BE | VET_TIME_HESDIRK2
  bool t2_ok{};                 // a valid FSAL slope K1 is stored
  int t2_solve{};               // what the next ImplicitSolve does (VET_T2S_*)
  bool t2_fail{};               // the last stage solve was not admissible
  int t2_dbg_fail{};            // debug time2_dbg_fail: fail stage 1 at this cycle
  void Time2Stash();
  // Ebnd (diagnostic, ebnd_on set by a pgen): radiation energy put in through the
  // physical x1 faces (inner minus outer) over the accepted steps. A hesdirk2 step
  // weights the face flux L = sum A F1 of K1, K2, K3 with the tableau weights of
  // U^{n+1}; a backward-Euler step uses L of its solve.
  bool ebnd_on = false;
  Real ebnd_acc = 0.0;            // accumulated energy (code units)
  // ebnd_split: inner (in) and outer (leaving) faces accumulated apart, ebnd_acc their
  // difference
  bool ebnd_split = false;
  Real ebnd_acc_in = 0.0, ebnd_acc_out = 0.0;
  // with ebnd_split also the esrc deposit (chat/c) sum esrc V with the same tableau
  // weights (K1 carries the previous step's last solve, FSAL)
  Real ebnd_acc_src = 0.0;
  Real ebnd_s1 = 0.0, ebnd_s2 = 0.0;
  // esrc sum of the last accepted solve (in the next step's K1); a restart takes it from
  // the drive trailer
  Real ebnd_s3 = 0.0;
  bool ebnd_s3_ok = false;
  Real EbndEsrc();
  void EbndStart();
  void EbndStage1();
  void EbndAccept(Real dt, bool t2step);
  // implicit_hold: cells a pgen holds (hold_m = 1) are Dirichlet rows of every implicit
  // solve (E' = held E, no increment, kick or work); the write-back leaves them to the
  // pgen
  bool impl_hold = false;
  DvceArray4D<Real> hold_m;
  Real t2_nstep{}, t2_nbe{}, t2_nfall{}; // stage steps, BE steps, fallbacks (counters)
  // time2_auto_be = N (0 off): after a stage fallback the next steps take the
  // backward-Euler redo path at once; stages are retried after 1, 2, 4, ... N skipped
  // steps (reset on success). The skip-run state is in the restart file.
  int t2_abe = 16;              // N
  int t2_abe_len = 0;           // current skip-run length (0 = stage mode)
  int t2_abe_left = 0;          // skipped steps left in the current run
  bool t2_abe_now = false;      // this cycle skips the stages
  // time_scheme = be runs as hesdirk2 with every step a time2_auto_be skip step
  // (Time2AutoBeBegin)
  bool t2_beonly = false;
  // the next step's hydro CFL dt is taken again after the radiation
  // update, on the gas it wrote (a backward-Euler kick can exceed the pre-coupling CFL)
  Real t2_abe_n = 0.0;          // steps taken by time2_auto_be (counter)
  Real t2_abe_ntry = 0.0;       // stage retries after a skip run
  Real t2_abe_nok = 0.0;        // of which succeeded
  // hesdirk2 wall-time split (rank 0, diagnostic): clean stage steps, failed attempts,
  // their redo, auto_be steps, other BE steps, discarded dt_retry attempts
  Real t2w_clean = 0.0, t2w_fail = 0.0, t2w_redo = 0.0;
  Real t2w_abe = 0.0, t2w_be = 0.0, t2w_retry = 0.0;
  bool Time2AutoBeBegin();
  void Time2AutoBeEnd(int kind, Real tw_a, Real tw_b);
  Real t2_dtprev{};             // dt of the previous step (vet_sc extrapolation)
  bool t2_vprev{};              // vet_prev holds the tensor of the previous step
  // time2_tableau: 0 hesdirk2 (Heun + H-ESDIRK2), 1 trbdf2 (TR-BDF2, RK2 hydro with
  // c2 = 2 - sqrt 2). Time2FormStage: stage 1 t2inc = t2_a21 dt K1; stage 2 start =
  // t2_w0 Y2 + t2_w1 U^n, t2inc = dt (t2_c31 K1 + t2_c32 K2); diagonal g = 1 - 1/sqrt 2.
  int t2_tab = 0;
  Real t2_a21 = 0.0, t2_w0 = 0.5, t2_w1 = 0.5, t2_c31 = 0.0, t2_c32 = 0.0;
  Real t2_lin_save{};           // implicit_lin_tol put back after a stage solve
  Real t2_nclip{};              // vet_sc cells whose extrapolated tensor was clipped
  DvceArray5D<Real> t2k1;       // (m,VET_T2_NK,k,j,i) the FSAL slope K1 (restart state)
  DvceArray5D<Real> t2k2;       // (m,VET_T2_NK,k,j,i) the stage-2 slope K2
  DvceArray5D<Real> t2inc;      // (m,VET_T2_NK,k,j,i) old vector - stage start state
  DvceArray4D<Real> t2f1, t2f2, t2f3;  // face fluxes of U^n (Heun average)
  bool pred2_ok{};
  Real pred2_dt{};
  DvceArray5D<Real> vet_prev;   // vet_sc: start-of-step tensor of the previous step
  void Time2Init(ParameterInput *pin);
  bool Time2Active();           // the next step runs the stages (else BE)
  void Time2FormStage(int stage, Real dt);
  // sdirk2 with Heun hydro: the stage-1 solve on the gas at t^n + g dt, the hydro's
  // stage-2 flux at t^{n+1} (t2hs: the gas U* of Heun stage 1)
  void Time2GasTime(Driver *pdrive, int which);
  DvceArray5D<Real> t2hs;
  int t2_gt = 0;   // last Time2GasTime step of this stage loop
  bool t2_gts = false;  // time2_gtime_slope: its eint by the limited exchange slope
  bool t2_gtp = true;   // time2_gtime_pos: the positive flux state, u0 restored
  void Time2Restore(Driver *pdrive);

  bool t2_inprobe = false;      // inside a helper solve (Time2RotCorrect)
  bool t2_x2old = false;        // a restart's time2_vet (time2_vet_col) = extrap2
  static constexpr Real kPredTcap = 1.1;   // first iterate: T within this factor of T^n
  // implicit_pos_floor: the E floor creates no energy.
  // The raise is paid by the cell's gas above its floor, then by the radiation above
  // e_floor of the same radial column, nearest cells first (donor slopes and
  // |F| <= c E kept). The column is in one MeshBlock: local, exact, layout invariant.
  // VET_POS_FLR_MV counts moves, VET_POS_FLR_UN columns with nothing left.
  // the gas raise the cell's E cannot cover also comes from the radial column's radiation
  DvceArray4D<Real> pgc_def;      // (m,k,j,i) uncovered gas raise x volume
  DvceArray4D<Real> pfc_def;      // (m,k,j,i) uncovered floor energy of the cell
  // old files (version < 4) carry time2_vet_col = extrap2 channels: read past
  static constexpr int VET_X2_NCH = 24;
  void Time2Report();
  int Time2RstNchWant();        // restart channels this run keeps (0 unless hesdirk2)
  int Time2RstNch();            // ...and writes (0 unless a slope is stored)
  bool Time2RstAccept(const int nfile, const int ver);   // a file's channels fit (sets
                                                         // t2_nfv, t2_rstd3)
  int t2_nfv = 0;               // conservative-ray intensity channels of the file
  int t2_rcl = -1;              // closure channels of the file read (closure_form =
                                // iterate, RstNcl; -1: this run's own layout)
  int t2_rskip = 0;             // a file's ring channels read with the ring off
  bool t2_rstd3 = false;        // version-1 file: the fkin block is d only (no K1_d)
  int t2_rstv = 5;              // clean file version: <= 3 zero predictor / extrap2
                                // channels and the ring in CLN_N ghosted slabs per
                                // record, 4 the cubed sphere's sc rays' cln_fvi ghosted
  void Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0 = 0, int w0 = 0);
  void Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb);

  // ---- off-diagonal Eddington terms and the closure lag (multi-D implicit)
  // <hydro>/rad_signal_speed reads opac(VET_OP_T) and the closure tensor of the last VET
  // step: they travel in the restart (kVetRssRstMagic), else are filled from the initial
  // state before the first dt (RssPrime).
  //! closure tensor the signal speed reads: 0 E/3, 4 vet_cell (D_dd), 5 cln_s (D, clean);
  //! 1 - 3 were removed closures
  int RssMode() const;
  int RssRstMode() const;   // RssMode + 16 (extended restart state)
  int RssRstNch(int mode) const;  // channels written: opac(VET_OP_T) + the tensor's
  //! array and channel of restart channel n of a mode
  void RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c);
  void RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode, int m0 = 0);
  void RssPrime(bool restart);
  DvceArray5D<Real> rss_stage;   // (nmb, nch, k, j, i) read from the restart file
  int rss_stage_mode = -1;       // mode of rss_stage, -1 = nothing staged

  // implicit_vimp_ewf: the lagged face enthalpy work E_f (1 + D_f) dv_f is taken as
  // (1 + D_f)(E_L dv_L + E_R dv_R)/2, each kick with its own cell's E
  bool vimp_now{};              // on for this step (the positivity fallback drops it)
  // cubed sphere: 3-component exchange (scalar, face-normal pair) of iw slots c0..c0+2
  DvceArray5D<Real> d3w;
  MeshBoundaryValuesCC *pbval_d3 = nullptr;
  // the implicit_hr_recon limiter ghosts (hr quad map across a seam)
  MeshBoundaryValuesCC *pbval_sg = nullptr;
  // time2_face_fct: Zalesak positivity limiter on the face fluxes of the explicit stage
  // increment t2inc. Each stage solve stores t2pf(m,d,k,j,i) = (chat/c) F_applied with
  // its slope; the face part is G_f = A_f (c1 t2pf1 + c2 t2pf2), the rest cell-local.
  // A cell with E + dE <= 0 gets R = (1 - 1e-12) Q/P (Q = V E + V dE + P - I + I(R),
  // P/I its out/inflow); a face takes its donor's R and scales the E flux and face F
  // slope on both sides (energy and momentum exact). Rounds until no cell is newly
  // limited, R halo-exchanged; the rest (Q <= 0, gas eint) goes to posfix.
  DvceArray5D<Real> t2pf1, t2pf2;   // (m,3,k,j,i) (chat/c) F_applied of K1 / K2's solve
  bool t2pf1_ok = false, t2pf2_ok = false;
  // t2pf is stored with the limiter or with clean staggered (BoundaryEnergyAcc): cln_lacc
  // (m,4,k,j) = energy through the physical x1 faces so far, inner (0, + inward), outer
  // (1, + outward), and the old vector's face part (2, 3); cln_ldt = the step's dt
  // (0: none yet, the history reads the faces' F)
  bool t2_pfst = false;
  DvceArray4D<Real> cln_lacc;
  Real cln_ldt = 0.0;
  void BoundaryEnergyAcc(int stage, Real c1, Real c2, Real w, Real dt);
  Real t2_ffct_calls = 0.0, t2_ffct_cells = 0.0, t2_ffct_left = 0.0;
  Real t2_ffct_rmax = 0.0, t2_ffct_de = 0.0, t2_ffct_nop = 0.0, t2_ffct_fill = 0.0;
  // t2_ffct_emin = th: the limiter keeps th V E^n in a limited cell; counted: cells
  // left below th E^n, min E/E^n of limited cells
  Real t2_ffct_emin = 0.0, t2_ffct_nemp = 0.0, t2_ffct_qmin = 1.0e300;
  void Time2FaceFct(int stage, Real c1, Real c2);
  bool am_vec = false;             // <hydro>/am_vector: the vector form
  DvceArray4D<Real> amk1, amk2;   // spin slopes (FSAL / stage 1)

  // ---- Anderson acceleration of the Picard map, implicit_accel (ImplicitAccelApply)
  DvceArray5D<Real> aa_xc;      // (m,nc,k,j,i) x_k, the scaled state entering the pass
  DvceArray5D<Real> aa_fc;      // (m,nc,k,j,i) g_k = G(x_k) - x_k, scaled
  DvceArray5D<Real> aa_xp, aa_fp;  // previous pass' pair, for the differences
  int aa_nh{};                  // history columns in use
  int aa_head{};                // ring head (column the next difference overwrites)
  bool aa_hasp{};               // aa_xp/aa_fp hold a usable previous pass
  Real aa_fnp{};                // ||g|| of the previous pass (< 0: none)

  // evolved variables
  DvceArray5D<Real> u0;         // (E, F1, F2, F3)
  DvceArray5D<Real> u1;         // state at the start of the substep
  DvceArray5D<Real> coarse_u0;  // on 2x coarser grid (for SMR/AMR)
  DvceFaceFld5D<Real> uflx;     // fluxes on cell faces
  DvceArray5D<Real> opac;       // rho*kappa per cell, VET_NOPAC components

  // implicit-solve counters (VET_CNT_*), Reals so that long-run sums cannot wrap
  DualArray1D<Real> cnt;

  // boundary communication
  MeshBoundaryValuesCC *pbval_u;

  // timestep bookkeeping
  Real dtnew{};   // implicit_cfl*min(dx)/c (implicit_cfl > 0), else FLT_MAX
  Real dt_sub{};  // dt of the stage solve being taken
  int nsub{};     // always 1 (no sub-cycling)

  RadiationVetTaskIDs id;

  // functions
  void AssembleRadVetTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  //! dt_sub = dt_mesh (nsub = 1); returns nsub
  int SetSubsteps(Real dt_mesh);
  //! print the implicit-solve counters (called from the destructor, rank 0)
  void ReportCounters();
  //! the two opacity tables (opacity = table), set by the pgen
  void SetOpacityTables(const DvceArray2D<Real> &kr, const DvceArray2D<Real> &kp,
                        const DvceArray1D<Real> &lT, const DvceArray1D<Real> &lD,
                        const int nT, const int nD);
  //! reference acceleration of force_reference = wb_arad (caller keeps it alive)
  void SetForceReference(const DvceArray4D<Real> &a);
  void SetX1maxBathColumns(const DvceArray3D<Real> &eb);
  void SetX2FaceFlux(const DvceArray3D<Real> &flo, const DvceArray3D<Real> &fhi);
  //! print v_max*tau_max/chat once, from the first Opacity task with a filled array
  // ---- implicit transport (rad_vet_params_implicit.cpp)
  //! read the implicit-solver parameters, allocate its arrays
  void ImplicitInit(ParameterInput *pin);
  // IC radiation start, fresh starts only (rad_vet_steady.cpp)
  int icr_mode = 0;             // ic_rad: 0 none, 1 fs, 2 relax, 3 fs_relax (IcRadInit)
  // icr_fmeas (set by a pgen): one more IcRadRun pass with the force on the gas; its net
  // x1 momentum deposit dm - dt rho a_ref goes to icr_dm (dt icr_dm_dt), the gas is
  // then restored
  bool icr_fmeas = false;
  DvceArray4D<Real> icr_dm;
  Real icr_dm_dt = 0.0;
  //! spherical-polar wedge: fatal on an unsupported configuration (after ImplicitInit)
  void SphericalCheck(ParameterInput *pin);
  // spherical-polar or cubed-sphere mesh: the implicit kernels take the Coordinates face
  // areas, volumes and centre distances
  bool sph_geom = false;
  // cubed sphere (sets sph_geom too)
  bool cs_geom = false;
  // cubed-sphere transverse operator (CubedInit):
  //   vetbcs     mb_bcs with BoundaryFlag::panel read as block (a seam is an open face)
  //   cs_seam   per MeshBlock: its inner_x2, outer_x2, inner_x3, outer_x3 face is a seam
  //   csg2/3    seam geometry (VetCsSeamGeom) of the x2 / x3 seam faces,
  //             (m, {dth, angm, sfoot, smid}, along-seam cell incl. ghosts, side lo/hi)
  DualArray2D<BoundaryFlag> vetbcs;
  DualArray2D<int> cs_seam;
  DvceArray4D<Real> csg2, csg3;
  // effective transverse geometry of the cs operator for kernels that read Coordinates
  // (ImplicitVimpBuild): two-point x2/x3 distances and the seam areas; x1 copied
  DvceFaceFld4D<Real> cs_dxf_eff{"vetcsdxf", 1, 1, 1, 1};
  DvceFaceFld4D<Real> cs_area_eff{"vetcsarea", 1, 1, 1, 1};
  // planar unit normals of the x2 faces (m, j = 0..n2, c) and x3 faces (m, k = 0..n3, c),
  // along +x2 / +x3 (VetCsLat0)
  DvceArray3D<Real> cs_fn2, cs_fn3;
  // the same normals at cell centres (m, j or k incl. ghosts, c): lateral axes of the cs
  // half-range ratios
  DvceArray3D<Real> cs_cn2, cs_cn3;
  void CubedInit();
  //! print the timers and positivity totals of the implicit transport
  void ImplicitReport();
  // dbg_energy_audit (output only): per solve, the change of gas energy + (c/chat) E
  // between the marks of EAuditMark, by cell class (rest, floor gas, reset), with the
  // boundary fluxes; rank 0 prints EAUDIT at mark 5 (mark 6: a discarded attempt)
  // acl_cnt_d: StageTimers clips over active cells (0 E raised, 1 |F| scaled)
  DvceArray1D<Real> acl_cnt_d;
  void EventTotals(Real *v);
  //! VET: read <rad_vet>/vet_*, check the mesh, allocate
  void ScCartInit(ParameterInput *pin);
  //! VET: short-characteristics formal solution and the uniaxial tensor (chi, n)
  void ScCartBuild();
  //! VET, vet_tensor = full: guarded D = K/J (and its ghosts) for the solve
  void ScCartTensor();
  //! VET: cost line at the end of the run
  void ScCartReport();
  //! VET, several MeshBlocks: allocate the banded buffers and neighbour tables
  void ScCartInitMB(ParameterInput *pin);
  //! VET, several MeshBlocks: exact sweep with per-layer band exchange
  void ScCartSweepMB();
  //! VET: free the multi-block state
  void ScCartFree();
  // ---- scheme = clean: the implicit stage solve (rad_vet_tasks.cpp)
  TaskStatus StageSolve(Driver *d, int stage);
  // the step's radiation (rad_vet_imex.cpp), called by Driver::Execute
  void RunStage(Driver *pdrive, const int stage);
  bool Time2Stage(Driver *pdrive, const int stage);
  void Time2Accept(Driver *pdrive);
  void BackwardEuler(Driver *pdrive);
  // the steps of StageSolve, in order (rad_vet_tasks.cpp)
  void SolveBegin();
  void StageClosure();
  SolvePass PassSetup();
  void FirstPass();
  void PassClosure(SolvePass &p);
  void PassRows(const SolvePass &p);
  bool PassSolve(SolvePass &p);
  void PassUpdate(SolvePass &p);
  bool PassSettled(SolvePass &p);
  void SolveEnd(SolvePass &p);
  void VerdictSolve(const SolvePass &p);
  void UpdateState(const SolvePass &p);
  void VerdictState(const SolvePass &p);
  void InitArrays();
  void CheckKeys(ParameterInput *pin);
  void CsGeometry();
  void SteadyStart(Driver *pdrive);
  void InitFaceFluxes();
  void RestartFaces();
  void OldVector();
  void BuildClosure();
  void PredictState(const bool use_k1);
  void RayMoments(DvceArray5D<Real> dst, const bool solve = true);
  void ClosureExtrap();
  bool ExtrapSx() const {return cln_extrap == 3 && t2_tab != 0;}
  void ClosureAfterSolve();
  void ClosureIterate();
  void StageOpacity();
  void TraceRays();
  void ClosureFromRays();
  void ClosureFromSc();
  void RayCut();
  void RayFaceSums();
  const SolveArgs *cln_gdx = nullptr;   // SolveArgs in device memory (slot 3)
  void RingSets(const Real tc, const int cyc);
  void RingMoments();
  void RingShift(const int64_t b);
  int64_t RingBlk(const Real dt) const;
  bool RingTwo() const;
  bool FhBuilt() const;
  void RingPreShift();   // dt_retry: shift before the checkpoint
  void RingJump();
  void RingSolveSet();
  void ClosureFill();
  void ClosureEddington();
  void FaceCoefficients();
  void AdvectiveCoefficients();
  void FirstIterate();
  // (beta: CLW_BT, CLW_RE only; lat false: CLW_L holds E_k's lateral sums)
  void RightHandSide(const bool beta = false, const bool lat = true);
  void PreconditionerSetup();
  void Precondition(const int src, const int dst);
  MgLevel GetMgLevel(const int l);
  void LineSolve(const MgLevel &lv, const int src, const int dst,
                 const int col, const int lim = -1);
  void CoarseSetup();
  void RingMatesInit(const int d, RingMates &rm);
  void VCycle(const int l, const int src, const int dst);
  void ApplyOperator(const int src, const int dst);
  int SolveBiCGStab();
  Real NewtonUpdate();
  Real VelocityUpdate();
  void CopyIterate();
  int CountNegative();
  int CountUnrealizable();
  void FaceBracket();
  void FluxSolve(const int e, const bool tie = true, const bool lat = true);
  void FluxMemoryKeep(const bool save);
  // ray history a failed stage's builds change: kept at the stage-loop start, put back
  // with U^n (Time2Stash, Time2Restore) for the backward-Euler redo
  struct RayKeep {
    DvceArray5D<Real> v[6];
    DvceArray2D<Real> fvd;
    std::vector<DvceArray2D<Real>> fvdm;
    std::vector<Real> fvkey;
    std::vector<int64_t> fvblk;
    int64_t rb = 0, rpb = 0, rbld = 0;
    Real rpt = 0.0, rs = 0.0, rg = 0.0, ro = 0.0;
    int rns = 0, fvn = 0, rcut = -1, rcut_n = 0;
    bool rok = false, rfix = false, rlast = false, rjump = false, rfdone = false,
         rhist = false;
  } cln_rk;
  void RayHistory(const bool save);
  void FluxMemoryInc(const bool two, const Real c1, const Real c2, const Real w1,
                     const Real w0);
  void FluxMemoryRetry(const bool save);   // dt_retry: d to / from cln_rstd
  int RingFields(int *fm, const bool file = false) const;   // ring
  // time2_face_fct's two scratch slots of iw (clean: the line-solve slots)
  int Time2FctSlot(const int q) const;
  void Time2FctHalo(const int c);
  int RingRstN() const;            // ring channels in the restart (version 4)
  // vet_rotate_avg_ema: a fill takes 2 sets; cln_rr = EMA offsets e, aux J and D11 and
  // the set count n; cln_rpd = Delta, a_ch and the next set's J and D11 (RingEma)
  int RingK() const {return vra_ema ? 2 : vra_k;}
  int RingNr(const int nf) const {return vra_ema ? nf + 3 : (vra_k - 1)*nf;}
  int RingNd(const int nf) const {return vra_ema ? nf + 3 : nf;}
  void RingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb, const int m0,
                   const bool pk, const bool have);
  // closure components in the restart (extrap history, ring), RstSlot
  // order; the cell layout keeps the first CLN_N - 1, CLN_JR is staggered only
  int RstNq() const;
  // ray intensity channels of the restart: the cs sc rays' (cln_fvi), and their last
  // build when a reflecting top mirrors it (cln_scmem: its top shell)
  int RstNfv() const;
  int RstNcl() const;   // closure_form = iterate: the closure in place (cln_s)
  void RingUnion();
  void FaceHalf(const DvceArray2D<Real> &dv, const int o, const bool zero,
                const int fo);
  DvceArray5D<Real> &RstFv() {return (cln_scmem && !cln_cssc) ? cln_top : cln_fvi;}
  void FaceDrag(const bool first, const bool vel = false, const bool lat = true);
  void FluxCsShared();
  void FluxCone(const int e, const int kd);
  void FluxBalanceProject(const int d, const int slot, const DvceArray4D<Real> &rs,
                          const DvceArray4D<Real> &rw, const bool cnt);
  void FluxTie(const int e, const int kd, const bool pt = false);
  void CsFaceTie(const int nc, const int c0, const int c1 = 0, const int c2 = 0);
  Real FluxSweepChange();
  void FluxSchur();
  void FluxSchurTile(const int gr[3], const int gs[3], const bool per[3],
                     const int ngd);
  void CsLateralStress(const int e, const bool st = false, const int part = 0);
  void WriteBack();
  void WbFct(const SolveArgs &x, const SolveArgs *px, const DvceArray5D<Real> &fp,
             const bool pfo, const DvceArray5D<Real> &pf);
  SolveArgs MakeArgs();
  Real MeshBlock0Dx(const int d) const;
  const void *DevCopy(const void *p, const size_t n, const int s);
  template <class T> const T *SolveDev(const T &x, const int s) {
    return static_cast<const T *>(DevCopy(&x, sizeof(T), s));
  }
  const SolveArgs *DevArgs(const SolveArgs &x, const int s);
  template <typename T>
  void GhostFillT(DvceArray5D<T> a, const int c0, const int nc, const bool phys,
                  const int part, const int dep = 0,
                  const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  void GhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                 const bool phys = true, const int part = 0, const int dep = 0,
                 const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  // dep: ghost layers its readers need (0 = all ng); fm[f]: components (bit n for c0 + n)
  // the ghosts of face family f take (3 at edges and corners, HaloInit)
  template <typename T> struct GSeg {
    DvceArray5D<T> a; int c0, nc; bool phys; int dep = 0;
    std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u};
  };
  template <typename T>
  void GhostExchange(const std::vector<GSeg<T>> &sg, const int part);
  template <typename T>
  void GhostPhys(DvceArray5D<T> a, const int c0, const int nc);
  void GhostDefer(DvceArray5D<Real> a, const int c0, const int nc,
                  const int dep = 0);
  void GhostFill(const WorkVec &a, const int c0, const int nc,   // (in one array)
                 const int part = 0, const int dep = 0,
                 const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  void GhostDefer(const WorkVec &a, const int c0, const int nc, const int dep = 0);
  void GhostFlush();
  std::array<Real, 2> GlobalDots(WorkVec w, const int nc, const int n,
                                 const std::array<int, 4> ab);
  Real GlobalReduce(Real v, const bool max);
  Real SumCells(const DvceArray4D<Real> &a);
  void ClosureSphere();

  Real cln_c0 = 0.1;            // closure_w_c0: C1 below which w = 0
  Real cln_c1 = 0.5;            // closure_w_c1: C1 above which w = 1
  Real cln_tol = 1.0e-10;       // solver_tol: relative residual of the linear solve
  // clean_init_cnorm: the steady start's per-cell test (0 off)
  bool cln_insteady = false;    // in SteadyStart's passes
  int cln_maxit = 800;          // solver_maxit: BiCGStab iterations
  std::map<std::array<int, 7>, DvceArray1D<int>> cln_pcrl;   // PCR lines of a colour
  int cln_optn = 3;             // gas_opac_passes: refreshed in passes < N
  // the operator's rows (CLN_ST per active cell, PreconditionerSetup)
  DvceArray5D<Real> cln_st;
  //   (cubed sphere: CLN_CSY coefficients of the lateral sums)
  DvceArray5D<Real> cln_csy;
  bool cln_csyok = false;       //   cln_csy is of the current closure (no rebuild)
  bool cln_ysum = false;        // clean_stag_lat_sums: the sums Y off the cubed sphere
  DvceArray5D<Real> cln_yl;     //   (in their own array: CLW_L is the F solve's scratch)
  int cln_fsweeps = 1;          // clean_stag_fsweeps: line-solve sweeps of the F solve
  bool cln_fschur = false;      // clean_stag_fschur: the F lines exact across MeshBlocks
  int cln_maxnl = 50;           // newton_maxit: Newton iterations
  Real cln_tol_nl = 1.0e-11;    // newton_tol: Newton tolerance (relative dE, dT, dv)
  bool cln_c1b = false;         // newton_closure_budget
  bool cln_bmn = false;         // newton_maxnorm
  bool cln_wbxl = false;        // update_exchange_local
  Real cln_flrun = 0.0;         // the last write-back's implicit_pos_floor creation
  int cln_flrid = -1;           // its cycle and t2s, 4 ncycle + t2_solve
  bool cln_wbmc = false;        // update_mcap
  int cln_wbfdg = 0;            // update_floor_diag (diagnostic)
  int cln_dbginf = -1;          // dbg_inf (test): the cycle of the injected inf
  Real cln_raymin = 1.0e-3;     // closure_ray_min: rays describe a cell when J >= this E
  // clean_ray_cut (RayCut): no rays below the shell where every deeper cell is
  // diffusive (|grad S|/(chi S) < ray_cut_r, none thinner than
  // vet_col_lat_cut_thin), less clean_ray_cut_margin optical depths, at
  // tau_top >= clean_ray_cut_taumin; rises only after ray_cut_wait builds
  bool cln_rcut_on = false;
  int cln_rcut_wait = 10, cln_rcut = -1, cln_rcut_n = 0;
  Real cln_rcut_r = 0.003;
  ParameterInput *pin_rcut_ = nullptr;
  bool cln_k1ok = false;        // K1 = f(U^n) valid (an accepted solve wrote it): the
                                // predictor and the budget's time error also in BE steps
  static constexpr int cln_maxpic = 12;   // passes that update the lagged slopes
  bool cln_gflux = false;       // clean_stag_gflux: rays' face flux as flux factors
  bool cln_rfdone = false;      // its fill taken
  bool cln_rmdone = false;      // its start taken
  void FaceDonor(const bool first, const bool lat = true, const bool hold = false);
  void WorstCell(const Real dnl, const Real dnl0);
  bool cln_freal = false;       // clean_stag_freal: realisable sb, restart faces
  // clean_stag_xdonor: cross flux of the iterate as a donor-cell coefficient of its face
  bool cln_xd_on = false;
  int cln_xqs = -1;   // (3 slots: qL >= 0, qR <= 0, the lagged explicit flux; CLF_XC: C(E_k),
                    //  C's last change)
  int cln_lza = 0;          // clean_cs_lz_abs
  bool cln_lzbook = false;  // (the L_z ledger: cubed sphere, or 3-D spherical polar)
  bool cln_lzt = false;     // lz_telescope
  bool cln_latwb = false;   // flux_balance_lat
  bool cln_thinf = false;   // face_thin
  DvceArray1D<Real> cln_mclc;    // (the Newton passes' momentum cap: the terms it moved)
  bool cln_fwb = false;     // clean_stag_fwb
  DvceArray3D<Real> cln_fsjt[3];   // (layout tile: the F join per tile, FluxTileArray)
  std::array<RingMates, 3> cln_fwr;   // (fwb: the lines' MeshBlocks)
  DvceArray1D<Real> cln_fwbc;    // (fwb: the faces, their momentum sums |dp|, x1 dp)
  bool cln_fcone = false;   // flux_cone
  bool cln_fcgas = false;    // flux_cone_gas
  bool cln_fclz = false;    // flux_cone_lz
  bool cln_fcemx = false;   // flux_cone_emax
  bool cln_lzsp = false;    // lz_spin
  DvceArray5D<Real> cln_fcx;   // (lz_spin: fcone's cut withheld, by face)
  DvceArray2D<Real> cln_cx[2];   // (lz_spin, Cartesian: the cell centres x, y)
  bool cln_fwbring = false;   // clean_stag_fwb_ring
  DvceArray1D<Real> cln_fconc;   // (fcone: the faces scaled, their momentum sum)
  int cln_xdfz = 0;         // face_xdonor_freeze
  bool cln_xdh = false;     // (its hold in this pass)
  // clean L_z ledger (clean_cs_lz_abs = 2; VetLzKey): cln_lzl the L_z the gas does not
  // take, booked as the radiation's; its sums as the gas momentum's slopes K1 (cln_lzs1),
  // K2 (cln_lzs2), stage-start increment (cln_lzc) and stage-1 part pending (cln_lzp).
  // 3-D spherical polar: the L_z that crossed the radial boundaries.
  Real cln_lzl = 0.0, cln_lzs1 = 0.0, cln_lzs2 = 0.0, cln_lzc = 0.0, cln_lzp = 0.0;
  ParameterInput *cln_lzpin = nullptr;
  void LzBook(const Real s, const bool start);
  Real LzTotal();   // gas + radiation + ledger L_z
  bool cln_gpic = false;   // a stalled solve's passes without the gas Newton slope
  Real cln_ngnsw = 0.0;   // solves switched so far
  int cln_fdg = -1;             // staggered bracket: face slot of its lagged flux DG
  bool cln_opt = false;         // clean_stag_opac_update: opacities at the iterate's T
  bool cln_consist = false;     // closure_consistency: bracket C = F_ray - f(Ur), cell
  bool cln_closit = false;      // closure_form: predict | iterate
  int cln_extrap = 0;           // closure_form = extrap (1) | extrap2 (2)
                                // | extrap1 (3, cell default)
  bool cln_edd = false;         // closure_form = eddington (D = I/3, no rays)
  bool cln_sphere = false;      // closure_form = sphere (spherical polar)
  Real cln_sph_r = 0.0;         // closure_sphere_radius
  static constexpr Real cln_trtol = 1.0e-13;   // residual stop of the exchange root
  int cln_wo_nc = 0;            // the slots of cln_w.o (BindWorkVecs)
  static constexpr int cln_hmax = 16;   // most variables of one ghost exchange
  bool cln_budget = false;      // newton_conv = budget (staggered default): StageSolve
  bool cln_bon = false;         // this solve takes it (Newton passes with a valid K1)
  Real cln_tolb = 0.0;          // its tolerance of the pass (kappa x time error)
  Real cln_bkap = 1.0;          // newton_conv_kappa
  Real cln_bmax = 1.0e-5;       // newton_conv_max
  DvceArray5D<Real> cln_wci;    // E before a backward-Euler solve's last pass
  Real cln_ewr = -1.0, cln_eta = 0.1;   // Eisenstat-Walker: last |r0|_S, forcing
  int cln_verbose = 0;          // verbose: 1 = also the per-solve diagnostics
  int cln_pcsweep = 0;          // solver_sweep: +1 a 2nd red pass, +2 a post-sweep
  bool cln_steady = false;      // steady_init = steady: radiation relaxed before step 1
  int cln_steady_n = 50;        // steady_n: most pre-solve passes
  Real cln_steady_dtfac = 1.0e6, cln_steady_tol = 1.0e-8;  // steady_dtfac, _tol
  std::vector<DvceArray5D<Real>> cln_mg;          // coarse levels of the multigrid
  std::vector<std::array<int, 3>> cln_mgn;        //   and their cells along x1..x3
  std::array<int, 3> cln_mgt = {0, 0, 0};       //   its tiles' cells (layout_tile)
  int cln_mgnl = 2;             //   most coarse levels (solver_mg_levels)
  bool cln_fsolve = false;      // staggered: the face F of a separate F solve
  bool cln_fkin = false;        // staggered fsolve: streaming share of the rays' kinetic
                                //   momentum flux
  bool cln_dok = true;          // staggered fkin: the last solve converged (memory d)
  bool cln_fwork = false;       // staggered fsolve: also in the Newton loop (the work W)
  bool cln_fsk = false;         //   CLF_N holds the F solve's F of the iterate
  bool cln_upw = false;         // staggered: donor-cell advective flux, lagged rest
  bool cln_split = false;       // staggered: energy flux apart from the face momentum
  int cln_nr1 = 0, cln_nr2 = 0;  // staggered: cells E < 0, faces |F| > c E_donor
  Real cln_fe = 0.0, cln_fv = 0.0;  // round-off levels of the Newton changes dE, dv
  bool cln_feedback = false;    // gas coupling active
  bool cln_init = false;        // arrays allocated
  bool cln_faces_set = false;   // starting face fluxes set
  bool cln_rstold = false;      // restart from a file without the clean state
  bool cln_havec = false;       // a closure has been built
  bool cln_cip = false;         // a closure is in place (closure_form = iterate: the
                                // last iterate's, restart state)
  // opac holds predicted opacities during the formal solution
  bool cln_opsaved = false;
  Real cln_tc = -1.0e300;       // time of the closure build
  Real cln_dtb = 0.0;           // dt of the closure build
  DvceArray5D<Real> cln_s;      // closure (CLN_*)
  bool cln_cssc = false;        // the cubed sphere, clean_rays = sc (vet_gd sweep)
  bool cln_rfh = false;         // cs cell half ranges in the ring (CLN_HO)
  DvceArray5D<Real> cln_fvi;    // the sc rays' intensities (m, d, k, j, i)
  DvceArray2D<Real> cln_fvd;    // the swept directions (x, y, z, weight; RingSets)
  int cln_fvn = 0;              // their number
  std::vector<Real> cln_fvkey;  // the set (block key) in each slot of a rotating union
  // ring: the rotating union swept as the current and the next set only, the
  // K - 1 older sets carried as offsets of their ray moments (RingMoments)
  int cln_ring = 0;             // ring = live (1: + the next set) | jump (2)
  bool cln_rfuse = false;       // ring_fuse: live ring's two sets in one sweep
  bool cln_rfuse_said = false;
  bool cln_ring_named = false;  // ring named in the input
  bool cln_rlast = false;       // jump: the build closes its block (the next set swept)
  bool cln_rjump = false;       // jump: cln_rnew holds the next block's moments
  DvceArray5D<Real> cln_rnew;   // jump: moments of block b+1's window at the boundary
  bool cln_rfix = false;        // no build of this cycle shifts (RingPreShift)
  bool cln_rok = false;         // the ring holds the offsets of the sets before cln_rb
  bool cln_rhist = false;       // the restart file's EMA history is converged
  int64_t cln_rb = 0;           // the current set of the ring (rotation block)
  int64_t cln_rpb = 0;          // the block of the pending Delta
  Real cln_rpt = 0.0;           // its time (the block boundary)
  Real cln_rs = 0.0;            // the cross-fade fraction of the build
  int cln_rns = 0;              // the swept slots of the build
  int64_t cln_rbld = 0;         // the block of the build
  std::vector<int64_t> cln_fvblk;   // the set in each slot
  std::vector<DvceArray2D<Real>> cln_fvdm;   // the directions with only slot q weighted
  DvceArray5D<Real> cln_fh;     // the sc rays' face half-range sums (FaceHalf)
  DvceArray5D<Real> cln_gdgl;   // their ghost layers (RayFaceSums, scratch)
  // the slots' moments g = M/J (CLN_J: J); a build's scratch
  DvceArray5D<Real> cln_rl;
  DvceArray5D<Real> cln_rr;     // offsets g_{b-j} - g_b, J_{b-j}/J_b (< 0: none), K-1
  DvceArray5D<Real> cln_rpd;    // the pending Delta g_{b+1} - g_b, J_{b+1}/J_b
  bool cln_rpfree = false;      // its memory is scratch of the build (RayFaceSums)
  Real cln_rguard_sum = 0.0, cln_roff_sum = 0.0;   // guard falls, off-boundary shifts
  DvceArray5D<Real> cln_v;      // closure fill, bracket, E < 0 flags (m,3,k,j,i)
  DvceArray5D<Real> cln_vf;     // cubed sphere: the F solve's cells, faces (m,18,k,j,i)
  DvceArray5D<Real> cln_opb;    // opacities of the stage, kept during the closure build
  // closure_form = extrap(2): ray moments (CLN_* layout, not normalised by J) of the
  // converged states at t^n (cln_mn), t^{n-1} (cln_mp), the end of the last step solved
  // (cln_mc) and the stage-1 solution (cln_m1), with their times and validity
  DvceArray5D<Real> cln_mn, cln_mp, cln_mc, cln_m1;
  Real cln_tn = -1.0e300, cln_tp = -1.0e300, cln_tm = -1.0e300, cln_t1 = -1.0e300;
  bool cln_have_n = false, cln_have_p = false, cln_have_c = false, cln_have_1 = false;
  bool cln_cand_be = false;     // the step that made cln_mc ended with backward Euler
  int cln_neint = 0;            // cells of the write-back left with gas eint < 0
  int cln_nguard = 0;           // cells whose extrapolated moments fell back to t^n
  Real cln_nguard_sum = 0.0, cln_next_sum = 0.0;   // guard falls, extrapolated closures
  WorkVec cln_w;                 // work vectors (CLW_*)
  FaceVars cln_fc1, cln_fc2, cln_fc3;   // face variables (CLF_*)
  DvceArray5D<Real> cln_fref;   // force_reference = clean: reference force (m,3,k,j,i)
  // staggered: cell flux of the faces' energy flux Phi (m,3,k,j,i; output vet_phi); u0 F
  // holds the face momenta's
  DvceArray5D<Real> cln_phic;
  DvceArray5D<Real> cln_rstd;   // staggered fkin: the restart's d and its K1 (m,6,k,j,i)
  bool cln_fref_set = false;    // cln_fref holds the force of the first coupled solve
  bool cln_mdcap = false;       // wb_arad_md: the write-back captures cln_fref (steady)
  bool cln_vzero = false;       // the iterates keep the gas momentum (steady start)
  // cubed sphere: face areas, centroid distances across the faces, normals and centroids
  DvceArray4D<Real> cln_ga[3], cln_gd[3];
  CsGeoArrays cln_gn;
  DvceArray3D<Real> cln_tg;   // spherical polar: SolveArgs::tg
  // vet_gd_theta_seam_sym: SolveArgs' geometry with the theta-periodic seam's ghost rows
  // their periodic images' (SeamImage): vol, area x1f / x2f, dxface x2f / x3f,
  // x2v, tg
  DvceArray4D<Real> cln_sgv, cln_sga0, cln_sga1, cln_sgd1, cln_sgd2;
  DvceArray2D<Real> cln_sgx;
  DvceArray3D<Real> cln_sgtg;
  void SeamImage();
  struct SolveHalo {            // ghost exchange of nc components: scratch array and
    DvceArray5D<Real> a, ca;    // its coarse dummy
  };
  std::map<int, SolveHalo> cln_halo;   // by nc
  DvceArray1D<char> cln_xdev;   // device copies of kernel arguments (DevCopy)
  // exchange object of every clean ghost fill (a fill sends only its own variables);
  // fills deferred to the next exchange (GhostDefer)
  std::unique_ptr<MeshBoundaryValuesCC> cln_hpb;
  std::vector<GSeg<Real>> cln_hpend;
  DvceArray5D<Real> cln_hbuf;   // storage of every scratch array (<= 32)
  // clean_halo_direct: cell lists read once off cln_hpb's exchange (HaloInit);
  // on-rank ghosts by one copy, one message per neighbour rank, receives posted first
  struct SolveHaloList {
    int state = 0;              // 0 lists to make at first use, 1 on, -1 while made
    int n1 = 0, n2 = 0, n3 = 0, nl = 0, ns = 0, nr = 0, par = 0;
    DvceArray1D<int> ld, ls, sc, rc;   // on-rank ghost <- its source; sent; received
    std::vector<int> rank, soff, slen, roff, rlen;
    // lists ordered by ghost depth 1..ng; per message and depth d the entries of
    // depth <= d (slend/rlend[s ng + d - 1], nld[d - 1]); rd, ldp, sdp: depth of each
    // received / on-rank / sent entry; dcur: the exchange's depth
    int ng = 0, dcur = 0;
    std::vector<int> slend, rlend, nld;
    DvceArray1D<int> rd, ldp, sdp;
    // device buffers on HIP too (GPU-aware MPI; host staging is copy bound on MI300A)
    DvceArray1D<Real> sbuf, rbuf;      // rbuf: two halves, alternate exchanges
    // face family of each on-rank, sent and received ghost (hsf, hrf on the host); per
    // mask set each entry's buffer offset and the message lengths by depth
    DvceArray1D<int> lf, sf, rf;
    std::vector<int> hsf, hrf;
    struct Mask {
      DvceArray1D<int> so, ro;
      std::vector<int> slm, rlm;
    };
    std::map<std::array<uint32_t, 5>, Mask> mask;
#if MPI_PARALLEL_ENABLED
    MPI_Comm comm = MPI_COMM_NULL;
    std::vector<MPI_Request> rq, sq;
#endif
  };
  SolveHaloList cln_hdl;
  void HaloInit(const int n3, const int n2, const int n1);
  template <typename T> void HaloPost(const std::vector<GSeg<T>> &sg);
  template <typename T> void HaloFinish(const std::vector<GSeg<T>> &sg);
  template <typename T>
  std::array<uint32_t, 5> HaloKey(const std::vector<GSeg<T>> &sg);
  SolveHaloList::Mask &HaloMask(const std::array<uint32_t, 5> &key);
  Real cln_nsolve = 0.0, cln_nit = 0.0, cln_nfail = 0.0, cln_nneg = 0.0;

  bool vray_ready{};            // VetRayInit has run
  // vet_scatter: electron scattering, the opacity split (otab.scat, opac VET_OP_S) and
  // VetScatterMix
  bool vscat = false;
  void VetRayInit();           // checks, buffers, the direction set (first call)
  void VetRayReport();

  // ---- the source (ln chi, ln S) of the vet_gd rays (rays/sc_source.cpp)
  // Athena++ <mesh>/x2rat, x3rat named != 1 (AthenaK grids are uniform): VetRayInit
  // refuses them
  std::string vray_rat_keys;
  int vray_icut = 0;           // first shell of the sweep (RayCut)
  Real vray_ncall = 0.0, vray_nclamp = 0.0;
  DvceArray5D<Real> vray_cs;   // (m, 2, k, j, i): ln chi, ln S (ghosts exchanged)
  DvceArray5D<Real> vray_cs_c;
  DvceArray1D<Real> vray_cnt;  // (1): gd: lateral reads clamped; Cartesian: I < 0 clipped
  MeshBoundaryValuesCC *pbval_vs = nullptr;   // vray_cs
  void VetRaySource();   // ln chi, ln S of a vet_gd build
  void VetRayExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                      MeshBoundaryValuesCC *pb);

  // ---- vet_gd (rays/sc_sph_sweep.cpp): Eddington tensor from an SC sweep along globally
  // fixed (Cartesian) HEALPix directions; sp wedge / cubed sphere
  bool vgd_on = false;         // closure = sc on spherical polar, cubed sphere
  // HEALPix nside of a set that names none (12 nside^2 directions)
  static constexpr int vgd_nside = 3;
  int vgd_n = 0;
  int vgd_ls = 0;              // vet_gd_ls = 6/8/10/12: level-symmetric LQ_N set instead
  double vgd_zper = 0.5*M_PI;  // z-rotation period of the set (ScInit)
  int vgd_gl_nmu = 0;          // vet_gd_gl_nmu > 0: product set GL(nmu in n_z) x nphi
  int vgd_gl_nphi = 12;        // vet_gd_gl_nphi (C_nphi about z)
  bool vgd_gl_stag = true;     // vet_gd_gl_stagger: alternate rings half a step
  int vgd_rot = 0;             // vet_gd_rotate_every: z-rotation every N cycles
  // vet_gd_rotate_dt > 0: rotate in physical time, block
  // b = floor((t + dt/2)/dt_rot); steps
  // land on block boundaries (RotAlignDt)
  Real vgd_rot_dt = 0.0;
  // vet_gd_rotate_dt = auto: vgd_rot_auto_fac x the first dt, recorded in pin_rot_
  bool vgd_rot_auto = false;
  // the unnamed default resolves from the uncapped first CFL dt (Mesh::dt_uncapped)
  bool vgd_rot_unc = false;
  Real vgd_rot_auto_fac = 20.0;
  ParameterInput *pin_rot_ = nullptr;
  double vgd_rot_nalign = 0.0;   // steps shortened by RotAlignDt (diagnostic)
  int64_t vgd_rot_lastblk = INT64_MIN;
  void RotNewBlock();
  // a Cartesian rotating set without the ring, the first step of a new rotation block: a
  // correction solve (Time2RotCorrect) of size dt_c = kT2G dt from U^n with the new set
  // turns K1 into K1 + (I - dt_c A_new)^{-1} (f_new - f_old)(U^n)
  bool t2_rotcorr = false;
  double vgd_rot_ncorr = 0.0, vgd_rot_ncfail = 0.0, vgd_rot_ncpass = 0.0;
  bool Time2RotCorrect(Driver *pdrive, Real dt);
  Real vgd_alpha = -1.0;       // current z-angle (-1: tables not built)
  // vet_gd_rotate_axis = z | 3d (default z on sp, 3d on the cubed sphere). 3d: the set of
  // block b is R_b n, R_b a 3-D rotation about a quasi-random axis (ScRot3); vgd_alpha
  // is then the block key, not an angle
  bool vgd_rot3 = false;
  void ScRot3(const Real key, double r[9]) const;
  std::vector<double> ScRotate(const Real alpha) const;
  // vet_rotate_avg = K >= 2 (rays/sc_cart_rotavg.cpp): the closure is the mean over the
  // last K
  // rotated direction sets. A closure registers its normalised fields (ratios to J); the
  // mean is kept as offsets d_j = g_{b-j} - g_b, each from two sweeps of the same state:
  // used = g_b(now) + sum_j d_j / K. Needs a rotating closure.
  static constexpr int VET_RA_NF = 16;
  struct RotAvgTab {           // registered fields (POD: captured by the kernels)
    int nf = 0;
    int arr[VET_RA_NF]{}, comp[VET_RA_NF]{}, kind[VET_RA_NF]{};
    Real lo[VET_RA_NF]{}, hi[VET_RA_NF]{};
    int farr = -1, fcomp = 0;  // validity flag field: invalid where < 0 (-1: none)
    int jf = -1;               // J field (index into the table) of kind-3 fields
    int jw = -1;               // the ratio field weighting the sets (-1: equal)
  };
  int vra_k = 0;               // vet_rotate_avg (K; 0/1 off)
  // vet_rotate_avg_ema = N: the K-window mean replaced by an exponential moving average
  // over the rotated sets, weight a = 1/N per block, in the same offset form (one slot
  // e = avg - g_b(now), r = Jbar/J_b, + an aux slot; vra_k = 3)
  bool vra_ema = false;
  Real vra_ema_n = 0.0;
  void RotAvgEmaInit(const int ilo);
  void RotAvgEmaShift(const int ds);
  void RotAvgEmaChange(const int ilo);
  static constexpr Real kEmaEps = 0.05;   // the change test's relative change
  bool vra_ema_chstep = false; // vet_rotate_avg_ema_change_form = step | linear
  // curved meshes (clean ring): a = a_ch alone (no 1/N floor, no 1/(n + 1)
  // count): a static state keeps its closure
  bool vra_ema_wb = false;
  void RotAvgEmaApply(const int ilo);
  // vet_gd_rotate_xfade (needs vet_gd_rotate_dt > 0): clean ring live, the closure
  // cross-faded from the current to the next set over each block (RingSets); both
  // sets' tables on the device (ScXfTables)
  bool vgd_xf = false;
  // vgd_xff: both sets of a cross-fade build in one sweep (one launch per shell, one
  // message per partner), bitwise the two sweeps. The second set keeps its own
  // intensities (vgd_i_xb), shell and halo lists. Falls back where ScXffOk fails.
  bool vgd_xff = false;
  bool vgd_xff_now = false;      // this sweep carries the next set too (ScSweep)
  bool vgd_xff_prep = false;     // ScSweep: the set's lists only, no passes
  DvceArray1D<int> vgd_shl_xb[2];
  int vgd_shn_xb[2] = {0, 0};
  Real vgd_shl_alpha_xb = -3.0;
  DvceArray1D<char> vgd_bbufx;     // fused shell closure (device copy)
  DvceArray1D<int> vgd_cslst_xb;   // cs: the other set's pass lists
  DvceArray1D<Real> vgd_xff_sb, vgd_xff_rb;   // merged messages (sp)
  double vgd_xff_n = 0.0;          // fused builds (diagnostic)
  double vgd_xff_nmsg = 0.0;       // merged exchanges (diagnostic)
  double vgd_xff_nfb = 0.0;        // builds that fell back to two sweeps (diagnostic)
  bool ScXffOk();
  void ScXffSwap();
  void ScEs2Make();   // second stream of the cubed sphere's early send
  void ScGbSwap() {   // the compact global band's kept lists go with their set
    std::swap(vgd_gbcf, vgd_gbcf_xb);
    std::swap(vgd_gbcb, vgd_gbcb_xb);
    std::swap(vgd_gbrec, vgd_gbrec_xb);
  }
  bool ScHlXff(const int i, const int ws);
  bool ScGbHaloX(const int i, const int ws);   // vgd_raw: both sets, one message
  // early send (ScCsSweep): the shell's send cells first on a second high-priority
  // stream; the host exchanges while the shell kernel (skipping them) runs
  DevExeSpace vgd_es2;
  bool vgd_es2_made = false;
  DvceArray4D<int> vgd_est[2];     // per set (m, kb, jb, d): stamp of the early cells
  int vgd_esid = 0;                // stamp of the current early shell
  Real vgd_xf_alt = -9.0;        // angle of the tables in the second device set
  void ScXfTables(const Real alpha);
  bool vra_on = false;
  RotAvgTab vra_tab;
  DvceArray5D<Real> vra_a0, vra_a1;   // closure's field arrays (set per call)
  DvceArray3D<Real> vra_c;     // column field (m, k, j) at the top cell (arr 2)
  DvceArray5D<Real> vra_d;     // (m, (K-1)(nf+1), k, j, i) interior: offsets + valid
  int64_t vra_gen = 0;        // bumped by every writer of vra_d
  int vra_head = 0;            // ring slot of d_1 (the previous set)
  int vra_state = 0;           // 0: no history yet (init at the next build), 1: ready
  int vra_pend = 0;            // 1: init slots raw, 3: Delta of the coming boundary held
  bool vra_raw = false;        // a snapshot sweep: the closure leaves its fields raw
  double vra_nsweep = 0.0, vra_nrot = 0.0, vra_ninit = 0.0;
  Real vra_time = 0.0;
  void RotAvgSnap(const int slot, const int ilo);
  void RotAvgUpdate(const int ilo);
  void RotAvgApply(const int ilo);
  void RotAvgParse(ParameterInput *pin);
  void RotAvgDelta(const int slot, const int ilo);
  void RotAvgShift(const int ds);
  bool RotLastStep() const;    // the step is the last of its rotation block
  // history in the restart file (kVetRaRstMagic, rays/sc_cart_rotavg.cpp)
  static constexpr int VET_RA_RST_NH = 11;   // int64 header entries (10)
  int RotAvgNf() const;
  int RotAvgRstNslab(const int ncell) const;
  int RotAvgRstI0() const;     // first written radial cell
  void RotAvgRstHdr(int64_t *h) const;
  void RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb, const int nsc);
  bool RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src, const int nmb);
  bool RotAvgRstCheck();
  int vra_rst_nf = -1;         // nf of a history loaded from the file (-1: none/checked)
  int64_t vra_blk = INT64_MIN; // block of the last build
  bool vra_last = false;       // this build is on the last step of its block
  DvceArray5D<Real> vra_bk_vc;              // vet_sc: vet_cell kept over the next-set
                                            // sweep of a last build
  void ScSweepPost();       // sweep, moments, ScPost
  bool vgd_tbath = false;   // the top (x1max Marshak) bath enters the sweep
  std::vector<double> vgd_base;   // unrotated set (x, y, z, w)
  DvceArray2D<Real> vgd_dir;   // (d, 0..3): n_x, n_y, n_z, weight (sum 1)
  // (m, d, k, j, i) intensities, ghosts lagged, as a ragged per-shell band: shell i keeps
  // only the lateral band it needs (depth wi(i) <= vgd_w) and only on sides whose
  // neighbour is off-rank. Index space (m, d, k, j, i) with k, j in the band space of
  // depth vgd_w; one slab per (m, i), d fastest, then k, then j. An access outside the
  // slab goes to a dummy element and is counted (ScRag::oos, must stay 0).
  struct ScRag {
    DvceArray1D<Real> d;
    Kokkos::View<int64_t**, DevMemSpace> off;  // (m, i): first element of slab (m, i)
    DvceArray1D<int> wi;                       // (i): band depth of shell i
    DvceArray2D<int> sb;                       // (m, 4): band on side k-, k+, j-, j+
    DvceArray1D<Real> oos;                     // (0): out-of-slab accesses
    int nmb = 0, n = 0, nx2 = 0, nx3 = 0, w = 0, c1 = 0;
    int64_t dummy = 0;
    Kokkos::View<int64_t**, DevMemSpace> offc;   // (m, i): compact slabs (vet_gd_band_ring)
    int64_t nth = 0;                             // thick shells' values (ScRagThick)
    // is (k, j) of shell i of block m stored?  (band index space of depth w)
    KOKKOS_INLINE_FUNCTION
    bool Has(const int m, const int k, const int j, const int i) const {
      const int wl = wi(i);
      const int k0 = w - sb(m,0)*wl, k1 = w + nx3 + sb(m,1)*wl;
      const int j0 = w - sb(m,2)*wl, j1 = w + nx2 + sb(m,3)*wl;
      return (k >= k0) && (k < k1) && (j >= j0) && (j < j1);
    }
    KOKKOS_INLINE_FUNCTION
    Real &operator()(const int m, const int v, const int k, const int j,
                     const int i) const {
      const int wl = wi(i);
      const int bk0 = sb(m,0)*wl, bj0 = sb(m,2)*wl;
      const int c3 = nx3 + bk0 + sb(m,1)*wl, c2 = nx2 + bj0 + sb(m,3)*wl;
      const int kk = k - (w - bk0), jj = j - (w - bj0);
      if (kk < 0 || kk >= c3 || jj < 0 || jj >= c2) {
        Kokkos::atomic_add(&oos(0), 1.0);
        return d(dummy);
      }
      return d(off(m,i) + v + static_cast<int64_t>(n)*(kk + static_cast<int64_t>(c3)*jj));
    }
    int extent_int(const int r) const {
      return (r == 0) ? nmb : (r == 1) ? n : (r == 2) ? (nx3 + 2*w) : (r == 3) ?
             (nx2 + 2*w) : c1;
    }
    size_t extent(const int r) const {return static_cast<size_t>(extent_int(r));}
  };
  using ScIView = ScRag;
  ScIView vgd_i;
  ScIView vgd_i_xb;   // the other cross-fade set's intensities
  void ScRagAlloc(ScRag &a);
  int vgd_rag_cut = 0;   // shells below: one slab for both fused sets (ScRagThick)
  bool vgd_rag_said = false;
  void ScRagThick();
  // memory arena (with vet_gd): the intensity arrays and the solve scratch dead across a
  // ScBuild in one allocation
  bool marena_on = false, marena_sealed = false;
  DvceArray1D<Real> marena;   // intensities, then (first checkpoint) solve scratch
  int64_t marena_sw = 0;      // Reals of the intensity arrays bound so far
  bool marena_bound = false;  // MemArenaBind done
  int64_t marena_hyd = 0;     // Reals of the hydro's fluxes at its start
  void MemArenaHydro();       // zeroes those (Hydro::Fluxes)
  bool marena_k2 = false;     // t2k2, t2pf2 in it (zeroed in stage 1's WriteBack)
  void MemArenaBind();        // the solver's and the hydro's scratch arrays over it
  bool MemArenaIn(const Real *p) const {
    return p >= marena.data() && p < marena.data() + marena.extent(0);
  }
  void MemArenaInit(ParameterInput *pin, int nmb, int niw, int n3, int n2, int n1);
  bool MemArenaRays(ScRag &a, int64_t n);
  void MemArenaGrow(int64_t need, bool keep = true);
  void MemArenaSeal();
  void BindWorkVecs();   // cln_w.o, .s at the first solve (mem_arena)
  // vet_gd_band_ring (MPI band with exact halo lists, or cs compact lists): the ragged
  // array keeps each shell's interior only and lends a shell one of VGD_RING_K slabs of
  // its full band depth while the sweep needs its band (outward pass: interior copied in,
  // the saved inward-exchange band put back). Only off(m, i) / wi(i) switch, so every
  // value read is the dense band's (bitwise).
  static constexpr int VGD_RING_K = 3;
  int vgd_rg_k = VGD_RING_K;          // slots of this allocation
  struct ScRingOp {
    int kind = 0;      // 1 save (inward band -> store), 2 release, 3 acquire
    int i = -1, slot = -1;
    bool copyin = false, restore = false;
  };
  bool vgd_ring = false;              // ring layout allocated (both arrays)
  bool vgd_ring_par = true;   // vet_gd_band_ring
  bool vgd_ring_act = false;          // the ring runs in this sweep
  bool vgd_rg_two = false;            // this sweep carries a second array (vgd_i_xb)
  bool vgd_rg_xb = false;             // the other xfade set's array
  std::vector<int64_t> vgd_rg_sn2;    // (i) the other set's saved count
  bool vgd_rg_lastx = false;          // the last halo used an exact list
  int vgd_rg_lastx_e = -1;
  std::vector<int64_t> vgd_rg_soff;   // (s nmb + m): slot offsets (host)
  std::vector<int> vgd_rg_wfull;      // (i): a shell's band depth when lent a slot
  std::vector<int> vgd_rg_state;      // (i): -1 compact, else slot
  std::vector<int> vgd_rg_sown;       // (s): shell holding slot s, -1 free
  std::vector<int64_t> vgd_rg_so, vgd_rg_sn;   // (i): save store offset, count/array
  std::vector<int> vgd_rg_sk;         // (i): save kind: 1 exact list, 2 dense band, 3 raw
  std::vector<int> vgd_rg_se;         // (i): exact list record (3: the raw key)
  // vgd_raw: band positions (as exact lists) each raw exchange writes, keyed as vgd_gbcf
  // (ScRawRec)
  std::map<int, DvceArray1D<int>> vgd_gbrec, vgd_gbrec_xb;
  int ScRawKey(const int i, const bool inw, const int ws) const;
  DvceArray1D<int> ScRawRec(const int key);
  int64_t vgd_rg_stop = 0;            // store fill this sweep
  DvceArray1D<Real> vgd_rg_store;     // saved inward band values (main, then twin)
  Kokkos::View<int64_t**, DevMemSpace> vgd_rg_dsoff;   // (s, m)
  std::vector<ScRingOp> vgd_rg_pend;
  Real vgd_rg_nk = 0.0, vgd_rg_ndense = 0.0;
  void ScRingFlush();
  void ScRingSweepBegin(const bool two);
  void ScRingSweepEnd();
  int ScRingFreeSlot(const int avoid) const;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_wall;    // (m, k, j)
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_map;    // (m, k, j, d)
  DvceArray4D<Real> vgd_mr;    // (m, k, j, d): n . r_hat of the value stored there
  // vet_gd_wall_interp: wall ghosts from the 3 nearest same-branch directions
  // (m, k, j, d, 3)
  bool vgd_wint = false;
  // vet_gd_wall_sym (needs _wall_interp): the set turned so a mirror plane is the wedge
  // mid-plane, phi-wall ghosts linear in azimuth on their ring (ScMirPick)
  bool vgd_wsym = false;
  bool vgd_wsx = false;       // and its phi wall map inexact: both-frame seam faces
  Real vgd_wsw = 0.0;         // phi period W
  DvceArray1D<int> vgd_mir;   // (n): mirror image of each direction
  // vet_gd_theta_seam_sym (needs _wall_interp): clean theta seam faces in both frames
  bool vgd_wst = false;
  DvceArray1D<int> vgd_mirz;  // (n): the equatorial mirror image of each direction
  ScLG vgd_lg;   // the vet_gd sweep geometry from the global index
  Kokkos::View<int *****, LayoutWrapper, DevMemSpace> vgd_m3;
  DvceArray5D<Real> vgd_w3;
  // vet_gd_theta_wall = periodic | mirror: mirror makes the theta walls specular; a theta
  // ghost at depth g reads the interior row at depth g across the wall, direction
  // (n.r, n.theta, n.phi) -> (n.r, -n.theta, n.phi), from up to 3 same-branch directions
  // (vgd_wall = 2, vgd_map / vgd_m3 / vgd_w3); ln chi, ln S from the mirror rows
  // (ScThMirCs)
  bool vgd_thmir = false;
  bool cln_scmem = false;      // clean sc rays, reflecting top: its mirror reads cln_top
  DvceArray5D<Real> cln_top;   // its last top shell: slot s, direction d at (s*nt + d/n1,
                               // i = d%n1), nt = ceil(n/n1) (ScTopCopy)
  void ScTopCopy();
  int cln_rq = 0;               // ring, sc rays: slot swept (RingSets)
  int cln_rqm = 0;              // its sp top-mirror slot (blk mod 2)
  void ScThMirCs(DvceArray5D<Real> &a, const int i0, const int i1);
  // vet_gd_theta_wall = vacuum (1) | extrap (2) at a physical theta wall (flag 3 / 4):
  // vacuum: I = 0 entering at the wall (ghost reads 0); extrap: ghosts read the edge
  // row's intensity (zero gradient); ln chi, ln S from the edge row (ScThVacCs)
  int vgd_thvac = 0;
  void ScThVacCs(DvceArray5D<Real> &a, const int i0, const int i1);
  MeshBoundaryValuesCC *pbval_gd = nullptr;
  // exact per-shell lateral halo: vgd_i and vgd_cs carry a ghost band of vgd_w cells,
  // filled shell by shell from the 8 lateral neighbours (on-rank copies + MPI)
  int vgd_w = 0;
  // clean closure on a wedge: band of any depth keyed on the global cell
  // (ScGlobalBandInit); the sweep reads an on-rank owner (vgd_own) directly, a remote
  // one from the band, and applies the wall map at the read
  bool vgd_raw = false;
  DvceArray4D<int> vgd_own;   // (m, band k, band j, 3): owner's local index (-1 remote),
                              // its band k, j
  DvceArray5D<Real> vgd_cs;    // (m, 2, k, j, i): ln chi, ln S with the vgd_w band
  DvceArray1D<int> vgd_hloc;   // (8 nmb): local index of the slot's neighbour, -1 remote
  std::vector<int> vgd_hrank, vgd_hlid;   // (8 nmb): remote neighbour's rank, lid
  // flat message buffers, one piece per remote (block, slot), contiguous per rank pair in
  // the receiver's (block, slot) order; offsets per band depth ws in units of nv*ni:
  // (ws, 8 m + o), -1 for an on-rank slot
  DvceArray1D<Real> vgd_sbuf, vgd_rbuf;
  DvceArray2D<int> vgd_soff, vgd_roff;
  std::vector<int> vgd_prk;                  // partner ranks
  std::vector<std::vector<int>> vgd_pdsp, vgd_pscnt, vgd_prdsp, vgd_prcnt;   // (ws, p)
  int vgd_maxcnt = 0;
  bool vgd_hmpi = false;
  Real vgd_thalo = 0.0, vgd_tmpi = 0.0, vgd_ttab = 0.0;
  // halo diagnostics (rank 0; clamps all ranks)
  Real vgd_tpost = 0.0, vgd_nexch = 0.0, vgd_nbyte = 0.0, vgd_nclamp_all = 0.0;
  Real vgd_nrb = 0.0, vgd_trb = 0.0;
  MeshBoundaryValuesCC *pbval_hr = nullptr;
  void ScHaloInit();
  // vet_gd on the cubed sphere (rays/sc_cs_sweep.cpp): a block's band is its own panel's
  // index space with the four neighbour panels' strips unfolded at the edges; each band
  // position holds one global cell, filled per shell by exact keyed copies. The upwind
  // point is found in Cartesian (panel by the largest normal component, (xi, eta) by the
  // gnomonic chart); bilinear reads, linear across a seam strip, spherical barycentric
  // weights at a cube vertex.
  bool vgd_csg = false;               // vet_gd on the cubed sphere
  int vgd_csn = 0;                    // cells per panel edge
  DvceArray4D<Real> vgd_csr;          // (m, k, j, 3): band cell unit vector, 0 none
  std::vector<double> vgd_csr_h;      // same on the host (double), (m, k, j, 3)
  DvceArray2D<int> vgd_csblk;         // (m, 3): panel, global j0, global k0
  DvceArray3D<int> vgd_cscm;          // (p, q, 7): q (a,t) -> p unfolded (J,K) affine
  DvceArray2D<int> vgd_csnb;          // (p, 4): panel across side -x2, +x2, -x3, +x3
  DvceArray2D<int> vgd_csloc;         // (e, 6): on-rank copy dst (m,k,j) <- src (m,k,j)
  DvceArray2D<int> vgd_cssl, vgd_csrl;   // (e, 3): MPI send src / recv dst (m,k,j)
  std::vector<int> vgd_csprk;         // partner ranks
  std::vector<int> vgd_csps, vgd_cspn, vgd_cspr, vgd_cspm;   // send/recv displ, count
  // message buffers in pinned host memory (pack/unpack kernels access them): MPI never
  // sees a GPU pointer here
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_cssb, vgd_csrb;
  int vgd_csnloc = 0, vgd_csnsnd = 0, vgd_csnrcv = 0;
  void ScCsInit();
  void ScCsSource();
  void ScCsSweep();
  void ScCsTab(const std::vector<double> &d);   // into vgd_tabh (walls off, mr)
  template <class V>
  void ScCsHalo(V &a, const int nv, const int i0, const int i1, const bool loc,
                const int ws = -1);
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_csown;   // (m, k, j, 3) owner
  // vet_gd_halo_compact > 0 on cs: chart-aware exact lists, per shell and pass only the
  // (remote band cell, direction) values the sweep reads (ScCsCompact: a read-mask
  // kernel replays the upwind stencils; masks sent to the owners once per set). Every
  // remote read is checked against the mask: fatal on a miss.
  struct ScCsCl {
    Real alpha = -9.0;
    int64_t use = 0;
    std::vector<int64_t> so, ro;     // ((pass*n1 + s)*(np + 1) + p): list offsets
    Kokkos::View<int *, LayoutWrapper, DevMemSpace> si, ri;   // e*n + d per value
    // (e_rcv, d, s): 1 in, 2 out, 0 not read; 2 bits per value, value
    // t = (e_rcv*n + d)*n1 + s in word t/16, bits 2 (t%16)
    Kokkos::View<uint32_t *, LayoutWrapper, DevMemSpace> msk;
  };
  bool vgd_cscomp = false;
  // vet_gd_cs_list_pad: read-mask stride per receive entry (n*nx1, or rounded up to 16);
  // list-build sub-timers (rank 0): mask kernel, download + check, mark exchange, host
  // lists, upload
  int64_t vgd_csper = 0;
  double vgd_cstb_k = 0.0, vgd_cstb_d = 0.0, vgd_cstb_m = 0.0, vgd_cstb_h = 0.0,
         vgd_cstb_u = 0.0;
  // persistent host buffers of the read-mask exchange: own marks (sent), the receivers'
  // (received)
  std::vector<uint32_t> vgd_csmw_h, vgd_csmw_r;
  // the cache holds max(2, vet_rotate_avg) sets (vgd_csncl): each set's lists made once
  static constexpr int VET_VGD_CSCL_MAX = 16;
  ScCsCl vgd_cscl[VET_VGD_CSCL_MAX];
  // read masks of two entries (the set in use and the next); the others keep lists only
  Kokkos::View<uint32_t *, LayoutWrapper, DevMemSpace> vgd_csmb[2];
  int vgd_csncl = 2;
  int vgd_cscur = 0;
  int64_t vgd_csuse = 0;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_csrix;   // (m, k, j): e_rcv / -1
  Kokkos::View<int *, LayoutWrapper, DevMemSpace> vgd_csmiss;   // reads not in the mask
  double vgd_csnbuild = 0.0, vgd_cstbuild = 0.0, vgd_csvsent = 0.0, vgd_csvdense = 0.0;
  void ScCsCompact();
  void ScCsHaloC(const int pass, const int i);
  bool ScCsHaloCX(const int pass, const int i, const int eb);
  DvceArray1D<int> vgd_cslst;               // per-pass branch lists
  // vet_gd_halo_compact = 0 (off), 1 (pass branch), 2 (+ ray reach): compact messages
  int vgd_hcomp = 0;
  bool vgd_capped = false, vgd_allow_clamp = false;
  Real vgd_nclamp_seen = 0.0;
  bool vgd_hinw = true;
  int vgd_nb2 = 1, vgd_nb3 = 1;
  // several MeshBlocks along x1: the levels swept in turn
  int vgd_nb1 = 1;
  std::vector<int> vgd_x1nb, vgd_x1lv;   // (2 m + 0/1): gid below / above, -1 none; lx1
  DvceArray1D<int> vgd_x1f;   // (m): 1 innermost, 2 outermost level, 4 reflecting top,
                              // + 8 lx1
  DvceArray1D<int> vgd_x1mo;  // blocks by level (the shell list order)
  std::vector<int> vgd_shlv[2];   // (level): first code of the level in vgd_shl[pass]
  DvceArray2D<Real> vgd_x1b, vgd_x1r;    // (m, v): one shell to / from the x1 neighbour
  void ScX1Hand(const bool up);
  void ScX1Pass(const bool up);
  DvceArray4D<int> vgd_rmap;   // (m, k, j, d): mirror direction at a reflecting top
  void ScTopMirror();
  DvceArray3D<Real> vgd_isc;   // (m, k, j): flux ratio of a reflecting inner x1
  void ScInnerMirror();
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_lxy;   // (m, 0/1): lx2, lx3
  DvceArray1D<Real> vgd_csb, vgd_crb;
  // the next shell's compact mask is built while this shell's messages travel; 2 slots
  // of flags/positions (vgd_hfs/hfr) and partner bounds (vgd_hpb, pinned), tagged
  // (i, inw, ws, icut, sweep id)
  DvceArray1D<int> vgd_hfs[2], vgd_hfr[2];
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hpb[2];
  int vgd_htag[2][5] = {{-1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1} };
  int vgd_hlast = 1, vgd_hsweep = 0;
  int vgd_hnext[3] = {-1, 0, 0};
  void ScHcPrep(const int slot, const int i0, const bool inw0, const int ws);
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_pbd;   // (ws, 2 (np+1))
  std::vector<double> vgd_r1v, vgd_r1f;
  template <class V>
  bool ScHaloCompact(V &a, const int nv, const int i0, const int ws);
  // compact-halo lists per (pass, shell), built from the masks at the first exchange of a
  // mask set (direction set, cut, depth) and reused: gather pack, MPI, scatter (+ zeros
  // on the other pass's entries) instead of flag / scan / pack / expand / dense unpack
  struct ScHlEntry {
    DvceArray1D<int> sl;          // (ns): packed band address of each sent value
    DvceArray1D<int> dl, sr;      // (nr): packed band address, compact source position
    DvceArray1D<int> dw, sw;      // (nw), (3 nw): wall_interp ghosts: address, sources
    std::vector<int> sp, rp;      // (np + 1): partner bounds of the compact messages
  std::vector<int> dsto;        // (np): start of my values in each partner's receive
                                // layout (batch-made exact lists)
    Real alpha = -2.0;
    int scut = -1, ws = -1, ns = 0, nr = 0, nw = 0;
    int64_t gen = -1;             // build id (> 0) while valid
    bool nofit = false;           // over the budget for (alpha, scut, ws)
    size_t bytes = 0;
    // vet_gd_halo_exact: only the values the sweep reads (sr / sw = -1: a read ghost with
    // no source in the mask, written 0)
    bool exact = false;
  };
  static constexpr int vgd_hl_mb = 4096;   // MB per rank: budget of the halo lists
  // the shell kernel runs over the (m, k, j, d) of the pass's branch only (ScSweep)
  DvceArray1D<int> vgd_shl[2];
  int vgd_shn[2] = {0, 0};
  Real vgd_shl_alpha = -3.0;
  std::vector<ScHlEntry> vgd_hle;
  size_t vgd_hl_bytes = 0, vgd_hl_cap = 0;   // cap: bytes (set at the first build)
  int64_t vgd_hl_gen = 0;
  Real vgd_hl_nmade = 0.0, vgd_hl_nuse = 0.0, vgd_hl_ndense = 0.0;
  DvceArray1D<int> vgd_hl_fn, vgd_hl_fw;    // build scratch: dest flags -> positions
  // per band array and pass: build id of the entry whose exchange last wrote each shell's
  // band (0: never written, all zero; -1: unknown)
  std::map<const void *, std::array<std::vector<int64_t>, 2>> vgd_hl_gw;
  bool ScHlBuild(const int e, const int i, const bool inw, const int ws,
                 const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                 const std::vector<int> &pbh);
  DvceArray4D<Real> vgd_trig;   // (m, k, j, 4): sin, cos of theta, phi (shell kernel)
  // vet_gd on a stretched theta and/or the full sphere (rays/sc_sph_globalband.cpp).
  // vgd_sth: band-row theta from the node table vgd_thn (extended chart: beyond a pole
  // theta continues as (-theta, phi + pi)), upwind rows linear between nodes. vgd_pol:
  // halo slots beyond a pole take the pole partner (lx2, lx3 + nb3/2) with rows reversed
  // (vgd_hpol 1 north, 2 south); reads beyond the stored band are clamped (counted).
  bool vgd_sth = false, vgd_pol = false;
  static constexpr Real vgd_bcolat = 30.0;   // (deg) full sphere's polar band zone
  DvceArray2D<Real> vgd_thn;              // (m, c2w) band-row theta (extended chart)
  std::vector<double> vgd_thn_h;          // host copy, m-major
  DvceArray1D<int> vgd_hpol;              // (8 nmb) pole slot flag
  DvceArray2D<int> vgd_wsd;               // (2, c1) per-shell halo depth: inward, outward
  // vet_gd_band_global (full sphere): a decomposition-independent band. Every band
  // position is an exact copy of its owner's value keyed on the global cell (the cs lists
  // vgd_cs*, ScCsHalo; any depth), per shell only within its halo depth; the upwind
  // stencil is clamped to the cell's own window (k +- wu, j +- wu), not the block's band
  bool vgd_gband = false;
  DvceArray2D<Real> vgd_phc;              // (m, band k): band-column phi, global index
  std::vector<double> vgd_phc_h;          // host copy, m-major
  DvceArray1D<int> vgd_gdps, vgd_gdpr;    // depth of every send / receive entry
  std::vector<int> vgd_gdnl;              // (ws): on-rank entries with depth <= ws
  std::vector<std::vector<int>> vgd_gdsn, vgd_gdrn;   // (ws)(partner): send / recv counts
  void ScGlobalBandInit();
  // vet_gd_halo_compact > 0 with vet_gd_band_global: per shell and pass only the values
  // ScHaloKeep passes at the source cell's global (J, K); rows in the polar zone
  // (colatitude < vet_gd_band_colat) + the band depth keep their whole branch
  bool vgd_gbc = false;
  DvceArray2D<int> vgd_gbjs, vgd_gbjr;    // (entry, 2): source cell global J, K
  DvceArray1D<Real> vgd_gbth, vgd_gbdt;   // (J): node theta, theta cell width
  DvceArray1D<int> vgd_gbcap;             // (J): 1 = polar zone (branch only)
  DvceArray1D<int> vgd_gbfs, vgd_gbfr;    // flags -> positions (entry*n + v; +1)
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_gbpb;   // segment bounds
  // vgd_raw: scanned flags and bounds per (shell, pass, depth) (ScGbHaloC); _xb: the
  // other cross-fade set's
  std::map<int, std::array<DvceArray1D<int>, 2>> vgd_gbcf, vgd_gbcf_xb;
  std::map<int, std::vector<int>> vgd_gbcb, vgd_gbcb_xb;
  Real vgd_gbdmin = 0.0;
  template <class V>
  bool ScGbHaloC(V &a, const int i, const int ws);
  void ScThetaInit();                  // node table and checks
  double ScThetaDmin() const;          // narrowest lateral spacing for the band
  int ScPoleSlot(const int o, const int lx2) const;  // 0, 1 north, 2 south
  Kokkos::View<int64_t*, DevMemSpace> vgd_scan_sums;   // ScScan2 chunk sums
  // vet_gd_halo_exact: the compact-halo lists cut to the band values the sweep reads
  // (marked from its ray geometry), senders told by one request per list build
  bool vgd_hx_on = false;
  DvceArray1D<int> vgd_hx_mk;              // (m, v, k, j) dense band: read marks
  int vgd_hx_stamp = 0;                    // marks of the current build
  Kokkos::View<int, DevMemSpace> vgd_hx_oc;   // marked ghosts outside the region
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_pin, vgd_hx_pin2;
  DvceArray1D<int> vgd_hx_need, vgd_hx_nc, vgd_hx_ns;   // build scratch
  DvceArray1D<int> vgd_hx_t0, vgd_hx_t1;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_hb;   // bounds, counts
  Real vgd_hx_nmade = 0.0, vgd_hx_nout = 0.0, vgd_hx_nfall = 0.0, vgd_hx_nun = 0.0;
  Real vgd_hx_vmask = 0.0, vgd_hx_vexact = 0.0;   // received values: mask lists, exact
  std::vector<ScHlEntry> vgd_hxe;         // (2 n1): exact lists per (pass, shell)
  std::vector<ScHlEntry> vgd_hxe_alt;     // the other xfade set's lists
  bool ScHxKey(const int e, const int ws) const;   // record e matches the masks' key
  bool ScHxCut(const int e, const ScHlEntry &X) const;   // ray-cut move test
  void ScHxBuild2(const int e, const int i, const int ws);
  // exact lists of many shells in a few collectives (ScHxBatch): scratch pools,
  // per-depth partner / piece tables, pinned bounds
  DvceArray1D<char> vgd_hxb_tmp;
  int vgd_hxb_g = 0, vgd_hxb_npc = 0;
  size_t vgd_hxb_rcap = 0, vgd_hxb_scap = 0;
  DvceArray1D<int> vgd_hxb_nd, vgd_hxb_f0, vgd_hxb_f1, vgd_hxb_rq, vgd_hxb_so, vgd_hxb_oc;
  Kokkos::View<int**, LayoutWrapper, DevMemSpace> vgd_hxb_pst, vgd_hxb_pcs, vgd_hxb_pci,
                                                  vgd_hxb_pds;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hxb_hb, vgd_hxb_sp;
  Real vgd_hxb_t = 0.0, vgd_hxb_ngroup = 0.0, vgd_hxb_nrealloc = 0.0;
  void ScHxBatch(const std::vector<std::array<int, 4>> &jobs);
  // vet_gd_halo_ipc: per-shell exchange of the batch-made exact lists over node-local
  // CUDA IPC, driven by the GPU (the host never waits in the sweep)
  bool vgd_ipc_on = false;
  int vgd_ipc_state = 0;        // 0 not built, 1 on, -1 not possible (MPI path)
  Real *vgd_ipc_rbuf = nullptr;
  long long *vgd_ipc_sig = nullptr;
  int *vgd_ipc_err = nullptr;
  size_t vgd_ipc_cs = 0, vgd_ipc_half = 0;
  long long vgd_ipc_seq = 0;
  Real vgd_ipc_n = 0.0;
  DualArray1D<unsigned long long> vgd_ipc_prb, vgd_ipc_psig;
  DualArray1D<long long> vgd_ipc_phalf, vgd_ipc_pcs;
  DualArray1D<int> vgd_ipc_pslot;
  std::vector<void *> vgd_ipc_open;
  void ScIpcInit(const std::vector<int> &prk, size_t cs, size_t half);
  void ScIpcCheck();
  // pending list-path exchange of the overlap (ScHlBegin / ScHlEnd)
  struct ScHlPend {
    bool on = false, xe = false;                // xe: X is an exact (vgd_hxe) entry
    bool ipc = false;                           // the IPC path
    // IPC: begin + end as one kernel (pack, flags, wait, scatter). now: the end follows
    // the begin; done: the scatter ran; cnt: teams done packing
    bool now = false, done = false;
    Kokkos::View<int, DevMemSpace> cnt;
    long long seq = 0;
    int ex = -1, px = 0, i = -1;
    std::array<std::vector<int64_t>, 2> *gwa = nullptr;
    ScRag a;
#if MPI_PARALLEL_ENABLED
    std::vector<MPI_Request> req;
#endif
    Kokkos::Timer tq;
  };
  ScHlPend vgd_hlp;
  template <class V>
  bool ScHlBegin(V &a, const int nv, const int i0, const int ws);
  void ScHlEnd();
  std::vector<int> vgd_wsh;    // (i): band depth the shell's data needs (<= vgd_w)
  std::vector<int> vgd_wsi, vgd_wso;   // (i): same per pass (inward, outward)
  template <class V>
  void ScHalo(V &a, const int nv, const int i0, const int i1, const int ws,
              const bool mapd);
  // multi-shell (ln chi, ln S) halo in chunks of shells that fit the per-shell buffers
  void ScHaloCs(DvceArray5D<Real> &a, const int i0, const int i1);
  Real vgd_tsrc = 0.0, vgd_tswp = 0.0;
  void ScInit();
  Real RotAngle(const int cyc) const;
  bool ScRotOn() const {return vgd_rot > 0 || vgd_rot_dt > 0.0 || vgd_rot_auto;}
  int64_t RotBlock(const int cyc) const;       // rotation block
  Real RotAngleBlk(const int64_t b) const;
  double RotPhase(const int64_t b) const;   // block -> fraction of period
  Real RotAlignDt(const Real time, const Real dt);
  void ScTables(const Real alpha);
  // host part of the tables in persistent pinned buffers; the next rotation's tables
  // prefetched (host thread + upload into a second device set on vgd_tex)
  struct ScTabHost {
    Kokkos::View<Real **, LayoutWrapper, Kokkos::SharedHostPinnedSpace> dir;
    Kokkos::View<int ***, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w;
    Kokkos::View<int ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mp;
    Kokkos::View<Real ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mr;
    Kokkos::View<int *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> m3;
    Kokkos::View<Real *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w3;
    std::vector<double> mt;
    Real alpha = -9.0;
    int nwall = 0, nexact = 0, nmap = 0, nphx = 0;   // nphx: inexact phi-only wall maps
    int nthx = 0;                                     // nthx: theta-only wall maps
  };
  struct ScTabDev {
    DvceArray2D<Real> dir;
    Kokkos::View<int ***, LayoutWrapper, DevMemSpace> wall;
    Kokkos::View<int ****, LayoutWrapper, DevMemSpace> map;
    DvceArray4D<Real> mr;
    Kokkos::View<int *****, LayoutWrapper, DevMemSpace> m3;
    DvceArray5D<Real> w3;
  };
  ScTabHost vgd_tabh;
  ScTabDev vgd_tabd;
  std::future<void> vgd_tabf;
  int vgd_tab_state = 0;       // 0 idle, 1 host thread running, 3 uploaded (vgd_tabd)
  Real vgd_tab_nswap = 0.0;
  DevExeSpace vgd_tex;
  void ScTabCompute(const Real alpha, const int nthr);
  void ScTabAlloc(const bool dev);
  void ScTabPrefetch(const int cyc);
  void ScTabPoll();
  void ScBuild();
  void ScSweep();
  // shell body closure in device memory, called through a pointer (ScClosureToDevice)
  DvceArray1D<char> vgd_bbuf0, vgd_bbuf1;
#if defined(KOKKOS_ENABLE_CUDA)
  DvceArray1D<char> vgd_bbuf2, vgd_bbuf3;   // CUDA gx-false closures
#endif
  int vgd_scut = 0;            // vray_icut the running sweep uses (sync: = vray_icut)
  DevExeSpace vgd_cur;         // instance the sweep/halo use
#if MPI_PARALLEL_ENABLED
  MPI_Comm vgd_comm = MPI_COMM_WORLD;
#endif
  void ScSource();
  void ScPost();
  void VetRaySourceEnd();

  // ...in "vet_before_stagen"
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "vet_stagen"
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus StageTimers(Driver *d, int stage);
  TaskStatus Opacity(Driver *d, int stage);
  TaskStatus CalculateFluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus GasUpdate(Driver *d, int stage);
  TaskStatus RestrictU(Driver *d, int stage);
  TaskStatus SendU(Driver *d, int stage);
  TaskStatus RecvU(Driver *d, int stage);
  TaskStatus ApplyPhysicalBCs(Driver *d, int stage);
  TaskStatus Prolongate(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  // ...in "vet_after_stagen"
  TaskStatus ClearSend(Driver *d, int stage);
  TaskStatus ClearRecv(Driver *d, int stage);

 private:
  MeshBlockPack *pmy_pack;
};

} // namespace radvet
#endif // RAD_VET_RAD_VET_HPP_
