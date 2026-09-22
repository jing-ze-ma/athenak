#ifndef RAD_M1_RAD_M1_HPP_
#define RAD_M1_RAD_M1_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1.hpp
//! \brief definitions for the grey photon two-moment (M1) radiation module, <rad_m1>.
//! See docs/dev/rad_m1_design.md.  MILESTONE 1a: explicit transport only -- the closure,
//! the HLL flux with thick_flux = none, the PD-ARS time integration with a null source,
//! sub-cycling and the cell-centred plumbing.  The implicit matter coupling (design
//! sect. 4) is milestone 1b and enters at the marked Coupling task in both stages.
//!
//! MILESTONE 1b adds: analytic (constant / power-law / user-hook) opacities stored as
//! rho*kappa per cell, the two thick-limit fluxes (ap_hll, scaled), and the implicit
//! local matter coupling to <hydro>, called inside BOTH PD-ARS stages.
//!
//! MILESTONE 1c adds: the advective enthalpy-flux split of the E equation for moving
//! media (<rad_m1>/advect_split, design sect. 3 "Moving fluid"), the unified form of the
//! ap_hll E-flux (<rad_m1>/ap_form, kept as an option -- the 1b pair measures better,
//! design sect. 11), and the O(v/c) source truncation used as the failing control of
//! T4/T4b (<rad_m1>/source_form).  Still absent: implicit transport (stage 3) and any
//! geometry but Cartesian (stage 4).

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations
class Driver;

namespace radm1 {

// indices of the evolved moments in u0(m,n,k,j,i)
constexpr int M1_E  = 0;
constexpr int M1_F1 = 1;
constexpr int M1_F2 = 2;
constexpr int M1_F3 = 3;
constexpr int M1_NVAR = 4;

// number of explicit stages of the PD-ARS (Chu et al. 2019) integrator used per substep.
// With S == 0 it reduces to Heun / SSP-RK2.
constexpr int M1_NSTAGE = 2;

// components of the per-cell opacity array, all stored as rho*kappa (inverse length)
constexpr int M1_OP_P = 0;   // rho kappa_P, Planck (emission) mean
constexpr int M1_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int M1_OP_T = 2;   // rho (kappa_F + kappa_s), the TRANSPORT opacity: what the
                             // flux source relaxes F with, and what tau_face is built of
constexpr int M1_NOPAC = 3;

// <rad_m1>/opacity
constexpr int M1_OPAC_CONST    = 0;
constexpr int M1_OPAC_POWERLAW = 1;
constexpr int M1_OPAC_USER     = 2;
constexpr int M1_OPAC_TABLE    = 3;   // stellar Rosseland + Planck tables (milestone 2a)

// <rad_m1>/force_reference: what the momentum coupling SUBTRACTS from the radiative
// force before handing it to the gas.  See RadiationM1::arad_ref and the bookkeeping
// note at the head of rad_m1_coupling.cpp.
constexpr int M1_FREF_NONE    = 0;
constexpr int M1_FREF_WB_ARAD = 1;

//----------------------------------------------------------------------------------------
//! \struct M1OpacTab
//! \brief the two stellar opacity tables, as device Views held BY VALUE on the module so
//! that the whole struct can be captured by value in a device lambda (a host pointer
//! cannot be dereferenced there).  Both tables share ONE (log10 T, log10 rho) grid -- the
//! format src/pgen/box_convection.cpp's ReadOpacityTable produces -- and both are read
//! with diffusion/conduction.hpp's RosselandTable, unchanged.
//!
//! The lookup is in KELVIN and g/cm^3: `tunit` converts the EOS's CODE temperature to
//! kelvin (T[K] = T_code*tunit, the same direction Conduction uses with
//! Units::temperature_cgs()), `dunit` the code density to g/cm^3, and `kunit` the
//! tabulated cm^2/g back to code units.

struct M1OpacTab {
  DvceArray2D<Real> kr;        // log10 kappa_Rosseland [cm^2/g], (iT, iD)
  DvceArray2D<Real> kp;        // log10 kappa_Planck    [cm^2/g], same grid
  DvceArray1D<Real> lT, lD;    // log10 T[K], log10 rho[g/cm^3]
  int nT = 0, nD = 0;
  Real tunit = 1.0;
  Real dunit = 1.0;
  Real kunit = 1.0;
};

// <rad_m1>/reconstruct, as a plain int for the device (same order as the code-wide
// ReconstructionMethod enum, so the two can be compared).  ppm4 and above read a
// 5-cell stencil per face state and need <mesh>/nghost >= 3.
constexpr int M1_RECON_DC    = 0;
constexpr int M1_RECON_PLM   = 1;
constexpr int M1_RECON_PPM4  = 2;
constexpr int M1_RECON_PPMX  = 3;
constexpr int M1_RECON_WENOZ = 4;

// entries of the implicit-solve diagnostic counter (see RadiationM1::cnt)
constexpr int M1_CNT_NSOLVE = 0;   // cells handed to the bracketed Newton
constexpr int M1_CNT_ITSUM  = 1;   // total iterations over those cells
constexpr int M1_CNT_ITMAX  = 2;   // largest iteration count of any one cell
constexpr int M1_CNT_NFAIL  = 3;   // cells that hit maxit without reaching the tolerance
constexpr int M1_CNT_NBRAK  = 4;   // cells for which no bracket [T_lo, T_hi] was found
constexpr int M1_NCNT = 5;

// bracketed-Newton controls (design sect. 4a)
constexpr int  M1_MAXIT = 100;
constexpr Real M1_RTOL  = 1.0e-12;

//----------------------------------------------------------------------------------------
//! \struct RadiationM1TaskIDs
//! \brief container to hold TaskIDs of all rad_m1 tasks

struct RadiationM1TaskIDs {
  TaskID irecv;
  TaskID copyu;
  TaskID closure;
  TaskID flux;
  TaskID sendf;
  TaskID recvf;
  TaskID opac;
  TaskID update;
  TaskID coupl;
  TaskID c2p;
  TaskID restu;
  TaskID sendu;
  TaskID recvu;
  TaskID bcs;
  TaskID prol;
  TaskID csend;
  TaskID crecv;
};

//----------------------------------------------------------------------------------------
//! \class RadiationM1

class RadiationM1 {
 public:
  RadiationM1(MeshBlockPack *ppack, ParameterInput *pin);
  ~RadiationM1();

  // parameters read from <rad_m1>
  Real c_light;        // speed of light in code units (required)
  Real chat_over_c;    // reduced speed of light factor (default 1)
  Real chat;           // = chat_over_c * c_light
  Real cfl_rad;        // CFL number for the radiation substep (default 0.4)
  Real e_floor;        // floor on E
  bool subcycle;       // sub-cycle the module inside the mesh timestep
  int chi_kind;        // M1_CHI_*: levermore (closure = m1) | minerbo | kershaw
  bool eddington;      // closure = eddington (chi = 1/3); default is the M1 closure
  std::string thick_flux_str;   // none | ap_hll | scaled
  int thick_flux;               // M1_THICK_*
  Real scaled_pref;             // thick_flux = scaled: tau_c = scaled_pref*tau_face
  std::string ap_form_str;      // unified | alpha2
  int ap_form;                  // M1_APFORM_*
  bool advect_split;            // the advective enthalpy-flux split of the E equation
  bool split_vel_recon;         // <rad_m1>/split_vel = recon: the velocity that builds
                                // the enthalpy flux A = v E + v.P is RECONSTRUCTED to
                                // the face with the same method as (E, f_i).  = cell
                                // reverts to the 1c-A behaviour (each side's own cell
                                // velocity), which is first order at the face.
  bool source_ovc;              // TRUNCATE the source to O(v/c) -- the control of T4/T4b,
                                // NOT a production option: it drops the beta^2 E and
                                // beta.P.beta pieces of E0 and uses F instead of F0 in
                                // the beta.g work term (Skinner & Ostriker type form)
  // MILESTONE 1d, design sect. 13.  <rad_m1>/f_source = cell | wb.  `cell` (the 1c
  // behaviour) divides the centred 2 dx radiation-pressure gradient by the cell's OWN
  // rho*kappa_F, which is wrong by O(1) in the two cells straddling an opacity jump;
  // `wb` uses the Bloch et al. (2021) eq. 18 trapezoidal interface source instead, i.e.
  // an effective opacity (1/2)[(rho kappa)_{i-1/2} + (rho kappa)_{i+1/2}] per direction
  // with the same arithmetic face mean the thick-limit flux uses.  The two agree exactly
  // wherever rho*kappa is linear across the cell, so only a jump can see the difference.
  bool f_source_wb;

  // opacities, per unit mass, in code units (design sect. 4).  kappa = kappa0 for
  // opacity = const, kappa0 (rho/rho_ref)^opac_a (T/t_ref)^opac_b for powerlaw.
  int opacity_type;             // M1_OPAC_*
  Real kappa_p, kappa_e, kappa_f, kappa_s;
  Real opac_rho_ref, opac_t_ref, opac_a, opac_b;
  bool opac_zero;               // all four opacities are identically zero
  M1OpacTab otab;               // opacity = table: filled by SetOpacityTables

  // <rad_m1>/force_reference (milestone 2a).  With M1_FREF_WB_ARAD the momentum handed
  // to the gas is the RESIDUAL between the instantaneous radiative force and the
  // reference acceleration arad_ref(z) that an external well-balanced scheme already
  // carries (box_convection's Phi_eff).  arad_ref is a per-cell device array filled by
  // the problem generator through SetForceReference.
  int force_ref;
  DvceArray4D<Real> arad_ref;   // (m,k,j,i), a reference x1 acceleration, code units

  // RSLA start-up check (design sect. 2): v_max*tau_max/chat, evaluated once, on the
  // first filled opacity array.  rsla_vmax <= 0 means "measure max |v| over the domain".
  Real rsla_warn, rsla_vmax;
  bool rsla_force, rsla_done;

  // matter coupling (design sect. 4)
  bool coupling;       // run the implicit local solve at all
  bool gas_feedback;   // write the solve's reaction back into hydro u0.  With this off
                       // the medium is a PRESCRIBED static background: the radiation
                       // still feels the opacity, the gas does not move or heat.  That
                       // is what the static-medium tests T3/T3b/T3c mean by "no hydro".
  Real arad;           // radiation constant in CODE units: the equilibrium energy
                       // density is arad*T^4 with T the EOS's code temperature

  // ---- DEBUG SWITCHES.  All default to the production value, so that leaving them out
  // of an input file is bitwise inert; they exist to take the coupling apart when a
  // configuration is unstable and it has to be decided WHICH piece drives it.
  bool opac_freeze;    // <rad_m1>/opac_freeze: fill the opacity array ONCE, at the first
                       // call, and never again.  Kills the kappa-mechanism feedback
                       // loop (a compression that raises kappa raises the trapped
                       // energy) while leaving transport and coupling otherwise intact.
  bool opac_frozen;    // internal: the first fill has happened
  bool dbg_gas_force;  // <rad_m1>/dbg_gas_force = false: the radiation field still
                       // relaxes its flux, but the gas receives NO momentum and no work
  bool dbg_gas_heat;   // <rad_m1>/dbg_gas_heat = false: skip the energy exchange (a)
                       // entirely, so gas and radiation exchange no energy
  int dbg_tensor;            // DIAGNOSTIC (VET scaffolding, 3i): 0 none | 1 frozen |
                       // 2 tilt | 3 tau: where the Eddington tensor of the multi-D
                       // implicit solve comes from (see ImplicitSolve step (b))
  Real dbg_tensor_tilt;      // amplitude [rad] of the prescribed tilt, dbg_tensor = tilt
  bool dbg_tensor_init;      // frozen / tilt: the stored tensor exists
  Real dbg_trans_memory;     // DIAGNOSTIC (3b phase G): weight of the F^n memory term
                       // of the transverse face flux; 1 = backward Euler (bitwise)
  bool dbg_gas_force_trans;  // 3b phase E, transport = implicit only:
                       // <rad_m1>/dbg_gas_force_trans = false keeps the x1 radiative
                       // force (and with it the well-balanced reference) but hands the
                       // gas NO x2/x3 momentum and no transverse work -- which is the
                       // one piece of the multi-D coupling that has no hydrostatic
                       // reference to be measured against

  // true when this module's dtnew must be folded into Mesh::NewTimeStep, i.e. when
  // sub-cycling is off OR when no other module sets the mesh timestep
  bool sets_mesh_dt;

  // ---- MILESTONE 3a: IMPLICIT x1 TRANSPORT (docs/dev/rad_m1_implicit_design.md).
  // Everything below is inert with the default <rad_m1>/transport = explicit; see
  // rad_m1_implicit.hpp / rad_m1_implicit.cpp.
  int transport;                // M1_TRANSPORT_*
  int nstage;                   // stages per substep: M1_NSTAGE explicit, 1 implicit
  Real impl_cfl;                // <rad_m1>/implicit_cfl: dtnew = impl_cfl*dx/c (<=0 off)
  Real impl_tol;                // Picard tolerance on max(|dE|/E, |dT|/T)
  int impl_maxit;               // maximum Picard iterations per solve
  bool impl_opac_update;        // re-evaluate the opacities inside the Picard loop
  bool impl_allow_multid;       // run a 2-D/3-D set of INDEPENDENT x1 columns
  Real marshak_q;               // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min, ibc_x1max;     // M1_IBC_*
  Real iflux_x1min, iflux_x1max;  // the imposed face flux of M1_IBC_FLUX
  Real iebath_x1min, iebath_x1max;  // M1_IBC_MARSHAK: the INCIDENT bath
                                // energy density, F_f = +-c q (E - E_bath);
                                // 0 is the plain free surface
  DvceArray4D<Real> f0x1;       // face-normal comoving flux on x1 faces, (m,k,j,i);
                                // PERSISTENT state, carried by the restart file
  DvceArray4D<Real> f0x1n;      // its start-of-step copy (not restarted)
  DvceArray5D<Real> iw;         // per-cell work array of the solve, M1_NIW components
  Real impl_nstep, impl_itsum, impl_itmax, impl_nfail;   // host-side Picard counters

  // ---- MILESTONE 3a2 (the limits of 3a).  See rad_m1_implicit.hpp for what each
  // constant means and docs/dev/rad_m1_implicit_design.md sect. 7 for the measurements.
  int impl_flux;                // M1_IFLUX_*: the spatial form of the implicit E-flux
  int impl_recon;               // M1_IRECON_*: dc, or plm by deferred correction
  Real impl_recon_w;            // implicit_recon = plm_dc: the weight the deferred
                                // correction carries.  <= 0 (the default) means the
                                // automatic 1/(1 + chat dt/dx), which is what keeps the
                                // deferred-correction fixed-point iteration a contraction
                                // at every CFL; 1.0 is the pure plm flux, which diverges
                                // above CFL ~ 1.
  int impl_part;                // M1_IPART_*: how a multi-block column is solved
  // ---- MILESTONE 3c: implicit_flux = blend, the smooth per-face convex blend of the
  // central and the berthon face flux.  All inert unless implicit_flux = blend.
  int impl_blend;               // M1_IBLEND_*: which lagged weight w_f is used
  int impl_blend_fmode;         // M1_IBFM_*: max or mean of the two cells' reduced flux
  int impl_blend_mode;          // M1_IBMODE_*: blend the whole flux, or add the HLL
                                // dissipation alone on top of the full central flux
  Real impl_blend_tau0;         // w = exp(-(tau_face/tau0)^2)
  Real impl_blend_flo;          // smoothstep lower edge in the reduced flux
  Real impl_blend_fhi;          // smoothstep upper edge in the reduced flux
  DvceArray5D<Real> ifw;        // per-x1-FACE work array, M1_NIFW components
  // the partitioned (gathered) line solve, LIMIT 4.  Every rank that owns a piece of a
  // column sends its (a,b,c,r) rows to the column's ROOT rank, which runs the identical
  // serial Thomas sweep and sends the solution back.
  // MILESTONE 3b, LIMIT 4: every rank that owns a piece of an x1 column sends the
  // assembled rows (a,b,c,r) of that column to the ROOT block (the one with the lowest
  // x1 logical location), which runs the IDENTICAL serial Thomas sweep on the whole
  // column and sends the solution back.  The arithmetic order is literally unchanged, so
  // 1, 2 and 4 MeshBlocks along x1 agree to the last bit by construction.
  int part_nblk;                // MeshBlocks along x1 of one column (1 = no partition)
  int part_nroot;               // x1 stacks whose root block is on THIS rank
  int part_nx1g;                // part_nblk*nx1, rows of one gathered column
  int part_nlay;                // ghost layers the x1 halo of the solve exchanges
  int part_nqa;                 // quantities of halo A (see ImplicitX1Halo)
  DualArray1D<int> part_pos;    // (nmb) position of block m in its x1 stack
  DualArray1D<int> part_slot;   // (nmb) gather slot of m's root, or -1 if it is remote
  DualArray1D<int> part_nbr;    // (2*nmb) local index of m's x1 neighbours (lo,hi), or -1
  DvceArray5D<Real> part_sys;   // (nroot,6,nk,nj,nx1g) gathered rows + Thomas scratch
  DvceArray2D<Real> part_sbuf;  // (nmb, 4*nk*nj*nx1) MPI staging of a member block
  DvceArray2D<Real> part_rbuf;  // (nroot*nblk, 4*nk*nj*nx1) MPI staging of a root
  HostArray2D<Real> part_sbuf_h, part_rbuf_h;
  DvceArray3D<Real> part_hbuf;  // (nmb,4,nq*nlay*nk*nj) x1 halo: send lo/hi, recv lo/hi
  HostArray3D<Real> part_hbuf_h;
  std::vector<int> part_mrank;  // (nroot*nblk) rank of each member of a rooted stack
  std::vector<int> part_mgid;   // (nroot*nblk) global id of each member
  std::vector<int> part_rootgid;  // (nmb) global id of block m's root
  std::vector<int> part_rootrank;  // (nmb) rank of block m's root
  std::vector<int> part_nbrrank;  // (2*nmb) rank of m's x1 neighbours, or -1
  std::vector<int> part_nbrgid;   // (2*nmb) global id of m's x1 neighbours, or -1
  bool part_any_mpi;            // true if any stack or x1 halo crosses a rank
  // the imposed-flux boundary hand-off, LIMIT 3.  See ImplicitSolve.
  bool impl_recon_freeze;       // implicit_recon_lag = step: evaluate the deferred
                                // correction once per step, not once per Picard pass
  int impl_recon_npass;         // implicit_recon_npass: FREEZE the plm deferred
                                // correction (and with it the limiter's choice) after
                                // this many Picard passes, to break the limit cycle the
                                // limiter otherwise drives.  <= 0 = never freeze (3a2)
  Real impl_res_floor;          // scale of the Picard convergence test: 0 = the pure
                                // relative change of 3a, x > 0 = |dE|/max(E, x*max(E))
  bool impl_bmom_half;          // give a physical boundary face HALF its flux to the one
                                // interior cell (the cell-averaged share) instead of all

  // ---- MILESTONE 3b phase B: transport = implicit, the TRANSVERSE (x2/x3) couplings.
  // All inert with transport = explicit | implicit_x1.
  int impl_solver;              // M1_ISOLV_*: line_jacobi | bicgstab
  Real impl_lin_tol;            // <rad_m1>/implicit_lin_tol: the max-norm residual of the
                                // FULL 7-point linear system, relative to the right-hand
                                // side, that the line-Jacobi outer iteration must reach.
                                // It is tested SEPARATELY from the Picard residual: a
                                // lagged transverse coupling can stall and look
                                // converged.
  Real impl_linsum, impl_linmax;  // statistics of the final linear residual per solve
  bool trans_on;                // transport == implicit AND the mesh is multi-D
  bool trans_x3;                // ...and three-dimensional
  DvceArray4D<Real> f0x2, f0x3;   // face-normal comoving fluxes on the x2 / x3 faces,
                                // (m,k,j,i).  PERSISTENT state, carried by the restart
                                // file next to f0x1; allocated only when trans_on.
  DvceArray4D<Real> f0x2n, f0x3n;  // their start-of-step copies (not restarted)
  DvceArray5D<Real> thw;        // (m,M1_NHALO_T,k,j,i) scratch the transverse halo
                                // exchanges through the module's ordinary CC boundary
                                // machinery
  DvceArray5D<Real> thw_c;      // its (unused) coarse buffer; SMR/AMR is a fatal
  MeshBoundaryValuesCC *pbval_th;  // the exchange object of thw
  DvceArray5D<Real> thq;        // (m,M1_NHALO_Q,k,j,i) the NARROW transverse halo: the
                                // components a Picard pass moves once the closure is
                                // frozen (implicit_closure_lag = step)
  DvceArray5D<Real> thq_c;      // its (unused) coarse buffer
  MeshBoundaryValuesCC *pbval_tq;  // the exchange object of thq
  bool halo_shell;              // the deep interior of the scratch halo arrays may be
                                // skipped (same-level neighbours only, no cubed sphere,
                                // no polar boundary)

  // ---- MILESTONE 3b phase C: implicit_solver = bicgstab.  All null under line_jacobi.
  bool bicg_on;                 // impl_solver == M1_ISOLV_BICGSTAB and trans_on
  int impl_lin_maxit;           // <rad_m1>/implicit_lin_maxit, the BiCGStab iteration cap
  DvceArray5D<Real> krw;        // (m,1,k,j,i) scratch: ONE Krylov vector, exchanged with
                                // the same cell-centred machinery thw uses
  DvceArray5D<Real> krw_c;      // its (unused) coarse buffer
  MeshBoundaryValuesCC *pbval_kr;  // the exchange object of krw
  Real bcg_nsolve, bcg_itsum, bcg_itmax;  // inner-iteration statistics
  Real bcg_nbreak, bcg_nfall, bcg_nred;   // breakdowns, line-Jacobi fallbacks, reductions

  // ---- MILESTONE 3b phase D: the OFF-DIAGONAL Eddington terms and the closure lag.
  // All inert with transport = explicit | implicit_x1 and on a 1-D mesh.
  int impl_offdiag;             // M1_OD_*: lagged | operator | none
  int od_now;                   // the mode this STEP is running: impl_offdiag, or
                                // M1_OD_NONE after an operator step produced a
                                // non-positive E (the positivity fallback)
  Real od_nfall;                // how often that fallback fired
  Real od_emin;                 // the smallest E the linear solve produced in the run
  Real impl_crelax;             // <rad_m1>/implicit_closure_relax: the weight w of
                                // chi^{k+1} = (1-w) chi^k + w chi(f^{k+1}) between Picard
                                // passes (1 = no relaxation), and the same for n
  bool impl_crelax_thin;        // relax only where the cell is optically THIN
                                // (theta = 1/(1 + c dt rho kappa_t) > 1/2)
  bool impl_clag_step;          // <rad_m1>/implicit_closure_lag = step: freeze chi and n
                                // at the START-of-step state for the whole step, so each
                                // step is one linear solve plus the T nonlinearity

  // ---- the TRANSVERSE realizability limiter, <rad_m1>/implicit_trans_limit.
  // Default `none` = not allocated and not referenced, so every earlier configuration
  // is bitwise unchanged.  See ImplicitTransTheta.
  int impl_tlim;                // M1_TLIM_*: none | lp
  Real impl_tfmax;              // <rad_m1>/implicit_trans_fmax, the reduced-flux cap
  DvceArray4D<Real> thx2, thx3;  // theta_f = 1/(1 + c dt (kt_f + klim_f)) on the x2/x3
                                // faces, the ONE number every use of the transverse
                                // theta reads
  DvceArray4D<Real> klx2, klx3;  // the lagged limiter opacity klim_f itself, which
                                // follows implicit_closure_lag

  // ---- MILESTONE 3g: the GAS-RADIATION energy coupling of the implicit solve.
  // Both options default to false; nothing below is then allocated or referenced, so
  // every earlier configuration is bitwise unchanged.  See ImplicitSolve steps (c)/(f).
  bool impl_gas_newton;         // <rad_m1>/implicit_gas_newton: update T' from the
                                // ELIMINATED (Schur) relation instead of re-running the
                                // bracketed root find in every Picard pass
  bool impl_eos_cache;          // <rad_m1>/implicit_eos_cache: serve e(T), c_v from a
                                // per-cell 1-D Hermite in ln T built once per step
  int impl_ecnt;                // <rad_m1>/implicit_eos_cache_nt, the half-width of the
                                // cached window in TABLE temperature cells
  bool impl_eccheck;            // <rad_m1>/implicit_eos_cache_check: one TRUE-table
                                // evaluation per cell at the end of the step, which
                                // measures the cache error and corrects T'
  int iw_gas;                   // first iw component of the M1_NIW_GAS block (see
                                // rad_m1_implicit.hpp); < 0 when neither option is on
  int impl_nec;                 // components of `ecache`
  DvceArray5D<Real> ecache;     // (m,nec,k,j,i) the frozen-density e(T) cache
  Real newt_nfb;                // Newton fallbacks to the bracketed root find, whole run
  Real gas_ncell;               // cell-passes of the gas solve, whole run (the scale the
                                // two counters above and below are read against)
  Real ec_nmiss;                // EOS-cache misses, whole run
  Real ec_emax;                 // max |e_cache - e_table|/e_table at the end of a step
  Real ec_tmax;                 // max relative error of the exchanged energy q =
                                // SRCR - SRCB E', measured against the true table

  // ---- MILESTONE 3e: ANDERSON ACCELERATION of the Picard map, <rad_m1>/implicit_accel.
  // Default `none` = nothing below is allocated and nothing is called, so every earlier
  // configuration is bitwise unchanged.  See ImplicitAccelApply.
  int impl_accel;               // M1_IACC_*: none | anderson
  int impl_line_solver;         // <rad_m1>/implicit_line_solver: 0 = thomas (one
                                // thread per x1 column), 1 = pcr (one team per column,
                                // parallel cyclic reduction in team scratch)
  int impl_pcr_team;            // <rad_m1>/implicit_pcr_team: team size of the pcr
                                // solve on a GPU (0 = next power of 2 >= nx1, max 256)
  bool impl_pcr_check;          // <rad_m1>/implicit_pcr_check: run BOTH line solvers
                                // every call and record max|dx|/max|x| (diagnostic)
  Real pcr_chk_max, pcr_chk_n;  // the check's running max and call count
  DvceArray4D<Real> pcr_chk;    // the other solver's answer (allocated on first use)
  int impl_and_m;               // <rad_m1>/implicit_anderson_m, the history depth
  Real impl_and_beta;           // <rad_m1>/implicit_anderson_beta, the mixing parameter
  int impl_and_start;           // <rad_m1>/implicit_anderson_start, the first (0-based)
                                // Picard pass that is accelerated; earlier passes are
                                // plain Picard and only feed the history
  int aa_nc;                    // components of the fixed-point vector: 4 in multi-D
                                // (E,F1,F2,F3), 2 on an x1-only solve (E,F1)
  DvceArray5D<Real> aa_xc;      // (m,nc,k,j,i) x_k, the SCALED state entering the pass
  DvceArray5D<Real> aa_fc;      // (m,nc,k,j,i) g_k = G(x_k) - x_k, scaled
  DvceArray5D<Real> aa_xp, aa_fp;  // the previous pass' pair, for the differences
  DvceArray6D<Real> aa_dx;      // (m,q,nc,k,j,i) history of dX_q = x_{q+1} - x_q
  DvceArray6D<Real> aa_df;      // ...and of dG_q = g_{q+1} - g_q; a RING of impl_and_m
  DvceArray4D<Real> aa_sc;      // (m,k,j,i) the per-cell scale max(E^n, e_floor), fixed
                                // over the step, that makes E and F/c commensurate
  int aa_nh;                    // history columns in use
  int aa_head;                  // ring head (the column the next difference overwrites)
  bool aa_hasp;                 // aa_xp/aa_fp hold a usable previous pass
  Real aa_fnp;                  // ||g|| of the previous pass (< 0: none)
  Real aa_nacc, aa_nrst;        // accelerated passes and history restarts, whole run

  // evolved variables
  DvceArray5D<Real> u0;         // (E, F1, F2, F3)
  DvceArray5D<Real> u1;         // state at the start of the substep
  DvceArray5D<Real> coarse_u0;  // on 2x coarser grid (for SMR/AMR)
  DvceFaceFld5D<Real> uflx;     // fluxes on cell faces
  DvceArray5D<Real> opac;       // rho*kappa per cell, M1_NOPAC components
  DvceArray5D<Real> ugas1;      // hydro (IEN, IM1, IM2, IM3) at the start of the substep

  // implicit-solve diagnostics, accumulated on the device and printed at the end of the
  // run (M1_CNT_*).  Reals rather than ints so that the sums cannot wrap on a long run.
  DualArray1D<Real> cnt;

  // boundary communication
  MeshBoundaryValuesCC *pbval_u;

  // timestep bookkeeping
  Real dtnew;     // cfl_rad*min(dx)/chat, over this pack
  Real dt_sub;    // the dt of the substep currently being taken
  int nsub;       // number of substeps in the current mesh step

  // PD-ARS explicit-stage weights, u0 = gam0*u0 + gam1*u1 - beta*dt_sub*div(F)
  Real gam0[M1_NSTAGE], gam1[M1_NSTAGE], beta[M1_NSTAGE];

  // reconstruction method
  ReconstructionMethod recon_method;
  int recon_code;               // M1_RECON_*, the device-side copy of recon_method
  std::string recon_str;

  RadiationM1TaskIDs id;

  // functions
  void AssembleRadM1Tasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  //! set dt_sub and nsub for a mesh step of length dt_mesh; returns nsub
  int SetSubsteps(Real dt_mesh);
  //! print the implicit-solve counters (called from the destructor, rank 0)
  void ReportCounters();
  //! hand the module its two opacity tables (opacity = table); called by the pgen
  void SetOpacityTables(const DvceArray2D<Real> &kr, const DvceArray2D<Real> &kp,
                        const DvceArray1D<Real> &lT, const DvceArray1D<Real> &lD,
                        const int nT, const int nD);
  //! the per-cell reference acceleration of force_reference = wb_arad; the caller owns
  //! the array and must keep it alive for the run
  void SetForceReference(const DvceArray4D<Real> &a);
  //! design sect. 2: evaluate and print v_max*tau_max/chat.  Runs once, from the first
  //! Opacity task that has a filled array (tau needs rho*kappa, which the ctor has not).
  void RSLACheck();
  // ---- milestone 3a (rad_m1_implicit.cpp)
  //! read the implicit-solver parameters and allocate its arrays; a no-op in
  //! transport = explicit
  void ImplicitInit(ParameterInput *pin);
  //! print the Picard statistics of the implicit solver (from the destructor, rank 0)
  void ImplicitReport();
  //! let a problem generator name the x1 boundary types of the implicit solve
  void SetImplicitX1BC(int lo_type, Real lo_flux, int hi_type, Real hi_flux);
  //! the whole backward-Euler step, in place of the explicit stage chain
  TaskStatus ImplicitSolve(Driver *d, int stage);
  //! milestone 3b LIMIT 4: build the x1 stack topology of the gathered line solve
  void ImplicitPartitionInit();
  //! exchange the x1 ghost layers the multi-block solve needs; `eponly` picks the
  //! one-quantity set (E of the new iterate) instead of the six lagged ones
  void ImplicitX1Halo(bool eponly);
  //! gather the assembled rows onto the roots, run Thomas there, scatter back
  void ImplicitGatherSolve();
  //! milestone 3b phase B: exchange the LAGGED quantities the x2/x3 faces need with
  //! ALL six neighbours (periodic wrap and MPI included), through the module's ordinary
  //! cell-centred boundary machinery on the scratch array thw.  `nq` is how many of the
  //! M1HaloCompT list go: M1_NHALO_T (all of them), M1_NHALO_Q (the components a Picard
  //! pass can move once the closure is frozen) or 1 (E of the new iterate alone)
  void ImplicitTransverseHalo(int nq);
  //! the common core of the two halo exchanges: pack `nq` components into the scratch
  //! array that matches that width, run the ordinary CC exchange, unpack.  `c0 >= 0`
  //! carries ONE named component of iw (the Krylov vector) instead of the M1HaloCompT
  //! list
  void ImplicitHaloExchange(int nq, int c0);
  //! copy the halo components between iw and the scratch array `sc` over the HALO SHELL
  //! (`topack` picks the direction); see the definition for why the deep interior is
  //! neither sent nor received
  void ImplicitHaloCopy(DvceArray5D<Real> &sc, int nq, int c0, bool topack);
  //! milestone 3b phase B: the lagged transverse (x2/x3) divergence of one Picard pass,
  //! split into its diagonal part (M1_IW_TDIA) and its right-hand side (M1_IW_TRHS)
  void ImplicitTransverseTerms(bool first);
  //! the TRANSVERSE realizability limiter: fill thx2/thx3 (and, when newk, the lagged
  //! limiter opacity klx2/klx3) for the current pass.  Inert under
  //! implicit_trans_limit = none, where thx2/thx3 are not even allocated.
  void ImplicitTransTheta(bool newk);
  //! milestone 3e: snapshot the SCALED state entering a Picard pass (x_k) into aa_xc.
  //! Called at the top of every pass when implicit_accel = anderson.
  void ImplicitAccelSave();
  //! milestone 3e: at the end of Picard pass `it` the iw state is G(x_k).  Form the
  //! fixed-point residual g_k, push (dX, dG) onto the history, solve the small least-
  //! squares problem and overwrite the iw state with the accelerated iterate x_{k+1}
  //! (realizability re-applied).  The statistics go into aa_nacc / aa_nrst.
  void ImplicitAccelApply(int it);
  //! milestone 3b phase C: the x1 line solve of the assembled rows (M1_IW_TA..TR ->
  //! M1_IW_S2), Thomas or cyclic Thomas, gathered over the stack when part_nblk > 1.
  //! It IS the preconditioner of the BiCGStab wrapper and the line-Jacobi pass itself.
  void ImplicitTridiagSolve();
  //! implicit_line_solver = pcr: the same line solve by parallel cyclic reduction,
  //! one Kokkos team per x1 column (GPU)
  void ImplicitPCRSolve();
  //! the Thomas / cyclic-Thomas sweep (implicit_line_solver = thomas)
  void ImplicitThomasSolve();
  //! milestone 3b phase C: exchange ONE component of iw with all six neighbours
  void ImplicitKrylovHalo(int comp);
  //! milestone 3b phase C: y = A x with the frozen 7-point operator (one halo of x)
  void ImplicitApplyOp(int xc, int yc);
  //! milestone 3b phase D: accumulate sgn * L_off(x) into the component yc, with
  //! L_off the cell-row contribution of the OFF-DIAGONAL Eddington terms of every face
  //! equation, evaluated at the component xc with the closure of this Picard pass
  //! FROZEN (so it is linear in x).  No communication: the caller has already put x in
  //! the ghost zones (ImplicitApplyOp) or x IS the iterate, whose halo is current.
  void ImplicitOffDiagOp(int xc, int yc, Real sgn);
  //! milestone 3b phase C: z = M^{-1} r, the x1 line solve applied to an arbitrary
  //! right-hand side (it overwrites M1_IW_TR and M1_IW_S1..S3)
  void ImplicitPrecond(int rc, int zc);
  //! milestone 3b phase C: solve the frozen 7-point system by BiCGStab; the answer is
  //! left in M1_IW_S2, exactly where the line-Jacobi pass leaves it.  Returns the number
  //! of inner iterations taken.
  int ImplicitBiCGStab(Real rhsmax);

  // ...in "m1_before_stagen"
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "m1_stagen"
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus ApplyClosureLimits(Driver *d, int stage);
  TaskStatus Opacity(Driver *d, int stage);
  TaskStatus CalculateFluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus Update(Driver *d, int stage);
  TaskStatus Coupling(Driver *d, int stage);
  TaskStatus HydroConToPrim(Driver *d, int stage);
  TaskStatus RestrictU(Driver *d, int stage);
  TaskStatus SendU(Driver *d, int stage);
  TaskStatus RecvU(Driver *d, int stage);
  TaskStatus ApplyPhysicalBCs(Driver *d, int stage);
  TaskStatus Prolongate(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  // ...in "m1_after_stagen"
  TaskStatus ClearSend(Driver *d, int stage);
  TaskStatus ClearRecv(Driver *d, int stage);

 private:
  MeshBlockPack *pmy_pack;
};

} // namespace radm1
#endif // RAD_M1_RAD_M1_HPP_
