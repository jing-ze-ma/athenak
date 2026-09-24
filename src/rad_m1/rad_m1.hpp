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
// components of RadiationM1::vet_cell (closure = vet_sc)
constexpr int M1_VET_CHX = 0;   // extinction rho (kappa_F + kappa_s) [1/length]
constexpr int M1_VET_SRC = 1;   // source function in E units (4 pi S / c)
constexpr int M1_VET_J   = 2;   // E of the formal solution (4 pi J / c)
constexpr int M1_VET_K11 = 3;   // K_ab in E units: 11 22 33 12 13 23
constexpr int M1_VET_H1  = 9;   // F_a / c of the formal solution: 1 2 3
constexpr int M1_VET_CHI = 12;  // the uniaxial projection: chi, n1, n2, n3
constexpr int M1_VET_N1  = 13;
constexpr int M1_VET_D11 = 16;  // vet_tensor = full: the GUARDED D = K/J handed to the
                                // solve, 11 22 33 12 13 23, ghosts filled (periodic)
constexpr int M1_VET_GD  = 22;  // full: |D_guarded - D_raw| (max norm) of the cell
constexpr int M1_VET_NC  = 23;
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

struct VetMBState;   // rad_m1_vet.cpp: the multi-block short-characteristics sweep

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
  // ---- VET by SHORT CHARACTERISTICS (rad_m1_vet.cpp), <rad_m1>/closure = vet_sc: the
  // Eddington tensor of the multi-D implicit solve is K/J of a grey formal solution for
  // the specific intensity, computed ONCE per hydro step from the start-of-step state
  bool vet_sc;               // closure = vet_sc (default false)
  int vet_nmu, vet_nphi, vet_nray;  // mu nodes per hemisphere, azimuths, rays
  bool vet_milne;            // DIAGNOSTIC: source = exact grey Milne S(tau), gate 3
  bool vet_x1per;            // vet_x1_periodic (default false): periodic x1 sweep
  int vet_x1npass;           // ...and its number of passes (runs_3r_radwave)
  bool vet_axis_flux;        // uniaxial axis: the SC flux (true) or the principal axis
  bool vet_full;             // vet_tensor = full: the solve reads all six D_ab = K_ab/J
  Real vet_eig_min;          // vet_tensor = full: eigenvalue floor of the guarded D
  Real vet_nguard, vet_ncell;  // full: cell-calls the guard changed D / all cell-calls
  Real vet_guard_max;        // full: largest |D_guarded - D_raw| (max norm) of the run
  Real vet_d11_min;          // full: smallest guarded D_11 of the run (x1 line diagonal)
  // DIAGNOSTIC <rad_m1>/dbg_opac_patch = f (default 1 = off): all three opacities are
  // multiplied by f inside the box [x1lo,x1hi] x [x2lo,x2hi] -- a horizontally
  // inhomogeneous absorber that casts a shadow (full-tensor VET test)
  Real dbg_opac_patch, dbg_opac_x1lo, dbg_opac_x1hi, dbg_opac_x2lo, dbg_opac_x2hi;
  int vet_dump_every;        // column dump every N calls (0 = the first call only)
  std::string vet_dump;      // column dump file stem ("" = no dump)
  Real vet_time, vet_itime, vet_ncall;  // SC seconds, implicit-solve seconds, calls
  DvceArray2D<Real> vet_ang;   // (nray, 4): mu1, mu2, mu3, weight (sum of weights = 1)
  DvceArray5D<Real> vet_cell;  // (nmb, M1_VET_NC, k, j, i): see M1_VET_* below
  DvceArray5D<Real> vet_ipl;   // (nmb, 2, nray, k, j): intensity of the last two layers
  // several MeshBlocks / MPI ranks (rad_m1_vet.cpp): the halo-banded plane buffers, the
  // neighbour tables and the exchange state of the exact layer-by-layer sweep; null on
  // one MeshBlock (the single-block sweep above is then run unchanged)
  VetMBState *vet_mbs;
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
  // implicit_bc_advect: a MARSHAK x1 end also carries the advected enthalpy flux
  // A E (outflow: the end cell's E; inflow: the bath).  Default off (bitwise).
  bool impl_bc_advect;
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
  int impl_enth;                // M1_IENTH_*: face value of the enthalpy flux a E
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
  int impl_bcg_sync;            // <rad_m1>/implicit_bcg_sync: 0 = the original loop
                                // (5 blocking reductions per its), 1 = fused reductions
                                // (3, the DEFAULT), 2 = 1 plus alpha kept on the device
                                // (2 on 1 rank)
  Kokkos::View<Real, DevMemSpace> bcg_rvd;  // rhat.v on the device (sync level 2)
  // ---- GPU cost of the Krylov iteration (bench/m1_fast_0923).  Every key defaults OFF;
  // off, nothing below is allocated or referenced and the path is bitwise HEAD.
  bool impl_halo_direct;        // <rad_m1>/implicit_halo_direct: fill the ghost zones of
                                // the implicit exchanges by ONE on-rank copy kernel
  bool halo_direct_on;          // ...and the mesh allows it (same level, all neighbours
                                // on this rank on EVERY rank, no seam/pole)
  DualArray2D<int> hd_src;      // (m, 27 directions): local source MeshBlock or -1
  bool impl_odc;                // <rad_m1>/implicit_od_cache: the off-diagonal Eddington
                                // operator from a per-cell cache of sum_e d_e P_de, and
                                // ONE fused 7-point + off-diagonal kernel (bitwise)
  DvceArray5D<Real> odc;        // (m,3,k,j,i) that cache
  bool impl_stencil;            // <rad_m1>/implicit_op_stencil: the operator of each
                                // Picard pass as a 19-point stencil, built once
  DvceArray5D<Real> ost;        // (m,19,k,j,i) that stencil
  bool st_edges;                // ...and whether any edge coefficient is non-zero
  bool impl_vfold;              // <rad_m1>/implicit_vimp_fold (tests_m1/runs_4a_accel):
                                // the implicit_vimp operator part folded into the stored
                                // stencil (x2/x3 +-1 into slots 3-6, the +-2 neighbours
                                // into slots 19-24) instead of M1VimpRow per apply
  int impl_onep;                // <rad_m1>/implicit_one_pass: check period N (0 = off;
                                // default 8 for be, see ImplicitInit)
  Real impl_onep_s;             // <rad_m1>/implicit_one_pass_safety (default 3)
  Real onep_qa[3], onep_qb[3];  // last two measured contractions (be, stage 1, stage 2)
  Real onep_cnt[3];             // solves since the last measurement
  Real impl_onep_n, impl_onep_nchk;  // accepted after one pass / measurements (counters)
  bool ew_tight;                // this pass: no Eisenstat-Walker loosening
  int impl_pord;                // <rad_m1>/implicit_predictor_order: 1 or 2 (default 2
                                // for be, see ImplicitInit)
  bool impl_opsplit;            // <rad_m1>/implicit_op_split_red: the stencil apply as a
                                // plain par_for, then a separate read-only reduction
                                // kernel (instead of the fused 4-value reduction)
  bool impl_fastk;              // <rad_m1>/implicit_fast_kernels: bitwise-exact shortcuts
  bool impl_odskip;             // ...its Eddington part: D_ab = 0 off the diagonal, so
                                // the stencil build and the right-hand side skip the od
                                // terms
  bool impl_prec_float;         // <rad_m1>/implicit_precond_float: the fast path's
                                // line solves in float
  int impl_kfuse;               // <rad_m1>/implicit_krylov_fuse: 1 = the preconditioner
                                // reads/writes the Krylov vectors directly and carries
                                // the p / s updates (bitwise); 2 = 1 plus the (rhat,v)
                                // and (t,s),(t,t) reductions inside the operator kernel;
                                // 3 = 2 with two blocking reductions per iteration
  int impl_prec;                // <rad_m1>/implicit_precond: 0 = x1 line Jacobi (the
                                // original), 1 = symmetric red-black transverse line
                                // Gauss-Seidel, block-local (no communication), 2 = its
                                // forward half (red, black)
  // ---- multi-rank Krylov (tests_m1/runs_3w_krylov, rad_m1_krylov.cpp).  Both keys
  // default OFF; off, nothing below is allocated and the path is bitwise the old one.
  bool impl_kpipe;              // <rad_m1>/implicit_krylov_pipe: pipelined (Cools-
                                // Vanroose) BiCGStab, reductions hidden behind the
                                // preconditioner + operator (needs krylov_fuse = 3)
  DvceArray5D<Real> kpw;        // (m,4,k,j,i) its extra vectors phat, s, shat, rhat
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> kp_h;  // its 2 x 6 reduction slots
  bool impl_halo_mpi;           // <rad_m1>/implicit_halo_mpi: implicit_halo_direct on
                                // several ranks: on-rank copy kernel + ONE message per
                                // neighbour rank (pack / unpack kernels)
  int hm_state;                 // 0 = not built, 1 = on, -1 = not possible on this mesh
  int hm_nsend, hm_nrecv;       // entries (cells) sent / received per component
  int hm_nq;                    // components the buffers are sized for
  std::vector<int> hm_rank, hm_soff, hm_slen, hm_roff, hm_rlen;  // per neighbour rank
  DualArray1D<int> hm_sm, hm_skji, hm_sseg, hm_rm, hm_rkji, hm_rseg;  // per entry
  DualArray1D<int> hm_segs;     // (4*nseg) soff, slen, roff, rlen
  DvceArray1D<Real> hm_sbuf, hm_rbuf;
  void *hm_comm;                // MPI_Comm* (the dup'ed communicator), opaque here
  void ImplicitHaloMPIInit();
  void ImplicitHaloMPI(int nq, int c0);
  void ImplicitHaloMPIPost(int nq, int c0);    // receives, pack, on-rank copy
  void ImplicitHaloMPIFinish(int nq, int c0);  // sends, wait, unpack
  // ---- implicit_halo_overlap (tests_m1/runs_3y_halo_overlap): the Krylov halo of x
  // overlapped with the operator on the interior cells; read only when named, default
  // off (then nothing below runs).  Needs implicit_halo_mpi and implicit_op_stencil.
  bool impl_halo_ovl;
  // <rad_m1>/implicit_halo_ovl_faces (tests_m1/runs_4l_sync; read only when named,
  // default false): the shell of the overlapped operator is only as deep as the faces
  // whose ghosts arrive by MPI (hm_face, the union over the pack, set by
  // ImplicitHaloMPIInit: x1-, x1+, x2-, x2+, x3-, x3+); the other faces go with the
  // interior, whose ghosts (on-rank copy, physical: coefficient 0) are in place
  bool impl_ovl_faces;
  int hm_face[6];
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> ho_h;  // 2 x 4 partial sums
  void ImplicitHaloOp(int xc, int yc, int red, Real *out);
  void ImplicitStencilOpPart(int xc, int yc, int red, int part, int w, Real *hs);
  int ImplicitBiCGStabPipe(Real rhsmax);
  // ---- launch and host-sync cuts (tests_m1/runs_4k_launch, rad_m1_launch.cpp).  Every
  // key defaults OFF and is read only when named; off, nothing below is allocated.
  int impl_kdev;                // <rad_m1>/implicit_krylov_dev: K > 0 = BiCGStab scalars
                                // and convergence test on the device, host status read
                                // every K iterations (one rank, halo_direct, stencil)
  int impl_kdev_halo;           // <rad_m1>/implicit_krylov_dev_halo: 0 = the operator
                                // reads ghosts from the neighbour block, 1 = halo kernel
  int kdev_last[3];             // inner iterations of the previous device solve, per
                                // slot (Picard pass 0, 1, >= 2)
  int kdev_slot;                // the slot of this solve (set by ImplicitSolve)
  int kdev_x1p;                 // -1 = not known, 1 = no x1 neighbour anywhere
  Real kdev_nchk, kdev_nq;      // host status reads, queued iterations (statistics)
  DvceArray1D<Real> kdv;        // (M1_KD_SIZE) the device scalars
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> kdh;  // their host copy
  bool ImplicitKrylovDevOK();
  void ImplicitPrecondXD(int zc, int upd);
  void ImplicitOpXD(int xc, int yc, int red);
  int ImplicitBiCGStabDev(Real rhsmax);
  // <rad_m1>/implicit_op_check = K (tests_m1/gates; read only when named, default 0 =
  // off, then nothing below runs): at the first |K| solves, every operator variant of
  // the configuration applied to the same pseudo-random vector and compared
  // (rad_m1_opcheck.cpp); K > 0 fatal on a mismatch, K < 0 report only
  int impl_opchk, opchk_n;
  Real impl_opchk_tol;          // <rad_m1>/implicit_op_check_tol (default 1e-12)
  void ImplicitOpCheck();
  // ---- the Picard pass count (bench/m1_picard_0923).  Since bench/m1_defaults_0923
  // lres_test = false, conv_est = true and lin_ew_max = 1e-2 (not predictor) are the
  // DEFAULTS for closure = eddington | vet_sc | tau (the OFF settings stay the defaults
  // for m1 | minerbo | kershaw); the OFF settings reproduce the earlier path bitwise.
  int impl_plog;                // <rad_m1>/implicit_picard_log: print one line per
                                // Picard pass for the first N solves (rank 0)
  Real bcg_r0rel;               // max|b - A x0|/max|b| of the last BiCGStab call
  bool impl_lres_test;          // <rad_m1>/implicit_lres_test (default false*): require
                                // the pass-to-pass transverse change lresid < lin_tol
  bool impl_conv_est;           // <rad_m1>/implicit_conv_est (default true*): stop when
                                // the linear-rate estimate of the REMAINING change is
                                // below implicit_tol
  Real impl_ew_max;             // <rad_m1>/implicit_lin_ew_max: Eisenstat-Walker forcing
                                // of the inner tolerance (0 = off, the fixed lin_tol;
                                // default 1e-2*, or 0 when implicit_bcg_sync = 0)
  Real impl_ew_gam;             // <rad_m1>/implicit_lin_ew_gamma
  Real impl_lin_cnorm;          // <rad_m1>/implicit_lin_cnorm: inner test on the per-cell
                                // relative E error max|r_i|/(s_i E_i) (0 = off)
  Real ew_fprev, ew_etaprev;    // max|r0| and eta of the previous pass of this step
                                // (* = for eddington | vet_sc | tau; m1-type closures
                                // default to the opposite, the pre-0923 path)
  bool impl_pred;               // <rad_m1>/implicit_predictor = step (default step*):
                                // start the Picard loop from E^n + (dt/dt_prev) dE_prev
                                // (and T likewise); none = cold start
  bool pred_ok;                 // a previous-step increment is stored
  Real pred_dt;                 // the dt of that step
  DvceArray5D<Real> ipred;      // (m,3,k,j,i): dE_prev, dT_prev, T^n of this step

  // ---- <rad_m1>/time_scheme = hesdirk2 (docs/dev/rad_m1_time2_design.md sect. 5,
  // tests_m1/runs_3x_hesdirk2): two implicit stage solves inside the Heun stages.  With
  // time_scheme = be (default) nothing below is allocated or read.
  int time_scheme;              // M1_TIME_BE | M1_TIME_HESDIRK2
  bool t2_ok;                   // a valid FSAL slope K1 is stored
  int t2_afmode;                // time2_enth_vel: 0 old | 1 start (plm a_f from the
                                // stage-start a, the rest as a face mean) | 2 central
                                // (a_f central, E_f plm), in stage solves with vimp
  int t2_solve;                 // what the next ImplicitSolve does (M1_T2S_*)
  bool t2_fail;                 // the last stage solve was not admissible
  int t2_dbg_fail;              // DEBUG time2_dbg_fail: fail stage 1 at this cycle
  Real t2_nstep, t2_nbe, t2_nfall;  // stage steps, BE steps, fallbacks (counters)
  Real t2_dtprev;               // the dt of the previous step (vet_sc extrapolation)
  bool t2_vprev;                // vet_prev holds the tensor of the previous step
  bool t2_vext;                 // time2_vet_extrap (default false: D^n)
  Real t2_lin_tol;              // time2_lin_tol: implicit_lin_tol of the stage solves (<= 0:
                                // implicit_lin_tol)
  Real t2_lin_save;             // the implicit_lin_tol put back after a stage solve
  bool t2_p2s1;                // time2_stage2_pred = stage1: the stage-2 predictor is
                                // d1/2 + the history of (d2 - d1/2) (runs_5f_h2fast)
  Real t2_nclip;                // vet_sc cells whose extrapolated tensor was clipped
  DvceArray5D<Real> t2k1;       // (m,M1_T2_NK,k,j,i) the FSAL slope K1 (restart state)
  DvceArray5D<Real> t2k2;       // (m,M1_T2_NK,k,j,i) the stage-2 slope K2
  DvceArray5D<Real> t2inc;      // (m,M1_T2_NK,k,j,i) old vector - stage start state
  DvceArray4D<Real> t2f1, t2f2, t2f3;  // the face fluxes of U^n (Heun average)
  DvceArray5D<Real> ipred2;     // the predictor increment of the stage-2 solve
  bool pred2_ok;
  Real pred2_dt;
  DvceArray5D<Real> vet_prev;   // vet_sc: the start-of-step tensor of the previous step
  DvceArray5D<Real> vet_now;    // vet_sc: the start-of-step tensor of this step
  DvceArray5D<Real> vet_opac;   // vet_sc: the stage-start opacities, saved over the
                                // formal solution at U^n
  void Time2Init(ParameterInput *pin);
  bool Time2Active();           // the next step runs the stages (else BE)
  void Time2FormStage(int stage, Real dt);
  void Time2Restore(Driver *pdrive);
  void Time2VetStart();        // vet_sc: formal solution at U^n, then extrapolate
  void Time2VetExtrapolate();
  void Time2Report();
  int Time2RstNchWant();        // restart channels this run keeps (0 unless hesdirk2)
  int Time2RstNch();            // ...and writes (0 unless a slope is stored)
  void Time2RstPack(DvceArray5D<Real> &a, int nmb);
  void Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb);

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
  // <rad_m1>/implicit_vimp (tests_m1/runs_3v_vimplicit): the gas velocity of the
  // enthalpy flux implicit through the radiative force of the solve (Newton form)
  bool impl_vimp;               // the switch (default false: nothing below is touched)
  Real impl_vimp_jscale;        // DIAGNOSTIC implicit_vimp_jscale: scale of P (1)
  bool vimp_now;                // on for this step (the positivity fallback drops it)
  int iw_vimp;                  // first iw component of the M1_NIW_VIMP block, or -1
  Real vimp_nfall;              // how often the positivity fallback dropped it
  Real vimp_emin;               // the smallest E the linear solve produced with it on
  DvceArray5D<Real> vmw, vmw_c;   // exchange scratch of the M1_NVIMP_X components
  MeshBoundaryValuesCC *pbval_vm;  // ...and its exchange object
  void ImplicitVimpBuild();     // the Jacobian of a(v') E' for this Picard pass
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
                                // thread per x1 column), 1 = pcr (the DEFAULT: one team
                                // per column, parallel cyclic reduction in team scratch;
                                // not used by the gathered stack sweep, part_nblk > 1)
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
  //! STAGE S1 (rad_m1_sph.cpp): the configuration the spherical-polar wedge supports;
  //! fatal on anything else.  Called after ImplicitInit when sph_geom is set.
  void SphericalS1Check(ParameterInput *pin);
  // true on a spherical-polar mesh (a wedge clear of the poles, S1): the implicit
  // kernels then take the face areas, cell volumes and centre-to-centre distances of
  // Coordinates instead of the uniform mb_size.dx1..3 (appended as overwrites, so
  // the Cartesian arithmetic is untouched)
  bool sph_geom = false;
  // STAGE S2 (rad_m1_sph.cpp): sph_geom with a closure that is not Eddington (m1,
  // minerbo, kershaw): the radial faces take the integrating factor (M1SphDrr) and every
  // face equation the lagged curvature M1SphCurv; false keeps the S1 arithmetic
  bool sph_q = false;
  //! print the Picard statistics of the implicit solver (from the destructor, rank 0)
  void ImplicitReport();
  //! VET (rad_m1_vet.cpp): read <rad_m1>/vet_*, check the mesh, allocate
  void VetInit(ParameterInput *pin);
  //! VET: the short-characteristics formal solution and the uniaxial tensor (chi, n)
  void VetShortChar();
  //! VET, vet_tensor = full: the guarded D = K/J (and its ghosts) for the solve
  void VetFullTensor();
  //! VET: write the column of j = js (m = 0) of the last formal solution to `fname`
  void VetDump(const std::string &fname);
  //! VET: cost line at the end of the run
  void VetReport();
  //! VET, several MeshBlocks: allocate the banded buffers and neighbour tables
  void VetMBInit(ParameterInput *pin);
  //! VET, several MeshBlocks: the exact sweep with per-layer band exchange
  void VetSweepMB(bool lagged);
  //! VET, several MeshBlocks: the sweep(s) of one call (exact, or the lag diagnostic)
  void VetMBSweeps();
  //! VET: free the multi-block state (destructor)
  void VetFree();
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
  //! implicit_krylov_fuse: the pcr line solve with the right-hand side read from `rc`
  //! (upd = 0), or made on the fly as the BiCGStab p (upd = 1) / s (upd = 2) update,
  //! and the answer written to `zc`.  `col` >= 0 solves only the columns of that
  //! red-black colour (implicit_precond = rbgs); `sub` >= 0 then subtracts the
  //! transverse 5-point coupling to component `sub` from the right-hand side.
  void ImplicitPCRSolveX(int rc, int zc, int upd, Real c1, Real c2, int col, int sub);
  //! implicit_od_cache: fill odc with sum_{e!=d} d_e P_de(x) of component xc
  void ImplicitODCache(int xc);
  //! the preconditioner of the fast path: z = M^{-1} r (see implicit_precond)
  void ImplicitPrecondX(int rc, int zc, int upd, Real c1, Real c2);
  //! implicit_od_cache: y += sgn L_off(x) (with7 = false) or y = A x (with7 = true) from
  //! the cache; red > 0 also returns reductions in out[4] (see the definition)
  void ImplicitOffDiagOpC(int xc, int yc, Real sgn, bool with7, int red, Real *out);
  //! implicit_krylov_fuse = 3: BiCGStab with TWO blocking reductions per iteration
  int ImplicitBiCGStabTwo(Real rhsmax);
  //! implicit_halo_direct: the neighbour table, and the one-kernel ghost fill
  void ImplicitHaloDirectInit();
  //! implicit_op_stencil: build the 19-point operator of the pass / apply it
  void ImplicitStencilBuild();
  void ImplicitStencilOp(int xc, int yc, int red, Real *out);
  //! the fast-path operator product: stencil or od cache
  void ImplicitOpX(int xc, int yc, int red, Real *out);
  void ImplicitHaloDirect(int nq, int c0);
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
  //! right-hand side (it overwrites M1_IW_TR and M1_IW_S1..S3).  rc < 0: the caller
  //! has already staged r in M1_IW_TR.
  void ImplicitPrecond(int rc, int zc);
  //! milestone 3b phase C: solve the frozen 7-point system by BiCGStab; the answer is
  //! left in M1_IW_S2, exactly where the line-Jacobi pass leaves it.  Returns the number
  //! of inner iterations taken.
  int ImplicitBiCGStab(Real rhsmax);
  //! the same BiCGStab with fused reductions and vector updates (implicit_bcg_sync > 0)
  int ImplicitBiCGStabFused(Real rhsmax);
  void ImplicitBiCGStabEnd(int nit, bool fell_back);  // fallback / output / statistics
  void ImplicitPicardLog(int it, int nin, Real resid, Real lresid, bool srct);

  // ---- <rad_m1>/closure = tau (rad_m1_tau.cpp): the Eddington tensor of the multi-D
  // implicit solve from the COLUMN OPTICAL DEPTH measured from the top of the domain
  // (x1max): chi = exact grey plane-parallel K/J(tau) (Hopf fit), axis along -grad tau.
  // Never reads the local flux.  Rebuilt once per hydro step from the start-of-step
  // transport opacity; x1 stacks of MeshBlocks across MPI ranks are supported.
  bool tau_closure;            // closure = tau (default false)
  bool tau_ready;              // TauClosureInit has run
  int tau_nblk;                // MeshBlocks per x1 stack (uniform mesh)
  bool tau_any_mpi;            // some stack member of a local block is on another rank
  Real tau_time, tau_ncall;    // seconds in TauClosureBuild (fenced), calls
  DualArray1D<int> tau_pos;    // (m): x1 position of the block in its stack
  DualArray1D<int> tau_loc;    // (m*nblk + q): local index of stack member q, -1 remote
  std::vector<int> tau_mrank, tau_mlid;  // (m*nblk + q): rank / local id of member q
  DvceArray3D<Real> tau_sum;   // (m,k,j): the block's own column integral of rho kappa_t
  DvceArray4D<Real> tau_rcv;   // (m,q,k,j): column integrals received from member q
  HostArray3D<Real> tau_sum_h;
  HostArray4D<Real> tau_rcv_h;
  DvceArray4D<Real> tau_col;   // (m,k,j,i): column optical depth at the cell centre
  DvceArray5D<Real> tau_ten;   // (m,4,k,j,i): chi, n1, n2, n3 read by ImplicitSolve (b)
  //! closure = tau: check the mesh, build the stack topology, allocate (first call)
  void TauClosureInit();
  //! closure = tau: fill tau_ten from the transport opacity M1_IW_KT of the work array
  //! (its x2/x3 ghost layers must be current: called after the transverse halo)
  void TauClosureBuild();
  //! closure = tau: cost line at the end of the run
  void TauClosureReport();

  // ---- <rad_m1>/closure = vet_col (rad_m1_vetcol.cpp, design stage S5, option D): the
  // Eddington factor f_K = K/J of a 1-D formal solution per RADIAL COLUMN (spherical
  // impact-parameter rays on the sp wedge, Gauss rays in plane-parallel on a Cartesian
  // mesh) with the column's own extinction and source, handed to the implicit solve as a
  // FIXED uniaxial tensor about r_hat (or the M1 flux axis) through the tau closure's
  // tau_ten (tau_closure is set as well; only the tensor build differs).
  bool vet_col;                // closure = vet_col (default false)
  bool vcol_sph;               // spherical rays (sp) or plane-parallel (Cartesian)
  bool vcol_axis_flux;         // vet_col_axis = flux: n = the M1 cell flux direction
  int vcol_nc, vcol_np, vcol_nmu, vcol_every, vcol_nray;
  int vcol_dump_every;         // vet_col_dump_every (0 = off)
  bool vcol_built;             // a tensor exists (vet_col_every > 1 skips builds)
  std::string vcol_dump;       // vet_col_dump: file prefix of the column dump
  Real vcol_time, vcol_ncall, vcol_nskip;
  DvceArray2D<Real> vcol_seg;  // (ray, shell l): path length from shell l to l+1 / top
  DvceArray2D<Real> vcol_mu;   // (ray, shell): mu of the ray at the shell
  DvceArray2D<Real> vcol_w;    // (ray, shell): hemisphere quadrature weight (sum 1)
  DvceArray2D<Real> vcol_ray;  // (ray, 0..3): L (first shell), type, seg0, a_t | mu_face
  DvceArray1D<int> vcol_klast; // (shell): index of the last ray active at the shell
  DvceArray2D<Real> vcol_buf;  // (ray, column): the running intensity of each ray
  DvceArray5D<Real> vcol_mom;  // (m,5,k,j,i): J, H, K, S, chi of the solution (dump)
  std::vector<double> vcol_rc; // shell radii (host, for the dump)
  void VetColInit();           // checks, ray tables, buffers (first TauClosureInit)
  void VetColBuild();          // the formal solution -> tau_ten (chi, n)
  void VetColReport();
  void VetColDumpColumn(int ncall);

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
