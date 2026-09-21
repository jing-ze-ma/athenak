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
