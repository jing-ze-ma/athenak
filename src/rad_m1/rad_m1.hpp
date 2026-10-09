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
#include <thread>  // NOLINT(build/c++11)
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
constexpr int M1_VET_HP1 = 23;  // hrup-1009: sum w I max(mu_a, 0), a = 1 2 3 (half-range)
constexpr int M1_VET_NC  = 26;
// tau_ten slots of the LATERAL correction (vet_col_lat; the cs CS2 interface): dD_ab in
// the mesh basis, order rr, r-a, r-b, aa, ab, bb (sp: a = theta, b = phi)
constexpr int M1_TT_LAT0 = 4;
constexpr int M1_TT_NLAT = 6;
constexpr int M1_VGD_NMAX = 432;   // vet_gd: at most 12 x 6^2 directions
constexpr int M1_OP_P = 0;  // rho kappa_P, Planck (emission) mean
constexpr int M1_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int M1_OP_T = 2;   // rho (kappa_F + kappa_s), the TRANSPORT opacity: what the
                             // flux source relaxes F with, and what tau_face is built of
constexpr int M1_NOPAC = 3;
// vet_scatter (vet-scatter-1006): a 4th component, allocated ONLY when the key is on
// (opac then has M1_NOPAC + 1 components): rho kappa_e, the electron-scattering part of
// M1_OP_T.  M1_OP_T stays the TOTAL extinction rho (kappa_a + kappa_e) = rho kappa_R.
constexpr int M1_OP_S = 3;

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
  // vet_scatter (vet-scatter-1006): the split of the tabulated TOTAL kappa_R into true
  // absorption kappa_a = kappa_R - kappa_e and electron scattering kappa_e (Ma+2026
  // sect. 2.2.4).  scat 0 = no split (os = 0, bitwise the old lookup), 1 = constant
  // kappa_e = kes, 2 = kappa_e from the EOS electron fraction, kxe x_e/mu (M1ScatterEos,
  // rad_m1_opacity.hpp).  kes, kfl, kxe in CODE opacity per unit mass (times kunit).
  int scat = 0;
  Real kes = 0.0;    // constant kappa_e, vet_kappa_es (0.34 cm^2/g = 0.2(1+X), X = 0.7)
  Real kfl = 0.0;    // Ma's floor: kappa_e <= kappa_R - kfl, vet_kappa_floor (1e-5)
  Real kxe = 0.0;    // sigma_T/m_u = 0.40062 cm^2/g: kappa_e = kxe n_e/(n_tot mu)
};

//----------------------------------------------------------------------------------------
//! \fn VetScatterMix
//! \brief <rad_m1>/vet_scatter (vet-scatter-1006, read only when named, default off =
//! bitwise): the grey source S of the VET formal solutions (vet_sc, vet_col,
//! vet_col_lat, vet_gd) in E units, for a cell of TOTAL extinction chi = rho kappa_R
//! (M1_OP_T / M1_IW_KT, unchanged by the split), Planck absorption rkp = rho kappa_P
//! and electron scattering rkes = rho kappa_e (opac M1_OP_S, so chi - rkes = rho kappa_a
//! is the true absorption of the split), from the thermal source b and the mean
//! intensity j (both in E units):
//!   form 1 (`ma`, Jiang 2021 ApJS 253 49 eq. 6 as used by Ma+2026 sect. 2.2.4):
//!     dI/ds = rho k_s (J - I) + rho k_a (B - I) + rho (k_P - k_a)(B - J),
//!     k_a + k_s = k_R  ->  S = J + (rho k_P/chi) (B - J), NOT capped (the split of k_R
//!     into k_a and k_e drops out exactly); clipped at S >= 0.
//!   form 2 (`absorption`): eq. 6 without the k_dP term, emission by true absorption:
//!     S = J + (rho k_a/chi)(B - J) = [(chi - rkes) B + rkes J]/chi.
//! b, j: vet_col_source = gas: b = a T^4, j = E (start of step / Picard iterate);
//! relaxed: b = a T*^4 after the local backward-Euler exchange and j = the relaxed E*
//! (vet_scatter_j = relaxed, default: the consistent end-of-step pair: S - E* =
//! (E* - E)/(c dt chi), which is ~0 in a stiff cell -- every form then collapses to
//! S ~ E*), or j = the CURRENT E (vet_scatter_j = current; with form ma a startup FATAL:
//! kappa_P/kappa_R amplifies the B*-E^n lag, Picard NON-CONV on the BSG).
//! The default source eps B + (1 - eps) E with eps = min(rho k_P/chi, 1) is NOT computed
//! here: the callers keep their own expression untouched when the key is off.
KOKKOS_INLINE_FUNCTION
Real VetScatterMix(const int form, const Real rkp, const Real chi, const Real rkes,
                   const Real b, const Real e) {
  if (form == 1) {
    const Real ep = rkp/chi;
    return fmax(e + ep*(b - e), 0.0);
  }
  return ((chi - rkes)*b + rkes*e)/chi;
}

//! \fn M1ScatterGuard
//! \brief the guarded kappa_e of the split (Ma+2026 sect. 2.2.4), per unit mass: kt is
//! the TOTAL kappa_R, ke the electron-scattering opacity; where kt < ke + kfl the
//! scattering part is kt - kfl (>= 0), so that kappa_a = kt - kappa_e >= kfl
KOKKOS_INLINE_FUNCTION
Real M1ScatterGuard(const Real kt, const Real ke, const Real kfl) {
  return (kt < ke + kfl) ? fmax(kt - kfl, 0.0) : ke;
}

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
  // OPTIONAL explicit energy source of the radiation (erg/cm^3/s, code units), (m,k,j,i).
  // Unallocated and esrc_on = false by default (every path bitwise).  A problem generator
  // that allocates and fills it sets esrc_on = true (transport = implicit only);
  // ImplicitSolve then adds dt_solve*(chat/c)*esrc to the OLD vector of E and the same
  // rate to the hesdirk2 slope, so a constant esrc raises E by exactly dt*esrc per step
  // under be and hesdirk2 alike (stage retries redo the step from U^n with be).
  // he_star_m1: the frozen MLT flux.
  DvceArray4D<Real> esrc;
  bool esrc_on = false;
  // rad-beam-1008 <rad_m1>/vet_source_noesrc (default false = bitwise): the vet_col /
  // vet_col_lat / vet_gd source reads E^n WITHOUT the explicit deposit dt (chat/c) esrc
  // that the implicit solve's old vector E^n carries (the frozen MLT scaffold): the
  // deposit is transported away inside the step by the solve, but a LOCAL relaxed
  // source keeps it in the cell (AG Car A: S up to 10 E at 0.985-0.994 R_ph, L_fs 3 L)
  // fixbundle-1009 F1: every reader (incl. the vet_col inner-boundary intensity) takes
  // the physical E^n = u0 directly (be only; dt-free, right for nsub > 1)
  bool vsrc_noes = false;

  // m1-mhd (docs/dev/m1_mhd_0927.md): the fluid is <hydro> or <mhd> (FluidRef,
  // m1_fluid.hpp).  Under MHD the conserved energy carries |B|^2/2, which every kernel
  // that forms the gas internal energy subtracts: emag0 = |B|^2/2 of the CURRENT field
  // (b0), emag1 = of the U^n field (b1), both built from the faces with the C2P's
  // average (EmagBuild).  Hydro: 1x1x1x1 dummies, never read (the `if (fl_mhd)`
  // branches are dead), so hydro is bitwise.
  bool fl_on = false;           // a fluid exists (<hydro> or <mhd>)
  bool fl_mhd = false;          // it is <mhd>
  bool emag_hold = false;       // Opacity must not rebuild emag0 (the U^n swap)
  DvceArray4D<Real> emag0, emag1;
  void EmagBuild(bool un);      // un = false: emag0 from b0; true: emag1 from b1
  // after a host-side write of the fluid's u0 (coupling write-back, restores): restrict,
  // exchange (u only), physical BCs, prolongate, ConToPrim.  Hydro: the hydro task calls
  // verbatim.  MHD: the u-only exchange through pbval_u -- MHD::InitRecv/ClearRecv(-1)
  // would post and then wait on B receives that nobody sends here (deadlock on > 1 rank).
  void FluidRefresh(Driver *pdrive);

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
  bool vet_bc_bath;          // a MARSHAK x1 end gives the SC rays its bath (m1-h2div)
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
  int dbg_opac_part;            // DIAGNOSTIC: 0 all, 1 only flux (KT), 2 only P/E
  bool impl_opac_newton;        // implicit_opac_newton: d(rho kappa_T)/dT in the x1 rows
  // implicit_opac_newton_guard (m1-opn-guard, default 0.5; <= 0 = off): a face's Newton
  // term is dropped from a row whose diagonal it would take below guard x its value
  // without it; opn_nskip_d counts the dropped (row, face) pairs on the device,
  // opn_nskip the all-rank total (summed at the report)
  Real impl_opn_guard = 0.0;
  int impl_opn_guard_mode = 2;   // implicit_opac_newton_guard_mode (bits 1 rhs, 2 off)
  // opacity cliffs (opacnewt-cliff-1007, both read only when named, 0 = off = bitwise):
  // with s = d ln kappa_T/d ln T of the one-sided Newton difference, a cell with |s| >
  // implicit_opac_newton_slope_off drops its Newton term (ktd = 0, the face stays
  // Picard); implicit_opac_newton_slope_max damps it, ktd x min(1, slope_max/|s|)
  Real impl_opn_soff = 0.0, impl_opn_smax = 0.0;
  DvceArray1D<Real> opn_nskip_d;
  Real opn_nskip = 0.0;
  // m1-positivity (every key read only when named; absent = off = bitwise the old code):
  //  implicit_g0_limit = w > 0: on the FIRST Picard pass the lagged net-absorption term
  //    g0 of the face-flux equations is clipped to |chat dt g0| <= w (max(E^k, E_old) +
  //    (chat/c) max(e_gas,0)) (the energy a cell can exchange in one step)
  //  implicit_g0_exchange: from the second pass on g0 = -(SRCR - SRCB E^k)/(chat dt), the
  //    exchange of the linearised source row (= the pointwise g0 at the fixed point)
  //  implicit_pos_gas: energy-conserving positivity limiter of the written-back gas eint
  //  implicit_pos_floor: the E floor takes the added energy from the gas (conserving)
  // pos_cnt_d: device counters (M1_POS_*), pos_cnt the all-rank totals at the report
  Real impl_g0_lim = 0.0;
  bool impl_g0_exch = false;    // implicit_g0_exchange (m1-positivity, see ImplicitSolve)
  bool impl_pos_gas = false;
  bool impl_pos_floor = false;
  Real impl_pos_gas_frac = 1.0e-3;
  DvceArray1D<Real> pos_cnt_d;
  Real pos_cnt[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // fixbundle-1009 F3 <rad_m1>/implicit_pos_floor_solve (default false = bitwise): the
  // T-solve / accept no longer raise the iterate to e_floor (EP = max(S2, e_floor) made
  // energy that implicit_pos_floor never saw, and the face fluxes rebuilt from the raised
  // EP no longer matched the solve's): EP = max(S2, 1e-6 e_floor), and the write-back
  // floor raises E to e_floor with the energy taken from the gas (FLR_DE) or counted
  // (FLR_UN); what the 1e-6 guard itself raised, (c/ch)(1 + SRCB)(EP - S2), joins it
  bool impl_pos_floor_s2 = false;
  DvceArray4D<Real> ktd;        // d(rho kappa_T)/dT at the iterate (implicit_opac_newton)
  bool impl_allow_multid;       // run a 2-D/3-D set of INDEPENDENT x1 columns
  Real marshak_q;               // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min, ibc_x1max;     // M1_IBC_*
  Real iflux_x1min, iflux_x1max;  // the imposed face flux of M1_IBC_FLUX
  Real iebath_x1min, iebath_x1max;  // M1_IBC_MARSHAK: the INCIDENT bath
  // implicit_marshak_face = linear (m1-sp-order2, spherical-polar wedge only): the
  // Marshak face flux takes the face E extrapolated from r^2 E of the two end cells
  // (M1SphMarshakCoef in the row, limiter remainder deferred) instead of the end cell's
  // E (the default `cell`, first order).  Read only when named.  The outer face keeps
  // the cell form under vet_col_surface_q, whose q = H(face)/J(top cell) is built for
  // the cell E.
  bool impl_mface_lin = false;
  // implicit_face_weight = distance (default on sp) | equal: how the x1 face values of
  // the transport opacity rho kappa_T, the velocity v1 and G0 are formed from the two
  // adjacent cells in the x1 face equation, the E rows, the vimp rows and the momentum
  // deposit.  equal: 0.5 (q_L + q_R).  distance (spherical-polar wedge only; a Cartesian
  // grid is uniform): linear interpolation to the face from the two cell midpoints,
  // q_f = (dr_R q_L + dr_L q_R)/(dr_L + dr_R), dr the radial cell widths; identical to
  // equal (bitwise) where the two widths agree to 1e-12.  Default distance on the wedge,
  // equal elsewhere; a restart whose input predates the key keeps equal (rad_m1_implicit
  // .cpp, where it is read).
  bool impl_face_wdist = false;
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
  // implicit_res_dmin / implicit_res_rmax (read only when named; 0 = off): cells with
  // rho < dmin or r > rmax are left out of the Picard convergence norm (resid, lresid);
  // they are still solved every pass.  res_exmax = the largest excluded residual on the
  // last pass of any step, res_nex = steps whose excluded residual was >= implicit_tol
  // (steps that would not have stopped there with the full norm).
  Real impl_res_dmin = 0.0, impl_res_rmax = 0.0;
  bool impl_real_couple = false;      // implicit_realisable_coupling (gd_physfix_1005)
  // implicit_thin_freeze (m1-picard-aa; default 1e-2 on fresh runs, 0 = off): cells with
  // c dt rho kappa_P < this at the start of the solve keep their start-of-step opacities
  // for every pass and take the frozen-opacity gas-T find
  Real impl_thin_frz = 0.0;
  DvceArray4D<Real> thin_frz;   // 1 = frozen cell of this solve (implicit_thin_freeze)
  // fixbundle-1009 F2 <rad_m1>/implicit_face_opac_n (default false = bitwise): the
  // step-frozen x1 face data of the blend (alpha, tau_f weight) from rho kappa_T at T^n
  // (snapshot ktn, taken before the pass-0 opacity update), not at the predicted T
  bool impl_face_ktn = false;
  DvceArray4D<Real> ktn;
  Real impl_stall_fac = 0.5;    // implicit_stall_fac (implicit_gas_newton_switch test)
  int impl_gn_sw = 0, impl_gn_sw_min = 20;   // implicit_gas_newton_switch (window, min)
  Real impl_gn_nsw = 0.0;       // solves switched off the gas Newton update
  // implicit_resid_fatal (read only when named; 0 = off): stop the run when a Picard
  // solve ends non-converged with resid above this value or non-finite (a DIVERGED
  // solve, not the usual 1e-7 stall), before its state is written back and spreads
  Real impl_res_fatal = 0.0;
  // implicit_resid_fatal_masked (resmask_1003; read only when the res mask is on, default
  // 1.0 then, 0 = off): stop the run when the largest residual of the EXCLUDED cells at
  // the last pass of a solve exceeds this value or is non-finite (a masked cell that
  // diverges is otherwise invisible to the stopping test and to implicit_resid_fatal)
  Real impl_res_fatal_masked = 0.0;
  bool impl_res_mask = false;
  Real impl_res_exmax = 0.0, impl_res_nex = 0.0;
  // implicit_timers = N (m1-fast4; default 0 = off): fenced host timers of the M1 stage
  // tasks and of the parts of ImplicitSolve, accumulated over cycles >= N and printed
  // by ImplicitReport.  The fences change the timing a little; only a diagnostic.
  int tmr_c0 = 0;
  bool tmr_on = false;
  double tmr_acc[12] = {}, tmr_cnt[8] = {}, tmr_last = 0.0;
  Kokkos::Timer tmr_t;
  void TmrMark(int c);
  // implicit_mr_every = k (m1-fast4, rad_m1_mr.cpp; default 1 = off): MULTI-RATE
  // radiation.  The hydro takes every step alone; every k-th step (the first after k/2)
  // the implicit radiation + coupling is advanced over the whole window Delta by an
  // operator-split, L-stable, stiffly accurate SDIRK2 (two stage solves), placed at the
  // window centres so that the composition is Strang (second order).
  int mr_every = 1;
  bool mr_on = false;          // an MR radiation step is being taken (Delta = mr_dt)
  bool mr_first = true;        // the first window is half long (Strang)
  int mr_cnt = 0;              // hydro steps since the last radiation step
  Real mr_acc = 0.0;           // time since the last radiation step
  Real mr_dt = 0.0;            // Delta of the current / last radiation step
  Real mr_nr = 0.0, mr_nfall = 0.0;   // radiation steps, BE fallbacks (counters)
  int mr_kc = 0, mr_kn = 0;    // window length (steps) of the last / next window
  Real mr_theta = 0.0;         // implicit_mr_theta: the adaptive guard (0 = fixed k)
  bool mr_peq = true;          // implicit_mr_peq: stage-A local-equilibrium start
  Real mr_thmax = 0.0, mr_ksum = 0.0;   // largest theta seen, sum of window lengths
  // vet_sc_every = N (m1-fast4, rad_m1_mr.cpp; default 1): the vet_sc formal solution
  // on every N-th stage-1 (or multi-rate) tensor build only; in between the tensor is
  // extrapolated linearly in time from the last two formal solutions (O(dt^2))
  int vsc_every = 1;
  int vsc_cnt = 0, vsc_nb = 0;  // builds requested, formal solutions stored (0-2)
  Real vsc_t0 = 0.0, vsc_t1 = 0.0, vsc_tn = 0.0;   // their times, this one
  Real vsc_nskip = 0.0, vsc_nclip = 0.0;          // extrapolated builds, clipped cells
  DvceArray5D<Real> vsc_d0, vsc_d1;               // (m,M1_T2_NVET,k,j,i)
  bool VscSkip();
  void VscStore();
  // coupling_split (ke-dt-0926, rad_m1_mr.cpp; read only when named, default = off):
  // how the implicit radiation + coupling R is placed relative to the Heun hydro H.
  //   0 default  hesdirk2 in the Heun stages (or implicit_mr_every)
  //   1 strang   H(dt) then R(dt) by the multi-rate SDIRK2 every step (k = 1)
  //   2 alternate  H R on even steps, R H on odd steps
  //   3 mix      coupling_mix_h2 hesdirk2 steps, then coupling_mix_s strang steps
  // The step counter travels in the multi-rate restart header (mr_cnt).
  int csplit = 0;
  int csplit_nh2 = 1, csplit_ns = 2;
  bool dbg_hydro_off = false;   // DIAGNOSTIC (ke-dt-0926): the Driver skips hydro stages
  // DIAGNOSTIC <rad_m1>/dbg_cell_lo, dbg_cell_hi (read only when named; default off):
  // the implicit write-back prints the per-cell energy budget of active x1 cells
  // lo..hi (offsets from is) of the first (k, j) column of MeshBlock 0
  int dbg_cell_lo = -1, dbg_cell_hi = -1;
  // force_reference_work (ke-dt-0926; read only when named, force_reference = wb_arad
  // only): `full` (default) = the solve hands the gas the work of the FULL force and the
  // WB source none (rad_m1_coupling.cpp header; its internal-energy cancellation is
  // O(dt^2) per step); `split` = the WB source gives the gas the work of rho arad_ref at
  // its own stage (the pgen reads fref_wsplit), the solve gives only the residual's
  // work, and the radiation pays both
  bool fref_wsplit = false;
  bool fref_wsplit_ok = false;  // set by a pgen whose WB source gives the reference work
                                // (BEFORE SetForceReference, which resolves `auto`)
  bool fref_wsplit_auto = false;  // force_reference_work = auto, not yet resolved
  // force_reference_work = split with time_scheme = be (fref-split-cons-1006): the WB
  // source gives the gas rho a_ref v at each hydro stage, and the Heun combination keeps
  // sum_s c_s beta_s dt (rho a_ref v)_s of it; the be solve, which runs once after the
  // whole step, paid v' dt rho a_ref at the post-kick velocity instead, an O(dt^2) per
  // step energy leak (-1.8e-3 L_top at dt 88 s in the BSG).  The pgen now also adds its
  // increment to fref_wacc (acc = gam0_s acc + increment, the same linear combination
  // the gas energy undergoes), and the be write-back makes E pay exactly fref_wacc
  // (times dt_solve/dt_mesh).  <rad_m1>/force_reference_work_pay = post (read only when
  // named) restores the old post-kick payment.  hesdirk2 and multi-rate keep theirs.
  DvceArray4D<Real> fref_wacc;
  bool fref_pay_post = false;
  int fref_hnst = 0;
  Real fref_hgam0[4] = {0.0, 0.0, 0.0, 0.0}, fref_hbeta[4] = {0.0, 0.0, 0.0, 0.0};
  bool FrefWaccOn() const;
  // gam0 of the hydro stage whose beta dt is bdt (FATAL if none matches)
  Real FrefWaccGam0(Real bdt, Real dt_mesh) const;
  // the accumulator, (re)allocated to the shape of arad_ref
  DvceArray4D<Real> FrefWacc();
  int mr_nsub = 1;              // implicit_mr_nsub: R(Delta) as nsub steps of Delta/nsub
  int mr_tab = 0;               // implicit_mr_tab: 0 sdirk2, 1 trbdf2 (DIAGNOSTIC)
  bool mr_k0ok = false;         // t2k1 holds f(Y_0) of this R (the last R's final slope)
  int CSplitStep();             // per step: 0 = default path, 1 = H then R, 2 = R then H
  bool F4RstHdr() const {return (mr_every > 1) || (vsc_every > 1) || (csplit > 0);}
  bool MRActive() const {return mr_every > 1;}
  void MRInit(ParameterInput *pin);
  void MRStep(Driver *pdrive, Real dt);
  void MRSolve(Driver *pdrive, Real dlt);
  void MRReport();

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
  Real impl_blend_xthin = 0.0;    // hrup-1009: transparency override X0 (0 off)
  Real impl_blend_r0 = 1.5;       // hrup-1009: implicit_blend = knudsen, R0
  Real impl_blend_alpha = 1.0;   // hrup-1009: implicit_blend = idort, x = alpha tau_f
  Real impl_blend_tau0;         // w = exp(-(tau_face/tau0)^2)
  Real impl_blend_flo;          // smoothstep lower edge in the reduced flux
  Real impl_blend_fhi;          // smoothstep upper edge in the reduced flux
  DvceArray5D<Real> ifw;        // per-x1-FACE work array, M1_NIFW components
  // blendall-1009: <rad_m1>/implicit_flux_faces = all.  berthon | blend on EVERY face of
  // the multi-D solve: Cartesian x1/x2/x3 and the sp LATERAL x2/x3 faces (the sp radial
  // faces keep the sp-blend-1008 path).  ifw2/ifw3 hold M1_IFW_AL/HCL/HCR of the x2/x3
  // faces (DG unused); the central part of a transverse face keeps 1 - AL through
  // thx2/thx3.  Off (the default): nothing allocated, every kernel bitwise unchanged.
  bool impl_flux_all = false;   // the key
  bool blat_on = false;         // transverse faces carry the berthon/blend coefficients
  bool bvec_x1 = false;         // Cartesian x1 faces take the multi-D beam vector
  bool impl_beam_fs = false;    // implicit_flux_beam = fs: |f| = |H|/J of vet_sc
  // hrup-1009: implicit_flux_beam = halfrange.  Every face flux of the all-faces blend
  // is (1 - w) F_central + w c (r+_L E_L + r-_R E_R), r+- = H+-/J the HALF-RANGE ratios
  // along the face normal (vet_sc / vet_gd rays; otherwise the isotropic + beam model of
  // the lagged closure), w = exp(-(tau_f/tau0)^2) (implicit_blend = tau) or 1 (berthon)
  bool impl_beam_hr = false;
  bool impl_hr_model = false;   // implicit_hr_model (no ray data: closure model)
  bool hr_q = false;            // the per-column outer q from h+_1 (vcol_q)
  DvceArray5D<Real> vgd_hr, vgd_hrt;   // (m, 6, k, j, i): r+_r r-_r r+_t r-_t r+_p r-_p
  DvceArray1D<Real> vgd_hrm;           // shell mean of one twin component
  void VetGdHalfRange(const int stage);
  DvceArray5D<Real> ifw2, ifw3;
  void ImplicitLatFaceCoef();
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
  // rad-beam-1008: implicit_recon_dgpass (default false): with the face coefficients
  // frozen for the step (implicit_recon_lag = step) the plm_dc correction is re-made on
  // every Picard pass from the iterate's E (plm in E only, the frozen face f), instead of
  // once from E^n (an explicit anti-diffusion that is unstable: pulse amplitude 2.7-3.5)
  bool impl_recon_dgpass = false;
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
  Real bcg_nkeep;               // fallbacks that kept the Krylov iterate (fallback best)
  int impl_bcg_maxrst;          // <rad_m1>/implicit_bcg_max_restarts (default 2)
  bool impl_bcg_keep;           // <rad_m1>/implicit_bcg_fallback = best
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
  bool impl_onep_auto;          // <rad_m1>/implicit_one_pass_auto (default true): switch
                                // one_pass off per kind of solve when it does not pay
  int impl_onep_awin;           // ..._auto_window: check periods without an acceptance
                                // before it is switched off (default 2)
  int impl_onep_arep;           // ..._auto_reprobe: switched-off solves before a re-probe
                                // of one check period (default 64)
  Real onep_off[3];             // 1: one_pass switched off for this kind (be, s1, s2)
  Real onep_actr[3];            // on: eligible solves since the last one-pass acceptance;
                                // off: switched-off solves since it was switched off
  Real onep_prb[3];             // 1: on as a re-probe, not yet confirmed by an acceptance
  Real impl_onep_ndis, impl_onep_nprb, impl_onep_nren, impl_onep_nfpr;  // switch-offs /
                                // re-probes / confirmed / failed (counters, this run)
  bool ew_tight;                // this pass: no Eisenstat-Walker loosening
  int impl_pord;                // <rad_m1>/implicit_predictor_order: 1 or 2 (default 2
                                // for be, see ImplicitInit)
  bool impl_opteam;             // <rad_m1>/implicit_op_team_red (m1-fast3): the fused
                                // stencil + dots as one-cell-per-thread teams of 256
                                // with a team reduction into opt_part, then a small
                                // reduction over the teams (round-off vs the fused one)
  DvceArray2D<Real> opt_part;   // (teams, 4) partial sums of implicit_op_team_red
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
                                // forward half (red, black), 3 = mg (rbgs_fwd + block-
                                // local semicoarsened levels, rad_m1_precond.cpp), 4 =
                                // mg_gc (mg + the global band coarse space)
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
  // ---- <rad_m1>/implicit_halo_ipc (m1-perf-0928; default false): the off-rank part of
  // implicit_halo_mpi between ranks on the SAME node without MPI and without a host
  // wait.  The pack kernel writes straight into the neighbour's receive buffer (CUDA
  // IPC, double-buffered by the parity of the exchange count), a 1-thread kernel sets
  // the neighbour's arrival flag; the receiver's 1-thread kernel spins on its own flag,
  // unpacks, and acknowledges.  Same numbers as the MPI path (bitwise).
  bool impl_halo_ipc;
  int hi_state;                 // 0 = not built, 1 = on, -1 = not possible (MPI path)
  long long hi_seq;             // exchanges so far (identical on every rank)
  long long hi_nchk;            // exchanges since the last check of the error flag
  Real *hi_rbuf;                // own receive buffer, 2 parities x hi_half Reals
  long long *hi_sig;            // own flags: [0,32) arrivals, [32,64) acks, by segment
  int *hi_err;                  // device error flag (a spin timed out)
  size_t hi_half;
  DualArray1D<unsigned long long> hi_prb, hi_psig;  // peer rbuf / sig, per segment
  DualArray1D<long long> hi_phalf;  // peer parity-half size, per segment
  DualArray1D<int> hi_proff;    // peer's receive offset (entries) for my data
  DualArray1D<int> hi_pslot;    // my segment index in the peer's list
  void ImplicitHaloIPCInit();
  void ImplicitHaloIPCPost(int nq, int c0);
  void ImplicitHaloIPCFinish(int nq, int c0);
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
  void ImplicitDetOpRed(int xc, int yc, int red, Real *out);   // implicit_det_reduce
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
  // implicit_precond = mg (rad_m1_precond.cpp): block-local (x2,x3) semicoarsening
  int mg_nlev;                  // <rad_m1>/implicit_mg_levels: levels incl. the fine one
  bool mg_halo;                 // <rad_m1>/implicit_mg_halo: fine residual with the halo
  std::vector<DvceArray5D<Real>> mgc;   // level l >= 1: (m,9,nk_l,nj_l,nx1)
  std::vector<int> mg_nj, mg_nk;        // per level (level 0 = the MeshBlock)
  void ImplicitMGBuild();
  void ImplicitMGApply(int rc, int zc, int upd, Real c1, Real c2);
  // implicit_precond = mg_gc (rad_m1_precond.cpp, tests_m1/runs_5p_coarse2): a GLOBAL
  // coarse space of per-x1-layer (x2,x3)-band means in front of mg; used only under it
  bool gc_on = false;
  int gc_b2 = 1, gc_b3 = 1;     // <rad_m1>/implicit_gc_bands2 / _bands3
  int gc_nb = 1, gc_n1 = 0;     // bands, global x1 cells
  int gc_n2 = 0, gc_n3 = 0;     // global x2 / x3 cells
  int gc_nlb2 = 1, gc_nlb = 1;  // band slots per block: in x2, in all
  int gc_kl = 0;                // half bandwidth of the coarse matrix
  int gc_nrc = 1;               // row chunks of the band reductions
  bool gc_per2 = false, gc_per3 = false;   // periodic x2 / x3 mesh boundaries
  bool gc_ok = false;           // the coarse factorisation of this pass is usable
  DvceArray5D<Real> gcw;        // (m,5,k,j,i) nbands = 1: row sums of A per x1 offset
  DvceArray1D<Real> gc_pd;      // (m,k,lb2,q,i) partial band sums over x2 (device)
  DvceArray1D<int> gc_off;      // (m*3) global fine offsets x1, x2, x3 of each block
  DvceArray1D<int> gc_binfo;    // (m*nlb*5) band id (or -1), j0, j1, k0, k1 (local)
  std::vector<int> gc_binfo_h, gc_off_h, gc_piv;
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_part;  // (m,lb,q,c,i) band sums
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_x;     // (gi*nb + b) x_g
  std::vector<Real> gc_lu, gc_g;  // banded LU of the coarse matrix, right-hand side
  void ImplicitGCInit();
  void ImplicitGCSum(int nq);
  void ImplicitGCBuild();
  void ImplicitGCPre(int rc, int upd, Real c1, Real c2);
  void ImplicitGCAdd(int zc);
  // implicit_precond = mg_gf (rad_m1_precond.cpp, tests_m1/runs_5t_fast5box): a GLOBAL
  // coarse space of the lowest (x2,x3) Fourier modes of each x1 layer in front of mg
  bool gf_on = false;
  int gf_k = 2;                 // <rad_m1>/implicit_gf_modes: |k2|, |k3| <= K
  int gf_nb2 = 1, gf_nb3 = 1, gf_nm = 1;   // 1-D basis sizes, modes per layer
  int gf_nkc = 1, gf_kcw = 1;   // k chunks of the layer reductions, their width
  bool gf_ok = false;           // the coarse space is set up (per-mode flags: gf_okd)
  DvceArray1D<Real> gf_u, gf_v;     // 1-D bases (a, gj) / (c, gk)
  DvceArray1D<Real> gf_cud, gf_cvd; // their correlation sums (a, d + 2), d = -2..2
  DvceArray1D<Real> gf_nrm;         // 1 / |phi_q|^2
  DvceArray1D<Real> gf_part;        // (m, kc, q, i) partial layer sums
  DvceArray1D<Real> gf_g;           // (gi*25 + s) layer sums of the stencil slots
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gf_gh;   // (gi*nm + q) multi-rank
  DvceArray1D<Real> gf_hg, gf_hx;   // ((c*n2+gj)*n1+gi): sum_a u_a(gj) (gn, x)(a, c)
  DvceArray1D<Real> gf_gn, gf_x;    // (gi*nm + q): P^T r / |phi|^2, x_g
  DvceArray1D<Real> gf_ainv;        // ((q*n1 + c)*n1 + gi): inverse of the mode-q matrix
  DvceArray1D<int> gf_okd;          // (q) the mode-q factorisation is usable
  void ImplicitGFInit();
  void ImplicitGFBuild();
  void ImplicitGFPre(int rc, int upd, Real c1, Real c2);
  void ImplicitGFAdd(int zc);
  bool impl_rho_direct = false;  // <rad_m1>/implicit_bcg_rho_direct (krylov_fuse = 3)
  int impl_dump_cyc;            // <rad_m1>/implicit_dump_op (debug, rad_m1_precond.cpp)
  bool impl_dump_done;
  void ImplicitDumpOp();
  void ImplicitOpCheck();
  // ---- the Picard pass count (bench/m1_picard_0923).  Since bench/m1_defaults_0923
  // lres_test = false, conv_est = true and lin_ew_max = 1e-2 (not predictor) are the
  // DEFAULTS for closure = eddington | vet_sc | tau (the OFF settings stay the defaults
  // for m1 | minerbo | kershaw); the OFF settings reproduce the earlier path bitwise.
  int impl_dtrace;              // <rad_m1>/implicit_det_trace (read only when named):
                                // bitwise hashes of the solve's arrays for the first N
                                // solves (run-to-run reproducibility diagnostic)
  void DetTrace(const char *tag);
  void DetTraceScalars(const char *tag, int n, const Real *v);
  bool impl_tsolve_opac;        // <rad_m1>/implicit_tsolve_opac (nc_cure_1002)
  int impl_tsolve_opac_start;   // implicit_tsolve_opac_start: first Picard pass
  int impl_tsolve_opac_mode = 0;      // implicit_tsolve_opac_mode (tsolve_root_fix_1003)
                                     // bit 1 nearest root, bit 2 slope guard, 0 = off
  Real impl_tsolve_opac_smin = -4.0;  // implicit_tsolve_opac_slope_min (bit 2 threshold)
  int impl_strace = 0;          // <rad_m1>/implicit_stall_trace (DEBUG, he_estall_1003)
  int impl_strace_p0 = 40;      // ...from this pass, the worst cell fixed (per rank)
  int impl_strace_max = 30;     // ...in at most this many solves
  int impl_strace_n = 0;        // solves traced so far
  int ncd_hot_loc = -1;         // set by ImplicitNCDump: its hottest in-block neighbour
  int impl_ncdump;              // <rad_m1>/implicit_nc_dump (DEBUG, default 0 = off):
                                // the local state of each rank's worst cell in the
                                // last N passes of a solve that hits implicit_maxit
  int impl_plog;                // <rad_m1>/implicit_picard_log: print one line per
                                // Picard pass for the first N solves (rank 0)
  Real bcg_r0rel;               // max|b - A x0|/max|b| of the last BiCGStab call
  bool impl_det;                // <rad_m1>/implicit_det_reduce: fixed-order Krylov sums
  DvceArray1D<Real> det_part;   // its level-1 partials and result (4*1024 + 4)
  bool impl_lres_test;          // <rad_m1>/implicit_lres_test (default false*): require
                                // the pass-to-pass transverse change lresid < lin_tol
  bool impl_src_stable = false;  // <rad_m1>/implicit_src_stable (stall_1002)
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
  int t2_dbg_adm;               // DEBUG dbg_t2_admiss: report the first N non-admissible
  int t2_dbg_adm_n;             // stages (which quantity <= 0, where); 0 = off
  Real t2_nstep, t2_nbe, t2_nfall;  // stage steps, BE steps, fallbacks (counters)
  DvceArray5D<Real> dbrow;      // DEBUG dbg_t2_admiss: the sp E row by term
  Real t2_dtprev;               // the dt of the previous step (vet_sc extrapolation)
  bool t2_vprev;                // vet_prev holds the tensor of the previous step
  bool t2_vext;                 // time2_vet_extrap (default false: D^n)
  // time2_tableau (ke-dt-0926; read only when named): 0 = hesdirk2 (Heun + H-ESDIRK2),
  // 1 = trbdf2 (TR-BDF2 implicit, stage order 2, with the RK2 hydro of c2 = 2 - sqrt 2;
  // the Driver sets the matching gam0/gam1/beta).  Stage coefficients of Time2FormStage:
  // stage 1 t2inc = t2_a21 dt K1; stage 2 radiation start = t2_w0 Y2 + t2_w1 U^n,
  // t2inc = dt (t2_c31 K1 + t2_c32 K2); the diagonal is g = 1 - 1/sqrt 2 in both.
  int t2_tab = 0;
  Real t2_a21 = 0.0, t2_w0 = 0.5, t2_w1 = 0.5, t2_c31 = 0.0, t2_c32 = 0.0;
  Real t2_lin_tol;              // time2_lin_tol: stage solves' implicit_lin_tol (<= 0:
                                // implicit_lin_tol x t2_lin_fac)
  Real t2_onep_s;               // time2_one_pass_safety: implicit_one_pass_safety of the
                                // stage solves (0: implicit_one_pass_safety)
  Real t2_lin_fac;              // time2_lin_tol_fac (default 1)
  Real t2_lin_save;             // the implicit_lin_tol put back after a stage solve
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
  // closure = vet_col under hesdirk2 (m1-sp-order2b): time2_vet_col = lag (D^n in both
  // stages, first order in time) | predict (ONE build per step at the predicted state
  // U^n + dt K1, at t^{n+1} to O(dt^2)) | rebuild (predict for stage 1, a second build
  // at the stage-1 solution Y1 for stage 2) | extrap (DIAGNOSTIC: D^n + r (D^n -
  // D^{n-1}), history in memory only, not in the restart)
  int t2_vcmode = 0;
  bool t2_fvnew = false;        // time2_vstage (rad_m1_time2.cpp)
  bool t2_wk = false;           // the last assembled E rows carried the gas work
  void ImplicitWorkRow(bool row);   // time2_vstage: the gas work in the E row
  bool t2_vcprev = false;       // extrap: vcol_prev holds the previous step's tensor
  Real t2_vcnfb = 0.0;          // predict/rebuild builds that fell back to U^n
  DvceArray4D<Real> vcol_prev;  // extrap: (m,k,j,i) f_K of the previous step
  DvceArray3D<Real> vcol_qprev; // extrap: (m,k,j) q of the previous step
  void Time2VetColAt(int which);    // vet_col build at 1 = U^n + dt K1, 2 = Y1
  void Time2VetColSaveY1();         // end of the stage-1 solve: E, T of Y1 (rebuild)
  void Time2VetColExtrap();
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
  // <rad_m1>/implicit_closure_thin_relax = C (runs_5c_thinstab; 0 = off, nothing is
  // allocated or read): under implicit_closure_lag = step, the (chi, n) a cell starts a
  // step with is relaxed from the one it used in the previous step with a weight
  // w = min(1, 2/(1 + G^2)), G = C max(chi', b/f)/tau_c + f chi'/min(chi, 1-chi), the
  // bound of the lagged-closure gain.  ctr_mem holds (chi, n1, n2, n3) of the last step.
  Real impl_ctrelax;
  bool ctr_init;
  DvceArray5D<Real> ctr_mem;
  // <hydro>/rad_signal_speed (hydro_newdt.cpp) reads opac(M1_OP_T) and the closure
  // tensor of the LAST M1 step, which a restart used to lose: the restarted run then
  // took one gas-speed-only step and its dt departed from the straight run's.  The
  // inputs travel in the restart file (kM1RssRstMagic) and are put back, or, on a fresh
  // start or a file without them, opac is filled from the initial state, before the
  // first dt (RssPrime, called from Driver::Initialize).
  //! the closure tensor the signal speed reads: 0 E/3, 1 M1Chi(F/cE), 2 tau_ten,
  //! 3 vet_cell (chi, n), 4 vet_cell (D_dd); the mode selection of hydro_newdt.cpp
  int RssMode() const;
  int RssRstNch(int mode) const;  // channels written: opac(M1_OP_T) + the tensor's
  //! the array and channel holding restart channel n of a given mode
  void RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c);
  void RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode);
  void RssPrime(bool restart);
  DvceArray5D<Real> rss_stage;   // (nmb, nch, k, j, i) read from the restart file
  int rss_stage_mode = -1;       // mode of rss_stage, -1 = nothing staged
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
  int impl_eccheck_every;       // <rad_m1>/implicit_eos_cache_check_every (m1-fast3):
                                // the check on every N-th cycle only (1 = every step)
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
  Real flr_ne, flr_ng;          // cell-solves whose solved E <= e_floor (floored) and
                                // whose written-back gas eint <= 0 (all ranks)
  Real vimp_emin;               // the smallest E the linear solve produced with it on
  DvceArray5D<Real> vmw, vmw_c;   // exchange scratch of the M1_NVIMP_X components
  MeshBoundaryValuesCC *pbval_vm;  // ...and its exchange object
  void ImplicitVimpBuild();     // the Jacobian of a(v') E' for this Picard pass
  int iw_gas;                   // first iw component of the M1_NIW_GAS block (see
                                // rad_m1_implicit.hpp); < 0 when neither option is on
  int impl_nec;                 // components of `ecache`
  DvceArray5D<Real> ecache;     // (m,nec,k,j,i) the frozen-density e(T) cache
  bool pin_report_newton_fb = false;   // <rad_m1>/report_newton_fb: print fallback cells
  DvceArray4D<Real> nfb_buf;   // report_newton_fb: the fallback counts copied to the host
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
  // STAGE CS0 (rad_m1_sph.cpp, m1-cs-implicit): the cubed sphere.  sph_geom is set there
  // too (x1 = r and the Coordinates areas, volumes and face distances exist on both
  // meshes), and cs_geom restricts it to transport = implicit, closure = eddington,
  // time_scheme = be, bicgstab + rbgs_fwd; the transverse faces are the plain two-point
  // form (no skew 1/sin, no cross term: stage CS1)
  bool cs_geom = false;
  // STAGE CS1 (rad_m1_sph.cpp CubedS1Init): the transverse operator on the cubed sphere.
  //   m1bcs     mb_bcs with BoundaryFlag::panel read as block, so that a panel seam is an
  //             OPEN face of the implicit operator (CS0 closed it); the kernels take it
  //             in place of mb_bcs only on cs (elsewhere the View they read is mb_bcs).
  //   cs_seam   per MeshBlock: its inner_x2, outer_x2, inner_x3, outer_x3 face is a seam.
  //   csg2/3    the canonical seam geometry (M1CsSeamGeom) of the x2 / x3 seam faces,
  //             (m, {dth, angm, sfoot, smid}, along-seam cell incl. ghosts, side lo/hi).
  // The implicit scratch exchanges (pbval_th, pbval_tq, pbval_kr) run without the
  // along-seam resample, so a seam ghost is the neighbour's own (mirror) cell.
  DualArray2D<BoundaryFlag> m1bcs;
  DualArray2D<int> cs_seam;
  DvceArray4D<Real> csg2, csg3;
  // STAGE CS3: the effective transverse geometry of the cs operator, for kernels that
  // read Coordinates through the sp branches (ImplicitVimpBuild): dxface with the x2/x3
  // two-point distances (centre arc x sin of the face angle; r angm at a panel seam)
  // and area with the canonical seam areas; x1 parts are copies
  DvceFaceFld4D<Real> cs_dxf_eff{"m1csdxf", 1, 1, 1, 1};
  DvceFaceFld4D<Real> cs_area_eff{"m1csarea", 1, 1, 1, 1};
  void CubedS1Init();
  //! the seam average of the stored face state f0x2/f0x3 (stage CS1, C5)
  void CubedSeamFaceAverage();
  int cs_seam_avg_n = 0;   // seam-face-average calls (diagnostic)
  bool cs_seam_avg_on = true;
  Real cs_seam_dmax = 0.0;   // max relative change of an x2 face value by the average
  ParameterInput *pin_cs_ = nullptr;
  //! print the Picard statistics of the implicit solver (from the destructor, rank 0)
  void ImplicitReport();
  // fixbundle-1009 <rad_m1>/dbg_energy_tally (default false, output only): rank 0 prints
  // at the start and end of every be solve the global sums of the gas total energy, E,
  // the reference-work accumulator and esrc (x volume): line ETALLY
  bool dbg_etally = false;
  void DbgEnergyTally(Real *o);
  // fixbundle-1009 F4: ApplyClosureLimits clips over ACTIVE cells (0 = E raised to
  // e_floor, 1 = |F| scaled to c E), every stage, and the event-log totals (collective)
  DvceArray1D<Real> acl_cnt_d;
  void EventTotals(Real *v);
  void OnePassAuto(const int t, const bool on, const bool one);
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
  void T2AdmissDebug(DvceArray5D<Real> uh, DvceArray5D<Real> u0_,
                     DvceArray5D<Real> t2i_, const Real cl, const Real ch,
                     const bool hh, const bool gq, const int t2s, const int it);
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
  // red_only (m1-fast3): y is already in yc (the overlap parts): the reductions only
  void ImplicitStencilOp(int xc, int yc, int red, Real *out, bool red_only = false);
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
  int ImplicitNCDump(int it, int tag, bool gasx, int igb, int igr, int igf, int igy,
                     int floc = -1);

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
  // vet_col_order2 (m1-sp-order2, default false): sp quadratic mu quadrature, S and chi
  // linear in r along the segments, the core rays from E and F at the inner face (the
  // mirrored incoming intensity at a reflecting inner x1)
  bool vcol_o2 = false;
  // vet_col_fk_min (m1-sp-order2b): the lower bound of f_K = K/J.  1/3 (the old clamp,
  // inherited from vet_sc's eigen/flux projection) or 0 (sp default: f_K < 1/3 is
  // physical behind a radiating front and D = diag(f, (1-f)/2, (1-f)/2) stays
  // realizable for 0 <= f <= 1)
  Real vcol_fkmin = 1.0/3.0;
  // vet_col_source (he-presn-m1 0929): the thermal source of the formal solution.
  // `relaxed` (default): a T^4 at the temperature T* the cell reaches after a LOCAL
  // backward-Euler exchange over the step (VcolRelaxedSource, rad_m1_vetcol.cpp); `gas`:
  // the start-of-step gas temperature as before (bitwise the old build)
  bool vcol_srelax = true;
  // vet_scatter (vet-scatter-1006, read only when named; default off = bitwise): ONE
  // switch for electron scattering: the opacity split (otab.scat, opac M1_OP_S) and the
  // closure source VetScatterMix.  vscat_form 1 = ma (Jiang 2021 eq. 6), 2 = absorption
  bool vscat = false;
  int vscat_form = 1;
  // vet_scatter_j (relaxed source mode only): false = current: J = E^n / the Picard
  // iterate next to B* (absorption form only); true = relaxed (the input default):
  // J = E*, the end-of-step pair of B*.  Set in the ctor whenever vet_scatter is on.
  bool vscat_jrel = false;
  // vet_col with a REFLECTING outer x1 (m1-sp-order2b): the incoming intensity at the
  // top face is the mirror of the outgoing one, I_in = b/(1 - a) per ray (b the outgoing
  // intensity of a vacuum-top sweep, a the ray's round-trip transmission), in a second
  // sweep; a is capped at vet_col_reflect_amax
  bool vcol_rtop = false;
  // vet_col_surface_face (m1-sp-order2b, sp): with vet_col_surface_q and
  // implicit_marshak_face = linear, the outer face flux is c q E_f with E_f the linear
  // face E and q = H(face)/J_f, J_f the formal J extrapolated to the face in the same
  // way (r^2 J linear from the top two shells): F = c H(face) (E/J)_face, second order
  // (the cell form c H(face) (E/J)_top-cell is first order when E/J varies).  vcol_jca,
  // vcol_jcb: the extrapolation weights of the top face (M1SphMarshakCoef), vcol_jcap
  // the limiter's 2 (r_c/r_f)^2
  bool vcol_sqf = false;
  Real vcol_jca = 0.0, vcol_jcb = 0.0, vcol_jcap = 0.0;
  Real vcol_amax = 0.99;
  DvceArray2D<Real> vcol_abuf; // (ray, column): transmission, then the top intensity
  int vcol_dump_every;         // vet_col_dump_every (0 = off)
  bool vcol_built;             // a tensor exists (vet_col_every > 1 skips builds)
  std::string vcol_dump;       // vet_col_dump: file prefix of the column dump
  Real vcol_time, vcol_ncall, vcol_nskip;
  DvceArray2D<Real> vcol_seg;  // (ray, shell l): path length from shell l to l+1 / top
  DvceArray2D<Real> vcol_mu;   // (ray, shell): mu of the ray at the shell
  DvceArray2D<Real> vcol_w;    // (ray, shell): hemisphere quadrature weight (sum 1)
  DvceArray2D<Real> vcol_ray;  // (ray, 0..3): L (first shell), type, seg0, a_t | mu_face
  // vet_col_order2 (sp): the mean over the segment of g(u) = (r(u) - r_a)/(r_b - r_a),
  // u the fraction of the path from the inner end r_a (1/2 radial, ~1/3 from a tangent
  // point): S and chi are taken linear in r, not in the path, along every segment
  DvceArray2D<Real> vcol_gb;   // (ray, shell l): the segment from shell l to l+1 / top
  DvceArray1D<Real> vcol_gb0;  // (ray): the first segment (core: r_in -> shell 0;
                               // sub-ray: tangent point -> its shell)
  DvceArray1D<int> vcol_klast; // (shell): index of the last ray active at the shell
  DvceArray3D<Real> vcol_sgt;  // (shell l, ray, 2): vcol_seg and vcol_gb transposed, so
                               // the team sweep's threads (rays) read them coalesced
  DvceArray2D<Real> vcol_buf;  // (ray, column): the running intensity of each ray
  DvceArray5D<Real> vcol_mom;  // (m,5,k,j,i): J, H, K, S, chi of the solution (dump)
  std::vector<double> vcol_rc; // shell radii (host, for the dump)
  // vet_col_surface_q (runs_5e): the outer-x1 Marshak q per column from the formal
  // solution, q = H(top face)/J(top cell), lagged like the tensor; vet_col_team: the
  // build as one team per column (rays over the team's threads, fixed-order moments)
  bool vcol_sq, vcol_team;
  Real vcol_qmin, vcol_qmax;
  int vcol_lc;                 // shells per chunk of the team build (scratch budget)
  int vcol_ts, vcol_lcin;      // vet_col_team_size (0 = AUTO), vet_col_chunk (0 = auto)
  DvceArray3D<Real> vcol_q;    // (m, k, j): the column's Marshak q at the top face
  DvceArray1D<Real> vcol_muf;  // (ray): mu at the top face
  DvceArray1D<Real> vcol_wf;   // (ray): hemisphere weight at the top face
  void VetColInit();           // checks, ray tables, buffers (first TauClosureInit)
  void VetColBuild();          // the formal solution -> tau_ten (chi, n)
  void VetColSpread(int ncall);   // vet_col_spread (CS2, gate T-S6)
  int vcol_spread = 0;
  void VetColBuildTeam(bool dmp);   // the same, one team per column (vet_col_team)
  void VetColReport();
  void VetColDumpColumn(int ncall);

  // ---- <rad_m1>/vet_col_lat (m1-vetcol-lat, rad_m1_vetlat.cpp; design
  // rt_design_1003/SP_SC_G0.md sect. 4, plan C1): the LATERAL CORRECTION of vet_col from
  // a few-angle short-characteristics sweep on the mesh (local mu x psi about r_hat).
  // Every direction is swept twice: through the 3-D field and through the column's own
  // laterally homogeneous field (the TWIN, the same arithmetic with the own column read
  // at every lateral point); their moment difference is the correction, exactly 0 on a
  // laterally uniform state.  Stored in tau_ten slots M1_TT_LAT0.. (the mesh basis):
  // dD_rr is folded into the implicit diagonal (chi of slot 0), D_r,lat enters the face
  // equations as a LAGGED Picard term (M1SphLat).  Read only when named; sp wedge only
  // (cs after CS2).
  bool vlat_on = false;        // vet_col_lat
  // vet_col_lat_offdiag: what D_r,lat does in the solve.  0 none, 1 lagged (a Picard
  // term, C1), 2 operator (default, C2: in the Krylov operator, M1SphLat of the
  // Krylov vector, through the stored stencil and VetLatOp)
  int vlat_odm = 2;
  bool vlat_now = false;       // the term is on in this step (positivity fallback)
  Real vlat_nfall = 0.0;       // steps dropped to `none` by the positivity guard
  Real vlat_odmax = 0.0;       // max |D_r,lat| of the last sweep (all ranks)
  // the operator form is active: on in this step and a non-zero D_r,lat (a zero one
  // leaves the operator, its stencil and the right-hand side untouched, bitwise)
  bool VlatOp() const {return vlat_now && (vlat_odm == 2) && (vlat_odmax > 0.0);}
  int vlat_nmu = 4, vlat_npsi = 4, vlat_every = 1, vlat_iinit = 3;
  Real vlat_taucut = 30.0;     // vet_col_lat_taucut: the sweep starts where every column
                               // has tau_top >= this (global shell index vlat_icut)
  int vlat_icut = 0;
  Real vlat_taumax = 0.0;      // vet_col_lat_taumax: depth taper (0 off)
  Real vlat_taumin = 0.0;      // vet_col_lat_taumin: thin-top cap (0 off)
  DvceArray4D<Real> vlat_wt;   // (m, k, j, i): the taper weight
  bool vlat_ready = false;     // VetLatInit done
  int vlat_nbuild = 0;         // builds done (the first one iterates vlat_iinit times)
  Real vlat_time = 0.0, vlat_ncall = 0.0, vlat_nclamp = 0.0;
  std::string vlat_dump;       // vet_col_lat_dump: per-rank dump prefix (gate c)
  int vlat_dump_every = 0;     // vet_col_lat_dump_every (0: the first build only)
  DvceArray5D<Real> vlat_i;    // (m, 2 nmu npsi, k, j, i): 3-D intensities
  DvceArray5D<Real> vlat_t;    // (m, 2 nmu npsi, k, j, i): twin intensities
  DvceArray5D<Real> vlat_t_c;
  DvceArray5D<Real> vlat_d;    // (m, 2 nmu npsi, k, j, i): 3-D minus twin, ghosts
                               // lagged (the lateral part of the inflow)
  DvceArray5D<Real> vlat_d_c;
  DvceArray5D<Real> vlat_cs;   // (m, 2, k, j, i): ln chi, ln S (ghosts exchanged)
  DvceArray5D<Real> vlat_cs_c;
  DvceArray5D<Real> vlat_lx;    // (m, 2, k, j, i): D_rt, D_rp for their ghost exchange
  DvceArray5D<Real> vlat_tt_c;  // its coarse twin
  DvceArray3D<Real> vlat_geo;  // (shell l, 2 a + dir, 0..3): per-shell ray geometry
  DvceArray1D<Real> vlat_mu;   // (a): Gauss nodes on [0, 1], descending
  DvceArray1D<Real> vlat_w;    // (a): hemisphere weights, sum 1/2 (J = sum w I)
  DvceArray1D<Real> vlat_cnt;  // (1): lateral reads clamped to the ghost band
  std::vector<double> vlat_nodes;   // the mu nodes (host)
  int vlat_geo_icut = -1;       // the first shell vlat_geo was built for
  MeshBoundaryValuesCC *pbval_vl = nullptr;   // vlat_i
  MeshBoundaryValuesCC *pbval_vs = nullptr;   // vlat_cs
  MeshBoundaryValuesCC *pbval_vt = nullptr;   // tau_ten (all slots)
  void VetLatInit();
  void VetLatBuild();
  void VetLatSweep(const int stage);   // 0 twin, 1 3-D, 2 moments
  void VetLatTTGhosts();
  void VetLatOdMax();
  void VetLatExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                      MeshBoundaryValuesCC *pb);
  void VetLatDump();
  // vet_col_lat_offdiag = operator: y += sgn L_lat(x) (the direct form, legacy
  // operator path and the right-hand side), and the same terms added to the stored
  // 19-point stencil of the pass (ImplicitStencilBuild)
  void VetLatOp(int xc, int yc, Real sgn);
  void VetLatStencilAdd();

  // ---- <rad_m1>/vet_gd (m1-vet-gd, rad_m1_vetgd.cpp; design
  // rt_design_1003/SC_PROPER_SP.md, candidate b): the Eddington tensor from an SC sweep
  // along GLOBALLY FIXED (Cartesian) HEALPix directions, handed over through the
  // vet_col_lat interface (sets vlat_on; the vet_col_lat keys apply).  Read only when
  // named; sp wedge only.
  bool vgd_on = false;         // vet_gd
  int vgd_nside = 3;           // vet_gd_nside: 12 nside^2 directions (3 -> 108)
  bool vgd_tan = true;         // vet_gd_tangential: the lagged M1SphTan term
  int vgd_iter = 1;            // vet_gd_iter: sweeps per ordinary build
  int vgd_smooth = 0;          // vet_gd_smooth: lateral 1-2-1 passes on dD_rr (diag.)
  DvceArray4D<Real> vgd_sm;    // (m, k, j, i) scratch of VetGdSmooth
  bool vgd_bandx = false;      // vet_gd_band_exit: shortened segments for reads
                               // beyond the ghost band instead of clamping
  Real vgd_taumin = 0.0;       // vet_gd_thin_taumin > 0: lateral parts tapered to 0
                               // from tau_top = 10 taumin to taumin (far thin top)
  bool vgd_replace = false;    // vet_gd_replace: no vet_col sweep (D = I/3 below the
                               // first gd shell, q from the gd top cell)
  int vgd_n = 0;
  int vgd_ls = 0;              // vet_gd_ls = 6/8/10/12: level-symmetric LQ_N set instead
  int vgd_gl_nmu = 0;          // vet_gd_gl_nmu > 0: product set GL(nmu in n_z) x nphi
  int vgd_gl_nphi = 12;        // vet_gd_gl_nphi (C_nphi about z)
  bool vgd_gl_stag = true;     // vet_gd_gl_stagger: alternate rings half a step
  int vgd_rot = 0;             // vet_gd_rotate_every: z-rotation every N cycles
  Real vgd_alpha = -1.0;       // the current z-angle (-1: tables not built)
  std::vector<double> vgd_base;   // the unrotated set (x, y, z, w)
  DvceArray2D<Real> vgd_dir;   // (d, 0..3): n_x, n_y, n_z, weight (sum 1)
  // (m, d, k, j, i): intensities, ghosts lagged.  accel-1009: LayoutLeft (i slowest), so
  // that one shell of the per-shell sweep and halo is one contiguous ~nd*nk*nj slab (the
  // LayoutRight array spread every shell over all of its GBs: strided, TLB-bound); the
  // values and every arithmetic are unchanged (bitwise)
  // accel-1009 (memory): RAGGED per-shell band.  Shell i stores only the lateral band
  // its own data needs (wi(i) = vgd_wsh, <= vgd_w; the dense array carried the deepest
  // band of all shells for every shell: 70 % of the device memory of AG Car A), and a
  // block only on the sides where a neighbour (edge or corner) is NOT on this rank (an
  // on-rank neighbour is read in its interior, never through the band).  Index space
  // unchanged: (m, d, k, j, i) with k, j in the band index space of depth vgd_w; one
  // slab per (m, i), d fastest, then k, then j.  An access outside
  // the shell's slab (never, by the band construction) goes to one dummy element and is
  // counted (VgdRag::oos, reported at the end): the values are bitwise those of the
  // dense array as long as the count is 0.
  struct VgdRag {
    DvceArray1D<Real> d;
    Kokkos::View<int64_t**, DevMemSpace> off;  // (m, i): first element of slab (m, i)
    DvceArray1D<int> wi;                       // (i): band depth of shell i
    DvceArray2D<int> sb;                       // (m, 4): band on side k-, k+, j-, j+
    DvceArray1D<Real> oos;                     // (0): out-of-slab accesses
    int nmb = 0, n = 0, nx2 = 0, nx3 = 0, w = 0, c1 = 0;
    int64_t dummy = 0;
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
  using VgdIView = VgdRag;
  VgdIView vgd_i;
  void VgdRagAlloc(VgdRag &a);
  DvceArray5D<Real> vgd_i_c;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_wall;    // (m, k, j)
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_map;    // (m, k, j, d)
  DvceArray4D<Real> vgd_mr;    // (m, k, j, d): n . r_hat of the value stored there
  // vet_gd_wall_interp (read only when named, default false): wall ghosts from the 3
  // nearest same-branch directions with weights (m, k, j, d, 3)
  bool vgd_wint = false;
  // vet_gd_twin (read only when named, default false): uniform-state twin subtraction
  bool vgd_twin = false, vgd_noq = false;
  bool vgd_twfull = false;    // vet_gd_twin_full: LAT0 -= full twin LAT0 (no shell mean)
  DvceArray1D<Real> vgd_twm, vgd_twm2;
  DvceArray5D<Real> vgd_twl, vgd_cs0;
  // vet_gd_twin_fuse (accel-1009): the twin's own source (ln chi, ln S of the shell
  // means) and intensities, swept in the main sweep's kernels; vgd_sw2 = this sweep
  // carries the twin; second compact halo buffers for it
  bool vgd_twfuse = false, vgd_sw2 = false;
  bool vgd_twdet = false;      // vet_gd_twin_det: fixed-order (reproducible) shell means
  bool vgd_twseq = false;      // vet_gd_twin_lowmem: fused twin swept first, one array
  DvceArray5D<Real> vgd_cst;
  VgdIView vgd_itw;
  DvceArray1D<Real> vgd_csb2, vgd_crb2, vgd_rbuf2;
  void VetGdTwinFusedMoments();
  Real vgd_ttwin = 0.0;
  void VetGdTwin(const int stage);
  Kokkos::View<int *****, LayoutWrapper, DevMemSpace> vgd_m3;
  DvceArray5D<Real> vgd_w3;
  MeshBoundaryValuesCC *pbval_gd = nullptr;
  // EXACT per-shell lateral halo (decomposition-invariant sweep): vgd_i and vgd_cs carry
  // a lateral ghost band of vgd_w cells (the deepest reach of an upwind point), filled
  // shell by shell from the 8 lateral neighbours (on-rank copies + MPI)
  int vgd_w = 0;
  DvceArray5D<Real> vgd_cs;    // (m, 2, k, j, i): ln chi, ln S with the vgd_w band
  DvceArray1D<int> vgd_hloc;   // (8 nmb): local index of the slot's neighbour, -1 remote
  std::vector<int> vgd_hrank, vgd_hlid;   // (8 nmb): the remote neighbour's rank, lid
  DvceArray1D<Real> vgd_sbuf, vgd_rbuf;   // flat message buffers, one piece per remote
  // (block, slot), the pieces to/from one rank contiguous (one message per rank pair)
  // in the order of the RECEIVER's (block, slot); offsets per band depth ws in units of
  // nv*ni: (ws, 8 m + o), -1 for an on-rank slot
  DvceArray2D<int> vgd_soff, vgd_roff;
  std::vector<int> vgd_prk;                  // partner ranks
  std::vector<std::vector<int>> vgd_pdsp, vgd_pscnt, vgd_prdsp, vgd_prcnt;   // (ws, p)
  int vgd_maxcnt = 0;
  bool vgd_hmpi = false;
  Real vgd_thalo = 0.0, vgd_tmpi = 0.0, vgd_ttab = 0.0;
  // halo diagnostics (rank 0; clamps all ranks)
  Real vgd_tpost = 0.0, vgd_nexch = 0.0, vgd_nbyte = 0.0, vgd_nclamp_all = 0.0;
  bool vgd_time_halo = false;
  // vet_gd_rebuild_every = k > 0 (read only when named; 0 = the tensor lagged to the
  // start of the step): rebuild the gd tensor inside the implicit solve at Picard passes
  // k, 2k, ... from the iterate's E (M1_IW_EP) and T (M1_IW_TP)
  int vgd_rbe = 0;
  // vet_gd_tan_operator (read only when named, default true): the tangential cross terms
  // (M1SphTan) join D_r,lat in the implicit operator under vet_col_lat_offdiag = operator
  bool vgd_tanop = true;
  // thin-top cap shape (read only when named): ramp decades, smoothstep, parts mask
  Real vgd_tdec = 1.0;
  bool vgd_tsmooth = false;
  int vgd_tparts = 3;
  int vgd_seam = 0;             // vet_gd_seam_mask (TEST)
  bool vlat_src_ep = false;
  DvceArray4D<Real> vgd_fk0;   // (m, k, j, i): slot 0 before the fold
  Real vgd_nrb = 0.0, vgd_trb = 0.0;
  void VetGdIterRebuild();
  void VetGdMms();
  void VetGdRealDiag();
  int vgd_rdiag = -1, vgd_rdcnt = 0;   // env VGD_TIME_HALO=1: fence before the halo timer
  void VetGdHaloInit();
  // vet_gd_halo_compact = 0 (off), 1 (pass branch), 2 (+ ray reach): compact messages
  int vgd_hcomp = 0;
  bool vgd_capped = false, vgd_allow_clamp = false;
  Real vgd_nclamp_seen = 0.0;
  bool vgd_hinw = true;
  int vgd_nb2 = 1, vgd_nb3 = 1;
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_lxy;   // (m, 0/1): lx2, lx3
  DvceArray1D<Real> vgd_csb, vgd_crb;
  // vet_gd_halo_pipe (default true, bitwise): the compact mask of the NEXT shell is built
  // while this shell's messages travel; 2 slots of flags/positions (vgd_hfs/hfr) and of
  // the partner boundaries in pinned host memory (vgd_hpb), tagged (i, inw, ws, icut,
  // sweep id)
  DvceArray1D<int> vgd_hfs[2], vgd_hfr[2];
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hpb[2];
  int vgd_htag[2][5] = {{-1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1} };
  int vgd_hlast = 1, vgd_hsweep = 0;
  // vet_gd_halo_cache_mb (accel-1009; default 0 = off): byte budget (MB per rank) of a
  // cache of the compact-halo masks per (pass, shell), valid while the direction set,
  // the cut and the depth are unchanged (one prep per shell per rotation window instead
  // of one per shell per sweep); bitwise the same masks
  struct VgdHcEntry {
    DvceArray1D<int> fs, fr;
    Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> pb;
    Real alpha = -2.0;
    int scut = -1, ws = -1;
  };
  int vgd_hc_mb = 0;
  std::vector<VgdHcEntry> vgd_hce;
  size_t vgd_hc_bytes = 0;
  Real vgd_hc_nmade = 0.0;
  int VetGdHcGet(const int i, const bool inw, const int ws);
  int vgd_hnext[3] = {-1, 0, 0};
  bool vgd_hpipe = true;
  void VetGdHcPrep(const int slot, const int i0, const bool inw0, const int ws);
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_pbd;   // (ws, 2 (np+1))
  std::vector<double> vgd_r1v, vgd_r1f;
  template <class V>
  void VetGdHaloCompact(V &a, const int nv, const int i0, const int ws, V *b = nullptr);
  std::vector<int> vgd_wsh;    // (i): band depth the shell's data needs (<= vgd_w)
  std::vector<int> vgd_wsi, vgd_wso;   // (i): the same per pass (inward, outward)
  template <class V>
  void VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                 const bool mapd, V *b = nullptr);
  Real vgd_tsrc = 0.0, vgd_tswp = 0.0, vgd_texc = 0.0, vgd_tmom = 0.0;
  void VetGdInit();
  Real VetGdAngle(const int cyc) const;
  void VetGdTables(const Real alpha);
  void VetGdBuild();
  void VetGdSweep();
  void VetGdWall(const int i0, const int i1);
  void VetGdMoments();
  void VetGdSmooth();
  // vet_gd_async (GD_ASYNC.md): the sweep of build n runs on its own device instance
  // and host thread while the rest of cycle n runs; the D of cycle n is the moments of
  // build n-1 (one cycle lag).  The pending build's source (vlat_cs) and cycle travel
  // in the restart file (kM1VgdRstMagic).  vgd_cur: the instance the sweep/halo use.
  bool vgd_async = false;
  bool vgd_ainl = false;       // the sweep inline at launch (host backends, or env
                               // VGD_ASYNC_INLINE=1): same numbers, no overlap
  bool vgd_apend = false;      // a swept build whose moments are not yet taken
  bool vgd_afly = false;       // the helper thread is running
  int vgd_acyc = -1;           // ncycle of the pending build (its z-angle)
  int vgd_pcut = -1;           // vlat_icut of the pending build (async)
  int vgd_scut = 0;            // vlat_icut the running sweep uses (sync: = vlat_icut)
  std::thread vgd_athr;
  DevExeSpace vgd_cur, vgd_ex;
  Real vgd_tjoin = 0.0, vgd_nasync = 0.0;
  DvceArray5D<Real> vgd_rst;   // (nmb, 2, k, j, i): vlat_cs of the pending build, staged
  int vgd_rst_cyc = -1;        // from the restart file (-1: nothing staged)
#if MPI_PARALLEL_ENABLED
  MPI_Comm vgd_comm = MPI_COMM_WORLD;
#endif
  void VetGdSource();
  void VetGdPost();
  void VetGdBuildAsync();
  void VetGdAsyncJoin();
  int VetGdAsyncRstNch() const {return (vgd_async && vgd_apend) ? 2 : 0;}
  void VetGdAsyncRstPack(DvceArray5D<Real> &dst, int nmb);
  void VetLatCsFinish();

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
