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

#include <array>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <future>
#include <iostream>
#include <thread>  // NOLINT(build/c++11)
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
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

//! unify-sc-gd-1012: ONE canonical key name per concept for both formal solutions
//! (vet_sc, vet_col / vet_gd), the older closure-specific spellings accepted as aliases.
//! Returns the name to read in <rad_m1>: the canonical key if the input has it, else the
//! first alias it has, else the canonical key (a GetOrAdd* on it then records the
//! resolved default).  An alias is read under its own name, so an input (or a restart
//! file) that names only the old spelling behaves exactly as before.  The canonical key
//! wins over an alias given with a different value (printed).
inline std::string M1KeyName(ParameterInput *pin, const std::string &canon,
                             const std::vector<std::string> &alias) {
  const bool hc = pin->DoesParameterExist("rad_m1", canon);
  for (const auto &a : alias) {
    if (!pin->DoesParameterExist("rad_m1", a)) {continue;}
    if (!hc) {return a;}
    if (pin->GetString("rad_m1", a).compare(pin->GetString("rad_m1", canon)) != 0 &&
        global_variable::my_rank == 0) {
      std::cout << "<rad_m1> NOTE: " << canon << " = " << pin->GetString("rad_m1", canon)
                << " overrides its old spelling " << a << " = "
                << pin->GetString("rad_m1", a) << std::endl;
    }
  }
  return canon;
}

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
constexpr int M1_VET_NL  = 26;  // thinsw-1010: sum w |I - S| (C1 non-locality), fs only
constexpr int M1_VET_EPS = 27;  // thinsw-1010: thermal fraction min(kappa_P/chi, 1) (fs)
constexpr int M1_VET_NC  = 28;
// tau_ten slots of the LATERAL correction (vet_col_lat; the cs CS2 interface): dD_ab in
// the mesh basis, order rr, r-a, r-b, aa, ab, bb (sp: a = theta, b = phi)
constexpr int M1_TT_LAT0 = 4;
constexpr int M1_TT_NLAT = 6;
// gdthin-1011: the extra vgd_hr components of the thin switch on vet_gd (VetGdHalfRange)
constexpr int M1_VGD_HRC1 = 6;     // C1 = sum w |I - S| / J (twin-corrected as r+-)
constexpr int M1_VGD_HREP = 7;     // eps = min(kappa_P/chi, 1) (1 where chi = 0)
constexpr int M1_VGD_HRJ  = 8;     // J (E units; twin-corrected by J_tw/<J_tw>)
constexpr int M1_VGD_HRNX = 9;
// vgdlat-1012 vet_gd_symref: the per-shell quantities of the 1-D reference (vgd_srm)
constexpr int M1_VGD_SR_F = 0;    // K/J
constexpr int M1_VGD_SR_HP = 1;   // H+_r/J
constexpr int M1_VGD_SR_HM = 2;   // H-_r/J (<= 0)
constexpr int M1_VGD_SR_HL = 3;   // H+_theta/J = H+_phi/J = -H-_theta/J
constexpr int M1_VGD_SR_C1 = 4;   // int |I - S|/J
constexpr int M1_VGD_SR_J = 5;    // J (E units)
constexpr int M1_VGD_SR_N = 6;
#ifdef M1_VGD_NMAX_BUILD
constexpr int M1_VGD_NMAX = M1_VGD_NMAX_BUILD;   // rotavg-1012: cmake -D Athena_M1_VGD_NMAX
#else
constexpr int M1_VGD_NMAX = 432;   // vet_gd: at most 432 directions.  prodcand-1011: kept at 432
#endif
                                    // (sizes the per-thread Real buf[M1_VGD_NMAX] in VetGdWall; 2048
                                    // would spill on A100).  The 4x-GL angular study (sptests P2)
                                    // needs a local build with 2048.
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
  // <rad_m1>/opac_abs_rho_max (default 0 = off): code densities below amask_rho get no
  // absorption (kappa_P = kappa_E = 0) and the transport opacity amask_kes (code units
  // per mass, = opac_abs_kappa_s [cm^2/g] x kunit); see rad_m1.cpp
  Real amask_rho = 0.0;
  Real amask_kes = 0.0;
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

  // fin4-dtretry (driver/dt_retry.cpp): <time>/dt_retry.  dtr_on: the driver checkpoints
  // the cycle and takes a failed one again with a smaller dt; then a solve that would
  // stop the run (DIVERGED) or end NON-CONVERGED outside a hesdirk2 stage (whose own redo
  // is the backward-Euler step) raises dtr_fail with the reason instead.
  bool dtr_on = false;
  bool dtr_fail = false;
  int dtr_why = 0;           // 1 NON-CONVERGED, 2 DIVERGED, 3 DIVERGED (masked), 4 eint<=0, 5 time2_reject
  Real dtr_resid = 0.0;
  // t2befb-1020: armed by DtRetry (driver/dt_retry.cpp) once <time>/dt_retry_t2_be cuts
  // (default 2) did not cure a time2_reject (the cycle is then retaken at its own dt), or
  // when the previous cycle took that redo: a non-admissible hesdirk2 stage is then redone
  // with backward Euler (as time2_reject = be) instead of raising another retry
  bool dtr_t2be = false;
  bool dtr_t2hit = false;    // t2befb-1020: this attempt took that backward-Euler redo
  template <class V> void DtRetryVisit(V &v);   // rad_m1_dtretry_visit.hpp (generated)
  void DtRetryQuiesce(bool restoring);          // no host thread writes while restoring
  void DtRetryAfterRestore();                   // async table prefetch state made safe

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
  // icrad-1011: under time_scheme = hesdirk2 the BACKWARD-EULER solves (the BESTORE
  // start step, the fallbacks, the ic_rad passes) read the polluted E^n + dt esrc too
  // (F1 was be-only).  vet_source_noesrc_h2be (default true; read only when named):
  // noes also for those solves; the stage solves keep iw EN (predict puts P there)
  bool vsrc_noes_h2be = true;
  bool VsrcNoes() const;

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
  // unify-sc-gd-1012 (#2): the x1 boundary intensity of the formal solutions, one rule
  // for vet_sc and vet_gd (VetBcInit, rad_m1.cpp).  vet_bc_inner (alias vet_gd_inner on
  // vet_gd): the entering intensity at a MARSHAK or FLUX inner x1 face, 0 diff (the
  // diffusion intensity), 1 emit (the face's own emission: the Marshak bath, or the
  // black sphere 4 F/c of a flux face), 2 auto (emit at a Marshak face, diff at a flux
  // face).  vet_bc_reflect: a REFLECTING x1 end, vet_bc_mirror = the entering rays
  // carry the mirrored leaving intensity (vet_sc; vet_gd: the vet_col column mirrors,
  // the gd sweep keeps diff at x1min and refuses a reflecting x1max).
  int vet_bc_inner = 0;
  bool vet_bc_mirror = false;   // mirror at a reflecting x1min (vet_bc_reflect = mirror)
  bool vet_bc_mirror_hi = false;   // ... at a reflecting x1max (mirror or auto, vet_sc)
  void VetBcInit(ParameterInput *pin, const bool gd);
  bool VetBcEmit() const;    // the inner rule at the current implicit_bc_x1min is emit
  Real VetBathLo() const;    // vet_sc: the entering inner intensity of emit (-1: diff)
  bool vet_bc_face = false;  // vet_bc_diff = face: inner diffusion I at the face (vet_sc,
                             // vet_gd, vet_col_lat, vet_col core; gdface-1011)
  bool vet_seg_split = false;  // vet_seg = split (segfix-1011): the SC segment u -> 0
  bool vet_seg_emit = false;   // vet_seg = emit: split, limiter only where S_u > S_0
                               // as two sub-segments split at the face
  bool vet_seg_jump = false;   // vet_seg = jump (coreL-1012): emit + the contrast rule
  // the segment mode VetSegSplit takes: 0 split, 1 emit, 2 jump (rad_m1_vetseg.hpp)
  int VetSegMode() const {return vet_seg_emit ? (vet_seg_jump ? 2 : 1) : 0;}
  bool vet_x1per;            // vet_x1_periodic (default false): periodic x1 sweep
  int vet_x1npass;           // ...and its number of passes (runs_3r_radwave)
  int vet_path = 0;          // open-od-1013: 1 = the x1 step subdivided at the row nodes
  int vet_interp = 0;        // unify-sc-gd-1012 (#6): upwind stencil 0 lin, 1 log, 2 corner
  int vet_mirror_npass = 4;  // unify-sc-gd-1012: passes of a cavity mirrored at both x1
                             // ends (vet_bc_reflect; vet_x1_npass when named)
  bool vet_x1src = false;    // vet_x1_periodic_inflow = source (predfix-1011)
  // bubbles-1012: a NON-periodic x2 (x3) mesh boundary (outflow, reflect, user, ...) is
  // a VACUUM face for the formal solution: a ray whose foot leaves the box through it
  // enters with I = 0 (no incoming radiation; single MeshBlock sweep only).  Periodic
  // x2 (x3) keeps the wrap, bit for bit.
  bool vet_x2open = false, vet_x3open = false;
  // bubbles-1012: the vet_sc angle set and its ROTATION on the vet_gd clock (the same
  // keys vet_gd_rotate_dt | _t0 | _align | _stretch | _slope and the same fresh default
  // auto = 20 x the first dt; VetGdBlock / VetGdAngleBlk / VetGdAlignDt / VetGdRotStep).
  // Double-Gauss: the azimuths about x1 turn by (2 pi/nphi) q; ring: the in-plane angles
  // turn by 0.8 (q - 1/2) (2 pi/N) (never past +-x2, where an x1 sweep has no ray);
  // q = alpha/(pi/2) in [0, 1) the low-discrepancy fraction of the block.
  int vet_ring_n = 0, vet_ring_g = 1;   // vet_ring, vet_ring_groups (0: double-Gauss)
  bool vsc_rot = false;                 // vet_sc rotates (physical clock)
  Real vsc_beta = 0.0;                  // the current rotation of the vet_sc set
  void VetBuildAngles(const Real beta);
  void VetRotInit(ParameterInput *pin);   // the shared rotation keys (rad_m1_vetgd.cpp)
  void VetRotUpdate();
  Real VetRotBeta(const int64_t blk) const;   // rotavg-1012: the set of block blk
  void VetRotSet(const Real beta);
  void VetScRotAvgSetup();     // rotavg-1012: vet_sc's averaged fields (vet_cell)
  void VetShortCharSweep();    // the build with the current set
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
  // picard-1014 implicit_tol_rel = kappa (default 0 = off): a solve that starts from the
  // step predictor stops at max(implicit_tol, kappa r_0), r_0 = its first change
  Real impl_tol_rel = 0.0;
  Real impl_tol_rel_max = 1.0e-4;   // implicit_tol_rel_max: cap on tol_eff
  Real impl_tolr_sum = 0.0, impl_tolr_n = 0.0;
  // picres-1017 implicit_final_resid_log (read only when named; default false = off).
  // DIAGNOSTIC: after every converged Picard solve ONE more pass is assembled at the
  // accepted iterate (the "shadow" pass: the next pass up to its linear solve) and the
  // nonlinear residual r = b(x) - A(x) x of that pass is evaluated per row and logged
  // (PICRES lines, rank 0) and summed (ImplicitReport); the state the shadow pass
  // touched is then put back, so the step continues exactly as without it.  The stop
  // test is NOT changed.  Needs implicit_solver = bicgstab.
  bool impl_frlog = false;
  bool impl_fr_sh = false;      // the current Picard pass is the shadow pass
  int impl_fr_how = 0;          // stop of the converged pass: 1 plain, 2 conv_est, 3 onep
  int impl_fr_it = 0;           // its number of passes
  Real impl_fr_tol = 0.0;       // its tol_eff
  Real impl_fr_chg = 0.0;       // its last change (the Picard resid)
  Real impl_fr_chp = 0.0;       // the change of the pass before (-1: none)
  double impl_fr_t0 = 0.0;      // wall time at the start of the snapshot
  double impl_fr_tsnap = 0.0;   // wall time of the snapshot copy
  Real impl_fr_sum[32] = {};    // summary accumulators (ImplicitFinalResid)
  DvceArray5D<Real> fr_iw, fr_opac, fr_ifw, fr_ifw2, fr_ifw3;
  DvceArray4D<Real> fr_ktd, fr_f0x1n, fr_f0x2n, fr_f0x3n, fr_cs2f1, fr_cs2f2, fr_cs2f3;
  DvceArray4D<Real> fr_f0x2, fr_f0x3, fr_thx2, fr_thx3;
  DvceArray1D<Real> fr_nsk;
  // picstop-1017 (Picard stop, OPT-IN; every key read only when named, default off =
  // the stop of picres-1017 bit for bit):
  //  implicit_conv_est_guard = K > 0: the contraction-estimate stop (implicit_conv_est)
  //    accepts only a pass whose own change is below K tol_eff;
  //  implicit_conv_est_nofloor = true: q of that stop is the ratio of the last two
  //    changes over the cells whose new E is above e_floor (a cell dropped to the floor
  //    makes a change ~1e10 that fakes q ~ 1e-10, picres_1017/NOTE.md);
  //  implicit_final_resid_stop = thr > 0 (and implicit_final_resid_stop_rel = s: the
  //    threshold is max(thr, s tol_eff)): a solve the change test accepts is accepted
  //    only if the scaled nonlinear residual max(rc, rT) of the accepted iterate over
  //    the cells above e_floor (ImplicitFinalResid, the shadow pass of
  //    implicit_final_resid_log) is <= the threshold; otherwise that assembly IS the
  //    next Picard pass (its linear solve follows) and the loop goes on, up to
  //    implicit_maxit (a check that fails with no pass left: NON-CONVERGED, the usual
  //    dt_retry / stage-redo handling).  Needs implicit_solver = bicgstab.
  //  Experiment 2 (picstop_1017, picres_1017/NOTE.md): nofloor brings the AG Car tail to
  //  the level of a tol 1e-10 reference; the residual stop is NOT recommended (in ~2-4 %
  //  of solves the change stalls at 1e-14 with rc > 1e-4: maxit, NC, dt retries).
  Real impl_ce_guard = 0.0;
  bool impl_ce_nofl = false;
  bool impl_frs_on = false;
  Real impl_frs = 0.0, impl_frs_rel = 0.0;
  int impl_fr_nchk = 0;           // residual checks of the current solve
  int impl_fr_it1 = 0;            // its passes at the first check
  // [0] checks [1] failed + continued [2] failed, no pass left (NC) [3] solves accepted
  // after >= 1 failed check [4] passes taken after a failed check [5] conv_est stops
  // refused by the guard [6] by nofloor (the plain q would have stopped) [7] nofloor
  // stops the plain q would have refused [8] logging time (fenced, s)
  Real impl_frs_sum[12] = {};
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
  // picard-1014 implicit_gas_newton_u4: emission-dominated cells take the gas Newton
  // step in u = T^4 (T' = (T^4 + 4 T^3 dT)^(1/4)); implicit_opac_lag = pred: the
  // opacities are evaluated at pass 0 (the predicted T) and frozen for the solve
  bool impl_gn_u4 = false;
  bool impl_opac_lagp = false;
  int impl_opac_lagn = 1;       // implicit_opac_lag_passes: updates in passes < N
  Real impl_ol_nlag = 0.0;      // solves run with the opacity lag (report)
  bool impl_allow_multid;       // run a 2-D/3-D set of INDEPENDENT x1 columns
  Real marshak_q;               // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min, ibc_x1max;     // M1_IBC_*
  // bubbles-1012: the transverse (x2, x3) physical faces of the implicit solve,
  // <rad_m1>/implicit_bc_x2min|x2max|x3min|x3max = auto (default: periodic on a
  // periodic mesh, else REFLECT = zero flux, as always) | reflect | marshak (the vacuum
  // free surface F_f = c marshak_q E of the boundary cell, outward; the formal solution
  // enters there with I = 0).  Read only when named.
  int ibc_x2min = 2, ibc_x2max = 2, ibc_x3min = 2, ibc_x3max = 2;
  Real iflux_x1min, iflux_x1max;  // the imposed face flux of M1_IBC_FLUX
  Real iebath_x1min, iebath_x1max;  // M1_IBC_MARSHAK: the INCIDENT bath
  // PER-COLUMN INFLOW columns of the outer x1 Marshak face (accretor-rhd-1008 S5, the
  // stream window's own radiation): a problem generator hands E_in(m,k,j) (nmb, n3, n2)
  // through SetX1maxBathColumns.  Columns with E_in > 0 are an optically thick INFLOW
  // face: comoving face flux 0 (no Marshak term) and the lab-frame enthalpy inflow
  // a E_in = v1 (1 + chi) E_in of implicit_bc_advect, i.e. F = (4/3) v E_in in the thick
  // limit; columns with E_in = 0 keep the Marshak face with the scalar bath (must be 0).
  // The face kernel stores the lab-frame advective energy flux through the outer face in
  // ibadv_x1max(m,k,j) (the end face's f0x1 is the comoving flux only).
  // Default: not set, bitwise.  (A first version used E_in as a Marshak bath, c q (E -
  // E_in): that shines an isotropic 65 kK wall into the thin ambient, c/(4 (4/3) v) ~450x
  // the advected radiation, and diverged in the first step, rhd/s5/smoke.)
  bool iebath_hi_col_on = false;
  DvceArray3D<Real> iebath_hi_col;
  // acc3d-1017: <rad_m1>/implicit_bc_x2min|max = flux: the x2 physical face carries an
  // IMPOSED radiative flux F_theta(m,k,i) (lab frame, + = toward larger theta), set by the
  // problem generator through SetX2FaceFlux (an equilibrium lateral flux, e.g. the von
  // Zeipel flux of a tilted stellar photosphere at a wedge wall); explicit in the solve
  bool ifx2_on = false;
  DvceArray3D<Real> ifx2_lo, ifx2_hi;
  DvceArray3D<Real> ibadv_x1max;
  // littests-1011 (implicit_beam_n > 0, vet_sc only): up to 2 collimated beams entering
  // through the x1min face on an x2 patch, carried by one SC ray (2-D: the +-mu3 mirror
  // pair).  ibm_i(m,q,k,j): the entering intensity of beam ray q at column (k,j) (local,
  // active index); ibm_f(m,k,j): its x1 face flux c sum w mu1 I, which replaces the
  // Marshak row of every column it is > 0 on (an imposed-flux face).
  int ibm_n = 0, ibm_nr = 0;
  int ibm_ray[4] = {-1, -1, -1, -1};
  DvceArray4D<Real> ibm_i;
  DvceArray3D<Real> ibm_f;
  void ImplBeamInit(ParameterInput *pin);
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
  int impl_blend_xthin_mode = 0;       // xthinfix-1009: 0 all, 1 beam, 2 beam_kn
  Real impl_blend_xthin_wmin = 0.0;    // xthinfix-1009: threshold of mode = beam
  Real impl_blend_xthin_r0 = 1.5;      // xthinfix-1009: Knudsen R0 of mode = beam_kn
  Real impl_blend_r0 = 1.5;       // hrup-1009: implicit_blend = knudsen, R0
  // thinsw-1010: <rad_m1>/implicit_thin_switch (M1_THSW_*, rad_m1_implicit.hpp): a
  // dt-free thin/thick weight frozen per solve (formal-solution |H|/J, its non-locality
  // C1, the E^n Knudsen number); implicit_thin_corr adds the lagged formal-solution flux
  int impl_thsw = 0;
  int impl_thsw_corr = 0;
  bool thsw_c1 = false;           // the sweep accumulates M1_VET_NL
  Real thsw_h0 = 0.3, thsw_h1 = 0.6, thsw_kn0 = 1.0, thsw_knp = 4.0;
  Real thsw_fsa = 0.0, thsw_fsj = 0.0, thsw_qp = 8.0;
  // thinfix-1011: implicit_thin_recon = js (M1ThinJs: half-range ratios from the face
  // mean of the formal-solution H+-), implicit_thin_bq = lin (outer Marshak q from the
  // face-extrapolated H+ of the top cell), the js factor limit
  bool thsw_js = false, thsw_bqlin = false;
  bool thsw_jse = false;   // js weighted by the thermal fraction eps of the face
  bool thsw_hctau = false;  // coreL-1012 implicit_thin_hcface = tau (M1ThinFace)
  bool thsw_jsmean = false;  // jsfix-1012 implicit_thin_js_face = mean (M1ThinJs)
  bool thsw_hctau_lat = false;  // ... and on the sp wedge's transverse faces (gd = all)
  Real thsw_jslim = 2.0;
  // unify-1012: implicit_thin_hcsf = jh: hc/hcs take the sf donor scale Jh as the cell's
  // formal-solution scale on faces with sf ratios (vet_sc); default j = the old J
  // (jh opt-in: +16 % in vacuum crossing beams, NOTE-unify)
  bool thsw_hcjh = false;
  bool thsw_bqgd = false;   // unify-1012: thin_bq lin also on vet_gd (opt-in)
  int thsw_jsgd = 1;   // unify-1012: js/jse on vet_gd, bit 1 x1, bit 2 lateral
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
  // cbfix-1011 implicit_hr_face = sf (default with vet_sc): the half-range face ratios
  // from the face mean of the formal-solution half-range fluxes over the DONOR SCALE
  // Jh = max(J, the outgoing face SC fluxes) wherever the cell's own J does not carry
  // the SC flux through its faces (M1HrSf); cell = the hrup-1009 cell ratios (opt-out)
  bool impl_hr_sf = false;
  int impl_hr_sfgd = 0;   // unify-sc-gd-1012 (#4) implicit_hr_face_gd: 1 x1, 3 all faces
  // sptests-1011 P2 (opt-in, vet_gd only): implicit_hr_excess = iso gives the part of the
  // lagged E^n that exceeds the formal-solution J of the cell the ISOTROPIC half-range
  // shape (+-1/4) instead of the formal solution's (M1HrExcess); none = off (default)
  bool impl_hr_xs = false;
  Real impl_hr_xs_cfl = 10.0;   // ... only where c dt/dx1 of the cell >= this
  int iw_hsf = -1;              // first iw component of the M1_NIW_HSF block
  double hsf_ncell = 0.0, hsf_nsf = 0.0;   // cell-solves, of which Jh > J (rank 0)
  bool hr_q = false;            // the per-column outer q from h+_1 (vcol_q)
  DvceArray5D<Real> vgd_hr, vgd_hrt;   // (m, 6, k, j, i): r+_r r-_r r+_t r-_t r+_p r-_p
  DvceArray1D<Real> vgd_hrm;           // shell mean of one twin component
  // gdthin-1011: with implicit_thin_switch = fs* or an implicit_thin_corr (vgd_thx) vgd_hr
  // carries M1_VGD_HRNX components: + C1 = sum w |I - S|/J, eps, J (E units)
  bool vgd_thx = false;
  bool vgd_hr_tw = false;              // stage 0 of the fused twin: S from vgd_cst
  DvceArray1D<Real> vgd_hrm6;          // vgdrot-1011: the 6 shell means at once
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_pin_hr, vgd_pin_tw, vgd_pin_cnt;
  void VetGdHalfRange(const int stage);
  DvceArray5D<Real> ifw2, ifw3;
  void ImplicitLatFaceCoef(bool first);
  void ImplicitLatFaceCoefCs(bool first);   // cs-s4-1017
  void CbDebugDump();
  void ImplicitHrSfBuild();   // cbfix-1011: the M1_NIW_HSF block (implicit_hr_face=sf)   // cbfix-1011: env M1_CBDBG=<prefix> per-step face/budget dump
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
  // moveq-1011 (moving media, NOTE-moveq.md): the radiation momentum of the enthalpy flux
  //   implicit_rad_inertia  the F0 face equation's dropped -dA/dt: after the solve
  //                         each face takes -theta_f dA_f (dA = A' - A^n, A = a E),
  //                         so the deposit dt k F0/c hands the gas -(1-theta) dA/
  //                         (c chat), plus the dt v_f g0_f/c the face equation
  //                         removed from the radiation
  //   implicit_f0_init      cold-start face flux: zero (old) | lab (F0 = F_pgen - A) |
  //                         auto (lab on faces with A != 0, else zero: static bitwise)
  //   implicit_nofb_work    gas_feedback = false: E still pays the work v.dm (fixed v)
  bool impl_rad_inertia = false;
  int impl_f0_init = 0;          // 0 zero, 1 lab, 2 auto
  bool impl_f0_init_done = false;
  bool impl_nofb_work = false;
  bool nofb_wrow = false;         // the last assembled E rows carried v.dm (nofb)
  // consrc-1011 <rad_m1>/implicit_gas_src = matched (default: the term-by-term deposit
  // above, bitwise) | conservative: the gas takes minus the radiation's own source of the
  // solve, dm = -(F'_c - F^old_c)/(c chat) - (1/2) sum_faces Z_f, with F_c the STORED
  // cell flux and Z_f the face's transport (and non-momentum-equation) part, written by
  // the face kernels of the solve (rad_m1_implicit.cpp, NOTE-consrc.md)
  bool impl_gas_src_cons = false;
  bool impl_gas_src_full = false;      // conservative_full: no drag-only HLL exception
  DvceArray4D<Real> gcz1, gcz2, gcz3;   // Z_f on the x1 / x2 / x3 faces (as f0x1..3)
  DvceArray5D<Real> csk1, csk2, csinc;  // hesdirk2: slopes / increment of the cell F (3)
  bool csk1_ok = false;                 // csk1 holds the slope that pairs with t2k1
  bool csinc_ok = false;                // csinc was formed for the coming stage solve
  bool gcs_first = true;                // the first conservative solve of this run
  // rotavg-rst-1012: csk1 (+ csk1_ok, gcs_first) in the restart file (kM1CskRstMagic):
  // the header whenever a conservative solve has run, the 3 slabs when csk1 holds a slope
  bool CskRstHdr() const {return impl_gas_src_cons && !gcs_first;}
  int CskRstNch() const {
    return (CskRstHdr() && csk1_ok && csk1.extent_int(0) > 0) ? 3 : 0;
  }
  // consfix-1011, implicit_gas_src = conservative2 (NOTE-consfix.md): conservative, plus
  //  (a) the face equation is the lab-F one, F0' = theta (F0^n - dA_f - ...): dA_c =
  //      A'_c - A^old_c the whole change of the cell's enthalpy flux (A^old = F^old -
  //      F0^old face mean), put into F0^n INSIDE the solve (ImplicitCs2Pass), so the E
  //      rows transport the F0 the faces end on.  The drag then hands the gas (1 - theta)
  //      of it (thick: all, the radiation inertia; thin: the lab flux is kept);
  //  (b) implicit_vimp predicts the velocity the conservative write-back gives
  //      (ImplicitCs2Vimp): the E rows' enthalpy flux and the stored F share one v;
  //  (c) opt-in, implicit_cs2_shift = true (default false since consfix2): the hydro's
  //      velocity change since the stored F as F0's EXPLICIT tendency -dA/dt
  //      (ImplicitCs2Shift at the solve start; the lab F is kept).  Off, it is part of
  //      the solve's dA (which gives clean order 2 on the T ladders).
  //  (d) implicit_cs2_hll = drag (default, consfix3): the AP-HLL (thin) part of a face
  //      has no face equation; there the gas pays only its drag, as matched, not dA.
  // iw components behind the other blocks: A^old (3), dA (3, exchanged), om_bar (3)
  bool impl_gas_src_cs2 = false;
  int iw_cs2 = -1;
  bool cs2_dav = false;                 // the dA components hold a pass's prediction
  bool cs2_shift = false;               // implicit_cs2_shift: (c) on (default off)
  bool cs2_hll_drag = true;             // implicit_cs2_hll = drag (default) | inertia
  bool cs2_clip_rad = false;            // cs2clip-1020 implicit_cs2_clip = rad (named)
  bool cs2_vimp = true;                 // implicit_cs2_vimp (open-ot-1013): (b) on
  int cs2_dbg_n = 0;                    // implicit_cs2_dbg bit 32: lines printed
  int cs2_vnp = 0;                      // implicit_cs2_vimp_npass (fresh default 3)
  int cs2_dbg = 0;                      // implicit_cs2_dbg (diagnostic): 1 no face split,
                                        // 2 no vimp prediction, 4 no hydro shift,
                                        // 8 no x1 split, 16 no x2/x3 split
  bool Cs2Active();
  // open-ot-1013 time2_enth_be: this (non-stage) solve takes the stage a_f rule
  bool T2EnthBe() const;
  // h2order2-1011 <rad_m1>/implicit_rad_inertia_solve (read only when named; default
  // false = the after-solve split, bitwise): the -theta dA of implicit_rad_inertia
  // (dA_c = a_c (E'_c - E^ov_c), a frozen for the solve) goes into the face equation
  // INSIDE the Picard passes (the F0^n shift of ImplicitCs2Pass, Picard-lagged by one
  // pass; the Picard stop then takes no contraction estimate), so the E rows transport
  // the F0 the faces end on and the ESDIRK slope of F0 carries it.  The after-solve
  // split leaves the E update and the stored F0 inconsistent by theta a dE: a first-
  // order (in dt) splitting error ~ v^2, suppressed in the stiff limit (rswv / NHvk0
  // floor, NOTE-h2order2).  matched / conservative only (conservative2 has it).
  bool impl_rinr_solve = false;
  bool RinrActive();
  void ImplicitRinrDeltaA();
  void ImplicitCs2DeltaA(bool usedv);
  void ImplicitCs2Vimp();
  void ImplicitCs2Shift();
  void ImplicitCs2Pass(bool restore);
  DvceArray4D<Real> cs2f1, cs2f2, cs2f3;   // the dA face means applied to F0^n (x1..x3)

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
  // hmstage-1019: HIP builds hand MPI host copies of the two buffers (see
  // ImplicitHaloMPIFinish); empty on CUDA / host builds
  Kokkos::View<Real*, Kokkos::HostSpace> hm_shost, hm_rhost;
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
  DvceArray1D<Real> gf_y;           // accel-1014: (q*n1 + gi) the layer sums of gf_solve
  DvceArray1D<Real> gf_el;          // accel-1014: banded LU of each mode (gf_inv)
  DvceArray1D<int> gf_pv;           // accel-1014: its pivots
  int acc_gfs = 1;                  // accel-1014: parallel gf_solve (env ACC_GFS=0: old)
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
  // unify-sc-gd-1012 (#10) implicit_predictor_frozen: the predictor runs wherever the
  // face coefficients of the solve do not read the iterate (sp: vfix; any mesh:
  // halfrange + implicit_thin_switch), not on the sp wedge only
  bool impl_pred_frz = false;
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
  bool t2_enth_be;              // open-ot-1013 time2_enth_be: the non-stage (backward-
                                // Euler) solves of a hesdirk2 run take the stage solves'
                                // a_f rule (time2_enth_vel = central only)
  int t2_solve;                 // what the next ImplicitSolve does (M1_T2S_*)
  bool t2_fail;                 // the last stage solve was not admissible
  int t2_dbg_fail;              // DEBUG time2_dbg_fail: fail stage 1 at this cycle
  int t2_dbg_adm;               // DEBUG dbg_t2_admiss: report the first N non-admissible
  int t2_dbg_adm_n;             // stages (which quantity <= 0, where); 0 = off
  // t2retry-1020: <rad_m1>/time2_reject = be (default) | retry (read only when named): a
  // non-admissible hesdirk2 stage raises the dt_retry flag (reason time2_reject) instead
  // of the backward-Euler redo, so the cycle is taken again at a smaller dt with hesdirk2
  int t2_reject = 0;            // 0 be, 1 retry
  Real t2_nstep, t2_nbe, t2_nfall;  // stage steps, BE steps, fallbacks (counters)
  DvceArray5D<Real> dbrow;      // DEBUG dbg_t2_admiss: the sp E row by term
  Real t2_dtprev;               // the dt of the previous step (vet_sc extrapolation)
  bool t2_vprev;                // vet_prev holds the tensor of the previous step
  bool t2_vext;                 // time2_vet_extrap (default false: D^n)
  int t2_vsmode = 0;            // time2_vet_sc: 0 lag, 1 predict, 2 rebuild
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
  // time2_start (h2order-1011; read only when named, default be = bitwise): how a step
  // without a stored slope K1 (the first step, a restart without the slope) starts.
  // be: the step is Lie split (Heun hydro, then one backward-Euler solve), first order
  // in time for that step.  slope: a backward-Euler PROBE solve of dt_p =
  // time2_start_frac x dt from U^n gives K1 = (Y_p - U^n)/dt_p = f_I(Y_p), U^n is put
  // back, and the step runs the hesdirk2 stages.  Default frac 1: K1 = f_I(U^n) + O(dt)
  // is enough for order 2 (one step), and where U^n is off the stiff manifold (a pulse
  // with F = 0, a moving gas in a radiation field at rest) Y_p is the relaxed state, so
  // (1 - g) dt K1 stays a bounded move toward it.  A small frac makes K1 ~ jump/dt_p and
  // the stage-1 old vector overshoots by dt/dt_p (h2order_1011: rsw, rlat floors).
  int t2_start = 0;             // 0 be, 1 slope
  // implicit_vimp_bc (h2order-1011, ImplicitVimpBuild; default true on fresh runs with
  // implicit_vimp since vimpbcdef-1011, false on a restart whose file lacks it): the implicit_bc_advect enthalpy flux of a Marshak x1 end face takes
  // the stage's radiative kick dv of its cell, as the interior faces do under vimp
  bool impl_vimp_bc = false;
  Real t2_start_frac = 1.0;     // time2_start_frac
  Real t2_nprobe = 0.0;         // probe solves taken (counter)
  bool t2_inprobe = false;      // inside the probe solve (ImplicitF0Init moves t2f*)
  bool Time2StartProbe(Driver *pdrive, Real dt);
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
  // time2_vet_col_source (vcpred-1010; read only when named, default relaxed = bitwise):
  // state = the predict/rebuild build takes its VET source from the predicted state
  // itself (vet_col_source = gas there), not relaxed by a further dt (rad_m1_time2.cpp)
  bool t2_vcsrc = false;
  bool t2_vcsbl = false;       // sporder-1012: time2_vet_col_source = blend
  void Time2VetSourceInit(ParameterInput *pin, const bool act, const bool col);
  // time2_pred_cap (predfix-1011): the predicted P = U^n + theta dt K1 of a predict /
  // rebuild build takes per cell theta = min(1, cap/r), r the relative change of E and
  // of the gas internal energy (0 = off, the unlimited predictor)
  Real t2_pcap = 0.0;
  Real t2_pnlim = 0.0, t2_pncell = 0.0;   // limited / all cell-builds (Time2Report)
  int t2_dbg_pred = 0, t2_dbg_pred_n = 0; // DEBUG dbg_t2_pred = N: the first N builds
  DvceArray4D<Real> t2pth;      // (m,k,j,i) theta of the current build
  DvceArray1D<Real> vet_dbgclip;  // DEBUG dbg_t2_pred: sweep intensities clipped at 0
  DvceArray4D<Real> vcol_prev;  // extrap: (m,k,j,i) f_K of the previous step
  DvceArray3D<Real> vcol_qprev; // extrap: (m,k,j) q of the previous step
  void Time2VetColAt(int which);    // vet_col build at 1 = U^n + dt K1, 2 = Y1
  void Time2VetColSaveY1();         // end of the stage-1 solve: E, T of Y1 (rebuild)
  void Time2PredTheta(const Real dt, const Real gdt, const bool dbg);  // predfix-1011
  void Time2PredDebugSc();
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
  // cs-open-1018 (cubed sphere only): a 3-component exchange (scalar, face-normal pair)
  // of iw slots c0..c0+2, for the dA / s vectors of the conservative sources
  DvceArray5D<Real> d3w, d3w_c;
  MeshBoundaryValuesCC *pbval_d3 = nullptr;
  void ImplicitCsPairHalo(int c0);
  void ImplicitCsPairHaloD3();
  // cs-open-1018: the implicit_hr_recon limiter ghosts (hr quad map across a seam)
  void ImplicitMusclSigHalo();
  MeshBoundaryValuesCC *pbval_sg = nullptr;
  DvceArray5D<Real> sgw, sgw_c;
  // cs-open-1018 <rad_m1>/implicit_am_cs (cubed sphere; default = <hydro>/am_transport):
  // the axial angular-momentum-conservative radiation-gas exchange (rad_m1_implicit.cpp,
  // after the write-back); amw = the pre-write-back gas momentum and cell F (6)
  bool impl_am_cs = false;
  int impl_am_cs_diag = 0;   // amdef-1019 diagnostic (read only when named): solves to print
  DvceArray5D<Real> amw;
  bool am_vec = false;             // + <hydro>/am_vector: the vector form
  DvceArray4D<Real> amk1, amk2, aminc;   // spin slopes (FSAL / stage 1) and increment
  void ImplicitVimpBuild();     // the Jacobian of a(v') E' for this Picard pass
  // <rad_m1>/implicit_hr_recon = dc | plm (xthinfix-1009 Fix B, rad_m1_implicit.hpp)
  bool impl_muscl = false;      // plm half-range face values (default dc: untouched)
  bool impl_muscl_led = false;  // implicit_hr_recon = led: sign-preserving linearisation
  bool impl_muscl_hq = false;   // cs-s7-1019 implicit_hr_recon_q = h: plm of H+- = h+- E
  bool impl_muscl_seam1 = false;  // cs-s7-1019 implicit_hr_recon_seam = onesided (cs)
  bool muscl_now = false;       // on for this step
  int muscl_nbuild = 0;         // limiter builds this step (then frozen)
  int impl_muscl_nfresh = 2;    // implicit_hr_recon_fresh: passes with a fresh limiter
  Real impl_muscl_damp = 0.0;   // implicit_hr_damp: beta of the consistent damping
  Real impl_muscl_qs = 0.0;     // implicit_hr_recon_qs: quasi-steady threshold eps
  bool impl_muscl_kill = true;  // implicit_hr_pos = kill (default) | floor | bound
  bool impl_muscl_bound = false;  // implicit_hr_pos = bound: per-face bound test
  Real muscl_nbound = 0.0;       // cell-axes switched to dc by the bound test
  void ImplicitMusclBound();     // the bound test of the solved iterate
  Real impl_muscl_qsrel = 0.0;  // implicit_hr_recon_qs_rel: vacuum-side threshold delta
  bool impl_muscl_qsstep = false;  // implicit_hr_recon_qs_mode = step (operator fixed)
  bool impl_muscl_pred = false;    // implicit_hr_recon_pred: dc predictor per step
  Real muscl_pred_dc = 0.0, muscl_pred_hi = 0.0;   // passes of the two phases
  int iw_muscl = -1;            // first iw component of the M1_NIW_MUSCL block
  Real muscl_nfall = 0.0;       // how often the positivity fallback dropped it
  Real muscl_emin = 1.0e300;    // the smallest E the linear solve produced with it on
  void ImplicitMusclBuild();    // the limiter and the row of this Picard pass
  void ImplicitMusclKill();     // positivity: sig = 0 around the cells with E <= 0
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
  void SetX1maxBathColumns(const DvceArray3D<Real> &eb);
  void SetX2FaceFlux(const DvceArray3D<Real> &flo, const DvceArray3D<Real> &fhi);
  //! design sect. 2: evaluate and print v_max*tau_max/chat.  Runs once, from the first
  //! Opacity task that has a filled array (tau needs rho*kappa, which the ctor has not).
  void RSLACheck();
  // ---- milestone 3a (rad_m1_implicit.cpp)
  //! read the implicit-solver parameters and allocate its arrays; a no-op in
  //! transport = explicit
  void ImplicitInit(ParameterInput *pin);
  // icrad-1011 (rad_m1_icrad.cpp): opt-in IC radiation start, fresh starts only
  int icr_mode = 0;             // <rad_m1>/ic_rad: 0 none, 1 fs, 2 relax, 3 fs_relax
                                // (default 3 when eligible, IcRadInit)
  int icr_fs_iter = 3;          // ic_rad_fs_iter: E := J passes
  int icr_relax_n = 40;         // ic_rad_relax_n: max radiation-only steps
  Real icr_relax_tol = 1.0e-4;  // ic_rad_relax_tol: stop at |dL_top/L_top| below this
  Real icr_relax_dtfac = 1.0;   // ic_rad_relax_dtfac: dt of the relaxation steps / dt
  int icr_diag = 0;             // ic_rad_diag = N: E/J shell profiles (icrad_ej.txt)
  int icr_ndiag_done = 0, icr_nrx = 0;
  bool icr_fs_pending = false, icr_fs_only = false, icr_active = false;
  bool icr_warned = false;
  Real icr_tfs = 0.0, icr_trx = 0.0;
  std::string icr_tag;
  // dic-1013 (set by a pgen, default false = nothing): at the end of IcRadRun, ONE more
  // pass at the run's dt with the radiative force on the gas, whose net momentum
  // deposit dm - dt rho a_ref per cell (x1 component) is kept in icr_dm (icr_dm_dt = the
  // dt it took); the gas momentum and kinetic energy are then restored (v as before)
  bool icr_fmeas = false;
  DvceArray4D<Real> icr_dm;
  Real icr_dm_dt = 0.0;
  void IcRadInit(ParameterInput *pin);
  void IcRadRun(Driver *pdrive);
  bool IcRadHook();
  void IcRadJ(DvceArray5D<Real> &jc);
  void IcRadProfile(const std::string &tag, DvceArray5D<Real> &jc);
  void IcRadShellSum(DvceArray5D<Real> &a, const int c0, const int nc,
                     std::vector<Real> &out);
  Real IcRadLtop();
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
  // cs-s4-1017: the planar unit normals (Cartesian) of the x2 faces (m, j = 0..n2, c)
  // and x3 faces (m, k = 0..n3, c) of each block, oriented along +x2 / +x3 (M1CsLat0)
  DvceArray3D<Real> cs_fn2, cs_fn3;
  // cs-s4-1017: the same normals at the CELL CENTRES (m, j or k incl. ghosts, c): the
  // lateral axes of the cs half-range ratios (VetGdHalfRange on cs)
  DvceArray3D<Real> cs_cn2, cs_cn3;
  void CubedS1Init();
  //! cs-s4-1017 test only (cs_vec_check): the seam halo of the tangential vector pairs
  void CubedVecCheck(ParameterInput *pin);
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
  void ImplicitX1Halo(bool eponly, int comp = -1);
  void ImplicitRadInertia();
  void ImplicitNofbWorkRow();
  void ImplicitF0Init();
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
  void ImplicitFinalResidSave(bool restore);   // picres-1017 (implicit_final_resid_log)
  // returns 0: go on (implicit_final_resid_stop failed, this pass continues), 1: accept,
  // 2: failed with no pass left (non-converged); 1 and 2 put the saved state back
  int ImplicitFinalResid(bool gnw, int igb, int igr, int igy, int it);
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
  // vetsrc-1020: vet_source = relaxed_esc (vcol_srelax true as well): the relaxed source
  // with the cell's photon escape (VetRelaxedSource qe, VetEscapeRate; rad_m1_vetsrc.hpp)
  bool vcol_sesc = false;
  bool vcol_sbl = false;       // sporder-1012: this build blends the relaxed source
  Real vcol_sx0 = 1.0;         // sporder-1012: time2_vet_col_source_x0 (blend)
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
  // vgd-latuni-1017: the Athena++ lateral ratio keys <mesh>/x2rat, x3rat named != 1 in
  // the input (AthenaK ignores them: its theta/phi grid is uniform); VetLatInit refuses
  // them under vet_col_lat / vet_gd, so a ported deck cannot think it stretched theta
  std::string vlat_rat_keys;
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
  Real vlat_cut_thin = 0.0;    // vet_col_lat_cut_thin (unify-sc-gd-1012): no cut above a
                               // cell thinner than this (0 off)
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
  // cs-s4-1017: the cubed-sphere forms (M1CsLat0 / M1CsLatF, Cartesian D_r,lat)
  void VetLatOpCs(int xc, int yc, Real sgn);
  void VetLatStencilAddCs();

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
  // vgdrot2-1012: vet_gd_rotate_dt > 0 (read only when named): the rotation in PHYSICAL
  // time instead (block b of the step = floor((t + dt/2 - t0)/dt_rot), the same angle
  // at the same t on every dt), vet_gd_rotate_t0 its origin; vet_gd_rotate_align (default
  // true): Mesh::NewTimeStep lands a step on every block boundary (VetGdAlignDt)
  Real vgd_rot_dt = 0.0, vgd_rot_t0 = 0.0;
  Real vgd_rot_stretch = 0.0;   // vet_gd_rotate_stretch: VetGdAlignDt may lengthen dt by <= this
  // defaults-1012: vet_gd_rotate_dt = auto (the fresh vet_gd default): resolved by the
  // first VetGdAlignDt to vgd_rot_auto_fac x that dt and recorded in pin_rot_
  bool vgd_rot_auto = false;
  Real vgd_rot_auto_fac = 20.0;
  ParameterInput *pin_rot_ = nullptr;
  bool vgd_rot_align = false;
  double vgd_rot_nalign = 0.0;   // steps shortened by VetGdAlignDt (diagnostic)
  // vet_gd_rotate_slope (with vet_gd_rotate_dt; probe | be | keep): the hesdirk2 FSAL
  // slope K1 of the first step of a new rotation block belongs to the OLD set (an O(dt)
  // stage-1 error per jump = first order); probe: a time2_start=slope probe solve of the
  // new set gives K1 (VetGdRotStep), be: that step is backward Euler, keep (default): as
  // before.  probe / be measured WORSE than keep (vgdrot2 NOTE: the BE probe slope is not
  // f(U^n) in stiff cells; be is Lie split): opt-in diagnostics only
  int vgd_rot_slope = 0;
  int64_t vgd_rot_lastblk = INT64_MIN;
  bool t2_rotprobe = false;
  double vgd_rot_nprobe = 0.0;
  void VetGdRotStep();
  // rotslope-1012, vet_gd_rotate_slope = correct (3): on the first step of a new block a
  // CORRECTION solve (Time2RotCorrect) of size dt_c = vet_gd_rotate_slope_frac x dt from
  // U^n with the new set and the old vector U^n - dt_c K1 turns the stored slope into
  // K1 + (I - dt_c A_new)^{-1} (f_new - f_old)(U^n): the closure jump of the slope in
  // the non-stiff cells, and keep's K1 in the stiff ones (where keep is 2nd order).
  // vet_gd_rotate_slope_reuse (default true): the stage-1 solve of that step keeps the
  // correction solve's closure build (the same set, P = U^n + theta dt K1) instead of a
  // second sweep
  bool t2_rotcorr = false;
  bool t2_vcreuse = false;
  bool vgd_rot_creuse = true;
  Real vgd_rot_cfrac = 0.0;
  double vgd_rot_ncorr = 0.0, vgd_rot_ncfail = 0.0, vgd_rot_ncpass = 0.0;
  bool Time2RotCorrect(Driver *pdrive, Real dt);
  Real vgd_alpha = -1.0;       // the current z-angle (-1: tables not built)
  // cs-s3b-1017 (USER decision 1): <rad_m1>/vet_gd_rotate_axis = z | 3d (default z on
  // the sp wedge: bitwise; 3d on the cubed sphere).  3d: the set of rotation block b is
  // R_b n with R_b a 3-D rotation about a quasi-random axis (VgdRot3); vgd_alpha is then
  // the block KEY (zig-zag code of b, or of b mod K with vet_rotate_avg = K), not an
  // angle
  bool vgd_rot3 = false;
  void VgdRot3(const Real key, double r[9]) const;
  // rotavg-1012: <rad_m1>/vet_rotate_avg = K >= 2 (read only when named; default 0 = off,
  // bitwise): the closure handed to the moment step is the mean over the LAST K rotated
  // direction sets, not the current set's alone (rad_m1_rotavg.cpp).  Closure-agnostic:
  // a closure registers its NORMALISED hand-over fields (ratios to J: the Eddington
  // tensor parts, the half-range ratios, C1/J; J itself as a ratio), and the average is
  // kept as per-set OFFSETS d_j = g_{b-j} - g_b of each older set from the current one,
  // every difference taken from two sweeps of the SAME state (the old set re-swept once
  // at the first build of a new block): used = g_b(now) + sum_j d_j / K.  The current set
  // is live (no lag of the closure); only the ray-pattern differences are carried
  // between blocks.  On a fresh start or a restart the K-1 older sets are swept once at
  // the first build.  Needs a rotating closure (vet_gd: vet_gd_rotate_dt / _every).
  static constexpr int M1_RA_NF = 16;
  struct RotAvgTab {           // the registered fields (POD: captured by the kernels)
    int nf = 0;
    int arr[M1_RA_NF], comp[M1_RA_NF], kind[M1_RA_NF];
    Real lo[M1_RA_NF], hi[M1_RA_NF];
    int farr = -1, fcomp = 0;  // validity flag field: invalid where < 0 (-1: none)
    int jf = -1;               // the J field (index into the table) of kind-3 fields
    int jw = -1;               // the ratio field weighting the sets (-1: equal)
  };
  int vra_k = 0;               // vet_rotate_avg (K; 0/1 off)
  bool vra_eq = false;         // vet_rotate_avg_weight = equal (default j: J-weighted)
  bool vra_on = false;
  RotAvgTab vra_tab;
  DvceArray5D<Real> vra_a0, vra_a1;   // the closure's field arrays (set per call)
  DvceArray3D<Real> vra_c;     // a column field (m, k, j) at the top cell (arr 2)
  DvceArray5D<Real> vra_d;     // (m, (K-1)(nf+1), k, j, i) interior: offsets + valid
  int vra_head = 0;            // the ring slot of d_1 (the previous set)
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
  // rotavg-rst-1012: the history in the restart file (kM1RaRstMagic, rad_m1_rotavg.cpp)
  static constexpr int M1_RA_RST_NH = 11;   // int64 header entries (rotavg-rst-1012: 10)
  int RotAvgNf() const;
  int RotAvgRstNslab(const int ncell) const;
  int RotAvgRstI0() const;     // rotavg-rst2-1013: the first written radial cell
  void RotAvgRstHdr(int64_t *h) const;
  void RotAvgRstPack(DvceArray5D<Real> &dst, const int nmb);
  void RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb, const int nsc);   // mem-1014
  bool RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src, const int nmb);
  bool RotAvgRstCheck();
  int vra_rst_nf = -1;         // nf of a history loaded from the file (-1: none/checked)
  int64_t vra_blk = INT64_MIN; // the block of the last build
  bool vra_last = false;       // this build is on the last step of its block
  DvceArray5D<Real> vra_bk_tt, vra_bk_hr;   // vet_gd: the hand-over fields kept over
  DvceArray3D<Real> vra_bk_q;               // the next-set sweep of a last build
  DvceArray5D<Real> vra_bk_vc;              // vet_sc: the same (vet_cell)
  void VetGdRotAvgSetup();     // vet_gd: registers tau_ten LAT + vgd_hr
  void VetGdSweepPost();       // twin, sweep(s), fused moments, VetGdPost
  // sptests-1011: vet_gd_inner = emit (read only when named; default diff = bitwise): with
  // the first gd shell at the physical inner face (icut = 0) the outward rays start from
  // the inner boundary's own emission (flux: black sphere I = 4 F_in/c; marshak: the bath;
  // vet_gd_beam: one set direction from a patch) instead of the diffusion intensity
  bool vgd_emit = false;
  bool vgd_tbath = false;   // vgdbath-1011: the top (x1max Marshak) bath enters the sweep
  bool vgd_beam = false;       // vet_gd_beam (needs emit)
  int vgd_nbeam = 0;           // 1, or 2 with the vet_gd_beam2_* keys
  int vgd_bdir[2] = {-1, -1};  // each beam's set direction
  Real vgd_bib[2] = {0.0, 0.0};    // its intensity e_b / w_d
  Real vgd_bfp[2] = {0.0, 0.0};    // patch fraction of the inner sphere (twin inflow)
  Real vgd_bth[2][2] = {{0.0, 0.0}, {0.0, 0.0}};   // the patches (theta, phi ranges)
  Real vgd_bph[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  Real vgd_bang[2][2] = {{0.0, 0.0}, {0.0, 0.0}};  // target directions (theta, phi)
  Real vgd_be[2] = {1.0, 1.0};     // vet_gd_beam_e: E of each beam
  DvceArray3D<Real> vgd_bfx;   // (m, k, j): the matching inner face flux of the M1 solve
  bool vgd_tws = false;        // the unfused twin's sweep is running (its inflow)
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
  // mem2-1014 vet_gd_band_ring (default true where it applies: MPI band with the exact
  // halo lists, no overlap / async / fused shells; memcut-1020: on the cubed sphere with
  // the cs compact lists, vet_gd_halo_compact > 0; read when named; env VGD_RING=0 off):
  // the ragged array keeps every shell's INTERIOR only (compact slab, depth 0) and lends
  // a shell one of VGD_RING_K wide slabs (its full band depth wi, the same rectangle as
  // before) only while the sweep needs its band: inward pass from just before the shell
  // is computed until the next shell has read it; outward pass from just before the
  // shell is computed (its interior copied in, and the band values of its INWARD
  // exchange, saved per shell after that exchange, put back: the turned rays read them)
  // until the next shell has read it.  The accessor (VgdRag::operator()) is unchanged:
  // only the device tables off(m, i) / wi(i) are switched (in stream order), so the
  // shell kernel's code and every value it reads are those of the dense-band array
  // (bitwise).  One extra small kernel per shell (copies), launched while the shell's
  // halo messages travel.
  static constexpr int VGD_RING_K = 3;
  struct VgdRingOp {
    int kind = 0;      // 1 save (inward band -> store), 2 release, 3 acquire
    int i = -1, slot = -1;
    bool copyin = false, restore = false;
  };
  bool vgd_ring = false;              // the ring layout is allocated (both arrays)
  bool vgd_ring_named = false, vgd_ring_par = true;
  bool vgd_ring_act = false;          // the ring runs in this sweep
  bool vgd_ring_poison = false;       // env VGD_RING_POISON=1: NaN in lent slabs (test)
  bool vgd_rg_two = false;            // this sweep carries vgd_itw too
  bool vgd_rg_lastx = false;          // the last halo used an exact list
  int vgd_rg_lastx_e = -1;
  std::vector<int64_t> vgd_rg_offc;   // (m c1 + i): compact slab offsets (host)
  std::vector<int64_t> vgd_rg_soff;   // (s nmb + m): slot offsets (host)
  std::vector<int> vgd_rg_wfull;      // (i): a shell's band depth when lent a slot
  std::vector<int> vgd_rg_state;      // (i): -1 compact, else slot
  std::vector<int> vgd_rg_sown;       // (s): the shell holding slot s, -1 free
  std::vector<int64_t> vgd_rg_so, vgd_rg_sn;   // (i): save store offset, count/array
  std::vector<int> vgd_rg_sk;         // (i): save kind: 1 exact list, 2 dense band
  std::vector<int> vgd_rg_se;         // (i): the exact list record of the save
  int64_t vgd_rg_stop = 0;            // store fill this sweep
  DvceArray1D<Real> vgd_rg_store;     // saved inward band values (main, then twin)
  Kokkos::View<int64_t**, DevMemSpace> vgd_rg_doffc;   // (m, i)
  Kokkos::View<int64_t**, DevMemSpace> vgd_rg_dsoff;   // (s, m)
  std::vector<VgdRingOp> vgd_rg_pend;
  Real vgd_rg_nk = 0.0, vgd_rg_ndense = 0.0;
  void VetGdRingAlloc();
  void VetGdRingFlush();
  void VetGdRingSweepBegin(const bool two);
  void VetGdRingSweepEnd();
  int VetGdRingFreeSlot(const int avoid) const;
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
  // twinfix-1012: vet_gd_twin_state = trans (the twin's shell state from the shell's mean
  // TRANSMISSION, not its mean opacity; vgd_twtr) and vet_gd_twin_weight = j (the twin
  // subtraction weighted per cell by the main/twin J agreement; vgd_twwt, vgd_twwq)
  bool vgd_twtr = false, vgd_twwt = false;
  Real vgd_twwq = 0.5;
  DvceArray1D<Real> vgd_twtb;  // (5*c1): sum t, sum (1-t), sum S(1-t), min, max ln chi
  void VetGdTwinTrans();
  DvceArray5D<Real> vgd_twl, vgd_cs0;
  // vgdlat-1012 vet_gd_symref = col (read only when named, default none; needs
  // vet_gd_twin): the twin's laterally symmetric reference (its shell mean) replaced by
  // the 1-D spherical p-ray solution of the same shell-mean state (VetGdSymRef)
  bool vgd_symref = false, vgd_sr_ok = false;
  int vgd_sr_nc = 64, vgd_sr_ns = 4;
  std::vector<Real> vgd_sr_chi, vgd_sr_s, vgd_sr_x1v, vgd_sr_x1f;
  DvceArray1D<Real> vgd_srm;   // (q c1 + i), q = M1_VGD_SR_*
  double vgd_tsr = 0.0, vgd_tsrj = 0.0;
  bool vgd_sr_pend = false;
  std::vector<double> vgd_sr_h;
  std::future<std::vector<double>> vgd_sr_fut;
  void VetGdSymRef();
  void VetGdSymRefJoin();
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
  void VetGdMmsCs();   // cs-s4-1017 (rad_m1_vetgd_cs.cpp)
  // cs-s4-1017: the cs half-range data at the start of a solve (|H|/J, seam ghosts)
  void VetGdCsHrPrep();
  MeshBoundaryValuesCC *pbval_hr = nullptr;
  DvceArray5D<Real> vgd_hr_c;
  void VetGdRealDiag();
  int vgd_rdiag = -1, vgd_rdcnt = 0;   // env VGD_TIME_HALO=1: fence before the halo timer
  void VetGdHaloInit();
  // cs-vetgd-1017 (CS_VETGD.md sect. 4A, stage S3; rad_m1_vetgd_cs.cpp): vet_gd on the
  // cubed sphere.  The band of a block is its OWN panel's index space (as on sp), and
  // the strips of the four neighbour panels are stored UNFOLDED next to the panel edges
  // (a neighbour cell at depth d from the seam sits where the mirror image of the own
  // cell at depth d would continue the index line; equal N, so the pairing is one to
  // one).  Every band position holds one global cell (or none: the corner quadrants
  // beyond two edges, depths >= N), filled per shell by exact keyed copies (on-rank +
  // MPI).  The upwind point is found in Cartesian, its panel by the largest normal
  // component and its (xi, eta) by the gnomonic chart; reads are bilinear in that
  // panel's native indices, the seam strip (half a cell on either side) linear across
  // the seam with each side's own tangential coordinate, the cube vertex triangle by
  // spherical barycentric weights of its 3 cells.  No walls, no direction maps.
  bool vgd_csg = false;               // vet_gd on the cubed sphere
  int vgd_csn = 0;                    // cells per panel edge
  DvceArray4D<Real> vgd_csr;          // (m, k, j, 3): band cell unit vector, 0 none
  std::vector<double> vgd_csr_h;      // the same on the host (double), (m, k, j, 3)
  DvceArray2D<int> vgd_csblk;         // (m, 3): panel, global j0, global k0
  DvceArray3D<int> vgd_cscm;          // (p, q, 7): q (a,t) -> p unfolded (J,K) affine
  DvceArray2D<int> vgd_csnb;          // (p, 4): panel across side -x2, +x2, -x3, +x3
  DvceArray2D<int> vgd_csloc;         // (e, 6): on-rank copy dst (m,k,j) <- src (m,k,j)
  DvceArray2D<int> vgd_cssl, vgd_csrl;   // (e, 3): MPI send src / recv dst (m,k,j)
  std::vector<int> vgd_csprk;         // partner ranks
  std::vector<int> vgd_csps, vgd_cspn, vgd_cspr, vgd_cspm;   // send/recv displ, count
  // message buffers.  xchg-1021: in PINNED HOST memory (the pack/unpack kernels write/
  // read them directly), so MPI never sees a GPU pointer here: no GPUDirect RDMA and
  // no UCX cuda staging for the vet_gd band (pure copies: bitwise the device buffers)
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_cssb, vgd_csrb;
  int vgd_csnloc = 0, vgd_csnsnd = 0, vgd_csnrcv = 0;
  void VetGdCsInit();
  void VetGdCsSource();
  void VetGdCsSweep();
  void VetGdCsMoments(const int part);   // 0: D_rr -> LAT0, 1: surface q
  void VetGdCsTab(const std::vector<double> &d);   // into vgd_tabh (walls off, mr)
  void VetGdCsDump(const char *env);
  template <class V>
  void VetGdCsHalo(V &a, const int nv, const int i0, const int i1, const bool loc);
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_csown;   // (m, k, j, 3) owner
  // cs-s3b-1017 (USER decision 2): vet_gd_halo_compact > 0 on cs -- the chart-aware
  // exact lists: per shell and pass only the (remote band cell, direction) values the
  // sweep reads (VetGdCsCompact: a read-mask kernel that replays the sweep's upwind
  // stencils for the direction set, the receivers' masks sent to the owners once per set,
  // both sides' lists in the same order).  Bitwise the dense band copies (the same values
  // are read; every remote read is checked against the mask in the sweep: FATAL on a
  // miss)
  struct VgdCsCl {
    Real alpha = -9.0;
    int64_t use = 0;
    std::vector<int64_t> so, ro;     // ((pass*n1 + s)*(np + 1) + p): list offsets
    Kokkos::View<int *, LayoutWrapper, DevMemSpace> si, ri;   // e*n + d per value
    // (e_rcv, d, s): 1 in, 2 out, 0 not read; memcut-1020: 2 bits per value, 16 values
    // per word (value t = (e_rcv*n + d)*n1 + s in word t/16, bits 2 (t%16))
    Kokkos::View<uint32_t *, LayoutWrapper, DevMemSpace> msk;
  };
  bool vgd_cscomp = false;
  VgdCsCl vgd_cscl[2];
  int vgd_cscur = 0;
  int64_t vgd_csuse = 0;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_csrix;   // (m, k, j): e_rcv / -1
  Kokkos::View<int *, LayoutWrapper, DevMemSpace> vgd_csmiss;   // reads not in the mask
  double vgd_csnbuild = 0.0, vgd_cstbuild = 0.0, vgd_csvsent = 0.0, vgd_csvdense = 0.0;
  void VetGdCsCompact();
  void VetGdCsHaloC(const int pass, const int i);
  // vet_gd_halo_compact = 0 (off), 1 (pass branch), 2 (+ ray reach): compact messages
  int vgd_hcomp = 0;
  bool vgd_capped = false, vgd_allow_clamp = false;
  Real vgd_nclamp_seen = 0.0;
  Real vgd_npath = 0.0;   // open-od-1013: vet_path subdivided steps (ALL ranks, all sweeps)
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
  bool VetGdHaloCompact(V &a, const int nv, const int i0, const int ws, V *b = nullptr);
  // vet_gd_halo_list_mb (vgdspeed-1009; read only when named; default 4096, 0 = off):
  // MB budget per rank of the compact-halo LISTS per (pass, shell), built from the masks
  // at the first exchange of a mask set (direction set, cut, depth) and used by every
  // later sweep of the set: gather pack, MPI, scatter (+ zeros on the other pass's
  // entries) instead of flag / scan / pack / expand / dense unpack.  Bitwise: the band
  // gets the same values (see VetGdHaloCompact for when the dense zeros are implied).
  struct VgdHlEntry {
    DvceArray1D<int> sl;          // (ns): packed band address of each sent value
    DvceArray1D<int> dl, sr;      // (nr): packed band address, compact source position
    DvceArray1D<int> dw, sw;      // (nw), (3 nw): wall_interp ghosts: address, sources
    std::vector<int> sp, rp;      // (np + 1): partner bounds of the compact messages
    std::vector<int> dsto;        // vgdrot-1011 (np): start of my values in each
                                  // partner's receive layout (batch-made exact lists)
    Real alpha = -2.0;
    int scut = -1, ws = -1, ns = 0, nr = 0, nw = 0;
    int64_t gen = -1;             // build id (> 0) while valid
    bool nofit = false;           // over the budget for (alpha, scut, ws)
    size_t bytes = 0;
    // vgdfuse-1010 vet_gd_halo_exact: the lists hold only the values the sweep READS
    // (sr / sw = -1: a read ghost whose source is not in the mask: written 0)
    bool exact = false;
  };
  int vgd_hl_mb = 4096;
  // vet_gd_shell_list (vgdspeed-1009; read only when named; default true; bitwise): the
  // shell kernel runs over the (m, k, j, d) of the pass's branch only (VetGdSweep)
  bool vgd_shl_on = true;
  DvceArray1D<int> vgd_shl[2];
  int vgd_shn[2] = {0, 0};
  Real vgd_shl_alpha = -3.0;
  std::vector<VgdHlEntry> vgd_hle;
  size_t vgd_hl_bytes = 0, vgd_hl_cap = 0;   // cap: bytes (set at the first build)
  int64_t vgd_hl_gen = 0;
  Real vgd_hl_nmade = 0.0, vgd_hl_nuse = 0.0, vgd_hl_ndense = 0.0;
  DvceArray1D<int> vgd_hl_fn, vgd_hl_fw;    // build scratch: dest flags -> positions
  // per band array (data pointer) and pass: the build id of the entry whose exchange
  // last wrote each shell's band (0 = never written: all zero; -1 = unknown)
  std::map<const void *, std::array<std::vector<int64_t>, 2>> vgd_hl_gw;
  bool VetGdHlBuild(const int e, const int i, const bool inw, const int ws,
                    const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                    const std::vector<int> &pbh);
  // vgdfuse-1010 (all bitwise unless stated)
  // vet_gd_overlap = G (read only when named; default 0 = off): split each shell launch
  // into the items that read no remote band of the previous shell (run while that
  // shell's list-path messages travel) and the rest (run after the scatter).  The split
  // lists are made per group of G shells (an item is "remote" if it is in any shell of
  // the group) when the direction set changes; vet_gd_overlap_mb caps their memory.
  int vgd_ovl_g = 0, vgd_ovl_mb = 1024;
  DvceArray4D<Real> vgd_trig;   // (m, k, j, 4): sin, cos of theta, phi (shell kernel)
  DvceArray1D<int> vgd_ovl[2];             // (ngroup cap): per group remote, then local
  DvceArray1D<int> vgd_ovl_f;              // scratch: flags -> positions
  std::vector<int> vgd_ovl_nr[2];          // per group: number of remote items
  std::vector<Real> vgd_ovl_alpha[2];      // per group: direction set of the lists
  int vgd_ovl_cap[2] = {0, 0};
  int vgd_ovl_cut = -1;                    // the cut the groups were made for
  bool vgd_ovl_nofit = false;
  Real vgd_ovl_nsplit = 0.0, vgd_ovl_nfull = 0.0, vgd_ovl_nrem = 0.0, vgd_ovl_nall = 0.0;
  // vet_gd_fuse_shells = H (read only when named; default 1): the lateral exchange only
  // after every H-th shell of a pass (and its last); the other shells read the band of
  // the last exchange that wrote it.  NOT bitwise for H > 1 (a lagged inflow)
  int vgd_fuse_h = 1;
  Real vgd_fuse_nskip = 0.0;
  // env VGD_SCAN_SPLIT=1: the paired mask / list scans in two launches (as before 1010)
  bool vgd_scan2 = true;
  Kokkos::View<int64_t*, DevMemSpace> vgd_scan_sums;   // VgdScan2 chunk sums
  // vet_gd_halo_exact (read only when named; default false; bitwise): the compact-halo
  // lists of each (pass, shell) cut down to the band values the sweep actually reads
  // (marked from the sweep's own ray geometry, VetGdSweep), the senders told by one
  // request message per list build; no zeros of the other pass's entries needed.
  bool vgd_hx_on = false;
  int vgd_hx_mode = 0;     // 1: masked (bitwise to the dense path), 2: mask-free
  bool vgd_hx_named = false;
  DvceArray1D<int> vgd_hx_mk;              // (m, v, k, j) dense band: read marks
  int vgd_hx_stamp = 0;                    // the marks of the current build
  Kokkos::View<int, DevMemSpace> vgd_hx_oc;   // marked ghosts outside the region
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_pin, vgd_hx_pin2;
  DvceArray1D<int> vgd_hx_need, vgd_hx_nc, vgd_hx_ns;   // build scratch
  DvceArray1D<int> vgd_hx_t0, vgd_hx_t1;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_hb;   // bounds, counts
  int vgd_hx_tag[4] = {-1, -1, -1, -1};    // the marks: pass, shell, scut, ws
  Real vgd_hx_alpha = -7.0;
  Real vgd_hx_nmade = 0.0, vgd_hx_nout = 0.0, vgd_hx_nfall = 0.0, vgd_hx_nun = 0.0;
  Real vgd_hx_vmask = 0.0, vgd_hx_vexact = 0.0;   // received values: mask lists, exact
  std::vector<VgdHlEntry> vgd_hxe;         // (2 n1): the exact lists per (pass, shell)
  bool VetGdHxKey(const int e, const int ws) const;   // record e matches the masks' key
  void VetGdHxBuild2(const int e, const int i, const int ws);
  // vgdrot-1011: the exact lists of many shells in a few collectives (VetGdHxBatch);
  // per-slot scratch pools, the per-depth partner / piece tables, pinned bounds
  bool vgd_hxb_on = true;
  bool vgd_hxb_cub = true;      // env VGD_HXB_CUB=0: the VgdScan2 scans
  DvceArray1D<char> vgd_hxb_tmp;
  int vgd_hxb_g = 0, vgd_hxb_npc = 0;
  size_t vgd_hxb_rcap = 0, vgd_hxb_scap = 0;
  DvceArray1D<int> vgd_hxb_nd, vgd_hxb_f0, vgd_hxb_f1, vgd_hxb_rq, vgd_hxb_so, vgd_hxb_oc;
  Kokkos::View<int**, LayoutWrapper, DevMemSpace> vgd_hxb_pst, vgd_hxb_pcs, vgd_hxb_pci,
                                                  vgd_hxb_pds;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hxb_hb, vgd_hxb_sp;
  Real vgd_hxb_t = 0.0, vgd_hxb_ngroup = 0.0, vgd_hxb_nrealloc = 0.0;
  void VetGdHxBatch(const std::vector<std::array<int, 4>> &jobs);
  // vgdrot-1011 vet_gd_halo_ipc (read only when named; default false; bitwise): the
  // per-shell exchange of the batch-made exact lists over node-local CUDA IPC, driven by
  // the GPU (spin / flag kernels on the sweep's stream; the host never waits in the
  // sweep).  env VGD_IPC=0/1 overrides.
  bool vgd_ipc_on = false;
  bool vgd_ipc_tfs = false;
  bool vgd_merge = true;         // bsgaccel-1011: env VGD_MERGE=0: twin in own messages
  DvceArray1D<Real> vgd_csbm, vgd_crbm;   // bsgaccel-1011: merged send / receive
  bool vgd_ipc_lean = true;      // bsgaccel-1011: env VGD_IPC_LEAN=0: the ack kernel k1
  int vgd_ipc_mode = 1;          // env VGD_IPC=2: local pack + peer copies
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
  void VetGdIpcInit();
  void VetGdIpcCheck();
  void VetGdHxBuild(const int e, const int i, const bool inw, const int ws,
                    const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                    const std::vector<int> &sp, const std::vector<int> &rp);
  // the pending list-path exchange of the overlap (VetGdHlBegin / VetGdHlEnd)
  struct VgdHlPend {
    bool on = false, two = false, xe = false;   // xe: X is an exact (vgd_hxe) entry
    bool ipc = false;                           // vgdrot-1011: the IPC path
    bool mrg = false;                           // bsgaccel-1011: merged twin messages
    bool nccl = false;                          // accel-1014: NCCL (stream-ordered) path
    long long seq = 0;
    int ex = -1, px = 0, i = -1;
    std::array<std::vector<int64_t>, 2> *gwa = nullptr, *gwb = nullptr;
    VgdRag a, b;
#if MPI_PARALLEL_ENABLED
    std::vector<MPI_Request> req;
#endif
    Kokkos::Timer tq;
  };
  VgdHlPend vgd_hlp;
  // accel-1014 (compile -DM1_VGD_NCCL, env ACC_NCCL=1): the per-shell list exchange of
  // the MPI path as stream-ordered NCCL send/recv on the sweep's stream (no host fence and
  // no MPI wait per shell; the same values reach the same band cells)
  int acc_nccl = -1;            // -1 not decided, 0 off, 1 on
  void *acc_ncomm = nullptr;    // ncclComm_t
  bool VgdNcclOn();
  template <class V>
  bool VetGdHlBegin(V &a, const int nv, const int i0, const int ws, V *b);
  void VetGdHlEnd();
  std::vector<int> vgd_wsh;    // (i): band depth the shell's data needs (<= vgd_w)
  std::vector<int> vgd_wsi, vgd_wso;   // (i): the same per pass (inward, outward)
  template <class V>
  void VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                 const bool mapd, V *b = nullptr);
  // mem-1014: the multi-shell (ln chi, ln S) halo in chunks of shells that fit the
  // per-shell message buffers (the same values; the buffers sized for one shell of all
  // directions instead of all shells of the 2 source fields)
  void VetGdHaloCs(DvceArray5D<Real> &a, const int i0, const int i1);
  Real vgd_tsrc = 0.0, vgd_tswp = 0.0, vgd_texc = 0.0, vgd_tmom = 0.0;
  void VetGdInit();
  Real VetGdAngle(const int cyc) const;
  bool VgdRotOn() const {return vgd_rot > 0 || vgd_rot_dt > 0.0 || vgd_rot_auto;}
  int64_t VetGdBlock(const int cyc) const;       // vgdrot2-1012: the rotation block
  Real VetGdAngleBlk(const int64_t b) const;
  double RotPhase(const int64_t b) const;   // rotavg-1012: block -> fraction of period
  Real VetGdAlignDt(const Real time, const Real dt);
  void VetGdTables(const Real alpha);
  // vgdrot-1011: the tables' host part in persistent pinned buffers (no per-rotation
  // mirrors), and the next rotation's tables prefetched (host thread + upload into a
  // second device set on vgd_tex, swapped in at the rotation); env VGD_TAB_PREFETCH=0 off
  struct VgdTabHost {
    Kokkos::View<Real **, LayoutWrapper, Kokkos::SharedHostPinnedSpace> dir;
    Kokkos::View<int ***, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w;
    Kokkos::View<int ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mp;
    Kokkos::View<Real ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mr;
    Kokkos::View<int *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> m3;
    Kokkos::View<Real *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w3;
    std::vector<double> mt;
    Real alpha = -9.0;
    int nwall = 0, nexact = 0, nmap = 0;
  };
  struct VgdTabDev {
    DvceArray2D<Real> dir;
    Kokkos::View<int ***, LayoutWrapper, DevMemSpace> wall;
    Kokkos::View<int ****, LayoutWrapper, DevMemSpace> map;
    DvceArray4D<Real> mr;
    Kokkos::View<int *****, LayoutWrapper, DevMemSpace> m3;
    DvceArray5D<Real> w3;
  };
  VgdTabHost vgd_tabh;
  VgdTabDev vgd_tabd;
  std::future<void> vgd_tabf;
  bool vgd_tab_pre = true;
  int vgd_tab_state = 0;       // 0 idle, 1 host thread running, 3 uploaded (vgd_tabd)
  Real vgd_tab_nswap = 0.0;
  DevExeSpace vgd_tex;
  void VgdTabCompute(const Real alpha, const int nthr);
  void VgdTabAlloc(const bool dev);
  void VetGdTabPrefetch(const int cyc);
  void VetGdTabPoll();
  void VetGdBeamInit(const std::vector<double> &d);   // sptests-1011: vet_gd_beam
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
